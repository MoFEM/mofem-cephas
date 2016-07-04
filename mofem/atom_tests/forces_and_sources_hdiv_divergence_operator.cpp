/* This file is part of MoFEM.
 * MoFEM is free software: you can redistribute it and/or modify it under
 * the terms of the GNU Lesser General Public License as published by the
 * Free Software Foundation, either version 3 of the License, or (at your
 * option) any later version.
 *
 * MoFEM is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
 * License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with MoFEM. If not, see <http://www.gnu.org/licenses/>. */

#include <MoFEM.hpp>

//#include <DirichletBC.hpp>
//#include <PostProcOnRefMesh.hpp>

#include <Projection10NodeCoordsOnField.hpp>

#include <boost/iostreams/tee.hpp>
#include <boost/iostreams/stream.hpp>
#include <fstream>
#include <iostream>

#include<moab/Skinner.hpp>

namespace bio = boost::iostreams;
using bio::tee_device;
using bio::stream;

using namespace MoFEM;

static char help[] = "...\n\n";

int main(int argc, char *argv[]) {

  ErrorCode rval;
  PetscErrorCode ierr;

  PetscInitialize(&argc,&argv,(char *)0,help);

  moab::Core mb_instance;
  Interface& moab = mb_instance;
  int rank;
  MPI_Comm_rank(PETSC_COMM_WORLD,&rank);

  //create one tet
  double tet_coords[] = {
    0,0,0,
    0.5,0,0,
    0,0.5,0,
    0,0,0.5
  };

  EntityHandle nodes[4];
  for(int nn = 0;nn<4;nn++) {
    rval = moab.create_vertex(&tet_coords[3*nn],nodes[nn]); CHKERRQ_MOAB(rval);
  }

  EntityHandle tet;
  rval = moab.create_element(MBTET,nodes,4,tet); CHKERRQ_MOAB(rval);

  ParallelComm* pcomm = ParallelComm::get_pcomm(&moab,MYPCOMM_INDEX);
  if(pcomm == NULL) pcomm =  new ParallelComm(&moab,PETSC_COMM_WORLD);

  //create MoFEM (Joseph) database
  MoFEM::Core core(moab);
  FieldInterface& m_field = core;
  //meshset consisting all entities in mesh
  EntityHandle root_set = moab.get_root_set();

  //set entitities bit level
  BitRefLevel bit_level0;
  bit_level0.set(0);
  EntityHandle meshset_level0;
  rval = moab.create_meshset(MESHSET_SET,meshset_level0); CHKERRQ_MOAB(rval);
  ierr = m_field.seed_ref_level_3D(0,bit_level0); CHKERRQ(ierr);

  //fields
  ierr = m_field.add_field("HDIV",HDIV,1); CHKERRQ(ierr);
  //add entities to field
  ierr = m_field.add_ents_to_field_by_TETs(root_set,"HDIV"); CHKERRQ(ierr);
  //set app. order
  int order = 4;
  ierr = m_field.set_field_order(root_set,MBTET,"HDIV",order); CHKERRQ(ierr);
  ierr = m_field.set_field_order(root_set,MBTRI,"HDIV",order); CHKERRQ(ierr);
  //build field
  ierr = m_field.build_fields(); CHKERRQ(ierr);

  //finite elements
  ierr = m_field.add_finite_element("TET_FE"); CHKERRQ(ierr);
  ierr = m_field.add_finite_element("SKIN_FE"); CHKERRQ(ierr);

  //Define rows/cols and element data
  ierr = m_field.modify_finite_element_add_field_row("TET_FE","HDIV"); CHKERRQ(ierr);
  ierr = m_field.modify_finite_element_add_field_col("TET_FE","HDIV"); CHKERRQ(ierr);
  ierr = m_field.modify_finite_element_add_field_data("TET_FE","HDIV"); CHKERRQ(ierr);
  ierr = m_field.modify_finite_element_add_field_row("SKIN_FE","HDIV"); CHKERRQ(ierr);
  ierr = m_field.modify_finite_element_add_field_col("SKIN_FE","HDIV"); CHKERRQ(ierr);
  ierr = m_field.modify_finite_element_add_field_data("SKIN_FE","HDIV"); CHKERRQ(ierr);
  //add entities to finite element
  ierr = m_field.add_ents_to_finite_element_by_TETs(root_set,"TET_FE"); CHKERRQ(ierr);
  Range tets;
  ierr = m_field.get_entities_by_type_and_ref_level(BitRefLevel().set(0),BitRefLevel().set(),MBTET,tets); CHKERRQ(ierr);
  Skinner skin(&moab);
  Range skin_faces; // skin faces from 3d ents
  rval = skin.find_skin(0,tets,false,skin_faces); CHKERR_MOAB(rval);
  ierr = m_field.add_ents_to_finite_element_by_TRIs(skin_faces,"SKIN_FE"); CHKERRQ(ierr);

  //build finite elemnts
  ierr = m_field.build_finite_elements(); CHKERRQ(ierr);

  //build adjacencies
  ierr = m_field.build_adjacencies(bit_level0); CHKERRQ(ierr);

  //problem
  ierr = m_field.add_problem("TEST_PROBLEM"); CHKERRQ(ierr);
  //set finite elements for problem
  ierr = m_field.modify_problem_add_finite_element("TEST_PROBLEM","TET_FE"); CHKERRQ(ierr);
  ierr = m_field.modify_problem_add_finite_element("TEST_PROBLEM","SKIN_FE"); CHKERRQ(ierr);
  //set refinment level for problem
  ierr = m_field.modify_problem_ref_level_add_bit("TEST_PROBLEM",bit_level0); CHKERRQ(ierr);
  //build problem
  ierr = m_field.build_problems(); CHKERRQ(ierr);

  //mesh partitioning
  //partition
  ierr = m_field.partition_simple_problem("TEST_PROBLEM"); CHKERRQ(ierr);
  ierr = m_field.partition_finite_elements("TEST_PROBLEM"); CHKERRQ(ierr);
  //what are ghost nodes, see Petsc Manual
  ierr = m_field.partition_ghost_dofs("TEST_PROBLEM"); CHKERRQ(ierr);

  struct OpTetDivergence: public VolumeElementForcesAndSourcesCore::UserDataOperator {

    FieldData &dIv;
    OpTetDivergence(FieldData &div):
      VolumeElementForcesAndSourcesCore::UserDataOperator("HDIV",UserDataOperator::OPROW),dIv(div) {}

    PetscErrorCode doWork(
      int side,
      EntityType type,
      DataForcesAndSurcesCore::EntData &data) {
      PetscFunctionBegin;

      if(type != MBTRI && type != MBTET) PetscFunctionReturn(0);


      PetscErrorCode ierr;

      if(data.getFieldData().size()==0) PetscFunctionReturn(0);

      //cout << "type " << type << " side " << side << std::endl;

      int nb_gauss_pts = data.getDiffHdivN().size1();
      int nb_dofs = data.getFieldData().size();

      VectorDouble div_vec;
      div_vec.resize(nb_dofs,0);

      int gg = 0;
      for(;gg<nb_gauss_pts;gg++) {
        ierr = getDivergenceMatrixOperator_Hdiv(side,type,data,gg,div_vec); CHKERRQ(ierr);
        //cout << std::fixed << div_vec << std::endl;
        unsigned int dd = 0;
        for(;dd<div_vec.size();dd++) {
          double w = getGaussPts()(3,gg)*getVolume();
          if(getHoGaussPtsDetJac().size()>0) {
            w *= getHoGaussPtsDetJac()[gg]; ///< higher order geometry
          }
          dIv += div_vec[dd]*w;
        }
        //cout << std::fixed << data.getDiffHdivN(gg) << std::endl;
      }

      //cout << std::fixed << data.getDiffHdivN() << std::endl;
      //cout << std::endl;


      PetscFunctionReturn(0);
    }

  };

  struct MyFE: public VolumeElementForcesAndSourcesCore {

    MyFE(FieldInterface &m_field): VolumeElementForcesAndSourcesCore(m_field) {}
    int getRule(int order) { return order; }; //order/2; };

  };

  struct MyTriFE: public FaceElementForcesAndSourcesCore {

    MyTriFE(FieldInterface &m_field): FaceElementForcesAndSourcesCore(m_field) {}
    int getRule(int order) { return order; };//2*order; }; //order/2; };

  };

  struct OpFacesFluxes: public FaceElementForcesAndSourcesCore::UserDataOperator {

    double &dIv;
    OpFacesFluxes(double &div):
      FaceElementForcesAndSourcesCore::UserDataOperator("HDIV",UserDataOperator::OPROW),
      dIv(div) {}

    PetscErrorCode doWork(
      int side,
      EntityType type,
      DataForcesAndSurcesCore::EntData &data
    ) {
      PetscFunctionBegin;

      if(type != MBTRI) PetscFunctionReturn(0);

      int nb_gauss_pts = data.getHdivN().size1();
      int nb_dofs = data.getFieldData().size();

      int gg = 0;
      for(;gg<nb_gauss_pts;gg++) {
        int dd = 0;
        for(;dd<nb_dofs;dd++) {
          double area;
          VectorDouble n;
          if(getNormals_at_GaussPt().size1() == (unsigned int)nb_gauss_pts) {
            n = getNormals_at_GaussPt(gg);
            area = norm_2(getNormals_at_GaussPt(gg))*0.5;
          } else {
            n = getNormal();
            area = getArea();
          }
          n /= norm_2(n);
          dIv +=
          ( n[0]*data.getHdivN(gg)(dd,0) +
          n[1]*data.getHdivN(gg)(dd,1) +
          n[2]*data.getHdivN(gg)(dd,2) )
          *getGaussPts()(2,gg)*area;
        }
        //cout << getNormal() << std::endl;
      }

      PetscFunctionReturn(0);
    }

  };

  double divergence_vol = 0;
  double divergence_skin = 0;


  MyFE tet_fe(m_field);
  tet_fe.getOpPtrVector().push_back(new OpTetDivergence(divergence_vol));

  MyTriFE skin_fe(m_field);
  skin_fe.getOpPtrVector().push_back(new OpFacesFluxes(divergence_skin));

  ierr = m_field.loop_finite_elements("TEST_PROBLEM","TET_FE",tet_fe);  CHKERRQ(ierr);
  ierr = m_field.loop_finite_elements("TEST_PROBLEM","SKIN_FE",skin_fe);  CHKERRQ(ierr);

  std::cout.precision(12);

  std::cout << "divergence_vol " << divergence_vol << std::endl;
  std::cout << "divergence_skin " << divergence_skin << std::endl;

  const double eps = 1e-6;
  if(fabs(divergence_vol-1.-1./3)>eps) {
    SETERRQ2(
      PETSC_COMM_SELF,
      MOFEM_ATOM_TEST_INVALID,
      "invalid divergence_vol = %6.4e, should be %6.4e\n",
      divergence_vol,1+1./3.
    );
  }
  if(fabs(divergence_skin-1.-1./3)>eps) {
     SETERRQ2(PETSC_COMM_SELF,MOFEM_ATOM_TEST_INVALID,"invalid fluxes = %6.4e, should be %6.4e\n",
	divergence_skin,1+1./3.);
  }

  ierr = m_field.add_field("MESH_NODE_POSITIONS",H1,3); CHKERRQ(ierr);

  ierr = m_field.add_ents_to_field_by_TETs(0,"MESH_NODE_POSITIONS"); CHKERRQ(ierr);
  ierr = m_field.set_field_order(0,MBVERTEX,"MESH_NODE_POSITIONS",1); CHKERRQ(ierr);
  ierr = m_field.set_field_order(0,MBEDGE,"MESH_NODE_POSITIONS",1); CHKERRQ(ierr);
  ierr = m_field.set_field_order(0,MBTRI,"MESH_NODE_POSITIONS",1); CHKERRQ(ierr);
  ierr = m_field.set_field_order(0,MBTET,"MESH_NODE_POSITIONS",1); CHKERRQ(ierr);

  ierr = m_field.modify_finite_element_add_field_data("TET_FE","MESH_NODE_POSITIONS"); CHKERRQ(ierr);
  ierr = m_field.modify_finite_element_add_field_data("SKIN_FE","MESH_NODE_POSITIONS"); CHKERRQ(ierr);

  ierr = m_field.build_fields(); CHKERRQ(ierr);
  ierr = m_field.build_finite_elements(); CHKERRQ(ierr);
  ierr = m_field.build_adjacencies(bit_level0); CHKERRQ(ierr);
  ierr = m_field.build_problems(); CHKERRQ(ierr);

  //mesh partitioning
  ierr = m_field.partition_simple_problem("TEST_PROBLEM"); CHKERRQ(ierr);
  ierr = m_field.partition_finite_elements("TEST_PROBLEM"); CHKERRQ(ierr);
  ierr = m_field.partition_ghost_dofs("TEST_PROBLEM"); CHKERRQ(ierr);

  for(_IT_GET_DOFS_FIELD_BY_NAME_AND_TYPE_FOR_LOOP_(m_field,"MESH_NODE_POSITIONS",MBVERTEX,dof)) {
    EntityHandle vert = (*dof)->getEnt();
    double coords[3];
    rval = moab.get_coords(&vert,1,coords); CHKERR_MOAB(rval);
    coords[0] *= 2;
    coords[1] *= 4;
    coords[2] *= 0.5;

    (*dof)->getFieldData() = coords[(*dof)->getDofCoeffIdx()];
  }

  divergence_vol = 0;
  divergence_skin = 0;

  ierr = m_field.loop_finite_elements("TEST_PROBLEM","TET_FE",tet_fe);  CHKERRQ(ierr);
  ierr = m_field.loop_finite_elements("TEST_PROBLEM","SKIN_FE",skin_fe);  CHKERRQ(ierr);

  std::cout.precision(12);

  std::cout << "divergence_vol " << divergence_vol << std::endl;
  std::cout << "divergence_skin " << divergence_skin << std::endl;

  if(fabs(divergence_skin-divergence_vol)>eps) {
     SETERRQ2(PETSC_COMM_SELF,MOFEM_ATOM_TEST_INVALID,"invalid surface flux or divergence or both\n",
	divergence_skin,divergence_vol);
  }


  ierr = PetscFinalize(); CHKERRQ(ierr);
}
