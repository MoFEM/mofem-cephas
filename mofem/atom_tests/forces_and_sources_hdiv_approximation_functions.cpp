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

namespace bio = boost::iostreams;
using bio::tee_device;
using bio::stream;

using namespace MoFEM;

static char help[] = "...\n\n";

int main(int argc, char *argv[]) {

  
  

  PetscInitialize(&argc,&argv,(char *)0,help);

  try {

  moab::Core mb_instance;
  moab::Interface& moab = mb_instance;
  int rank;
  MPI_Comm_rank(PETSC_COMM_WORLD,&rank);

  // create one tet
  double tet_coords[] = {
    0,0,0,
    2,0,0,
    0,2,0,
    0,0,2
  };
  EntityHandle nodes[4];
  for(int nn = 0;nn<4;nn++) {
    rval = moab.create_vertex(&tet_coords[3*nn],nodes[nn]); CHKERRG(rval);
  }
  EntityHandle tet;
  rval = moab.create_element(MBTET,nodes,4,tet); CHKERRG(rval);

  ParallelComm* pcomm = ParallelComm::get_pcomm(&moab,MYPCOMM_INDEX);
  if(pcomm == NULL) pcomm =  new ParallelComm(&moab,PETSC_COMM_WORLD);

  //Create MoFEM (Joseph) database
  MoFEM::Core core(moab);
  MoFEM::Interface& m_field = core;

  //set entitities bit level
  BitRefLevel bit_level0;
  bit_level0.set(0);
  EntityHandle meshset_level0;
  rval = moab.create_meshset(MESHSET_SET,meshset_level0); CHKERRG(rval);
  ierr = m_field.getInterface<BitRefManager>()->setBitRefLevelByDim(0,3,bit_level0); CHKERRG(ierr);

  //Fields
  ierr = m_field.add_field("HDIV",HDIV,AINSWORTH_LEGENDRE_BASE,1); CHKERRG(ierr);

  //FE TET
  ierr = m_field.add_finite_element("HDIV_TET_FE"); CHKERRG(ierr);
  //Define rows/cols and element data
  ierr = m_field.modify_finite_element_add_field_row("HDIV_TET_FE","HDIV"); CHKERRG(ierr);
  ierr = m_field.modify_finite_element_add_field_col("HDIV_TET_FE","HDIV"); CHKERRG(ierr);
  ierr = m_field.modify_finite_element_add_field_data("HDIV_TET_FE","HDIV"); CHKERRG(ierr);

  //FE TRI
  ierr = m_field.add_finite_element("HDIV_TRI_FE"); CHKERRG(ierr);
  //Define rows/cols and element data
  ierr = m_field.modify_finite_element_add_field_row("HDIV_TRI_FE","HDIV"); CHKERRG(ierr);
  ierr = m_field.modify_finite_element_add_field_col("HDIV_TRI_FE","HDIV"); CHKERRG(ierr);
  ierr = m_field.modify_finite_element_add_field_data("HDIV_TRI_FE","HDIV"); CHKERRG(ierr);

  //Problem
  ierr = m_field.add_problem("TEST_PROBLEM"); CHKERRG(ierr);

  //set finite elements for problem
  ierr = m_field.modify_problem_add_finite_element("TEST_PROBLEM","HDIV_TET_FE"); CHKERRG(ierr);
  ierr = m_field.modify_problem_add_finite_element("TEST_PROBLEM","HDIV_TRI_FE"); CHKERRG(ierr);
  //set refinement level for problem
  ierr = m_field.modify_problem_ref_level_add_bit("TEST_PROBLEM",bit_level0); CHKERRG(ierr);


  //meshset consisting all entities in mesh
  EntityHandle root_set = moab.get_root_set();
  //add entities to field
  ierr = m_field.add_ents_to_field_by_type(root_set,MBTET,"HDIV"); CHKERRG(ierr);
  //add entities to finite element
  ierr = m_field.add_ents_to_finite_element_by_type(root_set,MBTET,"HDIV_TET_FE"); CHKERRG(ierr);

  Range tets;
  ierr = m_field.getInterface<BitRefManager>()->getEntitiesByTypeAndRefLevel(
    BitRefLevel().set(0),BitRefLevel().set(),MBTET,tets
  ); CHKERRG(ierr);
  Skinner skin(&moab);
  Range skin_faces; // skin faces from 3d ents
  rval = skin.find_skin(0,tets,false,skin_faces); CHKERRG(rval);
  ierr = m_field.add_ents_to_finite_element_by_type(skin_faces,MBTRI,"HDIV_TRI_FE"); CHKERRG(ierr);

  //set app. order
  int order = 4;
  ierr = m_field.set_field_order(root_set,MBTET,"HDIV",order); CHKERRG(ierr);
  ierr = m_field.set_field_order(root_set,MBTRI,"HDIV",order); CHKERRG(ierr);

  /****/
  //build database
  //build field
  ierr = m_field.build_fields(); CHKERRG(ierr);
  //build finite elemnts
  ierr = m_field.build_finite_elements(); CHKERRG(ierr);
  //build adjacencies
  ierr = m_field.build_adjacencies(bit_level0); CHKERRG(ierr);
  //build problem

  ProblemsManager *prb_mng_ptr;
  ierr = m_field.getInterface(prb_mng_ptr); CHKERRG(ierr);
  ierr = prb_mng_ptr->buildProblem("TEST_PROBLEM",true); CHKERRG(ierr);

  /****/
  //mesh partitioning
  //partition
  ierr = prb_mng_ptr->partitionSimpleProblem("TEST_PROBLEM"); CHKERRG(ierr);
  ierr = prb_mng_ptr->partitionFiniteElements("TEST_PROBLEM"); CHKERRG(ierr);
  //what are ghost nodes, see Petsc Manual
  ierr = prb_mng_ptr->partitionGhostDofs("TEST_PROBLEM"); CHKERRG(ierr);

  typedef tee_device<std::ostream, std::ofstream> TeeDevice;
  typedef stream<TeeDevice> TeeStream;

  std::ofstream ofs("forces_and_sources_hdiv_approximation_functions.txt");
  TeeDevice my_tee(std::cout, ofs);
  TeeStream my_split(my_tee);

  struct OpPrintingHdivApproximationFunctions: public VolumeElementForcesAndSourcesCore::UserDataOperator {

    TeeStream &mySplit;
    OpPrintingHdivApproximationFunctions(TeeStream &my_split):
    VolumeElementForcesAndSourcesCore::UserDataOperator("HDIV",UserDataOperator::OPROW),
    mySplit(my_split) {
    }

    MoFEMErrorCode doWork(
      int side,
      EntityType type,
      DataForcesAndSourcesCore::EntData &data
    ) {
      MoFEMFunctionBeginHot;

      if(data.getFieldData().size()==0) MoFEMFunctionReturnHot(0);

      mySplit << std::endl << "type " << type << " side " << side << std::endl;
      mySplit.precision(5);

      const double eps = 1e-6;
      for(unsigned int dd = 0;dd<data.getHdivN().data().size();dd++) {
        if(fabs(data.getHdivN().data()[dd])<eps) data.getHdivN().data()[dd] = 0;
      }
      for(unsigned int dd = 0;dd<data.getDiffHdivN().data().size();dd++) {
        if(fabs(data.getDiffHdivN().data()[dd])<eps) data.getDiffHdivN().data()[dd] = 0;
      }

      mySplit << std::fixed << data.getHdivN() << std::endl;
      mySplit << std::fixed << data.getDiffHdivN() << std::endl;

      MoFEMFunctionReturnHot(0);
    }

  };

  struct MyFE: public VolumeElementForcesAndSourcesCore {

    MyFE(MoFEM::Interface &m_field): VolumeElementForcesAndSourcesCore(m_field) {}
    int getRule(int order) { return 1; };

  };

  struct OpFacePrintingHdivApproximationFunctions: public FaceElementForcesAndSourcesCore::UserDataOperator {

    TeeStream &mySplit;
    OpFacePrintingHdivApproximationFunctions(TeeStream &my_split):
    FaceElementForcesAndSourcesCore::UserDataOperator("HDIV",UserDataOperator::OPROW),
    mySplit(my_split) {
    }

    MoFEMErrorCode doWork(
      int side,
      EntityType type,
      DataForcesAndSourcesCore::EntData &data
    ) {
      MoFEMFunctionBeginHot;

      if(data.getFieldData().size()==0) MoFEMFunctionReturnHot(0);

      mySplit << std::endl << "type " << type << " side " << side << std::endl;
      mySplit.precision(5);

      const double eps = 1e-6;
      for(unsigned int dd = 0;dd<data.getHdivN().data().size();dd++) {
        if(fabs(data.getHdivN().data()[dd])<eps) data.getHdivN().data()[dd] = 0;
      }
      for(unsigned int dd = 0;dd<data.getDiffHdivN().data().size();dd++) {
        if(fabs(data.getDiffHdivN().data()[dd])<eps) data.getDiffHdivN().data()[dd] = 0;
      }

      mySplit << std::fixed << data.getHdivN() << std::endl;
      // mySplit << std::fixed << data.getDiffHdivN() << std::endl;

      MoFEMFunctionReturnHot(0);
    }

  };

  struct MyFaceFE: public FaceElementForcesAndSourcesCore {

    MyFaceFE(MoFEM::Interface &m_field):
    FaceElementForcesAndSourcesCore(m_field) {}
    int getRule(int order) { return 1; };

  };


  MyFE tet_fe(m_field);
  MyFaceFE tri_fe(m_field);

  tet_fe.getOpPtrVector().push_back(new OpPrintingHdivApproximationFunctions(my_split));
  tri_fe.getOpPtrVector().push_back(new OpFacePrintingHdivApproximationFunctions(my_split));

  ierr = m_field.loop_finite_elements("TEST_PROBLEM","HDIV_TET_FE",tet_fe);  CHKERRG(ierr);
  ierr = m_field.loop_finite_elements("TEST_PROBLEM","HDIV_TRI_FE",tri_fe);  CHKERRG(ierr);

  } catch (MoFEMException const &e) {
    SETERRQ(PETSC_COMM_SELF,e.errorCode,e.errorMessage);
  }

  ierr = PetscFinalize(); CHKERRG(ierr);
}
