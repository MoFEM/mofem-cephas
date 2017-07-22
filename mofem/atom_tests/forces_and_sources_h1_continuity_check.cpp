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

static const double face_coords[4][9] = {
  { 0,0,0, 1,0,0, 0,0,1 },
  { 1,0,0, 0,1,0, 0,0,1 },
  { 0,0,0, 0,1,0, 0,0,1 },
  { 0,0,0, 1,0,0, 0,1,0 }
};

static const double edge_coords[6][6] = {
  {0,0,0, 1,0,0},
  {1,0,0, 0,1,0},
  {0,1,0, 0,0,0},
  {0,0,0, 0,0,1},
  {1,0,0, 0,0,1},
  {0,1,0, 0,0,1}
};

int main(int argc, char *argv[]) {

  
  

  PetscInitialize(&argc,&argv,(char *)0,help);

  try {

  moab::Core mb_instance;
  moab::Interface& moab = mb_instance;
  int rank;
  MPI_Comm_rank(PETSC_COMM_WORLD,&rank);

  PetscBool flg = PETSC_TRUE;
  char mesh_file_name[255];
  #if PETSC_VERSION_GE(3,6,4)
  ierr = PetscOptionsGetString(PETSC_NULL,"","-my_file",mesh_file_name,255,&flg); CHKERRQ(ierr);
  #else
  ierr = PetscOptionsGetString(PETSC_NULL,PETSC_NULL,"-my_file",mesh_file_name,255,&flg); CHKERRQ(ierr);
  #endif
  if(flg != PETSC_TRUE) {
    SETERRQ(PETSC_COMM_SELF,1,"*** ERROR -my_file (MESH FILE NEEDED)");
  }
  const char *option;
  option = "";
  rval = moab.load_file(mesh_file_name, 0, option); CHKERRQ_MOAB(rval);

  ParallelComm* pcomm = ParallelComm::get_pcomm(&moab,MYPCOMM_INDEX);
  if(pcomm == NULL) pcomm =  new ParallelComm(&moab,PETSC_COMM_WORLD);

  //Create MoFEM (Joseph) database
  MoFEM::Core core(moab);
  MoFEM::Interface& m_field = core;

  //set entitities bit level
  BitRefLevel bit_level0;
  bit_level0.set(0);
  EntityHandle meshset_level0;
  rval = moab.create_meshset(MESHSET_SET,meshset_level0); CHKERRQ_MOAB(rval);
  ierr = m_field.seed_ref_level_3D(0,bit_level0); CHKERRQ(ierr);

  //Fields
  ierr = m_field.add_field("MESH_NODE_POSITIONS",H1,AINSWORTH_LEGENDRE_BASE,3); CHKERRQ(ierr);
  ierr = m_field.add_field("H1",H1,AINSWORTH_LEGENDRE_BASE,1); CHKERRQ(ierr);

  //FE
  ierr = m_field.add_finite_element("TET_FE"); CHKERRQ(ierr);
  ierr = m_field.add_finite_element("TRI_FE"); CHKERRQ(ierr);
  ierr = m_field.add_finite_element("SKIN_FE"); CHKERRQ(ierr);
  ierr = m_field.add_finite_element("EDGE_FE"); CHKERRQ(ierr);

  //Define rows/cols and element data
  ierr = m_field.modify_finite_element_add_field_row("TET_FE","H1"); CHKERRQ(ierr);
  ierr = m_field.modify_finite_element_add_field_col("TET_FE","H1"); CHKERRQ(ierr);
  ierr = m_field.modify_finite_element_add_field_data("TET_FE","H1"); CHKERRQ(ierr);
  ierr = m_field.modify_finite_element_add_field_data("TET_FE","MESH_NODE_POSITIONS"); CHKERRQ(ierr);

  ierr = m_field.modify_finite_element_add_field_row("SKIN_FE","H1"); CHKERRQ(ierr);
  ierr = m_field.modify_finite_element_add_field_col("SKIN_FE","H1"); CHKERRQ(ierr);
  ierr = m_field.modify_finite_element_add_field_data("SKIN_FE","H1"); CHKERRQ(ierr);
  ierr = m_field.modify_finite_element_add_field_data("SKIN_FE","MESH_NODE_POSITIONS"); CHKERRQ(ierr);

  ierr = m_field.modify_finite_element_add_field_row("TRI_FE","H1"); CHKERRQ(ierr);
  ierr = m_field.modify_finite_element_add_field_col("TRI_FE","H1"); CHKERRQ(ierr);
  ierr = m_field.modify_finite_element_add_field_data("TRI_FE","H1"); CHKERRQ(ierr);
  ierr = m_field.modify_finite_element_add_field_data("TRI_FE","MESH_NODE_POSITIONS"); CHKERRQ(ierr);

  ierr = m_field.modify_finite_element_add_field_row("EDGE_FE","H1"); CHKERRQ(ierr);
  ierr = m_field.modify_finite_element_add_field_col("EDGE_FE","H1"); CHKERRQ(ierr);
  ierr = m_field.modify_finite_element_add_field_data("EDGE_FE","H1"); CHKERRQ(ierr);
  ierr = m_field.modify_finite_element_add_field_data("EDGE_FE","MESH_NODE_POSITIONS"); CHKERRQ(ierr);

  //Problem
  ierr = m_field.add_problem("TEST_PROBLEM"); CHKERRQ(ierr);

  //set finite elements for problem
  ierr = m_field.modify_problem_add_finite_element("TEST_PROBLEM","TET_FE"); CHKERRQ(ierr);
  ierr = m_field.modify_problem_add_finite_element("TEST_PROBLEM","SKIN_FE"); CHKERRQ(ierr);
  ierr = m_field.modify_problem_add_finite_element("TEST_PROBLEM","TRI_FE"); CHKERRQ(ierr);
  ierr = m_field.modify_problem_add_finite_element("TEST_PROBLEM","EDGE_FE"); CHKERRQ(ierr);

  //set refinement level for problem
  ierr = m_field.modify_problem_ref_level_add_bit("TEST_PROBLEM",bit_level0); CHKERRQ(ierr);

  //meshset consisting all entities in mesh
  EntityHandle root_set = moab.get_root_set();
  //add entities to field
  ierr = m_field.add_ents_to_field_by_type(root_set,MBTET,"H1"); CHKERRQ(ierr);

  //add entities to finite element
  ierr = m_field.add_ents_to_finite_element_by_type(root_set,MBTET,"TET_FE"); CHKERRQ(ierr);

  Range tets;
  ierr = m_field.get_entities_by_type_and_ref_level(BitRefLevel().set(0),BitRefLevel().set(),MBTET,tets); CHKERRQ(ierr);
  Skinner skin(&moab);
  Range skin_faces; // skin faces from 3d ents
  rval = skin.find_skin(0,tets,false,skin_faces); CHKERRQ_MOAB(rval);
  ierr = m_field.add_ents_to_finite_element_by_type(skin_faces,MBTRI,"SKIN_FE"); CHKERRQ(ierr);

  Range faces;
  ierr = m_field.get_entities_by_type_and_ref_level(BitRefLevel().set(0),BitRefLevel().set(),MBTRI,faces); CHKERRQ(ierr);
  faces = subtract(faces,skin_faces);
  ierr = m_field.add_ents_to_finite_element_by_type(faces,MBTRI,"TRI_FE"); CHKERRQ(ierr);

  Range edges;
  ierr = moab.get_adjacencies(faces,1,false,edges,moab::Interface::UNION); CHKERRQ(ierr);
  ierr = m_field.add_ents_to_finite_element_by_type(edges,MBEDGE,"EDGE_FE"); CHKERRQ(ierr);
  Range skin_edges;
  ierr = moab.get_adjacencies(skin_faces,1,false,skin_edges,moab::Interface::UNION); CHKERRQ(ierr);
  ierr = m_field.add_ents_to_finite_element_by_type(edges,MBEDGE,"EDGE_FE"); CHKERRQ(ierr);

  //set app. order
  int order = 4;
  ierr = m_field.set_field_order(root_set,MBTET,"H1",order); CHKERRQ(ierr);
  ierr = m_field.set_field_order(root_set,MBTRI,"H1",order); CHKERRQ(ierr);
  ierr = m_field.set_field_order(root_set,MBEDGE,"H1",order); CHKERRQ(ierr);


  ierr = m_field.add_ents_to_field_by_type(0,MBTET,"MESH_NODE_POSITIONS"); CHKERRQ(ierr);
  ierr = m_field.set_field_order(0,MBVERTEX,"MESH_NODE_POSITIONS",1); CHKERRQ(ierr);
  ierr = m_field.set_field_order(0,MBEDGE,"MESH_NODE_POSITIONS",2); CHKERRQ(ierr);
  ierr = m_field.set_field_order(0,MBTRI,"MESH_NODE_POSITIONS",2); CHKERRQ(ierr);
  ierr = m_field.set_field_order(0,MBTET,"MESH_NODE_POSITIONS",2); CHKERRQ(ierr);

  /****/
  //build database
  //build field
  ierr = m_field.build_fields(); CHKERRQ(ierr);
  //build finite elemnts
  ierr = m_field.build_finite_elements(); CHKERRQ(ierr);
  //build adjacencies
  ierr = m_field.build_adjacencies(bit_level0); CHKERRQ(ierr);


  ProblemsManager *prb_mng_ptr;
  ierr = m_field.query_interface(prb_mng_ptr); CHKERRQ(ierr);
  //build problem
  ierr = prb_mng_ptr->buildProblem("TEST_PROBLEM",true); CHKERRQ(ierr);

  // project geometry form 10 node tets on higher order approx. functions
  Projection10NodeCoordsOnField ent_method(m_field,"MESH_NODE_POSITIONS");
  ierr = m_field.loop_dofs("MESH_NODE_POSITIONS",ent_method); CHKERRQ(ierr);

  /****/
  //mesh partitioning
  //partition
  ierr = prb_mng_ptr->partitionSimpleProblem("TEST_PROBLEM"); CHKERRQ(ierr);
  ierr = prb_mng_ptr->partitionFiniteElements("TEST_PROBLEM"); CHKERRQ(ierr);
  //what are ghost nodes, see Petsc Manual
  ierr = prb_mng_ptr->partitionGhostDofs("TEST_PROBLEM"); CHKERRQ(ierr);

  Vec v;
  ierr = m_field.VecCreateGhost("TEST_PROBLEM",ROW,&v);
  ierr = VecSetRandom(v,PETSC_NULL); CHKERRQ(ierr);
  ierr = m_field.set_local_ghost_vector("TEST_PROBLEM",ROW,v,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
  ierr = VecDestroy(&v); CHKERRQ(ierr);


  typedef tee_device<std::ostream, std::ofstream> TeeDevice;
  typedef stream<TeeDevice> TeeStream;
  std::ofstream ofs("forces_and_sources_hdiv_continuity_check.txt");
  TeeDevice my_tee(std::cout, ofs);
  TeeStream my_split(my_tee);

  struct OpTetFluxes: public VolumeElementForcesAndSourcesCore::UserDataOperator {

    MoFEM::Interface &mField;
    Tag tH;

    OpTetFluxes(MoFEM::Interface &m_field,Tag th):
      VolumeElementForcesAndSourcesCore::UserDataOperator("H1",UserDataOperator::OPROW),
      mField(m_field),tH(th) {}

    PetscErrorCode doWork(
      int side,
      EntityType type,
      DataForcesAndSurcesCore::EntData &data) {
      PetscFunctionBegin;

      

      if(data.getFieldData().size()==0) PetscFunctionReturn(0);

      if(type == MBTRI) {

        boost::shared_ptr<const NumeredEntFiniteElement> mofem_fe = getNumeredEntFiniteElementPtr();
        SideNumber_multiIndex &side_table = mofem_fe->getSideNumberTable();
        EntityHandle face = side_table.get<1>().find(boost::make_tuple(type,side))->get()->ent;
        int sense = side_table.get<1>().find(boost::make_tuple(type,side))->get()->sense;

        double t = 0;
        int nb_dofs = data.getN().size2();
        for(int dd = 0;dd<nb_dofs;dd++) {
            t += data.getN(side)(dd)*data.getFieldData()[dd];
        }

        double *t_ptr;
        rval = mField.get_moab().tag_get_by_ptr(tH,&face,1,(const void **)&t_ptr); CHKERRQ_MOAB(rval);
        *t_ptr += sense*t;

      }

      if(type == MBEDGE) {

        boost::shared_ptr<const NumeredEntFiniteElement> mofem_fe = getNumeredEntFiniteElementPtr();
        SideNumber_multiIndex &side_table = mofem_fe->getSideNumberTable();
        EntityHandle edge = side_table.get<1>().find(boost::make_tuple(type,side))->get()->ent;
        Range adj_tets;
        rval = mField.get_moab().get_adjacencies(&edge,1,3,false,adj_tets,moab::Interface::UNION);
        const int nb_adj_tets = adj_tets.size();

        double t = 0;
        int nb_dofs = data.getN().size2();
        for(int dd = 0;dd<nb_dofs;dd++) {
          t += data.getN(4+side)(dd)*data.getFieldData()[dd];
        }

        double *t_ptr;
        rval = mField.get_moab().tag_get_by_ptr(tH,&edge,1,(const void **)&t_ptr); CHKERRQ_MOAB(rval);
        *t_ptr += t/nb_adj_tets;

      }

      PetscFunctionReturn(0);
    }

  };

  struct MyTetFE: public VolumeElementForcesAndSourcesCore {

    MyTetFE(MoFEM::Interface &m_field):
    VolumeElementForcesAndSourcesCore(m_field) {}
    int getRule(int order) { return -1; };

    MatrixDouble N_tri;
    PetscErrorCode setGaussPts(int order) {
      
      PetscFunctionBegin;

      try {

        N_tri.resize(1,3);
        ierr = ShapeMBTRI(&N_tri(0,0),G_TRI_X1,G_TRI_Y1,1); CHKERRQ(ierr);

        gaussPts.resize(4,4+6);
        int ff = 0;
        for(;ff<4;ff++) {
          int dd = 0;
          for(;dd<3;dd++) {
            gaussPts(dd,ff) = cblas_ddot(3,&N_tri(0,0),1,&face_coords[ff][dd],3);
          }
          gaussPts(3,ff) = G_TRI_W1[0];
        }

        int ee = 0;
        for(;ee<6;ee++) {
          int dd = 0;
          for(;dd<3;dd++) {
            gaussPts(dd,4+ee) = (edge_coords[ee][0+dd]+edge_coords[ee][3+dd])*0.5;
          }
          gaussPts(3,4+ee) = 1;
        }

        //std::cerr << gaussPts << std::endl;

      } catch (std::exception& ex) {
        std::ostringstream ss;
        ss << "thorw in method: " << ex.what() << " at line " << __LINE__ << " in file " << __FILE__;
        SETERRQ(PETSC_COMM_SELF,MOFEM_STD_EXCEPTION_THROW,ss.str().c_str());
      }

      PetscFunctionReturn(0);
    }


  };

  struct OpFacesSkinFluxes: public FaceElementForcesAndSourcesCore::UserDataOperator {

    MoFEM::Interface &mField;
    Tag tH1,tH2;
    TeeStream &mySplit;

    OpFacesSkinFluxes(MoFEM::Interface &m_field,Tag th1,Tag th2,TeeStream &my_split):
      FaceElementForcesAndSourcesCore::UserDataOperator("H1",UserDataOperator::OPROW),
      mField(m_field),
      tH1(th1),
      tH2(th2),
      mySplit(my_split) {}

    PetscErrorCode doWork(
      int side,
      EntityType type,
      DataForcesAndSurcesCore::EntData &data) {
      PetscFunctionBegin;

      
      if(type != MBTRI) PetscFunctionReturn(0);
      EntityHandle face = getNumeredEntFiniteElementPtr()->getEnt();

      double *t_ptr;
      rval = mField.get_moab().tag_get_by_ptr(tH1,&face,1,(const void **)&t_ptr); CHKERRQ_MOAB(rval);
      double *tn_ptr;
      rval = mField.get_moab().tag_get_by_ptr(tH2,&face,1,(const void **)&tn_ptr); CHKERRQ_MOAB(rval);

      *tn_ptr = *t_ptr;

      int dd = 0;
      int nb_dofs = data.getN().size2();
      for(;dd<nb_dofs;dd++) {
        double val = data.getFieldData()[dd];
        *tn_ptr += -data.getN()(0,dd)*val;
      }

      const double eps = 1e-8;
      if(fabs(*tn_ptr)>eps) {
        SETERRQ1(
          PETSC_COMM_SELF,
          MOFEM_ATOM_TEST_INVALID,
          "H1 continuity failed %6.4e",
          *tn_ptr
        );
      }

      mySplit.precision(5);
      mySplit << face << " " << /*std::fixed <<*/ fabs(*tn_ptr) << std::endl;

      PetscFunctionReturn(0);
    }

  };

  struct OpFacesFluxes: public FaceElementForcesAndSourcesCore::UserDataOperator {

    MoFEM::Interface &mField;
    Tag tH1,tH2;
    TeeStream &mySplit;

    OpFacesFluxes(MoFEM::Interface &m_field,Tag _th1,Tag _th2,TeeStream &my_split):
      FaceElementForcesAndSourcesCore::UserDataOperator("H1",UserDataOperator::OPROW),
      mField(m_field),tH1(_th1),tH2(_th2),mySplit(my_split) {}

    PetscErrorCode doWork(
      int side,
      EntityType type,
      DataForcesAndSurcesCore::EntData &data) {
      PetscFunctionBegin;

      

      if(type != MBTRI) PetscFunctionReturn(0);

      EntityHandle face = getNumeredEntFiniteElementPtr()->getEnt();

      double *t_ptr;
      rval = mField.get_moab().tag_get_by_ptr(tH1,&face,1,(const void **)&t_ptr); CHKERRQ_MOAB(rval);
      double *tn_ptr;
      rval = mField.get_moab().tag_get_by_ptr(tH2,&face,1,(const void **)&tn_ptr); CHKERRQ_MOAB(rval);

      *tn_ptr = *t_ptr;

      const double eps = 1e-8;
      if(fabs(*tn_ptr)>eps) {
        SETERRQ1(
          PETSC_COMM_SELF,
          MOFEM_ATOM_TEST_INVALID,
          "H1 continuity failed %6.4e",
          *tn_ptr
        );
      }

      mySplit.precision(5);

      mySplit << face << " " << /*std::fixed <<*/ fabs(*tn_ptr) << std::endl;

      PetscFunctionReturn(0);
    }

  };

  struct MyTriFE: public FaceElementForcesAndSourcesCore {

    MyTriFE(MoFEM::Interface &m_field): FaceElementForcesAndSourcesCore(m_field) {}
    int getRule(int order) { return -1; };

    PetscErrorCode setGaussPts(int order) {
      PetscFunctionBegin;

      gaussPts.resize(3,1);
      gaussPts(0,0) = G_TRI_X1[0];
      gaussPts(1,0) = G_TRI_Y1[0];
      gaussPts(2,0) = G_TRI_W1[0];

      PetscFunctionReturn(0);
    }

  };

  struct OpEdgesFluxes: public EdgeElementForcesAndSurcesCore::UserDataOperator {

    MoFEM::Interface &mField;
    Tag tH1,tH2;
    TeeStream &mySplit;

    OpEdgesFluxes(MoFEM::Interface &m_field,Tag _th1,Tag _th2,TeeStream &my_split):
      EdgeElementForcesAndSurcesCore::UserDataOperator("H1",UserDataOperator::OPROW),
      mField(m_field),tH1(_th1),tH2(_th2),mySplit(my_split) {}

    PetscErrorCode doWork(
      int side,
      EntityType type,
      DataForcesAndSurcesCore::EntData &data) {
      PetscFunctionBegin;

      

      if(type != MBEDGE) PetscFunctionReturn(0);

      EntityHandle edge = getNumeredEntFiniteElementPtr()->getEnt();

      double *t_ptr;
      rval = mField.get_moab().tag_get_by_ptr(tH1,&edge,1,(const void **)&t_ptr); CHKERRQ_MOAB(rval);
      double *tn_ptr;
      rval = mField.get_moab().tag_get_by_ptr(tH2,&edge,1,(const void **)&tn_ptr); CHKERRQ_MOAB(rval);

      *tn_ptr = *t_ptr;

      double tn = 0;
      int nb_dofs = data.getN().size2();
      int dd = 0;
      for(;dd<nb_dofs;dd++) {
        double val = data.getFieldData()[dd];
        tn += data.getN()(0,dd)*val;
      }

      // mySplit << *tn_ptr << "  " << tn << " " << getLength() << endl;
      *tn_ptr -= tn;

      // mySplit << getTangetAtGaussPts() << " " << getDirection() << endl;

      // cerr << t_ptr[0] << " " << t_ptr[1] << " " << t_ptr[2] << endl;

      // const double eps = 1e-8;
      // if(fabs(*tn_ptr)>eps) {
      //   SETERRQ1(
      //     PETSC_COMM_SELF,
      //     MOFEM_ATOM_TEST_INVALID,
      //     "H1 continuity failed %6.4e",
      //     *tn_ptr
      //   );
      // }

      mySplit.precision(5);

      mySplit << edge << " " << /*std::fixed <<*/ fabs(*tn_ptr) << std::endl;

      PetscFunctionReturn(0);
    }

  };

  struct MyEdgeFE: public EdgeElementForcesAndSurcesCore {

    MyEdgeFE(MoFEM::Interface &m_field):
    EdgeElementForcesAndSurcesCore(m_field) {}
    int getRule(int order) { return -1; };

    PetscErrorCode setGaussPts(int order) {
      PetscFunctionBegin;

      gaussPts.resize(2,1);
      gaussPts(0,0) = 0.5;
      gaussPts(1,0) = 1;

      PetscFunctionReturn(0);
    }

  };



  MyTetFE tet_fe(m_field);
  MyTriFE tri_fe(m_field);
  MyTriFE skin_fe(m_field);
  MyEdgeFE edge_fe(m_field);

  Tag th1;
  double def_val[] = {0,0,0};
  rval = moab.tag_get_handle("T",3,MB_TYPE_DOUBLE,th1,MB_TAG_CREAT|MB_TAG_SPARSE,&def_val); CHKERRQ_MOAB(rval);
  tet_fe.getOpPtrVector().push_back(new OpTetFluxes(m_field,th1));

  Tag th2;
  rval = moab.tag_get_handle("TN",1,MB_TYPE_DOUBLE,th2,MB_TAG_CREAT|MB_TAG_SPARSE,&def_val); CHKERRQ_MOAB(rval);
  tri_fe.getOpPtrVector().push_back(new OpFacesFluxes(m_field,th1,th2,my_split));
  skin_fe.getOpPtrVector().push_back(new OpFacesSkinFluxes(m_field,th1,th2,my_split));
  edge_fe.getOpPtrVector().push_back(new OpEdgesFluxes(m_field,th1,th2,my_split));

  for(Range::iterator fit = faces.begin();fit!=faces.end();fit++) {
    rval = moab.tag_set_data(th1,&*fit,1,&def_val); CHKERRQ_MOAB(rval);
    rval = moab.tag_set_data(th2,&*fit,1,&def_val); CHKERRQ_MOAB(rval);
  }

  ierr = m_field.loop_finite_elements("TEST_PROBLEM","TET_FE",tet_fe);  CHKERRQ(ierr);
  my_split << "internal\n";
  ierr = m_field.loop_finite_elements("TEST_PROBLEM","TRI_FE",tri_fe);  CHKERRQ(ierr);
  my_split << "skin\n";
  ierr = m_field.loop_finite_elements("TEST_PROBLEM","SKIN_FE",skin_fe);  CHKERRQ(ierr);
  my_split << "edges\n";
  ierr = m_field.loop_finite_elements("TEST_PROBLEM","EDGE_FE",edge_fe);  CHKERRQ(ierr);


  EntityHandle meshset;
  rval = moab.create_meshset(MESHSET_SET,meshset); CHKERRQ_MOAB(rval);
  ierr = m_field.get_entities_by_type_and_ref_level(BitRefLevel().set(0),BitRefLevel().set(),MBTRI,meshset); CHKERRQ(ierr);
  rval = moab.write_file("out.vtk","VTK","",&meshset,1); CHKERRQ_MOAB(rval);


  } catch (MoFEMException const &e) {
    SETERRQ(PETSC_COMM_SELF,e.errorCode,e.errorMessage);
  }

  ierr = PetscFinalize(); CHKERRQ(ierr);
}
