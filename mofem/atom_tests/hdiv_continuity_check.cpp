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

int main(int argc, char *argv[]) {




  PetscInitialize(&argc,&argv,(char *)0,help);

  try {

  enum bases {
    HDIV_AINSWORTH,
    HDIV_DEMKOWICZ,
    LASTOP
  };

  const char *list[] = {
    "hdiv_ainsworth",
    "hdiv_demkowicz"
  };



  PetscBool flg;
  PetscInt choise_value = HDIV_AINSWORTH;
  ierr = PetscOptionsGetEList(
    PETSC_NULL,NULL,"-base",list,LASTOP,&choise_value,&flg
  ); CHKERRG(ierr);
  if(flg != PETSC_TRUE) {
    SETERRQ(PETSC_COMM_SELF,MOFEM_IMPOSIBLE_CASE,"base not set");
  }

  moab::Core mb_instance;
  moab::Interface& moab = mb_instance;
  int rank;
  MPI_Comm_rank(PETSC_COMM_WORLD,&rank);

  char mesh_file_name[255];
  #if PETSC_VERSION_GE(3,6,4)
  ierr = PetscOptionsGetString(
    PETSC_NULL,"","-my_file",mesh_file_name,255,&flg
  ); CHKERRG(ierr);
  #else
  ierr = PetscOptionsGetString(
    PETSC_NULL,PETSC_NULL,"-my_file",mesh_file_name,255,&flg
  ); CHKERRG(ierr);
  #endif
  if(flg != PETSC_TRUE) {
    SETERRQ(PETSC_COMM_SELF,1,"*** ERROR -my_file (MESH FILE NEEDED)");
  }
  const char *option;
  option = "";
  rval = moab.load_file(mesh_file_name, 0, option); CHKERRG(rval);

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
  ierr = m_field.add_field("MESH_NODE_POSITIONS",H1,AINSWORTH_LEGENDRE_BASE,3); CHKERRG(ierr);
  switch (choise_value) {
    case HDIV_AINSWORTH:
    ierr = m_field.add_field("HDIV",HDIV,AINSWORTH_LEGENDRE_BASE,1); CHKERRG(ierr);
    break;
    case HDIV_DEMKOWICZ:
    ierr = m_field.add_field("HDIV",HDIV,DEMKOWICZ_JACOBI_BASE,1); CHKERRG(ierr);
    break;
  }

  //FE
  ierr = m_field.add_finite_element("TET_FE"); CHKERRG(ierr);
  ierr = m_field.add_finite_element("TRI_FE"); CHKERRG(ierr);
  ierr = m_field.add_finite_element("SKIN_FE"); CHKERRG(ierr);

  //Define rows/cols and element data
  ierr = m_field.modify_finite_element_add_field_row("TET_FE","HDIV"); CHKERRG(ierr);
  ierr = m_field.modify_finite_element_add_field_col("TET_FE","HDIV"); CHKERRG(ierr);
  ierr = m_field.modify_finite_element_add_field_data("TET_FE","HDIV"); CHKERRG(ierr);
  ierr = m_field.modify_finite_element_add_field_data("TET_FE","MESH_NODE_POSITIONS"); CHKERRG(ierr);

  ierr = m_field.modify_finite_element_add_field_row("SKIN_FE","HDIV"); CHKERRG(ierr);
  ierr = m_field.modify_finite_element_add_field_col("SKIN_FE","HDIV"); CHKERRG(ierr);
  ierr = m_field.modify_finite_element_add_field_data("SKIN_FE","HDIV"); CHKERRG(ierr);
  ierr = m_field.modify_finite_element_add_field_data("SKIN_FE","MESH_NODE_POSITIONS"); CHKERRG(ierr);

  ierr = m_field.modify_finite_element_add_field_row("TRI_FE","HDIV"); CHKERRG(ierr);
  ierr = m_field.modify_finite_element_add_field_col("TRI_FE","HDIV"); CHKERRG(ierr);
  ierr = m_field.modify_finite_element_add_field_data("TRI_FE","HDIV"); CHKERRG(ierr);
  ierr = m_field.modify_finite_element_add_field_data("TRI_FE","MESH_NODE_POSITIONS"); CHKERRG(ierr);

  //Problem
  ierr = m_field.add_problem("TEST_PROBLEM"); CHKERRG(ierr);

  //set finite elements for problem
  ierr = m_field.modify_problem_add_finite_element("TEST_PROBLEM","TET_FE"); CHKERRG(ierr);
  ierr = m_field.modify_problem_add_finite_element("TEST_PROBLEM","SKIN_FE"); CHKERRG(ierr);
  ierr = m_field.modify_problem_add_finite_element("TEST_PROBLEM","TRI_FE"); CHKERRG(ierr);

  //set refinement level for problem
  ierr = m_field.modify_problem_ref_level_add_bit("TEST_PROBLEM",bit_level0); CHKERRG(ierr);

  //meshset consisting all entities in mesh
  EntityHandle root_set = moab.get_root_set();
  //add entities to field
  ierr = m_field.add_ents_to_field_by_type(root_set,MBTET,"HDIV"); CHKERRG(ierr);

  //add entities to finite element
  ierr = m_field.add_ents_to_finite_element_by_type(root_set,MBTET,"TET_FE"); CHKERRG(ierr);

  Range tets;
  ierr = m_field.getInterface<BitRefManager>()->getEntitiesByTypeAndRefLevel(BitRefLevel().set(0),BitRefLevel().set(),MBTET,tets); CHKERRG(ierr);
  Skinner skin(&moab);
  Range skin_faces; // skin faces from 3d ents
  rval = skin.find_skin(0,tets,false,skin_faces); CHKERRG(rval);
  ierr = m_field.add_ents_to_finite_element_by_type(skin_faces,MBTRI,"SKIN_FE"); CHKERRG(ierr);

  Range faces;
  ierr = m_field.getInterface<BitRefManager>()->getEntitiesByTypeAndRefLevel(BitRefLevel().set(0),BitRefLevel().set(),MBTRI,faces); CHKERRG(ierr);
  faces = subtract(faces,skin_faces);
  ierr = m_field.add_ents_to_finite_element_by_type(faces,MBTRI,"TRI_FE"); CHKERRG(ierr);

  //set app. order
  int order = 4;
  ierr = m_field.set_field_order(root_set,MBTET,"HDIV",order); CHKERRG(ierr);
  ierr = m_field.set_field_order(root_set,MBTRI,"HDIV",order); CHKERRG(ierr);

  ierr = m_field.add_ents_to_field_by_type(0,MBTET,"MESH_NODE_POSITIONS"); CHKERRG(ierr);
  ierr = m_field.set_field_order(0,MBVERTEX,"MESH_NODE_POSITIONS",1); CHKERRG(ierr);
  ierr = m_field.set_field_order(0,MBEDGE,"MESH_NODE_POSITIONS",2); CHKERRG(ierr);
  ierr = m_field.set_field_order(0,MBTRI,"MESH_NODE_POSITIONS",2); CHKERRG(ierr);
  ierr = m_field.set_field_order(0,MBTET,"MESH_NODE_POSITIONS",2); CHKERRG(ierr);

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

  //project geometry form 10 node tets on higher order approx. functions
  Projection10NodeCoordsOnField ent_method(m_field,"MESH_NODE_POSITIONS");
  ierr = m_field.loop_dofs("MESH_NODE_POSITIONS",ent_method); CHKERRG(ierr);

  /****/
  //mesh partitioning
  //partition
  ierr = prb_mng_ptr->partitionSimpleProblem("TEST_PROBLEM"); CHKERRG(ierr);
  ierr = prb_mng_ptr->partitionFiniteElements("TEST_PROBLEM"); CHKERRG(ierr);
  //what are ghost nodes, see Petsc Manual
  ierr = prb_mng_ptr->partitionGhostDofs("TEST_PROBLEM"); CHKERRG(ierr);

  Vec v;
  ierr = m_field.getInterface<VecManager>()->vecCreateGhost("TEST_PROBLEM",ROW,&v);
  ierr = VecSetRandom(v,PETSC_NULL); CHKERRG(ierr);
  ierr = m_field.getInterface<VecManager>()->setLocalGhostVector("TEST_PROBLEM",ROW,v,INSERT_VALUES,SCATTER_REVERSE); CHKERRG(ierr);
  ierr = VecDestroy(&v); CHKERRG(ierr);


  typedef tee_device<std::ostream, std::ofstream> TeeDevice;
  typedef stream<TeeDevice> TeeStream;
  std::ofstream ofs("forces_and_sources_hdiv_continuity_check.txt");
  TeeDevice my_tee(std::cout, ofs);
  TeeStream my_split(my_tee);

  struct OpTetFluxes: public VolumeElementForcesAndSourcesCore::UserDataOperator {

    MoFEM::Interface &m_field;
    Tag tH;

    OpTetFluxes(MoFEM::Interface &m_field,Tag _th):
      VolumeElementForcesAndSourcesCore::UserDataOperator("HDIV",UserDataOperator::OPROW),
      m_field(m_field),tH(_th) {}

    MoFEMErrorCode doWork(
      int side,
      EntityType type,
      DataForcesAndSourcesCore::EntData &data) {
      MoFEMFunctionBeginHot;



      if(data.getFieldData().size()==0) MoFEMFunctionReturnHot(0);

      if(type == MBTRI) {

        boost::shared_ptr<const NumeredEntFiniteElement> mofem_fe = getNumeredEntFiniteElementPtr();
        SideNumber_multiIndex &side_table = mofem_fe->getSideNumberTable();
        EntityHandle face = side_table.get<1>().find(boost::make_tuple(type,side))->get()->ent;
        int sense = side_table.get<1>().find(boost::make_tuple(type,side))->get()->sense;

        // cerr << data.getHcurlN() << endl;


        VectorDouble t(3,0);
        int dd = 0;
        int nb_dofs = data.getHdivN().size2()/3;
        for(;dd<nb_dofs;dd++) {
          int ddd = 0;
          for(;ddd<3;ddd++) {
            t(ddd) += data.getHdivN(side)(dd,ddd)*data.getFieldData()[dd];
          }
        }

        double *t_ptr;
        rval = m_field.get_moab().tag_get_by_ptr(tH,&face,1,(const void **)&t_ptr); CHKERRG(rval);
        dd = 0;
        for(;dd<3;dd++) {
          t_ptr[dd] += sense*t[dd];
        }

      }

      MoFEMFunctionReturnHot(0);
    }

  };

  struct MyTetFE: public VolumeElementForcesAndSourcesCore {

    MyTetFE(MoFEM::Interface &m_field):
    VolumeElementForcesAndSourcesCore(m_field) {}
    int getRule(int order) { return -1; };

    MatrixDouble N_tri;
    MoFEMErrorCode setGaussPts(int order) {

      MoFEMFunctionBeginHot;

      try {

        N_tri.resize(1,3);
        ierr = ShapeMBTRI(&N_tri(0,0),G_TRI_X1,G_TRI_Y1,1); CHKERRG(ierr);

        gaussPts.resize(4,4);
        int ff = 0;
        for(;ff<4;ff++) {
          int dd = 0;
          for(;dd<3;dd++) {
            gaussPts(dd,ff) = cblas_ddot(3,&N_tri(0,0),1,&face_coords[ff][dd],3);
          }
          gaussPts(3,ff) = G_TRI_W1[0];
        }

        //std::cerr << gaussPts << std::endl;

      } catch (std::exception& ex) {
        std::ostringstream ss;
        ss << "thorw in method: " << ex.what() << " at line " << __LINE__ << " in file " << __FILE__;
        SETERRQ(PETSC_COMM_SELF,MOFEM_STD_EXCEPTION_THROW,ss.str().c_str());
      }

      MoFEMFunctionReturnHot(0);
    }


  };

  struct OpFacesSkinFluxes: public FaceElementForcesAndSourcesCore::UserDataOperator {

    MoFEM::Interface &m_field;
    Tag tH1,tH2;
    TeeStream &mySplit;

    OpFacesSkinFluxes(MoFEM::Interface &m_field,Tag _th1,Tag _th2,TeeStream &my_split):
      FaceElementForcesAndSourcesCore::UserDataOperator("HDIV",UserDataOperator::OPROW),
      m_field(m_field),tH1(_th1),tH2(_th2),mySplit(my_split) {}

    MoFEMErrorCode doWork(
      int side,
      EntityType type,
      DataForcesAndSourcesCore::EntData &data) {
      MoFEMFunctionBeginHot;



      if(type != MBTRI) MoFEMFunctionReturnHot(0);

      EntityHandle face = getNumeredEntFiniteElementPtr()->getEnt();

      double *t_ptr;
      rval = m_field.get_moab().tag_get_by_ptr(tH1,&face,1,(const void **)&t_ptr); CHKERRG(rval);
      double *tn_ptr;
      rval = m_field.get_moab().tag_get_by_ptr(tH2,&face,1,(const void **)&tn_ptr); CHKERRG(rval);

      *tn_ptr = getNormalsAtGaussPt()(0,0)*t_ptr[0]+getNormalsAtGaussPt()(0,1)*t_ptr[1]+getNormalsAtGaussPt()(0,2)*t_ptr[2];

      int nb_dofs = data.getHdivN().size2()/3;
      int dd = 0;
      for(;dd<nb_dofs;dd++) {
        *tn_ptr +=
        -getNormalsAtGaussPt()(0,0)*data.getHdivN()(0,3*dd+0)*data.getFieldData()[dd]
        -getNormalsAtGaussPt()(0,1)*data.getHdivN()(0,3*dd+1)*data.getFieldData()[dd]
        -getNormalsAtGaussPt()(0,2)*data.getHdivN()(0,3*dd+2)*data.getFieldData()[dd];
      }

      const double eps = 1e-8;
      if(fabs(*tn_ptr)>eps) {
        SETERRQ1(
          PETSC_COMM_SELF,
          MOFEM_ATOM_TEST_INVALID,
          "HDiv continuity failed %6.4e",
          *tn_ptr
        );
      }

      mySplit.precision(5);


      mySplit << face << " " /*<< std::fixed*/ << fabs(*tn_ptr) << std::endl;

      MoFEMFunctionReturnHot(0);
    }

  };

  struct OpFacesFluxes: public FaceElementForcesAndSourcesCore::UserDataOperator {

    MoFEM::Interface &m_field;
    Tag tH1,tH2;
    TeeStream &mySplit;

    OpFacesFluxes(MoFEM::Interface &m_field,Tag _th1,Tag _th2,TeeStream &my_split):
      FaceElementForcesAndSourcesCore::UserDataOperator("HDIV",UserDataOperator::OPROW),
      m_field(m_field),tH1(_th1),tH2(_th2),mySplit(my_split) {}

    MoFEMErrorCode doWork(
      int side,
      EntityType type,
      DataForcesAndSourcesCore::EntData &data) {
      MoFEMFunctionBeginHot;



      if(type != MBTRI) MoFEMFunctionReturnHot(0);

      EntityHandle face = getNumeredEntFiniteElementPtr()->getEnt();

      double *t_ptr;
      rval = m_field.get_moab().tag_get_by_ptr(tH1,&face,1,(const void **)&t_ptr); CHKERRG(rval);
      double *tn_ptr;
      rval = m_field.get_moab().tag_get_by_ptr(tH2,&face,1,(const void **)&tn_ptr); CHKERRG(rval);

      *tn_ptr = getNormalsAtGaussPt()(0,0)*t_ptr[0]+getNormalsAtGaussPt()(0,1)*t_ptr[1]+getNormalsAtGaussPt()(0,2)*t_ptr[2];

      const double eps = 1e-8;
      if(fabs(*tn_ptr)>eps) {
        SETERRQ1(
          PETSC_COMM_SELF,
          MOFEM_ATOM_TEST_INVALID,
          "HDiv continuity failed %6.4e",
          *tn_ptr
        );
      }

      mySplit.precision(5);

      mySplit << face << " " /*<< std::fixed*/ << fabs(*tn_ptr) << std::endl;

      MoFEMFunctionReturnHot(0);
    }

  };

  struct MyTriFE: public FaceElementForcesAndSourcesCore {

    MyTriFE(MoFEM::Interface &m_field): FaceElementForcesAndSourcesCore(m_field) {}
    int getRule(int order) { return -1; };

    MoFEMErrorCode setGaussPts(int order) {
      MoFEMFunctionBeginHot;

      gaussPts.resize(3,1);
      gaussPts(0,0) = G_TRI_X1[0];
      gaussPts(1,0) = G_TRI_Y1[0];
      gaussPts(2,0) = G_TRI_W1[0];

      MoFEMFunctionReturnHot(0);
    }

  };

  MyTetFE tet_fe(m_field);
  MyTriFE tri_fe(m_field);
  MyTriFE skin_fe(m_field);


  Tag th1;
  double def_val[] = {0,0,0};
  rval = moab.tag_get_handle("T",3,MB_TYPE_DOUBLE,th1,MB_TAG_CREAT|MB_TAG_SPARSE,&def_val); CHKERRG(rval);
  tet_fe.getOpPtrVector().push_back(new OpTetFluxes(m_field,th1));

  Tag th2;
  rval = moab.tag_get_handle("TN",1,MB_TYPE_DOUBLE,th2,MB_TAG_CREAT|MB_TAG_SPARSE,&def_val); CHKERRG(rval);
  tri_fe.getOpPtrVector().push_back(new OpFacesFluxes(m_field,th1,th2,my_split));
  skin_fe.getOpPtrVector().push_back(new OpFacesSkinFluxes(m_field,th1,th2,my_split));

  for(Range::iterator fit = faces.begin();fit!=faces.end();fit++) {
    rval = moab.tag_set_data(th1,&*fit,1,&def_val); CHKERRG(rval);
    rval = moab.tag_set_data(th2,&*fit,1,&def_val); CHKERRG(rval);
  }

  ierr = m_field.loop_finite_elements("TEST_PROBLEM","TET_FE",tet_fe);  CHKERRG(ierr);
  my_split << "intrnal\n";
  ierr = m_field.loop_finite_elements("TEST_PROBLEM","TRI_FE",tri_fe);  CHKERRG(ierr);
  my_split << "skin\n";
  ierr = m_field.loop_finite_elements("TEST_PROBLEM","SKIN_FE",skin_fe);  CHKERRG(ierr);

  EntityHandle meshset;
  rval = moab.create_meshset(MESHSET_SET,meshset); CHKERRG(rval);
  ierr = m_field.getInterface<BitRefManager>()->getEntitiesByTypeAndRefLevel(BitRefLevel().set(0),BitRefLevel().set(),MBTRI,meshset); CHKERRG(ierr);
  rval = moab.write_file("out.vtk","VTK","",&meshset,1); CHKERRG(rval);


  } catch (MoFEMException const &e) {
    SETERRQ(PETSC_COMM_SELF,e.errorCode,e.errorMessage);
  }

  ierr = PetscFinalize(); CHKERRG(ierr);
}
