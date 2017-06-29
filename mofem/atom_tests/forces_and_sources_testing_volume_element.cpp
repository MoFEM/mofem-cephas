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

  ErrorCode rval;
  PetscErrorCode ierr;

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
    SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"*** ERROR -my_file (MESH FILE NEEDED)");
  }

  //Create MoFEM (Joseph) database
  MoFEM::Core core(moab);
  MoFEM::Interface& m_field = core;

  ParallelComm* pcomm = ParallelComm::get_pcomm(&moab,MYPCOMM_INDEX);
  if(pcomm == NULL) pcomm =  new ParallelComm(&moab,PETSC_COMM_WORLD);

  const char *option;
  option = "";//"PARALLEL=BCAST;";//;DEBUG_IO";
  BARRIER_RANK_START(pcomm)
  rval = moab.load_file(mesh_file_name, 0, option); CHKERRQ_MOAB(rval);
  BARRIER_RANK_END(pcomm)

  //set entitities bit level
  BitRefLevel bit_level0;
  bit_level0.set(0);
  EntityHandle meshset_level0;
  rval = moab.create_meshset(MESHSET_SET,meshset_level0); CHKERRQ_MOAB(rval);
  ierr = m_field.seed_ref_level_3D(0,bit_level0); CHKERRQ(ierr);

  //Fields
  ierr = m_field.add_field("FIELD1",H1,AINSWORTH_LEGENDRE_BASE,1); CHKERRQ(ierr);
  ierr = m_field.add_field("FIELD2",H1,AINSWORTH_LEGENDRE_BASE,3); CHKERRQ(ierr);
  ierr = m_field.add_field("FIELD3",NOFIELD,AINSWORTH_LEGENDRE_BASE,3); CHKERRQ(ierr);
  ierr = m_field.add_field("MESH_NODE_POSITIONS",H1,AINSWORTH_LEGENDRE_BASE,3); CHKERRQ(ierr);

  {
    // Creating and adding no field entities.
    const double coords[] = {0,0,0};
    EntityHandle no_field_vertex;
    rval = m_field.get_moab().create_vertex(coords,no_field_vertex); CHKERRQ_MOAB(rval);
    Range range_no_field_vertex;
    range_no_field_vertex.insert(no_field_vertex);
    ierr = m_field.seed_ref_level(range_no_field_vertex,BitRefLevel().set()); CHKERRQ(ierr);
    EntityHandle meshset = m_field.get_field_meshset("FIELD3");
    rval = m_field.get_moab().add_entities(meshset,range_no_field_vertex); CHKERRQ_MOAB(rval);
  }

  //FE
  ierr = m_field.add_finite_element("TEST_FE1"); CHKERRQ(ierr);
  ierr = m_field.add_finite_element("TEST_FE2"); CHKERRQ(ierr);

  //Define rows/cols and element data
  ierr = m_field.modify_finite_element_add_field_row("TEST_FE1","FIELD1"); CHKERRQ(ierr);
  ierr = m_field.modify_finite_element_add_field_col("TEST_FE1","FIELD2"); CHKERRQ(ierr);
  ierr = m_field.modify_finite_element_add_field_data("TEST_FE1","FIELD1"); CHKERRQ(ierr);
  ierr = m_field.modify_finite_element_add_field_data("TEST_FE1","FIELD2"); CHKERRQ(ierr);
  ierr = m_field.modify_finite_element_add_field_data("TEST_FE1","MESH_NODE_POSITIONS"); CHKERRQ(ierr);

  ierr = m_field.modify_finite_element_add_field_row("TEST_FE2","FIELD3"); CHKERRQ(ierr);
  ierr = m_field.modify_finite_element_add_field_col("TEST_FE2","FIELD1"); CHKERRQ(ierr);
  ierr = m_field.modify_finite_element_add_field_data("TEST_FE2","FIELD1"); CHKERRQ(ierr);
  ierr = m_field.modify_finite_element_add_field_data("TEST_FE2","FIELD3"); CHKERRQ(ierr);
  ierr = m_field.modify_finite_element_add_field_data("TEST_FE2","MESH_NODE_POSITIONS"); CHKERRQ(ierr);

  //Problem
  ierr = m_field.add_problem("TEST_PROBLEM"); CHKERRQ(ierr);

  //set finite elements for problem
  ierr = m_field.modify_problem_add_finite_element("TEST_PROBLEM","TEST_FE1"); CHKERRQ(ierr);
  ierr = m_field.modify_problem_add_finite_element("TEST_PROBLEM","TEST_FE2"); CHKERRQ(ierr);
  //set refinement level for problem
  ierr = m_field.modify_problem_ref_level_add_bit("TEST_PROBLEM",bit_level0); CHKERRQ(ierr);


  //meshset consisting all entities in mesh
  EntityHandle root_set = moab.get_root_set();
  //add entities to field
  ierr = m_field.add_ents_to_field_by_type(root_set,MBTET,"FIELD1"); CHKERRQ(ierr);
  ierr = m_field.add_ents_to_field_by_type(root_set,MBTET,"FIELD2"); CHKERRQ(ierr);
  ierr = m_field.add_ents_to_field_by_type(root_set,MBTET,"MESH_NODE_POSITIONS"); CHKERRQ(ierr);

  //add entities to finite element
  ierr = m_field.add_ents_to_finite_element_by_type(root_set,MBTET,"TEST_FE1"); CHKERRQ(ierr);
  ierr = m_field.add_ents_to_finite_element_by_type(root_set,MBTET,"TEST_FE2"); CHKERRQ(ierr);

  //set app. order
  //see Hierarchic Finite Element Bases on Unstructured Tetrahedral Meshes (Mark Ainsworth & Joe Coyle)
  int order = 4;
  ierr = m_field.set_field_order(root_set,MBTET,"FIELD1",order); CHKERRQ(ierr);
  ierr = m_field.set_field_order(root_set,MBTRI,"FIELD1",order); CHKERRQ(ierr);
  ierr = m_field.set_field_order(root_set,MBEDGE,"FIELD1",order); CHKERRQ(ierr);
  ierr = m_field.set_field_order(root_set,MBVERTEX,"FIELD1",1); CHKERRQ(ierr);
  ierr = m_field.set_field_order(root_set,MBTET,"FIELD2",order); CHKERRQ(ierr);
  ierr = m_field.set_field_order(root_set,MBTRI,"FIELD2",order); CHKERRQ(ierr);
  ierr = m_field.set_field_order(root_set,MBEDGE,"FIELD2",order); CHKERRQ(ierr);
  ierr = m_field.set_field_order(root_set,MBVERTEX,"FIELD2",1); CHKERRQ(ierr);

  ierr = m_field.set_field_order(root_set,MBTET,"MESH_NODE_POSITIONS",2); CHKERRQ(ierr);
  ierr = m_field.set_field_order(root_set,MBTRI,"MESH_NODE_POSITIONS",2); CHKERRQ(ierr);
  ierr = m_field.set_field_order(root_set,MBEDGE,"MESH_NODE_POSITIONS",2); CHKERRQ(ierr);
  ierr = m_field.set_field_order(root_set,MBVERTEX,"MESH_NODE_POSITIONS",1); CHKERRQ(ierr);

  /****/
  //build database
  //build field
  ierr = m_field.build_fields(); CHKERRQ(ierr);
  //build finite elemnts
  ierr = m_field.build_finite_elements(); CHKERRQ(ierr);
  //build adjacencies
  ierr = m_field.build_adjacencies(bit_level0); CHKERRQ(ierr);
  //build problem
  ProblemsManager *prb_mng_ptr;
  ierr = m_field.query_interface(prb_mng_ptr); CHKERRQ(ierr);
    const Problem_multiIndex *problems_ptr;
  ierr = prb_mng_ptr->buildProblem("TEST_PROBLEM",false); CHKERRQ(ierr);

  /****/
  //mesh partitioning
  //partition
  ierr = prb_mng_ptr->partitionSimpleProblem("TEST_PROBLEM"); CHKERRQ(ierr);
  ierr = prb_mng_ptr->partitionFiniteElements("TEST_PROBLEM"); CHKERRQ(ierr);
  //what are ghost nodes, see Petsc Manual
  ierr = prb_mng_ptr->partitionGhostDofs("TEST_PROBLEM"); CHKERRQ(ierr);

  //set from positions of 10 node tets
  Projection10NodeCoordsOnField ent_method(m_field,"MESH_NODE_POSITIONS");
  ierr = m_field.loop_dofs("MESH_NODE_POSITIONS",ent_method); CHKERRQ(ierr);


  typedef tee_device<std::ostream, std::ofstream> TeeDevice;
  typedef stream<TeeDevice> TeeStream;

  std::ofstream ofs("forces_and_sources_testing_volume_element.txt");
  TeeDevice my_tee(std::cout, ofs);
  TeeStream my_split(my_tee);

  struct MyOp1: public VolumeElementForcesAndSourcesCore::UserDataOperator {

    TeeStream &my_split;
    MyOp1(TeeStream &_my_split,char type):
      VolumeElementForcesAndSourcesCore::UserDataOperator("FIELD1","FIELD2",type),
      my_split(_my_split) {}

    PetscErrorCode doWork(
      int side,
      EntityType type,
      DataForcesAndSurcesCore::EntData &data) {
      PetscFunctionBegin;
      my_split << "NH1" << std::endl;
      my_split << "side: " << side << " type: " << type << std::endl;
      my_split << data << std::endl;
      my_split << std::setprecision(3) << getVolume() << std::endl;
      my_split << std::setprecision(3) << getCoords() << std::endl;
      my_split << std::setprecision(3) << getCoordsAtGaussPts() << std::endl;
      PetscFunctionReturn(0);
    }

    PetscErrorCode doWork(
      int row_side,int col_side,
      EntityType row_type,EntityType col_type,
      DataForcesAndSurcesCore::EntData &row_data,
      DataForcesAndSurcesCore::EntData &col_data) {
      PetscFunctionBegin;
      my_split << "NH1NH1" << std::endl;
      my_split << "row side: " << row_side << " row_type: " << row_type << std::endl;
      my_split << row_data << std::endl;
      my_split << "NH1NH1" << std::endl;
      my_split << "col side: " << col_side << " col_type: " << col_type << std::endl;
      my_split << col_data << std::endl;

      PetscErrorCode ierr;
      VectorInt row_indices,col_indices;
      ierr = getPorblemRowIndices("FIELD1",row_type,row_side,row_indices); CHKERRQ(ierr);
      ierr = getPorblemColIndices("FIELD2",col_type,col_side,col_indices); CHKERRQ(ierr);

      if(row_indices.size()!=row_data.getIndices().size()) {
        std::cerr << row_indices << std::endl;
        std::cerr << row_data.getIndices() << std::endl;
        SETERRQ2(
          PETSC_COMM_SELF,
          MOFEM_DATA_INCONSISTENCY,
          "row inconsistency %d != %d",
          row_indices.size(),
          row_data.getIndices().size()
        );
      }

      if(col_indices.size()!=col_data.getIndices().size()) {
        SETERRQ2(
          PETSC_COMM_SELF,
          MOFEM_DATA_INCONSISTENCY,
          "row inconsistency %d != %d",
          col_indices.size(),
          col_data.getIndices().size()
        );
      }

      for(unsigned int rr = 0;rr<row_indices.size();rr++) {
        if(row_indices[rr] != row_data.getIndices()[rr]) {
          std::cerr << row_indices << std::endl;
          std::cerr << row_data.getIndices() << std::endl;
          SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"row inconsistency");
        }
      }

      for(unsigned int cc = 0;cc<col_indices.size();cc++) {
        if(col_indices[cc] != col_data.getIndices()[cc]) {
          std::cerr << col_indices << std::endl;
          std::cerr << col_data.getIndices() << std::endl;
          SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"row inconsistency");
        }
      }

      PetscFunctionReturn(0);
    }

  };

  struct MyOp2: public VolumeElementForcesAndSourcesCore::UserDataOperator {

    TeeStream &my_split;
    MyOp2(TeeStream &_my_split,char type):
      VolumeElementForcesAndSourcesCore::UserDataOperator("FIELD3","FIELD1",type),
      my_split(_my_split) {}

    PetscErrorCode doWork(
      int side,
      EntityType type,
      DataForcesAndSurcesCore::EntData &data) {
      PetscFunctionBegin;

      if(type != MBENTITYSET) PetscFunctionReturn(0);

      my_split << "NOFIELD" << std::endl;
      my_split << "side: " << side << " type: " << type << std::endl;
      my_split << data << std::endl;
      PetscFunctionReturn(0);
    }

    PetscErrorCode doWork(
      int row_side,int col_side,
      EntityType row_type,EntityType col_type,
      DataForcesAndSurcesCore::EntData &row_data,
      DataForcesAndSurcesCore::EntData &col_data) {
      PetscFunctionBegin;

      unSetSymm();

      if(row_type != MBENTITYSET) PetscFunctionReturn(0);

      my_split << "NOFILEDH1" << std::endl;
      my_split << "row side: " << row_side << " row_type: " << row_type << std::endl;
      my_split << row_data << std::endl;
      my_split << "col side: " << col_side << " col_type: " << col_type << std::endl;
      my_split << col_data << std::endl;

      PetscFunctionReturn(0);
    }

  };

  VolumeElementForcesAndSourcesCore fe1(m_field);
  fe1.getOpPtrVector().push_back(new MyOp1(my_split,ForcesAndSurcesCore::UserDataOperator::OPROW));
  fe1.getOpPtrVector().push_back(new MyOp1(my_split,ForcesAndSurcesCore::UserDataOperator::OPROWCOL));

  VolumeElementForcesAndSourcesCore fe2(m_field);
  fe2.getOpPtrVector().push_back(new MyOp2(my_split,ForcesAndSurcesCore::UserDataOperator::OPROW));
  fe2.getOpPtrVector().push_back(new MyOp2(my_split,ForcesAndSurcesCore::UserDataOperator::OPROWCOL));

  ierr = m_field.loop_finite_elements("TEST_PROBLEM","TEST_FE1",fe1);  CHKERRQ(ierr);
  ierr = m_field.loop_finite_elements("TEST_PROBLEM","TEST_FE2",fe2);  CHKERRQ(ierr);


  } catch (MoFEMException const &e) {
    SETERRQ(PETSC_COMM_SELF,e.errorCode,e.errorMessage);
  }

  ierr = PetscFinalize(); CHKERRQ(ierr);

  return 0;

}
