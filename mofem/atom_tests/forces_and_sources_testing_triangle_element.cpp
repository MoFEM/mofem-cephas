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
#include <Projection10NodeCoordsOnField.hpp>

#include <moab/Skinner.hpp>

#include <boost/iostreams/tee.hpp>
#include <boost/iostreams/stream.hpp>
#include <fstream>
#include <iostream>

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

  PetscBool flg = PETSC_TRUE;
  char mesh_file_name[255];
  #if PETSC_VERSION_GE(3,6,4)
  ierr = PetscOptionsGetString(PETSC_NULL,"","-my_file",mesh_file_name,255,&flg); CHKERRQ(ierr);
  #else
  ierr = PetscOptionsGetString(PETSC_NULL,"-my_file",mesh_file_name,255,&flg); CHKERRQ(ierr);
  #endif
  if(flg != PETSC_TRUE) {
    SETERRQ(PETSC_COMM_SELF,1,"*** ERROR -my_file (MESH FILE NEEDED)");
  }

  //Create MoFEM (Joseph) database
  MoFEM::Core core(moab);
  FieldInterface& m_field = core;

  ParallelComm* pcomm = ParallelComm::get_pcomm(&moab,MYPCOMM_INDEX);
  if(pcomm == NULL) pcomm =  new ParallelComm(&moab,PETSC_COMM_WORLD);

  const char *option;
  option = "";//"PARALLEL=BCAST;";//;DEBUG_IO";
  rval = moab.load_file(mesh_file_name, 0, option); CHKERRQ_MOAB(rval);

  //set entitities bit level
  BitRefLevel bit_level0;
  bit_level0.set(0);
  EntityHandle meshset_level0;
  rval = moab.create_meshset(MESHSET_SET,meshset_level0); CHKERRQ_MOAB(rval);
  ierr = m_field.seed_ref_level_3D(0,bit_level0); CHKERRQ(ierr);

  //Fields
  ierr = m_field.add_field("FIELD1",H1,3); CHKERRQ(ierr);
  ierr = m_field.add_field("FIELD2",NOFIELD,3); CHKERRQ(ierr);

  {
    // Creating and adding no field entities.
    const double coords[] = {0,0,0};
    EntityHandle no_field_vertex;
    rval = m_field.get_moab().create_vertex(coords,no_field_vertex); CHKERRQ_MOAB(rval);
    Range range_no_field_vertex;
    range_no_field_vertex.insert(no_field_vertex);
    ierr = m_field.seed_ref_level(range_no_field_vertex,BitRefLevel().set()); CHKERRQ(ierr);
    EntityHandle meshset = m_field.get_field_meshset("FIELD2");
    rval = m_field.get_moab().add_entities(meshset,range_no_field_vertex); CHKERRQ_MOAB(rval);
  }

  //FE
  ierr = m_field.add_finite_element("TEST_FE1"); CHKERRQ(ierr);
  ierr = m_field.add_finite_element("TEST_FE2"); CHKERRQ(ierr);

  //Define rows/cols and element data
  ierr = m_field.modify_finite_element_add_field_row("TEST_FE1","FIELD1"); CHKERRQ(ierr);
  ierr = m_field.modify_finite_element_add_field_col("TEST_FE1","FIELD1"); CHKERRQ(ierr);
  ierr = m_field.modify_finite_element_add_field_data("TEST_FE1","FIELD1"); CHKERRQ(ierr);

  ierr = m_field.modify_finite_element_add_field_row("TEST_FE2","FIELD2"); CHKERRQ(ierr);
  ierr = m_field.modify_finite_element_add_field_col("TEST_FE2","FIELD1"); CHKERRQ(ierr);
  ierr = m_field.modify_finite_element_add_field_data("TEST_FE2","FIELD1"); CHKERRQ(ierr);
  ierr = m_field.modify_finite_element_add_field_data("TEST_FE2","FIELD2"); CHKERRQ(ierr);

  //Problem
  ierr = m_field.add_problem("TEST_PROBLEM"); CHKERRQ(ierr);

  //set finite elements for problem
  ierr = m_field.modify_problem_add_finite_element("TEST_PROBLEM","TEST_FE1"); CHKERRQ(ierr);
  ierr = m_field.modify_problem_add_finite_element("TEST_PROBLEM","TEST_FE2"); CHKERRQ(ierr);
  //set refinment level for problem
  ierr = m_field.modify_problem_ref_level_add_bit("TEST_PROBLEM",bit_level0); CHKERRQ(ierr);

  //meshset consisting all entities in mesh
  EntityHandle root_set = moab.get_root_set();
  //add entities to field
  ierr = m_field.add_ents_to_field_by_TETs(root_set,"FIELD1"); CHKERRQ(ierr);
  //add entities to finite element
  Range tets;
  rval = moab.get_entities_by_type(0,MBTET,tets,false); CHKERRQ_MOAB(rval);
  Skinner skin(&m_field.get_moab());
  Range tets_skin;
  rval = skin.find_skin(0,tets,false,tets_skin); CHKERR_MOAB(rval);
  ierr = m_field.add_ents_to_finite_element_by_TRIs(tets_skin,"TEST_FE1"); CHKERRQ(ierr);
  ierr = m_field.add_ents_to_finite_element_by_TRIs(tets_skin,"TEST_FE2"); CHKERRQ(ierr);

  //set app. order
  //see Hierarchic Finite Element Bases on Unstructured Tetrahedral Meshes (Mark Ainsworth & Joe Coyle)
  int order = 3;
  ierr = m_field.set_field_order(root_set,MBTET,"FIELD1",order); CHKERRQ(ierr);
  ierr = m_field.set_field_order(root_set,MBTRI,"FIELD1",order); CHKERRQ(ierr);
  ierr = m_field.set_field_order(root_set,MBEDGE,"FIELD1",order); CHKERRQ(ierr);
  ierr = m_field.set_field_order(root_set,MBVERTEX,"FIELD1",1); CHKERRQ(ierr);

  /****/
  //build database
  //build field
  ierr = m_field.build_fields(); CHKERRQ(ierr);
  //set FIELD1 from positions of 10 node tets
  Projection10NodeCoordsOnField ent_method(m_field,"FIELD1");
  ierr = m_field.loop_dofs("FIELD1",ent_method); CHKERRQ(ierr);
  //build finite elemnts
  ierr = m_field.build_finite_elements(); CHKERRQ(ierr);
  //build adjacencies
  ierr = m_field.build_adjacencies(bit_level0); CHKERRQ(ierr);
  //build problem
  ierr = m_field.build_problems(); CHKERRQ(ierr);

  /****/
  //mesh partitioning
  //partition
  ierr = m_field.partition_simple_problem("TEST_PROBLEM"); CHKERRQ(ierr);
  ierr = m_field.partition_finite_elements("TEST_PROBLEM"); CHKERRQ(ierr);
  //what are ghost nodes, see Petsc Manual
  ierr = m_field.partition_ghost_dofs("TEST_PROBLEM"); CHKERRQ(ierr);

  typedef tee_device<ostream, ofstream> TeeDevice;
  typedef stream<TeeDevice> TeeStream;

  ofstream ofs("forces_and_sources_testing_triangle_element.txt");
  TeeDevice my_tee(cout, ofs);
  TeeStream my_split(my_tee);

  struct MyOp1: public FaceElementForcesAndSourcesCore::UserDataOperator {

    TeeStream &my_split;
    MyOp1(TeeStream &_my_split,const char type):
      FaceElementForcesAndSourcesCore::UserDataOperator("FIELD1","FIELD1",type),
      my_split(_my_split) {}

    PetscErrorCode doWork(
      int side,
      EntityType type,
      DataForcesAndSurcesCore::EntData &data) {
      PetscFunctionBegin;

      const double eps = 1e-4;
      for(
        ublas::unbounded_array<double>::iterator it = getNormal().data().begin();
        it!=getNormal().data().end();it++
      ) {
        *it = fabs(*it)<eps ? 0.0 : *it;
      }

      my_split << "NH1" << endl;
      my_split << "side: " << side << " type: " << type << endl;
      my_split << "data: " << data << endl;
      my_split << setprecision(3) << getCoords() << endl;
      my_split << setprecision(3) << getCoordsAtGaussPts() << endl;
      my_split << setprecision(3) << getArea() << endl;
      my_split << setprecision(3) << getNormal() << endl;
      my_split << setprecision(3) << getHoCoordsAtGaussPts() << endl;
      my_split << setprecision(3) << getNormals_at_GaussPt() << endl;
      my_split << setprecision(3) << getTangent1_at_GaussPt() << endl;
      my_split << setprecision(3) << getTangent2_at_GaussPt() << endl;
      my_split << endl;
      PetscFunctionReturn(0);
    }

    PetscErrorCode doWork(
      int row_side,int col_side,
      EntityType row_type,EntityType col_type,
      DataForcesAndSurcesCore::EntData &row_data,
      DataForcesAndSurcesCore::EntData &col_data
    ) {

      PetscFunctionBegin;
      my_split << "NH1NH1" << endl;
      my_split << "row side: " << row_side << " row_type: " << row_type << endl;
      my_split << row_data << endl;
      my_split << "NH1NH1" << endl;
      my_split << "col side: " << col_side << " col_type: " << col_type << endl;
      my_split << row_data << endl;


      PetscErrorCode ierr;
      VectorInt row_indices,col_indices;
      ierr = getPorblemRowIndices("FIELD1",row_type,row_side,row_indices); CHKERRQ(ierr);
      ierr = getPorblemColIndices("FIELD1",col_type,col_side,col_indices); CHKERRQ(ierr);

      if(row_indices.size()!=row_data.getIndices().size()) {
        SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"row inconsistency");
      }

      if(col_indices.size()!=col_data.getIndices().size()) {
        SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"col inconsistency");
      }

      for(unsigned int rr = 0;rr<row_indices.size();rr++) {
        if(row_indices[rr] != row_data.getIndices()[rr]) {
          cerr << row_indices << endl;
          cerr << row_data.getIndices() << endl;
          SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"row inconsistency");
        }
      }

      for(unsigned int cc = 0;cc<col_indices.size();cc++) {
        if(col_indices[cc] != col_data.getIndices()[cc]) {
          cerr << col_indices << endl;
          cerr << col_data.getIndices() << endl;
          SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"row inconsistency");
        }
      }

      my_split << row_data << endl;

      PetscFunctionReturn(0);
    }

  };

  struct MyOp2: public FaceElementForcesAndSourcesCore::UserDataOperator {

    TeeStream &my_split;
    MyOp2(TeeStream &_my_split,const char type):
    FaceElementForcesAndSourcesCore::UserDataOperator("FIELD2","FIELD1",type),
      my_split(_my_split) {}

    PetscErrorCode doWork(
      int side,
      EntityType type,
      DataForcesAndSurcesCore::EntData &data) {
      PetscFunctionBegin;

      if(type != MBENTITYSET) PetscFunctionReturn(0);

      my_split << "NOFIELD" << endl;
      my_split << "side: " << side << " type: " << type << endl;
      my_split << data << endl;
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

      my_split << "NOFILEDH1" << endl;
      my_split << "row side: " << row_side << " row_type: " << row_type << endl;
      my_split << row_data << endl;
      my_split << "col side: " << col_side << " col_type: " << col_type << endl;
      my_split << col_data << endl;

      PetscFunctionReturn(0);
    }

  };

  FaceElementForcesAndSourcesCore fe1(m_field);
  fe1.getOpPtrVector().push_back(new MyOp1(my_split,ForcesAndSurcesCore::UserDataOperator::OPROW));
  fe1.getOpPtrVector().push_back(new MyOp1(my_split,ForcesAndSurcesCore::UserDataOperator::OPROWCOL));

  FaceElementForcesAndSourcesCore fe2(m_field);
  fe2.getOpPtrVector().push_back(new MyOp2(my_split,ForcesAndSurcesCore::UserDataOperator::OPROW));
  fe2.getOpPtrVector().push_back(new MyOp2(my_split,ForcesAndSurcesCore::UserDataOperator::OPROWCOL));

  ierr = m_field.loop_finite_elements("TEST_PROBLEM","TEST_FE1",fe1);  CHKERRQ(ierr);
  ierr = m_field.loop_finite_elements("TEST_PROBLEM","TEST_FE2",fe2);  CHKERRQ(ierr);

  ierr = PetscFinalize(); CHKERRQ(ierr);

  return 0;

}
