/** \file forces_and_sources_getting_mult_H1_H1_atom_test.cpp
  \brief Atom test verifying forces and sources operator on H1 approx. space

*/


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
  ierr = m_field.add_field("FIELD1",H1,1); CHKERRQ(ierr);
  ierr = m_field.add_field("FIELD2",H1,3); CHKERRQ(ierr);

  //FE
  ierr = m_field.add_finite_element("TEST_FE"); CHKERRQ(ierr);

  //Define rows/cols and element data
  ierr = m_field.modify_finite_element_add_field_row("TEST_FE","FIELD1"); CHKERRQ(ierr);
  ierr = m_field.modify_finite_element_add_field_col("TEST_FE","FIELD2"); CHKERRQ(ierr);
  ierr = m_field.modify_finite_element_add_field_data("TEST_FE","FIELD1"); CHKERRQ(ierr);
  ierr = m_field.modify_finite_element_add_field_data("TEST_FE","FIELD2"); CHKERRQ(ierr);

  //Problem
  ierr = m_field.add_problem("TEST_PROBLEM"); CHKERRQ(ierr);

  //set finite elements for problem
  ierr = m_field.modify_problem_add_finite_element("TEST_PROBLEM","TEST_FE"); CHKERRQ(ierr);
  //set refinment level for problem
  ierr = m_field.modify_problem_ref_level_add_bit("TEST_PROBLEM",bit_level0); CHKERRQ(ierr);


  //meshset consisting all entities in mesh
  EntityHandle root_set = moab.get_root_set();
  //add entities to field
  ierr = m_field.add_ents_to_field_by_TETs(root_set,"FIELD1"); CHKERRQ(ierr);
  ierr = m_field.add_ents_to_field_by_TETs(root_set,"FIELD2"); CHKERRQ(ierr);
  //add entities to finite element
  ierr = m_field.add_ents_to_finite_element_by_TETs(root_set,"TEST_FE"); CHKERRQ(ierr);


  //set app. order
  //see Hierarchic Finite Element Bases on Unstructured Tetrahedral Meshes (Mark Ainsworth & Joe Coyle)
  int order = 5;
  ierr = m_field.set_field_order(root_set,MBTET,"FIELD1",order); CHKERRQ(ierr);
  ierr = m_field.set_field_order(root_set,MBTRI,"FIELD1",order); CHKERRQ(ierr);
  ierr = m_field.set_field_order(root_set,MBEDGE,"FIELD1",order); CHKERRQ(ierr);
  ierr = m_field.set_field_order(root_set,MBVERTEX,"FIELD1",1); CHKERRQ(ierr);
  ierr = m_field.set_field_order(root_set,MBTET,"FIELD2",order); CHKERRQ(ierr);
  ierr = m_field.set_field_order(root_set,MBTRI,"FIELD2",order); CHKERRQ(ierr);
  ierr = m_field.set_field_order(root_set,MBEDGE,"FIELD2",order); CHKERRQ(ierr);
  ierr = m_field.set_field_order(root_set,MBVERTEX,"FIELD2",1); CHKERRQ(ierr);

  /****/
  //build database
  //build field
  ierr = m_field.build_fields(); CHKERRQ(ierr);
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

  struct ForcesAndSurcesCore_TestFE: public ForcesAndSurcesCore {

    ErrorCode rval;
    PetscErrorCode ierr;

    typedef tee_device<ostream, ofstream> TeeDevice;
    typedef stream<TeeDevice> TeeStream;

    struct my_mult_H1_H1: public DataOperator {

      ofstream ofs;
      TeeDevice my_tee;
      TeeStream my_split;

      my_mult_H1_H1():
      ofs("forces_and_sources_getting_mult_H1_H1_atom_test.txt"),
      my_tee(cout, ofs),my_split(my_tee
      ) {};

      ~my_mult_H1_H1() {
        my_split.close();
      }

      ublas::matrix<FieldData> NN;

      PetscErrorCode doWork(
        int row_side,int col_side,
        EntityType row_type,EntityType col_type,
        DataForcesAndSurcesCore::EntData &row_data,
        DataForcesAndSurcesCore::EntData &col_data) {
          PetscFunctionBegin;

          int nb_row_dofs = row_data.getN().size2();
          int nb_col_dofs = col_data.getN().size2();

          my_split << row_side << " " << col_side << " " << row_type << " " << col_type << endl;
          my_split << "nb_row_dofs " << nb_row_dofs << " nb_col_dofs " << nb_col_dofs << endl;
          NN.resize(nb_row_dofs,nb_col_dofs);


          my_split.precision(2);
          my_split << row_data.getN() << endl;
          my_split << col_data.getN() << endl;

          for(unsigned int gg = 0;gg<row_data.getN().size1();gg++) {

            bzero(&*NN.data().begin(),nb_row_dofs*nb_col_dofs*sizeof(FieldData));

            cblas_dger(CblasRowMajor,
              nb_row_dofs,nb_col_dofs,
              1,&row_data.getN()(gg,0),1,&col_data.getN()(gg,0),1,
              &*NN.data().begin(),nb_col_dofs);

              my_split << "gg " << gg << " : ";
              my_split.precision(2);
              //my_split << NN << endl;
              my_split << NN - outer_prod(row_data.getN(gg),col_data.getN(gg)) << endl;
              if(row_type != MBVERTEX) {
                my_split << row_data.getDiffN(gg) << endl;
              }

              if(row_type == MBVERTEX) {
                my_split << row_data.getDiffN() << endl;
              } else {
                typedef ublas::array_adaptor<FieldData> storage_t;
                storage_t st(nb_row_dofs*3,&row_data.getDiffN()(gg,0));
                ublas::matrix<FieldData,ublas::row_major,storage_t> digNatGaussPt(nb_row_dofs,3,st);
                my_split << endl << digNatGaussPt << endl;
              }

            }

            my_split << endl;

            PetscFunctionReturn(0);
          }

        };

    my_mult_H1_H1 op;

    ForcesAndSurcesCore_TestFE(FieldInterface &_m_field):
      ForcesAndSurcesCore(_m_field),data_row(MBTET),data_col(MBTET) {};

    PetscErrorCode preProcess() {
      PetscFunctionBegin;
      PetscFunctionReturn(0);
    }

    DataForcesAndSurcesCore data_row,data_col;

    PetscErrorCode operator()() {
      PetscFunctionBegin;

      ierr = getSpacesAndBaseOnEntities(data_row); CHKERRQ(ierr);
      ierr = getSpacesAndBaseOnEntities(data_col); CHKERRQ(ierr);

      ierr = getEdgesSense(data_row); CHKERRQ(ierr);
      ierr = getTrisSense(data_row); CHKERRQ(ierr);
      ierr = getEdgesSense(data_col); CHKERRQ(ierr);
      ierr = getTrisSense(data_col); CHKERRQ(ierr);

      ierr = getEdgesDataOrder(data_row,H1); CHKERRQ(ierr);
      ierr = getEdgesDataOrder(data_col,H1); CHKERRQ(ierr);
      ierr = getTrisDataOrder(data_row,H1); CHKERRQ(ierr);
      ierr = getTrisDataOrder(data_col,H1); CHKERRQ(ierr);
      ierr = getTetsDataOrder(data_row,H1); CHKERRQ(ierr);
      ierr = getTetsDataOrder(data_col,H1); CHKERRQ(ierr);
      ierr = getRowNodesIndices(data_row,"FIELD1"); CHKERRQ(ierr);
      ierr = getColNodesIndices(data_row,"FIELD2"); CHKERRQ(ierr);
      ierr = getEdgesRowIndices(data_row,"FIELD1"); CHKERRQ(ierr);
      ierr = getEdgesColIndices(data_col,"FIELD2"); CHKERRQ(ierr);
      ierr = getTrisRowIndices(data_row,"FIELD1"); CHKERRQ(ierr);
      ierr = getTrisColIndices(data_col,"FIELD2"); CHKERRQ(ierr);
      ierr = getTetsRowIndices(data_row,"FIELD1"); CHKERRQ(ierr);
      ierr = getTetsColIndices(data_col,"FIELD2"); CHKERRQ(ierr);
      ierr = getFaceTriNodes(data_row); CHKERRQ(ierr);
      ierr = getFaceTriNodes(data_col); CHKERRQ(ierr);

      data_row.dataOnEntities[MBVERTEX][0].getN().resize(4,4,false);
      ierr = ShapeMBTET(
        &*data_row.dataOnEntities[MBVERTEX][0].getN().data().begin(),G_TET_X4,G_TET_Y4,G_TET_Z4,4
      ); CHKERRQ(ierr);
      data_col.dataOnEntities[MBVERTEX][0].getN() = data_row.dataOnEntities[MBVERTEX][0].getN();
      ierr = shapeTETFunctions_H1(data_row,G_TET_X4,G_TET_Y4,G_TET_Z4,4,Legendre_polynomials); CHKERRQ(ierr);
      ierr = shapeTETFunctions_H1(data_col,G_TET_X4,G_TET_Y4,G_TET_Z4,4,Legendre_polynomials); CHKERRQ(ierr);

      try {
	ierr = op.opLhs(data_row,data_col,true); CHKERRQ(ierr);
      } catch (exception& ex) {
	ostringstream ss;
	ss << "thorw in method: " << ex.what() << " at line " << __LINE__ << " in file " << __FILE__ << endl;
	SETERRQ(PETSC_COMM_SELF,1,ss.str().c_str());
      }

      PetscFunctionReturn(0);
    }

    PetscErrorCode postProcess() {
      PetscFunctionBegin;

      PetscFunctionReturn(0);
    }


  };

  ForcesAndSurcesCore_TestFE fe1(m_field);
  ierr = m_field.loop_finite_elements("TEST_PROBLEM","TEST_FE",fe1);  CHKERRQ(ierr);

  ierr = PetscFinalize(); CHKERRQ(ierr);

  return 0;

}
