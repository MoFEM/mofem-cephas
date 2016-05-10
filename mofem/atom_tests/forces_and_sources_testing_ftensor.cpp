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

#define FTENSOR_DEBUG

#include <MoFEM.hpp>
#include <Projection10NodeCoordsOnField.hpp>

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
    SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"*** ERROR -my_file (MESH FILE NEEDED)");
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
  ierr = m_field.add_field("FIELD1",H1,3); CHKERRQ(ierr);

  //FE
  ierr = m_field.add_finite_element("TEST_FE1"); CHKERRQ(ierr);

  //Define rows/cols and element data
  ierr = m_field.modify_finite_element_add_field_col("TEST_FE1","FIELD1"); CHKERRQ(ierr);
  ierr = m_field.modify_finite_element_add_field_data("TEST_FE1","FIELD1"); CHKERRQ(ierr);

  //Problem
  ierr = m_field.add_problem("TEST_PROBLEM"); CHKERRQ(ierr);

  //set finite elements for problem
  ierr = m_field.modify_problem_add_finite_element("TEST_PROBLEM","TEST_FE1"); CHKERRQ(ierr);
  //set refinment level for problem
  ierr = m_field.modify_problem_ref_level_add_bit("TEST_PROBLEM",bit_level0); CHKERRQ(ierr);


  //meshset consisting all entities in mesh
  EntityHandle root_set = moab.get_root_set();
  //add entities to field
  ierr = m_field.add_ents_to_field_by_TETs(root_set,"FIELD1"); CHKERRQ(ierr);

  //add entities to finite element
  ierr = m_field.add_ents_to_finite_element_by_TETs(root_set,"TEST_FE1"); CHKERRQ(ierr);

  //set app. order
  //see Hierarchic Finite Element Bases on Unstructured Tetrahedral Meshes (Mark Ainsworth & Joe Coyle)
  int order = 4;
  ierr = m_field.set_field_order(root_set,MBTET,"FIELD1",order); CHKERRQ(ierr);
  ierr = m_field.set_field_order(root_set,MBTRI,"FIELD1",order); CHKERRQ(ierr);
  ierr = m_field.set_field_order(root_set,MBEDGE,"FIELD1",order); CHKERRQ(ierr);
  ierr = m_field.set_field_order(root_set,MBVERTEX,"FIELD1",1); CHKERRQ(ierr);

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

  typedef tee_device<std::ostream, std::ofstream> TeeDevice;
  typedef stream<TeeDevice> TeeStream;

  std::ofstream ofs("forces_and_sources_testing_ftensor.txt");
  TeeDevice my_tee(std::cout, ofs);
  TeeStream my_split(my_tee);

  struct MyOp1: public VolumeElementForcesAndSourcesCore::UserDataOperator {

    TeeStream &my_split;
    MyOp1(TeeStream &_my_split,char type):
    VolumeElementForcesAndSourcesCore::UserDataOperator("FIELD1","FIELD1",type),
    my_split(_my_split) {}


    PetscErrorCode doWork(
      int side,
      EntityType type,
      DataForcesAndSurcesCore::EntData &data) {
      PetscFunctionBegin;

      FTensor::Tensor2<double,3,3,row_major> t2;
      const int nb_gauss_pts = data.getN().size1();
      const int nb_base_functions = data.getN().size2();

      FTensor::Tensor0<double*> base_function = data.getFTensorN();
      Tensor1<double*,3> diff_base = data.getFTensorDiffN<3>();

      for(int gg = 0;gg!=nb_gauss_pts;gg++) {

        for(int bb = 0;bb!=nb_base_functions;bb++) {

          double base_val = base_function;
          // my_split << base_val << endl;
          if(base_val!=data.getN(gg)[bb]) {
            SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"Data inconsistency");
          }
          ++base_function;

          for(int dd = 0;dd<3;dd++) {
            if(diff_base(dd)!=data.getDiffN()(gg,3*bb+dd)) {
              SETERRQ2(
                PETSC_COMM_SELF,
                MOFEM_DATA_INCONSISTENCY,
                "Data inconsistency gg = %d bb = %d",
                gg,
                bb
              );
            }
          }

          FTensor::Index<'I',3> I;
          FTensor::Index<'J',3> J;

          t2(I,J) = diff_base(I)*diff_base(J);

          MatrixAdaptor mat = MatrixAdaptor(3,3,ublas::shallow_array_adaptor<double>(9,&t2(0,0)));

          for(int II = 0;II!=3;II++) {
            for(int JJ = 0;JJ!=3;JJ++) {
              if(t2(II,JJ)-mat(II,JJ)!=0) {
                SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"Data inconsistency");
              }
              //my_split << t2(II,JJ)-mat(II,JJ) << " ";
            }
            //my_split << endl;
          }
          //my_split << endl;
          //
          // double t0;
          // t0 = t2(I,I);
          // my_split << "trace " << t0 << endl;

          ++diff_base;

        }
      }
      PetscFunctionReturn(0);
    }

    PetscErrorCode doWork(
      int row_side,int col_side,
      EntityType row_type,EntityType col_type,
      DataForcesAndSurcesCore::EntData &row_data,
      DataForcesAndSurcesCore::EntData &col_data
    ) {
      PetscFunctionBegin;
      PetscErrorCode ierr;

      const int nb_gauss_pts = row_data.getN().size1();
      const int nb_base_functions_row = col_data.getN().size2();
      const int nb_base_functions_col = col_data.getN().size2();

      FTensor::Number<0> N0;
      FTensor::Number<1> N1;
      FTensor::Number<2> N2;

      FTensor::Index<'I',3> I;
      FTensor::Index<'J',3> J;

      for(int gg = 0;gg!=nb_gauss_pts;gg++) {

        double vol = getVolume();
        double weight = getGaussPts()(3,gg);

        for(int bb_row = 0;bb_row!=nb_base_functions_row;bb_row++) {


        }

      }

      PetscFunctionReturn(0);
    }

  };


  VolumeElementForcesAndSourcesCore fe1(m_field);
  fe1.getOpPtrVector().push_back(new MyOp1(my_split,ForcesAndSurcesCore::UserDataOperator::OPROW));
  fe1.getOpPtrVector().push_back(new MyOp1(my_split,ForcesAndSurcesCore::UserDataOperator::OPROWCOL));

  ierr = m_field.loop_finite_elements("TEST_PROBLEM","TEST_FE1",fe1);  CHKERRQ(ierr);

  ierr = PetscFinalize(); CHKERRQ(ierr);

  return 0;

}
