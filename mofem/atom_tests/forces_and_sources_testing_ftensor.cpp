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
using namespace MoFEM;

namespace bio = boost::iostreams;
using bio::tee_device;
using bio::stream;


static char help[] = "...\n\n";

int main(int argc, char *argv[]) {

  
  

  MoFEM::Core::Initialize(&argc,&argv,(char *)0,help);

  try {

  moab::Core mb_instance;
  moab::Interface& moab = mb_instance;
  int rank;
  MPI_Comm_rank(PETSC_COMM_WORLD,&rank);

  PetscBool flg = PETSC_TRUE;
  char mesh_file_name[255];
  #if PETSC_VERSION_GE(3,6,4)
  ierr = PetscOptionsGetString(PETSC_NULL,"","-my_file",mesh_file_name,255,&flg); CHKERRG(ierr);
  #else
  ierr = PetscOptionsGetString(PETSC_NULL,PETSC_NULL,"-my_file",mesh_file_name,255,&flg); CHKERRG(ierr);
  #endif
  if(flg != PETSC_TRUE) {
    SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"*** ERROR -my_file (MESH FILE NEEDED)");
  }

  //Create MoFEM (Joseph) database
  MoFEM::Core core(moab);
  MoFEM::Interface& m_field = core;

  const char *option;
  option = "";
  rval = moab.load_file(mesh_file_name, 0, option); CHKERRG(rval);

  //set entitities bit level
  BitRefLevel bit_level0;
  bit_level0.set(0);
  EntityHandle meshset_level0;
  rval = moab.create_meshset(MESHSET_SET,meshset_level0); CHKERRG(rval);
  ierr = m_field.getInterface<BitRefManager>()->setBitRefLevelByDim(0,3,bit_level0); CHKERRG(ierr);

  //Fields
  ierr = m_field.add_field("FIELD1",H1,AINSWORTH_LEGENDRE_BASE,3); CHKERRG(ierr);
  ierr = m_field.add_field("FIELD2",H1,AINSWORTH_LEGENDRE_BASE,1); CHKERRG(ierr);


  //FE
  ierr = m_field.add_finite_element("TEST_FE1"); CHKERRG(ierr);

  //Define rows/cols and element data
  ierr = m_field.modify_finite_element_add_field_row("TEST_FE1","FIELD1"); CHKERRG(ierr);
  ierr = m_field.modify_finite_element_add_field_col("TEST_FE1","FIELD1"); CHKERRG(ierr);
  ierr = m_field.modify_finite_element_add_field_data("TEST_FE1","FIELD1"); CHKERRG(ierr);
  ierr = m_field.modify_finite_element_add_field_data("TEST_FE1","FIELD2"); CHKERRG(ierr);

  //Problem
  ierr = m_field.add_problem("TEST_PROBLEM"); CHKERRG(ierr);

  //set finite elements for problem
  ierr = m_field.modify_problem_add_finite_element("TEST_PROBLEM","TEST_FE1"); CHKERRG(ierr);
  //set refinement level for problem
  ierr = m_field.modify_problem_ref_level_add_bit("TEST_PROBLEM",bit_level0); CHKERRG(ierr);


  //meshset consisting all entities in mesh
  EntityHandle root_set = moab.get_root_set();
  //add entities to field
  ierr = m_field.add_ents_to_field_by_type(root_set,MBTET,"FIELD1"); CHKERRG(ierr);
  ierr = m_field.add_ents_to_field_by_type(root_set,MBTET,"FIELD2"); CHKERRG(ierr);


  //add entities to finite element
  ierr = m_field.add_ents_to_finite_element_by_type(root_set,MBTET,"TEST_FE1"); CHKERRG(ierr);

  //set app. order
  //see Hierarchic Finite Element Bases on Unstructured Tetrahedral Meshes (Mark Ainsworth & Joe Coyle)
  int order = 4;
  ierr = m_field.set_field_order(root_set,MBTET,"FIELD1",order); CHKERRG(ierr);
  ierr = m_field.set_field_order(root_set,MBTRI,"FIELD1",order); CHKERRG(ierr);
  ierr = m_field.set_field_order(root_set,MBEDGE,"FIELD1",order); CHKERRG(ierr);
  ierr = m_field.set_field_order(root_set,MBVERTEX,"FIELD1",1); CHKERRG(ierr);

  ierr = m_field.set_field_order(root_set,MBTET,"FIELD2",order); CHKERRG(ierr);
  ierr = m_field.set_field_order(root_set,MBTRI,"FIELD2",order); CHKERRG(ierr);
  ierr = m_field.set_field_order(root_set,MBEDGE,"FIELD2",order); CHKERRG(ierr);
  ierr = m_field.set_field_order(root_set,MBVERTEX,"FIELD2",1); CHKERRG(ierr);

  /****/
  //build database
  //build field
  ierr = m_field.build_fields(); CHKERRG(ierr);
  //build finite elemnts
  ierr = m_field.build_finite_elements(); CHKERRG(ierr);
  //build adjacencies
  ierr = m_field.build_adjacencies(bit_level0); CHKERRG(ierr);


  ProblemsManager *prb_mng_ptr;
  ierr = m_field.getInterface(prb_mng_ptr); CHKERRG(ierr);
  //build problem
  ierr = prb_mng_ptr->buildProblem("TEST_PROBLEM",true); CHKERRG(ierr);
  ierr = prb_mng_ptr->partitionSimpleProblem("TEST_PROBLEM"); CHKERRG(ierr);
  ierr = prb_mng_ptr->partitionFiniteElements("TEST_PROBLEM"); CHKERRG(ierr);
  ierr = prb_mng_ptr->partitionGhostDofs("TEST_PROBLEM"); CHKERRG(ierr);

  typedef tee_device<std::ostream, std::ofstream> TeeDevice;
  typedef stream<TeeDevice> TeeStream;

  std::ofstream ofs("forces_and_sources_testing_ftensor.txt");
  TeeDevice my_tee(std::cout, ofs);
  TeeStream my_split(my_tee);

  struct MyOp1: public VolumeElementForcesAndSourcesCore::UserDataOperator {

    boost::shared_ptr<MatrixDouble> field1ValuesDataPtr;
    boost::shared_ptr<VectorDouble> field2ValuesDataPtr;
    boost::shared_ptr<MatrixDouble> grad1ValuesDataPtr;
    boost::shared_ptr<MatrixDouble> grad2ValuesDataPtr;
    TeeStream &my_split;

    MyOp1(
      boost::shared_ptr<MatrixDouble> field1_values_data_ptr,
      boost::shared_ptr<VectorDouble> field2_values_data_ptr,
      boost::shared_ptr<MatrixDouble> grad1_values_data_ptr,
      boost::shared_ptr<MatrixDouble> grad2_values_data_ptr,
      TeeStream &_my_split,
      char type
    ):
    VolumeElementForcesAndSourcesCore::UserDataOperator("FIELD1","FIELD2",type),
    field1ValuesDataPtr(field1_values_data_ptr),
    field2ValuesDataPtr(field2_values_data_ptr),
    grad1ValuesDataPtr(grad1_values_data_ptr),
    grad2ValuesDataPtr(grad2_values_data_ptr),
    my_split(_my_split) {}


    MoFEMErrorCode doWork(
      int side,
      EntityType type,
      DataForcesAndSourcesCore::EntData &data) {
      MoFEMFunctionBeginHot;

      FTensor::Tensor2<double, 3, 3> t2;
      const int nb_gauss_pts = data.getN().size1();
      const int nb_base_functions = data.getN().size2();

      FTensor::Index<'I',3> I;
      FTensor::Index<'J',3> J;

      auto base_function = data.getFTensor0N();
      auto diff_base = data.getFTensor1DiffN<3>();
      // FTensor::Tensor1<double*,3> field_values = getFTensor1FromMat<3>(field1ValuesDataPtr);
      // #ifdef WITH_ADOL_C
      // adouble val;
      // FTensor::Tensor0<adouble*> adouble_t0(&val);
      // // cerr << endl;
      // // cerr << base_function << endl;
      // adouble_t0<<=base_function;
      // adouble_t0>>=base_function;
      // // cerr << adouble_t0 << " " << base_function << endl;
      // if(adouble_t0 - base_function != 0) {
      //   SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"Data inconsistency");
      // }
      // FTensor::Tensor1<adouble,3> adouble_t1;
      // adouble_t1(I)<<=diff_base(I);
      // adouble_t1(I)>>=diff_base(I);
      // for(int II = 0;II!=3;II++) {
      //   // cerr << adouble_t1(II) << " " << diff_base(II) << endl;
      //   if(adouble_t1(II)-diff_base(II)!=0) {
      //     SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"Data inconsistency");
      //   }
      // }
      // FTensor::Tensor2<adouble,3,3> adouble_t2;
      // FTensor::Tensor2<double,3,3> double_t2;
      // double_t2(I,J) = diff_base(I)*diff_base(J);
      // adouble_t2(I,J)<<=double_t2(I,J);
      // adouble_t2(I,J)>>=double_t2(I,J);
      // cerr << endl;
      // for(int II = 0;II!=3;II++) {
      //   for(int JJ = 0;JJ!=3;JJ++) {
      //     cerr << adouble_t2(II,JJ) << " " << double_t2(II,JJ) << endl;
      //     if(adouble_t2(II,JJ)-double_t2(II,JJ)!=0) {
      //       SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"Data inconsistency");
      //     }
      //   }
      // }
      // adouble_t2(0,0)<<=1;
      // adouble_t1(0)<<=1;
      // double a;
      // adouble_t2(0,0)>>=a;
      // adouble_t1(0)>>=a;
      // #endif

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

          t2(I,J) = diff_base(I)*diff_base(J);

          MatrixAdaptor mat = MatrixAdaptor(3,3,ublas::shallow_array_adaptor<double>(9,&t2(0,0)));

          for(int II = 0;II!=3;II++) {
            for(int JJ = 0;JJ!=3;JJ++) {
              if(t2(II,JJ)-mat(II,JJ)!=0) {
                SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"Data inconsistency");
              }
            }
          }


          ++diff_base;

        }

        // VectorAdaptor vec = VectorAdaptor(3,ublas::shallow_array_adaptor<double>(3,&field_values(0)));
        // my_split << vec << endl;

        // ++field_values;

      }
      MoFEMFunctionReturnHot(0);
    }

    MoFEMErrorCode doWork(
      int row_side,int col_side,
      EntityType row_type,EntityType col_type,
      DataForcesAndSourcesCore::EntData &row_data,
      DataForcesAndSourcesCore::EntData &col_data
    ) {
      MoFEMFunctionBeginHot;
      //

      const int nb_gauss_pts = row_data.getN().size1();
      const int nb_base_functions_row = row_data.getN().size2();
      const int nb_base_functions_col = col_data.getN().size2();

      FTensor::Number<0> N0;
      FTensor::Number<1> N1;
      FTensor::Number<2> N2;

      FTensor::Index<'I',3> I;
      FTensor::Index<'J',3> J;
      FTensor::Index<'K',3> K;

      FTensor::Tensor2<double,3,3> t2;
      FTensor::Tensor3<double,3,3,3> t3;

      for(int br = 0;br!=nb_base_functions_row;br++) {
        auto base_row = row_data.getFTensor0N(br);
        auto diff_base_row = row_data.getFTensor1DiffN<3>(br);

        for(int bc = 0;bc!=nb_base_functions_col;bc++) {
          auto base_col = col_data.getFTensor0N(bc);
          auto diff_base_col= row_data.getFTensor1DiffN<3>(bc);

          auto field1_values = getFTensor1FromMat<3>(*field1ValuesDataPtr);
          auto field2_values = getFTensor0FromVec(*field2ValuesDataPtr);
          auto grad1_values = getFTensor2FromMat<3,3>(*grad1ValuesDataPtr);
          auto grad2_values = getFTensor1FromMat<3>(*grad2ValuesDataPtr);

          for(int gg = 0;gg!=nb_gauss_pts;gg++) {
            // This make no sense (just do some calculations)
            // FIXME: Some stuff can be calculated only in loop by integration pts for efficiency (I don't care for purpose of this test)
            t2(I,J) = diff_base_row(I)*diff_base_col(J)*base_row*base_col;
            t3(I,J,K) = grad1_values(I,J)*field1_values(K)*field2_values;
            ++field1_values;
            ++field2_values;
            ++grad1_values;
            ++grad2_values;
          }
          ++base_col;
          ++diff_base_col;

        }
        ++base_row;
        ++diff_base_row;

      }

      MoFEMFunctionReturnHot(0);
    }

  };

  boost::shared_ptr<MatrixDouble> values1_at_gauss_pts_ptr =
      boost::shared_ptr<MatrixDouble>(new MatrixDouble());
  boost::shared_ptr<VectorDouble> values2_at_gauss_pts_ptr =
      boost::shared_ptr<VectorDouble>(new VectorDouble());
  boost::shared_ptr<MatrixDouble> grad1_at_gauss_pts_ptr =
      boost::shared_ptr<MatrixDouble>(new MatrixDouble());
  boost::shared_ptr<MatrixDouble> grad2_at_gauss_pts_ptr =
      boost::shared_ptr<MatrixDouble>(new MatrixDouble());

  VolumeElementForcesAndSourcesCore fe1(m_field);

  fe1.getOpPtrVector().push_back(
      new OpCalculateVectorFieldValues<3>("FIELD1", values1_at_gauss_pts_ptr));
  fe1.getOpPtrVector().push_back(
      new OpCalculateScalarFieldValues("FIELD2", values2_at_gauss_pts_ptr));
  fe1.getOpPtrVector().push_back(new OpCalculateVectorFieldGradient<3, 3>(
      "FIELD1", grad1_at_gauss_pts_ptr));
  fe1.getOpPtrVector().push_back(
      new OpCalculateScalarFieldGradient<3>("FIELD2", grad2_at_gauss_pts_ptr));

  fe1.getOpPtrVector().push_back(
    new MyOp1(
      values1_at_gauss_pts_ptr,
      values2_at_gauss_pts_ptr,
      grad1_at_gauss_pts_ptr,
      grad2_at_gauss_pts_ptr,
      my_split,
      ForcesAndSourcesCore::UserDataOperator::OPROW
    )
  );
  fe1.getOpPtrVector().push_back(
    new MyOp1(
      values1_at_gauss_pts_ptr,
      values2_at_gauss_pts_ptr,
      grad1_at_gauss_pts_ptr,
      grad2_at_gauss_pts_ptr,
      my_split,
      ForcesAndSourcesCore::UserDataOperator::OPROWCOL
    )
  );

  ierr = m_field.loop_finite_elements("TEST_PROBLEM","TEST_FE1",fe1);  CHKERRG(ierr);


  } catch (MoFEMException const &e) {
    SETERRQ(PETSC_COMM_SELF,e.errorCode,e.errorMessage);
  }

  ierr = MoFEM::Core::Finalize(); CHKERRG(ierr);

  return 0;

}
