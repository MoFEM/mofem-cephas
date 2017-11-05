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

struct TestFE: public MoFEM::FaceElementForcesAndSourcesCore {
  TestFE(MoFEM::Interface &m_field):
  MoFEM::FaceElementForcesAndSourcesCore(m_field) {
  }
  int getRule(int order) { return 2*order; }
};

struct ApproxFunctions {
  static FTensor::Tensor1<double,3> fUn(const double x,const double y) {
    return
    FTensor::Tensor1<double,3>(
      pow(x,5)+pow(y,5) + pow(x,3)*pow(y,2),
      pow(x,5)+pow(y,5) + pow(x,2)*pow(y,3),
      0
    );
  }
  static FTensor::Tensor2<double,3,2> diffFun(const double x,const double y) {
    return
    FTensor::Tensor2<double,3,2>(
      5*pow(x,4)+3*pow(x,2)*pow(y,2), 5*pow(y,4)+2*pow(x,3)*pow(y,1),
      5*pow(x,4)+2*pow(x,1)*pow(y,3), 5*pow(y,4)+3*pow(x,2)*pow(y,2),
      0,0
    );
  }
};

struct OpAssembleMatAndVec: public MoFEM::FaceElementForcesAndSourcesCore::UserDataOperator {
  Mat A;
  Vec F;
  OpAssembleMatAndVec(Mat a,Vec f):
  MoFEM::FaceElementForcesAndSourcesCore::UserDataOperator("FIELD1","FIELD1",OPROW|OPROWCOL),
  A(a),
  F(f) {
    sYmm = false; // FIXME problem is symmetric, should use that
  }

  FTensor::Index<'i',3> i;

  VectorDouble nF;
  MoFEMErrorCode	doWork(
    int side, EntityType type, DataForcesAndSourcesCore::EntData &data
  ) {

    MoFEMFunctionBeginHot;
    const int nb_dofs = data.getIndices().size();
    if(nb_dofs==0) MoFEMFunctionReturnHot(0);
    const int nb_gauss_pts = data.getN().size1();
    nF.resize(nb_dofs,false);
    nF.clear();
    FTensor::Tensor1<double*,3> t_base_fun = data.getFTensor1HcurlN<3>();
    for(int gg = 0;gg!=nb_gauss_pts;gg++) {
      const double val = getArea()*getGaussPts()(2,gg);
      for(int bb = 0;bb!=nb_dofs;bb++) {
        const double x = getCoordsAtGaussPts()(gg,0);
        const double y = getCoordsAtGaussPts()(gg,1);
        nF(bb) += val*t_base_fun(i)*ApproxFunctions::fUn(x,y)(i);
        ++t_base_fun;
      }
    }
    ierr = VecSetValues(
      F,nb_dofs,&*data.getIndices().data().begin(),&*nF.data().begin(),ADD_VALUES
    ); CHKERRG(ierr);
    MoFEMFunctionReturnHot(0);
  }

  MatrixDouble nA;
  MoFEMErrorCode doWork(
    int row_side,
    int col_side,
    EntityType row_type,
    EntityType col_type,
    DataForcesAndSourcesCore::EntData &row_data,
    DataForcesAndSourcesCore::EntData &col_data
  ) {

    MoFEMFunctionBeginHot;
    const int nb_dofs_row = row_data.getIndices().size();
    if(nb_dofs_row==0) MoFEMFunctionReturnHot(0);
    const int nb_dofs_col = col_data.getIndices().size();
    if(nb_dofs_col==0) MoFEMFunctionReturnHot(0);
    nA.resize(nb_dofs_row,nb_dofs_col,false);
    nA.clear();
    const int nb_gauss_pts = row_data.getN().size1();
    FTensor::Tensor1<double*,3> t_base_row = row_data.getFTensor1HcurlN<3>();
    for(int gg = 0;gg!=nb_gauss_pts;gg++) {
      const double val = getArea()*getGaussPts()(2,gg);
      for(int rr = 0;rr!=nb_dofs_row;rr++) {
        FTensor::Tensor1<double*,3> t_base_col = col_data.getFTensor1HcurlN<3>(gg,0);
        for(int cc = 0;cc!=nb_dofs_col;cc++) {
          nA(rr,cc) += val*t_base_row(i)*t_base_col(i);
          ++t_base_col;
        }
        ++t_base_row;
      }
    }
    ierr = MatSetValues(
      A,
      nb_dofs_row,&*row_data.getIndices().begin(),
      nb_dofs_col,&*col_data.getIndices().begin(),
      &*nA.data().begin(),
      ADD_VALUES
    ); CHKERRG(ierr);
    MoFEMFunctionReturnHot(0);
  }
};

struct OpValsDiffVals: public MoFEM::FaceElementForcesAndSourcesCore::UserDataOperator {
  MatrixDouble &vAls;
  MatrixDouble &diffVals;
  OpValsDiffVals(MatrixDouble &vals,MatrixDouble &diff_vals):
  MoFEM::FaceElementForcesAndSourcesCore::UserDataOperator("FIELD1",OPROW),
  vAls(vals),
  diffVals(diff_vals) {
  }

  FTensor::Index<'i',3> i;
  FTensor::Index<'j',2> j;

  MoFEMErrorCode	doWork(
    int side, EntityType type, DataForcesAndSourcesCore::EntData &data
  ) {
    //
    MoFEMFunctionBeginHot;
    const int nb_dofs = data.getIndices().size();
    if(nb_dofs==0) MoFEMFunctionReturnHot(0);
    const int nb_gauss_pts = data.getN().size1();
    if(type == MBEDGE && side == 0) {
      vAls.resize(3,nb_gauss_pts,false);
      diffVals.resize(6,nb_gauss_pts,false);
      vAls.clear();
      diffVals.clear();
    }
    FTensor::Tensor1<double*,3> t_vals = getTensor1FormData<3>(vAls);
    FTensor::Tensor2<double*,3,2> t_diff_vals = getTensor2FormData<3,2>(diffVals);
    FTensor::Tensor1<double*,3> t_base_fun = data.getFTensor1HcurlN<3>();
    FTensor::Tensor2<double*,3,2> t_diff_base_fun = data.getFTensor2DiffHcurlN<3,2>();
    for(int gg = 0;gg!=nb_gauss_pts;gg++) {
      FTensor::Tensor0<double*> t_data = data.getFTensor0FieldData();
      for(int bb = 0;bb!=nb_dofs;bb++) {
        t_vals(i) += t_base_fun(i)*t_data;
        t_diff_vals(i,j) += t_diff_base_fun(i,j)*t_data;
        ++t_base_fun;
        ++t_diff_base_fun;
        ++t_data;
      }
      ++t_vals;
      ++t_diff_vals;
    }
    MoFEMFunctionReturnHot(0);
  }

};

struct OpCheckValsDiffVals: public MoFEM::FaceElementForcesAndSourcesCore::UserDataOperator {
  MatrixDouble &vAls;
  MatrixDouble &diffVals;
  OpCheckValsDiffVals(MatrixDouble &vals,MatrixDouble &diff_vals):
  MoFEM::FaceElementForcesAndSourcesCore::UserDataOperator("FIELD1",OPROW),
  vAls(vals),
  diffVals(diff_vals) {
  }

  FTensor::Index<'i',3> i;
  FTensor::Index<'j',2> j;

  MoFEMErrorCode	doWork(
    int side, EntityType type, DataForcesAndSourcesCore::EntData &data
  ) {
    //
    MoFEMFunctionBeginHot;
    const double eps = 1e-6;
    if(type == MBEDGE && side == 0) {
      const int nb_gauss_pts = data.getN().size1();
      FTensor::Tensor1<double*,3> t_vals = getTensor1FormData<3>(vAls);
      FTensor::Tensor2<double*,3,2> t_diff_vals = getTensor2FormData<3,2>(diffVals);
      for(int gg = 0;gg!=nb_gauss_pts;gg++) {
        const double x = getCoordsAtGaussPts()(gg,0);
        const double y = getCoordsAtGaussPts()(gg,1);
        FTensor::Tensor1<double,3> delta_val;
        delta_val(i) = t_vals(i)-ApproxFunctions::fUn(x,y)(i);
        FTensor::Tensor2<double,3,2> delta_diff_val;
        delta_diff_val(i,j) = t_diff_vals(i,j)-ApproxFunctions::diffFun(x,y)(i,j);
        double err_val = sqrt(delta_val(i)*delta_val(i));
        double err_diff_val = sqrt(delta_diff_val(i,j)*delta_diff_val(i,j));
        // cerr << "val " << err_val << " ";
        // cerr << "diff " << err_diff_val << endl;
        // cerr << delta_diff_val(0,0) << " " << delta_diff_val(0,1) << endl
        // << delta_diff_val(1,0) << " " << delta_diff_val(1,1) << endl;
        // cerr << sqrt(t_diff_vals(i,j)*t_diff_vals(i,j)) << endl;
        if(err_val>eps) {
          SETERRQ(PETSC_COMM_SELF,MOFEM_ATOM_TEST_INVALID,"Wrong value");
        }
        if(err_diff_val>eps) {
          SETERRQ(PETSC_COMM_SELF,MOFEM_ATOM_TEST_INVALID,"Wrong derivative of value");
        }
        ++t_vals;
        ++t_diff_vals;
      }
    }
    MoFEMFunctionReturnHot(0);
  }

};


int main(int argc, char *argv[]) {




  PetscInitialize(&argc,&argv,(char *)0,help);

  try {

    moab::Core mb_instance;
    moab::Interface& moab = mb_instance;
    int rank;
    MPI_Comm_rank(PETSC_COMM_WORLD,&rank);

    //Read mesh to MOAB
    rval = moab.load_file("rectangle.h5m", 0, ""); CHKERRG(rval);
    ParallelComm* pcomm = ParallelComm::get_pcomm(&moab,MYPCOMM_INDEX);
    if(pcomm == NULL) pcomm =  new ParallelComm(&moab,PETSC_COMM_WORLD);

    //Create MoFEM instance
    MoFEM::Core core(moab);
    MoFEM::Interface& m_field = core;

    //set entitities bit level
    BitRefLevel bit_level0;
    bit_level0.set(0);
    ierr = m_field.getInterface<BitRefManager>()->setBitRefLevelByDim(0,2,bit_level0); CHKERRG(ierr);

    // Declare elements
    ierr = m_field.add_field("FIELD1",HCURL,AINSWORTH_LEGENDRE_BASE,1); CHKERRG(ierr);
    ierr = m_field.add_finite_element("TEST_FE1"); CHKERRG(ierr);
    //Define rows/cols and element data
    ierr = m_field.modify_finite_element_add_field_row("TEST_FE1","FIELD1"); CHKERRG(ierr);
    ierr = m_field.modify_finite_element_add_field_col("TEST_FE1","FIELD1"); CHKERRG(ierr);
    ierr = m_field.modify_finite_element_add_field_data("TEST_FE1","FIELD1"); CHKERRG(ierr);
    //Problem
    ierr = m_field.add_problem("TEST_PROBLEM"); CHKERRG(ierr);
    //set finite elements for problem
    ierr = m_field.modify_problem_add_finite_element("TEST_PROBLEM","TEST_FE1"); CHKERRG(ierr);
    //set refinement level for problem
    ierr = m_field.modify_problem_ref_level_add_bit("TEST_PROBLEM",bit_level0); CHKERRG(ierr);

    // Add entities
    ierr = m_field.add_ents_to_field_by_type(0,MBTRI,"FIELD1"); CHKERRG(ierr);
    // Set order
    int order = 6;
    ierr = m_field.set_field_order(0,MBTRI,"FIELD1",order); CHKERRG(ierr);
    ierr = m_field.set_field_order(0,MBEDGE,"FIELD1",order); CHKERRG(ierr);
    ierr = m_field.add_ents_to_finite_element_by_type(0,MBTRI,"TEST_FE1"); CHKERRG(ierr);

    // Build database
    ierr = m_field.build_fields(); CHKERRG(ierr);
    //build finite elemnts
    ierr = m_field.build_finite_elements(); CHKERRG(ierr);
    //build adjacencies
    ierr = m_field.build_adjacencies(bit_level0); CHKERRG(ierr);

    //build problem
    ProblemsManager *prb_mng_ptr;
    ierr = m_field.getInterface(prb_mng_ptr); CHKERRG(ierr);
    ierr = prb_mng_ptr->buildProblem("TEST_PROBLEM",true); CHKERRG(ierr);
    // Partition
    ierr = prb_mng_ptr->partitionSimpleProblem("TEST_PROBLEM"); CHKERRG(ierr);
    ierr = prb_mng_ptr->partitionFiniteElements("TEST_PROBLEM"); CHKERRG(ierr);
    ierr = prb_mng_ptr->partitionGhostDofs("TEST_PROBLEM"); CHKERRG(ierr);

    // Create matrices
    Mat A;
    ierr = m_field.MatCreateMPIAIJWithArrays("TEST_PROBLEM",&A); CHKERRG(ierr);
    Vec F;
    ierr = m_field.getInterface<VecManager>()->vecCreateGhost("TEST_PROBLEM",ROW,&F); CHKERRG(ierr);
    Vec D;
    ierr = m_field.getInterface<VecManager>()->vecCreateGhost("TEST_PROBLEM",COL,&D); CHKERRG(ierr);

    {
      TestFE fe(m_field);
      fe.getOpPtrVector().push_back(new OpAssembleMatAndVec(A,F));
      ierr = m_field.loop_finite_elements("TEST_PROBLEM","TEST_FE1",fe); CHKERRG(ierr);
      ierr = VecAssemblyBegin(F); CHKERRG(ierr);
      ierr = VecAssemblyEnd(F); CHKERRG(ierr);
      ierr = MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY); CHKERRG(ierr);
      ierr = MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY); CHKERRG(ierr);
    }

    // Solver problem
    KSP solver;
    ierr = KSPCreate(PETSC_COMM_WORLD,&solver); CHKERRG(ierr);
    ierr = KSPSetOperators(solver,A,A); CHKERRG(ierr);
    ierr = KSPSetFromOptions(solver); CHKERRG(ierr);
    ierr = KSPSetUp(solver); CHKERRG(ierr);
    ierr = KSPSolve(solver,F,D); CHKERRG(ierr);
    ierr = VecGhostUpdateBegin(D,INSERT_VALUES,SCATTER_FORWARD); CHKERRG(ierr);
    ierr = VecGhostUpdateEnd(D,INSERT_VALUES,SCATTER_FORWARD); CHKERRG(ierr);
    ierr = m_field.getInterface<VecManager>()->setLocalGhostVector(
      "TEST_PROBLEM",COL,D,INSERT_VALUES,SCATTER_REVERSE
    ); CHKERRG(ierr);
    ierr = KSPDestroy(&solver); CHKERRG(ierr);

    {
      TestFE fe(m_field);
      MatrixDouble inv_jac,vals,diff_vals;
      fe.getOpPtrVector().push_back(new OpCalculateInvJacForFace(inv_jac));
      fe.getOpPtrVector().push_back(new OpSetInvJacHcurlFace(inv_jac));
      fe.getOpPtrVector().push_back(new OpValsDiffVals(vals,diff_vals));
      fe.getOpPtrVector().push_back(new OpCheckValsDiffVals(vals,diff_vals));
      ierr = m_field.loop_finite_elements("TEST_PROBLEM","TEST_FE1",fe); CHKERRG(ierr);
    }

    ierr = VecDestroy(&F); CHKERRG(ierr);
    ierr = VecDestroy(&D); CHKERRG(ierr);
    ierr = MatDestroy(&A); CHKERRG(ierr);

  } catch (MoFEMException const &e) {
    SETERRQ(PETSC_COMM_SELF,e.errorCode,e.errorMessage);
  }

  ierr = PetscFinalize(); CHKERRG(ierr);
}
