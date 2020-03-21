/**
 * \file scalar_check_approximation_2d.cpp
 * \example scalar_check_approximation_2d.cpp
 *
 * Approximate polynomial in 2D and check direvatives
 *
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

using namespace MoFEM;

static char help[] = "...\n\n";

using FaceEle = MoFEM::FaceElementForcesAndSourcesCoreSwitch<
    FaceElementForcesAndSourcesCore::NO_HO_GEOMETRY |
    FaceElementForcesAndSourcesCore::NO_CONTRAVARIANT_TRANSFORM_HDIV |
    FaceElementForcesAndSourcesCore::NO_COVARIANT_TRANSFORM_HCURL>;

using FaceEleOp = FaceEle::UserDataOperator;
using EntData = DataForcesAndSourcesCore::EntData;

struct ApproxFunctions {
  static double fUn(const double x, const double y) {
    return pow(x, 4) + pow(y, 4) + pow(x, 2) * pow(y, 2) + x * y + x + y;
  }

  static FTensor::Tensor1<double, 2> diffFun(const double x, const double y) {
    return FTensor::Tensor1<double, 2>(
        4 * pow(x, 3) + 2 * x * pow(y, 2) + y + 1.,
        4 * pow(y, 3) + 2 * y * pow(x, 2) + x + 1.);
  }
};

struct OpAssembleMat : public FaceEleOp {
  OpAssembleMat() : FaceEleOp("FIELD1", "FIELD1", OPROWCOL) {
    sYmm = false; // FIXME problem is symmetric, should use that
  }

  MatrixDouble nA;
  MoFEMErrorCode doWork(int row_side, int col_side, EntityType row_type,
                        EntityType col_type, EntData &row_data,
                        EntData &col_data) {

    MoFEMFunctionBegin;
    const int nb_dofs_row = row_data.getIndices().size();
    if (nb_dofs_row == 0)
      MoFEMFunctionReturnHot(0);
    const int nb_dofs_col = col_data.getIndices().size();
    if (nb_dofs_col == 0)
      MoFEMFunctionReturnHot(0);
    nA.resize(nb_dofs_row, nb_dofs_col, false);
    nA.clear();
    const int nb_gauss_pts = row_data.getN().size1();
    auto t_base_row = row_data.getFTensor0N();
    for (int gg = 0; gg != nb_gauss_pts; gg++) {
      const double val = getArea() * getGaussPts()(2, gg);
      for (int rr = 0; rr != nb_dofs_row; rr++) {
        auto t_base_col = col_data.getFTensor0N(gg, 0);
        for (int cc = 0; cc != nb_dofs_col; cc++) {
          nA(rr, cc) += val * t_base_row * t_base_col;
          ++t_base_col;
        }
        ++t_base_row;
      }
    }
    CHKERR MatSetValues(getFEMethod()->ksp_B, row_data, col_data,
                        &*nA.data().begin(), ADD_VALUES);
    MoFEMFunctionReturn(0);
  }
};

struct OpAssembleVec : public FaceEleOp {
  OpAssembleVec() : FaceEleOp("FIELD1", "FIELD1", OPROW) {}

  FTensor::Index<'i', 3> i;

  VectorDouble nF;
  MoFEMErrorCode doWork(int side, EntityType type,
                        DataForcesAndSourcesCore::EntData &data) {

    MoFEMFunctionBegin;
    const int nb_dofs = data.getIndices().size();
    if (nb_dofs == 0)
      MoFEMFunctionReturnHot(0);
    const int nb_gauss_pts = data.getN().size1();
    nF.resize(nb_dofs, false);
    nF.clear();
    auto t_base_fun = data.getFTensor0N();
    for (int gg = 0; gg != nb_gauss_pts; gg++) {
      const double val = getArea() * getGaussPts()(2, gg);
      for (int bb = 0; bb != nb_dofs; bb++) {
        const double x = getCoordsAtGaussPts()(gg, 0);
        const double y = getCoordsAtGaussPts()(gg, 1);
        nF(bb) += val * t_base_fun * ApproxFunctions::fUn(x, y);
        ++t_base_fun;
      }
    }
    CHKERR VecSetValues(getFEMethod()->ksp_f, data, &*nF.data().begin(),
                        ADD_VALUES);
    MoFEMFunctionReturn(0);
  }
};

struct OpValsDiffVals : public FaceEleOp {
  VectorDouble &vAls;
  MatrixDouble &diffVals;
  const bool checkGradients;
  OpValsDiffVals(VectorDouble &vals, MatrixDouble &diff_vals, bool check_grads)
      : FaceEleOp("FIELD1", OPROW), vAls(vals), diffVals(diff_vals),
        checkGradients(check_grads) {}

  FTensor::Index<'i', 2> i;

  MoFEMErrorCode doWork(int side, EntityType type,
                        DataForcesAndSourcesCore::EntData &data) {
    MoFEMFunctionBegin;
    const int nb_gauss_pts = getGaussPts().size2();
    if (type == MBVERTEX) {
      vAls.resize(nb_gauss_pts, false);
      diffVals.resize(2, nb_gauss_pts, false);
      vAls.clear();
      diffVals.clear();
    }
    const int nb_dofs = data.getIndices().size();
    if (nb_dofs) {
      auto t_vals = getFTensor0FromVec(vAls);
      auto t_base_fun = data.getFTensor0N();
      for (int gg = 0; gg != nb_gauss_pts; gg++) {
        auto t_data = data.getFTensor0FieldData();
        for (int bb = 0; bb != nb_dofs; bb++) {
          t_vals += t_base_fun * t_data;
          ++t_base_fun;
          ++t_data;
        }
        ++t_vals;
      }

      if (checkGradients) {
        auto t_diff_vals = getFTensor1FromMat<2>(diffVals);
        auto t_diff_base_fun = data.getFTensor1DiffN<2>();
        for (int gg = 0; gg != nb_gauss_pts; gg++) {
          auto t_data = data.getFTensor0FieldData();
          for (int bb = 0; bb != nb_dofs; bb++) {
            t_diff_vals(i) += t_diff_base_fun(i) * t_data;
            ++t_diff_base_fun;
            ++t_data;
          }
          ++t_diff_vals;
        }
      }
    }
    MoFEMFunctionReturn(0);
  }
};

struct OpCheckValsDiffVals : public FaceEleOp {
  VectorDouble &vAls;
  MatrixDouble &diffVals;
  boost::shared_ptr<VectorDouble> ptrVals;
  boost::shared_ptr<MatrixDouble> ptrDiffVals;
  const bool checkGradients;

  OpCheckValsDiffVals(VectorDouble &vals, MatrixDouble &diff_vals,
                      boost::shared_ptr<VectorDouble> &ptr_vals,
                      boost::shared_ptr<MatrixDouble> &ptr_diff_vals,
                      bool check_grads)
      : FaceEleOp("FIELD1", OPROW), vAls(vals), diffVals(diff_vals),
        ptrVals(ptr_vals), ptrDiffVals(ptr_diff_vals),
        checkGradients(check_grads) {
    std::fill(&doEntities[MBEDGE], &doEntities[MBMAXTYPE], false);
  }

  FTensor::Index<'i', 2> i;

  MoFEMErrorCode doWork(int side, EntityType type,
                        DataForcesAndSourcesCore::EntData &data) {
    MoFEMFunctionBegin;
    const double eps = 1e-6;
    const int nb_gauss_pts = data.getN().size1();

    auto t_vals = getFTensor0FromVec(vAls);

    auto t_ptr_vals = getFTensor0FromVec(*ptrVals);

    for (int gg = 0; gg != nb_gauss_pts; gg++) {
      const double x = getCoordsAtGaussPts()(gg, 0);
      const double y = getCoordsAtGaussPts()(gg, 1);

      // Check approximation
      const double delta_val = t_vals - ApproxFunctions::fUn(x, y);

      double err_val = sqrt(delta_val * delta_val);
      if (err_val > eps)
        SETERRQ1(PETSC_COMM_SELF, MOFEM_ATOM_TEST_INVALID, "Wrong value %4.3e",
                 err_val);

      // Check H1 user data operators
      err_val = t_vals - t_ptr_vals;
      if (err_val > eps)
        SETERRQ1(PETSC_COMM_SELF, MOFEM_ATOM_TEST_INVALID,
                 "Wrong value from operator %4.3e", err_val);

      ++t_vals;
      ++t_ptr_vals;
    }

    if (checkGradients) {

      auto t_diff_vals = getFTensor1FromMat<2>(diffVals);
      auto t_ptr_diff_vals = getFTensor1FromMat<2>(*ptrDiffVals);

      for (int gg = 0; gg != nb_gauss_pts; gg++) {
        const double x = getCoordsAtGaussPts()(gg, 0);
        const double y = getCoordsAtGaussPts()(gg, 1);

        // Check approximation
        FTensor::Tensor1<double, 2> t_delta_diff_val;
        t_delta_diff_val(i) =
            t_diff_vals(i) - ApproxFunctions::diffFun(x, y)(i);

        double err_diff_val = sqrt(t_delta_diff_val(i) * t_delta_diff_val(i));
        if (err_diff_val > eps)
          SETERRQ1(PETSC_COMM_SELF, MOFEM_ATOM_TEST_INVALID,
                   "Wrong derivative of value %4.3e", err_diff_val);

        t_delta_diff_val(i) = t_diff_vals(i) - t_ptr_diff_vals(i);
        err_diff_val = sqrt(t_delta_diff_val(i) * t_delta_diff_val(i));
        if (err_diff_val > eps)
          SETERRQ1(PETSC_COMM_SELF, MOFEM_ATOM_TEST_INVALID,
                   "Wrong direvatives from operator %4.3e", err_diff_val);

        ++t_diff_vals;
        ++t_ptr_diff_vals;
      }
    }

    MoFEMFunctionReturn(0);
  }
};

int main(int argc, char *argv[]) {

  MoFEM::Core::Initialize(&argc, &argv, (char *)0, help);

  try {

    DMType dm_name = "DMMOFEM";
    CHKERR DMRegister_MoFEM(dm_name);

    moab::Core mb_instance;
    moab::Interface &moab = mb_instance;

    // Create MoFEM instance
    MoFEM::Core core(moab);
    MoFEM::Interface &m_field = core;

    Simple *simple_interface = m_field.getInterface<Simple>();
    PipelineManager *pipeline_mng = m_field.getInterface<PipelineManager>();
    CHKERR simple_interface->getOptions();
    CHKERR simple_interface->loadFile("", "rectangle.h5m");

    // Declare elements
    enum bases { AINSWORTH, DEMKOWICZ, BERNSTEIN, LASBASETOP };
    const char *list_bases[] = {"ainsworth", "demkowicz", "bernstein"};
    PetscBool flg;
    PetscInt choice_base_value = AINSWORTH;
    CHKERR PetscOptionsGetEList(PETSC_NULL, NULL, "-base", list_bases,
                                LASBASETOP, &choice_base_value, &flg);

    if (flg != PETSC_TRUE)
      SETERRQ(PETSC_COMM_SELF, MOFEM_IMPOSIBLE_CASE, "base not set");
    FieldApproximationBase base = AINSWORTH_LEGENDRE_BASE;
    if (choice_base_value == AINSWORTH)
      base = AINSWORTH_LEGENDRE_BASE;
    else if (choice_base_value == DEMKOWICZ)
      base = DEMKOWICZ_JACOBI_BASE;
    else if (choice_base_value == BERNSTEIN)
      base = AINSWORTH_BERNSTEIN_BEZIER_BASE;

    enum spaces { H1SPACE, L2SPACE, LASBASETSPACE };
    const char *list_spaces[] = {"h1", "l2"};
    PetscInt choice_space_value = H1SPACE;
    CHKERR PetscOptionsGetEList(PETSC_NULL, NULL, "-space", list_spaces,
                                LASBASETSPACE, &choice_space_value, &flg);
    if (flg != PETSC_TRUE)
      SETERRQ(PETSC_COMM_SELF, MOFEM_IMPOSIBLE_CASE, "space not set");
    FieldSpace space = H1;
    if (choice_space_value == H1SPACE)
      space = H1;
    else if (choice_space_value == L2SPACE)
      space = L2;

    CHKERR simple_interface->addDomainField("FIELD1", space, base, 1);
    constexpr int order = 5;
    CHKERR simple_interface->setFieldOrder("FIELD1", order);
    CHKERR simple_interface->setUp();
    auto dm = simple_interface->getDM();

    VectorDouble vals;
    MatrixDouble jac(2, 2), inv_jac(2, 2), diff_vals;

    auto assemble_matrices_and_vectors = [&]() {
      MoFEMFunctionBegin;
      pipeline_mng->getOpDomainRhsPipeline().push_back(new OpAssembleVec());

      pipeline_mng->getOpDomainLhsPipeline().push_back(
          new OpCalculateInvJacForFace(inv_jac));
      pipeline_mng->getOpDomainLhsPipeline().push_back(new OpAssembleMat());

      auto integration_rule = [](int, int, int p_data) { return 2 * p_data; };
      CHKERR pipeline_mng->setDomainRhsIntegrationRule(integration_rule);
      CHKERR pipeline_mng->setDomainLhsIntegrationRule(integration_rule);

      MoFEMFunctionReturn(0);
    };

    auto solve_problem = [&] {
      MoFEMFunctionBegin;
      auto solver = pipeline_mng->createKSP();
      CHKERR KSPSetFromOptions(solver);
      CHKERR KSPSetUp(solver);

      auto D = smartCreateDMVector(dm);
      auto F = smartVectorDuplicate(D);

      CHKERR KSPSolve(solver, F, D);
      CHKERR VecGhostUpdateBegin(D, INSERT_VALUES, SCATTER_FORWARD);
      CHKERR VecGhostUpdateEnd(D, INSERT_VALUES, SCATTER_FORWARD);
      CHKERR DMoFEMMeshToLocalVector(dm, D, INSERT_VALUES, SCATTER_REVERSE);
      MoFEMFunctionReturn(0);
    };

    auto check_solution = [&] {
      MoFEMFunctionBegin;

      boost::shared_ptr<VectorDouble> ptr_values(new VectorDouble());
      boost::shared_ptr<MatrixDouble> ptr_diff_vals(new MatrixDouble());

      pipeline_mng->getDomainLhsFE().reset();
      pipeline_mng->getOpDomainRhsPipeline().clear();

      pipeline_mng->getOpDomainRhsPipeline().push_back(
          new OpCalculateInvJacForFace(inv_jac));
      pipeline_mng->getOpDomainRhsPipeline().push_back(
          new OpSetInvJacH1ForFace(inv_jac));
      pipeline_mng->getOpDomainRhsPipeline().push_back(
          new OpValsDiffVals(vals, diff_vals, choice_space_value == H1SPACE));
      pipeline_mng->getOpDomainRhsPipeline().push_back(
          new OpCalculateScalarFieldValues("FIELD1", ptr_values));
      if (choice_space_value == H1SPACE) {
        pipeline_mng->getOpDomainRhsPipeline().push_back(
            new OpCalculateScalarFieldGradient<2>("FIELD1", ptr_diff_vals));
      }
      pipeline_mng->getOpDomainRhsPipeline().push_back(
          new OpCheckValsDiffVals(vals, diff_vals, ptr_values, ptr_diff_vals,
                                  choice_space_value == H1SPACE));

      CHKERR pipeline_mng->loopFiniteElements();

      MoFEMFunctionReturn(0);
    };

    CHKERR assemble_matrices_and_vectors();
    CHKERR solve_problem();
    CHKERR check_solution();
  }
  CATCH_ERRORS;

  CHKERR MoFEM::Core::Finalize();
}
