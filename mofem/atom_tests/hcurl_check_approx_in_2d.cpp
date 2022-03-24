/**
 * \file hcurl_check_approx_in_2d
 * \example hcurl_check_approx_in_2d.cpp
 *
 * Approximate polynomial in 2D and check derivatives
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

using FaceEle = MoFEM::FaceElementForcesAndSourcesCore;

using FaceEleOp = FaceEle::UserDataOperator;

constexpr double a0 = 0.0;
constexpr double a1 = 2.0;
constexpr double a2 = -15.0 * a0;
constexpr double a3 = -20.0 / 6 * a1;
constexpr double a4 = 15.0 * a0;
constexpr double a5 = a1;
constexpr double a6 = -a0;

struct ApproxFunctions {
  static FTensor::Tensor1<double, 3> fUn(const double x, const double y) {
    return FTensor::Tensor1<double, 3>(
        6 * a6 * std::pow(x, 5) * std::pow(y, 0) +
            5 * a5 * std::pow(x, 4) * std::pow(y, 1) +
            4 * a4 * std::pow(x, 3) * std::pow(y, 2) +
            3 * a3 * std::pow(x, 2) * std::pow(y, 3) +
            2 * a2 * std::pow(x, 1) * std::pow(y, 4) +
            1 * a1 * std::pow(x, 0) * std::pow(y, 5),

        1 * a5 * std::pow(x, 5) * std::pow(y, 0) +
            2 * a4 * std::pow(x, 4) * std::pow(y, 1) +
            3 * a3 * std::pow(x, 3) * std::pow(y, 2) +
            4 * a2 * std::pow(x, 2) * std::pow(y, 3) +
            5 * a1 * std::pow(x, 1) * std::pow(y, 4) +
            6 * a0 * std::pow(x, 0) * std::pow(y, 5),

        0.);
  }

  static FTensor::Tensor2<double, 3, 2> diffFun(const double x,
                                                const double y) {
    return FTensor::Tensor2<double, 3, 2>(
        30 * a6 * pow(x, 4) * pow(y, 0) + 20 * a5 * pow(x, 3) * pow(y, 1) +
            12 * a4 * pow(x, 2) * pow(y, 2) + 6 * a3 * pow(x, 1) * pow(y, 3) +
            2 * a2 * pow(x, 0) * pow(y, 4),

        5 * a5 * pow(x, 4) * pow(y, 0) + 8 * a4 * pow(x, 3) * pow(y, 1) +
            9 * a3 * pow(x, 2) * pow(y, 2) + 8 * a2 * pow(x, 1) * pow(y, 3) +
            5 * a1 * pow(x, 0) * pow(y, 4),

        5 * a5 * pow(x, 4) * pow(y, 0) + 8 * a4 * pow(x, 3) * pow(y, 1) +
            9 * a3 * pow(x, 2) * pow(y, 2) + 8 * a2 * pow(x, 1) * pow(y, 3) +
            5 * a1 * pow(x, 0) * pow(y, 4),

        2 * a4 * pow(x, 4) * pow(y, 0) + 6 * a3 * pow(x, 3) * pow(y, 1) +
            12 * a2 * pow(x, 2) * pow(y, 2) + 20 * a1 * pow(x, 1) * pow(y, 3) +
            30 * a0 * pow(x, 0) * pow(y, 4),

        0., 0.);
  }
};

struct OpAssembleMat : public FaceEleOp {
  OpAssembleMat() : FaceEleOp("FIELD1", "FIELD1", OPROWCOL) {
    sYmm = false; // FIXME problem is symmetric, should use that
  }

  FTensor::Index<'i', 3> i;

  MatrixDouble nA;
  MoFEMErrorCode doWork(int row_side, int col_side, EntityType row_type,
                        EntityType col_type,
                        EntitiesFieldData::EntData &row_data,
                        EntitiesFieldData::EntData &col_data) {

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
    auto t_base_row = row_data.getFTensor1N<3>();
    for (int gg = 0; gg != nb_gauss_pts; gg++) {
      const double val = getArea() * getGaussPts()(2, gg);
      for (int rr = 0; rr != nb_dofs_row; rr++) {
        auto t_base_col = col_data.getFTensor1N<3>(gg, 0);
        for (int cc = 0; cc != nb_dofs_col; cc++) {
          nA(rr, cc) += val * t_base_row(i) * t_base_col(i);
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
                        EntitiesFieldData::EntData &data) {

    MoFEMFunctionBegin;
    const int nb_dofs = data.getIndices().size();
    if (nb_dofs == 0)
      MoFEMFunctionReturnHot(0);
    const int nb_gauss_pts = data.getN().size1();
    nF.resize(nb_dofs, false);
    nF.clear();
    auto t_base_fun = data.getFTensor1N<3>();
    for (int gg = 0; gg != nb_gauss_pts; gg++) {
      const double val = getArea() * getGaussPts()(2, gg);
      for (int bb = 0; bb != nb_dofs; bb++) {
        const double x = getCoordsAtGaussPts()(gg, 0);
        const double y = getCoordsAtGaussPts()(gg, 1);
        nF(bb) += val * t_base_fun(i) * ApproxFunctions::fUn(x, y)(i);
        ++t_base_fun;
      }
    }
    CHKERR VecSetValues(getFEMethod()->ksp_f, data, &*nF.data().begin(),
                        ADD_VALUES);
    MoFEMFunctionReturn(0);
  }
};

struct OpValsDiffVals : public FaceEleOp {
  MatrixDouble &vAls;
  MatrixDouble &diffVals;
  OpValsDiffVals(MatrixDouble &vals, MatrixDouble &diff_vals)
      : FaceEleOp("FIELD1", OPROW), vAls(vals), diffVals(diff_vals) {}

  FTensor::Index<'i', 3> i;
  FTensor::Index<'j', 2> j;

  MoFEMErrorCode doWork(int side, EntityType type,
                        EntitiesFieldData::EntData &data) {
    MoFEMFunctionBegin;
    const int nb_dofs = data.getIndices().size();
    if (nb_dofs == 0)
      MoFEMFunctionReturnHot(0);
    const int nb_gauss_pts = data.getN().size1();
    if (type == MBEDGE && side == 0) {
      vAls.resize(3, nb_gauss_pts, false);
      diffVals.resize(6, nb_gauss_pts, false);
      vAls.clear();
      diffVals.clear();
    }
    auto t_vals = getFTensor1FromMat<3>(vAls);
    auto t_diff_vals = getFTensor2FromMat<3, 2>(diffVals);
    auto t_base_fun = data.getFTensor1N<3>();
    auto t_diff_base_fun = data.getFTensor2DiffN<3, 2>();
    for (int gg = 0; gg != nb_gauss_pts; gg++) {
      auto t_data = data.getFTensor0FieldData();
      for (int bb = 0; bb != nb_dofs; bb++) {
        t_vals(i) += t_base_fun(i) * t_data;
        t_diff_vals(i, j) += t_diff_base_fun(i, j) * t_data;
        ++t_base_fun;
        ++t_diff_base_fun;
        ++t_data;
      }
      ++t_vals;
      ++t_diff_vals;
    }
    MoFEMFunctionReturn(0);
  }
};

struct OpCheckValsDiffVals : public FaceEleOp {
  MatrixDouble &vAls;
  MatrixDouble &diffVals;
  boost::shared_ptr<MatrixDouble> ptrVals;
  boost::shared_ptr<VectorDouble> ptrDiv;

  OpCheckValsDiffVals(MatrixDouble &vals, MatrixDouble &diff_vals,
                      boost::shared_ptr<MatrixDouble> ptr_vals,
                      boost::shared_ptr<VectorDouble> ptr_div)
      : FaceEleOp("FIELD1", OPROW), vAls(vals), diffVals(diff_vals),
        ptrVals(ptr_vals), ptrDiv(ptr_div) {}

  FTensor::Index<'i', 3> i;
  FTensor::Index<'j', 2> j;

  MoFEMErrorCode doWork(int side, EntityType type,
                        EntitiesFieldData::EntData &data) {
    MoFEMFunctionBegin;
    const double eps = 1e-6;
    if (type == MBEDGE && side == 0) {
      const int nb_gauss_pts = data.getN().size1();

      auto t_vals = getFTensor1FromMat<3>(vAls);
      auto t_diff_vals = getFTensor2FromMat<3, 2>(diffVals);
      auto t_vals_from_op = getFTensor1FromMat<3>(*ptrVals);
      auto t_div_from_op = getFTensor0FromVec(*ptrDiv);

      for (int gg = 0; gg != nb_gauss_pts; gg++) {
        const double x = getCoordsAtGaussPts()(gg, 0);
        const double y = getCoordsAtGaussPts()(gg, 1);

        // Check approximation
        FTensor::Tensor1<double, 3> delta_val;
        delta_val(i) = t_vals(i) - ApproxFunctions::fUn(x, y)(i);
        FTensor::Tensor2<double, 3, 2> delta_diff_val;
        delta_diff_val(i, j) =
            t_diff_vals(i, j) - ApproxFunctions::diffFun(x, y)(i, j);
        double err_val = sqrt(delta_val(i) * delta_val(i));
        double err_diff_val = sqrt(delta_diff_val(i, j) * delta_diff_val(i, j));
        if (err_val > eps)
          SETERRQ1(PETSC_COMM_SELF, MOFEM_ATOM_TEST_INVALID,
                   "Wrong value %4.3e", err_val);

        if (err_diff_val > eps)
          SETERRQ1(PETSC_COMM_SELF, MOFEM_ATOM_TEST_INVALID,
                   "Wrong derivative of value %4.3e", err_diff_val);

        // Check HDiv user data operators
        delta_val(i) = t_vals(i) - t_vals_from_op(i);
        err_val = sqrt(delta_val(i) * delta_val(i));
        if (err_val > eps)
          SETERRQ1(PETSC_COMM_SELF, MOFEM_ATOM_TEST_INVALID,
                   "Wrong value from operator %4.3e", err_val);

        double div = t_diff_vals(0, 0) + t_diff_vals(1, 1);
        double err_div = div - t_div_from_op;
        if (err_div > eps)
          SETERRQ3(PETSC_COMM_SELF, MOFEM_ATOM_TEST_INVALID,
                   "Wrong divergence from operator %4.3e (%4.3e != %4.3e)",
                   err_div, div, t_div_from_op);

        ++t_vals;
        ++t_diff_vals;
        ++t_vals_from_op;
        ++t_div_from_op;
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
    CHKERR simple_interface->loadFile("", "rectangle_tri.h5m");

    // Declare elements
    enum bases { AINSWORTH, DEMKOWICZ, LASBASETOP };
    const char *list_bases[] = {"ainsworth", "demkowicz"};
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

    CHKERR simple_interface->addDomainField("FIELD1", HCURL, base, 1);
    constexpr int order = 5;
    CHKERR simple_interface->setFieldOrder("FIELD1", order);
    CHKERR simple_interface->setUp();
    auto dm = simple_interface->getDM();

    MatrixDouble vals, diff_vals;

    auto assemble_matrices_and_vectors = [&]() {
      MoFEMFunctionBegin;
      auto jac_ptr = boost::make_shared<MatrixDouble>();
      auto inv_jac_ptr = boost::make_shared<MatrixDouble>();
      auto det_ptr = boost::make_shared<VectorDouble>();

      pipeline_mng->getOpDomainRhsPipeline().push_back(
          new OpCalculateHOJacForFace(jac_ptr));
      pipeline_mng->getOpDomainRhsPipeline().push_back(
          new OpInvertMatrix<2>(jac_ptr, det_ptr, inv_jac_ptr));
      pipeline_mng->getOpDomainRhsPipeline().push_back(
          new OpMakeHdivFromHcurl());
      pipeline_mng->getOpDomainRhsPipeline().push_back(
          new OpSetContravariantPiolaTransformOnFace2D(jac_ptr));
      pipeline_mng->getOpDomainRhsPipeline().push_back(new OpAssembleVec());

      pipeline_mng->getOpDomainLhsPipeline().push_back(
          new OpCalculateHOJacForFace(jac_ptr));
      pipeline_mng->getOpDomainLhsPipeline().push_back(
          new OpInvertMatrix<2>(jac_ptr, det_ptr, inv_jac_ptr));
      pipeline_mng->getOpDomainLhsPipeline().push_back(
          new OpMakeHdivFromHcurl());
      pipeline_mng->getOpDomainLhsPipeline().push_back(
          new OpSetContravariantPiolaTransformOnFace2D(jac_ptr));
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

      boost::shared_ptr<MatrixDouble> ptr_values(new MatrixDouble());
      boost::shared_ptr<VectorDouble> ptr_divergence(new VectorDouble());

      pipeline_mng->getOpDomainLhsPipeline().clear();
      pipeline_mng->getOpDomainRhsPipeline().clear();

      auto jac_ptr = boost::make_shared<MatrixDouble>();
      auto inv_jac_ptr = boost::make_shared<MatrixDouble>();
      auto det_ptr = boost::make_shared<VectorDouble>();

      pipeline_mng->getOpDomainRhsPipeline().push_back(
          new OpCalculateHOJacForFace(jac_ptr));
      pipeline_mng->getOpDomainRhsPipeline().push_back(
          new OpInvertMatrix<2>(jac_ptr, det_ptr, inv_jac_ptr));
      pipeline_mng->getOpDomainRhsPipeline().push_back(
          new OpMakeHdivFromHcurl());
      pipeline_mng->getOpDomainRhsPipeline().push_back(
          new OpSetContravariantPiolaTransformOnFace2D(jac_ptr));
      pipeline_mng->getOpDomainRhsPipeline().push_back(
          new OpSetInvJacHcurlFace(inv_jac_ptr));
      pipeline_mng->getOpDomainRhsPipeline().push_back(
          new OpValsDiffVals(vals, diff_vals));
      pipeline_mng->getOpDomainRhsPipeline().push_back(
          new OpCalculateHVecVectorField<3>("FIELD1", ptr_values));
      pipeline_mng->getOpDomainRhsPipeline().push_back(
          new OpCalculateHdivVectorDivergence<3, 2>("FIELD1", ptr_divergence));
      pipeline_mng->getOpDomainRhsPipeline().push_back(
          new OpCheckValsDiffVals(vals, diff_vals, ptr_values, ptr_divergence));

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
