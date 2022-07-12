/**
 * \file hcurl_check_approx_in_2d
 * \example hcurl_check_approx_in_2d.cpp
 *
 * Approximate polynomial in 2D and check derivatives
 *
 */

/* MIT License
 *
 * Copyright (c) 2022
 * 
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 * 
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 * 
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */

#include <MoFEM.hpp>

using namespace MoFEM;

static char help[] = "...\n\n";

using FaceEle = MoFEM::FaceElementForcesAndSourcesCore;
using FaceEleOp = FaceEle::UserDataOperator;

constexpr int BASE_DIM = 3;
constexpr int SPACE_DIM = 2;

/**
 * @brief  OPerator to integrate mass matrix for least square approximation
 *
 */
using OpDomainMass = FormsIntegrators<FaceEleOp>::Assembly<PETSC>::BiLinearForm<
    GAUSS>::OpMass<BASE_DIM, BASE_DIM>;

/**
 * @brief Operator to integrate the right hand side matrix for the problem
 *
 */
using OpDomainSource = FormsIntegrators<FaceEleOp>::Assembly<PETSC>::LinearForm<
    GAUSS>::OpSource<BASE_DIM, BASE_DIM>;

constexpr double a0 = 0.0;
constexpr double a1 = 2.0;
constexpr double a2 = -15.0 * a0;
constexpr double a3 = -20.0 / 6 * a1;
constexpr double a4 = 15.0 * a0;
constexpr double a5 = a1;
constexpr double a6 = -a0;
struct ApproxFunctions {

  static FTensor::Tensor1<double, BASE_DIM> fUn(const double x, const double y,
                                                double z) {
    return FTensor::Tensor1<double, BASE_DIM>(
        // x
        6 * a6 * std::pow(x, 5) * std::pow(y, 0) +
            5 * a5 * std::pow(x, 4) * std::pow(y, 1) +
            4 * a4 * std::pow(x, 3) * std::pow(y, 2) +
            3 * a3 * std::pow(x, 2) * std::pow(y, 3) +
            2 * a2 * std::pow(x, 1) * std::pow(y, 4) +
            1 * a1 * std::pow(x, 0) * std::pow(y, 5),
        // y
        1 * a5 * std::pow(x, 5) * std::pow(y, 0) +
            2 * a4 * std::pow(x, 4) * std::pow(y, 1) +
            3 * a3 * std::pow(x, 3) * std::pow(y, 2) +
            4 * a2 * std::pow(x, 2) * std::pow(y, 3) +
            5 * a1 * std::pow(x, 1) * std::pow(y, 4) +
            6 * a0 * std::pow(x, 0) * std::pow(y, 5),

        // z
        0.);
  }

  static FTensor::Tensor2<double, BASE_DIM, SPACE_DIM> diffFun(const double x,
                                                               const double y) {
    return FTensor::Tensor2<double, BASE_DIM, SPACE_DIM>(
        // x,x
        30 * a6 * pow(x, 4) * pow(y, 0) + 20 * a5 * pow(x, 3) * pow(y, 1) +
            12 * a4 * pow(x, 2) * pow(y, 2) + 6 * a3 * pow(x, 1) * pow(y, 3) +
            2 * a2 * pow(x, 0) * pow(y, 4),
        // x,y
        5 * a5 * pow(x, 4) * pow(y, 0) + 8 * a4 * pow(x, 3) * pow(y, 1) +
            9 * a3 * pow(x, 2) * pow(y, 2) + 8 * a2 * pow(x, 1) * pow(y, 3) +
            5 * a1 * pow(x, 0) * pow(y, 4),
        // y,x
        5 * a5 * pow(x, 4) * pow(y, 0) + 8 * a4 * pow(x, 3) * pow(y, 1) +
            9 * a3 * pow(x, 2) * pow(y, 2) + 8 * a2 * pow(x, 1) * pow(y, 3) +
            5 * a1 * pow(x, 0) * pow(y, 4),
        // y,y
        2 * a4 * pow(x, 4) * pow(y, 0) + 6 * a3 * pow(x, 3) * pow(y, 1) +
            12 * a2 * pow(x, 2) * pow(y, 2) + 20 * a1 * pow(x, 1) * pow(y, 3) +
            30 * a0 * pow(x, 0) * pow(y, 4),
        // z
        0., 0.);
  }

  static FTensor::Tensor3<double, BASE_DIM, SPACE_DIM, SPACE_DIM>
  diff2Fun(const double x, const double y) {
    return FTensor::Tensor3<double, BASE_DIM, SPACE_DIM, SPACE_DIM>(
        // x,xx 0/000

        30 * 4 * a6 * pow(x, 3) * pow(y, 0) +
            20 * 3 * a5 * pow(x, 2) * pow(y, 1) +
            12 * 2 * a4 * pow(x, 1) * pow(y, 2) +
            6 * 1 * a3 * pow(x, 0) * pow(y, 3),

        // x,xy 1/001

        20 * 1 * a5 * pow(x, 3) * pow(y, 0) +
            12 * 2 * a4 * pow(x, 2) * pow(y, 2) +
            6 * 3 * a3 * pow(x, 1) * pow(y, 2) +
            2 * 4 * a2 * pow(x, 0) * pow(y, 3),

        // x,yx 2/010

        5 * 4 * a5 * pow(x, 3) * pow(y, 0) +
            8 * 3 * a4 * pow(x, 2) * pow(y, 1) +
            9 * 2 * a3 * pow(x, 1) * pow(y, 2) +
            8 * 1 * a2 * pow(x, 0) * pow(y, 3),

        // x,yy 3/011

        8 * 1 * a4 * pow(x, 3) * pow(y, 0) +
            9 * 2 * a3 * pow(x, 2) * pow(y, 1) +
            8 * 3 * a2 * pow(x, 1) * pow(y, 2) +
            5 * 4 * a1 * pow(x, 0) * pow(y, 3),

        // y,xx 4/100

        5 * 4 * a5 * pow(x, 3) * pow(y, 0) +
            8 * 3 * a4 * pow(x, 2) * pow(y, 1) +
            9 * 2 * a3 * pow(x, 1) * pow(y, 2) +
            8 * 1 * a2 * pow(x, 0) * pow(y, 3),

        // y,xy 5/101

        8 * 1 * a4 * pow(x, 3) * pow(y, 0) +
            9 * 2 * a3 * pow(x, 2) * pow(y, 1) +
            8 * 3 * a2 * pow(x, 1) * pow(y, 2) +
            5 * 4 * a1 * pow(x, 0) * pow(y, 3),

        // y,yx 6/110

        2 * 4 * a4 * pow(x, 3) * pow(y, 0) +
            6 * 3 * a3 * pow(x, 2) * pow(y, 1) +
            12 * 2 * a2 * pow(x, 1) * pow(y, 2) +
            20 * 1 * a1 * pow(x, 0) * pow(y, 3),

        // y,yy 7/111

        6 * 1 * a3 * pow(x, 3) * pow(y, 0) +
            12 * 2 * a2 * pow(x, 2) * pow(y, 1) +
            20 * 3 * a1 * pow(x, 1) * pow(y, 2) +
            30 * 4 * a0 * pow(x, 0) * pow(y, 3),

        // z,xx 8/200
        0.,
        // z,xy 9/201
        0.,
        // z,yx 10/210
        0.,
        // z,yy 11/211
        0.);
  }
};

struct OpCheckValsDiffVals : public FaceEleOp {
  boost::shared_ptr<MatrixDouble> ptrVals;
  boost::shared_ptr<VectorDouble> ptrDiv;
  boost::shared_ptr<MatrixDouble> ptrGrad;
  boost::shared_ptr<MatrixDouble> ptrHess;

  OpCheckValsDiffVals(boost::shared_ptr<MatrixDouble> ptr_vals,
                      boost::shared_ptr<VectorDouble> ptr_div,
                      boost::shared_ptr<MatrixDouble> ptr_grad,
                      boost::shared_ptr<MatrixDouble> ptr_hess)
      : FaceEleOp("FIELD1", OPROW), ptrVals(ptr_vals), ptrDiv(ptr_div),
        ptrGrad(ptr_grad), ptrHess(ptr_hess) {}

  FTensor::Index<'i', 3> i;
  FTensor::Index<'j', 2> j;
  FTensor::Index<'k', 2> k;

  MoFEMErrorCode doWork(int side, EntityType type,
                        EntitiesFieldData::EntData &data) {
    MoFEMFunctionBegin;
    const double eps = 1e-6;
    if (type == MBEDGE && side == 0) {
      const int nb_gauss_pts = data.getN().size1();

      auto t_vals_from_op = getFTensor1FromMat<3>(*ptrVals);
      auto t_div_from_op = getFTensor0FromVec(*ptrDiv);
      auto t_grad_from_op = getFTensor2FromMat<3, 2>(*ptrGrad);
      auto t_hess_from_op = getFTensor3FromMat<3, 2, 2>(*ptrHess);

      for (int gg = 0; gg != nb_gauss_pts; gg++) {
        const double x = getCoordsAtGaussPts()(gg, 0);
        const double y = getCoordsAtGaussPts()(gg, 1);

        // Check approximation
        FTensor::Tensor1<double, 3> delta_val;
        delta_val(i) = t_vals_from_op(i) - ApproxFunctions::fUn(x, y, 0)(i);
        double err_val = sqrt(delta_val(i) * delta_val(i));
        if (err_val > eps)
          SETERRQ1(PETSC_COMM_SELF, MOFEM_ATOM_TEST_INVALID,
                   "Wrong value %4.3e", err_val);

        FTensor::Tensor2<double, 3, 2> delta_diff_val;
        delta_diff_val(i, j) =
            t_grad_from_op(i, j) - ApproxFunctions::diffFun(x, y)(i, j);
        double err_diff_val = sqrt(delta_diff_val(i, j) * delta_diff_val(i, j));
        if (err_diff_val > eps)
          SETERRQ1(PETSC_COMM_SELF, MOFEM_ATOM_TEST_INVALID,
                   "Wrong derivative of value %4.3e", err_diff_val);

        double div = t_grad_from_op(0, 0) + t_grad_from_op(1, 1);
        double err_div = div - t_div_from_op;
        if (err_div > eps)
          SETERRQ3(PETSC_COMM_SELF, MOFEM_ATOM_TEST_INVALID,
                   "Wrong divergence from operator %4.3e (%4.3e != %4.3e)",
                   err_div, div, t_div_from_op);

        FTensor::Tensor3<double, 3, 2, 2> delta_diff2_val;
        delta_diff2_val(i, j, k) =
            t_hess_from_op(i, j, k) - ApproxFunctions::diff2Fun(x, y)(i, j, k);
        double hess_diff_error =
            sqrt(delta_diff2_val(i, j, k) * delta_diff2_val(i, j, k));
        if (hess_diff_error > eps)
          SETERRQ1(PETSC_COMM_SELF, MOFEM_ATOM_TEST_INVALID,
                   "Wrong hessian from operator %4.3e", hess_diff_error);

        ++t_vals_from_op;
        ++t_div_from_op;
        ++t_grad_from_op;
        ++t_hess_from_op;
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
          new OpCalculateHOJac<SPACE_DIM>(jac_ptr));
      pipeline_mng->getOpDomainRhsPipeline().push_back(
          new OpInvertMatrix<SPACE_DIM>(jac_ptr, det_ptr, inv_jac_ptr));
      pipeline_mng->getOpDomainRhsPipeline().push_back(
          new OpMakeHdivFromHcurl());
      pipeline_mng->getOpDomainRhsPipeline().push_back(
          new OpSetContravariantPiolaTransformOnFace2D(jac_ptr));
      pipeline_mng->getOpDomainRhsPipeline().push_back(
          new OpDomainSource("FIELD1", ApproxFunctions::fUn));

      pipeline_mng->getOpDomainLhsPipeline().push_back(
          new OpCalculateHOJac<SPACE_DIM>(jac_ptr));
      pipeline_mng->getOpDomainLhsPipeline().push_back(
          new OpInvertMatrix<SPACE_DIM>(jac_ptr, det_ptr, inv_jac_ptr));
      pipeline_mng->getOpDomainLhsPipeline().push_back(
          new OpMakeHdivFromHcurl());
      pipeline_mng->getOpDomainLhsPipeline().push_back(
          new OpSetContravariantPiolaTransformOnFace2D(jac_ptr));
      pipeline_mng->getOpDomainLhsPipeline().push_back(

          new OpDomainMass("FIELD1", "FIELD1",
                           [](double, double, double) { return 1; })

      );

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

      pipeline_mng->getOpDomainLhsPipeline().clear();
      pipeline_mng->getOpDomainRhsPipeline().clear();

      auto ptr_values = boost::make_shared<MatrixDouble>();
      auto ptr_divergence = boost::make_shared<VectorDouble>();
      auto ptr_grad = boost::make_shared<MatrixDouble>();
      auto ptr_hessian = boost::make_shared<MatrixDouble>();

      auto jac_ptr = boost::make_shared<MatrixDouble>();
      auto inv_jac_ptr = boost::make_shared<MatrixDouble>();
      auto det_ptr = boost::make_shared<VectorDouble>();

      // Change H-curl to H-div in 2D, and apply Piola transform
      pipeline_mng->getOpDomainRhsPipeline().push_back(
          new OpCalculateHOJac<SPACE_DIM>(jac_ptr));
      pipeline_mng->getOpDomainRhsPipeline().push_back(
          new OpInvertMatrix<SPACE_DIM>(jac_ptr, det_ptr, inv_jac_ptr));
      pipeline_mng->getOpDomainRhsPipeline().push_back(
          new OpMakeHdivFromHcurl());
      pipeline_mng->getOpDomainRhsPipeline().push_back(
          new OpSetContravariantPiolaTransformOnFace2D(jac_ptr));
      pipeline_mng->getOpDomainRhsPipeline().push_back(
          new OpSetInvJacHcurlFace(inv_jac_ptr));

      // Evaluate base function second derivative
      auto base_mass = boost::make_shared<MatrixDouble>();
      auto data_l2 = boost::make_shared<EntitiesFieldData>(MBENTITYSET);

      pipeline_mng->getOpDomainRhsPipeline().push_back(
          new OpBaseDerivativesMass<BASE_DIM>(base_mass, data_l2, base, L2));
      pipeline_mng->getOpDomainRhsPipeline().push_back(
          new OpBaseDerivativesSetHOInvJacobian<SPACE_DIM>(data_l2,
                                                           inv_jac_ptr));
      pipeline_mng->getOpDomainRhsPipeline().push_back(
          new OpBaseDerivativesNext<BASE_DIM>(BaseDerivatives::SecondDerivative,
                                              base_mass, data_l2, base, HCURL));

      // Calculate field values at integration points
      pipeline_mng->getOpDomainRhsPipeline().push_back(
          new OpCalculateHVecVectorField<BASE_DIM>("FIELD1", ptr_values));
      // Gradient
      pipeline_mng->getOpDomainRhsPipeline().push_back(
          new OpCalculateHVecVectorGradient<BASE_DIM, SPACE_DIM>("FIELD1",
                                                                 ptr_grad));
      // Hessian
      pipeline_mng->getOpDomainRhsPipeline().push_back(
          new OpCalculateHdivVectorDivergence<BASE_DIM, SPACE_DIM>(
              "FIELD1", ptr_divergence)); 
      pipeline_mng->getOpDomainRhsPipeline().push_back(
          new OpCalculateHVecVectorHessian<BASE_DIM, SPACE_DIM>("FIELD1",
                                                                ptr_hessian));

      pipeline_mng->getOpDomainRhsPipeline().push_back(new OpCheckValsDiffVals(
          ptr_values, ptr_divergence, ptr_grad, ptr_hessian));

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
