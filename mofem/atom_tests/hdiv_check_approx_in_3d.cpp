/**
 * \file hdiv_check_approx_in_3d
 * \example hdiv_check_approx_in_3d.cpp
 *
 * Approximate vectorial polynomial in 3D and check derivatives
 *
 */

#include <MoFEM.hpp>

using namespace MoFEM;

static char help[] = "...\n\n";

constexpr int BASE_DIM = 3;
constexpr int SPACE_DIM = 3;

using Ele = MoFEM::PipelineManager::ElementsAndOpsByDim<SPACE_DIM>::DomainEle;
using EleOp = Ele::UserDataOperator;

/**
 * @brief  OPerator to integrate mass matrix for least square approximation
 *
 */
using OpDomainMass = FormsIntegrators<EleOp>::Assembly<PETSC>::BiLinearForm<
    GAUSS>::OpMass<BASE_DIM, SPACE_DIM>;

/**
 * @brief Operator to integrate the right hand side matrix for the problem
 *
 */
using OpDomainSource = FormsIntegrators<EleOp>::Assembly<PETSC>::LinearForm<
    GAUSS>::OpSource<BASE_DIM, SPACE_DIM>;

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
        x + y*y + x*x*x,
        // y
        y + z*x + y*y*y,
        // z
        z + x*x + z*z*z);
  }

  static FTensor::Tensor2<double, BASE_DIM, SPACE_DIM>
  diffFun(const double x, const double y, const double z) {
    return FTensor::Tensor2<double, BASE_DIM, SPACE_DIM>(
        // x,x
        1. + 3 * x * x,
        // x,y
        2. * y,
        // x, z
        0.,
        // y,x
        z,
        // y,y
        1. + 3 * y * y,
        // y,z
        x,
        // z
        2 * x, 0., 1. + 3 * z * z);
  }

}; // namespace ApproxFunctions

struct OpCheckValsDiffVals : public EleOp {
  boost::shared_ptr<MatrixDouble> ptrVals;
  boost::shared_ptr<MatrixDouble> ptrGrad;

  OpCheckValsDiffVals(boost::shared_ptr<MatrixDouble> ptr_vals,
                      boost::shared_ptr<MatrixDouble> ptr_grad)
      : EleOp(NOSPACE, OPSPACE), ptrVals(ptr_vals), ptrGrad(ptr_grad) {}

  FTENSOR_INDEX(SPACE_DIM, i);
  FTENSOR_INDEX(SPACE_DIM, j);
  FTENSOR_INDEX(SPACE_DIM, k);

  MoFEMErrorCode doWork(int side, EntityType type,
                        EntitiesFieldData::EntData &data) {
    MoFEMFunctionBegin;

    const double eps = 1e-6;
    const int nb_gauss_pts = getGaussPts().size2();

    auto t_vals_from_op = getFTensor1FromMat<SPACE_DIM>(*ptrVals);
    auto t_grad_from_op = getFTensor2FromMat<SPACE_DIM, SPACE_DIM>(*ptrGrad);

    for (int gg = 0; gg != nb_gauss_pts; gg++) {
      const double x = getCoordsAtGaussPts()(gg, 0);
      const double y = getCoordsAtGaussPts()(gg, 1);
      const double z = getCoordsAtGaussPts()(gg, 2);

      // Check approximation
      FTensor::Tensor1<double, SPACE_DIM> delta_val;
      delta_val(i) = t_vals_from_op(i) - ApproxFunctions::fUn(x, y, z)(i);
      double err_val = sqrt(delta_val(i) * delta_val(i));
      MOFEM_LOG("SELF", Sev::verbose) << "Approximation error: " << err_val;

      if (err_val > eps)
        SETERRQ1(PETSC_COMM_SELF, MOFEM_ATOM_TEST_INVALID, "Wrong value %4.3e",
                 err_val);


      FTensor::Tensor2<double, SPACE_DIM, SPACE_DIM> delta_diff_val;
      delta_diff_val(i, j) =
          t_grad_from_op(i, j) - ApproxFunctions::diffFun(x, y, z)(i, j);
      double err_diff_val = sqrt(delta_diff_val(i, j) * delta_diff_val(i, j));
      if (err_diff_val > eps)
        SETERRQ1(PETSC_COMM_SELF, MOFEM_ATOM_TEST_INVALID,
                 "Wrong derivative of value %4.3e", err_diff_val);

      ++t_vals_from_op;
      ++t_grad_from_op;
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

    Simple *simple = m_field.getInterface<Simple>();
    PipelineManager *pipeline_mng = m_field.getInterface<PipelineManager>();
    CHKERR simple->getOptions();
    CHKERR simple->loadFile();

    // Declare elements
    enum bases { AINSWORTH, DEMKOWICZ, LASBASETOP };
    const char *list_bases[] = {"ainsworth", "demkowicz"};
    PetscBool flg;
    PetscInt choice_base_value = AINSWORTH;
    CHKERR PetscOptionsGetEList(PETSC_NULL, NULL, "-base", list_bases,
                                LASBASETOP, &choice_base_value, &flg);

    if (flg != PETSC_TRUE)
      SETERRQ(PETSC_COMM_SELF, MOFEM_IMPOSSIBLE_CASE, "base not set");
    FieldApproximationBase base = AINSWORTH_LEGENDRE_BASE;
    if (choice_base_value == AINSWORTH)
      base = AINSWORTH_LEGENDRE_BASE;
    else if (choice_base_value == DEMKOWICZ)
      base = DEMKOWICZ_JACOBI_BASE;

    AinsworthOrderHooks::broken_nbfacetri_edge_hdiv = [](int p) { return p; };
    AinsworthOrderHooks::broken_nbfacetri_face_hdiv = [](int p) { return p; };
    AinsworthOrderHooks::broken_nbvolumetet_edge_hdiv = [](int p) { return p; };
    AinsworthOrderHooks::broken_nbvolumetet_face_hdiv = [](int p) { return p; };
    AinsworthOrderHooks::broken_nbvolumetet_volume_hdiv = [](int p) {
      return p;
    };

    CHKERR simple->addDomainField("FIELD1", HDIV, base, 1);
    constexpr int order = 4;
    CHKERR simple->setFieldOrder("FIELD1", order);
    CHKERR simple->setUp();
    auto dm = simple->getDM();

    MatrixDouble vals, diff_vals;

    auto assemble_matrices_and_vectors = [&]() {
      MoFEMFunctionBegin;

      CHKERR AddHOOps<SPACE_DIM, SPACE_DIM, SPACE_DIM>::add(
          pipeline_mng->getOpDomainRhsPipeline(), {HDIV});
      pipeline_mng->getOpDomainRhsPipeline().push_back(
          new OpDomainSource("FIELD1", ApproxFunctions::fUn));

      CHKERR AddHOOps<SPACE_DIM, SPACE_DIM, SPACE_DIM>::add(
          pipeline_mng->getOpDomainLhsPipeline(), {HDIV});
      pipeline_mng->getOpDomainLhsPipeline().push_back(

          new OpDomainMass("FIELD1", "FIELD1",
                           [](double x, double, double) { return 1; })

      );

      auto integration_rule = [](int, int, int p_data) {
        return 2 * p_data + 1;
      };
      CHKERR pipeline_mng->setDomainRhsIntegrationRule(integration_rule);
      CHKERR pipeline_mng->setDomainLhsIntegrationRule(integration_rule);

      MoFEMFunctionReturn(0);
    };

    auto solve_problem = [&] {
      MoFEMFunctionBegin;
      auto solver = pipeline_mng->createKSP();
      CHKERR KSPSetFromOptions(solver);
      CHKERR KSPSetUp(solver);

      auto D = createDMVector(dm);
      auto F = vectorDuplicate(D);

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
      auto ptr_grad = boost::make_shared<MatrixDouble>();

      // Change H-curl to H-div in 2D, and apply Piola transform
      CHKERR AddHOOps<SPACE_DIM, SPACE_DIM, SPACE_DIM>::add(
          pipeline_mng->getOpDomainRhsPipeline(), {HDIV});

      // Calculate field values at integration points
      pipeline_mng->getOpDomainRhsPipeline().push_back(
          new OpCalculateHVecVectorField<BASE_DIM>("FIELD1", ptr_values));
      // Gradient
      pipeline_mng->getOpDomainRhsPipeline().push_back(
          new OpCalculateHVecVectorGradient<BASE_DIM, SPACE_DIM>("FIELD1",
                                                                 ptr_grad));
      pipeline_mng->getOpDomainRhsPipeline().push_back(
          new OpCheckValsDiffVals(ptr_values, ptr_grad));

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
