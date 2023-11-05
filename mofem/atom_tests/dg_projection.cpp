/**
 * \example dg_projection.cpp
 *
 * Testing DG projection operators
 *
 */

#include <MoFEM.hpp>

using namespace MoFEM;

static char help[] = "...\n\n";

constexpr char FIELD_NAME[] = "U";
constexpr int BASE_DIM = 1;
constexpr int FIELD_DIM = 1;
constexpr int SPACE_DIM = 2;
constexpr int order = 2;

template <int DIM> struct ElementsAndOps {};

template <> struct ElementsAndOps<2> {
  using DomainEle = PipelineManager::FaceEle;
  using DomainEleOp = DomainEle::UserDataOperator;
};

using DomainEle = ElementsAndOps<SPACE_DIM>::DomainEle; ///< Finite elenent type
using DomainEleOp =
    DomainEle::UserDataOperator;            ///< Finire element operator type
using EntData = EntitiesFieldData::EntData; ///< Data on entities

/**
 * @brief Function to approximate
 *
 */
auto fun = [](const double x, const double y, const double z) {
  return x + y + x * x + y * y;
};

/**
 * @brief  OPerator to integrate mass matrix for least square approximation
 *
 */
using OpDomainMass = FormsIntegrators<DomainEleOp>::Assembly<
    PETSC>::BiLinearForm<GAUSS>::OpMass<BASE_DIM, FIELD_DIM>;

/**
 * @brief Operator to integrate the right hand side matrix for the problem
 *
 */
using OpDomainSource = FormsIntegrators<DomainEleOp>::Assembly<
    PETSC>::LinearForm<GAUSS>::OpSource<BASE_DIM, FIELD_DIM>;

struct AtomTest {

  AtomTest(MoFEM::Interface &m_field) : mField(m_field) {}

  MoFEMErrorCode runProblem();

private:
  MoFEM::Interface &mField;
  Simple *simpleInterface;

  MoFEMErrorCode readMesh();
  MoFEMErrorCode setupProblem();
  MoFEMErrorCode assembleSystem();
  MoFEMErrorCode solveSystem();
  MoFEMErrorCode checkResults();

  /**
   * @brief Collected data use d by operator to evaluate errors for the test
   *
   */
  struct CommonData {
    boost::shared_ptr<MatrixDouble> invJacPtr;
    boost::shared_ptr<VectorDouble> approxVals;
    boost::shared_ptr<MatrixDouble> approxGradVals;
    boost::shared_ptr<MatrixDouble> approxHessianVals;
    SmartPetscObj<Vec> L2Vec;
  };

  /**
   * @brief Operator to evaluate errors
   *
   */
  struct OpError;
};

/**
 * @brief Operator to evaluate errors
 *
 */
struct AtomTest::OpError : public DomainEleOp {
  boost::shared_ptr<CommonData> commonDataPtr;
  OpError(boost::shared_ptr<MatrixDouble> data_ptr)
      : DomainEleOp(NOSPACE, OPSPACE), dataPtr(data_ptr) {}

  MoFEMErrorCode doWork(int side, EntityType type, EntData &data) {
    MoFEMFunctionBegin;

    const int nb_integration_pts = getGaussPts().size2();
    auto t_val = getFTensor1FromMat<1>(
        *(dataPtr)); // get function approximation at gauss pts
    auto t_coords = getFTensor1CoordsAtGaussPts(); // get coordinates of
                                                   // integration points

    for (int gg = 0; gg != nb_integration_pts; ++gg) {

      // Calculate errors

      double diff = t_val(0) - fun(t_coords(0), t_coords(1), t_coords(2));
      constexpr double eps = 1e-8;
      if (std::abs(diff) > eps) {
        MOFEM_LOG("SELF", Sev::error) << "Wrong function value " << diff;
        SETERRQ(PETSC_COMM_SELF, MOFEM_ATOM_TEST_INVALID,
                "Wrong function value");
      }

      // move data to next integration point
      ++t_val;
      ++t_coords;
    }

    MOFEM_LOG("SELF", Sev::noisy) << "All is OK";

    MoFEMFunctionReturn(0);
  }

  private:
    boost::shared_ptr<MatrixDouble> dataPtr;
};

//! [Run programme]
MoFEMErrorCode AtomTest::runProblem() {
  MoFEMFunctionBegin;
  CHKERR readMesh();
  CHKERR setupProblem();
  CHKERR assembleSystem();
  CHKERR solveSystem();
  CHKERR checkResults();
  MoFEMFunctionReturn(0);
}
//! [Run programme]

//! [Read mesh]
MoFEMErrorCode AtomTest::readMesh() {
  MoFEMFunctionBegin;

  CHKERR mField.getInterface(simpleInterface);
  CHKERR simpleInterface->getOptions();
  CHKERR simpleInterface->loadFile();

  MoFEMFunctionReturn(0);
}
//! [Read mesh]

//! [Set up problem]
MoFEMErrorCode AtomTest::setupProblem() {
  MoFEMFunctionBegin;
  // Add field
  CHKERR simpleInterface->addDomainField(FIELD_NAME, H1,
                                         AINSWORTH_LEGENDRE_BASE, FIELD_DIM);
  CHKERR simpleInterface->setFieldOrder(FIELD_NAME, order);
  CHKERR simpleInterface->setUp();

  MoFEMFunctionReturn(0);
}
//! [Set up problem]

//! [Push operators to pipeline]
MoFEMErrorCode AtomTest::assembleSystem() {
  MoFEMFunctionBegin;

  auto rule = [](int, int, int p) -> int { return 2 * p; };

  PipelineManager *pipeline_mng = mField.getInterface<PipelineManager>();
  CHKERR pipeline_mng->setDomainLhsIntegrationRule(rule);
  CHKERR pipeline_mng->setDomainRhsIntegrationRule(rule);

  auto beta = [](const double, const double, const double) { return 1; };
  pipeline_mng->getOpDomainLhsPipeline().push_back(
      new OpDomainMass(FIELD_NAME, FIELD_NAME, beta));
  pipeline_mng->getOpDomainRhsPipeline().push_back(
      new OpDomainSource(FIELD_NAME, fun));

  MoFEMFunctionReturn(0);
}
//! [Push operators to pipeline]

//! [Solve]
MoFEMErrorCode AtomTest::solveSystem() {
  MoFEMFunctionBegin;
  PipelineManager *pipeline_mng = mField.getInterface<PipelineManager>();

  MOFEM_LOG("WORLD", Sev::inform) << "Solve problem";

  auto solver = pipeline_mng->createKSP();
  CHKERR KSPSetFromOptions(solver);
  CHKERR KSPSetUp(solver);

  auto dm = simpleInterface->getDM();
  auto D = createDMVector(dm);
  auto F = vectorDuplicate(D);

  CHKERR KSPSolve(solver, F, D);
  CHKERR VecGhostUpdateBegin(D, INSERT_VALUES, SCATTER_FORWARD);
  CHKERR VecGhostUpdateEnd(D, INSERT_VALUES, SCATTER_FORWARD);
  CHKERR DMoFEMMeshToLocalVector(dm, D, INSERT_VALUES, SCATTER_REVERSE);
  MoFEMFunctionReturn(0);
}

//! [Check results]
MoFEMErrorCode AtomTest::checkResults() {
  MoFEMFunctionBegin;
  auto simple = mField.getInterface<Simple>();
  auto pipeline_mng = mField.getInterface<PipelineManager>();
  pipeline_mng->getDomainLhsFE().reset();
  pipeline_mng->getDomainRhsFE().reset();
  pipeline_mng->getOpDomainRhsPipeline().clear();

  auto rule = [](int, int, int p) -> int { return 2 * p + 1; };
  CHKERR pipeline_mng->setDomainRhsIntegrationRule(
      rule); // set integration rule

  auto entity_data_l2 = boost::make_shared<EntitiesFieldData>(
      MBENTITYSET); // entity data shared between
                    // physical and post proc
                    // elements
  auto mass_ptr = boost::make_shared<MatrixDouble>(); // integrated mass matrix
                                                      // of post proc element
  auto coeffs_ptr =
      boost::make_shared<MatrixDouble>(); // vector of coeffs shared between
                                          // physical and post proc elements
  auto data_ptr =
      boost::make_shared<MatrixDouble>(); // data stored at integration points
                                          // of the physical element and
                                          // evaluated at integration points of
                                          // the post proc element

  auto op_this =
      new OpLoopThis<DomainEle>(mField, simple->getDomainFEName(), Sev::noisy);
  pipeline_mng->getOpDomainRhsPipeline().push_back(op_this); // 1
  pipeline_mng->getOpDomainRhsPipeline().push_back(new OpDGProjectionEvaluation(
      data_ptr, coeffs_ptr, entity_data_l2, AINSWORTH_LEGENDRE_BASE,
      L2)); // 5
  pipeline_mng->getOpDomainRhsPipeline().push_back(
      new OpError(data_ptr)); // 6

  auto fe_physics_ptr = op_this->getThisFEPtr();
  fe_physics_ptr->getRuleHook = [](int, int, int p) { return 2 * p; };

  fe_physics_ptr->getOpPtrVector().push_back(new OpDGProjectionMassMatrix(
      order, mass_ptr, entity_data_l2, AINSWORTH_LEGENDRE_BASE, L2)); // 2
  fe_physics_ptr->getOpPtrVector().push_back(
      new OpCalculateVectorFieldValues<FIELD_DIM>(FIELD_NAME,
                                                  data_ptr)); // 3
  fe_physics_ptr->getOpPtrVector().push_back(
      new OpDGProjectionCoefficients(data_ptr, coeffs_ptr, mass_ptr,
                                     entity_data_l2, AINSWORTH_LEGENDRE_BASE,
                                     L2)); // 4

  CHKERR pipeline_mng->loopFiniteElements();

  MoFEMFunctionReturn(0);
}
//! [Check results]

int main(int argc, char *argv[]) {

  // Initialisation of MoFEM/PETSc and MOAB data structures
  MoFEM::Core::Initialize(&argc, &argv, NULL, help);

  try {

    //! [Register MoFEM discrete manager in PETSc]
    DMType dm_name = "DMMOFEM";
    CHKERR DMRegister_MoFEM(dm_name);
    //! [Register MoFEM discrete manager in PETSc

    //! [Create MoAB]
    moab::Core mb_instance;              ///< mesh database
    moab::Interface &moab = mb_instance; ///< mesh database interface
    //! [Create MoAB]

    //! [Create MoFEM]
    MoFEM::Core core(moab);           ///< finite element database
    MoFEM::Interface &m_field = core; ///< finite element database insterface
    //! [Create MoFEM]

    //! [AtomTest]
    AtomTest ex(m_field);
    CHKERR ex.runProblem();
    //! [AtomTest]
  }
  CATCH_ERRORS;

  CHKERR MoFEM::Core::Finalize();
}
