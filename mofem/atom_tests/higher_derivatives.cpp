/**
 * \example higher_derivatives.cpp
 *
 * Testing higher derivatives of base functions
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

constexpr char FIELD_NAME[] = "U";
constexpr int BASE_DIM = 1;
constexpr int FIELD_DIM = 1;
constexpr int SPACE_DIM = 2;

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
  return x * x + y * y + x * y + pow(x, 3) + pow(y, 3) + pow(x, 4) + pow(y, 4);
};

/**
 * @brief Function derivative
 *
 */
auto diff_fun = [](const double x, const double y, const double z) {
  return FTensor::Tensor1<double, SPACE_DIM>{
      2 * x + y + 3 * pow(x, 2) + 4 * pow(x, 3),
      2 * y + x + 3 * pow(y, 2) + 4 * pow(y, 3)};
};

/**
 * @brief Function second derivative
 *
 */
auto diff2_fun = [](const double x, const double y, const double z) {
  return FTensor::Tensor2<double, SPACE_DIM, SPACE_DIM>{
      2 + 6 * x + 12 * pow(x, 2), 1.,

      1., 2 + 6 * y + 12 * pow(y, 2)};
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
  OpError(boost::shared_ptr<CommonData> &common_data_ptr)
      : DomainEleOp(FIELD_NAME, OPROW), commonDataPtr(common_data_ptr) {
    std::fill(doEntities.begin(), doEntities.end(), false);
    doEntities[MBVERTEX] = true;
  }
  MoFEMErrorCode doWork(int side, EntityType type, EntData &data) {
    MoFEMFunctionBegin;

    const int nb_integration_pts = getGaussPts().size2();
    auto t_w = getFTensor0IntegrationWeight(); // ger integration weights
    auto t_val = getFTensor0FromVec(*(
        commonDataPtr->approxVals)); // get function approximation at gauss pts
    auto t_grad_val = getFTensor1FromMat<SPACE_DIM>(
        *(commonDataPtr
              ->approxGradVals)); // get gradient of approximation at gauss pts
    auto t_hessian_val = getFTensor2FromMat<SPACE_DIM, SPACE_DIM>(
        *(commonDataPtr)->approxHessianVals); // get hessian of approximation
                                              // at integration pts

    auto t_inv_jac = getFTensor2FromMat<SPACE_DIM, SPACE_DIM>(
        *(commonDataPtr->invJacPtr)); // get inverse of element jacobian
    auto t_coords = getFTensor1CoordsAtGaussPts(); // get coordinates of
                                                   // integration points

    // Indices used for tensor operations
    FTensor::Index<'i', 2> i;
    FTensor::Index<'j', 2> j;
    FTensor::Index<'k', 2> k;
    FTensor::Index<'l', 2> l;

    const double volume = getMeasure(); // get finite element area

    auto t_row_base = data.getFTensor0N();
    auto t_diff_row_base = data.getFTensor1DiffN<2>();

    std::array<double, 3> error = {0, 0,
                                   0}; // array for storing operator errors

    for (int gg = 0; gg != nb_integration_pts; ++gg) {

      const double alpha = t_w * volume;

      // Calculate errors

      double diff = t_val - fun(t_coords(0), t_coords(1), t_coords(2));
      error[0] += alpha * pow(diff, 2);
      auto t_diff_grad = diff_fun(t_coords(0), t_coords(1), t_coords(2));
      t_diff_grad(i) -= t_grad_val(j) * t_inv_jac(j, i);
      error[1] += alpha * t_diff_grad(i) *
                  t_diff_grad(i); // note push forward derivatives

      FTensor::Tensor2<double, SPACE_DIM, SPACE_DIM> t_hessian_push;
      t_hessian_push(i, j) =
          (t_hessian_val(k, l) * t_inv_jac(k, i)) * t_inv_jac(l, j);

      MOFEM_LOG("SELF", Sev::noisy) << "t_hessian_val " << t_hessian_push;

      // hessian expected to have symmetry
      if (std::abs(t_hessian_push(0, 1) - t_hessian_push(1, 0)) >
          std::numeric_limits<float>::epsilon()) {
        MOFEM_LOG("SELF", Sev::error) << "t_hessian_push " << t_hessian_push;
        MOFEM_LOG("SELF", Sev::error) << "t_hessian_val " << t_hessian_val;
        SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
                "Hessian should be symmetric");
      }

      auto t_diff_hessian = diff2_fun(t_coords(0), t_coords(1), t_coords(2));
      t_diff_hessian(i, j) -= t_hessian_push(i, j);
      error[2] = t_diff_hessian(i, j) * t_diff_hessian(i, j);

      // move data to next integration point
      ++t_w;
      ++t_val;
      ++t_grad_val;
      ++t_hessian_val;
      ++t_inv_jac;
      ++t_coords;
    }

    // assemble error vector
    std::array<int, 3> index = {0, 1, 2};
    CHKERR VecSetValues(commonDataPtr->L2Vec, 3, index.data(), error.data(),
                        ADD_VALUES);

    MoFEMFunctionReturn(0);
  }
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
  constexpr int order = 4;
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
  auto D = smartCreateDMVector(dm);
  auto F = smartVectorDuplicate(D);

  CHKERR KSPSolve(solver, F, D);
  CHKERR VecGhostUpdateBegin(D, INSERT_VALUES, SCATTER_FORWARD);
  CHKERR VecGhostUpdateEnd(D, INSERT_VALUES, SCATTER_FORWARD);
  CHKERR DMoFEMMeshToLocalVector(dm, D, INSERT_VALUES, SCATTER_REVERSE);
  MoFEMFunctionReturn(0);
}

//! [Check results]
MoFEMErrorCode AtomTest::checkResults() {
  MoFEMFunctionBegin;
  PipelineManager *pipeline_mng = mField.getInterface<PipelineManager>();
  pipeline_mng->getDomainLhsFE().reset();
  pipeline_mng->getDomainRhsFE().reset();
  pipeline_mng->getOpDomainRhsPipeline().clear();

  auto rule = [](int, int, int p) -> int { return 2 * p; };
  CHKERR pipeline_mng->setDomainRhsIntegrationRule(
      rule); // set integration rule

  // create data structures for operator
  auto common_data_ptr = boost::make_shared<CommonData>();
  common_data_ptr->L2Vec = createSmartVectorMPI(
      mField.get_comm(), (!mField.get_comm_rank()) ? 3 : 0, 3);
  common_data_ptr->approxVals = boost::make_shared<VectorDouble>();
  common_data_ptr->approxGradVals = boost::make_shared<MatrixDouble>();
  common_data_ptr->approxHessianVals = boost::make_shared<MatrixDouble>();

  // create data strutires for evaluation of higher order derivatives
  auto base_mass = boost::make_shared<MatrixDouble>();
  auto data_l2 = boost::make_shared<EntitiesFieldData>(MBENTITYSET);
  auto jac_ptr = boost::make_shared<MatrixDouble>();
  auto det_ptr = boost::make_shared<VectorDouble>();
  auto inv_jac_ptr = boost::make_shared<MatrixDouble>();
  common_data_ptr->invJacPtr = inv_jac_ptr;

  // calculate jacobian at integration points
  pipeline_mng->getOpDomainRhsPipeline().push_back(
      new OpCalculateHOJacForFace(jac_ptr));
  // calculate incerse of jacobian at integration points
  pipeline_mng->getOpDomainRhsPipeline().push_back(
      new OpInvertMatrix<SPACE_DIM>(jac_ptr, det_ptr, inv_jac_ptr));

  // calculate value of function at integration points
  pipeline_mng->getOpDomainRhsPipeline().push_back(
      new OpCalculateScalarFieldValues(FIELD_NAME,
                                       common_data_ptr->approxVals));

  // calculate gradient at integration points
  pipeline_mng->getOpDomainRhsPipeline().push_back(
      new OpCalculateScalarFieldGradient<SPACE_DIM>(
          FIELD_NAME, common_data_ptr->approxGradVals));

  // calculate mass matrix to project derivatives
  pipeline_mng->getOpDomainRhsPipeline().push_back(
      new OpBaseDerivativesMass<BASE_DIM>(base_mass, data_l2,
                                          AINSWORTH_LEGENDRE_BASE, L2));
  // calculate second derivative of base functions, i.e. hessian
  pipeline_mng->getOpDomainRhsPipeline().push_back(
      new OpBaseDerivativesNext<BASE_DIM>(2, base_mass, data_l2,
                                          AINSWORTH_LEGENDRE_BASE, H1));
  // calculate third derivative
  pipeline_mng->getOpDomainRhsPipeline().push_back(
      new OpBaseDerivativesNext<BASE_DIM>(3, base_mass, data_l2,
                                          AINSWORTH_LEGENDRE_BASE, H1));
  // calculate forth derivative
  pipeline_mng->getOpDomainRhsPipeline().push_back(
      new OpBaseDerivativesNext<BASE_DIM>(4, base_mass, data_l2,
                                          AINSWORTH_LEGENDRE_BASE, H1));

  // calculate hessian at integration points
  pipeline_mng->getOpDomainRhsPipeline().push_back(
      new OpCalculateScalarFieldHessian<SPACE_DIM>(
          FIELD_NAME, common_data_ptr->approxHessianVals));

  //FIXME: Note third and forth derivative is calculated but not properly tested

  // push operator to evaluate errors
  pipeline_mng->getOpDomainRhsPipeline().push_back(
      new OpError(common_data_ptr));

  // zero error vector and iterate over all elements on the mesh
  CHKERR VecZeroEntries(common_data_ptr->L2Vec);
  CHKERR pipeline_mng->loopFiniteElements();

  // assemble error vector
  CHKERR VecAssemblyBegin(common_data_ptr->L2Vec);
  CHKERR VecAssemblyEnd(common_data_ptr->L2Vec);

  // check if errors are small and print results

  constexpr double eps = 1e-8;

  const double *array;
  CHKERR VecGetArrayRead(common_data_ptr->L2Vec, &array);
  MOFEM_LOG_C("WORLD", Sev::inform,
              "Error %6.4e Diff Error %6.4e Diff2 Error %6.4e\n",
              std::sqrt(array[0]), std::sqrt(array[1]), std::sqrt(array[2]));
  if (std::sqrt(array[0]) > eps)
    SETERRQ(PETSC_COMM_SELF, MOFEM_ATOM_TEST_INVALID, "Wrong function value");
  if (std::sqrt(array[1]) > eps)
    SETERRQ(PETSC_COMM_SELF, MOFEM_ATOM_TEST_INVALID,
            "Wrong function first direcative");
  if (std::sqrt(array[2]) > eps)
    SETERRQ(PETSC_COMM_SELF, MOFEM_ATOM_TEST_INVALID,
            "Wrong function second direcative");

  CHKERR VecRestoreArrayRead(common_data_ptr->L2Vec, &array);

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
