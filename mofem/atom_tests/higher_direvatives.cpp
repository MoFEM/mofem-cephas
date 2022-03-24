/**
 * \example higer_direvatives.cpp
 *
 * Testing higher direvatives of base functions
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
constexpr int FIELD_DIM = 1;
constexpr int SPACE_DIM = 2;

template <int DIM> struct ElementsAndOps {};

template <> struct ElementsAndOps<2> {
  using DomainEle = PipelineManager::FaceEle;
  using DomainEleOp = DomainEle::UserDataOperator;
  using DomianParentEle = FaceElementForcesAndSourcesCoreOnChildParentSwitch<0>;
};

template <> struct ElementsAndOps<3> {
  using DomainEle = VolumeElementForcesAndSourcesCore;
  using DomainEleOp = DomainEle::UserDataOperator;
};

using DomainEle = ElementsAndOps<SPACE_DIM>::DomainEle;
using DomainParentEle = ElementsAndOps<SPACE_DIM>::DomianParentEle;
using DomainEleOp = DomainEle::UserDataOperator;
using EntData = EntitiesFieldData::EntData;

auto fun = [](const double x, const double y, const double z) {
  return x * x + y * y + x * y + pow(x, 3) + pow(y, 3) + pow(x, 4) + pow(y, 4);
};

auto diff_fun = [](const double x, const double y, const double z) {
  return FTensor::Tensor1<double, 2>{2 * x + y + 3 * pow(x, 2) + 4 * pow(x, 3),
                                     2 * y + x + 3 * pow(y, 2) + 4 * pow(y, 3)};
};

auto diff2_fun = [](const double x, const double y, const double z) {
  return FTensor::Tensor2<double, 2, 2>{2 + 6 * x + 12 * pow(x, 2), 1.,

                                        1., 2 + 6 * y + 12 * pow(y, 2)};
};

using OpDomainMass = FormsIntegrators<DomainEleOp>::Assembly<
    PETSC>::BiLinearForm<GAUSS>::OpMass<1, FIELD_DIM>;
using OpDomainSource = FormsIntegrators<DomainEleOp>::Assembly<
    PETSC>::LinearForm<GAUSS>::OpSource<1, FIELD_DIM>;

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

  struct CommonData {
    boost::shared_ptr<MatrixDouble> invJacPtr;
    boost::shared_ptr<VectorDouble> approxVals;
    boost::shared_ptr<MatrixDouble> approxGradVals;
    boost::shared_ptr<MatrixDouble> approxHessianVals;
    SmartPetscObj<Vec> L2Vec;
    SmartPetscObj<Vec> resVec;
  };

  struct OpError;
};

/**
 * @brief Operator to evaluate errors
 * 
 */
struct AtomTest::OpError : public DomainEleOp {
  boost::shared_ptr<CommonData> commonDataPtr;
  OpError(boost::shared_ptr<CommonData> &common_data_ptr)
      : DomainEleOp(FIELD_NAME, OPROW), commonDataPtr(common_data_ptr) {}
  MoFEMErrorCode doWork(int side, EntityType type, EntData &data) {
    MoFEMFunctionBegin;

    if (const size_t nb_dofs = data.getIndices().size()) {

      const int nb_integration_pts = getGaussPts().size2();
      auto t_w = getFTensor0IntegrationWeight();
      auto t_val = getFTensor0FromVec(*(commonDataPtr->approxVals));
      auto t_grad_val = getFTensor1FromMat<2>(*(commonDataPtr->approxGradVals));
      auto t_hessian_val =
          getFTensor2FromMat<2, 2>(*(commonDataPtr)->approxHessianVals);
      auto t_inv_jac = getFTensor2FromMat<2, 2>(*(commonDataPtr->invJacPtr));
      auto t_coords = getFTensor1CoordsAtGaussPts();

      VectorDouble nf(nb_dofs, false);
      nf.clear();

      FTensor::Index<'i', 2> i;
      FTensor::Index<'j', 2> j;
      FTensor::Index<'k', 2> k;
      FTensor::Index<'l', 2> l;
      const double volume = getMeasure();

      auto t_row_base = data.getFTensor0N();
      auto t_diff_row_base = data.getFTensor1DiffN<2>();

      std::array<double, 3> error = {0, 0, 0};
      for (int gg = 0; gg != nb_integration_pts; ++gg) {

        const double alpha = t_w * volume;

        double diff = t_val - fun(t_coords(0), t_coords(1), t_coords(2));
        error[0] += alpha * pow(diff, 2);
        auto t_diff_grad = diff_fun(t_coords(0), t_coords(1), t_coords(2));
        t_diff_grad(i) -= t_grad_val(j) * t_inv_jac(j, i);
        error[1] += alpha * t_diff_grad(i) * t_diff_grad(i);

        FTensor::Tensor2<double, 2, 2> t_hessian_push;
        t_hessian_push(i, j) =
            (t_hessian_val(k, l) * t_inv_jac(k, i)) * t_inv_jac(l, j);

        MOFEM_LOG("SELF", Sev::noisy) << "t_hessian_val " << t_hessian_push;

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

        for (size_t r = 0; r != nb_dofs; ++r) {
          nf[r] += alpha * t_row_base * diff;
          ++t_row_base;
          ++t_diff_row_base;
        }

        ++t_w;
        ++t_val;
        ++t_grad_val;
        ++t_hessian_val;
        ++t_inv_jac;
        ++t_coords;
      }

      std::array<int, 3> index = {0, 1, 2};
      CHKERR VecSetValues(commonDataPtr->L2Vec, 3, index.data(), error.data(),
                          ADD_VALUES);
      CHKERR VecSetValues(commonDataPtr->resVec, data, &nf[0], ADD_VALUES);
    }

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

  auto rule = [](int, int, int p) -> int { return 2 * p + 1; };
  CHKERR pipeline_mng->setDomainRhsIntegrationRule(rule);

  auto common_data_ptr = boost::make_shared<CommonData>();
  common_data_ptr->resVec = smartCreateDMVector(simpleInterface->getDM());
  common_data_ptr->L2Vec = createSmartVectorMPI(
      mField.get_comm(), (!mField.get_comm_rank()) ? 3 : 0, 3);
  common_data_ptr->approxVals = boost::make_shared<VectorDouble>();
  common_data_ptr->approxGradVals = boost::make_shared<MatrixDouble>();
  common_data_ptr->approxHessianVals = boost::make_shared<MatrixDouble>();

  auto base_mass = boost::make_shared<MatrixDouble>();
  auto data_l2 = boost::make_shared<EntitiesFieldData>(MBENTITYSET);
  auto jac_ptr = boost::make_shared<MatrixDouble>();
  auto det_ptr = boost::make_shared<VectorDouble>();
  auto inv_jac_ptr = boost::make_shared<MatrixDouble>();
  common_data_ptr->invJacPtr = inv_jac_ptr;

  pipeline_mng->getOpDomainRhsPipeline().push_back(
      new OpCalculateHOJacForFace(jac_ptr));
  pipeline_mng->getOpDomainRhsPipeline().push_back(
      new OpInvertMatrix<2>(jac_ptr, det_ptr, inv_jac_ptr));

  pipeline_mng->getOpDomainRhsPipeline().push_back(
      new OpCalculateScalarFieldValues(FIELD_NAME,
                                       common_data_ptr->approxVals));
  pipeline_mng->getOpDomainRhsPipeline().push_back(
      new OpCalculateScalarFieldGradient<2>(FIELD_NAME,
                                            common_data_ptr->approxGradVals));

  // calculate mass matrix to project direvatives
  pipeline_mng->getOpDomainRhsPipeline().push_back(new OpBaseDerivativesMass<1>(
      base_mass, data_l2, AINSWORTH_LEGENDRE_BASE, L2));
  // calculate second direvative of base functions, i.e. hessian
  pipeline_mng->getOpDomainRhsPipeline().push_back(new OpBaseDerivativesNext<1>(
      2, base_mass, data_l2, AINSWORTH_LEGENDRE_BASE, H1));
  // calculate third direvative
  pipeline_mng->getOpDomainRhsPipeline().push_back(new OpBaseDerivativesNext<1>(
      3, base_mass, data_l2, AINSWORTH_LEGENDRE_BASE, H1));
  // calculate forth direvative
  pipeline_mng->getOpDomainRhsPipeline().push_back(new OpBaseDerivativesNext<1>(
      4, base_mass, data_l2, AINSWORTH_LEGENDRE_BASE, H1));

  pipeline_mng->getOpDomainRhsPipeline().push_back(
      new OpCalculateScalarFieldHessian<2>(FIELD_NAME,
                                           common_data_ptr->approxHessianVals));

  pipeline_mng->getOpDomainRhsPipeline().push_back(
      new OpError(common_data_ptr));

  CHKERR VecZeroEntries(common_data_ptr->L2Vec);
  CHKERR VecZeroEntries(common_data_ptr->resVec);

  CHKERR pipeline_mng->loopFiniteElements();

  CHKERR VecAssemblyBegin(common_data_ptr->L2Vec);
  CHKERR VecAssemblyEnd(common_data_ptr->L2Vec);
  CHKERR VecAssemblyBegin(common_data_ptr->resVec);
  CHKERR VecAssemblyEnd(common_data_ptr->resVec);

  constexpr double eps = 1e-8;

  double nrm2;
  CHKERR VecNorm(common_data_ptr->resVec, NORM_2, &nrm2);
  const double *array;
  CHKERR VecGetArrayRead(common_data_ptr->L2Vec, &array);
  MOFEM_LOG_C("WORLD", Sev::inform,
              "Error %6.4e Diff Error %6.4e Diff2 Error %6.4e Vec norm %6.4e\n",
              std::sqrt(array[0]), std::sqrt(array[1]), std::sqrt(array[2]),
              nrm2);
  if (std::sqrt(array[0]) > eps)
    SETERRQ(PETSC_COMM_SELF, MOFEM_ATOM_TEST_INVALID, "Wrong function value");
  if (std::sqrt(array[1]) > eps)
    SETERRQ(PETSC_COMM_SELF, MOFEM_ATOM_TEST_INVALID,
            "Wrong function first direcative");
  if (std::sqrt(array[2]) > eps)
    SETERRQ(PETSC_COMM_SELF, MOFEM_ATOM_TEST_INVALID,
            "Wrong function second direcative");

  CHKERR VecRestoreArrayRead(common_data_ptr->L2Vec, &array);

  if (nrm2 > eps)
    SETERRQ1(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
             "Not converged solution err = %6.4e", nrm2);
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
