/**
 * \file child_and_parent.cpp
 * \example child_and_parent.cpp
 *
 * Testing projection child and parent.
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
};
 
template <> struct ElementsAndOps<3> {
  using DomainEle = VolumeElementForcesAndSourcesCore;
  using DomainEleOp = DomainEle::UserDataOperator;
};
 
using DomainEle = ElementsAndOps<SPACE_DIM>::DomainEle;
using DomainEleOp = DomainEle::UserDataOperator;
using EntData = DataForcesAndSourcesCore::EntData;
 
template <int FIELD_DIM> struct ApproxFieldFunction;
 
template <> struct ApproxFieldFunction<1> {
  double operator()(const double x, const double y, const double z) {
    return sin(x * 10.) * cos(y * 10.);
  }
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
 
  static ApproxFieldFunction<FIELD_DIM> approxFunction;
 
  MoFEMErrorCode readMesh();
  MoFEMErrorCode setupProblem();
  MoFEMErrorCode setIntegrationRules();
  MoFEMErrorCode createCommonData();
  MoFEMErrorCode boundaryCondition();
  MoFEMErrorCode assembleSystem();
  MoFEMErrorCode solveSystem();
  MoFEMErrorCode outputResults();
  MoFEMErrorCode checkResults();
 
  struct CommonData {
    boost::shared_ptr<VectorDouble> approxVals;
    SmartPetscObj<Vec> L2Vec;
    SmartPetscObj<Vec> resVec;
  };
  boost::shared_ptr<CommonData> commonDataPtr;
 
  template <int FIELD_DIM> struct OpError;
};
 
ApproxFieldFunction<FIELD_DIM> AtomTest::approxFunction =
    ApproxFieldFunction<FIELD_DIM>();
 
template <> struct AtomTest::OpError<1> : public DomainEleOp {
  boost::shared_ptr<CommonData> commonDataPtr;
  OpError(boost::shared_ptr<CommonData> &common_data_ptr)
      : DomainEleOp(FIELD_NAME, OPROW), commonDataPtr(common_data_ptr) {}
  MoFEMErrorCode doWork(int side, EntityType type, EntData &data) {
    MoFEMFunctionBegin;
 
    if (const size_t nb_dofs = data.getIndices().size()) {
 
      const int nb_integration_pts = getGaussPts().size2();
      auto t_w = getFTensor0IntegrationWeight();
      auto t_val = getFTensor0FromVec(*(commonDataPtr->approxVals));
      auto t_coords = getFTensor1CoordsAtGaussPts();
 
      VectorDouble nf(nb_dofs, false);
      nf.clear();
 
      FTensor::Index<'i', 3> i;
      const double volume = getMeasure();
 
      auto t_row_base = data.getFTensor0N();
      double error = 0;
      for (int gg = 0; gg != nb_integration_pts; ++gg) {
 
        const double alpha = t_w * volume;
        double diff = t_val - AtomTest::approxFunction(t_coords(0), t_coords(1),
                                                      t_coords(2));
        error += alpha * pow(diff, 2);
 
        for (size_t r = 0; r != nb_dofs; ++r) {
          nf[r] += alpha * t_row_base * diff;
          ++t_row_base;
        }
 
        ++t_w;
        ++t_val;
        ++t_coords;
      }
 
      const int index = 0;
      CHKERR VecSetValue(commonDataPtr->L2Vec, index, error, ADD_VALUES);
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
  CHKERR setIntegrationRules();
  CHKERR createCommonData();
  CHKERR boundaryCondition();
  CHKERR assembleSystem();
  CHKERR solveSystem();
  CHKERR outputResults();
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

  MOFEM_LOG("WORLD", Sev::verbose) << "Dim " << simpleInterface->getDim();

  auto bit_level0 = simpleInterface->getBitRefLevel();

  auto &moab = mField.get_moab();

  auto refine_mesh = [&](auto bit_level1) {
    MoFEMFunctionBegin;

    auto refine = mField.getInterface<MeshRefinement>();

    auto meshset_level0_ptr = get_temp_meshset_ptr(moab);
    CHKERR mField.getInterface<BitRefManager>()->getEntitiesByRefLevel(
        bit_level0, BitRefLevel().set(), *meshset_level0_ptr);

    // random mesh refinement
    auto meshset_ref_edges_ptr = get_temp_meshset_ptr(moab);
    Range edges_to_refine;
    CHKERR moab.get_entities_by_type(*meshset_level0_ptr, MBEDGE,
                                     edges_to_refine);
    int ii = 0;
    for (Range::iterator eit = edges_to_refine.begin();
         eit != edges_to_refine.end(); eit++, ii++) {
      int numb = ii % 2;
      if (numb == 0) {
        CHKERR moab.add_entities(*meshset_ref_edges_ptr, &*eit, 1);
      }
    }
    CHKERR refine->addVerticesInTheMiddleOfEdges(
        *meshset_ref_edges_ptr, bit_level1, false, QUIET, 10000);
    if (simpleInterface->getDim() == 3) {
      CHKERR refine->refineTets(*meshset_level0_ptr, bit_level1, QUIET);
    } else if (simpleInterface->getDim() == 2) {
      CHKERR refine->refineTris(*meshset_level0_ptr, bit_level1, QUIET);
    } else {
      SETERRQ(PETSC_COMM_WORLD, MOFEM_ATOM_TEST_INVALID,
              "Dimension not handled by test");
    }

    MoFEMFunctionReturn(0);
  };

  BitRefLevel bit_level1;
  bit_level1.set(1);
  CHKERR refine_mesh(bit_level1);
  simpleInterface->getBitRefLevel() = bit_level0 | bit_level1;

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
 
//! [Set integration rule]
MoFEMErrorCode AtomTest::setIntegrationRules() {
  MoFEMFunctionBegin;
 
  auto rule = [](int, int, int p) -> int { return 2 * p; };
 
  PipelineManager *pipeline_mng = mField.getInterface<PipelineManager>();
  CHKERR pipeline_mng->setDomainLhsIntegrationRule(rule);
  CHKERR pipeline_mng->setDomainRhsIntegrationRule(rule);
 
  MoFEMFunctionReturn(0);
}
//! [Set integration rule]
 
//! [Create common data]
MoFEMErrorCode AtomTest::createCommonData() {
  MoFEMFunctionBegin;
  commonDataPtr = boost::make_shared<CommonData>();
  commonDataPtr->resVec = smartCreateDMVector(simpleInterface->getDM());
  commonDataPtr->L2Vec = createSmartVectorMPI(
      mField.get_comm(), (!mField.get_comm_rank()) ? 1 : 0, 1);
  commonDataPtr->approxVals = boost::make_shared<VectorDouble>();
  MoFEMFunctionReturn(0);
}
//! [Create common data]
 
//! [Boundary condition]
MoFEMErrorCode AtomTest::boundaryCondition() { return 0; }
//! [Boundary condition]
 
//! [Push operators to pipeline]
MoFEMErrorCode AtomTest::assembleSystem() {
  MoFEMFunctionBegin;
  PipelineManager *pipeline_mng = mField.getInterface<PipelineManager>();
  auto beta = [](const double, const double, const double) { return 1; };
  pipeline_mng->getOpDomainLhsPipeline().push_back(
      new OpDomainMass(FIELD_NAME, FIELD_NAME, beta));
  pipeline_mng->getOpDomainRhsPipeline().push_back(
      new OpDomainSource(FIELD_NAME, approxFunction));
  MoFEMFunctionReturn(0);
}
//! [Push operators to pipeline]
 
//! [Solve]
MoFEMErrorCode AtomTest::solveSystem() {
  MoFEMFunctionBegin;
  PipelineManager *pipeline_mng = mField.getInterface<PipelineManager>();
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
 
//! [Solve]
MoFEMErrorCode AtomTest::outputResults() {
  MoFEMFunctionBegin;
  MoFEMFunctionReturn(0);
}
//! [Postprocess results]
 
//! [Check results]
MoFEMErrorCode AtomTest::checkResults() {
  MoFEMFunctionBegin;
  PipelineManager *pipeline_mng = mField.getInterface<PipelineManager>();
  pipeline_mng->getDomainLhsFE().reset();
  pipeline_mng->getDomainRhsFE().reset();
  pipeline_mng->getOpDomainRhsPipeline().clear();
  pipeline_mng->getOpDomainRhsPipeline().push_back(
      new OpCalculateScalarFieldValues(FIELD_NAME, commonDataPtr->approxVals));
  pipeline_mng->getOpDomainRhsPipeline().push_back(
      new OpError<FIELD_DIM>(commonDataPtr));
  CHKERR pipeline_mng->loopFiniteElements();
  CHKERR VecAssemblyBegin(commonDataPtr->L2Vec);
  CHKERR VecAssemblyEnd(commonDataPtr->L2Vec);
  CHKERR VecAssemblyBegin(commonDataPtr->resVec);
  CHKERR VecAssemblyEnd(commonDataPtr->resVec);
  double nrm2;
  CHKERR VecNorm(commonDataPtr->resVec, NORM_2, &nrm2);
  const double *array;
  CHKERR VecGetArrayRead(commonDataPtr->L2Vec, &array);
  if (mField.get_comm_rank() == 0)
    PetscPrintf(PETSC_COMM_SELF, "Error %6.4e Vec norm %6.4e\n", std::sqrt(array[0]),
                nrm2);
  CHKERR VecRestoreArrayRead(commonDataPtr->L2Vec, &array);
  constexpr double eps = 1e-8;
  if (nrm2 > eps)
    SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
            "Not converged solution");
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
