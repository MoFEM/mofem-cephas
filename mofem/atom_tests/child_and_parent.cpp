/**
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
  using DomianParentEle = FaceElementForcesAndSourcesCoreOnChildParent;
};

template <> struct ElementsAndOps<3> {
  using DomainEle = VolumeElementForcesAndSourcesCore;
  using DomainEleOp = DomainEle::UserDataOperator;
};

using DomainEle = ElementsAndOps<SPACE_DIM>::DomainEle;
using DomainParentEle = ElementsAndOps<SPACE_DIM>::DomianParentEle;
using DomainEleOp = DomainEle::UserDataOperator;
using EntData = EntitiesFieldData::EntData;

template <int FIELD_DIM> struct ApproxFieldFunction;

template <> struct ApproxFieldFunction<1> {
  double operator()(const double x, const double y, const double z) {
    return x * x + y * y + x * y + pow(x, 3) + pow(y, 3) + pow(x, 4) +
           pow(y, 4);
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
  MoFEMErrorCode assembleSystem();
  MoFEMErrorCode solveSystem();
  MoFEMErrorCode
  checkResults(boost::function<bool(FEMethod *fe_method_ptr)> test_bit);
  MoFEMErrorCode refineResults();
  struct CommonData {
    boost::shared_ptr<VectorDouble> approxVals;
    SmartPetscObj<Vec> L2Vec;
    SmartPetscObj<Vec> resVec;
  };

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
  CHKERR assembleSystem();
  CHKERR solveSystem();

  auto test_bit_child = [](FEMethod *fe_ptr) {
    const auto &bit = fe_ptr->numeredEntFiniteElementPtr->getBitRefLevel();
    MOFEM_LOG("SELF", Sev::noisy) << bit << " " << bit.test(0);
    return bit.test(1);
  };

  CHKERR checkResults(test_bit_child);
  CHKERR refineResults();
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
    CHKERR refine->addVerticesInTheMiddleOfEdges(*meshset_ref_edges_ptr,
                                                 bit_level1, false, VERBOSE);
    if (simpleInterface->getDim() == 3) {
      CHKERR refine->refineTets(*meshset_level0_ptr, bit_level1, VERBOSE);
    } else if (simpleInterface->getDim() == 2) {
      CHKERR refine->refineTris(*meshset_level0_ptr, bit_level1, VERBOSE);
    } else {
      SETERRQ(PETSC_COMM_WORLD, MOFEM_ATOM_TEST_INVALID,
              "Dimension not handled by test");
    }

    MoFEMFunctionReturn(0);
  };

  BitRefLevel bit_level1;
  bit_level1.set(1);
  CHKERR refine_mesh(bit_level1);
  simpleInterface->getBitRefLevel() = BitRefLevel().set();
  simpleInterface->getBitRefLevelMask() = BitRefLevel().set();

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

  CHKERR mField.getInterface<ProblemsManager>()->removeDofsOnEntities(
      simpleInterface->getProblemName(), FIELD_NAME, BitRefLevel().set(0),
      BitRefLevel().set(0));

  MoFEMFunctionReturn(0);
}
//! [Set up problem]

boost::shared_ptr<DomainEle> domainChildLhs, domainChildRhs;

//! [Push operators to pipeline]
MoFEMErrorCode AtomTest::assembleSystem() {
  MoFEMFunctionBegin;
  PipelineManager *pipeline_mng = mField.getInterface<PipelineManager>();

  auto rule = [](int, int, int p) -> int { return 2 * p; };

  CHKERR pipeline_mng->setDomainLhsIntegrationRule(rule);
  CHKERR pipeline_mng->setDomainRhsIntegrationRule(rule);

  auto test_bit_parent = [](FEMethod *fe_ptr) {
    const auto &bit = fe_ptr->numeredEntFiniteElementPtr->getBitRefLevel();
    MOFEM_LOG("SELF", Sev::noisy) << bit << " " << bit.test(0);
    return bit.test(0);
  };

  pipeline_mng->getDomainLhsFE()->exeTestHook = test_bit_parent;
  pipeline_mng->getDomainRhsFE()->exeTestHook = test_bit_parent;

  auto beta = [](const double, const double, const double) { return 1; };

  // Make aliased shared pointer, and create child element
  domainChildLhs = boost::make_shared<DomainEle>(mField);
  domainChildLhs->getRuleHook = rule;
  domainChildLhs->getOpPtrVector().push_back(
      new OpDomainMass(FIELD_NAME, FIELD_NAME, beta));
      
  domainChildRhs = boost::make_shared<DomainEle>(mField);
  domainChildLhs->getRuleHook = rule;
  domainChildRhs->getOpPtrVector().push_back(
      new OpDomainSource(FIELD_NAME, approxFunction));

  auto parent_op_lhs = new DomainEleOp(NOSPACE, DomainEleOp::OPSPACE);
  parent_op_lhs->doWorkRhsHook = [&](DataOperator *op_ptr, int side,
                                     EntityType type,
                                     EntitiesFieldData::EntData &data) {
    auto domain_op = static_cast<DomainEleOp *>(op_ptr);
    MoFEMFunctionBegin;

    MOFEM_LOG("SELF", Sev::noisy) << "LHS Pipeline FE";

    if (!domainChildLhs)
      SETERRQ(PETSC_COMM_SELF, MOFEM_ATOM_TEST_INVALID, "FE not allocated");

    auto &bit =
        domain_op->getFEMethod()->numeredEntFiniteElementPtr->getBitRefLevel();
    if (bit == BitRefLevel().set(0)) {
      CHKERR domain_op->loopChildren(domain_op->getFEName(),
                                     domainChildLhs.get(), VERBOSE, Sev::noisy);
    } else {
      CHKERR domain_op->loopThis(domain_op->getFEName(), domainChildLhs.get(),
                                 VERBOSE, Sev::noisy);
    }
    MoFEMFunctionReturn(0);
  };

  auto parent_op_rhs = new DomainEleOp(NOSPACE, DomainEleOp::OPSPACE);
  parent_op_rhs->doWorkRhsHook = [&](DataOperator *op_ptr, int side,
                                     EntityType type,
                                     EntitiesFieldData::EntData &data) {
    auto domain_op = static_cast<DomainEleOp *>(op_ptr);
    MoFEMFunctionBegin;

    MOFEM_LOG("SELF", Sev::noisy) << "RHS Pipeline FE";

    if (!domainChildRhs)
      SETERRQ(PETSC_COMM_SELF, MOFEM_ATOM_TEST_INVALID, "FE not allocated");

    auto &bit =
        domain_op->getFEMethod()->numeredEntFiniteElementPtr->getBitRefLevel();
    if (bit == BitRefLevel().set(0)) {
      CHKERR domain_op->loopChildren(domain_op->getFEName(),
                                     domainChildRhs.get(), VERBOSE, Sev::noisy);
    } else if ((bit & BitRefLevel().set(0)).any()) {
      CHKERR domain_op->loopThis(domain_op->getFEName(), domainChildRhs.get(),
                                 VERBOSE, Sev::noisy);
    }
    MoFEMFunctionReturn(0);
  };

  pipeline_mng->getOpDomainLhsPipeline().push_back(parent_op_lhs);
  pipeline_mng->getOpDomainRhsPipeline().push_back(parent_op_rhs);

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

//! [Solve]
MoFEMErrorCode AtomTest::refineResults() {
  MoFEMFunctionBegin;

  auto &moab = mField.get_moab();

  auto bit_level0 = BitRefLevel().set(0);
  auto bit_level1 = BitRefLevel().set(1);
  auto bit_level2 = BitRefLevel().set(2);

  auto refine_mesh = [&]() {
    MoFEMFunctionBegin;

    auto refine = mField.getInterface<MeshRefinement>();

    auto meshset_level1_ptr = get_temp_meshset_ptr(moab);
    CHKERR mField.getInterface<BitRefManager>()->getEntitiesByDimAndRefLevel(
        bit_level1, BitRefLevel().set(), simpleInterface->getDim(),
        *meshset_level1_ptr);

    // random mesh refinement
    auto meshset_ref_edges_ptr = get_temp_meshset_ptr(moab);
    Range edges_to_refine;
    CHKERR mField.getInterface<BitRefManager>()->getEntitiesByTypeAndRefLevel(
        bit_level1, BitRefLevel().set(), MBEDGE, edges_to_refine);

    CHKERR refine->addVerticesInTheMiddleOfEdges(edges_to_refine, bit_level2,
                                                 VERBOSE);
    if (simpleInterface->getDim() == 3) {
      CHKERR refine->refineTets(*meshset_level1_ptr, bit_level2, VERBOSE);
    } else if (simpleInterface->getDim() == 2) {
      CHKERR refine->refineTris(*meshset_level1_ptr, bit_level2, VERBOSE);
    } else {
      SETERRQ(PETSC_COMM_WORLD, MOFEM_ATOM_TEST_INVALID,
              "Dimension not handled by test");
    }

    Range meshsets;
    CHKERR moab.get_entities_by_type(0, MBENTITYSET, meshsets, true);
    for (auto m : meshsets) {
      CHKERR mField.getInterface<BitRefManager>()
          ->updateMeshsetByEntitiesChildren(m, bit_level2, m, MBMAXTYPE, false);
    }

    MoFEMFunctionReturn(0);
  };

  CHKERR refine_mesh();

  simpleInterface->getBitRefLevel() = bit_level1 | bit_level2;
  simpleInterface->getBitRefLevelMask() = BitRefLevel().set();

  CHKERR simpleInterface->reSetUp();

  CHKERR mField.getInterface<ProblemsManager>()->removeDofsOnEntities(
      simpleInterface->getProblemName(), FIELD_NAME, BitRefLevel().set(),
      bit_level0 | bit_level1);

  auto project_data = [&]() {
    MoFEMFunctionBegin;

    PipelineManager *pipeline_mng = mField.getInterface<PipelineManager>();

    pipeline_mng->getDomainLhsFE().reset();
    pipeline_mng->getDomainRhsFE().reset();

    auto rule = [](int, int, int p) -> int { return 2 * p; };

    CHKERR pipeline_mng->setDomainLhsIntegrationRule(rule);
    CHKERR pipeline_mng->setDomainRhsIntegrationRule(rule);

    auto test_bit_ref = [](FEMethod *fe_ptr) {
      const auto &bit = fe_ptr->numeredEntFiniteElementPtr->getBitRefLevel();
      MOFEM_LOG("SELF", Sev::noisy) << "ref : " << bit << " " << bit.test(2);
      return bit.test(2);
    };

    pipeline_mng->getDomainLhsFE()->exeTestHook = test_bit_ref;
    pipeline_mng->getDomainRhsFE()->exeTestHook = test_bit_ref;

    auto beta = [](const double, const double, const double) { return 1; };
    auto field_vals_ptr = boost::make_shared<VectorDouble>();

    auto domainParentRhs = boost::make_shared<DomainParentEle>(mField);
    domainParentRhs->getOpPtrVector().push_back(
        new OpCalculateScalarFieldValues(FIELD_NAME, field_vals_ptr));

    pipeline_mng->getOpDomainLhsPipeline().push_back(
        new OpDomainMass(FIELD_NAME, FIELD_NAME, beta));
    pipeline_mng->getOpDomainRhsPipeline().push_back(
        new OpRunParent(domainParentRhs, bit_level2, bit_level2,
                        domainParentRhs, bit_level2, BitRefLevel().set()));

    using OpDomainTimesScalarField = FormsIntegrators<DomainEleOp>::Assembly<
        PETSC>::LinearForm<GAUSS>::OpBaseTimesScalarField<1>;
    pipeline_mng->getOpDomainRhsPipeline().push_back(
        new OpDomainTimesScalarField(FIELD_NAME, field_vals_ptr, beta));

    struct OpCheckGaussCoords : public DomainEleOp {
      OpCheckGaussCoords() : DomainEleOp(NOSPACE, DomainEleOp::OPSPACE) {}

      MoFEMErrorCode doWork(int side, EntityType type,
                            EntitiesFieldData::EntData &data) {
        MoFEMFunctionBegin;

        MatrixDouble parent_coords;

        DomainParentEle parent_fe(getPtrFE()->mField);
        auto domain_op = new DomainEleOp(NOSPACE, DomainEleOp::OPSPACE);
        domain_op->doWorkRhsHook =
            [&](DataOperator *op_ptr, int side, EntityType type,
                EntitiesFieldData::EntData &data) {
              MoFEMFunctionBegin;
              parent_coords =
                  static_cast<DomainEleOp *>(op_ptr)->getCoordsAtGaussPts();
              MoFEMFunctionReturn(0);
            };
        parent_fe.getOpPtrVector().push_back(domain_op);

        CHKERR loopParent(getFEName(), &parent_fe);

        MatrixDouble child_coords = getCoordsAtGaussPts();
        child_coords -= parent_coords;

        MOFEM_LOG("SELF", Sev::noisy) << "Corrds diffs" << child_coords;

        double n = 0;
        for (auto d : child_coords.data())
          n += d * d;

        if (sqrt(n) > 1e-12)
          SETERRQ1(PETSC_COMM_SELF, MOFEM_ATOM_TEST_INVALID,
                   "Parent and child global coords at integration points are "
                   "diffrent norm = %3.2e",
                   sqrt(n));

        MoFEMFunctionReturn(0);
      }
    };

    pipeline_mng->getOpDomainRhsPipeline().push_back(new OpCheckGaussCoords());

    CHKERR solveSystem();

    simpleInterface->getBitRefLevel() = bit_level2;
    simpleInterface->getBitRefLevelMask() = BitRefLevel().set();
    CHKERR simpleInterface->reSetUp();

    CHKERR checkResults([](FEMethod *fe_ptr) { return true; });

    MoFEMFunctionReturn(0);
  };

  CHKERR project_data();

  MoFEMFunctionReturn(0);
}
//! [Postprocess results]

//! [Check results]
MoFEMErrorCode AtomTest::checkResults(
    boost::function<bool(FEMethod *fe_method_ptr)> test_bit) {
  MoFEMFunctionBegin;
  PipelineManager *pipeline_mng = mField.getInterface<PipelineManager>();
  pipeline_mng->getDomainLhsFE().reset();
  pipeline_mng->getDomainRhsFE().reset();

  auto rule = [](int, int, int p) -> int { return 2 * p + 1; };
  CHKERR pipeline_mng->setDomainRhsIntegrationRule(rule);
  pipeline_mng->getDomainRhsFE()->exeTestHook = test_bit;

  auto common_data_ptr = boost::make_shared<CommonData>();
  common_data_ptr->resVec = smartCreateDMVector(simpleInterface->getDM());
  common_data_ptr->L2Vec = createSmartVectorMPI(
      mField.get_comm(), (!mField.get_comm_rank()) ? 1 : 0, 1);
  common_data_ptr->approxVals = boost::make_shared<VectorDouble>();

  pipeline_mng->getOpDomainRhsPipeline().push_back(
      new OpCalculateScalarFieldValues(FIELD_NAME,
      common_data_ptr->approxVals));
  pipeline_mng->getOpDomainRhsPipeline().push_back(
      new OpError<FIELD_DIM>(common_data_ptr));

  CHKERR VecZeroEntries(common_data_ptr->L2Vec);
  CHKERR VecZeroEntries(common_data_ptr->resVec);

  CHKERR pipeline_mng->loopFiniteElements();

  CHKERR VecAssemblyBegin(common_data_ptr->L2Vec);
  CHKERR VecAssemblyEnd(common_data_ptr->L2Vec);
  CHKERR VecAssemblyBegin(common_data_ptr->resVec);
  CHKERR VecAssemblyEnd(common_data_ptr->resVec);
  double nrm2;
  CHKERR VecNorm(common_data_ptr->resVec, NORM_2, &nrm2);
  const double *array;
  CHKERR VecGetArrayRead(common_data_ptr->L2Vec, &array);
  MOFEM_LOG_C("WORLD", Sev::inform, "Error %6.4e Vec norm %6.4e\n",
              std::sqrt(array[0]), nrm2);
  CHKERR VecRestoreArrayRead(common_data_ptr->L2Vec, &array);
  
  constexpr double eps = 1e-8;
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
