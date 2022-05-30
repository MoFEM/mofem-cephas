/**
 * \example hanging_node_approx.cpp
 *
 * Tetsing approximation with hanging nodes.
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
constexpr int nb_ref_levels = 3; ///< Three levels of refinement

template <int DIM> struct ElementsAndOps {};

template <> struct ElementsAndOps<2> {
  using DomainEle = PipelineManager::FaceEle;
  using DomainEleOp = DomainEle::UserDataOperator;
  using DomianParentEle = FaceElementForcesAndSourcesCoreOnChildParent;
  using SkeletonEle = PipelineManager::EdgeEle;
  using SkeletonEleOp = SkeletonEle::UserDataOperator;
  using SkeletonParentEle = EdgeElementForcesAndSourcesCoreOnChildParent;
};

template <> struct ElementsAndOps<3> {
  using DomainEle = VolumeElementForcesAndSourcesCore;
  using DomainEleOp = DomainEle::UserDataOperator;
};

using DomainEle = ElementsAndOps<SPACE_DIM>::DomainEle;
using DomainParentEle = ElementsAndOps<SPACE_DIM>::DomianParentEle;
using DomainEleOp = DomainEle::UserDataOperator;
using SkeletonEle = ElementsAndOps<SPACE_DIM>::SkeletonEle;
using SkeletonEleOp = SkeletonEle::UserDataOperator;
using SkeletonParentEle = ElementsAndOps<SPACE_DIM>::SkeletonParentEle;

using EntData = EntitiesFieldData::EntData;

template <int FIELD_DIM> struct ApproxFieldFunction;
template <int FIELD_DIM> struct ApproxFieldFunctionDerivative;

/**
 * @brief third order polynomial used for testing
 *
 */
template <> struct ApproxFieldFunction<1> {
  auto operator()(const double x, const double y, const double z) {
    return x * x + y * y + x * y * y + x * x * y;
  }
};

/**
 * @brief third order polynomial used for testing
 *
 */
template <> struct ApproxFieldFunctionDerivative<1> {
  auto operator()(const double x, const double y, const double z) {
    // x * x + y * y + x * y * y + x * x * y

    return FTensor::Tensor1<double, SPACE_DIM>{2 * x + y * y + 2 * x * y,
                                               2 * y + 2 * x * y + x * x};
  }
};

/**
 * @brief evaluate mass matrix
 *
 */
using OpDomainMass = FormsIntegrators<DomainEleOp>::Assembly<
    PETSC>::BiLinearForm<GAUSS>::OpMass<1, FIELD_DIM>;

/**
 * @brief evaluate source, i.e. rhs vector
 *
 */
using OpDomainSource = FormsIntegrators<DomainEleOp>::Assembly<
    PETSC>::LinearForm<GAUSS>::OpSource<1, FIELD_DIM>;

/**
 * @brief set bit
 *
 */
auto bit = [](auto l) { return BitRefLevel().set(l); };

/**
 * @brief set bit to marker
 *
 * Marker is used to mark field entities on skin on which we have hanging nodes
 */
auto marker = [](auto l) {
  return BitRefLevel().set(BITREFLEVEL_SIZE - 1 - l);
};

/**
 * @brief set levels of projection operators, which project field data from
 * parent entities, to child, up to to level, i.e. last mesh refinement.
 *
 */
template <typename PARENT_FE>
auto set_parent_dofs(MoFEM::Interface &m_field,
                     boost::shared_ptr<FEMethod> &fe_top,
                     ForcesAndSourcesCore::UserDataOperator::OpType op,
                     int verbosity, LogManager::SeverityLevel sev) {

  auto jac_ptr = boost::make_shared<MatrixDouble>();
  auto inv_jac_ptr = boost::make_shared<MatrixDouble>();
  auto det_ptr = boost::make_shared<VectorDouble>();

  BitRefLevel bit_marker;
  for (auto l = 1; l <= nb_ref_levels; ++l)
    bit_marker |= marker(l);

  boost::function<void(boost::shared_ptr<ForcesAndSourcesCore>, int)>
      add_parent_level =
          [&](boost::shared_ptr<ForcesAndSourcesCore> parent_fe_pt, int level) {
            if (level > 0) {

              auto fe_ptr_current = boost::shared_ptr<ForcesAndSourcesCore>(
                  new PARENT_FE(m_field));
              if (op == DomainEleOp::OPSPACE) {
                fe_ptr_current->getOpPtrVector().push_back(
                    new OpCalculateHOJacForFace(jac_ptr));
                fe_ptr_current->getOpPtrVector().push_back(
                    new OpInvertMatrix<2>(jac_ptr, det_ptr, inv_jac_ptr));
                fe_ptr_current->getOpPtrVector().push_back(
                    new OpSetInvJacH1ForFace(inv_jac_ptr));
              }

              add_parent_level(
                  boost::dynamic_pointer_cast<ForcesAndSourcesCore>(
                      fe_ptr_current),
                  level - 1);

              if (op == DomainEleOp::OPSPACE) {

                parent_fe_pt->getOpPtrVector().push_back(

                    new OpAddParentEntData(

                        H1, op, fe_ptr_current,

                        BitRefLevel().set(), bit(0).flip(),

                        bit_marker, BitRefLevel().set(),

                        verbosity, sev));

              } else {

                parent_fe_pt->getOpPtrVector().push_back(

                    new OpAddParentEntData(

                        FIELD_NAME, op, fe_ptr_current,

                        BitRefLevel().set(), bit(0).flip(),

                        bit_marker, BitRefLevel().set(),

                        verbosity, sev));
              }
            }
          };

  add_parent_level(boost::dynamic_pointer_cast<ForcesAndSourcesCore>(fe_top),
                   nb_ref_levels);
};

/**
 * @brief lambda function used to select elements on which finite element
 * pipelines are executed.
 *
 * @note childs elements on pipeline, retrive data from parents using operators
 * pushed by \ref set_parent_dofs
 *
 */
auto test_bit_child = [](FEMethod *fe_ptr) {
  return fe_ptr->numeredEntFiniteElementPtr->getBitRefLevel().test(
      nb_ref_levels);
};

struct AtomTest {

  AtomTest(MoFEM::Interface &m_field) : mField(m_field) {}

  MoFEMErrorCode runProblem();

private:
  MoFEM::Interface &mField;
  Simple *simpleInterface;

  static ApproxFieldFunction<FIELD_DIM> approxFunction;
  static ApproxFieldFunctionDerivative<FIELD_DIM> divApproxFunction;

  /**
   * @brief red mesh and randomly refine three times
   *
   * @return MoFEMErrorCode
   */
  MoFEMErrorCode readMesh();

  /**
   * @brief add field, and set up problem
   *
   * @return MoFEMErrorCode
   */
  MoFEMErrorCode setupProblem();
  MoFEMErrorCode assembleSystem();
  MoFEMErrorCode solveSystem();
  MoFEMErrorCode checkResults();
  MoFEMErrorCode printResults();
  struct CommonData {
    boost::shared_ptr<VectorDouble> approxVals;
    boost::shared_ptr<MatrixDouble> divApproxVals;
    SmartPetscObj<Vec> L2Vec;
    SmartPetscObj<Vec> resVec;
  };

  template <int FIELD_DIM> struct OpError;

  template <int FIELD_DIM> struct OpErrorSkel;
};

ApproxFieldFunction<FIELD_DIM> AtomTest::approxFunction =
    ApproxFieldFunction<FIELD_DIM>();
ApproxFieldFunctionDerivative<FIELD_DIM> AtomTest::divApproxFunction =
    ApproxFieldFunctionDerivative<FIELD_DIM>();
template <> struct AtomTest::OpError<1> : public DomainEleOp {
  boost::shared_ptr<CommonData> commonDataPtr;
  OpError(boost::shared_ptr<CommonData> &common_data_ptr)
      : DomainEleOp(FIELD_NAME, OPROW), commonDataPtr(common_data_ptr) {}
  MoFEMErrorCode doWork(int side, EntityType type, EntData &data) {
    MoFEMFunctionBegin;

    if (const size_t nb_dofs = data.getIndices().size()) {

      FTensor::Index<'i', SPACE_DIM> i;

      const int nb_integration_pts = getGaussPts().size2();
      auto t_w = getFTensor0IntegrationWeight();
      auto t_val = getFTensor0FromVec(*(commonDataPtr->approxVals));
      auto t_grad_val =
          getFTensor1FromMat<SPACE_DIM>(*(commonDataPtr->divApproxVals));
      auto t_coords = getFTensor1CoordsAtGaussPts();

      VectorDouble nf(nb_dofs, false);
      nf.clear();

      const double volume = getMeasure();

      auto t_row_base = data.getFTensor0N();
      double error = 0;
      for (int gg = 0; gg != nb_integration_pts; ++gg) {

        const double alpha = t_w * volume;
        double diff = t_val - AtomTest::approxFunction(t_coords(0), t_coords(1),
                                                       t_coords(2));

        auto t_grad_diff =
            AtomTest::divApproxFunction(t_coords(0), t_coords(1), t_coords(2));
        t_grad_diff(i) -= t_grad_val(i);

        error += alpha * (pow(diff, 2) + t_grad_diff(i) * t_grad_diff(i));

        for (size_t r = 0; r != nb_dofs; ++r) {
          nf[r] += alpha * t_row_base * diff;
          ++t_row_base;
        }

        ++t_w;
        ++t_val;
        ++t_grad_val;
        ++t_coords;
      }

      const int index = 0;
      CHKERR VecSetValue(commonDataPtr->L2Vec, index, error, ADD_VALUES);
      CHKERR VecSetValues(commonDataPtr->resVec, data, &nf[0], ADD_VALUES);
    }

    MoFEMFunctionReturn(0);
  }
};

template <> struct AtomTest::OpErrorSkel<1> : public SkeletonEleOp {
  boost::shared_ptr<CommonData> commonDataPtr;
  OpErrorSkel(boost::shared_ptr<CommonData> &common_data_ptr)
      : SkeletonEleOp(H1, OPSPACE), commonDataPtr(common_data_ptr) {}
  MoFEMErrorCode doWork(int side, EntityType type, EntData &data) {
    MoFEMFunctionBegin;

    FTensor::Index<'i', SPACE_DIM> i;

    const int nb_integration_pts = getGaussPts().size2();
    auto t_w = getFTensor0IntegrationWeight();
    auto t_val = getFTensor0FromVec(*(commonDataPtr->approxVals));
    auto t_coords = getFTensor1CoordsAtGaussPts();

    const double volume = getMeasure();

    double error2 = 0;
    for (int gg = 0; gg != nb_integration_pts; ++gg) {

      const double alpha = t_w * volume;
      double diff = t_val - AtomTest::approxFunction(t_coords(0), t_coords(1),
                                                     t_coords(2));
      error2 += alpha * (pow(diff, 2));

      ++t_w;
      ++t_val;
      ++t_coords;
    }

    constexpr double eps = 1e-8;
    if (sqrt(error2) > eps)
      SETERRQ1(PETSC_COMM_SELF, MOFEM_ATOM_TEST_INVALID,
               "Error on boundary = %6.4e", sqrt(error2));

    MOFEM_LOG("SELF", Sev::noisy) << "Skeleton error " << sqrt(error2);

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
  CHKERR printResults();
  MoFEMFunctionReturn(0);
}
//! [Run programme]

//! [Read mesh]
MoFEMErrorCode AtomTest::readMesh() {
  BitRefManager *bit_mng = mField.getInterface<BitRefManager>();
  ParallelComm *pcomm =
      ParallelComm::get_pcomm(&mField.get_moab(), MYPCOMM_INDEX);
  Skinner skin(&mField.get_moab());
  MoFEMFunctionBegin;

  CHKERR mField.getInterface(simpleInterface);
  CHKERR simpleInterface->getOptions();
  CHKERR simpleInterface->loadFile();

  MOFEM_LOG("WORLD", Sev::verbose) << "Dim " << simpleInterface->getDim();

  auto &moab = mField.get_moab();

  Range level0_ents;
  CHKERR mField.getInterface<BitRefManager>()->getEntitiesByDimAndRefLevel(
      bit(0), BitRefLevel().set(), SPACE_DIM, level0_ents);
  Range level0_skin;
  CHKERR skin.find_skin(0, level0_ents, false, level0_skin);
  CHKERR pcomm->filter_pstatus(level0_skin,
                               PSTATUS_SHARED | PSTATUS_MULTISHARED,
                               PSTATUS_NOT, -1, nullptr);

  auto refine_mesh = [&](auto l) {
    MoFEMFunctionBegin;

    auto refine = mField.getInterface<MeshRefinement>();

    auto meshset_level0_ptr = get_temp_meshset_ptr(moab);
    CHKERR bit_mng->getEntitiesByDimAndRefLevel(bit(l - 1), BitRefLevel().set(),
                                                SPACE_DIM, *meshset_level0_ptr);

    // random mesh refinement
    auto meshset_ref_edges_ptr = get_temp_meshset_ptr(moab);

    Range els;
    CHKERR moab.get_entities_by_dimension(*meshset_level0_ptr, SPACE_DIM, els);
    CHKERR bit_mng->filterEntitiesByRefLevel(bit(l - 1), bit(l - 1), els);

    Range ele_to_refine;

    if (l == 1) {
      int ii = 0;
      for (auto t : els) {
        if ((ii % 2)) {
          ele_to_refine.insert(t);
          std::vector<EntityHandle> adj_edges;
          CHKERR mField.get_moab().get_adjacencies(&t, 1, SPACE_DIM - 1, false,
                                                   adj_edges);
          CHKERR moab.add_entities(*meshset_ref_edges_ptr, &*adj_edges.begin(),
                                   adj_edges.size());
        }
        ++ii;
      }
    } else {
      Range level_skin;
      CHKERR skin.find_skin(0, els, false, level_skin);
      CHKERR pcomm->filter_pstatus(level_skin,
                                   PSTATUS_SHARED | PSTATUS_MULTISHARED,
                                   PSTATUS_NOT, -1, nullptr);
      level_skin = subtract(level_skin, level0_skin);
      Range adj;
      CHKERR mField.get_moab().get_adjacencies(level_skin, SPACE_DIM, false,
                                               adj, moab::Interface::UNION);
      els = subtract(els, adj);
      ele_to_refine.merge(els);
      Range adj_edges;
      CHKERR mField.get_moab().get_adjacencies(
          els, SPACE_DIM - 1, false, adj_edges, moab::Interface::UNION);
      CHKERR moab.add_entities(*meshset_ref_edges_ptr, adj_edges);
    }

    CHKERR refine->addVerticesInTheMiddleOfEdges(*meshset_ref_edges_ptr, bit(l),
                                                 false, VERBOSE);
    CHKERR refine->refineTrisHangingNodes(*meshset_level0_ptr, bit(l), VERBOSE);
    CHKERR bit_mng->updateRangeByChildren(level0_skin, level0_skin);

    CHKERR bit_mng->writeBitLevelByDim(
        bit(l), BitRefLevel().set(), SPACE_DIM,
        (boost::lexical_cast<std::string>(l) + "_ref_mesh.vtk").c_str(), "VTK",
        "");
    CHKERR bit_mng->writeBitLevelByDim(
        bit(l), bit(l), MBTRI,
        (boost::lexical_cast<std::string>(l) + "_only_ref_mesh.vtk").c_str(),
        "VTK", "");

    MoFEMFunctionReturn(0);
  };

  auto mark_skins = [&](auto l, auto m) {
    MoFEMFunctionBegin;
    Range ents;
    CHKERR bit_mng->getEntitiesByDimAndRefLevel(bit(l), bit(l), SPACE_DIM,
                                                ents);
    Range level_skin;
    CHKERR skin.find_skin(0, ents, false, level_skin);
    CHKERR pcomm->filter_pstatus(level_skin,
                                 PSTATUS_SHARED | PSTATUS_MULTISHARED,
                                 PSTATUS_NOT, -1, nullptr);
    level_skin = subtract(level_skin, level0_skin);
    CHKERR mField.get_moab().get_adjacencies(level_skin, 0, false, level_skin,
                                             moab::Interface::UNION);
    CHKERR bit_mng->addBitRefLevel(level_skin, marker(m));
    MoFEMFunctionReturn(0);
  };

  BitRefLevel bit_sum;
  for (auto l = 0; l != nb_ref_levels; ++l) {
    CHKERR refine_mesh(l + 1);
    CHKERR mark_skins(l, l + 1);
    CHKERR mark_skins(l + 1, l + 1);
    bit_sum |= bit(l);
  }
  bit_sum |= bit(nb_ref_levels);

  simpleInterface->getBitRefLevel() = bit_sum;
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
  CHKERR simpleInterface->addSkeletonField(FIELD_NAME, H1,
                                           AINSWORTH_LEGENDRE_BASE, FIELD_DIM);

  constexpr int order = 3;
  CHKERR simpleInterface->setFieldOrder(FIELD_NAME, order);

  // Simple interface will resolve adjacency to DOFs of parent of the element.
  // Using that information MAtrixManager  allocate appropriately size of
  // matrix.
  simpleInterface->getParentAdjacencies() = true;
  BitRefLevel bit_marker;
  for (auto l = 1; l <= nb_ref_levels; ++l)
    bit_marker |= marker(l);
  simpleInterface->getBitAdjEnt() = bit_marker;

  CHKERR simpleInterface->setUp();

  BitRefManager *bit_mng = mField.getInterface<BitRefManager>();
  ProblemsManager *prb_mng = mField.getInterface<ProblemsManager>();

  // remove obsolete DOFs from problem

  for (int l = 0; l != nb_ref_levels; ++l) {
    CHKERR prb_mng->removeDofsOnEntities(simpleInterface->getProblemName(),
                                         FIELD_NAME, bit(l), bit(l));
    CHKERR prb_mng->removeDofsOnEntities(simpleInterface->getProblemName(),
                                         FIELD_NAME, marker(l + 1),
                                         bit(l).flip());
  }
  MoFEMFunctionReturn(0);
}
//! [Set up problem]

//! [Push operators to pipeline]
MoFEMErrorCode AtomTest::assembleSystem() {
  MoFEMFunctionBegin;
  PipelineManager *pipeline_mng = mField.getInterface<PipelineManager>();

  auto rule = [](int, int, int p) -> int { return 2 * p + 1; };

  CHKERR pipeline_mng->setDomainLhsIntegrationRule(rule);
  CHKERR pipeline_mng->setDomainRhsIntegrationRule(rule);

  pipeline_mng->getDomainLhsFE()->exeTestHook = test_bit_child;
  pipeline_mng->getDomainRhsFE()->exeTestHook = test_bit_child;

  auto beta = [](const double, const double, const double) { return 1; };
  set_parent_dofs<DomainParentEle>(mField, pipeline_mng->getDomainLhsFE(),
                                   DomainEleOp::OPSPACE, QUIET, Sev::noisy);
  set_parent_dofs<DomainParentEle>(mField, pipeline_mng->getDomainLhsFE(),
                                   DomainEleOp::OPROW, QUIET, Sev::noisy);
  set_parent_dofs<DomainParentEle>(mField, pipeline_mng->getDomainLhsFE(),
                                   DomainEleOp::OPCOL, QUIET, Sev::noisy);
  pipeline_mng->getOpDomainLhsPipeline().push_back(
      new OpDomainMass(FIELD_NAME, FIELD_NAME, beta));

  auto field_op_row = new ForcesAndSourcesCore::UserDataOperator(
      FIELD_NAME, DomainEleOp::OPROW);
  field_op_row->doWorkRhsHook = [](DataOperator *op_ptr, int side,
                                   EntityType type,
                                   EntitiesFieldData::EntData &data) {
    MoFEMFunctionBegin;
    if (type == MBENTITYSET) {

      MOFEM_LOG("SELF", Sev::verbose)
          << "ROW: side/type: " << side << "/" << CN::EntityTypeName(type)
          << " op space/base: " << FieldSpaceNames[data.getSpace()] << "/"
          << ApproximationBaseNames[data.getBase()] << " DOFs "
          << data.getIndices() << " nb base functions " << data.getN().size2()
          << " nb base functions integration points " << data.getN().size1();

      auto get_bool = [](auto fe, auto bit) {
        return (bit & fe->getBitRefLevel()).any();
      };

      for (auto &field_ent : data.getFieldEntities()) {
        MOFEM_LOG("SELF", Sev::verbose)
            << "\t" << CN::EntityTypeName(field_ent->getEntType());
      }
    }
    MoFEMFunctionReturn(0);
  };

  set_parent_dofs<DomainParentEle>(
      mField, pipeline_mng->getDomainRhsFE(), DomainEleOp::OPSPACE, VERBOSE,
      Sev::verbose);
  set_parent_dofs<DomainParentEle>(mField, pipeline_mng->getDomainRhsFE(),
                                   DomainEleOp::OPROW, VERBOSE, Sev::noisy);
  pipeline_mng->getOpDomainRhsPipeline().push_back(field_op_row);

  pipeline_mng->getOpDomainRhsPipeline().push_back(
      new OpDomainSource(FIELD_NAME, approxFunction));

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
  pipeline_mng->getSkeletonRhsFE().reset();

  auto rule = [](int, int, int p) -> int { return 2 * p + 1; };
  CHKERR pipeline_mng->setDomainRhsIntegrationRule(rule);
  CHKERR pipeline_mng->setSkeletonRhsIntegrationRule(rule);
  pipeline_mng->getDomainRhsFE()->exeTestHook = test_bit_child;
  pipeline_mng->getSkeletonRhsFE()->exeTestHook = test_bit_child;

  auto common_data_ptr = boost::make_shared<CommonData>();
  common_data_ptr->resVec = smartCreateDMVector(simpleInterface->getDM());
  common_data_ptr->L2Vec = createSmartVectorMPI(
      mField.get_comm(), (!mField.get_comm_rank()) ? 1 : 0, 1);
  common_data_ptr->approxVals = boost::make_shared<VectorDouble>();
  common_data_ptr->divApproxVals = boost::make_shared<MatrixDouble>();

  auto jac_ptr = boost::make_shared<MatrixDouble>();
  auto inv_jac_ptr = boost::make_shared<MatrixDouble>();
  auto det_ptr = boost::make_shared<VectorDouble>();

  pipeline_mng->getOpDomainRhsPipeline().push_back(
      new OpCalculateHOJacForFace(jac_ptr));
  pipeline_mng->getOpDomainRhsPipeline().push_back(
      new OpInvertMatrix<2>(jac_ptr, det_ptr, inv_jac_ptr));
  pipeline_mng->getOpDomainRhsPipeline().push_back(
      new OpSetInvJacH1ForFace(inv_jac_ptr));

  set_parent_dofs<DomainParentEle>(mField, pipeline_mng->getDomainRhsFE(),
                                   DomainEleOp::OPSPACE, QUIET, Sev::noisy);
  set_parent_dofs<DomainParentEle>(mField, pipeline_mng->getDomainRhsFE(),
                                   DomainEleOp::OPROW, VERBOSE, Sev::noisy);

  pipeline_mng->getOpDomainRhsPipeline().push_back(
      new OpCalculateScalarFieldGradient<SPACE_DIM>(
          FIELD_NAME, common_data_ptr->divApproxVals));
  pipeline_mng->getOpDomainRhsPipeline().push_back(
      new OpCalculateScalarFieldValues(FIELD_NAME,
                                       common_data_ptr->approxVals));

  pipeline_mng->getOpDomainRhsPipeline().push_back(
      new OpError<FIELD_DIM>(common_data_ptr));

  set_parent_dofs<SkeletonParentEle>(mField, pipeline_mng->getSkeletonRhsFE(),
                                     SkeletonEleOp::OPSPACE, QUIET, Sev::noisy);
  set_parent_dofs<SkeletonParentEle>(mField, pipeline_mng->getSkeletonRhsFE(),
                                     SkeletonEleOp::OPROW, VERBOSE, Sev::noisy);
  pipeline_mng->getOpSkeletonRhsPipeline().push_back(
      new OpCalculateScalarFieldValues(FIELD_NAME,
                                       common_data_ptr->approxVals));
  pipeline_mng->getOpSkeletonRhsPipeline().push_back(
      new OpErrorSkel<FIELD_DIM>(common_data_ptr));

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

  constexpr double eps = 1e-8;
  if (nrm2 > eps)
    SETERRQ1(PETSC_COMM_SELF, MOFEM_ATOM_TEST_INVALID,
             "Not converged solution err = %6.4e", nrm2);
  if (std::sqrt(array[0]) > eps)
    SETERRQ1(PETSC_COMM_SELF, MOFEM_ATOM_TEST_INVALID,
             "Error in approximation err = %6.4e", std::sqrt(array[0]));

  CHKERR VecRestoreArrayRead(common_data_ptr->L2Vec, &array);

  MoFEMFunctionReturn(0);
}
//! [Check results]

MoFEMErrorCode AtomTest::printResults() {
  MoFEMFunctionBegin;
  PipelineManager *pipeline_mng = mField.getInterface<PipelineManager>();
  pipeline_mng->getDomainLhsFE().reset();
  pipeline_mng->getDomainRhsFE().reset();

  auto rule = [](int, int, int p) -> int { return -1; };
  CHKERR pipeline_mng->setDomainRhsIntegrationRule(rule);

  static_cast<ForcesAndSourcesCore *>(pipeline_mng->getDomainRhsFE().get())
      ->setRuleHook = [](ForcesAndSourcesCore *fe_raw_ptr, int order_row,
                         int order_col, int order_data) -> MoFEMErrorCode {
    MoFEMFunctionBeginHot;
    fe_raw_ptr->gaussPts.resize(3, 3);
    fe_raw_ptr->gaussPts(0, 0) = 0;
    fe_raw_ptr->gaussPts(1, 0) = 0;
    fe_raw_ptr->gaussPts(2, 0) = 0;
    fe_raw_ptr->gaussPts(0, 1) = 1;
    fe_raw_ptr->gaussPts(1, 1) = 0;
    fe_raw_ptr->gaussPts(2, 1) = 0;
    fe_raw_ptr->gaussPts(0, 2) = 0;
    fe_raw_ptr->gaussPts(1, 2) = 1;
    fe_raw_ptr->gaussPts(2, 2) = 0;
    MoFEMFunctionReturnHot(0);
  };

  auto field_op_row = new ForcesAndSourcesCore::UserDataOperator(
      FIELD_NAME, DomainEleOp::OPROW);

  auto approx_vals = boost::make_shared<VectorDouble>();

  auto &moab = mField.get_moab();
  Tag th;
  double def_val[] = {0};
  CHKERR moab.tag_get_handle("FIELD", 1, MB_TYPE_DOUBLE, th,
                             MB_TAG_CREAT | MB_TAG_SPARSE, &def_val);

  field_op_row->doWorkRhsHook = [&](DataOperator *base_op_ptr, int side,
                                    EntityType type,
                                    EntitiesFieldData::EntData &data) {
    MoFEMFunctionBegin;
    if (type == MBVERTEX) {
      auto op_ptr =
          static_cast<FaceElementForcesAndSourcesCore::UserDataOperator *>(
              base_op_ptr);
      auto t_field = getFTensor0FromVec(*approx_vals);
      auto nb_gauss_pts = op_ptr->getGaussPts().size2();
      if (nb_gauss_pts != 3)
        SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
                "Should be three guass pts.");
      auto conn = op_ptr->getConn();
      for (auto gg = 0; gg != nb_gauss_pts; ++gg) {
        const double v = t_field;
        CHKERR moab.tag_set_data(th, &conn[gg], 1, &v);
        ++t_field;
      }
    }
    MoFEMFunctionReturn(0);
  };

  set_parent_dofs<DomainParentEle>(mField, pipeline_mng->getDomainRhsFE(),
                                   DomainEleOp::OPSPACE, VERBOSE, Sev::noisy);
  set_parent_dofs<DomainParentEle>(mField, pipeline_mng->getDomainRhsFE(),
                                   DomainEleOp::OPROW, VERBOSE, Sev::noisy);
  pipeline_mng->getOpDomainRhsPipeline().push_back(
      new OpCalculateScalarFieldValues(FIELD_NAME, approx_vals));
  pipeline_mng->getOpDomainRhsPipeline().push_back(field_op_row);
  CHKERR pipeline_mng->loopFiniteElements();

  CHKERR mField.getInterface<BitRefManager>()->writeBitLevelByType(
      bit(nb_ref_levels), BitRefLevel().set(), MBTRI, "out.vtk", "VTK", "");

  MoFEMFunctionReturn(0);
}

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
