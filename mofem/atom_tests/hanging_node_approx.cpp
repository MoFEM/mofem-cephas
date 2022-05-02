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
  MoFEMErrorCode checkResults();

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

  auto refine_mesh = [&](auto bit0, auto bit1) {
    MoFEMFunctionBegin;

    auto refine = mField.getInterface<MeshRefinement>();

    auto meshset_level0_ptr = get_temp_meshset_ptr(moab);
    CHKERR mField.getInterface<BitRefManager>()->getEntitiesByRefLevel(
        bit0, BitRefLevel().set(), *meshset_level0_ptr);

    // random mesh refinement
    auto meshset_ref_edges_ptr = get_temp_meshset_ptr(moab);

    Range eles_to_refine;
    CHKERR moab.get_entities_by_dimension(
        *meshset_level0_ptr, simpleInterface->getDim(), eles_to_refine);
    int ii = 0;
    for (auto t : eles_to_refine) {
      if (ii % 2) {
        std::vector<EntityHandle> adj_edges;
        CHKERR mField.get_moab().get_adjacencies(&t, 1, 1, false, adj_edges);
        CHKERR moab.add_entities(*meshset_ref_edges_ptr, &*adj_edges.begin(),
                                 adj_edges.size());
      }
      ++ii;
    }

    CHKERR refine->addVerticesInTheMiddleOfEdges(*meshset_ref_edges_ptr, bit1,
                                                 false, VERBOSE);

    if (simpleInterface->getDim() == 3) {
      CHKERR refine->refineTetsHangingNodes(*meshset_level0_ptr, bit1, VERBOSE);
    } else if (simpleInterface->getDim() == 2) {
      CHKERR refine->refineTrisHangingNodes(*meshset_level0_ptr, bit1, VERBOSE);
    } else {
      SETERRQ(PETSC_COMM_WORLD, MOFEM_ATOM_TEST_INVALID,
              "Dimension not handled by test");
    }

    MoFEMFunctionReturn(0);
  };

  BitRefLevel bit_level1;
  bit_level1.set(1);
  CHKERR refine_mesh(bit_level0, bit_level1);
  simpleInterface->getBitRefLevel() = BitRefLevel().set();
  simpleInterface->getBitRefLevelMask() = BitRefLevel().set();

  BitRefManager *bit_mng = mField.getInterface<BitRefManager>();
  CHKERR bit_mng->writeBitLevelByType(bit_level1, BitRefLevel().set(), MBTRI,
                                      "out_ref_mesh.vtk", "VTK", "");

  auto bit_mark = BitRefLevel().set(2);
  CHKERR OpAddParentEntData::markHangingSkinParents(
      mField, 2, bit_level0, BitRefLevel().set(), bit_level1,
      BitRefLevel().set(), bit_mark, true, "parent_ref_skin.vtk");
  CHKERR OpAddParentEntData::markHangingSkinChildren(
      mField, bit_level1, bit_level1, bit_mark, BitRefLevel().set(),
      "children_ref_skin.vtk");

  MoFEMFunctionReturn(0);
}
//! [Read mesh]

//! [Set up problem]
MoFEMErrorCode AtomTest::setupProblem() {
  MoFEMFunctionBegin;
  // Add field
  CHKERR simpleInterface->addDomainField(FIELD_NAME, H1,
                                         AINSWORTH_LEGENDRE_BASE, FIELD_DIM);

  ElementAdjacencyFunct get_parent_adjacency =
      [&](moab::Interface &moab, const Field &field, const EntFiniteElement &fe,
          std::vector<EntityHandle> &adjacency) -> MoFEMErrorCode {
    MoFEMFunctionBegin;

    CHKERR DefaultElementAdjacency::defaultFace(moab, field, fe, adjacency);

    auto basic_entity_data_ptr = fe.getBasicDataPtr();
    auto th_parent_handle = basic_entity_data_ptr->th_RefParentHandle;

    using GetParent = boost::function<MoFEMErrorCode(
        EntityHandle fe, std::vector<EntityHandle> & parents)>;
    GetParent get_parent = [&](EntityHandle fe,
                               std::vector<EntityHandle> &parents) {
      MoFEMFunctionBegin;
      EntityHandle fe_parent;

      CHKERR moab.tag_get_data(th_parent_handle, &fe, 1, &fe_parent);
      auto parent_type = type_from_handle(fe_parent);
      auto back_type = type_from_handle(fe);
      if (fe_parent != 0 && fe_parent != fe && parent_type == back_type) {
        parents.push_back(fe_parent);
        CHKERR get_parent(parents.back(), parents);
      }
      MoFEMFunctionReturn(0);
    };

    if (field.getSpace() != NOFIELD) {

      const auto fe_name = fe.getName();
      const auto parent_ent = fe.getParentEnt();

      std::vector<EntityHandle> parents;
      parents.reserve(BITREFLEVEL_SIZE);

      if (parent_ent && parent_ent != fe.getEnt()) {
        CHKERR get_parent(fe.getEnt(), parents);
        for (auto fe_ent : parents) {
          switch (field.getSpace()) {
          case H1:
            CHKERR moab.get_adjacencies(&fe_ent, 1, 0, false, adjacency,
                                        moab::Interface::UNION);
          case HCURL:
            CHKERR moab.get_adjacencies(&fe_ent, 1, 1, false, adjacency,
                                        moab::Interface::UNION);
          case HDIV:
          case L2:
            adjacency.push_back(fe_ent);
            break;
          default:
            SETERRQ(PETSC_COMM_SELF, MOFEM_NOT_IMPLEMENTED,
                    "this field is not implemented for face finite element");
          }
        }

        std::sort(adjacency.begin(), adjacency.end());
        auto it = std::unique(adjacency.begin(), adjacency.end());
        adjacency.resize(std::distance(adjacency.begin(), it));

        for (auto e : adjacency) {
          auto side_table = fe.getSideNumberTable();
          if (side_table.find(e) == side_table.end())
            const_cast<SideNumber_multiIndex &>(fe.getSideNumberTable())
                .insert(
                    boost::shared_ptr<SideNumber>(new SideNumber(e, -1, 0, 0)));
        }
      }
    }

    MoFEMFunctionReturn(0);
  };

  constexpr int order = 4;
  CHKERR simpleInterface->setFieldOrder(FIELD_NAME, order);
  // CHKERR simpleInterface->setUp();

  CHKERR simpleInterface->defineFiniteElements();

  CHKERR mField.modify_finite_element_adjacency_table(
      simpleInterface->getDomainFEName(), MBTRI, get_parent_adjacency);
  CHKERR mField.modify_finite_element_adjacency_table(
      simpleInterface->getDomainFEName(), MBQUAD, get_parent_adjacency);

  CHKERR simpleInterface->defineProblem(PETSC_TRUE);
  CHKERR simpleInterface->buildFields();
  CHKERR simpleInterface->buildFiniteElements();
  CHKERR simpleInterface->buildProblem();

  auto bit_l0 = BitRefLevel().set(0);
  auto bit_l1 = BitRefLevel().set(1);
  auto bit_mark = BitRefLevel().set(2);

  BitRefManager *bit_mng = mField.getInterface<BitRefManager>();
  ProblemsManager *prb_mng = mField.getInterface<ProblemsManager>();

  CHKERR bit_mng->writeBitLevelByType(bit_mark, bit_l0 | bit_l1 | bit_mark,
                                      MBEDGE, "l0_ents_edges.vtk", "VTK", "");
  CHKERR bit_mng->writeBitLevelByType(bit_mark, bit_l1 | bit_mark, MBEDGE,
                                      "l1_ents_edges.vtk", "VTK", "");
  CHKERR bit_mng->writeBitLevelByType(bit_mark, bit_l1 | bit_mark, MBVERTEX,
                                      "l1_ents_verts.vtk", "VTK", "");

  CHKERR prb_mng->removeDofsOnEntities(simpleInterface->getProblemName(),
                                       FIELD_NAME, bit_l0, bit_l0);
  CHKERR prb_mng->removeDofsOnEntities(simpleInterface->getProblemName(),
                                       FIELD_NAME, bit_mark, bit_l1 | bit_mark);

  MoFEMFunctionReturn(0);
}
//! [Set up problem]

//! [Push operators to pipeline]
MoFEMErrorCode AtomTest::assembleSystem() {
  MoFEMFunctionBegin;
  PipelineManager *pipeline_mng = mField.getInterface<PipelineManager>();

  auto l0 = BitRefLevel().set(0);
  auto l1 = BitRefLevel().set(1);
  auto marker = BitRefLevel().set(2);

  auto rule = [](int, int, int p) -> int { return 2 * p; };

  CHKERR pipeline_mng->setDomainLhsIntegrationRule(rule);
  CHKERR pipeline_mng->setDomainRhsIntegrationRule(rule);

  auto test_bit_child = [](FEMethod *fe_ptr) {
    return fe_ptr->numeredEntFiniteElementPtr->getBitRefLevel().test(1);
  };

  pipeline_mng->getDomainLhsFE()->exeTestHook = test_bit_child;
  pipeline_mng->getDomainRhsFE()->exeTestHook = test_bit_child;

  auto parent_fe_ptr = boost::make_shared<DomainParentEle>(mField);

  auto set_parent_dofs = [&](auto &pipeline, auto op, auto verbosity,
                             auto sev) {
    pipeline.push_back(

        new OpAddParentEntData(

            FIELD_NAME, op, parent_fe_ptr,

            // level 1 elements
            l1, l1,

            // marked level 1
            marker, l0 | l1 | marker,

            verbosity, sev)

    );
  };

  auto beta = [](const double, const double, const double) { return 1; };
  set_parent_dofs(pipeline_mng->getOpDomainLhsPipeline(), DomainEleOp::OPROW,
                  QUIET, Sev::noisy);
  set_parent_dofs(pipeline_mng->getOpDomainLhsPipeline(), DomainEleOp::OPCOL,
                  QUIET, Sev::noisy);
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
    }
    MoFEMFunctionReturn(0);
  };

  auto field_op_col = new ForcesAndSourcesCore::UserDataOperator(
      FIELD_NAME, DomainEleOp::OPCOL);
  field_op_col->doWorkRhsHook = [](DataOperator *op_ptr, int side,
                                   EntityType type,
                                   EntitiesFieldData::EntData &data) {
    MoFEMFunctionBegin;
    if (type == MBENTITYSET) {
      MOFEM_LOG("SELF", Sev::verbose)
          << "COL: side/type: " << side << "/" << CN::EntityTypeName(type)
          << " op space/base: " << FieldSpaceNames[data.getSpace()] << "/"
          << ApproximationBaseNames[data.getBase()] << " DOFs "
          << data.getIndices() << " nb base functions " << data.getN().size2()
          << " nb base functions integration points " << data.getN().size1();
    }
    MoFEMFunctionReturn(0);
  };

  set_parent_dofs(pipeline_mng->getOpDomainRhsPipeline(), DomainEleOp::OPROW,
                  VERBOSE, Sev::noisy);
  // set_parent_dofs(pipeline_mng->getOpDomainRhsPipeline(), DomainEleOp::OPCOL,
  //                 VERBOSE, Sev::noisy);
  // pipeline_mng->getOpDomainRhsPipeline().push_back(field_op_row);
  // pipeline_mng->getOpDomainRhsPipeline().push_back(field_op_col);

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

  auto rule = [](int, int, int p) -> int { return 2 * p + 1; };
  CHKERR pipeline_mng->setDomainRhsIntegrationRule(rule);

  auto test_bit_child = [](FEMethod *fe_ptr) {
    return fe_ptr->numeredEntFiniteElementPtr->getBitRefLevel().test(1);
  };
  pipeline_mng->getDomainRhsFE()->exeTestHook = test_bit_child;

  auto common_data_ptr = boost::make_shared<CommonData>();
  common_data_ptr->resVec = smartCreateDMVector(simpleInterface->getDM());
  common_data_ptr->L2Vec = createSmartVectorMPI(
      mField.get_comm(), (!mField.get_comm_rank()) ? 1 : 0, 1);
  common_data_ptr->approxVals = boost::make_shared<VectorDouble>();


  auto parent_fe_ptr = boost::make_shared<DomainParentEle>(mField);

  auto l0 = BitRefLevel().set(0);
  auto l1 = BitRefLevel().set(1);
  auto marker = BitRefLevel().set(2);

  auto set_parent_dofs = [&](auto &pipeline, auto op, auto verbosity,
                             auto sev) {
    pipeline.push_back(

        new OpAddParentEntData(

            FIELD_NAME, op, parent_fe_ptr,

            // level 1 elements
            l1, l1,

            // marked level 1
            marker, l0 | l1 | marker,

            verbosity, sev)

    );
  };

  set_parent_dofs(pipeline_mng->getOpDomainRhsPipeline(), DomainEleOp::OPROW,
                  VERBOSE, Sev::noisy);
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
