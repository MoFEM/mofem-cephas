#include <stdlib.h>
#include <MoFEM.hpp>
#include <nonlinear_poisson_2d.hpp>

using namespace MoFEM;
using namespace NonlinearPoissonOps;

constexpr int SPACE_DIM = 2;
static char help[] = "...\n\n";

inline double sqr(const double x) { return x * x; }

inline double cube(const double x) { return x * x * x; }

struct NonlinearPoisson {
public:
  NonlinearPoisson(MoFEM::Interface &m_field);

  // Declaration of the main function to run analysis
  MoFEMErrorCode runProgram();

private:
  // Declaration of other main functions called in runProgram()
  MoFEMErrorCode readMesh();
  MoFEMErrorCode setupProblem();
  MoFEMErrorCode setIntegrationRules();
  MoFEMErrorCode boundaryCondition();
  MoFEMErrorCode assembleSystem();
  MoFEMErrorCode solveSystem();
  MoFEMErrorCode outputResults();
  MoFEMErrorCode checkResults();

  // Function to calculate the Source term
  static double sourceTermFunction(const double x, const double y,
                                   const double z) {

    return 2 * M_PI * M_PI *
           (cos(M_PI * x) * cos(M_PI * y) +
            cube(cos(M_PI * x)) * cube(cos(M_PI * y)) -
            cos(M_PI * x) * cos(M_PI * y) *
                (sqr(sin(M_PI * x)) * sqr(cos(M_PI * y)) +
                 sqr(sin(M_PI * y)) * sqr(cos(M_PI * x))));
  }

  // Function to calculate the Boundary term
  static double boundaryFunction(const double x, const double y,
                                 const double z) {
    return -cos(M_PI * x) *
           cos(M_PI * y); // here should put the negative of the proper formula
  }

  // Main interfaces
  MoFEM::Interface &mField;
  Simple *simpleInterface;

  // Field name and approximation order
  std::string domainField = "POTENTIAL";
  int order;

  // PETSc vector for storing norms
  SmartPetscObj<Vec> errorVec, exactVec;
  int atomTest = 0;
  enum NORMS { NORM = 0, LAST_NORM };
  enum EX_NORMS { EX_NORM = 0, LAST_EX_NORM };

  // Object to mark boundary entities for the assembling of domain elements
  boost::shared_ptr<std::vector<unsigned char>> boundaryMarker;

  using PostProcEle = PostProcBrokenMeshInMoab<DomainEle>;

  // Object needed for postprocessing
  boost::shared_ptr<PostProcEle> postProc;

  // Boundary entities marked for fieldsplit (block) solver - optional
  Range boundaryEntitiesForFieldsplit;
};

NonlinearPoisson::NonlinearPoisson(MoFEM::Interface &m_field)
    : mField(m_field) {}

MoFEMErrorCode NonlinearPoisson::runProgram() {
  MoFEMFunctionBegin;

  CHKERR readMesh();
  CHKERR setupProblem();
  CHKERR setIntegrationRules();
  CHKERR boundaryCondition();
  CHKERR assembleSystem();
  CHKERR solveSystem();
  CHKERR outputResults();
  CHKERR checkResults();

  MoFEMFunctionReturn(0);
}

MoFEMErrorCode NonlinearPoisson::readMesh() {
  MoFEMFunctionBegin;

  CHKERR mField.getInterface(simpleInterface);
  CHKERR simpleInterface->getOptions();
  CHKERR simpleInterface->loadFile();

  MoFEMFunctionReturn(0);
}

MoFEMErrorCode NonlinearPoisson::setupProblem() {
  MoFEMFunctionBegin;

  Range domain_ents;
  CHKERR mField.get_moab().get_entities_by_dimension(0, SPACE_DIM, domain_ents,
                                                     true);
  auto get_ents_by_dim = [&](const auto dim) {
    if (dim == SPACE_DIM) {
      return domain_ents;
    } else {
      Range ents;
      if (dim == 0)
        CHKERR mField.get_moab().get_connectivity(domain_ents, ents, true);
      else
        CHKERR mField.get_moab().get_entities_by_dimension(0, dim, ents, true);
      return ents;
    }
  };
  // Select base
  auto get_base = [&]() {
    auto domain_ents = get_ents_by_dim(SPACE_DIM);
    if (domain_ents.empty())
      CHK_THROW_MESSAGE(MOFEM_NOT_FOUND, "Empty mesh");
    const auto type = type_from_handle(domain_ents[0]);
    switch (type) {
    case MBQUAD:
      return DEMKOWICZ_JACOBI_BASE;
    case MBHEX:
      return DEMKOWICZ_JACOBI_BASE;
    case MBTRI:
      return AINSWORTH_LEGENDRE_BASE;
    case MBTET:
      return AINSWORTH_LEGENDRE_BASE;
    default:
      CHK_THROW_MESSAGE(MOFEM_NOT_FOUND, "Element type not handled");
    }
    return NOBASE;
  };
  auto base = get_base();
  CHKERR simpleInterface->addDomainField(domainField, H1, base, 1);
  CHKERR simpleInterface->addBoundaryField(domainField, H1, base, 1);

  order = 2;

  CHKERR PetscOptionsGetInt(PETSC_NULL, "", "-order", &order, PETSC_NULL);
  CHKERR PetscOptionsGetInt(PETSC_NULL, "", "-atom_test", &atomTest,
                            PETSC_NULL);
  CHKERR simpleInterface->setFieldOrder(domainField, order);

  CHKERR simpleInterface->setUp();

  MoFEMFunctionReturn(0);
}

MoFEMErrorCode NonlinearPoisson::setIntegrationRules() {
  MoFEMFunctionBegin;

  auto domain_rule_lhs = [](int, int, int p) -> int { return 2 * p - 1; };
  auto domain_rule_rhs = [](int, int, int p) -> int { return 2 * p - 1; };

  auto boundary_rule_lhs = [](int, int, int p) -> int { return 2 * p; };
  auto boundary_rule_rhs = [](int, int, int p) -> int { return 2 * p; };

  auto pipeline_mng = mField.getInterface<PipelineManager>();
  CHKERR pipeline_mng->setDomainRhsIntegrationRule(domain_rule_lhs);
  CHKERR pipeline_mng->setDomainLhsIntegrationRule(domain_rule_rhs);
  CHKERR pipeline_mng->setBoundaryLhsIntegrationRule(boundary_rule_lhs);
  CHKERR pipeline_mng->setBoundaryRhsIntegrationRule(boundary_rule_rhs);

  MoFEMFunctionReturn(0);
}

MoFEMErrorCode NonlinearPoisson::boundaryCondition() {
  MoFEMFunctionBegin;

  // Get boundary edges marked in block named "BOUNDARY_CONDITION"
  auto get_ents_on_mesh_skin = [&]() {
    Range boundary_entities;
    for (_IT_CUBITMESHSETS_BY_SET_TYPE_FOR_LOOP_(mField, BLOCKSET, it)) {
      std::string entity_name = it->getName();
      if (entity_name.compare(0, 18, "BOUNDARY_CONDITION") == 0) {
        CHKERR it->getMeshsetIdEntitiesByDimension(mField.get_moab(), 1,
                                                   boundary_entities, true);
      }
    }
    // Add vertices to boundary entities
    Range boundary_vertices;
    CHKERR mField.get_moab().get_connectivity(boundary_entities,
                                              boundary_vertices, true);
    boundary_entities.merge(boundary_vertices);

    // Store entities for fieldsplit (block) solver
    boundaryEntitiesForFieldsplit = boundary_entities;

    return boundary_entities;
  };

  auto mark_boundary_dofs = [&](Range &&skin_edges) {
    auto problem_manager = mField.getInterface<ProblemsManager>();
    auto marker_ptr = boost::make_shared<std::vector<unsigned char>>();
    problem_manager->markDofs(simpleInterface->getProblemName(), ROW,
                              ProblemsManager::OR, skin_edges, *marker_ptr);
    return marker_ptr;
  };

  // Get global local vector of marked DOFs. Is global, since is set for all
  // DOFs on processor. Is local since only DOFs on processor are in the
  // vector. To access DOFs use local indices.
  boundaryMarker = mark_boundary_dofs(get_ents_on_mesh_skin());

  MoFEMFunctionReturn(0);
}

MoFEMErrorCode NonlinearPoisson::assembleSystem() {
  MoFEMFunctionBegin;

  auto pipeline_mng = mField.getInterface<PipelineManager>();
  CHKERR AddHOOps<2, 2, 2>::add(pipeline_mng->getOpDomainLhsPipeline(), {H1});
  CHKERR AddHOOps<2, 2, 2>::add(pipeline_mng->getOpDomainRhsPipeline(), {H1});

  auto add_domain_lhs_ops = [&](auto &pipeline) {
    pipeline.push_back(new OpSetBc(domainField, true, boundaryMarker));
    auto data_u_at_gauss_pts = boost::make_shared<VectorDouble>();
    auto grad_u_at_gauss_pts = boost::make_shared<MatrixDouble>();
    pipeline.push_back(
        new OpCalculateScalarFieldValues(domainField, data_u_at_gauss_pts));
    pipeline.push_back(new OpCalculateScalarFieldGradient<2>(
        domainField, grad_u_at_gauss_pts));
    pipeline.push_back(new OpDomainLhs(
        domainField, domainField, data_u_at_gauss_pts, grad_u_at_gauss_pts));
    pipeline.push_back(new OpUnSetBc(domainField));
  };

  auto add_domain_rhs_ops = [&](auto &pipeline) {
    pipeline.push_back(new OpSetBc(domainField, true, boundaryMarker));
    auto data_u_at_gauss_pts = boost::make_shared<VectorDouble>();
    auto grad_u_at_gauss_pts = boost::make_shared<MatrixDouble>();
    pipeline.push_back(
        new OpCalculateScalarFieldValues(domainField, data_u_at_gauss_pts));
    pipeline.push_back(new OpCalculateScalarFieldGradient<2>(
        domainField, grad_u_at_gauss_pts));
    pipeline.push_back(new OpDomainRhs(domainField, sourceTermFunction,
                                       data_u_at_gauss_pts,
                                       grad_u_at_gauss_pts));
    pipeline.push_back(new OpUnSetBc(domainField));
  };

  auto add_boundary_lhs_ops = [&](auto &pipeline) {
    pipeline.push_back(new OpSetBc(domainField, false, boundaryMarker));
    pipeline.push_back(new OpBoundaryLhs(
        domainField, domainField,
        [](const double, const double, const double) { return 1; }));
    pipeline.push_back(new OpUnSetBc(domainField));
  };

  auto add_boundary_rhs_ops = [&](auto &pipeline) {
    pipeline.push_back(new OpSetBc(domainField, false, boundaryMarker));
    auto u_at_gauss_pts = boost::make_shared<VectorDouble>();
    pipeline.push_back(
        new OpCalculateScalarFieldValues(domainField, u_at_gauss_pts));
    pipeline.push_back(new OpBoundaryRhs(
        domainField, u_at_gauss_pts,
        [](const double, const double, const double) { return 1; }));
    pipeline.push_back(new OpBoundaryRhsSource(domainField, boundaryFunction));
    pipeline.push_back(new OpUnSetBc(domainField));
  };

  add_domain_lhs_ops(pipeline_mng->getOpDomainLhsPipeline());
  add_domain_rhs_ops(pipeline_mng->getOpDomainRhsPipeline());

  add_boundary_lhs_ops(pipeline_mng->getOpBoundaryLhsPipeline());
  add_boundary_rhs_ops(pipeline_mng->getOpBoundaryRhsPipeline());

  MoFEMFunctionReturn(0);
}

MoFEMErrorCode NonlinearPoisson::solveSystem() {
  MoFEMFunctionBegin;

  auto *simple = mField.getInterface<Simple>();

  auto set_fieldsplit_preconditioner = [&](auto snes) {
    MoFEMFunctionBeginHot;
    KSP ksp;
    CHKERR SNESGetKSP(snes, &ksp);
    PC pc;
    CHKERR KSPGetPC(ksp, &pc);
    PetscBool is_pcfs = PETSC_FALSE;
    PetscObjectTypeCompare((PetscObject)pc, PCFIELDSPLIT, &is_pcfs);

    // Set up FIELDSPLIT, only when option used -pc_type fieldsplit
    if (is_pcfs == PETSC_TRUE) {
      auto name_prb = simple->getProblemName();
      SmartPetscObj<IS> is_all_bc;
      CHKERR mField.getInterface<ISManager>()->isCreateProblemFieldAndRank(
          name_prb, ROW, domainField, 0, 1, is_all_bc,
          &boundaryEntitiesForFieldsplit);
      int is_all_bc_size;
      CHKERR ISGetSize(is_all_bc, &is_all_bc_size);
      MOFEM_LOG("EXAMPLE", Sev::inform)
          << "Field split block size " << is_all_bc_size;
      CHKERR PCFieldSplitSetIS(pc, PETSC_NULL,
                               is_all_bc); // boundary block
    }
    MoFEMFunctionReturnHot(0);
  };

  // Create global RHS and solution vectors
  auto dm = simple->getDM();
  SmartPetscObj<Vec> global_rhs, global_solution;
  CHKERR DMCreateGlobalVector_MoFEM(dm, global_rhs);
  global_solution = vectorDuplicate(global_rhs);

  // Create nonlinear solver (SNES)
  auto pipeline_mng = mField.getInterface<PipelineManager>();
  auto solver = pipeline_mng->createSNES();
  CHKERR SNESSetFromOptions(solver);
  CHKERR set_fieldsplit_preconditioner(solver);
  CHKERR SNESSetUp(solver);

  // Solve the system
  CHKERR SNESSolve(solver, global_rhs, global_solution);

  // Scatter result data on the mesh
  CHKERR DMoFEMMeshToGlobalVector(dm, global_solution, INSERT_VALUES,
                                  SCATTER_REVERSE);

  MoFEMFunctionReturn(0);
}

MoFEMErrorCode NonlinearPoisson::outputResults() {
  MoFEMFunctionBegin;

  auto post_proc = boost::make_shared<PostProcEle>(mField);

  auto u_ptr = boost::make_shared<VectorDouble>();
  post_proc->getOpPtrVector().push_back(
      new OpCalculateScalarFieldValues(domainField, u_ptr));

  using OpPPMap = OpPostProcMapInMoab<2, 2>;

  post_proc->getOpPtrVector().push_back(

      new OpPPMap(post_proc->getPostProcMesh(), post_proc->getMapGaussPts(),

                  {{domainField, u_ptr}},

                  {}, {}, {}

                  )

  );

  auto *simple = mField.getInterface<Simple>();
  auto dm = simple->getDM();
  CHKERR DMoFEMLoopFiniteElements(dm, simple->getDomainFEName(), post_proc);

  CHKERR post_proc->writeFile("out_result.h5m");

  MoFEMFunctionReturn(0);
}

//! [Check]
MoFEMErrorCode NonlinearPoisson::checkResults() {
  MoFEMFunctionBegin;
  auto check_result_fe_ptr = boost::make_shared<DomainEle>(mField);
  auto errorVec =
      createVectorMPI(mField.get_comm(),
                      (mField.get_comm_rank() == 0) ? LAST_NORM : 0, LAST_NORM);
  auto exactVec = createVectorMPI(
      mField.get_comm(), (mField.get_comm_rank() == 0) ? LAST_EX_NORM : 0,
      LAST_EX_NORM);
  CHK_THROW_MESSAGE((AddHOOps<SPACE_DIM, SPACE_DIM, SPACE_DIM>::add(
                        check_result_fe_ptr->getOpPtrVector(), {H1})),
                    "Apply transform");
  check_result_fe_ptr->getRuleHook = [](int, int, int p) { return p; };
  auto analyticalFunction = [&](double x, double y, double z) {
    return cos(M_PI * x) * cos(M_PI * y);
  };
  auto u_ptr = boost::make_shared<VectorDouble>();
  check_result_fe_ptr->getOpPtrVector().push_back(
      new OpCalculateScalarFieldValues(domainField, u_ptr));
  auto mValFuncPtr = boost::make_shared<VectorDouble>();
  check_result_fe_ptr->getOpPtrVector().push_back(
      new OpGetTensor0fromFunc(mValFuncPtr, analyticalFunction));
  check_result_fe_ptr->getOpPtrVector().push_back(
      new OpCalcNormL2Tensor0(u_ptr, errorVec, NORM, mValFuncPtr));
  check_result_fe_ptr->getOpPtrVector().push_back(
      new OpCalcNormL2Tensor0(u_ptr, exactVec, EX_NORM));
  CHKERR VecZeroEntries(errorVec);
  CHKERR VecZeroEntries(exactVec);
  CHKERR DMoFEMLoopFiniteElements(simpleInterface->getDM(),
                                  simpleInterface->getDomainFEName(),
                                  check_result_fe_ptr);
  CHKERR VecAssemblyBegin(errorVec);
  CHKERR VecAssemblyEnd(errorVec);
  CHKERR VecAssemblyBegin(exactVec);
  CHKERR VecAssemblyEnd(exactVec);
  MOFEM_LOG_CHANNEL("SELF"); // Clear channel from old tags
  // print norm in general log
  if (mField.get_comm_rank() == 0) {
    const double *L2_norms;
    const double *Ex_norms;
    CHKERR VecGetArrayRead(errorVec, &L2_norms);
    CHKERR VecGetArrayRead(exactVec, &Ex_norms);
    MOFEM_TAG_AND_LOG("SELF", Sev::inform, "Errors")
        << "L2_NORM: " << std::sqrt(L2_norms[NORM]);
    MOFEM_TAG_AND_LOG("SELF", Sev::inform, "Errors")
        << "NORMALISED_ERROR: "
        << (std::sqrt(L2_norms[NORM]) / std::sqrt(Ex_norms[EX_NORM]));
    CHKERR VecRestoreArrayRead(errorVec, &L2_norms);
    CHKERR VecRestoreArrayRead(exactVec, &Ex_norms);
  }
  // compare norm for ctest
  if (atomTest && !mField.get_comm_rank()) {
    const double *L2_norms;
    const double *Ex_norms;
    CHKERR VecGetArrayRead(errorVec, &L2_norms);
    CHKERR VecGetArrayRead(exactVec, &Ex_norms);
    double ref_percentage = 4.4e-04;
    double normalised_error;
    switch (atomTest) {
    case 1: // 2D
      normalised_error = std::sqrt(L2_norms[0]) / std::sqrt(Ex_norms[0]);
      break;
    default:
      SETERRQ1(PETSC_COMM_SELF, MOFEM_INVALID_DATA,
               "atom test %d does not exist", atomTest);
    }
    if (normalised_error > ref_percentage) {
      SETERRQ3(PETSC_COMM_SELF, MOFEM_ATOM_TEST_INVALID,
               "atom test %d failed! Calculated Norm %3.16e is greater than "
               "reference Norm %3.16e",
               atomTest, normalised_error, ref_percentage);
    }
    CHKERR VecRestoreArrayRead(errorVec, &L2_norms);
    CHKERR VecRestoreArrayRead(exactVec, &Ex_norms);
  }
  MoFEMFunctionReturn(0);
}
//! [Check]

int main(int argc, char *argv[]) {

  // Initialisation of MoFEM/PETSc and MOAB data structures
  const char param_file[] = "param_file.petsc";
  MoFEM::Core::Initialize(&argc, &argv, param_file, help);

  auto core_log = logging::core::get();
  core_log->add_sink(
      LogManager::createSink(LogManager::getStrmWorld(), "EXAMPLE"));
  LogManager::setLog("EXAMPLE");
  MOFEM_LOG_TAG("EXAMPLE", "example")

  // Error handling
  try {
    // Register MoFEM discrete manager in PETSc
    DMType dm_name = "DMMOFEM";
    CHKERR DMRegister_MoFEM(dm_name);

    // Create MOAB instance
    moab::Core mb_instance;              // mesh database
    moab::Interface &moab = mb_instance; // mesh database interface

    // Create MoFEM instance
    MoFEM::Core core(moab);           // finite element database
    MoFEM::Interface &m_field = core; // finite element interface

    // Run the main analysis
    NonlinearPoisson poisson_problem(m_field);
    CHKERR poisson_problem.runProgram();
  }
  CATCH_ERRORS;

  // Finish work: cleaning memory, getting statistics, etc.
  MoFEM::Core::Finalize();

  return 0;
}