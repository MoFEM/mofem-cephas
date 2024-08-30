#include <stdlib.h>
#include <MoFEM.hpp>
#include <nonlinear_poisson_2d.hpp>

using namespace MoFEM;
using namespace NonlinearPoissonOps;

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

  // Discrete Manager and nonliner SNES solver using SmartPetscObj
  SmartPetscObj<DM> dM;
  SmartPetscObj<SNES> snesSolver;

  // Field name and approximation order
  std::string domainField = "POTENTIAL";
  int order;

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

  readMesh();
  setupProblem();
  setIntegrationRules();
  boundaryCondition();
  assembleSystem();
  solveSystem();
  outputResults();

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

  CHKERR simpleInterface->addDomainField(domainField, H1,
                                         AINSWORTH_BERNSTEIN_BEZIER_BASE, 1);
  CHKERR simpleInterface->addBoundaryField(domainField, H1,
                                           AINSWORTH_BERNSTEIN_BEZIER_BASE, 1);

  order = 4;

  CHKERR PetscOptionsGetInt(PETSC_NULL, "", "-order", &order, PETSC_NULL);
  CHKERR simpleInterface->setFieldOrder(domainField, order);

  CHKERR simpleInterface->setUp();

  MoFEMFunctionReturn(0);
}

MoFEMErrorCode NonlinearPoisson::setIntegrationRules() {
  MoFEMFunctionBegin;

  auto domain_rule_lhs = [](int, int, int p) -> int { return 2 * (p - 1); };
  auto domain_rule_rhs = [](int, int, int p) -> int { return 2 * (p - 1); };

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