#include <poisson_2d_nonhomogeneous.hpp>

using namespace MoFEM;
using namespace Poisson2DNonhomogeneousOperators;

constexpr int SPACE_DIM = 2;
using PostProcFaceEle =
    PostProcBrokenMeshInMoab<FaceElementForcesAndSourcesCore>;
static char help[] = "...\n\n";
struct Poisson2DNonhomogeneous {
public:
  Poisson2DNonhomogeneous(MoFEM::Interface &m_field);

  // Declaration of the main function to run analysis
  MoFEMErrorCode runProgram();

private:
  // Declaration of other main functions called in runProgram()
  MoFEMErrorCode readMesh();
  MoFEMErrorCode setupProblem();
  MoFEMErrorCode boundaryCondition();
  MoFEMErrorCode assembleSystem();
  MoFEMErrorCode setIntegrationRules();
  MoFEMErrorCode solveSystem();
  MoFEMErrorCode outputResults();
  MoFEMErrorCode checkResults();

  // Function to calculate the Source term
  static double sourceTermFunction(const double x, const double y,
                                   const double z) {
    return 2. * M_PI * M_PI * cos(x * M_PI) * cos(y * M_PI);
  }
  // Function to calculate the Boundary term
  static double boundaryFunction(const double x, const double y,
                                 const double z) {
    return cos(x * M_PI) * cos(y * M_PI);
    // return 0;
  }

  // Main interfaces
  MoFEM::Interface &mField;
  Simple *simpleInterface;

  // Field name and approximation order
  std::string domainField;
  int oRder;

  // Object to mark boundary entities for the assembling of domain elements
  boost::shared_ptr<std::vector<unsigned char>> boundaryMarker;

  // Boundary entities marked for fieldsplit (block) solver - optional
  Range boundaryEntitiesForFieldsplit;
  int atom_test = 0;
  enum NORMS { NORM = 0, LAST_NORM };
};

Poisson2DNonhomogeneous::Poisson2DNonhomogeneous(MoFEM::Interface &m_field)
    : domainField("U"), mField(m_field) {}
//! [Read mesh]
MoFEMErrorCode Poisson2DNonhomogeneous::readMesh() {
  MoFEMFunctionBegin;

  CHKERR mField.getInterface(simpleInterface);
  CHKERR simpleInterface->getOptions();
  CHKERR simpleInterface->loadFile();

  MoFEMFunctionReturn(0);
}
//! [Read mesh]

//! [Setup problem]
MoFEMErrorCode Poisson2DNonhomogeneous::setupProblem() {
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
  int oRder = 3;
  CHKERR PetscOptionsGetInt(PETSC_NULL, "", "-atom_test", &atom_test,
                            PETSC_NULL);
  CHKERR PetscOptionsGetInt(PETSC_NULL, "", "-order", &oRder, PETSC_NULL);
  CHKERR simpleInterface->setFieldOrder(domainField, oRder);

  CHKERR simpleInterface->setUp();

  MoFEMFunctionReturn(0);
}
//! [Setup problem]

//! [Boundary condition]
MoFEMErrorCode Poisson2DNonhomogeneous::boundaryCondition() {
  MoFEMFunctionBegin;

  auto bc_mng = mField.getInterface<BcManager>();
  CHKERR bc_mng->pushMarkDOFsOnEntities(simpleInterface->getProblemName(),
                                        "BOUNDARY_CONDITION", domainField, 0, 1,
                                        true);
  // merge markers from all blocksets "BOUNDARY_CONDITION"
  boundaryMarker = bc_mng->getMergedBlocksMarker({"BOUNDARY_CONDITION"});
  // get entities on blocksets "BOUNDARY_CONDITION"
  boundaryEntitiesForFieldsplit =
      bc_mng->getMergedBlocksRange({"BOUNDARY_CONDITION"});

  MoFEMFunctionReturn(0);
}
//! [Boundary condition]

//! [Assemble system]
MoFEMErrorCode Poisson2DNonhomogeneous::assembleSystem() {
  MoFEMFunctionBegin;

  auto pipeline_mng = mField.getInterface<PipelineManager>();

  { // Push operators to the Pipeline that is responsible for calculating LHS of
    // domain elements
    CHKERR AddHOOps<SPACE_DIM, SPACE_DIM, SPACE_DIM>::add(
        pipeline_mng->getOpDomainLhsPipeline(), {H1});
    pipeline_mng->getOpDomainLhsPipeline().push_back(
        new OpSetBc(domainField, true, boundaryMarker));
    pipeline_mng->getOpDomainLhsPipeline().push_back(
        new OpDomainLhs(domainField, domainField));
    pipeline_mng->getOpDomainLhsPipeline().push_back(
        new OpUnSetBc(domainField));
  }

  { // Push operators to the Pipeline that is responsible for calculating RHS of
    // domain elements
    pipeline_mng->getOpDomainRhsPipeline().push_back(
        new OpSetBc(domainField, true, boundaryMarker));
    pipeline_mng->getOpDomainRhsPipeline().push_back(
        new OpDomainRhs(domainField, sourceTermFunction));
    pipeline_mng->getOpDomainRhsPipeline().push_back(
        new OpUnSetBc(domainField));
  }

  { // Push operators to the Pipeline that is responsible for calculating LHS of
    // boundary elements
    pipeline_mng->getOpBoundaryLhsPipeline().push_back(
        new OpSetBc(domainField, false, boundaryMarker));
    pipeline_mng->getOpBoundaryLhsPipeline().push_back(
        new OpBoundaryLhs(domainField, domainField));
    pipeline_mng->getOpBoundaryLhsPipeline().push_back(
        new OpUnSetBc(domainField));
  }

  { // Push operators to the Pipeline that is responsible for calculating RHS of
    // boundary elements
    pipeline_mng->getOpBoundaryRhsPipeline().push_back(
        new OpSetBc(domainField, false, boundaryMarker));
    pipeline_mng->getOpBoundaryRhsPipeline().push_back(
        new OpBoundaryRhs(domainField, boundaryFunction));
    pipeline_mng->getOpBoundaryRhsPipeline().push_back(
        new OpUnSetBc(domainField));
  }

  MoFEMFunctionReturn(0);
}
//! [Assemble system]

//! [Set integration rules]
MoFEMErrorCode Poisson2DNonhomogeneous::setIntegrationRules() {
  MoFEMFunctionBegin;

  auto pipeline_mng = mField.getInterface<PipelineManager>();

  auto domain_rule_lhs = [](int, int, int p) -> int { return 2 * (p - 1); };
  auto domain_rule_rhs = [](int, int, int p) -> int { return 2 * (p - 1); };
  CHKERR pipeline_mng->setDomainLhsIntegrationRule(domain_rule_lhs);
  CHKERR pipeline_mng->setDomainRhsIntegrationRule(domain_rule_rhs);

  auto boundary_rule_lhs = [](int, int, int p) -> int { return 2 * p; };
  auto boundary_rule_rhs = [](int, int, int p) -> int { return 2 * p; };
  CHKERR pipeline_mng->setBoundaryLhsIntegrationRule(boundary_rule_lhs);
  CHKERR pipeline_mng->setBoundaryLhsIntegrationRule(boundary_rule_rhs);

  MoFEMFunctionReturn(0);
}
//! [Set integration rules]

//! [Solve system]
MoFEMErrorCode Poisson2DNonhomogeneous::solveSystem() {
  MoFEMFunctionBegin;

  auto pipeline_mng = mField.getInterface<PipelineManager>();

  auto ksp_solver = pipeline_mng->createKSP();
  CHKERR KSPSetFromOptions(ksp_solver);

  // Create RHS and solution vectors
  auto dm = simpleInterface->getDM();
  auto F = createDMVector(dm);
  auto D = vectorDuplicate(F);

  // Setup fieldsplit (block) solver - optional: yes/no
  if (1) {
    PC pc;
    CHKERR KSPGetPC(ksp_solver, &pc);
    PetscBool is_pcfs = PETSC_FALSE;
    PetscObjectTypeCompare((PetscObject)pc, PCFIELDSPLIT, &is_pcfs);

    // Set up FIELDSPLIT, only when user set -pc_type fieldsplit
    // Identify the index for boundary entities, remaining will be for domain
    // Then split the fields for boundary and domain for solving
    if (is_pcfs == PETSC_TRUE) {
      const MoFEM::Problem *problem_ptr;
      CHKERR DMMoFEMGetProblemPtr(dm, &problem_ptr);
      IS is_boundary;
      CHKERR mField.getInterface<ISManager>()->isCreateProblemFieldAndRank(
          problem_ptr->getName(), ROW, domainField, 0, 1, &is_boundary,
          &boundaryEntitiesForFieldsplit);
      // CHKERR ISView(is_boundary, PETSC_VIEWER_STDOUT_SELF);
      CHKERR PCFieldSplitSetIS(pc, NULL, is_boundary);
      CHKERR ISDestroy(&is_boundary);
    }
  }

  CHKERR KSPSetUp(ksp_solver);

  // Solve the system
  CHKERR KSPSolve(ksp_solver, F, D);

  // Scatter result data on the mesh
  CHKERR VecGhostUpdateBegin(D, INSERT_VALUES, SCATTER_FORWARD);
  CHKERR VecGhostUpdateEnd(D, INSERT_VALUES, SCATTER_FORWARD);
  CHKERR DMoFEMMeshToLocalVector(dm, D, INSERT_VALUES, SCATTER_REVERSE);

  MoFEMFunctionReturn(0);
}
//! [Solve system]

//! [Output results]
MoFEMErrorCode Poisson2DNonhomogeneous::outputResults() {
  MoFEMFunctionBegin;

  auto pipeline_mng = mField.getInterface<PipelineManager>();
  pipeline_mng->getDomainLhsFE().reset();
  pipeline_mng->getBoundaryLhsFE().reset();
  pipeline_mng->getDomainRhsFE().reset();
  pipeline_mng->getBoundaryRhsFE().reset();

  auto d_ptr = boost::make_shared<VectorDouble>();
  auto l_ptr = boost::make_shared<VectorDouble>();

  using OpPPMap = OpPostProcMapInMoab<SPACE_DIM, SPACE_DIM>;

  auto post_proc_fe = boost::make_shared<PostProcFaceEle>(mField);

  post_proc_fe->getOpPtrVector().push_back(
      new OpCalculateScalarFieldValues(domainField, d_ptr));
  post_proc_fe->getOpPtrVector().push_back(new OpPPMap(
      post_proc_fe->getPostProcMesh(), post_proc_fe->getMapGaussPts(),
      {{domainField, d_ptr}}, {}, {}, {}));
  pipeline_mng->getDomainRhsFE() = post_proc_fe;

  CHKERR pipeline_mng->loopFiniteElements();
  CHKERR post_proc_fe->writeFile("out_result.h5m");

  MoFEMFunctionReturn(0);
}
//! [Output results]

//! [Check]
MoFEMErrorCode Poisson2DNonhomogeneous::checkResults() {
  MoFEMFunctionBegin;

  auto check_result_fe_ptr = boost::make_shared<FaceEle>(mField);
  auto petscVec =
      createVectorMPI(mField.get_comm(),
                      (mField.get_comm_rank() == 0) ? LAST_NORM : 0, LAST_NORM);

  CHK_THROW_MESSAGE((AddHOOps<SPACE_DIM, SPACE_DIM, SPACE_DIM>::add(
                        check_result_fe_ptr->getOpPtrVector(), {H1})),
                    "Apply transform");

  check_result_fe_ptr->getRuleHook = [](int, int, int p) { return 2 * p; };
  auto analyticalFunction = [&](double x, double y, double z) {
    return cos(x * M_PI) * cos(y * M_PI);
  };

  auto u_ptr = boost::make_shared<VectorDouble>();

  check_result_fe_ptr->getOpPtrVector().push_back(
      new OpCalculateScalarFieldValues(domainField, u_ptr));
  auto mValFuncPtr = boost::make_shared<VectorDouble>();
  check_result_fe_ptr->getOpPtrVector().push_back(
      new OpGetTensor0fromFunc(mValFuncPtr, analyticalFunction));
  check_result_fe_ptr->getOpPtrVector().push_back(
      new OpCalcNormL2Tensor0(u_ptr, petscVec, NORM, mValFuncPtr));
  CHKERR VecZeroEntries(petscVec);
  CHKERR DMoFEMLoopFiniteElements(simpleInterface->getDM(),
                                  simpleInterface->getDomainFEName(),
                                  check_result_fe_ptr);
  CHKERR VecAssemblyBegin(petscVec);
  CHKERR VecAssemblyEnd(petscVec);
  MOFEM_LOG_CHANNEL("SELF"); // Clear channel from old tags
  if (mField.get_comm_rank() == 0) {
    const double *norms;
    CHKERR VecGetArrayRead(petscVec, &norms);
    MOFEM_TAG_AND_LOG("SELF", Sev::inform, "Errors")
        << "NORM: " << std::sqrt(norms[NORM]);
    CHKERR VecRestoreArrayRead(petscVec, &norms);
  }
  if (atom_test && !mField.get_comm_rank()) {
    const double *t_ptr;
    CHKERR VecGetArrayRead(petscVec, &t_ptr);
    double ref_norm = 1.4e-04;
    double cal_norm;
    switch (atom_test) {
    case 1: // 2D
      cal_norm = sqrt(t_ptr[0]);
      break;
    default:
      SETERRQ1(PETSC_COMM_SELF, MOFEM_INVALID_DATA,
               "atom test %d does not exist", atom_test);
    }
    if (cal_norm > ref_norm) {
      SETERRQ3(PETSC_COMM_SELF, MOFEM_ATOM_TEST_INVALID,
               "atom test %d failed! Calculated Norm %3.16e is greater than "
               "reference Norm %3.16e",
               atom_test, cal_norm, ref_norm);
    }
    CHKERR VecRestoreArrayRead(petscVec, &t_ptr);
  }
  MoFEMFunctionReturn(0);
}
//! [Check]

//! [Run program]
MoFEMErrorCode Poisson2DNonhomogeneous::runProgram() {
  MoFEMFunctionBegin;

  CHKERR readMesh();
  CHKERR setupProblem();
  CHKERR boundaryCondition();
  CHKERR assembleSystem();
  CHKERR setIntegrationRules();
  CHKERR solveSystem();
  CHKERR outputResults();
  CHKERR checkResults();

  MoFEMFunctionReturn(0);
}
//! [Run program]

int main(int argc, char *argv[]) {

  // Initialisation of MoFEM/PETSc and MOAB data structures
  const char param_file[] = "param_file.petsc";
  MoFEM::Core::Initialize(&argc, &argv, param_file, help);

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
    Poisson2DNonhomogeneous poisson_problem(m_field);
    CHKERR poisson_problem.runProgram();
  }
  CATCH_ERRORS;

  // Finish work: cleaning memory, getting statistics, etc.
  MoFEM::Core::Finalize();

  return 0;
}
