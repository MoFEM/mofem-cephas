/**
 * \file poisson_2d_homogeneous.cpp
 * \example poisson_2d_homogeneous.cpp
 *
 * Solution of poisson equation. Direct implementation of User Data Operators
 * for teaching proposes.
 *
 * \note In practical application we suggest use form integrators to generalise
 * and simplify code. However, here we like to expose user to ways how to
 * implement data operator from scratch.
 */

constexpr auto field_name = "U";

constexpr int SPACE_DIM =
    EXECUTABLE_DIMENSION; //< Space dimension of problem, mesh

#include <poisson_2d_homogeneous.hpp>

using namespace MoFEM;
using namespace Poisson2DHomogeneousOperators;

using PostProcFaceEle = PostProcBrokenMeshInMoab<DomainEle>;

static char help[] = "...\n\n";

struct Poisson2DHomogeneous {
public:
  Poisson2DHomogeneous(MoFEM::Interface &m_field);

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

  // MoFEM interfaces
  Simple *simpleInterface;
  // Field name
  MoFEM::Interface &mField;
  // Approximation order
  int oRder;
  // Function to calculate the Source term
  static double sourceTermFunction(const double x, const double y,
                                   const double z) {
    return 2. * M_PI * M_PI * sin(M_PI * x) * sin(M_PI * y);
  }
  // PETSc vector for storing norms
  SmartPetscObj<Vec> petscVec;
  int atom_test = 0;
  enum NORMS { NORM = 0, LAST_NORM };
};

Poisson2DHomogeneous::Poisson2DHomogeneous(MoFEM::Interface &m_field)
    : mField(m_field) {}

//! [Read mesh]
MoFEMErrorCode Poisson2DHomogeneous::readMesh() {
  MoFEMFunctionBegin;

  CHKERR mField.getInterface(simpleInterface);
  CHKERR simpleInterface->getOptions();
  CHKERR simpleInterface->loadFile();

  MoFEMFunctionReturn(0);
}
//! [Read mesh]

//! [Setup problem]
MoFEMErrorCode Poisson2DHomogeneous::setupProblem() {
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
  CHKERR simpleInterface->addDomainField(field_name, H1, base, 1);

  int oRder = 3;
  CHKERR PetscOptionsGetInt(PETSC_NULL, "", "-order", &oRder, PETSC_NULL);
  CHKERR PetscOptionsGetInt(PETSC_NULL, "", "-atom_test", &atom_test,
                            PETSC_NULL);
  CHKERR simpleInterface->setFieldOrder(field_name, oRder);

  CHKERR simpleInterface->setUp();

  MoFEMFunctionReturn(0);
}
//! [Setup problem]

//! [Boundary condition]
MoFEMErrorCode Poisson2DHomogeneous::boundaryCondition() {
  MoFEMFunctionBegin;

  auto bc_mng = mField.getInterface<BcManager>();

  // Remove BCs from blockset name "BOUNDARY_CONDITION" or SETU, note that you
  // can use regular expression to put list of blocksets;
  CHKERR bc_mng->removeBlockDOFsOnEntities<BcScalarMeshsetType<BLOCKSET>>(
      simpleInterface->getProblemName(), "(BOUNDARY_CONDITION|SETU)",
      std::string(field_name), true);

  MoFEMFunctionReturn(0);
}
//! [Boundary condition]

//! [Assemble system]
MoFEMErrorCode Poisson2DHomogeneous::assembleSystem() {
  MoFEMFunctionBegin;

  auto pipeline_mng = mField.getInterface<PipelineManager>();

  { // Push operators to the Pipeline that is responsible for calculating LHS
    CHKERR AddHOOps<SPACE_DIM, SPACE_DIM, SPACE_DIM>::add(
        pipeline_mng->getOpDomainLhsPipeline(), {H1});
    pipeline_mng->getOpDomainLhsPipeline().push_back(
        new OpDomainLhsMatrixK(field_name, field_name));
  }

  { // Push operators to the Pipeline that is responsible for calculating RHS

    auto set_values_to_bc_dofs = [&](auto &fe) {
      auto get_bc_hook = [&]() {
        EssentialPreProc<TemperatureCubitBcData> hook(mField, fe, {});
        return hook;
      };
      fe->preProcessHook = get_bc_hook();
    };

    // you can skip that if boundary condition is prescribing zero
    auto calculate_residual_from_set_values_on_bc = [&](auto &pipeline) {
      using DomainEle =
          PipelineManager::ElementsAndOpsByDim<SPACE_DIM>::DomainEle;
      using DomainEleOp = DomainEle::UserDataOperator;
      using OpInternal = FormsIntegrators<DomainEleOp>::Assembly<
          PETSC>::LinearForm<GAUSS>::OpGradTimesTensor<1, 1, SPACE_DIM>;

      auto grad_u_vals_ptr = boost::make_shared<MatrixDouble>();
      pipeline_mng->getOpDomainRhsPipeline().push_back(
          new OpCalculateScalarFieldGradient<SPACE_DIM>(field_name,
                                                        grad_u_vals_ptr));
      pipeline_mng->getOpDomainRhsPipeline().push_back(
          new OpInternal(field_name, grad_u_vals_ptr,
                         [](double, double, double) constexpr { return -1; }));
    };

    CHKERR AddHOOps<SPACE_DIM, SPACE_DIM, SPACE_DIM>::add(
        pipeline_mng->getOpDomainRhsPipeline(), {H1});
    set_values_to_bc_dofs(pipeline_mng->getDomainRhsFE());
    calculate_residual_from_set_values_on_bc(
        pipeline_mng->getOpDomainRhsPipeline());
    pipeline_mng->getOpDomainRhsPipeline().push_back(
        new OpDomainRhsVectorF(field_name, sourceTermFunction));
  }

  MoFEMFunctionReturn(0);
}
//! [Assemble system]

//! [Set integration rules]
MoFEMErrorCode Poisson2DHomogeneous::setIntegrationRules() {
  MoFEMFunctionBegin;

  auto rule_lhs = [](int, int, int p) -> int { return 2 * (p - 1); };
  auto rule_rhs = [](int, int, int p) -> int { return p; };

  auto pipeline_mng = mField.getInterface<PipelineManager>();
  CHKERR pipeline_mng->setDomainLhsIntegrationRule(rule_lhs);
  CHKERR pipeline_mng->setDomainRhsIntegrationRule(rule_rhs);

  MoFEMFunctionReturn(0);
}
//! [Set integration rules]

//! [Solve system]
MoFEMErrorCode Poisson2DHomogeneous::solveSystem() {
  MoFEMFunctionBegin;

  auto pipeline_mng = mField.getInterface<PipelineManager>();

  auto ksp_solver = pipeline_mng->createKSP();
  CHKERR KSPSetFromOptions(ksp_solver);
  CHKERR KSPSetUp(ksp_solver);

  // Create RHS and solution vectors
  auto dm = simpleInterface->getDM();
  auto F = createDMVector(dm);
  auto D = vectorDuplicate(F);

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
MoFEMErrorCode Poisson2DHomogeneous::outputResults() {
  MoFEMFunctionBegin;

  auto pipeline_mng = mField.getInterface<PipelineManager>();
  pipeline_mng->getDomainLhsFE().reset();

  auto post_proc_fe = boost::make_shared<PostProcFaceEle>(mField);
  CHKERR AddHOOps<SPACE_DIM, SPACE_DIM, SPACE_DIM>::add(
      post_proc_fe->getOpPtrVector(), {H1});

  auto u_ptr = boost::make_shared<VectorDouble>();
  auto grad_u_ptr = boost::make_shared<MatrixDouble>();
  post_proc_fe->getOpPtrVector().push_back(
      new OpCalculateScalarFieldValues(field_name, u_ptr));

  post_proc_fe->getOpPtrVector().push_back(
      new OpCalculateScalarFieldGradient<SPACE_DIM>(field_name, grad_u_ptr));

  using OpPPMap = OpPostProcMapInMoab<SPACE_DIM, SPACE_DIM>;

  post_proc_fe->getOpPtrVector().push_back(

      new OpPPMap(post_proc_fe->getPostProcMesh(),
                  post_proc_fe->getMapGaussPts(),

                  OpPPMap::DataMapVec{{"U", u_ptr}},

                  OpPPMap::DataMapMat{{"GRAD_U", grad_u_ptr}},

                  OpPPMap::DataMapMat{},

                  OpPPMap::DataMapMat{}

                  )

  );

  pipeline_mng->getDomainRhsFE() = post_proc_fe;
  CHKERR pipeline_mng->loopFiniteElements();
  CHKERR post_proc_fe->writeFile("out_result.h5m");

  MoFEMFunctionReturn(0);
}
//! [Output results]

//! [Check]
MoFEMErrorCode Poisson2DHomogeneous::checkResults() {
  MoFEMFunctionBegin;

  auto check_result_fe_ptr = boost::make_shared<DomainEle>(mField);
  auto petscVec =
      createVectorMPI(mField.get_comm(),
                      (mField.get_comm_rank() == 0) ? LAST_NORM : 0, LAST_NORM);

  CHK_THROW_MESSAGE((AddHOOps<SPACE_DIM, SPACE_DIM, SPACE_DIM>::add(
                        check_result_fe_ptr->getOpPtrVector(), {H1})),
                    "Apply transform");

  check_result_fe_ptr->getRuleHook = [](int, int, int p) { return p; };
  auto analyticalFunction = [&](double x, double y, double z) {
    return sin(M_PI * x) * sin(M_PI * y);
  };

  auto u_ptr = boost::make_shared<VectorDouble>();

  check_result_fe_ptr->getOpPtrVector().push_back(
      new OpCalculateScalarFieldValues(field_name, u_ptr));
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
  // print norm in general log
  if (mField.get_comm_rank() == 0) {
    const double *norms;
    CHKERR VecGetArrayRead(petscVec, &norms);
    MOFEM_TAG_AND_LOG("SELF", Sev::inform, "Errors")
        << "NORM: " << std::sqrt(norms[NORM]);
    CHKERR VecRestoreArrayRead(petscVec, &norms);
  }
  // compare norm for ctest
  if (atom_test && !mField.get_comm_rank()) {
    const double *t_ptr;
    CHKERR VecGetArrayRead(petscVec, &t_ptr);
    double ref_norm = 2.2e-04;
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
MoFEMErrorCode Poisson2DHomogeneous::runProgram() {
  MoFEMFunctionBegin;

  CHKERR readMesh();
  CHKERR setupProblem();
  CHKERR boundaryCondition();
  CHKERR setIntegrationRules();
  CHKERR assembleSystem();
  CHKERR solveSystem();
  CHKERR outputResults();
  CHKERR checkResults();

  MoFEMFunctionReturn(0);
}
//! [Run program]

//! [Main]
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
    Poisson2DHomogeneous poisson_problem(m_field);
    CHKERR poisson_problem.runProgram();
  }
  CATCH_ERRORS;

  // Finish work: cleaning memory, getting statistics, etc.
  MoFEM::Core::Finalize();

  return 0;
}
//! [Main]