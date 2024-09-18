/**
 * \file mixed_nonlinear_poisson.cpp
 * \example mixed_nonlinear_poisson.cpp
 *
 * MixedNonlinearPoisson intended to show how to solve mixed formulation of the
 * Dirichlet problem for the Poisson equation using error indicators and
 * p-adaptivity
 *
 */
#ifndef EXECUTABLE_DIMENSION
  #define EXECUTABLE_DIMENSION 2
#endif

#include <stdlib.h>
#include <MoFEM.hpp>
#include <mixed_nonlinear_poisson.hpp>

using namespace MoFEM;
using namespace MixedNonlinearPoissonOps;

static char help[] = "...\n\n";

inline double sqr(double x) { return x * x; }

inline double cube(const double x) { return x * x * x; }

struct MixedNonlinearPoisson {

  MixedNonlinearPoisson(MoFEM::Interface &m_field) : mField(m_field) {}
  MoFEMErrorCode runProblem();

private:
  MoFEM::Interface &mField;
  Simple *simpleInterface;

  int initOrder;

  //! [Source function]
  static double sourceFunction(const double x, const double y, const double z) {
    return 2 * M_PI * M_PI *
           (cos(M_PI * x) * cos(M_PI * y) +
            cube(cos(M_PI * x)) * cube(cos(M_PI * y)) -
            cos(M_PI * x) * cos(M_PI * y) *
                (sqr(sin(M_PI * x)) * sqr(cos(M_PI * y)) +
                 sqr(sin(M_PI * y)) * sqr(cos(M_PI * x))));
  }
  //! [Source function]

  //! [Boundary function]
  // Function to calculate the Boundary term
  static double boundaryFunction(const double x, const double y,
                                 const double z) {
    return -cos(M_PI * x) *
           cos(M_PI * y); // here should put the negative of the proper formula
  }
  //! [Boundary function]

  MoFEMErrorCode readMesh();
  MoFEMErrorCode setupProblem();
  MoFEMErrorCode setIntegrationRules();
  MoFEMErrorCode assembleSystem();
  MoFEMErrorCode solveSystem();
  MoFEMErrorCode outputResults(int iter_num = 0);
};

//! [Run programme]
MoFEMErrorCode MixedNonlinearPoisson::runProblem() {
  MoFEMFunctionBegin;
  CHKERR readMesh();
  CHKERR setupProblem();
  CHKERR setIntegrationRules();
  CHKERR assembleSystem();
  CHKERR solveSystem();
  CHKERR outputResults();
  MoFEMFunctionReturn(0);
}
//! [Run programme]

//! [Read mesh]
MoFEMErrorCode MixedNonlinearPoisson::readMesh() {
  MoFEMFunctionBegin;
  CHKERR mField.getInterface(simpleInterface);
  CHKERR simpleInterface->getOptions();
  CHKERR simpleInterface->loadFile();
  MoFEMFunctionReturn(0);
}
//! [Read mesh]

//! [Set up problem]
MoFEMErrorCode MixedNonlinearPoisson::setupProblem() {
  MoFEMFunctionBegin;

  CHKERR simpleInterface->addDomainField("U", L2, AINSWORTH_LEGENDRE_BASE, 1);
  CHKERR simpleInterface->addBoundaryField("U", L2, AINSWORTH_LEGENDRE_BASE, 1);

  int nb_quads = 0;
  CHKERR mField.get_moab().get_number_entities_by_type(0, MBQUAD, nb_quads);
  auto base = AINSWORTH_LEGENDRE_BASE;
  if (nb_quads) {
    // AINSWORTH_LEGENDRE_BASE is not implemented for HDIV/HCURL space on quads
    base = DEMKOWICZ_JACOBI_BASE;
  }

  // Note that in 2D case HDIV and HCURL spaces are isomorphic, and therefore
  // only base for HCURL has been implemented in 2D. Base vectors for HDIV space
  // are be obtained after rotation of HCURL base vectors by a right angle
  CHKERR simpleInterface->addDomainField("Q", HCURL, base, 1);
  CHKERR simpleInterface->addBoundaryField("Q", HCURL, base, 1);

  initOrder = 2;
  CHKERR PetscOptionsGetInt(PETSC_NULL, "", "-order", &initOrder, PETSC_NULL);
  CHKERR simpleInterface->setFieldOrder("U", initOrder);
  CHKERR simpleInterface->setFieldOrder("Q", initOrder + 1);
  CHKERR simpleInterface->setUp();

  MoFEMFunctionReturn(0);
}
//! [Set up problem]

//! [Set integration rule]
MoFEMErrorCode MixedNonlinearPoisson::setIntegrationRules() {
  MoFEMFunctionBegin;

  auto rule = [](int, int, int p) -> int { return 2 * p; };

  PipelineManager *pipeline_mng = mField.getInterface<PipelineManager>();
  CHKERR pipeline_mng->setDomainLhsIntegrationRule(rule);
  CHKERR pipeline_mng->setDomainRhsIntegrationRule(rule);
  // CHKERR pipeline_mng->setBoundaryLhsIntegrationRule(rule);
  CHKERR pipeline_mng->setBoundaryRhsIntegrationRule(rule);

  MoFEMFunctionReturn(0);
}
//! [Set integration rule]

//! [Assemble system]
MoFEMErrorCode MixedNonlinearPoisson::assembleSystem() {
  MoFEMFunctionBegin;
  PipelineManager *pipeline_mng = mField.getInterface<PipelineManager>();
  pipeline_mng->getDomainLhsFE().reset();
  pipeline_mng->getDomainRhsFE().reset();
  pipeline_mng->getOpDomainRhsPipeline().clear();
  pipeline_mng->getOpDomainLhsPipeline().clear();

  CHKERR AddHOOps<SPACE_DIM, SPACE_DIM, SPACE_DIM>::add(
      pipeline_mng->getOpDomainLhsPipeline(), {HDIV});
  CHKERR AddHOOps<SPACE_DIM, SPACE_DIM, SPACE_DIM>::add(
      pipeline_mng->getOpDomainRhsPipeline(), {HDIV});

  auto add_domain_lhs_ops = [&](auto &pipeline) {
    auto data_u_at_gauss_pts = boost::make_shared<VectorDouble>();
    auto data_q_at_gauss_pts = boost::make_shared<MatrixDouble>();

    pipeline.push_back(
        new OpCalculateScalarFieldValues("U", data_u_at_gauss_pts));
    pipeline.push_back(
        new OpCalculateHVecVectorField<3>("Q", data_q_at_gauss_pts));

    auto unity = []() { return 1; };
    pipeline.push_back(new OpHdivU("Q", "U", unity, true));
    pipeline.push_back(
        new OpDomainQULhs("Q", "U", data_u_at_gauss_pts, data_q_at_gauss_pts));
    pipeline.push_back(new OpDomainQQLhs("Q", "Q", data_u_at_gauss_pts));
  };

  auto add_domain_rhs_ops = [&](auto &pipeline) {
    auto data_u_at_gauss_pts = boost::make_shared<VectorDouble>();
    auto data_divq_at_gauss_pts = boost::make_shared<VectorDouble>();
    auto data_q_at_gauss_pts = boost::make_shared<MatrixDouble>();
    auto source = [&](const double x, const double y, const double z) {
      return sourceFunction(x, y, z);
    };
    auto unity = [](const double, const double, const double) constexpr {
      return 1.;
    };

    pipeline.push_back(
        new OpCalculateScalarFieldValues("U", data_u_at_gauss_pts));
    pipeline.push_back(new OpCalculateHdivVectorDivergence<3, SPACE_DIM>(
        "Q", data_divq_at_gauss_pts));
    pipeline.push_back(
        new OpCalculateHVecVectorField<3>("Q", data_q_at_gauss_pts));

    pipeline.push_back(new OpBaseTimesScalarRhs("U", data_divq_at_gauss_pts, unity));
    pipeline.push_back(new OpDomainSource("U", source));
    pipeline.push_back(
        new OpHdivUHdivRhs("Q", data_u_at_gauss_pts, data_q_at_gauss_pts));
    pipeline.push_back(new OpHDivTimesScalarRhs("Q", data_u_at_gauss_pts, unity));
  };

  auto add_boundary_rhs_ops = [&](auto &pipeline) {
    // pipeline.push_back(new OpSetBc("Q", false, boundaryMarker));
    pipeline.push_back(new OpBoundaryRhsSource("Q", boundaryFunction));
    // pipeline.push_back(new OpUnSetBc("Q"));
  };

  add_domain_lhs_ops(pipeline_mng->getOpDomainLhsPipeline());
  add_domain_rhs_ops(pipeline_mng->getOpDomainRhsPipeline());
  add_boundary_rhs_ops(pipeline_mng->getOpBoundaryRhsPipeline());

  MoFEMFunctionReturn(0);
}
//! [Assemble system]

//! [Solve]
MoFEMErrorCode MixedNonlinearPoisson::solveSystem() {
  MoFEMFunctionBegin;
  auto dm = simpleInterface->getDM();
  SmartPetscObj<Vec> global_rhs, global_solution;
  CHKERR DMCreateGlobalVector_MoFEM(dm, global_rhs);
  global_solution = vectorDuplicate(global_rhs);

  // Create nonlinear solver (SNES)
  auto pipeline_mng = mField.getInterface<PipelineManager>();
  auto solver = pipeline_mng->createSNES();
  CHKERR SNESSetFromOptions(solver);
  CHKERR SNESSetUp(solver);
  // Solve the system
  CHKERR SNESSolve(solver, global_rhs, global_solution);
  CHKERR VecGhostUpdateBegin(global_solution, INSERT_VALUES, SCATTER_FORWARD);
  CHKERR VecGhostUpdateEnd(global_solution, INSERT_VALUES, SCATTER_FORWARD);
  // Scatter result data on the mesh
  CHKERR DMoFEMMeshToGlobalVector(dm, global_solution, INSERT_VALUES,
                                  SCATTER_REVERSE);
  MoFEMFunctionReturn(0);
}
//! [Solve]

//! [Output results]
MoFEMErrorCode MixedNonlinearPoisson::outputResults(int iter_num) {
  MoFEMFunctionBegin;
  PipelineManager *pipeline_mng = mField.getInterface<PipelineManager>();
  pipeline_mng->getDomainLhsFE().reset();

  auto post_proc_fe = boost::make_shared<PostProcEle>(mField);

  CHKERR AddHOOps<SPACE_DIM, SPACE_DIM, SPACE_DIM>::add(
      post_proc_fe->getOpPtrVector(), {HDIV});

  auto u_ptr = boost::make_shared<VectorDouble>();
  auto flux_ptr = boost::make_shared<MatrixDouble>();
  post_proc_fe->getOpPtrVector().push_back(
      new OpCalculateScalarFieldValues("U", u_ptr));
  post_proc_fe->getOpPtrVector().push_back(
      new OpCalculateHVecVectorField<3,SPACE_DIM>("Q", flux_ptr));

  using OpPPMap = OpPostProcMapInMoab<SPACE_DIM, SPACE_DIM>;

  post_proc_fe->getOpPtrVector().push_back(

      new OpPPMap(post_proc_fe->getPostProcMesh(),
                  post_proc_fe->getMapGaussPts(),

                  OpPPMap::DataMapVec{{"U", u_ptr}},

                  OpPPMap::DataMapMat{{"Q", flux_ptr}},

                  OpPPMap::DataMapMat{},

                  OpPPMap::DataMapMat{}

                  )

  );

  // pipeline_mng->getDomainRhsFE() = post_proc_fe;
  // CHKERR pipeline_mng->loopFiniteElements();

  // std::ostringstream strm;
  // strm << "out_" << iter_num << ".h5m";
  // CHKERR post_proc_fe->writeFile(strm.str().c_str());
  auto *simple = mField.getInterface<Simple>();
  auto dm = simple->getDM();
  CHKERR DMoFEMLoopFiniteElements(dm, simple->getDomainFEName(), post_proc_fe);

  CHKERR post_proc_fe->writeFile("out_result.h5m");
  MoFEMFunctionReturn(0);
}
//! [Output results]

int main(int argc, char *argv[]) {
  // Initialisation of MoFEM/PETSc and MOAB data structures
  const char param_file[] = "param_file.petsc";
  MoFEM::Core::Initialize(&argc, &argv, param_file, help);

  // Add logging channel for example problem
  auto core_log = logging::core::get();
  core_log->add_sink(
      LogManager::createSink(LogManager::getStrmWorld(), "EXAMPLE"));
  LogManager::setLog("EXAMPLE");
  MOFEM_LOG_TAG("EXAMPLE", "MixedNonlinearPoisson");

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
    MoFEM::Interface &m_field = core; ///< finite element database interface
    //! [Create MoFEM]

    //! [MixedNonlinearPoisson]
    MixedNonlinearPoisson ex(m_field);
    CHKERR ex.runProblem();
    //! [MixedNonlinearPoisson]
  }
  CATCH_ERRORS;

  CHKERR MoFEM::Core::Finalize();
}