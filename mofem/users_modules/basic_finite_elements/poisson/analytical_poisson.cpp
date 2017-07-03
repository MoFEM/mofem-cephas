/**
 * \file analytical_poisson.cpp
 * \ingroup mofem_simple_interface
 * \example analytical_poisson.cpp
 *
 * For more information and detailed explain of this
 * example see \ref poisson_tut1
 *
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

#include <BasicFiniteElements.hpp>
#include <PoissonOperators.hpp>
#include <AuxPoissonFunctions.hpp>

static char help[] = "...\n\n";

/**
 * \brief Function
 *
 * This is prescribed exact function. If this function is given by polynomial
 * order of "p" and order of approximation is "p" or higher, solution of
 * finite element method is exact (with machine precision).
 *
 *  \f[
 *  u = 1+x^2+y^2+z^3
 *  \f]
 *
 */
struct ExactFunction {
  double operator()(const double x,const double y,const double z) const {
    return 1+x*x+y*y+z*z*z;
  }
};

/**
 * \brief Exact gradient
 */
struct ExactFunctionGrad {
  FTensor::Tensor1<double,3> operator()(const double x,const double y,const double z) const {
    FTensor::Tensor1<double,3> grad;
    grad(0) = 2*x;
    grad(1) = 2*y;
    grad(2) = 3*z*z;
    return grad;
  }
};

/**
 * \brief Laplacian of function.
 *
 * This is Laplacian of \f$u\f$, it is calculated using formula
 * \f[
 * \nabla^2 u(x,y,z) = \nabla \cdot \nabla u
 * \frac{\partial^2 u}{\partial x^2}+
 * \frac{\partial^2 u}{\partial y^2}+
 * \frac{\partial^2 u}{\partial z^2}
 * \f]
 *
 */
struct ExactLaplacianFunction {
  double operator()(const double x,const double y,const double z) const {
    return 4+6*z;
  }
};

int main(int argc, char *argv[]) {

  ErrorCode rval;
  PetscErrorCode ierr;

  // Initialize PETSc
  PetscInitialize(&argc,&argv,(char *)0,help);
  // Create MoAB database
  moab::Core moab_core;                   // create database
  moab::Interface& moab = moab_core;      // create interface to database

  // Get command line options
  int order = 3;  // default approximation order
  PetscBool flg_test = PETSC_FALSE; // true check if error is numerical error
  ierr = PetscOptionsBegin(PETSC_COMM_WORLD,"", "Poisson's problem options","none"); CHKERRQ(ierr);
  // Set approximation order
  ierr = PetscOptionsInt("-order","approximation order","",order,&order,PETSC_NULL); CHKERRQ(ierr);
  // Set testing (used by CTest)
  ierr = PetscOptionsBool("-test","if true is ctest","",flg_test,&flg_test,PETSC_NULL); CHKERRQ(ierr);
  ierr = PetscOptionsEnd(); CHKERRQ(ierr);

  try {

    // Create MoFEM database and link it to MoAB
    MoFEM::Core mofem_core(moab);                      // create database
    MoFEM::Interface& m_field = mofem_core;            // create interface to database
    // Register DM Manager
    ierr = DMRegister_MoFEM("DMMOFEM"); CHKERRQ(ierr); // register MoFEM DM in PETSc

    // Create vector to store approximation global error
    Vec global_error;
    ierr = PoissonExample::AuxFunctions(m_field).createGhostVec(&global_error); CHKERRQ(ierr);

    // First we crate elements, implementation of elements is problem independent,
    // we do not know yet what fields are present in the problem, or
    // even we do not decided yet what approximation base or spaces we
    // are going to use. Implementation of element is free from
    // those constrains and can be used in different context.

    // Elements used by KSP & DM to assemble system of equations
    boost::shared_ptr<ForcesAndSurcesCore> domain_lhs_fe;     ///< Volume element for the matrix
    boost::shared_ptr<ForcesAndSurcesCore> boundary_lhs_fe;   ///< Boundary element for the matrix
    boost::shared_ptr<ForcesAndSurcesCore> domain_rhs_fe;     ///< Volume element to assemble vector
    boost::shared_ptr<ForcesAndSurcesCore> boundary_rhs_fe;   ///< Volume element to assemble vector
    boost::shared_ptr<ForcesAndSurcesCore> domain_error;      ///< Volume element evaluate error
    boost::shared_ptr<ForcesAndSurcesCore> post_proc_volume;  ///< Volume element to Post-process results
    boost::shared_ptr<ForcesAndSurcesCore> null;              ///< Null element do nothing
    {
      // Add problem specific operators the generic finite elements to calculate matrices and vectors.
      ierr = PoissonExample::CreateFiniteElementes(m_field).createFEToAssmbleMatrceAndVector(
        ExactFunction(),ExactLaplacianFunction(),
        domain_lhs_fe,boundary_lhs_fe,domain_rhs_fe,boundary_rhs_fe
      ); CHKERRQ(ierr);
      // Add problem specific operators the generic finite elements to calculate error on elements and global error
      // in H1 norm
      ierr = PoissonExample::CreateFiniteElementes(m_field).createFEToEvaluateError(
        ExactFunction(),ExactFunctionGrad(),global_error,domain_error
      ); CHKERRQ(ierr);
      // Post-process results
      ierr = PoissonExample::CreateFiniteElementes(m_field).creatFEToPostProcessResults(post_proc_volume); CHKERRQ(ierr);
    }

    // Get simple interface is simplified version enabling quick and
    // easy construction of problem.
    Simple *simple_interface;
    // Query interface and get pointer to Simple interface
    ierr = m_field.query_interface(simple_interface); CHKERRQ(ierr);

    // Build problem with simple interface
    {

      // Get options for simple interface from command line
      ierr = simple_interface->getOptions(); CHKERRQ(ierr);
      // Load mesh file to database
      ierr = simple_interface->loadFile(); CHKERRQ(ierr);

      // Add field on domain and boundary. Field is declared by space and base and rank. space
      // can be H1. Hcurl, Hdiv and L2, base can be AINSWORTH_LEGENDRE_BASE, DEMKOWICZ_JACOBI_BASE and more,
      // where rank is number of coefficients for dof.
      //
      // Simple interface assumes that there are four types of field; 1) defined
      // on body domain, 2) fields defined on body boundary, 3) skeleton field defined
      // on finite element skeleton. Finally data field is defined on body domain. Data field
      // is not solved but used for post-process or to keep material parameters, etc. Here
      // we using it to calculate approximation error on elements.

      // Add domain filed "U" in space H1 and Legendre base, Ainsworth recipe is used
      // to construct base functions.
      ierr = simple_interface->addDomainField("U",H1,AINSWORTH_LEGENDRE_BASE,1); CHKERRQ(ierr);
      // Add Lagrange multiplier field on body boundary
      ierr = simple_interface->addBoundaryField("L",H1,AINSWORTH_LEGENDRE_BASE,1); CHKERRQ(ierr);
      // Add error (data) field, we need only L2 norm. Later order is set to 0, so this
      // is piecewise discontinuous constant approx., i.e. 1 DOF for element. You can use
      // more DOFs and collate moments of error to drive anisotropic h/p-adaptivity, however
      // this is beyond this example.
      ierr = simple_interface->addDataField("ERROR",L2,AINSWORTH_LEGENDRE_BASE,1); CHKERRQ(ierr);

      // Set fields order domain and boundary fields.
      ierr = simple_interface->setFieldOrder("U",order); CHKERRQ(ierr); // to approximate function
      ierr = simple_interface->setFieldOrder("L",order); CHKERRQ(ierr); // to Lagrange multiplayers
      ierr = simple_interface->setFieldOrder("ERROR",0); CHKERRQ(ierr); // approximation order for error

      // Setup problem. At that point database is constructed, i.e. fields, finite elements and
      // problem data structures with local and global indexing.
      ierr = simple_interface->setUp(); CHKERRQ(ierr);

    }

    // Get access to PETSC-MoFEM DM manager.
    // or more derails see <http://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/DM/index.html>
    // Form that point internal MoFEM data structures are linked with PETSc interface. If
    // DM functions contains string MoFEM is is MoFEM specific DM interface function,
    // otherwise DM function part of standard PETSc interface.

    DM dm;
    // Get dm
    ierr = simple_interface->getDM(&dm); CHKERRQ(ierr);

    // Set KSP context for DM. At that point only elements are added to DM operators.
    // Calculations of matrices and vectors is executed by KSP solver. This part
    // of the code makes connection between implementation of finite elements and
    // data operators with finite element declarations in DM manager. From that
    // point DM takes responsibility for executing elements, calculations of
    // matrices and vectors, and solution of the problem.
    {
      // Set operators for KSP solver
      ierr = DMMoFEMKSPSetComputeOperators(
        dm,simple_interface->getDomainFEName(),domain_lhs_fe,null,null
      ); CHKERRQ(ierr);
      ierr = DMMoFEMKSPSetComputeOperators(
        dm,simple_interface->getBoundaryFEName(),boundary_lhs_fe,null,null
      ); CHKERRQ(ierr);
      // Set calculation of the right hand side vetor for KSP solver
      ierr = DMMoFEMKSPSetComputeRHS(
        dm,simple_interface->getDomainFEName(),domain_rhs_fe,null,null
      ); CHKERRQ(ierr);
      ierr = DMMoFEMKSPSetComputeRHS(
        dm,simple_interface->getBoundaryFEName(),boundary_rhs_fe,null,null
      ); CHKERRQ(ierr);
    }

    // Solve problem, only PETEc interface here.
    {

      // Create the right hand side vector and vector of unknowns
      Vec F,D;
      ierr = DMCreateGlobalVector(dm,&F); CHKERRQ(ierr);
      // Create unknown vector by creating duplicate copy of F vector. only
      // structure is duplicated no values.
      ierr = VecDuplicate(F,&D); CHKERRQ(ierr);

      // Create solver and link it to DM
      KSP solver;
      ierr = KSPCreate(PETSC_COMM_WORLD,&solver); CHKERRQ(ierr);
      ierr = KSPSetFromOptions(solver); CHKERRQ(ierr);
      ierr = KSPSetDM(solver,dm); CHKERRQ(ierr);
      // Set-up solver, is type of solver and pre-conditioners
      ierr = KSPSetUp(solver); CHKERRQ(ierr);
      // At solution process, KSP solver using DM creates matrices, Calculate
      // values of the left hand side and the right hand side vector. then
      // solves system of equations. Results are stored in vector D.
      ierr = KSPSolve(solver,F,D); CHKERRQ(ierr);

      // Scatter solution on the mesh. Stores unknown vector on field on the mesh.
      ierr = DMoFEMMeshToGlobalVector(dm,D,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);

      // Clean data. Solver and vector are not needed any more.
      ierr = KSPDestroy(&solver); CHKERRQ(ierr);
      ierr = VecDestroy(&D); CHKERRQ(ierr);
      ierr = VecDestroy(&F); CHKERRQ(ierr);
    }

    // Calculate error
    {
      // Loop over all elements in mesh, and run users operators on each element.
      ierr = DMoFEMLoopFiniteElements(dm,simple_interface->getDomainFEName(),domain_error); CHKERRQ(ierr);
      ierr = PoissonExample::AuxFunctions(m_field).assembleGhostVector(global_error); CHKERRQ(ierr);
      ierr = PoissonExample::AuxFunctions(m_field).printError(global_error); CHKERRQ(ierr);
      if(flg_test == PETSC_TRUE) {
        ierr = PoissonExample::AuxFunctions(m_field).testError(global_error); CHKERRQ(ierr);
      }
    }

    {
      // Loop over all elements in the mesh and for each execute post_proc_volume
      // element and operators on it.
      ierr = DMoFEMLoopFiniteElements(dm,simple_interface->getDomainFEName(),post_proc_volume); CHKERRQ(ierr);
      // Write results
      ierr = boost::static_pointer_cast<PostProcVolumeOnRefinedMesh>(post_proc_volume)->writeFile("out_vol.h5m"); CHKERRQ(ierr);
    }

    // Destroy DM, no longer needed.
    ierr = DMDestroy(&dm); CHKERRQ(ierr);

    // Destroy ghost vector
    ierr = VecDestroy(&global_error); CHKERRQ(ierr);

  } catch (MoFEMException const &e) {
    SETERRQ(PETSC_COMM_SELF,e.errorCode,e.errorMessage);
  }

  // finish work cleaning memory, getting statistics, etc.
  ierr = PetscFinalize(); CHKERRQ(ierr);

  return 0;

}
