/**
 * \file analytical_poisson.cpp
 * \ingroup mofem_simple_interface
 * \example analytical_poisson.cpp
 *
 * This is similar test to one form
 * <https://fenicsproject.org/pub/tutorial/html/._ftut1004.html#ch:fundamentals>,
 * but exploiting MoFEM capabilities.
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
#include <PoissonOperators.hpp>
using namespace PoissonOperators;

static char help[] = "...\n\n";

/**
 * \brief Function
 *
 * This is exact function. Finite element method is used to find This
 * function in the body volume. If this function is given by polynomial
 * order "p" and order of approximation is "p" or higher, solution of
 * finite element method is exact (with numerical precision)
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
 * \brief Laplacian of function
 *
 * This is laplacian of \f$u\f$
 *
 */
struct ExactLaplacianFunction {
  double operator()(const double x,const double y,const double z) const {
    return 0+2+2+3*2*z;
  }
};

/**
 * \brief Set integration rule to volume elements
 *
 * This rule is used to integrate \f$\nabla v \cdot \nabla u\f$, thus
 * if approximation field and testing field are polynomial order "p",
 * then rule for exact integration is 2*(p-1)
 *
 */
struct VolRule {
  int operator()(int,int,int p) const {
    // cerr << p << endl;
    return 2*p;
  }
};

/**
 * \brief Set integration rule to boundary elements
 *
 * This is uses to integrate values on the face. Is used to integrate
 * \f$(\mathbf{n} \cdot \lambda) u\f$, where Lagrange multiplayer
 * is order "p_row" and approximate function is order "p_col".
 *
 */
struct FaceRule {
  int operator()(int p_row,int p_col,int p_data) const {
    // cerr << p_row << " " << p_col << " " << p_data << endl;
    return p_row+p_col;
  }
};

int main(int argc, char *argv[]) {

  ErrorCode rval;
  PetscErrorCode ierr;

  // initialize petsc
  PetscInitialize(&argc,&argv,(char *)0,help);

  // Get command line options

  int order = 3;  // default approximation order
  PetscBool flg_test = PETSC_FALSE; // true check if error is numerical error
  ierr = PetscOptionsBegin(
    PETSC_COMM_WORLD,"",
    "Poisson's problem options","none"
  ); CHKERRQ(ierr);
  // Set approximation order
  ierr = PetscOptionsInt(
    "-order","approximation order","",order,&order,PETSC_NULL
  ); CHKERRQ(ierr);
  // Set testing (used by ctest)
  ierr = PetscOptionsBool(
    "-test","if true is ctest","",flg_test,&flg_test,PETSC_NULL
  ); CHKERRQ(ierr);
  ierr = PetscOptionsEnd(); CHKERRQ(ierr);

  try {

    // Create MoAB database
    moab::Core moab_core; // create database
    moab::Interface& moab = moab_core; // craete interface to database

    // Create MoFEM database and link it to MoAB
    MoFEM::Core mofem_core(moab); // create database
    MoFEM::Interface& m_field = mofem_core; // create interface to database

    // Register DM Manager
    DMType dm_name = "DMMOFEM"; // name of new DM manager
    ierr = DMRegister_MoFEM(dm_name); CHKERRQ(ierr); // register MoFEM DM in PETSc

    // Get simple interface. MoFEM::Interface is build to solve complex General
    // problems. Simple interface is simplified version enabling quick and
    // easy construction of problem.
    Simple *simple_interface;
    // Query interface and get pointer to Simple interface
    ierr = m_field.query_interface(simple_interface); CHKERRQ(ierr);

    // Build problem with simple ineterface
    {

      // Get options for simple interface from command line
      ierr = simple_interface->getOptions(); CHKERRQ(ierr);
      // Load mesh file to database
      ierr = simple_interface->loadFile(); CHKERRQ(ierr);

      // Add field on domain and boundary. Field is declared by space and base and rank. space
      // can be H1. Hcurl, Hdiv and L2, base can be AINSWORTH_LEGENDRE_BASE, DEMKOWICZ_JACOBI_BASE and more,
      // where rank is number of coefficients for dof.
      //
      // Simple interface assumes that there are four types of field, woman feels, defined
      // on body domain, boundary fields defined on body boundary, skeleton field defined
      // on finite element skeleton. Finally data field is defined on body domain. Data field
      // is not solved but used for Post-process or to keep material parameters, etc.
      //
      ierr = simple_interface->addDomainField("U",H1,AINSWORTH_LEGENDRE_BASE,1); CHKERRQ(ierr);
      ierr = simple_interface->addBoundaryField("L",H1,AINSWORTH_LEGENDRE_BASE,1); CHKERRQ(ierr);
      // Add error (data) field
      ierr = simple_interface->addDataField("ERROR",L2,AINSWORTH_LEGENDRE_BASE,1); CHKERRQ(ierr);

      // Set fields order domain and boundary fields.
      ierr = simple_interface->setFieldOrder("U",order); CHKERRQ(ierr); // to approximate function
      ierr = simple_interface->setFieldOrder("L",order); CHKERRQ(ierr); // to Lagrange multiplayers
      ierr = simple_interface->setFieldOrder("ERROR",0); CHKERRQ(ierr); // approximation order for error

      // Setup problem. At that point database is constructed.
      ierr = simple_interface->setUp(); CHKERRQ(ierr);

    }

    // Create finite elements and data operators on entities
    boost::shared_ptr<FEMethod> null;                       // This is used to pass null element
    boost::shared_ptr<ForcesAndSurcesCore> domain_lhs_fe;   // Domain element used to calculate matrix
    boost::shared_ptr<ForcesAndSurcesCore> boundary_lhs_fe; // Domain element to evaluate the right hand side
    boost::shared_ptr<ForcesAndSurcesCore> domain_rhs_fe;   // Boundary element to evaluate boundary constrains matrix
    boost::shared_ptr<ForcesAndSurcesCore> boundary_rhs_fe; // Boundary element to evaluate constrains vector
    {
      // Create elements element instances
      domain_lhs_fe = boost::shared_ptr<ForcesAndSurcesCore>(new VolumeElementForcesAndSourcesCore(m_field));
      boundary_lhs_fe = boost::shared_ptr<ForcesAndSurcesCore>(new FaceElementForcesAndSourcesCore(m_field));
      domain_rhs_fe = boost::shared_ptr<ForcesAndSurcesCore>(new VolumeElementForcesAndSourcesCore(m_field));
      boundary_rhs_fe = boost::shared_ptr<ForcesAndSurcesCore>(new FaceElementForcesAndSourcesCore(m_field));
      // Set integration rule to elements instances
      domain_lhs_fe->getRuleHook = VolRule();
      domain_rhs_fe->getRuleHook = VolRule();
      boundary_lhs_fe->getRuleHook = FaceRule();
      boundary_rhs_fe->getRuleHook = FaceRule();

      // Ass operators to element instances
      // Add operator grad-grad for calualte matrix
      domain_lhs_fe->getOpPtrVector().push_back(new OpGradGrad());
      // Add operator to calculate source terms
      domain_rhs_fe->getOpPtrVector().push_back(new OpVF(ExactLaplacianFunction()));
      // Add operator calculating constrains matrix
      boundary_lhs_fe->getOpPtrVector().push_back(new OpLU(true));
      // Add operator calculating constrains values
      boundary_rhs_fe->getOpPtrVector().push_back(new OpLg(ExactFunction()));
    }

    // Get access to PETSC-MoFEM DM manager. F
    // or more derails see <http://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/DM/index.html>
    DM dm;
    // Get dm
    ierr = simple_interface->getDM(&dm); CHKERRQ(ierr);

    // Set KSP context for DM. At that point only elements are added to DM operators.
    // Calculations of matrices and vectors is executed by KSP solver.
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

    // Solve problem
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
      // Create ghost vector to assemble errors from all element on distributed mesh.
      // Ghost vector has size 1, where one element is owned by processor 0, other processor
      // have one ghost element of zero element at processor 0.
      int ghosts[] = { 0 };
      Vec global_error;
      ierr = VecCreateGhost(
        PETSC_COMM_WORLD,m_field.get_comm_rank()==0?1:0,1,m_field.get_comm_rank()==0?0:1,ghosts,&global_error
      ); CHKERRQ(ierr);
      // Create finite element instance to calualte error
      boost::shared_ptr<ForcesAndSurcesCore> domain_error;
      domain_error = boost::shared_ptr<ForcesAndSurcesCore>(new VolumeElementForcesAndSourcesCore(m_field));
      // Set integration rule
      domain_error->getRuleHook = VolRule();
      // Crate shared vector storing values of field "u" on integration points on element. element
      // is local and is used to exchange data between operators.
      boost::shared_ptr<VectorDouble> values_at_integation_ptrs = boost::make_shared<VectorDouble>();
      // Add default operator to calculate field values at integration points
      domain_error->getOpPtrVector().push_back(new OpCalculateScalarFieldValues("U",values_at_integation_ptrs));
      // Add operator to integrate error element by element.
      domain_error->getOpPtrVector().push_back(new OpErrorL2(ExactFunction(),values_at_integation_ptrs,global_error));
      // Loop over all elements in mesh, and run users operators on each element.
      ierr = DMoFEMLoopFiniteElements(dm,simple_interface->getDomainFEName(),domain_error); CHKERRQ(ierr);
      // Assemble vector with globe error.
      ierr = VecAssemblyEnd(global_error); CHKERRQ(ierr);
      ierr = VecGhostUpdateBegin(global_error,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
      ierr = VecGhostUpdateEnd(global_error,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
      // Print error
      if(m_field.get_comm_rank()==0) {
        double *e;
        ierr = VecGetArray(global_error,&e); CHKERRQ(ierr);
        cout << "Global errror " << e[0] << endl;
        if(flg_test) {
          // Check if error is zero, otherwise throw error
          if(e[0]>1e-10||e[0]!=e[0]) {
            SETERRQ(PETSC_COMM_SELF,MOFEM_ATOM_TEST_INVALID,"Test failed, error too big");
          }
        }
        ierr = VecRestoreArray(global_error,&e); CHKERRQ(ierr);
      }
      // Destroy ghost vector
      ierr = VecDestroy(&global_error); CHKERRQ(ierr);
    }

    // Post-process results. This is standard element, with functionality
    // enabling refining mesh for post-processing. In addition in PostProcOnRefMesh.hpp
    // are implanted set of  users operators to post-processing fields. Here
    // using simplified mechanism for post-processing finite element, we
    // add operators to save data from field on mesh tags for preview
    // visualization.
    {
      PostProcVolumeOnRefinedMesh post_proc_volume(m_field);
      // Add operators to the elements, starting with some generic
      ierr = post_proc_volume.generateReferenceElementMesh(); CHKERRQ(ierr);
      ierr = post_proc_volume.addFieldValuesPostProc("U"); CHKERRQ(ierr);
      ierr = post_proc_volume.addFieldValuesPostProc("ERROR"); CHKERRQ(ierr);
      ierr = post_proc_volume.addFieldValuesGradientPostProc("U"); CHKERRQ(ierr);
      // Loop over all elements in the mesh and for each execute post_proc_volume
      // element and operators on it.
      ierr = DMoFEMLoopFiniteElements(
        dm,simple_interface->getDomainFEName().c_str(),&post_proc_volume
      ); CHKERRQ(ierr);
      ierr = post_proc_volume.writeFile("out_vol.h5m"); CHKERRQ(ierr);
    }

    // Destroy DM, no longer needed.
    ierr = DMDestroy(&dm); CHKERRQ(ierr);

  } catch (MoFEMException const &e) {
    SETERRQ(PETSC_COMM_SELF,e.errorCode,e.errorMessage);
  }

  // finish work cleaning memory, getting statistics, etc.
  ierr = PetscFinalize(); CHKERRQ(ierr);

  return 0;

}
