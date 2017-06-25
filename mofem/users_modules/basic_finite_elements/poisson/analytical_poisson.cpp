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
 * \brief Set integration rule
 */
struct IntRule {
  int operator()(
    int order_row,int order_col,int order_data
  ) const { return 2*(order_data-1); }
};

/**
 * \brief Function
 */
struct exactFunction {
  double operator()(const double x,const double y,const double z) const {
    return 1+x*x+y*y+z*z*z;
  }
};

/**
 * \brief Laplacian of function
 */
struct exactLaplacianFunction {
  double operator()(const double x,const double y,const double z) const {
    return 0+2+2+3*2*z;
  }
};

/**
 * \brief Set integration rule to volume elements
 */
struct VolRule { int operator()(int order_row,int order_col,int order_data) const { return 2*(order_data-1); } };

/**
 * \brief Set integration rule to boundary elements
 */
struct FaceRule { int operator()(int order_row,int order_col,int order_data) const { return 2*order_data; } };

int main(int argc, char *argv[]) {

  ErrorCode rval;
  PetscErrorCode ierr;

  // initialize petsc
  PetscInitialize(&argc,&argv,(char *)0,help);

  int order = 3;  //< approximation order
  ierr = PetscOptionsBegin(
    PETSC_COMM_WORLD,"",
    "Poisson's problem options","none"
  ); CHKERRQ(ierr);
  ierr = PetscOptionsInt(
    "-order","approximation order","",order,&order,PETSC_NULL
  ); CHKERRQ(ierr);
  ierr = PetscOptionsEnd(); CHKERRQ(ierr);

  try {

    // Create MoAB database
    moab::Core moab_core;
    moab::Interface& moab = moab_core;

    // Create MoFEM database and link it to MoAB
    MoFEM::Core mofem_core(moab);
    MoFEM::Interface& m_field = mofem_core;

    // Register DM Manager
    DMType dm_name = "DMMOFEM";
    ierr = DMRegister_MoFEM(dm_name); CHKERRQ(ierr);

    // Simple interface
    Simple *simple_interface;
    ierr = m_field.query_interface(simple_interface); CHKERRQ(ierr);

    // Build problem
    {

      // Get options from command line
      ierr = simple_interface->getOptions(); CHKERRQ(ierr);
      // Load mesh file
      ierr = simple_interface->loadFile(); CHKERRQ(ierr);
      // Add fields
      ierr = simple_interface->addDomainField("U",H1,AINSWORTH_LEGENDRE_BASE,1); CHKERRQ(ierr);
      ierr = simple_interface->addBoundaryField("L",HDIV,DEMKOWICZ_JACOBI_BASE,1); CHKERRQ(ierr);
      // ierr = simple_interface->addBoundaryField("L",L2,AINSWORTH_LEGENDRE_BASE,1); CHKERRQ(ierr);
      ierr = simple_interface->addDataField("ERROR",L2,AINSWORTH_LEGENDRE_BASE,1); CHKERRQ(ierr);
      // Set fields order
      ierr = simple_interface->setFieldOrder("U",order); CHKERRQ(ierr);
      ierr = simple_interface->setFieldOrder("L",order-1); CHKERRQ(ierr);
      ierr = simple_interface->setFieldOrder("ERROR",0); CHKERRQ(ierr);
      // Setup problem
      ierr = simple_interface->setUp(); CHKERRQ(ierr);

    }

    DM dm;
    // Get dm
    ierr = simple_interface->getDM(&dm); CHKERRQ(ierr);

    // Create finite elements and data operators on entities
    boost::shared_ptr<FEMethod> null;
    boost::shared_ptr<ForcesAndSurcesCore> domain_lhs_fe;
    boost::shared_ptr<ForcesAndSurcesCore> boundary_lhs_fe;
    boost::shared_ptr<ForcesAndSurcesCore> domain_rhs_fe;
    boost::shared_ptr<ForcesAndSurcesCore> boundary_rhs_fe;
    {
      // Create elements
      domain_lhs_fe = boost::shared_ptr<ForcesAndSurcesCore>(new VolumeElementForcesAndSourcesCore(m_field));
      boundary_lhs_fe = boost::shared_ptr<ForcesAndSurcesCore>(new FaceElementForcesAndSourcesCore(m_field));
      domain_rhs_fe = boost::shared_ptr<ForcesAndSurcesCore>(new VolumeElementForcesAndSourcesCore(m_field));
      boundary_rhs_fe = boost::shared_ptr<ForcesAndSurcesCore>(new FaceElementForcesAndSourcesCore(m_field));
      // Set integration rule
      domain_lhs_fe->getRuleHook = VolRule();
      domain_rhs_fe->getRuleHook = VolRule();
      boundary_lhs_fe->getRuleHook = FaceRule();
      boundary_rhs_fe->getRuleHook = FaceRule();
      // Set operator to the volume element
      domain_lhs_fe->getOpPtrVector().push_back(new OpGradGrad());
      domain_rhs_fe->getOpPtrVector().push_back(new OpVF(exactLaplacianFunction()));
      // Set operator to the boundary element
      boundary_lhs_fe->getOpPtrVector().push_back(new OpLUHdiv(true));
      boundary_rhs_fe->getOpPtrVector().push_back(new OpLgHdiv(exactFunction()));
      // boost::static_pointer_cast<ForcesAndSurcesCore>(boundary_lhs_fe)->
      // getOpPtrVector().push_back(new OpLU(true));
      // boost::static_pointer_cast<ForcesAndSurcesCore>(boundary_rhs_fe)->
      // getOpPtrVector().push_back(new OpLg(exactFunction()));
    }

    // Set KSP context
    {
      // Set operators for KSP solver
      ierr = DMMoFEMKSPSetComputeOperators(
        dm,simple_interface->getDomainFEName(),domain_lhs_fe,null,null
      ); CHKERRQ(ierr);
      ierr = DMMoFEMKSPSetComputeOperators(
        dm,simple_interface->getBoundaryFEName(),boundary_lhs_fe,null,null
      ); CHKERRQ(ierr);
      // Set calculation of rhs for KSP solver
      ierr = DMMoFEMKSPSetComputeRHS(
        dm,simple_interface->getDomainFEName(),domain_rhs_fe,null,null
      ); CHKERRQ(ierr);
      ierr = DMMoFEMKSPSetComputeRHS(
        dm,simple_interface->getBoundaryFEName(),boundary_rhs_fe,null,null
      ); CHKERRQ(ierr);
    }

    // Solve problem
    {
      Vec F,D;
      ierr = DMCreateGlobalVector_MoFEM(dm,&F); CHKERRQ(ierr);
      ierr = VecDuplicate(F,&D); CHKERRQ(ierr);

      // Create solver
      KSP solver;
      ierr = KSPCreate(PETSC_COMM_WORLD,&solver); CHKERRQ(ierr);
      // ierr = KSPSetDMActive(solver,PETSC_TRUE); CHKERRQ(ierr);
      ierr = KSPSetFromOptions(solver); CHKERRQ(ierr);
      // ierr = KSPSetOperators(solver,A,A); CHKERRQ(ierr);
      ierr = KSPSetDM(solver,dm); CHKERRQ(ierr);

      ierr = KSPSetUp(solver); CHKERRQ(ierr);
      ierr = KSPSolve(solver,F,D); CHKERRQ(ierr);

      // Scatter solution on the mesh
      ierr = DMoFEMMeshToGlobalVector(dm,D,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);

      // Clean data
      ierr = KSPDestroy(&solver); CHKERRQ(ierr);
      ierr = VecDestroy(&D); CHKERRQ(ierr);
      ierr = VecDestroy(&F); CHKERRQ(ierr);
    }

    // Calculate error
    {
      boost::shared_ptr<ForcesAndSurcesCore> domain_error;
      domain_error = boost::shared_ptr<ForcesAndSurcesCore>(new VolumeElementForcesAndSourcesCore(m_field));
      domain_error->getRuleHook = VolRule();
      boost::shared_ptr<VectorDouble> values_at_integation_ptrs = boost::make_shared<VectorDouble>();
      domain_error->getOpPtrVector().push_back(new OpCalculateScalarFieldValues("U",values_at_integation_ptrs));
      domain_error->getOpPtrVector().push_back(new OpErrorL2(exactFunction(),values_at_integation_ptrs));
      ierr = DMoFEMLoopFiniteElements(dm,simple_interface->getDomainFEName(),domain_error); CHKERRQ(ierr);
    }

    // Post-process results
    {
      PostProcVolumeOnRefinedMesh post_proc_volume(m_field);
      // Add operators to the elements, starting with some generic
      ierr = post_proc_volume.generateReferenceElementMesh(); CHKERRQ(ierr);
      ierr = post_proc_volume.addFieldValuesPostProc("U"); CHKERRQ(ierr);
      ierr = post_proc_volume.addFieldValuesPostProc("ERROR"); CHKERRQ(ierr);
      ierr = post_proc_volume.addFieldValuesGradientPostProc("U"); CHKERRQ(ierr);
      ierr = DMoFEMLoopFiniteElements(
        dm,simple_interface->getDomainFEName().c_str(),&post_proc_volume
      ); CHKERRQ(ierr);
      ierr = post_proc_volume.writeFile("out_vol.h5m"); CHKERRQ(ierr);
    }

    // Destroy DM and problem with it
    ierr = DMDestroy(&dm); CHKERRQ(ierr);

  } catch (MoFEMException const &e) {
    SETERRQ(PETSC_COMM_SELF,e.errorCode,e.errorMessage);
  }

  // finish work cleaning memory, getting statistics, etc.
  ierr = PetscFinalize(); CHKERRQ(ierr);

  return 0;

}
