/**
 * \file analytical_poisson_field_split.cpp
 * \ingroup mofem_simple_interface
 * \example analytical_poisson_field_split.cpp
 *
 * For more information and detailed explain of this
 * example see \ref poisson_tut3
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
static const bool debug = false;

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

struct OpS: public FaceElementForcesAndSourcesCore::UserDataOperator {

  OpS(const bool beta = 1):
  FaceElementForcesAndSourcesCore::UserDataOperator("U","U",OPROWCOL,true),
  bEta(beta) {
  }

  /**
   * \brief Do calculations for give operator
   * @param  row_side row side number (local number) of entity on element
   * @param  col_side column side number (local number) of entity on element
   * @param  row_type type of row entity MBVERTEX, MBEDGE, MBTRI or MBTET
   * @param  col_type type of column entity MBVERTEX, MBEDGE, MBTRI or MBTET
   * @param  row_data data for row
   * @param  col_data data for column
   * @return          error code
   */
  PetscErrorCode doWork(
    int row_side,int col_side,
    EntityType row_type,EntityType col_type,
    DataForcesAndSourcesCore::EntData &row_data,DataForcesAndSourcesCore::EntData &col_data
  ) {
    MoFEMFunctionBeginHot;
    // get number of dofs on row
    nbRows = row_data.getIndices().size();
    // if no dofs on row, exit that work, nothing to do here
    if(!nbRows) MoFEMFunctionReturnHot(0);
    // get number of dofs on column
    nbCols = col_data.getIndices().size();
    // if no dofs on Columbia, exit nothing to do here
    if(!nbCols) MoFEMFunctionReturnHot(0);
    // get number of integration points
    nbIntegrationPts = getGaussPts().size2();
    // check if entity block is on matrix diagonal
    if(
      row_side==col_side&&
      row_type==col_type
    ) {
      isDiag = true; // yes, it is
    } else {
      isDiag = false;
    }
    // integrate local matrix for entity block
    ierr = iNtegrate(row_data,col_data); CHKERRQ(ierr);
    // assemble local matrix
    ierr = aSsemble(row_data,col_data); CHKERRQ(ierr);
    MoFEMFunctionReturnHot(0);
  }

private:

  const double bEta;

    ///< error code

  int nbRows;           ///< number of dofs on rows
  int nbCols;           ///< number if dof on column
  int nbIntegrationPts; ///< number of integration points
  bool isDiag;          ///< true if this block is on diagonal

  FTensor::Index<'i',3> i;  ///< summit Index
  MatrixDouble locMat;      ///< local entity block matrix

  /**
   * \brief Integrate grad-grad operator
   * @param  row_data row data (consist base functions on row entity)
   * @param  col_data column data (consist base functions on column entity)
   * @return          error code
   */
  inline PetscErrorCode iNtegrate(
    DataForcesAndSourcesCore::EntData &row_data,DataForcesAndSourcesCore::EntData &col_data
  ) {
    MoFEMFunctionBeginHot;
    // set size of local entity bock
    locMat.resize(nbRows,nbCols,false);
    // clear matrix
    locMat.clear();
    // get element area
    double area = getArea()*bEta;
    // get integration weights
    FTensor::Tensor0<double*> t_w = getFTensor0IntegrationWeight();
    // get base function gradient on rows
    FTensor::Tensor0<double*> t_row_base = row_data.getFTensor0N();
    // loop over integration points
    for(int gg = 0;gg!=nbIntegrationPts;gg++) {
      // take into account Jacobean
      const double alpha = t_w*area;
      // take fist element to local matrix
      FTensor::Tensor0<double*> a(&*locMat.data().begin());
      // loop over rows base functions
      for(int rr = 0;rr!=nbRows;rr++) {
        // get column base functions gradient at gauss point gg
        FTensor::Tensor0<double*> t_col_base = col_data.getFTensor0N(gg,0);
        // loop over columns
        for(int cc = 0;cc!=nbCols;cc++) {
          // calculate element of local matrix
          a += alpha*t_row_base*t_col_base;
          ++t_col_base; // move to another gradient of base function on column
          ++a;  // move to another element of local matrix in column
        }
        ++t_row_base; // move to another element of gradient of base function on row
      }
      ++t_w; // move to another integration weight
    }
    MoFEMFunctionReturnHot(0);
  }

  /**
   * \brief Assemble local entity block matrix
   * @param  row_data row data (consist base functions on row entity)
   * @param  col_data column data (consist base functions on column entity)
   * @return          error code
   */
  inline PetscErrorCode aSsemble(
    DataForcesAndSourcesCore::EntData &row_data,DataForcesAndSourcesCore::EntData &col_data
  ) {
    MoFEMFunctionBeginHot;
    // get pointer to first global index on row
    const int* row_indices = &*row_data.getIndices().data().begin();
    // get pointer to first global index on column
    const int* col_indices = &*col_data.getIndices().data().begin();
    Mat B = getFEMethod()->ksp_B!=PETSC_NULL? getFEMethod()->ksp_B : getFEMethod()->snes_B;
    // assemble local matrix
    ierr = MatSetValues(
      B, nbRows,row_indices,nbCols,col_indices,&*locMat.data().begin(),ADD_VALUES
    ); CHKERRQ(ierr);
    if(!isDiag) {
      // if not diagonal term and since global matrix is symmetric assemble
      // transpose term.
      locMat = trans(locMat);
      ierr = MatSetValues(
        B, nbCols,col_indices,nbRows,row_indices,&*locMat.data().begin(),ADD_VALUES
      ); CHKERRQ(ierr);
    }
    MoFEMFunctionReturnHot(0);
  }

};

int main(int argc, char *argv[]) {

  // 
  

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
    boost::shared_ptr<ForcesAndSourcesCore> domain_lhs_fe;     ///< Volume element for the matrix
    boost::shared_ptr<ForcesAndSourcesCore> boundary_lhs_fe;   ///< Boundary element for the matrix
    boost::shared_ptr<ForcesAndSourcesCore> domain_rhs_fe;     ///< Volume element to assemble vector
    boost::shared_ptr<ForcesAndSourcesCore> boundary_rhs_fe;   ///< Volume element to assemble vector
    boost::shared_ptr<ForcesAndSourcesCore> domain_error;      ///< Volume element evaluate error
    boost::shared_ptr<ForcesAndSourcesCore> post_proc_volume;  ///< Volume element to Post-process results
    boost::shared_ptr<ForcesAndSourcesCore> null;              ///< Null element do nothing
    boost::shared_ptr<ForcesAndSourcesCore> boundary_penalty_lhs_fe;
    {
      // Add problem specific operators the generic finite elements to calculate matrices and vectors.
      ierr = PoissonExample::CreateFiniteElements(m_field).createFEToAssembleMatrixAndVector(
        ExactFunction(),ExactLaplacianFunction(),
        domain_lhs_fe,boundary_lhs_fe,domain_rhs_fe,boundary_rhs_fe,false
      ); CHKERRQ(ierr);
      // Add problem specific operators the generic finite elements to calculate error on elements and global error
      // in H1 norm
      ierr = PoissonExample::CreateFiniteElements(m_field).createFEToEvaluateError(
        ExactFunction(),ExactFunctionGrad(),global_error,domain_error
      ); CHKERRQ(ierr);
      // Post-process results
      ierr = PoissonExample::CreateFiniteElements(m_field).creatFEToPostProcessResults(post_proc_volume); CHKERRQ(ierr);

      const double beta = 1;
      boundary_penalty_lhs_fe = boost::shared_ptr<ForcesAndSourcesCore>(new FaceElementForcesAndSourcesCore(m_field));
      boundary_penalty_lhs_fe->getRuleHook = PoissonExample::FaceRule();
      boundary_penalty_lhs_fe->getOpPtrVector().push_back(new OpS(beta));
      boundary_rhs_fe->getOpPtrVector().push_back(new PoissonExample::Op_g(ExactFunction(),"U",beta));

    }

    // Get simple interface is simplified version enabling quick and
    // easy construction of problem.
    Simple *simple_interface;
    // Query interface and get pointer to Simple interface
    ierr = m_field.getInterface(simple_interface); CHKERRQ(ierr);

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
      ierr = simple_interface->setFieldOrder("L",order); CHKERRQ(ierr); // to Lagrange multipliers
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

    // Solve problem, only PETEc interface here.
    {

      // Create the right hand side vector and vector of unknowns
      Vec F,D;
      ierr = DMCreateGlobalVector(dm,&F); CHKERRQ(ierr);
      // Create unknown vector by creating duplicate copy of F vector. only
      // structure is duplicated no values.
      ierr = VecDuplicate(F,&D); CHKERRQ(ierr);

      DM dm_sub_KK,dm_sub_LU;
      ublas::matrix<Mat> nested_matrices(2,2);
      ublas::vector<IS> nested_is(2);

      ierr = DMCreate(PETSC_COMM_WORLD,&dm_sub_KK);CHKERRQ(ierr);
      ierr = DMSetType(dm_sub_KK,"DMMOFEM");CHKERRQ(ierr);
      ierr = DMMoFEMCreateSubDM(dm_sub_KK,dm,"SUB_KK"); CHKERRQ(ierr);
      ierr = DMSetFromOptions(dm_sub_KK); CHKERRQ(ierr);
      ierr = DMMoFEMSetSquareProblem(dm_sub_KK,PETSC_TRUE); CHKERRQ(ierr);
      ierr = DMMoFEMAddSubFieldRow(dm_sub_KK,"U"); CHKERRQ(ierr);
      ierr = DMMoFEMAddSubFieldCol(dm_sub_KK,"U"); CHKERRQ(ierr);
      ierr = DMMoFEMAddElement(dm_sub_KK,simple_interface->getDomainFEName().c_str()); CHKERRQ(ierr);
      ierr = DMMoFEMAddElement(dm_sub_KK,simple_interface->getBoundaryFEName().c_str()); CHKERRQ(ierr);
      ierr = DMSetUp(dm_sub_KK); CHKERRQ(ierr);
      ierr = DMMoFEMGetSubRowIS(dm_sub_KK,&nested_is[0]); CHKERRQ(ierr);
      ierr = DMCreateMatrix(dm_sub_KK,&nested_matrices(0,0)); CHKERRQ(ierr);
      domain_lhs_fe->ksp_B = nested_matrices(0,0);
      ierr = DMoFEMLoopFiniteElements(dm_sub_KK,simple_interface->getDomainFEName(),domain_lhs_fe); CHKERRQ(ierr);
      boundary_penalty_lhs_fe->ksp_B = nested_matrices(0,0);
      ierr = DMoFEMLoopFiniteElements(dm_sub_KK,simple_interface->getBoundaryFEName(),boundary_penalty_lhs_fe); CHKERRQ(ierr);
      ierr = MatAssemblyBegin(nested_matrices(0,0),MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
      ierr = MatAssemblyEnd(nested_matrices(0,0),MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
      ierr = DMDestroy(&dm_sub_KK); CHKERRQ(ierr);
      //
      ierr = DMCreate(PETSC_COMM_WORLD,&dm_sub_LU);CHKERRQ(ierr);
      ierr = DMSetType(dm_sub_LU,"DMMOFEM");CHKERRQ(ierr);
      ierr = DMMoFEMCreateSubDM(dm_sub_LU,dm,"SUB_LU"); CHKERRQ(ierr);
      ierr = DMSetFromOptions(dm_sub_LU); CHKERRQ(ierr);
      ierr = DMMoFEMSetSquareProblem(dm_sub_LU,PETSC_FALSE); CHKERRQ(ierr);
      ierr = DMMoFEMAddSubFieldRow(dm_sub_LU,"L"); CHKERRQ(ierr);
      ierr = DMMoFEMAddSubFieldCol(dm_sub_LU,"U"); CHKERRQ(ierr);
      ierr = DMMoFEMAddElement(dm_sub_LU,simple_interface->getBoundaryFEName().c_str()); CHKERRQ(ierr);
      ierr = DMSetUp(dm_sub_LU); CHKERRQ(ierr);
      ierr = DMMoFEMGetSubRowIS(dm_sub_LU,&nested_is[1]); CHKERRQ(ierr);
      ierr = DMCreateMatrix(dm_sub_LU,&nested_matrices(1,0)); CHKERRQ(ierr);
      boundary_lhs_fe->ksp_B = nested_matrices(1,0);
      ierr = DMoFEMLoopFiniteElements(dm_sub_LU,simple_interface->getBoundaryFEName(),boundary_lhs_fe); CHKERRQ(ierr);
      ierr = MatAssemblyBegin(nested_matrices(1,0),MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
      ierr = MatAssemblyEnd(nested_matrices(1,0),MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
      // ierr = MatCreateTranspose(nested_matrices(1,0),&nested_matrices(0,1)); CHKERRQ(ierr);
      ierr = MatTranspose(nested_matrices(1,0),MAT_INITIAL_MATRIX,&nested_matrices(0,1)); CHKERRQ(ierr);
      ierr = DMDestroy(&dm_sub_LU); CHKERRQ(ierr);

      domain_rhs_fe->ksp_f = F;
      ierr = DMoFEMLoopFiniteElements(dm,simple_interface->getDomainFEName(),domain_rhs_fe); CHKERRQ(ierr);
      boundary_rhs_fe->ksp_f = F;
      ierr = DMoFEMLoopFiniteElements(dm,simple_interface->getBoundaryFEName(),boundary_rhs_fe); CHKERRQ(ierr);
      ierr = VecAssemblyBegin(F); CHKERRQ(ierr);
      ierr = VecAssemblyEnd(F); CHKERRQ(ierr);

      Mat A;
      nested_matrices(1,1) = PETSC_NULL;

      if(debug) {
        MatType type;
        MatGetType(nested_matrices(0,0),&type);
        cerr << "K " << type << endl;
        MatGetType(nested_matrices(1,0),&type);
        cerr << "C " << type << endl;
        MatGetType(nested_matrices(0,1),&type);
        cerr << "CT " << type << endl;
        std::string wait;
        cerr << "UU" << endl;
        MatView(nested_matrices(0,0),PETSC_VIEWER_DRAW_WORLD);
        std::cin >> wait;
        cerr << "LU" << endl;
        MatView(nested_matrices(1,0),PETSC_VIEWER_DRAW_WORLD);
        std::cin >> wait;
        cerr << "UL" << endl;
        MatView(nested_matrices(0,1),PETSC_VIEWER_DRAW_WORLD);
        std::cin >> wait;
      }

      ierr = MatCreateNest(
        PETSC_COMM_WORLD,
        2,&nested_is[0],2,&nested_is[0],&nested_matrices(0,0),&A
      ); CHKERRQ(ierr);

      // Create solver and link it to DM
      KSP solver;
      ierr = KSPCreate(PETSC_COMM_WORLD,&solver); CHKERRQ(ierr);
      ierr = KSPSetFromOptions(solver); CHKERRQ(ierr);
      // Set operators
      ierr = KSPSetOperators(solver,A,A); CHKERRQ(ierr);
      PC pc;
      ierr = KSPGetPC(solver,&pc); CHKERRQ(ierr);
      PetscBool is_pcfs = PETSC_FALSE;
      PetscObjectTypeCompare((PetscObject)pc,PCFIELDSPLIT,&is_pcfs);
      if(is_pcfs) {
        ierr = PCFieldSplitSetIS(pc,NULL,nested_is[0]); CHKERRQ(ierr);
        ierr = PCFieldSplitSetIS(pc,NULL,nested_is[1]); CHKERRQ(ierr);
      } else {
        SETERRQ(
          PETSC_COMM_WORLD,
          MOFEM_DATA_INCONSISTENCY,
          "Only works with pre-conditioner PCFIELDSPLIT"
        );
      }
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
      for(int i = 0;i!=2;i++) {
        ierr = ISDestroy(&nested_is[i]); CHKERRQ(ierr);
        for(int j = 0;j!=2;j++) {
          ierr = MatDestroy(&nested_matrices(i,j)); CHKERRQ(ierr);
        }
      }
      ierr = MatDestroy(&A); CHKERRQ(ierr);
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
