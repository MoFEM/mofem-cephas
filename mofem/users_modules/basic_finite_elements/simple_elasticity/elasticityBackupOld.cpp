/** \file elasticity.cpp
 * \ingroup nonlinear_elastic_elem
 * \example elasticity.cpp

 The example shows how to solve the linear elastic problem. An example can read
 file with temperature field, then thermal stresses are included.

 What example can do:
 - take into account temperature field, i.e. calculate thermal stresses and deformation
 - stationary and time depend field is considered
 - take into account gravitational body forces
 - take in account fluid pressure
 - can work with higher order geometry definition
 - works on distributed meshes
 - multi-grid solver where each grid level is approximation order level
 - each mesh block can have different material parameters and approximation order

See example how code can be used \cite jordi:2017,
 \image html SquelaDamExampleByJordi.png "Example what you can do with this code. Analysis of the arch dam of Susqueda, located in Catalonia (Spain)" width=800px

 This is an example of application code; it does not show how elements are implemented. Example presents how to:
 - read mesh
 - set-up problem
 - run finite elements on the problem
 - assemble matrices and vectors
 - solve the linear system of equations
 - save results


 If you like to see how to implement finite elements, material, are other parts of the code, look here;
 - Hooke material, see \ref Hooke
 - Thermal-stress assembly, see \ref  ThermalElement
 - Body forces element, see \ref BodyFroceConstantField
 - Fluid pressure element, see \ref FluidPressure
 - The general implementation of an element for arbitrary Lagrangian-Eulerian
 aformulation for a nonlinear elastic problem is here \ref
 NonlinearElasticElement. Here we limit ourselves to Hooke equation and fix
 mesh, so the problem becomes linear. Not that elastic element is implemented
 with automatic differentiation.

*/

/* MoFEM is free software: you can redistribute it and/or modify it under
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

#define TOL 1.e-6

#include <BasicFiniteElements.hpp>
using namespace MoFEM;

#include <boost/program_options.hpp>
using namespace std;
namespace po = boost::program_options;

#include <Hooke.hpp>
#include <ElasticityNewOperator.hpp>
#include "../poisson/src/PoissonOperators.hpp"
#include "../poisson/src/AuxPoissonFunctions.hpp"

using namespace boost::numeric;

static char help[] =
  "-my_block_config set block data\n"
  "\n";

//const double young_modulus = 1;
//const double poisson_ratio = 0.0;

struct BlockOptionData {
  int oRder;
  double yOung;
  double pOisson;
  double initTemp;
  BlockOptionData():
    oRder(-1),
    yOung(-1),
    pOisson(-2),
    initTemp(0) {}
};

struct OpK: public VolumeElementForcesAndSourcesCore::UserDataOperator {
 
  MatrixDouble rowB;
  MatrixDouble colB;
  MatrixDouble CB;
  MatrixDouble K, transK;

  MatrixDouble D;
  double yOung;
  double pOisson;
  double coefficient;

  OpK(bool symm = true):
    VolumeElementForcesAndSourcesCore::UserDataOperator("U","U",OPROWCOL,symm){
      
      pOisson = 0.1;
      yOung = 10;
      coefficient = yOung / ( (1 + pOisson) * (1 - 2 * pOisson));
      D.resize(6, 6 ,false);
      D.clear();  

      D.insert_element (0, 0, 1 - pOisson);  
      D.insert_element (1, 1, 1 - pOisson);  
      D.insert_element (2, 2, 1 - pOisson);  
      D.insert_element (3, 3,(1 - 2* pOisson));  
      D.insert_element (4, 4,(1 - 2* pOisson));  
      D.insert_element (5, 5,(1 - 2* pOisson));  

      D.insert_element (0, 1, pOisson);  
      D.insert_element (0, 2, pOisson);  
      D.insert_element (1, 0, pOisson);  

      D.insert_element (1, 2, pOisson);  
      D.insert_element (2, 0, pOisson);  
      D.insert_element (2, 1, pOisson);  

      D *= coefficient;
    } 




 PetscErrorCode makeB(const MatrixAdaptor& diffN,
		       MatrixDouble& B 
		       ){

   PetscFunctionBegin;
   unsigned int nb_dofs = diffN.size1();
   B.resize(6,3*nb_dofs,false);
   B.clear();
   for(unsigned int dd = 0;dd<nb_dofs;dd++) {
     const double diff[] = {
       diffN(dd,0), diffN(dd,1), diffN(dd,2)
     };
     const int dd3 = 3*dd;
     for(int rr = 0;rr<3;rr++) {
       B(rr,dd3+rr) = diff[rr];
     }
     // gamma_xy
     B(3,dd3+0) = diff[1];
     B(3,dd3+1) = diff[0];
     // gamma_yz
     B(4,dd3+1) = diff[2];
     B(4,dd3+2) = diff[1];
     // gamma_xz
     B(5,dd3+0) = diff[2];
     B(5,dd3+2) = diff[0];
   }

      // int counter = 0;
      // for(int i = 0; i != 6; ++i){
      // 	for(unsigned int j = 0; j != 3*nb_dofs; ++j){
      // 	  if (counter == i) printf("%.2f ", B(i,j)  );
      // 	  else {
      // 	    printf("\n%.2f ",B(i,j));
      // 	    counter ++;
      // 	  }
      // 	}
      // }
       
   PetscFunctionReturn(0);

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
      PetscFunctionBegin;
      // get number of dofs on row
      nbRows = row_data.getIndices().size();
      // if no dofs on row, exit that work, nothing to do here
      if(!nbRows) PetscFunctionReturn(0);
      // get number of dofs on column
      nbCols = col_data.getIndices().size();
      // if no dofs on Columbia, exit nothing to do here
      if(!nbCols) PetscFunctionReturn(0);
      // get number of integration points
      nbIntegrationPts = getGaussPts().size2();
      // chekk if entity block is on matrix diagonal
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
      // asseble local matrix
      ierr = aSsemble(row_data,col_data); CHKERRQ(ierr);
      PetscFunctionReturn(0);
    }

  protected:

      ///< error code

    int nbRows;           ///< number of dofs on rows
    int nbCols;           ///< number if dof on column
    int nbIntegrationPts; ///< number of integration points
    bool isDiag;          ///< true if this block is on diagonal
    
    /**
     * \brief Integrate grad-grad operator
     * @param  row_data row data (consist base functions on row entity)
     * @param  col_data column data (consist base functions on column entity)
     * @return          error code
     */
   virtual PetscErrorCode iNtegrate(
      DataForcesAndSourcesCore::EntData &row_data,DataForcesAndSourcesCore::EntData &col_data
    ) {
      PetscFunctionBegin;

      int nb_dofs_row = row_data.getFieldData().size();
      if(nb_dofs_row==0) PetscFunctionReturn(0);
      int nb_dofs_col = col_data.getFieldData().size();
      if(nb_dofs_col==0) PetscFunctionReturn(0);
      
      K.resize(nb_dofs_row,nb_dofs_col,false);
      K.clear();
      CB.resize(6,nb_dofs_col,false);
      

      for(int gg = 0;gg!=nbIntegrationPts;gg++) {
	// get element volume
	// get integration weight
	double val = getVolume()*getGaussPts()(3,gg);

	const MatrixAdaptor &diffN_row = row_data.getDiffN(gg,nb_dofs_row/3);
	const MatrixAdaptor &diffN_col = col_data.getDiffN(gg,nb_dofs_col/3);
	//     printf("\nGauss Point %d at row:\n", gg);
	ierr = makeB(diffN_row,rowB); CHKERRQ(ierr);
	//printf("\nGauss Point %d at col:\n", gg);
	ierr = makeB(diffN_col,colB); CHKERRQ(ierr);
     
	
	noalias(CB) = prod(D,colB);
	noalias(K) += val*prod(trans(rowB),CB);
	// printf("\nGauss point %d:\n", gg);
	// for(unsigned int i = 0; i != D.size2(); ++i){
	//   for(unsigned int j = 0; j != D.size1(); ++j){
	//     printf("%e ", D(i,j) );
	//   }
	//   printf("\n");
	// }

      }
      PetscFunctionReturn(0);
    }
  
    /**
     * \brief Assemble local entity block matrix
     * @param  row_data row data (consist base functions on row entity)
     * @param  col_data column data (consist base functions on column entity)
     * @return          error code
     */
    virtual PetscErrorCode aSsemble(
      DataForcesAndSourcesCore::EntData &row_data,DataForcesAndSourcesCore::EntData &col_data
    ) {
      PetscFunctionBegin;
      // get pointer to first global index on row
      const int* row_indices = &*row_data.getIndices().data().begin();
      // get pointer to first global index on column
      const int* col_indices = &*col_data.getIndices().data().begin();
      Mat B = getFEMethod()->ksp_B!=PETSC_NULL? getFEMethod()->ksp_B : getFEMethod()->snes_B;
      // assemble local matrix
      ierr = MatSetValues(
        B, nbRows,row_indices,nbCols,col_indices,&*K.data().begin(),ADD_VALUES
      ); CHKERRQ(ierr);
      if(!isDiag&&sYmm) {
        // if not diagonal term and since global matrix is symmetric assemble
        // transpose term.
        K = trans(K);
        ierr = MatSetValues(
          B,nbCols,col_indices,nbRows,row_indices,&*K.data().begin(),ADD_VALUES
        ); CHKERRQ(ierr);
      }
      PetscFunctionReturn(0);
    }
    
  };



 struct ApplyDirichletBc: public MoFEM::FEMethod {

  Range fixFaces,fixNodes;

   ApplyDirichletBc(const Range& fix_faces,const Range& fix_nodes):
     MoFEM::FEMethod(),
    fixFaces(fix_faces),
    fixNodes(fix_nodes) {
  }
  
  PetscErrorCode postProcess() {
    
    PetscFunctionBegin;
    std::set<int> set_fix_dofs;

    cerr << fixFaces << endl;
    cerr << fixNodes << endl;

    //
    for(_IT_NUMEREDDOF_ROW_FOR_LOOP_(problemPtr,dit)) {
      if(dit->get()->getDofCoeffIdx()==2) {
	if(fixFaces.find(dit->get()->getEnt())!=fixFaces.end()) {
       set_fix_dofs.insert(dit->get()->getPetscGlobalDofIdx());
	}
      }
      
      if(fixNodes.find(dit->get()->getEnt())!=fixNodes.end()) {
	set_fix_dofs.insert(dit->get()->getPetscGlobalDofIdx());
      }
      
    }
    
    std::vector<int> fix_dofs(set_fix_dofs.size());
    
    std::copy(set_fix_dofs.begin(),set_fix_dofs.end(),fix_dofs.begin());
    
    ierr = MatAssemblyBegin(ksp_B,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
    ierr = MatAssemblyEnd(ksp_B,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
    ierr = VecAssemblyBegin(ksp_f); CHKERRQ(ierr);
    ierr = VecAssemblyEnd(ksp_f); CHKERRQ(ierr);
    

    Vec x;
    
    ierr = VecDuplicate(ksp_f,&x); CHKERRQ(ierr);
    ierr = VecZeroEntries(x); CHKERRQ(ierr);
    ierr = MatZeroRowsColumns(ksp_B,fix_dofs.size(),&fix_dofs[0],1,x,ksp_f); CHKERRQ(ierr);
    ierr = VecDestroy(&x); CHKERRQ(ierr);
    
    PetscFunctionReturn(0);
    
  }
  
};
 


struct OpPressure: MoFEM::FaceElementForcesAndSourcesCore::UserDataOperator {

double pressureVal;

OpPressure(const double pressure_val = 1):

MoFEM::FaceElementForcesAndSourcesCore::UserDataOperator("U",OPROW),
pressureVal(pressure_val) {
}
  //  virtual ~OpPressure();

VectorDouble nF;

FTensor::Index<'i',3> i;

PetscErrorCode doWork(int side, EntityType type, DataForcesAndSourcesCore::EntData &data) {

PetscFunctionBegin;

const int nb_dofs = data.getIndices().size();
if(nb_dofs==0) PetscFunctionReturn(0);

nF.resize(nb_dofs,false);
nF.clear();

const int nb_gauss_pts = data.getN().size1();
FTensor::Tensor1<double*,3> t_normal = getTensor1Normal();
FTensor::Tensor0<double*> t_base = data.getFTensor0N();

for(int gg = 0;gg!=nb_gauss_pts;gg++) {

double w = 0.5*getGaussPts()(2,gg);
FTensor::Tensor1<double*,3> t_nf(&nF[0],&nF[1],&nF[2],3);

for(int bb = 0;bb!=nb_dofs/3;bb++) {
t_nf(i) += (w*pressureVal*t_base)*t_normal(i);
++t_nf;
++t_base;
}

}

cerr << nF << endl;

ierr = VecSetValues(

getFEMethod()->ksp_f,nb_dofs,&data.getIndices()[0],&nF[0],ADD_VALUES

); CHKERRQ(ierr);

PetscFunctionReturn(0);

}

};

int main(int argc, char *argv[]) {

  // Initialize PETSCc
  PetscInitialize(&argc,&argv,(char *)0,help);

  try {

  // Create mesh databse
  moab::Core mb_instance;
  moab::Interface& moab = mb_instance;

// Get command line options
  int order = 3;  // default approximation order
  PetscBool flg_test = PETSC_FALSE; // true check if error is numerical error
  ierr = PetscOptionsBegin(PETSC_COMM_WORLD,"", "SimpleElasticProblem","none"); CHKERRQ(ierr);
  // Set approximation order
  ierr = PetscOptionsInt("-order","approximation order","",order,&order,PETSC_NULL); CHKERRQ(ierr);
  // Set testing (used by CTest)
  ierr = PetscOptionsBool("-test","if true is ctest","",flg_test,&flg_test,PETSC_NULL); CHKERRQ(ierr);
  ierr = PetscOptionsEnd(); CHKERRQ(ierr);


  // Create MoFEM databse and link it to MoAB
  MoFEM::Core core(moab);
  MoFEM::Interface& m_field = core;

  ierr = DMRegister_MoFEM("DMMOFEM"); CHKERRQ(ierr);

  Vec global_error;
  ierr = PoissonExample::AuxFunctions(m_field).createGhostVec(&global_error); CHKERRQ(ierr);
  
  boost::shared_ptr<ForcesAndSourcesCore> domain_lhs_fe;     ///< Volume element for the matrix
  boost::shared_ptr<ForcesAndSourcesCore> boundary_lhs_fe;   ///< Boundary element for the matrix
  boost::shared_ptr<ForcesAndSourcesCore> domain_rhs_fe;     ///< Volume element to assemble vector
  boost::shared_ptr<ForcesAndSourcesCore> boundary_rhs_fe;   ///< Volume element to assemble vector
  boost::shared_ptr<ForcesAndSourcesCore> domain_error;      ///< Volume element evaluate error
  boost::shared_ptr<ForcesAndSourcesCore> post_proc_volume;  ///< Volume element to Post-process results
  boost::shared_ptr<ForcesAndSourcesCore> null; 
 


  Simple *simple_interface;
  ierr = m_field.query_interface(simple_interface); CHKERRQ(ierr);

  ierr = simple_interface->getOptions(); CHKERRQ(ierr);
  ierr = simple_interface->loadFile(); CHKERRQ(ierr);


  ierr = simple_interface->addDomainField("U",H1,AINSWORTH_LEGENDRE_BASE,3); CHKERRQ(ierr);
  ierr = simple_interface->addDataField("ERROR",L2,AINSWORTH_LEGENDRE_BASE,3); CHKERRQ(ierr);
  ierr = simple_interface->setFieldOrder("U",order); CHKERRQ(ierr); // to approximate function
  ierr = simple_interface->setFieldOrder("ERROR",0); CHKERRQ(ierr); 


  Range fix_faces,pressure_faces,fix_nodes;


  for(_IT_CUBITMESHSETS_BY_BCDATA_TYPE_FOR_LOOP_(m_field,BLOCKSET,bit)) {
    EntityHandle meshset = bit->getMeshset();
    int id = bit->getMeshsetId();


    if(id == 1) {//brick-faces
      
      rval = m_field.get_moab().get_entities_by_dimension(meshset,2,fix_faces,true); CHKERRQ_MOAB(rval);

      Range adj_ents;
      rval = m_field.get_moab().get_adjacencies(fix_faces,0,false,adj_ents,moab::Interface::UNION); CHKERRQ_MOAB(rval);
      rval = m_field.get_moab().get_adjacencies(fix_faces,1,false,adj_ents,moab::Interface::UNION); CHKERRQ_MOAB(rval);
      fix_faces.merge(adj_ents);
      
    } else if(id == 2) {//node(s)

      rval = m_field.get_moab().get_entities_by_dimension(meshset,0,fix_nodes,true); CHKERRQ_MOAB(rval);
      
    } else if(id == 3) {//brick pressure faces

      rval = m_field.get_moab().get_entities_by_dimension(meshset,2,pressure_faces,true); CHKERRQ_MOAB(rval);
    
    } else {
      
      SETERRQ(PETSC_COMM_WORLD,MOFEM_DATA_INCONSISTENCY,"Unknown blockset");

    }

  }


ierr = simple_interface->defineFiniteElements(); CHKERRQ(ierr);

// Add pressure element

ierr = m_field.add_finite_element("PRESSURE"); CHKERRQ(ierr);
ierr = m_field.modify_finite_element_add_field_row("PRESSURE","U"); CHKERRQ(ierr);
ierr = m_field.modify_finite_element_add_field_col("PRESSURE","U"); CHKERRQ(ierr);
ierr = m_field.modify_finite_element_add_field_data("PRESSURE","U"); CHKERRQ(ierr);


ierr = simple_interface->defineProblem(); CHKERRQ(ierr);


 DM dm;
 ierr = simple_interface->getDM(&dm); CHKERRQ(ierr);
 ierr = DMMoFEMAddElement(dm,"PRESSURE"); CHKERRQ(ierr);
 ierr = DMMoFEMSetIsPartitioned(dm, PETSC_TRUE); CHKERRQ(ierr);

 ierr =simple_interface->buildFields(); CHKERRQ(ierr);
 ierr = simple_interface->buildFiniteElements(); CHKERRQ(ierr);
 
 ierr = m_field.add_ents_to_finite_element_by_dim(simple_interface->getMeshset(),simple_interface->getDim(),simple_interface->getDomainFEName(),true); CHKERRQ(ierr);

 ierr = m_field.build_finite_elements(simple_interface->getDomainFEName()); CHKERRQ(ierr);
 ierr = m_field.add_ents_to_finite_element_by_dim(pressure_faces,2,"PRESSURE"); CHKERRQ(ierr);
 ierr = m_field.build_finite_elements("PRESSURE",&pressure_faces); CHKERRQ(ierr);
 
 ierr = simple_interface->buildProblem(); CHKERRQ(ierr);


 boost::shared_ptr<FEMethod> nullFE;
 
 boost::shared_ptr<VolumeElementForcesAndSourcesCore> elastic_fe(new VolumeElementForcesAndSourcesCore(m_field));
 
 elastic_fe->getOpPtrVector().push_back(new OpK());
 
 // push operators to elastic_fe
 boost::shared_ptr<FaceElementForcesAndSourcesCore> pressure_fe(new FaceElementForcesAndSourcesCore(m_field));
 
 boost::shared_ptr<FEMethod> fix_dofs_fe(new ApplyDirichletBc(fix_faces, fix_nodes));
 
 pressure_fe->getOpPtrVector().push_back(new OpPressure());

 /*
  {
  post_proc_volume = boost::shared_ptr<ForcesAndSourcesCore>(new PostProcVolumeOnRefinedMesh(m_field));
      // Add operators to the elements, starting with some generic
  ierr = boost::static_pointer_cast<PostProcVolumeOnRefinedMesh>(post_proc_volume)->
    generateReferenceElementMesh(); CHKERRQ(ierr);
  ierr = boost::static_pointer_cast<PostProcVolumeOnRefinedMesh>(post_proc_volume)->
    addFieldValuesPostProc("U"); CHKERRQ(ierr);
  //ierr = boost::static_pointer_cast<PostProcVolumeOnRefinedMesh>(post_proc_volume)->
  // addFieldValuesPostProc("ERROR"); CHKERRQ(ierr);
  ierr = boost::static_pointer_cast<PostProcVolumeOnRefinedMesh>(post_proc_volume)->
    addFieldValuesGradientPostProc("U"); CHKERRQ(ierr);
  }  
 */
  
// Set operators for KSP solver
ierr = DMMoFEMKSPSetComputeOperators(
  dm,simple_interface->getDomainFEName(),elastic_fe,nullFE,nullFE
); CHKERRQ(ierr);
//ierr = DMMoFEMKSPSetComputeOperators(
//dm,simple_interface->getBoundaryFEName(), fix_dofs_fe,nullFE,nullFE
//); CHKERRQ(ierr);
// Set calculation of the right hand side vetor for KSP solver
ierr = DMMoFEMKSPSetComputeRHS(
  dm,"PRESSURE",pressure_fe,nullFE,nullFE
); CHKERRQ(ierr);
 

 Mat A;
 Vec x,f;
 
 ierr = DMCreateMatrix(dm,&A); CHKERRQ(ierr);
 ierr = DMCreateGlobalVector(dm,&x); CHKERRQ(ierr);
 ierr = VecDuplicate(x,&f); CHKERRQ(ierr);
 
 
 fix_dofs_fe->ksp_B = A;
 fix_dofs_fe->ksp_f = f;
 elastic_fe->ksp_B = A;
 pressure_fe->ksp_f = f;



// IGNATIOS: I do not know how the function with the nullFE objects can be called
//ierr = DMoFEMLoopFiniteElements(dm,simple_interface->getDomainFEName(),elastic_fe,nullFE,nullFE); CHKERRQ(ierr);
//ierr = DMoFEMLoopFiniteElements(dm,"PRESSURE",pressure_fe,nullFE,nullFE); CHKERRQ(ierr);
//ierr = DMoFEMLoopFiniteElements(dm,DM_NO_ELEMENT,nullFE,nullFE,fix_dofs_fe); CHKERRQ(ierr);

//IGNATIOS: invoking function without nulls
ierr = DMoFEMLoopFiniteElements(dm,simple_interface->getDomainFEName(),elastic_fe); CHKERRQ(ierr);
ierr = DMoFEMLoopFiniteElements(dm,"PRESSURE",pressure_fe); CHKERRQ(ierr);
//This is done because only post processor is implemented in the ApplyDirichletBc struct
 ierr = DMoFEMPostProcessFiniteElements(dm, fix_dofs_fe.get()); 
 //This would be problematic!
//ierr = DMoFEMLoopFiniteElements(dm,"DIRICHLETBC",fix_dofs_fe); CHKERRQ(ierr);

// MatView(A, PETSC_VIEWER_DRAW_WORLD);
// std::string wait;
// std::cin>>wait;
 MatView(A, PETSC_VIEWER_STDOUT_WORLD);

KSP solver;

ierr = KSPCreate(PETSC_COMM_WORLD,&solver); CHKERRQ(ierr);
ierr = KSPSetFromOptions(solver); CHKERRQ(ierr);
 ierr = KSPSetDM(solver,dm);
ierr = KSPSetOperators(solver,A,A); CHKERRQ(ierr);
ierr = KSPSetUp(solver); CHKERRQ(ierr);
ierr = KSPSolve(solver,f,x); CHKERRQ(ierr);

 VecView(x,PETSC_VIEWER_STDOUT_WORLD);


 // ierr = DMoFEMLoopFiniteElements(dm,simple_interface->getDomainFEName(),post_proc_volume); CHKERRQ(ierr);
 // Write results
 //ierr = boost::static_pointer_cast<PostProcVolumeOnRefinedMesh>(post_proc_volume)->writeFile("out_vol.h5m"); CHKERRQ(ierr);



ierr = MatDestroy(&A); CHKERRQ(ierr);
ierr = VecDestroy(&x); CHKERRQ(ierr);
ierr = VecDestroy(&f); CHKERRQ(ierr);
ierr = DMDestroy(&dm); CHKERRQ(ierr);

//This is a good reference for the future

/*
std::map<DofIdx,FieldData> mapZeroRows;

 const Problem *problemPtr;
 ierr = m_field.get_problem("SimpleProblem",&problemPtr); CHKERRQ(ierr);

 // for(_IT_CUBITMESHSETS_FOR_LOOP_(m_field,cit)){

 Range eNts;

 // get_root_set: is the mesh itself  
EntityHandle cubit_meshset = moab.get_root_set();
 rval = m_field.get_moab().get_entities_by_handle( cubit_meshset,eNts,true); CHKERRQ_MOAB(rval);

 for(std::vector<std::string>::iterator fit = simple_interface->getDomainFields()->begin();fit!=simple_interface->getDomainFields()->end();fit++) {
   for(Range::iterator eit = eNts.begin();eit!=eNts.end();++eit) {
  int counter = 0;//used to track dof directions
        for(_IT_NUMEREDDOF_ROW_BY_NAME_ENT_PART_FOR_LOOP_(problemPtr,*fit,*eit,m_field.get_comm_rank(),dof_ptr)) {
          NumeredDofEntity *dof = dof_ptr->get();
	  EntityHandle ent = (*dof).getEnt();
	  printf("\n %s \n", ent.getEntityType() );
	  double coords[3];
	  m_field.get_moab().get_coords(&ent,1,coords);
	  
	  if(fabs(coords[1]) <= TOL){
	    printf("type  of entity is  X coord %e and Z coord %e \n", coords[0], coords[2] );
	  }

	  ++counter;
	  if(fabs(coords[0]) <= TOL && 
	     fabs(coords[1]) <= TOL && 
	     fabs(coords[2]) <= TOL && 
	     (counter == 1 || counter == 3)   )  mapZeroRows[dof->getPetscGlobalDofIdx()] = 0;

	  if(counter == 2){
          mapZeroRows[dof->getPetscGlobalDofIdx()] = 0;
	  }

        }
   }
   }*/

  } catch (MoFEMException const &e) {
    SETERRQ(PETSC_COMM_SELF,e.errorCode,e.errorMessage);
  }

// finish work cleaning memory, getting statistics, etc
ierr =  PetscFinalize(); CHKERRQ(ierr);
  
  return 0;
}
