/* Copyright (C) 2013, Lukasz Kaczmarczyk (likask AT wp.pl)
 * --------------------------------------------------------------
 * FIXME: DESCRIPTION
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

#include <DirichletBC.hpp>
#include <Projection10NodeCoordsOnField.hpp>

#include <petsctime.h>

#include <SurfacePressure.hpp>
#include <FEMethod_LowLevelStudent.hpp>
#include <FEMethod_UpLevelStudent.hpp>

#include <PostProcVertexMethod.hpp>
#include <PostProcDisplacementAndStrainOnRefindedMesh.hpp>

#include <ElasticFEMethod.hpp>
#include "ElasticFEMethodTransIso.hpp"

#include "ElasticFE_RVELagrange_Disp.hpp"
#include "ElasticFE_RVELagrange_Disp_Multi_Rhs.hpp"
#include "ElasticFE_RVELagrange_Homogenized_Stress_Disp.hpp"
#include "RVEVolume.hpp"

using namespace ObosleteUsersModules;
#include "MaterialConstitutiveMatrix_FirstOrderDerivative.hpp"
#include "MaterialConstitutiveMatrix_SecondOrderDerivative.hpp"
#include "Trans_Iso_Rhs_r_PSFEM.hpp"
#include "Trans_Iso_Rhs_rs_PSFEM.hpp"

#include <FE2_ElasticFEMethod.hpp>

#include <FE2_Rhs_r_PSFEM.hpp>
#include <FE2_Rhs_rs_PSFEM.hpp>

#include <Reliability_SurfacePressure.hpp>

#include <FE2_Macro_Solver.hpp>

#include <FE2_PostProcStressForReliability.hpp>

using namespace boost::numeric;

//======================================================
// Declaration for reliability analysis
extern "C" {
#include <gm_rule.h>
#include <ltqnorm.h>
}

#include <iostream>
#include <fstream>
#include <ctime>
//#include <vector>
#include <new>
#include <ctype.h>

#include <boost/random.hpp>
#include <boost/math/distributions/normal.hpp>
#include <boost/math/distributions/lognormal.hpp>
#include <boost/math/distributions/exponential.hpp>
#include <boost/math/distributions/extreme_value.hpp>
#include <boost/math/distributions/weibull.hpp>
#include <boost/math/distributions/uniform.hpp>
#include <boost/math/distributions/gamma.hpp>

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/triangular.hpp>
#include <boost/numeric/ublas/io.hpp>

#include <cholesky.hpp>
#include <MatrixInverse.hpp> // download from http://proteowizard.sourceforge.net/dox/_matrix_inverse_8hpp.html

#include <ImportProbData.hpp>
#include <NatafTransformation.hpp>
#include <LimitStateFunction.hpp>

//======================================================

/*******************************************************************************
 *                                                                             *
 *             FUNCTIONS FOR STRUCTURAL RELIABILITY ANLYSIS                    *
 *                                                                             *
/*******************************************************************************/

//------------------------------------------------------------------------------
// Construct data structure <Stochastic_Model> for collecting data representing
//   statistical information of inputs including probability distribution,
//   correlation matrix and etc.

struct Stochastic_Model {
  int    num_vars;                       // Number of variables
  double dist_type;                      // Distribution type index
  double transf_type;                    // Type of joint distribution
  double R0_method;                      // Method for computation of the modified Nataf correlation matrix
  int flag_sens;                         // Flag for computation of sensitivities w.r.t. parameters
  vector<string> NameVars;               // Name of random variables
  ublas::matrix<double> correlation;     // Correlation matrix
  ublas::matrix<double> marg;            // Marginal distribution for each random variable
  ublas::matrix<double> mod_correlation; // modified correlation matrix
  ublas::matrix<double> Lo;              // Chelosky decomposition
  ublas::matrix<double> inv_Lo;          // inverse of matrix Lo
  ublas::matrix<double> MatStrength;     // Material strength
};

//------------------------------------------------------------------------------
// Construct data structure <Reliability_Options> to define calculation options
//

struct Reliability_Options {
  int echo_flag;            // Program interactive mode, 0: silent mode
  // FORM analysis options
  int istep_max;            // Maximum number of interations allowed in the search algorithm
  double e1;                // Tolerance on how close design point is to limit-state surface
  double e2;                // Tolerance on how accurately the gradient points towards the origin
  double step_code;         // 0: step size by Armijo rule, otherwise: given value is the step size
  int Recorded_u;           // 0: u-vector not recorded at all iterations, 1: u-vector recorded at all iterations
  int Recorded_x;           // 0: x-vector not recorded at all iterations, 1: x-vector recorded at all iterations
  int Recorded_beta;        // 0: beta not recorded at all iterations, 1: recorded at all iterations
  string grad_G;            // "PSFEM": perturbation based SFEM, "DDM": direct differentiation, 'ADM': automatic differentiation
};

//------------------------------------------------------------------------------
// Construct data structure <LSF_Options> to define limit-state function options
//

struct LSF_Options {
  string evaluator; // Type of limit-state function evaluator,
                    //   "basic", the LSF is defined by means of an analytical expression
                    //   "failure", the LSF is defined by using failure criterion
  int flag_sens;    // Flag for sensitivity computation w.r.t. specified parameters of the LSF, 0: no, 1: yes
};


//------------------------------------------------------------------------------
// To determine search direction
//

void search_dir(double val_G,
                ublas::vector<double> grad_G,
                ublas::vector<double> u,
                ublas::vector<double> &u_dir) {
  // Determin direction cosine vector
  boost::numeric::ublas::vector<double> alpha;
  alpha = -grad_G/norm_2(grad_G);
  
  // Compute direction
  u_dir = ((val_G/norm_2(grad_G)) + inner_prod(alpha, u))*alpha - u;
}


ErrorCode rval;
PetscErrorCode ierr;

static char help[] = "...\n\n";

const char* args[] = {
  "_r_Em", "_r_NUm",
  "_r_NUp",     "_r_NUpz",      "_r_Ep",    "_r_Ez",    "_r_Gzp",     // 1st order
  "_rs_EmEm", "_rs_NUmNUm",                                             // 2nd order
  "_rs_NUpNUp", "_rs_NUpzNUpz", "_rs_EpEp", "_rs_EzEz", "_rs_GzpGzp",
};

int nvars = 7;    // number of variables
int nders = 14;   // number of partial derivatives (firsr- and second- order)
vector<string> stochastic_fields(args, args + 14);


int main(int argc, char *argv[]) {
  
  clock_t start_time, finish_time;
  double total_time;
  start_time = clock();
  
  
  //============================================================================
  //
  //  A. Micro (RVE) Problem
  //
  //============================================================================
  
  PetscInitialize(&argc,&argv,(char *)0,help);
  
  moab::Core mb_instance_RVE;
  Interface& moab_RVE = mb_instance_RVE;
  int rank;
  MPI_Comm_rank(PETSC_COMM_WORLD,&rank);
  
  // ===========================================================================
  //
  //  A.I. READ MESH DATA AND FEA COMPUTATION PARAMETERS FROM FILE
  //
  // ===========================================================================
  
  /*****************************************************************************
   *
   * Read parameters from line command
   *
   ****************************************************************************/
  
  PetscBool flg = PETSC_TRUE;
  char mesh_file_name[255];
  ierr = PetscOptionsGetString(PETSC_NULL,"-my_file_RVE",mesh_file_name,255,&flg); CHKERRQ(ierr);
  if(flg != PETSC_TRUE) {
    SETERRQ(PETSC_COMM_SELF,1,"*** ERROR -my_file_RVE (MESH FILE NEEDED)");
  }
  PetscInt order;
  ierr = PetscOptionsGetInt(PETSC_NULL,"-my_order",&order,&flg); CHKERRQ(ierr);
  if(flg != PETSC_TRUE) {
    order = 1;
  }
  
 // char outName[PETSC_MAX_PATH_LEN]="out.vtk";
 // ierr = PetscOptionsGetString(PETSC_NULL,"-my_out",outName,sizeof(outName),&flg); CHKERRQ(ierr);
  
  //char outName2[PETSC_MAX_PATH_LEN]="out_post_proc.vtk";
  //ierr = PetscOptionsGetString(PETSC_NULL,"-my_post_out",outName2,sizeof(outName2),&flg); CHKERRQ(ierr);

  /*****************************************************************************
   *
   * Transfer mesh data to MOAB database
   *
   ****************************************************************************/
  
  const char *option;
  option = "";//"PARALLEL=BCAST;";//;DEBUG_IO";
  rval = moab_RVE.load_file(mesh_file_name, 0, option); CHKERR_PETSC(rval);
  ParallelComm* pcomm = ParallelComm::get_pcomm(&moab_RVE,MYPCOMM_INDEX);
  if(pcomm == NULL) pcomm =  new ParallelComm(&moab_RVE,PETSC_COMM_WORLD);
  
  //Create MoFEM (Joseph) database
  MoFEM::Core core_RVE(moab_RVE);
  FieldInterface& m_field_RVE = core_RVE;
  
  
  /*****************************************************************************
   *
   * Get fibre direction information from the potential-flow calculation
   *
   ****************************************************************************/
  
  Tag th_phi;
  //    double def_val  = 0;
  rval = moab_RVE.tag_get_handle("PHI",1,MB_TYPE_DOUBLE,th_phi); CHKERR_PETSC(rval);
	
  Tag th_meshset_info;
  int def_meshset_info[2] = {0,0};
  rval = moab_RVE.tag_get_handle("MESHSET_INFO",2,MB_TYPE_INTEGER,th_meshset_info,MB_TAG_CREAT|MB_TAG_SPARSE,&def_meshset_info);
  
  int meshset_data[2];
  EntityHandle root = moab_RVE.get_root_set();
  rval = moab_RVE.tag_get_data(th_meshset_info,&root,1,meshset_data); CHKERR_PETSC(rval);
  
  vector<BitRefLevel> bit_levels;
  bit_levels.push_back(BitRefLevel().set(meshset_data[0]-1));
  
  ierr = m_field_RVE.build_fields(); CHKERRQ(ierr);
  ierr = m_field_RVE.build_finite_elements(); CHKERRQ(ierr);
  ierr = m_field_RVE.build_adjacencies(bit_levels.back()); CHKERRQ(ierr);
  ierr = m_field_RVE.build_problems(); CHKERRQ(ierr);

  /*****************************************************************************
   *
   * Group element into various mesh-set
   * meshset_level0: all element
   * meshset_Matrix: element set for matrix
   * meshset_Inclusion: element set for inclusion/fibre
   *
   ****************************************************************************/
  
  EntityHandle out_meshset;
  rval = moab_RVE.create_meshset(MESHSET_SET,out_meshset); CHKERR_PETSC(rval);
  //    ierr = m_field_RVE.get_problem_finite_elements_entities("POTENTIAL_PROBLEM","POTENTIAL_ELEM",out_meshset); CHKERRQ(ierr);
  ierr = m_field_RVE.get_entities_by_ref_level(bit_levels.back(),BitRefLevel().set(),out_meshset); CHKERRQ(ierr);
  Range LatestRefinedTets;
  rval = moab_RVE.get_entities_by_type(out_meshset, MBTET,LatestRefinedTets,true); CHKERR_PETSC(rval);
  
  Range LatestRefinedPrisms;
  rval = moab_RVE.get_entities_by_type(out_meshset, MBPRISM,LatestRefinedPrisms,true); CHKERR_PETSC(rval);
	
  cout<<"No of Prisms/Interfaces = "<<LatestRefinedPrisms.size()<<endl;
	
  BitRefLevel problem_bit_level_RVE = bit_levels.back();
  
  EntityHandle meshset_Matrix, meshset_Fibre;
  rval = moab_RVE.create_meshset(MESHSET_SET,meshset_Matrix); CHKERR_PETSC(rval);
  rval = moab_RVE.create_meshset(MESHSET_SET,meshset_Fibre); CHKERR_PETSC(rval);
  
  // Select fibre reinforcement related elements into meshset of <meshset_Fibre>
  ///Getting No. of Fibres to be used for Potential Flow Problem
  int noOfFibres=0;
  for(_IT_CUBITMESHSETS_BY_BCDATA_TYPE_FOR_LOOP_(m_field_RVE,BLOCKSET|UNKNOWNCUBITNAME,it)) {
    
    std::size_t found=it->get_name().find("PotentialFlow");
    if (found==std::string::npos) continue;
    noOfFibres += 1;
  }
  cout<<"No. of Fibres for Potential Flow : "<<noOfFibres<<endl;

//  vector<int> fibreList(noOfFibres,0);
//  for (int aa=0; aa<noOfFibres; aa++) {
//    fibreList[aa] = aa + 1;
//  }
  
  vector<Range> RangeFibre(noOfFibres);
  vector<EntityHandle> fibre_meshset(noOfFibres);
  
  for (int ii=0; ii<noOfFibres; ii++) {
    ostringstream sss;
    sss << "POTENTIAL_ELEM" << ii+1;
    for(_IT_GET_FES_BY_NAME_FOR_LOOP_(m_field_RVE, sss.str().c_str() ,it)){
      RangeFibre[ii].insert(it->get_ent());
      rval = moab_RVE.create_meshset(MESHSET_SET,fibre_meshset[ii]); CHKERR_PETSC(rval);
      rval = moab_RVE.add_entities(fibre_meshset[ii],RangeFibre[ii]); CHKERR_PETSC(rval);
      rval = moab_RVE.unite_meshset(meshset_Fibre,fibre_meshset[ii]); CHKERR_PETSC(rval);
    }
  }
  rval = moab_RVE.write_file("meshset_Fibre.vtk","VTK","",&meshset_Fibre,1); CHKERR_PETSC(rval);
  
  // Select matrix related elements into meshset of <meshset_Matrix>
  for(_IT_CUBITMESHSETS_BY_SET_TYPE_FOR_LOOP_(m_field_RVE,BLOCKSET,it)) {
    
    if(it->get_name() == "MAT_ELASTIC_1") {
      Range TetsInBlock;
      rval = moab_RVE.get_entities_by_type(it->meshset, MBTET,TetsInBlock,true); CHKERR_PETSC(rval);
      Range block_rope_bit_level = intersect(LatestRefinedTets,TetsInBlock);
      
      cout<<"=============  TetsInBlock  "<< TetsInBlock.size() <<endl;
      
      rval = moab_RVE.add_entities(meshset_Matrix,block_rope_bit_level);CHKERR_PETSC(rval);
      
    }
  }

  ierr = m_field_RVE.seed_finite_elements(meshset_Matrix); CHKERRQ(ierr);
  
  Range prims_on_problem_bit_level_RVE;
  ierr = m_field_RVE.get_entities_by_type_and_ref_level(problem_bit_level_RVE,BitRefLevel().set(),MBPRISM,prims_on_problem_bit_level_RVE); CHKERRQ(ierr);
  //to create meshset from range
  EntityHandle meshset_prims_on_problem_bit_level_RVE;
  rval = moab_RVE.create_meshset(MESHSET_SET,meshset_prims_on_problem_bit_level_RVE); CHKERR_PETSC(rval);
  rval = moab_RVE.add_entities(meshset_prims_on_problem_bit_level_RVE,prims_on_problem_bit_level_RVE); CHKERR_PETSC(rval);
  ierr = m_field_RVE.seed_ref_level_MESHSET(meshset_prims_on_problem_bit_level_RVE,BitRefLevel().set()); CHKERRQ(ierr);

  
  // ===========================================================================
  //
  // A.II. DEFINE PROBLEM
  //
  // ===========================================================================
  
  /*****************************************************************************
   *
   * Add field
   *  (1) Deterministic fields
   *  (2) Stochastic fields
   *       (total 14 field for 1st and 2nd order stochastic PSFEM)
   *
   ****************************************************************************/
  // Deterministic fields
  int field_rank = 3;
  ierr = m_field_RVE.add_field("DISP_RVE",H1,field_rank); CHKERRQ(ierr);
  ierr = m_field_RVE.add_field("Lagrange_mul_disp",H1,field_rank); CHKERRQ(ierr);
  
  // Stochastic fields for perturbation method
  for(int ii=0; ii < nders; ii++ )
  {
    ostringstream ss_field;
    ss_field << "DISP_RVE" << stochastic_fields[ii];
    cout<<ss_field.str().c_str()<<endl;
    ierr = m_field_RVE.add_field(ss_field.str().c_str(),H1,field_rank,MF_ZERO); CHKERRQ(ierr);
  }
  
  
  /*****************************************************************************
   *
   * Create finite element for the defined fields
   *
   ****************************************************************************/
  ierr = m_field_RVE.add_finite_element("ELASTIC_FE_RVE"); CHKERRQ(ierr);
  ierr = m_field_RVE.add_finite_element("TRAN_ISO_FE_RVE"); CHKERRQ(ierr);
  ierr = m_field_RVE.add_finite_element("Lagrange_FE"); CHKERRQ(ierr);
  
  
  /*****************************************************************************
   *
   * set field data which finite element use
   * set field row which finite element use
   *
   * [K][U] = [F]                                           Zeroth-order problem
   * [K][D_r U] = - [D_r K][U]                               First-order problem
   * [K][H_rs U] = - [H_rs K][U] - 2[D_r K][D_s U]          Second-order problem
   *
   ****************************************************************************/
  //Define rows/cols and element data
  ierr = m_field_RVE.modify_finite_element_add_field_row("ELASTIC_FE_RVE","DISP_RVE"); CHKERRQ(ierr);
  ierr = m_field_RVE.modify_finite_element_add_field_col("ELASTIC_FE_RVE","DISP_RVE"); CHKERRQ(ierr);
  ierr = m_field_RVE.modify_finite_element_add_field_data("ELASTIC_FE_RVE","DISP_RVE"); CHKERRQ(ierr);
  
  // Stochastic
  for(int ii=0; ii < nvars; ii++ ) {
    ostringstream ss_field;
    ss_field << "DISP_RVE" << stochastic_fields[ii];
    cout<<ss_field.str().c_str()<<endl;
    ierr = m_field_RVE.modify_finite_element_add_field_data("ELASTIC_FE_RVE",ss_field.str().c_str()); CHKERRQ(ierr);
  }
  
  //FE Transverse Isotropic
  ierr = m_field_RVE.modify_finite_element_add_field_row("TRAN_ISO_FE_RVE","DISP_RVE"); CHKERRQ(ierr);
  ierr = m_field_RVE.modify_finite_element_add_field_col("TRAN_ISO_FE_RVE","DISP_RVE"); CHKERRQ(ierr);
  ierr = m_field_RVE.modify_finite_element_add_field_data("TRAN_ISO_FE_RVE","DISP_RVE"); CHKERRQ(ierr);
  ierr = m_field_RVE.modify_finite_element_add_field_data("TRAN_ISO_FE_RVE","POTENTIAL_FIELD"); CHKERRQ(ierr);
  
  for(int ii=0; ii < nvars; ii++ ) {
    ostringstream ss_field;
    ss_field << "DISP_RVE" << stochastic_fields[ii];
    cout<<ss_field.str().c_str()<<endl;
    ierr = m_field_RVE.modify_finite_element_add_field_data("TRAN_ISO_FE_RVE",ss_field.str().c_str()); CHKERRQ(ierr);
  }
  
  //C and CT
  //======================================================================================================
  //C row as Lagrange_mul_disp and col as DISPLACEMENT
  ierr = m_field_RVE.modify_finite_element_add_field_row("Lagrange_FE","Lagrange_mul_disp"); CHKERRQ(ierr);
  ierr = m_field_RVE.modify_finite_element_add_field_col("Lagrange_FE","DISP_RVE"); CHKERRQ(ierr);
  
  //CT col as Lagrange_mul_disp and row as DISPLACEMENT
  ierr = m_field_RVE.modify_finite_element_add_field_col("Lagrange_FE","Lagrange_mul_disp"); CHKERRQ(ierr);
  ierr = m_field_RVE.modify_finite_element_add_field_row("Lagrange_FE","DISP_RVE"); CHKERRQ(ierr);
  
  //data
  ierr = m_field_RVE.modify_finite_element_add_field_data("Lagrange_FE","Lagrange_mul_disp"); CHKERRQ(ierr);
  ierr = m_field_RVE.modify_finite_element_add_field_data("Lagrange_FE","DISP_RVE"); CHKERRQ(ierr);
  //======================================================================================================
  
  
  //define problems
  ierr = m_field_RVE.add_problem("ELASTIC_PROBLEM_RVE"); CHKERRQ(ierr);
  
  
  //set finite elements for problem
  ierr = m_field_RVE.modify_problem_add_finite_element("ELASTIC_PROBLEM_RVE","ELASTIC_FE_RVE"); CHKERRQ(ierr);
  ierr = m_field_RVE.modify_problem_add_finite_element("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE"); CHKERRQ(ierr);
  ierr = m_field_RVE.modify_problem_add_finite_element("ELASTIC_PROBLEM_RVE","Lagrange_FE"); CHKERRQ(ierr);
  
  
  //set refinment level for problem
  ierr = m_field_RVE.modify_problem_ref_level_add_bit("ELASTIC_PROBLEM_RVE",problem_bit_level_RVE); CHKERRQ(ierr); // problem_bit_level
  
  // ===========================================================================
  //
  // A.III. DECLARE PROBLEM
  //
  // ===========================================================================
  
  
  /*****************************************************************************
   *
   * Add entitities (by tets) to the field
   *
   ****************************************************************************/
  ierr = m_field_RVE.add_ents_to_field_by_TETs(0,"DISP_RVE"); CHKERRQ(ierr);
  for(int ii=0; ii < nders; ii++ ) {
    ostringstream ss_field;
    ss_field << "DISP_RVE" << stochastic_fields[ii];
    ierr = m_field_RVE.add_ents_to_field_by_TETs(0,ss_field.str().c_str()); CHKERRQ(ierr);
  }
  
  
  /*****************************************************************************
   *
   * Add finite elements entities
   *
   ****************************************************************************/
  ierr = m_field_RVE.add_ents_to_finite_element_by_TETs(meshset_Matrix,"ELASTIC_FE_RVE",true); CHKERRQ(ierr);
  ierr = m_field_RVE.add_ents_to_finite_element_by_TETs(meshset_Fibre,"TRAN_ISO_FE_RVE",true); CHKERRQ(ierr);

  Range SurfacesFaces;
  ierr = m_field_RVE.get_cubit_msId_entities_by_dimension(103,SIDESET,2,SurfacesFaces,true); CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"number of SideSet 103 = %d\n",SurfacesFaces.size()); CHKERRQ(ierr);
  ierr = m_field_RVE.add_ents_to_finite_element_by_TRIs(SurfacesFaces,"Lagrange_FE"); CHKERRQ(ierr);
  
  
  //to create meshset from range
  EntityHandle BoundFacesMeshset;
  rval = moab_RVE.create_meshset(MESHSET_SET,BoundFacesMeshset); CHKERR_PETSC(rval);
  rval = moab_RVE.add_entities(BoundFacesMeshset,SurfacesFaces); CHKERR_PETSC(rval);
  ierr = m_field_RVE.seed_ref_level_MESHSET(BoundFacesMeshset,BitRefLevel().set()); CHKERRQ(ierr);
  ierr = m_field_RVE.add_ents_to_field_by_TRIs(BoundFacesMeshset,"Lagrange_mul_disp",2); CHKERRQ(ierr);
  
  /*****************************************************************************
   *
   * Set applied order
   * See reference for detals:
   *   Ainsworth M. and Coyle J. (2003) Hierarchic finite element bases on
   unstructured tetrahedral meshes. IJNME, 58(14). pp.2103-2130.
   ****************************************************************************/
  ierr = m_field_RVE.set_field_order(0,MBTET,"DISP_RVE",order); CHKERRQ(ierr);
  ierr = m_field_RVE.set_field_order(0,MBTRI,"DISP_RVE",order); CHKERRQ(ierr);
  ierr = m_field_RVE.set_field_order(0,MBEDGE,"DISP_RVE",order); CHKERRQ(ierr);
  ierr = m_field_RVE.set_field_order(0,MBVERTEX,"DISP_RVE",1); CHKERRQ(ierr);
  
  ierr = m_field_RVE.set_field_order(0,MBTRI,"Lagrange_mul_disp",order); CHKERRQ(ierr);
  ierr = m_field_RVE.set_field_order(0,MBEDGE,"Lagrange_mul_disp",order); CHKERRQ(ierr);
  ierr = m_field_RVE.set_field_order(0,MBVERTEX,"Lagrange_mul_disp",1); CHKERRQ(ierr);
  
  int order_st=order;
  for(int ii=0; ii < nders; ii++ ) {
    ostringstream ss_field;
    ss_field << "DISP_RVE" << stochastic_fields[ii];
    //    cout<<ss_field.str().c_str()<<endl;
    ierr = m_field_RVE.set_field_order(0,MBTET,ss_field.str().c_str(),order_st); CHKERRQ(ierr);
    ierr = m_field_RVE.set_field_order(0,MBTRI,ss_field.str().c_str(),order_st); CHKERRQ(ierr);
    ierr = m_field_RVE.set_field_order(0,MBEDGE,ss_field.str().c_str(),order_st); CHKERRQ(ierr);
    ierr = m_field_RVE.set_field_order(0,MBVERTEX,ss_field.str().c_str(),1); CHKERRQ(ierr);
  }
  
  
  
  // ===========================================================================
  //
  //  A.IV. BUILD DATABASE
  //
  // ===========================================================================
  
  //build field
  ierr = m_field_RVE.build_fields(); CHKERRQ(ierr);
  
  //build finite elemnts
  ierr = m_field_RVE.build_finite_elements(); CHKERRQ(ierr);
  
  //build adjacencies
  ierr = m_field_RVE.build_adjacencies(problem_bit_level_RVE); CHKERRQ(ierr);
  
  //build problem
  ierr = m_field_RVE.build_problems(); CHKERRQ(ierr);
  
  
  
  // ===========================================================================
  //
  //  A.V. MESH PARTITION
  //
  // ===========================================================================
  
  //partition
  ierr = m_field_RVE.partition_problem("ELASTIC_PROBLEM_RVE"); CHKERRQ(ierr);
  ierr = m_field_RVE.partition_finite_elements("ELASTIC_PROBLEM_RVE"); CHKERRQ(ierr);
  //what are ghost nodes, see Petsc Manual
  ierr = m_field_RVE.partition_ghost_dofs("ELASTIC_PROBLEM_RVE"); CHKERRQ(ierr);
  
  //print bcs
  ierr = m_field_RVE.print_cubit_displacement_set(); CHKERRQ(ierr);
  ierr = m_field_RVE.print_cubit_force_set(); CHKERRQ(ierr);
  //print block sets with materials
  ierr = m_field_RVE.print_cubit_materials_set(); CHKERRQ(ierr);
  
  
  
  //============================================================================
  //
  //  B. Macro Problem
  //
  //============================================================================
  
  // ===========================================================================
  //
  // B. I. READ MESH DATA AND FINITE ELEMENT ANALYSIS CONTROL PARAMETERS FROM FILE
  //
  // ===========================================================================
  moab::Core mb_instance_Macro;
  Interface& moab_Macro = mb_instance_Macro;
  ParallelComm* pcomm_Macro = ParallelComm::get_pcomm(&moab_Macro,MYPCOMM_INDEX);
  if(pcomm == NULL) pcomm_Macro =  new ParallelComm(&moab_Macro,PETSC_COMM_WORLD);
  
  /*****************************************************************************
   *
   * Read parameters from line command
   *
   ****************************************************************************/
  ierr = PetscOptionsGetString(PETSC_NULL,"-my_file_macro",mesh_file_name,255,&flg); CHKERRQ(ierr);
  if(flg != PETSC_TRUE) {
    SETERRQ(PETSC_COMM_SELF,1,"*** ERROR -my_file_macro (MESH FILE NEEDED)");
  }
  
  /*****************************************************************************
   *
   * Transfer mesh data to MOAB database
   *
   ****************************************************************************/
  //option = "PARALLEL=BCAST_DELETE;"
  //"PARTITION=GEOM_DIMENSION,PARTITION_VAL=3,PARTITION_DISTRIBUTE";//;DEBUG_IO";
  rval = moab_Macro.load_file(mesh_file_name, 0, option); CHKERR_PETSC(rval);
  
  //Create MoFEM (Joseph) database
  MoFEM::Core core_Macro(moab_Macro);
  FieldInterface& m_field_Macro = core_Macro;
  
  //ref meshset ref level 0
  ierr = m_field_Macro.seed_ref_level_3D(0,0); CHKERRQ(ierr);
  
  // stl::bitset see for more details
  BitRefLevel bit_level0_Macro;
  bit_level0_Macro.set(0);
  
  //EntityHandle meshset_level0_Macro;
  //rval = moab_Macro.create_meshset(MESHSET_SET,meshset_level0_Macro); CHKERR_PETSC(rval);
  ierr = m_field_Macro.seed_ref_level_3D(0,bit_level0_Macro); CHKERRQ(ierr);
  //ierr = m_field_Macro.get_entities_by_ref_level(bit_level0_Macro,BitRefLevel().set(),meshset_level0_Macro); CHKERRQ(ierr);
  
  EntityHandle meshset_macro_Elastic,meshset_Reliability;
  rval = moab_Macro.create_meshset(MESHSET_SET,meshset_macro_Elastic); CHKERR_PETSC(rval);
  rval = moab_Macro.create_meshset(MESHSET_SET,meshset_Reliability); CHKERR_PETSC(rval);
  // select elements into group
  for(_IT_CUBITMESHSETS_BY_SET_TYPE_FOR_LOOP_(m_field_Macro,BLOCKSET,it)){
    
    if(it->get_name() == "MAT_ELASTIC") {
      Range TetsInBlock;
      rval = moab_Macro.get_entities_by_type(it->meshset, MBTET,TetsInBlock,true); CHKERR_PETSC(rval);
      Range block_rope_bit_level = intersect(LatestRefinedTets,TetsInBlock);
      
      cout<<"=============  TetsInBlock  "<< TetsInBlock.size() <<endl;
      
      rval = moab_Macro.add_entities(meshset_macro_Elastic,block_rope_bit_level);CHKERR_PETSC(rval);
      
    }
  }
  ierr = m_field_Macro.seed_finite_elements(meshset_macro_Elastic); CHKERRQ(ierr);
  
  // select reliability calculation related elements into element-set
  for(_IT_CUBITMESHSETS_BY_SET_TYPE_FOR_LOOP_(m_field_Macro,BLOCKSET,it)){
    
    if(it->get_name() == "RELIABILITY") {
      Range TetsInBlock;
      rval = moab_Macro.get_entities_by_type(it->meshset, MBTET,TetsInBlock,true); CHKERR_PETSC(rval);
      Range block_rope_bit_level = intersect(LatestRefinedTets,TetsInBlock);
      
      cout<<"=============  TetsInBlock  "<< TetsInBlock.size() <<endl;
      
      rval = moab_Macro.add_entities(meshset_Reliability,block_rope_bit_level);CHKERR_PETSC(rval);
      
    }
  }
  ierr = m_field_Macro.seed_finite_elements(meshset_Reliability); CHKERRQ(ierr);
  
  
  
  // ===========================================================================
  //
  // B.II. DEFINE PROBLEM
  //
  // ===========================================================================
  

  /*****************************************************************************
   *
   * Add stochastic field
   * (total 14 field for 1st and 2nd order stochastic PSFEM)
   *
   ****************************************************************************/
  ierr = m_field_Macro.add_field("DISP_MACRO",H1,3,MF_ZERO); CHKERRQ(ierr);
  ierr = m_field_Macro.add_field("MESH_NODE_POSITIONS",H1,3,MF_ZERO); CHKERRQ(ierr);
  
  // Stochastic fields for perturbation methods at macroscale
  for(int ii=0; ii < nders; ii++ ) {
    ostringstream ss_field;
    ss_field << "DISP_MACRO" << stochastic_fields[ii];
    //cout<<ss_field.str().c_str()<<endl;
    ierr = m_field_Macro.add_field(ss_field.str().c_str(),H1,field_rank,MF_ZERO); CHKERRQ(ierr);
  }

  
  /*****************************************************************************
   *
   * Create finite element for the defined fields
   *
   ****************************************************************************/
  //ierr = m_field_Macro.add_finite_element("ELASTIC_FE_MACRO",MF_ZERO); CHKERRQ(ierr);
  ierr = m_field_Macro.add_finite_element("ELASTIC_FE_MACRO"); CHKERRQ(ierr);
  ierr = m_field_Macro.add_finite_element("ELASTIC_FE_MACRO_REL"); CHKERRQ(ierr);
 

  /*****************************************************************************
   *
   * set field data which finite element use
   *
   ****************************************************************************/
  //Define rows/cols and element data
  ierr = m_field_Macro.modify_finite_element_add_field_row("ELASTIC_FE_MACRO","DISP_MACRO"); CHKERRQ(ierr);
  ierr = m_field_Macro.modify_finite_element_add_field_col("ELASTIC_FE_MACRO","DISP_MACRO"); CHKERRQ(ierr);
  ierr = m_field_Macro.modify_finite_element_add_field_data("ELASTIC_FE_MACRO","DISP_MACRO"); CHKERRQ(ierr);
  ierr = m_field_Macro.modify_finite_element_add_field_data("ELASTIC_FE_MACRO","MESH_NODE_POSITIONS"); CHKERRQ(ierr);
  
  for(int ii=0; ii < nvars; ii++ ) {
    ostringstream ss_field;
    ss_field << "DISP_MACRO" << stochastic_fields[ii];
    //    cout<<ss_field.str().c_str()<<endl;
    ierr = m_field_Macro.modify_finite_element_add_field_data("ELASTIC_FE_MACRO",ss_field.str().c_str()); CHKERRQ(ierr);
  }
  
  
  //Define rows/cols and element data
  ierr = m_field_Macro.modify_finite_element_add_field_row("ELASTIC_FE_MACRO_REL","DISP_MACRO"); CHKERRQ(ierr);
  ierr = m_field_Macro.modify_finite_element_add_field_col("ELASTIC_FE_MACRO_REL","DISP_MACRO"); CHKERRQ(ierr);
  ierr = m_field_Macro.modify_finite_element_add_field_data("ELASTIC_FE_MACRO_REL","DISP_MACRO"); CHKERRQ(ierr);
  ierr = m_field_Macro.modify_finite_element_add_field_data("ELASTIC_FE_MACRO_REL","MESH_NODE_POSITIONS"); CHKERRQ(ierr);
  
  for(int ii=0; ii < nvars; ii++ ) {
    ostringstream ss_field;
    ss_field << "DISP_MACRO" << stochastic_fields[ii];
    //    cout<<ss_field.str().c_str()<<endl;
    ierr = m_field_Macro.modify_finite_element_add_field_data("ELASTIC_FE_MACRO_REL",ss_field.str().c_str()); CHKERRQ(ierr);
  }
  
  //define problems
  ierr = m_field_Macro.add_problem("ELASTIC_PROBLEM_MACRO"); CHKERRQ(ierr);
  
  //set finite elements for problem
  ierr = m_field_Macro.modify_problem_add_finite_element("ELASTIC_PROBLEM_MACRO","ELASTIC_FE_MACRO"); CHKERRQ(ierr);
  ierr = m_field_Macro.modify_problem_add_finite_element("ELASTIC_PROBLEM_MACRO","ELASTIC_FE_MACRO_REL"); CHKERRQ(ierr);
  
  //set refinment level for problem
  ierr = m_field_Macro.modify_problem_ref_level_add_bit("ELASTIC_PROBLEM_MACRO",bit_level0_Macro); CHKERRQ(ierr);
  

  // ===========================================================================
  //
  // B.III. DECLARE PROBLEM
  //
  // ===========================================================================

  /*****************************************************************************
   *
   * Add entitities (by tets) to the field
   *
   ****************************************************************************/
  ierr = m_field_Macro.add_ents_to_field_by_TETs(0,"DISP_MACRO"); CHKERRQ(ierr);
  ierr = m_field_Macro.add_ents_to_field_by_TETs(0,"MESH_NODE_POSITIONS"); CHKERRQ(ierr);
  
  for(int ii=0; ii < nders; ii++ ) {
    ostringstream ss_field;
    ss_field << "DISP_MACRO" << stochastic_fields[ii];
    ierr = m_field_Macro.add_ents_to_field_by_TETs(0,ss_field.str().c_str()); CHKERRQ(ierr);
  }
  

  /*****************************************************************************
   *
   * Add finite elements entities
   *
   ****************************************************************************/
  //ierr = m_field_Macro.add_ents_to_finite_element_EntType_by_bit_ref(bit_level0_Macro,"ELASTIC_FE_MACRO",MBTET); CHKERRQ(ierr);
  ierr = m_field_Macro.add_ents_to_finite_element_by_TETs(meshset_macro_Elastic,"ELASTIC_FE_MACRO",true); CHKERRQ(ierr);
  ierr = m_field_Macro.add_ents_to_finite_element_by_TETs(meshset_Reliability,"ELASTIC_FE_MACRO_REL",true); CHKERRQ(ierr);
  

  /*****************************************************************************
   *
   * Set applied order
   * See reference for detals:
   *   Ainsworth M. and Coyle J. (2003) Hierarchic finite element bases on
   *      unstructured tetrahedral meshes. IJNME, 58(14). pp.2103-2130.
   ****************************************************************************/
  ierr = m_field_Macro.set_field_order(0,MBTET,"DISP_MACRO",order); CHKERRQ(ierr);
  ierr = m_field_Macro.set_field_order(0,MBTRI,"DISP_MACRO",order); CHKERRQ(ierr);
  ierr = m_field_Macro.set_field_order(0,MBEDGE,"DISP_MACRO",order); CHKERRQ(ierr);
  ierr = m_field_Macro.set_field_order(0,MBVERTEX,"DISP_MACRO",1); CHKERRQ(ierr);
  
  for(int ii=0; ii < nders; ii++ ) {
    ostringstream ss_field;
    ss_field << "DISP_MACRO" << stochastic_fields[ii];
    //    cout<<ss_field.str().c_str()<<endl;
    ierr = m_field_Macro.set_field_order(0,MBTET,ss_field.str().c_str(),order_st); CHKERRQ(ierr);
    ierr = m_field_Macro.set_field_order(0,MBTRI,ss_field.str().c_str(),order_st); CHKERRQ(ierr);
    ierr = m_field_Macro.set_field_order(0,MBEDGE,ss_field.str().c_str(),order_st); CHKERRQ(ierr);
    ierr = m_field_Macro.set_field_order(0,MBVERTEX,ss_field.str().c_str(),1); CHKERRQ(ierr);
  }
  
  //
  ierr = m_field_Macro.set_field_order(0,MBTET,"MESH_NODE_POSITIONS",2); CHKERRQ(ierr);
  ierr = m_field_Macro.set_field_order(0,MBTRI,"MESH_NODE_POSITIONS",2); CHKERRQ(ierr);
  ierr = m_field_Macro.set_field_order(0,MBEDGE,"MESH_NODE_POSITIONS",2); CHKERRQ(ierr);
  ierr = m_field_Macro.set_field_order(0,MBVERTEX,"MESH_NODE_POSITIONS",1); CHKERRQ(ierr);
 
 //***************************************************************************
   //Calculate Dmat for each Guass point here
  //Calculate_RVE_Dmat_TransIso_Disp calculate_rve_dmat_TransIso(m_field_Macro);
  
  //ierr = calculate_rve_dmat_TransIso.addElasticElements("DISP_MACRO"); CHKERRQ(ierr);
  //ierr = m_field_Macro.modify_finite_element_add_field_data("ELASTIC_FE_MACRO","Wt"); CHKERRQ(ierr);
  //ierr = m_field_Macro.modify_problem_add_finite_element("ELASTIC_PROBLEM_MACRO","ELASTIC_FE_MACRO"); CHKERRQ(ierr);
 //***************************************************************************

  //ierr = MetaNeummanForces::addNeumannBCElements(m_field_Macro,"ELASTIC_PROB","DISP_MACORO"); CHKERRQ(ierr);
  ierr = MetaNeummanForces::addNeumannBCElements(m_field_Macro,"DISP_MACRO"); CHKERRQ(ierr);
  ierr = m_field_Macro.modify_problem_add_finite_element("ELASTIC_PROBLEM_MACRO","FORCE_FE"); CHKERRQ(ierr);
  

  // ===========================================================================
  //
  //  B.IV. BUILD DATABASE
  //
  // ===========================================================================
  
  //build field
  ierr = m_field_Macro.build_fields(); CHKERRQ(ierr);
  
  //build finite elemnts
  ierr = m_field_Macro.build_finite_elements(); CHKERRQ(ierr);
  
  //build adjacencies
  ierr = m_field_Macro.build_adjacencies(bit_level0_Macro); CHKERRQ(ierr);
  
  //build problem
  ierr = m_field_Macro.build_problems(); CHKERRQ(ierr);
  

  // ===========================================================================
  //
  //  B.V. MESH PARTITION
  //
  // ===========================================================================
  
  //partition
  ierr = m_field_Macro.partition_problem("ELASTIC_PROBLEM_MACRO"); CHKERRQ(ierr);
  ierr = m_field_Macro.partition_finite_elements("ELASTIC_PROBLEM_MACRO"); CHKERRQ(ierr);
  ierr = m_field_Macro.partition_ghost_dofs("ELASTIC_PROBLEM_MACRO"); CHKERRQ(ierr);
  
  //m_field_Macro.list_dofs_by_field_name("DISP_MACORO",true);
  
  //print bcs
  ierr = m_field_Macro.print_cubit_displacement_set(); CHKERRQ(ierr);
  ierr = m_field_Macro.print_cubit_force_set(); CHKERRQ(ierr);
  //print block sets with materials
  ierr = m_field_Macro.print_cubit_materials_set(); CHKERRQ(ierr);

 
  // ===========================================================================
  //
  //  C. RELIABILITY ANLYSIS
  //
  // ===========================================================================
  
  cout<<"\n\n";
  cout<<"///////////////////////////////////////////////////////////////////\n";
  cout<<"//                                                               //\n";
  cout<<"//           Reliability calculation starts from here!           //\n";
  cout<<"//                                                               //\n";
  cout<<"/////////////////////////////////////////////////////////////////\n\n";
  
  //////////////////////////////////////////////////////////////////////////////
  //                                                                          //
  //        STEP 1: read inputs from data file                                //
  //                a) probability data                                       //
  //                 .1) marginal distribution type                           //
  //                 .2) parameters                                           //
  //                 .3) correlation matrix                                   //
  //                b) limit state function                                   //
  //                 .1) indicator of how limit-state function is given       //
  //                     1:                                                   //
  //                c) analysis option                                        //
  //                                                                          //
  //////////////////////////////////////////////////////////////////////////////

  double beta;                    // Reliability index
  
  /*
   *  Set reliability calculation options
   */
  Reliability_Options ReliabOpt;
  
  ReliabOpt.echo_flag     = 1;       // Program interactive mode, 0: silent mode
  // For FORM
  ReliabOpt.istep_max     = 2000;    // Maximum number of interation allowed in the search algorithm
  ReliabOpt.e1            = 0.001;   // Tolerance on how close design point is to limit-state surface
  ReliabOpt.e2            = 0.001;   // Tolerance on how accurately the gradient points towards the origin
  ReliabOpt.step_code     = 0.025;       // 0: step size by Armijo rule, otherwise: given value is the step size
  ReliabOpt.Recorded_u    = 1;       // 0: u-vector not recorded at all iterations, 1: u-vector recorded at all iterations
  ReliabOpt.Recorded_x    = 1;       // 0: x-vector not recorded at all iterations, 1: x-vector recorded at all iterations
  ReliabOpt.Recorded_beta = 1;
  ReliabOpt.grad_G        = "PSFEM"; // "PSFEM": perturbation, "DDM": direct differentiation, 'ADM': automatic differentiation
  
  int    echo_flag = ReliabOpt.echo_flag;
  double e1        = ReliabOpt.e1;
  double e2        = ReliabOpt.e2;
  int    istep_max = ReliabOpt.istep_max;
  double step_code = ReliabOpt.step_code;
  int    beta_flag = ReliabOpt.Recorded_beta;
  
  /*
   *  Read inputs' statistical properties from file to insert into <probdata>
   */
  Stochastic_Model probdata;
  
  // Import data from textile file
  ImportProbData readprobdata;
  ierr = readprobdata.ProbdataFileIn(); CHKERRQ(ierr);
  probdata.marg        = readprobdata.MargProb;
  probdata.correlation = readprobdata.CorrMat;
  probdata.num_vars    = readprobdata.NumVars;
  probdata.MatStrength = readprobdata.MatStrength;
  probdata.NameVars    = readprobdata.NameVars;
  
  
  //////////////////////////////////////////////////////////////////////////////
  //                                                                          //
  //        STEP 2: use probability transformation to obtain y in             //
  //                standard normal distribution space                        //
  //                y = norminv(F(x))                                         //
  //                                                                          //
  //////////////////////////////////////////////////////////////////////////////
  
  NatafTransformation my_nataf_transformation;
  ierr = my_nataf_transformation.ModCorrMat_Empirical(probdata.num_vars,probdata.marg,probdata.correlation,probdata.mod_correlation); CHKERRQ(ierr);
  
  /*
   * Perform Cholesky decomposition for the modified correlation matrix
   *    A = LL'
   */
  ublas::triangular_matrix<double, ublas::lower> Lmat(probdata.num_vars,probdata.num_vars);
  cholesky_decompose(probdata.mod_correlation,Lmat);
  probdata.Lo.resize(probdata.num_vars,probdata.num_vars);
  probdata.Lo = Lmat;
  
  //cout<<"Cholesky decomposed matrix: "<<Lmat<<endl;
  
  // Compute the inverse of Lo
  probdata.inv_Lo.resize(probdata.num_vars,probdata.num_vars);
  bool singular = false;
  probdata.inv_Lo = gjinverse(probdata.Lo,singular);
  
  
  
  //////////////////////////////////////////////////////////////////////////////
  //                                                                          //
  //        STEP 3: select an initial checking point, x                       //
  //                                                                          //
  //////////////////////////////////////////////////////////////////////////////
  
  ublas::vector<double> x(probdata.num_vars);
  for (int i=0; i<probdata.num_vars; i++) {
    x(i) = probdata.marg(i,3);
  }
  
  ublas::vector<double> u;
  ierr = my_nataf_transformation.x_to_u(x,probdata.num_vars,probdata.marg,probdata.inv_Lo,u); CHKERRQ(ierr);
  
  
  
  // ===========================================================================
  //
  //  C. POST-PROCESSING
  //
  // ===========================================================================
  // Declare matrix for stroring RVE constitutive matrix & its derivatives
  ublas::matrix<double> Dmat;
  ublas::matrix<double> Dmat_r_Em, Dmat_r_NUm, Dmat_r_Ep, Dmat_r_Ez;
  ublas::matrix<double> Dmat_r_NUp, Dmat_r_NUpz, Dmat_r_Gzp;
  
  ublas::matrix<double> Dmat_rs_EmEm, Dmat_rs_NUmNUm, Dmat_rs_EpEp, Dmat_rs_EzEz;
  ublas::matrix<double> Dmat_rs_NUpNUp, Dmat_rs_NUpzNUpz, Dmat_rs_GzpGzp;
  

  /*
   //Reading and writing binary files
   if(pcomm->rank()==0) {
   int fd;
   PetscViewer view_out;
   PetscViewerBinaryOpen(PETSC_COMM_WORLD,"input.dat",FILE_MODE_WRITE,&view_out);
   PetscViewerBinaryGetDescriptor(view_out,&fd);
   PetscBinaryWrite(fd,&Dmat(0,0),36,PETSC_DOUBLE,PETSC_FALSE);
   PetscViewerDestroy(&view_out);
   }
   
   ublas::matrix<FieldData> Dmat1;
   Dmat1.resize(6,6); Dmat1.clear();
   cout<< "Dmat1 Before Reading= "<<Dmat1<<endl;
   if(pcomm->rank()==0) {
   int fd;
   PetscViewer view_in;
   PetscViewerBinaryOpen(PETSC_COMM_WORLD,"input.dat",FILE_MODE_READ,&view_in);
   PetscViewerBinaryGetDescriptor(view_in,&fd);
   PetscBinaryRead(fd,&Dmat1(0,0),36,PETSC_DOUBLE);
   PetscViewerDestroy(&view_in);
   }
   cout<< "Dmat1 After Reading= "<<Dmat1<<endl;
   

  // ===========================================================================
  //
  //  B.VII. OUTPUT
  //
  // ===========================================================================
  
  
  ProjectionFieldOn10NodeTet ent_method_on_10nodeTet(m_field_Macro,"DISP_MACRO",true,false,"DISP_MACRO");
  ierr = m_field_Macro.loop_dofs("DISP_MACRO",ent_method_on_10nodeTet); CHKERRQ(ierr);
  ent_method_on_10nodeTet.set_nodes = false;
  ierr = m_field_Macro.loop_dofs("DISP_MACRO",ent_method_on_10nodeTet); CHKERRQ(ierr);
  
  if(pcomm->rank()==0) {
    EntityHandle out_meshset;
    rval = moab_Macro.create_meshset(MESHSET_SET,out_meshset); CHKERR_PETSC(rval);
    ierr = m_field_Macro.get_problem_finite_elements_entities("ELASTIC_PROBLEM_MACRO","ELASTIC_FE_MACRO",out_meshset); CHKERRQ(ierr);
    rval = moab_Macro.write_file("out.vtk","VTK","",&out_meshset,1); CHKERR_PETSC(rval);
    rval = moab_Macro.delete_entities(&out_meshset,1); CHKERR_PETSC(rval);
  }
  */
  
  
  
  //////////////////////////////////////////////////////////////////////////////
  //                                                                          //
  //        STEP 4: Perform iterative loop to find design point               //
  //                                                                          //
  //////////////////////////////////////////////////////////////////////////////
  
  /*
   *  Set parameters for the iterative loop
   */
  int    istep = 1;     // Initialize iterative counter
  int    conv_flag = 0; // Convergence is achieved when this flag is set to 1
  
  double val_G, val_G0;                                         // LSF for given inputs
  double step_size;                                             // Step size
  ublas::vector<double> grad_g(probdata.num_vars); // Gradient of LSF in x space
  ublas::vector<double> grad_G(probdata.num_vars); // Gradient of LSF in u space
  ublas::vector<double> alpha(probdata.num_vars);  // Direction cosine vector
  ublas::matrix<double> dudx(probdata.num_vars,probdata.num_vars);
  ublas::matrix<double> inv_dudx(probdata.num_vars,probdata.num_vars);
  ublas::vector<double> u_dir(probdata.num_vars);  // Direction
  ublas::vector<double> u_new(probdata.num_vars);  // New trial of checking point
  
  //
  ublas::matrix<double> StressGP;
  ublas::matrix<double> StressGP_r_Em;
  ublas::matrix<double> StressGP_r_NUm;
  ublas::matrix<double> StressGP_r_Ep;
  ublas::matrix<double> StressGP_r_Ez;
  ublas::matrix<double> StressGP_r_NUp;
  ublas::matrix<double> StressGP_r_NUpz;
  ublas::matrix<double> StressGP_r_Gzp;
  
  /*
   *  Start iteration
   */
  
  clock_t rel_t1, rel_t2;
  double rel_calc_time;
  rel_t1 = clock();
  
  FE2_Macro_Solver Solve_FE2_Problem;
  LimitStateFunction TheLSF;
  cout<<beta_flag<<endl;
  
  ofstream BetaFile;
  if (beta_flag == 1) {
    BetaFile.open("//mnt//home//Dropbox//DURACOMP_Cal//009_MoFEM//04_ReliabilityAnalysis//Result_Beta.txt",ofstream::out);
  }
  
  do {
    
    cout<<"\n\n*************************************************\n*\n";
    cout<<"*    This is "<<istep<<" step!"<<"\n*\n";
    cout<<"*************************************************\n";
    
    if (echo_flag) {
      cout<<"-------------------------------- \n";
      cout<<"Now carrying out iteration number: \t"<<istep<<endl;
    }
  
  
    ////////////////////////////////////////////////////////////////////////////
    //                                                                        //
    //        STEP 5: calculate Jacobian, J = dy/dx                           //
    //                                                                        //
    ////////////////////////////////////////////////////////////////////////////
  
    /*
     * Transformation from u to x space
     */
    double detj=1.0;
    x.resize(probdata.num_vars); x.clear();cout<<"u value: "<<u<<endl;
    ierr = my_nataf_transformation.u_to_x(u,probdata.num_vars,probdata.marg,probdata.Lo,x,detj); CHKERRQ(ierr);
    cout<<"Determinant of Jacobian: "<<detj<<"\t x value: "<<x<<endl;
    
    // Jacobian
    dudx.resize(probdata.num_vars,probdata.num_vars); dudx.clear();
    ierr = my_nataf_transformation.Jacobian_u_x(x,u,probdata.num_vars,probdata.marg,probdata.Lo,probdata.inv_Lo,dudx); CHKERRQ(ierr);
    
    inv_dudx.resize(probdata.num_vars,probdata.num_vars); inv_dudx.clear();
    inv_dudx = gjinverse(dudx,singular);
    
    /*
     * Evaluate limit-state function and its gradient
     */
    
    ierr = Solve_FE2_Problem.Calculate_RVEDmat(m_field_RVE, nvars, nders, stochastic_fields,x,probdata.num_vars,probdata.NameVars); CHKERRQ(ierr);
    ierr = Solve_FE2_Problem.Macro_FE_Solver(m_field_Macro, nvars, nders, stochastic_fields); CHKERRQ(ierr);
    
    Dmat             = Solve_FE2_Problem.Dmat; cout<<"Dmat: "<<Dmat<<"\n\n";
    Dmat_r_Em        = Solve_FE2_Problem.Dmat_r_Em; //cout<<"Dmat_r_Em: "<<Dmat_r_Em<<"\n\n";
    Dmat_r_NUm       = Solve_FE2_Problem.Dmat_r_NUm;
    Dmat_r_Ep        = Solve_FE2_Problem.Dmat_r_Ep;
    Dmat_r_Ez        = Solve_FE2_Problem.Dmat_r_Ez;
    Dmat_r_NUp       = Solve_FE2_Problem.Dmat_r_NUp;
    Dmat_r_NUpz      = Solve_FE2_Problem.Dmat_r_NUpz;
    Dmat_r_Gzp       = Solve_FE2_Problem.Dmat_r_Gzp;
    
    Dmat_rs_EmEm     = Solve_FE2_Problem.Dmat_rs_EmEm;
    Dmat_rs_NUmNUm   = Solve_FE2_Problem.Dmat_rs_NUmNUm;
    Dmat_rs_EpEp     = Solve_FE2_Problem.Dmat_rs_EpEp;
    Dmat_rs_EzEz     = Solve_FE2_Problem.Dmat_rs_EzEz;
    Dmat_rs_NUpNUp   = Solve_FE2_Problem.Dmat_rs_NUpNUp;
    Dmat_rs_NUpzNUpz = Solve_FE2_Problem.Dmat_rs_NUpzNUpz;
    Dmat_rs_GzpGzp   = Solve_FE2_Problem.Dmat_rs_GzpGzp;//cout<<"Dmat_rs_GzpGzp: "<<Dmat_rs_GzpGzp<<"\n\n";
    
    // Calculate the zeroth-order stress at Gauss points for specific element(s)
    FE2_PostProcStressForReliability_Zeroth Calc_Stress(m_field_Macro,"DISP_MACRO",Dmat);
    ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_FE_MACRO_REL",Calc_Stress);  CHKERRQ(ierr);
    StressGP.clear(); StressGP = Calc_Stress.StressGP; cout<<"Stress at GP: "<<StressGP<<endl;
  
    // Calculate the first-order partial derivative stress at Gauss points for specific element(s)
    // with respect to Em
    FE2_PostProcStressForReliability_First Calc_Stress_Em(m_field_Macro,"DISP_MACRO","DISP_MACRO_r_Em",Dmat,Dmat_r_Em);
    ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_FE_MACRO_REL",Calc_Stress_Em);  CHKERRQ(ierr);
    StressGP_r_Em.clear(); StressGP_r_Em = Calc_Stress_Em.StressGP_r;
    
    // with respect to NUm
    FE2_PostProcStressForReliability_First Calc_Stress_NUm(m_field_Macro,"DISP_MACRO","DISP_MACRO_r_NUm",Dmat,Dmat_r_NUm);
    ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_FE_MACRO_REL",Calc_Stress_NUm);  CHKERRQ(ierr);
    StressGP_r_NUm.clear(); StressGP_r_NUm = Calc_Stress_NUm.StressGP_r;
    
    
    // with respect to NUp
    FE2_PostProcStressForReliability_First Calc_Stress_NUp(m_field_Macro,"DISP_MACRO","DISP_MACRO_r_NUp",Dmat,Dmat_r_NUp);
    ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_FE_MACRO_REL",Calc_Stress_NUp);  CHKERRQ(ierr);
    StressGP_r_NUp.clear(); StressGP_r_NUp = Calc_Stress_NUp.StressGP_r;
    
    // with respect to NUpz
    FE2_PostProcStressForReliability_First Calc_Stress_NUpz(m_field_Macro,"DISP_MACRO","DISP_MACRO_r_NUpz",Dmat,Dmat_r_NUpz);
    ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_FE_MACRO_REL",Calc_Stress_NUpz);  CHKERRQ(ierr);
    StressGP_r_NUpz.clear(); StressGP_r_NUpz = Calc_Stress_NUpz.StressGP_r;
    
    // with respect to Ep
    FE2_PostProcStressForReliability_First Calc_Stress_Ep(m_field_Macro,"DISP_MACRO","DISP_MACRO_r_Ep",Dmat,Dmat_r_Ep);
    ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_FE_MACRO_REL",Calc_Stress_Ep);  CHKERRQ(ierr);
    StressGP_r_Ep.clear(); StressGP_r_Ep = Calc_Stress_Ep.StressGP_r;
    
    // with respect to Ez
    FE2_PostProcStressForReliability_First Calc_Stress_Ez(m_field_Macro,"DISP_MACRO","DISP_MACRO_r_Ez",Dmat,Dmat_r_Ez);
    ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_FE_MACRO_REL",Calc_Stress_Ez);  CHKERRQ(ierr);
    StressGP_r_Ez.clear(); StressGP_r_Ez = Calc_Stress_Ez.StressGP_r;
    
    // with respect to Gzp
    FE2_PostProcStressForReliability_First Calc_Stress_Gzp(m_field_Macro,"DISP_MACRO","DISP_MACRO_r_Gzp",Dmat,Dmat_r_Gzp);
    ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_FE_MACRO_REL",Calc_Stress_Gzp);  CHKERRQ(ierr);
    StressGP_r_Gzp.clear(); StressGP_r_Gzp = Calc_Stress_Gzp.StressGP_r;
    
    // Evaluate LSF and its gradient
    grad_g.resize(probdata.num_vars); grad_g.clear();
    //ierr = TheLSF.gfun(x,val_G,grad_g); CHKERRQ(ierr);
    ierr = TheLSF.gfun_ply_TW(x,probdata.NameVars,probdata.MatStrength,
                              StressGP,//probdata.MatStrength,StressGP,
                              StressGP_r_Em,StressGP_r_NUm,
                              StressGP_r_NUp,StressGP_r_NUpz,
                              StressGP_r_Ep,StressGP_r_Ez,StressGP_r_Gzp,
                              val_G,grad_g); CHKERRQ(ierr);
    
    grad_G.resize(probdata.num_vars); grad_G.clear();
    grad_G = prod(grad_g,inv_dudx);
    
    cout<<"LSF value is \t"<<val_G<<endl;
    cout<<"Gradient of LSF is \t"<<grad_g<<"\t"<<grad_G<<endl;
    
    /*
     *  Set scale parameter G0 and inform about structural response
     */
    if (istep == 1) {
      val_G0 = val_G;
      if (echo_flag) {
        cout<<"Value of limit-state function in the first step: \t"<<val_G<<endl;
      }
    }
    
    ////////////////////////////////////////////////////////////////////////////
    //                                                                        //
    //        STEP 6: compute the direction cosine                            //
    //                a) dg/dy = inv(J) dG/dx                                 //
    //                b) alpha_i = dg/dy_i / norm(dg/dy)                      //
    //                c) reliability index estimate: beta = sqrt(norm(y))     //
    //                                                                        //
    ////////////////////////////////////////////////////////////////////////////
  
    /*
     * Compute direction cosine alpha vector
     */
    alpha.resize(probdata.num_vars); alpha.clear();
    alpha = -grad_G/norm_2(grad_G);
    cout<<"Direction cosine "<<alpha<<endl;
    
    /*
     * Check convergence
     */
    if (((abs(val_G/val_G0)<e1) && (norm_2(u-inner_prod(alpha,u)*alpha))) || (istep == istep_max)) {
      conv_flag = 1;
    }
    
    ////////////////////////////////////////////////////////////////////////////
    //                                                                        //
    //        STEP 7: convergence check and compute new trial point:          //
    //                7.1 convergence check                                   //
    //                    (a) design point                                    //
    //                    (b) estimate of reliability index                   //
    //                7.2 new trial point                                     //
    //                           y_n = -alpha[beta + g/norm(dg/dy)]           //
    //                                                                        //
    ////////////////////////////////////////////////////////////////////////////
    
    /*
     * Take a step if convergence is not achieved
     */
    if (conv_flag == 0) {
      // Determine search direction
      u_dir.resize(probdata.num_vars); u_dir.clear();
      search_dir(val_G,grad_G,u,u_dir);
      // Determine step size
      if (step_code != 0) {
        step_size = step_code;
      }
      else {
        cout<<"\n!!!!!!!!!!!!!!!!!!!!!!!!!!!WARNING!!!!!!!!!!!!!!!!!!!!!!!!!!!";
        cout<<"\nArmijo rule will be used to determine step size for setting new trial point, \n";
        cout<<"but it isn't available yet in current version!!! \n";
        cout<<"Re-assign a nonzero value to step_code.\n";
        cout<<"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"<<endl;
        conv_flag = 1;
      }
      
      cout<<"\n\nReliability index estimate at "<<istep<<" is: "<<((val_G/norm_2(grad_G)) + inner_prod(alpha, u));
      cout<<"\t"<<norm_2(u)<<endl;
      //cout<<"design point: "<<u<<endl;
      
      // Write reliability index in file
      if (beta_flag == 1) {
        BetaFile<<setprecision(15)<<inner_prod(alpha,u)<<"\t"<<val_G<<"\n";
      }
      
      
      // Determin new trial point
      u_new.resize(probdata.num_vars); u_new.clear();
      u_new = u + step_size*u_dir;  // when step_size is 1, it is HLRF search algorithm
      // Prepare for a new round in the loop
      u.resize(probdata.num_vars); u.clear();
      u = u_new;
      istep = istep + 1;
    }
  
    rel_t2 = clock();
    rel_calc_time  = (double)(rel_t2 - rel_t1)/CLOCKS_PER_SEC;
    rel_t1 = rel_t2;
    cout<<"Elapsed time at this step is: "<<rel_calc_time<<" seconds.\n";
    
  
  } while (conv_flag == 0); // end of while
  
  // Close beta value writting file
  if (beta_flag == 1) { BetaFile.close(); }
  
  //////////////////////////////////////////////////////////////////////////////
  //                                                                          //
  //        STEP 8: Calculate the estimate of reliability index               //
  //                                                                          //
  //////////////////////////////////////////////////////////////////////////////

  beta = inner_prod(alpha,u);
  
  
  //////////////////////////////////////////////////////////////////////////////
  //                                                                          //
  //                               FINISH                                     //
  //                                                                          //
  //////////////////////////////////////////////////////////////////////////////
  
  finish_time = clock();
  total_time  = (double)(finish_time - start_time)/CLOCKS_PER_SEC;
  
  cout<<"\n\n*************************************************\n*\n";
  
  
  if (istep == istep_max) {
    cout<<"*  The maximum number of iteration is reached.\n";
  }
  
  cout<<"*  Optimal reliability index (beta) is: "<<beta<<endl;
  cout<<"*  Elapsed time is "<<total_time<<" seconds.\n";
  cout<<"*  The program finishes !!! \n*\n";
  cout<<"*************************************************"<<endl;

  ierr = PetscFinalize(); CHKERRQ(ierr);
  
}
