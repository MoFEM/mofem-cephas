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
#include <memory>
#include <MoFEM.hpp>
using namespace MoFEM;

#include <MethodForForceScaling.hpp>

#include <DirichletBC.hpp>
#include <Projection10NodeCoordsOnField.hpp>

#include <petsctime.h>

#include <SurfacePressure.hpp>
#include <FEMethod_LowLevelStudent.hpp>
#include <FEMethod_UpLevelStudent.hpp>

#include <PostProcOnRefMesh.hpp>
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
#include "MaterialConstitutiveMatrix_RuleOfMixture.hpp"
#include "MaterialConstitutiveMatrix_MoriTanaka.hpp"
#include "Trans_Iso_Rhs_r_PSFEM.hpp"
#include "Trans_Iso_Rhs_rs_PSFEM.hpp"
#include <FE2_ElasticFEMethod.hpp>

#include <FE2_Rhs_r_PSFEM.hpp>
#include <FE2_Rhs_rs_PSFEM.hpp>

#include <Reliability_SurfacePressure.hpp>

#include <Reliability_Functions.hpp>
#include <FE2_PostProcStressForReliability.hpp>

#include <FE2_Macro_Solver.hpp>

#include <Single_FE_Rhs_r_PSFEM.hpp>
#include <Single_FE_Solver_Laminate.hpp>

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


ErrorCode rval;
PetscErrorCode ierr;

//#include <Reliability_Methods.hpp>


static char help[] = "...\n\n";

const char* args[] = {
  "_r"                                                // 1st order
};

int nvars = 1;    // number of variables
int nders = 1;   // number of partial derivatives (firsr- and second- order)
vector<string> stochastic_fields(args, args + 1);

int main(int argc, char *argv[]) {
  
  clock_t start_time, finish_time;
  double total_time;
  start_time = clock();
  
  
  //============================================================================
  //
  //  B. Macro Problem
  //
  //============================================================================
  
  PetscInitialize(&argc,&argv,(char *)0,help);
  
  int rank;
  MPI_Comm_rank(PETSC_COMM_WORLD,&rank);
  
  moab::Core mb_instance_Macro;
  Interface& moab_Macro = mb_instance_Macro;
  
  // ===========================================================================
  //
  // B. I. READ MESH DATA AND FINITE ELEMENT ANALYSIS CONTROL PARAMETERS FROM FILE
  //
  // ===========================================================================
  
  PetscBool flg = PETSC_TRUE;
  char mesh_file_name[255];
  
  //ParallelComm* pcomm_Macro = ParallelComm::get_pcomm(&moab_Macro,MYPCOMM_INDEX);
  //if(pcomm == NULL) pcomm_Macro =  new ParallelComm(&moab_Macro,PETSC_COMM_WORLD);
  
  /*****************************************************************************
   *
   * Read parameters from line command
   *
   ****************************************************************************/
  ierr = PetscOptionsGetString(PETSC_NULL,"-my_file_macro",mesh_file_name,255,&flg); CHKERRQ(ierr);
  if(flg != PETSC_TRUE) {
    SETERRQ(PETSC_COMM_SELF,1,"*** ERROR -my_file_macro (MESH FILE NEEDED)");
  }
  
  PetscInt order_Macro;
  ierr = PetscOptionsGetInt(PETSC_NULL,"-my_order_Macro",&order_Macro,&flg); CHKERRQ(ierr);
  if(flg != PETSC_TRUE) {
    order_Macro = 1;
  }
  
  PetscInt NO_Layers;
  ierr = PetscOptionsGetInt(PETSC_NULL,"-my_NO_Layers",&NO_Layers,&flg); CHKERRQ(ierr);
  if(flg != PETSC_TRUE) {
    NO_Layers = 1;
  }
  cout<<"\n\nNumber of layers: "<<NO_Layers<<endl;
  
  /*****************************************************************************
   *
   * Transfer mesh data to MOAB database
   *
   ****************************************************************************/
  
  const char *option;
  option = "";//"PARALLEL=BCAST;";//;DEBUG_IO";
  rval = moab_Macro.load_file(mesh_file_name, 0, option); CHKERR_PETSC(rval);
  ParallelComm* pcomm = ParallelComm::get_pcomm(&moab_Macro,MYPCOMM_INDEX);
  if(pcomm == NULL) pcomm =  new ParallelComm(&moab_Macro,PETSC_COMM_WORLD);
  
  //Create MoFEM (Joseph) database
  MoFEM::Core core_Macro(moab_Macro);
  FieldInterface& m_field_Macro = core_Macro;

  //ref meshset ref level 0
  ierr = m_field_Macro.seed_ref_level_3D(0,0); CHKERRQ(ierr);

  // stl::bitset see for more details
  BitRefLevel bit_level0_Macro;
  bit_level0_Macro.set(0);

  EntityHandle meshset_level0;
  rval = moab_Macro.create_meshset(MESHSET_SET,meshset_level0); CHKERR_PETSC(rval);
  ierr = m_field_Macro.seed_ref_level_3D(0,bit_level0_Macro); CHKERRQ(ierr);
  ierr = m_field_Macro.get_entities_by_ref_level(bit_level0_Macro,BitRefLevel().set(),meshset_level0); CHKERRQ(ierr);
  
  Range TetsInBlock_1st_Ply, TetsInBlock_2nd_Ply, TetsInBlock_3rd_Ply, TetsInBlock_4th_Ply, TetsInBlock_Reliability;
  for(_IT_CUBITMESHSETS_BY_SET_TYPE_FOR_LOOP_(m_field_Macro,BLOCKSET,it)){
    if(it->get_name() == "MAT_ELASTIC_TRANSISO_First") {
      rval = moab_Macro.get_entities_by_type(it->meshset, MBTET,TetsInBlock_1st_Ply,true); CHKERR_PETSC(rval);
    }
    if(it->get_name() == "MAT_ELASTIC_TRANSISO_Second") {
      rval = moab_Macro.get_entities_by_type(it->meshset, MBTET,TetsInBlock_2nd_Ply,true); CHKERR_PETSC(rval);
    }
    if(it->get_name() == "MAT_ELASTIC_TRANSISO_Third") {
      rval = moab_Macro.get_entities_by_type(it->meshset, MBTET,TetsInBlock_3rd_Ply,true); CHKERR_PETSC(rval);
    }
    if(it->get_name() == "MAT_ELASTIC_TRANSISO_Fourth") {
      rval = moab_Macro.get_entities_by_type(it->meshset, MBTET,TetsInBlock_4th_Ply,true); CHKERR_PETSC(rval);
    }
    if(it->get_name() == "RELIABILITY") {
      rval = moab_Macro.get_entities_by_type(it->meshset, MBTET,TetsInBlock_Reliability,true); CHKERR_PETSC(rval);
    }
  }
  
  
  // ===========================================================================
  //
  // B.II. DEFINE PROBLEM
  //
  // ===========================================================================
  
  int field_rank=3;

  /*****************************************************************************
   *
   * Add stochastic field
   * (total 14 field for 1st and 2nd order stochastic PSFEM)
   *
   ****************************************************************************/
  ierr = m_field_Macro.add_field("DISP_MACRO",H1,3,MF_ZERO); CHKERRQ(ierr);
  ierr = m_field_Macro.add_field("MESH_NODE_POSITIONS",H1,3,MF_ZERO); CHKERRQ(ierr);
  
  // Stochastic fields for perturbation methods at macroscale
  for(int ii=0; ii < nvars; ii++ ) {
    ostringstream ss_field;
    ss_field << "DISP_MACRO" << stochastic_fields[ii];
    ierr = m_field_Macro.add_field(ss_field.str().c_str(),H1,field_rank,MF_ZERO); CHKERRQ(ierr);
  }
  
  /*****************************************************************************
   *
   * Create finite element for the defined fields
   *
   ****************************************************************************/
  // First layer
  ierr = m_field_Macro.add_finite_element("ELASTIC_1st_Ply"); CHKERRQ(ierr);
  // Second layer
  if (NO_Layers > 1) {
    ierr = m_field_Macro.add_finite_element("ELASTIC_2nd_Ply"); CHKERRQ(ierr);
  }
  // Third layer
  if (NO_Layers > 2) {
    ierr = m_field_Macro.add_finite_element("ELASTIC_3rd_Ply"); CHKERRQ(ierr);
  }
  // Fourth layer
  if (NO_Layers > 3) {
    ierr = m_field_Macro.add_finite_element("ELASTIC_4th_Ply"); CHKERRQ(ierr);
  }
  ierr = m_field_Macro.add_finite_element("ELASTIC_FE_MACRO_REL"); CHKERRQ(ierr);
 

  /*****************************************************************************
   *
   * set field data which finite element use
   *
   ****************************************************************************/
  // First layer
  //Define rows/cols and element data
  ierr = m_field_Macro.modify_finite_element_add_field_row("ELASTIC_1st_Ply","DISP_MACRO"); CHKERRQ(ierr);
  ierr = m_field_Macro.modify_finite_element_add_field_col("ELASTIC_1st_Ply","DISP_MACRO"); CHKERRQ(ierr);
  ierr = m_field_Macro.modify_finite_element_add_field_data("ELASTIC_1st_Ply","DISP_MACRO"); CHKERRQ(ierr);
  ierr = m_field_Macro.modify_finite_element_add_field_data("ELASTIC_1st_Ply","MESH_NODE_POSITIONS"); CHKERRQ(ierr);
  
  for(int ii=0; ii < nvars; ii++ ) {
    ostringstream ss_field;
    ss_field << "DISP_MACRO" << stochastic_fields[ii];
    ierr = m_field_Macro.modify_finite_element_add_field_data("ELASTIC_1st_Ply",ss_field.str().c_str()); CHKERRQ(ierr);
  }
//  ierr = m_field_Macro.modify_finite_element_add_field_data("ELASTIC_1st_Ply","DISP_MACRO_r_F"); CHKERRQ(ierr);
//  ierr = m_field_Macro.modify_finite_element_add_field_data("ELASTIC_1st_Ply","DISP_MACRO_r_Theta"); CHKERRQ(ierr);
//  ierr = m_field_Macro.modify_finite_element_add_field_data("ELASTIC_1st_Ply","DISP_MACRO_r_Theta_1st_Ply"); CHKERRQ(ierr);
  
  
  // Second layer
  if (NO_Layers > 1) {
    //Define rows/cols and element data
    ierr = m_field_Macro.modify_finite_element_add_field_row("ELASTIC_2nd_Ply","DISP_MACRO"); CHKERRQ(ierr);
    ierr = m_field_Macro.modify_finite_element_add_field_col("ELASTIC_2nd_Ply","DISP_MACRO"); CHKERRQ(ierr);
    ierr = m_field_Macro.modify_finite_element_add_field_data("ELASTIC_2nd_Ply","DISP_MACRO"); CHKERRQ(ierr);
    ierr = m_field_Macro.modify_finite_element_add_field_data("ELASTIC_2nd_Ply","MESH_NODE_POSITIONS"); CHKERRQ(ierr);
    
//    ierr = m_field_Macro.modify_finite_element_add_field_data("ELASTIC_1st_Ply","DISP_MACRO_r_Theta_2nd_Ply"); CHKERRQ(ierr);
    
    for(int ii=0; ii < nvars; ii++ ) {
      ostringstream ss_field;
      ss_field << "DISP_MACRO" << stochastic_fields[ii];
      ierr = m_field_Macro.modify_finite_element_add_field_data("ELASTIC_2nd_Ply",ss_field.str().c_str()); CHKERRQ(ierr);
    }
//    ierr = m_field_Macro.modify_finite_element_add_field_data("ELASTIC_2nd_Ply","DISP_MACRO_r_F"); CHKERRQ(ierr);
//    ierr = m_field_Macro.modify_finite_element_add_field_data("ELASTIC_2nd_Ply","DISP_MACRO_r_Theta"); CHKERRQ(ierr);
//    ierr = m_field_Macro.modify_finite_element_add_field_data("ELASTIC_2nd_Ply","DISP_MACRO_r_Theta_1st_Ply"); CHKERRQ(ierr);
//    ierr = m_field_Macro.modify_finite_element_add_field_data("ELASTIC_2nd_Ply","DISP_MACRO_r_Theta_2nd_Ply"); CHKERRQ(ierr);
  }
  // Third layer
  if (NO_Layers > 2) {
    //Define rows/cols and element data
    ierr = m_field_Macro.modify_finite_element_add_field_row("ELASTIC_3rd_Ply","DISP_MACRO"); CHKERRQ(ierr);
    ierr = m_field_Macro.modify_finite_element_add_field_col("ELASTIC_3rd_Ply","DISP_MACRO"); CHKERRQ(ierr);
    ierr = m_field_Macro.modify_finite_element_add_field_data("ELASTIC_3rd_Ply","DISP_MACRO"); CHKERRQ(ierr);
    ierr = m_field_Macro.modify_finite_element_add_field_data("ELASTIC_3rd_Ply","MESH_NODE_POSITIONS"); CHKERRQ(ierr);
    
//    ierr = m_field_Macro.modify_finite_element_add_field_data("ELASTIC_1st_Ply","DISP_MACRO_r_Theta_3rd_Ply"); CHKERRQ(ierr);
//    ierr = m_field_Macro.modify_finite_element_add_field_data("ELASTIC_2nd_Ply","DISP_MACRO_r_Theta_3rd_Ply"); CHKERRQ(ierr);
    
    for(int ii=0; ii < nvars; ii++ ) {
      ostringstream ss_field;
      ss_field << "DISP_MACRO" << stochastic_fields[ii];
      ierr = m_field_Macro.modify_finite_element_add_field_data("ELASTIC_3rd_Ply",ss_field.str().c_str()); CHKERRQ(ierr);
    }
//    ierr = m_field_Macro.modify_finite_element_add_field_data("ELASTIC_3rd_Ply","DISP_MACRO_r_F"); CHKERRQ(ierr);
//    ierr = m_field_Macro.modify_finite_element_add_field_data("ELASTIC_3rd_Ply","DISP_MACRO_r_Theta"); CHKERRQ(ierr);
//    ierr = m_field_Macro.modify_finite_element_add_field_data("ELASTIC_3rd_Ply","DISP_MACRO_r_Theta_1st_Ply"); CHKERRQ(ierr);
//    ierr = m_field_Macro.modify_finite_element_add_field_data("ELASTIC_3rd_Ply","DISP_MACRO_r_Theta_2nd_Ply"); CHKERRQ(ierr);
//    ierr = m_field_Macro.modify_finite_element_add_field_data("ELASTIC_3rd_Ply","DISP_MACRO_r_Theta_3rd_Ply"); CHKERRQ(ierr);
  }
  // Fourth layer
  if (NO_Layers > 3) {
    //Define rows/cols and element data
    ierr = m_field_Macro.modify_finite_element_add_field_row("ELASTIC_4th_Ply","DISP_MACRO"); CHKERRQ(ierr);
    ierr = m_field_Macro.modify_finite_element_add_field_col("ELASTIC_4th_Ply","DISP_MACRO"); CHKERRQ(ierr);
    ierr = m_field_Macro.modify_finite_element_add_field_data("ELASTIC_4th_Ply","DISP_MACRO"); CHKERRQ(ierr);
    ierr = m_field_Macro.modify_finite_element_add_field_data("ELASTIC_4th_Ply","MESH_NODE_POSITIONS"); CHKERRQ(ierr);
    
//    ierr = m_field_Macro.modify_finite_element_add_field_data("ELASTIC_1st_Ply","DISP_MACRO_r_Theta_4th_Ply"); CHKERRQ(ierr);
//    ierr = m_field_Macro.modify_finite_element_add_field_data("ELASTIC_2nd_Ply","DISP_MACRO_r_Theta_4th_Ply"); CHKERRQ(ierr);
//    ierr = m_field_Macro.modify_finite_element_add_field_data("ELASTIC_3rd_Ply","DISP_MACRO_r_Theta_4th_Ply"); CHKERRQ(ierr);
    
    for(int ii=0; ii < nvars; ii++ ) {
      ostringstream ss_field;
      ss_field << "DISP_MACRO" << stochastic_fields[ii];
      ierr = m_field_Macro.modify_finite_element_add_field_data("ELASTIC_4th_Ply",ss_field.str().c_str()); CHKERRQ(ierr);
    }
//    ierr = m_field_Macro.modify_finite_element_add_field_data("ELASTIC_4th_Ply","DISP_MACRO_r_F"); CHKERRQ(ierr);
//    ierr = m_field_Macro.modify_finite_element_add_field_data("ELASTIC_4th_Ply","DISP_MACRO_r_Theta"); CHKERRQ(ierr);
//    ierr = m_field_Macro.modify_finite_element_add_field_data("ELASTIC_4th_Ply","DISP_MACRO_r_Theta_1st_Ply"); CHKERRQ(ierr);
//    ierr = m_field_Macro.modify_finite_element_add_field_data("ELASTIC_4th_Ply","DISP_MACRO_r_Theta_2nd_Ply"); CHKERRQ(ierr);
//    ierr = m_field_Macro.modify_finite_element_add_field_data("ELASTIC_4th_Ply","DISP_MACRO_r_Theta_3rd_Ply"); CHKERRQ(ierr);
//    ierr = m_field_Macro.modify_finite_element_add_field_data("ELASTIC_4th_Ply","DISP_MACRO_r_Theta_4th_Ply"); CHKERRQ(ierr);
  }

  //Define rows/cols and element data
  ierr = m_field_Macro.modify_finite_element_add_field_row("ELASTIC_FE_MACRO_REL","DISP_MACRO"); CHKERRQ(ierr);
  ierr = m_field_Macro.modify_finite_element_add_field_col("ELASTIC_FE_MACRO_REL","DISP_MACRO"); CHKERRQ(ierr);
  ierr = m_field_Macro.modify_finite_element_add_field_data("ELASTIC_FE_MACRO_REL","DISP_MACRO"); CHKERRQ(ierr);
  ierr = m_field_Macro.modify_finite_element_add_field_data("ELASTIC_FE_MACRO_REL","MESH_NODE_POSITIONS"); CHKERRQ(ierr);
  
  for(int ii=0; ii < nvars; ii++ ) {
    ostringstream ss_field;
    ss_field << "DISP_MACRO" << stochastic_fields[ii];
    ierr = m_field_Macro.modify_finite_element_add_field_data("ELASTIC_FE_MACRO_REL",ss_field.str().c_str()); CHKERRQ(ierr);
  }
//  ierr = m_field_Macro.modify_finite_element_add_field_data("ELASTIC_FE_MACRO_REL","DISP_MACRO_r_F"); CHKERRQ(ierr);
//  ierr = m_field_Macro.modify_finite_element_add_field_data("ELASTIC_FE_MACRO_REL","DISP_MACRO_r_Theta"); CHKERRQ(ierr);
//  ierr = m_field_Macro.modify_finite_element_add_field_data("ELASTIC_FE_MACRO_REL","DISP_MACRO_r_Theta_1st_Ply"); CHKERRQ(ierr);
//  ierr = m_field_Macro.modify_finite_element_add_field_data("ELASTIC_FE_MACRO_REL","DISP_MACRO_r_Theta_2nd_Ply"); CHKERRQ(ierr);
//  ierr = m_field_Macro.modify_finite_element_add_field_data("ELASTIC_FE_MACRO_REL","DISP_MACRO_r_Theta_3rd_Ply"); CHKERRQ(ierr);
//  ierr = m_field_Macro.modify_finite_element_add_field_data("ELASTIC_FE_MACRO_REL","DISP_MACRO_r_Theta_4th_Ply"); CHKERRQ(ierr);
  
  
  //define problems
  ierr = m_field_Macro.add_problem("ELASTIC_PROBLEM_MACRO"); CHKERRQ(ierr);
  
  //set finite elements for problem
  // First layer
  ierr = m_field_Macro.modify_problem_add_finite_element("ELASTIC_PROBLEM_MACRO","ELASTIC_1st_Ply"); CHKERRQ(ierr);
  // Second layer
  if (NO_Layers > 1) {
    ierr = m_field_Macro.modify_problem_add_finite_element("ELASTIC_PROBLEM_MACRO","ELASTIC_2nd_Ply"); CHKERRQ(ierr);
  }
  // Third layer
  if (NO_Layers > 2) {
    ierr = m_field_Macro.modify_problem_add_finite_element("ELASTIC_PROBLEM_MACRO","ELASTIC_3rd_Ply"); CHKERRQ(ierr);
  }
  // Fourth layer
  if (NO_Layers > 3) {
    ierr = m_field_Macro.modify_problem_add_finite_element("ELASTIC_PROBLEM_MACRO","ELASTIC_4th_Ply"); CHKERRQ(ierr);
  }
  
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
  
  for(int ii=0; ii < nvars; ii++ ) {
    ostringstream ss_field;
    ss_field << "DISP_MACRO" << stochastic_fields[ii];
    ierr = m_field_Macro.add_ents_to_field_by_TETs(0,ss_field.str().c_str()); CHKERRQ(ierr);
  }

  /*****************************************************************************
   *
   * Add finite elements entities
   *
   ****************************************************************************/
  ierr = m_field_Macro.add_ents_to_finite_element_by_TETs(TetsInBlock_1st_Ply, "ELASTIC_1st_Ply"); CHKERRQ(ierr);
  if (NO_Layers > 1) {
    ierr = m_field_Macro.add_ents_to_finite_element_by_TETs(TetsInBlock_2nd_Ply,"ELASTIC_2nd_Ply"); CHKERRQ(ierr);
  }
  if (NO_Layers > 2) {
    ierr = m_field_Macro.add_ents_to_finite_element_by_TETs(TetsInBlock_3rd_Ply,"ELASTIC_3rd_Ply"); CHKERRQ(ierr);
  }
  if (NO_Layers > 3) {
    ierr = m_field_Macro.add_ents_to_finite_element_by_TETs(TetsInBlock_4th_Ply,"ELASTIC_4th_Ply"); CHKERRQ(ierr);
  }
  ierr = m_field_Macro.add_ents_to_finite_element_by_TETs(TetsInBlock_Reliability,"ELASTIC_FE_MACRO_REL"); CHKERRQ(ierr);

  /*****************************************************************************
   *
   * Set applied order
   * See reference for detals:
   *   Ainsworth M. and Coyle J. (2003) Hierarchic finite element bases on
   *      unstructured tetrahedral meshes. IJNME, 58(14). pp.2103-2130.
   ****************************************************************************/
  ierr = m_field_Macro.set_field_order(0,MBTET,   "DISP_MACRO",order_Macro); CHKERRQ(ierr);
  ierr = m_field_Macro.set_field_order(0,MBTRI,   "DISP_MACRO",order_Macro); CHKERRQ(ierr);
  ierr = m_field_Macro.set_field_order(0,MBEDGE,  "DISP_MACRO",order_Macro); CHKERRQ(ierr);
  ierr = m_field_Macro.set_field_order(0,MBVERTEX,"DISP_MACRO",1); CHKERRQ(ierr);
  
  for(int ii=0; ii < nvars; ii++ ) {
    ostringstream ss_field;
    ss_field << "DISP_MACRO" << stochastic_fields[ii];
    ierr = m_field_Macro.set_field_order(0,MBTET,ss_field.str().c_str(),order_Macro); CHKERRQ(ierr);
    ierr = m_field_Macro.set_field_order(0,MBTRI,ss_field.str().c_str(),order_Macro); CHKERRQ(ierr);
    ierr = m_field_Macro.set_field_order(0,MBEDGE,ss_field.str().c_str(),order_Macro); CHKERRQ(ierr);
    ierr = m_field_Macro.set_field_order(0,MBVERTEX,ss_field.str().c_str(),1); CHKERRQ(ierr);
  }
  
  //
  ierr = m_field_Macro.set_field_order(0,MBTET,   "MESH_NODE_POSITIONS",2); CHKERRQ(ierr);
  ierr = m_field_Macro.set_field_order(0,MBTRI,   "MESH_NODE_POSITIONS",2); CHKERRQ(ierr);
  ierr = m_field_Macro.set_field_order(0,MBEDGE,  "MESH_NODE_POSITIONS",2); CHKERRQ(ierr);
  ierr = m_field_Macro.set_field_order(0,MBVERTEX,"MESH_NODE_POSITIONS",1); CHKERRQ(ierr);
 
 //***************************************************************************
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
  cout<<"///////////////////////////////////////////////////////////////////\n\n";
  
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
  ReliabOpt.step_code     = 1.0;       // 0: step size by Armijo rule, otherwise: given value is the step size
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
  probdata.PlyAngle    = readprobdata.PlyAngle;
  probdata.ExaminedPly = readprobdata.ExaminedLayer;
  int FailureCriterion; // 1: Tsai-Wu, 2: Tsai-Hill
  string NameOfFailureCriterion;
  FailureCriterion = readprobdata.FailureCriterion;
  
  //////////////////////////////////////////////////////////////////////////////
  //                                                                          //
  //        STEP 2: use probability transformation to obtain y in             //
  //                standard normal distribution space                        //
  //                y = norminv(F(x))                                         //
  //                                                                          //
  //////////////////////////////////////////////////////////////////////////////
  NatafTransformation my_nataf_transformation;
  ierr = my_nataf_transformation.ModCorrMat_Empirical(probdata.num_vars,
                                                      probdata.marg,
                                                      probdata.correlation,
                                                      probdata.mod_correlation); CHKERRQ(ierr);

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
  //  C. SOLVING FE EQUATION
  //
  // ===========================================================================

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
  ublas::matrix<double> StressGP(3,3);           StressGP.clear();
  ublas::vector<ublas::matrix<double> > TheStress(probdata.num_vars+1);
  
  double theta_angle;
  ublas::vector<double> PlyAngle_new;
  PlyAngle_new = probdata.PlyAngle;
  cout<<"\nAngle "<<probdata.PlyAngle<<endl;
  /*
   *  Start iteration
   */
  
  clock_t rel_t1, rel_t2;
  double rel_calc_time;
  rel_t1 = clock();
  
  
  //cout<<beta_flag<<endl;
  
  ofstream BetaFile;
  if (beta_flag == 1) {
    //BetaFile.open("~/Dropbox/DURACOMP_Cal/009_MoFEM/04_ReliabilityAnalysis/Result_Beta.txt",ofstream::out);
    BetaFile.open("Result_Beta.txt",ofstream::out);

  }
  
  LSF_Composite_Lamina TheLSF;
  Single_FE_Solver_Laminate Solve_FE_Problem;

  
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
    
    //
    // Transformation from u to x space
    //
    double detj = 1.0;
    x.resize(probdata.num_vars); x.clear();
    ierr = my_nataf_transformation.u_to_x(u,probdata.num_vars,probdata.marg,probdata.Lo,x,detj); CHKERRQ(ierr);
    cout<<"\nx value: "<<x<<endl;
    // ----------------
    // Update ply angle
    // -----------------
    for (unsigned i=1; i<=x.size();i++) {
      if (probdata.NameVars[i].compare(0,11,"orientation") == 0) {
        cout<<"\nAngle "<<x(i-1)<<endl;
        for (int j=0; j<probdata.PlyAngle.size(); j++) {
          cout << "The original ply angle "<<probdata.PlyAngle(j)<<"\t delta x: "<<x(i-1)<<endl;
          PlyAngle_new(j) = probdata.PlyAngle(j) + x(i-1);
          cout << "The modified ply angle "<<PlyAngle_new(j)<<endl;
        }
      }
      else if (probdata.NameVars[i].compare(0,6,"theta1") == 0) {
        PlyAngle_new(0) = x(i-1);
      }
      else if (probdata.NameVars[i].compare(0,6,"theta2") == 0) {
        PlyAngle_new(1) = x(i-1);
      }
      else if (probdata.NameVars[i].compare(0,6,"theta3") == 0) {
        PlyAngle_new(2) = x(i-1);
      }
      else if (probdata.NameVars[i].compare(0,6,"theta4") == 0) {
        PlyAngle_new(3) = x(i-1);
      }
    }
    //cout << "\n\nThe modified ply angle "<<PlyAngle_new<<endl;
    
    // Jacobian
    dudx.resize(probdata.num_vars,probdata.num_vars); dudx.clear();
    ierr = my_nataf_transformation.Jacobian_u_x(x,u,probdata.num_vars,probdata.marg,probdata.Lo,probdata.inv_Lo,dudx); CHKERRQ(ierr);
    
    inv_dudx.resize(probdata.num_vars,probdata.num_vars); inv_dudx.clear();
    inv_dudx = gjinverse(dudx,singular);
    
    //
    // Evaluate limit-state function and its gradient
    //
    cout<<"\n\nStart to run FE"<<endl;
    ierr = Solve_FE_Problem.Macro_FE_REL(m_field_Macro, nvars, nders,
                                          stochastic_fields, x, probdata.num_vars,
                                          probdata.NameVars,PlyAngle_new,
                                          probdata.ExaminedPly,
                                          NO_Layers,
                                          TheStress); CHKERRQ(ierr);

    StressGP.clear(); StressGP = TheStress(0); TheStress(0).clear();
    //cout<<"Stress at GP in xyz: "<<Calc_Stress.StressGP<<endl;
    cout<<"Stress at GP in 123: "<<StressGP<<endl;
    
    
    // Evaluate LSF and its gradient
    grad_g.resize(probdata.num_vars); grad_g.clear();
    //ierr = TheLSF->gfun(x,val_G,grad_g); CHKERRQ(ierr);
    
    switch (FailureCriterion) {
//      case 12013: {
//        //
//        // Maximum stress theory - fibre failure
//        //
//        NameOfFailureCriterion = "Maximum stress theory - Fibre failure";
//        ierr = TheLSF.gfun_ply_MS_LD(x,probdata.NameVars,probdata.MatStrength,
//                                     StressGP,
//                                     StressGP_r_Em,StressGP_r_NUm,
//                                     StressGP_r_NUp,StressGP_r_NUpz,
//                                     StressGP_r_Ep,StressGP_r_Ez,StressGP_r_Gzp,
//                                     StressGP_r_Ef,StressGP_r_NUf,
//                                     StressGP_r_F,StressGP_r_Theta,
//                                     StressGP_r_Theta_1,StressGP_r_Theta_2,
//                                     StressGP_r_Theta_3,StressGP_r_Theta_4,
//                                     val_G,grad_g); CHKERRQ(ierr);
//        break;
//      }
//      case 22013: {
//        //
//        // Maximum stress theory - matrix failure
//        //
//        NameOfFailureCriterion = "Maximum stress theory - Matrix failure";
//        ierr = TheLSF.gfun_ply_MS_TD(x,probdata.NameVars,probdata.MatStrength,
//                                     StressGP,
//                                     StressGP_r_Em,StressGP_r_NUm,
//                                     StressGP_r_NUp,StressGP_r_NUpz,
//                                     StressGP_r_Ep,StressGP_r_Ez,StressGP_r_Gzp,
//                                     StressGP_r_Ef,StressGP_r_NUf,
//                                     StressGP_r_F,StressGP_r_Theta,
//                                     StressGP_r_Theta_1,StressGP_r_Theta_2,
//                                     StressGP_r_Theta_3,StressGP_r_Theta_4,
//                                     val_G,grad_g); CHKERRQ(ierr);
//        break;
//      }
//      case 22014: {
//        //
//        // Maximum stress theory - shear failure
//        //
//        NameOfFailureCriterion = "Maximum stress theory - Shear failure";
//        ierr = TheLSF.gfun_ply_MS_Shear(x,probdata.NameVars,probdata.MatStrength,
//                                        StressGP,
//                                        StressGP_r_Em,StressGP_r_NUm,
//                                        StressGP_r_NUp,StressGP_r_NUpz,
//                                        StressGP_r_Ep,StressGP_r_Ez,StressGP_r_Gzp,
//                                        StressGP_r_Ef,StressGP_r_NUf,
//                                        StressGP_r_F,StressGP_r_Theta,
//                                        StressGP_r_Theta_1,StressGP_r_Theta_2,
//                                        StressGP_r_Theta_3,StressGP_r_Theta_4,
//                                        val_G,grad_g); CHKERRQ(ierr);
//        break;
//      }
//      case 13033: {
//        //
//        // Hashin failure theory - fibre failure
//        //
//        NameOfFailureCriterion = "Hashin failure theory - Fibre failure";
//        ierr = TheLSF.gfun_ply_HF(x,probdata.NameVars,probdata.MatStrength,
//                                  StressGP,
//                                  StressGP_r_Em,StressGP_r_NUm,
//                                  StressGP_r_NUp,StressGP_r_NUpz,
//                                  StressGP_r_Ep,StressGP_r_Ez,StressGP_r_Gzp,
//                                  StressGP_r_Ef,StressGP_r_NUf,
//                                  StressGP_r_F,StressGP_r_Theta,
//                                  StressGP_r_Theta_1,StressGP_r_Theta_2,
//                                  StressGP_r_Theta_3,StressGP_r_Theta_4,
//                                  val_G,grad_g); CHKERRQ(ierr);
//        
//        break;
//      }
//      case 23033: {
//        //
//        // Hashin failure theory - matrix failure
//        //
//        NameOfFailureCriterion = "Hashin failure theory - Matrix failure";
//        ierr = TheLSF.gfun_ply_HM(x,probdata.NameVars,probdata.MatStrength,
//                                  StressGP,
//                                  StressGP_r_Em,StressGP_r_NUm,
//                                  StressGP_r_NUp,StressGP_r_NUpz,
//                                  StressGP_r_Ep,StressGP_r_Ez,StressGP_r_Gzp,
//                                  StressGP_r_Ef,StressGP_r_NUf,
//                                  StressGP_r_F,StressGP_r_Theta,
//                                  StressGP_r_Theta_1,StressGP_r_Theta_2,
//                                  StressGP_r_Theta_3,StressGP_r_Theta_4,
//                                  val_G,grad_g); CHKERRQ(ierr);
//        break;
//      }
//      case 42050: {
//        //
//        // Tsai-Wu failure criteria
//        //
//        NameOfFailureCriterion = "Tsai-Wu - 2D stress state";
//        ierr = TheLSF.gfun_ply_Tsai_Wu_2D(x,probdata.NameVars,probdata.MatStrength,
//                                          StressGP,
//                                          StressGP_r_Em,StressGP_r_NUm,
//                                          StressGP_r_NUp,StressGP_r_NUpz,
//                                          StressGP_r_Ep,StressGP_r_Ez,StressGP_r_Gzp,
//                                          StressGP_r_Ef,StressGP_r_NUf,
//                                          StressGP_r_F,StressGP_r_Theta,
//                                          StressGP_r_Theta_1,StressGP_r_Theta_2,
//                                          StressGP_r_Theta_3,StressGP_r_Theta_4,
//                                          val_G,grad_g); CHKERRQ(ierr);
//        break;
//      }
      case 43050: {
        //
        // Tsai-Wu failure criteria
        //
        NameOfFailureCriterion = "Tsai-Wu - 3D stress state";
        ierr = TheLSF.gfun_ply_Tsai_Wu(x,probdata.NameVars,probdata.MatStrength,
                                       StressGP,TheStress,val_G,grad_g); CHKERRQ(ierr);
        break;
      }
//      case 44050: {
//        //
//        // Tsai-Wu failure criteria
//        //
//        NameOfFailureCriterion = "Tsai-Wu-Christensen - 3D stress state";
//        ierr = TheLSF.gfun_ply_Tsai_Wu_Christensen(x,probdata.NameVars,probdata.MatStrength,
//                                                   StressGP,
//                                                   StressGP_r_Em,StressGP_r_NUm,
//                                                   StressGP_r_NUp,StressGP_r_NUpz,
//                                                   StressGP_r_Ep,StressGP_r_Ez,StressGP_r_Gzp,
//                                                   StressGP_r_Ef,StressGP_r_NUf,
//                                                   StressGP_r_F,StressGP_r_Theta,
//                                                   StressGP_r_Theta_1,StressGP_r_Theta_2,
//                                                   StressGP_r_Theta_3,StressGP_r_Theta_4,
//                                                   val_G,grad_g); CHKERRQ(ierr);
//        break;
//      }
//      case 42060: {
//        //
//        // Tsai-Hill failure criteria
//        //
//        NameOfFailureCriterion = "Tsai-Hill - 2D stress-state";
//        ierr = TheLSF.gfun_ply_Tsai_Hill_2D(x,probdata.NameVars,probdata.MatStrength,
//                                            StressGP,
//                                            StressGP_r_Em,StressGP_r_NUm,
//                                            StressGP_r_NUp,StressGP_r_NUpz,
//                                            StressGP_r_Ep,StressGP_r_Ez,StressGP_r_Gzp,
//                                            StressGP_r_Ef,StressGP_r_NUf,
//                                            StressGP_r_F,StressGP_r_Theta,
//                                            StressGP_r_Theta_1,StressGP_r_Theta_2,
//                                            StressGP_r_Theta_3,StressGP_r_Theta_4,
//                                            val_G,grad_g); CHKERRQ(ierr);
//        break;
//      }
//      case 43060: {
//        //
//        // Tsai-Hill failure criteria
//        //
//        NameOfFailureCriterion = "Tsai-Hill - 3D stress-state";
//        ierr = TheLSF.gfun_ply_Tsai_Hill(x,probdata.NameVars,probdata.MatStrength,
//                                         StressGP,
//                                         StressGP_r_Em,StressGP_r_NUm,
//                                         StressGP_r_NUp,StressGP_r_NUpz,
//                                         StressGP_r_Ep,StressGP_r_Ez,StressGP_r_Gzp,
//                                         StressGP_r_Ef,StressGP_r_NUf,
//                                         StressGP_r_F,StressGP_r_Theta,
//                                         StressGP_r_Theta_1,StressGP_r_Theta_2,
//                                         StressGP_r_Theta_3,StressGP_r_Theta_4,
//                                         val_G,grad_g); CHKERRQ(ierr);
//        break;
//      }
//      case 13073: {
//        //
//        // Richard Christen: Fibre controlled failure
//        //
//        NameOfFailureCriterion = "Christensen - Fibre controlled failure";
//        ierr = TheLSF.gfun_ply_RCF(x,probdata.NameVars,probdata.MatStrength,
//                                   StressGP,
//                                   StressGP_r_Em,StressGP_r_NUm,
//                                   StressGP_r_NUp,StressGP_r_NUpz,
//                                   StressGP_r_Ep,StressGP_r_Ez,StressGP_r_Gzp,
//                                   StressGP_r_Ef,StressGP_r_NUf,
//                                   StressGP_r_F,StressGP_r_Theta,
//                                   StressGP_r_Theta_1,StressGP_r_Theta_2,
//                                   StressGP_r_Theta_3,StressGP_r_Theta_4,
//                                   val_G,grad_g); CHKERRQ(ierr);
//        break;
//      }
//      case 23073: {
//        //
//        // Richard Christen: Matrix controlled failure
//        //
//        NameOfFailureCriterion = "Christensen - Matrix controlled failure";
//        ierr = TheLSF.gfun_ply_RCM(x,probdata.NameVars,probdata.MatStrength,
//                                   StressGP,
//                                   StressGP_r_Em,StressGP_r_NUm,
//                                   StressGP_r_NUp,StressGP_r_NUpz,
//                                   StressGP_r_Ep,StressGP_r_Ez,StressGP_r_Gzp,
//                                   StressGP_r_Ef,StressGP_r_NUf,
//                                   StressGP_r_F,StressGP_r_Theta,
//                                   StressGP_r_Theta_1,StressGP_r_Theta_2,
//                                   StressGP_r_Theta_3,StressGP_r_Theta_4,
//                                   val_G,grad_g); CHKERRQ(ierr);
//        break;
//      }
//      case 42080: {
//        //
//        // Hoffman failure theory
//        //
//        NameOfFailureCriterion = "Hoffman failure theory - 2D stress states";
//        ierr = TheLSF.gfun_ply_Hoffman_2D(x,probdata.NameVars,probdata.MatStrength,
//                                          StressGP,
//                                          StressGP_r_Em,StressGP_r_NUm,
//                                          StressGP_r_NUp,StressGP_r_NUpz,
//                                          StressGP_r_Ep,StressGP_r_Ez,StressGP_r_Gzp,
//                                          StressGP_r_Ef,StressGP_r_NUf,
//                                          StressGP_r_F,StressGP_r_Theta,
//                                          StressGP_r_Theta_1,StressGP_r_Theta_2,
//                                          StressGP_r_Theta_3,StressGP_r_Theta_4,
//                                          val_G,grad_g); CHKERRQ(ierr);
//        break;
//      }
//      case 43080: {
//        //
//        // Hoffman failure theory
//        //
//        NameOfFailureCriterion = "Hoffman failure theory - 3D stress states";
//        ierr = TheLSF.gfun_ply_Hoffman(x,probdata.NameVars,probdata.MatStrength,
//                                       StressGP,
//                                       StressGP_r_Em,StressGP_r_NUm,
//                                       StressGP_r_NUp,StressGP_r_NUpz,
//                                       StressGP_r_Ep,StressGP_r_Ez,StressGP_r_Gzp,
//                                       StressGP_r_Ef,StressGP_r_NUf,
//                                       StressGP_r_F,StressGP_r_Theta,
//                                       StressGP_r_Theta_1,StressGP_r_Theta_2,
//                                       StressGP_r_Theta_3,StressGP_r_Theta_4,
//                                       val_G,grad_g); CHKERRQ(ierr);
//        break;
//      }
        //default: {}
    }

    grad_G.resize(probdata.num_vars); grad_G.clear();
    grad_G = prod(grad_g,inv_dudx);
    
    cout<<"LSF value is \t"<<val_G<<endl;
    cout<<"Gradient of LSF is \t"<<grad_g<<"\t"<<grad_G<<endl;
    
    //
    // Set scale parameter G0 and inform about structural response
    //
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
    
    //
    // Compute direction cosine alpha vector
    //
    alpha.resize(probdata.num_vars); alpha.clear();
    alpha = -grad_G/norm_2(grad_G);
    cout<<"Direction cosine "<<alpha<<endl;
    
    //
    // Check convergence
    //
    if (((abs(val_G/val_G0)<e1) && (norm_2(u-inner_prod(alpha,u)*alpha)<e2)) || (istep == istep_max)) {
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
    
    //
    // Take a step if convergence is not achieved
    //
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
      cout<<"\t"<<norm_2(u)<<"\t"<<inner_prod(alpha,u)<<endl;
      //cout<<"design point: "<<u<<endl;
      
      // Write reliability index in file
      if (beta_flag == 1) {
        BetaFile<<setprecision(15)<<inner_prod(alpha,u)<<"\t"<<val_G;
        for (int i=0;i<probdata.num_vars;i++) {
          BetaFile<<"\t"<<x(i);
        }
        BetaFile<<"\n";
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
  } while (conv_flag == 0);
  
  
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
  } else {
    cout<<"*  The number of iterations is: "<<istep<<".\n";
  }
  
  
  cout<<"*  Optimal reliability index (beta) is: "<<beta<<endl;
  cout<<"*  The failure criterion used is: "<<NameOfFailureCriterion<<endl;
  cout<<"*  Elapsed time is "<<total_time<<" seconds.\n";
  cout<<"*  The program finishes !!! \n*\n";
  cout<<"*************************************************"<<endl;

  ierr = PetscFinalize(); CHKERRQ(ierr);
  return 0;
  
}
