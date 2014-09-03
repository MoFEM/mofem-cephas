/* Copyright (C) 2014, 
 *   Xiao-Yi Zhou (xiaoyi.zhou AT newcastle.ac.uk)
 *   Zahur Ullah (Zahur.Ullah@glasgow.ac.uk)
 * --------------------------------------------------------------
 * This routine conducts finite element implementation for stochastic multiscale
 * homogenization problem under constant traction boundary condition by perturbation 
 * technique for single-phase isotropic material which has two independent 
 * material properties of Young's modulus and Poisson's ratio.
 *
 * HISTORY
 *
 * 2014.08.31 (first version)
 *
 * REFERENCES
 * 1. Kleiber M. and Hien T. D. (1992) The stochastic finite element method - 
 *      Basic perturbation technique and computer implementation. John Wiley & 
 *      Sons.
 * 2. Kaczamarczyk L., Pearce C. J. and Bicanic N. (2008) Scale transition and 
 *      enforcement of RVE boundary conditions in second-order computational
 *      homogenization. International Journal for Numerical Methods in 
 *      Engineering, 74(3) p506-522. 
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

#include "FieldInterface.hpp"
#include "FieldCore.hpp"
#include "FEMethod_UpLevelStudent.hpp"
#include "cholesky.hpp"
#include <petscksp.h>

#include "ElasticFEMethod.hpp"

#include "K_rYoungFEMethod.hpp"
#include "K_rPoissonFEMethod.hpp"

#include "K_rsYoungFEMethod.hpp"
#include "K_rYoungPoissonFEMethod.hpp"
#include "K_rsPoissonFEMethod.hpp"

#include "ElasticFE_RVELagrange_Traction.hpp"
#include "ElasticFE_RVELagrange_Homogenized_Stress_Traction.hpp"
#include "ElasticFE_RVELagrange_RigidBodyTranslation.hpp"
#include "ElasticFE_RVELagrange_RigidBodyRotation.hpp"
#include "RVEVolume.hpp"

#include "PostProcVertexMethod.hpp"
#include "PostProcDisplacementAndStrainOnRefindedMesh.hpp"

using namespace MoFEM;

ErrorCode rval;
PetscErrorCode ierr;

static char help[] = "...\n\n";

// =============================================================================
//
//  SETTING FOR OUTPUTS
//
// =============================================================================
PetscErrorCode write_soltion(FieldInterface &mField,const string out_file, const string out_ref_file) {
  PetscFunctionBegin;

  PostProcVertexMethod ent_method(mField.get_moab(),"DISPLACEMENT");
  ierr = mField.loop_dofs("STOCHASIC_PROBLEM","DISPLACEMENT",ROW,ent_method); CHKERRQ(ierr);
  // 1st-order partial derivative of displacement with respect to Young's modulus
  PostProcVertexMethod ent_method_r_E(mField.get_moab(),"DISP_r_E");
  ierr = mField.loop_dofs("DISP_r_E",ent_method_r_E); CHKERRQ(ierr);
  // 1st-order partial derivative of displacement with respect to Poisson's ratio
  PostProcVertexMethod ent_method_r_P(mField.get_moab(),"DISP_r_P");
  ierr = mField.loop_dofs("DISP_r_P",ent_method_r_P); CHKERRQ(ierr);

  // 2nd-order partial derivative of displacement with respect to Young's modulus
  PostProcVertexMethod ent_method_rs_EE(mField.get_moab(),"DISP_rs_EE");
  ierr = mField.loop_dofs("DISP_rs_EE",ent_method_rs_EE); CHKERRQ(ierr);
  // 2nd-order partial derivative of displacement with respect to Poisson's ratio
  PostProcVertexMethod ent_method_rs_PP(mField.get_moab(),"DISP_rs_PP");
  ierr = mField.loop_dofs("DISP_rs_PP",ent_method_rs_PP); CHKERRQ(ierr);
  // 2nd-order partial derivative of displacement with respect to Young's modulus
  // and Poisson's ratio 
  PostProcVertexMethod ent_method_rs_EP(mField.get_moab(),"DISP_rs_EP");
  ierr = mField.loop_dofs("DISP_rs_EP",ent_method_rs_EP); CHKERRQ(ierr);
  
  ParallelComm* pcomm = ParallelComm::get_pcomm(&mField.get_moab(),MYPCOMM_INDEX);
  
  if(pcomm->rank()==0) {
    EntityHandle out_meshset;
    rval = mField.get_moab().create_meshset(MESHSET_SET,out_meshset); CHKERR_PETSC(rval);
    ierr = mField.problem_get_FE("STOCHASIC_PROBLEM","ELASTIC",out_meshset); CHKERRQ(ierr);
    rval = mField.get_moab().write_file(out_file.c_str(),"VTK","",&out_meshset,1); CHKERR_PETSC(rval);
    rval = mField.get_moab().delete_entities(&out_meshset,1); CHKERR_PETSC(rval);
  }
  
  PetscFunctionReturn(0);
}


//==============================================================================
//
// MAIN ROUTINE
// 
//==============================================================================
int main(int argc, char *argv[]) {
  
  PetscInitialize(&argc,&argv,(char *)0,help);
  
  Core mb_instance;
  Interface& moab = mb_instance;
  int rank;
  MPI_Comm_rank(PETSC_COMM_WORLD,&rank);

  // ===========================================================================
  //
  //  I. READ MESH DATA AND FINITE ELEMENT ANALYSIS CONTROL PARAMETERS FROM FILE
  //
  // ===========================================================================  

  /*****************************************************************************
   *
   * Read parameters from line command
   *
   ****************************************************************************/ 
  PetscBool flg = PETSC_TRUE;
  char mesh_file_name[255];
  ierr = PetscOptionsGetString(PETSC_NULL,"-my_file",mesh_file_name,255,&flg); CHKERRQ(ierr);
  if(flg != PETSC_TRUE) {
    SETERRQ(PETSC_COMM_SELF,1,"*** ERROR -my_file (MESH FILE NEEDED)");
  }
  PetscInt order;
  ierr = PetscOptionsGetInt(PETSC_NULL,"-my_order",&order,&flg); CHKERRQ(ierr);
  if(flg != PETSC_TRUE) {
    order = 5;
  }
  
  /*****************************************************************************
   *
   * Read Applied strain on the RVE (vector of length 6) 
   * strain=[xx, yy, zz, xy, xz, zy]^T
   *
   ****************************************************************************/
  double myapplied_strain[6];
  int nmax=6;
  ierr = PetscOptionsGetRealArray(PETSC_NULL,"-myapplied_strain",myapplied_strain,&nmax,&flg); CHKERRQ(ierr);
  ublas::vector<FieldData> applied_strain;
  applied_strain.resize(6);
  cblas_dcopy(6, &myapplied_strain[0], 1, &applied_strain(0), 1);
  //    cout<<"applied_strain ="<<applied_strain<<endl;
  
  /*****************************************************************************
   *
   * Transfer mesh data to MOAB database
   *
   ****************************************************************************/
  const char *option;
  option = "";//"PARALLEL=BCAST;";//;DEBUG_IO";
  rval = moab.load_file(mesh_file_name, 0, option); CHKERR_PETSC(rval);
  ParallelComm* pcomm = ParallelComm::get_pcomm(&moab,MYPCOMM_INDEX);
  if(pcomm == NULL) pcomm =  new ParallelComm(&moab,PETSC_COMM_WORLD);
  
  //We need that for code profiling
  PetscLogDouble t1,t2;
  PetscLogDouble v1,v2;
  ierr = PetscTime(&v1); CHKERRQ(ierr);
  ierr = PetscGetCPUTime(&t1); CHKERRQ(ierr);
  
  //Create MoFEM (Joseph) database
  FieldCore core(moab);
  FieldInterface& mField = core;
  
  //ref meshset ref level 0
  ierr = mField.seed_ref_level_3D(0,0); CHKERRQ(ierr);
  
  /*****************************************************************************
   *
   * Select element into various mesh-set
   * meshset_level0: all element
   *
   ****************************************************************************/  
  // stl::bitset see for more details
  BitRefLevel bit_level0;
  bit_level0.set(0);
  EntityHandle meshset_level0;
  rval = moab.create_meshset(MESHSET_SET,meshset_level0); CHKERR_PETSC(rval);
  ierr = mField.seed_ref_level_3D(0,bit_level0); CHKERRQ(ierr);
  ierr = mField.get_entities_by_ref_level(bit_level0,BitRefLevel().set(),meshset_level0); CHKERRQ(ierr);

  // ===========================================================================
  //
  // II. DEFINE PROBLEM
  //
  // ===========================================================================
  
  //Fields
  int field_rank=3;
  ierr = mField.add_field("DISPLACEMENT",H1,field_rank); CHKERRQ(ierr);
  ierr = mField.add_field("Lagrange_mul_disp",NOFIELD,6); CHKERRQ(ierr);
  // Control 3 rigid body translations in x, y and z axis
  ierr = mField.add_field("Lagrange_mul_disp_rigid_trans",NOFIELD,3); CHKERRQ(ierr);   
  // Control 3 rigid body rotations about x, y and z axis
  ierr = mField.add_field("Lagrange_mul_disp_rigid_rotation",NOFIELD,3); CHKERRQ(ierr); 
  
  /*****************************************************************************
   *
   * Add stochastic field
   *
   ****************************************************************************/
  // 1st-order derivative of displacment
  ierr = mField.add_field("DISP_r_E",H1,field_rank,MF_ZERO); CHKERRQ(ierr);
  ierr = mField.add_field("DISP_r_P",H1,field_rank,MF_ZERO); CHKERRQ(ierr);
  // 2nd-order derivative of displacement
  ierr = mField.add_field("DISP_rs_EE",H1,field_rank,MF_ZERO); CHKERRQ(ierr);
  ierr = mField.add_field("DISP_rs_EP",H1,field_rank,MF_ZERO); CHKERRQ(ierr);
  ierr = mField.add_field("DISP_rs_PP",H1,field_rank,MF_ZERO); CHKERRQ(ierr);

  /*****************************************************************************
   *
   * Create finite element for the defined fields
   *
   ****************************************************************************/
  ierr = mField.add_finite_element("ELASTIC"); CHKERRQ(ierr);
  ierr = mField.add_finite_element("Lagrange_elem"); CHKERRQ(ierr);
  ierr = mField.add_finite_element("Lagrange_elem_rigid_trans"); CHKERRQ(ierr);
  ierr = mField.add_finite_element("Lagrange_elem_rigid_rotation"); CHKERRQ(ierr);

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
  //Define rows/cols and element data for Elastic element
  ierr = mField.modify_finite_element_add_field_row("ELASTIC","DISPLACEMENT"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_col("ELASTIC","DISPLACEMENT"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_data("ELASTIC","DISPLACEMENT"); CHKERRQ(ierr);
  // required for solving zeroth-order problem
  ierr = mField.modify_finite_element_add_field_data("ELASTIC","DISPLACEMENT"); CHKERRQ(ierr);
  // required for solving first problem
  ierr = mField.modify_finite_element_add_field_data("ELASTIC","DISP_r_E"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_data("ELASTIC","DISP_r_P"); CHKERRQ(ierr);
  // required for solving second-order problem
  ierr = mField.modify_finite_element_add_field_data("ELASTIC","DISP_rs_EE"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_data("ELASTIC","DISP_rs_EP"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_data("ELASTIC","DISP_rs_PP"); CHKERRQ(ierr);

  
  //Define rows/cols and element data for C and CT (for lagrange multipliers)
  //============================================================================
  //C row as Lagrange_mul_disp and col as DISPLACEMENT
  ierr = mField.modify_finite_element_add_field_row("Lagrange_elem","Lagrange_mul_disp"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_col("Lagrange_elem","DISPLACEMENT"); CHKERRQ(ierr);
  
  //CT col as Lagrange_mul_disp and row as DISPLACEMENT
  ierr = mField.modify_finite_element_add_field_col("Lagrange_elem","Lagrange_mul_disp"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_row("Lagrange_elem","DISPLACEMENT"); CHKERRQ(ierr);
  
  //As for stress we need both displacement and temprature (Lukasz)
  ierr = mField.modify_finite_element_add_field_data("Lagrange_elem","Lagrange_mul_disp"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_data("Lagrange_elem","DISPLACEMENT"); CHKERRQ(ierr);
  
  //============================================================================ 
  
  //Define rows/cols and element data for C1 and C1T (for lagrange multipliers to contol the rigid body translations)
  //============================================================================
  //C1 row as Lagrange_elem_rigid_trans and col as DISPLACEMENT
  ierr = mField.modify_finite_element_add_field_row("Lagrange_elem_rigid_trans","Lagrange_mul_disp_rigid_trans"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_col("Lagrange_elem_rigid_trans","DISPLACEMENT"); CHKERRQ(ierr);
  
  //C1T col as Lagrange_elem_rigid_trans and row as DISPLACEMENT
  ierr = mField.modify_finite_element_add_field_col("Lagrange_elem_rigid_trans","Lagrange_mul_disp_rigid_trans"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_row("Lagrange_elem_rigid_trans","DISPLACEMENT"); CHKERRQ(ierr);
  
  //As for stress we need both displacement and temprature (Lukasz)
  ierr = mField.modify_finite_element_add_field_data("Lagrange_elem_rigid_trans","Lagrange_mul_disp_rigid_trans"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_data("Lagrange_elem_rigid_trans","DISPLACEMENT"); CHKERRQ(ierr);
  //============================================================================
  
  
  //Define rows/cols and element data for C2 and C2T (for lagrange multipliers to contol the rigid body rotations)
  //============================================================================
  //C2 row as Lagrange_elem_rigid_trans and col as DISPLACEMENT
  ierr = mField.modify_finite_element_add_field_row("Lagrange_elem_rigid_rotation","Lagrange_mul_disp_rigid_rotation"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_col("Lagrange_elem_rigid_rotation","DISPLACEMENT"); CHKERRQ(ierr);
  
  //C2T col as Lagrange_elem_rigid_trans and row as DISPLACEMENT
  ierr = mField.modify_finite_element_add_field_col("Lagrange_elem_rigid_rotation","Lagrange_mul_disp_rigid_rotation"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_row("Lagrange_elem_rigid_rotation","DISPLACEMENT"); CHKERRQ(ierr);
  
  //As for stress we need both displacement and temprature (Lukasz)
  ierr = mField.modify_finite_element_add_field_data("Lagrange_elem_rigid_rotation","Lagrange_mul_disp_rigid_rotation"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_data("Lagrange_elem_rigid_rotation","DISPLACEMENT"); CHKERRQ(ierr);
  //============================================================================ 
  
  //define problems
  ierr = mField.add_problem("STOCHASIC_PROBLEM"); CHKERRQ(ierr);
  
  
  //set finite elements for problem
  ierr = mField.modify_problem_add_finite_element("STOCHASIC_PROBLEM","ELASTIC"); CHKERRQ(ierr);
  ierr = mField.modify_problem_add_finite_element("STOCHASIC_PROBLEM","Lagrange_elem"); CHKERRQ(ierr);
  ierr = mField.modify_problem_add_finite_element("STOCHASIC_PROBLEM","Lagrange_elem_rigid_trans"); CHKERRQ(ierr);
  ierr = mField.modify_problem_add_finite_element("STOCHASIC_PROBLEM","Lagrange_elem_rigid_rotation"); CHKERRQ(ierr);
  
  
  //set refinment level for problem
  ierr = mField.modify_problem_ref_level_add_bit("STOCHASIC_PROBLEM",bit_level0); CHKERRQ(ierr);
  
  // ===========================================================================
  //
  // III. DECLARE PROBLEM
  //
  // ===========================================================================

  /*****************************************************************************
   *
   * Add entitities (by tets) to the field
   *
   ****************************************************************************/
  // Zeroth order
  ierr = mField.add_ents_to_field_by_TETs(0,"DISPLACEMENT",2); CHKERRQ(ierr);
  // First order
  ierr = mField.add_ents_to_field_by_TETs(0,"DISP_r_E"); CHKERRQ(ierr);
  ierr = mField.add_ents_to_field_by_TETs(0,"DISP_r_P"); CHKERRQ(ierr);
  // Second order
  ierr = mField.add_ents_to_field_by_TETs(0,"DISP_rs_EE"); CHKERRQ(ierr);
  ierr = mField.add_ents_to_field_by_TETs(0,"DISP_rs_EP"); CHKERRQ(ierr);
  ierr = mField.add_ents_to_field_by_TETs(0,"DISP_rs_PP"); CHKERRQ(ierr);

  /*****************************************************************************
   *
   * Add finite elements entities
   *
   ****************************************************************************/ 
  ierr = mField.add_ents_to_finite_element_EntType_by_bit_ref(bit_level0,"ELASTIC",MBTET); CHKERRQ(ierr);
  Range SurfacesFaces;
  //Add finite element to lagrange element for rigid body translation
  ierr = mField.get_Cubit_msId_entities_by_dimension(103,SIDESET,2,SurfacesFaces,true); CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"number of SideSet 103 = %d\n",SurfacesFaces.size()); CHKERRQ(ierr);
  
  //to create meshset from range
  EntityHandle BoundFacesMeshset;
  rval = moab.create_meshset(MESHSET_SET,BoundFacesMeshset); CHKERR_PETSC(rval);
	rval = moab.add_entities(BoundFacesMeshset,SurfacesFaces); CHKERR_PETSC(rval);
  ierr = mField.seed_ref_level_MESHSET(BoundFacesMeshset,BitRefLevel().set()); CHKERRQ(ierr);

  ierr = mField.add_ents_to_finite_element_by_TRIs(SurfacesFaces,"Lagrange_elem"); CHKERRQ(ierr);
  ierr = mField.add_ents_to_finite_element_by_TRIs(SurfacesFaces,"Lagrange_elem_rigid_trans"); CHKERRQ(ierr);
  ierr = mField.add_ents_to_finite_element_by_TRIs(SurfacesFaces,"Lagrange_elem_rigid_rotation"); CHKERRQ(ierr);
  
  /*****************************************************************************
   *
   * Set applied order
   * See reference for detals:
   *   Ainsworth M. and Coyle J. (2003) Hierarchic finite element bases on 
          unstructured tetrahedral meshes. IJNME, 58(14). pp.2103-2130.
   ****************************************************************************/  
  //see Hierarchic Finite Element Bases on Unstructured Tetrahedral Meshes (Mark Ainsworth & Joe Coyle)
  //int order = 5;
  ierr = mField.set_field_order(0,MBTET,"DISPLACEMENT",order); CHKERRQ(ierr);
  ierr = mField.set_field_order(0,MBTRI,"DISPLACEMENT",order); CHKERRQ(ierr);
  ierr = mField.set_field_order(0,MBEDGE,"DISPLACEMENT",order); CHKERRQ(ierr);
  ierr = mField.set_field_order(0,MBVERTEX,"DISPLACEMENT",1); CHKERRQ(ierr);
  
  int order_st=order;  
  // 1st order-1
  ierr = mField.set_field_order(0,MBTET,"DISP_r_E",order_st); CHKERRQ(ierr);
  ierr = mField.set_field_order(0,MBTRI,"DISP_r_E",order_st); CHKERRQ(ierr);
  ierr = mField.set_field_order(0,MBEDGE,"DISP_r_E",order_st); CHKERRQ(ierr);
  ierr = mField.set_field_order(0,MBVERTEX,"DISP_r_E",1); CHKERRQ(ierr);

  // 1st order-2
  ierr = mField.set_field_order(0,MBTET,"DISP_r_P",order_st); CHKERRQ(ierr);
  ierr = mField.set_field_order(0,MBTRI,"DISP_r_P",order_st); CHKERRQ(ierr);
  ierr = mField.set_field_order(0,MBEDGE,"DISP_r_P",order_st); CHKERRQ(ierr);
  ierr = mField.set_field_order(0,MBVERTEX,"DISP_r_P",1); CHKERRQ(ierr);
  
  // 2nd order-1
  ierr = mField.set_field_order(0,MBTET,"DISP_rs_EE",order_st); CHKERRQ(ierr);
  ierr = mField.set_field_order(0,MBTRI,"DISP_rs_EE",order_st); CHKERRQ(ierr);
  ierr = mField.set_field_order(0,MBEDGE,"DISP_rs_EE",order_st); CHKERRQ(ierr);
  ierr = mField.set_field_order(0,MBVERTEX,"DISP_rs_EE",1); CHKERRQ(ierr);

  // 2nd order-2
  ierr = mField.set_field_order(0,MBTET,"DISP_rs_EP",order_st); CHKERRQ(ierr);
  ierr = mField.set_field_order(0,MBTRI,"DISP_rs_EP",order_st); CHKERRQ(ierr);
  ierr = mField.set_field_order(0,MBEDGE,"DISP_rs_EP",order_st); CHKERRQ(ierr);
  ierr = mField.set_field_order(0,MBVERTEX,"DISP_rs_EP",1); CHKERRQ(ierr);

  // 2nd order-3
  ierr = mField.set_field_order(0,MBTET,"DISP_rs_PP",order_st); CHKERRQ(ierr);
  ierr = mField.set_field_order(0,MBTRI,"DISP_rs_PP",order_st); CHKERRQ(ierr);
  ierr = mField.set_field_order(0,MBEDGE,"DISP_rs_PP",order_st); CHKERRQ(ierr);
  ierr = mField.set_field_order(0,MBVERTEX,"DISP_rs_PP",1); CHKERRQ(ierr);


  // ===========================================================================
  //
  //  IV. BUILD DATABASE
  //
  // ===========================================================================
  
  //build field
  ierr = mField.build_fields(); CHKERRQ(ierr);
  
  //build finite elemnts
  ierr = mField.build_finite_elements(); CHKERRQ(ierr);
  
  //build adjacencies
  ierr = mField.build_adjacencies(bit_level0); CHKERRQ(ierr);
  
  //build problem
  ierr = mField.build_problems(); CHKERRQ(ierr);
  
  
  // ===========================================================================
  //
  //  V. MESH PARTITION
  //
  // ===========================================================================
  
  //partition
  ierr = mField.partition_problem("STOCHASIC_PROBLEM"); CHKERRQ(ierr);
  ierr = mField.partition_finite_elements("STOCHASIC_PROBLEM"); CHKERRQ(ierr);
  //what are ghost nodes, see Petsc Manual
  ierr = mField.partition_ghost_dofs("STOCHASIC_PROBLEM"); CHKERRQ(ierr);
  
  //print bcs
  ierr = mField.print_cubit_displacement_set(); CHKERRQ(ierr);
  ierr = mField.print_cubit_force_set(); CHKERRQ(ierr);
  //print block sets with materials
  ierr = mField.print_cubit_materials_set(); CHKERRQ(ierr);
  
  //    for(_IT_GET_DOFS_FIELD_BY_NAME_AND_TYPE_FOR_LOOP_(mField,"Lagrange_mul_disp",MBEDGE,dof)) {
  //        cerr << *dof << endl;
  //
  //    }


  // ===========================================================================
  //
  //  VI. SOLUTION PHASE
  //
  // =========================================================================== 

  /*****************************************************************************
   *
   *  0. PREPARATION FOR PROCESSING SOLVE
   *
   ****************************************************************************/  
  //create verctor and matrices for the stochastic problem inlcuiding zero, first and second order
  Vec F,dF,ddF;
  Vec D;                                   // 0th solution of nodal displacement
  Vec dD_E,dD_P;                           // 1st solution of nodal displacement 
  Vec ddD_EE,ddD_EP,ddD_PP;                // 2nd solution of nodal displacement

  ierr = mField.VecCreateGhost("STOCHASIC_PROBLEM",ROW,&F); CHKERRQ(ierr);
  ierr = mField.VecCreateGhost("STOCHASIC_PROBLEM",ROW,&dF); CHKERRQ(ierr);
  ierr = mField.VecCreateGhost("STOCHASIC_PROBLEM",ROW,&ddF); CHKERRQ(ierr);
  // zeroth-order problem  
  ierr = mField.VecCreateGhost("STOCHASIC_PROBLEM",COL,&D); CHKERRQ(ierr);
  // first-order problem
  ierr = mField.VecCreateGhost("STOCHASIC_PROBLEM",COL,&dD_E); CHKERRQ(ierr);
  ierr = mField.VecCreateGhost("STOCHASIC_PROBLEM",COL,&dD_P); CHKERRQ(ierr);
  // second-order problem
  ierr = mField.VecCreateGhost("STOCHASIC_PROBLEM",COL,&ddD_EE); CHKERRQ(ierr);
  ierr = mField.VecCreateGhost("STOCHASIC_PROBLEM",COL,&ddD_EP); CHKERRQ(ierr);
  ierr = mField.VecCreateGhost("STOCHASIC_PROBLEM",COL,&ddD_PP); CHKERRQ(ierr);

  /*****************************************************************************
   *
   *  1. Assembling global stiffness matrix K 
   *     and external force vector F
   ****************************************************************************/    
  Mat Aij;
  ierr = mField.MatCreateMPIAIJWithArrays("STOCHASIC_PROBLEM",&Aij); CHKERRQ(ierr);
  
  struct MyElasticFEMethod: public ElasticFEMethod {
    MyElasticFEMethod(FieldInterface& _mField,
                      Mat &_Aij,Vec &_D,Vec& _F,double _lambda,double _mu):
    ElasticFEMethod(_mField,_Aij,_D,_F,_lambda,_mu) {};
    
    PetscErrorCode Fint(Vec F_int) {
      PetscFunctionBegin;
      ierr = ElasticFEMethod::Fint(); CHKERRQ(ierr);
      for(int rr = 0;rr<row_mat;rr++) {
        if(RowGlob[rr].size()!=f_int[rr].size()) SETERRQ(PETSC_COMM_SELF,1,"data inconsistency");
        if(RowGlob[rr].size()==0) continue;
        f_int[rr] *= -1; //This is not SNES we solve K*D = -RES
        ierr = VecSetValues(F_int,RowGlob[rr].size(),&(RowGlob[rr])[0],&(f_int[rr].data()[0]),ADD_VALUES); CHKERRQ(ierr);
      }
      PetscFunctionReturn(0);
    }
  };
  
  //Assemble F and Aij
  const double young_modulus = 1;
  const double poisson_ratio = 0.0;
  MyElasticFEMethod my_fe(mField,Aij,D,F,LAMBDA(young_modulus,poisson_ratio),MU(young_modulus,poisson_ratio));
  
  ElasticFE_RVELagrange_Traction MyFE_RVELagrange(mField,Aij,D,F,applied_strain,"DISPLACEMENT","Lagrange_mul_disp",field_rank);
  ElasticFE_RVELagrange_RigidBodyTranslation MyFE_RVELagrangeRigidBodyTrans(mField,Aij,D,F,applied_strain,"DISPLACEMENT","Lagrange_mul_disp",field_rank,"Lagrange_mul_disp_rigid_trans");
  ElasticFE_RVELagrange_RigidBodyRotation MyFE_RVELagrangeRigidBodyRotation(mField,Aij,D,F,applied_strain,"DISPLACEMENT","Lagrange_mul_disp",field_rank);
  
  
  ierr = VecZeroEntries(F); CHKERRQ(ierr);
  ierr = VecGhostUpdateBegin(F,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = VecGhostUpdateEnd(F,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = MatZeroEntries(Aij); CHKERRQ(ierr);
  
  ierr = mField.loop_finite_elements("STOCHASIC_PROBLEM","Lagrange_elem",MyFE_RVELagrange);  CHKERRQ(ierr);
  ierr = mField.loop_finite_elements("STOCHASIC_PROBLEM","ELASTIC",my_fe);  CHKERRQ(ierr);
  ierr = mField.loop_finite_elements("STOCHASIC_PROBLEM","Lagrange_elem_rigid_trans",MyFE_RVELagrangeRigidBodyTrans);  CHKERRQ(ierr);
  ierr = mField.loop_finite_elements("STOCHASIC_PROBLEM","Lagrange_elem_rigid_rotation",MyFE_RVELagrangeRigidBodyRotation);  CHKERRQ(ierr);
  
  ierr = VecGhostUpdateBegin(F,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
  ierr = VecGhostUpdateEnd(F,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
  ierr = VecAssemblyBegin(F); CHKERRQ(ierr);
  ierr = VecAssemblyEnd(F); CHKERRQ(ierr);
  ierr = MatAssemblyBegin(Aij,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
  ierr = MatAssemblyEnd(Aij,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
  
  PetscSynchronizedFlush(PETSC_COMM_WORLD);
  //ierr = VecView(F,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);
  //ierr = MatView(Aij,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);


////Matrix View
//MatView(Aij,PETSC_VIEWER_DRAW_WORLD);//PETSC_VIEWER_STDOUT_WORLD);
//std::string wait;
//std::cin >> wait;

  /*****************************************************************************
   *
   *  2. Get the volume of RVE 
   *
   ****************************************************************************/ 
  double RVE_volume;    RVE_volume=0.0;  //RVE volume for full RVE We need this for stress calculation
  Vec RVE_volume_Vec;
  ierr = VecCreateMPI(PETSC_COMM_WORLD, 1, pcomm->size(), &RVE_volume_Vec);  CHKERRQ(ierr);
  ierr = VecZeroEntries(RVE_volume_Vec); CHKERRQ(ierr);
  
  RVEVolume MyRVEVol(mField,Aij,D,F,LAMBDA(young_modulus,poisson_ratio),MU(young_modulus,poisson_ratio), RVE_volume_Vec);
  ierr = mField.loop_finite_elements("STOCHASIC_PROBLEM","ELASTIC",MyRVEVol);  CHKERRQ(ierr);
  //    ierr = VecView(RVE_volume_Vec,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);
  ierr = VecSum(RVE_volume_Vec, &RVE_volume);  CHKERRQ(ierr);
    cout<<"Final RVE_volume = "<< RVE_volume <<endl;
  
  /*****************************************************************************
   *
   *  3. SOLVE THE ZEROTH-ORDER FINITE ELEMENT EQUILIBRIUM EQUATION
   *     [K][U] = [F]
   *
   ****************************************************************************/

  //----------------------------------------------------------------------------
  // 3.1 Solving the equation to get nodal displacement
  //----------------------------------------------------------------------------
  KSP solver;
  ierr = KSPCreate(PETSC_COMM_WORLD,&solver); CHKERRQ(ierr);
  ierr = KSPSetOperators(solver,Aij,Aij,SAME_NONZERO_PATTERN); CHKERRQ(ierr);
  ierr = KSPSetFromOptions(solver); CHKERRQ(ierr);
  ierr = KSPSetUp(solver); CHKERRQ(ierr);
  
  ierr = KSPSolve(solver,F,D); CHKERRQ(ierr);
  ierr = VecGhostUpdateBegin(D,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = VecGhostUpdateEnd(D,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  
  //Save data on mesh
  ierr = mField.set_global_VecCreateGhost("STOCHASIC_PROBLEM",ROW,D,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
//  ierr = VecView(D,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);
  
  //----------------------------------------------------------------------------
  // 3.2 Calculating zeroth-order homogenized stress using volume averaging theorem
  //----------------------------------------------------------------------------
  //create a vector for 6 components of homogenized stress
  Vec Stress_Homo;
  ierr = VecCreateMPI(PETSC_COMM_WORLD, 6, 6*pcomm->size(), &Stress_Homo);  CHKERRQ(ierr);
  ierr = VecZeroEntries(Stress_Homo); CHKERRQ(ierr);
  
  ////    ierr = VecView(D,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);
  ElasticFE_RVELagrange_Homogenized_Stress_Traction MyFE_RVEHomoStressTraction(mField,Aij,D,F,&RVE_volume,applied_strain, Stress_Homo,"DISPLACEMENT","Lagrange_mul_disp",field_rank);
  ierr = mField.loop_finite_elements("STOCHASIC_PROBLEM","Lagrange_elem",MyFE_RVEHomoStressTraction);  CHKERRQ(ierr);
  
//  if(pcomm->rank()) cout<< " Stress_Homo =  "<<endl;
//  ierr = VecView(Stress_Homo,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);
  
  if(pcomm->rank()==0){
    PetscScalar    *avec;
    VecGetArray(Stress_Homo, &avec);
    
    cout<< "\nStress_Homo = \n\n";
    for(int ii=0; ii<6; ii++){
      cout.precision(15);
      cout <<*avec<<endl;
      avec++;
    }
  }
  cout<< "\n\n";
  
  /*****************************************************************************
   *
   *  4. SOLVE THE FIRST-ORDER FINITE ELEMENT EQUILIBRIUM EQUATION
   *     [K][U_r] = -[K_r][U}  due to both young and poisson
   *
   ****************************************************************************/

  // ==================================
  // 4.1 due to Young's modulus
  // ===================================

  //----------------------------------------------------------------------------
  // a. Solving the equation to get nodal displacement
  //----------------------------------------------------------------------------
  ierr = VecZeroEntries(dF); CHKERRQ(ierr);
  ierr = VecGhostUpdateBegin(dF,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = VecGhostUpdateEnd(dF,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  
  K_rYoungFEMethod my_fe_k_r_E(mField,Aij,D,dF,LAMBDA(young_modulus,poisson_ratio),MU(young_modulus,poisson_ratio));
  ierr = mField.loop_finite_elements("STOCHASIC_PROBLEM","ELASTIC",my_fe_k_r_E);  CHKERRQ(ierr);
  
  ierr = VecGhostUpdateBegin(dF,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
  ierr = VecGhostUpdateEnd(dF,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
  ierr = VecAssemblyBegin(dF); CHKERRQ(ierr);
  ierr = VecAssemblyEnd(dF); CHKERRQ(ierr);
//  ierr = VecView(dF,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);
  
  ierr = KSPSolve(solver,dF,dD_E); CHKERRQ(ierr);
  ierr = VecGhostUpdateBegin(dD_E,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = VecGhostUpdateEnd(dD_E,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);

  ierr = mField.set_other_global_VecCreateGhost("STOCHASIC_PROBLEM","DISPLACEMENT","DISP_r_E",ROW,dD_E,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
  ierr = mField.set_other_global_VecCreateGhost("STOCHASIC_PROBLEM","Lagrange_mul_disp","Lagrange_mul_disp",ROW,dD_E,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
//  ierr = VecView(dD,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);

  
  //----------------------------------------------------------------------------
  // b. Calculating first-order homogenized stress
  //---------------------------------------------------------------------------- 
  Vec Stress_Homo_r_E;
  ierr = VecCreateMPI(PETSC_COMM_WORLD, 6, 6*pcomm->size(), &Stress_Homo_r_E);  CHKERRQ(ierr);
  ierr = VecZeroEntries(Stress_Homo_r_E); CHKERRQ(ierr);
  
  //    ierr = VecView(D,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);
  ElasticFE_RVELagrange_Homogenized_Stress_Traction MyFE_RVEHomoStressTraction_r_E(mField,Aij,D,F,&RVE_volume, applied_strain, Stress_Homo_r_E,"DISP_r_E","Lagrange_mul_disp",field_rank);
  ierr = mField.loop_finite_elements("STOCHASIC_PROBLEM","Lagrange_elem",MyFE_RVEHomoStressTraction_r_E);  CHKERRQ(ierr);
  
  if(pcomm->rank()==0){
    PetscScalar    *avec_r_E;
    VecGetArray(Stress_Homo_r_E, &avec_r_E);
    
    cout<< "\nStress_Homo_r_E = \n\n";
    for(int ii=0; ii<6; ii++){
      cout <<*avec_r_E<<endl; ;
      avec_r_E++;
    }
  }
  cout<< "\n\n";

  // ==================================
  // 4.2 due to Poisson's ratio
  // ===================================

  //----------------------------------------------------------------------------
  // a. Solving the equation to get nodal displacement
  //----------------------------------------------------------------------------
  ierr = VecZeroEntries(dF); CHKERRQ(ierr);
  ierr = VecGhostUpdateBegin(dF,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = VecGhostUpdateEnd(dF,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  
  K_rPoissonFEMethod my_fe_k_r_P(mField,Aij,D,dF,LAMBDA(young_modulus,poisson_ratio),MU(young_modulus,poisson_ratio));
  ierr = mField.loop_finite_elements("STOCHASIC_PROBLEM","ELASTIC",my_fe_k_r_P);  CHKERRQ(ierr);
  
  ierr = VecGhostUpdateBegin(dF,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
  ierr = VecGhostUpdateEnd(dF,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
  ierr = VecAssemblyBegin(dF); CHKERRQ(ierr);
  ierr = VecAssemblyEnd(dF); CHKERRQ(ierr);
//  ierr = VecView(dF,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);
  
  ierr = KSPSolve(solver,dF,dD_P); CHKERRQ(ierr);
  ierr = VecGhostUpdateBegin(dD_P,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = VecGhostUpdateEnd(dD_P,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);

  ierr = mField.set_other_global_VecCreateGhost("STOCHASIC_PROBLEM","DISPLACEMENT","DISP_r_P",ROW,dD_P,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
  ierr = mField.set_other_global_VecCreateGhost("STOCHASIC_PROBLEM","Lagrange_mul_disp","Lagrange_mul_disp",ROW,dD_P,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
//  ierr = VecView(dD,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);

  
  //----------------------------------------------------------------------------
  // b. Calculating first-order homogenized stress
  //---------------------------------------------------------------------------- 
  Vec Stress_Homo_r_P;
  ierr = VecCreateMPI(PETSC_COMM_WORLD, 6, 6*pcomm->size(), &Stress_Homo_r_P);  CHKERRQ(ierr);
  ierr = VecZeroEntries(Stress_Homo_r_P); CHKERRQ(ierr);
  
  //    ierr = VecView(D,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);
  ElasticFE_RVELagrange_Homogenized_Stress_Traction MyFE_RVEHomoStressTraction_r_P(mField,Aij,D,F,&RVE_volume, applied_strain, Stress_Homo_r_P,"DISP_r_P","Lagrange_mul_disp",field_rank);
  ierr = mField.loop_finite_elements("STOCHASIC_PROBLEM","Lagrange_elem",MyFE_RVEHomoStressTraction_r_P);  CHKERRQ(ierr);
  
  if(pcomm->rank()==0){
    PetscScalar    *avec_r_P;
    VecGetArray(Stress_Homo_r_P, &avec_r_P);
    
    cout<< "\nStress_Homo_r_P = \n\n";
    for(int ii=0; ii<6; ii++){
      cout <<*avec_r_P<<endl; ;
      avec_r_P++;
    }
  }
  cout<< "\n\n";

  /*****************************************************************************
   *
   *  5. SOLVE THE SECOND-ORDER FINITE ELEMENT EQUILIBRIUM EQUATION
   *     [K][U_rs] = -[K_rs][U]-2[K_r][U_s]
   *
   ****************************************************************************/

  // ==================================
  // 5.1 due to Young's modulus
  // ===================================
  //----------------------------------------------------------------------------
  // a. Solving the equation to get nodal displacement
  //----------------------------------------------------------------------------
  ierr = VecZeroEntries(ddF); CHKERRQ(ierr);
  ierr = VecGhostUpdateBegin(ddF,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = VecGhostUpdateEnd(ddF,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  
  K_rsYoungFEMethod my_fe_k_rs_EE(mField,Aij,D,ddF,LAMBDA(young_modulus,poisson_ratio),MU(young_modulus,poisson_ratio),"DISP_r_E");
  ierr = mField.loop_finite_elements("STOCHASIC_PROBLEM","ELASTIC",my_fe_k_rs_EE);  CHKERRQ(ierr);

  ierr = VecGhostUpdateBegin(ddF,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
  ierr = VecGhostUpdateEnd(ddF,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
  ierr = VecAssemblyBegin(ddF); CHKERRQ(ierr);
  ierr = VecAssemblyEnd(ddF); CHKERRQ(ierr);
//  ierr = VecView(ddF,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);
  
  ierr = KSPSolve(solver,ddF,ddD_EE); CHKERRQ(ierr);
  ierr = VecGhostUpdateBegin(ddD_EE,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = VecGhostUpdateEnd(ddD_EE,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  
  ierr = mField.set_other_global_VecCreateGhost("STOCHASIC_PROBLEM","DISPLACEMENT","DISP_rs_EE",ROW,ddD_EE,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
  ierr = mField.set_other_global_VecCreateGhost("STOCHASIC_PROBLEM","Lagrange_mul_disp","Lagrange_mul_disp",ROW,ddD_EE,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
//  ierr = VecView(dD,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);
  
  
  //----------------------------------------------------------------------------
  // b. Calculating second-order homogenized stress
  //----------------------------------------------------------------------------
  Vec Stress_Homo_rs_EE;
  ierr = VecCreateMPI(PETSC_COMM_WORLD, 6, 6*pcomm->size(), &Stress_Homo_rs_EE);  CHKERRQ(ierr);
  ierr = VecZeroEntries(Stress_Homo_rs_EE); CHKERRQ(ierr);
  
  //    ierr = VecView(D,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);
  ElasticFE_RVELagrange_Homogenized_Stress_Traction MyFE_RVEHomoStressTraction_rs_EE(mField,Aij,D,F,&RVE_volume, applied_strain, Stress_Homo_rs_EE,"DISP_rs_EE","Lagrange_mul_disp",field_rank);
  ierr = mField.loop_finite_elements("STOCHASIC_PROBLEM","Lagrange_elem",MyFE_RVEHomoStressTraction_rs_EE);  CHKERRQ(ierr);
  
  if(pcomm->rank()==0){
    PetscScalar    *avec_rs_EE;
    VecGetArray(Stress_Homo_rs_EE, &avec_rs_EE);
    
    cout<< "\nStress_Homo_rs_EE = \n\n";
    for(int ii=0; ii<6; ii++){
      cout <<*avec_rs_EE<<endl; ;
      avec_rs_EE++;
    }
  }
  cout<< "\n\n";

  // ==============================================
  // 5.2 due to Young's modulus and Poisson's ratio
  // ==============================================
  //----------------------------------------------------------------------------
  // a. Solving the equation to get nodal displacement
  //----------------------------------------------------------------------------
  ierr = VecZeroEntries(ddF); CHKERRQ(ierr);
  ierr = VecGhostUpdateBegin(ddF,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = VecGhostUpdateEnd(ddF,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  
  K_rYoungPoissonFEMethod my_fe_k_rs_EP(mField,Aij,D,ddF,LAMBDA(young_modulus,poisson_ratio),MU(young_modulus,poisson_ratio),"DISP_r_P");
  ierr = mField.loop_finite_elements("STOCHASIC_PROBLEM","ELASTIC",my_fe_k_rs_EP);  CHKERRQ(ierr);
  
  ierr = VecGhostUpdateBegin(ddF,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
  ierr = VecGhostUpdateEnd(ddF,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
  ierr = VecAssemblyBegin(ddF); CHKERRQ(ierr);
  ierr = VecAssemblyEnd(ddF); CHKERRQ(ierr);
//  ierr = VecView(ddF,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);
  
  ierr = KSPSolve(solver,ddF,ddD_EP); CHKERRQ(ierr);
  ierr = VecGhostUpdateBegin(ddD_EP,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = VecGhostUpdateEnd(ddD_EP,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = mField.set_other_global_VecCreateGhost("STOCHASIC_PROBLEM","DISPLACEMENT","DISP_rs_EP",ROW,ddD_EP,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
  ierr = mField.set_other_global_VecCreateGhost("STOCHASIC_PROBLEM","Lagrange_mul_disp","Lagrange_mul_disp",ROW,ddD_EP,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
//  ierr = VecView(ddD,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);

  //----------------------------------------------------------------------------
  // b. Calculating second-order homogenized stress
  //----------------------------------------------------------------------------
  
  Vec Stress_Homo_rs_EP;
  ierr = VecCreateMPI(PETSC_COMM_WORLD, 6, 6*pcomm->size(), &Stress_Homo_rs_EP);  CHKERRQ(ierr);
  ierr = VecZeroEntries(Stress_Homo_rs_EP); CHKERRQ(ierr);
  
  //    ierr = VecView(D,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);
  ElasticFE_RVELagrange_Homogenized_Stress_Traction MyFE_RVEHomoStressTraction_rs_EP(mField,Aij,D,F,&RVE_volume, applied_strain, Stress_Homo_rs_EP,"DISP_rs_EP","Lagrange_mul_disp",field_rank);
  ierr = mField.loop_finite_elements("STOCHASIC_PROBLEM","Lagrange_elem",MyFE_RVEHomoStressTraction_rs_EP);  CHKERRQ(ierr);
  
  if(pcomm->rank()==0){
    PetscScalar    *avec_rs_EP;
    VecGetArray(Stress_Homo_rs_EP, &avec_rs_EP);
    
    cout<< "\nStress_Homo_rs_EP = \n\n";
    for(int ii=0; ii<6; ii++){
      cout <<*avec_rs_EP<<endl; ;
      avec_rs_EP++;
    }
  }
  cout<< "\n\n";


  // ==================================
  // 5.3 due to Poisson's ratio
  // ===================================
  //----------------------------------------------------------------------------
  // a. Solving the equation to get nodal displacement
  //----------------------------------------------------------------------------
  ierr = VecZeroEntries(ddF); CHKERRQ(ierr);
  ierr = VecGhostUpdateBegin(ddF,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = VecGhostUpdateEnd(ddF,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  
  K_rsPoissonFEMethod my_fe_k_rs_PP(mField,Aij,D,ddF,LAMBDA(young_modulus,poisson_ratio),MU(young_modulus,poisson_ratio),"DISP_r_P");
  ierr = mField.loop_finite_elements("STOCHASIC_PROBLEM","ELASTIC",my_fe_k_rs_PP);  CHKERRQ(ierr);

  ierr = VecGhostUpdateBegin(ddF,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
  ierr = VecGhostUpdateEnd(ddF,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
  ierr = VecAssemblyBegin(ddF); CHKERRQ(ierr);
  ierr = VecAssemblyEnd(ddF); CHKERRQ(ierr);
//  ierr = VecView(ddF,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);
  
  ierr = KSPSolve(solver,ddF,ddD_PP); CHKERRQ(ierr);
  ierr = VecGhostUpdateBegin(ddD_PP,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = VecGhostUpdateEnd(ddD_PP,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  
  ierr = mField.set_other_global_VecCreateGhost("STOCHASIC_PROBLEM","DISPLACEMENT","DISP_rs_PP",ROW,ddD_PP,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
  ierr = mField.set_other_global_VecCreateGhost("STOCHASIC_PROBLEM","Lagrange_mul_disp","Lagrange_mul_disp",ROW,ddD_PP,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
//  ierr = VecView(dD,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);
  
  
  //----------------------------------------------------------------------------
  // b. Calculating second-order homogenized stress
  //----------------------------------------------------------------------------
  Vec Stress_Homo_rs_PP;
  ierr = VecCreateMPI(PETSC_COMM_WORLD, 6, 6*pcomm->size(), &Stress_Homo_rs_PP);  CHKERRQ(ierr);
  ierr = VecZeroEntries(Stress_Homo_rs_PP); CHKERRQ(ierr);
  
  //    ierr = VecView(D,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);
  ElasticFE_RVELagrange_Homogenized_Stress_Traction MyFE_RVEHomoStressTraction_rs_PP(mField,Aij,D,F,&RVE_volume, applied_strain, Stress_Homo_rs_PP,"DISP_rs_PP","Lagrange_mul_disp",field_rank);
  ierr = mField.loop_finite_elements("STOCHASIC_PROBLEM","Lagrange_elem",MyFE_RVEHomoStressTraction_rs_PP);  CHKERRQ(ierr);
  
  if(pcomm->rank()==0){
    PetscScalar    *avec_rs_PP;
    VecGetArray(Stress_Homo_rs_PP, &avec_rs_PP);
    
    cout<< "\nStress_Homo_rs_PP = \n\n";
    for(int ii=0; ii<6; ii++){
      cout.precision(15);
      cout <<*avec_rs_PP<<endl; ;
      avec_rs_PP++;
    }
  }
  cout<< "\n\n";

  // ===========================================================================
  //
  //  VII. OUTPUT
  //
  // =========================================================================== 

  //Save data on mesh
  ierr = write_soltion(mField,"out.vtk","out_post_proc.vtk");   CHKERRQ(ierr);
  
  
  //Destroy matrices
  ierr = VecDestroy(&F); CHKERRQ(ierr);
  ierr = VecDestroy(&dF); CHKERRQ(ierr);
  ierr = VecDestroy(&ddF); CHKERRQ(ierr);
  
  
  ierr = VecDestroy(&D); CHKERRQ(ierr);
  ierr = VecDestroy(&dD_E); CHKERRQ(ierr);
  ierr = VecDestroy(&dD_P); CHKERRQ(ierr);
  ierr = VecDestroy(&ddD_EE); CHKERRQ(ierr);
  ierr = VecDestroy(&ddD_EP); CHKERRQ(ierr);
  ierr = VecDestroy(&ddD_PP); CHKERRQ(ierr);
  
  ierr = MatDestroy(&Aij); CHKERRQ(ierr);
  ierr = KSPDestroy(&solver); CHKERRQ(ierr);
  
  
  ierr = PetscTime(&v2);CHKERRQ(ierr);
  ierr = PetscGetCPUTime(&t2);CHKERRQ(ierr);
  
  PetscSynchronizedPrintf(PETSC_COMM_WORLD,"Total Rank %d Time = %f CPU Time = %f\n",pcomm->rank(),v2-v1,t2-t1);
  PetscSynchronizedFlush(PETSC_COMM_WORLD);
  
  PetscFinalize();
  
}

