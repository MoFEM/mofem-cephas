/* Copyright (C) 2014, 
 *   Xiao-Yi Zhou (xiaoyi.zhou AT newcastle.ac.uk)
 *   Zahur Ullah (Zahur.Ullah@glasgow.ac.uk)
 * --------------------------------------------------------------
 * This routine conducts finite element implementation for stochastic multiscale
 * homogenization problem under linear displacement boundary condition by 
 * perturbation technique for for two-phase composite material consisting of 
 * matrix and inclusion of isotropic material, hence it has four independent 
 * material parameters, which are Young's modulus and Poisson's ration of matrix, 
 * and Young's modulus and Poisson's ratio of fibre/inclusion, to characterize  
 * the mechanical properties of this composite.
 *
 * HISTORY
 *
 * 2014.09.01 (first version)
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

#include "Iso_Rhs_r_PSFEM.hpp"
#include "Iso_Rhs_rs_PSFEM.hpp"

#include "ElasticFE_RVELagrange_Disp.hpp"
#include "ElasticFE_RVELagrange_Homogenized_Stress_Disp.hpp"
#include "RVEVolume.hpp"

#include "PostProcVertexMethod.hpp"
#include "PostProcDisplacementAndStrainOnRefindedMesh.hpp"

using namespace MoFEM;

ErrorCode rval;
PetscErrorCode ierr;

static char help[] = "...\n\n";

// =============================================================================
//  SETTING FOR OUTPUTS
// =============================================================================
PetscErrorCode write_soltion(FieldInterface &mField,const string out_file, const string out_ref_file) {
  PetscFunctionBegin;
  
  PostProcVertexMethod ent_method(mField.get_moab(),"DISPLACEMENT");
  ierr = mField.loop_dofs("STOCHASIC_PROBLEM","DISPLACEMENT",ROW,ent_method); CHKERRQ(ierr);

  const char* args[] = {"_r_Em",    "_r_Pm",    "_r_Ef",    "_r_Pf",    // 1st order
    "_rs_EmEm", "_rs_EmPm", "_rs_EmEf", "_rs_EmPf", // 2nd order
    "_rs_PmPm", "_rs_PmEf", "_rs_PmPf",
    "_rs_EfEf", "_rs_EfPf",
    "_rs_PfPf"};
  vector<string> stochastic_fields(args, args+14);
  
  for(int ii=0; ii<14; ii++ )
  {
    ostringstream ss_field;
    ss_field << "DISP" << stochastic_fields[ii];
    
    PostProcVertexMethod ent_method_r(mField.get_moab(),ss_field.str().c_str());
    ierr = mField.loop_dofs(ss_field.str().c_str(),ent_method_r); CHKERRQ(ierr);
  }
  
  ParallelComm* pcomm = ParallelComm::get_pcomm(&mField.get_moab(),MYPCOMM_INDEX);
  if(pcomm->rank()==0) {
    EntityHandle out_meshset;
    rval = mField.get_moab().create_meshset(MESHSET_SET,out_meshset); CHKERR_PETSC(rval);
    ierr = mField.problem_get_FE("STOCHASIC_PROBLEM","K_Matrix",out_meshset); CHKERRQ(ierr);
    ierr = mField.problem_get_FE("STOCHASIC_PROBLEM","K_Inclusion",out_meshset); CHKERRQ(ierr);
    rval = mField.get_moab().write_file(out_file.c_str(),"VTK","",&out_meshset,1); CHKERR_PETSC(rval);
    rval = mField.get_moab().delete_entities(&out_meshset,1); CHKERR_PETSC(rval);
  }
  PetscFunctionReturn(0);
}


int main(int argc, char *argv[]) {
  
  PetscInitialize(&argc,&argv,(char *)0,help);
  
  Core mb_instance;
  Interface& moab = mb_instance;
  int rank;
  MPI_Comm_rank(PETSC_COMM_WORLD,&rank);
  // ===========================================================================
  //  I. READ MESH DATA AND FINITE ELEMENT ANALYSIS CONTROL PARAMETERS FROM FILE
  // ===========================================================================

  /*****************************************************************************
   * Read parameters from line command
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
    order = 1;
  }
  
  /*****************************************************************************
   * Read Applied strain on the RVE (vector of length 6)
   * strain=[xx, yy, zz, xy, xz, zy]^T
   ****************************************************************************/
  double myapplied_strain[6];
  int nmax=6;
  ierr = PetscOptionsGetRealArray(PETSC_NULL,"-myapplied_strain",myapplied_strain,&nmax,&flg); CHKERRQ(ierr);
  ublas::vector<FieldData> applied_strain;
  applied_strain.resize(6);
  cblas_dcopy(6, &myapplied_strain[0], 1, &applied_strain(0), 1);
  cout<<"applied_strain ="<<applied_strain<<endl;
  
  /*****************************************************************************
   * Transfer mesh data to MOAB database
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
 
  /*****************************************************************************
   * Seting nodal coordinates on the surface to make sure they are periodic
   ****************************************************************************/
  
  Range SurTrisNeg, SurTrisPos;
  ierr = mField.get_Cubit_msId_entities_by_dimension(101,SIDESET,2,SurTrisNeg,true); CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"number of SideSet 101 = %d\n",SurTrisNeg.size()); CHKERRQ(ierr);
  ierr = mField.get_Cubit_msId_entities_by_dimension(102,SIDESET,2,SurTrisPos,true); CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"number of SideSet 102 = %d\n",SurTrisPos.size()); CHKERRQ(ierr);
  
  Range SurNodesNeg,SurNodesPos;
  rval = moab.get_connectivity(SurTrisNeg,SurNodesNeg,true); CHKERR_PETSC(rval);
  cout<<" All nodes on negative surfaces " << SurNodesNeg.size()<<endl;
  rval = moab.get_connectivity(SurTrisPos,SurNodesPos,true); CHKERR_PETSC(rval);
  cout<<" All nodes on positive surfaces " << SurNodesPos.size()<<endl;
  
  double roundfact=1000.0;   double coords_nodes[3];
  //Populating the Multi-index container with nodes on -ve faces
  for(Range::iterator nit = SurNodesNeg.begin(); nit!=SurNodesNeg.end();  nit++) {
    rval = moab.get_coords(&*nit,1,coords_nodes);  CHKERR_PETSC(rval);
    //round values to 3 disimal places
    if(coords_nodes[0]>=0) coords_nodes[0]=double(int(coords_nodes[0]*roundfact+0.5))/roundfact;  else coords_nodes[0]=double(int(coords_nodes[0]*roundfact-0.5))/roundfact;
    if(coords_nodes[1]>=0) coords_nodes[1]=double(int(coords_nodes[1]*roundfact+0.5))/roundfact;  else coords_nodes[1]=double(int(coords_nodes[1]*roundfact-0.5))/roundfact;
    if(coords_nodes[2]>=0) coords_nodes[2]=double(int(coords_nodes[2]*roundfact+0.5))/roundfact;  else coords_nodes[2]=double(int(coords_nodes[2]*roundfact-0.5))/roundfact;
    rval = moab.set_coords(&*nit,1,coords_nodes);  CHKERR_PETSC(rval);
    //      cout<<"   coords_nodes[0]= "<<coords_nodes[0] << "   coords_nodes[1]= "<< coords_nodes[1] << "   coords_nodes[2]= "<< coords_nodes[2] <<endl;
  }
  
  ///Populating the Multi-index container with nodes on +ve faces
  for(Range::iterator nit = SurNodesPos.begin(); nit!=SurNodesPos.end();  nit++) {
    rval = moab.get_coords(&*nit,1,coords_nodes);  CHKERR_PETSC(rval);
    //round values to 3 disimal places
    if(coords_nodes[0]>=0) coords_nodes[0]=double(int(coords_nodes[0]*roundfact+0.5))/roundfact;  else coords_nodes[0]=double(int(coords_nodes[0]*roundfact-0.5))/roundfact;
    if(coords_nodes[1]>=0) coords_nodes[1]=double(int(coords_nodes[1]*roundfact+0.5))/roundfact;  else coords_nodes[1]=double(int(coords_nodes[1]*roundfact-0.5))/roundfact;
    if(coords_nodes[2]>=0) coords_nodes[2]=double(int(coords_nodes[2]*roundfact+0.5))/roundfact;  else coords_nodes[2]=double(int(coords_nodes[2]*roundfact-0.5))/roundfact;
    rval = moab.set_coords(&*nit,1,coords_nodes);  CHKERR_PETSC(rval);
    //      cout<<"   coords_nodes[0]= "<<coords_nodes[0] << "   coords_nodes[1]= "<< coords_nodes[1] << "   coords_nodes[2]= "<< coords_nodes[2] <<endl;
  }

  //----------------------------------------------------------------------------
  //ref meshset ref level 0
  ierr = mField.seed_ref_level_3D(0,0); CHKERRQ(ierr);
  /*****************************************************************************
   *
   * Select element into various mesh-set
   * meshset_level0: all element
   * meshset_Matrix: element set for matrix
   * meshset_Inclusion: element set for inclusion/fibre
   *
   ****************************************************************************/  
  // stl::bitset see for more details
  BitRefLevel bit_level0;
  bit_level0.set(0);
  EntityHandle meshset_level0;
  rval = moab.create_meshset(MESHSET_SET,meshset_level0); CHKERR_PETSC(rval);
  ierr = mField.seed_ref_level_3D(0,bit_level0); CHKERRQ(ierr);
  ierr = mField.get_entities_by_ref_level(bit_level0,BitRefLevel().set(),meshset_level0); CHKERRQ(ierr);
 
  //
  EntityHandle meshset_Matrix, meshset_Inclusion;
  rval = moab.create_meshset(MESHSET_SET,meshset_Matrix); CHKERR_PETSC(rval);
  rval = moab.create_meshset(MESHSET_SET,meshset_Inclusion); CHKERR_PETSC(rval);

  // Element set for matrix
  Range TetsInBlock;
  for(_IT_CUBITMESHSETS_BY_SET_TYPE_FOR_LOOP_(mField,BLOCKSET,it)){
     if(it->get_Cubit_name() == "MAT_ELASTIC") {
	rval = moab.get_entities_by_type(it->meshset, MBTET,TetsInBlock,true); CHKERR_PETSC(rval);
		}
	}
  cout<<"TetsInBlock "<<TetsInBlock<<endl;
  rval = moab.add_entities(meshset_Matrix,TetsInBlock);CHKERR_PETSC(rval);
  
  // Element set for inclusion/fibre
  Range TetsInBlock_stiff_inc;
  for(_IT_CUBITMESHSETS_BY_SET_TYPE_FOR_LOOP_(mField,BLOCKSET,it)){
	if(it->get_Cubit_name() == "MAT_ELASTIC_Stiff_Inclusion") {
	  rval = moab.get_entities_by_type(it->meshset, MBTET,TetsInBlock_stiff_inc,true); CHKERR_PETSC(rval);
		}
	}
  cout<<"TetsInBlock "<<TetsInBlock_stiff_inc<<endl;
  rval = moab.add_entities(meshset_Inclusion,TetsInBlock_stiff_inc);CHKERR_PETSC(rval);

  // ===========================================================================
  // II. DEFINE PROBLEM
  // ===========================================================================
  
  //Fields
  int field_rank=3;
  ierr = mField.add_field("DISPLACEMENT",H1,field_rank); CHKERRQ(ierr);
  ierr = mField.add_field("Lagrange_mul_disp",H1,field_rank); CHKERRQ(ierr);
  
  /*****************************************************************************
   * Add stochastic field
   * (total 14 field for 1st and 2nd order stochastic PSFEM)
   ****************************************************************************/
  const char* args[] = {"_r_Em",    "_r_Pm",    "_r_Ef",    "_r_Pf",    // 1st order
                        "_rs_EmEm", "_rs_EmPm", "_rs_EmEf", "_rs_EmPf", // 2nd order
                                    "_rs_PmPm", "_rs_PmEf", "_rs_PmPf",
                                                "_rs_EfEf", "_rs_EfPf",
                                                            "_rs_PfPf",
    "_r_Em1",    "_r_Pm1",    "_r_Ef1",    "_r_Pf1",    // 1st order
    "_rs_EmEm1", "_rs_EmPm1", "_rs_EmEf1", "_rs_EmPf1", // 2nd order
    "_rs_PmPm1", "_rs_PmEf1", "_rs_PmPf1",
    "_rs_EfEf1", "_rs_EfPf1",
    "_rs_PfPf1",
   
    "_r_Em2",    "_r_Pm2",    "_r_Ef2",    "_r_Pf2",    // 1st order
    "_rs_EmEm2", "_rs_EmPm2", "_rs_EmEf2", "_rs_EmPf2", // 2nd order
    "_rs_PmPm2", "_rs_PmEf2", "_rs_PmPf2",
    "_rs_EfEf2", "_rs_EfPf2",
    "_rs_PfPf2"};
  
  
  vector<string> stochastic_fields(args, args+42);
  
  for(int ii=0; ii<3*14; ii++ )
  {
    ostringstream ss_field;
    ss_field << "DISP" << stochastic_fields[ii];
//    cout<<ss_field.str().c_str()<<endl;
    ierr = mField.add_field(ss_field.str().c_str(),H1,field_rank,MF_ZERO); CHKERRQ(ierr);
  }
  
//  /*****************************************************************************
//   * Create finite element for the defined fields
//   ****************************************************************************/
//  ierr = mField.add_finite_element("K_Matrix"); CHKERRQ(ierr);
//  ierr = mField.add_finite_element("K_Inclusion"); CHKERRQ(ierr);
//  ierr = mField.add_finite_element("Lagrange_elem"); CHKERRQ(ierr);
//
//   /*****************************************************************************
//   *
//   * set field data which finite element use
//   * set field row which finite element use
//   *
//   * [K][U] = [F]                                           Zeroth-order problem
//   * [K][D_r U] = - [D_r K][U]                               First-order problem
//   * [K][H_rs U] = - [H_rs K][U] - 2[D_r K][D_s U]          Second-order problem 
//   *
//   ****************************************************************************/  
//  // Define rows/cols and element data for K_Matrix
//  ierr = mField.modify_finite_element_add_field_row("K_Matrix","DISPLACEMENT"); CHKERRQ(ierr);
//  ierr = mField.modify_finite_element_add_field_col("K_Matrix","DISPLACEMENT"); CHKERRQ(ierr);
//  // required for solving first and second-order problem
//  ierr = mField.modify_finite_element_add_field_data("K_Matrix","DISPLACEMENT"); CHKERRQ(ierr);
//  // required for solving first problem
//  
//  for(int ii=0; ii<4; ii++ )
//  {
//    ostringstream ss_field;
//    ss_field << "DISP" << stochastic_fields[ii];
////    cout<<ss_field.str().c_str()<<endl;
//    ierr = mField.modify_finite_element_add_field_data("K_Matrix",ss_field.str().c_str()); CHKERRQ(ierr);
//  }
//
// // Define rows/cols and element data for K_Inclusion
//  ierr = mField.modify_finite_element_add_field_row("K_Inclusion","DISPLACEMENT"); CHKERRQ(ierr);
//  ierr = mField.modify_finite_element_add_field_col("K_Inclusion","DISPLACEMENT"); CHKERRQ(ierr);
//  // required for first and second-order problem
//  ierr = mField.modify_finite_element_add_field_data("K_Inclusion","DISPLACEMENT"); CHKERRQ(ierr);
//  // required for solving first problem
//  
//  for(int ii=0; ii<4; ii++ )
//  {
//    ostringstream ss_field;
//    ss_field << "DISP" << stochastic_fields[ii];
////    cout<<ss_field.str().c_str()<<endl;
//    ierr = mField.modify_finite_element_add_field_data("K_Inclusion",ss_field.str().c_str()); CHKERRQ(ierr);
//  }
//  
//  // Row and column and data for C and CT matrices
//  //======================================================================================================
//  //C row as Lagrange_mul_disp and col as DISPLACEMENT
//  ierr = mField.modify_finite_element_add_field_row("Lagrange_elem","Lagrange_mul_disp"); CHKERRQ(ierr);
//  ierr = mField.modify_finite_element_add_field_col("Lagrange_elem","DISPLACEMENT"); CHKERRQ(ierr);
//  
//  //CT col as Lagrange_mul_disp and row as DISPLACEMENT
//  ierr = mField.modify_finite_element_add_field_col("Lagrange_elem","Lagrange_mul_disp"); CHKERRQ(ierr);
//  ierr = mField.modify_finite_element_add_field_row("Lagrange_elem","DISPLACEMENT"); CHKERRQ(ierr);
//  
//  //As for stress we need both displacement and temprature (Lukasz)
//  ierr = mField.modify_finite_element_add_field_data("Lagrange_elem","Lagrange_mul_disp"); CHKERRQ(ierr);
//  ierr = mField.modify_finite_element_add_field_data("Lagrange_elem","DISPLACEMENT"); CHKERRQ(ierr);
//  
//  //======================================================================================================
//  //define problems
//  ierr = mField.add_problem("STOCHASIC_PROBLEM"); CHKERRQ(ierr);
//  
//  //set finite elements for problem
//  ierr = mField.modify_problem_add_finite_element("STOCHASIC_PROBLEM","K_Matrix"); CHKERRQ(ierr);
//  ierr = mField.modify_problem_add_finite_element("STOCHASIC_PROBLEM","K_Inclusion"); CHKERRQ(ierr);
//  ierr = mField.modify_problem_add_finite_element("STOCHASIC_PROBLEM","Lagrange_elem"); CHKERRQ(ierr);
//  
//  //set refinment level for problem
//  ierr = mField.modify_problem_ref_level_add_bit("STOCHASIC_PROBLEM",bit_level0); CHKERRQ(ierr);
//
//  // ===========================================================================
//  // III. DECLARE PROBLEM
//  // ===========================================================================
//
//  /*****************************************************************************
//   * Add entitities (by tets) to the field
//   ****************************************************************************/
//  // Zeroth order
//  ierr = mField.add_ents_to_field_by_TETs(0,"DISPLACEMENT",2); CHKERRQ(ierr);
//  
//  for(int ii=0; ii<14; ii++ )
//  {
//    ostringstream ss_field;
//    ss_field << "DISP" << stochastic_fields[ii];
////    cout<<ss_field.str().c_str()<<endl;
//    ierr = mField.add_ents_to_field_by_TETs(0,ss_field.str().c_str()); CHKERRQ(ierr);
//  }
//
//  /*****************************************************************************
//   * Add finite elements entities
//   ****************************************************************************/
//  ierr = mField.add_ents_to_finite_element_by_TETs(meshset_Matrix,"K_Matrix",true); CHKERRQ(ierr);
//  ierr = mField.add_ents_to_finite_element_by_TETs(meshset_Inclusion,"K_Inclusion",true); CHKERRQ(ierr);
//
//  Range SurfacesFaces;
//  ierr = mField.get_Cubit_msId_entities_by_dimension(103,SIDESET,2,SurfacesFaces,true); CHKERRQ(ierr);
//  ierr = PetscPrintf(PETSC_COMM_WORLD,"number of SideSet 103 = %d\n",SurfacesFaces.size()); CHKERRQ(ierr);
//  ierr = mField.add_ents_to_finite_element_by_TRIs(SurfacesFaces,"Lagrange_elem"); CHKERRQ(ierr);
//  
//  //to create meshset from range
//  EntityHandle BoundFacesMeshset;
//  rval = moab.create_meshset(MESHSET_SET,BoundFacesMeshset); CHKERR_PETSC(rval);
//	rval = moab.add_entities(BoundFacesMeshset,SurfacesFaces); CHKERR_PETSC(rval);
//  ierr = mField.seed_ref_level_MESHSET(BoundFacesMeshset,BitRefLevel().set()); CHKERRQ(ierr);
//  ierr = mField.add_ents_to_field_by_TRIs(BoundFacesMeshset,"Lagrange_mul_disp",2); CHKERRQ(ierr);
//
//  /*****************************************************************************
//   * Set applied order
//   ****************************************************************************/
//// int order = 5;
//  ierr = mField.set_field_order(0,MBTET,"DISPLACEMENT",order); CHKERRQ(ierr);
//  ierr = mField.set_field_order(0,MBTRI,"DISPLACEMENT",order); CHKERRQ(ierr);
//  ierr = mField.set_field_order(0,MBEDGE,"DISPLACEMENT",order); CHKERRQ(ierr);
//  ierr = mField.set_field_order(0,MBVERTEX,"DISPLACEMENT",1); CHKERRQ(ierr);
//  
//  int order_st=order;
//  for(int ii=0; ii<14; ii++ )
//  {
//    ostringstream ss_field;
//    ss_field << "DISP" << stochastic_fields[ii];
////    cout<<ss_field.str().c_str()<<endl;
//    ierr = mField.set_field_order(0,MBTET,ss_field.str().c_str(),order_st); CHKERRQ(ierr);
//    ierr = mField.set_field_order(0,MBTRI,ss_field.str().c_str(),order_st); CHKERRQ(ierr);
//    ierr = mField.set_field_order(0,MBEDGE,ss_field.str().c_str(),order_st); CHKERRQ(ierr);
//    ierr = mField.set_field_order(0,MBVERTEX,ss_field.str().c_str(),1); CHKERRQ(ierr);
//  }
//
//  // --------------
//  ierr = mField.set_field_order(0,MBTRI,"Lagrange_mul_disp",order); CHKERRQ(ierr);
//  ierr = mField.set_field_order(0,MBEDGE,"Lagrange_mul_disp",order); CHKERRQ(ierr);
//  ierr = mField.set_field_order(0,MBVERTEX,"Lagrange_mul_disp",1); CHKERRQ(ierr); 
//
//  // ===========================================================================
//  //  IV. BUILD DATABASE
//  // ===========================================================================
//  
//  //build field
//  ierr = mField.build_fields(); CHKERRQ(ierr);
//  
//  //build finite elemnts
//  ierr = mField.build_finite_elements(); CHKERRQ(ierr);
//  
//  //build adjacencies
//  ierr = mField.build_adjacencies(bit_level0); CHKERRQ(ierr);
//  
//  //build problem
//  ierr = mField.build_problems(); CHKERRQ(ierr);
//  
//  // ===========================================================================
//  //  V. MESH PARTITION
//  // ===========================================================================
//  
//  //partition
//  ierr = mField.partition_problem("STOCHASIC_PROBLEM"); CHKERRQ(ierr);
//  ierr = mField.partition_finite_elements("STOCHASIC_PROBLEM"); CHKERRQ(ierr);
//  //what are ghost nodes, see Petsc Manual
//  ierr = mField.partition_ghost_dofs("STOCHASIC_PROBLEM"); CHKERRQ(ierr);
//  
//  //print block sets with materials
//  ierr = mField.print_cubit_materials_set(); CHKERRQ(ierr);
//
//  // ===========================================================================
//  //  VI. SOLUTION PHASE
//  // ===========================================================================
//
//  /*****************************************************************************
//   *  0. PREPARATION FOR PROCESSING SOLVE
//   ****************************************************************************/
//  Vec F,dF,ddF; // vectors
//  Vec D,dD,ddD;
//  
//  ierr = mField.VecCreateGhost("STOCHASIC_PROBLEM",ROW,&F); CHKERRQ(ierr);
//  ierr = mField.VecCreateGhost("STOCHASIC_PROBLEM",ROW,&dF); CHKERRQ(ierr);
//  ierr = mField.VecCreateGhost("STOCHASIC_PROBLEM",ROW,&ddF); CHKERRQ(ierr);
//  ierr = mField.VecCreateGhost("STOCHASIC_PROBLEM",COL,&D); CHKERRQ(ierr);
//  ierr = mField.VecCreateGhost("STOCHASIC_PROBLEM",COL,&dD); CHKERRQ(ierr);
//  ierr = mField.VecCreateGhost("STOCHASIC_PROBLEM",COL,&ddD); CHKERRQ(ierr);
//
//  /*****************************************************************************
//   *  1. Assembling global stiffness matrix K 
//   *     and external force vector F
//   ****************************************************************************/  
//  Mat Aij;
//  ierr = mField.MatCreateMPIAIJWithArrays("STOCHASIC_PROBLEM",&Aij); CHKERRQ(ierr);
//  
//  struct MyElasticFEMethod: public ElasticFEMethod {
//    MyElasticFEMethod(FieldInterface& _mField, Mat &_Aij,Vec &_D,Vec& _F,double _lambda,double _mu):
//    ElasticFEMethod(_mField,_Aij,_D,_F,_lambda,_mu) {};
//    
//    PetscErrorCode Fint(Vec F_int) {
//      PetscFunctionBegin;
//      ierr = ElasticFEMethod::Fint(); CHKERRQ(ierr);
//      for(int rr = 0;rr<row_mat;rr++) {
//        if(RowGlob[rr].size()!=f_int[rr].size()) SETERRQ(PETSC_COMM_SELF,1,"data inconsistency");
//        if(RowGlob[rr].size()==0) continue;
//        f_int[rr] *= -1; //This is not SNES we solve K*D = -RES
//        ierr = VecSetValues(F_int,RowGlob[rr].size(),&(RowGlob[rr])[0],&(f_int[rr].data()[0]),ADD_VALUES); CHKERRQ(ierr);
//      }
//      PetscFunctionReturn(0);
//    }
//    
//  };
//  
//cout<<"solution 1"<<endl;
//  //Assemble F and Aij
//  const double young_modulus = 1;
//  const double poisson_ratio = 0.0;
//  MyElasticFEMethod my_fe(mField,Aij,D,F,LAMBDA(young_modulus,poisson_ratio),MU(young_modulus,poisson_ratio));
//  ElasticFE_RVELagrange_Disp MyFE_RVELagrange(mField,Aij,D,F,applied_strain,"DISPLACEMENT","Lagrange_mul_disp",field_rank);
//
//  ierr = VecZeroEntries(F); CHKERRQ(ierr);
//  ierr = VecGhostUpdateBegin(F,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
//  ierr = VecGhostUpdateEnd(F,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
//  ierr = MatZeroEntries(Aij); CHKERRQ(ierr);
//  
//  ierr = mField.loop_finite_elements("STOCHASIC_PROBLEM","K_Matrix",my_fe);  CHKERRQ(ierr);
//  ierr = mField.loop_finite_elements("STOCHASIC_PROBLEM","K_Inclusion",my_fe);  CHKERRQ(ierr);
//  ierr = mField.loop_finite_elements("STOCHASIC_PROBLEM","Lagrange_elem",MyFE_RVELagrange);  CHKERRQ(ierr);
//  
//  ierr = VecGhostUpdateBegin(F,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
//  ierr = VecGhostUpdateEnd(F,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
//  ierr = VecAssemblyBegin(F); CHKERRQ(ierr);
//  ierr = VecAssemblyEnd(F); CHKERRQ(ierr);
//  ierr = MatAssemblyBegin(Aij,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
//  ierr = MatAssemblyEnd(Aij,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
//  
//  PetscSynchronizedFlush(PETSC_COMM_WORLD);
//  //ierr = VecView(F,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);
//  //ierr = MatView(Aij,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);
//  
//////Matrix View
////MatView(Aij,PETSC_VIEWER_DRAW_WORLD);//PETSC_VIEWER_STDOUT_WORLD);
////std::string wait;
////std::cin >> wait;
//
//cout<<"solution"<<endl;
//  /*****************************************************************************
//   *  2. Get the volume of RVE
//   ****************************************************************************/
//  double RVE_volume;    RVE_volume=0.0;  //RVE volume for full RVE We need this for stress calculation
//  Vec RVE_volume_Vec;
//  ierr = VecCreateMPI(PETSC_COMM_WORLD, 1, pcomm->size(), &RVE_volume_Vec);  CHKERRQ(ierr);
//  ierr = VecZeroEntries(RVE_volume_Vec); CHKERRQ(ierr);
//  
//  RVEVolume MyRVEVol(mField,Aij,D,F,LAMBDA(young_modulus,poisson_ratio),MU(young_modulus,poisson_ratio), RVE_volume_Vec);
//  ierr = mField.loop_finite_elements("STOCHASIC_PROBLEM","K_Matrix",MyRVEVol);  CHKERRQ(ierr);
//  ierr = mField.loop_finite_elements("STOCHASIC_PROBLEM","K_Inclusion",MyRVEVol);  CHKERRQ(ierr);
//  //    ierr = VecView(RVE_volume_Vec,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);
//  ierr = VecSum(RVE_volume_Vec, &RVE_volume);  CHKERRQ(ierr);
//    cout<<"Final RVE_volume = "<< RVE_volume <<endl;
//  
//  /*****************************************************************************
//   *  3. SOLVE THE ZEROTH-ORDER FINITE ELEMENT EQUILIBRIUM EQUATION
//   *     [K][U] = [F]
//  ****************************************************************************/
//
//  //----------------------------------------------------------------------------
//  // 3.1 Solving the equation to get nodal displacement
//  //----------------------------------------------------------------------------
//  //Solver
//  KSP solver;
//  ierr = KSPCreate(PETSC_COMM_WORLD,&solver); CHKERRQ(ierr);
//  ierr = KSPSetOperators(solver,Aij,Aij,SAME_NONZERO_PATTERN); CHKERRQ(ierr);
//  ierr = KSPSetFromOptions(solver); CHKERRQ(ierr);
//  ierr = KSPSetUp(solver); CHKERRQ(ierr);
//
////  ierr = VecView(F,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);
//  ierr = KSPSolve(solver,F,D); CHKERRQ(ierr);
//  ierr = VecGhostUpdateBegin(D,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
//  ierr = VecGhostUpdateEnd(D,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
//  ierr = mField.set_global_VecCreateGhost("STOCHASIC_PROBLEM",ROW,D,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
////  ierr = VecView(D,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);
//  
//  //----------------------------------------------------------------------------
//  // 3.2 Calculating zeroth-order homogenized stress using volume averaging theorem
//  //----------------------------------------------------------------------------  
//  //create a vector for 6 components of homogenized stress
//  Vec Stress_Homo;
//  ierr = VecCreateMPI(PETSC_COMM_WORLD, 6, 6*pcomm->size(), &Stress_Homo);  CHKERRQ(ierr);
//  ierr = VecZeroEntries(Stress_Homo); CHKERRQ(ierr);
//  
//  //    ierr = VecView(D,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);
//  ElasticFE_RVELagrange_Homogenized_Stress_Disp MyFE_RVEHomoStressDisp(mField,Aij,D,F,&RVE_volume, applied_strain, Stress_Homo,"DISPLACEMENT","Lagrange_mul_disp",field_rank);
//  ierr = mField.loop_finite_elements("STOCHASIC_PROBLEM","Lagrange_elem",MyFE_RVEHomoStressDisp);  CHKERRQ(ierr);
//  
//  if(pcomm->rank()==0){
//    PetscScalar    *avec;
//    VecGetArray(Stress_Homo, &avec);
//    
//    cout<< "\nStress_Homo = \n\n";
//    for(int ii=0; ii<6; ii++){
//      cout.precision(15);
//      cout <<*avec<<endl;
//      avec++;
//    }
//  }
//  cout<< "\n\n";
//  
//   /*****************************************************************************
//   *  SOLVE THE FIRST-ORDER AND SECOND-ORDER FINITE ELEMENT EQUILIBRIUM EQUATION
//   *  1st order-[K][U_r] = -[K_r][U}  due to both young and poisson
//   *  2nd order-[K][U_rs] = -[K_rs][U]-2[K_r][U_s]
//   *
//   ****************************************************************************/
//  for(int ii=0; ii<14; ii++){
//    ierr = VecZeroEntries(dF); CHKERRQ(ierr);
//    ierr = VecGhostUpdateBegin(dF,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
//    ierr = VecGhostUpdateEnd(dF,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
//    
//    if (ii==0){//due to Young's modulus of matrix
//      K_rYoungFEMethod my_fe_k_r_Em(mField,Aij,D,dF,LAMBDA(young_modulus,poisson_ratio),MU(young_modulus,poisson_ratio));
//      ierr = mField.loop_finite_elements("STOCHASIC_PROBLEM","K_Matrix",my_fe_k_r_Em);  CHKERRQ(ierr);}
//    else if (ii==1){//due to Poisson's ratio of matrix
//      K_rPoissonFEMethod my_fe_k_r_Pm(mField,Aij,D,dF,LAMBDA(young_modulus,poisson_ratio),MU(young_modulus,poisson_ratio));
//      ierr = mField.loop_finite_elements("STOCHASIC_PROBLEM","K_Matrix",my_fe_k_r_Pm);  CHKERRQ(ierr);}
//    else if (ii==2){//due to Young's modulus of fibre/inclusion
//      K_rYoungFEMethod my_fe_k_r_Ef(mField,Aij,D,dF,LAMBDA(young_modulus,poisson_ratio),MU(young_modulus,poisson_ratio));
//      ierr = mField.loop_finite_elements("STOCHASIC_PROBLEM","K_Inclusion",my_fe_k_r_Ef);  CHKERRQ(ierr);}
//    else if (ii==3){//due to Poisson's ratio of fibre/inclusion
//      K_rPoissonFEMethod my_fe_k_r_Pf(mField,Aij,D,dF,LAMBDA(young_modulus,poisson_ratio),MU(young_modulus,poisson_ratio));
//      ierr = mField.loop_finite_elements("STOCHASIC_PROBLEM","K_Inclusion",my_fe_k_r_Pf);  CHKERRQ(ierr);}
//    else if (ii==4){//due to Young's modulus of matrix
//      K_rsYoungFEMethod my_fe_k_rs_EmEm(mField,Aij,D,dF,LAMBDA(young_modulus,poisson_ratio),MU(young_modulus,poisson_ratio),"DISP_r_Em");
//      ierr = mField.loop_finite_elements("STOCHASIC_PROBLEM","K_Matrix",my_fe_k_rs_EmEm);  CHKERRQ(ierr);}
//    else if (ii==5){//due to Young's modulus and Poisson's ratio of matrix
//      K_rYoungPoissonFEMethod my_fe_k_rs_EmPm(mField,Aij,D,dF,LAMBDA(young_modulus,poisson_ratio),MU(young_modulus,poisson_ratio),"DISP_r_Pm");
//      ierr = mField.loop_finite_elements("STOCHASIC_PROBLEM","K_Matrix",my_fe_k_rs_EmPm);  CHKERRQ(ierr);}
//    else if (ii==6){//due to Young's modulus of matrix and Young's modulus of fibre/inclusion
//      K_rs_EmEPf_FEMethod my_fe_k_rs_EmEf(mField,Aij,D,dF,LAMBDA(young_modulus,poisson_ratio),MU(young_modulus,poisson_ratio),"DISP_r_Ef");
//      ierr = mField.loop_finite_elements("STOCHASIC_PROBLEM","K_Matrix",my_fe_k_rs_EmEf);  CHKERRQ(ierr);}
//    else if (ii==7){//due to Young's modulus of matrix and Poisson's ratio of fibre/inclusion
//      K_rs_EmEPf_FEMethod my_fe_k_rs_EmPf(mField,Aij,D,dF,LAMBDA(young_modulus,poisson_ratio),MU(young_modulus,poisson_ratio),"DISP_r_Pf");
//      ierr = mField.loop_finite_elements("STOCHASIC_PROBLEM","K_Matrix",my_fe_k_rs_EmPf);  CHKERRQ(ierr);}
//    else if (ii==8){//due to Poisson's ratio of matrix
//      K_rsPoissonFEMethod my_fe_k_rs_PmPm(mField,Aij,D,dF,LAMBDA(young_modulus,poisson_ratio),MU(young_modulus,poisson_ratio),"DISP_r_Pm");
//      ierr = mField.loop_finite_elements("STOCHASIC_PROBLEM","K_Matrix",my_fe_k_rs_PmPm);  CHKERRQ(ierr);}
//    else if (ii==9){//due to Poisson's ratio of matrix and Young's modulus of fibre/inclusion
//      K_rs_PmEPf_FEMethod my_fe_k_rs_PmEf(mField,Aij,D,dF,LAMBDA(young_modulus,poisson_ratio),MU(young_modulus,poisson_ratio),"DISP_r_Ef");
//      ierr = mField.loop_finite_elements("STOCHASIC_PROBLEM","K_Matrix",my_fe_k_rs_PmEf);  CHKERRQ(ierr);}
//    else if (ii==10){//due to Poisson's ratio of matrix and Poisson's ratio of fibre/inclusion
//      K_rs_PmEPf_FEMethod my_fe_k_rs_PmPf(mField,Aij,D,dF,LAMBDA(young_modulus,poisson_ratio),MU(young_modulus,poisson_ratio),"DISP_r_Pf");
//      ierr = mField.loop_finite_elements("STOCHASIC_PROBLEM","K_Matrix",my_fe_k_rs_PmPf);  CHKERRQ(ierr);}
//    else if (ii==11){//due to Young's modulus of fibre/inclusion
//      K_rsYoungFEMethod my_fe_k_rs_EfEf(mField,Aij,D,dF,LAMBDA(young_modulus,poisson_ratio),MU(young_modulus,poisson_ratio),"DISP_r_Ef");
//      ierr = mField.loop_finite_elements("STOCHASIC_PROBLEM","K_Inclusion",my_fe_k_rs_EfEf);  CHKERRQ(ierr);}
//    else if (ii==12){//due to Young's modulus and Poisson's ratio of fibre/inclusion
//      K_rYoungPoissonFEMethod my_fe_k_rs_EfPf(mField,Aij,D,dF,LAMBDA(young_modulus,poisson_ratio),MU(young_modulus,poisson_ratio),"DISP_r_Pf");
//      ierr = mField.loop_finite_elements("STOCHASIC_PROBLEM","K_Inclusion",my_fe_k_rs_EfPf);  CHKERRQ(ierr);}
//    else if (ii==13){//due to Poisson's ratio of fibre/inclusion
//      K_rsPoissonFEMethod my_fe_k_rs_PfPf(mField,Aij,D,dF,LAMBDA(young_modulus,poisson_ratio),MU(young_modulus,poisson_ratio),"DISP_r_Pf");
//      ierr = mField.loop_finite_elements("STOCHASIC_PROBLEM","K_Inclusion",my_fe_k_rs_PfPf);  CHKERRQ(ierr);
//    }
//
//    ierr = VecGhostUpdateBegin(dF,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
//    ierr = VecGhostUpdateEnd(dF,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
//    ierr = VecAssemblyBegin(dF); CHKERRQ(ierr);
//    ierr = VecAssemblyEnd(dF); CHKERRQ(ierr);
//    
//    ierr = KSPSolve(solver,dF,dD); CHKERRQ(ierr);
//    ierr = VecGhostUpdateBegin(dD,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
//    ierr = VecGhostUpdateEnd(dD,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
//    
//    ostringstream ss_field;
//    ss_field << "DISP" << stochastic_fields[ii];
//    cout<<ss_field.str().c_str()<<endl;
//    ierr = mField.set_other_global_VecCreateGhost("STOCHASIC_PROBLEM","DISPLACEMENT",ss_field.str().c_str(),ROW,dD,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
//    ierr = mField.set_other_global_VecCreateGhost("STOCHASIC_PROBLEM","Lagrange_mul_disp","Lagrange_mul_disp",ROW,dD,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
//    
//    /*****************************************************************************
//     Calculating homogenized stress
//     ****************************************************************************/
//    Vec Stress_Homo_r;
//    ierr = VecCreateMPI(PETSC_COMM_WORLD, 6, 6*pcomm->size(), &Stress_Homo_r);  CHKERRQ(ierr);
//    ierr = VecZeroEntries(Stress_Homo_r); CHKERRQ(ierr);
//    
//    ElasticFE_RVELagrange_Homogenized_Stress_Disp MyFE_RVEHomoStressDisp_r(mField,Aij,D,F,&RVE_volume, applied_strain, Stress_Homo_r,ss_field.str().c_str(),"Lagrange_mul_disp",field_rank);
//    ierr = mField.loop_finite_elements("STOCHASIC_PROBLEM","Lagrange_elem",MyFE_RVEHomoStressDisp_r);  CHKERRQ(ierr);
//    
//    if(pcomm->rank()==0){
//      PetscScalar    *avec_r;
//      VecGetArray(Stress_Homo_r, &avec_r);
//      
//      cout<< "\nStress_Homo_r_Em = \n\n";
//      for(int ii=0; ii<6; ii++){
//        cout.precision(15);
//        cout <<*avec_r<<endl;
//        avec_r++;
//      }
//    }
//    cout<< "\n\n";
//  }
//  
//  // ===========================================================================
//  //  VII. OUTPUT
//  // ===========================================================================
//  //Save data on mesh
//  ierr = write_soltion(mField,"out.vtk","out_post_proc.vtk");   CHKERRQ(ierr);
//
//  ofstream TheFile;
//  TheFile.open("Result.txt",ofstream::out); 
//  if(pcomm->rank()==0){
//    PetscScalar    *avec;
//    VecGetArray(Stress_Homo, &avec);
//    
//    for(int ii=0; ii<6; ii++){
//      TheFile<<setprecision(15)<<*avec<<'\n';
//      avec++;
//    }
//  }
//  TheFile.close();
//
//  // ===========================================================================
//  //  VIII. FINISH
//  // ===========================================================================
//  //Destroy matrices
//  ierr = VecDestroy(&F); CHKERRQ(ierr);
//  ierr = VecDestroy(&dF); CHKERRQ(ierr);
//
//  ierr = VecDestroy(&D); CHKERRQ(ierr);
//  ierr = VecDestroy(&dD); CHKERRQ(ierr);
//
//  ierr = MatDestroy(&Aij); CHKERRQ(ierr);
//  ierr = KSPDestroy(&solver); CHKERRQ(ierr);
//
//  ierr = PetscTime(&v2);CHKERRQ(ierr);
//  ierr = PetscGetCPUTime(&t2);CHKERRQ(ierr);
  
  PetscSynchronizedPrintf(PETSC_COMM_WORLD,"Total Rank %d Time = %f CPU Time = %f\n",pcomm->rank(),v2-v1,t2-t1);
  PetscSynchronizedFlush(PETSC_COMM_WORLD);
  PetscFinalize();
}

