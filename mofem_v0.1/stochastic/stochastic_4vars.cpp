/* Copyright (C) 2013, Xiao-Yi Zhou (xiaoyi.zhou@ncl.ac.uk)
 * --------------------------------------------------------------
 * This routine performs SFEA for simple elastic problem for two-phase isotropic
 * material which has two independent material properties of Young's modulus and
 * Poisson's ratio for each constituent.
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

/* HISTORY
 *
 * 2014.08.21 (first version)
 *
 */

#include "common.hpp"

#include "FieldInterface.hpp"
#include "FieldCore.hpp"

#include "FEMethod_UpLevelStudent.hpp"
#include "ElasticFEMethod.hpp"

#include "K_rPoissonFEMethod.hpp"
#include "K_rYoungFEMethod.hpp"
#include "K_rsPoissonFEMethod.hpp"
#include "K_rsYoungFEMethod.hpp"
#include "K_rYoungPoissonFEMethod.hpp"


#include "K_rs_EmEPf_FEMethod.hpp"

#include "K_rs_PmEPf_FEMethod.hpp"

#include "ForcesAndSurcesCore.hpp"
#include "SnesCtx.hpp"
#include "TsCtx.hpp"

#ifdef __cplusplus
extern "C" {
#endif
#include<cblas.h>
#include<lapack_wrap.h>
#ifdef __cplusplus
}
#endif

#include "SurfacePressure.hpp"
#include "NodalForce.hpp"
#include "FluidPressure.hpp"
#include "BodyForce.hpp"

#include "ThermalStressElement.hpp"

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/vector.hpp>

#include "Projection10NodeCoordsOnField.hpp"
#include "PostProcVertexMethod.hpp"
#include "PostProcDisplacementAndStrainOnRefindedMesh.hpp"

#include <petscksp.h>

using namespace boost::numeric;
using namespace MoFEM;

ErrorCode rval;
PetscErrorCode ierr;

static char help[] = "...\n\n";

const double young_modulus = 1;
const double poisson_ratio = 0.0;
  // ===========================================================================
  //
  //  SETTING FOR OUTPUTS
  //
  // ===========================================================================
PetscErrorCode write_soltion(FieldInterface &mField,const string out_file, const string out_ref_file) {
  PetscFunctionBegin;
  
  PostProcVertexMethod ent_method(mField.get_moab(),"DISPLACEMENT");
  ierr = mField.loop_dofs("STOCHASIC_PROBLEM","DISPLACEMENT",ROW,ent_method); CHKERRQ(ierr);

  PostProcVertexMethod ent_method_disp_r_Em(mField.get_moab(),"DISP_r_Em");
  ierr = mField.loop_dofs("DISP_r_Em",ent_method_disp_r_Em); CHKERRQ(ierr);
  PostProcVertexMethod ent_method_disp_r_Pm(mField.get_moab(),"DISP_r_Pm");
  ierr = mField.loop_dofs("DISP_r_Pm",ent_method_disp_r_Pm); CHKERRQ(ierr);
  PostProcVertexMethod ent_method_disp_r_Ef(mField.get_moab(),"DISP_r_Ef");
  ierr = mField.loop_dofs("DISP_r_Ef",ent_method_disp_r_Ef); CHKERRQ(ierr);
  PostProcVertexMethod ent_method_disp_r_Pf(mField.get_moab(),"DISP_r_Pf");
  ierr = mField.loop_dofs("DISP_r_Pf",ent_method_disp_r_Pf); CHKERRQ(ierr);
  
  PostProcVertexMethod ent_method_disp_rs_EmEm(mField.get_moab(),"DISP_rs_EmEm");
  ierr = mField.loop_dofs("DISP_rs_EmEm",ent_method_disp_rs_EmEm); CHKERRQ(ierr);
  PostProcVertexMethod ent_method_disp_rs_EmPm(mField.get_moab(),"DISP_rs_EmPm");
  ierr = mField.loop_dofs("DISP_rs_EmPm",ent_method_disp_rs_EmPm); CHKERRQ(ierr);
  PostProcVertexMethod ent_method_disp_rs_EmEf(mField.get_moab(),"DISP_rs_EmEf");
  ierr = mField.loop_dofs("DISP_rs_EmEf",ent_method_disp_rs_EmEf); CHKERRQ(ierr);
  PostProcVertexMethod ent_method_disp_rs_EmPf(mField.get_moab(),"DISP_rs_EmPf");
  ierr = mField.loop_dofs("DISP_rs_EmPf",ent_method_disp_rs_EmPf); CHKERRQ(ierr);  

  PostProcVertexMethod ent_method_disp_rs_PmPm(mField.get_moab(),"DISP_rs_PmPm");
  ierr = mField.loop_dofs("DISP_rs_PmPm",ent_method_disp_rs_PmPm); CHKERRQ(ierr);
  PostProcVertexMethod ent_method_disp_rs_PmEf(mField.get_moab(),"DISP_rs_PmEf");
  ierr = mField.loop_dofs("DISP_rs_PmEf",ent_method_disp_rs_PmEf); CHKERRQ(ierr);
  PostProcVertexMethod ent_method_disp_rs_PmPf(mField.get_moab(),"DISP_rs_PmPf");
  ierr = mField.loop_dofs("DISP_rs_PmPf",ent_method_disp_rs_PmPf); CHKERRQ(ierr);

  PostProcVertexMethod ent_method_disp_rs_EfEf(mField.get_moab(),"DISP_rs_EfEf");
  ierr = mField.loop_dofs("DISP_rs_EfEf",ent_method_disp_rs_EfEf); CHKERRQ(ierr);
  PostProcVertexMethod ent_method_disp_rs_EfPf(mField.get_moab(),"DISP_rs_EfPf");
  ierr = mField.loop_dofs("DISP_rs_EfPf",ent_method_disp_rs_EfPf); CHKERRQ(ierr);

  PostProcVertexMethod ent_method_disp_rs_PfPf(mField.get_moab(),"DISP_rs_PfPf");
  ierr = mField.loop_dofs("DISP_rs_PfPf",ent_method_disp_rs_PfPf); CHKERRQ(ierr);

  ParallelComm* pcomm = ParallelComm::get_pcomm(&mField.get_moab(),MYPCOMM_INDEX);
  
  if(pcomm->rank()==0) {
    EntityHandle out_meshset;
    rval = mField.get_moab().create_meshset(MESHSET_SET,out_meshset); CHKERR_PETSC(rval);
    ierr = mField.problem_get_FE("STOCHASIC_PROBLEM","K",out_meshset); CHKERRQ(ierr);
    rval = mField.get_moab().write_file(out_file.c_str(),"VTK","",&out_meshset,1); CHKERR_PETSC(rval);
    rval = mField.get_moab().delete_entities(&out_meshset,1); CHKERR_PETSC(rval);
  }
  
  //  if(pcomm->rank()==0) {
  //
  //    PostProcDisplacemenysAndStarinAndElasticLinearStressOnRefMesh fe_post_proc_method(
  //      mField,"DISPLACEMENT",LAMBDA(young_modulus,poisson_ratio),MU(young_modulus,poisson_ratio));
  //    fe_post_proc_method.do_broadcast = false;
  //    ierr = mField.loop_finite_elements("STOCHASIC_PROBLEM","ELASTIC",fe_post_proc_method,0,pcomm->size());  CHKERRQ(ierr);
  //    rval = fe_post_proc_method.moab_post_proc.write_file(out_ref_file.c_str(),"VTK",""); CHKERR_PETSC(rval);
  //  }
  
  PetscFunctionReturn(0);
}

int main(int argc, char *argv[]) {
  
  PetscInitialize(&argc,&argv,(char *)0,help);
  
  Core mb_instance;
  Interface& moab = mb_instance;
  int rank;
  MPI_Comm_rank(PETSC_COMM_WORLD,&rank);
  ParallelComm* pcomm = ParallelComm::get_pcomm(&moab,MYPCOMM_INDEX);
  if(pcomm == NULL) pcomm =  new ParallelComm(&moab,PETSC_COMM_WORLD);
  // ===========================================================================
  //
  //  I. READ MESH DATA AND FINITE ELEMENT ANALYSIS CONTROL PARAMETERS FROM FILE
  //
  // ===========================================================================
  //Reade parameters from line command
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
  
  //Read mesh to MOAB
  const char *option;
  //option = "PARALLEL=BCAST_DELETE;"
  //"PARTITION=GEOM_DIMENSION,PARTITION_VAL=3,PARTITION_DISTRIBUTE";//;DEBUG_IO";
  option = "";
  rval = moab.load_file(mesh_file_name, 0, option); CHKERR_PETSC(rval);
  
  //Create MoFEM (Joseph) database
  FieldCore core(moab);
  FieldInterface& mField = core;
  
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
  // all elements
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
  //
  // II. DEFINE PROBLEM
  //
  // ===========================================================================
  //Define problem
  
  //Fields
  ierr = mField.add_field("DISPLACEMENT",H1,3,MF_ZERO); CHKERRQ(ierr);
  ierr = mField.add_field("MESH_NODE_POSITIONS",H1,3,MF_ZERO); CHKERRQ(ierr);
  
  //Adding stochastic field
  // 1st order
  ierr = mField.add_field("DISP_r_Pm",H1,3,MF_ZERO); CHKERRQ(ierr);
  ierr = mField.add_field("DISP_r_Pf",H1,3,MF_ZERO); CHKERRQ(ierr);
  ierr = mField.add_field("DISP_r_Em",H1,3,MF_ZERO); CHKERRQ(ierr);
  ierr = mField.add_field("DISP_r_Ef",H1,3,MF_ZERO); CHKERRQ(ierr);
  
  // 2nd order

  ierr = mField.add_field("DISP_rs_EmEm",H1,3,MF_ZERO); CHKERRQ(ierr);
  ierr = mField.add_field("DISP_rs_EmPm",H1,3,MF_ZERO); CHKERRQ(ierr);
  ierr = mField.add_field("DISP_rs_EmEf",H1,3,MF_ZERO); CHKERRQ(ierr);
  ierr = mField.add_field("DISP_rs_EmPf",H1,3,MF_ZERO); CHKERRQ(ierr);
  
  //ierr = mField.add_field("DISP_rs_PmEm",H1,3,MF_ZERO); CHKERRQ(ierr); // TBD
  ierr = mField.add_field("DISP_rs_PmPm",H1,3,MF_ZERO); CHKERRQ(ierr);
  ierr = mField.add_field("DISP_rs_PmEf",H1,3,MF_ZERO); CHKERRQ(ierr);
  ierr = mField.add_field("DISP_rs_PmPf",H1,3,MF_ZERO); CHKERRQ(ierr);

  
  //ierr = mField.add_field("DISP_rs_EfEm",H1,3,MF_ZERO); CHKERRQ(ierr); // TBD
  //ierr = mField.add_field("DISP_rs_EfPm",H1,3,MF_ZERO); CHKERRQ(ierr); // TBD
  ierr = mField.add_field("DISP_rs_EfEf",H1,3,MF_ZERO); CHKERRQ(ierr); 
  ierr = mField.add_field("DISP_rs_EfPf",H1,3,MF_ZERO); CHKERRQ(ierr);

  
  //ierr = mField.add_field("DISP_rs_PfEm",H1,3,MF_ZERO); CHKERRQ(ierr); // TBD
  //ierr = mField.add_field("DISP_rs_PfPm",H1,3,MF_ZERO); CHKERRQ(ierr); // TBD
  //ierr = mField.add_field("DISP_rs_PfEf",H1,3,MF_ZERO); CHKERRQ(ierr); // TBD
  ierr = mField.add_field("DISP_rs_PfPf",H1,3,MF_ZERO); CHKERRQ(ierr);
 
  /*****************************************************************************
   *
   * Create finite element
   *
   ****************************************************************************/
  // FE
  ierr = mField.add_finite_element("K",MF_ZERO); CHKERRQ(ierr);
  ierr = mField.add_finite_element("K_Matrix",MF_ZERO); CHKERRQ(ierr);
  ierr = mField.add_finite_element("K_Inclusion",MF_ZERO); CHKERRQ(ierr);
  
  /*****************************************************************************
   *
   * set field data which finite element use
   * set field row which finite element use
   *
   * [K][U] = [F]                                           Zeroth-order problem
   * [K][D_b U] = - [D_b K][U]                               First-order problem
   * [K][H_bb U] = - [H_bb K][U] - 2[D_b K][D_b U]          Second-order problem 
   *
   ****************************************************************************/
  // Define rows/cols and element data
  ierr = mField.modify_finite_element_add_field_row("K","DISPLACEMENT"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_col("K","DISPLACEMENT"); CHKERRQ(ierr);
  // required for solving zeroth-order problem
  ierr = mField.modify_finite_element_add_field_data("K","DISPLACEMENT"); CHKERRQ(ierr);

  // Define rows/cols and element data for K_Matrix
  ierr = mField.modify_finite_element_add_field_row("K_Matrix","DISPLACEMENT"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_col("K_Matrix","DISPLACEMENT"); CHKERRQ(ierr);
  // required for solving first and second-order problem
  ierr = mField.modify_finite_element_add_field_data("K_Matrix","DISPLACEMENT"); CHKERRQ(ierr);
  // required for solving second-order problem
  ierr = mField.modify_finite_element_add_field_data("K_Matrix","DISP_r_Em"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_data("K_Matrix","DISP_r_Pm"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_data("K_Matrix","DISP_r_Ef"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_data("K_Matrix","DISP_r_Pf"); CHKERRQ(ierr);


  // Define rows/cols and element data for K_Inclusion
  ierr = mField.modify_finite_element_add_field_row("K_Inclusion","DISPLACEMENT"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_col("K_Inclusion","DISPLACEMENT"); CHKERRQ(ierr);
  // required for first and second-order problem
  ierr = mField.modify_finite_element_add_field_data("K_Inclusion","DISPLACEMENT"); CHKERRQ(ierr);
  // required for second-order problem
  ierr = mField.modify_finite_element_add_field_data("K_Inclusion","DISP_r_Em"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_data("K_Inclusion","DISP_r_Pm"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_data("K_Inclusion","DISP_r_Ef"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_data("K_Inclusion","DISP_r_Pf"); CHKERRQ(ierr);
   
  //
  ierr = mField.modify_finite_element_add_field_data("K","MESH_NODE_POSITIONS"); CHKERRQ(ierr);
  
  // Define problem
  ierr = mField.add_problem("STOCHASIC_PROBLEM"); CHKERRQ(ierr);
  
  // Set finite elements for problem
  ierr = mField.modify_problem_add_finite_element("STOCHASIC_PROBLEM","K"); CHKERRQ(ierr);
  ierr = mField.modify_problem_add_finite_element("STOCHASIC_PROBLEM","K_Matrix"); CHKERRQ(ierr);
  ierr = mField.modify_problem_add_finite_element("STOCHASIC_PROBLEM","K_Inclusion"); CHKERRQ(ierr);
  
  // Set refinment level for problem
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
  ierr = mField.add_ents_to_field_by_TETs(0,"DISPLACEMENT"); CHKERRQ(ierr);
  // First order
  ierr = mField.add_ents_to_field_by_TETs(0,"DISP_r_Em"); CHKERRQ(ierr);
  ierr = mField.add_ents_to_field_by_TETs(0,"DISP_r_Pm"); CHKERRQ(ierr);
  ierr = mField.add_ents_to_field_by_TETs(0,"DISP_r_Ef"); CHKERRQ(ierr);
  ierr = mField.add_ents_to_field_by_TETs(0,"DISP_r_Pf"); CHKERRQ(ierr);
  
  // Second order -1
  ierr = mField.add_ents_to_field_by_TETs(0,"DISP_rs_EmEm"); CHKERRQ(ierr);
  ierr = mField.add_ents_to_field_by_TETs(0,"DISP_rs_EmPm"); CHKERRQ(ierr);
  ierr = mField.add_ents_to_field_by_TETs(0,"DISP_rs_EmEf"); CHKERRQ(ierr);
  ierr = mField.add_ents_to_field_by_TETs(0,"DISP_rs_EmPf"); CHKERRQ(ierr);

  // Second order -2
  
//  ierr = mField.add_ents_to_field_by_TETs(0,"DISP_rs_PmEm"); CHKERRQ(ierr); // TBD
  ierr = mField.add_ents_to_field_by_TETs(0,"DISP_rs_PmPm"); CHKERRQ(ierr);
  ierr = mField.add_ents_to_field_by_TETs(0,"DISP_rs_PmEf"); CHKERRQ(ierr);
  ierr = mField.add_ents_to_field_by_TETs(0,"DISP_rs_PmPf"); CHKERRQ(ierr);

  // Second order -3
  
//  ierr = mField.add_ents_to_field_by_TETs(0,"DISP_rs_EfEm"); CHKERRQ(ierr); // TBD
//  ierr = mField.add_ents_to_field_by_TETs(0,"DISP_rs_EfPm"); CHKERRQ(ierr); // TBD
  ierr = mField.add_ents_to_field_by_TETs(0,"DISP_rs_EfEf"); CHKERRQ(ierr);
  ierr = mField.add_ents_to_field_by_TETs(0,"DISP_rs_EfPf"); CHKERRQ(ierr);

  // Second order -4
//  ierr = mField.add_ents_to_field_by_TETs(0,"DISP_rs_PfEm"); CHKERRQ(ierr); // TBD
//  ierr = mField.add_ents_to_field_by_TETs(0,"DISP_rs_PfPm"); CHKERRQ(ierr); // TBD
//  ierr = mField.add_ents_to_field_by_TETs(0,"DISP_rs_PfEf"); CHKERRQ(ierr); // TBD
  ierr = mField.add_ents_to_field_by_TETs(0,"DISP_rs_PfPf"); CHKERRQ(ierr);

  ierr = mField.add_ents_to_field_by_TETs(0,"MESH_NODE_POSITIONS"); CHKERRQ(ierr);
  
  /*****************************************************************************
   *
   * Add finite elements entities
   *
   ****************************************************************************/ 
  ierr = mField.add_ents_to_finite_element_EntType_by_bit_ref(bit_level0,"K",MBTET); CHKERRQ(ierr);
  ierr = mField.add_ents_to_finite_element_by_TETs(meshset_Matrix,"K_Matrix",true); CHKERRQ(ierr);
  ierr = mField.add_ents_to_finite_element_by_TETs(meshset_Inclusion,"K_Inclusion",true); CHKERRQ(ierr);

  /*****************************************************************************
   *
   * Set applied order
   * See reference for detals:
   *   Ainsworth M. and Coyle J. (2003) Hierarchic finite element bases on 
          unstructured tetrahedral meshes. IJNME, 58(14). pp.2103-2130.
   ****************************************************************************/
  // int order = 5;
  ierr = mField.set_field_order(0,MBTET,"DISPLACEMENT",order); CHKERRQ(ierr);
  ierr = mField.set_field_order(0,MBTRI,"DISPLACEMENT",order); CHKERRQ(ierr);
  ierr = mField.set_field_order(0,MBEDGE,"DISPLACEMENT",order); CHKERRQ(ierr);
  ierr = mField.set_field_order(0,MBVERTEX,"DISPLACEMENT",1); CHKERRQ(ierr);
  
  //
  int order_st=order;
  ierr = mField.set_field_order(0,MBTET,"DISP_r_Em",order_st); CHKERRQ(ierr);
  ierr = mField.set_field_order(0,MBTRI,"DISP_r_Em",order_st); CHKERRQ(ierr);
  ierr = mField.set_field_order(0,MBEDGE,"DISP_r_Em",order_st); CHKERRQ(ierr);
  ierr = mField.set_field_order(0,MBVERTEX,"DISP_r_Em",1); CHKERRQ(ierr);
  
  ierr = mField.set_field_order(0,MBTET,"DISP_r_Pm",order_st); CHKERRQ(ierr);
  ierr = mField.set_field_order(0,MBTRI,"DISP_r_Pm",order_st); CHKERRQ(ierr);
  ierr = mField.set_field_order(0,MBEDGE,"DISP_r_Pm",order_st); CHKERRQ(ierr);
  ierr = mField.set_field_order(0,MBVERTEX,"DISP_r_Pm",1); CHKERRQ(ierr);

  ierr = mField.set_field_order(0,MBTET,"DISP_r_Ef",order_st); CHKERRQ(ierr);
  ierr = mField.set_field_order(0,MBTRI,"DISP_r_Ef",order_st); CHKERRQ(ierr);
  ierr = mField.set_field_order(0,MBEDGE,"DISP_r_Ef",order_st); CHKERRQ(ierr);
  ierr = mField.set_field_order(0,MBVERTEX,"DISP_r_Ef",1); CHKERRQ(ierr);
  
  ierr = mField.set_field_order(0,MBTET,"DISP_r_Pf",order_st); CHKERRQ(ierr);
  ierr = mField.set_field_order(0,MBTRI,"DISP_r_Pf",order_st); CHKERRQ(ierr);
  ierr = mField.set_field_order(0,MBEDGE,"DISP_r_Pf",order_st); CHKERRQ(ierr);
  ierr = mField.set_field_order(0,MBVERTEX,"DISP_r_Pf",1); CHKERRQ(ierr);
  
  // 2nd order -1
  ierr = mField.set_field_order(0,MBTET,"DISP_rs_EmEm",order_st); CHKERRQ(ierr);
  ierr = mField.set_field_order(0,MBTRI,"DISP_rs_EmEm",order_st); CHKERRQ(ierr);
  ierr = mField.set_field_order(0,MBEDGE,"DISP_rs_EmEm",order_st); CHKERRQ(ierr);
  ierr = mField.set_field_order(0,MBVERTEX,"DISP_rs_EmEm",1); CHKERRQ(ierr);
  
  ierr = mField.set_field_order(0,MBTET,"DISP_rs_EmPm",order_st); CHKERRQ(ierr);
  ierr = mField.set_field_order(0,MBTRI,"DISP_rs_EmPm",order_st); CHKERRQ(ierr);
  ierr = mField.set_field_order(0,MBEDGE,"DISP_rs_EmPm",order_st); CHKERRQ(ierr);
  ierr = mField.set_field_order(0,MBVERTEX,"DISP_rs_EmPm",1); CHKERRQ(ierr);
  
  ierr = mField.set_field_order(0,MBTET,"DISP_rs_EmEf",order_st); CHKERRQ(ierr);
  ierr = mField.set_field_order(0,MBTRI,"DISP_rs_EmEf",order_st); CHKERRQ(ierr);
  ierr = mField.set_field_order(0,MBEDGE,"DISP_rs_EmEf",order_st); CHKERRQ(ierr);
  ierr = mField.set_field_order(0,MBVERTEX,"DISP_rs_EmEf",1); CHKERRQ(ierr);
  
  ierr = mField.set_field_order(0,MBTET,"DISP_rs_EmPf",order_st); CHKERRQ(ierr);
  ierr = mField.set_field_order(0,MBTRI,"DISP_rs_EmPf",order_st); CHKERRQ(ierr);
  ierr = mField.set_field_order(0,MBEDGE,"DISP_rs_EmPf",order_st); CHKERRQ(ierr);
  ierr = mField.set_field_order(0,MBVERTEX,"DISP_rs_EmPf",1); CHKERRQ(ierr);

  // 2nd order -2
  /*
  ierr = mField.set_field_order(0,MBTET,"DISP_rs_PmEm",order_st); CHKERRQ(ierr);
  ierr = mField.set_field_order(0,MBTRI,"DISP_rs_PmEm",order_st); CHKERRQ(ierr);
  ierr = mField.set_field_order(0,MBEDGE,"DISP_rs_PmEm",order_st); CHKERRQ(ierr);
  ierr = mField.set_field_order(0,MBVERTEX,"DISP_rs_PmEm",1); CHKERRQ(ierr);
  */
  ierr = mField.set_field_order(0,MBTET,"DISP_rs_PmPm",order_st); CHKERRQ(ierr);
  ierr = mField.set_field_order(0,MBTRI,"DISP_rs_PmPm",order_st); CHKERRQ(ierr);
  ierr = mField.set_field_order(0,MBEDGE,"DISP_rs_PmPm",order_st); CHKERRQ(ierr);
  ierr = mField.set_field_order(0,MBVERTEX,"DISP_rs_PmPm",1); CHKERRQ(ierr);

  ierr = mField.set_field_order(0,MBTET,"DISP_rs_PmEf",order_st); CHKERRQ(ierr);
  ierr = mField.set_field_order(0,MBTRI,"DISP_rs_PmEf",order_st); CHKERRQ(ierr);
  ierr = mField.set_field_order(0,MBEDGE,"DISP_rs_PmEf",order_st); CHKERRQ(ierr);
  ierr = mField.set_field_order(0,MBVERTEX,"DISP_rs_PmEf",1); CHKERRQ(ierr);

  ierr = mField.set_field_order(0,MBTET,"DISP_rs_PmPf",order_st); CHKERRQ(ierr);
  ierr = mField.set_field_order(0,MBTRI,"DISP_rs_PmPf",order_st); CHKERRQ(ierr);
  ierr = mField.set_field_order(0,MBEDGE,"DISP_rs_PmPf",order_st); CHKERRQ(ierr);
  ierr = mField.set_field_order(0,MBVERTEX,"DISP_rs_PmPf",1); CHKERRQ(ierr);

  // 2nd order -3
  /*
  ierr = mField.set_field_order(0,MBTET,"DISP_rs_EfEm",order_st); CHKERRQ(ierr);
  ierr = mField.set_field_order(0,MBTRI,"DISP_rs_EfEm",order_st); CHKERRQ(ierr);
  ierr = mField.set_field_order(0,MBEDGE,"DISP_rs_EfEm",order_st); CHKERRQ(ierr);
  ierr = mField.set_field_order(0,MBVERTEX,"DISP_rs_EfEm",1); CHKERRQ(ierr);

  ierr = mField.set_field_order(0,MBTET,"DISP_rs_EfPm",order_st); CHKERRQ(ierr);
  ierr = mField.set_field_order(0,MBTRI,"DISP_rs_EfPm",order_st); CHKERRQ(ierr);
  ierr = mField.set_field_order(0,MBEDGE,"DISP_rs_EfPm",order_st); CHKERRQ(ierr);
  ierr = mField.set_field_order(0,MBVERTEX,"DISP_rs_EfPm",1); CHKERRQ(ierr);
  */

  ierr = mField.set_field_order(0,MBTET,"DISP_rs_EfEf",order_st); CHKERRQ(ierr);
  ierr = mField.set_field_order(0,MBTRI,"DISP_rs_EfEf",order_st); CHKERRQ(ierr);
  ierr = mField.set_field_order(0,MBEDGE,"DISP_rs_EfEf",order_st); CHKERRQ(ierr);
  ierr = mField.set_field_order(0,MBVERTEX,"DISP_rs_EfEf",1); CHKERRQ(ierr);


  ierr = mField.set_field_order(0,MBTET,"DISP_rs_EfPf",order_st); CHKERRQ(ierr);
  ierr = mField.set_field_order(0,MBTRI,"DISP_rs_EfPf",order_st); CHKERRQ(ierr);
  ierr = mField.set_field_order(0,MBEDGE,"DISP_rs_EfPf",order_st); CHKERRQ(ierr);
  ierr = mField.set_field_order(0,MBVERTEX,"DISP_rs_EfPf",1); CHKERRQ(ierr);

  // 2nd order-4
  /*
  ierr = mField.set_field_order(0,MBTET,"DISP_rs_PfEm",order_st); CHKERRQ(ierr);
  ierr = mField.set_field_order(0,MBTRI,"DISP_rs_PfEm",order_st); CHKERRQ(ierr);
  ierr = mField.set_field_order(0,MBEDGE,"DISP_rs_PfEm",order_st); CHKERRQ(ierr);
  ierr = mField.set_field_order(0,MBVERTEX,"DISP_rs_PfEm",1); CHKERRQ(ierr);

  ierr = mField.set_field_order(0,MBTET,"DISP_rs_PfPm",order_st); CHKERRQ(ierr);
  ierr = mField.set_field_order(0,MBTRI,"DISP_rs_PfPm",order_st); CHKERRQ(ierr);
  ierr = mField.set_field_order(0,MBEDGE,"DISP_rs_PfPm",order_st); CHKERRQ(ierr);
  ierr = mField.set_field_order(0,MBVERTEX,"DISP_rs_PfPm",1); CHKERRQ(ierr);

  ierr = mField.set_field_order(0,MBTET,"DISP_rs_PfEf",order_st); CHKERRQ(ierr);
  ierr = mField.set_field_order(0,MBTRI,"DISP_rs_PfEf",order_st); CHKERRQ(ierr);
  ierr = mField.set_field_order(0,MBEDGE,"DISP_rs_PfEf",order_st); CHKERRQ(ierr);
  ierr = mField.set_field_order(0,MBVERTEX,"DISP_rs_PfEf",1); CHKERRQ(ierr);
  */
  ierr = mField.set_field_order(0,MBTET,"DISP_rs_PfPf",order_st); CHKERRQ(ierr);
  ierr = mField.set_field_order(0,MBTRI,"DISP_rs_PfPf",order_st); CHKERRQ(ierr);
  ierr = mField.set_field_order(0,MBEDGE,"DISP_rs_PfPf",order_st); CHKERRQ(ierr);
  ierr = mField.set_field_order(0,MBVERTEX,"DISP_rs_PfPf",1); CHKERRQ(ierr);
  
  /*
  //
  ierr = mField.set_field_order(0,MBTET,"MESH_NODE_POSITIONS",2); CHKERRQ(ierr);
  ierr = mField.set_field_order(0,MBTRI,"MESH_NODE_POSITIONS",2); CHKERRQ(ierr);
  ierr = mField.set_field_order(0,MBEDGE,"MESH_NODE_POSITIONS",2); CHKERRQ(ierr);
  ierr = mField.set_field_order(0,MBVERTEX,"MESH_NODE_POSITIONS",1); CHKERRQ(ierr);
  */
  ierr = MetaNeummanForces::addNeumannBCElements(mField,"STOCHASIC_PROBLEM","DISPLACEMENT"); CHKERRQ(ierr);
  ierr = MetaNodalForces::addNodalForceElement(mField,"STOCHASIC_PROBLEM","DISPLACEMENT"); CHKERRQ(ierr);
  // ===========================================================================
  //
  //  IV. BUILD DATABASE
  //
  // ===========================================================================
  
  // build field
  ierr = mField.build_fields(); CHKERRQ(ierr);
  
  // build finite elemnts
  ierr = mField.build_finite_elements(); CHKERRQ(ierr);
  
  // build adjacencies
  ierr = mField.build_adjacencies(bit_level0); CHKERRQ(ierr);
  
  // build problem
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
  
  //mField.list_dofs_by_field_name("DISPLACEMENT",true);
  
  //print bcs
  ierr = mField.print_cubit_displacement_set(); CHKERRQ(ierr);
  ierr = mField.print_cubit_force_set(); CHKERRQ(ierr);
  //print block sets with materials
  ierr = mField.print_cubit_materials_set(); CHKERRQ(ierr);

  // ===========================================================================
  //
  //  V. SOLUTION PHASE
  //
  // ===========================================================================
  
  /*****************************************************************************
   *
   *  0. PREPARATION FOR PROCESSING SOLVE
   *
   ****************************************************************************/

  Vec F,dF,ddF; // External force vector or right hand side of finite element equilibrium equation
  Vec D;                                   // 0th solution of nodal displacement
  Vec dD_Em,dD_Pm,dD_Ef,dD_Pf;             // 1st solution of nodal displacement 
  Vec ddD_EmEm,ddD_EmPm,ddD_EmEf,ddD_EmPf; // 2nd solution of nodal displacement
  Vec ddD_PmEm,ddD_PmPm,ddD_PmEf,ddD_PmPf; // 2nd solution of nodal displacement
  Vec ddD_EfEm,ddD_EfPm,ddD_EfEf,ddD_EfPf; // 2nd solution of nodal displacement
  Vec ddD_PfEm,ddD_PfPm,ddD_PfEf,ddD_PfPf; // 2nd solution of nodal displacement
  
  ierr = mField.VecCreateGhost("STOCHASIC_PROBLEM",ROW,&F); CHKERRQ(ierr);
  ierr = mField.VecCreateGhost("STOCHASIC_PROBLEM",ROW,&dF); CHKERRQ(ierr);
  ierr = mField.VecCreateGhost("STOCHASIC_PROBLEM",ROW,&ddF); CHKERRQ(ierr);
  
  ierr = mField.VecCreateGhost("STOCHASIC_PROBLEM",COL,&D); CHKERRQ(ierr);
  // first-order problem
  ierr = mField.VecCreateGhost("STOCHASIC_PROBLEM",COL,&dD_Em); CHKERRQ(ierr);
  ierr = mField.VecCreateGhost("STOCHASIC_PROBLEM",COL,&dD_Pm); CHKERRQ(ierr);
  ierr = mField.VecCreateGhost("STOCHASIC_PROBLEM",COL,&dD_Ef); CHKERRQ(ierr);
  ierr = mField.VecCreateGhost("STOCHASIC_PROBLEM",COL,&dD_Pf); CHKERRQ(ierr);
  // second-order problem
  ierr = mField.VecCreateGhost("STOCHASIC_PROBLEM",COL,&ddD_EmEm); CHKERRQ(ierr);
  ierr = mField.VecCreateGhost("STOCHASIC_PROBLEM",COL,&ddD_EmPm); CHKERRQ(ierr);
  ierr = mField.VecCreateGhost("STOCHASIC_PROBLEM",COL,&ddD_EmEf); CHKERRQ(ierr);
  ierr = mField.VecCreateGhost("STOCHASIC_PROBLEM",COL,&ddD_EmPf); CHKERRQ(ierr);  
  
  ierr = mField.VecCreateGhost("STOCHASIC_PROBLEM",COL,&ddD_PmPm); CHKERRQ(ierr);
  ierr = mField.VecCreateGhost("STOCHASIC_PROBLEM",COL,&ddD_PmEf); CHKERRQ(ierr);
  ierr = mField.VecCreateGhost("STOCHASIC_PROBLEM",COL,&ddD_PmPf); CHKERRQ(ierr);

  ierr = mField.VecCreateGhost("STOCHASIC_PROBLEM",COL,&ddD_EfEf); CHKERRQ(ierr);
  ierr = mField.VecCreateGhost("STOCHASIC_PROBLEM",COL,&ddD_EfPf); CHKERRQ(ierr); 

  ierr = mField.VecCreateGhost("STOCHASIC_PROBLEM",COL,&ddD_PfPf); CHKERRQ(ierr);  

  /*****************************************************************************
   *
   *  1. Assembling global stiffness matrix K 
   *     and external force vector F
   ****************************************************************************/
  Mat Aij;
  ierr = mField.MatCreateMPIAIJWithArrays("STOCHASIC_PROBLEM",&Aij); CHKERRQ(ierr);
  
  struct MyElasticFEMethod: public ElasticFEMethod {
    MyElasticFEMethod(FieldInterface& _mField,Mat &_Aij,Vec &_D,Vec& _F,double _lambda,double _mu):
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
  
  Projection10NodeCoordsOnField ent_method_material(mField,"MESH_NODE_POSITIONS");
  ierr = mField.loop_dofs("MESH_NODE_POSITIONS",ent_method_material); CHKERRQ(ierr);
  
  // Definte boundary condition
  DisplacementBCFEMethodPreAndPostProc my_dirichlet_bc(mField,"DISPLACEMENT",Aij,D,F);
  MyElasticFEMethod my_fe(mField,Aij,D,F,LAMBDA(young_modulus,poisson_ratio),MU(young_modulus,poisson_ratio));
  
  ierr = VecZeroEntries(F); CHKERRQ(ierr);
  ierr = VecGhostUpdateBegin(F,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = VecGhostUpdateEnd(F,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = MatZeroEntries(Aij); CHKERRQ(ierr);
  
  //preproc
  ierr = mField.problem_basic_method_preProcess("STOCHASIC_PROBLEM",my_dirichlet_bc); CHKERRQ(ierr);
  
  // loop over elementss for assembling K
  ierr = mField.loop_finite_elements("STOCHASIC_PROBLEM","K",my_fe);  CHKERRQ(ierr);
  
  // forces and preassures on surface
  boost::ptr_map<string,NeummanForcesSurface> neumann_forces;
  ierr = MetaNeummanForces::setNeumannFiniteElementOperators(mField,neumann_forces,F,"DISPLACEMENT"); CHKERRQ(ierr);
  boost::ptr_map<string,NeummanForcesSurface>::iterator mit = neumann_forces.begin();
  for(;mit!=neumann_forces.end();mit++) {
    ierr = mField.loop_finite_elements("STOCHASIC_PROBLEM",mit->first,mit->second->getLoopFe()); CHKERRQ(ierr);
  }
  // noadl forces
  boost::ptr_map<string,NodalForce> nodal_forces;
  ierr = MetaNodalForces::setNodalForceElementOperators(mField,nodal_forces,F,"DISPLACEMENT"); CHKERRQ(ierr);
  boost::ptr_map<string,NodalForce>::iterator fit = nodal_forces.begin();
  for(;fit!=nodal_forces.end();fit++) {
    ierr = mField.loop_finite_elements("STOCHASIC_PROBLEM",fit->first,fit->second->getLoopFe()); CHKERRQ(ierr);
  }
  // postproc
  ierr = mField.problem_basic_method_postProcess("STOCHASIC_PROBLEM",my_dirichlet_bc); CHKERRQ(ierr);
  
  // set matrix possitives define and symetric for cholesky and icc preceonditionser
  ierr = MatSetOption(Aij,MAT_SPD,PETSC_TRUE); CHKERRQ(ierr);
  
  ierr = VecGhostUpdateBegin(F,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
  ierr = VecGhostUpdateEnd(F,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
  ierr = VecAssemblyBegin(F); CHKERRQ(ierr);
  ierr = VecAssemblyEnd(F); CHKERRQ(ierr);
  
  PetscSynchronizedFlush(PETSC_COMM_WORLD);
  
  //  //Matrix View
  //  MatView(Aij,PETSC_VIEWER_DRAW_WORLD);//PETSC_VIEWER_STDOUT_WORLD);
  //  std::string wait;
  //  std::cin >> wait;

  KSP solver;
  ierr = KSPCreate(PETSC_COMM_WORLD,&solver); CHKERRQ(ierr);
  ierr = KSPSetOperators(solver,Aij,Aij,SAME_NONZERO_PATTERN); CHKERRQ(ierr);
  ierr = KSPSetFromOptions(solver); CHKERRQ(ierr);
  ierr = KSPSetUp(solver); CHKERRQ(ierr);
  
  /*****************************************************************************
   *
   *  2. SOLVE THE ZEROTH-ORDER FINITE ELEMENT EQUILIBRIUM EQUATION
   *     [K][U] = [F]
   *
   ****************************************************************************/
  
  ierr = KSPSolve(solver,F,D); CHKERRQ(ierr);
  ierr = VecGhostUpdateBegin(D,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = VecGhostUpdateEnd(D,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = mField.set_global_VecCreateGhost("STOCHASIC_PROBLEM",ROW,D,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
  
  //  ierr = VecView(F,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);
  //  ierr = VecView(D,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);
  
  
  /*****************************************************************************
   *
   *  3. SOLVE THE FIRST-ORDER FINITE ELEMENT EQUILIBRIUM EQUATION
   *     [K][U_r] = -[K_r][U}  due to both young and poisson
   *
   ****************************************************************************/
  // ==================================
  // 3.1 due to Young modulus of matrix
  // ==================================

  // 3.1.1
  ierr = VecZeroEntries(dF); CHKERRQ(ierr);
  ierr = VecGhostUpdateBegin(dF,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = VecGhostUpdateEnd(dF,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  
  // 3.1.2 define/apply boundary condition
  DisplacementBCFEMethodPreAndPostProc my_dirichlet_bc_dEm(mField,"DISPLACEMENT",Aij,dD_Em,dF);
  ierr = mField.problem_basic_method_preProcess("STOCHASIC_PROBLEM",my_dirichlet_bc_dEm); CHKERRQ(ierr);
  
  // 3.1.3 construct the first-order derivative of global stiffness matrix
  K_rYoungFEMethod my_fe_kr_Em(mField,Aij,D,dF,LAMBDA(young_modulus,poisson_ratio),MU(young_modulus,poisson_ratio));
  ierr = mField.loop_finite_elements("STOCHASIC_PROBLEM","K_Matrix",my_fe_kr_Em);  CHKERRQ(ierr);
  
  // 3.1.4
  ierr = mField.problem_basic_method_postProcess("STOCHASIC_PROBLEM",my_dirichlet_bc_dEm); CHKERRQ(ierr);
  // 3.1.5
  ierr = VecGhostUpdateBegin(dF,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
  ierr = VecGhostUpdateEnd(dF,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
  ierr = VecAssemblyBegin(dF); CHKERRQ(ierr);
  ierr = VecAssemblyEnd(dF); CHKERRQ(ierr);
  //  ierr = VecView(dF,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);
  
  // 3.1.6 solve the finite element equation to get answer for
  //       the first order partial derivative of nodal displacement
  ierr = KSPSolve(solver,dF,dD_Em); CHKERRQ(ierr);
  ierr = VecGhostUpdateBegin(dD_Em,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = VecGhostUpdateEnd(dD_Em,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  // 3.1.7
//  ierr = VecView(dD_Em,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);
  ierr = mField.set_other_global_VecCreateGhost("STOCHASIC_PROBLEM","DISPLACEMENT","DISP_r_Em",ROW,dD_Em,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);

  // ====================================
  // 3.2 due to Poisson's ratio of matrix
  // ====================================
  // 3.2.1
  ierr = VecZeroEntries(dF); CHKERRQ(ierr);
  ierr = VecGhostUpdateBegin(dF,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = VecGhostUpdateEnd(dF,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  
  // 3.2.2 define/apply boundary condition
  DisplacementBCFEMethodPreAndPostProc my_dirichlet_bc_dPm(mField,"DISPLACEMENT",Aij,dD_Pm,dF);
  ierr = mField.problem_basic_method_preProcess("STOCHASIC_PROBLEM",my_dirichlet_bc_dPm); CHKERRQ(ierr);
  
  // 3.2.3 construct the first-order derivative of global stiffness matrix
  K_rPoissonFEMethod my_fe_kr_Pm(mField,Aij,D,dF,LAMBDA(young_modulus,poisson_ratio),MU(young_modulus,poisson_ratio));
  ierr = mField.loop_finite_elements("STOCHASIC_PROBLEM","K_Matrix",my_fe_kr_Pm);  CHKERRQ(ierr);
  
  // 3.2.4
  ierr = mField.problem_basic_method_postProcess("STOCHASIC_PROBLEM",my_dirichlet_bc_dPm); CHKERRQ(ierr);
  // 3.2.5
  ierr = VecGhostUpdateBegin(dF,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
  ierr = VecGhostUpdateEnd(dF,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
  ierr = VecAssemblyBegin(dF); CHKERRQ(ierr);
  ierr = VecAssemblyEnd(dF); CHKERRQ(ierr);
//  ierr = VecView(dF,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);
  
  // 3.2.6 solve the finite element equation to get answer for
  //       the first order partial derivative of nodal displacement
  ierr = KSPSolve(solver,dF,dD_Pm); CHKERRQ(ierr);
  ierr = VecGhostUpdateBegin(dD_Pm,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = VecGhostUpdateEnd(dD_Pm,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  // 3.2.7
//  ierr = VecView(dD_P,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);
  ierr = mField.set_other_global_VecCreateGhost("STOCHASIC_PROBLEM","DISPLACEMENT","DISP_r_Pm",ROW,dD_Pm,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
  
  // ===========================================
  // 3.3 due to Young modulus of inclusion/fibre
  // ===========================================

  // 3.3.1
  ierr = VecZeroEntries(dF); CHKERRQ(ierr);
  ierr = VecGhostUpdateBegin(dF,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = VecGhostUpdateEnd(dF,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  
  // 3.3.2 define/apply boundary condition
  DisplacementBCFEMethodPreAndPostProc my_dirichlet_bc_dEf(mField,"DISPLACEMENT",Aij,dD_Ef,dF);
  ierr = mField.problem_basic_method_preProcess("STOCHASIC_PROBLEM",my_dirichlet_bc_dEf); CHKERRQ(ierr);
  
  // 3.3.3 construct the first-order derivative of global stiffness matrix
  K_rYoungFEMethod my_fe_kr_Ef(mField,Aij,D,dF,LAMBDA(young_modulus,poisson_ratio),MU(young_modulus,poisson_ratio));
  ierr = mField.loop_finite_elements("STOCHASIC_PROBLEM","K_Inclusion",my_fe_kr_Ef);  CHKERRQ(ierr);
  
  // 3.3.4
  ierr = mField.problem_basic_method_postProcess("STOCHASIC_PROBLEM",my_dirichlet_bc_dEf); CHKERRQ(ierr);
  // 3.3.5
  ierr = VecGhostUpdateBegin(dF,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
  ierr = VecGhostUpdateEnd(dF,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
  ierr = VecAssemblyBegin(dF); CHKERRQ(ierr);
  ierr = VecAssemblyEnd(dF); CHKERRQ(ierr);
  //  ierr = VecView(dF,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);
  
  // 3.3.6 solve the finite element equation to get answer for
  //       the first order partial derivative of nodal displacement
  ierr = KSPSolve(solver,dF,dD_Ef); CHKERRQ(ierr);
  ierr = VecGhostUpdateBegin(dD_Ef,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = VecGhostUpdateEnd(dD_Ef,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  // 3.3.7
//  ierr = VecView(dD_Em,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);
  ierr = mField.set_other_global_VecCreateGhost("STOCHASIC_PROBLEM","DISPLACEMENT","DISP_r_Ef",ROW,dD_Ef,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);

  // =============================================
  // 3.4 due to Poisson's ratio of inclusion/fiber
  // =============================================
  // 3.4.1
  ierr = VecZeroEntries(dF); CHKERRQ(ierr);
  ierr = VecGhostUpdateBegin(dF,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = VecGhostUpdateEnd(dF,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  
  // 3.4.2 define/apply boundary condition
  DisplacementBCFEMethodPreAndPostProc my_dirichlet_bc_dPf(mField,"DISPLACEMENT",Aij,dD_Pf,dF);
  ierr = mField.problem_basic_method_preProcess("STOCHASIC_PROBLEM",my_dirichlet_bc_dPf); CHKERRQ(ierr);
  
  // 3.4.3 construct the first-order derivative of global stiffness matrix
  K_rPoissonFEMethod my_fe_kr_Pf(mField,Aij,D,dF,LAMBDA(young_modulus,poisson_ratio),MU(young_modulus,poisson_ratio));
  ierr = mField.loop_finite_elements("STOCHASIC_PROBLEM","K_Inclusion",my_fe_kr_Pf);  CHKERRQ(ierr);
  
  // 3.4.4
  ierr = mField.problem_basic_method_postProcess("STOCHASIC_PROBLEM",my_dirichlet_bc_dPf); CHKERRQ(ierr);
  // 3.4.5
  ierr = VecGhostUpdateBegin(dF,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
  ierr = VecGhostUpdateEnd(dF,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
  ierr = VecAssemblyBegin(dF); CHKERRQ(ierr);
  ierr = VecAssemblyEnd(dF); CHKERRQ(ierr);
//  ierr = VecView(dF,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);
  
  // 3.4.6 solve the finite element equation to get answer for
  //       the first order partial derivative of nodal displacement
  ierr = KSPSolve(solver,dF,dD_Pf); CHKERRQ(ierr);
  ierr = VecGhostUpdateBegin(dD_Pf,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = VecGhostUpdateEnd(dD_Pf,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  // 3.4.7
//  ierr = VecView(dD_P,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);
  ierr = mField.set_other_global_VecCreateGhost("STOCHASIC_PROBLEM","DISPLACEMENT","DISP_r_Pf",ROW,dD_Pf,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr); 
  /*****************************************************************************
   *
   *  4. SOLVE THE SECOND-ORDER FINITE ELEMENT EQUILIBRIUM EQUATION
   *     [K][U_rs] = -[K_rs][U]-2[K_r][U_s]
   *
   ****************************************************************************/
  
  // ==================================
  // 4.1.1 due to Em Em
  // ==================================
  // Vec ddD_EmEm,ddD_EmPm,ddD_EmEf,ddD_EmPf; // 2nd solution of nodal displacement
  // a.
  ierr = VecZeroEntries(ddF); CHKERRQ(ierr);
  ierr = VecGhostUpdateBegin(ddF,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = VecGhostUpdateEnd(ddF,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);

  // b. define/apply boundary condition
  DisplacementBCFEMethodPreAndPostProc my_dirichlet_bc_dEmEm(mField,"DISPLACEMENT",Aij,ddD_EmEm,ddF);
  K_rsYoungFEMethod my_fe_k_rs_EmEm(mField,Aij,D,ddF,LAMBDA(young_modulus,poisson_ratio),MU(young_modulus,poisson_ratio),"DISP_r_Em");

  // c. construct the first-order derivative of global stiffness matrix
  ierr = mField.problem_basic_method_preProcess("STOCHASIC_PROBLEM",my_dirichlet_bc_dEmEm); CHKERRQ(ierr);
  ierr = mField.loop_finite_elements("STOCHASIC_PROBLEM","K_Matrix",my_fe_k_rs_EmEm);  CHKERRQ(ierr);
  
  // d.
  ierr = mField.problem_basic_method_postProcess("STOCHASIC_PROBLEM",my_dirichlet_bc_dEmEm); CHKERRQ(ierr);
  
  // e.
  ierr = VecGhostUpdateBegin(ddF,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
  ierr = VecGhostUpdateEnd(ddF,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
  ierr = VecAssemblyBegin(ddF); CHKERRQ(ierr);
  ierr = VecAssemblyEnd(ddF); CHKERRQ(ierr);
  //  ierr = VecView(ddF,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);
  
  // f. solve the finite element equation to get answer for
  //     the second order partial derivative of nodal displacement
  ierr = KSPSolve(solver,ddF,ddD_EmEm); CHKERRQ(ierr);
  ierr = VecGhostUpdateBegin(ddD_EmEm,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = VecGhostUpdateEnd(ddD_EmEm,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  //  ierr = VecView(ddD,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);

  ierr = mField.set_other_global_VecCreateGhost("STOCHASIC_PROBLEM","DISPLACEMENT","DISP_rs_EmEm",ROW,ddD_EmEm,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
  

  // ==================================
  // 4.1.2 due to Em Pm
  // ==================================
  // a.
  ierr = VecZeroEntries(ddF); CHKERRQ(ierr);
  ierr = VecGhostUpdateBegin(ddF,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = VecGhostUpdateEnd(ddF,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);

  // b. define/apply boundary condition
  DisplacementBCFEMethodPreAndPostProc my_dirichlet_bc_dEmPm(mField,"DISPLACEMENT",Aij,ddD_EmPm,ddF);
  K_rYoungPoissonFEMethod my_fe_k_rs_EmPm(mField,Aij,D,ddF,LAMBDA(young_modulus,poisson_ratio),MU(young_modulus,poisson_ratio),"DISP_r_Pm");

  // c. construct the first-order derivative of global stiffness matrix
  ierr = mField.problem_basic_method_preProcess("STOCHASIC_PROBLEM",my_dirichlet_bc_dEmPm); CHKERRQ(ierr);
  ierr = mField.loop_finite_elements("STOCHASIC_PROBLEM","K_Matrix",my_fe_k_rs_EmPm);  CHKERRQ(ierr);

  // d.
  ierr = mField.problem_basic_method_postProcess("STOCHASIC_PROBLEM",my_dirichlet_bc_dEmPm); CHKERRQ(ierr);

  // e.
  ierr = VecGhostUpdateBegin(ddF,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
  ierr = VecGhostUpdateEnd(ddF,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
  ierr = VecAssemblyBegin(ddF); CHKERRQ(ierr);
  ierr = VecAssemblyEnd(ddF); CHKERRQ(ierr);
  //  ierr = VecView(ddF,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);

  // f. solve the finite element equation to get answer for
  //     the second order partial derivative of nodal displacement
  ierr = KSPSolve(solver,ddF,ddD_EmPm); CHKERRQ(ierr);
  ierr = VecGhostUpdateBegin(ddD_EmPm,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = VecGhostUpdateEnd(ddD_EmPm,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  //  ierr = VecView(ddD,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);

  ierr = mField.set_other_global_VecCreateGhost("STOCHASIC_PROBLEM","DISPLACEMENT","DISP_rs_EmPm",ROW,ddD_EmPm,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
  

  // ==================================
  // 4.1.3 due to Em Ef
  // ==================================
  // a.
  ierr = VecZeroEntries(ddF); CHKERRQ(ierr);
  ierr = VecGhostUpdateBegin(ddF,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = VecGhostUpdateEnd(ddF,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  
  // b. define/apply boundary condition
  DisplacementBCFEMethodPreAndPostProc my_dirichlet_bc_dEmEf(mField,"DISPLACEMENT",Aij,ddD_EmEf,ddF);
  K_rs_EmEPf_FEMethod my_fe_k_rs_EmEf(mField,Aij,D,ddF,LAMBDA(young_modulus,poisson_ratio),MU(young_modulus,poisson_ratio),"DISP_r_Ef");
  
  // c. construct the first-order derivative of global stiffness matrix
  ierr = mField.problem_basic_method_preProcess("STOCHASIC_PROBLEM",my_dirichlet_bc_dEmEf); CHKERRQ(ierr);
  ierr = mField.loop_finite_elements("STOCHASIC_PROBLEM","K_Matrix",my_fe_k_rs_EmEf);  CHKERRQ(ierr);
  
  // d.
  ierr = mField.problem_basic_method_postProcess("STOCHASIC_PROBLEM",my_dirichlet_bc_dEmEf); CHKERRQ(ierr);
  
  // e.
  ierr = VecGhostUpdateBegin(ddF,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
  ierr = VecGhostUpdateEnd(ddF,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
  ierr = VecAssemblyBegin(ddF); CHKERRQ(ierr);
  ierr = VecAssemblyEnd(ddF); CHKERRQ(ierr);
  //  ierr = VecView(ddF,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);
  
  // f. solve the finite element equation to get answer for
  //     the second order partial derivative of nodal displacement
  ierr = KSPSolve(solver,ddF,ddD_EmEf); CHKERRQ(ierr);
  ierr = VecGhostUpdateBegin(ddD_EmEf,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = VecGhostUpdateEnd(ddD_EmEf,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  //  ierr = VecView(ddD,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);

  ierr = mField.set_other_global_VecCreateGhost("STOCHASIC_PROBLEM","DISPLACEMENT","DISP_rs_EmEf",ROW,ddD_EmEf,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);  


  // ==================================
  // 4.1.4 due to Em Pf
  // ==================================
  // a.
  ierr = VecZeroEntries(ddF); CHKERRQ(ierr);
  ierr = VecGhostUpdateBegin(ddF,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = VecGhostUpdateEnd(ddF,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);

  // b. define/apply boundary condition
  DisplacementBCFEMethodPreAndPostProc my_dirichlet_bc_dEmPf(mField,"DISPLACEMENT",Aij,ddD_EmPf,ddF);
  K_rs_EmEPf_FEMethod my_fe_k_rs_EmPf(mField,Aij,D,ddF,LAMBDA(young_modulus,poisson_ratio),MU(young_modulus,poisson_ratio),"DISP_r_Pf");
 
  // c. construct the first-order derivative of global stiffness matrix
  ierr = mField.problem_basic_method_preProcess("STOCHASIC_PROBLEM",my_dirichlet_bc_dEmPf); CHKERRQ(ierr);
  ierr = mField.loop_finite_elements("STOCHASIC_PROBLEM","K_Matrix",my_fe_k_rs_EmPf);  CHKERRQ(ierr);

  // d.
  ierr = mField.problem_basic_method_postProcess("STOCHASIC_PROBLEM",my_dirichlet_bc_dEmPf); CHKERRQ(ierr);

  // e.
  ierr = VecGhostUpdateBegin(ddF,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
  ierr = VecGhostUpdateEnd(ddF,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
  ierr = VecAssemblyBegin(ddF); CHKERRQ(ierr);
  ierr = VecAssemblyEnd(ddF); CHKERRQ(ierr);
  //  ierr = VecView(ddF,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);

  // f. solve the finite element equation to get answer for
  //     the second order partial derivative of nodal displacement
  ierr = KSPSolve(solver,ddF,ddD_EmPf); CHKERRQ(ierr);
  ierr = VecGhostUpdateBegin(ddD_EmPf,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = VecGhostUpdateEnd(ddD_EmPf,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  //  ierr = VecView(ddD,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);

  ierr = mField.set_other_global_VecCreateGhost("STOCHASIC_PROBLEM","DISPLACEMENT","DISP_rs_EmPf",ROW,ddD_EmPf,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
  
  // ==================================
  // 4.2.1 due to Pm Em
  // ==================================
  
  
  // ==================================
  // 4.2.2 due to Pm Pm
  // ==================================
  // a.
  ierr = VecZeroEntries(ddF); CHKERRQ(ierr);
  ierr = VecGhostUpdateBegin(ddF,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = VecGhostUpdateEnd(ddF,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);

  // b. define/apply boundary condition
  DisplacementBCFEMethodPreAndPostProc my_dirichlet_bc_dPmPm(mField,"DISPLACEMENT",Aij,ddD_PmPm,ddF);
  K_rsPoissonFEMethod my_fe_k_rs_PmPm(mField,Aij,D,ddF,LAMBDA(young_modulus,poisson_ratio),MU(young_modulus,poisson_ratio),"DISP_r_Pm");

  // c. construct the first-order derivative of global stiffness matrix
  ierr = mField.problem_basic_method_preProcess("STOCHASIC_PROBLEM",my_dirichlet_bc_dPmPm); CHKERRQ(ierr);
  ierr = mField.loop_finite_elements("STOCHASIC_PROBLEM","K_Matrix",my_fe_k_rs_PmPm);  CHKERRQ(ierr);

  // d.
  ierr = mField.problem_basic_method_postProcess("STOCHASIC_PROBLEM",my_dirichlet_bc_dPmPm); CHKERRQ(ierr);

  // e.
  ierr = VecGhostUpdateBegin(ddF,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
  ierr = VecGhostUpdateEnd(ddF,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
  ierr = VecAssemblyBegin(ddF); CHKERRQ(ierr);
  ierr = VecAssemblyEnd(ddF); CHKERRQ(ierr);
  //  ierr = VecView(ddF,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);

  // f. solve the finite element equation to get answer for
  //     the second order partial derivative of nodal displacement
  ierr = KSPSolve(solver,ddF,ddD_PmPm); CHKERRQ(ierr);
  ierr = VecGhostUpdateBegin(ddD_PmPm,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = VecGhostUpdateEnd(ddD_PmPm,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  //  ierr = VecView(ddD,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);

  ierr = mField.set_other_global_VecCreateGhost("STOCHASIC_PROBLEM","DISPLACEMENT","DISP_rs_PmPm",ROW,ddD_PmPm,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);


  // ==================================
  // 4.2.3 due to Pm Ef
  // ==================================
  // a.
  ierr = VecZeroEntries(ddF); CHKERRQ(ierr);
  ierr = VecGhostUpdateBegin(ddF,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = VecGhostUpdateEnd(ddF,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);

  // b. define/apply boundary condition
  DisplacementBCFEMethodPreAndPostProc my_dirichlet_bc_dPmEf(mField,"DISPLACEMENT",Aij,ddD_PmEf,ddF);
  K_rs_PmEPf_FEMethod my_fe_k_rs_PmEf(mField,Aij,D,ddF,LAMBDA(young_modulus,poisson_ratio),MU(young_modulus,poisson_ratio),"DISP_r_Ef");

  // c. construct the first-order derivative of global stiffness matrix
  ierr = mField.problem_basic_method_preProcess("STOCHASIC_PROBLEM",my_dirichlet_bc_dPmEf); CHKERRQ(ierr);
  ierr = mField.loop_finite_elements("STOCHASIC_PROBLEM","K_Matrix",my_fe_k_rs_PmEf);  CHKERRQ(ierr);

  // d.
  ierr = mField.problem_basic_method_postProcess("STOCHASIC_PROBLEM",my_dirichlet_bc_dPmEf); CHKERRQ(ierr);

  // e.
  ierr = VecGhostUpdateBegin(ddF,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
  ierr = VecGhostUpdateEnd(ddF,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
  ierr = VecAssemblyBegin(ddF); CHKERRQ(ierr);
  ierr = VecAssemblyEnd(ddF); CHKERRQ(ierr);
  //  ierr = VecView(ddF,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);

  // f. solve the finite element equation to get answer for
  //     the second order partial derivative of nodal displacement
  ierr = KSPSolve(solver,ddF,ddD_PmEf); CHKERRQ(ierr);
  ierr = VecGhostUpdateBegin(ddD_PmEf,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = VecGhostUpdateEnd(ddD_PmEf,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  //  ierr = VecView(ddD,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);

  ierr = mField.set_other_global_VecCreateGhost("STOCHASIC_PROBLEM","DISPLACEMENT","DISP_rs_PmEf",ROW,ddD_PmEf,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);


  // ==================================
  // 4.2.4 due to Pm Pf
  // ==================================
  // a.
  ierr = VecZeroEntries(ddF); CHKERRQ(ierr);
  ierr = VecGhostUpdateBegin(ddF,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = VecGhostUpdateEnd(ddF,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);

  // b. define/apply boundary condition
  DisplacementBCFEMethodPreAndPostProc my_dirichlet_bc_dPmPf(mField,"DISPLACEMENT",Aij,ddD_PmPf,ddF);
  K_rs_PmEPf_FEMethod my_fe_k_rs_PmPf(mField,Aij,D,ddF,LAMBDA(young_modulus,poisson_ratio),MU(young_modulus,poisson_ratio),"DISP_r_Pf");

  // c. construct the first-order derivative of global stiffness matrix
  ierr = mField.problem_basic_method_preProcess("STOCHASIC_PROBLEM",my_dirichlet_bc_dPmPf); CHKERRQ(ierr);
  ierr = mField.loop_finite_elements("STOCHASIC_PROBLEM","K_Matrix",my_fe_k_rs_PmPf);  CHKERRQ(ierr);

  // d.
  ierr = mField.problem_basic_method_postProcess("STOCHASIC_PROBLEM",my_dirichlet_bc_dPmPf); CHKERRQ(ierr);

  // e.
  ierr = VecGhostUpdateBegin(ddF,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
  ierr = VecGhostUpdateEnd(ddF,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
  ierr = VecAssemblyBegin(ddF); CHKERRQ(ierr);
  ierr = VecAssemblyEnd(ddF); CHKERRQ(ierr);
  //  ierr = VecView(ddF,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);

  // f. solve the finite element equation to get answer for
  //     the second order partial derivative of nodal displacement
  ierr = KSPSolve(solver,ddF,ddD_PmPf); CHKERRQ(ierr);
  ierr = VecGhostUpdateBegin(ddD_PmPf,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = VecGhostUpdateEnd(ddD_PmPf,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  //  ierr = VecView(ddD,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);

  ierr = mField.set_other_global_VecCreateGhost("STOCHASIC_PROBLEM","DISPLACEMENT","DISP_rs_PmPf",ROW,ddD_PmPf,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);


  // ==================================
  // 4.3.1 due to Ef Em
  // ==================================


  // ==================================
  // 4.3.2 due to Ef Pm
  // ==================================


  // ==================================
  // 4.3.3 due to Ef Ef
  // ==================================
  // a.
  ierr = VecZeroEntries(ddF); CHKERRQ(ierr);
  ierr = VecGhostUpdateBegin(ddF,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = VecGhostUpdateEnd(ddF,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);

  // b. define/apply boundary condition
  DisplacementBCFEMethodPreAndPostProc my_dirichlet_bc_dEfEf(mField,"DISPLACEMENT",Aij,ddD_EfEf,ddF);
  K_rsYoungFEMethod my_fe_k_rs_EfEf(mField,Aij,D,ddF,LAMBDA(young_modulus,poisson_ratio),MU(young_modulus,poisson_ratio),"DISP_r_Ef");

  // c. construct the first-order derivative of global stiffness matrix
  ierr = mField.problem_basic_method_preProcess("STOCHASIC_PROBLEM",my_dirichlet_bc_dEfEf); CHKERRQ(ierr);
  ierr = mField.loop_finite_elements("STOCHASIC_PROBLEM","K_Inclusion",my_fe_k_rs_EfEf);  CHKERRQ(ierr);
  
  // d.
  ierr = mField.problem_basic_method_postProcess("STOCHASIC_PROBLEM",my_dirichlet_bc_dEfEf); CHKERRQ(ierr);

  // e.
  ierr = VecGhostUpdateBegin(ddF,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
  ierr = VecGhostUpdateEnd(ddF,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
  ierr = VecAssemblyBegin(ddF); CHKERRQ(ierr);
  ierr = VecAssemblyEnd(ddF); CHKERRQ(ierr);
  //  ierr = VecView(ddF,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);

  // f. solve the finite element equation to get answer for
  //     the second order partial derivative of nodal displacement
  ierr = KSPSolve(solver,ddF,ddD_EfEf); CHKERRQ(ierr);
  ierr = VecGhostUpdateBegin(ddD_EfEf,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = VecGhostUpdateEnd(ddD_EfEf,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  //  ierr = VecView(ddD,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);

  ierr = mField.set_other_global_VecCreateGhost("STOCHASIC_PROBLEM","DISPLACEMENT","DISP_rs_EfEf",ROW,ddD_EfEf,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
  

  // ==================================
  // 4.3.4 due to Ef Pf
  // ==================================
  // a.
  ierr = VecZeroEntries(ddF); CHKERRQ(ierr);
  ierr = VecGhostUpdateBegin(ddF,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = VecGhostUpdateEnd(ddF,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  
  // b. define/apply boundary condition
  DisplacementBCFEMethodPreAndPostProc my_dirichlet_bc_dEfPf(mField,"DISPLACEMENT",Aij,ddD_EfPf,ddF);
  K_rYoungPoissonFEMethod my_fe_k_rs_EfPf(mField,Aij,D,ddF,LAMBDA(young_modulus,poisson_ratio),MU(young_modulus,poisson_ratio),"DISP_r_Pf");

  // c. construct the first-order derivative of global stiffness matrix
  ierr = mField.problem_basic_method_preProcess("STOCHASIC_PROBLEM",my_dirichlet_bc_dEfPf); CHKERRQ(ierr);
  ierr = mField.loop_finite_elements("STOCHASIC_PROBLEM","K_Inclusion",my_fe_k_rs_EfPf);  CHKERRQ(ierr);

  // d.
  ierr = mField.problem_basic_method_postProcess("STOCHASIC_PROBLEM",my_dirichlet_bc_dEfPf); CHKERRQ(ierr);

  // e.
  ierr = VecGhostUpdateBegin(ddF,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
  ierr = VecGhostUpdateEnd(ddF,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
  ierr = VecAssemblyBegin(ddF); CHKERRQ(ierr);
  ierr = VecAssemblyEnd(ddF); CHKERRQ(ierr);
  //  ierr = VecView(ddF,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);

  // f. solve the finite element equation to get answer for
  //     the second order partial derivative of nodal displacement
  ierr = KSPSolve(solver,ddF,ddD_EfPf); CHKERRQ(ierr);
  ierr = VecGhostUpdateBegin(ddD_EfPf,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = VecGhostUpdateEnd(ddD_EfPf,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  //  ierr = VecView(ddD,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);

  ierr = mField.set_other_global_VecCreateGhost("STOCHASIC_PROBLEM","DISPLACEMENT","DISP_rs_EfPf",ROW,ddD_EfPf,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
  

  // ==================================
  // 4.4.1 due to Pf Em
  // ==================================


  // ==================================
  // 4.4.2 due to Pf Pm
  // ==================================


  // ==================================
  // 4.4.3 due to Pf Ef
  // ==================================


  cout<<"Solution phase Step 4.3.4 \n";
  // ==================================
  // 4.4.4 due to Pf Pf
  // ==================================
  // a.
  ierr = VecZeroEntries(ddF); CHKERRQ(ierr);
  ierr = VecGhostUpdateBegin(ddF,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = VecGhostUpdateEnd(ddF,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  cout<<"Solution phase substep a \n";
  // b. define/apply boundary condition
  DisplacementBCFEMethodPreAndPostProc my_dirichlet_bc_dPfPf(mField,"DISPLACEMENT",Aij,ddD_PfPf,ddF);
  K_rsPoissonFEMethod my_fe_k_rs_PfPf(mField,Aij,D,ddF,LAMBDA(young_modulus,poisson_ratio),MU(young_modulus,poisson_ratio),"DISP_r_Pf");

  // c. construct the first-order derivative of global stiffness matrix
  ierr = mField.problem_basic_method_preProcess("STOCHASIC_PROBLEM",my_dirichlet_bc_dPfPf); CHKERRQ(ierr);
  ierr = mField.loop_finite_elements("STOCHASIC_PROBLEM","K_Inclusion",my_fe_k_rs_PfPf);  CHKERRQ(ierr);

  // d.
  ierr = mField.problem_basic_method_postProcess("STOCHASIC_PROBLEM",my_dirichlet_bc_dPfPf); CHKERRQ(ierr);

  // e.
  ierr = VecGhostUpdateBegin(ddF,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
  ierr = VecGhostUpdateEnd(ddF,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
  ierr = VecAssemblyBegin(ddF); CHKERRQ(ierr);
  ierr = VecAssemblyEnd(ddF); CHKERRQ(ierr);
  //  ierr = VecView(ddF,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);

  // f. solve the finite element equation to get answer for
  //     the second order partial derivative of nodal displacement
  ierr = KSPSolve(solver,ddF,ddD_PfPf); CHKERRQ(ierr);
  ierr = VecGhostUpdateBegin(ddD_PfPf,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = VecGhostUpdateEnd(ddD_PfPf,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  //  ierr = VecView(ddD,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);

  ierr = mField.set_other_global_VecCreateGhost("STOCHASIC_PROBLEM","DISPLACEMENT","DISP_rs_PfPf",ROW,ddD_PfPf,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);



  // ===========================================================================
  //
  //  VI. OUTPUT
  //
  // ===========================================================================
  //Save data on mesh
  ierr = write_soltion(mField,"out.vtk","out_post_proc.vtk");   CHKERRQ(ierr);
  ierr = KSPDestroy(&solver); CHKERRQ(ierr);
  
   //Destroy matrices
  ierr = VecDestroy(&F); CHKERRQ(ierr);
  ierr = VecDestroy(&D); CHKERRQ(ierr);
  ierr = MatDestroy(&Aij); CHKERRQ(ierr);
  
  PetscSynchronizedFlush(PETSC_COMM_WORLD);
  
  PetscFinalize();
  
}

