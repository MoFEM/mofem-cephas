/* Copyright (C) 2015, Zahur Ullah (Zahur.Ullah AT glasgow.ac.uk)
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

#include <SurfacePressure.hpp>
#include <FEMethod_LowLevelStudent.hpp>
#include <FEMethod_UpLevelStudent.hpp>

#include <PostProcOnRefMesh.hpp>
#include <PostProcVertexMethod.hpp>
#include <PostProcDisplacementAndStrainOnRefindedMesh.hpp>

#include <ElasticFEMethod.hpp>  

using namespace ObosleteUsersModules;

#include <ElasticFEMethod_Matrix.hpp>
#include <ElasticFEMethod_Dmat_input.hpp>

#include "ElasticFE_RVELagrange_Disp.hpp"
#include "ElasticFE_RVELagrange_Disp_Multi_Rhs.hpp"
#include "RVEVolume.hpp"
#include "ElasticFE_RVELagrange_Homogenized_Stress_Disp.hpp"
#include <Calculate_RVE_Dmat_TransIso_Disp.hpp>

#include <Reliability_SurfacePressure.hpp>
#include <FE2_PostProcStressForReliability.hpp>

using namespace boost::numeric;


static char help[] = "...\n\n";

int main(int argc, char *argv[]) {
  ErrorCode rval;
  PetscErrorCode ierr;
  PetscInitialize(&argc,&argv,(char *)0,help);
  PetscBool flg = PETSC_TRUE;
  const char *option;

  //============================================================================
  //
  //  A. Micro (RVE) Problem
  //
  //============================================================================
  
  moab::Core mb_instance_RVE;
  Interface& moab_RVE = mb_instance_RVE;
  
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
  
  char mesh_file_name_RVE[255];
  ierr = PetscOptionsGetString(PETSC_NULL,"-my_file_RVE",mesh_file_name_RVE,255,&flg); CHKERRQ(ierr);
  if(flg != PETSC_TRUE) {
    SETERRQ(PETSC_COMM_SELF,1,"*** ERROR -my_file_RVE (MESH FILE NEEDED)");
  }
  
  char outName[PETSC_MAX_PATH_LEN]="out.vtk";
  ierr = PetscOptionsGetString(PETSC_NULL,"-my_out",outName,sizeof(outName),&flg); CHKERRQ(ierr);
  
  char outName2[PETSC_MAX_PATH_LEN]="out_post_proc.vtk";
  ierr = PetscOptionsGetString(PETSC_NULL,"-my_post_out",outName2,sizeof(outName2),&flg); CHKERRQ(ierr);

  PetscInt order_RVE;
  ierr = PetscOptionsGetInt(PETSC_NULL,"-my_order_RVE",&order_RVE,&flg); CHKERRQ(ierr);
  if(flg != PETSC_TRUE) {
    order_RVE = 1;
  }
  
  ParallelComm* pcomm_RVE = ParallelComm::get_pcomm(&moab_RVE,MYPCOMM_INDEX);
  if(pcomm_RVE == NULL) pcomm_RVE =  new ParallelComm(&moab_RVE,PETSC_COMM_SELF);
  
  
  /*****************************************************************************
   *
   * Transfer mesh data to MOAB database
   *
   ****************************************************************************/
  
  option = "";//"PARALLEL=BCAST;";//;DEBUG_IO";
  rval = moab_RVE.load_file(mesh_file_name_RVE, 0, option); CHKERR_PETSC(rval);
  
  //Create MoFEM (Joseph) database
  MoFEM::Core core_RVE(moab_RVE, PETSC_COMM_SELF);
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
  
  //    const clock_t begin_time = clock();
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
  
  EntityHandle meshset_Elastic, meshset_Trans_ISO;
  rval = moab_RVE.create_meshset(MESHSET_SET,meshset_Elastic); CHKERR_PETSC(rval);
  rval = moab_RVE.create_meshset(MESHSET_SET,meshset_Trans_ISO); CHKERR_PETSC(rval);
	
	///Getting No. of Fibres to be used for Potential Flow Problem
	int noOfFibres=0;
	for(_IT_CUBITMESHSETS_BY_BCDATA_TYPE_FOR_LOOP_(m_field_RVE,BLOCKSET|UNKNOWNCUBITNAME,it)) {
		
		std::size_t found=it->get_name().find("PotentialFlow");
		if (found==std::string::npos) continue;
		noOfFibres += 1;
	}
	cout<<"No. of Fibres for Potential Flow : "<<noOfFibres<<endl;
	
	vector<int> fibreList(noOfFibres,0);
	for (int aa=0; aa<noOfFibres; aa++) {
		fibreList[aa] = aa + 1;
	}
  
	Range RangeFibre[noOfFibres];
	EntityHandle fibre_meshset[noOfFibres];
	
	for (int ii=0; ii<noOfFibres; ii++) {
		ostringstream sss;
		sss << "POTENTIAL_ELEM" << ii+1;
		for(_IT_GET_FES_BY_NAME_FOR_LOOP_(m_field_RVE, sss.str().c_str() ,it)){
			RangeFibre[ii].insert(it->get_ent());
			rval = moab_RVE.create_meshset(MESHSET_SET,fibre_meshset[ii]); CHKERR_PETSC(rval);
			rval = moab_RVE.add_entities(fibre_meshset[ii],RangeFibre[ii]); CHKERR_PETSC(rval);
			rval = moab_RVE.unite_meshset(meshset_Trans_ISO,fibre_meshset[ii]); CHKERR_PETSC(rval);
		}
	}
  
  rval = moab_RVE.write_file("meshset_Trans_ISO.vtk","VTK","",&meshset_Trans_ISO,1); CHKERR_PETSC(rval);
  
	for(_IT_CUBITMESHSETS_BY_SET_TYPE_FOR_LOOP_(m_field_RVE,BLOCKSET,it)){
    
		if(it->get_name() == "MAT_ELASTIC_1") {
			Range TetsInBlock;
			rval = moab_RVE.get_entities_by_type(it->meshset, MBTET,TetsInBlock,true); CHKERR_PETSC(rval);
			Range block_rope_bit_level = intersect(LatestRefinedTets,TetsInBlock);
      
      cout<<"=============  TetsInBlock  "<< TetsInBlock.size() <<endl;
      
			rval = moab_RVE.add_entities(meshset_Elastic,block_rope_bit_level);CHKERR_PETSC(rval);
      
		}
	}
  ierr = m_field_RVE.seed_finite_elements(meshset_Elastic); CHKERRQ(ierr);
  
  Range prims_on_problem_bit_level_RVE;
	ierr = m_field_RVE.get_entities_by_type_and_ref_level(problem_bit_level_RVE,BitRefLevel().set(),MBPRISM,prims_on_problem_bit_level_RVE); CHKERRQ(ierr);
  //to create meshset from range
  EntityHandle meshset_prims_on_problem_bit_level_RVE;
  rval = moab_RVE.create_meshset(MESHSET_SET,meshset_prims_on_problem_bit_level_RVE); CHKERR_PETSC(rval);
	rval = moab_RVE.add_entities(meshset_prims_on_problem_bit_level_RVE,prims_on_problem_bit_level_RVE); CHKERR_PETSC(rval);
  ierr = m_field_RVE.seed_ref_level_MESHSET(meshset_prims_on_problem_bit_level_RVE,BitRefLevel().set()); CHKERRQ(ierr);

//  //set entitities bit level
//  BitRefLevel bit_level0_RVE;
//  bit_level0_RVE.set(0);
//  EntityHandle meshset_level0_RVE;
//  rval = moab_RVE.create_meshset(MESHSET_SET,meshset_level0_RVE); CHKERR_PETSC(rval);
//  ierr = m_field_RVE.seed_ref_level_3D(0,bit_level0_RVE); CHKERRQ(ierr);
  
  
  // ===========================================================================
  //
  // A.II. DEFINE PROBLEM
  //
  // ===========================================================================
  
  /*****************************************************************************
   *
   * Add field
   *  (1) Deterministic fields
   *
   ****************************************************************************/
  int field_rank=3;
  ierr = m_field_RVE.add_field("DISP_RVE",H1,3); CHKERRQ(ierr);
  ierr = m_field_RVE.add_field("Lagrange_mul_disp",H1,field_rank); CHKERRQ(ierr);
  
  
  //FE
  /*Here the Micro problem is divided in two parts (Matix(isotropic) and fibres(transverse-isotropic))
   Here we need to do [  Dmat_matrix=(1-w)*Dmat_matrix0  and Dmat_fibres=Dmat_fibres0]
   So degradaiton is happening only in matrix material and as we don't know about the acutal degredation
   in individual mechanical parameters [E, v] dut to limitted experimental data, so it is assumed that degradtion 
   happens in Dmat_material
  */
  /*****************************************************************************
   *
   * Create finite element for the defined fields
   *
   ****************************************************************************/
  ierr = m_field_RVE.add_finite_element("ELASTIC_FE_RVE"); CHKERRQ(ierr); //Matrix
  ierr = m_field_RVE.add_finite_element("TRAN_ISO_FE_RVE"); CHKERRQ(ierr); //Inclusion
  ierr = m_field_RVE.add_finite_element("Lagrange_FE"); CHKERRQ(ierr);
  
  //Define rows/cols and element data for ELASTIC_FE_RVE
  ierr = m_field_RVE.modify_finite_element_add_field_row("ELASTIC_FE_RVE","DISP_RVE"); CHKERRQ(ierr);
  ierr = m_field_RVE.modify_finite_element_add_field_col("ELASTIC_FE_RVE","DISP_RVE"); CHKERRQ(ierr);
  ierr = m_field_RVE.modify_finite_element_add_field_data("ELASTIC_FE_RVE","DISP_RVE"); CHKERRQ(ierr);
  
  //Define rows/cols and element data for TRAN_ISO_FE
  ierr = m_field_RVE.modify_finite_element_add_field_row("TRAN_ISO_FE_RVE","DISP_RVE"); CHKERRQ(ierr);
  ierr = m_field_RVE.modify_finite_element_add_field_col("TRAN_ISO_FE_RVE","DISP_RVE"); CHKERRQ(ierr);
  ierr = m_field_RVE.modify_finite_element_add_field_data("TRAN_ISO_FE_RVE","DISP_RVE"); CHKERRQ(ierr);
  ierr = m_field_RVE.modify_finite_element_add_field_data("TRAN_ISO_FE_RVE","POTENTIAL_FIELD"); CHKERRQ(ierr);

  //C row as Lagrange_mul_disp and col as DISPLACEMENT
  ierr = m_field_RVE.modify_finite_element_add_field_row("Lagrange_FE","Lagrange_mul_disp"); CHKERRQ(ierr);
  ierr = m_field_RVE.modify_finite_element_add_field_col("Lagrange_FE","DISP_RVE"); CHKERRQ(ierr);
  
  //CT col as Lagrange_mul_disp and row as DISPLACEMENT
  ierr = m_field_RVE.modify_finite_element_add_field_col("Lagrange_FE","Lagrange_mul_disp"); CHKERRQ(ierr);
  ierr = m_field_RVE.modify_finite_element_add_field_row("Lagrange_FE","DISP_RVE"); CHKERRQ(ierr);
  
  //data
  ierr = m_field_RVE.modify_finite_element_add_field_data("Lagrange_FE","Lagrange_mul_disp"); CHKERRQ(ierr);
  ierr = m_field_RVE.modify_finite_element_add_field_data("Lagrange_FE","DISP_RVE"); CHKERRQ(ierr);
  
  //define problems
  ierr = m_field_RVE.add_problem("ELASTIC_PROBLEM_RVE"); CHKERRQ(ierr);
  
  //set finite elements for problem
  ierr = m_field_RVE.modify_problem_add_finite_element("ELASTIC_PROBLEM_RVE","ELASTIC_FE_RVE"); CHKERRQ(ierr);
  ierr = m_field_RVE.modify_problem_add_finite_element("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE"); CHKERRQ(ierr);
  ierr = m_field_RVE.modify_problem_add_finite_element("ELASTIC_PROBLEM_RVE","Lagrange_FE"); CHKERRQ(ierr);
  
  
  //set refinment level for problem
  ierr = m_field_RVE.modify_problem_ref_level_add_bit("ELASTIC_PROBLEM_RVE",problem_bit_level_RVE); CHKERRQ(ierr);
  
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
  
  
  /*****************************************************************************
   *
   * Add finite elements entities
   *
   ****************************************************************************/
  
  ierr = m_field_RVE.add_ents_to_finite_element_by_TETs(meshset_Elastic,"ELASTIC_FE_RVE",true); CHKERRQ(ierr);
  ierr = m_field_RVE.add_ents_to_finite_element_by_TETs(meshset_Trans_ISO,"TRAN_ISO_FE_RVE",true); CHKERRQ(ierr);
  
  
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
   *   unstructured tetrahedral meshes. IJNME, 58(14). pp.2103-2130.
   ****************************************************************************/

  ierr = m_field_RVE.set_field_order(0,MBTET,"DISP_RVE",order_RVE); CHKERRQ(ierr);
  ierr = m_field_RVE.set_field_order(0,MBTRI,"DISP_RVE",order_RVE); CHKERRQ(ierr);
  ierr = m_field_RVE.set_field_order(0,MBEDGE,"DISP_RVE",order_RVE); CHKERRQ(ierr);
  ierr = m_field_RVE.set_field_order(0,MBVERTEX,"DISP_RVE",1); CHKERRQ(ierr);
  
  ierr = m_field_RVE.set_field_order(0,MBTRI,"Lagrange_mul_disp",order_RVE); CHKERRQ(ierr);
  ierr = m_field_RVE.set_field_order(0,MBEDGE,"Lagrange_mul_disp",order_RVE); CHKERRQ(ierr);
  ierr = m_field_RVE.set_field_order(0,MBVERTEX,"Lagrange_mul_disp",1); CHKERRQ(ierr);
  
  
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
  
  ParallelComm* pcomm = ParallelComm::get_pcomm(&moab_Macro,MYPCOMM_INDEX);
  if(pcomm == NULL) pcomm =  new ParallelComm(&moab_Macro,PETSC_COMM_WORLD);
  
  int rank;
  MPI_Comm_rank(PETSC_COMM_WORLD,&rank);
  
  /*****************************************************************************
   *
   * Read parameters from line command
   *
   ****************************************************************************/
  char mesh_file_name_Macro[255];
  ierr = PetscOptionsGetString(PETSC_NULL,"-my_file_Macro",mesh_file_name_Macro,255,&flg); CHKERRQ(ierr);
  if(flg != PETSC_TRUE) {
    SETERRQ(PETSC_COMM_SELF,1,"*** ERROR -my_file_Macro (MESH FILE NEEDED)");
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
  
  /*****************************************************************************
   *
   * Transfer mesh data to MOAB database
   *
   ****************************************************************************/
  option = "";//"PARALLEL=BCAST;";//;DEBUG_IO";
  rval = moab_Macro.load_file(mesh_file_name_Macro, 0, option); CHKERR_PETSC(rval);

  //Create MoFEM (Joseph) database
  MoFEM::Core core_Macro(moab_Macro);
  FieldInterface& m_field_Macro = core_Macro;

  //set entitities bit level
  BitRefLevel bit_level0_Macro;
  bit_level0_Macro.set(0);
  
  EntityHandle meshset_level0;
  rval = moab_Macro.create_meshset(MESHSET_SET,meshset_level0); CHKERR_PETSC(rval);
  ierr = m_field_Macro.seed_ref_level_3D(0,bit_level0_Macro); CHKERRQ(ierr);
  ierr = m_field_Macro.get_entities_by_ref_level(bit_level0_Macro,BitRefLevel().set(),meshset_level0); CHKERRQ(ierr);
  
  Range TetsInBlock_1st_Ply, TetsInBlock_2nd_Ply, TetsInBlock_3rd_Ply, TetsInBlock_4th_Ply, TetsInBlock_Reliability;
  for(_IT_CUBITMESHSETS_BY_SET_TYPE_FOR_LOOP_(m_field_Macro,BLOCKSET,it)){
    if(it->get_name() == "MAT_ELASTIC_First") {
      rval = moab_Macro.get_entities_by_type(it->meshset, MBTET,TetsInBlock_1st_Ply,true); CHKERR_PETSC(rval);
    }
    if(it->get_name() == "MAT_ELASTIC_Second") {
      rval = moab_Macro.get_entities_by_type(it->meshset, MBTET,TetsInBlock_2nd_Ply,true); CHKERR_PETSC(rval);
    }
    if(it->get_name() == "MAT_ELASTIC_Third") {
      rval = moab_Macro.get_entities_by_type(it->meshset, MBTET,TetsInBlock_3rd_Ply,true); CHKERR_PETSC(rval);
    }
    if(it->get_name() == "MAT_ELASTIC_Fourth") {
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

  //Fields
  ierr = m_field_Macro.add_field("DISP_MACRO",H1,3); CHKERRQ(ierr);

  //Problem
  ierr = m_field_Macro.add_problem("ELASTIC_PROBLEM_MACRO"); CHKERRQ(ierr);

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
  //meshset consisting all entities in mesh
  //EntityHandle root_set = moab_Macro.get_root_set();
  //add entities to field
  ierr = m_field_Macro.add_ents_to_field_by_TETs(0,"DISP_MACRO"); CHKERRQ(ierr);

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
  
  if(!(m_field_Macro.check_field("MESH_NODE_POSITIONS"))){
    ierr = m_field_Macro.add_field("MESH_NODE_POSITIONS",H1,3); CHKERRQ(ierr);
    ierr = m_field_Macro.add_ents_to_field_by_TETs(0,"MESH_NODE_POSITIONS"); CHKERRQ(ierr);
    ierr = m_field_Macro.set_field_order(0,MBTET,"MESH_NODE_POSITIONS",2); CHKERRQ(ierr);
    ierr = m_field_Macro.set_field_order(0,MBTRI,"MESH_NODE_POSITIONS",2); CHKERRQ(ierr);
    ierr = m_field_Macro.set_field_order(0,MBEDGE,"MESH_NODE_POSITIONS",2); CHKERRQ(ierr);
    ierr = m_field_Macro.set_field_order(0,MBVERTEX,"MESH_NODE_POSITIONS",1); CHKERRQ(ierr);
  }

  const char* args_new[] = {"ELASTIC_FE_MACRO_REL","ELASTIC_FE_MACRO", "ELASTIC_2nd_Ply","ELASTIC_3rd_Ply", "ELASTIC_4th_Ply"};
  vector<string> FE_Group_Name(args_new, args_new + 5);
  
  //Calculate Dmat for each Guass point here
  Calculate_RVE_Dmat_TransIso_Disp calculate_rve_dmat_TransIso(m_field_Macro);
  //ierr = calculate_rve_dmat_TransIso.addElasticElements("DISP_MACRO"); CHKERRQ(ierr);
  //ierr = m_field_Macro.modify_finite_element_add_field_data("ELASTIC_FE_MACRO","Wt"); CHKERRQ(ierr);
  for (int ife = 0; ife <=NO_Layers; ife++) {
    ierr = m_field_Macro.add_finite_element(FE_Group_Name[ife],MF_ZERO); CHKERRQ(ierr);
    ierr = m_field_Macro.modify_finite_element_add_field_row(FE_Group_Name[ife],"DISP_MACRO"); CHKERRQ(ierr);
    ierr = m_field_Macro.modify_finite_element_add_field_col(FE_Group_Name[ife],"DISP_MACRO"); CHKERRQ(ierr);
    ierr = m_field_Macro.modify_finite_element_add_field_data(FE_Group_Name[ife],"DISP_MACRO"); CHKERRQ(ierr);
    
    ierr = m_field_Macro.modify_finite_element_add_field_data(FE_Group_Name[ife],"Wt"); CHKERRQ(ierr);
    
    if(m_field_Macro.check_field("MESH_NODE_POSITIONS")) {
      ierr = m_field_Macro.modify_finite_element_add_field_data(FE_Group_Name[ife],"MESH_NODE_POSITIONS"); CHKERRQ(ierr);
    }
  }
  
  /*****************************************************************************
   *
   * Add finite elements entities
   *
   ****************************************************************************/
  ierr = m_field_Macro.add_ents_to_finite_element_by_TETs(TetsInBlock_1st_Ply, "ELASTIC_FE_MACRO"); CHKERRQ(ierr);
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
  * Set finite elements for problem
  *
  *****************************************************************************/
  ierr = m_field_Macro.modify_problem_add_finite_element("ELASTIC_PROBLEM_MACRO","ELASTIC_FE_MACRO_REL"); CHKERRQ(ierr);
  // First layer
  ierr = m_field_Macro.modify_problem_add_finite_element("ELASTIC_PROBLEM_MACRO","ELASTIC_FE_MACRO"); CHKERRQ(ierr);
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

  //****************************************************************************
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

  Projection10NodeCoordsOnField ent_method_material(m_field_Macro,"MESH_NODE_POSITIONS");
  ierr = m_field_Macro.loop_dofs("MESH_NODE_POSITIONS",ent_method_material); CHKERRQ(ierr);

  // ===========================================================================
  //
  //  B.V. MESH PARTITION
  //
  // ===========================================================================
  
  //partition
  ierr = m_field_Macro.partition_problem("ELASTIC_PROBLEM_MACRO"); CHKERRQ(ierr);
  ierr = m_field_Macro.partition_finite_elements("ELASTIC_PROBLEM_MACRO"); CHKERRQ(ierr);
  ierr = m_field_Macro.partition_ghost_dofs("ELASTIC_PROBLEM_MACRO"); CHKERRQ(ierr);

  // ===========================================================================
  //
  //  C. SOLVING FE EQUATION
  //
  // ===========================================================================
  
  Vec F; Vec D; Mat A;
  ierr = m_field_Macro.VecCreateGhost("ELASTIC_PROBLEM_MACRO",ROW,&F); CHKERRQ(ierr);
  ierr = VecDuplicate(F,&D); CHKERRQ(ierr);
  ierr = m_field_Macro.MatCreateMPIAIJWithArrays("ELASTIC_PROBLEM_MACRO",&A); CHKERRQ(ierr);


  ierr = VecZeroEntries(D); CHKERRQ(ierr);
  ierr = VecGhostUpdateBegin(D,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = VecGhostUpdateEnd(D,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = VecZeroEntries(F); CHKERRQ(ierr);
  ierr = VecGhostUpdateBegin(F,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = VecGhostUpdateEnd(F,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  
  //External forces, This vector is assemble only once (as this is not function of Dmat)
  boost::ptr_map<string,NeummanForcesSurface> neumann_forces;
  ierr = MetaNeummanForces::setNeumannFiniteElementOperators(m_field_Macro,neumann_forces,F,"DISP_MACRO"); CHKERRQ(ierr);
  boost::ptr_map<string,NeummanForcesSurface>::iterator mit = neumann_forces.begin();
  for(;mit!=neumann_forces.end();mit++) {
    ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO",mit->first,mit->second->getLoopFe()); CHKERRQ(ierr);
  }
  
//  ierr = VecView(F,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);
  
  ierr = calculate_rve_dmat_TransIso.setRVE_DmatRhsOperators(m_field_RVE, "DISP_MACRO","Wt"); CHKERRQ(ierr);
  Vec Fint;
  ierr = VecDuplicate(F,&Fint); CHKERRQ(ierr);

  PostPocOnRefinedMesh post_proc(m_field_Macro);
  ierr = post_proc.generateReferenceElementMesh(); CHKERRQ(ierr);
  ierr = post_proc.addFieldValuesPostProc("DISP_MACRO"); CHKERRQ(ierr);

  //read time series and do thermo elastci analysis
  SeriesRecorder *recorder_ptr;
  ierr = m_field_Macro.query_interface(recorder_ptr); CHKERRQ(ierr);
  int count=0;
  if( recorder_ptr->check_series("Wt_SERIES") ) {
    cout<<"============== Wt_SERIES exists =============== "<<endl;
    for(_IT_SERIES_STEPS_BY_NAME_FOR_LOOP_(recorder_ptr,"Wt_SERIES",sit)) {
      cout<<"\n**************************************"<<endl;
      cout<<"\n*      Time step: "<<count<<endl;
      cout<<"\n**************************************"<<endl;
      if(count==0) {
        PetscPrintf(PETSC_COMM_WORLD,"Process step %d\n",sit->get_step_number());
        ierr = recorder_ptr->load_series_data("Wt_SERIES",sit->get_step_number()); CHKERRQ(ierr);
        
        //Here we use ElasticFEMethod_Dmat_input, so will multiply Fint with -1
        DisplacementBCFEMethodPreAndPostProc my_dirichlet_bc(m_field_Macro,"DISP_MACRO",A,D,F);
        
        ElasticFEMethod_Dmat_input my_fe(m_field_Macro,A,D,Fint,0.0,0.0,calculate_rve_dmat_TransIso.commonData.Dmat_RVE,"DISP_MACRO");
        
        //preproc
        ierr = m_field_Macro.problem_basic_method_preProcess("ELASTIC_PROBLEM_MACRO",my_dirichlet_bc); CHKERRQ(ierr);
        
        // We need to assemble matrix A and internal force vector Fint at each
        // time step as these depends on Dmat, which will change at each time step
        ierr = MatZeroEntries(A); CHKERRQ(ierr);
        ierr = VecZeroEntries(Fint); CHKERRQ(ierr);
        ierr = VecGhostUpdateBegin(Fint,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
        ierr = VecGhostUpdateEnd(Fint,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
        
        ierr = VecZeroEntries(D); CHKERRQ(ierr);
        ierr = VecGhostUpdateBegin(D,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
        ierr = VecGhostUpdateEnd(D,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
        ierr = m_field_Macro.set_global_ghost_vector("ELASTIC_PROBLEM_MACRO",ROW,D,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
        
        // ------------------------
        // calculate Dmat for all Gauss points in the macro-mesh
        // ------------------------
        cout<<"\n\nSet 01"<<endl;
        ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_FE_MACRO",calculate_rve_dmat_TransIso.getLoopFeRhs()); CHKERRQ(ierr);
        ierr = VecScale(Fint,-1); CHKERRQ(ierr); //Multiply Fint with -1 (Fint=-Fint)
        ierr = VecAXPY(Fint,1,F); CHKERRQ(ierr); //Fint=Fint+F
        cout<<"\n\nSet 02"<<endl;
        //loop over macro elemnts to assemble A matrix and Fint vector
        ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_FE_MACRO",my_fe);  CHKERRQ(ierr);
        
        my_dirichlet_bc.snes_B=A;
        my_dirichlet_bc.snes_x = D;
        my_dirichlet_bc.snes_f = Fint;
        
        //postproc
        ierr = m_field_Macro.problem_basic_method_postProcess("ELASTIC_PROBLEM_MACRO",my_dirichlet_bc); CHKERRQ(ierr);
        
        //Solver
        KSP solver;
        ierr = KSPCreate(PETSC_COMM_WORLD,&solver); CHKERRQ(ierr);
        ierr = KSPSetOperators(solver,A,A); CHKERRQ(ierr);
        ierr = KSPSetFromOptions(solver); CHKERRQ(ierr);
        ierr = KSPSetUp(solver); CHKERRQ(ierr);
        
        ierr = KSPSolve(solver,Fint,D); CHKERRQ(ierr);
        ierr = VecGhostUpdateBegin(D,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
        ierr = VecGhostUpdateEnd(D,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
        
        
        //Save data on mesh
        ierr = m_field_Macro.set_local_ghost_vector("ELASTIC_PROBLEM_MACRO",ROW,D,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
        ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_FE_MACRO",post_proc); CHKERRQ(ierr);
        ostringstream o1;
        o1 << "FE2_out_" << sit->step_number << ".h5m";
        rval = post_proc.postProcMesh.write_file(o1.str().c_str(),"MOAB","PARALLEL=WRITE_PART"); CHKERR_PETSC(rval);
        
        if(count==100){
          //save the solution file for subsequent strain calculation analysis
          ierr = m_field_Macro.set_global_ghost_vector("ELASTIC_PROBLEM_MACRO",ROW,D,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
          rval = moab_Macro.write_file("FE2_solution_100.h5m"); CHKERR_PETSC(rval);
        }
        
        ierr = KSPDestroy(&solver); CHKERRQ(ierr);
        PetscPrintf(PETSC_COMM_WORLD,"End of step %d\n",sit->get_step_number());
      }
      count++;
    }
  }
  
  ierr = VecDestroy(&F); CHKERRQ(ierr);
  ierr = VecDestroy(&Fint); CHKERRQ(ierr);
  ierr = VecDestroy(&D); CHKERRQ(ierr);
  ierr = MatDestroy(&A); CHKERRQ(ierr);
  ierr = PetscFinalize(); CHKERRQ(ierr);



}


