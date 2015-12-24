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

#include <PostProcOnRefMesh.hpp>
#include <PostProcHookStresses.hpp>
#include <PostProcHookStressesLaminates.hpp>

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

using namespace boost::numeric;

ErrorCode rval;
PetscErrorCode ierr;

static char help[] = "...\n\n";

//const double young_modulus = 1;
//const double poisson_ratio = 0.0;

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
  
  PetscInitialize(&argc,&argv,(char *)0,help);
  
  int rank;
  MPI_Comm_rank(PETSC_COMM_WORLD,&rank);
  
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
  
  PetscBool flg = PETSC_TRUE;
  char mesh_file_name[255];
  
  ierr = PetscOptionsGetString(PETSC_NULL,"-my_file_RVE",mesh_file_name,255,&flg); CHKERRQ(ierr);
  if(flg != PETSC_TRUE) {
    SETERRQ(PETSC_COMM_SELF,1,"*** ERROR -my_file_RVE (MESH FILE NEEDED)");
  }
  
  PetscInt order_RVE;
  ierr = PetscOptionsGetInt(PETSC_NULL,"-my_order_RVE",&order_RVE,&flg); CHKERRQ(ierr);
  if(flg != PETSC_TRUE) {
    order_RVE = 1;
  }

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
  
  for(_IT_CUBITMESHSETS_BY_SET_TYPE_FOR_LOOP_(m_field_RVE,BLOCKSET,it)){
    
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
  int field_rank=3;
  ierr = m_field_RVE.add_field("DISP_RVE",H1,field_rank); CHKERRQ(ierr);
  ierr = m_field_RVE.add_field("Lagrange_mul_disp",H1,field_rank); CHKERRQ(ierr);
  
  // Stochastic fields for perturbation method
  for(int ii=0; ii < nders; ii++ ) {
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
  ierr = m_field_RVE.set_field_order(0,MBTET,"DISP_RVE",order_RVE); CHKERRQ(ierr);
  ierr = m_field_RVE.set_field_order(0,MBTRI,"DISP_RVE",order_RVE); CHKERRQ(ierr);
  ierr = m_field_RVE.set_field_order(0,MBEDGE,"DISP_RVE",order_RVE); CHKERRQ(ierr);
  ierr = m_field_RVE.set_field_order(0,MBVERTEX,"DISP_RVE",1); CHKERRQ(ierr);
  
  ierr = m_field_RVE.set_field_order(0,MBTRI,"Lagrange_mul_disp",order_RVE); CHKERRQ(ierr);
  ierr = m_field_RVE.set_field_order(0,MBEDGE,"Lagrange_mul_disp",order_RVE); CHKERRQ(ierr);
  ierr = m_field_RVE.set_field_order(0,MBVERTEX,"Lagrange_mul_disp",1); CHKERRQ(ierr);
  
  int order_st = order_RVE;
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
  
  //ParallelComm* pcomm_macro = ParallelComm::get_pcomm(&moab_Macro,MYPCOMM_INDEX);
  //if(pcomm == NULL) pcomm_macro =  new ParallelComm(&moab_Macro,PETSC_COMM_WORLD);
  
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
  
  /*****************************************************************************
   *
   * Transfer mesh data to MOAB database
   *
   ****************************************************************************/
  rval = moab_Macro.load_file(mesh_file_name, 0, option); CHKERR_PETSC(rval);
  
  //Create MoFEM (Joseph) database
  MoFEM::Core core_Macro(moab_Macro);
  FieldInterface& m_field_Macro = core_Macro;
  
  //ref meshset ref level 0
  ierr = m_field_Macro.seed_ref_level_3D(0,0); CHKERRQ(ierr);
  
  // stl::bitset see for more details
  BitRefLevel bit_level0_Macro;
  bit_level0_Macro.set(0);
  
  //EntityHandle meshset_level0;
  //rval = moab_Macro.create_meshset(MESHSET_SET,meshset_level0); CHKERR_PETSC(rval);
  ierr = m_field_Macro.seed_ref_level_3D(0,bit_level0_Macro); CHKERRQ(ierr);
  //ierr = m_field_Macro.get_entities_by_ref_level(bit_level0_Macro,BitRefLevel().set(),meshset_level0); CHKERRQ(ierr);
  
  
  EntityHandle meshset_macro_90deg, meshset_macro_pos25deg, meshset_macro_neg25deg;
  rval = moab_Macro.create_meshset(MESHSET_SET,meshset_macro_90deg); CHKERR_PETSC(rval);
  rval = moab_Macro.create_meshset(MESHSET_SET,meshset_macro_pos25deg); CHKERR_PETSC(rval);
  rval = moab_Macro.create_meshset(MESHSET_SET,meshset_macro_neg25deg); CHKERR_PETSC(rval);
  
  for(_IT_CUBITMESHSETS_BY_SET_TYPE_FOR_LOOP_(m_field_Macro,BLOCKSET,it)){
    if(it->get_name() == "MAT_ELASTIC_90deg") {
      Range TetsInBlock;
      rval = moab_Macro.get_entities_by_type(it->meshset, MBTET,TetsInBlock,true); CHKERR_PETSC(rval);
      Range block_rope_bit_level = intersect(LatestRefinedTets,TetsInBlock);
      rval = moab_Macro.add_entities(meshset_macro_90deg,block_rope_bit_level);CHKERR_PETSC(rval);
      
    }
  }
  ierr = m_field_Macro.seed_finite_elements(meshset_macro_90deg); CHKERRQ(ierr);
  
  for(_IT_CUBITMESHSETS_BY_SET_TYPE_FOR_LOOP_(m_field_Macro,BLOCKSET,it)){
    if(it->get_name() == "MAT_ELASTIC_Pos25deg") {
      Range TetsInBlock;
      rval = moab_Macro.get_entities_by_type(it->meshset, MBTET,TetsInBlock,true); CHKERR_PETSC(rval);
      Range block_rope_bit_level = intersect(LatestRefinedTets,TetsInBlock);
      rval = moab_Macro.add_entities(meshset_macro_pos25deg,block_rope_bit_level);CHKERR_PETSC(rval);
    }
  }
  ierr = m_field_Macro.seed_finite_elements(meshset_macro_pos25deg); CHKERRQ(ierr);
  
  for(_IT_CUBITMESHSETS_BY_SET_TYPE_FOR_LOOP_(m_field_Macro,BLOCKSET,it)){
    if(it->get_name() == "MAT_ELASTIC_Neg25deg") {
      Range TetsInBlock;
      rval = moab_Macro.get_entities_by_type(it->meshset, MBTET,TetsInBlock,true); CHKERR_PETSC(rval);
      Range block_rope_bit_level = intersect(LatestRefinedTets,TetsInBlock);
      rval = moab_Macro.add_entities(meshset_macro_neg25deg,block_rope_bit_level);CHKERR_PETSC(rval);
    }
  }
  ierr = m_field_Macro.seed_finite_elements(meshset_macro_neg25deg); CHKERRQ(ierr);
  
  //cout<<"===============================   TetsInBlock_90deg "<<TetsInBlock_90deg<<endl;
  //cout<<"===============================   TetsInBlock_pos25deg "<<TetsInBlock_pos25deg<<endl;
  //cout<<"===============================   TetsInBlock_neg25deg "<<TetsInBlock_neg25deg<<endl;
  
  
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
  for(int ii=0; ii < nvars; ii++ ) { // nders
    ostringstream ss_field;
    ss_field << "DISP_MACRO" << stochastic_fields[ii];
    ierr = m_field_Macro.add_field(ss_field.str().c_str(),H1,field_rank,MF_ZERO); CHKERRQ(ierr);
  }
  
  
  /*****************************************************************************
   *
   * Create finite element for the defined fields
   *
   ****************************************************************************/
  ierr = m_field_Macro.add_finite_element("ELASTIC_90deg"); CHKERRQ(ierr);
  ierr = m_field_Macro.add_finite_element("ELASTIC_pos25deg"); CHKERRQ(ierr);
  ierr = m_field_Macro.add_finite_element("ELASTIC_neg25deg"); CHKERRQ(ierr);
  ierr = m_field_Macro.add_finite_element("ELASTIC_0_90_post_process"); CHKERRQ(ierr);
  
  /*****************************************************************************
   *
   * set field data which finite element use
   *
   ****************************************************************************/
  //Define rows/cols and element data
  ierr = m_field_Macro.modify_finite_element_add_field_row("ELASTIC_90deg","DISP_MACRO"); CHKERRQ(ierr);
  ierr = m_field_Macro.modify_finite_element_add_field_col("ELASTIC_90deg","DISP_MACRO"); CHKERRQ(ierr);
  ierr = m_field_Macro.modify_finite_element_add_field_data("ELASTIC_90deg","DISP_MACRO"); CHKERRQ(ierr);
  ierr = m_field_Macro.modify_finite_element_add_field_data("ELASTIC_90deg","MESH_NODE_POSITIONS"); CHKERRQ(ierr);
  
  for(int ii=0; ii < nvars; ii++ ) {
    ostringstream ss_field;
    ss_field << "DISP_MACRO" << stochastic_fields[ii];
    ierr = m_field_Macro.modify_finite_element_add_field_data("ELASTIC_90deg",ss_field.str().c_str()); CHKERRQ(ierr);
  }
  
  //Define rows/cols and element data
  ierr = m_field_Macro.modify_finite_element_add_field_row("ELASTIC_pos25deg","DISP_MACRO"); CHKERRQ(ierr);
  ierr = m_field_Macro.modify_finite_element_add_field_col("ELASTIC_pos25deg","DISP_MACRO"); CHKERRQ(ierr);
  ierr = m_field_Macro.modify_finite_element_add_field_data("ELASTIC_pos25deg","DISP_MACRO"); CHKERRQ(ierr);
  ierr = m_field_Macro.modify_finite_element_add_field_data("ELASTIC_pos25deg","MESH_NODE_POSITIONS"); CHKERRQ(ierr);
  
  for(int ii=0; ii < nvars; ii++ ) {
    ostringstream ss_field;
    ss_field << "DISP_MACRO" << stochastic_fields[ii];
    ierr = m_field_Macro.modify_finite_element_add_field_data("ELASTIC_pos25deg",ss_field.str().c_str()); CHKERRQ(ierr);
  }
  
  //Define rows/cols and element data
  ierr = m_field_Macro.modify_finite_element_add_field_row("ELASTIC_neg25deg","DISP_MACRO"); CHKERRQ(ierr);
  ierr = m_field_Macro.modify_finite_element_add_field_col("ELASTIC_neg25deg","DISP_MACRO"); CHKERRQ(ierr);
  ierr = m_field_Macro.modify_finite_element_add_field_data("ELASTIC_neg25deg","DISP_MACRO"); CHKERRQ(ierr);
  ierr = m_field_Macro.modify_finite_element_add_field_data("ELASTIC_neg25deg","MESH_NODE_POSITIONS"); CHKERRQ(ierr);
  
  for(int ii=0; ii < nvars; ii++ ) {
    ostringstream ss_field;
    ss_field << "DISP_MACRO" << stochastic_fields[ii];
    ierr = m_field_Macro.modify_finite_element_add_field_data("ELASTIC_neg25deg",ss_field.str().c_str()); CHKERRQ(ierr);
  }
  
  //Define rows/cols and element data
  ierr = m_field_Macro.modify_finite_element_add_field_row("ELASTIC_0_90_post_process","DISP_MACRO"); CHKERRQ(ierr);
  ierr = m_field_Macro.modify_finite_element_add_field_col("ELASTIC_0_90_post_process","DISP_MACRO"); CHKERRQ(ierr);
  ierr = m_field_Macro.modify_finite_element_add_field_data("ELASTIC_0_90_post_process","DISP_MACRO"); CHKERRQ(ierr);
  ierr = m_field_Macro.modify_finite_element_add_field_data("ELASTIC_0_90_post_process","MESH_NODE_POSITIONS"); CHKERRQ(ierr);
  
  //define problems
  ierr = m_field_Macro.add_problem("ELASTIC_PROBLEM_MACRO"); CHKERRQ(ierr);
  
  //set finite elements for problem
  ierr = m_field_Macro.modify_problem_add_finite_element("ELASTIC_PROBLEM_MACRO","ELASTIC_90deg"); CHKERRQ(ierr);
  ierr = m_field_Macro.modify_problem_add_finite_element("ELASTIC_PROBLEM_MACRO","ELASTIC_pos25deg"); CHKERRQ(ierr);
  ierr = m_field_Macro.modify_problem_add_finite_element("ELASTIC_PROBLEM_MACRO","ELASTIC_neg25deg"); CHKERRQ(ierr);
  ierr = m_field_Macro.modify_problem_add_finite_element("ELASTIC_PROBLEM_MACRO","ELASTIC_0_90_post_process"); CHKERRQ(ierr);
  
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
  
  for(int ii=0; ii < nvars; ii++ ) { // nders
    ostringstream ss_field;
    ss_field << "DISP_MACRO" << stochastic_fields[ii];
    ierr = m_field_Macro.add_ents_to_field_by_TETs(0,ss_field.str().c_str()); CHKERRQ(ierr);
  }
  
  
  /*****************************************************************************
   *
   * Add finite elements entities
   *
   ****************************************************************************/
  ierr = m_field_Macro.add_ents_to_finite_element_by_TETs(meshset_macro_90deg, "ELASTIC_90deg",true); CHKERRQ(ierr);
  ierr = m_field_Macro.add_ents_to_finite_element_by_TETs(meshset_macro_pos25deg,"ELASTIC_pos25deg",true); CHKERRQ(ierr);
  ierr = m_field_Macro.add_ents_to_finite_element_by_TETs(meshset_macro_neg25deg,"ELASTIC_neg25deg",true); CHKERRQ(ierr);
  ierr = m_field_Macro.add_ents_to_finite_element_by_TETs(0, "ELASTIC_0_90_post_process"); CHKERRQ(ierr);
  
  
  /*****************************************************************************
   *
   * Set applied order
   * See reference for detals:
   *   Ainsworth M. and Coyle J. (2003) Hierarchic finite element bases on
   *      unstructured tetrahedral meshes. IJNME, 58(14). pp.2103-2130.
   ****************************************************************************/
  ierr = m_field_Macro.set_field_order(0,MBTET,"DISP_MACRO",order_Macro); CHKERRQ(ierr);
  ierr = m_field_Macro.set_field_order(0,MBTRI,"DISP_MACRO",order_Macro); CHKERRQ(ierr);
  ierr = m_field_Macro.set_field_order(0,MBEDGE,"DISP_MACRO",order_Macro); CHKERRQ(ierr);
  ierr = m_field_Macro.set_field_order(0,MBVERTEX,"DISP_MACRO",1); CHKERRQ(ierr);
  
  order_st = order_Macro;
  for(int ii=0; ii < nvars; ii++ ) { //nders
    ostringstream ss_field;
    ss_field << "DISP_MACRO" << stochastic_fields[ii];
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
  
  //mField.list_dofs_by_field_name("DISP",true);
  
  //print bcs
  ierr = m_field_Macro.print_cubit_displacement_set(); CHKERRQ(ierr);
  ierr = m_field_Macro.print_cubit_force_set(); CHKERRQ(ierr);
  //print block sets with materials
  ierr = m_field_Macro.print_cubit_materials_set(); CHKERRQ(ierr);
  
 
  // ===========================================================================
  //
  //  C. SOLUTION PHASE
  //
  // ===========================================================================
  
  
  FE2_Macro_Solver_Laminate Solve_FE2_Problem;
  ierr = Solve_FE2_Problem.Calculate_RVEDmat_PSFE(m_field_RVE, nvars, nders, stochastic_fields); CHKERRQ(ierr);
  ierr = Solve_FE2_Problem.Macro_FE_Solver_Laminate(m_field_Macro, nvars, nders, stochastic_fields); CHKERRQ(ierr);
  
  PetscFinalize();
  
}

