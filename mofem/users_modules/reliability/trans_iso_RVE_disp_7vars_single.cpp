/* Copyright (C) 2014,
 *   Zahur Ullah (Zahur.Ullah AT glasgow.ac.uk)
 *   Xiao-Yi Zhou (xiaoyi.zhou AT newcastle.ac.uk)
 * --------------------------------------------------------------
 * This routine conducts finite element implementation for stochastic multiscale
 * homogenization problem under linear displacement boundary condition by
 * perturbation technique for two-phase composites consisting of isotropic
 * material based matrix and transversely isotropic material based reinforcement
 * /inclusion/fibre, it thuse has seven independent material parameters, which
 * are Young's modulus and Poisson's ration of matrix, and Young's modulus in
 * p- and z-direction, Poisson's ratio in p- and z-direction, and shear modulus
 * in z-direction of reinforcement/fibre/inclusion, to characterize
 * the mechanical properties of this composite.
 *
 * HISTORY
 *
 * 2014.09.12 (first version)
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

#include <MoFEM.hpp>
using namespace MoFEM;

#include <Projection10NodeCoordsOnField.hpp>
#include <petsctime.h>

#include <FEMethod_LowLevelStudent.hpp>
#include <FEMethod_UpLevelStudent.hpp>

#include <PostProcVertexMethod.hpp>
#include <PostProcDisplacementAndStrainOnRefindedMesh.hpp>

#include "ElasticFEMethod.hpp"
#include "ElasticFEMethodTransIso.hpp"

using namespace ObosleteUsersModules;

#include "MaterialConstitutiveMatrix_FirstOrderDerivative.hpp"
#include "MaterialConstitutiveMatrix_SecondOrderDerivative.hpp"

#include "Trans_Iso_Rhs_r_PSFEM.hpp"
#include "Trans_Iso_Rhs_rs_PSFEM.hpp"

#include "ElasticFE_RVELagrange_Disp.hpp"
#include "ElasticFE_RVELagrange_Homogenized_Stress_Disp.hpp"
#include "RVEVolume.hpp"

ErrorCode rval;
PetscErrorCode ierr;

static char help[] = "...\n\n";

const char* args[] = {
    "_r_Em", "_r_Pm",
    "_r_NUp",     "_r_NUpz",      "_r_Ep",    "_r_Ez",    "_r_Gzp",     // 1st order
    "_rs_EmEm", "_rs_PmPm",                                             // 2nd order
    "_rs_NUpNUp", "_rs_NUpzNUpz", "_rs_EpEp", "_rs_EzEz", "_rs_GzpGzp",
    };

double nvars = 7;    // number of variables
double nders = 14;   // number of partial derivatives (firsr- and second- order)
vector<string> stochastic_fields(args, args + 14);

// =============================================================================
//
//  SETTING FOR OUTPUTS
//
// =============================================================================

PetscErrorCode write_soltion(FieldInterface &mField,const string out_file, const string out_ref_file) {
  PetscFunctionBegin;

  PostProcVertexMethod ent_method(mField.get_moab(),"DISPLACEMENT");
  ierr = mField.loop_dofs("STOCHASIC_PROBLEM","DISPLACEMENT",ROW,ent_method); CHKERRQ(ierr);

//  const char* args[] = {"_r_Em",    "_r_Pm",    "_r_Ef",    "_r_Pf",    // 1st order
//    "_rs_EmEm", "_rs_EmPm", "_rs_EmEf", "_rs_EmPf", // 2nd order
//    "_rs_PmPm", "_rs_PmEf", "_rs_PmPf",
//    "_rs_EfEf", "_rs_EfPf",
//    "_rs_PfPf"};

//  const char* args[] = {
//    "_r_NUp",     "_r_NUpz",      "_r_Ep",    "_r_Ez",    "_r_Gzp",     // 1st order
//    "_rs_NUpNUp", "_rs_NUpzNUpz", "_rs_EpEp", "_rs_EzEz", "_rs_GzpGzp", // 2nd order
//    };

//  vector<string> stochastic_fields(args, args+10);

  for(int ii=0; ii < nders; ii++ )
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
    ierr = mField.get_problem_finite_elements_entities("STOCHASIC_PROBLEM","ELASTIC",out_meshset); CHKERRQ(ierr);
    ierr = mField.get_problem_finite_elements_entities("STOCHASIC_PROBLEM","TRAN_ISOTROPIC_ELASTIC",out_meshset); CHKERRQ(ierr);
    ierr = mField.get_problem_finite_elements_entities("STOCHASIC_PROBLEM","INTERFACE",out_meshset); CHKERRQ(ierr);
    rval = mField.get_moab().write_file(out_file.c_str(),"VTK","",&out_meshset,1); CHKERR_PETSC(rval);
    rval = mField.get_moab().delete_entities(&out_meshset,1); CHKERR_PETSC(rval);
  }
  PetscFunctionReturn(0);
}


// =============================================================================
//
//  MAIN ROUTINE
//
// =============================================================================

int main(int argc, char *argv[]) {

  PetscInitialize(&argc,&argv,(char *)0,help);

  moab::Core mb_instance;
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
    order = 1;
  }

  char outName[PETSC_MAX_PATH_LEN]="out.vtk";
  ierr = PetscOptionsGetString(PETSC_NULL,"-my_out",outName,sizeof(outName),&flg); CHKERRQ(ierr);

  char outName2[PETSC_MAX_PATH_LEN]="out_post_proc.vtk";
  ierr = PetscOptionsGetString(PETSC_NULL,"-my_post_out",outName2,sizeof(outName2),&flg); CHKERRQ(ierr);

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
  MoFEM::Core core(moab);
  FieldInterface& mField = core;

  /*****************************************************************************
   *
   * Get fibre direction information from the potential-flow calculation
   *
   ****************************************************************************/
  Tag th_phi;
  //    double def_val  = 0;
  rval = moab.tag_get_handle("PHI",1,MB_TYPE_DOUBLE,th_phi); CHKERR_PETSC(rval);

  Tag th_meshset_info;
  int def_meshset_info[2] = {0,0};
  rval = moab.tag_get_handle("MESHSET_INFO",2,MB_TYPE_INTEGER,th_meshset_info,MB_TAG_CREAT|MB_TAG_SPARSE,&def_meshset_info);

  int meshset_data[2];
  EntityHandle root = moab.get_root_set();
  rval = moab.tag_get_data(th_meshset_info,&root,1,meshset_data); CHKERR_PETSC(rval);

  vector<BitRefLevel> bit_levels;
  bit_levels.push_back(BitRefLevel().set(meshset_data[0]-1));

  //    const clock_t begin_time = clock();
  ierr = mField.build_fields(); CHKERRQ(ierr);
  ierr = mField.build_finite_elements(); CHKERRQ(ierr);
  ierr = mField.build_adjacencies(bit_levels.back()); CHKERRQ(ierr);
  ierr = mField.build_problems(); CHKERRQ(ierr);

  //    std::cout << float( clock () - begin_time ) /  CLOCKS_PER_SEC<<endl;

	//Build FE
   /*****************************************************************************
   *
   * Select element into various mesh-set
   * meshset_level0: all element
   * meshset_Matrix: element set for matrix
   * meshset_Inclusion: element set for inclusion/fibre
   *
   ****************************************************************************/
  EntityHandle out_meshset;
  rval = moab.create_meshset(MESHSET_SET,out_meshset); CHKERR_PETSC(rval);
  //    ierr = mField.problem_get_FE("POTENTIAL_PROBLEM","POTENTIAL_ELEM",out_meshset); CHKERRQ(ierr);
  ierr = mField.get_entities_by_ref_level(bit_levels.back(),BitRefLevel().set(),out_meshset); CHKERRQ(ierr);
  Range LatestRefinedTets;
  rval = moab.get_entities_by_type(out_meshset, MBTET,LatestRefinedTets,true); CHKERR_PETSC(rval);

  Range LatestRefinedPrisms;
  rval = moab.get_entities_by_type(out_meshset, MBPRISM,LatestRefinedPrisms,true); CHKERR_PETSC(rval);

  cout<<"No of Prisms/Interfaces = "<<LatestRefinedPrisms.size()<<endl;

  BitRefLevel problem_bit_level = bit_levels.back();


  EntityHandle meshset_Elastic, meshset_Trans_ISO;
  rval = moab.create_meshset(MESHSET_SET,meshset_Elastic); CHKERR_PETSC(rval);
  rval = moab.create_meshset(MESHSET_SET,meshset_Trans_ISO); CHKERR_PETSC(rval);

  ///Getting No. of Fibres to be used for Potential Flow Problem
  int noOfFibres=0;
  for(_IT_CUBITMESHSETS_BY_BCDATA_TYPE_FOR_LOOP_(mField,BLOCKSET|UNKNOWNCUBITNAME,it)) {
     std::size_t found=it->get_name().find("PotentialFlow");
     if (found==std::string::npos) continue;
	noOfFibres += 1;
  }
  cout<<"No. of Fibres for Potential Flow : "<<noOfFibres<<endl;

  vector<int> fibreList(noOfFibres,0);
  for (int aa=0; aa<noOfFibres; aa++) {
     fibreList[aa] = aa + 1;
     }

  vector<Range> RangeFibre(noOfFibres);
  vector<EntityHandle> fibre_meshset(noOfFibres);

  for (int ii=0; ii<noOfFibres; ii++) {
    ostringstream sss;
    sss << "POTENTIAL_ELEM" << ii+1;
    for(_IT_GET_FES_BY_NAME_FOR_LOOP_(mField, sss.str().c_str() ,it)){
      RangeFibre[ii].insert(it->get_ent());
      rval = moab.create_meshset(MESHSET_SET,fibre_meshset[ii]); CHKERR_PETSC(rval);
      rval = moab.add_entities(fibre_meshset[ii],RangeFibre[ii]); CHKERR_PETSC(rval);
      rval = moab.unite_meshset(meshset_Trans_ISO,fibre_meshset[ii]); CHKERR_PETSC(rval);
    }
  }
  rval = moab.write_file("meshset_Trans_ISO.vtk","VTK","",&meshset_Trans_ISO,1); CHKERR_PETSC(rval);

  // Element set for matrix
  for(_IT_CUBITMESHSETS_BY_SET_TYPE_FOR_LOOP_(mField,BLOCKSET,it)){
    if(it->get_name() == "MAT_ELASTIC_1") {
      Range TetsInBlock;
      rval = moab.get_entities_by_type(it->meshset, MBTET,TetsInBlock,true); CHKERR_PETSC(rval);
      Range block_rope_bit_level = intersect(LatestRefinedTets,TetsInBlock);
      cout<<"=============  TetsInBlock  "<< TetsInBlock.size() <<endl;
      rval = moab.add_entities(meshset_Elastic,block_rope_bit_level);CHKERR_PETSC(rval);
    }
  }
  ierr = mField.seed_finite_elements(meshset_Elastic); CHKERRQ(ierr);

  Range prims_on_problem_bit_level;
  ierr = mField.get_entities_by_type_and_ref_level(problem_bit_level,BitRefLevel().set(),MBPRISM,prims_on_problem_bit_level); CHKERRQ(ierr);

  // to create meshset from range
  EntityHandle meshset_prims_on_problem_bit_level;
  rval = moab.create_meshset(MESHSET_SET,meshset_prims_on_problem_bit_level); CHKERR_PETSC(rval);
  rval = moab.add_entities(meshset_prims_on_problem_bit_level,prims_on_problem_bit_level); CHKERR_PETSC(rval);
  ierr = mField.seed_ref_level_MESHSET(meshset_prims_on_problem_bit_level,BitRefLevel().set()); CHKERRQ(ierr);

  // ===========================================================================
  //
  // II. DEFINE PROBLEM
  //
  // ===========================================================================

  //Fields
  int field_rank=3;
  ierr = mField.add_field("DISPLACEMENT",H1,field_rank); CHKERRQ(ierr);
  ierr = mField.add_field("Lagrange_mul_disp",H1,field_rank); CHKERRQ(ierr);

  /*****************************************************************************
   *
   * Add stochastic field
   * (total 14 field for 1st and 2nd order stochastic PSFEM)
   *
   ****************************************************************************/
//  const char* args[] = {
//    "_r_NUp",     "_r_NUpz",      "_r_Ep",    "_r_Ez",    "_r_Gzp",     // 1st order
//    "_rs_NUpNUp", "_rs_NUpzNUpz", "_rs_EpEp", "_rs_EzEz", "_rs_GzpGzp", // 2nd order
//    };

//  vector<string> stochastic_fields(args, args+10);

  for(int ii=0; ii < nders; ii++ )
  {
    ostringstream ss_field;
    ss_field << "DISP" << stochastic_fields[ii];
    //cout<<ss_field.str().c_str()<<endl;
    ierr = mField.add_field(ss_field.str().c_str(),H1,field_rank,MF_ZERO); CHKERRQ(ierr);
  }

//  ierr = mField.add_field("DISP_r_nup",H1,field_rank,MF_ZERO); CHKERRQ(ierr);
//  ierr = mField.add_field("DISP_rs_nup",H1,field_rank,MF_ZERO); CHKERRQ(ierr);

//  ierr = mField.add_field("DISP_r_Ep",H1,field_rank,MF_ZERO); CHKERRQ(ierr);
//  ierr = mField.add_field("DISP_rs_Ep",H1,field_rank,MF_ZERO); CHKERRQ(ierr);

  /*****************************************************************************
   *
   * Create finite element for the defined fields
   *
   ****************************************************************************/
  ierr = mField.add_finite_element("ELASTIC"); CHKERRQ(ierr);
  ierr = mField.add_finite_element("TRAN_ISOTROPIC_ELASTIC"); CHKERRQ(ierr);
  ierr = mField.add_finite_element("Lagrange_elem"); CHKERRQ(ierr);

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
  ierr = mField.modify_finite_element_add_field_row("ELASTIC","DISPLACEMENT"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_col("ELASTIC","DISPLACEMENT"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_data("ELASTIC","DISPLACEMENT"); CHKERRQ(ierr);

  //adding stochastic field to ELASTIC element
//  ierr = mField.modify_finite_element_add_field_data("ELASTIC","DISP_r_nup"); CHKERRQ(ierr);
//  ierr = mField.modify_finite_element_add_field_data("ELASTIC","DISP_r_Ep"); CHKERRQ(ierr);

  for(int ii=0; ii < nvars; ii++ )
  {
    ostringstream ss_field;
    ss_field << "DISP" << stochastic_fields[ii];
//    cout<<ss_field.str().c_str()<<endl;
    ierr = mField.modify_finite_element_add_field_data("ELASTIC",ss_field.str().c_str()); CHKERRQ(ierr);
  }

  //FE Transverse Isotropic
  ierr = mField.modify_finite_element_add_field_row("TRAN_ISOTROPIC_ELASTIC","DISPLACEMENT"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_col("TRAN_ISOTROPIC_ELASTIC","DISPLACEMENT"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_data("TRAN_ISOTROPIC_ELASTIC","DISPLACEMENT"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_data("TRAN_ISOTROPIC_ELASTIC","POTENTIAL_FIELD"); CHKERRQ(ierr);

  //adding stochastic field to TRAN_ISOTROPIC_ELASTIC element
//  ierr = mField.modify_finite_element_add_field_data("TRAN_ISOTROPIC_ELASTIC","DISP_r_nup"); CHKERRQ(ierr);
//  ierr = mField.modify_finite_element_add_field_data("TRAN_ISOTROPIC_ELASTIC","DISP_r_Ep"); CHKERRQ(ierr);

  for(int ii=0; ii < nvars; ii++ )
  {
    ostringstream ss_field;
    ss_field << "DISP" << stochastic_fields[ii];
//    cout<<ss_field.str().c_str()<<endl;
    ierr = mField.modify_finite_element_add_field_data("TRAN_ISOTROPIC_ELASTIC",ss_field.str().c_str()); CHKERRQ(ierr);
  }


  // Row and column and data for C and CT matrices
  //======================================================================================================
  //C row as Lagrange_mul_disp and col as DISPLACEMENT
  ierr = mField.modify_finite_element_add_field_row("Lagrange_elem","Lagrange_mul_disp"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_col("Lagrange_elem","DISPLACEMENT"); CHKERRQ(ierr);

  //CT col as Lagrange_mul_disp and row as DISPLACEMENT
  ierr = mField.modify_finite_element_add_field_col("Lagrange_elem","Lagrange_mul_disp"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_row("Lagrange_elem","DISPLACEMENT"); CHKERRQ(ierr);

  //As for stress we need both displacement and temprature (Lukasz)
  ierr = mField.modify_finite_element_add_field_data("Lagrange_elem","Lagrange_mul_disp"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_data("Lagrange_elem","DISPLACEMENT"); CHKERRQ(ierr);
  //======================================================================================================

  //define problems
  ierr = mField.add_problem("STOCHASIC_PROBLEM"); CHKERRQ(ierr);

  //set finite elements for problem
  ierr = mField.modify_problem_add_finite_element("STOCHASIC_PROBLEM","ELASTIC"); CHKERRQ(ierr);
  ierr = mField.modify_problem_add_finite_element("STOCHASIC_PROBLEM","TRAN_ISOTROPIC_ELASTIC"); CHKERRQ(ierr);
  ierr = mField.modify_problem_add_finite_element("STOCHASIC_PROBLEM","Lagrange_elem"); CHKERRQ(ierr);


  //set refinment level for problem
  ierr = mField.modify_problem_ref_level_add_bit("STOCHASIC_PROBLEM",problem_bit_level); CHKERRQ(ierr);

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
  ierr = mField.add_ents_to_field_by_TETs(0,"DISPLACEMENT"); CHKERRQ(ierr);

  for(int ii=0; ii < nders; ii++ )
  {
    ostringstream ss_field;
    ss_field << "DISP" << stochastic_fields[ii];
    ierr = mField.add_ents_to_field_by_TETs(0,ss_field.str().c_str()); CHKERRQ(ierr);
  }

  /*****************************************************************************
   *
   * Add finite elements entities
   *
   ****************************************************************************/
  ierr = mField.add_ents_to_finite_element_by_TETs(meshset_Elastic,"ELASTIC",true); CHKERRQ(ierr);
  ierr = mField.add_ents_to_finite_element_by_TETs(meshset_Trans_ISO,"TRAN_ISOTROPIC_ELASTIC",true); CHKERRQ(ierr);

  Range SurfacesFaces;
  ierr = mField.get_cubit_msId_entities_by_dimension(103,SIDESET,2,SurfacesFaces,true); CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"number of SIDESET 103 = %d\n",SurfacesFaces.size()); CHKERRQ(ierr);
  ierr = mField.add_ents_to_finite_element_by_TRIs(SurfacesFaces,"Lagrange_elem"); CHKERRQ(ierr);

  //to create meshset from range
  EntityHandle BoundFacesMeshset;
  rval = moab.create_meshset(MESHSET_SET,BoundFacesMeshset); CHKERR_PETSC(rval);
	rval = moab.add_entities(BoundFacesMeshset,SurfacesFaces); CHKERR_PETSC(rval);
  ierr = mField.seed_ref_level_MESHSET(BoundFacesMeshset,BitRefLevel().set()); CHKERRQ(ierr);

  ierr = mField.add_ents_to_field_by_TRIs(BoundFacesMeshset,"Lagrange_mul_disp",2); CHKERRQ(ierr);

  /*****************************************************************************
   *
   * Set applied order
   * See reference for detals:
   *   Ainsworth M. and Coyle J. (2003) Hierarchic finite element bases on
          unstructured tetrahedral meshes. IJNME, 58(14). pp.2103-2130.
   ****************************************************************************/
  //see Hierarchic Finite Element Bases on Unstructured Tetrahedral Meshes (Mark Ainsworth & Joe Coyle)
  ierr = mField.set_field_order(0,MBTET,"DISPLACEMENT",order); CHKERRQ(ierr);
  ierr = mField.set_field_order(0,MBTRI,"DISPLACEMENT",order); CHKERRQ(ierr);
  ierr = mField.set_field_order(0,MBEDGE,"DISPLACEMENT",order); CHKERRQ(ierr);
  ierr = mField.set_field_order(0,MBVERTEX,"DISPLACEMENT",1); CHKERRQ(ierr);

  ierr = mField.set_field_order(0,MBTRI,"Lagrange_mul_disp",order); CHKERRQ(ierr);
  ierr = mField.set_field_order(0,MBEDGE,"Lagrange_mul_disp",order); CHKERRQ(ierr);
  ierr = mField.set_field_order(0,MBVERTEX,"Lagrange_mul_disp",1); CHKERRQ(ierr);

  int order_st=order;
//  ierr = mField.set_field_order(0,MBTET,"DISP_r_nup",order_st); CHKERRQ(ierr);
//  ierr = mField.set_field_order(0,MBTRI,"DISP_r_nup",order_st); CHKERRQ(ierr);
//  ierr = mField.set_field_order(0,MBEDGE,"DISP_r_nup",order_st); CHKERRQ(ierr);
//  ierr = mField.set_field_order(0,MBVERTEX,"DISP_r_nup",1); CHKERRQ(ierr);
//
//  ierr = mField.set_field_order(0,MBTET,"DISP_rs_nup",order_st); CHKERRQ(ierr);
//  ierr = mField.set_field_order(0,MBTRI,"DISP_rs_nup",order_st); CHKERRQ(ierr);
//  ierr = mField.set_field_order(0,MBEDGE,"DISP_rs_nup",order_st); CHKERRQ(ierr);
//  ierr = mField.set_field_order(0,MBVERTEX,"DISP_rs_nup",1); CHKERRQ(ierr);

//  ierr = mField.set_field_order(0,MBTET,"DISP_r_Ep",order_st); CHKERRQ(ierr);
//  ierr = mField.set_field_order(0,MBTRI,"DISP_r_Ep",order_st); CHKERRQ(ierr);
//  ierr = mField.set_field_order(0,MBEDGE,"DISP_r_Ep",order_st); CHKERRQ(ierr);
//  ierr = mField.set_field_order(0,MBVERTEX,"DISP_r_Ep",1); CHKERRQ(ierr);
//
//  ierr = mField.set_field_order(0,MBTET,"DISP_rs_Ep",order_st); CHKERRQ(ierr);
//  ierr = mField.set_field_order(0,MBTRI,"DISP_rs_Ep",order_st); CHKERRQ(ierr);
//  ierr = mField.set_field_order(0,MBEDGE,"DISP_rs_Ep",order_st); CHKERRQ(ierr);
//  ierr = mField.set_field_order(0,MBVERTEX,"DISP_rs_Ep",1); CHKERRQ(ierr);

  for(int ii=0; ii < nders; ii++ )
  {
    ostringstream ss_field;
    ss_field << "DISP" << stochastic_fields[ii];
//    cout<<ss_field.str().c_str()<<endl;
    ierr = mField.set_field_order(0,MBTET,ss_field.str().c_str(),order_st); CHKERRQ(ierr);
    ierr = mField.set_field_order(0,MBTRI,ss_field.str().c_str(),order_st); CHKERRQ(ierr);
    ierr = mField.set_field_order(0,MBEDGE,ss_field.str().c_str(),order_st); CHKERRQ(ierr);
    ierr = mField.set_field_order(0,MBVERTEX,ss_field.str().c_str(),1); CHKERRQ(ierr);
  }


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
  ierr = mField.build_adjacencies(problem_bit_level); CHKERRQ(ierr);

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
  //create matrices
  Vec F,dF,ddF,D,dD,ddD;

  ierr = mField.VecCreateGhost("STOCHASIC_PROBLEM",ROW,&F); CHKERRQ(ierr);
  ierr = mField.VecCreateGhost("STOCHASIC_PROBLEM",ROW,&dF); CHKERRQ(ierr);
  ierr = mField.VecCreateGhost("STOCHASIC_PROBLEM",ROW,&ddF); CHKERRQ(ierr);

  ierr = mField.VecCreateGhost("STOCHASIC_PROBLEM",COL,&D); CHKERRQ(ierr);
  ierr = mField.VecCreateGhost("STOCHASIC_PROBLEM",COL,&dD); CHKERRQ(ierr);
  ierr = mField.VecCreateGhost("STOCHASIC_PROBLEM",COL,&ddD); CHKERRQ(ierr);

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


  for(_IT_GET_DOFS_FIELD_BY_NAME_FOR_LOOP_(mField,"POTENTIAL_FIELD",dof_ptr)) {
    if(dof_ptr->get_ent_type()!=MBVERTEX) continue;
    EntityHandle ent = dof_ptr->get_ent();
    double &fval = dof_ptr->get_FieldData();
    double phi;
    rval = moab.tag_get_data(th_phi,&ent,1,&phi); CHKERR_PETSC(rval);
    fval = phi;
  }

  //Assemble F and Aij
  double YoungModulus;
  double PoissonRatio;

  for(_IT_CUBITMESHSETS_BY_SET_TYPE_FOR_LOOP_(mField,BLOCKSET,it))
  {
    cout << endl << *it << endl;
    //Get block name
    string name = it->get_name();
    if (name.compare(0,11,"MAT_ELASTIC") == 0)
    {
      Mat_Elastic mydata;
      ierr = it->get_attribute_data_structure(mydata); CHKERRQ(ierr);
      cout << mydata;
      YoungModulus=mydata.data.Young;
      PoissonRatio=mydata.data.Poisson;
    }
  }
  // cout<<"the value"<<YoungModulus<<"\t"<<PoissonRatio<<endl;
  MyElasticFEMethod MyFE(mField,Aij,D,F,LAMBDA(YoungModulus,PoissonRatio),MU(YoungModulus,PoissonRatio));
  TranIsotropicFibreDirRotElasticFEMethod MyTIsotFE(mField,Aij,D,F);
  ElasticFE_RVELagrange_Disp MyFE_RVELagrange(mField,Aij,D,F,applied_strain,"DISPLACEMENT","Lagrange_mul_disp",field_rank);
  
  ierr = VecZeroEntries(F); CHKERRQ(ierr);
  ierr = VecGhostUpdateBegin(F,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = VecGhostUpdateEnd(F,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = MatZeroEntries(Aij); CHKERRQ(ierr);

  ierr = mField.loop_finite_elements("STOCHASIC_PROBLEM","ELASTIC",MyFE);  CHKERRQ(ierr);
  //PetscSynchronizedFlush(PETSC_COMM_WORLD,PETSC_STDOUT);
	ierr = mField.loop_finite_elements("STOCHASIC_PROBLEM","TRAN_ISOTROPIC_ELASTIC",MyTIsotFE);  CHKERRQ(ierr);
	PetscSynchronizedFlush(PETSC_COMM_WORLD,PETSC_STDOUT);
  ierr = mField.loop_finite_elements("STOCHASIC_PROBLEM","Lagrange_elem",MyFE_RVELagrange);  CHKERRQ(ierr);
  PetscSynchronizedFlush(PETSC_COMM_WORLD,PETSC_STDOUT);


  ierr = VecGhostUpdateBegin(F,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
  ierr = VecGhostUpdateEnd(F,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
  ierr = VecAssemblyBegin(F); CHKERRQ(ierr);
  ierr = VecAssemblyEnd(F); CHKERRQ(ierr);
  ierr = MatAssemblyBegin(Aij,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
  ierr = MatAssemblyEnd(Aij,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);

  //    //Matrix View
  //    MatView(Aij,PETSC_VIEWER_DRAW_WORLD);//PETSC_VIEWER_STDOUT_WORLD);
  //    std::string wait;
  //    std::cin >> wait;

  /*****************************************************************************
   *
   *  2. Get the volume of RVE
   *
   ****************************************************************************/
  double RVE_volume;    RVE_volume=0.0;  //RVE volume for full RVE We need this for stress calculation
  Vec RVE_volume_Vec;
  ierr = VecCreateMPI(PETSC_COMM_WORLD, 1, pcomm->size(), &RVE_volume_Vec);  CHKERRQ(ierr);
  ierr = VecZeroEntries(RVE_volume_Vec); CHKERRQ(ierr);

  RVEVolume MyRVEVol(mField,Aij,D,F,LAMBDA(YoungModulus,PoissonRatio),MU(YoungModulus,PoissonRatio), RVE_volume_Vec);
  RVEVolumeTrans MyRVEVolTrans(mField,Aij,D,F, RVE_volume_Vec);

  ierr = mField.loop_finite_elements("STOCHASIC_PROBLEM","ELASTIC",MyRVEVol);  CHKERRQ(ierr);
  ierr = mField.loop_finite_elements("STOCHASIC_PROBLEM","TRAN_ISOTROPIC_ELASTIC",MyRVEVolTrans);  CHKERRQ(ierr);

  ierr = VecView(RVE_volume_Vec,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);
  ierr = VecSum(RVE_volume_Vec, &RVE_volume);  CHKERRQ(ierr);
  cout<<"Final RVE_volume = "<< RVE_volume <<endl;
  cout<<"Actual RVE_volume = "<< 3*0.3*0.78<<endl;  //Lx=3, Ly=0.3; Lz=0.78

  /*****************************************************************************
   *
   *  3. SOLVE THE ZEROTH-ORDER FINITE ELEMENT EQUILIBRIUM EQUATION
   *     [K][U] = [F]
   *
   ****************************************************************************/

  //----------------------------------------------------------------------------
  // 3.1 Solving the equation to get nodal displacement
  //----------------------------------------------------------------------------
  //Solver
  KSP solver;
  ierr = KSPCreate(PETSC_COMM_WORLD,&solver); CHKERRQ(ierr);
  ierr = KSPSetOperators(solver,Aij,Aij); CHKERRQ(ierr);
  ierr = KSPSetFromOptions(solver); CHKERRQ(ierr);
  ierr = KSPSetUp(solver); CHKERRQ(ierr);

  ierr = KSPSolve(solver,F,D); CHKERRQ(ierr);
  ierr = VecGhostUpdateBegin(D,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = VecGhostUpdateEnd(D,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);

  //Save data on mesh
  //ierr = mField.set_global_VecCreateGhost("STOCHASIC_PROBLEM",ROW,D,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
  ierr = mField.set_global_ghost_vector("STOCHASIC_PROBLEM",ROW,D,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
  //ierr = VecView(F,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);

  //----------------------------------------------------------------------------
  // 3.2 Calculating zeroth-order homogenized stress using volume averaging theorem
  //----------------------------------------------------------------------------
  //create a vector for 6 components of homogenized stress
  Vec Stress_Homo;
  ierr = VecCreateMPI(PETSC_COMM_WORLD, 6, 6*pcomm->size(), &Stress_Homo);  CHKERRQ(ierr);
  ierr = VecZeroEntries(Stress_Homo); CHKERRQ(ierr);

  //    ierr = VecView(D,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);
  ElasticFE_RVELagrange_Homogenized_Stress_Disp MyFE_RVEHomoStressDisp(mField,Aij,D,F,&RVE_volume, applied_strain, Stress_Homo,"DISPLACEMENT","Lagrange_mul_disp",field_rank);
  ierr = mField.loop_finite_elements("STOCHASIC_PROBLEM","Lagrange_elem",MyFE_RVEHomoStressDisp);  CHKERRQ(ierr);

//  if(pcomm->rank()) cout<< " Stress_Homo =  "<<endl;
//  ierr = VecView(Stress_Homo,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);
// write the result in file
  ofstream TheFile;

  if (applied_strain[0]==1) {
	TheFile.open("//mnt//home//Dropbox//DURACOMP_Cal//009_MoFEM//03_EEP_Geometry_Uncertainties//Results//Result_TI_1.txt",ofstream::out);
  }
  else if (applied_strain[1]==1) {
	TheFile.open("//mnt//home//Dropbox//DURACOMP_Cal//009_MoFEM//03_EEP_Geometry_Uncertainties//Results//Result_TI_2.txt",ofstream::out);
  }
  else if (applied_strain[2]==1) {
	TheFile.open("//mnt//home//Dropbox//DURACOMP_Cal//009_MoFEM//03_EEP_Geometry_Uncertainties//Results//Result_TI_3.txt",ofstream::out);
  }
  else if (applied_strain[3]==1) {
	TheFile.open("//mnt//home//Dropbox//DURACOMP_Cal//009_MoFEM//03_EEP_Geometry_Uncertainties//Results//Result_TI_4.txt",ofstream::out);
  }
  else if (applied_strain[4]==1) {
	TheFile.open("//mnt//home//Dropbox//DURACOMP_Cal//009_MoFEM//03_EEP_Geometry_Uncertainties//Results//Result_TI_5.txt",ofstream::out);
  }
  else if (applied_strain[5]==1) {
	TheFile.open("//mnt//home//Dropbox//DURACOMP_Cal//009_MoFEM//03_EEP_Geometry_Uncertainties//Results//Result_TI_6.txt",ofstream::out);
  }

  if(pcomm->rank()==0){
    PetscScalar    *avec;
    VecGetArray(Stress_Homo, &avec);

    cout<< "\nStress_Homo =\n";
    for(int ii=0; ii<6; ii++){
      cout.precision(15);
      cout <<*avec<<endl;
	  TheFile<<setprecision(15)<<*avec<<'\n';
      avec++;
    }
  }
    cout<< "\n\n";

   /*****************************************************************************
   *  4. SOLVE THE FIRST-ORDER AND SECOND-ORDER FINITE ELEMENT EQUILIBRIUM EQUATION
   *     1st order-[K][U_r] = -[K_r][U}
   *     2nd order-[K][U_rs] = -[K_rs][U]-2[K_r][U_s]
   *
   ****************************************************************************/

  //----------------------------------------------------------------------------
  // a. Solving the equation to get nodal displacement
  //----------------------------------------------------------------------------

  for(int ii=0; ii < nders; ii++){
    //
    ierr = VecZeroEntries(dF); CHKERRQ(ierr);
    ierr = VecGhostUpdateBegin(dF,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
    ierr = VecGhostUpdateEnd(dF,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
    //
    ierr = VecZeroEntries(ddF); CHKERRQ(ierr);
    ierr = VecGhostUpdateBegin(ddF,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
    ierr = VecGhostUpdateEnd(ddF,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
    //
    if (ii == 0){ // due to Young's modulus of matrix - isotropic
      Trans_Iso_Rhs_r_PSFEM my_fe_k_r_Em(mField,Aij,D,dF,"DISPLACEMENT","Young","isotropic","matrix");
      ierr = mField.loop_finite_elements("STOCHASIC_PROBLEM","ELASTIC",my_fe_k_r_Em);  CHKERRQ(ierr);
      }
    else if (ii == 1){ // due to Poisson's ratio in p-direction of fibre
      cout<<"Poisson \t";
      Trans_Iso_Rhs_r_PSFEM my_fe_k_r_Pm(mField,Aij,D,dF,"DISPLACEMENT","Poisson","isotropic","matrix");
      ierr = mField.loop_finite_elements("STOCHASIC_PROBLEM","ELASTIC",my_fe_k_r_Pm);  CHKERRQ(ierr);
      }
    else if (ii == 2){ // due to Poisson's ratio in p-direction of fibre
      Trans_Iso_Rhs_r_PSFEM my_fe_k_r_NUp(mField,Aij,D,dF,"DISPLACEMENT","PoissonP", "transversely_isotropic", "reinforcement");
      ierr = mField.loop_finite_elements("STOCHASIC_PROBLEM","TRAN_ISOTROPIC_ELASTIC",my_fe_k_r_NUp);  CHKERRQ(ierr);
      }
    else if (ii == 3){ // due to Poisson's ratio in z-direction of fibre
      Trans_Iso_Rhs_r_PSFEM my_fe_k_r_NUpz(mField,Aij,D,dF,"DISPLACEMENT","PoissonZ", "transversely_isotropic", "reinforcement");
      ierr = mField.loop_finite_elements("STOCHASIC_PROBLEM","TRAN_ISOTROPIC_ELASTIC",my_fe_k_r_NUpz);  CHKERRQ(ierr);
      }
    else if (ii == 4){ // due to Young's modulus in p-direction of fibre
      Trans_Iso_Rhs_r_PSFEM my_fe_k_r_Ep(mField,Aij,D,dF,"DISPLACEMENT","YoungP", "transversely_isotropic", "reinforcement");
      ierr = mField.loop_finite_elements("STOCHASIC_PROBLEM","TRAN_ISOTROPIC_ELASTIC",my_fe_k_r_Ep);  CHKERRQ(ierr);
      }
    else if (ii == 5){ // due to Young's modulus in z-direction of fibre
      Trans_Iso_Rhs_r_PSFEM my_fe_k_r_Ez(mField,Aij,D,dF,"DISPLACEMENT","YoungZ", "transversely_isotropic", "reinforcement");
      ierr = mField.loop_finite_elements("STOCHASIC_PROBLEM","TRAN_ISOTROPIC_ELASTIC",my_fe_k_r_Ez);  CHKERRQ(ierr);
      }
    else if (ii == 6){ // due to shear modulus in z-direction of fibre
      Trans_Iso_Rhs_r_PSFEM my_fe_k_r_Gzp(mField,Aij,D,dF,"DISPLACEMENT","ShearZP", "transversely_isotropic", "reinforcement");
      ierr = mField.loop_finite_elements("STOCHASIC_PROBLEM","TRAN_ISOTROPIC_ELASTIC",my_fe_k_r_Gzp);  CHKERRQ(ierr);
      }
    else if (ii == 7){ // 2nd order derivative due to Young's modulus of matrix - isotropic
      Trans_Iso_Rhs_rs_PSFEM my_fe_k_rs_EmEm(mField,Aij,D,ddF,"DISPLACEMENT","DISP_r_Em","Young","Young", "isotropic", "matrix");
      ierr = mField.loop_finite_elements("STOCHASIC_PROBLEM","ELASTIC",my_fe_k_rs_EmEm);  CHKERRQ(ierr);
      }
    else if (ii == 8){ // 2nd order derivative due to Poisson's ratio of matrix - isotropic
      Trans_Iso_Rhs_rs_PSFEM my_fe_k_rs_PmPm(mField,Aij,D,ddF,"DISPLACEMENT","DISP_r_Pm","Poisson","Poisson", "isotropic", "matrix");
      ierr = mField.loop_finite_elements("STOCHASIC_PROBLEM","ELASTIC",my_fe_k_rs_PmPm);  CHKERRQ(ierr);
      }
    else if (ii == 9){ // 2nd order derivative due to Poisson's ratio in p-direction of fibre
      Trans_Iso_Rhs_rs_PSFEM my_fe_k_rs_NUpNUp(mField,Aij,D,ddF,"DISPLACEMENT","DISP_r_NUp","PoissonP","PoissonP", "transversely_isotropic", "reinforcement");
      ierr = mField.loop_finite_elements("STOCHASIC_PROBLEM","TRAN_ISOTROPIC_ELASTIC",my_fe_k_rs_NUpNUp);  CHKERRQ(ierr);
      }
    else if (ii == 10){ // 2nd order derivative due to Poisson's ratio in z-direction of fibre
      Trans_Iso_Rhs_rs_PSFEM my_fe_k_rs_NUpzNUpz(mField,Aij,D,ddF,"DISPLACEMENT","DISP_r_NUpz","PoissonZ","PoissonZ", "transversely_isotropic", "reinforcement");
      ierr = mField.loop_finite_elements("STOCHASIC_PROBLEM", "TRAN_ISOTROPIC_ELASTIC", my_fe_k_rs_NUpzNUpz);  CHKERRQ(ierr);
      }
    else if (ii == 11){ // 2nd order derivative due to Young's modulus in p-direction of fibre
      Trans_Iso_Rhs_rs_PSFEM my_fe_k_rs_EpEp(mField,Aij,D,ddF,"DISPLACEMENT","DISP_r_Ep","YoungP","YoungP", "transversely_isotropic", "reinforcement");
      ierr = mField.loop_finite_elements("STOCHASIC_PROBLEM","TRAN_ISOTROPIC_ELASTIC",my_fe_k_rs_EpEp);  CHKERRQ(ierr);
      }
    else if (ii == 12){ // 2nd order derivative due to Young's modulus in z-direction of fibre
      Trans_Iso_Rhs_rs_PSFEM my_fe_k_rs_EzEz(mField,Aij,D,ddF,"DISPLACEMENT","DISP_r_Ez","YoungZ","YoungZ", "transversely_isotropic", "reinforcement");
      ierr = mField.loop_finite_elements("STOCHASIC_PROBLEM","TRAN_ISOTROPIC_ELASTIC",my_fe_k_rs_EzEz);  CHKERRQ(ierr);
      }
    else if (ii == 13){ // 2nd order derivative due to shear modulus in z-direction of fibre
      Trans_Iso_Rhs_rs_PSFEM my_fe_k_rs_GzpGzp(mField,Aij,D,ddF,"DISPLACEMENT","DISP_r_Gzp","ShearZP","ShearZP", "transversely_isotropic", "reinforcement");
      ierr = mField.loop_finite_elements("STOCHASIC_PROBLEM","TRAN_ISOTROPIC_ELASTIC",my_fe_k_rs_GzpGzp);  CHKERRQ(ierr);
      }

    ostringstream ss_field;
    ss_field << "DISP" << stochastic_fields[ii];
    if (ii < nvars){ // solution for first-order problem
       ierr = VecGhostUpdateBegin(dF,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
       ierr = VecGhostUpdateEnd(dF,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
       ierr = VecAssemblyBegin(dF); CHKERRQ(ierr);
       ierr = VecAssemblyEnd(dF); CHKERRQ(ierr);

       ierr = KSPSolve(solver,dF,dD); CHKERRQ(ierr);
       ierr = VecGhostUpdateBegin(dD,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
       ierr = VecGhostUpdateEnd(dD,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
       // ierr = VecView(dD,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);
       ierr = mField.set_other_global_ghost_vector("STOCHASIC_PROBLEM","DISPLACEMENT",ss_field.str().c_str(),ROW,dD,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
       ierr = mField.set_other_global_ghost_vector("STOCHASIC_PROBLEM","Lagrange_mul_disp","Lagrange_mul_disp",ROW,dD,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
       }
    else { // solution for second-order problem
       ierr = VecGhostUpdateBegin(ddF,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
       ierr = VecGhostUpdateEnd(ddF,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
       ierr = VecAssemblyBegin(ddF); CHKERRQ(ierr);
       ierr = VecAssemblyEnd(ddF); CHKERRQ(ierr);

       ierr = KSPSolve(solver,ddF,ddD); CHKERRQ(ierr);
       ierr = VecGhostUpdateBegin(ddD,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
       ierr = VecGhostUpdateEnd(ddD,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
       // ierr = VecView(dD,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);
       ierr = mField.set_other_global_ghost_vector("STOCHASIC_PROBLEM","DISPLACEMENT",ss_field.str().c_str(),ROW,ddD,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
       ierr = mField.set_other_global_ghost_vector("STOCHASIC_PROBLEM","Lagrange_mul_disp","Lagrange_mul_disp",ROW,ddD,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
       }

  //----------------------------------------------------------------------------
  // b. Calculating first-order homogenized stress
  //----------------------------------------------------------------------------
  Vec Stress_Homo_r;
  ierr = VecCreateMPI(PETSC_COMM_WORLD, 6, 6*pcomm->size(), &Stress_Homo_r);  CHKERRQ(ierr);
  ierr = VecZeroEntries(Stress_Homo_r); CHKERRQ(ierr);

  ElasticFE_RVELagrange_Homogenized_Stress_Disp MyFE_RVEHomoStressDisp_r(mField,Aij,D,F,&RVE_volume, applied_strain, Stress_Homo_r,"DISPLACEMENT","Lagrange_mul_disp",field_rank);
    ierr = mField.loop_finite_elements("STOCHASIC_PROBLEM","Lagrange_elem",MyFE_RVEHomoStressDisp_r);  CHKERRQ(ierr);


    if(pcomm->rank()==0){
      PetscScalar    *avec_r;
      VecGetArray(Stress_Homo_r, &avec_r);

      cout<<"\n"<<ss_field.str().c_str()<<" =\n"<<endl;
      //cout<< "\n"<<ss_field<<" = \n\n";
      for(int ii=0; ii<6; ii++){
        cout.precision(15);
        cout<<*avec_r<<endl;
		// write result to output file
		TheFile<<setprecision(15)<<*avec_r<<'\n';

        avec_r++;
      }
    }
    cout<< "\n\n";
  }

 TheFile.close();

  // ===========================================================================
  //
  //  VII. OUTPUT
  //
  // ===========================================================================
  // Save data on mesh
  ierr = write_soltion(mField,"out.vtk","out_post_proc.vtk");   CHKERRQ(ierr);

  /* ofstream TheFile;
  TheFile.open("Result.txt",ofstream::out);
  if(pcomm->rank()==0){
    PetscScalar    *avec;
    VecGetArray(Stress_Homo, &avec);

    for(int ii=0; ii<6; ii++){
      TheFile<<setprecision(15)<<*avec<<'\n';
      avec++;
    }
  }
  TheFile.close();*/

  // ===========================================================================
  //
  //  VIII. FINISH
  //
  // ===========================================================================
  // Destroy matrices
  ierr = VecDestroy(&F); CHKERRQ(ierr);
  ierr = VecDestroy(&dF); CHKERRQ(ierr);
  ierr = VecDestroy(&ddF); CHKERRQ(ierr);


  ierr = VecDestroy(&D); CHKERRQ(ierr);
  ierr = VecDestroy(&dD); CHKERRQ(ierr);
  ierr = VecDestroy(&ddD); CHKERRQ(ierr);

  ierr = MatDestroy(&Aij); CHKERRQ(ierr);
  ierr = KSPDestroy(&solver); CHKERRQ(ierr);
  ierr = VecDestroy(&RVE_volume_Vec); CHKERRQ(ierr);
  ierr = VecDestroy(&Stress_Homo); CHKERRQ(ierr);


  ierr = PetscTime(&v2);CHKERRQ(ierr);
  ierr = PetscGetCPUTime(&t2);CHKERRQ(ierr);

  PetscSynchronizedPrintf(PETSC_COMM_WORLD,"Total Rank %d Time = %f CPU Time = %f\n",pcomm->rank(),v2-v1,t2-t1);
  PetscSynchronizedFlush(PETSC_COMM_WORLD,PETSC_STDOUT);

  PetscFinalize();

}
