/* Copyright (C) 2014,
 *   Zahur Ullah (Zahur.Ullah AT glasgow.ac.uk)
 *   Xiao-Yi Zhou (xiaoyi.zhou AT newcastle.ac.uk)
 * --------------------------------------------------------------
 * This routine conducts finite element implementation for stochastic multiscale
 * homogenization problem under periodic displacement and anti-periodic boundary
 * condition by perturbation technique for two-phase composites consisting of
 * isotropic material based matrix and transversely isotropic material based
 * reinforcement/inclusion/fibre, it thuse has seven independent material
 * parameters, which are Young's modulus and Poisson's ration of matrix, and
 * Young's modulus in p- and z-direction, Poisson's ratio in p- and z-direction,
 * and shear modulus in z-direction of reinforcement/fibre/inclusion, to
 * characterize the mechanical properties of this composite.
 *
 * HISTORY
 *
 * 2014.09.14 (first version)
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

#include "ElasticFE_RVELagrange_Periodic.hpp"
#include "ElasticFE_RVELagrange_Periodic_Multi_Rhs.hpp"
#include "ElasticFE_RVELagrange_RigidBodyTranslation.hpp"
#include "ElasticFE_RVELagrange_Homogenized_Stress_Periodic.hpp"
#include "RVEVolume.hpp"

ErrorCode rval;
PetscErrorCode ierr;

static char help[] = "...\n\n";

const char* args[] = {
    "_r_Em", "_r_Pm",
    "_r_NUp",     "_r_NUpz",      "_r_Ep",    "_r_Ez",    "_r_Gzp",     // 1st order
    "_rs_EmEm", "_rs_PmPm",                                             // 2nd order
  "_rs_NUpNUp", "_rs_NUpzNUpz", "_rs_EpEp", "_rs_EzEz", "_rs_GzpGzp"};
//    "_rs_NUpNUpz", "_rs_NUpEp",  "_rs_NUpEz",  "_rs_NUpGzp",  "_rs_NUpEm",  "_rs_NUpPm", // mix
//                   "_rs_NUpzEp", "_rs_NUpzEz", "_rs_NUpzGzp", "_rs_NUpzEm", "_rs_NUpzPm",
//                                 "_rs_EpEz",   "_rs_EpGzp",   "_rs_EpEm",   "_rs_EpPm",
//                                               "_rs_EzGzp",   "_rs_EzEm",   "_rs_EzPm",
//                                                              "_rs_GzpEm",  "_rs_GzpPm",
//                                                                             "_rs_EmPm",
//    };

double nvars = 7;    // number of variables
double nders = 14;   // number of partial derivatives (firsr- and second- order)
vector<string> stochastic_fields(args, args + 14);

//==============================================================================
//
// Define class and multindex container to store data for traiangles on the
// boundary of the RVE (it cannot be defined within main)
//
//==============================================================================

struct Face_CenPos_Handle
{
    double xcoord, ycoord, zcoord;
    const EntityHandle  Tri_Hand;
    Face_CenPos_Handle(double _xcoord, double _ycoord,  double _zcoord,  const EntityHandle _Tri_Hand):xcoord(_xcoord),
    ycoord(_ycoord), zcoord(_zcoord), Tri_Hand(_Tri_Hand) {}

};


struct xcoord_tag {}; //tags to used in the multindex container
struct ycoord_tag {};
struct zcoord_tag {};
struct Tri_Hand_tag {};
struct Composite_xyzcoord {};

typedef multi_index_container<
Face_CenPos_Handle,
indexed_by<
ordered_non_unique<
tag<xcoord_tag>, member<Face_CenPos_Handle,double,&Face_CenPos_Handle::xcoord> >,

ordered_non_unique<
tag<ycoord_tag>, member<Face_CenPos_Handle,double,&Face_CenPos_Handle::ycoord> >,

ordered_non_unique<
tag<zcoord_tag>, member<Face_CenPos_Handle,double,&Face_CenPos_Handle::zcoord> >,

ordered_unique<
tag<Tri_Hand_tag>, member<Face_CenPos_Handle,const EntityHandle,&Face_CenPos_Handle::Tri_Hand> >,

ordered_unique<
tag<Composite_xyzcoord>,
composite_key<
Face_CenPos_Handle,
member<Face_CenPos_Handle,double,&Face_CenPos_Handle::xcoord>,
member<Face_CenPos_Handle,double,&Face_CenPos_Handle::ycoord>,
member<Face_CenPos_Handle,double,&Face_CenPos_Handle::zcoord> > >
> > Face_CenPos_Handle_multiIndex;

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
  
  //============================================================================
  //
  //  I. READ MESH DATA AND FINITE ELEMENT ANALYSIS CONTROL PARAMETERS FROM FILE
  //
  //============================================================================
  
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
  cout<<"applied_strain ="<<applied_strain<<endl;
  
  
  /*****************************************************************************
   *
   * Transfer mesh data to MOAB database
   *
   ****************************************************************************/
  // Read mesh to MOAB
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
  
  for(_IT_CUBITMESHSETS_BY_SET_TYPE_FOR_LOOP_(mField,BLOCKSET,it)){
    
    if(it->get_name() == "MAT_ELASTIC_1") {
      Range TetsInBlock;
      rval = moab.get_entities_by_type(it->meshset, MBTET,TetsInBlock,true); CHKERR_PETSC(rval);
      Range block_rope_bit_level = intersect(LatestRefinedTets,TetsInBlock);
      rval = moab.add_entities(meshset_Elastic,block_rope_bit_level);CHKERR_PETSC(rval);
    }
  }
  ierr = mField.seed_finite_elements(meshset_Elastic); CHKERRQ(ierr);
  
  Range prims_on_problem_bit_level;
  ierr = mField.get_entities_by_type_and_ref_level(problem_bit_level,BitRefLevel().set(),MBPRISM,prims_on_problem_bit_level); CHKERRQ(ierr);
  //to create meshset from range
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
  ierr = mField.add_field("Lagrange_mul_disp",H1,field_rank); CHKERRQ(ierr);  //For lagrange multipliers to control the periodic motion
  ierr = mField.add_field("Lagrange_mul_disp_rigid_trans",NOFIELD,3); CHKERRQ(ierr);  //To control the rigid body motion (3 Traslations)
  
  /*****************************************************************************
   *
   * Add stochastic field
   * (total 14 field for 1st and 2nd order stochastic PSFEM)
   *
   ****************************************************************************/
  for(int ii=0; ii < nders; ii++ ) {
    ostringstream ss_field;
    ss_field << "DISP" << stochastic_fields[ii];
    cout<<"Here we are: \t"<<ii<<"\t"<<ss_field.str().c_str()<<endl;
    ierr = mField.add_field(ss_field.str().c_str(),H1,field_rank,MF_ZERO); CHKERRQ(ierr);
  }
  
  
  /*****************************************************************************
   *
   * Create finite element for the defined fields
   *
   ****************************************************************************/
  //FE
  ierr = mField.add_finite_element("ELASTIC"); CHKERRQ(ierr);
  ierr = mField.add_finite_element("TRAN_ISOTROPIC_ELASTIC"); CHKERRQ(ierr);
  ierr = mField.add_finite_element("Lagrange_elem"); CHKERRQ(ierr);
  ierr = mField.add_finite_element("Lagrange_elem_rigid_trans"); CHKERRQ(ierr);//For rigid body control
  
  
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
  for(int ii=0; ii < nvars; ii++ ) {
    ostringstream ss_field;
    ss_field << "DISP" << stochastic_fields[ii];
    ierr = mField.modify_finite_element_add_field_data("ELASTIC",ss_field.str().c_str()); CHKERRQ(ierr);
  }
  
  //FE Transverse Isotropic
  ierr = mField.modify_finite_element_add_field_row("TRAN_ISOTROPIC_ELASTIC","DISPLACEMENT"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_col("TRAN_ISOTROPIC_ELASTIC","DISPLACEMENT"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_data("TRAN_ISOTROPIC_ELASTIC","DISPLACEMENT"); CHKERRQ(ierr);
  
  
  //adding stochastic field to TRAN_ISOTROPIC_ELASTIC element
  for(int ii=0; ii < nvars; ii++ )
  {
    ostringstream ss_field;
    ss_field << "DISP" << stochastic_fields[ii];
    //    cout<<ss_field.str().c_str()<<endl;
    ierr = mField.modify_finite_element_add_field_data("TRAN_ISOTROPIC_ELASTIC",ss_field.str().c_str()); CHKERRQ(ierr);
  }
  
  ierr = mField.modify_finite_element_add_field_data("TRAN_ISOTROPIC_ELASTIC","POTENTIAL_FIELD"); CHKERRQ(ierr);
  
  //Define rows/cols and element data for C and CT (for lagrange multipliers)
  //==========================================================================
  //C row as Lagrange_mul_disp and col as DISPLACEMENT
  ierr = mField.modify_finite_element_add_field_row("Lagrange_elem","Lagrange_mul_disp"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_col("Lagrange_elem","DISPLACEMENT"); CHKERRQ(ierr);
  
  //CT col as Lagrange_mul_disp and row as DISPLACEMENT
  ierr = mField.modify_finite_element_add_field_col("Lagrange_elem","Lagrange_mul_disp"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_row("Lagrange_elem","DISPLACEMENT"); CHKERRQ(ierr);
  
  //As for stress we need both displacement and temprature (Lukasz)
  ierr = mField.modify_finite_element_add_field_data("Lagrange_elem","Lagrange_mul_disp"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_data("Lagrange_elem","DISPLACEMENT"); CHKERRQ(ierr);
  //==========================================================================
  
  
  //Define rows/cols and element data for C1 and C1T (for lagrange multipliers to contol the rigid body motion)
  //==========================================================================
  //C1 row as Lagrange_mul_disp and col as DISPLACEMENT
  ierr = mField.modify_finite_element_add_field_row("Lagrange_elem_rigid_trans","Lagrange_mul_disp_rigid_trans"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_col("Lagrange_elem_rigid_trans","DISPLACEMENT"); CHKERRQ(ierr);
  
  //C1T col as Lagrange_mul_disp and row as DISPLACEMENT
  ierr = mField.modify_finite_element_add_field_col("Lagrange_elem_rigid_trans","Lagrange_mul_disp_rigid_trans"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_row("Lagrange_elem_rigid_trans","DISPLACEMENT"); CHKERRQ(ierr);
  
  //As for stress we need both displacement and temprature (Lukasz)
  ierr = mField.modify_finite_element_add_field_data("Lagrange_elem_rigid_trans","Lagrange_mul_disp_rigid_trans"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_data("Lagrange_elem_rigid_trans","DISPLACEMENT"); CHKERRQ(ierr);
  //==========================================================================
  
  
  //define problems
  ierr = mField.add_problem("STOCHASIC_PROBLEM"); CHKERRQ(ierr);
  
  //set finite elements for problem
  ierr = mField.modify_problem_add_finite_element("STOCHASIC_PROBLEM","ELASTIC"); CHKERRQ(ierr);
  ierr = mField.modify_problem_add_finite_element("STOCHASIC_PROBLEM","TRAN_ISOTROPIC_ELASTIC"); CHKERRQ(ierr);
  ierr = mField.modify_problem_add_finite_element("STOCHASIC_PROBLEM","Lagrange_elem"); CHKERRQ(ierr);
  ierr = mField.modify_problem_add_finite_element("STOCHASIC_PROBLEM","Lagrange_elem_rigid_trans"); CHKERRQ(ierr);
  
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
  
  
  for(int ii=0; ii < nders; ii++ ) {
    ostringstream ss_field;
    ss_field << "DISP" << stochastic_fields[ii];
    cout<<"Here we are: \t"<<ii<<"\t"<<ss_field.str().c_str()<<endl;
    ierr = mField.add_ents_to_field_by_TETs(0,ss_field.str().c_str()); CHKERRQ(ierr);
  }
  
  /*****************************************************************************
   *
   * Add finite elements entities
   *
   ****************************************************************************/
  ierr = mField.add_ents_to_finite_element_by_TETs(meshset_Elastic,"ELASTIC",true); CHKERRQ(ierr);
  ierr = mField.add_ents_to_finite_element_by_TETs(meshset_Trans_ISO,"TRAN_ISOTROPIC_ELASTIC",true); CHKERRQ(ierr);
  
  
  //Add finite element to lagrange element for rigid body translation
  Range Tris_NewWholeMesh, Tri_OldNewSurf, SurfacesFaces;
  ierr = mField.get_entities_by_type_and_ref_level(bit_levels.back(),BitRefLevel().set(),MBTRI,Tris_NewWholeMesh); CHKERRQ(ierr);
  ierr = mField.get_cubit_msId_entities_by_dimension(103,SIDESET,2,Tri_OldNewSurf,true); CHKERRQ(ierr);
  SurfacesFaces = intersect(Tris_NewWholeMesh,Tri_OldNewSurf);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"number of SIDESET 103 = %d\n",SurfacesFaces.size()); CHKERRQ(ierr);
  
  
  //to create meshset from range
  EntityHandle BoundFacesMeshset;
  rval = moab.create_meshset(MESHSET_SET,BoundFacesMeshset); CHKERR_PETSC(rval);
  rval = moab.add_entities(BoundFacesMeshset,SurfacesFaces); CHKERR_PETSC(rval);
  ierr = mField.seed_ref_level_MESHSET(BoundFacesMeshset,BitRefLevel().set()); CHKERRQ(ierr);
  
  ierr = mField.add_ents_to_finite_element_by_TRIs(SurfacesFaces,"Lagrange_elem_rigid_trans"); CHKERRQ(ierr);
  
  
  //=======================================================================================================
  //Add Periodic Prisims Between Triangles on -ve and +ve faces to implement periodic bounary conditions
  //=======================================================================================================
  
  //Populating the Multi-index container with -ve triangles
  Range Tri_OldNewSurfNeg, SurTrisNeg;
  ierr = mField.get_cubit_msId_entities_by_dimension(101,SIDESET,2,Tri_OldNewSurfNeg,true); CHKERRQ(ierr);
  SurTrisNeg = intersect(Tris_NewWholeMesh,Tri_OldNewSurfNeg);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"number of SIDESET 101 = %d\n",SurTrisNeg.size()); CHKERRQ(ierr);
  
  Face_CenPos_Handle_multiIndex Face_CenPos_Handle_varNeg, Face_CenPos_Handle_varPos;
  double TriCen[3], coords_Tri[9];
  
  double roundfact=1000.0;
  for(Range::iterator it = SurTrisNeg.begin(); it!=SurTrisNeg.end();  it++) {
    //        cout<<"count1 ="<<count1<<endl;
    const EntityHandle* conn_face;  int num_nodes_Tri;
    
    //get nodes attached to the face
    rval = moab.get_connectivity(*it,conn_face,num_nodes_Tri,true); CHKERR_PETSC(rval);
    //get nodal coordinates
    rval = moab.get_coords(conn_face,num_nodes_Tri,coords_Tri); CHKERR_PETSC(rval);
    
    //Find triangle centriod
    TriCen[0]= (coords_Tri[0]+coords_Tri[3]+coords_Tri[6])/3.0;
    TriCen[1]= (coords_Tri[1]+coords_Tri[4]+coords_Tri[7])/3.0;
    TriCen[2]= (coords_Tri[2]+coords_Tri[5]+coords_Tri[8])/3.0;
    
    //round values to 3 disimal places
    if(TriCen[0]>=0) TriCen[0]=double(int(TriCen[0]*roundfact+0.5))/roundfact;  else TriCen[0]=double(int(TriCen[0]*roundfact-0.5))/roundfact; //-ve and +ve value
    if(TriCen[1]>=0) TriCen[1]=double(int(TriCen[1]*roundfact+0.5))/roundfact;  else TriCen[1]=double(int(TriCen[1]*roundfact-0.5))/roundfact;
    if(TriCen[2]>=0) TriCen[2]=double(int(TriCen[2]*roundfact+0.5))/roundfact;  else TriCen[2]=double(int(TriCen[2]*roundfact-0.5))/roundfact;
    //        cout<<"   TriCen[0]= "<<TriCen[0] << "   TriCen[1]= "<< TriCen[1] << "   TriCen[2]= "<< TriCen[2] <<endl;
    //fill the multi-index container with centriod coordinates and triangle handles
    Face_CenPos_Handle_varNeg.insert(Face_CenPos_Handle(TriCen[0], TriCen[1], TriCen[2], *it));
  }
  
  //    typedef Face_CenPos_Handle_multiIndex::index<Tri_Hand_tag>::type::iterator Tri_Hand_iterator;
  //    Tri_Hand_iterator Tri_Neg;
  //    count1=1;
  //    for(Range::iterator it = SurTrisNeg.begin(); it!=SurTrisNeg.end();  it++) {
  //        Tri_Neg=Face_CenPos_Handle_varNeg.get<Tri_Hand_tag>().find(*it);
  //        cout<<"it= "<<*it <<"     Tri_Neg->xcoord= "<<Tri_Neg->xcoord << "   Tri_Neg->ycoord "<< Tri_Neg->ycoord << "   Tri_Neg->zcoord= "<< Tri_Neg->zcoord <<endl;
  //    }
  
  
  //Populating the Multi-index container with +ve triangles
  Range Tri_OldNewSurfPos, SurTrisPos;
  ierr = mField.get_cubit_msId_entities_by_dimension(102,SIDESET,2,Tri_OldNewSurfPos,true); CHKERRQ(ierr);
  SurTrisPos = intersect(Tris_NewWholeMesh,Tri_OldNewSurfPos);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"number of SIDESET 102 = %d\n",SurTrisPos.size()); CHKERRQ(ierr);
  
  for(Range::iterator it = SurTrisPos.begin(); it!=SurTrisPos.end();  it++) {
    const EntityHandle* conn_face;  int num_nodes_Tri;
    
    //get nodes attached to the face
    rval = moab.get_connectivity(*it,conn_face,num_nodes_Tri,true); CHKERR_PETSC(rval);
    //get nodal coordinates
    rval = moab.get_coords(conn_face,num_nodes_Tri,coords_Tri); CHKERR_PETSC(rval);
    
    //Find triangle centriod
    TriCen[0]= (coords_Tri[0]+coords_Tri[3]+coords_Tri[6])/3.0;
    TriCen[1]= (coords_Tri[1]+coords_Tri[4]+coords_Tri[7])/3.0;
    TriCen[2]= (coords_Tri[2]+coords_Tri[5]+coords_Tri[8])/3.0;
    
    //round values to 3 disimal places
    if(TriCen[0]>=0) TriCen[0]=double(int(TriCen[0]*roundfact+0.5))/roundfact;  else TriCen[0]=double(int(TriCen[0]*roundfact-0.5))/roundfact;
    if(TriCen[1]>=0) TriCen[1]=double(int(TriCen[1]*roundfact+0.5))/roundfact;  else TriCen[1]=double(int(TriCen[1]*roundfact-0.5))/roundfact;
    if(TriCen[2]>=0) TriCen[2]=double(int(TriCen[2]*roundfact+0.5))/roundfact;  else TriCen[2]=double(int(TriCen[2]*roundfact-0.5))/roundfact;
    //        cout<<"TriCen[0]= "<<TriCen[0] << "   TriCen[1]= "<< TriCen[1] << "   TriCen[2]= "<< TriCen[2] <<endl;
    
    //fill the multi-index container with centriod coordinates and triangle handles
    Face_CenPos_Handle_varPos.insert(Face_CenPos_Handle(TriCen[0], TriCen[1], TriCen[2], *it));
  }
  
  
  //    typedef Face_CenPos_Handle_multiIndex::index<Tri_Hand_tag>::type::iterator Tri_Hand_iterator;
  //    Tri_Hand_iterator Tri_Pos;
  //    count1=1;
  //    for(Range::iterator it = SurTrisPos.begin(); it!=SurTrisPos.end();  it++) {
  //        Tri_Pos=Face_CenPos_Handle_varPos.get<Tri_Hand_tag>().find(*it);
  //        cout<<"it= "<<*it <<"     Tri_Pos->xcoord= "<<Tri_Pos->xcoord << "   Tri_Pos->ycoord "<< Tri_Pos->ycoord << "   Tri_Pos->zcoord= "<< Tri_Pos->zcoord <<endl;
  //    }
  
  
  
  //    //Find minimum and maximum X, Y and Z coordinates of the RVE (using multi-index container)
  double XcoordMin, YcoordMin, ZcoordMin, XcoordMax, YcoordMax, ZcoordMax;
  typedef Face_CenPos_Handle_multiIndex::index<xcoord_tag>::type::iterator Tri_Xcoord_iterator;
  typedef Face_CenPos_Handle_multiIndex::index<ycoord_tag>::type::iterator Tri_Ycoord_iterator;
  typedef Face_CenPos_Handle_multiIndex::index<zcoord_tag>::type::iterator Tri_Zcoord_iterator;
  Tri_Xcoord_iterator XcoordMin_it, XcoordMax_it;
  Tri_Ycoord_iterator YcoordMin_it, YcoordMax_it;
  Tri_Zcoord_iterator ZcoordMin_it, ZcoordMax_it;
  
  //XcoordMax_it-- because .end() will point iterator after the data range but .begin() will point the iteratore to the first value of range
  XcoordMin_it=Face_CenPos_Handle_varNeg.get<xcoord_tag>().begin();                  XcoordMin=XcoordMin_it->xcoord;
  XcoordMax_it=Face_CenPos_Handle_varPos.get<xcoord_tag>().end();    XcoordMax_it--; XcoordMax=XcoordMax_it->xcoord;
  YcoordMin_it=Face_CenPos_Handle_varNeg.get<ycoord_tag>().begin();                  YcoordMin=YcoordMin_it->ycoord;
  YcoordMax_it=Face_CenPos_Handle_varPos.get<ycoord_tag>().end();    YcoordMax_it--; YcoordMax=YcoordMax_it->ycoord;
  ZcoordMin_it=Face_CenPos_Handle_varNeg.get<zcoord_tag>().begin();                  ZcoordMin=ZcoordMin_it->zcoord;
  ZcoordMax_it=Face_CenPos_Handle_varPos.get<zcoord_tag>().end();    ZcoordMax_it--; ZcoordMax=ZcoordMax_it->zcoord;
  
  cout<<"XcoordMin "<<XcoordMin << "      XcoordMax "<<XcoordMax <<endl;
  cout<<"YcoordMin "<<YcoordMin << "      YcoordMax "<<YcoordMax <<endl;
  cout<<"ZcoordMin "<<ZcoordMin << "      ZcoordMax "<<ZcoordMax <<endl;
  
  //    XcoordMin =-1.5;      XcoordMax =1.5;
  //    YcoordMin =-0.15;      YcoordMax =0.15;
  //    ZcoordMin =-0.195;      ZcoordMax =0.585;
  
  //Creating Prisims between triangles on -ve and +ve faces
  typedef Face_CenPos_Handle_multiIndex::index<Tri_Hand_tag>::type::iterator Tri_Hand_iterator;
  Tri_Hand_iterator Tri_Neg;
  typedef Face_CenPos_Handle_multiIndex::index<Composite_xyzcoord>::type::iterator xyzcoord_iterator;
  xyzcoord_iterator Tri_Pos;
  Range PrismRange;
  double XPos, YPos, ZPos;
  //int count=0;
  
  //loop over -ve triangles to create prisims elemenet between +ve and -ve triangles
  //    count1=1;
  for(Range::iterator it = SurTrisNeg.begin(); it!=SurTrisNeg.end();  it++) {
    //        cout<<"count1 ="<<count1<<endl;  count1++;
    Tri_Neg=Face_CenPos_Handle_varNeg.get<Tri_Hand_tag>().find(*it);
    //        cout<<"Tri_Neg->xcoord= "<<Tri_Neg->xcoord << "   Tri_Neg->ycoord "<< Tri_Neg->ycoord << "   Tri_Neg->zcoord= "<< Tri_Neg->zcoord <<endl;
    //corresponding +ve triangle
    if(Tri_Neg->xcoord==XcoordMin){XPos=XcoordMax;         YPos=Tri_Neg->ycoord;  ZPos=Tri_Neg->zcoord;};
    if(Tri_Neg->ycoord==YcoordMin){XPos=Tri_Neg->xcoord;   YPos=YcoordMax;        ZPos=Tri_Neg->zcoord;};
    if(Tri_Neg->zcoord==ZcoordMin){XPos=Tri_Neg->xcoord;   YPos=Tri_Neg->ycoord;  ZPos=ZcoordMax;      };
    
    //        cout<<"Tri_Neg->xcoord= "<<Tri_Neg->xcoord << "   Tri_Neg->ycoord "<< Tri_Neg->ycoord << "   Tri_Neg->zcoord= "<< Tri_Neg->zcoord <<endl;
    //        cout<<"XPos= "<<XPos << "   YPos "<< YPos << "   ZPos= "<< ZPos <<endl;
    Tri_Pos=Face_CenPos_Handle_varPos.get<Composite_xyzcoord>().find(boost::make_tuple(XPos, YPos, ZPos));
    //        cout<<"Tri_Pos->xcoord= "<<Tri_Pos->xcoord << "   Tri_Pos->ycoord "<< Tri_Pos->ycoord << "   Tri_Pos->zcoord= "<< Tri_Pos->zcoord <<endl;
    
    //+ve and -ve nodes and their coords (+ve and -ve tiangles nodes can have matching problems, which can produce twisted prism)
    EntityHandle PrismNodes[6];
    vector<EntityHandle> TriNodesNeg, TriNodesPos;
    double CoordNodeNeg[9], CoordNodePos[9];
    rval = moab.get_connectivity(&(Tri_Neg->Tri_Hand),1,TriNodesNeg,true); CHKERR_PETSC(rval);
    rval = moab.get_connectivity(&(Tri_Pos->Tri_Hand),1,TriNodesPos,true); CHKERR_PETSC(rval);
    rval = moab.get_coords(&TriNodesNeg[0],3,CoordNodeNeg);  CHKERR_THROW(rval);
    rval = moab.get_coords(&TriNodesPos[0],3,CoordNodePos);  CHKERR_THROW(rval);
    for(int ii=0; ii<3; ii++){
      PrismNodes[ii]=TriNodesNeg[ii];
    }
    //        for(int ii=0; ii<3; ii++){
    //            cout<<"xcoord= "<<CoordNodeNeg[3*ii] << "   ycoord= "<< CoordNodeNeg[3*ii+1] << "   zcoord= "<< CoordNodeNeg[3*ii+2] <<endl;
    //        }
    //        for(int ii=0; ii<3; ii++){
    //            cout<<"xcoord= "<<CoordNodePos[3*ii] << "   ycoord= "<< CoordNodePos[3*ii+1] << "   zcoord= "<< CoordNodePos[3*ii+2] <<endl;
    //        }
    
    //Match exact nodes to each other to avoide the problem of twisted prisms
    double XNodeNeg, YNodeNeg, ZNodeNeg, XNodePos, YNodePos, ZNodePos;
    for(int ii=0; ii<3; ii++){
      if(Tri_Neg->xcoord==XcoordMin){XNodeNeg=XcoordMax;          YNodeNeg=CoordNodeNeg[3*ii+1];   ZNodeNeg=CoordNodeNeg[3*ii+2];};
      if(Tri_Neg->ycoord==YcoordMin){XNodeNeg=CoordNodeNeg[3*ii]; YNodeNeg=YcoordMax;              ZNodeNeg=CoordNodeNeg[3*ii+2];};
      if(Tri_Neg->zcoord==ZcoordMin){XNodeNeg=CoordNodeNeg[3*ii]; YNodeNeg=CoordNodeNeg[3*ii+1];   ZNodeNeg=ZcoordMax;};
      for(int jj=0; jj<3; jj++){
        XNodePos=CoordNodePos[3*jj]; YNodePos=CoordNodePos[3*jj+1]; ZNodePos=CoordNodePos[3*jj+2];
        
        if(XNodeNeg==XNodePos  &&  YNodeNeg==YNodePos  &&  ZNodeNeg==ZNodePos){
          PrismNodes[3+ii]=TriNodesPos[jj];
          break;
        }
      }
    }
    //prism nodes and their coordinates
    double CoordNodesPrisms[18];
    rval = moab.get_coords(PrismNodes,6,CoordNodesPrisms);  CHKERR_THROW(rval);
    //        for(int ii=0; ii<6; ii++){
    //            cout<<"xcoord= "<<CoordNodesPrisms[3*ii] << "   ycoord= "<< CoordNodesPrisms[3*ii+1] << "   zcoord= "<< CoordNodesPrisms[3*ii+2] <<endl;
    //        }
    //        cout<<endl<<endl;
    //insertion of individula prism element and its addition to range PrismRange
    EntityHandle PeriodicPrism;
    rval = moab.create_element(MBPRISM,PrismNodes,6,PeriodicPrism); CHKERR_PETSC(rval);
    PrismRange.insert(PeriodicPrism);
    
    //        //to see individual prisms
    //        Range Prism1;
    //        Prism1.insert(PeriodicPrism);
    //        EntityHandle out_meshset1;
    //        rval = moab.create_meshset(MESHSET_SET,out_meshset1); CHKERR_PETSC(rval);
    //        rval = moab.add_entities(out_meshset1,Prism1); CHKERR_PETSC(rval);
    //        ostringstream sss;
    //        sss << "Prism" << count << ".vtk"; count++;
    //        rval = moab.write_file(sss.str().c_str(),"VTK","",&out_meshset1,1); CHKERR_PETSC(rval);
  }
  
  
  //cout<<"PrismRange "<<PrismRange<<endl;
  //Saving prisms in interface.vtk
  EntityHandle out_meshset1;
  rval = moab.create_meshset(MESHSET_SET,out_meshset1); CHKERR_PETSC(rval);
  rval = moab.add_entities(out_meshset1,PrismRange); CHKERR_PETSC(rval);
  rval = moab.write_file("Prisms.vtk","VTK","",&out_meshset1,1); CHKERR_PETSC(rval);
  cout << "Prisms.vtk output" <<endl;
  
  //Adding Prisims to Element Lagrange_elem (to loop over these prisims)
  EntityHandle PrismRangeMeshset;
  rval = moab.create_meshset(MESHSET_SET,PrismRangeMeshset); CHKERR_PETSC(rval);
  rval = moab.add_entities(PrismRangeMeshset,PrismRange); CHKERR_PETSC(rval);
  //    cout << PrismRange <<endl;
  ierr = mField.seed_ref_level_3D(PrismRangeMeshset,problem_bit_level); CHKERRQ(ierr);
  //    mField.seed_finite_elements(PrismRange);
  ierr = mField.get_entities_by_ref_level(bit_levels.back(),BitRefLevel().set(),out_meshset); CHKERRQ(ierr);
  ierr = mField.add_ents_to_finite_element_by_PRISMs(PrismRange, "Lagrange_elem"); CHKERRQ(ierr);
  
  
  //Adding only -ve surfaces to the field Lagrange_mul_disp (in periodic boundary conditions size of C (3M/2 x 3N))
  //to create meshset from range  SurTrisNeg
  EntityHandle SurTrisNegMeshset;
  rval = moab.create_meshset(MESHSET_SET,SurTrisNegMeshset); CHKERR_PETSC(rval);
  rval = moab.add_entities(SurTrisNegMeshset,SurTrisNeg); CHKERR_PETSC(rval);
  ierr = mField.add_ents_to_field_by_TRIs(SurTrisNegMeshset,"Lagrange_mul_disp",2); CHKERRQ(ierr);
  
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
  Vec dF,ddF,D,dD,ddD;
  
  vector<Vec> F(6);
  ierr = mField.VecCreateGhost("STOCHASIC_PROBLEM",ROW,&F[0]); CHKERRQ(ierr);
  for(int ii = 1;ii<6;ii++) {
    ierr = VecDuplicate(F[0],&F[ii]); CHKERRQ(ierr);
  }
  
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
  //    double YoungModulusP;
  //    double PoissonRatioP;
  //    double YoungModulusZ;
  //    double PoissonRatioPZ;
  //    double ShearModulusZP;
  double YoungModulus = 3500;
  double PoissonRatio = 0.3;
  double alpha;
  
  /*****************************************************************************
   *
   *  2. Get the volume of RVE
   *
   ****************************************************************************/
  double RVE_volume;    RVE_volume=0.0;  //RVE volume for full RVE We need this for stress calculation
  Vec RVE_volume_Vec;
  ierr = VecCreateMPI(PETSC_COMM_WORLD, 1, pcomm->size(), &RVE_volume_Vec);  CHKERRQ(ierr);
  ierr = VecZeroEntries(RVE_volume_Vec); CHKERRQ(ierr);
  
  RVEVolume MyRVEVol(mField,Aij,D,F[0],LAMBDA(YoungModulus,PoissonRatio),MU(YoungModulus,PoissonRatio), RVE_volume_Vec);
  RVEVolumeTrans MyRVEVolTrans(mField,Aij,D,F[0], RVE_volume_Vec);
  
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
  
  for(_IT_CUBITMESHSETS_BY_SET_TYPE_FOR_LOOP_(mField,BLOCKSET,it)) {
    cout << endl << *it << endl;
    
    //Get block name
    string name = it->get_name();
    
    //        if (name.compare(0,20,"MAT_ELASTIC_TRANSISO") == 0)
    //        {
    //            Mat_Elastic_TransIso mydata;
    //            ierr = it->get_attribute_data_structure(mydata); CHKERRQ(ierr);
    //            cout << mydata;
    //            YoungModulusP=mydata.data.Youngp;
    //            YoungModulusZ=mydata.data.Youngz;
    //            PoissonRatioP=mydata.data.Poissonp;
    //            PoissonRatioPZ=mydata.data.Poissonpz;
    //            if (mydata.data.Shearzp!=0) {
    //                ShearModulusZP=mydata.data.Shearzp;
    //            }else{
    //                ShearModulusZP=YoungModulusZ/(2*(1+PoissonRatioPZ));}
    //        }
    if (name.compare(0,11,"MAT_ELASTIC") == 0)
    {
      Mat_Elastic mydata;
      ierr = it->get_attribute_data_structure(mydata); CHKERRQ(ierr);
      cout << mydata;
      YoungModulus=mydata.data.Young;
      PoissonRatio=mydata.data.Poisson;
    }
    else if (name.compare(0,10,"MAT_INTERF") == 0)
    {
      Mat_Interf mydata;
      ierr = it->get_attribute_data_structure(mydata); CHKERRQ(ierr);
      cout << mydata;
      alpha = mydata.data.alpha;
    }
  }
  
  //    alpha = 500;
  cout<<"alpha   = "<<alpha<<endl;
  
  applied_strain.resize(1.5*field_rank+1.5); applied_strain.clear();
  
  MyElasticFEMethod MyFE(mField,Aij,D,F[0],LAMBDA(YoungModulus,PoissonRatio),MU(YoungModulus,PoissonRatio));
  TranIsotropicFibreDirRotElasticFEMethod MyTIsotFE(mField,Aij,D,F[0]);
  
  // ElasticFE_RVELagrange_Periodic MyFE_RVELagrangePeriodic(mField,Aij,D,F,applied_strain,"DISPLACEMENT","Lagrange_mul_disp",field_rank);
  ElasticFE_RVELagrange_Periodic_Multi_Rhs MyFE_RVELagrangePeriodic(mField,Aij,D,F[0],F[1],F[2],F[3],F[4],F[5],applied_strain,"DISPLACEMENT","Lagrange_mul_disp",field_rank);
  ElasticFE_RVELagrange_RigidBodyTranslation MyFE_RVELagrangeRigidBodyTrans(mField,Aij,D,F[0],applied_strain,"DISPLACEMENT","Lagrange_mul_disp",field_rank,"Lagrange_mul_disp_rigid_trans");
  
  for(int ii = 0; ii<6; ii++) {
    ierr = VecZeroEntries(F[ii]); CHKERRQ(ierr);
    ierr = VecGhostUpdateBegin(F[ii],INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
    ierr = VecGhostUpdateEnd(F[ii],INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  }
  ierr = MatZeroEntries(Aij); CHKERRQ(ierr);
  ierr = VecZeroEntries(D); CHKERRQ(ierr);
  ierr = mField.set_global_ghost_vector("STOCHASIC_PROBLEM",ROW,D,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
  
  ierr = mField.loop_finite_elements("STOCHASIC_PROBLEM","ELASTIC",MyFE);  CHKERRQ(ierr);
  PetscSynchronizedFlush(PETSC_COMM_WORLD,PETSC_STDOUT);
  ierr = mField.loop_finite_elements("STOCHASIC_PROBLEM","TRAN_ISOTROPIC_ELASTIC",MyTIsotFE);  CHKERRQ(ierr);
  PetscSynchronizedFlush(PETSC_COMM_WORLD,PETSC_STDOUT);
  ierr = mField.loop_finite_elements("STOCHASIC_PROBLEM","Lagrange_elem",MyFE_RVELagrangePeriodic);  CHKERRQ(ierr);
  PetscSynchronizedFlush(PETSC_COMM_WORLD,PETSC_STDOUT);
  ierr = mField.loop_finite_elements("STOCHASIC_PROBLEM","Lagrange_elem_rigid_trans",MyFE_RVELagrangeRigidBodyTrans);  CHKERRQ(ierr);
  PetscSynchronizedFlush(PETSC_COMM_WORLD,PETSC_STDOUT);
  
  for(int ii = 0; ii<6; ii++) {
    ierr = VecGhostUpdateBegin(F[ii],ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
    ierr = VecGhostUpdateEnd(F[ii],ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
    ierr = VecAssemblyBegin(F[ii]); CHKERRQ(ierr);
    ierr = VecAssemblyEnd(F[ii]); CHKERRQ(ierr);
  }
  
  ierr = MatAssemblyBegin(Aij,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
  ierr = MatAssemblyEnd(Aij,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
  
  //    //Matrix View
  //    MatView(Aij,PETSC_VIEWER_DRAW_WORLD);//PETSC_VIEWER_STDOUT_WORLD);
  //    std::string wait;
  //    std::cin >> wait;
  
  
  //----------------------------------------------------------------------------
  // 3.1 Solving the equation to get nodal displacement
  //----------------------------------------------------------------------------
  KSP solver;
  ierr = KSPCreate(PETSC_COMM_WORLD,&solver); CHKERRQ(ierr);
  ierr = KSPSetOperators(solver,Aij,Aij); CHKERRQ(ierr);
  ierr = KSPSetFromOptions(solver); CHKERRQ(ierr);
  ierr = KSPSetUp(solver); CHKERRQ(ierr);
  
  //create a vector for 6 components of homogenized stress
  Vec Stress_Homo, Stress_Homo_r;
  ierr = VecCreateMPI(PETSC_COMM_WORLD, 6, 6*pcomm->size(), &Stress_Homo);  CHKERRQ(ierr);
  ierr = VecCreateMPI(PETSC_COMM_WORLD, 6, 6*pcomm->size(), &Stress_Homo_r);  CHKERRQ(ierr);
  
  ublas::matrix<double> Dmat(6,6); Dmat.clear();
  ublas::vector<ublas::matrix<double> > Dmat_r(nvars);
  ublas::vector<ublas::matrix<double> > Dmat_rs(nvars);
  for(int irv=0; irv<nvars; irv++) {
    Dmat_r(irv).resize(6,6);   Dmat_r(irv).clear();
    Dmat_rs(irv).resize(6,6); Dmat_rs(irv).clear();
  }
  
  for(int ics = 0; ics<6; ics++) {
    ierr = VecZeroEntries(D); CHKERRQ(ierr);
    ierr = KSPSolve(solver,F[ics],D); CHKERRQ(ierr);
    ierr = VecGhostUpdateBegin(D,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
    ierr = VecGhostUpdateEnd(D,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
    //Save data on mesh
    ierr = mField.set_global_VecCreateGhost("STOCHASIC_PROBLEM",ROW,D,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
    //    ierr = VecView(D,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);
    
    //----------------------------------------------------------------------------
    // 3.2 Calculating zeroth-order homogenized stress using volume averaging theorem
    //----------------------------------------------------------------------------
    //create a vector for 6 components of homogenized stress
    ierr = VecCreateMPI(PETSC_COMM_WORLD, 6, 6*pcomm->size(), &Stress_Homo);  CHKERRQ(ierr);
    ierr = VecZeroEntries(Stress_Homo); CHKERRQ(ierr);
    
    //    ierr = VecView(D,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);
    ElasticFE_RVELagrange_Homogenized_Stress_Periodic MyFE_RVEHomoStressPeriodic(mField,Aij,D,F[ics],&RVE_volume,applied_strain,Stress_Homo,"DISPLACEMENT","Lagrange_mul_disp",field_rank);
    ierr = mField.loop_finite_elements("STOCHASIC_PROBLEM","Lagrange_elem",MyFE_RVEHomoStressPeriodic);  CHKERRQ(ierr);
    
    //    if(pcomm->rank()) cout<< " Stress_Homo =  "<<endl;
    //    ierr = VecView(Stress_Homo,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);
    if(pcomm->rank()==0){
      PetscScalar    *avec;
      VecGetArray(Stress_Homo, &avec);
      
      cout<< "\nStress_Homo =\n";
      for(int kk=0; kk<6; kk++){
        Dmat(kk,ics) = *avec;
        //cout.precision(15);
        //cout <<*avec<<endl;
        avec++;
      }
      VecRestoreArray(Stress_Homo,&avec);
    }
    //cout<< "\n\n";
    
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
      ierr = VecZeroEntries(dD); CHKERRQ(ierr);
      ierr = VecGhostUpdateBegin(dF,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
      ierr = VecGhostUpdateEnd(dF,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
      //
      ierr = VecZeroEntries(ddF); CHKERRQ(ierr);
      ierr = VecZeroEntries(ddD); CHKERRQ(ierr);
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
      else if (ii == 14){ // 2nd order derivative due to shear modulus in z-direction of fibre
        Trans_Iso_Rhs_rs_PSFEM my_fe_k_rs_NUpNUpz(mField,Aij,D,ddF,"DISPLACEMENT","DISP_r_NUpz","PoissonP","PoissonZ", "transversely_isotropic", "reinforcement");
        ierr = mField.loop_finite_elements("STOCHASIC_PROBLEM","TRAN_ISOTROPIC_ELASTIC",my_fe_k_rs_NUpNUpz);  CHKERRQ(ierr);
      }
      else if (ii == 15){ // 2nd order derivative due to shear modulus in z-direction of fibre
        Trans_Iso_Rhs_rs_PSFEM my_fe_k_rs_NUpEp(mField,Aij,D,ddF,"DISPLACEMENT","DISP_r_Ep","PoissonP","YoungP", "transversely_isotropic", "reinforcement");
        ierr = mField.loop_finite_elements("STOCHASIC_PROBLEM","TRAN_ISOTROPIC_ELASTIC",my_fe_k_rs_NUpEp);  CHKERRQ(ierr);
      }
      else if (ii == 16){ // 2nd order derivative due to shear modulus in z-direction of fibre
        Trans_Iso_Rhs_rs_PSFEM my_fe_k_rs_NUpEz(mField,Aij,D,ddF,"DISPLACEMENT","DISP_r_Ez","PoissonP","YoungZ", "transversely_isotropic", "reinforcement");
        ierr = mField.loop_finite_elements("STOCHASIC_PROBLEM","TRAN_ISOTROPIC_ELASTIC",my_fe_k_rs_NUpEz);  CHKERRQ(ierr);
      }
      else if (ii == 17){ // 2nd order derivative due to shear modulus in z-direction of fibre
        Trans_Iso_Rhs_rs_PSFEM my_fe_k_rs_NUpGzp(mField,Aij,D,ddF,"DISPLACEMENT","DISP_r_Gzp","PoissonP","ShearZP", "transversely_isotropic", "reinforcement");
        ierr = mField.loop_finite_elements("STOCHASIC_PROBLEM","TRAN_ISOTROPIC_ELASTIC",my_fe_k_rs_NUpGzp);  CHKERRQ(ierr);
      }
      else if (ii == 18){ // 2nd order derivative due to shear modulus in z-direction of fibre
        Trans_Iso_Rhs_rs_PSFEM my_fe_k_rs_NUpEm(mField,Aij,D,ddF,"DISPLACEMENT","DISP_r_Em","PoissonP","YoungM", "transversely_isotropic", "reinforcement");
        ierr = mField.loop_finite_elements("STOCHASIC_PROBLEM","TRAN_ISOTROPIC_ELASTIC",my_fe_k_rs_NUpEm);  CHKERRQ(ierr);
      }
      else if (ii == 19){ // 2nd order derivative due to shear modulus in z-direction of fibre
        Trans_Iso_Rhs_rs_PSFEM my_fe_k_rs_NUpPm(mField,Aij,D,ddF,"DISPLACEMENT","DISP_r_Pm","PoissonP","PoissonM", "transversely_isotropic", "reinforcement");
        ierr = mField.loop_finite_elements("STOCHASIC_PROBLEM","TRAN_ISOTROPIC_ELASTIC",my_fe_k_rs_NUpPm);  CHKERRQ(ierr);
      }
      else if (ii == 20){ // 2nd order derivative due to shear modulus in z-direction of fibre
        Trans_Iso_Rhs_rs_PSFEM my_fe_k_rs_NUpzEp(mField,Aij,D,ddF,"DISPLACEMENT","DISP_r_Ep","PoissonZ","YoungP", "transversely_isotropic", "reinforcement");
        ierr = mField.loop_finite_elements("STOCHASIC_PROBLEM","TRAN_ISOTROPIC_ELASTIC",my_fe_k_rs_NUpzEp);  CHKERRQ(ierr);
      }
      else if (ii == 21){ // 2nd order derivative due to shear modulus in z-direction of fibre
        Trans_Iso_Rhs_rs_PSFEM my_fe_k_rs_NUpzEz(mField,Aij,D,ddF,"DISPLACEMENT","DISP_r_Ez","PoissonZ","YoungZ", "transversely_isotropic", "reinforcement");
        ierr = mField.loop_finite_elements("STOCHASIC_PROBLEM","TRAN_ISOTROPIC_ELASTIC",my_fe_k_rs_NUpzEz);  CHKERRQ(ierr);
      }
      else if (ii == 22){ // 2nd order derivative due to shear modulus in z-direction of fibre
        Trans_Iso_Rhs_rs_PSFEM my_fe_k_rs_NUpzGzp(mField,Aij,D,ddF,"DISPLACEMENT","DISP_r_Gzp","PoissonZ","ShearZP", "transversely_isotropic", "reinforcement");
        ierr = mField.loop_finite_elements("STOCHASIC_PROBLEM","TRAN_ISOTROPIC_ELASTIC",my_fe_k_rs_NUpzGzp);  CHKERRQ(ierr);
      }
      else if (ii == 23){ // 2nd order derivative due to shear modulus in z-direction of fibre
        Trans_Iso_Rhs_rs_PSFEM my_fe_k_rs_NUpzEm(mField,Aij,D,ddF,"DISPLACEMENT","DISP_r_Em","PoissonZ","YoungM", "transversely_isotropic", "reinforcement");
        ierr = mField.loop_finite_elements("STOCHASIC_PROBLEM","TRAN_ISOTROPIC_ELASTIC",my_fe_k_rs_NUpzEm);  CHKERRQ(ierr);
      }
      else if (ii == 24){ // 2nd order derivative due to shear modulus in z-direction of fibre
        Trans_Iso_Rhs_rs_PSFEM my_fe_k_rs_NUpzPm(mField,Aij,D,ddF,"DISPLACEMENT","DISP_r_Pm","PoissonZ","PoissonM", "transversely_isotropic", "reinforcement");
        ierr = mField.loop_finite_elements("STOCHASIC_PROBLEM","TRAN_ISOTROPIC_ELASTIC",my_fe_k_rs_NUpzPm);  CHKERRQ(ierr);
      }
      else if (ii == 25){ // 2nd order derivative due to shear modulus in z-direction of fibre
        Trans_Iso_Rhs_rs_PSFEM my_fe_k_rs_EpEz(mField,Aij,D,ddF,"DISPLACEMENT","DISP_r_Ez","YoungP","YoungZ", "transversely_isotropic", "reinforcement");
        ierr = mField.loop_finite_elements("STOCHASIC_PROBLEM","TRAN_ISOTROPIC_ELASTIC",my_fe_k_rs_EpEz);  CHKERRQ(ierr);
      }
      else if (ii == 26){ // 2nd order derivative due to shear modulus in z-direction of fibre
        Trans_Iso_Rhs_rs_PSFEM my_fe_k_rs_EpGzp(mField,Aij,D,ddF,"DISPLACEMENT","DISP_r_Gzp","YoungP","ShearZP", "transversely_isotropic", "reinforcement");
        ierr = mField.loop_finite_elements("STOCHASIC_PROBLEM","TRAN_ISOTROPIC_ELASTIC",my_fe_k_rs_EpGzp);  CHKERRQ(ierr);
      }
      else if (ii == 27){ // 2nd order derivative due to shear modulus in z-direction of fibre
        Trans_Iso_Rhs_rs_PSFEM my_fe_k_rs_EpEm(mField,Aij,D,ddF,"DISPLACEMENT","DISP_r_Em","YoungP","YoungM", "transversely_isotropic", "reinforcement");
        ierr = mField.loop_finite_elements("STOCHASIC_PROBLEM","TRAN_ISOTROPIC_ELASTIC",my_fe_k_rs_EpEm);  CHKERRQ(ierr);
      }
      else if (ii == 28){ // 2nd order derivative due to shear modulus in z-direction of fibre
        Trans_Iso_Rhs_rs_PSFEM my_fe_k_rs_EpPm(mField,Aij,D,ddF,"DISPLACEMENT","DISP_r_Pm","YoungP","PoissonM", "transversely_isotropic", "reinforcement");
        ierr = mField.loop_finite_elements("STOCHASIC_PROBLEM","TRAN_ISOTROPIC_ELASTIC",my_fe_k_rs_EpPm);  CHKERRQ(ierr);
      }
      else if (ii == 29){ // 2nd order derivative due to shear modulus in z-direction of fibre
        Trans_Iso_Rhs_rs_PSFEM my_fe_k_rs_EzGzp(mField,Aij,D,ddF,"DISPLACEMENT","DISP_r_Gzp","YoungZ","ShearZP", "transversely_isotropic", "reinforcement");
        ierr = mField.loop_finite_elements("STOCHASIC_PROBLEM","TRAN_ISOTROPIC_ELASTIC",my_fe_k_rs_EzGzp);  CHKERRQ(ierr);
      }
      else if (ii == 30){ // 2nd order derivative due to shear modulus in z-direction of fibre
        Trans_Iso_Rhs_rs_PSFEM my_fe_k_rs_EzEm(mField,Aij,D,ddF,"DISPLACEMENT","DISP_r_Em","YoungZ","YoungM", "transversely_isotropic", "reinforcement");
        ierr = mField.loop_finite_elements("STOCHASIC_PROBLEM","TRAN_ISOTROPIC_ELASTIC",my_fe_k_rs_EzEm);  CHKERRQ(ierr);
      }
      else if (ii == 31){ // 2nd order derivative due to shear modulus in z-direction of fibre
        Trans_Iso_Rhs_rs_PSFEM my_fe_k_rs_EzPm(mField,Aij,D,ddF,"DISPLACEMENT","DISP_r_Pm","YoungZ","PoissonM", "transversely_isotropic", "reinforcement");
        ierr = mField.loop_finite_elements("STOCHASIC_PROBLEM","TRAN_ISOTROPIC_ELASTIC",my_fe_k_rs_EzPm);  CHKERRQ(ierr);
      }
      else if (ii == 32){ // 2nd order derivative due to shear modulus in z-direction of fibre
        Trans_Iso_Rhs_rs_PSFEM my_fe_k_rs_GzpEm(mField,Aij,D,ddF,"DISPLACEMENT","DISP_r_Em","ShearZP","YoungM", "transversely_isotropic", "reinforcement");
        ierr = mField.loop_finite_elements("STOCHASIC_PROBLEM","TRAN_ISOTROPIC_ELASTIC",my_fe_k_rs_GzpEm);  CHKERRQ(ierr);
      }
      else if (ii == 33){ // 2nd order derivative due to shear modulus in z-direction of fibre
        Trans_Iso_Rhs_rs_PSFEM my_fe_k_rs_GzpPm(mField,Aij,D,ddF,"DISPLACEMENT","DISP_r_Pm","ShearZP","PoissonM", "transversely_isotropic", "reinforcement");
        ierr = mField.loop_finite_elements("STOCHASIC_PROBLEM","TRAN_ISOTROPIC_ELASTIC",my_fe_k_rs_GzpPm);  CHKERRQ(ierr);
      }
      else if (ii == 34){ // 2nd order derivative due to shear modulus in z-direction of fibre
        Trans_Iso_Rhs_rs_PSFEM my_fe_k_rs_EmPm(mField,Aij,D,ddF,"DISPLACEMENT","DISP_r_Pm","Young","Poisson", "isotropic", "matrix");
        ierr = mField.loop_finite_elements("STOCHASIC_PROBLEM","ELASTIC",my_fe_k_rs_EmPm);  CHKERRQ(ierr);
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
      ierr = VecCreateMPI(PETSC_COMM_WORLD, 6, 6*pcomm->size(), &Stress_Homo_r);  CHKERRQ(ierr);
      ierr = VecZeroEntries(Stress_Homo_r); CHKERRQ(ierr);
      
      ElasticFE_RVELagrange_Homogenized_Stress_Periodic MyFE_RVEHomoStressPeriodic_r(mField,Aij,dD,dF,&RVE_volume, applied_strain, Stress_Homo_r,"DISPLACEMENT","Lagrange_mul_disp",field_rank);
      ierr = mField.loop_finite_elements("STOCHASIC_PROBLEM","Lagrange_elem",MyFE_RVEHomoStressPeriodic_r);  CHKERRQ(ierr);
      
      
      if(pcomm->rank()==0){
        PetscScalar    *avec_r;
        VecGetArray(Stress_Homo_r, &avec_r);
        
        cout<<"\n"<<ss_field.str().c_str()<<" =\n"<<endl;
        //cout<< "\n"<<ss_field<<" = \n\n";
        for(int kk=0; kk<6; kk++){
          //cout.precision(15);
          //cout<<*avec_r<<endl;
          if (ii < nvars){
            Dmat_r(ii)(kk,ics) = *avec_r;
          }
          else {
            Dmat_rs(ii-nvars)(kk,ics) = *avec_r;
          }
          avec_r++;
        }
        VecRestoreArray(Stress_Homo_r,&avec_r);
      }
      //cout<< "\n";
    }
    
  }
  
  
  // ===========================================================================
  //
  //  VII. OUTPUT
  //
  // ===========================================================================
  
  PostProcVertexMethod ent_method(moab);
  ierr = mField.loop_dofs("STOCHASIC_PROBLEM","DISPLACEMENT",ROW,ent_method); CHKERRQ(ierr);
  
  if(pcomm->rank()==0) {
    EntityHandle out_meshset;
    rval = moab.create_meshset(MESHSET_SET,out_meshset); CHKERR_PETSC(rval);
    ierr = mField.get_problem_finite_elements_entities("STOCHASIC_PROBLEM","ELASTIC",out_meshset); CHKERRQ(ierr);
    ierr = mField.get_problem_finite_elements_entities("STOCHASIC_PROBLEM","TRAN_ISOTROPIC_ELASTIC",out_meshset); CHKERRQ(ierr);
    ierr = mField.get_problem_finite_elements_entities("STOCHASIC_PROBLEM","INTERFACE",out_meshset); CHKERRQ(ierr);
    rval = moab.write_file(outName,"VTK","",&out_meshset,1); CHKERR_PETSC(rval);
    rval = moab.delete_entities(&out_meshset,1); CHKERR_PETSC(rval);
  }
  
  //    TranIso_PostProc_FibreDirRot_OnRefMesh fe_post_proc_method( mField, LAMBDA(YoungModulusP,PoissonRatioP),MU(YoungModulusP,PoissonRatioP), YoungModulusP,YoungModulusZ,PoissonRatioP,PoissonRatioPZ,ShearModulusZP);
  //
  //    ierr = mField.loop_finite_elements("STOCHASIC_PROBLEM","ELASTIC",fe_post_proc_method);  CHKERRQ(ierr);
  //    ierr = mField.loop_finite_elements("STOCHASIC_PROBLEM","TRAN_ISOTROPIC_ELASTIC",fe_post_proc_method);  CHKERRQ(ierr);
  //
  //    PetscSynchronizedFlush(PETSC_COMM_WORLD);
  //    if(pcomm->rank()==0) {
  //        rval = fe_post_proc_method.moab_post_proc.write_file(outName2,"VTK",""); CHKERR_PETSC(rval);
  //    }
  //
  //    PostProcCohesiveForces fe_post_proc_prisms(mField,YoungModulus*alpha);
  //    ierr = mField.loop_finite_elements("STOCHASIC_PROBLEM","INTERFACE",fe_post_proc_prisms);  CHKERRQ(ierr);
  //    PetscSynchronizedFlush(PETSC_COMM_WORLD);
  //    if(pcomm->rank()==0) {
  //        rval = fe_post_proc_prisms.moab_post_proc.write_file("out_post_proc_prisms.vtk","VTK",""); CHKERR_PETSC(rval);
  //    }
  
  ofstream TheFile;
  TheFile.open("//mnt//home//Dropbox//DURACOMP_Cal//009_MoFEM//00_Post_Processing//Result_PSFE_Periodic.txt",ofstream::out);
  
  // Zero-order derivative of Dmat
  for (int irow=0; irow < 6; irow++) {
    for (int icol=0; icol < 6; icol++) {
      TheFile<<setprecision(15)<<Dmat(irow,icol)<<"\t";
    }
    TheFile<<"\n";
  }
  // First-order derivative of Dmat
  for (int irv = 0; irv < nvars; irv++) {
    for (int irow=0; irow < 6; irow++) {
      for (int icol=0; icol < 6; icol++) {
        TheFile<<setprecision(15)<<Dmat_r(irv)(irow,icol)<<"\t";
      }
      TheFile<<"\n";
    }
  }
  // Second-order derivative of Dmat
  for (int irv = 0; irv < nvars; irv++) {
    for (int irow=0; irow < 6; irow++) {
      for (int icol=0; icol < 6; icol++) {
        TheFile<<setprecision(15)<<Dmat_rs(irv)(irow,icol)<<"\t";
      }
      TheFile<<"\n";
    }
  }
  TheFile.close();

  
  // ===========================================================================
  //
  //  VIII. FINISH
  //
  // ===========================================================================
  //detroy matrices
  for(int ii = 0;ii<6;ii++) {
    ierr = VecDestroy(&F[ii]); CHKERRQ(ierr);
  }
  
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
