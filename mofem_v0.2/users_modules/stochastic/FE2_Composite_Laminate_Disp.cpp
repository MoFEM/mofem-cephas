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
#include <FE2_ElasticFEMethod.hpp>

#include "ElasticFE_RVELagrange_Disp_Multi_Rhs.hpp"
#include "ElasticFE_RVELagrange_Homogenized_Stress_Disp.hpp"
#include "RVEVolume.hpp"

using namespace boost::numeric;
using namespace ObosleteUsersModules;

ErrorCode rval;
PetscErrorCode ierr;

static char help[] = "...\n\n";

const double young_modulus = 1;
const double poisson_ratio = 0.0;


PetscErrorCode Dmat_Transfermation(double theta, ublas::matrix<FieldData> &Dmat_in, ublas::matrix<FieldData> &Dmat_out) {
  PetscFunctionBegin;
  
  double l1, l2, l3, m1, m2, m3, n1, n2, n3;
  l1=cos(theta);
  m1=sin(theta);
  n1=0.0;
  l2=-sin(theta);
  m2=cos(theta);
  n2=0.0;
  l3=0.0;
  m3=0.0;
  n3=1.0;
  
  //  cout<<"l1 = "<<l1<<endl;
  //  cout<<"m1 = "<<m1<<endl;
  //  cout<<"l2 = "<<l2<<endl;
  //  cout<<"m2 = "<<m2<<endl;
  
  ublas::matrix<FieldData> T_strain;    T_strain.resize(6,6);
  T_strain(0,0)=l1*l1;    T_strain(0,1)=m1*m1;   T_strain(0,2)=n1*n1;    T_strain(0,3)=l1*m1;        T_strain(0,4)=l1*n1;        T_strain(0,5)=m1*n1;
  T_strain(1,0)=l2*l2;    T_strain(1,1)=m2*m2;   T_strain(1,2)=n2*n2;    T_strain(1,3)=l2*m2;        T_strain(1,4)=l2*n2;        T_strain(1,5)=m2*n2;
  T_strain(2,0)=l3*l3;    T_strain(2,1)=m3*m3;   T_strain(2,2)=n3*n3;    T_strain(2,3)=l3*m3;        T_strain(2,4)=l3*n3;        T_strain(2,5)=m3*n3;
  T_strain(3,0)=2*l1*l2;  T_strain(3,1)=2*m1*m2; T_strain(3,2)=2*n1*n2;  T_strain(3,3)=l1*m2+m1*l2;  T_strain(3,4)=l1*n2+n1*l2;  T_strain(3,5)=m1*n2+n1*m2;
  T_strain(4,0)=2*l1*l3;  T_strain(4,1)=2*m1*m3; T_strain(4,2)=2*n1*n3;  T_strain(4,3)=l1*m3+m1*l3;  T_strain(4,4)=l1*n3+n1*l3;  T_strain(4,5)=m1*n3+n1*m3;
  T_strain(5,0)=2*l2*l3;  T_strain(5,1)=2*m2*m3; T_strain(5,2)=2*n2*n3;  T_strain(5,3)=l2*m3+m2*l3;  T_strain(5,4)=l2*n3+n2*l3;  T_strain(5,5)=m2*n3+n2*m3;
//  cout<<"\n\nT_strain = "<<T_strain<<endl;
  
  ublas::matrix<FieldData> Mat1=prod(Dmat_in,T_strain);
  Dmat_out=prod(trans(T_strain), Mat1);
//  cout<<"\n\n Dmat_out = "<<Dmat_out<<endl;

  PetscFunctionReturn(0);
}

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
  
  vector<Range> RangeFibre(noOfFibres);
  vector<EntityHandle> fibre_meshset(noOfFibres);
  
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
  //FE Transverse Isotropic
  ierr = m_field_RVE.modify_finite_element_add_field_row("TRAN_ISO_FE_RVE","DISP_RVE"); CHKERRQ(ierr);
  ierr = m_field_RVE.modify_finite_element_add_field_col("TRAN_ISO_FE_RVE","DISP_RVE"); CHKERRQ(ierr);
  ierr = m_field_RVE.modify_finite_element_add_field_data("TRAN_ISO_FE_RVE","DISP_RVE"); CHKERRQ(ierr);
  ierr = m_field_RVE.modify_finite_element_add_field_data("TRAN_ISO_FE_RVE","POTENTIAL_FIELD"); CHKERRQ(ierr);
  
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
   unstructured tetrahedral meshes. IJNME, 58(14). pp.2103-2130.
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
  
  
  // ===========================================================================
  //
  //  A.VI. SOLUTION PHASE
  //
  // ===========================================================================
  
  
  /*****************************************************************************
   *
   *  0. PREPARATION FOR PROCESSING SOLVE
   *
   ****************************************************************************/
  Vec F1,F2,F3,F4,F5,F6,D1,D2,D3,D4,D5,D6;
  ierr = m_field_RVE.VecCreateGhost("ELASTIC_PROBLEM_RVE",ROW,&F1); CHKERRQ(ierr);
  ierr = m_field_RVE.VecCreateGhost("ELASTIC_PROBLEM_RVE",ROW,&F2); CHKERRQ(ierr);
  ierr = m_field_RVE.VecCreateGhost("ELASTIC_PROBLEM_RVE",ROW,&F3); CHKERRQ(ierr);
  ierr = m_field_RVE.VecCreateGhost("ELASTIC_PROBLEM_RVE",ROW,&F4); CHKERRQ(ierr);
  ierr = m_field_RVE.VecCreateGhost("ELASTIC_PROBLEM_RVE",ROW,&F5); CHKERRQ(ierr);
  ierr = m_field_RVE.VecCreateGhost("ELASTIC_PROBLEM_RVE",ROW,&F6); CHKERRQ(ierr);
  
  ierr = m_field_RVE.VecCreateGhost("ELASTIC_PROBLEM_RVE",COL,&D1); CHKERRQ(ierr);
  ierr = m_field_RVE.VecCreateGhost("ELASTIC_PROBLEM_RVE",COL,&D2); CHKERRQ(ierr);
  ierr = m_field_RVE.VecCreateGhost("ELASTIC_PROBLEM_RVE",COL,&D3); CHKERRQ(ierr);
  ierr = m_field_RVE.VecCreateGhost("ELASTIC_PROBLEM_RVE",COL,&D4); CHKERRQ(ierr);
  ierr = m_field_RVE.VecCreateGhost("ELASTIC_PROBLEM_RVE",COL,&D5); CHKERRQ(ierr);
  ierr = m_field_RVE.VecCreateGhost("ELASTIC_PROBLEM_RVE",COL,&D6); CHKERRQ(ierr);
  
  
  /*****************************************************************************
   *
   *  1. Assembling global stiffness matrix K
   *     and external force vector F
   ****************************************************************************/
  Mat Aij;
  ierr = m_field_RVE.MatCreateMPIAIJWithArrays("ELASTIC_PROBLEM_RVE",&Aij); CHKERRQ(ierr);
  
  struct MyElasticFEMethod: public ElasticFEMethod {
    MyElasticFEMethod(FieldInterface& _m_field,
                      Mat& _Aij,Vec& _D,Vec& _F,double _lambda,double _mu, string _field_name = "DISPLACEMENT"):
    ElasticFEMethod(_m_field,_Aij,_D,_F,_lambda,_mu,_field_name) {};
    
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
  //const double young_modulus = 1;
  //const double poisson_ratio = 0.0;
  
  double YoungModulus = 1;
  double PoissonRatio = 0.0;
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
  
  RVEVolume MyRVEVol(m_field_RVE,Aij,D1,F1,LAMBDA(YoungModulus,PoissonRatio),MU(YoungModulus,PoissonRatio), RVE_volume_Vec);
  RVEVolumeTrans MyRVEVolTrans(m_field_RVE,Aij,D1,F1, RVE_volume_Vec);
  
  ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","ELASTIC_FE_RVE",MyRVEVol);  CHKERRQ(ierr);
  ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",MyRVEVolTrans);  CHKERRQ(ierr);
  //    ierr = VecView(RVE_volume_Vec,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);
  ierr = VecSum(RVE_volume_Vec, &RVE_volume);  CHKERRQ(ierr);
  cout<<"Final RVE_volume = "<< RVE_volume <<endl;
  
  for(_IT_CUBITMESHSETS_BY_SET_TYPE_FOR_LOOP_(m_field_RVE,BLOCKSET,it))
  {
    cout << endl << *it << endl;
    
    //Get block name
    string name = it->get_name();
    
    if (name.compare(0,20,"MAT_ELASTIC_TRANSISO") == 0) {
      Mat_Elastic_TransIso mydata;
      ierr = it->get_attribute_data_structure(mydata); CHKERRQ(ierr);
      /*double Youngp, Youngz, Poissonp, Poissonpz, Shearzp;
       Youngp    = mydata.data.Youngp;
       Youngz    = mydata.data.Youngz;
       Poissonp  = mydata.data.Poissonp;
       Poissonpz = mydata.data.Poissonpz;
       Shearzp   = mydata.data.Shearzp;
       
       mydata.data.Youngp    = Youngp;
       mydata.data.Youngz    = Youngz;
       mydata.data.Poissonp  = Poissonp;
       mydata.data.Poissonpz = Poissonpz;
       mydata.data.Shearzp   = Shearzp - 0.1;
       ierr = it->set_attribute_data_structure(mydata); CHKERRQ(ierr);*/
      cout << mydata;
    }
    if (name.compare(0,13,"MAT_ELASTIC_1") == 0) {
      Mat_Elastic mydata;
      ierr = it->get_attribute_data_structure(mydata); CHKERRQ(ierr);
      YoungModulus = mydata.data.Young;
      PoissonRatio = mydata.data.Poisson;
      //mydata.data.Young = YoungModulus;
      //mydata.data.Poisson = PoissonRatio - 0.001;
      //ierr = it->set_attribute_data_structure(mydata); CHKERRQ(ierr);
      cout << mydata;
    }
  }
  
  cout<<YoungModulus<<"\t"<<PoissonRatio<<endl;
  
  ublas::vector<FieldData> applied_strain;  //it is not used in the calculation, it is required by ElasticFE_RVELagrange_Disp as input
  applied_strain.resize(1.5*field_rank+1.5); applied_strain.clear();
  
  MyElasticFEMethod my_fe(m_field_RVE,Aij,D1,F1,LAMBDA(YoungModulus,PoissonRatio),MU(YoungModulus,PoissonRatio),"DISP_RVE");
  TranIsotropicFibreDirRotElasticFEMethod MyTIsotFE(m_field_RVE,Aij,D1,F1,"DISP_RVE");
  ElasticFE_RVELagrange_Disp_Multi_Rhs MyFE_RVELagrange(m_field_RVE,Aij,D1,F1,F2,F3,F4,F5,F6,applied_strain,"DISP_RVE","Lagrange_mul_disp",field_rank);
  
  cout<<"After ElasticFE_RVELagrange_Disp_Multi_Rhs = "<<endl;
  
  ierr = VecZeroEntries(F1); CHKERRQ(ierr);
  ierr = VecZeroEntries(F2); CHKERRQ(ierr);
  ierr = VecZeroEntries(F3); CHKERRQ(ierr);
  ierr = VecZeroEntries(F4); CHKERRQ(ierr);
  ierr = VecZeroEntries(F5); CHKERRQ(ierr);
  ierr = VecZeroEntries(F6); CHKERRQ(ierr);
  ierr = MatZeroEntries(Aij); CHKERRQ(ierr);
  ierr = VecZeroEntries(D1); CHKERRQ(ierr);
  ierr = m_field_RVE.set_global_VecCreateGhost("ELASTIC_PROBLEM_RVE",ROW,D1,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
  
  ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","ELASTIC_FE_RVE",my_fe);  CHKERRQ(ierr);
  ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",MyTIsotFE);  CHKERRQ(ierr);
  ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","Lagrange_FE",MyFE_RVELagrange);  CHKERRQ(ierr);
  
  
  ierr = MatAssemblyBegin(Aij,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
  ierr = MatAssemblyEnd(Aij,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
  
  ierr = VecAssemblyBegin(F1); CHKERRQ(ierr);
  ierr = VecAssemblyEnd(F1); CHKERRQ(ierr);
  
  ierr = VecAssemblyBegin(F2); CHKERRQ(ierr);
  ierr = VecAssemblyEnd(F2); CHKERRQ(ierr);
  
  ierr = VecAssemblyBegin(F3); CHKERRQ(ierr);
  ierr = VecAssemblyEnd(F3); CHKERRQ(ierr);
  
  ierr = VecAssemblyBegin(F4); CHKERRQ(ierr);
  ierr = VecAssemblyEnd(F4); CHKERRQ(ierr);
  
  ierr = VecAssemblyBegin(F5); CHKERRQ(ierr);
  ierr = VecAssemblyEnd(F5); CHKERRQ(ierr);
  
  ierr = VecAssemblyBegin(F6); CHKERRQ(ierr);
  ierr = VecAssemblyEnd(F6); CHKERRQ(ierr);
  
  
  /*****************************************************************************
   *
   *  3. SOLVE THE FINITE ELEMENT EQUILIBRIUM EQUATION
   *     [K][U] = [F]
   *
   ****************************************************************************/
  
  
  //Solver
  KSP solver;
  ierr = KSPCreate(PETSC_COMM_WORLD,&solver); CHKERRQ(ierr);
  ierr = KSPSetOperators(solver,Aij,Aij); CHKERRQ(ierr);
  ierr = KSPSetFromOptions(solver); CHKERRQ(ierr);
  ierr = KSPSetUp(solver); CHKERRQ(ierr);
  
  //----------------------------------------------------------------------------
  // 3.1 Solving the equation to get nodal displacement
  //     case 1: applied macro strain: [1 0 0 0 0 0]^T
  //----------------------------------------------------------------------------
  //solve for F2 and D2
  ierr = KSPSolve(solver,F1,D1); CHKERRQ(ierr);
  ierr = m_field_RVE.set_global_ghost_vector("ELASTIC_PROBLEM_RVE",ROW,D1,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
  
  ublas::matrix<FieldData> Dnmat;
  Dnmat.resize(6,6); Dnmat.clear();
  
  //create a vector for 6 components of homogenized stress
  Vec Stress_Homo;
  if(pcomm->rank()==0) {
    VecCreateGhost(PETSC_COMM_WORLD,6,6,0,PETSC_NULL,&Stress_Homo);
  } else {
    int ghost[] = {0,1,2,3,4,5};
    VecCreateGhost(PETSC_COMM_WORLD,0,6,6,ghost,&Stress_Homo);
    
  }
  
  {
    ierr = VecZeroEntries(Stress_Homo); CHKERRQ(ierr);
    
    ElasticFE_RVELagrange_Homogenized_Stress_Disp MyFE_RVEHomoStressDisp(m_field_RVE,Aij,D1,F1,&RVE_volume, applied_strain, Stress_Homo,"DISP_RVE","Lagrange_mul_disp",field_rank);
    ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","Lagrange_FE",MyFE_RVEHomoStressDisp);  CHKERRQ(ierr);
    
    PetscScalar *avec;
    VecGetArray(Stress_Homo, &avec);
    for(int ii=0; ii<6; ii++){
      Dnmat(ii,0)=*avec;
      avec++;
    }
    
    if(pcomm->rank()==0){
      cout<< "\nStress_Homo = \n\n";
      cout.precision(15);
      for(int ii=0; ii<6; ii++){
        cout <<Dnmat(ii,0)<<endl;
      }
    }
  }
  //----------------------------------------------------------------------------
  // 3.1 Solving the equation to get nodal displacement
  //     case 2: applied macro strain: [0 1 0 0 0 0]^T
  //----------------------------------------------------------------------------
  //solve for F2 and D2
  ierr = KSPSolve(solver,F2,D2); CHKERRQ(ierr);
  ierr = VecGhostUpdateBegin(D2,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = VecGhostUpdateEnd(D2,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = m_field_RVE.set_global_ghost_vector("ELASTIC_PROBLEM_RVE",ROW,D2,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
  {
    // Extract homogenized stress
    ierr = VecZeroEntries(Stress_Homo); CHKERRQ(ierr);
    
    ElasticFE_RVELagrange_Homogenized_Stress_Disp MyFE_RVEHomoStressDisp(m_field_RVE,Aij,D2,F2,&RVE_volume, applied_strain, Stress_Homo,"DISP_RVE","Lagrange_mul_disp",field_rank);
    ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","Lagrange_FE",MyFE_RVEHomoStressDisp);  CHKERRQ(ierr);
    
    PetscScalar *avec;
    VecGetArray(Stress_Homo, &avec);
    for(int ii=0; ii<6; ii++){
      Dnmat(ii,1)=*avec;
      avec++;
    }
    
    if(pcomm->rank()==0){
      cout<< "\nStress_Homo = \n\n";
      for(int ii=0; ii<6; ii++){
        cout <<Dnmat(ii,1)<<endl;
      }
    }
  }
  //----------------------------------------------------------------------------
  // 3.1 Solving the equation to get nodal displacement
  //     case 3: applied macro strain: [0 0 1 0 0 0]^T
  //----------------------------------------------------------------------------
  //solve for F3 and D3
  ierr = KSPSolve(solver,F3,D3); CHKERRQ(ierr);
  ierr = VecGhostUpdateBegin(D3,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = VecGhostUpdateEnd(D3,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = m_field_RVE.set_global_ghost_vector("ELASTIC_PROBLEM_RVE",ROW,D3,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
  {
    // Extract homogenized stress
    ierr = VecZeroEntries(Stress_Homo); CHKERRQ(ierr);
    ElasticFE_RVELagrange_Homogenized_Stress_Disp MyFE_RVEHomoStressDisp(m_field_RVE,Aij,D3,F3,&RVE_volume, applied_strain, Stress_Homo,"DISP_RVE","Lagrange_mul_disp",field_rank);
    ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","Lagrange_FE",MyFE_RVEHomoStressDisp);  CHKERRQ(ierr);
    
    PetscScalar *avec;
    VecGetArray(Stress_Homo, &avec);
    for(int ii=0; ii<6; ii++){
      Dnmat(ii,2)=*avec;
      avec++;
    }
    
    if(pcomm->rank()==0){
      cout<< "\nStress_Homo = \n\n";
      for(int ii=0; ii<6; ii++){
        cout <<Dnmat(ii,2)<<endl;
      }
    }
  }
  
  //----------------------------------------------------------------------------
  // 3.1 Solving the equation to get nodal displacement
  //     case 4: applied macro strain: [0 0 0 1 0 0]^T
  //----------------------------------------------------------------------------
  //solve for F4 and D4
  ierr = KSPSolve(solver,F4,D4); CHKERRQ(ierr);
  ierr = VecGhostUpdateBegin(D4,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = VecGhostUpdateEnd(D4,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = m_field_RVE.set_global_ghost_vector("ELASTIC_PROBLEM_RVE",ROW,D4,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
  {
    // Extract homogenized stress
    ierr = VecZeroEntries(Stress_Homo); CHKERRQ(ierr);
    ElasticFE_RVELagrange_Homogenized_Stress_Disp MyFE_RVEHomoStressDisp(m_field_RVE,Aij,D4,F4,&RVE_volume, applied_strain, Stress_Homo,"DISP_RVE","Lagrange_mul_disp",field_rank);
    ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","Lagrange_FE",MyFE_RVEHomoStressDisp);  CHKERRQ(ierr);
    
    PetscScalar *avec;
    VecGetArray(Stress_Homo, &avec);
    for(int ii=0; ii<6; ii++){
      Dnmat(ii,3)=*avec;
      avec++;
    }
    
    if(pcomm->rank()==0){
      cout<< "\nStress_Homo = \n\n";
      for(int ii=0; ii<6; ii++){
        cout <<Dnmat(ii,3)<<endl;
      }
    }
  }
  
  //----------------------------------------------------------------------------
  // 3.1 Solving the equation to get nodal displacement
  //     case 5: applied macro strain: [0 0 0 0 1 0]^T
  //----------------------------------------------------------------------------
  //solve for F5 and D5
  ierr = KSPSolve(solver,F5,D5); CHKERRQ(ierr);
  ierr = VecGhostUpdateBegin(D5,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = VecGhostUpdateEnd(D5,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = m_field_RVE.set_global_ghost_vector("ELASTIC_PROBLEM_RVE",ROW,D5,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
  {
    // Extract homogenized stress
    ierr = VecZeroEntries(Stress_Homo); CHKERRQ(ierr);
    ElasticFE_RVELagrange_Homogenized_Stress_Disp MyFE_RVEHomoStressDisp(m_field_RVE,Aij,D5,F5,&RVE_volume, applied_strain, Stress_Homo,"DISP_RVE","Lagrange_mul_disp",field_rank);
    ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","Lagrange_FE",MyFE_RVEHomoStressDisp);  CHKERRQ(ierr);
    
    PetscScalar *avec;
    VecGetArray(Stress_Homo, &avec);
    for(int ii=0; ii<6; ii++){
      Dnmat(ii,4)=*avec;
      avec++;
    }
    
    if(pcomm->rank()==0){
      cout<< "\nStress_Homo = \n\n";
      for(int ii=0; ii<6; ii++){
        cout <<Dnmat(ii,4)<<endl;
      }
    }
  }
  
  //----------------------------------------------------------------------------
  // 3.1 Solving the equation to get nodal displacement
  //     case 6: applied macro strain: [0 0 0 0 0 1]^T
  //----------------------------------------------------------------------------
  //solve for F6 and D6
  ierr = KSPSolve(solver,F6,D6); CHKERRQ(ierr);
  ierr = VecGhostUpdateBegin(D6,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = VecGhostUpdateEnd(D6,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = m_field_RVE.set_global_ghost_vector("ELASTIC_PROBLEM_RVE",ROW,D6,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
  {
    // Extract homogenized stress
    ierr = VecZeroEntries(Stress_Homo); CHKERRQ(ierr);
    ElasticFE_RVELagrange_Homogenized_Stress_Disp MyFE_RVEHomoStressDisp(m_field_RVE,Aij,D6,F6,&RVE_volume, applied_strain, Stress_Homo,"DISP_RVE","Lagrange_mul_disp",field_rank);
    ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","Lagrange_FE",MyFE_RVEHomoStressDisp);  CHKERRQ(ierr);
    
    //VecGhostUpdateBegin(Stress_Homo,INSERT_VALUES,SCATTER_FORWARD);
    //VecGhostUpdateEnd(Stress_Homo,INSERT_VALUES,SCATTER_FORWARD);
    PetscScalar *avec;
    VecGetArray(Stress_Homo, &avec);
    for(int ii=0; ii<6; ii++){
      Dnmat(ii,5)=*avec;
      avec++;
    }
    
    if(pcomm->rank()==0){
      cout<< "\nStress_Homo = \n\n";
      for(int ii=0; ii<6; ii++){
        cout <<Dnmat(ii,5)<<endl;
      }
      cout<< "\n\n";
      cout<< "Dnmat = "<<Dnmat<<endl;
    }
  }
  
  
  /*//Reading and writing binary files
  if(pcomm->rank()==0){
    int fd;
    PetscViewer view_out;
    PetscViewerBinaryOpen(PETSC_COMM_WORLD,"input.dat",FILE_MODE_WRITE,&view_out);
    PetscViewerBinaryGetDescriptor(view_out,&fd);
    PetscBinaryWrite(fd,&Dnmat(0,0),36,PETSC_DOUBLE,PETSC_FALSE);
    PetscViewerDestroy(&view_out);
  }
  //  Vec F1,F2,F3,F4,F5,F6,D1,D2,D3,D4,D5,D6;
  
  
  ublas::matrix<FieldData> Dmat1;
  Dmat1.resize(6,6); Dmat1.clear();
  cout<< "Dmat1 Before Reading= "<<Dmat1<<endl;
  if(pcomm->rank()==0){
    int fd;
    PetscViewer view_in;
    PetscViewerBinaryOpen(PETSC_COMM_WORLD,"input.dat",FILE_MODE_READ,&view_in);
    PetscViewerBinaryGetDescriptor(view_in,&fd);
    PetscBinaryRead(fd,&Dmat1(0,0),36,PETSC_DOUBLE);
    PetscViewerDestroy(&view_in);
  }
  cout<< "Dmat1 After Reading= "<<Dmat1<<endl;*/
  
  
  //Destroy matrices and vectors
  ierr = VecDestroy(&F1); CHKERRQ(ierr);
  ierr = VecDestroy(&F2); CHKERRQ(ierr);
  ierr = VecDestroy(&F3); CHKERRQ(ierr);
  ierr = VecDestroy(&F4); CHKERRQ(ierr);
  ierr = VecDestroy(&F5); CHKERRQ(ierr);
  ierr = VecDestroy(&F6); CHKERRQ(ierr);
  
  ierr = VecDestroy(&D1); CHKERRQ(ierr);
  ierr = VecDestroy(&D2); CHKERRQ(ierr);
  ierr = VecDestroy(&D3); CHKERRQ(ierr);
  ierr = VecDestroy(&D4); CHKERRQ(ierr);
  ierr = VecDestroy(&D5); CHKERRQ(ierr);
  ierr = VecDestroy(&D6); CHKERRQ(ierr);
  
  ierr = MatDestroy(&Aij); CHKERRQ(ierr);
  ierr = KSPDestroy(&solver); CHKERRQ(ierr);
  
  
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
  

  /*****************************************************************************
   *
   * Read parameters from line command
   *
   ****************************************************************************/
  //PetscBool flg = PETSC_TRUE;
  //char mesh_file_name[255];
  
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
  //const char *option;
  //option = "";//"PARALLEL=BCAST;";//;DEBUG_IO";
  
  rval = moab_Macro.load_file(mesh_file_name, 0, option); CHKERR_PETSC(rval);
  
  //ParallelComm* pcomm_macro = ParallelComm::get_pcomm(&moab_Macro,MYPCOMM_INDEX);
  //if(pcomm_macro == NULL) pcomm_macro =  new ParallelComm(&moab_Macro,PETSC_COMM_WORLD);
  
  
  //Create MoFEM (Joseph) database
  MoFEM::Core core_Macro(moab_Macro);
  FieldInterface& m_field_Macro = core_Macro;
  
  //ref meshset ref level 0
  ierr = m_field_Macro.seed_ref_level_3D(0,0); CHKERRQ(ierr);
  
  // stl::bitset see for more details
  BitRefLevel bit_level0;
  bit_level0.set(0);
  
  EntityHandle meshset_level0;
  rval = moab_Macro.create_meshset(MESHSET_SET,meshset_level0); CHKERR_PETSC(rval);
  ierr = m_field_Macro.seed_ref_level_3D(0,bit_level0); CHKERRQ(ierr);
  ierr = m_field_Macro.get_entities_by_ref_level(bit_level0,BitRefLevel().set(),meshset_level0); CHKERRQ(ierr);
  
  /*****************************************************************************
   *
   * Read the saved Dmat mechancial (from the computational homgenisaiton of the 0deg RVE)
   *
   ****************************************************************************/
  ublas::matrix<FieldData> Dmat_0deg;
  Dmat_0deg.resize(6,6);  Dmat_0deg.clear();

  cout<<"Dmat_0deg Before =  "<<Dmat_0deg<<endl;
  Dmat_0deg = Dnmat;

  //Find Dmat_90deg by rotating Dmat_0deg at an angle of +90 deg about the z-axis
  double theta; theta=90*(M_PI/180.0); //rotation angle about the Z-axis
  ublas::matrix<FieldData> Dmat_90deg;  Dmat_90deg.resize(6,6);   Dmat_90deg.clear();
  ierr = Dmat_Transfermation(theta, Dmat_0deg, Dmat_90deg); CHKERRQ(ierr);
  cout<<"\n\nDmat_90deg = "<<Dmat_90deg<<endl;

  //Find Dmat_pos25deg by rotating Dmat_pos25deg at an angle of +25 deg about the z-axis
  theta=25*(M_PI/180.0); //rotation angle about the Z-axis
  ublas::matrix<FieldData> Dmat_pos25deg;  Dmat_pos25deg.resize(6,6);   Dmat_pos25deg.clear();
  ierr = Dmat_Transfermation(theta, Dmat_0deg, Dmat_pos25deg); CHKERRQ(ierr);
  cout<<"\n\n Dmat_pos25deg = "<<Dmat_pos25deg<<endl;

  //Find Dmat_neg25deg by rotating Dmat_pos25deg at an angle of +25 deg about the z-axis
  theta=-25*(M_PI/180.0); //rotation angle about the Z-axis
  ublas::matrix<FieldData> Dmat_neg25deg;  Dmat_neg25deg.resize(6,6);   Dmat_neg25deg.clear();
  ierr = Dmat_Transfermation(theta, Dmat_0deg, Dmat_neg25deg); CHKERRQ(ierr);
  cout<<"\n\n Dmat_neg25deg = "<<Dmat_neg25deg<<endl;

  
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
  
  
  /*****************************************************************************
   *
   * Create finite element for the defined fields
   *
   ****************************************************************************/
  ierr = m_field_Macro.add_finite_element("ELASTIC_90deg",MF_ZERO); CHKERRQ(ierr);
  ierr = m_field_Macro.add_finite_element("ELASTIC_pos25deg",MF_ZERO); CHKERRQ(ierr);
  ierr = m_field_Macro.add_finite_element("ELASTIC_neg25deg",MF_ZERO); CHKERRQ(ierr);
  ierr = m_field_Macro.add_finite_element("ELASTIC_0_90_post_process",MF_ZERO); CHKERRQ(ierr);

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

  //Define rows/cols and element data
  ierr = m_field_Macro.modify_finite_element_add_field_row("ELASTIC_pos25deg","DISP_MACRO"); CHKERRQ(ierr);
  ierr = m_field_Macro.modify_finite_element_add_field_col("ELASTIC_pos25deg","DISP_MACRO"); CHKERRQ(ierr);
  ierr = m_field_Macro.modify_finite_element_add_field_data("ELASTIC_pos25deg","DISP_MACRO"); CHKERRQ(ierr);
  ierr = m_field_Macro.modify_finite_element_add_field_data("ELASTIC_pos25deg","MESH_NODE_POSITIONS"); CHKERRQ(ierr);
  
  //Define rows/cols and element data
  ierr = m_field_Macro.modify_finite_element_add_field_row("ELASTIC_neg25deg","DISP_MACRO"); CHKERRQ(ierr);
  ierr = m_field_Macro.modify_finite_element_add_field_col("ELASTIC_neg25deg","DISP_MACRO"); CHKERRQ(ierr);
  ierr = m_field_Macro.modify_finite_element_add_field_data("ELASTIC_neg25deg","DISP_MACRO"); CHKERRQ(ierr);
  ierr = m_field_Macro.modify_finite_element_add_field_data("ELASTIC_neg25deg","MESH_NODE_POSITIONS"); CHKERRQ(ierr);
  
  
  //Define rows/cols and element data
  ierr = m_field_Macro.modify_finite_element_add_field_row("ELASTIC_0_90_post_process","DISP_MACRO"); CHKERRQ(ierr);
  ierr = m_field_Macro.modify_finite_element_add_field_col("ELASTIC_0_90_post_process","DISP_MACRO"); CHKERRQ(ierr);
  ierr = m_field_Macro.modify_finite_element_add_field_data("ELASTIC_0_90_post_process","DISP_MACRO"); CHKERRQ(ierr);
  ierr = m_field_Macro.modify_finite_element_add_field_data("ELASTIC_0_90_post_process","MESH_NODE_POSITIONS"); CHKERRQ(ierr);

  //define problems
  ierr = m_field_Macro.add_problem("ELASTIC_PROB_MACRO"); CHKERRQ(ierr);
  
  //set finite elements for problem
  ierr = m_field_Macro.modify_problem_add_finite_element("ELASTIC_PROB_MACRO","ELASTIC_90deg"); CHKERRQ(ierr);
  ierr = m_field_Macro.modify_problem_add_finite_element("ELASTIC_PROB_MACRO","ELASTIC_pos25deg"); CHKERRQ(ierr);
  ierr = m_field_Macro.modify_problem_add_finite_element("ELASTIC_PROB_MACRO","ELASTIC_neg25deg"); CHKERRQ(ierr);
  ierr = m_field_Macro.modify_problem_add_finite_element("ELASTIC_PROB_MACRO","ELASTIC_0_90_post_process"); CHKERRQ(ierr);

  //set refinment level for problem
  ierr = m_field_Macro.modify_problem_ref_level_add_bit("ELASTIC_PROB_MACRO",bit_level0); CHKERRQ(ierr);
  
  
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

  Range TetsInBlock_90deg, TetsInBlock_pos25deg, TetsInBlock_neg25deg;
  for(_IT_CUBITMESHSETS_BY_SET_TYPE_FOR_LOOP_(m_field_Macro,BLOCKSET,it)){
    if(it->get_name() == "MAT_ELASTIC_90deg") {
      rval = moab_Macro.get_entities_by_type(it->meshset, MBTET,TetsInBlock_90deg,true); CHKERR_PETSC(rval);
    }
    if(it->get_name() == "MAT_ELASTIC_Pos25deg") {
      rval = moab_Macro.get_entities_by_type(it->meshset, MBTET,TetsInBlock_pos25deg,true); CHKERR_PETSC(rval);
    }
    if(it->get_name() == "MAT_ELASTIC_Neg25deg") {
      rval = moab_Macro.get_entities_by_type(it->meshset, MBTET,TetsInBlock_neg25deg,true); CHKERR_PETSC(rval);
    }
  }
  cout<<"===============================   TetsInBlock_90deg "<<TetsInBlock_90deg<<endl;
  cout<<"===============================   TetsInBlock_pos25deg "<<TetsInBlock_pos25deg<<endl;
  cout<<"===============================   TetsInBlock_neg25deg "<<TetsInBlock_neg25deg<<endl;
  
  
  /*****************************************************************************
   *
   * Add finite elements entities
   *
   ****************************************************************************/
  ierr = m_field_Macro.add_ents_to_finite_element_by_TETs(TetsInBlock_90deg, "ELASTIC_90deg"); CHKERRQ(ierr);
  ierr = m_field_Macro.add_ents_to_finite_element_by_TETs(TetsInBlock_pos25deg,"ELASTIC_pos25deg"); CHKERRQ(ierr);
  ierr = m_field_Macro.add_ents_to_finite_element_by_TETs(TetsInBlock_neg25deg,"ELASTIC_neg25deg"); CHKERRQ(ierr);
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
  //
  ierr = m_field_Macro.set_field_order(0,MBTET,"MESH_NODE_POSITIONS",2); CHKERRQ(ierr);
  ierr = m_field_Macro.set_field_order(0,MBTRI,"MESH_NODE_POSITIONS",2); CHKERRQ(ierr);
  ierr = m_field_Macro.set_field_order(0,MBEDGE,"MESH_NODE_POSITIONS",2); CHKERRQ(ierr);
  ierr = m_field_Macro.set_field_order(0,MBVERTEX,"MESH_NODE_POSITIONS",1); CHKERRQ(ierr);
  
  ierr = MetaNeummanForces::addNeumannBCElements(m_field_Macro,"DISP_MACRO"); CHKERRQ(ierr);
  ierr = m_field_Macro.modify_problem_add_finite_element("ELASTIC_PROB_MACRO","FORCE_FE"); CHKERRQ(ierr);
  
  
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
  ierr = m_field_Macro.build_adjacencies(bit_level0); CHKERRQ(ierr);
  
  //build problem
  ierr = m_field_Macro.build_problems(); CHKERRQ(ierr);
  
  
  // ===========================================================================
  //
  //  B.V. MESH PARTITION
  //
  // ===========================================================================
  
  //partition
  ierr = m_field_Macro.partition_problem("ELASTIC_PROB_MACRO"); CHKERRQ(ierr);
  //PetscBarrier(PETSC_NULL);
  ierr = m_field_Macro.partition_finite_elements("ELASTIC_PROB_MACRO"); CHKERRQ(ierr);
  //what are ghost nodes, see Petsc Manual
  ierr = m_field_Macro.partition_ghost_dofs("ELASTIC_PROB_MACRO"); CHKERRQ(ierr);
  
  //mField.list_dofs_by_field_name("DISP",true);
  
  //print bcs
  ierr = m_field_Macro.print_cubit_displacement_set(); CHKERRQ(ierr);
  ierr = m_field_Macro.print_cubit_force_set(); CHKERRQ(ierr);
  //print block sets with materials
  ierr = m_field_Macro.print_cubit_materials_set(); CHKERRQ(ierr);
  
 
  // ===========================================================================
  //
  //  B.VI. SOLUTION PHASE
  //
  // ===========================================================================
  
  /*****************************************************************************
   *
   *  0. PREPARATION FOR PROCESSING SOLVE
   *
   ****************************************************************************/
  //create matrices
  Vec F,D;
  ierr = m_field_Macro.VecCreateGhost("ELASTIC_PROB_MACRO",ROW,&F); CHKERRQ(ierr);
  ierr = m_field_Macro.VecCreateGhost("ELASTIC_PROB_MACRO",COL,&D); CHKERRQ(ierr);
  
  
  /*****************************************************************************
   *
   *  1. Assembling global stiffness matrix K
   *     and external force vector F
   ****************************************************************************/
  Mat A;
  ierr = m_field_Macro.MatCreateMPIAIJWithArrays("ELASTIC_PROB_MACRO",&A); CHKERRQ(ierr);
  
  
  struct MyElasticFEMethod_Macro: public FE2_ElasticFEMethod {
    MyElasticFEMethod_Macro(FieldInterface& _mField,Mat _A,Vec _D,Vec& _F, ublas::matrix<FieldData> _Dmat,string _field_name):
    FE2_ElasticFEMethod(_mField,_A,_D,_F, _Dmat, _field_name) {};
    
    PetscErrorCode Fint(Vec F_int) {
      PetscFunctionBegin;
      ierr = FE2_ElasticFEMethod::Fint(); CHKERRQ(ierr);
      for(int rr = 0;rr<row_mat;rr++) {
        if(RowGlob[rr].size()!=f_int[rr].size()) SETERRQ(PETSC_COMM_SELF,1,"data inconsistency");
        if(RowGlob[rr].size()==0) continue;
        f_int[rr] *= -1; //This is not SNES we solve K*D = -RES
        ierr = VecSetValues(F_int,RowGlob[rr].size(),&(RowGlob[rr])[0],&(f_int[rr].data()[0]),ADD_VALUES); CHKERRQ(ierr);
      }
      PetscFunctionReturn(0);
    }
    
  };
  
  
  Projection10NodeCoordsOnField ent_method_material(m_field_Macro,"MESH_NODE_POSITIONS");
  ierr = m_field_Macro.loop_dofs("MESH_NODE_POSITIONS",ent_method_material); CHKERRQ(ierr);
  
  //Assemble F and A
  DisplacementBCFEMethodPreAndPostProc my_dirichlet_bc(m_field_Macro,"DISP_MACRO",A,D,F);
  //preproc
  ierr = m_field_Macro.problem_basic_method_preProcess("ELASTIC_PROB_MACRO",my_dirichlet_bc); CHKERRQ(ierr);
  ierr = m_field_Macro.set_global_ghost_vector("ELASTIC_PROB_MACRO",ROW,D,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
//  ierr = VecView(D,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);
 
  
  MyElasticFEMethod_Macro my_fe_90deg   (m_field_Macro,A,D,F,Dmat_90deg,   "DISP_MACRO");
  MyElasticFEMethod_Macro my_fe_pos25deg(m_field_Macro,A,D,F,Dmat_pos25deg,"DISP_MACRO");
  MyElasticFEMethod_Macro my_fe_neg25deg(m_field_Macro,A,D,F,Dmat_neg25deg,"DISP_MACRO");

  
  ierr = VecZeroEntries(F); CHKERRQ(ierr);
  ierr = VecGhostUpdateBegin(F,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = VecGhostUpdateEnd(F,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = MatZeroEntries(A); CHKERRQ(ierr);
  
  //loop elems
  //PetscBarrier(PETSC_NULL);
  ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROB_MACRO","ELASTIC_90deg",my_fe_90deg);  CHKERRQ(ierr);
  ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROB_MACRO","ELASTIC_pos25deg",my_fe_pos25deg);  CHKERRQ(ierr);
  ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROB_MACRO","ELASTIC_neg25deg",my_fe_neg25deg);  CHKERRQ(ierr);

  //forces and preassures on surface
  boost::ptr_map<string,NeummanForcesSurface> neumann_forces;
  ierr = MetaNeummanForces::setNeumannFiniteElementOperators(m_field_Macro,neumann_forces,F,"DISP_MACRO"); CHKERRQ(ierr);
  boost::ptr_map<string,NeummanForcesSurface>::iterator mit = neumann_forces.begin();
  for(;mit!=neumann_forces.end();mit++) {
    ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROB_MACRO",mit->first,mit->second->getLoopFe()); CHKERRQ(ierr);
  }
  
  //postproc
  ierr = m_field_Macro.problem_basic_method_postProcess("ELASTIC_PROB_MACRO",my_dirichlet_bc); CHKERRQ(ierr);
  
  
  //set matrix possitives define and symetric for cholesky and icc preceonditionser
  ierr = MatSetOption(A,MAT_SPD,PETSC_TRUE); CHKERRQ(ierr);
  
  ierr = VecGhostUpdateBegin(F,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
  ierr = VecGhostUpdateEnd(F,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
  ierr = VecAssemblyBegin(F); CHKERRQ(ierr);
  ierr = VecAssemblyEnd(F); CHKERRQ(ierr);
//  ierr = VecView(F,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);

//  Matrix View
//  MatView(A,PETSC_VIEWER_DRAW_WORLD);//PETSC_VIEWER_STDOUT_WORLD);
//  std::string wait;
//  std::cin >> wait;

  
  /*****************************************************************************
   *
   *  3. SOLVE THE FINITE ELEMENT EQUILIBRIUM EQUATION
   *     [K][U] = [F]
   *
   ****************************************************************************/
  //Solver
  KSP solver_Macro;
  ierr = KSPCreate(PETSC_COMM_WORLD,&solver_Macro); CHKERRQ(ierr);
  ierr = KSPSetOperators(solver_Macro,A,A); CHKERRQ(ierr);
  ierr = KSPSetFromOptions(solver_Macro); CHKERRQ(ierr);
  ierr = KSPSetUp(solver_Macro); CHKERRQ(ierr);
  
  // elastic analys
  ierr = KSPSolve(solver_Macro,F,D); CHKERRQ(ierr);
  ierr = VecGhostUpdateBegin(D,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = VecGhostUpdateEnd(D,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = m_field_Macro.problem_basic_method_preProcess("ELASTIC_PROB_MACRO",my_dirichlet_bc); CHKERRQ(ierr);
//  ierr = VecView(D,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);

  
  // ===========================================================================
  //
  //  B.VII. OUTPUT
  //
  // ===========================================================================
  //Save data on mesh
  ierr = m_field_Macro.set_global_ghost_vector("ELASTIC_PROB_MACRO",ROW,D,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);

  PostPocOnRefinedMesh post_proc(m_field_Macro);
  ierr = post_proc.generateReferenceElementMesh(); CHKERRQ(ierr);
  ierr = post_proc.addFieldValuesPostProc("DISP_MACRO"); CHKERRQ(ierr);
  ierr = post_proc.addFieldValuesPostProc("MESH_NODE_POSITIONS"); CHKERRQ(ierr);
  ierr = post_proc.addFieldValuesGradientPostProc("DISP_MACRO"); CHKERRQ(ierr);
 
  //add postpocessing for sresses
  post_proc.getOpPtrVector().push_back(
    new PostPorcStressLaminates(
      m_field_Macro,
      post_proc.postProcMesh,
      post_proc.mapGaussPts,
      "DISP_MACRO",
      post_proc.commonData,Dmat_neg25deg)
  );

  
  ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROB_MACRO","ELASTIC_neg25deg",post_proc); CHKERRQ(ierr);
  rval = post_proc.postProcMesh.write_file("out_macro.h5m","MOAB","PARALLEL=WRITE_PART"); CHKERR_PETSC(rval);

  
  // ===========================================================================
  //
  //  B.VIII. FINISH
  //
  // ===========================================================================
  //Destroy matrices
  ierr = VecDestroy(&F); CHKERRQ(ierr);
  ierr = VecDestroy(&D); CHKERRQ(ierr);
  ierr = MatDestroy(&A); CHKERRQ(ierr);
  ierr = KSPDestroy(&solver); CHKERRQ(ierr);
  PetscFinalize();
  
}

