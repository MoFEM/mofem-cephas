/* Copyright (C) 2014, 
 *   Xiao-Yi Zhou (xiaoyi.zhou AT newcastle.ac.uk)
 *   Zahur Ullah (Zahur.Ullah@glasgow.ac.uk)
 * --------------------------------------------------------------
 * This routine conducts finite element implementation for stochastic multiscale
 * homogenization problem under periodic boundary condition by perturbation 
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
#include "K_r_PSFEMethod.hpp"
#include "K_rs_PSFEMethod.hpp"

#include "K_r_Young_PSFEMethod.hpp"
#include "K_r_Poisson_PSFEMethod.hpp"

#include "K_rs_Young_PSFEMethod.hpp"
#include "K_rs_EP_PSFEMethod.hpp"
#include "K_rs_Poisson_PSFEMethod.hpp"

#include "ElasticFE_RVELagrange_Periodic.hpp"
#include "ElasticFE_RVELagrange_RigidBodyTranslation.hpp"
#include "ElasticFE_RVELagrange_RigidBodyRotation.hpp"
#include "ElasticFE_RVELagrange_Homogenized_Stress_Periodic.hpp"
#include "RVEVolume.hpp"

#include "PostProcVertexMethod.hpp"
#include "PostProcDisplacementAndStrainOnRefindedMesh.hpp"

using namespace MoFEM;

ErrorCode rval;
PetscErrorCode ierr;
static char help[] = "...\n\n";


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
  // For lagrange multipliers to control the periodic motion
  ierr = mField.add_field("Lagrange_mul_disp",H1,field_rank); CHKERRQ(ierr);
  // To control the rigid body motion (3 Traslations and 3 rotations)  
  ierr = mField.add_field("Lagrange_mul_disp_rigid_trans",NOFIELD,3); CHKERRQ(ierr);
  
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
  
  
  //Define rows/cols and element data for C1 and C1T (for lagrange multipliers to contol the rigid body motion)
  //============================================================================
  //C1 row as Lagrange_mul_disp and col as DISPLACEMENT
  ierr = mField.modify_finite_element_add_field_row("Lagrange_elem_rigid_trans","Lagrange_mul_disp_rigid_trans"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_col("Lagrange_elem_rigid_trans","DISPLACEMENT"); CHKERRQ(ierr);
  
  //C1T col as Lagrange_mul_disp and row as DISPLACEMENT
  ierr = mField.modify_finite_element_add_field_col("Lagrange_elem_rigid_trans","Lagrange_mul_disp_rigid_trans"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_row("Lagrange_elem_rigid_trans","DISPLACEMENT"); CHKERRQ(ierr);
  
  //As for stress we need both displacement and temprature (Lukasz)
  ierr = mField.modify_finite_element_add_field_data("Lagrange_elem_rigid_trans","Lagrange_mul_disp_rigid_trans"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_data("Lagrange_elem_rigid_trans","DISPLACEMENT"); CHKERRQ(ierr);
  //============================================================================
  
  
  //define problems
  ierr = mField.add_problem("STOCHASIC_PROBLEM"); CHKERRQ(ierr);
  
  
  //set finite elements for problem
  ierr = mField.modify_problem_add_finite_element("STOCHASIC_PROBLEM","ELASTIC"); CHKERRQ(ierr);
  ierr = mField.modify_problem_add_finite_element("STOCHASIC_PROBLEM","Lagrange_elem"); CHKERRQ(ierr);
  ierr = mField.modify_problem_add_finite_element("STOCHASIC_PROBLEM","Lagrange_elem_rigid_trans"); CHKERRQ(ierr);
  
  
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
  
  //Add finite element to lagrange element for rigid body translation
  Range SurfacesFaces;
  ierr = mField.get_Cubit_msId_entities_by_dimension(103,SIDESET,2,SurfacesFaces,true); CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"number of SideSet 103 = %d\n",SurfacesFaces.size()); CHKERRQ(ierr);
  
  //to create meshset from range
  EntityHandle BoundFacesMeshset;
  rval = moab.create_meshset(MESHSET_SET,BoundFacesMeshset); CHKERR_PETSC(rval);
	rval = moab.add_entities(BoundFacesMeshset,SurfacesFaces); CHKERR_PETSC(rval);
  ierr = mField.seed_ref_level_MESHSET(BoundFacesMeshset,BitRefLevel().set()); CHKERRQ(ierr);

  ierr = mField.add_ents_to_finite_element_by_TRIs(SurfacesFaces,"Lagrange_elem_rigid_trans"); CHKERRQ(ierr);
  
  
  //============================================================================
  //
  // Add Periodic Prisims Between Triangles on -ve and +ve faces to 
  // implement periodic bounary conditions
  //
  //============================================================================
  
  //Populating the Multi-index container with -ve triangles
  Range SurTrisNeg;
  ierr = mField.get_Cubit_msId_entities_by_dimension(101,SIDESET,2,SurTrisNeg,true); CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"number of SideSet 101 = %d\n",SurTrisNeg.size()); CHKERRQ(ierr);
  Face_CenPos_Handle_multiIndex Face_CenPos_Handle_varNeg, Face_CenPos_Handle_varPos;
  double TriCen[3], coords_Tri[9];
  for(Range::iterator it = SurTrisNeg.begin(); it!=SurTrisNeg.end();  it++) {
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
    if(TriCen[0]>=0) TriCen[0]=double(int(TriCen[0]*1000.0+0.5))/1000.0;  else TriCen[0]=double(int(TriCen[0]*1000.0-0.5))/1000.0; //-ve and +ve value
    if(TriCen[1]>=0) TriCen[1]=double(int(TriCen[1]*1000.0+0.5))/1000.0;  else TriCen[1]=double(int(TriCen[1]*1000.0-0.5))/1000.0;
    if(TriCen[2]>=0) TriCen[2]=double(int(TriCen[2]*1000.0+0.5))/1000.0;  else TriCen[2]=double(int(TriCen[2]*1000.0-0.5))/1000.0;
    
    //        cout<<"\n\n\nTriCen[0]= "<<TriCen[0] << "   TriCen[1]= "<< TriCen[1] << "   TriCen[2]= "<< TriCen[2] <<endl;
    //fill the multi-index container with centriod coordinates and triangle handles
    Face_CenPos_Handle_varNeg.insert(Face_CenPos_Handle(TriCen[0], TriCen[1], TriCen[2], *it));
    //        for(int ii=0; ii<3; ii++) cout<<"TriCen "<<TriCen[ii]<<endl;
    //        cout<<endl<<endl;
  }
  
  //    double aaa;
  //    aaa=0.5011;
  //    cout<<"\n\n\n\n\nfloor(aaa+0.5) = "<<double(int(aaa*1000.0+0.5))/1000.0<<endl<<endl<<endl<<endl;
  //    aaa=-0.5011;
  //    cout<<"\n\n\n\n\nfloor(aaa+0.5) = "<<double(int(aaa*1000.0-0.5))/1000.0<<endl<<endl<<endl<<endl;
  
  
  //Populating the Multi-index container with +ve triangles
  Range SurTrisPos;
  ierr = mField.get_Cubit_msId_entities_by_dimension(102,SIDESET,2,SurTrisPos,true); CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"number of SideSet 102 = %d\n",SurTrisPos.size()); CHKERRQ(ierr);
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
    if(TriCen[0]>=0) TriCen[0]=double(int(TriCen[0]*1000.0+0.5))/1000.0;  else TriCen[0]=double(int(TriCen[0]*1000.0-0.5))/1000.0;
    if(TriCen[1]>=0) TriCen[1]=double(int(TriCen[1]*1000.0+0.5))/1000.0;  else TriCen[1]=double(int(TriCen[1]*1000.0-0.5))/1000.0;
    if(TriCen[2]>=0) TriCen[2]=double(int(TriCen[2]*1000.0+0.5))/1000.0;  else TriCen[2]=double(int(TriCen[2]*1000.0-0.5))/1000.0;
    //        cout<<"\n\n\nTriCen[0]= "<<TriCen[0] << "   TriCen[1]= "<< TriCen[1] << "   TriCen[2]= "<< TriCen[2] <<endl;
    
    //fill the multi-index container with centriod coordinates and triangle handles
    Face_CenPos_Handle_varPos.insert(Face_CenPos_Handle(TriCen[0], TriCen[1], TriCen[2], *it));
  }
  
  //Find minimum and maximum X, Y and Z coordinates of the RVE (using multi-index container)
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
  //    cout<<"XcoordMin "<<XcoordMin << "      XcoordMax "<<XcoordMax <<endl;
  //    cout<<"YcoordMin "<<YcoordMin << "      YcoordMax "<<YcoordMax <<endl;
  //    cout<<"ZcoordMin "<<ZcoordMin << "      ZcoordMax "<<ZcoordMax <<endl;
  
  
  //Creating Prisims between triangles on -ve and +ve faces
  typedef Face_CenPos_Handle_multiIndex::index<Tri_Hand_tag>::type::iterator Tri_Hand_iterator;
  Tri_Hand_iterator Tri_Neg;
  typedef Face_CenPos_Handle_multiIndex::index<Composite_xyzcoord>::type::iterator xyzcoord_iterator;
  xyzcoord_iterator Tri_Pos;
  Range PrismRange;
  double XPos, YPos, ZPos;
  //int count=0;
  
  //loop over -ve triangles to create prisims elemenet between +ve and -ve triangles
  for(Range::iterator it = SurTrisNeg.begin(); it!=SurTrisNeg.end();  it++) {
    
    Tri_Neg=Face_CenPos_Handle_varNeg.get<Tri_Hand_tag>().find(*it);
    //cout<<"xcoord= "<<Tri_iit->xcoord << "   ycoord= "<< Tri_iit->ycoord << "   ycoord= "<< Tri_iit->zcoord <<endl;
    
    //corresponding +ve triangle
    if(Tri_Neg->xcoord==XcoordMin){XPos=XcoordMax;              YPos=Tri_Neg->ycoord;  ZPos=Tri_Neg->zcoord;};
    if(Tri_Neg->ycoord==YcoordMin){XPos=YPos=Tri_Neg->xcoord;   YPos=YcoordMax;        ZPos=Tri_Neg->zcoord;};
    if(Tri_Neg->zcoord==ZcoordMin){XPos=YPos=Tri_Neg->xcoord;   YPos=Tri_Neg->ycoord;  ZPos=ZcoordMax;      };
    Tri_Pos=Face_CenPos_Handle_varPos.get<Composite_xyzcoord>().find(boost::make_tuple(XPos, YPos, ZPos));
    //        cout<<"Tri_Neg->xcoord= "<<Tri_Neg->xcoord << "   Tri_Neg->ycoord "<< Tri_Neg->ycoord << "   Tri_Neg->zcoord= "<< Tri_Neg->zcoord <<endl;
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
        //Round nodal coordinates to 3 dicimal places only for comparison
        //round values to 3 disimal places
        if(XNodeNeg>=0) XNodeNeg=double(int(XNodeNeg*1000.0+0.5))/1000.0;   else XNodeNeg=double(int(XNodeNeg*1000.0-0.5))/1000.0;
        if(YNodeNeg>=0) YNodeNeg=double(int(YNodeNeg*1000.0+0.5))/1000.0;   else YNodeNeg=double(int(YNodeNeg*1000.0-0.5))/1000.0;
        if(ZNodeNeg>=0) ZNodeNeg=double(int(ZNodeNeg*1000.0+0.5))/1000.0;   else ZNodeNeg=double(int(ZNodeNeg*1000.0-0.5))/1000.0;
        
        XNodePos=CoordNodePos[3*jj]; YNodePos=CoordNodePos[3*jj+1]; ZNodePos=CoordNodePos[3*jj+2];
        if(XNodePos>=0) XNodePos=double(int(XNodePos*1000.0+0.5))/1000.0;   else XNodePos=double(int(XNodePos*1000.0-0.5))/1000.0;
        if(YNodePos>=0) YNodePos=double(int(YNodePos*1000.0+0.5))/1000.0;   else YNodePos=double(int(YNodePos*1000.0-0.5))/1000.0;
        if(ZNodePos>=0) ZNodePos=double(int(ZNodePos*1000.0+0.5))/1000.0;   else ZNodePos=double(int(ZNodePos*1000.0-0.5))/1000.0;
        
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
  //    //Saving prisms in interface.vtk
  //      EntityHandle out_meshset1;
  //      rval = moab.create_meshset(MESHSET_SET,out_meshset1); CHKERR_PETSC(rval);
  //      rval = moab.add_entities(out_meshset1,PrismRange); CHKERR_PETSC(rval);
  //      rval = moab.write_file("Prisms.vtk","VTK","",&out_meshset1,1); CHKERR_PETSC(rval);
  
  
  //Adding Prisims to Element Lagrange_elem (to loop over these prisims)
  EntityHandle PrismRangeMeshset;
  rval = moab.create_meshset(MESHSET_SET,PrismRangeMeshset); CHKERR_PETSC(rval);
  rval = moab.add_entities(PrismRangeMeshset,PrismRange); CHKERR_PETSC(rval);
  ierr = mField.seed_ref_level_3D(PrismRangeMeshset,bit_level0); CHKERRQ(ierr);
  ierr = mField.get_entities_by_ref_level(bit_level0,BitRefLevel().set(),meshset_level0); CHKERRQ(ierr);
  //ierr = mField.add_ents_to_finite_element_EntType_by_bit_ref(bit_level0,"Lagrange_elem",MBPRISM); CHKERRQ(ierr);
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
  
  ierr = mField.set_field_order(0,MBTET,"DISPLACEMENT",order); CHKERRQ(ierr);
  ierr = mField.set_field_order(0,MBTRI,"DISPLACEMENT",order); CHKERRQ(ierr);
  ierr = mField.set_field_order(0,MBEDGE,"DISPLACEMENT",order); CHKERRQ(ierr);
  ierr = mField.set_field_order(0,MBVERTEX,"DISPLACEMENT",1); CHKERRQ(ierr);
  
  //ierr = mField.set_field_order(0,MBPRISM,"Lagrange_mul_disp",order); CHKERRQ(ierr);
  ierr = mField.set_field_order(0,MBTRI,"Lagrange_mul_disp",order); CHKERRQ(ierr);
  ierr = mField.set_field_order(0,MBEDGE,"Lagrange_mul_disp",order); CHKERRQ(ierr);
  ierr = mField.set_field_order(0,MBVERTEX,"Lagrange_mul_disp",1); CHKERRQ(ierr);
  
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
  
  
  //    //Matrix View
  //    MatView(Aij,PETSC_VIEWER_DRAW_WORLD);//PETSC_VIEWER_STDOUT_WORLD);
  //    std::string wait;
  //    std::cin >> wait;
  
  
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
  MyElasticFEMethod MyFE(mField,Aij,D,F,LAMBDA(young_modulus,poisson_ratio),MU(young_modulus,poisson_ratio));
  ElasticFE_RVELagrange_Periodic MyFE_RVELagrangePeriodic(mField,Aij,D,F,applied_strain,"DISPLACEMENT","Lagrange_mul_disp",field_rank);
  ElasticFE_RVELagrange_RigidBodyTranslation MyFE_RVELagrangeRigidBodyTrans(mField,Aij,D,F,applied_strain,"DISPLACEMENT","Lagrange_mul_disp",field_rank,"Lagrange_mul_disp_rigid_trans");
  
  
  ierr = VecZeroEntries(F); CHKERRQ(ierr);
  ierr = VecGhostUpdateBegin(F,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = VecGhostUpdateEnd(F,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = MatZeroEntries(Aij); CHKERRQ(ierr);
  
  ierr = mField.loop_finite_elements("STOCHASIC_PROBLEM","ELASTIC",MyFE);  CHKERRQ(ierr);
  ierr = mField.loop_finite_elements("STOCHASIC_PROBLEM","Lagrange_elem",MyFE_RVELagrangePeriodic);  CHKERRQ(ierr);
  ierr = mField.loop_finite_elements("STOCHASIC_PROBLEM","Lagrange_elem_rigid_trans",MyFE_RVELagrangeRigidBodyTrans);  CHKERRQ(ierr);
  
  
  ierr = VecGhostUpdateBegin(F,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
  ierr = VecGhostUpdateEnd(F,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
  ierr = VecAssemblyBegin(F); CHKERRQ(ierr);
  ierr = VecAssemblyEnd(F); CHKERRQ(ierr);
  ierr = MatAssemblyBegin(Aij,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
  ierr = MatAssemblyEnd(Aij,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
  
  PetscSynchronizedFlush(PETSC_COMM_WORLD);
  //ierr = VecView(F,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);
  //ierr = MatView(Aij,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);
  
  
  //    ierr = VecView(F,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);
  //
  //
  //Matrix View
  //    ierr = MatView(Aij,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);
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
  
  //  ierr = VecView(F,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);
  //     ierr = VecView(D,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);
  //
  
  //Save data on mesh
  ierr = mField.set_global_VecCreateGhost("STOCHASIC_PROBLEM",ROW,D,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
  //ierr = VecView(F,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);

  //----------------------------------------------------------------------------
  // 3.2 Calculating zeroth-order homogenized stress using volume averaging theorem
  //----------------------------------------------------------------------------

  //create a vector for 6 components of homogenized stress
  Vec Stress_Homo;
  ierr = VecCreateMPI(PETSC_COMM_WORLD, 6, 6*pcomm->size(), &Stress_Homo);  CHKERRQ(ierr);
  ierr = VecZeroEntries(Stress_Homo); CHKERRQ(ierr);
  
  //    ierr = VecView(D,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);
  ElasticFE_RVELagrange_Homogenized_Stress_Periodic MyFE_RVEHomoStressPeriodic(mField,Aij,D,F,&RVE_volume,applied_strain,Stress_Homo,"DISPLACEMENT","Lagrange_mul_disp",field_rank);
  ierr = mField.loop_finite_elements("STOCHASIC_PROBLEM","Lagrange_elem",MyFE_RVEHomoStressPeriodic);  CHKERRQ(ierr);
  
//  if(pcomm->rank()) cout<< " Stress_Homo =  "<<endl;
//  ierr = VecView(Stress_Homo,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);
  
  if(pcomm->rank()==0){
    PetscScalar    *avec;
    VecGetArray(Stress_Homo, &avec);
    
    cout<< "\nStress_Homo = \n\n";
    for(int ii=0; ii<6; ii++){
      cout <<*avec<<endl; ;
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
  
  K_r_Young_PSFEMethod my_fe_k_r_E(mField,Aij,D,dF,LAMBDA(young_modulus,poisson_ratio),MU(young_modulus,poisson_ratio));
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
  ElasticFE_RVELagrange_Homogenized_Stress_Periodic MyFE_RVEHomoStressPeriodic_r_E(mField,Aij,D,F,&RVE_volume, applied_strain, Stress_Homo_r_E,"DISP_r_E","Lagrange_mul_disp",field_rank);
  ierr = mField.loop_finite_elements("STOCHASIC_PROBLEM","Lagrange_elem",MyFE_RVEHomoStressPeriodic_r_E);  CHKERRQ(ierr);
  
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
  
  K_r_Poisson_PSFEMethod my_fe_k_r_P(mField,Aij,D,dF,LAMBDA(young_modulus,poisson_ratio),MU(young_modulus,poisson_ratio));
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
  ElasticFE_RVELagrange_Homogenized_Stress_Periodic MyFE_RVEHomoStressPeriodic_r_P(mField,Aij,D,F,&RVE_volume, applied_strain, Stress_Homo_r_P,"DISP_r_P","Lagrange_mul_disp",field_rank);
  ierr = mField.loop_finite_elements("STOCHASIC_PROBLEM","Lagrange_elem",MyFE_RVEHomoStressPeriodic_r_P);  CHKERRQ(ierr);
  
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
  
  K_rs_Young_PSFEMethod my_fe_k_rs_EE(mField,Aij,D,ddF,LAMBDA(young_modulus,poisson_ratio),MU(young_modulus,poisson_ratio),"DISP_r_E");
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
//  ierr = VecView(ddD,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);

  //----------------------------------------------------------------------------
  // b. Calculating second-order homogenized stress
  //----------------------------------------------------------------------------
  
  Vec Stress_Homo_rs_EE;
  ierr = VecCreateMPI(PETSC_COMM_WORLD, 6, 6*pcomm->size(), &Stress_Homo_rs_EE);  CHKERRQ(ierr);
  ierr = VecZeroEntries(Stress_Homo_rs_EE); CHKERRQ(ierr);
  
  //    ierr = VecView(D,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);
  ElasticFE_RVELagrange_Homogenized_Stress_Periodic MyFE_RVEHomoStressPeriodic_rs_EE(mField,Aij,D,F,&RVE_volume, applied_strain, Stress_Homo_rs_EE,"DISP_rs_EE","Lagrange_mul_disp",field_rank);
  ierr = mField.loop_finite_elements("STOCHASIC_PROBLEM","Lagrange_elem",MyFE_RVEHomoStressPeriodic_rs_EE);  CHKERRQ(ierr);
  
  if(pcomm->rank()==0){
    PetscScalar    *avec_r_EE;
    VecGetArray(Stress_Homo_rs_EE, &avec_r_EE);
    
    cout<< "\nStress_Homo_rs_EE = \n\n";
    for(int ii=0; ii<6; ii++){
      cout <<*avec_r_EE<<endl; ;
      avec_r_EE++;
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
  
  K_rs_EP_PSFEMethod my_fe_k_rs_EP(mField,Aij,D,ddF,LAMBDA(young_modulus,poisson_ratio),MU(young_modulus,poisson_ratio),"DISP_r_P");
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
  ElasticFE_RVELagrange_Homogenized_Stress_Periodic MyFE_RVEHomoStressPeriodic_rs_EP(mField,Aij,D,F,&RVE_volume, applied_strain, Stress_Homo_rs_EP,"DISP_rs_EP","Lagrange_mul_disp",field_rank);
  ierr = mField.loop_finite_elements("STOCHASIC_PROBLEM","Lagrange_elem",MyFE_RVEHomoStressPeriodic_rs_EP);  CHKERRQ(ierr);
  
  if(pcomm->rank()==0){
    PetscScalar    *avec_r_EP;
    VecGetArray(Stress_Homo_rs_EP, &avec_r_EP);
    
    cout<< "\nStress_Homo_rs_EP = \n\n";
    for(int ii=0; ii<6; ii++){
      cout <<*avec_r_EP<<endl; ;
      avec_r_EP++;
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
  
  K_rs_Poisson_PSFEMethod my_fe_k_rs_PP(mField,Aij,D,ddF,LAMBDA(young_modulus,poisson_ratio),MU(young_modulus,poisson_ratio),"DISP_r_P");
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
//  ierr = VecView(ddD_PP,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);

  //----------------------------------------------------------------------------
  // b. Calculating second-order homogenized stress
  //----------------------------------------------------------------------------
  
  Vec Stress_Homo_rs_PP;
  ierr = VecCreateMPI(PETSC_COMM_WORLD, 6, 6*pcomm->size(), &Stress_Homo_rs_PP);  CHKERRQ(ierr);
  ierr = VecZeroEntries(Stress_Homo_rs_PP); CHKERRQ(ierr);
  
  //    ierr = VecView(D,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);
  ElasticFE_RVELagrange_Homogenized_Stress_Periodic MyFE_RVEHomoStressPeriodic_rs_PP(mField,Aij,D,F,&RVE_volume, applied_strain, Stress_Homo_rs_PP,"DISP_rs_PP","Lagrange_mul_disp",field_rank);
  ierr = mField.loop_finite_elements("STOCHASIC_PROBLEM","Lagrange_elem",MyFE_RVEHomoStressPeriodic_rs_PP);  CHKERRQ(ierr);
  
  if(pcomm->rank()==0){
    PetscScalar    *avec_r_PP;
    VecGetArray(Stress_Homo_rs_PP, &avec_r_PP);
    
    cout<< "\nStress_Homo_rs_PP = \n\n";
    for(int ii=0; ii<6; ii++){
      cout <<*avec_r_PP<<endl; ;
      avec_r_PP++;
    }
  }
  cout<< "\n\n";  
  
  
  // ===========================================================================
  //
  //  VII. FINISH
  //
  // ===========================================================================   
  
  
  //Destroy matrices
  ierr = VecDestroy(&F); CHKERRQ(ierr);
  ierr = VecDestroy(&D); CHKERRQ(ierr);
  ierr = MatDestroy(&Aij); CHKERRQ(ierr);
  ierr = KSPDestroy(&solver); CHKERRQ(ierr);
  
  
  ierr = PetscTime(&v2);CHKERRQ(ierr);
  ierr = PetscGetCPUTime(&t2);CHKERRQ(ierr);
  
  PetscSynchronizedPrintf(PETSC_COMM_WORLD,"Total Rank %d Time = %f CPU Time = %f\n",pcomm->rank(),v2-v1,t2-t1);
  PetscSynchronizedFlush(PETSC_COMM_WORLD);
  
  PetscFinalize();
  
}

