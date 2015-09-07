/* Copyright (C) 2014, Zahur Ullah (Zahur.Ullah AT glasgow.ac.uk)
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

#include <Projection10NodeCoordsOnField.hpp>
#include <petsctime.h>
#include <FEMethod_LowLevelStudent.hpp>
#include <FEMethod_UpLevelStudent.hpp>
#include <PostProcVertexMethod.hpp>
#include <PostProcDisplacementAndStrainOnRefindedMesh.hpp>

#include <ThermalElement.hpp>
#include <MoistureTransportElement.hpp>
#include "ElasticFE_RVELagrange_Periodic.hpp"
#include "ElasticFE_RVELagrange_Periodic_Multi_Rhs.hpp"
#include "ElasticFE_RVELagrange_RigidBodyTranslation.hpp"
#include "RVEVolume.hpp"
#include "ElasticFE_RVELagrange_Homogenized_Stress_Periodic.hpp"


ErrorCode rval;
PetscErrorCode ierr;

static char help[] = "...\n\n";

//=================================================================================================================================
//Define class and multindex container to store data for traiangles on the boundary of the RVE (it cannot be defined within main)
//=================================================================================================================================

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

//============================================================================================
//============================================================================================

int main(int argc, char *argv[]) {
  
  PetscInitialize(&argc,&argv,(char *)0,help);
  
  moab::Core mb_instance;
  Interface& moab = mb_instance;
  int rank;
  MPI_Comm_rank(PETSC_COMM_WORLD,&rank);
  
  //Reade parameters from line command
  PetscBool flg = PETSC_TRUE;
  char mesh_file_name[255];
  ierr = PetscOptionsGetString(PETSC_NULL,"-my_file",mesh_file_name,255,&flg); CHKERRQ(ierr);
  if(flg != PETSC_TRUE) {
    SETERRQ(PETSC_COMM_SELF,1,"*** ERROR -my_file (MESH FILE NEEDED)");
  }
  
  char output_file_Dmat[255];
  ierr = PetscOptionsGetString(PETSC_NULL,"-my_output_file_Dmat",output_file_Dmat,255,&flg); CHKERRQ(ierr);
  if(flg != PETSC_TRUE) {
    SETERRQ(PETSC_COMM_SELF,1,"*** ERROR -my_file (NAME OF DMAT OUTPUT FILE IS NEEDED)");
  }

  
  PetscInt order;
  ierr = PetscOptionsGetInt(PETSC_NULL,"-my_order",&order,&flg); CHKERRQ(ierr);
  if(flg != PETSC_TRUE) {
    order = 1;
  }
  
  //Read mesh to MOAB
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
  
  //ref meshset ref level 0
  ierr = mField.seed_ref_level_3D(0,0); CHKERRQ(ierr);
  
  // stl::bitset see for more details
  BitRefLevel bit_level0;
  bit_level0.set(0);
  EntityHandle meshset_level0;
  rval = moab.create_meshset(MESHSET_SET,meshset_level0); CHKERR_PETSC(rval);
  ierr = mField.seed_ref_level_3D(0,bit_level0); CHKERRQ(ierr);
  ierr = mField.get_entities_by_ref_level(bit_level0,BitRefLevel().set(),meshset_level0); CHKERRQ(ierr);
  
  /***/
  //Define problem
  
  //Field
  int field_rank=1;
  ierr = mField.add_field("CONC",H1,field_rank); CHKERRQ(ierr);
  ierr = mField.add_field("LAGRANGE_MUL_FIELD",H1,field_rank); CHKERRQ(ierr);
  ierr = mField.add_field("LAGRANGE_MUL_FIELD_rigid_trans",NOFIELD,field_rank); CHKERRQ(ierr);
  
  //Problem
  ierr = mField.add_problem("MOISTURE_PROBLEM"); CHKERRQ(ierr);
  

  //meshset consisting all entities in mesh
  EntityHandle root_set = moab.get_root_set();
  //add entities to field
  ierr = mField.add_ents_to_field_by_TETs(root_set,"CONC"); CHKERRQ(ierr);

  //FE
  MoistureTransportElement moisture_elements(mField);
  ierr = moisture_elements.addDiffusionElement("MOISTURE_PROBLEM","CONC"); CHKERRQ(ierr);
  ierr = mField.add_finite_element("LAGRANGE_FE"); CHKERRQ(ierr);
  ierr = mField.add_finite_element("LAGRANGE_FE_rigid_trans"); CHKERRQ(ierr);//For rigid body control

  //Define Rows, cols and data for matrix C and CT for
  //============================================================================================================
  //C row as LAGRANGE_MUL_FIELD and col as CONC
  ierr = mField.modify_finite_element_add_field_row("LAGRANGE_FE","LAGRANGE_MUL_FIELD"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_col("LAGRANGE_FE","CONC"); CHKERRQ(ierr);
  
  //CT col as Lagrange_mul_disp and row as DISPLACEMENT
  ierr = mField.modify_finite_element_add_field_col("LAGRANGE_FE","LAGRANGE_MUL_FIELD"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_row("LAGRANGE_FE","CONC"); CHKERRQ(ierr);
  
  //Data
  ierr = mField.modify_finite_element_add_field_data("LAGRANGE_FE","LAGRANGE_MUL_FIELD"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_data("LAGRANGE_FE","CONC"); CHKERRQ(ierr);
  //============================================================================================================

  //Define rows/cols and element data for C1 and C1T (for lagrange multipliers to contol the rigid body motion)
  //============================================================================================================
  //C1 row as Lagrange_mul_disp and col as DISPLACEMENT
  ierr = mField.modify_finite_element_add_field_row("LAGRANGE_FE_rigid_trans","LAGRANGE_MUL_FIELD_rigid_trans"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_col("LAGRANGE_FE_rigid_trans","CONC"); CHKERRQ(ierr);
  
  //C1T col as Lagrange_mul_disp and row as DISPLACEMENT
  ierr = mField.modify_finite_element_add_field_col("LAGRANGE_FE_rigid_trans","LAGRANGE_MUL_FIELD_rigid_trans"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_row("LAGRANGE_FE_rigid_trans","CONC"); CHKERRQ(ierr);
  
  //Data
  ierr = mField.modify_finite_element_add_field_data("LAGRANGE_FE_rigid_trans","LAGRANGE_MUL_FIELD_rigid_trans"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_data("LAGRANGE_FE_rigid_trans","CONC"); CHKERRQ(ierr);
  //============================================================================================================

  
  //set finite elements for problem
  ierr = mField.modify_problem_add_finite_element("MOISTURE_PROBLEM","LAGRANGE_FE"); CHKERRQ(ierr);
  ierr = mField.modify_problem_add_finite_element("MOISTURE_PROBLEM","LAGRANGE_FE_rigid_trans"); CHKERRQ(ierr);
  
  //set refinment level for problem
  ierr = mField.modify_problem_ref_level_add_bit("MOISTURE_PROBLEM",bit_level0); CHKERRQ(ierr);
  
  /***/
  //Declare problem
  
  //add entitities (by tets) to the field
  ierr = mField.add_ents_to_field_by_TETs(0,"CONC"); CHKERRQ(ierr);
  
  Range SurfacesFaces;
  ierr = mField.get_cubit_msId_entities_by_dimension(103,SIDESET,2,SurfacesFaces,true); CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"number of SideSet 103 = %d\n",SurfacesFaces.size()); CHKERRQ(ierr);
  
  //to create meshset from range
  EntityHandle BoundFacesMeshset;
  rval = moab.create_meshset(MESHSET_SET,BoundFacesMeshset); CHKERR_PETSC(rval);
	rval = moab.add_entities(BoundFacesMeshset,SurfacesFaces); CHKERR_PETSC(rval);
  ierr = mField.seed_ref_level_MESHSET(BoundFacesMeshset,BitRefLevel().set()); CHKERRQ(ierr);
  ierr = mField.add_ents_to_finite_element_by_TRIs(SurfacesFaces,"LAGRANGE_FE_rigid_trans"); CHKERRQ(ierr);

  //=======================================================================================================
  //Add Periodic Prisims Between Triangles on -ve and +ve faces to implement periodic bounary conditions
  //=======================================================================================================
  
  //Populating the Multi-index container with -ve triangles
  Range SurTrisNeg;
  ierr = mField.get_cubit_msId_entities_by_dimension(101,SIDESET,2,SurTrisNeg,true); CHKERRQ(ierr);
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
  ierr = mField.get_cubit_msId_entities_by_dimension(102,SIDESET,2,SurTrisPos,true); CHKERRQ(ierr);
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
  //  ierr = mField.get_entities_by_ref_level(bit_level,BitRefLevel().set(),meshset_level0_RVE); CHKERRQ(ierr);
  //ierr = mField.add_ents_to_finite_element_EntType_by_bit_ref(problem_bit_level_RVE,"Lagrange_elem",MBPRISM); CHKERRQ(ierr);
  ierr = mField.add_ents_to_finite_element_by_PRISMs(PrismRange, "LAGRANGE_FE"); CHKERRQ(ierr);
  
  //Adding only -ve surfaces to the field Lagrange_mul_disp (in periodic boundary conditions size of C (3M/2 x 3N))
  //to create meshset from range  SurTrisNeg
  EntityHandle SurTrisNegMeshset;
  rval = moab.create_meshset(MESHSET_SET,SurTrisNegMeshset); CHKERR_PETSC(rval);
	rval = moab.add_entities(SurTrisNegMeshset,SurTrisNeg); CHKERR_PETSC(rval);
  ierr = mField.add_ents_to_field_by_TRIs(SurTrisNegMeshset,"LAGRANGE_MUL_FIELD",2); CHKERRQ(ierr);
  //=======================================================================================================
  //=======================================================================================================
  

  //set app. order
  //see Hierarchic Finite Element Bases on Unstructured Tetrahedral Meshes (Mark Ainsworth & Joe Coyle)
  //int order = 5;
  ierr = mField.set_field_order(root_set,MBTET,"CONC",order); CHKERRQ(ierr);
  ierr = mField.set_field_order(root_set,MBTRI,"CONC",order); CHKERRQ(ierr);
  ierr = mField.set_field_order(root_set,MBEDGE,"CONC",order); CHKERRQ(ierr);
  ierr = mField.set_field_order(root_set,MBVERTEX,"CONC",1); CHKERRQ(ierr);

  ierr = mField.set_field_order(0,MBTRI,"LAGRANGE_MUL_FIELD",order); CHKERRQ(ierr);
  ierr = mField.set_field_order(0,MBEDGE,"LAGRANGE_MUL_FIELD",order); CHKERRQ(ierr);
  ierr = mField.set_field_order(0,MBVERTEX,"LAGRANGE_MUL_FIELD",1); CHKERRQ(ierr);
  
  
  
  //create these elements to calculate separte volumes for each fibres and matrix.
  //======================================================================================
  //FE
  ierr = mField.add_finite_element("FE_MATRXIX"); CHKERRQ(ierr);
  ierr = mField.add_finite_element("FE_FIBRES"); CHKERRQ(ierr);
  
  //Row and columns
  ierr = mField.modify_finite_element_add_field_row("FE_MATRXIX","CONC"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_col("FE_MATRXIX","CONC"); CHKERRQ(ierr);
  
  ierr = mField.modify_finite_element_add_field_row("FE_FIBRES","CONC"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_col("FE_FIBRES","CONC"); CHKERRQ(ierr);
  
  ierr = mField.modify_finite_element_add_field_data("FE_MATRXIX","CONC"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_data("FE_FIBRES","CONC"); CHKERRQ(ierr);

  //set finite elements for problem
  ierr = mField.modify_problem_add_finite_element("MOISTURE_PROBLEM","FE_MATRXIX"); CHKERRQ(ierr);
  ierr = mField.modify_problem_add_finite_element("MOISTURE_PROBLEM","FE_FIBRES"); CHKERRQ(ierr);
  
  //add finite elements entities
  EntityHandle meshset_Elastic, meshset_Elastic_stiff_inclusion;
  rval = moab.create_meshset(MESHSET_SET,meshset_Elastic); CHKERR_PETSC(rval);
  rval = moab.create_meshset(MESHSET_SET,meshset_Elastic_stiff_inclusion); CHKERR_PETSC(rval);
  
  
  double density_matrix, density_fibres;
  Range TetsInBlock, TetsInBlock_stiff_inc;
  for(_IT_CUBITMESHSETS_BY_SET_TYPE_FOR_LOOP_(mField,BLOCKSET,it)){
		if(it->get_name() == "MAT_MOISTURE_MATRIX") {
			rval = moab.get_entities_by_type(it->meshset, MBTET,TetsInBlock,true); CHKERR_PETSC(rval);
      
      Mat_Moisture mydata;
      ierr = it->get_attribute_data_structure(mydata); CHKERRQ(ierr);
      density_matrix=mydata.data.Density;

		}
    if(it->get_name() == "MAT_MOISTURE_FIBRES") {
			rval = moab.get_entities_by_type(it->meshset, MBTET,TetsInBlock_stiff_inc,true); CHKERR_PETSC(rval);
      
      Mat_Moisture mydata;
      ierr = it->get_attribute_data_structure(mydata); CHKERRQ(ierr);
      density_fibres=mydata.data.Density;

		}
	}
  cout<<"TetsInBlock "<<TetsInBlock<<endl;
  rval = moab.add_entities(meshset_Elastic,TetsInBlock);CHKERR_PETSC(rval);
  
  cout<<"TetsInBlock "<<TetsInBlock_stiff_inc<<endl;
  rval = moab.add_entities(meshset_Elastic_stiff_inclusion,TetsInBlock_stiff_inc);CHKERR_PETSC(rval);
  
  ierr = mField.add_ents_to_finite_element_by_TETs(meshset_Elastic,"FE_MATRXIX",true); CHKERRQ(ierr);
  ierr = mField.add_ents_to_finite_element_by_TETs(meshset_Elastic_stiff_inclusion,"FE_FIBRES",true); CHKERRQ(ierr);
  
  
  //======================================================================================

  
  /****/
  //build database
  
  //build field
  ierr = mField.build_fields(); CHKERRQ(ierr);
  
  //build finite elemnts
  ierr = mField.build_finite_elements(); CHKERRQ(ierr);
  
  //build adjacencies
  ierr = mField.build_adjacencies(bit_level0); CHKERRQ(ierr);
  
  //build problem
  ierr = mField.build_problems(); CHKERRQ(ierr);
  
  
  /****/
  //mesh partitioning

  //partition
  ierr = mField.partition_problem("MOISTURE_PROBLEM"); CHKERRQ(ierr);
  ierr = mField.partition_finite_elements("MOISTURE_PROBLEM"); CHKERRQ(ierr);
  //what are ghost nodes, see Petsc Manual
  ierr = mField.partition_ghost_dofs("MOISTURE_PROBLEM"); CHKERRQ(ierr);
  
  //print block sets with materials
  ierr = mField.print_cubit_materials_set(); CHKERRQ(ierr);
  
  //create matrices (here F, D and Aij are matrices for the full problem)
  Vec F1,F2,F3,C1;
  ierr = mField.VecCreateGhost("MOISTURE_PROBLEM",ROW,&F1); CHKERRQ(ierr);
  ierr = mField.VecCreateGhost("MOISTURE_PROBLEM",ROW,&F2); CHKERRQ(ierr);
  ierr = mField.VecCreateGhost("MOISTURE_PROBLEM",ROW,&F3); CHKERRQ(ierr);
  
  ierr = mField.VecCreateGhost("MOISTURE_PROBLEM",COL,&C1); CHKERRQ(ierr);
  Mat A;
  ierr = mField.MatCreateMPIAIJWithArrays("MOISTURE_PROBLEM",&A); CHKERRQ(ierr);

  
  ublas::vector<FieldData> applied_strain;  //it is not used in the calculation, it is required by ElasticFE_RVELagrange_Disp as input
  applied_strain.resize(1.5*field_rank+1.5); applied_strain.clear();
  
  ierr = moisture_elements.setThermalFiniteElementLhsOperators("CONC",A); CHKERRQ(ierr);
  ElasticFE_RVELagrange_Periodic_Multi_Rhs MyFE_RVELagrange(mField,A,C1,F1,F2,F3,applied_strain,"CONC","LAGRANGE_MUL_FIELD",field_rank);
  ElasticFE_RVELagrange_RigidBodyTranslation MyFE_RVELagrangeRigidBodyTrans(mField,A,C1,F1,applied_strain,"CONC","LAGRANGE_MUL_FIELD",field_rank,"LAGRANGE_MUL_FIELD_rigid_trans");

  ierr = VecZeroEntries(F1); CHKERRQ(ierr);
  ierr = VecGhostUpdateBegin(F1,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = VecGhostUpdateEnd(F1,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  
  ierr = VecZeroEntries(F2); CHKERRQ(ierr);
  ierr = VecGhostUpdateBegin(F2,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = VecGhostUpdateEnd(F2,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  
  ierr = VecZeroEntries(F1); CHKERRQ(ierr);
  ierr = VecGhostUpdateBegin(F3,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = VecGhostUpdateEnd(F3,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);

  ierr = MatZeroEntries(A); CHKERRQ(ierr);
  
  ierr = mField.loop_finite_elements("MOISTURE_PROBLEM","DIFFUSION_FE",moisture_elements.getLoopFeLhs()); CHKERRQ(ierr);
  ierr = mField.loop_finite_elements("MOISTURE_PROBLEM","LAGRANGE_FE",MyFE_RVELagrange);  CHKERRQ(ierr);
  ierr = mField.loop_finite_elements("MOISTURE_PROBLEM","LAGRANGE_FE_rigid_trans",MyFE_RVELagrangeRigidBodyTrans);  CHKERRQ(ierr);

  ierr = VecGhostUpdateBegin(F1,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
  ierr = VecGhostUpdateEnd(F1,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
  ierr = VecAssemblyBegin(F1); CHKERRQ(ierr);
  ierr = VecAssemblyEnd(F1); CHKERRQ(ierr);

  ierr = VecGhostUpdateBegin(F2,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
  ierr = VecGhostUpdateEnd(F2,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
  ierr = VecAssemblyBegin(F2); CHKERRQ(ierr);
  ierr = VecAssemblyEnd(F2); CHKERRQ(ierr);

  ierr = VecGhostUpdateBegin(F3,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
  ierr = VecGhostUpdateEnd(F3,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
  ierr = VecAssemblyBegin(F3); CHKERRQ(ierr);
  ierr = VecAssemblyEnd(F3); CHKERRQ(ierr);

  ierr = MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
  ierr = MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);

  //ierr = VecView(F,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);
  //ierr = MatView(A,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);
  
//  //Matrix View
//  MatView(A,PETSC_VIEWER_DRAW_WORLD);//PETSC_VIEWER_STDOUT_WORLD);
//  std::string wait;
//  std::cin >> wait;

  
//Calculation of Homogenized stress
//=========================================================================================================================
  const double young_modulus = 1; //dummy values
  const double poisson_ratio = 0.0;
  
  double RVE_volume;    RVE_volume=0.0;  //RVE volume for full RVE We need this for stress calculation
  Vec RVE_volume_Vec;
  ierr = VecCreateMPI(PETSC_COMM_WORLD, 1, pcomm->size(), &RVE_volume_Vec);  CHKERRQ(ierr);
  ierr = VecZeroEntries(RVE_volume_Vec); CHKERRQ(ierr);
  
  RVEVolume MyRVEVol(mField,A,C1,F1,LAMBDA(young_modulus,poisson_ratio),MU(young_modulus,poisson_ratio), RVE_volume_Vec);
//  ierr = mField.loop_finite_elements("MOISTURE_PROBLEM","DIFFUSION_FE",MyRVEVol);  CHKERRQ(ierr);
  
  //  ierr = VecView(RVE_volume_Vec,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);
//  ierr = VecSum(RVE_volume_Vec, &RVE_volume);  CHKERRQ(ierr);
//  cout<<"Final RVE_volume = "<< RVE_volume <<endl;
  
  double volume_matrix, volume_fibres;
  ierr = mField.loop_finite_elements("MOISTURE_PROBLEM","FE_MATRXIX",MyRVEVol);  CHKERRQ(ierr);
  ierr = VecSum(RVE_volume_Vec, &volume_matrix);  CHKERRQ(ierr);
  cout<<"Matrix volume = "<< volume_matrix <<endl;

  ierr = mField.loop_finite_elements("MOISTURE_PROBLEM","FE_FIBRES",MyRVEVol);  CHKERRQ(ierr);
  ierr = VecSum(RVE_volume_Vec, &RVE_volume);  CHKERRQ(ierr);
  cout<<"Final RVE_volume = "<< RVE_volume <<endl;
  
  volume_fibres=RVE_volume-volume_matrix;
  cout<<"Fibres volume = "<< volume_fibres <<endl;
  
  cout<<"density_matrix = "<< density_matrix <<endl;
  cout<<"density_fibres = "<< density_fibres <<endl;
  double RVE_density;
  RVE_density=(volume_matrix/RVE_volume)*density_matrix +  (volume_fibres/RVE_volume)*density_fibres;
  
////=========================================================================================================================

  //Solver
  KSP solver;
  ierr = KSPCreate(PETSC_COMM_WORLD,&solver); CHKERRQ(ierr);
  ierr = KSPSetOperators(solver,A,A); CHKERRQ(ierr);
  ierr = KSPSetFromOptions(solver); CHKERRQ(ierr);
  ierr = KSPSetUp(solver); CHKERRQ(ierr);
  
  ierr = KSPSolve(solver,F1,C1); CHKERRQ(ierr);
  ierr = VecGhostUpdateBegin(C1,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = VecGhostUpdateEnd(C1,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = mField.set_global_ghost_vector("MOISTURE_PROBLEM",ROW,C1,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
  
  ublas::matrix<FieldData> Dmat;
  Dmat.resize(3,3); Dmat.clear();
  //=============================================================================================================
  // homogenised stress for strian [1 0 0]^T
  //=============================================================================================================
  //create a vector for 3 components of homogenized stress
  Vec Stress_Homo;
  if(pcomm->rank()==0) {
    VecCreateGhost(PETSC_COMM_WORLD,3,3,0,PETSC_NULL,&Stress_Homo);
  } else {
    int ghost[] = {0,1,2};
    VecCreateGhost(PETSC_COMM_WORLD,0,3,3,ghost,&Stress_Homo);
    
  }

  {
    ierr = VecZeroEntries(Stress_Homo); CHKERRQ(ierr);
    ElasticFE_RVELagrange_Homogenized_Stress_Periodic MyFE_RVEHomoStressPeriodic(mField,A,C1,F1,&RVE_volume, applied_strain, Stress_Homo,"CONC","LAGRANGE_MUL_FIELD",field_rank);
    ierr = mField.loop_finite_elements("MOISTURE_PROBLEM","LAGRANGE_FE",MyFE_RVEHomoStressPeriodic);  CHKERRQ(ierr);
    VecGhostUpdateBegin(Stress_Homo,INSERT_VALUES,SCATTER_FORWARD);
    VecGhostUpdateEnd(Stress_Homo,INSERT_VALUES,SCATTER_FORWARD);
    PetscScalar *avec;
    VecGetArray(Stress_Homo, &avec);
    for(int ii=0; ii<3; ii++){
      Dmat(ii,0)=*avec;
      avec++;
    }
    
    if(pcomm->rank()==0){
      cout<< "\nStress_Homo = \n\n";
      for(int ii=0; ii<3; ii++){
        cout <<Dmat(ii,0)<<endl;
      }
    }
  }

  //=============================================================================================================
  // homogenised stress for strian [0 1 0]^T
  //=============================================================================================================
  ierr = KSPSolve(solver,F2,C1); CHKERRQ(ierr);
  ierr = VecGhostUpdateBegin(C1,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = VecGhostUpdateEnd(C1,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = mField.set_global_ghost_vector("MOISTURE_PROBLEM",ROW,C1,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
  {
    ierr = VecZeroEntries(Stress_Homo); CHKERRQ(ierr);
    ElasticFE_RVELagrange_Homogenized_Stress_Periodic MyFE_RVEHomoStressPeriodic(mField,A,C1,F2,&RVE_volume, applied_strain, Stress_Homo,"CONC","LAGRANGE_MUL_FIELD",field_rank);
    ierr = mField.loop_finite_elements("MOISTURE_PROBLEM","LAGRANGE_FE",MyFE_RVEHomoStressPeriodic);  CHKERRQ(ierr);
    VecGhostUpdateBegin(Stress_Homo,INSERT_VALUES,SCATTER_FORWARD);
    VecGhostUpdateEnd(Stress_Homo,INSERT_VALUES,SCATTER_FORWARD);
    PetscScalar *avec;
    VecGetArray(Stress_Homo, &avec);
    for(int ii=0; ii<3; ii++){
      Dmat(ii,1)=*avec;
      avec++;
    }
    
    if(pcomm->rank()==0){
      cout<< "\nStress_Homo = \n\n";
      for(int ii=0; ii<3; ii++){
        cout <<Dmat(ii,1)<<endl;
      }
    }
  }

  
  //=============================================================================================================
  // homogenised stress for strian [0 0 1]^T
  //=============================================================================================================
  ierr = KSPSolve(solver,F3,C1); CHKERRQ(ierr);
  ierr = VecGhostUpdateBegin(C1,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = VecGhostUpdateEnd(C1,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = mField.set_global_ghost_vector("MOISTURE_PROBLEM",ROW,C1,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);

  {
    ierr = VecZeroEntries(Stress_Homo); CHKERRQ(ierr);
    ElasticFE_RVELagrange_Homogenized_Stress_Periodic MyFE_RVEHomoStressPeriodic(mField,A,C1,F3,&RVE_volume, applied_strain, Stress_Homo,"CONC","LAGRANGE_MUL_FIELD",field_rank);
    ierr = mField.loop_finite_elements("MOISTURE_PROBLEM","LAGRANGE_FE",MyFE_RVEHomoStressPeriodic);  CHKERRQ(ierr);
    VecGhostUpdateBegin(Stress_Homo,INSERT_VALUES,SCATTER_FORWARD);
    VecGhostUpdateEnd(Stress_Homo,INSERT_VALUES,SCATTER_FORWARD);
    PetscScalar *avec;
    VecGetArray(Stress_Homo, &avec);
    for(int ii=0; ii<3; ii++){
      Dmat(ii,2)=*avec;
      avec++;
    }
    
    if(pcomm->rank()==0){
      cout<< "\nStress_Homo = \n\n";
      for(int ii=0; ii<3; ii++){
        cout <<Dmat(ii,2)<<endl;
      }
    }
    
  }
  //====================================================================================================
  //Writing Dmat in to binary file for use in the unsteady diffusion problem
  cout <<"Conductivity matrix="<< Dmat<< endl;
  Dmat=Dmat/RVE_density; //Here we will save Dmat/RVE_density [mm2/s]
  
  cout <<"Effective desnsity RVE ="<< RVE_density<< endl;
  cout <<"Diffusivity matix ="<< Dmat<< endl;
  
  if(pcomm->rank()==0){
    int fd;
    PetscViewer view_out;
    PetscViewerBinaryOpen(PETSC_COMM_WORLD,output_file_Dmat,FILE_MODE_WRITE,&view_out);
    PetscViewerBinaryGetDescriptor(view_out,&fd);
    PetscBinaryWrite(fd,&Dmat(0,0),9,PETSC_DOUBLE,PETSC_FALSE);
    PetscViewerDestroy(&view_out);
  }
  //====================================================================================================


  ProjectionFieldOn10NodeTet ent_method_on_10nodeTet(mField,"CONC",true,false,"CONC");
  ierr = mField.loop_dofs("CONC",ent_method_on_10nodeTet); CHKERRQ(ierr);
  ent_method_on_10nodeTet.set_nodes = false;
  ierr = mField.loop_dofs("CONC",ent_method_on_10nodeTet); CHKERRQ(ierr);
  
  if(pcomm->rank()==0) {
    EntityHandle out_meshset;
    rval = moab.create_meshset(MESHSET_SET,out_meshset); CHKERR_PETSC(rval);
    ierr = mField.get_problem_finite_elements_entities("MOISTURE_PROBLEM","DIFFUSION_FE",out_meshset); CHKERRQ(ierr);
    rval = moab.write_file("out.vtk","VTK","",&out_meshset,1); CHKERR_PETSC(rval);
    rval = moab.delete_entities(&out_meshset,1); CHKERR_PETSC(rval);
  }
  
   //Destroy matrices
  ierr = VecDestroy(&F1); CHKERRQ(ierr);
  ierr = VecDestroy(&F2); CHKERRQ(ierr);
  ierr = VecDestroy(&F2); CHKERRQ(ierr);

  ierr = VecDestroy(&C1); CHKERRQ(ierr);

  ierr = MatDestroy(&A); CHKERRQ(ierr);
  ierr = KSPDestroy(&solver); CHKERRQ(ierr);
  
  
  ierr = PetscTime(&v2);CHKERRQ(ierr);
  ierr = PetscGetCPUTime(&t2);CHKERRQ(ierr);
  
  PetscSynchronizedPrintf(PETSC_COMM_WORLD,"Total Rank %d Time = %f CPU Time = %f\n",pcomm->rank(),v2-v1,t2-t1);
  
  PetscFinalize();
  
}

