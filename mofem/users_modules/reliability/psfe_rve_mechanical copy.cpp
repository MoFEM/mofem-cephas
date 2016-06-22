/** \file rve_mechanical.cpp
 * \brief Calculates stiffness matrix for elastic RVE.

 Three types of boundary conditions are implemented, i.e.
 HOMOBCDISP, HOMOBCPERIODIC, HOMOBCTRAC, NITSCHE.

 NITSHCE method allow to apply periodic boundary conditions
 to arbitrary convex shape RVE.

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

#include <boost/numeric/ublas/vector_proxy.hpp>

#include <MethodForForceScaling.hpp>
#include <TimeForceScale.hpp>

#include <BCs_RVELagrange_Disp.hpp>
#include <BCs_RVELagrange_Trac.hpp>
#include <BCs_RVELagrange_Trac_Rigid_Rot.hpp>
#include <BCs_RVELagrange_Trac_Rigid_Trans.hpp>
#include <BCs_RVELagrange_Periodic.hpp>

#ifndef WITH_ADOL_C
  #error "MoFEM need to be compiled with ADOL-C"
#endif

#include <adolc/adolc.h>
#include <NonLinearElasticElement.hpp>
#include <Hooke.hpp>
#include <boost/numeric/ublas/symmetric.hpp>
#include <SmallTransverselyIsotropic.hpp>
#include <VolumeCalculation.hpp>

#include <boost/ptr_container/ptr_map.hpp>
#include <PostProcOnRefMesh.hpp>
#include <PostProcStresses.hpp>

extern "C" {
  void tetcircumcenter_tp(double a[3],double b[3],double c[3], double d[3],
    double circumcenter[3],double *xi,double *eta,double *zeta);
  void tricircumcenter3d_tp(double a[3],double b[3],double c[3],
    double circumcenter[3],double *xi,double *eta);
  //#include <triangle_ncc_rule.h>
}

#include <NitscheMethod.hpp>
#include <moab/AdaptiveKDTree.hpp>
#include <NitschePeriodicMethod.hpp>

#include <algorithm>

ErrorCode rval;
PetscErrorCode ierr;

static char help[] = "...\n\n";

enum HomoBCTypes {
  HOMOBCDISP,
  HOMOBCPERIODIC,
  HOMOBCTRAC,
  NITSCHE,
  NOHOMOBC
};

const char *homo_bc_names[] = {
  "disp",
  "periodic",
  "trac",
  "nitsche",
  "nohomobc"
};

//=================================================================================================================================
//Define class and multindex container to store data for traiangles on the boundary of the RVE (it cannot be defined within main)
//=================================================================================================================================

struct Face_CenPos_Handle {
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
    SETERRQ(PETSC_COMM_SELF,MOFEM_NOT_FOUND,"*** ERROR -my_file (MESH FILE NEEDED)");
  }

  PetscInt order;
  ierr = PetscOptionsGetInt(PETSC_NULL,"-my_order",&order,&flg); CHKERRQ(ierr);
  if(flg != PETSC_TRUE) {
    order = 1;
  }

  PetscInt choise_value = NOHOMOBC;
  ierr = PetscOptionsGetEList(
    NULL,"-my_bc_type",homo_bc_names,NOHOMOBC,&choise_value,&flg
  ); CHKERRQ(ierr);
  if(flg != PETSC_TRUE) {
    SETERRQ(PETSC_COMM_SELF,MOFEM_IMPOSIBLE_CASE,"*** Boundary conditions not set");
  }

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
  FieldInterface& m_field = core;

  vector<BitRefLevel> bit_levels;
  {
    Tag th_meshset_info;
    int def_meshset_info[2] = {0,0};
    rval = moab.tag_get_handle(
      "MESHSET_INFO",2,MB_TYPE_INTEGER,th_meshset_info,MB_TAG_CREAT|MB_TAG_SPARSE,&def_meshset_info
    );
    int meshset_data[2];
    EntityHandle root = moab.get_root_set();
    rval = moab.tag_get_data(th_meshset_info,&root,1,meshset_data); CHKERR_PETSC(rval);
    if(meshset_data[0]==0) {
      meshset_data[0] = 1;
      rval = moab.tag_set_data(th_meshset_info,&root,1,meshset_data); CHKERR_PETSC(rval);

    }
    bit_levels.push_back(BitRefLevel().set(meshset_data[0]-1));
  }
  BitRefLevel problem_bit_level = bit_levels.back();
  ierr = m_field.seed_ref_level_3D(0,problem_bit_level); CHKERRQ(ierr);

  //    const clock_t begin_time = clock();
  ierr = m_field.build_fields(); CHKERRQ(ierr);
  ierr = m_field.build_finite_elements(); CHKERRQ(ierr);

  Range preriodic_prisms;
  if(choise_value == HOMOBCPERIODIC) {
    //FIXME: Naming convention is not consistent in this section of code
    //=======================================================================================================
    //Add Periodic Prisims Between Triangles on -ve and +ve faces to implement periodic bounary conditions
    //=======================================================================================================
    //Populating the Multi-index container with -ve triangles
    Range SurTrisNeg;
    ierr = m_field.get_cubit_msId_entities_by_dimension(101,SIDESET,2,SurTrisNeg,true); CHKERRQ(ierr);
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
    ierr = m_field.get_cubit_msId_entities_by_dimension(102,SIDESET,2,SurTrisPos,true); CHKERRQ(ierr);
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

    /*double LxRVE, LyRVE, LzRVE;
    LxRVE=XcoordMax-XcoordMin;
    LyRVE=YcoordMax-YcoordMin;
    LzRVE=ZcoordMax-ZcoordMin;*/

    //Creating Prisims between triangles on -ve and +ve faces
    typedef Face_CenPos_Handle_multiIndex::index<Tri_Hand_tag>::type::iterator Tri_Hand_iterator;
    Tri_Hand_iterator Tri_Neg;
    typedef Face_CenPos_Handle_multiIndex::index<Composite_xyzcoord>::type::iterator xyzcoord_iterator;
    xyzcoord_iterator Tri_Pos;
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
      //insertion of individula prism element and its addition to range preriodic_prisms
      EntityHandle PeriodicPrism;
      rval = moab.create_element(MBPRISM,PrismNodes,6,PeriodicPrism); CHKERR_PETSC(rval);
      preriodic_prisms.insert(PeriodicPrism);

    }

    //insertion of individual prism element and its addition to range preriodic_prisms
    ierr = m_field.seed_ref_level(preriodic_prisms,problem_bit_level); CHKERRQ(ierr);

  }

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
  //    ierr = m_field.get_problem_finite_elements_entities("POTENTIAL_PROBLEM","POTENTIAL_ELEM",out_meshset); CHKERRQ(ierr);
  ierr = m_field.get_entities_by_ref_level(bit_levels.back(),BitRefLevel().set(),out_meshset); CHKERRQ(ierr);
  Range LatestRefinedTets;
  rval = moab.get_entities_by_type(out_meshset, MBTET,LatestRefinedTets,true); CHKERR_PETSC(rval);
  Range LatestRefinedPrisms;
  rval = moab.get_entities_by_type(out_meshset, MBPRISM,LatestRefinedPrisms,true); CHKERR_PETSC(rval);

	Range prims_on_problem_bit_level;
	ierr = m_field.get_entities_by_type_and_ref_level(
    problem_bit_level,BitRefLevel().set(),MBPRISM,prims_on_problem_bit_level
  ); CHKERRQ(ierr);
  //to create meshset from range
  EntityHandle meshset_prims_on_problem_bit_level;
  rval = moab.create_meshset(MESHSET_SET,meshset_prims_on_problem_bit_level); CHKERR_PETSC(rval);
	rval = moab.add_entities(meshset_prims_on_problem_bit_level,prims_on_problem_bit_level); CHKERR_PETSC(rval);
  ierr = m_field.seed_ref_level_MESHSET(meshset_prims_on_problem_bit_level,BitRefLevel().set()); CHKERRQ(ierr);

  // ===========================================================================
  //
  // II. DEFINE PROBLEM
  //
  // ===========================================================================
  
  //Fields
  int field_rank=3;
  ierr = m_field.add_field("DISPLACEMENT",H1,field_rank); CHKERRQ(ierr);
  ierr = m_field.add_ents_to_field_by_TETs(0,"DISPLACEMENT"); CHKERRQ(ierr);
  ierr = m_field.add_field("MESH_NODE_POSITIONS",H1,field_rank,MF_ZERO); CHKERRQ(ierr);
  ierr = m_field.add_ents_to_field_by_TETs(0,"MESH_NODE_POSITIONS"); CHKERRQ(ierr);

  //add entitities (by tets) to the field
  if(choise_value == HOMOBCDISP) {
    ierr = m_field.add_field("LAGRANGE_MUL_DISP",H1,field_rank); CHKERRQ(ierr);
    for(_IT_CUBITMESHSETS_BY_SET_TYPE_FOR_LOOP_(m_field,SIDESET,it)) {
      if(it->get_name().compare(0,12,"AllBoundSurf") == 0 || it->get_msId() == 103) {
        Range tris;
        rval = m_field.get_moab().get_entities_by_type(it->meshset,MBTRI,tris,true); CHKERR_PETSC(rval);
        ierr = m_field.add_ents_to_field_by_TRIs(tris,"LAGRANGE_MUL_DISP",2); CHKERRQ(ierr);
      }
    }
  }

  if(choise_value == HOMOBCPERIODIC) {
    ierr = m_field.add_field("LAGRANGE_MUL_PERIODIC",H1,3); CHKERRQ(ierr);
    //Control 3 rigid body translations in x, y and z axis
    ierr = m_field.add_field("LAGRANGE_MUL_RIGID_TRANS",NOFIELD,3); CHKERRQ(ierr);
    // Setting up dummy no-field vertex
    EntityHandle no_field_vertex;
    {
      const double coords[] = {0,0,0};
      rval = m_field.get_moab().create_vertex(coords,no_field_vertex); CHKERR_PETSC(rval);
      Range range_no_field_vertex;
      range_no_field_vertex.insert(no_field_vertex);
      ierr = m_field.seed_ref_level(range_no_field_vertex,BitRefLevel().set()); CHKERRQ(ierr);
      EntityHandle meshset;
      meshset = m_field.get_field_meshset("LAGRANGE_MUL_RIGID_TRANS");
      rval = m_field.get_moab().add_entities(meshset,range_no_field_vertex); CHKERR_PETSC(rval);
    }
    Range surface_negative;
    ierr = m_field.get_cubit_msId_entities_by_dimension(101,SIDESET,2,surface_negative,true); CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_WORLD,"number of SideSet 101 = %d\n",surface_negative.size()); CHKERRQ(ierr);
    ierr = m_field.add_ents_to_field_by_TRIs(surface_negative,"LAGRANGE_MUL_PERIODIC",2); CHKERRQ(ierr);
  }

  if(choise_value == HOMOBCTRAC) {
    ierr = m_field.add_field("LAGRANGE_MUL_TRAC",NOFIELD,6); CHKERRQ(ierr);
    //Control 3 rigid body translations in x, y and z axis
    ierr = m_field.add_field("LAGRANGE_MUL_RIGID_TRANS",NOFIELD,3); CHKERRQ(ierr);
    //Controla 3 rigid body rotations about x, y and z axis
    ierr = m_field.add_field("LAGRANGE_MUL_RIGID_ROT",NOFIELD,3); CHKERRQ(ierr);
    EntityHandle no_field_vertex;
    {
      const double coords[] = {0,0,0};
      rval = m_field.get_moab().create_vertex(coords,no_field_vertex); CHKERR_PETSC(rval);
      Range range_no_field_vertex;
      range_no_field_vertex.insert(no_field_vertex);
      ierr = m_field.seed_ref_level(range_no_field_vertex,BitRefLevel().set()); CHKERRQ(ierr);
      EntityHandle meshset;
      meshset = m_field.get_field_meshset("LAGRANGE_MUL_TRAC");
      rval = m_field.get_moab().add_entities(meshset,range_no_field_vertex); CHKERR_PETSC(rval);
      meshset = m_field.get_field_meshset("LAGRANGE_MUL_RIGID_TRANS");
      rval = m_field.get_moab().add_entities(meshset,range_no_field_vertex); CHKERR_PETSC(rval);
      meshset = m_field.get_field_meshset("LAGRANGE_MUL_RIGID_ROT");
      rval = m_field.get_moab().add_entities(meshset,range_no_field_vertex); CHKERR_PETSC(rval);
    }
  }

  /*****************************************************************************
   *
   * Set applied order
   * See reference for detals:
   *   Ainsworth M. and Coyle J. (2003) Hierarchic finite element bases on
   unstructured tetrahedral meshes. IJNME, 58(14). pp.2103-2130.
   ****************************************************************************/
  ierr = m_field.set_field_order(0,MBTET,"DISPLACEMENT",order); CHKERRQ(ierr);
  ierr = m_field.set_field_order(0,MBTRI,"DISPLACEMENT",order); CHKERRQ(ierr);
  ierr = m_field.set_field_order(0,MBEDGE,"DISPLACEMENT",order); CHKERRQ(ierr);
  ierr = m_field.set_field_order(0,MBVERTEX,"DISPLACEMENT",1); CHKERRQ(ierr);

  ierr = m_field.set_field_order(0,MBTET,"MESH_NODE_POSITIONS",2); CHKERRQ(ierr);
  ierr = m_field.set_field_order(0,MBTRI,"MESH_NODE_POSITIONS",2); CHKERRQ(ierr);
  ierr = m_field.set_field_order(0,MBEDGE,"MESH_NODE_POSITIONS",2); CHKERRQ(ierr);
  ierr = m_field.set_field_order(0,MBVERTEX,"MESH_NODE_POSITIONS",1); CHKERRQ(ierr);

  PetscBool fo_boundary = PETSC_FALSE;
  ierr = PetscOptionsGetBool(PETSC_NULL,"-my_fo_boundary",&fo_boundary,PETSC_NULL); CHKERRQ(ierr);
  if(fo_boundary) {
    for(_IT_CUBITMESHSETS_BY_SET_TYPE_FOR_LOOP_(m_field,SIDESET,it)) {
      if(it->get_name().compare(0,12,"AllBoundSurf") == 0 || it->get_msId() == 103) {
        Range tris;
        rval = m_field.get_moab().get_entities_by_type(it->meshset,MBTRI,tris,true); CHKERR_PETSC(rval);
        Range tris_edges;
        rval = moab.get_adjacencies(tris,1,false,tris_edges,Interface::UNION); CHKERR_PETSC(rval);
        ierr = m_field.set_field_order(tris,"DISPLACEMENT",1); CHKERRQ(ierr);
        ierr = m_field.set_field_order(tris_edges,"DISPLACEMENT",1); CHKERRQ(ierr);
      }
    }
  }

  if(choise_value == HOMOBCDISP) {
    ierr = m_field.set_field_order(0,MBTRI,"LAGRANGE_MUL_DISP",order); CHKERRQ(ierr);
    ierr = m_field.set_field_order(0,MBEDGE,"LAGRANGE_MUL_DISP",order); CHKERRQ(ierr);
    ierr = m_field.set_field_order(0,MBVERTEX,"LAGRANGE_MUL_DISP",1); CHKERRQ(ierr);
  }

  if(choise_value == HOMOBCPERIODIC) {
    ierr = m_field.set_field_order(0,MBTRI,"LAGRANGE_MUL_PERIODIC",order); CHKERRQ(ierr);
    ierr = m_field.set_field_order(0,MBEDGE,"LAGRANGE_MUL_PERIODIC",order); CHKERRQ(ierr);
    ierr = m_field.set_field_order(0,MBVERTEX,"LAGRANGE_MUL_PERIODIC",1); CHKERRQ(ierr);
  }

  // ===========================================================================
  //
  //  IV. BUILD DATABASE
  //
  // ===========================================================================
  //build field
  ierr = m_field.build_fields(); CHKERRQ(ierr);

  Projection10NodeCoordsOnField ent_method_material(m_field,"MESH_NODE_POSITIONS");
  ierr = m_field.loop_dofs("MESH_NODE_POSITIONS",ent_method_material); CHKERRQ(ierr);

  Hooke<adouble> hooke_adouble;
  Hooke<double> hooke_double;

  NonlinearElasticElement iso_elastic(m_field,1);
  for(_IT_CUBITMESHSETS_BY_BCDATA_TYPE_FOR_LOOP_(m_field,BLOCKSET|MAT_ELASTICSET,it)) {
		if(it->get_name() != "MAT_ELASTIC_1") continue;
    Mat_Elastic mydata;
    ierr = it->get_attribute_data_structure(mydata); CHKERRQ(ierr);
    int id = it->get_msId();
    EntityHandle meshset = it->get_meshset();
    rval = m_field.get_moab().get_entities_by_type(
      meshset,MBTET,iso_elastic.setOfBlocks[id].tEts,true
    ); CHKERR_PETSC(rval);
    iso_elastic.setOfBlocks[id].iD = id;
    iso_elastic.setOfBlocks[id].E = mydata.data.Young;
    iso_elastic.setOfBlocks[id].PoissonRatio = mydata.data.Poisson;
    iso_elastic.setOfBlocks[id].materialDoublePtr = &hooke_double;
    iso_elastic.setOfBlocks[id].materialAdoublePtr = &hooke_adouble;
    ierr = m_field.seed_finite_elements(iso_elastic.setOfBlocks[id].tEts); CHKERRQ(ierr);
  }
  ierr = iso_elastic.addElement("ELASTIC","DISPLACEMENT"); CHKERRQ(ierr);
  ierr = iso_elastic.setOperators("DISPLACEMENT","MESH_NODE_POSITIONS",false,true); CHKERRQ(ierr);
  if(m_field.check_field("POTENTIAL_FIELD")) {
    ierr = m_field.modify_finite_element_add_field_data("ELASTIC","POTENTIAL_FIELD"); CHKERRQ(ierr);
  }

  NonlinearElasticElement trans_elastic(m_field,2);
  trans_elastic.commonData.spatialPositions = "DISPLACEMENT";
  trans_elastic.commonData.meshPositions = "MESH_NODE_POSITIONS";
  boost::ptr_map<int,SmallStrainTranverslyIsotropicADouble *> tranversly_isotropic_adouble_ptr_map;
  boost::ptr_map<int,SmallStrainTranverslyIsotropicDouble *> tranversly_isotropic_double_ptr_map;
  bool trans_iso_blocks = false;
  for(_IT_CUBITMESHSETS_BY_SET_TYPE_FOR_LOOP_(m_field,BLOCKSET,it)) {
    //Get block name
    string name = it->get_name();
    if (name.compare(0,20,"MAT_ELASTIC_TRANSISO") == 0) {
      trans_iso_blocks = true;
      int id = it->get_msId();
      Mat_Elastic_TransIso mydata;
      ierr = it->get_attribute_data_structure(mydata); CHKERRQ(ierr);
      tranversly_isotropic_adouble_ptr_map[id] = new SmallStrainTranverslyIsotropicADouble();
      tranversly_isotropic_double_ptr_map[id] = new SmallStrainTranverslyIsotropicDouble();
      //nu_p, nu_pz, E_p, E_z, G_zp
      tranversly_isotropic_adouble_ptr_map.at(id)->E_p = mydata.data.Youngp;
      tranversly_isotropic_double_ptr_map.at(id)->E_p = mydata.data.Youngp;
      tranversly_isotropic_adouble_ptr_map.at(id)->E_z = mydata.data.Youngz;
      tranversly_isotropic_double_ptr_map.at(id)->E_z = mydata.data.Youngz;
      tranversly_isotropic_adouble_ptr_map.at(id)->nu_p = mydata.data.Poissonp;
      tranversly_isotropic_double_ptr_map.at(id)->nu_p = mydata.data.Poissonp;
      tranversly_isotropic_adouble_ptr_map.at(id)->nu_pz = mydata.data.Poissonpz;
      tranversly_isotropic_double_ptr_map.at(id)->nu_pz = mydata.data.Poissonpz;
      double shear_zp;
      if(mydata.data.Shearzp!=0) {
        shear_zp = mydata.data.Shearzp;
      } else {
        shear_zp = mydata.data.Youngz/(2*(1+mydata.data.Poissonpz));
      }
      tranversly_isotropic_adouble_ptr_map.at(it->get_msId())->G_zp = shear_zp;
      tranversly_isotropic_double_ptr_map.at(it->get_msId())->G_zp = shear_zp;
      //get tets from block where material is defined
      EntityHandle meshset = it->get_meshset();
      rval = m_field.get_moab().get_entities_by_type(
        meshset,MBTET,trans_elastic.setOfBlocks[id].tEts,true
      ); CHKERR_PETSC(rval);
      //adding material to nonlinear class
      trans_elastic.setOfBlocks[id].iD = id;
      //note that material parameters are defined internally in material model
      trans_elastic.setOfBlocks[id].E = 0; // this is not working for this material
      trans_elastic.setOfBlocks[id].PoissonRatio = 0; // this is not working for this material
      trans_elastic.setOfBlocks[id].materialDoublePtr = tranversly_isotropic_double_ptr_map.at(id);
      trans_elastic.setOfBlocks[id].materialAdoublePtr = tranversly_isotropic_adouble_ptr_map.at(id);
      ierr = m_field.seed_finite_elements(trans_elastic.setOfBlocks[id].tEts); CHKERRQ(ierr);
    }
  }
  if(trans_iso_blocks) {
    ierr = m_field.add_finite_element("TRAN_ISOTROPIC_ELASTIC"); CHKERRQ(ierr);
    ierr = m_field.add_finite_element("TRAN_ISOTROPIC_ELASTIC",MF_ZERO); CHKERRQ(ierr);
    ierr = m_field.modify_finite_element_add_field_row("TRAN_ISOTROPIC_ELASTIC","DISPLACEMENT"); CHKERRQ(ierr);
    ierr = m_field.modify_finite_element_add_field_col("TRAN_ISOTROPIC_ELASTIC","DISPLACEMENT"); CHKERRQ(ierr);
    ierr = m_field.modify_finite_element_add_field_data("TRAN_ISOTROPIC_ELASTIC","DISPLACEMENT"); CHKERRQ(ierr);
    ierr = m_field.modify_finite_element_add_field_data("TRAN_ISOTROPIC_ELASTIC","POTENTIAL_FIELD"); CHKERRQ(ierr);
    if(m_field.check_field("MESH_NODE_POSITIONS")) {
      ierr = m_field.modify_finite_element_add_field_data("TRAN_ISOTROPIC_ELASTIC","MESH_NODE_POSITIONS"); CHKERRQ(ierr);
    }
    for(
      map<int,NonlinearElasticElement::BlockData>::iterator sit = trans_elastic.setOfBlocks.begin();
      sit!=trans_elastic.setOfBlocks.end();sit++
    ) {
      ierr = m_field.add_ents_to_finite_element_by_TETs(sit->second.tEts,"TRAN_ISOTROPIC_ELASTIC"); CHKERRQ(ierr);
    }
  }
  if(trans_iso_blocks) {
    //Rhs
    trans_elastic.feRhs.getOpPtrVector().push_back(
      new NonlinearElasticElement::OpGetCommonDataAtGaussPts("DISPLACEMENT",trans_elastic.commonData)
    );
    trans_elastic.feRhs.getOpPtrVector().push_back(
      new NonlinearElasticElement::OpGetCommonDataAtGaussPts("POTENTIAL_FIELD",trans_elastic.commonData)
    );
    if(m_field.check_field("MESH_NODE_POSITIONS")) {
      trans_elastic.feRhs.getOpPtrVector().push_back(
        new NonlinearElasticElement::OpGetCommonDataAtGaussPts("MESH_NODE_POSITIONS",trans_elastic.commonData)
      );
    }
    map<int,NonlinearElasticElement::BlockData>::iterator sit = trans_elastic.setOfBlocks.begin();
    for(;sit!=trans_elastic.setOfBlocks.end();sit++) {
      trans_elastic.feRhs.getOpPtrVector().push_back(
        new NonlinearElasticElement::OpJacobianPiolaKirchhoffStress(
          "DISPLACEMENT",sit->second,trans_elastic.commonData,2,false,false,true
        )
      );
      trans_elastic.feRhs.getOpPtrVector().push_back(
        new NonlinearElasticElement::OpRhsPiolaKirchhoff(
          "DISPLACEMENT",sit->second,trans_elastic.commonData
        )
      );
    }

    //Lhs
    trans_elastic.feLhs.getOpPtrVector().push_back(
      new NonlinearElasticElement::OpGetCommonDataAtGaussPts("DISPLACEMENT",trans_elastic.commonData)
    );
    trans_elastic.feLhs.getOpPtrVector().push_back(
      new NonlinearElasticElement::OpGetCommonDataAtGaussPts("POTENTIAL_FIELD",trans_elastic.commonData)
    );
    if(m_field.check_field("MESH_NODE_POSITIONS")) {
      trans_elastic.feLhs.getOpPtrVector().push_back(
        new NonlinearElasticElement::OpGetCommonDataAtGaussPts("MESH_NODE_POSITIONS",trans_elastic.commonData)
      );
    }
    sit = trans_elastic.setOfBlocks.begin();
    for(;sit!=trans_elastic.setOfBlocks.end();sit++) {
      trans_elastic.feLhs.getOpPtrVector().push_back(
        new NonlinearElasticElement::OpJacobianPiolaKirchhoffStress(
          "DISPLACEMENT",sit->second,trans_elastic.commonData,2,true,false,true
        )
      );
      trans_elastic.feLhs.getOpPtrVector().push_back(
        new NonlinearElasticElement::OpLhsPiolaKirchhoff_dx(
          "DISPLACEMENT","DISPLACEMENT",sit->second,trans_elastic.commonData
        )
      );
    }
  }

  BCs_RVELagrange_Disp lagrangian_element_disp(m_field);
  if(choise_value == HOMOBCDISP) {
    lagrangian_element_disp.addLagrangiangElement(
      "LAGRANGE_ELEM","DISPLACEMENT","LAGRANGE_MUL_DISP","MESH_NODE_POSITIONS"
    );
  }

  BCs_RVELagrange_Trac lagrangian_element_trac(m_field);
  BCs_RVELagrange_Trac_Rigid_Trans lagrangian_element_rigid_body_trans(m_field);
  BCs_RVELagrange_Trac_Rigid_Rot lagrangian_element_rigid_body_rot(m_field);
  if(choise_value == HOMOBCTRAC) {
    lagrangian_element_trac.addLagrangiangElement(
      "LAGRANGE_ELEM","DISPLACEMENT","LAGRANGE_MUL_TRAC","MESH_NODE_POSITIONS"
    );
    lagrangian_element_trac.addLagrangiangElement(
      "LAGRANGE_ELEM_TRANS","DISPLACEMENT","LAGRANGE_MUL_RIGID_TRANS","MESH_NODE_POSITIONS"
    );
    lagrangian_element_trac.addLagrangiangElement(
      "LAGRANGE_ELEM_ROT","DISPLACEMENT","LAGRANGE_MUL_RIGID_ROT","MESH_NODE_POSITIONS"
    );
  }

  BCs_RVELagrange_Periodic lagrangian_element_periodic(m_field);
  if(choise_value == HOMOBCPERIODIC) {
    lagrangian_element_periodic.addLagrangiangElement(
      "LAGRANGE_ELEM","DISPLACEMENT","LAGRANGE_MUL_PERIODIC","MESH_NODE_POSITIONS",preriodic_prisms
    );
    lagrangian_element_trac.addLagrangiangElement(
      "LAGRANGE_ELEM_TRANS","DISPLACEMENT","LAGRANGE_MUL_RIGID_TRANS","MESH_NODE_POSITIONS"
    );
  }

  struct MinMaxNodes {
    enum MINAMX { C0,MAXLAST };
    EntityHandle entMinMax[MAXLAST];
    ublas::vector<int> rowIndices;
    VectorDouble rHs[6];
    MinMaxNodes() {
      rowIndices.resize(3*MAXLAST);
      for(int ii = 0;ii<6;ii++) {
        rHs[ii].resize(3*MAXLAST);
      }
    }
  };
  MinMaxNodes minMaxNodes;

  if(choise_value == NITSCHE) { // Condensed traction/periodc BC
    ierr = m_field.add_finite_element("SURFACE_ELEMENTS"); CHKERRQ(ierr);
    ierr = m_field.modify_finite_element_add_field_row("SURFACE_ELEMENTS","DISPLACEMENT"); CHKERRQ(ierr);
    ierr = m_field.modify_finite_element_add_field_col("SURFACE_ELEMENTS","DISPLACEMENT"); CHKERRQ(ierr);
    ierr = m_field.modify_finite_element_add_field_data("SURFACE_ELEMENTS","DISPLACEMENT"); CHKERRQ(ierr);
    ierr = m_field.modify_finite_element_add_field_data("SURFACE_ELEMENTS","MESH_NODE_POSITIONS"); CHKERRQ(ierr);
    EntityHandle condensed_traction_element_meshset;
    rval = moab.create_meshset(MESHSET_TRACK_OWNER,condensed_traction_element_meshset); CHKERRQ(ierr);
    Range nodes;
    for(_IT_CUBITMESHSETS_BY_SET_TYPE_FOR_LOOP_(m_field,SIDESET,it)) {
      if(it->get_name().compare(0,12,"AllBoundSurf") == 0 || it->get_msId() == 103) {
        Range tris;
        rval = m_field.get_moab().get_entities_by_type(it->meshset,MBTRI,tris,true); CHKERR_PETSC(rval);
        ierr = m_field.add_ents_to_finite_element_by_TRIs(tris,"SURFACE_ELEMENTS"); CHKERRQ(ierr);
        Range tris_nodes;
        rval = moab.get_connectivity(tris,nodes,true); CHKERR_PETSC(rval);
        nodes.merge(tris_nodes);
      }
    }

    {
      ublas::vector<double> x,y,z;
      x.resize(nodes.size(),false);
      y.resize(nodes.size(),false);
      z.resize(nodes.size(),false);
      rval = moab.get_coords(nodes,&x[0],&y[0],&z[0]); CHKERR_PETSC(rval);
      int bc_nb = 0;
      for(int sx = -1; sx<=+1; sx+=2) {
        for(int sy = -1; sy<=+1; sy+=2) {
          for(int sz = -1; sz<=+1; sz+=2) {
            if(bc_nb == MinMaxNodes::MAXLAST) break;
            ublas::vector<double> dist_up_right;
            dist_up_right.resize(x.size(),false);
            dist_up_right.clear();
            for(unsigned int nn = 0;nn<x.size();nn++) {
              if(
                ((sx*x[nn])>0)&&
                ((sy*y[nn])>0)&&
                ((sz*z[nn])>0)
              ) {
                dist_up_right[nn] = sx*x[nn]+sy*y[nn]+sz*z[nn];
              }
            }
            ublas::vector<double>::iterator dist_it;
            dist_it = max_element(dist_up_right.begin(),dist_up_right.end());
            minMaxNodes.entMinMax[bc_nb++] = nodes[distance(dist_up_right.begin(),dist_it)];
          }
        }
      }
    }

  }

  //build finite elements
  ierr = m_field.build_finite_elements(); CHKERRQ(ierr);
  //build adjacencies
  ierr = m_field.build_adjacencies(problem_bit_level); CHKERRQ(ierr);

  //define problems
  ierr = m_field.add_problem("ELASTIC_MECHANICS"); CHKERRQ(ierr);
  //set finite elements for problem
  ierr = m_field.modify_problem_add_finite_element("ELASTIC_MECHANICS","ELASTIC"); CHKERRQ(ierr);
  if(trans_iso_blocks) {
    ierr = m_field.modify_problem_add_finite_element(
      "ELASTIC_MECHANICS","TRAN_ISOTROPIC_ELASTIC"
    ); CHKERRQ(ierr);
  }
  if(choise_value == HOMOBCDISP) {
    ierr = m_field.modify_problem_add_finite_element("ELASTIC_MECHANICS","LAGRANGE_ELEM"); CHKERRQ(ierr);
  }
  if(choise_value == HOMOBCTRAC) {
    ierr = m_field.modify_problem_add_finite_element("ELASTIC_MECHANICS","LAGRANGE_ELEM"); CHKERRQ(ierr);
    ierr = m_field.modify_problem_add_finite_element("ELASTIC_MECHANICS","LAGRANGE_ELEM_TRANS"); CHKERRQ(ierr);
    ierr = m_field.modify_problem_add_finite_element("ELASTIC_MECHANICS","LAGRANGE_ELEM_ROT"); CHKERRQ(ierr);
  }
  if(choise_value == HOMOBCPERIODIC) {
    ierr = m_field.modify_problem_add_finite_element("ELASTIC_MECHANICS","LAGRANGE_ELEM"); CHKERRQ(ierr);
    ierr = m_field.modify_problem_add_finite_element("ELASTIC_MECHANICS","LAGRANGE_ELEM_TRANS"); CHKERRQ(ierr);
  }
  if(choise_value == NITSCHE) {
    ierr = m_field.modify_problem_add_finite_element("ELASTIC_MECHANICS","SURFACE_ELEMENTS"); CHKERRQ(ierr);
  }

  //set refinment level for problem
  ierr = m_field.modify_problem_ref_level_add_bit("ELASTIC_MECHANICS",problem_bit_level); CHKERRQ(ierr);
  //build problem
  ierr = m_field.build_problems(); CHKERRQ(ierr);

  
  // ===========================================================================
  //
  //  V. MESH PARTITION
  //
  // ===========================================================================
  
  //partition
  ierr = m_field.partition_problem("ELASTIC_MECHANICS"); CHKERRQ(ierr);
  ierr = m_field.partition_finite_elements(
    "ELASTIC_MECHANICS",false,0,m_field.getCommSize() // build elements on all procs
  ); CHKERRQ(ierr);
  //what are ghost nodes, see Petsc Manual
  ierr = m_field.partition_ghost_dofs("ELASTIC_MECHANICS"); CHKERRQ(ierr);

  
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
  Vec D;
  vector<Vec> F(6);
  ierr = m_field.VecCreateGhost("ELASTIC_MECHANICS",ROW,&F[0]); CHKERRQ(ierr);
  for(int ii = 1;ii<6;ii++) {
    ierr = VecDuplicate(F[0],&F[ii]); CHKERRQ(ierr);
  }
  ierr = m_field.VecCreateGhost("ELASTIC_MECHANICS",COL,&D); CHKERRQ(ierr);

  /*****************************************************************************
   *
   *  1. Assembling global stiffness matrix K
   *     and external force vector F
   ****************************************************************************/
  Mat Aij;
  ierr = m_field.MatCreateMPIAIJWithArrays("ELASTIC_MECHANICS",&Aij); CHKERRQ(ierr);
  ierr = MatSetOption(Aij,MAT_STRUCTURALLY_SYMMETRIC,PETSC_TRUE); CHKERRQ(ierr);
  ierr = MatSetOption(Aij,MAT_NEW_NONZERO_LOCATIONS,PETSC_TRUE); CHKERRQ(ierr);
  ierr = MatSetOption(Aij,MAT_KEEP_NONZERO_PATTERN,PETSC_TRUE); CHKERRQ(ierr);
  ierr = MatSetOption(Aij,MAT_USE_INODES,PETSC_TRUE); CHKERRQ(ierr);
  ierr = MatSetOption(Aij,MAT_NEW_NONZERO_LOCATION_ERR,PETSC_FALSE); CHKERRQ(ierr);
  ierr = MatSetOption(Aij,MAT_KEEP_NONZERO_PATTERN,PETSC_FALSE); CHKERRQ(ierr);

  /*{
    ierr = MatView(Aij,PETSC_VIEWER_DRAW_SELF); CHKERRQ(ierr);
    std::string wait;
    std::cin >> wait;
  }*/

  ierr = VecZeroEntries(D); CHKERRQ(ierr);
  ierr = VecGhostUpdateBegin(D,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = VecGhostUpdateEnd(D,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = m_field.set_global_ghost_vector(
    "ELASTIC_MECHANICS",ROW,D,INSERT_VALUES,SCATTER_REVERSE
  ); CHKERRQ(ierr);
  for(int ii = 0;ii<6;ii++) {
    ierr = VecZeroEntries(F[ii]); CHKERRQ(ierr);
    ierr = VecGhostUpdateBegin(F[ii],INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
    ierr = VecGhostUpdateEnd(F[ii],INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  }
  ierr = MatZeroEntries(Aij); CHKERRQ(ierr);

  NitscheMethod::BlockData nitsche_block_data;
  NitscheMethod::CommonData nitsche_common_data;
  PeriodicNitscheConstrains::CommonData periodic_common_data;
  PeriodicNitscheConstrains::MyNitscheVolume nitsche_element_iso(
    m_field,nitsche_block_data,nitsche_common_data,periodic_common_data
  );
  PeriodicNitscheConstrains::MyNitscheVolume nitsche_element_trans(
    m_field,nitsche_block_data,nitsche_common_data,periodic_common_data
  );

  if(choise_value == NITSCHE) {

    nitsche_block_data.faceElemName = "SURFACE_ELEMENTS";
    for(_IT_CUBITMESHSETS_BY_SET_TYPE_FOR_LOOP_(m_field,SIDESET,it)) {
      if(it->get_name().compare(0,12,"AllBoundSurf") == 0 || it->get_msId() == 103) {
        Range tris;
        rval = m_field.get_moab().get_entities_by_type(it->meshset,MBTRI,tris,true); CHKERR_PETSC(rval);
        nitsche_block_data.fAces.merge(tris);
      }
    }

    for(_IT_CUBITMESHSETS_BY_SET_TYPE_FOR_LOOP_(m_field,SIDESET,it)) {
      if(it->get_name().compare(0,12,"AllBoundSurf") == 0 || it->get_msId() == 103) {
        Range tris;
        rval = m_field.get_moab().get_entities_by_type(it->meshset,MBTRI,tris,true); CHKERR_PETSC(rval);
        periodic_common_data.skinFaces.merge(tris);
      }
    }

    nitsche_block_data.gamma = 1e-7;
    nitsche_block_data.phi = -1;
    periodic_common_data.eps = 0;
    ierr = PetscOptionsGetReal(
      PETSC_NULL,"-my_gamma",&nitsche_block_data.gamma,PETSC_NULL
    ); CHKERRQ(ierr);
    ierr = PetscOptionsGetReal(
      PETSC_NULL,"-my_phi",&nitsche_block_data.phi,PETSC_NULL
    ); CHKERRQ(ierr);
    ierr = PetscOptionsGetReal(
      PETSC_NULL,"-my_eps",&periodic_common_data.eps,PETSC_NULL
    ); CHKERRQ(ierr);
    ierr = PetscPrintf(
      PETSC_COMM_WORLD,
      "Nitsche method gamma = %4.2e phi = %2.1f eps = %4.2e\n",
      nitsche_block_data.gamma,nitsche_block_data.phi,periodic_common_data.eps
    ); CHKERRQ(ierr);

    for(
      map<int,NonlinearElasticElement::BlockData>::iterator mit = iso_elastic.setOfBlocks.begin();
      mit!=iso_elastic.setOfBlocks.end();
      mit++
    ) {
      NonlinearElasticElement::CommonData &elastic_common_data = iso_elastic.commonData;
      NonlinearElasticElement::BlockData &elastic_block_data = mit->second;
      nitsche_element_iso.snes_B = Aij;

      nitsche_element_iso.getOpPtrVector().push_back(
        new NonlinearElasticElement::OpGetCommonDataAtGaussPts(
          "DISPLACEMENT",elastic_common_data
        )
      );
      nitsche_element_iso.getOpPtrVector().push_back(
        new NonlinearElasticElement::OpGetCommonDataAtGaussPts(
          "MESH_NODE_POSITIONS",elastic_common_data
        )
      );
      nitsche_element_iso.getOpPtrVector().push_back(
        new NonlinearElasticElement::OpJacobianPiolaKirchhoffStress(
          "DISPLACEMENT",elastic_block_data,elastic_common_data,1,true,false,true
        )
      );
      nitsche_element_iso.getOpPtrVector().push_back(
        new PeriodicNitscheConstrains::OpLhsPeriodicNormal(
          "DISPLACEMENT",nitsche_block_data,nitsche_common_data,
          elastic_block_data,elastic_common_data,
          periodic_common_data
        )
      );
      nitsche_element_iso.getOpPtrVector().push_back(
        new PeriodicNitscheConstrains::OpRhsPeriodicNormal(
          "DISPLACEMENT",nitsche_block_data,nitsche_common_data,
          elastic_block_data,elastic_common_data,
          periodic_common_data,
          F
        )
      );
      nitsche_element_iso.getOpPtrVector().push_back(
        new NitscheMethod::OpLhsNormal(
          "DISPLACEMENT",nitsche_block_data,nitsche_common_data,
          elastic_block_data,elastic_common_data,true
        )
      );

      // this is to get data on opposite element
      nitsche_element_iso.periodicVolume.getOpPtrVector().clear();
      nitsche_element_iso.periodicVolume.getOpPtrVector().push_back(
        new NonlinearElasticElement::OpGetCommonDataAtGaussPts(
          "DISPLACEMENT",elastic_common_data
        )
      );
      nitsche_element_iso.periodicVolume.getOpPtrVector().push_back(
        new NonlinearElasticElement::OpGetCommonDataAtGaussPts(
          "MESH_NODE_POSITIONS",elastic_common_data
        )
      );
      nitsche_element_iso.periodicVolume.getOpPtrVector().push_back(
        new NonlinearElasticElement::OpJacobianPiolaKirchhoffStress(
          "DISPLACEMENT",elastic_block_data,elastic_common_data,1,true,false,true
        )
      );
      nitsche_element_iso.periodicVolume.getOpPtrVector().push_back(
        new PeriodicNitscheConstrains::OpGetVolumeData(
          elastic_common_data,
          periodic_common_data
        )
      );
      periodic_common_data.volumeElemName = "ELASTIC";

      nitsche_element_iso.addToRule = 1;
      ierr = m_field.loop_finite_elements("ELASTIC_MECHANICS","ELASTIC",nitsche_element_iso); CHKERRQ(ierr);
    }

    for(
      map<int,NonlinearElasticElement::BlockData>::iterator mit = trans_elastic.setOfBlocks.begin();
      mit!=trans_elastic.setOfBlocks.end();
      mit++
    ) {
      NonlinearElasticElement::CommonData &elastic_common_data = trans_elastic.commonData;
      NonlinearElasticElement::BlockData &elastic_block_data = mit->second;
      nitsche_element_trans.snes_B = Aij;

      nitsche_element_trans.getOpPtrVector().push_back(
        new NonlinearElasticElement::OpGetCommonDataAtGaussPts(
          "DISPLACEMENT",elastic_common_data
        )
      );
      nitsche_element_trans.getOpPtrVector().push_back(
        new NonlinearElasticElement::OpGetCommonDataAtGaussPts(
          "MESH_NODE_POSITIONS",elastic_common_data
        )
      );
      nitsche_element_trans.getOpPtrVector().push_back(
        new NonlinearElasticElement::OpGetCommonDataAtGaussPts(
          "POTENTIAL_FIELD",elastic_common_data
        )
      );
      nitsche_element_trans.getOpPtrVector().push_back(
        new NonlinearElasticElement::OpJacobianPiolaKirchhoffStress(
          "DISPLACEMENT",elastic_block_data,elastic_common_data,2,true,false,true
        )
      );
      nitsche_element_trans.getOpPtrVector().push_back(
        new PeriodicNitscheConstrains::OpLhsPeriodicNormal(
          "DISPLACEMENT",nitsche_block_data,nitsche_common_data,
          elastic_block_data,elastic_common_data,
          periodic_common_data
        )
      );
      nitsche_element_trans.getOpPtrVector().push_back(
        new PeriodicNitscheConstrains::OpRhsPeriodicNormal(
          "DISPLACEMENT",nitsche_block_data,nitsche_common_data,
          elastic_block_data,elastic_common_data,
          periodic_common_data,
          F
        )
      );
      nitsche_element_trans.getOpPtrVector().push_back(
        new NitscheMethod::OpLhsNormal(
          "DISPLACEMENT",nitsche_block_data,nitsche_common_data,
          elastic_block_data,elastic_common_data,true
        )
      );

      // this is to get data on opposite element
      nitsche_element_trans.periodicVolume.getOpPtrVector().clear();
      nitsche_element_trans.periodicVolume.getOpPtrVector().push_back(
        new NonlinearElasticElement::OpGetCommonDataAtGaussPts(
          "DISPLACEMENT",elastic_common_data
        )
      );
      nitsche_element_trans.periodicVolume.getOpPtrVector().push_back(
        new NonlinearElasticElement::OpGetCommonDataAtGaussPts(
          "POTENTIAL_FIELD",elastic_common_data
        )
      );
      nitsche_element_trans.periodicVolume.getOpPtrVector().push_back(
        new NonlinearElasticElement::OpGetCommonDataAtGaussPts(
          "MESH_NODE_POSITIONS",elastic_common_data
        )
      );
      nitsche_element_trans.periodicVolume.getOpPtrVector().push_back(
        new NonlinearElasticElement::OpJacobianPiolaKirchhoffStress(
          "DISPLACEMENT",elastic_block_data,elastic_common_data,2,true,false,true
        )
      );
      nitsche_element_trans.periodicVolume.getOpPtrVector().push_back(
        new PeriodicNitscheConstrains::OpGetVolumeData(
          elastic_common_data,
          periodic_common_data
        )
      );
      periodic_common_data.volumeElemName = "TRAN_ISOTROPIC_ELASTIC";

      nitsche_element_trans.addToRule = 1;
      ierr = m_field.loop_finite_elements(
        "ELASTIC_MECHANICS","TRAN_ISOTROPIC_ELASTIC",nitsche_element_trans
      );  CHKERRQ(ierr);
    }


  }

  /*****************************************************************************
   *
   *  2. Get the volume of RVE
   *
   ****************************************************************************/
  Vec volume_vec;
  int volume_vec_ghost[] = { 0 };
  ierr = VecCreateGhost(
    PETSC_COMM_WORLD,(!m_field.getCommRank())?1:0,1,1,volume_vec_ghost,&volume_vec
  );  CHKERRQ(ierr);
  ierr = VecZeroEntries(volume_vec); CHKERRQ(ierr);

  iso_elastic.getLoopFeLhs().getOpPtrVector().push_back(new VolumeCalculation("DISPLACEMENT",volume_vec));
  trans_elastic.getLoopFeLhs().getOpPtrVector().push_back(new VolumeCalculation("DISPLACEMENT",volume_vec));

  /*****************************************************************************
   *
   *  3. SOLVE THE ZEROTH-ORDER FINITE ELEMENT EQUILIBRIUM EQUATION
   *     [K][U] = [F]
   *
   ****************************************************************************/
  //iso_elastic element matrix
  iso_elastic.getLoopFeLhs().snes_x = D;
  iso_elastic.getLoopFeLhs().snes_B = Aij;
  trans_elastic.getLoopFeLhs().snes_x = D;
  trans_elastic.getLoopFeLhs().snes_B = Aij;
  ierr = m_field.loop_finite_elements("ELASTIC_MECHANICS","ELASTIC",iso_elastic.getLoopFeLhs());  CHKERRQ(ierr);
  ierr = m_field.loop_finite_elements("ELASTIC_MECHANICS","TRAN_ISOTROPIC_ELASTIC",trans_elastic.getLoopFeLhs());  CHKERRQ(ierr);

  if(choise_value == HOMOBCDISP) {
    lagrangian_element_disp.setRVEBCsOperators("DISPLACEMENT","LAGRANGE_MUL_DISP","MESH_NODE_POSITIONS",Aij,F);
    ierr = m_field.loop_finite_elements("ELASTIC_MECHANICS","LAGRANGE_ELEM",lagrangian_element_disp.getLoopFeRVEBCLhs()); CHKERRQ(ierr);
    ierr = m_field.loop_finite_elements("ELASTIC_MECHANICS","LAGRANGE_ELEM",lagrangian_element_disp.getLoopFeRVEBCRhs()); CHKERRQ(ierr);
  }
  if(choise_value == HOMOBCTRAC) {
    lagrangian_element_trac.setRVEBCsOperators("DISPLACEMENT","LAGRANGE_MUL_TRAC","MESH_NODE_POSITIONS",Aij,F);
    ierr = m_field.loop_finite_elements("ELASTIC_MECHANICS","LAGRANGE_ELEM",lagrangian_element_trac.getLoopFeRVEBCLhs()); CHKERRQ(ierr);
    ierr = m_field.loop_finite_elements("ELASTIC_MECHANICS","LAGRANGE_ELEM",lagrangian_element_trac.getLoopFeRVEBCRhs()); CHKERRQ(ierr);
    lagrangian_element_rigid_body_trans.setRVEBCsRigidBodyTranOperators(
      "DISPLACEMENT","LAGRANGE_MUL_RIGID_TRANS",Aij,lagrangian_element_trac.setOfRVEBC
    );
    ierr = m_field.loop_finite_elements(
      "ELASTIC_MECHANICS","LAGRANGE_ELEM_TRANS",lagrangian_element_rigid_body_trans.getLoopFeRVEBCLhs()
    ); CHKERRQ(ierr);
    lagrangian_element_rigid_body_rot.setRVEBCsRigidBodyRotOperators(
      "DISPLACEMENT","LAGRANGE_MUL_RIGID_ROT",Aij,lagrangian_element_trac.setOfRVEBC
    );
    ierr = m_field.loop_finite_elements(
      "ELASTIC_MECHANICS","LAGRANGE_ELEM_ROT",lagrangian_element_rigid_body_rot.getLoopFeRVEBCLhs()
    ); CHKERRQ(ierr);
  }
  if(choise_value == HOMOBCPERIODIC) {
    lagrangian_element_periodic.setRVEBCsOperators("DISPLACEMENT","LAGRANGE_MUL_PERIODIC","MESH_NODE_POSITIONS",Aij,F);
    ierr = m_field.loop_finite_elements(
      "ELASTIC_MECHANICS","LAGRANGE_ELEM",lagrangian_element_periodic.getLoopFeRVEBCLhs()
    ); CHKERRQ(ierr);
    ierr = m_field.loop_finite_elements(
      "ELASTIC_MECHANICS","LAGRANGE_ELEM",lagrangian_element_periodic.getLoopFeRVEBCRhs()
    ); CHKERRQ(ierr);
    lagrangian_element_rigid_body_trans.setRVEBCsRigidBodyTranOperators(
      "DISPLACEMENT","LAGRANGE_MUL_RIGID_TRANS",Aij,lagrangian_element_periodic.setOfRVEBC
    );
    ierr = m_field.loop_finite_elements(
      "ELASTIC_MECHANICS","LAGRANGE_ELEM_TRANS",lagrangian_element_rigid_body_trans.getLoopFeRVEBCLhs()
    ); CHKERRQ(ierr);
  }

  for(int ii = 0;ii<6;ii++) {
    ierr = VecGhostUpdateBegin(F[ii],ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
    ierr = VecGhostUpdateEnd(F[ii],ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
    ierr = VecAssemblyBegin(F[ii]); CHKERRQ(ierr);
    ierr = VecAssemblyEnd(F[ii]); CHKERRQ(ierr);
  }
  ierr = MatAssemblyBegin(Aij,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
  ierr = MatAssemblyEnd(Aij,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);

  {
    //ierr = MatView(Aij,PETSC_VIEWER_DRAW_SELF); CHKERRQ(ierr);
    //std::string wait;
    //std::cin >> wait;
  }

  if(choise_value == NITSCHE) {
    if(periodic_common_data.eps==0) {
      const MoFEMProblem *problem_ptr;
      ierr = m_field.get_problem("ELASTIC_MECHANICS",&problem_ptr); CHKERRQ(ierr);
      for(int nn = 0;nn!=MinMaxNodes::MAXLAST;nn++) {
        for(
          _IT_NUMEREDDOFMOFEMENTITY_ROW_BY_ENT_FOR_LOOP_(
            problem_ptr,minMaxNodes.entMinMax[nn],dof
          )
        ) {
          minMaxNodes.rowIndices[3*nn+dof->get_dof_coeff_idx()]
          = dof->get_petsc_gloabl_dof_idx();
        }
      }
      VectorDouble coords;
      int nb_bcs = 3*MinMaxNodes::MAXLAST;
      coords.resize(nb_bcs,false);
      rval = moab.get_coords(&minMaxNodes.entMinMax[0],MinMaxNodes::MAXLAST,&coords[0]);
      VectorDouble strain;
      strain.resize(6,false);
      MatrixDouble mat_d;
      for(int ii = 0;ii<6;ii++) {
        minMaxNodes.rHs[ii].clear();
        strain.clear();
        strain[ii] = 1;
        for(int nn = 0;nn<MinMaxNodes::MAXLAST;nn++) {
          PeriodicNitscheConstrains::OpRhsPeriodicNormal::calcualteDMatrix(
            VectorAdaptor(3,ublas::shallow_array_adaptor<double>(3,&coords[3*nn])),
            mat_d
          );
          VectorAdaptor rhs(3,ublas::shallow_array_adaptor<double>(3,&minMaxNodes.rHs[ii][3*nn]));
          noalias(rhs) = prod(mat_d,strain);
        }
      }
      for(int ii = 0;ii<6;ii++) {
        ierr = VecSetValues(
          F[ii],nb_bcs,&minMaxNodes.rowIndices[0],&minMaxNodes.rHs[ii][0],INSERT_VALUES
        ); CHKERRQ(ierr);
        ierr = VecAssemblyBegin(F[ii]); CHKERRQ(ierr);
        ierr = VecAssemblyEnd(F[ii]); CHKERRQ(ierr);
      }
      ierr = MatZeroRows(
        Aij,nb_bcs,&minMaxNodes.rowIndices[0],1,PETSC_NULL,PETSC_NULL
      ); CHKERRQ(ierr);
    }
  }

  ierr = VecAssemblyBegin(volume_vec); CHKERRQ(ierr);
  ierr = VecAssemblyEnd(volume_vec); CHKERRQ(ierr);
  double rve_volume;
  ierr = VecSum(volume_vec,&rve_volume);  CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"RVE Volume %3.2g\n",rve_volume); CHKERRQ(ierr);

  ierr = VecDestroy(&volume_vec);

  //----------------------------------------------------------------------------
  // 3.1 Solving the equation to get nodal displacement
  //----------------------------------------------------------------------------
  // Solver
  KSP solver;
  ierr = KSPCreate(PETSC_COMM_WORLD,&solver); CHKERRQ(ierr);
  ierr = KSPSetOperators(solver,Aij,Aij); CHKERRQ(ierr);
  ierr = KSPSetFromOptions(solver); CHKERRQ(ierr);
  ierr = KSPSetUp(solver); CHKERRQ(ierr);

  ublas::matrix<double> Dmat;
  Dmat.resize(6,6);
  Dmat.clear();

  //create a vector for 6 components of homogenized stress
  Vec stress_homo;
  int stress_homo_ghost[] = { 0,1,2,3,4,5,6 };
  NonlinearElasticElement::MyVolumeFE ave_stress_iso(m_field);
  NonlinearElasticElement::MyVolumeFE ave_stress_trans(m_field);
  PetscBool stress_by_boundary_integral = PETSC_FALSE;
  ierr = VecCreateGhost(
    PETSC_COMM_WORLD,(!m_field.getCommRank())?6:0,6,6,stress_homo_ghost,&stress_homo
  );  CHKERRQ(ierr);
  switch(choise_value) {
    case HOMOBCDISP:
    lagrangian_element_disp.setRVEBCsHomoStressOperators(
      "DISPLACEMENT","LAGRANGE_MUL_DISP","MESH_NODE_POSITIONS",stress_homo
    );
    break;
    case HOMOBCTRAC:
    lagrangian_element_trac.setRVEBCsHomoStressOperators(
      "DISPLACEMENT","LAGRANGE_MUL_TRAC","MESH_NODE_POSITIONS",stress_homo
    );
    break;
    case HOMOBCPERIODIC:
    lagrangian_element_periodic.setRVEBCsHomoStressOperators(
      "DISPLACEMENT","LAGRANGE_MUL_PERIODIC","MESH_NODE_POSITIONS",stress_homo
    );
    break;
    case NITSCHE:
    if(stress_by_boundary_integral) {
      nitsche_element_iso.getOpPtrVector().clear();
      nitsche_element_trans.getOpPtrVector().clear();
      nitsche_element_iso.periodicVolume.getOpPtrVector().clear();
      nitsche_element_trans.periodicVolume.getOpPtrVector().clear();
      for(
        map<int,NonlinearElasticElement::BlockData>::iterator mit = iso_elastic.setOfBlocks.begin();
        mit!=iso_elastic.setOfBlocks.end();
        mit++
      ) {
        NonlinearElasticElement::CommonData &elastic_common_data = iso_elastic.commonData;
        NonlinearElasticElement::BlockData &elastic_block_data = mit->second;
        nitsche_element_iso.getOpPtrVector().push_back(
          new NonlinearElasticElement::OpGetCommonDataAtGaussPts(
            "DISPLACEMENT",elastic_common_data
          )
        );
        nitsche_element_iso.getOpPtrVector().push_back(
          new NonlinearElasticElement::OpGetCommonDataAtGaussPts(
            "MESH_NODE_POSITIONS",elastic_common_data
          )
        );
        nitsche_element_iso.getOpPtrVector().push_back(
          new NonlinearElasticElement::OpJacobianPiolaKirchhoffStress(
            "DISPLACEMENT",elastic_block_data,elastic_common_data,1,false,false,true
          )
        );
        nitsche_element_iso.getOpPtrVector().push_back(
          new PeriodicNitscheConstrains::OpCalculateHomogenisedStressSurfaceIntegral(
            "DISPLACEMENT",nitsche_block_data,nitsche_common_data,
            elastic_block_data,elastic_common_data,stress_homo
          )
        );
      }
      for(
        map<int,NonlinearElasticElement::BlockData>::iterator mit = trans_elastic.setOfBlocks.begin();
        mit!=trans_elastic.setOfBlocks.end();
        mit++
      ) {
        NonlinearElasticElement::CommonData &elastic_common_data = trans_elastic.commonData;
        NonlinearElasticElement::BlockData &elastic_block_data = mit->second;
        nitsche_element_trans.getOpPtrVector().push_back(
          new NonlinearElasticElement::OpGetCommonDataAtGaussPts(
            "DISPLACEMENT",elastic_common_data
          )
        );
        nitsche_element_trans.getOpPtrVector().push_back(
          new NonlinearElasticElement::OpGetCommonDataAtGaussPts(
            "MESH_NODE_POSITIONS",elastic_common_data
          )
        );
        nitsche_element_trans.getOpPtrVector().push_back(
          new NonlinearElasticElement::OpGetCommonDataAtGaussPts(
            "POTENTIAL_FIELD",elastic_common_data
          )
        );
        nitsche_element_trans.getOpPtrVector().push_back(
          new NonlinearElasticElement::OpJacobianPiolaKirchhoffStress(
            "DISPLACEMENT",elastic_block_data,elastic_common_data,2,false,false,true
          )
        );
        nitsche_element_trans.getOpPtrVector().push_back(
          new PeriodicNitscheConstrains::OpCalculateHomogenisedStressSurfaceIntegral(
            "DISPLACEMENT",nitsche_block_data,nitsche_common_data,
            elastic_block_data,elastic_common_data,stress_homo
          )
        );
      }
    } else {
      for(
        map<int,NonlinearElasticElement::BlockData>::iterator mit = iso_elastic.setOfBlocks.begin();
        mit!=iso_elastic.setOfBlocks.end();
        mit++
      ) {
        NonlinearElasticElement::CommonData &elastic_common_data = iso_elastic.commonData;
        NonlinearElasticElement::BlockData &elastic_block_data = mit->second;
        ave_stress_iso.getOpPtrVector().push_back(
          new NonlinearElasticElement::OpGetCommonDataAtGaussPts(
            "DISPLACEMENT",elastic_common_data
          )
        );
        ave_stress_iso.getOpPtrVector().push_back(
          new NonlinearElasticElement::OpGetCommonDataAtGaussPts(
            "MESH_NODE_POSITIONS",elastic_common_data
          )
        );
        ave_stress_iso.getOpPtrVector().push_back(
          new NonlinearElasticElement::OpJacobianPiolaKirchhoffStress(
            "DISPLACEMENT",elastic_block_data,elastic_common_data,1,false,false,true
          )
        );
        ave_stress_iso.getOpPtrVector().push_back(
          new PeriodicNitscheConstrains::OpCalculateHomogenisedStressVolumeIntegral(
            "DISPLACEMENT",elastic_block_data,elastic_common_data,stress_homo
          )
        );
      }
      for(
        map<int,NonlinearElasticElement::BlockData>::iterator mit = trans_elastic.setOfBlocks.begin();
        mit!=trans_elastic.setOfBlocks.end();
        mit++
      ) {
        NonlinearElasticElement::CommonData &elastic_common_data = trans_elastic.commonData;
        NonlinearElasticElement::BlockData &elastic_block_data = mit->second;
        ave_stress_trans.getOpPtrVector().push_back(
          new NonlinearElasticElement::OpGetCommonDataAtGaussPts(
            "DISPLACEMENT",elastic_common_data
          )
        );
        ave_stress_trans.getOpPtrVector().push_back(
          new NonlinearElasticElement::OpGetCommonDataAtGaussPts(
            "MESH_NODE_POSITIONS",elastic_common_data
          )
        );
        ave_stress_trans.getOpPtrVector().push_back(
          new NonlinearElasticElement::OpGetCommonDataAtGaussPts(
            "POTENTIAL_FIELD",elastic_common_data
          )
        );
        ave_stress_trans.getOpPtrVector().push_back(
          new NonlinearElasticElement::OpJacobianPiolaKirchhoffStress(
            "DISPLACEMENT",elastic_block_data,elastic_common_data,2,false,false,true
          )
        );
        ave_stress_trans.getOpPtrVector().push_back(
          new PeriodicNitscheConstrains::OpCalculateHomogenisedStressVolumeIntegral(
            "DISPLACEMENT",elastic_block_data,elastic_common_data,stress_homo
          )
        );
      }
    }
    break;
    case NOHOMOBC:
    SETERRQ(PETSC_COMM_SELF,MOFEM_IMPOSIBLE_CASE,"*** Boundary conditions not set");
  }

  //----------------------------------------------------------------------------
  // 3.2 Calculating zeroth-order homogenized stress using volume averaging theorem
  //----------------------------------------------------------------------------
  struct MyPostProc: public PostProcVolumeOnRefinedMesh {

    bool doPreProcess;
    bool doPostProcess;

    MyPostProc(FieldInterface &m_field):
    PostProcVolumeOnRefinedMesh(m_field),
    doPreProcess(true),
    doPostProcess(true)
    {}

    void setDoPreProcess() { doPreProcess = true; }
    void unSetDoPreProcess() { doPreProcess = false; }
    void setDoPostProcess() { doPostProcess = true; }
    void unSetDoPostProcess() { doPostProcess = false; }

    PetscErrorCode ierr;

    PetscErrorCode preProcess() {
      PetscFunctionBegin;
      if(doPreProcess) {
        ierr = PostProcVolumeOnRefinedMesh::preProcess(); CHKERRQ(ierr);
      }
      PetscFunctionReturn(0);
    }
    PetscErrorCode postProcess() {
      PetscFunctionBegin;
      if(doPostProcess) {
        ierr = PostProcVolumeOnRefinedMesh::postProcess(); CHKERRQ(ierr);
      }
      PetscFunctionReturn(0);
    }
  };

  MyPostProc post_proc(m_field);

  ierr = post_proc.generateReferenceElementMesh(); CHKERRQ(ierr);
  ierr = post_proc.addFieldValuesPostProc("DISPLACEMENT"); CHKERRQ(ierr);
  ierr = post_proc.addFieldValuesGradientPostProc("DISPLACEMENT"); CHKERRQ(ierr);
  ierr = post_proc.addFieldValuesPostProc("MESH_NODE_POSITIONS"); CHKERRQ(ierr);
  ierr = post_proc.addFieldValuesGradientPostProc("MESH_NODE_POSITIONS"); CHKERRQ(ierr);
  if(trans_iso_blocks) {
    ierr = post_proc.addFieldValuesGradientPostProc("POTENTIAL_FIELD"); CHKERRQ(ierr);
  }
  for(
    map<int,NonlinearElasticElement::BlockData>::iterator sit = iso_elastic.setOfBlocks.begin();
    sit != iso_elastic.setOfBlocks.end(); sit++
  ) {
    post_proc.getOpPtrVector().push_back(
      new PostPorcStress(
        post_proc.postProcMesh,
        post_proc.mapGaussPts,
        "DISPLACEMENT",
        sit->second,
        post_proc.commonData,
        false
      )
    );
  }
  for(
    map<int,NonlinearElasticElement::BlockData>::iterator sit = trans_elastic.setOfBlocks.begin();
    sit != trans_elastic.setOfBlocks.end(); sit++
  ) {
    post_proc.getOpPtrVector().push_back(
      new PostPorcStress(
        post_proc.postProcMesh,
        post_proc.mapGaussPts,
        "DISPLACEMENT",
        sit->second,
        post_proc.commonData
      )
    );
  }

  PetscScalar *avec;
  ierr = VecGetArray(stress_homo,&avec); CHKERRQ(ierr);
  for(int ii = 0;ii<6;ii++) {
    ierr = VecZeroEntries(D); CHKERRQ(ierr);
    ierr = KSPSolve(solver,F[ii],D); CHKERRQ(ierr);
    ierr = VecGhostUpdateBegin(D,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
    ierr = VecGhostUpdateEnd(D,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
    ierr = m_field.set_global_ghost_vector(
      "ELASTIC_MECHANICS",ROW,D,INSERT_VALUES,SCATTER_REVERSE
    ); CHKERRQ(ierr);
    post_proc.setDoPreProcess();
    post_proc.unSetDoPostProcess();
    ierr = m_field.loop_finite_elements(
      "ELASTIC_MECHANICS","ELASTIC",post_proc
    ); CHKERRQ(ierr);
    post_proc.unSetDoPreProcess();
    post_proc.setDoPostProcess();
    ierr = m_field.loop_finite_elements(
      "ELASTIC_MECHANICS","TRAN_ISOTROPIC_ELASTIC",post_proc
    ); CHKERRQ(ierr);
    {
      ostringstream sss;
      sss << "mode_" << homo_bc_names[choise_value] << "_" << ii << ".h5m";
      rval = post_proc.postProcMesh.write_file(
        sss.str().c_str(),"MOAB","PARALLEL=WRITE_PART"
      ); CHKERR_PETSC(rval);
    }
    ierr = VecZeroEntries(stress_homo); CHKERRQ(ierr);
    if(choise_value == HOMOBCDISP) {
      ierr = m_field.loop_finite_elements(
        "ELASTIC_MECHANICS","LAGRANGE_ELEM",lagrangian_element_disp.getLoopFeRVEBCStress()
      ); CHKERRQ(ierr);
    }
    if(choise_value == HOMOBCTRAC) {
      ierr = m_field.loop_finite_elements(
        "ELASTIC_MECHANICS","LAGRANGE_ELEM",lagrangian_element_trac.getLoopFeRVEBCStress()
      ); CHKERRQ(ierr);
    }
    if(choise_value == HOMOBCPERIODIC) {
      ierr = m_field.loop_finite_elements(
        "ELASTIC_MECHANICS","LAGRANGE_ELEM",lagrangian_element_periodic.getLoopFeRVEBCStress()
      ); CHKERRQ(ierr);
    }
    if(choise_value == NITSCHE) {
      if(stress_by_boundary_integral) {
        nitsche_element_iso.addToRule = 1;
        nitsche_element_trans.addToRule = 1;
        ierr = m_field.loop_finite_elements(
          "ELASTIC_MECHANICS","ELASTIC",nitsche_element_iso
        );  CHKERRQ(ierr);
        ierr = m_field.loop_finite_elements(
          "ELASTIC_MECHANICS","TRAN_ISOTROPIC_ELASTIC",nitsche_element_trans
        );  CHKERRQ(ierr);
      } else {
        ave_stress_iso.addToRule = 1;
        ave_stress_trans.addToRule = 1;
        ierr = m_field.loop_finite_elements(
          "ELASTIC_MECHANICS","ELASTIC",ave_stress_iso
        );  CHKERRQ(ierr);
        ierr = m_field.loop_finite_elements(
          "ELASTIC_MECHANICS","TRAN_ISOTROPIC_ELASTIC",ave_stress_trans
        );  CHKERRQ(ierr);
      }
    }
    ierr = PetscOptionsGetReal(
      PETSC_NULL,"-my_rve_volume",&rve_volume,PETSC_NULL
    ); CHKERRQ(ierr);
    ierr = VecAssemblyBegin(stress_homo); CHKERRQ(ierr);
    ierr = VecAssemblyEnd(stress_homo); CHKERRQ(ierr);
    ierr = VecScale(stress_homo,1.0/rve_volume); CHKERRQ(ierr);
    ierr = VecGhostUpdateBegin(stress_homo,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
    ierr = VecGhostUpdateEnd(stress_homo,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
    for(int jj=0; jj<6; jj++) {
      Dmat(jj,ii) = avec[jj];
    }
  }
  ierr = VecRestoreArray(stress_homo,&avec); CHKERRQ(ierr);

  
  // ===========================================================================
  //
  //  VII. OUTPUT
  //
  // ===========================================================================
  PetscPrintf(
    PETSC_COMM_WORLD,"\nHomogenised Stiffens Matrix = \n\n"
  );

  for(int ii=0; ii<6; ii++) {
    PetscPrintf(
      PETSC_COMM_WORLD,
      "stress %d\t\t%8.5e\t\t%8.5e\t\t%8.5e\t\t%8.5e\t\t%8.5e\t\t%8.5e\n",
      ii,Dmat(ii,0),Dmat(ii,1),Dmat(ii,2),Dmat(ii,3),Dmat(ii,4),Dmat(ii,5)
    );
  }

  
  // ===========================================================================
  //
  //  VIII. FINISH
  //
  // ===========================================================================
  //detroy matrices
  for(int ii = 0;ii<6;ii++) {
    ierr = VecDestroy(&F[ii]); CHKERRQ(ierr);
  }
  ierr = VecDestroy(&D); CHKERRQ(ierr);
  ierr = MatDestroy(&Aij); CHKERRQ(ierr);
  ierr = KSPDestroy(&solver); CHKERRQ(ierr);
  ierr = VecDestroy(&stress_homo); CHKERRQ(ierr);

  ierr = PetscTime(&v2);CHKERRQ(ierr);
  ierr = PetscGetCPUTime(&t2);CHKERRQ(ierr);

  PetscSynchronizedPrintf(PETSC_COMM_WORLD,"Total Rank %d Time = %f CPU Time = %f\n",pcomm->rank(),v2-v1,t2-t1);

  PetscFinalize();
}
