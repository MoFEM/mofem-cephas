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

#include "FieldInterface.hpp"
#include "FieldCore.hpp"
#include "FEMethod_UpLevelStudent.hpp"
#include "cholesky.hpp"
#include <petscksp.h>

#include "ElasticFEMethodTransIso.hpp"
#include "PostProcVertexMethod.hpp"
#include "PostProcDisplacementAndStrainOnRefindedMesh.hpp"

#include "ElasticFE_RVELagrange_Periodic.hpp"
#include "ElasticFE_RVELagrange_RigidBodyTranslation.hpp"
#include "ElasticFE_RVELagrange_RigidBodyRotation.hpp"
#include "ElasticFE_RVELagrange_Homogenized_Stress_Periodic.hpp"
#include "RVEVolume.hpp"

#include "ElasticFEMethodInterface.hpp"

using namespace MoFEM;

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
    
    Core mb_instance;
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
    
    PetscInt order;
    ierr = PetscOptionsGetInt(PETSC_NULL,"-my_order",&order,&flg); CHKERRQ(ierr);
    if(flg != PETSC_TRUE) {
        order = 1;
    }
    
    char outName[PETSC_MAX_PATH_LEN]="out.vtk";
    ierr = PetscOptionsGetString(PETSC_NULL,"-my_out",outName,sizeof(outName),&flg); CHKERRQ(ierr);
    
    char outName2[PETSC_MAX_PATH_LEN]="out_post_proc.vtk";
    ierr = PetscOptionsGetString(PETSC_NULL,"-my_post_out",outName2,sizeof(outName2),&flg); CHKERRQ(ierr);
    
    //Applied strain on the RVE (vector of length 6) strain=[xx, yy, zz, xy, xz, zy]^T
    double myapplied_strain[6];
    int nmax=6;
    ierr = PetscOptionsGetRealArray(PETSC_NULL,"-myapplied_strain",myapplied_strain,&nmax,&flg); CHKERRQ(ierr);
    ublas::vector<FieldData> applied_strain;
    applied_strain.resize(6);
    cblas_dcopy(6, &myapplied_strain[0], 1, &applied_strain(0), 1);
    //    cout<<"applied_strain ="<<applied_strain<<endl;

    
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
    FieldCore core(moab);
    FieldInterface& mField = core;
    
	
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
    ierr = mField.build_finiteElementsPtr(); CHKERRQ(ierr);
    ierr = mField.build_adjacencies(bit_levels.back()); CHKERRQ(ierr);
    ierr = mField.build_problems(); CHKERRQ(ierr);
    
    //    std::cout << float( clock () - begin_time ) /  CLOCKS_PER_SEC<<endl;
	
	//Build FE
    
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
		
		std::size_t found=it->get_Cubit_name().find("PotentialFlow");
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
		for(_IT_GET_FES_BY_NAME_FOR_LOOP_(mField, sss.str().c_str() ,it)){
			RangeFibre[ii].insert(it->get_ent());
			rval = moab.create_meshset(MESHSET_SET,fibre_meshset[ii]); CHKERR_PETSC(rval);
			rval = moab.add_entities(fibre_meshset[ii],RangeFibre[ii]); CHKERR_PETSC(rval);
			rval = moab.unite_meshset(meshset_Trans_ISO,fibre_meshset[ii]); CHKERR_PETSC(rval);
		}
	}
    
    rval = moab.write_file("meshset_Trans_ISO.vtk","VTK","",&meshset_Trans_ISO,1); CHKERR_PETSC(rval);
    
	for(_IT_CUBITMESHSETS_BY_SET_TYPE_FOR_LOOP_(mField,BLOCKSET,it)){
        
		if(it->get_Cubit_name() == "MAT_ELASTIC_1") {
			Range TetsInBlock;
			rval = moab.get_entities_by_type(it->meshset, MBTET,TetsInBlock,true); CHKERR_PETSC(rval);
			Range block_rope_bit_level = intersect(LatestRefinedTets,TetsInBlock);
			rval = moab.add_entities(meshset_Elastic,block_rope_bit_level);CHKERR_PETSC(rval);
		}
	}
    
	Range prims_on_problem_bit_level;
	ierr = mField.get_entities_by_type_and_ref_level(problem_bit_level,BitRefLevel().set(),MBPRISM,prims_on_problem_bit_level); CHKERRQ(ierr);
  //to create meshset from range
  EntityHandle meshset_prims_on_problem_bit_level;
  rval = moab.create_meshset(MESHSET_SET,meshset_prims_on_problem_bit_level); CHKERR_PETSC(rval);
	rval = moab.add_entities(meshset_prims_on_problem_bit_level,prims_on_problem_bit_level); CHKERR_PETSC(rval);
  ierr = mField.seed_ref_level_MESHSET(meshset_prims_on_problem_bit_level,BitRefLevel().set()); CHKERRQ(ierr);
	
    /***/
    //Define problem
    
    //Fields
    ierr = mField.add_field("DISPLACEMENT",H1,3); CHKERRQ(ierr);
    ierr = mField.add_field("Lagrange_mul_disp",H1,3); CHKERRQ(ierr);  //For lagrange multipliers to control the periodic motion
    ierr = mField.add_field("Lagrange_mul_disp_rigid_trans",NOFIELD,3); CHKERRQ(ierr);  //To control the rigid body motion (3 Traslations and 3 rotations)
    ierr = mField.add_field("Lagrange_mul_disp_rigid_rotation",NOFIELD,3); CHKERRQ(ierr); //Controla 3 rigid body rotations about x, y and z axis

    
    //FE
    ierr = mField.add_finite_element("ELASTIC"); CHKERRQ(ierr);
    ierr = mField.add_finite_element("TRAN_ISOTROPIC_ELASTIC"); CHKERRQ(ierr);
    ierr = mField.add_finite_element("INTERFACE"); CHKERRQ(ierr);
    ierr = mField.add_finite_element("Lagrange_elem"); CHKERRQ(ierr);
    ierr = mField.add_finite_element("Lagrange_elem_rigid_trans"); CHKERRQ(ierr);//For rigid body control
    ierr = mField.add_finite_element("Lagrange_elem_rigid_rotation"); CHKERRQ(ierr);

    //Define rows/cols and element data
    ierr = mField.modify_finite_element_add_field_row("ELASTIC","DISPLACEMENT"); CHKERRQ(ierr);
    ierr = mField.modify_finite_element_add_field_col("ELASTIC","DISPLACEMENT"); CHKERRQ(ierr);
    ierr = mField.modify_finite_element_add_field_data("ELASTIC","DISPLACEMENT"); CHKERRQ(ierr);
    //FE Transverse Isotropic
    ierr = mField.modify_finite_element_add_field_row("TRAN_ISOTROPIC_ELASTIC","DISPLACEMENT"); CHKERRQ(ierr);
    ierr = mField.modify_finite_element_add_field_col("TRAN_ISOTROPIC_ELASTIC","DISPLACEMENT"); CHKERRQ(ierr);
    ierr = mField.modify_finite_element_add_field_data("TRAN_ISOTROPIC_ELASTIC","DISPLACEMENT"); CHKERRQ(ierr);
    //FE Interface
    ierr = mField.modify_finite_element_add_field_row("INTERFACE","DISPLACEMENT"); CHKERRQ(ierr);
    ierr = mField.modify_finite_element_add_field_col("INTERFACE","DISPLACEMENT"); CHKERRQ(ierr);
    ierr = mField.modify_finite_element_add_field_data("INTERFACE","DISPLACEMENT"); CHKERRQ(ierr);
    ierr = mField.modify_finite_element_add_field_data("TRAN_ISOTROPIC_ELASTIC","POTENTIAL_FIELD"); CHKERRQ(ierr);
    
    
    //Define rows/cols and element data for C and CT (for lagrange multipliers)
    //============================================================================================================
    //C row as Lagrange_mul_disp and col as DISPLACEMENT
    ierr = mField.modify_finite_element_add_field_row("Lagrange_elem","Lagrange_mul_disp"); CHKERRQ(ierr);
    ierr = mField.modify_finite_element_add_field_col("Lagrange_elem","DISPLACEMENT"); CHKERRQ(ierr);
    
    //CT col as Lagrange_mul_disp and row as DISPLACEMENT
    ierr = mField.modify_finite_element_add_field_col("Lagrange_elem","Lagrange_mul_disp"); CHKERRQ(ierr);
    ierr = mField.modify_finite_element_add_field_row("Lagrange_elem","DISPLACEMENT"); CHKERRQ(ierr);
    
    //As for stress we need both displacement and temprature (Lukasz)
    ierr = mField.modify_finite_element_add_field_data("Lagrange_elem","Lagrange_mul_disp"); CHKERRQ(ierr);
    ierr = mField.modify_finite_element_add_field_data("Lagrange_elem","DISPLACEMENT"); CHKERRQ(ierr);
    //============================================================================================================
    
    
    //Define rows/cols and element data for C1 and C1T (for lagrange multipliers to contol the rigid body motion)
    //============================================================================================================
    //C1 row as Lagrange_mul_disp and col as DISPLACEMENT
    ierr = mField.modify_finite_element_add_field_row("Lagrange_elem_rigid_trans","Lagrange_mul_disp_rigid_trans"); CHKERRQ(ierr);
    ierr = mField.modify_finite_element_add_field_col("Lagrange_elem_rigid_trans","DISPLACEMENT"); CHKERRQ(ierr);
    
    //C1T col as Lagrange_mul_disp and row as DISPLACEMENT
    ierr = mField.modify_finite_element_add_field_col("Lagrange_elem_rigid_trans","Lagrange_mul_disp_rigid_trans"); CHKERRQ(ierr);
    ierr = mField.modify_finite_element_add_field_row("Lagrange_elem_rigid_trans","DISPLACEMENT"); CHKERRQ(ierr);
    
    //As for stress we need both displacement and temprature (Lukasz)
    ierr = mField.modify_finite_element_add_field_data("Lagrange_elem_rigid_trans","Lagrange_mul_disp_rigid_trans"); CHKERRQ(ierr);
    ierr = mField.modify_finite_element_add_field_data("Lagrange_elem_rigid_trans","DISPLACEMENT"); CHKERRQ(ierr);
    //============================================================================================================

    
    //Define rows/cols and element data for C2 and C2T (for lagrange multipliers to contol the rigid body rotations)
    //============================================================================================================
    //C2 row as Lagrange_elem_rigid_trans and col as DISPLACEMENT
    ierr = mField.modify_finite_element_add_field_row("Lagrange_elem_rigid_rotation","Lagrange_mul_disp_rigid_rotation"); CHKERRQ(ierr);
    ierr = mField.modify_finite_element_add_field_col("Lagrange_elem_rigid_rotation","DISPLACEMENT"); CHKERRQ(ierr);
    
    //C2T col as Lagrange_elem_rigid_trans and row as DISPLACEMENT
    ierr = mField.modify_finite_element_add_field_col("Lagrange_elem_rigid_rotation","Lagrange_mul_disp_rigid_rotation"); CHKERRQ(ierr);
    ierr = mField.modify_finite_element_add_field_row("Lagrange_elem_rigid_rotation","DISPLACEMENT"); CHKERRQ(ierr);
    
    //As for stress we need both displacement and temprature (Lukasz)
    ierr = mField.modify_finite_element_add_field_data("Lagrange_elem_rigid_rotation","Lagrange_mul_disp_rigid_rotation"); CHKERRQ(ierr);
    ierr = mField.modify_finite_element_add_field_data("Lagrange_elem_rigid_rotation","DISPLACEMENT"); CHKERRQ(ierr);
    //============================================================================================================

     //define problems
    ierr = mField.add_problem("ELASTIC_MECHANICS"); CHKERRQ(ierr);
    
    //set finite elements for problem
    ierr = mField.modify_problem_add_finite_element("ELASTIC_MECHANICS","ELASTIC"); CHKERRQ(ierr);
    ierr = mField.modify_problem_add_finite_element("ELASTIC_MECHANICS","INTERFACE"); CHKERRQ(ierr);
    ierr = mField.modify_problem_add_finite_element("ELASTIC_MECHANICS","TRAN_ISOTROPIC_ELASTIC"); CHKERRQ(ierr);
    ierr = mField.modify_problem_add_finite_element("ELASTIC_MECHANICS","Lagrange_elem"); CHKERRQ(ierr);
    ierr = mField.modify_problem_add_finite_element("ELASTIC_MECHANICS","Lagrange_elem_rigid_trans"); CHKERRQ(ierr);
    ierr = mField.modify_problem_add_finite_element("ELASTIC_MECHANICS","Lagrange_elem_rigid_rotation"); CHKERRQ(ierr);

    //set refinment level for problem
    ierr = mField.modify_problem_ref_level_add_bit("ELASTIC_MECHANICS",problem_bit_level); CHKERRQ(ierr);
    
    /***/
    //Declare problem
    
    //add entitities (by tets) to the field
    ierr = mField.add_ents_to_field_by_TETs(0,"DISPLACEMENT"); CHKERRQ(ierr);
    
    //add finite elements entities
    ierr = mField.add_ents_to_finite_element_by_TETs(meshset_Elastic,"ELASTIC",true); CHKERRQ(ierr);
    ierr = mField.add_ents_to_finite_element_by_TETs(meshset_Trans_ISO,"TRAN_ISOTROPIC_ELASTIC",true); CHKERRQ(ierr);
    
    
    ierr = mField.add_ents_to_finite_element_EntType_by_bit_ref(problem_bit_level,"INTERFACE",MBPRISM); CHKERRQ(ierr);
    
    
    //Add finite element to lagrange element for rigid body translation
    Range Tris_NewWholeMesh, Tri_OldNewSurf, SurfacesFaces;
    ierr = mField.get_entities_by_type_and_ref_level(bit_levels.back(),BitRefLevel().set(),MBTRI,Tris_NewWholeMesh); CHKERRQ(ierr);
    ierr = mField.get_Cubit_msId_entities_by_dimension(103,SIDESET,2,Tri_OldNewSurf,true); CHKERRQ(ierr);
    SurfacesFaces = intersect(Tris_NewWholeMesh,Tri_OldNewSurf);
    ierr = PetscPrintf(PETSC_COMM_WORLD,"number of SIDESET 103 = %d\n",SurfacesFaces.size()); CHKERRQ(ierr);
  
  //to create meshset from range
  EntityHandle BoundFacesMeshset;
  rval = moab.create_meshset(MESHSET_SET,BoundFacesMeshset); CHKERR_PETSC(rval);
	rval = moab.add_entities(BoundFacesMeshset,SurfacesFaces); CHKERR_PETSC(rval);
  ierr = mField.seed_ref_level_MESHSET(BoundFacesMeshset,BitRefLevel().set()); CHKERRQ(ierr);

  ierr = mField.add_ents_to_finite_element_by_TRIs(SurfacesFaces,"Lagrange_elem_rigid_trans"); CHKERRQ(ierr);
    ierr = mField.add_ents_to_finite_element_by_TRIs(SurfacesFaces,"Lagrange_elem_rigid_rotation"); CHKERRQ(ierr);

    
    //=======================================================================================================
    //Add Periodic Prisims Between Triangles on -ve and +ve faces to implement periodic bounary conditions
    //=======================================================================================================

    //Populating the Multi-index container with -ve triangles
    Range Tri_OldNewSurfNeg, SurTrisNeg;
    ierr = mField.get_Cubit_msId_entities_by_dimension(101,SIDESET,2,Tri_OldNewSurfNeg,true); CHKERRQ(ierr);
    SurTrisNeg = intersect(Tris_NewWholeMesh,Tri_OldNewSurfNeg);
    ierr = PetscPrintf(PETSC_COMM_WORLD,"number of SIDESET 101 = %d\n",SurTrisNeg.size()); CHKERRQ(ierr);

    Face_CenPos_Handle_multiIndex Face_CenPos_Handle_varNeg, Face_CenPos_Handle_varPos;
    double TriCen[3], coords_Tri[9];
    
    double roundfact=1000.0;
    int count1=1;
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
//        if(count1>2264){
//            cout<<"   TriCen[0]= "<<TriCen[0] << "   TriCen[1]= "<< TriCen[1] << "   TriCen[2]= "<< TriCen[2] <<endl;
//        }
        count1++;
        //fill the multi-index container with centriod coordinates and triangle handles
        Face_CenPos_Handle_varNeg.insert(Face_CenPos_Handle(TriCen[0], TriCen[1], TriCen[2], *it));
//        for(int ii=0; ii<3; ii++) cout<<"TriCen "<<TriCen[ii]<<endl;
//        cout<<endl<<endl;
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
    ierr = mField.get_Cubit_msId_entities_by_dimension(102,SIDESET,2,Tri_OldNewSurfPos,true); CHKERRQ(ierr);
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
//    typedef Face_CenPos_Handle_multiIndex::index<xcoord_tag>::type::iterator Tri_Xcoord_iterator;
//    typedef Face_CenPos_Handle_multiIndex::index<ycoord_tag>::type::iterator Tri_Ycoord_iterator;
//    typedef Face_CenPos_Handle_multiIndex::index<zcoord_tag>::type::iterator Tri_Zcoord_iterator;
//    Tri_Xcoord_iterator XcoordMin_it, XcoordMax_it;
//    Tri_Ycoord_iterator YcoordMin_it, YcoordMax_it;
//    Tri_Zcoord_iterator ZcoordMin_it, ZcoordMax_it;
//    
//    //XcoordMax_it-- because .end() will point iterator after the data range but .begin() will point the iteratore to the first value of range
//    XcoordMin_it=Face_CenPos_Handle_varNeg.get<xcoord_tag>().begin();                  XcoordMin=XcoordMin_it->xcoord;
//    XcoordMax_it=Face_CenPos_Handle_varPos.get<xcoord_tag>().end();    XcoordMax_it--; XcoordMax=XcoordMax_it->xcoord;
//    YcoordMin_it=Face_CenPos_Handle_varNeg.get<ycoord_tag>().begin();                  YcoordMin=YcoordMin_it->ycoord;
//    YcoordMax_it=Face_CenPos_Handle_varPos.get<ycoord_tag>().end();    YcoordMax_it--; YcoordMax=YcoordMax_it->ycoord;
//    ZcoordMin_it=Face_CenPos_Handle_varNeg.get<zcoord_tag>().begin();                  ZcoordMin=ZcoordMin_it->zcoord;
//    ZcoordMax_it=Face_CenPos_Handle_varPos.get<zcoord_tag>().end();    ZcoordMax_it--; ZcoordMax=ZcoordMax_it->zcoord;
//
//    cout<<"XcoordMin "<<XcoordMin << "      XcoordMax "<<XcoordMax <<endl;
//    cout<<"YcoordMin "<<YcoordMin << "      YcoordMax "<<YcoordMax <<endl;
//    cout<<"ZcoordMin "<<ZcoordMin << "      ZcoordMax "<<ZcoordMax <<endl;
    
    XcoordMin =-1.5;      XcoordMax =1.5;
    YcoordMin =-0.15;      YcoordMax =0.15;
    ZcoordMin =-0.195;      ZcoordMax =0.585;

    //Creating Prisims between triangles on -ve and +ve faces
    typedef Face_CenPos_Handle_multiIndex::index<Tri_Hand_tag>::type::iterator Tri_Hand_iterator;
    Tri_Hand_iterator Tri_Neg;
    typedef Face_CenPos_Handle_multiIndex::index<Composite_xyzcoord>::type::iterator xyzcoord_iterator;
    xyzcoord_iterator Tri_Pos;
    Range PrismRange;
    double XPos, YPos, ZPos;
    int count=0;
    
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
                //Round nodal coordinates to 3 dicimal places only for comparison
                //round values to 3 disimal places
                if(XNodeNeg>=0) XNodeNeg=double(int(XNodeNeg*roundfact+0.5))/roundfact;   else XNodeNeg=double(int(XNodeNeg*roundfact-0.5))/roundfact;
                if(YNodeNeg>=0) YNodeNeg=double(int(YNodeNeg*roundfact+0.5))/roundfact;   else YNodeNeg=double(int(YNodeNeg*roundfact-0.5))/roundfact;
                if(ZNodeNeg>=0) ZNodeNeg=double(int(ZNodeNeg*roundfact+0.5))/roundfact;   else ZNodeNeg=double(int(ZNodeNeg*roundfact-0.5))/roundfact;
                
                XNodePos=CoordNodePos[3*jj]; YNodePos=CoordNodePos[3*jj+1]; ZNodePos=CoordNodePos[3*jj+2];
                if(XNodePos>=0) XNodePos=double(int(XNodePos*roundfact+0.5))/roundfact;   else XNodePos=double(int(XNodePos*roundfact-0.5))/roundfact;
                if(YNodePos>=0) YNodePos=double(int(YNodePos*roundfact+0.5))/roundfact;   else YNodePos=double(int(YNodePos*roundfact-0.5))/roundfact;
                if(ZNodePos>=0) ZNodePos=double(int(ZNodePos*roundfact+0.5))/roundfact;   else ZNodePos=double(int(ZNodePos*roundfact-0.5))/roundfact;
                
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

    //set app. order
    //see Hierarchic Finite Element Bases on Unstructured Tetrahedral Meshes (Mark Ainsworth & Joe Coyle)
    ierr = mField.set_field_order(0,MBTET,"DISPLACEMENT",order); CHKERRQ(ierr);
    ierr = mField.set_field_order(0,MBTRI,"DISPLACEMENT",order); CHKERRQ(ierr);
    ierr = mField.set_field_order(0,MBEDGE,"DISPLACEMENT",order); CHKERRQ(ierr);
    ierr = mField.set_field_order(0,MBVERTEX,"DISPLACEMENT",1); CHKERRQ(ierr);
    
    ierr = mField.set_field_order(0,MBTRI,"Lagrange_mul_disp",order); CHKERRQ(ierr);
    ierr = mField.set_field_order(0,MBEDGE,"Lagrange_mul_disp",order); CHKERRQ(ierr);
    ierr = mField.set_field_order(0,MBVERTEX,"Lagrange_mul_disp",1); CHKERRQ(ierr);

    /****/
    //build database
    
    //build field
    ierr = mField.build_fields(); CHKERRQ(ierr);
    
    //build finite elemnts
    ierr = mField.build_finiteElementsPtr(); CHKERRQ(ierr);
    
    //build adjacencies
    ierr = mField.build_adjacencies(problem_bit_level); CHKERRQ(ierr);
    
    //build problem
    ierr = mField.build_problems(); CHKERRQ(ierr);
    
    /****/
    //mesh partitioning
    
    //partition
    ierr = mField.partition_problem("ELASTIC_MECHANICS"); CHKERRQ(ierr);
    ierr = mField.partition_finiteElementsPtr("ELASTIC_MECHANICS"); CHKERRQ(ierr);
    //what are ghost nodes, see Petsc Manual
    ierr = mField.partition_ghost_dofs("ELASTIC_MECHANICS"); CHKERRQ(ierr);
    
    //create matrices
    Vec F,D;
    ierr = mField.VecCreateGhost("ELASTIC_MECHANICS",ROW,&F); CHKERRQ(ierr);
    ierr = mField.VecCreateGhost("ELASTIC_MECHANICS",COL,&D); CHKERRQ(ierr);
    
    Mat Aij;
    ierr = mField.MatCreateMPIAIJWithArrays("ELASTIC_MECHANICS",&Aij); CHKERRQ(ierr);

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
    double YoungModulusP;
    double PoissonRatioP;
    double YoungModulusZ;
    double PoissonRatioPZ;
    double ShearModulusZP;
    double YoungModulus;
    double PoissonRatio;
    double alpha;
    
    for(_IT_CUBITMESHSETS_BY_SET_TYPE_FOR_LOOP_(mField,BLOCKSET,it))
    {
        cout << endl << *it << endl;
        
        //Get block name
        string name = it->get_Cubit_name();
        
        if (name.compare(0,20,"MAT_ELASTIC_TRANSISO") == 0)
        {
            Mat_Elastic_TransIso mydata;
            ierr = it->get_attribute_data_structure(mydata); CHKERRQ(ierr);
            cout << mydata;
            YoungModulusP=mydata.data.Youngp;
            YoungModulusZ=mydata.data.Youngz;
            PoissonRatioP=mydata.data.Poissonp;
            PoissonRatioPZ=mydata.data.Poissonpz;
            if (mydata.data.Shearzp!=0) {
                ShearModulusZP=mydata.data.Shearzp;
            }else{
                ShearModulusZP=YoungModulusZ/(2*(1+PoissonRatioPZ));}
        }
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
    
//    alpha = 50;
//    cout<<"alpha   = "<<alpha<<endl;
	
    InterfaceFEMethod IntMyFE(mField,Aij,D,F,YoungModulus*alpha);
    MyElasticFEMethod MyFE(mField,Aij,D,F,LAMBDA(YoungModulus,PoissonRatio),MU(YoungModulus,PoissonRatio));
    TranIsotropicFibreDirRotElasticFEMethod MyTIsotFE(mField,Aij,D,F);
    
    ElasticFE_RVELagrange_Periodic MyFE_RVELagrangePeriodic(mField,Aij,D,F,applied_strain);
    ElasticFE_RVELagrange_RigidBodyTranslation MyFE_RVELagrangeRigidBodyTrans(mField,Aij,D,F,applied_strain);
    ElasticFE_RVELagrange_RigidBodyRotation MyFE_RVELagrangeRigidBodyRotation(mField,Aij,D,F,applied_strain);

    ierr = VecZeroEntries(F); CHKERRQ(ierr);
    ierr = VecGhostUpdateBegin(F,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
    ierr = VecGhostUpdateEnd(F,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
    ierr = MatZeroEntries(Aij); CHKERRQ(ierr);
    
    ierr = mField.loop_finiteElementsPtr("ELASTIC_MECHANICS","ELASTIC",MyFE);  CHKERRQ(ierr);
    PetscSynchronizedFlush(PETSC_COMM_WORLD);
	ierr = mField.loop_finiteElementsPtr("ELASTIC_MECHANICS","TRAN_ISOTROPIC_ELASTIC",MyTIsotFE);  CHKERRQ(ierr);
	PetscSynchronizedFlush(PETSC_COMM_WORLD);
    ierr = mField.loop_finiteElementsPtr("ELASTIC_MECHANICS","INTERFACE",IntMyFE);  CHKERRQ(ierr);
    PetscSynchronizedFlush(PETSC_COMM_WORLD);
    ierr = mField.loop_finiteElementsPtr("ELASTIC_MECHANICS","Lagrange_elem",MyFE_RVELagrangePeriodic);  CHKERRQ(ierr);
    PetscSynchronizedFlush(PETSC_COMM_WORLD);
    ierr = mField.loop_finiteElementsPtr("ELASTIC_MECHANICS","Lagrange_elem_rigid_trans",MyFE_RVELagrangeRigidBodyTrans);  CHKERRQ(ierr);
    PetscSynchronizedFlush(PETSC_COMM_WORLD);
    ierr = mField.loop_finiteElementsPtr("ELASTIC_MECHANICS","Lagrange_elem_rigid_rotation",MyFE_RVELagrangeRigidBodyRotation);  CHKERRQ(ierr);
    PetscSynchronizedFlush(PETSC_COMM_WORLD);

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
    
    //Solver
    KSP solver;
    ierr = KSPCreate(PETSC_COMM_WORLD,&solver); CHKERRQ(ierr);
    ierr = KSPSetOperators(solver,Aij,Aij,SAME_NONZERO_PATTERN); CHKERRQ(ierr);
    ierr = KSPSetFromOptions(solver); CHKERRQ(ierr);
    ierr = KSPSetUp(solver); CHKERRQ(ierr);
    
    ierr = KSPSolve(solver,F,D); CHKERRQ(ierr);
    ierr = VecGhostUpdateBegin(D,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
    ierr = VecGhostUpdateEnd(D,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
    
    //Save data on mesh
    ierr = mField.set_global_VecCreateGhost("ELASTIC_MECHANICS",ROW,D,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
//    ierr = VecView(D,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);
    
    
    //Calculation of Homogenized stress
    //====================================================================================================================================
    double RVE_volume;    RVE_volume=0.0;  //RVE volume for full RVE We need this for stress calculation
    Vec RVE_volume_Vec;
    ierr = VecCreateMPI(PETSC_COMM_WORLD, 1, pcomm->size(), &RVE_volume_Vec);  CHKERRQ(ierr);
    ierr = VecZeroEntries(RVE_volume_Vec); CHKERRQ(ierr);
    
    RVEVolume MyRVEVol(mField,Aij,D,F,LAMBDA(YoungModulus,PoissonRatio),MU(YoungModulus,PoissonRatio), RVE_volume_Vec);
    RVEVolumeTrans MyRVEVolTrans(mField,Aij,D,F, RVE_volume_Vec);

    ierr = mField.loop_finiteElementsPtr("ELASTIC_MECHANICS","ELASTIC",MyRVEVol);  CHKERRQ(ierr);
    ierr = mField.loop_finiteElementsPtr("ELASTIC_MECHANICS","TRAN_ISOTROPIC_ELASTIC",MyRVEVolTrans);  CHKERRQ(ierr);
    
    ierr = VecView(RVE_volume_Vec,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);
    ierr = VecSum(RVE_volume_Vec, &RVE_volume);  CHKERRQ(ierr);
    cout<<"Final RVE_volume = "<< RVE_volume <<endl;
    cout<<"Actual RVE_volume = "<< 3*0.3*0.78<<endl;  //Lx=3, Ly=0.3; Lz=0.78
    
    //create a vector for 6 components of homogenized stress
    Vec Stress_Homo;
    ierr = VecCreateMPI(PETSC_COMM_WORLD, 6, 6*pcomm->size(), &Stress_Homo);  CHKERRQ(ierr);
    ierr = VecZeroEntries(Stress_Homo); CHKERRQ(ierr);

//    ierr = VecView(D,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);
    ElasticFE_RVELagrange_Homogenized_Stress_Periodic MyFE_RVEHomoStressPeriodic(mField,Aij,D,F,&RVE_volume,applied_strain, Stress_Homo);
    ierr = mField.loop_finiteElementsPtr("ELASTIC_MECHANICS","Lagrange_elem",MyFE_RVEHomoStressPeriodic);  CHKERRQ(ierr);

    if(pcomm->rank()) cout<< " Stress_Homo =  "<<endl;
    ierr = VecView(Stress_Homo,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);

    //====================================================================================================================================
 

     PostProcVertexMethod ent_method(moab);
    ierr = mField.loop_dofs("ELASTIC_MECHANICS","DISPLACEMENT",ROW,ent_method); CHKERRQ(ierr);
    
    if(pcomm->rank()==0) {
        EntityHandle out_meshset;
        rval = moab.create_meshset(MESHSET_SET,out_meshset); CHKERR_PETSC(rval);
        ierr = mField.problem_get_FE("ELASTIC_MECHANICS","ELASTIC",out_meshset); CHKERRQ(ierr);
        ierr = mField.problem_get_FE("ELASTIC_MECHANICS","TRAN_ISOTROPIC_ELASTIC",out_meshset); CHKERRQ(ierr);
        ierr = mField.problem_get_FE("ELASTIC_MECHANICS","INTERFACE",out_meshset); CHKERRQ(ierr);
        rval = moab.write_file(outName,"VTK","",&out_meshset,1); CHKERR_PETSC(rval);
        rval = moab.delete_entities(&out_meshset,1); CHKERR_PETSC(rval);
    }
    
//    TranIso_PostProc_FibreDirRot_OnRefMesh fe_post_proc_method( mField, LAMBDA(YoungModulusP,PoissonRatioP),MU(YoungModulusP,PoissonRatioP), YoungModulusP,YoungModulusZ,PoissonRatioP,PoissonRatioPZ,ShearModulusZP);
//
//    ierr = mField.loop_finiteElementsPtr("ELASTIC_MECHANICS","ELASTIC",fe_post_proc_method);  CHKERRQ(ierr);
//    ierr = mField.loop_finiteElementsPtr("ELASTIC_MECHANICS","TRAN_ISOTROPIC_ELASTIC",fe_post_proc_method);  CHKERRQ(ierr);
//    
//    PetscSynchronizedFlush(PETSC_COMM_WORLD);
//    if(pcomm->rank()==0) {
//        rval = fe_post_proc_method.moab_post_proc.write_file(outName2,"VTK",""); CHKERR_PETSC(rval);
//    }
//    
//    PostProcCohesiveForces fe_post_proc_prisms(mField,YoungModulus*alpha);
//    ierr = mField.loop_finiteElementsPtr("ELASTIC_MECHANICS","INTERFACE",fe_post_proc_prisms);  CHKERRQ(ierr);
//    PetscSynchronizedFlush(PETSC_COMM_WORLD);
//    if(pcomm->rank()==0) {
//        rval = fe_post_proc_prisms.moab_post_proc.write_file("out_post_proc_prisms.vtk","VTK",""); CHKERR_PETSC(rval);
//    }
    
    
    //detroy matrices
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