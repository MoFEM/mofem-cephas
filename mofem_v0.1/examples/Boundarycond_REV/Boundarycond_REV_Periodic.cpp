/* Copyright (C) 2013, Lukasz Kaczmarczyk (likask AT wp.pl)
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

#include "ElasticFEMethod.hpp"
#include "ElasticFE_RVELagrange.hpp"
#include "PostProcVertexMethod.hpp"
#include "PostProcDisplacementAndStrainOnRefindedMesh.hpp"

using namespace MoFEM;

ErrorCode rval;
PetscErrorCode ierr;
static char help[] = "...\n\n";


//=================================================================================================================================
//Define class and multindex container to store data for traiangles on the boundary of the RVE (it cannot be defined within main)
//=================================================================================================================================

    struct Face_CenPos_Handle
    {
        double xcoord;
        double ycoord, zcoord;
        EntityHandle  Tri_Hand;
        Face_CenPos_Handle(double _xcoord, double _ycoord,  double _zcoord,  EntityHandle _Tri_Hand):xcoord(_xcoord),
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
            tag<Tri_Hand_tag>, member<Face_CenPos_Handle,EntityHandle,&Face_CenPos_Handle::Tri_Hand> >,

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
        order = 5;
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
    FieldCore core(moab);
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
    
    //Fields
    ierr = mField.add_field("DISPLACEMENT",H1,3); CHKERRQ(ierr);
    ierr = mField.add_field("Lagrange_mul_disp",H1,3); CHKERRQ(ierr);
    
    //FE
    ierr = mField.add_finite_element("ELASTIC"); CHKERRQ(ierr);
    ierr = mField.add_finite_element("Lagrange_elem"); CHKERRQ(ierr);
    
    //Define rows/cols and element data
    ierr = mField.modify_finite_element_add_field_row("ELASTIC","DISPLACEMENT"); CHKERRQ(ierr);
    ierr = mField.modify_finite_element_add_field_col("ELASTIC","DISPLACEMENT"); CHKERRQ(ierr);
    ierr = mField.modify_finite_element_add_field_data("ELASTIC","DISPLACEMENT"); CHKERRQ(ierr);
    
    //C row as Lagrange_mul_disp and col as DISPLACEMENT
    ierr = mField.modify_finite_element_add_field_row("Lagrange_elem","Lagrange_mul_disp"); CHKERRQ(ierr);
    ierr = mField.modify_finite_element_add_field_col("Lagrange_elem","DISPLACEMENT"); CHKERRQ(ierr);
    
    //CT col as Lagrange_mul_disp and row as DISPLACEMENT
    ierr = mField.modify_finite_element_add_field_col("Lagrange_elem","Lagrange_mul_disp"); CHKERRQ(ierr);
    ierr = mField.modify_finite_element_add_field_row("Lagrange_elem","DISPLACEMENT"); CHKERRQ(ierr);
    
    //As for stress we need both displacement and temprature (Lukasz)
    ierr = mField.modify_finite_element_add_field_data("Lagrange_elem","Lagrange_mul_disp"); CHKERRQ(ierr);
    ierr = mField.modify_finite_element_add_field_data("Lagrange_elem","DISPLACEMENT"); CHKERRQ(ierr);
    
    
    //define problems
    ierr = mField.add_problem("ELASTIC_MECHANICS"); CHKERRQ(ierr);
    
    
    //set finite elements for problem
    ierr = mField.modify_problem_add_finite_element("ELASTIC_MECHANICS","ELASTIC"); CHKERRQ(ierr);
    ierr = mField.modify_problem_add_finite_element("ELASTIC_MECHANICS","Lagrange_elem"); CHKERRQ(ierr);
    
    
    //set refinment level for problem
    ierr = mField.modify_problem_ref_level_add_bit("ELASTIC_MECHANICS",bit_level0); CHKERRQ(ierr);
    
    /***/
    //Declare problem
    
    //add entitities (by tets) to the field
    ierr = mField.add_ents_to_field_by_TETs(0,"DISPLACEMENT",2); CHKERRQ(ierr);
    
    
    //add finite elements entities
    ierr = mField.add_ents_to_finite_element_EntType_by_bit_ref(bit_level0,"ELASTIC",MBTET); CHKERRQ(ierr);

    
    

    Range SurTris;
    ierr = mField.get_Cubit_msId_entities_by_dimension(101,SideSet,2,SurTris,true); CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_WORLD,"number of SideSet 101 = %d\n",SurTris.size()); CHKERRQ(ierr);
   
    map<EntityHandle,EntityHandle> bc_map[SurTris.size()/2];
    Face_CenPos_Handle_multiIndex Face_CenPos_Handle_var;

    
    double TriCen[3], coords_Tri[9];
    
    for(Range::iterator it = SurTris.begin(); it!=SurTris.end();  it++) {
        const EntityHandle* conn_face;  int num_nodes_Tri;
        
        //get nodes attached to the face
        rval = moab.get_connectivity(*it,conn_face,num_nodes_Tri,true); CHKERR_PETSC(rval);
        //get nodal coordinates
        rval = moab.get_coords(conn_face,num_nodes_Tri,coords_Tri); CHKERR_PETSC(rval);
        
        //Find triangle centriod
        TriCen[0]= (coords_Tri[0]+coords_Tri[3]+coords_Tri[6])/3.0;
        TriCen[1]= (coords_Tri[1]+coords_Tri[4]+coords_Tri[7])/3.0;
        TriCen[2]= (coords_Tri[2]+coords_Tri[5]+coords_Tri[8])/3.0;
        
        //fill the multi-index container with centriod coordinates and triangle handles
        Face_CenPos_Handle_var.insert(Face_CenPos_Handle(TriCen[0], TriCen[1], TriCen[2], *it));
//        for(int ii=0; ii<3; ii++) cout<<"TriCen "<<TriCen[ii]<<endl;
//        cout<<endl<<endl;
    }
    
    
    typedef Face_CenPos_Handle_multiIndex::index<xcoord_tag>::type::iterator xcoord_iterator;
    xcoord_iterator xiit_L, hi_xiit_L, xiit_R, hi_xiit_R;
    
    //iterators for Left face
    xiit_L = Face_CenPos_Handle_var.get<xcoord_tag>().lower_bound(0.0);
    hi_xiit_L = Face_CenPos_Handle_var.get<xcoord_tag>().upper_bound(0.0);
    
    //iterators for Right face
    xiit_R = Face_CenPos_Handle_var.get<xcoord_tag>().lower_bound(1.0);
    hi_xiit_R= Face_CenPos_Handle_var.get<xcoord_tag>().upper_bound(1.0);

//    int count=1;
//    for(;xiit_L!=hi_xiit_L;xiit_L++) {
//        cout<<"count  "<<count<<endl;  count++;
//    }
//
//    count=1;
//    for(;xiit_R!=hi_xiit_R;xiit_R++) {
//        cout<<"count  "<<count<<endl;  count++;
//    }

    typedef Face_CenPos_Handle_multiIndex::index<Composite_xyzcoord>::type::iterator xyzcoord_iterator;
    xyzcoord_iterator iit;
    
    
    for(;xiit_L!=hi_xiit_L;xiit_L++) {
        iit=Face_CenPos_Handle_var.get<Composite_xyzcoord>().find( boost::make_tuple(xiit_L->xcoord, xiit_L->ycoord, xiit_L->zcoord) );
//        cout<<"xcoord= "<<xiit_L->xcoord << "   ycoord= "<< xiit_L->ycoord << "   ycoord= "<< xiit_L->zcoord <<endl;
//        cout<<"xcoord= "<<iit->xcoord << "   ycoord= "<< iit->ycoord << "   ycoord= "<< iit->zcoord <<endl<<endl<<endl;
        rval = moab.add_adjacencies(xiit_L->Tri_Hand, &(iit->Tri_Hand), 1, true); CHKERR_PETSC(rval);
        
        
        
    }
   
    
    

    
    
    
    
    
    
    
//    for(int ii=0; ii<)
//    
//    
//    
//    ierr = mField.seed_finite_elements(SurfacesFaces); CHKERRQ(ierr);
//    ierr = mField.add_ents_to_finite_element_by_TRIs(SurfacesFaces,"Lagrange_elem"); CHKERRQ(ierr);
//    
//    
//    //to create meshset from range
//    EntityHandle BoundFacesMeshset;
//    rval = moab.create_meshset(MESHSET_SET,BoundFacesMeshset); CHKERR_PETSC(rval);
//	rval = moab.add_entities(BoundFacesMeshset,SurfacesFaces); CHKERR_PETSC(rval);
//    
//    ierr = mField.add_ents_to_field_by_TRIs(BoundFacesMeshset,"Lagrange_mul_disp",2); CHKERRQ(ierr);
//    
//    //set app. order
//    //see Hierarchic Finite Element Bases on Unstructured Tetrahedral Meshes (Mark Ainsworth & Joe Coyle)
//    //int order = 5;
//    ierr = mField.set_field_order(0,MBTET,"DISPLACEMENT",order); CHKERRQ(ierr);
//    ierr = mField.set_field_order(0,MBTRI,"DISPLACEMENT",order); CHKERRQ(ierr);
//    ierr = mField.set_field_order(0,MBEDGE,"DISPLACEMENT",order); CHKERRQ(ierr);
//    ierr = mField.set_field_order(0,MBVERTEX,"DISPLACEMENT",1); CHKERRQ(ierr);
//    
//    ierr = mField.set_field_order(0,MBTRI,"Lagrange_mul_disp",order); CHKERRQ(ierr);
//    ierr = mField.set_field_order(0,MBEDGE,"Lagrange_mul_disp",order); CHKERRQ(ierr);
//    ierr = mField.set_field_order(0,MBVERTEX,"Lagrange_mul_disp",1); CHKERRQ(ierr);
//    
//    /****/
//    //build database
//    
//    //build field
//    ierr = mField.build_fields(); CHKERRQ(ierr);
//    
//    //build finite elemnts
//    ierr = mField.build_finite_elements(); CHKERRQ(ierr);
//    
//    //build adjacencies
//    ierr = mField.build_adjacencies(bit_level0); CHKERRQ(ierr);
//    
//    //build problem
//    ierr = mField.build_problems(); CHKERRQ(ierr);
//    
//    
//    
//    /****/
//    //mesh partitioning
//    
//    //partition
//    ierr = mField.partition_problem("ELASTIC_MECHANICS"); CHKERRQ(ierr);
//    ierr = mField.partition_finite_elements("ELASTIC_MECHANICS"); CHKERRQ(ierr);
//    //what are ghost nodes, see Petsc Manual
//    ierr = mField.partition_ghost_dofs("ELASTIC_MECHANICS"); CHKERRQ(ierr);
//    
//    //print bcs
//    ierr = mField.printCubitDisplacementSet(); CHKERRQ(ierr);
//    ierr = mField.printCubitForceSet(); CHKERRQ(ierr);
//    //print block sets with materials
//    ierr = mField.printCubitMaterials(); CHKERRQ(ierr);
//    
//    //    for(_IT_GET_DOFS_FIELD_BY_NAME_AND_TYPE_FOR_LOOP_(mField,"Lagrange_mul_disp",MBEDGE,dof)) {
//    //        cerr << *dof << endl;
//    //
//    //    }
//    
//    //create matrices (here F, D and Aij are matrices for the full problem)
//    Vec F,D;
//    ierr = mField.VecCreateGhost("ELASTIC_MECHANICS",Row,&F); CHKERRQ(ierr);
//    ierr = mField.VecCreateGhost("ELASTIC_MECHANICS",Col,&D); CHKERRQ(ierr);
//    
//    Mat Aij;
//    ierr = mField.MatCreateMPIAIJWithArrays("ELASTIC_MECHANICS",&Aij); CHKERRQ(ierr);
//    
//    struct MyElasticFEMethod: public ElasticFEMethod {
//        MyElasticFEMethod(FieldInterface& _mField,BaseDirihletBC *_dirihlet_ptr,
//                          Mat &_Aij,Vec &_D,Vec& _F,double _lambda,double _mu):
//        ElasticFEMethod(_mField,_dirihlet_ptr,_Aij,_D,_F,_lambda,_mu) {};
//        
//        PetscErrorCode Fint(Vec F_int) {
//            PetscFunctionBegin;
//            ierr = ElasticFEMethod::Fint(); CHKERRQ(ierr);
//            for(int rr = 0;rr<row_mat;rr++) {
//                if(RowGlob[rr].size()!=f_int[rr].size()) SETERRQ(PETSC_COMM_SELF,1,"data inconsistency");
//                if(RowGlob[rr].size()==0) continue;
//                f_int[rr] *= -1; //This is not SNES we solve K*D = -RES
//                ierr = VecSetValues(F_int,RowGlob[rr].size(),&(RowGlob[rr])[0],&(f_int[rr].data()[0]),ADD_VALUES); CHKERRQ(ierr);
//            }
//            PetscFunctionReturn(0);
//        }
//        
//    };
//    
//    
//    CubitDisplacementDirihletBC myDirihletBC(mField,"ELASTIC_MECHANICS","DISPLACEMENT");
//    ierr = myDirihletBC.Init(); CHKERRQ(ierr);
//    
//    //Assemble F and Aij
//    const double YoungModulus = 1;
//    const double PoissonRatio = 0.0;
//    MyElasticFEMethod MyFE(mField,&myDirihletBC,Aij,D,F,LAMBDA(YoungModulus,PoissonRatio),MU(YoungModulus,PoissonRatio));
//    
//    ElasticFE_RVELagrange MyFE_RVELagrange(mField,&myDirihletBC,Aij,D,F);
//    
//    ierr = VecZeroEntries(F); CHKERRQ(ierr);
//    ierr = VecGhostUpdateBegin(F,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
//    ierr = VecGhostUpdateEnd(F,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
//    ierr = MatZeroEntries(Aij); CHKERRQ(ierr);
//    
//    ierr = mField.loop_finite_elements("ELASTIC_MECHANICS","ELASTIC",MyFE);  CHKERRQ(ierr);
//    ierr = mField.loop_finite_elements("ELASTIC_MECHANICS","Lagrange_elem",MyFE_RVELagrange);  CHKERRQ(ierr);
//    
//    ierr = VecGhostUpdateBegin(F,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
//    ierr = VecGhostUpdateEnd(F,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
//    ierr = VecAssemblyBegin(F); CHKERRQ(ierr);
//    ierr = VecAssemblyEnd(F); CHKERRQ(ierr);
//    ierr = MatAssemblyBegin(Aij,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
//    ierr = MatAssemblyEnd(Aij,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
//    
//    PetscSynchronizedFlush(PETSC_COMM_WORLD);
//    //ierr = VecView(F,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);
//    //ierr = MatView(Aij,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);
//    
//    
//    ////Matrix View
//    //MatView(Aij,PETSC_VIEWER_DRAW_WORLD);//PETSC_VIEWER_STDOUT_WORLD);
//    //std::string wait;
//    //std::cin >> wait;
//    
//    //Solver
//    KSP solver;
//    ierr = KSPCreate(PETSC_COMM_WORLD,&solver); CHKERRQ(ierr);
//    ierr = KSPSetOperators(solver,Aij,Aij,SAME_NONZERO_PATTERN); CHKERRQ(ierr);
//    ierr = KSPSetFromOptions(solver); CHKERRQ(ierr);
//    ierr = KSPSetUp(solver); CHKERRQ(ierr);
//    
//    ierr = KSPSolve(solver,F,D); CHKERRQ(ierr);
//    ierr = VecGhostUpdateBegin(D,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
//    ierr = VecGhostUpdateEnd(D,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
//    
//    //  ierr = VecView(F,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);
//    //     ierr = VecView(D,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);
//    //
//    
//    //Save data on mesh
//    ierr = mField.set_global_VecCreateGhost("ELASTIC_MECHANICS",Row,D,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
//    //ierr = VecView(F,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);
//    
//    PostProcVertexMethod ent_method(moab);
//    ierr = mField.loop_dofs("ELASTIC_MECHANICS","DISPLACEMENT",Row,ent_method); CHKERRQ(ierr);
//    
//    if(pcomm->rank()==0) {
//        EntityHandle out_meshset;
//        rval = moab.create_meshset(MESHSET_SET,out_meshset); CHKERR_PETSC(rval);
//        ierr = mField.problem_get_FE("ELASTIC_MECHANICS","ELASTIC",out_meshset); CHKERRQ(ierr);
//        rval = moab.write_file("out.vtk","VTK","",&out_meshset,1); CHKERR_PETSC(rval);
//        rval = moab.delete_entities(&out_meshset,1); CHKERR_PETSC(rval);
//    }
//    
//    //PostProcDisplacemenysAndStarinOnRefMesh fe_post_proc_method(moab);
//    PostProcDisplacemenysAndStarinAndElasticLinearStressOnRefMesh fe_post_proc_method(moab,"DISPLACEMENT",LAMBDA(YoungModulus,PoissonRatio),MU(YoungModulus,PoissonRatio));
//    ierr = mField.loop_finite_elements("ELASTIC_MECHANICS","ELASTIC",fe_post_proc_method);  CHKERRQ(ierr);
//    
//    if(pcomm->rank()==0) {
//        rval = fe_post_proc_method.moab_post_proc.write_file("out_post_proc.vtk","VTK",""); CHKERR_PETSC(rval);
//    }
//    
//    //End Disp
//    /*Range ents;*/
//    //ierr = mField.get_Cubit_msId_entities_by_dimension(1,NodeSet,0,ents,true); CHKERRQ(ierr);
//    //for(_IT_GET_DOFS_FIELD_BY_NAME_FOR_LOOP_(mField,"DISPLACEMENT",dit)) {
//    //if(find(ents.begin(),ents.end(),dit->get_ent())!=ents.end()) {
//    //PetscSynchronizedPrintf(PETSC_COMM_WORLD, "val = %6.7e\n",dit->get_FieldData());
//    //}
//    /*}*/
//    
//    //Support stresses
//    //Range ents;
//    //ierr = mField.get_Cubit_msId_entities_by_dimension(4,NodeSet,0,ents,true); CHKERRQ(ierr);
//    
//    //Destroy matrices
//    ierr = VecDestroy(&F); CHKERRQ(ierr);
//    ierr = VecDestroy(&D); CHKERRQ(ierr);
//    ierr = MatDestroy(&Aij); CHKERRQ(ierr);
//    ierr = KSPDestroy(&solver); CHKERRQ(ierr);
    
    
    ierr = PetscTime(&v2);CHKERRQ(ierr);
    ierr = PetscGetCPUTime(&t2);CHKERRQ(ierr);
    
    PetscSynchronizedPrintf(PETSC_COMM_WORLD,"Total Rank %d Time = %f CPU Time = %f\n",pcomm->rank(),v2-v1,t2-t1);
    PetscSynchronizedFlush(PETSC_COMM_WORLD);
    
    PetscFinalize();
    
}

