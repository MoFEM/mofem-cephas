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

#include "ElasticFEMethod_RVE_HomoStress.hpp"
#include "ElasticFEMethodTransIsoRVE_HomoStress.hpp"

#include "ElasticFEMethod_RVE_Try.hpp"
#include "ElasticFEMethodTransIsoRVE.hpp"
#include "PostProcVertexMethod.hpp"
#include "PostProcDisplacementAndStrainOnRefindedMesh.hpp"
#include "ElasticFEMethodForInterface.hpp"

using namespace MoFEM;

ErrorCode rval;
PetscErrorCode ierr;

static char help[] = "...\n\n";

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

    vector<BitRefLevel> bit_levels;
    bit_levels.push_back(BitRefLevel().set(4));

//    const clock_t begin_time = clock();
    ierr = mField.build_fields(); CHKERRQ(ierr);
    ierr = mField.build_finite_elements(); CHKERRQ(ierr);
    ierr = mField.build_adjacencies(bit_levels.back()); CHKERRQ(ierr);
    ierr = mField.build_problems(); CHKERRQ(ierr);

    Range RangeFibre1, RangeFibre2, RangeFibre3, RangeFibre4;
    
	for(_IT_GET_FES_BY_NAME_FOR_LOOP_(mField,"POTENTIAL_ELEM1",it)){
		RangeFibre1.insert(it->get_ent());
	}
	for(_IT_GET_FES_BY_NAME_FOR_LOOP_(mField,"POTENTIAL_ELEM2",it)){
		RangeFibre2.insert(it->get_ent());
	}
	for(_IT_GET_FES_BY_NAME_FOR_LOOP_(mField,"POTENTIAL_ELEM3",it)){
		RangeFibre3.insert(it->get_ent());
	}
	for(_IT_GET_FES_BY_NAME_FOR_LOOP_(mField,"POTENTIAL_ELEM4",it)){
		RangeFibre4.insert(it->get_ent());
	}
	
    EntityHandle fibre1_meshset, fibre2_meshset, fibre3_meshset, fibre4_meshset;
	
    rval = moab.create_meshset(MESHSET_SET,fibre1_meshset); CHKERR_PETSC(rval);
	rval = moab.create_meshset(MESHSET_SET,fibre2_meshset); CHKERR_PETSC(rval);
	rval = moab.create_meshset(MESHSET_SET,fibre3_meshset); CHKERR_PETSC(rval);
	rval = moab.create_meshset(MESHSET_SET,fibre4_meshset); CHKERR_PETSC(rval);

	rval = moab.add_entities(fibre1_meshset,RangeFibre1); CHKERR_PETSC(rval);
	rval = moab.add_entities(fibre2_meshset,RangeFibre2); CHKERR_PETSC(rval);
	rval = moab.add_entities(fibre3_meshset,RangeFibre3); CHKERR_PETSC(rval);
	rval = moab.add_entities(fibre4_meshset,RangeFibre4); CHKERR_PETSC(rval);
  
	rval = moab.write_file("fibre1.vtk","VTK","",&fibre1_meshset,1); CHKERR_PETSC(rval);
	rval = moab.write_file("fibre2.vtk","VTK","",&fibre2_meshset,1); CHKERR_PETSC(rval);
	rval = moab.write_file("fibre3.vtk","VTK","",&fibre3_meshset,1); CHKERR_PETSC(rval);
	rval = moab.write_file("fibre4.vtk","VTK","",&fibre4_meshset,1); CHKERR_PETSC(rval);

//    std::cout << float( clock () - begin_time ) /  CLOCKS_PER_SEC<<endl;
    
    EntityHandle out_meshset;
    rval = moab.create_meshset(MESHSET_SET,out_meshset); CHKERR_PETSC(rval);
//    ierr = mField.problem_get_FE("POTENTIAL_PROBLEM","POTENTIAL_ELEM",out_meshset); CHKERRQ(ierr);
    ierr = mField.refine_get_ents(bit_levels.back(),BitRefLevel().set(),MBTET,out_meshset); CHKERRQ(ierr);
    Range LatestRefinedTets;
    rval = moab.get_entities_by_type(out_meshset, MBTET,LatestRefinedTets,true); CHKERR_PETSC(rval);

    BitRefLevel problem_bit_level = bit_levels.back();
    
    EntityHandle meshset_Elastic, meshset_Trans_ISO;
    rval = moab.create_meshset(MESHSET_SET,meshset_Elastic); CHKERR_PETSC(rval);
    rval = moab.create_meshset(MESHSET_SET,meshset_Trans_ISO); CHKERR_PETSC(rval);

	rval = moab.unite_meshset(meshset_Trans_ISO,fibre1_meshset); CHKERR_PETSC(rval);
	rval = moab.unite_meshset(meshset_Trans_ISO,fibre2_meshset); CHKERR_PETSC(rval);
	rval = moab.unite_meshset(meshset_Trans_ISO,fibre3_meshset); CHKERR_PETSC(rval);
	rval = moab.unite_meshset(meshset_Trans_ISO,fibre4_meshset); CHKERR_PETSC(rval);
	
    rval = moab.write_file("meshset_Trans_ISO.vtk","VTK","",&meshset_Trans_ISO,1); CHKERR_PETSC(rval);

	for(_IT_CUBITMESHSETS_BY_SET_TYPE_FOR_LOOP_(mField,BlockSet,it)){
        
		if(it->get_Cubit_name() == "MAT_ELASTIC_1") {
			Range TetsInBlock;
			rval = moab.get_entities_by_type(it->meshset, MBTET,TetsInBlock,true); CHKERR_PETSC(rval);
			Range block_rope_bit_level = intersect(LatestRefinedTets,TetsInBlock);   // To find new block (including new nodes after inserting interfaces), here LatestRefinedTets is in the final bit_level, while TetsInBlock are in the start bit levle.
			rval = moab.add_entities(meshset_Elastic,block_rope_bit_level);CHKERR_PETSC(rval);
		}
	}
    
    rval = moab.write_file("meshset_Elastic_Zahur.vtk","VTK","",&meshset_Elastic,1); CHKERR_PETSC(rval);
  
    /***/
    //Define problem
    
    //Fields
    ierr = mField.add_field("DISPLACEMENT",H1,3); CHKERRQ(ierr);
    ierr = mField.add_field("DISPLACEMENT_FROM_APP_STRAIN",H1,3); CHKERRQ(ierr);
    
    //FE
    ierr = mField.add_finite_element("ELASTIC"); CHKERRQ(ierr);
    ierr = mField.add_finite_element("TRAN_ISOTROPIC_ELASTIC"); CHKERRQ(ierr);
    ierr = mField.add_finite_element("INTERFACE"); CHKERRQ(ierr);
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
    
    
    ierr = mField.modify_finite_element_add_field_data("TRAN_ISOTROPIC_ELASTIC","DISPLACEMENT_FROM_APP_STRAIN"); CHKERRQ(ierr);
    ierr = mField.modify_finite_element_add_field_data("ELASTIC","DISPLACEMENT_FROM_APP_STRAIN"); CHKERRQ(ierr);
    ierr = mField.modify_finite_element_add_field_data("INTERFACE","DISPLACEMENT_FROM_APP_STRAIN"); CHKERRQ(ierr);

    
    //define problems
    ierr = mField.add_problem("ELASTIC_MECHANICS"); CHKERRQ(ierr);
    
    //set finite elements for problem
    ierr = mField.modify_problem_add_finite_element("ELASTIC_MECHANICS","ELASTIC"); CHKERRQ(ierr);
    ierr = mField.modify_problem_add_finite_element("ELASTIC_MECHANICS","INTERFACE"); CHKERRQ(ierr);
    ierr = mField.modify_problem_add_finite_element("ELASTIC_MECHANICS","TRAN_ISOTROPIC_ELASTIC"); CHKERRQ(ierr);
    
    //set refinment level for problem
    ierr = mField.modify_problem_ref_level_add_bit("ELASTIC_MECHANICS",problem_bit_level); CHKERRQ(ierr);

    /***/
    //Declare problem
    
    //add entitities (by tets) to the field
    ierr = mField.add_ents_to_field_by_TETs(0,"DISPLACEMENT"); CHKERRQ(ierr);
    ierr = mField.add_ents_to_field_by_TETs(0,"DISPLACEMENT_FROM_APP_STRAIN"); CHKERRQ(ierr);

    //add finite elements entities
    ierr = mField.add_ents_to_finite_element_by_TETs(meshset_Elastic,"ELASTIC",true); CHKERRQ(ierr);
    ierr = mField.add_ents_to_finite_element_by_TETs(meshset_Trans_ISO,"TRAN_ISOTROPIC_ELASTIC",true); CHKERRQ(ierr);
    ierr = mField.add_ents_to_finite_element_EntType_by_bit_ref(problem_bit_level,"INTERFACE",MBPRISM); CHKERRQ(ierr);    
    //set app. order
    //see Hierarchic Finite Element Bases on Unstructured Tetrahedral Meshes (Mark Ainsworth & Joe Coyle)
    ierr = mField.set_field_order(0,MBTET,"DISPLACEMENT",order); CHKERRQ(ierr);
    ierr = mField.set_field_order(0,MBTRI,"DISPLACEMENT",order); CHKERRQ(ierr);
    ierr = mField.set_field_order(0,MBEDGE,"DISPLACEMENT",order); CHKERRQ(ierr);
    ierr = mField.set_field_order(0,MBVERTEX,"DISPLACEMENT",1); CHKERRQ(ierr);
    
    //see Hierarchic Finite Element Bases on Unstructured Tetrahedral Meshes (Mark Ainsworth & Joe Coyle)
    ierr = mField.set_field_order(0,MBTET,"DISPLACEMENT_FROM_APP_STRAIN",order); CHKERRQ(ierr);
    ierr = mField.set_field_order(0,MBTRI,"DISPLACEMENT_FROM_APP_STRAIN",order); CHKERRQ(ierr);
    ierr = mField.set_field_order(0,MBEDGE,"DISPLACEMENT_FROM_APP_STRAIN",order); CHKERRQ(ierr);
    ierr = mField.set_field_order(0,MBVERTEX,"DISPLACEMENT_FROM_APP_STRAIN",1); CHKERRQ(ierr);
    
    
    /****/
    //build database
    
    //build field
    ierr = mField.build_fields(); CHKERRQ(ierr);
    
    //build finite elemnts
    ierr = mField.build_finite_elements(); CHKERRQ(ierr);
    
    //build adjacencies
    ierr = mField.build_adjacencies(problem_bit_level); CHKERRQ(ierr);
    
    //build problem
    ierr = mField.build_problems(); CHKERRQ(ierr);
    
    /****/
    //mesh partitioning 
    
    //partition
    ierr = mField.partition_problem("ELASTIC_MECHANICS"); CHKERRQ(ierr);
    ierr = mField.partition_finite_elements("ELASTIC_MECHANICS"); CHKERRQ(ierr);
    //what are ghost nodes, see Petsc Manual
    ierr = mField.partition_ghost_dofs("ELASTIC_MECHANICS"); CHKERRQ(ierr);
    
    //create matrices
    Vec F,D,D_star,F_stress,coord_stress;
    ierr = mField.VecCreateGhost("ELASTIC_MECHANICS",Row,&F); CHKERRQ(ierr);
    ierr = mField.VecCreateGhost("ELASTIC_MECHANICS",Col,&D); CHKERRQ(ierr);
    ierr = mField.VecCreateGhost("ELASTIC_MECHANICS",Col,&D_star); CHKERRQ(ierr);
    ierr = mField.VecCreateGhost("ELASTIC_MECHANICS",Row,&F_stress); CHKERRQ(ierr);
    ierr = mField.VecCreateGhost("ELASTIC_MECHANICS",Row,&coord_stress); CHKERRQ(ierr);
    
    Mat Aij;
    ierr = mField.MatCreateMPIAIJWithArrays("ELASTIC_MECHANICS",&Aij); CHKERRQ(ierr);

    
       //Applied strain (specified by the user)
    ublas::matrix<double> strain_app;
    strain_app.resize(3,3);
    strain_app(0,0) = 0.01; strain_app(0,1)=0.0; strain_app(0,2)=0.0;
    strain_app(1,0) = 0.0; strain_app(1,1)=0.0; strain_app(1,2)=0.0;
    strain_app(2,0) = 0.0; strain_app(2,1)=0.0; strain_app(2,2)=0.0;
    
    //Apply the linear displacement to all nodes in the mesh
    Tag th_disp_1;
    double defaultval[3]={1,0,0};
    rval = moab.tag_get_handle("DISPLACEMENT_FROM_APP_STRAIN",3,MB_TYPE_DOUBLE,th_disp_1,MB_TAG_CREAT|MB_TAG_SPARSE,&defaultval); CHKERR_THROW(rval);
    
    EntityHandle node = 0;
    double coords[3], disp_applied[3];
    int count=0;
    for(_IT_GET_DOFS_FIELD_BY_NAME_FOR_LOOP_(mField,"DISPLACEMENT_FROM_APP_STRAIN",dof_ptr)) {
        if(dof_ptr->get_ent_type()!=MBVERTEX) continue;
        EntityHandle ent = dof_ptr->get_ent();
        int dof_rank = dof_ptr->get_dof_rank();
        double &fval = dof_ptr->get_FieldData();
        if(node!=ent) {    // to get coordinates only for 1 out of 3 ranks for disp field
            rval = moab.get_coords(&ent,1,coords); CHKERR_PETSC(rval);
            //cout<<"\n\n coord = " << coords[0]<<" "<< coords[1]<<" " << coords[2]<<" rank "<< dof_rank << "\n\n";
            for(int ii=0; ii<3; ii++) disp_applied[ii]=strain_app(ii,0)*coords[0] + strain_app(ii,1)*coords[1] + strain_app(ii,2)*coords[2];
            rval=moab.tag_set_data(th_disp_1,&ent,1,&disp_applied); CHKERR_PETSC(rval);
            node = ent;
        }
        fval = disp_applied[dof_rank];
    }
    
    //save the mesh to see in the paraview
    if(pcomm->rank()==0) {
        EntityHandle out_meshset;
        rval = moab.create_meshset(MESHSET_SET,out_meshset); CHKERR_PETSC(rval);
        //ierr = mField.problem_get_FE("ELASTIC_MECHANICS","ELASTIC",out_meshset); CHKERRQ(ierr);
        ierr = mField.refine_get_ents(problem_bit_level,BitRefLevel().set(),MBTET,out_meshset); CHKERRQ(ierr);
        rval = moab.write_file("Zahur_out_disp_from_app_strain.vtk","VTK","",&out_meshset,1); CHKERR_PETSC(rval);
        rval = moab.delete_entities(&out_meshset,1); CHKERR_PETSC(rval);
    }

    
    CubitDisplacementDirihletBC myDirihletBC(mField,"ELASTIC_MECHANICS","DISPLACEMENT");
    ierr = myDirihletBC.Init(); CHKERRQ(ierr);
    
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
    
    for(_IT_CUBITMESHSETS_BY_SET_TYPE_FOR_LOOP_(mField,BlockSet,it))
    {
        cout << endl << *it << endl;
        
        //Get block name
        string name = it->get_Cubit_name();
        
        if (name.compare(0,12,"MAT_TRANSISO") == 0)
        {
            Mat_TransIso mydata;
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
        
    InterfaceFEMethod IntMyFE(mField,&myDirihletBC,Aij,D,F,YoungModulus*alpha);
    ElasticFEMethod_RVE_Try MyFE(mField,&myDirihletBC,Aij,D,F,LAMBDA(YoungModulus,PoissonRatio),MU(YoungModulus,PoissonRatio),D_star,strain_app);
    
    ElasticFEMethodTransIsoRVE MyTIsotFE(mField,&myDirihletBC,Aij,D,F,LAMBDA(YoungModulusP,PoissonRatioP),MU(YoungModulusP,PoissonRatioP),YoungModulusP,YoungModulusZ,PoissonRatioP,PoissonRatioPZ,ShearModulusZP);

    ierr = VecZeroEntries(F); CHKERRQ(ierr);
    ierr = VecGhostUpdateBegin(F,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
    ierr = VecGhostUpdateEnd(F,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
    ierr = MatZeroEntries(Aij); CHKERRQ(ierr);
    
    ierr = mField.loop_finite_elements("ELASTIC_MECHANICS","INTERFACE",IntMyFE);  CHKERRQ(ierr);
    PetscSynchronizedFlush(PETSC_COMM_WORLD);
    ierr = mField.loop_finite_elements("ELASTIC_MECHANICS","ELASTIC",MyFE);  CHKERRQ(ierr);
    PetscSynchronizedFlush(PETSC_COMM_WORLD);
    ierr = mField.loop_finite_elements("ELASTIC_MECHANICS","TRAN_ISOTROPIC_ELASTIC",MyTIsotFE);  CHKERRQ(ierr);
    PetscSynchronizedFlush(PETSC_COMM_WORLD);
    
    ierr = VecGhostUpdateBegin(F,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
    ierr = VecGhostUpdateEnd(F,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
    ierr = VecAssemblyBegin(F); CHKERRQ(ierr);
    ierr = VecAssemblyEnd(F); CHKERRQ(ierr);
    ierr = MatAssemblyBegin(Aij,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
    ierr = MatAssemblyEnd(Aij,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
  

    //Matrix View
    //MatView(Aij,PETSC_VIEWER_DRAW_WORLD);//PETSC_VIEWER_STDOUT_WORLD);
    //std::string wait;
    //std::cin >> wait;
    
    //Solver
    KSP solver;
    ierr = KSPCreate(PETSC_COMM_WORLD,&solver); CHKERRQ(ierr);
    ierr = KSPSetOperators(solver,Aij,Aij,SAME_NONZERO_PATTERN); CHKERRQ(ierr);
    ierr = KSPSetFromOptions(solver); CHKERRQ(ierr);
    ierr = KSPSetUp(solver); CHKERRQ(ierr);
    
    ierr = KSPSolve(solver,F,D); CHKERRQ(ierr);
    ierr = VecGhostUpdateBegin(D,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
    ierr = VecGhostUpdateEnd(D,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
    
    ierr = VecAXPY(D,1,D_star); CHKERRQ(ierr);   // D=D+1*D_star (total displacement U_total=U_start+U_fluctuation)
    //Save data on mesh
    ierr = mField.set_global_VecCreateGhost("ELASTIC_MECHANICS",Row,D,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
    //ierr = VecView(F,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);

    PostProcVertexMethod ent_method(moab);
    ierr = mField.loop_dofs("ELASTIC_MECHANICS","DISPLACEMENT",Row,ent_method); CHKERRQ(ierr);
    
    if(pcomm->rank()==0) {
        EntityHandle out_meshset;
        rval = moab.create_meshset(MESHSET_SET,out_meshset); CHKERR_PETSC(rval);
        ierr = mField.problem_get_FE("ELASTIC_MECHANICS","ELASTIC",out_meshset); CHKERRQ(ierr);
        ierr = mField.problem_get_FE("ELASTIC_MECHANICS","TRAN_ISOTROPIC_ELASTIC",out_meshset); CHKERRQ(ierr);
        ierr = mField.problem_get_FE("ELASTIC_MECHANICS","INTERFACE",out_meshset); CHKERRQ(ierr);
        rval = moab.write_file(outName,"VTK","",&out_meshset,1); CHKERR_PETSC(rval);
        rval = moab.delete_entities(&out_meshset,1); CHKERR_PETSC(rval);
    }
    
    
    ierr = VecZeroEntries(F_stress); CHKERRQ(ierr);
    ierr = VecGhostUpdateBegin(F_stress,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
    ierr = VecGhostUpdateEnd(F_stress,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);

    double RVE_volume_Elastic, RVE_volume_TransIso, RVE_volume_total; RVE_volume_Elastic=0.0; RVE_volume_TransIso=0.0;
    ublas::matrix<double> Sigma_homo;
    
    ElasticFEMethod_RVE_HomoStress MyRVEFEStress(mField,&myDirihletBC,Aij,D,F,LAMBDA(YoungModulus,PoissonRatio),MU(YoungModulus,PoissonRatio), moab, F_stress, coord_stress, &RVE_volume_Elastic);
    
    ElasticFEMethodTransIsoRVE_HomoStress MyRVETransIsoStress(mField,&myDirihletBC,Aij,D,F,LAMBDA(YoungModulusP,PoissonRatioP),MU(YoungModulusP,PoissonRatioP),YoungModulusP,YoungModulusZ,PoissonRatioP,PoissonRatioPZ,ShearModulusZP, moab, F_stress, coord_stress,&RVE_volume_TransIso);
    
//    cout<<"\n RVE_volume_Elastic Before "<<RVE_volume_Elastic<<endl;
//    cout<<"\n RVE_volume_TransIso Before "<<RVE_volume_TransIso<<endl;
  
    ierr = mField.loop_finite_elements("ELASTIC_MECHANICS","ELASTIC",MyRVEFEStress);  CHKERRQ(ierr);
    ierr = mField.loop_finite_elements("ELASTIC_MECHANICS","TRAN_ISOTROPIC_ELASTIC",MyRVETransIsoStress);  CHKERRQ(ierr);

    ierr = VecGhostUpdateBegin(F_stress,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
    ierr = VecGhostUpdateEnd(F_stress,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
    ierr = VecAssemblyBegin(F_stress); CHKERRQ(ierr);
    ierr = VecAssemblyEnd(F_stress); CHKERRQ(ierr);

//    cout<<"\n RVE_volume_Elastic After "<<RVE_volume_Elastic<<endl;
//    cout<<"\n RVE_volume_TransIso After "<<RVE_volume_TransIso<<endl;
    RVE_volume_total=RVE_volume_Elastic+RVE_volume_TransIso;
    
    
//    cout<<"\n\n\n\n\n\n\n F_stress \n\n\n\n\n\n\\n\n\n\n\n";
//    ierr = VecView(F_stress,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);
//    cout<<"\n\n\n\n\n\n\n coord_stress \n\n\n\n\n\n\\n\n\n\n\n";
//    ierr = VecView(coord_stress,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);


    ublas::matrix< double > F_stress_mat;
    ublas::matrix< double > coord_stress_mat;
    Range nodes;
    rval = moab.get_entities_by_type(0,MBVERTEX,nodes,true); CHKERR_PETSC(rval);
    F_stress_mat.resize(3,nodes.size()); coord_stress_mat.resize(nodes.size(),3); Sigma_homo.resize(3,3);
    F_stress_mat.clear();  coord_stress_mat.clear(); Sigma_homo.clear();
    
    int array_indices[3];
    double array_output[3];
    
    for(int ii=0; ii<nodes.size(); ii++){
        array_indices[0]=3*ii;  array_indices[1]=3*ii+1;  array_indices[2]=3*ii+2;
        ierr = VecGetValues(F_stress,3,array_indices,array_output); CHKERRQ(ierr);
        F_stress_mat(0,ii)=array_output[0]; F_stress_mat(1,ii)=array_output[1],  F_stress_mat(2,ii)=array_output[2];
        
        ierr = VecGetValues(coord_stress,3,array_indices,array_output); CHKERRQ(ierr);
        coord_stress_mat(ii,0)=array_output[0]; coord_stress_mat(ii,1)=array_output[1], coord_stress_mat(ii,2)=array_output[2];
    }
    
    
    // Homogenized stress
    cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,
                Sigma_homo.size1(),Sigma_homo.size2(),F_stress_mat.size2(),
                1.,&*F_stress_mat.data().begin(),F_stress_mat.size2(),
                &*coord_stress_mat.data().begin(),coord_stress_mat.size2(),
                0.,&*Sigma_homo.data().begin(),Sigma_homo.size2());
    
    Sigma_homo=(1.0/RVE_volume_total)*(Sigma_homo+trans(Sigma_homo));
    
    //if(pcomm->rank()==1){
    //            cout<<"coord_stress_mat "<<coord_stress_mat<<endl;
    //            cout<<"F_stress_mat "<<F_stress_mat<<endl;  (1.0/RVE_volume)*
    cout<<"\n\n\nSigma_bar "<<Sigma_homo<<endl<<endl<<endl;

    
    
    
    
    
 //    TranIso_PostProc_FibreDirRot_OnRefMesh fe_post_proc_method( moab, LAMBDA(YoungModulusP,PoissonRatioP),MU(YoungModulusP,PoissonRatioP), YoungModulusP,YoungModulusZ,PoissonRatioP,PoissonRatioPZ,ShearModulusZP);
//    
////    ierr = mField.loop_finite_elements("ELASTIC_MECHANICS","ELASTIC",fe_post_proc_method);  CHKERRQ(ierr);
//    ierr = mField.loop_finite_elements("ELASTIC_MECHANICS","TRAN_ISOTROPIC_ELASTIC",fe_post_proc_method);  CHKERRQ(ierr);
//    
//    PetscSynchronizedFlush(PETSC_COMM_WORLD);
//    if(pcomm->rank()==0) {
//        rval = fe_post_proc_method.moab_post_proc.write_file(outName2,"VTK",""); CHKERR_PETSC(rval);
//    }
//    
//    PostProcCohesiveForces fe_post_proc_prisms(mField,YoungModulus*alpha);
//    ierr = mField.loop_finite_elements("ELASTIC_MECHANICS","INTERFACE",fe_post_proc_prisms);  CHKERRQ(ierr);
//    PetscSynchronizedFlush(PETSC_COMM_WORLD);
//    if(pcomm->rank()==0) {
//        rval = fe_post_proc_prisms.moab_post_proc.write_file("out_post_proc_prisms.vtk","VTK",""); CHKERR_PETSC(rval);
//    }
//    
    
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
