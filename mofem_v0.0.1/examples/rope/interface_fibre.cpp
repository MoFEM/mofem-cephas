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

#include "ElasticFEMethodTransIso.hpp"
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
        order = 3;
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
    ierr = PetscGetTime(&v1); CHKERRQ(ierr);
    ierr = PetscGetCPUTime(&t1); CHKERRQ(ierr);
    
    //Create MoFEM (Joseph) database
    FieldCore core(moab);
    FieldInterface& mField = core;
    
    //ref meshset ref level 0
    ierr = mField.seed_ref_level_3D(0,0); CHKERRQ(ierr);
    
    EntityHandle meshset_BlockSet1; //Dirihlet BC is there
    ierr = mField.get_msId_meshset(1,BlockSet,meshset_BlockSet1); CHKERRQ(ierr);
    EntityHandle meshset_BlockSet2; //Dirihlet BC is there
    ierr = mField.get_msId_meshset(2,BlockSet,meshset_BlockSet2); CHKERRQ(ierr);
    
    //Interface meshset 4
    EntityHandle meshset_interface;
    ierr = mField.get_msId_meshset(4,SideSet,meshset_interface); CHKERRQ(ierr);
    ierr = mField.get_msId_3dENTS_sides(meshset_interface,true); CHKERRQ(ierr);
    
    // stl::bitset see for more details
    BitRefLevel bit_level_interface;
    bit_level_interface.set(0);
    ierr = mField.get_msId_3dENTS_split_sides(0,bit_level_interface,meshset_interface,true,true); CHKERRQ(ierr);
    EntityHandle meshset_level_interface;
    rval = moab.create_meshset(MESHSET_SET,meshset_level_interface); CHKERR_PETSC(rval);
    ierr = mField.refine_get_ents(bit_level_interface,BitRefLevel().set(),meshset_level_interface); CHKERRQ(ierr);
    
    ierr = mField.refine_get_childern(meshset_BlockSet1,bit_level_interface,meshset_BlockSet1,MBTET,true); CHKERRQ(ierr);
    ierr = mField.refine_get_childern(meshset_BlockSet2,bit_level_interface,meshset_BlockSet2,MBTET,true); CHKERRQ(ierr);
    
    for(_IT_CUBITMESHSETS_FOR_LOOP_(mField,cubit_it)) {
        EntityHandle cubit_meshset = cubit_it->meshset; 
        ierr = mField.refine_get_childern(cubit_meshset,meshset_level_interface,cubit_meshset,MBTRI,true); CHKERRQ(ierr);
    }    
    
    EntityHandle temp_meshset;
    rval = moab.create_meshset(MESHSET_SET,temp_meshset); CHKERR_PETSC(rval);
    ierr = mField.refine_get_ents(bit_level_interface,BitRefLevel().set(),temp_meshset); CHKERRQ(ierr);
    
    //Interface meshset 5
    EntityHandle meshset_interface1;
    ierr = mField.get_msId_meshset(5,SideSet,meshset_interface1); CHKERRQ(ierr);
    ierr = mField.get_msId_3dENTS_sides(meshset_interface1,true); CHKERRQ(ierr);
    
    BitRefLevel bit_level_interface1;
    bit_level_interface1.set(1);
    ierr = mField.get_msId_3dENTS_split_sides(0,bit_level_interface1,meshset_interface1,true,true); CHKERRQ(ierr);
    EntityHandle meshset_level_interface1;
    rval = moab.create_meshset(MESHSET_SET,meshset_level_interface1); CHKERR_PETSC(rval);
    ierr = mField.refine_get_ents(bit_level_interface1,BitRefLevel().set(),meshset_level_interface1); CHKERRQ(ierr);
    
    ierr = mField.refine_get_childern(meshset_BlockSet1,bit_level_interface1,meshset_BlockSet1,MBTET,true); CHKERRQ(ierr);
    ierr = mField.refine_get_childern(meshset_BlockSet2,bit_level_interface1,meshset_BlockSet2,MBTET,true); CHKERRQ(ierr);
    
    for(_IT_CUBITMESHSETS_FOR_LOOP_(mField,cubit_it)) {
        EntityHandle cubit_meshset = cubit_it->meshset; 
        ierr = mField.refine_get_childern(cubit_meshset,meshset_level_interface1,cubit_meshset,MBTRI,true); CHKERRQ(ierr);
    }
    
    //Interface meshset 6
    EntityHandle meshset_interface2;
    ierr = mField.get_msId_meshset(6,SideSet,meshset_interface2); CHKERRQ(ierr);
    ierr = mField.get_msId_3dENTS_sides(meshset_interface2,true); CHKERRQ(ierr);
    
    BitRefLevel bit_level_interface2;
    bit_level_interface2.set(2);
    ierr = mField.get_msId_3dENTS_split_sides(0,bit_level_interface2,meshset_interface2,true,true); CHKERRQ(ierr);
    EntityHandle meshset_level_interface2;
    rval = moab.create_meshset(MESHSET_SET,meshset_level_interface2); CHKERR_PETSC(rval);
    ierr = mField.refine_get_ents(bit_level_interface2,BitRefLevel().set(),meshset_level_interface2); CHKERRQ(ierr);
    
    ierr = mField.refine_get_childern(meshset_BlockSet1,bit_level_interface2,meshset_BlockSet1,MBTET,true); CHKERRQ(ierr);
    ierr = mField.refine_get_childern(meshset_BlockSet2,bit_level_interface2,meshset_BlockSet2,MBTET,true); CHKERRQ(ierr);
    
    for(_IT_CUBITMESHSETS_FOR_LOOP_(mField,cubit_it)) {
        EntityHandle cubit_meshset = cubit_it->meshset; 
        ierr = mField.refine_get_childern(cubit_meshset,meshset_level_interface2,cubit_meshset,MBTRI,true); CHKERRQ(ierr);
    }
    
    // stl::bitset see for more details
    BitRefLevel bit_level0;
    bit_level0.set(1);
    EntityHandle meshset_level0;
    rval = moab.create_meshset(MESHSET_SET,meshset_level0); CHKERR_PETSC(rval);
    ierr = mField.seed_ref_level_3D(meshset_level_interface2,bit_level0); CHKERRQ(ierr);
    ierr = mField.refine_get_ents(bit_level0,BitRefLevel().set(),meshset_level0); CHKERRQ(ierr);
    
    Range TETSFormLevel0;
    Range TETsFromBlockSet1onmeshset_level0;
    rval = moab.get_entities_by_type(meshset_level0,MBTET,TETSFormLevel0,true); CHKERR_PETSC(rval);
    rval = moab.get_entities_by_type(meshset_BlockSet1,MBTET,TETsFromBlockSet1onmeshset_level0,true); CHKERR_PETSC(rval);
    TETsFromBlockSet1onmeshset_level0 = intersect(TETsFromBlockSet1onmeshset_level0,TETSFormLevel0);
    EntityHandle meshset_BlockSet1OnLevel0;
    rval = moab.create_meshset(MESHSET_SET,meshset_BlockSet1OnLevel0); CHKERR_PETSC(rval);
    rval = moab.add_entities(meshset_BlockSet1OnLevel0,TETsFromBlockSet1onmeshset_level0); CHKERR_PETSC(rval);
    
    
    Range TETsFromBlockSet2onmeshset_level0;
    rval = moab.get_entities_by_type(meshset_BlockSet2,MBTET,TETsFromBlockSet2onmeshset_level0,true); CHKERR_PETSC(rval);
    EntityHandle meshset_BlockSet2OnLevel0;
    rval = moab.create_meshset(MESHSET_SET,meshset_BlockSet2OnLevel0); CHKERR_PETSC(rval);
    TETsFromBlockSet2onmeshset_level0 = intersect(TETsFromBlockSet2onmeshset_level0,TETSFormLevel0);
    rval = moab.add_entities(meshset_BlockSet2OnLevel0,TETsFromBlockSet2onmeshset_level0); CHKERR_PETSC(rval);
    
    BitRefLevel problem_bit_level = bit_level0;
        
    /***/
    //Define problem
    
    //Fields
    ierr = mField.add_field("DISPLACEMENT",H1,3); CHKERRQ(ierr);
    
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
    
    //add finite elements entities
    ierr = mField.add_ents_to_finite_element_by_TETs(meshset_BlockSet1OnLevel0,"ELASTIC",true); CHKERRQ(ierr);
    ierr = mField.add_ents_to_finite_element_by_TETs(meshset_BlockSet2OnLevel0,"TRAN_ISOTROPIC_ELASTIC",true); CHKERRQ(ierr);
    ierr = mField.add_ents_to_finite_element_EntType_by_bit_ref(problem_bit_level,"INTERFACE",MBPRISM); CHKERRQ(ierr);    
    //set app. order
    //see Hierarchic Finite Element Bases on Unstructured Tetrahedral Meshes (Mark Ainsworth & Joe Coyle)
    ierr = mField.set_field_order(0,MBTET,"DISPLACEMENT",order); CHKERRQ(ierr);
    ierr = mField.set_field_order(0,MBTRI,"DISPLACEMENT",order); CHKERRQ(ierr);
    ierr = mField.set_field_order(0,MBEDGE,"DISPLACEMENT",order); CHKERRQ(ierr);
    ierr = mField.set_field_order(0,MBVERTEX,"DISPLACEMENT",1); CHKERRQ(ierr);
    
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
    Vec F,D;
    ierr = mField.VecCreateGhost("ELASTIC_MECHANICS",Row,&F); CHKERRQ(ierr);
    ierr = mField.VecCreateGhost("ELASTIC_MECHANICS",Col,&D); CHKERRQ(ierr);

    Mat Aij;
    ierr = mField.MatCreateMPIAIJWithArrays("ELASTIC_MECHANICS",&Aij); CHKERRQ(ierr);

    struct MyElasticFEMethod: public ElasticFEMethod {
        MyElasticFEMethod(FieldInterface& _mField,BaseDirihletBC *_dirihlet_ptr,
                          Mat &_Aij,Vec &_D,Vec& _F,double _lambda,double _mu): 
        ElasticFEMethod(_mField,_dirihlet_ptr,_Aij,_D,_F,_lambda,_mu) {};
        
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
    
    CubitDisplacementDirihletBC myDirihletBC(mField,"ELASTIC_MECHANICS","DISPLACEMENT");
    ierr = myDirihletBC.Init(); CHKERRQ(ierr);
    
    Tag th_phi;
    double def_val = 0;
    rval = moab.tag_get_handle("PHI",1,MB_TYPE_DOUBLE,th_phi,MB_TAG_CREAT|MB_TAG_SPARSE,&def_val); CHKERR_PETSC(rval);
    
    for(_IT_GET_DOFS_FIELD_BY_NAME_FOR_LOOP_(mField,"POTENTIAL_FIELD",dof_ptr)) {
        if(dof_ptr->get_ent_type()!=MBVERTEX) continue;
        EntityHandle ent = dof_ptr->get_ent();
        double &fval = dof_ptr->get_FieldData();
        double phi;
        rval = moab.tag_get_data(th_phi,&ent,1,&phi); CHKERR_PETSC(rval);  
        fval = phi;
    }
    
    //Assemble F and Aij
    //    vector<double> attributes;
    double YoungModulusP;
    double PoissonRatioP;
    double YoungModulusZ;
    double PoissonRatioPZ;
    double ShearModulusZP;
    double YoungModulus;
    double PoissonRatio;
    const double alpha = 0.05;
    
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
    }
        
    InterfaceFEMethod IntMyFE(mField,&myDirihletBC,Aij,D,F,YoungModulus*alpha);
    MyElasticFEMethod MyFE(mField,&myDirihletBC,Aij,D,F,LAMBDA(YoungModulus,PoissonRatio),MU(YoungModulus,PoissonRatio));
    TranIsotropicFibreDirRotElasticFEMethod MyTIsotFE(mField,&myDirihletBC,Aij,D,F,LAMBDA(YoungModulusP,PoissonRatioP),MU(YoungModulusP,PoissonRatioP),YoungModulusP,YoungModulusZ,PoissonRatioP,PoissonRatioPZ,ShearModulusZP);
    
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
    
    TranIso_PostProc_FibreDirRot_OnRefMesh fe_post_proc_method( moab, LAMBDA(YoungModulusP,PoissonRatioP),MU(YoungModulusP,PoissonRatioP), YoungModulusP,YoungModulusZ,PoissonRatioP,PoissonRatioPZ,ShearModulusZP);
    
//    ierr = mField.loop_finite_elements("ELASTIC_MECHANICS","ELASTIC",fe_post_proc_method);  CHKERRQ(ierr);
    ierr = mField.loop_finite_elements("ELASTIC_MECHANICS","TRAN_ISOTROPIC_ELASTIC",fe_post_proc_method);  CHKERRQ(ierr);
    
    PetscSynchronizedFlush(PETSC_COMM_WORLD);
    if(pcomm->rank()==0) {
        rval = fe_post_proc_method.moab_post_proc.write_file(outName2,"VTK",""); CHKERR_PETSC(rval);
    }
    
    PostProcCohesiveForces fe_post_proc_prisms(mField,YoungModulus*alpha);
    ierr = mField.loop_finite_elements("ELASTIC_MECHANICS","INTERFACE",fe_post_proc_prisms);  CHKERRQ(ierr);
    PetscSynchronizedFlush(PETSC_COMM_WORLD);
    if(pcomm->rank()==0) {
        rval = fe_post_proc_prisms.moab_post_proc.write_file("out_post_proc_prisms.vtk","VTK",""); CHKERR_PETSC(rval);
    }
    
    
    //detroy matrices
    ierr = VecDestroy(&F); CHKERRQ(ierr);
    ierr = VecDestroy(&D); CHKERRQ(ierr);
    ierr = MatDestroy(&Aij); CHKERRQ(ierr);
    ierr = KSPDestroy(&solver); CHKERRQ(ierr);
    
    
    ierr = PetscGetTime(&v2);CHKERRQ(ierr);
    ierr = PetscGetCPUTime(&t2);CHKERRQ(ierr);
    
    PetscSynchronizedPrintf(PETSC_COMM_WORLD,"Total Rank %d Time = %f CPU Time = %f\n",pcomm->rank(),v2-v1,t2-t1);
    PetscSynchronizedFlush(PETSC_COMM_WORLD);
    
    PetscFinalize();
    
}
