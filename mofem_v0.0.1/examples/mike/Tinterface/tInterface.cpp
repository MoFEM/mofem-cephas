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

#include "moabField.hpp"
#include "moabField_Core.hpp"
#include "moabFEMethod_UpLevelStudent.hpp"
#include "cholesky.hpp"
#include <petscksp.h>

#include "ElasticFEMethod.hpp"
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
    moabField_Core core(moab);
    moabField& mField = core;
    
    //ref meshset ref level 0
    ierr = mField.seed_ref_level_3D(0,0); CHKERRQ(ierr);
        
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
    ierr = mField.refine_get_ents(bit_level_interface,meshset_level_interface); CHKERRQ(ierr);
    
    //update BC for refined (with interface) mesh
    EntityHandle meshset_SideSet1; //Dirihlet BC is there
    ierr = mField.get_msId_meshset(1,SideSet,meshset_SideSet1); CHKERRQ(ierr);
    ierr = mField.refine_get_childern(meshset_SideSet1,bit_level_interface,meshset_SideSet1,MBTRI,true); CHKERRQ(ierr);
    EntityHandle meshset_SideSet2; //Dirihlet BC is there
    ierr = mField.get_msId_meshset(2,SideSet,meshset_SideSet2); CHKERRQ(ierr);
    ierr = mField.refine_get_childern(meshset_SideSet2,bit_level_interface,meshset_SideSet2,MBTRI,true); CHKERRQ(ierr);
//    EntityHandle meshset_BlockSet1; //Dirihlet BC is there
//    ierr = mField.get_msId_meshset(1,BlockSet,meshset_BlockSet1); CHKERRQ(ierr);
//    ierr = mField.refine_get_childern(meshset_BlockSet1,bit_level_interface,meshset_BlockSet1,MBTET,true); CHKERRQ(ierr);
//    EntityHandle meshset_BlockSet2; //Dirihlet BC is there
//    ierr = mField.get_msId_meshset(2,BlockSet,meshset_BlockSet2); CHKERRQ(ierr);
//    ierr = mField.refine_get_childern(meshset_BlockSet2,bit_level_interface,meshset_BlockSet2,MBTET,true); CHKERRQ(ierr);
//    EntityHandle meshset_BlockSet3; //Dirihlet BC is there
//    ierr = mField.get_msId_meshset(3,BlockSet,meshset_BlockSet3); CHKERRQ(ierr);
//    ierr = mField.refine_get_childern(meshset_BlockSet3,bit_level_interface,meshset_BlockSet3,MBTET,true); CHKERRQ(ierr);
    

    //Interface meshset 5
    EntityHandle meshset_interface1;
    ierr = mField.get_msId_meshset(5,SideSet,meshset_interface1); CHKERRQ(ierr);
    ierr = mField.get_msId_3dENTS_sides(meshset_interface1,true); CHKERRQ(ierr);

    BitRefLevel bit_level_interface1;
    bit_level_interface1.set(1);
    ierr = mField.get_msId_3dENTS_split_sides(0,bit_level_interface1,meshset_interface1,true,true); CHKERRQ(ierr);
    EntityHandle meshset_level_interface1;
    rval = moab.create_meshset(MESHSET_SET,meshset_level_interface1); CHKERR_PETSC(rval);
    ierr = mField.refine_get_ents(bit_level_interface1,meshset_level_interface1); CHKERRQ(ierr);
    
    ierr = mField.refine_get_childern(meshset_SideSet1, bit_level_interface1,meshset_SideSet1,MBTRI,true); CHKERRQ(ierr);
    ierr = mField.refine_get_childern(meshset_SideSet2, bit_level_interface1,meshset_SideSet2,MBTRI,true); CHKERRQ(ierr);
//    ierr = mField.refine_get_childern(meshset_BlockSet1,bit_level_interface1,meshset_BlockSet1,MBTET,true); CHKERRQ(ierr);
//    ierr = mField.refine_get_childern(meshset_BlockSet2,bit_level_interface1,meshset_BlockSet2,MBTET,true); CHKERRQ(ierr);
//    ierr = mField.refine_get_childern(meshset_BlockSet3,bit_level_interface1,meshset_BlockSet3,MBTET,true); CHKERRQ(ierr);

    //Interface meshset 6
    EntityHandle meshset_interface2;
    ierr = mField.get_msId_meshset(6,SideSet,meshset_interface2); CHKERRQ(ierr);
    ierr = mField.get_msId_3dENTS_sides(meshset_interface2,true); CHKERRQ(ierr);
    
    BitRefLevel bit_level_interface2;
    bit_level_interface1.set(2);
    ierr = mField.get_msId_3dENTS_split_sides(0,bit_level_interface2,meshset_interface2,true,true); CHKERRQ(ierr);
    EntityHandle meshset_level_interface2;
    rval = moab.create_meshset(MESHSET_SET,meshset_level_interface2); CHKERR_PETSC(rval);
    ierr = mField.refine_get_ents(bit_level_interface2,meshset_level_interface2); CHKERRQ(ierr);

    ierr = mField.refine_get_childern(meshset_SideSet1, bit_level_interface2,meshset_SideSet1,MBTRI,true); CHKERRQ(ierr);
    ierr = mField.refine_get_childern(meshset_SideSet2, bit_level_interface2,meshset_SideSet2,MBTRI,true); CHKERRQ(ierr);
//    ierr = mField.refine_get_childern(meshset_BlockSet1,bit_level_interface2,meshset_BlockSet1,MBTET,true); CHKERRQ(ierr);
//    ierr = mField.refine_get_childern(meshset_BlockSet2,bit_level_interface2,meshset_BlockSet2,MBTET,true); CHKERRQ(ierr);
//    ierr = mField.refine_get_childern(meshset_BlockSet3,bit_level_interface2,meshset_BlockSet3,MBTET,true); CHKERRQ(ierr);
    
//    rval = moab.write_file("AAA.vtk","VTK","",&meshset_level_interface,1); CHKERR_PETSC(rval);
    
    // stl::bitset see for more details
    BitRefLevel bit_level0;
    bit_level0.set(3);
    EntityHandle meshset_level0;
    rval = moab.create_meshset(MESHSET_SET,meshset_level0); CHKERR_PETSC(rval);
    ierr = mField.seed_ref_level_3D(meshset_level_interface2,bit_level0); CHKERRQ(ierr);
    ierr = mField.refine_get_ents(bit_level0,meshset_level0); CHKERRQ(ierr);
    
//    Range TETSFormLevel0;
//    rval = moab.get_entities_by_type(meshset_level0,MBTET,TETSFormLevel0,true); CHKERR_PETSC(rval);
//    Range TETsFromBlockSet1onmeshset_level0;
//    rval = moab.get_entities_by_type(meshset_BlockSet1,MBTET,TETsFromBlockSet1onmeshset_level0,true); CHKERR_PETSC(rval);
//    TETsFromBlockSet1onmeshset_level0 = intersect(TETsFromBlockSet1onmeshset_level0,TETSFormLevel0);
//    EntityHandle meshset_BlockSet1OnLevel0;
//    rval = moab.create_meshset(MESHSET_SET,meshset_BlockSet1OnLevel0); CHKERR_PETSC(rval); CHKERR_PETSC(rval);
//    rval = moab.add_entities(meshset_BlockSet1OnLevel0,TETsFromBlockSet1onmeshset_level0); CHKERR_PETSC(rval);
//    
//    Range TETsFromBlockSet2onmeshset_level0;
//    rval = moab.get_entities_by_type(meshset_BlockSet2,MBTET,TETsFromBlockSet2onmeshset_level0,true); CHKERR_PETSC(rval);
//    EntityHandle meshset_BlockSet2OnLevel0;
//    rval = moab.create_meshset(MESHSET_SET,meshset_BlockSet2OnLevel0); CHKERR_PETSC(rval); CHKERR_PETSC(rval);
//    TETsFromBlockSet2onmeshset_level0 = intersect(TETsFromBlockSet2onmeshset_level0,TETSFormLevel0);
//    rval = moab.add_entities(meshset_BlockSet2OnLevel0,TETsFromBlockSet2onmeshset_level0); CHKERR_PETSC(rval);
    
//    Range TETsFromBlockSet3onmeshset_level0;
//    rval = moab.get_entities_by_type(meshset_BlockSet3,MBTET,TETsFromBlockSet3onmeshset_level0,true); CHKERR_PETSC(rval);
//    EntityHandle meshset_BlockSet3OnLevel0;
//    rval = moab.create_meshset(MESHSET_SET,meshset_BlockSet3OnLevel0); CHKERR_PETSC(rval); CHKERR_PETSC(rval);
//    TETsFromBlockSet3onmeshset_level0 = intersect(TETsFromBlockSet3onmeshset_level0,TETSFormLevel0);
//    rval = moab.add_entities(meshset_BlockSet3OnLevel0,TETsFromBlockSet3onmeshset_level0); CHKERR_PETSC(rval);

    //rval = moab.write_file("AAA.vtk","VTK","",&meshset_BlockSet1,1); CHKERR_PETSC(rval);
    
    BitRefLevel problem_bit_level = bit_level0;
    
    /***/
    //Define problem
    
    //Fields
    ierr = mField.add_field("DISPLACEMENT",H1,3); CHKERRQ(ierr);
    
    //FE
    ierr = mField.add_finite_element("ELASTIC"); CHKERRQ(ierr);
    ierr = mField.add_finite_element("INTERFACE"); CHKERRQ(ierr);

    //Define rows/cols and element data
    ierr = mField.modify_finite_element_add_field_row("ELASTIC","DISPLACEMENT"); CHKERRQ(ierr);
    ierr = mField.modify_finite_element_add_field_col("ELASTIC","DISPLACEMENT"); CHKERRQ(ierr);
    ierr = mField.modify_finite_element_add_field_data("ELASTIC","DISPLACEMENT"); CHKERRQ(ierr);
    //FE Interface
    ierr = mField.modify_finite_element_add_field_row("INTERFACE","DISPLACEMENT"); CHKERRQ(ierr);
    ierr = mField.modify_finite_element_add_field_col("INTERFACE","DISPLACEMENT"); CHKERRQ(ierr);
    ierr = mField.modify_finite_element_add_field_data("INTERFACE","DISPLACEMENT"); CHKERRQ(ierr);
    
    //define problems
    ierr = mField.add_problem("ELASTIC_MECHANICS"); CHKERRQ(ierr);
    
    //set finite elements for problem
    ierr = mField.modify_problem_add_finite_element("ELASTIC_MECHANICS","ELASTIC"); CHKERRQ(ierr);
    ierr = mField.modify_problem_add_finite_element("ELASTIC_MECHANICS","INTERFACE"); CHKERRQ(ierr);

    
    //set refinment level for problem
    ierr = mField.modify_problem_ref_level_add_bit("ELASTIC_MECHANICS",problem_bit_level); CHKERRQ(ierr);
    
    /***/
    //Declare problem
    
    //add entitities (by tets) to the field
    ierr = mField.add_ents_to_field_by_TETs(0,"DISPLACEMENT"); CHKERRQ(ierr);

    //add finite elements entities
//    ierr = mField.add_ents_to_finite_element_by_TETs(meshset_BlockSet2OnLevel0,"ELASTIC",true); CHKERRQ(ierr);
//    ierr = mField.add_ents_to_finite_element_by_TETs(meshset_BlockSet1OnLevel0,"ELASTIC",true); CHKERRQ(ierr);
//    ierr = mField.add_ents_to_finite_element_by_TETs(meshset_BlockSet3OnLevel0,"TRAN_ISOTROPIC_ELASTIC",true); CHKERRQ(ierr);
    ierr = mField.add_ents_to_finite_element_EntType_by_bit_ref(problem_bit_level,"ELASTIC",MBTET); CHKERRQ(ierr);
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
    Vec F;
    ierr = mField.VecCreateGhost("ELASTIC_MECHANICS",Row,&F); CHKERRQ(ierr);
    Mat Aij;
    ierr = mField.MatCreateMPIAIJWithArrays("ELASTIC_MECHANICS",&Aij); CHKERRQ(ierr);
    
    //Get SideSet 1 and SideSet 2 defined in CUBIT
    Range SideSet1,SideSet2,SideSet3,SideSet4,SideSet5,SideSet6;
    ierr = mField.get_Cubit_msId_entities_by_dimension(1,SideSet,2,SideSet1,true); CHKERRQ(ierr);
    ierr = mField.get_Cubit_msId_entities_by_dimension(2,SideSet,2,SideSet2,true); CHKERRQ(ierr);
    ierr = mField.get_Cubit_msId_entities_by_dimension(3,SideSet,2,SideSet3,true); CHKERRQ(ierr);
    ierr = mField.get_Cubit_msId_entities_by_dimension(4,SideSet,2,SideSet4,true); CHKERRQ(ierr);
    ierr = mField.get_Cubit_msId_entities_by_dimension(5,SideSet,2,SideSet5,true); CHKERRQ(ierr);
    ierr = mField.get_Cubit_msId_entities_by_dimension(6,SideSet,2,SideSet6,true); CHKERRQ(ierr);

    PetscPrintf(PETSC_COMM_WORLD,"Nb. faces in SideSet 1 : %u\n",SideSet1.size());
    PetscPrintf(PETSC_COMM_WORLD,"Nb. faces in SideSet 2 : %u\n",SideSet2.size());
    PetscPrintf(PETSC_COMM_WORLD,"Nb. faces in SideSet 3 : %u\n",SideSet3.size());
    PetscPrintf(PETSC_COMM_WORLD,"Nb. faces in SideSet 4 : %u\n",SideSet4.size());
    PetscPrintf(PETSC_COMM_WORLD,"Nb. faces in SideSet 5 : %u\n",SideSet5.size());
    PetscPrintf(PETSC_COMM_WORLD,"Nb. faces in SideSet 6 : %u\n",SideSet6.size());
    
    //Assemble F and Aij
    const double YoungModulus = 200000;
    const double PoissonRatio = 0.3;
    const double alpha = 0.05;
            
    struct MyElasticFEMethod: public ElasticFEMethod {
        MyElasticFEMethod(Interface& _moab,BaseDirihletBC *_dirihlet_ptr,
                          Mat &_Aij,Vec& _F,
                          double _lambda,double _mu,Range &_SideSet1,Range &_SideSet2): 
        ElasticFEMethod(_moab,_dirihlet_ptr,_Aij,_F,_lambda,_mu,
                        _SideSet1,_SideSet2) {};
        
        /// Set Neumann Boundary Conditions on SideSet2
        PetscErrorCode NeumannBC() {
            PetscFunctionBegin;
            ublas::vector<FieldData,ublas::bounded_array<double,3> > traction(3);
            //Set Direction of Traction On SideSet2
            traction[0] = 0; //X
            traction[1] = 0; //Y 
            traction[2] = +100; //Z
            //ElasticFEMethod::NeumannBC(...) function calulating external forces (see file ElasticFEMethod.hpp)
            ierr = ElasticFEMethod::NeumannBC(F,traction,SideSet2); CHKERRQ(ierr);
            PetscFunctionReturn(0);
        }
        
        PetscErrorCode operator()() {
            PetscFunctionBegin;
            ierr = OpStudentStart_TET(g_NTET); CHKERRQ(ierr);
            ierr = GetMatrices(); CHKERRQ(ierr);
            
            //Dirihlet Boundary Condition
            ierr = SetDirihletBC_to_ElementIndicies(); CHKERRQ(ierr);
            if(Diagonal!=PETSC_NULL) {
                if(DirihletBC.size()>0) {
                    DirihletBCDiagVal.resize(DirihletBC.size());
                    fill(DirihletBCDiagVal.begin(),DirihletBCDiagVal.end(),1);
                    ierr = VecSetValues(Diagonal,DirihletBC.size(),&(DirihletBC[0]),&DirihletBCDiagVal[0],INSERT_VALUES); CHKERRQ(ierr);
                }
            }
            
            //Assembly Aij and F
            ierr = RhsAndLhs(); CHKERRQ(ierr);
            
            //Neumann Boundary Conditions
            ierr = NeumannBC(); CHKERRQ(ierr);
            
            ierr = OpStudentEnd(); CHKERRQ(ierr);
            PetscFunctionReturn(0);
        }
        
    };

    //**************************** Function to Set Boundary Conditions ******************************//  

    struct MyDirihletBC: public BaseDirihletBC {
        Range& SideSet1;
        Range SideSet1_;
        
        // Constructor
        MyDirihletBC(Interface &moab,Range& _SideSet1): 
        BaseDirihletBC(),SideSet1(_SideSet1){
            
            //Add to SideSet1_ nodes,edges, and faces, where dirihilet boundary conditions are applied.
            //Note that SideSet1 consist only faces in this particular example.
            ErrorCode rval;
            Range SideSet1Edges,SideSet1Nodes;
            
            rval = moab.get_adjacencies(SideSet1,1,false,SideSet1Edges,Interface::UNION); CHKERR_THROW(rval);
            rval = moab.get_connectivity(SideSet1,SideSet1Nodes,true); CHKERR_THROW(rval);
            SideSet1_.insert(SideSet1.begin(),SideSet1.end());
            SideSet1_.insert(SideSet1Edges.begin(),SideSet1Edges.end());
            SideSet1_.insert(SideSet1Nodes.begin(),SideSet1Nodes.end());
            
        }
        
        //This method is called insiaid finite element loop to apply boundary conditions
        PetscErrorCode SetDirihletBC_to_ElementIndicies(
                                                        moabField::FEMethod *fe_method_ptr,string field_name,
                                                        vector<vector<DofIdx> > &RowGlob,vector<vector<DofIdx> > &ColGlob,vector<DofIdx>& DirihletBC) {
            PetscFunctionBegin;
            
            ierr = InternalClassBCSet(fe_method_ptr,RowGlob,ColGlob,DirihletBC,field_name,SideSet1_,fixed_x|fixed_y|fixed_z,false); CHKERRQ(ierr);
            
            PetscFunctionReturn(0);
            
            
        }
        
    private:
        //Only to use in this auxiliary function
        enum bc_type { fixed_x = 1,fixed_y = 1<<1, fixed_z = 1<<2 };
        PetscErrorCode InternalClassBCSet(
                                          moabField::FEMethod *fe_method_ptr,vector<vector<DofIdx> > &RowGlob,vector<vector<DofIdx> > &ColGlob,vector<DofIdx>& DirihletBC,
                                          string field_name,Range &SideSet,unsigned int bc = fixed_x|fixed_y|fixed_z,bool zero_bc = true) {
            PetscFunctionBegin;
            //Dirihlet form SideSet1
            if(zero_bc) DirihletBC.resize(0);
            Range::iterator siit1 = SideSet.begin();
            for(;siit1!=SideSet.end();siit1++) {
                FENumeredDofMoFEMEntity_multiIndex::index<MoABEnt_mi_tag>::type::iterator riit = fe_method_ptr->row_multiIndex->get<MoABEnt_mi_tag>().lower_bound(*siit1);
                FENumeredDofMoFEMEntity_multiIndex::index<MoABEnt_mi_tag>::type::iterator hi_riit = fe_method_ptr->row_multiIndex->get<MoABEnt_mi_tag>().upper_bound(*siit1);
                for(;riit!=hi_riit;riit++) {
                    if(riit->get_name()!=field_name) continue;
                    unsigned int my_bc = 0;
                    switch(riit->get_dof_rank()) {
                        case 0: my_bc = fixed_x; break;
                        case 1: my_bc = fixed_y; break;
                        case 2: my_bc = fixed_z; break;
                        default:
                            SETERRQ(PETSC_COMM_SELF,1,"not implemented");
                    }
                    if((bc&my_bc) == 0) continue;
                    // all fixed
                    // if some ranks are selected then we could apply BC in particular direction
                    DirihletBC.push_back(riit->get_petsc_gloabl_dof_idx());
                    for(unsigned int cc = 0;cc<ColGlob.size();cc++) {
                        vector<DofIdx>::iterator it = find(ColGlob[cc].begin(),ColGlob[cc].end(),riit->get_petsc_gloabl_dof_idx());
                        if( it!=ColGlob[cc].end() ) *it = -1; // of idx is set -1 column is not assembled
                    }
                    for(unsigned int rr = 0;rr<RowGlob.size();rr++) {
                        vector<DofIdx>::iterator it = find(RowGlob[rr].begin(),RowGlob[rr].end(),riit->get_petsc_gloabl_dof_idx());
                        if( it!=RowGlob[rr].end() ) *it = -1; // of idx is set -1 row is not assembled
                    }
                }
            }
            PetscFunctionReturn(0);
        }
        
    };
    
    MyDirihletBC myDirihletBC(moab,SideSet1);
    MyElasticFEMethod MyFE(moab,&myDirihletBC,Aij,F,LAMBDA(YoungModulus,PoissonRatio),MU(YoungModulus,PoissonRatio),SideSet1,SideSet2);
    Range dummy;
    InterfaceFEMethod IntMyFE(moab,&myDirihletBC,Aij,F,YoungModulus*alpha,SideSet1,SideSet2,dummy);

    ierr = VecZeroEntries(F); CHKERRQ(ierr);
    ierr = VecGhostUpdateBegin(F,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
    ierr = VecGhostUpdateEnd(F,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
    ierr = MatZeroEntries(Aij); CHKERRQ(ierr);
    
    ierr = mField.loop_finite_elements("ELASTIC_MECHANICS","INTERFACE",IntMyFE);  CHKERRQ(ierr);
    PetscSynchronizedFlush(PETSC_COMM_WORLD);
    
    ierr = mField.loop_finite_elements("ELASTIC_MECHANICS","ELASTIC",MyFE);  CHKERRQ(ierr);
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
    
    Vec D;
    ierr = VecDuplicate(F,&D); CHKERRQ(ierr);
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
        ierr = mField.problem_get_FE("ELASTIC_MECHANICS","INTERFACE",out_meshset); CHKERRQ(ierr);
        rval = moab.write_file(outName,"VTK","",&out_meshset,1); CHKERR_PETSC(rval);
        rval = moab.delete_entities(&out_meshset,1); CHKERR_PETSC(rval);
    }
    
    PostProcDisplacemenysAndStarinAndElasticLinearStressOnRefMesh fe_post_proc_method( moab,  LAMBDA(YoungModulus,PoissonRatio),MU(YoungModulus,PoissonRatio));
    ierr = mField.loop_finite_elements("ELASTIC_MECHANICS","ELASTIC",fe_post_proc_method);  CHKERRQ(ierr);
    
    PetscSynchronizedFlush(PETSC_COMM_WORLD);
    if(pcomm->rank()==0) {
        rval = fe_post_proc_method.moab_post_proc.write_file(outName2,"VTK",""); CHKERR_PETSC(rval);
    }

    PostProcCohesiveForces fe_post_proc_prisms(moab,YoungModulus*alpha);
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
