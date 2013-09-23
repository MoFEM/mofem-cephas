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
#include "PostProcVertexMethod.hpp"

using namespace MoFEM;

ErrorCode rval;
PetscErrorCode ierr;

static char help[] = "-my_file _mesh_name_ (I can read h5m, cub, vtk, etc.) \n\n";

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

  // stl::bitset see for more details
  BitRefLevel bit_level0;
  bit_level0.set(0);
  EntityHandle meshset_level0;
  rval = moab.create_meshset(MESHSET_SET,meshset_level0); CHKERR_PETSC(rval);
  ierr = mField.seed_ref_level_3D(0,bit_level0); CHKERRQ(ierr);
  ierr = mField.refine_get_ents(bit_level0,meshset_level0); CHKERRQ(ierr);

  /***/
  //Define problem

  //Fields
  ierr = mField.add_field("DISPLACEMENT",H1,3); CHKERRQ(ierr);
  ierr = mField.add_field("TEMPERATURE",H1,1); CHKERRQ(ierr);

  //FE
  ierr = mField.add_finite_element("ELASTIC"); CHKERRQ(ierr);

  //Define rows/cols and element data
  ierr = mField.modify_finite_element_add_field_row("ELASTIC","DISPLACEMENT"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_col("ELASTIC","DISPLACEMENT"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_data("ELASTIC","DISPLACEMENT"); CHKERRQ(ierr);
    
  ierr = mField.modify_finite_element_add_field_data("ELASTIC","TEMPERATURE"); CHKERRQ(ierr);

  //define problems
  ierr = mField.add_problem("ELASTIC_MECHANICS"); CHKERRQ(ierr);

  //set finite elements for problem
  ierr = mField.modify_problem_add_finite_element("ELASTIC_MECHANICS","ELASTIC"); CHKERRQ(ierr);

  //set refinment level for problem
  ierr = mField.modify_problem_ref_level_add_bit("ELASTIC_MECHANICS",bit_level0); CHKERRQ(ierr);

  /***/
  //Declare problem

  //add entitities (by tets) to the field
  ierr = mField.add_ents_to_field_by_TETs(0,"DISPLACEMENT"); CHKERRQ(ierr);
  ierr = mField.add_ents_to_field_by_TETs(0,"TEMPERATURE"); CHKERRQ(ierr);

  //add finite elements entities
  ierr = mField.add_ents_to_finite_element_EntType_by_bit_ref(bit_level0,"ELASTIC",MBTET); CHKERRQ(ierr);

  //set app. order
  //see Hierarchic Finite Element Bases on Unstructured Tetrahedral Meshes (Mark Ainsworth & Joe Coyle)
  ierr = mField.set_field_order(0,MBTET,"DISPLACEMENT",1); CHKERRQ(ierr);
  ierr = mField.set_field_order(0,MBTRI,"DISPLACEMENT",1); CHKERRQ(ierr);
  ierr = mField.set_field_order(0,MBEDGE,"DISPLACEMENT",1); CHKERRQ(ierr);
  ierr = mField.set_field_order(0,MBVERTEX,"DISPLACEMENT",1); CHKERRQ(ierr);
    
    ierr = mField.set_field_order(0,MBTET,"TEMPERATURE",1); CHKERRQ(ierr);
    ierr = mField.set_field_order(0,MBTRI,"TEMPERATURE",1); CHKERRQ(ierr);
    ierr = mField.set_field_order(0,MBEDGE,"TEMPERATURE",1); CHKERRQ(ierr);
    ierr = mField.set_field_order(0,MBVERTEX,"TEMPERATURE",1); CHKERRQ(ierr);

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

    
    
//============================= New iterator to to get initial temprature and display it on the paraviewo mesh  =========================
    Tag th_temperature;
    double defaultval=0;
    rval = moab.tag_get_handle("TEMPERATURE",1,MB_TYPE_DOUBLE,th_temperature,MB_TAG_CREAT|MB_TAG_SPARSE,&defaultval); CHKERR_THROW(rval);
    EntityHandle node = 0;
    double coords[3];
    for(_IT_GET_DOFS_MOABFIELD_BY_NAME_FOR_LOOP_(mField,"TEMPERATURE",dof_ptr)) {
        if(dof_ptr->get_ent_type()!=MBVERTEX) continue;

        EntityHandle ent = dof_ptr->get_ent();
        double &fval = dof_ptr->get_FieldData();
        double tempcoord;
        if(node!=ent) {
            rval = moab.get_coords(&ent,1,coords); CHKERR_PETSC(rval);
            tempcoord = coords[0];
            fval = tempcoord;
            rval=moab.tag_set_data(th_temperature,&ent,1,&tempcoord); CHKERR_PETSC(rval);
            node = ent;
        }
    }

    
    FILE *fout_ptr;
    fout_ptr = fopen ("IT_GET_DOFS_MOABFIELD_BY_NAME_FOR_LOOP.txt", "w");
    fprintf(fout_ptr, "=====================     Display Data of the DOFs     =============================\n");
    
    for(_IT_GET_DOFS_MOABFIELD_BY_NAME_FOR_LOOP_(mField,"TEMPERATURE",dof_ptr)) {
            ostringstream ss;
            ss << *dof_ptr;
            PetscSynchronizedFPrintf(PETSC_COMM_WORLD, fout_ptr, "%s\n\n",ss.str().c_str());
    }
    fclose (fout_ptr);

    
    
    

    
    
    
    
    
    
    
    
    
    

  /****/
  //mesh partitioning 

  //partition
  ierr = mField.partition_problem("ELASTIC_MECHANICS"); CHKERRQ(ierr);
  ierr = mField.partition_finite_elements("ELASTIC_MECHANICS"); CHKERRQ(ierr);
  //what are ghost nodes, see Petsc Manual
  ierr = mField.partition_ghost_dofs("ELASTIC_MECHANICS"); CHKERRQ(ierr);

    struct MyFiniteElement_For_Printing: public moabField::FEMethod {
        MyFiniteElement_For_Printing() {};
        
        PetscErrorCode preProcess() {
            PetscFunctionBegin;
            PetscPrintf(PETSC_COMM_WORLD,"Start Loop\n");
            
            PetscFunctionReturn(0);
        }
        PetscErrorCode operator()() {
            PetscFunctionBegin;
            
//            cout<<"\n\n";
//            cout<<"FE name"<<fe_name.c_str(); 
//            cout<<"\n\n";
            
                 //if(fe_ptr->get_ent_id()==1){
//                cout<<"\n\n";
//                PetscSynchronizedPrintf(PETSC_COMM_WORLD,"FE name %s proc %u\n",fe_name.c_str(),fe_ptr->get_part());
//                PetscSynchronizedPrintf(PETSC_COMM_WORLD,"Entity ID %u\n",fe_ptr->get_ent_id());
//                PetscSynchronizedPrintf(PETSC_COMM_WORLD,"Row DOFs %u  Col DOFs %u  Data %u\n",fe_ptr->get_nb_dofs_row(), fe_ptr->get_nb_dofs_col(), fe_ptr->get_nb_dofs_data());
//                cout<<"\n\n";
//                ostringstream ss;
//                ss << *fe_ptr;
//                PetscSynchronizedPrintf(PETSC_COMM_WORLD,"%s",ss.str().c_str()  );
//                cout<<"\n\n";
//                
//                
//                                
//                cout<<"\n\n==============Contents of container row_multiIndex==================\n";
//                for(FENumeredDofMoFEMEntity_multiIndex::iterator it=row_multiIndex->begin();it!=row_multiIndex->end();++it){
//                    cout<<"Unique ID ="<<it->get_unique_id()<<"\tEntity ="<<it->get_ent()<<"\tName ="<<it->get_name()<<endl;
//                    //need to understand the use of the composite key here. 
//                    //pair<FENumeredDofMoFEMEntity_multiIndex::iterator,FENumeredDofMoFEMEntity_multiIndex::iterator> it1;                    
//                }
//
//                cout<<"\n\n==============Contents of container col_multiIndex==================\n";
//                for(FENumeredDofMoFEMEntity_multiIndex::iterator it=col_multiIndex->begin();it!=col_multiIndex->end();++it){
//                    cout<<"Unique ID ="<<it->get_unique_id()<<"\tEntity ="<<it->get_ent()<<"\tName ="<<it->get_name()<<endl;
//                }
//                
//                cout<<"\n\n==============Contents of container data_multiIndex==================\n";
//                for(FEDofMoFEMEntity_multiIndex::iterator it=data_multiIndex->begin();it!=data_multiIndex->end();++it){
//                    cout<<"Unique ID ="<<it->get_unique_id()<<"\tEntity ="<<it->get_ent()<<"\tName ="<<it->get_name()<<endl;
//                }
//
//                cout<<"\n\n==============Contents of container moabfields==================\n";
//                for(MoFEMField_multiIndex::iterator it=moabfields->begin();it!=moabfields->end();++it){
//                    cout<<"ID ="<<it->get_id()<<"\tmeshset ="<<it->meshset<<"\tName ="<<it->get_name()<<"\tSpace ="<<it->get_space()<<endl;
//                }
 
                
            //}
                
                
                
            //cout<<"\n\n\n\n";
            
            
//            //cout<<"\n\n\n\n\n\n";
//            for(_IT_GET_FEROW_BY_SIDE_DOFS_FOR_LOOP_(this,"DISPLACEMENT",MBVERTEX,1,it)) {
//                ostringstream ss;
//                ss << it->get_ent();
//                
//                cerr << ss << endl;
//                PetscSynchronizedPrintf(PETSC_COMM_WORLD, "%s\n\n",ss.str().c_str());
//            }
//            //cout<<"\n\n\n\n\n\n";
            
            

//            cout<<"\n\n\n\n\n";
//            for(_IT_GET_FEROW_BY_TYPE_DOFS_FOR_LOOP_(this,"TEMPERATURE",MBTET, dof_ptr)) {
//                
//                cout<<"-------------"<<endl;
//                cout<<dof_ptr-> get_FieldData();
//                //ostringstream ss;
//                //ss << *dof_ptr;
//                //PetscSynchronizedPrintf(PETSC_COMM_WORLD, "%s\n\n",ss.str().c_str());
//            }
//            cout<<"\n\n\n\n\n";
//
            
    
            
            
            //ostringstream ss;
            //ss << *fe_ptr;
            //PetscSynchronizedPrintf(PETSC_COMM_WORLD,"%s",ss.str().c_str());
            
            
            PetscFunctionReturn(0);
        }
        PetscErrorCode postProcess() {
            PetscFunctionBegin;
            PetscSynchronizedFlush(PETSC_COMM_WORLD);
            PetscPrintf(PETSC_COMM_WORLD,"End Loop\n");
            PetscFunctionReturn(0);
        }

    };

    MyFiniteElement_For_Printing myFE;
    ierr = mField.loop_finite_elements ("ELASTIC_MECHANICS","ELASTIC",myFE); CHKERRQ(ierr);
    
//Save data on mesh

    PostProcVertexMethod ent_method(moab);
    ierr = mField.loop_dofs("ELASTIC_MECHANICS","DISPLACEMENT",Row,ent_method); CHKERRQ(ierr);

    if(pcomm->rank()==0) {
        EntityHandle out_meshset;
        rval = moab.create_meshset(MESHSET_SET,out_meshset); CHKERR_PETSC(rval);
        ierr = mField.problem_get_FE("ELASTIC_MECHANICS","ELASTIC",out_meshset); CHKERRQ(ierr);
        rval = moab.write_file("out.vtk","VTK","",&out_meshset,1); CHKERR_PETSC(rval);
        rval = moab.delete_entities(&out_meshset,1); CHKERR_PETSC(rval);
    }
  
  ierr = PetscGetTime(&v2);CHKERRQ(ierr);
  ierr = PetscGetCPUTime(&t2);CHKERRQ(ierr);

  PetscSynchronizedPrintf(PETSC_COMM_WORLD,"Total Rank %d Time = %f CPU Time = %f\n",pcomm->rank(),v2-v1,t2-t1);
  PetscSynchronizedFlush(PETSC_COMM_WORLD);

  PetscFinalize();

}

