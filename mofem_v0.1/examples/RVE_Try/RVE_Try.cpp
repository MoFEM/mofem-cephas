/* Copyright (C) 2013, Zahur Ullah (Zahur.Ullah AT glasgow.ac.uk)
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

#include "ElasticFEMethod_RVE_Try.hpp"
#include "ElasticFEMethod_RVE_CalStress.hpp"
#include "PostProcVertexMethod.hpp"
#include "PostProcDisplacementAndStrainOnRefindedMesh.hpp"

using namespace MoFEM;

ErrorCode rval;
PetscErrorCode ierr;

static char help[] = "...\n\n";

int main(int argc, char *argv[]) {

  PetscInitialize(&argc,&argv,(char *)0,help);

  Core mb_instance;
  Interface& moab = mb_instance;
  int rank, num_processors;
  MPI_Comm_rank(PETSC_COMM_WORLD,&rank);
  MPI_Comm_size(PETSC_COMM_WORLD,&num_processors);
    
  //cout<<"/n/n Number of processors "<<num_processors<<endl<<endl;

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
  ierr = mField.refine_get_ents(bit_level0,BitRefLevel().set(),meshset_level0); CHKERRQ(ierr);

  /***/
  //Define problem

  //Fields
  ierr = mField.add_field("DISPLACEMENT",H1,3); CHKERRQ(ierr);
  ierr = mField.add_field("DISPLACEMENT_FROM_APP_STRAIN",H1,3); CHKERRQ(ierr);

  //FE
  ierr = mField.add_finite_element("ELASTIC"); CHKERRQ(ierr);

  //Define rows/cols and element data
  ierr = mField.modify_finite_element_add_field_row("ELASTIC","DISPLACEMENT"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_col("ELASTIC","DISPLACEMENT"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_data("ELASTIC","DISPLACEMENT"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_data("ELASTIC","DISPLACEMENT_FROM_APP_STRAIN"); CHKERRQ(ierr);

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
  ierr = mField.add_ents_to_field_by_TETs(0,"DISPLACEMENT_FROM_APP_STRAIN"); CHKERRQ(ierr);

  //add finite elements entities
  ierr = mField.add_ents_to_finite_element_EntType_by_bit_ref(bit_level0,"ELASTIC",MBTET); CHKERRQ(ierr);

  //set app. order
  //see Hierarchic Finite Element Bases on Unstructured Tetrahedral Meshes (Mark Ainsworth & Joe Coyle)
  //int order = 5;
  ierr = mField.set_field_order(0,MBTET,"DISPLACEMENT",order); CHKERRQ(ierr);
  ierr = mField.set_field_order(0,MBTRI,"DISPLACEMENT",order); CHKERRQ(ierr);
  ierr = mField.set_field_order(0,MBEDGE,"DISPLACEMENT",order); CHKERRQ(ierr);
  ierr = mField.set_field_order(0,MBVERTEX,"DISPLACEMENT",1); CHKERRQ(ierr);

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
  ierr = mField.build_adjacencies(bit_level0); CHKERRQ(ierr);

  //build problem
  ierr = mField.build_problems(); CHKERRQ(ierr);

    
    
   /****/
  //mesh partitioning 

  //partition
  ierr = mField.partition_problem("ELASTIC_MECHANICS"); CHKERRQ(ierr);
  ierr = mField.partition_finite_elements("ELASTIC_MECHANICS"); CHKERRQ(ierr);
  //what are ghost nodes, see Petsc Manual
  ierr = mField.partition_ghost_dofs("ELASTIC_MECHANICS"); CHKERRQ(ierr);

  //print bcs
  ierr = mField.printCubitDisplacementSet(); CHKERRQ(ierr);
  ierr = mField.printCubitForceSet(); CHKERRQ(ierr);
  //print block sets with materials
  ierr = mField.printCubitMaterials(); CHKERRQ(ierr);

  //create matrices
  Vec F,D,D_star,F_stress, coord_stress, F_stress_col;
  ierr = mField.VecCreateGhost("ELASTIC_MECHANICS",Row,&F); CHKERRQ(ierr);
  ierr = mField.VecCreateGhost("ELASTIC_MECHANICS",Col,&D); CHKERRQ(ierr);
  ierr = mField.VecCreateGhost("ELASTIC_MECHANICS",Col,&D_star); CHKERRQ(ierr);
  ierr = mField.VecCreateGhost("ELASTIC_MECHANICS",Row,&F_stress); CHKERRQ(ierr);
  ierr = mField.VecCreateGhost("ELASTIC_MECHANICS",Row,&coord_stress); CHKERRQ(ierr);
  ierr = mField.VecCreateGhost("ELASTIC_MECHANICS",Row,&F_stress_col); CHKERRQ(ierr);
    
   Mat Aij;
   ierr = mField.MatCreateMPIAIJWithArrays("ELASTIC_MECHANICS",&Aij); CHKERRQ(ierr);

   //Applied strain (specified by the user)
    ublas::matrix<double> strain_app;
    strain_app.resize(3,3);
    strain_app(0,0) = 0.01; strain_app(0,1)=0.0; strain_app(0,2)=0.0;
    strain_app(1,0) = 0.0;  strain_app(1,1)=0.0; strain_app(1,2)=0.0;
    strain_app(2,0) = 0.0;  strain_app(2,1)=0.0; strain_app(2,2)=0.0;
    
    //cout<<"\n\n strain_app "<<strain_app[0][0]<<"\n\n";
    //cout<<"\n\n strain_app "<<strain_app[1][1]<<"\n\n";
    
    
    //Apply the linear displacement to all nodes in the mesh
    Tag th_disp_1;
    double defaultval[3]={1,0,0};
    rval = moab.tag_get_handle("DISPLACEMENT_FROM_APP_STRAIN",3,MB_TYPE_DOUBLE,th_disp_1,MB_TAG_CREAT|MB_TAG_SPARSE,&defaultval); CHKERR_THROW(rval);  //_FROM_APP_STRAIN
    
    EntityHandle node = 0;
    double coords[3], disp_applied[3];
    int count=0;
    for(_IT_GET_DOFS_FIELD_BY_NAME_FOR_LOOP_(mField,"DISPLACEMENT_FROM_APP_STRAIN",dof_ptr)) {
        //cout<<"\n get_EntDofIdx "<<dof_ptr->get_EntDofIdx()<<"\n";
        //cout<<"\n dof_ptr->get_petsc_gloabl_dof_idx() "<<dof_ptr->get_petsc_gloabl_dof_idx()<<"\n";
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
        ierr = mField.refine_get_ents(bit_level0,BitRefLevel().set(),MBTET,out_meshset); CHKERRQ(ierr);
        rval = moab.write_file("Zahur_out_disp1.vtk","VTK","",&out_meshset,1); CHKERR_PETSC(rval);
        rval = moab.delete_entities(&out_meshset,1); CHKERR_PETSC(rval);
    }
    
    
  CubitDisplacementDirihletBC myDirihletBC(mField,"ELASTIC_MECHANICS","DISPLACEMENT");
  ierr = myDirihletBC.Init(); CHKERRQ(ierr);

  //Assemble F and Aij
  const double YoungModulus = 1;
  const double PoissonRatio = 0.0;
  ElasticFEMethod_RVE_Try MyFE(mField,&myDirihletBC,Aij,D,F,LAMBDA(YoungModulus,PoissonRatio),MU(YoungModulus,PoissonRatio),D_star,strain_app);
  
  ierr = VecZeroEntries(F); CHKERRQ(ierr);
  ierr = VecGhostUpdateBegin(F,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = VecGhostUpdateEnd(F,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = MatZeroEntries(Aij); CHKERRQ(ierr);
  ierr = VecZeroEntries(D_star); CHKERRQ(ierr);
  ierr = VecGhostUpdateBegin(D_star,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = VecGhostUpdateEnd(D_star,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);

  ierr = mField.loop_finite_elements("ELASTIC_MECHANICS","ELASTIC",MyFE);  CHKERRQ(ierr);
    
    
  ierr = VecGhostUpdateBegin(F,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
  ierr = VecGhostUpdateEnd(F,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
  ierr = VecAssemblyBegin(F); CHKERRQ(ierr);
  ierr = VecAssemblyEnd(F); CHKERRQ(ierr);
  ierr = MatAssemblyBegin(Aij,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
  ierr = MatAssemblyEnd(Aij,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
  PetscSynchronizedFlush(PETSC_COMM_WORLD);
  
//    ierr = VecView(D_star,  PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);
//    ierr = VecView(F,  PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);
 
  //Solver
  KSP solver;
  ierr = KSPCreate(PETSC_COMM_WORLD,&solver); CHKERRQ(ierr);
  ierr = KSPSetOperators(solver,Aij,Aij,SAME_NONZERO_PATTERN); CHKERRQ(ierr);
  ierr = KSPSetFromOptions(solver); CHKERRQ(ierr);
  ierr = KSPSetUp(solver); CHKERRQ(ierr);

  ierr = KSPSolve(solver,F,D); CHKERRQ(ierr);
  ierr = VecGhostUpdateBegin(D,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = VecGhostUpdateEnd(D,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);

  ierr = VecView(D, PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);

//    ierr = VecView(D_star,  PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);
    ierr = VecAXPY(D,1,D_star); CHKERRQ(ierr);   // D=D+1*D_star (total displacement start+fluctuation)
    ierr = VecView(D,  PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);
    
  //Save data on mesh
  ierr = mField.set_global_VecCreateGhost("ELASTIC_MECHANICS",Row,D,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
//  ierr = VecView(D,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);

    
  // calculate the homogenised stress
    
  ierr = VecZeroEntries(F_stress); CHKERRQ(ierr);
  ierr = VecZeroEntries(coord_stress); CHKERRQ(ierr);

   double RVE_volume;
   ublas::matrix<double> Sigma_homo;
    
   ElasticFEMethod_RVE_CalStress MyRVEStress(mField,&myDirihletBC,Aij,D,F,LAMBDA(YoungModulus,PoissonRatio),MU(YoungModulus,PoissonRatio), moab, F_stress, coord_stress, RVE_volume, Sigma_homo);
    ierr = mField.loop_finite_elements("ELASTIC_MECHANICS","ELASTIC",MyRVEStress);  CHKERRQ(ierr);
    
    
    PostProcVertexMethod ent_method(moab);
  ierr = mField.loop_dofs("ELASTIC_MECHANICS","DISPLACEMENT",Row,ent_method); CHKERRQ(ierr);
  if(pcomm->rank()==0) {
    EntityHandle out_meshset;
    rval = moab.create_meshset(MESHSET_SET,out_meshset); CHKERR_PETSC(rval);
    ierr = mField.problem_get_FE("ELASTIC_MECHANICS","ELASTIC",out_meshset); CHKERRQ(ierr);
    rval = moab.write_file("out.vtk","VTK","",&out_meshset,1); CHKERR_PETSC(rval);
    rval = moab.delete_entities(&out_meshset,1); CHKERR_PETSC(rval);
  }

    
//  PostProcDisplacemenysAndStarinOnRefMesh fe_post_proc_method(moab);
//  PostProcDisplacemenysAndStarinAndElasticLinearStressOnRefMesh fe_post_proc_method(moab,LAMBDA(YoungModulus,PoissonRatio),MU(YoungModulus,PoissonRatio));
//  ierr = mField.loop_finite_elements("ELASTIC_MECHANICS","ELASTIC",fe_post_proc_method);  CHKERRQ(ierr);
//
//  if(pcomm->rank()==0) {
//    rval = fe_post_proc_method.moab_post_proc.write_file("out_post_proc.vtk","VTK",""); CHKERR_PETSC(rval);
//  }

  //End Disp
  /*Range ents;*/
  //ierr = mField.get_Cubit_msId_entities_by_dimension(1,NodeSet,0,ents,true); CHKERRQ(ierr);
  //for(_IT_GET_DOFS_FIELD_BY_NAME_FOR_LOOP_(mField,"DISPLACEMENT",dit)) { 
    //if(find(ents.begin(),ents.end(),dit->get_ent())!=ents.end()) {
      //PetscSynchronizedPrintf(PETSC_COMM_WORLD, "val = %6.7e\n",dit->get_FieldData());
    //}
  /*}*/

  //Support stresses
  //Range ents;
  //ierr = mField.get_Cubit_msId_entities_by_dimension(4,NodeSet,0,ents,true); CHKERRQ(ierr);
      
  //Destroy matrices
    
    
//    //Matrix and Vector View
//    ierr = MatView(Aij,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);//PETSC_VIEWER_DRAW_WORLD);
//    ierr = VecView(D,  PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);
//    ierr = VecView(D_star,  PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);

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



