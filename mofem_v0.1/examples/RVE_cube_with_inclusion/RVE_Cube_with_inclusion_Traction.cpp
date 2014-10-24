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

#include "ElasticFEMethod.hpp"
#include "ElasticFE_RVELagrange_Traction.hpp"
#include "ElasticFE_RVELagrange_Homogenized_Stress_Traction.hpp"
#include "ElasticFE_RVELagrange_RigidBodyTranslation.hpp"
#include "ElasticFE_RVELagrange_RigidBodyRotation.hpp"
#include "RVEVolume.hpp"

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
  
  
  //=======================================================================================================
  //Seting nodal coordinates on the surface to make sure they are periodic
  //=======================================================================================================
  
  Range SurTrisNeg, SurTrisPos;
  ierr = mField.get_Cubit_msId_entities_by_dimension(101,SIDESET,2,SurTrisNeg,true); CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"number of SideSet 101 = %d\n",SurTrisNeg.size()); CHKERRQ(ierr);
  ierr = mField.get_Cubit_msId_entities_by_dimension(102,SIDESET,2,SurTrisPos,true); CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"number of SideSet 102 = %d\n",SurTrisPos.size()); CHKERRQ(ierr);
  
  Range SurNodesNeg,SurNodesPos;
  rval = moab.get_connectivity(SurTrisNeg,SurNodesNeg,true); CHKERR_PETSC(rval);
  cout<<" All nodes on negative surfaces " << SurNodesNeg.size()<<endl;
  rval = moab.get_connectivity(SurTrisPos,SurNodesPos,true); CHKERR_PETSC(rval);
  cout<<" All nodes on positive surfaces " << SurNodesPos.size()<<endl;
  
  
  double roundfact=1000.0;   double coords_nodes[3];
  //Populating the Multi-index container with nodes on -ve faces
  for(Range::iterator nit = SurNodesNeg.begin(); nit!=SurNodesNeg.end();  nit++) {
    rval = moab.get_coords(&*nit,1,coords_nodes);  CHKERR_PETSC(rval);
    //round values to 3 disimal places
    if(coords_nodes[0]>=0) coords_nodes[0]=double(int(coords_nodes[0]*roundfact+0.5))/roundfact;  else coords_nodes[0]=double(int(coords_nodes[0]*roundfact-0.5))/roundfact;
    if(coords_nodes[1]>=0) coords_nodes[1]=double(int(coords_nodes[1]*roundfact+0.5))/roundfact;  else coords_nodes[1]=double(int(coords_nodes[1]*roundfact-0.5))/roundfact;
    if(coords_nodes[2]>=0) coords_nodes[2]=double(int(coords_nodes[2]*roundfact+0.5))/roundfact;  else coords_nodes[2]=double(int(coords_nodes[2]*roundfact-0.5))/roundfact;
    rval = moab.set_coords(&*nit,1,coords_nodes);  CHKERR_PETSC(rval);
    //      cout<<"   coords_nodes[0]= "<<coords_nodes[0] << "   coords_nodes[1]= "<< coords_nodes[1] << "   coords_nodes[2]= "<< coords_nodes[2] <<endl;
  }
  
  ///Populating the Multi-index container with nodes on +ve faces
  for(Range::iterator nit = SurNodesPos.begin(); nit!=SurNodesPos.end();  nit++) {
    rval = moab.get_coords(&*nit,1,coords_nodes);  CHKERR_PETSC(rval);
    //round values to 3 disimal places
    if(coords_nodes[0]>=0) coords_nodes[0]=double(int(coords_nodes[0]*roundfact+0.5))/roundfact;  else coords_nodes[0]=double(int(coords_nodes[0]*roundfact-0.5))/roundfact;
    if(coords_nodes[1]>=0) coords_nodes[1]=double(int(coords_nodes[1]*roundfact+0.5))/roundfact;  else coords_nodes[1]=double(int(coords_nodes[1]*roundfact-0.5))/roundfact;
    if(coords_nodes[2]>=0) coords_nodes[2]=double(int(coords_nodes[2]*roundfact+0.5))/roundfact;  else coords_nodes[2]=double(int(coords_nodes[2]*roundfact-0.5))/roundfact;
    rval = moab.set_coords(&*nit,1,coords_nodes);  CHKERR_PETSC(rval);
    //      cout<<"   coords_nodes[0]= "<<coords_nodes[0] << "   coords_nodes[1]= "<< coords_nodes[1] << "   coords_nodes[2]= "<< coords_nodes[2] <<endl;
  }
  //=======================================================================================================

  
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
  int field_rank=3;
  ierr = mField.add_field("DISPLACEMENT",H1,field_rank); CHKERRQ(ierr);
  ierr = mField.add_field("Lagrange_mul_disp",NOFIELD,6); CHKERRQ(ierr);
  ierr = mField.add_field("Lagrange_mul_disp_rigid_trans",NOFIELD,3); CHKERRQ(ierr);   //Control 3 rigid body translations in x, y and z axis
  ierr = mField.add_field("Lagrange_mul_disp_rigid_rotation",NOFIELD,3); CHKERRQ(ierr); //Controla 3 rigid body rotations about x, y and z axis
  
    
    //FE
    ierr = mField.add_finite_element("ELASTIC"); CHKERRQ(ierr);
    ierr = mField.add_finite_element("ELASTIC_Inc"); CHKERRQ(ierr);
    ierr = mField.add_finite_element("Lagrange_elem"); CHKERRQ(ierr);
    ierr = mField.add_finite_element("Lagrange_elem_rigid_trans"); CHKERRQ(ierr);
    ierr = mField.add_finite_element("Lagrange_elem_rigid_rotation"); CHKERRQ(ierr);
    
    
    //Define rows/cols and element data
    ierr = mField.modify_finite_element_add_field_row("ELASTIC","DISPLACEMENT"); CHKERRQ(ierr);
    ierr = mField.modify_finite_element_add_field_col("ELASTIC","DISPLACEMENT"); CHKERRQ(ierr);
    ierr = mField.modify_finite_element_add_field_data("ELASTIC","DISPLACEMENT"); CHKERRQ(ierr);
    
    //Define rows/cols and element data for ELASTIC_Inclusion
    ierr = mField.modify_finite_element_add_field_row("ELASTIC_Inc","DISPLACEMENT"); CHKERRQ(ierr);
    ierr = mField.modify_finite_element_add_field_col("ELASTIC_Inc","DISPLACEMENT"); CHKERRQ(ierr);
    ierr = mField.modify_finite_element_add_field_data("ELASTIC_Inc","DISPLACEMENT"); CHKERRQ(ierr);
    
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
    
    
    //Define rows/cols and element data for C1 and C1T (for lagrange multipliers to contol the rigid body translations)
    //============================================================================================================
    //C1 row as Lagrange_elem_rigid_trans and col as DISPLACEMENT
    ierr = mField.modify_finite_element_add_field_row("Lagrange_elem_rigid_trans","Lagrange_mul_disp_rigid_trans"); CHKERRQ(ierr);
    ierr = mField.modify_finite_element_add_field_col("Lagrange_elem_rigid_trans","DISPLACEMENT"); CHKERRQ(ierr);
    
    //C1T col as Lagrange_elem_rigid_trans and row as DISPLACEMENT
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
    ierr = mField.modify_problem_add_finite_element("ELASTIC_MECHANICS","ELASTIC_Inc"); CHKERRQ(ierr);
    ierr = mField.modify_problem_add_finite_element("ELASTIC_MECHANICS","ELASTIC"); CHKERRQ(ierr);
    ierr = mField.modify_problem_add_finite_element("ELASTIC_MECHANICS","Lagrange_elem"); CHKERRQ(ierr);
    ierr = mField.modify_problem_add_finite_element("ELASTIC_MECHANICS","Lagrange_elem_rigid_trans"); CHKERRQ(ierr);
    ierr = mField.modify_problem_add_finite_element("ELASTIC_MECHANICS","Lagrange_elem_rigid_rotation"); CHKERRQ(ierr);
    
    
    //set refinment level for problem
    ierr = mField.modify_problem_ref_level_add_bit("ELASTIC_MECHANICS",bit_level0); CHKERRQ(ierr);
    
    /***/
    //Declare problem
    
    //add entitities (by tets) to the field
    ierr = mField.add_ents_to_field_by_TETs(0,"DISPLACEMENT",2); CHKERRQ(ierr);
    
    
    //add finite elements entities
    EntityHandle meshset_Elastic, meshset_Elastic_stiff_inclusion;
    rval = moab.create_meshset(MESHSET_SET,meshset_Elastic); CHKERR_PETSC(rval);
    rval = moab.create_meshset(MESHSET_SET,meshset_Elastic_stiff_inclusion); CHKERR_PETSC(rval);
    
    Range TetsInBlock;
    for(_IT_CUBITMESHSETS_BY_SET_TYPE_FOR_LOOP_(mField,BLOCKSET,it)){
        
		if(it->get_Cubit_name() == "MAT_ELASTIC") {
			rval = moab.get_entities_by_type(it->meshset, MBTET,TetsInBlock,true); CHKERR_PETSC(rval);
		}
	}
    //    cout<<"TetsInBlock "<<TetsInBlock<<endl;
    rval = moab.add_entities(meshset_Elastic,TetsInBlock);CHKERR_PETSC(rval);
    
    Range TetsInBlock_stiff_inc;
    for(_IT_CUBITMESHSETS_BY_SET_TYPE_FOR_LOOP_(mField,BLOCKSET,it)){
		if(it->get_Cubit_name() == "MAT_ELASTIC_Stiff_Inclusion") {
			rval = moab.get_entities_by_type(it->meshset, MBTET,TetsInBlock_stiff_inc,true); CHKERR_PETSC(rval);
		}
	}
    //    cout<<"TetsInBlock "<<TetsInBlock_stiff_inc<<endl;
    rval = moab.add_entities(meshset_Elastic_stiff_inclusion,TetsInBlock_stiff_inc);CHKERR_PETSC(rval);
    
    ierr = mField.add_ents_to_finite_element_by_TETs(meshset_Elastic,"ELASTIC",true); CHKERRQ(ierr);
    ierr = mField.add_ents_to_finite_element_by_TETs(meshset_Elastic_stiff_inclusion,"ELASTIC_Inc",true); CHKERRQ(ierr);

    
    //add finite elements entities
    Range SurfacesFaces;
    ierr = mField.get_Cubit_msId_entities_by_dimension(103,SIDESET,2,SurfacesFaces,true); CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_WORLD,"number of SideSet 103 = %d\n",SurfacesFaces.size()); CHKERRQ(ierr);
  
    //to create meshset from range
    EntityHandle BoundFacesMeshset;
    rval = moab.create_meshset(MESHSET_SET,BoundFacesMeshset); CHKERR_PETSC(rval);
    rval = moab.add_entities(BoundFacesMeshset,SurfacesFaces); CHKERR_PETSC(rval);
    ierr = mField.seed_ref_level_MESHSET(BoundFacesMeshset,BitRefLevel().set()); CHKERRQ(ierr);

    ierr = mField.add_ents_to_finite_element_by_TRIs(SurfacesFaces,"Lagrange_elem"); CHKERRQ(ierr);
    ierr = mField.add_ents_to_finite_element_by_TRIs(SurfacesFaces,"Lagrange_elem_rigid_trans"); CHKERRQ(ierr);
    ierr = mField.add_ents_to_finite_element_by_TRIs(SurfacesFaces,"Lagrange_elem_rigid_rotation"); CHKERRQ(ierr);
    
    
    //see Hierarchic Finite Element Bases on Unstructured Tetrahedral Meshes (Mark Ainsworth & Joe Coyle)
    //int order = 5;
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
    ierr = mField.print_cubit_displacement_set(); CHKERRQ(ierr);
    ierr = mField.print_cubit_force_set(); CHKERRQ(ierr);
    //print block sets with materials
    ierr = mField.print_cubit_materials_set(); CHKERRQ(ierr);
    
    //    for(_IT_GET_DOFS_FIELD_BY_NAME_AND_TYPE_FOR_LOOP_(mField,"Lagrange_mul_disp",MBEDGE,dof)) {
    //        cerr << *dof << endl;
    //
    //    }
    
    //create matrices (here F, D and Aij are matrices for the full problem)
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
    
  
    //Assemble F and Aij
    const double YoungModulus = 1;
    const double PoissonRatio = 0.0;
    MyElasticFEMethod MyFE(mField,Aij,D,F,LAMBDA(YoungModulus,PoissonRatio),MU(YoungModulus,PoissonRatio));
  ElasticFE_RVELagrange_Traction MyFE_RVELagrange(mField,Aij,D,F,applied_strain,"DISPLACEMENT","Lagrange_mul_disp",field_rank);
  
  ElasticFE_RVELagrange_RigidBodyTranslation MyFE_RVELagrangeRigidBodyTrans(mField,Aij,D,F,applied_strain,"DISPLACEMENT","Lagrange_mul_disp",field_rank,"Lagrange_mul_disp_rigid_trans");
  ElasticFE_RVELagrange_RigidBodyRotation MyFE_RVELagrangeRigidBodyRotation(mField,Aij,D,F,applied_strain,"DISPLACEMENT","Lagrange_mul_disp",field_rank);
    
    ierr = VecZeroEntries(F); CHKERRQ(ierr);
    ierr = VecGhostUpdateBegin(F,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
    ierr = VecGhostUpdateEnd(F,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
    ierr = MatZeroEntries(Aij); CHKERRQ(ierr);
    
    ierr = mField.loop_finite_elements("ELASTIC_MECHANICS","ELASTIC",MyFE);  CHKERRQ(ierr);
    ierr = mField.loop_finite_elements("ELASTIC_MECHANICS","ELASTIC_Inc",MyFE);  CHKERRQ(ierr);
    ierr = mField.loop_finite_elements("ELASTIC_MECHANICS","Lagrange_elem",MyFE_RVELagrange);  CHKERRQ(ierr);
    ierr = mField.loop_finite_elements("ELASTIC_MECHANICS","Lagrange_elem_rigid_trans",MyFE_RVELagrangeRigidBodyTrans);  CHKERRQ(ierr);
    ierr = mField.loop_finite_elements("ELASTIC_MECHANICS","Lagrange_elem_rigid_rotation",MyFE_RVELagrangeRigidBodyRotation);  CHKERRQ(ierr);
  
    ierr = VecGhostUpdateBegin(F,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
    ierr = VecGhostUpdateEnd(F,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
    ierr = VecAssemblyBegin(F); CHKERRQ(ierr);
    ierr = VecAssemblyEnd(F); CHKERRQ(ierr);
    ierr = MatAssemblyBegin(Aij,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
    ierr = MatAssemblyEnd(Aij,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
    
    PetscSynchronizedFlush(PETSC_COMM_WORLD);
    //ierr = VecView(F,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);
    //ierr = MatView(Aij,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);
    
    
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
    //ierr = VecView(F,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);
    
    
    //Calculation of Homogenized stress
    //============================================================================================================================================
    double RVE_volume;    RVE_volume=0.0;  //RVE volume for full RVE We need this for stress calculation
    Vec RVE_volume_Vec;
    ierr = VecCreateMPI(PETSC_COMM_WORLD, 1, pcomm->size(), &RVE_volume_Vec);  CHKERRQ(ierr);
    ierr = VecZeroEntries(RVE_volume_Vec); CHKERRQ(ierr);
    
    RVEVolume MyRVEVol(mField,Aij,D,F,LAMBDA(YoungModulus,PoissonRatio),MU(YoungModulus,PoissonRatio), RVE_volume_Vec);
    ierr = mField.loop_finite_elements("ELASTIC_MECHANICS","ELASTIC",MyRVEVol);  CHKERRQ(ierr);
    ierr = mField.loop_finite_elements("ELASTIC_MECHANICS","ELASTIC_Inc",MyRVEVol);  CHKERRQ(ierr);

    ierr = VecView(RVE_volume_Vec,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);
    ierr = VecSum(RVE_volume_Vec, &RVE_volume);  CHKERRQ(ierr);
    cout<<"Final RVE_volume = "<< RVE_volume <<endl;
    
    //create a vector for 6 components of homogenized stress
    Vec Stress_Homo;
    ierr = VecCreateMPI(PETSC_COMM_WORLD, 6, 6*pcomm->size(), &Stress_Homo);  CHKERRQ(ierr);
    ierr = VecZeroEntries(Stress_Homo); CHKERRQ(ierr);
    
    ////    ierr = VecView(D,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);
  ElasticFE_RVELagrange_Homogenized_Stress_Traction MyFE_RVEHomoStressTraction(mField,Aij,D,F,&RVE_volume,applied_strain, Stress_Homo,"DISPLACEMENT","Lagrange_mul_disp",field_rank);
    ierr = mField.loop_finite_elements("ELASTIC_MECHANICS","Lagrange_elem",MyFE_RVEHomoStressTraction);  CHKERRQ(ierr);
    
//    if(pcomm->rank()) cout<< " Stress_Homo =  "<<endl;
//    ierr = VecView(Stress_Homo,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);
  
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

  
    //=============================================================================================================================================

    PostProcVertexMethod ent_method(moab);
    ierr = mField.loop_dofs("ELASTIC_MECHANICS","DISPLACEMENT",ROW,ent_method); CHKERRQ(ierr);
    
    if(pcomm->rank()==0) {
        EntityHandle out_meshset;
        rval = moab.create_meshset(MESHSET_SET,out_meshset); CHKERR_PETSC(rval);
        ierr = mField.problem_get_FE("ELASTIC_MECHANICS","ELASTIC",out_meshset); CHKERRQ(ierr);
        ierr = mField.problem_get_FE("ELASTIC_MECHANICS","ELASTIC_Inc",out_meshset); CHKERRQ(ierr);
        rval = moab.write_file("out.vtk","VTK","",&out_meshset,1); CHKERR_PETSC(rval);
        rval = moab.delete_entities(&out_meshset,1); CHKERR_PETSC(rval);
    }
    
    //PostProcDisplacemenysAndStarinOnRefMesh fe_post_proc_method(moab);
    PostProcDisplacemenysAndStarinAndElasticLinearStressOnRefMesh fe_post_proc_method(mField,"DISPLACEMENT",LAMBDA(YoungModulus,PoissonRatio),MU(YoungModulus,PoissonRatio));
    ierr = mField.loop_finite_elements("ELASTIC_MECHANICS","ELASTIC",fe_post_proc_method);  CHKERRQ(ierr);
    
    if(pcomm->rank()==0) {
        rval = fe_post_proc_method.moab_post_proc.write_file("out_post_proc.vtk","VTK",""); CHKERR_PETSC(rval);
    }
    
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

