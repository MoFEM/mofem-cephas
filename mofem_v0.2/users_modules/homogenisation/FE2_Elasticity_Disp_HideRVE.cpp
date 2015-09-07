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

#include <MoFEM.hpp>
using namespace MoFEM;

#include <DirichletBC.hpp>

#include <Projection10NodeCoordsOnField.hpp>
#include <petsctime.h>

#include <SurfacePressure.hpp>

#include <FEMethod_LowLevelStudent.hpp>
#include <FEMethod_UpLevelStudent.hpp>

#include <PostProcVertexMethod.hpp>
#include <PostProcDisplacementAndStrainOnRefindedMesh.hpp>

#include <ElasticFEMethod.hpp>
#include "ElasticFE_RVELagrange_Disp.hpp"
#include <FE2_ElasticFEMethod.hpp>


#include "ElasticFE_RVELagrange_Disp_Multi_Rhs.hpp"
#include "ElasticFE_RVELagrange_Homogenized_Stress_Disp.hpp"
#include "RVEVolume.hpp"
#include "CalculateRVE_Dmat_Iso_Disp1.hpp"

using namespace boost::numeric;
using namespace ObosleteUsersModules;

ErrorCode rval;
PetscErrorCode ierr;

static char help[] = "...\n\n";

const double young_modulus = 1;
const double poisson_ratio = 0.0;


int main(int argc, char *argv[]) {
  
  
  //=============================================================================================================
  //  Micro (RVE) Problme
  //=============================================================================================================

  PetscInitialize(&argc,&argv,(char *)0,help);

  moab::Core mb_instance_RVE;
  Interface& moab_RVE = mb_instance_RVE;
  int rank;
  MPI_Comm_rank(PETSC_COMM_WORLD,&rank);
  
  //Reade parameters from line command
  PetscBool flg = PETSC_TRUE;
  char mesh_file_name[255];
  ierr = PetscOptionsGetString(PETSC_NULL,"-my_file_RVE",mesh_file_name,255,&flg); CHKERRQ(ierr);
  if(flg != PETSC_TRUE) {
    SETERRQ(PETSC_COMM_SELF,1,"*** ERROR -my_file_RVE (MESH FILE NEEDED)");
  }
  PetscInt order;
  ierr = PetscOptionsGetInt(PETSC_NULL,"-my_order",&order,&flg); CHKERRQ(ierr);
  if(flg != PETSC_TRUE) {
    order = 1;
  }

  //Read mesh to MOAB
  const char *option;
  option = "";//"PARALLEL=BCAST;";//;DEBUG_IO";
  rval = moab_RVE.load_file(mesh_file_name, 0, option); CHKERR_PETSC(rval);
  ParallelComm* pcomm = ParallelComm::get_pcomm(&moab_RVE,MYPCOMM_INDEX);
  if(pcomm == NULL) pcomm =  new ParallelComm(&moab_RVE,PETSC_COMM_WORLD);

  
  //Create MoFEM (Joseph) database
  MoFEM::Core core_RVE(moab_RVE);
  FieldInterface& mField_RVE = core_RVE;
  
  //ref meshset ref level 0
  ierr = mField_RVE.seed_ref_level_3D(0,0); CHKERRQ(ierr);
  
  // stl::bitset see for more details
  BitRefLevel bit_level0;
  bit_level0.set(0);
  EntityHandle meshset_level0;
  rval = moab_RVE.create_meshset(MESHSET_SET,meshset_level0); CHKERR_PETSC(rval);
  ierr = mField_RVE.seed_ref_level_3D(0,bit_level0); CHKERRQ(ierr);
  ierr = mField_RVE.get_entities_by_ref_level(bit_level0,BitRefLevel().set(),meshset_level0); CHKERRQ(ierr);
  

  /***/
  //Define problem
  
  //Fields
  int field_rank=3;
  ierr = mField_RVE.add_field("DISPLACEMENT",H1,field_rank); CHKERRQ(ierr);
  ierr = mField_RVE.add_field("Lagrange_mul_disp",H1,field_rank); CHKERRQ(ierr);
  
  //FE
  ierr = mField_RVE.add_finite_element("ELASTIC"); CHKERRQ(ierr);
  ierr = mField_RVE.add_finite_element("Lagrange_elem"); CHKERRQ(ierr);
  
  //Define rows/cols and element data
  ierr = mField_RVE.modify_finite_element_add_field_row("ELASTIC","DISPLACEMENT"); CHKERRQ(ierr);
  ierr = mField_RVE.modify_finite_element_add_field_col("ELASTIC","DISPLACEMENT"); CHKERRQ(ierr);
  ierr = mField_RVE.modify_finite_element_add_field_data("ELASTIC","DISPLACEMENT"); CHKERRQ(ierr);
  
  //C and CT
  //======================================================================================================
  //C row as Lagrange_mul_disp and col as DISPLACEMENT
  ierr = mField_RVE.modify_finite_element_add_field_row("Lagrange_elem","Lagrange_mul_disp"); CHKERRQ(ierr);
  ierr = mField_RVE.modify_finite_element_add_field_col("Lagrange_elem","DISPLACEMENT"); CHKERRQ(ierr);
  
  //CT col as Lagrange_mul_disp and row as DISPLACEMENT
  ierr = mField_RVE.modify_finite_element_add_field_col("Lagrange_elem","Lagrange_mul_disp"); CHKERRQ(ierr);
  ierr = mField_RVE.modify_finite_element_add_field_row("Lagrange_elem","DISPLACEMENT"); CHKERRQ(ierr);
  
  //data
  ierr = mField_RVE.modify_finite_element_add_field_data("Lagrange_elem","Lagrange_mul_disp"); CHKERRQ(ierr);
  ierr = mField_RVE.modify_finite_element_add_field_data("Lagrange_elem","DISPLACEMENT"); CHKERRQ(ierr);
  //======================================================================================================

  
  //define problems
  ierr = mField_RVE.add_problem("ELASTIC_MECHANICS"); CHKERRQ(ierr);
  
  
  //set finite elements for problem
  ierr = mField_RVE.modify_problem_add_finite_element("ELASTIC_MECHANICS","ELASTIC"); CHKERRQ(ierr);
  ierr = mField_RVE.modify_problem_add_finite_element("ELASTIC_MECHANICS","Lagrange_elem"); CHKERRQ(ierr);
  
  
  //set refinment level for problem
  ierr = mField_RVE.modify_problem_ref_level_add_bit("ELASTIC_MECHANICS",bit_level0); CHKERRQ(ierr);
  
  /***/
  //Declare problem
  
  //add entitities (by tets) to the field
  ierr = mField_RVE.add_ents_to_field_by_TETs(0,"DISPLACEMENT"); CHKERRQ(ierr);

  
  //add finite elements entities
  ierr = mField_RVE.add_ents_to_finite_element_EntType_by_bit_ref(bit_level0,"ELASTIC",MBTET); CHKERRQ(ierr);
  Range SurfacesFaces;
  ierr = mField_RVE.get_cubit_msId_entities_by_dimension(103,SIDESET,2,SurfacesFaces,true); CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"number of SideSet 103 = %d\n",SurfacesFaces.size()); CHKERRQ(ierr);
  ierr = mField_RVE.add_ents_to_finite_element_by_TRIs(SurfacesFaces,"Lagrange_elem"); CHKERRQ(ierr);
  

  //to create meshset from range
  EntityHandle BoundFacesMeshset;
  rval = moab_RVE.create_meshset(MESHSET_SET,BoundFacesMeshset); CHKERR_PETSC(rval);
	rval = moab_RVE.add_entities(BoundFacesMeshset,SurfacesFaces); CHKERR_PETSC(rval);
  ierr = mField_RVE.seed_ref_level_MESHSET(BoundFacesMeshset,BitRefLevel().set()); CHKERRQ(ierr);
  ierr = mField_RVE.add_ents_to_field_by_TRIs(BoundFacesMeshset,"Lagrange_mul_disp",2); CHKERRQ(ierr);
  
  //set app. order
  //see Hierarchic Finite Element Bases on Unstructured Tetrahedral Meshes (Mark Ainsworth & Joe Coyle)
  //int order = 5;
  ierr = mField_RVE.set_field_order(0,MBTET,"DISPLACEMENT",order); CHKERRQ(ierr);
  ierr = mField_RVE.set_field_order(0,MBTRI,"DISPLACEMENT",order); CHKERRQ(ierr);
  ierr = mField_RVE.set_field_order(0,MBEDGE,"DISPLACEMENT",order); CHKERRQ(ierr);
  ierr = mField_RVE.set_field_order(0,MBVERTEX,"DISPLACEMENT",1); CHKERRQ(ierr);
  
  ierr = mField_RVE.set_field_order(0,MBTRI,"Lagrange_mul_disp",order); CHKERRQ(ierr);
  ierr = mField_RVE.set_field_order(0,MBEDGE,"Lagrange_mul_disp",order); CHKERRQ(ierr);
  ierr = mField_RVE.set_field_order(0,MBVERTEX,"Lagrange_mul_disp",1); CHKERRQ(ierr);
  
  /****/
  //build database
  
  //build field
  ierr = mField_RVE.build_fields(); CHKERRQ(ierr);
  
  //build finite elemnts
  ierr = mField_RVE.build_finite_elements(); CHKERRQ(ierr);
  
  //build adjacencies
  ierr = mField_RVE.build_adjacencies(bit_level0); CHKERRQ(ierr);
  
  //build problem
  ierr = mField_RVE.build_problems(); CHKERRQ(ierr);

  
  /****/
  //mesh partitioning
  
  //partition
  ierr = mField_RVE.partition_problem("ELASTIC_MECHANICS"); CHKERRQ(ierr);
  ierr = mField_RVE.partition_finite_elements("ELASTIC_MECHANICS"); CHKERRQ(ierr);
  //what are ghost nodes, see Petsc Manual
  ierr = mField_RVE.partition_ghost_dofs("ELASTIC_MECHANICS"); CHKERRQ(ierr);
  
  //print bcs
  ierr = mField_RVE.print_cubit_displacement_set(); CHKERRQ(ierr);
  ierr = mField_RVE.print_cubit_force_set(); CHKERRQ(ierr);
  //print block sets with materials
  ierr = mField_RVE.print_cubit_materials_set(); CHKERRQ(ierr);

  
  
  ////=============================================================================================================
  ////=============================================================================================================
  //Here i hide all the RVE code

  Calculate_RVE_Dmat_iso_Disp1 calculate_rve_dmat_iso_disp;
  calculate_rve_dmat_iso_disp.Calculate_RVEDmat(mField_RVE);
  cout <<"\n\n\n\n\n Dmat in main program "<<calculate_rve_dmat_iso_disp.Dmat<<endl<<endl;
  
  ////=============================================================================================================
  ////=============================================================================================================
  
  
//=============================================================================================================
//  Macro Problme
//=============================================================================================================

  moab::Core mb_instance_Macro;
  Interface& moab_Macro = mb_instance_Macro;
  ParallelComm* pcomm_Macro = ParallelComm::get_pcomm(&moab_Macro,MYPCOMM_INDEX);
  if(pcomm == NULL) pcomm_Macro =  new ParallelComm(&moab_Macro,PETSC_COMM_WORLD);

  //Reade parameters from line command
  ierr = PetscOptionsGetString(PETSC_NULL,"-my_file_macro",mesh_file_name,255,&flg); CHKERRQ(ierr);
  if(flg != PETSC_TRUE) {
    SETERRQ(PETSC_COMM_SELF,1,"*** ERROR -my_file_macro (MESH FILE NEEDED)");
  }
  
  //Read mesh to MOAB
  //option = "PARALLEL=BCAST_DELETE;"
      //"PARTITION=GEOM_DIMENSION,PARTITION_VAL=3,PARTITION_DISTRIBUTE";//;DEBUG_IO";
  rval = moab_Macro.load_file(mesh_file_name, 0, option); CHKERR_PETSC(rval);

  //Create MoFEM (Joseph) database
  MoFEM::Core core_Macro(moab_Macro);
  FieldInterface& mField_Macro = core_Macro;

  //ref meshset ref level 0
  ierr = mField_Macro.seed_ref_level_3D(0,0); CHKERRQ(ierr);

  // stl::bitset see for more details
  BitRefLevel bit_level0_Macro;
  bit_level0_Macro.set(0);
  EntityHandle meshset_level0_Macro;
  rval = moab_Macro.create_meshset(MESHSET_SET,meshset_level0_Macro); CHKERR_PETSC(rval);
  ierr = mField_Macro.seed_ref_level_3D(0,bit_level0_Macro); CHKERRQ(ierr);
  ierr = mField_Macro.get_entities_by_ref_level(bit_level0_Macro,BitRefLevel().set(),meshset_level0_Macro); CHKERRQ(ierr);

  //Define problem

  //Fields
  ierr = mField_Macro.add_field("DISP_MACORO",H1,3,MF_ZERO); CHKERRQ(ierr);
  ierr = mField_Macro.add_field("MESH_NODE_POSITIONS",H1,3,MF_ZERO); CHKERRQ(ierr);

  //FE
  ierr = mField_Macro.add_finite_element("ELASTIC_MACRO",MF_ZERO); CHKERRQ(ierr);

  //Define rows/cols and element data
  ierr = mField_Macro.modify_finite_element_add_field_row("ELASTIC_MACRO","DISP_MACORO"); CHKERRQ(ierr);
  ierr = mField_Macro.modify_finite_element_add_field_col("ELASTIC_MACRO","DISP_MACORO"); CHKERRQ(ierr);
  ierr = mField_Macro.modify_finite_element_add_field_data("ELASTIC_MACRO","DISP_MACORO"); CHKERRQ(ierr);
  ierr = mField_Macro.modify_finite_element_add_field_data("ELASTIC_MACRO","MESH_NODE_POSITIONS"); CHKERRQ(ierr);

  //define problems
  ierr = mField_Macro.add_problem("ELASTIC_PROB"); CHKERRQ(ierr);

  //set finite elements for problem
  ierr = mField_Macro.modify_problem_add_finite_element("ELASTIC_PROB","ELASTIC_MACRO"); CHKERRQ(ierr);

  //set refinment level for problem
  ierr = mField_Macro.modify_problem_ref_level_add_bit("ELASTIC_PROB",bit_level0_Macro); CHKERRQ(ierr);

  //Declare problem

  //add entitities (by tets) to the field
  ierr = mField_Macro.add_ents_to_field_by_TETs(0,"DISP_MACORO"); CHKERRQ(ierr);
  ierr = mField_Macro.add_ents_to_field_by_TETs(0,"MESH_NODE_POSITIONS"); CHKERRQ(ierr);

  //add finite elements entities
  ierr = mField_Macro.add_ents_to_finite_element_EntType_by_bit_ref(bit_level0_Macro,"ELASTIC_MACRO",MBTET); CHKERRQ(ierr);

  //set app. order
  //see Hierarchic Finite Element Bases on Unstructured Tetrahedral Meshes (Mark Ainsworth & Joe Coyle)
  //int order = 5;
  ierr = mField_Macro.set_field_order(0,MBTET,"DISP_MACORO",order); CHKERRQ(ierr);
  ierr = mField_Macro.set_field_order(0,MBTRI,"DISP_MACORO",order); CHKERRQ(ierr);
  ierr = mField_Macro.set_field_order(0,MBEDGE,"DISP_MACORO",order); CHKERRQ(ierr);
  ierr = mField_Macro.set_field_order(0,MBVERTEX,"DISP_MACORO",1); CHKERRQ(ierr);
  //
  ierr = mField_Macro.set_field_order(0,MBTET,"MESH_NODE_POSITIONS",2); CHKERRQ(ierr);
  ierr = mField_Macro.set_field_order(0,MBTRI,"MESH_NODE_POSITIONS",2); CHKERRQ(ierr);
  ierr = mField_Macro.set_field_order(0,MBEDGE,"MESH_NODE_POSITIONS",2); CHKERRQ(ierr);
  ierr = mField_Macro.set_field_order(0,MBVERTEX,"MESH_NODE_POSITIONS",1); CHKERRQ(ierr);

  ierr = MetaNeummanForces::addNeumannBCElements(mField_Macro,"DISP_MACORO"); CHKERRQ(ierr);
  ierr = mField_Macro.modify_problem_add_finite_element("ELASTIC_PROB","FORCE_FE"); CHKERRQ(ierr);
  
  //build database

  //build field
  ierr = mField_Macro.build_fields(); CHKERRQ(ierr);

  //build finite elemnts
  ierr = mField_Macro.build_finite_elements(); CHKERRQ(ierr);

  //build adjacencies
  ierr = mField_Macro.build_adjacencies(bit_level0); CHKERRQ(ierr);

  //build problem
  ierr = mField_Macro.build_problems(); CHKERRQ(ierr);

  //mesh partitioning 

  //partition
  ierr = mField_Macro.partition_problem("ELASTIC_PROB"); CHKERRQ(ierr);
  //PetscBarrier(PETSC_NULL);
  ierr = mField_Macro.partition_finite_elements("ELASTIC_PROB"); CHKERRQ(ierr);
  //what are ghost nodes, see Petsc Manual
  ierr = mField_Macro.partition_ghost_dofs("ELASTIC_PROB"); CHKERRQ(ierr);

  //mField_Macro.list_dofs_by_field_name("DISP_MACORO",true);

  //print bcs
  ierr = mField_Macro.print_cubit_displacement_set(); CHKERRQ(ierr);
  ierr = mField_Macro.print_cubit_force_set(); CHKERRQ(ierr);
  //print block sets with materials
  ierr = mField_Macro.print_cubit_materials_set(); CHKERRQ(ierr);

  //create matrices
  Vec F,D;
  ierr = mField_Macro.VecCreateGhost("ELASTIC_PROB",ROW,&F); CHKERRQ(ierr);
  ierr = mField_Macro.VecCreateGhost("ELASTIC_PROB",COL,&D); CHKERRQ(ierr);

  Mat A;
  ierr = mField_Macro.MatCreateMPIAIJWithArrays("ELASTIC_PROB",&A); CHKERRQ(ierr);

  //Matrix View
  //MatView(A,PETSC_VIEWER_DRAW_WORLD);//PETSC_VIEWER_STDOUT_WORLD);
  //std::string wait;
  //std::cin >> wait;


  struct MyElasticFEMethod_Macro: public FE2_ElasticFEMethod {
    MyElasticFEMethod_Macro(FieldInterface& _mField_Macro,Mat _A,Vec _D,Vec& _F, ublas::matrix<FieldData> _Dmat,string _field_name):
    FE2_ElasticFEMethod(_mField_Macro,_A,_D,_F, _Dmat, _field_name) {};
    
    PetscErrorCode Fint(Vec F_int) {
      PetscFunctionBegin;
      ierr = FE2_ElasticFEMethod::Fint(); CHKERRQ(ierr);
      for(int rr = 0;rr<row_mat;rr++) {
        if(RowGlob[rr].size()!=f_int[rr].size()) SETERRQ(PETSC_COMM_SELF,1,"data inconsistency");
        if(RowGlob[rr].size()==0) continue;
        f_int[rr] *= -1; //This is not SNES we solve K*D = -RES
        ierr = VecSetValues(F_int,RowGlob[rr].size(),&(RowGlob[rr])[0],&(f_int[rr].data()[0]),ADD_VALUES); CHKERRQ(ierr);
      }
      PetscFunctionReturn(0);
    }
    
  };
  
  
  Projection10NodeCoordsOnField ent_method_material_Macro(mField_Macro,"MESH_NODE_POSITIONS");
  ierr = mField_Macro.loop_dofs("MESH_NODE_POSITIONS",ent_method_material_Macro); CHKERRQ(ierr);

  //Assemble F and A
  DisplacementBCFEMethodPreAndPostProc my_dirichlet_bc(mField_Macro,"DISP_MACORO",A,D,F);
  MyElasticFEMethod_Macro my_fe_Macro(mField_Macro,A,D,F,calculate_rve_dmat_iso_disp.Dmat,"DISP_MACORO");

  ierr = VecZeroEntries(F); CHKERRQ(ierr);
  ierr = VecGhostUpdateBegin(F,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = VecGhostUpdateEnd(F,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = MatZeroEntries(A); CHKERRQ(ierr);
  
  //preproc
  ierr = mField_Macro.problem_basic_method_preProcess("ELASTIC_PROB",my_dirichlet_bc); CHKERRQ(ierr);
  //loop elems
  //PetscBarrier(PETSC_NULL);
  ierr = mField_Macro.loop_finite_elements("ELASTIC_PROB","ELASTIC_MACRO",my_fe_Macro);  CHKERRQ(ierr);

  //forces and preassures on surface
  boost::ptr_map<string,NeummanForcesSurface> neumann_forces;
  ierr = MetaNeummanForces::setNeumannFiniteElementOperators(mField_Macro,neumann_forces,F,"DISP_MACORO"); CHKERRQ(ierr);
  boost::ptr_map<string,NeummanForcesSurface>::iterator mit = neumann_forces.begin();
  for(;mit!=neumann_forces.end();mit++) {
    ierr = mField_Macro.loop_finite_elements("ELASTIC_PROB",mit->first,mit->second->getLoopFe()); CHKERRQ(ierr);
  }

  //postproc
  ierr = mField_Macro.problem_basic_method_postProcess("ELASTIC_PROB",my_dirichlet_bc); CHKERRQ(ierr);

  //set matrix possitives define and symetric for cholesky and icc preceonditionser
  ierr = MatSetOption(A,MAT_SPD,PETSC_TRUE); CHKERRQ(ierr);

  ierr = VecGhostUpdateBegin(F,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
  ierr = VecGhostUpdateEnd(F,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
  ierr = VecAssemblyBegin(F); CHKERRQ(ierr);
  ierr = VecAssemblyEnd(F); CHKERRQ(ierr);

  //PetscSynchronizedFlush(PETSC_COMM_WORLD,PETSC_STDOUT);

  //Solver
  KSP solver_Macro;
  ierr = KSPCreate(PETSC_COMM_WORLD,&solver_Macro); CHKERRQ(ierr);
  ierr = KSPSetOperators(solver_Macro,A,A); CHKERRQ(ierr);
  ierr = KSPSetFromOptions(solver_Macro); CHKERRQ(ierr);
  ierr = KSPSetUp(solver_Macro); CHKERRQ(ierr);

//  MatView(A,PETSC_VIEWER_DRAW_WORLD);//PETSC_VIEWER_STDOUT_WORLD);

  // elastic analys
  ierr = KSPSolve(solver_Macro,F,D); CHKERRQ(ierr);
  ierr = VecGhostUpdateBegin(D,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = VecGhostUpdateEnd(D,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  
//  ierr = VecView(D,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);

  //Save data on mesh
  ierr = mField_Macro.set_global_ghost_vector("ELASTIC_PROB",ROW,D,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
  
  ProjectionFieldOn10NodeTet ent_method_on_10nodeTet(mField_Macro,"DISP_MACORO",true,false,"DISP_MACORO");
  ierr = mField_Macro.loop_dofs("DISP_MACORO",ent_method_on_10nodeTet); CHKERRQ(ierr);
  ent_method_on_10nodeTet.set_nodes = false;
  ierr = mField_Macro.loop_dofs("DISP_MACORO",ent_method_on_10nodeTet); CHKERRQ(ierr);
  
  if(pcomm->rank()==0) {
    EntityHandle out_meshset;
    rval = moab_Macro.create_meshset(MESHSET_SET,out_meshset); CHKERR_PETSC(rval);
    ierr = mField_Macro.get_problem_finite_elements_entities("ELASTIC_PROB","ELASTIC_MACRO",out_meshset); CHKERRQ(ierr);
    rval = moab_Macro.write_file("out.vtk","VTK","",&out_meshset,1); CHKERR_PETSC(rval);
    rval = moab_Macro.delete_entities(&out_meshset,1); CHKERR_PETSC(rval);
  }

  
  //Destroy matrices
  ierr = VecDestroy(&F); CHKERRQ(ierr);
  ierr = VecDestroy(&D); CHKERRQ(ierr);
  ierr = MatDestroy(&A); CHKERRQ(ierr);
  ierr = KSPDestroy(&solver_Macro); CHKERRQ(ierr);
  PetscFinalize();

}

