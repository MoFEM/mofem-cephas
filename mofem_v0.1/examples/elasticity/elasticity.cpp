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

#include "common.hpp"

#include "FieldInterface.hpp"
#include "FieldCore.hpp"

#include "FEMethod_UpLevelStudent.hpp"
#include "ElasticFEMethod.hpp"

#include "ForcesAndSurcesCore.hpp"
#include "SnesCtx.hpp"
#include "TsCtx.hpp"

#ifdef __cplusplus
extern "C" {
#endif
  #include<cblas.h>
  #include<lapack_wrap.h>
#ifdef __cplusplus
}
#endif

#include "SurfacePressure.hpp"
#include "NodalForce.hpp"
#include "FluidPressure.hpp"
#include "BodyForce.hpp"

#include "ThermalStressElement.hpp"

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/vector.hpp>

#include "Projection10NodeCoordsOnField.hpp"
#include "PostProcVertexMethod.hpp"
#include "PostProcDisplacementAndStrainOnRefindedMesh.hpp"

#include <petscksp.h>

using namespace boost::numeric;
using namespace MoFEM;

ErrorCode rval;
PetscErrorCode ierr;

static char help[] = "...\n\n";

const double young_modulus = 1;
const double poisson_ratio = 0.0;

PetscErrorCode write_soltion(FieldInterface &mField,const string out_file, const string out_ref_file) {
  PetscFunctionBegin;

  PostProcVertexMethod ent_method(mField.get_moab(),"DISPLACEMENT");
  ierr = mField.loop_dofs("ELASTIC_PROB","DISPLACEMENT",ROW,ent_method); CHKERRQ(ierr);
 
  ParallelComm* pcomm = ParallelComm::get_pcomm(&mField.get_moab(),MYPCOMM_INDEX);

  if(pcomm->rank()==0) {
    EntityHandle out_meshset;
    rval = mField.get_moab().create_meshset(MESHSET_SET,out_meshset); CHKERR_PETSC(rval);
    ierr = mField.problem_get_FE("ELASTIC_PROB","ELASTIC",out_meshset); CHKERRQ(ierr);
    rval = mField.get_moab().write_file(out_file.c_str(),"VTK","",&out_meshset,1); CHKERR_PETSC(rval);
    rval = mField.get_moab().delete_entities(&out_meshset,1); CHKERR_PETSC(rval);
  }
 
  if(pcomm->rank()==0) {

    PostProcDisplacemenysAndStarinAndElasticLinearStressOnRefMesh fe_post_proc_method(
      mField,"DISPLACEMENT",LAMBDA(young_modulus,poisson_ratio),MU(young_modulus,poisson_ratio));
    fe_post_proc_method.do_broadcast = false;
    ierr = mField.loop_finite_elements("ELASTIC_PROB","ELASTIC",fe_post_proc_method,0,pcomm->size());  CHKERRQ(ierr);
    rval = fe_post_proc_method.moab_post_proc.write_file(out_ref_file.c_str(),"VTK",""); CHKERR_PETSC(rval);
  }

  PetscFunctionReturn(0);
}

int main(int argc, char *argv[]) {

  PetscInitialize(&argc,&argv,(char *)0,help);

  Core mb_instance;
  Interface& moab = mb_instance;
  int rank;
  MPI_Comm_rank(PETSC_COMM_WORLD,&rank);
  ParallelComm* pcomm = ParallelComm::get_pcomm(&moab,MYPCOMM_INDEX);
  if(pcomm == NULL) pcomm =  new ParallelComm(&moab,PETSC_COMM_WORLD);

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
  //option = "PARALLEL=BCAST_DELETE;"
      //"PARTITION=GEOM_DIMENSION,PARTITION_VAL=3,PARTITION_DISTRIBUTE";//;DEBUG_IO";
  option = "";
  rval = moab.load_file(mesh_file_name, 0, option); CHKERR_PETSC(rval); 

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

  //Define problem

  //Fields
  ierr = mField.add_field("DISPLACEMENT",H1,3,MF_ZERO); CHKERRQ(ierr);
  ierr = mField.add_field("MESH_NODE_POSITIONS",H1,3,MF_ZERO); CHKERRQ(ierr);

  //FE
  ierr = mField.add_finite_element("ELASTIC",MF_ZERO); CHKERRQ(ierr);

  //Define rows/cols and element data
  ierr = mField.modify_finite_element_add_field_row("ELASTIC","DISPLACEMENT"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_col("ELASTIC","DISPLACEMENT"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_data("ELASTIC","DISPLACEMENT"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_data("ELASTIC","MESH_NODE_POSITIONS"); CHKERRQ(ierr);

  //define problems
  ierr = mField.add_problem("ELASTIC_PROB"); CHKERRQ(ierr);

  //set finite elements for problem
  ierr = mField.modify_problem_add_finite_element("ELASTIC_PROB","ELASTIC"); CHKERRQ(ierr);

  //set refinment level for problem
  ierr = mField.modify_problem_ref_level_add_bit("ELASTIC_PROB",bit_level0); CHKERRQ(ierr);

  //Declare problem

  //add entitities (by tets) to the field
  ierr = mField.add_ents_to_field_by_TETs(0,"DISPLACEMENT"); CHKERRQ(ierr);
  ierr = mField.add_ents_to_field_by_TETs(0,"MESH_NODE_POSITIONS"); CHKERRQ(ierr);

  //add finite elements entities
  ierr = mField.add_ents_to_finite_element_EntType_by_bit_ref(bit_level0,"ELASTIC",MBTET); CHKERRQ(ierr);

  //set app. order
  //see Hierarchic Finite Element Bases on Unstructured Tetrahedral Meshes (Mark Ainsworth & Joe Coyle)
  //int order = 5;
  ierr = mField.set_field_order(0,MBTET,"DISPLACEMENT",order); CHKERRQ(ierr);
  ierr = mField.set_field_order(0,MBTRI,"DISPLACEMENT",order); CHKERRQ(ierr);
  ierr = mField.set_field_order(0,MBEDGE,"DISPLACEMENT",order); CHKERRQ(ierr);
  ierr = mField.set_field_order(0,MBVERTEX,"DISPLACEMENT",1); CHKERRQ(ierr);
  //
  ierr = mField.set_field_order(0,MBTET,"MESH_NODE_POSITIONS",2); CHKERRQ(ierr);
  ierr = mField.set_field_order(0,MBTRI,"MESH_NODE_POSITIONS",2); CHKERRQ(ierr);
  ierr = mField.set_field_order(0,MBEDGE,"MESH_NODE_POSITIONS",2); CHKERRQ(ierr);
  ierr = mField.set_field_order(0,MBVERTEX,"MESH_NODE_POSITIONS",1); CHKERRQ(ierr);

  ierr = MetaNeummanForces::addNeumannBCElements(mField,"ELASTIC_PROB","DISPLACEMENT"); CHKERRQ(ierr);
  ierr = MetaNodalForces::addNodalForceElement(mField,"ELASTIC_PROB","DISPLACEMENT"); CHKERRQ(ierr);

  ierr = mField.add_finite_element("BODY_FORCE"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_row("BODY_FORCE","DISPLACEMENT"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_col("BODY_FORCE","DISPLACEMENT"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_data("BODY_FORCE","DISPLACEMENT"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_data("BODY_FORCE","MESH_NODE_POSITIONS"); CHKERRQ(ierr);
  ierr = mField.modify_problem_add_finite_element("ELASTIC_PROB","BODY_FORCE"); CHKERRQ(ierr);
  for(_IT_CUBITMESHSETS_BY_BCDATA_TYPE_FOR_LOOP_(mField,BLOCKSET|BLOCK_BODYFORCESSET,it)) {
    Range tets;
    rval = mField.get_moab().get_entities_by_type(it->meshset,MBTET,tets,true); CHKERR_PETSC(rval);
    ierr = mField.add_ents_to_finite_element_by_TETs(tets,"BODY_FORCE"); CHKERRQ(ierr);
  }

  //define fluid pressure finite elements
  FluidPressure fluid_pressure_fe(mField);
  fluid_pressure_fe.addNeumannFluidPressureBCElements("ELASTIC_PROB","DISPLACEMENT");

  //define elements for thermo elasticity if themperature field avelible
  ThermalStressElement thermal_stress_elem(mField);
  if(mField.check_field("TEMP")) {
    ierr = thermal_stress_elem.addThermalSterssElement("ELASTIC_PROB","ELASTIC","DISPLACEMENT","TEMP"); CHKERRQ(ierr);
  }

  //build database

  //build field
  ierr = mField.build_fields(); CHKERRQ(ierr);

  //build finite elemnts
  ierr = mField.build_finite_elements(); CHKERRQ(ierr);

  //build adjacencies
  ierr = mField.build_adjacencies(bit_level0); CHKERRQ(ierr);

  //build problem
  ierr = mField.build_problems(); CHKERRQ(ierr);

  //mesh partitioning 

  //partition
  ierr = mField.partition_problem("ELASTIC_PROB"); CHKERRQ(ierr);
  ierr = mField.partition_finite_elements("ELASTIC_PROB"); CHKERRQ(ierr);
  //what are ghost nodes, see Petsc Manual
  ierr = mField.partition_ghost_dofs("ELASTIC_PROB"); CHKERRQ(ierr);

  //mField.list_dofs_by_field_name("DISPLACEMENT",true);

  //print bcs
  ierr = mField.print_cubit_displacement_set(); CHKERRQ(ierr);
  ierr = mField.print_cubit_force_set(); CHKERRQ(ierr);
  //print block sets with materials
  ierr = mField.print_cubit_materials_set(); CHKERRQ(ierr);

  //create matrices
  Vec F,D;
  ierr = mField.VecCreateGhost("ELASTIC_PROB",ROW,&F); CHKERRQ(ierr);
  ierr = mField.VecCreateGhost("ELASTIC_PROB",COL,&D); CHKERRQ(ierr);

  Mat Aij;
  ierr = mField.MatCreateMPIAIJWithArrays("ELASTIC_PROB",&Aij); CHKERRQ(ierr);

  struct MyElasticFEMethod: public ElasticFEMethod {
    MyElasticFEMethod(FieldInterface& _mField,Mat &_Aij,Vec &_D,Vec& _F,double _lambda,double _mu): 
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

  Projection10NodeCoordsOnField ent_method_material(mField,"MESH_NODE_POSITIONS");
  ierr = mField.loop_dofs("MESH_NODE_POSITIONS",ent_method_material); CHKERRQ(ierr);

  //Assemble F and Aij
  DisplacementBCFEMethodPreAndPostProc my_dirichlet_bc(mField,"DISPLACEMENT",Aij,D,F);
  MyElasticFEMethod my_fe(mField,Aij,D,F,LAMBDA(young_modulus,poisson_ratio),MU(young_modulus,poisson_ratio));

  ierr = VecZeroEntries(F); CHKERRQ(ierr);
  ierr = VecGhostUpdateBegin(F,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = VecGhostUpdateEnd(F,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = MatZeroEntries(Aij); CHKERRQ(ierr);
  
  //preproc
  ierr = mField.problem_basic_method_preProcess("ELASTIC_PROB",my_dirichlet_bc); CHKERRQ(ierr);
  //loop elems
  ierr = mField.loop_finite_elements("ELASTIC_PROB","ELASTIC",my_fe);  CHKERRQ(ierr);
  //forces and preassures on surface
  boost::ptr_map<string,NeummanForcesSurface> neumann_forces;
  ierr = MetaNeummanForces::setNeumannFiniteElementOperators(mField,neumann_forces,F,"DISPLACEMENT"); CHKERRQ(ierr);
  boost::ptr_map<string,NeummanForcesSurface>::iterator mit = neumann_forces.begin();
  for(;mit!=neumann_forces.end();mit++) {
    ierr = mField.loop_finite_elements("ELASTIC_PROB",mit->first,mit->second->getLoopFe()); CHKERRQ(ierr);
  }
  //noadl forces
  boost::ptr_map<string,NodalForce> nodal_forces;
  ierr = MetaNodalForces::setNodalForceElementOperators(mField,nodal_forces,F,"DISPLACEMENT"); CHKERRQ(ierr);
  boost::ptr_map<string,NodalForce>::iterator fit = nodal_forces.begin();
  for(;fit!=nodal_forces.end();fit++) {
    ierr = mField.loop_finite_elements("ELASTIC_PROB",fit->first,fit->second->getLoopFe()); CHKERRQ(ierr);
  }
  //body forces
  BodyFroceConstantField body_forces_methods(mField);
  for(_IT_CUBITMESHSETS_BY_BCDATA_TYPE_FOR_LOOP_(mField,BLOCKSET|BLOCK_BODYFORCESSET,it)) {
    ierr = body_forces_methods.addBlock("DISPLACEMENT",F,it->get_msId()); CHKERRQ(ierr);
  }
  ierr = mField.loop_finite_elements("ELASTIC_PROB","BODY_FORCE",body_forces_methods.getLoopFe()); CHKERRQ(ierr);
  //fluid pressure
  ierr = fluid_pressure_fe.setNeumannFluidPressureFiniteElementOperators("DISPLACEMENT",F,false,true); CHKERRQ(ierr);
  ierr = mField.loop_finite_elements("ELASTIC_PROB","FLUID_PRESSURE_FE",fluid_pressure_fe.getLoopFe()); CHKERRQ(ierr);
  //postproc
  ierr = mField.problem_basic_method_postProcess("ELASTIC_PROB",my_dirichlet_bc); CHKERRQ(ierr);

  //set matrix possitives define and symetric for cholesky and icc preceonditionser
  ierr = MatSetOption(Aij,MAT_SPD,PETSC_TRUE); CHKERRQ(ierr);

  ierr = VecGhostUpdateBegin(F,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
  ierr = VecGhostUpdateEnd(F,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
  ierr = VecAssemblyBegin(F); CHKERRQ(ierr);
  ierr = VecAssemblyEnd(F); CHKERRQ(ierr);

  PetscSynchronizedFlush(PETSC_COMM_WORLD);

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

  if(mField.check_field("TEMP")) {
    
    //read time series and do thermo elastci analysis
    Vec F_thermal;
    ierr = VecDuplicate(F,&F_thermal); CHKERRQ(ierr);
    ierr = thermal_stress_elem.setThermalStressRhsOperators("DISPLACEMENT","TEMP",F_thermal); CHKERRQ(ierr);

    if( mField.check_series("THEMP_SERIES") ) {

      for(_IT_SERIES_STEPS_BY_NAME_FOR_LOOP_(mField,"THEMP_SERIES",sit)) {

	PetscPrintf(PETSC_COMM_WORLD,"Process step %d\n",sit->get_step_number());
	ierr = mField.load_series_data("THEMP_SERIES",sit->get_step_number()); CHKERRQ(ierr);
	ierr = VecZeroEntries(F_thermal); CHKERRQ(ierr);
	ierr = VecGhostUpdateBegin(F_thermal,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
	ierr = VecGhostUpdateEnd(F_thermal,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);

	ierr = mField.loop_finite_elements("ELASTIC_PROB","ELASTIC",thermal_stress_elem.getLoopThermalStressRhs()); CHKERRQ(ierr);
	ierr = VecGhostUpdateBegin(F_thermal,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
	ierr = VecGhostUpdateEnd(F_thermal,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
	ierr = VecAssemblyBegin(F_thermal); CHKERRQ(ierr);
	ierr = VecAssemblyEnd(F_thermal); CHKERRQ(ierr);

	PetscReal nrm_F;
	ierr = VecNorm(F,NORM_2,&nrm_F); CHKERRQ(ierr);
	PetscPrintf(PETSC_COMM_WORLD,"norm2 F = %6.4e\n",nrm_F);
  
	PetscReal nrm_F_thremal;
	ierr = VecNorm(F_thermal,NORM_2,&nrm_F_thremal); CHKERRQ(ierr);
	PetscPrintf(PETSC_COMM_WORLD,"norm2 F_thernal = %6.4e\n",nrm_F_thremal);

	ierr = VecScale(F_thermal,-1); CHKERRQ(ierr);
	ierr = VecAXPY(F_thermal,1,F); CHKERRQ(ierr);
	ierr = KSPSolve(solver,F_thermal,D); CHKERRQ(ierr);
	ierr = VecGhostUpdateBegin(D,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
	ierr = VecGhostUpdateEnd(D,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);

	//Save data on mesh
	ierr = mField.set_global_VecCreateGhost("ELASTIC_PROB",ROW,D,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);

	ostringstream o1;
	o1 << "out_" << sit->step_number << ".vtk";
	ostringstream o2;
	o2 << "out_post_proc_" << sit->step_number << ".vtk";
	ierr = write_soltion(mField,o1.str(),o2.str());   CHKERRQ(ierr);

      }

    } else {

      ierr = VecZeroEntries(F_thermal); CHKERRQ(ierr);
      ierr = VecGhostUpdateBegin(F_thermal,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
      ierr = VecGhostUpdateEnd(F_thermal,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);

      ierr = mField.loop_finite_elements("ELASTIC_PROB","ELASTIC",thermal_stress_elem.getLoopThermalStressRhs()); CHKERRQ(ierr);
      ierr = VecGhostUpdateBegin(F_thermal,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
      ierr = VecGhostUpdateEnd(F_thermal,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
      ierr = VecAssemblyBegin(F_thermal); CHKERRQ(ierr);
      ierr = VecAssemblyEnd(F_thermal); CHKERRQ(ierr);
 
      PetscReal nrm_F;
      ierr = VecNorm(F,NORM_2,&nrm_F); CHKERRQ(ierr);
      PetscPrintf(PETSC_COMM_WORLD,"norm2 F = %6.4e\n",nrm_F);
  
      PetscReal nrm_F_thremal;
      ierr = VecNorm(F_thermal,NORM_2,&nrm_F_thremal); CHKERRQ(ierr);
      PetscPrintf(PETSC_COMM_WORLD,"norm2 F_thernal = %6.4e\n",nrm_F_thremal);

      ierr = VecScale(F_thermal,-1); CHKERRQ(ierr);
      ierr = VecAXPY(F_thermal,1,F); CHKERRQ(ierr);
      ierr = KSPSolve(solver,F_thermal,D); CHKERRQ(ierr);
      ierr = VecGhostUpdateBegin(D,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
      ierr = VecGhostUpdateEnd(D,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);

      //Save data on mesh
      ierr = mField.set_global_VecCreateGhost("ELASTIC_PROB",ROW,D,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
      ierr = write_soltion(mField,"out.vtk","out_post_proc.vtk");   CHKERRQ(ierr);

    }

    ierr = VecDestroy(&F_thermal); CHKERRQ(ierr);

  } else {

    // elastic analys
  
    ierr = KSPSolve(solver,F,D); CHKERRQ(ierr);
    ierr = VecGhostUpdateBegin(D,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
    ierr = VecGhostUpdateEnd(D,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  
    //Save data on mesh
    ierr = mField.set_global_VecCreateGhost("ELASTIC_PROB",ROW,D,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
    ierr = write_soltion(mField,"out.vtk","out_post_proc.vtk");   CHKERRQ(ierr);

  }
      
  //Destroy matrices
  ierr = VecDestroy(&F); CHKERRQ(ierr);
  ierr = VecDestroy(&D); CHKERRQ(ierr);
  ierr = MatDestroy(&Aij); CHKERRQ(ierr);
  ierr = KSPDestroy(&solver); CHKERRQ(ierr);

  PetscSynchronizedFlush(PETSC_COMM_WORLD);

  PetscFinalize();

}

