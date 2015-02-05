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

#include <boost/numeric/ublas/vector_proxy.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/vector.hpp>

#include <SurfacePressure.hpp>
#include <NodalForce.hpp>
#include <FluidPressure.hpp>
#include <BodyForce.hpp>
#include <ThermalStressElement.hpp>

#include <FEMethod_LowLevelStudent.hpp>
#include <FEMethod_UpLevelStudent.hpp>

#include <PotsProcOnRefMesh.hpp>
#include <ElasticFEMethod.hpp>
#include <PostProcHookStresses.hpp>

using namespace boost::numeric;
using namespace ObosleteUsersModules;

ErrorCode rval;
PetscErrorCode ierr;

static char help[] = "...\n\n";

const double young_modulus = 1;
const double poisson_ratio = 0.0;

int main(int argc, char *argv[]) {

  PetscInitialize(&argc,&argv,(char *)0,help);

  moab::Core mb_instance;
  Interface& moab = mb_instance;
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

  PetscBool is_partitioned = PETSC_FALSE;
  ierr = PetscOptionsGetBool(PETSC_NULL,"-my_is_partitioned",&is_partitioned,&flg); CHKERRQ(ierr);
    
  if(is_partitioned == PETSC_TRUE) {
    //Read mesh to MOAB
    const char *option;
    option = "PARALLEL=BCAST_DELETE;PARALLEL_RESOLVE_SHARED_ENTS;PARTITION=PARALLEL_PARTITION;";
    rval = moab.load_file(mesh_file_name, 0, option); CHKERR_PETSC(rval); 
    rval = pcomm->resolve_shared_ents(0,3,0); CHKERR_PETSC(rval);
    rval = pcomm->resolve_shared_ents(0,3,1); CHKERR_PETSC(rval);
    rval = pcomm->resolve_shared_ents(0,3,2); CHKERR_PETSC(rval);
  } else {
    const char *option;
    option = "";
    rval = moab.load_file(mesh_file_name, 0, option); CHKERR_PETSC(rval); 
  }

  //Create MoFEM (Joseph) database
  MoFEM::Core core(moab);
  FieldInterface& m_field = core;

  // stl::bitset see for more details
  BitRefLevel bit_level0;
  bit_level0.set(0);
  ierr = m_field.seed_ref_level_3D(0,bit_level0); CHKERRQ(ierr);
  Range meshset_level0;
  ierr = m_field.get_entities_by_ref_level(bit_level0,BitRefLevel().set(),meshset_level0); CHKERRQ(ierr);
  PetscSynchronizedPrintf(PETSC_COMM_WORLD,"meshset_level0 %d\n",meshset_level0.size());
  PetscSynchronizedFlush(PETSC_COMM_WORLD,PETSC_STDOUT); 

  //Define problem

  //Fields
  ierr = m_field.add_field("DISPLACEMENT",H1,3,MF_ZERO); CHKERRQ(ierr);
  ierr = m_field.add_field("MESH_NODE_POSITIONS",H1,3,MF_ZERO); CHKERRQ(ierr);

  //FE
  ierr = m_field.add_finite_element("ELASTIC",MF_ZERO); CHKERRQ(ierr);

  //Define rows/cols and element data
  ierr = m_field.modify_finite_element_add_field_row("ELASTIC","DISPLACEMENT"); CHKERRQ(ierr);
  ierr = m_field.modify_finite_element_add_field_col("ELASTIC","DISPLACEMENT"); CHKERRQ(ierr);
  ierr = m_field.modify_finite_element_add_field_data("ELASTIC","DISPLACEMENT"); CHKERRQ(ierr);
  ierr = m_field.modify_finite_element_add_field_data("ELASTIC","MESH_NODE_POSITIONS"); CHKERRQ(ierr);

  //define problems
  ierr = m_field.add_problem("ELASTIC_PROB"); CHKERRQ(ierr);

  //set finite elements for problem
  ierr = m_field.modify_problem_add_finite_element("ELASTIC_PROB","ELASTIC"); CHKERRQ(ierr);

  //set refinment level for problem
  ierr = m_field.modify_problem_ref_level_add_bit("ELASTIC_PROB",bit_level0); CHKERRQ(ierr);

  //Declare problem

  //add entitities (by tets) to the field
  ierr = m_field.add_ents_to_field_by_TETs(0,"DISPLACEMENT",2); CHKERRQ(ierr);
  ierr = m_field.add_ents_to_field_by_TETs(0,"MESH_NODE_POSITIONS",2); CHKERRQ(ierr);

  //add finite elements entities
  ierr = m_field.add_ents_to_finite_element_EntType_by_bit_ref(bit_level0,"ELASTIC",MBTET); CHKERRQ(ierr);

  //set app. order
  //see Hierarchic Finite Element Bases on Unstructured Tetrahedral Meshes (Mark Ainsworth & Joe Coyle)
  ierr = m_field.set_field_order(0,MBTET,"DISPLACEMENT",order); CHKERRQ(ierr);
  ierr = m_field.set_field_order(0,MBTRI,"DISPLACEMENT",order); CHKERRQ(ierr);
  ierr = m_field.set_field_order(0,MBEDGE,"DISPLACEMENT",order); CHKERRQ(ierr);
  ierr = m_field.set_field_order(0,MBVERTEX,"DISPLACEMENT",1); CHKERRQ(ierr);
  ierr = m_field.set_field_order(0,MBTET,"MESH_NODE_POSITIONS",2); CHKERRQ(ierr);
  ierr = m_field.set_field_order(0,MBTRI,"MESH_NODE_POSITIONS",2); CHKERRQ(ierr);
  ierr = m_field.set_field_order(0,MBEDGE,"MESH_NODE_POSITIONS",2); CHKERRQ(ierr);
  ierr = m_field.set_field_order(0,MBVERTEX,"MESH_NODE_POSITIONS",1); CHKERRQ(ierr);

  ierr = MetaNeummanForces::addNeumannBCElements(m_field,"ELASTIC_PROB","DISPLACEMENT"); CHKERRQ(ierr);
  ierr = MetaNodalForces::addNodalForceElement(m_field,"ELASTIC_PROB","DISPLACEMENT"); CHKERRQ(ierr);

  ierr = m_field.add_finite_element("BODY_FORCE"); CHKERRQ(ierr);
  ierr = m_field.modify_finite_element_add_field_row("BODY_FORCE","DISPLACEMENT"); CHKERRQ(ierr);
  ierr = m_field.modify_finite_element_add_field_col("BODY_FORCE","DISPLACEMENT"); CHKERRQ(ierr);
  ierr = m_field.modify_finite_element_add_field_data("BODY_FORCE","DISPLACEMENT"); CHKERRQ(ierr);
  ierr = m_field.modify_finite_element_add_field_data("BODY_FORCE","MESH_NODE_POSITIONS"); CHKERRQ(ierr);
  ierr = m_field.modify_problem_add_finite_element("ELASTIC_PROB","BODY_FORCE"); CHKERRQ(ierr);
  for(_IT_CUBITMESHSETS_BY_BCDATA_TYPE_FOR_LOOP_(m_field,BLOCKSET|BODYFORCESSET,it)) {
    Range tets;
    rval = m_field.get_moab().get_entities_by_type(it->meshset,MBTET,tets,true); CHKERR_PETSC(rval);
    ierr = m_field.add_ents_to_finite_element_by_TETs(tets,"BODY_FORCE"); CHKERRQ(ierr);
  }

  //define fluid pressure finite elements
  FluidPressure fluid_pressure_fe(m_field);
  fluid_pressure_fe.addNeumannFluidPressureBCElements("ELASTIC_PROB","DISPLACEMENT");

  //define elements for thermo elasticity if themperature field avelible
  ThermalStressElement thermal_stress_elem(m_field);
  if(m_field.check_field("TEMP")) {
    ierr = thermal_stress_elem.addThermalSterssElement("ELASTIC_PROB","ELASTIC","DISPLACEMENT","TEMP"); CHKERRQ(ierr);
  }

  //build database

  //build field
  ierr = m_field.build_fields(); CHKERRQ(ierr);

  //build finite elemnts
  ierr = m_field.build_finite_elements(); CHKERRQ(ierr);

  //build adjacencies
  ierr = m_field.build_adjacencies(bit_level0); CHKERRQ(ierr);

  //build problem
  //ierr = m_field.build_problems(); CHKERRQ(ierr);
  if(is_partitioned) {
    ierr = m_field.build_partitioned_problems(1); CHKERRQ(ierr);
    ierr = m_field.partition_finite_elements("ELASTIC_PROB",true,0,pcomm->size(),1); CHKERRQ(ierr);
  } else {
    ierr = m_field.build_problems(); CHKERRQ(ierr);
    ierr = m_field.partition_problem("ELASTIC_PROB"); CHKERRQ(ierr);
    ierr = m_field.partition_finite_elements("ELASTIC_PROB"); CHKERRQ(ierr);
  }

  ierr = m_field.partition_ghost_dofs("ELASTIC_PROB"); CHKERRQ(ierr);

  //m_field.list_dofs_by_field_name("DISPLACEMENT",true);
  
  //print bcs
  ierr = m_field.print_cubit_displacement_set(); CHKERRQ(ierr);
  ierr = m_field.print_cubit_force_set(); CHKERRQ(ierr);
  //print block sets with materials
  ierr = m_field.print_cubit_materials_set(); CHKERRQ(ierr);

  //create matrices
  Vec F,D;
  ierr = m_field.VecCreateGhost("ELASTIC_PROB",ROW,&F); CHKERRQ(ierr);
  ierr = m_field.VecCreateGhost("ELASTIC_PROB",COL,&D); CHKERRQ(ierr);

  Mat Aij;
  ierr = m_field.MatCreateMPIAIJWithArrays("ELASTIC_PROB",&Aij); CHKERRQ(ierr);

  //Matrix View
  //MatView(Aij,PETSC_VIEWER_DRAW_WORLD);//PETSC_VIEWER_STDOUT_WORLD);
  //std::string wait;
  //std::cin >> wait;

  struct MyElasticFEMethod: public ElasticFEMethod {
    MyElasticFEMethod(FieldInterface& _m_field,Mat _Aij,Vec _D,Vec& _F,double _lambda,double _mu): 
      ElasticFEMethod(_m_field,_Aij,_D,_F,_lambda,_mu) {};

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

  Projection10NodeCoordsOnField ent_method_material(m_field,"MESH_NODE_POSITIONS");
  ierr = m_field.loop_dofs("MESH_NODE_POSITIONS",ent_method_material); CHKERRQ(ierr);

  //Assemble F and Aij
  DisplacementBCFEMethodPreAndPostProc my_dirichlet_bc(m_field,"DISPLACEMENT",Aij,D,F);
  MyElasticFEMethod my_fe(m_field,Aij,D,F,LAMBDA(young_modulus,poisson_ratio),MU(young_modulus,poisson_ratio));

  ierr = VecZeroEntries(F); CHKERRQ(ierr);
  ierr = VecGhostUpdateBegin(F,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = VecGhostUpdateEnd(F,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = MatZeroEntries(Aij); CHKERRQ(ierr);
  
  //preproc
  ierr = m_field.problem_basic_method_preProcess("ELASTIC_PROB",my_dirichlet_bc); CHKERRQ(ierr);
  //loop elems
  //PetscBarrier(PETSC_NULL);
  ierr = m_field.loop_finite_elements("ELASTIC_PROB","ELASTIC",my_fe);  CHKERRQ(ierr);
  //forces and preassures on surface
  boost::ptr_map<string,NeummanForcesSurface> neumann_forces;
  ierr = MetaNeummanForces::setNeumannFiniteElementOperators(m_field,neumann_forces,F,"DISPLACEMENT"); CHKERRQ(ierr);
  boost::ptr_map<string,NeummanForcesSurface>::iterator mit = neumann_forces.begin();
  for(;mit!=neumann_forces.end();mit++) {
    ierr = m_field.loop_finite_elements("ELASTIC_PROB",mit->first,mit->second->getLoopFe()); CHKERRQ(ierr);
  }
  //noadl forces
  boost::ptr_map<string,NodalForce> nodal_forces;
  ierr = MetaNodalForces::setNodalForceElementOperators(m_field,nodal_forces,F,"DISPLACEMENT"); CHKERRQ(ierr);
  boost::ptr_map<string,NodalForce>::iterator fit = nodal_forces.begin();
  for(;fit!=nodal_forces.end();fit++) {
    ierr = m_field.loop_finite_elements("ELASTIC_PROB",fit->first,fit->second->getLoopFe()); CHKERRQ(ierr);
  }
  //body forces
  BodyFroceConstantField body_forces_methods(m_field);
  for(_IT_CUBITMESHSETS_BY_BCDATA_TYPE_FOR_LOOP_(m_field,BLOCKSET|BODYFORCESSET,it)) {
    ierr = body_forces_methods.addBlock("DISPLACEMENT",F,it->get_msId()); CHKERRQ(ierr);
  }
  ierr = m_field.loop_finite_elements("ELASTIC_PROB","BODY_FORCE",body_forces_methods.getLoopFe()); CHKERRQ(ierr);
  //fluid pressure
  ierr = fluid_pressure_fe.setNeumannFluidPressureFiniteElementOperators("DISPLACEMENT",F,false,true); CHKERRQ(ierr);
  ierr = m_field.loop_finite_elements("ELASTIC_PROB","FLUID_PRESSURE_FE",fluid_pressure_fe.getLoopFe()); CHKERRQ(ierr);
  //postproc
  ierr = m_field.problem_basic_method_postProcess("ELASTIC_PROB",my_dirichlet_bc); CHKERRQ(ierr);

  //Matrix View
  /*MatView(Aij,PETSC_VIEWER_STDOUT_WORLD);
  MatView(Aij,PETSC_VIEWER_DRAW_WORLD);//PETSC_VIEWER_STDOUT_WORLD);
  std::string wait;
  std::cin >> wait;*/

  //set matrix possitives define and symetric for cholesky and icc preceonditionser
  ierr = MatSetOption(Aij,MAT_SPD,PETSC_TRUE); CHKERRQ(ierr);

  ierr = VecGhostUpdateBegin(F,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
  ierr = VecGhostUpdateEnd(F,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
  ierr = VecAssemblyBegin(F); CHKERRQ(ierr);
  ierr = VecAssemblyEnd(F); CHKERRQ(ierr);

  //PetscSynchronizedFlush(PETSC_COMM_WORLD,PETSC_STDOUT);

  //Solver
  KSP solver;
  ierr = KSPCreate(PETSC_COMM_WORLD,&solver); CHKERRQ(ierr);
  ierr = KSPSetOperators(solver,Aij,Aij); CHKERRQ(ierr);
  ierr = KSPSetFromOptions(solver); CHKERRQ(ierr);
  ierr = KSPSetUp(solver); CHKERRQ(ierr);

  {
    Tag th;
    int def_marker = 0;
    rval = moab.tag_get_handle("BLOCKID",1,MB_TYPE_INTEGER,th,MB_TAG_CREAT|MB_TAG_SPARSE,&def_marker); CHKERR_THROW(rval); 

    for(_IT_CUBITMESHSETS_BY_SET_TYPE_FOR_LOOP_(m_field,BLOCKSET,sit)) {

      int id = sit->get_msId();

      Range tets;
      rval = moab.get_entities_by_type(sit->meshset,MBTET,tets,true); CHKERR_PETSC(rval);
      Range::iterator it = tets.begin();
      for(;it!=tets.end();it++) {
	//cerr << id << endl;
	rval = moab.tag_set_data(th,&*it,1,&id); CHKERR_PETSC(rval);
      }

    }
  }

  PostPocOnRefinedMesh post_proc(m_field);
  ierr = post_proc.generateRefereneElemenMesh(); CHKERRQ(ierr);
  ierr = post_proc.addFieldValuesPostProc("DISPLACEMENT"); CHKERRQ(ierr);
  ierr = post_proc.addFieldValuesPostProc("MESH_NODE_POSITIONS"); CHKERRQ(ierr);
  ierr = post_proc.addFieldValuesGradientPostProc("DISPLACEMENT"); CHKERRQ(ierr);
  //add postpocessing for sresses
  post_proc.get_op_to_do_Rhs().push_back(
	  new PostPorcStress(
	    m_field,
	    post_proc.postProcMesh,
	    post_proc.mapGaussPts,
	    "DISPLACEMENT",
	    post_proc.commonData));

  if(m_field.check_field("TEMP")) {
    
    //read time series and do thermo elastci analysis
    Vec F_thermal;
    ierr = VecDuplicate(F,&F_thermal); CHKERRQ(ierr);
    ierr = thermal_stress_elem.setThermalStressRhsOperators("DISPLACEMENT","TEMP",F_thermal); CHKERRQ(ierr);

    SeriesRecorder *recorder_ptr;
    ierr = m_field.query_interface(recorder_ptr); CHKERRQ(ierr);
    if( recorder_ptr->check_series("THEMP_SERIES") ) {

      for(_IT_SERIES_STEPS_BY_NAME_FOR_LOOP_(recorder_ptr,"THEMP_SERIES",sit)) {

	PetscPrintf(PETSC_COMM_WORLD,"Process step %d\n",sit->get_step_number());
	ierr = recorder_ptr->load_series_data("THEMP_SERIES",sit->get_step_number()); CHKERRQ(ierr);
	ierr = VecZeroEntries(F_thermal); CHKERRQ(ierr);
	ierr = VecGhostUpdateBegin(F_thermal,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
	ierr = VecGhostUpdateEnd(F_thermal,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);

	ierr = m_field.loop_finite_elements("ELASTIC_PROB","ELASTIC",thermal_stress_elem.getLoopThermalStressRhs()); CHKERRQ(ierr);
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

        my_dirichlet_bc.snes_x = D;
        my_dirichlet_bc.snes_f = F_thermal;
	ierr = m_field.problem_basic_method_postProcess("ELASTIC_PROB",my_dirichlet_bc); CHKERRQ(ierr);

	ierr = KSPSolve(solver,F_thermal,D); CHKERRQ(ierr);
	ierr = VecGhostUpdateBegin(D,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
	ierr = VecGhostUpdateEnd(D,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);

	//Save data on mesh
	ierr = m_field.set_local_VecCreateGhost("ELASTIC_PROB",ROW,D,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
	ierr = m_field.loop_finite_elements("ELASTIC_PROB","ELASTIC",post_proc); CHKERRQ(ierr);
	ostringstream o1;
	o1 << "out_" << sit->step_number << ".h5m";
	rval = post_proc.postProcMesh.write_file(o1.str().c_str(),"MOAB","PARALLEL=WRITE_PART"); CHKERR_PETSC(rval);

      }

    } else {

      ierr = VecZeroEntries(F_thermal); CHKERRQ(ierr);
      ierr = VecGhostUpdateBegin(F_thermal,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
      ierr = VecGhostUpdateEnd(F_thermal,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);

      ierr = m_field.loop_finite_elements("ELASTIC_PROB","ELASTIC",thermal_stress_elem.getLoopThermalStressRhs()); CHKERRQ(ierr);
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

      my_dirichlet_bc.snes_x = D;
      my_dirichlet_bc.snes_f = F_thermal;
      ierr = m_field.problem_basic_method_postProcess("ELASTIC_PROB",my_dirichlet_bc); CHKERRQ(ierr);

      ierr = KSPSolve(solver,F_thermal,D); CHKERRQ(ierr);
      ierr = VecGhostUpdateBegin(D,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
      ierr = VecGhostUpdateEnd(D,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);

      //Save data on mesh
      ierr = m_field.set_local_VecCreateGhost("ELASTIC_PROB",ROW,D,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
      ierr = m_field.loop_finite_elements("ELASTIC_PROB","ELASTIC",post_proc); CHKERRQ(ierr);
      rval = post_proc.postProcMesh.write_file("out.h5m","MOAB","PARALLEL=WRITE_PART"); CHKERR_PETSC(rval);

    }

    ierr = VecDestroy(&F_thermal); CHKERRQ(ierr);

  } else {

    // elastic analys
  
    ierr = KSPSolve(solver,F,D); CHKERRQ(ierr);
    ierr = VecGhostUpdateBegin(D,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
    ierr = VecGhostUpdateEnd(D,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  
    //Save data on mesh
    ierr = m_field.set_local_VecCreateGhost("ELASTIC_PROB",COL,D,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
    ierr = m_field.loop_finite_elements("ELASTIC_PROB","ELASTIC",post_proc); CHKERRQ(ierr);
    rval = post_proc.postProcMesh.write_file("out.h5m","MOAB","PARALLEL=WRITE_PART"); CHKERR_PETSC(rval);

  }
      
  //Destroy matrices
  ierr = VecDestroy(&F); CHKERRQ(ierr);
  ierr = VecDestroy(&D); CHKERRQ(ierr);
  ierr = MatDestroy(&Aij); CHKERRQ(ierr);
  ierr = KSPDestroy(&solver); CHKERRQ(ierr);

  PetscFinalize();

}

