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

#ifdef WITH_TETGEN

#include <tetgen.h>
#ifdef REAL
  #undef REAL
#endif

#endif

#include <MoFEM.hpp>
using namespace MoFEM;

#include <DirichletBC.hpp>

#include <Projection10NodeCoordsOnField.hpp>

#include <boost/numeric/ublas/vector_proxy.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/vector.hpp>

#include <MethodForForceScaling.hpp>
#include <SurfacePressure.hpp>
#include <NodalForce.hpp>
#include <EdgeForce.hpp>

#include <FEMethod_LowLevelStudent.hpp>
#include <FEMethod_UpLevelStudent.hpp>
#include <FEMethod_SurfaceConstrains.hpp>

#include <PostProcVertexMethod.hpp>
#include <PostProcDisplacementAndStrainOnRefindedMesh.hpp>

extern "C" {
  #include <complex_for_lazy.h>
}

#include <ArcLengthTools.hpp>
#include <ConstrainMatrixCtx.hpp>
#include <FEMethod_ComplexForLazy.hpp>
#include <FEMethod_DriverComplexForLazy.hpp>

#include <SurfacePressureComplexForLazy.hpp>
#include <PostProcNonLinearElasticityStresseOnRefindedMesh.hpp>

#include <moab/Skinner.hpp>
#include <moab/AdaptiveKDTree.hpp>

#include <petsctime.h>

#include <ComplexConstArea.hpp>
#include <FaceSplittingTool.hpp>
#include <ConfigurationalFractureMechanics.hpp>
#include <MainCrackFunction.hpp>

#include <adolc/adolc.h>
#include <NonLinearElasticElement.hpp>
#include <Hooke.hpp>

#include <PostProcOnRefMesh.hpp>
#include <PostProcStresses.hpp>

using namespace ObosleteUsersModules;

PetscErrorCode main_spatial_solution(FieldInterface& m_field,ConfigurationalFractureMechanics& conf_prob) {
  PetscFunctionBegin;

  ErrorCode rval;
  PetscErrorCode ierr;

  ierr = conf_prob.set_material_fire_wall(m_field); CHKERRQ(ierr);

  Tag th_my_ref_level;
  rval = m_field.get_moab().tag_get_handle("_MY_REFINMENT_LEVEL",th_my_ref_level);
  const EntityHandle root_meshset = m_field.get_moab().get_root_set();
  BitRefLevel *ptr_bit_level0;
  rval = m_field.get_moab().tag_get_by_ptr(th_my_ref_level,&root_meshset,1,(const void**)&ptr_bit_level0); CHKERR_PETSC(rval);
  BitRefLevel& bit_level0 = *ptr_bit_level0;

  EntityHandle meshset_level0;
  rval = m_field.get_moab().create_meshset(MESHSET_SET,meshset_level0); CHKERR_PETSC(rval);
  ierr = m_field.get_entities_by_ref_level(bit_level0,BitRefLevel().set(),meshset_level0); CHKERRQ(ierr);

  ierr = conf_prob.spatial_problem_definition(m_field); CHKERRQ(ierr);

  //add finite elements entities
  ierr = m_field.add_ents_to_finite_element_EntType_by_bit_ref(bit_level0,"ELASTIC",MBTET); CHKERRQ(ierr);
  //set refinment level for problem
  ierr = m_field.modify_problem_ref_level_set_bit("ELASTIC_MECHANICS",bit_level0); CHKERRQ(ierr);

  //build field
  ierr = m_field.build_fields(); CHKERRQ(ierr);
  //build finite elemnts
  ierr = m_field.build_finite_elements(); CHKERRQ(ierr);
  //build adjacencies
  ierr = m_field.build_adjacencies(bit_level0); CHKERRQ(ierr);
  //build problem
  ierr = m_field.build_problems(); CHKERRQ(ierr);

  //partition problems
  ierr = conf_prob.spatial_partition_problems(m_field); CHKERRQ(ierr);

  //solve problem
  ierr = conf_prob.set_spatial_positions(m_field); CHKERRQ(ierr);
  ierr = conf_prob.set_material_positions(m_field); CHKERRQ(ierr);

  SNES snes;
  ierr = SNESCreate(PETSC_COMM_WORLD,&snes); CHKERRQ(ierr);
  ierr = SNESSetFromOptions(snes); CHKERRQ(ierr);
  ierr = conf_prob.solve_spatial_problem(m_field,&snes); CHKERRQ(ierr);
  ierr = SNESDestroy(&snes); CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

PetscErrorCode main_material_forces(FieldInterface& m_field,ConfigurationalFractureMechanics& conf_prob) {
  PetscFunctionBegin;

  ErrorCode rval;
  PetscErrorCode ierr;

  ierr = conf_prob.set_material_fire_wall(m_field); CHKERRQ(ierr);

  //ref meshset ref level 0
  Tag th_my_ref_level;
  rval = m_field.get_moab().tag_get_handle("_MY_REFINMENT_LEVEL",th_my_ref_level); CHKERR_PETSC(rval);
  const EntityHandle root_meshset = m_field.get_moab().get_root_set();
  BitRefLevel *ptr_bit_level0;
  rval = m_field.get_moab().tag_get_by_ptr(th_my_ref_level,&root_meshset,1,(const void**)&ptr_bit_level0); CHKERR_PETSC(rval);
  BitRefLevel& bit_level0 = *ptr_bit_level0;

  ierr = conf_prob.spatial_problem_definition(m_field); CHKERRQ(ierr);
  ierr = conf_prob.constrains_problem_definition(m_field); CHKERRQ(ierr);
  ierr = conf_prob.material_problem_definition(m_field); CHKERRQ(ierr);
  ierr = conf_prob.constrains_crack_front_problem_definition(m_field,"MATERIAL_MECHANICS"); CHKERRQ(ierr);

  //add finite elements entities
  ierr = m_field.add_ents_to_finite_element_EntType_by_bit_ref(bit_level0,"MATERIAL",MBTET); CHKERRQ(ierr);

  //set refinment level for problem
  ierr = m_field.modify_problem_ref_level_set_bit("MATERIAL_MECHANICS",bit_level0); CHKERRQ(ierr);

  //set refinment level for problem
  ierr = m_field.modify_problem_ref_level_set_bit("CCT_ALL_MATRIX",bit_level0); CHKERRQ(ierr);
  ierr = m_field.modify_problem_ref_level_set_bit("C_ALL_MATRIX",bit_level0); CHKERRQ(ierr);
  ierr = m_field.modify_problem_ref_level_set_bit("C_CRACKFRONT_MATRIX",bit_level0); CHKERRQ(ierr);
  ierr = m_field.modify_problem_ref_level_set_bit("CTC_CRACKFRONT_MATRIX",bit_level0); CHKERRQ(ierr);

  //build field
  ierr = m_field.build_fields(); CHKERRQ(ierr);
  //build finite elemnts
  ierr = m_field.build_finite_elements(); CHKERRQ(ierr);
  //build adjacencies
  ierr = m_field.build_adjacencies(bit_level0); CHKERRQ(ierr);
  //build problem
  ierr = m_field.build_problems(); CHKERRQ(ierr);

  //partition problems
  ierr = conf_prob.spatial_partition_problems(m_field); CHKERRQ(ierr);
  ierr = conf_prob.material_partition_problems(m_field); CHKERRQ(ierr);
  ierr = conf_prob.constrains_partition_problems(m_field,"MATERIAL_MECHANICS"); CHKERRQ(ierr);
  ierr = conf_prob.crackfront_partition_problems(m_field,"MATERIAL_MECHANICS"); CHKERRQ(ierr);

  //caculate material forces
  ierr = conf_prob.set_material_positions(m_field); CHKERRQ(ierr);
  ierr = conf_prob.front_projection_data(m_field,"MATERIAL_MECHANICS"); CHKERRQ(ierr);
  ierr = conf_prob.surface_projection_data(m_field,"MATERIAL_MECHANICS"); CHKERRQ(ierr);
  ierr = conf_prob.calculate_material_forces(m_field,"MATERIAL_MECHANICS","MATERIAL"); CHKERRQ(ierr);
  ierr = conf_prob.calculate_griffith_foces(m_field,"MATERIAL_MECHANICS"); CHKERRQ(ierr);
  ierr = conf_prob.project_force_vector(m_field,"MATERIAL_MECHANICS"); CHKERRQ(ierr);
  ierr = conf_prob.calculate_griffith_g(m_field,"MATERIAL_MECHANICS"); CHKERRQ(ierr);
  ierr = conf_prob.delete_surface_projection_data(m_field); CHKERRQ(ierr);
  ierr = conf_prob.delete_front_projection_data(m_field); CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

PetscErrorCode main_arc_length_setup(FieldInterface& m_field,ConfigurationalFractureMechanics& conf_prob) {
  PetscFunctionBegin;

  ErrorCode rval;
  PetscErrorCode ierr;
  //PetscBool flg = PETSC_TRUE;

  ierr = conf_prob.set_material_fire_wall(m_field); CHKERRQ(ierr);

  //ref meshset ref level 0
  Tag th_my_ref_level;
  rval = m_field.get_moab().tag_get_handle("_MY_REFINMENT_LEVEL",th_my_ref_level); CHKERR_PETSC(rval);
  const EntityHandle root_meshset = m_field.get_moab().get_root_set();
  BitRefLevel *ptr_bit_level0;
  rval = m_field.get_moab().tag_get_by_ptr(th_my_ref_level,&root_meshset,1,(const void**)&ptr_bit_level0); CHKERR_PETSC(rval);
  BitRefLevel& bit_level0 = *ptr_bit_level0;

  ierr = conf_prob.spatial_problem_definition(m_field); CHKERRQ(ierr);
  ierr = conf_prob.constrains_problem_definition(m_field); CHKERRQ(ierr);
  ierr = conf_prob.material_problem_definition(m_field); CHKERRQ(ierr);
  ierr = conf_prob.coupled_problem_definition(m_field); CHKERRQ(ierr);
  ierr = conf_prob.constrains_crack_front_problem_definition(m_field,"COUPLED_PROBLEM"); CHKERRQ(ierr);
  ierr = conf_prob.arclength_problem_definition(m_field); CHKERRQ(ierr);

  //add finite elements entities
  ierr = m_field.add_ents_to_finite_element_EntType_by_bit_ref(bit_level0,"ELASTIC_COUPLED",MBTET); CHKERRQ(ierr);
  ierr = m_field.add_ents_to_finite_element_EntType_by_bit_ref(bit_level0,"MESH_SMOOTHER",MBTET); CHKERRQ(ierr);
  ierr = m_field.add_ents_to_finite_element_EntType_by_bit_ref(bit_level0,"ELASTIC",MBTET); CHKERRQ(ierr);

  Range crack_front_edges;
  ierr = m_field.get_cubit_msId_entities_by_dimension(201,SIDESET,1,crack_front_edges,true); CHKERRQ(ierr);
  Range crack_front_nodes;
  rval = m_field.get_moab().get_connectivity(crack_front_edges,crack_front_nodes,true); CHKERR_PETSC(rval);
  Range crack_front_tets;
  rval = m_field.get_moab().get_adjacencies(crack_front_nodes,3,false,crack_front_tets,Interface::UNION); CHKERR_PETSC(rval);

  ierr = m_field.seed_finite_elements(crack_front_tets); CHKERRQ(ierr);
  ierr = m_field.add_ents_to_finite_element_by_TETs(crack_front_tets,"MATERIAL"); CHKERRQ(ierr);
  ierr = m_field.add_ents_to_finite_element_by_TETs(crack_front_tets,"MATERIAL_COUPLED"); CHKERRQ(ierr);

  //set refinment level for problem
  ierr = m_field.modify_problem_ref_level_set_bit("ELASTIC_MECHANICS",bit_level0); CHKERRQ(ierr);
  ierr = m_field.modify_problem_ref_level_set_bit("MATERIAL_MECHANICS",bit_level0); CHKERRQ(ierr);
  ierr = m_field.modify_problem_ref_level_set_bit("COUPLED_PROBLEM",bit_level0); CHKERRQ(ierr);
  ierr = m_field.modify_problem_ref_level_set_bit("CCT_ALL_MATRIX",bit_level0); CHKERRQ(ierr);
  ierr = m_field.modify_problem_ref_level_set_bit("C_ALL_MATRIX",bit_level0); CHKERRQ(ierr);
  ierr = m_field.modify_problem_ref_level_set_bit("C_CRACKFRONT_MATRIX",bit_level0); CHKERRQ(ierr);
  ierr = m_field.modify_problem_ref_level_set_bit("CTC_CRACKFRONT_MATRIX",bit_level0); CHKERRQ(ierr);

  //build field
  ierr = m_field.build_fields(); CHKERRQ(ierr);
  //build finite elemnts
  ierr = m_field.build_finite_elements(); CHKERRQ(ierr);
  //build adjacencies
  ierr = m_field.build_adjacencies(bit_level0); CHKERRQ(ierr);
  //build problem
  ierr = m_field.build_problems(); CHKERRQ(ierr);

  //partition problems
  ierr = conf_prob.spatial_partition_problems(m_field); CHKERRQ(ierr);
  ierr = conf_prob.material_partition_problems(m_field); CHKERRQ(ierr);
  ierr = conf_prob.coupled_partition_problems(m_field); CHKERRQ(ierr);
  ierr = conf_prob.constrains_partition_problems(m_field,"COUPLED_PROBLEM"); CHKERRQ(ierr);
  ierr = conf_prob.crackfront_partition_problems(m_field,"COUPLED_PROBLEM"); CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

PetscErrorCode main_arc_length_restart(FieldInterface& m_field,ConfigurationalFractureMechanics& conf_prob) {
  PetscFunctionBegin;

  ErrorCode rval;
  PetscErrorCode ierr;

  Tag th_griffith_force;
  rval = m_field.get_moab().tag_get_handle("GRIFFITH_FORCE",th_griffith_force); CHKERR_PETSC(rval);
  rval = m_field.get_moab().tag_delete(th_griffith_force); CHKERR_PETSC(rval);
  Tag th_freez;
  rval = m_field.get_moab().tag_get_handle("FROZEN_NODE",th_freez); CHKERR_PETSC(rval);
  rval = m_field.get_moab().tag_delete(th_freez); CHKERR_PETSC(rval);

  conf_prob.material_FirelWall->operator[](ConfigurationalFractureMechanics::FW_spatial_problem_definition) = 0;
  conf_prob.material_FirelWall->operator[](ConfigurationalFractureMechanics::FW_material_problem_definition) = 0;
  conf_prob.material_FirelWall->operator[](ConfigurationalFractureMechanics::FW_coupled_problem_definition) = 0;
  conf_prob.material_FirelWall->operator[](ConfigurationalFractureMechanics::FW_constrains_problem_definition) = 0;
  conf_prob.material_FirelWall->operator[](ConfigurationalFractureMechanics::FW_constrains_crack_front_problem_definition) = 0;
  conf_prob.material_FirelWall->operator[](ConfigurationalFractureMechanics::FW_set_spatial_positions) = 0;
  conf_prob.material_FirelWall->operator[](ConfigurationalFractureMechanics::FW_set_material_positions) = 0;

  ierr = main_arc_length_setup(m_field,conf_prob); CHKERRQ(ierr);
  ierr = m_field.check_number_of_ents_in_ents_field(); CHKERRQ(ierr);
  ierr = m_field.check_number_of_ents_in_ents_finite_element(); CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

PetscErrorCode main_rescale_load_factor(FieldInterface& m_field,ConfigurationalFractureMechanics& conf_prob) {
  PetscFunctionBegin;

  ErrorCode rval;
  PetscErrorCode ierr;

  double gc;
  PetscBool flg;
  ierr = PetscOptionsGetReal(PETSC_NULL,"-my_gc",&gc,&flg); CHKERRQ(ierr);
  if(flg != PETSC_TRUE) {
      SETERRQ(PETSC_COMM_SELF,1,"*** ERROR -my_gc (what is the fracture energy ?)");
  }

  const EntityHandle root_meshset = m_field.get_moab().get_root_set();

  Tag th_t_val;
  rval = m_field.get_moab().tag_get_handle("_LoadFactor_Scale_",th_t_val); CHKERR_PETSC(rval);
  double *load_factor_ptr;
  rval = m_field.get_moab().tag_get_by_ptr(th_t_val,&root_meshset,1,(const void**)&load_factor_ptr); CHKERR_THROW(rval);
  double& load_factor = *load_factor_ptr;

  double max_j = conf_prob.max_j;
  double a = fabs(max_j)/pow(load_factor,2);
  double new_load_factor = copysign(sqrt(gc/a),load_factor);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"\ncoefficient a = %6.4e\n",a); CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"new load factor value = %6.4e\n\n",new_load_factor); CHKERRQ(ierr);
  load_factor = new_load_factor;
  SNES snes;
  //solve spatial problem
  ierr = SNESCreate(PETSC_COMM_WORLD,&snes); CHKERRQ(ierr);
  ierr = SNESSetFromOptions(snes); CHKERRQ(ierr);
  ierr = conf_prob.solve_spatial_problem(m_field,&snes,false); CHKERRQ(ierr);
  ierr = SNESDestroy(&snes); CHKERRQ(ierr);
  ierr = conf_prob.calculate_material_forces(m_field,"COUPLED_PROBLEM","MATERIAL_COUPLED"); CHKERRQ(ierr);
  ierr = conf_prob.front_projection_data(m_field,"COUPLED_PROBLEM"); CHKERRQ(ierr);
  ierr = conf_prob.surface_projection_data(m_field,"COUPLED_PROBLEM"); CHKERRQ(ierr);
  ierr = conf_prob.project_force_vector(m_field,"COUPLED_PROBLEM"); CHKERRQ(ierr);
  ierr = conf_prob.calculate_griffith_foces(m_field,"COUPLED_PROBLEM"); CHKERRQ(ierr);
  ierr = conf_prob.calculate_griffith_g(m_field,"COUPLED_PROBLEM"); CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

PetscErrorCode main_arc_length_solve(FieldInterface& m_field,ConfigurationalFractureMechanics& conf_prob) {
  PetscFunctionBegin;

  ErrorCode rval;
  PetscErrorCode ierr;

  PetscBool flg = PETSC_TRUE;

  ParallelComm* pcomm = ParallelComm::get_pcomm(&m_field.get_moab(),MYPCOMM_INDEX);

  //ref meshset ref level 0
  Tag th_my_ref_level;
  rval = m_field.get_moab().tag_get_handle("_MY_REFINMENT_LEVEL",th_my_ref_level); CHKERR_PETSC(rval);
  const EntityHandle root_meshset = m_field.get_moab().get_root_set();
  BitRefLevel *ptr_bit_level0;
  rval = m_field.get_moab().tag_get_by_ptr(th_my_ref_level,&root_meshset,1,(const void**)&ptr_bit_level0); CHKERR_PETSC(rval);
  BitRefLevel& bit_level0 = *ptr_bit_level0;

  //caculate material forces
  ierr = conf_prob.set_material_positions(m_field); CHKERRQ(ierr);

  double da_0 = 0;
  ierr = PetscOptionsGetReal(PETSC_NULL,"-my_da",&da_0,&flg); CHKERRQ(ierr);
  if(flg != PETSC_TRUE) {
    SETERRQ(PETSC_COMM_SELF,1,"*** ERROR -my_da (what is the crack area increment ?)");
  }
  double da = da_0;

  int def_nb_load_steps = 0;
  Tag th_nb_load_steps;
  rval = m_field.get_moab().tag_get_handle("_NB_LOAD_STEPS",1,MB_TYPE_INTEGER,
    th_nb_load_steps,MB_TAG_CREAT|MB_TAG_SPARSE,&def_nb_load_steps);
  int *ptr_nb_load_steps;
  rval = m_field.get_moab().tag_get_by_ptr(th_nb_load_steps,&root_meshset,1,(const void**)&ptr_nb_load_steps); CHKERR_PETSC(rval);
  int &step = *ptr_nb_load_steps;
  if(step!=0) {
    step++;
  }

  PetscInt nb_load_steps = 0;
  ierr = PetscOptionsGetInt(PETSC_NULL,"-my_load_steps",&nb_load_steps,&flg); CHKERRQ(ierr);
  if(flg != PETSC_TRUE) {
    SETERRQ(PETSC_COMM_SELF,1,"*** ERROR -my_load_steps (what is the number of load_steps ?)");
  }

  if(step == 0) {
    //ierr = conf_prob.set_material_positions(m_field); CHKERRQ(ierr);
  }

  Tag th_t_val;
  rval = m_field.get_moab().tag_get_handle("_LoadFactor_Scale_",th_t_val); CHKERR_PETSC(rval);
  double *load_factor_ptr;
  rval = m_field.get_moab().tag_get_by_ptr(th_t_val,&root_meshset,1,(const void**)&load_factor_ptr); CHKERR_THROW(rval);
  double& load_factor = *load_factor_ptr;

  FaceSplittingTools face_splitting_tools(m_field);

  double fraction_treshold0 = 1e-1;
  ierr = PetscOptionsGetReal(PETSC_NULL,"-my_fraction_treshold",&fraction_treshold0,&flg); CHKERRQ(ierr);

  for(int aa = 0;step<nb_load_steps;step++,aa++) {

    PetscLogDouble t1,t2;
    PetscLogDouble v1,v2;
    ierr = PetscTime(&v1); CHKERRQ(ierr);
    ierr = PetscGetCPUTime(&t1); CHKERRQ(ierr);

    ierr = PetscPrintf(PETSC_COMM_WORLD,"\n\n** number of step = %D\n",step); CHKERRQ(ierr);

    if(aa == 0) {
      ierr = conf_prob.calculate_material_forces(m_field,"COUPLED_PROBLEM","MATERIAL_COUPLED"); CHKERRQ(ierr);
      ierr = conf_prob.front_projection_data(m_field,"COUPLED_PROBLEM"); CHKERRQ(ierr);
      ierr = conf_prob.surface_projection_data(m_field,"COUPLED_PROBLEM"); CHKERRQ(ierr);
      ierr = conf_prob.project_force_vector(m_field,"COUPLED_PROBLEM"); CHKERRQ(ierr);
      ierr = conf_prob.calculate_griffith_foces(m_field,"COUPLED_PROBLEM"); CHKERRQ(ierr);
      ierr = conf_prob.calculate_griffith_g(m_field,"COUPLED_PROBLEM"); CHKERRQ(ierr);
      ierr = conf_prob.delete_surface_projection_data(m_field); CHKERRQ(ierr);
      ierr = conf_prob.delete_front_projection_data(m_field); CHKERRQ(ierr);
    }

    double gc;
    PetscBool flg;
    ierr = PetscOptionsGetReal(PETSC_NULL,"-my_gc",&gc,&flg); CHKERRQ(ierr);
    if(flg != PETSC_TRUE) {
      SETERRQ(PETSC_COMM_SELF,1,"*** ERROR -my_gc (what is the fracture energy ?)");
    }

    //calculate initial load factor
    if(step == 0) {
      ierr = main_rescale_load_factor(m_field,conf_prob); CHKERRQ(ierr);
    }

    ierr = PetscPrintf(PETSC_COMM_WORLD,"* da = %6.4e\n",da); CHKERRQ(ierr);
    SNES snes;
    ierr = SNESCreate(PETSC_COMM_WORLD,&snes); CHKERRQ(ierr);
    ierr = SNESSetFromOptions(snes); CHKERRQ(ierr);

    Vec D0;
    ierr = m_field.VecCreateGhost("COUPLED_PROBLEM",COL,&D0); CHKERRQ(ierr);
    ierr = m_field.set_local_ghost_vector("COUPLED_PROBLEM",COL,D0,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
    double load_factor0 = load_factor;

    int nb_sub_steps;
    ierr = PetscOptionsGetInt("","-my_nb_sub_steps",&nb_sub_steps,&flg); CHKERRQ(ierr);
    if(flg != PETSC_TRUE) {
      nb_sub_steps = 20;
    }

    int its_d;
    ierr = PetscOptionsGetInt("","-my_its_d",&its_d,&flg); CHKERRQ(ierr);
    double penalty0 = 0;
    ierr = PetscOptionsGetReal(PETSC_NULL,"-my_front_penalty",&penalty0,&flg); CHKERRQ(ierr);

    double *load_factor_ptr;
    rval = m_field.get_moab().tag_get_by_ptr(th_t_val,&root_meshset,1,(const void**)&load_factor_ptr); CHKERR_THROW(rval);
    double& load_factor = *load_factor_ptr;

    double fraction_treshold = fraction_treshold0;
    bool at_least_one_step_converged = false;
    conf_prob.freeze_all_but_one = false;
    double penalty = (aa == 0) ? 0 : penalty0;
    double _da_ = (aa == 0) ? 0 : da;

    int ii = 0;
    for(;ii<nb_sub_steps;ii++) {
      ierr = PetscPrintf(PETSC_COMM_WORLD,"\n* number of substeps = %D _da_ = %6.4e\n\n",ii,_da_); CHKERRQ(ierr);

      SNESConvergedReason reason0;
      ierr = conf_prob.solve_coupled_problem(m_field,&snes,_da_,penalty,fraction_treshold,reason0); CHKERRQ(ierr);
      if(reason0 < 0) {
        fraction_treshold *= 0.25;
      }
      SNESConvergedReason reason;
      ierr = SNESGetConvergedReason(snes,&reason); CHKERRQ(ierr);
      if(conf_prob.total_its == 0) break;
      if(reason > 0) {
        if(da > 0) {
          if(aa > 0 && ii == 0) {
            if(flg != PETSC_TRUE) {
              its_d = 5;
            }
            double gamma = 0.5,reduction = 1;
            reduction = pow((double)its_d/(double)(conf_prob.total_its+1),gamma);
            const double max_da_reduction = 10;
            if(reduction<1 || da < max_da_reduction*da_0) {
              ierr = PetscPrintf(PETSC_COMM_WORLD,"\n* change of da = %6.4e\n\n",reduction); CHKERRQ(ierr);
              da = fmin(da*reduction,max_da_reduction*da_0);
            }
          }
        }
        at_least_one_step_converged = true;
        conf_prob.freeze_all_but_one = false;
        ierr = m_field.set_local_ghost_vector("COUPLED_PROBLEM",COL,D0,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
        load_factor0 = load_factor;
        ierr = conf_prob.calculate_material_forces(m_field,"COUPLED_PROBLEM","MATERIAL_COUPLED"); CHKERRQ(ierr);
        ierr = conf_prob.front_projection_data(m_field,"COUPLED_PROBLEM"); CHKERRQ(ierr);
        ierr = conf_prob.surface_projection_data(m_field,"COUPLED_PROBLEM"); CHKERRQ(ierr);
        ierr = conf_prob.project_force_vector(m_field,"COUPLED_PROBLEM"); CHKERRQ(ierr);
        ierr = conf_prob.calculate_griffith_foces(m_field,"COUPLED_PROBLEM"); CHKERRQ(ierr);
        ierr = conf_prob.calculate_griffith_g(m_field,"COUPLED_PROBLEM"); CHKERRQ(ierr);
        ierr = conf_prob.delete_surface_projection_data(m_field); CHKERRQ(ierr);
        ierr = conf_prob.delete_front_projection_data(m_field); CHKERRQ(ierr);
        _da_ = 0;
      } else {
        ierr = PetscPrintf(PETSC_COMM_WORLD,"* reset unknowns vector\n"); CHKERRQ(ierr);
        ierr = m_field.set_global_ghost_vector("COUPLED_PROBLEM",COL,D0,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
        load_factor = load_factor0;
        ierr = PetscPrintf(PETSC_COMM_WORLD,"* failed to converge, recalculate spatial positions\n"); CHKERRQ(ierr);
        {
          SNES snes_spatial;
          ierr = SNESCreate(PETSC_COMM_WORLD,&snes_spatial); CHKERRQ(ierr);
          ierr = SNESSetFromOptions(snes_spatial); CHKERRQ(ierr);
          ierr = conf_prob.solve_spatial_problem(m_field,&snes_spatial,false); CHKERRQ(ierr);
          ierr = SNESDestroy(&snes_spatial); CHKERRQ(ierr);
          //ierr = main_rescale_load_factor(m_field,conf_prob); CHKERRQ(ierr);
          ierr = conf_prob.calculate_material_forces(m_field,"COUPLED_PROBLEM","MATERIAL_COUPLED"); CHKERRQ(ierr);
          ierr = conf_prob.front_projection_data(m_field,"COUPLED_PROBLEM"); CHKERRQ(ierr);
          ierr = conf_prob.surface_projection_data(m_field,"COUPLED_PROBLEM"); CHKERRQ(ierr);
          ierr = conf_prob.project_force_vector(m_field,"COUPLED_PROBLEM"); CHKERRQ(ierr);
          ierr = conf_prob.calculate_griffith_foces(m_field,"COUPLED_PROBLEM"); CHKERRQ(ierr);
          ierr = conf_prob.calculate_griffith_g(m_field,"COUPLED_PROBLEM"); CHKERRQ(ierr);
        }
        if(!at_least_one_step_converged) {
          if(da > 0 && aa > 0) {
            da = 0.5*da;
            _da_ = da;
            ierr = PetscPrintf(PETSC_COMM_WORLD,"* failed to converge, set da = %6.4e ( 0.5 )\n",_da_); CHKERRQ(ierr);
          }
        }
        ///PetscAttachDebugger();

        if(_da_ == 0) {
          if(conf_prob.freeze_all_but_one) {
            SETERRQ(PETSC_COMM_SELF,1,"* unable to converge");
          } else {
            ierr = PetscPrintf(PETSC_COMM_WORLD,"* freez all but one\n"); CHKERRQ(ierr);
            conf_prob.freeze_all_but_one = true;
          }
        }
      }

    }
    ierr = VecDestroy(&D0); CHKERRQ(ierr);
    ierr = SNESDestroy(&snes); CHKERRQ(ierr);

    // Reduce step size, if many time need to freeze and unfreeze nodes
    double reduction = pow(2./(fmax(ii+1,3)-1),0.25);
    ierr = PetscPrintf(PETSC_COMM_WORLD,"\n* change of da = %6.4e (nb. sub-steps %d) \n\n",reduction,ii); CHKERRQ(ierr);
    da *= reduction;

    // Set max reduction to 0.1 of original value of da
    da = fmax(da,1e-1*da_0);

    ierr = conf_prob.set_coordinates_from_material_solution(m_field,false); CHKERRQ(ierr);

    ierr = PetscPrintf(PETSC_COMM_WORLD,
      "load_path: %4D Area %6.4e Lambda %6.4e Energy %6.4e\n",
      step,conf_prob.aRea,conf_prob.lambda,conf_prob.energy
    ); CHKERRQ(ierr);

    {

      Hooke<adouble> hooke_adouble;
      Hooke<double> hooke_double;

      PostPocOnRefinedMesh post_proc(m_field);
      ierr = post_proc.generateReferenceElementMesh(); CHKERRQ(ierr);
      ierr = post_proc.addFieldValuesPostProc("SPATIAL_POSITION"); CHKERRQ(ierr);
      ierr = post_proc.addFieldValuesPostProc("MESH_NODE_POSITIONS"); CHKERRQ(ierr);
      ierr = post_proc.addFieldValuesGradientPostProc("SPATIAL_POSITION"); CHKERRQ(ierr);

      map<int,NonlinearElasticElement::BlockData> set_of_blocks;
      for(_IT_CUBITMESHSETS_BY_BCDATA_TYPE_FOR_LOOP_(m_field,BLOCKSET|MAT_ELASTICSET,it)) {
        Mat_Elastic mydata;
        ierr = it->get_attribute_data_structure(mydata); CHKERRQ(ierr);
        int id = it->get_msId();
        EntityHandle meshset = it->get_meshset();
        rval = m_field.get_moab().get_entities_by_type(meshset,MBTET,set_of_blocks[id].tEts,true); CHKERR_PETSC(rval);
        set_of_blocks[id].iD = id;
        set_of_blocks[id].E = mydata.data.Young;
        set_of_blocks[id].PoissonRatio = mydata.data.Poisson;
        set_of_blocks[id].materialDoublePtr = &hooke_double;
        set_of_blocks[id].materialAdoublePtr = &hooke_adouble;
      }

      map<int,NonlinearElasticElement::BlockData>::iterator sit = set_of_blocks.begin();
      for (; sit != set_of_blocks.end(); sit++) {
        post_proc.getOpPtrVector().push_back(
          new PostPorcStress(
            post_proc.postProcMesh,
            post_proc.mapGaussPts,
            "SPATIAL_POSITION",
            sit->second,
            post_proc.commonData
          )
        );

      }

      ierr = m_field.loop_finite_elements("ELASTIC_MECHANICS","ELASTIC",post_proc); CHKERRQ(ierr);
      ostringstream ss;
      ss << "out_load_step_" << step << ".h5m";
      rval = post_proc.postProcMesh.write_file(ss.str().c_str(),"MOAB","PARALLEL=WRITE_PART"); CHKERR_PETSC(rval);

    }

    if(pcomm->rank()==0) {
      /*EntityHandle out_meshset;
      rval = m_field.get_moab().create_meshset(MESHSET_SET,out_meshset); CHKERR_PETSC(rval);
      ierr = m_field.get_problem_finite_elements_entities("COUPLED_PROBLEM","MATERIAL_COUPLED",out_meshset); CHKERRQ(ierr);
      ostringstream ss;
      ss << "out_load_step_" << step << ".vtk";
      rval = m_field.get_moab().write_file(ss.str().c_str(),"VTK","",&out_meshset,1); CHKERR_PETSC(rval);
      rval = m_field.get_moab().delete_entities(&out_meshset,1); CHKERR_PETSC(rval);*/

      {
        Range SurfacesFaces;
        ierr = m_field.get_cubit_msId_entities_by_dimension(102,SIDESET,2,SurfacesFaces,true); CHKERRQ(ierr);
        Range CrackSurfacesFaces;
        ierr = m_field.get_cubit_msId_entities_by_dimension(200,SIDESET,2,CrackSurfacesFaces,true); CHKERRQ(ierr);
        for(_IT_CUBITMESHSETS_BY_SET_TYPE_FOR_LOOP_(m_field,SIDESET,it)) {
          int msId = it->get_msId();
          if((msId < 10200)||(msId >= 10300)) continue;
          Range SurfacesFaces_msId;
          ierr = m_field.get_cubit_msId_entities_by_dimension(msId,SIDESET,2,SurfacesFaces_msId,true); CHKERRQ(ierr);
          SurfacesFaces.insert(SurfacesFaces_msId.begin(),SurfacesFaces_msId.end());
        }

        Range level_tris;
        ierr = m_field.get_entities_by_type_and_ref_level(bit_level0,BitRefLevel().set(),MBTRI,level_tris); CHKERRQ(ierr);
        SurfacesFaces = intersect(SurfacesFaces,level_tris);
        CrackSurfacesFaces = intersect(CrackSurfacesFaces,level_tris);

        EntityHandle out_meshset;
        rval = m_field.get_moab().create_meshset(MESHSET_SET,out_meshset); CHKERR_PETSC(rval);
        rval = m_field.get_moab().add_entities(out_meshset,CrackSurfacesFaces); CHKERR_PETSC(rval);
        ostringstream ss1;
        ss1 << "out_crack_surface_" << step << ".vtk";
        rval = m_field.get_moab().write_file(ss1.str().c_str(),"VTK","",&out_meshset,1); CHKERR_PETSC(rval);
        rval = m_field.get_moab().add_entities(out_meshset,SurfacesFaces); CHKERR_PETSC(rval);
        ostringstream ss2;
        ss2 << "out_surface_" << step << ".vtk";
        rval = m_field.get_moab().write_file(ss2.str().c_str(),"VTK","",&out_meshset,1); CHKERR_PETSC(rval);
        rval = m_field.get_moab().delete_entities(&out_meshset,1); CHKERR_PETSC(rval);
      }

      ostringstream ss3;
      ss3 << "restart_" << step << ".h5m";
      rval = m_field.get_moab().write_file(ss3.str().c_str()); CHKERR_PETSC(rval);

      const MoFEMProblem *problemPtr;
      ierr = m_field.get_problem("COUPLED_PROBLEM",&problemPtr); CHKERRQ(ierr);

      for(_IT_CUBITMESHSETS_BY_BCDATA_TYPE_FOR_LOOP_(m_field,BLOCKSET|UNKNOWNCUBITNAME,it)) {
        if(it->get_name() != "LoadPath") continue;

        Range nodes;
        rval = m_field.get_moab().get_entities_by_type(it->meshset,MBVERTEX,nodes,true); CHKERR_PETSC(rval);
        for(Range::iterator nit = nodes.begin();nit!=nodes.end();nit++) {
          double coords[3];
          rval = m_field.get_moab().get_coords(&*nit,1,coords); CHKERR_PETSC(rval);
          for(_IT_NUMEREDDOFMOFEMENTITY_ROW_BY_ENT_FOR_LOOP_(problemPtr,*nit,dof)) {
            if(dof->get_name()!="SPATIAL_POSITION") continue;
            ierr = PetscPrintf(PETSC_COMM_WORLD,
              "load_path_disp ent %ld dim %d "
              "coords ( %8.6f %8.6f %8.6f ) "
              "val %6.4e Lambda %6.4e\n",
              dof->get_ent(),dof->get_dof_rank(),
              coords[0],coords[1],coords[2],
              dof->get_FieldData()-coords[dof->get_dof_rank()],
              conf_prob.lambda); CHKERRQ(ierr);
            }
          }

        }
      }

    ierr = PetscTime(&v2);CHKERRQ(ierr);
    ierr = PetscGetCPUTime(&t2);CHKERRQ(ierr);
    PetscSynchronizedPrintf(PETSC_COMM_WORLD,"Step Time Rank %d Time = %f CPU Time = %f\n",pcomm->rank(),v2-v1,t2-t1);
    PetscSynchronizedFlush(PETSC_COMM_WORLD,PETSC_STDOUT);

    #ifdef WITH_TETGEN


    PetscBool flg_do_split = PETSC_TRUE;
    ierr = PetscOptionsGetBool(
      PETSC_NULL,"-my_do_split",&flg_do_split,PETSC_NULL
    ); CHKERRQ(ierr);
    PetscBool flg_do_tetgen = PETSC_TRUE;
    ierr = PetscOptionsGetBool(
      PETSC_NULL,"-my_do_tetgen",&flg_do_tetgen,PETSC_NULL
    ); CHKERRQ(ierr);

    bool do_tetgen = true;
    if(flg_do_split) {
      Range edges_to_cat;
      ierr = face_splitting_tools.getCornerEdges(edges_to_cat,10); CHKERRQ(ierr);
      if(edges_to_cat.size()>0) {
        Range new_nodes;
        ierr = face_splitting_tools.propagateBySplit(new_nodes,edges_to_cat,10); CHKERRQ(ierr);
        ierr = face_splitting_tools.cornerProblem(new_nodes,10); CHKERRQ(ierr);
        //if(!new_nodes.empty()) do_tetgen = false;
      }
    }

    if(do_tetgen) {
      if(flg_do_tetgen) {
        face_splitting_tools.moabTetGenMap.clear();
        face_splitting_tools.tetGenMoabMap.clear();
        face_splitting_tools.tetGenData.clear();
        vector<string> switches1;
        if(pcomm->rank() == 0) {
          switches1.push_back("rp178sqRS0JVV");
          ierr = face_splitting_tools.rebuildMeshWithTetGen(switches1,0); CHKERRQ(ierr);
        } else {
          switches1.push_back("rp178sqRS0JQ");
          ierr = face_splitting_tools.rebuildMeshWithTetGen(switches1,0); CHKERRQ(ierr);
        }
      }
    }

    bit_level0 = BitRefLevel().set(face_splitting_tools.meshIntefaceBitLevels.back());
    //retart analysis
    ierr = main_arc_length_restart(m_field,conf_prob); CHKERRQ(ierr);
    //project and set coords
    conf_prob.material_FirelWall->operator[](ConfigurationalFractureMechanics::FW_set_spatial_positions) = 0;
    conf_prob.material_FirelWall->operator[](ConfigurationalFractureMechanics::FW_set_material_positions) = 0;
    ierr = conf_prob.set_spatial_positions(m_field); CHKERRQ(ierr);
    ierr = conf_prob.set_material_positions(m_field); CHKERRQ(ierr);
    //solve spatial problem
    ierr = SNESCreate(PETSC_COMM_WORLD,&snes); CHKERRQ(ierr);
    ierr = SNESSetFromOptions(snes); CHKERRQ(ierr);
    ierr = conf_prob.solve_spatial_problem(m_field,&snes,false); CHKERRQ(ierr);
    ierr = SNESDestroy(&snes); CHKERRQ(ierr);
    //calculate Griffth forces
    ierr = conf_prob.calculate_material_forces(m_field,"COUPLED_PROBLEM","MATERIAL_COUPLED"); CHKERRQ(ierr);
    ierr = conf_prob.front_projection_data(m_field,"COUPLED_PROBLEM"); CHKERRQ(ierr);
    ierr = conf_prob.surface_projection_data(m_field,"COUPLED_PROBLEM"); CHKERRQ(ierr);
    ierr = conf_prob.project_force_vector(m_field,"COUPLED_PROBLEM"); CHKERRQ(ierr);
    ierr = conf_prob.calculate_griffith_foces(m_field,"COUPLED_PROBLEM"); CHKERRQ(ierr);
    ierr = conf_prob.calculate_griffith_g(m_field,"COUPLED_PROBLEM"); CHKERRQ(ierr);
    ierr = conf_prob.delete_surface_projection_data(m_field); CHKERRQ(ierr);
    ierr = conf_prob.delete_front_projection_data(m_field); CHKERRQ(ierr);
    #endif


  }

  PetscFunctionReturn(0);
}
