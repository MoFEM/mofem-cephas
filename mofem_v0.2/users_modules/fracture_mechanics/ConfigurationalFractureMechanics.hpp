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

#ifndef __CONFIGURATIONAL_MECHANICS_HPP__
#define __CONFIGURATIONAL_MECHANICS_HPP__

using namespace ObosleteUsersModules;

struct ConfigurationalFractureMechanics {

  Tag th_MaterialFireWall;
  typedef bitset<17> Material_FirelWall_def;
  Material_FirelWall_def *material_FirelWall;

  enum FirWall {
    FW_add_crack = 1,
    FW_refine_near_crack_tip,
    FW_set_load_factor,
    FW_spatial_problem_definition,
    FW_material_problem_definition,
    FW_coupled_problem_definition,
    FW_constrains_problem_definition,
    FW_constrains_crack_front_problem_definition,
    FW_set_spatial_positions,
    FW_set_material_positions,
    FW_arc_lenhghat_definition,
    FW_thermal_field
  };

  ConstrainMatrixCtx *projSurfaceCtx,*projFrontCtx;

  BitRefLevel *ptr_bit_level0;
  ConfigurationalFractureMechanics(FieldInterface& m_field): projSurfaceCtx(NULL),projFrontCtx(NULL) {

    ErrorCode rval;
    Tag th_my_ref_level;
    rval = m_field.get_moab().tag_get_handle("_MY_REFINMENT_LEVEL",th_my_ref_level); CHKERR_THROW(rval);
    const EntityHandle root_meshset = m_field.get_moab().get_root_set();
    rval = m_field.get_moab().tag_get_by_ptr(th_my_ref_level,&root_meshset,1,(const void**)&ptr_bit_level0); CHKERR_THROW(rval);

    fe_post_proc_stresses_method = NULL;

  };

  ~ConfigurationalFractureMechanics() {
    if(fe_post_proc_stresses_method!=NULL) {
      delete fe_post_proc_stresses_method;
    }
  }

  PetscErrorCode set_material_fire_wall(FieldInterface& m_field);
  PetscErrorCode thermal_field(FieldInterface& m_field);
  PetscErrorCode spatial_problem_definition(FieldInterface& m_field);
  PetscErrorCode material_problem_definition(FieldInterface& m_field);
  PetscErrorCode coupled_problem_definition(FieldInterface& m_field);
  PetscErrorCode arclength_problem_definition(FieldInterface& m_field);
  PetscErrorCode constrains_problem_definition(FieldInterface& m_field);
  PetscErrorCode constrains_crack_front_problem_definition(FieldInterface& m_field,string problem);
  PetscErrorCode spatial_partition_problems(FieldInterface& m_field);
  PetscErrorCode material_partition_problems(FieldInterface& m_field);
  PetscErrorCode coupled_partition_problems(FieldInterface& m_field);
  PetscErrorCode constrains_partition_problems(FieldInterface& m_field,string problem);
  PetscErrorCode crackfront_partition_problems(FieldInterface& m_field,string problem);
  PetscErrorCode set_spatial_positions(FieldInterface& m_field);
  PetscErrorCode set_material_positions(FieldInterface& m_field);
  PetscErrorCode set_coordinates_from_material_solution(FieldInterface& m_field,bool only_crack_front = false);

  PostProcStressNonLinearElasticity *fe_post_proc_stresses_method;
  PetscErrorCode solve_spatial_problem(FieldInterface& m_field,SNES *snes,bool postproc = true);
  PetscErrorCode solve_material_problem(FieldInterface& m_field,SNES *snes);

  PetscErrorCode fix_all_but_one(FieldInterface& m_field,Range &fix_nodes,const double fraction_treshold);

  double aRea,lambda,energy;
  int nb_un_freez_nodes;
  bool freeze_all_but_one;
  int total_its;
  PetscErrorCode solve_coupled_problem(FieldInterface& m_field,SNES *snes,double da,const double fraction_treshold = 5e-2);

  PetscErrorCode calculate_material_forces(FieldInterface& m_field,string problem,string fe);
  PetscErrorCode surface_projection_data(FieldInterface& m_field,string problem);
  PetscErrorCode delete_surface_projection_data(FieldInterface& m_field);
  PetscErrorCode project_force_vector(FieldInterface& m_field,string problem);
  PetscErrorCode front_projection_data(FieldInterface& m_field,string problem);
  PetscErrorCode delete_front_projection_data(FieldInterface& m_field);
  PetscErrorCode griffith_force_vector(FieldInterface& m_field,string problem);

  PetscErrorCode project_form_th_projection_tag(FieldInterface& m_field,string problem,bool do_not_project = false);

  map<EntityHandle,double> map_ent_g,map_ent_j,map_ent_work;
  PetscScalar ave_g,min_g,max_g;
  PetscScalar ave_j,min_j,max_j;
  PetscErrorCode griffith_g(FieldInterface& m_field,string problem);

 struct FrontAreaArcLengthControl: public FEMethod {

    FieldInterface& mField;
    ConfigurationalFractureMechanics *conf_prob;
    ArcLengthCtx* arc_ptr;

    Vec ghostDiag;
    Range crackSurfacesFaces;
    Range crackFrontNodes;
    PetscInt *isIdx;
    IS isSurface;
    Vec surfaceDofs;
    VecScatter surfaceScatter;
    Vec lambdaVec;

    FrontAreaArcLengthControl(FieldInterface& _mField,ConfigurationalFractureMechanics *_conf_prob,ArcLengthCtx *_arc_ptr);
    ~FrontAreaArcLengthControl();

    double aRea,aRea0,lambda_int;
    double resSpatialNrm2,resCrackFrontNrm2;

    PetscErrorCode set_dlambda_to_x(Vec x,double dlambda);
    PetscErrorCode calculate_area();
    PetscErrorCode calculate_lambda_int();
    PetscErrorCode calculate_db();
    PetscErrorCode get_dlambda(Vec x);

    PetscErrorCode preProcess();
    PetscErrorCode operator()();
    PetscErrorCode postProcess();

  };

};

PetscErrorCode SNESMonitorSpatialAndSmoothing_FEMEthod(SNES snes,PetscInt its,PetscReal fgnorm,void *dummy);


PetscErrorCode main_spatial_solution(FieldInterface& m_field,ConfigurationalFractureMechanics& conf_prob);
PetscErrorCode main_material_forces(FieldInterface& m_field,ConfigurationalFractureMechanics& conf_prob);

//crack propagation

/** \brief rescale load factor, such that maximally stressed crack front node has griffithe energy equal to gc
  *
  */
PetscErrorCode main_rescale_load_factor(FieldInterface& m_field,ConfigurationalFractureMechanics& conf_prob);

PetscErrorCode main_arc_length_setup(FieldInterface& m_field,ConfigurationalFractureMechanics& conf_prob);
PetscErrorCode main_arc_length_restart(FieldInterface& m_field,ConfigurationalFractureMechanics& conf_prob);
PetscErrorCode main_arc_length_solve(FieldInterface& m_field,ConfigurationalFractureMechanics& conf_prob);


#endif //__CONFIGURATIONAL_MECHANICS_HPP__
