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

#ifndef __CONFIGURATIONAL_MECHANICS_HPP__
#define __CONFIGURATIONAL_MECHANICS_HPP__

#include "FieldInterface.hpp"

#include "FEMethod_DriverComplexForLazy.hpp"
#include "PostProcNonLinearElasticityStresseOnRefindedMesh.hpp"

using namespace MoFEM;

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

  matPROJ_ctx *projSurfaceCtx,*projFrontCtx;

  BitRefLevel *ptr_bit_level0;
  ConfigurationalFractureMechanics(FieldInterface& mField): projSurfaceCtx(NULL),projFrontCtx(NULL) {

    ErrorCode rval;
    Tag th_my_ref_level;
    rval = mField.get_moab().tag_get_handle("_MY_REFINMENT_LEVEL",th_my_ref_level); CHKERR_THROW(rval);
    const EntityHandle root_meshset = mField.get_moab().get_root_set();
    rval = mField.get_moab().tag_get_by_ptr(th_my_ref_level,&root_meshset,1,(const void**)&ptr_bit_level0); CHKERR_THROW(rval);

    fe_post_proc_stresses_method = NULL;

  };

  ~ConfigurationalFractureMechanics() {
    if(fe_post_proc_stresses_method!=NULL) {
      delete fe_post_proc_stresses_method;
    }
  }
  
  PetscErrorCode set_material_fire_wall(FieldInterface& mField);
  PetscErrorCode thermal_field(FieldInterface& mField);
  PetscErrorCode spatial_problem_definition(FieldInterface& mField); 
  PetscErrorCode material_problem_definition(FieldInterface& mField);
  PetscErrorCode mesh_smoothing_problem_definition(FieldInterface& mField);
  PetscErrorCode coupled_problem_definition(FieldInterface& mField);
  PetscErrorCode arclength_problem_definition(FieldInterface& mField);
  PetscErrorCode constrains_problem_definition(FieldInterface& mField);
  PetscErrorCode constrains_crack_front_problem_definition(FieldInterface& mField,string problem);
  PetscErrorCode spatial_partition_problems(FieldInterface& mField);
  PetscErrorCode material_partition_problems(FieldInterface& mField);
  PetscErrorCode coupled_partition_problems(FieldInterface& mField);
  PetscErrorCode constrains_partition_problems(FieldInterface& mField,string problem);
  PetscErrorCode crackfront_partition_problems(FieldInterface& mField,string problem);
  PetscErrorCode set_spatial_positions(FieldInterface& mField);
  PetscErrorCode set_material_positions(FieldInterface& mField);
  PetscErrorCode set_coordinates_from_material_solution(FieldInterface& mField);

  PostProcStressNonLinearElasticity *fe_post_proc_stresses_method;
  PetscErrorCode solve_spatial_problem(FieldInterface& mField,SNES *snes,bool postproc = true);
  PetscErrorCode solve_material_problem(FieldInterface& mField,SNES *snes);
  PetscErrorCode solve_mesh_smooting_problem(FieldInterface& mField,SNES *snes);

  double aRea,lambda,energy;
  int nb_un_freez_nodes;
  bool freeze_all_but_one;
  PetscErrorCode solve_coupled_problem(FieldInterface& mField,SNES *snes,const double da,const double fraction_treshold = 5e-2);

  PetscErrorCode calculate_material_forces(FieldInterface& mField,string problem,string fe);
  PetscErrorCode surface_projection_data(FieldInterface& mField,string problem);
  PetscErrorCode delete_surface_projection_data(FieldInterface& mField);
  PetscErrorCode project_force_vector(FieldInterface& mField,string problem);
  PetscErrorCode front_projection_data(FieldInterface& mField,string problem);
  PetscErrorCode delete_front_projection_data(FieldInterface& mField);
  PetscErrorCode griffith_force_vector(FieldInterface& mField,string problem);

  PetscErrorCode project_form_th_projection_tag(FieldInterface& mField,string problem,bool do_not_project = false);

  map<EntityHandle,double> map_ent_g,map_ent_j;
  PetscScalar ave_g,min_g,max_g;
  PetscScalar ave_j,min_j,max_j;
  PetscErrorCode griffith_g(FieldInterface& mField,string problem);

  struct CubitDisplacementDirihletBC_Coupled: public CubitDisplacementDirihletBC {
  
    Range& cornersNodes;

    bool fixAllSpatialDispacements;
    CubitDisplacementDirihletBC_Coupled(FieldInterface& _mField,const string _problem_name,
      Range &_cornersNodes): 
      CubitDisplacementDirihletBC(_mField,_problem_name,"None"),
      cornersNodes(_cornersNodes),fixAllSpatialDispacements(false) {}

    PetscErrorCode SetDirihletBC_to_ElementIndicies(
      FieldInterface::FEMethod *fe_method_ptr,vector<vector<DofIdx> > &RowGlobDofs,vector<vector<DofIdx> > &ColGlobDofs,vector<DofIdx>& DirihletBC);
    PetscErrorCode SetDirihletBC_to_ElementIndiciesRow(
      FieldInterface::FEMethod *fe_method_ptr,vector<vector<DofIdx> > &RowGlobDofs,vector<DofIdx>& DirihletBC);  
    PetscErrorCode SetDirihletBC_to_ElementIndiciesCol(
      FieldInterface::FEMethod *fe_method_ptr,vector<vector<DofIdx> > &ColGlobDofs,vector<DofIdx>& DirihletBC);
    PetscErrorCode SetDirihletBC_to_MatrixDiagonal(FieldInterface::FEMethod *fe_method_ptr,Mat Aij);
    PetscErrorCode SetDirihletBC_to_RHS(FieldInterface::FEMethod *fe_method_ptr,Vec F);
    PetscErrorCode SetDirihletBC_to_ElementIndiciesFace(FieldInterface::FEMethod *fe_method_ptr,
      vector<DofIdx>& DirihletBC,vector<DofIdx>& FaceNodeGlobalDofs,vector<vector<DofIdx> > &FaceEdgeGlobalDofs,vector<DofIdx> &FaceGlobalDofs);
  
  };

 struct ArcLengthElemFEMethod: public FieldInterface::FEMethod {

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

    ArcLengthElemFEMethod(FieldInterface& _mField,ConfigurationalFractureMechanics *_conf_prob,ArcLengthCtx *_arc_ptr);
    ~ArcLengthElemFEMethod();

    double aRea,aRea0,lambda_int;

    PetscErrorCode set_dlambda_to_x(Vec x,double dlambda);
    PetscErrorCode calulate_area();
    PetscErrorCode calulate_lambda_int();
    PetscErrorCode calulate_db();
    PetscErrorCode get_dlambda(Vec x);

    PetscErrorCode preProcess();
    PetscErrorCode operator()();
    PetscErrorCode postProcess();

  };

};

PetscErrorCode SNESMonitorSpatialAndSmoothing_FEMEthod(SNES snes,PetscInt its,PetscReal fgnorm,void *dummy);


PetscErrorCode main_spatial_solution(FieldInterface& mField,ConfigurationalFractureMechanics& conf_prob);
PetscErrorCode main_material_forces(FieldInterface& mField,ConfigurationalFractureMechanics& conf_prob);

//crack propagation 

/** \brief rescale load factor, such that maximally stressed crack front node has griffithe energy equal to gc
  *
  */
PetscErrorCode main_rescale_load_factor(FieldInterface& mField,ConfigurationalFractureMechanics& conf_prob);

PetscErrorCode main_arc_length_setup(FieldInterface& mField,ConfigurationalFractureMechanics& conf_prob);
PetscErrorCode main_arc_length_solve(FieldInterface& mField,ConfigurationalFractureMechanics& conf_prob,bool face_splitting = false);

//face splitting
PetscErrorCode main_face_splitting_restart(FieldInterface& mField,ConfigurationalFractureMechanics& conf_prob);


#endif //__CONFIGURATIONAL_MECHANICS_HPP__
