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
#include "SnesCtx.hpp"
#include "ArcLengthTools.hpp"

using namespace MoFEM;

struct ConfigurationalFractureMechanics {
 
  Tag th_MaterialFireWall;
  typedef bitset<16> Material_FirelWall_def;
  Material_FirelWall_def *material_FirelWall;

  enum FirWall {
    FW_spatial_problem_definition = 0,
    FW_material_problem_definition,
    FW_coupled_problem_definition,
    FW_constrains_problem_definition,
    FW_constrains_crack_front_problem_definition,
    FW_set_spatial_positions,
    FW_set_material_positions,
    FW_arc_lenhghat_definition
  };

  EntityHandle cornersNodesMeshset,surfacesFacesNodesMeshset,crackSurfacesFacesNodesMeshset,crackForntMeshset;
  matPROJ_ctx *projSurfaceCtx,*projFrontCtx;

  BitRefLevel *ptr_bit_level0;
  BitRefLevel bit_level0;
  ConfigurationalFractureMechanics(FieldInterface& mField): projSurfaceCtx(NULL),projFrontCtx(NULL) {

    ErrorCode rval;
    Tag th_my_ref_level;
    rval = mField.get_moab().tag_get_handle("_MY_REFINMENT_LEVEL",th_my_ref_level); CHKERR_THROW(rval);
    const EntityHandle root_meshset = mField.get_moab().get_root_set();
    rval = mField.get_moab().tag_get_by_ptr(th_my_ref_level,&root_meshset,1,(const void**)&ptr_bit_level0); CHKERR_THROW(rval);
    bit_level0 = *ptr_bit_level0;

  };
  
  PetscErrorCode set_material_fire_wall(FieldInterface& mField);
  PetscErrorCode spatial_problem_definition(FieldInterface& mField); 
  PetscErrorCode material_problem_definition(FieldInterface& mField);
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
  PetscErrorCode solve_spatial_problem(FieldInterface& mField,SNES *snes,double step_size);
  PetscErrorCode solve_material_problem(FieldInterface& mField,SNES *snes);

  double nrm2_front_equlibrium;
  PetscErrorCode solve_coupled_problem(FieldInterface& mField,SNES *snes,double step_size,double alpha3);

  PetscErrorCode calculate_spatial_residual(FieldInterface& mField);
  PetscErrorCode calculate_material_forces(FieldInterface& mField,string problem,string fe);
  PetscErrorCode surface_projection_data(FieldInterface& mField,string problem);
  PetscErrorCode project_force_vector(FieldInterface& mField,string problem);
  PetscErrorCode front_projection_data(FieldInterface& mField,string problem);
  PetscErrorCode griffith_force_vector(FieldInterface& mField,string problem);

  PetscScalar ave_g,min_g,max_g;
  PetscErrorCode griffith_g(FieldInterface& mField,string problem);

  struct CubitDisplacementDirihletBC_Coupled: public CubitDisplacementDirihletBC {
  
    CubitDisplacementDirihletBC_Coupled (FieldInterface& _mField,const string _problem_name): 
      CubitDisplacementDirihletBC(_mField,_problem_name,"None") {}
  
    PetscErrorCode SetDirihletBC_to_ElementIndiciesRow(
      FieldInterface::FEMethod *fe_method_ptr,vector<vector<DofIdx> > &RowGlobDofs,vector<DofIdx>& DirihletBC);  
    PetscErrorCode SetDirihletBC_to_ElementIndiciesCol(
      FieldInterface::FEMethod *fe_method_ptr,vector<vector<DofIdx> > &ColGlobDofs,vector<DofIdx>& DirihletBC);
  
  };

  struct ConstrainCrackForntEdges_FEMethod: public FieldInterface::FEMethod {

    FieldInterface& mField;
    ConfigurationalFractureMechanics *conf_prob;

    double alpha3;
    ublas::vector<DofIdx> rowDofs,colDofs;
    ublas::vector<FieldData> dofsX,coords;
    ublas::vector<FieldData> delta0,delta;
    ublas::vector<FieldData> f;
    ublas::matrix<FieldData> K;

    ConstrainCrackForntEdges_FEMethod(
      FieldInterface& _mField,ConfigurationalFractureMechanics *_conf_prob,double _alpha3): mField(_mField),conf_prob(_conf_prob),alpha3(_alpha3) {
      rowDofs.resize(6);
      colDofs.resize(6);
      dofsX.resize(6);
      coords.resize(6);
      delta0.resize(3);
      delta.resize(3);
      f.resize(6);
      K.resize(6,6);
    }
    

    PetscErrorCode preProcess();
    PetscErrorCode operator()();
    PetscErrorCode postProcess();

  };

  struct ArcLengthElemFEMethod: public FieldInterface::FEMethod {

    FieldInterface& mField;
    ConfigurationalFractureMechanics *conf_prob;
    ArcLengthCtx* arc_ptr;
    Vec GhostDiag;

    Range crackSurfacesFaces;
    PetscInt *isIdx;
    IS isSurface;
    Vec surfaceDofs;
    VecScatter surfaceScatter;

    ArcLengthElemFEMethod(FieldInterface& _mField,ConfigurationalFractureMechanics *_conf_prob,ArcLengthCtx *_arc_ptr);
    ~ArcLengthElemFEMethod();

    PetscErrorCode set_dlambda_to_x(Vec x,double dlambda);
    double Area;
    PetscErrorCode calulate_lambda_int(double &_lambda_int_);
    PetscErrorCode calulate_db();
    PetscErrorCode calulate_dx_and_dlambda(Vec x);

    double lambda_int;

    PetscErrorCode preProcess();
    PetscErrorCode operator()();
    PetscErrorCode postProcess();

  };

};

PetscErrorCode  SNESMonitorSpatialAndSmoothing_FEMEthod(SNES snes,PetscInt its,PetscReal fgnorm,void *dummy);


#endif //__CONFIGURATIONAL_MECHANICS_HPP__
