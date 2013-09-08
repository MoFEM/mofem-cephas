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

#ifndef __HOOK_MATERIAL_FORCES_HPP__
#define __HOOK_MATERIAL_FORCES_HPP__

#include "moabField.hpp"
#include "moabField_Core.hpp"
#include "moabFEMethod_UpLevelStudent.hpp"
#include "cholesky.hpp"
#include <petscksp.h>

#include "moabSnes.hpp"
#include "moabFEMethod_ComplexForLazy.hpp"
#include "moabFEMethod_DriverComplexForLazy.hpp"
#include "petscShellMATs_ConstrainsByMarkAinsworth.hpp"
#include "moabFEMethod_SurfaceConstrains.hpp"
#include "ElasticFEMethod.hpp"

#include "complex_for_lazy.h"

namespace MoFEM {

struct SetPositionsEntMethod: public moabField::EntMethod {
    ErrorCode rval;
    PetscErrorCode ierr;
    Interface& moab;

    EntityHandle node;
    double coords[3];

    SetPositionsEntMethod(Interface& _moab): EntMethod(),moab(_moab),node(no_handle) {}
    
    PetscErrorCode preProcess() {
      PetscFunctionBegin;
      PetscPrintf(PETSC_COMM_WORLD,"Start Set Positions\n");
      PetscFunctionReturn(0);
    } 
     
    PetscErrorCode operator()() {
      PetscFunctionBegin;
      if(dof_ptr->get_ent_type()!=MBVERTEX) PetscFunctionReturn(0);
      EntityHandle ent = dof_ptr->get_ent();
      int dof_rank = dof_ptr->get_dof_rank();
      double &fval = dof_ptr->get_FieldData();
      if(node!=ent) {
	rval = moab.get_coords(&ent,1,coords); CHKERR_PETSC(rval);
	node = ent;
      }
      fval = coords[dof_rank];
      PetscFunctionReturn(0);
    }

    PetscErrorCode postProcess() {
      PetscFunctionBegin;
      PetscPrintf(PETSC_COMM_WORLD,"End Set Positions\n");
      PetscFunctionReturn(0);
    }

};

struct Spatial_ElasticFEMethod: public FEMethod_DriverComplexForLazy_Spatial {

  Range& SideSet2;

  Spatial_ElasticFEMethod(Interface& _moab,BaseDirihletBC *_dirihlet_bc_method_ptr,double _lambda,double _mu,Range &_SideSet2,int _verbose = 0): 
      FEMethod_DriverComplexForLazy_Spatial(_moab,_dirihlet_bc_method_ptr,_lambda,_mu,_verbose), SideSet2(_SideSet2)  {

    set_PhysicalEquationNumber(hooke);
    //set_PhysicalEquationNumber(neohookean);

  }

  PetscErrorCode operator()() {
    PetscFunctionBegin;
    ierr = FEMethod_DriverComplexForLazy_Spatial::operator()(SideSet2); CHKERRQ(ierr);
    PetscFunctionReturn(0);
  }

};

struct MaterialForcesFEMethod: public FEMethod_DriverComplexForLazy_Material {

  Range& SideSet2;
  Vec F_MATERIAL;
  MaterialForcesFEMethod(Interface& _moab,BaseDirihletBC *_dirihlet_bc_method_ptr,double _lambda,double _mu,Vec _F_MATERIAL,Range& _SideSet2,int _verbose = 0): 
      FEMethod_DriverComplexForLazy_Material(_moab,_dirihlet_bc_method_ptr,_lambda,_mu,_verbose),F_MATERIAL(_F_MATERIAL),SideSet2(_SideSet2) {

    set_PhysicalEquationNumber(hooke);
    type_of_analysis = material_analysis;

  }

  PetscErrorCode preProcess() {
    PetscFunctionBegin;
    PetscFunctionReturn(0);
  }

  PetscErrorCode operator()() {
    PetscFunctionBegin;
    ierr = OpComplexForLazyStart(); CHKERRQ(ierr);
    ierr = GetIndices(RowGlobMaterial,ColGlobMaterial,material_field_name); CHKERRQ(ierr);
    ierr = GetData(dofs_x_edge_data,dofs_x_edge,
      dofs_x_face_data,dofs_x_face,
      dofs_x_volume,dofs_x,
      spatial_field_name); CHKERRQ(ierr);
    ierr = GetFint(); CHKERRQ(ierr);
    double t[] = { 0,0,t_val, 0,0,t_val, 0,0,t_val };
    ierr = CaluclateMaterialFext(F_MATERIAL,t,SideSet2); CHKERRQ(ierr);
    ierr = VecSetValues(F_MATERIAL,RowGlobMaterial[0].size(),&(RowGlobMaterial[0])[0],&(Fint_H.data())[0],ADD_VALUES); CHKERRQ(ierr);
    PetscFunctionReturn(0);
  }

  PetscErrorCode postProcess() {
    PetscFunctionBegin;
    PetscFunctionReturn(0);
  }

};

struct Material_ElasticFEMethod: public FEMethod_DriverComplexForLazy_Material {

  moabField& mField;
  matPROJ_ctx& proj_all_ctx;
  Range& SideSet2;
  bool init;
  Material_ElasticFEMethod(
      Interface& _moab,moabField& _mField,matPROJ_ctx &_proj_all_ctx,BaseDirihletBC *_dirihlet_bc_method_ptr,
      double _lambda,double _mu,Range &_SideSet2,int _verbose = 0): 
      FEMethod_DriverComplexForLazy_Material(_moab,_dirihlet_bc_method_ptr,_lambda,_mu,_verbose),mField(_mField),proj_all_ctx(_proj_all_ctx),
      SideSet2(_SideSet2),init(true)  {
    set_PhysicalEquationNumber(hooke);
    //set_PhysicalEquationNumber(neohookean);
  }

  Range CornersEdges,CornersNodes,SurfacesFaces;
  C_SURFACE_FEMethod *CFE_SURFACE;
  g_SURFACE_FEMethod *gFE_SURFACE;
  C_CORNER_FEMethod *CFE_CORNER;
  g_CORNER_FEMethod *gFE_CORNER;
  PetscErrorCode preProcess() {
    PetscFunctionBegin;

    if(init) {
      init = false;
      ierr = mField.get_Cubit_msId_entities_by_dimension(100,SideSet,1,CornersEdges,true); CHKERRQ(ierr);
      ierr = mField.get_Cubit_msId_entities_by_dimension(101,NodeSet,0,CornersNodes,true); CHKERRQ(ierr);
      ierr = mField.get_Cubit_msId_entities_by_dimension(102,SideSet,2,SurfacesFaces,true); CHKERRQ(ierr);
      Range CornersEdgesNodes;
      rval = moab.get_adjacencies(CornersEdges,0,false,CornersEdgesNodes,Interface::UNION); CHKERR_PETSC(rval);
      CornersNodes.insert(CornersEdgesNodes.begin(),CornersEdgesNodes.end());
      Range SideSet1;
      ierr = mField.get_Cubit_msId_entities_by_dimension(1,SideSet,2,SideSet1,true); CHKERRQ(ierr);
      Range SideSet1Nodes;
      rval = moab.get_connectivity(SideSet1,SideSet1Nodes,true); CHKERR_PETSC(rval);
      CornersNodes.insert(SideSet1Nodes.begin(),SideSet1Nodes.end());
      CFE_SURFACE = new C_SURFACE_FEMethod(moab,SurfacesFaces,proj_all_ctx.C);
      gFE_SURFACE = new g_SURFACE_FEMethod(moab,SurfacesFaces,proj_all_ctx.g);
      CFE_CORNER = new C_CORNER_FEMethod(moab,CornersNodes,proj_all_ctx.C);
      gFE_CORNER = new g_CORNER_FEMethod(moab,CornersNodes,proj_all_ctx.g);

      ierr = MatSetOption(proj_all_ctx.C,MAT_NEW_NONZERO_LOCATION_ERR,PETSC_TRUE); CHKERRQ(ierr);
      ierr = MatSetOption(proj_all_ctx.C,MAT_NEW_NONZERO_ALLOCATION_ERR,PETSC_TRUE); CHKERRQ(ierr);

    }

    ierr = MatZeroEntries(proj_all_ctx.C); CHKERRQ(ierr);
    ierr = mField.loop_finite_elements("C_ALL_MATRIX","C_SURFACE_ELEM",*CFE_SURFACE);  CHKERRQ(ierr);
    ierr = MatAssemblyBegin(proj_all_ctx.C,MAT_FLUSH_ASSEMBLY); CHKERRQ(ierr);
    ierr = MatAssemblyEnd(proj_all_ctx.C,MAT_FLUSH_ASSEMBLY); CHKERRQ(ierr);
    ierr = mField.loop_finite_elements("C_ALL_MATRIX","C_CORNER_ELEM",*CFE_CORNER);  CHKERRQ(ierr);
    ierr = MatAssemblyBegin(proj_all_ctx.C,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
    ierr = MatAssemblyEnd(proj_all_ctx.C,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
    //
    ierr = proj_all_ctx.RecalculateCTandCCT(); CHKERRQ(ierr);
    ierr = proj_all_ctx.RecalulateCTC(); CHKERRQ(ierr);

    //{
      //MatView(proj_all_ctx.C,PETSC_VIEWER_DRAW_WORLD);
      //std::string wait;
      //std::cin >> wait;
    //}

    ierr = VecZeroEntries(proj_all_ctx.g); CHKERRQ(ierr);
    ierr = mField.loop_finite_elements("C_ALL_MATRIX","C_SURFACE_ELEM",*gFE_SURFACE);  CHKERRQ(ierr);
    ierr = VecAssemblyBegin(proj_all_ctx.g); CHKERRQ(ierr);
    ierr = VecAssemblyEnd(proj_all_ctx.g); CHKERRQ(ierr);
    ierr = mField.loop_finite_elements("C_ALL_MATRIX","C_CORNER_ELEM",*gFE_CORNER);  CHKERRQ(ierr);
    ierr = VecAssemblyBegin(proj_all_ctx.g); CHKERRQ(ierr);
    ierr = VecAssemblyEnd(proj_all_ctx.g); CHKERRQ(ierr);
  
    PetscReal g_nrm2;
    ierr = VecNorm(proj_all_ctx.g, NORM_2,&g_nrm2); CHKERRQ(ierr);
    PetscPrintf(PETSC_COMM_WORLD,"\t g_nrm2 = %6.4e\n",g_nrm2);


    ierr = FEMethod_DriverComplexForLazy_Material::preProcess(); CHKERRQ(ierr);
    PetscFunctionReturn(0);
  }

  PetscErrorCode operator()() {
    PetscFunctionBegin;

    ierr = OpComplexForLazyStart(); CHKERRQ(ierr);
    ierr = GetIndicesMaterial(); CHKERRQ(ierr);
    ierr = GetData(dofs_x_edge_data,dofs_x_edge,
      dofs_x_face_data,dofs_x_face,
      dofs_x_volume,dofs_x,
      spatial_field_name); CHKERRQ(ierr);

    ierr = SetDirihletBC_to_ElementIndicies(); CHKERRQ(ierr);
    if(Diagonal!=PETSC_NULL) {
      if(DirihletBC.size()>0) {
	DirihletBCDiagVal.resize(DirihletBC.size());
	fill(DirihletBCDiagVal.begin(),DirihletBCDiagVal.end(),1);
	ierr = VecSetValues(Diagonal,DirihletBC.size(),&(DirihletBC[0]),&DirihletBCDiagVal[0],INSERT_VALUES); CHKERRQ(ierr);
      }
    }

    double t[] = { 0,0,t_val, 0,0,t_val, 0,0,t_val };

    switch(snes_ctx) {
      case ctx_SNESSetFunction: { 
	ierr = CalculateMaterialFint(snes_f); CHKERRQ(ierr);
	//ierr = CaluclateMaterialFext(snes_f,t,SideSet2) ; CHKERRQ(ierr);
      }
      break;
      case ctx_SNESSetJacobian:
	ierr = CalculateMaterialTangent(proj_all_ctx.K); CHKERRQ(ierr);
	//ierr = CalculateMaterialTangentExt(*snes_B,t,SideSet2); CHKERRQ(ierr);
	break;
      default:
	SETERRQ(PETSC_COMM_SELF,1,"not implemented");
    }


    PetscFunctionReturn(0);
  }

  PetscErrorCode postProcess() {
    PetscFunctionBegin;

    ierr = FEMethod_DriverComplexForLazy_Material::postProcess(); CHKERRQ(ierr);

    switch(snes_ctx) {
      case ctx_SNESSetFunction: { 
	ierr = VecGhostUpdateBegin(snes_f,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
	ierr = VecGhostUpdateEnd(snes_f,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
	ierr = VecAssemblyBegin(snes_f); CHKERRQ(ierr);
	ierr = VecAssemblyEnd(snes_f); CHKERRQ(ierr);
	Vec Qf;
	ierr = VecDuplicate(snes_f,&Qf); CHKERRQ(ierr);
	int M,N,m,n;
	ierr = MatGetSize(proj_all_ctx.K,&M,&N); CHKERRQ(ierr);
	ierr = MatGetLocalSize(proj_all_ctx.K,&m,&n); CHKERRQ(ierr);
	Mat Q,R;
	ierr = MatCreateShell(PETSC_COMM_WORLD,m,n,M,N,&proj_all_ctx,&Q); CHKERRQ(ierr);
	ierr = MatShellSetOperation(Q,MATOP_MULT,(void(*)(void))matQ_mult_shell); CHKERRQ(ierr);
	int C_M,C_N,C_m,C_n;
	ierr = MatGetSize(proj_all_ctx.C,&C_M,&C_N); CHKERRQ(ierr);
	ierr = MatGetLocalSize(proj_all_ctx.C,&C_m,&C_n); CHKERRQ(ierr);
	ierr = MatCreateShell(PETSC_COMM_WORLD,m,C_m,M,C_m,&proj_all_ctx,&R); CHKERRQ(ierr);
	ierr = MatShellSetOperation(R,MATOP_MULT,(void(*)(void))matR_mult_shell); CHKERRQ(ierr);
	//Qf = R*g
	ierr = proj_all_ctx.InitQorP(snes_f); CHKERRQ(ierr);
	ierr = MatMult(R,proj_all_ctx.g,Qf); CHKERRQ(ierr);
	ierr = VecScale(Qf,-1); CHKERRQ(ierr);
	PetscReal Rg_nrm2;
	ierr = VecNorm(Qf,NORM_2,&Rg_nrm2); CHKERRQ(ierr);
	PetscPrintf(PETSC_COMM_WORLD,"\t Rg_nrm2 = %6.4e\n",Rg_nrm2);
	//snes_f = snes_f+K*R*g < snes_f += K*Qf >
	ierr = MatMultAdd(proj_all_ctx.K,Qf,snes_f,snes_f); CHKERRQ(ierr);
	//Qf = QT*(f + K*R*Qf) < Qf = Q*snes_f >
	ierr = MatMult(Q,snes_f,Qf); CHKERRQ(ierr);
	//Qf = CT*g + QT*(f+K*R*g) < Qf += Qf >
	ierr = MatMultAdd(proj_all_ctx.CT,proj_all_ctx.g,Qf,Qf); CHKERRQ(ierr);
	//Swap vectors Qf <--> snes_f
	ierr = VecSwap(snes_f,Qf); CHKERRQ(ierr);
	ierr = VecDestroy(&Qf); CHKERRQ(ierr);
	ierr = MatDestroy(&Q); CHKERRQ(ierr);
	ierr = MatDestroy(&R); CHKERRQ(ierr);
      }
      break;
      case ctx_SNESSetJacobian:
	ierr = proj_all_ctx.InitQTKQ(); CHKERRQ(ierr);
	ierr = MatAssemblyBegin(proj_all_ctx.K,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
	ierr = MatAssemblyEnd(proj_all_ctx.K,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
	ierr = MatCopy(proj_all_ctx.K,*snes_B,DIFFERENT_NONZERO_PATTERN); CHKERRQ(ierr);
	ierr = MatAXPY(*snes_B,1.,proj_all_ctx.CTC,DIFFERENT_NONZERO_PATTERN); CHKERRQ(ierr);
	break;
      default:
	SETERRQ(PETSC_COMM_SELF,1,"not implemented");
    }


    PetscFunctionReturn(0);
  }

};


}

#endif //__HOOK_MATERIAL_FORCES_HPP__
