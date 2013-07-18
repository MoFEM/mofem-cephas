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

#include "moabField.hpp"
#include "moabField_Core.hpp"
#include <petscksp.h>

#include "moabSnes.hpp"
#include "PostProcDisplacementOnMesh.hpp"
#include "PostProcDisplacementAndStrainOnRefindedMesh.hpp"

#include "nonlinear_elasticity.hpp"

using namespace MoFEM;

ErrorCode rval;
PetscErrorCode ierr;

struct MyElasticFEMethod: public FEMethod_DriverComplexForLazy {

  Vec F_lambda,b,db;
  Range& SideSet1;
  Range& SideSet2;
  Range SideSet1_;
  Range& SideSetArcLenght;
  Range SideSetArcLenght_;

  MyElasticFEMethod(Interface& _moab,double _lambda,double _mu,
      Vec _F_lambda,Vec _b,Vec _db,
      Range &_SideSet1,Range &_SideSet2,Range _SideSetArcLenght,
      int _verbose = 0): 
      FEMethod_DriverComplexForLazy(_moab,_lambda,_mu,_verbose), 
      F_lambda(_F_lambda),b(_b),db(_db),
      SideSet1(_SideSet1),SideSet2(_SideSet2),SideSetArcLenght(_SideSetArcLenght)  {

    set_PhysicalEquationNumber(neohookean);

    Range SideSet1Edges,SideSet1Nodes;
    rval = moab.get_adjacencies(SideSet1,1,false,SideSet1Edges,Interface::UNION); CHKERR_THROW(rval);
    rval = moab.get_connectivity(SideSet1,SideSet1Nodes,true); CHKERR_THROW(rval);
    SideSet1_.insert(SideSet1.begin(),SideSet1.end());
    SideSet1_.insert(SideSet1Edges.begin(),SideSet1Edges.end());
    SideSet1_.insert(SideSet1Nodes.begin(),SideSet1Nodes.end());
    //
    rval = moab.get_connectivity(SideSetArcLenght,SideSetArcLenght_,true); CHKERR_THROW(rval);
  }

  PetscErrorCode preProcess() {
    PetscFunctionBegin;
    ierr = FEMethod_DriverComplexForLazy::preProcess(); CHKERRQ(ierr);
    switch(ctx) {
      case ctx_SNESSetFunction: { 
      	ierr = VecZeroEntries(F_lambda); CHKERRQ(ierr);
	ierr = VecGhostUpdateBegin(F_lambda,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
	ierr = VecGhostUpdateEnd(F_lambda,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
      	ierr = VecZeroEntries(b); CHKERRQ(ierr);
	ierr = VecGhostUpdateBegin(b,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
	ierr = VecGhostUpdateEnd(b,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
      }
      break;
      case ctx_SNESSetJacobian: {
      }
      break;
      default:
	SETERRQ(PETSC_COMM_SELF,1,"not implemented");
    }
    PetscFunctionReturn(0);
  }


  PetscErrorCode operator()() {
    PetscFunctionBegin;

    ierr = OpComplexForLazyStart(); CHKERRQ(ierr);
    ierr = GetIndices(); CHKERRQ(ierr);

    Range& DirihletSideSet = SideSet1_;
    Range& NeumannSideSet = SideSet2;

    ierr = ApplyDirihletBC(DirihletSideSet); CHKERRQ(ierr);
    if(Diagonal!=PETSC_NULL) {
	if(DirihletBC.size()>0) {
	  DirihletBCDiagVal.resize(DirihletBC.size());
	  fill(DirihletBCDiagVal.begin(),DirihletBCDiagVal.end(),1);
	  ierr = VecSetValues(Diagonal,DirihletBC.size(),&(DirihletBC[0]),&DirihletBCDiagVal[0],INSERT_VALUES); CHKERRQ(ierr);
	}
    }

    switch(ctx) {
      case ctx_SNESSetFunction: { 
	  ierr = CalculateFint(snes_f); CHKERRQ(ierr);
	  double t[] = { 0,0,t_val, 0,0,t_val, 0,0,t_val };
	  ierr = CaluclateFext(F_lambda,t,NeumannSideSet); CHKERRQ(ierr);
	}
	break;
      case ctx_SNESSetJacobian: {
	  ierr = CalculateTangent(*snes_B); CHKERRQ(ierr);
	  FENumeredDofMoFEMEntity_multiIndex::index<FieldName_mi_tag>::type::iterator dit,hi_dit;
	  dit = row_multiIndex->get<FieldName_mi_tag>().lower_bound("LAMBDA");
	  hi_dit = row_multiIndex->get<FieldName_mi_tag>().upper_bound("LAMBDA");
	  //only one LAMBDA
	  if(distance(dit,hi_dit)!=1) SETERRQ(PETSC_COMM_SELF,1,"data inconsistency");
	  double _lambda_ = dit->get_FieldData();
	  //PetscPrintf(PETSC_COMM_WORLD,"snes _lambda_ = %6.4e\n",_lambda_);  
	  double t[] = { 0,0,t_val, 0,0,t_val, 0,0,t_val };
	  cblas_dscal(6,_lambda_,t,1);
	  ierr = CalculateTangentExt(*snes_B,t,NeumannSideSet); CHKERRQ(ierr);
	}
	break;
      default:
	SETERRQ(PETSC_COMM_SELF,1,"not implemented");
    }

    switch(ctx) {
      case ctx_SNESSetFunction: { 
  
	FENumeredDofMoFEMEntity_multiIndex::index<FieldName_mi_tag>::type::iterator dit,hi_dit;
	dit = row_multiIndex->get<FieldName_mi_tag>().lower_bound("SPATIAL_POSITION");
	hi_dit = row_multiIndex->get<FieldName_mi_tag>().upper_bound("SPATIAL_POSITION");
	for(;dit!=hi_dit;dit++) {
	  //if(dit->get_ent_type()!=MBVERTEX) continue;
	  //if(find(SideSetArcLenght_.begin(),SideSetArcLenght_.end(),dit->get_ent())==SideSetArcLenght.end()) continue;
	  //(x0+dx)*(x0+dx) - s = 0
	  //x0*x0 + 2x0*dx + dx*dx - s = 0
	  ierr = VecSetValue(b,dit->get_petsc_gloabl_dof_idx(),dit->get_FieldData(),INSERT_VALUES); CHKERRQ(ierr);
	  ierr = VecSetValue(db,dit->get_petsc_gloabl_dof_idx(),2.*dit->get_FieldData(),INSERT_VALUES); CHKERRQ(ierr);
	}
      }
      break;
      default:
      break;
    }
 
    PetscFunctionReturn(0);
  }

  PetscErrorCode postProcess() {
    PetscFunctionBegin;
    switch(ctx) {
      case ctx_SNESSetFunction: { 
	ierr = VecGhostUpdateBegin(F_lambda,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
	ierr = VecGhostUpdateEnd(F_lambda,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
	ierr = VecAssemblyBegin(F_lambda); CHKERRQ(ierr);
	ierr = VecAssemblyEnd(F_lambda); CHKERRQ(ierr);
	ierr = VecAssemblyBegin(b); CHKERRQ(ierr);
	ierr = VecAssemblyEnd(b); CHKERRQ(ierr);
	//ierr = VecView(b,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);
	ierr = VecAssemblyBegin(db); CHKERRQ(ierr);
	ierr = VecAssemblyEnd(db); CHKERRQ(ierr);
      }
      break;
      case ctx_SNESSetJacobian: {
	ierr = VecAssemblyBegin(Diagonal); CHKERRQ(ierr);
	ierr = VecAssemblyEnd(Diagonal); CHKERRQ(ierr);
	ierr = MatDiagonalSet(*snes_B,Diagonal,ADD_VALUES); CHKERRQ(ierr);
	ierr = VecDestroy(&Diagonal); CHKERRQ(ierr);
      }
      break;
      default:
	SETERRQ(PETSC_COMM_SELF,1,"not implemented");
    }
    PetscFunctionReturn(0);
  }

};

struct ArcLenghtElemFEMethod: public moabField::FEMethod {
  Interface& moab;

  double s;
  PetscErrorCode set_s(double _s) { 
    PetscFunctionBegin;
    s = _s;
    PetscFunctionReturn(0);
  }

  Vec GhostLambda;
  Vec F_lambda,b;
  ArcLenghtElemFEMethod(Interface& _moab,Vec _F_lambda,Vec _b): FEMethod(),moab(_moab),F_lambda(_F_lambda),b(_b) {
    PetscInt ghosts[1] = { 0 };
    ParallelComm* pcomm = ParallelComm::get_pcomm(&moab,MYPCOMM_INDEX);
    if(pcomm->rank() == 0) {
      VecCreateGhost(PETSC_COMM_WORLD,1,1,1,ghosts,&GhostLambda);
    } else {
      VecCreateGhost(PETSC_COMM_WORLD,0,1,1,ghosts,&GhostLambda);
    }
  }

  double s0;
  double b_dot_x;
  PetscErrorCode preProcess() {
    PetscFunctionBegin;
    switch(ctx) {
      case ctx_SNESSetFunction: { 
	ierr = VecDot(snes_x,b,&b_dot_x); CHKERRQ(ierr);
	PetscInt iter;
	ierr = SNESGetIterationNumber(snes,&iter); CHKERRQ(ierr);
	if(iter == 0) {
	  ierr = VecDot(snes_x,b,&s0); CHKERRQ(ierr);
	}
      }
      break;
      case ctx_SNESSetJacobian: {
      }
      break;
      default:
	SETERRQ(PETSC_COMM_SELF,1,"not implemented");
    }
    PetscFunctionReturn(0);
  }

  PetscErrorCode operator()() {
    PetscFunctionBegin;

    FENumeredDofMoFEMEntity_multiIndex::index<FieldName_mi_tag>::type::iterator dit,hi_dit;
    dit = row_multiIndex->get<FieldName_mi_tag>().lower_bound("LAMBDA");
    hi_dit = row_multiIndex->get<FieldName_mi_tag>().upper_bound("LAMBDA");
    //only one LAMBDA
    if(distance(dit,hi_dit)!=1) SETERRQ(PETSC_COMM_SELF,1,"data inconsistency");

    switch(ctx) {
      case ctx_SNESSetFunction: {
	double res_lambda,lambda;
	PetscScalar *array;
	res_lambda = b_dot_x - (s0 + s);
	ierr = VecSetValue(snes_f,dit->get_petsc_gloabl_dof_idx(),res_lambda,ADD_VALUES); CHKERRQ(ierr);
	//
	ierr = VecGetArray(snes_x,&array); CHKERRQ(ierr);
	lambda = array[dit->get_petsc_local_dof_idx()];
	ierr = VecRestoreArray(snes_x,&array); CHKERRQ(ierr);
	ierr = VecSetValue(GhostLambda,0,lambda,INSERT_VALUES); CHKERRQ(ierr);
	//PetscPrintf(PETSC_COMM_WORLD,"snes res_lambda = %6.4e, b_dot_d - s0 = %6.4e lambda = %6.4e\n",res_lambda,b_dot_x-s0,lambda);  
      }
      break; 
      case ctx_SNESSetJacobian: {
	ierr = MatSetValue(*snes_B,dit->get_petsc_gloabl_dof_idx(),dit->get_petsc_gloabl_dof_idx(),1,ADD_VALUES); CHKERRQ(ierr);
      }
      break;
      default:
	SETERRQ(PETSC_COMM_SELF,1,"no implemented");
    }
    
    PetscFunctionReturn(0);
  }

  PetscErrorCode postProcess() {
    PetscFunctionBegin;
    switch(ctx) {
      case ctx_SNESSetFunction: { 
	ierr = VecAssemblyBegin(GhostLambda); CHKERRQ(ierr);
	ierr = VecAssemblyEnd(GhostLambda); CHKERRQ(ierr);
	ierr = VecGhostUpdateBegin(GhostLambda,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
	ierr = VecGhostUpdateEnd(GhostLambda,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
	double *lambda;
	ierr = VecGetArray(GhostLambda,&lambda); CHKERRQ(ierr);
	ierr = VecAXPY(snes_f,*lambda,F_lambda); CHKERRQ(ierr);
	PetscPrintf(PETSC_COMM_WORLD,"\tlambda = %6.4e\n",*lambda);  
	ierr = VecRestoreArray(GhostLambda,&lambda); CHKERRQ(ierr);
	ierr = VecGhostUpdateBegin(snes_f,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
	ierr = VecGhostUpdateEnd(snes_f,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
	ierr = VecAssemblyBegin(snes_f); CHKERRQ(ierr);
	ierr = VecAssemblyEnd(snes_f); CHKERRQ(ierr);
      }
      break;
      case ctx_SNESSetJacobian: {
	ierr = MatAssemblyBegin(*snes_B,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
	ierr = MatAssemblyEnd(*snes_B,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
	//Matrix View
	//MatView(*snes_B,PETSC_VIEWER_DRAW_WORLD);//PETSC_VIEWER_STDOUT_WORLD);
	//std::string wait;
	//std::cin >> wait;
      }
      break;
      default:
	SETERRQ(PETSC_COMM_SELF,1,"not implemented");
    }
    PetscFunctionReturn(0);
  }

};

struct MatShellCtx {
  moabField& mField;

  Mat Aij;
  Vec F_lambda,db;
  MatShellCtx(moabField& _mField,Mat _Aij,Vec _F_lambda,Vec _db): mField(_mField),Aij(_Aij),F_lambda(_F_lambda),db(_db) {};
  PetscErrorCode set_lambda(Vec ksp_x,double *lambda,ScatterMode scattermode) {
    PetscFunctionBegin;
    const MoFEMProblem *problem_ptr;
    ierr = mField.get_problems_database("ELASTIC_MECHANICS",&problem_ptr); CHKERRQ(ierr);
    //get problem dofs
    NumeredDofMoFEMEntity_multiIndex &numered_dofs_rows = const_cast<NumeredDofMoFEMEntity_multiIndex&>(problem_ptr->numered_dofs_rows);
    NumeredDofMoFEMEntity_multiIndex::index<FieldName_mi_tag>::type::iterator dit;
    dit = numered_dofs_rows.get<FieldName_mi_tag>().find("LAMBDA");
    DofIdx lambda_dof_index = dit->get_petsc_local_dof_idx();
    int part = dit->part;
    if(lambda_dof_index!=-1) {
      PetscScalar *array;
      ierr = VecGetArray(ksp_x,&array); CHKERRQ(ierr);
      switch(scattermode) {
	case SCATTER_FORWARD:
	  *lambda = array[lambda_dof_index];
	  break;
	case SCATTER_REVERSE:
	  array[lambda_dof_index] = *lambda;
	  break;
	default:
	  SETERRQ(PETSC_COMM_SELF,1,"not implemented");
      }
      ierr = VecRestoreArray(ksp_x,&array); CHKERRQ(ierr);
    } 
    MPI_Bcast(lambda,1,MPI_DOUBLE,part,PETSC_COMM_WORLD);
    PetscFunctionReturn(0);
  }
  ~MatShellCtx() { }
  friend PetscErrorCode arc_lenght_mult_shell(Mat A,Vec x,Vec f);
};
PetscErrorCode arc_lenght_mult_shell(Mat A,Vec x,Vec f) {
  PetscFunctionBegin;
  void *void_ctx;
  MatShellGetContext(A,&void_ctx);
  MatShellCtx *ctx = (MatShellCtx*)void_ctx;
  ierr = MatMult(ctx->Aij,x,f); CHKERRQ(ierr);
  double db_dot_x;
  ierr = VecDot(ctx->db,x,&db_dot_x); CHKERRQ(ierr);
  ierr = ctx->set_lambda(f,&db_dot_x,SCATTER_REVERSE); CHKERRQ(ierr);
  double lambda;
  ierr = ctx->set_lambda(x,&lambda,SCATTER_FORWARD); CHKERRQ(ierr);
  //PetscPrintf(PETSC_COMM_WORLD,"mat_mult lambda = %6.4e\n",lambda);
  ierr = VecAXPY(f,lambda,ctx->F_lambda); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

struct PCShellCtx {
  PC pc;
  Mat ShellAij,Aij;
  Vec x_lambda;
  PCShellCtx(Mat _ShellAij,Mat _Aij,Vec _x_lambda): 
    ShellAij(_ShellAij),Aij(_Aij),x_lambda(_x_lambda) {
    PCCreate(PETSC_COMM_WORLD,&pc);
  }
  ~PCShellCtx() {
    PCDestroy(&pc);
  }
  friend PetscErrorCode pc_apply_arc_length(PC pc,Vec pc_f,Vec pc_x);
  friend PetscErrorCode pc_setup_arc_length(PC pc);
};
PetscErrorCode pc_apply_arc_length(PC pc,Vec pc_f,Vec pc_x) {
  PetscFunctionBegin;
  void *void_ctx;
  ierr = PCShellGetContext(pc,&void_ctx); CHKERRQ(ierr);
  PCShellCtx *PCCtx = (PCShellCtx*)void_ctx;
  void *void_MatCtx;
  MatShellGetContext(PCCtx->ShellAij,&void_MatCtx);
  MatShellCtx *MatCtx = (MatShellCtx*)void_MatCtx;
  double res_lambda;
  ierr = MatCtx->set_lambda(pc_f,&res_lambda,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = PCApply(PCCtx->pc,pc_f,pc_x); CHKERRQ(ierr);
  ierr = PCApply(PCCtx->pc,MatCtx->F_lambda,PCCtx->x_lambda); CHKERRQ(ierr);
  // b \dot pc_x - res_lambda = 0
  // pc_x = x_int + lambda*x_lambda
  // b \dot x_int + lambda*(b \dot x_lambda) - res_lambda = 0
  // lambda = (res_lambda - b \dot x_int)/(b \cdot x_lambda)
  double db_dot_x_int,db_dot_x_lambda;
  ierr = VecDot(MatCtx->db,pc_x,&db_dot_x_int); CHKERRQ(ierr);
  ierr = VecDot(MatCtx->db,PCCtx->x_lambda,&db_dot_x_lambda); CHKERRQ(ierr);
  double dlambda;
  dlambda = (res_lambda - db_dot_x_int)/db_dot_x_lambda;
  if(dlambda != dlambda) SETERRQ(PETSC_COMM_SELF,1,"db \\dot x_lambda = 0\nCheck constrint vector, SideSet, ect.");
  //PetscPrintf(PETSC_COMM_WORLD,"pc dlambda = %6.4e\n",dlambda);
  ierr = VecAXPY(pc_x,dlambda,PCCtx->x_lambda); CHKERRQ(ierr);
  ierr = MatCtx->set_lambda(pc_x,&dlambda,SCATTER_REVERSE); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
PetscErrorCode pc_setup_arc_length(PC pc) {
  PetscFunctionBegin;
  void *void_ctx;
  ierr = PCShellGetContext(pc,&void_ctx); CHKERRQ(ierr);
  PCShellCtx *ctx = (PCShellCtx*)void_ctx;
  ierr = PCSetFromOptions(ctx->pc); CHKERRQ(ierr);
  MatStructure flag;
  ierr = PCGetOperators(pc,&ctx->ShellAij,&ctx->Aij,&flag); CHKERRQ(ierr);
  ierr = PCSetOperators(ctx->pc,ctx->ShellAij,ctx->Aij,flag); CHKERRQ(ierr);
  ierr = PCSetUp(ctx->pc); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

static char help[] = "\
-my_file mesh file name\n\
-my_sr reduction of step size\n\
-my_ms maximal number of steps\n\n";

int main(int argc, char *argv[]) {

  PetscInitialize(&argc,&argv,(char *)0,help);

  Core mb_instance;
  Interface& moab = mb_instance;
  int rank;
  MPI_Comm_rank(PETSC_COMM_WORLD,&rank);

  PetscBool flg = PETSC_TRUE;
  char mesh_file_name[255];
  ierr = PetscOptionsGetString(PETSC_NULL,"-my_file",mesh_file_name,255,&flg); CHKERRQ(ierr);
  if(flg != PETSC_TRUE) {
    SETERRQ(PETSC_COMM_SELF,1,"*** ERROR -my_file (MESH FILE NEEDED)");
  }

  PetscScalar step_size_reduction;
  ierr = PetscOptionsGetReal(PETSC_NULL,"-my_sr",&step_size_reduction,&flg); CHKERRQ(ierr);
  if(flg != PETSC_TRUE) {
    step_size_reduction = 1.;
  }

  PetscInt max_steps;
  ierr = PetscOptionsGetInt(PETSC_NULL,"-my_ms",&max_steps,&flg); CHKERRQ(ierr);
  if(flg != PETSC_TRUE) {
    max_steps = 5;
  }
 
  const char *option;
  option = "";//"PARALLEL=BCAST;";//;DEBUG_IO";
  rval = moab.load_file(mesh_file_name, 0, option); CHKERR(rval); 
  ParallelComm* pcomm = ParallelComm::get_pcomm(&moab,MYPCOMM_INDEX);
  if(pcomm == NULL) pcomm =  new ParallelComm(&moab,PETSC_COMM_WORLD);

  //data stored on mesh for restart
  Tag th_step_size,th_step;
  double def_step_size = 5e-2;
  moab.tag_get_handle("_STEPSIZE",1,MB_TYPE_DOUBLE,th_step_size,MB_TAG_CREAT|MB_TAG_MESH,&def_step_size); 
  int def_step = 1;
  moab.tag_get_handle("_STEP",1,MB_TYPE_INTEGER,th_step,MB_TAG_CREAT|MB_TAG_MESH,&def_step); 
  const void* tag_data_step_size[1];
  EntityHandle root = moab.get_root_set();
  rval = moab.tag_get_by_ptr(th_step_size,&root,1,tag_data_step_size); CHKERR(rval);
  double& step_size = *(double *)tag_data_step_size[0];
  const void* tag_data_step[1];
  rval = moab.tag_get_by_ptr(th_step,&root,1,tag_data_step); CHKERR(rval);
  int& step= *(int *)tag_data_step[0];
  //end of data stored for restart
  step_size *= step_size_reduction;

  PetscLogDouble t1,t2;
  PetscLogDouble v1,v2;
  ierr = PetscGetTime(&v1); CHKERRQ(ierr);
  ierr = PetscGetCPUTime(&t1); CHKERRQ(ierr);

  moabField_Core core(moab);
  moabField& mField = core;

  BitRefLevel bit_level0;
  bit_level0.set(0);

  if(step == 1) {
    Range CubitSideSets_meshsets;
    ierr = mField.get_CubitBCType_meshsets(SideSet,CubitSideSets_meshsets); CHKERRQ(ierr);

    //ref meshset ref level 0
    ierr = mField.seed_ref_level_3D(0,0); CHKERRQ(ierr);
    EntityHandle meshset_level0;
    rval = moab.create_meshset(MESHSET_SET,meshset_level0); CHKERR_PETSC(rval);
    ierr = mField.seed_ref_level_3D(0,bit_level0); CHKERRQ(ierr);
    ierr = mField.refine_get_ents(bit_level0,meshset_level0); CHKERRQ(ierr);

    //Fields
    ierr = mField.add_field("SPATIAL_POSITION",H1,3); CHKERRQ(ierr);
    ierr = mField.add_field("LAMBDA",NoField,1); CHKERRQ(ierr);

    //FE
    ierr = mField.add_finite_element("ELASTIC"); CHKERRQ(ierr);
    ierr = mField.add_finite_element("ARC_LENGHT"); CHKERRQ(ierr);

    //Define rows/cols and element data
    ierr = mField.modify_finite_element_add_field_row("ELASTIC","SPATIAL_POSITION"); CHKERRQ(ierr);
    ierr = mField.modify_finite_element_add_field_row("ELASTIC","LAMBDA"); CHKERRQ(ierr);
    ierr = mField.modify_finite_element_add_field_col("ELASTIC","SPATIAL_POSITION"); CHKERRQ(ierr);
    ierr = mField.modify_finite_element_add_field_col("ELASTIC","LAMBDA"); CHKERRQ(ierr);
    ierr = mField.modify_finite_element_add_field_data("ELASTIC","SPATIAL_POSITION"); CHKERRQ(ierr);
    ierr = mField.modify_finite_element_add_field_data("ELASTIC","LAMBDA"); CHKERRQ(ierr);

    //Define rows/cols and element data
    ierr = mField.modify_finite_element_add_field_row("ARC_LENGHT","LAMBDA"); CHKERRQ(ierr);
    ierr = mField.modify_finite_element_add_field_col("ARC_LENGHT","LAMBDA"); CHKERRQ(ierr);
    //elem data
    ierr = mField.modify_finite_element_add_field_data("ARC_LENGHT","LAMBDA"); CHKERRQ(ierr);

    //define problems
    ierr = mField.add_problem("ELASTIC_MECHANICS"); CHKERRQ(ierr);

    //set finite elements for problems
    ierr = mField.modify_problem_add_finite_element("ELASTIC_MECHANICS","ELASTIC"); CHKERRQ(ierr);
    ierr = mField.modify_problem_add_finite_element("ELASTIC_MECHANICS","ARC_LENGHT"); CHKERRQ(ierr);

    //set refinment level for problem
    ierr = mField.modify_problem_ref_level_add_bit("ELASTIC_MECHANICS",bit_level0); CHKERRQ(ierr);

    //add entitities (by tets) to the field
    ierr = mField.add_ents_to_field_by_TETs(0,"SPATIAL_POSITION"); CHKERRQ(ierr);

    //add finite elements entities
    ierr = mField.add_ents_to_finite_element_EntType_by_bit_ref(bit_level0,"ELASTIC",MBTET); CHKERRQ(ierr);
  
    //this entity will carray data for this finite element
    EntityHandle meshset_FE_ARC_LENGHT;
    rval = moab.create_meshset(MESHSET_SET,meshset_FE_ARC_LENGHT); CHKERR_PETSC(rval);
    //get LAMBDA field meshset
    EntityHandle meshset_field_LAMBDA = mField.get_field_meshset("LAMBDA");
    //add LAMBDA field meshset to finite element ARC_LENGHT
    rval = moab.add_entities(meshset_FE_ARC_LENGHT,&meshset_field_LAMBDA,1); CHKERR_PETSC(rval);
    //add finite element ARC_LENGHT meshset to refinment database (all ref bit leveles)
    ierr = mField.seed_ref_level_MESHSET(meshset_FE_ARC_LENGHT,BitRefLevel().set()); CHKERRQ(ierr);
    //finally add created meshset to the ARC_LENGHT finite element
    ierr = mField.add_ents_to_finite_element_by_MESHSET(meshset_FE_ARC_LENGHT,"ARC_LENGHT"); CHKERRQ(ierr);

    //set app. order
    ierr = mField.set_field_order(0,MBTET,"SPATIAL_POSITION",4); CHKERRQ(ierr);
    ierr = mField.set_field_order(0,MBTRI,"SPATIAL_POSITION",4); CHKERRQ(ierr);
    ierr = mField.set_field_order(0,MBEDGE,"SPATIAL_POSITION",4); CHKERRQ(ierr);
    ierr = mField.set_field_order(0,MBVERTEX,"SPATIAL_POSITION",1); CHKERRQ(ierr);
  }

  //build field
  ierr = mField.build_fields(); CHKERRQ(ierr);

  //build finite elemnts
  ierr = mField.build_finite_elements(); CHKERRQ(ierr);

  //build adjacencies
  ierr = mField.build_adjacencies(bit_level0); CHKERRQ(ierr);

  //build problem
  ierr = mField.build_problems(); CHKERRQ(ierr);

  //partition
  ierr = mField.partition_problems("ELASTIC_MECHANICS"); CHKERRQ(ierr);
  ierr = mField.partition_finite_elements("ELASTIC_MECHANICS"); CHKERRQ(ierr);
  ierr = mField.partition_ghost_dofs("ELASTIC_MECHANICS"); CHKERRQ(ierr);

  //create matrices
  Vec F;
  ierr = mField.VecCreateGhost("ELASTIC_MECHANICS",Row,&F); CHKERRQ(ierr);
  Mat Aij;
  ierr = mField.MatCreateMPIAIJWithArrays("ELASTIC_MECHANICS",&Aij); CHKERRQ(ierr);
  Vec F_lambda,b,db,x_lambda;
  ierr = mField.VecCreateGhost("ELASTIC_MECHANICS",Col,&F_lambda); CHKERRQ(ierr);
  ierr = mField.VecCreateGhost("ELASTIC_MECHANICS",Row,&b); CHKERRQ(ierr);
  ierr = mField.VecCreateGhost("ELASTIC_MECHANICS",Row,&db); CHKERRQ(ierr);
  ierr = mField.VecCreateGhost("ELASTIC_MECHANICS",Row,&x_lambda); CHKERRQ(ierr);

  PetscInt M,N;
  ierr = MatGetSize(Aij,&M,&N); CHKERRQ(ierr);
  PetscInt m,n;
  MatGetLocalSize(Aij,&m,&n);
  MatShellCtx* MatCtx = new MatShellCtx(mField,Aij,F_lambda,db);
  Mat ShellAij;
  ierr = MatCreateShell(PETSC_COMM_WORLD,m,n,M,N,(void*)MatCtx,&ShellAij); CHKERRQ(ierr);
  ierr = MatShellSetOperation(ShellAij,MATOP_MULT,(void(*)(void))arc_lenght_mult_shell); CHKERRQ(ierr);

  if(step==1) {
    SetPositionsEntMethod set_positions(moab);
    ierr = mField.loop_dofs("ELASTIC_MECHANICS","SPATIAL_POSITION",Row,set_positions); CHKERRQ(ierr);
  }

  Range SideSet1,SideSet2;
  ierr = mField.get_Cubit_msId_entities_by_dimension(1,SideSet,2,SideSet1,true); CHKERRQ(ierr);
  ierr = mField.get_Cubit_msId_entities_by_dimension(2,SideSet,2,SideSet2,true); CHKERRQ(ierr);
  PetscPrintf(PETSC_COMM_WORLD,"Nb. faces in SideSet 1 : %u\n",SideSet1.size());
  PetscPrintf(PETSC_COMM_WORLD,"Nb. faces in SideSet 2 : %u\n",SideSet2.size());

  const double YoungModulus = 1;
  const double PoissonRatio = 0.25;
  MyElasticFEMethod MyFE(moab,
    LAMBDA(YoungModulus,PoissonRatio),MU(YoungModulus,PoissonRatio),
    F_lambda,b,db,
    SideSet1,SideSet2,SideSet2);

  ArcLenghtElemFEMethod MyArcMethod(moab,F_lambda,b);

  moabSnesCtx SnesCtx(mField,"ELASTIC_MECHANICS");
  
  SNES snes;
  ierr = SNESCreate(PETSC_COMM_WORLD,&snes); CHKERRQ(ierr);
  ierr = SNESSetApplicationContext(snes,&SnesCtx); CHKERRQ(ierr);
  ierr = SNESSetFunction(snes,F,SnesRhs,&SnesCtx); CHKERRQ(ierr);
  ierr = SNESSetJacobian(snes,ShellAij,Aij,SnesMat,&SnesCtx); CHKERRQ(ierr);
  ierr = SNESSetFromOptions(snes); CHKERRQ(ierr);

  KSP ksp;
  ierr = SNESGetKSP(snes,&ksp); CHKERRQ(ierr);
  PC pc;
  ierr = KSPGetPC(ksp,&pc); CHKERRQ(ierr);
  PCShellCtx* PCCtx = new PCShellCtx(Aij,ShellAij,x_lambda);
  ierr = PCSetType(pc,PCSHELL); CHKERRQ(ierr);
  ierr = PCShellSetContext(pc,PCCtx); CHKERRQ(ierr);
  ierr = PCShellSetApply(pc,pc_apply_arc_length); CHKERRQ(ierr);
  ierr = PCShellSetSetUp(pc,pc_setup_arc_length); CHKERRQ(ierr);

  moabSnesCtx::loops_to_do_type& loops_to_do_Rhs = SnesCtx.get_loops_to_do_Rhs();
  loops_to_do_Rhs.push_back(moabSnesCtx::loop_pair_type("ELASTIC",&MyFE));
  loops_to_do_Rhs.push_back(moabSnesCtx::loop_pair_type("ARC_LENGHT",&MyArcMethod));
  moabSnesCtx::loops_to_do_type& loops_to_do_Mat = SnesCtx.get_loops_to_do_Mat();
  loops_to_do_Mat.push_back(moabSnesCtx::loop_pair_type("ELASTIC",&MyFE));
  loops_to_do_Mat.push_back(moabSnesCtx::loop_pair_type("ARC_LENGHT",&MyArcMethod));


  Vec D;
  ierr = VecDuplicate(F,&D); CHKERRQ(ierr);
  ierr = mField.set_local_VecCreateGhost("ELASTIC_MECHANICS",Row,D,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = VecGhostUpdateBegin(D,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = VecGhostUpdateEnd(D,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);

  ierr = MyFE.set_t_val(1e-4); CHKERRQ(ierr);

  int its_d = 4;
  double gamma = 0.5;
  for(;step<max_steps; step++) {
    ierr = MyArcMethod.set_s(step_size); CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_WORLD,"Load Setp %D dlambda = %6.4e\n",step,step_size); CHKERRQ(ierr);
    ierr = SNESSolve(snes,PETSC_NULL,D); CHKERRQ(ierr);
    int its;
    ierr = SNESGetIterationNumber(snes,&its); CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_WORLD,"number of Newton iterations = %D\n",its); CHKERRQ(ierr);
    step_size *= pow((double)its_d/(double)its,gamma);
    //
    //Save data on mesh
    ierr = mField.set_global_VecCreateGhost("ELASTIC_MECHANICS",Row,D,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
    //
    PostProcDisplacementsEntMethod ent_method(moab,"SPATIAL_POSITION");
    ierr = mField.loop_dofs("ELASTIC_MECHANICS","SPATIAL_POSITION",Row,ent_method); CHKERRQ(ierr);
    //
    if(step % 1 == 0) {
      if(pcomm->rank()==0) {
	rval = moab.write_file("restart.h5m"); CHKERR_PETSC(rval);
	//
	EntityHandle out_meshset;
	rval = moab.create_meshset(MESHSET_SET,out_meshset); CHKERR_PETSC(rval);
	ierr = mField.problem_get_FE("ELASTIC_MECHANICS","ELASTIC",out_meshset); CHKERRQ(ierr);
	ostringstream sss;
	sss << "out_" << step << ".vtk";
	rval = moab.write_file(sss.str().c_str(),"VTK","",&out_meshset,1); CHKERR_PETSC(rval);
	rval = moab.delete_entities(&out_meshset,1); CHKERR_PETSC(rval);
      }
    }
  }
  
  ierr = VecGhostUpdateBegin(D,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = VecGhostUpdateEnd(D,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);

  //Save data on mesh
  ierr = mField.set_global_VecCreateGhost("ELASTIC_MECHANICS",Row,D,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
  //ierr = VecView(F,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);

  PostProcDisplacementsEntMethod ent_method(moab,"SPATIAL_POSITION");
  ierr = mField.loop_dofs("ELASTIC_MECHANICS","SPATIAL_POSITION",Row,ent_method); CHKERRQ(ierr);

  if(pcomm->rank()==0) {
    EntityHandle out_meshset;
    rval = moab.create_meshset(MESHSET_SET,out_meshset); CHKERR_PETSC(rval);
    ierr = mField.problem_get_FE("ELASTIC_MECHANICS","ELASTIC",out_meshset); CHKERRQ(ierr);
    rval = moab.write_file("out.vtk","VTK","",&out_meshset,1); CHKERR_PETSC(rval);
    rval = moab.delete_entities(&out_meshset,1); CHKERR_PETSC(rval);
  }

  PostProcFieldsAndGradientOnRefMesh fe_post_proc_method(moab);
  ierr = mField.loop_finite_elements("ELASTIC_MECHANICS","ELASTIC",fe_post_proc_method);  CHKERRQ(ierr);
  PetscSynchronizedFlush(PETSC_COMM_WORLD);
  if(pcomm->rank()==0) {
    rval = fe_post_proc_method.moab_post_proc.write_file("out_post_proc.vtk","VTK",""); CHKERR_PETSC(rval);
  }

  //detroy matrices
  ierr = VecDestroy(&F); CHKERRQ(ierr);
  ierr = VecDestroy(&D); CHKERRQ(ierr);
  ierr = MatDestroy(&Aij); CHKERRQ(ierr);
  ierr = VecDestroy(&F_lambda); CHKERRQ(ierr);
  ierr = VecDestroy(&b); CHKERRQ(ierr);
  ierr = VecDestroy(&db); CHKERRQ(ierr);
  ierr = VecDestroy(&x_lambda); CHKERRQ(ierr);
  ierr = MatDestroy(&ShellAij); CHKERRQ(ierr);
  ierr = SNESDestroy(&snes); CHKERRQ(ierr);
  delete MatCtx;
  delete PCCtx;

  ierr = PetscGetTime(&v2);CHKERRQ(ierr);
  ierr = PetscGetCPUTime(&t2);CHKERRQ(ierr);

  PetscSynchronizedPrintf(PETSC_COMM_WORLD,"Total Rank %d Time = %f CPU Time = %f\n",pcomm->rank(),v2-v1,t2-t1);
  PetscSynchronizedFlush(PETSC_COMM_WORLD);

  ierr = PetscFinalize(); CHKERRQ(ierr);

  return 0;

}



