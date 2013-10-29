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

#ifndef __ARCLEGHTTOOLS_HPP__
#define __ARCLEGHTTOOLS_HPP__

namespace MoFEM {

struct ArcLenghtCtx_DataOnMesh {
  ErrorCode rval;
  PetscErrorCode ierr;

  Interface &moab;
  const void* tag_data_dlambda[1];
  ArcLenghtCtx_DataOnMesh(FieldInterface &mField): moab(mField.get_moab()) {
    Tag th_dlambda;
    double def_dlambda = 0;
    rval = moab.tag_get_handle("_DLAMBDA",1,MB_TYPE_DOUBLE,th_dlambda,MB_TAG_CREAT|MB_TAG_MESH,&def_dlambda);  
    if(rval==MB_ALREADY_ALLOCATED) rval = MB_SUCCESS;
    CHKERR_THROW(rval);
    EntityHandle root = moab.get_root_set();
    rval = moab.tag_get_by_ptr(th_dlambda,&root,1,tag_data_dlambda); CHKERR_THROW(rval);
  }
};



struct ArcLenghtCtx: public ArcLenghtCtx_DataOnMesh {

  ErrorCode rval;
  PetscErrorCode ierr;

  double& dlambda; //reference to moab data see ArcLenghtCtc_DataOnMesh and constructor ArcLenghtCtx
  PetscErrorCode set_dlambda(double _dlambda) {
    PetscFunctionBegin;
    dlambda = _dlambda;
    PetscFunctionReturn(0);
  } 

  double s,beta,alpha;
  PetscErrorCode set_s(double _s) { 
    PetscFunctionBegin;
    s = _s;
    ierr = PetscPrintf(PETSC_COMM_WORLD,"\tSet s = %6.4e\n",s); CHKERRQ(ierr);
    PetscFunctionReturn(0);
  }

  PetscErrorCode set_alpha_and_beta(double _alpha,double _beta) { 
    PetscFunctionBegin;
    alpha = _alpha;
    beta = _beta;
    ierr = PetscPrintf(PETSC_COMM_WORLD,"\tSet alpha = %6.4e beta = %6.4e\n",alpha,beta); CHKERRQ(ierr);
    PetscFunctionReturn(0);
  }


  double diag,dx2,F_lambda2,res_lambda,;
  Vec F_lambda,db,x_lambda,y_residual,x0,dx;
  ArcLenghtCtx(FieldInterface &mField,const string &problem_name): ArcLenghtCtx_DataOnMesh(mField), dlambda(*(double *)tag_data_dlambda[0]) {

    mField.VecCreateGhost(problem_name,Row,&F_lambda);
    VecSetOption(F_lambda,VEC_IGNORE_NEGATIVE_INDICES,PETSC_TRUE); 
    mField.VecCreateGhost(problem_name,Row,&db);
    mField.VecCreateGhost(problem_name,Row,&x_lambda);
    mField.VecCreateGhost(problem_name,Row,&y_residual);
    mField.VecCreateGhost(problem_name,Row,&x0);
    mField.VecCreateGhost(problem_name,Row,&dx);

  }

  ~ArcLenghtCtx() {
    VecDestroy(&F_lambda);
    VecDestroy(&db);
    VecDestroy(&x_lambda);
    VecDestroy(&y_residual);
    VecDestroy(&x0);
    VecDestroy(&dx);
  }


};

struct ArcLenghtSnesCtx: public SnesCtx {

  ErrorCode rval;
  PetscErrorCode ierr;

  ArcLenghtCtx* arc_ptr;
  ArcLenghtSnesCtx(FieldInterface &_mField,const string &_problem_name,ArcLenghtCtx* _arc_ptr):
    SnesCtx(_mField,_problem_name),arc_ptr(_arc_ptr) {}

};

struct MatShellCtx {

  ErrorCode rval;
  PetscErrorCode ierr;


  double scale_lambda;
  FieldInterface& mField;

  Mat Aij;
  ArcLenghtCtx* arc_ptr;
  MatShellCtx(FieldInterface& _mField,Mat _Aij,ArcLenghtCtx *_arc_ptr): scale_lambda(1),mField(_mField),Aij(_Aij),arc_ptr(_arc_ptr) {};
  PetscErrorCode set_lambda(Vec ksp_x,double *lambda,ScatterMode scattermode) {
    PetscFunctionBegin;
    const MoFEMProblem *problem_ptr;
    ierr = mField.get_problem("ELASTIC_MECHANICS",&problem_ptr); CHKERRQ(ierr);
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
  PetscErrorCode ierr;
  void *void_ctx;
  ierr = MatShellGetContext(A,&void_ctx); CHKERRQ(ierr);
  MatShellCtx *ctx = (MatShellCtx*)void_ctx;
  ierr = MatMult(ctx->Aij,x,f); CHKERRQ(ierr);
  double lambda;
  ierr = ctx->set_lambda(x,&lambda,SCATTER_FORWARD); CHKERRQ(ierr);
  double db_dot_x;
  ierr = VecDot(ctx->arc_ptr->db,x,&db_dot_x); CHKERRQ(ierr);
  double f_lambda;
  f_lambda = ctx->scale_lambda*(ctx->arc_ptr->diag*lambda + db_dot_x);
  ierr = ctx->set_lambda(f,&f_lambda,SCATTER_REVERSE); CHKERRQ(ierr);
  ierr = VecAXPY(f,-lambda,ctx->arc_ptr->F_lambda); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

struct PCShellCtx {
  PC pc;
  Mat ShellAij,Aij;
  ArcLenghtCtx* arc_ptr;
  PCShellCtx(Mat _ShellAij,Mat _Aij,ArcLenghtCtx* _arc_ptr): 
    ShellAij(_ShellAij),Aij(_Aij),arc_ptr(_arc_ptr) {
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
  PetscErrorCode ierr;
  void *void_ctx;
  ierr = PCShellGetContext(pc,&void_ctx); CHKERRQ(ierr);
  PCShellCtx *PCCtx = (PCShellCtx*)void_ctx;
  void *void_MatCtx;
  MatShellGetContext(PCCtx->ShellAij,&void_MatCtx);
  MatShellCtx *MatCtx = (MatShellCtx*)void_MatCtx;
  ierr = PCApply(PCCtx->pc,pc_f,pc_x); CHKERRQ(ierr);
  ierr = PCApply(PCCtx->pc,PCCtx->arc_ptr->F_lambda,PCCtx->arc_ptr->x_lambda); CHKERRQ(ierr);
  double db_dot_pc_x,db_dot_x_lambda;
  ierr = VecDot(PCCtx->arc_ptr->db,pc_x,&db_dot_pc_x); CHKERRQ(ierr);
  ierr = VecDot(PCCtx->arc_ptr->db,PCCtx->arc_ptr->x_lambda,&db_dot_x_lambda); CHKERRQ(ierr);
  double denominator = MatCtx->scale_lambda*PCCtx->arc_ptr->diag+MatCtx->scale_lambda*db_dot_x_lambda;
  double res_lambda;
  ierr = MatCtx->set_lambda(pc_f,&res_lambda,SCATTER_FORWARD); CHKERRQ(ierr);
  double ddlambda = (res_lambda - MatCtx->scale_lambda*db_dot_pc_x)/denominator;
  //cerr << res_lambda << " " << ddlambda << " " << db_dot_pc_x << " " << db_dot_x_lambda << " " << PCCtx->arc_ptr->diag << endl;
  if(ddlambda != ddlambda) SETERRQ(PETSC_COMM_SELF,1,"problem with constraint");
  ierr = VecAXPY(pc_x,ddlambda,PCCtx->arc_ptr->x_lambda); CHKERRQ(ierr);
  ierr = MatCtx->set_lambda(pc_x,&ddlambda,SCATTER_REVERSE); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode pc_setup_arc_length(PC pc) {
  PetscFunctionBegin;
  PetscErrorCode ierr;
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

}

#endif // __ARCLEGHTTOOLS_HPP__


