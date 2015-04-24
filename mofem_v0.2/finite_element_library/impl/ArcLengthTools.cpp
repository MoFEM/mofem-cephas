/** \file ArcLengthTools.hpp
 *
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

#include <ArcLengthTools.hpp>

//static ErrorCode rval;
static PetscErrorCode ierr;

namespace MoFEM {

PetscErrorCode ArcLengthCtx::setS(double s) { 
  PetscFunctionBegin;
  this->s = s;
  ierr = PetscPrintf(mField.get_comm(),"\tSet s = %6.4e\n",this->s); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode ArcLengthCtx::setAlphaBeta(double alpha,double beta) { 
  PetscFunctionBegin;
  this->alpha = alpha;
  this->beta = beta;
  ierr = PetscPrintf(mField.get_comm(),"\tSet alpha = %6.4e beta = %6.4e\n",this->alpha,this->beta); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

ArcLengthCtx::ArcLengthCtx(FieldInterface &m_field,const string &problem_name):
  mField(m_field),dlambda(0),diag(0),dx2(0),F_lambda2(0),res_lambda(0) {
  ierr = m_field.VecCreateGhost(problem_name,ROW,&F_lambda); CHKERRABORT(PETSC_COMM_WORLD,ierr);
  ierr = VecSetOption(F_lambda,VEC_IGNORE_NEGATIVE_INDICES,PETSC_TRUE); CHKERRABORT(PETSC_COMM_WORLD,ierr);
  ierr = m_field.VecCreateGhost(problem_name,ROW,&db); CHKERRABORT(PETSC_COMM_WORLD,ierr);
  ierr = m_field.VecCreateGhost(problem_name,ROW,&x_lambda); CHKERRABORT(PETSC_COMM_WORLD,ierr);
  ierr = m_field.VecCreateGhost(problem_name,ROW,&x0); CHKERRABORT(PETSC_COMM_WORLD,ierr);
  ierr = m_field.VecCreateGhost(problem_name,ROW,&dx); CHKERRABORT(PETSC_COMM_WORLD,ierr);
  const MoFEMProblem *problem_ptr;
  ierr = m_field.get_problem(problem_name,&problem_ptr); CHKERRABORT(PETSC_COMM_WORLD,ierr);
  NumeredDofMoFEMEntity_multiIndex& dofs_ptr_no_const 
	    = const_cast<NumeredDofMoFEMEntity_multiIndex&>(problem_ptr->numered_dofs_rows);
  NumeredDofMoFEMEntity_multiIndex::index<FieldName_mi_tag>::type::iterator hi_dit;
  dIt = dofs_ptr_no_const.get<FieldName_mi_tag>().lower_bound("LAMBDA");
  hi_dit = dofs_ptr_no_const.get<FieldName_mi_tag>().upper_bound("LAMBDA");
  if(distance(dIt,hi_dit)!=1) {
    PetscTraceBackErrorHandler(
	PETSC_COMM_WORLD,
	__LINE__,PETSC_FUNCTION_NAME,__FILE__,
	MOFEM_DATA_INCONSISTENCT,PETSC_ERROR_INITIAL,"can not find unique LAMBDA (load factor)",PETSC_NULL);
    PetscMPIAbortErrorHandler(PETSC_COMM_WORLD,
	__LINE__,PETSC_FUNCTION_NAME,__FILE__,
	MOFEM_DATA_INCONSISTENCT,PETSC_ERROR_INITIAL,"can not find unique LAMBDA (load factor)",PETSC_NULL);
  }
}

ArcLengthCtx::~ArcLengthCtx() {
  ierr = VecDestroy(&F_lambda); CHKERRABORT(PETSC_COMM_WORLD,ierr);
  ierr = VecDestroy(&db); CHKERRABORT(PETSC_COMM_WORLD,ierr);
  ierr = VecDestroy(&x_lambda); CHKERRABORT(PETSC_COMM_WORLD,ierr);
  ierr = VecDestroy(&x0); CHKERRABORT(PETSC_COMM_WORLD,ierr);
  ierr = VecDestroy(&dx); CHKERRABORT(PETSC_COMM_WORLD,ierr);
}

ArcLengthMatShell::ArcLengthMatShell(Mat aij,ArcLengthCtx *arc_ptr,string problem_name): 
    Aij(aij),arcPtr(arc_ptr),problemName(problem_name) {}

PetscErrorCode ArcLengthMatShell::set_lambda(Vec ksp_x,double *lambda,ScatterMode scattermode) {
  PetscFunctionBegin;
  const MoFEMProblem *problem_ptr;
  ierr = arcPtr->mField.get_problem(problemName,&problem_ptr); CHKERRQ(ierr);
  if(arcPtr->getPetscLocalDofIdx()!=-1) {
    PetscScalar *array;
    ierr = VecGetArray(ksp_x,&array); CHKERRQ(ierr);
    switch(scattermode) {
	case SCATTER_FORWARD:
	  *lambda = array[arcPtr->getPetscLocalDofIdx()];
	  break;
	case SCATTER_REVERSE:
	  array[arcPtr->getPetscLocalDofIdx()] = *lambda;
	  break;
	default:
	  SETERRQ(PETSC_COMM_SELF,1,"not implemented");
    }
    ierr = VecRestoreArray(ksp_x,&array); CHKERRQ(ierr);
  } 
  unsigned int part = arcPtr->getPart();
  ParallelComm* pcomm = ParallelComm::get_pcomm(&arcPtr->mField.get_moab(),MYPCOMM_INDEX);
  Vec lambda_ghost;
  if(pcomm->rank()==part) {
    ierr = VecCreateGhostWithArray(arcPtr->mField.get_comm(),1,1,0,PETSC_NULL,lambda,&lambda_ghost); CHKERRQ(ierr);
  } else {
    int one[] = {0};
    ierr = VecCreateGhostWithArray(arcPtr->mField.get_comm(),0,1,1,one,lambda,&lambda_ghost); CHKERRQ(ierr);
  }
  ierr = VecGhostUpdateBegin(lambda_ghost,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = VecGhostUpdateEnd(lambda_ghost,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = VecDestroy(&lambda_ghost); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

ArcLengthMatShell::~ArcLengthMatShell() {}

PCArcLengthCtx::PCArcLengthCtx(Mat shell_Aij,Mat _Aij,ArcLengthCtx* arc_ptr): 
  shellAij(shell_Aij),Aij(_Aij),arcPtr(arc_ptr) {
  ierr = PCCreate(arc_ptr->mField.get_comm(),&pC); CHKERRABORT(PETSC_COMM_WORLD,ierr);
}

PCArcLengthCtx::~PCArcLengthCtx() {
  ierr = PCDestroy(&pC); CHKERRABORT(PETSC_COMM_WORLD,ierr);
}

PrePostProcessForArcLength::PrePostProcessForArcLength(ArcLengthCtx *arcPtr):
  arcPtr(arcPtr) {}

PetscErrorCode PrePostProcessForArcLength::preProcess() {
  PetscFunctionBegin;
  switch(snes_ctx) {
    case CTX_SNESNONE:
    case CTX_SNESSETFUNCTION: {
      //F_lambda
      ierr = VecZeroEntries(arcPtr->F_lambda); CHKERRQ(ierr);
      ierr = VecGhostUpdateBegin(arcPtr->F_lambda,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
      ierr = VecGhostUpdateEnd(arcPtr->F_lambda,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
    }
    break;
    case CTX_SNESSETJACOBIAN: {
    }
    break;
    default:
      SETERRQ(PETSC_COMM_SELF,1,"not implemented");
  }
  
  PetscFunctionReturn(0);
}

PetscErrorCode PrePostProcessForArcLength::postProcess() {
  PetscFunctionBegin;
  switch(snes_ctx) {
    case CTX_SNESNONE: {
    }
    case CTX_SNESSETFUNCTION: {
      //F_lambda
      ierr = VecGhostUpdateBegin(arcPtr->F_lambda,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
      ierr = VecGhostUpdateEnd(arcPtr->F_lambda,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
      ierr = VecAssemblyBegin(arcPtr->F_lambda); CHKERRQ(ierr);
      ierr = VecAssemblyEnd(arcPtr->F_lambda); CHKERRQ(ierr);
      //F_lambda2
      ierr = VecDot(arcPtr->F_lambda,arcPtr->F_lambda,&arcPtr->F_lambda2); CHKERRQ(ierr);
      PetscPrintf(arcPtr->mField.get_comm(),"\tFlambda2 = %6.4e\n",arcPtr->F_lambda2);
    }
    break;
    case CTX_SNESSETJACOBIAN: {
    }
    break;
    default:
      SETERRQ(PETSC_COMM_SELF,1,"not implemented");
  }
  
  PetscFunctionReturn(0);
}

PetscErrorCode ArcLengthMatMultShellOp(Mat A,Vec x,Vec f) {
  PetscFunctionBegin;
  void *void_ctx;
  ierr = MatShellGetContext(A,&void_ctx); CHKERRQ(ierr);
  ArcLengthMatShell *ctx = (ArcLengthMatShell*)void_ctx;
  ierr = MatMult(ctx->Aij,x,f); CHKERRQ(ierr);
  double lambda;
  ierr = ctx->set_lambda(x,&lambda,SCATTER_FORWARD); CHKERRQ(ierr);
  double db_dot_x;
  ierr = VecDot(ctx->arcPtr->db,x,&db_dot_x); CHKERRQ(ierr);
  double f_lambda;
  f_lambda = ctx->arcPtr->diag*lambda + db_dot_x;
  ierr = ctx->set_lambda(f,&f_lambda,SCATTER_REVERSE); CHKERRQ(ierr);
  ierr = VecAXPY(f,lambda,ctx->arcPtr->F_lambda); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode PCApplyArcLength(PC pc,Vec pc_f,Vec pc_x) {
  PetscFunctionBegin;
  void *void_ctx;
  ierr = PCShellGetContext(pc,&void_ctx); CHKERRQ(ierr);
  PCArcLengthCtx *ctx = (PCArcLengthCtx*)void_ctx;
  void *void_MatCtx;
  MatShellGetContext(ctx->shellAij,&void_MatCtx);
  ArcLengthMatShell *MatCtx = (ArcLengthMatShell*)void_MatCtx;
  ierr = PCApply(ctx->pC,pc_f,pc_x); CHKERRQ(ierr);
  ierr = PCApply(ctx->pC,ctx->arcPtr->F_lambda,ctx->arcPtr->x_lambda); CHKERRQ(ierr);
  double db_dot_pc_x,db_dot_x_lambda;
  ierr = VecDot(ctx->arcPtr->db,pc_x,&db_dot_pc_x); CHKERRQ(ierr);
  ierr = VecDot(ctx->arcPtr->db,ctx->arcPtr->x_lambda,&db_dot_x_lambda); CHKERRQ(ierr);
  double denominator = ctx->arcPtr->diag+db_dot_x_lambda;
  double res_lambda;
  ierr = MatCtx->set_lambda(pc_f,&res_lambda,SCATTER_FORWARD); CHKERRQ(ierr);
  double ddlambda = (res_lambda - db_dot_pc_x)/denominator;
  if(ddlambda != ddlambda) {
    ostringstream ss;
    ss << "problem with ddlambda: " << res_lambda << " " << ddlambda << " " << db_dot_pc_x << " " << db_dot_x_lambda << " " << ctx->arcPtr->diag;
    //cerr << ss.str() << endl;
    SETERRQ(PETSC_COMM_SELF,1,ss.str().c_str());
  }
  ierr = VecAXPY(pc_x,ddlambda,ctx->arcPtr->x_lambda); CHKERRQ(ierr);
  ierr = MatCtx->set_lambda(pc_x,&ddlambda,SCATTER_REVERSE); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode PCSetupArcLength(PC pc) {
  PetscFunctionBegin;
  void *void_ctx;
  ierr = PCShellGetContext(pc,&void_ctx); CHKERRQ(ierr);
  PCArcLengthCtx *ctx = (PCArcLengthCtx*)void_ctx;
  ierr = PCSetFromOptions(ctx->pC); CHKERRQ(ierr);
  ierr = PCGetOperators(pc,&ctx->shellAij,&ctx->Aij); CHKERRQ(ierr);
  ierr = PCSetOperators(ctx->pC,ctx->shellAij,ctx->Aij); CHKERRQ(ierr);
  ierr = PCSetUp(ctx->pC); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

SphericalArcLengthControl::SphericalArcLengthControl(ArcLengthCtx *arc_ptr): 
  FEMethod(),arcPtr(arc_ptr) {
  PetscInt ghosts[1] = { 0 };
  ParallelComm* pcomm = ParallelComm::get_pcomm(&arc_ptr->mField.get_moab(),MYPCOMM_INDEX);
  if(pcomm->rank() == 0) {
    VecCreateGhost(arc_ptr->mField.get_comm(),1,1,0,ghosts,&ghostDiag);
  } else {
    VecCreateGhost(arc_ptr->mField.get_comm(),0,1,1,ghosts,&ghostDiag);
  }
}

SphericalArcLengthControl::~SphericalArcLengthControl() {
  VecDestroy(&ghostDiag);
}

PetscErrorCode SphericalArcLengthControl::preProcess() {
  PetscFunctionBegin;
  switch (ts_ctx) {
    case CTX_TSSETIFUNCTION: {
      snes_ctx = CTX_SNESSETFUNCTION;
      snes_f = ts_F;
      break;
    }
    case CTX_TSSETIJACOBIAN: {
	snes_ctx = CTX_SNESSETJACOBIAN;
	snes_B = ts_B;
	break;
    }
    default:
    break;
  }
  switch(snes_ctx) {
    case CTX_SNESSETFUNCTION: { 
	ierr = calculateDxAndDlambda(snes_x); CHKERRQ(ierr);
	ierr = calculateDb(); CHKERRQ(ierr);
    }
    break;
    case CTX_SNESSETJACOBIAN: {
    }
    break;
    default:
    break;
  }
  PetscFunctionReturn(0);
}

double SphericalArcLengthControl::calculateLambdaInt() {
  PetscFunctionBegin;
  return arcPtr->alpha*arcPtr->dx2 + pow(arcPtr->dlambda,2)*pow(arcPtr->beta,2)*arcPtr->F_lambda2;
  PetscFunctionReturn(0);
}

PetscErrorCode SphericalArcLengthControl::calculateDb() {
  PetscFunctionBegin;
  //db
  ierr = VecCopy(arcPtr->dx,arcPtr->db); CHKERRQ(ierr);
  ierr = VecScale(arcPtr->db,2*arcPtr->alpha); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode SphericalArcLengthControl::operator()() {
  PetscFunctionBegin;
  switch(snes_ctx) {
    case CTX_SNESSETFUNCTION: {
      arcPtr->res_lambda = calculateLambdaInt() - pow(arcPtr->s,2);
      ierr = VecSetValue(snes_f,arcPtr->getPetscGloablDofIdx(),arcPtr->res_lambda,ADD_VALUES); CHKERRQ(ierr);
      PetscPrintf(arcPtr->mField.get_comm(),"\tres_lambda = %6.4e\n",arcPtr->res_lambda);
    }
    break; 
    case CTX_SNESSETJACOBIAN: {
      double diag = 2*arcPtr->dlambda*pow(arcPtr->beta,2)*arcPtr->F_lambda2;
      ierr = VecSetValue(ghostDiag,0,diag,INSERT_VALUES); CHKERRQ(ierr);
      ierr = MatSetValue(snes_B,arcPtr->getPetscGloablDofIdx(),arcPtr->getPetscGloablDofIdx(),1,ADD_VALUES); CHKERRQ(ierr);
    }
    break;
    default:
    break;
  }	
  
  PetscFunctionReturn(0);
}

PetscErrorCode SphericalArcLengthControl::postProcess() {
  PetscFunctionBegin;
  switch(snes_ctx) {
    case CTX_SNESSETFUNCTION: { 
      PetscPrintf(arcPtr->mField.get_comm(),"\tlambda = %6.4e\n",arcPtr->getFieldData());  
    }
    break;
    case CTX_SNESSETJACOBIAN: {
      ierr = MatAssemblyBegin(snes_B,MAT_FLUSH_ASSEMBLY); CHKERRQ(ierr);
      ierr = MatAssemblyEnd(snes_B,MAT_FLUSH_ASSEMBLY); CHKERRQ(ierr);
      ierr = VecAssemblyBegin(ghostDiag); CHKERRQ(ierr);
      ierr = VecAssemblyEnd(ghostDiag); CHKERRQ(ierr);
      ierr = VecGhostUpdateBegin(ghostDiag,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
      ierr = VecGhostUpdateEnd(ghostDiag,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
      double *diag;
      ierr = VecGetArray(ghostDiag,&diag); CHKERRQ(ierr);
      arcPtr->diag = *diag;
      ierr = VecRestoreArray(ghostDiag,&diag); CHKERRQ(ierr);
      PetscPrintf(arcPtr->mField.get_comm(),"\tdiag = %6.4e\n",arcPtr->diag);
    }
    break;
    default:
    break;
  }
  PetscFunctionReturn(0);
}

PetscErrorCode SphericalArcLengthControl::calculateDxAndDlambda(Vec x) {
  PetscFunctionBegin;
  //dx
  ierr = VecCopy(x,arcPtr->dx); CHKERRQ(ierr);
  ierr = VecAXPY(arcPtr->dx,-1,arcPtr->x0); CHKERRQ(ierr);
  //dlambda
  if(arcPtr->getPetscLocalDofIdx()!=-1) {
    double *array;
    ierr = VecGetArray(arcPtr->dx,&array); CHKERRQ(ierr);
    arcPtr->dlambda = array[arcPtr->getPetscLocalDofIdx()];
    array[arcPtr->getPetscLocalDofIdx()] = 0;
    ierr = VecRestoreArray(arcPtr->dx,&array); CHKERRQ(ierr);
  }
  int part = arcPtr->getPart();
  MPI_Bcast(&(arcPtr->dlambda),1,MPI_DOUBLE,part,arcPtr->mField.get_comm());
  //dx2
  ierr = VecDot(arcPtr->dx,arcPtr->dx,&arcPtr->dx2); CHKERRQ(ierr);
  PetscPrintf(arcPtr->mField.get_comm(),"\tdlambda = %6.4e dx2 = %6.4e\n",arcPtr->dlambda,arcPtr->dx2);
  PetscFunctionReturn(0);
}

PetscErrorCode SphericalArcLengthControl::calculateInitDlambda(double *dlambda) {
  PetscFunctionBegin;
  *dlambda = sqrt(pow(arcPtr->s,2)/(pow(arcPtr->beta,2)*arcPtr->F_lambda2));
  if(!(*dlambda == *dlambda)) {
    ostringstream sss;
    sss << "s " << arcPtr->s << " " << arcPtr->beta << " " << arcPtr->F_lambda2;
    SETERRQ(PETSC_COMM_SELF,1,sss.str().c_str());
  }
  PetscFunctionReturn(0);
}

PetscErrorCode SphericalArcLengthControl::setDlambdaToX(Vec x,double dlambda) {
  PetscFunctionBegin;
  //check if locl dof idx is non zero, i.e. that lambda is acessible from this processor
  if(arcPtr->getPetscLocalDofIdx()!=-1) {
    double *array;
    ierr = VecGetArray(x,&array); CHKERRQ(ierr);
    double lambda_old = array[arcPtr->getPetscLocalDofIdx()];
    if(!(dlambda == dlambda)) {
	ostringstream sss;
	sss << "s " << arcPtr->s << " " << arcPtr->beta << " " << arcPtr->F_lambda2;
	SETERRQ(PETSC_COMM_SELF,1,sss.str().c_str());
    }
    array[arcPtr->getPetscLocalDofIdx()] = lambda_old + dlambda;
    PetscPrintf(arcPtr->mField.get_comm(),"\tlambda = %6.4e, %6.4e (%6.4e)\n",
	lambda_old, array[arcPtr->getPetscLocalDofIdx()], dlambda);
    ierr = VecRestoreArray(x,&array); CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}



}
