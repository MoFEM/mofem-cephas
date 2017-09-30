/** \file ArcLengthTools.hpp

 Implementation of Arc Length element

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

// ********************
// Arc-length ctx class

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

ArcLengthCtx::ArcLengthCtx(
  MoFEM::Interface& m_field,const std::string& problem_name,const std::string& field_name
):
mField(m_field),
dx2(0),
F_lambda2(0),
res_lambda(0) {
  ierr = m_field.query_interface<VecManager>()->vecCreateGhost(problem_name,ROW,&F_lambda); CHKERRABORT(PETSC_COMM_WORLD,ierr);
  ierr = VecSetOption(F_lambda,VEC_IGNORE_NEGATIVE_INDICES,PETSC_TRUE); CHKERRABORT(PETSC_COMM_WORLD,ierr);
  ierr = VecDuplicate(F_lambda,&db); CHKERRABORT(PETSC_COMM_WORLD,ierr);
  ierr = VecDuplicate(F_lambda,&xLambda); CHKERRABORT(PETSC_COMM_WORLD,ierr);
  ierr = VecDuplicate(F_lambda,&x0); CHKERRABORT(PETSC_COMM_WORLD,ierr);
  ierr = VecDuplicate(F_lambda,&dx); CHKERRABORT(PETSC_COMM_WORLD,ierr);

  const Problem *problem_ptr;
  ierr = m_field.get_problem(problem_name,&problem_ptr); CHKERRABORT(PETSC_COMM_WORLD,ierr);
  boost::shared_ptr<NumeredDofEntity_multiIndex> dofs_ptr_no_const = problem_ptr->getNumeredDofsRows();
  NumeredDofEntityByFieldName::iterator hi_dit;
  dIt = dofs_ptr_no_const->get<FieldName_mi_tag>().lower_bound(field_name);
  hi_dit = dofs_ptr_no_const->get<FieldName_mi_tag>().upper_bound(field_name);

  if(distance(dIt,hi_dit)!=1) {
    PetscTraceBackErrorHandler(
      PETSC_COMM_WORLD,
      __LINE__,PETSC_FUNCTION_NAME,__FILE__,
      MOFEM_DATA_INCONSISTENCY,PETSC_ERROR_INITIAL,
      "can not find unique LAMBDA (load factor)",PETSC_NULL
    );
    PetscMPIAbortErrorHandler(PETSC_COMM_WORLD,
      __LINE__,PETSC_FUNCTION_NAME,__FILE__,
      MOFEM_DATA_INCONSISTENCY,PETSC_ERROR_INITIAL,
      "can not find unique LAMBDA (load factor)",PETSC_NULL
    );
  }

  // ParallelComm* pcomm = ParallelComm::get_pcomm(&mField.get_moab(),MYPCOMM_INDEX);
  // BARRIER_RANK_START(pcomm)
  // cerr << **dIt << endl;
  // BARRIER_RANK_END(pcomm)

  if((unsigned int)mField.get_comm_rank()==(*dIt)->getPart()) {
    ierr = VecCreateGhostWithArray(
      mField.get_comm(),1,1,0,PETSC_NULL,&dLambda,&ghosTdLambda
    ); CHKERRABORT(PETSC_COMM_WORLD,ierr);
    ierr = VecCreateGhostWithArray(
      mField.get_comm(),1,1,0,PETSC_NULL,&dIag,&ghostDiag
    ); CHKERRABORT(PETSC_COMM_WORLD,ierr);
  } else {
    int one[] = {0};
    ierr = VecCreateGhostWithArray(
      mField.get_comm(),0,1,1,one,&dLambda,&ghosTdLambda
    ); CHKERRABORT(PETSC_COMM_WORLD,ierr);
    ierr = VecCreateGhostWithArray(
      mField.get_comm(),0,1,1,one,&dIag,&ghostDiag
    ); CHKERRABORT(PETSC_COMM_WORLD,ierr);
  }
  dLambda = 0;
  dIag = 0;
}

ArcLengthCtx::~ArcLengthCtx() {
  ierr = VecDestroy(&F_lambda); CHKERRABORT(PETSC_COMM_WORLD,ierr);
  ierr = VecDestroy(&db); CHKERRABORT(PETSC_COMM_WORLD,ierr);
  ierr = VecDestroy(&xLambda); CHKERRABORT(PETSC_COMM_WORLD,ierr);
  ierr = VecDestroy(&x0); CHKERRABORT(PETSC_COMM_WORLD,ierr);
  ierr = VecDestroy(&dx); CHKERRABORT(PETSC_COMM_WORLD,ierr);
  ierr = VecDestroy(&ghosTdLambda); CHKERRABORT(PETSC_COMM_WORLD,ierr);
  ierr = VecDestroy(&ghostDiag); CHKERRABORT(PETSC_COMM_WORLD,ierr);
}

// ***********************
// Arc-length shell matrix

ArcLengthMatShell::ArcLengthMatShell(Mat aij,ArcLengthCtx *arc_ptr_raw,string problem_name):
Aij(aij),arcPtrRaw(arc_ptr_raw),problemName(problem_name) {
  ierr = PetscObjectReference((PetscObject)aij); CHKERRABORT(PETSC_COMM_WORLD,ierr);
}

ArcLengthMatShell::ArcLengthMatShell(Mat aij,boost::shared_ptr<ArcLengthCtx> arc_ptr,string problem_name):
Aij(aij),arcPtrRaw(arc_ptr.get()),problemName(problem_name),arcPtr(arc_ptr) {
  ierr = PetscObjectReference((PetscObject)aij); CHKERRABORT(PETSC_COMM_WORLD,ierr);
}

ArcLengthMatShell::~ArcLengthMatShell() {
  ierr = MatDestroy(&Aij); CHKERRABORT(PETSC_COMM_WORLD,ierr);
}

PetscErrorCode ArcLengthMatShell::setLambda(Vec ksp_x,double *lambda,ScatterMode scattermode) {
  PetscFunctionBegin;

  int part = arcPtrRaw->getPart();
  int rank = arcPtrRaw->mField.get_comm_rank();

  Vec lambda_ghost;
  if(rank==part) {
    ierr = VecCreateGhostWithArray(
      arcPtrRaw->mField.get_comm(),1,1,0,PETSC_NULL,lambda,&lambda_ghost
    ); CHKERRQ(ierr);
  } else {
    int one[] = {0};
    ierr = VecCreateGhostWithArray(
      arcPtrRaw->mField.get_comm(),0,1,1,one,lambda,&lambda_ghost
    ); CHKERRQ(ierr);
  }

  switch(scattermode) {
    case SCATTER_FORWARD: {
      int idx = arcPtrRaw->getPetscGlobalDofIdx();
      if(part == rank) {
        ierr = VecGetValues(ksp_x,1,&idx,lambda); CHKERRQ(ierr);
      }
      ierr = VecGhostUpdateBegin(lambda_ghost,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
      ierr = VecGhostUpdateEnd(lambda_ghost,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
    }
    break;
    case SCATTER_REVERSE: {
      // ierr = VecGhostUpdateBegin(lambda_ghost,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
      // ierr = VecGhostUpdateEnd(lambda_ghost,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
      if(part == rank) {
        PetscScalar *array;
        ierr = VecGetArray(ksp_x,&array); CHKERRQ(ierr);
        array[arcPtrRaw->getPetscLocalDofIdx()] = *lambda;
        ierr = VecRestoreArray(ksp_x,&array); CHKERRQ(ierr);
      }
    }
    break;
    default:
    SETERRQ(PETSC_COMM_SELF,MOFEM_NOT_IMPLEMENTED,"not implemented");
  }

  ierr = VecDestroy(&lambda_ghost); CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

PetscErrorCode ArcLengthMatMultShellOp(Mat A,Vec x,Vec f) {
  PetscFunctionBegin;
  void *void_ctx;
  ierr = MatShellGetContext(A,&void_ctx); CHKERRQ(ierr);
  ArcLengthMatShell *ctx = (ArcLengthMatShell*)void_ctx;
  ierr = MatMult(ctx->Aij,x,f); CHKERRQ(ierr);
  double lambda;
  ierr = ctx->setLambda(x,&lambda,SCATTER_FORWARD); CHKERRQ(ierr);
  double db_dot_x;
  ierr = VecDot(ctx->arcPtrRaw->db,x,&db_dot_x); CHKERRQ(ierr);
  double f_lambda;
  f_lambda = ctx->arcPtrRaw->dIag*lambda + db_dot_x;
  ierr = ctx->setLambda(f,&f_lambda,SCATTER_REVERSE); CHKERRQ(ierr);
  ierr = VecAXPY(f,lambda,ctx->arcPtrRaw->F_lambda); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

// arc-length preconditioner

PCArcLengthCtx::PCArcLengthCtx(Mat shell_Aij,Mat aij,ArcLengthCtx* arc_ptr):
  shellAij(shell_Aij),Aij(aij),arcPtrRaw(arc_ptr) {
  ierr = PCCreate(arc_ptr->mField.get_comm(),&pC); CHKERRABORT(PETSC_COMM_WORLD,ierr);
}

PCArcLengthCtx::PCArcLengthCtx(Mat shell_Aij,Mat aij,boost::shared_ptr<ArcLengthCtx>& arc_ptr):
  shellAij(shell_Aij),Aij(aij),arcPtrRaw(arc_ptr.get()),arcPtr(arc_ptr) {
  ierr = PCCreate(arc_ptr->mField.get_comm(),&pC); CHKERRABORT(PETSC_COMM_WORLD,ierr);
}

PCArcLengthCtx::PCArcLengthCtx(PC pc,Mat shell_Aij,Mat aij,boost::shared_ptr<ArcLengthCtx>& arc_ptr):
  pC(pc),shellAij(shell_Aij),Aij(aij),arcPtrRaw(arc_ptr.get()),arcPtr(arc_ptr) {
  ierr = PetscObjectReference((PetscObject)pC); CHKERRABORT(PETSC_COMM_WORLD,ierr);
}

PCArcLengthCtx::~PCArcLengthCtx() {
  ierr = PCDestroy(&pC); CHKERRABORT(PETSC_COMM_WORLD,ierr);
}

PetscErrorCode PCApplyArcLength(PC pc,Vec pc_f,Vec pc_x) {
  PetscFunctionBegin;
  void *void_ctx;
  ierr = PCShellGetContext(pc,&void_ctx); CHKERRQ(ierr);
  PCArcLengthCtx *ctx = (PCArcLengthCtx*)void_ctx;
  void *void_MatCtx;
  MatShellGetContext(ctx->shellAij,&void_MatCtx);
  ArcLengthMatShell *mat_ctx = (ArcLengthMatShell*)void_MatCtx;
  ierr = PCApply(ctx->pC,pc_f,pc_x); CHKERRQ(ierr);
  ierr = PCApply(ctx->pC,ctx->arcPtrRaw->F_lambda,ctx->arcPtrRaw->xLambda); CHKERRQ(ierr);
  double db_dot_pc_x,db_dot_x_lambda;
  ierr = VecDot(ctx->arcPtrRaw->db,pc_x,&db_dot_pc_x); CHKERRQ(ierr);
  ierr = VecDot(ctx->arcPtrRaw->db,ctx->arcPtrRaw->xLambda,&db_dot_x_lambda); CHKERRQ(ierr);
  double denominator = ctx->arcPtrRaw->dIag+db_dot_x_lambda;
  double res_lambda;
  ierr = mat_ctx->setLambda(pc_f,&res_lambda,SCATTER_FORWARD); CHKERRQ(ierr);
  double ddlambda = (res_lambda - db_dot_pc_x)/denominator;
  // cerr << denominator << " " << res_lambda << " " << ddlambda << endl;
  if(ddlambda != ddlambda || denominator == 0) {
    double nrm2_pc_f,nrm2_db,nrm2_pc_x,nrm2_xLambda;
    ierr = VecNorm(pc_f,NORM_2,&nrm2_pc_f); CHKERRQ(ierr);
    ierr = VecNorm(ctx->arcPtrRaw->db,NORM_2,&nrm2_db); CHKERRQ(ierr);
    ierr = VecNorm(pc_x,NORM_2,&nrm2_pc_x); CHKERRQ(ierr);
    ierr = VecNorm(ctx->arcPtrRaw->xLambda,NORM_2,&nrm2_xLambda); CHKERRQ(ierr);
    std::ostringstream ss;
    ss
    << "problem with ddlambda=" << res_lambda
    << " res_lambda=" << res_lambda
    << " denominator=" << denominator
    << " ddlamnda=" << ddlambda
    << " db_dot_pc_x=" << db_dot_pc_x
    << " db_dot_x_lambda=" << db_dot_x_lambda
    << " diag=" << ctx->arcPtrRaw->dIag
    << " nrm2_db=" << nrm2_db
    << " nrm2_pc_f=" << nrm2_pc_f
    << " nrm2_pc_x=" << nrm2_pc_x
    << " nrm2_xLambda=" << nrm2_xLambda;
    SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,ss.str().c_str());
  }
  ierr = VecAXPY(pc_x,ddlambda,ctx->arcPtrRaw->xLambda); CHKERRQ(ierr);
  ierr = mat_ctx->setLambda(pc_x,&ddlambda,SCATTER_REVERSE); CHKERRQ(ierr);
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

// ***********************
// Zero F_lambda vector

ZeroFLmabda::ZeroFLmabda(
  boost::shared_ptr<ArcLengthCtx> arc_ptr
): arcPtr(arc_ptr) {
}

PetscErrorCode ZeroFLmabda::preProcess() {
  PetscFunctionBegin;
  switch(snes_ctx) {
    case CTX_SNESSETFUNCTION: {
      ierr = VecZeroEntries(arcPtr->F_lambda); CHKERRQ(ierr);
      ierr = VecGhostUpdateBegin(arcPtr->F_lambda,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
      ierr = VecGhostUpdateEnd(arcPtr->F_lambda,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
    }
    break;
    default:
    SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"Impossible case");
  }
  PetscFunctionReturn(0);
}

AssembleFLmabda::AssembleFLmabda(
  boost::shared_ptr<ArcLengthCtx> arc_ptr,
  boost::shared_ptr<DirichletDisplacementBc> bc
): arcPtr(arc_ptr),bC(bc) {
}

PetscErrorCode AssembleFLmabda::preProcess() {
  PetscFunctionBegin;
  PetscFunctionReturn(0);
}
PetscErrorCode AssembleFLmabda::operator()() {
  PetscFunctionBegin;
  PetscFunctionReturn(0);
}
PetscErrorCode AssembleFLmabda::postProcess() {
  PetscFunctionBegin;
  switch(snes_ctx) {
    case CTX_SNESSETFUNCTION: {
      //F_lambda
      ierr = VecAssemblyBegin(arcPtr->F_lambda); CHKERRQ(ierr);
      ierr = VecAssemblyEnd(arcPtr->F_lambda); CHKERRQ(ierr);
      ierr = VecGhostUpdateBegin(arcPtr->F_lambda,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
      ierr = VecGhostUpdateEnd(arcPtr->F_lambda,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
      if(bC) {
        for(
          std::vector<int>::iterator vit = bC->dofsIndices.begin();
          vit!=bC->dofsIndices.end();vit++
        ) {
          ierr = VecSetValue(arcPtr->F_lambda,*vit,0,INSERT_VALUES); CHKERRQ(ierr);
        }
        ierr = VecAssemblyBegin(arcPtr->F_lambda); CHKERRQ(ierr);
        ierr = VecAssemblyEnd(arcPtr->F_lambda); CHKERRQ(ierr);
      }
      ierr = VecDot(arcPtr->F_lambda,arcPtr->F_lambda,&arcPtr->F_lambda2); CHKERRQ(ierr);
      PetscPrintf(PETSC_COMM_WORLD,"\tF_lambda2 = %6.4e\n",arcPtr->F_lambda2);
      //add F_lambda
      ierr = VecAssemblyBegin(snes_f); CHKERRQ(ierr);
      ierr = VecAssemblyEnd(snes_f); CHKERRQ(ierr);
      double lambda = arcPtr->getFieldData();
      ierr = VecAXPY(snes_f,lambda,arcPtr->F_lambda); CHKERRQ(ierr);
      double fnorm;
      ierr = VecNorm(snes_f,NORM_2,&fnorm); CHKERRQ(ierr);
      PetscPrintf(PETSC_COMM_WORLD,"\tfnorm = %6.4e lambda = %6.4g\n",fnorm,lambda);
    }
    break;
    default:
    SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"Impossible case");
  }
  PetscFunctionReturn(0);
}

// ************************
// Simple arc-length method

SimpleArcLengthControl::SimpleArcLengthControl(boost::shared_ptr<ArcLengthCtx>& arc_ptr,const bool assemble):
FEMethod(),
arcPtr(arc_ptr),
aSsemble(assemble) {
}

SimpleArcLengthControl::~SimpleArcLengthControl() {
}

PetscErrorCode SimpleArcLengthControl::preProcess() {
  PetscFunctionBegin;
  switch(snes_ctx) {
    case CTX_SNESSETFUNCTION: {
      if(aSsemble) {
        ierr = VecAssemblyBegin(snes_f); CHKERRQ(ierr);
        ierr = VecAssemblyEnd(snes_f); CHKERRQ(ierr);
      }
      ierr = calculateDxAndDlambda(snes_x); CHKERRQ(ierr);
      ierr = calculateDb(); CHKERRQ(ierr);
    }
    break;
    case CTX_SNESSETJACOBIAN: {
      if(aSsemble) {
        ierr = MatAssemblyBegin(snes_B,MAT_FLUSH_ASSEMBLY); CHKERRQ(ierr);
        ierr = MatAssemblyEnd(snes_B,MAT_FLUSH_ASSEMBLY); CHKERRQ(ierr);
      }
    }
    break;
    default:
    break;
  }
  PetscFunctionReturn(0);
}

PetscErrorCode SimpleArcLengthControl::operator()() {
  PetscFunctionBegin;
  switch(snes_ctx) {
    case CTX_SNESSETFUNCTION: {
      arcPtr->res_lambda = calculateLambdaInt() - arcPtr->s;
      ierr = VecSetValue(
        snes_f,arcPtr->getPetscGlobalDofIdx(),arcPtr->res_lambda,ADD_VALUES
      ); CHKERRQ(ierr);
    }
    break;
    case CTX_SNESSETJACOBIAN: {
      arcPtr->dIag = arcPtr->beta;
      ierr = MatSetValue(
        snes_B,arcPtr->getPetscGlobalDofIdx(),arcPtr->getPetscGlobalDofIdx(),1,ADD_VALUES
      ); CHKERRQ(ierr);
    }
    break;
    default:
    break;
  }
  PetscFunctionReturn(0);
}

PetscErrorCode SimpleArcLengthControl::postProcess() {
  PetscFunctionBegin;
  switch(snes_ctx) {
    case CTX_SNESSETFUNCTION: {
      if(aSsemble) {
        ierr = VecAssemblyBegin(snes_f); CHKERRQ(ierr);
        ierr = VecAssemblyEnd(snes_f); CHKERRQ(ierr);
      }
    }
    break;
    case CTX_SNESSETJACOBIAN: {
      if(aSsemble) {
        ierr = MatAssemblyBegin(snes_B,MAT_FLUSH_ASSEMBLY); CHKERRQ(ierr);
        ierr = MatAssemblyEnd(snes_B,MAT_FLUSH_ASSEMBLY); CHKERRQ(ierr);
      }
      ierr = VecGhostUpdateBegin(arcPtr->ghostDiag,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
      ierr = VecGhostUpdateEnd(arcPtr->ghostDiag,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
    }
    break;
    default:
    break;
  }
  PetscFunctionReturn(0);
}

double SimpleArcLengthControl::calculateLambdaInt() {
  return arcPtr->beta*arcPtr->dLambda;
}

PetscErrorCode SimpleArcLengthControl::calculateDb() {
  PetscFunctionBegin;
  ierr = VecZeroEntries(arcPtr->db); CHKERRQ(ierr);
  ierr = VecGhostUpdateBegin(arcPtr->db,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = VecGhostUpdateEnd(arcPtr->db,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode SimpleArcLengthControl::calculateDxAndDlambda(Vec x) {
  PetscFunctionBegin;
  // Calculate dx
  ierr = VecCopy(x,arcPtr->dx); CHKERRQ(ierr);
  ierr = VecAXPY(arcPtr->dx,-1,arcPtr->x0); CHKERRQ(ierr);
  ierr = VecGhostUpdateBegin(arcPtr->x0,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = VecGhostUpdateEnd(arcPtr->x0,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = VecGhostUpdateBegin(arcPtr->dx,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = VecGhostUpdateEnd(arcPtr->dx,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  // Calculate dlambda
  if(arcPtr->getPetscLocalDofIdx()!=-1) {
    double *array;
    ierr = VecGetArray(arcPtr->dx,&array); CHKERRQ(ierr);
    arcPtr->dLambda = array[arcPtr->getPetscLocalDofIdx()];
    array[arcPtr->getPetscLocalDofIdx()] = 0;
    ierr = VecRestoreArray(arcPtr->dx,&array); CHKERRQ(ierr);
  }
  ierr = VecGhostUpdateBegin(arcPtr->ghosTdLambda,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = VecGhostUpdateEnd(arcPtr->ghosTdLambda,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  // Calculate dx2
  ierr = VecDot(arcPtr->dx,arcPtr->dx,&arcPtr->dx2); CHKERRQ(ierr);
  PetscPrintf(PETSC_COMM_WORLD,"\tdx2 = %6.4e\n",arcPtr->dx2);
  PetscFunctionReturn(0);
}

// ***************************
// Spherical arc-length contri

SphericalArcLengthControl::SphericalArcLengthControl(ArcLengthCtx *arc_ptr_raw):
FEMethod(),
arcPtrRaw(arc_ptr_raw) {
}

SphericalArcLengthControl::SphericalArcLengthControl(boost::shared_ptr<ArcLengthCtx>& arc_ptr):
FEMethod(),
arcPtrRaw(arc_ptr.get()),
arcPtr(arc_ptr) {
}

SphericalArcLengthControl::~SphericalArcLengthControl() {
}

PetscErrorCode SphericalArcLengthControl::preProcess() {
  PetscFunctionBegin;
  switch (ts_ctx) {
    case CTX_TSSETIFUNCTION: {
      snes_ctx = CTX_SNESSETFUNCTION;
      snes_x = ts_u;
      snes_f = ts_F;
      break;
    }
    case CTX_TSSETIJACOBIAN: {
      snes_ctx = CTX_SNESSETJACOBIAN;
      snes_x = ts_u;
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
  return arcPtrRaw->alpha*arcPtrRaw->dx2 + pow(arcPtrRaw->dLambda,2)*pow(arcPtrRaw->beta,2)*arcPtrRaw->F_lambda2;
}

PetscErrorCode SphericalArcLengthControl::calculateDb() {
  PetscFunctionBegin;
  ierr = VecCopy(arcPtrRaw->dx,arcPtrRaw->db); CHKERRQ(ierr);
  ierr = VecScale(arcPtrRaw->db,2*arcPtrRaw->alpha); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode SphericalArcLengthControl::operator()() {
  PetscFunctionBegin;
  switch(snes_ctx) {
    case CTX_SNESSETFUNCTION: {
      arcPtrRaw->res_lambda = calculateLambdaInt() - pow(arcPtrRaw->s,2);
      ierr = VecSetValue(
        snes_f,arcPtrRaw->getPetscGlobalDofIdx(),arcPtrRaw->res_lambda,ADD_VALUES
      ); CHKERRQ(ierr);
      PetscPrintf(arcPtrRaw->mField.get_comm(),"\tres_lambda = %6.4e\n",arcPtrRaw->res_lambda);
    }
    break;
    case CTX_SNESSETJACOBIAN: {
      arcPtrRaw->dIag = 2*arcPtrRaw->dLambda*pow(arcPtrRaw->beta,2)*arcPtrRaw->F_lambda2;
      ierr = MatSetValue(
        snes_B,arcPtrRaw->getPetscGlobalDofIdx(),arcPtrRaw->getPetscGlobalDofIdx(),1,ADD_VALUES
      ); CHKERRQ(ierr);
    }
    break;
    default:
    break;
  }

  PetscFunctionReturn(0);
}

PetscErrorCode SphericalArcLengthControl::postProcess() {
  PetscFunctionBegin;
  switch (ts_ctx) {
    case CTX_TSSETIFUNCTION: {
      snes_ctx = CTX_SNESSETFUNCTION;
      snes_x = ts_u;
      snes_f = ts_F;
      break;
    }
    case CTX_TSSETIJACOBIAN: {
      snes_ctx = CTX_SNESSETJACOBIAN;
      snes_x = ts_u;
      snes_B = ts_B;
      break;
    }
    default:
    break;
  }
  switch(snes_ctx) {
    case CTX_SNESSETFUNCTION: {
      PetscPrintf(arcPtrRaw->mField.get_comm(),"\tlambda = %6.4e\n",arcPtrRaw->getFieldData());
    }
    break;
    case CTX_SNESSETJACOBIAN: {
      // VecView(arcPtrRaw->ghostDiag,PETSC_VIEWER_STDOUT_WORLD);
      ierr = VecGhostUpdateBegin(arcPtrRaw->ghostDiag,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
      ierr = VecGhostUpdateEnd(arcPtrRaw->ghostDiag,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
      ierr = MatAssemblyBegin(snes_B,MAT_FLUSH_ASSEMBLY); CHKERRQ(ierr);
      ierr = MatAssemblyEnd(snes_B,MAT_FLUSH_ASSEMBLY); CHKERRQ(ierr);
      PetscPrintf(arcPtrRaw->mField.get_comm(),"\tdiag = %6.4e\n",arcPtrRaw->dIag);
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
  ierr = VecCopy(x,arcPtrRaw->dx); CHKERRQ(ierr);
  ierr = VecAXPY(arcPtrRaw->dx,-1,arcPtrRaw->x0); CHKERRQ(ierr);
  ierr = VecGhostUpdateBegin(arcPtrRaw->dx,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = VecGhostUpdateEnd(arcPtrRaw->dx,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  //dlambda
  if(arcPtrRaw->getPetscLocalDofIdx()!=-1) {
    double *array;
    ierr = VecGetArray(arcPtrRaw->dx,&array); CHKERRQ(ierr);
    arcPtrRaw->dLambda = array[arcPtrRaw->getPetscLocalDofIdx()];
    array[arcPtrRaw->getPetscLocalDofIdx()] = 0;
    ierr = VecRestoreArray(arcPtrRaw->dx,&array); CHKERRQ(ierr);
  }
  ierr = VecGhostUpdateBegin(arcPtrRaw->ghosTdLambda,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = VecGhostUpdateEnd(arcPtrRaw->ghosTdLambda,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  //dx2
  ierr = VecDot(arcPtrRaw->dx,arcPtrRaw->dx,&arcPtrRaw->dx2); CHKERRQ(ierr);
  PetscPrintf(
    arcPtrRaw->mField.get_comm(),"\tdlambda = %6.4e dx2 = %6.4e\n",arcPtrRaw->dLambda,arcPtrRaw->dx2
  );
  PetscFunctionReturn(0);
}

PetscErrorCode SphericalArcLengthControl::calculateInitDlambda(double *dlambda) {
  PetscFunctionBegin;
  *dlambda = sqrt(pow(arcPtrRaw->s,2)/(pow(arcPtrRaw->beta,2)*arcPtrRaw->F_lambda2));
  if(!(*dlambda == *dlambda)) {
    std::ostringstream sss;
    sss << "s " << arcPtrRaw->s << " " << arcPtrRaw->beta << " " << arcPtrRaw->F_lambda2;
    SETERRQ(PETSC_COMM_SELF,MOFEM_IMPOSIBLE_CASE,sss.str().c_str());
  }
  PetscFunctionReturn(0);
}

PetscErrorCode SphericalArcLengthControl::setDlambdaToX(Vec x,double dlambda) {
  PetscFunctionBegin;
  //check if local dof idx is non zero, i.e. that lambda is accessible from this processor
  if(arcPtrRaw->getPetscLocalDofIdx()!=-1) {
    double *array;
    ierr = VecGetArray(x,&array); CHKERRQ(ierr);
    double lambda_old = array[arcPtrRaw->getPetscLocalDofIdx()];
    if(!(dlambda == dlambda)) {
      std::ostringstream sss;
      sss << "s " << arcPtrRaw->s << " " << arcPtrRaw->beta << " " << arcPtrRaw->F_lambda2;
      SETERRQ(PETSC_COMM_SELF,1,sss.str().c_str());
    }
    array[arcPtrRaw->getPetscLocalDofIdx()] = lambda_old + dlambda;
    PetscPrintf(
      arcPtrRaw->mField.get_comm(),
      "\tlambda = %6.4e, %6.4e (%6.4e)\n",
      lambda_old,array[arcPtrRaw->getPetscLocalDofIdx()],dlambda
    );
    ierr = VecRestoreArray(x,&array); CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}
