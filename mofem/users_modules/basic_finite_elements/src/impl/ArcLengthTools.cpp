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

MoFEMErrorCode ArcLengthCtx::setS(double s) {
  MoFEMFunctionBeginHot;
  this->s = s;
  ierr = PetscPrintf(mField.get_comm(), "\tSet s = %6.4e\n", this->s);
  CHKERRG(ierr);
  MoFEMFunctionReturnHot(0);
}

MoFEMErrorCode ArcLengthCtx::setAlphaBeta(double alpha, double beta) {
  MoFEMFunctionBeginHot;
  this->alpha = alpha;
  this->beta = beta;
  ierr = PetscPrintf(mField.get_comm(), "\tSet alpha = %6.4e beta = %6.4e\n",
                     this->alpha, this->beta);
  CHKERRG(ierr);
  MoFEMFunctionReturnHot(0);
}

ArcLengthCtx::ArcLengthCtx(MoFEM::Interface &m_field,
                           const std::string &problem_name,
                           const std::string &field_name)
    : mField(m_field), dx2(0), F_lambda2(0), res_lambda(0) {
  ierr = m_field.getInterface<VecManager>()->vecCreateGhost(problem_name, ROW,
                                                            &F_lambda);
  CHKERRABORT(PETSC_COMM_WORLD, ierr);
  ierr = VecSetOption(F_lambda, VEC_IGNORE_NEGATIVE_INDICES, PETSC_TRUE);
  CHKERRABORT(PETSC_COMM_WORLD, ierr);
  ierr = VecDuplicate(F_lambda, &db);
  CHKERRABORT(PETSC_COMM_WORLD, ierr);
  ierr = VecDuplicate(F_lambda, &xLambda);
  CHKERRABORT(PETSC_COMM_WORLD, ierr);
  ierr = VecDuplicate(F_lambda, &x0);
  CHKERRABORT(PETSC_COMM_WORLD, ierr);
  ierr = VecDuplicate(F_lambda, &dx);
  CHKERRABORT(PETSC_COMM_WORLD, ierr);

  const Problem *problem_ptr;
  ierr = m_field.get_problem(problem_name, &problem_ptr);
  CHKERRABORT(PETSC_COMM_WORLD, ierr);
  boost::shared_ptr<NumeredDofEntity_multiIndex> dofs_ptr_no_const =
      problem_ptr->getNumeredDofsRows();
  NumeredDofEntityByFieldName::iterator hi_dit;
  dIt = dofs_ptr_no_const->get<FieldName_mi_tag>().lower_bound(field_name);
  hi_dit = dofs_ptr_no_const->get<FieldName_mi_tag>().upper_bound(field_name);

  if (distance(dIt, hi_dit) != 1) {
    SETERRABORT(PETSC_COMM_WORLD, MOFEM_DATA_INCONSISTENCY,
                "can not find unique LAMBDA (load factor)");
  }

  if ((unsigned int)mField.get_comm_rank() == (*dIt)->getPart()) {
    ierr = VecCreateGhostWithArray(mField.get_comm(), 1, 1, 0, PETSC_NULL,
                                   &dLambda, &ghosTdLambda);
    CHKERRABORT(PETSC_COMM_WORLD, ierr);
    ierr = VecCreateGhostWithArray(mField.get_comm(), 1, 1, 0, PETSC_NULL,
                                   &dIag, &ghostDiag);
    CHKERRABORT(PETSC_COMM_WORLD, ierr);
  } else {
    int one[] = {0};
    ierr = VecCreateGhostWithArray(mField.get_comm(), 0, 1, 1, one, &dLambda,
                                   &ghosTdLambda);
    CHKERRABORT(PETSC_COMM_WORLD, ierr);
    ierr = VecCreateGhostWithArray(mField.get_comm(), 0, 1, 1, one, &dIag,
                                   &ghostDiag);
    CHKERRABORT(PETSC_COMM_WORLD, ierr);
  }
  dLambda = 0;
  dIag = 0;
}

ArcLengthCtx::~ArcLengthCtx() {
  ierr = VecDestroy(&F_lambda);
  CHKERRABORT(PETSC_COMM_WORLD, ierr);
  ierr = VecDestroy(&db);
  CHKERRABORT(PETSC_COMM_WORLD, ierr);
  ierr = VecDestroy(&xLambda);
  CHKERRABORT(PETSC_COMM_WORLD, ierr);
  ierr = VecDestroy(&x0);
  CHKERRABORT(PETSC_COMM_WORLD, ierr);
  ierr = VecDestroy(&dx);
  CHKERRABORT(PETSC_COMM_WORLD, ierr);
  ierr = VecDestroy(&ghosTdLambda);
  CHKERRABORT(PETSC_COMM_WORLD, ierr);
  ierr = VecDestroy(&ghostDiag);
  CHKERRABORT(PETSC_COMM_WORLD, ierr);
}

// ***********************
// Arc-length shell matrix

ArcLengthMatShell::ArcLengthMatShell(Mat aij, ArcLengthCtx *arc_ptr_raw,
                                     string problem_name)
    : Aij(aij), arcPtrRaw(arc_ptr_raw), problemName(problem_name) {
  ierr = PetscObjectReference((PetscObject)aij);
  CHKERRABORT(PETSC_COMM_WORLD, ierr);
}

ArcLengthMatShell::ArcLengthMatShell(Mat aij,
                                     boost::shared_ptr<ArcLengthCtx> arc_ptr,
                                     string problem_name)
    : Aij(aij), arcPtrRaw(arc_ptr.get()), problemName(problem_name),
      arcPtr(arc_ptr) {
  ierr = PetscObjectReference((PetscObject)aij);
  CHKERRABORT(PETSC_COMM_WORLD, ierr);
}

ArcLengthMatShell::~ArcLengthMatShell() {
  ierr = MatDestroy(&Aij);
  CHKERRABORT(PETSC_COMM_WORLD, ierr);
}

MoFEMErrorCode ArcLengthMatShell::setLambda(Vec ksp_x, double *lambda,
                                            ScatterMode scattermode) {
  MoFEMFunctionBeginHot;

  int part = arcPtrRaw->getPart();
  int rank = arcPtrRaw->mField.get_comm_rank();

  Vec lambda_ghost;
  if (rank == part) {
    ierr = VecCreateGhostWithArray(arcPtrRaw->mField.get_comm(), 1, 1, 0,
                                   PETSC_NULL, lambda, &lambda_ghost);
    CHKERRG(ierr);
  } else {
    int one[] = {0};
    ierr = VecCreateGhostWithArray(arcPtrRaw->mField.get_comm(), 0, 1, 1, one,
                                   lambda, &lambda_ghost);
    CHKERRG(ierr);
  }

  switch (scattermode) {
  case SCATTER_FORWARD: {
    int idx = arcPtrRaw->getPetscGlobalDofIdx();
    if (part == rank) {
      ierr = VecGetValues(ksp_x, 1, &idx, lambda);
      CHKERRG(ierr);
    }
    ierr = VecGhostUpdateBegin(lambda_ghost, INSERT_VALUES, SCATTER_FORWARD);
    CHKERRG(ierr);
    ierr = VecGhostUpdateEnd(lambda_ghost, INSERT_VALUES, SCATTER_FORWARD);
    CHKERRG(ierr);
  } break;
  case SCATTER_REVERSE: {
    // ierr = VecGhostUpdateBegin(lambda_ghost,INSERT_VALUES,SCATTER_FORWARD);
    // CHKERRG(ierr);
    // ierr = VecGhostUpdateEnd(lambda_ghost,INSERT_VALUES,SCATTER_FORWARD);
    // CHKERRG(ierr);
    if (part == rank) {
      PetscScalar *array;
      ierr = VecGetArray(ksp_x, &array);
      CHKERRG(ierr);
      array[arcPtrRaw->getPetscLocalDofIdx()] = *lambda;
      ierr = VecRestoreArray(ksp_x, &array);
      CHKERRG(ierr);
    }
  } break;
  default:
    SETERRQ(PETSC_COMM_SELF, MOFEM_NOT_IMPLEMENTED, "not implemented");
  }

  ierr = VecDestroy(&lambda_ghost);
  CHKERRG(ierr);

  MoFEMFunctionReturnHot(0);
}

MoFEMErrorCode ArcLengthMatMultShellOp(Mat A, Vec x, Vec f) {
  MoFEMFunctionBeginHot;
  void *void_ctx;
  ierr = MatShellGetContext(A, &void_ctx);
  CHKERRG(ierr);
  ArcLengthMatShell *ctx = (ArcLengthMatShell *)void_ctx;
  ierr = MatMult(ctx->Aij, x, f);
  CHKERRG(ierr);
  double lambda;
  ierr = ctx->setLambda(x, &lambda, SCATTER_FORWARD);
  CHKERRG(ierr);
  double db_dot_x;
  ierr = VecDot(ctx->arcPtrRaw->db, x, &db_dot_x);
  CHKERRG(ierr);
  double f_lambda;
  f_lambda = ctx->arcPtrRaw->dIag * lambda + db_dot_x;
  ierr = ctx->setLambda(f, &f_lambda, SCATTER_REVERSE);
  CHKERRG(ierr);
  ierr = VecAXPY(f, lambda, ctx->arcPtrRaw->F_lambda);
  CHKERRG(ierr);
  MoFEMFunctionReturnHot(0);
}

// arc-length preconditioner

PCArcLengthCtx::PCArcLengthCtx(Mat shell_Aij, Mat aij, ArcLengthCtx *arc_ptr)
    : kSP(PETSC_NULL), pC(PETSC_NULL), shellAij(shell_Aij), Aij(aij),
      arcPtrRaw(arc_ptr) {
  ierr = PCCreate(PetscObjectComm((PetscObject)aij), &pC);
  CHKERRABORT(PETSC_COMM_WORLD, ierr);
  ierr = KSPCreate(PetscObjectComm((PetscObject)pC), &kSP);
  CHKERRABORT(PETSC_COMM_WORLD, ierr);
  ierr = KSPAppendOptionsPrefix(kSP, "arc_length_");
  CHKERRABORT(PETSC_COMM_WORLD, ierr);
}

PCArcLengthCtx::PCArcLengthCtx(Mat shell_Aij, Mat aij,
                               boost::shared_ptr<ArcLengthCtx> &arc_ptr)
    : kSP(PETSC_NULL), shellAij(shell_Aij), Aij(aij), arcPtrRaw(arc_ptr.get()),
      arcPtr(arc_ptr) {
  ierr = PCCreate(PetscObjectComm((PetscObject)aij), &pC);
  CHKERRABORT(PETSC_COMM_WORLD, ierr);
  ierr = KSPCreate(PetscObjectComm((PetscObject)pC), &kSP);
  CHKERRABORT(PETSC_COMM_WORLD, ierr);
  ierr = KSPAppendOptionsPrefix(kSP, "arc_length_");
  CHKERRABORT(PETSC_COMM_WORLD, ierr);
}

PCArcLengthCtx::PCArcLengthCtx(PC pc, Mat shell_Aij, Mat aij,
                               boost::shared_ptr<ArcLengthCtx> &arc_ptr)
    : kSP(PETSC_NULL), pC(pc), shellAij(shell_Aij), Aij(aij),
      arcPtrRaw(arc_ptr.get()), arcPtr(arc_ptr) {
  ierr = PetscObjectReference((PetscObject)pC);
  CHKERRABORT(PETSC_COMM_WORLD, ierr);
  ierr = KSPCreate(PetscObjectComm((PetscObject)pC), &kSP);
  CHKERRABORT(PETSC_COMM_WORLD, ierr);
  ierr = KSPAppendOptionsPrefix(kSP, "arc_length_");
  CHKERRABORT(PETSC_COMM_WORLD, ierr);
}

PCArcLengthCtx::~PCArcLengthCtx() {
  if (kSP != PETSC_NULL) {
    ierr = KSPDestroy(&kSP);
    CHKERRABORT(PETSC_COMM_WORLD, ierr);
  }
  if (pC != PETSC_NULL) {
    ierr = PCDestroy(&pC);
    CHKERRABORT(PETSC_COMM_WORLD, ierr);
  }
}

MoFEMErrorCode PCApplyArcLength(PC pc, Vec pc_f, Vec pc_x) {
  MoFEMFunctionBeginHot;
  void *void_ctx;
  ierr = PCShellGetContext(pc, &void_ctx);
  CHKERRG(ierr);
  PCArcLengthCtx *ctx = (PCArcLengthCtx *)void_ctx;
  void *void_MatCtx;
  MatShellGetContext(ctx->shellAij, &void_MatCtx);
  ArcLengthMatShell *mat_ctx = (ArcLengthMatShell *)void_MatCtx;
  ierr = KSPSetInitialGuessNonzero(ctx->kSP, PETSC_FALSE);
  CHKERRG(ierr);
  ierr = KSPSolve(ctx->kSP, pc_f, pc_x);
  CHKERRG(ierr);
  PetscBool same;
  PetscObjectTypeCompare((PetscObject)ctx->kSP, KSPPREONLY, &same);
  if (same != PETSC_TRUE) {
    ierr = KSPSetInitialGuessNonzero(ctx->kSP, PETSC_TRUE);
    CHKERRG(ierr);
  }
  ierr = KSPSolve(ctx->kSP, ctx->arcPtrRaw->F_lambda, ctx->arcPtrRaw->xLambda);
  CHKERRG(ierr);
  double db_dot_pc_x, db_dot_x_lambda;
  ierr = VecDot(ctx->arcPtrRaw->db, pc_x, &db_dot_pc_x);
  CHKERRG(ierr);
  ierr = VecDot(ctx->arcPtrRaw->db, ctx->arcPtrRaw->xLambda, &db_dot_x_lambda);
  CHKERRG(ierr);
  double denominator = ctx->arcPtrRaw->dIag + db_dot_x_lambda;
  double res_lambda;
  ierr = mat_ctx->setLambda(pc_f, &res_lambda, SCATTER_FORWARD);
  CHKERRG(ierr);
  double ddlambda = (res_lambda - db_dot_pc_x) / denominator;
  // cerr << denominator << " " << res_lambda << " " << ddlambda << endl;
  if (ddlambda != ddlambda || denominator == 0) {
    double nrm2_pc_f, nrm2_db, nrm2_pc_x, nrm2_xLambda;
    ierr = VecNorm(pc_f, NORM_2, &nrm2_pc_f);
    CHKERRG(ierr);
    ierr = VecNorm(ctx->arcPtrRaw->db, NORM_2, &nrm2_db);
    CHKERRG(ierr);
    ierr = VecNorm(pc_x, NORM_2, &nrm2_pc_x);
    CHKERRG(ierr);
    ierr = VecNorm(ctx->arcPtrRaw->xLambda, NORM_2, &nrm2_xLambda);
    CHKERRG(ierr);
    std::ostringstream ss;
    ss << "problem with ddlambda=" << res_lambda << " res_lambda=" << res_lambda
       << " denominator=" << denominator << " ddlamnda=" << ddlambda
       << " db_dot_pc_x=" << db_dot_pc_x
       << " db_dot_x_lambda=" << db_dot_x_lambda
       << " diag=" << ctx->arcPtrRaw->dIag << " nrm2_db=" << nrm2_db
       << " nrm2_pc_f=" << nrm2_pc_f << " nrm2_pc_x=" << nrm2_pc_x
       << " nrm2_xLambda=" << nrm2_xLambda;
    SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY, ss.str().c_str());
  }
  ierr = VecAXPY(pc_x, ddlambda, ctx->arcPtrRaw->xLambda);
  CHKERRG(ierr);
  ierr = mat_ctx->setLambda(pc_x, &ddlambda, SCATTER_REVERSE);
  CHKERRG(ierr);
  MoFEMFunctionReturnHot(0);
}

MoFEMErrorCode PCSetupArcLength(PC pc) {
  MoFEMFunctionBeginHot;
  void *void_ctx;
  ierr = PCShellGetContext(pc, &void_ctx);
  CHKERRG(ierr);
  PCArcLengthCtx *ctx = (PCArcLengthCtx *)void_ctx;
  ierr = PCSetFromOptions(ctx->pC);
  CHKERRG(ierr);
  ierr = PCGetOperators(pc, &ctx->shellAij, &ctx->Aij);
  CHKERRG(ierr);
  ierr = PCSetOperators(ctx->pC, ctx->Aij, ctx->Aij);
  CHKERRG(ierr);
  ierr = PCSetUp(ctx->pC);
  CHKERRG(ierr);
  // SetUp PC KSP solver
  ierr = KSPSetType(ctx->kSP, KSPPREONLY);
  CHKERRG(ierr);
  ierr = KSPSetTabLevel(ctx->kSP, 3);
  CHKERRG(ierr);
  ierr = KSPSetFromOptions(ctx->kSP);
  CHKERRG(ierr);
  ierr = KSPSetOperators(ctx->kSP, ctx->Aij, ctx->Aij);
  CHKERRG(ierr);
  ierr = KSPSetPC(ctx->kSP, ctx->pC);
  CHKERRG(ierr);
  ierr = KSPSetUp(ctx->kSP);
  CHKERRG(ierr);
  MoFEMFunctionReturnHot(0);
}

// ***********************
// Zero F_lambda vector

ZeroFLmabda::ZeroFLmabda(boost::shared_ptr<ArcLengthCtx> arc_ptr)
    : arcPtr(arc_ptr) {}

MoFEMErrorCode ZeroFLmabda::preProcess() {
  MoFEMFunctionBeginHot;
  switch (snes_ctx) {
  case CTX_SNESSETFUNCTION: {
    ierr = VecZeroEntries(arcPtr->F_lambda);
    CHKERRG(ierr);
    ierr =
        VecGhostUpdateBegin(arcPtr->F_lambda, INSERT_VALUES, SCATTER_FORWARD);
    CHKERRG(ierr);
    ierr = VecGhostUpdateEnd(arcPtr->F_lambda, INSERT_VALUES, SCATTER_FORWARD);
    CHKERRG(ierr);
  } break;
  default:
    SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY, "Impossible case");
  }
  MoFEMFunctionReturnHot(0);
}

AssembleFlambda::AssembleFlambda(boost::shared_ptr<ArcLengthCtx> arc_ptr,
                                 boost::shared_ptr<DirichletDisplacementBc> bc)
    : arcPtr(arc_ptr), bC(bc) {}

MoFEMErrorCode AssembleFlambda::preProcess() {
  MoFEMFunctionBeginHot;
  MoFEMFunctionReturnHot(0);
}
MoFEMErrorCode AssembleFlambda::operator()() {
  MoFEMFunctionBeginHot;
  MoFEMFunctionReturnHot(0);
}
MoFEMErrorCode AssembleFlambda::postProcess() {
  MoFEMFunctionBeginHot;
  switch (snes_ctx) {
  case CTX_SNESSETFUNCTION: {
    // F_lambda
    ierr = VecAssemblyBegin(arcPtr->F_lambda);
    CHKERRG(ierr);
    ierr = VecAssemblyEnd(arcPtr->F_lambda);
    CHKERRG(ierr);
    ierr = VecGhostUpdateBegin(arcPtr->F_lambda, ADD_VALUES, SCATTER_REVERSE);
    CHKERRG(ierr);
    ierr = VecGhostUpdateEnd(arcPtr->F_lambda, ADD_VALUES, SCATTER_REVERSE);
    CHKERRG(ierr);
    if (bC) {
      for (std::vector<int>::iterator vit = bC->dofsIndices.begin();
           vit != bC->dofsIndices.end(); vit++) {
        ierr = VecSetValue(arcPtr->F_lambda, *vit, 0, INSERT_VALUES);
        CHKERRG(ierr);
      }
      ierr = VecAssemblyBegin(arcPtr->F_lambda);
      CHKERRG(ierr);
      ierr = VecAssemblyEnd(arcPtr->F_lambda);
      CHKERRG(ierr);
    }
    ierr = VecDot(arcPtr->F_lambda, arcPtr->F_lambda, &arcPtr->F_lambda2);
    CHKERRG(ierr);
    PetscPrintf(PETSC_COMM_WORLD, "\tF_lambda2 = %6.4e\n", arcPtr->F_lambda2);
    // add F_lambda
    ierr = VecAssemblyBegin(snes_f);
    CHKERRG(ierr);
    ierr = VecAssemblyEnd(snes_f);
    CHKERRG(ierr);
    double lambda = arcPtr->getFieldData();
    ierr = VecAXPY(snes_f, lambda, arcPtr->F_lambda);
    CHKERRG(ierr);
    double fnorm;
    ierr = VecNorm(snes_f, NORM_2, &fnorm);
    CHKERRG(ierr);
    PetscPrintf(PETSC_COMM_WORLD, "\tfnorm = %6.4e lambda = %6.4g\n", fnorm,
                lambda);
  } break;
  default:
    SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY, "Impossible case");
  }
  MoFEMFunctionReturnHot(0);
}

// ************************
// Simple arc-length method

SimpleArcLengthControl::SimpleArcLengthControl(
    boost::shared_ptr<ArcLengthCtx> &arc_ptr, const bool assemble)
    : FEMethod(), arcPtr(arc_ptr), aSsemble(assemble) {}

SimpleArcLengthControl::~SimpleArcLengthControl() {}

MoFEMErrorCode SimpleArcLengthControl::preProcess() {
  MoFEMFunctionBeginHot;
  switch (snes_ctx) {
  case CTX_SNESSETFUNCTION: {
    if (aSsemble) {
      ierr = VecAssemblyBegin(snes_f);
      CHKERRG(ierr);
      ierr = VecAssemblyEnd(snes_f);
      CHKERRG(ierr);
    }
    ierr = calculateDxAndDlambda(snes_x);
    CHKERRG(ierr);
    ierr = calculateDb();
    CHKERRG(ierr);
  } break;
  case CTX_SNESSETJACOBIAN: {
    if (aSsemble) {
      ierr = MatAssemblyBegin(snes_B, MAT_FLUSH_ASSEMBLY);
      CHKERRG(ierr);
      ierr = MatAssemblyEnd(snes_B, MAT_FLUSH_ASSEMBLY);
      CHKERRG(ierr);
    }
  } break;
  default:
    break;
  }
  MoFEMFunctionReturnHot(0);
}

MoFEMErrorCode SimpleArcLengthControl::operator()() {
  MoFEMFunctionBeginHot;
  switch (snes_ctx) {
  case CTX_SNESSETFUNCTION: {
    arcPtr->res_lambda = calculateLambdaInt() - arcPtr->s;
    ierr = VecSetValue(snes_f, arcPtr->getPetscGlobalDofIdx(),
                       arcPtr->res_lambda, ADD_VALUES);
    CHKERRG(ierr);
  } break;
  case CTX_SNESSETJACOBIAN: {
    arcPtr->dIag = arcPtr->beta;
    ierr = MatSetValue(snes_B, arcPtr->getPetscGlobalDofIdx(),
                       arcPtr->getPetscGlobalDofIdx(), 1, ADD_VALUES);
    CHKERRG(ierr);
  } break;
  default:
    break;
  }
  MoFEMFunctionReturnHot(0);
}

MoFEMErrorCode SimpleArcLengthControl::postProcess() {
  MoFEMFunctionBeginHot;
  switch (snes_ctx) {
  case CTX_SNESSETFUNCTION: {
    if (aSsemble) {
      ierr = VecAssemblyBegin(snes_f);
      CHKERRG(ierr);
      ierr = VecAssemblyEnd(snes_f);
      CHKERRG(ierr);
    }
  } break;
  case CTX_SNESSETJACOBIAN: {
    if (aSsemble) {
      ierr = MatAssemblyBegin(snes_B, MAT_FLUSH_ASSEMBLY);
      CHKERRG(ierr);
      ierr = MatAssemblyEnd(snes_B, MAT_FLUSH_ASSEMBLY);
      CHKERRG(ierr);
    }
    ierr =
        VecGhostUpdateBegin(arcPtr->ghostDiag, INSERT_VALUES, SCATTER_FORWARD);
    CHKERRG(ierr);
    ierr = VecGhostUpdateEnd(arcPtr->ghostDiag, INSERT_VALUES, SCATTER_FORWARD);
    CHKERRG(ierr);
  } break;
  default:
    break;
  }
  MoFEMFunctionReturnHot(0);
}

double SimpleArcLengthControl::calculateLambdaInt() {
  return arcPtr->beta * arcPtr->dLambda;
}

MoFEMErrorCode SimpleArcLengthControl::calculateDb() {
  MoFEMFunctionBeginHot;
  ierr = VecZeroEntries(arcPtr->db);
  CHKERRG(ierr);
  ierr = VecGhostUpdateBegin(arcPtr->db, INSERT_VALUES, SCATTER_FORWARD);
  CHKERRG(ierr);
  ierr = VecGhostUpdateEnd(arcPtr->db, INSERT_VALUES, SCATTER_FORWARD);
  CHKERRG(ierr);
  MoFEMFunctionReturnHot(0);
}

MoFEMErrorCode SimpleArcLengthControl::calculateDxAndDlambda(Vec x) {
  MoFEMFunctionBeginHot;
  // Calculate dx
  ierr = VecCopy(x, arcPtr->dx);
  CHKERRG(ierr);
  ierr = VecAXPY(arcPtr->dx, -1, arcPtr->x0);
  CHKERRG(ierr);
  ierr = VecGhostUpdateBegin(arcPtr->x0, INSERT_VALUES, SCATTER_FORWARD);
  CHKERRG(ierr);
  ierr = VecGhostUpdateEnd(arcPtr->x0, INSERT_VALUES, SCATTER_FORWARD);
  CHKERRG(ierr);
  ierr = VecGhostUpdateBegin(arcPtr->dx, INSERT_VALUES, SCATTER_FORWARD);
  CHKERRG(ierr);
  ierr = VecGhostUpdateEnd(arcPtr->dx, INSERT_VALUES, SCATTER_FORWARD);
  CHKERRG(ierr);
  // Calculate dlambda
  if (arcPtr->getPetscLocalDofIdx() != -1) {
    double *array;
    ierr = VecGetArray(arcPtr->dx, &array);
    CHKERRG(ierr);
    arcPtr->dLambda = array[arcPtr->getPetscLocalDofIdx()];
    array[arcPtr->getPetscLocalDofIdx()] = 0;
    ierr = VecRestoreArray(arcPtr->dx, &array);
    CHKERRG(ierr);
  }
  ierr =
      VecGhostUpdateBegin(arcPtr->ghosTdLambda, INSERT_VALUES, SCATTER_FORWARD);
  CHKERRG(ierr);
  ierr =
      VecGhostUpdateEnd(arcPtr->ghosTdLambda, INSERT_VALUES, SCATTER_FORWARD);
  CHKERRG(ierr);
  // Calculate dx2
  ierr = VecDot(arcPtr->dx, arcPtr->dx, &arcPtr->dx2);
  CHKERRG(ierr);
  PetscPrintf(PETSC_COMM_WORLD, "\tdx2 = %6.4e\n", arcPtr->dx2);
  MoFEMFunctionReturnHot(0);
}

// ***************************
// Spherical arc-length contri

SphericalArcLengthControl::SphericalArcLengthControl(ArcLengthCtx *arc_ptr_raw)
    : FEMethod(), arcPtrRaw(arc_ptr_raw) {}

SphericalArcLengthControl::SphericalArcLengthControl(
    boost::shared_ptr<ArcLengthCtx> &arc_ptr)
    : FEMethod(), arcPtrRaw(arc_ptr.get()), arcPtr(arc_ptr) {}

SphericalArcLengthControl::~SphericalArcLengthControl() {}

MoFEMErrorCode SphericalArcLengthControl::preProcess() {
  MoFEMFunctionBeginHot;
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
  switch (snes_ctx) {
  case CTX_SNESSETFUNCTION: {
    ierr = calculateDxAndDlambda(snes_x);
    CHKERRG(ierr);
    ierr = calculateDb();
    CHKERRG(ierr);
  } break;
  case CTX_SNESSETJACOBIAN: {
  } break;
  default:
    break;
  }
  MoFEMFunctionReturnHot(0);
}

double SphericalArcLengthControl::calculateLambdaInt() {
  return arcPtrRaw->alpha * arcPtrRaw->dx2 +
         pow(arcPtrRaw->dLambda, 2) * pow(arcPtrRaw->beta, 2) *
             arcPtrRaw->F_lambda2;
}

MoFEMErrorCode SphericalArcLengthControl::calculateDb() {
  MoFEMFunctionBeginHot;
  ierr = VecCopy(arcPtrRaw->dx, arcPtrRaw->db);
  CHKERRG(ierr);
  ierr = VecScale(arcPtrRaw->db, 2 * arcPtrRaw->alpha);
  CHKERRG(ierr);
  MoFEMFunctionReturnHot(0);
}

MoFEMErrorCode SphericalArcLengthControl::operator()() {
  MoFEMFunctionBeginHot;
  switch (snes_ctx) {
  case CTX_SNESSETFUNCTION: {
    arcPtrRaw->res_lambda = calculateLambdaInt() - pow(arcPtrRaw->s, 2);
    ierr = VecSetValue(snes_f, arcPtrRaw->getPetscGlobalDofIdx(),
                       arcPtrRaw->res_lambda, ADD_VALUES);
    CHKERRG(ierr);
    PetscPrintf(arcPtrRaw->mField.get_comm(), "\tres_lambda = %6.4e\n",
                arcPtrRaw->res_lambda);
  } break;
  case CTX_SNESSETJACOBIAN: {
    arcPtrRaw->dIag =
        2 * arcPtrRaw->dLambda * pow(arcPtrRaw->beta, 2) * arcPtrRaw->F_lambda2;
    ierr = MatSetValue(snes_B, arcPtrRaw->getPetscGlobalDofIdx(),
                       arcPtrRaw->getPetscGlobalDofIdx(), 1, ADD_VALUES);
    CHKERRG(ierr);
  } break;
  default:
    break;
  }

  MoFEMFunctionReturnHot(0);
}

MoFEMErrorCode SphericalArcLengthControl::postProcess() {
  MoFEMFunctionBeginHot;
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
  switch (snes_ctx) {
  case CTX_SNESSETFUNCTION: {
    PetscPrintf(arcPtrRaw->mField.get_comm(), "\tlambda = %6.4e\n",
                arcPtrRaw->getFieldData());
  } break;
  case CTX_SNESSETJACOBIAN: {
    // VecView(arcPtrRaw->ghostDiag,PETSC_VIEWER_STDOUT_WORLD);
    ierr = VecGhostUpdateBegin(arcPtrRaw->ghostDiag, INSERT_VALUES,
                               SCATTER_FORWARD);
    CHKERRG(ierr);
    ierr =
        VecGhostUpdateEnd(arcPtrRaw->ghostDiag, INSERT_VALUES, SCATTER_FORWARD);
    CHKERRG(ierr);
    ierr = MatAssemblyBegin(snes_B, MAT_FLUSH_ASSEMBLY);
    CHKERRG(ierr);
    ierr = MatAssemblyEnd(snes_B, MAT_FLUSH_ASSEMBLY);
    CHKERRG(ierr);
    PetscPrintf(arcPtrRaw->mField.get_comm(), "\tdiag = %6.4e\n",
                arcPtrRaw->dIag);
  } break;
  default:
    break;
  }
  MoFEMFunctionReturnHot(0);
}

MoFEMErrorCode SphericalArcLengthControl::calculateDxAndDlambda(Vec x) {
  MoFEMFunctionBeginHot;
  // dx
  ierr = VecCopy(x, arcPtrRaw->dx);
  CHKERRG(ierr);
  ierr = VecAXPY(arcPtrRaw->dx, -1, arcPtrRaw->x0);
  CHKERRG(ierr);
  ierr = VecGhostUpdateBegin(arcPtrRaw->dx, INSERT_VALUES, SCATTER_FORWARD);
  CHKERRG(ierr);
  ierr = VecGhostUpdateEnd(arcPtrRaw->dx, INSERT_VALUES, SCATTER_FORWARD);
  CHKERRG(ierr);
  // dlambda
  if (arcPtrRaw->getPetscLocalDofIdx() != -1) {
    double *array;
    ierr = VecGetArray(arcPtrRaw->dx, &array);
    CHKERRG(ierr);
    arcPtrRaw->dLambda = array[arcPtrRaw->getPetscLocalDofIdx()];
    array[arcPtrRaw->getPetscLocalDofIdx()] = 0;
    ierr = VecRestoreArray(arcPtrRaw->dx, &array);
    CHKERRG(ierr);
  }
  ierr = VecGhostUpdateBegin(arcPtrRaw->ghosTdLambda, INSERT_VALUES,
                             SCATTER_FORWARD);
  CHKERRG(ierr);
  ierr = VecGhostUpdateEnd(arcPtrRaw->ghosTdLambda, INSERT_VALUES,
                           SCATTER_FORWARD);
  CHKERRG(ierr);
  // dx2
  ierr = VecDot(arcPtrRaw->dx, arcPtrRaw->dx, &arcPtrRaw->dx2);
  CHKERRG(ierr);
  PetscPrintf(arcPtrRaw->mField.get_comm(), "\tdlambda = %6.4e dx2 = %6.4e\n",
              arcPtrRaw->dLambda, arcPtrRaw->dx2);
  MoFEMFunctionReturnHot(0);
}

MoFEMErrorCode
SphericalArcLengthControl::calculateInitDlambda(double *dlambda) {
  MoFEMFunctionBeginHot;
  *dlambda = sqrt(pow(arcPtrRaw->s, 2) /
                  (pow(arcPtrRaw->beta, 2) * arcPtrRaw->F_lambda2));
  if (!(*dlambda == *dlambda)) {
    std::ostringstream sss;
    sss << "s " << arcPtrRaw->s << " " << arcPtrRaw->beta << " "
        << arcPtrRaw->F_lambda2;
    SETERRQ(PETSC_COMM_SELF, MOFEM_IMPOSIBLE_CASE, sss.str().c_str());
  }
  MoFEMFunctionReturnHot(0);
}

MoFEMErrorCode SphericalArcLengthControl::setDlambdaToX(Vec x, double dlambda) {
  MoFEMFunctionBeginHot;
  // check if local dof idx is non zero, i.e. that lambda is accessible from
  // this processor
  if (arcPtrRaw->getPetscLocalDofIdx() != -1) {
    double *array;
    ierr = VecGetArray(x, &array);
    CHKERRG(ierr);
    double lambda_old = array[arcPtrRaw->getPetscLocalDofIdx()];
    if (!(dlambda == dlambda)) {
      std::ostringstream sss;
      sss << "s " << arcPtrRaw->s << " " << arcPtrRaw->beta << " "
          << arcPtrRaw->F_lambda2;
      SETERRQ(PETSC_COMM_SELF, 1, sss.str().c_str());
    }
    array[arcPtrRaw->getPetscLocalDofIdx()] = lambda_old + dlambda;
    PetscPrintf(arcPtrRaw->mField.get_comm(),
                "\tlambda = %6.4e, %6.4e (%6.4e)\n", lambda_old,
                array[arcPtrRaw->getPetscLocalDofIdx()], dlambda);
    ierr = VecRestoreArray(x, &array);
    CHKERRG(ierr);
  }
  MoFEMFunctionReturnHot(0);
}
