/** \file LobattoPolynomial.cpp
 * \brief implementation of multi-grid solver for p- adaptivity
 *
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
 * License along with MoFEM. If not, see <http://www.gnu.org/licenses/>
*/

#include <version.h>
#include <config.h>
#include <definitions.h>
#include <Includes.hpp>

#include <base_functions.h>
#include <Common.hpp>
#include <UnknownInterface.hpp>
using namespace MoFEM;

#include <BaseFunction.hpp>
#include <LegendrePolynomial.hpp>
#include <LobattoPolynomial.hpp>

PetscErrorCode LobattoPolynomialCtx::queryInterface(
  const MOFEMuuid& uuid,MoFEM::UnknownInterface** iface
) {
  PetscErrorCode ierr;
  PetscFunctionBegin;
  *iface = NULL;
  if(uuid == IDD_LOBATTO_BASE_FUNCTION) {
    *iface = dynamic_cast<LobattoPolynomialCtx*>(this);
    PetscFunctionReturn(0);
  } else {
    SETERRQ(PETSC_COMM_WORLD,MOFEM_DATA_INCONSISTENCY,"wrong interference");
  }
  ierr = LegendrePolynomialCtx::queryInterface(uuid,iface); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode LobattoPolynomial::queryInterface(
  const MOFEMuuid& uuid,MoFEM::UnknownInterface** iface
) {
  PetscErrorCode ierr;
  PetscFunctionBegin;
  *iface = NULL;
  if(uuid == IDD_LOBATTO_BASE_FUNCTION) {
    *iface = dynamic_cast<LobattoPolynomial*>(this);
    PetscFunctionReturn(0);
  } else {
    SETERRQ(PETSC_COMM_WORLD,MOFEM_DATA_INCONSISTENCY,"wrong interference");
  }
  ierr = LegendrePolynomial::queryInterface(uuid,iface); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode LobattoPolynomial::getValue(
  ublas::matrix<double> &pts,
  boost::shared_ptr<BaseFunctionCtx> ctx_ptr
) {
  PetscErrorCode ierr;
  PetscFunctionBegin;
  MoFEM::UnknownInterface *iface;
  ierr = ctx_ptr->queryInterface(IDD_LOBATTO_BASE_FUNCTION,&iface); CHKERRQ(ierr);
  LobattoPolynomialCtx *ctx = reinterpret_cast<LobattoPolynomialCtx*>(iface);
  // Polynomial order start from 2nd order
  ctx->baseFunPtr->resize(pts.size2(),ctx->P+1-2,false);
  ctx->baseDiffFunPtr->resize(pts.size2(),ctx->dIm*(ctx->P+1-2),false);
  double *l = NULL;
  double *diff_l = NULL;
  for(unsigned int gg = 0;gg<pts.size2();gg++) {
    if(ctx->baseFunPtr) l = &((*ctx->baseFunPtr)(gg,0));
    if(ctx->baseDiffFunPtr) diff_l = &((*ctx->baseDiffFunPtr)(gg,0));
    ierr = (ctx->basePolynomials)(ctx->P,pts(0,gg),ctx->diffS,l,diff_l,ctx->dIm); CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}

PetscErrorCode KernelLobattoPolynomialCtx::queryInterface(
  const MOFEMuuid& uuid,MoFEM::UnknownInterface** iface
) {
  PetscErrorCode ierr;
  PetscFunctionBegin;
  *iface = NULL;
  if(uuid == IDD_KERNEL_BASE_FUNCTION) {
    *iface = dynamic_cast<KernelLobattoPolynomialCtx*>(this);
    PetscFunctionReturn(0);
  } else {
    SETERRQ(PETSC_COMM_WORLD,MOFEM_DATA_INCONSISTENCY,"wrong interference");
  }
  ierr = LegendrePolynomialCtx::queryInterface(uuid,iface); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode KernelLobattoPolynomial::queryInterface(
  const MOFEMuuid& uuid,MoFEM::UnknownInterface** iface
) {
  PetscErrorCode ierr;
  PetscFunctionBegin;
  *iface = NULL;
  if(uuid == IDD_KERNEL_BASE_FUNCTION) {
    *iface = dynamic_cast<KernelLobattoPolynomial*>(this);
    PetscFunctionReturn(0);
  } else {
    SETERRQ(PETSC_COMM_WORLD,MOFEM_DATA_INCONSISTENCY,"wrong interference");
  }
  ierr = LegendrePolynomial::queryInterface(uuid,iface); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode KernelLobattoPolynomial::getValue(
  ublas::matrix<double> &pts,
  boost::shared_ptr<BaseFunctionCtx> ctx_ptr
) {
  PetscErrorCode ierr;
  PetscFunctionBegin;
  MoFEM::UnknownInterface *iface;
  ierr = ctx_ptr->queryInterface(IDD_KERNEL_BASE_FUNCTION,&iface); CHKERRQ(ierr);
  LobattoPolynomialCtx *ctx = reinterpret_cast<LobattoPolynomialCtx*>(iface);
  ctx->baseFunPtr->resize(pts.size2(),ctx->P+1,false);
  ctx->baseDiffFunPtr->resize(pts.size2(),ctx->dIm*(ctx->P+1),false);
  double *l = NULL;
  double *diff_l = NULL;
  for(unsigned int gg = 0;gg<pts.size2();gg++) {
    if(ctx->baseFunPtr) l = &((*ctx->baseFunPtr)(gg,0));
    if(ctx->baseDiffFunPtr) diff_l = &((*ctx->baseDiffFunPtr)(gg,0));
    ierr = (ctx->basePolynomials)(ctx->P,pts(0,gg),ctx->diffS,l,diff_l,ctx->dIm); CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}
