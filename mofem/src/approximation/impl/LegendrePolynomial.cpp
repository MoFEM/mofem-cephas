/** \file LegendrePolynomial.cpp
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

PetscErrorCode LegendrePolynomialCtx::queryInterface(
  const MOFEMuuid& uuid,MoFEM::UnknownInterface** iface
) {
  PetscErrorCode ierr;
  PetscFunctionBegin;
  *iface = NULL;
  if(uuid == IDD_LEGENDRE_BASE_FUNCTION) {
    *iface = dynamic_cast<LegendrePolynomialCtx*>(this);
    PetscFunctionReturn(0);
  } else {
    SETERRQ(PETSC_COMM_WORLD,MOFEM_DATA_INCONSISTENCY,"wrong interference");
  }
  ierr = BaseFunctionCtx::queryInterface(uuid,iface); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode LegendrePolynomial::queryInterface(
  const MOFEMuuid& uuid,MoFEM::UnknownInterface** iface
) {
  PetscErrorCode ierr;
  PetscFunctionBegin;
  *iface = NULL;
  if(uuid == IDD_LEGENDRE_BASE_FUNCTION) {
    *iface = dynamic_cast<LegendrePolynomial*>(this);
    PetscFunctionReturn(0);
  } else {
    SETERRQ(PETSC_COMM_WORLD,MOFEM_DATA_INCONSISTENCY,"wrong interference");
  }
  ierr = BaseFunction::queryInterface(uuid,iface); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode LegendrePolynomial::getValue(
  ublas::matrix<double> &pTs,
  boost::shared_ptr<ublas::matrix<double> > baseFunPtr,
  boost::shared_ptr<ublas::matrix<double> > baseDiffFunPtr,
  boost::shared_ptr<BaseFunctionCtx> ctxPtr
) {
  PetscErrorCode ierr;
  PetscFunctionBegin;
  MoFEM::UnknownInterface *iface;
  ierr = ctxPtr->queryInterface(IDD_LEGENDRE_BASE_FUNCTION,&iface); CHKERRQ(ierr);
  LegendrePolynomialCtx *ctx = reinterpret_cast<LegendrePolynomialCtx*>(iface);
  baseFunPtr->resize(1,ctx->P+1,false);
  baseDiffFunPtr->resize(ctx->dIm,ctx->P+1,false);
  double *l = NULL;
  double *diff_l = NULL;
  if(baseFunPtr) l = &*baseFunPtr->data().begin();
  if(baseDiffFunPtr) diff_l = &*baseDiffFunPtr->data().begin();
  for(int gg = 0;gg<pTs.size2();gg++) {
    ierr = (ctx->base_polynomials)(ctx->P,pTs(0,gg),ctx->diffS,l,diff_l,ctx->dIm); CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}
