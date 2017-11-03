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

namespace MoFEM {

MoFEMErrorCode LegendrePolynomialCtx::query_interface(
  const MOFEMuuid& uuid,MoFEM::UnknownInterface** iface
) const {

  MoFEMFunctionBeginHot;
  *iface = NULL;
  if(uuid == IDD_LEGENDRE_BASE_FUNCTION) {
    *iface = const_cast<LegendrePolynomialCtx*>(this);
    MoFEMFunctionReturnHot(0);
  } else {
    SETERRQ(PETSC_COMM_WORLD,MOFEM_DATA_INCONSISTENCY,"wrong interference");
  }
  ierr = BaseFunctionCtx::query_interface(uuid,iface); CHKERRQ(ierr);
  MoFEMFunctionReturnHot(0);
}

MoFEMErrorCode LegendrePolynomial::query_interface(
  const MOFEMuuid& uuid,MoFEM::UnknownInterface** iface
) const {

  MoFEMFunctionBeginHot;
  *iface = NULL;
  if(uuid == IDD_LEGENDRE_BASE_FUNCTION) {
    *iface = const_cast<LegendrePolynomial*>(this);
    MoFEMFunctionReturnHot(0);
  } else {
    SETERRQ(PETSC_COMM_WORLD,MOFEM_DATA_INCONSISTENCY,"wrong interference");
  }
  ierr = BaseFunction::query_interface(uuid,iface); CHKERRQ(ierr);
  MoFEMFunctionReturnHot(0);
}

MoFEMErrorCode LegendrePolynomial::getValue(
  MatrixDouble &pts,
  boost::shared_ptr<BaseFunctionCtx> ctx_ptr
) {

  MoFEMFunctionBeginHot;
  MoFEM::UnknownInterface *iface;
  ierr = ctx_ptr->query_interface(IDD_LEGENDRE_BASE_FUNCTION,&iface); CHKERRQ(ierr);
  LegendrePolynomialCtx *ctx = reinterpret_cast<LegendrePolynomialCtx*>(iface);
  ctx->baseFunPtr->resize(pts.size2(),ctx->P+1,false);
  ctx->baseDiffFunPtr->resize(pts.size2(),ctx->dIm*(ctx->P+1),false);
  double *l = NULL;
  double *diff_l = NULL;
  for(unsigned int gg = 0;gg<pts.size2();gg++) {
    if(ctx->baseFunPtr) l = &((*ctx->baseFunPtr)(gg,0));
    if(ctx->baseDiffFunPtr) diff_l = &((*ctx->baseDiffFunPtr)(gg,0));
    ierr = (ctx->basePolynomialsType0)(ctx->P,pts(0,gg),ctx->diffS,l,diff_l,ctx->dIm); CHKERRQ(ierr);
  }
  MoFEMFunctionReturnHot(0);
}

}