/** \file legendrepolynomial.cpp
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
#include <JacobiPolynomial.hpp>

PetscErrorCode JacobiPolynomialCtx::queryInterface(
  const MOFEMuuid& uuid,MoFEM::UnknownInterface** iface
) {
  
  MoFEMFunctionBeginHot;
  *iface = NULL;
  if(uuid == IDD_JACOBI_BASE_FUNCTION) {
    *iface = static_cast<JacobiPolynomialCtx*>(this);
    MoFEMFunctionReturnHot(0);
  } else {
    SETERRQ(PETSC_COMM_WORLD,MOFEM_DATA_INCONSISTENCY,"wrong interference");
  }
  ierr = BaseFunctionCtx::queryInterface(uuid,iface); CHKERRQ(ierr);
  MoFEMFunctionReturnHot(0);
}

PetscErrorCode JacobiPolynomial::queryInterface(
  const MOFEMuuid& uuid,MoFEM::UnknownInterface** iface
) {
  
  MoFEMFunctionBeginHot;
  *iface = NULL;
  if(uuid == IDD_JACOBI_BASE_FUNCTION) {
    *iface = static_cast<JacobiPolynomial*>(this);
    MoFEMFunctionReturnHot(0);
  } else {
    SETERRQ(PETSC_COMM_WORLD,MOFEM_DATA_INCONSISTENCY,"wrong interference");
  }
  ierr = BaseFunction::queryInterface(uuid,iface); CHKERRQ(ierr);
  MoFEMFunctionReturnHot(0);
}

PetscErrorCode JacobiPolynomial::getValue(
  MatrixDouble &pts_x,
  MatrixDouble &pts_t,
  boost::shared_ptr<BaseFunctionCtx> ctx_ptr
) {
  
  MoFEMFunctionBeginHot;
  MoFEM::UnknownInterface *iface;
  ierr = ctx_ptr->queryInterface(IDD_JACOBI_BASE_FUNCTION,&iface); CHKERRQ(ierr);
  JacobiPolynomialCtx *ctx = reinterpret_cast<JacobiPolynomialCtx*>(iface);
  ctx->baseFunPtr->resize(pts_x.size2(),ctx->P+1,false);
  ctx->baseDiffFunPtr->resize(pts_x.size2(),ctx->dIm*(ctx->P+1),false);
  if(
    pts_x.size1()!=pts_t.size1() ||
    pts_x.size2()!=pts_t.size2()
  ) {
    SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"Inconsistent size of arguments");
  }
  double *l = NULL;
  double *diff_l = NULL;
  for(unsigned int gg = 0;gg<pts_x.size2();gg++) {
    if(ctx->baseFunPtr) l = &((*ctx->baseFunPtr)(gg,0));
    if(ctx->baseDiffFunPtr) diff_l = &((*ctx->baseDiffFunPtr)(gg,0));
    ierr = (ctx->basePolynomialsType1)(
      ctx->P,ctx->aLpha,
      pts_x(0,gg),pts_t(0,gg),
      ctx->diffX,ctx->diffT,
      l,diff_l,
      ctx->dIm
    ); CHKERRQ(ierr);
  }
  MoFEMFunctionReturnHot(0);
}
