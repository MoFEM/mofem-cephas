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

namespace MoFEM {

MoFEMErrorCode LobattoPolynomialCtx::query_interface(
    const MOFEMuuid &uuid, BaseFunctionUnknownInterface **iface) const {

  MoFEMFunctionBeginHot;
  *iface = NULL;
  if (uuid == IDD_LOBATTO_BASE_FUNCTION) {
    *iface = const_cast<LobattoPolynomialCtx *>(this);
    MoFEMFunctionReturnHot(0);
  } else {
    SETERRQ(PETSC_COMM_WORLD, MOFEM_DATA_INCONSISTENCY, "wrong interference");
  }
  ierr = LegendrePolynomialCtx::query_interface(uuid, iface);
  CHKERRG(ierr);
  MoFEMFunctionReturnHot(0);
}

MoFEMErrorCode
LobattoPolynomial::query_interface(const MOFEMuuid &uuid,
                                   BaseFunctionUnknownInterface **iface) const {

  MoFEMFunctionBeginHot;
  *iface = NULL;
  if (uuid == IDD_LOBATTO_BASE_FUNCTION) {
    *iface = const_cast<LobattoPolynomial *>(this);
    MoFEMFunctionReturnHot(0);
  } else {
    SETERRQ(PETSC_COMM_WORLD, MOFEM_DATA_INCONSISTENCY, "wrong interference");
  }
  ierr = LegendrePolynomial::query_interface(uuid, iface);
  CHKERRG(ierr);
  MoFEMFunctionReturnHot(0);
}

MoFEMErrorCode
LobattoPolynomial::getValue(MatrixDouble &pts,
                            boost::shared_ptr<BaseFunctionCtx> ctx_ptr) {

  MoFEMFunctionBeginHot;
  BaseFunctionUnknownInterface *iface;
  ierr = ctx_ptr->query_interface(IDD_LOBATTO_BASE_FUNCTION, &iface);
  CHKERRG(ierr);
  LobattoPolynomialCtx *ctx = reinterpret_cast<LobattoPolynomialCtx *>(iface);
  // Polynomial order start from 2nd order
  ctx->baseFunPtr->resize(pts.size2(), ctx->P + 1, false);
  ctx->baseDiffFunPtr->resize(pts.size2(), ctx->dIm * (ctx->P + 1), false);
  double *l = NULL;
  double *diff_l = NULL;
  for (unsigned int gg = 0; gg < pts.size2(); gg++) {
    if (ctx->baseFunPtr)
      l = &((*ctx->baseFunPtr)(gg, 0));
    if (ctx->baseDiffFunPtr)
      diff_l = &((*ctx->baseDiffFunPtr)(gg, 0));
    ierr = (ctx->basePolynomialsType0)(ctx->P, pts(0, gg), ctx->diffS, l,
                                       diff_l, ctx->dIm);
    CHKERRG(ierr);
  }
  MoFEMFunctionReturnHot(0);
}

MoFEMErrorCode KernelLobattoPolynomialCtx::query_interface(
    const MOFEMuuid &uuid, BaseFunctionUnknownInterface **iface) const {

  MoFEMFunctionBeginHot;
  *iface = NULL;
  if (uuid == IDD_KERNEL_BASE_FUNCTION) {
    *iface = const_cast<KernelLobattoPolynomialCtx *>(this);
    MoFEMFunctionReturnHot(0);
  } else {
    SETERRQ(PETSC_COMM_WORLD, MOFEM_DATA_INCONSISTENCY, "wrong interference");
  }
  ierr = LegendrePolynomialCtx::query_interface(uuid, iface);
  CHKERRG(ierr);
  MoFEMFunctionReturnHot(0);
}

MoFEMErrorCode KernelLobattoPolynomial::query_interface(
    const MOFEMuuid &uuid, BaseFunctionUnknownInterface **iface) const {

  MoFEMFunctionBeginHot;
  *iface = NULL;
  if (uuid == IDD_KERNEL_BASE_FUNCTION) {
    *iface = const_cast<KernelLobattoPolynomial *>(this);
    MoFEMFunctionReturnHot(0);
  } else {
    SETERRQ(PETSC_COMM_WORLD, MOFEM_DATA_INCONSISTENCY, "wrong interference");
  }
  ierr = LegendrePolynomial::query_interface(uuid, iface);
  CHKERRG(ierr);
  MoFEMFunctionReturnHot(0);
}

MoFEMErrorCode
KernelLobattoPolynomial::getValue(MatrixDouble &pts,
                                  boost::shared_ptr<BaseFunctionCtx> ctx_ptr) {

  MoFEMFunctionBeginHot;
  BaseFunctionUnknownInterface *iface;
  ierr = ctx_ptr->query_interface(IDD_KERNEL_BASE_FUNCTION, &iface);
  CHKERRG(ierr);
  LobattoPolynomialCtx *ctx = reinterpret_cast<LobattoPolynomialCtx *>(iface);
  ctx->baseFunPtr->resize(pts.size2(), ctx->P + 1, false);
  ctx->baseDiffFunPtr->resize(pts.size2(), ctx->dIm * (ctx->P + 1), false);
  double *l = NULL;
  double *diff_l = NULL;
  for (unsigned int gg = 0; gg < pts.size2(); gg++) {
    if (ctx->baseFunPtr)
      l = &((*ctx->baseFunPtr)(gg, 0));
    if (ctx->baseDiffFunPtr)
      diff_l = &((*ctx->baseDiffFunPtr)(gg, 0));
    ierr = (ctx->basePolynomialsType0)(ctx->P, pts(0, gg), ctx->diffS, l,
                                       diff_l, ctx->dIm);
    CHKERRG(ierr);
  }
  MoFEMFunctionReturnHot(0);
}

} // namespace MoFEM