/** \file LobattoPolynomial.cpp
 * \brief implementation of multi-grid solver for p- adaptivity
 */

/* MIT License
 *
 * Copyright (c) 2022
 * 
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 * 
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 * 
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */

namespace MoFEM {

MoFEMErrorCode
LobattoPolynomialCtx::query_interface(boost::typeindex::type_index type_index,
                                      UnknownInterface **iface) const {
  *iface = const_cast<LobattoPolynomialCtx *>(this);
  return 0;
}

MoFEMErrorCode
LobattoPolynomial::query_interface(boost::typeindex::type_index type_index,
                                   UnknownInterface **iface) const {
  *iface = const_cast<LobattoPolynomial *>(this);
  return 0;
}

MoFEMErrorCode
LobattoPolynomial::getValue(MatrixDouble &pts,
                            boost::shared_ptr<BaseFunctionCtx> ctx_ptr) {

  MoFEMFunctionBeginHot;
  auto ctx = ctx_ptr->getInterface<LobattoPolynomialCtx>();
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
    boost::typeindex::type_index type_index, UnknownInterface **iface) const {
  *iface = const_cast<KernelLobattoPolynomialCtx *>(this);
  return 0;
}

MoFEMErrorCode KernelLobattoPolynomial::query_interface(
    boost::typeindex::type_index type_index, UnknownInterface **iface) const {
  *iface = const_cast<KernelLobattoPolynomial *>(this);
  return 0;
}

MoFEMErrorCode
KernelLobattoPolynomial::getValue(MatrixDouble &pts,
                                  boost::shared_ptr<BaseFunctionCtx> ctx_ptr) {

  MoFEMFunctionBeginHot;
  BaseFunctionUnknownInterface *iface;
  auto ctx = ctx_ptr->getInterface<KernelLobattoPolynomialCtx>();
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