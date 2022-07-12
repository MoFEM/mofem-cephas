/** \file legendrepolynomial.cpp
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
JacobiPolynomialCtx::query_interface(boost::typeindex::type_index type_index,
                                     UnknownInterface **iface) const {
  *iface = const_cast<JacobiPolynomialCtx *>(this);
  return 0;
}

MoFEMErrorCode
JacobiPolynomial::query_interface(boost::typeindex::type_index type_index,
                                  UnknownInterface **iface) const {
  MoFEMFunctionBeginHot;
  *iface = const_cast<JacobiPolynomial *>(this);
  MoFEMFunctionReturnHot(0);
}

template <class TYPE>
static MoFEMErrorCode get_value(MatrixDouble &pts_x, MatrixDouble &pts_t,
                                TYPE *ctx) {
  MoFEMFunctionBeginHot;
  ctx->baseFunPtr->resize(pts_x.size2(), ctx->P + 1, false);
  ctx->baseDiffFunPtr->resize(pts_x.size2(), ctx->dIm * (ctx->P + 1), false);
  if (pts_x.size1() != pts_t.size1() || pts_x.size2() != pts_t.size2()) {
    SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
            "Inconsistent size of arguments");
  }
  double *l = NULL;
  double *diff_l = NULL;
  for (unsigned int gg = 0; gg < pts_x.size2(); gg++) {
    if (ctx->baseFunPtr)
      l = &((*ctx->baseFunPtr)(gg, 0));
    if (ctx->baseDiffFunPtr)
      diff_l = &((*ctx->baseDiffFunPtr)(gg, 0));
    ierr = (ctx->basePolynomialsType1)(ctx->P, ctx->aLpha, pts_x(0, gg),
                                       pts_t(0, gg), ctx->diffX, ctx->diffT, l,
                                       diff_l, ctx->dIm);
    CHKERRG(ierr);
  }
  MoFEMFunctionReturnHot(0);
}

MoFEMErrorCode
JacobiPolynomial::getValue(MatrixDouble &pts_x, MatrixDouble &pts_t,
                           boost::shared_ptr<BaseFunctionCtx> ctx_ptr) {
  MoFEMFunctionBegin;
  auto ctx = ctx_ptr->getInterface<JacobiPolynomialCtx>();
  CHKERR get_value(pts_x, pts_t, ctx);
  MoFEMFunctionReturn(0);
}

MoFEMErrorCode IntegratedJacobiPolynomialCtx::query_interface(
    boost::typeindex::type_index type_index, UnknownInterface **iface) const {
  *iface = const_cast<IntegratedJacobiPolynomialCtx *>(this);
  return 0;
}

MoFEMErrorCode IntegratedJacobiPolynomial::query_interface(
    boost::typeindex::type_index type_index, UnknownInterface **iface) const {
  *iface = const_cast<IntegratedJacobiPolynomial *>(this);
  return 0;
}

MoFEMErrorCode IntegratedJacobiPolynomial::getValue(
    MatrixDouble &pts_x, MatrixDouble &pts_t,
    boost::shared_ptr<BaseFunctionCtx> ctx_ptr) {
  MoFEMFunctionBegin;
  auto ctx = ctx_ptr->getInterface<IntegratedJacobiPolynomialCtx>();
  CHKERR get_value(pts_x, pts_t, ctx);
  MoFEMFunctionReturn(0);
}

} // namespace MoFEM