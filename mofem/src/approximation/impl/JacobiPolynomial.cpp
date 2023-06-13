/** \file legendrepolynomial.cpp
 * \brief implementation of multi-grid solver for p- adaptivity
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