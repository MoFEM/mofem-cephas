/** \file BaseFunction.cpp
 * \brief implementation of multi-grid solver for p- adaptivity
 */



namespace MoFEM {

MoFEMErrorCode
BaseFunctionCtx::query_interface(boost::typeindex::type_index type_index,
                                 UnknownInterface **iface) const {
  *iface = const_cast<BaseFunctionCtx *>(this);
  return 0;
}

MoFEMErrorCode
BaseFunction::query_interface(boost::typeindex::type_index type_index,
                              UnknownInterface **iface) const {
  *iface = const_cast<BaseFunction *>(this);
  return 0;
}

MoFEMErrorCode
BaseFunction::getValue(MatrixDouble &pts,
                       boost::shared_ptr<BaseFunctionCtx> ctx_ptr) {
  MoFEMFunctionBeginHot;
  SETERRQ(PETSC_COMM_SELF, MOFEM_NOT_IMPLEMENTED,
          "BaseFunction has not valid implementation of any shape function");
  MoFEMFunctionReturnHot(0);
}

MoFEMErrorCode
BaseFunction::getValue(MatrixDouble &pts_x, MatrixDouble &pts_t,
                       boost::shared_ptr<BaseFunctionCtx> ctx_ptr) {
  MoFEMFunctionBeginHot;
  SETERRQ(PETSC_COMM_SELF, MOFEM_NOT_IMPLEMENTED,
          "BaseFunction has not valid implementation of any shape function");
  MoFEMFunctionReturnHot(0);
}

} // namespace MoFEM