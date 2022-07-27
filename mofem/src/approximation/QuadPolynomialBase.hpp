/** \file QuadPolynomialBase.hpp
\brief Implementation of H1 base on a quad face

\todo Quad element can be integrated exploiting tonsorial product. Current
implementation do not take that opportunity. That can be viewed as a bug. 

*/



#ifndef __H1QUADPOLYNOMIAL_HPP__
#define __H1QUADPOLYNOMIAL_HPP__

namespace MoFEM {

/**
 * \brief Calculate base functions on triangle
 *
 * \ingroup mofem_base_functions
 */
struct QuadPolynomialBase : public BaseFunction {

  MoFEMErrorCode query_interface(boost::typeindex::type_index type_index, UnknownInterface **iface) const;

  QuadPolynomialBase() = default;
  ~QuadPolynomialBase() = default;

  MoFEMErrorCode getValue(MatrixDouble &pts,
                          boost::shared_ptr<BaseFunctionCtx> ctx_ptr);

private:
  EntPolynomialBaseCtx *cTx;

  MoFEMErrorCode getValueH1(MatrixDouble &pts);
  MoFEMErrorCode getValueL2(MatrixDouble &pts);
  MoFEMErrorCode getValueHcurl(MatrixDouble &pts);
  MoFEMErrorCode getValueHdiv(MatrixDouble &pts);

  MoFEMErrorCode getValueH1AinsworthBase(MatrixDouble &pts);
  MoFEMErrorCode getValueH1DemkowiczBase(MatrixDouble &pts);
  MoFEMErrorCode getValueL2DemkowiczBase(MatrixDouble &pts);
  MoFEMErrorCode getValueHcurlDemkowiczBase(MatrixDouble &pts);
  MoFEMErrorCode getValueHdivDemkowiczBase(MatrixDouble &pts);

  MatrixDouble faceFamily;
  MatrixDouble diffFaceFamily;

};

} // namespace MoFEM

#endif //__H1QUADPOLYNOMIAL_HPP__
