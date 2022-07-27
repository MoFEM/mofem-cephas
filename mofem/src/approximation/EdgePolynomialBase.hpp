/** \file EdgePolynomialBase.hpp
\brief Implementation of base on tetrahedral for H1 bases.

TODO:
\todo L2 base on edge

*/

#ifndef __EDGEPOLYNOMIALBASE_HPP__
#define __EDGEPOLYNOMIALBASE_HPP__

namespace MoFEM {

/**
 * \brief Calculate base functions on tetrahedral
 *
 * \ingroup mofem_base_functions
 */
struct EdgePolynomialBase : public BaseFunction {

  MoFEMErrorCode query_interface(boost::typeindex::type_index type_index,
                                 UnknownInterface **iface) const;

  EdgePolynomialBase();
  ~EdgePolynomialBase();

  MoFEMErrorCode getValue(MatrixDouble &pts,
                          boost::shared_ptr<BaseFunctionCtx> ctx_ptr);

private:
  EntPolynomialBaseCtx *cTx;

  VectorDouble L, diffL;

  MoFEMErrorCode getValueH1(MatrixDouble &pts);

  MoFEMErrorCode getValueH1AinsworthBase(MatrixDouble &pts);

  MoFEMErrorCode getValueH1DemkowiczBase(MatrixDouble &pts);

  MoFEMErrorCode getValueH1BernsteinBezierBase(MatrixDouble &pts);

  MoFEMErrorCode getValueL2(MatrixDouble &pts);

  MoFEMErrorCode getValueL2DemkowiczBase(MatrixDouble &pts);

  MoFEMErrorCode getValueHdiv(MatrixDouble &pts);

  MoFEMErrorCode getValueHcurl(MatrixDouble &pts);

  MoFEMErrorCode getValueHcurlAinsworthBase(MatrixDouble &pts);

  MoFEMErrorCode getValueHcurlDemkowiczBase(MatrixDouble &pts);
};

} // namespace MoFEM

#endif //__EDGEPOLYNOMIALBASE_HPP__
