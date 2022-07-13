/** \file HexPolynomialBase.hpp
\brief Implementation of Ainsworth-Coyle / Demkowicz or any other H1, Hcurl,
Hdiv and L2 base on tetrahedral

*/



#ifndef __HEXPOLYNOMIALBASE_HPP__
#define __HEXPOLYNOMIALBASE_HPP__

namespace MoFEM {

/**
 * \brief Calculate base functions on tetrahedral
 *
 * \ingroup mofem_base_functions
 */
struct HexPolynomialBase : public BaseFunction {

  MoFEMErrorCode query_interface(boost::typeindex::type_index type_index,
                                 UnknownInterface **iface) const;
  HexPolynomialBase() = default;
  ~HexPolynomialBase() = default;

  MoFEMErrorCode getValue(MatrixDouble &pts,
                          boost::shared_ptr<BaseFunctionCtx> ctx_ptr);

private:
  EntPolynomialBaseCtx *cTx;

  /*
   * @brief Get base functions for H1 space
   *
   * @param pts matrix of intergation pts
   * @return MoFEMErrorCode
   *
   * \note matrix of integration points on rows has local coordinates of finite
   * element on columns are integration pts.
   */
  MoFEMErrorCode getValueH1(MatrixDouble &pts);

  /**
   * @brief Get base functions for L2 space
   *
   * @param pts matrix of intergation pts
   * @return MoFEMErrorCode
   *
   * \note matrix of integration points on rows has local coordinates of finite
   * element on columns are integration pts.
   */
  MoFEMErrorCode getValueL2(MatrixDouble &pts);

  /**
   * @brief Get base functions for Hdiv space
   *
   * @param pts matrix of intergation pts
   * @return MoFEMErrorCode
   *
   * \note matrix of integration points on rows has local coordinates of finite
   * element on columns are integration pts.
   */
  MoFEMErrorCode getValueHdiv(MatrixDouble &pts);

  /**
   * @brief Get base functions for Hcurl space
   *
   * @param pts matrix of intergation pts
   * @return MoFEMErrorCode
   *
   * \note matrix of integration points on rows has local coordinates of finite
   * element on columns are integration pts.
   */
  MoFEMErrorCode getValueHcurl(MatrixDouble &pts);

private:
  MoFEMErrorCode getValueH1DemkowiczBase(MatrixDouble &pts);
  MoFEMErrorCode getValueL2DemkowiczBase(MatrixDouble &pts);
  MoFEMErrorCode getValueHdivDemkowiczBase(MatrixDouble &pts);
  MoFEMErrorCode getValueHcurlDemkowiczBase(MatrixDouble &pts);

  std::array<MatrixDouble, 6> faceFamily;
  std::array<MatrixDouble, 6> diffFaceFamily;

  MatrixDouble volFamily;
  MatrixDouble diffVolFamily;
};

} // namespace MoFEM

#endif //__HEXPOLYNOMIALBASE_HPP__
