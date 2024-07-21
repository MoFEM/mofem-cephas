/** \file TetPolynomialBase.hpp
\brief Implementation of Ainsworth-Coyle / Demkowicz or any other H1, Hcurl,
Hdiv and L2 base on tetrahedral

*/

#ifndef __TETPOLYNOMIALBASE_HPP__
#define __TETPOLYNOMIALBASE_HPP__

namespace MoFEM {

/**
 * \brief Calculate base functions on tetrahedral
 *
 * \ingroup mofem_base_functions
 */
struct TetPolynomialBase : public BaseFunction {

  MoFEMErrorCode query_interface(boost::typeindex::type_index type_index,
                                 UnknownInterface **iface) const;

  TetPolynomialBase(const void *ptr = nullptr);
  virtual ~TetPolynomialBase();

  template <int SPACE>
  static bool swichCacheBaseFace(FieldApproximationBase base, void *ptr);

  template <int SPACE>
  static bool swichCacheBaseInterior(FieldApproximationBase base, void *ptr);

  template <int SPACE>
  static bool swichCacheBrokenBaseInterior(FieldApproximationBase base,
                                           void *ptr);

  template <int SPACE>
  static void swichCacheBaseOn(FieldApproximationBase base,
                               std::vector<void *> v);

  template <int SPACE>
  static void swichCacheBaseOff(FieldApproximationBase base,
                                std::vector<void *> v);

  template <int SPACE> static void swichCacheBaseOn(std::vector<void *> v);

  template <int SPACE> static void swichCacheBaseOff(std::vector<void *> v);

  MoFEMErrorCode getValue(MatrixDouble &pts,
                          boost::shared_ptr<BaseFunctionCtx> ctx_ptr);

private:
  const void *vPtr;
  EntPolynomialBaseCtx *cTx;

  /**
   * @brief Get base functions for H1 space
   *
   * @param pts matrix of integration pts
   * @return MoFEMErrorCode
   *
   * \note matrix of integration points on rows has local coordinates of
   * finite element on columns are integration pts.
   */
  MoFEMErrorCode getValueH1(MatrixDouble &pts);

  /**
   * @brief Get base functions for L2 space
   *
   * @param pts matrix of integration pts
   * @return MoFEMErrorCode
   *
   * \note matrix of integration points on rows has local coordinates of finite
   * element on columns are integration pts.
   */
  MoFEMErrorCode getValueL2(MatrixDouble &pts);

  /**
   * @brief Get base functions for Hdiv space
   *
   * @param pts matrix of integration pts
   * @return MoFEMErrorCode
   *
   * \note matrix of integration points on rows has local coordinates of finite
   * element on columns are integration pts.
   */
  MoFEMErrorCode getValueHdiv(MatrixDouble &pts);

  /**
   * @brief Get base functions for Hcurl space
   *
   * @param pts matrix of integration pts
   * @return MoFEMErrorCode
   *
   * \note matrix of integration points on rows has local coordinates of finite
   * element on columns are integration pts.
   */
  MoFEMErrorCode getValueHcurl(MatrixDouble &pts);

private:
  MoFEMErrorCode getValueH1AinsworthBase(MatrixDouble &pts);
  MoFEMErrorCode getValueH1BernsteinBezierBase(MatrixDouble &pts);

  MoFEMErrorCode getValueL2AinsworthBase(MatrixDouble &pts);
  MoFEMErrorCode getValueL2BernsteinBezierBase(MatrixDouble &pts);

  MoFEMErrorCode getValueHdivAinsworthBaseImpl(

      MatrixDouble &pts,

      MatrixDouble &shape_functions, MatrixDouble &diff_shape_functions,

      int volume_order, std::array<int, 4> &faces_order,
      std::array<int, 3 * 4> &faces_nodes

  );
  MoFEMErrorCode getValueHdivAinsworthBase(MatrixDouble &pts);
  MoFEMErrorCode getValueHdivAinsworthBrokenBase(MatrixDouble &pts);
  MoFEMErrorCode getValueHcurlAinsworthBase(MatrixDouble &pts);

  MoFEMErrorCode getValueHdivDemkowiczBase(MatrixDouble &pts);
  MoFEMErrorCode getValueHdivDemkowiczBrokenBase(MatrixDouble &pts);
  MoFEMErrorCode getValueHcurlDemkowiczBase(MatrixDouble &pts);

  MatrixInt senseFaceAlpha;

  ublas::matrix<MatrixDouble> N_face_edge;
  ublas::vector<MatrixDouble> N_face_bubble;
  ublas::vector<MatrixDouble> N_volume_edge;
  ublas::vector<MatrixDouble> N_volume_face;
  MatrixDouble N_volume_bubble;

  ublas::matrix<MatrixDouble> diffN_face_edge;
  ublas::vector<MatrixDouble> diffN_face_bubble;
  ublas::vector<MatrixDouble> diffN_volume_edge;
  ublas::vector<MatrixDouble> diffN_volume_face;
  MatrixDouble diffN_volume_bubble;
};

template <>
bool TetPolynomialBase::swichCacheBaseFace<HDIV>(FieldApproximationBase base,
                                                 void *ptr);
template <>
bool TetPolynomialBase::swichCacheBaseInterior<HDIV>(
    FieldApproximationBase base, void *ptr);

template <>
bool TetPolynomialBase::swichCacheBrokenBaseInterior<HDIV>(
    FieldApproximationBase base, void *ptr);

template <>
void TetPolynomialBase::swichCacheBaseOn<HDIV>(FieldApproximationBase base,
                                               std::vector<void *> v);

template <>
void TetPolynomialBase::swichCacheBaseOff<HDIV>(FieldApproximationBase base,
                                                std::vector<void *> v);

template <>
void TetPolynomialBase::swichCacheBaseOn<HDIV>(std::vector<void *> v);

template <>
void TetPolynomialBase::swichCacheBaseOff<HDIV>(std::vector<void *> v);

} // namespace MoFEM

#endif //__TETPOLYNOMIALBASE_HPP__
