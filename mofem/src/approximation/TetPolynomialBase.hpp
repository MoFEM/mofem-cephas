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
  static bool switchCacheBaseFace(FieldApproximationBase base, void *ptr);

  template <int SPACE>
  static bool switchCacheBaseInterior(FieldApproximationBase base, void *ptr);

  template <int SPACE>
  static bool switchCacheBrokenBaseInterior(FieldApproximationBase base,
                                           void *ptr);

  template <int SPACE>
  static void switchCacheBaseOn(FieldApproximationBase base,
                               std::vector<void *> v);

  template <int SPACE>
  static void switchCacheBaseOff(FieldApproximationBase base,
                                std::vector<void *> v);

  template <int SPACE> static void switchCacheBaseOn(std::vector<void *> v);

  template <int SPACE> static void switchCacheBaseOff(std::vector<void *> v);

  MoFEMErrorCode getValue(MatrixDouble &pts,
                          boost::shared_ptr<BaseFunctionCtx> ctx_ptr);

  /**
   * @brief Set map of dof to side number
   *
   * That is used for broken space to establish connection between dofs in the
   * interior of element/entity and side of element/entity to which that dof is
   * associated. That depends on implementation of the base for given space, and
   * has to be implemented while implementing base function for given space.
   *
   * @param space
   * @param continuity
   * @param base
   * @param DofsSideMap
   * @return MoFEMErrorCode
   */
  static MoFEMErrorCode
  setDofsSideMap(const FieldSpace space, const FieldContinuity continuity,
                 const FieldApproximationBase base,
                 DofsSideMap &);

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

  /**
   * @brief Set the Dofs Side Map Hdiv object
   *
   * @param space
   * @param continuity
   * @param base
   * @param dofs_side_map
   * @return MoFEMErrorCode
   */
  static MoFEMErrorCode setDofsSideMapHdiv(const FieldContinuity continuity,
                                           const FieldApproximationBase base,
                                           DofsSideMap &dofs_side_map);

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
bool TetPolynomialBase::switchCacheBaseFace<HDIV>(FieldApproximationBase base,
                                                 void *ptr);
template <>
bool TetPolynomialBase::switchCacheBaseInterior<HDIV>(
    FieldApproximationBase base, void *ptr);

template <>
bool TetPolynomialBase::switchCacheBrokenBaseInterior<HDIV>(
    FieldApproximationBase base, void *ptr);

template <>
void TetPolynomialBase::switchCacheBaseOn<HDIV>(FieldApproximationBase base,
                                               std::vector<void *> v);

template <>
void TetPolynomialBase::switchCacheBaseOff<HDIV>(FieldApproximationBase base,
                                                std::vector<void *> v);

template <>
void TetPolynomialBase::switchCacheBaseOn<HDIV>(std::vector<void *> v);

template <>
void TetPolynomialBase::switchCacheBaseOff<HDIV>(std::vector<void *> v);

template <>
bool TetPolynomialBase::switchCacheBaseInterior<L2>(FieldApproximationBase base,
                                                   void *ptr);

template <>
void TetPolynomialBase::switchCacheBaseOn<L2>(FieldApproximationBase base,
                                               std::vector<void *> v);

template <>
void TetPolynomialBase::switchCacheBaseOff<L2>(FieldApproximationBase base,
                                              std::vector<void *> v);

template <>
void TetPolynomialBase::switchCacheBaseOn<L2>(std::vector<void *> v);

template <>
void TetPolynomialBase::switchCacheBaseOff<L2>(std::vector<void *> v);


} // namespace MoFEM

#endif //__TETPOLYNOMIALBASE_HPP__
