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

  TetPolynomialBase() = default;
  ~TetPolynomialBase() = default;

  static bool swichCacheHDivBaseFaceDemkowiczCache() {
    cacheHDivBaseFaceDemkowiczCache = !cacheHDivBaseFaceDemkowiczCache;
    return cacheHDivBaseFaceDemkowiczCache;
  };

  static bool swichCacheHdivBaseInteriorDemkowiczCache() {
    cacheHdivBaseInteriorDemkowiczCache = !cacheHdivBaseInteriorDemkowiczCache;
    return cacheHdivBaseInteriorDemkowiczCache;
  };

  MoFEMErrorCode getValue(MatrixDouble &pts,
                          boost::shared_ptr<BaseFunctionCtx> ctx_ptr);

private:
  EntPolynomialBaseCtx *cTx;

  /**
   * @brief Get base functions for H1 space
   *
   * @param pts matrix of integration pts
   * @return MoFEMErrorCode
   *
   * \note matrix of integration points on rows has local coordinates of finite
   * element on columns are integration pts.
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

private:
  MoFEMErrorCode getValueH1AinsworthBase(MatrixDouble &pts);
  MoFEMErrorCode getValueH1BernsteinBezierBase(MatrixDouble &pts);

  MoFEMErrorCode getValueL2AinsworthBase(MatrixDouble &pts);
  MoFEMErrorCode getValueL2BernsteinBezierBase(MatrixDouble &pts);

  MoFEMErrorCode getValueHdivAinsworthBase(MatrixDouble &pts);
  MoFEMErrorCode getValueHcurlAinsworthBase(MatrixDouble &pts);

  MoFEMErrorCode getValueHdivDemkowiczBase(MatrixDouble &pts);
  MoFEMErrorCode getValueHcurlDemkowiczBase(MatrixDouble &pts);

  MatrixInt senseFaceAlpha;

  struct BaseCacheItem {
    int order;
    int nb_gauss_pts;
    mutable MatrixDouble N;
    mutable MatrixDouble diffN;
  };

  using BaseCacheMI = boost::multi_index_container<
      BaseCacheItem,
      boost::multi_index::indexed_by<

          boost::multi_index::hashed_unique<

              composite_key<

                  BaseCacheItem,
                  member<BaseCacheItem, int, &BaseCacheItem::order>,
                  member<BaseCacheItem, int, &BaseCacheItem::nb_gauss_pts>>>>

      >;

  struct HDivBaseCacheItem {

    int order;
    int nb_gauss_pts;

    // Number of permeations for tetrahedron
    // That is P(3, 4) = 24
    
    int n0;
    int n1;
    int n2;

    mutable MatrixDouble N;
    mutable MatrixDouble diffN;
  };

  using HDivBaseFaceCacheMI = boost::multi_index_container<
      HDivBaseCacheItem,
      boost::multi_index::indexed_by<

          boost::multi_index::hashed_unique<

              composite_key<

                  HDivBaseCacheItem,
                  member<HDivBaseCacheItem, int, &HDivBaseCacheItem::order>,
                  member<HDivBaseCacheItem, int,
                         &HDivBaseCacheItem::nb_gauss_pts>,
                  member<HDivBaseCacheItem, int, &HDivBaseCacheItem::n0>,
                  member<HDivBaseCacheItem, int, &HDivBaseCacheItem::n1>,
                  member<HDivBaseCacheItem, int, &HDivBaseCacheItem::n2>>>>

      >;

  static bool cacheHDivBaseFaceDemkowiczCache;
  HDivBaseFaceCacheMI hDivBaseDemkowiczCache;
  static bool cacheHdivBaseInteriorDemkowiczCache;
  BaseCacheMI hdivBaseInteriorDemkowiczCache;

};

} // namespace MoFEM

#endif //__TETPOLYNOMIALBASE_HPP__
