/** \file Tools.hpp
 * \brief Tools interface
 *
 * Implementation of some useful and very often used methods * in MoFEM.
 *
 * \ingroup mofem_tools
 */

/* MoFEM is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
 * License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with MoFEM. If not, see <http://www.gnu.org/licenses/>
 */

#ifndef __TOOLS_HPP__
#define __TOOLS_HPP__

namespace MoFEM {

static const MOFEMuuid IDD_MOFEMTools = MOFEMuuid(BitIntefaceId(TOOLS));

/**
 * \brief Auxiliary tools
 * \nosubgrouping
 * \ingroup mofem_tools
 */
struct Tools : public UnknownInterface {

  MoFEMErrorCode query_interface(const MOFEMuuid &uuid,
                                 UnknownInterface **iface) const;

  MoFEM::Core &cOre;
  Tools(const MoFEM::Core &core) : cOre(const_cast<MoFEM::Core &>(core)) {}

  /** \name Computational */

  /**
   * \brief Calculate tetrahedron volume length quality
   * @param  coords tet coordinates
   * @return        Volume-length quality
   */
  static double volumeLengthQuality(const double *coords);

  /**
   * @brief Calculate volume of tetrahedron
   *
   * @param coords
   * @return double volume
   */
  static double tetVolume(const double *coords);

  /**
   * \brief calculate minimal quality of tetrahedra in range
   * @param  tets        range
   * @param  min_quality mimimal quality
   * @return             error code
   */
  MoFEMErrorCode minTetsQuality(const Range &tets, double &min_quality,
                                Tag th = nullptr,
                                boost::function<double(double, double)> f =
                                    [](double a, double b) -> double {
                                  return std::min(a, b);
                                });

  static constexpr double diffN_MBEDGE0x = diffN_MBEDGE0;
  static constexpr double diffN_MBEDGE1x = diffN_MBEDGE1;

  static constexpr std::array<double, 2> diffShapeFunMBEDGE = {diffN_MBEDGE0x,
                                                               diffN_MBEDGE1x};

  static constexpr double diffShapeFunMBTRI0x =
      diffN_MBTRI0x; ///< derivative of triangle shape function
  static constexpr double diffShapeFunMBTRI0y =
      diffN_MBTRI0y; ///< derivative of triangle shape function
  static constexpr double diffShapeFunMBTRI1x =
      diffN_MBTRI1x; ///< derivative of triangle shape function
  static constexpr double diffShapeFunMBTRI1y =
      diffN_MBTRI1y; ///< derivative of triangle shape function
  static constexpr double diffShapeFunMBTRI2x =
      diffN_MBTRI2x; ///< derivative of triangle shape function
  static constexpr double diffShapeFunMBTRI2y =
      diffN_MBTRI2y; ///< derivative of triangle shape function

  static constexpr std::array<double, 6> diffShapeFunMBTRI = {

      diffShapeFunMBTRI0x, diffShapeFunMBTRI0y,

      diffShapeFunMBTRI1x, diffShapeFunMBTRI1y,

      diffShapeFunMBTRI2x, diffShapeFunMBTRI2y};

  static inline double shapeFunMBTRI0(const double x, const double y);

  static inline double shapeFunMBTRI1(const double x, const double y);

  static inline double shapeFunMBTRI2(const double x, const double y);

  /**
   * @brief Calculate shape functions on triangle
   *
   * \note Template parameter is leading dimension of point coordinate arrays,
   * such that \f$ksi_{n+1} = ksi[n + LDB]\f$
   *
   * @tparam 1
   * @param shape shape functions
   * @param ksi pointer to first local coordinates
   * @param eta pointer to second local coordinates
   * @param nb number of points
   * @return MoFEMErrorCode
   */
  template <int LDB = 1>
  static MoFEMErrorCode shapeFunMBTRI(double *shape, const double *ksi,
                                      const double *eta, const int nb);

  static constexpr double diffShapeFunMBQUADAtCenter0x =
      diffN_MBQUAD0x(0.5); ///< derivative of quad shape function
  static constexpr double diffShapeFunMBQUADAtCenter0y =
      diffN_MBQUAD0y(0.5); ///< derivative of quad shape function
  static constexpr double diffShapeFunMBQUADAtCenter1x =
      diffN_MBQUAD1x(0.5); ///< derivative of quad shape function
  static constexpr double diffShapeFunMBQUADAtCenter1y =
      diffN_MBQUAD1y(0.5); ///< derivative of quad shape function
  static constexpr double diffShapeFunMBQUADAtCenter2x =
      diffN_MBQUAD2x(0.5); ///< derivative of quad shape function
  static constexpr double diffShapeFunMBQUADAtCenter2y =
      diffN_MBQUAD2y(0.5); ///< derivative of quad shape function
  static constexpr double diffShapeFunMBQUADAtCenter3x =
      diffN_MBQUAD3x(0.5); ///< derivative of quad shape function
  static constexpr double diffShapeFunMBQUADAtCenter3y =
      diffN_MBQUAD3y(0.5); ///< derivative of quad shape function

  static constexpr std::array<double, 8> diffShapeFunMBQUADAtCenter = {
      diffShapeFunMBQUADAtCenter0x, diffShapeFunMBQUADAtCenter0y,
      diffShapeFunMBQUADAtCenter1x, diffShapeFunMBQUADAtCenter1y,
      diffShapeFunMBQUADAtCenter2x, diffShapeFunMBQUADAtCenter2y,
      diffShapeFunMBQUADAtCenter3x, diffShapeFunMBQUADAtCenter3y};

  static constexpr double diffShapeFunMBTET0x =
      diffN_MBTET0x; ///< derivative of tetrahedral shape function
  static constexpr double diffShapeFunMBTET0y =
      diffN_MBTET0y; ///< derivative of tetrahedral shape function
  static constexpr double diffShapeFunMBTET0z =
      diffN_MBTET0z; ///< derivative of tetrahedral shape function
  static constexpr double diffShapeFunMBTET1x =
      diffN_MBTET1x; ///< derivative of tetrahedral shape function
  static constexpr double diffShapeFunMBTET1y =
      diffN_MBTET1y; ///< derivative of tetrahedral shape function
  static constexpr double diffShapeFunMBTET1z =
      diffN_MBTET1z; ///< derivative of tetrahedral shape function
  static constexpr double diffShapeFunMBTET2x =
      diffN_MBTET2x; ///< derivative of tetrahedral shape function
  static constexpr double diffShapeFunMBTET2y =
      diffN_MBTET2y; ///< derivative of tetrahedral shape function
  static constexpr double diffShapeFunMBTET2z =
      diffN_MBTET2z; ///< derivative of tetrahedral shape function
  static constexpr double diffShapeFunMBTET3x =
      diffN_MBTET3x; ///< derivative of tetrahedral shape function
  static constexpr double diffShapeFunMBTET3y =
      diffN_MBTET3y; ///< derivative of tetrahedral shape function
  static constexpr double diffShapeFunMBTET3z =
      diffN_MBTET3z; ///< derivative of tetrahedral shape function

  static constexpr std::array<double, 12> diffShapeFunMBTET = {

      diffShapeFunMBTET0x, diffShapeFunMBTET0y, diffShapeFunMBTET0z,

      diffShapeFunMBTET1x, diffShapeFunMBTET1y, diffShapeFunMBTET1z,

      diffShapeFunMBTET2x, diffShapeFunMBTET2y, diffShapeFunMBTET2z,

      diffShapeFunMBTET3x, diffShapeFunMBTET3y, diffShapeFunMBTET3z};

  static inline double shapeFunMBTET0(const double x, const double y,
                                      const double z);

  static inline double shapeFunMBTET1(const double x, const double y,
                                      const double z);

  static inline double shapeFunMBTET2(const double x, const double y,
                                      const double z);

  static inline double shapeFunMBTET3(const double x, const double y,
                                      const double z);

  static constexpr double shapeFunMBTET0At000 = N_MBTET0(0, 0, 0);
  static constexpr double shapeFunMBTET1At000 = N_MBTET1(0, 0, 0);
  static constexpr double shapeFunMBTET2At000 = N_MBTET2(0, 0, 0);
  static constexpr double shapeFunMBTET3At000 = N_MBTET3(0, 0, 0);

  static constexpr double shapeFunMBTET0AtOneThird =
      N_MBTET0(1. / 3., 1. / 3., 1. / 3.);
  static constexpr double shapeFunMBTET1AtOneThird =
      N_MBTET1(1. / 3., 1. / 3., 1. / 3.);
  static constexpr double shapeFunMBTET2AtOneThird =
      N_MBTET2(1. / 3., 1. / 3., 1. / 3.);
  static constexpr double shapeFunMBTET3AtOneThird =
      N_MBTET3(1. / 3., 1. / 3., 1. / 3.);

  /**
   * @brief Calculate shape functions on tetrahedron
   *
   * \note Template parameter is leading dimension of point coordinate arrays,
   * such that \f$ksi_{n+1} = ksi[n + LDB]\f$
   *
   * @tparam 1
   * @param shape shape functions
   * @param ksi pointer to first local coordinates
   * @param eta pointer to second local coordinates
   * @param zeta pointer to first third coordinates
   * @param nb number of points
   * @return MoFEMErrorCode
   */
  template <int LDB = 1>
  static MoFEMErrorCode shapeFunMBTET(double *shape, const double *ksi,
                                      const double *eta, const double *zeta,
                                      const double nb);

  /**
   * @brief Array of shape function at zero local point on reference element
   *
   */
  static constexpr std::array<double, 4> shapeFunMBTETAt000 = {
      shapeFunMBTET0At000, shapeFunMBTET1At000, shapeFunMBTET2At000,
      shapeFunMBTET3At000};

  /**
   * @brief Array of shape function at center on reference element
   *
   */
  static constexpr std::array<double, 4> shapeFunMBTETAtOneThird = {
      shapeFunMBTET0AtOneThird, shapeFunMBTET1AtOneThird,
      shapeFunMBTET2AtOneThird, shapeFunMBTET3AtOneThird};


  /**
   * @brief Get the Local Coordinates On Reference Four Node Tet object
   *
   * \code
   * MatrixDouble elem_coords(4, 3);
   * // Set nodal coordinates
   * MatrixDouble global_coords(5, 3);
   * // Set global coordinates
   * MatrixDouble local_coords(global_coords.size1(), 3);
   * CHKERR Tools::getLocalCoordinatesOnReferenceFourNodeTet(
   *     &elem_coords(0, 0), &global_coords(0, 0), global_coords.size1(),
   *     &local_coords(0, 0))
   * \endcode
   *
   * @param elem_coords Global element node coordinates
   * @param glob_coords Globale coordinates
   * @param nb_nodes Number of points
   * @param local_coords Result
   * @return MoFEMErrorCode
   */
  static MoFEMErrorCode getLocalCoordinatesOnReferenceFourNodeTet(
      const double *elem_coords, const double *glob_coords, const int nb_nodes,
      double *local_coords);

  /**
   * @brief Get the Tets With Quality
   *
   * @param out_tets
   * @param tets
   * @param th
   * @param f
   * @return MoFEMErrorCode
   */
  MoFEMErrorCode
  getTetsWithQuality(Range &out_tets, const Range &tets, Tag th = nullptr,
                     boost::function<bool(double)> f = [](double q) -> bool {
                       if (q <= 0)
                         return true;
                       else
                         return false;
                     });
  /**
   * @brief Write file with tetrahedral of given quality
   *
   * @param file_name
   * @param file_type
   * @param options
   * @param tets
   * @param th
   * @param f
   * @return MoFEMErrorCode
   */
  MoFEMErrorCode
  writeTetsWithQuality(const char *file_name, const char *file_type,
                       const char *options, const Range &tets, Tag th = nullptr,
                       boost::function<bool(double)> f = [](double q) -> bool {
                         if (q <= 0)
                           return true;
                         else
                           return false;
                       });

  /**
   * @brief Check of point is in tetrahedral
   *
   * @param tet_coords
   * @param global_coord
   * @param tol
   * @param result
   * @return MoFEMErrorCode
   */
  static MoFEMErrorCode checkIfPointIsInTet(const double tet_coords[],
                                            const double global_coord[],
                                            const double tol, bool &result);

  /**
   * @brief Get the Tri Normal objectGet triangle normal
   *
   * @param coords
   * @param normal
   * @return MoFEMErrorCode
   */
  static MoFEMErrorCode getTriNormal(const double *coords, double *normal);

  /**
   * @brief Get triangle normal
   *
   * @param tri
   * @param normal
   * @return MoFEMErrorCode
   */
  MoFEMErrorCode getTriNormal(const EntityHandle tri, double *normal) const;

  /**
   * @brief Get triangle area
   *
   * @param tri
   * @return double
   */
  double getTriArea(const EntityHandle tri) const;

  /**
   * @brief Get edge length
   *
   * @param edge_coords
   * @return double
   */
  static double getEdgeLength(const double *edge_coords);

  /**
   * @brief Get edge length
   *
   * @param edge
   * @return double
   */
  double getEdgeLength(const EntityHandle edge);

  enum SEGMENT_MIN_DISTANCE {
    SOLUTION_EXIST,
    SEGMENT_ONE_IS_POINT,
    SEGMENT_TWO_IS_POINT,
    SEGMENT_TWO_AND_TWO_ARE_POINT,
    NO_SOLUTION
  };

  /**
   * @brief Find closet point on the segment from the point
   *
   * @param w_ptr segment first vertex coordinate
   * @param v_ptr segment second vertex coordinate
   * @param p_ptr coordinate of point
   * @param t_ptr distance on the segment
   *
   * \note If t is outside bounds [ 0,-1 ] point is on the line point
   * beyond segment.
   *
   * \code
   *
   * double w[] = {-1, 0, 0};
   * double v[] = {1, 0, 0};
   * double p[] = {0, 1, 0};
   * double t;
   * CHKERR Toolas::minDistancePointFromOnSegment(w, v, p, &t);
   * double point_on_segment[3];
   * for (int i = 0; i != 3; ++i)
   *   point_on_segment[i] = w[i] + t * (v[i] - w[i]);
   *
   * \endcode
   *
   * @return SEGMENT_MIN_DISTANCE
   */
  static SEGMENT_MIN_DISTANCE
  minDistancePointFromOnSegment(const double *w_ptr, const double *v_ptr,
                                const double *p_ptr,
                                double *const t_ptr = nullptr);
  /**
   * @brief Find points on two segments in closest distance
   *
   * @param w_ptr
   * @param v_ptr
   * @param k_ptr
   * @param l_ptr
   * @param tvw_ptr
   * @param tlk_ptr
   * @return SEGMENT_MIN_DISTANCE
   *
   * \note If tvwk or tlk are outside bound [0,-1], it means that points
   * are on the lines beyond segments, respectively for segment vw and
   * lk.
   *
   */
  static SEGMENT_MIN_DISTANCE
  minDistanceFromSegments(const double *w_ptr, const double *v_ptr,
                          const double *k_ptr, const double *l_ptr,
                          double *const tvw_ptr = nullptr,
                          double *const tlk_ptr = nullptr);

  /**
   * @brief Find minimal distance to edges
   *
   * \note Finding only edges with have smaller distance than distance
   * set on the input by min_dist_ptr
   *
   * @param v_ptr point coordinates
   * @param nb nb points
   * @param edges range of edges
   * @param min_dist_ptr on return minimal distance, on input starting
   * distance
   * @param o_ptr coordinates of the point on edge
   * @param o_segments closest segments
   * @return MoFEMErrorCode
   */
  MoFEMErrorCode
  findMinDistanceFromTheEdges(const double *v_ptr, const int nb, Range edges,
                              double *min_dist_ptr, double *o_ptr = nullptr,
                              EntityHandle *o_segments = nullptr) const;

  /**@}*/

  /** \name Debugging */

  /**@{*/

  /** \brief Print all DOFs for which element of vector is not a number
   *
   */
  MoFEMErrorCode checkVectorForNotANumber(const Problem *prb_ptr,
                                          const RowColData row_or_col, Vec v);

  /**@}*/

  static MoFEMErrorCode outerProductOfEdgeIntegrationPtsForQuad(MatrixDouble &pts,
                                                              const int edge0,
                                                              const int edge1);
};

double Tools::shapeFunMBTRI0(const double x, const double y) {
  return N_MBTRI0(x, y);
}

double Tools::shapeFunMBTRI1(const double x, const double y) {
  return N_MBTRI1(x, y);
}

double Tools::shapeFunMBTRI2(const double x, const double y) {
  return N_MBTRI2(x, y);
}

template <int LDB>
MoFEMErrorCode Tools::shapeFunMBTRI(double *shape, const double *ksi,
                                    const double *eta, const int nb) {
  MoFEMFunctionBeginHot;
  for (int n = 0; n != nb; ++n) {
    shape[0] = shapeFunMBTRI0(*ksi, *eta);
    shape[1] = shapeFunMBTRI1(*ksi, *eta);
    shape[2] = shapeFunMBTRI2(*ksi, *eta);
    shape += 3;
    ksi += LDB;
    eta += LDB;
  }
  MoFEMFunctionReturnHot(0);
}

double Tools::shapeFunMBTET0(const double x, const double y, const double z) {
  return N_MBTET0(x, y, z);
}

double Tools::shapeFunMBTET1(const double x, const double y, const double z) {
  return N_MBTET1(x, y, z);
}

double Tools::shapeFunMBTET2(const double x, const double y, const double z) {
  return N_MBTET2(x, y, z);
}

double Tools::shapeFunMBTET3(const double x, const double y, const double z) {
  return N_MBTET3(x, y, z);
};

template <int LDB>
MoFEMErrorCode Tools::shapeFunMBTET(double *shape, const double *ksi,
                                    const double *eta, const double *zeta,
                                    const double nb) {
  MoFEMFunctionBeginHot;
  for (int n = 0; n != nb; ++n) {
    shape[0] = shapeFunMBTET0(*ksi, *eta, *zeta);
    shape[1] = shapeFunMBTET1(*ksi, *eta, *zeta);
    shape[2] = shapeFunMBTET2(*ksi, *eta, *zeta);
    shape[3] = shapeFunMBTET3(*ksi, *eta, *zeta);
    shape += 4;
    ksi += LDB;
    eta += LDB;
    zeta += LDB;
  }
  MoFEMFunctionReturnHot(0);
}

} // namespace MoFEM

#endif // __TOOLS_HPP__a

/**
 * \defgroup mofem_tools Tools interface
 * \brief Interface for tools
 *
 * \ingroup mofem
 */
