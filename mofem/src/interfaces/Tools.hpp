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

  /**@{*/

  template <class T> static inline double dEterminant(T &t) {
    return t(0, 0) * t(1, 1) * t(2, 2) + t(1, 0) * t(2, 1) * t(0, 2) +
           t(2, 0) * t(0, 1) * t(1, 2) - t(0, 0) * t(2, 1) * t(1, 2) -
           t(2, 0) * t(1, 1) * t(0, 2) - t(1, 0) * t(0, 1) * t(2, 2);
  }

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

  static constexpr double diffNMBTET0x =
      diffN_MBTET0x; ///< derivative of tetrahedral shape function
  static constexpr double diffNMBTET0y =
      diffN_MBTET0y; ///< derivative of tetrahedral shape function
  static constexpr double diffNMBTET0z =
      diffN_MBTET0z; ///< derivative of tetrahedral shape function
  static constexpr double diffNMBTET1x =
      diffN_MBTET1x; ///< derivative of tetrahedral shape function
  static constexpr double diffNMBTET1y =
      diffN_MBTET1y; ///< derivative of tetrahedral shape function
  static constexpr double diffNMBTET1z =
      diffN_MBTET1z; ///< derivative of tetrahedral shape function
  static constexpr double diffNMBTET2x =
      diffN_MBTET2x; ///< derivative of tetrahedral shape function
  static constexpr double diffNMBTET2y =
      diffN_MBTET2y; ///< derivative of tetrahedral shape function
  static constexpr double diffNMBTET2z =
      diffN_MBTET2z; ///< derivative of tetrahedral shape function
  static constexpr double diffNMBTET3x =
      diffN_MBTET3x; ///< derivative of tetrahedral shape function
  static constexpr double diffNMBTET3y =
      diffN_MBTET3y; ///< derivative of tetrahedral shape function
  static constexpr double diffNMBTET3z =
      diffN_MBTET3z; ///< derivative of tetrahedral shape function

  static constexpr std::array<double, 12> diffNMBTET = {

      diffNMBTET0x, diffNMBTET0y, diffNMBTET0z,

      diffNMBTET1x, diffNMBTET1y, diffNMBTET1z,

      diffNMBTET2x, diffNMBTET2y, diffNMBTET2z,

      diffNMBTET3x, diffNMBTET3y, diffNMBTET3z};

  static inline double nMBTET0(const double x, const double y, const double z) {
    return N_MBTET0(x, y, z);
  }

  static inline double nMBTET1(const double x, const double y, const double z) {
    return N_MBTET1(x, y, z);
  }

  static inline double nMBTET2(const double x, const double y, const double z) {
    return N_MBTET2(x, y, z);
  }

  static inline double nMBTET3(const double x, const double y, const double z) {
    return N_MBTET3(x, y, z);
    ;
  };

  static constexpr double nMBTET0At000 = N_MBTET0(0, 0, 0);
  static constexpr double nMBTET1At000 = N_MBTET1(0, 0, 0);
  static constexpr double nMBTET2At000 = N_MBTET2(0, 0, 0);
  static constexpr double nMBTET3At000 = N_MBTET3(0, 0, 0);

  template <int LDB = 1>
  static MoFEMErrorCode nMBTET(double *shape, const double *ksi,
                               const double *eta, const double *zeta,
                               const double nb) {
    MoFEMFunctionBeginHot;
    for (int n = 0; n != nb; ++n) {
      shape[0] = nMBTET0(*ksi, *eta, *zeta);
      shape[1] = nMBTET1(*ksi, *eta, *zeta);
      shape[2] = nMBTET2(*ksi, *eta, *zeta);
      shape[3] = nMBTET3(*ksi, *eta, *zeta);
      shape += 4;
      ksi += LDB;
      eta += LDB;
      zeta += LDB;
    }
    MoFEMFunctionReturnHot(0);
  }

  static constexpr std::array<double, 4> nMBTETAt000 = {
      nMBTET0At000, nMBTET1At000, nMBTET2At000, nMBTET3At000};

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
   * @param edges range of edges
   * @param min_dist_ptr on return minimal distance, on input starting
   * distance
   * @param o_ptr coordinates of the point on edge
   * @return MoFEMErrorCode
   */
  MoFEMErrorCode findMinDistanceFromTheEdges(const double *v_ptr, Range edges,
                                             double *min_dist_ptr,
                                             double *o_ptr) const;

  /**@}*/

  /** \name Debugging */

  /**@{*/

  /** \brief Print all DOFs for which element of vector is not a number
   *
   */
  MoFEMErrorCode checkVectorForNotANumber(const Problem *prb_ptr,
                                          const RowColData row_or_col, Vec v);

  /**@}*/
};

} // namespace MoFEM

#endif // __TOOLS_HPP__a

/**
 * \defgroup mofem_tools Tools interface
 * \brief Interface for tools
 *
 * \ingroup mofem
 */
