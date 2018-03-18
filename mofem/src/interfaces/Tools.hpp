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

  static const MOFEMuuid IDD_MOFEMTools = MOFEMuuid( BitIntefaceId(TOOLS) );

  /**
   * \brief Auxiliary tools
   * \nosubgrouping
   * \ingroup mofem_tools
   */
  struct Tools: public UnknownInterface {

    MoFEMErrorCode query_interface(const MOFEMuuid &uuid,
                                   UnknownInterface **iface) const;

    MoFEM::Core &cOre;
    Tools(const MoFEM::Core &core) : cOre(const_cast<MoFEM::Core &>(core)) {}

    /** \name Computational */

    /**@{*/

    template<class T>
    static inline double dEterminant(T &t) {
      return
      +t(0,0)*t(1,1)*t(2,2) + t(1,0)*t(2,1)*t(0,2)
      +t(2,0)*t(0,1)*t(1,2) - t(0,0)*t(2,1)*t(1,2)
      -t(2,0)*t(1,1)*t(0,2) - t(1,0)*t(0,1)*t(2,2);
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
                                  Tag th = NULL,
                                  boost::function<double(double, double)> f =
                                      [](double a, double b) -> double {
                                    return std::min(a, b);
                                  });

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
    getTetsWithQuality(Range &out_tets, const Range &tets, Tag th = NULL,
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
    MoFEMErrorCode writeTetsWithQuality(
        const char *file_name, const char *file_type, const char *options,
        const Range &tets, Tag th = NULL,
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

}

#endif // __TOOLS_HPP__a

/***************************************************************************//**
 * \defgroup mofem_tools Tools interface
 * \brief Interface for tools
 *
 * \ingroup mofem
 ******************************************************************************/
