/** \file Tools.hpp
 * \brief Tools interface
 *
 * Implementatiom of some useful and very often used methods * in MoFEM.
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
                                  Tag th = NULL);

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
