/** \file CoordSystemsManager.hpp
 * \brief Interface managing coordinate systems set to fields
 *
 * Managing coordinate systems. Each field have attached coordinate system, this
 * allow to pass information about base between independent module
 * implementations.
 *
 */

/*
 * MoFEM is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
 * License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with MoFEM. If not, see <http://www.gnu.org/licenses/>
 */

#ifndef __COORDSTEMSMANAGER_HPP__
#define __COORDSTEMSMANAGER_HPP__

#include "UnknownInterface.hpp"

namespace MoFEM {

static const MOFEMuuid IDD_MOFEMCoordsSystemsManager =
    MOFEMuuid(BitIntefaceId(COORDSSYSTEMMANAGER_INTERFACE));

struct CoordSystemsManager : public UnknownInterface {

  MoFEMErrorCode query_interface(const MOFEMuuid &uuid,
                                 UnknownInterface **iface) const;

  MoFEM::Core &cOre;
  CoordSystemsManager(const MoFEM::Core &core);

  /**
   * \brief Destructor
   */
  ~CoordSystemsManager();

  /**
   * \brief get tags handlers used on meshsets conating information about
   coordinate systems

   */
  MoFEMErrorCode getTags(int verb = -1);

  /**
   * \brief clear multi-index container
   * @return error code
   */
  MoFEMErrorCode clearMap();

  /**
   * \brier initialize container form data on mesh
   * @return error code
   */
  MoFEMErrorCode initialiseDatabaseFromMesh(int verb = 0);

  inline CoordSys_multiIndex &getCoordinateSystemsMultindex() {
    return coordinateSystems;
  }

  /** \brief Add coordinate system

    * \param cs_id see \ref CoordSystems for options
    * \param name unique name of coordinate system
    */
  MoFEMErrorCode addCoordinateSystem(const int cs_dim[], const std::string name,
                                     const enum MoFEMTypes bh = MF_EXCL);

  /** \brief Set coordinate system to field

    * \param name of field
    * \param name unique name of coordinate system
    *
    */
  MoFEMErrorCode setFieldCoordinateSystem(const std::string field_name,
                                          const std::string cs_name);

  /**
   * @brief Get the Coord Sys Ptr object
   * 
   * @param id 
   * @param cs_ptr 
   * @return MoFEMErrorCode 
   */
  MoFEMErrorCode getCoordSysPtr(const EntityHandle id,
                                boost::shared_ptr<CoordSys> &cs_ptr);

  /**
   * @brief Get the Coord Sys Ptr object
   * 
   * @param name 
   * @param cs_ptr 
   * @return MoFEMErrorCode 
   */
  MoFEMErrorCode getCoordSysPtr(const string name,
                                boost::shared_ptr<CoordSys> &cs_ptr);

  /**
   * @brief Get the Coord Sys Ptr object
   * 
   * @param name 
   * @return boost::shared_ptr<CoordSys> 
   */
  boost::shared_ptr<CoordSys> getCoordSysPtr(const string name);

  inline Tag get_th_CoordSysName() const { return th_CoordSysName; }

protected:
  Tag th_CoordSysName; ///< Name of coordinate system
  Tag th_CoordSysDim;  ///< Tag on cordinate sys meshset for dimension of
                       ///< coordinate system associated to fields

  /// Coordinate systems multi-index
  CoordSys_multiIndex coordinateSystems;
};

} // namespace MoFEM

#endif //__COORDYSTEMSINTERFACE_HPP__
