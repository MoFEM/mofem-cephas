/** \file BcManager.hpp
 * \brief Manage boundary conditions set to the problem
 * \ingroup bc_manager
 *
 */

/* MoFEM is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
 * License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with MoFEM. If not, see <http://www.gnu.org/licenses/>
 */

#ifndef __BCMANAGER_HPP__
#define __BCMANAGER_HPP__

#include "UnknownInterface.hpp"

namespace MoFEM {

static const MOFEMuuid IDD_MOFEMBcManager =
    MOFEMuuid(BitIntefaceId(BC_MANAGER));

/**
 * \brief Simple interface for fast problem set-up
 * \ingroup mofem_simple_interface
 */
struct BcManager : public UnknownInterface {

  MoFEMErrorCode query_interface(const MOFEMuuid &uuid,
                                 UnknownInterface **iface) const;

  BcManager(const MoFEM::Core &core);
  virtual ~BcManager() = default;

  /**
   * @brief Data structure storing bc markers and atributes
   *
   */
  struct BCs : boost::enable_shared_from_this<BCs> {
    Range bcEdges;
    std::vector<double> bcAttributes;
    std::vector<unsigned char> bcMarkers;
    inline auto getBcEdgesPtr() {
      return boost::shared_ptr<Range>(shared_from_this(), &bcEdges);
    }
    inline auto getBcMarkersPtr() {
      return boost::shared_ptr<std::vector<unsigned char>>(shared_from_this(),
                                                           &bcMarkers);
    }
  };

  /**
   * \brief get options
   * @return error code
   */
  MoFEMErrorCode getOptions();

  Range getAdjEnts(Range ents);

  /**
   * @brief Remove DOFs from problem
   *
   * @param problem_name
   * @param block_name
   * @param field_name
   * @param lo lowest coefficient
   * @param hi highest coefficient
   * @param get_low_dim_ents get lower dimension entities
   * @return MoFEMErrorCode
   */
  MoFEMErrorCode removeBlockDOFsOnEntities(const std::string problem_name,
                                           const std::string block_name,
                                           const std::string field_name, int lo,
                                           int hi,
                                           bool get_low_dim_ents = true);

  /**
   * @brief Mark block dofs
   *
   * @param problem_name
   * @param block_name
   * @param field_name
   * @param lo lowest coefficient
   * @param hi highest coefficient
   * @param get_low_dim_ents get lower dimension entities
   * @return MoFEMErrorCode
   */
  MoFEMErrorCode pushMarkDOFsOnEntities(const std::string problem_name,
                                        const std::string block_name,
                                        const std::string field_name, int lo,
                                        int hi, bool get_low_dim_ents = true);


  /**
   * @brief Get bc data and remove element
   * 
   * @param block_name 
   * @return boost::shared_ptr<BCs> 
   */
  boost::shared_ptr<BCs> popMarkDOFsOnEntities(const std::string block_name);

  /** \todo Add markers for standard BCs from cubit on Nodests and Sidesets used
   * bu cubit for displacements, forces, etc. Also add function to add blockset
   * by id and type.
   */

  using BcMapByBlockName = std::map<string, boost::shared_ptr<BCs>>;

  /**
   * @brief Get the bc structure object
   *
   * @param block_name
   * @return auto
   */
  inline auto getBcStructure(const std::string bc_id) {
    return bcMapByBlockName.at(bc_id);
  }

  /**
   * @brief Get the bc map
   *
   * @return auto
   */
  inline BcMapByBlockName &getBcMapByBlockName() { return bcMapByBlockName; }

private:
  MoFEM::Core &cOre;

  BcMapByBlockName bcMapByBlockName;
};

} // namespace MoFEM

#endif //__BCMANAGER_HPP__

/**
 * \defgroup bc_manager Manages boundary conditions
 * \brief Implementation manages boundary conditions
 *
 * \ingroup mofem
 **/
