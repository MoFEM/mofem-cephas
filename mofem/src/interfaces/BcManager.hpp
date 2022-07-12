/** \file BcManager.hpp
 * \brief Manage boundary conditions set to the problem
 * \ingroup bc_manager
 *
 */

#ifndef __BCMANAGER_HPP__
#define __BCMANAGER_HPP__

#include "UnknownInterface.hpp"

namespace MoFEM {

/**
 * \brief Simple interface for fast problem set-up
 * \ingroup mofem_simple_interface
 */
struct BcManager : public UnknownInterface {

  MoFEMErrorCode query_interface(boost::typeindex::type_index type_index,
                                 UnknownInterface **iface) const;

  BcManager(const MoFEM::Core &core);
  virtual ~BcManager() = default;

  /**
   * @brief Data structure storing bc markers and atributes
   *
   */
  struct BCs : boost::enable_shared_from_this<BCs> {
    Range bcEnts;
    std::vector<double> bcAttributes;
    std::vector<unsigned char> bcMarkers;

    /// \deprecated use getBcEntsPtr
    DEPRECATED inline auto getBcEdgesPtr() {
      return boost::shared_ptr<Range>(shared_from_this(), &bcEnts);
    }

    inline auto getBcEntsPtr() {
      return boost::shared_ptr<Range>(shared_from_this(), &bcEnts);
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
   * @brief Mark block DOFs
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

  using BcMarkerPtr = boost::shared_ptr<std::vector<char unsigned>>;
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

  /**
   * @brief Get the Merged Boundary Marker object
   *
   * @param bc_regex_vec boundary name regex vector
   * @return boundaryMarker
   */
  BcMarkerPtr getMergedBlocksMarker(std::vector<std::regex> bc_regex_vec);
  /**
   * @brief Get the Merged Boundary Marker object
   *
   * @param bc_names vector of boundary names
   * @return boundaryMarker
   */
  inline auto getMergedBlocksMarker(std::vector<string> bc_names) {
    std::vector<std::regex> reg_vec(bc_names.size());
    for (int i = 0; i != bc_names.size(); ++i) {
      auto full_name = std::string("(.*)_") + bc_names[i] + std::string("(.*)");
      reg_vec[i] = std::regex(full_name);
    }
    return getMergedBlocksMarker(reg_vec);
  }
  /**
   * @brief Get the Merged Blocks Marker object
   *
   * @param boundary_markers_ptr_vec vector of boundary markers to merge
   * @return BcMarkerPtr
   */
  BcMarkerPtr getMergedBlocksMarker(
      const std::vector<BcMarkerPtr> &boundary_markers_ptr_vec);
  /**
   * @brief check if given boundary condition name is in the map bc element
   *
   * @param bc element of the map
   * @param reg bc regex
   * @return auto
   */
  inline auto checkBlock(const std::pair<string, boost::shared_ptr<BCs>> &bc,
                         std::regex reg) {
    return std::regex_match(bc.first, reg);
  }
  /**
   * @brief check if given boundary condition name is in the map bc element
   *
   * @param bc element of the map
   * @param name bc name
   * @return auto
   */
  inline auto
  checkBlock(const std::pair<std::string, boost::shared_ptr<BCs>> &bc,
             std::string name) {
    auto full_name = std::string("(.*)_") + name + std::string("(.*)");
    return checkBlock(bc, std::regex(full_name));
  }

  /**
   * @brief Get block IS
   *
   * @param block_prefix  for hashmap
   * @param block_name    for hash map
   * @param field_name    for hash map and IS
   * @param problem_name  for IS
   * @param lo
   * @param hi
   * @param is_expand is to extend
   * @return SmartPetscObj<IS>
   */
  SmartPetscObj<IS>
  getBlockIS(const std::string block_prefix, const std::string block_name,
             const std::string field_name, const std::string problem_name,
             int lo, int hi, SmartPetscObj<IS> is_expand = SmartPetscObj<IS>());

  /**
   * @brief Get block IS
   *
   * @param problem_name
   * @param block_name
   * @param field_name
   * @param lo
   * @param hi
   * @param is_expand is to extend
   * @return SmartPetscObj<IS>
   */
  SmartPetscObj<IS>
  getBlockIS(const std::string problem_name, const std::string block_name,
             const std::string field_name, int lo, int hi,
             SmartPetscObj<IS> is_expand = SmartPetscObj<IS>());

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
