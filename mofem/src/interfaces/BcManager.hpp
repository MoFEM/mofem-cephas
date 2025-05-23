/** \file BcManager.hpp
 * \brief Manage boundary conditions set to the problem
 * \ingroup bc_manager
 *
 */

#ifndef __BCMANAGER_HPP__
#define __BCMANAGER_HPP__

#include "UnknownInterface.hpp"

namespace MoFEM {

template <CubitBC BC> struct BcMeshsetType {};
template <CubitBC BC> struct BcScalarMeshsetType {};
template <CubitBC BC> using BcTemperature = BcScalarMeshsetType<BC>;
template <CubitBC BC> struct BcDisplacementMeshsetType {};
template <CubitBC BC> struct BcForceMeshsetType {};
template <CubitBC BC> struct BcFluxMeshsetType {};

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

    using DofsView = std::vector<boost::weak_ptr<NumeredDofEntity>>;
    boost::shared_ptr<DofsView> dofsViewPtr;

    boost::shared_ptr<DisplacementCubitBcData> dispBcPtr;
    boost::shared_ptr<TemperatureCubitBcData> tempBcPtr;
    boost::shared_ptr<HeatFluxCubitBcData> heatFluxBcPtr;
    boost::shared_ptr<ForceCubitBcData> forceBcPtr;
    boost::shared_ptr<MPCsType> mpcPtr;

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

  // FIXME: add functions to couple the dofs
  MoFEMErrorCode addBlockDOFsToMPCs(const std::string problem_name,
                                    const std::string field_name,
                                    bool get_low_dim_ents = false,
                                    bool block_name_field_prefix = false,
                                    bool is_distributed_mesh = false);

  /**
   * @brief Mark side DOFs 
   * 
   * @param problem_name 
   * @param block_name 
   * @param field_name 
   * @param lo 
   * @param hi 
   * @return MoFEMErrorCode 
   */
  MoFEMErrorCode pushMarkSideDofs(const std::string problem_name,
                                  const std::string block_name,
                                  const std::string field_name, int bridge_dim,
                                  int lo, int hi);

  /**
   * @brief Remove side DOFs
   * 
   * @param problem_name 
   * @param block_name 
   * @param field_name 
   * @param lo 
   * @param hi 
   * @param is_distributed_mesh 
   * @return MoFEMErrorCode 
   */
  MoFEMErrorCode removeSideDOFs(const std::string problem_name,
                                const std::string block_name,
                                const std::string field_name, int bridge_dim,
                                int lo, int hi,
                                bool is_distributed_mesh = true);

  /**
   * @brief Remove DOFs from problem
   *
   * @param problem_name
   * @param block_name
   * @param field_name
   * @param lo lowest coefficient
   * @param hi highest coefficient
   * @param get_low_dim_ents get lower dimension entities
   * @param is_distributed_mesh distributed mesh
   * @return MoFEMErrorCode
   */
  MoFEMErrorCode removeBlockDOFsOnEntities(const std::string problem_name,
                                           const std::string block_name,
                                           const std::string field_name, int lo,
                                           int hi, bool get_low_dim_ents = true,
                                           bool is_distributed_mesh = true);

  /**
   * @brief Mark block DOFs
   *
   * @param problem_name
   * @param block_name
   * @param field_name
   * @param lo lowest coefficient
   * @param hi highest coefficient
   * @param get_low_dim_ents get lower dimension entities
   * field name
   * @return MoFEMErrorCode
   */
  MoFEMErrorCode pushMarkDOFsOnEntities(const std::string problem_name,
                                        const std::string block_name,
                                        const std::string field_name, int lo,
                                        int hi, bool get_low_dim_ents = true);

  /**
   * @brief Mark block DOFs
   *
   * @tparam BCSET
   * @param problem_name
   * @param field_name
   * @param get_low_dim_ents
   * @param is_distributed_mesh
   * @param block_name_field_prefix
   * @return MoFEMErrorCode
   */
  template <typename T>
  MoFEMErrorCode removeBlockDOFsOnEntities(const std::string problem_name,
                                           const std::string field_name,
                                           bool get_low_dim_ents = true,
                                           bool block_name_field_prefix = false,
                                           bool is_distributed_mesh = true);

  /**
   * @brief Mark block DOFs
   *
   * @param problem_name
   * @param field_name
   * @param get_low_dim_ents get lower dimension entities
   * @param block_name_field_prefix
   * @return MoFEMErrorCode
   */
  template <typename T>
  MoFEMErrorCode pushMarkDOFsOnEntities(const std::string problem_name,
                                        const std::string field_name,
                                        bool get_low_dim_ents = true,
                                        bool block_name_field_prefix = false);

  /**
   * @brief Mark block DOFs
   *
   * @tparam BCSET
   * @param problem_name
   * @param field_name
   * @param block_name
   * @param get_low_dim_ents
   * @param is_distributed_mesh
   * @return MoFEMErrorCode
   */
  template <typename T>
  MoFEMErrorCode removeBlockDOFsOnEntities(const std::string problem_name,
                                           const std::string block_name,
                                           const std::string field_name,
                                           bool get_low_dim_ents = true,
                                           bool is_distributed_mesh = true);

  /**
   * @brief Mark block DOFs
   *
   * @param problem_name
   * @param field_name
   * @param block_name
   * @param get_low_dim_ents get lower dimension entities
   * @return MoFEMErrorCode
   */
  template <typename T>
  MoFEMErrorCode pushMarkDOFsOnEntities(const std::string problem_name,
                                        const std::string block_name,
                                        const std::string field_name,
                                        bool get_low_dim_ents = true);

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
   * @brief Merge block ranges
   *
   * @param bc_regex_vec
   * @return Range
   */
  Range getMergedBlocksRange(std::vector<std::regex> bc_regex_vec);

  /**
   * @brief Merge block ranges
   *
   * @param bc_names
   * @return auto
   */
  inline auto getMergedBlocksRange(std::vector<string> bc_names) {
    std::vector<std::regex> reg_vec(bc_names.size());
    for (int i = 0; i != bc_names.size(); ++i) {
      auto full_name = std::string("(.*)_") + bc_names[i] + std::string("(.*)");
      reg_vec[i] = std::regex(full_name);
    }
    return getMergedBlocksRange(reg_vec);
  }

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

  /**
   * @brief Extract block name and block name form block id
   *
   * @param block_id
   * @param prb_name
   * @return std::pair<std::string, std::string>
   */
  static std::pair<std::string, std::string>
  extractStringFromBlockId(const std::string block_id,
                           const std::string prb_name);

private:
  MoFEM::Core &cOre;

  BcMapByBlockName bcMapByBlockName;
};

template <>
MoFEMErrorCode
BcManager::removeBlockDOFsOnEntities<BcMeshsetType<DISPLACEMENTSET>>(
    const std::string problem_name, const std::string field_name,
    bool get_low_dim_ents, bool block_name_field_prefix,
    bool is_distributed_mesh);

template <>
MoFEMErrorCode
BcManager::removeBlockDOFsOnEntities<BcMeshsetType<TEMPERATURESET>>(
    const std::string problem_name, const std::string field_name,
    bool get_low_dim_ents, bool block_name_field_prefix,
    bool is_distributed_mesh);

template <>
MoFEMErrorCode BcManager::removeBlockDOFsOnEntities<BcMeshsetType<HEATFLUXSET>>(
    const std::string problem_name, const std::string field_name,
    bool get_low_dim_ents, bool block_name_field_prefix,
    bool is_distributed_mesh);

template <>
MoFEMErrorCode
BcManager::removeBlockDOFsOnEntities<BcDisplacementMeshsetType<BLOCKSET>>(
    const std::string problem_name, const std::string field_name,
    bool get_low_dim_ents, bool block_name_field_prefix,
    bool is_distributed_mesh);

template <>
MoFEMErrorCode
BcManager::removeBlockDOFsOnEntities<BcScalarMeshsetType<BLOCKSET>>(
    const std::string problem_name, const std::string block_name,
    const std::string field_name, bool get_low_dim_ents,
    bool is_distributed_mesh);

template <>
MoFEMErrorCode
BcManager::pushMarkDOFsOnEntities<BcMeshsetType<DISPLACEMENTSET>>(
    const std::string problem_name, const std::string field_name,
    bool get_low_dim_ents, bool block_name_field_prefix);

template <>
MoFEMErrorCode BcManager::pushMarkDOFsOnEntities<BcMeshsetType<TEMPERATURESET>>(
    const std::string problem_name, const std::string field_name,
    bool get_low_dim_ents, bool block_name_field_prefix);

template <>
MoFEMErrorCode BcManager::pushMarkDOFsOnEntities<BcMeshsetType<HEATFLUXSET>>(
    const std::string problem_name, const std::string field_name,
    bool get_low_dim_ents, bool block_name_field_prefix);

template <>
MoFEMErrorCode
BcManager::pushMarkDOFsOnEntities<BcDisplacementMeshsetType<BLOCKSET>>(
    const std::string problem_name, const std::string field_name,
    bool get_low_dim_ents, bool block_name_field_prefix);

template <>
MoFEMErrorCode BcManager::pushMarkDOFsOnEntities<BcScalarMeshsetType<BLOCKSET>>(
    const std::string problem_name, const std::string field_name,
    const std::string block_name, bool get_low_dim_ents);

template <>
MoFEMErrorCode BcManager::pushMarkDOFsOnEntities<DisplacementCubitBcData>(
    const std::string problem_name, const std::string field_name,
    bool get_low_dim_ents, bool block_name_field_prefix);

template <>
MoFEMErrorCode BcManager::removeBlockDOFsOnEntities<DisplacementCubitBcData>(
    const std::string problem_name, const std::string field_name,
    bool get_low_dim_ents, bool block_name_field_prefix,
    bool is_distributed_mesh);

template <>
MoFEMErrorCode BcManager::pushMarkDOFsOnEntities<TemperatureCubitBcData>(
    const std::string problem_name, const std::string field_name,
    bool get_low_dim_ents, bool block_name_field_prefix);

template <>
MoFEMErrorCode BcManager::removeBlockDOFsOnEntities<TemperatureCubitBcData>(
    const std::string problem_name, const std::string field_name,
    bool get_low_dim_ents, bool block_name_field_prefix,
    bool is_distributed_mesh);

template <>
MoFEMErrorCode BcManager::pushMarkDOFsOnEntities<HeatFluxCubitBcData>(
    const std::string problem_name, const std::string field_name,
    bool get_low_dim_ents, bool block_name_field_prefix);

template <>
MoFEMErrorCode BcManager::removeBlockDOFsOnEntities<HeatFluxCubitBcData>(
    const std::string problem_name, const std::string field_name,
    bool get_low_dim_ents, bool block_name_field_prefix,
    bool is_distributed_mesh);

template <>
MoFEMErrorCode BcManager::pushMarkDOFsOnEntities<BcMeshsetType<FORCESET>>(
    const std::string problem_name, const std::string field_name,
    bool get_low_dim_ents, bool block_name_field_prefix);

template <>
MoFEMErrorCode BcManager::pushMarkDOFsOnEntities<BcForceMeshsetType<BLOCKSET>>(
    const std::string problem_name, const std::string field_name,
    bool get_low_dim_ents, bool block_name_field_prefix);

template <>
MoFEMErrorCode BcManager::pushMarkDOFsOnEntities<ForceCubitBcData>(
    const std::string problem_name, const std::string field_name,
    bool get_low_dim_ents, bool block_name_field_prefix);

template <>
MoFEMErrorCode BcManager::removeBlockDOFsOnEntities<BcMeshsetType<FORCESET>>(
    const std::string problem_name, const std::string field_name,
    bool get_low_dim_ents, bool block_name_field_prefix,
    bool is_distributed_mesh);

template <>
MoFEMErrorCode
BcManager::removeBlockDOFsOnEntities<BcForceMeshsetType<BLOCKSET>>(
    const std::string problem_name, const std::string field_name,
    bool get_low_dim_ents, bool block_name_field_prefix,
    bool is_distributed_mesh);

template <>
MoFEMErrorCode BcManager::removeBlockDOFsOnEntities<ForceCubitBcData>(
    const std::string problem_name, const std::string field_name,
    bool get_low_dim_ents, bool block_name_field_prefix,
    bool is_distributed_mesh);

} // namespace MoFEM

#endif //__BCMANAGER_HPP__

/**
 * \defgroup bc_manager Manages boundary conditions
 * \brief Implementation manages boundary conditions
 *
 * \ingroup mofem
 **/
