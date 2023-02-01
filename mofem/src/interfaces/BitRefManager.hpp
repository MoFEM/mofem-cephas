/** \file BitRefManager.hpp
 * \brief Interface managing BitRefLevels
 * \ingroup mofem_bit_ref
 *
 * Managing BitRef levels
 *
 */

#ifndef __BITREFMANAGER_HPP__
#define __BITREFMANAGER_HPP__

#include "UnknownInterface.hpp"

namespace MoFEM {

/**
 * \brief Managing BitRefLevels
 * \ingroup mofem_bit_ref
 * \nosubgrouping
 */
struct BitRefManager : public UnknownInterface {

  MoFEMErrorCode query_interface(boost::typeindex::type_index type_index,
                                 UnknownInterface **iface) const;

  MoFEM::Core &cOre;
  bool dEbug;

  BitRefManager(const MoFEM::Core &core);

  /**
   * \brief Destructor
   */
  virtual ~BitRefManager() = default;

  /** \name Setting and shifting bits */

  /**@{*/

  /**
   * \brief add entities to database and set bit ref level
   * \ingroup mofem_bit_ref
   *
   * This function set bit ref level, add entries to core database and
   * create ref finite elements. Finite elements are create of entities in
   function
   * argument, whereas all lower dimension entities are added as a field
   entities
   *
   *

   Example:\code
   EntityHandle meshset1; //contains ent1,ent2,ent3
   BitRefLevel myLevel0;
   myLevel0.set(0);
   m_field.getInterface<BitRefManager>()->setBitRefLevelByDim(meshset1,3,myLevel0);
   //refine meshset1 into meshset2 and get new ents which are ent4, ent5
   EntityHandle meshset2; //contains ent1,ent2,ent3,ent4,ent5
   BitRefLevel myLevel1;
   myLevel1.set(1);
   m_field.getInterface<BitRefManager>()->setBitRefLevelByDim(meshset2,3,myLevel1);
   \endcode

   * So entities 1,2,3 would be assigned to bit level 0 and 1 <br>
   * ent1[1,1,0,0,0,0,0], ent2[1,1,0,0,0,0,0], ent3[1,1,0,0,0,0,0], <br>
   * and entities 4 and 5 are assigned to bit level 1 only <br>
   * ent4[0,1,0,0,0,0,0], ent5[0,1,0,0,0,0,0] <br>


   * @param  ents      entities to set
   * @param  bit       bit refinement level
   * @param  only_tets only add entities on tetrahedral (obsolete need to be
   fixed)
   * @param  verb      verbosity level
   * @param adj_ents_ptr if pointer is given, it is used to get adj entities by
   dimension/type
   * @return           error code
   */
  MoFEMErrorCode setBitRefLevel(const Range &ents, const BitRefLevel bit,
                                const bool only_tets = true, int verb = 0,
                                Range *adj_ents_ptr = nullptr) const;

  /**
   * @brief add entities to database and set bit ref level
   * \ingroup mofem_bit_ref
   *
   * \note In THIS variant only entities in range are added and ref finite
   * elements reated.
   *
   * @param ents
   * @param bit
   * @param only_tets
   * @param verb
   * @return MoFEMErrorCode setBitRefLevel
   */
  MoFEMErrorCode setElementsBitRefLevel(const Range &ents,
                                        const BitRefLevel bit = BitRefLevel(),
                                        int verb = QUIET) const;

  /**
   * @brief add entities to database and set bit ref level
   * \ingroup mofem_bit_ref
   *
   * \note In THIS variant only entities in range are added. And DO NOT create
   * elements.
   *
   * @param ents
   * @param bit
   * @param only_tets
   * @param verb
   * @return MoFEMErrorCode setBitRefLevel
   */
  MoFEMErrorCode setEntitiesBitRefLevel(const Range &ents,
                                        const BitRefLevel bit = BitRefLevel(),
                                        int verb = QUIET) const;

  /**
   * @brief Set the bit ref level to entities in the field meshset
   * 
   * @param field_name 
   * @param bit 
   * @param verb 
   * @return MoFEMErrorCode 
   */
  MoFEMErrorCode
  setFieldEntitiesBitRefLevel(const std::string field_name,
                              const BitRefLevel bit = BitRefLevel(),
                              int verb = QUIET) const;

  /**
   * @brief Set the Bit Ref Level By Dim object
   *
   * \note In THIS variant only entities in range are added. And DO NOT create
   * elements.
   *
   * @param meshset
   * @param dim
   * @param bit
   * @param verb
   * @return MoFEMErrorCode
   */
  MoFEMErrorCode setBitRefLevelByDim(const EntityHandle meshset, const int dim,
                                     const BitRefLevel bit,
                                     int verb = QUIET) const;

  /**
   * @brief Set the Bit Ref Level By Type object
   * 
   * \note In THIS variant only entities in range are added. And DO NOT create
   * elements.
   * 
   * @param meshset 
   * @param type 
   * @param bit 
   * @param verb 
   * @return MoFEMErrorCode 
   */
  MoFEMErrorCode setBitRefLevelByType(const EntityHandle meshset,
                                      const EntityType type,
                                      const BitRefLevel bit,
                                      int verb = QUIET) const;

  /** brief add meshset and set bit ref level
   * \ingroup mofem_bit_ref
   *
   * \param EntityHandle MeshSet
   * \param BitRefLevel bitLevel
   */
  MoFEMErrorCode setBitLevelToMeshset(const EntityHandle meshset,
                                      const BitRefLevel bit,
                                      int verb = 0) const;

  /**
   * @brief Add entities which exist in MoAB database, and have set appropiate
   * BitRef level tag, to multi-indices in MoFEM.
   *
   * \note Every entity, used for create DoFS, or elements has to have set
   * BitRefLevel and be added to MoFEM database.
   *
   * \note This functions add lower dimension entities by calling
   * setEntitiesBitRefLevel
   *
   * @param type of entity
   * @param bit bit ref level
   * @param mask 
   * @param verb verbosity level
   * @return MoFEMErrorCode
   */
  MoFEMErrorCode
  addToDatabaseBitRefLevelByType(const EntityType type, const BitRefLevel bit,
                                 const BitRefLevel mask = BitRefLevel().set(),
                                 int verb = QUIET) const;

  /**
   * @brief Add entities which exist in MoAB database, and have set appropiate
   * BitRef level tag, to multi-indices in MoFEM.
   *
   * \note Every entity, used for create DoFS, or elements has to have set
   * BitRefLevel and be added to MoFEM database.
   *
   * \note This functions add lower dimension entities by calling
   * setEntitiesBitRefLevel
   *
   * @param dim dimension of entity
   * @param bit bit ref level
   * @param verb verbosity level
   * @return MoFEMErrorCode
   */
  MoFEMErrorCode
  addToDatabaseBitRefLevelByDim(const int dim, const BitRefLevel bit,
                                const BitRefLevel mask = BitRefLevel().set(),
                                int verb = QUIET) const;

  /**
   * @brief Process bit ref level by lambda function
   *
   * \note That not apply to type of MBENTITY. To avoid problems with problem
   * meshsets. 
   *
   * @param fun
   * @return MoFEMErrorCode
   */
  MoFEMErrorCode lambdaBitRefLevel(
      boost::function<void(EntityHandle ent, BitRefLevel &bit)> fun) const;

  /**
   * @brief Process bit ref level by lambda function
   * 
   * @param fun 
   * @return MoFEMErrorCode 
   */
  MoFEMErrorCode lambdaBitRefLevel(
      const Range &ents,
      boost::function<void(EntityHandle ent, BitRefLevel &bit)> fun) const;

  /**
   * \brief add bit ref level to ref entity
   * \ingroup mofem_bit_ref
   * @param  ents range of entities
   * @param  bit  bit ref level
   * @param  verb verbosity level
   * @return      error code
   */
  MoFEMErrorCode addBitRefLevel(const Range &ents, const BitRefLevel &bit,
                                int verb = QUIET) const;

  /**
   * @brief add bit ref level by dimension
   *
   * @param meshset
   * @param dim dimension of entities
   * @param bit added bit
   * @param verb verbosity level
   * @return MoFEMErrorCode
   */
  MoFEMErrorCode addBitRefLevelByDim(const EntityHandle meshset, const int dim,
                                     const BitRefLevel bit,
                                     int verb = QUIET) const;

  /**
   * \brief Set nth bit ref level
   *
   * \note This function modify bits only on entities in RefEntity_multiindex
   *
   * @param  ents entities to set bit ref level
   * @param  n    nth bit
   * @param  b    value to set
   * @return      error code
   */
  MoFEMErrorCode setNthBitRefLevel(const Range &ents, const int n, const bool b,
                                   int verb = QUIET) const;

  /**
   * \brief Set nth bit ref level to all entities in database
   * \ingroup mofem_bit_ref
   * @param  n    nth bit
   * @param  b    value to set
   * @return      error code
   */
  MoFEMErrorCode setNthBitRefLevel(const int n, const bool b,
                                   int verb = QUIET) const;

  /** \brief left shift bit ref level
   * \ingroup mofem_bit_ref
   * this results of deletion of entities on far left side
   *
   * \note Not implemented
   */
  MoFEMErrorCode shiftLeftBitRef(const int shift,
                                 const BitRefLevel mask = BitRefLevel().set(),
                                 int verb = DEFAULT_VERBOSITY) const;

  /** \brief right shift bit ref level
   * \ingroup mofem_bit_ref
   */
  MoFEMErrorCode shiftRightBitRef(const int shift,
                                  const BitRefLevel mask = BitRefLevel().set(),
                                  int verb = DEFAULT_VERBOSITY,
                                  MoFEMTypes mf = MF_ZERO) const;

  /**@}*/

  /** \name Entity handlers by bit ref level */

  /**@{*/

  /**
   * @brief filter entities by bit ref level
   *
   * @param bit
   * @param mask
   * @param ents
   * @param QUIET
   * @return MoFEMErrorCode
   */
  MoFEMErrorCode filterEntitiesByRefLevel(const BitRefLevel bit,
                                          const BitRefLevel mask, Range &ents,
                                          int verb = QUIET) const;

  /**\brief add all ents from ref level given by bit to meshset
   * \ingroup mofem_bit_ref
   *
   * \note Entities NOT have to be added to MoFEM database
   *
   * \param BitRefLevel bitLevel
   * \param BitRefLevel mask
   * \param EntityType type of entities
   * \retval EntityHandle meshset
   *
   */
  MoFEMErrorCode getEntitiesByTypeAndRefLevel(const BitRefLevel bit,
                                              const BitRefLevel mask,
                                              const EntityType type,
                                              const EntityHandle meshset,
                                              int verb = 0) const;

  /**\brief add all ents from ref level given by bit to meshset
   * \ingroup mofem_bit_ref
   *
   * \note Entities NOT have to be added to MoFEM database
   *
   * \param BitRefLevel bitLevel
   * \param BitRefLevel mask
   * \param EntityType type of entities
   * \retval ents
   *
   */
  MoFEMErrorCode getEntitiesByTypeAndRefLevel(const BitRefLevel bit,
                                              const BitRefLevel mask,
                                              const EntityType type,
                                              Range &ents, int verb = 0) const;

  /**\brief add all ents from ref level given by bit to meshset
   * \ingroup mofem_bit_ref
   *
   * \note Entities NOT have to be added to MoFEM database
   *
   * \param BitRefLevel bitLevel
   * \param BitRefLevel mask
   * \param EntityType dimension of entities
   * \retval ents
   *
   */
  MoFEMErrorCode getEntitiesByDimAndRefLevel(const BitRefLevel bit,
                                             const BitRefLevel mask,
                                             const int dim,
                                             const EntityHandle meshset,
                                             int verb = 0) const;

  /**\brief add all ents from ref level given by bit to meshset
   * \ingroup mofem_bit_ref
   *
   * \note Entities NOT have to be added to MoFEM database
   *
   * \param BitRefLevel bitLevel
   * \param BitRefLevel mask
   * \param EntityType dimension of entities
   * \retval ents
   *
   */
  MoFEMErrorCode getEntitiesByDimAndRefLevel(const BitRefLevel bit,
                                             const BitRefLevel mask,
                                             const int dim, Range &ents,
                                             int verb = 0) const;

  /**\brief add all ents from ref level given by bit to meshset
   * \ingroup mofem_bit_ref
   *
   * \note Entities NOT have to be added to MoFEM database
   *
   * \param BitRefLevel bitLevel
   * \param BitRefLevel mask
   * \param EntityHandle meshset
   *
   */
  MoFEMErrorCode getEntitiesByRefLevel(const BitRefLevel bit,
                                       const BitRefLevel mask,
                                       const EntityHandle meshset,
                                       const int verb = QUIET) const;

  /**\brief add all ents from ref level given by bit to meshset
   * \ingroup mofem_bit_ref
   *
   * \note Entities NOT have to be added to MoFEM database
   *
   * \param BitRefLevel bitLevel
   * \param BitRefLevel mask
   * \retval ents
   */
  MoFEMErrorCode getEntitiesByRefLevel(const BitRefLevel bit,
                                       const BitRefLevel mask, Range &ents,
                                       const int verb = QUIET) const;

  /**
   * \brief get entities by bit ref level and type of parent
   *
   * \note Entities have to be added to MoFEM database
   *
   * \param BitRefLevel bitLevel
   * \param BitRefLevel mask
   * @param  type of parent
   * @param  ents returned ents
   * @return      error code
   */
  MoFEMErrorCode getEntitiesByParentType(const BitRefLevel bit,
                                         const BitRefLevel mask,
                                         const EntityType type, Range &ents,
                                         const int verb = QUIET) const;

  /**
   * @brief Get all entities not in database
   *
   * @param ents
   * @return MoFEMErrorCode
   */
  MoFEMErrorCode getAllEntitiesNotInDatabase(Range &ents) const;

  /**
   * @brief \brief Get entities not in database
   *
   * @param ents
   * @return MoFEMErrorCode
   */
  MoFEMErrorCode filterEntitiesNotInDatabase(Range &ents) const;

  /**@}*/

  /** \name Get adjacencies bit ref level */

  /**@{*/

  /** \brief Get the adjacencies associated with a entity to entities of a
   * specified dimension.
   * \ingroup mofem_bit_ref
   *
   * bit ref level of adjacent entities is equal to bit ref level of adjacent
   * entities
   */
  virtual MoFEMErrorCode getAdjacenciesEquality(const EntityHandle from_entity,
                                                const int to_dimension,
                                                Range &adj_entities) const;

  /** \brief Get the adjacencies associated with a entity to entities of a
   * specified dimension.
   * \ingroup mofem_bit_ref
   *
   * bit ref level of adjacent entities is any of bit ref level of adjacent
   * entities
   */
  virtual MoFEMErrorCode getAdjacenciesAny(const EntityHandle from_entity,
                                           const int to_dimension,
                                           Range &adj_entities) const;

  /** \brief Get the adjacencies associated with a entity to entities of a
   * specified dimension.
   * \ingroup mofem_bit_ref
   *
   * bit ref level of adjacent entities is equal to bit ref level of adjacent
   * entities
   */
  virtual MoFEMErrorCode
  getAdjacencies(const Problem *problem_ptr, const EntityHandle *from_entities,
                 const int num_entities, const int to_dimension,
                 Range &adj_entities,
                 const int operation_type = moab::Interface::INTERSECT,
                 const int verb = 0) const;

  /** \brief Get the adjacencies associated with a entity to entities of a
   * specified dimension.
   * \ingroup mofem_bit_ref
   *
   * bit ref level of adjacent entities is equal to bit ref level of adjacent
   * entities
   */
  virtual MoFEMErrorCode
  getAdjacencies(const BitRefLevel bit, const EntityHandle *from_entities,
                 const int num_entities, const int to_dimension,
                 Range &adj_entities,
                 const int operation_type = moab::Interface::INTERSECT,
                 const int verb = 0) const;

  /**@}*/

  /** \name Update meshsets and ranges by children */

  /**@{*/

  /** \brief Get child entities form meshset containing parent entities
   * \ingroup mofem_bit_ref
   *
   * Search for refined entities of given type whose parent are entities in
   *the
   * parent meshset. It can be used for example to transfer information about
   * boundary conditions to refined mesh or split mesh by interface
   * elements. It is used by function refineMeshset, to update MESHSET
   *finite elements.
   *
   * \param parent meshset
   * \param parent_bit refinement level
   * \param mask of parent bit ref level
   * \param child_bit refinement level
   * \param mask of child bit ref level
   * \param child meshset where child entities are stored (if the child
   *meshset is set to be the parent meshset, the parent would be updated with
   *the refined entities)
   * \param child_type meshset is update only by entities of specified type. if
   *type is set to MBMAXTYPE all types are updated.
   * \param recursive if true parent meshset is searched recursively
   *
   **/
  MoFEMErrorCode updateMeshsetByEntitiesChildren(
      const EntityHandle parent, const BitRefLevel &parent_bit,
      const BitRefLevel &parent_mask, const BitRefLevel &child_bit,
      const BitRefLevel &child_mask, const EntityHandle child,
      EntityType child_type, const bool recursive = false, int verb = 0);

  /** \copydoc updateMeshsetByEntitiesChildren
   **/
  MoFEMErrorCode updateMeshsetByEntitiesChildren(const EntityHandle parent,
                                                 const BitRefLevel &child_bit,
                                                 const EntityHandle child,
                                                 EntityType child_type,
                                                 const bool recursive = false,
                                                 int verb = 0);

  /** \brief update fields meshesets by child entities
   * \ingroup mofem_update_meshsets_and_ranges
   *
   * \note This calls updateMeshsetByEntitiesChildren for all entity types.
   *
   */
  MoFEMErrorCode
  updateFieldMeshsetByEntitiesChildren(const BitRefLevel &child_bit,
                                       int verb = 0);

  /** \brief update field meshset by child entities
   * \ingroup mofem_update_meshsets_and_ranges
   */
  MoFEMErrorCode updateFieldMeshsetByEntitiesChildren(
      const std::string name, const BitRefLevel &child_bit, int verb = 0);

  /** \brief update finite element meshset by child entities
   * \ingroup mofem_update_meshsets_and_ranges
   */
  MoFEMErrorCode updateFiniteElementMeshsetByEntitiesChildren(
      const std::string name, const BitRefLevel &child_bit,
      const EntityType fe_ent_type, int verb = 0);

  /**
   * \brief Update range by childrens
   *
   * FIXME: NOT TESTED
   *
   * @param  parent parent range
   * @param  child  children range
   * @return        error code
   */
  MoFEMErrorCode updateRangeByChildren(const Range &parent, Range &child,
                                       MoFEMTypes bh = MF_ZERO);

  /**
   * \brief Update range by parents
   *
   * FIXME: NOT TESTED
   *
   * @param  child parent range
   * @param  parent  children range
   * @return        error code
   */
  MoFEMErrorCode updateRangeByParent(const Range &child_ents,
                                     Range &parent_ents,
                                     MoFEMTypes bh = MF_ZERO);

  /**
   * @deprecated use updateRangeByChildren
   */
  DEPRECATED inline MoFEMErrorCode
  updateRange(const Range &parent, Range &child, MoFEMTypes bh = MF_ZERO) {
    return updateRangeByChildren(parent, child, bh);
  }

  /**@}*/

  /** \name Writing files */

  /**@{*/

  /**
   * \brief Write bit ref level to file
   * @param  bit       bit ref level
   * @param  mask      mask of bit ref level
   * @param  dim       dimension
   * @param  file_name file name (see moab documentation)
   * @param  file_type file type (see moab documentation)
   * @param  options   file options (see moab documentation)
   * @return           error code
   */
  MoFEMErrorCode writeBitLevel(const BitRefLevel bit, const BitRefLevel mask,
                               const char *file_name, const char *file_type,
                               const char *options,
                               const bool check_for_empty = true) const;

  /**
   * \brief Write bit ref level to file
   * @param  bit       bit ref level
   * @param  mask      mask of bit ref level
   * @param  dim       dimension
   * @param  file_name file name (see moab documentation)
   * @param  file_type file type (see moab documentation)
   * @param  options   file options (see moab documentation)
   * @return           error code
   */
  MoFEMErrorCode writeBitLevelByDim(const BitRefLevel bit,
                                    const BitRefLevel mask, const int dim,
                                    const char *file_name,
                                    const char *file_type, const char *options,
                                    const bool check_for_empty = true) const;

  /**
   * \brief Write bit ref level to file
   * @param  bit       bit ref level
   * @param  mask      mask of bit ref level
   * @param  type      type of entity
   * @param  file_name file name (see moab documentation)
   * @param  file_type file type (see moab documentation)
   * @param  options   file options (see moab documentation)
   * @return           error code
   */
  MoFEMErrorCode writeBitLevelByType(const BitRefLevel bit,
                                     const BitRefLevel mask,
                                     const EntityType type,
                                     const char *file_name,
                                     const char *file_type, const char *options,
                                     const bool check_for_empty = true) const;

  /**
   * @brief Write ents not in database
   *
   * @param file_name
   * @param file_type for example "VTK"
   * @param options
   * @param check_for_empty
   * @return MoFEMErrorCode
   */
  MoFEMErrorCode
  writeEntitiesNotInDatabase(const char *file_name, const char *file_type,
                             const char *options,
                             const bool check_for_empty = true) const;

  /**
   * @brief Write all entities by bit levels and type
   * 
   * @param mask 
   * @param type 
   * @param file_name 
   * @param file_type 
   * @param options 
   * @return MoFEMErrorCode 
   */
  MoFEMErrorCode writeEntitiesAllBitLevelsByType(const BitRefLevel mask,
                                                 const EntityType type,
                                                 const char *file_name,
                                                 const char *file_type,
                                                 const char *options);

  /**@}*/

  /**@{*/

  /**
   * @brief Fix tag size when BITREFLEVEL_SIZE of core library is different than
   * file BITREFLEVEL_SIZE
   *
   * @return MoFEMErrorCode
   */
  static MoFEMErrorCode fixTagSize(moab::Interface &moab,
                                   bool *changed = nullptr);

  /**@}*/


  /** \name Get tag handles to data on the mesh */

  /**@{*/

  inline Tag get_th_RefParentHandle() const {
    return cOre.get_th_RefParentHandle();
  }
  inline Tag get_th_RefBitLevel() const { return cOre.get_th_RefBitLevel(); }

  /**@}*/

};
} // namespace MoFEM

#endif //__BITREFMANAGER_HPP__

/**
 * \defgroup mofem_bit_ref BitRefManager
 * \brief Managing BitRefLevels
 *
 * \ingroup mofem
 */
