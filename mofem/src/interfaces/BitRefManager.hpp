/** \file BitRefManager.hpp
 * \brief Interface managing BitRefLevels
 * \ingroup mofem_bit_ref
 *
 * Managing BitRef levels
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

#ifndef __BITREFMANAGER_HPP__
#define __BITREFMANAGER_HPP__

#include "UnknownInterface.hpp"

namespace MoFEM {

  static const MOFEMuuid IDD_MOFEMBitRefManager = MOFEMuuid( BitIntefaceId(BITREFMANAGER_INTERFACE) );

  /**
   * \brief Managing BitRefLevels
   * \ingroup mofem_bit_ref
   * \nosubgrouping
   */
  struct BitRefManager : public UnknownInterface {

    PetscErrorCode query_interface(const MOFEMuuid &uuid,
                                   UnknownInterface **iface) const;

    MoFEM::Core &cOre;
    bool dEbug;

    BitRefManager(const MoFEM::Core &core);

    /**
    * \brief Destructor
    */
    ~BitRefManager();

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
     * @return           error code
     */
    PetscErrorCode setBitRefLevel(const Range &ents, const BitRefLevel &bit,
                                  const bool only_tets = true,
                                  int verb = 0) const;

    PetscErrorCode setBitRefLevelByDim(const EntityHandle meshset,
                                       const int dim, const BitRefLevel &bit,
                                       int verb = 0) const;

    PetscErrorCode setBitRefLevelByType(const EntityHandle meshset,
                                        const EntityType type,
                                        const BitRefLevel &bit,
                                        int verb = 0) const;

    /** brief add meshset and set bit ref level
    * \ingroup mofem_bit_ref
     *
     * \param EntityHandle MeshSet
     * \param BitRefLevel bitLevel
     */
    PetscErrorCode setBitLevelToMeshset(const EntityHandle meshset,
                                        const BitRefLevel &bit,
                                        int verb = 0) const;

    /**
     * \brief add bit ref level to ref entity
     * \ingroup mofem_bit_ref
     * @param  ents range of entities
     * @param  bit  bit ref level
     * @param  verb verbosity level
     * @return      error code
     */
    PetscErrorCode addBitRefLevel(const Range &ents, const BitRefLevel &bit,
                                  int verb = 0) const;

    /**
     * \brief Set nth bith ref lecel
     * @param  ents entities to set bit ref level
     * @param  n    nth bit
     * @param  b    value to set
     * @return      error code
     */
    PetscErrorCode setNthBitRefLevel(const Range &ents, const int n,
                                     const bool b, int verb = 0) const;

    /**
     * \brief Set nth bith ref lecel
     * \ingroup mofem_bit_ref
     * @param  n    nth bit
     * @param  b    value to set
     * @return      error code
     */
    PetscErrorCode setNthBitRefLevel(const int n, const bool b,
                                     int verb = 0) const;

    /** \brief left shift bit ref level
      * \ingroup mofem_bit_ref
      * this results of deletion of entities on far left side
      * 
      * \note Not implemented
      */
    PetscErrorCode shiftLeftBitRef(const int shif,
                                   const BitRefLevel mask = BitRefLevel().set(),
                                   int verb = -1) const;

    /** \brief right shift bit ref level
      * \ingroup mofem_bit_ref
      */
    PetscErrorCode
    shiftRightBitRef(const int shift,
                     const BitRefLevel mask = BitRefLevel().set(),
                     int verb = -1) const;

    /**@}*/

    /** \name Entity hanlders by bit ref level */

    /**@{*/

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
    PetscErrorCode getEntitiesByTypeAndRefLevel(const BitRefLevel &bit,
                                                const BitRefLevel &mask,
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
    PetscErrorCode getEntitiesByTypeAndRefLevel(const BitRefLevel &bit,
                                                const BitRefLevel &mask,
                                                const EntityType type,
                                                Range &ents,
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
    PetscErrorCode getEntitiesByRefLevel(const BitRefLevel &bit,
                                         const BitRefLevel &mask,
                                         const EntityHandle meshset) const;

    /**\brief add all ents from ref level given by bit to meshset
     * \ingroup mofem_bit_ref
     *
     * \note Entities NOT have to be added to MoFEM database
     *
     * \param BitRefLevel bitLevel
     * \param BitRefLevel mask
     * \retval ents
     */
    PetscErrorCode getEntitiesByRefLevel(const BitRefLevel &bit,
                                         const BitRefLevel &mask,
                                         Range &ents) const;

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
    PetscErrorCode getEntitiesByParentType(const BitRefLevel &bit,
                                           const BitRefLevel &mask,
                                           const EntityType type,
                                           Range &ents) const;

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
    virtual PetscErrorCode
    getAdjacenciesEquality(const EntityHandle from_entiti,
                           const int to_dimension, Range &adj_entities) const;

    /** \brief Get the adjacencies associated with a entity to entities of a
     * specified dimension.
      * \ingroup mofem_bit_ref
      *
      * bit ref level of adjacent entities is any of bit ref level of adjacent
     * entities
      */
    virtual PetscErrorCode getAdjacenciesAny(const EntityHandle from_entiti,
                                             const int to_dimension,
                                             Range &adj_entities) const;

    /** \brief Get the adjacencies associated with a entity to entities of a
     * specified dimension.
      * \ingroup mofem_bit_ref
      *
      * bit ref level of adjacent entities is equal to bit ref level of adjacent
     * entities
      */
    virtual PetscErrorCode
    getAdjacencies(const Problem *problem_ptr,
                   const EntityHandle *from_entities, const int num_netities,
                   const int to_dimension, Range &adj_entities,
                   const int operation_type = moab::Interface::INTERSECT,
                   const int verb = 0) const;

    /** \brief Get the adjacencies associated with a entity to entities of a
     * specified dimension.
      * \ingroup mofem_bit_ref
      *
      * bit ref level of adjacent entities is equal to bit ref level of adjacent
     * entities
      */
    virtual PetscErrorCode
    getAdjacencies(const BitRefLevel &bit, const EntityHandle *from_entities,
                   const int num_netities, const int to_dimension,
                   Range &adj_entities,
                   const int operation_type = moab::Interface::INTERSECT,
                   const int verb = 0) const;

    /**@}*/

    /** \name Update meshsets and ranges by children */

    /**@{*/

    /** \brief Get child entities form meshset containing parent entities
      * \ingroup mofem_update_meshsets_and_ranges
      *
      * Search for refined entities of given type whose parent are entities in
      *the
      * parent meshset. It can be used for example to transfer information about
      * boundary conditions to refined mesh or split mesh by interface
      * elements. It is used by function refine_MESHSET, to update MESHSET
      *finite elements.
      *
      * \param parent meshset
      * \param child_bit refinement level
      * \param type of refined entity
      * \param child_type meshset where child entities are stored (if the child
      *meshset is set to be the parent meshset, the parent would be updated with
      *the refined entities)
      * \param recursive if true parent meshset is searched recursively
      *
     **/
    PetscErrorCode updateMeshsetByEntitiesChildren(const EntityHandle parent,
                                                   const BitRefLevel &child_bit,
                                                   const EntityHandle child,
                                                   EntityType child_type,
                                                   const bool recursive = false,
                                                   int verb = 0);

    /** \brief update fields meshesets by child entities
      * \ingroup mofem_update_meshsets_and_ranges
      *
      */
    PetscErrorCode
    updateFieldMeshsetByEntitiesChildren(const BitRefLevel &child_bit,
                                         int verb = 0);

    /** \brief update field mesheset by child entities
      * \ingroup mofem_update_meshsets_and_ranges
      */
    PetscErrorCode updateFieldMeshsetByEntitiesChildren(
        const std::string name, const BitRefLevel &child_bit, int verb = 0);

    /** \brief update finite element mesheset by child entities
     * \ingroup mofem_update_meshsets_and_ranges
     */
    PetscErrorCode updateFiniteElementMeshsetByEntitiesChildren(
        const std::string name, const BitRefLevel &child_bit,
        const EntityType fe_ent_type, int verb = 0);

    /**
     * \brief Update range by prents
     *
     * FIXME: NOT TESTED
     *
     * @param  parent parent range
     * @param  child  childeren range
     * @return        error code
     */
    PetscErrorCode updateRange(const Range &parent, Range &child);

    /**@}*/

    /** \name Writting files */

    /**@{*/

    /**
     * \brief write bit ref level to file
     * @param  bit       bit ref level
     * @param  mask      mask of bit ref level
     * @param  type      type of entity
     * @param  file_name file name (see moab documentation)
     * @param  file_type file type (see moab documentation)
     * @param  options   file options (see moab documentation)
     * @return           error code
     */
    PetscErrorCode
    writeBitLevelByType(const BitRefLevel &bit, const BitRefLevel &mask,
                        const EntityType type, const char *file_name,
                        const char *file_type, const char *options) const;

    /**@}*/
  };
}

#endif //__BITREFMANAGER_HPP__


/***************************************************************************//**
 * \defgroup mofem_bit_ref BitRefLevels manager
 * \brief Managing BitRefLevels
 *
 * \ingroup mofem
 ******************************************************************************/
