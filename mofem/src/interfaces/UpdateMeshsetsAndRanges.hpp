/** \file UpdateMeshsetsAndRanges.hpp
 * \brief Updata meshsets and ranges based on parent child relation
 *
 * Update meshsets and ranges based on parten child relation. If range or
 * meshset which children is given, add them to meshset or range to which
 * parent belongs.
 *
 * \ingroup mofem_update_meshsets_and_ranges
 */

/* MoFEM is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
 * License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with MoFEM. If not, see <http://www.gnu.org/licenses/>
*/

#ifndef __UPDATE_MESHSETS_AND_RANGES_HPP__
#define __UPDATE_MESHSETS_AND_RANGES_HPP__

namespace MoFEM {

  static const MOFEMuuid IDD_MOFEMUpdateMeshsetsAndRanges = MOFEMuuid(BitIntefaceId(UPDATEMESHSETSANDRANGES_INTERFACE) );

  /**
   * \brief Interface for updating child entities based on their relation to parent
   * \ingroup mofem_update_meshsets_and_ranges
   *
   */
  struct UpdateMeshsetsAndRanges: public UnknownInterface {

    PetscErrorCode queryInterface(const MOFEMuuid& uuid, UnknownInterface** iface);

    MoFEM::Core& cOre;
    UpdateMeshsetsAndRanges(const MoFEM::Core& core):
    cOre(const_cast<MoFEM::Core&>(core)) {
    }

    /** \brief Get child entities form meshset containing parent entities
      * \ingroup mofem_update_meshsets_and_ranges
      *
      * Search for refined entities of given type whose parent are entities in the
      * parent meshset. It can be used for example to transfer information about
      * boundary conditions to refined mesh or split mesh by interface
      * elements. It is used by function refine_MESHSET, to update MESHSET finite elements.
      *
      * \param parent meshset
      * \param child_bit refinement level
      * \param type of refined entity
      * \param child_type meshset where child entities are stored (if the child meshset is set to be the parent meshset, the parent would be updated with the refined entities)
      * \param recursive if true parent meshset is searched recursively
      *
     **/
    PetscErrorCode updateMeshsetByEntitiesChildren(
      const EntityHandle parent, const BitRefLevel &child_bit,const EntityHandle child, EntityType child_type,
      const bool recursive = false, int verb = 0
    );

    /** \brief update fields meshesets by child entities
      * \ingroup mofem_update_meshsets_and_ranges
      *
      */
    PetscErrorCode updateFieldMeshsetByEntitiesChildren(const BitRefLevel &child_bit,int verb = 0);

    /** \brief update field mesheset by child entities
      * \ingroup mofem_update_meshsets_and_ranges
      */
    PetscErrorCode updateFieldMeshsetByEntitiesChildren(const std::string name,const BitRefLevel &child_bit,int verb = 0);

    /** \brief update finite element mesheset by child entities
     * \ingroup mofem_update_meshsets_and_ranges
     */
    PetscErrorCode updateFiniteElementMeshsetByEntitiesChildren(
      const std::string name,const BitRefLevel &child_bit,const EntityType fe_ent_type,int verb = 0
    );

    /**
     * \brief Update range by prents
     *
     * FIXME: NOT TESTED
     *
     * @param  parent parent range
     * @param  child  childeren range
     * @return        error code
     */
    PetscErrorCode updateRange(const Range& parent,Range& child);

  };

}

#endif // __UPDATE_MESHSETS_AND_RANGES_HPP__


/***************************************************************************//**
 * \defgroup mofem_update_meshsets_and_ranges Updata meshsets and ranges based on child relation
 * \brief Interface for UpdateMeshsetsAndRanges
 *
 * \ingroup mofem
 ******************************************************************************/
