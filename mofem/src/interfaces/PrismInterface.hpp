/** \file PrismInterface.hpp
 * \brief MoFEM interface
 *
 * Low level data structures not used directly by user

 * \ingroup mofem_prism_interface
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

#ifndef __PRISMINTERFACE_HPP__
#define __PRISMINTERFACE_HPP__

#include "UnknownInterface.hpp"

namespace MoFEM {

static const MOFEMuuid IDD_MOFEMPrismInterface =
    MOFEMuuid(BitIntefaceId(PRISM_INTEFACE));

/** \brief Make interface on given faces and create flat prism in that space

  \todo FIXME Names of methods do not follow naming convention and are difficult
  to work with.

  \ingroup mofem_prism_interface

*/
struct PrismInterface : public UnknownInterface {

  MoFEMErrorCode query_interface(const MOFEMuuid &uuid,
                                 UnknownInterface **iface) const;

  MoFEM::Core &cOre;
  PrismInterface(const MoFEM::Core &core);

  /// destructor
  ~PrismInterface() {}

  /** \brief create two children meshsets in the meshset containing tetrahedral
   * on two sides of faces
   *
   * \param msId Id of meshset
   * \param cubit_bc_type type of meshset (NODESET, SIDESET or BLOCKSET and
   * more) \param mesh_bit_level add interface on bit level is bit_level =
   * BitRefLevel.set() then add interface on all bit levels \param recursive if
   * true parent meshset is searched recursively
   */
  MoFEMErrorCode getSides(const int msId, const CubitBCType cubit_bc_type,
                          const BitRefLevel mesh_bit_level,
                          const bool recursive, int verb = QUIET);

  /** \brief create two children meshsets in the meshset containing tetrahedral
   * on two sides of faces
   *
   * Get tets adj to faces. Take skin form tets and get edges from that skin.
   * Take skin form triangles (the face). Subtrac skin faces edges form skin
   * edges in order to get edges on the boundary of the face which is in the
   * volume of the body, but is not on the boundary.
   * Each child set has a child containing nodes which can be split and skin
   * edges. After that simply iterate under all tets on one side which are
   * adjacent to the face are found. Side tets are stored in to children
   * meshsets of the SIDESET meshset.
   */
  MoFEMErrorCode getSides(const EntityHandle sideset,
                          const BitRefLevel mesh_bit_level,
                          const bool recursive, int verb = QUIET);

  /** 
   * \brief Find if triangle has three nodes on internal surface skin
   * 
   * Internal surface skin is a set of edges in the body on the boundary of the
   * surface. This set of edges is called a surface front. If the surface the
   * face has three nodes on the surface front, none of the face nodes is split
   * and should be removed from the surface. Otherwise, such a triangle cannot
   * be split.
   *
   * @param  sideset        meshset with surface
   * @param  mesh_bit_level bit ref level of the volume mesh
   * @param  recursive      search in sub-meshsets
   * @param  faces_with_three_nodes_on_front returned faces
   * @param  verb           verbosity level
   *
   * @return error code
   */
  MoFEMErrorCode findIfTringleHasThreeNodesOnInternalSurfaceSkin(
      const EntityHandle sideset, const BitRefLevel mesh_bit_level,
      const bool recursive, Range &faces_with_three_nodes_on_front,
      int verb = QUIET);

  /**
   * \brief split nodes and other entities of tetrahedral in children sets and
   *add prism elements
   *
   * The all new entities (prisms, tets) are added to refinement level given by
   *bit \param meshset meshset to get entities from \param BitRefLevel new level
   *where refinement would be stored \param msId meshset ID imported from cubit
   * \param cubit_bc_type type of meshset (NODESET, SIDESET or BLOCKSET and
   *more) \param add_interface_entities meshset which contain the interface
   * \param recursive if true parent meshset is searched recursively
   *
   * Each inteface face has two tags,
   *    const int def_side[] = {0};
   *	rval = moab.tag_get_handle("INTERFACE_SIDE",1,MB_TYPE_INTEGER,
   *	  th_interface_side,MB_TAG_CREAT|MB_TAG_SPARSE,def_side);
   *CHKERRQ_MOAB(rval);
   *
   *  	const EntityHandle def_node[] = {0};
   *	rval = moab.tag_get_handle("SIDE_INTFACE_ELEMENT",1,MB_TYPE_HANDLE,
   *	  th_side_elem,MB_TAG_CREAT|MB_TAG_SPARSE,def_node); CHKERRQ_MOAB(rval);
   *
   * First tag inform about inteface side, second tag inform about side adjacent
   * inteface element.
   *
   */
  MoFEMErrorCode splitSides(const EntityHandle meshset, const BitRefLevel &bit,
                            const int msId, const CubitBCType cubit_bc_type,
                            const bool add_interface_entities,
                            const bool recursive = false, int verb = QUIET);

  /**
   * \brief split nodes and other entities of tetrahedral in children sets and
   * add prism elements
   *
   * The all new entities (prisms, tets) are added to refinement level given by
   * bit
   */
  MoFEMErrorCode splitSides(const EntityHandle meshset, const BitRefLevel &bit,
                            const EntityHandle,
                            const bool add_interface_entities,
                            const bool recursive = false, int verb = QUIET);

  /**
   * \brief split nodes and other entities of tetrahedrons in children sets and
   * add prism elements
   *
   * The all new entities (prisms, tets) are added to refinement level given by
   * bit
   *
   * \param meshset
   * \param Refinement bit level of new mesh
   * \param inhered_from_bit_level inhered nodes and other entities form this
   * bit level. \param add_interface_entities add prism elements at interface
   * \param recursive do meshesets in the meshset
   *
   * note inhered_from_bit_level is need to be specified to some meshsets
   * with interfaces. Some nodes on some refinement levels dividing edges but
   * not splitting faces. Inheriting those nodes will not split faces.
   *
   */
  MoFEMErrorCode splitSides(const EntityHandle meshset, const BitRefLevel &bit,
                            const BitRefLevel &inhered_from_bit_level,
                            const BitRefLevel &inhered_from_bit_level_mask,
                            const EntityHandle sideset,
                            const bool add_interface_entities,
                            const bool recursive = false, int verb = QUIET);
};

} // namespace MoFEM

/**
 * \defgroup mofem_prism_interface PrismInterface
 * \brief Make interface between faces
 *
 * Make interface between faces (surface) and put in between prism element if
 *needed.
 *
 * \ingroup mofem
 */

#endif // __PRISMINTERFACE_HPP__
