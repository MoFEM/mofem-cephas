/** 
 * \file PrismInterface.hpp
 * 
 * \brief MoFEM interface
 *
 * Insert prisms in the interface between two surfaces
 * 
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

/** 
 * \brief Create interface at given faces and add flat prisms in-between 
 * 
 * \ingroup mofem_prism_interface
*/
struct PrismInterface : public UnknownInterface {

  MoFEMErrorCode query_interface(const MOFEMuuid &uuid,
                                 UnknownInterface **iface) const;

  MoFEM::Core &cOre;
  PrismInterface(const MoFEM::Core &core);

  /// destructor
  ~PrismInterface() {}

  /** 
   * \brief Store tetrahedra from each side of the interface 
   * separately in two child meshsets of the parent meshset
   *
   * Additional third child meshset contains nodes which can be split
   * and skin edges
   *
   * \param msId Id of meshset
   * \param cubit_bc_type type of meshset (NODESET, SIDESET or BLOCKSET and
   * more)
   * \param mesh_bit_level interface is added on this bit level
   * \param recursive if true parent meshset is searched recursively
   * \param  verb verbosity level
   * 
   * \note if bit_level == BitRefLevel.set() then interface will be added 
   * on all bit levels
   * 
   */
  MoFEMErrorCode getSides(const int msId, const CubitBCType cubit_bc_type,
                          const BitRefLevel mesh_bit_level,
                          const bool recursive, int verb = QUIET);

  /**
   * \brief Store tetrahedra from each side of the interface
   * in two child meshsets of the parent meshset
   *
   * Additional third child meshset contains nodes which can be split
   * and skin edges
   *
   * \param sideset parent meshset with the surface
   * \param mesh_bit_level interface is added on this bit level
   * \param recursive if true parent meshset is searched recursively
   * \param  verb verbosity level
   *
   * \note if bit_level == BitRefLevel.set() then interface will be added
   * on all bit levels
   *
   * 1. Get tets adjacent to nodes of the interface meshset.
   * 2. Take skin faces from these tets and get edges from that skin.
   * 3. Take skin from triangles of the interface.
   * 4. Subtract edges of skin faces from skin of triangles in order to get
   * edges in the volume of the body, and not on the interface boundary.
   * 5. Iterate between all triangles of the interface and find adjacent tets
   * on each side of the interface
   *
   */
  MoFEMErrorCode getSides(const EntityHandle sideset,
                          const BitRefLevel mesh_bit_level,
                          const bool recursive, int verb = QUIET);

  /**
   * \brief Find triangles which have three nodes on internal surface skin
   *
   * Internal surface skin is a set of all edges on the boundary of a given
   * surface inside the body. This set of edges is also called the surface
   * front. If a triangle has three nodes on the surface front, none of these 
   * nodes can be split. Therefore, such a triangle cannot be split and
   * should be removed from the surface. 
   *
   * @param  sideset        meshset with surface
   * @param  mesh_bit_level bit ref level of the volume mesh
   * @param  recursive      if true search in sub-meshsets
   * @param  faces_with_three_nodes_on_front returned faces
   * @param  verb           verbosity level
   *
   * @return error code
   */
  MoFEMErrorCode findFacesWithThreeNodesOnInternalSurfaceSkin(
      const EntityHandle sideset, const BitRefLevel mesh_bit_level,
      const bool recursive, Range &faces_with_three_nodes_on_front,
      int verb = QUIET);

  /**
   * \brief Split nodes and other entities of tetrahedra on both sides
   * of the interface and add prism elements in-between
   *
   * \param meshset parent meshset to get entities from
   * \param bit bit ref level on which new entities will be stored
   * \param msId meshset ID imported from cubit
   * \param cubit_bc_type type of meshset (NODESET, SIDESET or BLOCKSET and
   *more)
   * \param add_interface_entities if true add prism elements at interface
   * \param recursive if true parent meshset is searched recursively
   *
   * \note Parent meshset must have three child meshsets: two with tetrahedra
   * from each side of the interface, third containing nodes which can be split
   * and skin edges
   *
   */
  MoFEMErrorCode splitSides(const EntityHandle meshset, const BitRefLevel &bit,
                            const int msId, const CubitBCType cubit_bc_type,
                            const bool add_interface_entities,
                            const bool recursive = false, int verb = QUIET);

  /**
   * \brief Split nodes and other entities of tetrahedra on both sides
   * of the interface and add prism elements in-between
   *
   * \param meshset parent meshset to get entities from
   * \param bit bit ref level on which new entities will be stored
   * \param msId meshset ID imported from cubit
   * \param cubit_bc_type type of meshset (NODESET, SIDESET or BLOCKSET and
   *more)
   * \param add_interface_entities if true add prism elements at interface
   * \param recursive if true parent meshset is searched recursively
   *
   * \note Parent meshset must have three child meshsets: two with tetrahedra
   * from each side of the interface, third containing nodes which can be split
   * and skin edges
   */
  MoFEMErrorCode splitSides(const EntityHandle meshset, const BitRefLevel &bit,
                            const EntityHandle sideset,
                            const bool add_interface_entities,
                            const bool recursive = false, int verb = QUIET);

  /**
   * \brief Split nodes and other entities of tetrahedrons in children sets and
   * add prism elements
   *
   * \param meshset parent meshset to get entities from
   * \param bit bit ref level on which new entities will be stored
   * \param inhered_from_bit_level inherit nodes and other entities form this
   * bit level
   * \param inhered_from_bit_level_mask corresponding mask
   * \param add_interface_entities if true add prism elements at interface
   * \param recursive if true parent meshset is searched recursively
   *
   * \note Parent meshset must have three child meshsets: two with tetrahedra
   * from each side of the interface, third containing nodes which can be split
   * and skin edges
   * \note inhered_from_bit_level needs to be specified for some meshsets
   * with interfaces. Some nodes on some refinement levels are dividing edges 
   * but not splitting faces. Inheriting those nodes will not split faces.
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
 * \brief Create interface between faces
 *
 * Create interface between given faces and add in-between prisms
 * if needed
 *
 * \ingroup mofem
 */

#endif // __PRISMINTERFACE_HPP__
