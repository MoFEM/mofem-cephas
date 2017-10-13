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

static const MOFEMuuid IDD_MOFEMPrismInterface = MOFEMuuid( BitIntefaceId(PRISM_INTEFACE) );

/** \brief Make interface on given faces and create flat prism in that space

  \todo FIXME Names of methods do not follow naming convention and are difficult
  to work with.

  \ingroup mofem_prism_interface

*/
struct PrismInterface: public UnknownInterface {

  PetscErrorCode queryInterface(const MOFEMuuid& uuid, UnknownInterface** iface);

  MoFEM::Core& cOre;
  PrismInterface(const MoFEM::Core &core);

  ///destructor
  ~PrismInterface() {}


  /** \brief create two children meshsets in the meshset containing tetrahedral on two sides of faces
    *
    * \param msId Id of meshset
    * \param cubit_bc_type type of meshset (NODESET, SIDESET or BLOCKSET and more)
    * \param mesh_bit_level add interface on bit level is bit_level = BitRefLevel.set() then add interface on all bit levels
    * \param recursive if true parent meshset is searched recursively
    */
  PetscErrorCode getSides(
    const int msId,
    const CubitBCType cubit_bc_type,
    const BitRefLevel mesh_bit_level,
    const bool recursive,int verb = -1
  );


  /** \brief create two children meshsets in the meshset containing tetrahedral on two sides of faces
   *
   * Get tets adj to faces. Take skin form tets and get edges from that skin. Take skin form triangles (the face).
   * Subtrac skin faces edges form skin edges in order to get eges on the boundary of the face which is in the
   * voulume of the body, but is not on the boundary.
   * Each child set has a child containing nodes which can be split and skin edges.
   * After that simply iterate under all tets on one side which are adjacent to the face are found.
   * Side tets are stored in to children meshsets of the SIDESET meshset.
   */
  PetscErrorCode getSides(
    const EntityHandle sideset,
    const BitRefLevel mesh_bit_level,
    const bool recursive,int verb = -1
  );

  /**
   * \brief Find if tringle has three nodes on internal surface skin
   *
   * Internal surface skin is a set of edges in interia of the body on boundary
   * of surface. This set of edges is called surface front. If surface face has three nodes on
   * surface front, non of the face nodes is split and should be removed from surface
   * if it is going to be split.
   *
   * @param  sideset        meshset with surface
   * @param  mesh_bit_level bit ref level of the volume mesh
   * @param  recursive      search in sub-meshsets
   * @param  faces_with_three_nodes_on_front returned faces
   * @param  verb           error code
   *
   * @return error code
   */
  PetscErrorCode findIfTringleHasThreeNodesOnInternalSurfaceSkin(
    const EntityHandle sideset,
    const BitRefLevel mesh_bit_level,
    const bool recursive,
    Range& faces_with_three_nodes_on_front,
    int verb = -1
  );

  /**
   * \brief split nodes and other entities of tetrahedral in children sets and add prism elements
   *
   * The all new entities (prisms, tets) are added to refinement level given by bit
   * \param meshset meshset to get entities from
   * \param BitRefLevel new level where refinement would be stored
   * \param msId meshset ID imported from cubit
   * \param cubit_bc_type type of meshset (NODESET, SIDESET or BLOCKSET and more)
   * \param add_intefece_entities meshset which contain the interface
   * \param recursive if true parent meshset is searched recursively
   *
   * Each inteface face has two tages,
   *    const int def_side[] = {0};
   *	rval = moab.tag_get_handle("INTERFACE_SIDE",1,MB_TYPE_INTEGER,
   *	  th_interface_side,MB_TAG_CREAT|MB_TAG_SPARSE,def_side); CHKERRQ_MOAB(rval);
   *
   *  	const EntityHandle def_node[] = {0};
   *	rval = moab.tag_get_handle("SIDE_INTFACE_ELEMENT",1,MB_TYPE_HANDLE,
   *	  th_side_elem,MB_TAG_CREAT|MB_TAG_SPARSE,def_node); CHKERRQ_MOAB(rval);
   *
   * First tag inform abot inteface side, second tag inform about side adjacent
   * inteface element.
   *
   */
  PetscErrorCode splitSides(
    const EntityHandle meshset,const BitRefLevel &bit,
    const int msId,const CubitBCType cubit_bc_type,
    const bool add_iterfece_entities,const bool recursive = false,int verb = -1
  );

  /**
   * \brief split nodes and other entities of tetrahedral in children sets and add prism elements
   *
   * The all new entities (prisms, tets) are added to refinement level given by bit
   */
  PetscErrorCode splitSides(
    const EntityHandle meshset,const BitRefLevel &bit,
    const EntityHandle sideset,const bool add_iterfece_entities,
    const bool recursive = false,int verb = -1
  );


  /**
   * \brief split nodes and other entities of tetrahedrons in children sets and add prism elements
   *
   * The all new entities (prisms, tets) are added to refinement level given by bit
   *
   * \param meshset
   * \param Refinement bit level of new mesh
   * \param inheret_from_bit_level inheret nodes and other entities form this bit level.
   * \param add_iterfece_entities add prism elements at interface
   * \param recuslsive do meshesets in the meshset
   *
   * note inheret_from_bit_level is need to be specified to some meshsets
   * with interfaces. Some nodes on some refinement levels dividing edges but
   * not splitting faces. Inheriting those nodes will not split faces.
   *
   */
  PetscErrorCode splitSides(
    const EntityHandle meshset,const BitRefLevel &bit,
    const BitRefLevel &inheret_from_bit_level,const BitRefLevel &inheret_from_bit_level_mask,
    const EntityHandle sideset,const bool add_iterfece_entities,const bool recursive = false,int verb = -1
  );

};

}

/***************************************************************************//**
 * \defgroup mofem_prism_interface Prism interface
 * \brief Make interface between faces
 *
 * Make interface between faces (surface) and put in between prism element
 *
 * \ingroup mofem
 ******************************************************************************/


#endif // __PRISMINTERFACE_HPP__
