/** \file MeshRefinement.hpp
 * \brief Interface managing boundary and blockset data on meshsets
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

#ifndef __MESHREFINE_HPP__
#define __MESHREFINE_HPP__

#include "UnknownInterface.hpp"

namespace MoFEM {

static const MOFEMuuid IDD_MOFEMMeshRefine = MOFEMuuid( BitIntefaceId(MESH_REFINE) );

/** \brief Mesh refinement interface

  Currently this class is abstraction to Core interface. In future should be
  outsourced as independent interface.

  \bug Not working on partitioned meshes
  \bug Need to be implemented as a stand alone interface not as a part of core
  structure which should be only basic database
  \bug If outsourced, class member functions should follow name convention

  */
struct MeshRefinement: public UnknownInterface {

  PetscErrorCode queryInterface(const MOFEMuuid& uuid, UnknownInterface** iface);

  Tag th_RefType;
  Tag th_RefParentHandle;
  Tag th_RefBitLevel;
  Tag th_RefBitLevel_Mask;
  Tag th_RefBitEdge;

  MoFEM::Core& cOre;
  MeshRefinement(const MoFEM::Core &core);

  ///destructor
  ~MeshRefinement() {}

  /**
   * \brief make vertices in the middle of edges in meshset and add them to refinement levels defined by bit
   *
   * Takes entities fromm meshsets and queried recursively (get entities from meshsets in meshsets, usually have to be used for CUBIT meshset).
   * If meshset does not contain any edges, get entities in dimension 3 and get edge adjacencies.
   *
   * \param EntityHandle meshset
   * \param BitRefLevel bitLevel
   * \param recursive If true, meshsets containing meshsets are queried recursively.  Returns the contents of meshsets, but not the meshsets themselves if true.
   */
  virtual PetscErrorCode add_verices_in_the_middel_of_edges(
    const EntityHandle meshset,const BitRefLevel &bit,const bool recursive = false,int verb = -1
  );

  /**
   * \brief make vertices in the middle of edges in meshset and add them to Refinement levels defined by bit
   *
   * Takes entities from meshsets and queried recursively (get entities from meshsets in meshsets, usually have to be used for CUBIT meshset).
   * If meshset does not contain any edges, get entities in dimension 3 and get edge adjacencies.
   *
   * \param Range consisting edges for refine
   * \param BitRefLevel bitLevel
   * \param recursive If true, meshsets containing meshsets are queried recursively.  Returns the contents of meshsets, but not the meshsets themselves if true.
   */
  virtual PetscErrorCode add_verices_in_the_middel_of_edges(const Range &edges,const BitRefLevel &bit,int verb = 0);

  /**\brief refine TET in the meshset
   *
   * \param EntityHandle meshset
   * \param BitRefLevel bitLevel
   * \param If TRUE, interface elements would be refined too
   */
  virtual PetscErrorCode refine_TET(const EntityHandle meshset,const BitRefLevel &bit,const bool respect_interface = false,int verb = 0);

  /**\brief refine TET in the meshset
   *
   * \param Range of tets to refine
   * \param BitRefLevel bitLevel
   * \param If TRUE, interface elements would be refined too
   */
  virtual PetscErrorCode refine_TET(const Range &tets,const BitRefLevel &bit,const bool respect_interface = false,int verb = 0);

  /**\brief refine PRISM in the meshset
   *
   * \param EntityHandle meshset
   * \param BitRefLevel bitLevel
   */
  virtual PetscErrorCode refine_PRISM(const EntityHandle meshset,const BitRefLevel &bit,int verb = 0);

  /**\brief refinem meshset, i.e. add child of refined entities to meshset
   *
   * \param EntityHandle meshset where to save the child refined entities
   * \param BitRefLevel bitLevel
   * \param recursive If true, meshsets containing meshsets are queried recursively.  Returns the contents of meshsets, but not the meshsets themselves if true.
   */
  virtual PetscErrorCode refine_MESHSET(
    const EntityHandle meshset,const BitRefLevel &bit,const bool recursive = false,int verb = 0
  );

};

DEPRECATED typedef MeshRefinement MeshRefinment;

}

#endif // __MESHREFINE_HPP__
