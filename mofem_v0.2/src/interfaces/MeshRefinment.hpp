/** \file MeshRefinment.hpp
 * \brief MoFEM interface 
 * 
 * Low level data structures not used directly by user
 *
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

#include "FieldUnknownInterface.hpp"

namespace MoFEM {

static const MOFEMuuid IDD_MOFEMMeshRefine = MOFEMuuid( BitIntefaceId(MESH_REFINE) );

/** \brief Mesh refinement interface

  Currently this class is abstraction to Core interface. In future should be
  outsourced as independent interface.
  
  \bug not working on partitioned mesha

  \bug need to be implemented as a stand alone interface not as a part of core
  structure which should be only basic database

  \bug if outsourced, class member functions should follow name convention

  */
struct MeshRefinment: public FieldUnknownInterface {

  ///destructor
  virtual ~MeshRefinment() {}

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
    const EntityHandle meshset,const BitRefLevel &bit,const bool recursive = false,int verb = -1) = 0;

  /**
   * \brief make vertices in the middle of edges in meshset and add them to refinment levels defined by bit
   *
   * Takes entities from meshsets and queried recursively (get entities from meshsets in meshsets, usually have to be used for CUBIT meshset).
   * If meshset does not contain any edges, get entities in dimension 3 and get edge adjacencies.
   *
   * \param Range consisting edges for refine 
   * \param BitRefLevel bitLevel
   * \param recursive If true, meshsets containing meshsets are queried recursively.  Returns the contents of meshsets, but not the meshsets themselves if true.
   */
  virtual PetscErrorCode add_verices_in_the_middel_of_edges(const Range &edges,const BitRefLevel &bit,int verb = -1) = 0;

  /**\brief refine TET in the meshset
   *
   * \param EntityHandle meshset
   * \param BitRefLevel bitLevel
   * \param If TRUE, interface elements would be refined too
   */
  virtual PetscErrorCode refine_TET(const EntityHandle meshset,const BitRefLevel &bit,const bool respect_interface = false) = 0;

  /**\brief refine TET in the meshset
   *
   * \param Range of tets to refine
   * \param BitRefLevel bitLevel
   * \param If TRUE, interface elements would be refined too
   */
  virtual PetscErrorCode refine_TET(const Range &tets,const BitRefLevel &bit,const bool respect_interface = false) = 0;

  /**\brief refine PRISM in the meshset
   *
   * \param EntityHandle meshset
   * \param BitRefLevel bitLevel
   */
  virtual PetscErrorCode refine_PRISM(const EntityHandle meshset,const BitRefLevel &bit,int verb = -1) = 0;

  /**\brief refinem meshset, i.e. add child of refined entities to meshset
   *
   * \param EntityHandle meshset where to save the child refined entities
   * \param BitRefLevel bitLevel
   * \param recursive If true, meshsets containing meshsets are queried recursively.  Returns the contents of meshsets, but not the meshsets themselves if true.
   */
  virtual PetscErrorCode refine_MESHSET(const EntityHandle meshset,const BitRefLevel &bit,
    const bool recursive = false,int verb = -1) = 0;

};

}

#endif // __MESHREFINE_HPP__




