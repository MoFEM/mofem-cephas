/** \file PrismsFromSurface.hpp
 * \brief PrismsFromSurface interface
 *
 * Create prisms from surface triangle elements
 */

/* MoFEM is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
 * License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with MoFEM. If not, see <http://www.gnu.org/licenses/>
*/

#ifndef __PRISMS_FORM_SURFACE_HPP__
#define __PRISMS_FORM_SURFACE_HPP__

namespace MoFEM {

static const MOFEMuuid IDD_MOFEMPrismsFromSurface = MOFEMuuid(BitIntefaceId(PRISMSFROMSURFACE_INTERFACE));

/** \brief merge node from two bit levels
  * \ingroup mofem
  */
struct PrismsFromSurfaceInterface: public UnknownInterface {

  PetscErrorCode queryInterface(const MOFEMuuid& uuid, UnknownInterface** iface);

  MoFEM::Core& cOre;
  PrismsFromSurfaceInterface(MoFEM::Core& core): cOre(core) {};

  map<EntityHandle,EntityHandle> createdVertices;

  /**
   * \brief Make prisms from triangles
   * @param  ents   Range of triangles
   * @param  prisms Returned range of prisms
   * @param  verb   Verbosity level
   * @return        Error code
   */
  PetscErrorCode createPrisms(const Range &ents,Range &prisms,int verb = -1);

  /**
   * \brief Seed prism entities by bit level
   * @param  prisms Range of entities
   * @param  bit    BitRefLevel
   * @param  verb   Verbosity level
   * @return        Error code
   */
  PetscErrorCode seedPrismsEntities(Range &prisms,const BitRefLevel &bit,int verb = -1);

  /**
   * \brief Make prisms by extruding top or bottom prisms
   * @param  prisms      Input prisms
   * @param  from_down  Use top or down face, if true from f3
   * @param  out_prisms  Returned prisms entities
   * @param  verb        Verbosity level
   * @return             Error code
   */
  PetscErrorCode createPrismsFromPrisms(const Range &prisms,bool from_down,Range &out_prisms,int verb = -1);


  /**
   * Set uniform thickness
   * @param  prisms   Range of prisms
   * @param  director3 Displacement of face 3
   * @param  director4 Displacement of face 4
   * @return
   */
  PetscErrorCode setThickness(const Range &prisms,const double director3[],const double director4[]);

};

}

#endif //__PRISMS_FORM_SURFACE_HPP__
