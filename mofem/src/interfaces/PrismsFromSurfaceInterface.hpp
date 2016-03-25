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
  PetscErrorCode createPrisms(const Range &ents,Range &prisms,int verb = -1);
  PetscErrorCode seedPrismsEntities(Range &prisms,const BitRefLevel &bit,int verb = -1);

};

}

#endif //__PRISMS_FORM_SURFACE_HPP__
