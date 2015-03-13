/** \file BitLevelCoupler.hpp
 * \brief BitLevelCoupler interface 
 * 

 * Is used to couple bit levels to enable easy and efficient projection between
 * levels. It is not assumed that print children relation between entities,
 * however if such relation exist is used coupling algorithm.


 *
 * MoFEM is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
 * License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with MoFEM. If not, see <http://www.gnu.org/licenses/>
*/

#ifndef __BITLEVELCOUPLER_HPP__
#define __BITLEVELCOUPLER_HPP__

namespace MoFEM {

static const MOFEMuuid IDD_MOFENBitLevelCoupler = MOFEMuuid( BitIntefaceId(BITLEVELCOUPLER_INTERFACE) );

struct BitLevelCouplerInterface: public FieldUnknownInterface {

  PetscErrorCode queryInterface(const MOFEMuuid& uuid, FieldUnknownInterface** iface);

  MoFEM::Core& cOre;
  BitLevelCouplerInterface(MoFEM::Core& core): cOre(core) {};

  PetscErrorCode buidlDirectAdjacencies(const BitRefLevel &parent_level,Range &children);
  PetscErrorCode buidlDirectAdjacencies(const BitRefLevel &parent_level,const BitRefLevel &children_level);


};

}

#endif //__BITLEVELCOUPLER_HPP__
