/** \file UnknownInterface.hpp
 * \brief MoFEM interface
 *
 * Low level data structures not used directly by user
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

#ifndef __MOFEMUNKNOWNINTERFACE_HPP__
#define __MOFEMUNKNOWNINTERFACE_HPP__

namespace MoFEM {

struct MOFEMuuid {

  MOFEMuuid() { memset( this, 0, sizeof(MOFEMuuid)); }
  MOFEMuuid(BitIntefaceId uuid) { uUId = uuid; }

  //! returns whether two uuid's are equal
  bool operator==(const MOFEMuuid& orig) const {
    return (uUId&orig.uUId) == orig.uUId;
  }

  //uuid
  BitIntefaceId uUId;

};

//! uuid for an unknown interface
//! this can be used to either return a default interface
//! or a NULL interface
static const MOFEMuuid IDD_MOFEMUnknown = MOFEMuuid( BitIntefaceId(UNKNOWNINTERFACE) );

/** \brief base class for all interface classes
  * \ingroup mofem
  */
struct UnknownInterface {
  virtual PetscErrorCode queryInterface (const MOFEMuuid& uuid, UnknownInterface** iface) = 0;
  virtual ~UnknownInterface() {}
};

}

#endif // __MOFEMUNKNOWNINTERFACE_HPP__
