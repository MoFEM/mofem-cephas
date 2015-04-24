/** \file FieldUnknownInterface.hpp
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

#ifndef __MOFEMUNKNOWNFIELD_HPP__
#define __MOFEMUNKNOWNFIELD_HPP__

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
static const MOFEMuuid IDD_MOFEMUnknown = MOFEMuuid( BitIntefaceId(FIELD_UNKNOWNINTERFACE) );

/** \brief base class for all interface classes
  * \ingroup mofem
  */
struct FieldUnknownInterface {
  virtual PetscErrorCode queryInterface (const MOFEMuuid& uuid, FieldUnknownInterface** iface) = 0;
  virtual ~FieldUnknownInterface() {}
};

}

#endif // __MOFEMUNKNOWNFIELD_HPP__
 
