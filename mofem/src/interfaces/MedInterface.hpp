/** \file MedInterface.hpp
 * \brief Med file interface interface
 *
 * Interface loading mesh and data on mesh directly to mofem & moab
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

#ifdef WITH_MED

#ifndef __MED_INTERFACE_HPP__
#define __MED_INTERFACE_HPP__

namespace MoFEM {

static const MOFEMuuid IDD_MOFEMMedInterface = MOFEMuuid( BitIntefaceId(MED_INTERFACE) );

/** \brief merge node from two bit levels
  * \ingroup mofem
  */
struct MedInterface: public UnknownInterface {

  PetscErrorCode queryInterface(const MOFEMuuid& uuid, UnknownInterface** iface);

  MoFEM::interface& mField;
  MedInterface(MoFEM::Interface& m_field):
  mField(m_field) {
  }

};

#endif //__MED_INTERFACE_HPP__
#endif //WITH_MED
