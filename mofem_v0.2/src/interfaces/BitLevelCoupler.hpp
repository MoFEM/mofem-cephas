/** \file BitLevelCoupler.hpp
 * \brief BitLevelCoupler interface 

 * Is used to couple bit levels to enable easy and efficient projection between
 * levels. It is not assumed that print children relation between entities,
 * however if such relation exist is used coupling algorithm.

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

/** \brief manage adjacencies between mesh bit levels
  * \ingroup mofem
  *
  */
struct BitLevelCouplerInterface: public FieldUnknownInterface {

  PetscErrorCode queryInterface(const MOFEMuuid& uuid, FieldUnknownInterface** iface);

  MoFEM::Core& cOre;
  BitLevelCouplerInterface(MoFEM::Core& core): cOre(core) {};

  /** \brief finding adjacencies between vertices and tetrahedrons. 
    *
    * Use boundary volume tree to find tetrahedral or other volume element. 
    */
  PetscErrorCode buidlAdjacenciesVerticesOnTets(const BitRefLevel &parent_level,Range &children,
    const double iter_tol = 1.0e-10,
    const double inside_tol = 1.0e-6,
    bool vertex_elements = false,int verb = 0);

  /** \brief finding adjacencies between vertices and faces, edges or other vertices
    *
    * Assumes that adjacencies between vertices and tetrahedrons are known 
    */
  PetscErrorCode buidlAdjacenciesVerticesOnFacesEdgesVolumes(
    const BitRefLevel &parent_level,Range &children,bool vertex_elements = true,const double inside_tol = 1e-6,
    bool throw_error = true,int verb = 0);

  /** \brief build adjacencies for edges, faces and volumes
    *
    * Assumes that adjacencies for vertices has been build
    */
  PetscErrorCode buidlAdjacenciesEdgesFacesVolumes(
    const BitRefLevel &parent_level,Range &children,bool elements = true,int verb = 0);


  private:

  PetscErrorCode chanegParent(RefMoFEMEntity_multiIndex::iterator it,EntityHandle parent,bool element);
  PetscErrorCode verifyParent(RefMoFEMEntity_multiIndex::iterator it,EntityHandle parent);



};

}

#endif //__BITLEVELCOUPLER_HPP__
