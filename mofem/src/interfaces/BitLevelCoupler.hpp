/** \file BitLevelCoupler.hpp
 * \brief BitLevelCoupler interface 

 * Is used to couple bit levels to enable easy and efficient projection between
 * levels. It is not assumed that print children relation between entities,
 * however if such relation exist is used coupling algorithm.
 */

/* MoFEM is distributed in the hope that it will be useful, but WITHOUT
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

static const MOFEMuuid IDD_MOFEMBitLevelCoupler = MOFEMuuid( BitIntefaceId(BITLEVELCOUPLER_INTERFACE) );

/** \brief Interface set parent for verrtices, edges, triangles and tetrahedrons. 
  * \ingroup mofem
  *
  */
struct BitLevelCouplerInterface: public FieldUnknownInterface {

  PetscErrorCode queryInterface(const MOFEMuuid& uuid, FieldUnknownInterface** iface);

  MoFEM::Core& cOre;
  bool vErify;	///< by defualt is switched off, swith it on to verify if existing parent is equal to parent set by interface

  BitLevelCouplerInterface(MoFEM::Core& core): cOre(core),vErify(false) {};

  /** \brief build adaptive kd-tree
    */
  PetscErrorCode buildTree(const BitRefLevel &parent_level,int verb = 0);

  /** \brief reset adaptive kd-tree
    */
  PetscErrorCode resetTree(const BitRefLevel &parent_level,int verb = 0);

  /** \brief get parent entity

    * Use kd-tree to find tetrahedral or other volume element. 
  
    \param coordinate
    \param parent returned parent entity    
    \param iter_tol tolerance for convergence of point search
    \param inside_tol tolerance for inside element calculation
    \param throw_error if parent can not be found
    \param verbose level 

    */
  PetscErrorCode getParent(const double *coords,EntityHandle &parent,
    bool tet_only = false,const double iter_tol = 1.0e-10,const double inside_tol = 1.0e-6,int verb = 0);

  /** \brief finding parents for vertices 
    *
    * Use kd-tree to find tetrahedral or other volume element. 
  
    \param parent_level bit level of parents

    \param children list of vertices for which parents are being set

    \param vertex_elements if true algorithm assumes that vertices elements are
    used. IF NOT SET AND SUCH ELEMENTS EXIST IT WILL RESULT IN UNPREDICTABLE
    BEHAVIOUR.

    \param iter_tol tolerance for convergence of point search
    
    \param inside_tol tolerance for inside element calculation

    \param throw_error if parent can not be found
  
    \param verbose level 
   
    */
  PetscErrorCode buidlAdjacenciesVerticesOnTets(const BitRefLevel &parent_level,Range &children,
    bool vertex_elements = false,
    const double iter_tol = 1.0e-10,
    const double inside_tol = 1.0e-6,
    bool throw_error = true,
    int verb = 0);

  /** \brief finding parents for edegs, faces and tets

    It assumes that parents for vertices are known. Run
    buidlAdjacenciesVerticesOnTets if parents for vertices are not set.

    \param parent_level bit level of parents

    \param children list of entities for which parents are being set

    \param vertex_elements if true algorithm assumes that vertices elements are
    used. IF NOT SET AND SUCH ELEMENTS EXIST IT WILL RESULT IN UNPREDICTABLE
    BEHAVIOUR.

    \param iter_tol tolerance for convergence of point search
    
    \param inside_tol tolerance for inside element calculation

    \param throw_error if parent can not be found
  
    \param verbose level 

    */
  PetscErrorCode buidlAdjacenciesEdgesFacesVolumes(
    const BitRefLevel &parent_level,Range &children,bool elements = true,int verb = 0);

  /** \brief reset parent entities

    This is needed for testing.
  
    */
  PetscErrorCode resetParents(Range &children,bool elements = true,int verb = 0);

  private:

  PetscErrorCode chanegParent(RefMoFEMEntity_multiIndex::iterator it,EntityHandle parent,bool element);
  PetscErrorCode verifyParent(RefMoFEMEntity_multiIndex::iterator it,EntityHandle parent);

  double cOords[12+3];
  double diffN[12],N[4];
  double locCoords[3];
  const EntityHandle *cOnn;

  PetscErrorCode getLocCoordsOnTet(EntityHandle tet,const double *glob_coords,int verb = 0);

  boost::scoped_ptr<AdaptiveKDTree> treePtr;

};

}

#endif //__BITLEVELCOUPLER_HPP__
