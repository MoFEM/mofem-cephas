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

static const MOFEMuuid IDD_MOFEMBitLevelCoupler =
    MOFEMuuid(BitIntefaceId(BITLEVELCOUPLER_INTERFACE));

/** \brief Interface set parent for vertices, edges, triangles and tetrahedrons.
 * \ingroup mofem
 *
 * FIXME: Not tested, slow, bugs
 *
 */
struct BitLevelCoupler : public UnknownInterface {

  MoFEMErrorCode query_interface(const MOFEMuuid &uuid,
                                 UnknownInterface **iface) const;

  MoFEM::Core &cOre;
  bool vErify; ///< by default is switched off, with it on to verify if existing
               ///< parent is equal to parent set by interface

  BitLevelCoupler(const MoFEM::Core &core)
      : cOre(const_cast<MoFEM::Core &>(core)), vErify(false) {}

  /** \brief build adaptive kd-tree
   */
  MoFEMErrorCode buildTree(const BitRefLevel &parent_level, int verb = 0);

  /** \brief reset adaptive kd-tree
   */
  MoFEMErrorCode resetTree(const BitRefLevel &parent_level, int verb = 0);

  /** \brief get parent entity

    * Use kd-tree to find tetrahedral or other volume element.

    \param coordinate
    \param parent returned parent entity
    \param iter_tol tolerance for convergence of point search
    \param inside_tol tolerance for inside element calculation
    \param throw_error if parent can not be found
    \param verbose level

    */
  MoFEMErrorCode getParent(const double *coords, EntityHandle &parent,
                           bool tet_only = false,
                           const double iter_tol = 1.0e-10,
                           const double inside_tol = 1.0e-6, int verb = 0);

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
  MoFEMErrorCode buildAdjacenciesVerticesOnTets(
      const BitRefLevel &parent_level, Range &children,
      bool vertex_elements = false, const double iter_tol = 1.0e-10,
      const double inside_tol = 1.0e-6, bool throw_error = true, int verb = 0);

  /** \brief finding parents for edegs, faces and tets

    It assumes that parents for vertices are known. Run
    buildAdjacenciesVerticesOnTets if parents for vertices are not set.

    \param parent_level bit level of parents

    \param children list of entities for which parents are being set

    \param vertex_elements if true algorithm assumes that vertices elements are
    used. IF NOT SET AND SUCH ELEMENTS EXIST IT WILL RESULT IN UNPREDICTABLE
    BEHAVIOR.

    \param iter_tol tolerance for convergence of point search

    \param inside_tol tolerance for inside element calculation

    \param throw_error if parent can not be found

    \param verbose level

    */
  MoFEMErrorCode
  buildAdjacenciesEdgesFacesVolumes(const BitRefLevel &parent_level,
                                    Range &children, bool elements = true,
                                    int verb = 0);

  /** \brief reset parent entities

    This is needed for testing.

    */
  MoFEMErrorCode resetParents(Range &children, bool elements = true,
                              int verb = 0);

  /**
   * \brief copy data from parents
   *
   * This not approximate date from, simply copy DOFs values from one mesh to
   * another. This is useful for special case of refinement, e.g. insertion of
   * interface, where entities on the interface are doubled to create
   * displacement jump. In general case use this function could lead to wrong
   * results.
   *
   * \note Move data only if type of child and parent is the same.
   *
   * @param  parent pointer to array of parent entities
   * @param  children pointer to array of children entities
   * @param  verify if true verifi consistency with database
   * @return     error code
   */
  MoFEMErrorCode
  copyFieldDataFromParentToChildren(const std::vector<EntityHandle> &parents,
                                    const std::vector<EntityHandle> &children,
                                    const bool verify = true);

  /**
   * \brief copy data from parents
   *
   * This not approximate date from, simply copy DOFs values from one mesh to
   * another. This is useful for special case of refinement, e.g. insertion of
   * interface, where entities on the interface are doubled to create
   * displacement jump.
   *
   * In general case use this function could lead to wrong results.
   *
   *
   * \note Move data only if type of child and parent is the same.
   *
   * @param  bit bit ref level
   * @param  verify if true verifi consistency with database
   * @return     error code
   */
  MoFEMErrorCode copyFieldDataFromParentToChildren(const BitRefLevel bit,
                                                   const BitRefLevel mask,
                                                   const bool verify = true);

private:
  MoFEMErrorCode chanegParent(RefEntity_multiIndex::iterator it,
                              EntityHandle parent);
  MoFEMErrorCode verifyParent(RefEntity_multiIndex::iterator it,
                              EntityHandle parent);

  double cOords[12 + 3];
  double diffN[12], N[4];
  double locCoords[3];
  const EntityHandle *cOnn;

  MoFEMErrorCode getLocCoordsOnTet(EntityHandle tet, const double *glob_coords,
                                   int verb = 0);

  boost::scoped_ptr<AdaptiveKDTree> treePtr;
};

} // namespace MoFEM

#endif //__BITLEVELCOUPLER_HPP__
