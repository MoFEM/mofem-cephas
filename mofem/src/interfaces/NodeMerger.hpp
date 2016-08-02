/** \file NodeMerger.hpp
 * \brief NodeMerger interface
 *
 * Node merger interface
 */

/* MoFEM is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
 * License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with MoFEM. If not, see <http://www.gnu.org/licenses/>
*/

#ifndef __NODE_MERGER_HPP__
#define __NODE_MERGER_HPP__

namespace MoFEM {


static const MOFEMuuid IDD_MOFEMNodeMerger = MOFEMuuid( BitIntefaceId(NODEMERGER_INTERFACE) );

/** \brief merge node from two bit levels
  * \ingroup mofem
  */
struct NodeMergerInterface: public UnknownInterface {

  PetscErrorCode queryInterface(const MOFEMuuid& uuid, UnknownInterface** iface);

  MoFEM::Core& cOre;
  NodeMergerInterface(MoFEM::Core& core):
  cOre(core),
  successMerge(false),
  errorIfNoCommonEdge(false) {
  }

  /**
   * \brief Return true if successful merge.
   * @return Error code
   */
  inline bool getSucessMerge() { return successMerge; }

  /**
   * \brief Set error if no common edge
   * @param  b If true send error if false no error
   */
  inline void setErrorIfNoCommonEdge(const bool b = true) {
    errorIfNoCommonEdge = b;
  }


  /** \brief merge nodes which sharing edge

    Father is sties, mother is merged.

    \param father node to which mother is merged to.
    \param mother merged node
    \param test only tets_ptr from range are changed
    \param only_if_improve_quality Do merge if that improve quality
    \param move father by fraction of edge length move=[0,1]

    Move node on the edge, 0 not move, 1 move to mother side, 0.5 will be in the
    middle.

    */
  PetscErrorCode mergeNodes(
    EntityHandle father,
    EntityHandle mother,
    BitRefLevel bit,
    Range *tets_ptr = NULL,
    const bool only_if_improve_quality = false,
    const double move = 0
  );

  /** \brief merge nodes which sharing edge

    Father is sties, mother is merged.

    \param father node to which mother is merged to.
    \param mother merged node
    \param tets_from_bit_ref_level only tetrahedrons from bit level are changed
    \param only_if_improve_quality Do merge if that improve quality
    \param move father by fraction of edge length move=[0,1]

    Move node on the edge, 0 not move, 1 move to mother side, 0.5 will be in the
    middle.

    */
  PetscErrorCode mergeNodes(
    EntityHandle father,
    EntityHandle mother,
    BitRefLevel bit,
    BitRefLevel tets_from_bit_ref_level,
    const bool only_if_improve_quality = false,
    const double move = 0
  );

private:

  bool successMerge; ///< True if marge is success
  bool errorIfNoCommonEdge; ///< Send error if no common edge

};

}

#endif //__NODE_MERGER_HPP__
