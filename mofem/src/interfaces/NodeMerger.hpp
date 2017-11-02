/** \file NodeMerger.hpp
 * \brief NodeMerger interface
 *
 * Node merger interface
 *
 * \ingroup mofem_node_merger
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
  *
  * \ingroup mofem_node_merger
  */
struct NodeMergerInterface: public UnknownInterface {

  PetscErrorCode query_interface(const MOFEMuuid& uuid, UnknownInterface** iface) const;

  MoFEM::Core& cOre;
  NodeMergerInterface(const MoFEM::Core& core):
  cOre(const_cast<MoFEM::Core&>(core)),
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

  /**
   * \brief calculated quality of tets adjacent to edge
   * @param  edge        edge handle
   * @param  tets_ptr    pointer to range of tets
   * @param  min_quality returned quality of tets
   * @param  th          handle to tag with nodes positions
   * @return             error code
   */
  PetscErrorCode edgeMinQuality(
    EntityHandle edge,const Range *tets_ptr,double &min_quality,Tag th = NULL
  );

  /** \brief merge nodes which sharing edge

    Father is sties, mother is merged.

    \param father node to which mother is merged to.
    \param mother merged node
    \param tetrahedra after merge
    \param test only tets_ptr from range are changed
    \param only_if_improve_quality Do merge if that improve quality
    \param move father by fraction of edge length move=[0,1]

    Move node on the edge, 0 not move, 1 move to mother side, 0.5 will be in the
    middle.

    */
  PetscErrorCode mergeNodes(
    EntityHandle father,
    EntityHandle mother,
    Range &out_tets,
    Range *tets_ptr = NULL,
    const bool only_if_improve_quality = false,
    const double move = 0,
    const int line_search = 0,
    Tag th = NULL,
    const int verb = 0
  );


  /** \brief merge nodes which sharing edge

    Father is sties, mother is merged.

    \param father node to which mother is merged to.
    \param mother merged node
    \param bit level of mesh merged nodes mesh
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
    const double move = 0,
    Tag th = NULL
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
    const double move = 0,
    Tag th = NULL
  );

  struct ParentChild {
    EntityHandle pArent;
    EntityHandle cHild;
    ParentChild(const EntityHandle parent,const EntityHandle child):
    pArent(parent),
    cHild(child) {
    }
  };

  typedef multi_index_container<
    ParentChild,
    indexed_by<
      hashed_unique<
        member<ParentChild,EntityHandle,&ParentChild::pArent>
      >,
      hashed_non_unique<
        member<ParentChild,EntityHandle,&ParentChild::cHild>
      >
    >
  > ParentChildMap;

  /**
   * \brief Get map of parent cand child
   * @return
   */
  inline ParentChildMap& getParentChildMap() {
    return parentChildMap;
  }

private:

  bool successMerge; ///< True if marge is success
  bool errorIfNoCommonEdge; ///< Send error if no common edge

  /**
   * \brief Calualte quality if nodes merged
   * @param  check_tests tets to check
   * @param  father      fisrt node of the edge
   * @param  mother      second node of the edge
   * @param  coords_move moved father node
   * @param  min_quality calculated quality
   * @return             error code
   */
  PetscErrorCode minQuality(
    Range &check_tests,
    EntityHandle father,
    EntityHandle mother,
    double *coords_move,
    double &min_quality,
    Tag th = NULL
  );

  /**
   * \brief Use bisection method to find point of edge collapse
   * @param  check_tests range of tets to check quality
   * @param  father      first node of the edge
   * @param  mother      second node of the edge
   * @param  line_search number of iterations
   * @param  coords_move node to move
   * @return             error code
   */
  PetscErrorCode lineSearch(
    Range &check_tests,
    EntityHandle father,
    EntityHandle mother,
    int line_search,
    double *coords_move,
    Tag th = NULL
  );

  ParentChildMap parentChildMap;

  struct FaceMap {
    EntityHandle e,n0,n1;
    FaceMap(const EntityHandle e,const EntityHandle n0,const EntityHandle n1):
    e(e),n0(n0),n1(n1) {
    }
  };

  typedef multi_index_container<
    FaceMap,
    indexed_by<
      hashed_unique<
        composite_key<
          FaceMap,
          member<FaceMap,EntityHandle,&FaceMap::n0>,
          member<FaceMap,EntityHandle,&FaceMap::n1>
        >
      >
    >
  > FaceMapIdx;

};

}

#endif //__NODE_MERGER_HPP__

/***************************************************************************//**
 * \defgroup mofem_node_merger Node merger
 * \brief Interface for node merger
 *
 * \ingroup mofem
 ******************************************************************************/
