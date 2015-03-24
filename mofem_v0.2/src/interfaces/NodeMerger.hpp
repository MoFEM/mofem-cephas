/** \file NodeMerger.hpp
 * \brief NodeMerger interface 
 * 
 * Low level data structures not used directly by user
 *
 * MoFEM is distributed in the hope that it will be useful, but WITHOUT
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


static const MOFEMuuid IDD_MOFENNodeMerger = MOFEMuuid( BitIntefaceId(NODEMERGER_INTERFACE) );

/** \brief merge node from two bit levels
  * \ingroup mofem
  */
struct NodeMergerInterface: public FieldUnknownInterface {

  PetscErrorCode queryInterface(const MOFEMuuid& uuid, FieldUnknownInterface** iface);

  MoFEM::Core& cOre;
  NodeMergerInterface(MoFEM::Core& core): cOre(core) {};

  /** \brief merge nodes which sharing edge

    I apologise that it could be traditional her. Father is sties, mother is merged.

    \param father node to which mother is merged to.
    \param mother merged node
    \param test only tets_ptr from range are changed

    */
  PetscErrorCode mergeNodes(EntityHandle father,EntityHandle mother,BitRefLevel bit,Range *tets_ptr = NULL);

  /** \brief merge nodes which sharing edge

    I apologise that it could be traditional her. Father is sties, mother is merged.

    \param father node to which mother is merged to.
    \param mother merged node
    \param tets_from_bit_ref_level only tets from bit level are changed

    */
  PetscErrorCode mergeNodes(EntityHandle father,EntityHandle mother,BitRefLevel bit,BitRefLevel tets_from_bit_ref_level);

};

}

#endif //__NODE_MERGER_HPP__
