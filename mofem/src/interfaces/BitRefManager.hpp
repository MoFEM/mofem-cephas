/** \file BitRefManager.hpp
 * \brief Interface managing BitRefLevels
 * \mofem_bit_ref
 *
 * Managing BitRef levels
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

#ifndef __BITREFMANAGER_HPP__
#define __BITREFMANAGER_HPP__

#include "UnknownInterface.hpp"

namespace MoFEM {

  static const MOFEMuuid IDD_MOFEMBitRefManager = MOFEMuuid( BitIntefaceId(BITREFMANAGER_INTERFACE) );

  /**
   * \brief Managing BitRefLevels
   * \mofem_bit_ref
   *
   */
  struct BitRefManager: public UnknownInterface {

    PetscErrorCode queryInterface(const MOFEMuuid& uuid, UnknownInterface** iface);

    MoFEM::Core& cOre;
    bool dEbug;

    BitRefManager(const MoFEM::Core& core);

    /**
    * \brief Destructor
    */
    ~BitRefManager();

    /**
     * \brief add entities to database and set bit ref level
     * \mofem_bit_ref
     *

     Example:\code
     EntityHandle meshset1; //contains ent1,ent2,ent3
     BitRefLevel myLevel0;
     myLevel0.set(0);
     m_field.query_interface<BitRefManager>()->setBitRefLevel(meshset1,myLevel0);
     //refine meshset1 into meshset2 and get new ents which are ent4, ent5
     EntityHandle meshset2; //contains ent1,ent2,ent3,ent4,ent5
     BitRefLevel myLevel1;
     myLevel1.set(1);
     m_field.query_interface<BitRefManager>()->setBitRefLevel(meshset2,myLevel1);
     \endcode

     * So entities 1,2,3 would be assigned to bit level 0 and 1 <br>
     * ent1[1,1,0,0,0,0,0], ent2[1,1,0,0,0,0,0], ent3[1,1,0,0,0,0,0], <br>
     * and entities 4 and 5 are assigned to bit level 1 only <br>
     * ent4[0,1,0,0,0,0,0], ent5[0,1,0,0,0,0,0] <br>


     * @param  ents      entities to set
     * @param  bit       bit refinement level
     * @param  only_tets only add entities on tetrahedral (obsolete need to be fixed)
     * @param  verb      verbosity level
     * @return           error code
     */
    PetscErrorCode setBitRefLevel(
      const Range &ents,const BitRefLevel &bit,const bool only_tets = true,int verb = 0
    ) const;

    PetscErrorCode setBitRefLevelByDim(
      const EntityHandle meshset,const int dim,const BitRefLevel &bit,int verb = 0
    ) const;

    PetscErrorCode setBitRefLevelByType(
      const EntityHandle meshset,const EntityType type,const BitRefLevel &bit,int verb = 0
    ) const;

    /** brief add meshset and set bit ref level
    * \mofem_bit_ref
     *
     * \param EntityHandle MeshSet
     * \param BitRefLevel bitLevel
     */
    PetscErrorCode setBitLevelToMeshset(
      const EntityHandle meshset,const BitRefLevel &bit,int verb
    ) const;

  };

}

#endif //__BITREFMANAGER_HPP__


/***************************************************************************//**
 * \defgroup mofem_bit_ref BitRefLevels manager
 * \brief Managing BitRefLevels
 *
 * \ingroup mofem
 ******************************************************************************/
