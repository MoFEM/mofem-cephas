/** \file MoABTools.hpp
 * \brief MoABTools interface
 *
 * Implementatiom of some useful and very often used methods * in MoFEM.
 *
 * \ingroup mofem_moab_tools
 */

/* MoFEM is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
 * License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with MoFEM. If not, see <http://www.gnu.org/licenses/>
*/

#ifndef __TOOLS_HPP__
#define __TOOLS_HPP__

namespace MoFEM {

  static const MOFEMuuid IDD_MOFEMTools = MOFEMuuid( BitIntefaceId(TOOLS) );

  /**
   * \brief Auxilairy tools
   * \nosubgrouping
   * \ingroup mofem_tools
   */
  struct Tools: public UnknownInterface {

    PetscErrorCode queryInterface(const MOFEMuuid& uuid, UnknownInterface** iface);

    MoFEM::Core& cOre;
    Tools(const MoFEM::Core& core):
    cOre(const_cast<MoFEM::Core&>(core)) {
    }

    /** \name Computational */

    /**@{*/

    template<class T>
    static inline double dEterminant(T &t) {
      return
      +t(0,0)*t(1,1)*t(2,2) + t(1,0)*t(2,1)*t(0,2)
      +t(2,0)*t(0,1)*t(1,2) - t(0,0)*t(2,1)*t(1,2)
      -t(2,0)*t(1,1)*t(0,2) - t(1,0)*t(0,1)*t(2,2);
    }

    /**
     * \brief calulate tetrahedron volume length quality
     * @param  coords tet coordinates
     * @return        qilaity
     */
    static double volumeLengthQuality(const double *coords);

    /**
     * \brief calculate minimal quality of tetrahedra in range
     * @param  tets        range
     * @param  min_quality mimimal quality
     * @return             error code
     */
    PetscErrorCode minTetsQuality(const Range& tets,double &min_quality);

    /**@}*/

    /** \name Entity hanlders by bit ref level */

    /**@{*/

    /**\brief add all ents from ref level given by bit to meshset
      * \ingroup mofem_tools
      *
      * \note Entities NOT have to be added to MoFEM database
      *
      * \param BitRefLevel bitLevel
      * \param BitRefLevel mask
      * \param EntityType type of entities
      * \retval EntityHandle meshset
      *
      */
    PetscErrorCode getEntitiesByTypeAndRefLevel(
      const BitRefLevel &bit,const BitRefLevel &mask,const EntityType type,const EntityHandle meshset,int verb = 0
    ) const;

    /**\brief add all ents from ref level given by bit to meshset
     * \ingroup mofem_tools
     *
     * \note Entities NOT have to be added to MoFEM database
     *
     * \param BitRefLevel bitLevel
     * \param BitRefLevel mask
     * \param EntityType type of entities
     * \retval ents
     *
     */
    PetscErrorCode getEntitiesByTypeAndRefLevel(
      const BitRefLevel &bit,const BitRefLevel &mask,const EntityType type,Range &ents,int verb = 0
    ) const;

    /**\brief add all ents from ref level given by bit to meshset
     * \ingroup mofem_tools
     *
     * \note Entities NOT have to be added to MoFEM database
     *
     * \param BitRefLevel bitLevel
     * \param BitRefLevel mask
     * \param EntityHandle meshset
     *
     */
    PetscErrorCode getEntitiesByRefLevel(
      const BitRefLevel &bit,const BitRefLevel &mask,const EntityHandle meshset
    ) const;

    /**\brief add all ents from ref level given by bit to meshset
     * \ingroup mofem_tools
     *
     * \note Entities NOT have to be added to MoFEM database
     *
     * \param BitRefLevel bitLevel
     * \param BitRefLevel mask
     * \retval ents
     */
    PetscErrorCode getEntitiesByRefLevel(
      const BitRefLevel &bit,const BitRefLevel &mask,Range &ents
    ) const;


    /**
     * \brief get entities by bit ref level and type of parent
     *
     * \note Entities have to be added to MoFEM database
     *
     * \param BitRefLevel bitLevel
     * \param BitRefLevel mask
     * @param  type of parent
     * @param  ents returned ents
     * @return      error code
     */
    PetscErrorCode getEntitiesByParentType(
      const BitRefLevel &bit,const BitRefLevel &mask,const EntityType type,Range &ents
    ) const;

    /**@}*/

    /** \name Get adjacencies bi bit ref level */

    /**@{*/

    /** \brief Get the adjacencies associated with a entity to entities of a specified dimension.
      * \ingroup mofem_tools
      *
      * bit ref level of adjacent entities is equal to bit ref level of adjacent entities
      */
    virtual PetscErrorCode getAdjacenciesEquality(const EntityHandle from_entiti,const int to_dimension,Range &adj_entities) const;

    /** \brief Get the adjacencies associated with a entity to entities of a specified dimension.
      * \ingroup mofem_tools
      *
      * bit ref level of adjacent entities is any of bit ref level of adjacent entities
      */
    virtual PetscErrorCode getAdjacenciesAny(const EntityHandle from_entiti,const int to_dimension,Range &adj_entities) const;

    /** \brief Get the adjacencies associated with a entity to entities of a specified dimension.
      * \ingroup mofem_tools
      *
      * bit ref level of adjacent entities is equal to bit ref level of adjacent entities
      */
    virtual PetscErrorCode getAdjacencies(
      const Problem *problem_ptr,
      const EntityHandle *from_entities,
      const int num_netities,
      const int to_dimension,
      Range &adj_entities,
      const int operation_type = moab::Interface::INTERSECT,
      const int verb = 0
    ) const;

    /** \brief Get the adjacencies associated with a entity to entities of a specified dimension.
      * \ingroup mofem_tools
      *
      * bit ref level of adjacent entities is equal to bit ref level of adjacent entities
      */
    virtual PetscErrorCode getAdjacencies(
      const BitRefLevel &bit,
      const EntityHandle *from_entities,
      const int num_netities,
      const int to_dimension,
      Range &adj_entities,
      const int operation_type = moab::Interface::INTERSECT,
      const int verb = 0
    ) const;

    /**@}*/

    /** \name Writting files */

    /**@{*/

    PetscErrorCode writeBitLevelByType(
      const BitRefLevel& bit,
      const BitRefLevel& mask,
      const EntityType type,
      const char * 	file_name,
      const char * 	file_type,
      const char * 	options
    ) const;

    /**@}*/

  };

}

#endif // __TOOLS_HPP__a

/***************************************************************************//**
 * \defgroup mofem_tools Tools interface
 * \brief Interface for tools
 *
 * \ingroup mofem
 ******************************************************************************/
