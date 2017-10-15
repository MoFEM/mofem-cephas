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

    #ifdef FTENSOR_HPP

    template<class T>
    static double volumeLengthQuality(T *coords) {
      T lrms = 0;
      for(int dd = 0;dd!=3;dd++) {
        lrms +=
        pow(coords[0*3+dd]-coords[1*3+dd],2)+
        pow(coords[0*3+dd]-coords[2*3+dd],2)+
        pow(coords[0*3+dd]-coords[3*3+dd],2)+
        pow(coords[1*3+dd]-coords[2*3+dd],2)+
        pow(coords[1*3+dd]-coords[3*3+dd],2)+
        pow(coords[2*3+dd]-coords[3*3+dd],2);
      }
      lrms = sqrt((1./6.)*lrms);
      T diff_n[12];
      ShapeDiffMBTET(diff_n);
      FTensor::Tensor1<T*,3> t_diff_n(&diff_n[0],&diff_n[1],&diff_n[2],3);
      FTensor::Tensor1<T*,3> t_coords(&coords[0],&coords[1],&coords[2],3);
      FTensor::Tensor2<T,3,3> jac;
      FTensor::Index<'i',3> i;
      FTensor::Index<'j',3> j;
      jac(i,j) = 0;
      for(int nn = 0;nn!=4;nn++) {
        jac(i,j) += t_coords(i)*t_diff_n(j);
        ++t_coords;
        ++t_diff_n;
      }
      T volume = dEterminant(jac)/6.;
      return 6.*sqrt(2.)*volume/pow(lrms,3);
    }

    #endif // FTENSOR_HPP

    /**@}*/

    /** \name Entity hanlders by bit ref level */

    /**@{*/

    /**\brief add all ents from ref level given by bit to meshset
      * \ingroup mofem_tools
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
     * \param BitRefLevel bitLevel
     * \param BitRefLevel mask
     * \retval ents
     */
    PetscErrorCode getEntitiesByRefLevel(
      const BitRefLevel &bit,const BitRefLevel &mask,Range &ents
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
