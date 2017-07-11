/** \file Simple.hpp
 * \brief Header file for simple interface
 * \ingroup mofem_simple_interface
 *
 * Make simplified interface, to speedup problem setup and analysts.
 * See discussion here
 * <a href=https://groups.google.com/d/msg/mofem-group/Vkc00aia4dU/o9RF3ZmPAAAJ>link to google groups</a>
 *
 */

/* MoFEM is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
 * License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with MoFEM. If not, see <http://www.gnu.org/licenses/>
*/

#ifndef __SIMPLE_HPP__
#define __SIMPLE_HPP__

#include "UnknownInterface.hpp"

namespace MoFEM {

  static const MOFEMuuid IDD_MOFEMSimple = MOFEMuuid( BitIntefaceId(SIMPLE_INTERFACE) );

  /**
   * \brief Simple interface for fast problem set-up
   * \ingroup mofem_simple_interface
   */
  struct Simple: public UnknownInterface {

    PetscErrorCode queryInterface(const MOFEMuuid& uuid, UnknownInterface** iface);

    MoFEM::Core& cOre;
    Simple(const MoFEM::Core& core);

    /**
    * \brief Destructor
    */
    ~Simple();

    /**
     * \brief get options
     * @return error code
     */
    PetscErrorCode getOptions();

    /**
     * \brief Load mesh file
     * @param  field_name file name
     * @return            error code
     */
    PetscErrorCode loadFile();

    /**
     * \brief Add field on domain
     * @param  name              name of the filed
     * @param  space             space (L2,H1,Hdiv,Hcurl)
     * @param  base              approximation base, see FieldApproximationBase
     * @param  nb_of_cooficients number of field coefficients
     * @param  tag_type          type of the tag MB_TAG_DENSE or MB_TAG_SPARSE (DENSE is faster and uses less memory, SPARSE is more flexible if you define field on subdomains)
     * @param  bh                if MF_EXCL throws error if field exits, MF_ZERO no error if field exist
     * @param  verb              verbosity leve
     * @return                   error code
     */
    PetscErrorCode addDomainField(
      const std::string& name,
      const FieldSpace space,
      const FieldApproximationBase base,
      const FieldCoefficientsNumber nb_of_cooficients,
      const TagType tag_type = MB_TAG_SPARSE,
      const enum MoFEMTypes bh = MF_ZERO,
      int verb = -1
    );

    /**
     * \brief Add field on boundary
     * @param  name              name of the filed
     * @param  space             space (L2,H1,Hdiv,Hcurl)
     * @param  base              approximation base, see FieldApproximationBase
     * @param  nb_of_cooficients number of field coefficients
     * @param  tag_type          type of the tag MB_TAG_DENSE or MB_TAG_SPARSE (DENSE is faster and uses less memory, SPARSE is more flexible if you define field on subdomains)
     * @param  bh                if MF_EXCL throws error if field exits, MF_ZERO no error if field exist
     * @param  verb              verbosity leve
     * @return                   error code
     */
    PetscErrorCode addBoundaryField(
      const std::string& name,
      const FieldSpace space,
      const FieldApproximationBase base,
      const FieldCoefficientsNumber nb_of_cooficients,
      const TagType tag_type = MB_TAG_SPARSE,
      const enum MoFEMTypes bh = MF_ZERO,
      int verb = -1
    );

    /**
     * \brief Add field on skeleton
     * @param  name              name of the filed
     * @param  space             space (L2,H1,Hdiv,Hcurl)
     * @param  base              approximation base, see FieldApproximationBase
     * @param  nb_of_cooficients number of field coefficients
     * @param  tag_type          type of the tag MB_TAG_DENSE or MB_TAG_SPARSE (DENSE is faster and uses less memory, SPARSE is more flexible if you define field on subdomains)
     * @param  bh                if MF_EXCL throws error if field exits, MF_ZERO no error if field exist
     * @param  verb              verbosity leve
     * @return                   error code
     */
    PetscErrorCode addSkeletionField(
      const std::string& name,
      const FieldSpace space,
      const FieldApproximationBase base,
      const FieldCoefficientsNumber nb_of_cooficients,
      const TagType tag_type = MB_TAG_SPARSE,
      const enum MoFEMTypes bh = MF_ZERO,
      int verb = -1
    );

    /**
     * \brief Add field on domain
     * @param  name              name of the filed
     * @param  space             space (L2,H1,Hdiv,Hcurl)
     * @param  base              approximation base, see FieldApproximationBase
     * @param  nb_of_cooficients number of field coefficients
     * @param  tag_type          type of the tag MB_TAG_DENSE or MB_TAG_SPARSE (DENSE is faster and uses less memory, SPARSE is more flexible if you define field on subdomains)
     * @param  bh                if MF_EXCL throws error if field exits, MF_ZERO no error if field exist
     * @param  verb              verbosity leve
     * @return                   error code
     */
    PetscErrorCode addDataField(
      const std::string& name,
      const FieldSpace space,
      const FieldApproximationBase base,
      const FieldCoefficientsNumber nb_of_cooficients,
      const TagType tag_type = MB_TAG_SPARSE,
      const enum MoFEMTypes bh = MF_ZERO,
      int verb = -1
    );


    /**
     * \brief Define finite elements
     * @return         Error code
     */
    PetscErrorCode defineFiniteElements();

    /**
     * \brief define problem
     * @return error code
     */
    PetscErrorCode defineProblem();

    /**
     * \brief Set field order
     * @param  std::field_name field name
     * @param  order           order
     * @param  range of entities to which order is set (If null it sat to all entities)
     * @return                 error code
     */
    PetscErrorCode setFieldOrder(const std::string field_name,const int order,const Range* ents = NULL);

    /**
     * \brief Build fields
     * @return error code
     */
    PetscErrorCode buildFields();

    /**
     * \brief Build finite elements
     * @return error code
     */
    PetscErrorCode buildFiniteElements();

    /**
     * \brief Build problem
     * @return error code
     */
    PetscErrorCode buildProblem();

    /**
     * \brief Setup problem
     * @return error code
     */
    PetscErrorCode setUp();

    /**
     * \brief Get DM
     * @param  dm discrete manager
     * @return    error code
     */
    PetscErrorCode getDM(DM *dm);

    inline int getDim() const { return dIm; }
    inline const std::string& getDomainFEName() { return domainFE; }
    inline const std::string& getBoundaryFEName() { return boundaryFE; }
    inline const std::string& getSkeletonFEName() { return skeletonFE; }

  private:

    const BitRefLevel bitLevel;

    PetscLogEvent USER_EVENT_SimpleLoadMesh;
    PetscLogEvent USER_EVENT_SimpleBuildFields;
    PetscLogEvent USER_EVENT_SimpleBuildFiniteElements;
    PetscLogEvent USER_EVENT_SimpleBuildProblem;
    PetscLogEvent USER_EVENT_SimpleKSPSolve;

    EntityHandle meshSet;                     ///< domain meshset
    EntityHandle boundaryMeshset;             ///< meshset with boundary
    std::string nameOfProblem;                ///< problem name
    std::vector<std::string> domainFields;    ///< domain fields
    std::vector<std::string> boundaryFields;  ///< boundary fields
    std::vector<std::string> skeletonFields;  ///< fields on the skeleton
    std::vector<std::string> dataFields;      ///< Data fields

    std::map<std::string,std::pair<int,Range> > fieldsOrder;    ///< fields order

    std::string domainFE;      ///< domain finite element
    std::string boundaryFE;     ///< boundary finite element
    std::string skeletonFE;     ///< skeleton finite element

    char meshFileName[255]; ///< mesh file name
    int dIm;                ///< dimension of problem

    DM dM;

  };

}

#endif // __SIMPLE_HPP__

/***************************************************************************//**
 * \defgroup mofem_simple_interface Simple interface
 * \brief Implementation of simple interface for fast problem set-up.
 *
 * \ingroup mofem
 ******************************************************************************/
