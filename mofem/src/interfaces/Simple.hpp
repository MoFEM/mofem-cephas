/** \file Simple.hpp
 * \brief Header file for simple interface
 * \ingroup mofem_simple_interface
 *
 * Make simplified interface, to speedup problem setup and analysts.
 * See discussion here
 * <a
 * href=https://groups.google.com/d/msg/mofem-group/Vkc00aia4dU/o9RF3ZmPAAAJ>link
 * to google groups</a>
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

static const MOFEMuuid IDD_MOFEMSimple =
    MOFEMuuid(BitIntefaceId(SIMPLE_INTERFACE));

/**
 * \brief Simple interface for fast problem set-up
 * \ingroup mofem_simple_interface
 */
struct Simple : public UnknownInterface {

  MoFEMErrorCode query_interface(const MOFEMuuid &uuid,
                                 UnknownInterface **iface) const;

  Simple(const MoFEM::Core &core);

  /**
   * \brief get options
   * @return error code
   */
  MoFEMErrorCode getOptions();

  /**
   * \brief Load mesh file
   * @param  field_name file name
   * @return            error code
   */
  MoFEMErrorCode
  loadFile(const std::string options = "PARALLEL=READ_PART;"
                                       "PARALLEL_RESOLVE_SHARED_ENTS;"
                                       "PARTITION=PARALLEL_PARTITION;");

  /**
   * \brief Add field on domain
   * @param  name              name of the filed
   * @param  space             space (L2,H1,Hdiv,Hcurl)
   * @param  base              approximation base, see FieldApproximationBase
   * @param  nb_of_coefficients number of field coefficients
   * @param  tag_type          type of the tag MB_TAG_DENSE or MB_TAG_SPARSE
   * (DENSE is faster and uses less memory, SPARSE is more flexible if you
   * define field on subdomains)
   * @param  bh                if MF_EXCL throws error if field exits, MF_ZERO
   * no error if field exist
   * @param  verb              verbosity level
   * @return                   error code
   */
  MoFEMErrorCode
  addDomainField(const std::string &name, const FieldSpace space,
                 const FieldApproximationBase base,
                 const FieldCoefficientsNumber nb_of_coefficients,
                 const TagType tag_type = MB_TAG_SPARSE,
                 const enum MoFEMTypes bh = MF_ZERO, int verb = -1);

  /**
   * \brief Add field on boundary
   * @param  name              name of the filed
   * @param  space             space (L2,H1,Hdiv,Hcurl)
   * @param  base              approximation base, see FieldApproximationBase
   * @param  nb_of_coefficients number of field coefficients
   * @param  tag_type          type of the tag MB_TAG_DENSE or MB_TAG_SPARSE
   * (DENSE is faster and uses less memory, SPARSE is more flexible if you
   * define field on subdomains)
   * @param  bh                if MF_EXCL throws error if field exits, MF_ZERO
   * no error if field exist
   * @param  verb              verbosity level
   * @return                   error code
   */
  MoFEMErrorCode
  addBoundaryField(const std::string &name, const FieldSpace space,
                   const FieldApproximationBase base,
                   const FieldCoefficientsNumber nb_of_coefficients,
                   const TagType tag_type = MB_TAG_SPARSE,
                   const enum MoFEMTypes bh = MF_ZERO, int verb = -1);

  /**
   * \brief Add field on skeleton
   * @param  name              name of the filed
   * @param  space             space (L2,H1,Hdiv,Hcurl)
   * @param  base              approximation base, see FieldApproximationBase
   * @param  nb_of_coefficients number of field coefficients
   * @param  tag_type          type of the tag MB_TAG_DENSE or MB_TAG_SPARSE
   * (DENSE is faster and uses less memory, SPARSE is more flexible if you
   * define field on subdomains)
   * @param  bh                if MF_EXCL throws error if field exits, MF_ZERO
   * no error if field exist
   * @param  verb              verbosity level
   * @return                   error code
   */
  MoFEMErrorCode
  addSkeletonField(const std::string &name, const FieldSpace space,
                   const FieldApproximationBase base,
                   const FieldCoefficientsNumber nb_of_coefficients,
                   const TagType tag_type = MB_TAG_SPARSE,
                   const enum MoFEMTypes bh = MF_ZERO, int verb = -1);

  /**
   * \brief Add field on domain
   * @param  name              name of the filed
   * @param  space             space (L2,H1,Hdiv,Hcurl)
   * @param  base              approximation base, see FieldApproximationBase
   * @param  nb_of_coefficients number of field coefficients
   * @param  tag_type          type of the tag MB_TAG_DENSE or MB_TAG_SPARSE
   * (DENSE is faster and uses less memory, SPARSE is more flexible if you
   * define field on subdomains)
   * @param  bh                if MF_EXCL throws error if field exits, MF_ZERO
   * no error if field exist
   * @param  verb              verbosity level
   * @return                   error code
   */
  MoFEMErrorCode addDataField(const std::string &name, const FieldSpace space,
                              const FieldApproximationBase base,
                              const FieldCoefficientsNumber nb_of_coefficients,
                              const TagType tag_type = MB_TAG_SPARSE,
                              const enum MoFEMTypes bh = MF_ZERO,
                              int verb = -1);

  /**
   * \brief Define finite elements
   * @return         Error code
   */
  MoFEMErrorCode defineFiniteElements();

  /**
   * \brief define problem
   * @return error code
   */
  MoFEMErrorCode defineProblem(const PetscBool is_partitioned = PETSC_TRUE);

  /**
   * \brief Set field order
   * @param  std::field_name field name
   * @param  order           order
   * @param  range of entities to which order is set (If null it sat to all
   * entities)
   * @return                 error code
   */
  MoFEMErrorCode setFieldOrder(const std::string field_name, const int order,
                               const Range *ents = NULL);

  /**
   * \brief Build fields
   * @return error code
   */
  MoFEMErrorCode buildFields();

  /**
   * \brief Build finite elements
   * @return error code
   */
  MoFEMErrorCode buildFiniteElements();

  /**
   * \brief Build problem
   * @return error code
   */
  MoFEMErrorCode buildProblem();

  /**
   * \brief Setup problem
   * @return error code
   */
  MoFEMErrorCode setUp(const PetscBool is_partitioned = PETSC_TRUE);

  /**
   * \brief Get DM
   * @param  dm discrete manager
   * @return    error code
   */
  MoFEMErrorCode getDM(DM *dm);

  /**
   * @brief Return smart DM object
   *
   * \code
   * {
   *    auto dm  = simple_interface->getDM();
   *
   *    // ...
   *
   *    // dm is automatically destroyed when out of the scope
   * }
   * \endcode
   *
   * @return SmartPetscObj<DM>
   */
  inline SmartPetscObj<DM> getDM() { return dM; }

  /**
   * @brief Get the problem dimensio
   *
   * Problem dimension is determined by highest dimension of entities on the
   * mesh.
   *
   * @return int
   */
  inline int getDim() const { return dIm; }

  /**
   * @brief Get the Domain FE Name
   *
   * @return const std::string
   */
  inline const std::string getDomainFEName() const { return domainFE; }

  /**
   * @brief Get the Boundary FE Name
   *
   * @return const std::string
   */
  inline const std::string getBoundaryFEName() const { return boundaryFE; }

  /**
   * @brief Get the Skeleton FE Name
   *
   * @return const std::string
   */
  inline const std::string getSkeletonFEName() const { return skeletonFE; }

  /**
   * @brief Get the Problem Name
   *
   * @return const std::string
   */
  inline const std::string getProblemName() const { return nameOfProblem; }

  /**
   * @brief Get the Domain FE Name
   *
   * @return std::string&
   */
  inline std::string &getDomainFEName() { return domainFE; }

  /**
   * @brief Get the Boundary FE Name
   *
   * @return std::string&
   */
  inline std::string &getBoundaryFEName() { return boundaryFE; }

  /**
   * @brief Get the Skeleton FE Name
   *
   * @return std::string&
   */
  inline std::string &getSkeletonFEName() { return skeletonFE; }

  /**
   * @brief Get the Problem Name
   *
   * @return std::string&
   */
  inline std::string &getProblemName() { return nameOfProblem; }

  /**
   * @brief Get the Other Finite Elements
   *
   * User can create finite elements using directly core interface and
   * and them to simple problem by this function
   *
   * @return std::vector<std::string>&
   */
  inline std::vector<std::string> &getOtherFiniteElements() { return otherFEs; }

private:
  MoFEM::Core &cOre;

  const BitRefLevel bitLevel; ///< BitRefLevel of the probelm

  PetscLogEvent MOFEM_EVENT_SimpleLoadMesh;
  PetscLogEvent MOFEM_EVENT_SimpleBuildFields;
  PetscLogEvent MOFEM_EVENT_SimpleBuildFiniteElements;
  PetscLogEvent MOFEM_EVENT_SimpleBuildProblem;
  PetscLogEvent MOFEM_EVENT_SimpleKSPSolve;

  EntityHandle meshSet;                       ///< domain meshset
  EntityHandle boundaryMeshset;               ///< meshset with boundary
  std::vector<std::string> domainFields;      ///< domain fields
  std::vector<std::string> boundaryFields;    ///< boundary fields
  std::vector<std::string> skeletonFields;    ///< fields on the skeleton
  std::vector<std::string> dataFields;        ///< Data fields
  std::vector<std::string> noFieldFields;     ///< NOFIELD field name
  std::vector<std::string> noFieldDataFields; ///< NOFIELD field name

  std::map<std::string, std::pair<int, Range>> fieldsOrder; ///< fields order

  std::string nameOfProblem; ///< problem name
  std::string domainFE;      ///< domain finite element
  std::string boundaryFE;    ///< boundary finite element
  std::string skeletonFE;    ///< skeleton finite element

  std::vector<std::string> otherFEs; ///< Other finite elements

  char meshFileName[255]; ///< mesh file name
  int dIm;                ///< dimension of problem

  SmartPetscObj<DM>
      dM; ///< Discrete manager (interface to PETSc/MoFEM functions)
};

} // namespace MoFEM

#endif // __SIMPLE_HPP__

/**
 * \defgroup mofem_simple_interface Simple interface
 * \brief Implementation of simple interface for fast problem set-up.
 *
 * \ingroup mofem
 **/
