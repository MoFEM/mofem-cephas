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

/**
 * \brief Simple interface for fast problem set-up
 * \ingroup mofem_simple_interface
 */
struct Simple : public UnknownInterface {

  MoFEMErrorCode query_interface(boost::typeindex::type_index type_index,
                                 UnknownInterface **iface) const;

  Simple(const MoFEM::Core &core);
  virtual ~Simple() = default;

  /**
   * \brief get options
   * @return error code
   */
  MoFEMErrorCode getOptions();

  /**
   * \brief Load mesh file
   * @param  options file load options
   * @param  mesh_file_name file name if not set default or set by command line
   * is used.
   * @return            error code
   */
  MoFEMErrorCode loadFile(const std::string options,
                          const std::string mesh_file_name);
  /**
   * \brief Load mesh file with parallel options if number of cores > 1
   * @param  mesh_file_name file name if not set default or set by command line
   * is used.
   * @return            error code
   */
  MoFEMErrorCode loadFile(const std::string mesh_file_name = "");
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
   * @brief Remove field form domain
   *
   * @param name
   * @return MoFEMErrorCode
   */
  MoFEMErrorCode removeDomainField(const std::string &name);

  /**
   * @brief Remove field form boundary
   *
   * @param name
   * @return MoFEMErrorCode
   */
  MoFEMErrorCode removeBoundaryField(const std::string &name);

  /**
   * @brief Remove field form skeleton
   *
   * @param name
   * @return MoFEMErrorCode
   */
  MoFEMErrorCode removeSkeletonField(const std::string &name);

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
   * @brief Set the skeleton adjacency object
   *
   * @param dim
   * @return MoFEMErrorCode
   */
  MoFEMErrorCode setSkeletonAdjacency(int dim = -1);

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
   * @brief Rebuild internal MoFEM data structures
   *
   * Call this function after you add field or remove it.
   *
   * \note If you add field, or remove it, finite element and problem needs to
   * be rebuild. However DM can remain the same.
   *
   * @return MoFEMErrorCode
   */
  MoFEMErrorCode reSetUp();

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
   * @brief Get the problem dimension
   *
   * Problem dimension is determined by highest dimension of entities on the
   * mesh.
   *
   * @return int
   */
  inline int getDim() const { return dIm; }

  /**
   * @brief Set the problem dimension
   *
   * @return int
   */
  void setDim(int dim) { dIm = dim; };

  /**
   * @brief Get the MeshSet object
   *
   * @return EntityHandle&
   */
  inline EntityHandle &getMeshSet() { return meshSet; }

  /**
   * @brief Get the BoundaryMeshSet object
   *
   * @return EntityHandle&
   */
  inline EntityHandle &getBoundaryMeshSet() { return boundaryMeshset; }

  /**
   * @brief Get the SkeletonMeshSet object
   *
   * @return EntityHandle&
   */
  inline EntityHandle &getSkeletonMeshSet() { return skeletonMeshset; }

  /**
   * @brief Get the BitRefLevel
   *
   * @return BitRefLevel
   */
  inline BitRefLevel &getBitRefLevel() { return bitLevel; }

  /**
   * @brief Get the BitRefLevel
   *
   * @return BitRefLevel
   */
  inline BitRefLevel &getBitRefLevelMask() { return bitLevelMask; }

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

  /**
   * @brief Delete dm
   *
   * @return MoFEMErrorCode
   */
  MoFEMErrorCode deleteDM();

  /**
   * @brief Delete finite elements
   *
   * @return MoFEMErrorCode
   */
  MoFEMErrorCode deleteFiniteElements();

  /**
   * @brief Get the Add Skeleton F E object
   *
   * If variable set to true, skeleton element is created regardless field on
   * skelton is added or not.
   *
   * @return true
   * @return false
   */
  bool &getAddSkeletonFE() { return addSkeletonFE; }

  /**
   * @brief Get the Add Boundary F E object
   *
   * If variable set to true, boundary element is created regardless field on
   * skelton is added or not.
   *
   * @return true
   * @return false
   */
  bool &getAddBoundaryFE() { return addBoundaryFE; }

  /**
   * @brief add empty block to problem
   *
   * MatrixManager assumes that all blocks, i.e. all fields combinations are non
   * zero. This is not always the case, to optimise code and reduce memory usage
   * user can specifi which blocks are empty.
   *
   * @param row_field row filed name
   * @param col_field col field name
   * @return MoFEMErrorCode
   */
  MoFEMErrorCode addFieldToEmptyFieldBlocks(const std::string row_field,
                                            const std::string col_field) const;

private:
  MoFEM::Core &cOre;

  BitRefLevel bitLevel;     ///< BitRefLevel of the probelm
  BitRefLevel bitLevelMask; ///< BitRefLevel of the probelm

  PetscLogEvent MOFEM_EVENT_SimpleSetUP;
  PetscLogEvent MOFEM_EVENT_SimpleLoadMesh;
  PetscLogEvent MOFEM_EVENT_SimpleBuildFields;
  PetscLogEvent MOFEM_EVENT_SimpleBuildFiniteElements;
  PetscLogEvent MOFEM_EVENT_SimpleBuildProblem;
  PetscLogEvent MOFEM_EVENT_SimpleKSPSolve;

  EntityHandle meshSet;         ///< domain meshset
  EntityHandle boundaryMeshset; ///< meshset with boundary
  EntityHandle skeletonMeshset; ///< skeleton meshset with boundary

  bool addSkeletonFE; ///< Add skeleton FE
  bool addBoundaryFE; ///< Add boundary FE

  std::vector<std::string> domainFields;      ///< domain fields
  std::vector<std::string> boundaryFields;    ///< boundary fields
  std::vector<std::string> skeletonFields;    ///< fields on the skeleton
  std::vector<std::string> dataFields;        ///< Data fields
  std::vector<std::string> noFieldFields;     ///< NOFIELD field name
  std::vector<std::string> noFieldDataFields; ///< NOFIELD field name

  std::multimap<std::string, std::pair<int, Range>>
      fieldsOrder; ///< fields order

  std::string nameOfProblem; ///< problem name
  std::string domainFE;      ///< domain finite element
  std::string boundaryFE;    ///< boundary finite element
  std::string skeletonFE;    ///< skeleton finite element

  std::vector<std::string> otherFEs; ///< Other finite elements

  char meshFileName[255]; ///< mesh file name
  int dIm;                ///< dimension of problem

  SmartPetscObj<DM>
      dM; ///< Discrete manager (interface to PETSc/MoFEM functions)

  template <int DIM = -1> MoFEMErrorCode setSkeletonAdjacency();
};

} // namespace MoFEM

#endif // __SIMPLE_HPP__

/**
 * \defgroup mofem_simple_interface Simple interface
 * \brief Implementation of simple interface for fast problem set-up.
 *
 * \ingroup mofem
 **/
