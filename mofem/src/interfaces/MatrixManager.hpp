/** \file MatrixManager.hpp
 * \brief Interface for creating matrices and managing matrices
 * \ingroup mofem_mat_interface
 *
 * Creating and managing matrices
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

#ifndef __MATMANAGER_HPP__
#define __MATMANAGER_HPP__

#include "UnknownInterface.hpp"

namespace MoFEM {

static const MOFEMuuid IDD_MOFEMMatrixManager =
    MOFEMuuid(BitIntefaceId(MATRIX_MANAGER_INTERFACE));

/**
 * \brief Matrix manager is used to build and partition problems
 * \ingroup mofem_mat_interface
 *
 */
struct MatrixManager : public UnknownInterface {

  MoFEMErrorCode query_interface(const MOFEMuuid &uuid,
                                 UnknownInterface **iface) const;

  MoFEM::Core &cOre;
  MatrixManager(const MoFEM::Core &core);

  /**
   * \brief Destructor
   */
  virtual ~MatrixManager() = default;

  /**
   * @brief Creates a MPI AIJ matrix using arrays that contain in standard CSR
   * format the local rows.
   *
   * <a
   href=https://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/Mat/MatCreateMPIAIJWithArrays.html>See
   PETSc for details</a>

   * @tparam Tag
   * @param name
   * @param Aij
   * @param verb
   * @return MoFEMErrorCode
   */
  template <class Tag>
  MoFEMErrorCode createMPIAIJWithArrays(const std::string name, Mat *Aij,
                                        int verb = QUIET) {
    static_assert(!std::is_same<Tag, Tag>::value, "not implemented");
    return 0;
  }

  /** \copydoc MoFEM::MatrixManager::createMPIAIJWithArrays
   */
  template <class Tag>
  MoFEMErrorCode createMPIAIJWithArrays(const std::string name,
                                        SmartPetscObj<Mat> &aij_ptr,
                                        int verb = QUIET) {
    MoFEMFunctionBegin;
    Mat aij;
    CHKERR createMPIAIJWithArrays<Tag>(name, &aij, verb);
    aij_ptr.reset(aij, false);
    MoFEMFunctionReturn(0);
  }

  /**
   * @brief Creates a MPI AIJ matrix using arrays that contain in standard CSR
   * format the local rows.
   *
   * <a
   href=https://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/Mat/MATMPIAIJ.html>See
   PETSc for details</a>

   * @tparam Tag
   * @param name
   * @param Aij
   * @param verb
   * @return MoFEMErrorCode
   */
  template <class Tag>
  MoFEMErrorCode createMPIAIJ(const std::string name, Mat *Aij,
                                        int verb = QUIET) {
    static_assert(!std::is_same<Tag, Tag>::value, "not implemented");
    return 0;
  }

  /** \copydoc MoFEM::MatrixManager::createMPIAIJ
   */
  template <class Tag>
  MoFEMErrorCode createMPIAIJ(const std::string name,
                                        SmartPetscObj<Mat> &aij_ptr,
                                        int verb = QUIET) {
    MoFEMFunctionBegin;
    Mat aij;
    CHKERR createMPIAIJ<Tag>(name, &aij, verb);
    aij_ptr.reset(aij, false);
    MoFEMFunctionReturn(0);
  }

  /**
   * @brief Creates a sparse matrix representing an adjacency list.
   *
   * The matrix
   * does not have numerical values associated with it, but is intended for
   * ordering (to reduce bandwidth etc) and partitioning.
   *
   * <a
   * href=https://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/Mat/MatCreateMPIAdj.html>See
   * PETSc for details</a>
   *
   * \note This matrix object does not support most matrix operations, include
   * MatSetValues(). You must NOT free the ii, values and jj arrays yourself.
   * PETSc will free them when the matrix is destroyed; you must allocate them
   * with PetscMalloc(). If you call from Fortran you need not create the arrays
   * with PetscMalloc(). Should not include the matrix diagonals.
   *
   * @tparam Tag
   * @param name
   * @param Adj
   * @param verb
   * @return MoFEMErrorCode
   */
  template <class Tag>
  MoFEMErrorCode createMPIAdjWithArrays(const std::string name, Mat *Adj,
                                        int verb = QUIET) {
    static_assert(!std::is_same<Tag, Tag>::value, "not implemented");
    return MOFEM_NOT_IMPLEMENTED;
  }

  /**
   * @brief Create sequencial matrix
   *
   * Creates a sparse matrix in AIJ (compressed row) format (the default
   * parallel PETSc format). For good matrix assembly performance the user
   * should preallocate the matrix storage by setting the parameter nz (or the
   * array nnz). By setting these parameters accurately, performance during
   * matrix assembly can be increased by more than a factor of 50.
   *
   * <a
   * href=https://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/Mat/MatCreateSeqAIJ.html>See
   * PETSc for details</a>
   *
   * @tparam Tag
   * @param name
   * @param Aij
   * @param i
   * @param j
   * @param v
   * @param verb
   * @return MoFEMErrorCode
   */
  template <class Tag>
  MoFEMErrorCode createSeqAIJWithArrays(const std::string name, Mat *Aij,
                                        int verb = QUIET) {
    static_assert(!std::is_same<Tag, Tag>::value, "not implemented");
    return 0;
  }

    /**
   * @brief Create sequencial matrix
   *
   * Creates a sparse matrix in AIJ (compressed row) format (the default
   * parallel PETSc format). For good matrix assembly performance the user
   * should preallocate the matrix storage by setting the parameter nz (or the
   * array nnz). By setting these parameters accurately, performance during
   * matrix assembly can be increased by more than a factor of 50.
   *
   * <a
   * href=https://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/Mat/MatCreateSeqAIJ.html>See
   * PETSc for details</a>
   *
   * @tparam Tag
   * @param name
   * @param Aij (SmartPetscObj)
   * @param i
   * @param j
   * @param v
   * @param verb
   * @return MoFEMErrorCode
   */
  template <class Tag>
  MoFEMErrorCode createSeqAIJWithArrays(const std::string name,
                                        SmartPetscObj<Mat> &aij_ptr,
                                        int verb = QUIET) {
    MoFEMFunctionBegin;
    Mat aij;
    CHKERR createSeqAIJWithArrays<Tag>(name, &aij, verb);
    aij_ptr.reset(aij, false);
    MoFEMFunctionReturn(0);
  }

  /** \brief check if matrix fill in correspond to finite element indices

  This is used to check consistency of code. If problem is notices with
  additional non-zero elements in matrix, this function can help detect
  problem. Should be used as a part of atom tests

  * @param problem_name
  * @param row print info at particular row
  * @param col print info at particular col
  * @return MoFEMErrorCode

  */
  MoFEMErrorCode checkMatrixFillIn(const std::string problem_name,
                                   int row_print, int col_print, Mat A,
                                   int verb = QUIET);

  /** \brief check if matrix fill in correspond to finite element indices

  This is used to check consistency of code. If problem is notices with
  additional non-zero elements in matrix, this function can help detect
  problem. Should be used as a part of atom tests

  * @tparam Tag
  * @param problem_name
  * @param row print info at particular row
  * @param col print info at particular col
  * @return MoFEMErrorCode

  */
  template <class Tag>
  MoFEMErrorCode
  checkMPIAIJWithArraysMatrixFillIn(const std::string problem_name,
                                    int row_print, int col_print,
                                    int verb = QUIET) {
    static_assert(!std::is_same<Tag, Tag>::value, "not implemented");
    return 0;
  }

  /** \brief check if matrix fill in correspond to finite element indices

  This is used to check consistency of code. If problem is notices with
  additional non-zero elements in matrix, this function can help detect
  problem. Should be used as a part of atom tests

  * @tparam Tag
  * @param problem_name
  * @param row print info at particular row
  * @param col print info at particular col
  * @return MoFEMErrorCode

  */
  template <class Tag>
  MoFEMErrorCode checkMPIAIJMatrixFillIn(const std::string problem_name,
                                         int row_print, int col_print,
                                         int verb = QUIET) {
    static_assert(!std::is_same<Tag, Tag>::value, "not implemented");
    return 0;
  }

private:
  PetscLogEvent MOFEM_EVENT_createMPIAIJ;
  PetscLogEvent MOFEM_EVENT_createMPIAIJWithArrays;
  PetscLogEvent MOFEM_EVENT_createMPIAdjWithArrays;
  PetscLogEvent MOFEM_EVENT_createSeqAIJWithArrays;
  PetscLogEvent MOFEM_EVENT_checkMatrixFillIn;
};

template <>
MoFEMErrorCode MatrixManager::createMPIAIJWithArrays<PetscGlobalIdx_mi_tag>(
    const std::string name, Mat *Aij, int verb);

template <>
MoFEMErrorCode
MatrixManager::createMPIAIJ<PetscGlobalIdx_mi_tag>(const std::string name,
                                                   Mat *Aij, int verb);

template <>
MoFEMErrorCode
MatrixManager::createMPIAdjWithArrays<Idx_mi_tag>(const std::string name,
                                                  Mat *Adj, int verb);

template <>
MoFEMErrorCode MatrixManager::createSeqAIJWithArrays<PetscLocalIdx_mi_tag>(
    const std::string name, Mat *Aij, int verb);

template <>
MoFEMErrorCode
MatrixManager::checkMPIAIJWithArraysMatrixFillIn<PetscGlobalIdx_mi_tag>(
    const std::string problem_name, int row_print, int col_print, int verb);

template <>
MoFEMErrorCode MatrixManager::checkMPIAIJMatrixFillIn<PetscGlobalIdx_mi_tag>(
    const std::string problem_name, int row_print, int col_print, int verb);

} // namespace MoFEM

#endif // __MATMANAGER_HPP__a

/**
 * \defgroup mofem_mat_interface Matrix Manager
 * \brief Creating and managing matrices
 *
 * \ingroup mofem
 */
