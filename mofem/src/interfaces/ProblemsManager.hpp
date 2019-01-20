/** \file ProblemsManager.hpp
 * \brief Interface managing problems
 * \ingroup mofem_problems_manager
 *
 * Managing problems, build and partitioning.
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

#ifndef __PROBLEMSMANAGER_HPP__
#define __PROBLEMSMANAGER_HPP__

#include "UnknownInterface.hpp"

namespace MoFEM {

static const MOFEMuuid IDD_MOFEMProblemsManager =
    MOFEMuuid(BitIntefaceId(PROBLEMSMANAGER_INTERFACE));

/**
 * \brief Problem manager is used to build and partition problems
 * \mofem_problems_manager
 *
 */
struct ProblemsManager : public UnknownInterface {

  MoFEMErrorCode query_interface(const MOFEMuuid &uuid,
                                 UnknownInterface **iface) const;

  MoFEM::Core &cOre;
  ProblemsManager(const MoFEM::Core &core);

  /**
   * \brief Destructor
   */
  ~ProblemsManager();

  PetscBool buildProblemFromFields; ///< If set to true, problem is build from
  /// DOFs in fields, not from DOFs on elements

  PetscBool synchroniseProblemEntities;

  MoFEMErrorCode getOptions();

  /**
   * \brief Set partition tag to each finite element in the problem
   * \ingroup mofem_problems_manager
   *
   * This will use one of the mesh partitioning programs available from PETSc
   * See
   * <http://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/Mat/MatPartitioningType.html>
   *
   * @param  ents        Entities to partition
   * @param  dim         Dimension of entities to partition
   * @param  adj_dim     Adjacency dimension
   * @param  n_parts     Number of partitions
   * @param  verb        Verbosity level
   * @return             Error code
   */
  MoFEMErrorCode partitionMesh(const Range &ents, const int dim,
                               const int adj_dim, const int n_parts,
                               Tag *th_vertex_weights = nullptr,
                               Tag *th_edge_weights = nullptr,
                               Tag *th_part_weights = nullptr,
                               int verb = VERBOSE, const bool debug = false);

  /** \brief build problem data structures
   * \ingroup mofem_problems_manager
   *
   * \note If square_matrix is set to true, that indicate that problem is
   * structurally
   * symmetric, i.e. rows and columns have the same dofs and are indexed in the
   * same
   * way.
   *
   * @param  name          problem name
   * @param  square_matrix structurally symmetric problem
   * @param  verb          verbosity level
   * @return               error code
   *
   */
  MoFEMErrorCode buildProblem(const std::string name, const bool square_matrix,
                              int verb = VERBOSE);

  /** \brief build problem data structures
   * \ingroup mofem_problems_manager
   *
   * \note If square_matrix is set to true, that indicate that problem is
   * structurally
   * symmetric, i.e. rows and columns have the same dofs and are indexed in the
   * same
   * way.
   *
   * @param  problem pointer
   * @param  square_matrix structurally symmetric problem
   * @param  verb          verbosity level
   * @return               error code
   *
   */
  MoFEMErrorCode buildProblem(Problem *problem_ptr, const bool square_matrix,
                              int verb = VERBOSE);

  /** \brief build problem data structures, assuming that mesh is distributed
   (collective)
   * \ingroup mofem_problems_manager

   Mesh is distributed, that means that each processor keeps only own part of
   the mesh and shared entities.

   Collective - need to be run on all processors in communicator, i.e. each
   function has to call this function.

   */
  MoFEMErrorCode buildProblemOnDistributedMesh(const std::string name,
                                               const bool square_matrix,
                                               int verb = VERBOSE);

  /** \brief build problem data structures, assuming that mesh is distributed
   (collective)
   * \ingroup mofem_problems_manager

   Mesh is distributed, that means that each processor keeps only own part of
   the mesh and shared entities.

   Collective - need to be run on all processors in communicator, i.e. each
   function has to call this function.

   */
  MoFEMErrorCode buildProblemOnDistributedMesh(Problem *problem_ptr,
                                               const bool square_matrix = true,
                                               int verb = VERBOSE);

  /**
   * \brief build sub problem
   * @param  out_name problem
   * @param  fields_row  vector of fields composing problem
   * @param  fields_col  vector of fields composing problem
   * @param  main_problem main problem
   * @return              error code
   */
  MoFEMErrorCode buildSubProblem(const std::string out_name,
                                 const std::vector<std::string> &fields_row,
                                 const std::vector<std::string> &fields_col,
                                 const std::string main_problem,
                                 const bool square_matrix = true,
                                 int verb = VERBOSE);

  /**
   * \brief build composite problem
   * @param  out_name         name of build problem
   * @param  add_row_problems vector of add row problems
   * @param  add_col_problems vector of add col problems
   * @param  square_matrix    true if structurally squared matrix
   * @param  verb             verbosity level
   * @return                  error code
   */
  MoFEMErrorCode
  buildCompsedProblem(const std::string out_name,
                      const std::vector<std::string> add_row_problems,
                      const std::vector<std::string> add_col_problems,
                      const bool square_matrix = true, int verb = 1);

  /**
   * \brief build indexing and partition problem inheriting indexing and
   * partitioning from two other problems
   * \ingroup mofem_problems_manager
   *
   * \param name problem name
   * \param problem_for_rows problem used to index rows
   * \param copy_rows just copy rows dofs
   * \param problem_for_cols problem used to index cols
   * \param copy_cols just copy cols dofs
   *
   * If copy_rows/copy_cols is set to false only partition is copied between
   * problems.
   *
   */
  MoFEMErrorCode inheritPartition(const std::string name,
                                  const std::string problem_for_rows,
                                  bool copy_rows,
                                  const std::string problem_for_cols,
                                  bool copy_cols, int verb = VERBOSE);

  /** \brief partition problem dofs
   * \ingroup mofem_problems_manager
   *
   * \param name problem name
   */
  MoFEMErrorCode partitionSimpleProblem(const std::string name,
                                        int verb = VERBOSE);

  /** \brief partition problem dofs (collective)
   * \ingroup mofem_problems_manager
   *
   * \param name problem name
   */
  MoFEMErrorCode partitionProblem(const std::string name, int verb = VERBOSE);

  MoFEMErrorCode printPartitionedProblem(const Problem *problem_ptr,
                                         int verb = VERBOSE);

  MoFEMErrorCode debugPartitionedProblem(const Problem *problem_ptr,
                                         int verb = VERBOSE);

  /** \brief partition finite elements
   * \ingroup mofem_problems_manager
   *
   * Function which partition finite elements based on dofs partitioning.<br>
   * In addition it sets information about local row and cols dofs at given
   * element on partition.
   *
   * \param name problem name
   */
  MoFEMErrorCode partitionFiniteElements(const std::string name,
                                         bool part_from_moab = false,
                                         int low_proc = -1, int hi_proc = -1,
                                         int verb = VERBOSE);

  /** \brief determine ghost nodes
   * \ingroup mofem_problems_manager
   * \param name problem name
   *
   * DOFs are ghost dofs if are used by elements on given partitition, but not
   * owned by that partitition.
   *
   */
  MoFEMErrorCode partitionGhostDofs(const std::string name, int verb = VERBOSE);

  /** \brief determine ghost nodes on distributed meshes
   * \ingroup mofem_problems_manager
   * \param name problem name
   *
   * It is very similar for partitionGhostDofs, however this explits that mesh
   * is
   * distributed.
   *
   * DOFs are ghosted if are shered but not owned by partition.
   *
   */
  MoFEMErrorCode partitionGhostDofsOnDistributedMesh(const std::string name,
                                                     int verb = VERBOSE);

  /**
   * \create meshset problem finite elements
   */
  MoFEMErrorCode getFEMeshset(const std::string prb_name,
                              const std::string fe_name,
                              EntityHandle *meshset) const;

  /**
   * \brief Get layout of elements in the problem
   * \ingroup mofem_problems_manager
   *
   * In layout is stored information how many elements is on each processor, for
   * more information look int petsc documentation
   * <http://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/IS/PetscLayoutCreate.html#PetscLayoutCreate>
   *
   * @param  name    problem name
   * @param  fe_name finite elements
   * @param  layout  layout
   * @param  verb    verbosity level
   * @return         error code
   */
  MoFEMErrorCode getProblemElementsLayout(const std::string name,
                                          const std::string fe_name,
                                          PetscLayout *layout) const;

  /**
   * @brief Remove DOFs from problem
   *
   * Remove DOFs from problem which are on entities on the given range and given
   * field name. On the finite element level, DOFs can be still accessed however
   * local PETSc indices and global PETSc indices are marked with the index -1.
   *
   * @note If the index is marked -1 it is not assembled and dropped by
   * VecSetValues and MatSetValues.
   *
   * @param problem_name name of the problem
   * @param field_name name of the field
   * @param ents entities on which DOFs are removed
   * @param verb
   * @return MoFEMErrorCode
   */
  MoFEMErrorCode removeDofsOnEntities(const std::string problem_name,
                                      const std::string field_name,
                                      const Range ents, int verb = VERBOSE);

private:
  PetscLogEvent MOFEM_EVENT_ProblemsManager;
};
} // namespace MoFEM

#endif //__PROBLEMSMANAGER_HPP__

/**
 * \defgroup mofem_problems_manager ProblemsManager
 * \brief Adding and managing problems
 *
 * \ingroup mofem
 */
