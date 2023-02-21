/** \file ProblemsManager.hpp
 * \brief Interface managing problems
 * \ingroup mofem_problems_manager
 *
 * Managing problems, build and partitioning.
 *
 */

#ifndef __PROBLEMSMANAGER_HPP__
#define __PROBLEMSMANAGER_HPP__

#include "UnknownInterface.hpp"

namespace MoFEM {

/**
 * \brief Problem manager is used to build and partition problems
 * \ingroup mofem_problems_manager
 *
 */
struct ProblemsManager : public UnknownInterface {

  MoFEMErrorCode query_interface(boost::typeindex::type_index type_index,
                                 UnknownInterface **iface) const;

  MoFEM::Core &cOre;
  ProblemsManager(const MoFEM::Core &core);
  virtual ~ProblemsManager() = default;

  PetscBool buildProblemFromFields; ///< If set to true, problem is build from
  /// DOFs in fields, not from DOFs on elements

  PetscBool synchroniseProblemEntities;

  MoFEMErrorCode getOptions();

  /**
   * \brief Set partition tag to each finite element in the problem
   * \deprecated Use CommInterface
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
  DEPRECATED MoFEMErrorCode partitionMesh(const Range &ents, const int dim,
                                          const int adj_dim, const int n_parts,
                                          Tag *th_vertex_weights = nullptr,
                                          Tag *th_edge_weights = nullptr,
                                          Tag *th_part_weights = nullptr,
                                          int verb = VERBOSE,
                                          const bool debug = false);

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
   * \ingroup mofem_problems_manager
   * 
   * @param  out_name problem
   * @param  fields_row  vector of fields composing problem
   * @param  fields_col  vector of fields composing problem
   * @param  main_problem main problem
   * @return              error code
   */
  MoFEMErrorCode buildSubProblem(
      const std::string out_name, const std::vector<std::string> &fields_row,
      const std::vector<std::string> &fields_col,
      const std::string main_problem, const bool square_matrix = true,
      const map<std::string, boost::shared_ptr<Range>> *entityMapRow = nullptr,
      const map<std::string, boost::shared_ptr<Range>> *entityMapCol = nullptr,
      int verb = VERBOSE);

  /**
   * \brief build composite problem
   * \ingroup mofem_problems_manager
   * 
   * @param  out_name         name of build problem
   * @param  add_row_problems vector of add row problems
   * @param  add_col_problems vector of add col problems
   * @param  square_matrix    true if structurally squared matrix
   * @param  verb             verbosity level
   * @return                  error code
   */
  MoFEMErrorCode
  buildComposedProblem(const std::string out_name,
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
   * @param name 
   * @param part_from_moab 
   * @param low_proc 
   * @param hi_proc 
   * @param verb 
   * @return MoFEMErrorCode 
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
   * \create add entities of finite element in the problem
   * \ingroup mofem_problems_manager
   *
   * \note Meshset entity has to be crated
   * 
   */

  /**
   * @brief create add entities of finite element in the problem
   *
   * @note Meshset entity has to be crated
   * 
   * @param prb_name name of the problem
   * @param fe_name name of the entity
   * @param meshset pointer meshset handle
   * @return MoFEMErrorCode 
   */
  MoFEMErrorCode getFEMeshset(const std::string prb_name,
                              const std::string &fe_name,
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
                                          const std::string &fe_name,
                                          PetscLayout *layout) const;

  /**
   * @brief Remove DOFs from problem
   * @ingroup mofem_problems_manager
   *
   * Remove DOFs from problem which are on entities on the given range and given
   * field name. On the finite element level, DOFs can be still accessed however
   * local PETSc indices and global PETSc indices are marked with the index -1.
   *
   * @note If the index is marked -1 it is not assembled and dropped by
   * VecSetValues and MatSetValues.
   *
   * @todo Not yet implemented update for AO maps and IS ranges if removed
   * entities in composite problem or sub-problem
   *
   * @param problem_name name of the problem
   * @param field_name name of the field
   * @param ents entities on which DOFs are removed
   * @param lo_coeff low dof coefficient (rank)
   * @param hi_coeff high dof coefficient (rank)
   * @param verb  verbosity level
   * @param debug to debug and seek for inconsistencies set to true
   * @return MoFEMErrorCode
   */
  MoFEMErrorCode removeDofsOnEntities(
      const std::string problem_name, const std::string field_name,
      const Range ents, const int lo_coeff = 0,
      const int hi_coeff = MAX_DOFS_ON_ENTITY, const int lo_order = 0,
      const int hi_order = 100, int verb = VERBOSE, const bool debug = false);

  /**
   * @copydoc removeDofsOnEntities
   *
   * \note Use this function for nondistributed meshes
   */
  MoFEMErrorCode removeDofsOnEntitiesNotDistributed(
      const std::string problem_name, const std::string field_name,
      const Range ents, const int lo_coeff = 0,
      const int hi_coeff = MAX_DOFS_ON_ENTITY, const int lo_order = 0,
      const int hi_order = 100, int verb = VERBOSE, const bool debug = false);

  /**
   * @brief Remove DOFs from problem by bit ref level
   * @ingroup mofem_problems_manager
   *
   * See for more detail other implementation for removeDofsOnEntities.
   *
   * @param problem_name name of the problem
   * @param field_name name of the field
   * @param bit_ref_level bit ref level on which DOFs are removed
   * @param bit_ref_mask bit ref mask on which DOFs are removed
   * @param ents_ptr filter entities with given bit and mask
   * @param lo_coeff low dof coefficient (rank)
   * @param hi_coeff high dof coefficient (rank)
   * @param verb  verbosity level
   * @param debug to debug and seek for inconsistencies set to true
   * @return MoFEMErrorCode
   */
  MoFEMErrorCode removeDofsOnEntities(
      const std::string problem_name, const std::string field_name,
      const BitRefLevel bit_ref_level, const BitRefLevel bit_ref_mask,
      Range *ents_ptr = nullptr, const int lo_coeff = 0,
      const int hi_coeff = MAX_DOFS_ON_ENTITY, const int lo_order = 0,
      const int hi_order = 100, int verb = VERBOSE, const bool debug = false);

  /**
   * @copydoc removeDofsOnEntities
   *
   * \note Use this function for nondistributed meshes
   */
  MoFEMErrorCode removeDofsOnEntitiesNotDistributed(
      const std::string problem_name, const std::string field_name,
      const BitRefLevel bit_ref_level, const BitRefLevel bit_ref_mask,
      Range *ents_ptr = nullptr, const int lo_coeff = 0,
      const int hi_coeff = MAX_DOFS_ON_ENTITY, const int lo_order = 0,
      const int hi_order = 100, int verb = VERBOSE, const bool debug = false);

  enum MarkOP { OR, AND };

  /**
   * @brief Create vector with marked indices 
   * 
   * Vector with local DOFs marked by entities 
   * 
   * 
   * @param problem_name 
   * @param row 
   * @param ents 
   * @param marker 
   * @return MoFEMErrorCode 
   */
  MoFEMErrorCode markDofs(const std::string problem_name, RowColData rc,
                          const enum MarkOP op, const Range ents,
                          std::vector<unsigned char> &marker) const;

  /**
   * @deprecated use function with MarkOP
   */
  inline DEPRECATED MoFEMErrorCode
  markDofs(const std::string problem_name, RowColData rc, const Range ents,
           std::vector<unsigned char> &marker) const {
    return markDofs(problem_name, rc, MarkOP::OR, ents, marker);
  }

  /**
   * @brief Mark DOFs
   *
   * @param problem_name
   * @param rc
   * @param field_name
   * @param lo
   * @param hi
   * @param op
   * @param marker
   * @return MoFEMErrorCode
   */
  MoFEMErrorCode modifyMarkDofs(const std::string problem_name, RowColData rc,
                                const std::string field_name, const int lo,
                                const int hi, const enum MarkOP op,
                                const unsigned char c,
                                std::vector<unsigned char> &marker) const;

  /**
   * @brief add empty block to problem
   *
   * MatrixManager assumes that all blocks, i.e. all fields combinations are non
   * zero. This is not always the case, to optimise code and reduce memory usage
   * user can specifi which blocks are empty.
   *
   * @param problem_name problem name
   * @param row_field row filed name
   * @param col_field col field name
   * @return MoFEMErrorCode
   */
  MoFEMErrorCode addFieldToEmptyFieldBlocks(const std::string problem_name,
                                            const std::string row_field,
                                            const std::string col_field) const;

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
