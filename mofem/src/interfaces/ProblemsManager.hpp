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

  static const MOFEMuuid IDD_MOFEMProblemsManager = MOFEMuuid( BitIntefaceId(PROBLEMSMANAGER_INTERFACE) );

  /**
   * \brief Problem manager is used to build and partition problems
   * \mofem_problems_manager
   *
   */
  struct ProblemsManager: public UnknownInterface {

    PetscErrorCode queryInterface(const MOFEMuuid& uuid, UnknownInterface** iface);

    MoFEM::Core& cOre;
    ProblemsManager(const MoFEM::Core& core);

    /**
    * \brief Destructor
    */
    ~ProblemsManager();

    /**
     * \brief Set partition tag to each finite element in the problem
     * \ingroup mofem_problems_manager
     *
     * This will use one of the mesh partitioning programs available from PETSc
     * See <http://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/Mat/MatPartitioningType.html>
     *
     * @param  ents        Entities to partition
     * @param  dim         Dimension of entities to partition
     * @param  adj_dim     Adjacency dimension
     * @param  n_parts     Number of partitions
     * @param  verb        Verbosity level
     * @return             Error code
     */
    PetscErrorCode partitionMesh(
      const Range &ents,const int dim,const int adj_dim,const int n_parts,int verb = 1
    );

    /** \brief build problem data structures
     * \ingroup mofem_problems_manager
     *
     * \note If square_matrix is set to true, that indicate that problem is structurally
     * symmetric, i.e. rows and columns have the same dofs and are indexed in the same
     * way.
     *
     * @param  name          problem name
     * @param  square_matrix structurally symmetric problem
     * @param  verb          verbosity level
     * @return               error code
     *
     */
     PetscErrorCode buildProblem(const std::string &name,const bool square_matrix,int verb = 1);

     /** \brief build problem data structures
      * \ingroup mofem_problems_manager
      *
      * \note If square_matrix is set to true, that indicate that problem is structurally
      * symmetric, i.e. rows and columns have the same dofs and are indexed in the same
      * way.
      *
      * @param  problem pointer
      * @param  square_matrix structurally symmetric problem
      * @param  verb          verbosity level
      * @return               error code
      *
      */
    PetscErrorCode buildProblem(Problem *problem_ptr,const bool square_matrix,int verb = 1);

    /** \brief build problem data structures, assuming that mesh is distributed (collective)
     * \ingroup mofem_problems_manager

     Mesh is distributed, that means that each processor keeps only own part of
     the mesh and shared entities.

     Collective - need to be run on all processors in communicator, i.e. each
     function has to call this function.

     */
    PetscErrorCode buildProblemOnDistributedMesh(
      const std::string &name,const bool square_matrix,int verb = 1
    );

    /** \brief build problem data structures, assuming that mesh is distributed (collective)
     * \ingroup mofem_problems_manager

     Mesh is distributed, that means that each processor keeps only own part of
     the mesh and shared entities.

     Collective - need to be run on all processors in communicator, i.e. each
     function has to call this function.

     */
    PetscErrorCode buildProblemOnDistributedMesh(
      Problem *problem_ptr,
      const bool square_matrix = true,
      int verb = 1
    );

    /**
     * \brief build sub problem
     * @param  out_name problem
     * @param  fields_row  vector of fields composing problem
     * @param  fields_col  vector of fields composing problem
     * @param  main_problem main problem
     * @return              error code
     */
    PetscErrorCode buildSubProblem(
      const std::string &out_name,
      const std::vector<std::string> &fields_row,
      const std::vector<std::string> &fields_col,
      const std::string &main_problem,
      const bool square_matrix = true,
      int verb = 1
    );

    /**
     * \brief build composite problem
     * @param  out_name         name of build problem
     * @param  add_row_problems vector of add row problems
     * @param  add_col_problems vector of add col problems
     * @param  square_matrix    true if structurally squared matrix
     * @param  verb             verbosity level
     * @return                  error code
     */
    PetscErrorCode buildCompsedProblem(
      const std::string &out_name,
      const std::vector<std::string> add_row_problems,
      const std::vector<std::string> add_col_problems,
      const bool square_matrix = true,
      int verb = 1
    );

    /**
      * \brief build indexing and partition problem inheriting indexing and partitioning from two other problems
      * \ingroup mofem_problems_manager
      *
      * \param name problem name
      * \param problem_for_rows problem used to index rows
      * \param copy_rows just copy rows dofs
      * \param problem_for_cols problem used to index cols
      * \param copy_cols just copy cols dofs
      *
      * If copy_rows/copy_cols is set to false only partition is copied between problems.
      *
      */
    PetscErrorCode inheretPartition(
      const std::string &name,
      const std::string &problem_for_rows,
      bool copy_rows,
      const std::string &problem_for_cols,
      bool copy_cols,
      int verb = 1
    );

    /** \brief partition problem dofs
     * \ingroup mofem_problems_manager
     *
     * \param name problem name
     */
    PetscErrorCode partitionSimpleProblem(const std::string &name,int verb = 1);

    /** \brief partition problem dofs (collective)
     * \ingroup mofem_problems_manager
     *
     * \param name problem name
     */
    PetscErrorCode partitionProblem(const std::string &name,int verb = 1);

    PetscErrorCode printPartitionedProblem(
      const Problem *problem_ptr,int verb = 1
    );

    PetscErrorCode debugPartitionedProblem(
      const Problem *problem_ptr,int verb = 1
    );

    /** \brief partition finite elements
     * \ingroup mofem_problems_manager
     *
     * Function which partition finite elements based on dofs partitioning.<br>
     * In addition it sets information about local row and cols dofs at given element on partition.
     *
     * \param name problem name
     */
    PetscErrorCode partitionFiniteElements(const std::string &name,
      bool part_from_moab = false,int low_proc = -1,int hi_proc = -1,int verb = 1
    );

    /** \brief determine ghost nodes
     * \ingroup mofem_problems_manager
     *
     * \param name problem name
     */
    PetscErrorCode partitionGhostDofs(const std::string &name,int verb = 1);


  private:

    PetscLogEvent USER_EVENT_ProblemsManager;


  };



}

#endif //__PROBLEMSMANAGER_HPP__


/***************************************************************************//**
 * \defgroup mofem_problems_manager Problems manager
 * \brief Adding and managing problems
 *
 * \ingroup mofem
 ******************************************************************************/
