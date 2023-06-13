/** \file DeprecatedCoreInterface.hpp
 * \brief Deprecated interface functions
 */

#ifndef __INTERFACE_DEPRECATED_HPP__
#define __INTERFACE_DEPRECATED_HPP__

/** \brief name space of MoFEM library functions and classes
 */
namespace MoFEM {

/**
 * \brief Deprecated interface functions
 * \nosubgrouping
 */
struct DeprecatedCoreInterface : public CoreInterface {

//   /** \name Interfaces */

//   /**@}*/

  /** \name Seed entities */

  /**@{*/

  /** \deprecated use BitRefManager
   * \brief seed 2D entities (Triangles entities only) in the meshset and their
   * adjacencies (only TRIs adjacencies) in a particular BitRefLevel \todo
   * Should be outsourced to separate interface, i.e. BitLevelManager
   *
   * \param EntityHandle MeshSet
   * \param BitRefLevel bitLevel
   *
   */
  DEPRECATED virtual MoFEMErrorCode
  seed_ref_level_2D(const EntityHandle meshset, const BitRefLevel &bit,
                    int verb = -1);

  /** \deprecated use BitRefManager
  * \brief seed 2D entities in the meshset and their adjacencies (only TETs
  adjacencies) in a particular BitRefLevel
  * \todo Should be outsourced to separate interface, i.e. BitLevelManager
  *
  * \param EntityHandle MeshSet
  * \param BitRefLevel bitLevel
  *
  * \brief Example:\code
  EntityHandle meshset1; //contains ent1,ent2,ent3
  BitRefLevel myLevel0;
  myLevel0.set(0);
  seed_ref_level_3D(meshset1,myLevel0);
  //refine meshset1 into meshset2 and get new ents which are ent4, ent5
  EntityHandle meshset2; //contains ent1,ent2,ent3,ent4,ent5
  BitRefLevel myLevel1;
  myLevel1.set(1);
  seed_ref_level_3D(meshset2,myLevel1);  \endcode

  * So entities 1,2,3 would be assigned to bit level 0 and 1 <br>
  * ent1[1,1,0,0,0,0,0], ent2[1,1,0,0,0,0,0], ent3[1,1,0,0,0,0,0], <br>
  * and entities 4 and 5 are assigned to bit level 1 only <br>
  * ent4[0,1,0,0,0,0,0], ent5[0,1,0,0,0,0,0] <br>
  *
  */
  DEPRECATED MoFEMErrorCode seed_ref_level_3D(const EntityHandle meshset,
                                              const BitRefLevel &bit,
                                              int verb = -1);

  /** \deprecated use BitRefManager
   * \brief seed entities in the range and their adjacencies in a particular
   * BitRefLevel \todo Should be outsourced to separate interface, i.e.
   * BitLevelManager
   */
  DEPRECATED MoFEMErrorCode seed_ref_level(const Range &ents,
                                           const BitRefLevel &bit,
                                           const bool only_tets = true,
                                           int verb = -1);

  /**
   * \brief Set partition tag to each finite element in the problem
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
  DEPRECATED MoFEMErrorCode partition_mesh(const Range &ents, const int dim,
                                           const int adj_dim, const int n_parts,
                                           int verb = -1);

  /**@}*/

  /** \name Synchronize */

  /**@{*/

  /** synchronize entity range on processors (collective)

  collective - need tu be run on all processors in communicator

  @deprecated Use Comm Interface
  */
  DEPRECATED MoFEMErrorCode synchronise_entities(Range &ent,
                                                 int verb = DEFAULT_VERBOSITY);

  /** synchronize entity range on processors (collective)
  * \ingroup mofem_field

  collective - need tu be run on all processors in communicator

  \param name field
  \param verbose level

  \deprecated Use CommInterface
  */
  DEPRECATED MoFEMErrorCode synchronise_field_entities(
      const std::string &name, int verb = DEFAULT_VERBOSITY);

  /**@}*/

};

} // namespace MoFEM

#endif // __INTERFACE_DEPRECATED_HPP__
