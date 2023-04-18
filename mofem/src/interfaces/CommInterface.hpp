/** \file CommInterface.hpp
 * \brief Interface for communication functions
 * \ingroup mofem_comm
 *
 * Functions used to communicate, share entities, share data, etc.
 *
 */

#ifndef __COMMINTERFACE_HPP__
#define __COMMINTERFACE_HPP__

#include "UnknownInterface.hpp"

namespace MoFEM {

/**
 * \brief Managing BitRefLevels
 * \ingroup mofem_bit_ref
 * \nosubgrouping
 */
struct CommInterface : public UnknownInterface {

  MoFEMErrorCode query_interface(boost::typeindex::type_index type_index,
                                 UnknownInterface **iface) const;

  MoFEM::Core &cOre;
  bool dEbug;

  CommInterface(const MoFEM::Core &core);

  /**
   * \brief Destructor
   */
  ~CommInterface() = default;

  /** \name Make elemnts multishared */

  /**@{*/

  /**
  * \brief resolve shared entities for finite elements in the problem
  * \ingroup mofem_problems

  * @param  problem_ptr  problem pointer
  * @param  fe_name     finite element name
  * @param  verb        verbosity level
  * @return             error code
  *
  * This allows for tag reduction or tag exchange, f.e.

  \code
  CHKERR m_field.resolveSharedFiniteElements(problem_ptr,"SHELL_ELEMENT");
  Tag th;
  CHKERR mField.get_moab().tag_get_handle("ADAPT_ORDER",th);
  CHKERR ParallelComm* pcomm =
  ParallelComm::get_pcomm(&mField.get_moab(),MYPCOMM_INDEX);
  CHKERR pcomm->exchange_tags(th,prisms);
  \endcode

  *
  */
  MoFEMErrorCode resolveSharedFiniteElements(const Problem *problem_ptr,
                                             const std::string &fe_name,
                                             int verb = DEFAULT_VERBOSITY);

  /**
  * \brief resolve shared entities for finite elements in the problem
  * \ingroup mofem_problems

  * @param  name        problem name
  * @param  fe_name     finite element name
  * @param  verb        verbosity level
  * @return             error code
  *
  * This allows for tag reduction or tag exchange, f.e.

  \code
  CHKERR m_field.resolveSharedFiniteElements(problem_ptr,"SHELL_ELEMENT");
  Tag th;
  CHKERR mField.get_moab().tag_get_handle("ADAPT_ORDER",th);
  ParallelComm* pcomm =
  ParallelComm::get_pcomm(&mField.get_moab(),MYPCOMM_INDEX);
  // CHKERR pcomm->reduce_tags(th,MPI_SUM,prisms);
  CHKERR pcomm->exchange_tags(th,prisms);
  \endcode

  *
  */
  MoFEMErrorCode resolveSharedFiniteElements(const std::string name,
                                             const std::string &fe_name,
                                             int verb = DEFAULT_VERBOSITY);

  /** \name Make entities multishared */

  /**
   * @brief make entities from proc 0 shared on all proc
   *
   * \note collective - need tu be run on all processors in communicator
   *
   * @param entities
   * @param num_entities
   * @param my_proc default proc id to share from
   * @param verb
   * @return MoFEMErrorCode
   */
  MoFEMErrorCode makeEntitiesMultishared(const EntityHandle *entities,
                                         const int num_entities,
                                         const int owner_proc = 0,
                                         int verb = DEFAULT_VERBOSITY);

  /**
   * @brief make entities from proc 0 shared on all proc
   *
   * \note collective - need tu be run on all processors in communicator
   *
   * @param entities
   * @param my_proc default proc id to share from
   * @param verb
   * @return MoFEMErrorCode
   */
  MoFEMErrorCode makeEntitiesMultishared(Range &entities,
                                         const int owner_proc = 0,
                                         int verb = DEFAULT_VERBOSITY);

  /**
   * @brief make field entities multi shared
   *
   * \note collective - need tu be run on all processors in communicator
   *
   * @param field_name
   * @param owner_proc
   * @param verb
   * @return MoFEMErrorCode
   */
  MoFEMErrorCode makeFieldEntitiesMultishared(const std::string field_name,
                                              const int owner_proc = 0,
                                              int verb = DEFAULT_VERBOSITY);

  /**
   * @brief Exchange field data
   *
   * Exchange field for all shared and ghosted entities. This function should be
   * called collectively over the communicator for this ParallelComm. If the
   * entities vector is empty, all shared entities participate in the exchange.
   * If a proc has no owned entities this function must still be called since it
   * is collective.
   *
   * \note collective - need tu be run on all processors in communicator
   *
   * \todo It is not working if field has entities diffrent than vertices.
   *
   * @param verb
   * @param field_name @return MoFEMErrorCode
   */
  MoFEMErrorCode exchangeFieldData(const std::string field_name,
                                   int verb = DEFAULT_VERBOSITY);

  /**@}*/

  /** \name Synchronize entities (Following functions in future will be
   * deprecated) */

  /**@{*/

  /** synchronize entity range on processors (collective)
   *
   * \note collective - need tu be run on all processors in communicator
   *
   */

  /**
   * @brief synchronize entity range on processors (collective)
   * 
   * \note collective - need tu be run on all processors in communicator
   * 
   * @param ent ents to send and received
   * @param received_ents pointer to map with received entities
   * @param verb 
   * @return MoFEMErrorCode 
   */
  MoFEMErrorCode synchroniseEntities(Range &ent,
                                     std::map<int, Range> *received_ents,
                                     int verb = DEFAULT_VERBOSITY);

   /**
   * @brief synchronize entity range on processors (collective)
   * 
   * \note collective - need tu be run on all processors in communicator
   * 
   * @param ent ents to send and received
   * @param verb 
   * @return MoFEMErrorCode 
   */ 
  MoFEMErrorCode synchroniseEntities(Range &ent, int verb = DEFAULT_VERBOSITY);

  /** synchronize entity range on processors (collective)
   * \ingroup mofem_field
   *
   *  \note collective - need tu be run on all processors in communicator
   *
   *  \param name field
   *  \param verbose level
   *
   */
  MoFEMErrorCode synchroniseFieldEntities(const std::string name,
                                          int verb = DEFAULT_VERBOSITY);

  /**
   * @brief Synchronise parent entities
   *
   *  \note collective - need tu be run on all processors in communicator
   *
   * Exchange parent entity tag and bitref of entity. Note thar parent handle
   * can be different on each processor.
   *
   * @param ent
   * @param verb
   * @return MoFEMErrorCode
   */
  MoFEMErrorCode resolveParentEntities(const Range &ent,
                                       int verb = DEFAULT_VERBOSITY);

  /**@}*/

  /**@*/

  /** \name Read load and boradcoast */

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


  /**@}*/

};
} // namespace MoFEM

#endif //__COMMINTERFACE_HPP__

/**
 * \defgroup mofem_comm Comm intrface
 * \brief Comm interface
 *
 * \ingroup mofem
 */