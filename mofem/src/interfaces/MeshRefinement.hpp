/** \file MeshRefinement.hpp
 * \brief Interface for mesh refinement
 *
 * \ingroup mofem_refiner
 */

#ifndef __MESHREFINE_HPP__
#define __MESHREFINE_HPP__

namespace MoFEM {

/** \brief Mesh refinement interface

  Currently this class is abstraction to Core interface. In future should be
  outsourced as independent interface.

  \bug Not working on partitioned meshes
  \bug Need to be implemented as a stand alone interface not as a part of core
  structure which should be only basic database
  \bug If outsourced, class member functions should follow name convention
  \bug Spelling mistakes will be corrected with names fix to follow name
  convetion

  \ingroup mofem_refiner
  */
struct MeshRefinement : public UnknownInterface {

  MoFEMErrorCode query_interface(boost::typeindex::type_index type_index,
                                 UnknownInterface **iface) const;

  MoFEM::Core &cOre;
  MeshRefinement(const MoFEM::Core &core);

  virtual ~MeshRefinement() = default;

  /**
   * \brief make vertices in the middle of edges in meshset and add them to
   * refinement levels defined by bit
   *
   * Takes entities fromm meshsets and queried recursively (get entities from
   * meshsets in meshsets, usually have to be used for CUBIT meshset).
   * If meshset does not contain any edges, get entities in dimension 3 and get
   * edge adjacencies.
   *
   * \param EntityHandle meshset
   * \param BitRefLevel bitLevel
   * \param recursive If true, meshsets containing meshsets are queried
   * recursively.  Returns the contents of meshsets, but not the meshsets
   * themselves if true.
   */
  MoFEMErrorCode addVerticesInTheMiddleOfEdges(const EntityHandle meshset,
                                               const BitRefLevel &bit,
                                               const bool recursive = false,
                                               int verb = QUIET,
                                               EntityHandle start_v = 0);

  /**
   * \brief make vertices in the middle of edges in meshset and add them to
   * Refinement levels defined by bit
   *
   * Takes entities from meshsets and queried recursively (get entities from
   * meshsets in meshsets, usually have to be used for CUBIT meshset).
   * If meshset does not contain any edges, get entities in dimension 3 and get
   * edge adjacencies.
   *
   * \param Range consisting edges for refine
   * \param BitRefLevel bitLevel
   * \param recursive If true, meshsets containing meshsets are queried
   * recursively.  Returns the contents of meshsets, but not the meshsets
   * themselves if true.
   */
  MoFEMErrorCode addVerticesInTheMiddleOfEdges(const Range &edges,
                                               const BitRefLevel &bit,
                                               int verb = QUIET,
                                               EntityHandle start_v = 0);

  /**\brief refine TET in the meshset
   *
   * \param EntityHandle meshset
   * \param BitRefLevel bitLevel
   * \param verb verbosity level
   */
  MoFEMErrorCode refineTets(const EntityHandle meshset, const BitRefLevel &bit,
                            int verb = QUIET, const bool debug = false);

  /**\brief refine TET in the meshset
   *
   * \param Range of tets to refine
   * \param BitRefLevel bitLevel
   * \param BitRefLevel bitLevel
   * \param verb verbosity level
   */
  MoFEMErrorCode refineTets(const Range &tets, const BitRefLevel &bit,
                            int verb = QUIET, const bool debug = false);

  /**\brief refine TET in the meshset
   *
   * \param Range of tets to refine
   * \param BitRefLevel bitLevel
   * \param BitRefLevel bitLevel
   * \param verb verbosity level
   */
  MoFEMErrorCode refineTetsHangingNodes(const Range &tets,
                                        const BitRefLevel &bit,
                                        int verb = QUIET,
                                        const bool debug = false);

  /**\brief refine TET in the meshset
   *
   * \param Range of tets to refine
   * \param BitRefLevel bitLevel
   * \param BitRefLevel bitLevel
   * \param verb verbosity level
   */
  MoFEMErrorCode refineTetsHangingNodes(const EntityHandle meshset,
                                        const BitRefLevel &bit,
                                        int verb = QUIET,
                                        const bool debug = false);

  /**\brief refine PRISM in the meshset
   *
   * \param EntityHandle meshset
   * \param BitRefLevel bitLevel
   */
  MoFEMErrorCode refinePrisms(const EntityHandle meshset,
                              const BitRefLevel &bit, int verb = QUIET);

  /**\brief refine meshset, i.e. add child of refined entities to meshset
   *
   * \param EntityHandle meshset where to save the child refined entities
   * \param BitRefLevel bitLevel
   * \param recursive If true, meshsets containing meshsets are queried
   * recursively.  Returns the contents of meshsets, but not the meshsets
   * themselves if true.
   */
  MoFEMErrorCode refineMeshset(const EntityHandle meshset,
                               const BitRefLevel &bit,
                               const bool recursive = false, int verb = QUIET);

  /**\brief refine triangles in the meshset
   *
   * \param EntityHandle meshset
   * \param BitRefLevel bitLevel
   * \param verb verbosity level
   */
  MoFEMErrorCode refineTris(const EntityHandle meshset, const BitRefLevel &bit,
                            int verb = QUIET, const bool debug = false);

  /**\brief refine TRI in the meshset
   *
   * \param meshset of entities to refine
   * \param BitRefLevel bit level of created entities
   * \param verb verbosity level
   */
  MoFEMErrorCode refineTris(const Range &tris, const BitRefLevel &bit,
                            int verb = QUIET, const bool debug = false);

  /**\brief refine TRI in the meshset
   *
   * \param Range of entities to refine
   * \param BitRefLevel bit level of created entities
   * \param verb verbosity level
   */
  MoFEMErrorCode refineTrisHangingNodes(const EntityHandle meshset,
                                        const BitRefLevel &bit,
                                        int verb = QUIET,
                                        const bool debug = false);

  /**\brief refine TRI in the meshset
   *
   * \param Range of entities to refine
   * \param BitRefLevel bit level of created entities
   * \param verb verbosity level
   */
  MoFEMErrorCode refineTrisHangingNodes(const Range &tris,
                                        const BitRefLevel &bit,
                                        int verb = QUIET,
                                        const bool debug = false);

private:
  struct SetParent {
    map<EntityHandle, EntityHandle> parentsToChange;
    MoFEMErrorCode operator()(const EntityHandle ent, const EntityHandle parent,
                              const RefEntity_multiIndex *ref_ents_ptr,
                              MoFEM::Core &cOre);

    MoFEMErrorCode operator()(const RefEntity_multiIndex *ref_ents_ptr);
  };


  /**
   * @brief Functions setting edges for refinemnt on enetity level 
   * 
   */
  using SetEdgeBitsFun = boost::function<

      MoFEMErrorCode(moab::Interface &moab,
                     RefEntity_multiIndex_view_by_ordered_parent_entity
                         &ref_parent_ents_view,
                     EntityHandle tet, BitRefEdges &parent_edges_bit,
                     EntityHandle *edge_new_nodes, int *split_edges

                     )>;

  /**\brief refine TET in the meshset
   *
   * \param Range of tets to refine
   * \param BitRefLevel bitLevel
   * \param verb verbosity level
   */
  MoFEMErrorCode refineTets(const Range &tets, const BitRefLevel &bit,
                            SetEdgeBitsFun set_edge_bits, int verb,
                            const bool debug);

  MoFEMErrorCode refineTris(const Range &tris, const BitRefLevel &bit,
                            SetEdgeBitsFun set_edge_bits, int verb,
                            const bool debug);
};

} // namespace MoFEM

#endif // __MESHREFINE_HPP__

/**
 * \defgroup mofem_refiner MeshRefinement
 * \brief Refine mesh by splitting edges
 *
 */
