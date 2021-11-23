/* \file EntityRefine.hpp
 * --------------------------------------------------------------
 * \brief Tetrahedral refinement algorithm
 *
 * It is based on \cite ruprecht1998scheme
 *
 * \todo tet Refinement should be rewritten, with better error control, and more
 * clear refinement patterns.
 */

#ifndef __ENTITYREFINE_HPP__
#define __ENTITYREFINE_HPP__

namespace MoFEM {

// TETS
void tet_type_6(moab::Interface &moab, const EntityHandle *conn,
                const EntityHandle *edge_new_nodes,
                EntityHandle *new_tets_conn);
int tet_type_5(moab::Interface &moab, const EntityHandle *conn,
               const EntityHandle *edge_new_nodes, EntityHandle *new_tets_conn);
int tet_type_4(const EntityHandle *conn, const int *split_edges,
               const EntityHandle *edge_new_nodes, EntityHandle *new_tets_conn);
int tet_type_3(const EntityHandle *conn, const int *split_edges,
               const EntityHandle *edge_new_nodes, EntityHandle *new_tets_conn);
int tet_type_2(const EntityHandle *conn, const int *split_edges,
               const EntityHandle *edge_new_nodes, EntityHandle *new_tets_conn);
void tet_type_1(const EntityHandle *conn, const int split_edge,
                const EntityHandle edge_new_node, EntityHandle *new_tets_conn);

// TRIS
MoFEMErrorCode tri_type_3(const EntityHandle *conn,
                          const BitRefEdges split_edges,
                          const EntityHandle *edge_new_nodes,
                          EntityHandle *new_tris_conn);

// PRISM
MoFEMErrorCode prism_type_1(const EntityHandle *conn,
                            const BitRefEdges split_edges,
                            const EntityHandle *edge_new_nodes,
                            EntityHandle *new_prism_conn);
MoFEMErrorCode prism_type_2(const EntityHandle *conn,
                            const BitRefEdges split_edges,
                            const EntityHandle *edge_new_nodes,
                            EntityHandle *new_prism_conn);
MoFEMErrorCode prism_type_3(const EntityHandle *conn,
                            const BitRefEdges split_edges,
                            const EntityHandle *edge_new_nodes,
                            EntityHandle *new_prism_conn);

// Quad
MoFEMErrorCode quad_split_all_edges(const EntityHandle *conn,
                                    const EntityHandle *edge_new_nodes,
                                    EntityHandle *new_quad_conn);

} // namespace MoFEM

#endif //__ENTITYREFINE_HPP__
