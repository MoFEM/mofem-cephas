/* \file EntityRefine.cpp
 * \brief Tetrahedral refinement algorithm

 It is based on \cite ruprecht1998scheme

 \todo tet refinement should be rewritten, with better error control, and more
 clear refinement patterns.

*/

/* This file is part of MoFEM.
 * MoFEM is free software: you can redistribute it and/or modify it under
 * the terms of the GNU Lesser General Public License as published by the
 * Free Software Foundation, either version 3 of the License, or (at your
 * option) any later version.
 *
 * MoFEM is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
 * License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with MoFEM. If not, see <http://www.gnu.org/licenses/>. */

// A scheme for Edge-based Adaptive Tetrahedral Subdivision; Delft Ruprecht

namespace MoFEM {

// TET
static constexpr int edges_conn[] = {0, 1, 1, 2, 2, 0, 0, 3, 1, 3, 2, 3};
static constexpr int oposite_edge[] = {5, 3, 4, 1, 2, 0};
static constexpr int edge_permutations[6][6] = {
    {0, 1, 2, 3, 4, 5}, {1, 2, 0, 4, 5, 3}, {2, 0, 1, 5, 3, 4},
    {3, 4, 0, 2, 5, 1}, {4, 5, 1, 0, 3, 2}, {5, 3, 2, 1, 4, 0}};
static constexpr int edge_mirror_cross[6] = {0, 3, 4, 1, 2, 5};
static constexpr int edge_mirror_vertical[6] = {0, 4, 3, 2, 1, 5};
static constexpr int cyclic_node_rotate_face_3[3][4] = {
    {3, 1, 0, 2}, {0, 1, 2, 3}, {2, 1, 3, 0}}; // 2,0,1
static constexpr int cyclic_edge_rotate_face_3[3][6] = {
    {4, 0, 3, 5, 1, 2}, {0, 1, 2, 3, 4, 5}, {1, 4, 5, 2, 0, 3}};
static constexpr char edge_bits_mark[] = {1, 2, 4, 8, 16, 32};

void tet_type_6(moab::Interface &moab, const EntityHandle *conn,
                const EntityHandle *edge_new_nodes,
                EntityHandle *new_tets_conn) {
  // 0:01 - 4 1:12-5 2:20-6 3:03-7 4:13-8 5:23-9
  // TET0
  new_tets_conn[0 * 4 + 0] = conn[0];
  new_tets_conn[0 * 4 + 1] = edge_new_nodes[0];
  new_tets_conn[0 * 4 + 2] = edge_new_nodes[2];
  new_tets_conn[0 * 4 + 3] = edge_new_nodes[3];
  // TET1
  new_tets_conn[1 * 4 + 0] = conn[1];
  new_tets_conn[1 * 4 + 1] = edge_new_nodes[0];
  new_tets_conn[1 * 4 + 2] = edge_new_nodes[4];
  new_tets_conn[1 * 4 + 3] = edge_new_nodes[1];
  // TET2
  new_tets_conn[2 * 4 + 0] = conn[2];
  new_tets_conn[2 * 4 + 1] = edge_new_nodes[1];
  new_tets_conn[2 * 4 + 2] = edge_new_nodes[5];
  new_tets_conn[2 * 4 + 3] = edge_new_nodes[2];
  // TET3
  new_tets_conn[3 * 4 + 0] = conn[3];
  new_tets_conn[3 * 4 + 1] = edge_new_nodes[3];
  new_tets_conn[3 * 4 + 2] = edge_new_nodes[5];
  new_tets_conn[3 * 4 + 3] = edge_new_nodes[4];
  double coords[6 * 3];
  moab.get_coords(edge_new_nodes, 6, coords);
  cblas_daxpy(3, -1, &coords[4 * 3], 1, &coords[2 * 3], 1);
  cblas_daxpy(3, -1, &coords[3 * 3], 1, &coords[1 * 3], 1);
  cblas_daxpy(3, -1, &coords[5 * 3], 1, &coords[0 * 3], 1);
  double L[3] = {cblas_dnrm2(3, &coords[2 * 3], 1),
                 cblas_dnrm2(3, &coords[1 * 3], 1),
                 cblas_dnrm2(3, &coords[0 * 3], 1)};
  // VARIANT 1 - diag 4-2
  if (L[0] <= L[1] && L[0] <= L[2]) {
    // TET4
    new_tets_conn[4 * 4 + 0] = edge_new_nodes[4];
    new_tets_conn[4 * 4 + 1] = edge_new_nodes[3];
    new_tets_conn[4 * 4 + 2] = edge_new_nodes[2];
    new_tets_conn[4 * 4 + 3] = edge_new_nodes[0];
    // TET5
    new_tets_conn[5 * 4 + 0] = edge_new_nodes[4];
    new_tets_conn[5 * 4 + 1] = edge_new_nodes[2];
    new_tets_conn[5 * 4 + 2] = edge_new_nodes[3];
    new_tets_conn[5 * 4 + 3] = edge_new_nodes[5];
    // TET6
    new_tets_conn[6 * 4 + 0] = edge_new_nodes[4];
    new_tets_conn[6 * 4 + 1] = edge_new_nodes[2];
    new_tets_conn[6 * 4 + 2] = edge_new_nodes[1];
    new_tets_conn[6 * 4 + 3] = edge_new_nodes[0];
    // TET7
    new_tets_conn[7 * 4 + 0] = edge_new_nodes[4];
    new_tets_conn[7 * 4 + 1] = edge_new_nodes[1];
    new_tets_conn[7 * 4 + 2] = edge_new_nodes[2];
    new_tets_conn[7 * 4 + 3] = edge_new_nodes[5];
    return;
  }
  // VARIANT 2 - diag 3-1
  if (L[1] <= L[0] && L[1] <= L[2]) {
    // TET4
    new_tets_conn[4 * 4 + 0] = edge_new_nodes[4];
    new_tets_conn[4 * 4 + 1] = edge_new_nodes[3];
    new_tets_conn[4 * 4 + 2] = edge_new_nodes[1];
    new_tets_conn[4 * 4 + 3] = edge_new_nodes[0];
    // TET5
    new_tets_conn[5 * 4 + 0] = edge_new_nodes[4];
    new_tets_conn[5 * 4 + 1] = edge_new_nodes[1];
    new_tets_conn[5 * 4 + 2] = edge_new_nodes[3];
    new_tets_conn[5 * 4 + 3] = edge_new_nodes[5];
    // TET6
    new_tets_conn[6 * 4 + 0] = edge_new_nodes[1];
    new_tets_conn[6 * 4 + 1] = edge_new_nodes[3];
    new_tets_conn[6 * 4 + 2] = edge_new_nodes[2];
    new_tets_conn[6 * 4 + 3] = edge_new_nodes[0];
    // TET7
    new_tets_conn[7 * 4 + 0] = edge_new_nodes[1];
    new_tets_conn[7 * 4 + 1] = edge_new_nodes[2];
    new_tets_conn[7 * 4 + 2] = edge_new_nodes[3];
    new_tets_conn[7 * 4 + 3] = edge_new_nodes[5];
    return;
  }
  // VARIANT 3 - diag 5-0
  // TET4
  new_tets_conn[4 * 4 + 0] = edge_new_nodes[5];
  new_tets_conn[4 * 4 + 1] = edge_new_nodes[2];
  new_tets_conn[4 * 4 + 2] = edge_new_nodes[0];
  new_tets_conn[4 * 4 + 3] = edge_new_nodes[3];
  // TET5
  new_tets_conn[5 * 4 + 0] = edge_new_nodes[5];
  new_tets_conn[5 * 4 + 1] = edge_new_nodes[0];
  new_tets_conn[5 * 4 + 2] = edge_new_nodes[2];
  new_tets_conn[5 * 4 + 3] = edge_new_nodes[1];
  // TET6
  new_tets_conn[6 * 4 + 0] = edge_new_nodes[5];
  new_tets_conn[6 * 4 + 1] = edge_new_nodes[0];
  new_tets_conn[6 * 4 + 2] = edge_new_nodes[4];
  new_tets_conn[6 * 4 + 3] = edge_new_nodes[3];
  // TET7
  new_tets_conn[7 * 4 + 0] = edge_new_nodes[5];
  new_tets_conn[7 * 4 + 1] = edge_new_nodes[4];
  new_tets_conn[7 * 4 + 2] = edge_new_nodes[0];
  new_tets_conn[7 * 4 + 3] = edge_new_nodes[1];
}
int tet_type_5(moab::Interface &moab, const EntityHandle *conn,
               const EntityHandle *edge_new_nodes,
               EntityHandle *new_tets_conn) {
  int free_edge = -1;
  for (int ee = 0; ee < 6; ee++) {
    if (edge_new_nodes[ee] == no_handle) {
      free_edge = ee;
      break;
    }
  }
  int edge0 = oposite_edge[free_edge];
  EntityHandle conn_[] = {
      conn[edges_conn[edge0 * 2 + 0]], conn[edges_conn[edge0 * 2 + 1]],
      conn[edges_conn[free_edge * 2 + 0]], conn[edges_conn[free_edge * 2 + 1]]};
  const int *edges_ = &edge_permutations[edge0][0];
  EntityHandle edge_new_nodes_[6];
  for (int ee = 0; ee < 6; ee++)
    edge_new_nodes_[ee] = edge_new_nodes[edges_[ee]];
  bool free_edge_swappped = false;
  if (conn_[3] < conn_[2]) {
    free_edge_swappped = true;
    EntityHandle conn__2_ = conn_[2];
    conn_[2] = conn_[3];
    conn_[3] = conn__2_;
  }
  assert(conn_[0] != no_handle);
  assert(conn_[1] != no_handle);
  assert(conn_[2] != no_handle);
  assert(conn_[3] != no_handle);
  assert(edge_new_nodes_[0] != no_handle);
  assert(edge_new_nodes_[1] != no_handle);
  assert(edge_new_nodes_[2] != no_handle);
  assert(edge_new_nodes_[3] != no_handle);
  assert(edge_new_nodes_[4] != no_handle);
  // TET0
  new_tets_conn[0 * 4 + 0] = edge_new_nodes_[4];
  new_tets_conn[0 * 4 + 1] = edge_new_nodes_[0];
  new_tets_conn[0 * 4 + 2] = edge_new_nodes_[1];
  new_tets_conn[0 * 4 + 3] = conn_[1];
  // TET1
  new_tets_conn[1 * 4 + 0] = edge_new_nodes_[2];
  new_tets_conn[1 * 4 + 1] = edge_new_nodes_[0];
  new_tets_conn[1 * 4 + 2] = edge_new_nodes_[3];
  new_tets_conn[1 * 4 + 3] = conn_[0];
  // TET4
  new_tets_conn[2 * 4 + 0] = conn_[2];
  new_tets_conn[2 * 4 + 1] = edge_new_nodes_[3];
  new_tets_conn[2 * 4 + 2] = edge_new_nodes_[4];
  if (free_edge_swappped) {
    new_tets_conn[2 * 4 + 1] = edge_new_nodes_[1];
    new_tets_conn[2 * 4 + 2] = edge_new_nodes_[2];
  }
  new_tets_conn[2 * 4 + 3] = conn_[3];
  double coords[6 * 3];
  moab.get_coords(edge_new_nodes_, 6, coords);
  cblas_daxpy(3, -1, &coords[4 * 3], 1, &coords[2 * 3], 1);
  cblas_daxpy(3, -1, &coords[3 * 3], 1, &coords[1 * 3], 1);
  double L[2] = {cblas_dnrm2(3, &coords[2 * 3], 1),
                 cblas_dnrm2(3, &coords[1 * 3], 1)};
  if (L[1] <= L[0]) {
    // VARIANT 1 diag 4-2
    // TET2
    new_tets_conn[3 * 4 + 0] = edge_new_nodes_[4];
    new_tets_conn[3 * 4 + 1] = edge_new_nodes_[3];
    new_tets_conn[3 * 4 + 2] = edge_new_nodes_[1];
    new_tets_conn[3 * 4 + 3] = edge_new_nodes_[0];
    // TET3
    new_tets_conn[4 * 4 + 0] = edge_new_nodes_[2];
    new_tets_conn[4 * 4 + 1] = edge_new_nodes_[1];
    new_tets_conn[4 * 4 + 2] = edge_new_nodes_[3];
    new_tets_conn[4 * 4 + 3] = edge_new_nodes_[0];
    // TET5
    new_tets_conn[5 * 4 + 0] = edge_new_nodes_[4];
    new_tets_conn[5 * 4 + 1] = edge_new_nodes_[1];
    new_tets_conn[5 * 4 + 2] = edge_new_nodes_[3];
    new_tets_conn[5 * 4 + 3] = conn_[2];
    // TET6
    new_tets_conn[6 * 4 + 0] = edge_new_nodes_[1];
    new_tets_conn[6 * 4 + 1] = edge_new_nodes_[2];
    new_tets_conn[6 * 4 + 2] = edge_new_nodes_[3];
    new_tets_conn[6 * 4 + 3] = conn_[2];
    return 1;
  }
  // VARIANT 2 diag 1-3
  // TET2
  new_tets_conn[3 * 4 + 0] = edge_new_nodes_[4];
  new_tets_conn[3 * 4 + 1] = edge_new_nodes_[3];
  new_tets_conn[3 * 4 + 2] = edge_new_nodes_[2];
  new_tets_conn[3 * 4 + 3] = edge_new_nodes_[0];
  // TET3
  new_tets_conn[4 * 4 + 0] = edge_new_nodes_[2];
  new_tets_conn[4 * 4 + 1] = edge_new_nodes_[1];
  new_tets_conn[4 * 4 + 2] = edge_new_nodes_[4];
  new_tets_conn[4 * 4 + 3] = edge_new_nodes_[0];
  // TET5
  new_tets_conn[5 * 4 + 0] = edge_new_nodes_[4];
  new_tets_conn[5 * 4 + 1] = edge_new_nodes_[1];
  new_tets_conn[5 * 4 + 2] = edge_new_nodes_[2];
  new_tets_conn[5 * 4 + 3] = conn_[2];
  // TET6
  new_tets_conn[6 * 4 + 0] = edge_new_nodes_[4];
  new_tets_conn[6 * 4 + 1] = edge_new_nodes_[2];
  new_tets_conn[6 * 4 + 2] = edge_new_nodes_[3];
  new_tets_conn[6 * 4 + 3] = conn_[2];
  return 0;
}
int tet_type_4(const EntityHandle *conn, const int *split_edges,
               const EntityHandle *edge_new_nodes,
               EntityHandle *new_tets_conn) {
  char mach_pattern = 0;
  for (int ee = 0; ee < 4; ee++)
    mach_pattern |= edge_bits_mark[split_edges[ee]];
  EntityHandle conn_[4] = {no_handle, no_handle, no_handle, no_handle};
  int edges_[6];
  int type = -1;
  for (int ee = 0; ee < 6; ee++) {
    char pattern0, pattern1;
    pattern0 = edge_bits_mark[edge_permutations[ee][0]] |
               edge_bits_mark[edge_permutations[ee][1]] |
               edge_bits_mark[edge_permutations[ee][4]] |
               edge_bits_mark[edge_permutations[ee][3]];
    pattern1 = edge_bits_mark[edge_permutations[ee][1]] |
               edge_bits_mark[edge_permutations[ee][4]] |
               edge_bits_mark[edge_permutations[ee][2]] |
               edge_bits_mark[edge_permutations[ee][3]];
    if (pattern0 == mach_pattern || pattern1 == mach_pattern) {
      int free_edge = oposite_edge[ee];
      conn_[0] = conn[edges_conn[ee * 2 + 0]];
      conn_[1] = conn[edges_conn[ee * 2 + 1]];
      conn_[2] = conn[edges_conn[free_edge * 2 + 0]];
      conn_[3] = conn[edges_conn[free_edge * 2 + 1]];
      for (int EE = 0; EE < 6; EE++)
        edges_[EE] = edge_permutations[ee][EE];
      if (pattern0 == mach_pattern)
        type = 0;
      else if (pattern1 == mach_pattern)
        type = 1;
      // printf("no mirror\n");
      break;
    }
    pattern0 = edge_bits_mark[edge_permutations[ee][edge_mirror_cross[0]]] |
               edge_bits_mark[edge_permutations[ee][edge_mirror_cross[1]]] |
               edge_bits_mark[edge_permutations[ee][edge_mirror_cross[4]]] |
               edge_bits_mark[edge_permutations[ee][edge_mirror_cross[3]]];
    pattern1 = edge_bits_mark[edge_permutations[ee][1]] |
               edge_bits_mark[edge_permutations[ee][4]] |
               edge_bits_mark[edge_permutations[ee][2]] |
               edge_bits_mark[edge_permutations[ee][3]];
    if (pattern0 == mach_pattern || pattern1 == mach_pattern) {
      int free_edge = oposite_edge[ee];
      conn_[0] = conn[edges_conn[ee * 2 + 1]];
      conn_[1] = conn[edges_conn[ee * 2 + 0]];
      conn_[2] = conn[edges_conn[free_edge * 2 + 1]];
      conn_[3] = conn[edges_conn[free_edge * 2 + 0]];
      for (int EE = 0; EE < 6; EE++)
        edges_[EE] = edge_permutations[ee][edge_mirror_cross[EE]];
      if (pattern0 == mach_pattern)
        type = 0;
      else if (pattern1 == mach_pattern)
        type = 1;
      // printf("mirror\n");
      break;
    }
  }
  assert(type != -1);
  EntityHandle edge_new_nodes_[6];
  for (int ee = 0; ee < 6; ee++)
    edge_new_nodes_[ee] = edge_new_nodes[edges_[ee]];
  if (type == 0) {
    assert(edge_new_nodes_[0] != no_handle);
    assert(edge_new_nodes_[1] != no_handle);
    assert(edge_new_nodes_[4] != no_handle);
    assert(edge_new_nodes_[3] != no_handle);
    bool free_edge_swappped5 = false;
    if (conn_[3] < conn_[2]) {
      free_edge_swappped5 = true;
    }
    bool free_edge_swappped2 = false;
    if (conn_[0] < conn_[2]) {
      free_edge_swappped2 = true;
    }
    // TET0
    new_tets_conn[0 * 4 + 0] = conn_[1];
    new_tets_conn[0 * 4 + 1] = edge_new_nodes_[1];
    new_tets_conn[0 * 4 + 2] = edge_new_nodes_[0];
    new_tets_conn[0 * 4 + 3] = edge_new_nodes_[4];
    if (free_edge_swappped5 && (!free_edge_swappped2)) {
      // TET1
      new_tets_conn[1 * 4 + 0] = edge_new_nodes_[1];
      new_tets_conn[1 * 4 + 1] = edge_new_nodes_[3];
      new_tets_conn[1 * 4 + 2] = edge_new_nodes_[4];
      new_tets_conn[1 * 4 + 3] = conn_[3];
      // TET2
      new_tets_conn[2 * 4 + 0] = conn_[3];
      new_tets_conn[2 * 4 + 1] = edge_new_nodes_[1];
      new_tets_conn[2 * 4 + 2] = edge_new_nodes_[3];
      new_tets_conn[2 * 4 + 3] = conn_[2];
      // TET3
      new_tets_conn[3 * 4 + 0] = edge_new_nodes_[0];
      new_tets_conn[3 * 4 + 1] = edge_new_nodes_[3];
      new_tets_conn[3 * 4 + 2] = edge_new_nodes_[4];
      new_tets_conn[3 * 4 + 3] = edge_new_nodes_[1];
      // TET4
      new_tets_conn[4 * 4 + 0] = conn_[2];
      new_tets_conn[4 * 4 + 1] = edge_new_nodes_[0];
      new_tets_conn[4 * 4 + 2] = edge_new_nodes_[3];
      new_tets_conn[4 * 4 + 3] = conn_[0];
      // TET5
      new_tets_conn[5 * 4 + 0] = edge_new_nodes_[1];
      new_tets_conn[5 * 4 + 1] = edge_new_nodes_[0];
      new_tets_conn[5 * 4 + 2] = edge_new_nodes_[3];
      new_tets_conn[5 * 4 + 3] = conn_[2];
      return 2;
    } else if (free_edge_swappped2 && (!free_edge_swappped5)) {
      // TET1
      new_tets_conn[1 * 4 + 0] = edge_new_nodes_[3];
      new_tets_conn[1 * 4 + 1] = edge_new_nodes_[1];
      new_tets_conn[1 * 4 + 2] = edge_new_nodes_[0];
      new_tets_conn[1 * 4 + 3] = conn_[0];
      // TET2
      new_tets_conn[2 * 4 + 0] = edge_new_nodes_[3];
      new_tets_conn[2 * 4 + 1] = edge_new_nodes_[1];
      new_tets_conn[2 * 4 + 2] = conn_[0];
      new_tets_conn[2 * 4 + 3] = conn_[2];
      // TET3
      new_tets_conn[3 * 4 + 0] = edge_new_nodes_[0];
      new_tets_conn[3 * 4 + 1] = edge_new_nodes_[3];
      new_tets_conn[3 * 4 + 2] = edge_new_nodes_[4];
      new_tets_conn[3 * 4 + 3] = edge_new_nodes_[1];
      // TET4
      new_tets_conn[4 * 4 + 0] = conn_[2];
      new_tets_conn[4 * 4 + 1] = edge_new_nodes_[3];
      new_tets_conn[4 * 4 + 2] = edge_new_nodes_[4];
      new_tets_conn[4 * 4 + 3] = conn_[3];
      // TET5
      new_tets_conn[5 * 4 + 0] = edge_new_nodes_[1];
      new_tets_conn[5 * 4 + 1] = edge_new_nodes_[3];
      new_tets_conn[5 * 4 + 2] = edge_new_nodes_[4];
      new_tets_conn[5 * 4 + 3] = conn_[2];
      return 3;
    } else if (free_edge_swappped2 && free_edge_swappped5) {
      // TET1
      new_tets_conn[1 * 4 + 0] = edge_new_nodes_[3];
      new_tets_conn[1 * 4 + 1] = edge_new_nodes_[1];
      new_tets_conn[1 * 4 + 2] = edge_new_nodes_[0];
      new_tets_conn[1 * 4 + 3] = conn_[0];
      // TET2
      new_tets_conn[2 * 4 + 0] = edge_new_nodes_[3];
      new_tets_conn[2 * 4 + 1] = edge_new_nodes_[1];
      new_tets_conn[2 * 4 + 2] = conn_[0];
      new_tets_conn[2 * 4 + 3] = conn_[2];
      // TET3
      new_tets_conn[3 * 4 + 0] = edge_new_nodes_[0];
      new_tets_conn[3 * 4 + 1] = edge_new_nodes_[3];
      new_tets_conn[3 * 4 + 2] = edge_new_nodes_[4];
      new_tets_conn[3 * 4 + 3] = edge_new_nodes_[1];
      // TET4
      new_tets_conn[4 * 4 + 0] = edge_new_nodes_[1];
      new_tets_conn[4 * 4 + 1] = edge_new_nodes_[3];
      new_tets_conn[4 * 4 + 2] = edge_new_nodes_[4];
      new_tets_conn[4 * 4 + 3] = conn_[3];
      // TET5
      new_tets_conn[5 * 4 + 0] = conn_[3];
      new_tets_conn[5 * 4 + 1] = edge_new_nodes_[1];
      new_tets_conn[5 * 4 + 2] = edge_new_nodes_[3];
      new_tets_conn[5 * 4 + 3] = conn_[2];
      return 4;
    } else {
      // TET1
      new_tets_conn[1 * 4 + 0] = edge_new_nodes_[0];
      new_tets_conn[1 * 4 + 1] = edge_new_nodes_[3];
      new_tets_conn[1 * 4 + 2] = edge_new_nodes_[4];
      new_tets_conn[1 * 4 + 3] = conn_[2];
      // TET2
      new_tets_conn[2 * 4 + 0] = conn_[2];
      new_tets_conn[2 * 4 + 1] = edge_new_nodes_[1];
      new_tets_conn[2 * 4 + 2] = edge_new_nodes_[4];
      new_tets_conn[2 * 4 + 3] = edge_new_nodes_[0];
      // TET3
      new_tets_conn[3 * 4 + 0] = conn_[2];
      new_tets_conn[3 * 4 + 1] = edge_new_nodes_[0];
      new_tets_conn[3 * 4 + 2] = edge_new_nodes_[3];
      new_tets_conn[3 * 4 + 3] = conn_[0];
      // TET4
      new_tets_conn[4 * 4 + 0] = conn_[2];
      new_tets_conn[4 * 4 + 1] = edge_new_nodes_[3];
      new_tets_conn[4 * 4 + 2] = edge_new_nodes_[4];
      new_tets_conn[4 * 4 + 3] = conn_[3];
    }
  } else if (type == 1) {
    assert(edge_new_nodes_[1] != no_handle);
    assert(edge_new_nodes_[2] != no_handle);
    assert(edge_new_nodes_[3] != no_handle);
    assert(edge_new_nodes_[4] != no_handle);
    bool free_edge_swappped5 = false;
    if (conn_[3] < conn_[2]) {
      free_edge_swappped5 = true;
    }
    bool free_edge_swappped0 = false;
    if (conn_[0] > conn_[1]) {
      free_edge_swappped0 = true;
    }
    if (free_edge_swappped0 && (!free_edge_swappped5)) {
      // TET0
      new_tets_conn[0 * 4 + 0] = edge_new_nodes_[3];
      new_tets_conn[0 * 4 + 1] = edge_new_nodes_[2];
      new_tets_conn[0 * 4 + 2] = conn_[1];
      new_tets_conn[0 * 4 + 3] = conn_[0];
      // TET1
      new_tets_conn[1 * 4 + 0] = edge_new_nodes_[3];
      new_tets_conn[1 * 4 + 1] = edge_new_nodes_[1];
      new_tets_conn[1 * 4 + 2] = edge_new_nodes_[4];
      new_tets_conn[1 * 4 + 3] = conn_[1];
      // TET2
      new_tets_conn[2 * 4 + 0] = edge_new_nodes_[3];
      new_tets_conn[2 * 4 + 1] = edge_new_nodes_[2];
      new_tets_conn[2 * 4 + 2] = edge_new_nodes_[1];
      new_tets_conn[2 * 4 + 3] = conn_[1];
      // TET3
      new_tets_conn[3 * 4 + 0] = edge_new_nodes_[1];
      new_tets_conn[3 * 4 + 1] = edge_new_nodes_[2];
      new_tets_conn[3 * 4 + 2] = edge_new_nodes_[3];
      new_tets_conn[3 * 4 + 3] = conn_[2];
      // TET4
      new_tets_conn[4 * 4 + 0] = edge_new_nodes_[1];
      new_tets_conn[4 * 4 + 1] = edge_new_nodes_[3];
      new_tets_conn[4 * 4 + 2] = edge_new_nodes_[4];
      new_tets_conn[4 * 4 + 3] = conn_[2];
      // TET5
      new_tets_conn[5 * 4 + 0] = edge_new_nodes_[3];
      new_tets_conn[5 * 4 + 1] = edge_new_nodes_[4];
      new_tets_conn[5 * 4 + 2] = conn_[2];
      new_tets_conn[5 * 4 + 3] = conn_[3];
      return 5;
    } else if (free_edge_swappped5 && (!free_edge_swappped0)) {
      // TET0
      new_tets_conn[0 * 4 + 0] = edge_new_nodes_[4];
      new_tets_conn[0 * 4 + 1] = edge_new_nodes_[1];
      new_tets_conn[0 * 4 + 2] = conn_[1];
      new_tets_conn[0 * 4 + 3] = conn_[0];
      // TET1
      new_tets_conn[1 * 4 + 0] = edge_new_nodes_[3];
      new_tets_conn[1 * 4 + 1] = edge_new_nodes_[1];
      new_tets_conn[1 * 4 + 2] = edge_new_nodes_[4];
      new_tets_conn[1 * 4 + 3] = conn_[0];
      // TET2
      new_tets_conn[2 * 4 + 0] = edge_new_nodes_[2];
      new_tets_conn[2 * 4 + 1] = edge_new_nodes_[1];
      new_tets_conn[2 * 4 + 2] = edge_new_nodes_[3];
      new_tets_conn[2 * 4 + 3] = conn_[0];
      // TET3
      new_tets_conn[3 * 4 + 0] = edge_new_nodes_[1];
      new_tets_conn[3 * 4 + 1] = edge_new_nodes_[2];
      new_tets_conn[3 * 4 + 2] = edge_new_nodes_[3];
      new_tets_conn[3 * 4 + 3] = conn_[3];
      // TET4
      new_tets_conn[4 * 4 + 0] = edge_new_nodes_[1];
      new_tets_conn[4 * 4 + 1] = edge_new_nodes_[3];
      new_tets_conn[4 * 4 + 2] = edge_new_nodes_[4];
      new_tets_conn[4 * 4 + 3] = conn_[3];
      // TET5
      new_tets_conn[5 * 4 + 0] = edge_new_nodes_[1];
      new_tets_conn[5 * 4 + 1] = edge_new_nodes_[2];
      new_tets_conn[5 * 4 + 2] = conn_[3];
      new_tets_conn[5 * 4 + 3] = conn_[2];
      return 6;
    } else if (free_edge_swappped5 && free_edge_swappped0) {
      // TET0
      new_tets_conn[0 * 4 + 0] = edge_new_nodes_[3];
      new_tets_conn[0 * 4 + 1] = edge_new_nodes_[2];
      new_tets_conn[0 * 4 + 2] = conn_[1];
      new_tets_conn[0 * 4 + 3] = conn_[0];
      // TET1
      new_tets_conn[1 * 4 + 0] = edge_new_nodes_[3];
      new_tets_conn[1 * 4 + 1] = edge_new_nodes_[1];
      new_tets_conn[1 * 4 + 2] = edge_new_nodes_[4];
      new_tets_conn[1 * 4 + 3] = conn_[1];
      // TET2
      new_tets_conn[2 * 4 + 0] = edge_new_nodes_[3];
      new_tets_conn[2 * 4 + 1] = edge_new_nodes_[2];
      new_tets_conn[2 * 4 + 2] = edge_new_nodes_[1];
      new_tets_conn[2 * 4 + 3] = conn_[1];
      // TET3
      new_tets_conn[3 * 4 + 0] = edge_new_nodes_[1];
      new_tets_conn[3 * 4 + 1] = edge_new_nodes_[2];
      new_tets_conn[3 * 4 + 2] = edge_new_nodes_[3];
      new_tets_conn[3 * 4 + 3] = conn_[3];
      // TET4
      new_tets_conn[4 * 4 + 0] = edge_new_nodes_[1];
      new_tets_conn[4 * 4 + 1] = edge_new_nodes_[3];
      new_tets_conn[4 * 4 + 2] = edge_new_nodes_[4];
      new_tets_conn[4 * 4 + 3] = conn_[3];
      // TET5
      new_tets_conn[5 * 4 + 0] = edge_new_nodes_[1];
      new_tets_conn[5 * 4 + 1] = edge_new_nodes_[2];
      new_tets_conn[5 * 4 + 2] = conn_[3];
      new_tets_conn[5 * 4 + 3] = conn_[2];
      return 7;
    }
    // TET0
    new_tets_conn[0 * 4 + 0] = edge_new_nodes_[4];
    new_tets_conn[0 * 4 + 1] = edge_new_nodes_[1];
    new_tets_conn[0 * 4 + 2] = conn_[1];
    new_tets_conn[0 * 4 + 3] = conn_[0];
    // TET1
    new_tets_conn[1 * 4 + 0] = edge_new_nodes_[3];
    new_tets_conn[1 * 4 + 1] = edge_new_nodes_[1];
    new_tets_conn[1 * 4 + 2] = edge_new_nodes_[4];
    new_tets_conn[1 * 4 + 3] = conn_[0];
    // TET2
    new_tets_conn[2 * 4 + 0] = edge_new_nodes_[2];
    new_tets_conn[2 * 4 + 1] = edge_new_nodes_[1];
    new_tets_conn[2 * 4 + 2] = edge_new_nodes_[3];
    new_tets_conn[2 * 4 + 3] = conn_[0];
    // TET3
    new_tets_conn[3 * 4 + 0] = edge_new_nodes_[1];
    new_tets_conn[3 * 4 + 1] = edge_new_nodes_[2];
    new_tets_conn[3 * 4 + 2] = edge_new_nodes_[3];
    new_tets_conn[3 * 4 + 3] = conn_[2];
    // TET4
    new_tets_conn[4 * 4 + 0] = edge_new_nodes_[1];
    new_tets_conn[4 * 4 + 1] = edge_new_nodes_[3];
    new_tets_conn[4 * 4 + 2] = edge_new_nodes_[4];
    new_tets_conn[4 * 4 + 3] = conn_[2];
    // TET5
    new_tets_conn[5 * 4 + 0] = edge_new_nodes_[3];
    new_tets_conn[5 * 4 + 1] = edge_new_nodes_[4];
    new_tets_conn[5 * 4 + 2] = conn_[2];
    new_tets_conn[5 * 4 + 3] = conn_[3];
  }
  return type;
}
int tet_type_3(const EntityHandle *conn, const int *split_edges,
               const EntityHandle *edge_new_nodes,
               EntityHandle *new_tets_conn) {
  char mach_pattern = 0;
  for (int ee = 0; ee < 3; ee++)
    mach_pattern |= edge_bits_mark[split_edges[ee]];
  EntityHandle conn_[4];
  int edges_[6];
  int type = -1;
  bool is_rotated = false;
  for (int ee = 0; ee < 6; ee++) {
    char pattern0, pattern1, pattern2;
    pattern0 = edge_bits_mark[edge_permutations[ee][0]] |
               edge_bits_mark[edge_permutations[ee][1]] |
               edge_bits_mark[edge_permutations[ee][4]];
    pattern1 = edge_bits_mark[edge_permutations[ee][1]] |
               edge_bits_mark[edge_permutations[ee][4]] |
               edge_bits_mark[edge_permutations[ee][5]];
    pattern2 = edge_bits_mark[edge_permutations[ee][4]] |
               edge_bits_mark[edge_permutations[ee][1]] |
               edge_bits_mark[edge_permutations[ee][2]];
    if (pattern0 == mach_pattern || pattern1 == mach_pattern ||
        pattern2 == mach_pattern) {
      // printf("nothing\n");
      int free_edge = oposite_edge[ee];
      conn_[0] = conn[edges_conn[ee * 2 + 0]];
      conn_[1] = conn[edges_conn[ee * 2 + 1]];
      conn_[2] = conn[edges_conn[free_edge * 2 + 0]];
      conn_[3] = conn[edges_conn[free_edge * 2 + 1]];
      for (int EE = 0; EE < 6; EE++)
        edges_[EE] = edge_permutations[ee][EE];
      if (pattern0 == mach_pattern) {
        // printf("nothing\n");
        type = 0;
      } else if (pattern1 == mach_pattern)
        type = 1;
      else if (pattern2 == mach_pattern)
        type = 2;
      break;
    }
    pattern0 = edge_bits_mark[edge_permutations[ee][edge_mirror_cross[0]]] |
               edge_bits_mark[edge_permutations[ee][edge_mirror_cross[1]]] |
               edge_bits_mark[edge_permutations[ee][edge_mirror_cross[4]]];
    pattern1 = edge_bits_mark[edge_permutations[ee][edge_mirror_cross[1]]] |
               edge_bits_mark[edge_permutations[ee][edge_mirror_cross[4]]] |
               edge_bits_mark[edge_permutations[ee][edge_mirror_cross[5]]];
    pattern2 = edge_bits_mark[edge_permutations[ee][edge_mirror_cross[4]]] |
               edge_bits_mark[edge_permutations[ee][edge_mirror_cross[1]]] |
               edge_bits_mark[edge_permutations[ee][edge_mirror_cross[2]]];
    if (pattern0 == mach_pattern || pattern1 == mach_pattern ||
        pattern2 == mach_pattern) {
      // printf("edge_mirror_cross\n");
      int free_edge = oposite_edge[ee];
      conn_[0] = conn[edges_conn[ee * 2 + 1]];
      conn_[1] = conn[edges_conn[ee * 2 + 0]];
      conn_[2] = conn[edges_conn[free_edge * 2 + 1]];
      conn_[3] = conn[edges_conn[free_edge * 2 + 0]];
      for (int EE = 0; EE < 6; EE++)
        edges_[EE] = edge_permutations[ee][edge_mirror_cross[EE]];
      if (pattern0 == mach_pattern) {
        // printf("edge_mirror_cross\n");
        type = 0;
      } else if (pattern1 == mach_pattern)
        type = 1;
      else if (pattern2 == mach_pattern)
        type = 2;
      break;
    }
    pattern2 = edge_bits_mark[edge_permutations[ee][edge_mirror_vertical[4]]] |
               edge_bits_mark[edge_permutations[ee][edge_mirror_vertical[1]]] |
               edge_bits_mark[edge_permutations[ee][edge_mirror_vertical[2]]];
    if (pattern2 == mach_pattern) {
      is_rotated = true;
      // printf("edge_mirror_vertical\n");
      int free_edge = oposite_edge[ee];
      conn_[0] = conn[edges_conn[ee * 2 + 0]];
      conn_[1] = conn[edges_conn[ee * 2 + 1]];
      conn_[2] = conn[edges_conn[free_edge * 2 + 1]];
      conn_[3] = conn[edges_conn[free_edge * 2 + 0]];
      for (int EE = 0; EE < 6; EE++)
        edges_[EE] = edge_permutations[ee][edge_mirror_vertical[EE]];
      if (pattern0 == mach_pattern) {
        // printf("edge_mirror_vertical\n");
        type = 0;
      } else if (pattern1 == mach_pattern)
        type = 1;
      else if (pattern2 == mach_pattern)
        type = 2;
      break;
    }
    pattern2 =
        edge_bits_mark[edge_permutations
                           [ee][edge_mirror_cross[edge_mirror_vertical[4]]]] |
        edge_bits_mark[edge_permutations
                           [ee][edge_mirror_cross[edge_mirror_vertical[1]]]] |
        edge_bits_mark
            [edge_permutations[ee][edge_mirror_cross[edge_mirror_vertical[2]]]];
    if (pattern2 == mach_pattern) {
      is_rotated = true;
      int free_edge = oposite_edge[ee];
      conn_[0] = conn[edges_conn[ee * 2 + 1]];
      conn_[1] = conn[edges_conn[ee * 2 + 0]];
      conn_[2] = conn[edges_conn[free_edge * 2 + 0]];
      conn_[3] = conn[edges_conn[free_edge * 2 + 1]];
      for (int EE = 0; EE < 6; EE++)
        edges_[EE] =
            edge_permutations[ee][edge_mirror_cross[edge_mirror_vertical[EE]]];
      if (pattern0 == mach_pattern) {
        // printf("edge_mirror_cross|edge_mirror_vertical\n");
        type = 0;
      } else if (pattern1 == mach_pattern)
        type = 1;
      else if (pattern2 == mach_pattern)
        type = 2;
      break;
    }
  }
  assert(type != -1);
  EntityHandle edge_new_nodes_[6];
  for (int ee = 0; ee < 6; ee++)
    edge_new_nodes_[ee] = edge_new_nodes[edges_[ee]];
  if (type == 0) {
    EntityHandle conn__[4];
    EntityHandle edge_new_nodes__[6];
    bcopy(conn_, conn__, 4 * sizeof(EntityHandle));
    bcopy(edge_new_nodes_, edge_new_nodes__, 6 * sizeof(EntityHandle));
    for (int rotate_idx = 0; rotate_idx < 3; rotate_idx++) {
      // fprintf(stderr,"%d\n",rotate_idx);
      for (int ii = 0; ii < 4; ii++)
        conn_[ii] = conn__[cyclic_node_rotate_face_3[rotate_idx][ii]];
      for (int ee = 0; ee < 6; ee++)
        edge_new_nodes_[ee] =
            edge_new_nodes__[cyclic_edge_rotate_face_3[rotate_idx][ee]];
      if ((conn_[0] > conn_[2]) && (conn_[0] > conn_[3]))
        break;
    }
    assert(conn_[0] > conn_[2]);
    assert(conn_[0] > conn_[3]);
    assert(edge_new_nodes_[0] != no_handle);
    assert(edge_new_nodes_[1] != no_handle);
    assert(edge_new_nodes_[4] != no_handle);
    // TET0
    new_tets_conn[0 * 4 + 0] = edge_new_nodes_[0];
    new_tets_conn[0 * 4 + 1] = edge_new_nodes_[1];
    new_tets_conn[0 * 4 + 2] = edge_new_nodes_[4];
    new_tets_conn[0 * 4 + 3] = conn_[1];
    bool free_edge_swappped5 = false;
    if (conn_[3] < conn_[2]) {
      free_edge_swappped5 = true;
    }
    if (free_edge_swappped5) {
      // TET1
      new_tets_conn[1 * 4 + 0] = edge_new_nodes_[1];
      new_tets_conn[1 * 4 + 1] = edge_new_nodes_[0];
      new_tets_conn[1 * 4 + 2] = edge_new_nodes_[4];
      new_tets_conn[1 * 4 + 3] = conn_[3];
      // TET2
      new_tets_conn[2 * 4 + 0] = edge_new_nodes_[0];
      new_tets_conn[2 * 4 + 1] = edge_new_nodes_[1];
      new_tets_conn[2 * 4 + 2] = conn_[2];
      new_tets_conn[2 * 4 + 3] = conn_[3];
      // TET3
      new_tets_conn[3 * 4 + 0] = edge_new_nodes_[0];
      new_tets_conn[3 * 4 + 1] = conn_[0];
      new_tets_conn[3 * 4 + 2] = conn_[3];
      new_tets_conn[3 * 4 + 3] = conn_[2];
      // printf("free_edge_swappped5\n");
      return 4;
    }
    assert(conn_[3] > conn_[2]);
    // TET1
    new_tets_conn[1 * 4 + 0] = edge_new_nodes_[1];
    new_tets_conn[1 * 4 + 1] = edge_new_nodes_[0];
    new_tets_conn[1 * 4 + 2] = edge_new_nodes_[4];
    new_tets_conn[1 * 4 + 3] = conn_[2];
    // TET2
    new_tets_conn[2 * 4 + 0] = edge_new_nodes_[0];
    new_tets_conn[2 * 4 + 1] = edge_new_nodes_[4];
    new_tets_conn[2 * 4 + 2] = conn_[2];
    new_tets_conn[2 * 4 + 3] = conn_[3];
    // TET3
    new_tets_conn[3 * 4 + 0] = edge_new_nodes_[0];
    new_tets_conn[3 * 4 + 1] = conn_[0];
    new_tets_conn[3 * 4 + 2] = conn_[3];
    new_tets_conn[3 * 4 + 3] = conn_[2];
    // printf("no free_edge_swappped5\n");
  } else if (type == 1) {
    assert(edge_new_nodes_[1] != no_handle);
    assert(edge_new_nodes_[4] != no_handle);
    assert(edge_new_nodes_[5] != no_handle);
    // TET0
    new_tets_conn[0 * 4 + 0] = edge_new_nodes_[1];
    new_tets_conn[0 * 4 + 1] = edge_new_nodes_[4];
    new_tets_conn[0 * 4 + 2] = edge_new_nodes_[5];
    new_tets_conn[0 * 4 + 3] = conn_[0];
    // TET1
    new_tets_conn[1 * 4 + 0] = edge_new_nodes_[4];
    new_tets_conn[1 * 4 + 1] = edge_new_nodes_[5];
    new_tets_conn[1 * 4 + 2] = conn_[0];
    new_tets_conn[1 * 4 + 3] = conn_[3];
    // TET2
    new_tets_conn[2 * 4 + 0] = edge_new_nodes_[1];
    new_tets_conn[2 * 4 + 1] = edge_new_nodes_[5];
    new_tets_conn[2 * 4 + 2] = conn_[2];
    new_tets_conn[2 * 4 + 3] = conn_[0];
    // TET3
    new_tets_conn[3 * 4 + 0] = edge_new_nodes_[4];
    new_tets_conn[3 * 4 + 1] = edge_new_nodes_[1];
    new_tets_conn[3 * 4 + 2] = conn_[1];
    new_tets_conn[3 * 4 + 3] = conn_[0];
  } else if (type == 2) {
    assert(edge_new_nodes_[1] != no_handle);
    assert(edge_new_nodes_[4] != no_handle);
    assert(edge_new_nodes_[2] != no_handle);
    bool free_edge_swappped5 = false;
    if (conn_[3] < conn_[2]) {
      free_edge_swappped5 = true;
    }
    bool free_edge_swappped0 = false;
    if (conn_[0] > conn_[1]) {
      free_edge_swappped0 = true;
    }
    if (free_edge_swappped5 && (!free_edge_swappped0)) {
      // TET0
      new_tets_conn[0 * 4 + 0] = edge_new_nodes_[4];
      new_tets_conn[0 * 4 + 1] = edge_new_nodes_[1];
      if (is_rotated) {
        new_tets_conn[0 * 4 + 2] = conn_[0];
        new_tets_conn[0 * 4 + 3] = conn_[1];
      } else {
        new_tets_conn[0 * 4 + 2] = conn_[1];
        new_tets_conn[0 * 4 + 3] = conn_[0];
      }
      // TET1
      new_tets_conn[1 * 4 + 0] = edge_new_nodes_[1];
      new_tets_conn[1 * 4 + 1] = edge_new_nodes_[2];
      if (is_rotated) {
        new_tets_conn[1 * 4 + 2] = conn_[3];
        new_tets_conn[1 * 4 + 3] = conn_[0];
      } else {
        new_tets_conn[1 * 4 + 2] = conn_[0];
        new_tets_conn[1 * 4 + 3] = conn_[3];
      }
      // TET2
      new_tets_conn[2 * 4 + 0] = edge_new_nodes_[1];
      new_tets_conn[2 * 4 + 1] = edge_new_nodes_[2];
      if (is_rotated) {
        new_tets_conn[2 * 4 + 2] = conn_[2];
        new_tets_conn[2 * 4 + 3] = conn_[3];
      } else {
        new_tets_conn[2 * 4 + 2] = conn_[3];
        new_tets_conn[2 * 4 + 3] = conn_[2];
      }
      // TET3
      new_tets_conn[3 * 4 + 0] = edge_new_nodes_[4];
      new_tets_conn[3 * 4 + 1] = edge_new_nodes_[1];
      if (is_rotated) {
        new_tets_conn[3 * 4 + 2] = conn_[3];
        new_tets_conn[3 * 4 + 3] = conn_[0];
      } else {
        new_tets_conn[3 * 4 + 2] = conn_[0];
        new_tets_conn[3 * 4 + 3] = conn_[3];
      }
      return type;
    } else if (free_edge_swappped0 && (!free_edge_swappped5)) {
      // TET0
      new_tets_conn[0 * 4 + 0] = edge_new_nodes_[4];
      new_tets_conn[0 * 4 + 1] = edge_new_nodes_[1];
      if (is_rotated) {
        new_tets_conn[0 * 4 + 2] = edge_new_nodes_[2];
        new_tets_conn[0 * 4 + 3] = conn_[1];
      } else {
        new_tets_conn[0 * 4 + 2] = conn_[1];
        new_tets_conn[0 * 4 + 3] = edge_new_nodes_[2];
      }
      // TET1
      new_tets_conn[1 * 4 + 0] = edge_new_nodes_[4];
      new_tets_conn[1 * 4 + 1] = edge_new_nodes_[1];
      if (is_rotated) {
        new_tets_conn[1 * 4 + 2] = conn_[2];
        new_tets_conn[1 * 4 + 3] = edge_new_nodes_[2];
      } else {
        new_tets_conn[1 * 4 + 2] = edge_new_nodes_[2];
        new_tets_conn[1 * 4 + 3] = conn_[2];
      }
      // TET2
      new_tets_conn[2 * 4 + 0] = edge_new_nodes_[4];
      new_tets_conn[2 * 4 + 1] = edge_new_nodes_[2];
      if (is_rotated) {
        new_tets_conn[2 * 4 + 2] = conn_[0];
        new_tets_conn[2 * 4 + 3] = conn_[1];
      } else {
        new_tets_conn[2 * 4 + 2] = conn_[1];
        new_tets_conn[2 * 4 + 3] = conn_[0];
      }
      // TET3
      new_tets_conn[3 * 4 + 0] = edge_new_nodes_[4];
      new_tets_conn[3 * 4 + 1] = edge_new_nodes_[2];
      if (is_rotated) {
        new_tets_conn[3 * 4 + 2] = conn_[3];
        new_tets_conn[3 * 4 + 3] = conn_[0];
      } else {
        new_tets_conn[3 * 4 + 2] = conn_[0];
        new_tets_conn[3 * 4 + 3] = conn_[3];
      }
      // TET4
      new_tets_conn[4 * 4 + 0] = edge_new_nodes_[4];
      new_tets_conn[4 * 4 + 1] = edge_new_nodes_[2];
      if (is_rotated) {
        new_tets_conn[4 * 4 + 2] = conn_[2];
        new_tets_conn[4 * 4 + 3] = conn_[3];
      } else {
        new_tets_conn[4 * 4 + 2] = conn_[3];
        new_tets_conn[4 * 4 + 3] = conn_[2];
      }
      return 5;
    } else if (free_edge_swappped0 && free_edge_swappped5) {
      // TET0
      new_tets_conn[0 * 4 + 0] = edge_new_nodes_[4];
      new_tets_conn[0 * 4 + 1] = edge_new_nodes_[1];
      if (is_rotated) {
        new_tets_conn[0 * 4 + 2] = edge_new_nodes_[2];
        new_tets_conn[0 * 4 + 3] = conn_[1];
      } else {
        new_tets_conn[0 * 4 + 2] = conn_[1];
        new_tets_conn[0 * 4 + 3] = edge_new_nodes_[2];
      }
      // TET1
      new_tets_conn[1 * 4 + 0] = edge_new_nodes_[1];
      new_tets_conn[1 * 4 + 1] = edge_new_nodes_[2];
      if (is_rotated) {
        new_tets_conn[1 * 4 + 2] = conn_[2];
        new_tets_conn[1 * 4 + 3] = conn_[3];
      } else {
        new_tets_conn[1 * 4 + 2] = conn_[3];
        new_tets_conn[1 * 4 + 3] = conn_[2];
      }
      // TET2
      new_tets_conn[2 * 4 + 0] = edge_new_nodes_[4];
      new_tets_conn[2 * 4 + 1] = edge_new_nodes_[2];
      if (is_rotated) {
        new_tets_conn[2 * 4 + 2] = conn_[0];
        new_tets_conn[2 * 4 + 3] = conn_[1];
      } else {
        new_tets_conn[2 * 4 + 2] = conn_[1];
        new_tets_conn[2 * 4 + 3] = conn_[0];
      }
      // TET3
      new_tets_conn[3 * 4 + 0] = edge_new_nodes_[4];
      new_tets_conn[3 * 4 + 1] = edge_new_nodes_[2];
      if (is_rotated) {
        new_tets_conn[3 * 4 + 2] = conn_[3];
        new_tets_conn[3 * 4 + 3] = conn_[0];
      } else {
        new_tets_conn[3 * 4 + 2] = conn_[0];
        new_tets_conn[3 * 4 + 3] = conn_[3];
      }
      // TET4
      new_tets_conn[4 * 4 + 0] = edge_new_nodes_[4];
      new_tets_conn[4 * 4 + 1] = edge_new_nodes_[2];
      if (is_rotated) {
        new_tets_conn[4 * 4 + 2] = edge_new_nodes_[1];
        new_tets_conn[4 * 4 + 3] = conn_[3];
      } else {
        new_tets_conn[4 * 4 + 2] = conn_[3];
        new_tets_conn[4 * 4 + 3] = edge_new_nodes_[1];
      }
      return 6;
    }
    // TET0
    new_tets_conn[0 * 4 + 0] = edge_new_nodes_[4];
    new_tets_conn[0 * 4 + 1] = edge_new_nodes_[1];
    if (is_rotated) {
      new_tets_conn[0 * 4 + 2] = conn_[0];
      new_tets_conn[0 * 4 + 3] = conn_[1];
    } else {
      new_tets_conn[0 * 4 + 2] = conn_[1];
      new_tets_conn[0 * 4 + 3] = conn_[0];
    }
    // TET1
    new_tets_conn[1 * 4 + 0] = edge_new_nodes_[1];
    new_tets_conn[1 * 4 + 1] = edge_new_nodes_[2];
    if (is_rotated) {
      new_tets_conn[1 * 4 + 2] = conn_[2];
      new_tets_conn[1 * 4 + 3] = edge_new_nodes_[4];
    } else {
      new_tets_conn[1 * 4 + 2] = edge_new_nodes_[4];
      new_tets_conn[1 * 4 + 3] = conn_[2];
    }
    // TET2
    new_tets_conn[2 * 4 + 0] = edge_new_nodes_[4];
    new_tets_conn[2 * 4 + 1] = edge_new_nodes_[2];
    if (is_rotated) {
      new_tets_conn[2 * 4 + 2] = conn_[2];
      new_tets_conn[2 * 4 + 3] = conn_[3];
    } else {
      new_tets_conn[2 * 4 + 2] = conn_[3];
      new_tets_conn[2 * 4 + 3] = conn_[2];
    }
    // TET3
    new_tets_conn[3 * 4 + 0] = edge_new_nodes_[4];
    new_tets_conn[3 * 4 + 1] = edge_new_nodes_[2];
    if (is_rotated) {
      new_tets_conn[3 * 4 + 2] = conn_[0];
      new_tets_conn[3 * 4 + 3] = edge_new_nodes_[1];
    } else {
      new_tets_conn[3 * 4 + 2] = edge_new_nodes_[1];
      new_tets_conn[3 * 4 + 3] = conn_[0];
    }
    // TET4
    new_tets_conn[4 * 4 + 0] = edge_new_nodes_[4];
    new_tets_conn[4 * 4 + 1] = edge_new_nodes_[2];
    if (is_rotated) {
      new_tets_conn[4 * 4 + 2] = conn_[3];
      new_tets_conn[4 * 4 + 3] = conn_[0];
    } else {
      new_tets_conn[4 * 4 + 2] = conn_[0];
      new_tets_conn[4 * 4 + 3] = conn_[3];
    }
    return 7;
  }
  return type;
}
int tet_type_2(const EntityHandle *conn, const int *split_edges,
               const EntityHandle *edge_new_nodes,
               EntityHandle *new_tets_conn) {
  if (split_edges[0] == oposite_edge[split_edges[1]]) {
    // type 2b
    EntityHandle n0 = conn[edges_conn[2 * split_edges[0] + 0]];
    EntityHandle n1 = conn[edges_conn[2 * split_edges[0] + 1]];
    EntityHandle n2 = conn[edges_conn[2 * split_edges[1] + 0]];
    EntityHandle n3 = conn[edges_conn[2 * split_edges[1] + 1]];
    // TET0
    new_tets_conn[0 * 4 + 0] = edge_new_nodes[split_edges[0]];
    new_tets_conn[0 * 4 + 1] = edge_new_nodes[split_edges[1]];
    new_tets_conn[0 * 4 + 2] = n2;
    new_tets_conn[0 * 4 + 3] = n0;
    // TET1
    new_tets_conn[1 * 4 + 0] = edge_new_nodes[split_edges[0]];
    new_tets_conn[1 * 4 + 1] = edge_new_nodes[split_edges[1]];
    new_tets_conn[1 * 4 + 2] = n1;
    new_tets_conn[1 * 4 + 3] = n2;
    // TET3
    new_tets_conn[2 * 4 + 0] = edge_new_nodes[split_edges[0]];
    new_tets_conn[2 * 4 + 1] = edge_new_nodes[split_edges[1]];
    new_tets_conn[2 * 4 + 2] = n0;
    new_tets_conn[2 * 4 + 3] = n3;
    // TET4
    new_tets_conn[3 * 4 + 0] = edge_new_nodes[split_edges[0]];
    new_tets_conn[3 * 4 + 1] = edge_new_nodes[split_edges[1]];
    new_tets_conn[3 * 4 + 2] = n3;
    new_tets_conn[3 * 4 + 3] = n1;
    return 1;
  } else {
    int sub_type = 0;
    // type 2a
    const char mach_pattern =
        edge_bits_mark[split_edges[0]] | edge_bits_mark[split_edges[1]];
    EntityHandle conn_[4] = {no_handle, no_handle, no_handle, no_handle};
    int edges_[6];
    for (int ee = 0; ee < 6; ee++) {
      char pattern;
      pattern = edge_bits_mark[edge_permutations[ee][1]] |
                edge_bits_mark[edge_permutations[ee][4]];
      if (pattern == mach_pattern) {
        int free_edge = oposite_edge[ee];
        conn_[0] = conn[edges_conn[ee * 2 + 0]];
        conn_[1] = conn[edges_conn[ee * 2 + 1]];
        conn_[2] = conn[edges_conn[free_edge * 2 + 0]];
        conn_[3] = conn[edges_conn[free_edge * 2 + 1]];
        for (int EE = 0; EE < 6; EE++)
          edges_[EE] = edge_permutations[ee][EE];
        sub_type |= 4;
        break;
      }
      pattern = edge_bits_mark[edge_permutations[ee][edge_mirror_cross[1]]] |
                edge_bits_mark[edge_permutations[ee][edge_mirror_cross[4]]];
      if (pattern == mach_pattern) {
        int free_edge = oposite_edge[ee];
        conn_[0] = conn[edges_conn[ee * 2 + 1]];
        conn_[1] = conn[edges_conn[ee * 2 + 0]];
        conn_[2] = conn[edges_conn[free_edge * 2 + 1]];
        conn_[3] = conn[edges_conn[free_edge * 2 + 0]];
        for (int EE = 0; EE < 6; EE++)
          edges_[EE] = edge_permutations[ee][edge_mirror_cross[EE]];
        sub_type |= 8;
        break;
      }
    }
    EntityHandle edge_new_nodes_[6];
    for (int ee = 0; ee < 6; ee++)
      edge_new_nodes_[ee] = edge_new_nodes[edges_[ee]];
    assert(edge_new_nodes_[1] != no_handle);
    assert(edge_new_nodes_[4] != no_handle);
    // TET0
    new_tets_conn[0 * 4 + 0] = conn_[1];
    new_tets_conn[0 * 4 + 1] = edge_new_nodes_[4];
    new_tets_conn[0 * 4 + 2] = edge_new_nodes_[1];
    new_tets_conn[0 * 4 + 3] = conn_[0];
    bool free_edge_swappped5 = false;
    if (conn_[3] < conn_[2]) {
      free_edge_swappped5 = true;
      sub_type |= 16;
    }
    if (free_edge_swappped5) {
      // TET1
      new_tets_conn[1 * 4 + 0] = edge_new_nodes_[1];
      new_tets_conn[1 * 4 + 1] = conn_[2];
      new_tets_conn[1 * 4 + 2] = conn_[0];
      new_tets_conn[1 * 4 + 3] = conn_[3];
      // TET2
      new_tets_conn[2 * 4 + 0] = edge_new_nodes_[1];
      new_tets_conn[2 * 4 + 1] = edge_new_nodes_[4];
      new_tets_conn[2 * 4 + 2] = conn_[3];
      new_tets_conn[2 * 4 + 3] = conn_[0];
      return sub_type;
    }
    // TET1
    new_tets_conn[1 * 4 + 0] = edge_new_nodes_[4];
    new_tets_conn[1 * 4 + 1] = conn_[2];
    new_tets_conn[1 * 4 + 2] = conn_[0];
    new_tets_conn[1 * 4 + 3] = conn_[3];
    // TET2
    new_tets_conn[2 * 4 + 0] = edge_new_nodes_[1];
    new_tets_conn[2 * 4 + 1] = edge_new_nodes_[4];
    new_tets_conn[2 * 4 + 2] = conn_[2];
    new_tets_conn[2 * 4 + 3] = conn_[0];
    return sub_type;
  }
  return -1;
}
void tet_type_1(const EntityHandle *conn, const int split_edge,
                const EntityHandle edge_new_node, EntityHandle *new_tets_conn) {
  // TET0
  new_tets_conn[0 * 4 + 0] = edge_new_node;
  new_tets_conn[0 * 4 + 1] = conn[edges_conn[2 * split_edge + 0]];
  new_tets_conn[0 * 4 + 2] = conn[edges_conn[2 * oposite_edge[split_edge] + 1]];
  new_tets_conn[0 * 4 + 3] = conn[edges_conn[2 * oposite_edge[split_edge] + 0]];
  // TET1
  new_tets_conn[1 * 4 + 0] = edge_new_node;
  new_tets_conn[1 * 4 + 1] = conn[edges_conn[2 * split_edge + 1]];
  new_tets_conn[1 * 4 + 2] = conn[edges_conn[2 * oposite_edge[split_edge] + 0]];
  new_tets_conn[1 * 4 + 3] = conn[edges_conn[2 * oposite_edge[split_edge] + 1]];
}

// TRIS
MoFEMErrorCode tri_type_3(const EntityHandle *conn,
                          const BitRefEdges split_edges,
                          const EntityHandle *edge_new_nodes,
                          EntityHandle *new_tris_conn) {
  MoFEMFunctionBeginHot;
  for (int ee = 0; ee < 3; ee++) {
    if (!split_edges.test(ee)) {
      SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY, "data inconsistency");
    }
  }
  // TRI0
  new_tris_conn[0 * 3 + 0] = edge_new_nodes[0];
  new_tris_conn[0 * 3 + 1] = edge_new_nodes[2];
  new_tris_conn[0 * 3 + 2] = conn[0];
  // TRI1
  new_tris_conn[1 * 3 + 0] = edge_new_nodes[0];
  new_tris_conn[1 * 3 + 1] = edge_new_nodes[1];
  new_tris_conn[1 * 3 + 2] = conn[1];
  // TRI2
  new_tris_conn[2 * 3 + 0] = edge_new_nodes[1];
  new_tris_conn[2 * 3 + 1] = edge_new_nodes[2];
  new_tris_conn[2 * 3 + 2] = conn[2];
  // TRI3
  new_tris_conn[3 * 3 + 0] = edge_new_nodes[0];
  new_tris_conn[3 * 3 + 1] = edge_new_nodes[1];
  new_tris_conn[3 * 3 + 2] = edge_new_nodes[2];
  MoFEMFunctionReturnHot(0);
}

// PRISM
MoFEMErrorCode prism_type_1(const EntityHandle *conn,
                            const BitRefEdges split_edges,
                            const EntityHandle *edge_new_nodes,
                            EntityHandle *new_prism_conn) {
  MoFEMFunctionBeginHot;
  int ee = 0;
  for (; ee < 6; ee++) {
    if (split_edges.test(ee))
      break;
  }
  if (ee > 2)
    SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY, "data inconsistency");
  const int cycle_edges[3][6] = {
      {0, 1, 2, 3, 4, 5}, {1, 2, 0, 4, 5, 3}, {2, 0, 1, 5, 3, 4}};
  const int cycle_nodes[3][6] = {
      {0, 1, 2, 3, 4, 5}, {1, 2, 0, 4, 5, 3}, {2, 0, 1, 5, 3, 4}};
  if (edge_new_nodes[cycle_nodes[ee][0]] == no_handle)
    SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY, "data inconsistency");
  if (edge_new_nodes[cycle_nodes[ee][3]] == no_handle)
    SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY, "data inconsistency");
  // PRISM0
  new_prism_conn[0 * 6 + 0] = edge_new_nodes[cycle_edges[ee][0]];
  new_prism_conn[0 * 6 + 1] = conn[cycle_nodes[ee][2]];
  new_prism_conn[0 * 6 + 2] = conn[cycle_nodes[ee][0]];
  new_prism_conn[0 * 6 + 3] = edge_new_nodes[cycle_edges[ee][3]];
  new_prism_conn[0 * 6 + 4] = conn[cycle_nodes[ee][5]];
  new_prism_conn[0 * 6 + 5] = conn[cycle_nodes[ee][3]];
  // PRISM1
  new_prism_conn[1 * 6 + 0] = edge_new_nodes[cycle_edges[ee][0]];
  new_prism_conn[1 * 6 + 1] = conn[cycle_nodes[ee][2]];
  new_prism_conn[1 * 6 + 2] = conn[cycle_nodes[ee][1]];
  new_prism_conn[1 * 6 + 3] = edge_new_nodes[cycle_edges[ee][3]];
  new_prism_conn[1 * 6 + 4] = conn[cycle_nodes[ee][5]];
  new_prism_conn[1 * 6 + 5] = conn[cycle_nodes[ee][4]];
  MoFEMFunctionReturnHot(0);
}
MoFEMErrorCode prism_type_2(const EntityHandle *conn,
                            const BitRefEdges split_edges,
                            const EntityHandle *edge_new_nodes,
                            EntityHandle *new_prism_conn) {
  MoFEMFunctionBeginHot;
  const int cycle_edges[3][6] = {
      {0, 1, 2, 3, 4, 5}, {1, 2, 0, 4, 5, 3}, {2, 0, 1, 5, 3, 4}};
  const int cycle_nodes[3][6] = {
      {0, 1, 2, 3, 4, 5}, {1, 2, 0, 4, 5, 3}, {2, 0, 1, 5, 3, 4}};
  int ee = 0;
  for (; ee < 3; ee++) {
    BitRefEdges mach_pattern(0);
    mach_pattern.set(cycle_edges[ee][0]);
    mach_pattern.set(cycle_edges[ee][2]);
    mach_pattern.set(cycle_edges[ee][3]);
    mach_pattern.set(cycle_edges[ee][5]);
    if (mach_pattern == split_edges)
      break;
  }
  if (ee > 2)
    SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY, "data inconsistency");
  EntityHandle _conn_[6], _edge_new_nodes_[6];
  int nn = 0;
  for (; nn < 6; nn++) {
    _conn_[nn] = conn[cycle_nodes[ee][nn]];
    _edge_new_nodes_[nn] = edge_new_nodes[cycle_edges[ee][nn]];
  }
  if (_conn_[1] < _conn_[2]) {
    EntityHandle _conn____;
    _conn____ = _conn_[2];
    _conn_[2] = _conn_[1];
    _conn_[1] = _conn____;
    _conn____ = _conn_[5];
    _conn_[5] = _conn_[4];
    _conn_[4] = _conn____;
    _conn____ = _edge_new_nodes_[0];
    _edge_new_nodes_[0] = _edge_new_nodes_[2];
    _edge_new_nodes_[2] = _conn____;
    _conn____ = _edge_new_nodes_[6 - 3];
    _edge_new_nodes_[6 - 3] = _edge_new_nodes_[8 - 3];
    _edge_new_nodes_[8 - 3] = _conn____;
  }
  if (_edge_new_nodes_[0] == no_handle)
    SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY, "data inconsistency");
  if (_edge_new_nodes_[2] == no_handle)
    SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY, "data inconsistency");
  if (_edge_new_nodes_[6 - 3] == no_handle)
    SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY, "data inconsistency");
  if (_edge_new_nodes_[8 - 3] == no_handle)
    SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY, "data inconsistency");
  // PRIMS0
  new_prism_conn[0 * 6 + 0] = _conn_[0];
  new_prism_conn[0 * 6 + 1] = _edge_new_nodes_[0];
  new_prism_conn[0 * 6 + 2] = _edge_new_nodes_[2];
  new_prism_conn[0 * 6 + 3] = _conn_[3];
  new_prism_conn[0 * 6 + 4] = _edge_new_nodes_[6 - 3];
  new_prism_conn[0 * 6 + 5] = _edge_new_nodes_[8 - 3];
  // PRISM1
  new_prism_conn[1 * 6 + 0] = _conn_[1];
  new_prism_conn[1 * 6 + 1] = _conn_[2];
  new_prism_conn[1 * 6 + 2] = _edge_new_nodes_[0];
  new_prism_conn[1 * 6 + 3] = _conn_[4];
  new_prism_conn[1 * 6 + 4] = _conn_[5];
  new_prism_conn[1 * 6 + 5] = _edge_new_nodes_[6 - 3];
  // PRISM2
  new_prism_conn[2 * 6 + 0] = _conn_[2];
  new_prism_conn[2 * 6 + 1] = _edge_new_nodes_[0];
  new_prism_conn[2 * 6 + 2] = _edge_new_nodes_[2];
  new_prism_conn[2 * 6 + 3] = _conn_[5];
  new_prism_conn[2 * 6 + 4] = _edge_new_nodes_[6 - 3];
  new_prism_conn[2 * 6 + 5] = _edge_new_nodes_[8 - 3];
  MoFEMFunctionReturnHot(0);
}
MoFEMErrorCode prism_type_3(const EntityHandle *conn,
                            const BitRefEdges split_edges,
                            const EntityHandle *edge_new_nodes,
                            EntityHandle *new_prism_conn) {
  MoFEMFunctionBeginHot;
  int ee = 0;
  for (; ee < 6; ee++) {
    if (!split_edges.test(ee))
      SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY, "data inconsistency");
  }
  // PRISM0
  new_prism_conn[0 * 6 + 0] = edge_new_nodes[0];
  new_prism_conn[0 * 6 + 1] = edge_new_nodes[2];
  new_prism_conn[0 * 6 + 2] = conn[0];
  new_prism_conn[0 * 6 + 3] = edge_new_nodes[6 - 3];
  new_prism_conn[0 * 6 + 4] = edge_new_nodes[8 - 3];
  new_prism_conn[0 * 6 + 5] = conn[3];
  // PRISM1
  new_prism_conn[1 * 6 + 0] = edge_new_nodes[0];
  new_prism_conn[1 * 6 + 1] = edge_new_nodes[1];
  new_prism_conn[1 * 6 + 2] = conn[1];
  new_prism_conn[1 * 6 + 3] = edge_new_nodes[6 - 3];
  new_prism_conn[1 * 6 + 4] = edge_new_nodes[7 - 3];
  new_prism_conn[1 * 6 + 5] = conn[4];
  // PRISM2
  new_prism_conn[2 * 6 + 0] = edge_new_nodes[1];
  new_prism_conn[2 * 6 + 1] = edge_new_nodes[2];
  new_prism_conn[2 * 6 + 2] = conn[2];
  new_prism_conn[2 * 6 + 3] = edge_new_nodes[7 - 3];
  new_prism_conn[2 * 6 + 4] = edge_new_nodes[8 - 3];
  new_prism_conn[2 * 6 + 5] = conn[5];
  // PRISM3
  new_prism_conn[3 * 6 + 0] = edge_new_nodes[0];
  new_prism_conn[3 * 6 + 1] = edge_new_nodes[1];
  new_prism_conn[3 * 6 + 2] = edge_new_nodes[2];
  new_prism_conn[3 * 6 + 3] = edge_new_nodes[6 - 3];
  new_prism_conn[3 * 6 + 4] = edge_new_nodes[7 - 3];
  new_prism_conn[3 * 6 + 5] = edge_new_nodes[8 - 3];
  MoFEMFunctionReturnHot(0);
}

MoFEMErrorCode quad_split_all_edges(const EntityHandle *conn,
                                    const BitRefEdges split_edges,
                                    const EntityHandle *edge_new_nodes,
                                    EntityHandle *new_quad_conn) {




                                    }

} // namespace MoFEM
