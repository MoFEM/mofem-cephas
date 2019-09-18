/** \file edge_and_bubble_shape_functions_on_quad.cpp
  \example edge_and_bubble_shape_functions_on_quad.cpp
  \brief Testing edge and bubble shape functions on a quad face

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

#include <MoFEM.hpp>

using namespace MoFEM;
static char help[] = "...\n\n";

static double sum_matrix(MatrixDouble &m) {
  double s = 0;
  for (unsigned int ii = 0; ii < m.size1(); ii++) {
    for (unsigned int jj = 0; jj < m.size2(); jj++) {
      s += m(ii, jj);
    }
  }
  return s;
}

int main(int argc, char *argv[]) {

  MoFEM::Core::Initialize(&argc, &argv, (char *)0, help);
  try {
    moab::Core mb_instance;
    moab::Interface &moab = mb_instance;
    ParallelComm *pcomm = ParallelComm::get_pcomm(&moab, MYPCOMM_INDEX);
    if (pcomm == NULL)
      pcomm = new ParallelComm(&moab, PETSC_COMM_WORLD);
    MoFEM::Core core(moab);
    MoFEM::Interface &m_field = core;

    PetscBool flg = PETSC_TRUE;
    char mesh_file_name[255];
#if PETSC_VERSION_GE(3, 6, 4)
    CHKERR PetscOptionsGetString(PETSC_NULL, "", "-my_file", mesh_file_name,
                                 255, &flg);
#else
    CHKERR PetscOptionsGetString(PETSC_NULL, PETSC_NULL, "-my_file",
                                 mesh_file_name, 255, &flg);
#endif
    if (flg != PETSC_TRUE) {
      SETERRQ(PETSC_COMM_SELF, 1, "*** ERROR -my_file (MESH FILE NEEDED)");
    }

    const char *option;
    option = ""; //"PARALLEL=BCAST;";//;DEBUG_IO";
    CHKERR moab.load_file(mesh_file_name, 0, option);
    Range verts;
    CHKERR moab.get_entities_by_type(0, MBVERTEX, verts, true);

    MatrixDouble coords(verts.size(), 3);
    MatrixDouble N(verts.size(), 4);
    MatrixDouble diffN(verts.size(), 8);
    CHKERR moab.get_coords(verts, &coords(0, 0));
    for (int i = 0; i < verts.size(); ++i) {
      N(i, 0) = N_MBQUAD0(coords(i, 0), coords(i, 1));
      N(i, 1) = N_MBQUAD1(coords(i, 0), coords(i, 1));
      N(i, 2) = N_MBQUAD2(coords(i, 0), coords(i, 1));
      N(i, 3) = N_MBQUAD3(coords(i, 0), coords(i, 1));

      diffN(i, 0) = diffN_MBQUAD0x(coords(i, 1));
      diffN(i, 1) = diffN_MBQUAD0y(coords(i, 0));
      diffN(i, 2) = diffN_MBQUAD1x(coords(i, 1));
      diffN(i, 3) = diffN_MBQUAD1y(coords(i, 0));
      diffN(i, 4) = diffN_MBQUAD2x(coords(i, 1));
      diffN(i, 5) = diffN_MBQUAD2y(coords(i, 0));
      diffN(i, 6) = diffN_MBQUAD3x(coords(i, 1));
      diffN(i, 7) = diffN_MBQUAD3y(coords(i, 0));
    }

    /* BUBBLES */
    {
      std::array<int, 4> faces_nodes = {0, 1, 2, 3};
      int p = 7;
      int P = NBFACEQUAD_H1(p);
      MatrixDouble quad_bubbles(verts.size(), P);
      MatrixDouble quad_diff_bubbles(verts.size(), P * 2);
      double eps = 1e-8;
      double quad_bubbles_sum = 1.56816;
      double quad_diff_bubbles_sum = -7.57944;

      CHKERR H1_QuadShapeFunctions_MBQUAD(
          faces_nodes.data(), p, &*N.data().begin(), &*diffN.data().begin(),
          &*quad_bubbles.data().begin(), &*quad_diff_bubbles.data().begin(),
          verts.size(), Legendre_polynomials);

      if (fabs(quad_bubbles_sum - sum_matrix(quad_bubbles)) > eps) {
        SETERRQ(PETSC_COMM_WORLD, MOFEM_ATOM_TEST_INVALID, "wrong result");
      }
      if (fabs(quad_diff_bubbles_sum - sum_matrix(quad_diff_bubbles)) > eps) {
        SETERRQ(PETSC_COMM_WORLD, MOFEM_ATOM_TEST_INVALID, "wrong result");
      }
    }

    /* EDGES */
    {
      int sense[] = {1, -1, 1, -1};
      int p[] = {3, 4, 5, 6};
      int P[4];
      int ee = 0;
      for (; ee < 4; ee++) {
        P[ee] = NBEDGE_H1(p[ee]);
      }

      double *quad_edges_ptr[4];
      std::array<MatrixDouble, 4> quad_edges;
      for (auto ee : {0, 1, 2, 3}) {
        quad_edges[ee] = MatrixDouble(verts.size(), P[ee]);
        quad_edges_ptr[ee] = &quad_edges[ee](0, 0);
      }

      double *quad_diff_edges_ptr[4];
      std::array<MatrixDouble, 4> quad_diff_edges;
      for (auto ee : {0, 1, 2, 3}) {
        quad_diff_edges[ee] = MatrixDouble(verts.size(), P[ee] * 2);
        quad_diff_edges_ptr[ee] = &quad_diff_edges[ee](0, 0);
      }

      double eps = 1e-8;
      double quad_edges_sum[] = {6.3525, 4.38007416, 4.38007416, 4.849528992};
      double quad_diff_edges_sum[] = {-21.4775, 18.15242, 19.96258,
                                      -19.7395528};

      CHKERR H1_EdgeShapeFunctions_MBQUAD(
          sense, p, &*N.data().begin(), &*diffN.data().begin(), quad_edges_ptr,
          quad_diff_edges_ptr, verts.size(), Legendre_polynomials);

      for (int ee = 0; ee < 4; ++ee) {
        if (fabs(quad_edges_sum[ee] - sum_matrix(quad_edges[ee])) > eps) {
          SETERRQ(PETSC_COMM_WORLD, MOFEM_ATOM_TEST_INVALID, "wrong result");
        }
        if (fabs(quad_diff_edges_sum[ee] - sum_matrix(quad_diff_edges[ee])) >
            eps) {
          SETERRQ(PETSC_COMM_WORLD, MOFEM_ATOM_TEST_INVALID, "wrong result");
        }
      }
    }
  }
  CATCH_ERRORS;

  CHKERR MoFEM::Core::Finalize();
}
