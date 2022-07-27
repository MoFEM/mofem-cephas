/** \file edge_and_bubble_shape_functions_on_quad.cpp
  \example edge_and_bubble_shape_functions_on_quad.cpp
  \brief Testing edge and bubble shape functions on a quad face

*/

#include <MoFEM.hpp>

using namespace MoFEM;
static char help[] = "...\n\n";

static double sum_matrix(MatrixDouble &m) {
  double s = 0;
  for (unsigned int ii = 0; ii < m.size1(); ii++) {
    for (unsigned int jj = 0; jj < m.size2(); jj++) {
      s += std::abs(m(ii, jj));
    }
  }
  return s;
}

int main(int argc, char *argv[]) {

  MoFEM::Core::Initialize(&argc, &argv, (char *)0, help);
  try {
    moab::Core mb_instance;
    moab::Interface &moab = mb_instance;
    
    MoFEM::Core core(moab);

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
      SETERRQ(PETSC_COMM_SELF, 1, "Error -my_file (mesh file needed)");
    }

    const char *option;
    option = ""; 
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
      constexpr double eps = 1e-4;
      double quad_bubbles_sum = 3.275491e+01;
      double quad_diff_bubbles_sum = 5.838131e+02;

      CHKERR H1_QuadShapeFunctions_MBQUAD(
          faces_nodes.data(), p, &*N.data().begin(), &*diffN.data().begin(),
          &*quad_bubbles.data().begin(), &*quad_diff_bubbles.data().begin(),
          verts.size(), Legendre_polynomials);

      if (std::abs(quad_bubbles_sum - sum_matrix(quad_bubbles)) > eps) {
        SETERRQ1(PETSC_COMM_WORLD, MOFEM_ATOM_TEST_INVALID,
                 "wrong result = %8.6e", sum_matrix(quad_bubbles));
      }
      if (std::abs(quad_diff_bubbles_sum - sum_matrix(quad_diff_bubbles)) >
          eps) {
        SETERRQ1(PETSC_COMM_WORLD, MOFEM_ATOM_TEST_INVALID,
                 "wrong result = %8.6e", sum_matrix(quad_diff_bubbles));
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

      double eps = 1e-4;
      double quad_edges_sum[] = {1.237500e+01, 1.535050e+01, 1.781890e+01,
                                 2.015678e+01};
      double quad_diff_edges_sum[] = {8.470000e+01, 1.192510e+02, 1.580128e+02,
                                      1.976378e+02};

      CHKERR H1_EdgeShapeFunctions_MBQUAD(
          sense, p, &*N.data().begin(), &*diffN.data().begin(), quad_edges_ptr,
          quad_diff_edges_ptr, verts.size(), Legendre_polynomials);

      for (int ee = 0; ee < 4; ++ee) {
        if (std::abs(quad_edges_sum[ee] - sum_matrix(quad_edges[ee])) > eps)
          SETERRQ2(PETSC_COMM_WORLD, MOFEM_ATOM_TEST_INVALID,
                   "edge %d wrong result = %8.6e", ee,
                   sum_matrix(quad_edges[ee]));

        if (std::abs(quad_diff_edges_sum[ee] -
                     sum_matrix(quad_diff_edges[ee])) > eps)
          SETERRQ2(PETSC_COMM_WORLD, MOFEM_ATOM_TEST_INVALID,
                   "edge %d wrong result = %8.6e", ee,
                   sum_matrix(quad_diff_edges[ee]));
      }
    }
  }
  CATCH_ERRORS;

  CHKERR MoFEM::Core::Finalize();
}
