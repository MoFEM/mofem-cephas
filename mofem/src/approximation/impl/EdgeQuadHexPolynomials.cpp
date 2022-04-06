/** \file EdgeQuadHexPolynomials.cpp

  \brief Implementation of hierarchical Edge, Quad, and Hex shape functions of
  type H1, Hcurl, Hdiv
*/

using namespace MoFEM;

namespace DemkowiczHexAndQuad {

MoFEMErrorCode monom_ordering(int *perm, int p, int q, int r = 0) {
  MoFEMFunctionBeginHot;

  if (r > 0) {

    for (int m = 0; m != std::max(std::max(p, q), r) + 1; ++m) {

      const int i = std::min(m, p);
      const int j = std::min(m, q);
      const int k = std::min(m, r);

      if (i == m)
        for (int jj = 0; jj != j; ++jj) {
          for (int kk = 0; kk != k; ++kk) {
            *(perm++) = i;
            *(perm++) = jj;
            *(perm++) = kk;
          }
        }

      if (j == m)
        for (int ii = 0; ii != i; ++ii) {
          for (int kk = 0; kk != k; ++kk) {
            *(perm++) = ii;
            *(perm++) = j;
            *(perm++) = kk;
          }
        }

      if (k == m)
        for (int ii = 0; ii != i; ++ii) {
          for (int jj = 0; jj != j; ++jj) {
            *(perm++) = ii;
            *(perm++) = jj;
            *(perm++) = k;
          }
        }

      if (j == m || k == m)
        for (int ii = 0; ii != i; ++ii) {
          *(perm++) = ii;
          *(perm++) = j;
          *(perm++) = k;
        }

      if (i == m || k == m)
        for (int jj = 0; jj != j; ++jj) {
          *(perm++) = i;
          *(perm++) = jj;
          *(perm++) = k;
        }

      if (i == m || j == m)
        for (int kk = 0; kk != k; ++kk) {
          *(perm++) = i;
          *(perm++) = j;
          *(perm++) = kk;
        }

      *(perm++) = i;
      *(perm++) = j;
      *(perm++) = k;
    }
  } else {

    for (int m = 0; m != std::max(p, q) + 1; ++m) {

      const int i = std::min(m, p);
      const int j = std::min(m, q);

      if (j == m)
        for (int ii = 0; ii != i; ++ii) {
          *(perm++) = ii;
          *(perm++) = j;
          *(perm++) = 0;
        }

      if (i == m)
        for (int jj = 0; jj != j; ++jj) {
          *(perm++) = i;
          *(perm++) = jj;
          *(perm++) = 0;
        }

      *(perm++) = i;
      *(perm++) = j;
      *(perm++) = 0;
    }
  }
  MoFEMFunctionReturnHot(0);
}

static inline void get_ksi_hex(int shift, double *N, double *N_diff,
                               double ksi[3], double diff_ksi[3][3]) {

  constexpr std::array<size_t, 4> ksi_nodes[2][3] = {

      {{1, 2, 6, 5}, {3, 2, 6, 7}, {4, 5, 6, 7}},

      {{0, 3, 7, 4}, {0, 1, 5, 4}, {0, 1, 2, 3}}

  };

  for (size_t i = 0; i != 3; ++i) {
    for (auto n : ksi_nodes[0][i])
      ksi[i] += N[shift + n];
    for (auto n : ksi_nodes[1][i])
      ksi[i] -= N[shift + n];
  }

  for (size_t i = 0; i != 3; ++i) {
    for (auto d = 0; d != 3; ++d) {
      for (auto n : ksi_nodes[0][i]) {

        diff_ksi[i][d] += N_diff[3 * (shift + n) + d];
      }
    }
    for (auto n : ksi_nodes[1][i]) {
      for (auto d = 0; d != 3; ++d) {
        diff_ksi[i][d] -= N_diff[3 * (shift + n) + d];
      }
    }
  }
}

} // namespace DemkowiczHexAndQuad

/*
    0 *------------* 1
*/
MoFEMErrorCode MoFEM::DemkowiczHexAndQuad::H1_BubbleShapeFunctions_ONSEGMENT(
    int p, double *N, double *diffN, double *bubbleN, double *diff_bubbleN,
    int nb_integration_pts) {
  MoFEMFunctionBeginHot;
  const int nb_dofs = (p - 1);
  if (nb_dofs > 0) {
    constexpr int n0 = 0;
    constexpr int n1 = 1;
    double diff_mu = diffN[n1] - diffN[n0];
    for (int q = 0; q != nb_integration_pts; q++) {
      int shift = 2 * q;
      double mu = N[shift + n1] - N[shift + n0];
      double L[p + 2];
      double diffL[p + 2];
      CHKERR Lobatto_polynomials(p + 1, mu, &diff_mu, L, diffL, 1);
      int qd_shift = nb_dofs * q;
      for (int n = 0; n != nb_dofs; n++) {
        bubbleN[qd_shift + n] = L[n + 2];
        diff_bubbleN[qd_shift + n] = diffL[n + 2];
      }
    }
  }
  MoFEMFunctionReturnHot(0);
}

MoFEMErrorCode MoFEM::DemkowiczHexAndQuad::L2_ShapeFunctions_ONSEGMENT(
    int p, double *N, double *diffN, double *funN, double *funDiffN,
    int nb_integration_pts) {
  MoFEMFunctionBeginHot;
  const int nb_dofs = p + 1;
  if (nb_dofs > 0) {
    constexpr int n0 = 0;
    constexpr int n1 = 1;
    double diff_mu = diffN[n1] - diffN[n0];
    for (int q = 0; q != nb_integration_pts; q++) {
      int shift = 2 * q;
      double mu = N[shift + n1] - N[shift + n0];
      double L[p + 1];
      double diffL[p + 1];
      CHKERR Legendre_polynomials(p, mu, &diff_mu, L, diffL, 1);
      int qd_shift = (p + 1) * q;
      for (int n = 0; n <= p; n++) {
        funN[qd_shift + n] = L[n];
        funDiffN[qd_shift + n] = diffL[n];
      }
    }
  }
  MoFEMFunctionReturnHot(0);
}

/*
      Quads
 3-------2------2
 |              |       eta
 |              |       ^
 3              1       |
 |              |       |
 |              |       0-----  > ksi
 0-------0------1
*/

MoFEMErrorCode MoFEM::DemkowiczHexAndQuad::H1_EdgeShapeFunctions_ONQUAD(
    int *sense, int *p, double *N, double *diffN, double *edgeN[4],
    double *diff_edgeN[4], int nb_integration_pts) {
  MoFEMFunctionBeginHot;

  constexpr int n0 = 0;
  constexpr int n1 = 1;
  constexpr int n2 = 2;
  constexpr int n3 = 3;

  for (int q = 0; q != nb_integration_pts; q++) {
    const int shift = 4 * q;
    const double shape0 = N[shift + n0];
    const double shape1 = N[shift + n1];
    const double shape2 = N[shift + n2];
    const double shape3 = N[shift + n3];
    const double ksi01 = (shape1 + shape2 - shape0 - shape3) * sense[n0];
    const double ksi12 = (shape2 + shape3 - shape1 - shape0) * sense[n1];
    const double ksi23 = (shape3 + shape0 - shape2 - shape1) * sense[n2];
    const double ksi30 = (shape0 + shape1 - shape3 - shape2) * sense[n3];

    const double mu[] = {ksi01, ksi12, ksi23, ksi30};
    const double mu_const[] = {shape0 + shape1, shape1 + shape2,
                               shape2 + shape3, shape3 + shape0};

    const int diff_shift = 2 * shift;
    double diff_mu_const[4][2];
    double diff_mu[4][2];
    for (int d = 0; d != 2; d++) {
      const double diff_shape0 = diffN[diff_shift + 2 * n0 + d];
      const double diff_shape1 = diffN[diff_shift + 2 * n1 + d];
      const double diff_shape2 = diffN[diff_shift + 2 * n2 + d];
      const double diff_shape3 = diffN[diff_shift + 2 * n3 + d];
      diff_mu_const[0][d] = diff_shape0 + diff_shape1;
      diff_mu_const[1][d] = diff_shape1 + diff_shape2;
      diff_mu_const[2][d] = diff_shape2 + diff_shape3;
      diff_mu_const[3][d] = diff_shape3 + diff_shape0;

      const double diff_ksi01 =
          (diff_shape1 + diff_shape2 - diff_shape0 - diff_shape3) * sense[0];
      const double diff_ksi12 =
          (diff_shape2 + diff_shape3 - diff_shape1 - diff_shape0) * sense[1];
      const double diff_ksi23 =
          (diff_shape3 + diff_shape0 - diff_shape2 - diff_shape1) * sense[2];
      const double diff_ksi30 =
          (diff_shape0 + diff_shape1 - diff_shape3 - diff_shape2) * sense[3];
      diff_mu[0][d] = diff_ksi01;
      diff_mu[1][d] = diff_ksi12;
      diff_mu[2][d] = diff_ksi23;
      diff_mu[3][d] = diff_ksi30;
    }

    for (int e = 0; e != 4; e++) {
      const int nb_dofs = p[e] - 1;
      if (nb_dofs > 0) {
        double L[p[e] + 2];
        double diffL[2 * (p[e] + 2)];
        CHKERR Lobatto_polynomials(p[e] + 1, mu[e], diff_mu[e], L, diffL, 2);
        int qd_shift = (p[e] - 1) * q;
        for (int n = 0; n != p[e] - 1; ++n) {
          edgeN[e][qd_shift + n] = mu_const[e] * L[n + 2];
          for (int d = 0; d != 2; ++d) {
            diff_edgeN[e][2 * qd_shift + 2 * n + d] =
                mu_const[e] * diffL[d * (p[e] + 2) + n + 2] +
                diff_mu_const[e][d] * L[n + 2];
          }
        }
      }
    }
  }
  MoFEMFunctionReturnHot(0);
}

MoFEMErrorCode MoFEM::DemkowiczHexAndQuad::H1_FaceShapeFunctions_ONQUAD(
    int *faces_nodes, int *p, double *N, double *diffN, double *faceN,
    double *diff_faceN, int nb_integration_pts) {

  const int nb_dofs = (p[0] - 1) * (p[1] - 1);

  if (nb_dofs > 0) {

    const int n0 = faces_nodes[0];
    const int n1 = faces_nodes[1];
    const int n2 = faces_nodes[2];
    const int n3 = faces_nodes[3];

    MoFEMFunctionBeginHot;
    int permute[(p[0] - 1) * (p[1] - 1)][3];
    CHKERR ::DemkowiczHexAndQuad::monom_ordering(&permute[0][0], p[0] - 2,
                                                 p[1] - 2);
    for (int q = 0; q != nb_integration_pts; q++) {

      const int shift = 4 * q;
      const double shape0 = N[shift + n0];
      const double shape1 = N[shift + n1];
      const double shape2 = N[shift + n2];
      const double shape3 = N[shift + n3];
      const double ksi01 = (shape1 + shape2 - shape0 - shape3);
      const double ksi12 = (shape2 + shape3 - shape1 - shape0);

      const int diff_shift = 2 * shift;
      double diff_ksi01[2], diff_ksi12[2];
      for (int d = 0; d != 2; d++) {
        const double diff_shape0 = diffN[diff_shift + 2 * n0 + d];
        const double diff_shape1 = diffN[diff_shift + 2 * n1 + d];
        const double diff_shape2 = diffN[diff_shift + 2 * n2 + d];
        const double diff_shape3 = diffN[diff_shift + 2 * n3 + d];
        diff_ksi01[d] = (diff_shape1 + diff_shape2 - diff_shape0 - diff_shape3);
        diff_ksi12[d] = (diff_shape2 + diff_shape3 - diff_shape1 - diff_shape0);
      }

      double L01[p[0] + 2];
      double diffL01[2 * (p[0] + 2)];
      CHKERR Lobatto_polynomials(p[0] + 1, ksi01, diff_ksi01, L01, diffL01, 2);
      double L12[p[1] + 2];
      double diffL12[2 * (p[1] + 2)];
      CHKERR Lobatto_polynomials(p[1] + 1, ksi12, diff_ksi12, L12, diffL12, 2);

      int qd_shift = nb_dofs * q;
      for (int n = 0; n != nb_dofs; ++n) {
        int s1 = permute[n][0];
        int s2 = permute[n][1];
        faceN[qd_shift + n] = L01[s1 + 2] * L12[s2 + 2];
        for (int d = 0; d != 2; ++d) {
          diff_faceN[2 * (qd_shift + n) + d] =
              diffL01[d * (p[0] + 2) + s1 + 2] * L12[s2 + 2] +
              L01[s1 + 2] * diffL12[d * (p[1] + 2) + s2 + 2];
        }
      }
    }
  }
  MoFEMFunctionReturnHot(0);
}

MoFEMErrorCode MoFEM::DemkowiczHexAndQuad::L2_FaceShapeFunctions_ONQUAD(
    int *p, double *N, double *diffN, double *faceN, double *diff_faceN,
    int nb_integration_pts) {
  MoFEMFunctionBeginHot;
  const int nb_dofs = (p[0] + 1) * (p[1] + 1);
  if (nb_dofs > 0) {

    int permute[nb_dofs][3];
    CHKERR ::DemkowiczHexAndQuad::monom_ordering(&permute[0][0], p[0], p[1]);

    constexpr int n0 = 0;
    constexpr int n1 = 1;
    constexpr int n2 = 2;
    constexpr int n3 = 3;

    for (int q = 0; q != nb_integration_pts; q++) {
      const int shift = 4 * q;
      const double shape0 = N[shift + n0];
      const double shape1 = N[shift + n1];
      const double shape2 = N[shift + n2];
      const double shape3 = N[shift + n3];
      const double ksi01 = (shape1 + shape2 - shape0 - shape3);
      const double ksi12 = (shape2 + shape3 - shape1 - shape0);

      const int diff_shift = 2 * shift;
      double diff_ksi01[2], diff_ksi12[2];
      for (int d = 0; d != 2; d++) {
        const double diff_shape0 = diffN[diff_shift + 2 * n0 + d];
        const double diff_shape1 = diffN[diff_shift + 2 * n1 + d];
        const double diff_shape2 = diffN[diff_shift + 2 * n2 + d];
        const double diff_shape3 = diffN[diff_shift + 2 * n3 + d];
        diff_ksi01[d] = (diff_shape1 + diff_shape2 - diff_shape0 - diff_shape3);
        diff_ksi12[d] = (diff_shape2 + diff_shape3 - diff_shape1 - diff_shape0);
      }

      double L01[p[0] + 2];
      double diffL01[2 * (p[0] + 2)];
      CHKERR Legendre_polynomials(p[0] + 1, ksi01, diff_ksi01, L01, diffL01, 2);
      double L12[p[1] + 2];
      double diffL12[2 * (p[1] + 2)];
      CHKERR Legendre_polynomials(p[1] + 1, ksi12, diff_ksi12, L12, diffL12, 2);

      int qd_shift = nb_dofs * q;
      for (int n = 0; n != nb_dofs; ++n) {
        int s1 = permute[n][0];
        int s2 = permute[n][1];
        faceN[qd_shift + n] = L01[s1] * L12[s2];
        for (int d = 0; d != 2; ++d) {
          diff_faceN[2 * (qd_shift + n) + d] =
              diffL01[d * (p[0] + 2) + s1] * L12[s2] +
              L01[s1] * diffL12[d * (p[1] + 2) + s2];
        }
      }
    }
  }
  MoFEMFunctionReturnHot(0);
}

MoFEMErrorCode MoFEM::DemkowiczHexAndQuad::Hcurl_EdgeShapeFunctions_ONQUAD(
    int *sense, int *p, double *N, double *diffN, double *edgeN[4],
    double *diff_edgeN[4], int nb_integration_pts) {
  MoFEMFunctionBeginHot;

  constexpr int n0 = 0;
  constexpr int n1 = 1;
  constexpr int n2 = 2;
  constexpr int n3 = 3;

  for (int q = 0; q != nb_integration_pts; q++) {

    const int shift = 4 * q;
    const double shape0 = N[shift + n0];
    const double shape1 = N[shift + n1];
    const double shape2 = N[shift + n2];
    const double shape3 = N[shift + n3];
    const double ksi01 = (shape1 + shape2 - shape0 - shape3) * sense[n0];
    const double ksi12 = (shape2 + shape3 - shape1 - shape0) * sense[n1];
    const double ksi23 = (shape3 + shape0 - shape2 - shape1) * sense[n2];
    const double ksi30 = (shape0 + shape1 - shape3 - shape2) * sense[n3];

    const double mu[] = {ksi01, ksi12, ksi23, ksi30};
    const double mu_const[] = {shape0 + shape1, shape1 + shape2,
                               shape2 + shape3, shape3 + shape0};

    const int diff_shift = 2 * shift;
    double diff_mu_const[4][2];
    double diff_mu[4][2];
    for (int d = 0; d != 2; d++) {
      const double diff_shape0 = diffN[diff_shift + 2 * n0 + d];
      const double diff_shape1 = diffN[diff_shift + 2 * n1 + d];
      const double diff_shape2 = diffN[diff_shift + 2 * n2 + d];
      const double diff_shape3 = diffN[diff_shift + 2 * n3 + d];
      diff_mu_const[0][d] = diff_shape0 + diff_shape1;
      diff_mu_const[1][d] = diff_shape1 + diff_shape2;
      diff_mu_const[2][d] = diff_shape2 + diff_shape3;
      diff_mu_const[3][d] = diff_shape3 + diff_shape0;

      const double diff_ksi01 =
          (diff_shape1 + diff_shape2 - diff_shape0 - diff_shape3) * sense[n0];
      const double diff_ksi12 =
          (diff_shape2 + diff_shape3 - diff_shape1 - diff_shape0) * sense[n1];
      const double diff_ksi23 =
          (diff_shape3 + diff_shape0 - diff_shape2 - diff_shape1) * sense[n2];
      const double diff_ksi30 =
          (diff_shape0 + diff_shape1 - diff_shape3 - diff_shape2) * sense[n3];
      diff_mu[0][d] = diff_ksi01;
      diff_mu[1][d] = diff_ksi12;
      diff_mu[2][d] = diff_ksi23;
      diff_mu[3][d] = diff_ksi30;
    }

    for (int e = 0; e != 4; e++) {

      if (p[e] > 0) {
        double L[p[e]];
        double diffL[2 * p[e]];

        CHKERR Legendre_polynomials(p[e] - 1, mu[e], diff_mu[e], L, diffL, 2);

        int qd_shift = p[e] * q;
        double *t_n_ptr = &edgeN[e][3 * qd_shift];
        double *t_diff_n_ptr = &diff_edgeN[e][3 * 2 * qd_shift];
        auto t_n = getFTensor1FromPtr<3>(t_n_ptr);
        auto t_diff_n = getFTensor2FromPtr<3, 2>(t_diff_n_ptr);

        for (int n = 0; n != p[e]; ++n) {
          const double a = mu_const[e] * L[n];
          const double d_a[] = {
              diff_mu_const[e][0] * L[n] + mu_const[e] * diffL[0 * p[e] + n],
              diff_mu_const[e][1] * L[n] + mu_const[e] * diffL[1 * p[e] + n]};

          for (int d = 0; d != 2; ++d) {
            t_n(d) = (a / 2) * diff_mu[e][d];
            for (int j = 0; j != 2; ++j) {
              t_diff_n(d, j) = (d_a[j] / 2) * diff_mu[e][d];
            }
          }
          t_n(2) = 0;
          t_diff_n(2, 0) = 0;
          t_diff_n(2, 1) = 0;
          ++t_n;
          ++t_diff_n;
        }
      }
    }
  }
  MoFEMFunctionReturnHot(0);
}

MoFEMErrorCode MoFEM::DemkowiczHexAndQuad::Hcurl_FaceShapeFunctions_ONQUAD(
    int *face_nodes, int *p, double *N, double *diffN, double *faceN[],
    double *diff_faceN[], int nb_integration_pts) {
  MoFEMFunctionBeginHot;

  const int pq[2] = {p[0], p[1]};
  const int qp[2] = {p[1], p[0]};

  const int nb_dofs_fm0 = pq[0] * (qp[1] - 1);
  const int nb_dofs_fm1 = (pq[0] - 1) * qp[1];
  int permute_fm0[3 * nb_dofs_fm0];
  int permute_fm1[3 * nb_dofs_fm1];

  std::array<int *, 2> permute = {permute_fm0, permute_fm1};
  for (int fm = 0; fm != 2; ++fm) {
    const int pp = pq[fm];
    const int qq = qp[fm];
    CHKERR ::DemkowiczHexAndQuad::monom_ordering(permute[fm], qq - 1, pp - 2);
  }

  const int n0 = face_nodes[0];
  const int n1 = face_nodes[1];
  const int n2 = face_nodes[2];
  const int n3 = face_nodes[3];

  for (int q = 0; q != nb_integration_pts; q++) {

    const int shift = 4 * q;
    const double shape0 = N[shift + n0];
    const double shape1 = N[shift + n1];
    const double shape2 = N[shift + n2];
    const double shape3 = N[shift + n3];
    const double ksi01 = (shape1 + shape2 - shape0 - shape3);
    const double ksi12 = (shape2 + shape3 - shape1 - shape0);

    const int diff_shift = 2 * shift;
    double diff_ksi01[2], diff_ksi12[2];
    for (int d = 0; d != 2; d++) {
      const double diff_shape0 = diffN[diff_shift + 2 * n0 + d];
      const double diff_shape1 = diffN[diff_shift + 2 * n1 + d];
      const double diff_shape2 = diffN[diff_shift + 2 * n2 + d];
      const double diff_shape3 = diffN[diff_shift + 2 * n3 + d];
      diff_ksi01[d] = (diff_shape1 + diff_shape2 - diff_shape0 - diff_shape3);
      diff_ksi12[d] = (diff_shape2 + diff_shape3 - diff_shape1 - diff_shape0);
    }

    double mu_ksi_eta[2] = {ksi01, ksi12};
    double mu_eta_ksi[2] = {ksi12, ksi01};

    double *diff_ksi_eta[2] = {diff_ksi01, diff_ksi12};
    double *diff_eta_ksi[2] = {diff_ksi12, diff_ksi01};

    for (int family = 0; family != 2; family++) {
      const int pp = pq[family];
      const int qq = qp[family];
      const int nb_dofs = (pp - 1) * qq;

      if (nb_dofs > 0) {

        double zeta[pp + 2];
        double diff_zeta[2 * (pp + 2)];
        CHKERR Lobatto_polynomials(pp + 1, mu_ksi_eta[family],
                                   diff_ksi_eta[family], zeta, diff_zeta, 2);

        double eta[qq];
        double diff_eta[2 * qq];
        CHKERR Legendre_polynomials(qq - 1, mu_eta_ksi[family],
                                    diff_eta_ksi[family], eta, diff_eta, 2);

        const int qd_shift = nb_dofs * q;
        double *t_n_ptr = &faceN[family][3 * qd_shift];
        double *t_diff_n_ptr = &diff_faceN[family][6 * qd_shift];
        auto t_n = getFTensor1FromPtr<3>(t_n_ptr);
        auto t_diff_n = getFTensor2FromPtr<3, 2>(t_diff_n_ptr);

        for (int n = 0; n != nb_dofs; n++) {
          int i = permute[family][3 * n + 0];
          int j = permute[family][3 * n + 1];

          const double z = zeta[j + 2];
          const double e = eta[i];
          const double a = z * e;
          const double d_a[] = {

              diff_zeta[0 * (pp + 2) + j + 2] * e + z * diff_eta[0 * qq + i],

              diff_zeta[1 * (pp + 2) + j + 2] * e + z * diff_eta[1 * qq + i]};

          for (int d = 0; d != 2; ++d) {
            t_n(d) = a * diff_eta_ksi[family][d];
            for (int m = 0; m != 2; ++m) {
              t_diff_n(d, m) = d_a[m] * diff_eta_ksi[family][d];
            }
          }
          t_n(2) = 0;
          t_diff_n(2, 0) = 0;
          t_diff_n(2, 1) = 0;
          ++t_n;
          ++t_diff_n;
        }
      }
    }
  }
  MoFEMFunctionReturnHot(0);
}

MoFEMErrorCode MoFEM::DemkowiczHexAndQuad::Hdiv_FaceShapeFunctions_ONQUAD(
    int *face_nodes, int *p, double *N, double *diffN, double *faceN,
    double *diff_faceN, int nb_integration_pts) {
  MoFEMFunctionBeginHot;

  const int nb_dofs = (p[0] * p[1]);

  if (nb_dofs > 0) {

    int permute[nb_dofs][3];
    CHKERR ::DemkowiczHexAndQuad::monom_ordering(&permute[0][0], p[0] - 1,
                                                 p[1] - 1);

    const int n0 = face_nodes[0];
    const int n1 = face_nodes[1];
    const int n2 = face_nodes[2];
    const int n3 = face_nodes[3];

    FTensor::Index<'i', 3> i;
    FTensor::Index<'j', 3> j;
    FTensor::Index<'k', 3> k;

    auto t_n = getFTensor1FromPtr<3>(faceN);
    auto t_diff_n = getFTensor2FromPtr<3, 2>(diff_faceN);

    for (int q = 0; q != nb_integration_pts; q++) {

      const int shift = 4 * q;
      const double shape0 = N[shift + n0];
      const double shape1 = N[shift + n1];
      const double shape2 = N[shift + n2];
      const double shape3 = N[shift + n3];
      const double ksi01 = (shape1 + shape2 - shape0 - shape3);
      const double ksi12 = (shape2 + shape3 - shape1 - shape0);

      const int diff_shift = 2 * shift;
      FTensor::Tensor1<double, 3> t_diff_ksi01;
      FTensor::Tensor1<double, 3> t_diff_ksi12;
      for (int d = 0; d != 2; d++) {
        const double diff_shape0 = diffN[diff_shift + 2 * n0 + d];
        const double diff_shape1 = diffN[diff_shift + 2 * n1 + d];
        const double diff_shape2 = diffN[diff_shift + 2 * n2 + d];
        const double diff_shape3 = diffN[diff_shift + 2 * n3 + d];
        t_diff_ksi01(d + 1) =
            (diff_shape1 + diff_shape2 - diff_shape0 - diff_shape3);
        t_diff_ksi12(d + 1) =
            (diff_shape2 + diff_shape3 - diff_shape1 - diff_shape0);
      }
      t_diff_ksi01(0) = 0;
      t_diff_ksi12(0) = 0;

      FTensor::Tensor1<double, 3> t_cross;
      t_cross(i) =
          FTensor::levi_civita(i, j, k) * t_diff_ksi01(j) * t_diff_ksi12(k);

      double zeta[p[0] + 1];
      double diff_zeta[2 * (p[0] + 1)];
      CHKERR Legendre_polynomials(p[0], ksi01, &t_diff_ksi01(0), zeta,
                                  diff_zeta, 2);

      double eta[p[1] + 1];
      double diff_eta[2 * (p[1] + 1)];
      CHKERR Legendre_polynomials(p[1], ksi12, &t_diff_ksi12(0), eta, diff_eta,
                                  2);

      for (int n = 0; n != nb_dofs; ++n) {
        int ii = permute[n][0];
        int jj = permute[n][1];

        const double z = zeta[ii];
        const double e = eta[jj];
        const double ez = e * z;

        auto t_diff_zeta = FTensor::Tensor1<double, 2>{
            diff_zeta[ii], diff_zeta[(p[0] + 1) + ii]};
        auto t_diff_eta = FTensor::Tensor1<double, 2>{
            diff_eta[jj], diff_eta[(p[1] + 1) + jj]};

        FTensor::Index<'J', 2> J;
        t_n(i) = t_cross(i) * ez;
        t_diff_n(i, J) = t_cross(i) * (t_diff_zeta(J) * e + t_diff_eta(J) * z);

        ++t_n;
        ++t_diff_n;
      }
    }
  }
  MoFEMFunctionReturnHot(0);
}

/* Reference Hex and its cannonical vertex and edge numbering
               7 ---------10--------- 6
              /|                    /|
             / |                   / |
           11  |                  9  |             x3
           /   7                 /   |             |
          /    |                /    6             |    x2
        4 ----------8--------- 5     |             |   /
         |     |               |     |             |  /
         |    3 ----------2---------- 2            | /
         4    /                |    /              o -------- x1
         |   /                 5   /
         |  3                  |  1
         | /                   | /
         |/                    |/
        0 ---------0---------- 1

  Hex Face Cannonical numbering

        1. 0 1 2 3
        2. 0 1 5 4
        3. 1 2 6 5
        4. 3 2 6 7
        5. 0 3 7 4
        6. 4 5 6 7
*/

MoFEMErrorCode MoFEM::DemkowiczHexAndQuad::H1_EdgeShapeFunctions_ONHEX(
    int *sense, int *p, double *N, double *diff_N, double *edgeN[12],
    double *diff_edgeN[12], int nb_integration_pts) {
  MoFEMFunctionBeginHot;

  EntityType sub_type;
  int num_sub_ent_vertices;
  const short int *edges_conn[12];
  for (int e = 0; e != 12; ++e)
    edges_conn[e] =
        CN::SubEntityVertexIndices(MBHEX, 1, e, sub_type, num_sub_ent_vertices);

  const short int *face_conn[6];
  for (int f = 0; f != 6; ++f)
    face_conn[f] =
        CN::SubEntityVertexIndices(MBHEX, 2, f, sub_type, num_sub_ent_vertices);

  constexpr int quad_edge[12][2] = {{3, 1}, {0, 2}, {1, 3}, {2, 0},
                                    {4, 5}, {4, 5}, {4, 5}, {4, 5},
                                    {3, 1}, {0, 2}, {1, 3}, {2, 0}};

  for (int qq = 0; qq != nb_integration_pts; ++qq) {

    const int shift = 8 * qq;

    double ksi[12];
    double diff_ksi[12][3];
    for (int e = 0; e != 12; ++e) {

      ksi[e] = 0;
      for (int d = 0; d != 3; ++d)
        diff_ksi[e][d] = 0;
      for (int n = 0; n != 4; ++n) {
        const auto n1 = shift + face_conn[quad_edge[e][1]][n];
        const auto n0 = shift + face_conn[quad_edge[e][0]][n];
        ksi[e] += N[n1] - N[n0];
        const auto n03 = 3 * n0;
        const auto n13 = 3 * n1;
        for (int d = 0; d != 3; ++d)
          diff_ksi[e][d] += diff_N[n13 + d] - diff_N[n03 + d];
      }

      ksi[e] *= sense[e];
      for (int d = 0; d != 3; ++d)
        diff_ksi[e][d] *= sense[e];
    }

    double mu[12];
    double diff_mu[12][3];
    for (int e = 0; e != 12; ++e) {
      const auto n0 = shift + edges_conn[e][0];
      const auto n1 = shift + edges_conn[e][1];
      mu[e] = N[n0] + N[n1];
      const auto n03 = 3 * n0;
      const auto n13 = 3 * n1;
      for (int d = 0; d != 3; ++d) {
        diff_mu[e][d] = diff_N[n03 + d] + diff_N[n13 + d];
      }
    }

    for (int e = 0; e != 12; e++) {

      const int nb_dofs = (p[e] - 1);
      if (nb_dofs > 0) {

        double L[p[e] + 2];
        double diffL[3 * (p[e] + 2)];
        CHKERR Lobatto_polynomials(p[e] + 1, ksi[e], diff_ksi[e], L, diffL, 3);

        const int qd_shift = nb_dofs * qq;
        for (int n = 0; n != nb_dofs; n++) {
          edgeN[e][qd_shift + n] = mu[e] * L[n + 2];
          for (int d = 0; d != 3; ++d) {
            diff_edgeN[e][3 * (qd_shift + n) + d] =

                diff_mu[e][d] * L[n + 2]

                +

                mu[e] * diffL[d * (p[e] + 2) + n + 2];
          }
        }
      }
    }
  }
  MoFEMFunctionReturnHot(0);
}

MoFEMErrorCode MoFEM::DemkowiczHexAndQuad::H1_FaceShapeFunctions_ONHEX(
    int *face_nodes, int *face_nodes_order, int *p, double *N, double *diffN,
    double *faceN[6], double *diff_faceN[6], int nb_integration_pts) {
  MoFEMFunctionBeginHot;

  constexpr int opposite_face_node[6][4] = {

      {3, 2, 6, 7},

      {0, 3, 7, 4},

      {1, 0, 4, 5},

      {2, 1, 5, 6},

      {4, 7, 6, 5},

      {0, 1, 2, 3}

  };

  for (int face = 0; face != 6; face++) {

    if (p[face] > 1) {
      const int nb_dofs = (p[face] - 1) * (p[face] - 1);

      const int n0 = face_nodes[4 * face + 0];
      const int n1 = face_nodes[4 * face + 1];
      const int n2 = face_nodes[4 * face + 2];
      const int n3 = face_nodes[4 * face + 3];

      const int o0 = opposite_face_node[face][face_nodes_order[4 * face + 0]];
      const int o1 = opposite_face_node[face][face_nodes_order[4 * face + 1]];
      const int o2 = opposite_face_node[face][face_nodes_order[4 * face + 2]];
      const int o3 = opposite_face_node[face][face_nodes_order[4 * face + 3]];

      int permute[nb_dofs][3];
      CHKERR ::DemkowiczHexAndQuad::monom_ordering(&permute[0][0], p[face] - 2,
                                                   p[face] - 2);

      for (int qq = 0; qq != nb_integration_pts; qq++) {

        const int shift = 8 * qq;

        const double shape0 = N[shift + n0];
        const double shape1 = N[shift + n1];
        const double shape2 = N[shift + n2];
        const double shape3 = N[shift + n3];

        const double o_shape0 = N[shift + o0];
        const double o_shape1 = N[shift + o1];
        const double o_shape2 = N[shift + o2];
        const double o_shape3 = N[shift + o3];

        const double ksi01 = (shape1 + shape2 - shape0 - shape3) +
                             (o_shape1 + o_shape2 - o_shape0 - o_shape3);
        const double ksi12 = (shape2 + shape3 - shape1 - shape0) +
                             (o_shape2 + o_shape3 - o_shape1 - o_shape0);
        const double mu = shape1 + shape2 + shape0 + shape3;

        const int diff_shift = 3 * shift;
        double diff_ksi01[3], diff_ksi12[3], diff_mu[3];
        for (int d = 0; d != 3; d++) {
          const double diff_shape0 = diffN[diff_shift + 3 * n0 + d];
          const double diff_shape1 = diffN[diff_shift + 3 * n1 + d];
          const double diff_shape2 = diffN[diff_shift + 3 * n2 + d];
          const double diff_shape3 = diffN[diff_shift + 3 * n3 + d];
          const double o_diff_shape0 = diffN[diff_shift + 3 * o0 + d];
          const double o_diff_shape1 = diffN[diff_shift + 3 * o1 + d];
          const double o_diff_shape2 = diffN[diff_shift + 3 * o2 + d];
          const double o_diff_shape3 = diffN[diff_shift + 3 * o3 + d];

          diff_ksi01[d] =
              (diff_shape1 + diff_shape2 - diff_shape0 - diff_shape3) +
              (o_diff_shape1 + o_diff_shape2 - o_diff_shape0 - o_diff_shape3);
          diff_ksi12[d] =
              (diff_shape2 + diff_shape3 - diff_shape1 - diff_shape0) +
              (o_diff_shape2 + o_diff_shape3 - o_diff_shape1 - o_diff_shape0);
          diff_mu[d] = (diff_shape1 + diff_shape2 + diff_shape0 + diff_shape3);
        }

        double L01[p[face] + 2];
        double diffL01[3 * (p[face] + 2)];
        CHKERR Lobatto_polynomials(p[face] + 1, ksi01, diff_ksi01, L01, diffL01,
                                   3);
        double L12[p[face] + 2];
        double diffL12[3 * (p[face] + 2)];
        CHKERR Lobatto_polynomials(p[face] + 1, ksi12, diff_ksi12, L12, diffL12,
                                   3);

        int qd_shift = nb_dofs * qq;
        for (int n = 0; n != nb_dofs; ++n) {
          const int s1 = permute[n][0];
          const int s2 = permute[n][1];
          const double vol = L01[s1 + 2] * L12[s2 + 2];
          faceN[face][qd_shift + n] = vol * mu;
          for (int d = 0; d != 3; ++d) {
            diff_faceN[face][3 * (qd_shift + n) + d] =
                (diffL01[d * (p[face] + 2) + s1 + 2] * L12[s2 + 2] +
                 diffL12[d * (p[face] + 2) + s2 + 2] * L01[s1 + 2]) *
                    mu +
                vol * diff_mu[d];
          }
        }
      }
    }
  }
  MoFEMFunctionReturnHot(0);
}

MoFEMErrorCode MoFEM::DemkowiczHexAndQuad::H1_InteriorShapeFunctions_ONHEX(
    const int *p, double *N, double *N_diff, double *faceN, double *diff_faceN,
    int nb_integration_pts) {
  MoFEMFunctionBeginHot;

  const int nb_bases = (p[0] - 1) * (p[1] - 1) * (p[2] - 1);

  if (nb_bases > 0) {

    int permute[nb_bases][3];
    CHKERR ::DemkowiczHexAndQuad::monom_ordering(&permute[0][0], p[0] - 2,
                                                 p[1] - 2, p[2] - 2);

    double P0[p[0] + 2];
    double diffL0[3 * (p[0] + 2)];
    double P1[p[1] + 2];
    double diffL1[3 * (p[1] + 2)];
    double P2[p[2] + 2];
    double diffL2[3 * (p[2] + 2)];

    for (int qq = 0; qq != nb_integration_pts; ++qq) {

      const int shift = 8 * qq;
      double ksi[3] = {0, 0, 0};
      double diff_ksi[3][3] = {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}};
      ::DemkowiczHexAndQuad::get_ksi_hex(shift, N, N_diff, ksi, diff_ksi);

      double L0[p[0] + 2];
      double diffL0[3 * (p[0] + 2)];
      double L1[p[1] + 2];
      double diffL1[3 * (p[1] + 2)];
      double L2[p[2] + 2];
      double diffL2[3 * (p[2] + 2)];

      CHKERR Lobatto_polynomials(p[0] + 1, ksi[0], diff_ksi[0], L0, diffL0, 3);
      CHKERR Lobatto_polynomials(p[1] + 1, ksi[1], diff_ksi[1], L1, diffL1, 3);
      CHKERR Lobatto_polynomials(p[2] + 1, ksi[2], diff_ksi[2], L2, diffL2, 3);

      const int qd_shift = nb_bases * qq;
      for (int n = 0; n != nb_bases; ++n) {
        const int s1 = permute[n][0];
        const int s2 = permute[n][1];
        const int s3 = permute[n][2];

        const double l0l1 = L0[s1 + 2] * L1[s2 + 2];
        const double l0l2 = L0[s1 + 2] * L2[s3 + 2];
        const double l1l2 = L1[s2 + 2] * L2[s3 + 2];

        faceN[qd_shift + n] = l0l1 * L2[s3 + 2];
        for (int d = 0; d != 3; ++d) {
          diff_faceN[3 * (qd_shift + n) + d] =
              (diffL0[d * (p[0] + 2) + s1 + 2] * l1l2 +
               diffL1[d * (p[1] + 2) + s2 + 2] * l0l2 +
               diffL2[d * (p[2] + 2) + s3 + 2] * l0l1);
        }
      }
    }
  }
  MoFEMFunctionReturnHot(0);
}

MoFEMErrorCode MoFEM::DemkowiczHexAndQuad::L2_InteriorShapeFunctions_ONHEX(
    const int *p, double *N, double *N_diff, double *volN, double *diff_volN,
    int nb_integration_pts) {
  MoFEMFunctionBeginHot;

  const int nb_bases = (p[0] + 1) * (p[1] + 1) * (p[2] + 1);
  if (nb_bases > 0) {

    int permute[nb_bases][3];
    CHKERR ::DemkowiczHexAndQuad::monom_ordering(&permute[0][0], p[0], p[1],
                                                 p[2]);

    double P0[p[0] + 2];
    double diffL0[3 * (p[0] + 2)];
    double P1[p[1] + 2];
    double diffL1[3 * (p[1] + 2)];
    double P2[p[2] + 2];
    double diffL2[3 * (p[2] + 2)];

    for (int qq = 0; qq != nb_integration_pts; qq++) {

      const int shift = 8 * qq;
      double ksi[3] = {0, 0, 0};
      double diff_ksi[3][3] = {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}};
      ::DemkowiczHexAndQuad::get_ksi_hex(shift, N, N_diff, ksi, diff_ksi);

      CHKERR Legendre_polynomials(p[0] + 1, ksi[0], diff_ksi[0], P0, diffL0, 3);
      CHKERR Legendre_polynomials(p[1] + 1, ksi[1], diff_ksi[1], P1, diffL1, 3);
      CHKERR Legendre_polynomials(p[2] + 1, ksi[2], diff_ksi[2], P2, diffL2, 3);

      const int qd_shift = qq * nb_bases;
      for (int n = 0; n != nb_bases; ++n) {
        const int ii = permute[n][0];
        const int jj = permute[n][1];
        const int kk = permute[n][2];

        const double p1p2 = P1[jj] * P2[kk];
        const double p0p1 = P0[ii] * P1[jj];
        const double p0p2 = P0[ii] * P2[kk];

        volN[qd_shift + n] = p0p1 * P2[kk];
        for (int d = 0; d != 3; ++d) {
          diff_volN[3 * (qd_shift + n) + d] =
              p1p2 * diffL0[d * (p[0] + 2) + ii] +
              p0p2 * diffL1[d * (p[1] + 2) + jj] +
              p0p1 * diffL2[d * (p[2] + 2) + kk];
        }
      }
    }
  }
  MoFEMFunctionReturnHot(0);
}

MoFEMErrorCode MoFEM::DemkowiczHexAndQuad::Hcurl_EdgeShapeFunctions_ONHEX(
    int *sense, int *p, double *N, double *diff_N, double *edgeN[12],
    double *diff_edgeN[12], int nb_integration_pts) {

  FTensor::Index<'i', 3> i;
  FTensor::Index<'j', 3> j;

  MoFEMFunctionBeginHot;

  EntityType sub_type;
  int num_sub_ent_vertices;
  const short int *edges_conn[12];
  for (int e = 0; e != 12; ++e)
    edges_conn[e] =
        CN::SubEntityVertexIndices(MBHEX, 1, e, sub_type, num_sub_ent_vertices);

  const short int *face_conn[6];
  for (int f = 0; f != 6; ++f)
    face_conn[f] =
        CN::SubEntityVertexIndices(MBHEX, 2, f, sub_type, num_sub_ent_vertices);

  constexpr int quad_edge[12][2] = {{3, 1}, {0, 2}, {1, 3}, {2, 0},
                                    {4, 5}, {4, 5}, {4, 5}, {4, 5},
                                    {3, 1}, {0, 2}, {1, 3}, {2, 0}};

  for (int qq = 0; qq != nb_integration_pts; qq++) {

    double ksi[12];
    double diff_ksi[12 * 3];

    const int shift = 8 * qq;

    auto calc_ksi = [&]() {
      auto t_diff_ksi = getFTensor1FromPtr<3>(diff_ksi);
      for (int e = 0; e != 12; ++e) {
        if (p[e] > 0) {
          ksi[e] = 0;
          t_diff_ksi(i) = 0;
          for (int n = 0; n != 4; ++n) {
            const auto n1 = shift + face_conn[quad_edge[e][1]][n];
            const auto n0 = shift + face_conn[quad_edge[e][0]][n];
            ksi[e] += N[n1] - N[n0];
            const auto n03 = 3 * n0;
            const auto n13 = 3 * n1;
            for (int d = 0; d != 3; ++d)
              t_diff_ksi(d) += diff_N[n13 + d] - diff_N[n03 + d];
          }

          ksi[e] *= sense[e];
          t_diff_ksi(i) *= sense[e];

          ++t_diff_ksi;
        }
      }
    };

    double mu[12];
    double diff_mu[12][3];
    auto calc_mu = [&]() {
      for (int e = 0; e != 12; ++e) {
        if (p[e] > 0) {
          const auto n0 = shift + edges_conn[e][0];
          const auto n1 = shift + edges_conn[e][1];
          mu[e] = N[n0] + N[n1];
          const auto n03 = 3 * n0;
          const auto n13 = 3 * n1;
          for (int d = 0; d != 3; ++d) {
            diff_mu[e][d] = diff_N[n03 + d] + diff_N[n13 + d];
          }
        }
      }
    };

    auto calc_base = [&]() {
      auto t_diff_ksi = getFTensor1FromPtr<3>(diff_ksi);
      for (int ee = 0; ee != 12; ee++) {
        if (p[ee] > 0) {

          double L[p[ee]];
          double diffL[3 * (p[ee])];
          CHKERR Legendre_polynomials(p[ee] - 1, ksi[ee], &(t_diff_ksi(0)), L,
                                      diffL, 3);

          int qd_shift = p[ee] * qq;
          auto t_n = getFTensor1FromPtr<3>(&edgeN[ee][3 * qd_shift]);
          auto t_diff_n =
              getFTensor2FromPtr<3, 3>(&diff_edgeN[ee][9 * qd_shift]);

          for (int ii = 0; ii != p[ee]; ii++) {

            const double a = mu[ee] * L[ii];
            auto t_diff_a = FTensor::Tensor1<const double, 3>{

                diff_mu[ee][0] * L[ii] + diffL[0 * p[ee] + ii],

                diff_mu[ee][1] * L[ii] + diffL[1 * p[ee] + ii],

                diff_mu[ee][2] * L[ii] + diffL[2 * p[ee] + ii]

            };

            t_n(i) = (a / 2) * t_diff_ksi(i);
            t_diff_n(i, j) = (t_diff_a(j) / 2) * t_diff_ksi(i);

            ++t_n;
            ++t_diff_n;
          }

          ++t_diff_ksi;
        }
      }
    };

    calc_ksi();
    calc_mu();
    calc_base();
  }
  MoFEMFunctionReturnHot(0);
}

MoFEMErrorCode MoFEM::DemkowiczHexAndQuad::Hcurl_FaceShapeFunctions_ONHEX(
    int *face_nodes, int *face_nodes_order, int *p, double *N, double *diffN,
    double *faceN[6][2], double *diff_faceN[6][2], int nb_integration_pts) {
  MoFEMFunctionBeginHot;

  constexpr int opposite_face_node[6][4] = {

      {3, 2, 6, 7},

      {0, 3, 7, 4},

      {1, 0, 4, 5},

      {2, 1, 5, 6},

      {4, 7, 6, 5},

      {0, 1, 2, 3}

  };

  for (int face = 0; face != 6; face++) {
    if ((p[face] - 1) * p[face] > 0) {

      const int n0 = face_nodes[4 * face + 0];
      const int n1 = face_nodes[4 * face + 1];
      const int n2 = face_nodes[4 * face + 2];
      const int n3 = face_nodes[4 * face + 3];

      const int o0 = opposite_face_node[face][face_nodes_order[4 * face + 0]];
      const int o1 = opposite_face_node[face][face_nodes_order[4 * face + 1]];
      const int o2 = opposite_face_node[face][face_nodes_order[4 * face + 2]];
      const int o3 = opposite_face_node[face][face_nodes_order[4 * face + 3]];

      int permute[(p[face] - 1) * p[face]][3];
      CHKERR ::DemkowiczHexAndQuad::monom_ordering(&permute[0][0], p[face] - 1,
                                                   p[face] - 2);

      for (int q = 0; q != nb_integration_pts; ++q) {

        const int shift = 8 * q;

        const double shape0 = N[shift + n0];
        const double shape1 = N[shift + n1];
        const double shape2 = N[shift + n2];
        const double shape3 = N[shift + n3];

        const double o_shape0 = N[shift + o0];
        const double o_shape1 = N[shift + o1];
        const double o_shape2 = N[shift + o2];
        const double o_shape3 = N[shift + o3];

        const double ksi01 = (shape1 + shape2 - shape0 - shape3) +
                             (o_shape1 + o_shape2 - o_shape0 - o_shape3);
        const double ksi12 = (shape2 + shape3 - shape1 - shape0) +
                             (o_shape2 + o_shape3 - o_shape1 - o_shape0);
        const double mu = shape1 + shape2 + shape0 + shape3;

        const int diff_shift = 3 * shift;
        double diff_ksi01[3], diff_ksi12[3], diff_mu[3];
        for (int d = 0; d != 3; ++d) {
          const double diff_shape0 = diffN[diff_shift + 3 * n0 + d];
          const double diff_shape1 = diffN[diff_shift + 3 * n1 + d];
          const double diff_shape2 = diffN[diff_shift + 3 * n2 + d];
          const double diff_shape3 = diffN[diff_shift + 3 * n3 + d];
          const double o_diff_shape0 = diffN[diff_shift + 3 * o0 + d];
          const double o_diff_shape1 = diffN[diff_shift + 3 * o1 + d];
          const double o_diff_shape2 = diffN[diff_shift + 3 * o2 + d];
          const double o_diff_shape3 = diffN[diff_shift + 3 * o3 + d];

          diff_ksi01[d] =
              (diff_shape1 + diff_shape2 - diff_shape0 - diff_shape3) +
              (o_diff_shape1 + o_diff_shape2 - o_diff_shape0 - o_diff_shape3);
          diff_ksi12[d] =
              (diff_shape2 + diff_shape3 - diff_shape1 - diff_shape0) +
              (o_diff_shape2 + o_diff_shape3 - o_diff_shape1 - o_diff_shape0);
          diff_mu[d] = (diff_shape1 + diff_shape2 + diff_shape0 + diff_shape3);
        }

        int pq[2] = {p[face], p[face]};
        int qp[2] = {p[face], p[face]};
        double mu_ksi_eta[2] = {ksi01, ksi12};
        double mu_eta_ksi[2] = {ksi12, ksi01};
        double *diff_ksi_eta[2] = {diff_ksi01, diff_ksi12};
        double *diff_eta_ksi[2] = {diff_ksi12, diff_ksi01};

        for (int family = 0; family != 2; ++family) {

          const int pp = pq[family];
          const int qq = qp[family];
          const int nb_dofs = (pp - 1) * qq;

          if (nb_dofs > 0) {

            double zeta[pp + 2];
            double diff_zeta[3 * (pp + 2)];
            CHKERR Lobatto_polynomials(pp + 1, mu_ksi_eta[family],
                                       diff_ksi_eta[family], zeta, diff_zeta,
                                       3);

            double eta[qq];
            double diff_eta[3 * qq];
            CHKERR Legendre_polynomials(qq - 1, mu_eta_ksi[family],
                                        diff_eta_ksi[family], eta, diff_eta, 3);

            const int qd_shift = nb_dofs * q;
            double *t_n_ptr = &(faceN[face][family][3 * qd_shift]);
            double *t_diff_n_ptr = &(diff_faceN[face][family][9 * qd_shift]);
            auto t_n = getFTensor1FromPtr<3>(t_n_ptr);
            auto t_diff_n = getFTensor2FromPtr<3, 3>(t_diff_n_ptr);

            for (int n = 0; n != nb_dofs; n++) {
              int i = permute[n][0];
              int j = permute[n][1];

              const double z = zeta[j + 2];
              const double e = eta[i];
              const double ze = z * e;
              const double a = ze * mu;
              double d_a[3];

              for (int d = 0; d != 3; ++d)
                d_a[d] = (diff_zeta[d * (pp + 2) + j + 2] * e +
                          z * diff_eta[d * qq + i]) *
                             mu +
                         ze * diff_mu[d];

              for (int d = 0; d != 3; ++d) {
                t_n(d) = a * diff_eta_ksi[family][d];
                for (int m = 0; m != 3; ++m) {
                  t_diff_n(d, m) = d_a[m] * diff_eta_ksi[family][d];
                }
              }

              ++t_n;
              ++t_diff_n;
            }
          }
        }
      }
    }
  }

  MoFEMFunctionReturnHot(0);
}

MoFEMErrorCode MoFEM::DemkowiczHexAndQuad::Hcurl_InteriorShapeFunctions_ONHEX(
    int *p, double *N, double *diffN, double *volN[], double *diff_volN[],
    int nb_integration_pts) {
  MoFEMFunctionBeginHot;

  int pqr[3] = {p[0], p[1], p[2]};
  int qrp[3] = {p[1], p[2], p[0]};
  int rpq[3] = {p[2], p[0], p[1]};

  const int nb_dofs_fm0 = (p[0] - 1) * p[1] * p[2];
  const int nb_dofs_fm1 = p[0] * (p[1] - 1) * p[2];
  const int nb_dofs_fm2 = p[0] * p[1] * (p[2] - 1);
  int permute_fm0[3 * nb_dofs_fm0];
  int permute_fm1[3 * nb_dofs_fm1];
  int permute_fm2[3 * nb_dofs_fm2];

  std::array<int *, 3> permute = {&permute_fm0[0], &permute_fm1[0],
                                  &permute_fm2[0]};

  for (int fam = 0; fam != 3; ++fam) {
    const int qqq = qrp[fam];
    const int rrr = rpq[fam];
    const int ppp = pqr[fam];
    CHKERR ::DemkowiczHexAndQuad::monom_ordering(permute[fam], ppp - 2, qqq - 1,
                                                 rrr - 2);
  }

  for (int qq = 0; qq != nb_integration_pts; ++qq) {

    const int shift = 8 * qq;
    double ksi[3] = {0, 0, 0};
    double diff_ksi[3][3] = {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}};
    ::DemkowiczHexAndQuad::get_ksi_hex(shift, N, diffN, ksi, diff_ksi);

    double ksi_eta_gma[3] = {ksi[0], ksi[1], ksi[2]};
    double eta_gma_ksi[3] = {ksi[1], ksi[2], ksi[0]};
    double gma_ksi_eta[3] = {ksi[2], ksi[0], ksi[1]};

    double *diff_ksi_eta_gma[3] = {diff_ksi[0], diff_ksi[1], diff_ksi[2]};
    double *diff_eta_gma_ksi[3] = {diff_ksi[1], diff_ksi[2], diff_ksi[0]};
    double *diff_gma_ksi_eta[3] = {diff_ksi[2], diff_ksi[0], diff_ksi[1]};

    int pqr[3] = {p[0], p[1], p[2]};
    int qrp[3] = {p[1], p[2], p[0]};
    int rpq[3] = {p[2], p[0], p[1]};
    for (int fam = 0; fam != 3; ++fam) {

      const int qqq = qrp[fam];
      const int rrr = rpq[fam];
      const int ppp = pqr[fam];

      const int nb_dofs = (ppp - 1) * qqq * (rrr - 1);
      if (nb_dofs > 0) {

        double phi_j[ppp + 2];
        double diff_phi_j[3 * (ppp + 2)];

        CHKERR Lobatto_polynomials(ppp + 1, ksi_eta_gma[fam],
                                   diff_ksi_eta_gma[fam], phi_j, diff_phi_j, 3);

        double eta_i[qqq];
        double diff_eta_i[3 * qqq];

        CHKERR Legendre_polynomials(qqq - 1, eta_gma_ksi[fam],
                                    diff_eta_gma_ksi[fam], eta_i, diff_eta_i,
                                    3);

        double phi_k[rrr + 2];
        double diff_phi_k[3 * (rrr + 2)];

        CHKERR Lobatto_polynomials(rrr + 1, gma_ksi_eta[fam],
                                   diff_gma_ksi_eta[fam], phi_k, diff_phi_k, 3);

        int qd_shift = nb_dofs * qq;
        double *t_n_ptr = &volN[fam][3 * qd_shift];
        double *t_diff_n_ptr = &diff_volN[fam][3 * 3 * qd_shift];
        auto t_n = getFTensor1FromPtr<3>(t_n_ptr);
        auto t_diff_n = getFTensor2FromPtr<3, 3>(t_diff_n_ptr);

        int n = 0;
        for (; n != nb_dofs; n++) {
          int ii = permute[fam][3 * n + 0];
          int jj = permute[fam][3 * n + 1];
          int kk = permute[fam][3 * n + 2];

          const double p_k = phi_k[kk + 2];
          const double p_j = phi_j[jj + 2];
          const double e = eta_i[ii];
          const double a = p_j * p_k * e;

          const double d_a[] = {
              diff_phi_k[0 * (ppp + 2) + kk + 2] * p_j * e +
                  p_k * diff_phi_j[0 * (rrr + 2) + jj + 2] * e +
                  p_k * p_j * diff_eta_i[0 * qqq + ii],

              diff_phi_k[1 * (ppp + 2) + kk + 2] * p_j * e +
                  p_k * diff_phi_j[1 * (rrr + 2) + jj + 2] * e +
                  p_k * p_j * diff_eta_i[1 * qqq + ii],

              diff_phi_k[2 * (ppp + 2) + kk + 2] * p_j * e +
                  p_k * diff_phi_j[2 * (rrr + 2) + jj + 2] * e +
                  p_k * p_j * diff_eta_i[2 * qqq + ii]};

          for (int d = 0; d != 3; ++d) {
            t_n(d) = a * diff_eta_gma_ksi[fam][d];
            for (int m = 0; m != 3; ++m) {
              t_diff_n(d, m) = d_a[m] * diff_eta_gma_ksi[fam][d];
            }
          }
          ++t_n;
          ++t_diff_n;
        }
      }
    }
  }

  MoFEMFunctionReturnHot(0);
}

MoFEMErrorCode MoFEM::DemkowiczHexAndQuad::Hdiv_FaceShapeFunctions_ONHEX(
    int *face_nodes, int *face_nodes_order, int *p, double *N, double *diffN,
    double *faceN[6], double *div_faceN[6], int nb_integration_pts) {
  MoFEMFunctionBeginHot;

  FTensor::Index<'i', 3> i;
  FTensor::Index<'j', 3> j;
  FTensor::Index<'k', 3> k;

  constexpr int opposite_face_node[6][4] = {

      {3, 2, 6, 7},

      {0, 3, 7, 4},

      {1, 0, 4, 5},

      {2, 1, 5, 6},

      {4, 7, 6, 5},

      {0, 1, 2, 3}

  };

  for (int face = 0; face != 6; face++) {

    const int nb_dofs = (p[face] * p[face]);

    if (nb_dofs > 0) {

      auto t_n = getFTensor1FromPtr<3>(faceN[face]);
      auto t_diff_n = getFTensor2FromPtr<3, 3>(div_faceN[face]);

      const int n0 = face_nodes[4 * face + 0];
      const int n1 = face_nodes[4 * face + 1];
      const int n2 = face_nodes[4 * face + 2];
      const int n3 = face_nodes[4 * face + 3];

      const int o0 = opposite_face_node[face][face_nodes_order[4 * face + 0]];
      const int o1 = opposite_face_node[face][face_nodes_order[4 * face + 1]];
      const int o2 = opposite_face_node[face][face_nodes_order[4 * face + 2]];
      const int o3 = opposite_face_node[face][face_nodes_order[4 * face + 3]];

      int permute[nb_dofs][3];
      CHKERR ::DemkowiczHexAndQuad::monom_ordering(&permute[0][0], p[face] - 1,
                                                   p[face] - 1);

      for (int q = 0; q != nb_integration_pts; ++q) {

        const int shift = 8 * q;

        const double shape0 = N[shift + n0];
        const double shape1 = N[shift + n1];
        const double shape2 = N[shift + n2];
        const double shape3 = N[shift + n3];

        const double o_shape0 = N[shift + o0];
        const double o_shape1 = N[shift + o1];
        const double o_shape2 = N[shift + o2];
        const double o_shape3 = N[shift + o3];

        const double ksi01 = (shape1 + shape2 - shape0 - shape3) +
                             (o_shape1 + o_shape2 - o_shape0 - o_shape3);
        const double ksi12 = (shape2 + shape3 - shape1 - shape0) +
                             (o_shape2 + o_shape3 - o_shape1 - o_shape0);
        const double mu = shape1 + shape2 + shape0 + shape3;

        const int diff_shift = 3 * shift;
        FTensor::Tensor1<double, 3> t_diff_ksi01;
        FTensor::Tensor1<double, 3> t_diff_ksi12;
        FTensor::Tensor1<double, 3> t_diff_mu;

        for (int d = 0; d != 3; ++d) {
          const double diff_shape0 = diffN[diff_shift + 3 * n0 + d];
          const double diff_shape1 = diffN[diff_shift + 3 * n1 + d];
          const double diff_shape2 = diffN[diff_shift + 3 * n2 + d];
          const double diff_shape3 = diffN[diff_shift + 3 * n3 + d];
          const double o_diff_shape0 = diffN[diff_shift + 3 * o0 + d];
          const double o_diff_shape1 = diffN[diff_shift + 3 * o1 + d];
          const double o_diff_shape2 = diffN[diff_shift + 3 * o2 + d];
          const double o_diff_shape3 = diffN[diff_shift + 3 * o3 + d];
          t_diff_ksi01(d) =
              (diff_shape1 + diff_shape2 - diff_shape0 - diff_shape3) +
              (o_diff_shape1 + o_diff_shape2 - o_diff_shape0 - o_diff_shape3);
          t_diff_ksi12(d) =
              (diff_shape2 + diff_shape3 - diff_shape1 - diff_shape0) +
              (o_diff_shape2 + o_diff_shape3 - o_diff_shape1 - o_diff_shape0);
          t_diff_mu(d) =
              (diff_shape1 + diff_shape2 + diff_shape0 + diff_shape3);
        }

        FTensor::Tensor1<double, 3> t_cross;
        t_cross(i) =
            FTensor::levi_civita(i, j, k) * t_diff_ksi01(j) * t_diff_ksi12(k);

        double zeta[p[0] + 1];
        double diff_zeta[3 * (p[0] + 1)];
        CHKERR Legendre_polynomials(p[0], ksi01, &t_diff_ksi01(0), zeta,
                                    diff_zeta, 3);

        double eta[p[1] + 1];
        double diff_eta[3 * (p[1] + 1)];
        CHKERR Legendre_polynomials(p[1], ksi12, &t_diff_ksi12(0), eta,
                                    diff_eta, 3);

        for (int n = 0; n != nb_dofs; ++n) {
          int ii = permute[n][0];
          int jj = permute[n][1];

          const double z = zeta[ii];
          const double e = eta[jj];
          const double ez = e * z;

          auto t_diff_zeta = FTensor::Tensor1<double, 3>{
              diff_zeta[ii], diff_zeta[1 * (p[0] + 1) + ii],
              diff_zeta[2 * (p[0] + 1) + ii]};
          auto t_diff_eta = FTensor::Tensor1<double, 3>{
              diff_eta[jj], diff_eta[1 * (p[1] + 1) + jj],
              diff_eta[2 * (p[1] + 1) + jj]};

          t_n(i) = t_cross(i) * ez * mu;
          t_diff_n(i, j) =
              t_cross(i) * ((t_diff_zeta(j) * e + z * t_diff_eta(j)) * mu +
                            ez * t_diff_mu(j));

          ++t_n;
          ++t_diff_n;
        }
      }
    }
  }

  MoFEMFunctionReturnHot(0);
}

MoFEMErrorCode MoFEM::DemkowiczHexAndQuad::Hdiv_InteriorShapeFunctions_ONHEX(
    int *p, double *N, double *diffN, double *volN[3], double *diff_volN[3],
    int nb_integration_pts) {
  MoFEMFunctionBeginHot;

  int pqr[3] = {p[0], p[1], p[2]};
  int qrp[3] = {p[1], p[2], p[0]};
  int rpq[3] = {p[2], p[0], p[1]};

  int perm_fam0[3 * (pqr[0] - 1) * qrp[0] * rpq[0]];
  int perm_fam1[3 * (pqr[1] - 1) * qrp[1] * rpq[1]];
  int perm_fam2[3 * (pqr[2] - 1) * qrp[2] * rpq[2]];

  std::array<int *, 3> permute = {perm_fam0, perm_fam1, perm_fam2};
  for (int fam = 0; fam != 3; ++fam) {
    const int ppp = pqr[fam];
    const int qqq = qrp[fam];
    const int rrr = rpq[fam];
    CHKERR ::DemkowiczHexAndQuad::monom_ordering(permute[fam], ppp - 2, qqq - 1,
                                                 rrr - 1);
  }

  FTensor::Index<'i', 3> i;
  FTensor::Index<'j', 3> j;
  FTensor::Index<'k', 3> k;

  FTensor::Tensor1<double, 3> t_cross;

  //  = {
  // {1., 0., 0.}, {0., 1., 0.}, {0., 0., 1.}};

  for (int qq = 0; qq != nb_integration_pts; ++qq) {

    const int shift = 8 * qq;
    double ksi[3] = {0, 0, 0};
    double diff_ksi[3][3] = {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}};
    ::DemkowiczHexAndQuad::get_ksi_hex(shift, N, diffN, ksi, diff_ksi);

    double ksi_eta_gma[3] = {ksi[0], ksi[1], ksi[2]};
    double eta_gma_ksi[3] = {ksi[1], ksi[2], ksi[0]};
    double gma_ksi_eta[3] = {ksi[2], ksi[0], ksi[1]};

    double *diff_ksi_eta_gma[3] = {diff_ksi[0], diff_ksi[1], diff_ksi[2]};
    double *diff_eta_gma_ksi[3] = {diff_ksi[1], diff_ksi[2], diff_ksi[0]};
    double *diff_gma_ksi_eta[3] = {diff_ksi[2], diff_ksi[0], diff_ksi[1]};

    for (int fam = 0; fam != 3; ++fam) {

      const int ppp = pqr[fam];
      const int qqq = qrp[fam];
      const int rrr = rpq[fam];

      const int nb_dofs = (ppp - 1) * qqq * rrr;
      if (nb_dofs > 0) {

        FTensor::Tensor1<double, 3> t_e1{diff_eta_gma_ksi[fam][0],
                                         diff_eta_gma_ksi[fam][1],
                                         diff_eta_gma_ksi[fam][2]};
        FTensor::Tensor1<double, 3> t_e2{diff_gma_ksi_eta[fam][0],
                                         diff_gma_ksi_eta[fam][1],
                                         diff_gma_ksi_eta[fam][2]};

        t_cross(i) = FTensor::levi_civita(i, j, k) * t_e1(j) * t_e2(k);

        double eta_i[ppp + 2];
        double diff_eta_i[3 * (ppp + 2)];

        CHKERR Lobatto_polynomials(ppp + 1, ksi_eta_gma[fam],
                                   diff_ksi_eta_gma[fam], eta_i, diff_eta_i, 3);

        double phi_j[qqq];
        double diff_phi_j[3 * qqq];

        CHKERR Legendre_polynomials(qqq - 1, eta_gma_ksi[fam],
                                    diff_eta_gma_ksi[fam], phi_j, diff_phi_j,
                                    3);

        double phi_k[rrr];
        double diff_phi_k[3 * rrr];

        CHKERR Legendre_polynomials(rrr - 1, gma_ksi_eta[fam],
                                    diff_gma_ksi_eta[fam], phi_k, diff_phi_k,
                                    3);

        int qd_shift = nb_dofs * qq;
        double *t_n_ptr = &volN[fam][3 * qd_shift];
        double *t_diff_n_ptr = &diff_volN[fam][9 * qd_shift];
        auto t_n = getFTensor1FromPtr<3>(t_n_ptr);
        auto t_diff_n = getFTensor2FromPtr<3, 3>(t_diff_n_ptr);

        for (int n = 0; n != nb_dofs; n++) {
          int ii = permute[fam][3 * n + 0];
          int jj = permute[fam][3 * n + 1];
          int kk = permute[fam][3 * n + 2];

          const double e_i = eta_i[ii + 2];
          const double p_j = phi_j[jj];
          const double p_k = phi_k[kk];

          const double p_jk = p_j * p_k;
          const double ep_ij = e_i * p_j;
          const double ep_ik = e_i * p_k;

          const double a = e_i * p_jk;

          FTensor::Tensor1<double, 3> t_d_a;
          for (int d = 0; d != 3; ++d)
            t_d_a(d) = diff_eta_i[d * (ppp + 2) + ii + 2] * p_jk +
                       diff_phi_j[d * qqq + jj] * ep_ik +
                       diff_phi_k[d * rrr + kk] * ep_ij;

          t_n(i) = a * t_cross(i);
          t_diff_n(i, j) = t_cross(i) * t_d_a(j);

          ++t_n;
          ++t_diff_n;
        }
      }
    }
  }

  MoFEMFunctionReturnHot(0);
}
