/** \file EdgeQuadHexPolynomials.cpp

  \brief Implementation of hierarchical Edge, Quad, and Hex shape functions of
  type H1, Hcurl, Hdiv
*/

using namespace MoFEM;

struct RefHex {
  RefHex()
      : vertices{{0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}, {1.0, 1.0, 0.0},
                 {0.0, 1.0, 0.0}, {0.0, 0.0, 1.0}, {1.0, 0.0, 1.0},
                 {1.0, 1.0, 1.0}, {0.0, 1.0, 1.0}},
        faces{{0, 1, 2, 3}, {0, 1, 5, 4}, {1, 2, 6, 5},
              {3, 2, 6, 7}, {0, 3, 7, 4}, {4, 5, 6, 7}} {}

  void get_volume_coords(double (&Nq)[8], double (&volume_coords)[3]) {
    volume_coords[0] = 0.0;
    volume_coords[1] = 0.0;
    volume_coords[2] = 0.0;
    for (int vv = 0; vv < 8; vv++) {
      volume_coords[0] += Nq[vv] * vertices[vv][0];
      volume_coords[1] += Nq[vv] * vertices[vv][1];
      volume_coords[2] += Nq[vv] * vertices[vv][2];
    }
  }

  void get_volume_diff_coords(double (&volume_diff_coords)[3][3]) {
    double diff_coords[3][3] = {
        {1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {0.0, 0.0, 1.0}};
    for (int c1 = 0; c1 < 3; c1++)
      for (int c2 = 0; c2 < 3; c2++)
        volume_diff_coords[c1][c2] = diff_coords[c1][c2];
  }

  void get_edge_affines(double (&volume_coords)[3],
                        double (&edge_affines)[12][2]) {

    double ksi = volume_coords[0];
    double eta = volume_coords[1];
    double gma = volume_coords[2];

    double affines[12][2] = {
        {1.0 - eta, 1.0 - gma}, {0.0 + ksi, 1.0 - gma}, {0.0 + eta, 1.0 - gma},
        {1.0 - ksi, 1.0 - gma}, {1.0 - ksi, 1.0 - eta}, {0.0 + ksi, 1.0 - eta},
        {0.0 + ksi, 0.0 + eta}, {1.0 - ksi, 0.0 + eta}, {1.0 - eta, 0.0 + gma},
        {0.0 + ksi, 0.0 + gma}, {0.0 + eta, 0.0 + gma}, {1.0 - ksi, 0.0 + gma}};

    for (int ee = 0; ee < 12; ee++)
      for (int comp = 0; comp < 3; comp++)
        edge_affines[ee][comp] = affines[ee][comp];
  }

  void get_edge_diff_affines(double (&edge_diff_affines)[12][2][3]) {

    double diff_affines[3][2][3] = {{{-1.0, 0.0, 0.0}, {1.0, 0.0, 0.0}},
                                    {{0.0, -1.0, 0.0}, {0.0, 1.0, 0.0}},
                                    {{0.0, 0.0, -1.0}, {0.0, 0.0, 1.0}}};
    int II[12][2][2] = {{{1, 0}, {2, 0}}, {{0, 1}, {2, 0}}, {{1, 1}, {2, 0}},
                        {{0, 0}, {2, 0}}, {{0, 0}, {1, 0}}, {{0, 1}, {1, 0}},
                        {{0, 1}, {1, 1}}, {{0, 0}, {1, 1}}, {{1, 0}, {2, 1}},
                        {{0, 1}, {2, 1}}, {{1, 1}, {2, 1}}, {{0, 0}, {2, 1}}};
    for (int ee = 0; ee < 12; ee++) {
      int comp11 = II[ee][0][0];
      int comp12 = II[ee][0][1];
      int comp21 = II[ee][1][0];
      int comp22 = II[ee][1][1];
      for (int c = 0; c < 3; c++) {
        edge_diff_affines[ee][0][c] = diff_affines[comp11][comp12][c];
        edge_diff_affines[ee][1][c] = diff_affines[comp21][comp22][c];
      }
    }
  }
  void get_edge_coords(int *sense, double (&volume_coords)[3],
                       double (&edge_coords)[12]) {

    int free_edge_coords[12] = {0, 1, 0, 1, 2, 2, 2, 2, 0, 1, 0, 1};

    for (int ee = 0; ee < 12; ee++) {
      int cc = free_edge_coords[ee];
      edge_coords[ee] = (double)sense[ee] * volume_coords[cc] +
                        0.5 * (1.0 - (double)sense[ee]);
    }
  }
  void get_edge_diff_coords(int *sense, double (&edge_diff_coords)[12][3]) {

    double diff_coords[3][3] = {
        {1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {0.0, 0.0, 1.0}};
    int free_edge_coords[12] = {0, 1, 0, 1, 2, 2, 2, 2, 0, 1, 0, 1};
    for (int ee = 0; ee < 12; ee++) {
      for (int cc = 0; cc < 3; cc++) {
        int n = free_edge_coords[ee];
        edge_diff_coords[ee][cc] = (double)sense[ee] * diff_coords[n][cc];
      }
    }
  }
  void get_face_affines(double (&volume_coords)[3], double (&face_affines)[6]) {

    double ksi = volume_coords[0];
    double eta = volume_coords[1];
    double gma = volume_coords[2];

    double faceAffine[6] = {1.0 - gma, 1.0 - eta, 0.0 + ksi,
                            0.0 + eta, 1.0 - ksi, 0.0 + gma};
    for (int ff = 0; ff != 6; ff++)
      face_affines[ff] = faceAffine[ff];
  }

  void get_face_diff_affines(double (&face_diff_affines)[6][3]) {

    int free_face[6][2] = {{2, 0}, {1, 0}, {0, 1}, {1, 1}, {0, 0}, {2, 1}};

    double diff_affines[3][2][3] = {{{-1.0, 0.0, 0.0}, {1.0, 0.0, 0.0}},
                                    {{0.0, -1.0, 0.0}, {0.0, 1.0, 0.0}},
                                    {{0.0, 0.0, -1.0}, {0.0, 0.0, 1.0}}};
    for (int ff = 0; ff < 6; ff++) {
      int var_i = free_face[ff][0];
      int var_j = free_face[ff][1];
      for (int cc = 0; cc < 3; cc++) {
        face_diff_affines[ff][cc] = diff_affines[var_i][var_j][cc];
      }
    }
  }

  void get_face_coords(int *face_nodes[6], double (&Nq)[8],
                       double (&face_coords)[6][2]) {

    int par_face_nodes[6][8] = {
        {4, 5, 6, 7, 0, 1, 2, 3}, {3, 2, 1, 0, 7, 6, 5, 4},
        {1, 0, 3, 2, 5, 4, 7, 6}, {3, 2, 1, 0, 7, 6, 5, 4},
        {1, 0, 3, 2, 5, 4, 7, 6}, {4, 5, 6, 7, 0, 1, 2, 3}};

    int free_coords[6][2] = {{0, 1}, {0, 2}, {1, 2}, {0, 2}, {1, 2}, {0, 1}};

    for (int ff = 0; ff != 6; ff++) {
      int v0 = free_coords[ff][0];
      int v1 = free_coords[ff][1];
      face_coords[ff][0] = 0.0;
      face_coords[ff][1] = 0.0;
      for (int fv = 0; fv != 4; fv++) {
        int n0 = face_nodes[ff][fv];
        int n1 = par_face_nodes[ff][n0];
        int index = faces[ff][fv];
        double N = Nq[n0] + Nq[n1];
        face_coords[ff][0] += vertices[index][v0] * N;
        face_coords[ff][1] += vertices[index][v1] * N;
      }
    }
  }

  void get_face_diff_coords(int *face_nodes[6],
                            double (&face_diff_coords)[6][2][3]) {

    int I1[8] = {0, 1, 1, 0, 0, 1, 1, 0};
    int I2[8] = {0, 0, 1, 1, 0, 0, 1, 1};
    int I3[8] = {0, 0, 0, 0, 1, 1, 1, 1};

    double n1[8] = {-1.0, 1.0, 1.0, -1.0, -1.0, 1.0, 1.0, -1.0};
    double n2[8] = {-1.0, -1.0, 1.0, 1.0, -1.0, -1.0, 1.0, 1.0};
    double n3[8] = {-1.0, -1.0, -1.0, -1.0, 1.0, 1.0, 1.0, 1.0};

    double diff_N[8][3];

    for (int n = 0; n < 8; n++) {
      diff_N[n][0] = n1[n] * 0.25;
      diff_N[n][1] = n2[n] * 0.25;
      diff_N[n][2] = n3[n] * 0.25;
    }

    int par_face_nodes[6][8] = {
        {4, 5, 6, 7, 0, 1, 2, 3}, {3, 2, 1, 0, 7, 6, 5, 4},
        {1, 0, 3, 2, 5, 4, 7, 6}, {3, 2, 1, 0, 7, 6, 5, 4},
        {1, 0, 3, 2, 5, 4, 7, 6}, {4, 5, 6, 7, 0, 1, 2, 3}};

    int free_coords[6][2] = {{0, 1}, {0, 2}, {1, 2}, {0, 2}, {1, 2}, {0, 1}};

    for (int ff = 0; ff != 6; ff++) {
      int v0 = free_coords[ff][0];
      int v1 = free_coords[ff][1];
      face_diff_coords[ff][0][0] = 0.0;
      face_diff_coords[ff][0][1] = 0.0;
      face_diff_coords[ff][0][2] = 0.0;

      face_diff_coords[ff][1][0] = 0.0;
      face_diff_coords[ff][1][1] = 0.0;
      face_diff_coords[ff][1][2] = 0.0;
      for (int fv = 0; fv != 4; fv++) {
        int n0 = face_nodes[ff][fv];
        int n1 = par_face_nodes[ff][n0];
        int index = faces[ff][fv];
        double Nx = diff_N[n0][0] + diff_N[n1][0];
        double Ny = diff_N[n0][1] + diff_N[n1][1];
        double Nz = diff_N[n0][2] + diff_N[n1][2];
        face_diff_coords[ff][0][0] += vertices[index][v0] * Nx;
        face_diff_coords[ff][0][1] += vertices[index][v0] * Ny;
        face_diff_coords[ff][0][2] += vertices[index][v0] * Nz;

        face_diff_coords[ff][1][0] += vertices[index][v1] * Nx;
        face_diff_coords[ff][1][1] += vertices[index][v1] * Ny;
        face_diff_coords[ff][1][2] += vertices[index][v1] * Nz;
      }
    }
  }

  void Cross_product(double (&v1)[3], double (&v2)[3], double (&result)[3]) {

    result[0] = v1[1] * v2[2] - v1[2] * v2[1];
    result[1] = v1[2] * v2[0] - v1[0] * v2[2];
    result[2] = v1[0] * v2[1] - v1[1] * v2[0];
  }

private:
  double vertices[8][3];
  int faces[6][4];
};

MoFEMErrorCode MoFEM::DemkowiczHexAndQuad::MonomOrdering(int perm[][3], int p,
                                                         int q, int r) {
  MoFEMFunctionBeginHot;
  int n = 0;
  for (int m = 0; m != max(max(p, q), r) + 1; m++) {
    for (int i = 0; i != min(m, p) + 1; i++) {
      for (int j = 0; j != min(m, q) + 1; j++) {
        for (int k = 0; k != min(m, r) + 1; k++) {
          if (i == m || j == m || k == m) {
            perm[n][0] = i;
            perm[n][1] = j;
            perm[n][2] = k;
            n++;
          } else
            continue;
        }
      }
    }
  }
  MoFEMFunctionReturnHot(0);
}
// Auxilary functions

MoFEMErrorCode MoFEM::DemkowiczHexAndQuad::Legendre_polynomials01(int p,
                                                                  double s01,
                                                                  double *L) {

  MoFEMFunctionBeginHot;
  double s = 2.0 * s01 - 1.0;
  CHKERR Legendre_polynomials(p, s, NULL, L, NULL, 1);
  MoFEMFunctionReturnHot(0);
}

MoFEMErrorCode
MoFEM::DemkowiczHexAndQuad::Integrated_Legendre01(int p, double s01, double *L,
                                                  double *diffL) {
  MoFEMFunctionBeginHot;
  double l[p + 3];
  CHKERR Legendre_polynomials01(p, s01, l);
  if (p >= 2) {
    for (int i = 0; i != p - 1; i++) {
      double factor = 1.0 / (2 * (2.0 * (i + 1) + 1.0));
      L[i] = factor * (l[i + 2] - l[i]);
      diffL[i] = l[i + 1];
    }
  }
  MoFEMFunctionReturnHot(0);
}
/*

    0--------------1

*/
MoFEMErrorCode MoFEM::DemkowiczHexAndQuad::H1_BubbleShapeFunctions_ONSEGMENT(
    int p, double *N, double *diffN, double *bubbleN, double *diff_bubbleN,
    int nb_integration_pts) {
  MoFEMFunctionBeginHot;
  constexpr int n0 = 0;
  constexpr int n1 = 1;
  double diff_mu = diffN[n1] - diffN[n0];
  for (int q = 0; q != nb_integration_pts; q++) {
    int shift = 2 * q;
    double mu = N[shift + n1] - N[shift + n0];
    double L[p + 2];
    double diffL[p + 2];
    CHKERR Lobatto_polynomials(p + 1, mu, &diff_mu, L, diffL, 1);
    int qd_shift = (p - 1) * q;
    for (int n = 0; n != p - 1; n++) {
      bubbleN[qd_shift + n] = L[n + 2];
      diff_bubbleN[qd_shift + n] = diffL[n + 2];
    }
  }
  MoFEMFunctionReturnHot(0);
}

MoFEMErrorCode MoFEM::DemkowiczHexAndQuad::L2_ShapeFunctions_ONSEGMENT(
    int p, double *N, double *diffN, double *funN, double *funDiffN,
    int nb_integration_pts) {
  MoFEMFunctionBeginHot;
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
  MoFEMFunctionReturnHot(0);
}

MoFEMErrorCode MoFEM::DemkowiczHexAndQuad::H1_FaceShapeFunctions_ONQUAD(
    int *faces_nodes, int *p, double *N, double *diffN, double *faceN,
    double *diff_faceN, int nb_integration_pts) {

  const int n0 = faces_nodes[0];
  const int n1 = faces_nodes[1];
  const int n2 = faces_nodes[2];
  const int n3 = faces_nodes[3];

  MoFEMFunctionBeginHot;
  int permute[(p[0] - 1) * (p[1] - 1)][3];
  CHKERR MonomOrdering(permute, p[0] - 2, p[1] - 2);
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

    int qd_shift = (p[0] - 1) * (p[1] - 1) * q;
    for (int n = 0; n != (p[0] - 1) * (p[1] - 1); ++n) {
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
  MoFEMFunctionReturnHot(0);
}

MoFEMErrorCode MoFEM::DemkowiczHexAndQuad::L2_FaceShapeFunctions_ONQUAD(
    int *p, double *N, double *diffN, double *faceN, double *diff_faceN,
    int nb_integration_pts) {
  MoFEMFunctionBeginHot;
  int permute[(p[0] + 1) * (p[1] + 1)][3];
  CHKERR MonomOrdering(permute, p[0], p[1]);

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

    int qd_shift = (p[0] + 1) * (p[1] + 1) * q;
    for (int n = 0; n != (p[0] + 1) * (p[1] + 1); ++n) {
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

    int pp[4] = {p[0], p[1], p[2], p[3]};
    for (int e = 0; e != 4; e++) {

      double L[pp[e]];
      double diffL[2 * pp[e]];

      CHKERR Legendre_polynomials(pp[e] - 1, mu[e], diff_mu[e], L, diffL, 2);

      int qd_shift = pp[e] * q;
      double *t_n_ptr = &edgeN[e][3 * qd_shift];
      double *t_diff_n_ptr = &diff_edgeN[e][3 * 2 * qd_shift];
      auto t_n = getFTensor1FromPtr<3>(t_n_ptr);
      auto t_diff_n = getFTensor2FromPtr<3, 2>(t_diff_n_ptr);

      for (int n = 0; n != pp[e]; ++n) {
        const double a = mu_const[e] * L[n];
        const double d_a[] = {
            diff_mu_const[e][0] * L[n] + mu_const[e] * diffL[0 * pp[e] + n],
            diff_mu_const[e][1] * L[n] + mu_const[e] * diffL[1 * pp[e] + n]};

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
  MoFEMFunctionReturnHot(0);
}

MoFEMErrorCode MoFEM::DemkowiczHexAndQuad::Hcurl_FaceShapeFunctions_ONQUAD(
    int *face_nodes, int *p, double *N, double *diffN, double *faceN[],
    double *diff_faceN[], int nb_integration_pts) {
  MoFEMFunctionBeginHot;

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

    int pq[2] = {p[0], p[1]};
    int qp[2] = {p[1], p[0]};
    double mu_ksi_eta[2] = {ksi01, ksi12};
    double mu_eta_ksi[2] = {ksi12, ksi01};

    double *diff_mu_ksi_eta[2] = {diff_ksi01, diff_ksi12};
    double *diff_mu_eta_ksi[2] = {diff_ksi12, diff_ksi01};

    for (int family = 0; family != 2; family++) {
      const int pp = pq[family];
      const int qq = qp[family];

      double Phi[pp + 2];
      double diffPhi[2 * (pp + 2)];
      CHKERR Lobatto_polynomials(qq + 1, mu_ksi_eta[family],
                                 diff_mu_ksi_eta[family], Phi, diffPhi, 2);

      double E[qq];
      double diffE[2 * qq];
      CHKERR Legendre_polynomials(pp - 1, mu_eta_ksi[family],
                                  diff_mu_eta_ksi[family], E, diffE, 2);

      int permute[qq * (pp - 1)][3];
      CHKERR MonomOrdering(permute, qq - 1, pp - 2);

      const int qd_shift = (pp - 1) * qq * q;
      double *t_n_ptr = &faceN[family][3 * qd_shift];
      double *t_diff_n_ptr = &diff_faceN[family][3 * 2 * qd_shift];
      auto t_n = getFTensor1FromPtr<3>(t_n_ptr);
      auto t_diff_n = getFTensor2FromPtr<3, 2>(t_diff_n_ptr);

      for (int n = 0; n != (pp - 1) * qq; n++) {
        int i = permute[n][0];
        int j = permute[n][1];

        const double phi = Phi[j + 2];
        const double e = E[i];
        const double a = phi * e;
        const double d_a[] = {

            diffPhi[0 * (pp + 2) + j + 2] * e + phi * diffE[0 * qq + i],

            diffPhi[1 * (pp + 2) + j + 2] * e + phi * diffE[1 * qq + i]};

        for (int d = 0; d != 2; ++d) {
          t_n(d) = a * diff_mu_eta_ksi[family][d];
          for (int m = 0; m != 2; ++m) {
            t_diff_n(d, m) = d_a[m] * diff_mu_eta_ksi[family][d];
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
  MoFEMFunctionReturnHot(0);
}

MoFEMErrorCode MoFEM::DemkowiczHexAndQuad::Hdiv_FaceShapeFunctions_ONQUAD(
    int *p, double *N, double *faceN[], double *div_faceN[],
    int nb_integration_pts) {
  MoFEMFunctionBeginHot;

  // Here should be trace of H-div space from Hex

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
    int *sense, int *p, double *N, double *edgeN[12], double *diff_edgeN[12],
    int nb_integration_pts) {
  MoFEMFunctionBeginHot;
  RefHex ref_hex;

  for (int q = 0; q != nb_integration_pts; q++) {
    // General *******************************
    int shift = 8 * q;
    double quad_coords[3];
    double Nq[8];
    for (int vv = 0; vv < 8; vv++)
      Nq[vv] = N[shift + vv];

    ref_hex.get_volume_coords(Nq, quad_coords);
    // ***************************************

    double mu[12][2];
    ref_hex.get_edge_affines(quad_coords, mu);

    double diff_mu[12][2][3];
    ref_hex.get_edge_diff_affines(diff_mu);

    double ksi[12];
    ref_hex.get_edge_coords(sense, quad_coords, ksi);
    double diff_ksi[12][3];
    ref_hex.get_edge_diff_coords(sense, diff_ksi);

    int pp[12] = {p[0], p[1], p[0], p[1], p[2], p[2],
                  p[2], p[2], p[0], p[1], p[0], p[1]};

    for (int e = 0; e != 12; e++) {
      double L[pp[e] - 1];
      double diffL[pp[e] - 1];

      CHKERR Integrated_Legendre01(pp[e], ksi[e], L, diffL);
      int qd_shift = (pp[e] - 1) * q;

      for (int n = 0; n != pp[e] - 1; n++) {
        edgeN[e][qd_shift + n] = mu[e][0] * mu[e][1] * L[n];

        diff_edgeN[e][3 * (qd_shift + n) + 0] =
            mu[e][0] * mu[e][1] * diffL[n] * diff_ksi[e][0] +
            (mu[e][0] * diff_mu[e][0][0] + mu[e][1] * diff_mu[e][1][0]) * L[n];

        diff_edgeN[e][3 * (qd_shift + n) + 1] =
            mu[e][0] * mu[e][1] * diffL[n] * diff_ksi[e][1] +
            (mu[e][0] * diff_mu[e][0][1] + mu[e][1] * diff_mu[e][1][1]) * L[n];

        diff_edgeN[e][3 * (qd_shift + n) + 0] =
            mu[e][0] * mu[e][1] * diffL[n] * diff_ksi[e][2] +
            (mu[e][0] * diff_mu[e][0][2] + mu[e][1] * diff_mu[e][1][2]) * L[n];
      }
    }
  }
  MoFEMFunctionReturnHot(0);
}
// TODO
MoFEMErrorCode MoFEM::DemkowiczHexAndQuad::H1_FaceShapeFunctions_ONHEX(
    int *face_nodes[6], int *p, double *N, double *faceN[6],
    double *diff_faceN[6], int nb_integration_pts) {
  MoFEMFunctionBeginHot;

  RefHex ref_hex;

  for (int q = 0; q != nb_integration_pts; q++) {
    // general ******************************
    int shift = 8 * q;
    double quad_coords[3];
    double Nq[8];
    for (int vv = 0; vv < 8; vv++)
      Nq[vv] = N[shift + vv];

    ref_hex.get_volume_coords(Nq, quad_coords);
    // ****************************************

    double mu[6];
    ref_hex.get_face_affines(quad_coords, mu);
    double diff_mu[6][3];
    ref_hex.get_face_diff_affines(diff_mu);

    double ksi[6][2];
    ref_hex.get_face_coords(face_nodes, Nq, ksi);
    double diff_ksi[6][2][3];
    ref_hex.get_face_diff_coords(face_nodes, diff_ksi);

    int pp[6][2] = {{p[0], p[1]}, {p[0], p[2]}, {p[1], p[2]},
                    {p[0], p[2]}, {p[1], p[2]}, {p[0], p[1]}};

    for (int face = 0; face != 6; face++) {
      int pq[2] = {pp[face][0], pp[face][1]};

      double L0[pq[0] - 1];
      double diffL0[pq[0] - 1];
      double L1[pq[1] - 1];
      double diffL1[pq[1] - 1];

      double ss0 = ksi[face][0];
      double ss1 = ksi[face][1];

      CHKERR Integrated_Legendre01(pq[0], ss0, L0, diffL0);
      CHKERR Integrated_Legendre01(pq[1], ss1, L1, diffL1);

      int permute[(pq[0] - 1) * (pq[1] - 1)][3];
      CHKERR MonomOrdering(permute, pq[0] - 2, pq[1] - 2);

      int qd_shift = (pq[0] - 1) * (pq[1] - 1) * q;
      int n = 0;
      for (; n != (pq[0] - 1) * (pq[1] - 1); n++) {
        int s1 = permute[n][0];
        int s2 = permute[n][1];
        faceN[face][qd_shift + n] = mu[face] * L0[s1] * L1[s2];

        diff_faceN[face][3 * (qd_shift + n) + 0] =
            diff_mu[face][0] * L0[s1] * L1[s2] +
            mu[face] * diffL0[s1] * diff_ksi[face][0][0] * L1[s2] +
            mu[face] * L0[s1] * diffL1[s1] * diff_ksi[face][1][0];

        diff_faceN[face][3 * (qd_shift + n) + 1] =
            diff_mu[face][1] * L0[s1] * L1[s2] +
            mu[face] * diffL0[s1] * diff_ksi[face][0][1] * L1[s2] +
            mu[face] * L0[s1] * diffL1[s1] * diff_ksi[face][1][1];

        diff_faceN[face][3 * (qd_shift + n) + 2] =
            diff_mu[face][2] * L0[s1] * L1[s2] +
            mu[face] * diffL0[s1] * diff_ksi[face][0][2] * L1[s2] +
            mu[face] * L0[s1] * diffL1[s1] * diff_ksi[face][1][2];
      }

      // for (int s1 = 0; s1 != pq[0] - 1; s1++) {
      //   for (int s2 = 0; s2 != pq[1] - 1; s2++) {
      //     faceN[face][qd_shift + n] = mu[face] * L0[s1] * L1[s2];

      //     diff_faceN[face][3 * (qd_shift + n) + 0] =
      //         diff_mu[face][0] * L0[s1] * L1[s2] +
      //         mu[face] * diffL0[s1] * diff_ksi[face][0][0] * L1[s2] +
      //         mu[face] * L0[s1] * diffL1[s1] * diff_ksi[face][1][0];

      //     diff_faceN[face][3 * (qd_shift + n) + 1] =
      //         diff_mu[face][1] * L0[s1] * L1[s2] +
      //         mu[face] * diffL0[s1] * diff_ksi[face][0][1] * L1[s2] +
      //         mu[face] * L0[s1] * diffL1[s1] * diff_ksi[face][1][1];

      //     diff_faceN[face][3 * (qd_shift + n) + 2] =
      //         diff_mu[face][2] * L0[s1] * L1[s2] +
      //         mu[face] * diffL0[s1] * diff_ksi[face][0][2] * L1[s2] +
      //         mu[face] * L0[s1] * diffL1[s1] * diff_ksi[face][1][2];

      //     ++n;
      //   }
      // }
    }
  }
  MoFEMFunctionReturnHot(0);
}

MoFEMErrorCode MoFEM::DemkowiczHexAndQuad::H1_InteriorShapeFunctions_ONHEX(
    int *p, double *N, double *faceN, double *diff_faceN,
    int nb_integration_pts) {
  MoFEMFunctionBeginHot;
  RefHex ref_hex;

  int permute[(p[0] - 1) * (p[1] - 1) * (p[2] - 1)][3];
  CHKERR MonomOrdering(permute, p[0] - 2, p[1] - 2, p[2] - 2);

  for (int q = 0; q != nb_integration_pts; q++) {

    // general ******************************
    int shift = 8 * q;
    double ksi[3];
    double Nq[8];
    for (int vv = 0; vv < 8; vv++)
      Nq[vv] = N[shift + vv];

    ref_hex.get_volume_coords(Nq, ksi);
    // ****************************************

    double diff_ksi[3][3];

    ref_hex.get_volume_diff_coords(diff_ksi);

    double L0[p[0] - 1];
    double diffL0[p[0] - 1];
    double L1[p[1] - 1];
    double diffL1[p[1] - 1];
    double L2[p[2] - 1];
    double diffL2[p[2] - 1];

    CHKERR Integrated_Legendre01(p[0], ksi[0], L0, diffL0);
    CHKERR Integrated_Legendre01(p[1], ksi[1], L1, diffL1);
    CHKERR Integrated_Legendre01(p[2], ksi[2], L2, diffL2);

    // cout << "In Face H1" << endl;

    int qd_shift = (p[0] - 1) * (p[1] - 1) * (p[2] - 1) * q;
    int n = 0;
    for (; n != (p[0] - 1) * (p[1] - 1) * (p[2] - 1); n++) {
      int s1 = permute[n][0];
      int s2 = permute[n][1];
      int s3 = permute[n][2];

      faceN[qd_shift + n] = L0[s1] * L1[s2] * L2[s3];

      diff_faceN[3 * (qd_shift + n) + 0] =
          diffL0[s1] * diff_ksi[0][0] * L1[s2] * L2[s3] +
          L0[s1] * diffL1[s2] * diff_ksi[1][0] * L2[s3] +
          L0[s1] * L1[s2] * diffL2[s3] * diff_ksi[2][0];

      diff_faceN[3 * (qd_shift + n) + 1] =
          diffL0[s1] * diff_ksi[0][1] * L1[s2] * L2[s3] +
          L0[s1] * diffL1[s2] * diff_ksi[1][1] * L2[s3] +
          L0[s1] * L1[s2] * diffL2[s3] * diff_ksi[2][1];

      diff_faceN[3 * (qd_shift + n) + 2] =
          diffL0[s1] * diff_ksi[0][2] * L1[s2] * L2[s3] +
          L0[s1] * diffL1[s2] * diff_ksi[1][2] * L2[s3] +
          L0[s1] * L1[s2] * diffL2[s3] * diff_ksi[2][2];
    }
    // for (int s1 = 0; s1 != p[0] - 1; s1++) {
    //   for (int s2 = 0; s2 != p[1] - 1; s2++) {
    //     for (int s3 = 0; s3 < p[2] - 1; s3++) {
    //       faceN[qd_shift + n] = L0[s1] * L1[s2] * L2[s3];

    //       diff_faceN[3 * (qd_shift + n) + 0] =
    //           diffL0[s1] * diff_ksi[0][0] * L1[s2] * L2[s3] +
    //           L0[s1] * diffL1[s2] * diff_ksi[1][0] * L2[s3] +
    //           L0[s1] * L1[s2] * diffL2[s3] * diff_ksi[2][0];

    //       diff_faceN[3 * (qd_shift + n) + 1] =
    //           diffL0[s1] * diff_ksi[0][1] * L1[s2] * L2[s3] +
    //           L0[s1] * diffL1[s2] * diff_ksi[1][1] * L2[s3] +
    //           L0[s1] * L1[s2] * diffL2[s3] * diff_ksi[2][1];

    //       diff_faceN[3 * (qd_shift + n) + 2] =
    //           diffL0[s1] * diff_ksi[0][2] * L1[s2] * L2[s3] +
    //           L0[s1] * diffL1[s2] * diff_ksi[1][2] * L2[s3] +
    //           L0[s1] * L1[s2] * diffL2[s3] * diff_ksi[2][2];

    //       ++n;
    //     }
    //   }
    // }
  }
  MoFEMFunctionReturnHot(0);
}

MoFEMErrorCode MoFEM::DemkowiczHexAndQuad::L2_InteriorShapeFunctions_ONHEX(
    int *p, double *N, double *volN, int nb_integration_pts) {
  MoFEMFunctionBeginHot;
  RefHex ref_hex;
  int permute[p[0] * p[1] * p[2]][3];
  CHKERR MonomOrdering(permute, p[0] - 1, p[1] - 1, p[2] - 1);
  for (int qq = 0; qq != nb_integration_pts; qq++) {

    // general ******************************
    int shift = 8 * qq;
    double ksi[3];
    double Nq[8];
    for (int vv = 0; vv < 8; vv++)
      Nq[vv] = N[shift + vv];

    ref_hex.get_volume_coords(Nq, ksi);
    // ****************************************

    double P0[p[0]];
    double P1[p[1]];
    double P2[p[2]];

    int qd_shift = qq * p[0] * p[1] * p[2];

    CHKERR Legendre_polynomials01(p[0] - 1, ksi[0], P0);
    CHKERR Legendre_polynomials01(p[1] - 1, ksi[1], P1);
    CHKERR Legendre_polynomials01(p[2] - 1, ksi[2], P2);
    int n = 0;
    for (; n != p[0] * p[1] * p[2]; n++) {
      int ii = permute[n][0];
      int jj = permute[n][1];
      int kk = permute[n][2];

      volN[qd_shift + n] = P0[ii] * P1[jj] * P2[kk];
    }

    // for (int ii = 0; ii != p[0]; ii++) {
    //   for (int jj = 0; jj != p[1]; jj++) {
    //     for (int kk = 0; kk != p[2]; kk++) {
    //       volN[qd_shift + n] = P0[ii] * P1[jj] * P2[kk];
    //     }
    //   }
    // }
  }

  MoFEMFunctionReturnHot(0);
}

MoFEMErrorCode MoFEM::DemkowiczHexAndQuad::Hcurl_EdgeShapeFunctions_ONHEX(
    int *sense, int *p, double *N, double *edgeN[12], double *curl_edgeN[12],
    int nb_integration_pts) {
  MoFEMFunctionBeginHot;
  RefHex ref_hex;

  for (int qq = 0; qq != nb_integration_pts; qq++) {
    // general ******************************
    int shift = 8 * qq;
    double quad_coords[3];
    double Nq[8];
    for (int vv = 0; vv < 8; vv++)
      Nq[vv] = N[shift + vv];

    ref_hex.get_volume_coords(Nq, quad_coords);
    // ****************************************

    double mu[12][2];
    ref_hex.get_edge_affines(quad_coords, mu);

    double diff_mu[12][2][3];
    ref_hex.get_edge_diff_affines(diff_mu);

    double ksi[12];
    ref_hex.get_edge_coords(sense, quad_coords, ksi);
    double diff_ksi[12][3];
    ref_hex.get_edge_diff_coords(sense, diff_ksi);

    int pp[12] = {p[0], p[1], p[0], p[1], p[2], p[2],
                  p[2], p[2], p[0], p[1], p[0], p[1]};
    for (int ee = 0; ee != 12; ee++) {

      double L[pp[ee]];

      CHKERR Legendre_polynomials01(pp[ee] - 1, ksi[ee], L);
      int qd_shift = pp[ee] * qq;

      for (int ii = 0; ii != pp[ee]; ii++) {
        edgeN[ee][3 * (qd_shift + ii) + 0] =
            mu[ee][0] * mu[ee][1] * L[ii] * diff_ksi[ee][0];
        edgeN[ee][3 * (qd_shift + ii) + 1] =
            mu[ee][0] * mu[ee][1] * L[ii] * diff_ksi[ee][1];
        edgeN[ee][3 * (qd_shift + ii) + 2] =
            mu[ee][0] * mu[ee][1] * L[ii] * diff_ksi[ee][2];

        double E1[3] = {
            diff_mu[ee][0][0] * mu[ee][1] + mu[ee][0] * diff_mu[ee][1][0],
            diff_mu[ee][0][1] * mu[ee][1] + mu[ee][0] * diff_mu[ee][1][1],
            diff_mu[ee][0][2] * mu[ee][1] + mu[ee][0] * diff_mu[ee][1][2]};
        double E2[3] = {L[ii] * diff_ksi[ee][0], L[ii] * diff_ksi[ee][1],
                        L[ii] * diff_ksi[ee][2]};

        double E1_X_E2[3];
        ref_hex.Cross_product(E1, E2, E1_X_E2);

        curl_edgeN[ee][3 * (qd_shift + ii) + 0] = E1_X_E2[0];
        curl_edgeN[ee][3 * (qd_shift + ii) + 1] = E1_X_E2[1];
        curl_edgeN[ee][3 * (qd_shift + ii) + 2] = E1_X_E2[2];
      }
    }
  }
  MoFEMFunctionReturnHot(0);
}

MoFEMErrorCode MoFEM::DemkowiczHexAndQuad::Hcurl_FaceShapeFunctions_ONHEX(
    int *face_nodes[6], int *p, double *N, double *faceN[6][2],
    double *curl_faceN[6][2], int nb_integration_pts) {
  MoFEMFunctionBeginHot;
  RefHex ref_hex;
  for (int qq = 0; qq != nb_integration_pts; qq++) {

    for (int qq = 0; qq != nb_integration_pts; qq++) {
      // general ******************************
      int shift = 8 * qq;
      double quad_coords[3];
      double Nq[8];
      for (int vv = 0; vv < 8; vv++)
        Nq[vv] = N[shift + vv];

      ref_hex.get_volume_coords(Nq, quad_coords);
      // ****************************************

      double mu[6];
      ref_hex.get_face_affines(quad_coords, mu);
      double diff_mu[6][3];
      ref_hex.get_face_diff_affines(diff_mu);

      double ksi[6][2];
      ref_hex.get_face_coords(face_nodes, Nq, ksi);
      double diff_ksi[6][2][3];
      ref_hex.get_face_diff_coords(face_nodes, diff_ksi);

      int px[6][2] = {{p[0], p[1]}, {p[0], p[2]}, {p[1], p[2]},
                      {p[0], p[2]}, {p[1], p[2]}, {p[0], p[1]}};

      for (int ff = 0; ff != 6; ff++) {
        double ksi_eta[2] = {ksi[ff][0], ksi[ff][1]};
        double eta_ksi[2] = {ksi[ff][1], ksi[ff][0]};

        double *diff_ksi_eta[2] = {diff_ksi[ff][0], diff_ksi[ff][1]};
        double *diff_eta_ksi[2] = {diff_ksi[ff][1], diff_ksi[ff][0]};

        int pq[2] = {px[ff][0], px[ff][1]};
        int qp[2] = {px[ff][1], px[ff][0]};

        for (int fam = 0; fam != 2; fam++) {

          double Phi[pq[fam] - 1];
          double diffPhi[pq[fam] - 1];
          CHKERR Integrated_Legendre01(pq[fam], ksi_eta[fam], Phi, diffPhi);

          double E[qp[fam]];
          CHKERR Legendre_polynomials01(qp[fam] - 1, eta_ksi[fam], E);

          int permute[(pq[fam] - 1) * qp[fam]][3];
          CHKERR MonomOrdering(permute, qp[fam] - 1, pq[fam] - 2);

          int qd_shift = (pq[fam] - 1) * qp[fam] * qq;
          int n = 0;
          for (; n != (pq[fam] - 1) * qp[fam]; n++) {
            int i = permute[n][0];
            int j = permute[n][1];

            faceN[ff][fam][3 * (qd_shift + n) + 0] =
                mu[ff] * Phi[j] * E[i] * diff_eta_ksi[fam][0];
            faceN[ff][fam][3 * (qd_shift + n) + 1] =
                mu[ff] * Phi[j] * E[i] * diff_eta_ksi[fam][1];
            faceN[ff][fam][3 * (qd_shift + n) + 2] =
                mu[ff] * Phi[j] * E[i] * diff_eta_ksi[fam][2];

            double Ei[3] = {E[i] * diff_eta_ksi[fam][0],
                            E[i] * diff_eta_ksi[fam][1],
                            E[i] * diff_eta_ksi[fam][2]};

            double Eij[3] = {Phi[j] * Ei[0], Phi[j] * Ei[1], Phi[j] * Ei[2]};

            double diff_Phi[3] = {diffPhi[j] * diff_ksi_eta[fam][0],
                                  diffPhi[j] * diff_ksi_eta[fam][1],
                                  diffPhi[j] * diff_ksi_eta[fam][2]};
            double diff_muVec[3] = {diff_mu[ff][0], diff_mu[ff][1],
                                    diff_mu[ff][2]};

            double diff_mu_cross_Eij[3];
            ref_hex.Cross_product(diff_muVec, Eij, diff_mu_cross_Eij);
            double diff_Phi_cross_Ei[3];
            ref_hex.Cross_product(diff_Phi, Ei, diff_Phi_cross_Ei);

            curl_faceN[ff][fam][3 * (qd_shift + n) + 0] =
                mu[ff] * diff_Phi_cross_Ei[0] + diff_mu_cross_Eij[0];
            curl_faceN[ff][fam][3 * (qd_shift + n) + 1] =
                mu[ff] * diff_Phi_cross_Ei[1] + diff_mu_cross_Eij[1];
            curl_faceN[ff][fam][3 * (qd_shift + n) + 2] =
                mu[ff] * diff_Phi_cross_Ei[2] + diff_mu_cross_Eij[2];
          }
          // for (int i = 0; i != qp[fam]; i++) {
          //   for (int j = 0; j != pq[fam] - 1; j++) {
          //     faceN[ff][fam][3 * (qd_shift + n) + 0] =
          //         mu[ff] * Phi[j] * E[i] * diff_eta_ksi[fam][0];
          //     faceN[ff][fam][3 * (qd_shift + n) + 1] =
          //         mu[ff] * Phi[j] * E[i] * diff_eta_ksi[fam][1];
          //     faceN[ff][fam][3 * (qd_shift + n) + 2] =
          //         mu[ff] * Phi[j] * E[i] * diff_eta_ksi[fam][2];

          //     double Ei[3] = {E[i] * diff_eta_ksi[fam][0],
          //                     E[i] * diff_eta_ksi[fam][1],
          //                     E[i] * diff_eta_ksi[fam][2]};

          //     double Eij[3] = {Phi[j] * Ei[0], Phi[j] * Ei[1], Phi[j] *
          //     Ei[2]};

          //     double diff_Phi[3] = {diffPhi[j] * diff_ksi_eta[fam][0],
          //                           diffPhi[j] * diff_ksi_eta[fam][1],
          //                           diffPhi[j] * diff_ksi_eta[fam][2]};
          //     double diff_muVec[3] = {diff_mu[ff][0], diff_mu[ff][1],
          //                             diff_mu[ff][2]};

          //     double diff_mu_cross_Eij[3];
          //     ref_hex.Cross_product(diff_muVec, Eij, diff_mu_cross_Eij);
          //     double diff_Phi_cross_Ei[3];
          //     ref_hex.Cross_product(diff_Phi, Ei, diff_Phi_cross_Ei);

          //     curl_faceN[ff][fam][3 * (qd_shift + n) + 0] =
          //         mu[ff] * diff_Phi_cross_Ei[0] + diff_mu_cross_Eij[0];
          //     curl_faceN[ff][fam][3 * (qd_shift + n) + 1] =
          //         mu[ff] * diff_Phi_cross_Ei[1] + diff_mu_cross_Eij[1];
          //     curl_faceN[ff][fam][3 * (qd_shift + n) + 2] =
          //         mu[ff] * diff_Phi_cross_Ei[2] + diff_mu_cross_Eij[2];
          //     ++n;
          //   }
          // }
        }
      }
    }
  }
  MoFEMFunctionReturnHot(0);
}

MoFEMErrorCode MoFEM::DemkowiczHexAndQuad::Hcurl_InteriorShapeFunctions_ONHEX(
    int *p, double *N, double *volN[3], double *curl_volN[3],
    int nb_integration_pts) {
  MoFEMFunctionBeginHot;
  RefHex ref_hex;

  for (int qq = 0; qq < nb_integration_pts; qq++) {
    // general ******************************
    int shift = 8 * qq;
    double quad_coords[3];
    double Nq[8];
    for (int vv = 0; vv < 8; vv++)
      Nq[vv] = N[shift + vv];

    ref_hex.get_volume_coords(Nq, quad_coords);
    // ****************************************
    double *ksi = quad_coords;

    double diff_ksi[3][3];
    ref_hex.get_volume_diff_coords(diff_ksi);

    double ksi_eta_gma[3] = {ksi[0], ksi[1], ksi[2]};
    double eta_gma_ksi[3] = {ksi[1], ksi[2], ksi[0]};
    double gma_ksi_eta[3] = {ksi[2], ksi[0], ksi[1]};

    double *diff_ksi_eta_gma[3] = {diff_ksi[0], diff_ksi[1], diff_ksi[2]};
    double *diff_eta_gma_ksi[3] = {diff_ksi[1], diff_ksi[2], diff_ksi[0]};
    double *diff_gma_ksi_eta[3] = {diff_ksi[2], diff_ksi[0], diff_ksi[1]};

    int pqr[3] = {p[0], p[1], p[2]};
    int qrp[3] = {p[1], p[2], p[0]};
    int rpq[3] = {p[2], p[0], p[1]};
    for (int fam = 0; fam < 3; fam++) {

      int ppp = pqr[fam];
      double PhiJ[ppp - 1];
      double diffPhiJ[ppp - 1];
      CHKERR Integrated_Legendre01(ppp, ksi_eta_gma[fam], PhiJ, diffPhiJ);

      int qqq = qrp[fam];
      double EI[qqq];
      CHKERR Legendre_polynomials01(qqq - 1, eta_gma_ksi[fam], EI);

      int rrr = rpq[fam];

      double PhiK[rrr - 1];
      double diffPhiK[rrr - 1];
      CHKERR Integrated_Legendre01(rrr, gma_ksi_eta[fam], PhiK, diffPhiK);

      int permute[(ppp - 1) * qqq * (rrr - 1)][3];
      CHKERR MonomOrdering(permute, ppp - 2, qqq - 1, rrr - 2);

      int qd_shift = (ppp - 1) * qqq * (rrr - 1) * qq;
      int n = 0;
      for (; n != (ppp - 1) * qqq * (rrr - 1); n++) {
        int ii = permute[n][0];
        int jj = permute[n][1];
        int kk = permute[n][2];

        volN[fam][3 * (qd_shift + n) + 0] =
            PhiK[kk] * PhiJ[jj] * EI[ii] * diff_eta_gma_ksi[fam][0];
        volN[fam][3 * (qd_shift + n) + 1] =
            PhiK[kk] * PhiJ[jj] * EI[ii] * diff_eta_gma_ksi[fam][1];
        volN[fam][3 * (qd_shift + n) + 2] =
            PhiK[kk] * PhiJ[jj] * EI[ii] * diff_eta_gma_ksi[fam][2];

        double EEI[3] = {EI[ii] * diff_eta_gma_ksi[fam][0],
                         EI[ii] * diff_eta_gma_ksi[fam][1],
                         EI[ii] * diff_eta_gma_ksi[fam][2]};

        double EEIJ[3] = {PhiJ[jj] * EEI[0], PhiJ[jj] * EEI[1],
                          PhiJ[jj] * EEI[2]};

        double diff_PhiJ[3] = {diffPhiJ[jj] * diff_ksi_eta_gma[fam][0],
                               diffPhiJ[jj] * diff_ksi_eta_gma[fam][1],
                               diffPhiJ[jj] * diff_ksi_eta_gma[fam][2]};

        double diff_PhiK[3] = {diffPhiK[jj] * diff_gma_ksi_eta[fam][0],
                               diffPhiJ[jj] * diff_gma_ksi_eta[fam][1],
                               diffPhiJ[jj] * diff_gma_ksi_eta[fam][2]};

        double diff_PhiK_cross_EEIJ[3];
        ref_hex.Cross_product(diff_PhiK, EEIJ, diff_PhiK_cross_EEIJ);
        double diff_PhiJ_cross_EEI[3];
        ref_hex.Cross_product(diff_PhiJ, EEI, diff_PhiJ_cross_EEI);

        curl_volN[fam][3 * (qd_shift + n) + 0] =
            PhiK[kk] * diff_PhiJ_cross_EEI[0] + diff_PhiK_cross_EEIJ[0];
        curl_volN[fam][3 * (qd_shift + n) + 1] =
            PhiK[kk] * diff_PhiJ_cross_EEI[1] + diff_PhiK_cross_EEIJ[1];
        curl_volN[fam][3 * (qd_shift + n) + 2] =
            PhiK[kk] * diff_PhiJ_cross_EEI[2] + diff_PhiK_cross_EEIJ[2];
      }

      // for (int ii = 0; ii < qqq; ii++) {
      //   for (int jj = 0; jj < ppp - 1; jj++) {
      //     for (int kk = 0; kk < rrr - 1; kk++) {
      //       volN[fam][3 * (qd_shift + n) + 0] =
      //           PhiK[kk] * PhiJ[jj] * EI[ii] * diff_eta_gma_ksi[fam][0];
      //       volN[fam][3 * (qd_shift + n) + 1] =
      //           PhiK[kk] * PhiJ[jj] * EI[ii] * diff_eta_gma_ksi[fam][1];
      //       volN[fam][3 * (qd_shift + n) + 2] =
      //           PhiK[kk] * PhiJ[jj] * EI[ii] * diff_eta_gma_ksi[fam][2];

      //       double EEI[3] = {EI[ii] * diff_eta_gma_ksi[fam][0],
      //                        EI[ii] * diff_eta_gma_ksi[fam][1],
      //                        EI[ii] * diff_eta_gma_ksi[fam][2]};

      //       double EEIJ[3] = {PhiJ[jj] * EEI[0], PhiJ[jj] * EEI[1],
      //                         PhiJ[jj] * EEI[2]};

      //       double diff_PhiJ[3] = {diffPhiJ[jj] * diff_ksi_eta_gma[fam][0],
      //                              diffPhiJ[jj] * diff_ksi_eta_gma[fam][1],
      //                              diffPhiJ[jj] * diff_ksi_eta_gma[fam][2]};

      //       double diff_PhiK[3] = {diffPhiK[jj] * diff_gma_ksi_eta[fam][0],
      //                              diffPhiJ[jj] * diff_gma_ksi_eta[fam][1],
      //                              diffPhiJ[jj] * diff_gma_ksi_eta[fam][2]};

      //       double diff_PhiK_cross_EEIJ[3];
      //       ref_hex.Cross_product(diff_PhiK, EEIJ, diff_PhiK_cross_EEIJ);
      //       double diff_PhiJ_cross_EEI[3];
      //       ref_hex.Cross_product(diff_PhiJ, EEI, diff_PhiJ_cross_EEI);

      //       curl_volN[fam][3 * (qd_shift + n) + 0] =
      //           PhiK[kk] * diff_PhiJ_cross_EEI[0] + diff_PhiK_cross_EEIJ[0];
      //       curl_volN[fam][3 * (qd_shift + n) + 1] =
      //           PhiK[kk] * diff_PhiJ_cross_EEI[1] + diff_PhiK_cross_EEIJ[1];
      //       curl_volN[fam][3 * (qd_shift + n) + 2] =
      //           PhiK[kk] * diff_PhiJ_cross_EEI[2] + diff_PhiK_cross_EEIJ[2];

      //       n++;
      //     }
      //   }
      // }
    }
  }

  MoFEMFunctionReturnHot(0);
}

MoFEMErrorCode MoFEM::DemkowiczHexAndQuad::Hdiv_FaceShapeFunctions_ONHEX(
    int *face_nodes[6], int *p, double *N, double *faceN[6],
    double *div_faceN[6], int nb_integration_pts) {
  MoFEMFunctionBeginHot;
  RefHex ref_hex;

  int pp[6][2] = {{p[0], p[1]}, {p[0], p[2]}, {p[1], p[2]},
                  {p[0], p[2]}, {p[1], p[2]}, {p[0], p[1]}};

  for (int qq = 0; qq < 6; qq++) {
    // general ******************************
    int shift = 8 * qq;
    double quad_coords[3];
    double Nq[8];
    for (int vv = 0; vv < 8; vv++)
      Nq[vv] = N[shift + vv];

    ref_hex.get_volume_coords(Nq, quad_coords);
    // ****************************************

    double mu[6];
    ref_hex.get_face_affines(quad_coords, mu);
    double diff_mu[6][3];
    ref_hex.get_face_diff_affines(diff_mu);

    double ksi[6][2];
    ref_hex.get_face_coords(face_nodes, Nq, ksi);
    double diff_ksi[6][2][3];
    ref_hex.get_face_diff_coords(face_nodes, diff_ksi);

    for (int ff = 0; ff < 6; ff++) {
      double EI[pp[ff][0]];
      CHKERR Legendre_polynomials01(pp[ff][0] - 1, ksi[ff][0], EI);

      double EJ[pp[ff][1]];
      CHKERR Legendre_polynomials01(pp[ff][1] - 1, ksi[ff][1], EJ);

      int permute[pp[ff][0] * pp[ff][1]][3];
      CHKERR MonomOrdering(permute, pp[ff][0] - 1, pp[ff][1] - 1);

      int qd_shift = pp[ff][0] * pp[ff][1] * qq;
      int n = 0;
      for (; n != pp[ff][0] * pp[ff][1]; n++) {
        int ii = permute[n][0];
        int jj = permute[n][1];

        double EEI[3] = {EI[ii] * diff_ksi[ff][0][0],
                         EI[ii] * diff_ksi[ff][0][1],
                         EI[ii] * diff_ksi[ff][0][2]};
        double EEJ[3] = {EJ[jj] * diff_ksi[ff][1][0],
                         EJ[jj] * diff_ksi[ff][1][1],
                         EJ[jj] * diff_ksi[ff][1][2]};

        double VIJ_square[3];
        ref_hex.Cross_product(EEI, EEJ, VIJ_square);

        faceN[ff][3 * (qd_shift + n) + 0] = mu[ff] * VIJ_square[0];
        faceN[ff][3 * (qd_shift + n) + 1] = mu[ff] * VIJ_square[1];
        faceN[ff][3 * (qd_shift + n) + 2] = mu[ff] * VIJ_square[2];

        div_faceN[ff][qd_shift + n] = diff_mu[ff][0] * VIJ_square[0] +
                                      diff_mu[ff][1] * VIJ_square[1] +
                                      diff_mu[ff][2] * VIJ_square[2];
      }
      // for (int ii = 0; ii < pp[ff][0]; ii++) {
      //   for (int jj = 0; jj < pp[ff][1]; jj++) {
      //     double EEI[3] = {EI[ii] * diff_ksi[ff][0][0],
      //                      EI[ii] * diff_ksi[ff][0][1],
      //                      EI[ii] * diff_ksi[ff][0][2]};
      //     double EEJ[3] = {EJ[jj] * diff_ksi[ff][1][0],
      //                      EJ[jj] * diff_ksi[ff][1][1],
      //                      EJ[jj] * diff_ksi[ff][1][2]};

      //     double VIJ_square[3];
      //     ref_hex.Cross_product(EEI, EEJ, VIJ_square);

      //     faceN[ff][3 * (qd_shift + n) + 0] = mu[ff] * VIJ_square[0];
      //     faceN[ff][3 * (qd_shift + n) + 1] = mu[ff] * VIJ_square[1];
      //     faceN[ff][3 * (qd_shift + n) + 2] = mu[ff] * VIJ_square[2];

      //     div_faceN[ff][qd_shift + n] = diff_mu[ff][0] * VIJ_square[0] +
      //                                   diff_mu[ff][1] * VIJ_square[1] +
      //                                   diff_mu[ff][2] * VIJ_square[2];

      //     n++;
      //   }
      // }
    }
  }
  MoFEMFunctionReturnHot(0);
}

MoFEMErrorCode MoFEM::DemkowiczHexAndQuad::Hdiv_InteriorShapeFunctions_ONHEX(
    int *p, double *N, double *bubbleN[3], double *div_bubbleN[3],
    int nb_integration_pts) {
  MoFEMFunctionBeginHot;
  RefHex ref_hex;

  for (int qq = 0; qq < nb_integration_pts; qq++) {
    // general ******************************
    int shift = 8 * qq;
    double quad_coords[3];
    double Nq[8];
    for (int vv = 0; vv < 8; vv++)
      Nq[vv] = N[shift + vv];

    ref_hex.get_volume_coords(Nq, quad_coords);
    // ****************************************
    double *ksi = quad_coords;
    double diff_ksi[3][3];

    ref_hex.get_volume_diff_coords(diff_ksi);

    double ksi_eta_gma[3] = {ksi[0], ksi[1], ksi[2]};
    double eta_gma_ksi[3] = {ksi[1], ksi[2], ksi[0]};
    double gma_ksi_eta[3] = {ksi[2], ksi[0], ksi[1]};

    double *diff_ksi_eta_gma[3] = {diff_ksi[0], diff_ksi[1], diff_ksi[2]};
    double *diff_eta_gma_ksi[3] = {diff_ksi[1], diff_ksi[2], diff_ksi[0]};
    double *diff_gma_ksi_eta[3] = {diff_ksi[2], diff_ksi[0], diff_ksi[1]};

    int pqr[3] = {p[0], p[1], p[2]};
    int qrp[3] = {p[1], p[2], p[0]};
    int rpq[3] = {p[2], p[0], p[1]};

    for (int fam = 0; fam < 3; fam++) {
      double PhiK[pqr[fam] - 1];
      double diffPhiK[pqr[fam] - 1];
      CHKERR Integrated_Legendre01(pqr[fam], ksi_eta_gma[fam], PhiK, diffPhiK);

      double EI[qrp[fam]];
      CHKERR Legendre_polynomials01(qrp[fam] - 1, eta_gma_ksi[fam], EI);

      double EJ[rpq[fam]];
      CHKERR Legendre_polynomials01(rpq[fam] - 1, gma_ksi_eta[fam], EJ);

      int permute[(pqr[fam] - 1) * qrp[fam] * rpq[fam]][3];
      CHKERR MonomOrdering(permute, qrp[fam] - 1, rpq[fam] - 1, pqr[fam] - 2);

      int qd_shift = (pqr[fam] - 1) * qrp[fam] * rpq[fam] * qq;
      int n = 0;
      for (; n != (pqr[fam] - 1) * qrp[fam] * rpq[fam]; n++) {
        int ii = permute[n][0];
        int jj = permute[n][1];
        int kk = permute[n][2];

        double EEI[3] = {EI[ii] * diff_eta_gma_ksi[fam][0],
                         EI[ii] * diff_eta_gma_ksi[fam][1],
                         EI[ii] * diff_eta_gma_ksi[fam][2]};
        double EEJ[3] = {EJ[jj] * diff_gma_ksi_eta[fam][0],
                         EJ[jj] * diff_gma_ksi_eta[fam][1],
                         EJ[jj] * diff_gma_ksi_eta[fam][2]};

        double VIJ_square[3];
        ref_hex.Cross_product(EEI, EEJ, VIJ_square);

        bubbleN[fam][3 * (qd_shift + n) + 0] = PhiK[kk] * VIJ_square[0];
        bubbleN[fam][3 * (qd_shift + n) + 1] = PhiK[kk] * VIJ_square[1];
        bubbleN[fam][3 * (qd_shift + n) + 2] = PhiK[kk] * VIJ_square[2];

        div_bubbleN[fam][qd_shift + n] =
            diffPhiK[kk] * diff_ksi_eta_gma[fam][0] * VIJ_square[0] +
            diffPhiK[kk] * diff_ksi_eta_gma[fam][1] * VIJ_square[1] +
            diffPhiK[kk] * diff_ksi_eta_gma[fam][2] * VIJ_square[2];
      }
      // for (int ii = 0; ii < qrp[fam]; ii++) {
      //   for (int jj = 0; jj < rpq[fam]; jj++) {
      //     for (int kk = 0; kk < pqr[fam] - 1; kk++) {
      //       double EEI[3] = {EI[ii] * diff_eta_gma_ksi[fam][0],
      //                        EI[ii] * diff_eta_gma_ksi[fam][1],
      //                        EI[ii] * diff_eta_gma_ksi[fam][2]};
      //       double EEJ[3] = {EJ[jj] * diff_gma_ksi_eta[fam][0],
      //                        EJ[jj] * diff_gma_ksi_eta[fam][1],
      //                        EJ[jj] * diff_gma_ksi_eta[fam][2]};

      //       double VIJ_square[3];
      //       ref_hex.Cross_product(EEI, EEJ, VIJ_square);

      //       bubbleN[fam][3 * (qd_shift + n) + 0] = PhiK[kk] * VIJ_square[0];
      //       bubbleN[fam][3 * (qd_shift + n) + 1] = PhiK[kk] * VIJ_square[1];
      //       bubbleN[fam][3 * (qd_shift + n) + 2] = PhiK[kk] * VIJ_square[2];

      //       div_bubbleN[fam][qd_shift + n] =
      //           diffPhiK[kk] * diff_ksi_eta_gma[fam][0] * VIJ_square[0] +
      //           diffPhiK[kk] * diff_ksi_eta_gma[fam][1] * VIJ_square[1] +
      //           diffPhiK[kk] * diff_ksi_eta_gma[fam][2] * VIJ_square[2];
      //       n++;
      //     }
      //   }
      // }
    }
  }

  MoFEMFunctionReturnHot(0);
}
