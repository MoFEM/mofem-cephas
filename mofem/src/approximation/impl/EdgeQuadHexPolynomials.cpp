/** \file EdgeQuadHexPolynomials.cpp

  \brief Implementation of hierarchical Edge, Quad, and Hex shape functions of 
  type H1, Hcurl, Hdiv, Hdiv
*/

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

using namespace MoFEM;

MoFEMErrorCode MoFEM::Integrated_Legendre(int p, double s, double *L,
                                          double *diffL) {

  MoFEMFunctionBeginHot;

  double l[p + 1];
  CHKERR Legendre_polynomials(p, s, NULL, l, NULL, 1);
  if (p >= 2)
  {
    for (int i = 0; i != p-1; i++)
    {
      double factor = 1.0 / ((2.0 * (i+1) + 1.0));

      L[i] = factor * (l[i+2] - l[i]);
      diffL[i] = l[i+1];
    }
  }
  
  MoFEMFunctionReturnHot(0);
}

MoFEMErrorCode MoFEM::H1_EdgeShapeFunctions_ONQUAD(int *sense, int *p,
                                                   double *N, double
                                                   *edgeN[4], double
                                                   *diff_edgeN[4], int
                                                   nb_integration_pts) {
  MoFEMFunctionBeginHot;
    for (int q = 0; q != nb_integration_pts; q++) {
      int shift = 4 * q;
      double ksi = -N[shift + 0] + N[shift + 1] + N[shift + 2] - N[shift + 3]; 
      double eta = -N[shift + 0] - N[shift + 1] + N[shift + 2] + N[shift+ 3];

      // Affine coordinates of each eadge mu0 and mu1
      double mu_ksi0 = 1.0 - 0.5 * (ksi + 1.0);
      double mu_eta0 = 1.0 - 0.5 * (eta + 1.0);

      double mu_ksi1 = 0.5 * (ksi + 1.0);
      double mu_eta1 = 0.5 * (eta + 1.0);

      // constant affine coordinates of edges per cannonical numbering
      double mu[4] = {mu_eta0, mu_ksi1, mu_eta1, mu_ksi0};
      double diff_mu[4][2] = {{0.0, -0.5}, {0.5, 0.0}, {0.0, 0.5}, {-0.5,
      0.0}};

      // parametrization of each edge per cannonical numbering
      double s[4] = {ksi, eta, ksi, eta};
      double diff_s[4][2] = {{1.0, 0.0}, {0.0, 1.0}, {1.0, 0.0}, {0.0, 1.0}};

      for (int e = 0; e != 4; e++) {
        double L[p[e] - 1];
        double diffL[p[e] - 1];
        double ss = s[e] * sense[e];
        CHKERR Integrated_Legendre(p[e], ss, L, diffL);
        int qd_shift = (p[e]- 1) * q;
        for (int n = 0; n != p[e]-1; n++) {
          edgeN[e][qd_shift + n] = mu[e] * L[n];
          diff_edgeN[e][2 * qd_shift + 2 * n] =
              mu[e] * diffL[n] * sense[e] * diff_s[e][0] + diff_mu[e][0] *
              L[n];
          diff_edgeN[e][2 * qd_shift + 2 * n + 1] =
              mu[e] * diffL[n] * sense[e] * diff_s[e][1] + diff_mu[e][1] *
              L[n];
        }
      }
    }
  MoFEMFunctionReturnHot(0);
}

MoFEMErrorCode MoFEM::H1_FaceShapeFunctions_ONQUAD(int *p, double *N,
                                                   double *faceN,
                                                   double *diff_faceN,
                                                   int nb_integration_pts) {
  MoFEMFunctionBeginHot;
  for (int q = 0; q != nb_integration_pts; q++) {
    int shift = 4 * q;
    double ksi = -N[shift + 0] + N[shift + 1] + N[shift + 2] - N[shift + 3];
    double eta = -N[shift + 0] - N[shift + 1] + N[shift + 2] + N[shift + 3];

    double mu_ksi0 = 1.0 - 0.5 * (ksi + 1.0);
    double mu_eta0 = 1.0 - 0.5 * (eta + 1.0);

    double mu_ksi1 = 0.5 * (ksi + 1.0);
    double mu_eta1 = 0.5 * (eta + 1.0);
    
    double s[2] = {ksi, eta};
    double diff_s[2][2] = {{1.0, 0.0}, {0.0, 1.0}};
    double L0[p[0]-1];   double diffL0[p[0]-1];
    double L1[p[1]-1];   double diffL1[p[1]-1];
    CHKERR Integrated_Legendre(p[0], ksi, L0, diffL0);
    CHKERR Integrated_Legendre(p[1], eta, L1, diffL1);

    int qd_shift = (p[0] - 1) * (p[1] - 1) * q;
    int n = 0;
    for (int s1 = 0; s1 != p[0]-1; s1++){
      for (int s2 = 0; s2 != p[1]-1; s2++){
        faceN[qd_shift + n] = L0[s1] * L1[s2];
        diff_faceN[2 * qd_shift + 2 * n] = diffL0[s1] * L1[s2]*diff_s[0][0] + L0[s1] * diffL1[s2] * diff_s[1][0];
        diff_faceN[2 * qd_shift + 2 * n + 1] = diffL0[s1] * L1[s2] * diff_s[0][1] + L0[s1] * diffL1[s2] * diff_s[1][1];
        ++n;
      }
    }
  }
  MoFEMFunctionReturnHot(0);
}

MoFEMErrorCode MoFEM::L2_FaceShapeFunctions_ONQUAD(int *p, double *N, double *faceN,
                                                   double *diff_faceN, int nb_integration_pts) {
  MoFEMFunctionBeginHot;
  for (int q = 0; q != nb_integration_pts; q++) {
    int shift = 4 * q;
    double ksi = -N[shift + 0] + N[shift + 1] + N[shift + 2] - N[shift + 3];
    double eta = -N[shift + 0] - N[shift + 1] + N[shift + 2] + N[shift + 3];

    double L0[p[0] - 1];
    double L1[p[1] - 1];

    CHKERR Legendre_polynomials(p[0], ksi, NULL, L0, NULL, 1);
    CHKERR Legendre_polynomials(p[1], eta, NULL, L1, NULL, 1);

    int qd_shift = (p[0] - 1) * (p[1] - 1)  * q;
    int n = 0;
    for (int s1 = 0; s1 != p[0]; s1++) {
      for (int s2 = 0; s2 != p[1]; s2++) {
        faceN[qd_shift + n] = L0[s1] * L1[s2];
        ++n;
      }
    }
  }
  MoFEMFunctionReturnHot(0);
}

MoFEMErrorCode MoFEM::Hcurl_EdgeShapeFunctions_ONQUAD(int *sense, int *p,
                                                      double *N,
                                                      double *edgeN[],
                                                      double *curl_edgeN[],
                                                      int nb_integration_pts) {
  MoFEMFunctionBeginHot;
  for (int q = 0; q != nb_integration_pts; q++) {
    int shift = 4 * q;
    double ksi = -N[shift + 0] + N[shift + 1] + N[shift + 2] - N[shift + 3];
    double eta = -N[shift + 0] - N[shift + 1] + N[shift + 2] + N[shift + 3];

    // Affine coordinates of each eadge mu0 and mu1
    double mu_ksi0 = 1.0 - 0.5 * (ksi + 1.0);
    double mu_eta0 = 1.0 - 0.5 * (eta + 1.0);

    double mu_ksi1 = 0.5 * (ksi + 1.0);
    double mu_eta1 = 0.5 * (eta + 1.0);

    // affine coordinates parametrising of edges per cannonical numbering
    double mu[4] = {mu_eta0, mu_ksi1, mu_eta1, mu_ksi0};
    double diff_mu[4][2] = {{0.0, -0.5}, {0.5, 0.0}, {0.0, 0.5}, {-0.5, 0.0}};

    // parametrization of each edge per cannonical numbering
    double s[4] = {ksi, eta, ksi, eta};
    double diff_s[4][2] = {{1.0, 0.0}, {0.0, 1.0}, {1.0, 0.0}, {0.0, 1.0}};

    for (int e = 0; e != 4; e++) {
      double L[p[e] - 1];
      double ss = s[e] * sense[e];
      CHKERR Legendre_polynomials(p[0], ss, NULL, L, NULL, 1);
      int qd_shift = (p[e] - 1) * q;
      for (int n = 0; n != p[e]; n++) {
        edgeN[e][2 * (qd_shift + n) + 0] = mu[e] * L[n] * diff_s[e][0];
        edgeN[e][2 * (qd_shift + n) + 1] = mu[e] * L[n] * diff_s[e][1];

        double E1[2] = {diff_mu[e][0], diff_mu[e][1]};
        double E2[2] = {L[n] * diff_s[e][0], L[n] * diff_s[e][0]};
        curl_edgeN[e][qd_shift + n] = E1[0] * E2[1] - E1[1] * E2[0];
      }
    }
  }
  MoFEMFunctionReturnHot(0);
}

MoFEMErrorCode Hcurl_FaceShapeFunctions_ONQUAD(int *p, double *N,
                                               double *faceN[],
                                               double *curl_faceN[],
                                               int nb_integration_pts){
  MoFEMFunctionBeginHot;
  for (int q = 0; q != nb_integration_pts; q++) {
    int shift = 4 * q;

    double ksi = -N[shift + 0] + N[shift + 1] + N[shift + 2] - N[shift + 3];
    double eta = -N[shift + 0] - N[shift + 1] + N[shift + 2] + N[shift + 3];

    double s[2] = {ksi, eta};
    double diff_s[2][2] = {{1.0, 0.0}, {0.0, 1.0}};

    int sgn[2] = {-1, 1};

    for (int typ = 0; typ != 2; typ++)
    {
      int pp = p[typ];
      int qq = p[(typ + 1) % 2];

      double Phi[pp - 1];
      double diffPhi[pp - 1];
      CHKERR Integrated_Legendre(pp, s[typ], Phi, diffPhi);

      double E[qq];
      CHKERR Legendre_polynomials(qq - 1, s[typ], NULL, E, NULL, 1);

      int qd_shift = (pp - 1) * qq * q;
      int n = 0;
      for (int i = 0; i != qq - 1; i++) {
        for (int j = 0; j != pp - 2; i++) {
          faceN[typ][2 * (qd_shift + n) + 0] = Phi[j] * E[i] * diff_s[typ][0];
          faceN[typ][2 * (qd_shift + n) + 1] = Phi[j] * E[i] * diff_s[typ][1];
          curl_faceN[typ][qd_shift + n] = (double)sgn[typ] * diffPhi[j] * E[i];
          ++n;
        }
      }
    }
  }
  MoFEMFunctionReturnHot(0);
}
