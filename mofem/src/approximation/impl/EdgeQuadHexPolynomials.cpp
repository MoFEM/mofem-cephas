/** \file Hdiv.cpp

  \brief Implementation of hierarchical Edge, Quad, and Hex shape functions of 
  type H1, Hcurl, Hdiv, Hdiv
*/

using namespace MoFEM;

MoFEMErrorCode MoFEM::Integrated_Legendre(int p, double s, double *L,
                                          double *diffL) {

  MoFEMFunctionBeginHot;

  double l[p + 2];
  CHKERR Legendre_polynomials(p + 1, s, NULL, l, NULL, 1);
  if (p >= 0)
  {
    L[0] = 0.5*(s + 1.0);
    diffL[0] = 0.5 * l[0];
    if (p >= 1)
    {
      for (int i = 0; i != p; i++)
      {
        double factor = 1.0 / (2.0 * (2.0*(i+1) + 1.0));

        L[i+1] = factor * (l[i+2] - factor * l[i]);
        diffL[i+1] = 0.5 * l[i+1];
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

MoFEMErrorCode MoFEM::H1_EdgeShapeFunctions_ONQUAD(int *sense, int *p,
                                                   double *N, double *edgeN[],
                                                   double *diff_edgeN[],
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

    // constant affine coordinates of edges per cannonical numbering
    double mu[4] = {mu_eta0, mu_ksi1, mu_eta1, mu_ksi0};
    double diff_mu[4][2] = {{0, -0.5}, {0.5, 0.0}, {0.0, 0.5}, {-0.5, 0.0}};

    // parametrization of each edge per cannonical numbering
    double s[4] = {ksi, eta, ksi, eta};
    double diff_s[4][2] = {{1.0, 0.0}, {0.0, 1.0}, {1.0, 0.0}, {0.0, 1.0}};

    for (int e = 0; e != 4; e++) {
      double L[p[e] + 1];
      double diffL[p[e] + 1];
      double ss = s[e] * sense[e];
      CHKERR Integrated_Legendre(p[e], ss, L, diffL);
      int qd_shift = (p[e]+1) * q;
      for (int n = 0; n != p[e] + 1; n++) {
        edgeN[qd_shift + n][e] = mu[e] * L[n];
        diff_edgeN[2 * qd_shift + 2 * n][e] =
            mu[e] * diffL[n] * sense[e] * diff_s[e][0] + diff_mu[e][0] * L[n];
        diff_edgeN[2 * qd_shift + 2 * n + 1][e] =
            mu[e] * diffL[n] * sense[e] * diff_s[e][1] + diff_mu[e][1] * L[n];
      }
    }
  }

  MoFEMFunctionReturnHot(0);
}

MoFEMErrorCode MoFEM::H1_FaceShapeFunctions_ONQUAD(int *p,
                                                   double *N, double *edgeN[],
                                                   double *diff_edgeN[],
                                                   int nb_integration_pts) {
    MoFEMFunctionBeginHot;
    for (int q = 0; q != nb_integration_pts; q++) {
      int shift = 4 * q;
      double ksi = -N[shift + 0] + N[shift + 1] + N[shift + 2] - N[shift + 3];
      double eta = -N[shift + 0] - N[shift + 1] + N[shift + 2] + N[shift + 3];

    }

    MoFEMFunctionReturnHot(0);
  }
