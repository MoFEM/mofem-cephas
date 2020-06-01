/** \file EdgeQuadHexPolynomials.cpp

  \brief Implementation of hierarchical Edge, Quad, and Hex shape functions of 
  type H1, Hcurl, Hdiv
*/

/*
      Quads
 4-------3------3
 |              |       eta
 |              |       ^
 4              2       |
 |              |       |
 |              |       0-----  > ksi
 1-------1------2
*/

using namespace MoFEM;

// Auxilary functions

MoFEMErrorCode MoFEM::Integrated_Legendre(int p, double s, double *L,
                                          double *diffL) {

  MoFEMFunctionBeginHot;

  double l[p + 1];
  CHKERR Legendre_polynomials(p, s, NULL, l, NULL, 1);
  if (p >= 2) {
    for (int i = 0; i != p - 1; i++) {
      double factor = 1.0 / ((2.0 * (i + 1) + 1.0));

      L[i] = factor * (l[i + 2] - l[i]);
      diffL[i] = l[i + 1];
    }
  }

  MoFEMFunctionReturnHot(0);
}

MoFEMErrorCode MoFEM::Face_orientMat(int *face_nodes, double orientMat[2][2]) {

  MoFEMFunctionBeginHot;

  for (int i = 0; i != 2; i++)
    for (int j = 0; j != 2; j++)
      orientMat[i][j] = 0.0;

  int node0 = face_nodes[0]; int node1 = face_nodes[1];

  switch (node0) {
  case 0 : 
    switch (node1){
      case 1 :  //positive direction
        orientMat[0][0] = 1.0;
        orientMat[1][1] = 1.0; //
        break;
      case 3 :  // negative direction
        orientMat[0][1] = 1.0;
        orientMat[1][0] = 1.0; //
        break;
    }
    break;
  case 1 : 
    switch (node1){
    case 2 :  //positive direction
      orientMat[0][1] = 1.0;
      orientMat[1][0] = -1.0; //
      break;
      case 0 :  // negative direction
        orientMat[0][0] = -1.0;
        orientMat[1][1] = 1.0; //
        break;
    }
    break;
  case 2 : 
    switch (node1){
      case 3 :
        orientMat[0][0] = -1.0;
        orientMat[1][1] = -1.0;
        break;
      case 1 :
        orientMat[0][1] = -1.0;
        orientMat[1][0] = -1.0;
        break; 
    }
    break;
    case 3 :
      switch (node1){
        case 0 :
          orientMat[0][1] = -1.0;
          orientMat[1][0] = 1.0;
          break;
        case 2 :
          orientMat[0][0] = 1.0;
          orientMat[1][1] = -1.0;
          break;
      }
      break;
    default :
      break;
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
      double eta = -N[shift + 0] - N[shift + 1] + N[shift + 2] + N[shift + 3];

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

    // cout << "In Face H1" << endl;

    int qd_shift = (p[0] - 1) * (p[1] - 1) * q;
    int n = 0;
    for (int s1 = 0; s1 != p[0]-1; s1++){
      for (int s2 = 0; s2 != p[1]-1; s2++){
        faceN[qd_shift + n] = L0[s1] * L1[s2];
        diff_faceN[2 * (qd_shift + n) + 0] = diffL0[s1] * L1[s2]*diff_s[0][0] + L0[s1] * diffL1[s2] * diff_s[1][0];
        diff_faceN[2 * (qd_shift + n) + 1] = diffL0[s1] * L1[s2] * diff_s[0][1] + L0[s1] * diffL1[s2] * diff_s[1][1];
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

    double L0[p[0]];
    double L1[p[1]];

    CHKERR Legendre_polynomials(p[0] - 1, ksi, NULL, L0, NULL, 1);
    CHKERR Legendre_polynomials(p[1] - 1, eta, NULL, L1, NULL, 1);

    int qd_shift = p[0] * p[1]  * q;
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
                                                      double *edgeN[4],
                                                      double *curl_edgeN[4],
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
    int pp[4] = {p[0], p[1], p[0], p[1]};
    for (int e = 0; e != 4; e++) {
      double L[pp[e]];
      double ss = s[e] * sense[e];
      CHKERR Legendre_polynomials(pp[e] - 1, ss, NULL, L, NULL, 1);
      int qd_shift = pp[e] * q;
      for (int n = 0; n != pp[e]; n++) {
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

MoFEMErrorCode MoFEM::Hcurl_FaceShapeFunctions_ONQUAD(int *p, double *N,
                                                      double *faceN[2],
                                                      double *curl_faceN[2],
                                                      int nb_integration_pts) {
  MoFEMFunctionBeginHot;
  for (int q = 0; q != nb_integration_pts; q++) {
    int shift = 4 * q;

    double ksi = -N[shift + 0] + N[shift + 1] + N[shift + 2] - N[shift + 3];
    double eta = -N[shift + 0] - N[shift + 1] + N[shift + 2] + N[shift + 3];
    double s[2][2] = {{ksi, eta}, {eta, ksi}};
    double diff_s[2][2] = {{1.0, 0.0}, {0.0, 1.0}};

    double sgn[2] = {-1.0, 1.0};

    int pq[2] = {p[0], p[1]};
    int qp[2] = {p[1], p[0]};
    for (int typ = 0; typ != 2; typ++)
    {
      int pp = pq[typ];
      int qq = qp[typ];

      

      double Phi[pp - 1];
      double diffPhi[pp - 1];
      CHKERR Integrated_Legendre(pp, s[typ][1], Phi, diffPhi);

      

      double E[qq];
      CHKERR Legendre_polynomials(qq - 1, s[typ][0], NULL, E, NULL, 1);

      int qd_shift = (pp - 1) * qq * q;
      int n = 0;
      for (int i = 0; i != qq; i++) {
        for (int j = 0; j != pp - 1; j++) {
          faceN[typ][2 * (qd_shift + n) + 0] = Phi[j] * E[i] * diff_s[typ][0];
          faceN[typ][2 * (qd_shift + n) + 1] = Phi[j] * E[i] * diff_s[typ][1];
          curl_faceN[typ][qd_shift + n] = sgn[typ] * diffPhi[j] * E[i];
          ++n;
        }
      }
    }
  }
  MoFEMFunctionReturnHot(0);
}

MoFEMErrorCode MoFEM::Hdiv_EdgeShapeFunctions_ONQUAD(int *sense, int *p,
                                                     double *N, double *edgeN[],
                                                     double *div_edgeN[],
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
    int pp[4] = {p[0], p[1], p[0], p[1]};
    for (int e = 0; e != 4; e++) {
      double L[pp[e]];
      double ss = s[e] * sense[e];
      CHKERR Legendre_polynomials(pp[e] - 1, ss, NULL, L, NULL, 1);
      int qd_shift = pp[e] * q;
      for (int n = 0; n != pp[e]; n++) {
        edgeN[e][2 * (qd_shift + n) + 0] = mu[e] * L[n] * diff_s[e][1];
        edgeN[e][2 * (qd_shift + n) + 1] = -mu[e] * L[n] * diff_s[e][0];

        double E1[2] = {diff_mu[e][0], diff_mu[e][1]};
        double E2[2] = {L[n] * diff_s[e][0], L[n] * diff_s[e][0]};
        div_edgeN[e][qd_shift + n] = E1[0] * E2[1] - E1[1] * E2[0];
      }
    }
  }
  MoFEMFunctionReturnHot(0);
}

MoFEMErrorCode MoFEM::Hdiv_FaceShapeFunctions_ONQUAD(int *p, double *N,
                                                     double *faceN[],
                                                     double *div_faceN[],
                                                     int nb_integration_pts) {
  MoFEMFunctionBeginHot;
  for (int q = 0; q != nb_integration_pts; q++) {
    int shift = 4 * q;

    double ksi = -N[shift + 0] + N[shift + 1] + N[shift + 2] - N[shift + 3];
    double eta = -N[shift + 0] - N[shift + 1] + N[shift + 2] + N[shift + 3];
    double s[2][2] = {{ksi, eta}, {eta, ksi}};
    double diff_s[2][2] = {{1.0, 0.0}, {0.0, 1.0}};

    double sgn[2] = {-1.0, 1.0};

    int pq[2] = {p[0], p[1]};
    int qp[2] = {p[1], p[0]};
    for (int typ = 0; typ != 2; typ++) {
      int pp = pq[typ];
      int qq = qp[typ];

      double Phi[pp - 1];
      double diffPhi[pp - 1];
      CHKERR Integrated_Legendre(pp, s[typ][1], Phi, diffPhi);

      double E[qq];
      CHKERR Legendre_polynomials(qq - 1, s[typ][0], NULL, E, NULL, 1);

      int qd_shift = (pp - 1) * qq * q;
      int n = 0;
      for (int i = 0; i != qq; i++) {
        for (int j = 0; j != pp - 1; j++) {
          faceN[typ][2 * (qd_shift + n) + 0] = Phi[j] * E[i] * diff_s[typ][1];
          faceN[typ][2 * (qd_shift + n) + 1] = -Phi[j] * E[i] * diff_s[typ][0];
          div_faceN[typ][qd_shift + n] = sgn[typ] * diffPhi[j] * E[i];
          ++n;
        }
      }
    }
  }
  MoFEMFunctionReturnHot(0);
}
/* Reference Hex and its cannonical vertex and edge numbering
               8 ---------11--------- 7 
              /|                    /|
             / |                   / |
           12  |                 10  |             x3
           /   8                 /   |             |     
          /    |                /    7             |    x2
        5 ----------9--------- 6     |             |   /
         |     |               |     |             |  /
         |    4 ----------3---------- 3            | /
        5|    /                |    /              o -------- x1
         |   /                 6   /
         |  4                  |  2
         | /                   | / 
         |/                    |/
        1 ---------1---------- 2

  Hex Face Cannonical numbering

        1. 1 2 3 4
        2. 1 2 6 5
        3. 2 3 7 6
        4. 4 3 7 8
        5. 1 4 8 5
        6. 5 6 7 8
*/

MoFEMErrorCode MoFEM::H1_EdgeShapeFunctions_ONHEX(int        *sense, 
                                                  int         *p,   
                                                  double      *N,
                                                  double      *edgeN[12],
                                                  double      *diff_edgeN[12],
                                                   int         nb_integration_pts){
  MoFEMFunctionBeginHot;
  double vertices[8][3] =  {{-1., -1., -1.}, 
                            { 1., -1., -1.},
                            { 1.,  1., -1.},
                            {-1.,  1., -1.},
                            {-1., -1.,  1.},
                            { 1., -1.,  1.},
                            { 1.,  1.,  1.},
                            {-1.,  1.,  1.}};
  cout << "++++++++++++++++++++++++++++++++++++" << endl;
  for (int q = 0; q != nb_integration_pts; q++) {
    
    int shift = 8 * q;
    double ksi = 0.; double eta = 0.; double gma = 0.;
    for (int vv = 0; vv != 8; vv++)
    {
      ksi += N[shift + vv] * vertices[vv][0];
      eta += N[shift + vv] * vertices[vv][1];
      gma += N[shift + vv] * vertices[vv][2];
    }

    cout << "ksi : " << ksi << "\t";
    cout << "eta : " << eta << "\t";
    cout << "gma : " << gma << endl;

    double diff_ksi[3] = {1.0, 0.0, 0.0};
    double diff_eta[3] = {0.0, 1.0, 0.0};
    double diff_gma[3] = {0.0, 0.0, 1.0};

    // Affine coordinates of each eadge mu0 and mu1
    double mu_ksi0 = 1.0 - 0.5 * (ksi + 1.0); double diff_mu_ksi0[3] = {-0.5, 0.0, 0.0};
    double mu_eta0 = 1.0 - 0.5 * (eta + 1.0); double diff_mu_eta0[3] = {0.0, -0.5, 0.0};
    double mu_gma0 = 1.0 - 0.5 * (gma + 1.0); double diff_mu_gma0[3] = {0.0, 0.0, -0.5};

    double mu_ksi1 = 0.5 * (ksi + 1.0);  double diff_mu_ksi1[3] = {0.5, 0.0, 0.0};
    double mu_eta1 = 0.5 * (eta + 1.0);  double diff_mu_eta1[3] = {0.0, 0.5, 0.0};
    double mu_gma1 = 0.5 * (gma + 1.0);  double diff_mu_gma1[3] = {0.0, 0.0, 0.5};

    // constant affine coordinates of edges per cannonical numbering
    double mu[12][2] = {{mu_eta0, mu_gma0}, {mu_ksi1, mu_gma0}, {mu_eta1, mu_gma0}, {mu_ksi0, mu_gma0},
                        {mu_ksi0, mu_eta0}, {mu_ksi1, mu_eta0}, {mu_ksi1, mu_eta1}, {mu_ksi0, mu_eta1},
                        {mu_eta0, mu_gma1}, {mu_ksi1, mu_gma1}, {mu_eta1, mu_gma1}, {mu_ksi0, mu_gma1}};

    double *diff_mu[12][2] = {{diff_mu_eta0, diff_mu_gma0}, {diff_mu_ksi1, diff_mu_gma0}, {diff_mu_eta1, diff_mu_gma0}, {diff_mu_ksi0, diff_mu_gma0},
                              {diff_mu_ksi0, diff_mu_eta0}, {diff_mu_ksi1, diff_mu_eta0}, {diff_mu_ksi1, diff_mu_eta1}, {diff_mu_ksi0, diff_mu_eta1},
                              {diff_mu_eta0, diff_mu_gma1}, {diff_mu_ksi1, diff_mu_gma1}, {diff_mu_eta1, diff_mu_gma1}, {diff_mu_ksi0, diff_mu_gma1}};
                        
    

    // parametrization of each edge per cannonical numbering
    double s[12] = {ksi, eta, ksi, eta, 
                    gma, gma, gma, gma, 
                    ksi, eta, ksi, eta};
    double *diff_s[12] = {diff_ksi,  diff_eta, diff_ksi, diff_eta,
                            diff_gma,  diff_gma, diff_gma, diff_gma,
                            diff_ksi,  diff_eta, diff_ksi, diff_eta};

    int pp[12] = {p[0], p[1], p[0], p[1], p[2], p[2], p[2], p[2], p[0], p[1], p[0], p[1]};

    

    for (int e = 0; e != 12; e++) {
      double L[pp[e] - 1];
      double diffL[pp[e] - 1];
      double ss = s[e] * sense[e];
      
      CHKERR Integrated_Legendre(pp[e], ss, L, diffL);
      int qd_shift = (pp[e] - 1) * q;

      for (int n = 0; n != pp[e] - 1; n++) {
        edgeN[e][qd_shift + n] = mu[e][0] * mu[e][1] * L[n];
        
        diff_edgeN[e][3 * (qd_shift + n) + 0] = mu[e][0] * mu[e][1] * diffL[n] * sense[e] * diff_s[e][0] + 
                                                (mu[e][0] * diff_mu[e][1][0] + diff_mu[e][0][0] * mu[e][1]) * L[n];

        diff_edgeN[e][3 * (qd_shift + n) + 1] = mu[e][0] * mu[e][1] * diffL[n] * sense[e] * diff_s[e][1] + 
                                                (mu[e][0] * diff_mu[e][1][1] + diff_mu[e][0][1] * mu[e][1]) * L[n];

        diff_edgeN[e][3 * (qd_shift + n) + 2] = mu[e][0] * mu[e][1] * diffL[n] * sense[e] * diff_s[e][2] + 
                                                (mu[e][0] * diff_mu[e][1][2] + diff_mu[e][0][2] * mu[e][1]) * L[n];

      }
    }
  }
  MoFEMFunctionReturnHot(0);
}
// TODO
MoFEMErrorCode MoFEM::H1_FaceShapeFunctions_ONHEX(int *face_nodes[6], 
                                                  int *p, double *N,
                                                   double *faceN[6],
                                                   double *diff_faceN[6],
                                                   int nb_integration_pts) {
  MoFEMFunctionBeginHot;

  double vertices[8][3] = {{-1., -1., -1.}, 
                           { 1., -1., -1.}, 
                           { 1.,  1., -1.}, 
                           {-1.,  1., -1.},  
                           {-1., -1.,  1.}, 
                           { 1., -1.,  1.},
                           { 1.,  1.,  1.}, 
                           {-1.,  1.,  1.}};
  for (int q = 0; q != nb_integration_pts; q++) {
    
    int shift = 8 * q;
    double ksi = 0.;
    double eta = 0.;
    double gma = 0.;
    
    for (int vv = 0; vv != 8; vv++) {
      ksi += N[shift + vv] * vertices[vv][0];
      eta += N[shift + vv] * vertices[vv][1];
      gma += N[shift + vv] * vertices[vv][2];
    }

    double diff_ksi[3] = {1.0, 0.0, 0.0};
    double diff_eta[3] = {0.0, 1.0, 0.0};
    double diff_gma[3] = {0.0, 0.0, 1.0};

    // Affine coordinates of each eadge mu0 and mu1
    double mu_ksi0 = 1.0 - 0.5 * (ksi + 1.0); double diff_mu_ksi0[3] = {-0.5, 0.0, 0.0};
    double mu_eta0 = 1.0 - 0.5 * (eta + 1.0); double diff_mu_eta0[3] = {0.5, -0.5, 0.0};
    double mu_gma0 = 1.0 - 0.5 * (gma + 1.0); double diff_mu_gma0[3] = {0.5, 0.0, -0.5};

    double mu_ksi1 = 0.5 * (ksi + 1.0);  double diff_mu_ksi1[3] = {0.5, 0.0, 0.0};
    double mu_eta1 = 0.5 * (eta + 1.0);  double diff_mu_eta1[3] = {0.5, 0.5, 0.0};
    double mu_gma1 = 0.5 * (gma + 1.0);  double diff_mu_gma1[3] = {0.5, 0.0, 0.5};

    // constant affine, face normal coordinates of each face per cannonical numbering
    double mu[6] = {mu_gma0, mu_eta0, mu_ksi1, mu_eta1, mu_ksi0, mu_gma1};
    double *diff_mu[6] = {diff_mu_gma0, diff_mu_eta0, diff_mu_ksi1, diff_mu_eta1, diff_mu_ksi0, diff_mu_gma1};

    // face tangent, coordinates 
    double coords[6][2] = {{ksi, eta}, 
                            {ksi, gma}, 
                            {eta, gma}, 
                            {ksi, gma}, 
                            {eta, gma}, 
                            {ksi, eta}};

    // polynomial orders for each face

    int pp[6][2] = {{p[0], p[1]},
                    {p[0], p[2]},
                    {p[1], p[2]},
                    {p[0], p[2]},
                    {p[1], p[2]},
                    {p[0], p[1]}};

    int pos[6][2] = {{0, 1}, {0, 2}, {1, 2}, {0,2}, {1, 2}, {0, 1}};

    for (int face = 0; face != 6; face++) {

      int *nodes = face_nodes[face];
      double s[2] = {0.0, 0.0}; 
      double orient_mat[2][2];
      
      CHKERR Face_orientMat(nodes, orient_mat);
      // cout << "oMat : [" << orient_mat[0][0] << ", " << orient_mat[0][1] << endl;
      // cout << "        " << orient_mat[1][0] << ", " << orient_mat[1][1] << "]" << endl;
      for (int r = 0; r != 2; r++){
        for (int c = 0; c != 2; c++){
          s[r] += orient_mat[r][c] * coords[face][c];
        }
      }
      // cout << "coords[face] : [" << coords[face][0] << ", " << coords[face][1] << "]" << endl; 
      // cout << "           s : [" << s[0] << ", " << s[1] << "]" << endl; 
      double diff_s[2][3] = {{0., 0., 0.}, {0., 0., 0.}};

      diff_s[0][pos[face][0]] = orient_mat[0][0];
      diff_s[0][pos[face][1]] = orient_mat[0][1];

      diff_s[1][pos[face][0]] = orient_mat[1][0];
      diff_s[1][pos[face][1]] = orient_mat[1][1];

      int pq[2] = {pp[face][0], pp[face][1]};


      double L0[pq[0] - 1];          double diffL0[pq[0] - 1];
      double L1[pq[1] - 1];          double diffL1[pq[1] - 1];

      CHKERR Integrated_Legendre(pq[0], s[0], L0, diffL0);
      CHKERR Integrated_Legendre(pq[1], s[1], L1, diffL1);

      int qd_shift = (pq[0] - 1) * (pq[1] - 1) * q;
      int n = 0;
      // cout << "=============================================" << endl;
      for (int s1 = 0; s1 != pq[0] - 1; s1++) {
        for (int s2 = 0; s2 != pq[1] - 1; s2++) {
          faceN[face][qd_shift + n] = mu[face] * L0[s1] * L1[s2];
          // cout << "mu[ff] : " << mu[face] << endl;
          // cout << "L0[s1] : " << L0[s1] << endl;
          // cout << "L1[s2] : " << L1[s2] << endl;

          diff_faceN[face][3 * (qd_shift + n) + 0] = diff_mu[face][0] * L0[s1] * L1[s2] +
                                                mu[face] * diffL0[s1] * diff_s[0][0] * L1[s2] +
                                                mu[face] * L0[s1] * diffL1[s1] * diff_s[1][0]; 

          diff_faceN[face][3 * (qd_shift + n) + 1] = diff_mu[face][1] * L0[s1] * L1[s2] +
                                                mu[face] * diffL0[s1] * diff_s[0][1] * L1[s2] +
                                                mu[face] * L0[s1] * diffL1[s1] * diff_s[1][1];

          diff_faceN[face][3 * (qd_shift + n) + 2] = diff_mu[face][2] * L0[s1] * L1[s2] +
                                                mu[face] * diffL0[s1] * diff_s[0][2] * L1[s2] +
                                                mu[face] * L0[s1] * diffL1[s1] * diff_s[1][2];

          ++n;
        }
      }
    }
  }
  MoFEMFunctionReturnHot(0);
}

MoFEMErrorCode MoFEM::H1_InteriorShapeFunctions_ONHEX(int *p, 
                                                      double *N,
                                                      double *faceN,
                                                      double *diff_faceN,
                                                      int nb_integration_pts) {
  MoFEMFunctionBeginHot;
  double vertices[8][3] = {{-1., -1., -1.}, {1., -1., -1.}, {1., 1., -1.},
                           {-1., 1., -1.},  {-1., -1., 1.}, {1., -1., 1.},
                           {1., 1., 1.},    {-1., 1., 1.}};
  for (int q = 0; q != nb_integration_pts; q++) {
    int shift = 8 * q;
    double ksi = 0.;
    double eta = 0.;
    double gma = 0.;
    
    for (int vv = 0; vv != 8; vv++) {
      ksi += N[shift + vv] * vertices[vv][0];
      eta += N[shift + vv] * vertices[vv][1];
      gma += N[shift + vv] * vertices[vv][2];
    }

    double diff_ksi[3] = {1.0, 0.0, 0.0};
    double diff_eta[3] = {0.0, 1.0, 0.0};
    double diff_gma[3] = {0.0, 0.0, 1.0};

    // Affine coordinates of each eadge mu0 and mu1
    double mu_ksi0 = 1.0 - 0.5 * (ksi + 1.0); double diff_mu_ksi0[3] = {-0.5, 0.0, 0.0};
    double mu_eta0 = 1.0 - 0.5 * (eta + 1.0); double diff_mu_eta0[3] = {0.0, -0.5, 0.0};
    double mu_gma0 = 1.0 - 0.5 * (gma + 1.0); double diff_mu_gma0[3] = {0.0, 0.0, -0.5};

    double mu_ksi1 = 0.5 * (ksi + 1.0);  double diff_mu_ksi1[3] = {0.5, 0.0, 0.0};
    double mu_eta1 = 0.5 * (eta + 1.0);  double diff_mu_eta1[3] = {0.0, 0.5, 0.0};
    double mu_gma1 = 0.5 * (gma + 1.0);  double diff_mu_gma1[3] = {0.0, 0.0, 0.5};

    double s[3] = {ksi, eta, gma};
    double diff_s[3][3] = {{1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {0.0, 0.0, 1.0}};
    double L0[p[0] - 1];         double diffL0[p[0] - 1];
    double L1[p[1] - 1];         double diffL1[p[1] - 1];
    double L2[p[2] - 1];         double diffL2[p[2] - 1];

    CHKERR Integrated_Legendre(p[0], ksi, L0, diffL0);
    CHKERR Integrated_Legendre(p[1], eta, L1, diffL1);
    CHKERR Integrated_Legendre(p[2], gma, L2, diffL2);

    // cout << "In Face H1" << endl;

    int qd_shift = (p[0] - 1) * (p[1] - 1) * (p[2] - 1) * q;
    int n = 0;
    for (int s1 = 0; s1 != p[0] - 1; s1++) {
      for (int s2 = 0; s2 != p[1] - 1; s2++) {
        for (int s3 = 0; s3 < p[2] - 1; s3++)
        {
          faceN[qd_shift + n] = L0[s1] * L1[s2] * L2[s3];

          diff_faceN[3 * (qd_shift + n) + 0] = diffL0[s1] * diff_s[0][0] * L1[s2] * L2[s3] +
                                               L0[s1] * diffL1[s2] * diff_s[1][0] * L2[s3] +
                                               L0[s1] * L1[s2] * diffL2[s3] * diff_s[2][0];

          diff_faceN[3 * (qd_shift + n) + 1] = diffL0[s1] * diff_s[0][1] * L1[s2] * L2[s3] +
                                               L0[s1] * diffL1[s2] * diff_s[1][1] * L2[s3] +
                                               L0[s1] * L1[s2] * diffL2[s3] * diff_s[2][1];

          diff_faceN[3 * (qd_shift + n) + 2] = diffL0[s1] * diff_s[0][2] * L1[s2] * L2[s3] +
                                               L0[s1] * diffL1[s2] * diff_s[1][2] * L2[s3] +
                                               L0[s1] * L1[s2] * diffL2[s3] * diff_s[2][2];

          ++n;
        } 
      }
    }
  }
  MoFEMFunctionReturnHot(0);
}
                   