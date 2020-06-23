/** \file EdgeQuadHexPolynomials.cpp

  \brief Implementation of hierarchical Edge, Quad, and Hex shape functions of 
  type H1, Hcurl, Hdiv
*/

using namespace MoFEM;

// Auxilary functions

MoFEMErrorCode MoFEM::Legendre_polynomials01(int p, double s01, double *L) {

  MoFEMFunctionBeginHot;
  double s = 2.0 * s01 - 1.0;
  CHKERR Legendre_polynomials(p, s, NULL, L, NULL, 1);

  MoFEMFunctionReturnHot(0);
}

MoFEMErrorCode MoFEM::Integrated_Legendre01(int p, double s01, double *L,
                                          double *diffL) {

  MoFEMFunctionBeginHot;

  double l[p + 1];
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
MoFEMErrorCode MoFEM::H1_BubbleShapeFunctions_ONSEGMENT(int p, double *N,
                                                      double *bubbleN,
                                                      double *diff_bubbleN,
                                                      int nb_integration_pts) {
  MoFEMFunctionBeginHot;
  double coords[2] = {0., 1.}; // vertex coordinates of the reference element [0, 1]
  for (int q = 0; q != nb_integration_pts; q++) {
    int shift = 2 * q;
    
    double ksi = coords[0] * N[shift + 0] + coords[1] * N[shift + 1];

    double mu[2] = {1.0 - ksi, ksi};   // affine coordinates of the ref. element [0, 1]

    double diff_mu[2] = {-1.0, 1.0};  

    
    double L[p-1];
    double diffL[p-1];
    CHKERR Integrated_Legendre01(p, mu[1], L, diffL); // diffL = (d L)/(d ksi)
    int qd_shift = (p - 1) * q;
    for (int n = 0; n != p-1; n++){
      bubbleN[qd_shift + n] = L[n];
      diff_bubbleN[qd_shift + n] = diffL[n] * diff_mu[1]; // diff_bubbleN = (d bubbleN) / (d gp) 
    }
    
  }   
  MoFEMFunctionReturnHot(0);
}

MoFEMErrorCode MoFEM::L2_ShapeFunctions_ONSEGMENT(int p, double *N,
                                                  double *funN,
                                                  int nb_integration_pts) {
  MoFEMFunctionBeginHot;
  double coords[2] = {0., 1.}; // vertex coordinates of the reference element [0, 1]
  for (int q = 0; q != nb_integration_pts; q++) {
    int shift = 2 * q;
    double ksi = coords[0] * N[shift + 0] + coords[1] * N[shift + 1];

    double mu[2] = {1.0 - ksi, ksi};   // affine coordinates over the ref element [0, 1]
    
    double L[p];
    CHKERR Legendre_polynomials01(p - 1, ksi, L);
    int qd_shift = p * q;
    for (int n = 0; n != p; n++) {
      funN[qd_shift + n] = L[n] * mu[1];
    }
  }   
  MoFEMFunctionReturnHot(0);
}

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


MoFEMErrorCode MoFEM::H1_EdgeShapeFunctions_ONQUAD(int *sense, int *p,
                                                   double *N, double
                                                   *edgeN[4], double
                                                   *diff_edgeN[4], int
                                                   nb_integration_pts) {
  MoFEMFunctionBeginHot;
  double coords[4][2] = {{0.0, 0.0},
                      {1.0, 0.0},
                      {1.0, 1.0},
                      {0.0, 1.0}};
    for (int q = 0; q != nb_integration_pts; q++) {
      int shift = 4 * q;
      double ksi = 0.0; double eta = 0.0;
      for (int vv = 0; vv != 4; vv++)
      {
        ksi += coords[vv][0] * N[shift + vv];
        eta += coords[vv][1] * N[shift + vv]; 
      }
      

      // Affine coordinates of each eadge mu0 and mu1
      double mu_ksi[2] = {1.0 - ksi, ksi};
      double mu_eta[2] = {1.0 - eta, eta};

      double diff_mu_ksi[2][2] = {{-1.0, 0.0}, {1.0, 0.0}};
      double diff_mu_eta[2][2] = {{0.0, -1.0}, {0.0, 1.0}};

      // constant affine coordinates of edges per cannonical numbering
      double mu_const[4] = {mu_eta[0], mu_ksi[1], mu_eta[1], mu_ksi[0]};
      double *diff_mu_const[4] = {diff_mu_eta[0], diff_mu_ksi[1], diff_mu_eta[1], diff_mu_ksi[0]};

      // parametrization of each edge per cannonical numbering
      double *mu[4] = {mu_ksi, mu_eta, mu_ksi, mu_eta};
      double *diff_mu[4][2] = {{&*diff_mu_ksi[0], &*diff_mu_ksi[1]}, 
                               {&*diff_mu_eta[0], &*diff_mu_eta[1]}, 
                               {&*diff_mu_ksi[0], &*diff_mu_ksi[1]}, 
                               {&*diff_mu_eta[0], &*diff_mu_eta[1]}};

      for (int e = 0; e != 4; e++) {
        double L[p[e] - 1];
        double diffL[p[e] - 1];
        int index = (sense[e] + 1) / 2;
        double ss = mu[e][index];
        CHKERR Integrated_Legendre01(p[e], ss, L, diffL);
        int qd_shift = (p[e]- 1) * q;
        for (int n = 0; n != p[e]-1; n++) {
          edgeN[e][qd_shift + n] = mu_const[e] * L[n];
          diff_edgeN[e][2 * qd_shift + 2 * n] =
              mu_const[e] * diffL[n] * diff_mu[e][index][0] + diff_mu_const[e][0] *
              L[n];
          diff_edgeN[e][2 * qd_shift + 2 * n + 1] =
              mu_const[e] * diffL[n] * diff_mu[e][index][1] + diff_mu_const[e][1] *
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
  double coords[4][2] = {{0.0, 0.0},
                      {1.0, 0.0},
                      {1.0, 1.0},
                      {0.0, 1.0}};
    for (int q = 0; q != nb_integration_pts; q++) {
      int shift = 4 * q;
      double ksi = 0.0; double eta = 0.0;
      for (int vv = 0; vv != 4; vv++)
      {
        ksi += coords[vv][0] * N[shift + vv];
        eta += coords[vv][1] * N[shift + vv]; 
      }
      

      // Affine coordinates of each eadge mu0 and mu1
      double mu_ksi[2] = {1.0 - ksi, ksi};
      double mu_eta[2] = {1.0 - eta, eta};

      double diff_mu_ksi[2][2] = {{-1.0, 0.0}, {1.0, 0.0}};
      double diff_mu_eta[2][2] = {{0.0, -1.0}, {0.0, 1.0}};
    
    double mu[2] = {mu_ksi[1], mu_eta[1]};
    double *diff_mu[2] = {diff_mu_ksi[1], diff_mu_eta[1]};
  
    double L0[p[0]-1];   double diffL0[p[0]-1];
    double L1[p[1]-1];   double diffL1[p[1]-1];
    CHKERR Integrated_Legendre01(p[0], mu[0], L0, diffL0);
    CHKERR Integrated_Legendre01(p[1], mu[1], L1, diffL1);

    // cout << "In Face H1" << endl;

    int qd_shift = (p[0] - 1) * (p[1] - 1) * q;
    int n = 0;
    for (int s1 = 0; s1 != p[0]-1; s1++){
      for (int s2 = 0; s2 != p[1]-1; s2++){
        faceN[qd_shift + n] = L0[s1] * L1[s2];
        diff_faceN[2 * (qd_shift + n) + 0] = diffL0[s1] * L1[s2] * diff_mu[0][0] + L0[s1] * diffL1[s2] * diff_mu[1][0];
        diff_faceN[2 * (qd_shift + n) + 1] = diffL0[s1] * L1[s2] * diff_mu[0][1] + L0[s1] * diffL1[s2] * diff_mu[1][1];
        ++n;
      }
    }
  }
  MoFEMFunctionReturnHot(0);
}

MoFEMErrorCode MoFEM::L2_FaceShapeFunctions_ONQUAD(int *p, double *N, double *faceN,
                                                   double *diff_faceN, int nb_integration_pts) {
  MoFEMFunctionBeginHot;
  double coords[4][2] = {{0.0, 0.0},
                      {1.0, 0.0},
                      {1.0, 1.0},
                      {0.0, 1.0}};
    for (int q = 0; q != nb_integration_pts; q++) {
      int shift = 4 * q;
      double ksi = 0.0; double eta = 0.0;
      for (int vv = 0; vv != 4; vv++)
      {
        ksi += coords[vv][0] * N[shift + vv];
        eta += coords[vv][1] * N[shift + vv]; 
      }
      

      // Affine coordinates of each eadge mu0 and mu1
      double mu_ksi[2] = {1.0 - ksi, ksi};
      double mu_eta[2] = {1.0 - eta, eta};

      double diff_mu_ksi[2][2] = {{-1.0, 0.0}, {1.0, 0.0}};
      double diff_mu_eta[2][2] = {{0.0, -1.0}, {0.0, 1.0}};

    double L0[p[0]];
    double L1[p[1]];

    CHKERR Legendre_polynomials01(p[0] - 1, mu_ksi[1], L0);
    CHKERR Legendre_polynomials01(p[1] - 1, mu_eta[1], L1);

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
  double coords[4][2] = {{0.0, 0.0},
                         {1.0, 0.0},
                         {1.0, 1.0},
                         {0.0, 1.0}};
  for (int q = 0; q != nb_integration_pts; q++) {
    int shift = 4 * q;
    double ksi = 0.0; double eta = 0.0;
    for (int vv = 0; vv != 4; vv++)
    {
      ksi += coords[vv][0] * N[shift + vv];
      eta += coords[vv][1] * N[shift + vv]; 
    }
    

    // Affine coordinates of each eadge mu0 and mu1
    double mu_ksi[2] = {1.0 - ksi, ksi};
    double mu_eta[2] = {1.0 - eta, eta};

    double diff_mu_ksi[2][2] = {{-1.0, 0.0}, {1.0, 0.0}};
    double diff_mu_eta[2][2] = {{0.0, -1.0}, {0.0, 1.0}};

    // constant affine coordinates of edges per cannonical numbering
      double mu_const[4] = {mu_eta[0], mu_ksi[1], mu_eta[1], mu_ksi[0]};
      double *diff_mu_const[4] = {diff_mu_eta[0], diff_mu_ksi[1], diff_mu_eta[1], diff_mu_ksi[0]};

    // parametrization of each edge per cannonical numbering
    double *mu[4] = {mu_ksi, mu_eta, mu_ksi, mu_eta};
    double *diff_mu[4][2] = {{&*diff_mu_ksi[0], &*diff_mu_ksi[1]},
                             {&*diff_mu_eta[0], &*diff_mu_eta[1]},
                             {&*diff_mu_ksi[0], &*diff_mu_ksi[1]},
                             {&*diff_mu_eta[0], &*diff_mu_eta[1]}};

    int pp[4] = {p[0], p[1], p[0], p[1]};
    for (int e = 0; e != 4; e++) {
      int index = (sense[e] + 1) / 2;
      double L[pp[e]];
      double ss = mu[e][index];
      CHKERR Legendre_polynomials01(pp[e] - 1, ss, L);
      int qd_shift = pp[e] * q;
      for (int n = 0; n != pp[e]; n++) {
        edgeN[e][2 * (qd_shift + n) + 0] = mu_const[e] * L[n] * diff_mu[e][index][0];
        edgeN[e][2 * (qd_shift + n) + 1] = mu_const[e] * L[n] * diff_mu[e][index][1];

        double E1[2] = {diff_mu_const[e][0], diff_mu_const[e][1]};
        double E2[2] = {L[n] * diff_mu[e][index][0], L[n] * diff_mu[e][index][1]};
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
  double coords[4][2] = {{0.0, 0.0},
                      {1.0, 0.0},
                      {1.0, 1.0},
                      {0.0, 1.0}};
  for (int q = 0; q != nb_integration_pts; q++) {
    int shift = 4 * q;
    double ksi = 0.0; double eta = 0.0;
    for (int vv = 0; vv != 4; vv++)
    {
      ksi += coords[vv][0] * N[shift + vv];
      eta += coords[vv][1] * N[shift + vv]; 
    }
    

    // Affine coordinates of each eadge mu0 and mu1
    double mu_ksi[2] = {1.0 - ksi, ksi};
    double mu_eta[2] = {1.0 - eta, eta};

    

    double diff_mu_ksi[2][2] = {{-1.0, 0.0}, {1.0, 0.0}};
    double diff_mu_eta[2][2] = {{0.0, -1.0}, {0.0, 1.0}};

    int pq[2] = {p[0], p[1]};
    int qp[2] = {p[1], p[0]};
    double mu_ksi_eta[2] = {mu_ksi[1], mu_eta[1]};
    double mu_eta_ksi[2] = {mu_eta[1], mu_ksi[1]};

    double *diff_mu_ksi_eta[2] = {diff_mu_ksi[1], diff_mu_eta[1]};
    double *diff_mu_eta_ksi[2] = {diff_mu_eta[1], diff_mu_ksi[1]};
    double sgn[2] = {-1.0, 1.0};
    for (int typ = 0; typ != 2; typ++)
    {
      int pp = pq[typ];
      

      double Phi[pp - 1];
      double diffPhi[pp - 1];
      CHKERR Integrated_Legendre01(pp, mu_ksi_eta[typ], Phi, diffPhi);

      
      int qq = qp[typ];
      double E[qq];
      CHKERR Legendre_polynomials01(qq - 1, mu_eta_ksi[typ], E);

      int qd_shift = (pp - 1) * qq * q;
      int n = 0;
      for (int i = 0; i != qq; i++) {
        for (int j = 0; j != pp - 1; j++) {
          faceN[typ][2 * (qd_shift + n) + 0] = Phi[j] * E[i] * diff_mu_eta_ksi[typ][0];
          faceN[typ][2 * (qd_shift + n) + 1] = Phi[j] * E[i] * diff_mu_eta_ksi[typ][1];
          



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
  double coords[4][2] = {{0.0, 0.0},
                      {1.0, 0.0},
                      {1.0, 1.0},
                      {0.0, 1.0}};
  for (int q = 0; q != nb_integration_pts; q++) {
    int shift = 4 * q;
    double ksi = 0.0; double eta = 0.0;
    for (int vv = 0; vv != 4; vv++)
    {
      ksi += coords[vv][0] * N[shift + vv];
      eta += coords[vv][1] * N[shift + vv]; 
    }
    

    // Affine coordinates of each eadge mu0 and mu1
    double mu_ksi[2] = {1.0 - ksi, ksi};
    double mu_eta[2] = {1.0 - eta, eta};

    double diff_mu_ksi[2][2] = {{-1.0, 0.0}, {1.0, 0.0}};
    double diff_mu_eta[2][2] = {{0.0, -1.0}, {0.0, 1.0}};

    // constant affine coordinates of edges per cannonical numbering
      double mu_const[4] = {mu_eta[0], mu_ksi[1], mu_eta[1], mu_ksi[0]};
      double *diff_mu_const[4] = {diff_mu_eta[0], diff_mu_ksi[1], diff_mu_eta[1], diff_mu_ksi[0]};

    // parametrization of each edge per cannonical numbering
    double *mu[4] = {mu_ksi, mu_eta, mu_ksi, mu_eta};
    double *diff_mu[4][2] = {{&*diff_mu_ksi[0], &*diff_mu_ksi[1]},
                             {&*diff_mu_eta[0], &*diff_mu_eta[1]},
                             {&*diff_mu_ksi[0], &*diff_mu_ksi[1]},
                             {&*diff_mu_eta[0], &*diff_mu_eta[1]}};

    int pp[4] = {p[0], p[1], p[0], p[1]};
    for (int e = 0; e != 4; e++) {
      int index = (sense[e] + 1) / 2;
      double L[pp[e]];
      double ss = mu[e][index];
      CHKERR Legendre_polynomials01(pp[e] - 1, ss, L);
      int qd_shift = pp[e] * q;
      for (int n = 0; n != pp[e]; n++) {
        edgeN[e][2 * (qd_shift + n) + 0] =  mu_const[e] * L[n] * diff_mu[e][index][1];
        edgeN[e][2 * (qd_shift + n) + 1] = -mu_const[e] * L[n] * diff_mu[e][index][0];

        double E1[2] = {diff_mu_const[e][0], diff_mu_const[e][1]};
        double E2[2] = {L[n] * diff_mu[e][index][0], L[n] * diff_mu[e][index][1]};
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
  double coords[4][2] = {{0.0, 0.0},
                         {1.0, 0.0},
                         {1.0, 1.0},
                         {0.0, 1.0}};
  for (int q = 0; q != nb_integration_pts; q++) {
    int shift = 4 * q;
    double ksi = 0.0; double eta = 0.0;
    for (int vv = 0; vv != 4; vv++)
    {
      ksi += coords[vv][0] * N[shift + vv];
      eta += coords[vv][1] * N[shift + vv]; 
    }
    

    // Affine coordinates of each eadge mu0 and mu1
    double mu_ksi[2] = {1.0 - ksi, ksi};
    double mu_eta[2] = {1.0 - eta, eta};

    

    double diff_mu_ksi[2][2] = {{-1.0, 0.0}, {1.0, 0.0}};
    double diff_mu_eta[2][2] = {{0.0, -1.0}, {0.0, 1.0}};

    int pq[2] = {p[0], p[1]};
    int qp[2] = {p[1], p[0]};
    double mu_ksi_eta[2] = {mu_ksi[1], mu_eta[1]};
    double mu_eta_ksi[2] = {mu_eta[1], mu_ksi[1]};

    double *diff_mu_ksi_eta[2] = {diff_mu_ksi[1], diff_mu_eta[1]};
    double *diff_mu_eta_ksi[2] = {diff_mu_eta[1], diff_mu_ksi[1]};
    double sgn[2] = {-1.0, 1.0};
    for (int typ = 0; typ != 2; typ++){
      int pp = pq[typ];
      

      double Phi[pp - 1];
      double diffPhi[pp - 1];
      CHKERR Integrated_Legendre01(pp, mu_ksi_eta[typ], Phi, diffPhi);

      
      int qq = qp[typ];
      double E[qq];
      CHKERR Legendre_polynomials01(qq - 1, mu_eta_ksi[typ], E);

      int qd_shift = (pp - 1) * qq * q;
      int n = 0;
      for (int i = 0; i != qq; i++) {
        for (int j = 0; j != pp - 1; j++) {
          faceN[typ][2 * (qd_shift + n) + 0] = Phi[j] * E[i] * diff_mu_eta_ksi[typ][1];
          faceN[typ][2 * (qd_shift + n) + 1] = -Phi[j] * E[i] * diff_mu_eta_ksi[typ][0];

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
  RefHex ref_hex(N, nb_integration_pts);
  auto edge_coords     = ref_hex.get_edge_coords();
  auto edge_diff_coords = ref_hex.get_edge_diff_coords();
  auto edge_affine      = ref_hex.get_edge_affines();
  auto edge_diff_affine  = ref_hex.get_edge_diff_affines();

  for (int q = 0; q != nb_integration_pts; q++) {

    int pp[12] = {p[0], p[1], p[0], p[1], p[2], p[2], p[2], p[2], p[0], p[1], p[0], p[1]};

    for (int e = 0; e != 12; e++) {
      double L[pp[e] - 1];
      double diffL[pp[e] - 1];
      double ss = (double)sense[e] * edge_coords[q][e] + 0.5 * (1.0 - (double)sense[e]);
      
      CHKERR Integrated_Legendre01(pp[e], ss, L, diffL);
      int qd_shift = (pp[e] - 1) * q;

      for (int n = 0; n != pp[e] - 1; n++) {
        edgeN[e][qd_shift + n] = edge_affine[q][2 * e + 0] * edge_affine[q][2 * e + 1] * L[n];
        
        diff_edgeN[e][3 * (qd_shift + n) + 0] = edge_affine[q][2 * e + 0] * edge_affine[q][2 * e + 1] * 
                                                diffL[n] * sense[e] * edge_diff_coords[e][0] + 
                                                (edge_affine[q][2 * e + 0] * edge_diff_affine[e][0][0] + 
                                                 edge_affine[q][2 * e + 1] * edge_diff_affine[e][1][0]) * L[n];

        diff_edgeN[e][3 * (qd_shift + n) + 1] = edge_affine[q][2 * e + 0] * edge_affine[q][2 * e + 1] * 
                                                diffL[n] * sense[e] * edge_diff_coords[e][1] + 
                                                (edge_affine[q][2 * e + 0] * edge_diff_affine[e][0][1] + 
                                                 edge_affine[q][2 * e + 1] * edge_diff_affine[e][1][1]) * L[n];

        diff_edgeN[e][3 * (qd_shift + n) + 2] = edge_affine[q][2 * e + 0] * edge_affine[q][2 * e + 1] * 
                                                diffL[n] * sense[e] * edge_diff_coords[e][2] + 
                                                (edge_affine[q][2 * e + 0] * edge_diff_affine[e][0][2] + 
                                                 edge_affine[q][2 * e + 1] * edge_diff_affine[e][1][2]) * L[n];

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

  RefHex ref_hex(N, nb_integration_pts);
  auto face_coords = ref_hex.get_face_coords(face_nodes);
  auto face_diff_coords = ref_hex.get_face_diff_coords(face_nodes);
  auto face_affine = ref_hex.get_face_affines();
  auto face_diff_affine = ref_hex.get_face_diff_affines();
  for (int q = 0; q != nb_integration_pts; q++) {
        int pp[6][2] = {{p[0], p[1]},
                        {p[0], p[2]},
                        {p[1], p[2]},
                        {p[0], p[2]},
                        {p[1], p[2]},
                        {p[0], p[1]}};

    for (int face = 0; face != 6; face++) {
      int pq[2] = {pp[face][0], pp[face][1]};


      double L0[pq[0] - 1];          double diffL0[pq[0] - 1];
      double L1[pq[1] - 1];          double diffL1[pq[1] - 1];

      double ss0 = face_coords[q][2 * face + 0];
      double ss1 = face_coords[q][2 * face + 1];

      CHKERR Integrated_Legendre01(pq[0], ss0, L0, diffL0);
      CHKERR Integrated_Legendre01(pq[1], ss1, L1, diffL1);

      int qd_shift = (pq[0] - 1) * (pq[1] - 1) * q;
      int n = 0;
      for (int s1 = 0; s1 != pq[0] - 1; s1++) {
        for (int s2 = 0; s2 != pq[1] - 1; s2++) {
          faceN[face][qd_shift + n] = face_affine[q][face] * L0[s1] * L1[s2];

          diff_faceN[face][3 * (qd_shift + n) + 0] = face_diff_affine[face][0] * L0[s1] * L1[s2] +
                                                     face_affine[q][face] * diffL0[s1] * face_diff_coords[face][0][0] * L1[s2] +
                                                     face_affine[q][face] * L0[s1] * diffL1[s1] * face_diff_coords[face][1][0]; 

                                                 

          diff_faceN[face][3 * (qd_shift + n) + 1] = face_diff_affine[face][1] * L0[s1] * L1[s2] +
                                                     face_affine[q][face] * diffL0[s1] * face_diff_coords[face][0][1] * L1[s2] +
                                                     face_affine[q][face] * L0[s1] * diffL1[s1] * face_diff_coords[face][1][1]; 

          diff_faceN[face][3 * (qd_shift + n) + 2] = face_diff_affine[face][2] * L0[s1] * L1[s2] +
                                                     face_affine[q][face] * diffL0[s1] * face_diff_coords[face][0][2] * L1[s2] +
                                                     face_affine[q][face] * L0[s1] * diffL1[s1] * face_diff_coords[face][1][2]; 

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
  RefHex ref_hex(N, nb_integration_pts);
  auto volume_coords = ref_hex.get_integrationPts();
  auto volume_diff_coords = ref_hex.get_volume_diff_coords();
  for (int q = 0; q != nb_integration_pts; q++) {


  

    double ksi = volume_coords[q][0];
    double eta = volume_coords[q][1];
    double gma = volume_coords[q][2];

    double L0[p[0] - 1];         double diffL0[p[0] - 1];
    double L1[p[1] - 1];         double diffL1[p[1] - 1];
    double L2[p[2] - 1];         double diffL2[p[2] - 1];

    CHKERR Integrated_Legendre01(p[0], ksi, L0, diffL0);
    CHKERR Integrated_Legendre01(p[1], eta, L1, diffL1);
    CHKERR Integrated_Legendre01(p[2], gma, L2, diffL2);

    // cout << "In Face H1" << endl;

    int qd_shift = (p[0] - 1) * (p[1] - 1) * (p[2] - 1) * q;
    int n = 0;
    for (int s1 = 0; s1 != p[0] - 1; s1++) {
      for (int s2 = 0; s2 != p[1] - 1; s2++) {
        for (int s3 = 0; s3 < p[2] - 1; s3++)
        {
          faceN[qd_shift + n] = L0[s1] * L1[s2] * L2[s3];
                                                            
          diff_faceN[3 * (qd_shift + n) + 0] = diffL0[s1] * volume_diff_coords[0][0] * L1[s2] * L2[s3] +
                                               L0[s1] * diffL1[s2] * volume_diff_coords[1][0] * L2[s3] +
                                               L0[s1] * L1[s2] * diffL2[s3] * volume_diff_coords[2][0];

          diff_faceN[3 * (qd_shift + n) + 1] = diffL0[s1] * volume_diff_coords[0][1] * L1[s2] * L2[s3] +
                                               L0[s1] * diffL1[s2] * volume_diff_coords[1][1] * L2[s3] +
                                               L0[s1] * L1[s2] * diffL2[s3] * volume_diff_coords[2][1];

          diff_faceN[3 * (qd_shift + n) + 2] = diffL0[s1] * volume_diff_coords[0][2] * L1[s2] * L2[s3] +
                                               L0[s1] * diffL1[s2] * volume_diff_coords[1][2] * L2[s3] +
                                               L0[s1] * L1[s2] * diffL2[s3] * volume_diff_coords[2][2];

          ++n;
        } 
      }
    }
  }
  MoFEMFunctionReturnHot(0);
}

MoFEMErrorCode MoFEM::L2_InteriorShapeFunctions_ONHEX(int *p, double *N,
                                                      double *volN,
                                                      int nb_integration_pts){
  MoFEMFunctionBeginHot;
  RefHex ref_hex(N, nb_integration_pts);
  auto volume_coords = ref_hex.get_integrationPts();

  

  for(int qq = 0; qq != nb_integration_pts; qq++){
    double ksi = volume_coords[qq][0];
    double eta = volume_coords[qq][1];
    double gma = volume_coords[qq][2];

    double P0[p[0]];
    double P1[p[1]];
    double P2[p[2]];
   
    int qd_shift = qq * p[0] * p[1] * p[2];

    CHKERR Legendre_polynomials01(p[0]-1, ksi, P0);
    CHKERR Legendre_polynomials01(p[1]-1, eta, P1);
    CHKERR Legendre_polynomials01(p[2]-1, gma, P2);
    int n = 0;
    for(int ii = 0; ii != p[0]; ii++){
      for(int jj = 0; jj != p[1]; jj++){
        for(int kk = 0; kk != p[2]; kk++){
          volN[qd_shift + n] = P0[ii] * P1[jj] * P2[kk];
        }
      }
    }
  }

  MoFEMFunctionReturnHot(0);
}

MoFEMErrorCode MoFEM::Hcurl_EdgeShapeFunctions_ONHEX(int *sense, int *p,
                                                     double *N,
                                                     double *edgeN[12],
                                                     double *curl_edgeN[12],
                                                     int nb_integration_pts) {
  MoFEMFunctionBeginHot;
  RefHex ref_hex(N, nb_integration_pts);
  auto mu = ref_hex.get_edge_affines();
  auto diff_mu = ref_hex.get_edge_diff_affines();
  auto ksi = ref_hex.get_edge_coords();
  auto diff_ksi = ref_hex.get_edge_diff_coords();

  for (int qq = 0; qq != nb_integration_pts; qq++){
    int pp[12] = {p[0], p[1], p[0], p[1], p[2], p[2], p[2], p[2] ,p[0], p[1], p[0], p[1]};
    for (int ee = 0; ee != 12; ee++){

      double L[pp[ee]];
      double diff_ss[3] = {(double)sense[ee] * diff_ksi[ee][0],
                           (double)sense[ee] * diff_ksi[ee][1],
                           (double)sense[ee] * diff_ksi[ee][2]};
      double ss = (double)sense[ee] * ksi[qq][ee] + 0.5 * (1.0 - (double)sense[ee]);
      CHKERR Legendre_polynomials01(pp[ee] - 1, ss, L);
      int qd_shift = pp[ee] * qq;

      for (int ii = 0; ii != pp[ee]; ii++){
        edgeN[ee][3 * (qd_shift + ii) + 0] = mu[qq][2 * ee + 0] * mu[qq][2 * ee + 1] * L[ii] * diff_ss[0];
        edgeN[ee][3 * (qd_shift + ii) + 1] = mu[qq][2 * ee + 0] * mu[qq][2 * ee + 1] * L[ii] * diff_ss[1];
        edgeN[ee][3 * (qd_shift + ii) + 2] = mu[qq][2 * ee + 0] * mu[qq][2 * ee + 1] * L[ii] * diff_ss[2];

        double E1[3] = {diff_mu[ee][0][0] * mu[qq][2 * ee + 1] + mu[qq][2 * ee + 0] * diff_mu[ee][1][0],
                        diff_mu[ee][0][1] * mu[qq][2 * ee + 1] + mu[qq][2 * ee + 0] * diff_mu[ee][1][1],
                        diff_mu[ee][0][2] * mu[qq][2 * ee + 1] + mu[qq][2 * ee + 0] * diff_mu[ee][1][2]};
        double E2[3] = {L[ii] * diff_ss[0], L[ii] * diff_ss[1], L[ii] * diff_ss[2]};

        auto cross_product = ref_hex.Cross_product(E1, E2);
        curl_edgeN[ee][3 * (qd_shift + ii) + 0] = cross_product[0];
        curl_edgeN[ee][3 * (qd_shift + ii) + 1] = cross_product[1];
        curl_edgeN[ee][3 * (qd_shift + ii) + 2] = cross_product[2];
      }
      
    }
    
  }
  

  MoFEMFunctionReturnHot(0);
}

MoFEMErrorCode MoFEM::Hcurl_FaceShapeFunctions_ONHEX(int *face_nodes[6], 
                                                     int *p,
                                                     double *N,
                                                     double *faceN[6][2],
                                                     double *curl_faceN[6][2],
                                                     int nb_integration_pts) {
  MoFEMFunctionBeginHot;
  RefHex ref_hex(N, nb_integration_pts);
  auto ksi = ref_hex.get_face_coords(face_nodes);
  auto diff_ksi = ref_hex.get_face_diff_coords(face_nodes);
  auto mu = ref_hex.get_face_affines();
  auto diff_mu = ref_hex.get_face_diff_affines();
  for (int qq = 0; qq != nb_integration_pts; qq++){

    for(int ff = 0; ff != 6; ff++){
      double ksi_eta[2] = {ksi[qq][2 * ff + 0], ksi[qq][2 * ff + 1]};
      double eta_ksi[2] = {ksi[qq][2 * ff + 1], ksi[qq][2 * ff + 0]};

      double *diff_ksi_eta[2] = {diff_ksi[ff][0], diff_ksi[ff][1]};
      double *diff_eta_ksi[2] = {diff_ksi[ff][1], diff_ksi[ff][0]};

      int pq[2] = {p[0], p[1]};
      int qp[2] = {p[1], p[0]};
      
      for (int fam = 0; fam != 2; fam++){
        int pp = pq[fam];

        double Phi[pp - 1];
        double diffPhi[pp - 1];
        CHKERR Integrated_Legendre01(pp, ksi_eta[fam], Phi, diffPhi);

        int qq = qp[fam];
        double E[qq];
        CHKERR Legendre_polynomials01(qq - 1, eta_ksi[fam], E);

        int qd_shift = (pp - 1) * qq * qq;
        int n = 0;
        for (int i = 0; i != qq; i++) {
          for (int j = 0; j != pp - 1; j++) {
            faceN[ff][fam][3 * (qd_shift + n) + 0] = mu[qq][ff] * Phi[j] * E[i] * diff_eta_ksi[fam][0];
            faceN[ff][fam][3 * (qd_shift + n) + 0] = mu[qq][ff] * Phi[j] * E[i] * diff_eta_ksi[fam][1];
            faceN[ff][fam][3 * (qd_shift + n) + 0] = mu[qq][ff] * Phi[j] * E[i] * diff_eta_ksi[fam][2];

            double Ei[3] = {E[i] * diff_eta_ksi[fam][0],
                            E[i] * diff_eta_ksi[fam][1],
                            E[i] * diff_eta_ksi[fam][2]};
            double Eij[3] = {Phi[j] * Ei[0], Phi[j] * Ei[1], Phi[j] * Ei[2]};

            double diff_Phi[3] = {diffPhi[j] * diff_ksi_eta[fam][0],
                                  diffPhi[j] * diff_ksi_eta[fam][1],
                                  diffPhi[j] * diff_ksi_eta[fam][2]};
            double diff_muVec[3] = {diff_mu[ff][0], diff_mu[ff][1], diff_mu[ff][2]};
            auto diff_mu_cross_Eij = ref_hex.Cross_product(diff_muVec, Eij);
            auto diff_Phi_cross_Ei = ref_hex.Cross_product(diff_Phi, Ei);

            curl_faceN[ff][fam][3 * (qd_shift + n) + 0] = mu[qq][ff] * diff_Phi_cross_Ei[0] + diff_mu_cross_Eij[0];
            curl_faceN[ff][fam][3 * (qd_shift + n) + 0] = mu[qq][ff] * diff_Phi_cross_Ei[1] + diff_mu_cross_Eij[1];
            curl_faceN[ff][fam][3 * (qd_shift + n) + 0] = mu[qq][ff] * diff_Phi_cross_Ei[2] + diff_mu_cross_Eij[2];
            ++n;
          }
        }
      }
      

    }


  }
  MoFEMFunctionReturnHot(0);
}

MoFEMErrorCode MoFEM::Hcurl_InteriorShapeFunctions_ONHEX(int *p, 
                                                         double *N,
                                                         double *volN[3],
                                                         double *curl_volN[3],
                                                         int nb_integration_pts){
  MoFEMFunctionBeginHot;
  RefHex ref_hex(N, nb_integration_pts);
  auto ksi = ref_hex.get_integrationPts();
  auto diff_ksi = ref_hex.get_volume_diff_coords();

  for (int qq = 0; qq < nb_integration_pts; qq++){
    double ksi_eta_gma[3] = {ksi[qq][0], ksi[qq][1], ksi[qq][2]};
    double eta_gma_ksi[3] = {ksi[qq][1], ksi[qq][2], ksi[qq][0]};   
    double gma_ksi_eta[3] = {ksi[qq][2], ksi[qq][0], ksi[qq][1]};   

    double *diff_ksi_eta_gma[3] = {diff_ksi[0], diff_ksi[1], diff_ksi[2]};
    double *diff_eta_gma_ksi[3] = {diff_ksi[1], diff_ksi[2], diff_ksi[0]};
    double *diff_gma_ksi_eta[3] = {diff_ksi[2], diff_ksi[0], diff_ksi[1]};

    int pqr[3] = {p[0], p[1], p[2]};
    int qrp[3] = {p[1], p[2], p[0]};
    int rpq[3] = {p[2], p[0], p[1]};
    for (int fam = 0; fam < 3; fam++){

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

      int qd_shift = (ppp - 1) * qqq * (rrr - 1) * qq;
      int n = 0;
      for (int ii = 0; ii < qqq; ii++){
        for (int jj = 0; jj < ppp - 1; jj++){
          for (int kk = 0; kk < rrr - 1; kk++){
            volN[fam][3 * (qd_shift + n) + 0] = PhiK[kk] * PhiJ[jj] * EI[ii] * diff_eta_gma_ksi[fam][0];
            volN[fam][3 * (qd_shift + n) + 1] = PhiK[kk] * PhiJ[jj] * EI[ii] * diff_eta_gma_ksi[fam][1];
            volN[fam][3 * (qd_shift + n) + 2] = PhiK[kk] * PhiJ[jj] * EI[ii] * diff_eta_gma_ksi[fam][2];

            double EEI[3] = {EI[ii] * diff_eta_gma_ksi[fam][0],
                             EI[ii] * diff_eta_gma_ksi[fam][1],
                             EI[ii] * diff_eta_gma_ksi[fam][2]};

            double EEIJ[3] = {PhiJ[jj] * EEI[0], PhiJ[jj] * EEI[1], PhiJ[jj] * EEI[2]};

            double diff_PhiJ[3] = {diffPhiJ[jj] * diff_ksi_eta_gma[fam][0],
                                   diffPhiJ[jj] * diff_ksi_eta_gma[fam][1],
                                   diffPhiJ[jj] * diff_ksi_eta_gma[fam][2]};

            double diff_PhiK[3] = {diffPhiK[jj] * diff_gma_ksi_eta[fam][0],
                                   diffPhiJ[jj] * diff_gma_ksi_eta[fam][1],
                                   diffPhiJ[jj] * diff_gma_ksi_eta[fam][2]};

            auto diff_PhiK_cross_EEIJ = ref_hex.Cross_product(diff_PhiK, EEIJ);
            auto diff_PhiJ_cross_EEI = ref_hex.Cross_product(diff_PhiJ, EEI);

            curl_volN[fam][3 * (qd_shift + n) + 0] = PhiK[kk] * diff_PhiJ_cross_EEI[0] + diff_PhiK_cross_EEIJ[0];
            curl_volN[fam][3 * (qd_shift + n) + 1] = PhiK[kk] * diff_PhiJ_cross_EEI[1] + diff_PhiK_cross_EEIJ[1];
            curl_volN[fam][3 * (qd_shift + n) + 2] = PhiK[kk] * diff_PhiJ_cross_EEI[2] + diff_PhiK_cross_EEIJ[2];
            
            n++;
          }
          
        }
        
      }
      
    }
    
  }
  

  MoFEMFunctionReturnHot(0);
}

MoFEMErrorCode MoFEM::Hdiv_FaceShapeFunctions_ONHEX(int *face_nodes[6], 
                                                    int *p, 
                                                    double *N,
                                                    double *faceN[6],
                                                    double *div_faceN[6],
                                                    int nb_integration_pts){
  MoFEMFunctionBeginHot;
  RefHex ref_hex(N, nb_integration_pts);
  auto ksi = ref_hex.get_face_coords(face_nodes);
  auto diff_ksi = ref_hex.get_face_diff_coords(face_nodes);
  auto mu = ref_hex.get_face_affines();
  auto diff_mu = ref_hex.get_face_diff_affines();
  int pp[6][2] = {{p[0], p[1]},
                  {p[0], p[2]},
                  {p[1], p[2]},
                  {p[0], p[2]},
                  {p[1], p[2]},
                  {p[0], p[1]}};

  for (int qq = 0; qq < 6; qq++){
    for (int ff = 0; ff < 6; ff++) {
      double EI[pp[ff][0]];
      CHKERR Legendre_polynomials01(pp[ff][0] - 1, ksi[qq][2 * ff + 0], EI);

      double EJ[pp[ff][1]];
      CHKERR Legendre_polynomials01(pp[ff][1] - 1, ksi[qq][2 * ff + 1], EJ);
      int qd_shift = pp[ff][0] * pp[ff][1] * qq;
      int n = 0;
      for (int ii = 0; ii < pp[ff][0]; ii++){
        for (int jj = 0; jj < pp[ff][1]; jj++){
          double EEI[3] = {EI[ii] * diff_ksi[ff][0][0], EI[ii] * diff_ksi[ff][0][1], EI[ii] * diff_ksi[ff][0][2]};
          double EEJ[3] = {EJ[jj] * diff_ksi[ff][1][0], EJ[jj] * diff_ksi[ff][1][1], EJ[jj] * diff_ksi[ff][1][2]};
          auto VIJ_square = ref_hex.Cross_product(EEI, EEJ);

          faceN[ff][3 * (qd_shift + n) + 0] = mu[qq][ff] * VIJ_square[0];
          faceN[ff][3 * (qd_shift + n) + 1] = mu[qq][ff] * VIJ_square[1];
          faceN[ff][3 * (qd_shift + n) + 2] = mu[qq][ff] * VIJ_square[2];

          div_faceN[ff][qd_shift + n] = diff_mu[ff][0] * VIJ_square[0] +
                                         diff_mu[ff][1] * VIJ_square[1] +
                                         diff_mu[ff][2] * VIJ_square[2];

          n++;
        }
        
      }
      
    }
  }
  MoFEMFunctionReturnHot(0);

}

MoFEMErrorCode MoFEM::Hdiv_InteriorShapeFunctions_ONHEX(int *p, double *N,
                                                 double *bubbleN[3],
                                                 double *div_bubbleN[3],
                                                 int nb_integration_pts){
  MoFEMFunctionBeginHot;
  RefHex ref_hex(N, nb_integration_pts);
  auto ksi = ref_hex.get_integrationPts();
  auto diff_ksi = ref_hex.get_volume_diff_coords();

  for (int qq = 0; qq < nb_integration_pts; qq++) {
    double ksi_eta_gma[3] = {ksi[qq][0], ksi[qq][1], ksi[qq][2]};
    double eta_gma_ksi[3] = {ksi[qq][1], ksi[qq][2], ksi[qq][0]};
    double gma_ksi_eta[3] = {ksi[qq][2], ksi[qq][0], ksi[qq][1]};

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

      int qd_shift = pqr[fam] * (qrp[fam] - 1) * (rpq[fam] - 1) * qq;
      int n = 0;
      for (int ii = 0; ii < qrp[fam]; ii++){
        for (int jj = 0; jj < rpq[fam]; jj++){
          for (int kk = 0; kk < pqr[fam] - 1; kk++){
              double EEI[3] = {EI[ii] * diff_eta_gma_ksi[fam][0], EI[ii] * diff_eta_gma_ksi[fam][1], EI[ii] * diff_eta_gma_ksi[fam][2]};
              double EEJ[3] = {EJ[jj] * diff_gma_ksi_eta[fam][0], EJ[jj] * diff_gma_ksi_eta[fam][1], EJ[jj] * diff_gma_ksi_eta[fam][2]};

              auto VIJ_square = ref_hex.Cross_product(EEI, EEJ);

              bubbleN[fam][3 * (qd_shift + n) + 0] = PhiK[kk] * VIJ_square[0];
              bubbleN[fam][3 * (qd_shift + n) + 1] = PhiK[kk] * VIJ_square[2];
              bubbleN[fam][3 * (qd_shift + n) + 2] = PhiK[kk] * VIJ_square[2];

              div_bubbleN[fam][qd_shift + n] = diffPhiK[kk] * diff_ksi_eta_gma[fam][0] * VIJ_square[0] +
                                               diffPhiK[kk] * diff_ksi_eta_gma[fam][1] * VIJ_square[1] +
                                               diffPhiK[kk] * diff_ksi_eta_gma[fam][2] * VIJ_square[2];
              n++;
          }
          
        }
        
      }
       
    }
  }

    MoFEMFunctionReturnHot(0);
}
