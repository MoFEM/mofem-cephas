/** \file EdgeQuadHexPolynomials.cpp

  \brief Implementation of hierarchical Edge, Quad, and Hex shape functions of
  type H1, Hcurl, Hdiv
*/

using namespace MoFEM;

struct RefHex {
  RefHex()
      : vertices{{-1.0, -1.0, -1.0}, {1.0, -1.0, -1.0}, {1.0, 1.0, -1.0},
                 {-1.0, 1.0, -1.0},  {-1.0, -1.0, 1.0}, {1.0, -1.0, 1.0},
                 {1.0, 1.0, 1.0},    {-1.0, 1.0, 1.0}},
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

  void get_volume_diff_coords(double (&Nq_diff)[8][3],
                              double (&volume_diff_coords)[3][3]) {
    for (int c1 = 0; c1 < 3; c1++)
      for (int c2 = 0; c2 < 3; c2++)
        volume_diff_coords[c1][c2] = 0.0;

    for (int cc = 0; cc < 3; cc++) {
      for (int vv = 0; vv < 8; vv++) {
        volume_diff_coords[cc][0] += Nq_diff[vv][0] * vertices[vv][cc];
        volume_diff_coords[cc][1] += Nq_diff[vv][1] * vertices[vv][cc];
        volume_diff_coords[cc][2] += Nq_diff[vv][2] * vertices[vv][cc];
      }
    }
  }

  void get_edge_affines(double (&Nq)[8], double (&edge_affines)[12][2]) {

    double face0 = Nq[faces[0][0]] + Nq[faces[0][1]] + Nq[faces[0][2]] +
                   Nq[faces[0][3]]; // gmma_0
    double face1 = Nq[faces[1][0]] + Nq[faces[1][1]] + Nq[faces[1][2]] +
                   Nq[faces[1][3]]; // eta_0
    double face2 = Nq[faces[2][0]] + Nq[faces[2][1]] + Nq[faces[2][2]] +
                   Nq[faces[2][3]]; // ksi_1
    double face3 = Nq[faces[3][0]] + Nq[faces[3][1]] + Nq[faces[3][2]] +
                   Nq[faces[3][3]]; // eta_1
    double face4 = Nq[faces[4][0]] + Nq[faces[4][1]] + Nq[faces[4][2]] +
                   Nq[faces[4][3]]; // ksi_0
    double face5 = Nq[faces[5][0]] + Nq[faces[5][1]] + Nq[faces[5][2]] +
                   Nq[faces[5][3]]; // gmma_1

    double affines[12][2] = {{face1, face0}, {face2, face0}, {face3, face0},
                             {face4, face0}, {face1, face4}, {face1, face2},
                             {face3, face2}, {face3, face4}, {face1, face5},
                             {face2, face5}, {face3, face5}, {face4, face5}};

    for (int ee = 0; ee < 12; ee++)
      for (int comp = 0; comp < 3; comp++)
        edge_affines[ee][comp] = affines[ee][comp];
  }

  void get_edge_diff_affines(double (&Nq_diff)[8][3],
                             double (&edge_diff_affines)[12][2][3]) {

    double face0_diff[3] = {
        Nq_diff[faces[0][0]][0] + Nq_diff[faces[0][1]][0] +
            Nq_diff[faces[0][2]][0] + Nq_diff[faces[0][3]][0],
        Nq_diff[faces[0][0]][1] + Nq_diff[faces[0][1]][1] +
            Nq_diff[faces[0][2]][1] + Nq_diff[faces[0][3]][1],
        Nq_diff[faces[0][0]][2] + Nq_diff[faces[0][1]][2] +
            Nq_diff[faces[0][2]][2] + Nq_diff[faces[0][3]][2]}; // gma_0_diff

    double face1_diff[3] = {
        Nq_diff[faces[1][0]][0] + Nq_diff[faces[1][1]][0] +
            Nq_diff[faces[1][2]][0] + Nq_diff[faces[1][3]][0],
        Nq_diff[faces[1][0]][1] + Nq_diff[faces[1][1]][1] +
            Nq_diff[faces[1][2]][1] + Nq_diff[faces[1][3]][1],
        Nq_diff[faces[1][0]][2] + Nq_diff[faces[1][1]][2] +
            Nq_diff[faces[1][2]][2] + Nq_diff[faces[1][3]][2]}; // eta_0_diff

    double face2_diff[3] = {
        Nq_diff[faces[2][0]][0] + Nq_diff[faces[2][1]][0] +
            Nq_diff[faces[2][2]][0] + Nq_diff[faces[2][3]][0],
        Nq_diff[faces[2][0]][1] + Nq_diff[faces[2][1]][1] +
            Nq_diff[faces[2][2]][1] + Nq_diff[faces[2][3]][1],
        Nq_diff[faces[2][0]][2] + Nq_diff[faces[2][1]][2] +
            Nq_diff[faces[2][2]][2] + Nq_diff[faces[2][3]][2]}; // ksi_1_diff

    double face3_diff[3] = {
        Nq_diff[faces[3][0]][0] + Nq_diff[faces[3][1]][0] +
            Nq_diff[faces[3][2]][0] + Nq_diff[faces[3][3]][0],
        Nq_diff[faces[3][0]][1] + Nq_diff[faces[3][1]][1] +
            Nq_diff[faces[3][2]][1] + Nq_diff[faces[3][3]][1],
        Nq_diff[faces[3][0]][2] + Nq_diff[faces[3][1]][2] +
            Nq_diff[faces[3][2]][2] + Nq_diff[faces[3][3]][2]}; // eta_1_diff

    double face4_diff[3] = {
        Nq_diff[faces[4][0]][0] + Nq_diff[faces[4][1]][0] +
            Nq_diff[faces[4][2]][0] + Nq_diff[faces[4][3]][0],
        Nq_diff[faces[4][0]][1] + Nq_diff[faces[4][1]][1] +
            Nq_diff[faces[4][2]][1] + Nq_diff[faces[4][3]][1],
        Nq_diff[faces[4][0]][2] + Nq_diff[faces[4][1]][2] +
            Nq_diff[faces[4][2]][2] + Nq_diff[faces[4][3]][2]}; // ksi_0_diff

    double face5_diff[3] = {
        Nq_diff[faces[5][0]][0] + Nq_diff[faces[5][1]][0] +
            Nq_diff[faces[5][2]][0] + Nq_diff[faces[5][3]][0],
        Nq_diff[faces[5][0]][1] + Nq_diff[faces[5][1]][1] +
            Nq_diff[faces[5][2]][1] + Nq_diff[faces[5][3]][1],
        Nq_diff[faces[5][0]][2] + Nq_diff[faces[5][1]][2] +
            Nq_diff[faces[5][2]][2] + Nq_diff[faces[5][3]][2]}; // gmma_1_diff

    double *affines_diff[12][2] = {
        {face1_diff, face0_diff}, {face2_diff, face0_diff},
        {face3_diff, face0_diff}, {face4_diff, face0_diff},
        {face1_diff, face4_diff}, {face1_diff, face2_diff},
        {face3_diff, face2_diff}, {face3_diff, face4_diff},
        {face1_diff, face5_diff}, {face2_diff, face5_diff},
        {face3_diff, face5_diff}, {face4_diff, face5_diff}};

    for (int ee = 0; ee < 12; ee++) {
      for (int c = 0; c < 3; c++) {
        edge_diff_affines[ee][0][c] = affines_diff[ee][0][c];
        edge_diff_affines[ee][1][c] = affines_diff[ee][1][c];
      }
    }
  }

  void get_edge_coords(int *sense, double (&Nq)[8], double (&edge_coords)[12]) {

    double face0 = Nq[faces[0][0]] + Nq[faces[0][1]] + Nq[faces[0][2]] +
                   Nq[faces[0][3]]; // gmma_0
    double face1 = Nq[faces[1][0]] + Nq[faces[1][1]] + Nq[faces[1][2]] +
                   Nq[faces[1][3]]; // eta_0
    double face2 = Nq[faces[2][0]] + Nq[faces[2][1]] + Nq[faces[2][2]] +
                   Nq[faces[2][3]]; // ksi_1
    double face3 = Nq[faces[3][0]] + Nq[faces[3][1]] + Nq[faces[3][2]] +
                   Nq[faces[3][3]]; // eta_1
    double face4 = Nq[faces[4][0]] + Nq[faces[4][1]] + Nq[faces[4][2]] +
                   Nq[faces[4][3]]; // ksi_0
    double face5 = Nq[faces[5][0]] + Nq[faces[5][1]] + Nq[faces[5][2]] +
                   Nq[faces[5][3]]; // gmma_1

    double free_edge_coords[12] = {face2 - face4, face3 - face1, face4 - face2,
                                   face1 - face3, face5 - face0, face5 - face0,
                                   face5 - face0, face5 - face0, face2 - face4,
                                   face3 - face1, face4 - face2, face1 - face3};

    for (int ee = 0; ee < 12; ee++) {
      edge_coords[ee] = free_edge_coords[ee] * (double)sense[ee];
    }
  }

  void get_edge_diff_coords(int *sense, double (&Nq_diff)[8][3],
                            double (&edge_diff_coords)[12][3]) {

    double face0_diff[3] = {
        Nq_diff[faces[0][0]][0] + Nq_diff[faces[0][1]][0] +
            Nq_diff[faces[0][2]][0] + Nq_diff[faces[0][3]][0],
        Nq_diff[faces[0][0]][1] + Nq_diff[faces[0][1]][1] +
            Nq_diff[faces[0][2]][1] + Nq_diff[faces[0][3]][1],
        Nq_diff[faces[0][0]][2] + Nq_diff[faces[0][1]][2] +
            Nq_diff[faces[0][2]][2] + Nq_diff[faces[0][3]][2]}; // gma_0_diff

    double face1_diff[3] = {
        Nq_diff[faces[1][0]][0] + Nq_diff[faces[1][1]][0] +
            Nq_diff[faces[1][2]][0] + Nq_diff[faces[1][3]][0],
        Nq_diff[faces[1][0]][1] + Nq_diff[faces[1][1]][1] +
            Nq_diff[faces[1][2]][1] + Nq_diff[faces[1][3]][1],
        Nq_diff[faces[1][0]][2] + Nq_diff[faces[1][1]][2] +
            Nq_diff[faces[1][2]][2] + Nq_diff[faces[1][3]][2]};
    ; // eta_0_diff

    double face2_diff[3] = {
        Nq_diff[faces[2][0]][0] + Nq_diff[faces[2][1]][0] +
            Nq_diff[faces[2][2]][0] + Nq_diff[faces[2][3]][0],
        Nq_diff[faces[2][0]][1] + Nq_diff[faces[2][1]][1] +
            Nq_diff[faces[2][2]][1] + Nq_diff[faces[2][3]][1],
        Nq_diff[faces[2][0]][2] + Nq_diff[faces[2][1]][2] +
            Nq_diff[faces[2][2]][2] + Nq_diff[faces[2][3]][2]}; // ksi_1_diff

    double face3_diff[3] = {
        Nq_diff[faces[3][0]][0] + Nq_diff[faces[3][1]][0] +
            Nq_diff[faces[3][2]][0] + Nq_diff[faces[3][3]][0],
        Nq_diff[faces[3][0]][1] + Nq_diff[faces[3][1]][1] +
            Nq_diff[faces[3][2]][1] + Nq_diff[faces[3][3]][1],
        Nq_diff[faces[3][0]][2] + Nq_diff[faces[3][1]][2] +
            Nq_diff[faces[3][2]][2] + Nq_diff[faces[3][3]][2]}; // eta_1_diff

    double face4_diff[3] = {
        Nq_diff[faces[4][0]][0] + Nq_diff[faces[4][1]][0] +
            Nq_diff[faces[4][2]][0] + Nq_diff[faces[4][3]][0],
        Nq_diff[faces[4][0]][1] + Nq_diff[faces[4][1]][1] +
            Nq_diff[faces[4][2]][1] + Nq_diff[faces[4][3]][1],
        Nq_diff[faces[4][0]][2] + Nq_diff[faces[4][1]][2] +
            Nq_diff[faces[4][2]][2] + Nq_diff[faces[4][3]][2]}; // ksi_0_diff

    double face5_diff[3] = {
        Nq_diff[faces[5][0]][0] + Nq_diff[faces[5][1]][0] +
            Nq_diff[faces[5][2]][0] + Nq_diff[faces[5][3]][0],
        Nq_diff[faces[5][0]][1] + Nq_diff[faces[5][1]][1] +
            Nq_diff[faces[5][2]][1] + Nq_diff[faces[5][3]][1],
        Nq_diff[faces[5][0]][2] + Nq_diff[faces[5][1]][2] +
            Nq_diff[faces[5][2]][2] + Nq_diff[faces[5][3]][2]}; // gmma_1_diff

    double f0_diff[3] = {face2_diff[0] - face4_diff[0],
                         face2_diff[1] - face4_diff[1],
                         face2_diff[2] - face4_diff[2]};
    double f1_diff[3] = {face3_diff[0] - face1_diff[0],
                         face3_diff[1] - face1_diff[1],
                         face3_diff[2] - face1_diff[2]};
    double f2_diff[3] = {face4_diff[0] - face2_diff[0],
                         face4_diff[1] - face2_diff[1],
                         face4_diff[2] - face2_diff[2]};
    double f3_diff[3] = {face1_diff[0] - face3_diff[0],
                         face1_diff[1] - face3_diff[1],
                         face1_diff[2] - face3_diff[2]};
    double f4_diff[3] = {face5_diff[0] - face0_diff[0],
                         face5_diff[1] - face0_diff[1],
                         face5_diff[2] - face0_diff[2]};

    double *free_edge_coords[12] = {f0_diff, f1_diff, f2_diff, f3_diff,
                                    f4_diff, f4_diff, f4_diff, f4_diff,
                                    f0_diff, f1_diff, f2_diff, f3_diff};

    for (int ee = 0; ee < 12; ee++) {
      for (int cc = 0; cc < 3; cc++) {
        edge_diff_coords[ee][cc] = free_edge_coords[ee][cc] * (double)sense[ee];
      }
    }
  }

  void get_face_affines(double (&Nq)[8], double (&face_affines)[6]) {

    double face0 = Nq[faces[0][0]] + Nq[faces[0][1]] + Nq[faces[0][2]] +
                   Nq[faces[0][3]]; // gmma_0
    double face1 = Nq[faces[1][0]] + Nq[faces[1][1]] + Nq[faces[1][2]] +
                   Nq[faces[1][3]]; // eta_0
    double face2 = Nq[faces[2][0]] + Nq[faces[2][1]] + Nq[faces[2][2]] +
                   Nq[faces[2][3]]; // ksi_1
    double face3 = Nq[faces[3][0]] + Nq[faces[3][1]] + Nq[faces[3][2]] +
                   Nq[faces[3][3]]; // eta_1
    double face4 = Nq[faces[4][0]] + Nq[faces[4][1]] + Nq[faces[4][2]] +
                   Nq[faces[4][3]]; // ksi_0
    double face5 = Nq[faces[5][0]] + Nq[faces[5][1]] + Nq[faces[5][2]] +
                   Nq[faces[5][3]]; // gmma_1

    double faceAffine[6] = {face0, face1, face2, face3, face4, face5};
    for (int ff = 0; ff != 6; ff++)
      face_affines[ff] = faceAffine[ff];
  }

  void get_face_diff_affines(double (&Nq_diff)[8][3],
                             double (&face_diff_affines)[6][3]) {

    double face0_diff[3] = {
        Nq_diff[faces[0][0]][0] + Nq_diff[faces[0][1]][0] +
            Nq_diff[faces[0][2]][0] + Nq_diff[faces[0][3]][0],
        Nq_diff[faces[0][0]][1] + Nq_diff[faces[0][1]][1] +
            Nq_diff[faces[0][2]][1] + Nq_diff[faces[0][3]][1],
        Nq_diff[faces[0][0]][2] + Nq_diff[faces[0][1]][2] +
            Nq_diff[faces[0][2]][2] + Nq_diff[faces[0][3]][2]}; // gma_0_diff

    double face1_diff[3] = {
        Nq_diff[faces[1][0]][0] + Nq_diff[faces[1][1]][0] +
            Nq_diff[faces[1][2]][0] + Nq_diff[faces[1][3]][0],
        Nq_diff[faces[1][0]][1] + Nq_diff[faces[1][1]][1] +
            Nq_diff[faces[1][2]][1] + Nq_diff[faces[1][3]][1],
        Nq_diff[faces[1][0]][2] + Nq_diff[faces[1][1]][2] +
            Nq_diff[faces[1][2]][2] + Nq_diff[faces[1][3]][2]};
    ; // eta_0_diff

    double face2_diff[3] = {
        Nq_diff[faces[2][0]][0] + Nq_diff[faces[2][1]][0] +
            Nq_diff[faces[2][2]][0] + Nq_diff[faces[2][3]][0],
        Nq_diff[faces[2][0]][1] + Nq_diff[faces[2][1]][1] +
            Nq_diff[faces[2][2]][1] + Nq_diff[faces[2][3]][1],
        Nq_diff[faces[2][0]][2] + Nq_diff[faces[2][1]][2] +
            Nq_diff[faces[2][2]][2] + Nq_diff[faces[2][3]][2]}; // ksi_1_diff

    double face3_diff[3] = {
        Nq_diff[faces[3][0]][0] + Nq_diff[faces[3][1]][0] +
            Nq_diff[faces[3][2]][0] + Nq_diff[faces[3][3]][0],
        Nq_diff[faces[3][0]][1] + Nq_diff[faces[3][1]][1] +
            Nq_diff[faces[3][2]][1] + Nq_diff[faces[3][3]][1],
        Nq_diff[faces[3][0]][2] + Nq_diff[faces[3][1]][2] +
            Nq_diff[faces[3][2]][2] + Nq_diff[faces[3][3]][2]}; // eta_1_diff

    double face4_diff[3] = {
        Nq_diff[faces[4][0]][0] + Nq_diff[faces[4][1]][0] +
            Nq_diff[faces[4][2]][0] + Nq_diff[faces[4][3]][0],
        Nq_diff[faces[4][0]][1] + Nq_diff[faces[4][1]][1] +
            Nq_diff[faces[4][2]][1] + Nq_diff[faces[4][3]][1],
        Nq_diff[faces[4][0]][2] + Nq_diff[faces[4][1]][2] +
            Nq_diff[faces[4][2]][2] + Nq_diff[faces[4][3]][2]}; // ksi_0_diff

    double face5_diff[3] = {
        Nq_diff[faces[5][0]][0] + Nq_diff[faces[5][1]][0] +
            Nq_diff[faces[5][2]][0] + Nq_diff[faces[5][3]][0],
        Nq_diff[faces[5][0]][1] + Nq_diff[faces[5][1]][1] +
            Nq_diff[faces[5][2]][1] + Nq_diff[faces[5][3]][1],
        Nq_diff[faces[5][0]][2] + Nq_diff[faces[5][1]][2] +
            Nq_diff[faces[5][2]][2] + Nq_diff[faces[5][3]][2]}; // gmma_1_diff

    double *free_face[6] = {face0_diff, face1_diff, face2_diff,
                            face3_diff, face4_diff, face5_diff};

    for (int ff = 0; ff < 6; ff++) {
      for (int cc = 0; cc < 3; cc++) {
        face_diff_affines[ff][cc] = free_face[ff][cc];
      }
    }
  }

  void get_face_coords(int *face_nodes[6], double (&Nq)[8],
                       double (&face_coords)[6][2]) {

    int par_face_nodes[6][8] = {
        {4, 5, 6, 7, 0, 1, 2, 3}, {3, 2, 1, 0, 7, 6, 5, 4},
        {1, 0, 3, 2, 5, 4, 7, 6}, {3, 2, 1, 0, 7, 6, 5, 4},
        {1, 0, 3, 2, 5, 4, 7, 6}, {4, 5, 6, 7, 0, 1, 2, 3}};

    double face_vertices[4][2] = {
        {-1.0, -1.0}, {1.0, -1.0}, {1.0, 1.0}, {-1.0, 1.0}};

    for (int ff = 0; ff != 6; ff++) {
      int n00 = par_face_nodes[ff][0];
      int n01 = par_face_nodes[ff][4];

      int n10 = par_face_nodes[ff][1];
      int n11 = par_face_nodes[ff][5];

      int n20 = par_face_nodes[ff][2];
      int n21 = par_face_nodes[ff][6];

      int n30 = par_face_nodes[ff][3];
      int n31 = par_face_nodes[ff][7];

      double N0 = Nq[n00] + Nq[n01];
      double N1 = Nq[n10] + Nq[n11];
      double N2 = Nq[n20] + Nq[n21];
      double N3 = Nq[n30] + Nq[n31];

      int n0 = face_nodes[ff][0];
      int n1 = face_nodes[ff][1];
      int n2 = face_nodes[ff][2];
      int n3 = face_nodes[ff][3];

      for (int cc = 0; cc < 2; cc++) {
        face_coords[ff][cc] =
            N0 * face_vertices[n0][cc] + N1 * face_vertices[n1][cc] +
            N2 * face_vertices[n2][cc] + N3 * face_vertices[n3][cc];
      }
    }
  }

  void get_face_diff_coords(int *face_nodes[6], double (&Nq_diff)[8][3],
                            double (&face_diff_coords)[6][2][3]) {

    int par_face_nodes[6][8] = {
        {4, 5, 6, 7, 0, 1, 2, 3}, {3, 2, 1, 0, 7, 6, 5, 4},
        {1, 0, 3, 2, 5, 4, 7, 6}, {3, 2, 1, 0, 7, 6, 5, 4},
        {1, 0, 3, 2, 5, 4, 7, 6}, {4, 5, 6, 7, 0, 1, 2, 3}};

    double face_vertices[4][2] = {
        {-1.0, -1.0}, {1.0, -1.0}, {1.0, 1.0}, {-1.0, 1.0}};

    for (int ff = 0; ff != 6; ff++) {
      int n00 = par_face_nodes[ff][0];
      int n01 = par_face_nodes[ff][4];

      int n10 = par_face_nodes[ff][1];
      int n11 = par_face_nodes[ff][5];

      int n20 = par_face_nodes[ff][2];
      int n21 = par_face_nodes[ff][6];

      int n30 = par_face_nodes[ff][3];
      int n31 = par_face_nodes[ff][7];

      double N0_diff[3] = {Nq_diff[n00][0] + Nq_diff[n01][0],
                           Nq_diff[n00][1] + Nq_diff[n01][1],
                           Nq_diff[n00][2] + Nq_diff[n01][2]};
      double N1_diff[3] = {Nq_diff[n10][0] + Nq_diff[n11][0],
                           Nq_diff[n10][1] + Nq_diff[n11][1],
                           Nq_diff[n10][2] + Nq_diff[n11][2]};
      double N2_diff[3] = {Nq_diff[n20][0] + Nq_diff[n21][0],
                           Nq_diff[n20][1] + Nq_diff[n21][1],
                           Nq_diff[n20][2] + Nq_diff[n21][2]};
      double N3_diff[3] = {Nq_diff[n30][0] + Nq_diff[n31][0],
                           Nq_diff[n30][1] + Nq_diff[n31][1],
                           Nq_diff[n30][2] + Nq_diff[n31][2]};

      int n0 = face_nodes[ff][0];
      int n1 = face_nodes[ff][1];
      int n2 = face_nodes[ff][2];
      int n3 = face_nodes[ff][3];

      for (int cc = 0; cc < 2; cc++) {
        for (int dd = 0; dd < 3; dd++) {
          face_diff_coords[ff][cc][dd] = N0_diff[dd] * face_vertices[n0][cc] +
                                         N1_diff[dd] * face_vertices[n1][cc] +
                                         N2_diff[dd] * face_vertices[n2][cc] +
                                         N3_diff[dd] * face_vertices[n3][cc];
        }
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
          }
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
    for (int i = 2; i != p + 1; i++) {
      double factor = 1.0 / (2.0 * (2.0 * (double)i - 1.0));
      L[i - 2] = factor * (l[i] - l[i - 2]);
      diffL[i - 2] = l[i - 1];
    }
  }

  MoFEMFunctionReturnHot(0);
}

/*
    0 *------------* 1
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

namespace DemkowiczHexAndQuad {
static inline void get_ksi_hex(int shift, double *N, double *N_diff,
                               double ksi[3], double diff_ksi[3][2]) {

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
    for (auto n : ksi_nodes[0][i]) {
      for (auto d = 0; d != 3; ++d) {
        diff_ksi[i][d] += N_diff[3 * shift + 3 * n + d];
      }
    }
    for (auto n : ksi_nodes[1][i]) {
      for (auto d = 0; d != 3; ++d) {
        diff_ksi[i][d] -= N_diff[3 * shift + 3 * n + d];
      }
    }
  }
};
} // namespace DemkowiczHexAndQuad

MoFEMErrorCode MoFEM::DemkowiczHexAndQuad::H1_EdgeShapeFunctions_ONHEX(
    int *sense, int *p, double *N, double *N_diff, double *edgeN[12],
    double *diff_edgeN[12], int nb_integration_pts) {
  MoFEMFunctionBeginHot;
  RefHex ref_hex;

  for (int q = 0; q != nb_integration_pts; q++) {
    // General *******************************
    int shift = 8 * q;
    double quad_coords[3];
    double Nq[8];
    double Nq_diff[8][3];
    for (int vv = 0; vv < 8; vv++) {
      Nq[vv] = N[shift + vv];
      Nq_diff[vv][0] = N_diff[3 * (shift + vv) + 0];
      Nq_diff[vv][1] = N_diff[3 * (shift + vv) + 1];
      Nq_diff[vv][2] = N_diff[3 * (shift + vv) + 2];
    }

    double mu[12][2];
    ref_hex.get_edge_affines(Nq, mu);

    double diff_mu[12][2][3];
    ref_hex.get_edge_diff_affines(Nq_diff, diff_mu);

    double ksi[12];
    ref_hex.get_edge_coords(sense, Nq, ksi);
    double diff_ksi[12][3];
    ref_hex.get_edge_diff_coords(sense, Nq_diff, diff_ksi);

    int pp[12] = {p[0], p[1], p[0], p[1], p[2], p[2],
                  p[2], p[2], p[0], p[1], p[0], p[1]};

    for (int e = 0; e != 12; e++) {
      double L[pp[e] + 2];
      double diffL[3 * (pp[e] + 2)];

      CHKERR Lobatto_polynomials(pp[e] + 1, ksi[e], diff_ksi[e], L, diffL, 3);

      int qd_shift = (pp[e] - 1) * q;

      for (int n = 0; n != pp[e] - 1; n++) {
        edgeN[e][qd_shift + n] = mu[e][0] * mu[e][1] * L[n + 2];
        for (int d = 0; d != 3; ++d) {
          diff_edgeN[e][3 * (qd_shift + n) + d] =
              diff_mu[e][0][d] * mu[e][1] * L[n + 2] + mu[e][0] +
              diff_mu[e][1][d] + L[n + 2] +
              mu[e][0] * mu[e][1] * diffL[d * (p[e] + 2) + n + 2];
        }
      }
    }
  }
  MoFEMFunctionReturnHot(0);
}
// TODO
MoFEMErrorCode MoFEM::DemkowiczHexAndQuad::H1_FaceShapeFunctions_ONHEX(
    int *face_nodes[6], int *p, double *N, double *N_diff, double *faceN[6],
    double *diff_faceN[6], int nb_integration_pts) {
  MoFEMFunctionBeginHot;

  RefHex ref_hex;

  for (int q = 0; q != nb_integration_pts; q++) {
    // general ******************************
    int shift = 8 * q;
    double quad_coords[3];
    double Nq[8];
    double Nq_diff[8][3];
    for (int vv = 0; vv < 8; vv++) {
      Nq[vv] = N[shift + vv];
      Nq_diff[vv][0] = N_diff[3 * (shift + vv) + 0];
      Nq_diff[vv][1] = N_diff[3 * (shift + vv) + 1];
      Nq_diff[vv][2] = N_diff[3 * (shift + vv) + 2];
    }

    double mu[6];
    ref_hex.get_face_affines(Nq, mu);
    double diff_mu[6][3];
    ref_hex.get_face_diff_affines(Nq_diff, diff_mu);

    double ksi[6][2];
    ref_hex.get_face_coords(face_nodes, Nq, ksi);
    double diff_ksi[6][2][3];
    ref_hex.get_face_diff_coords(face_nodes, Nq_diff, diff_ksi);

    int pp[6][2] = {{p[0], p[1]}, {p[0], p[2]}, {p[1], p[2]},
                    {p[0], p[2]}, {p[1], p[2]}, {p[0], p[1]}};

    for (int face = 0; face != 6; face++) {
      int pq[2] = {pp[face][0], pp[face][1]};

      double ksi0 = ksi[face][0];
      double ksi1 = ksi[face][1];

      double diff_ksi0[3] = {diff_ksi[face][0][0], diff_ksi[face][0][1],
                             diff_ksi[face][0][2]};
      double diff_ksi1[3] = {diff_ksi[face][1][0], diff_ksi[face][1][1],
                             diff_ksi[face][1][2]};

      double L0[pq[0] + 1];
      double diffL0[3 * (pq[0] + 2)];
      CHKERR Lobatto_polynomials(pq[0] + 1, ksi0, diff_ksi0, L0, diffL0, 3);
      double L1[pq[1] + 2];
      double diffL1[3 * (pq[1] + 2)];
      CHKERR Lobatto_polynomials(pq[1] + 1, ksi1, diff_ksi1, L1, diffL1, 3);

      int permute[(pq[0] - 1) * (pq[1] - 1)][3];
      CHKERR MonomOrdering(permute, pq[0] - 2, pq[1] - 2);

      int qd_shift = (pq[0] - 1) * (pq[1] - 1) * q;
      int n = 0;
      for (; n != (pq[0] - 1) * (pq[1] - 1); n++) {
        int s1 = permute[n][0];
        int s2 = permute[n][1];
        faceN[face][qd_shift + n] = mu[face] * L0[s1 + 2] * L1[s2 + 2];
        for (int d = 0; d != 3; ++d) {
          diff_faceN[face][3 * (qd_shift + n) + d] =
              diff_mu[face][0] * L0[s1 + 2] * L1[s2 + 2] +
              mu[face] * diffL0[d * (pq[0] + 2) + s1 + 2] * L1[s2] +
              mu[face] * L0[s1 + 2] * diffL1[d * (pq[1] + 2) + s2 + 2];
        }
      }
    }
  }
  MoFEMFunctionReturnHot(0);
}

MoFEMErrorCode MoFEM::DemkowiczHexAndQuad::H1_InteriorShapeFunctions_ONHEX(
    int *p, double *N, double *N_diff, double *faceN, double *diff_faceN,
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
    double Nq_diff[8][3];
    for (int vv = 0; vv < 8; vv++) {
      Nq[vv] = N[shift + vv];
      Nq_diff[vv][0] = N_diff[3 * (shift + vv) + 0];
      Nq_diff[vv][1] = N_diff[3 * (shift + vv) + 1];
      Nq_diff[vv][2] = N_diff[3 * (shift + vv) + 2];
    }
    ref_hex.get_volume_coords(Nq, ksi);
    // ****************************************

    double diff_ksi[3][3];

    ref_hex.get_volume_diff_coords(Nq_diff, diff_ksi);

    double L0[p[0] + 2];
    double diffL0[3 * (p[0] + 2)];
    double L1[p[1] + 2];
    double diffL1[3 * (p[1] + 2)];
    double L2[p[2] + 2];
    double diffL2[3 * (p[2] + 2)];

    CHKERR Lobatto_polynomials(p[0] + 1, ksi[0], diff_ksi[0], L0, diffL0, 3);
    CHKERR Lobatto_polynomials(p[1] + 1, ksi[1], diff_ksi[1], L1, diffL1, 3);
    CHKERR Lobatto_polynomials(p[2] + 1, ksi[2], diff_ksi[2], L2, diffL2, 3);

    // cout << "In Face H1" << endl;

    int qd_shift = (p[0] - 1) * (p[1] - 1) * (p[2] - 1) * q;
    int n = 0;
    for (; n != (p[0] - 1) * (p[1] - 1) * (p[2] - 1); n++) {
      int s1 = permute[n][0];
      int s2 = permute[n][1];
      int s3 = permute[n][2];

      faceN[qd_shift + n] = L0[s1] * L1[s2] * L2[s3];
      for (int d = 0; d != 3; ++d) {
        diff_faceN[3 * (qd_shift + n) + 0] =
            diffL0[d * (p[0] + 2) + s1 + 2] * L1[s2 + 2] * L2[s3 + 2] +
            L0[s1 + 2] * diffL1[d * (p[1] + 2) + s2 + 2] * L2[s3 + 2] +
            L0[s1 + 2] * L1[s2 + 2] * diffL2[d * (p[2] + 2) + s3 + 2];
      }
    }
  }
  MoFEMFunctionReturnHot(0);
}

MoFEMErrorCode MoFEM::DemkowiczHexAndQuad::L2_InteriorShapeFunctions_ONHEX(
    const int *p, double *N, double *N_diff, double *volN, double *diff_volN,
    int nb_integration_pts) {
  MoFEMFunctionBeginHot;
  RefHex ref_hex;
  int permute[p[0] * p[1] * p[2]][3];
  CHKERR MonomOrdering(permute, p[0] - 1, p[1] - 1, p[2] - 1);

  for (int qq = 0; qq != nb_integration_pts; qq++) {

    int shift = 8 * qq;
    double ksi[3] = {0, 0, 0};
    double diff_ksi[3][2] = {{0, 0}, {0, 0}, {0, 0}};
    ::DemkowiczHexAndQuad::get_ksi_hex(shift, N, N_diff, ksi, diff_ksi);

    double P0[p[0]];
    double diffL0[3 * (p[0] + 2)];
    double P1[p[1]];
    double diffL1[3 * (p[1] + 2)];
    double P2[p[2]];
    double diffL2[3 * (p[2] + 2)];

    int qd_shift = qq * p[0] * p[1] * p[2];
    CHKERR Legendre_polynomials(p[0] + 1, ksi[0], diff_ksi[0], P0, diffL0, 3);

    CHKERR Legendre_polynomials(p[1] + 1, ksi[1], diff_ksi[1], P1, diffL1, 3);

    CHKERR Legendre_polynomials(p[2] + 1, ksi[2], diff_ksi[2], P2, diffL2, 3);

    int n = 0;
    for (; n != p[0] * p[1] * p[2]; n++) {
      const int ii = permute[n][0];
      const int jj = permute[n][1];
      const int kk = permute[n][2];
      volN[qd_shift + n] = P0[ii] * P1[jj] * P2[kk];
      for (int d = 0; d != 3; ++d) {
        diff_volN[3 * (qd_shift + n) + d] =
            diffL0[d * (p[0] + 2) + ii] * P1[jj] * P2[kk] +
            P0[ii] * diffL1[d * (p[1] + 2) + jj] * P2[kk] +
            P0[ii] * P1[jj] * diffL2[d * (p[2] + 2) + kk];
      }
    }
  }

  MoFEMFunctionReturnHot(0);
}

MoFEMErrorCode MoFEM::DemkowiczHexAndQuad::Hcurl_EdgeShapeFunctions_ONHEX(
    int *sense, int *p, double *N, double *N_diff, double *edgeN[12],
    double *diff_edgeN[12], int nb_integration_pts) {
  MoFEMFunctionBeginHot;
  RefHex ref_hex;

  for (int qq = 0; qq != nb_integration_pts; qq++) {
    // general ******************************
    int shift = 8 * qq;
    double Nq[8];
    double Nq_diff[8][3];
    for (int vv = 0; vv < 8; vv++) {
      Nq[vv] = N[shift + vv];
      Nq_diff[vv][0] = N_diff[3 * (shift + vv) + 0];
      Nq_diff[vv][1] = N_diff[3 * (shift + vv) + 1];
      Nq_diff[vv][2] = N_diff[3 * (shift + vv) + 2];
    }

    double mu[12][2];
    ref_hex.get_edge_affines(Nq, mu);

    double diff_mu[12][2][3];
    ref_hex.get_edge_diff_affines(Nq_diff, diff_mu);

    double ksi[12];
    ref_hex.get_edge_coords(sense, Nq, ksi);
    double diff_ksi[12][3];
    ref_hex.get_edge_diff_coords(sense, Nq_diff, diff_ksi);

    int pp[12] = {p[0], p[1], p[0], p[1], p[2], p[2],
                  p[2], p[2], p[0], p[1], p[0], p[1]};
    for (int ee = 0; ee != 12; ee++) {

      double L[pp[ee]];
      double diffL[3 * (pp[ee] + 2)];
      CHKERR Legendre_polynomials(pp[ee] - 1, ksi[ee], diff_ksi[ee], L, diffL,
                                  3);

      // CHKERR Legendre_polynomials01(pp[ee] - 1, ksi[ee], L);
      int qd_shift = pp[ee] * qq;
      double *t_n_ptr = &edgeN[ee][3 * qd_shift];
      double *t_diff_n_ptr = &diff_edgeN[ee][3 * 3 * qd_shift];
      auto t_n = getFTensor1FromPtr<3>(t_n_ptr);
      auto t_diff_n = getFTensor2FromPtr<3, 3>(t_diff_n_ptr);

      for (int ii = 0; ii != pp[ee]; ii++) {
        const double a = mu[ee][0] * mu[ee][1] * L[ii];

        const double d_a[] = {
            diff_mu[ee][0][0] * mu[ee][1] * L[ii] +
                mu[ee][0] * diff_mu[ee][1][0] * L[ii] +
                mu[ee][0] * mu[ee][1] * diffL[0 * pp[ee] + ii],

            diff_mu[ee][0][1] * mu[ee][1] * L[ii] +
                mu[ee][0] * diff_mu[ee][1][1] * L[ii] +
                mu[ee][0] * mu[ee][1] * diffL[1 * pp[ee] + ii],

            diff_mu[ee][0][2] * mu[ee][1] * L[ii] +
                mu[ee][0] * diff_mu[ee][1][2] * L[ii] +
                mu[ee][0] * mu[ee][1] * diffL[2 * pp[ee] + ii]};

        for (int d = 0; d != 3; ++d) {
          t_n(d) = 0.5 * a * diff_ksi[ee][d];
          for (int j = 0; j != 2; ++j) {
            t_diff_n(d, j) = 0.5 * d_a[j] * diff_ksi[ee][d];
          }
        }
        ++t_n;
        ++t_diff_n;
      }
    }
  }
  MoFEMFunctionReturnHot(0);
}

MoFEMErrorCode MoFEM::DemkowiczHexAndQuad::Hcurl_FaceShapeFunctions_ONHEX(
    int *face_nodes[6], int *p, double *N, double *N_diff, double *faceN[6][2],
    double *diff_faceN[6][2], int nb_integration_pts) {
  MoFEMFunctionBeginHot;
  RefHex ref_hex;
  for (int qq = 0; qq != nb_integration_pts; qq++) {
    // general ******************************
    int shift = 8 * qq;
    double quad_coords[3];
    double Nq[8];
    double Nq_diff[8][3];
    for (int vv = 0; vv < 8; vv++) {
      Nq[vv] = N[shift + vv];
      Nq_diff[vv][0] = N_diff[3 * (shift + vv) + 0];
      Nq_diff[vv][1] = N_diff[3 * (shift + vv) + 1];
      Nq_diff[vv][2] = N_diff[3 * (shift + vv) + 2];
    }

    double mu[6];
    ref_hex.get_face_affines(Nq, mu);
    double diff_mu[6][3];
    ref_hex.get_face_diff_affines(Nq_diff, diff_mu);

    double ksi[6][2];
    ref_hex.get_face_coords(face_nodes, Nq, ksi);
    double diff_ksi[6][2][3];
    ref_hex.get_face_diff_coords(face_nodes, Nq_diff, diff_ksi);

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

        double Phi[pq[fam] + 2];
        double diffPhi[pq[fam] + 2];

        CHKERR Lobatto_polynomials(pq[fam] + 1, ksi_eta[fam], diff_ksi_eta[fam],
                                   Phi, diffPhi, 3);

        double E[qp[fam]];
        double diffE[2 * qp[fam]];
        CHKERR Legendre_polynomials(qp[fam] - 1, eta_ksi[fam],
                                    diff_eta_ksi[fam], E, diffE, 3);

        int permute[(pq[fam] - 1) * qp[fam]][3];
        CHKERR MonomOrdering(permute, qp[fam] - 1, pq[fam] - 2);

        int qd_shift = (pq[fam] - 1) * qp[fam] * qq;
        double *t_n_ptr = &faceN[ff][fam][3 * qd_shift];
        double *t_diff_n_ptr = &diff_faceN[ff][fam][3 * 3 * qd_shift];
        auto t_n = getFTensor1FromPtr<3>(t_n_ptr);
        auto t_diff_n = getFTensor2FromPtr<3, 3>(t_diff_n_ptr);

        int n = 0;
        for (; n != (pq[fam] - 1) * qp[fam]; n++) {
          int i = permute[n][0];
          int j = permute[n][1];

          const double phi = Phi[j + 2];
          const double e = E[i];
          const double a = phi * e;
          const double d_a[] = {diffPhi[0 * (pq[fam] + 2) + j + 2] * e +
                                    phi * diffE[0 * qp[fam] + i],

                                diffPhi[1 * (pq[fam] + 2) + j + 2] * e +
                                    phi * diffE[1 * qp[fam] + i],

                                diffPhi[2 * (pq[fam] + 2) + j + 2] * e +
                                    phi * diffE[2 * qp[fam] + i]};

          for (int d = 0; d != 3; ++d) {
            t_n(d) = a * diff_eta_ksi[fam][d];
            for (int m = 0; m != 2; ++m) {
              t_diff_n(d, m) = d_a[m] * diff_eta_ksi[fam][d];
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

MoFEMErrorCode MoFEM::DemkowiczHexAndQuad::Hcurl_InteriorShapeFunctions_ONHEX(
    int *p, double *N, double *N_diff, double *volN[], double *diff_volN[],
    int nb_integration_pts) {
  MoFEMFunctionBeginHot;
  RefHex ref_hex;

  for (int qq = 0; qq < nb_integration_pts; qq++) {
    // general ******************************
    int shift = 8 * qq;
    double quad_coords[3];
    double Nq[8];
    double Nq_diff[8][3];
    for (int vv = 0; vv < 8; vv++) {
      Nq[vv] = N[shift + vv];
      Nq_diff[vv][0] = N_diff[3 * (shift + vv) + 0];
      Nq_diff[vv][1] = N_diff[3 * (shift + vv) + 1];
      Nq_diff[vv][2] = N_diff[3 * (shift + vv) + 2];
    }

    ref_hex.get_volume_coords(Nq, quad_coords);
    // ****************************************
    double *ksi = quad_coords;

    double diff_ksi[3][3];
    ref_hex.get_volume_diff_coords(Nq_diff, diff_ksi);

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
      double PhiJ[ppp + 2];
      double diffPhiJ[3 * (ppp + 2)];

      CHKERR Lobatto_polynomials(ppp + 1, ksi_eta_gma[fam],
                                 diff_ksi_eta_gma[fam], PhiJ, diffPhiJ, 3);

      int qqq = qrp[fam];
      double EI[qqq];
      double diffEI[3 * qqq];

      CHKERR Legendre_polynomials(qqq - 1, eta_gma_ksi[fam],
                                  diff_eta_gma_ksi[fam], EI, diffEI, 3);

      int rrr = rpq[fam];

      double PhiK[rrr + 2];
      double diffPhiK[3 * (rrr + 2)];

      CHKERR Lobatto_polynomials(rrr + 1, gma_ksi_eta[fam],
                                 diff_gma_ksi_eta[fam], PhiK, diffPhiK, 3);

      int permute[(ppp - 1) * qqq * (rrr - 1)][3];
      CHKERR MonomOrdering(permute, ppp - 2, qqq - 1, rrr - 2);

      int qd_shift = (ppp - 1) * qqq * (rrr - 1) * qq;
      double *t_n_ptr = &volN[fam][3 * qd_shift];
      double *t_diff_n_ptr = &diff_volN[fam][3 * 3 * qd_shift];
      auto t_n = getFTensor1FromPtr<3>(t_n_ptr);
      auto t_diff_n = getFTensor2FromPtr<3, 3>(t_diff_n_ptr);

      int n = 0;
      for (; n != (ppp - 1) * qqq * (rrr - 1); n++) {
        int ii = permute[n][0];
        int jj = permute[n][1];
        int kk = permute[n][2];

        const double phiK = PhiK[kk + 2];
        const double phiJ = PhiJ[jj + 2];
        const double e = EI[ii];
        const double a = phiJ * phiK * e;

        const double d_a[] = {diffPhiK[0 * (ppp + 2) + kk + 2] * phiJ * e +
                                  phiK * diffPhiJ[0 * (rrr + 2) + jj + 2] * e +
                                  phiK * phiJ * diffEI[0 * qqq + ii],

                              diffPhiK[1 * (ppp + 2) + kk + 2] * phiJ * e +
                                  phiK * diffPhiJ[1 * (rrr + 2) + jj + 2] * e +
                                  phiK * phiJ * diffEI[1 * qqq + ii],

                              diffPhiK[2 * (ppp + 2) + kk + 2] * phiJ * e +
                                  phiK * diffPhiJ[2 * (rrr + 2) + jj + 2] * e +
                                  phiK * phiJ * diffEI[2 * qqq + ii]};

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

  MoFEMFunctionReturnHot(0);
}

MoFEMErrorCode MoFEM::DemkowiczHexAndQuad::Hdiv_FaceShapeFunctions_ONHEX(
    int *face_nodes[6], int *p, double *N, double *N_diff, double *faceN[6],
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
    double Nq_diff[8][3];
    for (int vv = 0; vv < 8; vv++) {
      Nq[vv] = N[shift + vv];
      Nq_diff[vv][0] = N_diff[3 * (shift + vv) + 0];
      Nq_diff[vv][1] = N_diff[3 * (shift + vv) + 1];
      Nq_diff[vv][2] = N_diff[3 * (shift + vv) + 2];
    }

    double mu[6];
    ref_hex.get_face_affines(Nq, mu);
    double diff_mu[6][3];
    ref_hex.get_face_diff_affines(Nq_diff, diff_mu);

    double ksi[6][2];
    ref_hex.get_face_coords(face_nodes, Nq, ksi);
    double diff_ksi[6][2][3];
    ref_hex.get_face_diff_coords(face_nodes, Nq_diff, diff_ksi);

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
    int *p, double *N, double *N_diff, double *bubbleN[3],
    double *div_bubbleN[3], int nb_integration_pts) {
  MoFEMFunctionBeginHot;
  RefHex ref_hex;

  for (int qq = 0; qq < nb_integration_pts; qq++) {
    // general ******************************
    int shift = 8 * qq;
    double quad_coords[3];
    double Nq[8];
    double Nq_diff[8][3];
    for (int vv = 0; vv < 8; vv++) {
      Nq[vv] = N[shift + vv];
      Nq_diff[vv][0] = N_diff[3 * (shift + vv) + 0];
      Nq_diff[vv][1] = N_diff[3 * (shift + vv) + 1];
      Nq_diff[vv][2] = N_diff[3 * (shift + vv) + 2];
    }

    ref_hex.get_volume_coords(Nq, quad_coords);
    // ****************************************
    double *ksi = quad_coords;
    double diff_ksi[3][3];

    ref_hex.get_volume_diff_coords(Nq_diff, diff_ksi);

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
