/** \file EdgeQuadHexPolynomials.hpp

  \brief Implementation of base functions on 1D segment (Edge), .

  Based on Hierarchic Finite Element Bases on Unstructured QUADS (2D) 
  and HEXES (3D) Meshes

*/

#ifndef __EDGE_QUAD_HEX_HPP__
#define __EDGE_QUAD_HEX_HPP__

namespace MoFEM {



MoFEMErrorCode Legendre_polynomials01(int p, double s, double *L);

MoFEMErrorCode Integrated_Legendre01(int p, double s, double *L,
                                          double *diffL);

MoFEMErrorCode Face_orientMat(int *face_nodes, double orientMat[2][2]);

/*
Segment (1d) basis functions

       0-----------1
*/
MoFEMErrorCode H1_BubbleShapeFunctions_ONSEGMENT(int p, double *L,
                                                 double *bubbleN,
                                                 double *diff_bubbleN,
                                                int nb_integration_pts);

MoFEMErrorCode L2_ShapeFunctions_ONSEGMENT(int p, double *L, double *funN,
                                           int nb_integration_pts);

/*
      Quads
 3-------2------2
 |              |       y
 |              |       ^
 3              1       |
 |              |       |
 |              |       0-----  > x
 0-------0------1
*/

/**
* \brief H1 Edge base functions on Quad

Function generates hierarchical base of H1 comforting functions on a 2D
Quad.

* @param  sense            array of orientation of edges (take 1 or -1)
* @param  p               array of orders (in each direction) of base
functions
* @param  N                array vertex shape functions evaluated at each
integration point
* @param  diffN            derivatives of vertex shape functions
* @return  edgeN           base functions on edges at qd pts
* @return  diff_edgeN      derivatives of edge shape functions at qd pts
* @param  nb_integration_pts           number of integration points
* @param  base_polynomials polynomial base function (f.e. Legendre of
Integrated Legendre)
* @return                  error code
*/
MoFEMErrorCode H1_EdgeShapeFunctions_ONQUAD(int *sense, int *p, double *N,
                                            double *edgeN[4],
                                            double *diff_edgeN[4],
                                            int nb_integration_pts);

/**
* \brief H1 Face bubble functions on Quad.

Function generates hierarchical base of H1 comforting functions on a 2D quad.

* @param  sense            array of orientation of edges (take 1 or -1)
* @param  p               array of orders (in each direction) of base functions
* @param  N                array vertex shape functions evaluated at each
integration point
* @param  diffN            derivatives of vertex shape functions
* @return  edgeN           base functions on edges at qd pts
* @return  diff_edgeN      derivatives of edge shape functions at qd pts
* @param  nb_integration_pts           number of integration points
* @param  base_polynomials polynomial base function (f.e. Legendre of Integrated
Legendre)
* @return                  error code
*/
MoFEMErrorCode H1_FaceShapeFunctions_ONQUAD(int *p, double *N, double *faceN,
                                            double *diff_faceN, int nb_integration_pts);

/**
* \brief L2 Face base functions on Quad

Function generates hierarchical base of H1 comforting functions on a 2D quad.

* @param  p               array of orders (in each direction) of base functions
* @param  N                array vertex shape functions evaluated at each
integration point
* @param  diffN            derivatives of vertex shape functions
* @return  edgeN           base functions on edges at qd pts
* @return  diff_edgeN      derivatives of edge shape functions at qd pts
* @param  nb_integration_pts           number of integration points
* @param  base_polynomials polynomial base function (f.e. Legendre of Integrated
Legendre)
* @return                  error code
*/
MoFEMErrorCode L2_FaceShapeFunctions_ONQUAD(int *p, double *N, 
                                            double *face_buble, double *diff_face_bubble,
                                            int nb_integration_pts);

MoFEMErrorCode Hcurl_EdgeShapeFunctions_ONQUAD(int *sense, int *p, double *N,
                                               double *edgeN[4],
                                               double *curl_edgeN[4],
                                               int nb_integration_pts);

MoFEMErrorCode Hcurl_FaceShapeFunctions_ONQUAD(int *p, double *N,
                                               double *faceN[2],
                                               double *curl_faceN[2],
                                               int nb_integration_pts);

MoFEMErrorCode Hdiv_EdgeShapeFunctions_ONQUAD(int *sense, int *p, double *N,
                                               double *edgeN[],
                                               double *div_edgeN[],
                                               int nb_integration_pts);

MoFEMErrorCode Hdiv_FaceShapeFunctions_ONQUAD(int *p, double *N,
                                               double *faceN[],
                                               double *div_faceN[],
                                               int nb_integration_pts);

/* Reference Hex and its canonical vertex and edge numbering
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

  Hex Face Canonical numbering

        1. 1 2 3 4
        2. 1 2 6 5
        3. 2 3 7 6
        4. 4 3 7 8
        5. 1 4 8 5
        6. 5 6 7 8
*/

MoFEMErrorCode H1_EdgeShapeFunctions_ONHEX(int        *sense, 
                                           int         *p,   
                                           double      *N,
                                           double      *edgeN[12],
                                           double      *diff_edgeN[12],
                                           int         nb_integration_pts);

MoFEMErrorCode H1_FaceShapeFunctions_ONHEX(int        *face_nodes[6],
                                           int        *p, 
                                           double     *N, 
                                           double     *faceN[6],
                                           double     *diff_faceN[6],
                                           int        nb_integration_pts);

MoFEMErrorCode H1_InteriorShapeFunctions_ONHEX(int       *p, 
                                               double    *N,
                                               double    *faceN,
                                               double    *diff_faceN,
                                               int       nb_integration_pts);

MoFEMErrorCode L2_InteriorShapeFunctions_ONHEX(int        *p, 
                                           double     *N,
                                           double     *volN,
                                           int        nb_integration_pts);

MoFEMErrorCode Hcurl_EdgeShapeFunctions_ONHEX(int       *sense, 
                                              int       *p, 
                                              double    *N,
                                              double    *edgeN[12],
                                              double    *curl_edgeN[12],
                                              int       nb_integration_pts);

MoFEMErrorCode Hcurl_FaceShapeFunctions_ONHEX(int *face_nodes[6],
                                              int       *p, 
                                              double    *N,
                                              double    *faceN[6][2],
                                              double    *curl_faceN[6][2],
                                              int       nb_integration_pts);

MoFEMErrorCode Hcurl_InteriorShapeFunctions_ONHEX(int *p, 
                                              double *N,
                                              double *volN[3],
                                              double *curl_volN[3],
                                              int nb_integration_pts);

MoFEMErrorCode Hdiv_FaceShapeFunctions_ONHEX(int       *sense[6], 
                                             int       *p, 
                                             double    *N,
                                             double    *faceN[6],
                                             double    *div_faceN[6],
                                             int       nb_integration_pts);

MoFEMErrorCode Hdiv_InteriorShapeFunctions_ONHEX(int      *p, 
                                                 double   *N,
                                                 double   *bubleN[3],
                                                 double   *div_bubleN[3],
                                                 int      nb_integration_pts);

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
  void get_edge_coords(int *sense, double (&volume_coords)[3], double (&edge_coords)[12]) {

    int free_edge_coords[12] = {0, 1, 0, 1, 2, 2, 2, 2, 0, 1, 0, 1};

    for (int ee = 0; ee < 12; ee++) {
      int cc = free_edge_coords[ee];
      edge_coords[ee] =
          (double)sense[ee] * volume_coords[cc] + 0.5 * (1.0 - (double)sense[ee]);
    }
  }
  void get_edge_diff_coords(int *sense, double (&edge_diff_coords)[12][3]) {

    double diff_coords[3][3] = {
        {1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {0.0, 0.0, 1.0}};
    int free_edge_coords[12] = {0, 1, 0, 1, 2, 2, 2, 2, 0, 1, 0, 1};
    for (int ee = 0; ee < 12; ee++) {
      for (int cc = 0; cc < 3; cc++) {
        int n = free_edge_coords[ee];
        edge_diff_coords[ee][cc] = (double) sense[ee] * diff_coords[n][cc];
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

} // namespace MoFEM

#endif // __EDGE_QUAD_HEX_HPP__