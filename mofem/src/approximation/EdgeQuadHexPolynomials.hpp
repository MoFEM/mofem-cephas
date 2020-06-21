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

// MoFEMErrorCode L2_FaceShapeFunctions_ONHEX(int        *p, 
//                                            double     *N,
//                                            double     *face_buble,
//                                            double     *diff_face_bubble,
//                                            int        nb_integration_pts);

// MoFEMErrorCode Hcurl_EdgeShapeFunctions_ONHEX(int       *sense, 
//                                               int       *p, 
//                                               double    *N,
//                                               double    *edgeN[4],
//                                               double    *curl_edgeN[4],
//                                               int       nb_integration_pts);

// MoFEMErrorCode Hcurl_FaceShapeFunctions_ONHEX(int       *p, 
//                                                double    *N,
//                                                double    *faceN[2],
//                                                double    *curl_faceN[2],
//                                                int       nb_integration_pts);

// MoFEMErrorCode Hdiv_EdgeShapeFunctions_ONHEX(int       *sense, 
//                                              int       *p, 
//                                              double    *N,
//                                              double    *edgeN[],
//                                              double    *div_edgeN[],
//                                              int       nb_integration_pts);

// MoFEMErrorCode Hdiv_FaceShapeFunctions_ONHEX(int      *p, 
//                                              double   *N,
//                                              double   *faceN[],
//                                              double   *div_faceN[],
//                                              int      nb_integration_pts);

struct RefHex {
  RefHex(double *N, int nb_integration_pts)
      : vertices{{0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}, {1.0, 1.0, 0.0},
                 {0.0, 1.0, 0.0}, {0.0, 0.0, 1.0}, {1.0, 0.0, 1.0},
                 {1.0, 1.0, 1.0}, {0.0, 1.0, 1.0}},
        faces{{0, 1, 2, 3}, {0, 1, 5, 4}, {1, 2, 6, 5},
              {3, 2, 6, 7}, {0, 3, 7, 4}, {4, 5, 6, 7}},
        nbIntegrationPts(nb_integration_pts), vNshape(N) {
    intergration_pts = new double *[nb_integration_pts];
    for (int qq = 0; qq != nb_integration_pts; qq++) {
      int shift = 8 * qq;
      intergration_pts[qq] = new double[3];
      for (int vv = 0; vv < 8; vv++) {
        intergration_pts[qq][0] += N[shift + vv] * vertices[vv][0];
        intergration_pts[qq][1] += N[shift + vv] * vertices[vv][1];
        intergration_pts[qq][2] += N[shift + vv] * vertices[vv][2];
      }
    }
  }

  double **get_integrationPts() { return intergration_pts; }
  double **get_edge_affines() {
    double **edge_affine = 0;
    edge_affine = new double *[nbIntegrationPts];
    for (int qq = 0; qq != nbIntegrationPts; qq++) {
      edge_affine[qq] = new double[24];
      double ksi = intergration_pts[qq][0];
      double eta = intergration_pts[qq][1];
      double gma = intergration_pts[qq][2];

      double edgeAffines[24] = {
          1.0 - eta, 1.0 - gma, 0.0 + ksi, 1.0 - gma, 0.0 + eta, 1.0 - gma,
          1.0 - ksi, 1.0 - gma, 1.0 - ksi, 1.0 - eta, 0.0 + ksi, 1.0 - eta,
          0.0 + ksi, 0.0 + eta, 1.0 - ksi, 0.0 + eta, 1.0 - eta, 0.0 + gma,
          0.0 + ksi, 0.0 + gma, 0.0 + eta, 0.0 + gma, 1.0 - ksi, 0.0 + gma};
      for (int ee = 0; ee != 12; ee++) {
        edge_affine[qq][2 * ee + 0] = edgeAffines[2 * ee + 0];
        edge_affine[qq][2 * ee + 1] = edgeAffines[2 * ee + 1];
      }
    }
    return edge_affine;
  }
  double **get_edge_coords() {
    double **edge_coords = 0;
    edge_coords = new double *[nbIntegrationPts];
    int free_edge_coords[12] = {0, 1, 0, 1, 2, 2, 2, 2, 0, 1, 0, 1};
    for (int qq = 0; qq != nbIntegrationPts; qq++) {
      edge_coords[qq] = new double[12];
      for (int ee = 0; ee < 12; ee++) {
        int cc = free_edge_coords[ee];
        edge_coords[qq][ee] = intergration_pts[qq][cc];
      }
    }
    return edge_coords;
  }
  double **get_face_affines() {
    double **face_affine = 0;
    face_affine = new double *[nbIntegrationPts];
    for (int qq = 0; qq != nbIntegrationPts; qq++) {
      face_affine[qq] = new double[6];
      double ksi = intergration_pts[qq][0];
      double eta = intergration_pts[qq][1];
      double gma = intergration_pts[qq][2];

      double faceAffine[6] = {1.0 - gma, 1.0 - eta, 0.0 + ksi,
                              0.0 + eta, 1.0 - ksi, 0.0 + gma};
      for (int ff = 0; ff != 6; ff++)
        face_affine[qq][ff] = faceAffine[ff];
    }
    return face_affine;
  }

  double **get_face_coords(int face_nodes[6][4]) {
    double **face_coords = 0;
    int par_face_nodes[6][8] = {
        {4, 5, 6, 7, 0, 1, 2, 3}, {3, 2, 1, 0, 7, 6, 5, 4},
        {1, 0, 3, 2, 5, 4, 7, 6}, {3, 2, 1, 0, 7, 6, 5, 4},
        {1, 0, 3, 2, 5, 4, 7, 6}, {4, 5, 6, 7, 0, 1, 2, 3}};
    face_coords = new double *[nbIntegrationPts];
    int free_coords[6][2] = {{0, 1}, {0, 2}, {1, 2}, {0, 2}, {1, 2}, {0, 1}};
    for (int qq = 0; qq != nbIntegrationPts; qq++) {
      int quad_shift = 8 * qq;
      face_coords[qq] = new double[12];
      for (int ff = 0; ff != 6; ff++) {
        int v0 = free_coords[ff][0];
        int v1 = free_coords[ff][1];
        for (int fv = 0; fv != 4; fv++) {
          int n0 = face_nodes[ff][fv];
          int n1 = par_face_nodes[ff][n0];
          int index = faces[ff][fv];
          double N = vNshape[quad_shift + n0] + vNshape[quad_shift + n1];
          face_coords[qq][2 * ff + 0] += vertices[index][v0] * N;
          face_coords[qq][2 * ff + 1] += vertices[index][v1] * N;
        }
      }
    }
    return face_coords;
  }

private:
  double vertices[8][3];
  int faces[6][4];

  double **intergration_pts;
  int nbIntegrationPts;
  double *vNshape;
};

} // namespace MoFEM

#endif // __EDGE_QUAD_HEX_HPP__