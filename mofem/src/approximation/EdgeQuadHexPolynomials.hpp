/** \file EdgeQuadHexPolynomials.hpp

  \brief Implementation of base Demkowicz functions on quad and hex

  Based on Hierarchic Finite Element Bases on Unstructured QUADS (2D)
  and HEXES (3D) Meshes

*/

#ifndef __EDGE_QUAD_HEX_HPP__
#define __EDGE_QUAD_HEX_HPP__

namespace MoFEM {
namespace DemkowiczHexAndQuad {

MoFEMErrorCode Legendre_polynomials01(int p, double s, double *L);

MoFEMErrorCode Integrated_Legendre01(int p, double s, double *L, double *diffL);

MoFEMErrorCode Face_orientMat(int *face_nodes, double orientMat[2][2]);

MoFEMErrorCode MonomOrdering(int perm[][3], int p, int q, int r = 0);

/*
Segment (1d) basis functions

       0-----------1
*/
MoFEMErrorCode H1_BubbleShapeFunctions_ONSEGMENT(int p, double *N,
                                                 double *diffN, double *bubbleN,
                                                 double *diff_bubbleN,
                                                 int nb_integration_pts);

MoFEMErrorCode L2_ShapeFunctions_ONSEGMENT(int p, double *N, double *diffN,
                                           double *funN, double *funDiffN,
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
                                            double *diffN, double *edgeN[4],
                                            double *edgeDiffN[4],
                                            int nb_integration_pts);

/**
* \brief H1 Face bubble functions on Quad.

Function generates hierarchical base of H1 comforting functions on a 2D quad.

* @param  face_nodes       face nodes order
* @param  p                array of orders (in each direction) of base functions
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
MoFEMErrorCode H1_FaceShapeFunctions_ONQUAD(int *face_nodes, int *p, double *N,
                                            double *diffN, double *faceN,
                                            double *diff_faceN,
                                            int nb_integration_pts);

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
MoFEMErrorCode L2_FaceShapeFunctions_ONQUAD(int *p, double *N, double *diffN,
                                            double *face_buble,
                                            double *diff_face_bubble,
                                            int nb_integration_pts);

MoFEMErrorCode Hcurl_EdgeShapeFunctions_ONQUAD(int *sense, int *p, double *N,
                                               double *diffN, double *edgeN[4],
                                               double *curl_edgeN[4],
                                               int nb_integration_pts);

MoFEMErrorCode Hcurl_FaceShapeFunctions_ONQUAD(int *face_nodes, int *p,
                                               double *N, double *diffN,
                                               double *faceN[],
                                               double *diff_faceN[],
                                               int nb_integration_pts);

MoFEMErrorCode Hdiv_FaceShapeFunctions_ONQUAD(int *p, double *N,
                                              double *faceN[],
                                              double *div_faceN[],
                                              int nb_integration_pts);

/* Reference Hex and its canonical vertex and edge numbering
               7---------10----------6
              /|                    /|
             / |                   / |
           11  |                  9  |             x3
           /   7                 /   |             |
          /    |                /    6             |    x2
         4----------8----------5     |             |   /
         |     |               |     |             |  /
         |     3----------2----------2            | /
         |    /                |    /              o -------- x1
         4   /                 5   /
         |  3                  |  1
         | /                   | /
         |/                    |/
         0---------0-----------1

  Hex Face Canonical numbering

        1. 1 2 3 4
        2. 1 2 6 5
        3. 2 3 7 6
        4. 4 3 7 8
        5. 1 4 8 5
        6. 5 6 7 8
*/

MoFEMErrorCode H1_EdgeShapeFunctions_ONHEX(int *sense, int *p, double *N,
                                           double *N_diff, double *edgeN[12],
                                           double *diff_edgeN[12],
                                           int nb_integration_pts);

MoFEMErrorCode H1_FaceShapeFunctions_ONHEX(int *face_nodes[6], int *p,
                                           double *N, double *N_diff,
                                           double *faceN[6],
                                           double *diff_faceN[6],
                                           int nb_integration_pts);

MoFEMErrorCode H1_InteriorShapeFunctions_ONHEX(int *p, double *N,
                                               double *N_diff, double *faceN,
                                               double *diff_faceN,
                                               int nb_integration_pts);

MoFEMErrorCode L2_InteriorShapeFunctions_ONHEX(const int *p, double *N,
                                               double *N_diff, double *volN,
                                               double *diff_volN,
                                               int nb_integration_pts);

MoFEMErrorCode Hcurl_EdgeShapeFunctions_ONHEX(int *sense, int *p, double *N,
                                              double *N_diff, double *edgeN[12],
                                              double *diff_edgeN[12],
                                              int nb_integration_pts);

MoFEMErrorCode Hcurl_FaceShapeFunctions_ONHEX(int *face_nodes[6], int *p,
                                              double *N, double *N_diff,
                                              double *faceN[6][2],
                                              double *diff_faceN[6][2],
                                              int nb_integration_pts);

MoFEMErrorCode Hcurl_InteriorShapeFunctions_ONHEX(int *p, double *N,
                                                  double *N_diff,
                                                  double *volN[3],
                                                  double *diff_volN[3],
                                                  int nb_integration_pts);

MoFEMErrorCode Hdiv_FaceShapeFunctions_ONHEX(int *sense[6], int *p, double *N,
                                             double *N_diff, double *faceN[6],
                                             double *div_faceN[6],
                                             int nb_integration_pts);

MoFEMErrorCode Hdiv_InteriorShapeFunctions_ONHEX(int *p, double *N,
                                                 double *N_diff,
                                                 double *bubleN[3],
                                                 double *div_bubleN[3],
                                                 int nb_integration_pts);

} // namespace DemkowiczHexAndQuad

} // namespace MoFEM

#endif // __EDGE_QUAD_HEX_HPP__