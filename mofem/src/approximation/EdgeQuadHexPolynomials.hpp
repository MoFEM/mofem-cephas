/** \file EdgeQuadHexPolynomials.hpp

  \brief Implementation of base functions on 1D segment (Edge), .

  Based on Hierarchic Finite Element Bases on Unstructured QUADS (2D) 
  and HEXES (3D) Meshes

*/

#ifndef __EDGE_QUAD_HEX_HPP__
#define __EDGE_QUAD_HEX_HPP__

namespace MoFEM {

MoFEMErrorCode Integrated_Legendre(int p, double s, double *L,
                                          double *diffL);

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
  MoFEMErrorCode H1_EdgeShapeFunctions_ONQUAD(int *sense, int *p, double *N, double *edgeN[4],
                                 double *diff_edgeN[4], int nb_integration_pts);

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
                                               double *edgeN[],
                                               double *curl_edgeN[],
                                               int nb_integration_pts);

MoFEMErrorCode Hcurl_FaceShapeFunctions_ONQUAD(int *p, double *N,
                                               double *faceN[],
                                               double *curl_faceN[],
                                               int nb_integration_pts);

} // namespace MoFEM

#endif // __EDGE_QUAD_HEX_HPP__