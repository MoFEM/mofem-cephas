/** \file EdgeQuadHexPolynomials.hpp

  \brief Implementation of base functions on 1D segment (Edge), .

  Based on Hierarchic Finite Element Bases on Unstructured Tetrahedral
  Meshes, by \cite NME:NME847 and \cite fuentes2015orientation
  Shape functions for MBTRI/MBTET and HCurl space

*/

#ifndef __EDGE_QUAD_HEX_HPP__
#define __EDGE_QUAD_HEX_HPP__

namespace MoFEM {

/**
 * \brief H1 Edge bubble base functions on Edge

 Function generates hierarchical base of H1 comforting functions on a 1D segment
 (edge).

 * @param  sense            orientation of edge (takes 1 or -1)
 * @param  p                order of base functions being generated
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
MoFEMErrorCode H1_TP_EdgeShapeFunctions_MBEDGE(
    int sense, int p, double *N, double *diffN, double *edgeN[],
    double *diff_edgeN[], int nb_integration_pts,
    PetscErrorCode (*base_polynomials)(int p, double s, double *diff_s,
                                       double *L, double *diffL,
                                       const int dim));

/**
* \brief L2 Edge base functions on Edge

Function generates hierarchical base of L2 comforting functions on a 1D segment
(edge).

* @param  sense            orientation of edge (takes 1 or -1)
* @param  p                order of base functions
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
MoFEMErrorCode L2_TP_EdgeShapeFunctions_MBEDGE(
    int sense, int p, double *N, double *diffN, double *edgeN[],
    double *diff_edgeN[], int nb_integration_pts,
    PetscErrorCode (*base_polynomials)(int p, double s, double *diff_s,
                                       double *L, double *diffL,
                                       const int dim));

/**
* \brief H1 Edge base functions on Quad

Function generates hierarchical base of H1 comforting functions on a 2D Quad.

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
MoFEMErrorCode H1_TP_EdgeShapeFunctions_MBQUAD(
    int *sense, int *p, double *N, double *diffN, double *edgeN[],
    double *diff_edgeN[], int nb_integration_pts,
    PetscErrorCode (*base_polynomials)(int p, double s, double *diff_s,
                                       double *L, double *diffL,
                                       const int dim));

/**
* \brief H1 Face bubble functions on Quad

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
MoFEMErrorCode H1_TP_QuadShapeFunctions_MBQUAD(
    int *sense, int *p, double *N, double *diffN, double *edgeN[],
    double *diff_edgeN[], int nb_integration_pts,
    PetscErrorCode (*base_polynomials)(int p, double s, double *diff_s,
                                       double *L, double *diffL,
                                       const int dim));

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
MoFEMErrorCode L2_TP_QuadShapeFunctions_MBQUAD(
    int *p, double *N, double *diffN, double *edgeN[], double *diff_edgeN[],
    int nb_integration_pts,
    PetscErrorCode (*base_polynomials)(int p, double s, double *diff_s,
                                       double *L, double *diffL,
                                       const int dim));

} // namespace MoFEM

#endif // __EDGE_QUAD_HEX_HPP__