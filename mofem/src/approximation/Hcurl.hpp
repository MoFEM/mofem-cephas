/** \file Hcurl.hpp

  \brief Implementation of H-curl base function.

  Based on Hierarchic Finite Element Bases on Unstructured Tetrahedral
  Meshes, by \cite NME:NME847 and \cite fuentes2015orientation
  Shape functions for MBTRI/MBTET and HCurl space

*/

#ifndef __HCURL_HPP__
#define __HCURL_HPP__

namespace MoFEM {

/**
 * \brief Edge based H-curl base functions on tetrahedral

 Function generates hierarchical base of h-curl comforting functions on
 tetrahedral edge.  For more details see \cite ainsworth2011bernstein.

 On each tetrahedral's edge we have P+1 functions. See NBEDGE_AINSWORTH_HCURL

 * @param  sense            sense fo edge (i.e. unique orientation)
 * @param  p                array of oder for each edge
 * @param  N                array shape functions evaluated at each integration
 point
 * @param  diffN            derivatives of shape functions
 * @param  edgeN            base functions on edges
 * @param  diff_edgeN       derivatives of edge shape functions
 * @param  nb_integration_pts             number of integration points
 * @param  base_polynomials polynomial base function (f.e. Legendre of Lobatto)
 * @return                  error code
 */
MoFEMErrorCode Hcurl_Ainsworth_EdgeBaseFunctions_MBTET(
    int *sense, int *p, double *N, double *diffN, double *edgeN[],
    double *diff_edgeN[], int nb_integration_pts,
    PetscErrorCode (*base_polynomials)(int p, double s, double *diff_s,
                                       double *L, double *diffL,
                                       const int dim));

/**
 * \brief Edge based H-curl base functions on edge

 Function generates hierarchical base of h-curl comforting functions on
 tetrahedral edge.  For more details see \cite ainsworth2011bernstein.

 On each edge we have P+1 functions. See NBEDGE_AINSWORTH_HCURL

 * @param  sense            sense fo edge (i.e. unique orientation)
 * @param  p                array of oder for each edge
 * @param  N                array shape functions evaluated at each integration
 point
 * @param  diffN            derivatives of shape functions
 * @param  edgeN            base functions on edges
 * @param  diff_edgeN       derivatives of edge shape functions
 * @param  nb_integration_pts             number of integration points
 * @param  base_polynomials polynomial base function (f.e. Legendre of Lobatto)
 * @return                  error code
 */
MoFEMErrorCode Hcurl_Ainsworth_EdgeBaseFunctions_MBTET_ON_EDGE(
    int sense, int p, double *N, double *diffN, double *edgeN,
    double *diff_edgeN, int nb_integration_pts,
    PetscErrorCode (*base_polynomials)(int p, double s, double *diff_s,
                                       double *L, double *diffL,
                                       const int dim));

/**
 * \brief Edge based H-curl base functions on face

 Function generates hierarchical base of h-curl comforting functions on
 tetrahedral edge.  For more details see \cite ainsworth2011bernstein.

 On each face's edge we have P+1 functions. See NBEDGE_AINSWORTH_HCURL

 * @param  sense            sense fo edge (i.e. unique orientation)
 * @param  p                array of oder for each edge
 * @param  N                array shape functions evaluated at each integration
 point
 * @param  diffN            derivatives of shape functions
 * @param  edgeN            base functions on edges
 * @param  diff_edgeN       derivatives of edge shape functions
 * @param  nb_integration_pts             number of integration points
 * @param  base_polynomials polynomial base function (f.e. Legendre of Lobatto)
 * @return                  error code
 */
MoFEMErrorCode Hcurl_Ainsworth_EdgeBaseFunctions_MBTET_ON_FACE(
    int *sense, int *p, double *N, double *diffN, double *edgeN[],
    double *diff_edgeN[], int nb_integration_pts,
    PetscErrorCode (*base_polynomials)(int p, double s, double *diff_s,
                                       double *L, double *diffL,
                                       const int dim));

/** \brief Face edge base functions of Hcurl space on tetrahedral.

  On each edge we have (P-1) base functions, and each face has 3 edges and are 4
  faces on tetrahedral.

  See NBFACETRI_AINSWORTH_EDGE_HCURL

  * @param  face_nodes       array [4*3] of local indices of face nodes
  * @param  p                approximation order
  * @param  N                array shape functions evaluated at each integration
  point
  * @param  diffN            derivatives of nodal shape functions
  * @param  phi_f[4]         calculated shape functions for each face
  * @param  diff_phi_v[4]    derivatives of shape functions for each face
  * @param  nb_integration_pts             number of shape functions
  * @param  base_polynomials polynomial base function (f.e. Legendre of Lobatto)
  * @return                  error code

*/
MoFEMErrorCode Hcurl_Ainsworth_EdgeBasedFaceFunctions_MBTET(
    int *faces_nodes, int *p, double *N, double *diffN, double *phi_f_e[4][3],
    double *diff_phi_f_e[4][3], int nb_integration_pts,
    PetscErrorCode (*base_polynomials)(int p, double s, double *diff_s,
                                       double *L, double *diffL,
                                       const int dim));

/** \brief Face edge base functions of Hcurl space.

  On each edge we have (P-1) base functions, and each face has 3 edges and are 4
  faces on tetrahedral.

  See NBFACETRI_AINSWORTH_EDGE_HCURL

  * @param  face_nodes       array [4*3] of local indices of face nodes
  * @param  p                approximation order
  * @param  N                array shape functions evaluated at each integration
  point
  * @param  diffN            derivatives of nodal shape functions
  * @param  phi_f[4]         calculated shape functions for each face
  * @param  diff_phi_v[4]    derivatives of shape functions for each face
  * @param  nb_integration_pts             number of shape functions
  * @param  base_polynomials polynomial base function (f.e. Legendre of Lobatto)
  * @return                  error code

*/
MoFEMErrorCode Hcurl_Ainsworth_EdgeBasedFaceFunctions_MBTET_ON_FACE(
    int *faces_nodes, int p, double *N, double *diffN, double *phi_f_e[3],
    double *diff_phi_f_e[3], int nb_integration_pts,
    PetscErrorCode (*base_polynomials)(int p, double s, double *diff_s,
                                       double *L, double *diffL,
                                       const int dim));

/** \brief Face edge base functions of Hcurl space on face on tetrahedral.

  On each face we have P*(P-1) base functions and are 4 faces.

  See NBFACETRI_AINSWORTH_EDGE_HCURL

  * @param  face_nodes       array [4*3] of local indices of face nodes
  * @param  p                approximation order
  * @param  N                array shape functions evaluated at each integration
  point
  * @param  diffN            derivatives of nodal shape functions
  * @param  phi_f[4]         calculated shape functions for each face
  * @param  diff_phi_v[4]    derivatives of shape functions for each face
  * @param  nb_integration_pts             number of shape functions
  * @param  base_polynomials polynomial base function (f.e. Legendre of Lobatto)
  * @return                  error code

*/
MoFEMErrorCode Hcurl_Ainsworth_BubbleFaceFunctions_MBTET(
    int *faces_nodes, int *p, double *N, double *diffN, double *phi_f[4],
    double *diff_phi_f[4], int nb_integration_pts,
    PetscErrorCode (*base_polynomials)(int p, double s, double *diff_s,
                                       double *L, double *diffL,
                                       const int dim));

/** \brief Face edge base functions of Hcurl space on face.

  On each face we have P*(P-1) base functions and are 4 faces.

  See NBFACETRI_AINSWORTH_EDGE_HCURL

  * @param  face_nodes       array [4*3] of local indices of face nodes
  * @param  p                approximation order
  * @param  N                array shape functions evaluated at each integration
  point
  * @param  diffN            derivatives of nodal shape functions
  * @param  phi_f[4]         calculated shape functions for each face
  * @param  diff_phi_v[4]    derivatives of shape functions for each face
  * @param  nb_integration_pts             number of shape functions
  * @param  base_polynomials polynomial base function (f.e. Legendre of Lobatto)
  * @return                  error code

*/
MoFEMErrorCode Hcurl_Ainsworth_BubbleFaceFunctions_MBTET_ON_FACE(
    int *faces_nodes, int p, double *N, double *diffN, double *phi_f,
    double *diff_phi_f, int nb_integration_pts,
    PetscErrorCode (*base_polynomials)(int p, double s, double *diff_s,
                                       double *L, double *diffL,
                                       const int dim));

/** \brief Face base interior function

On each face we have P*(P-1)/2 and are 4 faces on Tetrahedral.

See NBFACETRI_AINSWORTH_FACE_HCURL

* @param  face_nodes       array [4*3] of local indices of face nodes
* @param  p                approximation order
* @param  N                nodal shape functions
* @param  diffN            derivatives of nodal shape functions
* @param  phi_v            calculated shape functions
* @param  diff_phi_v       derivatives of shape functions
* @param  nb_integration_pts             number of shape functions
* @param  base_polynomials polynomial base function (f.e. Legendre of Lobatto)
* @return                  error code


*/
MoFEMErrorCode Hcurl_Ainsworth_FaceInteriorFunctions_MBTET(
    int *faces_nodes, int p, double *N, double *diffN, double *phi_v,
    double *diff_phi_v, int nb_integration_pts,
    PetscErrorCode (*base_polynomials)(int p, double s, double *diff_s,
                                       double *L, double *diffL,
                                       const int dim));

/** \brief Volume interior function

On volume have (P-3)*(P-2)*(P-1)/2.

See NBVOLUMETET_AINSWORTH_TET_HCURL

* @param  p                approximation order
* @param  N                nodal shape functions
* @param  diffN            derivatives of nodal shape functions
* @param  phi_v            calculated shape functions
* @param  diff_phi_v       derivatives of shape functions
* @param  nb_integration_pts             number of shape functions
* @param  base_polynomials polynomial base function (f.e. Legendre of Lobatto)
* @return                  error code

*/
MoFEMErrorCode Hcurl_Ainsworth_VolumeInteriorFunctions_MBTET(
    int p, double *N, double *diffN, double *phi_v, double *diff_phi_v,
    int nb_integration_pts,
    PetscErrorCode (*base_polynomials)(int p, double s, double *diff_s,
                                       double *L, double *diffL,
                                       const int dim));

/** \brief Face H-curl functions

  Face H-curl functions are set of Eddge-Based Face functions
  (Hcurl_Ainsworth_EdgeBasedFaceFunctions_MBTET) and Bubble-Face functions
  (Hcurl_Ainsworth_BubbleFaceFunctions_MBTET).

  See NBVOLUMETET_AINSWORTH_FACE_HCURL

  * @param  face_nodes       array [4*3] of local indices of face nodes
  * @param  p                approximation order
  * @param  N                nodal shape functions
  * @param  diffN            derivatives of nodal shape functions
  * @param  phi_f[4]         calculated shape functions for each face
  * @param  diff_phi_v[4]    derivatives of shape functions for each face
  * @param  nb_integration_pts             number of shape functions
  * @param  base_polynomials polynomial base function (f.e. Legendre of Lobatto)
  * @return                  error code

*/
MoFEMErrorCode Hcurl_Ainsworth_FaceFunctions_MBTET(
    int *face_nodes, int *p, double *N, double *diffN, double *phi_f[4],
    double *diff_phi_f[4], int nb_integration_pts,
    PetscErrorCode (*base_polynomials)(int p, double s, double *diff_s,
                                       double *L, double *diffL,
                                       const int dim));

/** \brief Face H-curl functions

  Face H-curl functions are set of Eddge-Based Face functions
  (Hcurl_Ainsworth_EdgeBasedFaceFunctions_MBTET) and Bubble-Face functions
  (Hcurl_Ainsworth_BubbleFaceFunctions_MBTET).

  See NBVOLUMETET_AINSWORTH_FACE_HCURL

  * @param  face_nodes       array [4*3] of local indices of face nodes
  * @param  p                approximation order
  * @param  N                nodal shape functions
  * @param  diffN            derivatives of nodal shape functions
  * @param  phi_f[4]         calculated shape functions for each face
  * @param  diff_phi_v[4]    derivatives of shape functions for each face
  * @param  nb_integration_pts             number of shape functions
  * @param  base_polynomials polynomial base function (f.e. Legendre of Lobatto)
  * @return                  error code

*/
MoFEMErrorCode Hcurl_Ainsworth_FaceFunctions_MBTET_ON_FACE(
    int *faces_nodes, int p, double *N, double *diffN, double *phi_f,
    double *diff_phi_f, int nb_integration_pts,
    PetscErrorCode (*base_polynomials)(int p, double s, double *diff_s,
                                       double *L, double *diffL,
                                       const int dim));

/**
 \brief H-curl volume base functions

 Volume base functions are collection of iace interior base functions and
 volume interior base functions.

 * @param  p                approximation order
 * @param  N                nodal shape functions
 * @param  diffN            derivatives of nodal shape functions
 * @param  phi_v            calculated shape functions
 * @param  diff_phi_v       derivatives of shape functions
 * @param  nb_integration_pts             number of shape functions
 * @param  base_polynomials polynomial base function (f.e. Legendre of Lobatto)
 * @return                  error code

 */
MoFEMErrorCode Hcurl_Ainsworth_VolumeFunctions_MBTET(
    int p, double *N, double *diffN, double *phi_v, double *diff_phi_v,
    int nb_integration_pts,
    PetscErrorCode (*base_polynomials)(int p, double s, double *diff_s,
                                       double *L, double *diffL,
                                       const int dim));

/**
 * \brief Edge based H-curl base functions on tetrahedral

 Function generates hierarchical base of h-curl comforting functions on
 tetrahedral edge.  For more details see \cite ainsworth2011bernstein.

 * @param  sense            sense fo edge (i.e. unique orientation)
 * @param  p                array of oder for each edge
 * @param  n                array shape functions evaluated at each integration
 point
 * @param  diff_n           derivatives of shape functions
 * @param  phi              base functions on edges
 * @param  diff_phi         derivatives of edge shape functions
 * @param  nb_integration_pts             number of integration points
 * @return                  error code
 */
MoFEMErrorCode Hcurl_Demkowicz_EdgeBaseFunctions_MBTET(
    int *sense, int *p, double *n, double *diff_n, double *phi[],
    double *diff_phi[], int nb_integration_pts);

/**
 * \brief Edge based H-curl base functions on teriangle

 Function generates hierarchical base of h-curl comforting functions on
 tetrahedral edge.  For more details see \cite ainsworth2011bernstein.

 * @param  sense            sense fo edge (i.e. unique orientation)
 * @param  p                array of oder for each edge
 * @param  n                array shape functions evaluated at each integration
 point
 * @param  diff_n           derivatives of shape functions
 * @param  phi              base functions on edges
 * @param  diff_phi         derivatives of edge shape functions
 * @param  nb_integration_pts             number of integration points
 * @return                  error code
 */
MoFEMErrorCode Hcurl_Demkowicz_EdgeBaseFunctions_MBTRI(
    int *sense, int *p, double *n, double *diff_n, double *phi[],
    double *diff_phi[], int nb_integration_pts);

/**
 * \brief Edge based H-curl base functions on edge

 Function generates hierarchical base of h-curl comforting functions on
 tetrahedral edge.  For more details see \cite ainsworth2011bernstein.

 * @param  sense            sense fo edge (i.e. unique orientation)
 * @param  p                array of oder for each edge
 * @param  n                array shape functions evaluated at each integration
 point
 * @param  diff_n           derivatives of shape functions
 * @param  phi              base functions on edges
 * @param  diff_phi         derivatives of edge shape functions
 * @param  nb_integration_pts             number of integration points
 * @return                  error code
 */
MoFEMErrorCode Hcurl_Demkowicz_EdgeBaseFunctions_MBEDGE(
    int sense, int p, double *n, double *diff_n, double *phi,
    double *diff_phi, int nb_integration_pts);

/** \brief Face base interior function

* @param  face_nodes       array [4*3] of local indices of face nodes
* @param  p                approximation order
* @param  n                nodal shape functions
* @param  diff_n            derivatives of nodal shape functions
* @param  phi            calculated shape functions
* @param  diff_phi       derivatives of shape functions
* @param  nb_integration_pts             number of shape functions
* @return                  error code

*/
MoFEMErrorCode Hcurl_Demkowicz_FaceBaseFunctions_MBTET(
    int *faces_nodes, int *p, double *n, double *diff_n, double *phi[],
    double *diff_phi[], int nb_integration_pts);

/** \brief Face base interior function

* @param  face_nodes       array [4*3] of local indices of face nodes
* @param  p                approximation order
* @param  n                nodal shape functions
* @param  diff_n            derivatives of nodal shape functions
* @param  phi            calculated shape functions
* @param  diff_phi       derivatives of shape functions
* @param  nb_integration_pts             number of shape functions
* @return                  error code

*/
MoFEMErrorCode Hcurl_Demkowicz_FaceBaseFunctions_MBTRI(
    int *faces_nodes, int p, double *n, double *diff_n, double *phi,
    double *diff_phi, int nb_integration_pts);

/** \brief Volume base interior function

* @param  p                approximation order
* @param  n                nodal shape functions
* @param  diff_n            derivatives of nodal shape functions
* @param  phi            calculated shape functions
* @param  diff_phi       derivatives of shape functions
* @param  nb_integration_pts             number of shape functions
* @return                  error code


*/
MoFEMErrorCode Hcurl_Demkowicz_VolumeBaseFunctions_MBTET(
    int p, double *n, double *diff_n, double *phi,
    double *diff_phi, int nb_integration_pts);


} // namespace MoFEM

#endif // __HCURL_HPP__
