/** \file Hcurl.gpp

  Based on Hierarchic Finite Element Bases on Unstructured Tetrahedral
  Meshes, by Mark Ainsworth and Joe Coyle
  Shape functions for MBTRI/MBTET and HCurl space

*/

#ifndef __HCURL_HPP__
#define __HCURL_HPP__

namespace MoFEM {

/**
 * \brief Edge based H-curl base functions
 * @param  sense            sense fo edge (i.e. unique orientation)
 * @param  p                array of oder for each edge
 * @param  N                array shape functions evaluated at each integration point
 * @param  diffN            derivatives of shape functions
 * @param  edgeN            base functions on edges
 * @param  diff_edgeN       derivatives of edge shape functions
 * @param  GDIM             number of integration points
 * @param  base_polynomials polynomial base function (f.e. Legendre of Lobatto)
 * @return                  error code
 */
PetscErrorCode Hcurl_EdgeBaseFunctions_MBTET(
  int *sense,int *p,double *N,double *diffN,double *edgeN[],double *diff_edgeN[],int GDIM,
  PetscErrorCode (*base_polynomials)(int p,double s,double *diff_s,double *L,double *diffL,const int dim)
);

/** \brief Face edge base functions of Hcurl space.

  On each edge we have (P-1) base functions, and each face has 3 edges and are 4
  faces on tets.

  See NBEDGE_HCURL

  * @param  face_nodes       array [4*3] of local indices of face nodes
  * @param  p                approximation order
  * @param  N                array shape functions evaluated at each integration point
  * @param  diffN            derivatives of nodal shape functions
  * @param  phi_f[4]         calculated shape functions for each face
  * @param  diff_phi_v[4]    derivatives of shape functions for each face
  * @param  GDIM             number of shape functions
  * @param  base_polynomials polynomial base function (f.e. Legendre of Lobatto)
  * @return                  error code

*/
PetscErrorCode Hcurl_EdgeBasedFaceFunctions_MBTET(
  int *faces_nodes,int *p,double *N,double *diffN,double *phi_f_e[4][3],double *diff_phi_f_e[4][3],int GDIM,
  PetscErrorCode (*base_polynomials)(int p,double s,double *diff_s,double *L,double *diffL,const int dim)
);

/** \brief Face edge base functions of Hcurl space.

  On each face we have P*(P-1) base functions and are 4 faces.

  See NBFACETRI_EDGE_HCURL

  * @param  face_nodes       array [4*3] of local indices of face nodes
  * @param  p                approximation order
  * @param  N                array shape functions evaluated at each integration point
  * @param  diffN            derivatives of nodal shape functions
  * @param  phi_f[4]         calculated shape functions for each face
  * @param  diff_phi_v[4]    derivatives of shape functions for each face
  * @param  GDIM             number of shape functions
  * @param  base_polynomials polynomial base function (f.e. Legendre of Lobatto)
  * @return                  error code

*/
PetscErrorCode Hcurl_BubbleFaceFunctions_MBTET(
  int *faces_nodes,int *p,double *N,double *diffN,double *phi_f[4],double *diff_phi_f[4],int GDIM,
  PetscErrorCode (*base_polynomials)(int p,double s,double *diff_s,double *L,double *diffL,const int dim)
);

/** \brief Face base interior function

On each face we have P*(P-1)/2 and are 4 faces on Tetrahedral.

See NBVOLUMETET_FACE_HCURL

* @param  face_nodes       array [4*3] of local indices of face nodes
* @param  p                approximation order
* @param  N                nodal shape functions
* @param  diffN            derivatives of nodal shape functions
* @param  phi_v            calculated shape functions
* @param  diff_phi_v       derivatives of shape functions
* @param  GDIM             number of shape functions
* @param  base_polynomials polynomial base function (f.e. Legendre of Lobatto)
* @return                  error code


*/
PetscErrorCode Hcurl_FaceInteriorFunctions_MBTET(
  int *faces_nodes,int p,double *N,double *diffN,double *phi_v,double *diff_phi_v,int GDIM,
  PetscErrorCode (*base_polynomials)(int p,double s,double *diff_s,double *L,double *diffL,const int dim)
);

/** \brief Volume interior function

On volume have (P-3)*(P-2)*(P-1)/2.

See NBVOLUMETET_TET_HCURL

* @param  p                approximation order
* @param  N                nodal shape functions
* @param  diffN            derivatives of nodal shape functions
* @param  phi_v            calculated shape functions
* @param  diff_phi_v       derivatives of shape functions
* @param  GDIM             number of shape functions
* @param  base_polynomials polynomial base function (f.e. Legendre of Lobatto)
* @return                  error code

*/
PetscErrorCode Hcurl_VolumeInteriorFunctions_MBTET(
  int p,double *N,double *diffN,double *phi_v,double *diff_phi_v,int GDIM,
  PetscErrorCode (*base_polynomials)(int p,double s,double *diff_s,double *L,double *diffL,const int dim)
);

/** \brief Face H-curl functions

  Face H-curl functions are set of Eddge-Based Face functions
  (Hcurl_EdgeBasedFaceFunctions_MBTET) and Bubble-Face functions (Hcurl_BubbleFaceFunctions_MBTET).

  See NBVOLUMETET_FACE_HCURL

  * @param  face_nodes       array [4*3] of local indices of face nodes
  * @param  p                approximation order
  * @param  N                nodal shape functions
  * @param  diffN            derivatives of nodal shape functions
  * @param  phi_f[4]         calculated shape functions for each face
  * @param  diff_phi_v[4]    derivatives of shape functions for each face
  * @param  GDIM             number of shape functions
  * @param  base_polynomials polynomial base function (f.e. Legendre of Lobatto)
  * @return                  error code

*/
PetscErrorCode Hcurl_FaceFunctions_MBTET(
  int *faces_nodes,int *p,double *N,double *diffN,double *phi_f[4],double *diff_phi_f[4],int GDIM,
  PetscErrorCode (*base_polynomials)(int p,double s,double *diff_s,double *L,double *diffL,const int dim)
);

/**
 \brief H-curl volume base functions

 Volume base functions are collection of iace interior base functions and
 volume interior base functions.

 * @param  p                approximation order
 * @param  N                nodal shape functions
 * @param  diffN            derivatives of nodal shape functions
 * @param  phi_v            calculated shape functions
 * @param  diff_phi_v       derivatives of shape functions
 * @param  GDIM             number of shape functions
 * @param  base_polynomials polynomial base function (f.e. Legendre of Lobatto)
 * @return                  error code

 */
PetscErrorCode Hcurl_VolumeFunctions_MBTET(
  int p,double *N,double *diffN,double *phi_v,double *diff_phi_v,int GDIM,
  PetscErrorCode (*base_polynomials)(int p,double s,double *diff_s,double *L,double *diffL,const int dim)
);


}

#endif // __HCURL_HPP__
