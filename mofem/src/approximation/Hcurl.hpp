/** \file Hcurl.gpp

  Based on Hierarchic Finite Element Bases on Unstructured Tetrahedral
  Meshes, by Mark Ainsworth and Joe Coyle
  Shape functions for MBTRI/MBTET and HCurl space

*/

#ifndef __HCURL_HPP__
#define __HCURL_HPP__

namespace MoFEM {

PetscErrorCode Hcurl_EdgeBaseFunctions_MBTET(
  int *sense,int *p,double *N,double *diffN,double *edgeN[],double *diff_edgeN[],int GDIM,
  PetscErrorCode (*base_polynomials)(int p,double s,double *diff_s,double *L,double *diffL,const int dim)
);

/** \brief Face edge base functions of Hcurl space.

  On each edge we have (P-1) base functions, and each face has 3 edges and are 4
  faces on tets.

  See NBEDGE_HCURL

*/
PetscErrorCode Hcurl_EdgeBasedFaceFunctions_MBTET(
  int *faces_nodes,int *p,double *N,double *diffN,double *phi_f_e[4][3],double *diff_phi_f_e[4][3],int GDIM,
  PetscErrorCode (*base_polynomials)(int p,double s,double *diff_s,double *L,double *diffL,const int dim)
);

/** \brief Face edge base functions of Hcurl space.

  On each face we have P*(P-1) base functions and are 4 faces.

  See NBFACETRI_EDGE_HCURL

*/
PetscErrorCode Hcurl_BubbleFaceFunctions_MBTET(
  int *faces_nodes,int *p,double *N,double *diffN,double *phi_f[4],double *diff_phi_f[4],int GDIM,
  PetscErrorCode (*base_polynomials)(int p,double s,double *diff_s,double *L,double *diffL,const int dim)
);

/** \brief Face base interior function

On each face we have P*(P-1)/2 and are 4 faces on Tetrahedral.

See NBVOLUMETET_FACE_HCURL

*/
PetscErrorCode Hcurl_FaceInteriorFunctions_MBTET(
  int *faces_nodes,int p,double *N,double *diffN,double *phi_v,double *diff_phi_v,int GDIM,
  PetscErrorCode (*base_polynomials)(int p,double s,double *diff_s,double *L,double *diffL,const int dim)
);

/** \brief Volume interior function

On volume have (P-3)*(P-2)*(P-1)/2.

See NBVOLUMETET_TET_HCURL

*/
PetscErrorCode Hcurl_VolumeInteriorFunctions_MBTET(
  int p,double *N,double *diffN,double *phi_v,double *diff_phi_v,int GDIM,
  PetscErrorCode (*base_polynomials)(int p,double s,double *diff_s,double *L,double *diffL,const int dim)
);

}

#endif // __HCURL_HPP__
