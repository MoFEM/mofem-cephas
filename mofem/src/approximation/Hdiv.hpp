/** \file Hdiv.hpp

  \brief Implementation of H-curl base function.

  Based on Hierarchic Finite Element Bases on Unstructured Tetrahedral
  Meshes, by Mark Ainsworth and Joe Coyle
  Shape functions for MBTRI/MBTET and HCurl space

*/

#ifndef __HDIV_HPP__
#define __HDIV_HPP__

namespace MoFEM {

  /**
   * \brief Hdiv base functions, Edge-based face functions
   * @param  faces_nodes      Face nodes on tetrahedral
   * @param  p                Approximation order on faces
   * @param  N                Shape functions
   * @param  diffN            Derivatives of shape functions
   * @param  phi_f_e          Base functions (returned)
   * @param  diff_phi_f_e     Derivatives of base functions (returned)
   * @param  gdim             Number of integration pts
   * @param  base_polynomials base function (Legendre/Lobbatto-Gauss)
   * @return                  error code
   */
  PetscErrorCode Hdiv_EdgeFaceShapeFunctions_MBTET(
    int *faces_nodes,int *p,double *N,double *diffN,double *phi_f_e[4][3],double *diff_phi_f_e[4][3],int gdim,
    PetscErrorCode (*base_polynomials)(int p,double s,double *diff_s,double *L,double *diffL,const int dim)
  );

  /**
   * \brief Hdiv base functions, Edge-based face functions
   * @param  faces_nodes      Face nodes on face
   * @param  p                Approximation order on faces
   * @param  N                Shape functions
   * @param  diffN            Derivatives of shape functions
   * @param  phi_f_e          Base functions (returned)
   * @param  diff_phi_f_e     Derivatives of base functions (returned)
   * @param  gdim             Number of integration pts
   * @param  base_polynomials base function (Legendre/Lobbatto-Gauss)
   * @param  nb               Number of nodes on entity (4 if tet, 3 if triangle)
   * @return                  error code
   */
  PetscErrorCode Hdiv_EdgeFaceShapeFunctions_MBTET_ON_FACE(
    int *faces_nodes,int p,double *N,double *diffN,double *phi_f_e[3],double *diff_phi_f_e[3],int gdim,int nb,
    PetscErrorCode (*base_polynomials)(int p,double s,double *diff_s,double *L,double *diffL,const int dim)
  );

  /**
   * \brief Face bubble functions
   * @param  faces_nodes      Face nodes on tetrahedral
   * @param  p                Approx. order on faces
   * @param  N                Shape function
   * @param  diffN            Derivatives of shape functions
   * @param  phi_f            Base functions
   * @param  diff_phi_f       Derivatives of base functions
   * @param  gdim             Number of integration pts
   * @param  base_polynomials Base function (Legendre/Lobbatto-Gauss)
   * @return                  error code
   */
  PetscErrorCode Hdiv_FaceBubbleShapeFunctions_MBTET(
    int *faces_nodes,int *p,double *N,double *diffN,double *phi_f[],double *diff_phi_f[],int gdim,
    PetscErrorCode (*base_polynomials)(int p,double s,double *diff_s,double *L,double *diffL,const int dim)
  );

  /**
   * \brief Face bubble functions
   * @param  faces_nodes      Face nodes on face
   * @param  p                Approx. order on face
   * @param  N                Shape function
   * @param  diffN            Derivatives of shape functions
   * @param  phi_f            Base functions
   * @param  diff_phi_f       Derivatives of base functions
   * @param  gdim             Number of integration pts
   * @param  nb               Number of nodes on entity (4 if tet, 3 if triangle)
   * @param  base_polynomials Base function (Legendre/Lobbatto-Gauss)
   * @return                  error code
   */
  PetscErrorCode Hdiv_FaceBubbleShapeFunctions_MBTET_ON_FACE(
    int *faces_nodes,int p,double *N,double *diffN,double *phi_f,double *diff_phi_f,int gdim,int nb,
    PetscErrorCode (*base_polynomials)(int p,double s,double *diff_s,double *L,double *diffL,const int dim)
  );

  /**
   * \brief Hdiv base function, Edge-based interior (volume) functions
   * @param  p                volume approx. order
   * @param  N                shape functions
   * @param  diffN            derivatives of shape functions
   * @param  phi_v_e          base functions (return)
   * @param  diff_phi_v_e     derivatives of base functions (returned)
   * @param  gdim             number of integration points
   * @param  base_polynomials base function (Legendre/Lobbatto-Gauss)
   * @return                  error code
   */
  PetscErrorCode Hdiv_EdgeBasedVolumeShapeFunctions_MBTET(
    int p,
    double *N,
    double *diffN,
    double *phi_v_e[6],
    double *diff_phi_v_e[6],
    int gdim,
    PetscErrorCode (*base_polynomials)(
      int p,double s,double *diff_s,double *L,double *diffL,const int dim
    )
  );


}

#endif // __HDIV_HPP__
