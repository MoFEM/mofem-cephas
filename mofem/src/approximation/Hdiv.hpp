/** \file Hdiv.hpp

  \brief Implementation of H-curl base function.

  Based on Hierarchic Finite Element Bases on Unstructured Tetrahedral
  Meshes, by Mark Ainsworth and Joe Coyle
  Shape functions for MBTRI/MBTET and HCurl space

*/

#ifndef __HDIV_HPP__
#define __HDIV_HPP__

namespace MoFEM {

/** \brief Broken base Ainsworth subentries order change hooks.
 *
 * Hooks enabling change of DOFs numbers/polynomial order on subentries, when
 * base on element is constructed.
 *
 * Note that this functionality is global, that all functions are static.
 *
 **/
struct AinsworthOrderHooks {
  static boost::function<int(int)> broken_nbfacetri_edge_hdiv;
  static boost::function<int(int)> broken_nbfacetri_face_hdiv;
  static boost::function<int(int)> broken_nbvolumetet_edge_hdiv;
  static boost::function<int(int)> broken_nbvolumetet_face_hdiv;
  static boost::function<int(int)> broken_nbvolumetet_volume_hdiv;
};

/**
 * \brief Hdiv base functions, Edge-based face functions by Ainsworth \cite
 * NME:NME847
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
MoFEMErrorCode Hdiv_Ainsworth_EdgeFaceShapeFunctions_MBTET(
    int *faces_nodes, int *p, double *N, double *diffN, double *phi_f_e[4][3],
    double *diff_phi_f_e[4][3], int gdim,
    PetscErrorCode (*base_polynomials)(int p, double s, double *diff_s,
                                       double *L, double *diffL,
                                       const int dim));

/**
 * \brief Hdiv base functions, Edge-based face functions by Ainsworth \cite
 * NME:NME847
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
MoFEMErrorCode Hdiv_Ainsworth_EdgeFaceShapeFunctions_MBTET_ON_FACE(
    int *faces_nodes, int p, double *N, double *diffN, double *phi_f_e[3],
    double *diff_phi_f_e[3], int gdim, int nb,
    PetscErrorCode (*base_polynomials)(int p, double s, double *diff_s,
                                       double *L, double *diffL,
                                       const int dim));

/**
 * \brief Face bubble functions by Ainsworth \cite NME:NME847
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
MoFEMErrorCode Hdiv_Ainsworth_FaceBubbleShapeFunctions(
    int *faces_nodes, int *p, double *N, double *diffN, double *phi_f[],
    double *diff_phi_f[], int gdim,
    PetscErrorCode (*base_polynomials)(int p, double s, double *diff_s,
                                       double *L, double *diffL,
                                       const int dim));

/**
 * \brief Face bubble functions by Ainsworth \cite NME:NME847
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
MoFEMErrorCode Hdiv_Ainsworth_FaceBubbleShapeFunctions_ON_FACE(
    int *faces_nodes, int p, double *N, double *diffN, double *phi_f,
    double *diff_phi_f, int gdim, int nb,
    PetscErrorCode (*base_polynomials)(int p, double s, double *diff_s,
                                       double *L, double *diffL,
                                       const int dim));

/**
 * \brief Hdiv base function, Edge-based interior (volume) functions by
 * Ainsworth \cite NME:NME847
 * @param  p                volume approx. order
 * @param  N                shape functions
 * @param  diffN            derivatives of shape functions
 * @param  phi_v_e          base functions (return)
 * @param  diff_phi_v_e     derivatives of base functions (returned)
 * @param  gdim             number of integration points
 * @param  base_polynomials base function (Legendre/Lobbatto-Gauss)
 * @return                  error code
 */
MoFEMErrorCode Hdiv_Ainsworth_EdgeBasedVolumeShapeFunctions_MBTET(
    int p, double *N, double *diffN, double *phi_v_e[6],
    double *diff_phi_v_e[6], int gdim,
    PetscErrorCode (*base_polynomials)(int p, double s, double *diff_s,
                                       double *L, double *diffL,
                                       const int dim));

/** Hdiv Face-based interior functions by Ainsworth \cite NME:NME847
 * @param  p                Approximation order on face
 * @param  N                Shape functions on face
 * @param  diffN            Derivatives of shape functions of face
 * @param  phi_v_f          Base functions (returned)
 * @param  diff_phi_v_f     Derivatives of base functions (returned)
 * @param  gdim             Number of integration points
 * @param  base_polynomials Base function (Legendre/Lobbatto-Gauss)
 * @return                  Error code
 */
MoFEMErrorCode Hdiv_Ainsworth_FaceBasedVolumeShapeFunctions_MBTET(
    int p, double *N, double *diffN, double *phi_v_f[], double *diff_phi_v_f[],
    int gdim,
    PetscErrorCode (*base_polynomials)(int p, double s, double *diff_s,
                                       double *L, double *diffL,
                                       const int dim));

/**
 * \brief Interior bubble functions by Ainsworth \cite NME:NME847
 * @param  p                Volume order
 * @param  N                Shape functions
 * @param  diffN            Derivatives of shape functions
 * @param  phi_v            Base functions
 * @param  diff_phi_v       Derivatives of shape functions
 * @param  gdim             Number of integration points
 * @param  base_polynomials Base function (Legendre/Lobbatto-Gauss)
 * @return                  Error code
 */
MoFEMErrorCode Hdiv_Ainsworth_VolumeBubbleShapeFunctions_MBTET(
    int p, double *N, double *diffN, double *phi_v, double *diff_phi_v,
    int gdim,
    PetscErrorCode (*base_polynomials)(int p, double s, double *diff_s,
                                       double *L, double *diffL,
                                       const int dim));

/**
 * \brirf HDiv base finuctions on triangle by Demkowicz \cite
 * fuentes2015orientation
 *
 * @param  faces_nodes Nodes on face
 * @param  p           Approx. order
 * @param  N           Shape functions
 * @param  diffN       Derivatives of shape functions
 * @param  phi_f       Returned base functions on face
 * @param  diff_phi_f  Returned derivatives of base functions on face
 * @param  gdim        Number of integration points
 * @param  nb          nb is 4 for face on tetrahedral and 3 for face
 * @return             error code
 */
MoFEMErrorCode Hdiv_Demkowicz_Face_MBTET_ON_FACE(int *faces_nodes, int p,
                                                 double *N, double *diffN,
                                                 double *phi_f,
                                                 double *diff_phi_f, int gdim,
                                                 int nb);

/**
 * \brirf HDiv base finuctions in tetrahedral interior by Demkowicz \cite
 * fuentes2015orientation
 *
 * @param  p           Approximation order
 * @param  N           Shape functions
 * @param  diffN       Derivatives of base functions
 * @param  p_face      Max order on faces
 * @param  phi_f       Precalculated face base functions
 * @param  diff_phi_f  Precalculated derivatives of face base functions
 * @param  phi_v       Returned base functions in volume
 * @param  diff_phi_v  Returned derivatives of base functions in volume
 * @param  gdim        Number of integration points
 * @return             error code
 */
MoFEMErrorCode Hdiv_Demkowicz_Interior_MBTET(int p, double *N, double *diffN,
                                             int p_face[], double *phi_f[4],
                                             double *diff_phi_f[4],
                                             double *phi_v, double *diff_phi_v,
                                             int gdim);

} // namespace MoFEM

#endif // __HDIV_HPP__
