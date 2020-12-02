/** \file h1_hdiv_hcurl_l2.h
\brief Functions to approximate hierarchical spaces

\FIXME: Name Shape Functions is used, in that context is more appropriate
to use base functions. Need to be changed.

*/

/*
 * MoFEM is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
 * License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with MoFEM. If not, see <http://www.gnu.org/licenses/>
 */

#ifndef __H1_HDIV_HCURL_L2_H__
#define __H1_HDIV_HCURL_L2_H__

#ifndef DEPRECATED
  #define DEPRECATED
#endif

#ifdef __cplusplus
extern "C" {
#endif

// L2

// Number of dofs for L2 space

/**
 * @brief Number of base functions on tetrahedron for L2 space
 */
#define NBVOLUMETET_L2(P) ((P + 1) * (P + 2) * (P + 3) / 6)

/**
 * @brief Number of base functions on triangle for L2 space
 */
#define NBFACETRI_L2(P) ((P + 1) * (P + 2) / 2)

/**
 * @brief Number of base functions on edge fro L2 space
 * 
 */
#define NBEDGE_L2(P) (P + 1)

/**
 * @brief Number of base functions on edge for L2 (Demkowicz) space
 *
 */
#define NBEDGE_EXACT_L2(P) (P)

// H1

/**
 * @brief Numer of base function on edge for H1 space
 */
#define NBEDGE_H1(P) (((P) > 1) ? (P - 1) : 0)

/**
 * @brief Number of base function on triangle for H1 space
 */
#define NBFACETRI_H1(P) (((P) > 2) ? ((P - 2) * (P - 1) / 2) : 0)

/**
 * @brief Number of base functions on quad for H1 space
 */
#define NBFACEQUAD_H1(P) (((P) > 1) ? ((P - 1) * (P - 1)) : 0)

/**
 * @brief Number of base functions on quad for L2 space
 */
#define NBFACEQUAD_L2(P) (((P) > 1) ? ((P - 1) * (P - 1)) : 0)

/**
 * @brief Number of base functions on quad edge for Hcurl space
 */
#define NBEDGEQUAD_FULL_HCURL(P) (((P) > 0) ? P : 0)

/**
 * @brief Number of base functions on quad for Hcurl space
 */
#define NBFACE_TYP_QUAD_FULL_HCURL(P) (((P) > 1) ? P * (P - 1) : 0)
/**
 * @brief Number of base functions on quad for H1 space
 */
#define NBFACEQUAD_FULL_H1(P) (((P) > 1) ? ((P - 1) * (P - 1)) : 0)

/**
 * @brief Number of base functions on tetrahedron for H1 space
 */
#define NBVOLUMETET_H1(P) (((P) > 3) ? ((P - 3) * (P - 2) * (P - 1) / 6) : 0)

/**
 * @brief Number of base functions on prism for H1 space
 */
#define NBVOLUMEPRISM_H1(P) ((P > 4) ? ((P - 2) * (P - 3) * (P - 4) / 6) : 0)

// H curl

#define NBEDGE_AINSWORTH_HCURL(P) (((P) > 0) ? (P + 1) : 0)
#define NBFACETRI_AINSWORTH_EDGE_HCURL(P) (((P) > 1) ? P - 1 : 0)
#define NBFACETRI_AINSWORTH_FACE_HCURL(P) (((P) > 2) ? (P - 1) * (P - 2) : 0)
#define NBFACETRI_AINSWORTH_HCURL(P) ((P > 1) ? ((P)-1) * (P + 1) : 0)
#define NBVOLUMETET_AINSWORTH_FACE_HCURL(P)                                    \
  (((P) > 2) ? (2 * (P - 1) * (P - 2)) : 0)
#define NBVOLUMETET_AINSWORTH_TET_HCURL(P)                                     \
  (((P) > 3) ? ((P - 3) * (P - 2) * (P - 1) / 2) : 0)
#define NBVOLUMETET_AINSWORTH_HCURL(P)                                         \
  (((P) > 2) ? (P - 2) * (P - 1) * (P + 1) / 2 : 0)

#define NBEDGE_DEMKOWICZ_HCURL(P) (((P) > 0) ? (P) : 0)
#define NBFACETRI_DEMKOWICZ_HCURL(P) (((P) > 1) ? (P) * ((P)-1) : 0)
#define NBVOLUMETET_DEMKOWICZ_HCURL(P)                                         \
  (((P) > 2) ? ((P) * ((P)-1) * ((P)-2) / 2) : 0)

// H div

#define NBEDGE_HDIV(P) (0)
#define NBFACETRI_AINSWORTH_EDGE_HDIV(P) (((P) > 0) ? (P) : 0)
#define NBFACETRI_AINSWORTH_FACE_HDIV(P) (((P) > 2) ? (P - 1) * (P - 2) / 2 : 0)
#define NBFACETRI_AINSWORTH_HDIV(P) (((P) > 0) ? (P + 1) * (P + 2) / 2 : 0)
#define NBVOLUMETET_AINSWORTH_EDGE_HDIV(P) (((P) > 1) ? (P - 1) : 0)
#define NBVOLUMETET_AINSWORTH_FACE_HDIV(P) (((P) > 2) ? (P - 1) * (P - 2) : 0)
#define NBVOLUMETET_AINSWORTH_VOLUME_HDIV(P)                                   \
  (((P) > 3) ? (P - 3) * (P - 2) * (P - 1) / 2 : 0)
#define NBVOLUMETET_AINSWORTH_HDIV(P)                                          \
  (((P) > 1) ? (P - 1) * (P + 1) * (P + 2) / 2 : 0)
#define NBFACETRI_DEMKOWICZ_HDIV(P) ((P > 0) ? (P) * (P + 1) / 2 : 0)
#define NBVOLUMETET_DEMKOWICZ_HDIV(P)                                          \
  (((P) > 1) ? (P) * (P - 1) * (P + 1) / 2 : 0)

// Bubbles for H div space

/**
 * @brief Get base functions on triangle for L2 space
 * 
 * @param p polynomial order
 * @param N barycentric coordinates (shape functions) at integration points 
 * @param diffN direvatives of barycentric coordinates, i.e. direvatives of shape functions 
 * @param L2N values of L2 base at integration points
 * @param diff_L2N dirvatives of base functions at integration points
 * @param GDIM number of integration points
 * @param base_polynomials polynomial base used to construct L2 base on element 
 * @return PetscErrorCode 
 */
PetscErrorCode L2_Ainsworth_ShapeFunctions_MBTRI(
    int p, double *N, double *diffN, double *L2N, double *diff_L2N, int GDIM,
    PetscErrorCode (*base_polynomials)(int p, double s, double *diff_s,
                                       double *L, double *diffL,
                                       const int dim));
/**
 * @brief Get base functions on tetrahedron for L2 space
 * 
 * @param p polynomial order
 * @param N barycentric coordinates (shape functions) at integration points 
 * @param diffN direvatives of barycentric coordinates, i.e. direvatives of shape functions 
 * @param L2N values of L2 base at integration points
 * @param diff_L2N dirvatives of base functions at integration points
 * @param GDIM number of integration points
 * @param base_polynomials polynomial base used to construct L2 base on element 
 * @return PetscErrorCode 
 */
PetscErrorCode L2_Ainsworth_ShapeFunctions_MBTET(
    int p, double *N, double *diffN, double *L2N, double *diff_L2N, int GDIM,
    PetscErrorCode (*base_polynomials)(int p, double s, double *diff_s,
                                       double *L, double *diffL,
                                       const int dim));

/**
 * \deprecated Use L2_Ainsworth_ShapeFunctions_MBTRI 
 */
DEPRECATED PetscErrorCode L2_ShapeFunctions_MBTRI(
    int p, double *N, double *diffN, double *L2N, double *diff_L2N, int GDIM,
    PetscErrorCode (*base_polynomials)(int p, double s, double *diff_s,
                                       double *L, double *diffL,
                                       const int dim));
/**
 * \deprecated Use L2_Ainsworth_ShapeFunctions_MBTET
 */
DEPRECATED PetscErrorCode L2_ShapeFunctions_MBTET(
    int p, double *N, double *diffN, double *L2N, double *diff_L2N, int GDIM,
    PetscErrorCode (*base_polynomials)(int p, double s, double *diff_s,
                                       double *L, double *diffL,
                                       const int dim));

PetscErrorCode L2_VolumeShapeDiffMBTETinvJ(int base_p, int p,
                                           double *volume_diffN, double *invJac,
                                           double *volume_diffNinvJac,
                                           int GDIM);

/**
 * \brief H1_EdgeShapeFunctions_MBTRI
 *
 * \param sense of edges, it is array of integers dim 3 (3-edges of triangle)
 * \param p of edges
 */
PetscErrorCode H1_EdgeShapeFunctions_MBTRI(
    int *sense, int *p, double *N, double *diffN, double *edgeN[3],
    double *diff_edgeN[3], int GDIM,
    PetscErrorCode (*base_polynomials)(int p, double s, double *diff_s,
                                       double *L, double *diffL,
                                       const int dim));
PetscErrorCode H1_FaceShapeFunctions_MBTRI(
    const int *face_nodes, int p, double *N, double *diffN, double *faceN,
    double *diff_faceN, int GDIM,
    PetscErrorCode (*base_polynomials)(int p, double s, double *diff_s,
                                       double *L, double *diffL,
                                       const int dim));
PetscErrorCode H1_EdgeShapeFunctions_MBTET(
    int *sense, int *p, double *N, double *diffN, double *edgeN[],
    double *diff_edgeN[], int GDIM,
    PetscErrorCode (*base_polynomials)(int p, double s, double *diff_s,
                                       double *L, double *diffL,
                                       const int dim));
PetscErrorCode H1_FaceShapeFunctions_MBTET(
    int *faces_nodes, int *p, double *N, double *diffN, double *faceN[],
    double *diff_faceN[], int GDIM,
    PetscErrorCode (*base_polynomials)(int p, double s, double *diff_s,
                                       double *L, double *diffL,
                                       const int dim));
PetscErrorCode H1_VolumeShapeFunctions_MBTET(
    int p, double *N, double *diffN, double *volumeN, double *diff_volumeN,
    int GDIM,
    PetscErrorCode (*base_polynomials)(int p, double s, double *diff_s,
                                       double *L, double *diffL,
                                       const int dim));
PetscErrorCode H1_EdgeShapeDiffMBTETinvJ(int *base_p, int *p,
                                         double *edge_diffN[], double *invJac,
                                         double *edge_diffNinvJac[], int GDIM);
PetscErrorCode H1_FaceShapeDiffMBTETinvJ(int *base_p, int *p,
                                         double *face_diffN[], double *invJac,
                                         double *face_diffNinvJac[], int GDIM);
PetscErrorCode H1_VolumeShapeDiffMBTETinvJ(int base_p, int p,
                                           double *volume_diffN, double *invJac,
                                           double *volume_diffNinvJac,
                                           int GDIM);

PetscErrorCode H1_EdgeGradientOfDeformation_hierarchical(int p, double *diffN,
                                                        double *dofs,
                                                        double *F);
PetscErrorCode H1_FaceGradientOfDeformation_hierarchical(int p, double *diffN,
                                                        double *dofs,
                                                        double *F);
PetscErrorCode H1_VolumeGradientOfDeformation_hierarchical(int p, double *diffN,
                                                      double *dofs, double *F);

/**
 * \deprecated use H1_EdgeGradientOfDeformation_hierarchical
 */
DEPRECATED PetscErrorCode H1_EdgeGradientOfDeformation_hierachical(
    int p, double *diffN, double *dofs, double *F);

/**
 * \deprecated use H1_FaceGradientOfDeformation_hierarchical
 */
DEPRECATED PetscErrorCode H1_FaceGradientOfDeformation_hierachical(
    int p, double *diffN, double *dofs, double *F);

/**
 * \deprecated use H1_VolumeGradientOfDeformation_hierarchical
 */
DEPRECATED PetscErrorCode H1_VolumeGradientOfDeformation_hierachical(
    int p, double *diffN, double *dofs, double *F);

PetscErrorCode H1_QuadShapeFunctions_MBPRISM(
    int *faces_nodes, int *p, double *N, double *diffN, double *faceN[],
    double *diff_faceN[], int GDIM,
    PetscErrorCode (*base_polynomials)(int p, double s, double *diff_s,
                                       double *L, double *diffL,
                                       const int dim));
PetscErrorCode H1_VolumeShapeFunctions_MBPRISM(
    int p, double *N, double *diffN, double *volumeN, double *diff_volumeN,
    int GDIM,
    PetscErrorCode (*base_polynomials)(int p, double s, double *diff_s,
                                       double *L, double *diffL,
                                       const int dim));
PetscErrorCode H1_QuadShapeFunctions_MBQUAD(
    int *faces_nodes, int p, double *N, double *diffN, double *faceN,
    double *diff_faceN, int GDIM,
    PetscErrorCode (*base_polynomials)(int p, double s, double *diff_s,
                                       double *L, double *diffL,
                                       const int dim));
PetscErrorCode H1_EdgeShapeFunctions_MBQUAD(
    int *sense, int *p, double *N, double *diffN, double *edgeN[4],
    double *diff_edgeN[4], int GDIM,
    PetscErrorCode (*base_polynomials)(int p, double s, double *diff_s,
                                       double *L, double *diffL, const int dim));

// Hdiv and Hcurl are implemented and declared in other files



#ifdef __cplusplus
}
#endif

#endif // __H1_HDIV_HCURL_L2_H__
