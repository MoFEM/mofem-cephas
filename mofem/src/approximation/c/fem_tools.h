/** \file fem_tools.h
 * \brief Loose implementation of some useful functions
 *
 * FIXME: Implementation here is very unstructured, need cleaning and pruning
 *
 */

#ifndef __FEM_H__
#define __FEM_H__

#include <petscsys.h>

#ifdef __cplusplus
extern "C" {
#endif
#include <cblas.h>
#include <lapack_wrap.h>
#ifdef __cplusplus
}
#endif

#define LAMBDA(E, NU) (E * NU / ((1. + NU) * (1. - 2. * NU)))
#define MU(E, NU) (0.5 * E / (1. + NU))
#define DELTA(NU_P, NU_PZ, E_P, E_Z)                                           \
  (((1 + NU_P) * (1 - NU_P - 2 * NU_PZ * (NU_PZ * E_Z / E_P))) /               \
   (E_P * E_P * E_Z))

#define N_MBTET0(x, y, z) (1. - x - y - z) ///< tetrahedral shape function
#define N_MBTET1(x, y, z) (x)              ///< tetrahedral shape function
#define N_MBTET2(x, y, z) (y)              ///< tetrahedral shape function
#define N_MBTET3(x, y, z) (z)              ///< tetrahedral shape function
#define diffN_MBTET0x (-1.) ///< derivative of tetrahedral shape function
#define diffN_MBTET0y (-1.) ///< derivative of tetrahedral shape function
#define diffN_MBTET0z (-1.) ///< derivative of tetrahedral shape function
#define diffN_MBTET1x (1.)  ///< derivative of tetrahedral shape function
#define diffN_MBTET1y (0.)  ///< derivative of tetrahedral shape function
#define diffN_MBTET1z (0.)  ///< derivative of tetrahedral shape function
#define diffN_MBTET2x (0.)  ///< derivative of tetrahedral shape function
#define diffN_MBTET2y (1.)  ///< derivative of tetrahedral shape function
#define diffN_MBTET2z (0.)  ///< derivative of tetrahedral shape function
#define diffN_MBTET3x (0.)  ///< derivative of tetrahedral shape function
#define diffN_MBTET3y (0.)  ///< derivative of tetrahedral shape function
#define diffN_MBTET3z (1.)  ///< derivative of tetrahedral shape function

// MBTRI
#define N_MBTRI0(x, y) (1. - x - y) ///< triangle shape function
#define N_MBTRI1(x, y) (x)          ///< triangle shape function
#define N_MBTRI2(x, y) (y)          ///< triangle shape function
#define diffN_MBTRI0x (-1.)         ///< derivative of triangle shape function
#define diffN_MBTRI0y (-1.)         ///< derivative of triangle shape function
#define diffN_MBTRI1x (1.)          ///< derivative of triangle shape function
#define diffN_MBTRI1y (0.)          ///< derivative of triangle shape function
#define diffN_MBTRI2x (0.)          ///< derivative of triangle shape function
#define diffN_MBTRI2y (1.)          ///< derivative of triangle shape function

// MBQUAD
#define N_MBQUAD0(x, y) ((1. - x) * (1. - y)) ///< quad shape function
#define N_MBQUAD1(x, y) ((x) * (1. - y))         ///< quad shape function
#define N_MBQUAD2(x, y) ((x) * (y))             ///< quad shape function
#define N_MBQUAD3(x, y) ((1. - x) * (y))      ///< quad shape function
#define diffN_MBQUAD0x(y) (-(1. - y))
#define diffN_MBQUAD0y(x) (-(1. - x))
#define diffN_MBQUAD1x(y) ((1. - y))
#define diffN_MBQUAD1y(x) (-x)
#define diffN_MBQUAD2x(y) (y)
#define diffN_MBQUAD2y(x) (x)
#define diffN_MBQUAD3x(y) (-y)
#define diffN_MBQUAD3y(x) ((1. - x))

// MBHEX
#define N_MBHEX0(x, y, z) (N_MBQUAD0(x, y) * (1 - z))
#define N_MBHEX1(x, y, z) (N_MBQUAD1(x, y) * (1 - z))
#define N_MBHEX2(x, y, z) (N_MBQUAD2(x, y) * (1 - z))
#define N_MBHEX3(x, y, z) (N_MBQUAD3(x, y) * (1 - z))
#define N_MBHEX4(x, y, z) (N_MBQUAD0(x, y) * (z))
#define N_MBHEX5(x, y, z) (N_MBQUAD1(x, y) * (z))
#define N_MBHEX6(x, y, z) (N_MBQUAD2(x, y) * (z))
#define N_MBHEX7(x, y, z) (N_MBQUAD3(x, y) * (z))
#define diffN_MBHEX0x(y, z) (diffN_MBQUAD0x(y) * (1 - z))
#define diffN_MBHEX1x(y, z) (diffN_MBQUAD1x(y) * (1 - z))
#define diffN_MBHEX2x(y, z) (diffN_MBQUAD2x(y) * (1 - z))
#define diffN_MBHEX3x(y, z) (diffN_MBQUAD3x(y) * (1 - z))
#define diffN_MBHEX4x(y, z) (diffN_MBQUAD0x(y) * (z))
#define diffN_MBHEX5x(y, z) (diffN_MBQUAD1x(y) * (z))
#define diffN_MBHEX6x(y, z) (diffN_MBQUAD2x(y) * (z))
#define diffN_MBHEX7x(y, z) (diffN_MBQUAD3x(y) * (z))
#define diffN_MBHEX0y(x, z) (diffN_MBQUAD0y(x) * (1 - z))
#define diffN_MBHEX1y(x, z) (diffN_MBQUAD1y(x) * (1 - z))
#define diffN_MBHEX2y(x, z) (diffN_MBQUAD2y(x) * (1 - z))
#define diffN_MBHEX3y(x, z) (diffN_MBQUAD3y(x) * (1 - z))
#define diffN_MBHEX4y(x, z) (diffN_MBQUAD0y(x) * (z))
#define diffN_MBHEX5y(x, z) (diffN_MBQUAD1y(x) * (z))
#define diffN_MBHEX6y(x, z) (diffN_MBQUAD2y(x) * (z))
#define diffN_MBHEX7y(x, z) (diffN_MBQUAD3y(x) * (z))
#define diffN_MBHEX0z(x, y) (-N_MBQUAD0(x, y))
#define diffN_MBHEX1z(x, y) (-N_MBQUAD1(x, y))
#define diffN_MBHEX2z(x, y) (-N_MBQUAD2(x, y))
#define diffN_MBHEX3z(x, y) (-N_MBQUAD3(x, y))
#define diffN_MBHEX4z(x, y) (N_MBQUAD0(x, y))
#define diffN_MBHEX5z(x, y) (N_MBQUAD1(x, y))
#define diffN_MBHEX6z(x, y) (N_MBQUAD2(x, y))
#define diffN_MBHEX7z(x, y) (N_MBQUAD3(x, y))

// MBEDGE
#define N_MBEDGE0(x) (1. - (x)) ///< edge shape function
#define N_MBEDGE1(x) (x)        ///< edge shape function
#define diffN_MBEDGE0 (-1.)     ///< derivative of edge shape function
#define diffN_MBEDGE1 (1.)      ///< derivative of edge shape function

// MBTRIQ
#define N_MBTRIQ0(x, y) ((1. - x - y) * (2 * (1. - x - y) - 1.))
#define N_MBTRIQ1(x, y) (x * (2. * x - 1.))
#define N_MBTRIQ2(x, y) (y * (2. * y - 1.))
#define N_MBTRIQ3(x, y) (4. * (1. - x - y) * x)
#define N_MBTRIQ4(x, y) (4. * x * y)
#define N_MBTRIQ5(x, y) (4. * (1. - x - y) * y)
#define diffN_MBTRIQ0x(x, y) (x + y - 3. * (1. - x - y))
#define diffN_MBTRIQ0y(x, y) (x + y - 3. * (1. - x - y))
#define diffN_MBTRIQ1x(x, y) (-1. + 4. * x)
#define diffN_MBTRIQ1y(x, y) (0.)
#define diffN_MBTRIQ2x(x, y) (0.)
#define diffN_MBTRIQ2y(x, y) (-1. + 4. * y)
#define diffN_MBTRIQ3x(x, y) (4. * ((1. - x - y) - x))
#define diffN_MBTRIQ3y(x, y) (-4. * x)
#define diffN_MBTRIQ4x(x, y) (4. * y)
#define diffN_MBTRIQ4y(x, y) (4. * x)
#define diffN_MBTRIQ5x(x, y) (-4. * y)
#define diffN_MBTRIQ5y(x, y) (4. * ((1. - x - y) - y))

#ifdef __cplusplus
extern "C" {
#endif

/// print matric M
void print_mat(double *M, int m, int n);
/// print upper part of the symmetric matrix
void print_mat_sym_upper(double *M, int m, int n);
/// priint complex matrix
void print_mat_complex(__CLPK_doublecomplex *M, int m, int n);

/// \brief calculate shape functions on triangle
/// \param N shape function array
/// \param X array of Gauss X coordinates
/// \param Y array of Gauss Y coordinates
/// \param G_DIM number of Gauss points
PetscErrorCode ShapeMBTRI(double *N, const double *X, const double *Y,
                          const int G_DIM);
/// calculate derivatives of shape functions
PetscErrorCode ShapeDiffMBTRI(double *diffN);

/// calculate face normal
/// \param diffN derivatives of shape functions
/// \param coords is position of the nodes
/// \param normal vector
PetscErrorCode ShapeFaceNormalMBTRI(double *diffN, const double *coords,
                                    double *normal);
PetscErrorCode ShapeFaceBaseMBTRI(double *diffN, const double *coords,
                                  double *normal, double *s1, double *s2);

/// calculate derivative of normal in respect to nodal positions
PetscErrorCode ShapeFaceDiffNormalMBTRI(double *diffN, const double *coords,
                                        double *diff_normal);
/// calculate jacobioan
void ShapeJacMBTRI(double *diffN, const double *coords, double *Jac);
/// calculate derivatives of shape functions in space
void ShapeDiffMBTRIinvJ(double *diffN, double *invJac, double *diffNinvJac);
/// calculate shape functions
PetscErrorCode ShapeMBTET(double *N, const double *G_X, const double *G_Y,
                          const double *G_Z, int DIM);
/// calculate derivatives of shape functions
PetscErrorCode ShapeDiffMBTET(double *diffN);
/// determined of jacobian
double ShapeDetJacVolume(double *jac);
/// calculate jacobian
PetscErrorCode ShapeJacMBTET(double *diffN, const double *coords, double *jac);
// calculate inverse of jacobian
PetscErrorCode ShapeInvJacVolume(double *jac);
/// calculate TET volume
double ShapeVolumeMBTET(double *diffN, const double *coords);
/// calculate shape functions derivatives in space
PetscErrorCode ShapeDiffMBTETinvJ(double *diffN, double *invJac,
                                  double *diffNinvJac);

/// calculate spin matrix from vector
// \param spinOmega is a spin matrix
// \param vecOmega is a spin vector
PetscErrorCode Spin(double *spinOmega, double *vecOmega);

/// Compose complex matrix (3x3) from two real matrices
PetscErrorCode make_complex_matrix(double *reA, double *imA,
                                   __CLPK_doublecomplex *xA);
/// Complex normal
PetscErrorCode
Normal_hierarchical(int order_approx, int *order_edge_approx, int order,
                    int *order_edge, double *diffN, double *diffN_face,
                    double *diffN_edge[], double *dofs, double *dofs_edge[],
                    double *dofs_face, double *idofs, double *idofs_edge[],
                    double *idofs_face, __CLPK_doublecomplex *xnormal,
                    __CLPK_doublecomplex *s1, __CLPK_doublecomplex *s2, int gg);
PetscErrorCode Base_scale(__CLPK_doublecomplex *xnormal,
                          __CLPK_doublecomplex *xs1, __CLPK_doublecomplex *xs2);

/**
 * \brief calculate local coordinates for given global coordinates
 *
 * new version for multiple points need to be implemented
 */
PetscErrorCode ShapeMBTET_inverse(double *N, double *diffN,
                                  const double *elem_coords,
                                  const double *glob_coords,
                                  double *loc_coords);

/**
 * \brief calculate local coordinates of triangle element for given global
coordinates in 2D (Assume e.g. z=0) \f[ \left[\begin{array}{cc}
\frac{\partial N_{1}}{\partial\xi}x_{N_{1}}+\frac{\partial
N_{2}}{\partial\xi}x_{N_{2}}+\frac{\partial N_{3}}{\partial\xi}x_{N_{3}} &
\frac{\partial N_{1}}{\partial\eta}x_{N_{1}}+\frac{\partial
N_{2}}{\partial\eta}x_{N_{2}}+\frac{\partial N_{3}}{\partial\eta}x_{N_{3}}\\
\frac{\partial N_{1}}{\partial\xi}y_{N_{1}}+\frac{\partial
N_{2}}{\partial\xi}y_{N_{2}}+\frac{\partial N_{3}}{\partial\xi}y_{N_{3}} &
\frac{\partial N_{1}}{\partial\eta}y_{N_{1}}+\frac{\partial
N_{2}}{\partial\eta}y_{N_{2}}+\frac{\partial N_{3}}{\partial\eta}y_{N_{3}}
\end{array}\right]\left\{ \begin{array}{c}
\xi\\
\eta
\end{array}\right\} =\left\{ \begin{array}{c}
x_{gp}-\left(N_{1}x_{N_{1}}+N_{2}x_{N_{2}}+N_{3}x_{N_{3}}\right)\\
y_{gp}-\left(N_{1}y_{N_{1}}+N_{2}y_{N_{2}}+N_{3}y_{N_{3}}\right)
\end{array}\right\}
 \f]

 /// \param N shape function array
 /// \param diffN array of shape function derivative w.r.t to local coordinates
 /// \param elem_coords global coordinates of element
 /// \param glob_coords global coordinates of required point
 /// \param loc_coords  local coordinates of required point
 */
PetscErrorCode ShapeMBTRI_inverse(double *N, double *diffN,
                                  const double *elem_coords,
                                  const double *glob_coords,
                                  double *loc_coords);

/// calculate gradient of deformation
PetscErrorCode GradientOfDeformation(double *diffN, double *dofs, double *F);

// 2 Node edge
PetscErrorCode ShapeMBEDGE(double *N, const double *G_X, int DIM);
PetscErrorCode ShapeDiffMBEDGE(double *diffN);

// 10 Node Tet
PetscErrorCode ShapeMBTRIQ(double *N, const double *X, const double *Y,
                           const int G_DIM);
PetscErrorCode ShapeDiffMBTRIQ(double *diffN, const double *X, const double *Y,
                               const int G_DIM);
PetscErrorCode ShapeMBTETQ(double *N, const double x, const double y,
                           const double z);
PetscErrorCode ShapeDiffMBTETQ(double *diffN, const double x, const double y,
                               const double z);
PetscErrorCode ShapeMBTETQ_GAUSS(double *N, const double *X, const double *Y,
                                 const double *Z, const int G_DIM);
PetscErrorCode ShapeDiffMBTETQ_GAUSS(double *diffN, const double *X,
                                     const double *Y, const double *Z,
                                     const int G_DIM);
PetscErrorCode ShapeJacMBTETQ(const double *diffN, const double *coords,
                              double *Jac);
PetscErrorCode
ShapeMBTETQ_detJac_at_Gauss_Points(double *detJac_at_Gauss_Points,
                                   const double *diffN, const double *coords,
                                   int G_DIM);
double ShapeVolumeMBTETQ(const double *diffN, const double *coords, int G_DIM,
                         double *G_TET_W);
PetscErrorCode ShapeMBTETQ_inverse(double *N, double *diffN,
                                   const double *elem_coords,
                                   const double *glob_coords,
                                   double *loc_coords, const double eps);

// complex part
void ShapeDiffMBTETinvJ_complex(double *diffN, __CLPK_doublecomplex *invJac,
                                __CLPK_doublecomplex *diffNinvJac,
                                enum CBLAS_TRANSPOSE Trans);
PetscErrorCode ShapeFaceNormalMBTRI_complex(double *diffN,
                                            __CLPK_doublecomplex *xcoords,
                                            __CLPK_doublecomplex *xnormal);
PetscErrorCode MakeComplexTensor(double *reA, double *imA,
                                 __CLPK_doublecomplex *xA);
PetscErrorCode InvertComplexGradient(__CLPK_doublecomplex *xF);
PetscErrorCode DeterminantComplexGradient(__CLPK_doublecomplex *xF,
                                          __CLPK_doublecomplex *det_xF);

// integration
/// Compute weights and integration points for edge using Grundmann_Moeller rule
PetscErrorCode Grundmann_Moeller_integration_points_1D_EDGE(int rule,
                                                            double *G_TRI_X,
                                                            double *G_TRI_W);
/// Compute weights and integration points for 2D Triangle using
/// Grundmann_Moeller rule
PetscErrorCode Grundmann_Moeller_integration_points_2D_TRI(int rule,
                                                           double *G_TRI_X,
                                                           double *G_TRI_Y,
                                                           double *G_TRI_W);
/// Compute weights and integration points for 3D Tet using Grundmann_Moeller
/// rule
PetscErrorCode Grundmann_Moeller_integration_points_3D_TET(int rule,
                                                           double *G_TET_X,
                                                           double *G_TET_Y,
                                                           double *G_TET_Z,
                                                           double *G_TET_W);

// http://people.sc.fsu.edu/~jburkardt/datasets/quadrature_rules_tet/quadrature_rules_tet.html
// TRI
static const double G_TRI_X1[] = {3.3333333333333331e-01};
static const double G_TRI_Y1[] = {3.3333333333333331e-01};
static const double G_TRI_W1[] = {1.};
static const double G_TRI_X3[] = {0.5, 0., 0.5};
static const double G_TRI_Y3[] = {0., 0.5, 0.5};
static const double G_TRI_W3[] = {
    3.3333333333333331e-01, 3.3333333333333331e-01, 3.3333333333333331e-01};
static const double G_TRI_X4[] = {
    7.503111022260811058e-02, 1.785587282636164064e-01,
    2.800199154990741235e-01, 6.663902460147014262e-01};
static const double G_TRI_Y4[] = {
    2.800199154990741235e-01, 6.663902460147014262e-01,
    7.503111022260811058e-02, 1.785587282636164064e-01};
static const double G_TRI_W4[] = {
    1.8195861825602258066e-01, 3.1804138174397683647e-01,
    1.8195861825602258066e-01, 3.1804138174397683647e-01};
static const double G_TRI_X7[] = {
    0.333333333333333, 0.736712498968435, 0.736712498968435, 0.237932366472434,
    0.237932366472434, 0.025355134551932, 0.025355134551932};
static const double G_TRI_Y7[] = {
    0.333333333333333, 0.237932366472434, 0.025355134551932, 0.736712498968435,
    0.025355134551932, 0.736712498968435, 0.237932366472434};
static const double G_TRI_W7[] = {
    0.375000000000000, 0.104166666666667, 0.104166666666667, 0.104166666666667,
    0.104166666666667, 0.104166666666667, 0.104166666666667};
static const double G_TRI_X13[] = {
    0.333333333333333, 0.479308067841923, 0.260345966079038, 0.260345966079038,
    0.869739794195568, 0.065130102902216, 0.065130102902216, 0.638444188569809,
    0.638444188569809, 0.312865496004875, 0.312865496004875, 0.048690315425316,
    0.048690315425316};
static const double G_TRI_Y13[] = {
    0.333333333333333, 0.260345966079038, 0.479308067841923, 0.260345966079038,
    0.065130102902216, 0.869739794195568, 0.065130102902216, 0.312865496004875,
    0.048690315425316, 0.638444188569809, 0.048690315425316, 0.638444188569809,
    0.312865496004875};
static const double G_TRI_W13[] = {
    -0.149570044467670, 0.175615257433204, 0.175615257433204, 0.175615257433204,
    0.053347235608839,  0.053347235608839, 0.053347235608839, 0.077113760890257,
    0.077113760890257,  0.077113760890257, 0.077113760890257, 0.077113760890257,
    0.077113760890257};

static const double G_TRI_X19[] = {
    0.333333333333333, 0.797426985353087, 0.101286507323456, 0.101286507323456,
    0.059715871789770, 0.470142064105115, 0.470142064105115, 0.535795346449899,
    0.232102326775050, 0.232102326775050, 0.941038278231121, 0.029480860884440,
    0.029480860884440, 0.738416812340510, 0.738416812340510, 0.232102326775050,
    0.232102326775050, 0.029480860884440, 0.029480860884440};
static const double G_TRI_Y19[] = {
    0.333333333333333, 0.101286507323456, 0.797426985353087, 0.101286507323456,
    0.470142064105115, 0.059715871789770, 0.470142064105115, 0.232102326775050,
    0.535795346449899, 0.232102326775050, 0.029480860884440, 0.941038278231121,
    0.029480860884440, 0.232102326775050, 0.029480860884440, 0.738416812340510,
    0.029480860884440, 0.738416812340510, 0.232102326775050};
static const double G_TRI_W19[] = {
    9.71357962827961025E-002, 3.13347002271398278E-002,
    3.13347002271398278E-002, 3.13347002271398278E-002,
    7.78275410047754301E-002, 7.78275410047754301E-002,
    7.78275410047754301E-002, 7.96477389272090969E-002,
    7.96477389272090969E-002, 7.96477389272090969E-002,
    2.55776756586981006E-002, 2.55776756586981006E-002,
    2.55776756586981006E-002, 4.32835393772893970E-002,
    4.32835393772893970E-002, 4.32835393772893970E-002,
    4.32835393772893970E-002, 4.32835393772893970E-002,
    4.32835393772893970E-002};

static const double G_TRI_X28[] = {
    0.333333333333333,  0.948021718143423,  0.025989140928288,
    0.025989140928288,  0.811424994704155,  0.094287502647923,
    0.094287502647923,  0.010726449965571,  0.494636775017215,
    0.494636775017215,  0.585313234770972,  0.207343382614514,
    0.207343382614514,  0.122184388599019,  0.438907805700491,
    0.438907805700491,  0.677937654882590,  0.677937654882590,
    0.044841677589131,  0.044841677589131,  0.277220667528279,
    0.277220667528279,  0.858870281282636,  0.858870281282636,
    0.0000000000000000, 0.0000000000000000, 0.141129718717364,
    0.141129718717364};
static const double G_TRI_Y28[] = {
    0.333333333333333, 0.025989140928288, 0.948021718143423, 0.025989140928288,
    0.094287502647923, 0.811424994704155, 0.094287502647923, 0.494636775017215,
    0.010726449965571, 0.494636775017215, 0.207343382614514, 0.585313234770972,
    0.207343382614514, 0.438907805700491, 0.122184388599019, 0.438907805700491,
    0.044841677589131, 0.277220667528279, 0.677937654882590, 0.277220667528279,
    0.677937654882590, 0.044841677589131, 0.000000000000000, 0.141129718717364,
    0.858870281282636, 0.141129718717364, 0.858870281282636, 0.000000000000000};
static const double G_TRI_W28[] = {
    0.08797730116222190,  0.008744311553736190, 0.008744311553736190,
    0.008744311553736190, 0.03808157199393533,  0.03808157199393533,
    0.03808157199393533,  0.01885544805613125,  0.01885544805613125,
    0.01885544805613125,  0.07215969754474100,  0.07215969754474100,
    0.07215969754474100,  0.06932913870553720,  0.06932913870553720,
    0.06932913870553720,  0.04105631542928860,  0.04105631542928860,
    0.04105631542928860,  0.04105631542928860,  0.04105631542928860,
    0.04105631542928860,  0.007362383783300573, 0.007362383783300573,
    0.007362383783300573, 0.007362383783300573, 0.007362383783300573,
    0.007362383783300573};

static const double G_TRI_X37[] = {
    0.333333333333333, 0.950275662924106, 0.024862168537947, 0.024862168537947,
    0.171614914923835, 0.414192542538082, 0.414192542538082, 0.539412243677190,
    0.230293878161405, 0.230293878161405, 0.772160036676533, 0.113919981661734,
    0.113919981661734, 0.009085399949835, 0.495457300025082, 0.495457300025082,
    0.062277290305887, 0.468861354847056, 0.468861354847056, 0.022076289653624,
    0.022076289653624, 0.851306504174348, 0.851306504174348, 0.126617206172027,
    0.126617206172027, 0.018620522802521, 0.018620522802521, 0.689441970728591,
    0.689441970728591, 0.291937506468888, 0.291937506468888, 0.096506481292159,
    0.096506481292159, 0.635867859433873, 0.635867859433873, 0.267625659273968,
    0.267625659273968};
static const double G_TRI_Y37[] = {
    0.333333333333333, 0.024862168537947, 0.950275662924106, 0.024862168537947,
    0.414192542538082, 0.171614914923835, 0.414192542538082, 0.230293878161405,
    0.539412243677190, 0.230293878161405, 0.113919981661734, 0.772160036676533,
    0.113919981661734, 0.495457300025082, 0.009085399949835, 0.495457300025082,
    0.468861354847056, 0.062277290305887, 0.468861354847056, 0.851306504174348,
    0.126617206172027, 0.022076289653624, 0.126617206172027, 0.022076289653624,
    0.851306504174348, 0.689441970728591, 0.291937506468888, 0.018620522802521,
    0.291937506468888, 0.018620522802521, 0.689441970728591, 0.635867859433873,
    0.267625659273968, 0.096506481292159, 0.267625659273968, 0.096506481292159,
    0.635867859433873};
static const double G_TRI_W37[] = {
    0.051739766065744, 0.008007799555565, 0.008007799555565, 0.008007799555565,
    0.046868898981822, 0.046868898981822, 0.046868898981822, 0.046590940183976,
    0.046590940183976, 0.046590940183976, 0.031016943313796, 0.031016943313796,
    0.031016943313796, 0.010791612736631, 0.010791612736631, 0.010791612736631,
    0.032195534242432, 0.032195534242432, 0.032195534242432, 0.015445834210702,
    0.015445834210702, 0.015445834210702, 0.015445834210702, 0.015445834210702,
    0.015445834210702, 0.017822989923179, 0.017822989923179, 0.017822989923179,
    0.017822989923179, 0.017822989923179, 0.017822989923179, 0.037038683681385,
    0.037038683681385, 0.037038683681385, 0.037038683681385, 0.037038683681385,
    0.037038683681385};
static const double G_TRI_X286[] = {0.04347826086956522,
                                    0.1304347826086956,
                                    0.2173913043478261,
                                    0.3043478260869565,
                                    0.391304347826087,
                                    0.4782608695652174,
                                    0.5652173913043478,
                                    0.6521739130434783,
                                    0.7391304347826086,
                                    0.8260869565217391,
                                    0.9130434782608695,
                                    0.04347826086956522,
                                    0.1304347826086956,
                                    0.2173913043478261,
                                    0.3043478260869565,
                                    0.391304347826087,
                                    0.4782608695652174,
                                    0.5652173913043478,
                                    0.6521739130434783,
                                    0.7391304347826086,
                                    0.8260869565217391,
                                    0.04347826086956522,
                                    0.1304347826086956,
                                    0.2173913043478261,
                                    0.3043478260869565,
                                    0.391304347826087,
                                    0.4782608695652174,
                                    0.5652173913043478,
                                    0.6521739130434783,
                                    0.7391304347826086,
                                    0.04347826086956522,
                                    0.1304347826086956,
                                    0.2173913043478261,
                                    0.3043478260869565,
                                    0.391304347826087,
                                    0.4782608695652174,
                                    0.5652173913043478,
                                    0.6521739130434783,
                                    0.04347826086956522,
                                    0.1304347826086956,
                                    0.2173913043478261,
                                    0.3043478260869565,
                                    0.391304347826087,
                                    0.4782608695652174,
                                    0.5652173913043478,
                                    0.04347826086956522,
                                    0.1304347826086956,
                                    0.2173913043478261,
                                    0.3043478260869565,
                                    0.391304347826087,
                                    0.4782608695652174,
                                    0.04347826086956522,
                                    0.1304347826086956,
                                    0.2173913043478261,
                                    0.3043478260869565,
                                    0.391304347826087,
                                    0.04347826086956522,
                                    0.1304347826086956,
                                    0.2173913043478261,
                                    0.3043478260869565,
                                    0.04347826086956522,
                                    0.1304347826086956,
                                    0.2173913043478261,
                                    0.04347826086956522,
                                    0.1304347826086956,
                                    0.04347826086956522,
                                    0.04761904761904762,
                                    0.1428571428571428,
                                    0.2380952380952381,
                                    0.3333333333333333,
                                    0.4285714285714285,
                                    0.5238095238095238,
                                    0.6190476190476191,
                                    0.7142857142857143,
                                    0.8095238095238095,
                                    0.9047619047619048,
                                    0.04761904761904762,
                                    0.1428571428571428,
                                    0.2380952380952381,
                                    0.3333333333333333,
                                    0.4285714285714285,
                                    0.5238095238095238,
                                    0.6190476190476191,
                                    0.7142857142857143,
                                    0.8095238095238095,
                                    0.04761904761904762,
                                    0.1428571428571428,
                                    0.2380952380952381,
                                    0.3333333333333333,
                                    0.4285714285714285,
                                    0.5238095238095238,
                                    0.6190476190476191,
                                    0.7142857142857143,
                                    0.04761904761904762,
                                    0.1428571428571428,
                                    0.2380952380952381,
                                    0.3333333333333333,
                                    0.4285714285714285,
                                    0.5238095238095238,
                                    0.6190476190476191,
                                    0.04761904761904762,
                                    0.1428571428571428,
                                    0.2380952380952381,
                                    0.3333333333333333,
                                    0.4285714285714285,
                                    0.5238095238095238,
                                    0.04761904761904762,
                                    0.1428571428571428,
                                    0.2380952380952381,
                                    0.3333333333333333,
                                    0.4285714285714285,
                                    0.04761904761904762,
                                    0.1428571428571428,
                                    0.2380952380952381,
                                    0.3333333333333333,
                                    0.04761904761904762,
                                    0.1428571428571428,
                                    0.2380952380952381,
                                    0.04761904761904762,
                                    0.1428571428571428,
                                    0.04761904761904762,
                                    0.05263157894736842,
                                    0.1578947368421053,
                                    0.2631578947368421,
                                    0.3684210526315789,
                                    0.4736842105263158,
                                    0.5789473684210527,
                                    0.6842105263157895,
                                    0.7894736842105263,
                                    0.8947368421052632,
                                    0.05263157894736842,
                                    0.1578947368421053,
                                    0.2631578947368421,
                                    0.3684210526315789,
                                    0.4736842105263158,
                                    0.5789473684210527,
                                    0.6842105263157895,
                                    0.7894736842105263,
                                    0.05263157894736842,
                                    0.1578947368421053,
                                    0.2631578947368421,
                                    0.3684210526315789,
                                    0.4736842105263158,
                                    0.5789473684210527,
                                    0.6842105263157895,
                                    0.05263157894736842,
                                    0.1578947368421053,
                                    0.2631578947368421,
                                    0.3684210526315789,
                                    0.4736842105263158,
                                    0.5789473684210527,
                                    0.05263157894736842,
                                    0.1578947368421053,
                                    0.2631578947368421,
                                    0.3684210526315789,
                                    0.4736842105263158,
                                    0.05263157894736842,
                                    0.1578947368421053,
                                    0.2631578947368421,
                                    0.3684210526315789,
                                    0.05263157894736842,
                                    0.1578947368421053,
                                    0.2631578947368421,
                                    0.05263157894736842,
                                    0.1578947368421053,
                                    0.05263157894736842,
                                    0.05882352941176471,
                                    0.1764705882352941,
                                    0.2941176470588235,
                                    0.4117647058823529,
                                    0.5294117647058824,
                                    0.6470588235294118,
                                    0.7647058823529411,
                                    0.8823529411764706,
                                    0.05882352941176471,
                                    0.1764705882352941,
                                    0.2941176470588235,
                                    0.4117647058823529,
                                    0.5294117647058824,
                                    0.6470588235294118,
                                    0.7647058823529411,
                                    0.05882352941176471,
                                    0.1764705882352941,
                                    0.2941176470588235,
                                    0.4117647058823529,
                                    0.5294117647058824,
                                    0.6470588235294118,
                                    0.05882352941176471,
                                    0.1764705882352941,
                                    0.2941176470588235,
                                    0.4117647058823529,
                                    0.5294117647058824,
                                    0.05882352941176471,
                                    0.1764705882352941,
                                    0.2941176470588235,
                                    0.4117647058823529,
                                    0.05882352941176471,
                                    0.1764705882352941,
                                    0.2941176470588235,
                                    0.05882352941176471,
                                    0.1764705882352941,
                                    0.05882352941176471,
                                    0.06666666666666667,
                                    0.2,
                                    0.3333333333333333,
                                    0.4666666666666667,
                                    0.6,
                                    0.7333333333333333,
                                    0.8666666666666667,
                                    0.06666666666666667,
                                    0.2,
                                    0.3333333333333333,
                                    0.4666666666666667,
                                    0.6,
                                    0.7333333333333333,
                                    0.06666666666666667,
                                    0.2,
                                    0.3333333333333333,
                                    0.4666666666666667,
                                    0.6,
                                    0.06666666666666667,
                                    0.2,
                                    0.3333333333333333,
                                    0.4666666666666667,
                                    0.06666666666666667,
                                    0.2,
                                    0.3333333333333333,
                                    0.06666666666666667,
                                    0.2,
                                    0.06666666666666667,
                                    0.07692307692307693,
                                    0.2307692307692308,
                                    0.3846153846153846,
                                    0.5384615384615384,
                                    0.6923076923076923,
                                    0.8461538461538461,
                                    0.07692307692307693,
                                    0.2307692307692308,
                                    0.3846153846153846,
                                    0.5384615384615384,
                                    0.6923076923076923,
                                    0.07692307692307693,
                                    0.2307692307692308,
                                    0.3846153846153846,
                                    0.5384615384615384,
                                    0.07692307692307693,
                                    0.2307692307692308,
                                    0.3846153846153846,
                                    0.07692307692307693,
                                    0.2307692307692308,
                                    0.07692307692307693,
                                    0.09090909090909091,
                                    0.2727272727272727,
                                    0.4545454545454545,
                                    0.6363636363636364,
                                    0.8181818181818182,
                                    0.09090909090909091,
                                    0.2727272727272727,
                                    0.4545454545454545,
                                    0.6363636363636364,
                                    0.09090909090909091,
                                    0.2727272727272727,
                                    0.4545454545454545,
                                    0.09090909090909091,
                                    0.2727272727272727,
                                    0.09090909090909091,
                                    0.1111111111111111,
                                    0.3333333333333333,
                                    0.5555555555555556,
                                    0.7777777777777778,
                                    0.1111111111111111,
                                    0.3333333333333333,
                                    0.5555555555555556,
                                    0.1111111111111111,
                                    0.3333333333333333,
                                    0.1111111111111111,
                                    0.1428571428571428,
                                    0.4285714285714285,
                                    0.7142857142857143,
                                    0.1428571428571428,
                                    0.4285714285714285,
                                    0.1428571428571428,
                                    0.2,
                                    0.6,
                                    0.2,
                                    0.3333333333333333};
static const double G_TRI_Y286[] = {0.04347826086956522,
                                    0.04347826086956522,
                                    0.04347826086956522,
                                    0.04347826086956522,
                                    0.04347826086956522,
                                    0.04347826086956522,
                                    0.04347826086956522,
                                    0.04347826086956522,
                                    0.04347826086956522,
                                    0.04347826086956522,
                                    0.04347826086956522,
                                    0.1304347826086956,
                                    0.1304347826086956,
                                    0.1304347826086956,
                                    0.1304347826086956,
                                    0.1304347826086956,
                                    0.1304347826086956,
                                    0.1304347826086956,
                                    0.1304347826086956,
                                    0.1304347826086956,
                                    0.1304347826086956,
                                    0.2173913043478261,
                                    0.2173913043478261,
                                    0.2173913043478261,
                                    0.2173913043478261,
                                    0.2173913043478261,
                                    0.2173913043478261,
                                    0.2173913043478261,
                                    0.2173913043478261,
                                    0.2173913043478261,
                                    0.3043478260869565,
                                    0.3043478260869565,
                                    0.3043478260869565,
                                    0.3043478260869565,
                                    0.3043478260869565,
                                    0.3043478260869565,
                                    0.3043478260869565,
                                    0.3043478260869565,
                                    0.391304347826087,
                                    0.391304347826087,
                                    0.391304347826087,
                                    0.391304347826087,
                                    0.391304347826087,
                                    0.391304347826087,
                                    0.391304347826087,
                                    0.4782608695652174,
                                    0.4782608695652174,
                                    0.4782608695652174,
                                    0.4782608695652174,
                                    0.4782608695652174,
                                    0.4782608695652174,
                                    0.5652173913043478,
                                    0.5652173913043478,
                                    0.5652173913043478,
                                    0.5652173913043478,
                                    0.5652173913043478,
                                    0.6521739130434783,
                                    0.6521739130434783,
                                    0.6521739130434783,
                                    0.6521739130434783,
                                    0.7391304347826086,
                                    0.7391304347826086,
                                    0.7391304347826086,
                                    0.8260869565217391,
                                    0.8260869565217391,
                                    0.9130434782608695,
                                    0.04761904761904762,
                                    0.04761904761904762,
                                    0.04761904761904762,
                                    0.04761904761904762,
                                    0.04761904761904762,
                                    0.04761904761904762,
                                    0.04761904761904762,
                                    0.04761904761904762,
                                    0.04761904761904762,
                                    0.04761904761904762,
                                    0.1428571428571428,
                                    0.1428571428571428,
                                    0.1428571428571428,
                                    0.1428571428571428,
                                    0.1428571428571428,
                                    0.1428571428571428,
                                    0.1428571428571428,
                                    0.1428571428571428,
                                    0.1428571428571428,
                                    0.2380952380952381,
                                    0.2380952380952381,
                                    0.2380952380952381,
                                    0.2380952380952381,
                                    0.2380952380952381,
                                    0.2380952380952381,
                                    0.2380952380952381,
                                    0.2380952380952381,
                                    0.3333333333333333,
                                    0.3333333333333333,
                                    0.3333333333333333,
                                    0.3333333333333333,
                                    0.3333333333333333,
                                    0.3333333333333333,
                                    0.3333333333333333,
                                    0.4285714285714285,
                                    0.4285714285714285,
                                    0.4285714285714285,
                                    0.4285714285714285,
                                    0.4285714285714285,
                                    0.4285714285714285,
                                    0.5238095238095238,
                                    0.5238095238095238,
                                    0.5238095238095238,
                                    0.5238095238095238,
                                    0.5238095238095238,
                                    0.6190476190476191,
                                    0.6190476190476191,
                                    0.6190476190476191,
                                    0.6190476190476191,
                                    0.7142857142857143,
                                    0.7142857142857143,
                                    0.7142857142857143,
                                    0.8095238095238095,
                                    0.8095238095238095,
                                    0.9047619047619048,
                                    0.05263157894736842,
                                    0.05263157894736842,
                                    0.05263157894736842,
                                    0.05263157894736842,
                                    0.05263157894736842,
                                    0.05263157894736842,
                                    0.05263157894736842,
                                    0.05263157894736842,
                                    0.05263157894736842,
                                    0.1578947368421053,
                                    0.1578947368421053,
                                    0.1578947368421053,
                                    0.1578947368421053,
                                    0.1578947368421053,
                                    0.1578947368421053,
                                    0.1578947368421053,
                                    0.1578947368421053,
                                    0.2631578947368421,
                                    0.2631578947368421,
                                    0.2631578947368421,
                                    0.2631578947368421,
                                    0.2631578947368421,
                                    0.2631578947368421,
                                    0.2631578947368421,
                                    0.3684210526315789,
                                    0.3684210526315789,
                                    0.3684210526315789,
                                    0.3684210526315789,
                                    0.3684210526315789,
                                    0.3684210526315789,
                                    0.4736842105263158,
                                    0.4736842105263158,
                                    0.4736842105263158,
                                    0.4736842105263158,
                                    0.4736842105263158,
                                    0.5789473684210527,
                                    0.5789473684210527,
                                    0.5789473684210527,
                                    0.5789473684210527,
                                    0.6842105263157895,
                                    0.6842105263157895,
                                    0.6842105263157895,
                                    0.7894736842105263,
                                    0.7894736842105263,
                                    0.8947368421052632,
                                    0.05882352941176471,
                                    0.05882352941176471,
                                    0.05882352941176471,
                                    0.05882352941176471,
                                    0.05882352941176471,
                                    0.05882352941176471,
                                    0.05882352941176471,
                                    0.05882352941176471,
                                    0.1764705882352941,
                                    0.1764705882352941,
                                    0.1764705882352941,
                                    0.1764705882352941,
                                    0.1764705882352941,
                                    0.1764705882352941,
                                    0.1764705882352941,
                                    0.2941176470588235,
                                    0.2941176470588235,
                                    0.2941176470588235,
                                    0.2941176470588235,
                                    0.2941176470588235,
                                    0.2941176470588235,
                                    0.4117647058823529,
                                    0.4117647058823529,
                                    0.4117647058823529,
                                    0.4117647058823529,
                                    0.4117647058823529,
                                    0.5294117647058824,
                                    0.5294117647058824,
                                    0.5294117647058824,
                                    0.5294117647058824,
                                    0.6470588235294118,
                                    0.6470588235294118,
                                    0.6470588235294118,
                                    0.7647058823529411,
                                    0.7647058823529411,
                                    0.8823529411764706,
                                    0.06666666666666667,
                                    0.06666666666666667,
                                    0.06666666666666667,
                                    0.06666666666666667,
                                    0.06666666666666667,
                                    0.06666666666666667,
                                    0.06666666666666667,
                                    0.2,
                                    0.2,
                                    0.2,
                                    0.2,
                                    0.2,
                                    0.2,
                                    0.3333333333333333,
                                    0.3333333333333333,
                                    0.3333333333333333,
                                    0.3333333333333333,
                                    0.3333333333333333,
                                    0.4666666666666667,
                                    0.4666666666666667,
                                    0.4666666666666667,
                                    0.4666666666666667,
                                    0.6,
                                    0.6,
                                    0.6,
                                    0.7333333333333333,
                                    0.7333333333333333,
                                    0.8666666666666667,
                                    0.07692307692307693,
                                    0.07692307692307693,
                                    0.07692307692307693,
                                    0.07692307692307693,
                                    0.07692307692307693,
                                    0.07692307692307693,
                                    0.2307692307692308,
                                    0.2307692307692308,
                                    0.2307692307692308,
                                    0.2307692307692308,
                                    0.2307692307692308,
                                    0.3846153846153846,
                                    0.3846153846153846,
                                    0.3846153846153846,
                                    0.3846153846153846,
                                    0.5384615384615384,
                                    0.5384615384615384,
                                    0.5384615384615384,
                                    0.6923076923076923,
                                    0.6923076923076923,
                                    0.8461538461538461,
                                    0.09090909090909091,
                                    0.09090909090909091,
                                    0.09090909090909091,
                                    0.09090909090909091,
                                    0.09090909090909091,
                                    0.2727272727272727,
                                    0.2727272727272727,
                                    0.2727272727272727,
                                    0.2727272727272727,
                                    0.4545454545454545,
                                    0.4545454545454545,
                                    0.4545454545454545,
                                    0.6363636363636364,
                                    0.6363636363636364,
                                    0.8181818181818182,
                                    0.1111111111111111,
                                    0.1111111111111111,
                                    0.1111111111111111,
                                    0.1111111111111111,
                                    0.3333333333333333,
                                    0.3333333333333333,
                                    0.3333333333333333,
                                    0.5555555555555556,
                                    0.5555555555555556,
                                    0.7777777777777778,
                                    0.1428571428571428,
                                    0.1428571428571428,
                                    0.1428571428571428,
                                    0.4285714285714285,
                                    0.4285714285714285,
                                    0.7142857142857143,
                                    0.2,
                                    0.2,
                                    0.6,
                                    0.3333333333333333};
static const double G_TRI_W286[] = {
    2.912193380035668,      2.912193380035668,      2.912193380035668,
    2.912193380035668,      2.912193380035668,      2.912193380035668,
    2.912193380035668,      2.912193380035668,      2.912193380035668,
    2.912193380035668,      2.912193380035668,      2.912193380035668,
    2.912193380035668,      2.912193380035668,      2.912193380035668,
    2.912193380035668,      2.912193380035668,      2.912193380035668,
    2.912193380035668,      2.912193380035668,      2.912193380035668,
    2.912193380035668,      2.912193380035668,      2.912193380035668,
    2.912193380035668,      2.912193380035668,      2.912193380035668,
    2.912193380035668,      2.912193380035668,      2.912193380035668,
    2.912193380035668,      2.912193380035668,      2.912193380035668,
    2.912193380035668,      2.912193380035668,      2.912193380035668,
    2.912193380035668,      2.912193380035668,      2.912193380035668,
    2.912193380035668,      2.912193380035668,      2.912193380035668,
    2.912193380035668,      2.912193380035668,      2.912193380035668,
    2.912193380035668,      2.912193380035668,      2.912193380035668,
    2.912193380035668,      2.912193380035668,      2.912193380035668,
    2.912193380035668,      2.912193380035668,      2.912193380035668,
    2.912193380035668,      2.912193380035668,      2.912193380035668,
    2.912193380035668,      2.912193380035668,      2.912193380035668,
    2.912193380035668,      2.912193380035668,      2.912193380035668,
    2.912193380035668,      2.912193380035668,      2.912193380035668,
    -9.914451197589852,     -9.914451197589852,     -9.914451197589852,
    -9.914451197589852,     -9.914451197589852,     -9.914451197589852,
    -9.914451197589852,     -9.914451197589852,     -9.914451197589852,
    -9.914451197589852,     -9.914451197589852,     -9.914451197589852,
    -9.914451197589852,     -9.914451197589852,     -9.914451197589852,
    -9.914451197589852,     -9.914451197589852,     -9.914451197589852,
    -9.914451197589852,     -9.914451197589852,     -9.914451197589852,
    -9.914451197589852,     -9.914451197589852,     -9.914451197589852,
    -9.914451197589852,     -9.914451197589852,     -9.914451197589852,
    -9.914451197589852,     -9.914451197589852,     -9.914451197589852,
    -9.914451197589852,     -9.914451197589852,     -9.914451197589852,
    -9.914451197589852,     -9.914451197589852,     -9.914451197589852,
    -9.914451197589852,     -9.914451197589852,     -9.914451197589852,
    -9.914451197589852,     -9.914451197589852,     -9.914451197589852,
    -9.914451197589852,     -9.914451197589852,     -9.914451197589852,
    -9.914451197589852,     -9.914451197589852,     -9.914451197589852,
    -9.914451197589852,     -9.914451197589852,     -9.914451197589852,
    -9.914451197589852,     -9.914451197589852,     -9.914451197589852,
    -9.914451197589852,     13.33158527957992,      13.33158527957992,
    13.33158527957992,      13.33158527957992,      13.33158527957992,
    13.33158527957992,      13.33158527957992,      13.33158527957992,
    13.33158527957992,      13.33158527957992,      13.33158527957992,
    13.33158527957992,      13.33158527957992,      13.33158527957992,
    13.33158527957992,      13.33158527957992,      13.33158527957992,
    13.33158527957992,      13.33158527957992,      13.33158527957992,
    13.33158527957992,      13.33158527957992,      13.33158527957992,
    13.33158527957992,      13.33158527957992,      13.33158527957992,
    13.33158527957992,      13.33158527957992,      13.33158527957992,
    13.33158527957992,      13.33158527957992,      13.33158527957992,
    13.33158527957992,      13.33158527957992,      13.33158527957992,
    13.33158527957992,      13.33158527957992,      13.33158527957992,
    13.33158527957992,      13.33158527957992,      13.33158527957992,
    13.33158527957992,      13.33158527957992,      13.33158527957992,
    13.33158527957992,      -9.027792408986382,     -9.027792408986382,
    -9.027792408986382,     -9.027792408986382,     -9.027792408986382,
    -9.027792408986382,     -9.027792408986382,     -9.027792408986382,
    -9.027792408986382,     -9.027792408986382,     -9.027792408986382,
    -9.027792408986382,     -9.027792408986382,     -9.027792408986382,
    -9.027792408986382,     -9.027792408986382,     -9.027792408986382,
    -9.027792408986382,     -9.027792408986382,     -9.027792408986382,
    -9.027792408986382,     -9.027792408986382,     -9.027792408986382,
    -9.027792408986382,     -9.027792408986382,     -9.027792408986382,
    -9.027792408986382,     -9.027792408986382,     -9.027792408986382,
    -9.027792408986382,     -9.027792408986382,     -9.027792408986382,
    -9.027792408986382,     -9.027792408986382,     -9.027792408986382,
    -9.027792408986382,     3.258672079964582,      3.258672079964582,
    3.258672079964582,      3.258672079964582,      3.258672079964582,
    3.258672079964582,      3.258672079964582,      3.258672079964582,
    3.258672079964582,      3.258672079964582,      3.258672079964582,
    3.258672079964582,      3.258672079964582,      3.258672079964582,
    3.258672079964582,      3.258672079964582,      3.258672079964582,
    3.258672079964582,      3.258672079964582,      3.258672079964582,
    3.258672079964582,      3.258672079964582,      3.258672079964582,
    3.258672079964582,      3.258672079964582,      3.258672079964582,
    3.258672079964582,      3.258672079964582,      -0.6133639040302452,
    -0.6133639040302452,    -0.6133639040302452,    -0.6133639040302452,
    -0.6133639040302452,    -0.6133639040302452,    -0.6133639040302452,
    -0.6133639040302452,    -0.6133639040302452,    -0.6133639040302452,
    -0.6133639040302452,    -0.6133639040302452,    -0.6133639040302452,
    -0.6133639040302452,    -0.6133639040302452,    -0.6133639040302452,
    -0.6133639040302452,    -0.6133639040302452,    -0.6133639040302452,
    -0.6133639040302452,    -0.6133639040302452,    0.05511571669513555,
    0.05511571669513555,    0.05511571669513555,    0.05511571669513555,
    0.05511571669513555,    0.05511571669513555,    0.05511571669513555,
    0.05511571669513555,    0.05511571669513555,    0.05511571669513555,
    0.05511571669513555,    0.05511571669513555,    0.05511571669513555,
    0.05511571669513555,    0.05511571669513555,    -0.001979122382447095,
    -0.001979122382447095,  -0.001979122382447095,  -0.001979122382447095,
    -0.001979122382447095,  -0.001979122382447095,  -0.001979122382447095,
    -0.001979122382447095,  -0.001979122382447095,  -0.001979122382447095,
    2.02054621415273e-05,   2.02054621415273e-05,   2.02054621415273e-05,
    2.02054621415273e-05,   2.02054621415273e-05,   2.02054621415273e-05,
    -2.874940020535803e-08, -2.874940020535803e-08, -2.874940020535803e-08,
    8.829438425435718e-13};

// TET
static const double G_TET_X1[] = {0.25};
static const double G_TET_Y1[] = {0.25};
static const double G_TET_Z1[] = {0.25};
static const double G_TET_W1[] = {1.};
static const double G_TET_X4[] = {0.1757281246520584, 0.2445310270213291,
                                  0.5556470949048655, 0.0240937534217468};
static const double G_TET_Y4[] = {0.5656137776620919, 0.0501800797762026,
                                  0.1487681308666864, 0.2354380116950194};
static const double G_TET_Z4[] = {0.2180665126782654, 0.5635595064952189,
                                  0.0350112499848832, 0.1833627308416330};
static const double G_TET_W4[] = {0.25, 0.25, 0.25, 0.25};
static const double G_TET_X5[] = {0.25000000000000000, 0.50000000000000000,
                                  0.16666666666666667, 0.16666666666666667,
                                  0.16666666666666667};
static const double G_TET_Y5[] = {0.25000000000000000, 0.16666666666666667,
                                  0.50000000000000000, 0.16666666666666667,
                                  0.16666666666666667};
static const double G_TET_Z5[] = {0.25000000000000000, 0.16666666666666667,
                                  0.16666666666666667, 0.50000000000000000,
                                  0.16666666666666667};
static const double G_TET_W5[] = {-0.80000000000000000, 0.45000000000000000,
                                  0.45000000000000000, 0.45000000000000000,
                                  0.45000000000000000};
static const double G_TET_X10[] = {0.5684305841968444, 0.1438564719343852,
                                   0.1438564719343852, 0.1438564719343852,
                                   0.0000000000000000, 0.5000000000000000,
                                   0.5000000000000000, 0.5000000000000000,
                                   0.0000000000000000, 0.0000000000000000};
static const double G_TET_Y10[] = {0.1438564719343852, 0.1438564719343852,
                                   0.1438564719343852, 0.5684305841968444,
                                   0.5000000000000000, 0.0000000000000000,
                                   0.5000000000000000, 0.0000000000000000,
                                   0.5000000000000000, 0.0000000000000000};
static const double G_TET_Z10[] = {0.1438564719343852, 0.1438564719343852,
                                   0.5684305841968444, 0.1438564719343852,
                                   0.5000000000000000, 0.5000000000000000,
                                   0.0000000000000000, 0.0000000000000000,
                                   0.0000000000000000, 0.5000000000000000};
static const double G_TET_W10[] = {0.2177650698804054, 0.2177650698804054,
                                   0.2177650698804054, 0.2177650698804054,
                                   0.0214899534130631, 0.0214899534130631,
                                   0.0214899534130631, 0.0214899534130631,
                                   0.0214899534130631, 0.0214899534130631};

static const double G_TET_X45[] = {
    0.2500000000000000, 0.6175871903000830, 0.1274709365666390,
    0.1274709365666390, 0.1274709365666390, 0.9037635088221031,
    0.0320788303926323, 0.0320788303926323, 0.0320788303926323,
    0.4502229043567190, 0.0497770956432810, 0.0497770956432810,
    0.0497770956432810, 0.4502229043567190, 0.4502229043567190,
    0.3162695526014501, 0.1837304473985499, 0.1837304473985499,
    0.1837304473985499, 0.3162695526014501, 0.3162695526014501,
    0.0229177878448171, 0.2319010893971509, 0.2319010893971509,
    0.5132800333608811, 0.2319010893971509, 0.2319010893971509,
    0.2319010893971509, 0.0229177878448171, 0.5132800333608811,
    0.2319010893971509, 0.0229177878448171, 0.5132800333608811,
    0.7303134278075384, 0.0379700484718286, 0.0379700484718286,
    0.1937464752488044, 0.0379700484718286, 0.0379700484718286,
    0.0379700484718286, 0.7303134278075384, 0.1937464752488044,
    0.0379700484718286, 0.7303134278075384, 0.1937464752488044};
static const double G_TET_Y45[] = {
    0.2500000000000000, 0.1274709365666390, 0.1274709365666390,
    0.1274709365666390, 0.6175871903000830, 0.0320788303926323,
    0.0320788303926323, 0.0320788303926323, 0.9037635088221031,
    0.0497770956432810, 0.4502229043567190, 0.0497770956432810,
    0.4502229043567190, 0.0497770956432810, 0.4502229043567190,
    0.1837304473985499, 0.3162695526014501, 0.1837304473985499,
    0.3162695526014501, 0.1837304473985499, 0.3162695526014501,
    0.2319010893971509, 0.0229177878448171, 0.2319010893971509,
    0.2319010893971509, 0.5132800333608811, 0.2319010893971509,
    0.0229177878448171, 0.5132800333608811, 0.2319010893971509,
    0.5132800333608811, 0.2319010893971509, 0.0229177878448171,
    0.0379700484718286, 0.7303134278075384, 0.0379700484718286,
    0.0379700484718286, 0.1937464752488044, 0.0379700484718286,
    0.7303134278075384, 0.1937464752488044, 0.0379700484718286,
    0.1937464752488044, 0.0379700484718286, 0.7303134278075384};
static const double G_TET_Z45[] = {
    0.2500000000000000, 0.1274709365666390, 0.1274709365666390,
    0.6175871903000830, 0.1274709365666390, 0.0320788303926323,
    0.0320788303926323, 0.9037635088221031, 0.0320788303926323,
    0.0497770956432810, 0.0497770956432810, 0.4502229043567190,
    0.4502229043567190, 0.4502229043567190, 0.0497770956432810,
    0.1837304473985499, 0.1837304473985499, 0.3162695526014501,
    0.3162695526014501, 0.3162695526014501, 0.1837304473985499,
    0.2319010893971509, 0.2319010893971509, 0.0229177878448171,
    0.2319010893971509, 0.2319010893971509, 0.5132800333608811,
    0.5132800333608811, 0.2319010893971509, 0.0229177878448171,
    0.0229177878448171, 0.5132800333608811, 0.2319010893971509,
    0.0379700484718286, 0.0379700484718286, 0.7303134278075384,
    0.0379700484718286, 0.0379700484718286, 0.1937464752488044,
    0.1937464752488044, 0.0379700484718286, 0.7303134278075384,
    0.7303134278075384, 0.1937464752488044, 0.0379700484718286};
static const double G_TET_W45[] = {
    -0.2359620398477557, 0.0244878963560562, 0.0244878963560562,
    0.0244878963560562,  0.0244878963560562, 0.0039485206398261,
    0.0039485206398261,  0.0039485206398261, 0.0039485206398261,
    0.0263055529507371,  0.0263055529507371, 0.0263055529507371,
    0.0263055529507371,  0.0263055529507371, 0.0263055529507371,
    0.0829803830550589,  0.0829803830550589, 0.0829803830550589,
    0.0829803830550589,  0.0829803830550589, 0.0829803830550589,
    0.0254426245481023,  0.0254426245481023, 0.0254426245481023,
    0.0254426245481023,  0.0254426245481023, 0.0254426245481023,
    0.0254426245481023,  0.0254426245481023, 0.0254426245481023,
    0.0254426245481023,  0.0254426245481023, 0.0254426245481023,
    0.0134324384376852,  0.0134324384376852, 0.0134324384376852,
    0.0134324384376852,  0.0134324384376852, 0.0134324384376852,
    0.0134324384376852,  0.0134324384376852, 0.0134324384376852,
    0.0134324384376852,  0.0134324384376852, 0.0134324384376852};
static const double NC_TET_X84[] = {
    0.1000000000000000, 0.1000000000000000, 0.1000000000000000,
    0.7000000000000000, 0.1000000000000000, 0.1000000000000000,
    0.1000000000000000, 0.1000000000000000, 0.1000000000000000,
    0.1000000000000000, 0.6000000000000000, 0.6000000000000000,
    0.6000000000000000, 0.2000000000000000, 0.2000000000000000,
    0.2000000000000000, 0.1000000000000000, 0.1000000000000000,
    0.1000000000000000, 0.1000000000000000, 0.1000000000000000,
    0.1000000000000000, 0.5000000000000000, 0.5000000000000000,
    0.5000000000000000, 0.3000000000000000, 0.3000000000000000,
    0.3000000000000000, 0.2000000000000000, 0.2000000000000000,
    0.2000000000000000, 0.2000000000000000, 0.2000000000000000,
    0.2000000000000000, 0.1000000000000000, 0.1000000000000000,
    0.1000000000000000, 0.5000000000000000, 0.5000000000000000,
    0.5000000000000000, 0.4000000000000000, 0.4000000000000000,
    0.4000000000000000, 0.1000000000000000, 0.1000000000000000,
    0.1000000000000000, 0.4000000000000000, 0.4000000000000000,
    0.4000000000000000, 0.4000000000000000, 0.4000000000000000,
    0.4000000000000000, 0.3000000000000000, 0.3000000000000000,
    0.3000000000000000, 0.3000000000000000, 0.3000000000000000,
    0.3000000000000000, 0.2000000000000000, 0.2000000000000000,
    0.2000000000000000, 0.2000000000000000, 0.2000000000000000,
    0.2000000000000000, 0.1000000000000000, 0.1000000000000000,
    0.1000000000000000, 0.1000000000000000, 0.1000000000000000,
    0.1000000000000000, 0.2000000000000000, 0.2000000000000000,
    0.2000000000000000, 0.4000000000000000, 0.3000000000000000,
    0.3000000000000000, 0.3000000000000000, 0.1000000000000000,
    0.3000000000000000, 0.3000000000000000, 0.3000000000000000,
    0.2000000000000000, 0.2000000000000000, 0.2000000000000000};
static const double NC_TET_Y84[] = {
    0.1000000000000000, 0.1000000000000000, 0.7000000000000000,
    0.1000000000000000, 0.1000000000000000, 0.1000000000000000,
    0.6000000000000000, 0.6000000000000000, 0.2000000000000000,
    0.2000000000000000, 0.1000000000000000, 0.1000000000000000,
    0.2000000000000000, 0.1000000000000000, 0.1000000000000000,
    0.6000000000000000, 0.1000000000000000, 0.1000000000000000,
    0.5000000000000000, 0.5000000000000000, 0.3000000000000000,
    0.3000000000000000, 0.1000000000000000, 0.1000000000000000,
    0.3000000000000000, 0.1000000000000000, 0.1000000000000000,
    0.5000000000000000, 0.2000000000000000, 0.2000000000000000,
    0.1000000000000000, 0.1000000000000000, 0.5000000000000000,
    0.5000000000000000, 0.2000000000000000, 0.2000000000000000,
    0.5000000000000000, 0.2000000000000000, 0.2000000000000000,
    0.1000000000000000, 0.4000000000000000, 0.1000000000000000,
    0.1000000000000000, 0.4000000000000000, 0.4000000000000000,
    0.1000000000000000, 0.3000000000000000, 0.3000000000000000,
    0.2000000000000000, 0.2000000000000000, 0.1000000000000000,
    0.1000000000000000, 0.4000000000000000, 0.4000000000000000,
    0.2000000000000000, 0.2000000000000000, 0.1000000000000000,
    0.1000000000000000, 0.4000000000000000, 0.4000000000000000,
    0.3000000000000000, 0.3000000000000000, 0.1000000000000000,
    0.1000000000000000, 0.4000000000000000, 0.4000000000000000,
    0.3000000000000000, 0.3000000000000000, 0.2000000000000000,
    0.2000000000000000, 0.2000000000000000, 0.2000000000000000,
    0.4000000000000000, 0.2000000000000000, 0.3000000000000000,
    0.3000000000000000, 0.1000000000000000, 0.3000000000000000,
    0.3000000000000000, 0.2000000000000000, 0.2000000000000000,
    0.3000000000000000, 0.3000000000000000, 0.2000000000000000};
static const double NC_TET_Z84[] = {
    0.1000000000000000, 0.7000000000000000, 0.1000000000000000,
    0.1000000000000000, 0.6000000000000000, 0.2000000000000000,
    0.1000000000000000, 0.2000000000000000, 0.1000000000000000,
    0.6000000000000000, 0.1000000000000000, 0.2000000000000000,
    0.1000000000000000, 0.1000000000000000, 0.6000000000000000,
    0.1000000000000000, 0.5000000000000000, 0.3000000000000000,
    0.1000000000000000, 0.3000000000000000, 0.1000000000000000,
    0.5000000000000000, 0.1000000000000000, 0.3000000000000000,
    0.1000000000000000, 0.1000000000000000, 0.5000000000000000,
    0.1000000000000000, 0.1000000000000000, 0.5000000000000000,
    0.2000000000000000, 0.5000000000000000, 0.2000000000000000,
    0.1000000000000000, 0.2000000000000000, 0.5000000000000000,
    0.2000000000000000, 0.2000000000000000, 0.1000000000000000,
    0.2000000000000000, 0.1000000000000000, 0.4000000000000000,
    0.1000000000000000, 0.4000000000000000, 0.1000000000000000,
    0.4000000000000000, 0.2000000000000000, 0.1000000000000000,
    0.3000000000000000, 0.1000000000000000, 0.3000000000000000,
    0.2000000000000000, 0.2000000000000000, 0.1000000000000000,
    0.4000000000000000, 0.1000000000000000, 0.4000000000000000,
    0.2000000000000000, 0.3000000000000000, 0.1000000000000000,
    0.4000000000000000, 0.1000000000000000, 0.4000000000000000,
    0.3000000000000000, 0.3000000000000000, 0.2000000000000000,
    0.4000000000000000, 0.2000000000000000, 0.4000000000000000,
    0.3000000000000000, 0.2000000000000000, 0.4000000000000000,
    0.2000000000000000, 0.2000000000000000, 0.3000000000000000,
    0.1000000000000000, 0.3000000000000000, 0.3000000000000000,
    0.2000000000000000, 0.3000000000000000, 0.2000000000000000,
    0.3000000000000000, 0.2000000000000000, 0.3000000000000000};
static const double NC_TET_W84[] = {
    0.2843915343915344,  0.2843915343915344,  0.2843915343915344,
    0.2843915343915344,  -0.3882275132275133, -0.3882275132275133,
    -0.3882275132275133, -0.3882275132275133, -0.3882275132275133,
    -0.3882275132275133, -0.3882275132275133, -0.3882275132275133,
    -0.3882275132275133, -0.3882275132275133, -0.3882275132275133,
    -0.3882275132275133, 0.8776455026455027,  0.8776455026455027,
    0.8776455026455027,  0.8776455026455027,  0.8776455026455027,
    0.8776455026455027,  0.8776455026455027,  0.8776455026455027,
    0.8776455026455027,  0.8776455026455027,  0.8776455026455027,
    0.8776455026455027,  0.1236772486772487,  0.1236772486772487,
    0.1236772486772487,  0.1236772486772487,  0.1236772486772487,
    0.1236772486772487,  0.1236772486772487,  0.1236772486772487,
    0.1236772486772487,  0.1236772486772487,  0.1236772486772487,
    0.1236772486772487,  -0.8584656084656085, -0.8584656084656085,
    -0.8584656084656085, -0.8584656084656085, -0.8584656084656085,
    -0.8584656084656085, -0.2632275132275133, -0.2632275132275133,
    -0.2632275132275133, -0.2632275132275133, -0.2632275132275133,
    -0.2632275132275133, -0.2632275132275133, -0.2632275132275133,
    -0.2632275132275133, -0.2632275132275133, -0.2632275132275133,
    -0.2632275132275133, -0.2632275132275133, -0.2632275132275133,
    -0.2632275132275133, -0.2632275132275133, -0.2632275132275133,
    -0.2632275132275133, -0.2632275132275133, -0.2632275132275133,
    -0.2632275132275133, -0.2632275132275133, -0.2632275132275133,
    -0.2632275132275133, 0.0145502645502645,  0.0145502645502645,
    0.0145502645502645,  0.0145502645502645,  1.0165343915343916,
    1.0165343915343916,  1.0165343915343916,  1.0165343915343916,
    -0.0251322751322751, -0.0251322751322751, -0.0251322751322751,
    -0.0251322751322751, -0.0251322751322751, -0.0251322751322751};

#ifdef __cplusplus
}
#endif

#endif
