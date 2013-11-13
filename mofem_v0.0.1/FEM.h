/** \file FEM.h
 * \brief Core FieldInterface class for user interface
 * 
 * Copyright (C) 2013, Lukasz Kaczmarczyk (likask AT wp.pl) <br>
 * MoFEM is free software: you can redistribute it and/or modify it under
 * the terms of the GNU Lesser General Public License as published by the
 * Free Software Foundation, either version 3 of the License, or (at your
 * option) any later version.
 *
 * MoFEM is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
 * License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with MoFEM. If not, see <http://www.gnu.org/licenses/>
*/

#ifndef __FEM_H__
#define __FEM_H__

#include<stdlib.h>

#include<petscmat.h>
#include<petscao.h>
#include<petscbt.h>

#include<petscmat.h>
#include<petscao.h>
#include<petscbt.h>
#include<petsclog.h>


#ifdef __APPLE__
  #ifdef __CBLAS__
  #include<cblas.h>
  #else  
  #include <Accelerate/Accelerate.h>
  #endif
  #include<lapack_wrap.h>
#else 
  #include<cblas.h>
  #include<lapack_wrap.h>
#endif

#define LAMBDA(E,NU) (E*NU/((1.+NU)*(1.-2.*NU)))
#define MU(E,NU) (0.5*E/(1.+NU))
#define DELTA(NU_P,NU_PZ,E_P,E_Z) (((1+NU_P)*(1-NU_P-2*NU_PZ*(NU_PZ*E_Z/E_P)))/(E_P*E_P*E_Z))

#ifdef __cplusplus
extern "C" {
#endif

/// print matric M
void print_mat(double *M,int m,int n);
/// print upper part of the symmetric matrix
void print_mat_sym_upper(double *M,int m,int n);
/// priint complex matrix
void print_mat_complex(__CLPK_doublecomplex *M,int m,int n);

/// \brief calulate shape functions of trianangle
/// \param N shape function array
/// \param X array of Guass X coordinates 
/// \param Y array of Guass Y coordinates 
/// \param G_DIM number of Gauss points 
PetscErrorCode ShapeMBTRI(double *N,const double *X,const double *Y,const int G_DIM);
/// calulate direvatives of shape functions
PetscErrorCode ShapeDiffMBTRI(double *diffN);

/// calulate face nornam
/// \param diffN direvatives of shape functions
/// \param coords is position of the nodes
/// \param normal vector
PetscErrorCode ShapeFaceNormalMBTRI(double *diffN,const double *coords,double *normal);
PetscErrorCode ShapeFaceBaseMBTRI(
  double *diffN,const double *coords,
  double *normal,double *s1,double *s2);
/// calulate jacobioan 
void ShapeJacMBTRI(double *diffN,const double *coords,double *Jac);
/// calulate direvatives of shape functions in space
void ShapeDiffMBTRIinvJ(double *diffN,double *invJac,double *diffNinvJac);
/// caluate shape functions
PetscErrorCode ShapeMBTET(double *N,const double *G_X,const double *G_Y,const double *G_Z,int DIM);
/// calulare direvatives of shape functions
void ShapeDiffMBTET(double *diffN);
/// determinad of jacobian
double Shape_detJac(double *Jac);
/// calulate jacobian
PetscErrorCode ShapeJacMBTET(double *diffN,const double *coords,double *Jac);
// calulate inverse of jacobian
PetscErrorCode Shape_invJac(double *Jac);
/// calulate TET volume
double Shape_intVolumeMBTET(double *diffN,const double *coords);
/// calulate shape functions direvatives in space
PetscErrorCode ShapeDiffMBTETinvJ(double *diffN,double *invJac,double *diffNinvJac);

/// calulate spin matrix from vector
// \param spinOmega is a spin matrxi
// \param vecOmega is a spin vector
PetscErrorCode Spin(double *spinOmega,double *vecOmega);

/// Compose complex matrix (3x3) from two real matrices
PetscErrorCode make_complex_matrix(double *reA,double *imA,__CLPK_doublecomplex *xA);

/** 
 * \brief calculate local coordinates for given global coordinates
 *
 * new verison for multiple points need to be implemented
 */
PetscErrorCode ShapeMBTET_inverse(double *N,double *diffN,const double *elem_coords,const double *glob_coords,double *loc_coords);
/// calculate gradient of deformation
PetscErrorCode GradientOfDeformation(double *diffN,double *dofs,double *F);

/** 
 * \brief Calulate Lagrange approximation basis
 *
 * \param p is approximation order
 * \param s is is position [-1,1]
 * \param diff_s direvatives of shape functions
 * \param L appeoximation functions
 * \param diffL direvatives
 * \param dim dimension
 */
PetscErrorCode Lagrange_basis(int p,double s,double *diff_s,double *L,double *diffL,const int dim);

//complex part
void ShapeDiffMBTETinvJ_complex(double *diffN,__CLPK_doublecomplex *invJac,__CLPK_doublecomplex *diffNinvJac,const enum CBLAS_TRANSPOSE Trans);
PetscErrorCode ShapeFaceNormalMBTRI_complex(double *diffN,__CLPK_doublecomplex *xcoords,__CLPK_doublecomplex *xnormal);
PetscErrorCode MakeComplexTensor(double *reA,double *imA,__CLPK_doublecomplex *xA);
PetscErrorCode InvertComplexGradient(__CLPK_doublecomplex *xF);
PetscErrorCode DeterminantComplexGradient(__CLPK_doublecomplex *xF,__CLPK_doublecomplex *det_xF);

//http://people.sc.fsu.edu/~jburkardt/datasets/quadrature_rules_tet/quadrature_rules_tet.html
//TRI
static const double G_TRI_X1[] = {
3.3333333333333331e-01 
};
static const double G_TRI_Y1[] = {
3.3333333333333331e-01 
};
static const double G_TRI_W1[] = {
1. 
};
static const double G_TRI_X3[] = {
0.5, 0., 0.5 
};
static const double G_TRI_Y3[] = {
0., 0.5, 0.5 
};
static const double G_TRI_W3[] = {
3.3333333333333331e-01, 3.3333333333333331e-01, 3.3333333333333331e-01 
};
static const double G_TRI_X4[] = {
7.503111022260811058e-02, 1.785587282636164064e-01,
2.800199154990741235e-01, 6.663902460147014262e-01 
};
static const double G_TRI_Y4[] = {
2.800199154990741235e-01, 6.663902460147014262e-01,
7.503111022260811058e-02, 1.785587282636164064e-01 
};
static const double G_TRI_W4[] = {
1.8195861825602258066e-01, 3.1804138174397683647e-01,
1.8195861825602258066e-01, 3.1804138174397683647e-01
};
static const double G_TRI_X7[] = {
0.333333333333333, 0.736712498968435, 0.736712498968435,
0.237932366472434, 0.237932366472434, 0.025355134551932,
0.025355134551932
};
static const double G_TRI_Y7[] = {
0.333333333333333, 0.237932366472434, 0.025355134551932,
0.736712498968435, 0.025355134551932, 0.736712498968435,
0.237932366472434
};
static const double G_TRI_W7[] = {
  0.375000000000000, 0.104166666666667, 0.104166666666667, 0.104166666666667, 0.104166666666667, 0.104166666666667, 0.104166666666667
};
static const double G_TRI_X13[] = {
0.333333333333333, 0.479308067841923, 0.260345966079038, 0.260345966079038, 0.869739794195568,
0.065130102902216, 0.065130102902216, 0.638444188569809, 0.638444188569809, 0.312865496004875,
0.312865496004875, 0.048690315425316, 0.048690315425316
};
static const double G_TRI_Y13[] = {
0.333333333333333, 0.260345966079038, 0.479308067841923, 0.260345966079038, 0.065130102902216,
0.869739794195568, 0.065130102902216, 0.312865496004875, 0.048690315425316, 0.638444188569809,
0.048690315425316, 0.638444188569809, 0.312865496004875
};
static const double G_TRI_W13[] = {
 -0.149570044467670, 0.175615257433204, 0.175615257433204, 0.175615257433204, 0.053347235608839,
  0.053347235608839, 0.053347235608839, 0.077113760890257, 0.077113760890257, 0.077113760890257,
  0.077113760890257, 0.077113760890257, 0.077113760890257 };
//TET
static const double G_TET_X1[] = {
0.25 
};
static const double G_TET_Y1[] = {
0.25 
};
static const double G_TET_Z1[] = {
0.25 
};
static const double G_TET_W1[] = {
1.
};
static const double G_TET_X4[] = {
0.1757281246520584, 0.2445310270213291, 
0.5556470949048655, 0.0240937534217468 
};
static const double G_TET_Y4[] = {
0.5656137776620919, 0.0501800797762026, 
0.1487681308666864, 0.2354380116950194
};
static const double G_TET_Z4[] = {
0.2180665126782654, 0.5635595064952189, 
0.0350112499848832, 0.1833627308416330
};
static const double G_TET_W4[] = {
0.25,0.25,0.25,0.25
};
static const double G_TET_X5[] = {
0.25000000000000000, 0.50000000000000000, 
0.16666666666666667, 0.16666666666666667, 
0.16666666666666667 
};
static const double G_TET_Y5[] = {
0.25000000000000000, 0.16666666666666667, 
0.50000000000000000, 0.16666666666666667, 
0.16666666666666667
};
static const double G_TET_Z5[] = {
0.25000000000000000, 0.16666666666666667, 
0.16666666666666667, 0.50000000000000000, 
0.16666666666666667
};
static const double G_TET_W5[] = {
-0.80000000000000000, 0.45000000000000000, 
0.45000000000000000, 0.45000000000000000, 
0.45000000000000000 
};
static const double G_TET_X10[] = {
  0.5684305841968444, 0.1438564719343852, 0.1438564719343852, 0.1438564719343852,
  0.0000000000000000, 0.5000000000000000, 0.5000000000000000, 0.5000000000000000,
  0.0000000000000000, 0.0000000000000000
};
static const double G_TET_Y10[] = {
  0.1438564719343852, 0.1438564719343852, 0.1438564719343852, 0.5684305841968444, 
  0.5000000000000000, 0.0000000000000000, 0.5000000000000000, 0.0000000000000000, 
  0.5000000000000000, 0.0000000000000000
};
static const double G_TET_Z10[] = {
  0.1438564719343852, 0.1438564719343852, 0.5684305841968444, 0.1438564719343852,
  0.5000000000000000, 0.5000000000000000, 0.0000000000000000, 0.0000000000000000,
  0.0000000000000000, 0.5000000000000000
};
static const double G_TET_W10[] = {
  0.2177650698804054, 0.2177650698804054, 0.2177650698804054, 0.2177650698804054,
  0.0214899534130631, 0.0214899534130631, 0.0214899534130631, 0.0214899534130631,
  0.0214899534130631, 0.0214899534130631
};

static const double G_TET_X45[] = {
0.2500000000000000, 0.6175871903000830, 0.1274709365666390, 0.1274709365666390, 0.1274709365666390,
0.9037635088221031, 0.0320788303926323, 0.0320788303926323, 0.0320788303926323, 0.4502229043567190,
0.0497770956432810, 0.0497770956432810, 0.0497770956432810, 0.4502229043567190, 0.4502229043567190,
0.3162695526014501, 0.1837304473985499, 0.1837304473985499, 0.1837304473985499, 0.3162695526014501,
0.3162695526014501, 0.0229177878448171, 0.2319010893971509, 0.2319010893971509, 0.5132800333608811,
0.2319010893971509, 0.2319010893971509, 0.2319010893971509, 0.0229177878448171, 0.5132800333608811,
0.2319010893971509, 0.0229177878448171, 0.5132800333608811, 0.7303134278075384, 0.0379700484718286,
0.0379700484718286, 0.1937464752488044, 0.0379700484718286, 0.0379700484718286, 0.0379700484718286,
0.7303134278075384, 0.1937464752488044, 0.0379700484718286, 0.7303134278075384, 0.1937464752488044
};
static const double G_TET_Y45[] = {
0.2500000000000000, 0.1274709365666390, 0.1274709365666390, 0.1274709365666390, 0.6175871903000830,
0.0320788303926323, 0.0320788303926323, 0.0320788303926323, 0.9037635088221031, 0.0497770956432810,
0.4502229043567190, 0.0497770956432810, 0.4502229043567190, 0.0497770956432810, 0.4502229043567190,
0.1837304473985499, 0.3162695526014501, 0.1837304473985499, 0.3162695526014501, 0.1837304473985499,
0.3162695526014501, 0.2319010893971509, 0.0229177878448171, 0.2319010893971509, 0.2319010893971509,
0.5132800333608811, 0.2319010893971509, 0.0229177878448171, 0.5132800333608811, 0.2319010893971509,
0.5132800333608811, 0.2319010893971509, 0.0229177878448171, 0.0379700484718286, 0.7303134278075384,
0.0379700484718286, 0.0379700484718286, 0.1937464752488044, 0.0379700484718286, 0.7303134278075384,
0.1937464752488044, 0.0379700484718286, 0.1937464752488044, 0.0379700484718286, 0.7303134278075384
};
static const double G_TET_Z45[] = {
0.2500000000000000, 0.1274709365666390, 0.1274709365666390, 0.6175871903000830, 0.1274709365666390,
0.0320788303926323, 0.0320788303926323, 0.9037635088221031, 0.0320788303926323, 0.0497770956432810,
0.0497770956432810, 0.4502229043567190, 0.4502229043567190, 0.4502229043567190, 0.0497770956432810,
0.1837304473985499, 0.1837304473985499, 0.3162695526014501, 0.3162695526014501, 0.3162695526014501,
0.1837304473985499, 0.2319010893971509, 0.2319010893971509, 0.0229177878448171, 0.2319010893971509,
0.2319010893971509, 0.5132800333608811, 0.5132800333608811, 0.2319010893971509, 0.0229177878448171,
0.0229177878448171, 0.5132800333608811, 0.2319010893971509, 0.0379700484718286, 0.0379700484718286,
0.7303134278075384, 0.0379700484718286, 0.0379700484718286, 0.1937464752488044, 0.1937464752488044,
0.0379700484718286, 0.7303134278075384, 0.7303134278075384, 0.1937464752488044, 0.0379700484718286
};
static const double G_TET_W45[] = {
 -0.2359620398477557, 0.0244878963560562, 0.0244878963560562, 0.0244878963560562, 0.0244878963560562,
  0.0039485206398261, 0.0039485206398261, 0.0039485206398261, 0.0039485206398261, 0.0263055529507371,
  0.0263055529507371, 0.0263055529507371, 0.0263055529507371, 0.0263055529507371, 0.0263055529507371,
  0.0829803830550589, 0.0829803830550589, 0.0829803830550589, 0.0829803830550589, 0.0829803830550589,
  0.0829803830550589, 0.0254426245481023, 0.0254426245481023, 0.0254426245481023, 0.0254426245481023,
  0.0254426245481023, 0.0254426245481023, 0.0254426245481023, 0.0254426245481023, 0.0254426245481023,
  0.0254426245481023, 0.0254426245481023, 0.0254426245481023, 0.0134324384376852, 0.0134324384376852,
  0.0134324384376852, 0.0134324384376852, 0.0134324384376852, 0.0134324384376852, 0.0134324384376852,
  0.0134324384376852, 0.0134324384376852, 0.0134324384376852, 0.0134324384376852, 0.0134324384376852
};

#ifdef __cplusplus
}
#endif

#endif
