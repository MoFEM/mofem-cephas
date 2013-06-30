/* Copyright (C) 2009, Lukasz Kaczmarczyk (likask AT civil.gla.ac.uk)
 * --------------------------------------------------------------
 * FIXME: DESCRIPTION
 */

/* This file is part of mofem.
 * mofem is free software: you can redistribute it and/or modify it under
 * the terms of the GNU Lesser General Public License as published by the
 * Free Software Foundation, either version 3 of the License, or (at your
 * option) any later version.
 *
 * mofem is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
 * License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with mofem. If not, see <http://www.gnu.org/licenses/>. */


#ifndef __COMPLEX_FOR_LAZY_H__
#define __COMPLEX_FOR_LAZY_H__

#include "FEM.h"
#include "H1HdivHcurlL2.h"
#include "complex.h"

enum phisical_equation_volume { hooke = 0, stvenant_kirchhoff = 1,neohookean = 2,eberleinholzapfel1 = 3};

#ifdef __cplusplus
extern "C" {
#endif

void TakeIm(__CLPK_doublecomplex *xA,double *imA);
void TakeRe(__CLPK_doublecomplex *xA,double *reA);

void SpatialGradientOfDeformation(__CLPK_doublecomplex *xh,__CLPK_doublecomplex *inv_xH,__CLPK_doublecomplex *xF);
void CauchyGreenDeformation(__CLPK_doublecomplex *xF,__CLPK_doublecomplex *xC);
void PiolaKrihoff1_PullBack(__CLPK_doublecomplex *det_xH,__CLPK_doublecomplex *inv_xH,__CLPK_doublecomplex *xP,__CLPK_doublecomplex *xP_PullBack);
void ElshebyStress_PullBack(__CLPK_doublecomplex *det_xH,__CLPK_doublecomplex *inv_xH,__CLPK_doublecomplex *xStress,__CLPK_doublecomplex *xStress_PullBack);
void ComputeEneregyAtGauss_Point_hierarchical(
  int *order_edge,int *order_face,int order_volume,double alpha,double lambda,double mu,void *matctx,
  double *diffN,double *diffN_edge[],double *diffN_face[],double *diffN_volume,
  double *dofs_X,double *dofs_x_node,double *dofs_x_edge[],double *dofs_x_face[],double *dofs_x_volume,
  double *rePsi,double *reP,double *reSigma,double *reF,double *intPsi,
  int G_DIM,const double *G_W);

void Fint_Hh(double alpha,double lambda,double mu,void *matctx,double *diffN,double *dofs_X,double *dofs_x,double *dofs_iX,double *dofs_ix,double *Fint_H,double *Fint_h,double *Fint_iH,double *Fint_ih);
void Fext_h(double alpha,double *NTRI,double *tractions,double *Fext);
void Tangent_hh(double alpha,double eps,double lambda,double mu,void *matctx,double *diffN,double *dofs_X,double *dofs_x,double *K,double *Koff);
void Tangent_HH(double alpha,double eps,double lambda,double mu,void *matctx,double *diffN,double *dofs_X,double *dofs_x,double *K,double *Koff);
void Fint_Hh_hierarchical(int *order_edge,int *order_face,int order_volume,double alpha,double lambda,double mu,void *matctx,
  double *diffN,double *diffN_edge[],double *diffN_face[],double *diffN_volume,
  double *dofs_X,double *dofs_x_node,double *dofs_iX,double *dofs_ix_node,
  double *dofs_x_edge[],double *dofs_x_face[],double *dofs_x_volume,
  double *Fint_H,double *Fint_h,
  double *Fint_h_edge[],double *Fint_h_face[],double *Fint_h_volume,
  double *Fint_iH,double *Fint_ih,
  double *Fint_ih_edge[],double *Fint_ih_face[],double *Fint_ih_volume,
  int G_DIM,const double *G_W);
void Tangent_HH_hierachical(int *order_edge,int *order_face,int order_volume,double alpha,double eps,double lambda,double mu,void *matctx,
  double *diffN,double *diffN_edge[],double *diffN_face[],double *diffN_volume,
  double *dofs_X,double *dofs_x_node,double *dofs_x_edge[],double *dofs_x_face[],double *dofs_x_volume,
  double *K,double *Koff,double *Koff_edge[6],double *Koff_face[4],double *Koff_volume,int G_DIM,const double *G_W);
void Tangent_hh_hierachical(int *order_edge,int *order_face,int order_volume,double alpha,double eps,double lambda,double mu,void *matctx,
  double *diffN,double *diffN_edge[],double *diffN_face[],double *diffN_volume,
  double *dofs_X,double *dofs_x_node,double *dofs_x_edge[],double *dofs_x_face[],double *dofs_x_volume,
  double *K,double *Koff,double *K_edge[6],double *K_face[4],double *K_volume,int G_DIM,const double *G_W);
void Tangent_hh_hierachical_edge(int *order_edge,int *order_face,int order_volume,double alpha,double eps,double lambda,double mu,void *matctx,
  double *diffN,double *diffN_edge[],double *diffN_face[],double *diffN_volume,
  double *dofs_X,double *dofs_x_node,double *dofs_x_edge[],double *dofs_x_face[],double *dofs_x_volume,
  double *K[6],double *Koff[6],
  double *K_edge[6][6],double *K_face[4][6],double *K_volume[6],
  int G_DIM,const double *G_W);
void Tangent_hh_hierachical_face(int *order_edge,int *order_face,int order_volume,double alpha,double eps,double lambda,double mu,void *matctx,
  double *diffN,double *diffN_edge[],double *diffN_face[],double *diffN_volume,
  double *dofs_X,double *dofs_x_node,double *dofs_x_edge[],double *dofs_x_face[],double *dofs_x_volume,
  double *K[4],double *Koff[4],
  double *K_edge[6][4],double *K_face[4][4],double *K_volume[4],
  int G_DIM,const double *G_W);
void Tangent_hh_hierachical_volume(int *order_edge,int *order_face,int order_volume,double alpha,double eps,double lambda,double mu,void *matctx,
  double *diffN,double *diffN_edge[],double *diffN_face[],double *diffN_volume,
  double *dofs_X,double *dofs_x_node,double *dofs_x_edge[],double *dofs_x_face[],double *dofs_x_volume,
  double *K,double *Koff,double *K_edge[6],double *K_face[4],double *K_volume,int G_DIM,const double *G_W);

void EdgeGradientOfDeformation_hierachical(int p,double *diffN,double *dofs,double *F);
void FaceGradientOfDeformation_hierachical(int p,double *diffN,double *dofs,double *F);
void VolumeGradientOfDeformation_hierachical(int p,double *diffN,double *dofs,double *F);

//Fibres, by R. Eberlein, G.A. Holzapfel, and C.A.J. Schulze-Bauer. 
//An anisotropic model for annulus tissue and enhanced finite element analyses of intact lumbar disc bodies. 
//Computer Methods in Biomechanics and Biomedical Engineering, 4(3):209–229, 2001. ISSN 1025-5842.
typedef struct {
  enum phisical_equation_volume eq_solid;
  double k1,k2;
  //those are in material space
  double fibre_vector_a1[3];
  double fibre_vector_a2[3];
} ctx_EberleinHolzapfel1;


#ifdef __cplusplus
}
#endif

#endif // __COMPLEX_FOR_LAZY_H__
