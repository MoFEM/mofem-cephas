/* Copyright (C) 2009, Lukasz Kaczmarczyk (likask AT wp.pl)
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
#include "assert.h"

#ifdef __cplusplus
extern "C" {
#endif

enum phisical_equation_volume { hooke = 0, stvenant_kirchhoff = 1,neohookean = 2,eberleinholzapfel1 = 3};
void set_PhysicalEquationNumber(enum phisical_equation_volume eq);
enum phisical_equation_volume get_PhysicalEquationNumber();

PetscErrorCode StrainEnergy(double lambda,double mu,__CLPK_doublecomplex *xF,__CLPK_doublecomplex *xC,__CLPK_doublecomplex *xJ,__CLPK_doublecomplex *xPsi,void *ctx);
PetscErrorCode PiolaKirhoiff2(double lambda,double mu,__CLPK_doublecomplex *xF,__CLPK_doublecomplex *xC,__CLPK_doublecomplex *xJ,__CLPK_doublecomplex *xS,void *ctx);
PetscErrorCode PiolaKirhoiff1(double lambda,double mu,__CLPK_doublecomplex *xF,__CLPK_doublecomplex *xS,__CLPK_doublecomplex *xP);
PetscErrorCode CauchyStress(__CLPK_doublecomplex *xF,__CLPK_doublecomplex *xJ,__CLPK_doublecomplex *xP,__CLPK_doublecomplex *xCauchyStress);
PetscErrorCode ElshebyStress(__CLPK_doublecomplex *xPsi,__CLPK_doublecomplex *xF,__CLPK_doublecomplex *xP,__CLPK_doublecomplex *xSigma);

PetscErrorCode TakeIm(__CLPK_doublecomplex *xA,double *imA);
PetscErrorCode TakeRe(__CLPK_doublecomplex *xA,double *reA);
PetscErrorCode SpatialGradientOfDeformation(__CLPK_doublecomplex *xh,__CLPK_doublecomplex *inv_xH,__CLPK_doublecomplex *xF);
PetscErrorCode CauchyGreenDeformation(__CLPK_doublecomplex *xF,__CLPK_doublecomplex *xC);
PetscErrorCode PiolaKrihoff1_PullBack(__CLPK_doublecomplex *det_xH,__CLPK_doublecomplex *inv_xH,__CLPK_doublecomplex *xP,__CLPK_doublecomplex *xP_PullBack);
PetscErrorCode ElshebyStress_PullBack(__CLPK_doublecomplex *det_xH,__CLPK_doublecomplex *inv_xH,__CLPK_doublecomplex *xStress,__CLPK_doublecomplex *xStress_PullBack);

PetscErrorCode ThermalDeformationGradient(double alpha,__CLPK_doublecomplex xT,__CLPK_doublecomplex *xF);

PetscErrorCode Calulate_Stresses_at_GaussPoint(int *order_edge,int *order_face,int order_volume,double alpha,double lambda,double mu,void *matctx,
  double *diffN,double *diffN_edge[],double *diffN_face[],double *diffN_volume,
  double *dofs_X,double *dofs_x_node,double *dofs_x_edge[],double *dofs_x_face[],double *dofs_x_volume,
  //temperature
  double termal_expansion,
  double *N,double *N_edge[],double *N_face[],double *N_volume,
  int *order_T_edge,int *order_T_face,int order_T_volume,
  double *dofs_T,double *dofs_T_edge[],double *dofs_T_face[],double *dofs_T_volume,
  //
  double *Piola1Stress,double *CauhyStress,double *EshelbyStress,double *Psi,double *J,
  int gg);
PetscErrorCode Fint_Hh_hierarchical(
  int *order_edge,int *order_face,int order_volume,double alpha,double lambda,double mu,void *matctx,
  double *diffN,double *diffN_edge[],double *diffN_face[],double *diffN_volume,
  double *dofs_X,double *dofs_x_node,double *dofs_iX,double *dofs_ix_node,
  double *dofs_x_edge[],double *dofs_x_face[],double *dofs_x_volume,
  //temperature
  double termal_expansion,
  double *N,double *N_edge[],double *N_face[],double *N_volume,
  int *order_T_edge,int *order_T_face,int order_T_volume,
  double *dofs_T,double *dofs_T_edge[],double *dofs_T_face[],double *dofs_T_volume,
  //rhs
  double *Fint_H,double *Fint_h,
  double *Fint_h_edge[],double *Fint_h_face[],double *Fint_h_volume,
  double *Fint_iH,double *Fint_ih,
  double *Fint_ih_edge[],double *Fint_ih_face[],double *Fint_ih_volume,
  int G_DIM,const double *G_W);
PetscErrorCode Tangent_HH_hierachical(int *order_edge,int *order_face,int order_volume,double alpha,double eps,double lambda,double mu,void *matctx,
  double *diffN,double *diffN_edge[],double *diffN_face[],double *diffN_volume,
  double *dofs_X,double *dofs_x_node,double *dofs_x_edge[],double *dofs_x_face[],double *dofs_x_volume,
  //temperature
  double termal_expansion,
  double *N,double *N_edge[],double *N_face[],double *N_volume,
  int *order_T_edge,int *order_T_face,int order_T_volume,
  double *dofs_T,double *dofs_T_edge[],double *dofs_T_face[],double *dofs_T_volume,
  //
  double *K,double *Koff,double *Koff_edge[6],double *Koff_face[4],double *Koff_volume,int G_DIM,const double *G_W);
PetscErrorCode Tangent_hh_hierachical(int *order_edge,int *order_face,int order_volume,double alpha,double eps,double lambda,double mu,void *matctx,
  double *diffN,double *diffN_edge[],double *diffN_face[],double *diffN_volume,
  double *dofs_X,double *dofs_x_node,double *dofs_x_edge[],double *dofs_x_face[],double *dofs_x_volume,
  //temperature
  double termal_expansion,
  double *N,double *N_edge[],double *N_face[],double *N_volume,
  int *order_T_edge,int *order_T_face,int order_T_volume,
  double *dofs_T,double *dofs_T_edge[],double *dofs_T_face[],double *dofs_T_volume,
  //
  double *K,double *Koff,double *K_edge[6],double *K_face[4],double *K_volume,int G_DIM,const double *G_W);
PetscErrorCode Tangent_hh_hierachical_edge(int *order_edge,int *order_face,int order_volume,double alpha,double eps,double lambda,double mu,void *matctx,
  double *diffN,double *diffN_edge[],double *diffN_face[],double *diffN_volume,
  double *dofs_X,double *dofs_x_node,double *dofs_x_edge[],double *dofs_x_face[],double *dofs_x_volume,
  //temperature
  double termal_expansion,
  double *N,double *N_edge[],double *N_face[],double *N_volume,
  int *order_T_edge,int *order_T_face,int order_T_volume,
  double *dofs_T,double *dofs_T_edge[],double *dofs_T_face[],double *dofs_T_volume,
  //
  double *K[6],double *Koff[6],
  double *K_edge[6][6],double *K_face[4][6],double *K_volume[6],
  int G_DIM,const double *G_W);
PetscErrorCode Tangent_hh_hierachical_face(int *order_edge,int *order_face,int order_volume,double alpha,double eps,double lambda,double mu,void *matctx,
  double *diffN,double *diffN_edge[],double *diffN_face[],double *diffN_volume,
  double *dofs_X,double *dofs_x_node,double *dofs_x_edge[],double *dofs_x_face[],double *dofs_x_volume,
  //temperature
  double termal_expansion,
  double *N,double *N_edge[],double *N_face[],double *N_volume,
  int *order_T_edge,int *order_T_face,int order_T_volume,
  double *dofs_T,double *dofs_T_edge[],double *dofs_T_face[],double *dofs_T_volume,
  //
  double *K[4],double *Koff[4],
  double *K_edge[6][4],double *K_face[4][4],double *K_volume[4],
  int G_DIM,const double *G_W);
PetscErrorCode Tangent_hh_hierachical_volume(int *order_edge,int *order_face,int order_volume,double alpha,double eps,double lambda,double mu,void *matctx,
  double *diffN,double *diffN_edge[],double *diffN_face[],double *diffN_volume,
  double *dofs_X,double *dofs_x_node,double *dofs_x_edge[],double *dofs_x_face[],double *dofs_x_volume,
  //temperature
  double termal_expansion,
  double *N,double *N_edge[],double *N_face[],double *N_volume,
  int *order_T_edge,int *order_T_face,int order_T_volume,
  double *dofs_T,double *dofs_T_edge[],double *dofs_T_face[],double *dofs_T_volume,
  //
  double *K,double *Koff,double *K_edge[6],double *K_face[4],double *K_volume,int G_DIM,const double *G_W);
PetscErrorCode Traction_hierarchical(int order,int *order_edge,
  double *N,double *N_face,double *N_edge[],
  double *t,double *t_edge[],double *t_face,
  double *traction,int gg);
PetscErrorCode Fext_h_hierarchical(int order,int *order_edge,
  double *N,double *N_face,double *N_edge[],
  double *diffN,double *diffN_face,double *diffN_edge[],
  double *t,double *t_edge[],double *t_face,
  double *dofs_x,double *dofs_x_edge[],double *dofs_x_face,
  double *idofs_x,double *idofs_x_edge[],double *idofs_x_face,
  double *Fext,double *Fext_egde[],double *Fext_face,
  double *iFext,double *iFext_egde[],double *iFext_face,
  int g_dim,const double *g_w);
PetscErrorCode KExt_hh_hierarchical(double eps,int order,int *order_edge,
  double *N,double *N_face,double *N_edge[],
  double *diffN,double *diffN_face,double *diffN_edge[],
  double *t,double *t_edge[],double *t_face,
  double *dofs_x,double *dofs_x_edge[],double *dofs_x_face,
  double *KExt_hh,double* KExt_egdeh[3],double *KExt_faceh,
  int g_dim,const double *g_w);
PetscErrorCode KExt_hh_hierarchical_edge(double eps,int order,int *order_edge,
  double *N,double *N_face,double *N_edge[],
  double *diffN,double *diffN_face,double *diffN_edge[],
  double *t,double *t_edge[],double *t_face,
  double *dofs_x,double *dofs_x_edge[],double *dofs_x_face,
  double *Khext_edge[3],double *KExt_edgeegde[3][3],double *KExt_faceedge[3],
  int g_dim,const double *g_w);
PetscErrorCode KExt_hh_hierarchical_face(double eps,int order,int *order_edge,
  double *N,double *N_face,double *N_edge[],
  double *diffN,double *diffN_face,double *diffN_edge[],
  double *t,double *t_edge[],double *t_face,
  double *dofs_x,double *dofs_x_edge[],double *dofs_x_face,
  double *KExt_hface,double *KExt_egdeface[3],double *KExt_faceface,
  int g_dim,const double *g_w);
PetscErrorCode Fext_H(int order,int *order_edge,
  double *N,double *N_face,double *N_edge[],
  double *diffN,double *diffN_face,double *diffN_edge[],
  double *t,double *t_edge[],double *t_face,
  double *dofs_X,double *idofs_X,
  double *dofs_x,double *dofs_x_edge[],double *dofs_x_face,
  double *idofs_x,double *idofs_x_edge[],double *idofs_x_face,
  double *Fext,double *iFext,int g_dim,const double *g_w);
PetscErrorCode KExt_HH(double eps,int order,int *order_edge,
  double *N,double *N_face,double *N_edge[],
  double *diffN,double *diffN_face,double *diffN_edge[],
  double *t,double *t_edge[],double *t_face,
  double *dofs_X,
  double *dofs_x,double *dofs_x_edge[],double *dofs_x_face,
  double *KExt_HH,int g_dim,const double *g_w);

//quality
void set_qual_ver(int ver);
int get_qual_ver();
PetscErrorCode quality_volume_length_F(double alpha,double *alpha2,double gamma,double *diffN,
  double *coords_edges,double *dofs_X,double *dofs_x,double *dofs_iX,double *dofs_ix,double *quality0,double *quality,double *b,
  double *F,double *iF);
int quality_volume_length_K(double eps,double alpha,double *alpha2,double gamma,double *diffN,double *coords_edges,double *dofs_X,double *dofs_x,double *K,double *Koff);

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
