/** \file H1HdivHcurlL2.h
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

#ifndef __H1_H__
#define __H1_H__

#ifdef __cplusplus
extern "C" {
#endif

/// number of dofs for L2 space
#define NBVOLUME_L2(P) ((P+1)*(P+2)*(P+3)/6)
#define NBFACE_L2(P) ((P+1)*(P+2)/2)
#define NBEDGE_L2(P) (P+1)
/// number of dofs on edge for H1 space
#define NBEDGE_H1(P) ((P>0) ? (P-1) : 0)
/// number of dofs on face for H1 space
#define NBFACE_H1(P) ((P>1) ? ((P-2)*(P-1)/2) : 0)
/// number of dofs on volume for H1 space
#define NBVOLUME_H1(P) ((P>2) ? ((P-3)*(P-2)*(P-1)/6) : 0)
#define NBEDGE_Hcurl(P) (P+1)
#define NBFACE_Hcurl(P) ((P>0) ? (P-1)*(P+1) : 0)
#define NBVOLUME_Hcurl(P) ((P>1) ? (P-2)*(P-1)*(P+1)/2 : 0)
#define NBEDGE_Hdiv(P) (0)
#define NBFACE_Hdiv(P) ((P+1)*(P+2)/2)
#define NBVOLUME_Hdiv(P) ((P>1) ? (P-1)*(P+1)*(P+2)/2 : 0)

/** 
 * \brief Calulate Lagrange approximation basis
 *
 * \param p is approximation order
 * \param s is is position [-1,1]
 * \parem L appeoximation functions
 * \param diffL direvatives
 * \param dim dimension
 */
PetscErrorCode Lagrange_basis(int p,double s,double *diff_s,double *L,double *diffL,const int dim);

PetscErrorCode L2_FaceShapeFunctions_MBTRI(int p,double *N,double *diffN,double *L2N,double *diff_L2N,int GDIM);
PetscErrorCode L2_ShapeFunctions_MBTET(int p,double *N,double *diffN,double *L2N,double *diff_L2N,int GDIM);
PetscErrorCode L2_VolumeShapeDiffMBTETinvJ(int base_p,int p,double *volume_diffN,double *invJac,double *volume_diffNinvJac,int GDIM);

/** 
 * \brief H1_EdgeShapeFunctions_MBTRI
 * 
 * \param sense of edges, it is array of inegers dim 3 (3-egdes of triangle)
 * \param p of edges
 */
PetscErrorCode H1_EdgeShapeFunctions_MBTRI(int *sense,int *p,double *N,double *diffN,double *edgeN[3],double *diff_edgeN[3],int GDIM);
PetscErrorCode H1_FaceShapeFunctions_MBTRI(int *face_nodes,int p,double *N,double *diffN,double *faceN,double *diff_faceN,int GDIM);
PetscErrorCode H1_EdgeShapeFunctions_MBTET(int *sense,int *p,double *N,double *diffN,double *edgeN[],double *diff_edgeN[],int GDIM);
PetscErrorCode H1_FaceShapeFunctions_MBTET(int *faces_nodes,int *p,double *N,double *diffN,double *faceN[],double *diff_faceN[],int GDIM);
PetscErrorCode H1_VolumeShapeFunctions_MBTET(int p,double *N,double *diffN,double *volumeN,double *diff_volumeN,int GDIM);
PetscErrorCode H1_EdgeShapeDiffMBTETinvJ(int *base_p,int *p,double *edge_diffN[],double *invJac,double *edge_diffNinvJac[],int GDIM);
PetscErrorCode H1_FaceShapeDiffMBTETinvJ(int *base_p,int *p,double *face_diffN[],double *invJac,double *face_diffNinvJac[],int GDIM);
PetscErrorCode H1_VolumeShapeDiffMBTETinvJ(int base_p,int p,double *volume_diffN,double *invJac,double *volume_diffNinvJac,int GDIM);
PetscErrorCode H1_EdgeGradientOfDeformation_hierachical(int p,double *diffN,double *dofs,double *F);
PetscErrorCode H1_FaceGradientOfDeformation_hierachical(int p,double *diffN,double *dofs,double *F);
PetscErrorCode H1_VolumeGradientOfDeformation_hierachical(int p,double *diffN,double *dofs,double *F);

// Hdiv shape functions

PetscErrorCode Hdiv_EdgeFaceShapeFunctions_MBTET(int *faces_nodes,int *p,double *N,double *diffN,double *PHI_f_e[4][3],int GDIM);
PetscErrorCode Hdiv_FaceBubbleShapeFunctions_MBTET(int *faces_nodes,int *p,double *N,double *diffN,double *PHI_f[4],int GDIM);
PetscErrorCode Hdiv_EdgeBasedVolumeShapeFunctions_MBTET(int *sense,int *p,double *coords,double *N,double *PHI_v_e[6],int GDIM);
PetscErrorCode Hdiv_FaceBasedVolumeShapeFunctions_MBTET(int *faces_nodes,int *p,double *coords,double *N,double *PHI_v_f[4],int GDIM);
PetscErrorCode Hdiv_VolumeBubbleShapeFunctions_MBTET(int *sense,int p,double *N,double *PHI_v,int GDIM);

#ifdef __cplusplus
}
#endif

#endif
