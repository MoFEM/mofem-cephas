/** \file l2.c

  Based on Hierarchic Finite Element Bases on Unstructured Tetrahedral
  Meshes, by Mark Ainsworth and Joe Coyle
  Shape functions for MBTRI and H1 approximation

*/

/**
 * MoFEM is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
 * License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with MoFEM. If not, see <http://www.gnu.org/licenses/>
*/

#include <petscsys.h>
#include <cblas.h>

#include <definitions.h>
#include <fem_tools.h>
#include <base_functions.h>
#include <h1_hdiv_hcurl_l2.h>

PetscErrorCode L2_FaceShapeFunctions_MBTRI(int p,double *N,double *diffN,double *L2N,double *diff_L2N,int GDIM) {
  PetscFunctionBegin;
  PetscErrorCode ierr;
  int P = NBFACETRI_L2_AINSWORTH_COLE(p);
  if(P==0) PetscFunctionReturn(0);
  double diff_ksiL01[2],diff_ksiL20[2];
  int dd = 0;
  for(;dd<2;dd++) {
    diff_ksiL01[dd] = ( diffN[1*2 + dd] - diffN[0*2 + dd] );
    diff_ksiL20[dd] = ( diffN[0*2 + dd] - diffN[2*2 + dd] );
  }
  int ii = 0;
  for(;ii<GDIM;ii++) {
    int node_shift = ii*3;
    double ksiL01 = N[ node_shift+1 ] - N[ node_shift + 0];
    double ksiL20 = N[ node_shift+0 ] - N[ node_shift + 2];
    double L01[ p+1 ],L20[ p+1 ];
    double diffL01[ 2*(p+1) ],diffL20[ 2*(p+1) ];
    ierr = Legendre_polynomials(p,ksiL01,diff_ksiL01,L01,diffL01,2); CHKERRQ(ierr);
    ierr = Legendre_polynomials(p,ksiL20,diff_ksiL20,L20,diffL20,2); CHKERRQ(ierr);
    int shift = ii*P;
    int jj = 0;
    int oo = 0;
    for(;oo<=p;oo++) {
      int pp0 = 0;
      for(;pp0<=oo;pp0++) {
        int pp1 = oo-pp0;
        if(pp1>=0) {
          if(L2N!=NULL) {
            L2N[shift+jj] = L01[pp0]*L20[pp1];
          }
          if(diff_L2N!=NULL) {
            int dd = 0;
            for(;dd<2;dd++) {
              diff_L2N[2*shift+2*jj+dd] = diffL01[dd*(p+1)+pp0]*L20[pp1]+L01[pp0]*diffL20[dd*(p+1)+pp1];
            }
          }
          jj++;
        }
      }
    }
    if(jj!=P) SETERRQ1(PETSC_COMM_SELF,1,"wrong order %d",jj);
  }
  PetscFunctionReturn(0);
}
PetscErrorCode L2_ShapeFunctions_MBTET(int p,double *N,double *diffN,double *L2N,double *diff_L2N,int GDIM) {
  PetscFunctionBegin;
  PetscErrorCode ierr;
  int P = NBVOLUMETET_L2_AINSWORTH_COLE(p);
  if(P==0) PetscFunctionReturn(0);
  double diff_ksiL0[3],diff_ksiL1[3],diff_ksiL2[3];
  int dd = 0;
  for(;dd<3;dd++) {
    diff_ksiL0[dd] = ( diffN[1*3 + dd] - diffN[0*3 + dd] );
    diff_ksiL1[dd] = ( diffN[2*3 + dd] - diffN[0*3 + dd] );
    diff_ksiL2[dd] = ( diffN[3*3 + dd] - diffN[0*3 + dd] );
  }
  int ii = 0;
  for(;ii<GDIM;ii++) {
    int node_shift = ii*4;
    double ksiL0 = N[ node_shift+1 ] - N[ node_shift + 0];
    double ksiL1 = N[ node_shift+2 ] - N[ node_shift + 0];
    double ksiL2 = N[ node_shift+3 ] - N[ node_shift + 0];
    double L0[ p+1 ],L1[ p+1 ],L2[ p+1 ];
    double diffL0[ 3*(p+1) ],diffL1[ 3*(p+1) ],diffL2[ 3*(p+1) ];
    ierr = Legendre_polynomials(p,ksiL0,diff_ksiL0,L0,diffL0,3);  CHKERRQ(ierr);
    ierr = Legendre_polynomials(p,ksiL1,diff_ksiL1,L1,diffL1,3);  CHKERRQ(ierr);
    ierr = Legendre_polynomials(p,ksiL2,diff_ksiL2,L2,diffL2,3);  CHKERRQ(ierr);
    int shift = ii*P;
    int jj = 0;
    int oo = 0;
    for(;oo<=p;oo++) {
      int pp0 = 0;
      for(;pp0<=oo;pp0++) {
        int pp1 = 0;
        for(;(pp0+pp1)<=oo;pp1++) {
          int pp2 = oo - pp0 - pp1;
          if(pp2>=0) {
            if(L2N!=NULL) {
              L2N[shift+jj] = L0[pp0]*L1[pp1]*L2[pp2];
            }
            if(diff_L2N!=NULL) {
              int dd = 0;
              for(;dd<3;dd++) {
                diff_L2N[3*shift+3*jj+dd] =
                diffL0[dd*(p+1)+pp0]*L1[pp1]*L2[pp2]+L0[pp0]*diffL1[dd*(p+1)+pp1]
                *L2[pp2]+L0[pp0]*L1[pp1]*diffL2[dd*(p+1)+pp2];
              }
            }
            jj++;
          }
        }
      }
    }
    if(jj!=P) SETERRQ2(PETSC_COMM_SELF,1,"wrong order %d != %d",jj,P);
  }
  PetscFunctionReturn(0);
}
PetscErrorCode L2_VolumeShapeDiffMBTETinvJ(int base_p,int p,double *volume_diffN,double *invJac,double *volume_diffNinvJac,int GDIM) {
  PetscFunctionBegin;
  int ii,gg;
  for(ii = 0;ii<NBVOLUMETET_L2_AINSWORTH_COLE(p);ii++) {
    for(gg = 0;gg<GDIM;gg++) {
      int shift1 = NBVOLUMETET_L2_AINSWORTH_COLE(base_p)*gg;
      int shift2 = NBVOLUMETET_L2_AINSWORTH_COLE(p)*gg;
      cblas_dgemv(CblasRowMajor,CblasTrans,3,3,1.,
	invJac,3,&(volume_diffN)[3*shift1+3*ii],1,0.,&(volume_diffNinvJac)[3*shift2+3*ii],1);
  }}
  PetscFunctionReturn(0);
}
