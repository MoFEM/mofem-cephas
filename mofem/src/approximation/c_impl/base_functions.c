/** \file base_functions.c

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
#include <phg-quadrule/quad.h>

PetscErrorCode Legendre_polynomials(
  int p,double s,double *diff_s,double *L,double *diffL,const int dim
) {
  PetscFunctionBegin;
  if(dim < 1) SETERRQ(PETSC_COMM_SELF,1,"dim < 1");
  if(dim > 3) SETERRQ(PETSC_COMM_SELF,1,"dim > 3");
  if(p<0) SETERRQ(PETSC_COMM_SELF,1,"p < 0");
  L[0] = 1;
  if(diffL!=NULL) {
    diffL[0*(p+1)+0] = 0;
    if(dim >= 2) {
      diffL[1*(p+1)+0] = 0;
      if(dim == 3) {
        diffL[2*(p+1)+0] = 0;
      }
    }
  }
  if(p==0) PetscFunctionReturn(0);
  L[1] = s;
  if(diffL != NULL) {
    if(diff_s == NULL) {
      SETERRQ(PETSC_COMM_SELF,1,"diff_s == NULL");
    }
    diffL[0*(p+1)+1] = diff_s[0];
    if(dim >= 2) {
      diffL[1*(p+1)+1] = diff_s[1];
      if(dim == 3) {
        diffL[2*(p+1)+1] = diff_s[2];
      }
    }
  }
  if(p==1) PetscFunctionReturn(0);
  int l = 1;
  for(;l<p;l++) {
    double A = ( (2*(double)l+1)/((double)l+1) );
    double B = ( (double)l/((double)l+1) );
    L[l+1] = A*s*L[l] - B*L[l-1];
    if(diffL!=NULL) {
      if(diff_s==NULL) {
        SETERRQ(PETSC_COMM_SELF,1,"diff_s == NULL");
      }
      diffL[0*(p+1)+l+1] = A*(s*diffL[0*(p+1)+l] + diff_s[0]*L[l]) - B*diffL[0*(p+1)+l-1];
      if(dim >= 2) {
        diffL[1*(p+1)+l+1] = A*(s*diffL[1*(p+1)+l] + diff_s[1]*L[l]) - B*diffL[1*(p+1)+l-1];
        if(dim == 3) {
          diffL[2*(p+1)+l+1] = A*(s*diffL[2*(p+1)+l] + diff_s[2]*L[l]) - B*diffL[2*(p+1)+l-1];
        }
      }
    }
  }
  PetscFunctionReturn(0);
}

PetscErrorCode Lobatto_polynomials(
  int p,double s,double *diff_s,double *L,double *diffL,const int dim
) {
  PetscErrorCode ierr;
  PetscFunctionBegin;
  if(dim < 1) SETERRQ(PETSC_COMM_SELF,1,"dim < 1");
  if(dim > 3) SETERRQ(PETSC_COMM_SELF,1,"dim > 3");
  if(p<0) SETERRQ(PETSC_COMM_SELF,1,"p < 0");
  L[0] = 1;
  if(diffL!=NULL) {
    diffL[0*(p+1)+0] = 0;
    if(dim >= 2) {
      diffL[1*(p+1)+0] = 0;
      if(dim == 3) {
        diffL[2*(p+1)+0] = 0;
      }
    }
  }
  if(p==0) PetscFunctionReturn(0);
  L[1] = s;
  if(diffL != NULL) {
    if(diff_s == NULL) {
      SETERRQ(PETSC_COMM_SELF,1,"diff_s == NULL");
    }
    diffL[0*(p+1)+1] = diff_s[0];
    if(dim >= 2) {
      diffL[1*(p+1)+1] = diff_s[1];
      if(dim == 3) {
        diffL[2*(p+1)+1] = diff_s[2];
      }
    }
  }
  if(p==1) PetscFunctionReturn(0);
  double l[p+1];
  ierr = Legendre_polynomials(p,s,NULL,l,NULL,1); CHKERRQ(ierr);
  {
    // Derivatives
    int k = 1;
    for(;k<p;k++) {
      if(diffL!=NULL) {
        if(diff_s==NULL) {
          SETERRQ(PETSC_COMM_SELF,1,"diff_s == NULL");
        }
        diffL[0*(p+1)+k+1] = l[k+1]*diff_s[1];
        if(dim >= 2) {
          diffL[1*(p+1)+k+1] = l[k+1]*diff_s[1];
          if(dim == 3) {
            diffL[2*(p+1)+k+1] = l[k+1]*diff_s[1];
          }
        }
      }
    }
  }
  {
    // Functions
    bzero(&L[2],(p+1)*sizeof(double));
    int nb_gauss_pts = QUAD_1D_TABLE[p]->npoints;
    for(int gg = 0;gg<nb_gauss_pts;gg++) {
      double ksi = QUAD_1D_TABLE[p]->points[1]-QUAD_1D_TABLE[p]->points[0];
      ierr = Legendre_polynomials(
        p,ksi,NULL,l,NULL,1
      ); CHKERRQ(ierr);
      double w = 2*QUAD_2D_TABLE[p]->weights[gg];
      int k = 1;
      for(;k<p;k++) {
        L[k+1] += w*l[k+1];
      }
    }
  }
  PetscFunctionReturn(0);
}
