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
#include <cblas.h>

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
  if(p<2) SETERRQ(PETSC_COMM_SELF,1,"p < 2");
  p -= 2; // polynomial order starts from 2
  double l[p+2];
  ierr = Legendre_polynomials(p+1,s,NULL,l,NULL,1); CHKERRQ(ierr);
  {
    // Derivatives
    int k = 0;
    for(;k<=p;k++) {
      if(diffL!=NULL) {
        if(diff_s==NULL) {
          SETERRQ(PETSC_COMM_SELF,1,"diff_s == NULL");
        }
        double a = l[k+1];
        diffL[0*(p+1)+k] = a*diff_s[0];
        if(dim >= 2) {
          diffL[1*(p+1)+k] = a*diff_s[1];
          if(dim == 3) {
            diffL[2*(p+1)+k] = a*diff_s[2];
          }
        }
      }
    }
  }
  {
    // Functions
    bzero(L,(p+1)*sizeof(double));
    int nb_gauss_pts = QUAD_1D_TABLE[p+2]->npoints;
    double *points = QUAD_1D_TABLE[p+2]->points;
    double *weights = QUAD_1D_TABLE[p+2]->weights;
    s = s+1;
    int gg = 0;
    for(;gg!=nb_gauss_pts;gg++) {
      double ksi = points[2*gg+1];
      double zeta = s*ksi-1;
      ierr = Legendre_polynomials(p+1,zeta,NULL,l,NULL,1); CHKERRQ(ierr);
      double w = s*weights[gg];
      cblas_daxpy(p,w,&l[1],1,&L[0],1);
    }
  }
  {
    int k = 0;
    for(;k<=p;k++) {
      double a = sqrt(k+2-0.5);
      if(L!=NULL) L[k] *= a;
      if(diffL!=NULL) diffL[k] *= a;
    }
  }
  PetscFunctionReturn(0);
}

/// Definitions taken from Hermes2d code

/// kernel functions
#define phi0(x) (-2.0 * 1.22474487139158904909864203735)
#define phi1(x) (-2.0 * 1.58113883008418966599944677222 * (x))
#define phi2(x) (-1.0 / 2.0 * 1.87082869338697069279187436616 * (5 * (x) * (x) - 1))
#define phi3(x) (-1.0 / 2.0 * 2.12132034355964257320253308631 * (7 * (x) * (x) - 3) * (x))
#define phi4(x) (-1.0 / 4.0 * 2.34520787991171477728281505677 * (21 * (x) * (x) * (x) * (x) - 14 * (x) * (x) + 1))
#define phi5(x) (-1.0 / 4.0 * 2.54950975679639241501411205451 * ((33 * (x) * (x) - 30) * (x) * (x) + 5) * (x))
#define phi6(x) (-1.0 / 32.0 * 2.73861278752583056728484891400 * (((429 * (x) * (x) - 495) * (x) * (x) + 135) * (x) * (x) - 5))
#define phi7(x) (-1.0 / 32.0 * 2.91547594742265023543707643877 * (((715 * (x) * (x) - 1001) * (x) * (x) + 385) * (x) * (x) - 35) * (x))
#define phi8(x) (-1.0 / 64.0 * 3.08220700148448822512509619073 * ((((2431 * (x) * (x) - 4004) * (x) * (x) + 2002) * (x) * (x) - 308) * (x) * (x) + 7))
#define phi9(x) (-1.0 / 128.0 * 6.4807406984078603784382721642 * ((((4199 * (x) * (x) - 7956) * (x) * (x) + 4914) * (x) * (x) - 1092) * (x) * (x) + 63) * (x))

static double f_phi0(double x) { return phi0(x); }
static double f_phi1(double x) { return phi1(x); }
static double f_phi2(double x) { return phi2(x); }
static double f_phi3(double x) { return phi3(x); }
static double f_phi4(double x) { return phi4(x); }
static double f_phi5(double x) { return phi5(x); }
static double f_phi6(double x) { return phi6(x); }
static double f_phi7(double x) { return phi7(x); }
static double f_phi8(double x) { return phi8(x); }
static double f_phi9(double x) { return phi9(x); }

static double (*f_phi[])(double x) = {
  f_phi0, f_phi1, f_phi2, f_phi3, f_phi4, f_phi5, f_phi6, f_phi7, f_phi8, f_phi9
};

/// derivatives of kernel functions
#define phi0x(x) (0)
#define phi1x(x) (-2.0 * 1.58113883008418966599944677222)
#define phi2x(x) (-1.0 / 2.0 * 1.87082869338697069279187436616 * (10 * (x)))
#define phi3x(x) (-1.0 / 2.0 * 2.12132034355964257320253308631 * (21.0*(x)*(x)-3.0))
#define phi4x(x) (-1.0 / 4.0 * 2.34520787991171477728281505677 * ((84.0*(x)*(x)-28.0)*(x)))
#define phi5x(x) (-1.0 / 4.0 * 2.54950975679639241501411205451 * ((165.0*(x)*(x)-90.0)*(x)*(x)+5.0))
#define phi6x(x) (-1.0 / 32.0 * 2.73861278752583056728484891400 * (((2574.0*(x)*(x)-1980.0)*(x)*(x)+270.0)*(x)))
#define phi7x(x) (-1.0 / 32.0 * 2.91547594742265023543707643877 * (((5005.0*(x)*(x)-5005.0)*(x)*(x)+1155.0)*(x)*(x)-35.0))
#define phi8x(x) (-1.0 / 64.0 * 3.08220700148448822512509619073 * ((((19448.0*(x)*(x)-24024.0)*(x)*(x)+8008.0)*(x)*(x)-616.0)*(x)))
#define phi9x(x) (-1.0 / 128.0 * 6.4807406984078603784382721642 * ((((37791.0*(x)*(x)-55692.0)*(x)*(x)+24570.0)*(x)*(x)-3276.0)*(x)*(x)-63.0))

static double f_phi0x(double x) { return phi0x(x); }
static double f_phi1x(double x) { return phi1x(x); }
static double f_phi2x(double x) { return phi2x(x); }
static double f_phi3x(double x) { return phi3x(x); }
static double f_phi4x(double x) { return phi4x(x); }
static double f_phi5x(double x) { return phi5x(x); }
static double f_phi6x(double x) { return phi6x(x); }
static double f_phi7x(double x) { return phi7x(x); }
static double f_phi8x(double x) { return phi8x(x); }
static double f_phi9x(double x) { return phi9x(x); }

static double (*f_phix[])(double x) = {
  f_phi0x, f_phi1x, f_phi2x, f_phi3x, f_phi4x, f_phi5x, f_phi6x, f_phi7x, f_phi8x, f_phi9x
};

// /// Legendre polynomials
// #define Legendre0(x) (1.0)
// #define Legendre1(x) (x)
// #define Legendre2(x) (1.0 / 2.0 * (3 * (x) * (x) - 1))
// #define Legendre3(x) (1.0 / 2.0 * (5 * (x) * (x) - 3) * (x))
// #define Legendre4(x) (1.0 / 8.0 * ((35 * (x) * (x) - 30) * (x) * (x) + 3))
// #define Legendre5(x) (1.0 / 8.0 * ((63 * (x) * (x) - 70) * (x) * (x) + 15) * (x))
// #define Legendre6(x) (1.0 / 16.0 * (((231 * (x) * (x) - 315) * (x) * (x) + 105) * (x) * (x) - 5))
// #define Legendre7(x) (1.0 / 16.0 * (((429 * (x) * (x) - 693) * (x) * (x) + 315) * (x) * (x) - 35) * (x))
// #define Legendre8(x) (1.0 / 128.0 * ((((6435 * (x) * (x) - 12012) * (x) * (x) + 6930) * (x) * (x) - 1260) * (x) * (x) + 35))
// #define Legendre9(x) (1.0 / 128.0 * ((((12155 * (x) * (x) - 25740) * (x) * (x) + 18018) * (x) * (x) - 4620) * (x) * (x) + 315) * (x))
// #define Legendre10(x) (1.0 / 256.0 * (((((46189 * (x) * (x) - 109395) * (x) * (x) + 90090) * (x) * (x) - 30030) * (x) * (x) + 3465) * (x) * (x) - 63))
//
// /// derivatives of Legendre polynomials
// #define Legendre0x(x) (0.0)
// #define Legendre1x(x) (1.0)
// #define Legendre2x(x) (3.0 * (x))
// #define Legendre3x(x) (15.0 / 2.0 * (x) * (x) - 3.0 / 2.0)
// #define Legendre4x(x) (5.0 / 2.0 * (x) * (7.0 * (x) * (x) - 3.0))
// #define Legendre5x(x) ((315.0 / 8.0 * (x) * (x) - 105.0 / 4.0) * (x) * (x) + 15.0 / 8.0)
// #define Legendre6x(x) (21.0 / 8.0 * (x) * ((33.0 * (x) * (x) - 30.0) * (x) * (x) + 5.0))
// #define Legendre7x(x) (((3003.0 / 16.0 * (x) * (x) - 3465.0 / 16.0) * (x) * (x) + 945.0 / 16.0) * (x) * (x) - 35.0 / 16.0)
// #define Legendre8x(x) (9.0 / 16.0 * (x) * (((715.0 * (x) * (x) - 1001.0) * (x) * (x) + 385.0) * (x) * (x) - 35.0))
// #define Legendre9x(x) ((((109395.0 / 128.0 * (x) * (x) - 45045.0 / 32.0) * (x) * (x) + 45045.0 / 64.0) * (x) * (x) - 3465.0 / 32.0) * (x) * (x) + 315.0 / 128.0)
// #define Legendre10x(x) (2.0 / 256.0 * (x) * ((((230945.0 * (x) * (x) - 437580.0) * (x) * (x) + 270270.0) * (x) * (x) - 60060.0) * (x) * (x) + 3465.0))
//
// /// first two Lobatto shape functions
// #define l0(x) ((1.0 - (x)) * 0.5)
// #define l1(x) ((1.0 + (x)) * 0.5)
//
// #define l0l1(x) ((1.0 - (x)*(x)) * 0.25)
//
// /// other Lobatto shape functions
// #define l2(x)  (phi0(x) * l0l1(x))
// #define l3(x)  (phi1(x) * l0l1(x))
// #define l4(x)  (phi2(x) * l0l1(x))
// #define l5(x)  (phi3(x) * l0l1(x))
// #define l6(x)  (phi4(x) * l0l1(x))
// #define l7(x)  (phi5(x) * l0l1(x))
// #define l8(x)  (phi6(x) * l0l1(x))
// #define l9(x)  (phi7(x) * l0l1(x))
// #define l10(x) (phi8(x) * l0l1(x))
// #define l11(x) (phi9(x) * l0l1(x))
//
// /// derivatives of Lobatto functions
// #define dl0(x)  (-0.5)
// #define dl1(x)  (0.5)
// #define dl2(x)  (sqrt(3.0/2.0) * Legendre1(x))
// #define dl3(x)  (sqrt(5.0/2.0) * Legendre2(x))
// #define dl4(x)  (sqrt(7.0/2.0) * Legendre3(x))
// #define dl5(x)  (sqrt(9.0/2.0) * Legendre4(x))
// #define dl6(x)  (sqrt(11.0/2.0) * Legendre5(x))
// #define dl7(x)  (sqrt(13.0/2.0) * Legendre6(x))
// #define dl8(x)  (sqrt(15.0/2.0) * Legendre7(x))
// #define dl9(x)  (sqrt(17.0/2.0) * Legendre8(x))
// #define dl10(x) (sqrt(19.0/2.0) * Legendre9(x))
// #define dl11(x) (sqrt(21.0/2.0) * Legendre10(x))


PetscErrorCode LobattoKernel_polynomials(
  int p,double s,double *diff_s,double *L,double *diffL,const int dim
) {
  PetscFunctionBegin;
  if(dim < 1) SETERRQ(PETSC_COMM_SELF,1,"dim < 1");
  if(dim > 3) SETERRQ(PETSC_COMM_SELF,1,"dim > 3");
  if(p<0) SETERRQ(PETSC_COMM_SELF,1,"p < 0");
  if(p>9) SETERRQ(PETSC_COMM_SELF,1,"p > 9");
  if(L) {
    for(int l = 0;l!=p+1;l++) {
      L[l] = f_phi[l](s)*0.25;
    }
  }
  if(diffL!=NULL) {
    if(diff_s==NULL) {
      SETERRQ(PETSC_COMM_SELF,1,"diff_s == NULL");
    }
    for(int l = 0;l!=p+1;l++) {
      double a = f_phix[l](s)*0.25;
      diffL[0*(p+1)+l] = diff_s[0]*a;
      if(dim >= 2) {
        diffL[1*(p+1)+l] = diff_s[1]*a;
        if(dim == 3) {
          diffL[2*(p+1)+l] = diff_s[2]*a;
        }
      }
    }
  }
  PetscFunctionReturn(0);
}
