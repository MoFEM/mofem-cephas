/** \file base_functions.c

*/

#include <cblas.h>
#include <petscsys.h>
#include <phg-quadrule/quad.h>

#include <definitions.h>

#include <base_functions.h>

static PetscErrorCode ierr;

PetscErrorCode Legendre_polynomials(int p, double s, double *diff_s, double *L,
                                    double *diffL, const int dim) {
  MoFEMFunctionBeginHot;
#ifndef NDEBUG
  if (dim < 1)
    SETERRQ(PETSC_COMM_SELF, MOFEM_INVALID_DATA, "dim < 1");
  if (dim > 3)
    SETERRQ(PETSC_COMM_SELF, MOFEM_INVALID_DATA, "dim > 3");
  if (p < 0)
    SETERRQ(PETSC_COMM_SELF, MOFEM_INVALID_DATA, "p < 0");
#endif // NDEBUG

  L[0] = 1;
  if (diffL != NULL) {
    for (int d = 0; d != dim; ++d) {
      diffL[d * (p + 1) + 0] = 0;
    }
  }
  if (p == 0)
    MoFEMFunctionReturnHot(0);

  L[1] = s;
  if (diffL != NULL) {
#ifndef NDEBUG
    if (diff_s == NULL) {
      SETERRQ(PETSC_COMM_SELF, MOFEM_INVALID_DATA, "diff_s == NULL");
    }
#endif // NDEBUG
    for (int d = 0; d != dim; ++d) {
      diffL[d * (p + 1) + 1] = diff_s[d];
    }
  }
  if (p == 1)
    MoFEMFunctionReturnHot(0);

  int l = 1;
  for (; l < p; l++) {
    double A = ((2 * (double)l + 1) / ((double)l + 1));
    double B = ((double)l / ((double)l + 1));
    L[l + 1] = A * s * L[l] - B * L[l - 1];
    if (diffL != NULL) {
      for (int d = 0; d != dim; ++d) {
        diffL[d * (p + 1) + l + 1] =
            A * (diff_s[d] * L[l] + s * diffL[d * (p + 1) + l]) -
            B * diffL[d * (p + 1) + l - 1];
      }
    }
  }

  MoFEMFunctionReturnHot(0);
}

PetscErrorCode Jacobi_polynomials(int p, double alpha, double x, double t,
                                  double *diff_x, double *diff_t, double *L,
                                  double *diffL, const int dim) {
  MoFEMFunctionBeginHot;
#ifndef NDEBUG
  if (dim < 1)
    SETERRQ(PETSC_COMM_SELF, MOFEM_INVALID_DATA, "dim < 1");
  if (dim > 3)
    SETERRQ(PETSC_COMM_SELF, MOFEM_INVALID_DATA, "dim > 3");
  if (p < 0)
    SETERRQ(PETSC_COMM_SELF, MOFEM_INVALID_DATA, "p < 0");

  if (diffL != NULL) {
    if (diff_x == NULL) {
      SETERRQ(PETSC_COMM_SELF, MOFEM_INVALID_DATA, "diff_s == NULL");
    }
  }

#endif // NDEBUG
  L[0] = 1;
  if (diffL != NULL) {
    diffL[0 * (p + 1) + 0] = 0;
    if (dim >= 2) {
      diffL[1 * (p + 1) + 0] = 0;
      if (dim == 3) {
        diffL[2 * (p + 1) + 0] = 0;
      }
    }
  }
  if (p == 0)
    MoFEMFunctionReturnHot(0);
  L[1] = 2 * x - t + alpha * x;
  if (diffL != NULL) {
    int d = 0;
    for (; d < dim; ++d) {
      double d_t = (diff_t) ? diff_t[d] : 0;
      diffL[d * (p + 1) + 1] = (2 + alpha) * diff_x[d] - d_t;
    }
  }
  if (p == 1)
    MoFEMFunctionReturnHot(0);
  int l = 1;
  for (; l < p; l++) {
    int lp1 = l + 1;
    double a = 2 * lp1 * (lp1 + alpha) * (2 * lp1 + alpha - 2);
    double b = 2 * lp1 + alpha - 1;
    double c = (2 * lp1 + alpha) * (2 * lp1 + alpha - 2);
    double d = 2 * (lp1 + alpha - 1) * (lp1 - 1) * (2 * lp1 + alpha);
    double A = b * (c * (2 * x - t) + alpha * alpha * t) / a;
    double B = d * t * t / a;
    L[lp1] = A * L[l] - B * L[l - 1];
    if (diffL != NULL) {
      int z = 0;
      for (; z < dim; ++z) {
        double d_t = (diff_t) ? diff_t[z] : 0;
        double diffA =
            b * (c * (2 * diff_x[z] - d_t) + alpha * alpha * d_t) / a;
        double diffB = d * 2 * t * d_t / a;
        diffL[z * (p + 1) + lp1] = A * diffL[z * (p + 1) + l] -
                                   B * diffL[z * (p + 1) + l - 1] +
                                   diffA * L[l] - diffB * L[l - 1];
      }
    }
  }
  MoFEMFunctionReturnHot(0);
}

PetscErrorCode IntegratedJacobi_polynomials(int p, double alpha, double x,
                                            double t, double *diff_x,
                                            double *diff_t, double *L,
                                            double *diffL, const int dim) {
  MoFEMFunctionBeginHot;
#ifndef NDEBUG
  if (dim < 1)
    SETERRQ(PETSC_COMM_SELF, MOFEM_INVALID_DATA, "dim < 1");
  if (dim > 3)
    SETERRQ(PETSC_COMM_SELF, MOFEM_INVALID_DATA, "dim > 3");
  if (p < 1)
    SETERRQ(PETSC_COMM_SELF, MOFEM_INVALID_DATA, "p < 1");
#endif // NDEBUG
  L[0] = x;
  if (diffL != NULL) {
    int d = 0;
    for (; d != dim; ++d) {
      diffL[d * p + 0] = diff_x[d];
    }
  }
  if (p == 0)
    MoFEMFunctionReturnHot(0);
  double jacobi[(p + 1)];
  double diff_jacobi[(p + 1) * dim];
  ierr = Jacobi_polynomials(p, alpha, x, t, diff_x, diff_t, jacobi, diff_jacobi,
                            dim);
  CHKERRQ(ierr);
  int l = 1;
  for (; l < p; l++) {
    int i = l + 1;
    const double a = (i + alpha) / ((2 * i + alpha - 1) * (2 * i + alpha));
    const double b = alpha / ((2 * i + alpha - 2) * (2 * i + alpha));
    const double c = (i - 1) / ((2 * i + alpha - 2) * (2 * i + alpha - 1));
    L[l] = a * jacobi[i] + b * t * jacobi[i - 1] - c * t * t * jacobi[i - 2];
    if (diffL != NULL) {
      int d = 0;
      for (; d != dim; ++d) {
        diffL[d * p + l] = a * diff_jacobi[d * (p + 1) + i] +
                            b * (t * diff_jacobi[d * (p + 1) + i - 1] +
                                 diff_t[d] * jacobi[i - 1]) -
                            c * (t * t * diff_jacobi[d * (p + 1) + i - 2] +
                                 2 * t * diff_t[d] * jacobi[i - 2]);
      }
    }
  }
  MoFEMFunctionReturnHot(0);
}

PetscErrorCode Lobatto_polynomials(int p, double s, double *diff_s, double *L,
                                   double *diffL, const int dim) {

  MoFEMFunctionBeginHot;
#ifndef NDEBUG
  if (dim < 1)
    SETERRQ(PETSC_COMM_SELF, MOFEM_INVALID_DATA, "dim < 1");
  if (dim > 3)
    SETERRQ(PETSC_COMM_SELF, MOFEM_INVALID_DATA, "dim > 3");
  if (p < 2)
    SETERRQ(PETSC_COMM_SELF, MOFEM_INVALID_DATA, "p < 2");
#endif // NDEBUG
  double l[p + 1];
  ierr = Legendre_polynomials(p, s, NULL, l, NULL, 1);
  CHKERRQ(ierr);

  L[0] = 1;
  if (diffL != NULL) {
    for (int d = 0; d != dim; ++d) {
      diffL[d * (p + 1) + 0] = 0;
    }
  }
  L[1] = s;
  if (diffL != NULL) {
#ifndef NDEBUG
    if (diff_s == NULL) {
      SETERRQ(PETSC_COMM_SELF, MOFEM_INVALID_DATA, "diff_s == NULL");
    }
#endif // NDEBUG
    for (int d = 0; d != dim; ++d) {
      diffL[d * (p + 1) + 1] = diff_s[d];
    }
  }

  // Integrated Legendre
  for (int k = 2; k <= p; k++) {
    const double factor = 2 * (2 * k - 1);
    L[k] = 1.0 / factor * (l[k] - l[k - 2]);
  }

  if (diffL != NULL) {
    for (int k = 2; k <= p; k++) {
      double a = l[k - 1] / 2.;
      for (int d = 0; d != dim; ++d) {
        diffL[d * (p + 1) + k] = a * diff_s[d];
      }
    }
  }
  MoFEMFunctionReturnHot(0);
}

static double f_phi0(double x) { return LOBATTO_PHI0(x); }
static double f_phi1(double x) { return LOBATTO_PHI1(x); }
static double f_phi2(double x) { return LOBATTO_PHI2(x); }
static double f_phi3(double x) { return LOBATTO_PHI3(x); }
static double f_phi4(double x) { return LOBATTO_PHI4(x); }
static double f_phi5(double x) { return LOBATTO_PHI5(x); }
static double f_phi6(double x) { return LOBATTO_PHI6(x); }
static double f_phi7(double x) { return LOBATTO_PHI7(x); }
static double f_phi8(double x) { return LOBATTO_PHI8(x); }
static double f_phi9(double x) { return LOBATTO_PHI9(x); }

static double (*f_phi[])(double x) = {f_phi0, f_phi1, f_phi2, f_phi3, f_phi4,
                                      f_phi5, f_phi6, f_phi7, f_phi8, f_phi9};

static double f_phi0x(double x) { return LOBATTO_PHI0X(x); }
static double f_phi1x(double x) { return LOBATTO_PHI1X(x); }
static double f_phi2x(double x) { return LOBATTO_PHI2X(x); }
static double f_phi3x(double x) { return LOBATTO_PHI3X(x); }
static double f_phi4x(double x) { return LOBATTO_PHI4X(x); }
static double f_phi5x(double x) { return LOBATTO_PHI5X(x); }
static double f_phi6x(double x) { return LOBATTO_PHI6X(x); }
static double f_phi7x(double x) { return LOBATTO_PHI7X(x); }
static double f_phi8x(double x) { return LOBATTO_PHI8X(x); }
static double f_phi9x(double x) { return LOBATTO_PHI9X(x); }

static double (*f_phix[])(double x) = {f_phi0x, f_phi1x, f_phi2x, f_phi3x,
                                       f_phi4x, f_phi5x, f_phi6x, f_phi7x,
                                       f_phi8x, f_phi9x};

// /// Legendre polynomials
// #define Legendre0(x) (1.0)
// #define Legendre1(x) (x)
// #define Legendre2(x) (1.0 / 2.0 * (3 * (x) * (x) - 1))
// #define Legendre3(x) (1.0 / 2.0 * (5 * (x) * (x) - 3) * (x))
// #define Legendre4(x) (1.0 / 8.0 * ((35 * (x) * (x) - 30) * (x) * (x) + 3))
// #define Legendre5(x) (1.0 / 8.0 * ((63 * (x) * (x) - 70) * (x) * (x) + 15) *
// (x)) #define Legendre6(x) (1.0 / 16.0 * (((231 * (x) * (x) - 315) * (x) * (x)
// + 105) * (x) * (x) - 5)) #define Legendre7(x) (1.0 / 16.0 * (((429 * (x) *
// (x) - 693) * (x) * (x) + 315) * (x) * (x) - 35) * (x)) #define Legendre8(x)
// (1.0 / 128.0 * ((((6435 * (x) * (x) - 12012) * (x) * (x) + 6930) * (x) * (x)
// - 1260) * (x) * (x) + 35)) #define Legendre9(x) (1.0 / 128.0 * ((((12155 *
// (x) * (x) - 25740) * (x) * (x) + 18018) * (x) * (x) - 4620) * (x) * (x) +
// 315) * (x)) #define Legendre10(x) (1.0 / 256.0 * (((((46189 * (x) * (x) -
// 109395) * (x) * (x) + 90090) * (x) * (x) - 30030) * (x) * (x) + 3465) * (x) *
// (x) - 63))
//
// /// derivatives of Legendre polynomials
// #define Legendre0x(x) (0.0)
// #define Legendre1x(x) (1.0)
// #define Legendre2x(x) (3.0 * (x))
// #define Legendre3x(x) (15.0 / 2.0 * (x) * (x) - 3.0 / 2.0)
// #define Legendre4x(x) (5.0 / 2.0 * (x) * (7.0 * (x) * (x) - 3.0))
// #define Legendre5x(x) ((315.0 / 8.0 * (x) * (x) - 105.0 / 4.0) * (x) * (x)
// + 15.0 / 8.0) #define Legendre6x(x) (21.0 / 8.0 * (x) * ((33.0 * (x) * (x)
// - 30.0) * (x) * (x) + 5.0)) #define Legendre7x(x) (((3003.0 / 16.0 * (x) *
// (x) - 3465.0 / 16.0) * (x) * (x) + 945.0 / 16.0) * (x) * (x) - 35.0 / 16.0)
// #define Legendre8x(x) (9.0 / 16.0 * (x) * (((715.0 * (x) * (x) - 1001.0) *
// (x) * (x) + 385.0) * (x) * (x) - 35.0)) #define Legendre9x(x) ((((109395.0 /
// 128.0 * (x) * (x) - 45045.0 / 32.0) * (x) * (x) + 45045.0 / 64.0) * (x) * (x)
// - 3465.0 / 32.0) * (x) * (x) + 315.0 / 128.0) #define Legendre10x(x) (2.0 /
// 256.0 * (x) * ((((230945.0 * (x) * (x) - 437580.0) * (x) * (x) + 270270.0) *
// (x) * (x) - 60060.0) * (x) * (x) + 3465.0))
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

PetscErrorCode LobattoKernel_polynomials(int p, double s, double *diff_s,
                                         double *L, double *diffL,
                                         const int dim) {
  MoFEMFunctionBeginHot;
#ifndef NDEBUG
  if (dim < 1)
    SETERRQ(PETSC_COMM_SELF, MOFEM_INVALID_DATA, "dim < 1");
  if (dim > 3)
    SETERRQ(PETSC_COMM_SELF, MOFEM_INVALID_DATA, "dim > 3");
  if (p < 0)
    SETERRQ(PETSC_COMM_SELF, MOFEM_INVALID_DATA, "p < 0");
  if (p > 9)
    SETERRQ(PETSC_COMM_SELF, MOFEM_NOT_IMPLEMENTED,
            "Polynomial beyond order 9 is not implemented");
  if (diffL != NULL) {
    if (diff_s == NULL) {
      SETERRQ(PETSC_COMM_SELF, MOFEM_INVALID_DATA, "diff_s == NULL");
    }
  }

#endif // NDEBUG
  if (L) {
    int l = 0;
    for (; l != p + 1; l++) {
      L[l] = f_phi[l](s);
    }
  }
  if (diffL != NULL) {
    int l = 0;
    for (; l != p + 1; l++) {
      double a = f_phix[l](s);
      int d = 0;
      for (; d < dim; ++d) {
        diffL[d * (p + 1) + l] = diff_s[d] * a;
      }
    }
  }
  MoFEMFunctionReturnHot(0);
}
