/** \file l2.c

  Based on Hierarchic Finite Element Bases on Unstructured Tetrahedral
  Meshes, by Mark Ainsworth and Joe Coyle
  Shape functions for MBTRI and H1 approximation

*/

#include <petscsys.h>
#include <cblas.h>

#include <definitions.h>
#include <fem_tools.h>
#include <base_functions.h>
#include <h1_hdiv_hcurl_l2.h>

static PetscErrorCode ierr;

PetscErrorCode L2_Ainsworth_ShapeFunctions_MBTRI(
    int p, double *N, double *diffN, double *L2N, double *diff_L2N, int GDIM,
    PetscErrorCode (*base_polynomials)(int p, double s, double *diff_s,
                                       double *L, double *diffL,
                                       const int dim)) {
  MoFEMFunctionBeginHot;

  int P = NBFACETRI_L2(p);
  if (P == 0)
    MoFEMFunctionReturnHot(0);
  double diff_ksiL01[2], diff_ksiL20[2];
  int dd = 0;
  for (; dd < 2; dd++) {
    diff_ksiL01[dd] = (diffN[1 * 2 + dd] - diffN[0 * 2 + dd]);
    diff_ksiL20[dd] = (diffN[2 * 2 + dd] - diffN[0 * 2 + dd]);
  }
  int ii = 0;
  for (; ii != GDIM; ++ii) {
    int node_shift = ii * 3;
    double ksiL01 = N[node_shift + 1] - N[node_shift + 0];
    double ksiL20 = N[node_shift + 2] - N[node_shift + 0];
    double L01[p + 1], L20[p + 1];
    double diffL01[2 * (p + 1)], diffL20[2 * (p + 1)];
    ierr = base_polynomials(p, ksiL01, diff_ksiL01, L01, diffL01, 2);
    CHKERRQ(ierr);
    ierr = base_polynomials(p, ksiL20, diff_ksiL20, L20, diffL20, 2);
    CHKERRQ(ierr);
    int shift = ii * P;
    int jj = 0;
    int oo = 0;
    for (; oo <= p; oo++) {
      int pp0 = 0;
      for (; pp0 <= oo; pp0++) {
        int pp1 = oo - pp0;
        if (pp1 >= 0) {
          if (L2N != NULL) {
            L2N[shift + jj] = L01[pp0] * L20[pp1];
          }
          if (diff_L2N != NULL) {
            int dd = 0;
            for (; dd < 2; dd++) {
              diff_L2N[2 * shift + 2 * jj + dd] =
                  diffL01[dd * (p + 1) + pp0] * L20[pp1] +
                  L01[pp0] * diffL20[dd * (p + 1) + pp1];
            }
          }
          jj++;
        }
      }
    }
    if (jj != P)
      SETERRQ1(PETSC_COMM_SELF, 1, "wrong order %d", jj);
  }
  MoFEMFunctionReturnHot(0);
}
PetscErrorCode L2_Ainsworth_ShapeFunctions_MBTET(
    int p, double *N, double *diffN, double *L2N, double *diff_L2N, int GDIM,
    PetscErrorCode (*base_polynomials)(int p, double s, double *diff_s,
                                       double *L, double *diffL,
                                       const int dim)) {
  MoFEMFunctionBeginHot;

  int P = NBVOLUMETET_L2(p);
  if (P == 0)
    MoFEMFunctionReturnHot(0);
  double diff_ksiL0[3], diff_ksiL1[3], diff_ksiL2[3];
  int dd = 0;
  if (diffN != NULL) {
    for (; dd < 3; dd++) {
      diff_ksiL0[dd] = (diffN[1 * 3 + dd] - diffN[0 * 3 + dd]);
      diff_ksiL1[dd] = (diffN[2 * 3 + dd] - diffN[0 * 3 + dd]);
      diff_ksiL2[dd] = (diffN[3 * 3 + dd] - diffN[0 * 3 + dd]);
    }
  }
  int ii = 0;
  for (; ii != GDIM; ++ii) {
    int node_shift = ii * 4;
    double ksiL0 = N[node_shift + 1] - N[node_shift + 0];
    double ksiL1 = N[node_shift + 2] - N[node_shift + 0];
    double ksiL2 = N[node_shift + 3] - N[node_shift + 0];
    double L0[p + 1], L1[p + 1], L2[p + 1];
    double diffL0[3 * (p + 1)], diffL1[3 * (p + 1)], diffL2[3 * (p + 1)];
    if (diffN != NULL) {
      ierr = base_polynomials(p, ksiL0, diff_ksiL0, L0, diffL0, 3);
      CHKERRQ(ierr);
      ierr = base_polynomials(p, ksiL1, diff_ksiL1, L1, diffL1, 3);
      CHKERRQ(ierr);
      ierr = base_polynomials(p, ksiL2, diff_ksiL2, L2, diffL2, 3);
      CHKERRQ(ierr);
    } else {
      ierr = base_polynomials(p, ksiL0, NULL, L0, NULL, 3);
      CHKERRQ(ierr);
      ierr = base_polynomials(p, ksiL1, NULL, L1, NULL, 3);
      CHKERRQ(ierr);
      ierr = base_polynomials(p, ksiL2, NULL, L2, NULL, 3);
      CHKERRQ(ierr);
    }
    int shift = ii * P;
    int jj = 0;
    int oo = 0;
    for (; oo <= p; oo++) {
      int pp0 = 0;
      for (; pp0 <= oo; pp0++) {
        int pp1 = 0;
        for (; (pp0 + pp1) <= oo; pp1++) {
          int pp2 = oo - pp0 - pp1;
          if (pp2 >= 0) {
            if (L2N != NULL) {
              L2N[shift + jj] = L0[pp0] * L1[pp1] * L2[pp2];
            }
            if (diff_L2N != NULL) {
              int dd = 0;
              for (; dd < 3; dd++) {
                diff_L2N[3 * shift + 3 * jj + dd] =
                    diffL0[dd * (p + 1) + pp0] * L1[pp1] * L2[pp2] +
                    L0[pp0] * diffL1[dd * (p + 1) + pp1] * L2[pp2] +
                    L0[pp0] * L1[pp1] * diffL2[dd * (p + 1) + pp2];
              }
            }
            jj++;
          }
        }
      }
    }
    if (jj != P)
      SETERRQ2(PETSC_COMM_SELF, 1, "wrong order %d != %d", jj, P);
  }
  MoFEMFunctionReturnHot(0);
}
;