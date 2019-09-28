/** \file h1.c

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
#include <definitions.h>
#include <cblas.h>

#include <h1_hdiv_hcurl_l2.h>

static PetscErrorCode ierr;

PetscErrorCode H1_EdgeShapeFunctions_MBTRI(
    int *sense, int *p, double *N, double *diffN, double *edgeN[3],
    double *diff_edgeN[3], int GDIM,
    PetscErrorCode (*base_polynomials)(int p, double s, double *diff_s,
                                       double *L, double *diffL,
                                       const int dim)) {
  MoFEMFunctionBeginHot;

  double *edgeN01 = NULL, *edgeN12 = NULL, *edgeN20 = NULL;
  if (edgeN != NULL) {
    edgeN01 = edgeN[0];
    edgeN12 = edgeN[1];
    edgeN20 = edgeN[2];
  }
  double *diff_edgeN01 = NULL, *diff_edgeN12 = NULL, *diff_edgeN20 = NULL;
  if (diff_edgeN != NULL) {
    diff_edgeN01 = diff_edgeN[0];
    diff_edgeN12 = diff_edgeN[1];
    diff_edgeN20 = diff_edgeN[2];
  }
  int P[3];
  int ee = 0;
  for (; ee < 3; ee++) {
    P[ee] = NBEDGE_H1(p[ee]);
  }
  int dd = 0;
  double diff_ksi01[2], diff_ksi12[2], diff_ksi20[2];
  if (diff_edgeN != NULL) {
    for (; dd < 2; dd++) {
      diff_ksi01[dd] = (diffN[1 * 2 + dd] - diffN[0 * 2 + dd]) * sense[0];
      diff_ksi12[dd] = (diffN[2 * 2 + dd] - diffN[1 * 2 + dd]) * sense[1];
      diff_ksi20[dd] = (diffN[0 * 2 + dd] - diffN[2 * 2 + dd]) * sense[2];
    }
  }
  int ii = 0;
  for (; ii < GDIM; ii++) {
    int node_shift = ii * 3;
    double ksi01 = (N[node_shift + 1] - N[node_shift + 0]) * sense[0];
    double ksi12 = (N[node_shift + 2] - N[node_shift + 1]) * sense[1];
    double ksi20 = (N[node_shift + 0] - N[node_shift + 2]) * sense[2];
    double L01[p[0] + 1], L12[p[1] + 1], L20[p[2] + 1];
    double diffL01[2 * (p[0] + 1)], diffL12[2 * (p[1] + 1)],
        diffL20[2 * (p[2] + 1)];
    ierr = base_polynomials(p[0], ksi01, diff_ksi01, L01, diffL01, 2);
    CHKERRQ(ierr);
    ierr = base_polynomials(p[1], ksi12, diff_ksi12, L12, diffL12, 2);
    CHKERRQ(ierr);
    ierr = base_polynomials(p[2], ksi20, diff_ksi20, L20, diffL20, 2);
    CHKERRQ(ierr);
    int shift;
    if (edgeN != NULL) {
      // edge01
      shift = ii * (P[0]);
      cblas_dcopy(P[0], L01, 1, &edgeN01[shift], 1);
      cblas_dscal(P[0], N[node_shift + 0] * N[node_shift + 1], &edgeN01[shift],
                  1);
      // edge12
      shift = ii * (P[1]);
      cblas_dcopy(P[1], L12, 1, &edgeN12[shift], 1);
      cblas_dscal(P[1], N[node_shift + 1] * N[node_shift + 2], &edgeN12[shift],
                  1);
      // edge20
      shift = ii * (P[2]);
      cblas_dcopy(P[2], L20, 1, &edgeN20[shift], 1);
      cblas_dscal(P[2], N[node_shift + 0] * N[node_shift + 2], &edgeN20[shift],
                  1);
    }
    if (diff_edgeN != NULL) {
      if (P[0] > 0) {
        // edge01
        shift = ii * (P[0]);
        bzero(&diff_edgeN01[2 * shift], sizeof(double) * 2 * (P[0]));
        // diffX
        cblas_daxpy(P[0], N[node_shift + 0] * N[node_shift + 1],
                    &diffL01[0 * (p[0] + 1)], 1, &diff_edgeN01[2 * shift + 0],
                    2);
        cblas_daxpy(P[0],
                    diffN[2 * 0 + 0] * N[node_shift + 1] +
                        N[node_shift + 0] * diffN[2 * 1 + 0],
                    L01, 1, &diff_edgeN01[2 * shift + 0], 2);
        // diff  Y
        cblas_daxpy(P[0], N[node_shift + 0] * N[node_shift + 1],
                    &diffL01[1 * (p[0] + 1)], 1, &diff_edgeN01[2 * shift + 1],
                    2);
        cblas_daxpy(P[0],
                    diffN[2 * 0 + 1] * N[node_shift + 1] +
                        N[node_shift + 0] * diffN[2 * 1 + 1],
                    L01, 1, &diff_edgeN01[2 * shift + 1], 2);
      }
      if (P[1] > 0) {
        // edge12
        shift = ii * (P[1]);
        bzero(&diff_edgeN12[2 * shift], sizeof(double) * 2 * (P[1]));
        // diffX
        cblas_daxpy(P[1], N[node_shift + 1] * N[node_shift + 2],
                    &diffL12[0 * (p[1] + 1)], 1, &diff_edgeN12[2 * shift + 0],
                    2);
        cblas_daxpy(P[1],
                    diffN[2 * 1 + 0] * N[node_shift + 2] +
                        N[node_shift + 1] * diffN[2 * 2 + 0],
                    L12, 1, &diff_edgeN12[2 * shift + 0], 2);
        // diffY
        cblas_daxpy(P[1], N[node_shift + 1] * N[node_shift + 2],
                    &diffL12[1 * (p[1] + 1)], 1, &diff_edgeN12[2 * shift + 1],
                    2);
        cblas_daxpy(P[1],
                    diffN[2 * 1 + 1] * N[node_shift + 2] +
                        N[node_shift + 1] * diffN[2 * 2 + 1],
                    L12, 1, &diff_edgeN12[2 * shift + 1], 2);
      }
      if (P[2] > 0) {
        // edge20
        shift = ii * (P[2]);
        bzero(&diff_edgeN20[2 * shift], sizeof(double) * 2 * (P[2]));
        // diffX
        cblas_daxpy(P[2], N[node_shift + 0] * N[node_shift + 2],
                    &diffL20[0 * (p[2] + 1)], 1, &diff_edgeN20[2 * shift + 0],
                    2);
        cblas_daxpy(P[2],
                    diffN[2 * 0 + 0] * N[node_shift + 2] +
                        N[node_shift + 0] * diffN[2 * 2 + 0],
                    L20, 1, &diff_edgeN20[2 * shift + 0], 2);
        // diffY
        cblas_daxpy(P[2], N[node_shift + 0] * N[node_shift + 2],
                    &diffL20[1 * (p[2] + 1)], 1, &diff_edgeN20[2 * shift + 1],
                    2);
        cblas_daxpy(P[2],
                    diffN[2 * 0 + 1] * N[node_shift + 2] +
                        N[node_shift + 0] * diffN[2 * 2 + 1],
                    L20, 1, &diff_edgeN20[2 * shift + 1], 2);
      }
    }
  }
  MoFEMFunctionReturnHot(0);
}
PetscErrorCode H1_FaceShapeFunctions_MBTRI(
    const int *face_nodes, int p, double *N, double *diffN, double *faceN,
    double *diff_faceN, int GDIM,
    PetscErrorCode (*base_polynomials)(int p, double s, double *diff_s,
                                       double *L, double *diffL,
                                       const int dim)) {
  MoFEMFunctionBeginHot;

  int P = NBFACETRI_H1(p);
  if (P == 0)
    MoFEMFunctionReturnHot(0);
  double diff_ksiL0[2], diff_ksiL1[2];
  double *diff_ksi_faces[] = {diff_ksiL0, diff_ksiL1};
  int dd = 0;
  if (diff_faceN != NULL) {
    for (; dd < 2; dd++) {
      diff_ksi_faces[0][dd] =
          diffN[face_nodes[1] * 2 + dd] - diffN[face_nodes[0] * 2 + dd];
      diff_ksi_faces[1][dd] =
          diffN[face_nodes[2] * 2 + dd] - diffN[face_nodes[0] * 2 + dd];
    }
  }
  int ii = 0;
  for (; ii < GDIM; ii++) {
    int node_shift = ii * 3;
    double ksi_faces[2];
    ksi_faces[0] =
        N[node_shift + face_nodes[1]] - N[node_shift + face_nodes[0]];
    ksi_faces[1] =
        N[node_shift + face_nodes[2]] - N[node_shift + face_nodes[0]];
    double L0[p + 1], L1[p + 1];
    double diffL0[2 * (p + 1)], diffL1[2 * (p + 1)];
    ierr = base_polynomials(p, ksi_faces[0], diff_ksi_faces[0], L0, diffL0, 2);
    CHKERRQ(ierr);
    ierr = base_polynomials(p, ksi_faces[1], diff_ksi_faces[1], L1, diffL1, 2);
    CHKERRQ(ierr);
    double v = N[node_shift + face_nodes[0]] * N[node_shift + face_nodes[1]] *
               N[node_shift + face_nodes[2]];
    double v2[2] = {0, 0};
    if (diff_faceN != NULL) {
      double n1n2 =
          N[node_shift + face_nodes[1]] * N[node_shift + face_nodes[2]];
      double n0n2 =
          N[node_shift + face_nodes[0]] * N[node_shift + face_nodes[2]];
      double n0n1 =
          N[node_shift + face_nodes[0]] * N[node_shift + face_nodes[1]];
      dd = 0;
      for (; dd < 2; dd++) {
        v2[dd] = diffN[face_nodes[0] * 2 + dd] * n1n2 +
                 diffN[face_nodes[1] * 2 + dd] * n0n2 +
                 diffN[face_nodes[2] * 2 + dd] * n0n1;
      }
    }
    int shift = ii * P;
    int jj = 0;
    int oo = 0;
    for (; oo <= (p - 3); oo++) {
      int pp0 = 0;
      for (; pp0 <= oo; pp0++) {
        int pp1 = oo - pp0;
        if (pp1 >= 0) {
          if (faceN != NULL) {
            faceN[shift + jj] = v * L0[pp0] * L1[pp1];
          }
          if (diff_faceN != NULL) {
            dd = 0;
            for (; dd < 2; dd++) {
              double *val = &diff_faceN[2 * shift + 2 * jj + dd];
              *val = (L0[pp0] * diffL1[dd * (p + 1) + pp1] +
                      diffL0[dd * (p + 1) + pp0] * L1[pp1]) *
                     v;
              *val += L0[pp0] * L1[pp1] * v2[dd];
            }
          }
          jj++;
        }
      }
    }
    if (jj != P)
      SETERRQ2(PETSC_COMM_SELF, 1, "wrong order %d != %d", jj, P);
  }
  MoFEMFunctionReturnHot(0);
}
PetscErrorCode H1_EdgeShapeFunctions_MBTET(
    int *sense, int *p, double *N, double *diffN, double *edgeN[],
    double *diff_edgeN[], int GDIM,
    PetscErrorCode (*base_polynomials)(int p, double s, double *diff_s,
                                       double *L, double *diffL,
                                       const int dim)) {
  MoFEMFunctionBeginHot;

  int P[6];
  int ee = 0;
  for (; ee < 6; ee++)
    P[ee] = NBEDGE_H1(p[ee]);
  int edges_nodes[2 * 6] = {0, 1, 1, 2, 2, 0, 0, 3, 1, 3, 2, 3};
  double diff_ksi01[3], diff_ksi12[3], diff_ksi20[3];
  double diff_ksi03[3], diff_ksi13[3], diff_ksi23[3];
  double *edges_diff_ksi[6] = {diff_ksi01, diff_ksi12, diff_ksi20,
                               diff_ksi03, diff_ksi13, diff_ksi23};
  for (ee = 0; ee < 6; ee++) {
    int dd = 0;
    for (; dd < 3; dd++) {
      edges_diff_ksi[ee][dd] = (diffN[edges_nodes[2 * ee + 1] * 3 + dd] -
                                diffN[edges_nodes[2 * ee + 0] * 3 + dd]) *
                               sense[ee];
    }
  }
  int ii = 0;
  for (; ii < GDIM; ii++) {
    int node_shift = ii * 4;
    double edges_ksi[6];
    for (ee = 0; ee < 6; ee++) {
      edges_ksi[ee] = (N[node_shift + edges_nodes[2 * ee + 1]] -
                       N[node_shift + edges_nodes[2 * ee + 0]]) *
                      sense[ee];
    }
    for (ee = 0; ee < 6; ee++) {
      if (P[ee] == 0)
        continue;
      double L[p[ee] + 1], diffL[3 * (p[ee] + 1)];
      ierr = base_polynomials(p[ee], edges_ksi[ee], edges_diff_ksi[ee], L,
                              diffL, 3);
      CHKERRQ(ierr);
      double v = N[node_shift + edges_nodes[2 * ee + 0]] *
                 N[node_shift + edges_nodes[2 * ee + 1]];
      if (edgeN != NULL)
        if (edgeN[ee] != NULL) {
          int shift = ii * P[ee];
          cblas_dcopy(P[ee], L, 1, &edgeN[ee][shift], 1);
          cblas_dscal(P[ee],
                      N[node_shift + edges_nodes[2 * ee + 0]] *
                          N[node_shift + edges_nodes[2 * ee + 1]],
                      &edgeN[ee][shift], 1);
        }
      if (diff_edgeN != NULL)
        if (diff_edgeN[ee] != NULL) {
          int shift = ii * P[ee];
          bzero(&diff_edgeN[ee][3 * shift], sizeof(double) * 3 * P[ee]);
          int dd = 0;
          for (; dd < 3; dd++) {
            cblas_daxpy(P[ee], v, &diffL[dd * (p[ee] + 1)], 1,
                        &diff_edgeN[ee][3 * shift + dd], 3);
            cblas_daxpy(P[ee],
                        diffN[3 * edges_nodes[2 * ee + 0] + dd] *
                                N[node_shift + edges_nodes[2 * ee + 1]] +
                            N[node_shift + edges_nodes[2 * ee + 0]] *
                                diffN[3 * edges_nodes[2 * ee + 1] + dd],
                        L, 1, &diff_edgeN[ee][3 * shift + dd], 3);
          }
        }
    }
  }
  MoFEMFunctionReturnHot(0);
}
PetscErrorCode H1_FaceShapeFunctions_MBTET(
    int *faces_nodes, int *p, double *N, double *diffN, double *faceN[],
    double *diff_faceN[], int GDIM,
    PetscErrorCode (*base_polynomials)(int p, double s, double *diff_s,
                                       double *L, double *diffL,
                                       const int dim)) {
  MoFEMFunctionBeginHot;

  int P[4];
  int ff = 0;
  for (; ff < 4; ff++)
    P[ff] = NBFACETRI_H1(p[ff]);
  double diff_ksiL0F0[3], diff_ksiL1F0[3];
  double diff_ksiL0F1[3], diff_ksiL1F1[3];
  double diff_ksiL0F2[3], diff_ksiL1F2[3];
  double diff_ksiL0F3[3], diff_ksiL1F3[3];
  double *diff_ksi_faces[] = {diff_ksiL0F0, diff_ksiL1F0, diff_ksiL0F1,
                              diff_ksiL1F1, diff_ksiL0F2, diff_ksiL1F2,
                              diff_ksiL0F3, diff_ksiL1F3};
  int dd = 0;
  for (; dd < 3; dd++) {
    for (ff = 0; ff < 4; ff++) {
      diff_ksi_faces[ff * 2 + 0][dd] =
          (diffN[faces_nodes[3 * ff + 1] * 3 + dd] -
           diffN[faces_nodes[3 * ff + 0] * 3 + dd]);
      diff_ksi_faces[ff * 2 + 1][dd] =
          (diffN[faces_nodes[3 * ff + 2] * 3 + dd] -
           diffN[faces_nodes[3 * ff + 0] * 3 + dd]);
    }
  }
  int ii = 0;
  for (; ii < GDIM; ii++) {
    int node_shift = ii * 4;
    double ksi_faces[8];
    for (ff = 0; ff < 4; ff++) {
      ksi_faces[2 * ff + 0] = N[node_shift + faces_nodes[3 * ff + 1]] -
                              N[node_shift + faces_nodes[3 * ff + 0]];
      ksi_faces[2 * ff + 1] = N[node_shift + faces_nodes[3 * ff + 2]] -
                              N[node_shift + faces_nodes[3 * ff + 0]];
    }
    int shift;
    for (ff = 0; ff < 4; ff++) {
      if (P[ff] == 0)
        continue;
      double L0[p[ff] + 1], L1[p[ff] + 1];
      double diffL0[3 * (p[ff] + 1)], diffL1[3 * (p[ff] + 1)];
      ierr = base_polynomials(p[ff], ksi_faces[ff * 2 + 0],
                              diff_ksi_faces[ff * 2 + 0], L0, diffL0, 3);
      CHKERRQ(ierr);
      ierr = base_polynomials(p[ff], ksi_faces[ff * 2 + 1],
                              diff_ksi_faces[ff * 2 + 1], L1, diffL1, 3);
      CHKERRQ(ierr);
      double v = N[node_shift + faces_nodes[3 * ff + 0]] *
                 N[node_shift + faces_nodes[3 * ff + 1]] *
                 N[node_shift + faces_nodes[3 * ff + 2]];
      double v2[3] = {0, 0, 0};
      dd = 0;
      double n1n2 = N[node_shift + faces_nodes[3 * ff + 1]] *
                    N[node_shift + faces_nodes[3 * ff + 2]];
      double n0n2 = N[node_shift + faces_nodes[3 * ff + 0]] *
                    N[node_shift + faces_nodes[3 * ff + 2]];
      double n0n1 = N[node_shift + faces_nodes[3 * ff + 0]] *
                    N[node_shift + faces_nodes[3 * ff + 1]];
      for (; dd < 3; dd++) {
        v2[dd] = diffN[3 * faces_nodes[3 * ff + 0] + dd] * n1n2 +
                 diffN[3 * faces_nodes[3 * ff + 1] + dd] * n0n2 +
                 diffN[3 * faces_nodes[3 * ff + 2] + dd] * n0n1;
      }
      shift = ii * P[ff];
      int jj = 0;
      int oo = 0;
      for (; oo <= (p[ff] - 3); oo++) {
        int pp0 = 0;
        for (; pp0 <= oo; pp0++) {
          int pp1 = oo - pp0;
          if (pp1 >= 0) {
            if (faceN != NULL)
              if (faceN[ff] != NULL) {
                faceN[ff][shift + jj] = v * L0[pp0] * L1[pp1];
              }
            if (diff_faceN != NULL)
              if (diff_faceN[ff] != NULL) {
                dd = 0;
                double L0L1 = L0[pp0] * L1[pp1];
                for (; dd < 3; dd++) {
                  diff_faceN[ff][3 * shift + 3 * jj + dd] =
                      (L0[pp0] * diffL1[dd * (p[ff] + 1) + pp1] +
                       diffL0[dd * (p[ff] + 1) + pp0] * L1[pp1]) *
                      v;
                  diff_faceN[ff][3 * shift + 3 * jj + dd] += L0L1 * v2[dd];
                }
              }
            jj++;
          }
        }
      }
      if (jj != P[ff])
        SETERRQ2(PETSC_COMM_SELF, 1, "wrong order %d != %d", jj, P[ff]);
    }
  }
  MoFEMFunctionReturnHot(0);
}
PetscErrorCode H1_VolumeShapeFunctions_MBTET(
    int p, double *N, double *diffN, double *volumeN, double *diff_volumeN,
    int GDIM,
    PetscErrorCode (*base_polynomials)(int p, double s, double *diff_s,
                                       double *L, double *diffL,
                                       const int dim)) {
  MoFEMFunctionBeginHot;

  int P = NBVOLUMETET_H1(p);
  if (P == 0)
    MoFEMFunctionReturnHot(0);
  double diff_ksiL0[3], diff_ksiL1[3], diff_ksiL2[3];
  int dd = 0;
  for (; dd < 3; dd++) {
    diff_ksiL0[dd] = (diffN[1 * 3 + dd] - diffN[0 * 3 + dd]);
    diff_ksiL1[dd] = (diffN[2 * 3 + dd] - diffN[0 * 3 + dd]);
    diff_ksiL2[dd] = (diffN[3 * 3 + dd] - diffN[0 * 3 + dd]);
  }
  int ii = 0;
  for (; ii < GDIM; ii++) {
    int node_shift = ii * 4;
    double ksiL0 = N[node_shift + 1] - N[node_shift + 0];
    double ksiL1 = N[node_shift + 2] - N[node_shift + 0];
    double ksiL2 = N[node_shift + 3] - N[node_shift + 0];
    double L0[p + 1], L1[p + 1], L2[p + 1];
    double diffL0[3 * (p + 1)], diffL1[3 * (p + 1)], diffL2[3 * (p + 1)];
    ierr = base_polynomials(p, ksiL0, diff_ksiL0, L0, diffL0, 3);
    CHKERRQ(ierr);
    ierr = base_polynomials(p, ksiL1, diff_ksiL1, L1, diffL1, 3);
    CHKERRQ(ierr);
    ierr = base_polynomials(p, ksiL2, diff_ksiL2, L2, diffL2, 3);
    CHKERRQ(ierr);
    double v = N[node_shift + 0] * N[node_shift + 1] * N[node_shift + 2] *
               N[node_shift + 3];
    double v2[3] = {0, 0, 0};
    dd = 0;
    for (; dd < 3; dd++) {
      v2[dd] = diffN[3 * 0 + dd] * N[node_shift + 1] * N[node_shift + 2] *
                   N[node_shift + 3] +
               N[node_shift + 0] * diffN[3 * 1 + dd] * N[node_shift + 2] *
                   N[node_shift + 3] +
               N[node_shift + 0] * N[node_shift + 1] * diffN[3 * 2 + dd] *
                   N[node_shift + 3] +
               N[node_shift + 0] * N[node_shift + 1] * N[node_shift + 2] *
                   diffN[3 * 3 + dd];
    }
    int shift = ii * P;
    int jj = 0;
    int oo = 0;
    for (; oo <= (p - 4); oo++) {
      int pp0 = 0;
      for (; pp0 <= oo; pp0++) {
        int pp1 = 0;
        for (; (pp0 + pp1) <= oo; pp1++) {
          int pp2 = oo - pp0 - pp1;
          if (pp2 >= 0) {
            if (volumeN != NULL) {
              volumeN[shift + jj] = L0[pp0] * L1[pp1] * L2[pp2] * v;
            }
            if (diff_volumeN != NULL) {
              dd = 0;
              for (; dd < 3; dd++) {
                diff_volumeN[3 * shift + 3 * jj + dd] =
                    (diffL0[dd * (p + 1) + pp0] * L1[pp1] * L2[pp2] +
                     L0[pp0] * diffL1[dd * (p + 1) + pp1] * L2[pp2] +
                     L0[pp0] * L1[pp1] * diffL2[dd * (p + 1) + pp2]) *
                    v;
                diff_volumeN[3 * shift + 3 * jj + dd] +=
                    L0[pp0] * L1[pp1] * L2[pp2] * v2[dd];
              }
            }
            jj++;
          }
        }
      }
    }
    if (jj != P)
      SETERRQ1(PETSC_COMM_SELF, 1, "wrong order %d", jj);
  }
  MoFEMFunctionReturnHot(0);
}
PetscErrorCode H1_EdgeShapeDiffMBTETinvJ(int *base_p, int *p,
                                         double *edge_diffN[], double *invJac,
                                         double *edge_diffNinvJac[], int GDIM) {
  MoFEMFunctionBeginHot;
  int ee = 0, ii, gg;
  for (; ee < 6; ee++) {
    for (ii = 0; ii < NBEDGE_H1(p[ee]); ii++) {
      for (gg = 0; gg < GDIM; gg++) {
        int shift1 = NBEDGE_H1(base_p[ee]) * gg;
        int shift2 = NBEDGE_H1(p[ee]) * gg;
        cblas_dgemv(CblasRowMajor, CblasTrans, 3, 3, 1., invJac, 3,
                    &(edge_diffN[ee])[3 * shift1 + 3 * ii], 1, 0.,
                    &(edge_diffNinvJac[ee])[3 * shift2 + 3 * ii], 1);
      }
    }
  }
  MoFEMFunctionReturnHot(0);
}
PetscErrorCode H1_FaceShapeDiffMBTETinvJ(int *base_p, int *p,
                                         double *face_diffN[], double *invJac,
                                         double *face_diffNinvJac[], int GDIM) {
  MoFEMFunctionBeginHot;
  int ff = 0, ii, gg;
  for (; ff < 4; ff++) {
    for (ii = 0; ii < NBFACETRI_H1(p[ff]); ii++) {
      for (gg = 0; gg < GDIM; gg++) {
        int shift1 = NBFACETRI_H1(base_p[ff]) * gg;
        int shift2 = NBFACETRI_H1(p[ff]) * gg;
        cblas_dgemv(CblasRowMajor, CblasTrans, 3, 3, 1., invJac, 3,
                    &(face_diffN[ff])[3 * shift1 + 3 * ii], 1, 0.,
                    &(face_diffNinvJac[ff])[3 * shift2 + 3 * ii], 1);
      }
    }
  }
  MoFEMFunctionReturnHot(0);
}
PetscErrorCode H1_VolumeShapeDiffMBTETinvJ(int base_p, int p,
                                           double *volume_diffN, double *invJac,
                                           double *volume_diffNinvJac,
                                           int GDIM) {
  MoFEMFunctionBeginHot;
  int ii, gg;
  for (ii = 0; ii < NBVOLUMETET_H1(p); ii++) {
    for (gg = 0; gg < GDIM; gg++) {
      int shift1 = NBVOLUMETET_H1(base_p) * gg;
      int shift2 = NBVOLUMETET_H1(p) * gg;
      cblas_dgemv(CblasRowMajor, CblasTrans, 3, 3, 1., invJac, 3,
                  &(volume_diffN)[3 * shift1 + 3 * ii], 1, 0.,
                  &(volume_diffNinvJac)[3 * shift2 + 3 * ii], 1);
    }
  }
  MoFEMFunctionReturnHot(0);
}
PetscErrorCode H1_EdgeGradientOfDeformation_hierarchical(int p, double *diffN,
                                                        double *dofs,
                                                        double *F) {
  MoFEMFunctionBeginHot;
  int col, row = 0;
  for (; row < 3; row++)
    for (col = 0; col < 3; col++)
      F[3 * row + col] =
          cblas_ddot(NBEDGE_H1(p), &diffN[col], 3, &dofs[row], 3);
  MoFEMFunctionReturnHot(0);
}
PetscErrorCode H1_FaceGradientOfDeformation_hierarchical(int p, double *diffN,
                                                        double *dofs,
                                                        double *F) {
  MoFEMFunctionBeginHot;
  int col, row = 0;
  for (; row < 3; row++)
    for (col = 0; col < 3; col++)
      F[3 * row + col] =
          cblas_ddot(NBFACETRI_H1(p), &diffN[col], 3, &dofs[row], 3);
  MoFEMFunctionReturnHot(0);
}
PetscErrorCode H1_VolumeGradientOfDeformation_hierarchical(int p, double *diffN,
                                                          double *dofs,
                                                          double *F) {
  MoFEMFunctionBeginHot;
  int col, row = 0;
  for (; row < 3; row++)
    for (col = 0; col < 3; col++)
      F[3 * row + col] =
          cblas_ddot(NBVOLUMETET_H1(p), &diffN[col], 3, &dofs[row], 3);
  MoFEMFunctionReturnHot(0);
}
PetscErrorCode H1_QuadShapeFunctions_MBPRISM(
    int *faces_nodes, int *p, double *N, double *diffN, double *faceN[],
    double *diff_faceN[], int GDIM,
    PetscErrorCode (*base_polynomials)(int p, double s, double *diff_s,
                                       double *L, double *diffL,
                                       const int dim)) {
  MoFEMFunctionBeginHot;
  // TODO: return separately components of the tensor product between two edges 

  int P[3];
  int ff = 0;
  for (; ff < 3; ff++) {
    P[ff] = NBFACEQUAD_H1(p[ff]);
  }
  double ksi_faces[6];
  double diff_ksiL0F0[3], diff_ksiL3F0[3];
  double diff_ksiL0F1[3], diff_ksiL3F1[3];
  double diff_ksiL0F2[3], diff_ksiL3F2[3];
  double *diff_ksi_faces[] = {diff_ksiL0F0, diff_ksiL3F0, diff_ksiL0F1,
                              diff_ksiL3F1, diff_ksiL0F2, diff_ksiL3F2};
  int ii = 0;
  for (; ii < GDIM; ii++) {
    int node_shift = ii * 6;
    int node_diff_shift = 3 * node_shift;
    int ff = 0;
    for (; ff < 3; ff++) {
      if (P[ff] == 0)
        continue;
      int n0 = faces_nodes[4 * ff + 0];
      int n1 = faces_nodes[4 * ff + 1];
      int n2 = faces_nodes[4 * ff + 2];
      int n3 = faces_nodes[4 * ff + 3];
      int e0 = 2 * ff + 0;
      int e1 = 2 * ff + 1;
      ksi_faces[e0] = N[node_shift + n1] + N[node_shift + n2] -
                      N[node_shift + n0] - N[node_shift + n3];
      ksi_faces[e1] = N[node_shift + n2] + N[node_shift + n3] -
                      N[node_shift + n0] - N[node_shift + n1];
      double eps = 1e-14;
      if (ksi_faces[e0] < -(1 + eps) || ksi_faces[e0] > (1 + eps)) {
        SETERRQ(PETSC_COMM_SELF, 1, "Data inconsistency");
      }
      if (ksi_faces[e1] < -(1 + eps) || ksi_faces[e1] > (1 + eps)) {
        SETERRQ(PETSC_COMM_SELF, 1, "Data inconsistency");
      }
      int dd = 0;
      for (; dd < 3; dd++) {
        diff_ksi_faces[e0][dd] = diffN[node_diff_shift + 2 * n1 + dd] +
                                 diffN[node_diff_shift + 2 * n2 + dd] -
                                 diffN[node_diff_shift + 2 * n0 + dd] -
                                 diffN[node_diff_shift + 2 * n3 + dd];
        diff_ksi_faces[e1][dd] = diffN[node_diff_shift + 2 * n2 + dd] +
                                 diffN[node_diff_shift + 2 * n3 + dd] -
                                 diffN[node_diff_shift + 2 * n0 + dd] -
                                 diffN[node_diff_shift + 2 * n1 + dd];
      }
      double L0[p[ff] + 1], L1[p[ff] + 1];
      double diffL0[3 * (p[ff] + 1)], diffL1[3 * (p[ff] + 1)];
      ierr = base_polynomials(p[ff], ksi_faces[e0], diff_ksi_faces[e0], L0,
                              diffL0, 3);
      CHKERRQ(ierr);
      ierr = base_polynomials(p[ff], ksi_faces[e1], diff_ksi_faces[e1], L1,
                              diffL1, 3);
      CHKERRQ(ierr);
      double v = N[node_shift + n0] * N[node_shift + n2];
      double v2[3] = {0, 0, 0};
      dd = 0;
      for (; dd < 3; dd++) {
        v2[dd] = diffN[node_diff_shift + 3 * n0 + dd] * N[node_shift + n2] +
                 N[node_shift + n0] * diffN[node_diff_shift + 3 * n2 + dd];
      }
      int shift;
      shift = ii * P[ff];
      int jj = 0;
      int oo = 0;
      for (; oo <= (p[ff] - 4); oo++) {
        int pp0 = 0;
        for (; pp0 <= oo; pp0++) {
          int pp1 = oo - pp0;
          if (pp1 >= 0) {
            if (faceN != NULL) {
              if (faceN[ff] != NULL) {
                faceN[ff][shift + jj] = v * L0[pp0] * L1[pp1];
              }
            }
            if (diff_faceN != NULL) {
              if (diff_faceN[ff] != NULL) {
                dd = 0;
                for (; dd < 3; dd++) {
                  diff_faceN[ff][3 * shift + 3 * jj + dd] =
                      (L0[pp0] * diffL1[dd * (p[ff] + 1) + pp1] +
                       diffL0[dd * (p[ff] + 1) + pp0] * L1[pp1]) *
                      v;
                  diff_faceN[ff][3 * shift + 3 * jj + dd] +=
                      L0[pp0] * L1[pp1] * v2[dd];
                }
              }
            }
            jj++;
          }
        }
      }
      if (jj != P[ff])
        SETERRQ2(PETSC_COMM_SELF, 1, "wrong order %d != %d", jj, P[ff]);
    }
  }
  MoFEMFunctionReturnHot(0);
}
PetscErrorCode H1_VolumeShapeFunctions_MBPRISM(
    int p, double *N, double *diffN, double *volumeN, double *diff_volumeN,
    int GDIM,
    PetscErrorCode (*base_polynomials)(int p, double s, double *diff_s,
                                       double *L, double *diffL,
                                       const int dim)) {
  MoFEMFunctionBeginHot;

  int P = NBVOLUMEPRISM_H1(p);
  if (P == 0)
    MoFEMFunctionReturnHot(0);
  double diff_ksiL0[3], diff_ksiL1[3], diff_ksiL2[3];
  int ii = 0;
  for (; ii < GDIM; ii++) {
    int node_shift = ii * 6;
    int node_diff_shift = ii * 18;
    double ksiL0 = N[node_shift + 1] - N[node_shift + 0];
    double ksiL1 = N[node_shift + 2] - N[node_shift + 0];
    double ksiL2 =
        N[node_shift + 3] - N[node_shift + 0]; // 2*gauss_pts[2*GDIM+ii]-1.;
    int dd = 0;
    for (; dd < 3; dd++) {
      diff_ksiL0[dd] = (diffN[node_diff_shift + 1 * 3 + dd] -
                        diffN[node_diff_shift + 0 * 3 + dd]);
      diff_ksiL1[dd] = (diffN[node_diff_shift + 2 * 3 + dd] -
                        diffN[node_diff_shift + 0 * 3 + dd]);
      diff_ksiL2[dd] = 2;
    }
    double L0[p + 1], L1[p + 1], L2[p + 1];
    double diffL0[3 * (p + 1)], diffL1[3 * (p + 1)], diffL2[3 * (p + 1)];
    ierr = base_polynomials(p, ksiL0, diff_ksiL0, L0, diffL0, 3);
    CHKERRQ(ierr);
    ierr = base_polynomials(p, ksiL1, diff_ksiL1, L1, diffL1, 3);
    CHKERRQ(ierr);
    ierr = base_polynomials(p, ksiL2, diff_ksiL2, L2, diffL2, 3);
    CHKERRQ(ierr);
    double t0 = (N[node_shift + 2] + N[node_shift + 5]);
    double v = N[node_shift + 0] * N[node_shift + 4] * t0;
    double v2[3] = {0, 0, 0};
    dd = 0;
    for (; dd < 3; dd++) {
      v2[dd] = diffN[node_diff_shift + 3 * 0 + dd] * N[node_shift + 4] * t0 +
               N[node_shift + 0] * diffN[node_diff_shift + 3 * 4 + dd] * t0 +
               N[node_shift + 0] * N[node_shift + 4] *
                   (diffN[node_diff_shift + 3 * 2 + dd] +
                    diffN[node_diff_shift + 3 * 5 + dd]);
    }
    int shift = ii * P;
    int jj = 0;
    int oo = 0;
    for (; oo <= (p - 6); oo++) {
      int pp0 = 0;
      for (; pp0 <= oo; pp0++) {
        int pp1 = 0;
        for (; (pp0 + pp1) <= oo; pp1++) {
          int pp2 = oo - pp0 - pp1;
          if (pp2 >= 0) {
            if (volumeN != NULL) {
              volumeN[shift + jj] = L0[pp0] * L1[pp1] * L2[pp2] * v;
            }
            if (diff_volumeN != NULL) {
              dd = 0;
              for (; dd < 3; dd++) {
                diff_volumeN[3 * shift + 3 * jj + dd] =
                    (diffL0[dd * (p + 1) + pp0] * L1[pp1] * L2[pp2] +
                     L0[pp0] * diffL1[dd * (p + 1) + pp1] * L2[pp2] +
                     L0[pp0] * L1[pp1] * diffL2[dd * (p + 1) + pp2]) *
                    v;
                diff_volumeN[3 * shift + 3 * jj + dd] +=
                    L0[pp0] * L1[pp1] * L2[pp2] * v2[dd];
              }
            }
            jj++;
          }
        }
      }
    }
    if (jj != P)
      SETERRQ1(PETSC_COMM_SELF, 1, "wrong order %d", jj);
  }
  MoFEMFunctionReturnHot(0);
}
PetscErrorCode H1_QuadShapeFunctions_MBQUAD(
    int *faces_nodes, int p, double *N, double *diffN, double *faceN,
    double *diff_faceN, int GDIM,
    PetscErrorCode (*base_polynomials)(int p, double s, double *diff_s,
                                       double *L, double *diffL,
                                       const int dim)) {
  MoFEMFunctionBeginHot;
  // TODO: return separately components of the tensor product between two edges

  int P = NBFACEQUAD_H1(p);
  if (P == 0)
    MoFEMFunctionReturnHot(0);
  double ksi_faces[2];
  double diff_ksiL0F0[2], diff_ksiL1F0[2];
  double *diff_ksi_faces[] = {diff_ksiL0F0, diff_ksiL1F0}; 
  int ii = 0;
  for (; ii < GDIM; ii++) {
    int node_shift = ii * 4;
    int node_diff_shift = 2 * node_shift;

    int n0 = faces_nodes[0];
    int n1 = faces_nodes[1];
    int n2 = faces_nodes[2];
    int n3 = faces_nodes[3];
    int e0 = 0;
    int e1 = 1;
    ksi_faces[e0] = N[node_shift + n1] + N[node_shift + n2] -
                    N[node_shift + n0] - N[node_shift + n3];
    ksi_faces[e1] = N[node_shift + n2] + N[node_shift + n3] -
                    N[node_shift + n0] - N[node_shift + n1];
    double eps = 1e-14;
    if (ksi_faces[e0] < -(1 + eps) || ksi_faces[e0] > (1 + eps)) {
      SETERRQ(PETSC_COMM_SELF, 1, "Data inconsistency");
    }
    if (ksi_faces[e1] < -(1 + eps) || ksi_faces[e1] > (1 + eps)) {
      SETERRQ(PETSC_COMM_SELF, 1, "Data inconsistency");
    }
    int dd = 0;
    for (; dd < 2; dd++) {
      diff_ksi_faces[e0][dd] = diffN[node_diff_shift + 2 * n1 + dd] +
                               diffN[node_diff_shift + 2 * n2 + dd] -
                               diffN[node_diff_shift + 2 * n0 + dd] -
                               diffN[node_diff_shift + 2 * n3 + dd];
      diff_ksi_faces[e1][dd] = diffN[node_diff_shift + 2 * n2 + dd] +
                               diffN[node_diff_shift + 2 * n3 + dd] -
                               diffN[node_diff_shift + 2 * n0 + dd] -
                               diffN[node_diff_shift + 2 * n1 + dd];
    }
    double L0[p + 1], L1[p + 1];
    double diffL0[2 * (p + 1)], diffL1[2 * (p + 1)]; 
    ierr =
        base_polynomials(p, ksi_faces[e0], diff_ksi_faces[e0], L0, diffL0, 2); 
    CHKERRQ(ierr);
    ierr =
        base_polynomials(p, ksi_faces[e1], diff_ksi_faces[e1], L1, diffL1, 2); 
    CHKERRQ(ierr);
    double v = N[node_shift + n0] * N[node_shift + n2];
    double v2[] = {0, 0};
    dd = 0;
    for (; dd < 2; dd++) {
      v2[dd] = diffN[node_diff_shift + 2 * n0 + dd] * N[node_shift + n2] +
               N[node_shift + n0] * diffN[node_diff_shift + 2 * n2 + dd];
    }
    int shift;
    shift = ii * P;
    int jj = 0;
    int oo = 0;
    for (; oo <= (p - 4); oo++) {
      int pp0 = 0;
      for (; pp0 <= oo; pp0++) {
        int pp1 = oo - pp0;
        if (pp1 >= 0) {
          if (faceN != NULL) {
            faceN[shift + jj] = L0[pp0] * L1[pp1] * v;
          }
          if (diff_faceN != NULL) {
            dd = 0;
            for (; dd < 2; dd++) {
              diff_faceN[2 * shift + 2 * jj + dd] =
                  (L0[pp0] * diffL1[dd * (p + 1) + pp1] +
                   diffL0[dd * (p + 1) + pp0] * L1[pp1]) *
                  v;
              diff_faceN[2 * shift + 2 * jj + dd] += L0[pp0] * L1[pp1] * v2[dd];
            }
          }
          jj++;
        }
      }
    }
    if (jj != P)
      SETERRQ2(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY, "wrong order %d != %d", jj, P);
  }
  MoFEMFunctionReturnHot(0);
}

PetscErrorCode H1_EdgeShapeFunctions_MBQUAD(
    int *sense, int *p, double *N, double *diffN, double *edgeN[4],
    double *diff_edgeN[4], int GDIM,
    PetscErrorCode (*base_polynomials)(int p, double s, double *diff_s,
                                       double *L, double *diffL,
                                       const int dim)) {
  MoFEMFunctionBeginHot;

  double *edgeN01 = NULL, *edgeN12 = NULL, *edgeN23 = NULL, *edgeN30 = NULL;
  if (edgeN != NULL) {
    edgeN01 = edgeN[0];
    edgeN12 = edgeN[1];
    edgeN23 = edgeN[2];
    edgeN30 = edgeN[3];
  }
  double *diff_edgeN01 = NULL, *diff_edgeN12 = NULL, *diff_edgeN23 = NULL,
         *diff_edgeN30 = NULL;
  if (diff_edgeN != NULL) {
    diff_edgeN01 = diff_edgeN[0];
    diff_edgeN12 = diff_edgeN[1];
    diff_edgeN23 = diff_edgeN[2];
    diff_edgeN30 = diff_edgeN[3];
  }
  int P[4];
  int ee = 0;
  for (; ee < 4; ee++) 
    P[ee] = NBEDGE_H1(p[ee]);
  
  int n0 = 0;
  int n1 = 1;
  int n2 = 2;
  int n3 = 3;

  int ii = 0;
  for (; ii < GDIM; ii++) {
    int node_shift = ii * 4;
    int node_diff_shift = 2 * node_shift;
    double ksi01 = (N[node_shift + n1] - N[node_shift + n0]) * sense[n0];
    double ksi12 = (N[node_shift + n2] - N[node_shift + n1]) * sense[n1];
    double ksi23 = (N[node_shift + n3] - N[node_shift + n2]) * sense[n2];
    double ksi30 = (N[node_shift + n0] - N[node_shift + n3]) * sense[n3];

    double diff_ksi01[2], diff_ksi12[2], diff_ksi23[2], diff_ksi30[2];
    int dd = 0;
    for (; dd < 2; dd++) {
      diff_ksi01[dd] = (diffN[node_diff_shift + 2 * n1 + dd] -
                        diffN[node_diff_shift + 2 * n0 + dd]) *
                       sense[n0];
      diff_ksi12[dd] = (diffN[node_diff_shift + 2 * n2 + dd] -
                        diffN[node_diff_shift + 2 * n1 + dd]) *
                       sense[n1];
      diff_ksi23[dd] = (diffN[node_diff_shift + 2 * n3 + dd] -
                        diffN[node_diff_shift + 2 * n2 + dd]) *
                       sense[n2];
      diff_ksi30[dd] = (diffN[node_diff_shift + 2 * n0 + dd] -
                        diffN[node_diff_shift + 2 * n3 + dd]) *
                       sense[n3];
    }

    double L01[p[0] + 1], L12[p[1] + 1], L23[p[2] + 1], L30[p[3] + 1];
    double diffL01[2 * (p[0] + 1)], diffL12[2 * (p[1] + 1)],
        diffL23[2 * (p[2] + 1)], diffL30[2 * (p[3] + 1)];
    ierr = base_polynomials(p[0], ksi01, diff_ksi01, L01, diffL01, 2);
    CHKERRQ(ierr);
    ierr = base_polynomials(p[1], ksi12, diff_ksi12, L12, diffL12, 2);
    CHKERRQ(ierr);
    ierr = base_polynomials(p[2], ksi23, diff_ksi23, L23, diffL23, 2);
    CHKERRQ(ierr);
    ierr = base_polynomials(p[3], ksi30, diff_ksi30, L30, diffL30, 2);
    CHKERRQ(ierr);

    int shift;
    if (edgeN != NULL) {
      // edge01
      shift = ii * (P[0]);
      cblas_dcopy(P[0], L01, 1, &edgeN01[shift], 1);
      cblas_dscal(P[0], N[node_shift + n0] * N[node_shift + n1],
                  &edgeN01[shift], 1);
      // edge12
      shift = ii * (P[1]);
      cblas_dcopy(P[1], L12, 1, &edgeN12[shift], 1);
      cblas_dscal(P[1], N[node_shift + n1] * N[node_shift + n2],
                  &edgeN12[shift], 1);
      // edge23
      shift = ii * (P[2]);
      cblas_dcopy(P[2], L23, 1, &edgeN23[shift], 1);
      cblas_dscal(P[2], N[node_shift + n2] * N[node_shift + n3],
                  &edgeN23[shift], 1);
      // edge30
      shift = ii * (P[3]);
      cblas_dcopy(P[3], L30, 1, &edgeN30[shift], 1);
      cblas_dscal(P[3], N[node_shift + n3] * N[node_shift + n0],
                  &edgeN30[shift], 1);
    }
    if (diff_edgeN != NULL) {
      if (P[0] > 0) {
        // edge01
        shift = ii * (P[0]);

        bzero(&diff_edgeN01[2 * shift], sizeof(double) * 2 * (P[0]));
        // diffX
        cblas_daxpy(P[0], N[node_shift + n0] * N[node_shift + n1],
                    &diffL01[0 * (p[0] + 1)], 1, &diff_edgeN01[2 * shift + 0],
                    2);
        cblas_daxpy(P[0],
                    diffN[node_diff_shift + 2 * n0 + 0] * N[node_shift + n1] +
                        N[node_shift + n0] *
                            diffN[node_diff_shift + 2 * n1 + 0],
                    L01, 1, &diff_edgeN01[2 * shift + 0], 2);
        // diffY
        cblas_daxpy(P[0], N[node_shift + n0] * N[node_shift + n1],
                    &diffL01[1 * (p[0] + 1)], 1, &diff_edgeN01[2 * shift + 1],
                    2);
        cblas_daxpy(P[0],
                    diffN[node_diff_shift + 2 * n0 + 1] * N[node_shift + n1] +
                        N[node_shift + n0] *
                            diffN[node_diff_shift + 2 * n1 + 1],
                    L01, 1, &diff_edgeN01[2 * shift + 1], 2);
      }
      if (P[1] > 0) {
        // edge12
        shift = ii * (P[1]);
        bzero(&diff_edgeN12[2 * shift], sizeof(double) * 2 * (P[1]));
        // diffX
        cblas_daxpy(P[1], N[node_shift + n1] * N[node_shift + n2],
                    &diffL12[0 * (p[1] + 1)], 1, &diff_edgeN12[2 * shift + 0],
                    2);
        cblas_daxpy(P[1],
                    diffN[node_diff_shift + 2 * n1 + 0] * N[node_shift + n2] +
                        N[node_shift + n1] *
                            diffN[node_diff_shift + 2 * n2 + 0],
                    L12, 1, &diff_edgeN12[2 * shift + 0], 2);
        // diffY
        cblas_daxpy(P[1], N[node_shift + n1] * N[node_shift + n2],
                    &diffL12[1 * (p[1] + 1)], 1, &diff_edgeN12[2 * shift + 1],
                    2);
        cblas_daxpy(P[1],
                    diffN[node_diff_shift + 2 * n1 + 1] * N[node_shift + n2] +
                        N[node_shift + n1] *
                            diffN[node_diff_shift + 2 * n2 + 1],
                    L12, 1, &diff_edgeN12[2 * shift + 1], 2);
      }
      if (P[2] > 0) {
        // edge23
        shift = ii * (P[2]);
        bzero(&diff_edgeN23[2 * shift], sizeof(double) * 2 * (P[2]));
        // diffX
        cblas_daxpy(P[2], N[node_shift + n2] * N[node_shift + n3],
                    &diffL23[0 * (p[2] + 1)], 1, &diff_edgeN23[2 * shift + 0],
                    2);
        cblas_daxpy(P[2],
                    diffN[node_diff_shift + 2 * n2 + 0] * N[node_shift + n3] +
                        N[node_shift + n2] *
                            diffN[node_diff_shift + 2 * n3 + 0],
                    L23, 1, &diff_edgeN23[2 * shift + 0], 2);
        // diffY
        cblas_daxpy(P[2], N[node_shift + n2] * N[node_shift + n3],
                    &diffL23[1 * (p[2] + 1)], 1, &diff_edgeN23[2 * shift + 1],
                    2);
        cblas_daxpy(P[2],
                    diffN[node_diff_shift + 2 * n2 + 1] * N[node_shift + n3] +
                        N[node_shift + n2] *
                            diffN[node_diff_shift + 2 * n3 + 1],
                    L23, 1, &diff_edgeN23[2 * shift + 1], 2);
      }
      if (P[3] > 0) {
        // edge30
        shift = ii * (P[3]);
        bzero(&diff_edgeN30[2 * shift], sizeof(double) * 2 * (P[3]));
        // diffX
        cblas_daxpy(P[3], N[node_shift + n3] * N[node_shift + n0],
                    &diffL30[0 * (p[3] + 1)], 1, &diff_edgeN30[2 * shift + 0],
                    2);
        cblas_daxpy(P[3],
                    diffN[node_diff_shift + 2 * n3 + 0] * N[node_shift + n0] +
                        N[node_shift + n3] *
                            diffN[node_diff_shift + 2 * n0 + 0],
                    L30, 1, &diff_edgeN30[2 * shift + 0], 2);
        // diffY
        cblas_daxpy(P[3], N[node_shift + n3] * N[node_shift + n0],
                    &diffL30[1 * (p[3] + 1)], 1, &diff_edgeN30[2 * shift + 1],
                    2);
        cblas_daxpy(P[3],
                    diffN[node_diff_shift + 2 * n3 + 1] * N[node_shift + n0] +
                        N[node_shift + n3] *
                            diffN[node_diff_shift + 2 * n0 + 1],
                    L30, 1, &diff_edgeN30[2 * shift + 1], 2);
      }
    }
  }
  MoFEMFunctionReturnHot(0);
}

PetscErrorCode H1_EdgeGradientOfDeformation_hierachical(int p, double *diffN,
                                                        double *dofs,
                                                        double *F) {
  return H1_EdgeGradientOfDeformation_hierarchical(p, diffN, dofs, F);
}
PetscErrorCode H1_FaceGradientOfDeformation_hierachical(int p, double *diffN,
                                                        double *dofs,
                                                        double *F) {
  return H1_FaceGradientOfDeformation_hierarchical(p, diffN, dofs, F);
}
PetscErrorCode H1_VolumeGradientOfDeformation_hierachical(int p, double *diffN,
                                                          double *dofs,
                                                          double *F) {
  return H1_VolumeGradientOfDeformation_hierarchical(p, diffN, dofs, F);
}
