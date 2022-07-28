

// FIXME This is obsolete implementation, need to be rewritten

#ifdef __cplusplus
extern "C" {
#endif

#ifndef __lapack__
#define __lapack__

#define __CLPK_TYPES__
#ifdef __CLPK_TYPES__

#if defined(__LP64__) /* In LP64 match sizes with the 32 bit ABI */
typedef int __CLPK_integer;
typedef int __CLPK_logical;
typedef float __CLPK_real;
typedef double __CLPK_doublereal;
typedef __CLPK_logical (*__CLPK_L_fp)();
typedef int __CLPK_ftnlen;
#else
typedef long int __CLPK_integer;
typedef long int __CLPK_logical;
typedef float __CLPK_real;
typedef double __CLPK_doublereal;
typedef __CLPK_logical (*__CLPK_L_fp)();
typedef long int __CLPK_ftnlen;
#endif

typedef struct {
  __CLPK_real r, i;
} __CLPK_complex;
typedef struct {
  __CLPK_doublereal r, i;
} __CLPK_doublecomplex;

__CLPK_integer sgetrf_(__CLPK_integer *m, __CLPK_integer *n, float *a,
                       __CLPK_integer *lda, __CLPK_integer *ipiv,
                       __CLPK_integer *info);
#ifndef MOAB_dgetrf
__CLPK_integer dgetrf_(__CLPK_integer *m, __CLPK_integer *n,
                       __CLPK_doublereal *a, __CLPK_integer *lda,
                       __CLPK_integer *ipiv, __CLPK_integer *info);
#endif // MOAB_dgetrf
__CLPK_integer dgetrs_(char *trans, __CLPK_integer *n, __CLPK_integer *nrhs,
                       __CLPK_doublereal *a, __CLPK_integer *lda,
                       __CLPK_integer *ipiv, __CLPK_doublereal *b,
                       __CLPK_integer *ldb, __CLPK_integer *info);
__CLPK_integer dgesv_(__CLPK_integer *n, __CLPK_integer *nrhs,
                      __CLPK_doublereal *a, __CLPK_integer *lda,
                      __CLPK_integer *ipiv, __CLPK_doublereal *b,
                      __CLPK_integer *ldb, __CLPK_integer *info);
#ifndef MOAB_dgetri
__CLPK_integer dgetri_(__CLPK_integer *n, __CLPK_doublereal *a,
                       __CLPK_integer *lda, __CLPK_integer *ipiv,
                       __CLPK_doublereal *work, __CLPK_integer *lwork,
                       __CLPK_integer *info);
#endif // MOAB_dgetri
__CLPK_integer dpotrf_(char *uplo, __CLPK_integer *n, __CLPK_doublereal *a,
                       __CLPK_integer *lda, __CLPK_integer *info);
__CLPK_integer dpotrs_(char *uplo, __CLPK_integer *n, __CLPK_integer *nrhs,
                       __CLPK_doublereal *a, __CLPK_integer *lda,
                       __CLPK_doublereal *b, __CLPK_integer *ldb,
                       __CLPK_integer *info);
__CLPK_integer dposv_(char *uplo, __CLPK_integer *n, __CLPK_integer *nrhs,
                      __CLPK_doublereal *a, __CLPK_integer *lda,
                      __CLPK_doublereal *b, __CLPK_integer *ldb,
                      __CLPK_integer *info);
__CLPK_integer dsysv_(char *uplo, __CLPK_integer *n, __CLPK_integer *nrhs,
                      __CLPK_doublereal *a, __CLPK_integer *lda,
                      __CLPK_integer *ipiv, __CLPK_doublereal *b,
                      __CLPK_integer *ldb, __CLPK_doublereal *work,
                      __CLPK_integer *lwork, __CLPK_integer *info);
__CLPK_integer dpotri_(char *uplo, __CLPK_integer *n, __CLPK_doublereal *a,
                       __CLPK_integer *lda, __CLPK_integer *info);
__CLPK_integer dgesvd_(char *jobu, char *jobvt, __CLPK_integer *m,
                       __CLPK_integer *n, __CLPK_doublereal *a,
                       __CLPK_integer *lda, __CLPK_doublereal *s,
                       __CLPK_doublereal *u, __CLPK_integer *ldu,
                       __CLPK_doublereal *vt, __CLPK_integer *ldvt,
                       __CLPK_doublereal *work, __CLPK_integer *lwork,
                       __CLPK_integer *info);
__CLPK_integer sgesvd_(char *jobu, char *jobvt, __CLPK_integer *m,
                       __CLPK_integer *n, float *a, __CLPK_integer *lda,
                       float *s, float *u, __CLPK_integer *ldu, float *vt,
                       __CLPK_integer *ldvt, float *work, __CLPK_integer *lwork,
                       __CLPK_integer *info);
__CLPK_integer dsyev_(char *jobz, char *uplo, __CLPK_integer *n,
                      __CLPK_doublereal *a, __CLPK_integer *lda,
                      __CLPK_doublereal *w, __CLPK_doublereal *work,
                      __CLPK_integer *lwork, __CLPK_integer *info);
__CLPK_integer zheev_(char *jobz, char *uplo, __CLPK_integer *n,
                      __CLPK_doublecomplex *a, __CLPK_integer *lda,
                      __CLPK_doublereal *w, __CLPK_doublecomplex *work,
                      __CLPK_integer *lwork, __CLPK_doublereal *rwork,
                      __CLPK_integer *info);
__CLPK_integer zgeev_(char *jobvl, char *jobvr, __CLPK_integer *n,
                      __CLPK_doublecomplex *a, __CLPK_integer *lda,
                      __CLPK_doublecomplex *w, __CLPK_doublecomplex *vl,
                      __CLPK_integer *ldvl, __CLPK_doublecomplex *vr,
                      __CLPK_integer *ldvr, __CLPK_doublecomplex *work,
                      __CLPK_integer *lwork, __CLPK_doublereal *rwork,
                      __CLPK_integer *info);
__CLPK_integer dgelsy_(char *trans, __CLPK_integer *m, __CLPK_integer *n,
                       __CLPK_integer *nrhs, __CLPK_doublereal *a,
                       __CLPK_integer *lda, __CLPK_doublereal *b,
                       __CLPK_integer *ldb, __CLPK_doublereal *work,
                       __CLPK_integer *lwork, __CLPK_integer *info);
__CLPK_integer dgels_(char *trans, __CLPK_integer *m, __CLPK_integer *n,
                      __CLPK_integer *nrhs, __CLPK_doublereal *a,
                      __CLPK_integer *lda, __CLPK_doublereal *b,
                      __CLPK_integer *ldb, __CLPK_doublereal *work,
                      __CLPK_integer *lwork, __CLPK_integer *info);
__CLPK_integer dgesdd_(char *jobz, __CLPK_integer *m, __CLPK_integer *n,
                       __CLPK_doublereal *a, __CLPK_integer *lda,
                       __CLPK_doublereal *s, __CLPK_doublereal *u,
                       __CLPK_integer *ldu, __CLPK_doublereal *vt,
                       __CLPK_integer *ldvt, __CLPK_doublereal *work,
                       __CLPK_integer *lwork, __CLPK_integer *iwork,
                       __CLPK_integer *info);
#ifndef MOAB_dsyevd
__CLPK_integer dsyevd_(char *jobz, char *uplo, __CLPK_integer *n,
                       __CLPK_doublereal *a, __CLPK_integer *lda,
                       __CLPK_doublereal *w, __CLPK_doublereal *work,
                       __CLPK_integer *lwork, __CLPK_integer *iwork,
                       __CLPK_integer *liwork, __CLPK_integer *info);
#endif // MOAB_dsyevd
__CLPK_integer zgetri_(__CLPK_integer *n, __CLPK_doublecomplex *a,
                       __CLPK_integer *lda, __CLPK_integer *ipiv,
                       __CLPK_doublecomplex *work, __CLPK_integer *lwork,
                       __CLPK_integer *info);
__CLPK_integer zgetrf_(__CLPK_integer *m, __CLPK_integer *n,
                       __CLPK_doublecomplex *a, __CLPK_integer *lda,
                       __CLPK_integer *ipiv, __CLPK_integer *info);
__CLPK_integer zpotri_(char *uplo, __CLPK_integer *n, __CLPK_doublecomplex *a,
                       __CLPK_integer *lda, __CLPK_integer *info);
__CLPK_integer zpotrf_(char *uplo, __CLPK_integer *n, __CLPK_doublecomplex *a,
                       __CLPK_integer *lda, __CLPK_integer *info);

/// SUBROUTINE DPBTRS( UPLO, N, KD, NRHS, AB, LDAB, B, LDB, INFO )
__CLPK_integer dpbtrs_(char *uplo, __CLPK_integer *n, __CLPK_integer *kd,
                       __CLPK_integer *nrhs, __CLPK_doublereal *ab,
                       __CLPK_integer *ldab, __CLPK_doublereal *b,
                       __CLPK_integer *ldb, __CLPK_integer *info);

#endif

inline static __CLPK_integer lapack_sgetrf(__CLPK_integer m, __CLPK_integer n,
                                           float *a, __CLPK_integer lda,
                                           __CLPK_integer *ipiv) {
  __CLPK_integer info;
  sgetrf_(&m, &n, a, &lda, ipiv, &info);
  return info;
}

inline static __CLPK_integer lapack_dgetrf(__CLPK_integer m, __CLPK_integer n,
                                           __CLPK_doublereal *a,
                                           __CLPK_integer lda,
                                           __CLPK_integer *ipiv) {
  __CLPK_integer info;
  dgetrf_(&m, &n, a, &lda, ipiv, &info);
  return info;
}

inline static __CLPK_integer
lapack_dgetrs(char trans, __CLPK_integer n, __CLPK_integer nrhs,
              __CLPK_doublereal *a, __CLPK_integer lda, __CLPK_integer *ipiv,
              __CLPK_doublereal *b, __CLPK_integer ldb) {
  __CLPK_integer info;
  dgetrs_(&trans, &n, &nrhs, a, &lda, ipiv, b, &ldb, &info);
  return info;
}

inline static __CLPK_integer
lapack_dgesv(__CLPK_integer n, __CLPK_integer nrhs, __CLPK_doublereal *a,
             __CLPK_integer lda, __CLPK_integer *ipiv, __CLPK_doublereal *b,
             __CLPK_integer ldb) {
  __CLPK_integer info;
  dgesv_(&n, &nrhs, a, &lda, ipiv, b, &ldb, &info);
  return info;
}

inline static __CLPK_integer
lapack_dgetri(__CLPK_integer n, __CLPK_doublereal *a, __CLPK_integer lda,
              __CLPK_integer *ipiv, __CLPK_doublereal *work,
              __CLPK_integer lwork) {
  __CLPK_integer info;
  dgetri_(&n, a, &lda, ipiv, work, &lwork, &info);
  return info;
}

inline static __CLPK_integer lapack_dpotrf(char uplo, __CLPK_integer n,
                                           __CLPK_doublereal *a,
                                           __CLPK_integer lda) {
  __CLPK_integer info;
  dpotrf_(&uplo, &n, a, &lda, &info);
  return info;
}

inline static __CLPK_integer
lapack_dpotrs(char uplo, __CLPK_integer n, __CLPK_integer nrhs,
              __CLPK_doublereal *a, __CLPK_integer lda, __CLPK_doublereal *b,
              __CLPK_integer ldb) {
  __CLPK_integer info;
  dpotrs_(&uplo, &n, &nrhs, a, &lda, b, &ldb, &info);
  return info;
}

inline static __CLPK_integer
lapack_dposv(char uplo, __CLPK_integer n, __CLPK_integer nrhs,
             __CLPK_doublereal *a, __CLPK_integer lda, __CLPK_doublereal *b,
             __CLPK_integer ldb) {
  __CLPK_integer info;
  dposv_(&uplo, &n, &nrhs, a, &lda, b, &ldb, &info);
  return info;
}

inline static __CLPK_integer
lapack_dsysv(char uplo, __CLPK_integer n, __CLPK_integer nrhs,
             __CLPK_doublereal *a, __CLPK_integer lda, __CLPK_integer *ipiv,
             __CLPK_doublereal *b, __CLPK_integer ldb, __CLPK_doublereal *work,
             __CLPK_integer lwork) {
  __CLPK_integer info = 0;
  dsysv_(&uplo, &n, &nrhs, a, &lda, ipiv, b, &ldb, work, &lwork, &info);
  return info;
}

inline static __CLPK_integer lapack_dpotri(char uplo, __CLPK_integer n,
                                           __CLPK_doublereal *a,
                                           __CLPK_integer lda) {
  __CLPK_integer info;
  dpotri_(&uplo, &n, a, &lda, &info);
  return info;
}

inline static __CLPK_integer
lapack_dgesvd(char jobu, char jobvt, __CLPK_integer m, __CLPK_integer n,
              __CLPK_doublereal *a, __CLPK_integer lda, __CLPK_doublereal *s,
              __CLPK_doublereal *u, __CLPK_integer ldu, __CLPK_doublereal *vt,
              __CLPK_integer ldvt, __CLPK_doublereal *work,
              __CLPK_integer lwork) {
  __CLPK_integer info;
  dgesvd_(&jobu, &jobvt, &m, &n, a, &lda, s, u, &ldu, vt, &ldvt, work, &lwork,
          &info);
  return info;
}

inline static __CLPK_integer
lapack_sgesvd(char jobu, char jobvt, __CLPK_integer m, __CLPK_integer n,
              float *a, __CLPK_integer lda, float *s, float *u,
              __CLPK_integer ldu, float *vt, __CLPK_integer ldvt, float *work,
              __CLPK_integer lwork) {
  __CLPK_integer info;
  sgesvd_(&jobu, &jobvt, &m, &n, a, &lda, s, u, &ldu, vt, &ldvt, work, &lwork,
          &info);
  return info;
}

inline static __CLPK_integer
lapack_dsyev(char jobz, char uplo, __CLPK_integer n, __CLPK_doublereal *a,
             __CLPK_integer lda, __CLPK_doublereal *w, __CLPK_doublereal *work,
             __CLPK_integer lwork) {
  __CLPK_integer info;
  dsyev_(&jobz, &uplo, &n, a, &lda, w, work, &lwork, &info);
  return info;
}
inline static __CLPK_integer
lapack_zheev(char jobz, char uplo, __CLPK_integer n, __CLPK_doublecomplex *a,
             __CLPK_integer lda, __CLPK_doublereal *w,
             __CLPK_doublecomplex *work, __CLPK_integer lwork,
             __CLPK_doublereal *rwork) {
  __CLPK_integer info;
  zheev_(&jobz, &uplo, &n, a, &lda, w, work, &lwork, rwork, &info);
  return info;
}

inline static __CLPK_integer
lapack_zgeev(char jobvl, char jobvr, __CLPK_integer n, __CLPK_doublecomplex *a,
             __CLPK_integer lda, __CLPK_doublecomplex *w,
             __CLPK_doublecomplex *vl, __CLPK_integer ldvl,
             __CLPK_doublecomplex *vr, __CLPK_integer ldvr,
             __CLPK_doublecomplex *work, __CLPK_integer lwork,
             __CLPK_doublereal *rwork) {
  __CLPK_integer info;
  zgeev_(&jobvl, &jobvr, &n, a, &lda, w, vl, &ldvl, vr, &ldvr, work, &lwork,
         rwork, &info);
  return info;
}

inline static __CLPK_integer
lapack_dgels(char trans, __CLPK_integer m, __CLPK_integer n,
             __CLPK_integer nrhs, __CLPK_doublereal *a, __CLPK_integer lda,
             __CLPK_doublereal *b, __CLPK_integer ldb, __CLPK_doublereal *work,
             __CLPK_integer lwork) {
  __CLPK_integer info;
  dgels_(&trans, &m, &n, &nrhs, a, &lda, b, &ldb, work, &lwork, &info);
  return info;
}

inline static __CLPK_integer
lapack_dgesdd(char jobz, __CLPK_integer m, __CLPK_integer n,
              __CLPK_doublereal *a, __CLPK_integer lda, __CLPK_doublereal *s,
              __CLPK_doublereal *u, __CLPK_integer ldu, __CLPK_doublereal *vt,
              __CLPK_integer ldvt, __CLPK_doublereal *work,
              __CLPK_integer lwork, __CLPK_integer *iwork) {
  __CLPK_integer info;
  dgesdd_(&jobz, &m, &n, a, &lda, s, u, &ldu, vt, &ldvt, work, &lwork, iwork,
          &info);
  return info;
}

inline static __CLPK_integer
lapack_dsyevd(char jobz, char uplo, __CLPK_integer n, __CLPK_doublereal *a,
              __CLPK_integer lda, __CLPK_doublereal *w, __CLPK_doublereal *work,
              __CLPK_integer lwork, __CLPK_integer *iwork,
              __CLPK_integer liwork) {
  __CLPK_integer info;
  dsyevd_(&jobz, &uplo, &n, a, &lda, w, work, &lwork, iwork, &liwork, &info);
  return info;
}

inline static __CLPK_integer lapack_zgetrf(__CLPK_integer m, __CLPK_integer n,
                                           __CLPK_doublecomplex *a,
                                           __CLPK_integer lda,
                                           __CLPK_integer *ipiv) {
  __CLPK_integer info;
  zgetrf_(&m, &n, a, &lda, ipiv, &info);
  return info;
}

inline static __CLPK_integer
lapack_zgetri(__CLPK_integer n, __CLPK_doublecomplex *a, __CLPK_integer lda,
              __CLPK_integer *ipiv, __CLPK_doublecomplex *work,
              __CLPK_integer lwork) {
  __CLPK_integer info;
  zgetri_(&n, a, &lda, ipiv, work, &lwork, &info);
  return info;
}

inline static __CLPK_integer lapack_zpotri(char uplo, __CLPK_integer n,
                                           __CLPK_doublecomplex *a,
                                           __CLPK_integer lda) {
  __CLPK_integer info;
  zpotri_(&uplo, &n, a, &lda, &info);
  return info;
}
inline static __CLPK_integer lapack_zpotrf(char uplo, __CLPK_integer n,
                                           __CLPK_doublecomplex *a,
                                           __CLPK_integer lda) {
  __CLPK_integer info;
  zpotrf_(&uplo, &n, a, &lda, &info);
  return info;
}

inline static __CLPK_integer
lapack_dpbtrs(char uplo, __CLPK_integer n, __CLPK_integer kd,
              __CLPK_integer nrhs, __CLPK_doublereal *ab, __CLPK_integer ldab,
              __CLPK_doublereal *b, __CLPK_integer ldb) {
  __CLPK_integer info;
  dpbtrs_(&uplo, &n, &kd, &nrhs, ab, &ldab, b, &ldb, &info);
  return info;
}

#endif

#ifdef __cplusplus
}
#endif
