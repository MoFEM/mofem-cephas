/** \file DeprecatedPetsc.hpp
 * \brief Deprecated PETSc functions
 */



#ifndef __DEPRECATED_PETSC_HPP__
#define __DEPRECATED_PETSC_HPP__

namespace MoFEM {

#if PETSC_VERSION_GE(3, 7, 0)

/**
\deprecated Function is deprecated use function
http://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/Sys/PetscOptionsGetInt.html>
*/
inline PetscErrorCode PetscOptionsGetInt(const char pre[], const char name[],
                                         PetscInt *ivalue, PetscBool *set) {
  PetscErrorCode ierr;
  MoFEMFunctionBeginHot;
  ierr = ::PetscOptionsGetInt(PETSC_NULL, pre, name, ivalue, set);
  CHKERRQ(ierr);
  MoFEMFunctionReturnHot(0);
}

/**
\deprecated Function is deprecated use function
http://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/Sys/PetscOptionsGetReal.html>
*/
inline PetscErrorCode PetscOptionsGetReal(const char pre[], const char name[],
                                          PetscReal *dval, PetscBool *set) {
  PetscErrorCode ierr;
  MoFEMFunctionBeginHot;
  ierr = ::PetscOptionsGetReal(PETSC_NULL, pre, name, dval, set);
  CHKERRQ(ierr);
  MoFEMFunctionReturnHot(0);
}

/**
\deprecated Function is deprecated use function
http://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/Sys/PetscOptionsGetScalar.html>
*/
inline PetscErrorCode PetscOptionsGetScalar(const char pre[], const char name[],
                                            PetscScalar *dval, PetscBool *set) {
  PetscErrorCode ierr;
  MoFEMFunctionBeginHot;
  ierr = ::PetscOptionsGetScalar(PETSC_NULL, pre, name, dval, set);
  CHKERRQ(ierr);
  MoFEMFunctionReturnHot(0);
}

/**
\deprecated Function is deprecated use function
http://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/Sys/PetscOptionsGetString.html>
*/
inline PetscErrorCode PetscOptionsGetString(const char pre[], const char name[],
                                            char str[], size_t size,
                                            PetscBool *set) {
  MoFEMFunctionBeginHot;
  ierr = ::PetscOptionsGetString(PETSC_NULL, pre, name, str, size, set);
  CHKERRQ(ierr);
  MoFEMFunctionReturnHot(0);
}

/**
\deprecated Function is deprecated
usectiontp://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/Sys/PetscOptionsGetBool.html>
*/
inline PetscErrorCode PetscOptionsGetBool(const char pre[], const char name[],
                                          PetscBool *bval, PetscBool *set) {
  PetscErrorCode ierr;
  MoFEMFunctionBeginHot;
  ierr = ::PetscOptionsGetBool(PETSC_NULL, pre, name, bval, set);
  CHKERRQ(ierr);
  MoFEMFunctionReturnHot(0);
}

/**
\deprecated Function is deprecated
usectiontp://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/Sys/PetscOptionsGetRealArray.html>
*/
inline PetscErrorCode PetscOptionsGetRealArray(const char pre[],
                                               const char name[],
                                               PetscReal dval[], PetscInt *nmax,
                                               PetscBool *set) {
  PetscErrorCode ierr;
  MoFEMFunctionBeginHot;
  ierr = ::PetscOptionsGetRealArray(PETSC_NULL, pre, name, dval, nmax, set);
  CHKERRQ(ierr);
  MoFEMFunctionReturnHot(0);
}

/**
\deprecated Function is deprecated
usectiontp://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/Sys/PetscOptionsGetEList.html>
*/
inline PetscErrorCode PetscOptionsGetEList(const char pre[], const char name[],
                                           const char *const *list,
                                           PetscInt next, PetscInt *value,
                                           PetscBool *set) {
  PetscErrorCode ierr;
  MoFEMFunctionBeginHot;
  ierr = ::PetscOptionsGetEList(PETSC_NULL, pre, name, list, next, value, set);
  CHKERRQ(ierr);
  MoFEMFunctionReturnHot(0);
}

/**
\deprecated Function is deprecated
usectiontp://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/Sys/PetscOptionsGetIntArray.html>
*/
inline PetscErrorCode PetscOptionsGetIntArray(const char pre[],
                                              const char name[],
                                              PetscInt dvalue[], PetscInt *nmax,
                                              PetscBool *set) {
  PetscErrorCode ierr;
  MoFEMFunctionBeginHot;
  ierr = ::PetscOptionsGetIntArray(PETSC_NULL, pre, name, dvalue, nmax, set);
  CHKERRQ(ierr);
  MoFEMFunctionReturnHot(0);
}

/**
 *\deprecated Function is deprecated
 *usectiontp://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/Sys/PetscOptionsGetScalarArray.html>
 */
inline PetscErrorCode PetscOptionsGetScalarArray(const char pre[],
                                                 const char name[],
                                                 PetscScalar dvalue[],
                                                 PetscInt *nmax,
                                                 PetscBool *set) {
  PetscErrorCode ierr;
  MoFEMFunctionBeginHot;
  ierr = ::PetscOptionsGetScalarArray(PETSC_NULL, pre, name, dvalue, nmax, set);
  CHKERRQ(ierr);
  MoFEMFunctionReturnHot(0);
}

#else

inline PetscErrorCode PetscOptionsGetInt(PetscOptions *, const char pre[],
                                         const char name[], PetscInt *ivalue,
                                         PetscBool *set) {
  PetscErrorCode ierr;
  MoFEMFunctionBeginHot;
  ierr = ::PetscOptionsGetInt(pre, name, ivalue, set);
  CHKERRQ(ierr);
  MoFEMFunctionReturnHot(0);
}

inline PetscErrorCode PetscOptionsGetReal(PetscOptions *, const char pre[],
                                          const char name[], PetscReal *dval,
                                          PetscBool *set) {
  PetscErrorCode ierr;
  MoFEMFunctionBeginHot;
  ierr = ::PetscOptionsGetReal(pre, name, dval, set);
  CHKERRQ(ierr);
  MoFEMFunctionReturnHot(0);
}

inline PetscErrorCode PetscOptionsGetScalar(PetscOptions *, const char pre[],
                                            const char name[],
                                            PetscScalar *dval, PetscBool *set) {
  PetscErrorCode ierr;
  MoFEMFunctionBeginHot;
  ierr = ::PetscOptionsGetScalar(pre, name, dval, set);
  CHKERRQ(ierr);
  MoFEMFunctionReturnHot(0);
}

inline PetscErrorCode PetscOptionsGetString(PetscOptions *, const char pre[],
                                            const char name[], char str[],
                                            size_t size, PetscBool *set) {
  PetscErrorCode ierr;
  MoFEMFunctionBeginHot;
  ierr = ::PetscOptionsGetString(pre, name, str, size, set);
  CHKERRQ(ierr);
  MoFEMFunctionReturnHot(0);
}

inline PetscErrorCode PetscOptionsGetBool(PetscOptions *, const char pre[],
                                          const char name[], PetscBool *bval,
                                          PetscBool *set) {
  PetscErrorCode ierr;
  MoFEMFunctionBeginHot;
  ierr = ::PetscOptionsGetBool(pre, name, bval, set);
  CHKERRQ(ierr);
  MoFEMFunctionReturnHot(0);
}

inline PetscErrorCode PetscOptionsGetRealArray(PetscOptions *, const char pre[],
                                               const char name[],
                                               PetscReal dval[], PetscInt *nmax,
                                               PetscBool *set) {
  PetscErrorCode ierr;
  MoFEMFunctionBeginHot;
  ierr = ::PetscOptionsGetRealArray(pre, name, dval, nmax, set);
  CHKERRQ(ierr);
  MoFEMFunctionReturnHot(0);
}

inline PetscErrorCode PetscOptionsGetEList(PetscOptions *, const char pre[],
                                           const char name[],
                                           const char *const *list,
                                           PetscInt next, PetscInt *value,
                                           PetscBool *set) {
  PetscErrorCode ierr;
  MoFEMFunctionBeginHot;
  ierr = ::PetscOptionsGetEList(pre, name, list, next, value, set);
  CHKERRQ(ierr);
  MoFEMFunctionReturnHot(0);
}

inline PetscErrorCode PetscOptionsGetIntArray(PetscOptions options,
                                              const char pre[],
                                              const char name[],
                                              PetscInt dvalue[], PetscInt *nmax,
                                              PetscBool *set) {
  PetscErrorCode ierr;
  MoFEMFunctionBeginHot;
  ierr = ::PetscOptionsGetIntArray(pre, name, dvalue, nmax, set);
  CHKERRQ(ierr);
  MoFEMFunctionReturnHot(0);
}

inline PetscErrorCode
PetscOptionsGetScalarArray(PetscOptions options, const char pre[],
                           const char name[], PetscScalar dvalue[],
                           PetscInt *nmax, PetscBool *set) {
  PetscErrorCode ierr;
  MoFEMFunctionBeginHot;
  ierr = ::PetscOptionsGetScalarArray(pre, name, dvalue, nmax, set);
  CHKERRQ(ierr);
  MoFEMFunctionReturnHot(0);
}

#endif

} // namespace MoFEM

#endif //__DEPRECATED_PETSC_HPP__
