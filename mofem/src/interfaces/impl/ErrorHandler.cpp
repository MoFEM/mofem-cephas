
#if PETSC_VERSION_GE(3, 8, 0)
#include <petsc/private/petscimpl.h>
#endif

static PetscErrorCode mofem_error_handler(MPI_Comm comm, int line,
                                          const char *fun, const char *file,
                                          PetscErrorCode n, PetscErrorType p,
                                          const char *mess, void *ctx) {

  static int cnt = 1;
  MoFEMFunctionBeginHot;
  PetscVFPrintf = PetscVFPrintfDefault;
  MoFEM::LogManager::dummy_mofem_fd = PETSC_STDERR;

  int rank = 0;
  if (comm != PETSC_COMM_SELF)
    MPI_Comm_rank(comm, &rank);

  if (!rank) {

    if (p == PETSC_ERROR_INITIAL) {

      char petsc_version[255];
      PetscGetVersion(petsc_version, 255);
      MOFEM_LOG_CHANNEL("SELF");

      MOFEM_LOG_C("SELF", MoFEM::Sev::error, "%s",
                  "--------------------- MoFEM Error Message "
                  "---------------------------------------------------");

      MOFEM_LOG_C("SELF", MoFEM::Sev::error, "MoFEM version %d.%d.%d (%s %s)",
                  MoFEM_VERSION_MAJOR, MoFEM_VERSION_MINOR, MoFEM_VERSION_BUILD,
                  MOAB_VERSION_STRING, petsc_version);

      MOFEM_LOG_C("SELF", MoFEM::Sev::error, "MoFEM git commit id %s",
                  GIT_SHA1_NAME);

      MOFEM_LOG_C("SELF", MoFEM::Sev::error, "%s",
                  "See http://mofem.eng.gla.ac.uk/mofem/html/ "
                  "guidelines_bug_reporting.html for bug reporting.");

      MOFEM_LOG_C(
          "SELF", MoFEM::Sev::error, "%s",
          "Write to https://groups.google.com/forum/#!forum/mofem-group to "
          "seek help.");
    }

    if (n >= MOFEM_DATA_INCONSISTENCY) {

      if (p == PETSC_ERROR_INITIAL) {
        if (mess)
          MOFEM_LOG_C("SELF", MoFEM::Sev::error, "%s", mess);
      }
      MOFEM_LOG_C("SELF", MoFEM::Sev::error, "#%d %s() line %d in %s", cnt++,
                  fun, line, file);

    } else {
      PetscTraceBackErrorHandler(PETSC_COMM_SELF, line, fun, file, n, p, mess,
                                 ctx);
    }

    PetscBool ismain, isunknown;
    PetscStrncmp(fun, "main", 4, &ismain);
    PetscStrncmp(fun, "unknown", 7, &isunknown);
    if (ismain || isunknown) {

      if (n >= MOFEM_DATA_INCONSISTENCY) {
#if PETSC_VERSION_GE(3, 7, 0)
        PetscOptionsView(NULL, PETSC_VIEWER_STDERR_SELF);
#else
        PetscOptionsView(PETSC_VIEWER_STDERR_SELF);
#endif
      }

      // error_printf_highlight();
      MOFEM_LOG_C("SELF", MoFEM::Sev::error, "%s",
                  "-- MoFEM End of Error Message -- send entire error "
                  "message to mofem-group@googlegroups.com --");
      // error_printf_normal();
    }

  } else {

    /* do not print error messages since process 0 will print them, sleep before
     * aborting so will not accidentally kill process 0*/
    PetscSleep(10.0);
    abort();
  }

  MoFEMFunctionReturnHot(n);
}