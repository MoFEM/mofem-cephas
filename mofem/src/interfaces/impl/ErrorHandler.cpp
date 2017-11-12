/* MoFEM is free software: you can redistribute it and/or modify it under
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

static PetscErrorCode mofem_error_printf(const char format[], ...) {
  va_list Argp;
  static PetscBool PetscErrorPrintfCalled = PETSC_FALSE;

  /*
      This function does not call PetscFunctionBegin and PetscFunctionReturn()
    because it may be called by PetscStackView().

      This function does not do error checking because it is called by the error
    handlers.
  */

  if (!PetscErrorPrintfCalled) {
    PetscErrorPrintfCalled = PETSC_TRUE;

    /*
        On the SGI machines and Cray T3E, if errors are generated
      "simultaneously" by different processors, the messages are printed all
      jumbled up; to try to prevent this we have each processor wait based on
      their rank
    */
#if defined(PETSC_CAN_SLEEP_AFTER_ERROR)
    {
      PetscMPIInt rank;
      if (PetscGlobalRank > 8)
        rank = 8;
      else
        rank = PetscGlobalRank;
      PetscSleep((PetscReal)rank);
    }
#endif
  }

  PetscFPrintf(PETSC_COMM_SELF, PETSC_STDERR,
               "[%d]MoFEM ERROR: ", PetscGlobalRank);
  va_start(Argp, format);
  (*PetscVFPrintf)(PETSC_STDERR, format, Argp);
  va_end(Argp);
  return 0;
}

static void error_printf_highlight(void) {
#if defined(PETSC_HAVE_UNISTD_H) && defined(PETSC_USE_ISATTY)
  if (isatty(fileno(PETSC_STDERR)))
    fprintf(PETSC_STDERR, "\033[1;32m");
#endif
}

static void error_printf_normal(void) {
#if defined(PETSC_HAVE_UNISTD_H) && defined(PETSC_USE_ISATTY)
  if (isatty(fileno(PETSC_STDERR)))
    fprintf(PETSC_STDERR, "\033[0;39m\033[0;49m");
#endif
}

static PetscErrorCode mofem_error_handler(MPI_Comm comm, int line,
                                          const char *fun, const char *file,
                                          PetscErrorCode n, PetscErrorType p,
                                          const char *mess, void *ctx) {

  static int cnt = 1;
  MoFEMFunctionBeginHot;

  int rank = 0;
  if (comm != PETSC_COMM_SELF)
    MPI_Comm_rank(comm, &rank);

  if (!rank) {

    if (p == PETSC_ERROR_INITIAL) {
      error_printf_highlight();
      mofem_error_printf("--------------------- MoFEM Error Message ---------------------------------------------------\n");
      mofem_error_printf("MoFEM version %d.%d.%d\n", MoFEM_VERSION_MAJOR,
                         MoFEM_VERSION_MINOR, MoFEM_VERSION_BUILD);
      mofem_error_printf("MoFEM git commit id %s\n", GIT_SHA1_NAME);
      mofem_error_printf("See "
                         "http://mofem.eng.gla.ac.uk/mofem/html/"
                         "guidelines_bug_reporting.html for bug reporting.\n");
      mofem_error_printf("See "
                         "http://mofem.eng.gla.ac.uk/mofem/html/"
                         "faq_and_bugs.html for trouble shooting.\n");
      error_printf_normal();
    }

    if (n >= MOFEM_DATA_INCONSISTENCY) {
      if (p == PETSC_ERROR_INITIAL) {
        if (mess)
          mofem_error_printf("%s\n", mess);
      }
      mofem_error_printf("#%d %s() line %d in %s\n", cnt++, fun, line, file);
    } else {
      PetscTraceBackErrorHandler(PETSC_COMM_SELF, line, fun, file, n, p, mess,
                                 ctx);
    }

    PetscBool ismain, isunknown;
    PetscStrncmp(fun, "main", 4, &ismain);
    PetscStrncmp(fun, "unknown", 7, &isunknown);
    if (ismain || isunknown) {

      std::stringstream strs_version;
      strs_version << "MoFEM_version_" << MoFEM_VERSION_MAJOR << "."
                   << MoFEM_VERSION_MINOR << "." << MoFEM_VERSION_BUILD;

      error_printf_highlight();
      mofem_error_printf("-- MoFEM End of Error Message -- send entire error message to mofem-group@googlegroups.com --\n");
      error_printf_normal();
    }

  } else {

    /* do not print error messages since process 0 will print them, sleep before
     * aborting so will not accidentally kill process 0*/
    PetscSleep(10.0);
    abort();
  }


  MoFEMFunctionReturnHot(n);
}