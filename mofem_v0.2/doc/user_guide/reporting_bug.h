/*! \page guidelines_bug_reporting Guidelines for Bug Reports

The MoFEM maintenance e-mail, CMatGU <cmatgu@googlegroups.com>, is intended for users to

- report bugs,
- ask for clarification,
- ask for help in tracking down bugs, and
- request new features within MoFEM.

\section bug_reporting Guidelines for Bug Reports

The more information that you convey about a bug, the easier it will be for us to target the problem. We suggest providing the following information:

- Line command 
- Version and git commit id ( shown at beginning of analysis )
- Error trace or full log file 
- Mesh file ( if error is not big )
- Run ctest \code $MOFEM_BUILD_DIRECTORY/scripts/mofem_fast_check.sh \endcode

\section example_bug_repoty Example bug report

In this particular case a wrong name of mesh files was used to trigger error. Bug report for that case should look be as follows:

Line command:
\code
./elasticity -my_file aaa.cub -ksp_type fgmres -pc_type lu -pc_factor_mat_solver_package mumps -ksp_monitor
\endcode

Error Trace:
\code
[0]PETSC ERROR: --------------------- MoFEM Error Message---------------------------------------------------------------------------
[0]PETSC ERROR: MoFEM version 0.2.1
[0]PETSC ERROR: MoFEM git commit id b92aa7ece87bd423f082a11d120fd3b65c11c8f9
[0]PETSC ERROR: See http://userweb.eng.gla.ac.uk/lukasz.kaczmarczyk/MoFem/html/guidelines_bug_reporting.html for bug reporting.
[0]PETSC ERROR: See http://userweb.eng.gla.ac.uk/lukasz.kaczmarczyk/MoFem/html/faq_and_bugs.html for trouble shooting.
[0]PETSC ERROR: --------------------- Error Message --------------------------------------------------------------
[0]PETSC ERROR: Error code  7 at /mnt/home/MyBuild/mofem-cephas/mofem_v0.2/users_modules/elasticity/elasticity.cpp:112

[0]PETSC ERROR: See http://www.mcs.anl.gov/petsc/documentation/faq.html for trouble shooting.
[0]PETSC ERROR: Petsc Development GIT revision: v3.5.3-1524-gee900cc  GIT Date: 2015-01-31 17:44:15 -0600
[0]PETSC ERROR: /mofem_build/um_debug/elasticity/elasticity on a arch-linux2-c-debug named likask by root Fri Mar 20 14:59:28 2015
[0]PETSC ERROR: Configure options --download-blacs=1 --download-hypre=1 --download-metis=1 --download-moab --download-mumps=1 --download-netcdf=1 --download-parmetis=1 --download-ptscotch=1 --download-scalapack=1 --download-superlu_dist=1 --download-zoltan=1 --with-debugging=1 --with-hdf5-dir=/usr --with-hdf5=1 --with-mpi=1 --with-shared-libraries=1 1 = PETSC_ARCH=arch-linux2-c-debug
[0]PETSC ERROR: #1 main() line 112 in /mnt/home/MyBuild/mofem-cephas/mofem_v0.2/users_modules/elasticity/elasticity.cpp
[0]PETSC ERROR: PETSc Option Table entries:
[0]PETSC ERROR: -ksp_monitor
[0]PETSC ERROR: -ksp_type fgmres
[0]PETSC ERROR: -my_file aaaa.cub
[0]PETSC ERROR: -pc_factor_mat_solver_package mumps
[0]PETSC ERROR: -pc_type lu
[0]PETSC ERROR: ----------------End of Error Message -------send entire error message to petsc-maint@mcs.anl.gov----------
[0]PETSC ERROR: ----------MoFEM End of Error Message -------send entire error message to CMatGU <cmatgu@googlegroups.com> ----------
\endcode

*/
