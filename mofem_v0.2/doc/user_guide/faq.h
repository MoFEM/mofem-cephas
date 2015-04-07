/*! \page faqs Frequently Asked Questions

\section faq Frequently Asked Questions

\subsection partition_mesh How to partition mesh?

If problem is large, mesh can be partitioned for to save memory and improve efficiency. This can be done using native MoAB tool,
\code
mbpart -t -p PartKway  16 cylinder.cub cylinder_16parts.h5m
\endcode

\subsection running_multi_grid How to multi-grid solver via approximation orders?

Code is run using direct solver, i.e. \em Mumps on coarse level. Note that
loaded mesh is portioned and each processor only reads part of the mesh, i.e.
\em  -my_is_partitioned.

\code
mpirun -np 4 ./elasticity \
  -my_file dam_4parts.h5m -my_is_partitioned \
  -ksp_type gmres -ksp_max_it 1000  -ksp_atol 1e-13 -ksp_rtol 0  -ksp_monitor_lg_residualnorm  -ksp_final_residual -ksp_monitor -ksp_converged_reason 
  -my_order 1  -my_block_config block_congig.in   \
  -mofem_mg_verbose 1 -mofem_mg_coarse_order 1 -mofem_mg_levels 4  \
  -pc_type mg \
  -mg_coarse_ksp_type preonly -mg_coarse_pc_type lu -mg_coarse_pc_factor_mat_solver_package mumps \
  -pc_mg_smoothup 20 -pc_mg_smoothdown 20  -pc_mg_type multiplicative
\endcode

- Option \em -my_is_partitioned is set if mesh is partitioned using \em mbpart 
- Option \em -mofem_mg_coarse_order 1 set coarse level is for linear approximation, i.e. order 1. 
- Option \em -mofem_mg_levels 4 set number of multi-grid level. In that case maximal approx. order for some part of mesh is 4, thus 4 multi-grid levels.
- Option -pc_mg_smoothup 20 -pc_mg_smoothdown 20 set number of smoothing iterations, for more details look to PETSc manual.
- In line \code
-mg_coarse_ksp_type preonly -mg_coarse_pc_type lu -mg_coarse_pc_factor_mat_solver_package mumps
\endcode a direct solver for coarse mesh is set.


\subsection update_on_memory_stick Ho to update MoFEM on Live USB Stick?

MoFEM update on Live USB stick:
\code
$ mofem_update.sh
$ mofem_build.sh
\endcode

Following command run test verifying updated code:
\code
$ mofem_fast_check.sh
\endcode

If you run MoFEM update at University of Glasgow behind proxy server, set proxy
servers as follows:
\code
$ export http_proxy=http://wwwcache.gla.ac.uk:8080
$ export https_proxy=http://wwwcache.gla.ac.uk:8080
\endcode

\subsection ctest How to run ctest?

You can run tests and report results to MoFEM CDash web page. Form mofem user
modules build directory executing  run script
\code./bin/mofem_fast_check.sh\endcode
Results of test can be seen on <http://cdash.eng.gla.ac.uk/cdash/>. 

Note that test tests for MoFEM library and \em User \em Modules are run independently
and can be seen as a two different projects.

If you run test behind proxy server you can have to set \em http_proxy and \em
http_proxy environmental variables. For example if you run mofem at Glasgow
University, please do:
\code
$ export http_proxy=http://wwwcache.gla.ac.uk:8080
$ export https_proxy=http://wwwcache.gla.ac.uk:8080
\endcode

You can as well run ctest directly by simply executing command line:
\code 
ctest -V -D Experimental
\endcode
where option -V sets verbose version and all test output is printed on screen
and -D Experimental tels ctest to submit results to Experimental build on CDash
MoFEM server.


*/


