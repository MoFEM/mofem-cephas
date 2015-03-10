/*! \page elastic_dynamics Elastic dynamic (Usage example)

\section dynamic_elastic_bar Prismatic 3d bar

\subsection dynamic_elastic_bar_input_file Input files

All files for this example can be found in directory
\em users_modules/nonlinear_elasticity/examples/prismatic_bar

  - mesh input file, .cub, .h5m. .vtk and more.
  - time data history, f.e. history_data.in


\subsubsection journal Journal file

Mesh file with boundary conditions, material properties can be created with Cubit journal file with prismatic bar (rod.jou),
\code
reset
set duplicate block elements on

brick x 1 y 1 z 10

create Displacement  on surface 1 dof 1 dof 2 dof 3 fix 0
create pressure  on surface 12 magnitude 2
#create force  on surface 2 force value 1 direction   y 

block 1 volume all 
block 1 name 'MAT_ELASTIC'
block 1 attribute count 2
block 1 attribute index 1 1 #young modulus
block 1 attribute index 2 0 #poisson ratio

block 2 volume 1
block 2 name "BODY_FORCES")
block 2 attribute count 4
block 2 attribute index 1 0.01	#material density
block 2 attribute index 2 0.	#constant acceleration in x-direction
block 2 attribute index 3 0.	#constant acceleration in y-direction
block 2 attribute index 4 0.	#constant acceleration in z-direction

volume all scheme Tetmesh
volume all size auto factor 9
mesh volume all
\endcode 

\subsubsection time_history Time history data

Time data history file (rod_history.in),
\code
0.0 	0.0
0.1 	1.0
0.1001	0.0
100000	0.0
\endcode

where first column represents time and second column represent force multiplier.

\subsection execution Executing code

If bar material model is NeoHookean (note that uisng key word KIRCHOFF a St Venant Kirchoff material is used)
\code
mpirun -np 4 ../../nonlinear_dynamics \
  -my_file rod.cub \
  -ksp_type fgmres -pc_type lu -pc_factor_mat_solver_package superlu_dist -ksp_atol 1e-10 -ksp_rtol 1e-10 \
  -snes_monitor -snes_type newtonls -snes_linesearch_type basic -snes_max_it 100 -snes_atol 1e-7 -snes_rtol 1e-7 \
  -ts_monitor -ts_type alpha -ts_dt 0.001 -ts_final_time 0.9 \
  -my_output_prt -10 -my_max_post_proc_ref_level 0 \
  -my_disp_order 2 \
  -default_material NEOHOOKEAN -my_time_data_file rod_history.in 2>&1 | tee log
\endcode

If bar is made form Hooke material subjected to small strains and displacements,
\code
mpirun -np 4 ../../nonlinear_dynamics \
  -my_file rod.cub \
  -ksp_type fgmres -pc_type lu -pc_factor_mat_solver_package superlu_dist -ksp_atol 1e-10 -ksp_rtol 1e-10 \
  -snes_monitor -snes_type newtonls -snes_linesearch_type basic -snes_max_it 100 -snes_atol 1e-7 -snes_rtol 1e-7 \
  -ts_monitor -ts_type alpha -ts_dt 0.001 -ts_final_time 0.9 \
  -my_output_prt -10 -my_max_post_proc_ref_level 0 \
  -my_disp_order 2 \
  -snes_lag_jacobian -2 \
  -default_material HOOKE -is_linear -my_time_data_file rod_history.in 2>&1 | tee log
\endcode



*/
