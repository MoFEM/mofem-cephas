# Running code

~~~~
../../tools/mofem_part \
-my_file test_seepage.cub -output_file mesh.h5m -my_nparts 1 -dim 2 -adj_dim 1
~~~~

~~~~
mpirun -np 1 ./seepage_2d_1 -file_name mesh.h5m -ksp_monitor -order 2 -ts_max_steps 10 -ts_dt 0.1 -young_modulus 0.0009 -poisson_ratio 0.2 -conductivity 0.4 -biot_constant 1 -storage 1 -log_sl inform -field_eval_coords 0 0 
~~~~

mpirun -np 1 ./seepage_3d -file_name mesh.h5m -ksp_monitor -order 2 -ts_max_steps 10 -ts_dt 0.1 -young_modulus 0.0009 -poisson_ratio 0.2 -conductivity 0.4 -biot_constant 1 -storage 1 -log_sl inform -field_eval_coords 0 0 0