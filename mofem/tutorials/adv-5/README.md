# Running code

~~~~
../../tools/mofem_part \
-my_file test_seepage.cub -output_file mesh.h5m -my_nparts 1 -dim 2 -adj_dim 1
~~~~

~~~~
mpirun -np 1 ./seepage_2d -file_name mesh.h5m -ksp_monitor -order 2 -ts_max_steps 10 -ts_dt 0.1 -log_sl inform -field_eval_coords 0 0 
~~~~