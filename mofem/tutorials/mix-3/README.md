mofem_part -my_file surface.cub -output_file surface.h5m -my_nparts 2 -dim 2 -adj_dim 1

mpirun -np 2 ./mixed_nonlinear_poisson -file_name surface.h5m -order 2