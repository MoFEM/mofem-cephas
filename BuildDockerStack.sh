#!/bin/bash

mkdir -p empty

# Create environment
#docker build -t likask/mofem-spack-env \
#  -f Dockerfile-spack-env . 

# Build core and users modules
docker build -t likask/mofem-jupyter-um \
  -f Dockerfile-spack-build . 

# Install softmech module
docker build -t likask/mofem-jupyter-softmech \
  -f Dockerfile-spack-jupyter-softmech empty/ 

# Install other modules ...

# Create volume
docker build -t likask/mofem-spack-build \
  -f Dockerfile-spack-volume empty

# Make Jupyter
docker build -t likask/mofem-spack-jupyter \
  -f Dockerfile-spack-jupyter empty

  