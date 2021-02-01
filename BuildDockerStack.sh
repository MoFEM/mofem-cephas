#!/bin/bash

mkdir -p empty

export VERSION=0.10.2

# Create environment
docker build -t likask/mofem-spack-env -f Dockerfile-spack-env empty/
docker tag likask/mofem-spack-env:latest likask/mofem-spack-env:$VERSION

# Build core and users modules
docker build -t likask/mofem-spack-mofem -f Dockerfile-spack-mofem . 
docker tag likask/mofem-spack-mofem:latest ikask/mofem-spack-mofem:$VERSION

# Install softmech module
docker build -t likask/mofem-jupyter-softmech -f Dockerfile-spack-softmech empty/ 
docker tag likask/mofem-softmech:latest likask/mofem-softmech:$VERSION 
# Install other modules ...


# Tag last module
docker tag likask/mofem-jupyter-softmech:latest likask/mofem-intermidiate:latest

# Volume
docker build -t likask/mofem-spack-build -f Dockerfile-spack-volume empty
docker tag likask/mofem-spack-build:latest likask/mofem-spack-build:$VERSION 

# Jupyter
docker build -t likask/mofem-spack-jupyter -f Dockerfile-spack-jupyter empty
docker tag likask/likask/mofem-spack-jupyter:latest likask/likask/mofem-spack-jupyter:$VERSION 

  