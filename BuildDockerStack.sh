#!/bin/bash

mkdir -p empty

export VERSION=0.13.1

# Create environment
#docker build -t likask/mofem-spack-env -f Dockerfile-spack-env empty/
#docker tag likask/mofem-spack-env:latest likask/mofem-spack-env:$VERSION

Build core and users modules
docker build -t likask/mofem-spack-mofem -f Dockerfile-spack-mofem . 
docker tag likask/mofem-spack-mofem:latest likask/mofem-spack-mofem:$VERSION
docker tag likask/mofem-spack-mofem:latest likask/mofem-intermidiate:latest

Install softmech module
docker build -t likask/mofem-spack-softmech -f Dockerfile-spack-softmech empty/ 
docker tag likask/mofem-spack-softmech:latest likask/mofem-spack-softmech:$VERSION 
docker tag likask/mofem-spack-softmech:latest likask/mofem-intermidiate:latest

# Install other modules ...

Jupyter
docker build -t likask/mofem-spack-jupyter -f Dockerfile-spack-jupyter empty
docker tag likask/mofem-spack-jupyter:latest likask/mofem-spack-jupyter:$VERSION 
docker tag likask/mofem-spack-jupyter:latest likask/mofem-intermidiate:latest

# Hub
docker build -t likask/mofem-spack-jupyterhub -f Dockerfile-spack-jupyterhub jupyter
docker tag likask/mofem-spack-jupyterhub:latest likask/mofem-spack-jupyterhub:$VERSION
docker tag likask/mofem-spack-jupyterhub:latest likask/mofem-intermidiate:latest

# Volume
docker build -t likask/mofem-spack-build -f Dockerfile-spack-volume empty
docker tag likask/mofem-spack-build:latest likask/mofem-spack-build:$VERSION 



