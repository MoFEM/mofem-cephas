#!/bin/bash

mkdir -p empty

export ACCOUNT=agshvarts
export VERSION=mos-demo

# Create environment
docker build -t $ACCOUNT/mofem-spack-env -f Dockerfile-spack-env .
docker tag $ACCOUNT/mofem-spack-env:latest $ACCOUNT/mofem-spack-env:$VERSION

#Build core and users modules
docker build -t $ACCOUNT/mofem-spack-mofem -f Dockerfile-spack-mofem .
docker tag $ACCOUNT/mofem-spack-mofem:latest $ACCOUNT/mofem-spack-mofem:$VERSION
# docker tag $ACCOUNT/mofem-spack-mofem:latest $ACCOUNT/mofem-intermidiate:latest

#Install softmech module
#docker build -t likask/mofem-spack-softmech -f Dockerfile-spack-softmech empty/ 
#docker tag likask/mofem-spack-softmech:latest likask/mofem-spack-softmech:$VERSION 
#docker tag likask/mofem-spack-softmech:latest likask/mofem-intermidiate:latest

# Install other modules ...

#Jupyter
# docker build -t likask/mofem-spack-jupyter -f Dockerfile-spack-jupyter empty
# docker tag likask/mofem-spack-jupyter:latest likask/mofem-spack-jupyter:$VERSION 
# docker tag likask/mofem-spack-jupyter:latest likask/mofem-intermidiate:latest

# Labs indentation
#dodocker build -t likask/mofem-spack-jupyter-labs-indentation:VERSION-new -f Dockerfile-spack-jupyter-labs-indentation .

# Hub
# docker build -t likask/mofem-spack-jupyterhub -f Dockerfile-spack-jupyterhub jupyter
# docker tag likask/mofem-spack-jupyterhub:latest likask/mofem-spack-jupyterhub:$VERSION
# docker tag likask/mofem-spack-jupyterhub:latest likask/mofem-intermidiate:latest

# Labs indentation
#docker build -t likask/mofem-spack-jupyterhub-labs-indentation:VERSION-new -f Dockerfile-spack-jupyterhub-labs-indentation .

# Volume
#docker build -t likask/mofem-spack-build -f Dockerfile-spack-volume empty
#docker tag likask/mofem-spack-build:latest likask/mofem-spack-build:$VERSION 



