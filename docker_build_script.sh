#!/bin/bash

# Create volume
docker volume create --name mofem_build

# Pull docker conatiner
docker pull likask/ubuntu_mofem:latest

# Build mofem from scrach
docker run \
  -it --rm=true \
  -v mofem_build:/mofem_build \
  -v $PWD/mofem:/mofem \
  -v $HOME:$HOME \
  likask/ubuntu_mofem:v0.9 \
  /bin/bash /mofem/docker_build_mofem_script.sh
