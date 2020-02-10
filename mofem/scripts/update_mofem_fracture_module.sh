#!/bin/bash
# This script automises the update to a specific version of fracture module.
# Newest version of fracture module is installed as default.
# This script should be used only if MoFEM was installed using spack in advance.
# If MoFEM has not been installed, install it using install_mofem_user.sh or
# install_mofem_developer.sh before running this script
# The update should take few minutes.
#
# Usage:
#       1. Copy update_mofem_fracture_module.sh to the directory where new
#       version of fracture module will be installed.
# 
#       2. Run './update_mofem_fracture_module.sh fracture_version' from the 
#       command line where fracture_version is the version of choice.
#       If fracture_version is not specified, the newest version will be 
#       installed. 
#       For example, to install MoFEM Fracture Module version 0.9.60, 
#       run this from the command line: ./update_mofem_fracture_module.sh 0.9.60
#

echo "Start time: $(date +"%T")"

##############################
# INITIALISATION
##############################

# Set the version of fracture module to be installed
MOFEM_FM_NEWEST_VER=`spack info mofem-fracture-module | grep tag | awk '{print $1}' | head -n 1`
if [ $# -eq 0 ]
then
    echo "No version specified. Newest version will be installed."
    MOFEM_FM_VER=$MOFEM_FM_NEWEST_VER
else
    if spack info mofem-fracture-module | grep -Fq  "$1  "  
    then
      MOFEM_FM_VER=$1
    else
      echo -e "\nError: No matching version of MoFEM Fracture Module!\n"
      exit 1
    fi
fi
echo -e "\nInstalling MoFEM Fracture Module version $MOFEM_FM_VER ...\n"

# Setup installation directory
export MOFEM_FM_DIR="$PWD/fracture_module_v$MOFEM_FM_VER"
mkdir -p $MOFEM_FM_DIR
echo "Directory created: $MOFEM_FM_DIR"
  
# Check operating system
unameOut="$(uname -s)"
case "${unameOut}" in
    Linux*)     machine=Linux;;
    Darwin*)    machine=Mac;;
    CYGWIN*)    machine=Cygwin;;
    MINGW*)     machine=MinGw;;
    *)          machine="UNKNOWN:${unameOut}"
esac
  
# Get numbers of processor
if [ ${machine} = "Linux" ]
then
    NumberOfProcs=$(nproc)
elif [ ${machine} = "Mac" ]
then
    NumberOfProcs=$(($(sysctl -n hw.ncpu)/2))
fi
# max(NumberOfProcs, 1)
NumberOfProcs=$(( NumberOfProcs > 1 ? NumberOfProcs : 1 ))
  
echo "The number of processors is $NumberOfProcs"

# Locate home directory
cd $MOFEM_FM_DIR


# Download mirror
SPACK_MIRROR_DIR=$MOFEM_FM_DIR/mofem_mirror

if [ ! -d "$SPACK_MIRROR_DIR" ]; then
  if [ ! -f "$PWD/mirror.tgz" ]; then
    echo "Downloading spack mofem mirror ..."
    mkdir -p $SPACK_MIRROR_DIR && \
    curl -s -L https://bitbucket.org/likask/mofem-cephas/downloads/mirror_v0.9.1.tar.gz \
    | tar xzC $SPACK_MIRROR_DIR --strip 1
    echo -e "Done.\n"
  else 
    mkdir -p $SPACK_MIRROR_DIR && \
    tar xzf $PWD/mirror.tgz -C $SPACK_MIRROR_DIR  --strip 1
  fi
fi

# Add mirror
spack mirror remove mofem_mirror 2> /dev/null
spack mirror add mofem_mirror $SPACK_MIRROR_DIR

# Uninstall same version if already installed
yes | spack uninstall --all mofem-fracture-module@$MOFEM_FM_VER 2> /dev/null

# Install MoFEM Fracture Module
spack install mofem-fracture-module@$MOFEM_FM_VER build_type=Release

MOFEM_FM_SPACK_ID=`spack find -lvd mofem-fracture-module@$MOFEM_FM_VER | grep mofem-fracture-module | awk '{print $1}'| tail -n 1`
echo "mofem-fracture-module ID: $MOFEM_FM_SPACK_ID"

MOFEM_CEPHAS_VER=`spack spec mofem-fracture-module@$MOFEM_FM_VER | grep mofem-cephas | sed 's/\%.*//' | sed 's/^[^@]*@//g'`
echo "mofem-cephas version: $MOFEM_CEPHAS_VER"

# Activate fracture module
spack view --verbose symlink -i um_view mofem-cephas@$MOFEM_CEPHAS_VER
spack activate -v um_view mofem-fracture-module@$MOFEM_FM_VER

# Test fracture crack propagation
cd $MOFEM_FM_DIR/um_view/fracture_mechanics
echo "Current directory: $PWD"
./crack_propagation \
-my_file examples/analytical_bc/out_10.h5m \
-my_order 2 \
-my_ref 0 2>&1 | tee log

# Check the output message and finalise the installation
if tail -n 1 log | grep -q "Done rank = 0"
then
  echo -e "\nUpdate SUCCESSFUL!\n"
  
  # Export view and make view visible from any directory
  export PATH=$PWD/um_view/bin:$PATH 
  # Add PATH to .bashrc on Ubuntu or .bash_profile on Mac
  if [ ${machine} = "Linux" ]
  then
    echo "export PATH=$PWD/um_view/bin:\$PATH" >> ~/.bashrc
  elif [ ${machine} = "Mac" ]
  then
    echo "export PATH=$PWD/um_view/bin:\$PATH" >> ~/.bash_profile
  fi

  echo "Please check PATH in .bashrc (Ubuntu) or .bash_profile (macOS) and remove the old ones."
else
   echo -e "\nUpdate FAILED!\n"
fi

echo "End time: $(date +"%T")"