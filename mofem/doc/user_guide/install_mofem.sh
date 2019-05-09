#!/bin/bash
# This script automises installation of MoFEM user modules and fracture module
# Full description of the installation process using spack can be found at:
# http://mofem.eng.gla.ac.uk/mofem/html/install_spack.html
# The script should work for both Ubuntu and macOS platforms
# The installation will take two to three hours
#
# Usage:
#       1. Copy install_mofem.sh to the directory where MoFEM will be instaslled
#       2. Run install_mofem.sh from the command line
#

##############################
# INITIALISATION
##############################

# Setup installation directory
pwd
export MOFEM_INSTALL_DIR="$PWD/mofem_install"
mkdir -p $MOFEM_INSTALL_DIR

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
    NumberOfProcs=$(($(cat /proc/cpuinfo | grep -i "^processor" | wc -l)/2))
elif [ ${machine} = "Mac" ]
then
    NumberOfProcs=$(($(sysctl -n hw.ncpu)/2))
fi

if
    [ "$NumberOfProcs" -lt 1 ]
then 
    NumberOfProcs=1
fi

echo "There is $NumberOfProcs processor(s)"


##############################
### SPACK
##############################

echo -e "****************************\nInstalling SPACK...\n****************************"

# Locate home directory
cd ~

# Retrieve Spack for MoFEM
git clone --single-branch -b mofem https://github.com/likask/spack.git

# Initialise Spack environment variables:
. $HOME/spack/share/spack/setup-env.sh

# Add command to .bashrc
echo ". $HOME/spack/share/spack/setup-env.sh" >> ~/.bashrc

# Install packages required by Spack
spack bootstrap

echo -e "Done.\n"


##############################
### PREREQUISITES
##############################

echo -e "****************************\nInstalling PREREQUISITES...\n****************************"
echo -e "User password can be asked at some point. Please wait..."

# Install appropriate prerequisite packages
if [ ${machine} = "Linux" ]
then
    apt-get update \
    && apt-get install -y --no-install-recommends \
    autoconf \
    build-essential \
    ca-certificates \
    coreutils \
    curl \
    environment-modules \
    git \
    python \
    unzip \
    vim \ 
    gfortran

elif [ ${machine} = "Mac" ]
then
    ## Assuming Xcode already installed
    # xcode-select --install
    # sudo xcodebuild -license accept
    # /usr/bin/ruby -e "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install)"
    # brew install curl git gcc

    spack install -j $NumberOfProcs gcc
    spack load gcc
    spack compiler find
    spack install -j $NumberOfProcs curl
    spack load curl
fi

echo -e "Done.\n"
echo -e "No user password will be asked from now on.\nYou can come back and check in few hours."

##############################
### MoFEM USER MODULES
##############################

echo "****************************\nInstalling MoFEM USER MODULES...\n****************************"

# Locate MoFEM installation directory
cd $MOFEM_INSTALL_DIR

# Install user modules
spack install -j $NumberOfProcs mofem-users-modules

# Activate user modules
spack view --verbose symlink -i um_view mofem-cephas
spack activate -v um_view mofem-users-modules

# Test elasticity module
cd um_view/elasticity
./elasticity \
-my_file LShape.h5m \
-my_order 2 \
-ksp_type gmres \
-pc_type lu -pc_factor_mat_solver_package mumps \
-ksp_monitor 2>&1 | tee log

export DIR="$PWD/um_view"
echo -e "Done.\n"
echo "MoFEM user modules were sucessfully installed and tested in $DIR. "


##############################
### MoFEM FRACTURE MODULES
##############################

echo "****************************\nInstalling MoFEM FRACTURE MODULE...\n****************************"

# Locate MoFEM installation directory
cd $MOFEM_INSTALL_DIR

# Install fracture module
spack install  -j $NumberOfProcs mofem-fracture-module

# Activate fracture module
spack view --verbose symlink -i um_view_fracture mofem-cephas
spack activate -v um_view_fracture mofem-fracture-module

# Test fracture module
cd um_view_fracture/mofem_um_fracture_mechanics 
./crack_propagation \
-my_file examples/analytical_bc/out_10.h5m \
-my_max_post_proc_ref_level 0 \
-is_conservative_force false \
-my_order 2 \
-my_ref 0\
-my_geom_order 1 \
-my_ref_order 2 \
-material HOOKE \
-mofem_mg_verbose 1 \
-mofem_mg_coarse_order 1 \
-mofem_mg_levels 2 \
-analytical_alpha 0  2>&1 | tee log

export DIR="$PWD/um_view_fracture"
echo -e "Done.\n"
echo "MoFEM fracture module was succesfully installed and tested in $DIR... "