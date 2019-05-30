#!/bin/bash
# This script automises installation of MoFEM user version with basic modules
# and fracture module
# Full description of the installation process using spack can be found at:
# http://mofem.eng.gla.ac.uk/mofem/html/install_spack.html
# The script should work for both Ubuntu and macOS platforms
# The last known working platforms are Ubuntu (18.04) and macOS (10.13, 10.14)
# The installation may take two to three hours
#
# Usage:
#       1. Copy install_mofem_user.sh to the directory where MoFEM will be installed
#       2. Run install_mofem_user.sh from the command line
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
  
echo "The number of processors is $NumberOfProcs"
  
  
##############################
### PREREQUISITES
##############################
  
echo -e "\n****************************\nInstalling PREREQUISITES...\n****************************\n"
echo "User password can be asked at some point. Please wait..."
  
# Install appropriate prerequisite packages
if [ ${machine} = "Linux" ]
then
    echo -e "\nRunning in Linux\n"
    sudo apt-get update \
    && sudo apt-get install -y --no-install-recommends \
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
    echo -e "\nRunning in macOS\n"

    # Install Xcode
    xcode-select --install
    sudo xcodebuild -license accept
    /usr/bin/ruby -e "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install)"
    brew install curl git gcc
 
    # Install XQuartz
    if ! which 'xquartz' &>/dev/null
    then
        echo -e "\nXQuartz is not installed yet. Installing XQuartz ...\n"
        brew install caskroom/cask/brew-cask 2> /dev/null
        brew cask install xquartz
    else
        echo -e "\nXQuartz is already installed.\n"
    fi

fi
  
echo -e "\nFinished installing Prerequisites.\n"
echo "No user password will be asked from now on."
  
echo "Current directory: $PWD"
  
  
##############################
### SPACK
##############################
  
echo -e "\n****************************\nInstalling SPACK...\n****************************\n"
  
# Locate home directory
cd ~
echo "$PWD"
  
# Retrieve Spack for MoFEM
git clone https://github.com/likask/spack.git
  
# Initialise Spack environment variables:
. $HOME/spack/share/spack/setup-env.sh
  
# Add command to .bash_profile
echo ". $HOME/spack/share/spack/setup-env.sh" >> ~/.bash_profile
  
# Install packages required by Spack
spack bootstrap
  
echo -e "\nFinished installing Spack.\n"
  
echo "Current directory: $PWD"
  
##############################
### MoFEM USER MODULES
##############################
  
echo -e "\n********************************************************\n"
echo -e "Installing USER MODULE & FRACTURE MODULE..."
echo -e "\n********************************************************\n"
  
# Locate MoFEM installation directory
cd $MOFEM_INSTALL_DIR
echo "Current directory: $PWD"
  
# Install MoFEM packages
spack install  -j $NumberOfProcs mofem-fracture-module build_type=Release
  
# Activate fracture module
spack view --verbose symlink -i um_view mofem-fracture-module
 
# Export view and make view visible from any directory
export PATH=$PWD/um_view/bin:$PATH 
echo "export PATH=$PWD/um_view/bin:$PATH" >> ~/.bash_profile
 
echo -e "\nFinished installing MoFEM User Module and Fracture Module.\n"
 
 
# Test elasticity
cd $MOFEM_INSTALL_DIR/um_view/elasticity
echo "Current directory: $PWD"
./elasticity \
-my_file LShape.h5m \
-my_order 2 \
-ksp_type gmres \
-pc_type lu -pc_factor_mat_solver_package mumps \
-ksp_monitor 2>&1 | tee log
  
echo -e "\nFinished testing elasticity.\n"
 
# Test fracture crack propagation
cd $MOFEM_INSTALL_DIR/um_view/mofem_um_fracture_mechanics
echo "Current directory: $PWD"
./crack_propagation \
-my_file examples/analytical_bc/out_10.h5m \
-my_order 2 \
-my_ref 0 2>&1 | tee log
  
echo -e "\nFinished testing crack propagation.\n"