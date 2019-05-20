#!/bin/bash
# This script automises installation of MoFEM developer version (RelWithDebInfo)
# Full description of the installation process using spack can be found at:
# http://mofem.eng.gla.ac.uk/mofem/html/install_spack.html
# The script should work for both Ubuntu and macOS platforms
# The last known working platforms are Ubuntu (18.04) and macOS (10.13, 10.14)
# The installation may take two to three hours
#
# Usage:
#       1. Copy install_mofem_developer.sh to the directory where MoFEM will be installed
#       2. Run install_mofem_developer.sh from the command line
#
  
##############################
# INITIALISATION
##############################
  
# Setup installation directory
pwd
export MOFEM_INSTALL_DIR="$PWD/mofem_install"
#export MOFEM_INSTALL_DIR="$HOME"
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
    echo "Running in Linux"
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
    echo "Running in macOS"
    ## Assuming Xcode already installed
    xcode-select --install
    sudo xcodebuild -license accept
    /usr/bin/ruby -e "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install)"
    brew install curl git gcc
 
fi
  
echo -e "\nFinished installing Prerequisites.\n"
echo -e "\nNo user password will be asked from now on.\n"
  
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
  
# Add command to .profile
echo ". $HOME/spack/share/spack/setup-env.sh" >> ~/.profile
  
# Install packages required by Spack
spack bootstrap
  
echo -e "\nFinished installing Spack.\n"
  
echo "Current directory: $PWD"
  
########################################
### MoFEM CORE LIBRARY & USER MODULES
########################################
  
echo -e "\n********************************************************\n"
echo -e "Installing CORE LIBRARY ..."
echo -e "\n********************************************************\n"
  
# Locate MoFEM installation directory
cd $MOFEM_INSTALL_DIR
echo "Current directory: $PWD"
  
# Clone MoFEM core library
git clone -b develop --recurse-submodules https://bitbucket.org/likask/mofem-cephas.git mofem-cephas

# Installation of core library
mkdir lib
cd lib

spack install --only dependencies mofem-cephas 
spack setup mofem-cephas@develop copy_user_modules=False build_type=RelWithDebInfo

echo -e "\n----------------------------\n"
echo -e "CORE LIBRARY: spconfig ..."
echo -e "\n----------------------------n"

./spconfig.py -DMOFEM_BUILD_TESTS=ON $MOFEM_INSTALL_DIR/mofem-cephas/mofem

echo -e "\n----------------------------\n"
echo -e "CORE LIBRARY: make -j $NumberOfProcs ..."
echo -e "\n----------------------------\n"

make -j $NumberOfProcs

# Run the tests on core library
echo -e "\n----------------------------\n"
echo -e "CORE LIBRARY: ctest ..."
echo -e "\n----------------------------\n"

spack load cmake
ctest -D Experimental

# Install the library
echo -e "\n----------------------------\n"
echo -e "CORE LIBRARY: make install ..."
echo -e "\n----------------------------\n"

make install

echo -e "\nFinished installing and testing the Core Library.\n"


echo -e "\n********************************************************\n"
echo -e "Installing USER MODULES ..."
echo -e "\n********************************************************\n"

# Installation of user modules
cd $MOFEM_INSTALL_DIR

mkdir um
cd um/
spack view --verbose symlink -i um_view mofem-cephas@develop build_type=RelWithDebInfo 
export PATH=$PWD/um_view/bin:$PATH

mkdir build
cd build/
spack setup mofem-users-modules@develop copy_user_modules=False build_type=RelWithDebInfo \
    ^mofem-cephas@develop copy_user_modules=False build_type=RelWithDebInfo

echo -e "\n----------------------------\n"
echo -e "USER MODULE: spconfig ..."
echo -e "\n----------------------------n"

./spconfig.py -DMOFEM_UM_BUILD_TESTS=ON -DMOFEM_DIR=../um_view $MOFEM_INSTALL_DIR/mofem-cephas/mofem/users_modules

echo -e "\n----------------------------\n"
echo -e "USER MODULE: make -j $NumberOfProcs ..."
echo -e "\n----------------------------\n"

make -j $NumberOfProcs

# Run the tests on user modules
echo -e "\n----------------------------\n"
echo -e "USER MODULE: ctest ..."
echo -e "\n----------------------------\n"

spack load cmake
ctest -D Experimental

# Install the user module
echo -e "\n----------------------------\n"
echo -e "USER MODULE: make install ..."
echo -e "\n----------------------------\n"

make install

echo -e "\nFinished installing and testing the User Module.\n"