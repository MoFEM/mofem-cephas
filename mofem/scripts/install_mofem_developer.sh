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
# Note: Installation script changes .bash_profile. Inspect that file after installation.

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
    if ! which 'brew' &>/dev/null
    then
        xcode-select --install
        sudo xcodebuild -license accept
        /usr/bin/ruby -e "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install)"
    else 
        echo -e "\nHomebrew installed"
    fi
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
echo -e "\nNo user password will be asked from now on.\n"
  
echo "Current directory: $PWD"
  
  
##############################
### SPACK
##############################
  
echo -e "\n****************************\nInstalling SPACK...\n****************************\n"
  
# Locate home directory
cd $MOFEM_INSTALL_DIR
echo "$PWD"

SPACK_ROOT_DIR=$MOFEM_INSTALL_DIR/spack
SPACK_MIRROR_DIR=$MOFEM_INSTALL_DIR/mofem_mirror

# Retrieve Spack for MoFEM
if [ ! -d "$SPACK_ROOT_DIR" ]; then
  if [ ! -f "$PWD/spack.tgz" ]; then
    echo "Download spack"
    mkdir -p $SPACK_ROOT_DIR &&\
    curl -s -L https://api.github.com/repos/likask/spack/tarball/mofem \
    | tar xzC $SPACK_ROOT_DIR --strip 1
  else 
    mkdir -p $SPACK_ROOT_DIR &&\
    tar xzf $PWD/spack.tgz -C $SPACK_ROOT_DIR --strip 1
  fi
fi

# Download mirror
if [ ! -d "$SPACK_ROOT_DIR" ]; then
  if [ ! -f "$PWD/mirror.tgz" ]; then
    echo "Download spack mofem mirror"
    mkdir -p mofem_mirror && \
    curl -s -L https://bitbucket.org/likask/mofem-cephas/downloads/mirror_v0.9.0.tar.gz \
    | tar xzC $SPACK_MIRROR_DIR --strip 1
  else 
    mkdir -p $SPACK_MIRROR_DIR && \
    tar xzf $PWD/mirror.tgz -C $SPACK_MIRROR_DIR  --strip 1
  fi
fi
 
# Initialise Spack environment variables:
. $SPACK_ROOT_DIR/share/spack/setup-env.sh

# Add mirror
spack mirror add mofem_mirror $SPACK_MIRROR_DIR
  
# Add command to configuration file .bash_profile
echo ". $SPACK_ROOT_DIR/share/spack/setup-env.sh" >> ~/.bash_profile
  
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
if [ ! -d "$PWD/mofem-cephas" ]; then
  git clone -b develop --recurse-submodules https://bitbucket.org/likask/mofem-cephas.git mofem-cephas
else 
  echo -e "\nMoFEM source directory is found"
fi

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
echo "export PATH=$PWD/um_view/bin:$PATH" >> ~/.bash_profile

mkdir build
cd build/
spack setup mofem-users-modules@develop \
    copy_user_modules=False build_type=RelWithDebInfo \
    ^mofem-cephas@develop 

echo -e "\n----------------------------\n"
echo -e "USER MODULE: spconfig ..."
echo -e "\n----------------------------n"

./spconfig.py -DMOFEM_UM_BUILD_TESTS=ON -DMOFEM_DIR=../um_view \
    $MOFEM_INSTALL_DIR/mofem-cephas/mofem/users_modules

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