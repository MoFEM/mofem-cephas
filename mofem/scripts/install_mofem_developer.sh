#!/bin/bash
# This script automises installation of MoFEM developer version: Release & Debug
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
# Note: Installation script changes .bashrc on Ubuntu or .bash_profile on Mac.
# Please inspect the file after installation.

echo "Start time: $(date +"%T")"

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
    NumberOfProcs=$(nproc)
elif [ ${machine} = "Mac" ]
then
    NumberOfProcs=$(($(sysctl -n hw.ncpu)/2))
fi
# max(NumberOfProcs, 1)
NumberOfProcs=$(( NumberOfProcs > 1 ? NumberOfProcs : 1 ))
  
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
    echo "Downloading spack ..."
    mkdir -p $SPACK_ROOT_DIR &&\
    curl -s -L https://api.github.com/repos/likask/spack/tarball/mofem \
    | tar xzC $SPACK_ROOT_DIR --strip 1
    echo -e "Done.\n"
  else 
    mkdir -p $SPACK_ROOT_DIR &&\
    tar xzf $PWD/spack.tgz -C $SPACK_ROOT_DIR --strip 1
  fi

  # Initialise Spack environment variables:
  . $SPACK_ROOT_DIR/share/spack/setup-env.sh
  # Add command to configuration file .bashrc on Ubuntu or .bash_profile on Mac
  if [ ${machine} = "Linux" ]
  then
    echo ". $SPACK_ROOT_DIR/share/spack/setup-env.sh" >> ~/.bashrc
  elif [ ${machine} = "Mac" ]
  then
    echo ". $SPACK_ROOT_DIR/share/spack/setup-env.sh" >> ~/.bash_profile
  fi

  # Download mirror
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

  # Install packages required by Spack
  spack bootstrap
fi
 
echo -e "\nFinished installing Spack.\n"
  
echo "Current directory: $PWD"
  
########################################
### MoFEM CORE LIBRARY & USER MODULES
########################################
  
echo -e "\n********************************************************\n"
echo -e "Installing CORE LIBRARY - Release version ..."
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

# Clone MoFEM Fracture Module
cd $MOFEM_INSTALL_DIR/mofem-cephas/mofem/users_modules
git clone -b develop https://bitbucket.org/likask/mofem_um_fracture_mechanics.git

# Installation of core library
cd $MOFEM_INSTALL_DIR
mkdir -p lib_release
cd lib_release

spack install --only dependencies mofem-cephas
spack setup mofem-cephas@develop copy_user_modules=False build_type=Release

echo -e "\n----------------------------\n"
echo -e "CORE LIBRARY - Release version: spconfig ..."
echo -e "\n----------------------------\n"

./spconfig.py -DMOFEM_BUILD_TESTS=ON $MOFEM_INSTALL_DIR/mofem-cephas/mofem

echo -e "\n----------------------------\n"
echo -e "CORE LIBRARY - Release version: make -j $NumberOfProcs ..."
echo -e "\n----------------------------\n"

make -j $NumberOfProcs

# Run the tests on core library
echo -e "\n----------------------------\n"
echo -e "CORE LIBRARY - Release version: ctest ..."
echo -e "\n----------------------------\n"

spack load cmake
ctest -D Experimental

# Install the library
echo -e "\n----------------------------\n"
echo -e "CORE LIBRARY - Release version: make install ..."
echo -e "\n----------------------------\n"

make install

echo -e "\nFinished installing and testing the Core Library.\n"

echo -e "\n********************************************************\n"
echo -e "Installing USER MODULES - Release version ..."
echo -e "\n********************************************************\n"

# Installation of user modules
cd $MOFEM_INSTALL_DIR

mkdir -p um
cd um/
spack view --verbose symlink -i um_view mofem-cephas@develop build_type=Release 
export PATH=$PWD/um_view/bin:$PATH
# Add PATH to .bashrc on Ubuntu or .bash_profile on Mac
if [ ${machine} = "Linux" ]
then
  echo "export PATH=$PWD/um_view/bin:\$PATH" >> ~/.bashrc
elif [ ${machine} = "Mac" ]
then
  echo "export PATH=$PWD/um_view/bin:\$PATH" >> ~/.bash_profile
fi

mkdir -p build_release
cd build_release/

MOFEM_CEPHAS_SPACK_ID=`spack find -lv mofem-cephas | grep mofem-cephas@develop | grep build_type=Release | awk '{print $1}'| tail -n 1`
echo "mofem-cephas id for release: $MOFEM_CEPHAS_SPACK_ID"

spack setup mofem-users-modules@develop \
    copy_user_modules=False build_type=Release \
    ^/$MOFEM_CEPHAS_SPACK_ID

echo -e "\n----------------------------\n"
echo -e "USER MODULE - Release version: spconfig ..."
echo -e "\n----------------------------\n"

./spconfig.py -DMOFEM_UM_BUILD_TESTS=ON -DFM_VERSION_MAJOR=0 -DFM_VERSION_MINOR=0 -DFM_VERSION_BUILD=0 -DMOFEM_DIR=../um_view \
    $MOFEM_INSTALL_DIR/mofem-cephas/mofem/users_modules

echo -e "\n----------------------------\n"
echo -e "USER MODULE - Release version: make -j $NumberOfProcs ..."
echo -e "\n----------------------------\n"

make -j $NumberOfProcs

# Run the tests on user modules
echo -e "\n----------------------------\n"
echo -e "USER MODULE - Release version: ctest ..."
echo -e "\n----------------------------\n"

ctest -D Experimental

# Install the user module
echo -e "\n----------------------------\n"
echo -e "USER MODULE - Release version: make install ..."
echo -e "\n----------------------------\n"

make install

echo -e "\nFinished installing and testing the User Module - Release version.\n"


# ************************************************************************
# DEBUG VERSION
# ************************************************************************

echo -e "\n********************************************************\n"
echo -e "Installing MoFEM - DEBUG VERSION..."
echo -e "\n********************************************************\n"
  

########################################
### MoFEM CORE LIBRARY & USER MODULES
########################################

# Locate MoFEM installation directory
cd $MOFEM_INSTALL_DIR

# Installation of core library
mkdir -p lib_debug
cd lib_debug

spack setup mofem-cephas@develop copy_user_modules=False build_type=Debug

echo -e "\n----------------------------\n"
echo -e "CORE LIBRARY - Debug version: spconfig ..."
echo -e "\n----------------------------\n"

./spconfig.py -DMOFEM_BUILD_TESTS=ON -DCMAKE_EXPORT_COMPILE_COMMANDS=1 $MOFEM_INSTALL_DIR/mofem-cephas/mofem

echo -e "\n----------------------------\n"
echo -e "CORE LIBRARY - Debug version: make -j $NumberOfProcs ..."
echo -e "\n----------------------------\n"

make -j $NumberOfProcs

# Run the tests on core library
echo -e "\n----------------------------\n"
echo -e "CORE LIBRARY - Debug version: ctest ..."
echo -e "\n----------------------------\n"

ctest -D Experimental

# Install the library
echo -e "\n----------------------------\n"
echo -e "CORE LIBRARY - Debug version: make install ..."
echo -e "\n----------------------------\n"

make install

echo -e "\nFinished installing and testing the Core Library - Debug version.\n"

echo -e "\n********************************************************\n"
echo -e "Installing USER MODULES - Debug version ..."
echo -e "\n********************************************************\n"

# Installation of user modules
cd $MOFEM_INSTALL_DIR

# mkdir um
cd um/

spack view --verbose symlink -i um_view_debug mofem-cephas@develop build_type=Debug 

mkdir -p build_debug
cd build_debug/

MOFEM_CEPHAS_SPACK_ID=`spack find -lv mofem-cephas | grep mofem-cephas@develop | grep build_type=Debug | awk '{print $1}'| tail -n 1`
echo "mofem-cephas id for debug: $MOFEM_CEPHAS_SPACK_ID"

spack setup mofem-users-modules@develop \
    copy_user_modules=False build_type=Debug \
    ^/$MOFEM_CEPHAS_SPACK_ID

echo -e "\n----------------------------\n"
echo -e "USER MODULE: spconfig - Debug version ..."
echo -e "\n----------------------------\n"

./spconfig.py -DMOFEM_UM_BUILD_TESTS=ON -DCMAKE_EXPORT_COMPILE_COMMANDS=1  -DFM_VERSION_MAJOR=0 -DFM_VERSION_MINOR=0 -DFM_VERSION_BUILD=0 -DMOFEM_DIR=../um_view_debug \
    $MOFEM_INSTALL_DIR/mofem-cephas/mofem/users_modules

echo -e "\n----------------------------\n"
echo -e "USER MODULE - Debug version: make -j $NumberOfProcs ..."
echo -e "\n----------------------------\n"

make -j $NumberOfProcs

# Run the tests on user modules
echo -e "\n----------------------------\n"
echo -e "USER MODULE - Debug version: ctest ..."
echo -e "\n----------------------------\n"

ctest -D Experimental

# Install the user module
echo -e "\n----------------------------\n"
echo -e "USER MODULE - Debug version: make install ..."
echo -e "\n----------------------------\n"

make install

echo -e "\nFinished installing and testing the User Module - Debug version.\n"

echo "End time: $(date +"%T")"