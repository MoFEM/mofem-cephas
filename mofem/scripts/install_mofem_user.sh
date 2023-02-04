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
# Note: Installation script changes .bashrc on Ubuntu or .bash_profile on Mac.
# Please inspect the file after installation.

# Only for debugging
#set -x 

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
    pkgconf \
    cmake \
    git \
    python \
    python3-distutils \
    unzip \
    ssh \
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
        /bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"
    else 
        echo -e "\nHomebrew installed"
    fi

    brew install \
      curl \
      git \
      python \
      gcc@9 \
      cmake \
      autoconf \
      automake \
      libtool \
      doxygen \
      pkg-config

    # Install XQuartz
    if ! which 'xquartz' &>/dev/null
    then
        echo -e "\nXQuartz is not installed yet. Installing XQuartz ...\n"
        brew install --cask xquartz
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
  
cd $MOFEM_INSTALL_DIR
echo "$PWD"

SPACK_ROOT_DIR=$MOFEM_INSTALL_DIR/spack
SPACK_MIRROR_DIR=$MOFEM_INSTALL_DIR/mofem_mirror

# Retrieve Spack for MoFEM
if [ ! -d "$SPACK_ROOT_DIR" ]; then

  # Remove .spack directory in $HOME from previous installation (if any)
  if [ -d "$HOME/.spack" ]; then
    mv $HOME/.spack $HOME/.spack_old
  fi

  if [ ! -f "$PWD/spack.tgz" ]; then
    echo "Downloading spack ..."
    mkdir -p $SPACK_ROOT_DIR &&\
    curl -s -L https://api.github.com/repos/likask/spack/tarball/master \
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
    echo ". $SPACK_ROOT_DIR/share/spack/setup-env.sh" >> ~/.zshrc
  fi
  
  # Download mirror
  if [ ! -d "$SPACK_MIRROR_DIR" ]; then
    if [ ! -f "$PWD/mirror.tgz" ]; then
      # echo "Downloading mirror of spack packages for MoFEM..."
      # mkdir -p $SPACK_MIRROR_DIR && \
      # curl -s -L http://mofem.eng.gla.ac.uk/mofem/downloads/mirror_v0.16.tar.gz \
      # | tar xzC $SPACK_MIRROR_DIR --strip 1
      # echo -e "Done.\n"
    else 
      # mkdir -p $SPACK_MIRROR_DIR && \
      # tar xzf $PWD/mirror.tgz -C $SPACK_MIRROR_DIR  --strip 1
    fi
  fi
 
  # FIXME: We do not have mirror build for most recent version
  # Add mirror
  # spack mirror remove mofem_mirror 2> /dev/null
  # spack mirror add mofem_mirror $SPACK_MIRROR_DIR

  # Install packages required by Spack
  spack compiler find
  spack external find

  # # Set fortran compiler to version 9
  # # Set fortran compiler to version 9
  # if [ ${machine} = "Mac" ]
  # then
  #   sed 's/gfortran$/gfortran-9/g' $HOME/.spack/darwin/compilers.yaml
  # fi

fi
 
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
spack install --only dependencies mofem-cephas ^petsc+X
spack install --test all mofem-fracture-module build_type=Release ^petsc+X

# Activate fracture module
spack view --verbose symlink -i um_view mofem-fracture-module
 
# Export view and make view visible from any directory
export PATH=$PWD/um_view/bin:$PATH 

echo "End time: $(date +"%T")"