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
        /usr/bin/ruby -e "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install)"
    else 
        echo -e "\nHomebrew is already installed."
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
      xcode-select --install
      sudo xcodebuild -license accept
      /bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebreinstall/HEAD/install.sh)" 
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

  # Remove .spack directory in $HOME from previous installation (if any)
  if [ -d "$HOME/.spack" ]; then
    mv $HOME/.spack $HOME/.spack_old
  fi

  echo "Cloning spack ..."
  git clone -b master https://github.com/likask/spack.git
  echo -e "Done.\n"

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
      echo "Downloading mirror of spack packages for MoFEM..."
      mkdir -p $SPACK_MIRROR_DIR && \
      curl -s -L http://mofem.eng.gla.ac.uk/mofem/downloads/mirror_v0.16.tar.gz \
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
  spack compiler find
  spack external find

  # Set fortran compiler to version 9
  if [ ${machine} = "Mac" ]
  then
    sed 's/gfortran$/gfortran-9/g' $HOME/.spack/darwin/compilers.yaml
  fi

else
  spack external find  
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
lse 
  echo -e "\nMoFEM source directory is found"
fi

# Installation of core library
cd $MOFEM_INSTALL_DIR

echo -e "\n----------------------------\n"
echo -e "CORE LIBRARY - Install depenencies ..."
echo -e "\n----------------------------\n"

spack install --only dependencies mofem-cephas ^petsc+X

echo -e "\n----------------------------\n"
echo -e "CORE LIBRARY - Release version ..."
echo -e "\n----------------------------\n"

spack dev-build \
  --source-path $MOFEM_INSTALL_DIR/mofem-cephas \
  --keep-prefix \
  --test root \
  mofem-cephas@develop~copy_user_modules build_type=RelWithDebInfo ^petsc+X

echo -e "\n********************************************************\n"
echo -e "Installing USER MODULES - Release version ..."
echo -e "\n********************************************************\n"

# Get mofem-cephas spack hash for release version
TODAY=`date +%F` 
MOFEM_CEPHAS_HASH=`spack find -lv --start-date $TODAY | grep mofem-cephas@develop | grep RelWithDebInfo | awk '{print $1}'` 
echo "mofem-cephas id for release: $MOFEM_CEPHAS_HASH"

spack dev-build \
  --test root  \
  --source-path $MOFEM_INSTALL_DIR/mofem-cephas/mofem/users_modules \
  mofem-users-modules@develop build_type=RelWithDebInfo \
  ^/$MOFEM_CEPHAS_HASH

TODAY=`date +%F` 
MOFEM_UN_HASH=`spack find -lv --start-date $TODAY | grep mofem-users-modulesdevelop | grep RelWithDebInfo | awk '{print $1}'` 
echo "mofem-users-modules id for release: $MOFEM_UN_HASH"

# ************************************************************************
# DEBUG VERSION
# ************************************************************************

echo -e "\n********************************************************\n"
echo -e "Installing MoFEM - DEBUG VERSION..."
echo -e "\n********************************************************\n"
  

########################################
### MoFEM CORE LIBRARY & USER MODULES
########################################

echo -e "\n----------------------------\n"
echo -e "CORE LIBRARY - Debug version ..."
echo -e "\n----------------------------\n"

spack dev-build \
  --source-path $MOFEM_INSTALL_DIR/mofem-cephas \
  --keep-prefix \
  --test root \
  mofem-cephas@develop~copy_user_modules build_type=Debug ^petsc+X

echo -e "\n********************************************************\n"
echo -e "Installing USER MODULES - Debug version ..."
echo -e "\n********************************************************\n"

# Get mofem-cephas spack hash for debug version
TODAY=`date +%F` 
MOFEM_CEPHAS_HASH=`spack find -lv --start-date $TODAY | grep mofem-cephas@develop | grep Debug | awk '{print $1}'` 
echo "mofem-cephas id for debug: $MOFEM_CEPHAS_HASH"

echo -e "\n----------------------------\n"
echo -e "USER MODULE - Debug version ..."
echo -e "\n----------------------------\n"

spack dev-build \
  --test root  \
  --source-path $MOFEM_INSTALL_DIR/mofem-cephas/mofem/users_modules \
  mofem-users-modules@develop build_type=Debug \
  ^/$MOFEM_CEPHAS_HASH

TODAY=`date +%F` 
MOFEM_UN_HASH=`spack find -lv --start-date $TODAY | grep mofem-users-modulesdevelop | grep Debug | awk '{print $1}'` 
echo "mofem-users-modules id for debug: $MOFEM_UN_HASH"

echo -e "\nFinished installing and testing the User Module - Debug version.\n"
cd $MOFEM_INSTALL_DIR

echo "End time: $(date +"%T")"