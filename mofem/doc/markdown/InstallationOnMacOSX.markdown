## Installation on OS X ##

Installation of MoFEM on MacOS X can be made as well within docker see details
[go to this page](doc/markdown/InstalltionWithDocker.markdown). Installation in
docker is faster where prebuild enviroent is downloaded directly form
[docker hub](https://hub.docker.com/r/likask/ubuntu_mofem/).

If you have any problems, feedback or would like to suggest corrections,
please email [cmatgu@googlegroups.com](mailto:cmatgu@googlegroups.com).
Before you start you need to install
XCode [https://developer.apple.com/xcode/downloads/](https://developer.apple.com/xcode/downloads/)
and XQuartz
[http://xquartz.macosforge.org/landing/](http://xquartz.macosforge.org/landing/).
Had to run XCode standalone to say OK to license agreement.

###1. Open a terminal

Create mofem install directory:
~~~~~~
export MOFEM_INSTALL_DIR=$HOME/mofem_installation
mkdir $MOFEM_INSTALL_DIR
~~~~~~

###2. Install libraries using homebrew

~~~~~~
ruby -e "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install)"
brew update
brew tap homebrew/science
brew tap homebrew/versions
brew install cmake
brew install --with-x11 gnuplot
brew install --with-graphviz doxygen
brew install wget
brew install boost
brew install homebrew/versions/gcc49 --with-fortran
brew install gnu-sed
~~~~~~

Set the default fortran compiler
~~~~~~
ln -s /usr/local/bin/gfortran-4.9 /usr/local/bin/gfortran`
~~~~~~

###3. Install PETSc and other libraries

~~~~~~
# Change to your $MOFEM_INSTALL_DIR
cd $MOFEM_INSTALL_DIR

git clone https://bitbucket.org/petsc/petsc.git
cd $MOFEM_INSTALL_DIR/petsc

# Fix PETSc vetsion
export PETSC_VERSION=3.6.3
git checkout tags/v$PETSC_VERSION

# Configure and compile petsc:
wget ftp://ftp.unidata.ucar.edu/pub/netcdf/netcdf-4.3.3.1.tar.gz
./configure --with-debugging=0 --download-fblaslapack=1 --download-superlu_dist=1 --download-metis=1 --download-parmetis=1 --download-hypre=1 --download-mumps=1 --download-scalapack=1 --download-blacs=1 --download-moab=1 --download-hdf5=1 --download-netcdf=$MOFEM_INSTALL_DIR/petsc/netcdf-4.3.3.1.tar.gz  --download-openmpi=1 --download-mumps=1
make PETSC_DIR=$PWD PETSC_ARCH=arch-darwin-c-opt all
~~~~~~

Note: PETSc is compiled with debugging switch off for efficiency. If you
develop code is recommended that you compile PETSc with debugging flag on in
addition. You can have two versions of MoFEM compiled, for debugging and
development and other version for larger calculations.

###4. Compile ADOL-C

~~~~~~
# Change to your $MOFEM_INSTALL_DIR
cd $MOFEM_INSTALL_DIR

# Get library, configure make and install
wget http://www.coin-or.org/download/source/ADOL-C/ADOL-C-2.5.2.zip
tar -xzf ADOL-C-2.5.2.zip
cd ADOL-C-2.5.2
./configure --prefix=$MOFEM_INSTALL_DIR/local
make install
~~~~~~

###5. Install other libraries

####5.3 Boost 1.57

~~~~~~
cd $MOFEM_INSTALL_DIR
# Download boost from http://www.boost.org/users/history/version_1_57_0.html
# and place it in mofem install directory ($MOFEM_INSTALL_DIR)
tar -xvvzf boost_1_57_0.tar.gz
cd boost_1_57_0
./bootstrap.sh --prefix=$MOFEM_INSTALL_DIR/local
./b2 install
~~~~~~

Note: It will take some time to build boost.

####5.2 TetGen

~~~~~~
cd $MOFEM_INSTALL_DIR
wget https://bitbucket.org/likask/mofem-joseph/downloads/tetgen1.5.0.tgz
tar -xvvzf tetgen1.5.0.tgz
cd tetgen1.5.0
cmake .
make
cp libtet.a lib/
~~~~~~

###6. Clone source code and install MoFEM library

~~~~~~
# Change to your $MOFEM_INSTALL_DIR
cd $MOFEM_INSTALL_DIR

# Cloning MoFEM sourcecode:
git clone https://bitbucket.org/likask/mofem-cephas.git mofem-cephas

# Make a build directory
mkdir $MOFEM_INSTALL_DIR/lib
cd $MOFEM_INSTALL_DIR/lib

# Configuring and compiling code:
cmake -DCMAKE_Fortran_COMPILER=/usr/local/bin/gfortran -DCMAKE_BUILD_TYPE=Release -DCMAKE_C_FLAGS="-Wall"  -DCMAKE_CXX_FLAGS="-Wall -Wno-bind-to-temporary-copy -Wno-overloaded-virtual" -DPETSC_DIR=$MOFEM_INSTALL_DIR/petsc/ -DPETSC_ARCH=arch-darwin-c-opt -DMOAB_DIR=$MOFEM_INSTALL_DIR/petsc/arch-darwin-c-opt/ -DADOL-C_DIR=$MOFEM_INSTALL_DIR/local/ -DTETGEN_DIR=$MOFEM_INSTALL_DIR/tetgen1.5.0 -DCMAKE_INSTALL_PREFIX=$MOFEM_INSTALL_DIR/users_modules $MOFEM_INSTALL_DIR/mofem-cephas/mofem

# Building code (assuming that you have computer with 4 cores):
make -j4 install

# Testing and publishing results on MoFEM CDashTesting WebPage:
ctest -D Experimental
~~~~~~

###7. Configuration, compilation and testing user modules

Before you start this version, change directory to install directory
~~~~~~
cd $MOFEM_INSTALL_DIR/users_modules
~~~~~~
Some elements still using some obsolete implementation which is gradually
removed. At this stage you need to install "obsolete" user modules:
~~~~~~
git clone https://bitbucket.org/likask/mofem_um_obsolete users_modules/obsolete
~~~~~~
List of some additional users modules is available on the main page.

~~~~~~
# Configuration:
cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_C_FLAGS="-Wall"  -DCMAKE_CXX_FLAGS="-Wall -Wno-bind-to-temporary-copy -Wno-overloaded-virtual" -DCMAKE_EXE_LINKER_FLAGS="-L$MOFEM_INSTALL_DIR/local/lib" users_modules

# Build:
make -j4
~~~~~~

###8. Dynamic Linked Shared Libraries

When executing a binary there may be a `dyld` related error. This is because the required dynamic libraries haven not been linked to the executable.

There are two ways of solving this problem:

1. Hack

  Add this to your ~/.bashrc:
  ~~~~~
  export MOFEM_INSTALL_DIR=$HOME/mofem_installation
  export PATH=$PATH:$MOFEM_INSTALL_DIR/petsc/arch-darwin-c-opt/bin
  export DYLD_LIBRARY_PATH=$DYLD_LIBRARY_PATH:$MOFEM_INSTALL_DIR/local/lib64:$MOFEM_INSTALL_DIR/local/lib
  ~~~~~

  Note: you'll need to reload the current session for this to take effect

2. Linking

  If you don't want to-do the above hack you can instead link the required dylib files to your executable. This requires `install_name_tool` and `otool`, which should come as part of Xcode.

  `otool` can check a binary for dynamicallly linked libraries e.g. `otool -L arc_length_nonlinear_elasticity`. Check the printout to see which libraries are not linked correctly. That can be a case when two or more version of the same library is installed
  in the system.

  Then change the libraries that are not linked correclty e.g. `install_name_tool -change libboost_program_options.dylib $MOFEM_INSTALL_DIR/local/lib/libboost_program_options.dylib arc_length_nonlinear_elasticity`

###9. Testing

~~~~~~
ctest -D Experimental
~~~~~~

Note that results of the tests are publish on MoFEM CDashTesting web page. If you do not wish to publish results remove option ``-D Experimental``
