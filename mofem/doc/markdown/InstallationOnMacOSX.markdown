## Installation on OS X ## {#install_macosx}

If you have any problems, feedback or would like to suggest corrections,
please email [mofem-group@googlegroups.com](https://groups.google.com/forum/#!forum/mofem-group).
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
ln -s /usr/local/bin/gfortran-4.9 /usr/local/bin/gfortran
~~~~~~

###3. Install PETSc and other libraries

~~~~~~
# Change to your $MOFEM_INSTALL_DIR
cd $MOFEM_INSTALL_DIR

git clone https://bitbucket.org/petsc/petsc.git
cd $MOFEM_INSTALL_DIR/petsc

# Fix PETSc version
export PETSC_VERSION=3.8.4
git checkout tags/v$PETSC_VERSION

# Configure and compile petsc:
./configure --with-debugging=0 --download-fblaslapack=1 --download-superlu_dist=1 \
--download-metis=1 --download-parmetis=1 --download-hypre=1 --download-mumps=1 \
--download-scalapack=1 --download-blacs=1 --download-moab=1 --download-hdf5=1 \
--download-netcdf=1 --download-mumps=1 --download-openmpi=1 && \
make PETSC_DIR=$PWD PETSC_ARCH=arch-darwin-c-opt all

# Add path to petsc binaries, you can add that line to .bashrc
export PATH=$MOFEM_INSTALL_DIR/petsc/$PETSC_ARCH/bin:$PATH
~~~~~~

Note: PETSc is compiled with debugging switch off for efficiency. If you
develop code is recommended that you compile PETSc with debugging flag on in
addition. You can have two versions of MoFEM compiled, for debugging and
development and other version for larger calculations.

###4. Clone source code and install MoFEM library

~~~~~~
# Change to your $MOFEM_INSTALL_DIR
cd $MOFEM_INSTALL_DIR

# Cloning MoFEM sourcecode:
git clone https://bitbucket.org/likask/mofem-cephas.git mofem-cephas

# Make a build directory
mkdir $MOFEM_INSTALL_DIR/lib
cd $MOFEM_INSTALL_DIR/lib

# Configuring and compiling code:
cmake -DCMAKE_BUILD_TYPE=Release \
 -DCMAKE_C_FLAGS="-Wall" \
 -DCMAKE_CXX_FLAGS="-Wall -Wno-bind-to-temporary-copy -Wno-overloaded-virtual" \
 -DPETSC_DIR=$MOFEM_INSTALL_DIR/petsc/ -DPETSC_ARCH=arch-darwin-c-opt \
 -DMOAB_DIR=$MOFEM_INSTALL_DIR/petsc/arch-darwin-c-opt/  \
 -DWITH_ADOL-C=1 -DWITH_TETGEN=1 -DWITH_MED=1 \
 -DCMAKE_INSTALL_PREFIX=$MOFEM_INSTALL_DIR/users_modules \
 $MOFEM_INSTALL_DIR/mofem-cephas/mofem

# Building code (assuming that you have computer with 4 cores):
make -j4 install

# Testing and publishing results on MoFEM CDashTesting WebPage:
ctest -D Experimental
~~~~~~

###5. Configuration, compilation and testing user modules

Before you start this version, change directory to install directory
~~~~~~
cd $MOFEM_INSTALL_DIR/users_modules
~~~~~~
Some elements still using some obsolete implementation which is gradually
removed. At this stage you need to install "obsolete" for user modules:
~~~~~~
git clone https://bitbucket.org/likask/mofem_um_obsolete users_modules/obsolete
~~~~~~
List of some additional users modules is available on the main page.

~~~~~~
# Configuration:
cmake -DCMAKE_BUILD_TYPE=Release \
-DCMAKE_C_FLAGS="-Wall"  \
-DCMAKE_CXX_FLAGS="-Wall -Wno-bind-to-temporary-copy -Wno-overloaded-virtual" \
-DCMAKE_EXE_LINKER_FLAGS="-L$MOFEM_INSTALL_DIR/local/lib" users_modules

# Build:
make -j4
~~~~~~

###6. Testing

~~~~~~
ctest -D Experimental
~~~~~~

Note that results of the tests are publish on MoFEM CDashTesting web page. If you 
do not wish to publish results remove option ``-D Experimental``
