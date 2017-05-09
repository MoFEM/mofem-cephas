##Installation on Ubuntu## {#install_ubuntu}

###1. Install packages

Following lines install minimum number of packages need to work with MoFEM:
~~~~~~
apt-get update
apt-get install -y wget valgrind git g++ gdb m4 automake libsigsegv2 build-essential libibverbs-dev libblas-dev gfortran libatlas-dev libatlas-base-dev libhdf5-openmpi-dev libjpeg-dev graphviz doxygen cmake gnuplot pstack ca-certificates python libadolc-dev bison flex libx11-dev libboost-all-dev xauth xterm libtool
~~~~~~

###1. Open a terminal

Create mofem install director:
~~~~~~
export MOFEM_INSTALL_DIR=$HOME/mofem_installation
mkdir $MOFEM_INSTALL_DIR
~~~~~~

###2. Install PETSc and other libraries

~~~~~~
# Change to your $MOFEM_INSTALL_DIR
cd $MOFEM_INSTALL_DIR

# Colone PETSc repository:

git clone https://bitbucket.org/petsc/petsc.git
cd $MOFEM_INSTALL_DIR/petsc

# Fix PETSc version
export PETSC_VERSION=3.6.3
git checkout tags/v$PETSC_VERSION

# Configure and compile petsc:
wget ftp://ftp.unidata.ucar.edu/pub/netcdf/netcdf-4.3.3.1.tar.gz
./configure --with-mpi=1 --with-debugging=0 --download-superlu_dist=1 --download-metis=1 --download-parmetis=1 --download-hypre=1 --download-mumps=1 --download-scalapack=1 --download-zoltan=1 --download-blacs=1 --download-moab=1 --download-ptscotch=1 --with-hdf5=1 --with-hdf5-dir=/usr --download-netcdf=$MOFEM_INSTALL_DIR/petsc/netcdf-4.3.3.1.tar.gz --with-shared-libraries=1
make PETSC_DIR=$PWD PETSC_ARCH=arch-linux2-c-opt all
~~~~~~

Note: PETSc is compiled with debugging switch off for efficiency. If you
develop code is recommended that you compile PETSc with debugging flag on in
addition. You can have two versions of MoFEM compiled, for debugging and
development and other version for larger calculations.

###3. Install TetGen and other libraries

####3.1 TetGen

~~~~~~
cd $MOFEM_INSTALL_DIR
wget https://bitbucket.org/likask/mofem-joseph/downloads/tetgen1.5.0.tgz
tar -xvvzf tetgen1.5.0.tgz
cd tetgen1.5.0
cmake .
make
cp libtet.a lib/
~~~~~~

###4. Clone source code and install core MoFEM library

~~~~~~
# Change to your $MOFEM_INSTALL_DIR
cd $MOFEM_INSTALL_DIR

# Cloning MoFEM source code:
git clone https://bitbucket.org/likask/mofem-cephas.git mofem-cephas

# Make a build directory
mkdir $MOFEM_INSTALL_DIR/lib
cd $MOFEM_INSTALL_DIR/lib

# Configuring and compiling code:
cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_C_FLAGS="-Wall"  -DCMAKE_CXX_FLAGS="-Wall" -DPETSC_DIR=$MOFEM_INSTALL_DIR/petsc/ -DPETSC_ARCH=arch-linux2-c-opt -DMOAB_DIR=$MOFEM_INSTALL_DIR/petsc/arch-linux2-c-opt/ -DADOL-C_DIR=$MOFEM_INSTALL_DIR/local/ -DTETGEN_DIR=$MOFEM_INSTALL_DIR/tetgen1.5.0 -DCMAKE_INSTALL_PREFIX=$MOFEM_INSTALL_DIR/users_modules $MOFEM_INSTALL_DIR/mofem-cephas/mofem

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
removed. At this stage you need to install "obsolete" user modules:
~~~~~~
git clone https://bitbucket.org/likask/mofem_um_obsolete users_modules/obsolete
~~~~~~
List of some additional users modules is available on the main page.

~~~~~~
# Configuration:
cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_C_FLAGS="-Wall" -DCMAKE_CXX_FLAGS="-Wall"  users_modules

# Build:
make -j4

# Testing:
ctest -D Experimental
~~~~~~

Note that results of the test are publish on MoFEM CDashTesting web page. If you do not like publish results pleas remove option ``-D Experimental``
