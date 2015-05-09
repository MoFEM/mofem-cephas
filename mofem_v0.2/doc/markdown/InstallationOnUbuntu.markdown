##Installation on Ubuntu##


###1. Install packages

Following lines install minimum number of packages need to work with MoFEM:
~~~~~~
apt-get update
apt-get install -y wget valgrind git g++ gdb m4 automake libsigsegv2 build-essential libibverbs-dev libblas-dev gfortran libatlas-dev libatlas-base-dev libhdf5-openmpi-dev libjpeg-dev graphviz doxygen cmake gnuplot pstack ca-certificates python libadolc-dev bison flex libx11-dev libboost-all-dev xauth xterm
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
cd $HOME/$MOFEM_INSTALL_DIR

# Colone PETSc repository:
export PETSC_VERSION 3.5.3
git clone https://bitbucket.org/petsc/petsc.git
cd /opt/petsc
git checkout tags/v$PETSC_VERSION


# Configure and compile petsc:
wget ftp://ftp.unidata.ucar.edu/pub/netcdf/netcdf-4.3.3.1.tar.gz
./configure --with-mpi=1 --with-debugging=0 --download-superlu_dist=1 --download-metis=1 --download-parmetis=1 --download-hypre=1 --download-mumps=1 --download-scalapack=1 --download-zoltan=1 --download-blacs=1 --download-moab=1 --download-ptscotch=1 --with-hdf5=1 --with-hdf5-dir=/usr --download-netcdf=/opt/petsc/netcdf-4.3.3.1.tar.gz --with-shared-libraries=1
make PETSC_DIR=/opt/petsc PETSC_ARCH=arch-linux2-c-opt all

# Building code (assuming that you have computer with 4 cores):
make -j4 install

# Testing and publishing results on MoFEM CDashTesting WebPage:
ctest -D Experimental
~~~~~~

Note: PETSc is compiled with debugging switch off for efficiency. If you
develop code is recommened that you compile petsce with debugging flag on.

Note: Version pets3.5.3 is checkout. MoFEM should work with development
version of petsc if mofem is configured with `` -DCMAKE_C_FLAGS="-Wall
-DPETSC_DEV" -DCMAKE_CXX_FLAGS="-Wall -DPETSC_DEV"``.

###3. Install TetGem and other libraries

####3.1 TetGen

~~~~~~
cd $HOME/$MOFEM_INSTALL_DIR
wget https://bitbucket.org/likask/mofem-joseph/downloads/tetgen1.5.0.tgz
tar -xvvzf tetgen1.5.0.tgz
cd tetgen1.5.0
cmake .
make
cp libtet.a lib/
~~~~~~

###4. Clone sourcecode and install core MoFEM library

~~~~~~
# Change to your $MOFEM_INSTALL_DIR
cd $HOME/$MOFEM_INSTALL_DIR

# Cloning MoFEM sourcecode:
git clone https://bitbucket.org/likask/mofem-cephas.git mofem-cephas

# Make a build directory
mkdir $MOFEM_INSTALL_DIR/lib
cd $MOFEM_INSTALL_DIR/lib

# Configuring and compiling code:
cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_C_FLAGS="-Wall"  -DCMAKE_CXX_FLAGS="-Wall" -DPETSC_DIR=$MOFEM_INSTALL_DIR/petsc/ -DPETSC_ARCH=arch-darwin-c-opt -DMOAB_DIR=$MOFEM_INSTALL_DIR/petsc/arch-darwin-c-debug/ -DADOL-C_DIR=$MOFEM_INSTALL_DIR/local/ -DCMAKE_INSTALL_PREFIX=$MOFEM_INSTALL_DIR/user_modules $MOFEM_INSTALL_DIR/mofem-cephas/mofem_v0.2

# Building code (assuming that you have computer with 4 cores):
make -j4 install

# Testing and publishing results on MoFEM CDashTesting WebPage:
ctest -D Experimental
~~~~~~

###5. Configuration, compilation and testing user modules

~~~~~~
# Change to your $MOFEM_INSTALL_DIR
cd $HOME/$MOFEM_INSTALL_DIR/user_modules

# Configuration:
cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_C_FLAGS="-Wall" -DCMAKE_CXX_FLAGS="-Wall" $MOFEM_INSTALL_DIR/user_modules

# Build:
make -j4

# Testing:
cmake -D Experimental
~~~~~~

Note that results of the test are publish on MoFEM CDashTesting web page. If you do not like publish results pleas remove option ``-D Experimental``
