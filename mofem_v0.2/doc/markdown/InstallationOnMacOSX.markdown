## Installation on MacOSX ##

If you have any problems, feedback or would like to suggest corrections,
please email [cmatgu@googlegroups.com](mailto:cmatgu@googlegroups.com).

###1. Open a terminal

Create mofem install director:
~~~~~~
export MOFEM_INSTALL_DIR=$HOME/mofem_installation
mkdir $MOFEM_INSTALL_DIR
~~~~~~

###2. Install libraries using homebrew

~~~~~~
ruby -e "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install)"
brew update
brew tap homebrew/science
brew install cmake
brew install --with-x11 gnuplot
brew install --with-graphviz doxygen
brew install wget
~~~~~~

###3. Install PETSc and other liblaries

Change to your MOFEM_INSTALL_DIR
~~~~~~
cd $HOME/$MOFEM_INSTALL_DIR
~~~~~~

Colone PETSc repository:
~~~~~~
export PETSC_VERSION 3.5.3
git clone https://bitbucket.org/petsc/petsc.git
cd /opt/petsc
git checkout tags/v$PETSC_VERSION
~~~~~~

Note that version pets3.5.3 is checkout. MoFEM should work with development version of petsc if mofem is configured with `` -DCMAKE_C_FLAGS="-Wall -DPETSC_DEV"  -DCMAKE_CXX_FLAGS="-Wall -DPETSC_DEV"``.

Configure and compile petsc:
~~~~~~
wget ftp://ftp.unidata.ucar.edu/pub/netcdf/netcdf-4.3.3.1.tar.gz
./configure --with-debugging=0 --download-fblaslapack=1 --download-superlu_dist=1 --download-metis=1 --download-parmetis=1 --download-hypre=1 --download-mumps=1 --download-scalapack=1 --download-blacs=1 --download-moab=1 --download-hdf5=1 --download-netcdf=/opt/petsc/netcdf-4.3.3.1.tar.gz  --download-openmpi=1 --download-mumps=1 --download-ptscotch=1
make PETSC_DIR=$PWD PETSC_ARCH=arch-darwin-c-opt all
~~~~~~

Note that PETSc is compiled with debugging switch off for efficiency. If you
develop code is recommened that you compile petsce with debugging flag on.

###4. Compile ADOL-C

Change to your MOFEM_INSTALL_DIR
~~~~~~
cd $HOME/$MOFEM_INSTALL_DIR
~~~~~~

Get linraryu, configure make and install
~~~~~~
wget http://www.coin-or.org/download/source/ADOL-C/ADOL-C-2.5.2.zip
tar -xzf ADOL-C-2.5.2.zip
cd ADOL-C-2.5.2
./configure --prefix=$MOFEM_INSTALL_DIR/local
make install
~~~~~~

###5. Install TetGem and other libraries

####5.1 TetGen
~~~~~~
cd $HOME/$MOFEM_INSTALL_DIR
wget https://bitbucket.org/likask/mofem-joseph/downloads/tetgen1.5.0.tgz
tar -xvvzf tetgen1.5.0.tgz
cd tetgen1.5.0
cmake .
make
cp libtet.a lib/
~~~~~~

###6. Clone sourcecode and install MoFEM liblary

Change to your MOFEM_INSTALL_DIR
~~~~~~
cd $HOME/$MOFEM_INSTALL_DIR
~~~~~~

Cloning MoFEM sourcecode:
~~~~~~
git clone https://bitbucket.org/likask/mofem-cephas.git mofem-cephas
~~~~~~

Make a build directory
~~~~~~
mkdir $MOFEM_INSTALL_DIR/lib
cd $MOFEM_INSTALL_DIR/lib
~~~~~~


Configuring and compiling code:
~~~~~~
cmake -DCMAKE_Fortran_COMPILER=/usr/local/bin/gfortran -DCMAKE_BUILD_TYPE=Release -DCMAKE_C_FLAGS="-Wall"  -DCMAKE_CXX_FLAGS="-Wall  -Wno-bind-to-temporary-copy -DPETSC_DEV" -DPETSC_DIR=$MOFEM_INSTALL_DIR/petsc/ -DPETSC_ARCH=arch-darwin-c-opt -DMOAB_DIR=$MOFEM_INSTALL_DIR/petsc/arch-darwin-c-debug/ -DADOL-C_DIR=$MOFEM_INSTALL_DIR/local/ -DCMAKE_INSTALL_PREFIX=$MOFEM_INSTALL_DIR/user_modules $MOFEM_INSTALL_DIR/mofem-cephas/mofem_v0.2
~~~~~~

Building code (assuming that you have computer with 4 cores):
~~~~~~
make -j4 install
~~~~~~

Testing and publishing results on MoFEM CDashTesting WebPage:
~~~~~~
ctest -D Experimental
~~~~~~

###7. Configuration, compilation and testing user modules

Change to your MOFEM_INSTALL_DIR
~~~~~~
cd $HOME/$MOFEM_INSTALL_DIR/user_modules
~~~~~~

Configuration:
~~~~~~
cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_C_FLAGS="-Wall"  -DCMAKE_CXX_FLAGS="-Wall  -Wno-bind-to-temporary-copy -DPETSC_DEV" $MOFEM_INSTALL_DIR/user_modules
~~~~~~

Build:
~~~~~~
make -j4
~~~~~~

Testing:
~~~~~~
cmake -D Experimental
~~~~~~

Note that results of the test are publish on MoFEM CDashTesting web page. If you do not like publish results pleas remove option ``-D Experimental``
