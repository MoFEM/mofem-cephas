## Installation with Spack (Recommended for Linux & Mac OS X)## {#install_spack}

All that you need to know about [Spack](https://spack.io) and more you find
here [https://spack.io](https://spack.io). For short Spack is a package
manager for supercomputers, Linux, and Mac OS X. It is designed to make
installation of scientific packages as easy as possible.

## Quick installation snippet

~~~~~~
git clone -single-branch -b develop https://github.com/likask/spack.git
. spack/share/spack/setup-env.sh
spack install mofem-users-modules
spack view --verbose symlink  um_view mofem-cephas
export PATH=$HOME/um_view/bin:$PATH
~~~~~~

Installation can take some time, with MoFEM all dependent libraries are
installed like PETSc, MoAB, Mumps, SuperLU, Parmetis and many others. You
have to be patient.

## Installation Spack

Clone Spack from GitHub and you’re ready to go:
~~~~~~
git clone -single-branch -b develop https://github.com/likask/spack.git
. spack/share/spack/setup-env.sh
~~~~~~
Note that we used forked Spack repository on GitHub. Forked Spack repository
is under control of MoFEM developers. Moreover is cloned develop branch,
where we put most up to date changes in the core library and users modules.
Over time we will make pull-requests to official Spack repository and make
MoFEM installation with main Spack repository.

## Installation MoFEM

You can install MoFEM in two different ways;
- if you are going to develop or use users modules, i.e. use or write an application using MoFEM library. 
- if you are going to develop core MoFEM library 

If you do not know at that point what to do, you probably like to choose the first option.

Spack has excellent documentation, reading it is not essential to make
successful MoFEM installation, however sooner or later you will have to look
into it. For a start read
- [Basic usage](http://spack.readthedocs.io/en/latest/basic_usage.html)
- [Filesystem Views](http://spack.readthedocs.io/en/latest/workflows.html#filesystem-views)

### Quick MoFEM installation

##### Installation
Basic installation of users modules is as short us
~~~~~~
spack install mofem-users-modules
~~~~~~
Spack will install all dependencies including core MoFEM library, PETSc and
MoAB, i.e. MoFEM software ecosystem. Once basic users modules are installed
you can start to work. For example, you can link your home directory to your
installation.

##### Users modules

You can but not have to add some users modules for example fracture module
 ~~~~~~
spack install mofem-fracture-module
 ~~~~~~
or minimal surface equation tutorial module
~~~~~~
spack install mofem-minimal-surface-equation
~~~~~~

You can list all modules in spack 
~~~~~~
spack list mofem
~~~~~~

Not all modules are yet added to spack; if your module is not yet there, you
can install it by hand, that is by cloning appropriate users module
repository and complaining code after spack package view is created.

Each spack package can be installed with variants and with a specific
version, to see information about the package
~~~~~~
spack info mofem-fracture-module
~~~~~~
and if you like to install mofem-fracture-module for version  v0.9.38, you do
~~~~~~
spack info mofem-fracture-module@0.9.38
~~~~~~
or from a development branch
~~~~~~
spack info mofem-fracture-module@develop
~~~~~~

##### Package view

A filesystem view is a single directory tree that is the union of the
directory hierarchies of a number of installed packages; it is similar to the
directory hierarchy that might exist under */usr/local*. The files of the
view’s installed MoFEM packages are brought into the view by symbolic or hard
links, referencing the original Spack installation.
~~~~~~
cd $HOME
spack view --verbose symlink  um_view mofem-cephas
~~~~~~

Note that mofem package view has a directory hierarchy including *bin*,
*etc*, *lib* and other system directories. At the time you work with a
particular view (you can have more than one with different MoFEM versions)
you can add view *bin* directory to the shell search path
~~~~~~
export PATH=$HOME/um_view/bin
~~~~~~

To check if all is working you can run tests and submit results to CDash
~~~~~~
cd $HOME/um_view
make 
./bin/ctest -D Experimental
~~~~~~

Now you can start to work, use code from users modules, or develop your own
user module code, for example run elastic analys of L-shape beam
~~~~~~
cd $HOME/um_view/basic_finite_elements/elasticity
../../bin/mpirun -np 2 ./elasticity -my_file LShape.h5m -ksp_type gmres -pc_type lu -pc_factor_mat_solver_package mumps -ksp_monitor -my_order 2
../../bin/mbconvert out.h5m out.vtk
~~~~~~
and finally open VTK file in [ParaView](https://www.paraview.org). You can
install ParaView using Spack or use install binary for your native OS.

### Fine grain installation

You can have more control on version and core library and other dependent
mofem libraries, i.e. do the fine installation

~~~~~~
cd $HOME
spack install --verbose --test root mofem-cephas@0.8.1+slepc
spack install mofem-users-modules ^mofem-cephas@0.8.1+slepc
spack install mofem-minimal-surface-equation@develop ^mofem-cephas@0.8.1+slepc
spack install mofem-fracture-module@0.9.38 ^mofem-cephas@0.8.1+slepc
spack view --verbose symlink  um_view_0.8.1 mofem-users-modules^mofem-cephas@0.8.1+slepc
~~~~~~

With the above snippet, we install a specific version of core mofem library,
i.e. 0.8.1, then we extend mofem by two users modules. Also, MoFEM is
installed with the variant SLEPC, i.e. PETSc extension enabling parallel
calculation of eigenproblems. Module minimal-surface-equation is installed
form developemnt branch and fracture module is installed with the version
0.9.38. At the end we male link to mofem core library with dependencies like
PETSc and Moab and extensions like users modules.

Also, note that we install MoFEM core library, displaying verbose build
output while installing and running tests verifying the correctness of
installation.

To uninstall mofem-cephas@0.8.1 with all dependents libraries, you run spack command
~~~~~
spack uninstall --dependents  mofem-cephas@0.8.1+slepc
~~~~~

### Core lib developer installation

If you are going to develop MoFEM core library, it means that you are core
developer; you can install mofem directly from the source.

Create mofem_install folder in the home directory and clone MoFEMrespository
~~~~~
mkdir $HOME/mofem_install
cd $HOME/mofem_install
git clone -b develop https://bitbucket.org/likask/mofem-cephas.git mofem-cephas
cd $HOME/mofem_install/mofem-cephas
git submodule update --init mofem/users_modules
~~~~~

Create *build* directory for core liblary
~~~~~
mkdir $HOME/mofem_install/lib
cd $HOME/mofem_install/lib
~~~~~

Configure installation using Spack,
~~~~~
spack setup mofem-cephas@develop build_type=Debug copy_user_modules=False
~~~~~
The result will be a file spconfig.py in the top-level
*$HOME/mofem_install/lib* directory. It is a short script that calls CMake
with the dependencies and options determined by Spack — similar to what
happens in spack install but now written out in script form. From a
developer’s point of view, you can think of spconfig.py as a stand-in for the
cmake command.

Note variant *copy_user_modules=False* is set so users modules are not copied
to install directory by indicating that symbolic is created. That is useful
when you do changes both in the core library and basic users modules.

Run config file
~~~~~
./spconfig.py ../mofem-cephas/mofem/
~~~~~

Create a view of installed CMake package in Spack
~~~~~
spack view --verbose symlink cmake_view cmake
~~~~~
and now we can compile, and run tests and finally install code
~~~~~
make -j4 
./cmake_view/bin/ctest -D Experimental
make install
~~~~~

Once that is done a view to users modules is created
~~~~~
cd $HOME/mofem_install
spack view --verbose symlink  um_view mofem-cephas@develop build_type=Debug copy_user_modules=False
~~~~~

Before we compile users modules you fix the link to users modules source
code, which in this case is a soft link to a directory in core lib repository
~~~~~
cd $HOME/mofem_install/um_view
rmdir users_modules
ln -s ../mofem-cephas/mofem/users_modules
~~~~~
One can that is done you can compile users modules and run tests
~~~~~
./bin/cmake -DCMAKE_BUILD_TYPE=Debug -DSTAND_ALLONE_USERS_MODULES=0 users_modules
make -j4
ctest -D Experimental
~~~~~

### Adding user module

Look at [Spack creation tutorial](https://spack.readthedocs.io/en/latest/tutorial_packaging.html#packaging-tutorial)
and [Spack Package Build Systems](https://spack.readthedocs.io/en/latest/tutorial_buildsystems.html). 
You can look how we have created packages for *fracture module* or *minimal-surface-equation*,
located in *$HOME/spack/var/spack/repos/builtin/packages/mofem-fracture-module*
 and *$HOME/spack/var/spack/repos/builtin/packages/mofem-minimal-surface-equation*.
You can open package file
~~~~~
spack edit mofem-minimal-surface-equation
~~~~~

To create your package for the user module, you have to 
- create directory following naming convention
- copy existing mofem user module package 
- create package class following naming convention
- modify file changing names and locations appropriately for your user module


