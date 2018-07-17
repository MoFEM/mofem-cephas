## Installation with Spack (Recommended for Linux & Mac OS X)## {#install_spack}

All that you need to know about [Spack](https://spack.io) and more you find
here [https://spack.io](https://spack.io). For short Spack is a package
manager for supercomputers, Linux, and Mac OS X. It is designed to make
installation of scientific packages as easy as possible.

## Quick installation snippet

Befor you start see Spack [getting started](https://spack.readthedocs.io/en/v0.10.0/getting_started.html), to see all
prerequisites.
~~~~~~
git clone --single-branch -b mofem https://github.com/likask/spack.git
. spack/share/spack/setup-env.sh
. ${SPACK_ROOT}/share/spack/setup-env.sh
spack install mofem-users-modules
spack view --verbose symlink um_view mofem-cephas
spack activate -v um_view mofem-users-modules
export PATH=$HOME/um_view/bin:$PATH
~~~~~~

Installation can take some time, with MoFEM all dependent libraries are
installed like PETSc, MoAB, Mumps, SuperLU, Parmetis and many others. You
have to be patient.

## Installation Spack

Check if you have all [prerequisites](https://spack.readthedocs.io/en/v0.10.0/getting_started.html) and you are ready 
to start. Clone Spack from GitHub and you’re ready to go:
~~~~~~
git clone --single-branch -b mofem https://github.com/likask/spack.git
. spack/share/spack/setup-env.sh
. ${SPACK_ROOT}/share/spack/setup-env.sh
~~~~~~
Note that we used forked Spack repository on GitHub. Forked Spack repository
is under control of MoFEM developers, where we put most up to date changes in
the core library and users modules. 

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

##### Installation basic users modules

Basic installation of users modules is as short us
~~~~~~
spack install mofem-users-modules
spack view --verbose symlink  um_view mofem-cephas
spack activate -v um_view mofem-users-modules
~~~~~~
Spack will install all dependencies including core MoFEM library, PETSc and
MoAB, i.e. MoFEM software ecosystem. 

##### Adding more users modules

You can but not have to add some users modules for example fracture module
~~~~~~
spack install mofem-users-module+mofem-fracture-module
cd $HOME
spack view --verbose symlink um_view_foo mofem-cephas
spack activate -v um_view_foo mofem-users-module+mofem-fracture-module
 ~~~~~~
or minimal surface equation tutorial module
~~~~~~
spack deactivate -v um_view_foo mofem-users-module+mofem-fracture-module
spack install mofem-users-module+mofem-minimal-surface-equation
spack activate -v um_view_foo mofem-users-module+mofem-minimal-surface-equation
~~~~~~
You can as well install install multiple users modules, for example
~~~~~~
spack install mofem-users-module+mofem-minimal-surface-equation+mofem-fracture-module
cd $HOME
spack view --verbose symlink um_view_combined mofem-cephas
spack activate -v um_view_combined mofem-users-module+mofem-minimal-surface-equationmofem-fracture-module
~~~~~~
You can have only one activated * mofem-users-module* in view if you like to create a different one with a specific set of users modules, you need to deactivate the old one. Alternatively, you can create another view and in it activate desired user modules.

You can list all MoFEM packages in Spack 
~~~~~~
spack list mofem
~~~~~~
or extensions 
~~~~~~
spack extensions mofem-cephas
~~~~~~

Not all modules are yet added to spack; if your module is not yet there, you
can install it by hand, that is by cloning appropriate users module
repository and complaining code after spack package view is created.

##### Package view

A filesystem view is a single directory tree that is the union of the
directory hierarchies of a number of installed packages; it is similar to the
directory hierarchy that might exist under */usr/local*. The files of the
view’s installed MoFEM packages are broughtinto the view by symbolic or hard
links, referencing the original Spack installation.
~~~~~~
cd $HOME
spack view --verbose symlink um_view mofem-cephas
spack activate -v um_view mofem-users-modules
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
ctest -D Experimental
~~~~~~

Now you can start to work, use code from users modules, or develop your own
user module code, for example run elastic analys of L-shape beam
~~~~~~
cd $HOME/um_view/build/basic_finite_elements/elasticity
mpirun -np 2 ./elasticity -my_file LShape.h5m -ksp_type gmres -pc_type lu -pc_factor_mat_solver_package mumps -ksp_monitor -my_order 2
mbconvert out.h5m out.vtk
~~~~~~
and finally open VTK file in [ParaView](https://www.paraview.org). You can
install ParaView using Spack or use install binary for your native OS.

### Core lib developer installation

If you are going to develop MoFEM core library, it means that you are core
developer; you can install mofem directly from the source.

Create mofem_install folder in the home directory and clone MoFEMrespository
~~~~~
mkdir $HOME/mofem_install
cd $HOME/mofem_install
git clone -b develop --recurse-submodules https://bitbucket.org/likask/mofem-cephas.git mofem-cephas
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

We load cmake tools
~~~~~
spack load cmake
~~~~~
and run config file
~~~~~
./spconfig.py ../mofem-cephas/mofem/
~~~~~
next make code and run test
~~~~~
make -j4 
ctest -D Experimental
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
spack edit mofem-users-modules
spack edit mofem-minimal-surface-equation
~~~~~

To create your package for the user module, you have to 
- create directory following naming convention
- copy existing mofem user module package 
- create package class following naming convention
- modify file changing names and locations appropriately for your user module
- edit *mofem-users-modules* to add variant for created extension 


