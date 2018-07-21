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
spack bootstrap
~~~~~~
Having spack installed you can install a basic version of MoFEM
~~~~~~
spack install mofem-users-modules
spack view --verbose symlink um_view mofem-cephas
spack activate -v um_view mofem-users-modules
export PATH=$HOME/um_view/bin:$PATH
~~~~~~
Installation can take some time, with MoFEM all dependent libraries are
installed like PETSc, MoAB, Mumps, SuperLU, Parmetis and many others. You
have to be patient.

MoFEM is extendable by users modules, called in spack extensions. Available in spack extension can be seen by calling
~~~~~~
spack extensions mofem-cephas
~~~~~~

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

You can install all other prerequisites needed by spack by calling command
~~~~~~
spack bootstrap
~~~~~~

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

##### Installation basic users modules

Basic installation of users modules is as short us
~~~~~~
spack install mofem-users-modules
spack view --verbose symlink  um_view mofem-cephas
spack activate -v um_view mofem-users-modules
~~~~~~
Spack will install all dependencies including core MoFEM library, PETSc and
MoAB, i.e. MoFEM software ecosystem. 

Note that MoFEM package view has a directory hierarchy including *bin*,
*etc*, *lib* and other system directories. At the time you work with a
particular view (you can have more than one with different MoFEM versions)
you can add view *bin* directory to the shell search path
~~~~~~
export PATH=$HOME/um_view/bin
~~~~~~

Now you can start to work, use code from users modules, or develop your own
user module code, for example run elastic analysis of L-shape beam
~~~~~~
cd $HOME/um_view/elasticity
mpirun -np 2 ./elasticity \
-my_file LShape.h5m \
-my_order 2 \
-ksp_type gmres \
-pc_type lu -pc_factor_mat_solver_package mumps \
-ksp_monitor
mbconvert out.h5m out.vtk
~~~~~~
and finally open VTK file in [ParaView](https://www.paraview.org). You can
install ParaView using Spack or use install binary for your native OS.


##### Adding more users modules

You can but not have to add some users modules for example fracture module
~~~~~~
spack install mofem-fracture-module
cd $HOME
spack view --verbose symlink um_view_foo mofem-cephas
spack activate -v um_view_foo mofem-fracture-module
 ~~~~~~
or minimal surface equation tutorial module
~~~~~~
spack install mofem-minimal-surface-equation
spack activate -v um_view_foo mofem-minimal-surface-equation
~~~~~~

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
view’s installed MoFEM packages are brought in to the view by symbolic or hard
links, referencing the original Spack installation.
~~~~~~
cd $HOME
spack view --verbose symlink um_view mofem-cephas
spack activate -v um_view mofem-users-modules
~~~~~~

## For developers

You can develop code in MoFEM on different levels, work with development of user module, i.e. solving a particular problem using finite elements. Or contribute to necessary modules which are used for the majority of MoFEM users, or finally be a core MoFEM developer. Depending on what you are going to do you have several options to choose how to set up the environment for development with spack. Below we list some of them.

### Making code in view

This is an option for one who likes to play with the code, make changes and modifications and check how they work.  You can work with the code provided in modules or have add own directory which your module.
~~~~
spack install mofem-users-modules
cd $HOME
spack view --verbose symlink um_view mofem-cephas
cd um_view
mkdir build
cmake \
-DCMAKE_BUILD_TYPE=Debug \
-DSTAND_ALLONE_USERS_MODULES=YES \
-DEXTERNAL_MODULE_SOURCE_DIRS=../ext_users_modules \
../users_modules
~~~~

You can provide another directory which is pointing to your module which you (secretly) developing
~~~~
-DEXTERNAL_MODULE_SOURCE_DIRS=../ext_users_modules\;$PATH_TO_MY_SECRET_MODULE
~~~~

### Module developer

To develop module, you need install mofem-users-modules, clone from
repository which you like to work with, set up configuration and make the
code.
~~~~
spack install mofem-users-modules
cd $HOME
spack view --verbose symlink um_view mofem-cephas
mkdir $HOME/mod_developer
cd mod_developer/
git clone -b develop https://bitbucket.org/likask/mofem_um_minimal_surface_equation minimal_surface_equation
spack setup mofem-minimal-surface-equation@develop
./spconfig.py -DMOFEM_DIR=$HOME/um_view \
-DEXTERNAL_MODULE_SOURCE_DIRS=$HOME/mod_developer \
$HOME/um_view/users_modules
make -j4
~~~~
Note command 
~~~~
spack setup mofem-minimal-surface-equation@develop
~~~~
It creates spconfig.py file, which is python script which calls CMake with
appropriate paths. You run it to set configuration paths to other libraries.
For general notes on workflows how to build packages with Spack look
[here](https://spack.readthedocs.io/en/latest/workflows.html?highlight=spconfig.py#build-with-spack).

### Users module developer

To develop users modules procedure is similar to one shown above, one needs
to install mofem-cephas, clone source code, run configuration and finally
make code
~~~~
spack install mofem-cephas
mkdir $HOME/um_developer
cd $HOME/um_developer/
git clone -b develop https://likask@bitbucket.org/mofem/users-modules-cephas.git 
spack setup mofem-users-modules@develop
spack view --verbose symlink um_view mofem-cephas@0.8.7
./spconfig.py -DMOFEM_DIR=$HOME/um_view users-modules-cephas/ 
make -j4
~~~~

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
CMake command.

Note variant *copy_user_modules=False* is set so users modules are not copied
to install directory by indicating that symbolic is created. That is useful
when you do changes both in the core library and basic users modules. 

We load CMake tools
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
make install
~~~~~

Once that is done a view to users modules is created
~~~~~
cd $HOME/mofem_install
spack view --verbose symlink um_view mofem-cephas@develop build_type=Debug copy_user_modules=False
~~~~~

One can that is done you can compile users modules and run tests
~~~~~
cd um_view
mkdir build
cd build
cmake -DMOFEM_DIR=$HOME/mofem_install/um_view/ $HOME/mofem_install/mofem-cephas/mofem/users_modules
make -j4
~~~~~

## Adding MoFEM extension to Spack, i.e. user module

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

## Mirrors

Some sites may not have access to the internet for fetching packages. These sites will need a local repository of tarballs from which they can get their files. Spack has support for this with mirrors. Look to Spack documentation to learn more about mirrors, see [here](https://spack.readthedocs.io/en/latest/mirrors.html?highlight=mirror).

You need to create mirror first
~~~~~
spack mirror create -D mofem-fracture-module
~~~~~
as a result, a directory is created with are prerequisites need to Spack
installation of mofem-fracture-module. You can now move that directory to 
your secure server with spack package.

On a secure server, you add mirror directory to Spack, for example
~~~~~
spack mirror add local_filesystem file://$HOME/spack-mirror-2018-07-21
~~~~~
and with that at hand kick-start installation process described above.
