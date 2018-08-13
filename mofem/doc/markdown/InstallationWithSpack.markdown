Installation with Spack (Recommended for HPC, Linux & Mac OS X) {#install_spack}
==========================================================

All that you need to know about [Spack](https://spack.io) and more you find
here [https://spack.io](https://spack.io). For short Spack is a package
manager for supercomputers, Linux, and Mac OS X. It is designed to make
installation of scientific packages as easy as possible.

If you have any problems, feedback or would like to suggest corrections,
please email
[mofem-group@googlegroups.com](https://groups.google.com/forum/#!forum/mofem-group)..

[TOC]

# Quick installation snippet {#spack_quick}


If you have [git](https://www.atlassian.com/git/tutorials/what-is-git) and
[curl](https://en.wikipedia.org/wiki/CURL) and compilers you can proceed to
Spack installation. If you already develop some code, most likely you have
all that you need. If not look to for install
[prerequisites](@ref spack_prerequisites).

You can start by downloading Spack and installing essential, as follows
~~~~~~
git clone --single-branch -b mofem https://github.com/likask/spack.git
. spack/share/spack/setup-env.sh
spack bootstrap
. ${SPACK_ROOT}/share/spack/setup-env.sh
~~~~~~

Having spack installed you can install a basic version of MoFEM
~~~~~~
spack install mofem-users-modules
spack view --verbose symlink -i um_view mofem-cephas
spack activate -v um_view mofem-users-modules
export PATH=$HOME/um_view/bin:$PATH
~~~~~~
Installation can take some time, with MoFEM all dependent libraries are
installed like PETSc, MoAB, Mumps, SuperLU, Parmetis and many others. You
have to be patient.

MoFEM is extendable by users modules, called in spack extensions. Available
in spack extension can be seen by calling
~~~~~~
spack extensions mofem-cephas
~~~~~~

# Prerequisites {#spack_prerequisites}

Before you start see Spack [getting
started](https://spack.readthedocs.io/en/v0.10.0/getting_started.html), to
see all prerequisites. You need to have installed
[git](https://www.atlassian.com/git/tutorials/what-is-git) and
[curl](https://en.wikipedia.org/wiki/CURL). You will need as well C++
compilers, f.e. *gcc* or *clang*, and Fortran compiler, f.e. *gfrortran*. If
you are developing any code on your computer, you probably have probably most
of the software (if not all) installed.

For example on Ubuntu or Debian
like system you can install both running from command line
~~~~~
apt-get update \
&& apt-get install -y --no-install-recommends \
autoconf \
build-essential \
ca-certificates \
coreutils \
curl \
environment-modules \
git \
python \
unzip \
vim \ 
gfortran
~~~~~

If you work on OS X, and if have not done it already, install
build-essentials and [Homebrew](https://brew.sh) and then
~~~~~
xcode-select --install
sudo xcodebuild -license accept
/usr/bin/ruby -e "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install)"
brew install curl git 
~~~~~
For Mac OS X you can as well follow instrucion on Spack documentation, called
[Mixed
Toolchains](http://spack.readthedocs.io/en/latest/getting_started.html#mixed-toolchains)

Check if you have Fortran compiler, f.e. *gfortran* if not you can install
with homebrew or Ubuntu/Debian packages. 
~~~~~
apt-get install -y gfortran
~~~~~
or
~~~~~
brew install gfortran
~~~~~
respectively for Ubuntu and Mac OS X.

Alternatively, If your system does not provide any Fortran compiler or you
want to have the most recent gcc. Just before kick-starting MoFEM
installation command with Spack, do as follows
~~~~~
spack install gcc
spack load gcc
spack compiler find
spack install mofem-users-modules
~~~~~

Installing *gcc* with Spack is more resistant to future problems. For
example, updating codes on hombrew old *gfortran* could be updated, causing
problems with linking of installed packages. Adding *gcc* compiler with Spack
resistant to those type of problems.

# Installation Spack {#spack_spack}

Check if you have all [prerequisites](https://spack.readthedocs.io/en/v0.10.0/getting_started.html) and you are ready 
to start. Clone Spack from GitHub and you’re ready to go:
~~~~~~
git clone --single-branch -b mofem https://github.com/likask/spack.git
. spack/share/spack/setup-env.sh
~~~~~~
Note that we used forked Spack repository on GitHub. Forked Spack repository
is under control of MoFEM developers, where we put most up to date changes in
the core library and users modules. 

You can install all other prerequisites needed by spack by calling command
~~~~~~
spack bootstrap
. ${SPACK_ROOT}/share/spack/setup-env.sh
~~~~~~

# Installation of MoFEM {#spack_mofem}

You can install MoFEM in two different ways;
- if you are going to develop or use users modules, i.e. use or write an application using MoFEM library. 
- if you are going to develop core MoFEM library 

If you do not know at that point what to do, you probably like to choose the first option.

Spack has excellent documentation, reading it is not essential to make
successful MoFEM installation, however sooner or later you will have to look
into it. For a start read
- [Basic usage](http://spack.readthedocs.io/en/latest/basic_usage.html)
- [Filesystem Views](http://spack.readthedocs.io/en/latest/workflows.html#filesystem-views)

## Installation of basic users modules {#spack_basic}

Basic installation of users modules is as short us
~~~~~~
spack install mofem-users-modules
spack view --verbose symlink -i um_view mofem-cephas
spack activate -v um_view mofem-users-modules
~~~~~~
Spack will install all dependencies including core MoFEM library, PETSc and
MoAB, i.e. MoFEM software ecosystem. 

Note that MoFEM package view has a directory hierarchy including *bin*,
*etc*, *lib* and other system directories. At the time you work with a
particular view (you can have more than one with different MoFEM versions)
you can add view *bin* directory to the shell search path
~~~~~~
export PATH=$HOME/um_view/bin:$PATH
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

## Adding more users modules {#spack_add_more}

You can but not have to add some users modules for example fracture module
~~~~~~
spack install mofem-fracture-module
cd $HOME
spack view --verbose symlink -i um_view_foo mofem-cephas
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

## Package view {#spack_view}

A filesystem view is a single directory tree that is the union of the
directory hierarchies of a number of installed packages; it is similar to the
directory hierarchy that might exist under */usr/local*. The files of the
view’s installed MoFEM packages are brought in to the view by symbolic or hard
links, referencing the original Spack installation.
~~~~~~
cd $HOME
spack view --verbose symlink -i um_view mofem-cephas
spack activate -v um_view mofem-users-modules
~~~~~~

## Package install and running tests {#spack_running_tests}

You can install package and run tests, for example
~~~~~~
spack install --test=root -j 4 mofem-cephas
spack install --test=root -j 4 mofem-users-modules
spack install --test=root -j 4 mofem-fracture-module
~~~~~~

## Setting build type and compiler flags

During spack installation, you can set build type, as follows
~~~~~~
spack install mofem-cephas build_type="Debug"
~~~~~~
or compiler flaks, as follows
~~~~~~
spack install mofem-cephas cppflags="-march=native -O3"
~~~~~~

# Installation on specific servers {#spack_servers}

## RDB development server {#spack_rdb}

On RDB server we are running Ubuntu 12.04.5 LTS, which is very old and does
not have compilers which can compile the C++14 code required by MoFEM and
some dependent packages. Difficulties with the compiler are not only one. The
[curl](https://en.wikipedia.org/wiki/CURL) version on Ubuntu 12.04.5 LTS is
compiled with old
[SSL](https://en.wikipedia.org/wiki/Transport_Layer_Security) protocol and
fetching files from some servers does work. In order to solve those problems,
we will use a local mirror, which previously downloaded packages and install
a new version of GCC.

First step is to get Spack
~~~~~
git clone -b mofem https://github.com/likask/spack.git
. spack/share/spack/setup-env.sh
~~~~~
Next, we add local mirror to Spack,
~~~~~
spack mirror add local_filesystem file://opt/spack-mirror-2018-07-26 
~~~~~
and compile all prerequisites
~~~~~
spack bootstrap
. ${SPACK_ROOT}/share/spack/setup-env.sh
~~~~~

Now we install GCC package, load and add compiler, install newest version of
curl, if we need fetch something from external source, and finally build
MoFEM with all dependencies
~~~~~
spack install -j 8 gcc
spack load gcc
spack compiler find
spack install -j 8 curl
spack load curl
spack install -j 8 mofem-users-modules%gcc@8.1.0
~~~~~
All this take some time, so you can run this on
[screen](https://www.youtube.com/watch?v=3txYaF_IVZQ) terminal. You can
easily go for a walk to Kelvingrove Park and come back.

Now you can create a symlink to install directory including dependent
libraries, using commands below
~~~~
spack view symlink -i um_view mofem-cephas
spack activate -v um_view mofem-users-modules
~~~~

## Server Buckethead {#spack_buckedhead}

### Installation {#spack_buckedhead_installation}

Buckethead is a Linux cluster running Centos7. More information about
Buckethead you will find 
[here](http://malmsteen.eng.gla.ac.uk/wiki/hpc-wiki/index.php/Resources#Buckethead). 
Install Spack and load cluster modules
~~~~~
module load gcc/6.4.0
module load gridengine
git clone --single-branch -b mofem https://github.com/likask/spack.git
. spack/share/spack/setup-env.sh
~~~~~

We need to set up compilers since the location of standard gcc@6.4.0 libraries is in an unusual place when loaded by the module. In order to do that you have to edit file *.spack/linux/compilers.yaml* in your home directory
~~~~~
     1    compilers:
     2    - compiler:
     3        environment: {}
     4        extra_rpaths: []
     5        flags: {}
     6        modules: []
     7        operating_system: centos7
     8        paths:
     9          cc: /usr/bin/gcc
    10          cxx: /usr/bin/g++
    11          f77: null
    12          fc: null
    13        spec: gcc@4.8.5
    14        target: x86_64
    15    - compiler:
    16        environment: {}
    17        extra_rpaths: []
    18        flags: {}
    19        modules: []
    20        operating_system: centos7
    21        paths:
    22          cc: /software/compilers/gcc/6.4.0/bin/gcc
    23          cxx: /software/compilers/gcc/6.4.0/bin/g++
    24          f77: /software/compilers/gcc/6.4.0/bin/gfortran
    25          fc: /software/compilers/gcc/6.4.0/bin/gfortran
    26        spec: gcc@6.4.0
    27        target: x86_64
~~~~~
and substitute line 17 by 
~~~~~
    17    extra_rpaths:
          - /software/compilers/gcc/6.4.0/lib64
~~~~~

At that point, we can follow the standard installation procedure, as follows
~~~~~
spack bootstrap
spack install mofem-users-modules
spack view --verbose symlink -i um_view mofem-cephas
spack activate -v um_view mofem-users-modules
~~~~~
If needs you can add more users modules or compile them by yourself.
Installation can take some time since everything is installed from scratch.
You can consider to run it in
[screen](https://www.youtube.com/watch?v=3txYaF_IVZQ) terminal, and go for a
coffee. 

Now you can create a symlink to install directory including dependent
libraries, using commands below
~~~~
spack view symlink -i um_view mofem-cephas
spack activate -v um_view mofem-users-modules
~~~~

### Job file {#spack_buckedhead_job}

Create a script file with content as below and name it, for example, *job_spack*
~~~~~
#! /bin/bash

# The job's name
#$ -N MoFEM

# The queue in which to run the job 
#$ -q gcec.q

# File to which standard error should be directed
#$ -e ./stderr

# File to which standard output should be directed
#$ -o ./stdout

# E-mail address to which status updates should be sent
# N.B.: in an array job, a separate e-mail will be sent for each task!
#$ -M your.email@uni.ac.uk

# Events on which to send a status update
#$ -m beas

# Request for 1.0 GB of memory per task (needed on Miffy and Dusty)
#$ -l mem_tokens=1.0G

#$ -pe mpi 2 # where N is the number of processors required

# List of commands which do the actual work
echo "$NSLOTS received"
cat $PE_HOSTFILE

# Load compiler
module load gcc/6.4.0

# List of commands which do the actual work
cd $HOME/um_view/elasticity
$HOME/um_view/bin/mpirun -np $NSLOTS \
./elasticity \
-my_file LShape.h5m \
-ksp_type gmres -pc_type lu \
-pc_factor_mat_solver_package mumps \
-ksp_monitor \
-my_order 2 2>&1 | tee log
~~~~~
and run it as follows
~~~~~
qsub job_spack
~~~~~
Results of the analysis are located in $HOME/um_view/elasticity. 

# For developers {#spack_developers}

You can develop code in MoFEM on different levels, work with development of user module, i.e. solving a particular problem using finite elements. Or contribute to necessary modules which are used for the majority of MoFEM users, or finally be a core MoFEM developer. Depending on what you are going to do you have several options to choose how to set up the environment for development with spack. Below we list some of them.

## Making code in view {#spack_make_in_view}

This is an option for one who likes to play with the code, make changes and modifications and check how they work.  You can work with the code provided in modules or have add own directory which your module.
~~~~
spack install mofem-users-modules
cd $HOME
spack view --verbose symlink -i um_view mofem-cephas
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

## Module developer {#spack_module_developer}

To develop module, you need install mofem-users-modules, clone from
repository which you like to work with, set up configuration and make the
code.
~~~~
spack install mofem-users-modules
cd $HOME
spack view --verbose symlink -i um_view mofem-cephas
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

## Users module developer {#spack_users_modules_developer}

To develop users modules procedure is similar to one shown above, one needs
to install mofem-cephas, clone source code, run configuration and finally
make code
~~~~
spack install mofem-cephas
mkdir $HOME/um_developer
cd $HOME/um_developer/
git clone -b develop https://likask@bitbucket.org/mofem/users-modules-cephas.git 
spack setup mofem-users-modules@develop
spack view --verbose symlink -i um_view mofem-cephas
./spconfig.py -DMOFEM_DIR=$HOME/um_view users-modules-cephas/ 
make -j4
~~~~

## Core lib developer {#spack_core_lib_developer}

If you are going to develop MoFEM core library, it means that you are core
developer; you can install mofem directly from the source.

Create mofem_install folder in the home directory and clone MoFEMrespository
~~~~~
mkdir $HOME/mofem_install
cd $HOME/mofem_install
git clone -b develop --recurse-submodules https://bitbucket.org/likask/mofem-cephas.git mofem-cephas
~~~~~
and kick-start installation of the core library
~~~~~
cd $HOME/mofem_install
mkdir lib
cd lib/
spack setup mofem-cephas@develop copy_user_modules=False
./spconfig.py $HOME/mofem_install/mofem-cephas/mofem/
make -j4
ctest
make install
~~~~~
Next, install users modules
~~~~~
cd $HOME/mofem_install
mkdir um
cd um/
spack view --verbose symlink -i um_view mofem-cephas@develop
mkdir build 
cd build/
spack setup mofem-users-modules@develop copy_user_modules=False ^mofem-cephas@develop
./spconfig.py -DMOFEM_DIR=../um_view $HOME/mofem_install/mofem-cephas/mofem/users_modules
make -j4
ctest
make install
~~~~~

You can add extended users modules to
*$HOME/mofem_install/mofem-cephas/mofem/users_modules*. To make the included
in build process you can re-do steps
~~~~~
./spconfig.py -DMOFEM_DIR=../um_view $HOME/mofem_install/mofem-cephas/mofem/users_modules
make -j4
ctest
make install
~~~~~

Alternatively, you can add your users modules to independent folder and run
snipped below
~~~~~
./spconfig.py \
-DEXTERNAL_MODULE_SOURCE_DIRS=$PATH_TO_MY_SECRET_MODULE \
-DMOFEM_DIR=../um_view \
$HOME/mofem_install/mofem-cephas/mofem/users_modules
make -j4
ctest
make install
~~~~~


# Adding MoFEM extension to Spack {#spack_adding_package}

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

# Basic settings in config.yaml {#spack_config}

Spack’s basic configuration options are set in config.yaml. You can see the
default settings by *$HOME/.spack/config.yaml*, for example, you can set the
default number of jobs to use when running `make` in parallel. If set to 4,
for example, `spack install` will run `make -j4`. This can be done by adding
line 
~~~~~ 
build_jobs: 4 
~~~~~ 
For more details see
[here](https://spack.readthedocs.io/en/latest/config_yaml.html?highlight=-jobs#)

# Mirrors {#spack_mirrors}

Some sites may not have access to the internet for fetching packages. These
sites will need a local repository of tarballs from which they can get their
files. Spack has support for this with mirrors. Look to Spack documentation
to learn more about mirrors, see
[here](https://spack.readthedocs.io/en/latest/mirrors.html?highlight=mirror).

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

# Known issues {#spack_issues}

## Mac OS X and DYLD_LIBRARY_PATH {#spack_dyld_library_path}

On your Mac OS X you could have already installed some libraries, like HDF5.
Some packages installed with Spack can use those libraries if a path to them
is set in DYLD_LIBRARY_PATH. DYLD_LIBRARY_PATH is a path to dynamics
libraries.

Before you run the installation set, unset DYLD_LIBRARY_PATH, as follows
~~~~~
export DYLD_LIBRARY_PATH=
spack install mofem-users-modules
~~~~~

# FAQ {#spack_faq}

## How to get packages installed today?

Run command
~~~~~
spack find -p -L --start-date `date +%F`
~~~~~
and you will get 
~~~~~
-- darwin-elcapitan-x86_64 / clang@7.3.0-apple ------------------
ugi2gm7    mofem-cephas@develop     /spack/opt/spack/darwin-elcapitan-x86_64/clang-7.3.0-apple/mofem-cephas-develop-ugi2gm7lydyjwmiwneefhxq7ahvpnuzc
gu4uheh    mofem-users-modules@develop  /spack/opt/spack/darwin-elcapitan-x86_64/clang-7.3.0-apple/mofem-users-modules-develop-gu4uhehai3edtqow7kivuxgpw5lmvbxz
~~~~~

## How to check when and were packages were installed?

Run command
~~~~~
ls -lhd `spack find -lp  | grep mofem | awk '{print $3}'`
~~~~~

# Contact {#spack_contact}

Any problems with this installation, please contact us by [mofem-group@googlegroups.com](https://groups.google.com/forum/#!forum/mofem-group)
or on Slack [MoFEM Slack](https://mofem.slack.com/).