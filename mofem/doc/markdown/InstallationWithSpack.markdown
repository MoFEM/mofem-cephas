Installation with Spack (Recommended for Linux, macOS, HPC) {#install_spack}
==========================================================

Spack is "A flexible package manager that supports multiple versions,
configurations, platforms, and compilers." -
[https://spack.io](https://spack.io)

A community of Spack users and developers contribute to developing package
configurations for a range of platforms, including macOS and Linux. This creates
a consistent build setup for supported scientific packages.

MoFEM can be deployed and developed using Spack, which is the recommended way of
installing. If you have any problems, feedback or would like to suggest
corrections,
please email to
[mofem-group@googlegroups.com](https://groups.google.com/forum/#!forum/mofem-group).
It should be noted that Spack compiles packages from source, pre-compiled binaries are not
available.

[TOC]

# Installation using scripts {#installation_scripts}

The remaining of this article presents a step-by-step guidance with
comments/explanations on how to install MoFEM using Spack in different
platforms (Linux, macOS, HPC). However, those who prefer to automise the process of the
installation can use the following scripts for different purposes:

- For user only (just the binary files): [`install_mofem_user.sh`](scripts/install_mofem_user.sh)
- For developer (full source and binary files): [`install_mofem_developer.sh`](scripts/install_mofem_developer.sh)

After downloading the appropriate file of choice, one should copy the file to
the directory where MoFEM will be installed. Then simply run the file from the
terminal by executing the command, for example,
~~~~~~
./install_mofem_user.sh
~~~~~~ 
It is worth noting that running the scripts may require user password for sudo privileges.

# Prerequisites {#spack_prerequisites}

The installation of MoFEM requires
[git](https://www.atlassian.com/git/tutorials/what-is-git),
[curl](https://en.wikipedia.org/wiki/CURL) and C++ and Fortran compilers (GCC
and/or clang). They must also be available in your `PATH`.

System requirements: Minimum 8GB RAM. On Linux ensure this RAM is free.

##Ubuntu/Debian

Install the following packages:
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

## macOS

Xcode contains required compilers for macOS. The last known working version of
Xcode is 10.2 with clang 10.0.1; there may be issues with future/different
versions.

To install the latest Xcode for your release of macOS, follow the commands
below, or see Apple's Xcode
[download page](https://developer.apple.com/download/) (Apple login
required) for older versions.
~~~~~
xcode-select --install
sudo xcodebuild -license accept
~~~~~

Additional packages are required - install [homebrew](https://brew.sh) package
manager:
~~~~~
/usr/bin/ruby -e "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install)"
~~~~~

Install packages through `homebrew`:
~~~~~
brew install curl git gcc
~~~~~

Check the path of the Fortran compiler (shipped with `gcc`):
~~~~~
which gfortran
~~~~~

If it is not already in the `PATH`, you should add it there.

Note: recent releases of macOS stopped shipping a Fortran compiler and therefore
require [Mixed
Toolchains](http://spack.readthedocs.io/en/latest/getting_started.html#mixed-toolchains).
The installing of gfortran through homebrew is another way of solving this.

# Spack setup {#spack_setup}

Retrieve Spack for MoFEM:
~~~~~~
git clone --single-branch -b mofem https://github.com/likask/spack.git
~~~~~~

Initialise Spack's environment variables:
~~~~~~
. spack/share/spack/setup-env.sh
~~~~~~

Spack's environment variables will be lost when the terminal session is closed.
Consider adding the previous command to your `.bash_profile` or `.bashrc`, e.g.:
~~~~~~
echo ". $HOME/spack/share/spack/setup-env.sh" >> ~/.bash_profile
~~~~~~


Finally, install packages required by Spack:
~~~~~~
spack bootstrap
~~~~~~

Note that there are further instructions on [Spack usage and configuration](#spack_usage_config).

# User only {#spack_user}

The first time installing MoFEM takes considerable time (hours). Spack fetches
and compiles MoFEM's dependencies from source and this includes several large
packages. 

Those seeking to develop MoFEM see [developer's](#spack_developers)
instructions.

## Install basic users modules

MoFEM's basic users modules consist of tools for solving a range of common
problems, from elasticity to the Poisson equation. To install them run the
following command:

~~~~~~
spack install mofem-users-modules
~~~~~~

MoFEM can be extended by adding user modules. More modules are available as
extensions in Spack. This will show available extensions:
~~~~~~
spack extensions mofem-cephas
~~~~~~

The extensions can be installed using `spack install <extension>`.

To access the installed users modules, create an `um_view` directory. This
should be created in an appropriate directory:
~~~~~~
spack view --verbose symlink -i um_view mofem-cephas
spack activate -v um_view mofem-users-modules
~~~~~~

This filesystem view is a single directory tree that is the union of the
directory hierarchies of a number of installed packages; it is similar to the
directory hierarchy that might exist under `/usr/local`. The files of the view's
installed MoFEM packages are brought into the view by symbolic or hard links,
referencing the original Spack installation. Different 'views' can be created
depending on the MoFEM version you wish to access. To make the new 'view'
visible from any directory, add its `bin` directory to your `PATH`, e.g.: 

~~~~~~
export PATH=$PWD/um_view/bin:$PATH 
~~~~~~ 

Consider also adding this command to your `.bash_profile` or `.bashrc`, e.g.: 
~~~~~~ 
echo "export PATH=$PWD/um_view/bin:\$PATH" >> ~/.bash_profile
~~~~~~

## Test elasticity module

~~~~~~
cd um_view/elasticity

mpirun -np 2 ./elasticity \
-my_file LShape.h5m \
-my_order 2 \
-ksp_type gmres \
-pc_type lu -pc_factor_mat_solver_type mumps \
-ksp_monitor
mbconvert out.h5m out.vtk
~~~~~~

Open the output VTK file in [ParaView](https://www.paraview.org). You can
install ParaView using Spack, [directly from the
website](https://www.paraview.org/download/) or through your OS package manager.

## Install additional users modules

For example, the fracture module can be installed by:
~~~~~~
spack install mofem-fracture-module
cd $HOME
spack view --verbose symlink -i um_view_foo mofem-cephas
spack activate -v um_view_foo mofem-fracture-module
~~~~~~
or the minimal surface equation tutorial module:
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

Not all modules have been added to spack. If your module is not yet there, you
can install manually by cloning the appropriate users module.

## Running tests {#spack_running_tests}

You can install package and run tests, for example:
~~~~~~
spack install --test=root -j 4 mofem-cephas
spack install --test=root -j 4 mofem-users-modules
spack install --test=root -j 4 mofem-fracture-module
~~~~~~


# Developers {#spack_developers}

MoFEM can be developed in different ways:

1. [In-situ](#spack_modules_insitu)
2. [Specific module](#spack_specific_module)
3. [Basic users modules](#spack_basic_users_modules)
4. [Core libraries](#spack_core_libraries)

The developer installation requires knowing more about how Spack works. See [Spack usage and configuration](#spack_usage_config)
before proceeding with the installation. In particular, the instructions below
will use so-called *specs* to obtain the desired build configuration, see
[Spack manual page](https://spack.readthedocs.io/en/latest/basic_usage.html#specs-dependencies)
for more details. For example, the default Spack specifier for the build type
will be given explicitly here: `build_type=RelWithDebInfo`, however keep in mind
that two other build types can be specified in the same way:
`build_type=Release` or `build_type=Debug`, see also [Change the build_type](#spack_build_type).

<!-- The instructions below will use `mofem-cephas@develop` and
`mofem-users-modules@develop` to denote the name of an existing installation.
However, if multiple versions are installed then you will need to [specify the
version ID](#spack_mofem_package_versions) instead. You can also [change the build_type](#spack_build_type). -->

## 1. In-situ {#spack_modules_insitu}

This is the simplest method but also limited. You can work with the code
provided in modules or add your own directory with your module.
~~~~
spack install mofem-users-modules build_type=RelWithDebInfo
cd $HOME
spack view --verbose symlink -i um_view mofem-users-modules
export PATH=$PWD/um_view/bin:$PATH
cd um_view
mkdir build
cd build
cmake \
-DCMAKE_BUILD_TYPE=Debug \
-DSTAND_ALLONE_USERS_MODULES=YES \
-DEXTERNAL_MODULE_SOURCE_DIRS=../ext_users_modules \
../users_modules
make -j4
ctest
~~~~

You can also provide another directory which is pointing to a module that you are
developing privately:
~~~~
-DEXTERNAL_MODULE_SOURCE_DIRS=../ext_users_modules\;$PATH_TO_MY_SECRET_MODULE
~~~~

## 2. Specific module {#spack_specific_module}

To develop a module, you need install mofem-users-modules, clone from the
repository which you like to work with, set up configuration and build the
code.
~~~~
spack install mofem-users-modules 
cd $HOME
spack view --verbose symlink -i um_view mofem-users-modules
export PATH=$PWD/um_view/bin:$PATH
mkdir $HOME/mod_developer
cd mod_developer/
git clone -b develop https://bitbucket.org/likask/mofem_um_minimal_surface_equation minimal_surface_equation
spack setup mofem-minimal-surface-equation@develop build_type=RelWithDebInfo
./spconfig.py -DMOFEM_DIR=$HOME/um_view \
-DEXTERNAL_MODULE_SOURCE_DIRS=$HOME/mod_developer \
$HOME/um_view/users_modules
make -j4
~~~~

## 3. Basic users modules {#spack_basic_users_modules}

To develop the basic users modules, the procedure is similar to one shown above.
One needs
to install mofem-cephas, clone source code, run configuration and finally
make the code.
~~~~
spack install mofem-cephas 
mkdir $HOME/um_developer
cd $HOME/um_developer/
git clone -b develop https://likask@bitbucket.org/mofem/users-modules-cephas.git 
spack setup mofem-users-modules@develop build_type=RelWithDebInfo
spack view --verbose symlink -i um_view mofem-cephas
export PATH=$PWD/um_view/bin:$PATH
./spconfig.py -DMOFEM_DIR=$HOME/um_view users-modules-cephas/ 
make -j4
~~~~

## 4. Core libraries {#spack_core_libraries}

If you are going to develop MoFEM's core library, it means that you are a core
developer and you can install MoFEM directly from the source.

Create mofem_install folder in the home directory and clone MoFEM repository:
~~~~~
mkdir $HOME/mofem_install
cd $HOME/mofem_install
git clone -b develop --recurse-submodules https://bitbucket.org/likask/mofem-cephas.git mofem-cephas
~~~~~
and kick-start installation of the core library:
~~~~~
cd $HOME/mofem_install
mkdir lib
cd lib/
spack install --only dependencies mofem-cephas 
spack setup mofem-cephas@develop copy_user_modules=False build_type=RelWithDebInfo
./spconfig.py $HOME/mofem_install/mofem-cephas/mofem/
make -j4
ctest
make install
~~~~~
Note that in addition to `build_type` another specification of the build configuration (*spec*) was used: `copy_user_modules=False `. 

Next, install users modules
~~~~~
cd $HOME/mofem_install
mkdir um
cd um/
spack view --verbose symlink -i um_view mofem-cephas@develop
export PATH=$PWD/um_view/bin:$PATH
mkdir build 
cd build/
spack setup mofem-users-modules@develop ^mofem-cephas@develop copy_user_modules=False build_type=RelWithDebInfo
./spconfig.py -DMOFEM_DIR=../um_view $HOME/mofem_install/mofem-cephas/mofem/users_modules
make -j4
ctest
make install
~~~~~

In the `spack setup` command of the snippet above `^` is a *dependency spec*, i.e. a descriptor defining the dependency of the package that we are currently installing on another package (in this case, on the core library). Note that `^` is
followed by the package name and its own *specs* for a particular
version, i.e. *specs* are defined recursively, see
[Spack manual page](https://spack.readthedocs.io/en/latest/basic_usage.html#specs-dependencies)
for more details. For example, if the *Debug* configuration is needed for both the core library and for user modules, then
`build_type=Debug` must be specified for both packages in `spack setup` command above. Note that in the considered case *specs* for `mofem-cephas@develop` package must coincide with the ones provided during the installation of the core library discussed above.
Alternatively, a particular version of an already installed library can be
defined by its unique ID, which can be determined using the instructions
given below in [MoFEM package versions](#spack_mofem_package_versions). This unique ID
is to be provided after the *dependency spec* `^` with a `/` (*slash*)
preceding, e.g.:
~~~~~ 
spack setup mofem-users-modules@develop copy_user_modules=False build_type=Debug ^/yk45ivx 
~~~~~

<!-- 
You can add extended users modules to
`$HOME/mofem_install/mofem-cephas/mofem/users_modules`. To include this in
the build process:
~~~~~
./spconfig.py -DMOFEM_DIR=../um_view $HOME/mofem_install/mofem-cephas/mofem/users_modules
make -j4
ctest
make install
~~~~~ 
-->

Alternatively, you can add your users modules to an independent folder and run
the snippet below
~~~~~
./spconfig.py \
-DEXTERNAL_MODULE_SOURCE_DIRS=$PATH_TO_MY_SECRET_MODULE \
-DMOFEM_DIR=../um_view \
$HOME/mofem_install/mofem-cephas/mofem/users_modules
make -j4
ctest
make install
~~~~~

# Installation on specific servers {#spack_servers}

## RDB development server {#spack_rdb}

On RDB server we are running Ubuntu 12.04.5 LTS, which is very old and does
not have compilers which can compile the C++14 code required by MoFEM and
some dependent packages. There are also problems with
[curl](https://en.wikipedia.org/wiki/CURL) version on Ubuntu 12.04.5 LTS is
compiled with old
[SSL](https://en.wikipedia.org/wiki/Transport_Layer_Security) protocol and
fetching files from some servers does not work. In order to solve those problems,
we will use a local mirror, which has previously downloaded packages and install
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

## Server Buckethead {#spack_buckethead}

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

We need to set up compilers since the location of standard `gcc@6.4.0` libraries
is in an unusual place when loaded by the module. In order to do that you have
to edit file `.spack/linux/compilers.yaml` in your home directory
~~~~~
     1	compilers:
     2	- compiler:
     3	    environment: {}
     4	    extra_rpaths: []
     5	    flags: {}
     6	    modules: []
     7	    operating_system: centos7
     8	    paths:
     9	      cc: /usr/bin/gcc
    10	      cxx: /usr/bin/g++
    11	      f77: null
    12	      fc: null
    13	    spec: gcc@4.8.5
    14	    target: x86_64
    15	- compiler:
    16	    environment: {}
    17	    extra_rpaths:
    18	    - /software/compilers/gcc/6.4.0/lib64
    19	    flags: {}
    20	    modules: []
    21	    operating_system: centos7
    22	    paths:
    23	      cc: /software/compilers/gcc/6.4.0/bin/gcc
    24	      cxx: /software/compilers/gcc/6.4.0/bin/g++
    25	      f77: /software/compilers/gcc/6.4.0/bin/gfortran
    26	      fc: /software/compilers/gcc/6.4.0/bin/gfortran
    27	    spec: gcc@6.4.0
    28	    target: x86_64
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
-pc_factor_mat_solver_type mumps \
-ksp_monitor \
-my_order 2 2>&1 | tee log
~~~~~
and run it as follows
~~~~~
qsub job_spack
~~~~~
Results of the analysis are located in $HOME/um_view/elasticity. 

# Spack usage and configuration {#spack_usage_config}

## Spack configuration setup {#spack_config_setup}

The `spack setup` command generates a `spconfig.py` configuration file. This
script calls CMake and can be edited to change dependencie paths and options.
Instructions on how to build with Spack can be found
[here](https://spack.readthedocs.io/en/latest/workflows.html?highlight=spconfig.py#build-with-spack).

## Enabling tests

MoFEM comes with built-in tests but by default they are off. To enable the tests
modify `spconfig.py` before building.

For `mofem-cephas` or core libraries change `OFF` to `ON`:
~~~~~~
DMOFEM_BUILD_TESTS=OFF
~~~~~~
And the same change for `users_modules`:
~~~~~~
DMOFEM_UM_BUILD_TESTS=OFF
~~~~~~

The following command will do this. Note this is for `GNU sed` and macOS
uses `BSD sed`. To use `GNU sed` install it via `homebrew install gnu-sed`
and call it with `gsed` instead of `sed`. Linux defaults to `GNU sed`.
~~~~~~
-sed -i 's/DMOFEM_UM_BUILD_TESTS=OFF/DMOFEM_UM_BUILD_TESTS=ON/ ; s/DMOFEM_BUILD_TESTS=OFF/DMOFEM_BUILD_TESTS=ON/' spconfig.py
~~~~~~

## Setting build type and compiler flags {#spack_build_type}

The build type affects performance and the availability of debugging symbols.
Select an appropriate type. With spack during installation, you can set the
`build_type`:

* `Release`
* `Debug`
* `RelWithDebInfo`

By default Spack will use `RelWithDebInfo`. 

For example:
~~~~~~
spack install mofem-cephas build_type=Debug
~~~~~~

Or when when using `spack setup` for developers:
~~~~~~
spack setup mofem-cephas@develop copy_user_modules=False build_type=Debug
~~~~~~

An example compiler flag:
~~~~~~
spack install mofem-cephas cppflags="-march=native -O3"
~~~~~~
## MoFEM package versions {#spack_mofem_package_versions}

Spack is capable of managing multiple versions of MoFEM. It stores each version
under a unique ID. The following command displays installed MoFEM packages with
their unique ID, build type and other info:
~~~~~~
spack find -lvd | grep mofem
~~~~~~

An example output:
~~~~~~
yk45ivx    mofem-cephas@develop+adol-c build_type=Debug ~copy_user_modules+med~slepc+tetgen
mco5jh5    mofem-cephas@develop+adol-c build_type=RelWithDebInfo ~copy_user_modules+med~slepc+tetgen
7f4u3sj    mofem-cephas@develop+adol-c build_type=Release ~copy_user_modules+med~slepc+tetgen
jtolxoh    mofem-users-modules@develop build_type=Debug ~copy_user_modules
yk45ivx        ^mofem-cephas@develop+adol-c build_type=Debug ~copy_user_modules+med~slepc+tetgen
~~~~~~

Proceeding the name of the package is its ID. e.g. `yk45ivx`. This is important
when specifying which package version you wish to use.

## Adding MoFEM extension to Spack {#spack_adding_package}

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

## Basic settings in config.yaml {#spack_config}

Spack's basic configuration options are set in config.yaml. You can see the
default settings by *$HOME/.spack/config.yaml*, for example, you can set the
default number of jobs to use when running `make` in parallel. If set to 4,
for example, `spack install` will run `make -j4`. This can be done by adding line to *config.yaml*. 
Example file should look like this:
~~~~~ 
config:	
	build_jobs: 4 
~~~~~ 
For more details see
[here](https://spack.readthedocs.io/en/latest/config_yaml.html?highlight=-jobs#)

## Mirrors {#spack_mirrors}

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

## How to find packages that were only installed today?

Run command
~~~~~
ls -lhd `spack find -lp  | grep mofem | awk '\''{print $3}'\''`
~~~~~

# Contact {#spack_contact}

Any problems with this installation, please contact us by
[mofem-group@googlegroups.com](https://groups.google.com/forum/#!forum/mofem-group)
or on Slack [MoFEM Slack](https://mofem.slack.com/).
