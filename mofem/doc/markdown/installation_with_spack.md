Installation with Spack (Linux, macOS, HPC) {#install_spack}
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

You will copy the script to the directory where MoFEM will be installed and run
it from the terminal, for example
~~~~~~
./install_mofem_user.sh
~~~~~~ 

\note Running the scripts may require user password for sudo privileges.

## Scripts on secure server

If you are going to install MoFEM on a secure server, or in server without
access to the internet (in a train or during long flight), you can download
spack and spack mirror first, and then start the installation.

Download spack
~~~~~~
curl -L https://api.github.com/repos/likask/spack/tarball/mofem \
--output spack.tgz
~~~~~~
and download mirror
~~~~~~
curl -L  http://mofem.eng.gla.ac.uk/downloads/mirror_v0.16.tar.gz \
--output mirror.tgz
~~~~~~
and then you can install MoFEM
~~~~~
./install_mofem_user.sh
~~~~~

# Prerequisites {#spack_prerequisites}

The installation of MoFEM requires
[git](https://www.atlassian.com/git/tutorials/what-is-git),

System requirements: Minimum 8GB RAM. On Linux ensure this RAM is free.

##Ubuntu/Debian

Install the following packages:
~~~~~
apt-get update \
&& apt-get install -y --no-install-recommends \
build-essential \
ca-certificates \
coreutils \
environment-modules \
python \
python3-distutils \
gfortran \
curl \
git \
cmake \
doxygen \
vim 
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
/bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"
~~~~~

Install packages through `homebrew`:
~~~~~
brew install curl git python gcc@9 cmake autoconf automake libtool doxygen graphviz pkg-config
~~~~~

If it is not already in the `PATH`, you should add it there.

Note: recent releases of macOS stopped shipping a Fortran compiler and therefore
require [Mixed
Toolchains](http://spack.readthedocs.io/en/latest/getting_started.html#mixed-toolchains).
The installing of gfortran through homebrew is another way of solvingr this.

# Spack setup {#spack_setup}

Make `mofem_install` director
~~~~~~
mkdir -p $HOME/mofem_install
cd $HOME/mofem_install
~~~~~~

Retrieve Spack for MoFEM:
~~~~~~
git clone --single-branch -b develop_spack_v0.16.1 https://github.com/likask/spack.git
~~~~~~

Initialise Spack's environment variables:
~~~~~~
. spack/share/spack/setup-env.sh
~~~~~~

Download spack packages in the mirror necessary to install MoFEM
~~~~~~
mkdir -p mofem_mirror &&
curl -s -L  http://mofem.eng.gla.ac.uk/downloads/mirror_v0.16.tar.gz \
| tar xzC $PWD/mofem_mirror  --strip 1
spack mirror add mofem_mirror $PWD/mofem_mirror
~~~~~~
Above option is not needed. If you do not download and add a mirror, packages
will be downloaded from the internet. However, locations of libraries change,
or some server could be temporarily down. Using spack mirror making
installation resistant to those problems.

Spack's environment variables will be lost when the terminal session is closed.
Consider adding the previous command to your `.bash_profile` or `.bashrc`, e.g.:
~~~~~~
echo ". $HOME/spack/share/spack/setup-env.sh" >> ~/.bash_profile
~~~~~~

If you using Big Sur or Catalina, or newer Mac OX, you have to add config to
`.zshrc`, e.g.:
~~~~~~
echo ". $HOME/spack/share/spack/setup-env.sh" >> ~/.zshrc
~~~~~~

Finally, install packages required by Spack:
~~~~~~
spack compiler find
spack external find
~~~~~~

If you using system where gfortant v10 is installed, some packages like
openmpi will not compile on Mac OSX. You can check this running code
~~~~~
gfortran -v
~~~~~
as reusult if you get
~~~~~
gcc version 10.2.0 (GCC) 
~~~~~
t means that you have version 10. This is temporary problem and will be
fixed over the time, once various patches and fixes will be applied to those
libraries. In mean time you can fix that problem by editind `compilers.yaml`
in Mac OSX located in `~/.spack/darwin/compilers.yaml` to set version 9 of
gfortran compiler,
~~~~~~~
- compiler:
    spec: apple-clang@12.0.0
    paths:
      cc: /usr/bin/clang
      cxx: /usr/bin/clang++
      f77: /usr/local/bin/gfortran-9
      fc: /usr/local/bin/gfortran-9
    flags: {}
    operating_system: macos
    target: x86_64
    modules: []
    environment: {}
    extra_rpaths: []
~~~~~~~
Note that fortran compiler is set to version `9`, as follows
~~~~~~~
      f77: /usr/local/bin/gfortran-9
      fc: /usr/local/bin/gfortran-9
~~~~~~~

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
spack view --verbose symlink -i um_view mofem-users-modules
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

Consider also adding this command to your `.bash_profile` or `.bashrc`, or `.zshrc`, e.g.: 
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
mbconvert out_skin.h5m out_skin.vtk
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
spack install --test=root -j 1 mofem-cephas
spack install --test=root -j 1 mofem-users-modules
spack install --test=root -j 1 mofem-fracture-module
~~~~~~

# Developers {#spack_developers}

MoFEM can be developed in many different ways, here you can follow to step
installation for core and user modules developers:

1. [Core libraries](#spack_core_libraries)
2. [Install users modules](#spack_users_modules)

The developer installation requires knowing more about how Spack works. See [Spack usage and configuration](#spack_usage_config)
before proceeding with the installation. In particular, the instructions below
will use so-called *specs* to obtain the desired build configuration, see
[Spack manual page](https://spack.readthedocs.io/en/latest/basic_usage.html#specs-dependencies)
for more details. For example, the default Spack specifier for the build type
will be given explicitly here: `build_type=RelWithDebInfo`, however keep in mind
that two other build types can be specified in the same way:
`build_type=Release` or `build_type=Debug`, see also [Change the build_type](#spack_build_type).

You can skip the first installation method and jump to second if you are not
going to be the core MoFEM developer. However, to avoid many installation
paths, and avoid potential bugs which are difficult to reproduce, we
recommend following installation step 1, and step 2.

## 1. Core libraries {#spack_core_libraries}

If you are going to develop MoFEM's core library, it means that you are a core
developer and you can install MoFEM directly from the source.

Create mofem_install folder in the home directory and clone MoFEM repository:
~~~~~
mkdir $HOME/mofem_install
cd $HOME/mofem_install
git clone \
  -b develop \
  --recurse-submodules https://bitbucket.org/likask/mofem-cephas.git \
  mofem-cephas
~~~~~
and kick-start installation of the core library. First install all dependencies,
~~~~~
spack install --only dependencies mofem-cephas ^petsc+X
~~~~~
and then core MoFEM library
~~~~~
spack dev-build \
  --source-path $HOME/mofem_install/mofem-cephas \
  --keep-prefix \
  --test root \
  mofem-cephas@develop~copy_user_modules ^petsc+X
~~~~~
Note that in addition to `build_type` another specification of the build configuration (*spec*) was used: `copy_user_modules=False `. 

You can do partial install, if before of after some installation phase
using command line option `-b BEFORE` or `-u UNTIL`, for example if you like
to investigate build issues, you can do,
~~~~~
spack dev-build \
  -b build \
  --source-path $HOME/mofem_install/mofem-cephas \
  --keep-prefix \
  --test root \
  mofem-cephas@develop~copy_user_modules ^petsc+X
~~~~~
You can find build director in `$HOME/mofem_install/mofem-cephas`, its name
will start as `core-build-WithDebugInfo-`, for example `core-build-WithDebugInfo-vhv7opa`. 
Note that installation at that point is partial

If installation is successfully, by executing, 
~~~~~
spack find -lv mofem-cephas
~~~~~
you should see
~~~~~
==> 1 installed package
-- darwin-macos-haswell / apple-clang@12.0.0-gfortran9 ----------
vhv7opa mofem-cephas@develop+adol-c~copy_user_modules~ipo+med+slepc+tetgen build_type=RelWithDebInfo dev_path=/Users/lukaszkaczmarczyk/mofem_install/mofem-cephas install_id=0
==> 1 installed package
-- darwin-macos-haswell / apple-clang@12.0.0 ----------
~~~~~

Also in directory `$HOME/mofem_install/mofem-cephas` you will find
build directory `$HOME/mofem_install/mofem-cephas/core-build-WithDebugInfo-vhv7opa`. Note 
that suffix is matching first column when executed column `spack find -lv mofem-cephas`.

You can now start develop code, if you 
`cd $HOME/mofem_install/mofem-cephas/core-build-WithDebugInfo-vhv7opa`, and
can run 
~~~~~
make -j4
ctest -D Experimental
make install
~~~~~
and do typical developer work.

You can install simultaneously debugging version of code, as follows
~~~~~
spack dev-build \
  --source-path $HOME/mofem_install/mofem-cephas \
  --keep-prefix \
  --test root \
  mofem-cephas@develop~copy_user_modules build_type=Debug ^petsc+X
~~~~~ 
## 2. Install users modules {#spack_users_modules}

Run first `spack find -lv mofem-cephas`, and pick core library
~~~~~
==> 2 installed packages
-- linux-ubuntu18.04-zen / gcc@7.5.0 ----------------------------
nnnvprd mofem-cephas@develop+adol-c~copy_user_modules~ipo+med+slepc+tetgen build_type=Debug dev_path=/home/lukasz/mofem_install/mofem-cephas install_id=0
pa3httg mofem-cephas@develop+adol-c~copy_user_modules~ipo+med+slepc+tetgen build_type=RelWithDebInfo dev_path=/home/lukasz/mofem_install/mofem-cephas install_id=0
~~~~~
For example, if you like to intall users modules against build `RelWithDebInfo` pick
second row, and build `pa3httg`. Next, install users modules
~~~~~
spack dev-build \
  --test root  \
  --source-pat $HOME/mofem_install/mofem-cephas/mofem/users_modules \
  mofem-users-modules@develop build_type=RelWithDebInfo \
  ^/pa3httg
~~~~~
You can do partial install, if before of after some installation phase
using command line option `-b BEFORE` or `-u UNTIL`, for example if you like
to investigate build issues, you can do,
~~~~~
spack dev-build \
  -u configure \
  --test root  \
  --source-pat $HOME/mofem_install/mofem-cephas/mofem/users_modules \
  mofem-users-modules@develop build_type=RelWithDebInfo \
  ^/pa3httg
~~~~~
After that point, `spack` will configure environment and `cmake`, and beyond that
point you can run build and install
~~~~~
cd $HOME/mofem_install/mofem-cephas/mofem/users_modules/um-build-WithDebugInfo-c6ts56b
make -j4
make install
~~~~~
and use code.

Once installation is successfully, you can execute `spack find -lv mofem-users-modules`, and
as result you will get
~~~~~
==> 1 installed package
-- darwin-macos-haswell / apple-clang@12.0.0-gfortran9 ----------
c6ts56b mofem-users-modules@develop+copy_user_modules~ipo build_type=Debug dev_path=/Users/lukaszkaczmarczyk/mofem_install/mofem-cephas/mofem/users_modules install_id=0
~~~~~
Also, you will find director
`$HONE/mofem_install/mofem-cephas/mofem/users_modules/um-build-WithDebugInfo-c6ts56b`.
In that directory is build directory particular versions of
mofem-users-modules. In it you can do typical developer work, `make`, `ctest`,
and `make install`. You can add other modules to that directory if needed.

In the `spack dev-build` command of the snippet above `^` is a *dependency
spec*, i.e. a descriptor defining the dependency of the package that we are
currently installing on another package (in this case, on the core library).
Note that `^` is followed by the package name and its own *specs* for a
particular version, i.e. *specs* are defined recursively, see [Spack manual
page](https://spack.readthedocs.io/en/latest/basic_usage.html#specs-dependencies)
for more details. For example, if the *Debug* configuration is needed for
both the core library and for user modules, then `build_type=Debug` must be
specified for both packages in `spack dev-build` command above. Note that in the
considered case *specs* for `mofem-cephas@develop` package must coincide with
the ones provided during the installation of the core library discussed
above. Alternatively, a particular version of an already installed library
can be defined by its unique ID, which can be determined using the
instructions given below in [MoFEM package
versions](#spack_mofem_package_versions). This unique ID is to be provided
after the *dependency spec* `^` with a `/` (*slash*) preceding, e.g.:
~~~~~ 
spack spack dev-build \
  --test root  \
  --source-pat $HOME/mofem_install/mofem-cephas/mofem/users_modules \
  mofem-users-modules@develop~copy_user_modules build_type=Debug \
  ^/nnnvprd
~~~~~

# Installation on specific servers {#spack_servers}

## Server Buckethead {#spack_buckethead}

### Installation {#spack_buckedhead_installation}

Buckethead is a Linux cluster running Centos7. More information about
Buckethead you will find 
[here](http://malmsteen.eng.gla.ac.uk/wiki/hpc-wiki/index.php/Resources#Buckethead). 

Load cluster modules
~~~~~
module load mpi/openmpi/3.1.4
module load gcc/9.2.0 
module load gridengine
~~~~~

Download spack
~~~~~
mkdir -p spack &&\
curl -s -L https://api.github.com/repos/likask/spack/tarball/mofem \
| tar xzC $PWD/spack --strip 1
~~~~~

Download packages mirror
~~~~~
mkdir -p mofem_mirror &&\
curl -s -L  http://mofem.eng.gla.ac.uk/downloads/mirror_v0.9.2.tar.gz \
| tar xzC $PWD/mofem_mirror  --strip 1
~~~~~

Set spack environment and mirror
~~~~~
. spack/share/spack/setup-env.sh
spack mirror add mofem_mirror $PWD/mofem_mirror
spack extensions find
spack compiler find
~~~~~
##### Setup packages and compiler

Edit file `.spack/linux/compilers.yaml`
~~~~~
compilers:
- compiler:
    environment: {}
    extra_rpaths: []
    flags: {}
    modules: []
    operating_system: centos7
    paths:
      cc: /usr/bin/gcc
      cxx: /usr/bin/g++
      f77: null
      fc: null
    spec: gcc@4.8.5
    target: x86_64
- compiler:
    environment: {}
    extra_rpaths: 
        - /software/compilers/gcc/9.2.0/lib64
    flags: {}
    modules: []
    operating_system: centos7
    paths:
      cc: /software/compilers/gcc/9.2.0/bin/gcc
      cxx: /software/compilers/gcc/9.2.0/bin/g++
      f77: /software/compilers/gcc/9.2.0/bin/gfortran
      fc: /software/compilers/gcc/9.2.0/bin/gfortran
    spec: gcc@9.2.0
    target: x86_64
~~~~~

Note that most of the file can be created by command
~~~~~
spack compiler find
~~~~~
however you have to put line 
~~~~~
    extra_rpaths: 
        - /software/compilers/gcc/9.2.0/lib64
~~~~~
which a enable linking std lib c++ libraries.

#### Running job to compile code

At that point, we can follow the standard installation procedure, as follows
~~~~~
spack bootstrap
spack install -j 2 -v --only dependencies mofem-cephas%gcc@9.2.0 ^openmpi@3.1.4%gcc@9.2.0
spack install mofem-users-modules
~~~~~

##### Creating view to compiled MoFEM modules

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
#$ -e ./stdout_job_$JOB_ID

# File to which standard output should be directed
#$ -o ./stderr_job_$JOB_ID

# E-mail address to which status updates should be sent
# N.B.: in an array job, a separate e-mail will be sent for each task!
#$ -M your.email@glasgow.ac.uk

# Events on which to send a status update
#$ -m beas

# Request for 1.0 GB of memory per task (needed on Miffy and Dusty)
#$ -l mem_tokens=1.0G

#$ -pe mofem-* 2 # where N is the number of processors required

# List of commands which do the actual work
echo "$NSLOTS received"
cat $PE_HOSTFILE

# Load compiler
module load mpi/openmpi/3.1.4

# List of commands which do the actual work
cd $HOME/um_view/elasticity
mpirun -np $NSLOTS \
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

### Downloading mirror with prerequisites

You can download mirror with all necessary packages from MoFEM repository and
untar and unzip to director
~~~~~
mkdir -p mofem_mirror &&
curl -s -L  http://mofem.eng.gla.ac.uk/downloads/mirror_spack12.tar.gz \
| tar xzC $PWD/mofem_mirror  --strip 1
~~~~~
Note that packages are expanded to directory `mofem_mirror`, and mirror is
made for MoFEM version v0.9.0.

Once mirror os downloaded, you can add it to your spack 
~~~~~
spack mirror add mofem_mirror $PWD/mofem_mirror
~~~~~

### Making onw mirror 

You need to create mirror first. Some sites may not have access to the
internet for fetching packages. These sites will need a local repository of
tarballs from which they can get their files. Spack has support for this with
mirrors. Look to Spack documentation to learn more about mirrors, see
[here](https://spack.readthedocs.io/en/latest/mirrors.html?highlight=mirror).

~~~~~
spack mirror create -D mofem-fracture-module
~~~~~
as a result, a directory is created with are prerequisites need to Spack
installation of mofem-fracture-module. You can now move that directory to 
your secure server with spack package.

On a secure server, you add mirror directory to Spack, for example
~~~~~
spack mirror add local_filesystem /$HOME/spack-mirror-2018-07-21
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
[mofem-group@googlegroups.com](https://groups.google.com/forum/#!forum/mofem-group).
