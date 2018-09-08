Installation with Spack (Recommended for HPC, Linux & Mac OS X) {#install_spack}
==========================================================

Spack is "A flexible package manager that supports multiple versions, configurations, platforms, and compilers." - [https://spack.io](https://spack.io)

A community of Spack users and developers contribute to developing package configurations for a range of platforms, including macOS and Linux. This creates a consistant build setup for supported scientific packages.

MoFEM can be deployed and developed using Spack, and is the recommended way of installing. If you have any problems, feedback or would like to suggest corrections,
please email
[mofem-group@googlegroups.com](https://groups.google.com/forum/#!forum/mofem-group).

Note: Spack compiles packages from source, pre-compiled binaries are not available. 

[TOC]

# Prerequisites {#spack_prerequisites}

The insallation of MoFEM requires [git](https://www.atlassian.com/git/tutorials/what-is-git),
[curl](https://en.wikipedia.org/wiki/CURL) and C++ and Fortran compilers (GCC and/or clang). They must also be available in your PATH.

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

If not already done so, install Xcode and command-line tools:

~~~~~
xcode-select --install
sudo xcodebuild -license accept
~~~~~

Install [homebrew](https://brew.sh) package manager:
~~~~~
/usr/bin/ruby -e "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install)"
~~~~~

Install packages through hombrew:
~~~~~
brew install curl git gfortran
~~~~~

Check PATH of Fortran compiler:
~~~~~
which gfortran
~~~~~

If not there, then add to PATH.

Note: recent releases of macOS stopped shipping a Fortran compiler and so require [Mixed
Toolchains](http://spack.readthedocs.io/en/latest/getting_started.html#mixed-toolchains). The installing of gfortran through homebrew is one way of solving this.

# Spack setup {#spack_setup}

Retrieve Spack's configuration for MoFEM:
~~~~~~
git clone --single-branch -b mofem https://github.com/likask/spack.git
~~~~~~

Initialise Spack's environment variables:
~~~~~~
. spack/share/spack/setup-env.sh
~~~~~~

Spack's environment variables will be lost at the close of a terminal session. Consider adding to `.bashrc` or `.profile`

Install packages required by Spack:
~~~~~~
spack bootstrap
~~~~~~

# User only {#spack_user}

This will install MoFEM with all of its required dependencies from source using Spack. This will take considerable time > 1 hour. Those seeking to develop MoFEM see [developer's](#spack_developers) instructions.

## Install base users modules MoFEM

Install basic MoFEM with required dependences. This will take considerable time:
~~~~~~
spack install mofem-users-modules
~~~~~~

MoFEM is extendable by adding user modules. More modules are available as extensions in Spack. This will show available extensions:
~~~~~~
spack extensions mofem-cephas
~~~~~~

The extensions can be installed using `spack install <extension>`.

Create a `um_view` directory to conveniently access the installed users-modules. This should be created in an appropriate directory:
~~~~~~
spack view --verbose symlink -i um_view mofem-cephas
spack activate -v um_view mofem-users-modules
~~~~~~

This filesystem view is a single directory tree that is the union of the
directory hierarchies of a number of installed packages; it is similar to the
directory hierarchy that might exist under */usr/local*. The files of the
view’s installed MoFEM packages are brought into the view by symbolic or hard
links, referencing the original Spack installation. Different 'views' can be created depending on the MoFEM version you wish to access.
Add the new 'view' *bin* directory to your PATH. e.g.:
~~~~~~
export PATH=$PWD/um_view/bin:$PATH
~~~~~~

## Test elasticity module

~~~~~~
cd um_view/elasticity

mpirun -np 2 ./elasticity \
-my_file LShape.h5m \
-my_order 2 \
-ksp_type gmres \
-pc_type lu -pc_factor_mat_solver_package mumps \
-ksp_monitor
mbconvert out.h5m out.vtk
~~~~~~

Open the output VTK file in [ParaView](https://www.paraview.org). You can
install ParaView using Spack, [directly from the website](https://www.paraview.org/download/) or through your OS package manager.

## Install additional users modules

For example, the fracture module can be installed by:
~~~~~~
spack install mofem-fracture-module
cd $HOME
spack view --verbose symlink -i um_view_foo mofem-cephas
spack activate -v um_view_foo mofem-fracture-module
~~~~~~
or minimal surface equation tutorial module:
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

You can install package and run tests, for example
~~~~~~
spack install --test=root -j 4 mofem-cephas
spack install --test=root -j 4 mofem-users-modules
spack install --test=root -j 4 mofem-fracture-module
~~~~~~

## Setting build type and compiler flags

The build type affects performance and the availability of debugging symbols. Select an appropropriate type. With spack during installation, you can set the `build_type`:

* `Release`
* `Debug`
* `RelWithDebInfo`

By default Spack will use `RelWithDebInfo`. 

For example:
~~~~~~
spack install mofem-cephas build_type="Debug"
~~~~~~

An example compiler flat:
~~~~~~
spack install mofem-cephas cppflags="-march=native -O3"
~~~~~~

