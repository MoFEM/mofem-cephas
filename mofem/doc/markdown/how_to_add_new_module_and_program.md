How to add a new module and program {#how_to_add_new_module_and_program}
==========================================================

[TOC]

This tutorial assumes that you have installed MoFEM for developer using
the script provided in \ref installation. And therefore, you probably have paths
for the source code and the binary files (e.g. for release version) of the Core
Library and the Default User Module as follows

- Core Library
  - Source code: *$HOME/mofem_install/mofem-cephas/*
  - Binary files (build directory): *$HOME/mofem_install/lib_release/*
- Default User Module
  - Source code: *$HOME/mofem_install/mofem-cephas/mofem/users_modules/*
  - Binary files (build directory): *$HOME/mofem_install/um/build_release/*

# How to add a new module

## Module developed by someone else

To add a new MoFEM module that already developed, for example Solid Shell
Module, you may wish to follow these steps

- Locate Default User Module source code
```
cd $HOME/mofem_install/mofem-cephas/mofem/users_modules/
```

- Clone Solid Shell Module
```
git clone https://bitbucket.org/likask/mofem_um_solid_shell_prism_element.git
```

\note If you have issue with repo access for Solid Shell Module, please contact
us on [MoFEM Q&A](https://groups.google.com/forum/#!categories/mofem-group)

- Reconfigure build directory
```
cd $HOME/mofem_install/um/build_release/
./spconfig.py \
-DMOFEM_UM_BUILD_TESTS=ON \
-DFM_VERSION_MAJOR=0 \
-DFM_VERSION_MINOR=0 \
-DFM_VERSION_BUILD=0 \
-DMOFEM_DIR=../um_view $HOME/mofem_install/mofem-cephas/mofem/users_modules
```

- Compile codes for all modules
```
cd $HOME/mofem_install/um/build_release/
make -j4
```

- By now you should see all the executables generated
```
cd mofem_um_solid_shell_prism_element
```

- Run the unit test (if the module has one)
```
spack load cmake
ctest
```

If you did every thing correctly, you should see the paths for MoFEM
Solid Shell Module are as follows

- MoFEM Solid Shell Module
  - Source code: *$HOME/mofem_install/mofem-cephas/mofem/users_modules/mofem_um_solid_shell_prism_element*
  - Binary files (build directory): *$HOME/mofem_install/um/build_release/mofem_um_solid_shell_prism_element*

## New module for your own purpose

The quickest way to create your own purpose is to replicate a module that already
available in MoFEM. For example, you can replicate MoFEM Solid Shell Module
following general steps as follows

- Clone the module to somewhere outside of MoFEM
- Change file names and descriptions
- Be careful with `CMakeLists.txt`, `InstalledAddModule.cmake`. Make sure all
  the names and paths are consistent with the new module name of your choice
- Delete unnecessary source code files
- Once everything done, create a new repository of your own on
  [Github](https://github.com) or [Bitbucket](https://bitbucket.org)
- Upload your new files to the remote repository
- Follow steps in the previous section to clone and setup your new module in MoFEM.

\note If you have issues creating new module of your own, please contact us at [MoFEM
Q&A](https://groups.google.com/forum/#!categories/mofem-group). We are happy to help.


# How to add a new program

The quickest way to add a new program to an module that already exists in your
MoFEM installation directory is to replicate a program that already available in
the same module of your choice. For example, you can replicate `elasticity`
program in MoFEM Default User Module following general steps as follows

- Go to:
  $HOME/mofem_install/mofem-cephas/mofem/users_modules/basic_finite_elements/
- Copy `elasticity` and paste it in the same place
- Change file names and descriptions
- Be careful with `CMakeLists.txt`. Make sure all the names and paths are consistent with the program name of your choice
- Delete unnecessary source code files
- Once everything done, compile source code by
```
cd $HOME/mofem_install/um/build_release/basic_finite_elements/
make -j4
```
\note If you have issues creating new program of your own, please contact us at [MoFEM
Q&A](https://groups.google.com/forum/#!categories/mofem-group). We are happy to help.