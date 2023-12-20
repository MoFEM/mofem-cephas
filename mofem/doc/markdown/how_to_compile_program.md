How to compile a program {#how_to_compile_program}
==========================================================

[TOC]

This tutorial assumes that you have installed MoFEM with developer version using either instructions in \ref install_docker_jupyterhub
or the script provided in \ref installation. And therefore, you probably have paths
for the source code and the binary files (e.g. for `Release`  build type) of the
Core Library and the Basic User Modules as follows:

- Core Library
  - Source code: *$HOME/mofem_install/mofem-cephas/*
  - Binary files (build directory): *$HOME/mofem_install/mofem-cephas/core-build-Release-*
- Basic User Modules
  - Source code: *$HOME/mofem_install/mofem-cephas/mofem/users_modules/*
  - Binary files (build directory): *$HOME/mofem_install/mofem-cephas/mofem/users_modules/um-build-Release-*

# Compile the Core Library

You will need to compile the Core Library if you change any file related to C++ source codes that are located in 
```
  $HOME/mofem_install/mofem-cephas/
  $HOME/mofem_install/mofem-cephas/mofem/
```

You can compile the Core Library by running the following command lines 

```
  cd $HOME/mofem_install/mofem-cephas/core-build-Release-*
  make -j4 install
```
where `-j4` indicates you will use four processors to compile the code.

# Compile Basic Users Modules

You will need to compile Basic Users Modules if you change any file related
to C++ source codes that are located in
```
  $HOME/mofem_install/mofem-cephas/mofem/users_modules/
```
The changes include those when you modify the codes of Basic User Modules or when
you add a new module of your own and change its source code.

You can compile both the Basic User Modules and your own module (if any) by
running the following command lines:

```
  cd $HOME/mofem_install/mofem-cephas/mofem/users_modules/um-build-Release-*
  make -j4 install
```

\note You only need to compile Core Library before compiling Basic User Modules
if you have changes in the Core Library. If you only have changes in users modules
(either the Basic one or the one added by yourself), you just need to
compile the users modules.