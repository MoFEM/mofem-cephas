How to compile a program {#how_to_compile_program}
==========================================================

[TOC]

This tutorial assumes that you have installed MoFEM with developer version using
the script provided in \ref installation. And therefore, you probably have paths
for the source code and the binary files (e.g. for `RelWithDebInfo` -- Release
with Debug info -- build type) of the
Core Library and the Basic User Module as follows

- Core Library
  - Source code: *$HOME/mofem_install/mofem-cephas/*
  - Binary files (build directory): *$HOME/mofem_install/mofem-cephas/core-build-RelWithDebInfo-abcd1234*
- Basic User Module
  - Source code: *$HOME/mofem_install/mofem-cephas/mofem/users_modules/*
  - Binary files (build directory): *$HOME/mofem_install/mofem-cephas/mofem/users_modules/um-build-RelWithDebInfo-abcd1234*

# Compile the Core Library

You will need to compile the Core Library if you change any file related to C++ source codes that are located in 
```
  $HOME/mofem_install/mofem-cephas/
  $HOME/mofem_install/mofem-cephas/mofem/
```

You compile the Core Library by running the following command lines 

```
  cd $HOME/mofem_install/mofem-cephas/core-build-RelWithDebInfo-abcd1234
  make -j4
  make -j4 install
```
where `-j4` indicates you will use four processors to compile the codes.

\note Remember that compiling Core Library requires both `make` and `make install`



# Compile the Basic User Module

You will need to compile the Basic User Module if you change any file related
to C++ source codes that are located in
```
  $HOME/mofem_install/mofem-cephas/mofem/users_modules/
```
The changes include those when you modify codes for Basic User Module or when
you add a new module of your own and change its source code.

You can compile both the Basic User Module and your own module (if any) by
running the following command lines

```
  cd $HOME/mofem_install/mofem-cephas/mofem/users_modules/um-build-RelWithDebInfo-abcd1234
  make -j4
```

\note You only need to compile Core Library before compiling Basic User Module
if you have changes in Core Library. If you only have changes in user modules
(either the Basic one or the one added by yourself), you just need to
compile the Basic User Module.