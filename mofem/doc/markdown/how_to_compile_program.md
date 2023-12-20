How to compile a program {#how_to_compile_program}
==========================================================

[TOC]

This brief tutorial assumes that you have installed MoFEM with developer version using either instructions in \ref install_docker_jupyterhub or the script provided in \ref installation. Therefore, you should have paths
for the source code and the binary files (e.g. for the `Release`  build type) of the Core Library and Users Modules as follows:

- Core Library
  - Source code: `$HOME/mofem_install/mofem-cephas/`
  - Binary files (build directory): `$HOME/mofem_install/mofem-cephas/core-build-Release-*`
- Users Modules
  - Source code: `$HOME/mofem_install/mofem-cephas/mofem/users_modules/`
  - Binary files (build directory): `$HOME/mofem_install/mofem-cephas/mofem/users_modules/um-build-Release-*`

> In build directories described above, symbol `*` stands for the `hash` of the particular build, e.g. `5sehreo`

# Compile the Core Library

You will need to compile the Core Library if you change any file related to C++ source codes that are located in 
```
  $HOME/mofem_install/mofem-cephas/mofem/
```

You can compile the Core Library by running the following command lines:

```
  cd $HOME/mofem_install/mofem-cephas/core-build-Release-*
  make -j4 install
```
where `-j4` indicates that you will use four processors to compile the code.

> If you have several builds of the Core library, e.g. `Debug` and `Release`, you will need to substitute the symbol `*` in the command above by the `hash` of the build you want to recompile.

# Compile Users Modules

You will need to compile Users Modules if you change any file related
to C++ source codes that are located in
```
  $HOME/mofem_install/mofem-cephas/mofem/users_modules/
```
The changes include those when you modify the codes of Basic Users Modules or when
you add a new module of your own and change its source code.

You can compile both the Basic User Modules and your own module (if any) by
running the following command lines:

```
  cd $HOME/mofem_install/mofem-cephas/mofem/users_modules/um-build-Release-*
  make -j4 install
```

> If you have several builds of the Users Modules, e.g. `Debug` and `Release`, you will need to substitute the symbol `*` in the command above by the `hash` of the build you want to recompile.

\note You only need to compile Core Library before compiling User Modules
if you have changes in the Core Library. If you only have changes in Users Modules
(either the Basic one or the one added by yourself), you just need to
compile the Users Modules.