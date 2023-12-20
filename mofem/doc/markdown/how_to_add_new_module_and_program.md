How to add a new module and program {#how_to_add_new_module_and_program}
==========================================================

[TOC]

This tutorial assumes that you have installed MoFEM with developer version using either instructions in \ref install_docker_jupyterhub or the script provided in \ref installation. Therefore, you should have paths
for the source code and the binary files (e.g. for the `Release`  build type) of the Core Library and Users Modules as follows:

- Core Library
  - Source code: `$HOME/mofem_install/mofem-cephas/`
  - Binary files (build directory): `$HOME/mofem_install/mofem-cephas/core-build-Release-*`
- Users Modules
  - Source code: `$HOME/mofem_install/mofem-cephas/mofem/users_modules/`
  - Binary files (build directory): `$HOME/mofem_install/mofem-cephas/mofem/users_modules/um-build-Release-*`

> In build directories described above, symbol `*` stands for the `hash` of the particular build, e.g. `5sehreo`

# How to add a new module

## An existing module

To add a new MoFEM module that is already developed, for example [Solid Shell
Module](https://bitbucket.org/likask/mofem_um_solid_shell_prism_element/src/master/), you may wish to follow these steps

- Locate Users Modules source code
```
cd $HOME/mofem_install/mofem-cephas/mofem/users_modules/
```

- Clone Solid Shell Module
```
git clone https://bitbucket.org/likask/mofem_um_solid_shell_prism_element.git
```

- Compile codes for all modules
```
cd $HOME/mofem_install/mofem-cephas/mofem/users_modules/um-build-Release-*
make -j4
```

> If you have several builds of the Users Modules, e.g. `Debug` and `Release`, you will need to substitute the symbol `*` in the command above by the `hash` of the build you want to recompile.

- By now you should see all the executables generated in the folder within `um-build-Release-*`
```
cd mofem_um_solid_shell_prism_element
```

- Run the unit test (if the module has one)
```
spack load cmake
ctest
```

If you did everything correctly, you should see the paths for %MoFEM Solid Shell Module are as follows

- MoFEM Solid Shell Module
  - Source code: `$HOME/mofem_install/mofem-cephas/mofem/users_modules/mofem_um_solid_shell_prism_element`
  - Binary files (build directory): `$HOME/mofem_install/mofem-cephas/mofem/users_modules/um-build-Release-*/mofem_um_solid_shell_prism_element`

## New module for your own purpose

The quickest way to create your own module is to replicate a module that is already
available in MoFEM and then replace the existing source codes with your own. For
example, you can replicate MoFEM Solid Shell Module following general steps as follows:

- Clone the module to somewhere outside of MoFEM installation directory
- Change file names and descriptions
- Be careful with `CMakeLists.txt`, `InstalledAddModule.cmake`. Make sure all
  the names and paths are consistent with the new module name of your choice
- Delete unnecessary source code files
- Once everything done, create a new repository of your own on
  [Github](https://github.com) or [Bitbucket](https://bitbucket.org)
- Upload your new files to the remote repository
- Follow steps in the previous section to clone and setup your new module in
  MoFEM (How to add a new module - An existing module).

\note If you have issues creating new module of your own, please contact us at [MoFEM
Q&A](https://groups.google.com/forum/#!categories/mofem-group). 


# How to add a new program

The quickest way to add a new program to an existing module in your
MoFEM installation directory is to replicate a program that already available in
the same module of your choice. For example, you can replicate `vec-0`
program in MoFEM tutorials following general steps as follows

- Go to:
  `$HOME/mofem_install/mofem-cephas/mofem/users_modules/tutorials/`
- Copy tutorial, e.g. `vec-0` and paste it in the same place
- Change file names and descriptions
- Be careful with `CMakeLists.txt`. Make sure all the names and paths are consistent with the program name of your choice
- Delete unnecessary source code files
- Once everything done, compile source code by
```
cd $HOME/mofem_install/mofem-cephas/mofem/users_modules/um-build-Release-*/tutorials
make -j4
```
- You will see a new directory containing your new program

\note If you have issues creating a new program of your own, please contact us at [MoFEM
Q&A](https://groups.google.com/forum/#!categories/mofem-group). 