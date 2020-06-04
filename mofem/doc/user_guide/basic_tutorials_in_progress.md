Basic tutorials (in progress) {#basic_tutorials_in_progress}
==========================================================

[TOC]

The tutorials contain a number of basic programs that are built on top of each
other. Each of the tutorials has the following structure

1. **Introduction**: Intended learning outcome, problem being solved,
mathematical equations and derivations leading to finite element implementation
2. **Implementation**: Detailed explanation of the source code
3. **Results**: How to run the program, output, visualisation, interpretation
and comments, possible extensions.
4. **Plain program**: Full source code without extended comments

The source code and corresponding binary files of the tutorials are located at
the following directories

- Source code:
  *$HOME/mofem_install/mofem-cephas/mofem/users_modules/basic_finite_elements/tutorials*
- Binary files (build directory):
  *$HOME/mofem_install/um/build/basic_finite_elements/tutorials* 

Each tutorial in this page includes code name and keywords for quick reference
and search within the browser.

# Mesh generation

This part presents how to create a MoFEM-compatible mesh from open-source mesh
generators, e.g. Gmsh, Salome

- \ref basic_tutorials_mesh_generation_2d (Keywords: *block definition, config file, read_med*)
- \ref basic_tutorials_mesh_generation_3d (Keywords: *block definition, config file, read_med*)



# Scalar-field problems

In this part of the tutorials, you will be guided through the process of getting
familiar with MoFEM. Firstly, solving Poisson's equation and its variants. This
equation is of scalar-field problem meaning the field value that needs to be
solved has *only one* component. In these tutorials, you will:
- solve partial differential equations in MoFEM from the simplest case to more complicated ones
- get familiar with the concept of User Data Operator (UDO) and how to
implement new operators and \b push them to the main program
- get familiar with Simple Interface among other \ref basic_lessons1_interfaces
- expand the code to solve 3D problem from 2D one.
- be able to use PETSc linear solver ([KSP](https://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/KSP/index.html)),
nonlinear solver ([SNES](https://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/SNES/index.html)),
and time-stepping solver ([TS](https://www.mcs.anl.gov/petsc/documentation/tutorials/ECP19/ECP19_TS.pdf))
- get familiar with ([PCFIELDSPLIT](http://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/PC/PCFIELDSPLIT.html))
pre-conditioner to solve the problem with multiple fields.

## Linear and nonlinear Poisson's equation

Let's start to learn MoFEM by solving the linear Poisson's equation in 2D with
homogeneous boundary condition. At the end of this tutorial, there will be
instruction on how to effortlessly and quickly switch the implementation to
solve the equation in 3D with few changes in the code

- \ref basic_tutorials_poisson_homogeneous (Keywords: *Simple Interface, KSP solver, mofem_part, mbconvert, Paraview*)

Based on the implementation of the previous step, let's explore how to expand
the code to cover the non-homogeneous boundary condition

- \ref basic_tutorials_poisson_nonhomogeneous (Keywords: *Simple Interface, KSP solver, least-square approx.*)

Below is an alternative way to solve the equation with non-homogeneous boundary
condition using Lagrange multiplier which means making use of an additional field 

- \ref basic_tutorials_poisson_lagrange_multiplier (Keywords: *Simple Interface, PCFIELDSPLIT block solver*)

That was for linear problem using KSP solver. Now, let's try to solve a
nonlinear Poisson's equation by mean of SNES solver

- \ref basic_tutorials_poisson_nonlinear (Keywords: *Simple Interface, SNES solver*)

Let's concretise what you have learnt by solving a variant of the nonlinear
Poisson's equation, namely Minimal surface equation

- \ref basic_tutorials_minimal_surface_equation (Keywords: *Simple Interface, SNES solver*)

## Time-dependent problems

You have finished some Poisson's equations. By now, you have general idea how
MoFEM works and how to implement UDOs which is the crucial part in MoFEM.

Experience with all of the implementation for linear/nonlinear Poisson's equation
would be helpful as you progress to the next stage which involves solving
time-dependent heat equation, wave equation, and reaction-diffusion equation
using PETSc time-stepping solver
([TS](https://www.mcs.anl.gov/petsc/documentation/tutorials/ECP19/ECP19_TS.pdf))

- \ref basic_tutorials_heat_equation (Keywords: *Simple Interface, TS solver, implicit scheme*)

- \ref basic_tutorials_wave_equation (Keywords: *Simple Interface, TS solver, implicit scheme*)

- \ref basic_tutorials_reaction_diffusion_equation (Keywords: *Simple Interface, TS solver, Implicit-Explicit (IMEX) scheme*)


# Vector-field problems

This part of the tutorials contains sample programs that deals with vector-field
problems meaning the field value that needs to be solved has *three* components.

## Elasticity

This includes linear and nonlinear elasticity problems.

- \ref VEC-1: Linear elasticity
- \ref VEC-2: Nonlinear elasticity

# Complex and mixed field problems

- \ref MIX-1: Acoustic wave
- \ref MIX-2: Mixed problem

# Advanced topics

- \ref ADV-1: Plasticity
- \ref ADV-2: Contact
- \ref ADV-3: Fracture
- \ref ADV-3: Problem with various fields and finite element spaces
- \ref ADV-3: Problem with mixed dimension (1D, 2D, 3D)
