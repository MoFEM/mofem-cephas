# **MoFEM (JosePH)** #

**MoFEM (Mesh Oriented Finite Element Method)** is a new finite element analysis code tailored for the solution of multi-physics problems with arbitrary levels of approximation, different levels of mesh refinement and optimised for high-performance computing. MoFEM is the blend of the **[Boost](http://www.boost.org) MultiIndex** containers, **[MOAB](https://trac.mcs.anl.gov/projects/ITAPS/wiki/MOAB)** (Mesh Oriented Database) and **[PETSc](http://www.mcs.anl.gov/petsc/)** (Portable, Extensible Toolkit for Scientific Computation). MoFEM is developed in C++ and it is open-source software under the [GNU](http://www.gnu.org/licenses/) Lesser General Public License. The current version of MoFEM has full support for **[CUBIT](https://cubit.sandia.gov/)**/**[TRELIS](http://csimsoft.com/)** v13 for pre-processing and **[ParaView](http://www.paraview.org/)** v4.1 for post-processing. 

### Versions: ###

* v0.0.1 Start JosePH
* v0.1.1 MoFEM is stabile enough to be used by other outside group. It is working version with Linear Elasticity, Non Linear Elasticity, Damage at interface element, Brittle Crack propagation, Lienar Dynamics, Potential Flow, Mesh Smoothing and Potential Flow.
* v0.1.2 Solves several bugs. Higher order approximation of gemetry, projection from 10 Nodet Tetrahederal from/to hierarchical approximation basis. New atom tests and benchmarks and more.
* v0.1.3 Face splitting for crack propagation. Several other changes.
* v0.1.4 Add new mechanism for implementation of finite elements. Look to thermal examples to see details. Moreover efficiency improvments, consistent handling of dirihlet and numan boundary conditions.

[Go to Home Page](https://bitbucket.org/likask/mofem-joseph/wiki/Home)
