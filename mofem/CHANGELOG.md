### MoFEM v0.13.0

- Improvements in mesh refinement
- Storing parent ent in multi-index (makes initialisation data structures faster(
- Fixes on side loops on the skeleton
- Adding child and parent loops
- Fixes testing postprocessing
- Fixes in bit ref level in finite elements and problem dofs
- Add user data operator hooks
- Combine FieldBlas test in one programme. 
- Fix getting 2d Dg tensor from matrix
- Refactor FieldBals and add new functionality
- Improvements in getFTensor1DiffN
- Add variants for TS solver, in particular IMEX to pipelines.
- Improve getting agencies with on building elements level with vector not range. 
- Fixing parallel pertaining and scalability
- Remove obsolete code and other fixes.
- Other changes and fixes

### MoFEM v0.12.0

- Higher order geometry operator for volume
- New FTensor functions and operators
- Extensions of Form integrators
- Newton-Cotes integration rule for tets
- MoFEM webpage major update
- Fixes for installation pages and installations scripts
- Add L2 base for AINSWORTH_BERNSTEIN_BEZIER_BASE on tets
- Generalise approximation test to work on faces and volumes
- Add series of tests for approximation on tets for H1 and L2 spaces
- Improvements in MPI communicators
- Fixes for PETSc 3.14.0
- Improvements in error messaging
- Other fixes

----

### MoFEM v0.11.0

- Added functionality for matrix functions and their derivatives
- Fixes of form integrators
- Added new tensor operators
- Novel functionality of marking DOFs for boundary conditions
- Extended handling of node data from `med` and `rmed` files
- Implementation of quad elements
- Extension of Docker functionality
- Installation with new version of `spack`
- Support for `Jupyter` notebooks with `Docker`
- Various minor fixes and improvements

----

### MoFEM v0.10.0

- Improved memory and run-time efficiency
- Fixes and improvements in the mesh-cutting algorithm
- New method of enforcing boundary conditions
- Refactoring of the implementation of multi-indices
- Refactoring of logging
- Changes on web-pages
- New tutorials and changes for tutorials web-pages 
- Various minor bug fixes

----

### MoFEM v0.9.2

- New logging interface based on boost.log
- Fixes in PrismInterface
- Developments of mesh cutting algorithm
- Interface for side volume elements for contact prism elements
- H-div for contact elements
- Surface pressure ALE
- Contact element ALE
- Gauss point convection for contact (moderate slip)
- ALM frictionless contact
- Computation of the real contact area and active set of gauss points
- Rotational Dirichlet BCs
- 8 new lessons
- New tutorials
- Code refactoring and minor fixes