/*! \page coding_practice Coding practice


\section mofam_and_user_modules MoFEM and User Modules

All new finite elements or modification or finite elements are implemented in
users modules (UM). Only "reliable" and "proven" to work finite elements are
moved to MoFEM finite element library. Implementation to MoFEM library (not UM)
need to be always discussed, agreed and planed on CMatGU
<cmatgu@googlegroups.com>.  Work on UM is less controlled and contained to user
directory. It is planed that it will be single version of MoFEM library but
many UM forks. Some UM can be in in depended from MoFEM library repository.

 In default UM are build as a part of MoFEM
library, i.e. source files from MoFEM source directory are used in compilation
process. UM can be created as "stand alone" source, when out of source build is
created as follows
\code
cmake -DSTAND_ALLONE_USERS_MODULES ../user_modules
\endcode 

\section user_module Adding user module

In each module you can have two type directory data, 
  -# Simple
  -# Extended

User module is added to \em ModulesLists.cmake file using cmake command:
\code
add_subdirectory(my_new_module)
\endcode
Simple directory structure consist no subdirectories. f.e. elasticity. Extended (recommended)
data structure consist subdirectories, f.e. homogenisation, and follows pattern

\code
-> /atom_tests
-> /src <- hpp files
-> /src/impl <- cpp files form library
-> /meshes
-> /data
-> /doc
\endcode

Not all elements of module source tree are compulsory, however each user module
and new MoFEM functionality should have associated \em atom \em test verifying
implementation and each module should have README file or module documentation
in doc using Doxygen.

\section Source Code Style and Best Practices

MoFEM code should follow MoAB code style and best pratices listed here
<http://www.mcs.anl.gov/~fathom/moab-docs/html/styleguide.html>.

- Names:
  - Class names should be in the CamelBack style, e.g. EdgeMesh or VertexMesher.

  - Class member variables should be camelBack, e.g. EdgeMesh::schemeType; each
    member variable, e.g. int memberVariable, should have set/get functions
    void member_variable(int newval) and int member_variable(), respectively

  - Enumeration values should be all captitalized, with underscores avoided if
    possible (the enumeration name indicates the general purpose of the
    enumeration, so e.g. we use EQUAL, not EQUAL_MESH)

  - Each class header should be fully commented; that includes:
  
  - A \\file comment block at the top of the file; DO NOT include things like
    Author and Date blocks; this stuff is available from subversion if we
    really need to know.

  - Each function in both the public and private interfaces should be
    commented, INCLUDING ANY ARGUMENTS AND RETURN VALUES. See the MOAB classes
    for examples of how to format these comments. 

  - Developers should avoid using #include in header files, as they propagate
    dependencies more widely than necessary. 

  - Local variables and function argument have names with small letters, i.e.
    temperature_val, nodal_position. 


\section Making Repository Commits

As a general rule, developers should update frequently, and commit changes
often. However, the repository should always remain in a state where the code
can be compiled. Most of the time, the code should also successfully execute
"ctest" run from the top-level directory. If you commit code that violates this
principal, it should be your first priority to return the repository code to a
compilable state, and your second priority to make sure "ctest" runs without
errors.

Finished proration of work should be committed by pull-request to CDashTesting
branch. Commits to the CDashTesting branch should also come with a non-trivial,
useful, non-verbose log message. It is required that before pull request user
merge current CDashTesting branch and verify if all ctest run without fails. If
pull request is accepted it is user responsibility to verify results on
CDash server <http://cdash.eng.gla.ac.uk/cdash/>. The first priority will be to
eliminate compilation errors, completion warnings, failed tests and memory
leaks.

*/
