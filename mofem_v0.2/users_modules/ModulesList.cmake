# Obsolete, i.e. implementation is obsolete and need to be changed or replaced
# by alternative classes or functions. Use of oboslete implementation should be
# avoided. This code will be removed in futuer versions of MoFEM

include_directories(${UM_SOURCE_DIR}/obsolete)
include_directories(${UM_SOURCE_DIR}/obsolete/c)
include_directories(${UM_SOURCE_DIR}/obsolete/c_impl)

add_library(users_modules_obsolete
  ${UM_SOURCE_DIR}/obsolete/impl/FEMethod_SurfaceConstrains.cpp
  ${UM_SOURCE_DIR}/obsolete/impl/FEMethod_LowLevelStudent.cpp
  ${UM_SOURCE_DIR}/obsolete/impl/FEMethod_UpLevelStudent.cpp
)

add_library(complex_for_lazy_obsolete
  ${UM_SOURCE_DIR}/obsolete/c_impl/complex_for_lazy.c
  ${UM_SOURCE_DIR}/obsolete/impl/FEMethod_ComplexForLazy.cpp
  ${UM_SOURCE_DIR}/obsolete/impl/SurfacePressureComplexForLazy.cpp
)

# List of sources for user_modules libaries. This list
# is set by users in user modules subdirectories
set(UM_LIB_SOURCES "")

include(analytical_dirichlet_boundary_conditions/CMakeLists.txt)

# Users modules library, common for all programs
add_library(users_modules ${UM_LIB_SOURCES})

# Atom test user modules
add_subdirectory(atom_tests)
add_subdirectory(
  analytical_dirichlet_boundary_conditions/atom_tests
)

# Users modules
add_subdirectory(thermal)
add_subdirectory(ultraweak)
add_subdirectory(convective_mass_element)
add_subdirectory(ground_surface_temperature)
add_subdirectory(elasticity)
add_subdirectory(cohesive_interface)
add_subdirectory(nonlinear_elasticity)
add_subdirectory(fracture_mechanics)
add_subdirectory(helmholtz)
add_subdirectory(homogenisation)
add_subdirectory(moisture_transport)
