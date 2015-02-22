# Users modules

# Users modules library, common for all programs
#add_library(users_modules
#)

add_subdirectory(thermal)
add_subdirectory(ultraweak)
add_subdirectory(convective_mass_element)
add_subdirectory(ground_surface_temperature)
add_subdirectory(field_approximation)
add_subdirectory(analytical_dirihlet_boundary_conditions)

# Obsolete, i.e. implementation is obsolete and need to be changed or replaced
# by alternative classes or functions

include_directories(${UM_SOURCE_DIR}/obsolete)
include_directories(${UM_SOURCE_DIR}/obsolete/c)
include_directories(${UM_SOURCE_DIR}/obsolete/c_impl)

add_library(users_modules_obsolete
  ${UM_SOURCE_DIR}/obsolete/impl/ArcLengthTools.cpp
  ${UM_SOURCE_DIR}/obsolete/impl/FEMethod_SurfaceConstrains.cpp
  ${UM_SOURCE_DIR}/obsolete/impl/FEMethod_LowLevelStudent.cpp
  ${UM_SOURCE_DIR}/obsolete/impl/FEMethod_UpLevelStudent.cpp
)

add_library(complex_for_lazy_obsolete
  ${UM_SOURCE_DIR}/obsolete/c_impl/complex_for_lazy.c
  ${UM_SOURCE_DIR}/obsolete/impl/FEMethod_ComplexForLazy.cpp
  ${UM_SOURCE_DIR}/obsolete/impl/SurfacePressureComplexForLazy.cpp
)

add_subdirectory(elasticity)
add_subdirectory(cohesive_interface)
add_subdirectory(nonlinear_elasticity)
add_subdirectory(fracture_mechanics)
add_subdirectory(homogenisation)
add_subdirectory(moisture_transport)

# Atom tests for users modules

add_subdirectory(atom_tests)

