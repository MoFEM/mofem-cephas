# Users modules

add_subdirectory(thermal)
add_subdirectory(ultraweak)
add_subdirectory(convective_mass_element)
add_subdirectory(ground_surface_temperature)
add_subdirectory(field_approximation)

# Obsolete

include_directories(${UM_SOURCE_DIR}/obsolete)
include_directories(${UM_SOURCE_DIR}/obsolete/c)
include_directories(${UM_SOURCE_DIR}/obsolete/c_impl)

add_library(users_modules_obsolete
  ${UM_SOURCE_DIR}/obsolete/impl/ArcLengthTools.cpp
  ${UM_SOURCE_DIR}/obsolete/impl/FEMethod_SurfaceConstrains.cpp
  ${UM_SOURCE_DIR}/obsolete/impl/FEMethod_LowLevelStudent.cpp
  ${UM_SOURCE_DIR}/obsolete/impl/FEMethod_UpLevelStudent.cpp
  ${UM_SOURCE_DIR}/obsolete/impl/MatShellConstrainsByMarkAinsworth.cpp
)

add_library(complex_for_lazy_obsolete
  ${UM_SOURCE_DIR}/obsolete/c_impl/complex_for_lazy.c
  ${UM_SOURCE_DIR}/obsolete/impl/FEMethod_ComplexForLazy.cpp
  ${UM_SOURCE_DIR}/obsolete/impl/SurfacePressureComplexForLazy.cpp
)

add_subdirectory(analytical_dirihlet_boundary_conditions)
add_subdirectory(elasticity)
add_subdirectory(cohesive_interface)
add_subdirectory(nonlinear_elasticity)
add_subdirectory(fracture_mechanics)

add_subdirectory(helmholtz)

add_subdirectory(homogenisation)
add_subdirectory(moisture_transport)


#include_directories("${MoFEM_PROJECT_SOURCE_DIR}/users_modules/analytical_dirihlet_boundary_conditions/src")
#add_subdirectory(analytical_dirihlet_boundary_conditions)

#include_directories("${MoFEM_PROJECT_SOURCE_DIR}/users_modules/convective_mass_element/src")
#add_subdirectory(convective_mass_element)

add_subdirectory(atom_tests)

