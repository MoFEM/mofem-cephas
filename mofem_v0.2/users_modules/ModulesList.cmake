# Users modules

add_subdirectory(thermal)
add_subdirectory(ultraweak)

# Obsolete

include_directories(${UM_SOURCE_DIR}/obsolete)
include_directories(${UM_SOURCE_DIR}/obsolete/c)
include_directories(${UM_SOURCE_DIR}/obsolete/c_impl)

add_library(users_modules_obsolete
  ${UM_SOURCE_DIR}/obsolete/impl/ArcLengthTools.cpp
  ${UM_SOURCE_DIR}/obsolete/impl/FEMethod_ComplexForLazy.cpp
  ${UM_SOURCE_DIR}/obsolete/impl/FEMethod_SurfaceConstrains.cpp
  ${UM_SOURCE_DIR}/obsolete/impl/FEMethod_LowLevelStudent.cpp
  ${UM_SOURCE_DIR}/obsolete/impl/FEMethod_UpLevelStudent.cpp
  ${UM_SOURCE_DIR}/obsolete/impl/SurfacePressureComplexForLazy.cpp
  ${UM_SOURCE_DIR}/obsolete/impl/MatShellConstrainsByMarkAinsworth.cpp
)

add_library(complex_for_lazy_obsolete
  ${UM_SOURCE_DIR}/obsolete/c_impl/complex_for_lazy.c
)

add_subdirectory(atom_tests)
add_subdirectory(elasticity)
add_subdirectory(arc_length_nonlinear_elasticity)
add_subdirectory(fracture_mechanics)
add_subdirectory(homogenisation)

include_directories("${MoFEM_PROJECT_SOURCE_DIR}/users_modules/analytical_dirihlet_boundary_conditions/src")
add_subdirectory(analytical_dirihlet_boundary_conditions)

