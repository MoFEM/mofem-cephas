# Users modules

add_subdirectory(thermal)
add_subdirectory(ultraweak)

# Obsolete

include_directories(${MoFEM_PROJECT_SOURCE_DIR}/users_modules/obsolete)
include_directories(${MoFEM_PROJECT_SOURCE_DIR}/users_modules/obsolete/c)

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

add_subdirectory(elasticity)
add_subdirectory(arc_length_nonlinear_elasticity)
add_subdirectory(fracture_mechanics)
add_subdirectory(homogenisation)

