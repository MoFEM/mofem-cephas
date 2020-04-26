/**
 * @file LogManager.cpp
 * @brief Log and register warnings
 * 
 */

namespace MoFEM {

LogManager::LogManager(const MoFEM::Core &core) : cOre(core) {}

MoFEMErrorCode LogManager::getSubInterfaceOptions() { return getOptions(); }

MoFEMErrorCode LogManager::getOptions() {
  MoFEMFunctionBegin;
  CHKERR PetscOptionsBegin(PETSC_COMM_WORLD, "waning_",
                           "Warning interface options", "none");

  defaultVerbosityLevel = VERBOSE;
  CHKERR PetscOptionsInt("-level", "verbosity level", "", defaultVerbosityLevel,
                         &defaultVerbosityLevel, PETSC_NULL);

  ierr = PetscOptionsEnd();
  CHKERRG(ierr);
  MoFEMFunctionReturn(0);
}


} // MOFEM namespace