/** \file magnetostatic.cpp
  * \ingroup maxwell_element
  *
  * \brief Example implementation of magnetostaic problem
  *
  */

/* This file is part of MoFEM.
 * MoFEM is free software: you can redistribute it and/or modify it under
 * the terms of the GNU Lesser General Public License as published by the
 * Free Software Foundation, either version 3 of the License, or (at your
 * option) any later version.
 *
 * MoFEM is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
 * License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with MoFEM. If not, see <http://www.gnu.org/licenses/>. */

#include <BasicFiniteElements.hpp>
#include <MagneticElement.hpp>
using namespace MoFEM;

static char help[] = "...\n\n";

int main(int argc, char *argv[]) {

  
  

  PetscInitialize(&argc,&argv,(char *)0,help);

  moab::Core mb_instance;
  moab::Interface& moab = mb_instance;

  ParallelComm* pcomm = ParallelComm::get_pcomm(&moab,MYPCOMM_INDEX);
  if(pcomm == NULL) pcomm =  new ParallelComm(&moab,PETSC_COMM_WORLD);

  //Read parameters from line command
  PetscBool flg_file;
  char mesh_file_name[255];
  PetscInt order = 2;
  PetscBool is_partitioned = PETSC_FALSE;

  ierr = PetscOptionsBegin(
    PETSC_COMM_WORLD,"","Shell prism configure","none"
  ); CHKERRG(ierr);
  ierr = PetscOptionsString(
    "-my_file",
    "mesh file name","", "mesh.h5m",mesh_file_name, 255, &flg_file
  ); CHKERRG(ierr);
  ierr = PetscOptionsInt(
    "-my_order",
    "default approximation order",
    "",order,&order,PETSC_NULL
  ); CHKERRG(ierr);
  ierr = PetscOptionsBool(
    "-my_is_partitioned",
    "set if mesh is partitioned (this result that each process keeps only part of the mesh)",
    "",PETSC_FALSE,&is_partitioned,PETSC_NULL
  ); CHKERRG(ierr);
  ierr = PetscOptionsEnd(); CHKERRG(ierr);

  if(is_partitioned == PETSC_TRUE) {
    //Read mesh to MOAB
    const char *option;
    option = "PARALLEL=READ_PART;"
      "PARALLEL_RESOLVE_SHARED_ENTS;"
      "PARTITION=PARALLEL_PARTITION;";
    rval = moab.load_file(mesh_file_name, 0, option); CHKERRG(rval);
  } else {
    const char *option;
    option = "";
    rval = moab.load_file(mesh_file_name, 0, option); CHKERRG(rval);
  }

  //Create mofem interface
  MoFEM::Core core(moab);
  MoFEM::Interface& m_field = core;

  MagneticElement magnetic(m_field);
  magnetic.blockData.oRder = order;
  ierr = magnetic.getNaturalBc(); CHKERRG(ierr);
  ierr = magnetic.getEssentialBc(); CHKERRG(ierr);
  ierr = magnetic.createFields(); CHKERRG(ierr);
  ierr = magnetic.createElements(); CHKERRG(ierr);
  ierr = magnetic.createProblem(); CHKERRG(ierr);
  ierr = magnetic.solveProblem(); CHKERRG(ierr);
  ierr = magnetic.postProcessResults(); CHKERRG(ierr);
  ierr = magnetic.destroyProblem(); CHKERRG(ierr);

  ierr = PetscFinalize(); CHKERRG(ierr);

  return 0;
}
