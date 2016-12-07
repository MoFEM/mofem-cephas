/** \file reading_med.cpp

  \brief Testing interface for reading and writing med files

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

#include <MoFEM.hpp>
using namespace MoFEM;

static char help[] = "...\n\n";

int main(int argc, char *argv[]) {

  ErrorCode rval;
  PetscErrorCode ierr;

  PetscInitialize(&argc,&argv,(char *)0,help);

  moab::Core mb_instance;
  moab::Interface& moab = mb_instance;
  ParallelComm* pcomm = ParallelComm::get_pcomm(&moab,MYPCOMM_INDEX);
  if(pcomm == NULL) pcomm =  new ParallelComm(&moab,PETSC_COMM_WORLD);

  //Create MoFEM (Joseph) database
  MoFEM::Core core(moab);
  MoFEM::Interface& m_field = core;

  MedInterface *med_interface_ptr;
  ierr = m_field.query_interface(med_interface_ptr); CHKERRQ(ierr);
  ierr = med_interface_ptr->readMed(); CHKERRQ(ierr);
  // ierr = med_interface_ptr->medGetFieldNames(); CHKERRQ(ierr);

  MeshsetsManager *meshsets_interface_ptr;
  ierr = m_field.query_interface(meshsets_interface_ptr); CHKERRQ(ierr);
  ierr = meshsets_interface_ptr->setMeshsetFromFile(); CHKERRQ(ierr);

  for(
    CubitMeshSet_multiIndex::iterator cit = meshsets_interface_ptr->getBegin();
    cit!=meshsets_interface_ptr->getEnd(); cit++
  ) {
    std::cout << *cit << endl;
  }

  rval = moab.write_file("out.h5m"); CHKERRQ_MOAB(rval);

  ierr = PetscFinalize(); CHKERRQ(ierr);

  return 0;
}
