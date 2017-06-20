/** \file simple_interface.cpp
  * \example simple_interface.hpp
  * \brief Simple interface
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

  // initialize petsc
  PetscInitialize(&argc,&argv,(char *)0,help);

  try {

    // Create MoAB database
    moab::Core moab_core;
    moab::Interface& moab = moab_core;

    // Create MoFEM database and link it to MoAB
    MoFEM::Core mofem_core(moab);
    MoFEM::Interface& m_field = mofem_core;

    // Register DM Manager
    DMType dm_name = "DMMOFEM";
    ierr = DMRegister_MoFEM(dm_name); CHKERRQ(ierr);

    // Simple interface
    Simple *simple_interface;
    ierr = m_field.query_interface(simple_interface); CHKERRQ(ierr);
    {
      ierr = simple_interface->getOptions(); CHKERRQ(ierr);
      ierr = simple_interface->loadFile(); CHKERRQ(ierr);
      ierr = simple_interface->addDomainField("U",H1,AINSWORTH_LEGENDRE_BASE,1); CHKERRQ(ierr);
      ierr = simple_interface->addBoundaryField("L",H1,AINSWORTH_LEGENDRE_BASE,1); CHKERRQ(ierr);
      ierr = simple_interface->setFieldOrder("U",2); CHKERRQ(ierr);
      ierr = simple_interface->setFieldOrder("L",1); CHKERRQ(ierr);
      ierr = simple_interface->setUp(); CHKERRQ(ierr);
    }

  } catch (MoFEMException const &e) {
    SETERRQ(PETSC_COMM_SELF,e.errorCode,e.errorMessage);
  }

  // finish work cleaning memory, getting statistics, etc.
  ierr = PetscFinalize(); CHKERRQ(ierr);

  return 0;
}
