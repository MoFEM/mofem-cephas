/* Copyright (C) 2013, Lukasz Kaczmarczyk (likask AT wp.pl)
 * --------------------------------------------------------------
 * FIXME: DESCRIPTION
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

#include "FieldInterface.hpp"
#include "FieldCore.hpp"
#include "FEMethod_UpLevelStudent.hpp"
#include "PotentialFlowFEMethod.hpp"
#include "cholesky.hpp"
#include <petscksp.h>

using namespace MoFEM;

ErrorCode rval;
PetscErrorCode ierr;

static char help[] = "...\n\n";


int main(int argc, char *argv[]) {

  try {

  PetscInitialize(&argc,&argv,(char *)0,help);

  Core mb_instance;
  Interface& moab = mb_instance;
  int rank;
  MPI_Comm_rank(PETSC_COMM_WORLD,&rank);
  /*if(rank==0) {
    EntityHandle dummy_meshset;
    rval = moab.create_meshset(MESHSET_SET,dummy_meshset); CHKERR_PETSC(rval);
  }*/

  PetscBool flg = PETSC_TRUE;
  char mesh_file_name[255];
  ierr = PetscOptionsGetString(PETSC_NULL,"-my_file",mesh_file_name,255,&flg); CHKERRQ(ierr);
  if(flg != PETSC_TRUE) {
    SETERRQ(PETSC_COMM_SELF,1,"*** ERROR -my_file (MESH FILE NEEDED)");
  }

  const char *option;
  option = "";//"PARALLEL=BCAST;";//;DEBUG_IO";
  rval = moab.load_file(mesh_file_name, 0, option); CHKERR(rval); 
  ParallelComm* pcomm = ParallelComm::get_pcomm(&moab,MYPCOMM_INDEX);
  if(pcomm == NULL) pcomm =  new ParallelComm(&moab,PETSC_COMM_WORLD);

  FieldCore core(moab);
  FieldInterface& mField = core;

  BitRefLevel bit_level0;
  bit_level0.set(0);

  //build fields
  ierr = mField.build_fields(); CHKERRQ(ierr);
  //build finite elements
  ierr = mField.build_finite_elements(); CHKERRQ(ierr);
  //build adjacencies
  ierr = mField.build_adjacencies(bit_level0); CHKERRQ(ierr);
  //build problem
  ierr = mField.build_problems(); CHKERRQ(ierr);

  //partition problems
  ierr = mField.partition_problem("PRESSURE_PROBLEM"); CHKERRQ(ierr);
  ierr = mField.partition_finite_elements("PRESSURE_PROBLEM"); CHKERRQ(ierr);
  ierr = mField.partition_ghost_dofs("PRESSURE_PROBLEM"); CHKERRQ(ierr);

  Tag th_p;
  double def_val = 0;
  rval = moab.tag_get_handle("P",1,MB_TYPE_DOUBLE,th_p,MB_TAG_CREAT|MB_TAG_SPARSE,&def_val); CHKERR_PETSC(rval);
  for(_IT_GET_DOFS_FIELD_BY_NAME_FOR_LOOP_(mField,"PRESSURE_FIELD",dof)) {
    EntityHandle ent = dof->get_ent();
    double val = dof->get_FieldData();
    rval = moab.tag_set_data(th_p,&ent,1,&val); CHKERR_PETSC(rval);
  }

  SurfaceForces forces_from_pressures(mField);
  Vec total_forces;
  ierr = forces_from_pressures.create_totalForceVector(&total_forces); CHKERRQ(ierr);
  ierr = mField.loop_finite_elements("PRESSURE_PROBLEM","PRESSURE_ELEM",forces_from_pressures);  CHKERRQ(ierr);
  ierr = VecAssemblyBegin(total_forces); CHKERRQ(ierr);
  ierr = VecAssemblyEnd(total_forces); CHKERRQ(ierr);
  ierr = VecView(total_forces,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);
  ierr = VecDestroy(&total_forces); CHKERRQ(ierr);

  PetscFinalize();

  } catch (const char* msg) {
    SETERRQ(PETSC_COMM_SELF,1,msg);
  } catch (const std::exception& ex) {
    ostringstream ss;
    ss << "thorw in method: " << ex.what() << endl;
    SETERRQ(PETSC_COMM_SELF,1,ss.str().c_str());
  }

}
