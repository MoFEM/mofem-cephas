/** \file partition_mesh.cpp

  \brief Atom testing mesh partitioning

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

  //initialize petsc
  PetscInitialize(&argc,&argv,(char *)0,help);

  {

    PetscBool flg = PETSC_TRUE;
    char mesh_file_name[255];
    #if PETSC_VERSION_GE(3,6,4)
    ierr = PetscOptionsGetString(PETSC_NULL,"","-my_file",mesh_file_name,255,&flg); CHKERRQ(ierr);
    #else
    ierr = PetscOptionsGetString(PETSC_NULL,PETSC_NULL,"-my_file",mesh_file_name,255,&flg); CHKERRQ(ierr);
    #endif
    if(flg != PETSC_TRUE) {
      SETERRQ(PETSC_COMM_SELF,1,"*** ERROR -my_file (MESH FILE NEEDED)");
    }

    //read mesh and create moab and mofem datastrutures

    moab::Core mb_instance;
    moab::Interface& moab = mb_instance;
    {
      const char *option;
      option = "";
      // option = "PARALLEL=BCAST_DELETE;"
      // "PARALLEL_RESOLVE_SHARED_ENTS;"
      // "PARTITION=PARALLEL_PARTITION;";

      rval = moab.load_file(mesh_file_name, 0, option); CHKERRQ_MOAB(rval);
    }
    MoFEM::Core core(moab);
    MoFEM::Interface& m_field = core;

    EntityHandle root_set = moab.get_root_set();
    Range tets;
    moab.get_entities_by_type(root_set,MBTET,tets,false);
    ierr = m_field.partition_mesh(tets,3,2,2); CHKERRQ(ierr);

    EntityHandle meshset;
    rval = moab.create_meshset(MESHSET_SET,meshset); CHKERRQ_MOAB(rval);
    rval = moab.add_entities(meshset,tets); CHKERRQ_MOAB(rval);
    // // resolve shared entities
    // ParallelComm* pcomm = ParallelComm::get_pcomm(&moab,MYPCOMM_INDEX);
    if(!m_field.getCommRank()) {
      // rval = moab.write_file("partitioned_mesh.h5m","MOAB","",&meshset,1); CHKERRQ_MOAB(rval);
      rval = moab.write_file("partitioned_mesh.h5m"); CHKERRQ_MOAB(rval);
      // rval = moab.write_file("partitioned_mesh.h5m"); CHKERRQ_MOAB(rval);
    }

  }

  PetscBarrier(PETSC_NULL);

  moab::Core mb_instance2;
  moab::Interface& moab2 = mb_instance2;
  {
    const char *option =
    "DEBUG_IO;"
    "PARALLEL=BCAST_DELETE;"
    "PARALLEL_RESOLVE_SHARED_ENTS;"
    "PARTITION=PARALLEL_PARTITION;";
    rval = moab2.load_file("partitioned_mesh.h5m",0,option); CHKERRQ_MOAB(rval);
  }

  ierr = PetscFinalize(); CHKERRQ(ierr);

  return 0;

}
