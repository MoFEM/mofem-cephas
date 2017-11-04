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




  //initialize petsc
  PetscInitialize(&argc,&argv,(char *)0,help);

  try {

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

      rval = moab.load_file(mesh_file_name, 0, option); CHKERRG(rval);
    }
    MoFEM::Core core(moab);
    MoFEM::Interface& m_field = core;

    EntityHandle root_set = moab.get_root_set();
    Range tets;
    moab.get_entities_by_type(root_set,MBTET,tets,false);

    Tag th_vertex_weight;
    int def_val = 1;
    rval = moab.tag_get_handle(
      "VERTEX_WEIGHT",1,MB_TYPE_INTEGER,th_vertex_weight,MB_TAG_CREAT|MB_TAG_DENSE,&def_val
    ); CHKERRQ(ierr);

    ProblemsManager *prb_mng_ptr;
    ierr = m_field.getInterface(prb_mng_ptr); CHKERRQ(ierr);
    ierr = prb_mng_ptr->partitionMesh(tets,3,2,2,&th_vertex_weight,NULL,NULL); CHKERRQ(ierr);

    EntityHandle meshset;
    rval = moab.create_meshset(MESHSET_SET,meshset); CHKERRG(rval);
    rval = moab.add_entities(meshset,tets); CHKERRG(rval);
    // // resolve shared entities
    // ParallelComm* pcomm = ParallelComm::get_pcomm(&moab,MYPCOMM_INDEX);
    if(!m_field.get_comm_rank()) {
      // rval = moab.write_file("partitioned_mesh.h5m","MOAB","",&meshset,1); CHKERRG(rval);
      rval = moab.write_file("partitioned_mesh.h5m"); CHKERRG(rval);
      // rval = moab.write_file("partitioned_mesh.h5m"); CHKERRG(rval);
    }

  } catch (MoFEMException const &e) {
    SETERRQ(PETSC_COMM_SELF,e.errorCode,e.errorMessage);
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
    rval = moab2.load_file("partitioned_mesh.h5m",0,option); CHKERRG(rval);
  }

  ierr = PetscFinalize(); CHKERRQ(ierr);

  return 0;

}
