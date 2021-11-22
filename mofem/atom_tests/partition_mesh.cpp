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

  // initialize petsc
  MoFEM::Core::Initialize(&argc, &argv, (char *)0, help);

  try {

    PetscBool flg = PETSC_TRUE;
    char mesh_file_name[255];
#if PETSC_VERSION_GE(3, 6, 4)
    CHKERR PetscOptionsGetString(PETSC_NULL, "", "-my_file", mesh_file_name,
                                 255, &flg);
#else
    CHKERR PetscOptionsGetString(PETSC_NULL, PETSC_NULL, "-my_file",
                                 mesh_file_name, 255, &flg);
#endif
    if (flg != PETSC_TRUE) 
      SETERRQ(PETSC_COMM_SELF, 1, "*** ERROR -my_file (MESH FILE NEEDED)");

    // read mesh and create moab and mofem datastrutures

    moab::Core mb_instance;
    moab::Interface &moab = mb_instance;
    {
      const char *option;
      option = "";
      CHKERR moab.load_file(mesh_file_name, 0, option);
    }
    MoFEM::Core core(moab);
    MoFEM::Interface &m_field = core;

    EntityHandle root_set = moab.get_root_set();
    Range tets;
    moab.get_entities_by_type(root_set, MBTET, tets, false);

    Tag th_vertex_weight;
    int def_val = 1;
    CHKERR moab.tag_get_handle("VERTEX_WEIGHT", 1, MB_TYPE_INTEGER,
                               th_vertex_weight, MB_TAG_CREAT | MB_TAG_DENSE,
                               &def_val);

    ProblemsManager *prb_mng_ptr;
    CHKERR m_field.getInterface(prb_mng_ptr);
    CHKERR prb_mng_ptr->partitionMesh(tets, 3, 2, 2, &th_vertex_weight, NULL,
                                      NULL, VERBOSE, false);

    EntityHandle meshset;
    CHKERR moab.create_meshset(MESHSET_SET, meshset);
    CHKERR moab.add_entities(meshset, tets);
    if (!m_field.get_comm_rank()) {
      CHKERR moab.write_file("partitioned_mesh.h5m");
    }
  }
  CATCH_ERRORS;

  PetscBarrier(PETSC_NULL);

  moab::Core mb_instance2;
  moab::Interface &moab2 = mb_instance2;
  {
    const char *option = "DEBUG_IO;"
                         "PARALLEL=BCAST_DELETE;"
                         "PARALLEL_RESOLVE_SHARED_ENTS;"
                         "PARTITION=PARALLEL_PARTITION;";
    CHKERR moab2.load_file("partitioned_mesh.h5m", 0, option);
  }

  CHKERR MoFEM::Core::Finalize();

  return 0;
}
