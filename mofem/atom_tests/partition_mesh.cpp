/** \file partition_mesh.cpp
  \example partition_mesh.cpp 
  \brief Atom testing mesh partitioning

*/

/* MIT License
 *
 * Copyright (c) 2022
 * 
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 * 
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 * 
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */

#include <MoFEM.hpp>

using namespace MoFEM;

static char help[] = "...\n\n";

int main(int argc, char *argv[]) {

  // initialize petsc
  MoFEM::Core::Initialize(&argc, &argv, (char *)0, help);

  int nb_vertices;

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

    CHKERR moab.get_number_entities_by_dimension(0, 0, nb_vertices, true);

    EntityHandle root_set = moab.get_root_set();
    Range tets;
    moab.get_entities_by_type(root_set, MBTET, tets, false);

    Tag th_vertex_weight;
    int def_val = 1;
    CHKERR moab.tag_get_handle("VERTEX_WEIGHT", 1, MB_TYPE_INTEGER,
                               th_vertex_weight, MB_TAG_CREAT | MB_TAG_DENSE,
                               &def_val);

    CommInterface *comm_interafce_ptr = m_field.getInterface<CommInterface>();
    CHKERR comm_interafce_ptr->partitionMesh(
        tets, 3, 2, m_field.get_comm_size(), &th_vertex_weight, NULL, NULL,
        VERBOSE, false);

    if (!m_field.get_comm_rank()) {
      CHKERR moab.write_file("partitioned_mesh.h5m");
    }
  }
  CATCH_ERRORS;

  PetscBarrier(PETSC_NULL);

  try {

    moab::Core mb_instance2;
    moab::Interface &moab2 = mb_instance2;


    MoFEM::CoreTmp<1> core(moab2);
    MoFEM::Interface &m_field = core;

    // Register DM Manager
    DMType dm_name = "DMMOFEM";
    CHKERR DMRegister_MoFEM(dm_name);

    // Test build simple problem
    const char *option = "DEBUG_IO;"
                         "PARALLEL=READ_PART;"
                         "PARALLEL_RESOLVE_SHARED_ENTS;"
                         "PARTITION=PARALLEL_PARTITION;";

    CHKERR m_field.getInterface<Simple>()->getOptions();
    CHKERR m_field.getInterface<Simple>()->loadFile(option,
                                                    "partitioned_mesh.h5m");
    CHKERR m_field.getInterface<Simple>()->addDomainField(
        "U", H1, AINSWORTH_LEGENDRE_BASE, 1);
    CHKERR m_field.getInterface<Simple>()->setFieldOrder("U", 1);
    CHKERR m_field.getInterface<Simple>()->setUp();

    auto dm  = m_field.getInterface<Simple>()->getDM();

    const MoFEM::Problem *problem_ptr;
    CHKERR DMMoFEMGetProblemPtr(dm, &problem_ptr);

    if(problem_ptr->nbDofsRow != nb_vertices)
      SETERRQ(PETSC_COMM_WORLD, MOFEM_ATOM_TEST_INVALID,
              "Number of vertices and DOFs is inconstent");

    MOFEM_LOG_CHANNEL("WORLD");
    MOFEM_TAG_AND_LOG("WORLD", Sev::inform, "Atom test")
        << "All is good in this test";
  }
  CATCH_ERRORS;

  CHKERR MoFEM::Core::Finalize();

  return 0;
}
