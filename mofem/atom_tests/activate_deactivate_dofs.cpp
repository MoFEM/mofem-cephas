/** \file build_problems.cpp

  \brief Atom test for building composite problem

  \bug Not verifying what if partitioned mesh is loaded.

*/

#include <MoFEM.hpp>

using namespace MoFEM;

static char help[] = "...\n\n";

int main(int argc, char *argv[]) {
  MoFEM::Core::Initialize(&argc, &argv, (char *)0, help);

  try {

    moab::Core mb_instance;
    moab::Interface &moab = mb_instance;
    int rank;
    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

    PetscBool flg = PETSC_TRUE;
    char mesh_file_name[255];
#if PETSC_VERSION_GE(3, 6, 4)
    CHKERR PetscOptionsGetString(PETSC_NULL, "", "-my_file", mesh_file_name,
                                 255, &flg);
#else
    CHKERR PetscOptionsGetString(PETSC_NULL, PETSC_NULL, "-my_file",
                                 mesh_file_name, 255, &flg);
#endif
    if (flg != PETSC_TRUE) {
      SETERRQ(PETSC_COMM_SELF, 1, "*** ERROR -my_file (MESH FILE NEEDED)");
    }

    const char *option;
    option = "";
    CHKERR moab.load_file(mesh_file_name, 0, option);
    CHKERRG(rval);

    // Create MoFEM (Joseph) database
    MoFEM::Core core(moab);
    MoFEM::Interface &m_field = core;

    // set entitities bit level
    BitRefLevel bit_level0;
    bit_level0.set(0);
    EntityHandle meshset_level0;
    CHKERR moab.create_meshset(MESHSET_SET, meshset_level0);
    CHKERRG(rval);
    CHKERR m_field.getInterface<BitRefManager>()->setBitRefLevelByDim(
        0, 3, bit_level0);

    // Fields
    CHKERR m_field.add_field("F1", L2, AINSWORTH_LEGENDRE_BASE, 1);
    CHKERR m_field.add_field("F2", HDIV, AINSWORTH_LEGENDRE_BASE, 1);

    // meshset consisting all entities in mesh
    EntityHandle root_set = moab.get_root_set();
    // add entities to field
    CHKERR m_field.add_ents_to_field_by_type(root_set, MBTET, "F1");
    CHKERR m_field.add_ents_to_field_by_type(root_set, MBTET, "F2");

    // set app. order
    // see Hierarchic Finite Element Bases on Unstructured Tetrahedral Meshes
    // (Mark Ainsworth & Joe Coyle)
    int order = 4;
    CHKERR m_field.set_field_order(root_set, MBTET, "F1", order - 1, 2);
    CHKERR m_field.set_field_order(root_set, MBTET, "F2", order, 2);
    CHKERR m_field.set_field_order(root_set, MBTRI, "F2", order, 2);

    CHKERR m_field.build_fields();

    order = 2;
    CHKERR m_field.set_field_order(root_set, MBTET, "F1", order - 1, 2);
    CHKERR m_field.set_field_order(root_set, MBTET, "F2", order, 2);
    CHKERR m_field.set_field_order(root_set, MBTRI, "F2", order, 2);

    CHKERR m_field.clear_inactive_dofs();
    CHKERR m_field.build_fields();

    auto dofs_ptr = m_field.get_dofs();
    const unsigned int expected_size = 300;
    if (dofs_ptr->size() != expected_size) {
      SETERRQ2(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
               "Data inconsistency, should %d dofs but is %d", expected_size,
               dofs_ptr->size());
    }
  }
  CATCH_ERRORS;

  // finish work cleaning memory, getting statistics, ect.
  CHKERR MoFEM::Core::Finalize();

  return 0;
}
