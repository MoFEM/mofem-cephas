/** \file dm_mofem.cpp

  \brief Atom test for Data Manager Interface in MoFEM

*/

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
    if (flg != PETSC_TRUE) {
      SETERRQ(PETSC_COMM_SELF, 1, "*** ERROR -my_file (MESH FILE NEEDED)");
    }
    const char *option;
    option = "";

    // register new dm type, i.e. mofem
    DMType dm_name = "DMMOFEM";
    CHKERR DMRegister_MoFEM(dm_name);

    // Create dm instance
    DM dm;
    CHKERR DMCreate(PETSC_COMM_WORLD, &dm);
    CHKERR DMSetType(dm, dm_name);

    // read mesh and create moab and mofem data structures
    moab::Core mb_instance;
    moab::Interface &moab = mb_instance;
    CHKERR moab.load_file(mesh_file_name, 0, option);
    MoFEM::Core core(moab);
    MoFEM::Interface &m_field = core;

    EntityHandle root_set = moab.get_root_set();
    // add all entities to database, all of them will be used
    BitRefLevel bit_level0;
    bit_level0.set(0);
    CHKERR m_field.getInterface<BitRefManager>()->setBitRefLevelByDim(
        root_set, 3, bit_level0);
    // define & build field
    const int field_rank = 1;
    CHKERR m_field.add_field("FIELD", H1, AINSWORTH_LEGENDRE_BASE, field_rank);
    // add entities to field
    CHKERR m_field.add_ents_to_field_by_type(root_set, MBTET, "FIELD");
    // set app. order
    int order = 4;
    CHKERR m_field.set_field_order(root_set, MBTET, "FIELD", order);
    CHKERR m_field.set_field_order(root_set, MBTRI, "FIELD", order);
    CHKERR m_field.set_field_order(root_set, MBEDGE, "FIELD", order);
    CHKERR m_field.set_field_order(root_set, MBVERTEX, "FIELD", 1);
    // build data structures for fields
    CHKERR m_field.build_fields();

    // define & build finite elements
    CHKERR m_field.add_finite_element("FE");
    // Define rows/cols and element data
    CHKERR m_field.modify_finite_element_add_field_row("FE", "FIELD");
    CHKERR m_field.modify_finite_element_add_field_col("FE", "FIELD");
    CHKERR m_field.modify_finite_element_add_field_data("FE", "FIELD");
    // add entities to finite element
    CHKERR m_field.add_ents_to_finite_element_by_type(root_set, MBTET, "FE");
    // build finite elemnts
    CHKERR m_field.build_finite_elements();
    // build adjacencies
    CHKERR m_field.build_adjacencies(bit_level0);

    // set dm data structure which created mofem data structures
    CHKERR DMMoFEMCreateMoFEM(dm, &m_field, dm_name, bit_level0);
    CHKERR DMSetFromOptions(dm);
    CHKERR DMMoFEMAddElement(dm, "FE");
    CHKERR DMSetUp(dm);

    Mat m;
    Vec l, g;

    CHKERR DMCreateGlobalVector(dm, &g);
    CHKERR DMCreateLocalVector(dm, &l);
    CHKERR DMCreateMatrix(dm, &m);

    // glob loc
    CHKERR VecSet(g, 1.1);
    CHKERR VecGhostUpdateBegin(g, INSERT_VALUES, SCATTER_FORWARD);
    CHKERR VecGhostUpdateEnd(g, INSERT_VALUES, SCATTER_FORWARD);

    CHKERR DMGlobalToLocalBegin(dm, g, ADD_VALUES, l);
    CHKERR DMGlobalToLocalEnd(dm, g, ADD_VALUES, l);

    // loc glob
    CHKERR DMLocalToGlobalBegin(dm, l, ADD_VALUES, g);
    CHKERR DMLocalToGlobalEnd(dm, l, ADD_VALUES, g);

    PetscViewer viewer;
    CHKERR PetscViewerASCIIOpen(PETSC_COMM_WORLD, "dm_mofem.txt", &viewer);
    const double chop = 1e-4;
    CHKERR VecChop(g, chop);
    VecView(g, viewer);
    CHKERR PetscViewerDestroy(&viewer);

    CHKERR VecDestroy(&g);
    CHKERR VecDestroy(&l);
    CHKERR MatDestroy(&m);

    // destroy dm
    CHKERR DMDestroy(&dm);
  }
  CATCH_ERRORS;

  // finish work cleaning memory, getting statistics, ect.
  CHKERR MoFEM::Core::Finalize();

  return 0;
}
