/** \file dm_build_partitioned_mesh.cpp
  \example dm_partitioned_no_field.cpp
  \brief Testing problem for partitioned mesh with NOFIELD field

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
constexpr bool debug = false;

int main(int argc, char *argv[]) {
  // initialize petsc
  MoFEM::Core::Initialize(&argc, &argv, (char *)0, help);

  try {

    PetscBool flg = PETSC_TRUE;
    char mesh_file_name[255];
    CHKERR PetscOptionsGetString(PETSC_NULL, "", "-my_file", mesh_file_name,
                                 255, &flg);
    if (flg != PETSC_TRUE)
      SETERRQ(PETSC_COMM_SELF, 1, "*** ERROR -my_file (MESH FILE NEEDED)");

    // register new dm type, i.e. mofem
    DMType dm_name = "DMMOFEM";
    CHKERR DMRegister_MoFEM(dm_name);

    // create dm instance
    DM dm;
    CHKERR DMCreate(PETSC_COMM_WORLD, &dm);
    CHKERR DMSetType(dm, dm_name);

    // read mesh and create moab and mofem data structures
    moab::Core mb_instance;
    moab::Interface &moab = mb_instance;

    ParallelComm *pcomm = ParallelComm::get_pcomm(&moab, MYPCOMM_INDEX);
    auto moab_comm_wrap =
        boost::make_shared<WrapMPIComm>(PETSC_COMM_WORLD, false);
    if (pcomm == NULL)
      pcomm =
          new ParallelComm(&moab, moab_comm_wrap->get_comm(), MYPCOMM_INDEX);
          
    const std::string options = "PARALLEL=READ_PART;"
                                "PARALLEL_RESOLVE_SHARED_ENTS;"
                                "PARTITION=PARALLEL_PARTITION;";
    CHKERR moab.load_file(mesh_file_name, 0, options.c_str());

    MoFEM::Core core(moab, PETSC_COMM_WORLD);
    MoFEM::Interface &m_field = core;

    auto *bit_ref_ptr = m_field.getInterface<BitRefManager>();
    auto *comm_interface_ptr = m_field.getInterface<CommInterface>();

    EntityHandle root_set = moab.get_root_set();
    // add all entities to database, all of them will be used
    CHKERR bit_ref_ptr->setBitRefLevelByDim(root_set, 3, BitRefLevel().set(0));

    // add field
    CHKERR m_field.add_field("FIELD", H1, AINSWORTH_LEGENDRE_BASE, 1);
    CHKERR m_field.add_field("LAMBDA", NOFIELD, NOBASE, 3);

    // add entities to field
    CHKERR m_field.add_ents_to_field_by_type(root_set, MBTET, "FIELD");

    // set app. order
    CHKERR m_field.set_field_order(root_set, MBVERTEX, "FIELD", 1);

    // Create vertices for NOFILE
    std::array<double, 6> coords = {0, 0, 0, 0, 0, 0};
    CHKERR m_field.create_vertices_and_add_to_field("LAMBDA", coords.data(), 2);
    CHKERR comm_interface_ptr->makeFieldEntitiesMultishared("LAMBDA", 0,
                                                               NOISY);
    CHKERR bit_ref_ptr->setFieldEntitiesBitRefLevel("LAMBDA",
                                                    BitRefLevel().set(0));

    // build data structures for fields
    CHKERR m_field.build_fields();

    auto print_field_ents = [&](const std::string field_name) {
      auto *field_ents = m_field.get_field_ents();
      auto field_bit_number = m_field.get_field_bit_number(field_name);
      auto lo = field_ents->get<Unique_mi_tag>().lower_bound(
          FieldEntity::getLoBitNumberUId(field_bit_number));
      auto hi = field_ents->get<Unique_mi_tag>().upper_bound(
          FieldEntity::getHiBitNumberUId(field_bit_number));
      for (auto it = lo; it != hi; ++it) {
        std::ostringstream ss;
        ss << "Rank " << m_field.get_comm_rank() << " -> " << **it << endl;
        PetscSynchronizedPrintf(m_field.get_comm(), "%s", ss.str().c_str());
        PetscSynchronizedFlush(m_field.get_comm(), PETSC_STDOUT);
      }
    };

    print_field_ents("LAMBDA");

    if (m_field.get_comm_rank() == 0) {
      CHKERR m_field.getInterface<FieldBlas>()->setField(1, "LAMBDA");
      CHKERR m_field.getInterface<FieldBlas>()->setField(1, "FIELD");
    }

    CHKERR comm_interface_ptr->exchangeFieldData("LAMBDA");
    CHKERR comm_interface_ptr->exchangeFieldData("FIELD");

    auto check_exchanged_values = [&](const std::string field_name) {
      MoFEMFunctionBegin;
      if (m_field.get_comm_rank() != 0) {
        auto *field_ents = m_field.get_field_ents();
        auto field_bit_number = m_field.get_field_bit_number(field_name);
        auto lo = field_ents->get<Unique_mi_tag>().lower_bound(
            FieldEntity::getLoBitNumberUId(field_bit_number));
        auto hi = field_ents->get<Unique_mi_tag>().upper_bound(
            FieldEntity::getHiBitNumberUId(field_bit_number));
        for (auto it = lo; it != hi; ++it) {
          VectorDouble field_data = (*it)->getEntFieldData();
          for (auto v : field_data)
            if (v != 1)
              SETERRQ1(PETSC_COMM_SELF, MOFEM_ATOM_TEST_INVALID,
                       "Wrong value on field %4.3f", v);
        }
      }
      MoFEMFunctionReturn(0);
    };

    CHKERR check_exchanged_values("LAMBDA");

    // define & build finite elements
    CHKERR m_field.add_finite_element("FE");
    // Define rows/cols and element data
    CHKERR m_field.modify_finite_element_add_field_row("FE", "FIELD");
    CHKERR m_field.modify_finite_element_add_field_col("FE", "FIELD");
    CHKERR m_field.modify_finite_element_add_field_data("FE", "FIELD");
    CHKERR m_field.modify_finite_element_add_field_row("FE", "LAMBDA");
    CHKERR m_field.modify_finite_element_add_field_col("FE", "LAMBDA");
    CHKERR m_field.modify_finite_element_add_field_data("FE", "LAMBDA");

    // add entities to finite element
    CHKERR m_field.add_ents_to_finite_element_by_type(root_set, MBTET, "FE");
    // build finite elemnts
    CHKERR m_field.build_finite_elements();
    // build adjacencies
    CHKERR m_field.build_adjacencies(BitRefLevel().set());

    // set dm data structure which created mofem data structures
    CHKERR DMMoFEMCreateMoFEM(dm, &m_field, "TEST_PROBLEM",
                              BitRefLevel().set(0));
    CHKERR DMMoFEMSetSquareProblem(
        dm, PETSC_FALSE); // this is for testing (this problem has the same rows
                          // and cols)
    CHKERR DMSetFromOptions(dm);
    CHKERR DMMoFEMAddElement(dm, "FE");
    CHKERR DMSetUp(dm);

    Mat m;
    CHKERR DMCreateMatrix(dm, &m);

    CHKERR m_field.getInterface<MatrixManager>()
        ->checkMPIAIJWithArraysMatrixFillIn<PetscGlobalIdx_mi_tag>(
            "TEST_PROBLEM", -1, -1, 1);

    CHKERR MatDestroy(&m);

    // check if file can be saved
    if (debug)
      CHKERR moab.write_file("test_out.h5m", "MOAB", "PARALLEL=WRITE_PART");

    // destry dm
    CHKERR DMDestroy(&dm);
  }
  CATCH_ERRORS;

  // finish work cleaning memory, getting statistics, ect.
  CHKERR MoFEM::Core::Finalize();

  return 0;
}
