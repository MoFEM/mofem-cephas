/** \file build_composite_problems.cpp
  \example build_composite_problems.cpp
  \brief Atom test for building composite problems

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

    // read mesh and create moab and mofem data structures

    moab::Core mb_instance;
    moab::Interface &moab = mb_instance;
    const char *option;
    option = "";
    CHKERR moab.load_file(mesh_file_name, 0, option);


    MoFEM::Core core(moab);
    MoFEM::Interface &m_field = core;

    EntityHandle root_set = moab.get_root_set();
    Range tets;
    CHKERR moab.get_entities_by_type(root_set, MBTET, tets, false);

    ProblemsManager *prb_mng_ptr;
    CHKERR m_field.getInterface(prb_mng_ptr);
    CHKERR prb_mng_ptr->partitionMesh(tets, 3, 2, m_field.get_comm_size(), NULL,
                                      NULL, NULL);

    EntityHandle part_set;
    CHKERR moab.create_meshset(MESHSET_SET, part_set);

    ParallelComm *pcomm = ParallelComm::get_pcomm(&moab, MYPCOMM_INDEX);
    if (pcomm == NULL)
      SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
              "Communicator should be allocated");

    Tag part_tag = pcomm->part_tag();
    Range proc_ents;
    Range tagged_sets;
    CHKERR m_field.get_moab().get_entities_by_type_and_tag(
        0, MBENTITYSET, &part_tag, NULL, 1, tagged_sets,
        moab::Interface::UNION);
    for (Range::iterator mit = tagged_sets.begin(); mit != tagged_sets.end();
         mit++) {
      int part;
      CHKERR moab.tag_get_data(part_tag, &*mit, 1, &part);
      if (part == m_field.get_comm_rank()) {
        // pcomm->partition_sets().insert(*mit);
        CHKERR moab.get_entities_by_type(*mit, MBTET, proc_ents, true);
        CHKERR moab.add_entities(part_set, proc_ents);
      }
    }

    Skinner skin(&m_field.get_moab());
    Range tets_skin;
    CHKERR skin.find_skin(0, tets, false, tets_skin);
    Range proc_ents_skin[4];
    proc_ents_skin[3] = proc_ents;
    CHKERR skin.find_skin(0, proc_ents, false, proc_ents_skin[2]);
    proc_ents_skin[2] = subtract(proc_ents_skin[2], tets_skin);
    CHKERR moab.get_adjacencies(proc_ents_skin[2], 1, false, proc_ents_skin[1],
                                moab::Interface::UNION);
    CHKERR moab.get_connectivity(proc_ents_skin[1], proc_ents_skin[0], true);
    for (int dd = 0; dd != 3; dd++) {
      CHKERR moab.add_entities(part_set, proc_ents_skin[dd]);
    }

    if (0) {
      std::ostringstream file_skin;
      file_skin << "out_skin_" << m_field.get_comm_rank() << ".vtk";
      EntityHandle meshset_skin;
      CHKERR moab.create_meshset(MESHSET_SET, meshset_skin);
      CHKERR moab.add_entities(meshset_skin, proc_ents_skin[2]);
      CHKERR moab.add_entities(meshset_skin, proc_ents_skin[1]);
      CHKERR moab.add_entities(meshset_skin, proc_ents_skin[0]);
      CHKERR moab.write_file(file_skin.str().c_str(), "VTK", "", &meshset_skin,
                             1);
    }

    CHKERR pcomm->resolve_shared_ents(0, proc_ents, 3, -1, proc_ents_skin);
    Range owned_tets = proc_ents;

    if (0) {
      std::ostringstream file_owned;
      file_owned << "out_owned_" << m_field.get_comm_rank() << ".vtk";
      EntityHandle meshset_owned;
      CHKERR moab.create_meshset(MESHSET_SET, meshset_owned);
      CHKERR moab.add_entities(meshset_owned, owned_tets);
      CHKERR moab.write_file(file_owned.str().c_str(), "VTK", "",
                             &meshset_owned, 1);
    }
    Range shared_ents;
    // Get entities shared with all other processors
    CHKERR pcomm->get_shared_entities(-1, shared_ents);

    if (0) {
      std::ostringstream file_shared_owned;
      file_shared_owned << "out_shared_owned_" << m_field.get_comm_rank()
                        << ".vtk";
      EntityHandle meshset_shared_owned;
      CHKERR moab.create_meshset(MESHSET_SET, meshset_shared_owned);
      CHKERR moab.add_entities(meshset_shared_owned, shared_ents);
      CHKERR moab.write_file(file_shared_owned.str().c_str(), "VTK", "",
                             &meshset_shared_owned, 1);
    }

    // set entities bit level
    BitRefLevel bit_level0;
    bit_level0.set(0);
    CHKERR m_field.getInterface<BitRefManager>()->setBitRefLevelByDim(
        part_set, 3, bit_level0);

    // Fields
    CHKERR m_field.add_field("F1", H1, AINSWORTH_LEGENDRE_BASE, 1);
    CHKERR m_field.add_field("F2", H1, AINSWORTH_LEGENDRE_BASE, 1);
    // add entities to field
    CHKERR m_field.add_ents_to_field_by_type(part_set, MBTET, "F1");
    CHKERR m_field.add_ents_to_field_by_type(part_set, MBTET, "F2");
    int order = 4;
    CHKERR m_field.set_field_order(part_set, MBTET, "F1", order);
    CHKERR m_field.set_field_order(part_set, MBTRI, "F1", order);
    CHKERR m_field.set_field_order(part_set, MBEDGE, "F1", order);
    CHKERR m_field.set_field_order(part_set, MBVERTEX, "F1", 1);
    CHKERR m_field.set_field_order(part_set, MBTET, "F2", order);
    CHKERR m_field.set_field_order(part_set, MBTRI, "F2", order);
    CHKERR m_field.set_field_order(part_set, MBEDGE, "F2", order);
    CHKERR m_field.set_field_order(part_set, MBVERTEX, "F2", 1);
    CHKERR m_field.build_fields();

    // Elements
    CHKERR m_field.add_finite_element("E1");
    CHKERR m_field.add_finite_element("E2");
    CHKERR m_field.modify_finite_element_add_field_row("E1", "F1");
    CHKERR m_field.modify_finite_element_add_field_col("E1", "F1");
    CHKERR m_field.modify_finite_element_add_field_row("E2", "F2");
    CHKERR m_field.modify_finite_element_add_field_col("E2", "F2");
    CHKERR m_field.add_ents_to_finite_element_by_type(part_set, MBTET, "E1");
    CHKERR m_field.add_ents_to_finite_element_by_type(part_set, MBTET, "E2");
    CHKERR m_field.build_finite_elements();
    CHKERR m_field.build_adjacencies(bit_level0);

    // Problems
    CHKERR m_field.add_problem("P1");
    CHKERR m_field.add_problem("P2");
    // set refinement level for problem
    CHKERR m_field.modify_problem_ref_level_add_bit("P1", bit_level0);
    CHKERR m_field.modify_problem_ref_level_add_bit("P2", bit_level0);
    CHKERR m_field.modify_problem_add_finite_element("P1", "E1");
    CHKERR m_field.modify_problem_add_finite_element("P2", "E2");

    // Build problems
    CHKERR prb_mng_ptr->buildProblemOnDistributedMesh("P1", true, 1);
    CHKERR prb_mng_ptr->buildProblemOnDistributedMesh("P2", true, 1);
    CHKERR prb_mng_ptr->partitionFiniteElements("P1", true, 0,
                                                m_field.get_comm_size(), 1);
    CHKERR prb_mng_ptr->partitionGhostDofs("P1", 1);
    CHKERR prb_mng_ptr->partitionFiniteElements("P2", true, 0,
                                                m_field.get_comm_size(), 1);
    CHKERR prb_mng_ptr->partitionGhostDofs("P2", 1);

    if (0) {
      SmartPetscObj<Mat> m;
      CHKERR m_field.getInterface<MatrixManager>()
          ->createMPIAIJWithArrays<PetscGlobalIdx_mi_tag>("P1", m);
      MatView(m, PETSC_VIEWER_DRAW_WORLD);
      std::string wait;
      std::cin >> wait;
    }

    CHKERR m_field.getInterface<MatrixManager>()
        ->checkMPIAIJWithArraysMatrixFillIn<PetscGlobalIdx_mi_tag>("P1", -1, -1,
                                                                   1);
    CHKERR m_field.getInterface<MatrixManager>()
        ->checkMPIAIJWithArraysMatrixFillIn<PetscGlobalIdx_mi_tag>("P2", -1, -1,
                                                                   1);

    // register new dm type, i.e. mofem
    DMType dm_name = "MOFEM";
    CHKERR DMRegister_MoFEM(dm_name);
    // create dm instance
    auto dm = createSmartDM(m_field.get_comm(), dm_name);

    CHKERR DMMoFEMCreateMoFEM(dm, &m_field, "COMP", bit_level0);
    CHKERR DMSetFromOptions(dm);
    CHKERR DMMoFEMSetIsPartitioned(dm, PETSC_TRUE);
    CHKERR DMMoFEMAddElement(dm, "E1");
    CHKERR DMMoFEMAddElement(dm, "E2");
    CHKERR DMMoFEMAddRowCompositeProblem(dm, "P1");
    CHKERR DMMoFEMAddRowCompositeProblem(dm, "P2");
    CHKERR DMSetUp(dm);

    if (0) {
      SmartPetscObj<Mat> m;
      CHKERR m_field.getInterface<MatrixManager>()
          ->createMPIAIJWithArrays<PetscGlobalIdx_mi_tag>("COMP", m);
      MatView(m, PETSC_VIEWER_DRAW_WORLD);
      std::string wait;
      std::cin >> wait;
    }

    CHKERR m_field.getInterface<MatrixManager>()
        ->checkMPIAIJWithArraysMatrixFillIn<PetscGlobalIdx_mi_tag>("COMP", -1,
                                                                   -1, 1);

  }
  CATCH_ERRORS;

  // Finish work cleaning memory, getting statistics, etc.
  CHKERR MoFEM::Core::Finalize();

  return 0;
}
