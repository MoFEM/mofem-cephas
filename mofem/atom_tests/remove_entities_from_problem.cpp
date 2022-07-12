/** \file remove_entities_from_problem.cpp
  \example remove_entities_from_problem.cpp
  \brief Remove field entities from the problem

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
    if (!flg)
      SETERRQ(PETSC_COMM_WORLD, 1, "*** ERROR -my_file (MESH FILE NEEDED)");
    
    ParallelComm *pcomm = ParallelComm::get_pcomm(&moab, MYPCOMM_INDEX);
    if (pcomm == NULL)
      pcomm = new ParallelComm(&moab, PETSC_COMM_WORLD);
    const char *option;
    option = "PARALLEL=BCAST_DELETE;PARALLEL_RESOLVE_SHARED_ENTS;PARTITION="
             "PARALLEL_PARTITION;";
    CHKERR moab.load_file(mesh_file_name, 0, option);

    // Create MoFEM database
    MoFEM::Core core(moab);
    MoFEM::Interface &m_field = core;

    // set entities bit level
    const auto bit_level = BitRefLevel().set(0);
    CHKERR m_field.getInterface<BitRefManager>()->setBitRefLevelByDim(
        0, 3, bit_level);

    // Fields
    CHKERR m_field.add_field("F1", HCURL, DEMKOWICZ_JACOBI_BASE, 3);
    CHKERR m_field.add_field("F2", HDIV, DEMKOWICZ_JACOBI_BASE, 1);

    // meshset consisting all entities in mesh
    EntityHandle root_set = moab.get_root_set();
    // add entities to field
    CHKERR m_field.add_ents_to_field_by_type(root_set, MBTET, "F1");
    CHKERR m_field.add_ents_to_field_by_type(root_set, MBTET, "F2");

    // set app. order
    // see Hierarchic Finite Element Bases on Unstructured Tetrahedral Meshes
    // (Mark Ainsworth & Joe Coyle)
    int order = 2;
    CHKERR m_field.set_field_order(root_set, MBEDGE, "F1", order);
    CHKERR m_field.set_field_order(root_set, MBTRI, "F1", order);
    CHKERR m_field.set_field_order(root_set, MBTET, "F1", order);
    CHKERR m_field.set_field_order(root_set, MBTET, "F2", order);
    CHKERR m_field.set_field_order(root_set, MBTRI, "F2", order);

    CHKERR m_field.build_fields();

    // add elements
    CHKERR m_field.add_finite_element("E1");

    CHKERR m_field.modify_finite_element_add_field_row("E1", "F1");
    CHKERR m_field.modify_finite_element_add_field_row("E1", "F2");
    CHKERR m_field.modify_finite_element_add_field_col("E1", "F2");
    CHKERR m_field.modify_finite_element_add_field_col("E1", "F1");
    CHKERR m_field.modify_finite_element_add_field_data("E1", "F1");
    CHKERR m_field.modify_finite_element_add_field_data("E1", "F2");

    CHKERR m_field.add_ents_to_finite_element_by_type(root_set, MBTET, "E1");

    CHKERR m_field.build_finite_elements();
    CHKERR m_field.build_adjacencies(bit_level);

    // Problems
    CHKERR m_field.add_problem("P1");

    // set refinement level for problem
    CHKERR m_field.modify_problem_ref_level_add_bit("P1", bit_level);
    CHKERR m_field.modify_problem_add_finite_element("P1", "E1");

    // build problems
    ProblemsManager *prb_mng_ptr;
    CHKERR m_field.getInterface(prb_mng_ptr);
    CHKERR prb_mng_ptr->buildProblemOnDistributedMesh("P1", true);
    CHKERR prb_mng_ptr->partitionFiniteElements("P1");
    CHKERR prb_mng_ptr->partitionGhostDofs("P1");

    auto get_triangles_on_skin = [&](Range &tets_skin) {
      MoFEMFunctionBegin;
      Range tets;
      CHKERR m_field.get_moab().get_entities_by_type(root_set, MBTET, tets);
      Range tets_skin_part;
      Skinner skin(&m_field.get_moab());
      CHKERR skin.find_skin(0, tets, false, tets_skin_part);
      CHKERR pcomm->filter_pstatus(tets_skin_part,
                                   PSTATUS_SHARED | PSTATUS_MULTISHARED,
                                   PSTATUS_NOT, -1, &tets_skin);
      MoFEMFunctionReturn(0);
    };

    Range tets_skin;
    CHKERR get_triangles_on_skin(tets_skin);
    CHKERR prb_mng_ptr->removeDofsOnEntities("P1", "F1", tets_skin, 0, 1, NOISY,
                                             true);

    CHKERR m_field.getInterface<MatrixManager>()
        ->checkMPIAIJWithArraysMatrixFillIn<PetscGlobalIdx_mi_tag>("P1", -2, -2,
                                                                   0);
  }
  CATCH_ERRORS;

  // finish work cleaning memory, getting statistics, etc.
  CHKERR MoFEM::Core::Finalize();

  return 0;
}
