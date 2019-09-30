/** \file normals_at_integration_points_for_quad.cpp
  \example normals_at_integration_points_for_quad
  \brief Testing normals at integration points for a quad face

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


namespace bio = boost::iostreams;
using bio::stream;
using bio::tee_device;

using namespace MoFEM;
static char help[] = "...\n\n";
static int debug = 1;

using EntData = DataForcesAndSourcesCore::EntData;

static double sum_matrix(MatrixDouble &m) {
  double s = 0;
  for (unsigned int ii = 0; ii < m.size1(); ii++) {
    for (unsigned int jj = 0; jj < m.size2(); jj++) {
      s += m(ii, jj);
    }
  }
  return s;
}

int main(int argc, char *argv[]) {

  MoFEM::Core::Initialize(&argc, &argv, (char *)0, help);
  try {
    moab::Core mb_instance;
    moab::Interface &moab = mb_instance;
    ParallelComm *pcomm = ParallelComm::get_pcomm(&moab, MYPCOMM_INDEX);
    if (pcomm == NULL)
      pcomm = new ParallelComm(&moab, PETSC_COMM_WORLD);
    MoFEM::Core core(moab);
    MoFEM::Interface &m_field = core;

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
    option = ""; //"PARALLEL=BCAST;";//;DEBUG_IO";
    CHKERR moab.load_file(mesh_file_name, 0, option);

    // PrismsFromSurfaceInterface *prisms_from_surface_interface;
    // CHKERR m_field.getInterface(prisms_from_surface_interface);

    // Range tris;
    // CHKERR moab.get_entities_by_type(0, MBTRI, tris, false);
    // Range prisms;
    // CHKERR prisms_from_surface_interface->createPrisms(tris, prisms);
    // prisms_from_surface_interface->createdVertices.clear();

    // double d0[3] = {0, 0, 0};
    // double d1[3] = {0, 0, 1};
    // CHKERR prisms_from_surface_interface->setThickness(prisms, d0, d1);

    // Range faces, quads;
    // CHKERR moab.get_adjacencies(prisms, 2, true, faces, moab::Interface::UNION);
    // quads = faces.subset_by_type(MBQUAD);

    // EntityHandle meshset;
    // CHKERR moab.create_meshset(MESHSET_SET | MESHSET_TRACK_OWNER, meshset);
    // CHKERR moab.add_entities(meshset, prisms);

    BitRefLevel bit_level0;
    bit_level0.set(0);
    CHKERR m_field.getInterface<BitRefManager>()->setBitRefLevelByDim(
        0, 2, bit_level0);

    // Fields
    CHKERR m_field.add_field("FIELD1", H1, AINSWORTH_LEGENDRE_BASE, 1);
    CHKERR m_field.add_ents_to_field_by_type(0, MBQUAD, "FIELD1");

    int order = 1;

    CHKERR m_field.set_field_order(0, MBVERTEX, "FIELD1", 1);
    CHKERR m_field.set_field_order(0, MBEDGE, "FIELD1", order);
    CHKERR m_field.set_field_order(0, MBQUAD, "FIELD1", order);

    CHKERR m_field.build_fields();

    // FE
    CHKERR m_field.add_finite_element("TEST_FE1");

    // Define rows/cols and element data
    CHKERR m_field.modify_finite_element_add_field_row("TEST_FE1", "FIELD1");
    CHKERR m_field.modify_finite_element_add_field_col("TEST_FE1", "FIELD1");
    CHKERR m_field.modify_finite_element_add_field_data("TEST_FE1", "FIELD1");

    CHKERR m_field.add_ents_to_finite_element_by_type(0, MBQUAD,
                                                      "TEST_FE1");

    // build finite elemnts
    CHKERR m_field.build_finite_elements();
    // //build adjacencies
    CHKERR m_field.build_adjacencies(bit_level0);

    // Problem
    CHKERR m_field.add_problem("TEST_PROBLEM");

    // set finite elements for problem
    CHKERR m_field.modify_problem_add_finite_element("TEST_PROBLEM",
                                                     "TEST_FE1");
    // set refinement level for problem
    CHKERR m_field.modify_problem_ref_level_add_bit("TEST_PROBLEM", bit_level0);

    // build problem
    ProblemsManager *prb_mng_ptr;
    CHKERR m_field.getInterface(prb_mng_ptr);
    CHKERR prb_mng_ptr->buildProblem("TEST_PROBLEM", true);
    // partition
    CHKERR prb_mng_ptr->partitionSimpleProblem("TEST_PROBLEM");
    CHKERR prb_mng_ptr->partitionFiniteElements("TEST_PROBLEM");
    // what are ghost nodes, see Petsc Manual
    CHKERR prb_mng_ptr->partitionGhostDofs("TEST_PROBLEM");


    EntityHandle root_set = moab.get_root_set();
    CHKERR moab.write_file("sphere.vtk", "VTK", "", &root_set, 1);
  }
  CATCH_ERRORS;

  CHKERR MoFEM::Core::Finalize();
}
