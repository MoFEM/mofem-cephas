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

// Rounding
#define RND_EPS 1e-6
double roundn(double n) {

  // break n into fractional part (fract) and integral part (intp)
  double fract, intp;
  fract = modf(n, &intp);

  // case where n approximates zero, set n to "positive" zero
  if (std::abs(intp) == 0) {
    if (std::abs(fract) <= RND_EPS) {
      n = 0.000;
    }
  }

  return n;
}

int main(int argc, char *argv[]) {

  MoFEM::Core::Initialize(&argc, &argv, (char *)0, help);

  try {

    moab::Core mb_instance;
    moab::Interface &moab = mb_instance;
    int rank;
    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
    /*if(rank==0) {
      EntityHandle dummy_meshset;
      CHKERR moab.create_meshset(MESHSET_SET,dummy_meshset);
    }*/

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

    MoFEM::Core core(moab);
    MoFEM::Interface &m_field = core;

    // add filds
    CHKERR m_field.add_field("MESH_NODE_POSITIONS", H1, AINSWORTH_LEGENDRE_BASE,
                             3);

    // add finite elements
    CHKERR m_field.add_finite_element("TET_ELEM");
    CHKERR m_field.modify_finite_element_add_field_row("TET_ELEM",
                                                       "MESH_NODE_POSITIONS");
    CHKERR m_field.modify_finite_element_add_field_col("TET_ELEM",
                                                       "MESH_NODE_POSITIONS");
    CHKERR m_field.modify_finite_element_add_field_data("TET_ELEM",
                                                        "MESH_NODE_POSITIONS");

    // add problems
    // CHKERR m_field.add_problem("EDGE_PROJECTOR_PROBLEM");
    CHKERR m_field.add_problem("TET_PROBLEM");

    // define problems and finite elements
    CHKERR m_field.modify_problem_add_finite_element("TET_PROBLEM", "TET_ELEM");

    BitRefLevel bit_level0;
    bit_level0.set(0);
    CHKERR m_field.getInterface<BitRefManager>()->setBitRefLevelByDim(
        0, 3, bit_level0);
    Range tets;
    CHKERR moab.get_entities_by_type(0, MBTET, tets, true);
    Range edges;
    CHKERR moab.get_entities_by_type(0, MBEDGE, edges, true);
    CHKERR m_field.getInterface<BitRefManager>()->setElementsBitRefLevel(edges);

    // add ents to field and set app. order

    CHKERR m_field.add_ents_to_field_by_type(0, MBTET, "MESH_NODE_POSITIONS");
    CHKERR m_field.set_field_order(0, MBVERTEX, "MESH_NODE_POSITIONS", 1);
    CHKERR m_field.set_field_order(0, MBEDGE, "MESH_NODE_POSITIONS", 2);
    CHKERR m_field.set_field_order(0, MBTRI, "MESH_NODE_POSITIONS", 2);
    CHKERR m_field.set_field_order(0, MBTET, "MESH_NODE_POSITIONS", 2);

    // add finite elements entities
    CHKERR m_field.add_ents_to_finite_element_by_type(tets, MBTET, "TET_ELEM");

    // set problem level
    CHKERR m_field.modify_problem_ref_level_add_bit("TET_PROBLEM", bit_level0);

    // build fields
    CHKERR m_field.build_fields();
    // build finite elements
    CHKERR m_field.build_finite_elements();
    // build adjacencies
    CHKERR m_field.build_adjacencies(bit_level0);
    // build problem
    ProblemsManager *prb_mng_ptr;
    CHKERR m_field.getInterface(prb_mng_ptr);
    CHKERR prb_mng_ptr->buildProblem("TET_PROBLEM", false);

    // partition
    CHKERR prb_mng_ptr->partitionProblem("TET_PROBLEM");
    CHKERR prb_mng_ptr->partitionFiniteElements("TET_PROBLEM");
    // what are ghost nodes, see Petsc Manual
    CHKERR prb_mng_ptr->partitionGhostDofs("TET_PROBLEM");

    Projection10NodeCoordsOnField ent_method(m_field, "MESH_NODE_POSITIONS");
    CHKERR m_field.loop_dofs("MESH_NODE_POSITIONS", ent_method);

    // Open mesh_file_name.txt for writing
    std::ofstream myfile;
    myfile.open("10node_sphere.txt");

    // Output displacements
    std::cout << "<<<< Dofs (X-Translation, Y-Translation, Z-Translation) >>>>>"
              << std::endl;
    myfile << "<<<< Dofs (X-Translation, Y-Translation, Z-Translation) >>>>>"
           << std::endl;

    auto dofs_ptr = m_field.get_dofs();
    DofEntity_multiIndex_uid_view dofs_view;
    for (_IT_GET_DOFS_FIELD_BY_NAME_FOR_LOOP_(m_field, "MESH_NODE_POSITIONS",
                                              dof_ptr)) {
      dofs_view.insert(*dof_ptr);
    }
    for (DofEntity_multiIndex_uid_view::iterator dit = dofs_view.begin();
         dit != dofs_view.end(); dit++) {
      // if(dof_ptr->getEntType()!=MBEDGE) continue;

      if ((*dit)->getDofCoeffIdx() == 0) {
        // Round and truncate to 3 decimal places
        double fval = (*dit)->getFieldData();
        std::cout << boost::format("%.3lf") % roundn(fval) << "  ";
        myfile << boost::format("%.3lf") % roundn(fval) << "  ";
      }
      if ((*dit)->getDofCoeffIdx() == 1) {
        // Round and truncate to 3 decimal places
        double fval = (*dit)->getFieldData();
        std::cout << boost::format("%.3lf") % roundn(fval) << "  ";
        myfile << boost::format("%.3lf") % roundn(fval) << "  ";
      }
      if ((*dit)->getDofCoeffIdx() == 2) {
        // Round and truncate to 3 decimal places
        double fval = (*dit)->getFieldData();
        std::cout << boost::format("%.3lf") % roundn(fval) << std::endl;
        myfile << boost::format("%.3lf") % roundn(fval) << std::endl;
      }
    }
  }
  CATCH_ERRORS;

  MoFEM::Core::Finalize();
}
