/** \file  forces_and_sources_getting_higher_order_skin_normals_atom_tets.cpp

  \brief Atom test to calculate normals on faces approximated with HO geometry
  representation

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

    // Create MoFEM (Joseph) database
    MoFEM::Core core(moab);
    MoFEM::Interface &m_field = core;

    const char *option;
    option = ""; 
    CHKERR moab.load_file(mesh_file_name, 0, option);

    // set entitities bit level
    BitRefLevel bit_level0;
    bit_level0.set(0);
    EntityHandle meshset_level0;
    CHKERR moab.create_meshset(MESHSET_SET, meshset_level0);
    CHKERR m_field.getInterface<BitRefManager>()->setBitRefLevelByDim(
        0, 3, bit_level0);

    // Fields
    CHKERR m_field.add_field("FIELD1", H1, AINSWORTH_LEGENDRE_BASE, 3);

    // FE
    CHKERR m_field.add_finite_element("TEST_FE");

    // Define rows/cols and element data
    CHKERR m_field.modify_finite_element_add_field_row("TEST_FE", "FIELD1");
    CHKERR m_field.modify_finite_element_add_field_col("TEST_FE", "FIELD1");
    CHKERR m_field.modify_finite_element_add_field_data("TEST_FE", "FIELD1");

    // Problem
    CHKERR m_field.add_problem("TEST_PROBLEM");

    // set finite elements for problem
    CHKERR m_field.modify_problem_add_finite_element("TEST_PROBLEM", "TEST_FE");
    // set refinement level for problem
    CHKERR m_field.modify_problem_ref_level_add_bit("TEST_PROBLEM", bit_level0);

    // meshset consisting all entities in mesh
    EntityHandle root_set = moab.get_root_set();
    // add entities to field
    CHKERR m_field.add_ents_to_field_by_type(root_set, MBTET, "FIELD1");
    // add entities to finite element
    Range tets;
    CHKERR moab.get_entities_by_type(0, MBTET, tets, false);
    Skinner skin(&m_field.get_moab());
    Range tets_skin;
    CHKERR skin.find_skin(0, tets, false, tets_skin);
    CHKERR m_field.add_ents_to_finite_element_by_type(tets_skin, MBTRI,
                                                      "TEST_FE");

    // set app. order
    // see Hierarchic Finite Element Bases on Unstructured Tetrahedral Meshes
    // (Mark Ainsworth & Joe Coyle)
    int order = 3;
    CHKERR m_field.set_field_order(root_set, MBTET, "FIELD1", order);
    CHKERR m_field.set_field_order(root_set, MBTRI, "FIELD1", order);
    CHKERR m_field.set_field_order(root_set, MBEDGE, "FIELD1", order);
    CHKERR m_field.set_field_order(root_set, MBVERTEX, "FIELD1", 1);

    /****/
    // build database
    // build field
    CHKERR m_field.build_fields();
    // set FIELD1 from positions of 10 node tets
    Projection10NodeCoordsOnField ent_method(m_field, "FIELD1");
    CHKERR m_field.loop_dofs("FIELD1", ent_method);
    // build finite elemnts
    CHKERR m_field.build_finite_elements();
    // build adjacencies
    CHKERR m_field.build_adjacencies(bit_level0);
    // build problem
    ProblemsManager *prb_mng_ptr;
    CHKERR m_field.getInterface(prb_mng_ptr);
    CHKERR prb_mng_ptr->buildProblem("TEST_PROBLEM", true);

    /****/
    // mesh partitioning
    // partition
    CHKERR prb_mng_ptr->partitionSimpleProblem("TEST_PROBLEM");
    CHKERR prb_mng_ptr->partitionFiniteElements("TEST_PROBLEM");
    // what are ghost nodes, see Petsc Manual
    CHKERR prb_mng_ptr->partitionGhostDofs("TEST_PROBLEM");

    struct ForcesAndSourcesCore_TestFE
        : public FaceElementForcesAndSourcesCoreBase {

      typedef tee_device<std::ostream, std::ofstream> TeeDevice;
      typedef stream<TeeDevice> TeeStream;
      std::ofstream ofs;
      TeeDevice my_tee;
      TeeStream my_split;

      ForcesAndSourcesCore_TestFE(MoFEM::Interface &m_field)
          : FaceElementForcesAndSourcesCoreBase(m_field),
            ofs("forces_and_sources_getting_higher_order_skin_normals_atom."
                "txt"),
            my_tee(std::cout, ofs), my_split(my_tee){};

      MoFEMErrorCode preProcess() {
        MoFEMFunctionBeginHot;
        MoFEMFunctionReturnHot(0);
      }

      MoFEMErrorCode operator()() {
        MoFEMFunctionBegin;

        CHKERR opSwitch<0>();

        my_split.precision(3);
        my_split << "coords: " << coordsAtGaussPts << std::endl;
        my_split << "normals: " << normalsAtGaussPts << std::endl;
        my_split << "tangent1: " << tangentOneAtGaussPts << std::endl;
        my_split << "tangent2: " << tangentTwoAtGaussPts << std::endl;

        MoFEMFunctionReturn(0);
      }

      MoFEMErrorCode postProcess() {
        MoFEMFunctionBeginHot;
        MoFEMFunctionReturnHot(0);
      }
    };

    ForcesAndSourcesCore_TestFE fe1(m_field);
    fe1.meshPositionsFieldName = "FIELD1";

    fe1.getRuleHook = [](int, int, int approx_order) { return -1; };
    fe1.setRuleHook = [](ForcesAndSourcesCore *fe_ptr, int ro, int co, int ao) {
      MoFEMFunctionBegin;
      auto &gauss_pts = fe_ptr->gaussPts;
      gauss_pts.resize(3, 4, false);
      for (int gg = 0; gg < 4; gg++) {
        gauss_pts(0, gg) = G_TRI_X4[gg];
        gauss_pts(1, gg) = G_TRI_Y4[gg];
        gauss_pts(2, gg) = 0;
      }
      MoFEMFunctionReturn(0);
    };

    fe1.getOpPtrVector().push_back(new OpCalculateHOCoords("FIELD1"));
    fe1.getOpPtrVector().push_back(new OpGetHONormalsOnFace("FIELD1"));

    CHKERR m_field.loop_finite_elements("TEST_PROBLEM", "TEST_FE", fe1);
  }
  CATCH_ERRORS;

  CHKERR MoFEM::Core::Finalize();

  return 0;
}
