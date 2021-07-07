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
      SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
              "*** ERROR -my_file (MESH FILE NEEDED)");
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
    CHKERR m_field.add_field("FIELD1", H1, AINSWORTH_LEGENDRE_BASE, 1);
    CHKERR m_field.add_field("FIELD2", H1, AINSWORTH_LEGENDRE_BASE, 3);
    CHKERR m_field.add_field("FIELD3", NOFIELD, NOBASE, 3);
    CHKERR m_field.add_field("MESH_NODE_POSITIONS", H1, AINSWORTH_LEGENDRE_BASE,
                             3);

    {
      // Creating and adding no field entities.
      const double coords[] = {0, 0, 0};
      EntityHandle no_field_vertex;
      CHKERR m_field.get_moab().create_vertex(coords, no_field_vertex);
      Range range_no_field_vertex;
      range_no_field_vertex.insert(no_field_vertex);
      CHKERR m_field.getInterface<BitRefManager>()->setBitRefLevel(
          range_no_field_vertex, BitRefLevel().set());
      EntityHandle meshset = m_field.get_field_meshset("FIELD3");
      CHKERR m_field.get_moab().add_entities(meshset, range_no_field_vertex);
    }

    // FE
    CHKERR m_field.add_finite_element("TEST_FE1");
    CHKERR m_field.add_finite_element("TEST_FE2");

    // Define rows/cols and element data
    CHKERR m_field.modify_finite_element_add_field_row("TEST_FE1", "FIELD1");
    CHKERR m_field.modify_finite_element_add_field_col("TEST_FE1", "FIELD2");
    CHKERR m_field.modify_finite_element_add_field_data("TEST_FE1", "FIELD1");
    CHKERR m_field.modify_finite_element_add_field_data("TEST_FE1", "FIELD2");
    CHKERR m_field.modify_finite_element_add_field_data("TEST_FE1",
                                                        "MESH_NODE_POSITIONS");

    CHKERR m_field.modify_finite_element_add_field_row("TEST_FE2", "FIELD3");
    CHKERR m_field.modify_finite_element_add_field_col("TEST_FE2", "FIELD1");
    CHKERR m_field.modify_finite_element_add_field_data("TEST_FE2", "FIELD1");
    CHKERR m_field.modify_finite_element_add_field_data("TEST_FE2", "FIELD3");
    CHKERR m_field.modify_finite_element_add_field_data("TEST_FE2",
                                                        "MESH_NODE_POSITIONS");

    // Problem
    CHKERR m_field.add_problem("TEST_PROBLEM");

    // set finite elements for problem
    CHKERR m_field.modify_problem_add_finite_element("TEST_PROBLEM",
                                                     "TEST_FE1");
    CHKERR m_field.modify_problem_add_finite_element("TEST_PROBLEM",
                                                     "TEST_FE2");
    // set refinement level for problem
    CHKERR m_field.modify_problem_ref_level_add_bit("TEST_PROBLEM", bit_level0);

    // meshset consisting all entities in mesh
    EntityHandle root_set = moab.get_root_set();
    // add entities to field
    CHKERR m_field.add_ents_to_field_by_type(root_set, MBTET, "FIELD1");
    CHKERR m_field.add_ents_to_field_by_type(root_set, MBTET, "FIELD2");
    CHKERR m_field.add_ents_to_field_by_type(root_set, MBTET,
                                             "MESH_NODE_POSITIONS");

    // add entities to finite element
    CHKERR m_field.add_ents_to_finite_element_by_type(root_set, MBTET,
                                                      "TEST_FE1");
    CHKERR m_field.add_ents_to_finite_element_by_type(root_set, MBTET,
                                                      "TEST_FE2");

    // set app. order
    // see Hierarchic Finite Element Bases on Unstructured Tetrahedral Meshes
    // (Mark Ainsworth & Joe Coyle)
    int order = 4;
    CHKERR m_field.set_field_order(root_set, MBTET, "FIELD1", order);
    CHKERR m_field.set_field_order(root_set, MBTRI, "FIELD1", order);
    CHKERR m_field.set_field_order(root_set, MBEDGE, "FIELD1", order);
    CHKERR m_field.set_field_order(root_set, MBVERTEX, "FIELD1", 1);
    CHKERR m_field.set_field_order(root_set, MBTET, "FIELD2", order);
    CHKERR m_field.set_field_order(root_set, MBTRI, "FIELD2", order);
    CHKERR m_field.set_field_order(root_set, MBEDGE, "FIELD2", order);
    CHKERR m_field.set_field_order(root_set, MBVERTEX, "FIELD2", 1);

    CHKERR m_field.set_field_order(root_set, MBTET, "MESH_NODE_POSITIONS", 2);
    CHKERR m_field.set_field_order(root_set, MBTRI, "MESH_NODE_POSITIONS", 2);
    CHKERR m_field.set_field_order(root_set, MBEDGE, "MESH_NODE_POSITIONS", 2);
    CHKERR m_field.set_field_order(root_set, MBVERTEX, "MESH_NODE_POSITIONS",
                                   1);

    /****/
    // build database
    // build field
    CHKERR m_field.build_fields();
    // build finite elemnts
    CHKERR m_field.build_finite_elements();
    // build adjacencies
    CHKERR m_field.build_adjacencies(bit_level0);
    // build problem
    ProblemsManager *prb_mng_ptr;
    CHKERR m_field.getInterface(prb_mng_ptr);
    //  const Problem_multiIndex *problems_ptr;
    CHKERR prb_mng_ptr->buildProblem("TEST_PROBLEM", false);

    /****/
    // mesh partitioning
    // partition
    CHKERR prb_mng_ptr->partitionSimpleProblem("TEST_PROBLEM");
    CHKERR prb_mng_ptr->partitionFiniteElements("TEST_PROBLEM");
    // what are ghost nodes, see Petsc Manual
    CHKERR prb_mng_ptr->partitionGhostDofs("TEST_PROBLEM");

    // set from positions of 10 node tets
    Projection10NodeCoordsOnField ent_method(m_field, "MESH_NODE_POSITIONS");
    CHKERR m_field.loop_dofs("MESH_NODE_POSITIONS", ent_method);

    typedef tee_device<std::ostream, std::ofstream> TeeDevice;
    typedef stream<TeeDevice> TeeStream;

    std::ofstream ofs("forces_and_sources_testing_volume_element.txt");
    TeeDevice my_tee(std::cout, ofs);
    TeeStream my_split(my_tee);

    struct MyOp1 : public VolumeElementForcesAndSourcesCore::UserDataOperator {

      TeeStream &my_split;
      MyOp1(TeeStream &_my_split, char type)
          : VolumeElementForcesAndSourcesCore::UserDataOperator("FIELD1",
                                                                "FIELD2", type),
            my_split(_my_split) {}

      MoFEMErrorCode doWork(int side, EntityType type,
                            DataForcesAndSourcesCore::EntData &data) {
        MoFEMFunctionBeginHot;
        my_split << "NH1" << std::endl;
        my_split << "side: " << side << " type: " << type << std::endl;
        my_split << data << std::endl;
        my_split << std::setprecision(3) << getVolume() << std::endl;
        my_split << std::setprecision(3) << getCoords() << std::endl;
        my_split << std::setprecision(3) << getCoordsAtGaussPts() << std::endl;
        MoFEMFunctionReturnHot(0);
      }

      MoFEMErrorCode doWork(int row_side, int col_side, EntityType row_type,
                            EntityType col_type,
                            DataForcesAndSourcesCore::EntData &row_data,
                            DataForcesAndSourcesCore::EntData &col_data) {
        MoFEMFunctionBeginHot;
        my_split << "NH1NH1" << std::endl;
        my_split << "row side: " << row_side << " row_type: " << row_type
                 << std::endl;
        my_split << row_data << std::endl;
        my_split << "NH1NH1" << std::endl;
        my_split << "col side: " << col_side << " col_type: " << col_type
                 << std::endl;
        my_split << col_data << std::endl;

        VectorInt row_indices, col_indices;
        CHKERR getProblemRowIndices("FIELD1", row_type, row_side, row_indices);
        CHKERR getProblemColIndices("FIELD2", col_type, col_side, col_indices);

        if (row_indices.size() != row_data.getIndices().size()) {
          std::cerr << row_indices << std::endl;
          std::cerr << row_data.getIndices() << std::endl;
          SETERRQ2(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
                   "row inconsistency %d != %d", row_indices.size(),
                   row_data.getIndices().size());
        }

        if (col_indices.size() != col_data.getIndices().size()) {
          SETERRQ2(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
                   "row inconsistency %d != %d", col_indices.size(),
                   col_data.getIndices().size());
        }

        for (unsigned int rr = 0; rr < row_indices.size(); rr++) {
          if (row_indices[rr] != row_data.getIndices()[rr]) {
            std::cerr << row_indices << std::endl;
            std::cerr << row_data.getIndices() << std::endl;
            SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
                    "row inconsistency");
          }
        }

        for (unsigned int cc = 0; cc < col_indices.size(); cc++) {
          if (col_indices[cc] != col_data.getIndices()[cc]) {
            std::cerr << col_indices << std::endl;
            std::cerr << col_data.getIndices() << std::endl;
            SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
                    "row inconsistency");
          }
        }

        MoFEMFunctionReturnHot(0);
      }
    };

    struct MyOp2 : public VolumeElementForcesAndSourcesCore::UserDataOperator {

      TeeStream &my_split;
      MyOp2(TeeStream &_my_split, char type)
          : VolumeElementForcesAndSourcesCore::UserDataOperator("FIELD3",
                                                                "FIELD1", type),
            my_split(_my_split) {
        sYmm = 0;
      }

      MoFEMErrorCode doWork(int side, EntityType type,
                            DataForcesAndSourcesCore::EntData &data) {
        MoFEMFunctionBeginHot;

        if (type != MBENTITYSET)
          MoFEMFunctionReturnHot(0);

        my_split << "NOFIELD" << std::endl;
        my_split << "side: " << side << " type: " << type << std::endl;
        my_split << data << std::endl;
        MoFEMFunctionReturnHot(0);
      }

      MoFEMErrorCode doWork(int row_side, int col_side, EntityType row_type,
                            EntityType col_type,
                            DataForcesAndSourcesCore::EntData &row_data,
                            DataForcesAndSourcesCore::EntData &col_data) {
        MoFEMFunctionBeginHot;

        unSetSymm();

        if (row_type != MBENTITYSET)
          MoFEMFunctionReturnHot(0);

        my_split << "NOFILEDH1" << std::endl;
        my_split << "row side: " << row_side << " row_type: " << row_type
                 << std::endl;
        my_split << row_data << std::endl;
        my_split << "col side: " << col_side << " col_type: " << col_type
                 << std::endl;
        my_split << col_data << std::endl;

        MoFEMFunctionReturnHot(0);
      }
    };

    VolumeElementForcesAndSourcesCore fe1(m_field);
    fe1.meshPositionsFieldName = "none";

    fe1.getOpPtrVector().push_back(
        new MyOp1(my_split, ForcesAndSourcesCore::UserDataOperator::OPROW));
    fe1.getOpPtrVector().push_back(
        new MyOp1(my_split, ForcesAndSourcesCore::UserDataOperator::OPROWCOL));

    VolumeElementForcesAndSourcesCore fe2(m_field);
    fe2.getOpPtrVector().push_back(
        new MyOp2(my_split, ForcesAndSourcesCore::UserDataOperator::OPROW));
    fe2.getOpPtrVector().push_back(
        new MyOp2(my_split, ForcesAndSourcesCore::UserDataOperator::OPROWCOL));

    CHKERR m_field.loop_finite_elements("TEST_PROBLEM", "TEST_FE1", fe1);
    CHKERR m_field.loop_finite_elements("TEST_PROBLEM", "TEST_FE2", fe2);
  }
  CATCH_ERRORS;

  MoFEM::Core::Finalize();

  return 0;
}
