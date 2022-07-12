/** \file forces_and_sources_testing_flat_prism_element.cpp
 * \brief test for flat prism element
 * \example forces_and_sources_testing_flat_prism_element.cpp
 *
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

    // Create MoFEM (Joseph) database
    MoFEM::Core core(moab);
    MoFEM::Interface &m_field = core;
    PrismInterface *interface;
    CHKERR m_field.getInterface(interface);

    // set entitities bit level
    CHKERR m_field.getInterface<BitRefManager>()->setBitRefLevelByDim(
        0, 3, BitRefLevel().set(0));
    std::vector<BitRefLevel> bit_levels;
    bit_levels.push_back(BitRefLevel().set(0));

    int ll = 1;
    // for(_IT_CUBITMESHSETS_BY_BCDATA_TYPE_FOR_LOOP_(m_field,SIDESET|INTERFACESET,cit))
    // {
    for (_IT_CUBITMESHSETS_BY_SET_TYPE_FOR_LOOP_(m_field, SIDESET, cit)) {
      CHKERR PetscPrintf(PETSC_COMM_WORLD, "Insert Interface %d\n",
                         cit->getMeshsetId());
      EntityHandle cubit_meshset = cit->getMeshset();
      {
        // get tet enties form back bit_level
        EntityHandle ref_level_meshset = 0;
        CHKERR moab.create_meshset(MESHSET_SET, ref_level_meshset);
        ;
        CHKERR m_field.getInterface<BitRefManager>()
            ->getEntitiesByTypeAndRefLevel(bit_levels.back(),
                                           BitRefLevel().set(), MBTET,
                                           ref_level_meshset);
        CHKERR m_field.getInterface<BitRefManager>()
            ->getEntitiesByTypeAndRefLevel(bit_levels.back(),
                                           BitRefLevel().set(), MBPRISM,
                                           ref_level_meshset);
        Range ref_level_tets;
        CHKERR moab.get_entities_by_handle(ref_level_meshset, ref_level_tets,
                                           true);
        ;
        // get faces and test to split
        CHKERR interface->getSides(cubit_meshset, bit_levels.back(), true, 0);
        // set new bit level
        bit_levels.push_back(BitRefLevel().set(ll++));
        // split faces and
        CHKERR interface->splitSides(ref_level_meshset, bit_levels.back(),
                                     cubit_meshset, true, true, 0);
        // clean meshsets
        CHKERR moab.delete_entities(&ref_level_meshset, 1);
        ;
      }
      // update cubit meshsets
      for (_IT_CUBITMESHSETS_FOR_LOOP_(m_field, ciit)) {
        EntityHandle cubit_meshset = ciit->meshset;
        CHKERR m_field.getInterface<BitRefManager>()
            ->updateMeshsetByEntitiesChildren(cubit_meshset, bit_levels.back(),
                                              cubit_meshset, MBMAXTYPE, true);
      }
    }

    // Fields
    CHKERR m_field.add_field("FIELD1", H1, AINSWORTH_LEGENDRE_BASE, 3);
    CHKERR m_field.add_field("MESH_NODE_POSITIONS", H1, AINSWORTH_LEGENDRE_BASE,
                             3);
    CHKERR m_field.add_field("FIELD2", NOFIELD, NOBASE, 3);

    {
      // Creating and adding no field entities.
      const double coords[] = {0, 0, 0};
      EntityHandle no_field_vertex;
      CHKERR m_field.get_moab().create_vertex(coords, no_field_vertex);
      ;
      Range range_no_field_vertex;
      range_no_field_vertex.insert(no_field_vertex);
      CHKERR m_field.getInterface<BitRefManager>()->setBitRefLevel(
          range_no_field_vertex, BitRefLevel().set());
      EntityHandle meshset = m_field.get_field_meshset("FIELD2");
      CHKERR m_field.get_moab().add_entities(meshset, range_no_field_vertex);
      ;
    }

    // FE
    CHKERR m_field.add_finite_element("TEST_FE1");
    CHKERR m_field.add_finite_element("TEST_FE2");

    // Define rows/cols and element data
    CHKERR m_field.modify_finite_element_add_field_row("TEST_FE1", "FIELD1");
    CHKERR m_field.modify_finite_element_add_field_col("TEST_FE1", "FIELD1");
    CHKERR m_field.modify_finite_element_add_field_data("TEST_FE1", "FIELD1");
    CHKERR m_field.modify_finite_element_add_field_data("TEST_FE1",
                                                        "MESH_NODE_POSITIONS");

    CHKERR m_field.modify_finite_element_add_field_row("TEST_FE2", "FIELD1");
    // CHKERR m_field.modify_finite_element_add_field_row("TEST_FE2","FIELD2");
    // CHKERR m_field.modify_finite_element_add_field_col("TEST_FE2","FIELD1");
    CHKERR m_field.modify_finite_element_add_field_col("TEST_FE2", "FIELD2");
    CHKERR m_field.modify_finite_element_add_field_data("TEST_FE2", "FIELD1");
    CHKERR m_field.modify_finite_element_add_field_data("TEST_FE2", "FIELD2");

    // Problem
    CHKERR m_field.add_problem("TEST_PROBLEM");

    // set finite elements for problem
    CHKERR m_field.modify_problem_add_finite_element("TEST_PROBLEM",
                                                     "TEST_FE1");
    CHKERR m_field.modify_problem_add_finite_element("TEST_PROBLEM",
                                                     "TEST_FE2");
    // set refinement level for problem
    CHKERR m_field.modify_problem_ref_level_add_bit("TEST_PROBLEM",
                                                    bit_levels.back());

    // meshset consisting all entities in mesh
    EntityHandle root_set = moab.get_root_set();
    // add entities to field
    CHKERR m_field.add_ents_to_field_by_type(root_set, MBTET, "FIELD1");
    CHKERR m_field.add_ents_to_field_by_type(root_set, MBTET,
                                             "MESH_NODE_POSITIONS");
    // add entities to finite element
    CHKERR m_field.add_ents_to_finite_element_by_type(root_set, MBPRISM,
                                                      "TEST_FE1", 10);
    CHKERR m_field.add_ents_to_finite_element_by_type(root_set, MBPRISM,
                                                      "TEST_FE2", 10);

    // set app. order
    // see Hierarchic Finite Element Bases on Unstructured Tetrahedral Meshes
    // (Mark Ainsworth & Joe Coyle)
    int order = 3;
    CHKERR m_field.set_field_order(root_set, MBTET, "FIELD1", order);
    CHKERR m_field.set_field_order(root_set, MBTRI, "FIELD1", order);
    CHKERR m_field.set_field_order(root_set, MBEDGE, "FIELD1", order);
    CHKERR m_field.set_field_order(root_set, MBVERTEX, "FIELD1", 1);

    CHKERR m_field.set_field_order(root_set, MBTET, "MESH_NODE_POSITIONS", 2);
    CHKERR m_field.set_field_order(root_set, MBTRI, "MESH_NODE_POSITIONS", 2);
    CHKERR m_field.set_field_order(root_set, MBEDGE, "MESH_NODE_POSITIONS", 2);
    CHKERR m_field.set_field_order(root_set, MBVERTEX, "MESH_NODE_POSITIONS",
                                   1);

    /****/
    // build database
    // build field
    CHKERR m_field.build_fields();
    // set FIELD1 from positions of 10 node tets
    Projection10NodeCoordsOnField ent_method_field1(m_field, "FIELD1");
    CHKERR m_field.loop_dofs("FIELD1", ent_method_field1);
    Projection10NodeCoordsOnField ent_method_mesh_positions(
        m_field, "MESH_NODE_POSITIONS");
    CHKERR m_field.loop_dofs("MESH_NODE_POSITIONS", ent_method_mesh_positions);

    // build finite elemnts
    CHKERR m_field.build_finite_elements();
    // build adjacencies
    CHKERR m_field.build_adjacencies(bit_levels.back());
    // build problem
    ProblemsManager *prb_mng_ptr;
    CHKERR m_field.getInterface(prb_mng_ptr);
    CHKERR prb_mng_ptr->buildProblem("TEST_PROBLEM", false);

    /****/
    // mesh partitioning
    // partition
    CHKERR prb_mng_ptr->partitionSimpleProblem("TEST_PROBLEM");
    CHKERR prb_mng_ptr->partitionFiniteElements("TEST_PROBLEM");
    // what are ghost nodes, see Petsc Manual
    CHKERR prb_mng_ptr->partitionGhostDofs("TEST_PROBLEM");

    typedef tee_device<std::ostream, std::ofstream> TeeDevice;
    typedef stream<TeeDevice> TeeStream;

    std::ofstream ofs("forces_and_sources_testing_flat_prism_element.txt");
    TeeDevice my_tee(std::cout, ofs);
    TeeStream my_split(my_tee);

    struct MyOp
        : public FlatPrismElementForcesAndSourcesCore::UserDataOperator {

      TeeStream &mySplit;
      MyOp(TeeStream &mySplit, const char type)
          : FlatPrismElementForcesAndSourcesCore::UserDataOperator(
                "FIELD1", "FIELD1", type),
            mySplit(mySplit) {}

      MoFEMErrorCode doWork(int side, EntityType type,
                            EntitiesFieldData::EntData &data) {
        MoFEMFunctionBeginHot;

        if (data.getFieldData().empty())
          MoFEMFunctionReturnHot(0);

        const double eps = 1e-4;
        for (DoubleAllocator::iterator it = getNormal().data().begin();
             it != getNormal().data().end(); it++) {
          *it = fabs(*it) < eps ? 0.0 : *it;
        }
        for (DoubleAllocator::iterator it =
                 getNormalsAtGaussPtsF3().data().begin();
             it != getNormalsAtGaussPtsF3().data().end(); it++) {
          *it = fabs(*it) < eps ? 0.0 : *it;
        }
        for (DoubleAllocator::iterator it =
                 getTangent1AtGaussPtF3().data().begin();
             it != getTangent1AtGaussPtF3().data().end(); it++) {
          *it = fabs(*it) < eps ? 0.0 : *it;
        }
        for (DoubleAllocator::iterator it =
                 getTangent2AtGaussPtF3().data().begin();
             it != getTangent2AtGaussPtF3().data().end(); it++) {
          *it = fabs(*it) < eps ? 0.0 : *it;
        }
        for (DoubleAllocator::iterator it =
                 getNormalsAtGaussPtsF4().data().begin();
             it != getNormalsAtGaussPtsF4().data().end(); it++) {
          *it = fabs(*it) < eps ? 0.0 : *it;
        }
        for (DoubleAllocator::iterator it =
                 getTangent1AtGaussPtF4().data().begin();
             it != getTangent1AtGaussPtF4().data().end(); it++) {
          *it = fabs(*it) < eps ? 0.0 : *it;
        }
        for (DoubleAllocator::iterator it =
                 getTangent2AtGaussPtF4().data().begin();
             it != getTangent2AtGaussPtF4().data().end(); it++) {
          *it = fabs(*it) < eps ? 0.0 : *it;
        }

        mySplit << "NH1" << std::endl;
        mySplit << "side: " << side << " type: " << type << std::endl;
        mySplit << data << std::endl;
        mySplit << std::setprecision(3) << getCoords() << std::endl;
        mySplit << std::setprecision(3) << getCoordsAtGaussPts() << std::endl;
        mySplit << std::setprecision(3) << getArea(0) << std::endl;
        mySplit << std::setprecision(3) << getArea(1) << std::endl;
        mySplit << std::setprecision(3) << "normal F3 " << getNormalF3()
                << std::endl;
        mySplit << std::setprecision(3) << "normal F4 " << getNormalF4()
                << std::endl;
        mySplit << std::setprecision(3) << "normal at Gauss pt F3 "
                << getNormalsAtGaussPtsF3() << std::endl;
        mySplit << std::setprecision(3) << getTangent1AtGaussPtF3()
                << std::endl;
        mySplit << std::setprecision(3) << getTangent2AtGaussPtF3()
                << std::endl;
        mySplit << std::setprecision(3) << "normal at Gauss pt F4 "
                << getNormalsAtGaussPtsF4() << std::endl;
        mySplit << std::setprecision(3) << getTangent1AtGaussPtF4()
                << std::endl;
        mySplit << std::setprecision(3) << getTangent2AtGaussPtF4()
                << std::endl;
        MoFEMFunctionReturnHot(0);
      }

      MoFEMErrorCode doWork(int row_side, int col_side, EntityType row_type,
                            EntityType col_type,
                            EntitiesFieldData::EntData &row_data,
                            EntitiesFieldData::EntData &col_data) {
        MoFEMFunctionBeginHot;

        if (row_data.getFieldData().empty())
          MoFEMFunctionReturnHot(0);

        mySplit << "NH1NH1" << std::endl;
        mySplit << "row side: " << row_side << " row_type: " << row_type
                << std::endl;
        mySplit << row_data << std::endl;
        mySplit << "NH1NH1" << std::endl;
        mySplit << "col side: " << col_side << " col_type: " << col_type
                << std::endl;
        mySplit << row_data << std::endl;

        MoFEMFunctionReturnHot(0);
      }
    };

    struct MyOp2
        : public FlatPrismElementForcesAndSourcesCore::UserDataOperator {

      TeeStream &mySplit;
      MyOp2(TeeStream &my_split, const char type)
          : FlatPrismElementForcesAndSourcesCore::UserDataOperator(
                "FIELD1", "FIELD2", type),
            mySplit(my_split) {}

      MoFEMErrorCode doWork(int side, EntityType type,
                            EntitiesFieldData::EntData &data) {
        MoFEMFunctionBeginHot;

        if (type != MBENTITYSET)
          MoFEMFunctionReturnHot(0);

        mySplit << "NOFIELD" << std::endl;
        mySplit << "side: " << side << " type: " << type << std::endl;
        mySplit << data << std::endl;
        MoFEMFunctionReturnHot(0);
      }

      MoFEMErrorCode doWork(int row_side, int col_side, EntityType row_type,
                            EntityType col_type,
                            EntitiesFieldData::EntData &row_data,
                            EntitiesFieldData::EntData &col_data) {
        MoFEMFunctionBeginHot;

        unSetSymm();

        if (col_type != MBENTITYSET)
          MoFEMFunctionReturnHot(0);

        mySplit << "NOFILEDH1" << std::endl;
        mySplit << "row side: " << row_side << " row_type: " << row_type
                << std::endl;
        mySplit << row_data << std::endl;
        mySplit << "col side: " << col_side << " col_type: " << col_type
                << std::endl;
        mySplit << col_data << std::endl;

        MoFEMFunctionReturnHot(0);
      }
    };

    FlatPrismElementForcesAndSourcesCore fe1(m_field);
    fe1.getOpPtrVector().push_back(
        new MyOp(my_split, ForcesAndSourcesCore::UserDataOperator::OPROW));
    fe1.getOpPtrVector().push_back(
        new MyOp(my_split, ForcesAndSourcesCore::UserDataOperator::OPROWCOL));
    CHKERR m_field.loop_finite_elements("TEST_PROBLEM", "TEST_FE1", fe1);

    FlatPrismElementForcesAndSourcesCore fe2(m_field);
    fe2.getOpPtrVector().push_back(
        new MyOp2(my_split, ForcesAndSourcesCore::UserDataOperator::OPCOL));
    fe2.getOpPtrVector().push_back(
        new MyOp2(my_split, ForcesAndSourcesCore::UserDataOperator::OPROWCOL));
    CHKERR m_field.loop_finite_elements("TEST_PROBLEM", "TEST_FE2", fe2);
  }
  CATCH_ERRORS;

  CHKERR MoFEM::Core::Finalize();

  return 0;
}
