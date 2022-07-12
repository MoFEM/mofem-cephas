/** \file forces_and_sources_testing_edge_element.cpp
 * \example forces_and_sources_testing_edge_element.cpp
 * \brief Testing edge elements
 * nodesets.
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

    // Create MoFEM database
    MoFEM::Core core(moab);
    MoFEM::Interface &m_field = core;

    const char *option;
    option = ""; ;
    CHKERR moab.load_file(mesh_file_name, 0, option);

    // set ebturues bit level
    BitRefLevel bit_level0;
    bit_level0.set(0);
    CHKERR m_field.getInterface<BitRefManager>()->setBitRefLevelByDim(
        0, 3, bit_level0);

    enum bases { AINSWORTH, BERNSTEIN_BEZIER, LASBASETOP };
    const char *list_bases[] = {"ainsworth", "bernstein_bezier"};
    PetscInt choice_base_value = AINSWORTH;
    CHKERR PetscOptionsGetEList(PETSC_NULL, NULL, "-base", list_bases,
                                LASBASETOP, &choice_base_value, &flg);
    FieldApproximationBase base = NOBASE;
    if (choice_base_value == AINSWORTH) 
      base = AINSWORTH_LEGENDRE_BASE;
    else if (choice_base_value == BERNSTEIN_BEZIER) 
      base = AINSWORTH_BERNSTEIN_BEZIER_BASE;

    // Fields
    CHKERR m_field.add_field("FIELD1", H1, base, 3);
    CHKERR m_field.add_field("MESH_NODE_POSITIONS", H1, AINSWORTH_LEGENDRE_BASE,
                             3);

    // FE
    CHKERR m_field.add_finite_element("TEST_FE");
    

    // Define rows/cols and element data
    auto add_field_to_fe = [&m_field](const std::string field_name) {
      MoFEMFunctionBegin;
      CHKERR m_field.modify_finite_element_add_field_row("TEST_FE", field_name);
      CHKERR m_field.modify_finite_element_add_field_col("TEST_FE", field_name);
      CHKERR m_field.modify_finite_element_add_field_data("TEST_FE",
                                                          field_name);
      MoFEMFunctionReturn(0);
    };
    CHKERR add_field_to_fe("FIELD1");
    CHKERR m_field.modify_finite_element_add_field_data("TEST_FE",
                                                        "MESH_NODE_POSITIONS");

    // Problem
    CHKERR m_field.add_problem("TEST_PROBLEM");

    // set finite elements for problem
    CHKERR m_field.modify_problem_add_finite_element("TEST_PROBLEM", "TEST_FE");
    // set refinement level for problem
    CHKERR m_field.modify_problem_ref_level_add_bit("TEST_PROBLEM", bit_level0);
    

    // meshset consisting all entities in mesh
    EntityHandle root_set = moab.get_root_set();
    // add entities to field
    CHKERR m_field.add_ents_to_field_by_type(root_set, MBEDGE, "FIELD1");
    CHKERR m_field.add_ents_to_field_by_type(root_set, MBEDGE,
                                             "MESH_NODE_POSITIONS");
    CHKERR m_field.set_field_order(root_set, MBVERTEX, "MESH_NODE_POSITIONS",
                                   1);
    CHKERR m_field.set_field_order(root_set, MBEDGE, "MESH_NODE_POSITIONS", 1);
    CHKERR m_field.set_field_order(root_set, MBTRI, "MESH_NODE_POSITIONS", 1);
    CHKERR m_field.set_field_order(root_set, MBTET, "MESH_NODE_POSITIONS", 1);
    

    // add entities to finite element
    Range tets;
    CHKERR moab.get_entities_by_type(0, MBTET, tets, false);
    Skinner skin(&m_field.get_moab());
    Range tets_skin;
    CHKERR skin.find_skin(0, tets, false, tets_skin);
    Range tets_skin_edges;
    CHKERR moab.get_adjacencies(tets_skin, 1, false, tets_skin_edges,
                                moab::Interface::UNION);
    CHKERR m_field.add_ents_to_finite_element_by_type(tets_skin_edges, MBEDGE,
                                                      "TEST_FE");
    
    // set app. order
    // see Hierarchic Finite Element Bases on Unstructured Tetrahedral Meshes
    // (Mark Ainsworth & Joe Coyle)
    int order = 3;
    CHKERR m_field.set_field_order(root_set, MBEDGE, "FIELD1", order);
    CHKERR m_field.set_field_order(root_set, MBVERTEX, "FIELD1", 1);

    /****/
    // build database
    // build field
    CHKERR m_field.build_fields();
    
    // set FIELD1 from positions of 10 node tets
    if (base == AINSWORTH_LEGENDRE_BASE) {
      Projection10NodeCoordsOnField ent_method(m_field, "FIELD1");
      CHKERR m_field.loop_dofs("FIELD1", ent_method);
      
    }
    {
      Projection10NodeCoordsOnField ent_method(m_field, "MESH_NODE_POSITIONS");
      CHKERR m_field.loop_dofs("MESH_NODE_POSITIONS", ent_method);
      
    }
    // build finite elemnts
    CHKERR m_field.build_finite_elements();
    
    // build adjacencies
    CHKERR m_field.build_adjacencies(bit_level0);
    
    // build problem
    ProblemsManager *prb_mng_ptr;
    CHKERR m_field.getInterface(prb_mng_ptr);
    
    CHKERR prb_mng_ptr->buildProblem("TEST_PROBLEM", true);
    
    // mesh partitioning
    // partition
    CHKERR prb_mng_ptr->partitionSimpleProblem("TEST_PROBLEM");
    CHKERR prb_mng_ptr->partitionFiniteElements("TEST_PROBLEM");
    
    // what are ghost nodes, see Petsc Manual
    CHKERR prb_mng_ptr->partitionGhostDofs("TEST_PROBLEM");

    EdgeElementForcesAndSourcesCore fe1(m_field);

    typedef tee_device<std::ostream, std::ofstream> TeeDevice;
    typedef stream<TeeDevice> TeeStream;

    std::ofstream ofs(
        (std ::string("forces_and_sources_testing_edge_element_") +
         ApproximationBaseNames[base] + ".txt")
            .c_str());
    TeeDevice my_tee(std::cout, ofs);
    TeeStream my_split(my_tee);

    struct MyOp : public EdgeElementForcesAndSourcesCore::UserDataOperator {

      TeeStream &my_split;
      MyOp(TeeStream &_my_split, const char type)
          : EdgeElementForcesAndSourcesCore::UserDataOperator("FIELD1",
                                                              "FIELD1", type),
            my_split(_my_split) {}

      MoFEMErrorCode doWork(int side, EntityType type,
                            EntitiesFieldData::EntData &data) {
        MoFEMFunctionBegin;

        my_split << "NH1" << std::endl;
        my_split << "side: " << side << " type: " << type << std::endl;
        my_split << "data: " << data << std::endl;
        my_split << "coords: " << std::setprecision(3) << getCoords()
                 << std::endl;
        my_split << "getCoordsAtGaussPts: " << std::setprecision(3)
                 << getCoordsAtGaussPts() << std::endl;
        my_split << "length: " << std::setprecision(3) << getLength()
                 << std::endl;
        my_split << "direction: " << std::setprecision(3) << getDirection()
                 << std::endl;

        int nb_gauss_pts = data.getN().size1();
        for (int gg = 0; gg < nb_gauss_pts; gg++) {
          my_split << "tangent " << gg << " " << getTangetAtGaussPts()
                   << std::endl;
        }

        MoFEMFunctionReturn(0);
      }

      MoFEMErrorCode doWork(int row_side, int col_side, EntityType row_type,
                            EntityType col_type,
                            EntitiesFieldData::EntData &row_data,
                            EntitiesFieldData::EntData &col_data) {
        MoFEMFunctionBegin;
        my_split << "ROW NH1NH1" << std::endl;
        my_split << "row side: " << row_side << " row_type: " << row_type
                 << std::endl;
        my_split << row_data << std::endl;
        my_split << "COL NH1NH1" << std::endl;
        my_split << "col side: " << col_side << " col_type: " << col_type
                 << std::endl;
        my_split << col_data << std::endl;
        MoFEMFunctionReturn(0);
      }
    };

    fe1.getOpPtrVector().push_back(
        new OpGetHOTangentsOnEdge("MESH_NODE_POSITIONS"));
    fe1.getOpPtrVector().push_back(
        new MyOp(my_split, ForcesAndSourcesCore::UserDataOperator::OPROW));
    fe1.getOpPtrVector().push_back(
        new MyOp(my_split, ForcesAndSourcesCore::UserDataOperator::OPROWCOL));

    CHKERR m_field.loop_finite_elements("TEST_PROBLEM", "TEST_FE", fe1);
    }
  CATCH_ERRORS;

  MoFEM::Core::Finalize();

  return 0;
}
