

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
    CHKERR m_field.add_field("FIELD2", NOFIELD, NOBASE, 3);

    {
      // Creating and adding no field entities.
      const double coords[] = {0, 0, 0};
      EntityHandle no_field_vertex;
      CHKERR m_field.get_moab().create_vertex(coords, no_field_vertex);
      Range range_no_field_vertex;
      range_no_field_vertex.insert(no_field_vertex);
      CHKERR m_field.getInterface<BitRefManager>()->setBitRefLevel(
          range_no_field_vertex, BitRefLevel().set());
      EntityHandle meshset = m_field.get_field_meshset("FIELD2");
      CHKERR m_field.get_moab().add_entities(meshset, range_no_field_vertex);
    }

    // FE
    CHKERR m_field.add_finite_element("TEST_FE1");
    CHKERR m_field.add_finite_element("TEST_FE2");

    // Define rows/cols and element data
    CHKERR m_field.modify_finite_element_add_field_row("TEST_FE1", "FIELD1");
    CHKERR m_field.modify_finite_element_add_field_col("TEST_FE1", "FIELD1");
    CHKERR m_field.modify_finite_element_add_field_data("TEST_FE1", "FIELD1");

    CHKERR m_field.modify_finite_element_add_field_row("TEST_FE2", "FIELD2");
    CHKERR m_field.modify_finite_element_add_field_col("TEST_FE2", "FIELD1");
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
                                                      "TEST_FE1");
    CHKERR m_field.add_ents_to_finite_element_by_type(tets_skin, MBTRI,
                                                      "TEST_FE2");

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

    typedef tee_device<std::ostream, std::ofstream> TeeDevice;
    typedef stream<TeeDevice> TeeStream;

    std::ofstream ofs("forces_and_sources_testing_triangle_element.txt");
    TeeDevice my_tee(std::cout, ofs);
    TeeStream my_split(my_tee);

    struct MyOp1 : public FaceElementForcesAndSourcesCore::UserDataOperator {

      TeeStream &my_split;
      MyOp1(TeeStream &_my_split, const char type)
          : FaceElementForcesAndSourcesCore::UserDataOperator("FIELD1",
                                                              "FIELD1", type),
            my_split(_my_split) {}

      MoFEMErrorCode doWork(int side, EntityType type,
                            EntitiesFieldData::EntData &data) {
        MoFEMFunctionBeginHot;

        const double eps = 1e-4;
        for (DoubleAllocator::iterator it = getNormal().data().begin();
             it != getNormal().data().end(); it++) {
          *it = fabs(*it) < eps ? 0.0 : *it;
        }

        my_split << "NH1" << std::endl;
        my_split << "side: " << side << " type: " << type << std::endl;
        my_split << "data: " << data << std::endl;
        my_split << std::setprecision(3) << getCoords() << std::endl;
        my_split << std::setprecision(3) << getCoordsAtGaussPts() << std::endl;
        my_split << std::setprecision(3) << getArea() << std::endl;
        my_split << std::setprecision(3) << getNormal() << std::endl;
        my_split << std::setprecision(3) << getNormalsAtGaussPts() << std::endl;
        my_split << std::setprecision(3) << getTangent1AtGaussPts()
                 << std::endl;
        my_split << std::setprecision(3) << getTangent2AtGaussPts()
                 << std::endl;
        my_split << std::endl;
        MoFEMFunctionReturnHot(0);
      }

      MoFEMErrorCode doWork(int row_side, int col_side, EntityType row_type,
                            EntityType col_type,
                            EntitiesFieldData::EntData &row_data,
                            EntitiesFieldData::EntData &col_data) {

        MoFEMFunctionBeginHot;
        my_split << "NH1NH1" << std::endl;
        my_split << "row side: " << row_side << " row_type: " << row_type
                 << std::endl;
        my_split << row_data << std::endl;
        my_split << "NH1NH1" << std::endl;
        my_split << "col side: " << col_side << " col_type: " << col_type
                 << std::endl;
        my_split << row_data << std::endl;

        VectorInt row_indices, col_indices;
        CHKERR getProblemRowIndices("FIELD1", row_type, row_side, row_indices);
        CHKERR getProblemColIndices("FIELD1", col_type, col_side, col_indices);

        if (row_indices.size() != row_data.getIndices().size()) {
          SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
                  "row inconsistency");
        }

        if (col_indices.size() != col_data.getIndices().size()) {
          SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
                  "col inconsistency");
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

        my_split << row_data << std::endl;

        MoFEMFunctionReturnHot(0);
      }
    };

    struct MyOp2 : public FaceElementForcesAndSourcesCore::UserDataOperator {

      TeeStream &my_split;
      MyOp2(TeeStream &_my_split, const char type)
          : FaceElementForcesAndSourcesCore::UserDataOperator("FIELD2",
                                                              "FIELD1", type),
            my_split(_my_split) {
        sYmm = false;
       }

      MoFEMErrorCode doWork(int side, EntityType type,
                            EntitiesFieldData::EntData &data) {
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
                            EntitiesFieldData::EntData &row_data,
                            EntitiesFieldData::EntData &col_data) {
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

    FaceElementForcesAndSourcesCore fe1(m_field);
    fe1.getOpPtrVector().push_back(
        new MyOp1(my_split, ForcesAndSourcesCore::UserDataOperator::OPROW));
    fe1.getOpPtrVector().push_back(
        new MyOp1(my_split, ForcesAndSourcesCore::UserDataOperator::OPROWCOL));

    FaceElementForcesAndSourcesCore fe2(m_field);
    fe2.getOpPtrVector().push_back(
        new MyOp2(my_split, ForcesAndSourcesCore::UserDataOperator::OPROW));
    fe2.getOpPtrVector().push_back(
        new MyOp2(my_split, ForcesAndSourcesCore::UserDataOperator::OPROWCOL));

    CHKERR m_field.loop_finite_elements("TEST_PROBLEM", "TEST_FE1", fe1);
    CHKERR m_field.loop_finite_elements("TEST_PROBLEM", "TEST_FE2", fe2);
  }
  CATCH_ERRORS;

  CHKERR MoFEM::Core::Finalize();

  return 0;
}
