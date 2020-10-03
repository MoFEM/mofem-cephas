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

typedef tee_device<std::ostream, std::ofstream> TeeDevice;
typedef stream<TeeDevice> TeeStream;

using namespace MoFEM;

static char help[] = "...\n\n";

struct MyOp : public ContactPrismElementForcesAndSourcesCore::UserDataOperator {

  const char faceType;
  MyOp(const char type, const char face_type)
      : ContactPrismElementForcesAndSourcesCore::UserDataOperator(
            "FIELD1", "FIELD1", type, face_type),
        faceType(face_type) {}

  MoFEMErrorCode doWork(int side, EntityType type,
                        DataForcesAndSourcesCore::EntData &data) {
    MoFEMFunctionBeginHot;

    if (data.getFieldData().empty())
      MoFEMFunctionReturnHot(0);

    MOFEM_LOG("ATOM_TEST", Sev::inform) << "NH1";
    MOFEM_LOG("ATOM_TEST", Sev::inform)
        << "side: " << side << " type: " << type;
    MOFEM_LOG("ATOM_TEST", Sev::inform) << data;

    if (faceType == FACEMASTER) {
      MOFEM_LOG("ATOM_TEST", Sev::inform)
          << "coords Master " << std::scientific << std::setprecision(2)
          << getCoordsMaster();
      MOFEM_LOG("ATOM_TEST", Sev::inform)
          << "area Master " << std::scientific << std::setprecision(2)
          << getAreaMaster();
      MOFEM_LOG("ATOM_TEST", Sev::inform)
          << "normal Master " << std::scientific << std::setprecision(2)
          << getNormalMaster();
      MOFEM_LOG("ATOM_TEST", Sev::inform)
          << "coords at Gauss Pts Master " << std::scientific
          << std::setprecision(2) << getCoordsAtGaussPtsMaster();
    } else {
      MOFEM_LOG("ATOM_TEST", Sev::inform)
          << "coords Slave " << std::scientific << std::setprecision(2)
          << getCoordsSlave();
      MOFEM_LOG("ATOM_TEST", Sev::inform)
          << "area Slave " << std::scientific << std::setprecision(2)
          << getAreaSlave();
      MOFEM_LOG("ATOM_TEST", Sev::inform)
          << "normal Slave " << std::scientific << std::setprecision(2)
          << getNormalSlave();
      MOFEM_LOG("ATOM_TEST", Sev::inform)
          << "coords at Gauss Pts Slave " << std::scientific
          << std::setprecision(2) << getCoordsAtGaussPtsSlave();
    }
    MoFEMFunctionReturnHot(0);
  }

  MoFEMErrorCode doWork(int row_side, int col_side, EntityType row_type,
                        EntityType col_type,
                        DataForcesAndSourcesCore::EntData &row_data,
                        DataForcesAndSourcesCore::EntData &col_data) {
    MoFEMFunctionBeginHot;

    if (row_data.getFieldData().empty())
      MoFEMFunctionReturnHot(0);

    MOFEM_LOG("ATOM_TEST", Sev::inform) << "NH1NH1";
    MOFEM_LOG("ATOM_TEST", Sev::inform)
        << "row side: " << row_side << " row_type: " << row_type;
    MOFEM_LOG("ATOM_TEST", Sev::inform) << row_data;
    MOFEM_LOG("ATOM_TEST", Sev::inform) << "NH1NH1";
    MOFEM_LOG("ATOM_TEST", Sev::inform)
        << "col side: " << col_side << " col_type: " << col_type;
    MOFEM_LOG("ATOM_TEST", Sev::inform) << col_data;

    MoFEMFunctionReturnHot(0);
  }
};

struct CallingOp : public ForcesAndSourcesCore::UserDataOperator {

  CallingOp(const char type)
      : ForcesAndSourcesCore::UserDataOperator("FIELD1", "FIELD1", type) {}

  MoFEMErrorCode doWork(int side, EntityType type,
                        DataForcesAndSourcesCore::EntData &data) {
    MoFEMFunctionBeginHot;

    if (data.getFieldData().empty())
      MoFEMFunctionReturnHot(0);

    MOFEM_LOG("ATOM_TEST", Sev::inform) << "Calling Operator NH1";
    MOFEM_LOG("ATOM_TEST", Sev::inform)
        << "side: " << side << " type: " << type;
    MOFEM_LOG("ATOM_TEST", Sev::inform) << data;

    MoFEMFunctionReturnHot(0);
  }

  MoFEMErrorCode doWork(int row_side, int col_side, EntityType row_type,
                        EntityType col_type,
                        DataForcesAndSourcesCore::EntData &row_data,
                        DataForcesAndSourcesCore::EntData &col_data) {
    MoFEMFunctionBeginHot;

    if (row_data.getFieldData().empty())
      MoFEMFunctionReturnHot(0);

    MOFEM_LOG("ATOM_TEST", Sev::inform) << "Calling Operator NH1NH1";
    MOFEM_LOG("ATOM_TEST", Sev::inform)
        << "row side: " << row_side << " row_type: " << row_type;
    MOFEM_LOG("ATOM_TEST", Sev::inform) << row_data;
    MOFEM_LOG("ATOM_TEST", Sev::inform) << "NH1NH1";
    MOFEM_LOG("ATOM_TEST", Sev::inform)
        << "col side: " << col_side << " col_type: " << col_type;
    MOFEM_LOG("ATOM_TEST", Sev::inform) << col_data;

    MoFEMFunctionReturnHot(0);
  }
};

struct MyOp2
    : public ContactPrismElementForcesAndSourcesCore::UserDataOperator {

  MyOp2(const char type, const char face_type)
      : ContactPrismElementForcesAndSourcesCore::UserDataOperator(
            "FIELD1", "FIELD2", type, face_type) {}

  MoFEMErrorCode doWork(int side, EntityType type,
                        DataForcesAndSourcesCore::EntData &data) {
    MoFEMFunctionBeginHot;

    if (type != MBENTITYSET)
      MoFEMFunctionReturnHot(0);

    MOFEM_LOG("ATOM_TEST", Sev::inform) << "NPFIELD";
    MOFEM_LOG("ATOM_TEST", Sev::inform)
        << "side: " << side << " type: " << type;
    MOFEM_LOG("ATOM_TEST", Sev::inform) << data;
    MoFEMFunctionReturnHot(0);
  }

  MoFEMErrorCode doWork(int row_side, int col_side, EntityType row_type,
                        EntityType col_type,
                        DataForcesAndSourcesCore::EntData &row_data,
                        DataForcesAndSourcesCore::EntData &col_data) {
    MoFEMFunctionBeginHot;

    unSetSymm();

    if (col_type != MBENTITYSET)
      MoFEMFunctionReturnHot(0);

    MOFEM_LOG("ATOM_TEST", Sev::inform) << "NOFILEDH1";
    MOFEM_LOG("ATOM_TEST", Sev::inform)
        << "row side: " << row_side << " row_type: " << row_type;
    MOFEM_LOG("ATOM_TEST", Sev::inform) << row_data;
    MOFEM_LOG("ATOM_TEST", Sev::inform)
        << "col side: " << col_side << " col_type: " << col_type;
    MOFEM_LOG("ATOM_TEST", Sev::inform) << col_data;

    MoFEMFunctionReturnHot(0);
  }
};

int main(int argc, char *argv[]) {

  MoFEM::Core::Initialize(&argc, &argv, (char *)0, help);

  try {

    moab::Core mb_instance;
    moab::Interface &moab = mb_instance;
    ParallelComm *pcomm = ParallelComm::get_pcomm(&moab, MYPCOMM_INDEX);
    if (pcomm == NULL)
      pcomm = new ParallelComm(&moab, PETSC_COMM_WORLD);

    PetscBool flg = PETSC_TRUE;
    PetscInt is_hdiv = 0;
    char mesh_file_name[255];
#if PETSC_VERSION_GE(3, 6, 4)
    CHKERR PetscOptionsGetString(PETSC_NULL, "", "-my_file", mesh_file_name,
                                 255, &flg);
    CHKERR PetscOptionsGetInt(PETSC_NULL, "", "-is_hdiv", &is_hdiv, PETSC_NULL);
#else
    CHKERR PetscOptionsGetString(PETSC_NULL, PETSC_NULL, "-my_file",
                                 mesh_file_name, 255, &flg);
    CHKERR PetscOptionsGetInt(PETSC_NULL, PETSC_NULL, "-is_hdiv", &is_hdiv,
                              PETSC_NULL);
#endif

    if (flg != PETSC_TRUE)
      SETERRQ(PETSC_COMM_SELF, 1, "error -my_file (mesh file not given");

    const char *option;
    option = "";
    CHKERR moab.load_file(mesh_file_name, 0, option);

    // Create MoFEM (Joseph) database
    MoFEM::Core core(moab);
    MoFEM::Interface &m_field = core;
    PrismInterface *interface;
    CHKERR m_field.getInterface(interface);

    // set entities bit level
    CHKERR m_field.getInterface<BitRefManager>()->setBitRefLevelByDim(
        0, 3, BitRefLevel().set(0));
    std::vector<BitRefLevel> bit_levels;
    bit_levels.push_back(BitRefLevel().set(0));

    int ll = 1;
    for (_IT_CUBITMESHSETS_BY_SET_TYPE_FOR_LOOP_(m_field, SIDESET, cit)) {

      CHKERR PetscPrintf(PETSC_COMM_WORLD, "Insert Interface %d\n",
                         cit->getMeshsetId());
      EntityHandle cubit_meshset = cit->getMeshset();
      {
        // get tet enties form back bit_level
        EntityHandle ref_level_meshset = 0;
        CHKERR moab.create_meshset(MESHSET_SET, ref_level_meshset);
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
        // get faces and test to split
        CHKERR interface->getSides(cubit_meshset, bit_levels.back(), true, 0);
        // set new bit level
        bit_levels.push_back(BitRefLevel().set(ll++));
        // split faces and
        CHKERR interface->splitSides(ref_level_meshset, bit_levels.back(),
                                     cubit_meshset, true, true, 0);
        // clean meshsets
        CHKERR moab.delete_entities(&ref_level_meshset, 1);
      }
      // update cubit meshsets
      for (_IT_CUBITMESHSETS_FOR_LOOP_(m_field, ciit)) {
        EntityHandle cubit_meshset = ciit->meshset;
        CHKERR m_field.getInterface<BitRefManager>()
            ->updateMeshsetByEntitiesChildren(cubit_meshset, bit_levels.back(),
                                              cubit_meshset, MBVERTEX, true);
        CHKERR m_field.getInterface<BitRefManager>()
            ->updateMeshsetByEntitiesChildren(cubit_meshset, bit_levels.back(),
                                              cubit_meshset, MBEDGE, true);
        CHKERR m_field.getInterface<BitRefManager>()
            ->updateMeshsetByEntitiesChildren(cubit_meshset, bit_levels.back(),
                                              cubit_meshset, MBTRI, true);
        CHKERR m_field.getInterface<BitRefManager>()
            ->updateMeshsetByEntitiesChildren(cubit_meshset, bit_levels.back(),
                                              cubit_meshset, MBTET, true);
      }
    }

    // Fields
    if (is_hdiv)
      CHKERR m_field.add_field("FIELD1", HDIV, DEMKOWICZ_JACOBI_BASE, 3);
    else
      CHKERR m_field.add_field("FIELD1", H1, AINSWORTH_LEGENDRE_BASE, 3);

    CHKERR m_field.add_field("MESH_NODE_POSITIONS", H1, AINSWORTH_LEGENDRE_BASE,
                             3);
    CHKERR m_field.add_field("FIELD2", NOFIELD, NOBASE, 3);

    auto set_no_field_vertex = [&]() {
      MoFEMFunctionBegin;
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
      MoFEMFunctionReturn(0);
    };

    CHKERR set_no_field_vertex();

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
                                                      "TEST_FE1");
    CHKERR m_field.add_ents_to_finite_element_by_type(root_set, MBPRISM,
                                                      "TEST_FE2");

    // set app. order
    constexpr int order = 3;

    if (is_hdiv) {
      CHKERR m_field.set_field_order(root_set, MBTET, "FIELD1", 0);
      CHKERR m_field.set_field_order(root_set, MBTRI, "FIELD1", order);
    } else {
      CHKERR m_field.set_field_order(root_set, MBTET, "FIELD1", order);
      CHKERR m_field.set_field_order(root_set, MBTRI, "FIELD1", order);
      CHKERR m_field.set_field_order(root_set, MBEDGE, "FIELD1", order);
      CHKERR m_field.set_field_order(root_set, MBVERTEX, "FIELD1", 1);
    }

    CHKERR m_field.set_field_order(root_set, MBTET, "MESH_NODE_POSITIONS", 2);
    CHKERR m_field.set_field_order(root_set, MBTRI, "MESH_NODE_POSITIONS", 2);
    CHKERR m_field.set_field_order(root_set, MBEDGE, "MESH_NODE_POSITIONS", 2);
    CHKERR m_field.set_field_order(root_set, MBVERTEX, "MESH_NODE_POSITIONS",
                                   1);

    // build field
    CHKERR m_field.build_fields();

    if (!is_hdiv) {
      // set FIELD1 from positions of 10 node tets
      Projection10NodeCoordsOnField ent_method_field1(m_field, "FIELD1");
      CHKERR m_field.loop_dofs("FIELD1", ent_method_field1);
    }
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

    // mesh partitioning
    // partition
    CHKERR prb_mng_ptr->partitionSimpleProblem("TEST_PROBLEM");
    CHKERR prb_mng_ptr->partitionFiniteElements("TEST_PROBLEM");
    // what are ghost nodes, see Petsc Manual
    CHKERR prb_mng_ptr->partitionGhostDofs("TEST_PROBLEM");

    auto add_atop_logging = [is_hdiv]() {
      auto get_log_file_name = [is_hdiv]() {
        if (is_hdiv)
          return "forces_and_sources_testing_contact_prism_element_HDIV.txt";
        else
          return "forces_and_sources_testing_contact_prism_element.txt";
      };

      auto core_log = logging::core::get();
      core_log->add_sink(
          LogManager::createSink(LogManager::getStrmSelf(), "ATOM_TEST"));
      LogManager::setLog("ATOM_TEST");
      MOFEM_LOG_TAG("ATOM_TEST", "atom test");

      // Add log to file
      logging::add_file_log(keywords::file_name = get_log_file_name(),
                            keywords::channel = "ATOM_TEST");

    };

    add_atop_logging();

    using UMDataOp = ForcesAndSourcesCore::UserDataOperator;
    using ContactDataOp =
        ContactPrismElementForcesAndSourcesCore::UserDataOperator;

    ContactPrismElementForcesAndSourcesCore fe1(m_field);
    fe1.getOpPtrVector().push_back(
        new MyOp(UMDataOp::OPROW, ContactDataOp::FACEMASTER));
    fe1.getOpPtrVector().push_back(
        new MyOp(UMDataOp::OPROW, ContactDataOp::FACESLAVE));
    fe1.getOpPtrVector().push_back(
        new MyOp(UMDataOp::OPROWCOL, ContactDataOp::FACEMASTERMASTER));
    fe1.getOpPtrVector().push_back(
        new MyOp(UMDataOp::OPROWCOL, ContactDataOp::FACEMASTERSLAVE));
    fe1.getOpPtrVector().push_back(
        new MyOp(UMDataOp::OPROWCOL, ContactDataOp::FACESLAVEMASTER));
    fe1.getOpPtrVector().push_back(
        new MyOp(UMDataOp::OPROWCOL, ContactDataOp::FACESLAVESLAVE));
    fe1.getOpPtrVector().push_back(new CallingOp(UMDataOp::OPCOL));
    fe1.getOpPtrVector().push_back(new CallingOp(UMDataOp::OPROW));
    fe1.getOpPtrVector().push_back(new CallingOp(UMDataOp::OPROWCOL));
    CHKERR m_field.loop_finite_elements("TEST_PROBLEM", "TEST_FE1", fe1);

    ContactPrismElementForcesAndSourcesCore fe2(m_field);
    fe2.getOpPtrVector().push_back(
        new MyOp2(UMDataOp::OPCOL, ContactDataOp::FACEMASTER));
    fe2.getOpPtrVector().push_back(
        new MyOp2(UMDataOp::OPCOL, ContactDataOp::FACESLAVE));
    fe2.getOpPtrVector().push_back(
        new MyOp2(UMDataOp::OPROWCOL, ContactDataOp::FACEMASTERMASTER));
    fe2.getOpPtrVector().push_back(
        new MyOp2(UMDataOp::OPROWCOL, ContactDataOp::FACEMASTERSLAVE));
    fe2.getOpPtrVector().push_back(
        new MyOp2(UMDataOp::OPROWCOL, ContactDataOp::FACESLAVEMASTER));
    fe2.getOpPtrVector().push_back(
        new MyOp2(UMDataOp::OPROWCOL, ContactDataOp::FACESLAVESLAVE));
    CHKERR m_field.loop_finite_elements("TEST_PROBLEM", "TEST_FE2", fe2);
  }
  CATCH_ERRORS;

  MoFEM::Core::Finalize();

  return 0;
}
