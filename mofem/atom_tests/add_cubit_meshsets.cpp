/**
 * @file add_cubit_meshsets.cpp
 * @example add_cubit_meshsets.cpp
 * @brief Test and example setting cubit meshsets
 *
 */

#include <MoFEM.hpp>

using namespace MoFEM;

static char help[] = "...\n\n";

int main(int argc, char *argv[]) {

  MoFEM::Core::Initialize(&argc, &argv, (char *)0, help);

  auto core_log = logging::core::get();
  core_log->add_sink(
      LogManager::createSink(LogManager::getStrmSelf(), "ATOM_TEST"));
  LogManager::setLog("ATOM_TEST");

  try {

    moab::Core mb_instance;
    moab::Interface &moab = mb_instance;
    int rank;
    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

    // Create MoFEM (Joseph) database
    MoFEM::Core core(moab);
    MoFEM::Interface &m_field = core;

    MeshsetsManager *meshsets_manager_ptr;
    CHKERR m_field.getInterface(meshsets_manager_ptr);

    MOFEM_LOG_CHANNEL("ATOM_TEST")
    MOFEM_LOG_ATTRIBUTES("ATOM_TEST",
                         LogManager::BitLineID | LogManager::BitScope);
    MOFEM_LOG_TAG("ATOM_TEST", "atom test");

    MOFEM_LOG("ATOM_TEST", Sev::verbose) << "<<<< SIDESETs >>>>>";

    bool add_block_is_there = false;
    CHKERR meshsets_manager_ptr->addMeshset(SIDESET, 1002);
    {
      PressureCubitBcData mybc;
      strncpy(mybc.data.name, "Pressure", 8);
      mybc.data.flag1 = 0;
      mybc.data.flag2 = 0;
      mybc.data.value1 = 1;
      CHKERR meshsets_manager_ptr->setBcData(SIDESET, 1002, mybc);
    }
    for (_IT_CUBITMESHSETS_BY_SET_TYPE_FOR_LOOP_(m_field, SIDESET, it)) {
      if (it->getMeshsetId() != 1002)
        continue;
      add_block_is_there = true;
      PressureCubitBcData mydata;
      CHKERR it->getBcDataStructure(mydata);
      // Print data
      MOFEM_LOG("ATOM_TEST", Sev::inform) << mydata;
    }
    if (!add_block_is_there)
      SETERRQ(PETSC_COMM_WORLD, MOFEM_OPERATION_UNSUCCESSFUL,
              "no added block set");

    MOFEM_LOG("ATOM_TEST", Sev::inform) << "<<<< BLOCKSETs >>>>>";

    add_block_is_there = false;
    CHKERR meshsets_manager_ptr->addMeshset(BLOCKSET, 1000, "ADD_BLOCK_SET");
    std::vector<double> attr(3);
    attr[0] = 0;
    attr[1] = 1;
    attr[2] = 2;
    CHKERR meshsets_manager_ptr->setAtributes(BLOCKSET, 1000, attr);
    for (_IT_CUBITMESHSETS_BY_SET_TYPE_FOR_LOOP_(m_field, BLOCKSET, it)) {
      // Get block name
      std::string name = it->getName();
      if (name.compare(0, 13, "ADD_BLOCK_SET") == 0) {
        add_block_is_there = true;
        std::vector<double> attributes;
        it->getAttributes(attributes);
        if (attributes.size() != 3) {
          SETERRQ1(PETSC_COMM_WORLD, MOFEM_ATOM_TEST_INVALID,
                   "should be 3 attributes but is %d", attributes.size());
        }
        if (attributes[0] != 0 || attributes[1] != 1 || attributes[2] != 2) {
          SETERRQ(PETSC_COMM_WORLD, MOFEM_ATOM_TEST_INVALID,
                  "wrong values of attributes");
        }
      }
    }
    if (!add_block_is_there) {
      SETERRQ(PETSC_COMM_WORLD, MOFEM_OPERATION_UNSUCCESSFUL,
              "no added block set");
    }
    add_block_is_there = false;
    CHKERR meshsets_manager_ptr->addMeshset(BLOCKSET, 1001, "MAT_ELASTIC");
    {
      Mat_Elastic mydata;
      mydata.data.Young = 1;
      mydata.data.Poisson = 0.25;
      CHKERR meshsets_manager_ptr->setAtributesByDataStructure(BLOCKSET, 1001,
                                                               mydata);
    }
    for (_IT_CUBITMESHSETS_BY_SET_TYPE_FOR_LOOP_(m_field, BLOCKSET, it)) {
      if (it->getMeshsetId() != 1001)
        continue;
      // Get block name
      std::string name = it->getName();
      if (name.compare(0, 13, "MAT_ELASTIC") == 0 &&
          (it->getBcType() & CubitBCType(MAT_ELASTICSET)).any()) {
        add_block_is_there = true;
        Mat_Elastic mydata;
        CHKERR it->getAttributeDataStructure(mydata);
        // Print data
        MOFEM_LOG("ATOM_TEST", Sev::inform) << mydata;
        if (mydata.data.Young != 1 || mydata.data.Poisson != 0.25) {
          SETERRQ(PETSC_COMM_WORLD, MOFEM_ATOM_TEST_INVALID,
                  "wrong values of attributes");
        }
      }
    }
    if (!add_block_is_there) {
      SETERRQ(PETSC_COMM_WORLD, MOFEM_OPERATION_UNSUCCESSFUL,
              "no added block set");
    }

    MOFEM_LOG("ATOM_TEST", Sev::inform) << "<<<< NODESET >>>>>";

    CHKERR meshsets_manager_ptr->addMeshset(NODESET, 1010);
    DisplacementCubitBcData disp_bc;
    std::memcpy(disp_bc.data.name, "Displacement", 12);
    disp_bc.data.flag1 = 1;
    disp_bc.data.flag2 = 1;
    disp_bc.data.flag3 = 1;
    disp_bc.data.flag4 = 0;
    disp_bc.data.flag5 = 0;
    disp_bc.data.flag6 = 0;
    disp_bc.data.value1 = 0;
    disp_bc.data.value2 = 0;
    disp_bc.data.value3 = 0;
    disp_bc.data.value4 = 0;
    disp_bc.data.value5 = 0;
    disp_bc.data.value6 = 0;

    CHKERR meshsets_manager_ptr->setBcData(NODESET, 1010, disp_bc);

    for (_IT_CUBITMESHSETS_BY_BCDATA_TYPE_FOR_LOOP_(
             m_field, NODESET | DISPLACEMENTSET, it)) {
      DisplacementCubitBcData disp_data;
      CHKERR it->getBcDataStructure(disp_data);
      MOFEM_LOG("ATOM_TEST", Sev::inform) << disp_data;
    }

    MOFEM_LOG("ATOM_TEST", Sev::inform)
        << "<<<< ADD BLOCKSETs FROM CONFIG FILE >>>>>";

    CHKERR meshsets_manager_ptr->addMeshset(BLOCKSET, 1002, "ADD_BLOCK_SET");
    CHKERR meshsets_manager_ptr->addMeshset(BLOCKSET, 1003, "ADD_BLOCK_SET");
    CHKERR meshsets_manager_ptr->addMeshset(BLOCKSET, 1004, "ADD_BLOCK_SET");
    CHKERR meshsets_manager_ptr->addMeshset(BLOCKSET, 1005, "ADD_BLOCK_SET");
    CHKERR meshsets_manager_ptr->addMeshset(BLOCKSET, 1006, "ADD_BLOCK_SET");
    CHKERR meshsets_manager_ptr->addMeshset(BLOCKSET, 1007, "ADD_BLOCK_SET");
    CHKERR meshsets_manager_ptr->addMeshset(BLOCKSET, 1008, "ADD_BLOCK_SET");
    CHKERR meshsets_manager_ptr->addMeshset(BLOCKSET, 1009, "ADD_BLOCK_SET");

    CHKERR meshsets_manager_ptr->setMeshsetFromFile(
        /*"add_cubit_meshsets.in"*/);

    // List all meshsets
    MOFEM_LOG("ATOM_TEST", Sev::inform) << "Iterate blocksets";

    bool mat_elastic_trans_is_found = true;
    for (_IT_CUBITMESHSETS_FOR_LOOP_(m_field, it)) {
      MOFEM_LOG("ATOM_TEST", Sev::inform) << *it;
      if ((it->getBcType() & CubitBCType(BLOCKSET)).any()) {
        std::vector<double> attributes;
        it->getAttributes(attributes);
        std::ostringstream ss;
        ss << "Attr: ";
        for (unsigned int ii = 0; ii != attributes.size(); ii++) {
          ss << attributes[ii] << " ";
        }
        MOFEM_LOG("ATOM_TEST", Sev::inform) << ss.str();

        std::string block_name = it->getName();
        if (block_name.compare(0, block_name.size(), "MAT_ELASTIC_TRANS_ISO") ==
            0) {
          MOFEM_LOG("ATOM_TEST", Sev::inform) << "Mat Trans Iso block ";
          mat_elastic_trans_is_found = true;
          Mat_Elastic_TransIso mydata;
          CHKERR it->getAttributeDataStructure(mydata);
          // Print data
          MOFEM_LOG("ATOM_TEST", Sev::inform) << mydata;
          if (mydata.data.Youngp != 1)
            SETERRQ(PETSC_COMM_SELF, MOFEM_ATOM_TEST_INVALID, "Wrong value");
          if (mydata.data.Youngz != 2)
            SETERRQ(PETSC_COMM_SELF, MOFEM_ATOM_TEST_INVALID, "Wrong value");
          if (mydata.data.Poissonp != 3)
            SETERRQ(PETSC_COMM_SELF, MOFEM_ATOM_TEST_INVALID, "Wrong value");
          if (mydata.data.Poissonpz != 4)
            SETERRQ(PETSC_COMM_SELF, MOFEM_ATOM_TEST_INVALID, "Wrong value");
          if (mydata.data.Shearzp != 5)
            SETERRQ(PETSC_COMM_SELF, MOFEM_ATOM_TEST_INVALID, "Wrong value");
        } 
      }
      if(!mat_elastic_trans_is_found)
        SETERRQ(PETSC_COMM_SELF, MOFEM_ATOM_TEST_INVALID, "Block not found");

      if ((it->getBcType() & CubitBCType(MAT_ELASTICSET)).any()) {
        Mat_Elastic mydata;
        CHKERR it->getAttributeDataStructure(mydata);
        MOFEM_LOG("ATOM_TEST", Sev::inform) << "Mat elastic found ";
        MOFEM_LOG("ATOM_TEST", Sev::inform) << mydata;
      }
      if ((it->getBcType() & CubitBCType(MAT_THERMALSET)).any()) {
        Mat_Thermal mydata;
        CHKERR it->getAttributeDataStructure(mydata);
        MOFEM_LOG("ATOM_TEST", Sev::inform) << "Mat thermal found ";
        MOFEM_LOG("ATOM_TEST", Sev::inform) << mydata;
      }
      if ((it->getBcType() & CubitBCType(DISPLACEMENTSET)).any()) {
        DisplacementCubitBcData mydata;
        CHKERR it->getBcDataStructure(mydata);
        MOFEM_LOG("ATOM_TEST", Sev::inform) << mydata;
      }
      if ((it->getBcType() & CubitBCType(FORCESET)).any()) {
        ForceCubitBcData mydata;
        CHKERR it->getBcDataStructure(mydata);
        MOFEM_LOG("ATOM_TEST", Sev::inform) << mydata;
      }
      if ((it->getBcType() & CubitBCType(PRESSURESET)).any()) {
        PressureCubitBcData mydata;
        CHKERR it->getBcDataStructure(mydata);
        MOFEM_LOG("ATOM_TEST", Sev::inform) << mydata;
      }
      if ((it->getBcType() & CubitBCType(TEMPERATURESET)).any()) {
        TemperatureCubitBcData mydata;
        CHKERR it->getBcDataStructure(mydata);
        MOFEM_LOG("ATOM_TEST", Sev::inform) << mydata;
      }
      if ((it->getBcType() & CubitBCType(HEATFLUXSET)).any()) {
        HeatFluxCubitBcData mydata;
        CHKERR it->getBcDataStructure(mydata);
        MOFEM_LOG("ATOM_TEST", Sev::inform) << mydata;
      }
    }
  }
  CATCH_ERRORS;

  MoFEM::Core::Finalize();
}
