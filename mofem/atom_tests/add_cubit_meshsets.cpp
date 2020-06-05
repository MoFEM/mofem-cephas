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

    MOFEM_LOG_CHANNEL("WORLD")
    MOFEM_LOG_ATTRIBUTES("WORLD", LogManager::BitLineID | LogManager::BitScope);

    MOFEM_LOG("WORLD", Sev::verbose) << "<<<< SIDESETs >>>>>";

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
      MOFEM_LOG("WORLD", Sev::inform) << mydata;
    }
    if (!add_block_is_there) 
      SETERRQ(PETSC_COMM_WORLD, MOFEM_OPERATION_UNSUCCESSFUL,
              "no added block set");

    MOFEM_LOG("WORLD", Sev::inform) << "<<<< BLOCKSETs >>>>>";

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
        MOFEM_LOG("WORLD", Sev::inform) << mydata;
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

    CHKERR meshsets_manager_ptr->addMeshset(BLOCKSET, 1002, "ADD_BLOCK_SET");
    CHKERR meshsets_manager_ptr->addMeshset(BLOCKSET, 1003, "ADD_BLOCK_SET");
    CHKERR meshsets_manager_ptr->addMeshset(BLOCKSET, 1004, "ADD_BLOCK_SET");
    CHKERR meshsets_manager_ptr->addMeshset(BLOCKSET, 1005, "ADD_BLOCK_SET");
    CHKERR meshsets_manager_ptr->addMeshset(BLOCKSET, 1006, "ADD_BLOCK_SET");
    CHKERR meshsets_manager_ptr->addMeshset(BLOCKSET, 1007, "ADD_BLOCK_SET");
    CHKERR meshsets_manager_ptr->addMeshset(BLOCKSET, 1008, "ADD_BLOCK_SET");

    MOFEM_LOG("WORLD", Sev::inform)
        << "<<<< ADD BLOCKSETs FROM CONFIG FILE >>>>>";

    CHKERR meshsets_manager_ptr->setMeshsetFromFile(
        /*"add_cubit_meshsets.in"*/);

    // List all meshsets
    for (_IT_CUBITMESHSETS_FOR_LOOP_(m_field, it)) {
      MOFEM_LOG("WORLD", Sev::inform) << *it;
      if ((it->getBcType() & CubitBCType(BLOCKSET)).any()) {
        std::vector<double> attributes;
        it->getAttributes(attributes);
        std::ostringstream ss;
        ss << "Attr: ";
        for (unsigned int ii = 0; ii != attributes.size(); ii++) {
          ss << attributes[ii] << " ";
        }
        MOFEM_LOG("WORLD", Sev::inform) << ss.str(); 
      }
      if ((it->getBcType() & CubitBCType(MAT_ELASTICSET)).any()) {
        Mat_Elastic mydata;
        CHKERR it->getAttributeDataStructure(mydata);
        MOFEM_LOG("WORLD", Sev::inform) << "Mat elastic found ";
        MOFEM_LOG("WORLD", Sev::inform) << mydata;
      }
      if ((it->getBcType() & CubitBCType(MAT_THERMALSET)).any()) {
        Mat_Thermal mydata;
        CHKERR it->getAttributeDataStructure(mydata);
        MOFEM_LOG("WORLD", Sev::inform) << "Mat thermal found ";
        MOFEM_LOG("WORLD", Sev::inform) << mydata;
      }
      if ((it->getBcType() & CubitBCType(DISPLACEMENTSET)).any()) {
        DisplacementCubitBcData mydata;
        CHKERR it->getBcDataStructure(mydata);
        MOFEM_LOG("WORLD", Sev::inform) << mydata;
      }
      if ((it->getBcType() & CubitBCType(FORCESET)).any()) {
        ForceCubitBcData mydata;
        CHKERR it->getBcDataStructure(mydata);
        MOFEM_LOG("WORLD", Sev::inform) << mydata;
      }
      if ((it->getBcType() & CubitBCType(PRESSURESET)).any()) {
        PressureCubitBcData mydata;
        CHKERR it->getBcDataStructure(mydata);
        MOFEM_LOG("WORLD", Sev::inform) << mydata;
      }
      if ((it->getBcType() & CubitBCType(TEMPERATURESET)).any()) {
        TemperatureCubitBcData mydata;
        CHKERR it->getBcDataStructure(mydata);
        MOFEM_LOG("WORLD", Sev::inform) << mydata;
      }
      if ((it->getBcType() & CubitBCType(HEATFLUXSET)).any()) {
        HeatFluxCubitBcData mydata;
        CHKERR it->getBcDataStructure(mydata);
        MOFEM_LOG("WORLD", Sev::inform) << mydata;
      }
    }
  }
  CATCH_ERRORS;

  MoFEM::Core::Finalize();
}
