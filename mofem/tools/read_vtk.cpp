/** \file reading_vtk.cpp

  \brief Reading vtk files

*/



#include <MoFEM.hpp>
using namespace MoFEM;

struct VtkInterface {
    
    VtkInterface(MoFEM::Interface &m_field, moab::Interface &moab) : mField(m_field), mOab(moab) {}

    MoFEMErrorCode readVtk();

private:
    MoFEM::Interface &mField;
    moab::Interface &mOab;
    Simple *simpleInterface;

    MoFEMErrorCode readMesh();
    MoFEMErrorCode getBoundaryConditions();
    MoFEMErrorCode setBoundaryConditions();
    MoFEMErrorCode getMaterialProperties();
    MoFEMErrorCode writeOutput();


};

MoFEMErrorCode VtkInterface::readVtk() {
  MoFEMFunctionBegin;
  CHKERR readMesh();
  CHKERR getBoundaryConditions();
  CHKERR setBoundaryConditions();
  CHKERR getMaterialProperties();
  CHKERR writeOutput();
  MoFEMFunctionReturn(0);
}


MoFEMErrorCode VtkInterface::readMesh() {
  MoFEMFunctionBegin;

  CHKERR mField.getInterface(simpleInterface);
  CHKERR simpleInterface->getOptions();
  CHKERR simpleInterface->loadFile();

  MoFEMFunctionReturn(0);
}

MoFEMErrorCode VtkInterface::getBoundaryConditions() {
  MoFEMFunctionBegin;
  MeshsetsManager *meshsets_manager_ptr;
  CHKERR mField.getInterface(meshsets_manager_ptr);

  Tag bc_tag;
  CHKERR mOab.tag_get_handle("BOUNDARY_CONDITIONS", 1, MB_TYPE_INTEGER, bc_tag,
                             MB_TAG_CREAT | MB_TAG_SPARSE);

  std::vector<std::tuple<int, std::string, std::string>> conditions = {
      {1, "FIX_ALL", "FIX_ALL"},      {2, "FIX_X", "FIX_X"},
      {3, "FIX_Y", "FIX_Y"},          {4, "FIX_Z", "FIX_Z"},
      {5, "FIX_X_2", "DISP_X"},       {6, "FIX_Y_2", "DISP_Y"},
      {7, "FIX_Z_2", "DISP_Z"},       {8, "FORCE_X", "FORCE_X"},
      {9, "FORCE_Y", "FORCE_Y"},      {10, "FORCE_Z", "FORCE_Z"},
      {11, "VEL_X", "VEL_X"},         {12, "VEL_Y", "VEL_Y"},
      {13, "VEL_Z", "VEL_Z"},         {14, "ACCL_X", "ACCL_X"},
      {15, "ACCL_Y", "ACCL_Y"},       {16, "ACCL_Z", "ACCL_Z"},
      {17, "TEMP", "TEMP"},           {18, "PRESSURE", "PRESSURE"},
      {19, "HEAT_FLUX", "HEAT_FLUX"}, {20, "CONTACT", "CONTACT"}};

  for (auto &condition : conditions) {
    int desired_val = std::get<0>(condition);
    const void *desired_val_ptr = &desired_val;
    Range bc;
    CHKERR mOab.get_entities_by_type_and_tag(0, MBVERTEX, &bc_tag,
                                             &desired_val_ptr, 1, bc);

    if (!bc.empty()) {
      Range faces;
      CHKERR mOab.get_adjacencies(bc, 2, true, faces, moab::Interface::UNION);
      Range edges;
      CHKERR mOab.get_adjacencies(bc, 1, true, edges, moab::Interface::UNION);

      Range adj_nodes;
      CHKERR mOab.get_adjacencies(faces, 0, false, adj_nodes,
                                  moab::Interface::UNION);

      Range non_BC_nodes;
      non_BC_nodes = subtract(adj_nodes, bc);
      Range non_BC_faces;
      CHKERR mOab.get_adjacencies(non_BC_nodes, 2, false, non_BC_faces,
                                  moab::Interface::UNION);
      Range non_BC_edges;
      CHKERR mOab.get_adjacencies(non_BC_nodes, 1, false, non_BC_edges,
                                  moab::Interface::UNION);

      faces = subtract(faces, non_BC_faces);
      edges = subtract(edges, non_BC_edges);

      MOFEM_LOG_C("WORLD", Sev::inform, "Number of nodes assigned to %s = %d",
                  std::get<2>(condition).c_str(), bc.size());
      MOFEM_LOG_C("WORLD", Sev::inform, "Number of edges assigned to %s = %d",
                  std::get<2>(condition).c_str(), edges.size());
      MOFEM_LOG_C("WORLD", Sev::inform, "Number of faces assigned to %s = %d",
                  std::get<2>(condition).c_str(), faces.size());

      CHKERR meshsets_manager_ptr->addMeshset(BLOCKSET, desired_val,
                                              std::get<1>(condition));
      CHKERR meshsets_manager_ptr->addEntitiesToMeshset(BLOCKSET, desired_val,
                                                        bc);
      CHKERR meshsets_manager_ptr->addEntitiesToMeshset(BLOCKSET, desired_val,
                                                        edges);
      CHKERR meshsets_manager_ptr->addEntitiesToMeshset(BLOCKSET, desired_val,
                                                        faces);
    }
  }

  MoFEMFunctionReturn(0);
}

MoFEMErrorCode VtkInterface::setBoundaryConditions() {
    MoFEMFunctionBegin;
    // Add meshsets if config file provided
    MeshsetsManager *meshsets_interface_ptr;
    CHKERR mField.getInterface(meshsets_interface_ptr);
    CHKERR meshsets_interface_ptr->setMeshsetFromFile();

    MOFEM_LOG_CHANNEL("WORLD");
    MOFEM_LOG_TAG("WORLD", "read_vtk")
    MOFEM_LOG("WORLD", Sev::inform)
        << "Print all meshsets (old and added from meshsets "
           "configurational file)";
    for (auto cit = meshsets_interface_ptr->getBegin();
         cit != meshsets_interface_ptr->getEnd(); cit++)
      MOFEM_LOG("WORLD", Sev::inform) << *cit;

      MoFEMFunctionReturn(0);
}

MoFEMErrorCode VtkInterface::getMaterialProperties() {
  MoFEMFunctionBegin;
  MeshsetsManager *meshsets_manager_ptr;
  CHKERR mField.getInterface(meshsets_manager_ptr);

  Range ents;
  CHKERR mOab.get_entities_by_type(0, MBHEX, ents, true);
  CHKERR meshsets_manager_ptr->addMeshset(BLOCKSET, 100, "ADOLCMAT");
  CHKERR meshsets_manager_ptr->addEntitiesToMeshset(BLOCKSET, 100, ents);
  MoFEMFunctionReturn(0);
}

MoFEMErrorCode VtkInterface::writeOutput() {
  MoFEMFunctionBegin;
  CHKERR mOab.write_file("out_vtk.h5m");
  MoFEMFunctionReturn(0);
}

static char help[] = "...\n\n";

int main(int argc, char *argv[]) {
  MoFEM::Core::Initialize(&argc, &argv, (char *)0, help);

  try {

    moab::Core mb_instance;
    moab::Interface &moab = mb_instance;
    ParallelComm *pcomm = ParallelComm::get_pcomm(&moab, MYPCOMM_INDEX);
    if (pcomm == NULL)
      pcomm = new ParallelComm(&moab, PETSC_COMM_WORLD);

    // Create MoFEM database
    MoFEM::Core core(moab);
    MoFEM::Interface &m_field = core;

    char mesh_out_file[255] = "out.h5m";

    int time_step = 0;
    CHKERR PetscOptionsBegin(m_field.get_comm(), "", "Read vtk tool", "none");
    CHKERR PetscOptionsString("-output_file", "output mesh file name", "",
                              "out.h5m", mesh_out_file, 255, PETSC_NULL);
    ierr = PetscOptionsEnd();
    CHKERRQ(ierr);

    VtkInterface read_vtk(m_field, moab);
    CHKERR read_vtk.readVtk();




  }
  CATCH_ERRORS;

  MoFEM::Core::Finalize();

  return 0;
}