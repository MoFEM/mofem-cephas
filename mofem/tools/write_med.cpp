/** \file write_med.cpp

  \brief write med files

*/

#include <MoFEM.hpp>
using namespace MoFEM;

static char help[] = "...\n\n";

int main(int argc, char *argv[]) {
  MoFEM::Core::Initialize(&argc, &argv, (char *)0, help);

  try {

    moab::Core mb_instance;
    moab::Interface &moab = mb_instance;
    ParallelComm *pcomm = ParallelComm::get_pcomm(&moab, MYPCOMM_INDEX);
    if (pcomm == NULL)
      pcomm = new ParallelComm(&moab, PETSC_COMM_WORLD);

    char mesh_out_file[255] = "out.h5m";
    char mesh_file_name[255] = "mesh.h5m";
    char meshset_to_write[255] = "";
    PetscBool meshset_to_write_set = PETSC_FALSE;
    int time_step = 0;
    // Create MoFEM database
    MoFEM::Core core(moab, PETSC_COMM_WORLD, VERBOSE);
    MoFEM::Interface &m_field = core;
    PetscBool bit_ref_set = PETSC_FALSE;
    int bit_ref_level = 0;

    CHKERR PetscOptionsBegin(m_field.get_comm(), "", "Read MED tool", "none");
    CHKERR PetscOptionsGetString(PETSC_NULL, "", "-file_name", mesh_file_name,
                                 255, PETSC_NULL);
    CHKERR PetscOptionsInt("-bit_ref_level", "bit ref level", "", bit_ref_level,
                           &bit_ref_level, &bit_ref_set);
    CHKERR PetscOptionsGetString(PETSC_NULL, "", "-meshset_to_write",
                                 meshset_to_write, 255, &meshset_to_write_set);
    CHKERR PetscOptionsInt("-med_time_step", "time step", "", time_step,
                           &time_step, PETSC_NULL);
    CHKERR PetscOptionsString("-output_file", "output mesh file name", "",
                              "out.h5m", mesh_out_file, 255, PETSC_NULL);
    ierr = PetscOptionsEnd();
    CHKERRQ(ierr);

    const char *option = "";

    // load mesh
    CHKERR moab.load_file(mesh_file_name, 0, option);

    // read meshsets
    CHKERR m_field.getInterface<MeshsetsManager>()->readMeshsets();

    // define range to write
    Range entities_to_write;
    if (bit_ref_set) {
      // for writing mesh with crack
      BitRefLevel bit_level;
      // write mesh with crack
      bit_level.set(bit_ref_level);

      // get all entities in bit ref level
      CHKERR m_field.getInterface<BitRefManager>()->getEntitiesByRefLevel(
          bit_level, BitRefLevel().set(), entities_to_write);

      // remove prism entities
      entities_to_write = subtract(entities_to_write,
                                   entities_to_write.subset_by_type(MBPRISM));
      // remove quad entities
      entities_to_write =
          subtract(entities_to_write, entities_to_write.subset_by_type(MBQUAD));
      MOFEM_LOG("MEDWORLD", Sev::inform)
          << "Total number of entities in bit ref level " << bit_ref_level
          << ": " << entities_to_write.size() << std::endl;

      MOFEM_LOG("MEDWORLD", Sev::inform)
          << "Removing interior edges and face entities"<< std::endl;
      // get all subset ranges
      Range entities_3d = entities_to_write.subset_by_dimension(3);
      Range entities_2d = entities_to_write.subset_by_dimension(2);
      Range entities_1d = entities_to_write.subset_by_dimension(1);

      // find skin
      Range skin_2d;
      Range skin_1d;
      Skinner skin(&moab);
      CHKERR skin.find_skin(0, entities_3d, false, skin_2d);
      CHKERR moab.get_adjacencies(skin_2d, 1, false, skin_1d,
                                  moab::Interface::UNION);

      // remove interior entities from entities to write
      entities_2d = subtract(entities_2d, skin_2d);
      entities_to_write = subtract(entities_to_write, entities_2d);

      entities_1d = subtract(entities_1d, skin_1d);
      entities_to_write = subtract(entities_to_write, entities_1d);

      MOFEM_LOG("MEDWORLD", Sev::inform)
          << "Removing nodes with no connectivity"<< std::endl;
      // find nodes with no connectivity
      Range nodes;
      //moab.get_entities_by_type(0, MBVERTEX, nodes, true);
      nodes = entities_to_write.subset_by_type(MBVERTEX);
      // loop to find nodes with no adjacencies
      Range free_nodes;
      for (auto node : nodes) {
        Range adj_edges;
        EntityHandle node_handle = node;
        CHKERR moab.get_adjacencies(&node_handle, 1, 1, false, adj_edges,
                                    moab::Interface::UNION);
        if (adj_edges.empty()) {
          free_nodes.insert(node);
        }
      }
      // get coordinates of nodes with no connectivity
      for (auto node : free_nodes) {
        double coords[3];
        CHKERR moab.get_coords(&node, 1, coords);
      }

      // remove nodes with no connectivity
      entities_to_write = subtract(entities_to_write, free_nodes);

      // Check orientation of crack surfaces
      // Get Meshest 10000 and 10001
      EntityHandle meshset_10000;
      EntityHandle meshset_10001;

      CHKERR m_field.getInterface<MeshsetsManager>()->getMeshset(10000, SIDESET,
                                                                 meshset_10000);
      CHKERR m_field.getInterface<MeshsetsManager>()->getMeshset(10001, SIDESET,
                                                                 meshset_10001);

      Range entities_10000;
      Range entities_10001;
      moab.get_entities_by_handle(meshset_10000, entities_10000);
      moab.get_entities_by_handle(meshset_10001, entities_10001);
      
      auto check_face_orientation = [&moab](const Range &faces) {
        for (auto face : faces) {
          // check orientation of face
          // get adjacent tets
          Range adj_tets;
          CHKERR(moab.get_adjacencies(&face, 1, 3, false, adj_tets,
                                      moab::Interface::UNION));
          // remove other entities from adj_tets
          Range actual_tets;
          actual_tets = adj_tets.subset_by_type(MBTET);

          int side_number;
          int sense;
          int offset;
          for (auto tet : actual_tets) {
            CHKERR(moab.side_number(tet, face, side_number, sense, offset));
            if (sense == -1) {
              // get nodes of face
              std::vector<EntityHandle> nodes_face;
              CHKERR(moab.get_connectivity(&face, 1, nodes_face));
              std::vector<EntityHandle> face_nodes_new;
              // change orientation
              face_nodes_new = {nodes_face[1], nodes_face[0], nodes_face[2]};
              CHKERR(moab.set_connectivity(face, face_nodes_new.data(), 3));
              // recheck orientation
              CHKERR(moab.side_number(tet, face, side_number, sense, offset));
              if (sense == -1) {
                MOFEM_LOG("MEDWORLD", Sev::warning)
                    << "Face: " << face << " has wrong orientation";
              }
            }
          }
        }
      };

      check_face_orientation(entities_10000);
      check_face_orientation(entities_10001);

    } else if (meshset_to_write_set) {
      // loop meshsets to find entities to write
      auto &meshsets_idx =
          m_field.getInterface<MeshsetsManager>()->getMeshsetsMultindex();

      for (auto &meshset : meshsets_idx) {
        auto meshset_name = meshset.getName();
        auto trim_name = [&](std::string &name) {
          name.erase(std::remove(name.begin(), name.end(), ' '), name.end());
        };
        trim_name(meshset_name);

        if (meshset_to_write == meshset_name) {
          CHKERR moab.get_entities_by_handle(meshset.getMeshset(),
                                             entities_to_write, true);
          break;
        }
      }
      MOFEM_LOG("MEDWORLD", Sev::inform)
          << "Wrtiting all entitiies from meshset: " << meshset_to_write
          << std::endl;
    } else {
      // get all entities
      CHKERR moab.get_entities_by_handle(0, entities_to_write, true);
      MOFEM_LOG("MEDWORLD", Sev::inform)
          << "Wrtiting all entitiies" << std::endl;
    }
    MOFEM_LOG("MEDWORLD", Sev::inform)
        << "Number of entities to write: " << entities_to_write.size()
        << std::endl;
    // loop to print size by entity type
    for (EntityType ent_type = MBVERTEX; ent_type < MBMAXTYPE; ent_type++) {
      Range entities;
      moab.get_entities_by_type(0, ent_type, entities, true);
      Range ents_to_write;
      ents_to_write = intersect(entities, entities_to_write);

      if (ents_to_write.empty())
        continue;

      MOFEM_LOG("MEDWORLD", Sev::inform)
          << "Number of entities to write: " << ents_to_write.size()
          << " of type: " << moab::CN::EntityTypeName(ent_type)
          << " from total of: " << entities.size() << std::endl;
    }

    // make shared pointer to range
    boost::shared_ptr<Range> entities_to_write_ptr =
        boost::make_shared<Range>(entities_to_write);

    MedInterface *med_interface_ptr;
    CHKERR m_field.getInterface(med_interface_ptr);
    CHKERR med_interface_ptr->writeMed(entities_to_write_ptr);
  }
  CATCH_ERRORS;

  MoFEM::Core::Finalize();

  return 0;
}
