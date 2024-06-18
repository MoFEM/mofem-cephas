/** \file MedInterface.cpp
 * \brief Med file interface interface
 *
 * Interface loading mesh and data on mesh directly to mofem & moab
 *
 */

#ifdef WITH_MED

extern "C" {
#include <med.h>
}

#if (MED_MAJOR_NUM >= 3)
// To avoid too many ifdefs below we use defines for the bits of the
// API that did not change too much between MED2 and MED3. If we
// remove MED2 support at some point, please remove these defines and
// replace the symbols accordingly.
#define med_geometrie_element med_geometry_type
#define med_entite_maillage med_entity_type
#define med_type_champ med_field_type
#define MED_LECTURE MED_ACC_RDONLY
#define MED_LECTURE_AJOUT MED_ACC_RDEXT
#define MED_NOEUD MED_NODE
#define MED_MAILLE MED_CELL
#define MED_NOEUD_MAILLE MED_NODE_ELEMENT
#define MED_NOPFL MED_NO_PROFILE
#define MEDouvrir MEDfileOpen
#define MEDfermer MEDfileClose
#define MEDnChamp MEDfieldnComponent
#define MEDnValProfil MEDprofileSizeByName
#define MEDfichDesEcr MEDfileCommentWr
#define MED_ARETE MED_EDGE
#define MED_FACETTE MED_FACE
#else
#error "MED file is not MED_MAJOR_NUM == 3"
#endif

#include <MedInterface.hpp>
#include <unordered_set>

namespace MoFEM {

MoFEMErrorCode
MedInterface::query_interface(boost::typeindex::type_index type_index,
                              UnknownInterface **iface) const {
  *iface = const_cast<MedInterface *>(this);
  return 0;
}

MedInterface::MedInterface(const Core &core)
    : cOre(const_cast<Core &>(core)), flgFile(PETSC_FALSE) {

  if (!LogManager::checkIfChannelExist("MEDWORLD")) {
    auto core_log = logging::core::get();
    core_log->add_sink(
        LogManager::createSink(LogManager::getStrmWorld(), "MEDWORLD"));
    LogManager::setLog("MEDWORLD");
    MOFEM_LOG_TAG("MEDWORLD", "MED");
  }

  MOFEM_LOG("MEDWORLD", Sev::noisy) << "Mashset manager created";
}

MoFEMErrorCode MedInterface::getFileNameFromCommandLine(int verb) {
  Interface &m_field = cOre;
  char mesh_file_name[255];
  MoFEMFunctionBeginHot;
  ierr = PetscOptionsBegin(m_field.get_comm(), "", "MED Interface", "none");
  CHKERRG(ierr);
  ierr = PetscOptionsString("-med_file", "med file name", "", "mesh.med",
                            mesh_file_name, 255, &flgFile);
  CHKERRG(ierr);
  ierr = PetscOptionsEnd();
  CHKERRG(ierr);
  if (flgFile == PETSC_TRUE) {
    medFileName = std::string(mesh_file_name);
  } else {
    medFileName = std::string("mesh.med");
  }
  MoFEMFunctionReturnHot(0);
}

MoFEMErrorCode MedInterface::medGetFieldNames(const string &file, int verb) {
  Interface &m_field = cOre;
  MoFEMFunctionBegin;
  med_idt fid = MEDfileOpen(file.c_str(), MED_ACC_RDONLY);
  if (fid < 0) {
    SETERRQ1(m_field.get_comm(), MOFEM_OPERATION_UNSUCCESSFUL,
             "Unable to open file '%s'", file.c_str());
  }
  med_int num_fields = MEDnField(fid);
  for (int index = 0; index < num_fields; index++) {

    med_int num_comp = MEDfieldnComponent(fid, index + 1);
    if (num_comp <= 0) {
      SETERRQ(m_field.get_comm(), MOFEM_IMPOSSIBLE_CASE,
              "Could not get number of components for MED field");
    }

    char name[MED_NAME_SIZE + 1], mesh_name[MED_NAME_SIZE + 1];
    char dt_unit[MED_SNAME_SIZE + 1];
    std::vector<char> comp_name(num_comp * MED_SNAME_SIZE + 1);
    std::vector<char> comp_unit(num_comp * MED_SNAME_SIZE + 1);
    med_int num_steps = 0;
    med_type_champ type;
    med_bool local_mesh;
    if (MEDfieldInfo(fid, index + 1, name, mesh_name, &local_mesh, &type,
                     &comp_name[0], &comp_unit[0], dt_unit, &num_steps) < 0) {
      SETERRQ(m_field.get_comm(), MOFEM_OPERATION_UNSUCCESSFUL,
              "Could not get MED field info");
    }

    std::string field_name = std::string(name);
    fieldNames[field_name] = FieldData();
    fieldNames[field_name].fieldName = field_name;
    fieldNames[field_name].meshName = std::string(mesh_name);
    fieldNames[field_name].localMesh = (local_mesh == MED_TRUE);
    for (int ff = 0; ff != num_comp; ff++) {
      fieldNames[field_name].componentNames.push_back(
          std::string(&comp_name[ff * MED_SNAME_SIZE], MED_SNAME_SIZE));
      fieldNames[field_name].componentUnits.push_back(
          std::string(&comp_unit[ff * MED_SNAME_SIZE], MED_SNAME_SIZE));
    }
    fieldNames[field_name].dtUnit = std::string(&dt_unit[0]);
    fieldNames[field_name].ncSteps = num_steps;

    if (verb > 0)
      MOFEM_LOG("MEDWORLD", Sev::inform) << fieldNames[name];
  }
  if (MEDfileClose(fid) < 0) {
    SETERRQ1(m_field.get_comm(), MOFEM_OPERATION_UNSUCCESSFUL,
             "Unable to close file '%s'", file.c_str());
  }
  MoFEMFunctionReturn(0);
}

MoFEMErrorCode MedInterface::medGetFieldNames(int verb) {
  MoFEMFunctionBegin;
  if (medFileName.empty()) {
    CHKERR getFileNameFromCommandLine(verb);
  }
  CHKERR medGetFieldNames(medFileName, verb);
  MoFEMFunctionReturn(0);
}

MoFEMErrorCode MedInterface::readMed(const string &file, int verb) {
  Interface &m_field = cOre;
  MoFEMFunctionBegin;

  med_idt fid = MEDfileOpen(file.c_str(), MED_ACC_RDONLY);
  if (fid < 0) {
    SETERRQ1(m_field.get_comm(), MOFEM_OPERATION_UNSUCCESSFUL,
             "Unable to open file '%s'", file.c_str());
  }

  med_int v[3], vf[3];
  MEDlibraryNumVersion(&v[0], &v[1], &v[2]);
  MEDfileNumVersionRd(fid, &vf[0], &vf[1], &vf[2]);

  if (verb > 0)
    MOFEM_LOG_C("MEDWORLD", Sev::inform,
                "Reading MED file V%d.%d.%d using MED library V%d.%d.%d", vf[0],
                vf[1], vf[2], v[0], v[1], v[2]);

  if (vf[0] < 2 || (vf[0] == 2 && vf[1] < 2)) {
    SETERRQ(m_field.get_comm(), MOFEM_DATA_INCONSISTENCY,
            "Cannot read MED file older than V2.2");
  }

  for (int i = 0; i < MEDnMesh(fid); i++) {
    char mesh_name[MED_NAME_SIZE + 1], mesh_desc[MED_COMMENT_SIZE + 1];
    bzero(mesh_name, MED_NAME_SIZE);
    bzero(mesh_desc, MED_COMMENT_SIZE);
    med_int space_dim;
    med_mesh_type mesh_type;
    med_int mesh_dim, n_step;
    char dt_unit[MED_SNAME_SIZE + 1];
    char axis_name[3 * MED_SNAME_SIZE + 1], axis_unit[3 * MED_SNAME_SIZE + 1];
    med_sorting_type sorting_type;
    med_axis_type axis_type;
    if (MEDmeshInfo(fid, i + 1, mesh_name, &space_dim, &mesh_dim, &mesh_type,
                    mesh_desc, dt_unit, &sorting_type, &n_step, &axis_type,
                    axis_name, axis_unit) < 0) {
      SETERRQ(m_field.get_comm(), MOFEM_OPERATION_UNSUCCESSFUL,
              "Unable to read mesh information");
    }
    meshNames.push_back(std::string(mesh_name));
    if (verb > 0) {
      MOFEM_LOG_C("MEDWORLD", Sev::inform, "Check mesh %s nsteps %d", mesh_name,
                  n_step);
    }
  }

  std::map<int, Range> family_elem_map;
  std::map<string, Range> group_elem_map;

  for (unsigned int ii = 0; ii != meshNames.size(); ii++) {
    CHKERR readMesh(file, ii, family_elem_map, verb);
    CHKERR readFamily(file, ii, family_elem_map, group_elem_map, verb);
    CHKERR makeBlockSets(group_elem_map, verb);
  }

  if (MEDfileClose(fid) < 0) {
    SETERRQ1(m_field.get_comm(), MOFEM_OPERATION_UNSUCCESSFUL,
             "Unable to close file '%s'", file.c_str());
  }

  MoFEMFunctionReturn(0);
}

static std::vector<med_geometrie_element>
moab2med_element_type(const EntityType type) {

  std::vector<med_geometrie_element> types;

  switch (type) {
  case MBEDGE:
    types.push_back(MED_SEG2);
    types.push_back(MED_SEG3);
    break;
  case MBTRI:
    types.push_back(MED_TRIA3);
    types.push_back(MED_TRIA6);
    break;
  case MBQUAD:
    types.push_back(MED_QUAD4);
    break;
  case MBTET:
    types.push_back(MED_TETRA4);
    types.push_back(MED_TETRA10);
    break;
  case MBHEX:
    types.push_back(MED_HEXA8);
    break;
  case MBPRISM:
    types.push_back(MED_PENTA6);
    break;
  case MBPYRAMID:
    types.push_back(MED_PYRA5);
    break;
  case MBVERTEX:
    types.push_back(MED_POINT1);
  default:
    break;
  }

  return types;
}

// static EntityType
// med2moab_element_type(std::vector<med_geometrie_element> med_type) {

//   EntityType type;

//   switch (med_type) {
//   case MED_SEG2 || MED_SEG3:
//     type = MBEDGE;
//     break;
//   case MED_TRIA3 || MED_TRIA6:
//     type = MBTRI;
//     break;
//   case MED_QUAD4:
//     type = MBQUAD;
//     break;
//   case MED_TETRA4 || MED_TETRA10:
//     type = MBTET;
//     break;
//   case MED_HEXA8:
//     type = MBHEX;
//     break;
//   case MED_PENTA6: 
//     type = MBPRISM;
//     break;
//   case MED_PYRA5:
//     type = MBPYRAMID;
//     break;
//   case MED_POINT1:
//     type = MBVERTEX;
//   default:
//     break;
//   }
//   return type;
// }

MoFEMErrorCode MedInterface::readMesh(const string &file, const int index,
                                      std::map<int, Range> &family_elem_map,
                                      int verb) {

  Interface &m_field = cOre;
  MoFEMFunctionBegin;

  med_idt fid = MEDfileOpen(file.c_str(), MED_ACC_RDONLY);
  if (fid < 0) {
    SETERRQ1(m_field.get_comm(), MOFEM_OPERATION_UNSUCCESSFUL,
             "Unable to open file '%s'", file.c_str());
  }
  med_int v[3], vf[3];
  MEDlibraryNumVersion(&v[0], &v[1], &v[2]);
  MEDfileNumVersionRd(fid, &vf[0], &vf[1], &vf[2]);
  if (verb > 1)
    MOFEM_LOG_C("MEDWORLD", Sev::noisy,
                "Reading MED file V%d.%d.%d using MED library V%d.%d.%d", vf[0],
                vf[1], vf[2], v[0], v[1], v[2]);

  if (vf[0] < 2 || (vf[0] == 2 && vf[1] < 2)) {
    SETERRQ(m_field.get_comm(), MOFEM_DATA_INCONSISTENCY,
            "Cannot read MED file older than V2.2");
  }

  char mesh_name[MED_NAME_SIZE + 1], mesh_desc[MED_COMMENT_SIZE + 1];
  bzero(mesh_name, MED_NAME_SIZE + 1);
  bzero(mesh_desc, MED_COMMENT_SIZE + 1);
  med_int space_dim;
  med_mesh_type mesh_type;
  med_int mesh_dim, n_step;
  char dt_unit[MED_SNAME_SIZE + 1];
  char axis_name[3 * MED_SNAME_SIZE + 1], axis_unit[3 * MED_SNAME_SIZE + 1];
  med_sorting_type sorting_type;
  med_axis_type axis_type;
  if (MEDmeshInfo(fid, index + 1, mesh_name, &space_dim, &mesh_dim, &mesh_type,
                  mesh_desc, dt_unit, &sorting_type, &n_step, &axis_type,
                  axis_name, axis_unit) < 0) {
    SETERRQ(m_field.get_comm(), MOFEM_OPERATION_UNSUCCESSFUL,
            "Unable to read mesh information");
  }
  if (verb > 0)
    MOFEM_LOG_C("MEDWORLD", Sev::inform, "Reading mesh %s nsteps %d", mesh_name,
                n_step);

  switch (axis_type) {
  case MED_CARTESIAN:
    break;
  case MED_SPHERICAL:
  case MED_CYLINDRICAL:
  default:
    SETERRQ(m_field.get_comm(), MOFEM_NOT_IMPLEMENTED,
            "Curvilinear coordinate system implemented");
  }
  if (space_dim < 2) {
    SETERRQ1(m_field.get_comm(), MOFEM_NOT_IMPLEMENTED,
             "Not implemented for space dim %d", space_dim);
  }

  EntityHandle mesh_meshset;
  {
    MeshsetsManager *meshsets_manager_ptr;
    CHKERR m_field.getInterface(meshsets_manager_ptr);
    int max_id = 0;
    for (_IT_CUBITMESHSETS_BY_SET_TYPE_FOR_LOOP_(m_field, BLOCKSET, cit)) {
      max_id = (max_id < cit->getMeshsetId()) ? cit->getMeshsetId() : max_id;
    }
    max_id++;
    CHKERR meshsets_manager_ptr->addMeshset(BLOCKSET, max_id,
                                            std::string(mesh_name));
    CubitMeshSet_multiIndex::index<
        Composite_Cubit_msId_And_MeshsetType_mi_tag>::type::iterator cit;
    cit =
        meshsets_manager_ptr->getMeshsetsMultindex()
            .get<Composite_Cubit_msId_And_MeshsetType_mi_tag>()
            .find(boost::make_tuple(max_id, CubitBCType(BLOCKSET).to_ulong()));
    max_id++;
    mesh_meshset = cit->getMeshset();
    meshMeshsets.push_back(mesh_meshset);
  }

  med_bool change_of_coord, geo_transform;
  med_int num_nodes = MEDmeshnEntity(
      fid, mesh_name, MED_NO_DT, MED_NO_IT, MED_NODE, MED_NO_GEOTYPE,
      MED_COORDINATE, MED_NO_CMODE, &change_of_coord, &geo_transform);
  if (num_nodes < 0) {
    SETERRQ(m_field.get_comm(), MOFEM_OPERATION_UNSUCCESSFUL,
            "Could not read number of MED nodes");
  }
  if (num_nodes == 0) {
    SETERRQ(m_field.get_comm(), MOFEM_OPERATION_UNSUCCESSFUL,
            "No nodes in MED mesh");
  }
  if (verb > 0)
    MOFEM_LOG_C("MEDWORLD", Sev::inform, "Read number of nodes %d", num_nodes);

  std::vector<med_float> coord_med(space_dim * num_nodes);
  if (MEDmeshNodeCoordinateRd(fid, mesh_name, MED_NO_DT, MED_NO_IT,
                              MED_NO_INTERLACE, &coord_med[0]) < 0) {
    SETERRQ(m_field.get_comm(), MOFEM_OPERATION_UNSUCCESSFUL,
            "Could not read MED node coordinates");
  }

  // Add vertices to moab
  ReadUtilIface *iface;
  vector<double *> arrays_coord;
  EntityHandle startv;
  CHKERR m_field.get_moab().query_interface(iface);
  CHKERR iface->get_node_coords(3, num_nodes, 0, startv, arrays_coord);
  Range verts(startv, startv + num_nodes - 1);
  std::copy(&coord_med[0 * num_nodes], &coord_med[1 * num_nodes],
            arrays_coord[0]);
  std::copy(&coord_med[1 * num_nodes], &coord_med[2 * num_nodes],
            arrays_coord[1]);
  if (space_dim == 2) {
    std::fill(arrays_coord[2], &arrays_coord[2][num_nodes], 0.);
  } else {
    std::copy(&coord_med[2 * num_nodes], &coord_med[3 * num_nodes],
              arrays_coord[2]);
  }
  CHKERR m_field.get_moab().add_entities(mesh_meshset, verts);
  family_elem_map.clear();

  // get family for vertices
  {
    std::vector<med_int> nodes_tags(num_nodes, 0);
    if (MEDmeshEntityFamilyNumberRd(fid, mesh_name, MED_NO_DT, MED_NO_IT,
                                    MED_NODE, MED_NO_GEOTYPE,
                                    &nodes_tags[0]) < 0) {
      nodes_tags.clear();
      // SETERRQ(
      //   m_field.get_comm(),
      //   MOFEM_OPERATION_UNSUCCESSFUL,
      //   "No family number for elements: using 0 as default family number"
      // );
    }
    for (int i = 0; i < num_nodes; i++) {
      // cerr << verts[i] << " " /*<< ele_tags[j] << " "*/ << nodes_tags[i] <<
      // endl;
      family_elem_map[nodes_tags.empty() ? i : nodes_tags[i]].insert(verts[i]);
    }
  }

  // read elements (loop over all possible MSH element types)
  for (EntityType ent_type = MBVERTEX; ent_type < MBMAXTYPE; ent_type++) {

    auto types = moab2med_element_type(ent_type);

    for (auto type : types) {

      // get number of cells
      med_bool change_of_coord;
      med_bool geo_transform;
      med_int num_ele = MEDmeshnEntity(
          fid, mesh_name, MED_NO_DT, MED_NO_IT, MED_CELL, type,
          MED_CONNECTIVITY, MED_NODAL, &change_of_coord, &geo_transform);

      // get connectivity
      if (num_ele <= 0)
        continue;

      int num_nod_per_ele = type % 100;
      if (verb > 0)
        MOFEM_LOG("MEDWORLD", Sev::inform)
            << "Reading elements " << num_ele << " of type "
            << moab::CN::EntityTypeName(ent_type) << " number of nodes "
            << num_nod_per_ele;

      std::vector<med_int> conn_med(num_ele * num_nod_per_ele);
      if (MEDmeshElementConnectivityRd(fid, mesh_name, MED_NO_DT, MED_NO_IT,
                                       MED_CELL, type, MED_NODAL,
                                       MED_FULL_INTERLACE, &conn_med[0]) < 0) {
        SETERRQ(m_field.get_comm(), MOFEM_OPERATION_UNSUCCESSFUL,
                "Could not read MED elements");
      }

      // cerr << "type " << ent_type << " ";
      // cerr << "num_ele " << num_ele << " " << num_nod_per_ele << endl;;

      Range ents;

      if (ent_type != MBVERTEX) {
        EntityHandle *conn_moab;
        EntityHandle starte;
        CHKERR iface->get_element_connect(num_ele, num_nod_per_ele, ent_type, 0,
                                          starte, conn_moab);
        switch (ent_type) {
        // FIXME: Some connectivity could not work, need to run and test
        case MBTET: {
          int ii = 0;
          for (int ee = 0; ee != num_ele; ee++) {
            EntityHandle n[4];
            for (int nn = 0; nn != num_nod_per_ele; nn++) {
              // conn_moab[ii] = verts[conn_med[ii]-1];
              n[nn] = verts[conn_med[ii + nn] - 1];
            }

            std::array<int, 10> nodes_map{

                1, 0, 2, 3,

                4, 6, 5, 8,
                7, 9

            };

            for (int nn = 0; nn != num_nod_per_ele; nn++, ii++) {
              conn_moab[ii] = n[nodes_map[nn]];
            }
          }
        } break;
        default: {
          int ii = 0;
          for (int ee = 0; ee != num_ele; ee++) {
            for (int nn = 0; nn != num_nod_per_ele; nn++, ii++) {
              // cerr << conn_med[ii] << " ";
              conn_moab[ii] = verts[conn_med[ii] - 1];
            }
            // cerr << endl;
          }
        }
        }
        CHKERR iface->update_adjacencies(starte, num_ele, num_nod_per_ele,
                                         conn_moab);
        ents = Range(starte, starte + num_ele - 1);
      } else {
        // This is special case when in med vertices are defined as elements
        int ii = 0;
        std::vector<EntityHandle> conn_moab(num_ele * num_nod_per_ele);
        for (int ee = 0; ee != num_ele; ++ee)
          for (int nn = 0; nn != num_nod_per_ele; ++nn, ++ii)
            conn_moab[ii] = verts[conn_med[ii] - 1];
        ents.insert_list(conn_moab.begin(), conn_moab.end());
      }

      // Add elements to family meshset
      CHKERR m_field.get_moab().add_entities(mesh_meshset, ents);

      // get family for cells
      {
        std::vector<med_int> fam(num_ele, 0);
        if (MEDmeshEntityFamilyNumberRd(fid, mesh_name, MED_NO_DT, MED_NO_IT,
                                        MED_CELL, type, &fam[0]) < 0) {
          SETERRQ(m_field.get_comm(), MOFEM_OPERATION_UNSUCCESSFUL,
                  "No family number for elements: using 0 as default family "
                  "number");
        }
        // std::vector<med_int> ele_tags(num_ele);
        // if(MEDmeshEntityNumberRd(
        //   fid, mesh_name,MED_NO_DT,MED_NO_IT,MED_CELL,type,&ele_tags[0]) < 0
        // ) {
        //   ele_tags.clear();
        // }
        for (int j = 0; j < num_ele; j++) {
          // cerr << ents[j] << " " /*<< ele_tags[j] << " "*/ << fam[j] << endl;
          family_elem_map[fam[j]].insert(ents[j]);
        }
      }
    }
  }

  if (MEDfileClose(fid) < 0) {
    SETERRQ1(m_field.get_comm(), MOFEM_OPERATION_UNSUCCESSFUL,
             "Unable to close file '%s'", file.c_str());
  }

  MoFEMFunctionReturn(0);
}

MoFEMErrorCode
MedInterface::readFamily(const string &file, const int index,
                         const std::map<int, Range> &family_elem_map,
                         std::map<string, Range> &group_elem_map, int verb) {
  //
  Interface &m_field = cOre;
  MoFEMFunctionBegin;

  med_idt fid = MEDfileOpen(file.c_str(), MED_ACC_RDONLY);
  if (fid < 0) {
    SETERRQ1(m_field.get_comm(), MOFEM_OPERATION_UNSUCCESSFUL,
             "Unable to open file '%s'", file.c_str());
  }
  med_int v[3], vf[3];
  MEDlibraryNumVersion(&v[0], &v[1], &v[2]);
  MEDfileNumVersionRd(fid, &vf[0], &vf[1], &vf[2]);

  if (verb > 1) {
    MOFEM_LOG_C("MEDWORLD", Sev::noisy,
                "Reading MED file V%d.%d.%d using MED library V%d.%d.%d", vf[0],
                vf[1], vf[2], v[0], v[1], v[2]);
  }
  if (vf[0] < 2 || (vf[0] == 2 && vf[1] < 2)) {
    SETERRQ(m_field.get_comm(), MOFEM_DATA_INCONSISTENCY,
            "Cannot read MED file older than V2.2");
  }

  // clear groups
  group_elem_map.clear();

  med_int num_families = MEDnFamily(fid, meshNames[index].c_str());
  if (num_families < 0) {
    SETERRQ(m_field.get_comm(), MOFEM_OPERATION_UNSUCCESSFUL,
            "Could not read MED families");
  }
  for (int i = 0; i < num_families; i++) {
    med_int num_attrib =
        (vf[0] == 2)
            ? MEDnFamily23Attribute(fid, meshNames[index].c_str(), i + 1)
            : 0;
    med_int num_groups = MEDnFamilyGroup(fid, meshNames[index].c_str(), i + 1);
    if (num_attrib < 0 || num_groups < 0) {
      SETERRQ(m_field.get_comm(), MOFEM_OPERATION_UNSUCCESSFUL,
              "Could not read MED groups or attributes");
    }

    std::vector<med_int> attribId(num_attrib + 1);
    std::vector<med_int> attrib_val(num_attrib + 1);
    std::vector<char> attrib_des(MED_COMMENT_SIZE * num_attrib + 1);
    std::vector<char> group_names(MED_LNAME_SIZE * num_groups + 1);
    char family_name[MED_NAME_SIZE + 1];
    med_int family_num;

    if (vf[0] == 2) { // MED2 file
      if (MEDfamily23Info(fid, meshNames[index].c_str(), i + 1, family_name,
                          &attribId[0], &attrib_val[0], &attrib_des[0],
                          &family_num, &group_names[0]) < 0) {
        SETERRQ1(m_field.get_comm(), MOFEM_OPERATION_UNSUCCESSFUL,
                 "Could not read info for MED2 family %d", i + 1);
      }
    } else {
      if (MEDfamilyInfo(fid, meshNames[index].c_str(), i + 1, family_name,
                        &family_num, &group_names[0]) < 0) {
        SETERRQ1(m_field.get_comm(), MOFEM_OPERATION_UNSUCCESSFUL,
                 "Could not read info for MED3 family %d", i + 1);
      }
    }

    // cerr << family_name << " " << family_num  << " " << num_groups << endl;
    for (int g = 0; g != num_groups; g++) {
      std::string name =
          std::string(&group_names[MED_LNAME_SIZE * g], MED_LNAME_SIZE - 1);
      name.resize(NAME_TAG_SIZE - 1);
      if (family_elem_map.find(family_num) == family_elem_map.end()) {
        MOFEM_LOG_C(
            "MEDWORLD", Sev::warning,
            "Family %d not read, likely type of element is not added "
            "to moab database. Currently only triangle, quad, tetrahedral and "
            "hexahedral elements are read to moab database",
            family_num);
      } else {
        group_elem_map[name].merge(family_elem_map.at(family_num));
        // cerr << string(&group_names[MED_LNAME_SIZE*g]) << endl;
      }
    }
  }

  if (MEDfileClose(fid) < 0) {
    SETERRQ1(m_field.get_comm(), MOFEM_OPERATION_UNSUCCESSFUL,
             "Unable to close file '%s'", file.c_str());
  }

  MoFEMFunctionReturn(0);
}

MoFEMErrorCode
MedInterface::makeBlockSets(const std::map<string, Range> &group_elem_map,
                            int verb) {

  Interface &m_field = cOre;
  MoFEMFunctionBegin;
  MeshsetsManager *meshsets_manager_ptr;
  CHKERR m_field.getInterface(meshsets_manager_ptr);

  int max_id = 0;
  for (_IT_CUBITMESHSETS_BY_SET_TYPE_FOR_LOOP_(m_field, BLOCKSET, cit)) {
    max_id = (max_id < cit->getMeshsetId()) ? cit->getMeshsetId() : max_id;
  }
  max_id++;

  for (std::map<string, Range>::const_iterator git = group_elem_map.begin();
       git != group_elem_map.end(); git++) {
    CHKERR meshsets_manager_ptr->addMeshset(BLOCKSET, max_id, git->first);
    CubitMeshSet_multiIndex::index<
        Composite_Cubit_msId_And_MeshsetType_mi_tag>::type::iterator cit;
    cit =
        meshsets_manager_ptr->getMeshsetsMultindex()
            .get<Composite_Cubit_msId_And_MeshsetType_mi_tag>()
            .find(boost::make_tuple(max_id, CubitBCType(BLOCKSET).to_ulong()));
    EntityHandle meshsets = cit->getMeshset();
    if (!git->second.empty()) {
      CHKERR m_field.get_moab().add_entities(meshsets, git->second);
    }
    max_id++;
  }

  for (_IT_CUBITMESHSETS_BY_SET_TYPE_FOR_LOOP_(m_field, BLOCKSET, cit))
    MOFEM_LOG("MEDWORLD", Sev::verbose) << *cit;

  MoFEMFunctionReturn(0);
}

MoFEMErrorCode MedInterface::readMed(int verb) {
  MoFEMFunctionBegin;
  if (medFileName.empty()) {
    CHKERR getFileNameFromCommandLine(verb);
  }
  CHKERR readMed(medFileName, verb);
  MoFEMFunctionReturn(0);
}

MoFEMErrorCode MedInterface::getMeshsets(
    boost::shared_ptr<std::vector<const CubitMeshSets *>> &meshsets_ptr) {
  Interface &m_field = cOre;
  MoFEMFunctionBegin;

  auto &meshsets_idx =
      m_field.getInterface<MeshsetsManager>()->getMeshsetsMultindex();
  std::vector<const CubitMeshSets *> meshsets;

  for (auto &m : meshsets_idx) {
    meshsets.push_back(&m);
  }
    // Sort meshsets based on meshsetId
  std::sort(meshsets.begin(), meshsets.end(),
      [](const CubitMeshSets* a, const CubitMeshSets* b) {
          return a->getMeshsetId() < b->getMeshsetId();
      }
  );

  meshsets_ptr = boost::make_shared<std::vector<const CubitMeshSets *>>(meshsets);
  MoFEMFunctionReturn(0);
}

MoFEMErrorCode MedInterface::writeMed(boost::shared_ptr<Range> range_ptr, int verb) {
  MoFEMFunctionBegin;
  if (medFileName.empty()) {
    CHKERR getFileNameFromCommandLine(verb);
  }

  // Get all meshsets
  boost::shared_ptr<std::vector<const CubitMeshSets *>> meshsets_ptr;
  getMeshsets(meshsets_ptr);

  CHKERR writeMed(medFileName, meshsets_ptr, range_ptr, verb);
  MoFEMFunctionReturn(0);
}

MoFEMErrorCode MedInterface::writeMed(
    const string &file,
    boost::shared_ptr<std::vector<const CubitMeshSets *>> meshsets_ptr,
    boost::shared_ptr<Range> write_range_ptr, int verb) {
  Interface &m_field = cOre;
  // Get the MOAB instance from the field
  moab::Interface &moab = m_field.get_moab();

  MoFEMFunctionBeginHot;
  MOFEM_LOG("MEDWORLD", Sev::warning) << "WRITE_MED IS EXPERIMENTAL";
  // Open a med file with the specified version
  med_idt fid =
      MEDfileVersionOpen((char *)file.c_str(), MED_ACC_CREAT, MED_MAJOR_NUM,
                         MED_MINOR_NUM, MED_RELEASE_NUM);
  // Throw an error if the file could not be opened
  if (fid < 0) {
    SETERRQ(m_field.get_comm(), MOFEM_OPERATION_UNSUCCESSFUL,
            "Cannot create MED file");
    return 0;
  }

  // Write appropriate header
  CHKERR MEDfichDesEcr(fid, (char *)"MED file generated by MoFEM");

  // Generate empty mesh
  char dtUnit[MED_SNAME_SIZE + 1] = "";       // time unit
  char axisName[3 * MED_SNAME_SIZE + 1] = ""; // axis name
  char axisUnit[3 * MED_SNAME_SIZE + 1] = ""; // axis unit

  PetscBool is_cubit_meshset = PETSC_TRUE;
  int med_mesh_name_id = 0;

  // Create the mesh
  char mesh_name_char[255];
  std::string mesh_name = "Mesh";

  // Get mesh_name from command line 
  CHKERR PetscOptionsBegin(PETSC_COMM_WORLD, "", "MED mesh options", "");
  CHKERR PetscOptionsString("-med_mesh_name", "get med mesh name", "",
                     mesh_name.c_str(), mesh_name_char, 255, PETSC_NULL);
  ierr = PetscOptionsEnd();
  CHKERRG(ierr);

  mesh_name = mesh_name_char;
  
  // check if mesh_name is a meshset
  for (auto &m : *meshsets_ptr) {
    if (m->getName() == mesh_name) {
      med_mesh_name_id = m->getMeshsetId();
      is_cubit_meshset = PETSC_FALSE;
      break;
    }
  }

  // Get the maximum meshset id
  int max_id = 0;
  for (auto &m : *meshsets_ptr) {
    if (m->getMeshsetId() == med_mesh_name_id)
      mesh_name = m->getName();
    max_id = (max_id < m->getMeshsetId()) ? m->getMeshsetId() : max_id;
  }

  // Create a mesh
    CHKERR MEDmeshCr(fid, mesh_name.c_str(), 3, 3, MED_UNSTRUCTURED_MESH,
                     "Mesh created", dtUnit, MED_SORT_DTIT, MED_CARTESIAN,
                     axisName, axisUnit);

    // Create a map of meshset names to family ids and shared meshsets
    med_int family_id = 0;
    std::map<std::vector<std::string>, std::tuple<med_int, std::vector<int>>> shared_meshsets_map;
    std::map<EntityHandle, med_int> entityHandle_family_map;

    // initialise zero family id
    shared_meshsets_map[std::vector<string>()] =
        std::make_tuple(family_id, std::vector<int>());

    // function to get meshset names 
    auto get_set_name = [&](const CubitMeshSets *iit) {
      if (iit->getBcTypeULong() & BLOCKSET) {
        EntityHandle meshset = iit->getMeshset();

        std::string name = iit->getName();
        if (name == "NoNameSet") {
          name = "BLOCKSET_NoNameSet_";
          name += std::to_string(iit->getMeshsetId());
        }
        return name;
      } else if (iit->getBcTypeULong() & SIDESET ||
                 iit->getBcTypeULong() & NODESET) {
        EntityHandle meshset = iit->getMeshset();

        CubitBCType cubitBcType(iit->getBcTypeULong());
        auto test = iit->getBcType();

        // std::cout << iit->getBcTypeULong() << std::endl;

        unsigned jj = 0;
        string bc_type_name;
        while (1 << jj != LASTSET_BC) {
          const CubitBCType jj_bc_type = 1 << jj;
          if ((iit->getBcType() & jj_bc_type).any()) {

            bc_type_name += string(CubitBCNames[jj + 1]);
            bc_type_name += "_";
          }
          ++jj;
        }
        bc_type_name += std::to_string(iit->getMeshsetId());
        return bc_type_name;
      }
    };

    // loop over all entities in the write range
    for (auto &entity : *write_range_ptr) {
      // check if entity is shared with other meshsets
      std::vector<int> shared_meshsets;
      std::vector<string> shared_names;

      for (auto &other_meshset : *meshsets_ptr) {
        if (med_mesh_name_id == other_meshset->getMeshsetId())
          continue;

        Range other_entities;
        EntityHandle other_set = other_meshset->getMeshset();
        moab.get_entities_by_handle(other_set, other_entities);

        //   get entity type
        EntityType ent_type = moab.type_from_handle(entity);

        bool is_in_meshset = moab.contains_entities(other_set, &entity, 1);
        if (is_in_meshset) {
          if (ent_type == MBVERTEX) {
            // add shared meshset id to list
            // check if higher dimension entity is in meshset i.e not a nodeset
            bool is_in_higher_dim = false;
            Range entities_in_higher_dim;
            for (int i = 1; i < 4; i++) {
              moab.get_entities_by_dimension(other_set, i,
                                             entities_in_higher_dim);
              if (!entities_in_higher_dim.empty()) {
                is_in_higher_dim = true;
                break;
              }
            }
            if (is_in_higher_dim) {
              continue;
            }
            // check if name is already added to shared_names
            if (std::find(shared_names.begin(), shared_names.end(),
                          get_set_name(other_meshset)) == shared_names.end()) {
              shared_meshsets.push_back(other_meshset->getMeshsetId());
              shared_names.push_back(get_set_name(other_meshset));
            }
          } else {
            // add shared meshset id to list
            if (std::find(shared_names.begin(), shared_names.end(),
                          get_set_name(other_meshset)) == shared_names.end()) {
              shared_meshsets.push_back(other_meshset->getMeshsetId());
              shared_names.push_back(get_set_name(other_meshset));
            }
          }
        }
      }
      // check if shared meshset is already in map
      if (shared_meshsets_map.find(shared_names) != shared_meshsets_map.end()) {
        // assign family id to entity
      } else {
        // create new family id
        family_id++;
        std::tuple<med_int, std::vector<int>> family_tuple;
        family_tuple = std::make_tuple(family_id, shared_meshsets);
        shared_meshsets_map[shared_names] = family_tuple;
      }
      // assign family id to entity
      entityHandle_family_map[entity] =
          std::get<0>(shared_meshsets_map[shared_names]);
    }

    // loop to create families based on shared meshsets map
    for (auto &it : shared_meshsets_map) {
      // create family
      std::string family_name = "F_";
      std::tuple<med_int, std::vector<int>> family_tuple = it.second;
      family_name += std::to_string(std::get<0>(family_tuple));
      std::vector<std::string> shared_meshset_names = it.first;
      std::string group_name;
      for (auto &name : shared_meshset_names) {
        // get meshset name
        std::string meshset_name = name;
        meshset_name.resize(MED_LNAME_SIZE, ' ');
        group_name += meshset_name;
      }
      // create family
      CHKERR MEDfamilyCr(fid, mesh_name.c_str(), family_name.c_str(),
                         std::get<0>(family_tuple), shared_meshset_names.size(),
                         group_name.c_str());
      MOFEM_LOG("MEDWORLD", Sev::inform)
          << "Creating family " << family_name << " with id "
          << std::get<0>(family_tuple) << " and " << shared_meshset_names.size()
          << " groups with name " << group_name << std::endl;
    }

    // write nodes
    Range verts;
    moab.get_entities_by_type(0, MBVERTEX, verts);
    // Prepare arrays to hold the coordinates
    std::vector<med_float> coord_med(3 * verts.size());
    std::vector<med_int> fam;
    std::vector<med_int> tags;

    // For each vertex, get its coordinates
    for (Range::iterator it = verts.begin(); it != verts.end(); ++it) {
      double coords[3];
      moab.get_coords(&(*it), 1, coords);
      coord_med[3 * (it - verts.begin())] = coords[0];
      coord_med[3 * (it - verts.begin()) + 1] = coords[1];
      coord_med[3 * (it - verts.begin()) + 2] = coords[2];
      fam.push_back(entityHandle_family_map[*it]);
      tags.push_back(*it);
    }
    // Write the coordinates to the MED file
    CHKERR MEDmeshNodeWr(fid, mesh_name.c_str(), MED_NO_DT, MED_NO_IT,
                         MED_UNDEF_DT, MED_FULL_INTERLACE, verts.size(),
                         &coord_med[0], MED_FALSE, "", MED_TRUE, &tags[0],
                         MED_TRUE, &fam[0]);

    // loop to write elements to families
    for (EntityType ent_type = MBVERTEX; ent_type < MBMAXTYPE; ent_type++) {
      // get all entities of type ent_type
      Range entities;
      moab.get_entities_by_type(0, ent_type, entities, true);
      Range ents_to_write;
      ents_to_write = intersect(entities, *write_range_ptr);

      if (ents_to_write.empty())
        continue;

      // vectors to write
      std::vector<med_int> tag_number;
      std::vector<med_int> family_number;
      std::vector<med_int> connectivity;

      // loop over all entities
      for (auto &entity : ents_to_write) {
        if (ent_type != MBVERTEX) {
          // get family id for entity
          family_number.push_back(entityHandle_family_map[entity]);
          // get tag number for entity
          tag_number.push_back(entity);
          // get connectivity for entity
          std::vector<EntityHandle> conn;
          moab.get_connectivity(&entity, 1, conn);
          for (auto &c : conn) {
            connectivity.push_back(c);
          }
        } else {
          continue;
        }
      }

      auto get_med_element_type = [](EntityType ent_type) {
        med_geometrie_element type;
        switch (ent_type) {
        case MBHEX:
          type = MED_HEXA8;
          break;
        case MBTET:
          type = MED_TETRA4;
          break;
        case MBQUAD:
          type = MED_QUAD4;
          break;
        case MBTRI:
          type = MED_TRIA3;
          break;
        case MBEDGE:
          type = MED_SEG2;
          break;
        case MBPRISM:
          type = MED_PENTA6;
          break;
        default:
          type = MED_POINT1;
          break;
        }
        return type;
      };

      med_geometrie_element med_type = get_med_element_type(ent_type);
      if (ent_type == MBENTITYSET) {
        continue;
      }
      CHKERR MEDmeshElementWr(fid, mesh_name.c_str(), MED_NO_DT, MED_NO_IT, 0.,
                              MED_CELL, med_type, MED_NODAL, MED_FULL_INTERLACE,
                              family_number.size(), &connectivity[0], MED_FALSE,
                              nullptr, MED_TRUE, &tag_number[0], MED_TRUE,
                              &family_number[0]);

      MOFEM_LOG_C("MEDWORLD", Sev::inform,
                  "Writing %i elements of type %i (%s) ", family_number.size(),
                  med_type, moab::CN::EntityTypeName(ent_type));
    }

    // Close the MED file
    CHKERR MEDfermer(fid);
    MoFEMFunctionReturnHot(0);
}

MoFEMErrorCode MedInterface::readFields(const std::string &file_name,
                                        const std::string &field_name,
                                        const bool load_series,
                                        const int only_step, int verb) {

  Interface &m_field = cOre;
  MoFEMFunctionBeginHot;
  med_idt fid = MEDfileOpen((char *)file_name.c_str(), MED_LECTURE);
  if (fid < 0) {
    SETERRQ1(m_field.get_comm(), MOFEM_DATA_INCONSISTENCY,
             "Unable to open file '%s'", file_name.c_str());
  }

  med_int num_comp = MEDfieldnComponentByName(fid, field_name.c_str());
  if (num_comp <= 0) {
    SETERRQ(m_field.get_comm(), MOFEM_DATA_INCONSISTENCY,
            "Could not get number of components for MED field");
  }

  char meshName[MED_NAME_SIZE + 1];
  char dtUnit[MED_SNAME_SIZE + 1];
  std::vector<char> compName(num_comp * MED_SNAME_SIZE + 1);
  std::vector<char> compUnit(num_comp * MED_SNAME_SIZE + 1);
  med_int numSteps = 0;
  med_type_champ type;
  med_bool localMesh;
  if (MEDfieldInfoByName(fid, field_name.c_str(), meshName, &localMesh, &type,
                         &compName[0], &compUnit[0], dtUnit, &numSteps) < 0) {
    SETERRQ(m_field.get_comm(), MOFEM_DATA_INCONSISTENCY,
            "Could not get MED field info");
  }

  // Get meshset
  MeshsetsManager *meshsets_manager_ptr;
  CHKERR m_field.getInterface(meshsets_manager_ptr);
  const CubitMeshSets *cubit_meshset_ptr;
  CHKERR meshsets_manager_ptr->getCubitMeshsetPtr(meshName, &cubit_meshset_ptr);
  EntityHandle meshset = cubit_meshset_ptr->getMeshset();

  // Paraview can only plot field which have 1, 3, or 9 components. If field has
  // more that 9 comonents, it is stored on MOAB mesh but not viable in
  // ParaView.
  int num_comp_msh = (num_comp <= 1)   ? 1
                     : (num_comp <= 3) ? 3
                     : (num_comp <= 9) ? 9
                                       : num_comp;

  // Create tag to store nodal or cell values read form med file. Note that tag
  // has prefix MED to avoid collision with other tags.
  Tag th;
  std::string tag_name = "MED_" + field_name;
  {
    std::vector<double> def_val(num_comp_msh, 0);
    CHKERR m_field.get_moab().tag_get_handle(
        tag_name.c_str(), num_comp_msh, MB_TYPE_DOUBLE, th,
        MB_TAG_CREAT | MB_TAG_SPARSE, &def_val[0]);
  }

  // Warning! The ordering of the elements in the last two lists is
  // important: it should match the ordering of the MSH element types
  // (when elements are saved without tags, this governs the order
  // with which we implicitly index them in GModel::readMED)
  const med_entity_type entType[] = {MED_NODE, MED_CELL, MED_NODE_ELEMENT};
  const med_geometrie_element eleType[] = {
      MED_NONE,   MED_SEG2,   MED_TRIA3, MED_QUAD4,  MED_TETRA4,  MED_HEXA8,
      MED_PENTA6, MED_PYRA5,  MED_SEG3,  MED_TRIA6,  MED_QUAD9,   MED_TETRA10,
      MED_HEXA27, MED_POINT1, MED_QUAD8, MED_HEXA20, MED_PENTA15, MED_PYRA13};
  // const int nodesPerEle[] = {0, 2, 3, 4, 4, 8, 6, 5, 3, 6, 9, 10, 27, 1, 8,
  // 20, 15, 13};

  std::vector<std::pair<int, int>> pairs;
  for (unsigned int i = 0; i < sizeof(entType) / sizeof(entType[0]); i++) {
    for (unsigned int j = 0; j < sizeof(eleType) / sizeof(eleType[0]); j++) {
      if ((!i && !j) || j) {
        med_int n = numSteps;
        if (n > 0) {
          pairs.push_back(std::pair<int, int>(i, j));
          numSteps = std::max(numSteps, n);
        }
        if (!i && !j)
          break; // MED_NOEUD does not care about eleType
      }
    }
  }

  if (numSteps < 1 || pairs.empty()) {
    SETERRQ(m_field.get_comm(), MOFEM_DATA_INCONSISTENCY,
            "Nothing to import from MED file");
  }

  if (load_series) {
    SETERRQ(m_field.get_comm(), MOFEM_NOT_IMPLEMENTED, "Not yet implemented");
  }

  for (int step = (only_step == -1) ? 0 : only_step; step < numSteps; step++) {

    if (!load_series && only_step != step)
      continue;

    // cerr << only_step << " " << step << endl;

    // FIXME: in MED3 we might want to loop over all profiles instead
    // of relying of the default one

    // FIXME: MED3 allows to store multi-step meshes; we should

    for (unsigned int pair = 0; pair < pairs.size(); pair++) {

      // get step info
      med_entite_maillage ent = entType[pairs[pair].first];
      med_geometrie_element ele = eleType[pairs[pair].second];
      med_int numdt, numit, ngauss;
      med_float dt;
      if (MEDfieldComputingStepInfo(fid, field_name.c_str(), step + 1, &numdt,
                                    &numit, &dt) < 0) {
        SETERRQ(m_field.get_comm(), MOFEM_OPERATION_UNSUCCESSFUL,
                "Could not read step info");
      }

      char locName[MED_NAME_SIZE + 1], profileName[MED_NAME_SIZE + 1];

      // get number of values in the field (numVal takes the number of
      // Gauss points or the number of nodes per element into account,
      // but not the number of components)
      med_int profileSize;
      med_int numVal = MEDfieldnValueWithProfile(
          fid, field_name.c_str(), numdt, numit, ent, ele, 1,
          MED_COMPACT_STMODE, profileName, &profileSize, locName, &ngauss);
      numVal *= ngauss;

      if (numVal <= 0)
        continue;

      // int mult = 1;
      // if(ent == MED_NODE_ELEMENT) {
      //  mult = nodesPerEle[pairs[pair].second];
      //}
      // else if(ngauss != 1){
      //  mult = ngauss;
      //}

      // read field data
      std::vector<double> val(numVal * num_comp);
      if (MEDfieldValueWithProfileRd(fid, field_name.c_str(), numdt, numit, ent,
                                     ele, MED_COMPACT_STMODE, profileName,
                                     MED_FULL_INTERLACE, MED_ALL_CONSTITUENT,
                                     (unsigned char *)&val[0]) < 0) {
        SETERRQ(m_field.get_comm(), MOFEM_OPERATION_UNSUCCESSFUL,
                "Could not read field values");
      }

      // if (verb > 2) {
      //   // FIXME: This not looks ok for me
      //   cerr << ent << " " << ele << endl;
      //   cerr << string(meshName) << " : " << string(profileName) << " : "
      //        << string(locName) << " : " << profileSize << " : " << ngauss
      //        << endl;
      // }

      switch (ent) {
      case MED_CELL: {
        EntityType ent_type = MBMAXTYPE;
        switch (ele) {
        case MED_TETRA4:
        case MED_TETRA10:
          ent_type = MBTET;
          break;
        case MED_HEXA8:
          ent_type = MBHEX;
          break;
        default:
          MOFEM_LOG_C("MEDWORLD", Sev::warning,
                      "Not yet implemented for this cell %d", ele);
        }
        if (ent_type != MBMAXTYPE) {
          if (ngauss == 1) {
            Range ents;
            CHKERR m_field.get_moab().get_entities_by_type(meshset, ent_type,
                                                           ents, true);
            double e_vals[num_comp_msh];
            bzero(e_vals, sizeof(double) * num_comp_msh);
            std::vector<double>::iterator vit = val.begin();
            for (Range::iterator eit = ents.begin(); eit != ents.end(); eit++) {
              for (int ii = 0; ii != num_comp; ii++, vit++) {
                e_vals[ii] = *vit;
              }
              CHKERR m_field.get_moab().tag_set_data(th, &*eit, 1, e_vals);
            }
          } else {
            Range ents;
            CHKERR m_field.get_moab().get_entities_by_type(meshset, ent_type,
                                                           ents, true);
            if (ents.size() * ngauss * num_comp != val.size()) {
              SETERRQ(m_field.get_comm(), MOFEM_DATA_INCONSISTENCY,
                      "data inconsistency");
            }
            // FIXME simply average gauss values, far from perfect need fix
            double e_vals[num_comp_msh];
            std::vector<double>::iterator vit = val.begin();
            for (Range::iterator eit = ents.begin(); eit != ents.end(); eit++) {
              bzero(e_vals, sizeof(double) * num_comp_msh);
              for (int gg = 0; gg != ngauss; gg++) {
                for (int ii = 0; ii != num_comp; ii++, vit++) {
                  e_vals[ii] += *vit / ngauss;
                }
              }
              CHKERR m_field.get_moab().tag_set_data(th, &*eit, 1, e_vals);
            }
          }
        }
        // SETERRQ1(
        //   m_field.get_comm(),
        //   MOFEM_NOT_IMPLEMENTED,
        //   "Not implemented for more gauss pts ngauss = %d",
        //   ngauss
        // );
      } break;
      case MED_NODE:
      case MED_NODE_ELEMENT: {
        EntityType ent_type = MBVERTEX;
        Range ents;
        CHKERR m_field.get_moab().get_entities_by_type(meshset, ent_type, ents,
                                                       true);
        double e_vals[num_comp_msh];
        bzero(e_vals, sizeof(double) * num_comp_msh);
        std::vector<double>::iterator vit = val.begin();
        for (Range::iterator eit = ents.begin(); eit != ents.end(); eit++) {
          for (int ii = 0; ii != num_comp; ii++, vit++) {
            e_vals[ii] = *vit;
          }
          CHKERR m_field.get_moab().tag_set_data(th, &*eit, 1, e_vals);
        }
      } break;
      default:
        MOFEM_LOG_C("MEDWORLD", Sev::inform, "Entity type %d not implemented",
                    ent);
      }
    }
  }

  MoFEMFunctionReturnHot(0);
}

std::ostream &operator<<(std::ostream &os,
                         const MedInterface::FieldData &field_data) {
  os << "field name: " << field_data.fieldName;
  os << " mesh name: " << field_data.meshName;
  os << " local mesh: " << ((field_data.localMesh) ? "true" : "false");
  os << std::endl;
  // os << " field type: ";
  // switch (field_data.fieldType) {
  //   case MED_FLOAT64: os << "MED_FLOAT64"; break;
  //   case MED_INT32: os << "MED_INT32"; break;
  //   case MED_INT64: os << "MED_INT64"; break;
  //   case MED_INT: os << "MED_INT"; break;
  // };
  os << " componentNames:";
  for (unsigned int ff = 0; ff != field_data.componentNames.size(); ff++) {
    os << " " << field_data.componentNames[ff];
  }
  os << std::endl;
  os << " componentUnits:";
  for (unsigned int ff = 0; ff != field_data.componentUnits.size(); ff++) {
    os << " " << field_data.componentUnits[ff];
  }
  os << std::endl;
  os << " dtUnit: " << field_data.dtUnit << endl;
  os << " number of steps: " << field_data.ncSteps;
  return os;
}

  } // namespace MoFEM

#endif // WITH_MED
