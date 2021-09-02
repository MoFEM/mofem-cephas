/** \file MeshsetsManager.cpp
 * \brief Interface to manage meshsets which carrying information about boundary
 * conditions and material blocks
 *
 */

/**
 * The MoFEM package is copyrighted by Lukasz Kaczmarczyk.
 * It can be freely used for educational and research purposes
 * by other institutions. If you use this software pleas cite my work.
 *
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
 * License along with MoFEM. If not, see <http://www.gnu.org/licenses/>
 */

using namespace std;
namespace po = boost::program_options;

namespace MoFEM {

bool MeshsetsManager::brodcastMeshsets = true;

MoFEMErrorCode
MeshsetsManager::query_interface(boost::typeindex::type_index type_index, UnknownInterface **iface) const {
  *iface = const_cast<MeshsetsManager *>(this);
  return 0;
}

MeshsetsManager::MeshsetsManager(const Core &core)
    : cOre(const_cast<Core &>(core)) {

  if (!LogManager::checkIfChannelExist("MeshsetMngWorld")) {
    auto core_log = logging::core::get();

    core_log->add_sink(
        LogManager::createSink(LogManager::getStrmWorld(), "MeshsetMngWorld"));
    core_log->add_sink(
        LogManager::createSink(LogManager::getStrmSync(), "MeshsetMngSync"));
    core_log->add_sink(
        LogManager::createSink(LogManager::getStrmSelf(), "MeshsetMngSelf"));

    LogManager::setLog("MeshsetMngWorld");
    LogManager::setLog("MeshsetMngSync");
    LogManager::setLog("MeshsetMngSelf");

    MOFEM_LOG_TAG("MeshsetMngWorld", "MeshsetMng");
    MOFEM_LOG_TAG("MeshsetMngSync", "MeshsetMng");
    MOFEM_LOG_TAG("MeshsetMngSelf", "MeshsetMng");
  }

  MOFEM_LOG("MeshsetMngWorld", Sev::noisy) << "Mashset manager created";
}

MoFEMErrorCode MeshsetsManager::clearMap() {
  MoFEMFunctionBeginHot;
  cubitMeshsets.clear();
  MoFEMFunctionReturnHot(0);
}

MoFEMErrorCode MeshsetsManager::initialiseDatabaseFromMesh(int verb) {
  Interface &m_field = cOre;
  MoFEMFunctionBegin;
  CHKERR readMeshsets(verb);
  if (brodcastMeshsets)
    CHKERR broadcastMeshsets(verb);

  for (auto &m : cubitMeshsets) {
    MOFEM_LOG("MeshsetMngWorld", Sev::inform) << m;
  }

  MoFEMFunctionReturn(0);
}

MoFEMErrorCode MeshsetsManager::readMeshsets(int verb) {
  Interface &m_field = cOre;
  moab::Interface &moab = m_field.get_moab();
  MoFEMFunctionBegin;

  Range meshsets;
  CHKERR moab.get_entities_by_type(0, MBENTITYSET, meshsets, false);
  for (auto m : meshsets) {
    // check if meshset is cubit meshset
    CubitMeshSets block(moab, m);
    if ((block.cubitBcType & CubitBCType(NODESET | SIDESET | BLOCKSET)).any()) {
      auto p = cubitMeshsets.insert(block);
      if (!p.second)
        SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
                "meshset not inserted");
      MOFEM_LOG("MeshsetMngSelf", Sev::noisy) << "read " << block;
    }
  }

  MoFEMFunctionReturn(0);
}

MoFEMErrorCode MeshsetsManager::broadcastMeshsets(int verb) {
  Interface &m_field = cOre;
  moab::Interface &moab = m_field.get_moab();
  MoFEMFunctionBegin;

 ParallelComm *pcomm = ParallelComm::get_pcomm(&moab, MYPCOMM_INDEX);
  if (pcomm == NULL)
    SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
            "MOAB communicator not set");

  const double coords[] = {0, 0, 0};

  auto set_tags_dummy_node = [&](const EntityHandle dummy_node,
                                 const EntityHandle meshset) {
    MoFEMFunctionBegin;
    std::vector<Tag> tag_handles;
    CHKERR moab.tag_get_tags_on_entity(meshset, tag_handles);
    for (auto th : tag_handles) {
      void *data[1];
      int tag_size;
      CHKERR moab.tag_get_by_ptr(th, &meshset, 1, (const void **)data,
                                 &tag_size);
      CHKERR moab.tag_set_by_ptr(th, &dummy_node, 1, data, &tag_size);
    }
    MoFEMFunctionReturn(0);
  };

  for (int from_proc = 0; from_proc < pcomm->size(); ++from_proc) {

    Range r_dummy_nodes;

    if (from_proc == pcomm->rank()) {
      std::vector<EntityHandle> dummy_nodes(cubitMeshsets.size(), 0);
      int i = 0;
      for (auto &m : cubitMeshsets) {
        CHKERR moab.create_vertex(coords, dummy_nodes[i]);
        CHKERR set_tags_dummy_node(dummy_nodes[i], m.getMeshset());
        ++i;
      }
      r_dummy_nodes.insert_list(dummy_nodes.begin(), dummy_nodes.end());
    }

    CHKERR pcomm->broadcast_entities(from_proc, r_dummy_nodes, false, true);

    for (auto dummy_node : r_dummy_nodes) {
      EntityHandle m;
      CHKERR moab.create_meshset(MESHSET_SET, m);
      CHKERR set_tags_dummy_node(m, dummy_node);

      CubitMeshSets broadcast_block(moab, m);
      if ((broadcast_block.cubitBcType &
           CubitBCType(NODESET | SIDESET | BLOCKSET))
              .any()) {
        auto p = cubitMeshsets.insert(broadcast_block);
        if (!p.second) {
          CHKERR moab.delete_entities(&m, 1);
        } else {
          MOFEM_LOG("MeshsetMngSelf", Sev::noisy)
              << "broadcast " << broadcast_block;
        }
      }
    }

    CHKERR moab.delete_entities(r_dummy_nodes);
  }

  MoFEMFunctionReturn(0);
}

MoFEMErrorCode MeshsetsManager::getTags(int verb) {
  MoFEMFunctionBegin;
  Interface &m_field = cOre;
  moab::Interface &moab = m_field.get_moab();
  int default_val = -1;
  CHKERR moab.tag_get_handle(DIRICHLET_SET_TAG_NAME, 1, MB_TYPE_INTEGER, nsTag,
                             MB_TAG_SPARSE | MB_TAG_CREAT, &default_val);

  CHKERR moab.tag_get_handle(NEUMANN_SET_TAG_NAME, 1, MB_TYPE_INTEGER, ssTag,
                             MB_TAG_SPARSE | MB_TAG_CREAT, &default_val);

  const int def_bc_data_len = 0;
  std::string tag_name = std::string(DIRICHLET_SET_TAG_NAME) + "__BC_DATA";
  CHKERR moab.tag_get_handle(
      tag_name.c_str(), def_bc_data_len, MB_TYPE_OPAQUE, nsTag_data,
      MB_TAG_CREAT | MB_TAG_SPARSE | MB_TAG_BYTES | MB_TAG_VARLEN, NULL);

  tag_name = std::string(NEUMANN_SET_TAG_NAME) + "__BC_DATA";
  CHKERR moab.tag_get_handle(
      tag_name.c_str(), def_bc_data_len, MB_TYPE_OPAQUE, ssTag_data,
      MB_TAG_CREAT | MB_TAG_SPARSE | MB_TAG_BYTES | MB_TAG_VARLEN, NULL);

  CHKERR moab.tag_get_handle(MATERIAL_SET_TAG_NAME, 1, MB_TYPE_INTEGER, bhTag,
                             MB_TAG_SPARSE | MB_TAG_CREAT, &default_val);

  std::vector<unsigned int> def_uint_zero(3, 0);
  CHKERR moab.tag_get_handle(
      BLOCK_HEADER, 3 * sizeof(unsigned int), MB_TYPE_INTEGER, bhTag_header,
      MB_TAG_CREAT | MB_TAG_SPARSE | MB_TAG_BYTES, &def_uint_zero[0]);

  Tag block_attribs;
  int def_Block_Attributes_length = 0;
  CHKERR moab.tag_get_handle(
      BLOCK_ATTRIBUTES, def_Block_Attributes_length, MB_TYPE_DOUBLE,
      block_attribs, MB_TAG_CREAT | MB_TAG_SPARSE | MB_TAG_VARLEN, NULL);

  Tag entity_name_tag;
  CHKERR moab.tag_get_handle(NAME_TAG_NAME, NAME_TAG_SIZE, MB_TYPE_OPAQUE,
                             entity_name_tag, MB_TAG_SPARSE | MB_TAG_CREAT);

  MoFEMFunctionReturn(0);
}

MoFEMErrorCode MeshsetsManager::printDisplacementSet() const {
  DisplacementCubitBcData mydata;
  MoFEMFunctionBegin;
  CHKERR printBcSet(mydata, NODESET | mydata.tYpe.to_ulong());
  MoFEMFunctionReturn(0);
}

MoFEMErrorCode MeshsetsManager::printPressureSet() const {
  PressureCubitBcData mydata;
  MoFEMFunctionBegin;
  CHKERR printBcSet(mydata, SIDESET | mydata.tYpe.to_ulong());
  MoFEMFunctionReturn(0);
}

MoFEMErrorCode MeshsetsManager::printForceSet() const {
  ForceCubitBcData mydata;
  MoFEMFunctionBegin;
  CHKERR printBcSet(mydata, NODESET | mydata.tYpe.to_ulong());
  MoFEMFunctionReturn(0);
}

MoFEMErrorCode MeshsetsManager::printTemperatureSet() const {
  TemperatureCubitBcData mydata;
  MoFEMFunctionBegin;
  CHKERR printBcSet(mydata, NODESET | mydata.tYpe.to_ulong());
  MoFEMFunctionReturn(0);
}

MoFEMErrorCode MeshsetsManager::printHeatFluxSet() const {
  HeatFluxCubitBcData mydata;
  MoFEMFunctionBegin;
  CHKERR printBcSet(mydata, SIDESET | mydata.tYpe.to_ulong());
  MoFEMFunctionReturn(0);
}

MoFEMErrorCode MeshsetsManager::printMaterialsSet() const {
  MoFEMFunctionBegin;
  const Interface &m_field = cOre;
  const moab::Interface &moab = m_field.get_moab();
  for (_IT_CUBITMESHSETS_BY_BCDATA_TYPE_FOR_LOOP_(
           (*this), BLOCKSET | MAT_ELASTICSET, it)) {
    Mat_Elastic data;
    CHKERR it->getAttributeDataStructure(data);
    MOFEM_LOG("MeshsetMngWorld", Sev::inform) << *it;
    MOFEM_LOG("MeshsetMngWorld", Sev::inform) << data;
    Range tets;
    CHKERR moab.get_entities_by_dimension(it->meshset, 3, tets, true);
    MOFEM_LOG("MeshsetMngWorld", Sev::inform)
        << "MAT_ELATIC msId " << it->getMeshsetId() << " nb. volumes "
        << tets.size();
  }

  for (_IT_CUBITMESHSETS_BY_BCDATA_TYPE_FOR_LOOP_(
           m_field, BLOCKSET | MAT_THERMALSET, it)) {
    Mat_Thermal data;
    CHKERR it->getAttributeDataStructure(data);
    MOFEM_LOG("MeshsetMngWorld", Sev::inform) << *it;
    MOFEM_LOG("MeshsetMngWorld", Sev::inform) << data;
  }

  for (_IT_CUBITMESHSETS_BY_BCDATA_TYPE_FOR_LOOP_(
           m_field, BLOCKSET | MAT_MOISTURESET, it)) {
    Mat_Moisture data;
    CHKERR it->getAttributeDataStructure(data);
    MOFEM_LOG("MeshsetMngWorld", Sev::inform) << *it;
    MOFEM_LOG("MeshsetMngWorld", Sev::inform) << data;
  }
  MoFEMFunctionReturn(0);
}

bool MeshsetsManager::checkMeshset(const int ms_id,
                                   const CubitBCType cubit_bc_type) const {
  auto miit =
      cubitMeshsets.get<Composite_Cubit_msId_And_MeshSetType_mi_tag>().find(
          boost::make_tuple(ms_id, cubit_bc_type.to_ulong()));
  if (miit !=
      cubitMeshsets.get<Composite_Cubit_msId_And_MeshSetType_mi_tag>().end()) {
    return true;
  }
  return false;
}

bool MeshsetsManager::checkMeshset(const string name,
                                   int *const number_of_meshsets_ptr) const {
  auto miit = cubitMeshsets.get<CubitMeshSets_name>().lower_bound(name);
  auto hi_miit = cubitMeshsets.get<CubitMeshSets_name>().upper_bound(name);
  if (std::distance(miit, hi_miit) == 0) {
    return false;
  }
  if (number_of_meshsets_ptr) {
    *number_of_meshsets_ptr = std::distance(miit, hi_miit);
  }
  return true;
}

MoFEMErrorCode MeshsetsManager::addMeshset(const CubitBCType cubit_bc_type,
                                           const int ms_id,
                                           const std::string name) {
  Interface &m_field = cOre;
  moab::Interface &moab = m_field.get_moab();
  MoFEMFunctionBegin;
  if (checkMeshset(ms_id, cubit_bc_type)) {
    SETERRQ1(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
             "such cubit meshset is already there", ms_id);
  }

  CubitMeshSets cmeshset(moab, cubit_bc_type, ms_id);
  if ((cmeshset.cubitBcType & CubitBCType(NODESET | SIDESET | BLOCKSET))
          .any()) {
    auto p = cubitMeshsets.insert(cmeshset);
    if (!p.second) {
      SETERRQ(PETSC_COMM_SELF, MOFEM_OPERATION_UNSUCCESSFUL,
              "meshset not inserted");
    }
    if (name.size() > 0) {
      bool success =
          cubitMeshsets.modify(p.first, CubitMeshSets_change_name(moab, name));
      if (!success) {
        SETERRQ(PETSC_COMM_SELF, MOFEM_OPERATION_UNSUCCESSFUL,
                "name to cubit meshset can not be set");
      }
    }
  }

  MoFEMFunctionReturn(0);
}

MoFEMErrorCode
MeshsetsManager::addEntitiesToMeshset(const CubitBCType cubit_bc_type,
                                      const int ms_id, const Range &ents) {
  Interface &m_field = cOre;
  moab::Interface &moab = m_field.get_moab();
  MoFEMFunctionBegin;
  auto cit =
      cubitMeshsets.get<Composite_Cubit_msId_And_MeshSetType_mi_tag>().find(
          boost::make_tuple(ms_id, cubit_bc_type.to_ulong()));
  if (cit ==
      cubitMeshsets.get<Composite_Cubit_msId_And_MeshSetType_mi_tag>().end()) {
    SETERRQ1(m_field.get_comm(), MOFEM_DATA_INCONSISTENCY,
             "Cannot find Cubit meshset with id: %d", ms_id);
  }
  EntityHandle meshset = cit->getMeshset();
  CHKERR moab.add_entities(meshset, ents);

  MoFEMFunctionReturn(0);
}

MoFEMErrorCode
MeshsetsManager::addEntitiesToMeshset(const CubitBCType cubit_bc_type,
                                      const int ms_id, const EntityHandle *ents,
                                      const int nb_ents) {
  Interface &m_field = cOre;
  moab::Interface &moab = m_field.get_moab();
  MoFEMFunctionBegin;
  auto cit =
      cubitMeshsets.get<Composite_Cubit_msId_And_MeshSetType_mi_tag>().find(
          boost::make_tuple(ms_id, cubit_bc_type.to_ulong()));
  if (cit ==
      cubitMeshsets.get<Composite_Cubit_msId_And_MeshSetType_mi_tag>().end()) {
    SETERRQ1(m_field.get_comm(), MOFEM_DATA_INCONSISTENCY,
             "Cannot find Cubit meshset with id: %d", ms_id);
  }
  EntityHandle meshset = cit->getMeshset();
  CHKERR moab.add_entities(meshset, ents, nb_ents);

  MoFEMFunctionReturn(0);
}

MoFEMErrorCode
MeshsetsManager::setAtributes(const CubitBCType cubit_bc_type, const int ms_id,
                              const std::vector<double> &attributes,
                              const std::string name) {
  Interface &m_field = cOre;
  moab::Interface &moab = m_field.get_moab();
  MoFEMFunctionBegin;
  auto cit =
      cubitMeshsets.get<Composite_Cubit_msId_And_MeshSetType_mi_tag>().find(
          boost::make_tuple(ms_id, cubit_bc_type.to_ulong()));
  if (cit ==
      cubitMeshsets.get<Composite_Cubit_msId_And_MeshSetType_mi_tag>().end()) {
    SETERRQ1(m_field.get_comm(), MOFEM_DATA_INCONSISTENCY,
             "Cannot find Cubit meshset with id: %d", ms_id);
  }
  if (name.size() > 0) {
    bool success = cubitMeshsets.modify(cubitMeshsets.project<0>(cit),
                                        CubitMeshSets_change_name(moab, name));
    if (!success) {
      SETERRQ(m_field.get_comm(), MOFEM_OPERATION_UNSUCCESSFUL,
              "name to cubit meshset can not be set");
    }
  }
  bool success =
      cubitMeshsets.modify(cubitMeshsets.project<0>(cit),
                           CubitMeshSets_change_attributes(moab, attributes));
  if (!success)
    SETERRQ(PETSC_COMM_SELF, MOFEM_OPERATION_UNSUCCESSFUL,
            "modification unsuccessful");

  std::ostringstream ss;
  ss << "Block " << cit->getName();
  ss << " Add attr: ";
  for (unsigned int ii = 0; ii != attributes.size(); ii++) {
    ss << attributes[ii] << " ";
  }
  MOFEM_LOG("MeshsetMngSelf", Sev::noisy) << ss.str();

  MoFEMFunctionReturn(0);
}

MoFEMErrorCode MeshsetsManager::setAtributesByDataStructure(
    const CubitBCType cubit_bc_type, const int ms_id,
    const GenericAttributeData &data, const std::string name) {
  Interface &m_field = cOre;
  moab::Interface &moab = m_field.get_moab();
  MoFEMFunctionBegin;
  auto cit =
      cubitMeshsets.get<Composite_Cubit_msId_And_MeshSetType_mi_tag>().find(
          boost::make_tuple(ms_id, cubit_bc_type.to_ulong()));
  if (cit ==
      cubitMeshsets.get<Composite_Cubit_msId_And_MeshSetType_mi_tag>().end()) {
    SETERRQ1(m_field.get_comm(), MOFEM_DATA_INCONSISTENCY,
             "Cannot find Cubit meshset with id: %d", ms_id);
  }
  if (name.size() > 0) {
    bool success = cubitMeshsets.modify(cubitMeshsets.project<0>(cit),
                                        CubitMeshSets_change_name(moab, name));
    if (!success) {
      SETERRQ(PETSC_COMM_SELF, MOFEM_OPERATION_UNSUCCESSFUL,
              "name to cubit meshset can not be set");
    }
  }
  bool success = cubitMeshsets.modify(
      cubitMeshsets.project<0>(cit),
      CubitMeshSets_change_attributes_data_structure(moab, data));
  if (!success)
    SETERRQ(PETSC_COMM_SELF, MOFEM_OPERATION_UNSUCCESSFUL,
            "modification unsuccessful");
  MoFEMFunctionReturn(0);
}

MoFEMErrorCode MeshsetsManager::setBcData(const CubitBCType cubit_bc_type,
                                          const int ms_id,
                                          const GenericCubitBcData &data) {
  Interface &m_field = cOre;
  moab::Interface &moab = m_field.get_moab();
  MoFEMFunctionBegin;
  auto cit =
      cubitMeshsets.get<Composite_Cubit_msId_And_MeshSetType_mi_tag>().find(
          boost::make_tuple(ms_id, cubit_bc_type.to_ulong()));
  if (cit ==
      cubitMeshsets.get<Composite_Cubit_msId_And_MeshSetType_mi_tag>().end()) {
    SETERRQ1(m_field.get_comm(), MOFEM_DATA_INCONSISTENCY,
             "Cubit meshset with id is already there", ms_id);
  }
  bool success =
      cubitMeshsets.modify(cubitMeshsets.project<0>(cit),
                           CubitMeshSets_change_bc_data_structure(moab, data));
  if (!success)
    SETERRQ(PETSC_COMM_SELF, MOFEM_OPERATION_UNSUCCESSFUL,
            "modification unsuccessful");
  MoFEMFunctionReturn(0);
}

MoFEMErrorCode MeshsetsManager::deleteMeshset(const CubitBCType cubit_bc_type,
                                              const int ms_id,
                                              const MoFEMTypes bh) {
  Interface &m_field = cOre;
  moab::Interface &moab = m_field.get_moab();
  MoFEMFunctionBegin;
  auto miit =
      cubitMeshsets.get<Composite_Cubit_msId_And_MeshSetType_mi_tag>().find(
          boost::make_tuple(ms_id, cubit_bc_type.to_ulong()));
  if (miit ==
      cubitMeshsets.get<Composite_Cubit_msId_And_MeshSetType_mi_tag>().end()) {
    if (bh & MF_EXIST) {
      SETERRQ1(m_field.get_comm(), MOFEM_DATA_INCONSISTENCY,
               "meshset not found", ms_id);
    } else {
      MoFEMFunctionReturnHot(0);
    }
  }
  EntityHandle meshset = miit->getMeshset();
  cubitMeshsets.get<Composite_Cubit_msId_And_MeshSetType_mi_tag>().erase(miit);
  CHKERR moab.delete_entities(&meshset, 1);

  MoFEMFunctionReturn(0);
}

MoFEMErrorCode MeshsetsManager::getCubitMeshsetPtr(
    const int ms_id, const CubitBCType cubit_bc_type,
    const CubitMeshSets **cubit_meshset_ptr) const {
  Interface &m_field = cOre;
  MoFEMFunctionBegin;
  auto miit =
      cubitMeshsets.get<Composite_Cubit_msId_And_MeshSetType_mi_tag>().find(
          boost::make_tuple(ms_id, cubit_bc_type.to_ulong()));
  if (miit !=
      cubitMeshsets.get<Composite_Cubit_msId_And_MeshSetType_mi_tag>().end()) {
    *cubit_meshset_ptr = &*miit;
  } else {
    SETERRQ1(m_field.get_comm(), MOFEM_DATA_INCONSISTENCY,
             "msId = %d is not there", ms_id);
  }
  MoFEMFunctionReturn(0);
}

MoFEMErrorCode MeshsetsManager::getCubitMeshsetPtr(
    const string name, const CubitMeshSets **cubit_meshset_ptr) const {
  Interface &m_field = cOre;
  MoFEMFunctionBegin;
  auto miit = cubitMeshsets.get<CubitMeshSets_name>().lower_bound(name);
  auto hi_miit = cubitMeshsets.get<CubitMeshSets_name>().upper_bound(name);
  if (std::distance(miit, hi_miit) == 0) {
    SETERRQ1(m_field.get_comm(), MOFEM_DATA_INCONSISTENCY,
             "meshset name <%s> is not there", name.c_str());
  }
  if (std::distance(miit, hi_miit) > 1) {
    SETERRQ1(m_field.get_comm(), MOFEM_DATA_INCONSISTENCY,
             "more that one meshser of that name <%s>", name.c_str());
  }
  *cubit_meshset_ptr = &*miit;
  MoFEMFunctionReturn(0);
}

MoFEMErrorCode MeshsetsManager::getEntitiesByDimension(
    const int msId, const unsigned int cubit_bc_type, const int dimension,
    Range &entities, const bool recursive) const {
  Interface &m_field = cOre;
  moab::Interface &moab = m_field.get_moab();
  MoFEMFunctionBegin;
  auto miit =
      cubitMeshsets.get<Composite_Cubit_msId_And_MeshSetType_mi_tag>().find(
          boost::make_tuple(msId, cubit_bc_type));
  if (miit !=
      cubitMeshsets.get<Composite_Cubit_msId_And_MeshSetType_mi_tag>().end()) {
    CHKERR miit->getMeshsetIdEntitiesByDimension(moab, dimension, entities,
                                                 recursive);
  } else {
    SETERRQ1(m_field.get_comm(), MOFEM_DATA_INCONSISTENCY,
             "msId = %d is not there", msId);
  }
  MoFEMFunctionReturn(0);
}

MoFEMErrorCode MeshsetsManager::getEntitiesByDimension(
    const int ms_id, const unsigned int cubit_bc_type, Range &entities,
    const bool recursive) const {
  Interface &m_field = cOre;
  moab::Interface &moab = m_field.get_moab();
  MoFEMFunctionBegin;
  auto miit =
      cubitMeshsets.get<Composite_Cubit_msId_And_MeshSetType_mi_tag>().find(
          boost::make_tuple(ms_id, cubit_bc_type));
  if (miit !=
      cubitMeshsets.get<Composite_Cubit_msId_And_MeshSetType_mi_tag>().end()) {
    CHKERR miit->getMeshsetIdEntitiesByDimension(moab, entities, recursive);
  } else {
    SETERRQ1(m_field.get_comm(), MOFEM_DATA_INCONSISTENCY,
             "ms_id = %d is not there", ms_id);
  }
  MoFEMFunctionReturn(0);
}

MoFEMErrorCode MeshsetsManager::getMeshset(const int ms_id,
                                           const unsigned int cubit_bc_type,
                                           EntityHandle &meshset) const {
  Interface &m_field = cOre;
  MoFEMFunctionBegin;
  auto miit =
      cubitMeshsets.get<Composite_Cubit_msId_And_MeshSetType_mi_tag>().find(
          boost::make_tuple(ms_id, cubit_bc_type));
  if (miit !=
      cubitMeshsets.get<Composite_Cubit_msId_And_MeshSetType_mi_tag>().end()) {
    meshset = miit->meshset;
  } else {
    SETERRQ1(m_field.get_comm(), MOFEM_DATA_INCONSISTENCY,
             "ms_id = %d is not there", ms_id);
  }
  MoFEMFunctionReturn(0);
}

bool MeshsetsManager::checkIfMeshsetContainsEntities(
    const int ms_id, const unsigned int cubit_bc_type,
    const EntityHandle *entities, int num_entities, const int operation_type) {
  auto miit =
      cubitMeshsets.get<Composite_Cubit_msId_And_MeshSetType_mi_tag>().find(
          boost::make_tuple(ms_id, cubit_bc_type));
  if (miit !=
      cubitMeshsets.get<Composite_Cubit_msId_And_MeshSetType_mi_tag>().end()) {
    Interface &m_field = cOre;
    return m_field.get_moab().contains_entities(miit->meshset, entities,
                                                num_entities, operation_type);
  } else
    return false;
}

MoFEMErrorCode
MeshsetsManager::getMeshsetsByType(const unsigned int cubit_bc_type,
                                   Range &meshsets) const {
  MoFEMFunctionBegin;
  auto miit =
      cubitMeshsets.get<CubitMeshSets_mi_tag>().lower_bound(cubit_bc_type);
  auto hi_miit =
      cubitMeshsets.get<CubitMeshSets_mi_tag>().upper_bound(cubit_bc_type);
  for (; miit != hi_miit; miit++) {
    meshsets.insert(miit->meshset);
  }
  MoFEMFunctionReturn(0);
}

struct BlockData {

  EntityHandle cubitMeshset;

  int iD;
  string addType;
  string nAme;
  CubitBC bcType;

  // Materials
  Mat_Elastic matElastic;
  Mat_Elastic_TransIso matTransIso;
  Mat_Thermal matThermal;
  Mat_Interf matInterf;

  // BCs
  DisplacementCubitBcData dispBc;
  ForceCubitBcData forceBc;
  PressureCubitBcData pressureBc;
  TemperatureCubitBcData temperatureBc;
  HeatFluxCubitBcData heatFluxBc;
  CfgCubitBcData cfgBc;

  int numberOfAttributes;
  std::vector<double> aTtr;

  BlockData() : aTtr(10, 0) {
    std::memcpy(dispBc.data.name, "Displacement", 12);
    std::memcpy(forceBc.data.name, "Force", 5);
    std::memcpy(pressureBc.data.name, "Pressure", 8);
    std::memcpy(temperatureBc.data.name, "Temperature", 11);
    std::memcpy(heatFluxBc.data.name, "HeatFlux", 8);
    std::memcpy(cfgBc.data.name, "cfd_bc", 6);
  }
};

MoFEMErrorCode MeshsetsManager::setMeshsetFromFile(const string file_name,
                                                   bool clean_file_options) {
  Interface &m_field = cOre;
  MoFEMFunctionBegin;
  std::ifstream ini_file(file_name.c_str(), std::ifstream::in);
  po::variables_map vm;
  if (clean_file_options) {
    configFileOptionsPtr =
        boost::shared_ptr<boost::program_options::options_description>(
            new po::options_description());
  }

  auto add_block_attributes = [&](auto prefix, auto &block_lists, auto &it) {
    // remove spacec from blockset name
    configFileOptionsPtr->add_options()(
        (prefix + ".number_of_attributes").c_str(),
        po::value<int>(&block_lists[it->getMeshsetId()].numberOfAttributes)
            ->default_value(-1),
        "Number of blockset attribute");
    for (int ii = 1; ii <= 10; ii++) {
      std::string surfix = ".user" + boost::lexical_cast<std::string>(ii);
      configFileOptionsPtr->add_options()(
          (prefix + surfix).c_str(),
          po::value<double>(&block_lists[it->getMeshsetId()].aTtr[ii - 1])
              ->default_value(0.0),
          "Add block attribute");
    }
  };

  // Add blocks
  map<int, BlockData> block_lists;
  for (_IT_CUBITMESHSETS_BY_SET_TYPE_FOR_LOOP_(m_field, BLOCKSET, it)) {
    block_lists[it->getMeshsetId()].cubitMeshset = it->getMeshset();
    std::string prefix =
        "block_" + boost::lexical_cast<std::string>(it->getMeshsetId());
    configFileOptionsPtr->add_options()(
        (prefix + ".add").c_str(),
        po::value<string>(&block_lists[it->getMeshsetId()].addType)
            ->default_value("UNKNOWNSET"),
        "Add block set")(
        (prefix + ".id").c_str(),
        po::value<int>(&block_lists[it->getMeshsetId()].iD)->default_value(-1),
        "Id of meshset")(
        (prefix + ".name").c_str(),
        po::value<string>(&block_lists[it->getMeshsetId()].nAme)
            ->default_value(""),
        "Name of the meshset");

    // Block attributes
    add_block_attributes(prefix, block_lists, it);

    // Mat elastic
    {
      // double Young; 			///< Young's modulus
      // double Poisson; 		///< Poisson's ratio
      // double ThermalExpansion;	///< Thermal expansion
      configFileOptionsPtr->add_options()(
          (prefix + ".young").c_str(),
          po::value<double>(
              &block_lists[it->getMeshsetId()].matElastic.data.Young)
              ->default_value(-1),
          "Young modulus")(
          (prefix + ".poisson").c_str(),
          po::value<double>(
              &block_lists[it->getMeshsetId()].matElastic.data.Poisson)
              ->default_value(-2),
          "Poisson ratio")(
          (prefix + ".thermalexpansion").c_str(),
          po::value<double>(
              &block_lists[it->getMeshsetId()].matElastic.data.ThermalExpansion)
              ->default_value(-1),
          "Thermal expansion");
      // TODO Add users parameters
    }
    // Mat Trans Iso
    {
      // Young's modulus in xy plane (Ep)
      // Young's modulus in z-direction (Ez)
      // Poisson's ratio in xy plane (vp)
      // Poisson's ratio in z-direction (vpz)
      // Shear modulus in z-direction (Gzp)
      configFileOptionsPtr->add_options()(
          (prefix + ".Youngp").c_str(),
          po::value<double>(
              &block_lists[it->getMeshsetId()].matTransIso.data.Youngp)
              ->default_value(-1),
          "Youngp")(
          (prefix + ".Youngz").c_str(),
          po::value<double>(
              &block_lists[it->getMeshsetId()].matTransIso.data.Youngz)
              ->default_value(-1),
          "Youngz")(
          (prefix + ".Poissonp").c_str(),
          po::value<double>(
              &block_lists[it->getMeshsetId()].matTransIso.data.Poissonp)
              ->default_value(0),
          "Poissonp")(
          (prefix + ".Poissonpz").c_str(),
          po::value<double>(
              &block_lists[it->getMeshsetId()].matTransIso.data.Poissonpz)
              ->default_value(0),
          "Poissonpz")(
          (prefix + ".Shearzp").c_str(),
          po::value<double>(
              &block_lists[it->getMeshsetId()].matTransIso.data.Shearzp)
              ->default_value(-1),
          "Shearzp");
      // TODO Add users parameters
    }
    // Mat thermal
    {
      // double Conductivity; ///< Thermal conductivity
      // double HeatCapacity; ///< Heat Capacity
      configFileOptionsPtr->add_options()(
          (prefix + ".conductivity").c_str(),
          po::value<double>(
              &block_lists[it->getMeshsetId()].matThermal.data.Conductivity)
              ->default_value(-1),
          "Conductivity")(
          (prefix + ".capacity").c_str(),
          po::value<double>(
              &block_lists[it->getMeshsetId()].matThermal.data.HeatCapacity)
              ->default_value(-1),
          "Capacity");
      // TODO Add users parameters
    }
    // Mat interface
    {
      // double alpha; ///< Elastic modulus multiplier
      // double beta;  ///< Damage Coupling multiplier between normal and
      // shear (g=sqrt(gn^2 + beta(gt1^2 + gt2^2))) double ft;    ///< Maximum
      // stress of crack double Gf;    ///< Fracture Energy
      configFileOptionsPtr->add_options()(
          (prefix + ".interface_alpha").c_str(),
          po::value<double>(
              &block_lists[it->getMeshsetId()].matInterf.data.alpha)
              ->default_value(-1),
          "alpha")((prefix + ".interface_beta").c_str(),
                   po::value<double>(
                       &block_lists[it->getMeshsetId()].matInterf.data.beta)
                       ->default_value(-1),
                   "beta")(
          (prefix + ".interface_ft").c_str(),
          po::value<double>(&block_lists[it->getMeshsetId()].matInterf.data.ft)
              ->default_value(-1),
          "ft")(
          (prefix + ".interface_Gf").c_str(),
          po::value<double>(&block_lists[it->getMeshsetId()].matInterf.data.Gf)
              ->default_value(-1),
          "Gf");
      // TODO Add users parameters
    }

    // Displacement bc
    {
      // char flag1; //< Flag for X-Translation (0: N/A, 1: specified)
      // char flag2; //< Flag for Y-Translation (0: N/A, 1: specified)
      // char flag3; //< Flag for Z-Translation (0: N/A, 1: specified)
      // char flag4; //< Flag for X-Rotation (0: N/A, 1: specified)
      // char flag5; //< Flag for Y-Rotation (0: N/A, 1: specified)
      // char flag6; //< Flag for Z-Rotation (0: N/A, 1: specified)
      // double value1; //< Value for X-Translation
      // double value2; //< Value for Y-Translation
      // double value3; //< Value for Z-Translation
      // double value4; //< Value for X-Rotation
      // double value5; //< Value for Y-Rotation
      // double value6; //< Value for Z-Rotation
      configFileOptionsPtr->add_options()(
          (prefix + ".disp_flag1").c_str(),
          po::value<char>(&block_lists[it->getMeshsetId()].dispBc.data.flag1)
              ->default_value(0),
          "flag1")(
          (prefix + ".disp_flag2").c_str(),
          po::value<char>(&block_lists[it->getMeshsetId()].dispBc.data.flag2)
              ->default_value(0),
          "flag2")(
          (prefix + ".disp_flag3").c_str(),
          po::value<char>(&block_lists[it->getMeshsetId()].dispBc.data.flag3)
              ->default_value(0),
          "flag3")(
          (prefix + ".disp_flag4").c_str(),
          po::value<char>(&block_lists[it->getMeshsetId()].dispBc.data.flag4)
              ->default_value(0),
          "flag4")(
          (prefix + ".disp_flag5").c_str(),
          po::value<char>(&block_lists[it->getMeshsetId()].dispBc.data.flag5)
              ->default_value(0),
          "flag5")(
          (prefix + ".disp_flag6").c_str(),
          po::value<char>(&block_lists[it->getMeshsetId()].dispBc.data.flag6)
              ->default_value(0),
          "flag6")(
          (prefix + ".disp_ux").c_str(),
          po::value<double>(&block_lists[it->getMeshsetId()].dispBc.data.value1)
              ->default_value(0),
          "value1")(
          (prefix + ".disp_uy").c_str(),
          po::value<double>(&block_lists[it->getMeshsetId()].dispBc.data.value2)
              ->default_value(0),
          "value2")(
          (prefix + ".disp_uz").c_str(),
          po::value<double>(&block_lists[it->getMeshsetId()].dispBc.data.value3)
              ->default_value(0),
          "value3")(
          (prefix + ".disp_rx").c_str(),
          po::value<double>(&block_lists[it->getMeshsetId()].dispBc.data.value4)
              ->default_value(0),
          "value4")(
          (prefix + ".disp_ry").c_str(),
          po::value<double>(&block_lists[it->getMeshsetId()].dispBc.data.value5)
              ->default_value(0),
          "value5")(
          (prefix + ".disp_rz").c_str(),
          po::value<double>(&block_lists[it->getMeshsetId()].dispBc.data.value6)
              ->default_value(0),
          "value6");
    }
    // Force BC data
    {
      // char zero[3]; //< 3 zeros
      // double value1; //< Force magnitude
      // double value2; //< Moment magnitude
      // double value3; //< X-component of force direction vector
      // double value4; //< Y-component of force direction vector
      // double value5; //< Z-component of force direction vector
      // double value6; //< X-component of moment direction vector
      // double value7; //< Y-component of moment direction vector
      // double value8; //< Z-component of moment direction vector
      // char zero2; // 0
      configFileOptionsPtr->add_options()(
          (prefix + ".force_magnitude").c_str(),
          po::value<double>(
              &block_lists[it->getMeshsetId()].forceBc.data.value1)
              ->default_value(0),
          "value1")((prefix + ".moment_magnitude").c_str(),
                    po::value<double>(
                        &block_lists[it->getMeshsetId()].forceBc.data.value2)
                        ->default_value(0),
                    "value2")(
          (prefix + ".force_fx").c_str(),
          po::value<double>(
              &block_lists[it->getMeshsetId()].forceBc.data.value3)
              ->default_value(0),
          "value3")((prefix + ".force_fy").c_str(),
                    po::value<double>(
                        &block_lists[it->getMeshsetId()].forceBc.data.value4)
                        ->default_value(0),
                    "value4")(
          (prefix + ".force_fz").c_str(),
          po::value<double>(
              &block_lists[it->getMeshsetId()].forceBc.data.value5)
              ->default_value(0),
          "value5")((prefix + ".moment_mx").c_str(),
                    po::value<double>(
                        &block_lists[it->getMeshsetId()].forceBc.data.value6)
                        ->default_value(0),
                    "value6")(
          (prefix + ".moment_my").c_str(),
          po::value<double>(
              &block_lists[it->getMeshsetId()].forceBc.data.value7)
              ->default_value(0),
          "value7")((prefix + ".moment_mz").c_str(),
                    po::value<double>(
                        &block_lists[it->getMeshsetId()].forceBc.data.value8)
                        ->default_value(0),
                    "value8");
    }
    {
      // char name[11]; //< 11 characters for "Temperature"
      // char pre1; //< This is always zero
      // char pre2; //< 0: temperature is not applied on thin shells
      // (default); 1: temperature is applied on thin shells char flag1; //<
      // 0: N/A, 1: temperature value applied (not on thin shells) char flag2;
      // //< 0: N/A, 1: temperature applied on thin shell middle char flag3;
      // //< 0: N/A, 1: thin shell temperature gradient specified char flag4;
      // //< 0: N/A, 1: top thin shell temperature char flag5; //< 0: N/A, 1:
      // bottom thin shell temperature char flag6; //< This is always zero
      // double value1; //< Temperature (default case - no thin shells)
      // double value2; //< Temperature for middle of thin shells
      // double value3; //< Temperature gradient for thin shells
      // double value4; //< Temperature for top of thin shells
      // double value5; //< Temperature for bottom of thin shells
      // double value6; //< This is always zero, i.e. ignore
      configFileOptionsPtr->add_options()(
          (prefix + ".temperature_flag1").c_str(),
          po::value<char>(
              &block_lists[it->getMeshsetId()].temperatureBc.data.flag1)
              ->default_value(0),
          "flag1")(
          (prefix + ".temperature_t").c_str(),
          po::value<double>(
              &block_lists[it->getMeshsetId()].temperatureBc.data.value1)
              ->default_value(0),
          "value1");
      // TODO: Add more cases, see above
    }
    // Sideset
    {
      // char name[8];   //< 8 characters for "Pressure"
      // char zero;      //< This is always zero
      // char flag2;     //< 0: Pressure is interpreted as pure pressure 1:
      // pressure is interpreted as total force double value1;  //< Pressure
      // value
      configFileOptionsPtr->add_options()(
          (prefix + ".pressure_flag2").c_str(),
          po::value<char>(
              &block_lists[it->getMeshsetId()].pressureBc.data.flag2)
              ->default_value(0),
          "flag2")((prefix + ".pressure_magnitude").c_str(),
                   po::value<double>(
                       &block_lists[it->getMeshsetId()].pressureBc.data.value1)
                       ->default_value(0),
                   "value1");
    }
    {
      // char name[8]; //< 8 characters for "HeatFlux" (no space)
      // char pre1; //< This is always zero
      // char pre2; //< 0: heat flux is not applied on thin shells (default);
      // 1: heat flux is applied on thin shells char flag1; //< 0: N/A, 1:
      // normal heat flux case (i.e. single value, case without thin shells)
      // char flag2; //< 0: N/A, 1: Thin shell top heat flux specified
      // char flag3; //< 0: N/A, 1: Thin shell bottom heat flux specidied
      // double value1; //< Heat flux value for default case (no thin shells)
      // double value2; //< Heat flux (thin shell top)
      // double value3; //< Heat flux (thin shell bottom)
      configFileOptionsPtr->add_options()(
          (prefix + ".heatflux_flag1").c_str(),
          po::value<char>(
              &block_lists[it->getMeshsetId()].heatFluxBc.data.flag1)
              ->default_value(0),
          "flag1")((prefix + ".heatflux_magnitude").c_str(),
                   po::value<double>(
                       &block_lists[it->getMeshsetId()].heatFluxBc.data.value1)
                       ->default_value(0),
                   "value1");
    }
    // Interface set
    {
      configFileOptionsPtr->add_options()(
          (prefix + ".interface_type").c_str(),
          po::value<char>(&block_lists[it->getMeshsetId()].cfgBc.data.type)
              ->default_value(0),
          "type");
    }
  }

  map<int, BlockData> block_set_attributes;
  for (_IT_CUBITMESHSETS_BY_SET_TYPE_FOR_LOOP_(m_field, BLOCKSET, it)) {
    block_set_attributes[it->getMeshsetId()].cubitMeshset = it->getMeshset();
    block_set_attributes[it->getMeshsetId()].iD = it->getMeshsetId();
    block_set_attributes[it->getMeshsetId()].bcType = BLOCKSET;
    std::string block_name = it->getName();
    block_name.erase(
        std::remove_if(block_name.begin(), block_name.end(), ::isspace),
        block_name.end());
    block_set_attributes[it->getMeshsetId()].nAme = block_name;
    // Only blocks which have name
    if (block_name.compare("NoNameSet") != 0) {
      std::string prefix = "SET_ATTR_" + block_name;
      // Block attributes
      add_block_attributes(prefix, block_set_attributes, it);
    }
  }

  po::parsed_options parsed =
      parse_config_file(ini_file, *configFileOptionsPtr, true);
  store(parsed, vm);
  po::notify(vm);

  // Set type from name
  for (_IT_CUBITMESHSETS_BY_SET_TYPE_FOR_LOOP_(m_field, BLOCKSET, it)) {

    CubitBCType bc_type;
    unsigned jj = 0;
    while (1 << jj != LASTSET_BC) {
      if (string(CubitBCNames[jj + 1]) ==
          block_lists[it->getMeshsetId()].addType) {
        bc_type = 1 << jj;
      }
      ++jj;
    }
    if (bc_type.none()) {
      block_lists[it->getMeshsetId()].bcType = UNKNOWNSET;
      // Skip the bockset nothing is defined for it
      continue;
    }

    if (bc_type.to_ulong() == BLOCKSET)
      block_lists[it->getMeshsetId()].bcType = BLOCKSET;
    else if (bc_type.to_ulong() == NODESET)
      block_lists[it->getMeshsetId()].bcType = NODESET;
    else if (bc_type.to_ulong() == SIDESET)
      block_lists[it->getMeshsetId()].bcType = SIDESET;
    else {
      SETERRQ1(m_field.get_comm(), MOFEM_DATA_INCONSISTENCY,
               "Not yet implemented type %s\n",
               block_lists[it->getMeshsetId()].addType.c_str());
    }
    if (block_lists[it->getMeshsetId()].iD == -1) {
      SETERRQ1(m_field.get_comm(), MOFEM_DATA_INCONSISTENCY,
               "Unset iD number %d\n", block_lists[it->getMeshsetId()].iD);
    }
  }

  std::vector<std::string> additional_parameters;
  additional_parameters =
      collect_unrecognized(parsed.options, po::include_positional);
  for (std::vector<std::string>::iterator vit = additional_parameters.begin();
       vit != additional_parameters.end(); vit++) {
    MOFEM_LOG_C("MeshsetMngSelf", Sev::warning, "Unrecognized option %s",
                vit->c_str());
  }
  for (map<int, BlockData>::iterator mit = block_lists.begin();
       mit != block_lists.end(); mit++) {
    CubitMeshSet_multiIndex::iterator cubit_meshset_it =
        cubitMeshsets.find(mit->second.cubitMeshset);
    if (cubit_meshset_it == cubitMeshsets.end()) {
      SETERRQ(m_field.get_comm(), MOFEM_DATA_INCONSISTENCY,
              "Data inconsistency\n");
    }
    switch (mit->second.bcType) {
    case UNKNOWNSET:
      break;
    case BLOCKSET: {
      if ((CubitBCType(mit->second.bcType) & cubit_meshset_it->getBcType())
              .any() &&
          mit->second.iD == cubit_meshset_it->getMeshsetId()) {
        // Meshset is the same, only modification
      } else {
        CHKERR addMeshset(mit->second.bcType, mit->second.iD, mit->second.nAme);
        EntityHandle meshset = cubit_meshset_it->getMeshset();
        CHKERR addEntitiesToMeshset(mit->second.bcType, mit->second.iD,
                                    &meshset, 1);
      }
      // Add attributes
      CHKERR setAtributes(mit->second.bcType, mit->second.iD, mit->second.aTtr);
      // Add material elastic data if value are physical (i.e. Young > 0,
      // Poisson in (-1.0.5) and ThermalExpansion>0)
      if (mit->second.matElastic.data.Young != -1) {
        CHKERR setAtributesByDataStructure(mit->second.bcType, mit->second.iD,
                                           mit->second.matElastic);
      }
      if (mit->second.matTransIso.data.Youngp != -1) {
        CHKERR setAtributesByDataStructure(mit->second.bcType, mit->second.iD,
                                           mit->second.matTransIso);
      }
      if (mit->second.matThermal.data.Conductivity != -1) {
        CHKERR setAtributesByDataStructure(mit->second.bcType, mit->second.iD,
                                           mit->second.matThermal);
      }
      if (mit->second.matInterf.data.ft != -1) {
        CHKERR setAtributesByDataStructure(mit->second.bcType, mit->second.iD,
                                           mit->second.matInterf);
      }
    } break;
    case NODESET: {
      if ((CubitBCType(mit->second.bcType) & cubit_meshset_it->getBcType())
              .any() &&
          mit->second.iD == cubit_meshset_it->getMeshsetId()) {
        // Meshset is the same, only modification
      } else {
        CHKERR addMeshset(mit->second.bcType, mit->second.iD);
        EntityHandle meshset = cubit_meshset_it->getMeshset();
        CHKERR addEntitiesToMeshset(mit->second.bcType, mit->second.iD,
                                    &meshset, 1);
      }
      // Add displacement bc
      if (mit->second.dispBc.data.flag1 || mit->second.dispBc.data.flag2 ||
          mit->second.dispBc.data.flag3 || mit->second.dispBc.data.flag4 ||
          mit->second.dispBc.data.flag5 || mit->second.dispBc.data.flag6) {
        if (mit->second.dispBc.data.flag1 == '0')
          mit->second.dispBc.data.flag1 = 0;
        if (mit->second.dispBc.data.flag1 == 'N')
          mit->second.dispBc.data.flag1 = 0;
        if (mit->second.dispBc.data.flag1)
          mit->second.dispBc.data.flag1 = 1;
        if (mit->second.dispBc.data.flag2 == '0')
          mit->second.dispBc.data.flag2 = 0;
        if (mit->second.dispBc.data.flag2 == 'N')
          mit->second.dispBc.data.flag2 = 0;
        if (mit->second.dispBc.data.flag2)
          mit->second.dispBc.data.flag2 = 1;
        if (mit->second.dispBc.data.flag3 == '0')
          mit->second.dispBc.data.flag3 = 0;
        if (mit->second.dispBc.data.flag3 == 'N')
          mit->second.dispBc.data.flag3 = 0;
        if (mit->second.dispBc.data.flag3)
          mit->second.dispBc.data.flag3 = 1;
        if (mit->second.dispBc.data.flag4 == '0')
          mit->second.dispBc.data.flag4 = 0;
        if (mit->second.dispBc.data.flag4 == 'N')
          mit->second.dispBc.data.flag4 = 0;
        if (mit->second.dispBc.data.flag4)
          mit->second.dispBc.data.flag4 = 1;
        if (mit->second.dispBc.data.flag5 == '0')
          mit->second.dispBc.data.flag5 = 0;
        if (mit->second.dispBc.data.flag5 == 'N')
          mit->second.dispBc.data.flag5 = 0;
        if (mit->second.dispBc.data.flag5)
          mit->second.dispBc.data.flag5 = 1;
        if (mit->second.dispBc.data.flag6 == '0')
          mit->second.dispBc.data.flag6 = 0;
        if (mit->second.dispBc.data.flag6 == 'N')
          mit->second.dispBc.data.flag6 = 0;
        if (mit->second.dispBc.data.flag6)
          mit->second.dispBc.data.flag6 = 1;
        CHKERR setBcData(mit->second.bcType, mit->second.iD,
                         mit->second.dispBc);
      }
      if (mit->second.forceBc.data.value1 != 0 ||
          mit->second.forceBc.data.value2 != 0) {
        CHKERR setBcData(mit->second.bcType, mit->second.iD,
                         mit->second.forceBc);
      }
      // Add temperature boundary condition
      if (mit->second.temperatureBc.data.flag1) {
        if (mit->second.temperatureBc.data.flag1 == '0')
          mit->second.temperatureBc.data.flag1 = 0;
        if (mit->second.temperatureBc.data.flag1 == 'N')
          mit->second.temperatureBc.data.flag1 = 0;
        if (mit->second.temperatureBc.data.flag1)
          mit->second.temperatureBc.data.flag1 = 1;
        CHKERR setBcData(mit->second.bcType, mit->second.iD,
                         mit->second.temperatureBc);
      }
    } break;
    case SIDESET: {
      if ((CubitBCType(mit->second.bcType) & cubit_meshset_it->getBcType())
              .any() &&
          mit->second.iD == cubit_meshset_it->getMeshsetId()) {
        // Meshset is the same, only modification
      } else {
        CHKERR addMeshset(mit->second.bcType, mit->second.iD);
        EntityHandle meshset = cubit_meshset_it->getMeshset();
        CHKERR addEntitiesToMeshset(mit->second.bcType, mit->second.iD,
                                    &meshset, 1);
      }
      // Add pressure
      if (mit->second.pressureBc.data.value1 != 0) {
        if (mit->second.pressureBc.data.flag2 == '0')
          mit->second.pressureBc.data.flag2 = 0;
        if (mit->second.pressureBc.data.flag2 == 'N')
          mit->second.pressureBc.data.flag2 = 0;
        if (mit->second.pressureBc.data.flag2)
          mit->second.pressureBc.data.flag2 = 1;
        CHKERR setBcData(mit->second.bcType, mit->second.iD,
                         mit->second.pressureBc);
      }
      // Add heat flux
      if (mit->second.heatFluxBc.data.value1 != 0) {
        if (mit->second.heatFluxBc.data.flag1 == '0')
          mit->second.heatFluxBc.data.flag1 = 0;
        if (mit->second.heatFluxBc.data.flag1 == 'N')
          mit->second.heatFluxBc.data.flag1 = 0;
        if (mit->second.heatFluxBc.data.flag1)
          mit->second.heatFluxBc.data.flag1 = 1;
        CHKERR setBcData(mit->second.bcType, mit->second.iD,
                         mit->second.heatFluxBc);
      }
      // Add Interface
      if (mit->second.cfgBc.data.type != 0) {
        CHKERR setBcData(mit->second.bcType, mit->second.iD, mit->second.cfgBc);
      }
    } break;
    default:
      SETERRQ(m_field.get_comm(), MOFEM_DATA_INCONSISTENCY,
              "Not yet implemented type\n");
    }
  }

  for (auto set_attr : block_set_attributes) {
    // Add attributes
    if (set_attr.second.numberOfAttributes > 0) {
      MOFEM_LOG("MeshsetMngSelf", Sev::verbose)
          << "Set attributes to blockset " << set_attr.second.nAme;
      set_attr.second.aTtr.resize(set_attr.second.numberOfAttributes);
      CHKERR setAtributes(set_attr.second.bcType, set_attr.second.iD,
                          set_attr.second.aTtr);
    }
  }

  MoFEMFunctionReturn(0);
}

MoFEMErrorCode MeshsetsManager::setMeshsetFromFile() {
  Interface &m_field = cOre;
  PetscBool flg_file;
  char meshset_file_name[255] = "config_file.cfg";
  MoFEMFunctionBegin;
  CHKERR PetscOptionsBegin(m_field.get_comm(), "", "Set meshsets form file",
                           "none");
  CHKERR PetscOptionsString("-meshsets_config", "meshsets config  file name",
                            "", "add_cubit_meshsets.in", meshset_file_name, 255,
                            &flg_file);
  if (flg_file == PETSC_TRUE) {
    ifstream f(meshset_file_name);
    if (!f.good()) {
      SETERRQ1(m_field.get_comm(), MOFEM_DATA_INCONSISTENCY,
               "File configuring meshsets ( %s ) can not be open\n",
               meshset_file_name);
    }
    CHKERR setMeshsetFromFile(string(meshset_file_name));
  }
  ierr = PetscOptionsEnd();
  CHKERRG(ierr);
  MoFEMFunctionReturn(0);
}

MoFEMErrorCode MeshsetsManager::saveMeshsetToFile(
    const int ms_id, const unsigned int cubit_bc_type,
    const std::string file_name, const std::string file_type,
    const std::string options) const {

  MoFEMFunctionBegin;
  MoFEM::Interface &m_field = cOre;
  const CubitMeshSets *cubit_meshset_ptr;
  CHKERR getCubitMeshsetPtr(ms_id, cubit_bc_type, &cubit_meshset_ptr);
  EntityHandle meshset = cubit_meshset_ptr->getMeshset();
  CHKERR m_field.get_moab().write_file(file_name.c_str(), file_type.c_str(),
                                       options.c_str(), &meshset, 1);
  MoFEMFunctionReturn(0);
}

MoFEMErrorCode MeshsetsManager::saveMeshsetToFile(
    const int ms_id, const unsigned int cubit_bc_type, const int dim,
    const std::string file_name, const bool recursive,
    const std::string file_type, const std::string options) const {

  MoFEMFunctionBegin;
  MoFEM::Interface &m_field = cOre;
  moab::Interface &moab = m_field.get_moab();
  Range entities;
  CHKERR getEntitiesByDimension(ms_id, cubit_bc_type, dim, entities, recursive);
  EntityHandle meshset;
  CHKERR moab.create_meshset(MESHSET_SET, meshset);
  CHKERR moab.add_entities(meshset, entities);
  CHKERR moab.write_file(file_name.c_str(), file_type.c_str(), options.c_str(),
                         &meshset, 1);
  CHKERR moab.delete_entities(&meshset, 1);
  MoFEMFunctionReturn(0);
}

MoFEMErrorCode
MeshsetsManager::updateAllMeshsetsByEntitiesChildren(const BitRefLevel &bit) {
  MoFEMFunctionBegin;
  BitRefManager *bit_mng = cOre.getInterface<BitRefManager>();
  for (_IT_CUBITMESHSETS_FOR_LOOP_((*this), iit)) {
    EntityHandle meshset = iit->getMeshset();
    for (EntityType t = MBVERTEX; t != MBENTITYSET; ++t)
      CHKERR bit_mng->updateMeshsetByEntitiesChildren(meshset, bit, meshset, t,
                                                      true);
  }
  MoFEMFunctionReturn(0);
}

} // namespace MoFEM
