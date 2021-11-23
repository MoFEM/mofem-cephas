/** \file BCMultiIndices.cpp
 * \brief Structures for managing boundary conditions
 */

/* MoFEM is free software: you can redistribute it and/or modify it under
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

namespace MoFEM {

// moab base meshsets
MoFEMErrorCode CubitMeshSets::getTagsHandlers(moab::Interface &moab) {
  MoFEMFunctionBegin;
  CHKERR moab.tag_get_handle(DIRICHLET_SET_TAG_NAME, nsTag);
  CHKERR moab.tag_get_handle(NEUMANN_SET_TAG_NAME, ssTag);
  CHKERR moab.tag_get_handle(
      (std::string(DIRICHLET_SET_TAG_NAME) + "__BC_DATA").c_str(), nsTag_data);
  CHKERR moab.tag_get_handle(
      (std::string(NEUMANN_SET_TAG_NAME) + "__BC_DATA").c_str(), ssTag_data);
  CHKERR moab.tag_get_handle(MATERIAL_SET_TAG_NAME, bhTag);
  CHKERR moab.tag_get_handle(BLOCK_HEADER, bhTag_header);
  CHKERR moab.tag_get_handle(BLOCK_ATTRIBUTES, thBlockAttribs);
  CHKERR moab.tag_get_handle(NAME_TAG_NAME, entityNameTag);
  MoFEMFunctionReturn(0);
}
CubitMeshSets::CubitMeshSets(moab::Interface &moab, const EntityHandle meshset)
    : meshset(meshset), cubitBcType(UNKNOWNSET), msId(nullptr),
      tagBcData(nullptr), tagBcSize(0), tagBlockHeaderData(nullptr),
      tagBlockAttributes(nullptr), tagBlockAttributesSize(0), tagName(nullptr),
      meshsetsMask(NODESET | SIDESET | BLOCKSET) {

  CHKERR getTagsHandlers(moab);
  CHKERR moab.tag_get_tags_on_entity(meshset, tag_handles);
  for (auto tit = tag_handles.begin(); tit != tag_handles.end(); tit++) {
    if (*tit == nsTag || *tit == ssTag || *tit == bhTag) {
      CHKERR moab.tag_get_by_ptr(*tit, &meshset, 1, (const void **)&msId);
    }
    if (*tit == nsTag) {
      if (*msId != -1) {
        cubitBcType = NODESET;
      }
    }
    if (*tit == ssTag) {
      if (*msId != -1) {
        cubitBcType = SIDESET;
      }
    }
    if (*tit == bhTag) {
      if (*msId != -1) {
        cubitBcType = BLOCKSET;
      }
    }
    if ((*tit == nsTag_data) || (*tit == ssTag_data)) {
      CHKERR moab.tag_get_by_ptr(*tit, &meshset, 1, (const void **)&tagBcData,
                                 &tagBcSize);

      CHKERR getTypeFromBcData(cubitBcType);
    }
    if (*tit == bhTag_header) {
      int tag_length;
      CHKERR moab.tag_get_length(*tit, tag_length);
      CHKERR moab.tag_get_by_ptr(*tit, &meshset, 1,
                                 (const void **)&tagBlockHeaderData);
      if (tagBlockHeaderData[1] > 0)
        cubitBcType |= MATERIALSET;
    }
    if (*tit == thBlockAttribs) {
      CHKERR moab.tag_get_by_ptr(*tit, &meshset, 1,
                                 (const void **)&tagBlockAttributes,
                                 &tagBlockAttributesSize);
    }
    if (*tit == entityNameTag) {
      CHKERR moab.tag_get_by_ptr(entityNameTag, &meshset, 1,
                                 (const void **)&tagName);
      CHKERR getTypeFromName(cubitBcType);
    }
  }

  // If BC set has name, unset UNKNOWNNAME
  if (cubitBcType.to_ulong() &
      (DISPLACEMENTSET | FORCESET | PRESSURESET | VELOCITYSET |
       ACCELERATIONSET | TEMPERATURESET | HEATFLUXSET | INTERFACESET)) {
    if ((cubitBcType & CubitBCType(UNKNOWNNAME)).any()) {
      cubitBcType = cubitBcType & (~CubitBCType(UNKNOWNNAME));
    }
  }
}
CubitMeshSets::CubitMeshSets(moab::Interface &moab,
                             const CubitBCType cubit_bc_type, const int ms_id)
    : cubitBcType(cubit_bc_type), msId(nullptr), tagBcData(nullptr),
      tagBcSize(0), tagBlockHeaderData(nullptr), tagBlockAttributes(nullptr),
      tagBlockAttributesSize(0), tagName(nullptr),
      meshsetsMask(NODESET | SIDESET | BLOCKSET) {

  CHKERR getTagsHandlers(moab);
  CHKERR moab.create_meshset(MESHSET_SET, meshset);
  switch (cubitBcType.to_ulong()) {
  case NODESET:
    CHKERR moab.tag_set_data(nsTag, &meshset, 1, &ms_id);
    CHKERR moab.tag_get_by_ptr(nsTag, &meshset, 1, (const void **)&msId);
    break;
  case SIDESET:
    CHKERR moab.tag_set_data(ssTag, &meshset, 1, &ms_id);
    CHKERR moab.tag_get_by_ptr(ssTag, &meshset, 1, (const void **)&msId);
    break;
  case BLOCKSET:
    CHKERR moab.tag_set_data(bhTag, &meshset, 1, &ms_id);
    CHKERR moab.tag_get_by_ptr(bhTag, &meshset, 1, (const void **)&msId);
    break;
  default: {
    PetscTraceBackErrorHandler(PETSC_COMM_WORLD, __LINE__, PETSC_FUNCTION_NAME,
                               __FILE__, MOFEM_DATA_INCONSISTENCY,
                               PETSC_ERROR_INITIAL, "Unknow meshset type",
                               PETSC_NULL);
    PetscMPIAbortErrorHandler(PETSC_COMM_WORLD, __LINE__, PETSC_FUNCTION_NAME,
                              __FILE__, MOFEM_DATA_INCONSISTENCY,
                              PETSC_ERROR_INITIAL, "Unknown meshset type",
                              PETSC_NULL);
  }
  }
}
MoFEMErrorCode CubitMeshSets::getMeshsetIdEntitiesByDimension(
    moab::Interface &moab, const int dimension, Range &entities,
    const bool recursive) const {
  MoFEMFunctionBegin;
  rval =
      moab.get_entities_by_dimension(meshset, dimension, entities, recursive);
  if (rval != MB_SUCCESS) {
    std::ostringstream ss;
    ss << "bc set " << *this << std::endl;
    PetscPrintf(PETSC_COMM_WORLD, ss.str().c_str());
  }
  CHKERR rval;
  MoFEMFunctionReturn(0);
}
MoFEMErrorCode CubitMeshSets::getMeshsetIdEntitiesByDimension(
    moab::Interface &moab, Range &entities, const bool recursive) const {
  MoFEMFunctionBegin;
  if ((cubitBcType & CubitBCType(BLOCKSET)).any()) {
    if (tagBlockHeaderData != nullptr) {
      return getMeshsetIdEntitiesByDimension(moab, tagBlockHeaderData[2],
                                             entities, recursive);
    } else {
      SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY, "dimension unknown");
    }
  }
  if ((cubitBcType & CubitBCType(NODESET)).any()) {
    return getMeshsetIdEntitiesByDimension(moab, 0, entities, recursive);
  }
  MoFEMFunctionReturn(0);
}
MoFEMErrorCode CubitMeshSets::getMeshsetIdEntitiesByType(
    moab::Interface &moab, const EntityType type, Range &entities,
    const bool recursive) const {
  MoFEMFunctionBegin;
  rval = moab.get_entities_by_type(meshset, type, entities, recursive);
  if (rval != MB_SUCCESS) {
    std::ostringstream ss;
    ss << "bc set " << *this << std::endl;
    PetscPrintf(PETSC_COMM_WORLD, ss.str().c_str());
  }
  CHKERR rval;
  MoFEMFunctionReturn(0);
}

MoFEMErrorCode CubitMeshSets::getBcData(std::vector<char> &bc_data) const {
  MoFEMFunctionBegin;
  bc_data.resize(tagBcSize);
  copy(&tagBcData[0], &tagBcData[tagBcSize], bc_data.begin());
  MoFEMFunctionReturn(0);
}

MoFEMErrorCode CubitMeshSets::getBlockHeaderData(
    std::vector<unsigned int> &material_data) const {
  MoFEMFunctionBegin;
  material_data.resize(3);
  copy(&tagBlockHeaderData[0], &tagBlockHeaderData[3], material_data.begin());
  MoFEMFunctionReturn(0);
}

MoFEMErrorCode CubitMeshSets::printBlockHeaderData(std::ostream &os) const {
  MoFEMFunctionBegin;
  if (tagBlockHeaderData != nullptr) {
    std::vector<unsigned int> material_data;
    getBlockHeaderData(material_data);
    os << "block_header_data = ";
    std::vector<unsigned int>::iterator vit = material_data.begin();
    for (; vit != material_data.end(); vit++) {
      os << std::hex << (int)((unsigned int)*vit) << " ";
    }
    os << ": ";
    vit = material_data.begin();
    for (; vit != material_data.end(); vit++) {
      os << *vit;
    }
    os << std::endl;
  } else {
    os << "no block header data" << std::endl;
  }
  MoFEMFunctionReturn(0);
}

std::string CubitMeshSets::getName() const {
  if (tagName != nullptr) {
    return std::string(tagName);
  } else {
    return "NoNameSet";
  }
}

MoFEMErrorCode CubitMeshSets::printName(std::ostream &os) const {
  MoFEMFunctionBegin;
  std::string name = getName();
  os << std::endl;
  os << "Block name:  " << name << std::endl;
  MoFEMFunctionReturn(0);
}

MoFEMErrorCode
CubitMeshSets::getTypeFromBcData(const std::vector<char> &bc_data,
                                 CubitBCType &type) const {
  MoFEMFunctionBegin;

  // See CubitBCType in common.hpp
  if (bc_data.size() == 0) {
    MoFEMFunctionReturnHot(0);
  }

  if (strcmp(&bc_data[0], "Displacement") == 0)
    type |= DISPLACEMENTSET;
  else if (strcmp(&bc_data[0], "Force") == 0)
    type |= FORCESET;
  else if (strcmp(&bc_data[0], "Velocity") == 0)
    type |= VELOCITYSET;
  else if (strcmp(&bc_data[0], "Acceleration") == 0)
    type |= ACCELERATIONSET;
  else if (strcmp(&bc_data[0], "Temperature") == 0)
    type |= TEMPERATURESET;
  else if (strcmp(&bc_data[0], "Pressure") == 0)
    type |= PRESSURESET;
  else if (strcmp(&bc_data[0], "HeatFlux") == 0)
    type |= HEATFLUXSET;
  else if (strcmp(&bc_data[0], "cfd_bc") == 0)
    type |= INTERFACESET;
  else
    type |= UNKNOWNNAME;

  MoFEMFunctionReturn(0);
}
MoFEMErrorCode CubitMeshSets::getTypeFromBcData(CubitBCType &type) const {
  MoFEMFunctionBegin;
  std::vector<char> bc_data;
  CHKERR getBcData(bc_data);
  CHKERR getTypeFromBcData(bc_data, type);
  MoFEMFunctionReturn(0);
}

MoFEMErrorCode CubitMeshSets::printBcData(std::ostream &os) const {
  MoFEMFunctionBegin;
  std::vector<char> bc_data;
  CHKERR getBcData(bc_data);
  os << "bc_data = ";
  std::vector<char>::iterator vit = bc_data.begin();
  for (; vit != bc_data.end(); vit++) {
    os << std::hex << (int)((unsigned char)*vit) << " ";
  }
  os << ": ";
  vit = bc_data.begin();
  for (; vit != bc_data.end(); vit++) {
    os << *vit;
  }
  os << std::endl;
  MoFEMFunctionReturn(0);
}

MoFEMErrorCode
CubitMeshSets::getAttributes(std::vector<double> &attributes) const {
  MoFEMFunctionBegin;
  attributes.resize(tagBlockAttributesSize);
  if (tagBlockAttributesSize > 0) {
    copy(&tagBlockAttributes[0], &tagBlockAttributes[tagBlockAttributesSize],
         attributes.begin());
  }
  MoFEMFunctionReturn(0);
}

MoFEMErrorCode
CubitMeshSets::setAttributes(moab::Interface &moab,
                             const std::vector<double> &attributes) {

  MoFEMFunctionBegin;
  int tag_size[] = {(int)attributes.size()};
  void const *tag_data[] = {&*attributes.begin()};
  CHKERR moab.tag_set_by_ptr(thBlockAttribs, &meshset, 1, tag_data, tag_size);
  CHKERR moab.tag_get_by_ptr(thBlockAttribs, &meshset, 1,
                             (const void **)&tagBlockAttributes,
                             &tagBlockAttributesSize);
  MoFEMFunctionReturn(0);
}

MoFEMErrorCode CubitMeshSets::printAttributes(std::ostream &os) const {
  MoFEMFunctionBegin;
  std::vector<double> attributes;
  CHKERR getAttributes(attributes);
  os << std::endl;
  os << "Block attributes" << std::endl;
  os << "----------------" << std::endl;
  for (unsigned int ii = 0; ii < attributes.size(); ii++) {
    os << "attr. no: " << ii + 1 << "   value: " << attributes[ii] << std::endl;
  }
  os << std::endl;
  MoFEMFunctionReturn(0);
}

MoFEMErrorCode CubitMeshSets::getTypeFromName(const std::string &name,
                                              CubitBCType &type) const {
  MoFEMFunctionBegin;
  // See CubitBCType in common.hpp
  if (name.compare(0, 11, "MAT_ELASTIC") == 0) {
    type |= MAT_ELASTICSET;
  } else if (name.compare(0, 11, "MAT_THERMAL") == 0) {
    type |= MAT_THERMALSET;
  } else if (name.compare(0, 12, "MAT_MOISTURE") == 0) {
    type |= MAT_MOISTURESET;
  } else if (name.compare(0, 10, "MAT_INTERF") == 0) {
    type |= MAT_INTERFSET;
  } else if (name.compare(0, 11, "BODY_FORCES") == 0) {
    type |= BODYFORCESSET;
  }
  // To be extended as appropriate
  else {
    type |= UNKNOWNNAME;
  }
  MoFEMFunctionReturn(0);
}

MoFEMErrorCode CubitMeshSets::getTypeFromName(CubitBCType &type) const {
  MoFEMFunctionBegin;
  std::string name = getName();
  CHKERR getTypeFromName(name, type);
  MoFEMFunctionReturn(0);
}

std::ostream &operator<<(std::ostream &os, const CubitMeshSets &e) {
  // get name of cubit meshset
  std::ostringstream ss;
  unsigned jj = 0;
  while (1 << jj != LASTSET_BC) {
    const CubitBCType jj_bc_type = 1 << jj;
    if ((e.cubitBcType & jj_bc_type).any()) {
      string bc_type_name;
      ss << " " << string(CubitBCNames[jj + 1]);
    }
    ++jj;
  }

  // push data to stream
  os << "meshset " << e.meshset << " type" << ss.str();
  if (e.msId != nullptr)
    os << " msId " << *(e.msId);
  if (e.tagName != nullptr) {
    os << " name " << e.getName();
  }
  if (e.tagBlockHeaderData != nullptr) {
    os << " block header: ";
    os << " blockCol = " << e.tagBlockHeaderData[0];
    os << " blockMat = " << e.tagBlockHeaderData[1];
    os << " blockDimension = " << e.tagBlockHeaderData[2];
  }
  return os;
}

void CubitMeshSets_change_name::operator()(CubitMeshSets &e) {

  switch (e.cubitBcType.to_ulong()) {
  case BLOCKSET: {
    nAme.resize(NAME_TAG_SIZE);
    CHKERR mOab.tag_set_data(e.entityNameTag, &e.meshset, 1, nAme.c_str());
    CHKERR mOab.tag_get_by_ptr(e.entityNameTag, &e.meshset, 1,
                               (const void **)&e.tagName);

    CubitBCType type;
    CHKERR e.getTypeFromName(type);
    e.cubitBcType |= type;
  }; break;
  case NODESET:
  case SIDESET: {
    nAme.resize(NAME_TAG_SIZE);
    CHKERR mOab.tag_set_data(e.entityNameTag, &e.meshset, 1, nAme.c_str());
    CHKERR mOab.tag_get_by_ptr(e.entityNameTag, &e.meshset, 1,
                               (const void **)&e.tagName);
  }; break;
  default:
    THROW_MESSAGE("not implemented for this CubitBC type");
  }
}

void CubitMeshSets_change_add_bit_to_cubit_bc_type::operator()(
    CubitMeshSets &e) {
  e.cubitBcType |= bIt;
}

void CubitMeshSets_change_attributes::operator()(CubitMeshSets &e) {
  CHKERR e.setAttributes(mOab, aTtr);
}

void CubitMeshSets_change_attributes_data_structure::operator()(
    CubitMeshSets &e) {
  // Need to run this to set tag size in number of doubles, don;t know nothing
  // about structure
  int tag_size[] = {(int)(aTtr.getSizeOfData() / sizeof(double))};
  void const *tag_data[] = {aTtr.getDataPtr()};
  CHKERR mOab.tag_set_by_ptr(e.thBlockAttribs, &e.meshset, 1, tag_data,
                             tag_size);
  CHKERR mOab.tag_get_by_ptr(e.thBlockAttribs, &e.meshset, 1,
                             (const void **)&e.tagBlockAttributes,
                             &e.tagBlockAttributesSize);
  // Here I know about structure
  CHKERR e.setAttributeDataStructure(aTtr);
}

void CubitMeshSets_change_bc_data_structure::operator()(CubitMeshSets &e) {

  // Need to run this to set tag size, don;t know nothing about structure
  int tag_size[] = {(int)bcData.getSizeOfData()};
  void const *tag_data[] = {bcData.getDataPtr()};
  if ((e.cubitBcType & CubitBCType(NODESET)).any()) {
    CHKERR mOab.tag_set_by_ptr(e.nsTag_data, &e.meshset, 1, tag_data, tag_size);
    CHKERR mOab.tag_get_by_ptr(e.nsTag_data, &e.meshset, 1,
                               (const void **)&e.tagBcData, &e.tagBcSize);
  } else if ((e.cubitBcType & CubitBCType(SIDESET)).any()) {
    CHKERR mOab.tag_set_by_ptr(e.ssTag_data, &e.meshset, 1, tag_data, tag_size);
    CHKERR mOab.tag_get_by_ptr(e.ssTag_data, &e.meshset, 1,
                               (const void **)&e.tagBcData, &e.tagBcSize);
  } else {
    THROW_MESSAGE("You have to have NODESET or SIDESET to apply BC data on it");
  }
  // Here I know about structure
  CHKERR e.setBcDataStructure(bcData);
  // Get Type form BC data
  CHKERR e.getTypeFromBcData(e.cubitBcType);
}

} // namespace MoFEM
