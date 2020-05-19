/** \file CoordSystemsManager.cpp
 * \brief Interface managing coordinate systems set to fields
 *
 */

/*
 * MoFEM is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
 * License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with MoFEM. If not, see <http://www.gnu.org/licenses/>
 */

/** \file MeshsetsManager.cpp
 * \brief Interface to manage meshsets which carrying information about boundary
 * conditions and material blocks
 *
 */

/**
 * The MoFEM package is copyrighted by Lukasz Kaczmarczyk.
 * It can be freely used for educational and research purposes
 * by other institutions. If you use this softwre pleas cite my work.
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

namespace MoFEM {

MoFEMErrorCode
CoordSystemsManager::query_interface(const MOFEMuuid &uuid,
                                     UnknownInterface **iface) const {

  MoFEMFunctionBeginHot;
  *iface = NULL;
  if (uuid == IDD_MOFEMMeshsetsManager) {
    *iface = const_cast<CoordSystemsManager *>(this);
    MoFEMFunctionReturnHot(0);
  }
  SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY, "unknown interface");
  MoFEMFunctionReturnHot(0);
}

CoordSystemsManager::CoordSystemsManager(const Core &core)
    : cOre(const_cast<Core &>(core)) {}

CoordSystemsManager::~CoordSystemsManager() {}

MoFEMErrorCode CoordSystemsManager::getTags(int verb) {

  Interface &m_field = cOre;
  moab::Interface &moab = m_field.get_moab();
  MoFEMFunctionBegin;
  // Coordinate systems
  const int def_coord_sys_dim[] = {0, 0, 0, 0};
  CHKERR moab.tag_get_handle("_CoordSysDim", 4, MB_TYPE_INTEGER, th_CoordSysDim,
                             MB_TAG_CREAT | MB_TAG_SPARSE, &def_coord_sys_dim);
  const int def_val_len = 0;
  CHKERR moab.tag_get_handle(
      "_CoordSysName", def_val_len, MB_TYPE_OPAQUE, th_CoordSysName,
      MB_TAG_CREAT | MB_TAG_BYTES | MB_TAG_VARLEN | MB_TAG_SPARSE, NULL);
  MoFEMFunctionReturn(0);
}

MoFEMErrorCode CoordSystemsManager::clearMap() {
  MoFEMFunctionBeginHot;
  coordinateSystems.clear();
  MoFEMFunctionReturnHot(0);
}

MoFEMErrorCode CoordSystemsManager::initialiseDatabaseFromMesh(int verb) {
  Interface &m_field = cOre;
  moab::Interface &moab = m_field.get_moab();
  MoFEMFunctionBegin;

  Range meshsets;
  CHKERR moab.get_entities_by_type(0, MBENTITYSET, meshsets, false);

  // Iterate all mesh sets to find coordinate system
  for (auto meshset : meshsets) {
    const void *cs_name_ptr = nullptr;
    int cs_name_size = 0;
    rval = moab.tag_get_by_ptr(th_CoordSysName, &meshset, 1,
                               (const void **)&cs_name_ptr, &cs_name_size);
    if (rval == MB_SUCCESS && cs_name_size) {
      std::string cs_name(static_cast<const char *>(cs_name_ptr), cs_name_size);

      std::array<int, 4> dim;
      rval = moab.tag_get_data(th_CoordSysDim, &meshset, 1, dim.data());
      if (rval == MB_SUCCESS && (dim[0] + dim[1] + dim[2] + dim[3]) != 0) {

        std::pair<CoordSys_multiIndex::iterator, bool> p =
            coordinateSystems.insert(
                boost::make_shared<CoordSys>(moab, meshset));

        if (!p.second) {
          // Coordinate system is in database. Could be created by another
          // processor Check consistency of both coordinate systems with the
          // same name.
          if (((*p.first)->getDim(0) != dim[0]) ||
              ((*p.first)->getDim(1) != dim[1]) ||
              ((*p.first)->getDim(2) != dim[2]) ||
              ((*p.first)->getDim(3) != dim[3])) {
            SETERRQ2(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
                     "meshset to coord system not inserted "
                     "cs_name %s dim = %d",
                     cs_name.c_str(), dim[0] + dim[1] + dim[2] + dim[3]);
          } else {
            // Remove duplicate mesheset
            CHKERR moab.delete_entities(&meshset, 1);
          }
        }
      }
    }
  }

  MoFEMFunctionReturn(0);
}

MoFEMErrorCode CoordSystemsManager::addCoordinateSystem(
    const int cs_dim[], const std::string cs_name, const enum MoFEMTypes bh) {

  auto check_cs = [&](const std::string cs_name) {
    auto undefined_cs_it =
        coordinateSystems.get<CoordSysName_mi_tag>().find(cs_name);
    if (undefined_cs_it == coordinateSystems.get<CoordSysName_mi_tag>().end())
      return false;
    else
      return true;
  };

  MoFEMFunctionBegin;

  if (!check_cs(cs_name)) {

    Interface &m_field = cOre;
    moab::Interface &moab = m_field.get_moab();
    EntityHandle meshset;
    CHKERR moab.create_meshset(MESHSET_SET, meshset);
    CHKERR moab.tag_set_data(th_CoordSysDim, &meshset, 1, cs_dim);
    void const *cs_name_ptr[] = {cs_name.c_str()};
    int cs_name_size[1];
    cs_name_size[0] = cs_name.size();
    CHKERR moab.tag_set_by_ptr(th_CoordSysName, &meshset, 1, cs_name_ptr,
                               cs_name_size);
    auto p =
        coordinateSystems.insert(boost::make_shared<CoordSys>(moab, meshset));
    if (!p.second)
      SETERRQ1(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
               "MeshSet to coord system <%s> not inserted", cs_name.c_str());

  } else if (bh == MF_EXCL)

    SETERRQ1(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
             "MeshSet to coord system <%s> exist", cs_name.c_str());

  else {

    auto cs_ptr = getCoordSysPtr(cs_name);
    for (auto d : {0, 1, 2, 3})
      if (cs_ptr->getDim(d) != cs_dim[d])
        SETERRQ3(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
                 "Coord system <%s> has inconsistent dimension %d: %d != %d", d,
                 cs_ptr->getDim(d), cs_dim[d]);
  }

  MoFEMFunctionReturn(0);
}

MoFEMErrorCode
CoordSystemsManager::setFieldCoordinateSystem(const std::string field_name,
                                              const std::string cs_name) {

  Interface &m_field = cOre;
  const Field_multiIndex *fields_ptr;
  MoFEMFunctionBegin;
  CHKERR m_field.get_fields(&fields_ptr);

  // Find field
  auto field_it = fields_ptr->get<FieldName_mi_tag>().find(field_name);
  if (field_it == fields_ptr->get<FieldName_mi_tag>().end())
    SETERRQ1(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
             "Field < %s > not found", field_name.c_str());

  // Remove field from other coordinate system
  EntityHandle field_meshest = (*field_it)->getMeshset();
  for (auto &cs : coordinateSystems)
    CHKERR m_field.get_moab().remove_entities(cs->getMeshset(), &field_meshest,
                                              1);

  // Find coordinate system
  auto cs_it = coordinateSystems.get<CoordSysName_mi_tag>().find(cs_name);
  if (cs_it == coordinateSystems.get<CoordSysName_mi_tag>().end())
    SETERRQ1(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
             "Coord system < %s > not found", cs_name.c_str());

  // Add field to coordinate system
  CHKERR m_field.get_moab().add_entities((*cs_it)->getMeshset(), &field_meshest,
                                         1);

  int dim = 1;
  for (int alpha = 0; alpha < 4; alpha++) {
    if ((*cs_it)->getDim(alpha) > 0) {
      dim *= (*cs_it)->getDim(alpha);
    }
  }

  // Check consistency of field and coordinate system
  switch ((*field_it)->getSpace()) {
  case NOSPACE:
    SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY, "No space given");
  case H1:
    if ((*field_it)->getNbOfCoeffs() != dim) {
      SETERRQ2(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
               "dimension mismatch of field and coordinate system"
               "cs dim %d field rank %d",
               dim, (*field_it)->getNbOfCoeffs());
    }
    break;
  case HDIV:
  case HCURL:
    if (3 * (*field_it)->getNbOfCoeffs() != dim) {
      SETERRQ2(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
               "dimension mismatch of field and coordinate system"
               "cs dim %d field rank %d",
               dim, (*field_it)->getNbOfCoeffs());
    }
    break;
  case L2:
    if ((*field_it)->getNbOfCoeffs() != dim) {
      SETERRQ2(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
               "dimension mismatch of field and coordinate system"
               "cs dim %d field rank %d",
               dim, (*field_it)->getNbOfCoeffs());
    }
  case NOFIELD:
  case LASTSPACE:
    break;
  default:
    SETERRQ1(PETSC_COMM_SELF, MOFEM_NOT_IMPLEMENTED,
             "Not implemented for this space", (*field_it)->getSpace());
  }
  bool success = const_cast<Field_multiIndex *>(fields_ptr)
                     ->modify(fields_ptr->project<0>(field_it),
                              FieldChangeCoordinateSystem(*cs_it));
  if (!success)
    SETERRQ(PETSC_COMM_SELF, MOFEM_OPERATION_UNSUCCESSFUL,
            "modification unsuccessful");

  MoFEMFunctionReturn(0);
}

MoFEMErrorCode
CoordSystemsManager::getCoordSysPtr(const EntityHandle id,
                                    boost::shared_ptr<CoordSys> &cs_ptr) {
  MoFEMFunctionBeginHot;
  auto cs_it = coordinateSystems.get<Meshset_mi_tag>().find(id);
  if (cs_it == coordinateSystems.get<Meshset_mi_tag>().end())
    SETERRQ1(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
             "Unknown Coordinate System ms_id %lu", id);
  cs_ptr = *cs_it;
  MoFEMFunctionReturnHot(0);
}

MoFEMErrorCode
CoordSystemsManager::getCoordSysPtr(const string name,
                                    boost::shared_ptr<CoordSys> &cs_ptr) {
  MoFEMFunctionBeginHot;
  auto cs_it = coordinateSystems.get<CoordSysName_mi_tag>().find(name);
  if (cs_it == coordinateSystems.get<CoordSysName_mi_tag>().end()) {
    SETERRQ1(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
             "Unknown Coordinate System <%s>", name.c_str());
  }
  cs_ptr = *cs_it;
  MoFEMFunctionReturnHot(0);
}

boost::shared_ptr<CoordSys>
CoordSystemsManager::getCoordSysPtr(const string name) {
  MoFEMFunctionBeginHot;
  auto cs_it = coordinateSystems.get<CoordSysName_mi_tag>().find(name);
  if (cs_it == coordinateSystems.get<CoordSysName_mi_tag>().end())
    return boost::shared_ptr<CoordSys>();
  else
    return *cs_it;
  MoFEMFunctionReturnHot(0);
}

} // namespace MoFEM
