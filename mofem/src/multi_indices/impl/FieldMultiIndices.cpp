/** \file FieldMultiIndices.cpp
 * \brief Multi-index containers for fields
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

// Not partitioned
const bool Idx_mi_tag::IamNotPartitioned = true;

// This tag is used for partitioned problems
const bool PetscGlobalIdx_mi_tag::IamNotPartitioned = false;
const bool PetscLocalIdx_mi_tag::IamNotPartitioned = false;

// fields
Field::Field(const moab::Interface &moab, const EntityHandle meshset,
             const boost::shared_ptr<CoordSys> coord_sys_ptr)
    : moab(const_cast<moab::Interface &>(moab)), meshSet(meshset),
      coordSysPtr(coord_sys_ptr), tagId(NULL), tagSpaceData(NULL),
      tagNbCoeffData(NULL), tagName(NULL), tagNameSize(0) {

  auto get_tag_data_ptr = [&](const auto name, auto &tag_data) {
    MoFEMFunctionBegin;
    Tag th;
    CHKERR moab.tag_get_handle(name, th);
    CHKERR moab.tag_get_by_ptr(th, &meshset, 1, (const void **)&tag_data);
    MoFEMFunctionReturn(0);
  };

  // id
  ierr = get_tag_data_ptr("_FieldId", tagId);
  CHKERRABORT(PETSC_COMM_SELF, ierr);
  // space
  ierr = get_tag_data_ptr("_FieldSpace", tagSpaceData);
  CHKERRABORT(PETSC_COMM_SELF, ierr);

  // approx. base
  ierr = get_tag_data_ptr("_FieldBase", tagBaseData);
  CHKERRABORT(PETSC_COMM_SELF, ierr);

  auto get_all_tags = [&]() {
    MoFEMFunctionBegin;
    // name
    Tag th_field_name;
    CHKERR moab.tag_get_handle("_FieldName", th_field_name);
    CHKERR moab.tag_get_by_ptr(th_field_name, &meshSet, 1,
                               (const void **)&tagName, &tagNameSize);
    // name prefix
    Tag th_field_name_data_name_prefix;
    CHKERR moab.tag_get_handle("_FieldName_DataNamePrefix",
                               th_field_name_data_name_prefix);
    CHKERR moab.tag_get_by_ptr(th_field_name_data_name_prefix, &meshSet, 1,
                               (const void **)&tagNamePrefixData,
                               &tagNamePrefixSize);
    std::string name_data_prefix((char *)tagNamePrefixData, tagNamePrefixSize);
    // data
    std::string tag_data_name = name_data_prefix + getName();
    CHKERR moab.tag_get_handle(tag_data_name.c_str(), th_FieldData);
    std::string tag_data_name_verts = name_data_prefix + getName() + "V";
    CHKERR moab.tag_get_handle(tag_data_name_verts.c_str(), th_FieldDataVerts);
    CHKERR moab.tag_get_type(th_FieldDataVerts, th_FieldDataVertsType);
    // order
    std::string tag_approximation_order_name = "_App_Order_" + getName();
    CHKERR moab.tag_get_handle(tag_approximation_order_name.c_str(),
                               th_AppOrder);
    // rank
    Tag th_rank;
    std::string Tag_rank_name = "_Field_Rank_" + getName();
    CHKERR moab.tag_get_handle(Tag_rank_name.c_str(), th_rank);
    CHKERR moab.tag_get_by_ptr(th_rank, &meshSet, 1,
                               (const void **)&tagNbCoeffData);
    MoFEMFunctionReturn(0);
  };

  ierr = get_all_tags();
  CHKERRABORT(PETSC_COMM_SELF, ierr);

  bitNumber = getBitNumberCalculate();

  auto reset_entity_order_table = [&]() {
    for (int tt = 0; tt != MBMAXTYPE; ++tt) 
      forderTable[tt] = NULL;
  };

  auto set_entity_order_table = [&]() {
    switch (*tagBaseData) {
    case AINSWORTH_LEGENDRE_BASE:
    case AINSWORTH_LOBATTO_BASE:
      switch (*tagSpaceData) {
      case H1:
        forderTable[MBVERTEX] = fNBVERTEX_H1;
        forderTable[MBEDGE] = fNBEDGE_H1;
        forderTable[MBTRI] = fNBFACETRI_H1;
        forderTable[MBQUAD] = fNBFACEQUAD_H1;
        forderTable[MBTET] = fNBVOLUMETET_H1;
        forderTable[MBPRISM] = fNBVOLUMEPRISM_H1;
        break;
      case HCURL:
        forderTable[MBVERTEX] = fNBVERTEX_HCURL;
        forderTable[MBEDGE] = fNBEDGE_AINSWORTH_HCURL;
        forderTable[MBTRI] = fNBFACETRI_AINSWORTH_HCURL;
        forderTable[MBTET] = fNBVOLUMETET_AINSWORTH_HCURL;
        break;
      case HDIV:
        forderTable[MBVERTEX] = fNBVERTEX_HDIV;
        forderTable[MBEDGE] = fNBEDGE_HDIV;
        forderTable[MBTRI] = fNBFACETRI_AINSWORTH_HDIV;
        forderTable[MBTET] = fNBVOLUMETET_AINSWORTH_HDIV;
        break;
      case L2:
        forderTable[MBVERTEX] = fNBVERTEX_L2;
        forderTable[MBEDGE] = fNBEDGE_L2;
        forderTable[MBTRI] = fNBFACETRI_L2;
        forderTable[MBTET] = fNBVOLUMETET_L2;
        break;
      default:
        THROW_MESSAGE("unknown approximation space");
      }
      break;
    case AINSWORTH_BERNSTEIN_BEZIER_BASE:
      switch (*tagSpaceData) {
      case H1:
        forderTable[MBVERTEX] = fNBVERTEX_H1;
        forderTable[MBEDGE] = fNBEDGE_H1;
        forderTable[MBTRI] = fNBFACETRI_H1;
        forderTable[MBQUAD] = fNBFACEQUAD_H1;
        forderTable[MBTET] = fNBVOLUMETET_H1;
        forderTable[MBPRISM] = fNBVOLUMEPRISM_H1;
        break;
      case L2:
        forderTable[MBVERTEX] = fNBVERTEX_L2;
        forderTable[MBEDGE] = fNBEDGE_L2;
        forderTable[MBTRI] = fNBFACETRI_L2;
        forderTable[MBTET] = fNBVOLUMETET_L2;
        break;
      default:
        THROW_MESSAGE("unknown approximation space or not yet implemented");
      }
      break;
    case DEMKOWICZ_JACOBI_BASE:
      switch (*tagSpaceData) {
      case HCURL:
        forderTable[MBVERTEX] = fNBVERTEX_HCURL;
        forderTable[MBEDGE] = fNBEDGE_DEMKOWICZ_HCURL;
        forderTable[MBTRI] = fNBFACETRI_DEMKOWICZ_HCURL;
        forderTable[MBTET] = fNBVOLUMETET_DEMKOWICZ_HCURL;
        break;
      case HDIV:
        forderTable[MBVERTEX] = fNBVERTEX_HDIV;
        forderTable[MBEDGE] = fNBEDGE_HDIV;
        forderTable[MBTRI] = fNBFACETRI_DEMKOWICZ_HDIV;
        forderTable[MBTET] = fNBVOLUMETET_DEMKOWICZ_HDIV;
        break;
      default:
        THROW_MESSAGE("unknown approximation space or not yet implemented");
      }
      break;
    case USER_BASE:
      for (int ee = 0; ee < MBMAXTYPE; ee++) {
        forderTable[ee] = fNBENTITY_GENERIC;
      }
      break;
    default:
      if (*tagSpaceData != NOFIELD) {
        THROW_MESSAGE("unknown approximation base");
      } else {
        for (EntityType t = MBVERTEX; t < MBMAXTYPE; t++) 
          forderTable[t] = fNBENTITYSET_NOFIELD;
      }
    }
  };

  reset_entity_order_table();
  set_entity_order_table();
  ierr = rebuildDofsOrderMap();
  CHKERRABORT(PETSC_COMM_SELF, ierr);
};

MoFEMErrorCode Field::rebuildDofsOrderMap() const {
  MoFEMFunctionBegin;

  for (auto t = MBVERTEX; t != MBMAXTYPE; ++t) {

    int DD = 0;
    int nb_last_order_dofs = 0;
    const int rank = (*tagNbCoeffData);

    if (forderTable[t]) {

      for (int oo = 0; oo < MAX_DOFS_ON_ENTITY; ++oo) {

        const int nb_order_dofs = forderTable[t](oo);
        const int diff_oo = nb_order_dofs - nb_last_order_dofs;
        if (diff_oo >= 0) {

          if ((DD + rank * diff_oo) < MAX_DOFS_ON_ENTITY)
            for (int dd = 0; dd < diff_oo; ++dd)
              for (int rr = 0; rr != rank; ++rr, ++DD)
                dofOrderMap[t][DD] = oo;
          else
            break;

          nb_last_order_dofs = nb_order_dofs;

        } else {
          break;
        }
      }
    }

    std::fill(&dofOrderMap[t][DD], dofOrderMap[t].end(), -1);
  }

  MoFEMFunctionReturn(0);
}

std::ostream &operator<<(std::ostream &os, const Field &e) {
  os << "name " << e.getNameRef() << " BitFieldId " << e.getId().to_ulong()
     << " bit number " << e.getBitNumber() << " space "
     << FieldSpaceNames[e.getSpace()] << " approximation base "
     << ApproximationBaseNames[e.getApproxBase()] << " rank "
     << e.getNbOfCoeffs() << " meshset " << e.meshSet;
  return os;
}

// FieldEntityEntFiniteElementAdjacencyMap
FieldEntityEntFiniteElementAdjacencyMap::
    FieldEntityEntFiniteElementAdjacencyMap(
        const boost::shared_ptr<FieldEntity> &ent_field_ptr,
        const boost::shared_ptr<EntFiniteElement> &ent_fe_ptr)
    : byWhat(0), entFieldPtr(ent_field_ptr), entFePtr(ent_fe_ptr) {}

std::ostream &operator<<(std::ostream &os,
                         const FieldEntityEntFiniteElementAdjacencyMap &e) {
  os << "byWhat " << std::bitset<3>(e.byWhat) << " " << *e.entFieldPtr
     << std::endl
     << *e.entFePtr->sFePtr;
  return os;
}

MoFEMErrorCode test_moab(moab::Interface &moab, const EntityHandle ent) {
  MoFEMFunctionBeginHot;
  // tets type
  EntityType type = (EntityType)((ent & MB_TYPE_MASK) >> MB_ID_WIDTH);
  if (type != moab.type_from_handle(ent))
    SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
            "inconsistencies with type_from_handle");
  // tets id
  EntityID id = (EntityType)(ent & MB_ID_MASK);
  if (id != moab.id_from_handle(ent))
    SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
            "inconsistencies with id_from_handle");
  MoFEMFunctionReturnHot(0);
}

} // namespace MoFEM
