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
Field::Field(moab::Interface &moab, const EntityHandle meshset)
    : moab(moab), meshSet(meshset), tagId(NULL), tagSpaceData(NULL),
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

  // rank
  std::string Tag_rank_name = "_Field_Rank_" + getName();
  CHKERR moab.tag_get_handle(Tag_rank_name.c_str(), th_FieldRank);
  CHKERR moab.tag_get_by_ptr(th_FieldRank, &meshSet, 1,
                             (const void **)&tagNbCoeffData);

  auto get_all_tags = [&]() {
    MoFEMFunctionBegin;
    // order
    ApproximationOrder def_approx_order = -1;
    std::string tag_approximation_order_name = "_App_Order_" + getName();
    rval = moab.tag_get_handle(tag_approximation_order_name.c_str(), 1,
                               MB_TYPE_INTEGER, th_AppOrder,
                               MB_TAG_CREAT | MB_TAG_SPARSE, &def_approx_order);
    if (rval == MB_ALREADY_ALLOCATED)
      rval = MB_SUCCESS;
    MOAB_THROW(rval);

    // data
    std::string tag_data_name = name_data_prefix + getName();
    const int def_len = 0;
    rval = moab.tag_get_handle(
        tag_data_name.c_str(), def_len, MB_TYPE_DOUBLE, th_FieldData,
        MB_TAG_CREAT | MB_TAG_VARLEN | MB_TAG_SPARSE, NULL);
    if (rval == MB_ALREADY_ALLOCATED)
      rval = MB_SUCCESS;
    MOAB_THROW(rval);

    std::string tag_data_name_verts = name_data_prefix + getName() + "_V";
    rval = moab.tag_get_handle(tag_data_name_verts.c_str(), th_FieldDataVerts);
    if (rval == MB_SUCCESS)
      CHKERR moab.tag_get_type(th_FieldDataVerts, tagFieldDataVertsType);
    else {
      // Since vertex tag is not it mesh that tag is not dense, it is sparse,
      // sinc it is set to all vertices on the mesh. Is unlikely that mesh has
      // no vertices, then above assumption does not hold.
      tagFieldDataVertsType = MB_TAG_SPARSE;
      VectorDouble def_vert_data(*tagNbCoeffData);
      def_vert_data.clear();
      rval = moab.tag_get_handle(tag_data_name_verts.c_str(), *tagNbCoeffData,
                                 MB_TYPE_DOUBLE, th_FieldDataVerts,
                                 MB_TAG_CREAT | tagFieldDataVertsType,
                                 &*def_vert_data.begin());
      if (rval == MB_ALREADY_ALLOCATED)
        rval = MB_SUCCESS;
      MOAB_THROW(rval);
    }

    MoFEMFunctionReturn(0);
  };

  auto get_all_tags_deprecated = [&]() {
    MoFEMFunctionBegin;
    // order
    ApproximationOrder def_approx_order = -1;
    std::string tag_approximation_order_name = "_App_Order_" + getName();
    rval = moab.tag_get_handle(tag_approximation_order_name.c_str(), 1,
                               MB_TYPE_INTEGER, th_AppOrder,
                               MB_TAG_CREAT | MB_TAG_SPARSE, &def_approx_order);
    if (rval == MB_ALREADY_ALLOCATED)
      rval = MB_SUCCESS;
    MOAB_THROW(rval);

    // data
    std::string tag_data_name = name_data_prefix + getName();
    const int def_len = 0;
    rval = moab.tag_get_handle(
        tag_data_name.c_str(), def_len, MB_TYPE_DOUBLE, th_FieldData,
        MB_TAG_CREAT | MB_TAG_VARLEN | MB_TAG_SPARSE, NULL);
    if (rval == MB_ALREADY_ALLOCATED)
      rval = MB_SUCCESS;
    MOAB_THROW(rval);

    std::string tag_data_name_verts = name_data_prefix + getName() + "V";
    rval = moab.tag_get_handle(tag_data_name_verts.c_str(), th_FieldDataVerts);
    if (rval == MB_SUCCESS)
      CHKERR moab.tag_get_type(th_FieldDataVerts, tagFieldDataVertsType);
    else {
      // Since vertex tag is not it mesh that tag is not dense, it is sparse,
      // sinc it is set to all vertices on the mesh. Is unlikely that mesh has
      // no vertices, then above assumption does not hold.
      tagFieldDataVertsType = MB_TAG_SPARSE;
      VectorDouble def_vert_data(*tagNbCoeffData);
      def_vert_data.clear();
      rval = moab.tag_get_handle(tag_data_name_verts.c_str(), *tagNbCoeffData,
                                 MB_TYPE_DOUBLE, th_FieldDataVerts,
                                 MB_TAG_CREAT | tagFieldDataVertsType,
                                 &*def_vert_data.begin());
      if (rval == MB_ALREADY_ALLOCATED)
        rval = MB_SUCCESS;
      MOAB_THROW(rval);
    }

    MoFEMFunctionReturn(0);
  };

  Version file_ver;
  ierr = UnknownInterface::getFileVersion(moab, file_ver);
  CHK_THROW_MESSAGE(ierr, "Not known file version");
  if (file_ver.majorVersion >= 0 && file_ver.minorVersion >= 12 &&
      file_ver.buildVersion >= 1) {
    ierr = get_all_tags();
    CHKERRABORT(PETSC_COMM_SELF, ierr);
  } else {
    ierr = get_all_tags_deprecated();
    CHKERRABORT(PETSC_COMM_SELF, ierr);
  }

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
        forderTable[MBVERTEX] = [](int P) -> int { return (P > 0) ? 1 : 0; };
        forderTable[MBEDGE] = [](int P) -> int { return NBEDGE_H1(P); };
        forderTable[MBTRI] = [](int P) -> int { return NBFACETRI_H1(P); };
        forderTable[MBQUAD] = [](int P) -> int { return NBFACEQUAD_H1(P); };
        forderTable[MBTET] = [](int P) -> int { return NBVOLUMETET_H1(P); };
        forderTable[MBHEX] = [](int P) -> int { return NBVOLUMEHEX_H1(P); };
        forderTable[MBPRISM] = [](int P) -> int { return NBVOLUMEPRISM_H1(P); };
        break;
      case HCURL:
        forderTable[MBVERTEX] = [](int P) -> int {
          (void)P;
          return 0;
        };
        forderTable[MBEDGE] = [](int P) -> int {
          return NBEDGE_AINSWORTH_HCURL(P);
        };
        forderTable[MBTRI] = [](int P) -> int {
          return NBFACETRI_AINSWORTH_HCURL(P);
        };
        forderTable[MBTET] = [](int P) -> int {
          return NBVOLUMETET_AINSWORTH_HCURL(P);
        };
        break;
      case HDIV:
        forderTable[MBVERTEX] = [](int P) -> int {
          (void)P;
          return 0;
        };
        forderTable[MBEDGE] = [](int P) -> int {
          (void)P;
          return NBEDGE_HDIV(P);
        };
        forderTable[MBTRI] = [](int P) -> int {
          return NBFACETRI_AINSWORTH_HDIV(P);
        };
        forderTable[MBTET] = [](int P) -> int {
          return NBVOLUMETET_AINSWORTH_HDIV(P);
        };
        break;
      case L2:
        forderTable[MBVERTEX] = [](int P) -> int {
          (void)P;
          return 1;
        };
        forderTable[MBEDGE] = [](int P) -> int { return NBEDGE_L2(P); };
        forderTable[MBTRI] = [](int P) -> int { return NBFACETRI_L2(P); };
        forderTable[MBQUAD] = [](int P) -> int { return NBFACEQUAD_L2(P); };
        forderTable[MBTET] = [](int P) -> int { return NBVOLUMETET_L2(P); };
        forderTable[MBHEX] = [](int P) -> int { return NBVOLUMEHEX_L2(P); };
        break;
      default:
        THROW_MESSAGE("unknown approximation space");
      }
      break;
    case AINSWORTH_BERNSTEIN_BEZIER_BASE:
      switch (*tagSpaceData) {
      case H1:
        forderTable[MBVERTEX] = [](int P) -> int { return (P > 0) ? 1 : 0; };
        forderTable[MBEDGE] = [](int P) -> int { return NBEDGE_H1(P); };
        forderTable[MBTRI] = [](int P) -> int { return NBFACETRI_H1(P); };
        forderTable[MBQUAD] = [](int P) -> int { return NBFACEQUAD_H1(P); };
        forderTable[MBTET] = [](int P) -> int { return NBVOLUMETET_H1(P); };
        forderTable[MBPRISM] = [](int P) -> int { return NBVOLUMEPRISM_H1(P); };
        break;
      case L2:
        forderTable[MBVERTEX] = [](int P) -> int {
          (void)P;
          return 1;
        };
        forderTable[MBEDGE] = [](int P) -> int { return NBEDGE_L2(P); };
        forderTable[MBTRI] = [](int P) -> int { return NBFACETRI_L2(P); };
        forderTable[MBQUAD] = [](int P) -> int { return NBFACEQUAD_L2(P); };
        forderTable[MBTET] = [](int P) -> int { return NBVOLUMETET_L2(P); };
        forderTable[MBHEX] = [](int P) -> int { return NBVOLUMEHEX_L2(P); };
        break;
      default:
        THROW_MESSAGE("unknown approximation space or not yet implemented");
      }
      break;
    case DEMKOWICZ_JACOBI_BASE:
      switch (*tagSpaceData) {
      case H1:
        forderTable[MBVERTEX] = [](int P) -> int { return (P > 0) ? 1 : 0; };
        forderTable[MBEDGE] = [](int P) -> int { return NBEDGE_H1(P); };
        forderTable[MBTRI] = [](int P) -> int { return NBFACETRI_H1(P); };
        forderTable[MBQUAD] = [](int P) -> int { return NBFACEQUAD_H1(P); };
        forderTable[MBTET] = [](int P) -> int { return NBVOLUMETET_H1(P); };
        forderTable[MBHEX] = [](int P) -> int { return NBVOLUMEHEX_H1(P); };
        forderTable[MBPRISM] = [](int P) -> int { return NBVOLUMEPRISM_H1(P); };
        break;
      case HCURL:
        forderTable[MBVERTEX] = [](int P) -> int {
          (void)P;
          return 0;
        };
        forderTable[MBEDGE] = [](int P) -> int {
          return NBEDGE_DEMKOWICZ_HCURL(P);
        };
        forderTable[MBTRI] = [](int P) -> int {
          return NBFACETRI_DEMKOWICZ_HCURL(P);
        };
        forderTable[MBQUAD] = [](int P) -> int {
          return NBFACEQUAD_DEMKOWICZ_HCURL(P);
        };
        forderTable[MBTET] = [](int P) -> int {
          return NBVOLUMETET_DEMKOWICZ_HCURL(P);
        };
        break;
      case HDIV:
        forderTable[MBVERTEX] = [](int P) -> int {
          (void)P;
          return 0;
        };
        forderTable[MBEDGE] = [](int P) -> int {
          (void)P;
          return 0;
        };
        forderTable[MBTRI] = [](int P) -> int {
          return NBFACETRI_DEMKOWICZ_HDIV(P);
        };
        forderTable[MBQUAD] = [](int P) -> int {
          return NBFACEQUAD_DEMKOWICZ_HDIV(P);
        };
        forderTable[MBTET] = [](int P) -> int {
          return NBVOLUMETET_DEMKOWICZ_HDIV(P);
        };
        forderTable[MBHEX] = [](int P) -> int {
          return NBVOLUMEHEX_DEMKOWICZ_HDIV(P);
        };
        break;
      case L2:
        forderTable[MBVERTEX] = [](int P) -> int {
          (void)P;
          return 1;
        };
        forderTable[MBEDGE] = [](int P) -> int { return NBEDGE_L2(P); };
        forderTable[MBTRI] = [](int P) -> int { return NBFACETRI_L2(P); };
        forderTable[MBQUAD] = [](int P) -> int { return NBFACEQUAD_L2(P); };
        forderTable[MBTET] = [](int P) -> int { return NBVOLUMETET_L2(P); };
        forderTable[MBHEX] = [](int P) -> int { return NBVOLUMEHEX_L2(P); };
        break;
      default:
        THROW_MESSAGE("unknown approximation space or not yet implemented");
      }
      break;
    case USER_BASE:
      for (int ee = 0; ee < MBMAXTYPE; ee++) {
        forderTable[ee] = [](int P) -> int {
          (void)P;
          return 0;
        };
      }
      break;
    default:
      if (*tagSpaceData != NOFIELD) {
        THROW_MESSAGE("unknown approximation base");
      } else {
        for (EntityType t = MBVERTEX; t < MBMAXTYPE; t++)
          forderTable[t] = [](int P) -> int {
            (void)P;
            return 1;
          };
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
  os << e.getNameRef() << " field_id " << e.getId().to_ulong() << " space "
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
     << *e.entFePtr->getFiniteElementPtr();
  return os;
}

} // namespace MoFEM
