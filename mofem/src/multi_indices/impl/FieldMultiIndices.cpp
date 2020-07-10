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

constexpr const int FieldTmp<0, 0>::CoreValue;
constexpr const int FieldTmp<0, 0>::FieldValue;
constexpr const int FieldTmp<-1, -1>::CoreValue;
constexpr const int FieldTmp<-1, -1>::FieldValue;

// fields
Field::FieldTmp(const moab::Interface &moab, const EntityHandle meshset,
                const boost::shared_ptr<CoordSys> coord_sys_ptr)
    : moab(const_cast<moab::Interface &>(moab)), meshSet(meshset),
      coordSysPtr(coord_sys_ptr), tagId(NULL), tagSpaceData(NULL),
      tagNbCoeffData(NULL), tagName(NULL), tagNameSize(0),
      destructorCalled(false) {

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
        forderTable[MBVERTEX] = [](int P) -> int { return (P > 0) ? 1 : 0; };
        forderTable[MBEDGE] = [](int P) -> int { return NBEDGE_H1(P); };
        forderTable[MBTRI] = [](int P) -> int { return NBFACETRI_H1(P); };
        forderTable[MBQUAD] = [](int P) -> int { return NBFACEQUAD_H1(P); };
        forderTable[MBTET] = [](int P) -> int { return NBVOLUMETET_H1(P); };
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
        forderTable[MBTET] = [](int P) -> int { return NBVOLUMETET_L2(P); };
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
        forderTable[MBTET] = [](int P) -> int { return NBVOLUMETET_L2(P); };
        break;
      default:
        THROW_MESSAGE("unknown approximation space or not yet implemented");
      }
      break;
    case DEMKOWICZ_JACOBI_BASE:
      switch (*tagSpaceData) {
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
        forderTable[MBTET] = [](int P) -> int {
          return NBVOLUMETET_DEMKOWICZ_HDIV(P);
        };
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

FieldTmp<0, 0>::~FieldTmp() {
  if (!destructorCalled) 
    FieldEntityTmp<0, 0>::sFieldPtr.reset();
  
  destructorCalled = true;
}

FieldTmp<-1, -1>::~FieldTmp() { this->destructorCalled = true; }

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
