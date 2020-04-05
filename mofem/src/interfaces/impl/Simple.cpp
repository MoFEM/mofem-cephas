/** \file Simple.cpp
 * \brief Implementation of simple interface
 * \ingroup mofem_simple_interface
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

MoFEMErrorCode Simple::query_interface(const MOFEMuuid &uuid,
                                       UnknownInterface **iface) const {
  MoFEMFunctionBeginHot;
  *iface = NULL;
  if (uuid == IDD_MOFEMSimple) {
    *iface = const_cast<Simple *>(this);
    MoFEMFunctionReturnHot(0);
  }
  SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY, "unknown interface");
  MoFEMFunctionReturnHot(0);
}

template<int DIM>
MoFEMErrorCode Simple::setSkeletonAdjacency() {
  static_assert(DIM == 2 || DIM == 3, "not implemented");
  return MOFEM_NOT_IMPLEMENTED;
}

template<>
MoFEMErrorCode Simple::setSkeletonAdjacency<2>() {
  Interface &m_field = cOre;
  MoFEMFunctionBegin;

  auto defaultSkeletonEdge = [&](moab::Interface &moab, const Field &field,
                                 const EntFiniteElement &fe,
                                 Range &adjacency) -> MoFEMErrorCode {
    MoFEMFunctionBegin;

    CHKERR DefaultElementAdjacency::defaultEdge(moab, field, fe, adjacency);

    if (std::find(domainFields.begin(), domainFields.end(),
                  field.getName()) != domainFields.end()) {

      const EntityHandle fe_ent = fe.getEnt();
      Range bride_adjacency_edge, bride_adjacency;
      CHKERR moab.get_adjacencies(&fe_ent, 1, 2, false, bride_adjacency_edge);

      switch (field.getSpace()) {
      case H1:
        CHKERR moab.get_connectivity(bride_adjacency_edge, bride_adjacency,
                                     true);
      case L2:
      case HCURL:
      case HDIV:
        CHKERR moab.get_adjacencies(bride_adjacency_edge, 1, false,
                                    bride_adjacency, moab::Interface::UNION);
        bride_adjacency.merge(bride_adjacency_edge);
        break;
      case NOFIELD:
        break;
      default:
        SETERRQ(PETSC_COMM_SELF, MOFEM_NOT_IMPLEMENTED,
                "this field is not implemented for TRI finite element");
      }

      bride_adjacency = subtract(bride_adjacency, adjacency);

      for (auto e : bride_adjacency)
        const_cast<SideNumber_multiIndex &>(fe.getSideNumberTable())
            .insert(boost::shared_ptr<SideNumber>(new SideNumber(e, -1, 0, 0)));

      adjacency.merge(bride_adjacency);
    }

    MoFEMFunctionReturn(0);
  };

  CHKERR m_field.modify_finite_element_adjacency_table(skeletonFE, MBEDGE,
                                                       defaultSkeletonEdge);

  MoFEMFunctionReturn(0);
}

template<>
MoFEMErrorCode Simple::setSkeletonAdjacency<3>() {
  Interface &m_field = cOre;
  MoFEMFunctionBegin;

  auto defaultSkeletonEdge = [&](moab::Interface &moab, const Field &field,
                                 const EntFiniteElement &fe,
                                 Range &adjacency) -> MoFEMErrorCode {
    MoFEMFunctionBegin;

    CHKERR DefaultElementAdjacency::defaultFace(moab, field, fe, adjacency);

    if (std::find(domainFields.begin(), domainFields.end(),
                  field.getName()) != domainFields.end()) {

      const EntityHandle fe_ent = fe.getEnt();
      Range bride_adjacency_face, bride_adjacency;
      CHKERR moab.get_adjacencies(&fe_ent, 1, 3, false, bride_adjacency_face);

      switch (field.getSpace()) {
      case H1:
        CHKERR moab.get_connectivity(bride_adjacency_face, bride_adjacency,
                                     true);
      case HCURL:
        CHKERR moab.get_adjacencies(bride_adjacency_face, 1, false,
                                    bride_adjacency, moab::Interface::UNION);
      case L2:
      case HDIV:
        CHKERR moab.get_adjacencies(bride_adjacency_face, 2, false,
                                    bride_adjacency, moab::Interface::UNION);
        bride_adjacency.merge(bride_adjacency_face);
        break;
      case NOFIELD:
        break;
      default:
        SETERRQ(PETSC_COMM_SELF, MOFEM_NOT_IMPLEMENTED,
                "this field is not implemented for TRI finite element");
      }

      bride_adjacency = subtract(bride_adjacency, adjacency);

      for (auto e : bride_adjacency)
        const_cast<SideNumber_multiIndex &>(fe.getSideNumberTable())
            .insert(boost::shared_ptr<SideNumber>(new SideNumber(e, -1, 0, 0)));

      adjacency.merge(bride_adjacency);
    }

    MoFEMFunctionReturn(0);
  };

  CHKERR m_field.modify_finite_element_adjacency_table(skeletonFE, MBTRI,
                                                       defaultSkeletonEdge);
  CHKERR m_field.modify_finite_element_adjacency_table(skeletonFE, MBQUAD,
                                                       defaultSkeletonEdge);

  MoFEMFunctionReturn(0);
}

template<>
MoFEMErrorCode Simple::setSkeletonAdjacency<-1>() {
  MoFEMFunctionBegin;
  switch (getDim()) {
  case 1:
    THROW_MESSAGE("Not implemented");
  case 2:
    return setSkeletonAdjacency<2>();
  case 3:
    return setSkeletonAdjacency<3>(); 
  default:
    THROW_MESSAGE("Not implemented");
  }
  MoFEMFunctionReturn(0);
}

MoFEMErrorCode Simple::setSkeletonAdjacency(int dim) {
  MoFEMFunctionBegin;
  if(dim == -1)
    dim = getDim();
  switch(dim) {
  case 2:
    return setSkeletonAdjacency<2>();
  case 3:
    return setSkeletonAdjacency<3>();
  default:
    SETERRQ(PETSC_COMM_WORLD, MOFEM_NOT_IMPLEMENTED, "Not implemented");
  }
  MoFEMFunctionReturn(0);
}

Simple::Simple(const Core &core)
    : cOre(const_cast<Core &>(core)), bitLevel(BitRefLevel().set(0)),
      meshSet(0), boundaryMeshset(0), nameOfProblem("SimpleProblem"),
      domainFE("dFE"), boundaryFE("bFE"), skeletonFE("sFE"), dIm(-1) {
  PetscLogEventRegister("LoadMesh", 0, &MOFEM_EVENT_SimpleLoadMesh);
  PetscLogEventRegister("buildFields", 0, &MOFEM_EVENT_SimpleBuildFields);
  PetscLogEventRegister("buildFiniteElements", 0,
                        &MOFEM_EVENT_SimpleBuildFiniteElements);
  PetscLogEventRegister("SimpleSetUp", 0, &MOFEM_EVENT_SimpleBuildProblem);
  PetscLogEventRegister("SimpleKSPSolve", 0, &MOFEM_EVENT_SimpleKSPSolve);
  strcpy(meshFileName, "mesh.h5m");
}

MoFEMErrorCode Simple::getOptions() {
  PetscBool flg = PETSC_TRUE;
  MoFEMFunctionBeginHot;
  ierr = PetscOptionsBegin(PETSC_COMM_WORLD, "", "Simple interface options",
                           "none");
  CHKERRG(ierr);
  ierr = PetscOptionsString("-file_name", "file name", "", "mesh.h5m",
                            meshFileName, 255, &flg);
  CHKERRG(ierr);
  ierr = PetscOptionsEnd();
  CHKERRG(ierr);
  MoFEMFunctionReturnHot(0);
}

MoFEMErrorCode Simple::loadFile(const std::string options,
                                const std::string mesh_file_name) {
  Interface &m_field = cOre;
  MoFEMFunctionBegin;
  PetscLogEventBegin(MOFEM_EVENT_SimpleLoadMesh, 0, 0, 0, 0);

  if (!mesh_file_name.empty())
    strcpy(meshFileName, mesh_file_name.c_str());

  // This is a case of distributed mesh and algebra. In that case each processor
  // keep only part of the problem.
  if (m_field.get_comm_size() > 1)
    CHKERR m_field.get_moab().load_file(meshFileName, 0, options.c_str());
  else
    CHKERR m_field.get_moab().load_file(meshFileName, 0, "");
  CHKERR m_field.rebuild_database();
  // determine problem dimension
  if (dIm == -1) {
    int nb_ents_3d;
    CHKERR m_field.get_moab().get_number_entities_by_dimension(
        meshSet, 3, nb_ents_3d, true);
    if (nb_ents_3d > 0) {
      dIm = 3;
    } else {
      int nb_ents_2d;
      CHKERR m_field.get_moab().get_number_entities_by_dimension(
          meshSet, 2, nb_ents_2d, true);
      if (nb_ents_2d > 0) {
        dIm = 2;
      } else {
        dIm = 1;
      }
    }
  }
  Range ents;
  CHKERR m_field.get_moab().get_entities_by_dimension(meshSet, dIm, ents, true);
  CHKERR m_field.getInterface<BitRefManager>()->setBitRefLevel(ents, bitLevel,
                                                               false);
  ParallelComm *pcomm =
      ParallelComm::get_pcomm(&m_field.get_moab(), MYPCOMM_INDEX);
  if (pcomm == NULL)
    pcomm = new ParallelComm(&m_field.get_moab(), m_field.get_comm());
  PetscLogEventEnd(MOFEM_EVENT_SimpleLoadMesh, 0, 0, 0, 0);
  MoFEMFunctionReturn(0);
}

MoFEMErrorCode
Simple::addDomainField(const std::string &name, const FieldSpace space,
                       const FieldApproximationBase base,
                       const FieldCoefficientsNumber nb_of_cooficients,
                       const TagType tag_type, const enum MoFEMTypes bh,
                       int verb) {

  Interface &m_field = cOre;
  MoFEMFunctionBegin;
  CHKERR m_field.add_field(name, space, base, nb_of_cooficients, tag_type, bh,
                           verb);
  if (space == NOFIELD)
    noFieldFields.push_back(name);
  else
    domainFields.push_back(name);
  MoFEMFunctionReturn(0);
}

MoFEMErrorCode
Simple::addBoundaryField(const std::string &name, const FieldSpace space,
                         const FieldApproximationBase base,
                         const FieldCoefficientsNumber nb_of_cooficients,
                         const TagType tag_type, const enum MoFEMTypes bh,
                         int verb) {
  Interface &m_field = cOre;
  MoFEMFunctionBegin;
  CHKERR m_field.add_field(name, space, base, nb_of_cooficients, tag_type, bh,
                           verb);
  boundaryFields.push_back(name);
  if (space == NOFIELD)
    SETERRQ(
        PETSC_COMM_SELF, MOFEM_NOT_IMPLEMENTED,
        "NOFIELD space for boundary filed not implemented in Simple interface");
  MoFEMFunctionReturn(0);
}

MoFEMErrorCode
Simple::addSkeletonField(const std::string &name, const FieldSpace space,
                         const FieldApproximationBase base,
                         const FieldCoefficientsNumber nb_of_cooficients,
                         const TagType tag_type, const enum MoFEMTypes bh,
                         int verb) {

  Interface &m_field = cOre;
  MoFEMFunctionBegin;
  CHKERR m_field.add_field(name, space, base, nb_of_cooficients, tag_type, bh,
                           verb);
  skeletonFields.push_back(name);
  if (space == NOFIELD)
    SETERRQ(
        PETSC_COMM_SELF, MOFEM_NOT_IMPLEMENTED,
        "NOFIELD space for boundary filed not implemented in Simple interface");
  MoFEMFunctionReturn(0);
}

MoFEMErrorCode
Simple::addDataField(const std::string &name, const FieldSpace space,
                     const FieldApproximationBase base,
                     const FieldCoefficientsNumber nb_of_cooficients,
                     const TagType tag_type, const enum MoFEMTypes bh,
                     int verb) {

  Interface &m_field = cOre;
  MoFEMFunctionBegin;
  CHKERR m_field.add_field(name, space, base, nb_of_cooficients, tag_type, bh,
                           verb);
  if (space == NOFIELD)
    noFieldDataFields.push_back(name);
  else
    dataFields.push_back(name);
  MoFEMFunctionReturn(0);
}

MoFEMErrorCode Simple::removeDomainField(const std::string &name) {
  Interface &m_field = cOre;
  MoFEMFunctionBegin;

  auto remove_field_from_list = [&](auto &vec) {
    auto it = std::find(vec.begin(), vec.end(), name);
    if (it != vec.end())
      vec.erase(it);
  };

  remove_field_from_list(noFieldFields);
  remove_field_from_list(domainFields);

  MoFEMFunctionReturn(0);
}

MoFEMErrorCode Simple::removeBoundaryField(const std::string &name) {
  Interface &m_field = cOre;
  MoFEMFunctionBegin;

  auto remove_field_from_list = [&](auto &vec) {
    auto it = std::find(vec.begin(), vec.end(), name);
    if (it != vec.end())
      vec.erase(it);
  };

  remove_field_from_list(boundaryFields);

  MoFEMFunctionReturn(0);
}

MoFEMErrorCode Simple::removeSkeletonField(const std::string &name) {
  Interface &m_field = cOre;
  MoFEMFunctionBegin;

  auto remove_field_from_list = [&](auto &vec) {
    auto it = std::find(vec.begin(), vec.end(), name);
    if (it != vec.end())
      vec.erase(it);
  };

  remove_field_from_list(skeletonFields);

  MoFEMFunctionReturn(0);
}

MoFEMErrorCode Simple::defineFiniteElements() {
  Interface &m_field = cOre;
  MoFEMFunctionBeginHot;

  auto clear_rows_and_cols = [&](auto &fe_name) {
    MoFEMFunctionBeginHot;
    const FiniteElement_multiIndex *fe_ptr;
    CHKERR m_field.get_finite_elements(&fe_ptr);
    auto &fe_by_name = const_cast<FiniteElement_multiIndex *>(fe_ptr)
                           ->get<FiniteElement_name_mi_tag>();
    auto it_fe = fe_by_name.find(fe_name);
    if (it_fe != fe_by_name.end()) {

      if(!fe_by_name.modify(it_fe, FiniteElement_row_change_bit_reset()))
        SETERRQ(m_field.get_comm(), MOFEM_OPERATION_UNSUCCESSFUL,
                "modification unsuccessful");

      if(!fe_by_name.modify(it_fe, FiniteElement_col_change_bit_reset()))
        SETERRQ(m_field.get_comm(), MOFEM_OPERATION_UNSUCCESSFUL,
                "modification unsuccessful");
    }
    MoFEMFunctionReturnHot(0);
  };
  CHKERR clear_rows_and_cols(domainFE);
  CHKERR clear_rows_and_cols(boundaryFE);
  CHKERR clear_rows_and_cols(skeletonFE);

  // Define finite elements
  CHKERR m_field.add_finite_element(domainFE, MF_ZERO);

  auto add_fields = [&](auto &fe_name, auto &fields) {
    MoFEMFunctionBeginHot;
    for (auto &field : fields) {
      CHKERR m_field.modify_finite_element_add_field_row(fe_name, field);
      CHKERR m_field.modify_finite_element_add_field_col(fe_name, field);
      CHKERR m_field.modify_finite_element_add_field_data(fe_name, field);
    }
    MoFEMFunctionReturnHot(0);
  };

  auto add_data_fields = [&](auto &fe_name, auto &fields) {
    MoFEMFunctionBeginHot;
    for (auto &field : fields)
      CHKERR m_field.modify_finite_element_add_field_data(fe_name, field);
    MoFEMFunctionReturnHot(0);
  };

  CHKERR add_fields(domainFE, domainFields);
  CHKERR add_data_fields(domainFE, dataFields);
  CHKERR add_fields(domainFE, noFieldFields);
  CHKERR add_data_fields(domainFE, noFieldDataFields);

  if (!boundaryFields.empty()) {
    CHKERR m_field.add_finite_element(boundaryFE, MF_ZERO);
    CHKERR add_fields(boundaryFE, domainFields);
    CHKERR add_fields(boundaryFE, boundaryFields);
    CHKERR add_fields(boundaryFE, skeletonFields);
    CHKERR add_data_fields(boundaryFE, dataFields);
    CHKERR add_data_fields(boundaryFE, noFieldDataFields);
    CHKERR add_fields(boundaryFE, noFieldFields);
  }
  if (!skeletonFields.empty()) {
    CHKERR m_field.add_finite_element(skeletonFE, MF_ZERO);
    CHKERR add_fields(skeletonFE, domainFields);
    CHKERR add_fields(skeletonFE, boundaryFields);
    CHKERR add_fields(skeletonFE, skeletonFields);
    CHKERR add_data_fields(skeletonFE, dataFields);
    CHKERR add_data_fields(skeletonFE, noFieldDataFields);
    CHKERR add_fields(skeletonFE, noFieldFields);
  }
  MoFEMFunctionReturnHot(0);
}

MoFEMErrorCode Simple::defineProblem(const PetscBool is_partitioned) {
  Interface &m_field = cOre;
  MoFEMFunctionBegin;
  // Create dm instance
  dM = createSmartDM(m_field.get_comm(), "DMMOFEM");
  // set dm data structure which created mofem data structures
  CHKERR DMMoFEMCreateMoFEM(dM, &m_field, nameOfProblem.c_str(), bitLevel);
  CHKERR DMSetFromOptions(dM);
  CHKERR DMMoFEMAddElement(dM, domainFE.c_str());
  if (!boundaryFields.empty()) {
    CHKERR DMMoFEMAddElement(dM, boundaryFE.c_str());
  }
  if (!skeletonFields.empty()) {
    CHKERR DMMoFEMAddElement(dM, skeletonFE.c_str());
  }
  for (std::vector<std::string>::iterator fit = otherFEs.begin();
       fit != otherFEs.end(); ++fit) {
    CHKERR DMMoFEMAddElement(dM, fit->c_str());
  }
  CHKERR DMMoFEMSetIsPartitioned(dM, is_partitioned);
  MoFEMFunctionReturn(0);
}

MoFEMErrorCode Simple::setFieldOrder(const std::string field_name,
                                     const int order, const Range *ents) {
  MoFEMFunctionBeginHot;
  fieldsOrder.insert(
      {field_name, {order, ents == NULL ? Range() : Range(*ents)}});
  MoFEMFunctionReturnHot(0);
}

MoFEMErrorCode Simple::buildFields() {
  Interface &m_field = cOre;
  MoFEMFunctionBegin;
  PetscLogEventBegin(MOFEM_EVENT_SimpleBuildFields, 0, 0, 0, 0);

  auto get_skin = [&](auto meshset) {
    Range domain_ents;
    CHKERR m_field.get_moab().get_entities_by_dimension(meshset, dIm,
                                                        domain_ents, true);
    Skinner skin(&m_field.get_moab());
    Range domain_skin;
    CHKERR skin.find_skin(0, domain_ents, false, domain_skin);
    // filter not owned entities, those are not on boundary
    ParallelComm *pcomm =
        ParallelComm::get_pcomm(&m_field.get_moab(), MYPCOMM_INDEX);
    Range proc_domain_skin;
    CHKERR pcomm->filter_pstatus(domain_skin,
                                 PSTATUS_SHARED | PSTATUS_MULTISHARED,
                                 PSTATUS_NOT, -1, &proc_domain_skin);
    return proc_domain_skin;
  };

  auto create_boundary_meshset = [&](auto &&proc_domain_skin) {
    MoFEMFunctionBeginHot;
    // create boundary meshset
    if (boundaryMeshset != 0) {
      CHKERR m_field.get_moab().delete_entities(&boundaryMeshset, 1);
    }
    CHKERR m_field.get_moab().create_meshset(MESHSET_SET, boundaryMeshset);
    CHKERR m_field.get_moab().add_entities(boundaryMeshset, proc_domain_skin);
    for (int dd = 0; dd != dIm - 1; dd++) {
      Range adj;
      CHKERR m_field.get_moab().get_adjacencies(proc_domain_skin, dd, false,
                                                adj, moab::Interface::UNION);
      CHKERR m_field.get_moab().add_entities(boundaryMeshset, adj);
    }
    MoFEMFunctionReturnHot(0);
  };

  CHKERR create_boundary_meshset(get_skin(meshSet));

  auto comm_interface_ptr = m_field.getInterface<CommInterface>();
  auto bit_ref_ptr = m_field.getInterface<BitRefManager>();

  // Add entities to the fields
  auto add_ents_to_field = [&](auto meshset, auto dim, auto &fields) {
    MoFEMFunctionBeginHot;
    for (auto &field : fields) {
      CHKERR m_field.add_ents_to_field_by_dim(meshset, dim, field);
      CHKERR comm_interface_ptr->synchroniseFieldEntities(field, 0);
    }
    MoFEMFunctionReturnHot(0);
  };

  auto make_no_field_ents = [&](auto &fields) {
    MoFEMFunctionBeginHot;
    for (auto &field : fields) {
      std::array<double, 6> coords = {0, 0, 0, 0, 0, 0};
      CHKERR m_field.create_vertices_and_add_to_field(field, coords.data(), 2);
      CHKERR comm_interface_ptr->makeFieldEntitiesMultishared(field, 0);
      CHKERR bit_ref_ptr->setFieldEntitiesBitRefLevel(field, bitLevel);
    }
    MoFEMFunctionReturnHot(0);
  };

  CHKERR add_ents_to_field(meshSet, dIm, domainFields);
  CHKERR add_ents_to_field(meshSet, dIm, dataFields);
  CHKERR add_ents_to_field(boundaryMeshset, dIm - 1, boundaryFields);
  CHKERR add_ents_to_field(meshSet, dIm - 1, skeletonFields);
  CHKERR make_no_field_ents(noFieldFields);
  CHKERR make_no_field_ents(noFieldDataFields);

  // Set order
  auto set_order = [&](auto meshset, auto dim, auto &fields) {
    MoFEMFunctionBeginHot;
    for (auto &field : fields) {
      if (fieldsOrder.find(field) == fieldsOrder.end()) {
        SETERRQ1(PETSC_COMM_WORLD, MOFEM_INVALID_DATA,
                 "Order for field not set %s", field.c_str());
      }
      int dds = 0;
      const Field *field_ptr = m_field.get_field_structure(field);
      switch (field_ptr->getSpace()) {
      case L2:
        dds = dim;
        break;
      case HDIV:
        dds = 2;
        break;
      case HCURL:
        dds = 1;
        break;
      case H1:
        dds = 1;
        break;
      default:
        SETERRQ(PETSC_COMM_WORLD, MOFEM_DATA_INCONSISTENCY,
                "Glasgow we have a problem");
      }

      auto set_order = [&](auto field, auto &ents) {
        MoFEMFunctionBegin;

        auto range = fieldsOrder.equal_range(field);
        for (auto o = range.first; o != range.second; ++o) {
          if (!o->second.second.empty())
            ents = intersect(ents, o->second.second);
          CHKERR m_field.set_field_order(ents, field, o->second.first);
        }

        MoFEMFunctionReturn(0);
      };

      if (field_ptr->getSpace() == H1) {
        if (field_ptr->getApproxBase() == AINSWORTH_BERNSTEIN_BEZIER_BASE) {
          Range ents;
          CHKERR m_field.get_field_entities_by_dimension(field, 0, ents);
          CHKERR set_order(field, ents);
        } else {
          CHKERR m_field.set_field_order(meshSet, MBVERTEX, field, 1);
        }
      }
      for (int dd = dds; dd <= dim; dd++) {
        Range ents;
        CHKERR m_field.get_field_entities_by_dimension(field, dd, ents);
        CHKERR set_order(field, ents);
      }
    }
    MoFEMFunctionReturnHot(0);
  };

  CHKERR set_order(meshSet, dIm, domainFields);
  CHKERR set_order(meshSet, dIm, dataFields);
  CHKERR set_order(meshSet, dIm - 1, boundaryFields);
  CHKERR set_order(meshSet, dIm - 1, skeletonFields);

  // Build fields
  CHKERR m_field.build_fields();
  PetscLogEventEnd(MOFEM_EVENT_SimpleBuildFields, 0, 0, 0, 0);
  MoFEMFunctionReturn(0);
}

MoFEMErrorCode Simple::buildFiniteElements() {
  Interface &m_field = cOre;
  MoFEMFunctionBegin;
  PetscLogEventBegin(MOFEM_EVENT_SimpleBuildFiniteElements, 0, 0, 0, 0);
  // Add finite elements
  CHKERR m_field.add_ents_to_finite_element_by_dim(meshSet, dIm, domainFE,
                                                   true);
  CHKERR m_field.build_finite_elements(domainFE);
  if (!boundaryFields.empty()) {
    CHKERR m_field.add_ents_to_finite_element_by_dim(boundaryMeshset, dIm - 1,
                                                     boundaryFE, true);
    CHKERR m_field.build_finite_elements(boundaryFE);
  }
  if (!skeletonFields.empty()) {
    CHKERR m_field.add_ents_to_finite_element_by_dim(meshSet, dIm - 1,
                                                     skeletonFE, true);
    CHKERR m_field.build_finite_elements(skeletonFE);
  }
  for (std::vector<std::string>::iterator fit = otherFEs.begin();
       fit != otherFEs.end(); ++fit) {
    CHKERR m_field.build_finite_elements(*fit);
  }
  PetscLogEventEnd(MOFEM_EVENT_SimpleBuildFiniteElements, 0, 0, 0, 0);
  MoFEMFunctionReturn(0);
}

MoFEMErrorCode Simple::buildProblem() {
  Interface &m_field = cOre;
  MoFEMFunctionBegin;
  PetscLogEventBegin(MOFEM_EVENT_SimpleBuildProblem, 0, 0, 0, 0);
  CHKERR m_field.build_adjacencies(bitLevel);
  // Set problem by the DOFs on the fields rather that by adding DOFs on the
  // elements
  m_field.getInterface<ProblemsManager>()->buildProblemFromFields = PETSC_TRUE;
  CHKERR DMSetUp(dM);
  m_field.getInterface<ProblemsManager>()->buildProblemFromFields = PETSC_FALSE;
  PetscLogEventEnd(MOFEM_EVENT_SimpleBuildProblem, 0, 0, 0, 0);
  MoFEMFunctionReturn(0);
}

MoFEMErrorCode Simple::setUp(const PetscBool is_partitioned) {
  MoFEMFunctionBegin;
  CHKERR defineFiniteElements();
  if (!skeletonFields.empty())
    CHKERR setSkeletonAdjacency();
  CHKERR defineProblem(is_partitioned);
  CHKERR buildFields();
  CHKERR buildFiniteElements();
  CHKERR buildProblem();
  MoFEMFunctionReturn(0);
}

MoFEMErrorCode Simple::reSetUp() {
  Interface &m_field = cOre;
  MoFEMFunctionBegin;

  CHKERR defineFiniteElements();
  if (!skeletonFields.empty())
    CHKERR setSkeletonAdjacency();
  CHKERR buildFields();
  CHKERR buildFiniteElements();

  PetscLogEventBegin(MOFEM_EVENT_SimpleBuildProblem, 0, 0, 0, 0);
  CHKERR m_field.build_adjacencies(bitLevel);
  // Set problem by the DOFs on the fields rather that by adding DOFs on the
  // elements
  m_field.getInterface<ProblemsManager>()->buildProblemFromFields = PETSC_TRUE;
  CHKERR DMSetUp_MoFEM(dM);
  m_field.getInterface<ProblemsManager>()->buildProblemFromFields = PETSC_FALSE;
  PetscLogEventEnd(MOFEM_EVENT_SimpleBuildProblem, 0, 0, 0, 0);

  MoFEMFunctionReturn(0);
}

MoFEMErrorCode Simple::getDM(DM *dm) {
  MoFEMFunctionBegin;
  CHKERR PetscObjectReference(getPetscObject(dM.get()));
  *dm = dM.get();
  MoFEMFunctionReturn(0);
}

} // namespace MoFEM
