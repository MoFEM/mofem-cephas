/** \file Simple.cpp
 * \brief Implementation of simple interface
 * \ingroup mofem_simple_interface
 */

/* MIT License
 *
 * Copyright (c) 2022
 * 
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 * 
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 * 
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */
namespace MoFEM {

MoFEMErrorCode Simple::query_interface(boost::typeindex::type_index type_index,
                                       UnknownInterface **iface) const {
  *iface = const_cast<Simple *>(this);
  return 0;
}

template <int DIM> MoFEMErrorCode Simple::setSkeletonAdjacency() {
  static_assert(DIM == 2 || DIM == 3, "not implemented");
  return MOFEM_NOT_IMPLEMENTED;
}

template <> MoFEMErrorCode Simple::setSkeletonAdjacency<2>() {
  Interface &m_field = cOre;
  MoFEMFunctionBegin;

  auto defaultSkeletonEdge =
      [&](moab::Interface &moab, const Field &field, const EntFiniteElement &fe,
          std::vector<EntityHandle> &adjacency) -> MoFEMErrorCode {
    MoFEMFunctionBegin;

    CHKERR DefaultElementAdjacency::defaultEdge(moab, field, fe, adjacency);

    if (std::find(domainFields.begin(), domainFields.end(), field.getName()) !=
        domainFields.end()) {

      const EntityHandle fe_ent = fe.getEnt();
      std::vector<EntityHandle> bride_adjacency_edge;
      CHKERR moab.get_adjacencies(&fe_ent, 1, 2, false, bride_adjacency_edge);

      switch (field.getSpace()) {
      case H1:
        CHKERR moab.get_connectivity(&*bride_adjacency_edge.begin(),
                                     bride_adjacency_edge.size(), adjacency,
                                     true);
      case HCURL:
      case HDIV:
        CHKERR moab.get_adjacencies(&*bride_adjacency_edge.begin(),
                                    bride_adjacency_edge.size(), 1, false,
                                    adjacency, moab::Interface::UNION);
      case L2:
        adjacency.insert(adjacency.end(), bride_adjacency_edge.begin(),
                         bride_adjacency_edge.end());
        break;
      case NOFIELD:
        break;
      default:
        SETERRQ(PETSC_COMM_SELF, MOFEM_NOT_IMPLEMENTED,
                "this field is not implemented for TRI finite element");
      }

      std::sort(adjacency.begin(), adjacency.end());
      auto it = std::unique(adjacency.begin(), adjacency.end());

      std::vector<EntityHandle> new_adjacency(
          std::distance(adjacency.begin(), it));
      std::copy(adjacency.begin(), it, new_adjacency.begin());

      for (auto e : new_adjacency) {
        auto side_table = fe.getSideNumberTable();
        if (side_table.find(e) == side_table.end())
          const_cast<SideNumber_multiIndex &>(fe.getSideNumberTable())
              .insert(
                  boost::shared_ptr<SideNumber>(new SideNumber(e, -1, 0, 0)));
      }

      adjacency.swap(new_adjacency);
    }

    MoFEMFunctionReturn(0);
  };

  CHKERR m_field.modify_finite_element_adjacency_table(skeletonFE, MBEDGE,
                                                       defaultSkeletonEdge);

  MoFEMFunctionReturn(0);
}

template <> MoFEMErrorCode Simple::setSkeletonAdjacency<3>() {
  Interface &m_field = cOre;
  MoFEMFunctionBegin;

  auto defaultSkeletonEdge =
      [&](moab::Interface &moab, const Field &field, const EntFiniteElement &fe,
          std::vector<EntityHandle> &adjacency) -> MoFEMErrorCode {
    MoFEMFunctionBegin;

    CHKERR DefaultElementAdjacency::defaultFace(moab, field, fe, adjacency);

    if (std::find(domainFields.begin(), domainFields.end(), field.getName()) !=
        domainFields.end()) {

      const EntityHandle fe_ent = fe.getEnt();
      std::vector<EntityHandle> bride_adjacency_edge;
      CHKERR moab.get_adjacencies(&fe_ent, 1, 2, false, bride_adjacency_edge);

      switch (field.getSpace()) {
      case H1:
        CHKERR moab.get_connectivity(&*bride_adjacency_edge.begin(),
                                     bride_adjacency_edge.size(), adjacency,
                                     true);
      case HCURL:
        CHKERR moab.get_adjacencies(&*bride_adjacency_edge.begin(),
                                    bride_adjacency_edge.size(), 1, false,
                                    adjacency, moab::Interface::UNION);
      case HDIV:
        CHKERR moab.get_adjacencies(&*bride_adjacency_edge.begin(),
                                    bride_adjacency_edge.size(), 2, false,
                                    adjacency, moab::Interface::UNION);
      case L2:
        adjacency.insert(adjacency.end(), bride_adjacency_edge.begin(),
                         bride_adjacency_edge.end());
        break;
      case NOFIELD:
        break;
      default:
        SETERRQ(PETSC_COMM_SELF, MOFEM_NOT_IMPLEMENTED,
                "this field is not implemented for TRI finite element");
      }

      std::sort(adjacency.begin(), adjacency.end());
      auto it = std::unique(adjacency.begin(), adjacency.end());

      std::vector<EntityHandle> new_adjacency(
          std::distance(adjacency.begin(), it));
      std::copy(adjacency.begin(), it, new_adjacency.begin());

      for (auto e : new_adjacency) {
        auto side_table = fe.getSideNumberTable();
        if (side_table.find(e) == side_table.end())
          const_cast<SideNumber_multiIndex &>(fe.getSideNumberTable())
              .insert(
                  boost::shared_ptr<SideNumber>(new SideNumber(e, -1, 0, 0)));
      }

      adjacency.swap(new_adjacency);
    }

    MoFEMFunctionReturn(0);
  };

  CHKERR m_field.modify_finite_element_adjacency_table(skeletonFE, MBTRI,
                                                       defaultSkeletonEdge);
  CHKERR m_field.modify_finite_element_adjacency_table(skeletonFE, MBQUAD,
                                                       defaultSkeletonEdge);

  MoFEMFunctionReturn(0);
}

template <> MoFEMErrorCode Simple::setSkeletonAdjacency<-1>() {
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

template <int DIM> MoFEMErrorCode Simple::setParentAdjacency() {
  static_assert(DIM == 2 || DIM == 3, "not implemented");
  return MOFEM_NOT_IMPLEMENTED;
}

template <> MoFEMErrorCode Simple::setParentAdjacency<3>() {
  Interface &m_field = cOre;
  MoFEMFunctionBegin;

  parentAdjFunctionDim3 =
      boost::make_shared<ParentFiniteElementAdjacencyFunction<3>>(
          bitAdjParent, bitAdjParentMask, bitAdjEnt, bitAdjEntMask);
  parentAdjFunctionDim2 =
      boost::make_shared<ParentFiniteElementAdjacencyFunction<2>>(
          bitAdjParent, bitAdjParentMask, bitAdjEnt, bitAdjEntMask);

  CHKERR m_field.modify_finite_element_adjacency_table(domainFE, MBTET,
                                                       *parentAdjFunctionDim3);
  CHKERR m_field.modify_finite_element_adjacency_table(domainFE, MBHEX,
                                                       *parentAdjFunctionDim3);
  if (addBoundaryFE || !boundaryFields.empty()) {
    CHKERR m_field.modify_finite_element_adjacency_table(
        boundaryFE, MBTRI, *parentAdjFunctionDim2);
    CHKERR m_field.modify_finite_element_adjacency_table(
        boundaryFE, MBQUAD, *parentAdjFunctionDim2);
  }
  if (addSkeletonFE || !skeletonFields.empty()) {
    CHKERR m_field.modify_finite_element_adjacency_table(
        skeletonFE, MBTRI, *parentAdjFunctionDim2);
    CHKERR m_field.modify_finite_element_adjacency_table(
        skeletonFE, MBQUAD, *parentAdjFunctionDim2);
  }

  MoFEMFunctionReturn(0);
}

template <> MoFEMErrorCode Simple::setParentAdjacency<2>() {
  Interface &m_field = cOre;
  MoFEMFunctionBegin;

  parentAdjFunctionDim2 =
      boost::make_shared<ParentFiniteElementAdjacencyFunction<2>>(
          bitAdjParent, bitAdjParentMask, bitAdjEnt, bitAdjEntMask);
  parentAdjFunctionDim1 =
      boost::make_shared<ParentFiniteElementAdjacencyFunction<1>>(
          bitAdjParent, bitAdjParentMask, bitAdjEnt, bitAdjEntMask);

  CHKERR m_field.modify_finite_element_adjacency_table(domainFE, MBTRI,
                                                       *parentAdjFunctionDim2);
  CHKERR m_field.modify_finite_element_adjacency_table(domainFE, MBQUAD,
                                                       *parentAdjFunctionDim2);
  if (addBoundaryFE || !boundaryFields.empty())
    CHKERR m_field.modify_finite_element_adjacency_table(
        boundaryFE, MBEDGE, *parentAdjFunctionDim1);
  if (addSkeletonFE || !skeletonFields.empty())
    CHKERR m_field.modify_finite_element_adjacency_table(
        skeletonFE, MBEDGE, *parentAdjFunctionDim1);

  MoFEMFunctionReturn(0);
}

template <> MoFEMErrorCode Simple::setParentAdjacency<-1>() {
  MoFEMFunctionBegin;
  switch (getDim()) {
  case 1:
    THROW_MESSAGE("Not implemented");
  case 2:
    return setParentAdjacency<2>();
  case 3:
    return setParentAdjacency<3>();
  default:
    THROW_MESSAGE("Not implemented");
  }
  MoFEMFunctionReturn(0);
}

MoFEMErrorCode Simple::setSkeletonAdjacency(int dim) {
  MoFEMFunctionBegin;
  if (dim == -1)
    dim = getDim();
  switch (dim) {
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
      bitLevelMask(BitRefLevel().set()), meshSet(0), boundaryMeshset(0),
      skeletonMeshset(0), nameOfProblem("SimpleProblem"), domainFE("dFE"),
      boundaryFE("bFE"), skeletonFE("sFE"), dIm(-1), addSkeletonFE(false),
      addBoundaryFE(false), addParentAdjacencies(false),
      bitAdjParent(BitRefLevel().set()), bitAdjParentMask(BitRefLevel().set()),
      bitAdjEnt(BitRefLevel().set()), bitAdjEntMask(BitRefLevel().set()) {
  PetscLogEventRegister("SimpleSetUp", 0, &MOFEM_EVENT_SimpleSetUP);
  PetscLogEventRegister("SimpleLoadMesh", 0, &MOFEM_EVENT_SimpleLoadMesh);
  PetscLogEventRegister("SimpleBuildFields", 0, &MOFEM_EVENT_SimpleBuildFields);
  PetscLogEventRegister("SimpleBuildFiniteElements", 0,
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
  CHKERR m_field.get_moab().load_file(meshFileName, 0, options.c_str());
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

  if (!boundaryMeshset)
    CHKERR createBoundaryMeshset();
  if (!skeletonMeshset)
    CHKERR createSkeletonMeshset();
  if (addSkeletonFE)
    CHKERR exchangeGhostCells();

  if (bitLevel.any()) {
    Range ents;
    CHKERR m_field.get_moab().get_entities_by_dimension(meshSet, dIm, ents,
                                                        true);
    CHKERR m_field.getInterface<BitRefManager>()->setBitRefLevel(ents, bitLevel,
                                                                 false);
  } else {
    MOFEM_LOG("WORLD", Sev::warning) << "BitRefLevel is none and not set";
    CHKERR m_field.getInterface<BitRefManager>()->addToDatabaseBitRefLevelByDim(
        dIm, BitRefLevel().set(), BitRefLevel().set());
  }

  PetscLogEventEnd(MOFEM_EVENT_SimpleLoadMesh, 0, 0, 0, 0);
  MoFEMFunctionReturn(0);
}

MoFEMErrorCode Simple::loadFile(const std::string mesh_file_name) {
  MoFEMFunctionBegin;
  Interface &m_field = cOre;
  if (m_field.get_comm_size() == 1)
    CHKERR loadFile("", meshFileName);
  else
    CHKERR loadFile("PARALLEL=READ_PART;"
                    "PARALLEL_RESOLVE_SHARED_ENTS;"
                    "PARTITION=PARALLEL_PARTITION;",
                    meshFileName);
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
    auto fe_ptr = m_field.get_finite_elements();
    auto &fe_by_name = const_cast<FiniteElement_multiIndex *>(fe_ptr)
                           ->get<FiniteElement_name_mi_tag>();
    auto it_fe = fe_by_name.find(fe_name);
    if (it_fe != fe_by_name.end()) {

      if (!fe_by_name.modify(it_fe, FiniteElement_row_change_bit_reset()))
        SETERRQ(m_field.get_comm(), MOFEM_OPERATION_UNSUCCESSFUL,
                "modification unsuccessful");

      if (!fe_by_name.modify(it_fe, FiniteElement_col_change_bit_reset()))
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

  if (addBoundaryFE || !boundaryFields.empty()) {
    CHKERR m_field.add_finite_element(boundaryFE, MF_ZERO);
    CHKERR add_fields(boundaryFE, domainFields);
    if (!boundaryFields.empty())
      CHKERR add_fields(boundaryFE, boundaryFields);
    CHKERR add_data_fields(boundaryFE, dataFields);
    CHKERR add_data_fields(boundaryFE, noFieldDataFields);
    CHKERR add_fields(boundaryFE, noFieldFields);
  }
  if (addSkeletonFE || !skeletonFields.empty()) {
    CHKERR m_field.add_finite_element(skeletonFE, MF_ZERO);
    CHKERR add_fields(skeletonFE, domainFields);
    if (!skeletonFields.empty())
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
  CHKERR DMMoFEMCreateMoFEM(dM, &m_field, nameOfProblem.c_str(), bitLevel,
                            bitLevelMask);
  CHKERR DMSetFromOptions(dM);
  CHKERR DMMoFEMAddElement(dM, domainFE.c_str());
  if (addBoundaryFE || !boundaryFields.empty()) {
    CHKERR DMMoFEMAddElement(dM, boundaryFE.c_str());
  }
  if (addSkeletonFE || !skeletonFields.empty()) {
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
      CHKERR m_field.create_vertices_and_add_to_field(field, coords.data(), 1);
      CHKERR comm_interface_ptr->makeFieldEntitiesMultishared(field, 0);
      CHKERR bit_ref_ptr->setFieldEntitiesBitRefLevel(field, bitLevel);
    }
    MoFEMFunctionReturnHot(0);
  };

  CHKERR add_ents_to_field(meshSet, dIm, domainFields);
  CHKERR add_ents_to_field(meshSet, dIm, dataFields);
  CHKERR add_ents_to_field(boundaryMeshset, dIm - 1, boundaryFields);
  CHKERR add_ents_to_field(skeletonMeshset, dIm - 1, skeletonFields);

  std::set<std::string> nofield_fields;
  for (auto &f : noFieldFields)
    nofield_fields.insert(f);
  for (auto &f : noFieldDataFields)
    nofield_fields.insert(f);

  CHKERR make_no_field_ents(nofield_fields);

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
  if (addBoundaryFE || !boundaryFields.empty()) {
    CHKERR m_field.add_ents_to_finite_element_by_dim(boundaryMeshset, dIm - 1,
                                                     boundaryFE, true);
    CHKERR m_field.build_finite_elements(boundaryFE);
  }
  if (addSkeletonFE || !skeletonFields.empty()) {
    CHKERR m_field.add_ents_to_finite_element_by_dim(skeletonMeshset, dIm - 1,
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
  CHKERR m_field.build_adjacencies(bitLevel, bitLevelMask);
  // Set problem by the DOFs on the fields rather that by adding DOFs on the
  // elements
  m_field.getInterface<ProblemsManager>()->buildProblemFromFields = PETSC_TRUE;
  CHKERR DMSetUp_MoFEM(dM);
  m_field.getInterface<ProblemsManager>()->buildProblemFromFields = PETSC_FALSE;
  PetscLogEventEnd(MOFEM_EVENT_SimpleBuildProblem, 0, 0, 0, 0);
  MoFEMFunctionReturn(0);
}

MoFEMErrorCode Simple::setUp(const PetscBool is_partitioned) {
  Interface &m_field = cOre;
  MoFEMFunctionBegin;

  PetscLogEventBegin(MOFEM_EVENT_SimpleSetUP, 0, 0, 0, 0);

  CHKERR defineFiniteElements();

  if (addSkeletonFE || !skeletonFields.empty())
    CHKERR setSkeletonAdjacency();

  if (addParentAdjacencies)
    CHKERR setParentAdjacency();

  CHKERR defineProblem(is_partitioned);
  CHKERR buildFields();
  CHKERR buildFiniteElements();
  CHKERR buildProblem();

  PetscLogEventEnd(MOFEM_EVENT_SimpleSetUP, 0, 0, 0, 0);
  MoFEMFunctionReturn(0);
}

MoFEMErrorCode Simple::reSetUp(bool only_dm) {
  Interface &m_field = cOre;
  MoFEMFunctionBegin;
  PetscLogEventBegin(MOFEM_EVENT_SimpleBuildProblem, 0, 0, 0, 0);

  if (!only_dm) {
    CHKERR defineFiniteElements();
    CHKERR buildFields();
    if (addSkeletonFE || !skeletonFields.empty())
      CHKERR setSkeletonAdjacency();
    if (addParentAdjacencies)
      CHKERR setParentAdjacency();
    CHKERR buildFiniteElements();
  }

  CHKERR m_field.build_adjacencies(bitLevel, bitLevelMask);

  const Problem *problem_ptr;
  CHKERR DMMoFEMGetProblemPtr(dM, &problem_ptr);
  const auto problem_name = problem_ptr->getName();
  CHKERR m_field.modify_problem_ref_level_set_bit(problem_name, bitLevel);
  CHKERR m_field.modify_problem_mask_ref_level_set_bit(problem_name,
                                                       bitLevelMask);

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

/**
 * @brief Delete dm
 *
 * @return MoFEMErrorCode
 */
MoFEMErrorCode Simple::deleteDM() {
  MoFEMFunctionBegin;
  dM.reset();
  MoFEMFunctionReturn(0);
}

/**
 * @brief Delete finite elements
 *
 * @return MoFEMErrorCode
 */
MoFEMErrorCode Simple::deleteFiniteElements() {
  Interface &m_field = cOre;
  MoFEMFunctionBegin;
  for (auto fe : {domainFE, boundaryFE, skeletonFE}) {
    if (m_field.check_finite_element(fe)) {
      CHKERR m_field.delete_finite_element(fe);
    }
  }
  MoFEMFunctionReturn(0);
}

MoFEMErrorCode
Simple::addFieldToEmptyFieldBlocks(const std::string row_field,
                                   const std::string col_field) const {
  Interface &m_field = cOre;
  MoFEMFunctionBegin;
  CHKERR m_field.getInterface<ProblemsManager>()->addFieldToEmptyFieldBlocks(
      getProblemName(), row_field, col_field);
  MoFEMFunctionReturn(0);
}

MoFEMErrorCode Simple::createBoundaryMeshset() {
  Interface &m_field = cOre;
  MoFEMFunctionBegin;
  ParallelComm *pcomm =
      ParallelComm::get_pcomm(&m_field.get_moab(), MYPCOMM_INDEX);

  auto get_skin = [&](auto meshset) {
    // filter not owned entities, those are not on boundary

    Range domain_ents;
    CHKERR m_field.get_moab().get_entities_by_dimension(meshset, dIm,
                                                        domain_ents, true);
    CHKERR pcomm->filter_pstatus(domain_ents,
                                 PSTATUS_SHARED | PSTATUS_MULTISHARED,
                                 PSTATUS_NOT, -1, nullptr);

    Skinner skin(&m_field.get_moab());
    Range domain_skin;
    CHKERR skin.find_skin(0, domain_ents, false, domain_skin);
    CHKERR pcomm->filter_pstatus(domain_skin,
                                 PSTATUS_SHARED | PSTATUS_MULTISHARED,
                                 PSTATUS_NOT, -1, nullptr);
    return domain_skin;
  };

  auto create_boundary_meshset = [&](auto &&domain_skin) {
    MoFEMFunctionBeginHot;
    // create boundary meshset
    if (boundaryMeshset != 0) {
      MoFEMFunctionReturnHot(0);
    }
    CHKERR m_field.get_moab().create_meshset(MESHSET_SET, boundaryMeshset);
    CHKERR m_field.get_moab().add_entities(boundaryMeshset, domain_skin);
    for (int dd = 0; dd != dIm - 1; dd++) {
      Range adj;
      CHKERR m_field.get_moab().get_adjacencies(domain_skin, dd, false, adj,
                                                moab::Interface::UNION);
      CHKERR m_field.get_moab().add_entities(boundaryMeshset, adj);
    }
    MoFEMFunctionReturnHot(0);
  };

  CHKERR create_boundary_meshset(get_skin(meshSet));

  MoFEMFunctionReturn(0);
}

MoFEMErrorCode Simple::createSkeletonMeshset() {
  Interface &m_field = cOre;
  MoFEMFunctionBegin;

  ParallelComm *pcomm =
      ParallelComm::get_pcomm(&m_field.get_moab(), MYPCOMM_INDEX);

  auto create_skeleton_meshset = [&](auto meshset) {
    MoFEMFunctionBeginHot;
    // create boundary meshset
    if (skeletonMeshset != 0) {
      MoFEMFunctionReturnHot(0);
    }
    Range boundary_ents, skeleton_ents;
    CHKERR m_field.get_moab().get_entities_by_dimension(boundaryMeshset,
                                                        dIm - 1, boundary_ents);
    Range domain_ents;
    CHKERR m_field.get_moab().get_entities_by_dimension(meshset, dIm,
                                                        domain_ents, true);
    CHKERR m_field.get_moab().get_adjacencies(
        domain_ents, dIm - 1, false, skeleton_ents, moab::Interface::UNION);
    skeleton_ents = subtract(skeleton_ents, boundary_ents);
    CHKERR pcomm->filter_pstatus(skeleton_ents, PSTATUS_NOT_OWNED, PSTATUS_NOT,
                                 -1, nullptr);
    CHKERR m_field.get_moab().create_meshset(MESHSET_SET, skeletonMeshset);
    CHKERR m_field.get_moab().add_entities(skeletonMeshset, skeleton_ents);
    MoFEMFunctionReturnHot(0);
  };

  CHKERR create_skeleton_meshset(meshSet);

  MoFEMFunctionReturn(0);
}

MoFEMErrorCode Simple::exchangeGhostCells() {
  Interface &m_field = cOre;
  MoFEMFunctionBegin;
  MOFEM_LOG("WORLD", Sev::verbose) << "Exchange ghost cells";

  ParallelComm *pcomm =
      ParallelComm::get_pcomm(&m_field.get_moab(), MYPCOMM_INDEX);
  if (pcomm == NULL)
    pcomm = new ParallelComm(&m_field.get_moab(), m_field.get_comm());

  Range verts;
  CHKERR m_field.get_moab().get_entities_by_type(0, MBVERTEX, verts);
  CHKERR pcomm->exchange_ghost_cells(getDim(), getDim() - 1, 1,
                                     3 /**get all adjacent ghosted entities */,
                                     true, false, meshSet ? &meshSet : nullptr);

  Range shared;
  CHKERR m_field.get_moab().get_entities_by_dimension(0, dIm, shared);
  for (auto d = dIm - 1; d >= 1; --d) {
    CHKERR m_field.get_moab().get_adjacencies(shared, d, false, shared,
                                              moab::Interface::UNION);
  }
  CHKERR pcomm->filter_pstatus(shared, PSTATUS_SHARED | PSTATUS_MULTISHARED,
                               PSTATUS_OR, -1, &shared);
  Tag part_tag = pcomm->part_tag();
  CHKERR pcomm->exchange_tags(part_tag, shared);

  MoFEMFunctionReturn(0);
}

} // namespace MoFEM
