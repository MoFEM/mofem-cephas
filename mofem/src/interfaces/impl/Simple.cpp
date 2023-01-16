/** \file Simple.cpp
 * \brief Implementation of simple interface
 * \ingroup mofem_simple_interface
 */

namespace MoFEM {

MoFEMErrorCode Simple::query_interface(boost::typeindex::type_index type_index,
                                       UnknownInterface **iface) const {
  *iface = const_cast<Simple *>(this);
  return 0;
}

template <int DIM>
MoFEMErrorCode Simple::setSkeletonAdjacency(std::string fe_name) {
  static_assert(DIM == 2 || DIM == 3, "not implemented");
  return MOFEM_NOT_IMPLEMENTED;
}

template <>
MoFEMErrorCode Simple::setSkeletonAdjacency<2>(std::string fe_name) {
  Interface &m_field = cOre;
  MoFEMFunctionBegin;

  parentAdjSkeletonFunctionDim1 =
      boost::make_shared<ParentFiniteElementAdjacencyFunctionSkeleton<1>>(
          bitAdjParent, bitAdjParentMask, bitAdjEnt, bitAdjEntMask);
  CHKERR m_field.modify_finite_element_adjacency_table(
      fe_name, MBEDGE, *parentAdjSkeletonFunctionDim1);

  MoFEMFunctionReturn(0);
}

template <>
MoFEMErrorCode Simple::setSkeletonAdjacency<3>(std::string fe_name) {
  Interface &m_field = cOre;
  MoFEMFunctionBegin;

  parentAdjSkeletonFunctionDim2 =
      boost::make_shared<ParentFiniteElementAdjacencyFunctionSkeleton<2>>(
          bitAdjParent, bitAdjParentMask, bitAdjEnt, bitAdjEntMask);

  CHKERR m_field.modify_finite_element_adjacency_table(
      fe_name, MBTRI, *parentAdjSkeletonFunctionDim2);
  CHKERR m_field.modify_finite_element_adjacency_table(
      fe_name, MBQUAD, *parentAdjSkeletonFunctionDim2);

  MoFEMFunctionReturn(0);
}

template <>
MoFEMErrorCode Simple::setSkeletonAdjacency<-1>(std::string fe_name) {
  MoFEMFunctionBegin;

  switch (getDim()) {
  case 1:
    THROW_MESSAGE("Not implemented");
  case 2:
    return setSkeletonAdjacency<2>(fe_name);
  case 3:
    return setSkeletonAdjacency<3>(fe_name);
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

MoFEMErrorCode Simple::setSkeletonAdjacency(int dim, std::string fe_name) {
  MoFEMFunctionBegin;
  if (dim == -1)
    dim = getDim();

  if (fe_name.empty())
    fe_name = skeletonFE;

  switch (dim) {
  case 2:
    return setSkeletonAdjacency<2>(fe_name);
  case 3:
    return setSkeletonAdjacency<3>(fe_name);
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
  CHKERR DMMoFEMAddElement(dM, domainFE);
  if (addBoundaryFE || !boundaryFields.empty()) {
    CHKERR DMMoFEMAddElement(dM, boundaryFE);
  }
  if (addSkeletonFE || !skeletonFields.empty()) {
    CHKERR DMMoFEMAddElement(dM, skeletonFE);
  }
  CHKERR DMMoFEMAddElement(dM, otherFEs);
  CHKERR DMMoFEMSetIsPartitioned(dM, is_partitioned);
  MoFEMFunctionReturn(0);
}

MoFEMErrorCode Simple::setFieldOrder(const std::string field_name,
                                     const int order, const Range *ents) {
  MoFEMFunctionBeginHot;
  fieldsOrder.emplace_back(field_name, order,
                           ents == NULL ? Range() : Range(*ents));
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

  auto get_field_ptr = [&](auto &f) { return m_field.get_field_structure(f); };
  for (auto &t : fieldsOrder) {
    const auto f = std::get<0>(t);
    const auto order = std::get<1>(t);

    MOFEM_TAG_AND_LOG("WORLD", Sev::inform, "Simple")
        << "Set order to field " << f << " order " << order;
    if (!std::get<2>(t).empty()) {
      MOFEM_LOG_CHANNEL("SYNC");
      MOFEM_TAG_AND_LOG("SYNC", Sev::verbose, "Simple")
          << "To ents: " << std::endl
          << std::get<2>(t) << std::endl;
      MOFEM_LOG_SYNCHRONISE(m_field.get_comm());
    }

    if (std::get<2>(t).empty()) {
      auto f_ptr = get_field_ptr(f);

      if (f_ptr->getSpace() == H1) {
        if (f_ptr->getApproxBase() == AINSWORTH_BERNSTEIN_BEZIER_BASE) {
          CHKERR m_field.set_field_order(meshSet, MBVERTEX, f, order);
        } else {
          CHKERR m_field.set_field_order(meshSet, MBVERTEX, f, 1);
        }
      }

      for (auto d = 1; d <= dIm; ++d) {
        for (EntityType t = CN::TypeDimensionMap[d].first;
             t <= CN::TypeDimensionMap[d].second; ++t) {
          CHKERR m_field.set_field_order(meshSet, t, f, order);
        }
      }
    } else {
      CHKERR m_field.set_field_order(std::get<2>(t), f, order);
    }
  }
  MOFEM_LOG_CHANNEL("WORLD");
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
