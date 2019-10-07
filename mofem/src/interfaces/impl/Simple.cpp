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
Simple::~Simple() {}

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

MoFEMErrorCode Simple::loadFile(const std::string options) {
  Interface &m_field = cOre;
  MoFEMFunctionBegin;
  PetscLogEventBegin(MOFEM_EVENT_SimpleLoadMesh, 0, 0, 0, 0);
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
  dataFields.push_back(name);
  MoFEMFunctionReturn(0);
}

MoFEMErrorCode Simple::defineFiniteElements() {
  Interface &m_field = cOre;
  MoFEMFunctionBeginHot;
  // Define finite elements
  CHKERR m_field.add_finite_element(domainFE);
  for (unsigned int ff = 0; ff != domainFields.size(); ff++) {
    CHKERR m_field.modify_finite_element_add_field_row(domainFE,
                                                       domainFields[ff]);
    CHKERR m_field.modify_finite_element_add_field_col(domainFE,
                                                       domainFields[ff]);
    CHKERR m_field.modify_finite_element_add_field_data(domainFE,
                                                        domainFields[ff]);
  }
  for (unsigned int ff = 0; ff != dataFields.size(); ff++) {
    CHKERR m_field.modify_finite_element_add_field_data(domainFE,
                                                        dataFields[ff]);
  }
  if (!boundaryFields.empty()) {
    CHKERR m_field.add_finite_element(boundaryFE);
    for (unsigned int ff = 0; ff != domainFields.size(); ff++) {
      CHKERR m_field.modify_finite_element_add_field_row(boundaryFE,
                                                         domainFields[ff]);
      CHKERR m_field.modify_finite_element_add_field_col(boundaryFE,
                                                         domainFields[ff]);
      CHKERR m_field.modify_finite_element_add_field_data(boundaryFE,
                                                          domainFields[ff]);
    }
    for (unsigned int ff = 0; ff != boundaryFields.size(); ff++) {
      CHKERR m_field.modify_finite_element_add_field_row(boundaryFE,
                                                         boundaryFields[ff]);
      CHKERR m_field.modify_finite_element_add_field_col(boundaryFE,
                                                         boundaryFields[ff]);
      CHKERR m_field.modify_finite_element_add_field_data(boundaryFE,
                                                          boundaryFields[ff]);
    }
    for (unsigned int ff = 0; ff != skeletonFields.size(); ff++) {
      CHKERR m_field.modify_finite_element_add_field_row(boundaryFE,
                                                         skeletonFields[ff]);
      CHKERR m_field.modify_finite_element_add_field_col(boundaryFE,
                                                         skeletonFields[ff]);
      CHKERR m_field.modify_finite_element_add_field_data(boundaryFE,
                                                          skeletonFields[ff]);
    }
  }
  if (!skeletonFields.empty()) {
    CHKERR m_field.add_finite_element(skeletonFE);
    for (unsigned int ff = 0; ff != domainFields.size(); ff++) {
      CHKERR m_field.modify_finite_element_add_field_row(skeletonFE,
                                                         domainFields[ff]);
      CHKERR m_field.modify_finite_element_add_field_col(skeletonFE,
                                                         domainFields[ff]);
      CHKERR m_field.modify_finite_element_add_field_data(skeletonFE,
                                                          domainFields[ff]);
    }
    for (unsigned int ff = 0; ff != boundaryFields.size(); ff++) {
      CHKERR m_field.modify_finite_element_add_field_row(skeletonFE,
                                                         boundaryFields[ff]);
      CHKERR m_field.modify_finite_element_add_field_col(skeletonFE,
                                                         boundaryFields[ff]);
      CHKERR m_field.modify_finite_element_add_field_data(skeletonFE,
                                                          boundaryFields[ff]);
    }
    for (unsigned int ff = 0; ff != skeletonFields.size(); ff++) {
      CHKERR m_field.modify_finite_element_add_field_row(skeletonFE,
                                                         skeletonFields[ff]);
      CHKERR m_field.modify_finite_element_add_field_col(skeletonFE,
                                                         skeletonFields[ff]);
      CHKERR m_field.modify_finite_element_add_field_data(skeletonFE,
                                                          skeletonFields[ff]);
    }
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
  fieldsOrder[field_name] =
      std::pair<int, Range>(order, ents == NULL ? Range() : Range(*ents));
  MoFEMFunctionReturnHot(0);
}

MoFEMErrorCode Simple::buildFields() {
  Interface &m_field = cOre;
  MoFEMFunctionBegin;
  PetscLogEventBegin(MOFEM_EVENT_SimpleBuildFields, 0, 0, 0, 0);
  // take skin
  {
    Range domain_ents;
    CHKERR m_field.get_moab().get_entities_by_dimension(meshSet, dIm,
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
  }
  auto comm_interface_ptr = m_field.getInterface<CommInterface>();
  // Add entities to the fields
  for (unsigned int ff = 0; ff != domainFields.size(); ff++) {
    CHKERR m_field.add_ents_to_field_by_dim(meshSet, dIm, domainFields[ff]);
    CHKERR comm_interface_ptr->synchroniseFieldEntities(domainFields[ff], 0);
  }
  for (unsigned int ff = 0; ff != dataFields.size(); ff++) {
    CHKERR m_field.add_ents_to_field_by_dim(meshSet, dIm, dataFields[ff]);
    CHKERR comm_interface_ptr->synchroniseFieldEntities(dataFields[ff], 0);
  }
  for (unsigned int ff = 0; ff != boundaryFields.size(); ff++) {
    CHKERR m_field.add_ents_to_field_by_dim(boundaryMeshset, dIm - 1,
                                            boundaryFields[ff]);
    CHKERR comm_interface_ptr->synchroniseFieldEntities(boundaryFields[ff], 0);
  }
  for (unsigned int ff = 0; ff != skeletonFields.size(); ff++) {
    CHKERR m_field.add_ents_to_field_by_dim(meshSet, dIm - 1,
                                            skeletonFields[ff]);
    CHKERR comm_interface_ptr->synchroniseFieldEntities(skeletonFields[ff], 0);
  }
  // Set order
  for (unsigned int ff = 0; ff != domainFields.size(); ff++) {
    if (fieldsOrder.find(domainFields[ff]) == fieldsOrder.end()) {
      SETERRQ1(PETSC_COMM_WORLD, MOFEM_INVALID_DATA,
               "Order for field not set %s", domainFields[ff].c_str());
    }
    int dds = 0;
    const Field *field = m_field.get_field_structure(domainFields[ff]);
    switch (field->getSpace()) {
    case L2:
      dds = dIm;
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
              "Huston we have a problem");
    }
    if (field->getSpace() == H1) {
      CHKERR m_field.set_field_order(meshSet, MBVERTEX, domainFields[ff], 1);
    }
    for (int dd = dds; dd <= dIm; dd++) {
      Range ents;
      CHKERR m_field.get_field_entities_by_dimension(domainFields[ff], dd,
                                                     ents);
      if (!fieldsOrder.at(domainFields[ff]).second.empty()) {
        ents = intersect(ents, fieldsOrder.at(domainFields[ff]).second);
      }
      CHKERR m_field.set_field_order(ents, domainFields[ff],
                                     fieldsOrder.at(domainFields[ff]).first);
    }
  }
  // Set order to data fields
  for (unsigned int ff = 0; ff != dataFields.size(); ff++) {
    if (fieldsOrder.find(dataFields[ff]) == fieldsOrder.end()) {
      SETERRQ1(PETSC_COMM_WORLD, MOFEM_INVALID_DATA,
               "Order for field not set %s", dataFields[ff].c_str());
    }
    int dds = 0;
    const Field *field = m_field.get_field_structure(dataFields[ff]);
    switch (field->getSpace()) {
    case L2:
      dds = dIm;
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
              "Huston we have a problem");
    }
    if (field->getSpace() == H1) {
      CHKERR m_field.set_field_order(meshSet, MBVERTEX, dataFields[ff], 1);
    }
    for (int dd = dds; dd <= dIm; dd++) {
      Range ents;
      CHKERR m_field.get_field_entities_by_dimension(dataFields[ff], dd, ents);
      CHKERR m_field.set_field_order(ents, dataFields[ff],
                                     fieldsOrder.at(dataFields[ff]).first);
    }
  }
  // Set order to boundary
  for (unsigned int ff = 0; ff != boundaryFields.size(); ff++) {
    if (fieldsOrder.find(boundaryFields[ff]) == fieldsOrder.end()) {
      SETERRQ1(PETSC_COMM_WORLD, MOFEM_INVALID_DATA,
               "Order for field not set %s", boundaryFields[ff].c_str());
    }
    int dds = 0;
    const Field *field = m_field.get_field_structure(boundaryFields[ff]);
    switch (field->getSpace()) {
    case L2:
      dds = dIm - 1;
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
              "Huston we have a problem");
    }
    if (field->getSpace() == H1) {
      CHKERR m_field.set_field_order(meshSet, MBVERTEX, boundaryFields[ff], 1);
    }
    for (int dd = dds; dd <= dIm - 1; dd++) {
      Range ents;
      CHKERR m_field.get_field_entities_by_dimension(boundaryFields[ff], dd,
                                                     ents);
      CHKERR m_field.set_field_order(ents, boundaryFields[ff],
                                     fieldsOrder.at(boundaryFields[ff]).first);
    }
  }
  // Set order to skeleton
  for (unsigned int ff = 0; ff != skeletonFields.size(); ff++) {
    if (fieldsOrder.find(skeletonFields[ff]) == fieldsOrder.end()) {
      SETERRQ1(PETSC_COMM_WORLD, MOFEM_INVALID_DATA,
               "Order for field not set %s", skeletonFields[ff].c_str());
    }
    int dds = 0;
    const Field *field = m_field.get_field_structure(skeletonFields[ff]);
    switch (field->getSpace()) {
    case L2:
      dds = dIm - 1;
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
              "Huston we have a problem");
    }
    if (field->getSpace() == H1) {
      CHKERR m_field.set_field_order(meshSet, MBVERTEX, skeletonFields[ff], 1);
    }
    for (int dd = dds; dd <= dIm - 1; dd++) {
      Range ents;
      CHKERR m_field.get_field_entities_by_dimension(skeletonFields[ff], dd,
                                                     ents);
      CHKERR m_field.set_field_order(ents, skeletonFields[ff],
                                     fieldsOrder.at(skeletonFields[ff]).first);
    }
  }
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
  PetscLogEventEnd(MOFEM_EVENT_SimpleBuildProblem, 0, 0, 0, 0);
  MoFEMFunctionReturn(0);
}

MoFEMErrorCode Simple::setUp(const PetscBool is_partitioned) {
  MoFEMFunctionBegin;
  CHKERR defineFiniteElements();
  CHKERR defineProblem(is_partitioned);
  CHKERR buildFields();
  CHKERR buildFiniteElements();
  CHKERR buildProblem();
  MoFEMFunctionReturn(0);
}

MoFEMErrorCode Simple::getDM(DM *dm) {
  MoFEMFunctionBegin;
  CHKERR PetscObjectReference(getPetscObject(dM.get()));
  *dm = dM.get();
  MoFEMFunctionReturn(0);
}

} // namespace MoFEM
