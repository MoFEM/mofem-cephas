/** \file Projection10NodeCoordsOnField.cpp
 * \brief Project displacements/coordinates from 10 node tetrahedra on
 * hierarchical approximation base.
 */



namespace MoFEM {

Projection10NodeCoordsOnField::Projection10NodeCoordsOnField(
    Interface &m_field, std::string field_name, int verb)
    : mField(m_field), fieldName(field_name), vErbose(verb) {}

MoFEMErrorCode Projection10NodeCoordsOnField::preProcess() {
  MoFEMFunctionBeginHot;
  auto field_ptr = mField.get_field_structure(fieldName);
  if (field_ptr->getApproxBase() == AINSWORTH_BERNSTEIN_BEZIER_BASE) {
    MOFEM_TAG_AND_LOG("WORLD", Sev::warning, "Projection10NodeCoordsOnField")
        << "Only working well for first order AINSWORTH_BERNSTEIN_BEZIER_BASE!";
    MOFEM_LOG_CHANNEL("WORLD");
  }
  MoFEMFunctionReturnHot(0);
}

MoFEMErrorCode Projection10NodeCoordsOnField::operator()() {
  MoFEMFunctionBeginHot;

  if (dofPtr == NULL) {
    SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY, "data inconsistency");
  }
  if (dofPtr->getName() != fieldName)
    MoFEMFunctionReturnHot(0);
  if (dofPtr->getEntType() == MBVERTEX) {
    EntityHandle node = dofPtr->getEnt();
    coords.resize(3);
    CHKERR mField.get_moab().get_coords(&node, 1, &*coords.data().begin());
    dofPtr->getFieldData() = coords[dofPtr->getDofCoeffIdx()];
    if (vErbose > 0) {
      MOFEM_TAG_AND_LOG_C("SELF", Sev::noisy, "Projection10NodeCoordsOnField",
                          "val = %6.7e\n", dofPtr->getFieldData());
      MOFEM_LOG_CHANNEL("SELF");
    }

    MoFEMFunctionReturnHot(0);
  }
  if (dofPtr->getEntType() != MBEDGE) {
    MoFEMFunctionReturnHot(0);
  }
  if (dofPtr->getEntDofIdx() != dofPtr->getDofCoeffIdx()) {
    MoFEMFunctionReturnHot(0);
  }
  EntityHandle edge = dofPtr->getEnt();
  if (type_from_handle(edge) != MBEDGE) {
    SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
            "this method works only elements which are type of MBEDGE");
  }
  // coords
  int num_nodes;
  const EntityHandle *conn;
  CHKERR mField.get_moab().get_connectivity(edge, conn, num_nodes, false);
  if ((num_nodes != 2) && (num_nodes != 3)) {
    SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
            "this method works only 4 node and 10 node tets");
  }
  coords.resize(num_nodes * 3);
  CHKERR mField.get_moab().get_coords(conn, num_nodes, &*coords.data().begin());
  aveMidCoord.resize(3);
  midNodeCoord.resize(3);
  for (int dd = 0; dd < 3; dd++) {
    aveMidCoord[dd] = (coords[0 * 3 + dd] + coords[1 * 3 + dd]) * 0.5;
    if (num_nodes == 3) {
      midNodeCoord[dd] = coords[2 * 3 + dd];
    } else {
      midNodeCoord[dd] = aveMidCoord[dd];
    }
  }

  const auto base = dofPtr->getApproxBase();
  double edge_shape_function_val;
  switch (base) {
  case AINSWORTH_LEGENDRE_BASE:
    edge_shape_function_val = 0.25;
    break;
  case AINSWORTH_LOBATTO_BASE:
    edge_shape_function_val = 0.25 * LOBATTO_PHI0(0);
    break;
  case DEMKOWICZ_JACOBI_BASE: {
    double L[3];
    CHKERR Legendre_polynomials(2, 0, NULL, L, NULL, 1);
    edge_shape_function_val = 0.125 * LOBATTO_PHI0(0);
  }; break;
  default:
    SETERRQ(PETSC_COMM_SELF, MOFEM_NOT_IMPLEMENTED, "not yet implemented");
  }

  diffNodeCoord.resize(3);
  ublas::noalias(diffNodeCoord) = midNodeCoord - aveMidCoord;
  dOf.resize(3);
  ublas::noalias(dOf) = diffNodeCoord / edge_shape_function_val;
  if (dofPtr->getNbOfCoeffs() > 3) {
    SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
            "this method works only fields which are rank 3 or lower");
  }
  dofPtr->getFieldData() = dOf[dofPtr->getDofCoeffIdx()];

  MoFEMFunctionReturnHot(0);
}

MoFEMErrorCode Projection10NodeCoordsOnField::postProcess() {
  MoFEMFunctionBeginHot;
  MoFEMFunctionReturnHot(0);
}

ProjectionFieldOn10NodeTet::ProjectionFieldOn10NodeTet(Interface &m_field,
                                                       std::string _fieldName,
                                                       bool set_nodes,
                                                       bool on_coords,
                                                       std::string on_tag)
    : Projection10NodeCoordsOnField(m_field, _fieldName), setNodes(set_nodes),
      onCoords(on_coords), onTag(on_tag), maxApproximationOrder(20) {}

Tag th;
Field_multiIndex::index<FieldName_mi_tag>::type::iterator field_it;
VectorDouble L;
VectorDouble K;

MoFEMErrorCode ProjectionFieldOn10NodeTet::preProcess() {
  MoFEMFunctionBeginHot;
  if (!onCoords) {
    if (onTag == "NoNE") {
      SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
              "tag name not specified");
    }
    field_it = fieldsPtr->get<FieldName_mi_tag>().find(fieldName);
    if (field_it == fieldsPtr->get<FieldName_mi_tag>().end()) {
      SETERRQ1(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY, "field not found %s",
               fieldName.c_str());
    }
    int field_rank = (*field_it)->getNbOfCoeffs();
    VectorDouble def_VAL = ublas::zero_vector<double>(field_rank);
    CHKERR mField.get_moab().tag_get_handle(
        onTag.c_str(), field_rank, MB_TYPE_DOUBLE, th,
        MB_TAG_CREAT | MB_TAG_SPARSE, &*def_VAL.data().begin());
    MOAB_THROW(rval);
  }

  L.resize(maxApproximationOrder + 1);
  CHKERR Legendre_polynomials(maxApproximationOrder, 0., NULL,
                              &*L.data().begin(), NULL, 3);
  K.resize(10);
  CHKERR LobattoKernel_polynomials(9, 0., NULL, &*K.data().begin(), NULL, 3);

  MoFEMFunctionReturnHot(0);
}

MoFEMErrorCode ProjectionFieldOn10NodeTet::operator()() {
  MoFEMFunctionBeginHot;
  if (dofPtr->getName() != fieldName)
    MoFEMFunctionReturnHot(0);
  if (setNodes) {
    if (dofPtr->getEntType() == MBVERTEX) {
      EntityHandle node = dofPtr->getEnt();
      if (onCoords) {
        coords.resize(3);
        CHKERR mField.get_moab().get_coords(&node, 1, &*coords.data().begin());
        coords[dofPtr->getDofCoeffIdx()] = dofPtr->getFieldData();
        CHKERR mField.get_moab().set_coords(&node, 1, &*coords.data().begin());
      } else {
        int field_rank = (*field_it)->getNbOfCoeffs();
        if (field_rank != dofPtr->getNbOfCoeffs()) {
          SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
                  "data inconsistency");
        }
        double *tag_value;
        int tag_size;
        CHKERR mField.get_moab().tag_get_by_ptr(
            th, &node, 1, (const void **)&tag_value, &tag_size);
        if (tag_size != dofPtr->getNbOfCoeffs()) {
          SETERRQ(PETSC_COMM_SELF, 1, "data inconsistency");
        }
        tag_value[dofPtr->getDofCoeffIdx()] = dofPtr->getFieldData();
      }
    }
    MoFEMFunctionReturnHot(0);
  }
  if (dofPtr->getEntType() != MBEDGE) {
    MoFEMFunctionReturnHot(0);
  }
  EntityHandle edge = dofPtr->getEnt();
  if (type_from_handle(edge) != MBEDGE) {
    SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
            "this method works only elements which are type of MBEDGE");
  }

  int num_nodes;
  const EntityHandle *conn;
  CHKERR mField.get_moab().get_connectivity(edge, conn, num_nodes, false);
  if ((num_nodes != 2) && (num_nodes != 3)) {
    SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
            "this method works only 4 node and 10 node tets");
  }
  if (num_nodes == 2) {
    MoFEMFunctionReturnHot(0);
  }

  if (dofPtr->getDofOrder() >= maxApproximationOrder) {
    SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
            "too big approximation order, increase constant "
            "max_ApproximationOrder");
  }
  double approx_val = 0;
  FieldApproximationBase base = dofPtr->getApproxBase();
  switch (base) {
  case AINSWORTH_LEGENDRE_BASE:
    approx_val = 0.25 * L[dofPtr->getDofOrder() - 2] * dofPtr->getFieldData();
    break;
  case AINSWORTH_LOBATTO_BASE:
    approx_val = 0.25 * K[dofPtr->getDofOrder() - 2] * dofPtr->getFieldData();
    break;
  default:
    SETERRQ(PETSC_COMM_SELF, MOFEM_NOT_IMPLEMENTED, "not yet implemented");
  }

  if (onCoords) {
    coords.resize(num_nodes * 3);
    CHKERR mField.get_moab().get_coords(conn, num_nodes,
                                        &*coords.data().begin());
    if (dofPtr->getEntDofIdx() == dofPtr->getDofCoeffIdx()) {
      // add only one when higher order terms present
      double ave_mid = (coords[3 * 0 + dofPtr->getDofCoeffIdx()] +
                        coords[3 * 1 + dofPtr->getDofCoeffIdx()]) *
                       0.5;
      coords[2 * 3 + dofPtr->getDofCoeffIdx()] = ave_mid;
    }
    coords[2 * 3 + dofPtr->getDofCoeffIdx()] += approx_val;
    CHKERR mField.get_moab().set_coords(&conn[2], 1, &coords[3 * 2]);
  } else {
    int tag_size;
    double *tag_value[num_nodes];
    CHKERR mField.get_moab().tag_get_by_ptr(
        th, conn, num_nodes, (const void **)tag_value, &tag_size);
    if (tag_size != dofPtr->getNbOfCoeffs()) {
      SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY, "data inconsistency");
    }
    if (dofPtr->getEntDofIdx() == dofPtr->getDofCoeffIdx()) {
      // add only one when higher order terms present
      double ave_mid = (tag_value[0][dofPtr->getDofCoeffIdx()] +
                        tag_value[1][dofPtr->getDofCoeffIdx()]) *
                       0.5;
      tag_value[2][dofPtr->getDofCoeffIdx()] = ave_mid;
    }
    tag_value[2][dofPtr->getDofCoeffIdx()] += approx_val;
  }
  MoFEMFunctionReturnHot(0);
}
}; // namespace MoFEM