/** \file EntPolynomialBaseCtx.cpp
\brief Implementation of Ainsworth-Cole H1 base on tetrahedral

*/



using namespace MoFEM;

MoFEMErrorCode
EntPolynomialBaseCtx::query_interface(boost::typeindex::type_index type_index,
                                      UnknownInterface **iface) const {
  *iface = const_cast<EntPolynomialBaseCtx *>(this);
  return 0;
}

EntPolynomialBaseCtx::EntPolynomialBaseCtx(
    EntitiesFieldData &data, const FieldSpace space,
    const FieldApproximationBase base,
    const FieldApproximationBase copy_node_base)
    : dAta(data), sPace(space), bAse(base), copyNodeBase(copy_node_base) {
  ierr = setBase();
  CHKERRABORT(PETSC_COMM_WORLD, ierr);
}

EntPolynomialBaseCtx::EntPolynomialBaseCtx(
    EntitiesFieldData &data, const std::string field_name,
    const FieldSpace space, const FieldApproximationBase base,
    const FieldApproximationBase copy_node_base)
    : dAta(data), sPace(space), bAse(base), fieldName(field_name),
      copyNodeBase(copy_node_base) {
  ierr = setBase();
  CHKERRABORT(PETSC_COMM_WORLD, ierr);
}

MoFEMErrorCode EntPolynomialBaseCtx::setBase() {
  MoFEMFunctionBeginHot;
  switch (bAse) {
  case AINSWORTH_LEGENDRE_BASE:
    switch (sPace) {
    case NOSPACE:
    case NOFIELD:
      SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY, "Makes no sense");
    case H1:
    case HCURL:
    case HDIV:
    case L2:
      basePolynomialsType0 = Legendre_polynomials;
      break;
    default:
      SETERRQ1(PETSC_COMM_SELF, MOFEM_NOT_IMPLEMENTED,
               "Not implemented for this space", FieldSpaceNames[sPace]);
    }
    break;
  case AINSWORTH_LOBATTO_BASE:
    switch (sPace) {
    case NOSPACE:
    case NOFIELD:
      SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY, "Makes no sense");
    case H1:
    case HCURL:
    case HDIV:
    case L2:
      basePolynomialsType0 = LobattoKernel_polynomials;
      break;
    default:
      SETERRQ1(PETSC_COMM_SELF, MOFEM_NOT_IMPLEMENTED,
               "Not implemented for this space", FieldSpaceNames[sPace]);
    }
    break;
  case AINSWORTH_BERNSTEIN_BEZIER_BASE:
    if(fieldName.empty())
      SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY, "Field name not set");
    switch (sPace) {
    case NOSPACE:
    case NOFIELD:
      SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY, "Makes no sense");
    case H1:
    case L2:
      break;
    default:
      SETERRQ1(PETSC_COMM_SELF, MOFEM_NOT_IMPLEMENTED,
               "Not implemented for this space", FieldSpaceNames[sPace]);
    }
    break;
  case DEMKOWICZ_JACOBI_BASE:
    switch (sPace) {
    case NOSPACE:
    case NOFIELD:
      SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY, "Makes no sense");
    case H1:
    case HCURL:
    case HDIV:
    case L2:
      basePolynomialsType1 = Jacobi_polynomials;
      break;
    default:
      SETERRQ1(PETSC_COMM_SELF, MOFEM_NOT_IMPLEMENTED,
               "Not implemented for this space",
               FieldSpaceNames[sPace]);
    }
    break;
  case USER_BASE:
    break;
  default:
    SETERRQ1(PETSC_COMM_SELF, MOFEM_NOT_IMPLEMENTED,
             "Not implemented for this base <%s>",
             ApproximationBaseNames[bAse]);
  }
  MoFEMFunctionReturnHot(0);
}
