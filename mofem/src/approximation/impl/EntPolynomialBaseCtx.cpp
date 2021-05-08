/** \file EntPolynomialBaseCtx.cpp
\brief Implementation of Ainsworth-Cole H1 base on tetrahedral

*/

/* This file is part of MoFEM.
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
 * License along with MoFEM. If not, see <http://www.gnu.org/licenses/>. */

using namespace MoFEM;

MoFEMErrorCode EntPolynomialBaseCtx::query_interface(
    const MOFEMuuid &uuid, BaseFunctionUnknownInterface **iface) const {

  MoFEMFunctionBeginHot;
  *iface = NULL;
  if (uuid == IDD_TET_BASE_FUNCTION || uuid == IDD_TRI_BASE_FUNCTION ||
          uuid == IDD_EDGE_BASE_FUNCTION || uuid == IDD_QUAD_BASE_FUNCTION) {
    *iface = const_cast<EntPolynomialBaseCtx *>(this);
    MoFEMFunctionReturnHot(0);
  } else {
    SETERRQ(PETSC_COMM_WORLD, MOFEM_DATA_INCONSISTENCY, "wrong interference");
  }
  ierr = BaseFunctionCtx::query_interface(uuid, iface);
  CHKERRG(ierr);
  MoFEMFunctionReturnHot(0);
}

EntPolynomialBaseCtx::EntPolynomialBaseCtx(
    DataForcesAndSourcesCore &data, const FieldSpace space,
    const FieldApproximationBase base,
    const FieldApproximationBase copy_node_base)
    : dAta(data), sPace(space), bAse(base), copyNodeBase(copy_node_base) {
  ierr = setBase();
  CHKERRABORT(PETSC_COMM_WORLD, ierr);
}

EntPolynomialBaseCtx::EntPolynomialBaseCtx(
    DataForcesAndSourcesCore &data, const std::string field_name,
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
