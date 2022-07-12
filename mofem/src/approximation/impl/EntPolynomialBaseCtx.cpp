/** \file EntPolynomialBaseCtx.cpp
\brief Implementation of Ainsworth-Cole H1 base on tetrahedral

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
