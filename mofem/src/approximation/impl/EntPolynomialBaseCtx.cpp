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

#include <version.h>
#include <config.h>
#include <definitions.h>
#include <Includes.hpp>

#include <base_functions.h>
#include <fem_tools.h>
#include <h1_hdiv_hcurl_l2.h>
#include <Common.hpp>
#include <UnknownInterface.hpp>
using namespace MoFEM;

#include <TagMultiIndices.hpp>
#include <CoordSysMultiIndices.hpp>
#include <FieldMultiIndices.hpp>
#include <EntsMultiIndices.hpp>
#include <DofsMultiIndices.hpp>
#include <FEMMultiIndices.hpp>
#include <DataStructures.hpp>

#include <BaseFunction.hpp>
#include <EntPolynomialBaseCtx.hpp>

PetscErrorCode EntPolynomialBaseCtx::queryInterface(
  const MOFEMuuid& uuid,MoFEM::UnknownInterface** iface
) {
  PetscErrorCode ierr;
  PetscFunctionBegin;
  *iface = NULL;
  if(
    uuid == IDD_TET_BASE_FUNCTION ||
    uuid == IDD_TRI_BASE_FUNCTION
  ) {
    *iface = dynamic_cast<EntPolynomialBaseCtx*>(this);
    PetscFunctionReturn(0);
  } else {
    SETERRQ(PETSC_COMM_WORLD,MOFEM_DATA_INCONSISTENCY,"wrong interference");
  }
  ierr = BaseFunctionCtx::queryInterface(uuid,iface); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

EntPolynomialBaseCtx::EntPolynomialBaseCtx(
  DataForcesAndSurcesCore &data,
  const FieldSpace space,
  const FieldApproximationBase base,
  const FieldApproximationBase copy_node_base
):
dAta(data),
sPace(space),
bAse(base),
copyNodeBase(copy_node_base) {
  switch(bAse) {
    case AINSWORTH_COLE_BASE:
    basePolynomials = Legendre_polynomials;
    break;
    case LOBATTO_BASE:
    basePolynomials = Lobatto_polynomials;
    break;
    default:
    THROW_MESSAGE("Not implemented for this base")
  }
}

EntPolynomialBaseCtx::~EntPolynomialBaseCtx() {
}
