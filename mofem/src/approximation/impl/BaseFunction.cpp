/** \file BaseFunction.cpp
 * \brief implementation of multi-grid solver for p- adaptivity
 *
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
 * License along with MoFEM. If not, see <http://www.gnu.org/licenses/>
*/

#include <version.h>
#include <config.h>
#include <definitions.h>
#include <Includes.hpp>
#include <Common.hpp>
#include <UnknownInterface.hpp>
using namespace MoFEM;

#include <BaseFunction.hpp>

PetscErrorCode BaseFunctionCtx::queryInterface(
  const MOFEMuuid& uuid,MoFEM::UnknownInterface** iface
) {
  PetscFunctionBegin;
  *iface = NULL;
  if(uuid == IDD_UNKNOWN_BASE_FUNCTION) {
    *iface = static_cast<BaseFunctionCtx*>(this);
    PetscFunctionReturn(0);
  } else {
    SETERRQ(PETSC_COMM_WORLD,MOFEM_DATA_INCONSISTENCY,"wrong interference");
  }
  PetscFunctionReturn(0);
}

PetscErrorCode BaseFunction::queryInterface(
  const MOFEMuuid& uuid,MoFEM::UnknownInterface** iface
) {
  PetscFunctionBegin;
  *iface = NULL;
  if(uuid == IDD_UNKNOWN_BASE_FUNCTION) {
    *iface = static_cast<BaseFunction*>(this);
    PetscFunctionReturn(0);
  } else {
    SETERRQ(PETSC_COMM_WORLD,MOFEM_DATA_INCONSISTENCY,"wrong interference");
  }
  PetscFunctionReturn(0);
}

PetscErrorCode BaseFunction::getValue(
  ublas::matrix<double> &pts,
  boost::shared_ptr<BaseFunctionCtx> ctx_ptr
) {
  PetscFunctionBegin;
  SETERRQ(
    PETSC_COMM_SELF,
    MOFEM_NOT_IMPLEMENTED,
    "BaseFunction has not valid implementation of any shape function"
  );
  PetscFunctionReturn(0);
}
