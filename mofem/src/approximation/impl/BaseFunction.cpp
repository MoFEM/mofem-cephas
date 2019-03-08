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

namespace MoFEM {

MoFEMErrorCode BaseFunctionCtx::query_interface(
  const MOFEMuuid& uuid,MoFEM::UnknownInterface** iface
) const {
  MoFEMFunctionBeginHot;
  *iface = NULL;
  if(uuid == IDD_UNKNOWN_BASE_FUNCTION) {
    *iface = const_cast<BaseFunctionCtx*>(this);
    MoFEMFunctionReturnHot(0);
  } else {
    SETERRQ(PETSC_COMM_WORLD,MOFEM_DATA_INCONSISTENCY,"wrong interference");
  }
  MoFEMFunctionReturnHot(0);
}

MoFEMErrorCode BaseFunction::query_interface(
  const MOFEMuuid& uuid,MoFEM::UnknownInterface** iface
) const {
  MoFEMFunctionBeginHot;
  *iface = NULL;
  if(uuid == IDD_UNKNOWN_BASE_FUNCTION) {
    *iface = const_cast<BaseFunction*>(this);
    MoFEMFunctionReturnHot(0);
  } else {
    SETERRQ(PETSC_COMM_WORLD,MOFEM_DATA_INCONSISTENCY,"wrong interference");
  }
  MoFEMFunctionReturnHot(0);
}

MoFEMErrorCode BaseFunction::getValue(
  MatrixDouble &pts,
  boost::shared_ptr<BaseFunctionCtx> ctx_ptr
) {
  MoFEMFunctionBeginHot;
  SETERRQ(
    PETSC_COMM_SELF,
    MOFEM_NOT_IMPLEMENTED,
    "BaseFunction has not valid implementation of any shape function"
  );
  MoFEMFunctionReturnHot(0);
}

MoFEMErrorCode BaseFunction::getValue(
  MatrixDouble &pts_x,
  MatrixDouble &pts_t,
  boost::shared_ptr<BaseFunctionCtx> ctx_ptr
) {
  MoFEMFunctionBeginHot;
  SETERRQ(
    PETSC_COMM_SELF,
    MOFEM_NOT_IMPLEMENTED,
    "BaseFunction has not valid implementation of any shape function"
  );
  MoFEMFunctionReturnHot(0);
}

}