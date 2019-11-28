/** \file Basic.cpp
 * \brief Implementation of basic interface
 * \ingroup mofem_basic_interface
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

MoFEMErrorCode Basic::query_interface(const MOFEMuuid &uuid,
                                      UnknownInterface **iface) const {
  MoFEMFunctionBeginHot;
  *iface = NULL;
  if (uuid == IDD_MOFEMBasic) {
    *iface = const_cast<Basic *>(this);
    MoFEMFunctionReturnHot(0);
  }
  if (uuid == IDD_MOFEMSimple) {
    *iface = static_cast<Simple *>(const_cast<Basic *>(this));
    MoFEMFunctionReturnHot(0);
  }
  SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY, "unknown interface");
  MoFEMFunctionReturnHot(0);
}

MoFEMErrorCode Basic::loopFiniteElements() {
  MoFEMFunctionBegin;

  // Add element to calculate lhs of stiff part
  if (feDomainLhs)
    CHKERR DMoFEMLoopFiniteElements(getDM(), getDomainFEName(), feDomainLhs);
  if (feBoundaryLhs)
    CHKERR DMoFEMLoopFiniteElements(getDM(), getBoundaryFEName(),
                                    feBoundaryLhs);
  if (feSkeletonLhs)
    CHKERR DMoFEMLoopFiniteElements(getDM(), getSkeletonFEName(),
                                    feSkeletonLhs);

  // Add element to calculate rhs of stiff part
  if (feDomainRhs)
    CHKERR DMoFEMLoopFiniteElements(getDM(), getDomainFEName(), feDomainRhs);
  if (feBoundaryRhs)
    CHKERR DMoFEMLoopFiniteElements(getDM(), getBoundaryFEName(),
                                    feBoundaryRhs);
  if (feSkeletonRhs)
    CHKERR DMoFEMLoopFiniteElements(getDM(), getSkeletonFEName(),
                                    feSkeletonRhs);

  MoFEMFunctionReturn(0);
}

SmartPetscObj<KSP> Basic::createKSP() {
  Interface &m_field = cOre;

  boost::shared_ptr<KspCtx> snes_ctx(new KspCtx(m_field, getProblemName()));
  CHKERR DMMoFEMSetKspCtx(getDM(), snes_ctx);

  boost::shared_ptr<FEMethod> null;

  // Add element to calculate lhs of stiff part
  if (feDomainLhs)
    CHKERR DMMoFEMKSPSetComputeOperators(getDM(), getDomainFEName(),
                                         feDomainLhs, null, null);
  if (feBoundaryLhs)
    CHKERR DMMoFEMKSPSetComputeOperators(getDM(), getBoundaryFEName(),
                                         feBoundaryLhs, null, null);
  if (feSkeletonLhs)
    CHKERR DMMoFEMKSPSetComputeOperators(getDM(), getSkeletonFEName(),
                                         feSkeletonLhs, null, null);

  // Add element to calculate rhs of stiff part
  if (feDomainRhs)
    CHKERR DMMoFEMKSPSetComputeRHS(getDM(), getDomainFEName(), feDomainRhs,
                                   null, null);
  if (feBoundaryRhs)
    CHKERR DMMoFEMKSPSetComputeRHS(getDM(), getBoundaryFEName(), feBoundaryRhs,
                                   null, null);
  if (feSkeletonRhs)
    CHKERR DMMoFEMKSPSetComputeRHS(getDM(), getSkeletonFEName(), feSkeletonRhs,
                                   null, null);

  auto ksp = MoFEM::createKSP(m_field.get_comm());
  CHKERR KSPSetDM(ksp, getDM());
  return ksp;
}

SmartPetscObj<SNES> Basic::createSNES() {
  Interface &m_field = cOre;

  boost::shared_ptr<MoFEM::SnesCtx> snes_ctx(
      new SnesCtx(m_field, getProblemName()));
  CHKERR DMMoFEMSetSnesCtx(getDM(), snes_ctx);

  boost::shared_ptr<FEMethod> null;

  // Add element to calculate lhs of stiff part
  if (feDomainLhs)
    CHKERR DMMoFEMSNESSetJacobian(getDM(), getDomainFEName(), feDomainLhs, null,
                                  null);
  if (feBoundaryLhs)
    CHKERR DMMoFEMSNESSetJacobian(getDM(), getBoundaryFEName(), feBoundaryLhs,
                                  null, null);
  if (feSkeletonLhs)
    CHKERR DMMoFEMSNESSetJacobian(getDM(), getSkeletonFEName(), feSkeletonLhs,
                                  null, null);

  // Add element to calculate rhs of stiff part
  if (feDomainRhs)
    CHKERR DMMoFEMSNESSetFunction(getDM(), getDomainFEName(), feDomainRhs, null,
                                  null);
  if (feBoundaryRhs)
    CHKERR DMMoFEMSNESSetFunction(getDM(), getBoundaryFEName(), feBoundaryRhs,
                                  null, null);
  if (feSkeletonRhs)
    CHKERR DMMoFEMSNESSetFunction(getDM(), getSkeletonFEName(), feSkeletonRhs,
                                  null, null);

  auto snes = MoFEM::createSNES(m_field.get_comm());
  CHKERR SNESSetDM(snes, getDM());
  return snes;
}

SmartPetscObj<TS> Basic::createTS() {
  Interface &m_field = cOre;

  boost::shared_ptr<MoFEM::TsCtx> ts_ctx(new TsCtx(m_field, getProblemName()));
  CHKERR DMMoFEMSetTsCtx(getDM(), ts_ctx);

  boost::shared_ptr<FEMethod> null;

  // Add element to calculate lhs of stiff part
  if (feDomainLhs)
    CHKERR DMMoFEMTSSetIJacobian(getDM(), getDomainFEName(), feDomainLhs, null,
                                 null);
  if (feBoundaryLhs)
    CHKERR DMMoFEMTSSetIJacobian(getDM(), getBoundaryFEName(), feBoundaryLhs,
                                 null, null);
  if (feSkeletonLhs)
    CHKERR DMMoFEMTSSetIJacobian(getDM(), getSkeletonFEName(), feSkeletonLhs,
                                 null, null);

  // Add element to calculate rhs of stiff part
  if (feDomainRhs)
    CHKERR DMMoFEMTSSetIFunction(getDM(), getDomainFEName(), feDomainRhs, null,
                                 null);
  if (feBoundaryRhs)
    CHKERR DMMoFEMTSSetIFunction(getDM(), getBoundaryFEName(), feBoundaryRhs,
                                 null, null);
  if (feSkeletonRhs)
    CHKERR DMMoFEMTSSetIFunction(getDM(), getSkeletonFEName(), feSkeletonRhs,
                                 null, null);

  // Note: More cases for explit, and implicit time ingeration cases can be
  // implemented here.

  auto ts = MoFEM::createTS(m_field.get_comm());
  CHKERR TSSetDM(ts, getDM());
  return ts;
}

} // namespace MoFEM
