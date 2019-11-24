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

boost::shared_ptr<ForcesAndSourcesCore>
Basic::createDomainFEPipeline(boost::shared_ptr<ForcesAndSourcesCore> &fe,
                              const bool reset) {
  if (!fe || reset) {
    switch (getDim()) {
    case 2:
      fe = boost::make_shared<FaceElementForcesAndSourcesCore>(cOre);
      break;
    case 3:
      fe = boost::make_shared<VolumeElementForcesAndSourcesCore>(cOre);
      break;
    default:
      THROW_MESSAGE("Dimension not implemented Dim = " +
                    boost::lexical_cast<std::string>(getDim()));
    }
  }
  return fe;
}

boost::shared_ptr<ForcesAndSourcesCore>
Basic::createBoundaryFEPipeline(boost::shared_ptr<ForcesAndSourcesCore> &fe,
                                const bool reset) {
  if (!fe || reset) {
    switch (getDim()) {
    case 2:
      fe = boost::make_shared<EdgeElementForcesAndSourcesCore>(cOre);
      break;
    case 3:
      fe = boost::make_shared<FaceElementForcesAndSourcesCore>(cOre);
      break;
    default:
      THROW_MESSAGE("Dimension not implemented Dim = " +
                    boost::lexical_cast<std::string>(getDim()));
    }
  }
  return fe;
}

boost::ptr_vector<Basic::UserDataOperator> &
Basic::getOpDomainLhsPipeline(const bool reset) {
  return createDomainFEPipeline(feDomainLhs, reset)->getOpPtrVector();
}

boost::ptr_vector<Basic::UserDataOperator> &
Basic::getOpDomainRhsPipeline(const bool reset) {
  return createDomainFEPipeline(feDomainRhs, reset)->getOpPtrVector();
}

boost::ptr_vector<Basic::UserDataOperator> &
Basic::getOpBoundaryLhsPipeline(const bool reset) {
  return createBoundaryFEPipeline(feBcLhs, reset)->getOpPtrVector();
}

boost::ptr_vector<Basic::UserDataOperator> &
Basic::getOpBoundaryRhsPipeline(const bool reset) {
  return createBoundaryFEPipeline(feBcRhs, reset)->getOpPtrVector();
}

boost::ptr_vector<Basic::UserDataOperator> &
Basic::getOpSkeletonLhsPipeline(const bool reset) {
  return createBoundaryFEPipeline(feSkeletonLhs, reset)->getOpPtrVector();
}

boost::ptr_vector<Basic::UserDataOperator> &
Basic::getOpSkeletonRhsPipeline(const bool reset) {
  return createBoundaryFEPipeline(feSkeletonRhs, reset)->getOpPtrVector();
}

MoFEMErrorCode Basic::loopFiniteElements() {
  MoFEMFunctionBegin;

  // Add element to calculate lhs of stiff part
  if (feDomainLhs)
    CHKERR DMoFEMLoopFiniteElements(getDM(), getDomainFEName(), feDomainLhs);
  if (feBcLhs)
    CHKERR DMoFEMLoopFiniteElements(getDM(), getBoundaryFEName(), feBcLhs);
  if (feSkeletonLhs)
    CHKERR DMoFEMLoopFiniteElements(getDM(), getSkeletonFEName(),
                                    feSkeletonLhs);

  // Add element to calculate rhs of stiff part
  if (feDomainRhs)
    CHKERR DMoFEMLoopFiniteElements(getDM(), getDomainFEName(), feDomainRhs);
  if (feBcRhs)
    CHKERR DMoFEMLoopFiniteElements(getDM(), getBoundaryFEName(), feBcRhs);
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
  if (feBcLhs)
    CHKERR DMMoFEMKSPSetComputeOperators(getDM(), getBoundaryFEName(), feBcLhs,
                                         null, null);
  if (feSkeletonLhs)
    CHKERR DMMoFEMKSPSetComputeOperators(getDM(), getSkeletonFEName(),
                                         feSkeletonLhs, null, null);

  // Add element to calculate rhs of stiff part
  if (feDomainRhs)
    CHKERR DMMoFEMKSPSetComputeRHS(getDM(), getDomainFEName(), feDomainRhs,
                                   null, null);
  if (feBcRhs)
    CHKERR DMMoFEMKSPSetComputeRHS(getDM(), getBoundaryFEName(), feBcRhs, null,
                                   null);
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
  if (feBcLhs)
    CHKERR DMMoFEMSNESSetJacobian(getDM(), getBoundaryFEName(), feBcLhs, null,
                                  null);
  if (feSkeletonLhs)
    CHKERR DMMoFEMSNESSetJacobian(getDM(), getSkeletonFEName(), feSkeletonLhs,
                                  null, null);

  // Add element to calculate rhs of stiff part
  if (feDomainRhs)
    CHKERR DMMoFEMSNESSetFunction(getDM(), getDomainFEName(), feDomainRhs, null,
                                  null);
  if (feBcRhs)
    CHKERR DMMoFEMSNESSetFunction(getDM(), getBoundaryFEName(), feBcRhs, null,
                                  null);
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
  if (feBcLhs)
    CHKERR DMMoFEMTSSetIJacobian(getDM(), getBoundaryFEName(), feBcLhs, null,
                                 null);
  if (feSkeletonLhs)
    CHKERR DMMoFEMTSSetIJacobian(getDM(), getSkeletonFEName(), feSkeletonLhs,
                                 null, null);

  // Add element to calculate rhs of stiff part
  if (feDomainRhs)
    CHKERR DMMoFEMTSSetIFunction(getDM(), getDomainFEName(), feDomainRhs, null,
                                 null);
  if (feBcRhs)
    CHKERR DMMoFEMTSSetIFunction(getDM(), getBoundaryFEName(), feBcRhs, null,
                                 null);
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
