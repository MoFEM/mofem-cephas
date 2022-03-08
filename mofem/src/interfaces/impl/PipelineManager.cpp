/** \file PipelineManager.cpp
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

MoFEMErrorCode
PipelineManager::query_interface(boost::typeindex::type_index type_index,
                                 UnknownInterface **iface) const {
  *iface = const_cast<PipelineManager *>(this);
  return 0;
}

PipelineManager::PipelineManager(const MoFEM::Core &core)
    : cOre(const_cast<Core &>(core)) {}

MoFEMErrorCode PipelineManager::loopFiniteElements(SmartPetscObj<DM> dm) {
  MoFEMFunctionBegin;
  Interface &m_field = cOre;
  Simple *simple = m_field.getInterface<Simple>();
  if (!dm)
    dm = simple->getDM();

  // Add element to calculate lhs of stiff part
  if (feDomainLhs)
    CHKERR DMoFEMLoopFiniteElements(dm, simple->getDomainFEName(), feDomainLhs);
  if (feBoundaryLhs)
    CHKERR DMoFEMLoopFiniteElements(dm, simple->getBoundaryFEName(),
                                    feBoundaryLhs);
  if (feSkeletonLhs)
    CHKERR DMoFEMLoopFiniteElements(dm, simple->getSkeletonFEName(),
                                    feSkeletonLhs);

  // Add element to calculate rhs of stiff part
  if (feDomainRhs)
    CHKERR DMoFEMLoopFiniteElements(dm, simple->getDomainFEName(), feDomainRhs);
  if (feBoundaryRhs)
    CHKERR DMoFEMLoopFiniteElements(dm, simple->getBoundaryFEName(),
                                    feBoundaryRhs);
  if (feSkeletonRhs)
    CHKERR DMoFEMLoopFiniteElements(dm, simple->getSkeletonFEName(),
                                    feSkeletonRhs);

  MoFEMFunctionReturn(0);
}

SmartPetscObj<KSP> PipelineManager::createKSP(SmartPetscObj<DM> dm) {
  Interface &m_field = cOre;
  Simple *simple = m_field.getInterface<Simple>();

  auto copy_dm_struture = [&](auto simple_dm) {
    MPI_Comm comm;
    CHKERR PetscObjectGetComm(getPetscObject(simple_dm.get()), &comm);
    DMType type;
    CHKERR DMGetType(simple_dm, &type);
    dm = createSmartDM(comm, type);
    CHKERR DMMoFEMDuplicateDMCtx(simple_dm, dm);
    return dm;
  };

  if (!dm)
    dm = copy_dm_struture(simple->getDM());
  else
    dm = copy_dm_struture(dm);

  const MoFEM::Problem *prb_ptr;
  CHKERR DMMoFEMGetProblemPtr(dm, &prb_ptr);

  auto set_dm_section = [&](auto dm) {
    MoFEMFunctionBegin;
    auto section =
        m_field.getInterface<ISManager>()->sectionCreate(prb_ptr->getName());
    CHKERR DMSetDefaultSection(dm, section);
    CHKERR DMSetDefaultGlobalSection(dm, section);
    MoFEMFunctionReturn(0);
  };
  CHKERR set_dm_section(dm);

  boost::shared_ptr<FEMethod> null;

  smartGetDMKspCtx(dm)->clearLoops();

  // Add element to calculate lhs of stiff part
  if (feDomainLhs)
    CHKERR DMMoFEMKSPSetComputeOperators(dm, simple->getDomainFEName(),
                                         feDomainLhs, null, null);
  if (feBoundaryLhs)
    CHKERR DMMoFEMKSPSetComputeOperators(dm, simple->getBoundaryFEName(),
                                         feBoundaryLhs, null, null);
  if (feSkeletonLhs)
    CHKERR DMMoFEMKSPSetComputeOperators(dm, simple->getSkeletonFEName(),
                                         feSkeletonLhs, null, null);

  // Add element to calculate rhs of stiff part
  if (feDomainRhs)
    CHKERR DMMoFEMKSPSetComputeRHS(dm, simple->getDomainFEName(), feDomainRhs,
                                   null, null);
  if (feBoundaryRhs)
    CHKERR DMMoFEMKSPSetComputeRHS(dm, simple->getBoundaryFEName(),
                                   feBoundaryRhs, null, null);
  if (feSkeletonRhs)
    CHKERR DMMoFEMKSPSetComputeRHS(dm, simple->getSkeletonFEName(),
                                   feSkeletonRhs, null, null);

  auto ksp = MoFEM::createKSP(m_field.get_comm());
  CHKERR KSPSetDM(ksp, dm);
  return ksp;
}

SmartPetscObj<SNES> PipelineManager::createSNES(SmartPetscObj<DM> dm) {
  Interface &m_field = cOre;
  Simple *simple = m_field.getInterface<Simple>();

  auto copy_dm_struture = [&](auto simple_dm) {
    MPI_Comm comm;
    CHKERR PetscObjectGetComm(getPetscObject(simple_dm.get()), &comm);
    DMType type;
    CHKERR DMGetType(simple_dm, &type);
    dm = createSmartDM(comm, type);
    CHKERR DMMoFEMDuplicateDMCtx(simple_dm, dm);
    return dm;
  };

  if (!dm)
    dm = copy_dm_struture(simple->getDM());
  else
    dm = copy_dm_struture(dm);

  const MoFEM::Problem *prb_ptr;
  CHKERR DMMoFEMGetProblemPtr(dm, &prb_ptr);

  auto set_dm_section = [&](auto dm) {
    MoFEMFunctionBegin;
    auto section =
        m_field.getInterface<ISManager>()->sectionCreate(prb_ptr->getName());
    CHKERR DMSetDefaultSection(dm, section);
    CHKERR DMSetDefaultGlobalSection(dm, section);
    MoFEMFunctionReturn(0);
  };
  CHKERR set_dm_section(dm);

  smartGetDMSnesCtx(dm)->clearLoops();

  boost::shared_ptr<FEMethod> null;

  // Add element to calculate lhs of stiff part
  if (feDomainLhs)
    CHKERR DMMoFEMSNESSetJacobian(dm, simple->getDomainFEName(), feDomainLhs,
                                  null, null);
  if (feBoundaryLhs)
    CHKERR DMMoFEMSNESSetJacobian(dm, simple->getBoundaryFEName(),
                                  feBoundaryLhs, null, null);
  if (feSkeletonLhs)
    CHKERR DMMoFEMSNESSetJacobian(dm, simple->getSkeletonFEName(),
                                  feSkeletonLhs, null, null);

  // Add element to calculate rhs of stiff part
  if (feDomainRhs)
    CHKERR DMMoFEMSNESSetFunction(dm, simple->getDomainFEName(), feDomainRhs,
                                  null, null);
  if (feBoundaryRhs)
    CHKERR DMMoFEMSNESSetFunction(dm, simple->getBoundaryFEName(),
                                  feBoundaryRhs, null, null);
  if (feSkeletonRhs)
    CHKERR DMMoFEMSNESSetFunction(dm, simple->getSkeletonFEName(),
                                  feSkeletonRhs, null, null);

  auto snes = MoFEM::createSNES(m_field.get_comm());
  CHKERR SNESSetDM(snes, dm);
  return snes;
}

SmartPetscObj<TS> PipelineManager::createTS(const TSType type,
                                            SmartPetscObj<DM> dm) {
  switch (type) {
  case EX:
    return createTSEX(dm);
    break;
  case IM:
    return createTSIM(dm);
    break;
  case IM2:
    return createTSIM2(dm);
    break;
  case IMEX:
    return createTSIMEX(dm);
    break;
  default:
    CHK_THROW_MESSAGE(MOFEM_NOT_IMPLEMENTED,
                      "TS solver handling not implemented");
    break;
  }
  return SmartPetscObj<TS>();
}

SmartPetscObj<TS> PipelineManager::createTSEX(SmartPetscObj<DM> dm) {
  Interface &m_field = cOre;
  Simple *simple = m_field.getInterface<Simple>();

  auto copy_dm_struture = [&](auto simple_dm) {
    MPI_Comm comm;
    CHKERR PetscObjectGetComm(getPetscObject(simple_dm.get()), &comm);
    DMType type;
    CHKERR DMGetType(simple_dm, &type);
    dm = createSmartDM(comm, type);
    CHKERR DMMoFEMDuplicateDMCtx(simple_dm, dm);
    return dm;
  };

  if (!dm)
    dm = copy_dm_struture(simple->getDM());
  else
    dm = copy_dm_struture(dm);

  const MoFEM::Problem *prb_ptr;
  CHKERR DMMoFEMGetProblemPtr(dm, &prb_ptr);

  auto set_dm_section = [&](auto dm) {
    MoFEMFunctionBegin;
    auto section =
        m_field.getInterface<ISManager>()->sectionCreate(prb_ptr->getName());
    CHKERR DMSetDefaultSection(dm, section);
    CHKERR DMSetDefaultGlobalSection(dm, section);
    MoFEMFunctionReturn(0);
  };
  CHKERR set_dm_section(dm);

  boost::shared_ptr<FEMethod> null;

  smartGetDMTsCtx(dm)->clearLoops();

  // Add element to calculate rhs of stiff part
  if (feDomainRhs)
    CHKERR DMMoFEMTSSetRHSFunction(dm, simple->getDomainFEName(),
                                   feDomainExplicitRhs, null, null);
  if (feBoundaryRhs)
    CHKERR DMMoFEMTSSetRHSFunction(dm, simple->getBoundaryFEName(),
                                   feBoundaryExplicitRhs, null, null);
  if (feSkeletonExplicitRhs)
    CHKERR DMMoFEMTSSetRHSFunction(dm, simple->getSkeletonFEName(),
                                   feSkeletonExplicitRhs, null, null);

  // Note: More cases for explit, and implicit time ingeration cases can be
  // implemented here.

  auto ts = MoFEM::createTS(m_field.get_comm());
  CHKERR TSSetDM(ts, dm);
  return ts;
}

SmartPetscObj<TS> PipelineManager::createTSIM(SmartPetscObj<DM> dm) {
  Interface &m_field = cOre;
  Simple *simple = m_field.getInterface<Simple>();

  auto copy_dm_struture = [&](auto simple_dm) {
    MPI_Comm comm;
    CHKERR PetscObjectGetComm(getPetscObject(simple_dm.get()), &comm);
    DMType type;
    CHKERR DMGetType(simple_dm, &type);
    dm = createSmartDM(comm, type);
    CHKERR DMMoFEMDuplicateDMCtx(simple_dm, dm);
    return dm;
  };

  if (!dm)
    dm = copy_dm_struture(simple->getDM());
  else
    dm = copy_dm_struture(dm);

  const MoFEM::Problem *prb_ptr;
  CHKERR DMMoFEMGetProblemPtr(dm, &prb_ptr);

  auto set_dm_section = [&](auto dm) {
    MoFEMFunctionBegin;
    auto section =
        m_field.getInterface<ISManager>()->sectionCreate(prb_ptr->getName());
    CHKERR DMSetDefaultSection(dm, section);
    CHKERR DMSetDefaultGlobalSection(dm, section);
    MoFEMFunctionReturn(0);
  };
  CHKERR set_dm_section(dm);

  boost::shared_ptr<FEMethod> null;

  // Add element to calculate lhs of stiff part
  if (feDomainLhs)
    CHKERR DMMoFEMTSSetIJacobian(dm, simple->getDomainFEName(), feDomainLhs,
                                 null, null);
  if (feBoundaryLhs)
    CHKERR DMMoFEMTSSetIJacobian(dm, simple->getBoundaryFEName(), feBoundaryLhs,
                                 null, null);
  if (feSkeletonLhs)
    CHKERR DMMoFEMTSSetIJacobian(dm, simple->getSkeletonFEName(), feSkeletonLhs,
                                 null, null);

  // Add element to calculate rhs of stiff part
  if (feDomainRhs)
    CHKERR DMMoFEMTSSetIFunction(dm, simple->getDomainFEName(), feDomainRhs,
                                 null, null);
  if (feBoundaryRhs)
    CHKERR DMMoFEMTSSetIFunction(dm, simple->getBoundaryFEName(), feBoundaryRhs,
                                 null, null);
  if (feSkeletonRhs)
    CHKERR DMMoFEMTSSetIFunction(dm, simple->getSkeletonFEName(), feSkeletonRhs,
                                 null, null);

  // Note: More cases for explit, and implicit time ingeration cases can be
  // implemented here.

  auto ts = MoFEM::createTS(m_field.get_comm());
  CHKERR TSSetDM(ts, dm);
  return ts;
}

SmartPetscObj<TS> PipelineManager::createTSIM2(SmartPetscObj<DM> dm) {
  Interface &m_field = cOre;
  Simple *simple = m_field.getInterface<Simple>();

  auto copy_dm_struture = [&](auto simple_dm) {
    MPI_Comm comm;
    CHKERR PetscObjectGetComm(getPetscObject(simple_dm.get()), &comm);
    DMType type;
    CHKERR DMGetType(simple_dm, &type);
    dm = createSmartDM(comm, type);
    CHKERR DMMoFEMDuplicateDMCtx(simple_dm, dm);
    return dm;
  };

  if (!dm)
    dm = copy_dm_struture(simple->getDM());
  else
    dm = copy_dm_struture(dm);

  const MoFEM::Problem *prb_ptr;
  CHKERR DMMoFEMGetProblemPtr(dm, &prb_ptr);

  auto set_dm_section = [&](auto dm) {
    MoFEMFunctionBegin;
    auto section =
        m_field.getInterface<ISManager>()->sectionCreate(prb_ptr->getName());
    CHKERR DMSetDefaultSection(dm, section);
    CHKERR DMSetDefaultGlobalSection(dm, section);
    MoFEMFunctionReturn(0);
  };
  CHKERR set_dm_section(dm);

  boost::shared_ptr<FEMethod> null;

  // Add element to calculate lhs of stiff part
  if (feDomainLhs)
    CHKERR DMMoFEMTSSetI2Jacobian(dm, simple->getDomainFEName(), feDomainLhs,
                                  null, null);
  if (feBoundaryLhs)
    CHKERR DMMoFEMTSSetI2Jacobian(dm, simple->getBoundaryFEName(),
                                  feBoundaryLhs, null, null);
  if (feSkeletonLhs)
    CHKERR DMMoFEMTSSetI2Jacobian(dm, simple->getSkeletonFEName(),
                                  feSkeletonLhs, null, null);

  // Add element to calculate rhs of stiff part
  if (feDomainRhs)
    CHKERR DMMoFEMTSSetI2Function(dm, simple->getDomainFEName(), feDomainRhs,
                                  null, null);
  if (feBoundaryRhs)
    CHKERR DMMoFEMTSSetI2Function(dm, simple->getBoundaryFEName(),
                                  feBoundaryRhs, null, null);
  if (feSkeletonRhs)
    CHKERR DMMoFEMTSSetI2Function(dm, simple->getSkeletonFEName(),
                                  feSkeletonRhs, null, null);

  // Note: More cases for explit, and implicit time ingeration cases can be
  // implemented here.

  auto ts = MoFEM::createTS(m_field.get_comm());
  CHKERR TSSetDM(ts, dm);
  return ts;
}

SmartPetscObj<TS> PipelineManager::createTSIMEX(SmartPetscObj<DM> dm) {
  Interface &m_field = cOre;
  Simple *simple = m_field.getInterface<Simple>();

  auto copy_dm_struture = [&](auto simple_dm) {
    MPI_Comm comm;
    CHKERR PetscObjectGetComm(getPetscObject(simple_dm.get()), &comm);
    DMType type;
    CHKERR DMGetType(simple_dm, &type);
    dm = createSmartDM(comm, type);
    CHKERR DMMoFEMDuplicateDMCtx(simple_dm, dm);
    return dm;
  };

  if (!dm)
    dm = copy_dm_struture(simple->getDM());
  else
    dm = copy_dm_struture(dm);

  const MoFEM::Problem *prb_ptr;
  CHKERR DMMoFEMGetProblemPtr(dm, &prb_ptr);

  auto set_dm_section = [&](auto dm) {
    MoFEMFunctionBegin;
    auto section =
        m_field.getInterface<ISManager>()->sectionCreate(prb_ptr->getName());
    CHKERR DMSetDefaultSection(dm, section);
    CHKERR DMSetDefaultGlobalSection(dm, section);
    MoFEMFunctionReturn(0);
  };
  CHKERR set_dm_section(dm);

  boost::shared_ptr<FEMethod> null;

  // Add element to calculate lhs of stiff part
  if (feDomainLhs)
    CHKERR DMMoFEMTSSetIJacobian(dm, simple->getDomainFEName(), feDomainLhs,
                                 null, null);
  if (feBoundaryLhs)
    CHKERR DMMoFEMTSSetIJacobian(dm, simple->getBoundaryFEName(), feBoundaryLhs,
                                 null, null);
  if (feSkeletonLhs)
    CHKERR DMMoFEMTSSetIJacobian(dm, simple->getSkeletonFEName(), feSkeletonLhs,
                                 null, null);

  // Add element to calculate rhs of stiff part
  if (feDomainRhs)
    CHKERR DMMoFEMTSSetIFunction(dm, simple->getDomainFEName(), feDomainRhs,
                                 null, null);
  if (feBoundaryRhs)
    CHKERR DMMoFEMTSSetIFunction(dm, simple->getBoundaryFEName(), feBoundaryRhs,
                                 null, null);
  if (feSkeletonRhs)
    CHKERR DMMoFEMTSSetIFunction(dm, simple->getSkeletonFEName(), feSkeletonRhs,
                                 null, null);

  // Add element to calculate rhs of stiff part
  if (feDomainExplicitRhs)
    CHKERR DMMoFEMTSSetRHSFunction(dm, simple->getDomainFEName(),
                                   feDomainExplicitRhs, null, null);
  if (feBoundaryExplicitRhs)
    CHKERR DMMoFEMTSSetRHSFunction(dm, simple->getBoundaryFEName(),
                                   feBoundaryExplicitRhs, null, null);
  if (feSkeletonExplicitRhs)
    CHKERR DMMoFEMTSSetRHSFunction(dm, simple->getSkeletonFEName(),
                                   feSkeletonExplicitRhs, null, null);

  // Note: More cases for explit, and implicit time ingeration cases can be
  // implemented here.

  auto ts = MoFEM::createTS(m_field.get_comm());
  CHKERR TSSetDM(ts, dm);
  return ts;
}

} // namespace MoFEM
