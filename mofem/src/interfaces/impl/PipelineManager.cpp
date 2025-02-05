/** \file PipelineManager.cpp
 * \brief Implementation of basic interface
 * \ingroup mofem_basic_interface
 */


namespace MoFEM {

struct PipelineManager::MeshsetFE : public ForcesAndSourcesCore {
  using ForcesAndSourcesCore::ForcesAndSourcesCore;
  MoFEMErrorCode operator()() {
    MoFEMFunctionBegin;
    const auto type = numeredEntFiniteElementPtr->getEntType();
    if (type != MBENTITYSET) {
      SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
              "Expected entity set as a finite element");
    }
    CHKERR loopOverOperators();
    MoFEMFunctionReturn(0);
  }
  MoFEMErrorCode preProcess() { return 0; }
  MoFEMErrorCode postProcess() { return 0; }
};

boost::shared_ptr<FEMethod> &
PipelineManager::createMeshsetFEPipeline(boost::shared_ptr<FEMethod> &fe) {
  if (!fe)
    fe = boost::make_shared<MeshsetFE>(cOre);
  return fe;
}

boost::ptr_deque<PipelineManager::UserDataOperator> &
PipelineManager::getOpMeshsetRhsPipeline() {
  return boost::dynamic_pointer_cast<ForcesAndSourcesCore>(
             createMeshsetFEPipeline(feMeshsetRhs))
      ->getOpPtrVector();
}

boost::ptr_deque<PipelineManager::UserDataOperator> &
PipelineManager::getOpMeshsetLhsPipeline() {
  return boost::dynamic_pointer_cast<ForcesAndSourcesCore>(
             createMeshsetFEPipeline(feMeshsetLhs))
      ->getOpPtrVector();
}

boost::ptr_deque<PipelineManager::UserDataOperator> &
PipelineManager::getOpMeshsetExplicitRhsPipeline() {
  return boost::dynamic_pointer_cast<ForcesAndSourcesCore>(
             createMeshsetFEPipeline(feMeshsetExplicitRhs))
      ->getOpPtrVector();
}

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
  if (feMeshsetLhs)
    CHKERR DMoFEMLoopFiniteElements(dm, simple->getMeshsetFEName(),
                                    feMeshsetLhs);

  // Add element to calculate rhs of stiff part
  if (feDomainRhs)
    CHKERR DMoFEMLoopFiniteElements(dm, simple->getDomainFEName(), feDomainRhs);
  if (feBoundaryRhs)
    CHKERR DMoFEMLoopFiniteElements(dm, simple->getBoundaryFEName(),
                                    feBoundaryRhs);
  if (feSkeletonRhs)
    CHKERR DMoFEMLoopFiniteElements(dm, simple->getSkeletonFEName(),
                                    feSkeletonRhs);
  if (feMeshsetRhs)
    CHKERR DMoFEMLoopFiniteElements(dm, simple->getMeshsetFEName(),
                                    feMeshsetRhs);

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
    dm = createDM(comm, type);
    CHKERR DMMoFEMDuplicateDMCtx(simple_dm, dm);
    return dm;
  };

  if (!dm)
    dm = copy_dm_struture(simple->getDM());

  const MoFEM::Problem *prb_ptr;
  CHKERR DMMoFEMGetProblemPtr(dm, &prb_ptr);

  auto set_dm_section = [&](auto dm) {
    MoFEMFunctionBegin;
    auto section =
        m_field.getInterface<ISManager>()->sectionCreate(prb_ptr->getName());
    CHKERR DMSetSection(dm, section);
    MoFEMFunctionReturn(0);
  };
  CHKERR set_dm_section(dm);

  boost::shared_ptr<FEMethod> null;

  getDMKspCtx(dm)->clearLoops();

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
  if (feMeshsetLhs)
    CHKERR DMMoFEMKSPSetComputeOperators(dm, simple->getMeshsetFEName(),
                                         feMeshsetLhs, null, null);

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
  if (feMeshsetRhs)
    CHKERR DMMoFEMKSPSetComputeRHS(dm, simple->getMeshsetFEName(), feMeshsetRhs,
                                   null, null);

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
    dm = createDM(comm, type);
    CHKERR DMMoFEMDuplicateDMCtx(simple_dm, dm);
    return dm;
  };

  if (!dm)
    dm = copy_dm_struture(simple->getDM());

  const MoFEM::Problem *prb_ptr;
  CHKERR DMMoFEMGetProblemPtr(dm, &prb_ptr);

  auto set_dm_section = [&](auto dm) {
    MoFEMFunctionBegin;
    auto section =
        m_field.getInterface<ISManager>()->sectionCreate(prb_ptr->getName());
    CHKERR DMSetSection(dm, section);
    MoFEMFunctionReturn(0);
  };
  CHKERR set_dm_section(dm);

  getDMSnesCtx(dm)->clearLoops();

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
  if (feMeshsetLhs)
    CHKERR DMMoFEMSNESSetJacobian(dm, simple->getMeshsetFEName(), feMeshsetLhs,
                                  null, null);

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
  if (feMeshsetRhs)
    CHKERR DMMoFEMSNESSetFunction(dm, simple->getMeshsetFEName(), feMeshsetRhs,
                                  null, null);

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
    dm = createDM(comm, type);
    CHKERR DMMoFEMDuplicateDMCtx(simple_dm, dm);
    return dm;
  };

  if (!dm)
    dm = copy_dm_struture(simple->getDM());

  const MoFEM::Problem *prb_ptr;
  CHKERR DMMoFEMGetProblemPtr(dm, &prb_ptr);

  auto set_dm_section = [&](auto dm) {
    MoFEMFunctionBegin;
    auto section =
        m_field.getInterface<ISManager>()->sectionCreate(prb_ptr->getName());
    CHKERR DMSetSection(dm, section);
    MoFEMFunctionReturn(0);
  };
  CHKERR set_dm_section(dm);

  boost::shared_ptr<FEMethod> null;

  getDMTsCtx(dm)->clearLoops();

  // Add element to calculate rhs of slow part
  if (feDomainExplicitRhs)
    CHKERR DMMoFEMTSSetRHSFunction(dm, simple->getDomainFEName(),
                                   feDomainExplicitRhs, null, null);
  if (feBoundaryExplicitRhs)
    CHKERR DMMoFEMTSSetRHSFunction(dm, simple->getBoundaryFEName(),
                                   feBoundaryExplicitRhs, null, null);
  if (feSkeletonExplicitRhs)
    CHKERR DMMoFEMTSSetRHSFunction(dm, simple->getSkeletonFEName(),
                                   feSkeletonExplicitRhs, null, null);
  if (feMeshsetRhs)
    CHKERR DMMoFEMTSSetRHSFunction(dm, simple->getMeshsetFEName(),
                                   feMeshsetExplicitRhs, null, null);

  // Note: More cases for explict, and implicit time interation cases can be
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
    dm = createDM(comm, type);
    CHKERR DMMoFEMDuplicateDMCtx(simple_dm, dm);
    return dm;
  };

  if (!dm)
    dm = copy_dm_struture(simple->getDM());

  const MoFEM::Problem *prb_ptr;
  CHKERR DMMoFEMGetProblemPtr(dm, &prb_ptr);

  auto set_dm_section = [&](auto dm) {
    MoFEMFunctionBegin;
    auto section =
        m_field.getInterface<ISManager>()->sectionCreate(prb_ptr->getName());
    CHKERR DMSetSection(dm, section);
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
  if (feMeshsetLhs)
    CHKERR DMMoFEMTSSetIJacobian(dm, simple->getMeshsetFEName(), feMeshsetLhs,
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
  if (feMeshsetRhs)
    CHKERR DMMoFEMTSSetIFunction(dm, simple->getMeshsetFEName(), feMeshsetRhs,
                                 null, null);

  // Note: More cases for explict, and implicit time interation cases can be
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
    dm = createDM(comm, type);
    CHKERR DMMoFEMDuplicateDMCtx(simple_dm, dm);
    return dm;
  };

  if (!dm)
    dm = copy_dm_struture(simple->getDM());

  const MoFEM::Problem *prb_ptr;
  CHKERR DMMoFEMGetProblemPtr(dm, &prb_ptr);

  auto set_dm_section = [&](auto dm) {
    MoFEMFunctionBegin;
    auto section =
        m_field.getInterface<ISManager>()->sectionCreate(prb_ptr->getName());
    CHKERR DMSetSection(dm, section);
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
  if (feMeshsetLhs)
    CHKERR DMMoFEMTSSetI2Jacobian(dm, simple->getMeshsetFEName(), feMeshsetLhs,
                                  null, null);

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
  if (feMeshsetRhs)
    CHKERR DMMoFEMTSSetI2Function(dm, simple->getMeshsetFEName(), feMeshsetRhs,
                                  null, null);

  // Note: More cases for explict, and implicit time interation cases can be
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
    dm = createDM(comm, type);
    CHKERR DMMoFEMDuplicateDMCtx(simple_dm, dm);
    return dm;
  };

  if (!dm)
    dm = copy_dm_struture(simple->getDM());

  const MoFEM::Problem *prb_ptr;
  CHKERR DMMoFEMGetProblemPtr(dm, &prb_ptr);

  auto set_dm_section = [&](auto dm) {
    MoFEMFunctionBegin;
    auto section =
        m_field.getInterface<ISManager>()->sectionCreate(prb_ptr->getName());
    CHKERR DMSetSection(dm, section);
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
  if (feMeshsetLhs)
    CHKERR DMMoFEMTSSetIJacobian(dm, simple->getMeshsetFEName(), feMeshsetLhs,
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
  if (feMeshsetRhs)
    CHKERR DMMoFEMTSSetIFunction(dm, simple->getMeshsetFEName(), feMeshsetRhs,
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
  if (feMeshsetExplicitRhs)
    CHKERR DMMoFEMTSSetRHSFunction(dm, simple->getMeshsetFEName(),
                                   feMeshsetExplicitRhs, null, null);

  // Note: More cases for explict, and implicit time interation cases can be
  // implemented here.

  auto ts = MoFEM::createTS(m_field.get_comm());
  CHKERR TSSetDM(ts, dm);
  return ts;
}

} // namespace MoFEM
