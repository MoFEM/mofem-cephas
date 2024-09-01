/**
 * @file Electrostatics.cpp
 * \example Electrostatics.cpp
 *  */

#ifndef EXECUTABLE_DIMENSION
  #define EXECUTABLE_DIMENSION 3
#endif

#include <electrostatics.hpp>
static char help[] = "...\n\n";
struct Electrostatics {
public:
  Electrostatics(MoFEM::Interface &m_field);

  // Declaration of the main function to run analysis
  MoFEMErrorCode runProgram();

private:
  // Declaration of other main functions called in runProgram()
  MoFEMErrorCode readMesh();
  MoFEMErrorCode setupProblem();
  MoFEMErrorCode boundaryCondition();
  MoFEMErrorCode setIntegrationRules();
  MoFEMErrorCode assembleSystem();
  MoFEMErrorCode solveSystem();
  MoFEMErrorCode outputResults();
  MoFEMErrorCode getTotalEnergy();
  MoFEMErrorCode getElectrodeCharge();

  int oRder;
  MoFEM::Interface &mField;
  Simple *simpleInterface;
  std::string domainField;
  boost::shared_ptr<std::map<int, BlockData>> permBlockSetsPtr;
  boost::shared_ptr<std::map<int, BlockData>> intBlockSetsPtr;
  boost::shared_ptr<std::map<int, BlockData>> electrodeBlockSetsPtr;
  boost::shared_ptr<std::map<int, BlockData>> intrnlDomainBlockSetsPtr;
  boost::shared_ptr<DataAtIntegrationPts> commonDataPtr;

  boost::shared_ptr<ForcesAndSourcesCore> interface_rhs_fe;
  boost::shared_ptr<ForcesAndSourcesCore> electrode_rhs_fe;

  double ALPHA = 0.0;       // declaration for total charge on first electrode
  double BETA = 0.0;        // declaration for total charge on second electrode
  double totalEnergy = 0.0; // declaration for total energy
  SmartPetscObj<Vec> petscVec;       // petsc vector for the charge solution
  SmartPetscObj<Vec> petscVecEnergy; // petsc vector for the energy solution
  PetscBool out_skin = PETSC_FALSE;  //
  PetscBool out_volume = PETSC_FALSE;
  PetscBool is_partitioned = PETSC_FALSE;
  enum VecElements { ZERO = 0, ONE = 1, LAST_ELEMENT };
  int atom_test = 0;
};

Electrostatics::Electrostatics(MoFEM::Interface &m_field)
    : domainField("POTENTIAL"), mField(m_field) {}

//! [Read mesh]
MoFEMErrorCode Electrostatics::readMesh() {
  MoFEMFunctionBegin;
  CHKERR mField.getInterface(simpleInterface);
  CHKERR simpleInterface->getOptions();

  simpleInterface->getAddSkeletonFE() =
      true; // create lower dimensional element to ensure sharing of partitioed
            // entities at the boundaries.
  CHKERR simpleInterface->loadFile();
  MoFEMFunctionReturn(0);
}
//! [Read mesh]

//! [Setup problem]
MoFEMErrorCode Electrostatics::setupProblem() {
  MoFEMFunctionBegin;

  Range domain_ents;
  CHKERR mField.get_moab().get_entities_by_dimension(0, SPACE_DIM, domain_ents,
                                                     true);
  auto get_ents_by_dim = [&](const auto dim) {
    if (dim == SPACE_DIM) {
      return domain_ents;
    } else {
      Range ents;
      if (dim == 0)
        CHKERR mField.get_moab().get_connectivity(domain_ents, ents, true);
      else
        CHKERR mField.get_moab().get_entities_by_dimension(0, dim, ents, true);
      return ents;
    }
  };

  // Select base for the field based on the element type
  auto get_base = [&]() {
    auto domain_ents = get_ents_by_dim(SPACE_DIM);
    if (domain_ents.empty())
      CHK_THROW_MESSAGE(MOFEM_NOT_FOUND, "Empty mesh");
    const auto type = type_from_handle(domain_ents[0]);
    switch (type) {
    case MBQUAD:
      return DEMKOWICZ_JACOBI_BASE;
    case MBHEX:
      return DEMKOWICZ_JACOBI_BASE;
    case MBTRI:
      return AINSWORTH_LEGENDRE_BASE;
    case MBTET:
      return AINSWORTH_LEGENDRE_BASE;
    default:
      CHK_THROW_MESSAGE(MOFEM_NOT_FOUND, "Element type is not handled");
    }
    return NOBASE;
  };

  auto base = get_base();
  CHKERR simpleInterface->addDomainField(domainField, H1, base, 1);
  CHKERR simpleInterface->addBoundaryField(domainField, H1, base, 1);

  int oRder = 2; // default order

  CHKERR PetscOptionsGetInt(PETSC_NULL, "", "-order", &oRder, PETSC_NULL);
  CHKERR simpleInterface->setFieldOrder(domainField, oRder);
  CHKERR PetscOptionsGetInt(PETSC_NULL, "", "-atom_test", &atom_test,
                            PETSC_NULL);
  commonDataPtr = boost::make_shared<DataAtIntegrationPts>(mField);

  // gets the map of the permittivity attributes and the block sets
  permBlockSetsPtr = boost::make_shared<std::map<int, BlockData>>();
  Range mat_electr_ents; // range of entities with the permittivity
  for (_IT_CUBITMESHSETS_BY_SET_TYPE_FOR_LOOP_(mField, BLOCKSET, bit)) {
    if (bit->getName().compare(0, 12, "MAT_ELECTRIC") == 0) {
      const int id = bit->getMeshsetId();
      auto &block_data = (*permBlockSetsPtr)[id];

      CHKERR mField.get_moab().get_entities_by_dimension(
          bit->getMeshset(), SPACE_DIM, block_data.domainEnts, true);
      mat_electr_ents.merge(block_data.domainEnts);

      std::vector<double> attributes;
      bit->getAttributes(attributes);
      if (attributes.size() < 1) {
        SETERRQ1(PETSC_COMM_WORLD, MOFEM_DATA_INCONSISTENCY,
                 " At least one permittivity attributes should be given but "
                 "found %d",
                 attributes.size());
      }
      block_data.iD = id;                   // id of the block
      block_data.epsPermit = attributes[0]; // permittivity value of the block
    }
  }

  // gets the map of the charge attributes and the block sets
  intBlockSetsPtr = boost::make_shared<std::map<int, BlockData>>();
  Range int_electr_ents; // range of entities with the charge
  for (_IT_CUBITMESHSETS_BY_SET_TYPE_FOR_LOOP_(mField, BLOCKSET, bit)) {
    if (bit->getName().compare(0, 12, "INT_ELECTRIC") == 0) {
      const int id = bit->getMeshsetId();
      auto &block_data = (*intBlockSetsPtr)[id];

      CHKERR mField.get_moab().get_entities_by_dimension(
          bit->getMeshset(), SPACE_DIM - 1, block_data.interfaceEnts, true);
      int_electr_ents.merge(block_data.interfaceEnts);

      std::vector<double> attributes;
      bit->getAttributes(attributes);
      if (attributes.size() < 1) {
        SETERRQ1(PETSC_COMM_WORLD, MOFEM_DATA_INCONSISTENCY,
                 "At least one charge attributes should be given but found %d",
                 attributes.size());
      }

      block_data.iD = id;                       // id of the block
      block_data.chargeDensity = attributes[0]; // charge value of the block
    }
  }
  // gets the map of the electrode entity range in the  block sets
  electrodeBlockSetsPtr = boost::make_shared<std::map<int, BlockData>>();
  Range electrode_ents; // range of entities with the electrode
  int electrodeCount = 0;
  for (_IT_CUBITMESHSETS_BY_SET_TYPE_FOR_LOOP_(mField, BLOCKSET, bit)) {
    if (bit->getName().compare(0, 9, "ELECTRODE") == 0) {
      const int id = bit->getMeshsetId();
      auto &block_data = (*electrodeBlockSetsPtr)[id];
      ++electrodeCount;

      CHKERR mField.get_moab().get_entities_by_dimension(
          bit->getMeshset(), SPACE_DIM - 1, block_data.electrodeEnts, true);
      electrode_ents.merge(block_data.electrodeEnts);
      block_data.iD = id;
      auto print_range_on_procs = [&](const std::string &name, int meshsetId,
                                      const Range &range) {
        MOFEM_LOG("SYNC", Sev::inform)
            << name << " in meshID: " << id << " with range "
            << block_data.electrodeEnts << " on proc ["
            << mField.get_comm_rank() << "] \n";
        MOFEM_LOG_SYNCHRONISE(mField.get_comm());
      };
      print_range_on_procs(bit->getName(), id, electrode_ents);
      if (electrodeCount > 2) {
        SETERRQ(PETSC_COMM_SELF, MOFEM_INVALID_DATA,
                "Three or more electrode blocksets found");
        ;
      }
    }
  }

  // gets the map of the internal domain entity range to get the total energy
  intrnlDomainBlockSetsPtr = boost::make_shared<std::map<int, BlockData>>();
  Range internal_domain; // range of entities marked the internal domain
  for (_IT_CUBITMESHSETS_BY_SET_TYPE_FOR_LOOP_(mField, BLOCKSET, bit)) {
    if (bit->getName().compare(0, 10, "DOMAIN_INT") == 0) {
      const int id = bit->getMeshsetId();
      auto &block_data = (*intrnlDomainBlockSetsPtr)[id];

      CHKERR mField.get_moab().get_entities_by_dimension(
          bit->getMeshset(), SPACE_DIM, block_data.internalDomainEnts, true);
      internal_domain.merge(block_data.internalDomainEnts);
      block_data.iD = id;
    }
  }
  // sync entities
  CHKERR mField.getInterface<CommInterface>()->synchroniseEntities(
      mat_electr_ents);
  CHKERR mField.getInterface<CommInterface>()->synchroniseEntities(
      int_electr_ents);
  CHKERR mField.getInterface<CommInterface>()->synchroniseEntities(
      electrode_ents);
  CHKERR mField.getInterface<CommInterface>()->synchroniseEntities(
      internal_domain);

  CHKERR PetscOptionsGetBool(PETSC_NULL, "", "-out_skin", &out_skin,
                             PETSC_NULL);
  CHKERR PetscOptionsGetBool(PETSC_NULL, "", "-out_volume", &out_volume,
                             PETSC_NULL);
  CHKERR PetscOptionsGetBool(PETSC_NULL, "", "-is_partitioned", &is_partitioned,
                             PETSC_NULL);
  // get the skin entities
  Skinner skinner(&mField.get_moab());
  Range skin_tris;
  CHKERR skinner.find_skin(0, mat_electr_ents, false, skin_tris);
  Range proc_skin;
  ParallelComm *pcomm =
      ParallelComm::get_pcomm(&mField.get_moab(), MYPCOMM_INDEX);
  if (is_partitioned) {
    CHKERR pcomm->filter_pstatus(skin_tris,
                                 PSTATUS_SHARED | PSTATUS_MULTISHARED,
                                 PSTATUS_NOT, -1, &proc_skin);
  } else {
    proc_skin = skin_tris;
  }
  // add the skin entities to the field
  CHKERR mField.add_finite_element("SKIN", MF_ZERO);
  CHKERR mField.add_ents_to_finite_element_by_dim(proc_skin, SPACE_DIM - 1,
                                                  "SKIN");
  CHKERR mField.modify_finite_element_add_field_row("SKIN", domainField);
  CHKERR mField.modify_finite_element_add_field_col("SKIN", domainField);
  CHKERR mField.modify_finite_element_add_field_data("SKIN", domainField);
  // add the interface entities to the field
  CHKERR mField.add_finite_element("INTERFACE");
  CHKERR mField.modify_finite_element_add_field_row("INTERFACE", domainField);
  CHKERR mField.modify_finite_element_add_field_col("INTERFACE", domainField);
  CHKERR mField.modify_finite_element_add_field_data("INTERFACE", domainField);
  CHKERR mField.add_ents_to_finite_element_by_dim(int_electr_ents,
                                                  SPACE_DIM - 1, "INTERFACE");
  // add the electrode entities to the field
  CHKERR mField.add_finite_element("ELECTRODE");
  CHKERR mField.modify_finite_element_add_field_row("ELECTRODE", domainField);
  CHKERR mField.modify_finite_element_add_field_col("ELECTRODE", domainField);
  CHKERR mField.modify_finite_element_add_field_data("ELECTRODE", domainField);
  CHKERR mField.add_ents_to_finite_element_by_dim(electrode_ents, SPACE_DIM - 1,
                                                  "ELECTRODE");

  // sync field entities
  mField.getInterface<CommInterface>()->synchroniseFieldEntities(domainField);
  CHKERR simpleInterface->defineFiniteElements();
  CHKERR simpleInterface->defineProblem(PETSC_TRUE);
  CHKERR simpleInterface->buildFields();
  CHKERR simpleInterface->buildFiniteElements();
  CHKERR simpleInterface->buildProblem();
  CHKERR mField.build_finite_elements("SKIN");
  CHKERR mField.build_finite_elements("INTERFACE");
  CHKERR mField.build_finite_elements("ELECTRODE");

  CHKERR DMMoFEMAddElement(simpleInterface->getDM(), "SKIN");
  CHKERR DMMoFEMAddElement(simpleInterface->getDM(), "INTERFACE");
  CHKERR DMMoFEMAddElement(simpleInterface->getDM(), "ELECTRODE");

  DMType dm_name = "DMMOFEM";
  CHKERR DMRegister_MoFEM(dm_name);

  SmartPetscObj<DM> dm;
  dm = createDM(mField.get_comm(), dm_name);

  // create dm instance
  CHKERR DMSetType(dm, dm_name);

  CHKERR simpleInterface->buildProblem();

  // initialise petsc vector for required processor
  int local_size;
  if (mField.get_comm_rank() == 0) // get_comm_rank() gets processor number

    local_size = LAST_ELEMENT; // last element gives size of vector

  else
    // other processors (e.g. 1, 2, 3, etc.)
    local_size = 0; // local size of vector is zero on other processors

  petscVec = createVectorMPI(mField.get_comm(), local_size, LAST_ELEMENT);
  petscVecEnergy = createVectorMPI(mField.get_comm(), local_size, LAST_ELEMENT);
  MoFEMFunctionReturn(0);
}
//! [Setup problem]

//! [Boundary condition]
MoFEMErrorCode Electrostatics::boundaryCondition() {
  MoFEMFunctionBegin;

  auto bc_mng = mField.getInterface<BcManager>();

  // Remove_BCs_from_blockset name "BOUNDARY_CONDITION";
  CHKERR bc_mng->removeBlockDOFsOnEntities<BcScalarMeshsetType<BLOCKSET>>(
      simpleInterface->getProblemName(), "BOUNDARY_CONDITION",
      std::string(domainField), true);

  MoFEMFunctionReturn(0);
}
//! [Boundary condition]

//! [Set integration rules]
MoFEMErrorCode Electrostatics::setIntegrationRules() {
  MoFEMFunctionBegin;

  auto rule_lhs = [](int, int, int p) -> int { return 2 * p; };
  auto rule_rhs = [](int, int, int p) -> int { return p; };

  auto pipeline_mng = mField.getInterface<PipelineManager>();
  CHKERR pipeline_mng->setDomainLhsIntegrationRule(rule_lhs);
  CHKERR pipeline_mng->setDomainRhsIntegrationRule(rule_rhs);

  MoFEMFunctionReturn(0);
}
//! [Set integration rules]

//! [Assemble system]
MoFEMErrorCode Electrostatics::assembleSystem() {
  MoFEMFunctionBegin;

  auto pipeline_mng = mField.getInterface<PipelineManager>();

  commonDataPtr = boost::make_shared<DataAtIntegrationPts>(mField);
  auto add_base_ops = [&](auto &pipeline) {
    CHKERR AddHOOps<SPACE_DIM, SPACE_DIM, SPACE_DIM>::add(pipeline, {H1});

    pipeline.push_back(
        new OpBlockPermittivity(commonDataPtr, permBlockSetsPtr, domainField));
  };

  add_base_ops(pipeline_mng->getOpDomainLhsPipeline());
  auto epsilon = [&](const double, const double, const double) {
    return commonDataPtr->blockPermittivity;
  };

  { // Push operators to the Pipeline that is responsible for calculating LHS
    pipeline_mng->getOpDomainLhsPipeline().push_back(
        new OpDomainLhsMatrixK(domainField, domainField, epsilon));
  }

  { // Push operators to the Pipeline that is responsible for calculating RHS
    auto set_values_to_bc_dofs = [&](auto &fe) {
      auto get_bc_hook = [&]() {
        EssentialPreProc<TemperatureCubitBcData> hook(mField, fe, {});
        return hook;
      };
      fe->preProcessHook = get_bc_hook();
    };
    // Set essential BC
    auto calculate_residual_from_set_values_on_bc = [&](auto &pipeline) {
      using OpInternal =
          FormsIntegrators<DomainEleOp>::Assembly<PETSC>::LinearForm<
              GAUSS>::OpGradTimesTensor<BASE_DIM, FIELD_DIM, SPACE_DIM>;

      auto gradUValsPtr = boost::make_shared<MatrixDouble>();
      pipeline_mng->getOpDomainRhsPipeline().push_back(
          new OpCalculateScalarFieldGradient<SPACE_DIM>(domainField,
                                                        gradUValsPtr));

      add_base_ops(pipeline_mng->getOpDomainRhsPipeline());

      auto minus_epsilon = [&](double, double, double) constexpr {
        return -commonDataPtr->blockPermittivity;
      };
      pipeline_mng->getOpDomainRhsPipeline().push_back(
          new OpInternal(domainField, gradUValsPtr, minus_epsilon));
    };

    set_values_to_bc_dofs(pipeline_mng->getDomainRhsFE());
    calculate_residual_from_set_values_on_bc(
        pipeline_mng->getOpDomainRhsPipeline());
    auto bodySourceTerm = [&](const double, const double, const double) {
      return bodySource;
    };
    pipeline_mng->getOpDomainRhsPipeline().push_back(
        new OpBodySourceVectorb(domainField, bodySourceTerm));
  }

  interface_rhs_fe = boost::shared_ptr<ForcesAndSourcesCore>(
      new intElementForcesAndSourcesCore(mField));
  interface_rhs_fe->getRuleHook = [](int, int, int p) { return p; };

  {
    interface_rhs_fe->getOpPtrVector().push_back(
        new OpBlockChargeDensity(commonDataPtr, intBlockSetsPtr, domainField));

    auto sIgma = [&](const double, const double, const double) {
      return commonDataPtr->blockChrgDens;
    };

    interface_rhs_fe->getOpPtrVector().push_back(
        new OpInterfaceRhsVectorF(domainField, sIgma));
  }

  MoFEMFunctionReturn(0);
}
//! [Assemble system]

//! [Solve system]
MoFEMErrorCode Electrostatics::solveSystem() {
  MoFEMFunctionBegin;

  auto pipeline_mng = mField.getInterface<PipelineManager>();

  auto ksp_solver = pipeline_mng->createKSP();

  boost::shared_ptr<ForcesAndSourcesCore> null; ///< Null element does
  CHKERR DMMoFEMKSPSetComputeRHS(simpleInterface->getDM(), "INTERFACE",
                                 interface_rhs_fe, null, null);

  CHKERR KSPSetFromOptions(ksp_solver);
  CHKERR KSPSetUp(ksp_solver);

  // Create RHS and solution vectors
  auto dm = simpleInterface->getDM();
  auto F = createDMVector(dm);
  auto D = vectorDuplicate(F);
  // Solve the system
  CHKERR KSPSolve(ksp_solver, F, D);

  CHKERR VecGhostUpdateBegin(F, INSERT_VALUES, SCATTER_FORWARD);
  CHKERR VecGhostUpdateEnd(F, INSERT_VALUES, SCATTER_FORWARD);

  double fnorm;
  CHKERR VecNorm(F, NORM_2, &fnorm);
  CHKERR PetscPrintf(PETSC_COMM_WORLD, "F norm  = %9.8e\n", fnorm);

  // Scatter result data on the mesh
  CHKERR VecGhostUpdateBegin(D, INSERT_VALUES, SCATTER_FORWARD);
  CHKERR VecGhostUpdateEnd(D, INSERT_VALUES, SCATTER_FORWARD);

  double dnorm;
  CHKERR VecNorm(D, NORM_2, &dnorm);
  CHKERR PetscPrintf(PETSC_COMM_WORLD, "D norm  = %9.8e\n", dnorm);

  CHKERR DMoFEMMeshToLocalVector(dm, D, INSERT_VALUES, SCATTER_REVERSE);

  MoFEMFunctionReturn(0);
}
//! [Solve system]

//! [Output results]
MoFEMErrorCode Electrostatics::outputResults() {
  MoFEMFunctionBegin;
  auto pipeline_mng = mField.getInterface<PipelineManager>();
  pipeline_mng->getDomainRhsFE().reset();
  pipeline_mng->getBoundaryLhsFE().reset();
  pipeline_mng->getBoundaryRhsFE().reset(); //
  pipeline_mng->getDomainLhsFE().reset();
  auto post_proc_fe = boost::make_shared<PostProcEle>(mField);
  // add higher order operator
  CHKERR AddHOOps<SPACE_DIM, SPACE_DIM, SPACE_DIM>::add(
      post_proc_fe->getOpPtrVector(), {H1});

  auto uPtr = boost::make_shared<VectorDouble>();
  auto gradUPtr = boost::make_shared<MatrixDouble>();
  post_proc_fe->getOpPtrVector().push_back(
      new OpCalculateScalarFieldValues(domainField, uPtr));

  post_proc_fe->getOpPtrVector().push_back(
      new OpCalculateScalarFieldGradient<SPACE_DIM>(domainField, gradUPtr));

  using OpPPMap = OpPostProcMapInMoab<SPACE_DIM, SPACE_DIM>;

  auto negGradUPtr = boost::make_shared<MatrixDouble>();
  auto negGradUPtrTimesPerm = boost::make_shared<MatrixDouble>();
  auto energyDensityPtr = boost::make_shared<VectorDouble>();
  post_proc_fe->getOpPtrVector().push_back(
      new OpElectricField(negGradUPtr, gradUPtr));

  post_proc_fe->getOpPtrVector().push_back(
      new OpGradTimesPerm(domainField, negGradUPtr, negGradUPtrTimesPerm,
                          permBlockSetsPtr, commonDataPtr));
  post_proc_fe->getOpPtrVector().push_back(
      new OpEnergyDensity(domainField, negGradUPtr, energyDensityPtr,
                          permBlockSetsPtr, commonDataPtr));

  post_proc_fe->getOpPtrVector().push_back(

      new OpPPMap(post_proc_fe->getPostProcMesh(),
                  post_proc_fe->getMapGaussPts(),

                  OpPPMap::DataMapVec{{"POTENTIAL", uPtr},
                                      {"ENERGY_DENSITY", energyDensityPtr}},
                  OpPPMap::DataMapMat{
                      {"ELECTRIC_FIELD", negGradUPtr},
                      {"ELECTRIC_DISPLACEMENT", negGradUPtrTimesPerm},
                  },
                  OpPPMap::DataMapMat{},

                  OpPPMap::DataMapMat{})

  );

  pipeline_mng->getDomainRhsFE() = post_proc_fe;
  CHKERR pipeline_mng->loopFiniteElements();
  CHKERR post_proc_fe->writeFile("out.h5m");

  if (out_skin && SPACE_DIM == 3) {
    auto post_proc_skin = boost::make_shared<PostProcFaceEle>(mField);
    auto op_loop_side = new OpLoopSide<SideEle>(
        mField, simpleInterface->getDomainFEName(), SPACE_DIM);
    CHKERR AddHOOps<SPACE_DIM, SPACE_DIM, SPACE_DIM>::add(
        op_loop_side->getOpPtrVector(), {H1});

    auto uPtrSkin = boost::make_shared<VectorDouble>();
    auto gradUPtrSkin = boost::make_shared<MatrixDouble>();
    auto negGradUSkinPtrSkin = boost::make_shared<MatrixDouble>();
    auto elecDisplacementSkin = boost::make_shared<MatrixDouble>();
    auto energyDensityPtrSkin = boost::make_shared<VectorDouble>();

    op_loop_side->getOpPtrVector().push_back(
        new OpCalculateScalarFieldValues(domainField, uPtrSkin));
    op_loop_side->getOpPtrVector().push_back(
        new OpCalculateScalarFieldGradient<SPACE_DIM>(domainField,
                                                      gradUPtrSkin));
    op_loop_side->getOpPtrVector().push_back(
        new OpElectricField(negGradUSkinPtrSkin, gradUPtrSkin));
    op_loop_side->getOpPtrVector().push_back(new OpGradTimesPerm(
        domainField, negGradUSkinPtrSkin, elecDisplacementSkin,
        permBlockSetsPtr, commonDataPtr));
    op_loop_side->getOpPtrVector().push_back(new OpEnergyDensity(
        domainField, negGradUSkinPtrSkin, energyDensityPtrSkin,
        permBlockSetsPtr, commonDataPtr));

    // push op to boundary element
    post_proc_skin->getOpPtrVector().push_back(op_loop_side);

    post_proc_skin->getOpPtrVector().push_back(new OpPPMap(
        post_proc_skin->getPostProcMesh(), post_proc_skin->getMapGaussPts(),
        OpPPMap::DataMapVec{{"POTENTIAL", uPtrSkin},
                            {"ENERGY_DENSITY", energyDensityPtrSkin}},
        OpPPMap::DataMapMat{{"ELECTRIC_FIELD", negGradUSkinPtrSkin},
                            {"ELECTRIC_DISPLACEMENT", elecDisplacementSkin}},
        OpPPMap::DataMapMat{}, OpPPMap::DataMapMat{}));

    CHKERR DMoFEMLoopFiniteElements(simpleInterface->getDM(), "SKIN",
                                    post_proc_skin);
    CHKERR post_proc_skin->writeFile("out_skin.h5m");
  }
  MoFEMFunctionReturn(0);
}
//! [Output results]

//! [Get Total Energy]
MoFEMErrorCode Electrostatics::getTotalEnergy() {
  MoFEMFunctionBegin;
  auto pipeline_energy = mField.getInterface<PipelineManager>();
  pipeline_energy->getDomainRhsFE().reset();
  pipeline_energy->getDomainLhsFE().reset();
  pipeline_energy->getBoundaryLhsFE().reset();
  pipeline_energy->getBoundaryRhsFE().reset();
  pipeline_energy->getOpDomainRhsPipeline().clear();
  pipeline_energy->getOpDomainLhsPipeline().clear();

  CHKERR AddHOOps<SPACE_DIM, SPACE_DIM, SPACE_DIM>::add(
      pipeline_energy->getOpDomainRhsPipeline(), {H1});

  auto uPtr = boost::make_shared<VectorDouble>();
  auto gradUValsPtr = boost::make_shared<MatrixDouble>();
  auto negGradUPtr = boost::make_shared<MatrixDouble>();

  pipeline_energy->getOpDomainRhsPipeline().push_back(
      new OpCalculateScalarFieldValues(domainField, uPtr));
  pipeline_energy->getOpDomainRhsPipeline().push_back(
      new OpCalculateScalarFieldGradient<SPACE_DIM>(domainField, gradUValsPtr));
  pipeline_energy->getOpDomainRhsPipeline().push_back(
      new OpElectricField(negGradUPtr, gradUValsPtr));

  commonDataPtr = boost::make_shared<DataAtIntegrationPts>(mField);

  pipeline_energy->getOpDomainRhsPipeline().push_back(new OpTotalEnergy(
      domainField, gradUValsPtr, permBlockSetsPtr, intrnlDomainBlockSetsPtr,
      commonDataPtr, petscVecEnergy));

  pipeline_energy->loopFiniteElements();
  CHKERR VecAssemblyBegin(petscVecEnergy);
  CHKERR VecAssemblyEnd(petscVecEnergy);

  if (!mField.get_comm_rank()) {
    const double *array;

    CHKERR VecGetArrayRead(petscVecEnergy, &array);
    totalEnergy = array[ZERO];
    MOFEM_LOG_CHANNEL("SELF");
    MOFEM_LOG_C("SELF", Sev::inform, "Total Energy: %6.15f", totalEnergy);
    CHKERR VecRestoreArrayRead(petscVecEnergy, &array);
  }

  MoFEMFunctionReturn(0);
}
//! [Get Total Energy]

//! [Get Charges]
MoFEMErrorCode Electrostatics::getElectrodeCharge() {
  MoFEMFunctionBegin;
  auto op_loop_side = new OpLoopSide<SideEle>(
      mField, simpleInterface->getDomainFEName(), SPACE_DIM);

  CHKERR AddHOOps<SPACE_DIM, SPACE_DIM, SPACE_DIM>::add(
      op_loop_side->getOpPtrVector(), {H1});

  auto gradPtrCharge = boost::make_shared<MatrixDouble>();
  auto negGradUPtrCharge = boost::make_shared<MatrixDouble>();

  op_loop_side->getOpPtrVector().push_back(
      new OpCalculateScalarFieldGradient<SPACE_DIM>(domainField,
                                                    gradPtrCharge));

  op_loop_side->getOpPtrVector().push_back(
      new OpElectricField(negGradUPtrCharge, gradPtrCharge));
  auto dJump = boost::make_shared<MatrixDouble>();
  op_loop_side->getOpPtrVector().push_back(new OpElectricDispJump<SPACE_DIM>(
      domainField, negGradUPtrCharge, dJump, commonDataPtr, permBlockSetsPtr));

  electrode_rhs_fe = boost::shared_ptr<ForcesAndSourcesCore>(
      new intElementForcesAndSourcesCore(mField));
  electrode_rhs_fe->getRuleHook = [](int, int, int p) { return p; };

  // push all the operators in on the side to the electrode_rhs_fe
  electrode_rhs_fe->getOpPtrVector().push_back(op_loop_side);

  electrode_rhs_fe->getOpPtrVector().push_back(new OpElectrodeCharge<SPACE_DIM>(
      domainField, dJump, petscVec, electrodeBlockSetsPtr));
  CHKERR VecZeroEntries(petscVec);
  CHKERR DMoFEMLoopFiniteElementsUpAndLowRank(
      simpleInterface->getDM(), "ELECTRODE", electrode_rhs_fe,
      mField.get_comm_rank(), mField.get_comm_rank());
  CHKERR VecAssemblyBegin(petscVec);
  CHKERR VecAssemblyEnd(petscVec);

  if (!mField.get_comm_rank()) {
    const double *array;

    CHKERR(VecGetArrayRead(petscVec, &array));
    double ALPHA = array[0]; // Use explicit index instead of ZERO
    double BETA = array[1];  // Use explicit index instead of ONE
    MOFEM_LOG_CHANNEL("SELF");
    MOFEM_LOG_C("SELF", Sev::inform,
                "CHARGE_ELEC_1: %6.15f , CHARGE_ELEC_2: %6.15f", ALPHA, BETA);

    CHKERR(VecRestoreArrayRead(petscVec, &array));
  }

  if (atom_test && !mField.get_comm_rank()) {
    double cal_charge_elec1;
    double cal_charge_elec2;
    double cal_total_energy;
    const double *c_ptr, *e_ptr;

    // Get a pointer to the PETSc vector data
    CHKERR(VecGetArrayRead(petscVec, &c_ptr));
    CHKERR(VecGetArrayRead(petscVecEnergy, &e_ptr));

    // Expected charges at the electrodes
    double charge_elec1 = 50.0;
    double charge_elec2 = -50.0;
    // Expected total energy of the system
    double total_energy = 500.0;

    switch (atom_test) {
    case 1:                        // 2D & 3D test
      cal_charge_elec1 = c_ptr[0]; // Read the  charge at the first electrode
      cal_charge_elec2 = c_ptr[1]; // Read the charge at the second electrode
      cal_total_energy = e_ptr[0]; // Read the total energy of the system
      break;
    default:
      SETERRQ1(PETSC_COMM_SELF, MOFEM_INVALID_DATA,
               "atom test %d does not exist", atom_test);
    }

    // Validate the results
    if (std::abs(charge_elec1 - cal_charge_elec1) > 1e-10 ||
        std::abs(charge_elec2 - cal_charge_elec2) > 1e-10 ||
        std::abs(total_energy - cal_total_energy) > 1e-10) {
      SETERRQ1(
          PETSC_COMM_SELF, MOFEM_ATOM_TEST_INVALID,
          "atom test %d failed! Calculated values do not match expected values",
          atom_test);
    }

    CHKERR(VecRestoreArrayRead(petscVec,
                               &c_ptr)); // Restore the PETSc vector array
    CHKERR(VecRestoreArrayRead(petscVecEnergy, &e_ptr));
  }

  MoFEMFunctionReturn(0);
}
//! [Get Charges]

//! [Run program]
MoFEMErrorCode Electrostatics::runProgram() {
  MoFEMFunctionBegin;

  CHKERR readMesh();
  CHKERR setupProblem();
  CHKERR boundaryCondition();
  CHKERR setIntegrationRules();
  CHKERR assembleSystem();
  CHKERR solveSystem();
  CHKERR outputResults();
  CHKERR getTotalEnergy();
  CHKERR getElectrodeCharge();
  MoFEMFunctionReturn(0);
}
//! [Run program]

//! [Main]
int main(int argc, char *argv[]) {
  // Initialisation of MoFEM/PETSc and MOAB data structures
  const char param_file[] = "param_file.petsc";
  MoFEM::Core::Initialize(&argc, &argv, param_file, help);

  // Error handling
  try {
    // Register MoFEM discrete manager in PETSc
    DMType dm_name = "DMMOFEM";
    CHKERR DMRegister_MoFEM(dm_name);

    // Create MOAB instance
    moab::Core mb_instance;              // mesh database
    moab::Interface &moab = mb_instance; // mesh database interface

    // Create MoFEM instance
    MoFEM::Core core(moab);           // finite element database
    MoFEM::Interface &m_field = core; // finite element interface

    // Run the main analysis
    Electrostatics Electrostatics_problem(m_field);
    CHKERR Electrostatics_problem.runProgram();
  }
  CATCH_ERRORS;

  // Finish work: cleaning memory, getting statistics, etc.
  MoFEM::Core::Finalize();

  return 0;
}
//! [Main]
