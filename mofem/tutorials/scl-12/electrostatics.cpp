/**
 * @file Electrostatics.cpp
 * \example Electrostatics.cpp
 *  */

#ifndef EXECUTABLE_DIMENSION
  #define EXECUTABLE_DIMENSION 3
#endif

#include <electrostatics.hpp>

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
  MoFEMErrorCode assembleSystem();
  MoFEMErrorCode setIntegrationRules();
  MoFEMErrorCode solveSystem();
  MoFEMErrorCode outputResults();
  MoFEMErrorCode getTotalEnergy();
  MoFEMErrorCode getElectrodeCharge();

  MoFEM::Interface &mField;
  boost::shared_ptr<std::map<int, BlockData>> perm_block_sets_ptr;
  boost::shared_ptr<std::map<int, BlockData>> int_block_sets_ptr;
  boost::shared_ptr<std::map<int, BlockData>> electrode_block_sets_ptr;
  boost::shared_ptr<std::map<int, BlockData>> internal_domain_block_sets_ptr;

  Simple *simpleInterface;
  boost::shared_ptr<ForcesAndSourcesCore> interface_rhs_fe;
  boost::shared_ptr<ForcesAndSourcesCore> electrode_rhs_fe;
  boost::shared_ptr<DataAtIntegrationPts> common_data_ptr;
  std::string domainField;
  int oRder;
  double ALPHA = 0.0;       // declaration for total charge on first electrode
  double BETA = 0.0;        // declaration for total charge on second electrode
  double totalEnergy = 0.0; // declaration for total energy
  SmartPetscObj<Vec> petscVec;        // petsc vector for the charge solution
  SmartPetscObj<Vec> petscVec_energy; // petsc vector for the energy solution
  PetscBool out_skin = PETSC_FALSE;   //
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
      true; // creates lower-dimensional elements to make integration possible
            // across partitioned entities
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
      CHK_THROW_MESSAGE(MOFEM_NOT_FOUND, "Element type not handled");
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
  common_data_ptr = boost::make_shared<DataAtIntegrationPts>(mField);

  // gets the map of the permittivity attributes and the block sets
  perm_block_sets_ptr = boost::make_shared<std::map<int, BlockData>>();
  Range electrIcs; // range of entities with the permittivity
  for (_IT_CUBITMESHSETS_BY_SET_TYPE_FOR_LOOP_(mField, BLOCKSET, bit)) {
    if (bit->getName().compare(0, 12, "MAT_ELECTRIC") == 0) {
      const int id = bit->getMeshsetId();
      auto &block_data = (*perm_block_sets_ptr)[id];

      CHKERR mField.get_moab().get_entities_by_dimension(
          bit->getMeshset(), SPACE_DIM, block_data.domainEnts, true);
      electrIcs.merge(block_data.domainEnts);

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
  int_block_sets_ptr = boost::make_shared<std::map<int, BlockData>>();
  Range interfIcs; // range of entities with the charge
  for (_IT_CUBITMESHSETS_BY_SET_TYPE_FOR_LOOP_(mField, BLOCKSET, bit)) {
    if (bit->getName().compare(0, 12, "INT_ELECTRIC") == 0) {
      const int id = bit->getMeshsetId();
      auto &block_data = (*int_block_sets_ptr)[id];

      CHKERR mField.get_moab().get_entities_by_dimension(
          bit->getMeshset(), SPACE_DIM - 1, block_data.interfaceEnts, true);
      interfIcs.merge(block_data.interfaceEnts);

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
  electrode_block_sets_ptr = boost::make_shared<std::map<int, BlockData>>();
  Range ElectrodIcs; // range of entities with the electrode
  for (_IT_CUBITMESHSETS_BY_SET_TYPE_FOR_LOOP_(mField, BLOCKSET, bit)) {
    if (bit->getName().compare(0, 9, "ELECTRODE") == 0) {
      const int id = bit->getMeshsetId();
      auto &block_data = (*electrode_block_sets_ptr)[id];

      CHKERR mField.get_moab().get_entities_by_dimension(
          bit->getMeshset(), SPACE_DIM - 1, block_data.electrodeEnts, true);
      ElectrodIcs.merge(block_data.electrodeEnts);
      block_data.iD = id;
      auto print_range_on_procs = [&](const std::string &name, int meshsetId,
                                      const Range &range) {
        MOFEM_LOG("SYNC", Sev::inform)
            << name << " in meshID: " << id << " with range "
            << block_data.electrodeEnts << " on proc ["
            << mField.get_comm_rank() << "] \n";
        MOFEM_LOG_SYNCHRONISE(mField.get_comm());
      };
      print_range_on_procs(bit->getName(), id, ElectrodIcs);
    }
  }

  // gets the map of the internal domain entity range to get the total energy
  Range internal_domain; // range of entities marked the internal domain
  internal_domain_block_sets_ptr =
      boost::make_shared<std::map<int, BlockData>>();
  for (_IT_CUBITMESHSETS_BY_SET_TYPE_FOR_LOOP_(mField, BLOCKSET, bit)) {
    if (bit->getName().compare(0, 10, "DOMAIN_INT") == 0) {
      const int id = bit->getMeshsetId();
      auto &block_data = (*internal_domain_block_sets_ptr)[id];

      CHKERR mField.get_moab().get_entities_by_dimension(
          bit->getMeshset(), SPACE_DIM, block_data.internalDomainEnts, true);
      internal_domain.merge(block_data.internalDomainEnts);
      block_data.iD = id;
      auto print_range_on_procs = [&](const std::string &name, int meshsetId,
                                      const Range &range) {
        MOFEM_LOG("SYNC", Sev::inform)
            << name << " in meshID: " << id << " with range " << internal_domain
            << " on proc [" << mField.get_comm_rank() << "] \n";
        MOFEM_LOG_SYNCHRONISE(mField.get_comm());
      };
      print_range_on_procs(bit->getName(), id, internal_domain);
    }
  }
  // sync entities
  CHKERR mField.getInterface<CommInterface>()->synchroniseEntities(electrIcs);
  CHKERR mField.getInterface<CommInterface>()->synchroniseEntities(interfIcs);
  CHKERR mField.getInterface<CommInterface>()->synchroniseEntities(ElectrodIcs);
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
  CHKERR skinner.find_skin(0, electrIcs, false, skin_tris);
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
  CHKERR mField.add_ents_to_finite_element_by_dim(interfIcs, SPACE_DIM - 1,
                                                  "INTERFACE");
  // add the electrode entities to the field
  CHKERR mField.add_finite_element("ELECTRODE");
  CHKERR mField.modify_finite_element_add_field_row("ELECTRODE", domainField);
  CHKERR mField.modify_finite_element_add_field_col("ELECTRODE", domainField);
  CHKERR mField.modify_finite_element_add_field_data("ELECTRODE", domainField);
  CHKERR mField.add_ents_to_finite_element_by_dim(ElectrodIcs, SPACE_DIM - 1,
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
  petscVec_energy =
      createVectorMPI(mField.get_comm(), local_size, LAST_ELEMENT);
  // LAST_ELEMENT);
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

//! [Assemble system]
MoFEMErrorCode Electrostatics::assembleSystem() {
  MoFEMFunctionBegin;

  auto pipeline_mng = mField.getInterface<PipelineManager>();
  pipeline_mng->getDomainRhsFE().reset();
  pipeline_mng->getDomainLhsFE().reset();
  pipeline_mng->getBoundaryLhsFE().reset();
  pipeline_mng->getBoundaryRhsFE().reset();

  common_data_ptr = boost::make_shared<DataAtIntegrationPts>(mField);
  auto add_domain_lhs_ops = [&](auto &pipeline) {
    pipeline.push_back(new OpBlockPermittivity(
        common_data_ptr, perm_block_sets_ptr, domainField));
  };

  add_domain_lhs_ops(pipeline_mng->getOpDomainLhsPipeline());
  auto epsilon = [&](const double, const double, const double) {
    return common_data_ptr->blockPermittivity;
  };

  { // Push operators to the Pipeline that is responsible for calculating LHS
    CHKERR AddHOOps<SPACE_DIM, SPACE_DIM, SPACE_DIM>::add(
        pipeline_mng->getOpDomainLhsPipeline(), {H1});
    pipeline_mng->getOpDomainLhsPipeline().push_back(
        new OpDomainLhsMatrixK(domainField, domainField, epsilon));
  }

  { // Push operators to the Pipeline that is responsible for calculating LHS
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

      auto grad_u_vals_ptr = boost::make_shared<MatrixDouble>();
      pipeline_mng->getOpDomainRhsPipeline().push_back(
          new OpCalculateScalarFieldGradient<SPACE_DIM>(domainField,
                                                        grad_u_vals_ptr));
      add_domain_lhs_ops(pipeline_mng->getOpDomainRhsPipeline());
      auto minus_epsilon = [&](double, double, double) constexpr {
        return -common_data_ptr->blockPermittivity;
      };
      pipeline_mng->getOpDomainRhsPipeline().push_back(
          new OpInternal(domainField, grad_u_vals_ptr, minus_epsilon));
    };

    CHKERR AddHOOps<SPACE_DIM, SPACE_DIM, SPACE_DIM>::add(
        pipeline_mng->getOpDomainRhsPipeline(), {H1});
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

  {
    interface_rhs_fe->getOpPtrVector().push_back(new OpBlockChargeDensity(
        common_data_ptr, int_block_sets_ptr, domainField));

    auto sIgma = [&](const double, const double, const double) {
      return common_data_ptr->blockChrgDens;
    };

    interface_rhs_fe->getOpPtrVector().push_back(
        new OpInterfaceRhsVectorF(domainField, sIgma));
  }

  MoFEMFunctionReturn(0);
}
//! [Assemble system]

//! [Set integration rules]
MoFEMErrorCode Electrostatics::setIntegrationRules() {
  MoFEMFunctionBegin;

  auto rule_lhs = [](int, int, int p) -> int { return 2 * p + 1; };
  auto rule_rhs = [](int, int, int p) -> int { return p + 1; };

  auto pipeline_mng = mField.getInterface<PipelineManager>();
  CHKERR pipeline_mng->setDomainLhsIntegrationRule(rule_lhs);
  CHKERR pipeline_mng->setDomainRhsIntegrationRule(rule_rhs);

  MoFEMFunctionReturn(0);
}
//! [Set integration rules]

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

  auto u_ptr = boost::make_shared<VectorDouble>();
  auto grad_u_ptr = boost::make_shared<MatrixDouble>();
  post_proc_fe->getOpPtrVector().push_back(
      new OpCalculateScalarFieldValues(domainField, u_ptr));

  post_proc_fe->getOpPtrVector().push_back(
      new OpCalculateScalarFieldGradient<SPACE_DIM>(domainField, grad_u_ptr));

  using OpPPMap = OpPostProcMapInMoab<SPACE_DIM, SPACE_DIM>;

  auto neg_grad_u_ptr = boost::make_shared<MatrixDouble>();
  auto neg_grad_u_ptr_times_perm = boost::make_shared<MatrixDouble>();
  auto energy_density_ptr = boost::make_shared<MatrixDouble>();
  post_proc_fe->getOpPtrVector().push_back(
      new OpElectricField(neg_grad_u_ptr, grad_u_ptr));

  post_proc_fe->getOpPtrVector().push_back(new OpGradTimesperm(
      domainField, neg_grad_u_ptr, neg_grad_u_ptr_times_perm,
      perm_block_sets_ptr, common_data_ptr));
  post_proc_fe->getOpPtrVector().push_back(
      new OpEnergyDensity(domainField, neg_grad_u_ptr, energy_density_ptr,
                          perm_block_sets_ptr, common_data_ptr));

  post_proc_fe->getOpPtrVector().push_back(

      new OpPPMap(
          post_proc_fe->getPostProcMesh(), post_proc_fe->getMapGaussPts(),

          OpPPMap::DataMapVec{{"POTENTIAL", u_ptr}},

          OpPPMap::DataMapMat{{"ELECTRIC FIELD", neg_grad_u_ptr},
                              {"FLUX DENSITY", neg_grad_u_ptr_times_perm},
                              {"ENERGY DENSITY", energy_density_ptr}},
          OpPPMap::DataMapMat{},

          OpPPMap::DataMapMat{}

          )

  );

  pipeline_mng->getDomainRhsFE() = post_proc_fe;
  CHKERR pipeline_mng->loopFiniteElements();
  CHKERR post_proc_fe->writeFile("out_resultFD.h5m");

  if (out_skin) {
    auto post_proc_skin = boost::make_shared<PostProcFaceEle>(mField);
    auto op_loop_side = new OpLoopSide<SideEle>(
        mField, simpleInterface->getDomainFEName(), SPACE_DIM);

    auto jac_ptr = boost::make_shared<MatrixDouble>();
    auto inv_jac_ptr = boost::make_shared<MatrixDouble>();
    auto det_ptr = boost::make_shared<VectorDouble>();

    auto u_ptr_skin = boost::make_shared<VectorDouble>();
    auto grad_u_ptr_skin = boost::make_shared<MatrixDouble>();
    auto neg_grad_u_ptr_times_perm1 = boost::make_shared<MatrixDouble>();
    auto energy_density_ptr = boost::make_shared<MatrixDouble>();
    op_loop_side->getOpPtrVector().push_back(
        new OpCalculateHOJac<SPACE_DIM>(jac_ptr));
    op_loop_side->getOpPtrVector().push_back(
        new OpInvertMatrix<SPACE_DIM>(jac_ptr, det_ptr, inv_jac_ptr));
    op_loop_side->getOpPtrVector().push_back(
        new OpSetHOInvJacToScalarBases<SPACE_DIM>(H1, inv_jac_ptr));
    op_loop_side->getOpPtrVector().push_back(
        new OpCalculateScalarFieldValues(domainField, u_ptr_skin));
    op_loop_side->getOpPtrVector().push_back(
        new OpCalculateScalarFieldGradient<SPACE_DIM>(domainField,
                                                      grad_u_ptr_skin));
    op_loop_side->getOpPtrVector().push_back(
        new OpElectricField(neg_grad_u_ptr, grad_u_ptr_skin));
    op_loop_side->getOpPtrVector().push_back(new OpGradTimesperm(
        domainField, grad_u_ptr_skin, neg_grad_u_ptr_times_perm1,
        perm_block_sets_ptr, common_data_ptr));
    op_loop_side->getOpPtrVector().push_back(
        new OpEnergyDensity(domainField, grad_u_ptr_skin, energy_density_ptr,
                            perm_block_sets_ptr, common_data_ptr));

    // push op to boundary element
    post_proc_skin->getOpPtrVector().push_back(op_loop_side);

    post_proc_skin->getOpPtrVector().push_back(new OpPPMap(
        post_proc_skin->getPostProcMesh(), post_proc_skin->getMapGaussPts(),
        OpPPMap::DataMapVec{{"POTENTIAL", u_ptr_skin}},
        {{"ELECTRIC_FIELD", neg_grad_u_ptr}}, {}, {}));
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
  auto pipeline_total_energy = mField.getInterface<PipelineManager>();
  pipeline_total_energy->getDomainRhsFE().reset();
  pipeline_total_energy->getDomainLhsFE().reset();
  pipeline_total_energy->getBoundaryLhsFE().reset();
  pipeline_total_energy->getBoundaryRhsFE().reset();
  pipeline_total_energy->getOpDomainRhsPipeline().clear();
  pipeline_total_energy->getOpDomainLhsPipeline().clear();

  CHKERR AddHOOps<SPACE_DIM, SPACE_DIM, SPACE_DIM>::add(
      pipeline_total_energy->getOpDomainRhsPipeline(), {H1});

  auto grad_u_vals_ptr = boost::make_shared<MatrixDouble>();
  auto neg_grad_u_ptr = boost::make_shared<MatrixDouble>();
  auto u_ptr = boost::make_shared<VectorDouble>();
  pipeline_total_energy->getOpDomainRhsPipeline().push_back(
      new OpCalculateScalarFieldValues(domainField, u_ptr));
  pipeline_total_energy->getOpDomainRhsPipeline().push_back(
      new OpCalculateScalarFieldGradient<SPACE_DIM>(domainField,
                                                    grad_u_vals_ptr));
  pipeline_total_energy->getOpDomainRhsPipeline().push_back(
      new OpElectricField(neg_grad_u_ptr, grad_u_vals_ptr));

  common_data_ptr = boost::make_shared<DataAtIntegrationPts>(mField);

  pipeline_total_energy->getOpDomainRhsPipeline().push_back(
      new OpTotalEnergyDensity(
          domainField, grad_u_vals_ptr, perm_block_sets_ptr,
          internal_domain_block_sets_ptr, common_data_ptr, petscVec_energy));

  pipeline_total_energy->loopFiniteElements();
  CHKERR VecAssemblyBegin(petscVec_energy);
  CHKERR VecAssemblyEnd(petscVec_energy);

  if (!mField.get_comm_rank()) {
    const double *array;

    CHKERR VecGetArrayRead(petscVec_energy, &array);
    totalEnergy = array[ZERO];
    MOFEM_LOG_CHANNEL("SELF");
    MOFEM_LOG_C("SELF", Sev::inform, "Total Energy: %6.15f", totalEnergy);
    CHKERR VecRestoreArrayRead(petscVec_energy, &array);
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

  auto grad_grad_ptr = boost::make_shared<MatrixDouble>();
  auto neg_grad_u_ptr_alfa = boost::make_shared<MatrixDouble>();

  op_loop_side->getOpPtrVector().push_back(
      new OpCalculateScalarFieldGradient<SPACE_DIM>(domainField,
                                                    grad_grad_ptr));

  op_loop_side->getOpPtrVector().push_back(
      new OpElectricField(neg_grad_u_ptr_alfa, grad_grad_ptr));
  auto d_jump = boost::make_shared<MatrixDouble>();
  op_loop_side->getOpPtrVector().push_back(
      new OpELectricJump<SPACE_DIM>(domainField, neg_grad_u_ptr_alfa, d_jump,
                                    common_data_ptr, perm_block_sets_ptr));

  electrode_rhs_fe = boost::shared_ptr<ForcesAndSourcesCore>(
      new intElementForcesAndSourcesCore(mField));
  electrode_rhs_fe->getRuleHook = [](int, int, int p) { return 3 * p; };

  // push all the operators in on the side  to the electrode_rhs_fe
  electrode_rhs_fe->getOpPtrVector().push_back(op_loop_side);

  electrode_rhs_fe->getOpPtrVector().push_back(new OpAlpha<SPACE_DIM>(
      domainField, d_jump, petscVec, electrode_block_sets_ptr));
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
    CHKERR(VecGetArrayRead(petscVec_energy, &e_ptr));

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
    CHKERR(VecRestoreArrayRead(petscVec_energy, &e_ptr));
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
