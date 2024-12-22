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

  int oRder = 2;      // default order
  int geom_order = 1; // default gemoetric order
  MoFEM::Interface &mField;
  Simple *simpleInterface;
  std::string domainField;
  boost::shared_ptr<std::map<int, BlockData>> permBlockSetsPtr;
  boost::shared_ptr<std::map<int, BlockData>> intBlockSetsPtr;
  boost::shared_ptr<std::map<int, BlockData>> electrodeBlockSetsPtr;
  boost::shared_ptr<DataAtIntegrationPts> commonDataPtr;

  boost::shared_ptr<ForcesAndSourcesCore> interFaceRhsFe;
  boost::shared_ptr<ForcesAndSourcesCore> electrodeRhsFe;

  double aLpha = 0.0; // declaration for total charge on first electrode
  double bEta = 0.0;  // declaration for total charge on second electrode
  SmartPetscObj<Vec> petscVec;       // petsc vector for the charge solution
  SmartPetscObj<Vec> petscVecEnergy; // petsc vector for the energy solution
  PetscBool out_skin = PETSC_FALSE;  //
  PetscBool is_partitioned = PETSC_FALSE;
  enum VecElements { ZERO = 0, ONE = 1, LAST_ELEMENT };
  int atomTest = 0;
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

  CHKERR PetscOptionsGetInt(PETSC_NULL, "", "-order", &oRder, PETSC_NULL);

  CHKERR PetscOptionsGetInt(PETSC_NULL, "", "-geom_order", &geom_order,
                            PETSC_NULL);

  CHKERR PetscOptionsGetInt(PETSC_NULL, "", "-atom_test", &atomTest,
                            PETSC_NULL);
  CHKERR simpleInterface->setFieldOrder(domainField, oRder);
  CHKERR simpleInterface->addDataField("GEOMETRY", H1, base, SPACE_DIM);
  CHKERR simpleInterface->setFieldOrder("GEOMETRY", geom_order);

  auto project_ho_geometry = [&]() {
    Projection10NodeCoordsOnField ent_method(mField, "GEOMETRY");
    return mField.loop_dofs("GEOMETRY", ent_method);
  };

  CHKERR simpleInterface->setUp();
  CHKERR project_ho_geometry();

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

      block_data.iD = id;                       // id-> block ID
      block_data.chargeDensity = attributes[0]; // block charge attribute
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

  // sync entities
  CHKERR mField.getInterface<CommInterface>()->synchroniseEntities(
      mat_electr_ents);
  CHKERR mField.getInterface<CommInterface>()->synchroniseEntities(
      int_electr_ents);
  CHKERR mField.getInterface<CommInterface>()->synchroniseEntities(
      electrode_ents);

  CHKERR PetscOptionsGetBool(PETSC_NULL, "", "-out_skin", &out_skin,
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
  CHKERR mField.modify_finite_element_add_field_data("INTERFACE", "GEOMETRY");

  CHKERR mField.add_ents_to_finite_element_by_dim(int_electr_ents,
                                                  SPACE_DIM - 1, "INTERFACE");
  // add the electrode entities to the field
  CHKERR mField.add_finite_element("ELECTRODE");
  CHKERR mField.modify_finite_element_add_field_row("ELECTRODE", domainField);
  CHKERR mField.modify_finite_element_add_field_col("ELECTRODE", domainField);
  CHKERR mField.modify_finite_element_add_field_data("ELECTRODE", domainField);
  CHKERR mField.modify_finite_element_add_field_data("ELECTRODE", "GEOMETRY");

  CHKERR mField.add_ents_to_finite_element_by_dim(electrode_ents, SPACE_DIM - 1,
                                                  "ELECTRODE");

  // sync field entities
  mField.getInterface<CommInterface>()->synchroniseFieldEntities(domainField);
  mField.getInterface<CommInterface>()->synchroniseFieldEntities("GEOMETRY");
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

  auto rule_lhs = [this](int, int, int p) -> int { return 2 * p + geom_order; };
  auto rule_rhs = [this](int, int, int p) -> int { return 2 * p + geom_order; };

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

  auto add_domain_base_ops = [&](auto &pipeline) {
    CHKERR AddHOOps<SPACE_DIM, SPACE_DIM, SPACE_DIM>::add(pipeline, {H1},
                                                          "GEOMETRY");

    pipeline.push_back(
        new OpBlockPermittivity(commonDataPtr, permBlockSetsPtr, domainField));
  };

  add_domain_base_ops(pipeline_mng->getOpDomainLhsPipeline());
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

      add_domain_base_ops(pipeline_mng->getOpDomainRhsPipeline());

      auto grad_u_ptr = boost::make_shared<MatrixDouble>();
      pipeline_mng->getOpDomainRhsPipeline().push_back(
          new OpCalculateScalarFieldGradient<SPACE_DIM>(domainField,
                                                        grad_u_ptr));
      auto minus_epsilon = [&](double, double, double) constexpr {
        return -commonDataPtr->blockPermittivity;
      };
      pipeline_mng->getOpDomainRhsPipeline().push_back(
          new OpInternal(domainField, grad_u_ptr, minus_epsilon));
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

  interFaceRhsFe = boost::shared_ptr<ForcesAndSourcesCore>(
      new IntElementForcesAndSourcesCore(mField));
  interFaceRhsFe->getRuleHook = [this](int, int, int p) {
    return 2 * p + geom_order;
  };

  {

    CHKERR AddHOOps<SPACE_DIM - 1, SPACE_DIM, SPACE_DIM>::add(
        interFaceRhsFe->getOpPtrVector(), {NOSPACE}, "GEOMETRY");

    interFaceRhsFe->getOpPtrVector().push_back(
        new OpBlockChargeDensity(commonDataPtr, intBlockSetsPtr, domainField));

    auto sIgma = [&](const double, const double, const double) {
      return commonDataPtr->blockChrgDens;
    };

    interFaceRhsFe->getOpPtrVector().push_back(
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
  DM dm;
  CHKERR simpleInterface->getDM(&dm);

  CHKERR DMMoFEMKSPSetComputeRHS(dm, "INTERFACE", interFaceRhsFe, null, null);

  CHKERR KSPSetFromOptions(ksp_solver);

  // Create RHS and solution vectors
  auto F = createDMVector(dm);
  auto D = vectorDuplicate(F);
  // Solve the system
  CHKERR KSPSetUp(ksp_solver);
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

  // lamda function to calculate electric field
  auto calculate_e_field = [&](auto &pipeline) {
    auto u_ptr = boost::make_shared<VectorDouble>();
    auto x_ptr = boost::make_shared<MatrixDouble>();
    auto grad_u_ptr = boost::make_shared<MatrixDouble>();
    auto e_field_ptr = boost::make_shared<MatrixDouble>();
    // add higher order operator
    CHKERR AddHOOps<SPACE_DIM, SPACE_DIM, SPACE_DIM>::add(pipeline, {H1},
                                                          "GEOMETRY");
    // calculate field values
    pipeline.push_back(new OpCalculateScalarFieldValues(domainField, u_ptr));
    pipeline.push_back(
        new OpCalculateVectorFieldValues<SPACE_DIM>("GEOMETRY", x_ptr));

    // calculate gradient
    pipeline.push_back(
        new OpCalculateScalarFieldGradient<SPACE_DIM>(domainField, grad_u_ptr));
    // calculate electric field
    pipeline.push_back(new OpElectricField(e_field_ptr, grad_u_ptr));
    return boost::make_tuple(u_ptr, e_field_ptr, x_ptr);
  };

  auto [u_ptr, e_field_ptr, x_ptr] =
      calculate_e_field(post_proc_fe->getOpPtrVector());

  auto e_field_times_perm_ptr = boost::make_shared<MatrixDouble>();
  auto energy_density_ptr = boost::make_shared<VectorDouble>();

  post_proc_fe->getOpPtrVector().push_back(
      new OpGradTimesPerm(domainField, e_field_ptr, e_field_times_perm_ptr,
                          permBlockSetsPtr, commonDataPtr));
  post_proc_fe->getOpPtrVector().push_back(
      new OpEnergyDensity(domainField, e_field_ptr, energy_density_ptr,
                          permBlockSetsPtr, commonDataPtr));

  using OpPPMap = OpPostProcMapInMoab<SPACE_DIM, SPACE_DIM>;
  post_proc_fe->getOpPtrVector().push_back(new OpPPMap(
      post_proc_fe->getPostProcMesh(), post_proc_fe->getMapGaussPts(),

      OpPPMap::DataMapVec{{"POTENTIAL", u_ptr},
                          {"ENERGY_DENSITY", energy_density_ptr}},
      OpPPMap::DataMapMat{
          {"GEOMETRY", x_ptr},
          {"ELECTRIC_FIELD", e_field_ptr},
          {"ELECTRIC_DISPLACEMENT", e_field_times_perm_ptr},
      },
      OpPPMap::DataMapMat{},

      OpPPMap::DataMapMat{})

  );

  pipeline_mng->getDomainRhsFE() = post_proc_fe;
  CHKERR pipeline_mng->loopFiniteElements();
  CHKERR post_proc_fe->writeFile("out.h5m");

  if (out_skin && SPACE_DIM == 3) {

    auto post_proc_skin = boost::make_shared<PostProcFaceEle>(mField);
    auto op_loop_skin = new OpLoopSide<SideEle>(
        mField, simpleInterface->getDomainFEName(), SPACE_DIM);

    auto [u_ptr, e_field_ptr, x_ptr] =
        calculate_e_field(op_loop_skin->getOpPtrVector());

    op_loop_skin->getOpPtrVector().push_back(
        new OpGradTimesPerm(domainField, e_field_ptr, e_field_times_perm_ptr,
                            permBlockSetsPtr, commonDataPtr));
    op_loop_skin->getOpPtrVector().push_back(
        new OpEnergyDensity(domainField, e_field_ptr, energy_density_ptr,
                            permBlockSetsPtr, commonDataPtr));

    // push op to boundary element
    post_proc_skin->getOpPtrVector().push_back(op_loop_skin);

    post_proc_skin->getOpPtrVector().push_back(new OpPPMap(
        post_proc_skin->getPostProcMesh(), post_proc_skin->getMapGaussPts(),
        OpPPMap::DataMapVec{{"POTENTIAL", u_ptr},
                            {"ENERGY_DENSITY", energy_density_ptr}},
        OpPPMap::DataMapMat{{"ELECTRIC_FIELD", e_field_ptr},
                            {"GEOMETRY", x_ptr},
                            {"ELECTRIC_DISPLACEMENT", e_field_times_perm_ptr}},
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
  auto pip_energy = mField.getInterface<PipelineManager>();
  pip_energy->getDomainRhsFE().reset();
  pip_energy->getDomainLhsFE().reset();
  pip_energy->getBoundaryLhsFE().reset();
  pip_energy->getBoundaryRhsFE().reset();
  pip_energy->getOpDomainRhsPipeline().clear();
  pip_energy->getOpDomainLhsPipeline().clear();

  // gets the map of the internal domain entity range to get the total energy

  boost::shared_ptr<std::map<int, BlockData>> intrnlDomnBlckSetPtr =
      boost::make_shared<std::map<int, BlockData>>();
  Range internal_domain; // range of entities marked the internal domain
  for (_IT_CUBITMESHSETS_BY_SET_TYPE_FOR_LOOP_(mField, BLOCKSET, bit)) {
    if (bit->getName().compare(0, 10, "DOMAIN_INT") == 0) {
      const int id = bit->getMeshsetId();
      auto &block_data = (*intrnlDomnBlckSetPtr)[id];

      CHKERR mField.get_moab().get_entities_by_dimension(
          bit->getMeshset(), SPACE_DIM, block_data.internalDomainEnts, true);
      internal_domain.merge(block_data.internalDomainEnts);
      block_data.iD = id;
    }
  }
  CHKERR mField.getInterface<CommInterface>()->synchroniseEntities(
      internal_domain);

  CHKERR AddHOOps<SPACE_DIM, SPACE_DIM, SPACE_DIM>::add(
      pip_energy->getOpDomainRhsPipeline(), {H1}, "GEOMETRY");

  auto grad_u_ptr = boost::make_shared<MatrixDouble>();
  auto e_field_ptr = boost::make_shared<MatrixDouble>();

  pip_energy->getOpDomainRhsPipeline().push_back(
      new OpCalculateScalarFieldGradient<SPACE_DIM>(domainField, grad_u_ptr));
  pip_energy->getOpDomainRhsPipeline().push_back(
      new OpElectricField(e_field_ptr, grad_u_ptr));

  commonDataPtr = boost::make_shared<DataAtIntegrationPts>(mField);

  pip_energy->getOpDomainRhsPipeline().push_back(
      new OpTotalEnergy(domainField, grad_u_ptr, permBlockSetsPtr,
                        intrnlDomnBlckSetPtr, commonDataPtr, petscVecEnergy));

  pip_energy->loopFiniteElements();
  CHKERR VecAssemblyBegin(petscVecEnergy);
  CHKERR VecAssemblyEnd(petscVecEnergy);

  double total_energy = 0.0; // declaration for total energy
  if (!mField.get_comm_rank()) {
    const double *array;

    CHKERR VecGetArrayRead(petscVecEnergy, &array);
    total_energy = array[ZERO];
    MOFEM_LOG_CHANNEL("SELF");
    MOFEM_LOG_C("SELF", Sev::inform, "Total Energy: %6.15f", total_energy);
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
      op_loop_side->getOpPtrVector(), {H1}, "GEOMETRY");

  auto grad_u_ptr_charge = boost::make_shared<MatrixDouble>();
  auto e_ptr_charge = boost::make_shared<MatrixDouble>();

  op_loop_side->getOpPtrVector().push_back(
      new OpCalculateScalarFieldGradient<SPACE_DIM>(domainField,
                                                    grad_u_ptr_charge));

  op_loop_side->getOpPtrVector().push_back(
      new OpElectricField(e_ptr_charge, grad_u_ptr_charge));
  auto d_jump = boost::make_shared<MatrixDouble>();
  op_loop_side->getOpPtrVector().push_back(new OpElectricDispJump<SPACE_DIM>(
      domainField, e_ptr_charge, d_jump, commonDataPtr, permBlockSetsPtr));

  electrodeRhsFe = boost::shared_ptr<ForcesAndSourcesCore>(
      new IntElementForcesAndSourcesCore(mField));
  electrodeRhsFe->getRuleHook = [this](int, int, int p) {
    return 2 * p + geom_order;
  };

  // push all the operators in on the side to the electrodeRhsFe
  electrodeRhsFe->getOpPtrVector().push_back(op_loop_side);

  electrodeRhsFe->getOpPtrVector().push_back(new OpElectrodeCharge<SPACE_DIM>(
      domainField, d_jump, petscVec, electrodeBlockSetsPtr));
  CHKERR VecZeroEntries(petscVec);
  CHKERR DMoFEMLoopFiniteElementsUpAndLowRank(simpleInterface->getDM(),
                                              "ELECTRODE", electrodeRhsFe, 0,
                                              mField.get_comm_size());
  CHKERR VecAssemblyBegin(petscVec);
  CHKERR VecAssemblyEnd(petscVec);

  if (!mField.get_comm_rank()) {
    const double *array;

    CHKERR(VecGetArrayRead(petscVec, &array));
    double aLpha = array[0]; // Use explicit index instead of ZERO
    double bEta = array[1];  // Use explicit index instead of ONE
    MOFEM_LOG_CHANNEL("SELF");
    MOFEM_LOG_C("SELF", Sev::inform,
                "CHARGE_ELEC_1: %6.15f , CHARGE_ELEC_2: %6.15f", aLpha, bEta);

    CHKERR(VecRestoreArrayRead(petscVec, &array));
  }
  if (atomTest && !mField.get_comm_rank()) {
    double cal_charge_elec1;
    double cal_charge_elec2;
    double cal_total_energy;
    const double *c_ptr, *te_ptr;

    // Get a pointer to the PETSc vector data
    CHKERR(VecGetArrayRead(petscVec, &c_ptr));
    CHKERR(VecGetArrayRead(petscVecEnergy, &te_ptr));

    // Expected charges at the electrodes
    double ref_charge_elec1;
    double ref_charge_elec2;
    // Expected total energy of the system
    double ref_tot_energy;
    double tol;
    cal_charge_elec1 = c_ptr[0];  // Read charge at the first electrode
    cal_charge_elec2 = c_ptr[1];  // Read charge at the second electrode
    cal_total_energy = te_ptr[0]; // Read total energy of the system
    if (std::isnan(cal_charge_elec1) || std::isnan(cal_charge_elec2) ||
        std::isnan(cal_total_energy)) {
      SETERRQ(PETSC_COMM_SELF, MOFEM_ATOM_TEST_INVALID,
              "Atom test failed! NaN detected in calculated values.");
    }
    switch (atomTest) {
    case 1: // 2D & 3D test
      // Expected charges at the electrodes
      ref_charge_elec1 = 50.0;
      ref_charge_elec2 = -50.0;
      // Expected total energy of the system
      ref_tot_energy = 500.0;
      tol = 1e-10;
      break;
    case 2: // wavy 3D test
      ref_charge_elec1 = 10.00968352472943;
      ref_charge_elec2 = 0.0; // no electrode
      ref_tot_energy = 50.5978;
      tol = 1e-4;
      break;
    default:
      SETERRQ1(PETSC_COMM_SELF, MOFEM_INVALID_DATA,
               "atom test %d does not exist", atomTest);
    }

    // Validate the results
    if (std::abs(ref_charge_elec1 - cal_charge_elec1) > tol ||
        std::abs(ref_charge_elec2 - cal_charge_elec2) > tol ||
        std::abs(ref_tot_energy - cal_total_energy) > tol) {
      SETERRQ1(
          PETSC_COMM_SELF, MOFEM_ATOM_TEST_INVALID,
          "atom test %d failed! Calculated values do not match expected values",
          atomTest);
    }

    CHKERR(VecRestoreArrayRead(petscVec,
                               &c_ptr)); // Restore the PETSc vector array
    CHKERR(VecRestoreArrayRead(petscVecEnergy, &te_ptr));
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
