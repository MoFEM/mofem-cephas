/** \file MFrontInterface.cpp
 * \brief MFrontInterface
 *
 * MFrontInterface
 *
 */

#ifdef WITH_MGIS

  #include <MFrontInterface.hpp>

namespace MoFEM {

template struct MFrontInterface<TRIDIMENSIONAL>;
template struct MFrontInterface<AXISYMMETRICAL>;
template struct MFrontInterface<PLANESTRAIN>;

template <ModelHypothesis H>
MFrontInterface<H>::MFrontInterface(const MoFEM::Core &core)
    : cOre(const_cast<MoFEM::Core &>(core)) {

  isFiniteKinematics = true;
  saveGauss = PETSC_FALSE;
  saveVolume = PETSC_TRUE;
  testJacobian = PETSC_FALSE;
  randomFieldScale = 1.0;
  optionsPrefix = "mf_";
  monitorPtr = nullptr;
  fieldBase = NOBASE;

  // if (!LogManager::checkIfChannelExist("FieldEvaluatorWorld")) {
  //   auto core_log = logging::core::get();

  //   core_log->add_sink(LogManager::createSink(LogManager::getStrmWorld(),
  //                                             "FieldEvaluatorWorld"));
  //   core_log->add_sink(LogManager::createSink(LogManager::getStrmSync(),
  //                                             "FieldEvaluatorSync"));
  //   core_log->add_sink(LogManager::createSink(LogManager::getStrmSelf(),
  //                                             "FieldEvaluatorSelf"));

  //   LogManager::setLog("FieldEvaluatorWorld");
  //   LogManager::setLog("FieldEvaluatorSync");
  //   LogManager::setLog("FieldEvaluatorSelf");

  //   MOFEM_LOG_TAG("FieldEvaluatorWorld", "FieldEvaluator");
  //   MOFEM_LOG_TAG("FieldEvaluatorSync", "FieldEvaluator");
  //   MOFEM_LOG_TAG("FieldEvaluatorSelf", "FieldEvaluator");
  // }

  // MOFEM_LOG("FieldEvaluatorWorld", Sev::noisy) << "Field evaluator
  // intreface";
}

template <ModelHypothesis H>
MoFEMErrorCode
MFrontInterface<H>::query_interface(boost::typeindex::type_index type_index,
                                    UnknownInterface **iface) const {
  MoFEMFunctionBeginHot;
  *iface = const_cast<MFrontInterface<H> *>(this);
  MoFEMFunctionReturnHot(0);
}

template <ModelHypothesis H>
MoFEMErrorCode MFrontInterface<H>::getCommandLineParameters() {
  MoFEMFunctionBegin;
  // isQuasiStatic = PETSC_FALSE;
  // if (oRder == -1)
  //   oRder = 2;
  CHKERR PetscOptionsBegin(PETSC_COMM_WORLD, optionsPrefix.c_str(), "", "none");

  CHKERR PetscOptionsBool("-save_gauss", "save gauss pts (internal variables)",
                          "", saveGauss, &saveGauss, PETSC_NULL);
  CHKERR PetscOptionsBool("-save_volume", "save results on a volumetric mesh",
                          "", saveVolume, &saveVolume, PETSC_NULL);

  CHKERR PetscOptionsBool("-test_jacobian", "test Jacobian (LHS matrix)", "",
                          testJacobian, &testJacobian, PETSC_NULL);
  CHKERR PetscOptionsReal("-random_field_scale",
                          "scale for the finite difference jacobian", "",
                          randomFieldScale, &randomFieldScale, PETSC_NULL);

  if (saveGauss)
    moabGaussIntPtr = boost::shared_ptr<moab::Interface>(new moab::Core());

  MoFEM::Interface &m_field = cOre;
  commonDataPtr = boost::make_shared<CommonData>(m_field);

  commonDataPtr->setBlocks(DIM);
  commonDataPtr->createTags();

  commonDataPtr->mGradPtr = boost::make_shared<MatrixDouble>();
  commonDataPtr->mStressPtr = boost::make_shared<MatrixDouble>();
  commonDataPtr->mFullStrainPtr = boost::make_shared<MatrixDouble>();
  commonDataPtr->mFullStressPtr = boost::make_shared<MatrixDouble>();
  commonDataPtr->mDispPtr = boost::make_shared<MatrixDouble>();
  commonDataPtr->mPrevGradPtr = boost::make_shared<MatrixDouble>();
  commonDataPtr->mPrevStressPtr = boost::make_shared<MatrixDouble>();
  commonDataPtr->materialTangentPtr = boost::make_shared<MatrixDouble>();
  commonDataPtr->mFullTangentPtr = boost::make_shared<MatrixDouble>();
  commonDataPtr->internalVariablePtr = boost::make_shared<MatrixDouble>();

  if (commonDataPtr->setOfBlocksData.empty())
    MOFEM_LOG("WORLD", Sev::inform)
        << "No blocksets on the mesh have been provided for MFront (e.g. "
           "MFRONT_MAT_1)";

  auto check_lib_finite_strain = [&](const std::string &lib,
                                     const std::string &beh_name, bool &flag) {
    MoFEMFunctionBeginHot;

    ifstream f(lib.c_str());
    if (!f.good())
      MOFEM_LOG("WORLD", Sev::error)
          << "Problem with the behaviour path: " << lib;

    auto &lm = LibrariesManager::get();
    flag = bool(lm.getBehaviourType(lib, beh_name) == 2) &&
           (lm.getBehaviourKinematic(lib, beh_name) == 3);
    MoFEMFunctionReturnHot(0);
  };

  auto op = FiniteStrainBehaviourOptions{};
  op.stress_measure = FiniteStrainBehaviourOptions::PK1;
  op.tangent_operator = FiniteStrainBehaviourOptions::DPK1_DF;

  for (auto &block : commonDataPtr->setOfBlocksData) {
    const int &id = block.first;
    auto &lib_path = block.second.behaviourPath;
    auto &name = block.second.behaviourName;
    const string param_name = "-block_" + to_string(id);
    const string param_path = "-lib_path_" + to_string(id);
    const string param_from_blocks = "-my_params_" + to_string(id);
    PetscBool set_from_blocks = PETSC_FALSE;
    char char_name[255];
    PetscBool is_param;

    CHKERR PetscOptionsBool(param_from_blocks.c_str(),
                            "set parameters from blocks", "", set_from_blocks,
                            &set_from_blocks, PETSC_NULL);

    CHKERR PetscOptionsString(param_name.c_str(), "name of the behaviour", "",
                              "IsotropicLinearHardeningPlasticity", char_name,
                              255, &is_param);
    if (is_param)
      name = string(char_name);
    string default_lib_path =
        "src/libBehaviour." + string(DEFAULT_LIB_EXTENSION);
    CHKERR PetscOptionsString(
        param_path.c_str(), "path to the behaviour library", "",
        default_lib_path.c_str(), char_name, 255, &is_param);
    if (is_param)
      lib_path = string(char_name);
    auto &mgis_bv_ptr = block.second.mGisBehaviour;
    bool is_finite_strain = false;

    cout << "Loading behaviour: " << name << " from " << lib_path << endl;

    CHKERR check_lib_finite_strain(lib_path, name, is_finite_strain);

    mgis::behaviour::Hypothesis h;
    switch (H) {
    case TRIDIMENSIONAL:
      h = mgis::behaviour::Hypothesis::TRIDIMENSIONAL;
      CHKERR PetscPrintf(PETSC_COMM_WORLD,
                         "MFront material model: TRIDIMENSIONAL \n");
      break;
    case PLANESTRAIN:
      h = mgis::behaviour::Hypothesis::PLANESTRAIN;
      CHKERR PetscPrintf(PETSC_COMM_WORLD,
                         "MFront material model: PLANESTRAIN \n");
      break;
    case AXISYMMETRICAL:
      h = mgis::behaviour::Hypothesis::AXISYMMETRICAL;
      CHKERR PetscPrintf(PETSC_COMM_WORLD,
                         "MFront material model: AXISYMMETRICAL \n");
      break;
    default:
      break;
    }

    if (is_finite_strain) {
      mgis_bv_ptr = boost::make_shared<Behaviour>(load(op, lib_path, name, h));
      block.second.isFiniteStrain = true;
    } else
      mgis_bv_ptr = boost::make_shared<Behaviour>(load(lib_path, name, h));

    CHKERR block.second.setBlockBehaviourData(set_from_blocks);
    for (size_t dd = 0; dd < mgis_bv_ptr->mps.size(); ++dd) {
      double my_param = 0;
      PetscBool is_set = PETSC_FALSE;
      string param_cmd = "-param_" + to_string(id) + "_" + to_string(dd);
      CHKERR PetscOptionsScalar(param_cmd.c_str(), "parameter from cmd", "",
                                my_param, &my_param, &is_set);
      if (!is_set)
        continue;
      setMaterialProperty(block.second.behDataPtr->s0, dd, my_param);
      setMaterialProperty(block.second.behDataPtr->s1, dd, my_param);
    }

    int nb = 0;

    // PRINT PROPERLY WITH SHOWING WHAT WAS ASSIGNED BY THE USER!!!
    CHKERR PetscPrintf(PETSC_COMM_WORLD, "%s behaviour loaded on block %d. \n",
                       mgis_bv_ptr->behaviour.c_str(), block.first);
    if (is_finite_strain)
      CHKERR PetscPrintf(PETSC_COMM_WORLD, "Finite Strain Kinematics \n");
    else
      CHKERR PetscPrintf(PETSC_COMM_WORLD, "Small Strain Kinematics \n");

    CHKERR PetscPrintf(PETSC_COMM_WORLD, "Internal variables: \n");
    for (const auto &is : mgis_bv_ptr->isvs)
      CHKERR PetscPrintf(PETSC_COMM_WORLD, ": %s\n", is.name.c_str());
    CHKERR PetscPrintf(PETSC_COMM_WORLD, "External variables: \n");
    for (const auto &es : mgis_bv_ptr->esvs)
      CHKERR PetscPrintf(PETSC_COMM_WORLD, ": %s\n", es.name.c_str());

    auto it = block.second.behDataPtr->s0.material_properties.begin();
    CHKERR PetscPrintf(PETSC_COMM_WORLD, "Material properties: \n");
    for (const auto &mp : mgis_bv_ptr->mps)
      CHKERR PetscPrintf(PETSC_COMM_WORLD, "%d : %s = %g\n", nb++,
                         mp.name.c_str(), *it++);

    CHKERR PetscPrintf(PETSC_COMM_WORLD, "Real parameters: \n");
    for (auto &p : mgis_bv_ptr->params)
      CHKERR PetscPrintf(PETSC_COMM_WORLD, "%d : %s\n", nb++, p.c_str());
    CHKERR PetscPrintf(PETSC_COMM_WORLD, "Integer parameters: \n");
    for (auto &p : mgis_bv_ptr->iparams)
      CHKERR PetscPrintf(PETSC_COMM_WORLD, "%d : %s\n", nb++, p.c_str());
    CHKERR PetscPrintf(PETSC_COMM_WORLD, "Unsigned short parameters: \n");
    for (auto &p : mgis_bv_ptr->usparams)
      CHKERR PetscPrintf(PETSC_COMM_WORLD, "%d : %s\n", nb++, p.c_str());
  }

  ierr = PetscOptionsEnd();
  CHKERRQ(ierr);

  auto check_behaviours_kinematics = [&](bool &is_finite_kin) {
    MoFEMFunctionBeginHot;
    is_finite_kin =
        commonDataPtr->setOfBlocksData.begin()->second.isFiniteStrain;
    for (auto &block : commonDataPtr->setOfBlocksData) {
      if (block.second.isFiniteStrain != is_finite_kin)
        SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
                "All used MFront behaviours have to be of same kinematics "
                "(small or "
                "large strains)");
    }
    MoFEMFunctionReturnHot(0);
  };

  CHKERR check_behaviours_kinematics(isFiniteKinematics);

  auto get_base = [&]() {
    Range domain_ents;
    MoFEM::Interface &m_field = cOre;
    CHKERR m_field.get_moab().get_entities_by_dimension(0, DIM, domain_ents,
                                                       true);
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

  fieldBase = get_base();
  MOFEM_LOG("WORLD", Sev::inform)
      << "MFront approximation base " << ApproximationBaseNames[fieldBase];

  MoFEMFunctionReturn(0);
};

} // namespace MoFEM

#endif // WITH_MGIS
