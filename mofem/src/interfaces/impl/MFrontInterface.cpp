/** \file MFrontInterface.cpp
 * \brief MFrontInterface
 *
 * MFrontInterface
 *
 */

#ifdef WITH_MGIS

  #include <MFrontInterface.hpp>

  #include "MGIS/LibrariesManager.hxx"
  #include "MGIS/Behaviour/Integrate.hxx"

using namespace FTensor;
using namespace mgis::behaviour;

namespace MoFEM {

  #define VOIGT_VEC_SYMM_3D(VEC)                                               \
    VEC[0], inv_sqr2 *VEC[3], inv_sqr2 *VEC[4], VEC[1], inv_sqr2 *VEC[5], VEC[2]

  #define VOIGT_VEC_SYMM_2D(VEC) VEC[0], inv_sqr2 *VEC[3], VEC[1]

  #define VOIGT_VEC_SYMM_2D_FULL(VEC)                                          \
    VEC[0], inv_sqr2 *VEC[3], 0, VEC[1], 0, VEC[2]

  #define VOIGT_VEC_3D(VEC)                                                    \
    VEC[0], VEC[3], VEC[5], VEC[4], VEC[1], VEC[7], VEC[6], VEC[8], VEC[2]

  #define VOIGT_VEC_2D(VEC) VEC[0], VEC[3], VEC[4], VEC[1]

  #define VOIGT_VEC_2D_FULL(VEC)                                               \
    VEC[0], VEC[3], 0, VEC[4], VEC[1], 0, 0, 0, VEC[2]

MoFEMErrorCode MFrontInterfaceBase::BlockData::setBlockBehaviourData(
    bool set_params_from_blocks) {
  MoFEMFunctionBegin;
  if (mGisBehaviour) {

    auto &mgis_bv = *mGisBehaviour;

    sizeIntVar = getArraySize(mgis_bv.isvs, mgis_bv.hypothesis);
    sizeExtVar = getArraySize(mgis_bv.esvs, mgis_bv.hypothesis);
    sizeGradVar = getArraySize(mgis_bv.gradients, mgis_bv.hypothesis);
    sizeStressVar =
        getArraySize(mgis_bv.thermodynamic_forces, mgis_bv.hypothesis);

    behDataPtr = boost::make_shared<BehaviourData>(BehaviourData{mgis_bv});
    bView = make_view(*behDataPtr);
    const int total_number_of_params = mgis_bv.mps.size();
    // const int total_number_of_params = mgis_bv.mps.size() +
    // mgis_bv.params.size() + mgis_bv.iparams.size() +
    // mgis_bv.usparams.size();

    if (set_params_from_blocks) {

      if (params.size() < total_number_of_params)
        SETERRQ2(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
                 "Not enough parameters supplied for this block. We have %d "
                 "provided where %d are necessary for this block",
                 params.size(), total_number_of_params);

      for (int dd = 0; dd < total_number_of_params; ++dd) {
        setMaterialProperty(behDataPtr->s0, dd, params[dd]);
        setMaterialProperty(behDataPtr->s1, dd, params[dd]);
      }
    }

    if (isFiniteStrain) {
      behDataPtr->K[0] = 0; // no tangent
      behDataPtr->K[1] = 2; // PK1
      behDataPtr->K[2] = 2; // PK1
    } else {
      behDataPtr->K[0] = 0; // no tangent
      behDataPtr->K[1] = 0; // cauchy
    }

    for (auto &mb : {&behDataPtr->s0, &behDataPtr->s1}) {
      mb->dissipated_energy = dIssipation;
      mb->stored_energy = storedEnergy;
      setExternalStateVariable(*mb, 0, externalVariable);
    }
  }

  MoFEMFunctionReturn(0);
}

MoFEMErrorCode MFrontInterfaceBase::CommonData::setBlocks(int dim) {
  MoFEMFunctionBegin;
  string block_name = "MFRONT";
  for (_IT_CUBITMESHSETS_BY_SET_TYPE_FOR_LOOP_(mField, BLOCKSET, it)) {
    if (it->getName().compare(0, block_name.size(), block_name) == 0) {
      std::vector<double> block_data;
      // FIXME: TODO: maybe this should be set only from the command line!!!
      CHKERR it->getAttributes(block_data);
      const int id = it->getMeshsetId();
      EntityHandle meshset = it->getMeshset();
      CHKERR mField.get_moab().get_entities_by_dimension(
          meshset, dim, setOfBlocksData[id].eNts, true);
      for (auto ent : setOfBlocksData[id].eNts)
        blocksIDmap[ent] = id;

      setOfBlocksData[id].iD = id;
      setOfBlocksData[id].params.resize(block_data.size());

      for (int n = 0; n != block_data.size(); n++)
        setOfBlocksData[id].params[n] = block_data[n];
    }
  }

  MoFEMFunctionReturn(0);
}

MoFEMErrorCode MFrontInterfaceBase::CommonData::getInternalVar(
    const EntityHandle fe_ent, const int nb_gauss_pts, const int var_size,
    const int grad_size, const int stress_size, bool is_large_strain) {
  MoFEMFunctionBegin;

  auto mget_tag_data = [&](Tag &m_tag, boost::shared_ptr<MatrixDouble> &m_mat,
                           const int &m_size, bool is_def_grad = false) {
    MoFEMFunctionBeginHot;

    double *tag_data;
    int tag_size;
    rval = mField.get_moab().tag_get_by_ptr(
        m_tag, &fe_ent, 1, (const void **)&tag_data, &tag_size);

    if (rval != MB_SUCCESS || tag_size != m_size * nb_gauss_pts) {
      m_mat->resize(nb_gauss_pts, m_size, false);
      m_mat->clear();
      // initialize deformation gradient properly
      if (is_def_grad && is_large_strain)
        for (int gg = 0; gg != nb_gauss_pts; ++gg) {
          (*m_mat)(gg, 0) = 1;
          (*m_mat)(gg, 1) = 1;
          (*m_mat)(gg, 2) = 1;
        }
      void const *tag_data2[] = {&*m_mat->data().begin()};
      const int tag_size2 = m_mat->data().size();
      CHKERR mField.get_moab().tag_set_by_ptr(m_tag, &fe_ent, 1, tag_data2,
                                              &tag_size2);
    } else {
      MatrixAdaptor tag_vec = MatrixAdaptor(
          nb_gauss_pts, m_size,
          ublas::shallow_array_adaptor<double>(tag_size, tag_data));

      *m_mat = tag_vec;
    }

    MoFEMFunctionReturnHot(0);
  };

  CHKERR mget_tag_data(internalVariableTag, internalVariablePtr, var_size);
  CHKERR mget_tag_data(stressTag, mPrevStressPtr, grad_size);
  CHKERR mget_tag_data(gradientTag, mPrevGradPtr, stress_size, true);

  MoFEMFunctionReturn(0);
}

MoFEMErrorCode
MFrontInterfaceBase::CommonData::setInternalVar(const EntityHandle fe_ent) {
  MoFEMFunctionBegin;

  auto mset_tag_data = [&](Tag &m_tag, boost::shared_ptr<MatrixDouble> &m_mat) {
    MoFEMFunctionBeginHot;
    void const *tag_data[] = {&*m_mat->data().begin()};
    const int tag_size = m_mat->data().size();
    CHKERR mField.get_moab().tag_set_by_ptr(m_tag, &fe_ent, 1, tag_data,
                                            &tag_size);
    MoFEMFunctionReturnHot(0);
  };

  CHKERR mset_tag_data(internalVariableTag, internalVariablePtr);
  CHKERR mset_tag_data(stressTag, mPrevStressPtr);
  CHKERR mset_tag_data(gradientTag, mPrevGradPtr);

  MoFEMFunctionReturn(0);
}

MoFEMErrorCode MFrontInterfaceBase::CommonData::createTags() {
  MoFEMFunctionBegin;
  double def_val = 0.0;
  const int default_length = 0;
  CHKERR mField.get_moab().tag_get_handle(
      "_INTERNAL_VAR", default_length, MB_TYPE_DOUBLE, internalVariableTag,
      MB_TAG_CREAT | MB_TAG_VARLEN | MB_TAG_SPARSE, PETSC_NULL);
  CHKERR mField.get_moab().tag_get_handle(
      "_STRESS_TAG", default_length, MB_TYPE_DOUBLE, stressTag,
      MB_TAG_CREAT | MB_TAG_VARLEN | MB_TAG_SPARSE, PETSC_NULL);
  CHKERR mField.get_moab().tag_get_handle(
      "_GRAD_TAG", default_length, MB_TYPE_DOUBLE, gradientTag,
      MB_TAG_CREAT | MB_TAG_VARLEN | MB_TAG_SPARSE, PETSC_NULL);

  MoFEMFunctionReturn(0);
}

MoFEMErrorCode MFrontInterfaceBase::CommonData::clearTags() {
  MoFEMFunctionBegin;
  double zero = 0;
  for (auto &[id, data] : setOfBlocksData) {
    CHKERR mField.get_moab().tag_clear_data(internalVariableTag, data.eNts,
                                            &zero);
    CHKERR mField.get_moab().tag_clear_data(stressTag, data.eNts, &zero);
    CHKERR mField.get_moab().tag_clear_data(gradientTag, data.eNts, &zero);
  }
  MoFEMFunctionReturn(0);
}

template <ModelHypothesis MH, AssemblyType AT>
MFrontInterface<MH, AT>::MFrontInterface(const MoFEM::Core &core)
    : cOre(const_cast<MoFEM::Core &>(core)) {

  isFiniteKinematics = false;
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

template <ModelHypothesis MH, AssemblyType AT>
MoFEMErrorCode MFrontInterface<MH, AT>::query_interface(
    boost::typeindex::type_index type_index, UnknownInterface **iface) const {
  MoFEMFunctionBeginHot;
  *iface = const_cast<MFrontInterface<MH, AT> *>(this);
  MoFEMFunctionReturnHot(0);
}

template <ModelHypothesis MH, AssemblyType AT>
MoFEMErrorCode MFrontInterface<MH, AT>::getCommandLineParameters() {
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

    auto &lm = mgis::LibrariesManager::get();
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
                              "LinearElasticity", char_name, 255, &is_param);
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
    switch (MH) {
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
}

// template struct OpSaveStress<false, TRIDIMENSIONAL>;
// template struct OpSaveStress<true, TRIDIMENSIONAL>;

// template struct OpSaveStress<false, PLANESTRAIN>;
// template struct OpSaveStress<true, PLANESTRAIN>;

// template struct OpSaveStress<false, AXISYMMETRICAL>;
// template struct OpSaveStress<true, AXISYMMETRICAL>;

// template struct OpSaveGaussPts<TRIDIMENSIONAL>;
// template struct OpSaveGaussPts<PLANESTRAIN>;
// template struct OpSaveGaussPts<AXISYMMETRICAL>;

// template struct MFrontInterfaceBase::OpTangent<Tensor4Pack<3>,
// TRIDIMENSIONAL>; template struct MFrontInterfaceBase::OpTangent<DdgPack<3>,
// TRIDIMENSIONAL>;

// template struct MFrontInterfaceBase::OpTangent<Tensor4Pack<2>, PLANESTRAIN>;
// template struct MFrontInterfaceBase::OpTangent<DdgPack<2>, PLANESTRAIN>;

// template struct MFrontInterfaceBase::OpTangent<Tensor4Pack<2>,
// AXISYMMETRICAL>; template struct MFrontInterfaceBase::OpTangent<DdgPack<2>,
// AXISYMMETRICAL>;

// template <int DIM, ModelHypothesis MH>
// using OpTangentFiniteStrains =
//     struct MFrontInterfaceBase::OpTangent<Tensor4Pack<DIM>, MH>;

// template <int DIM, ModelHypothesis MH>
// using OpTangentSmallStrains =
//     struct MFrontInterfaceBase::OpTangent<DdgPack<DIM>, MH>;

template <int DIM>
using Tensor4Pack =
    FTensor::Tensor4<FTensor::PackPtr<double *, 1>, DIM, DIM, DIM, DIM>;

template <int DIM>
using Tensor2Pack = FTensor::Tensor2<FTensor::PackPtr<double *, 1>, DIM, DIM>;

template <int DIM>
using Tensor1Pack = FTensor::Tensor1<FTensor::PackPtr<double *, 1>, DIM>;

using Tensor1PackCoords = FTensor::Tensor1<FTensor::PackPtr<double *, 3>, 3>;

template <int DIM>
using DdgPack = FTensor::Ddg<FTensor::PackPtr<double *, 1>, DIM, DIM>;

constexpr double sqr2 = boost::math::constants::root_two<double>();
constexpr double inv_sqr2 = boost::math::constants::half_root_two<double>();

template <typename T> inline size_t get_paraview_size(T &vsize) {
  return vsize > 1 ? (vsize > 3 ? 9 : 3) : 1;
};

template <typename T> inline auto get_voigt_vec_symm(T &t_grad) {
  Tensor2_symmetric<double, 3> t_strain;
  Index<'i', 3> i;
  Index<'j', 3> j;
  t_strain(i, j) = (t_grad(i, j) || t_grad(j, i)) / 2;

  array<double, 9> vec_sym{t_strain(0, 0),
                           t_strain(1, 1),
                           t_strain(2, 2),
                           sqr2 * t_strain(0, 1),
                           sqr2 * t_strain(0, 2),
                           sqr2 * t_strain(1, 2),
                           0,
                           0,
                           0};
  return vec_sym;
};

template <typename T> inline auto get_voigt_vec_symm_plane_strain(T &t_grad) {
  Tensor2_symmetric<double, 2> t_strain;
  Index<'i', 2> i;
  Index<'j', 2> j;
  t_strain(i, j) = (t_grad(i, j) || t_grad(j, i)) / 2;

  array<double, 4> vec_sym{t_strain(0, 0), t_strain(1, 1), 0,
                           sqr2 * t_strain(0, 1)};
  return vec_sym;
};

template <typename T1, typename T2, typename T3>
inline auto get_voigt_vec_symm_axisymm(T1 &t_grad, T2 &t_disp, T3 &t_coords) {
  Tensor2_symmetric<double, 2> t_strain;
  Index<'i', 2> i;
  Index<'j', 2> j;
  t_strain(i, j) = (t_grad(i, j) || t_grad(j, i)) / 2;

  array<double, 5> vec_sym{t_strain(0, 0), t_strain(1, 1),
                           t_disp(0) / t_coords(0), sqr2 * t_strain(0, 1)};
  return vec_sym;
};

template <typename T> inline auto get_voigt_vec(T &t_grad) {
  Tensor2<double, 3, 3> F;
  Index<'i', 3> i;
  Index<'j', 3> j;
  F(i, j) = t_grad(i, j) + kronecker_delta(i, j);

  // double det;
  // CHKERR determinantTensor3by3(F, det);
  // if (det < 0)
  //   MOFEM_LOG("WORLD", Sev::error) << "NEGATIVE DET!!!" << det;

  array<double, 9> vec{F(0, 0), F(1, 1), F(2, 2), F(0, 1), F(1, 0),
                       F(0, 2), F(2, 0), F(1, 2), F(2, 1)};

  return vec;
};

template <typename T> inline auto get_voigt_vec_plane_strain(T &t_grad) {
  Tensor2<double, 2, 2> F;
  Index<'i', 2> i;
  Index<'j', 2> j;
  F(i, j) = t_grad(i, j) + kronecker_delta(i, j);

  array<double, 5> vec{F(0, 0), F(1, 1), 1, F(0, 1), F(1, 0)};

  return vec;
};

template <typename T1, typename T2, typename T3>
inline auto get_voigt_vec_axisymm(T1 &t_grad, T2 &t_disp, T3 &t_coords) {
  Tensor2<double, 2, 2> F;
  Index<'i', 2> i;
  Index<'j', 2> j;
  F(i, j) = t_grad(i, j) + kronecker_delta(i, j);

  array<double, 5> vec{F(0, 0), F(1, 1), 1. + t_disp(0) / t_coords(0), F(0, 1),
                       F(1, 0)};

  return vec;
};

template <typename T>
inline auto to_non_symm_3d(const Tensor2_symmetric<T, 3> &symm) {
  Tensor2<double, 3, 3> non_symm;
  Number<0> N0;
  Number<1> N1;
  Number<2> N2;
  non_symm(N0, N0) = symm(N0, N0);
  non_symm(N1, N1) = symm(N1, N1);
  non_symm(N2, N2) = symm(N2, N2);
  non_symm(N0, N1) = non_symm(N1, N0) = symm(N0, N1);
  non_symm(N0, N2) = non_symm(N2, N0) = symm(N0, N2);
  non_symm(N1, N2) = non_symm(N2, N1) = symm(N1, N2);
  return non_symm;
};

template <typename T>
inline auto to_non_symm_2d(const Tensor2_symmetric<T, 2> &symm) {
  Tensor2<double, 2, 2> non_symm;
  Number<0> N0;
  Number<1> N1;
  non_symm(N0, N0) = symm(N0, N0);
  non_symm(N1, N1) = symm(N1, N1);
  non_symm(N0, N1) = non_symm(N1, N0) = symm(N0, N1);
  return non_symm;
};

template <typename T1, typename T2>
inline MoFEMErrorCode get_tensor4_from_voigt(const T1 &K, T2 &D) {
  MoFEMFunctionBeginHot;

  Number<0> N0;
  Number<1> N1;
  Number<2> N2;

  if (std::is_same<T2, Tensor4Pack<3>>::value) { // 3D finite strain
    D(N0, N0, N0, N0) = K[0];
    D(N0, N0, N1, N1) = K[1];
    D(N0, N0, N2, N2) = K[2];
    D(N0, N0, N0, N1) = K[3];
    D(N0, N0, N1, N0) = K[4];
    D(N0, N0, N0, N2) = K[5];
    D(N0, N0, N2, N0) = K[6];
    D(N0, N0, N1, N2) = K[7];
    D(N0, N0, N2, N1) = K[8];
    D(N1, N1, N0, N0) = K[9];
    D(N1, N1, N1, N1) = K[10];
    D(N1, N1, N2, N2) = K[11];
    D(N1, N1, N0, N1) = K[12];
    D(N1, N1, N1, N0) = K[13];
    D(N1, N1, N0, N2) = K[14];
    D(N1, N1, N2, N0) = K[15];
    D(N1, N1, N1, N2) = K[16];
    D(N1, N1, N2, N1) = K[17];
    D(N2, N2, N0, N0) = K[18];
    D(N2, N2, N1, N1) = K[19];
    D(N2, N2, N2, N2) = K[20];
    D(N2, N2, N0, N1) = K[21];
    D(N2, N2, N1, N0) = K[22];
    D(N2, N2, N0, N2) = K[23];
    D(N2, N2, N2, N0) = K[24];
    D(N2, N2, N1, N2) = K[25];
    D(N2, N2, N2, N1) = K[26];
    D(N0, N1, N0, N0) = K[27];
    D(N0, N1, N1, N1) = K[28];
    D(N0, N1, N2, N2) = K[29];
    D(N0, N1, N0, N1) = K[30];
    D(N0, N1, N1, N0) = K[31];
    D(N0, N1, N0, N2) = K[32];
    D(N0, N1, N2, N0) = K[33];
    D(N0, N1, N1, N2) = K[34];
    D(N0, N1, N2, N1) = K[35];
    D(N1, N0, N0, N0) = K[36];
    D(N1, N0, N1, N1) = K[37];
    D(N1, N0, N2, N2) = K[38];
    D(N1, N0, N0, N1) = K[39];
    D(N1, N0, N1, N0) = K[40];
    D(N1, N0, N0, N2) = K[41];
    D(N1, N0, N2, N0) = K[42];
    D(N1, N0, N1, N2) = K[43];
    D(N1, N0, N2, N1) = K[44];
    D(N0, N2, N0, N0) = K[45];
    D(N0, N2, N1, N1) = K[46];
    D(N0, N2, N2, N2) = K[47];
    D(N0, N2, N0, N1) = K[48];
    D(N0, N2, N1, N0) = K[49];
    D(N0, N2, N0, N2) = K[50];
    D(N0, N2, N2, N0) = K[51];
    D(N0, N2, N1, N2) = K[52];
    D(N0, N2, N2, N1) = K[53];
    D(N2, N0, N0, N0) = K[54];
    D(N2, N0, N1, N1) = K[55];
    D(N2, N0, N2, N2) = K[56];
    D(N2, N0, N0, N1) = K[57];
    D(N2, N0, N1, N0) = K[58];
    D(N2, N0, N0, N2) = K[59];
    D(N2, N0, N2, N0) = K[60];
    D(N2, N0, N1, N2) = K[61];
    D(N2, N0, N2, N1) = K[62];
    D(N1, N2, N0, N0) = K[63];
    D(N1, N2, N1, N1) = K[64];
    D(N1, N2, N2, N2) = K[65];
    D(N1, N2, N0, N1) = K[66];
    D(N1, N2, N1, N0) = K[67];
    D(N1, N2, N0, N2) = K[68];
    D(N1, N2, N2, N0) = K[69];
    D(N1, N2, N1, N2) = K[70];
    D(N1, N2, N2, N1) = K[71];
    D(N2, N1, N0, N0) = K[72];
    D(N2, N1, N1, N1) = K[73];
    D(N2, N1, N2, N2) = K[74];
    D(N2, N1, N0, N1) = K[75];
    D(N2, N1, N1, N0) = K[76];
    D(N2, N1, N0, N2) = K[77];
    D(N2, N1, N2, N0) = K[78];
    D(N2, N1, N1, N2) = K[79];
    D(N2, N1, N2, N1) = K[80];
  }

  if (std::is_same<T2,
                   Tensor4Pack<2>>::value) { // plane strain, finite strain
    D(N0, N0, N0, N0) = K[0];
    D(N0, N0, N1, N1) = K[1];
    // D(N0, N0, N2, N2) = K[2];
    D(N0, N0, N0, N1) = K[3];
    D(N0, N0, N1, N0) = K[4];
    D(N1, N1, N0, N0) = K[5];
    D(N1, N1, N1, N1) = K[6];
    // D(N1, N1, N2, N2) = K[7];
    D(N1, N1, N0, N1) = K[8];
    D(N1, N1, N1, N0) = K[9];
    // D(N2, N2, N0, N0) = K[10];
    // D(N2, N2, N1, N1) = K[11];
    // D(N2, N2, N2, N2) = K[12];
    // D(N2, N2, N0, N1) = K[13];
    // D(N2, N2, N1, N0) = K[14];
    D(N0, N1, N0, N0) = K[15];
    D(N0, N1, N1, N1) = K[16];
    // D(N0, N1, N2, N2) = K[17];
    D(N0, N1, N0, N1) = K[18];
    D(N0, N1, N1, N0) = K[19];
    D(N1, N0, N0, N0) = K[20];
    D(N1, N0, N1, N1) = K[21];
    // D(N1, N0, N2, N2) = K[22];
    D(N1, N0, N0, N1) = K[23];
    D(N1, N0, N1, N0) = K[24];
  }

  if (std::is_same<T2, DdgPack<3>>::value) { // 3D small strain
    D(N0, N0, N0, N0) = K[0];
    D(N0, N0, N1, N1) = K[1];
    D(N0, N0, N2, N2) = K[2];

    D(N0, N0, N0, N1) = inv_sqr2 * K[3];
    D(N0, N0, N0, N2) = inv_sqr2 * K[4];
    D(N0, N0, N1, N2) = inv_sqr2 * K[5];

    D(N1, N1, N0, N0) = K[6];
    D(N1, N1, N1, N1) = K[7];
    D(N1, N1, N2, N2) = K[8];

    D(N1, N1, N0, N1) = inv_sqr2 * K[9];
    D(N1, N1, N0, N2) = inv_sqr2 * K[10];
    D(N1, N1, N1, N2) = inv_sqr2 * K[11];

    D(N2, N2, N0, N0) = K[12];
    D(N2, N2, N1, N1) = K[13];
    D(N2, N2, N2, N2) = K[14];

    D(N2, N2, N0, N1) = inv_sqr2 * K[15];
    D(N2, N2, N0, N2) = inv_sqr2 * K[16];
    D(N2, N2, N1, N2) = inv_sqr2 * K[17];

    D(N0, N1, N0, N0) = inv_sqr2 * K[18];
    D(N0, N1, N1, N1) = inv_sqr2 * K[19];
    D(N0, N1, N2, N2) = inv_sqr2 * K[20];

    D(N0, N1, N0, N1) = 0.5 * K[21];
    D(N0, N1, N0, N2) = 0.5 * K[22];
    D(N0, N1, N1, N2) = 0.5 * K[23];

    D(N0, N2, N0, N0) = inv_sqr2 * K[24];
    D(N0, N2, N1, N1) = inv_sqr2 * K[25];
    D(N0, N2, N2, N2) = inv_sqr2 * K[26];

    D(N0, N2, N0, N1) = 0.5 * K[27];
    D(N0, N2, N0, N2) = 0.5 * K[28];
    D(N0, N2, N1, N2) = 0.5 * K[29];

    D(N1, N2, N0, N0) = inv_sqr2 * K[30];
    D(N1, N2, N1, N1) = inv_sqr2 * K[31];
    D(N1, N2, N2, N2) = inv_sqr2 * K[32];

    D(N1, N2, N0, N1) = 0.5 * K[33];
    D(N1, N2, N0, N2) = 0.5 * K[34];
    D(N1, N2, N1, N2) = 0.5 * K[35];
  }

  if (std::is_same<T2, DdgPack<2>>::value) { // plane strain, small strain

    D(N0, N0, N0, N0) = K[0];
    D(N0, N0, N1, N1) = K[1];
    // D(N0, N0, N2, N2) = K[2];

    D(N0, N0, N0, N1) = inv_sqr2 * K[3];

    D(N1, N1, N0, N0) = K[4];
    D(N1, N1, N1, N1) = K[5];
    // D(N1, N1, N2, N2) = K[6];

    D(N1, N1, N0, N1) = inv_sqr2 * K[7];

    // D(N2, N2, N0, N0) = K[8];
    // D(N2, N2, N1, N1) = K[9];
    // D(N2, N2, N2, N2) = K[10];

    // D(N2, N2, N0, N1) = inv_sqr2 * K[11];

    D(N0, N1, N0, N0) = inv_sqr2 * K[12];
    D(N0, N1, N1, N1) = inv_sqr2 * K[13];

    // D(N0, N1, N2, N2) = inv_sqr2 * K[14];

    D(N0, N1, N0, N1) = 0.5 * K[15];
  }

  MoFEMFunctionReturnHot(0);
};

template <bool IS_LARGE_STRAIN, typename T1, typename T2>
inline MoFEMErrorCode get_full_tensor4_from_voigt(const T1 &K, T2 &D) {
  MoFEMFunctionBeginHot;

  Number<0> N0;
  Number<1> N1;
  Number<2> N2;

  if constexpr (IS_LARGE_STRAIN) { // 2D finite strain, full tensor
    D(N0, N0, N0, N0) = K[0];
    D(N0, N0, N1, N1) = K[1];
    D(N0, N0, N2, N2) = K[2];
    D(N0, N0, N0, N1) = K[3];
    D(N0, N0, N1, N0) = K[4];
    D(N1, N1, N0, N0) = K[5];
    D(N1, N1, N1, N1) = K[6];
    D(N1, N1, N2, N2) = K[7];
    D(N1, N1, N0, N1) = K[8];
    D(N1, N1, N1, N0) = K[9];
    D(N2, N2, N0, N0) = K[10];
    D(N2, N2, N1, N1) = K[11];
    D(N2, N2, N2, N2) = K[12];
    D(N2, N2, N0, N1) = K[13];
    D(N2, N2, N1, N0) = K[14];
    D(N0, N1, N0, N0) = K[15];
    D(N0, N1, N1, N1) = K[16];
    D(N0, N1, N2, N2) = K[17];
    D(N0, N1, N0, N1) = K[18];
    D(N0, N1, N1, N0) = K[19];
    D(N1, N0, N0, N0) = K[20];
    D(N1, N0, N1, N1) = K[21];
    D(N1, N0, N2, N2) = K[22];
    D(N1, N0, N0, N1) = K[23];
    D(N1, N0, N1, N0) = K[24];
  } else { // 2D small strain, full tensor

    D(N0, N0, N0, N0) = K[0];
    D(N0, N0, N1, N1) = K[1];
    D(N0, N0, N2, N2) = K[2];

    D(N0, N0, N0, N1) = inv_sqr2 * K[3];

    D(N1, N1, N0, N0) = K[4];
    D(N1, N1, N1, N1) = K[5];
    D(N1, N1, N2, N2) = K[6];

    D(N1, N1, N0, N1) = inv_sqr2 * K[7];

    D(N2, N2, N0, N0) = K[8];
    D(N2, N2, N1, N1) = K[9];
    D(N2, N2, N2, N2) = K[10];

    D(N2, N2, N0, N1) = inv_sqr2 * K[11];
    D(N2, N2, N1, N0) = D(N2, N2, N0, N1);

    D(N0, N1, N0, N0) = inv_sqr2 * K[12];
    D(N0, N1, N1, N1) = inv_sqr2 * K[13];

    D(N0, N1, N2, N2) = inv_sqr2 * K[14];
    D(N1, N0, N2, N2) = D(N0, N1, N2, N2);

    D(N0, N1, N0, N1) = 0.5 * K[15];
  }

  MoFEMFunctionReturnHot(0);
};

template <typename T> T get_tangent_tensor(MatrixDouble &mat);

template <>
Tensor4Pack<3> get_tangent_tensor<Tensor4Pack<3>>(MatrixDouble &mat) {
  return getFTensor4FromMat<3, 3, 3, 3>(mat);
}

template <>
Tensor4Pack<2> get_tangent_tensor<Tensor4Pack<2>>(MatrixDouble &mat) {
  return getFTensor4FromMat<2, 2, 2, 2>(mat);
}

template <> DdgPack<3> get_tangent_tensor<DdgPack<3>>(MatrixDouble &mat) {
  return getFTensor4DdgFromMat<3, 3>(mat);
}

template <> DdgPack<2> get_tangent_tensor<DdgPack<2>>(MatrixDouble &mat) {
  return getFTensor4DdgFromMat<2, 2>(mat);
}

template <bool IS_LARGE_STRAIN, ModelHypothesis MH>
MoFEMErrorCode
mgis_integration(size_t gg, Tensor2Pack<MFrontEleType<MH>::SPACE_DIM> &t_grad,
                 Tensor1Pack<MFrontEleType<MH>::SPACE_DIM> &t_disp,
                 Tensor1PackCoords &t_coords,
                 MFrontInterfaceBase::CommonData &common_data,
                 MFrontInterfaceBase::BlockData &block_data) {
  MoFEMFunctionBegin;

  static constexpr int DIM = MFrontEleType<MH>::SPACE_DIM;

  int check_integration;
  MatrixDouble &mat_int = *common_data.internalVariablePtr;
  MatrixDouble &mat_grad0 = *common_data.mPrevGradPtr;
  MatrixDouble &mat_stress0 = *common_data.mPrevStressPtr;

  int &size_of_vars = block_data.sizeIntVar;
  int &size_of_grad = block_data.sizeGradVar;
  int &size_of_stress = block_data.sizeStressVar;

  auto &mgis_bv = *block_data.mGisBehaviour;

  if constexpr (IS_LARGE_STRAIN) {
    if constexpr (DIM == 3)
      setGradient(block_data.behDataPtr->s1, 0, size_of_grad,
                  &*get_voigt_vec(t_grad).data());
    if constexpr (DIM == 2) {
      if constexpr (MH == AXISYMMETRICAL) {
        setGradient(block_data.behDataPtr->s1, 0, size_of_grad,
                    &*get_voigt_vec_axisymm(t_grad, t_disp, t_coords).data());
      } else
        setGradient(block_data.behDataPtr->s1, 0, size_of_grad,
                    &*get_voigt_vec_plane_strain(t_grad).data());
    }
  } else {
    if constexpr (DIM == 3)
      setGradient(block_data.behDataPtr->s1, 0, size_of_grad,
                  &*get_voigt_vec_symm(t_grad).data());
    if constexpr (DIM == 2) {
      if constexpr (MH == AXISYMMETRICAL)
        setGradient(
            block_data.behDataPtr->s1, 0, size_of_grad,
            &*get_voigt_vec_symm_axisymm(t_grad, t_disp, t_coords).data());
      else
        setGradient(block_data.behDataPtr->s1, 0, size_of_grad,
                    &*get_voigt_vec_symm_plane_strain(t_grad).data());
    }
  }

  auto grad0_vec =
      getVectorAdaptor(&mat_grad0.data()[gg * size_of_grad], size_of_grad);
  setGradient(block_data.behDataPtr->s0, 0, size_of_grad, &*grad0_vec.begin());

  auto stress0_vec = getVectorAdaptor(&mat_stress0.data()[gg * size_of_stress],
                                      size_of_stress);
  setThermodynamicForce(block_data.behDataPtr->s0, 0, size_of_stress,
                        &*stress0_vec.begin());

  if (size_of_vars) {
    auto internal_var =
        getVectorAdaptor(&mat_int.data()[gg * size_of_vars], size_of_vars);
    setInternalStateVariable(block_data.behDataPtr->s0, 0, size_of_vars,
                             &*internal_var.begin());
  }

  check_integration = mgis::behaviour::integrate(block_data.bView, mgis_bv);
  switch (check_integration) {
  case -1:
    MOFEM_LOG("WORLD", Sev::error) << "Mfront integration failed";
    break;
  case 0:
    MOFEM_LOG("WORLD", Sev::warning)
        << "Mfront integration succeeded but results are unreliable";
    break;
  case 1:
  default:
    break;
  }

  MoFEMFunctionReturn(0);
}

template <bool UPDATE, bool IS_LARGE_STRAIN, ModelHypothesis MH>
struct OpStressTmp : public MFrontEleType<MH>::DomainEleOp {
  static constexpr int DIM = MFrontEleType<MH>::SPACE_DIM;
  using DomainEleOp = typename MFrontEleType<MH>::DomainEleOp;

  OpStressTmp(
      const std::string field_name,
      boost::shared_ptr<MFrontInterfaceBase::CommonData> common_data_ptr,
      boost::shared_ptr<FEMethod> monitor_ptr)
      : DomainEleOp(field_name, DomainEleOp::OPROW),
        commonDataPtr(common_data_ptr), monitorPtr(monitor_ptr) {
    std::fill(&DomainEleOp::doEntities[MBEDGE],
              &DomainEleOp::doEntities[MBMAXTYPE], false);
  }
  MoFEMErrorCode doWork(int side, EntityType type, EntData &data);

private:
  boost::shared_ptr<MFrontInterfaceBase::CommonData> commonDataPtr;
  boost::shared_ptr<FEMethod> monitorPtr;
};

template <ModelHypothesis MH>
using OpUpdateVariablesFiniteStrains = OpStressTmp<true, true, MH>;

template <ModelHypothesis MH>
using OpUpdateVariablesSmallStrains = OpStressTmp<true, false, MH>;

template <ModelHypothesis MH>
using OpStressFiniteStrains = OpStressTmp<false, true, MH>;

template <ModelHypothesis MH>
using OpStressSmallStrains = OpStressTmp<false, false, MH>;

template <bool UPDATE, bool IS_LARGE_STRAIN, ModelHypothesis MH>
MoFEMErrorCode OpStressTmp<UPDATE, IS_LARGE_STRAIN, MH>::doWork(int side,
                                                                EntityType type,
                                                                EntData &data) {
  MoFEMFunctionBegin;

  Index<'i', 3> i;
  Index<'j', 3> j;

  Index<'I', 2> I;
  Index<'J', 2> J;

  const size_t nb_gauss_pts = commonDataPtr->mGradPtr->size2();
  auto fe_ent = DomainEleOp::getNumeredEntFiniteElementPtr()->getEnt();
  auto id = commonDataPtr->blocksIDmap.at(fe_ent);
  auto &dAta = commonDataPtr->setOfBlocksData.at(id);
  auto &mgis_bv = *dAta.mGisBehaviour;

  if (monitorPtr == nullptr)
    SETERRQ(PETSC_COMM_WORLD, MOFEM_NOT_FOUND,
            "Time Monitor (FEMethod) has not been set for MFrontInterface. "
            "Make sure to call setMonitorPtr(monitor_ptr) before calling "
            "setOperators().");

  dAta.setTag(MFrontInterfaceBase::DataTags::RHS);
  dAta.behDataPtr->dt = monitorPtr->ts_dt;
  dAta.bView.dt = monitorPtr->ts_dt;

  CHKERR commonDataPtr->getInternalVar(fe_ent, nb_gauss_pts, dAta.sizeIntVar,
                                       dAta.sizeGradVar, dAta.sizeStressVar,
                                       IS_LARGE_STRAIN);

  MatrixDouble &mat_int = *commonDataPtr->internalVariablePtr;
  MatrixDouble &mat_grad0 = *commonDataPtr->mPrevGradPtr;
  MatrixDouble &mat_stress0 = *commonDataPtr->mPrevStressPtr;

  auto t_grad = getFTensor2FromMat<DIM, DIM>(*(commonDataPtr->mGradPtr));
  auto t_disp = getFTensor1FromMat<DIM>(*(commonDataPtr->mDispPtr));
  auto t_coords = DomainEleOp::getFTensor1CoordsAtGaussPts();

  commonDataPtr->mStressPtr->resize(DIM * DIM, nb_gauss_pts);
  commonDataPtr->mStressPtr->clear();
  auto t_stress = getFTensor2FromMat<DIM, DIM>(*(commonDataPtr->mStressPtr));

  commonDataPtr->mFullStressPtr->resize(3 * 3, nb_gauss_pts);
  commonDataPtr->mFullStressPtr->clear();
  auto t_full_stress =
      getFTensor2FromMat<3, 3>(*(commonDataPtr->mFullStressPtr));

  for (size_t gg = 0; gg != nb_gauss_pts; ++gg) {

    CHKERR mgis_integration<IS_LARGE_STRAIN, MH>(gg, t_grad, t_disp, t_coords,
                                                 *commonDataPtr, dAta);

    if constexpr (DIM == 3) {
      Tensor2<double, 3, 3> t_force;
      if constexpr (IS_LARGE_STRAIN) {
        t_force = Tensor2<double, 3, 3>(
            VOIGT_VEC_3D(getThermodynamicForce(dAta.behDataPtr->s1, 0)));
      } else {
        t_force = to_non_symm_3d(Tensor2_symmetric<double, 3>(
            VOIGT_VEC_SYMM_3D(getThermodynamicForce(dAta.behDataPtr->s1, 0))));
      }
      t_stress(i, j) = t_force(i, j);
    } else if constexpr (DIM == 2) {
      Tensor2<double, 2, 2> t_force;
      Tensor2<double, 3, 3> t_full_force;
      if constexpr (IS_LARGE_STRAIN) {
        t_force = Tensor2<double, 2, 2>(
            VOIGT_VEC_2D(getThermodynamicForce(dAta.behDataPtr->s1, 0)));
        t_full_force = Tensor2<double, 3, 3>(
            VOIGT_VEC_2D_FULL(getThermodynamicForce(dAta.behDataPtr->s1, 0)));
      } else {
        t_force = to_non_symm_2d(Tensor2_symmetric<double, 2>(
            VOIGT_VEC_SYMM_2D(getThermodynamicForce(dAta.behDataPtr->s1, 0))));
        t_full_force =
            to_non_symm_3d(Tensor2_symmetric<double, 3>(VOIGT_VEC_SYMM_2D_FULL(
                getThermodynamicForce(dAta.behDataPtr->s1, 0))));
      }
      t_stress(I, J) = t_force(I, J);
      t_full_stress(i, j) = t_full_force(i, j);
    }

    if constexpr (UPDATE) {
      for (int dd = 0; dd != dAta.sizeIntVar; ++dd) {
        mat_int(gg, dd) = *getInternalStateVariable(dAta.behDataPtr->s1, dd);
      }
      for (int dd = 0; dd != dAta.sizeGradVar; ++dd) {
        mat_grad0(gg, dd) = *getGradient(dAta.behDataPtr->s1, dd);
      }
      for (int dd = 0; dd != dAta.sizeStressVar; ++dd) {
        mat_stress0(gg, dd) = *getThermodynamicForce(dAta.behDataPtr->s1, dd);
      }
    }

    ++t_stress;
    ++t_full_stress;
    ++t_grad;
    ++t_disp;
    ++t_coords;
  }

  if constexpr (UPDATE) {
    CHKERR commonDataPtr->setInternalVar(fe_ent);
    // mfront_dt_prop = mfront_dt * b_view.rdt;
  }

  MoFEMFunctionReturn(0);
}

template <typename T, ModelHypothesis MH>
struct OpTangent : public MFrontEleType<MH>::DomainEleOp {
  static constexpr int DIM = MFrontEleType<MH>::SPACE_DIM;
  using DomainEleOp = typename MFrontEleType<MH>::DomainEleOp;

  OpTangent(const std::string field_name,
            boost::shared_ptr<MFrontInterfaceBase::CommonData> common_data_ptr,
            boost::shared_ptr<FEMethod> monitor_ptr)
      : DomainEleOp(field_name, DomainEleOp::OPROW),
        commonDataPtr(common_data_ptr), monitorPtr(monitor_ptr) {
    std::fill(&DomainEleOp::doEntities[MBEDGE],
              &DomainEleOp::doEntities[MBMAXTYPE], false);
  }
  MoFEMErrorCode doWork(int side, EntityType type, EntData &data);

private:
  boost::shared_ptr<MFrontInterfaceBase::CommonData> commonDataPtr;
  boost::shared_ptr<FEMethod> monitorPtr;
};

template <ModelHypothesis MH>
using OpTangentFiniteStrains =
    struct OpTangent<Tensor4Pack<MFrontEleType<MH>::SPACE_DIM>, MH>;

template <ModelHypothesis MH>
using OpTangentSmallStrains =
    struct OpTangent<DdgPack<MFrontEleType<MH>::SPACE_DIM>, MH>;

template <typename T, ModelHypothesis MH>
MoFEMErrorCode OpTangent<T, MH>::doWork(int side, EntityType type,
                                        EntData &data) {
  MoFEMFunctionBegin;

  const size_t nb_gauss_pts = commonDataPtr->mGradPtr->size2();
  auto fe_ent = DomainEleOp::getNumeredEntFiniteElementPtr()->getEnt();
  auto id = commonDataPtr->blocksIDmap.at(fe_ent);
  auto &dAta = commonDataPtr->setOfBlocksData.at(id);
  auto &mgis_bv = *dAta.mGisBehaviour;

  if (monitorPtr == nullptr)
    SETERRQ(PETSC_COMM_WORLD, MOFEM_NOT_FOUND,
            "Time Monitor (FEMethod) has not been set for MFrontInterface. "
            "Make sure to call setMonitorPtr(monitor_ptr) before calling "
            "setOperators().");

  dAta.setTag(MFrontInterfaceBase::DataTags::LHS);
  dAta.behDataPtr->dt = monitorPtr->ts_dt;
  dAta.bView.dt = monitorPtr->ts_dt;

  constexpr bool IS_LARGE_STRAIN = std::is_same<T, Tensor4Pack<3>>::value ||
                                   std::is_same<T, Tensor4Pack<2>>::value;

  CHKERR commonDataPtr->getInternalVar(fe_ent, nb_gauss_pts, dAta.sizeIntVar,
                                       dAta.sizeGradVar, dAta.sizeStressVar,
                                       IS_LARGE_STRAIN);

  MatrixDouble &S_E = *(commonDataPtr->materialTangentPtr);
  MatrixDouble &F_E = *(commonDataPtr->mFullTangentPtr);

  size_t tens_size = 36;

  if constexpr (DIM == 2) {
    // plane strain
    if constexpr (IS_LARGE_STRAIN)
      tens_size = 16;
    else
      tens_size = 9;
  } else {
    if constexpr (IS_LARGE_STRAIN) // for finite strains
      tens_size = 81;
  }
  // FIXME plain strain
  S_E.resize(tens_size, nb_gauss_pts, false);
  auto D1 = get_tangent_tensor<T>(S_E);

  size_t full_tens_size = 81;
  F_E.resize(full_tens_size, nb_gauss_pts, false);
  F_E.clear();

  auto D2 = get_tangent_tensor<Tensor4Pack<3>>(F_E);

  auto t_grad = getFTensor2FromMat<DIM, DIM>(*(commonDataPtr->mGradPtr));
  auto t_disp = getFTensor1FromMat<DIM>(*(commonDataPtr->mDispPtr));
  auto t_coords = DomainEleOp::getFTensor1CoordsAtGaussPts();

  for (size_t gg = 0; gg != nb_gauss_pts; ++gg) {

    CHKERR mgis_integration<IS_LARGE_STRAIN, MH>(gg, t_grad, t_disp, t_coords,
                                                 *commonDataPtr, dAta);

    CHKERR get_tensor4_from_voigt(&*dAta.behDataPtr->K.begin(), D1);
    CHKERR get_full_tensor4_from_voigt<IS_LARGE_STRAIN>(
        &*dAta.behDataPtr->K.begin(), D2);

    ++D1;
    ++D2;
    ++t_grad;
    ++t_disp;
    ++t_coords;
  }

  MoFEMFunctionReturn(0);
}

template <AssemblyType AT>
struct OpAxisymmetricRhs
    : public FormsIntegrators<
          MFrontEleType<AXISYMMETRICAL>::DomainEleOp>::Assembly<AT>::OpBase {
  using DomainEleOp = typename MFrontEleType<AXISYMMETRICAL>::DomainEleOp;
  using OpBase = typename FormsIntegrators<DomainEleOp>::Assembly<AT>::OpBase;

  OpAxisymmetricRhs(
      const std::string field_name,
      boost::shared_ptr<MFrontInterfaceBase::CommonData> common_data_ptr)
      : OpBase(field_name, field_name, DomainEleOp::OPROW),
        commonDataPtr(common_data_ptr) {};

private:
  boost::shared_ptr<MFrontInterfaceBase::CommonData> commonDataPtr;

  MoFEMErrorCode iNtegrate(EntData &row_data);
};

template <AssemblyType AT>
MoFEMErrorCode OpAxisymmetricRhs<AT>::iNtegrate(EntData &row_data) {
  MoFEMFunctionBegin;

  auto fe_ent = OpBase::getNumeredEntFiniteElementPtr()->getEnt();
  auto id = commonDataPtr->blocksIDmap.at(fe_ent);
  auto &dAta = commonDataPtr->setOfBlocksData.at(id);

  // get element volume
  const double vol = OpBase::getMeasure();
  // get integration weights
  auto t_w = OpBase::getFTensor0IntegrationWeight();
  // get coordinate at integration points

  auto t_full_stress =
      getFTensor2FromMat<3, 3>(*(commonDataPtr->mFullStressPtr));

  // loop over integration points
  for (int gg = 0; gg != OpBase::nbIntegrationPts; gg++) {

    auto t_nf = OpBase::template getNf<2>();

    FTensor::Tensor0<double *> t_base(&row_data.getN()(gg, 0));

    // take into account Jacobian
    const double alpha = t_w * vol * 2. * M_PI;
    // loop over rows base functions
    for (int rr = 0; rr != OpBase::nbRows / 2; ++rr) {

      t_nf(0) += alpha * t_full_stress(2, 2) * t_base;

      ++t_base;
      ++t_nf;
    }

    ++t_full_stress;
    ++t_w; // move to another integration weight
  }

  MoFEMFunctionReturn(0);
}

template <AssemblyType AT>
struct OpAxisymmetricLhs
    : public FormsIntegrators<
          MFrontEleType<AXISYMMETRICAL>::DomainEleOp>::Assembly<AT>::OpBase {
  using DomainEleOp = typename MFrontEleType<AXISYMMETRICAL>::DomainEleOp;
  using OpBase = typename FormsIntegrators<DomainEleOp>::Assembly<AT>::OpBase;

  OpAxisymmetricLhs(
      const std::string field_name,
      boost::shared_ptr<MFrontInterfaceBase::CommonData> common_data_ptr)
      : OpBase(field_name, field_name, DomainEleOp::OPROWCOL),
        commonDataPtr(common_data_ptr) {};

private:
  boost::shared_ptr<MFrontInterfaceBase::CommonData> commonDataPtr;

  MoFEMErrorCode iNtegrate(EntData &row_data, EntData &col_data);
};

template <AssemblyType AT>
MoFEMErrorCode OpAxisymmetricLhs<AT>::iNtegrate(EntData &row_data,
                                                EntData &col_data) {
  MoFEMFunctionBegin;

  auto fe_ent = OpBase::getNumeredEntFiniteElementPtr()->getEnt();
  auto id = commonDataPtr->blocksIDmap.at(fe_ent);
  auto &dAta = commonDataPtr->setOfBlocksData.at(id);

  // get element volume
  const double vol = OpBase::getMeasure();
  // get integration weights
  auto t_w = OpBase::getFTensor0IntegrationWeight();
  // get base function gradient on rows
  auto t_row_base = row_data.getFTensor0N();
  // get derivatives of base functions on rows
  auto t_row_diff_base = row_data.getFTensor1DiffN<2>();
  // get coordinate at integration points
  auto t_coords = OpBase::getFTensor1CoordsAtGaussPts();

  Number<0> N0;
  Number<1> N1;
  Number<2> N2;

  auto t_D =
      getFTensor4FromMat<3, 3, 3, 3, 1>(*(commonDataPtr->mFullTangentPtr));

  // loop over integration points
  for (int gg = 0; gg != OpBase::nbIntegrationPts; gg++) {

    // Cylinder radius
    const double r_cylinder = t_coords(0);

    // take into account Jacobean
    const double alpha = t_w * vol * 2. * M_PI;

    // loop over rows base functions
    int rr = 0;
    for (; rr != OpBase::nbRows / 2; ++rr) {

      // get sub matrix for the row
      auto t_m = OpBase::template getLocMat<2>(2 * rr);

      // get derivatives of base functions for columns
      auto t_col_diff_base = col_data.getFTensor1DiffN<2>(gg, 0);
      // get base functions for columns
      auto t_col_base = col_data.getFTensor0N(gg, 0);

      // loop over columns
      for (int cc = 0; cc != OpBase::nbCols / 2; cc++) {

        t_m(0, 0) +=
            alpha * t_D(N0, N0, N2, N2) * t_col_base * t_row_diff_base(0);

        t_m(1, 0) +=
            alpha * t_D(N1, N1, N2, N2) * t_col_base * t_row_diff_base(1);

        t_m(0, 0) +=
            alpha * t_D(N2, N2, N0, N0) * t_col_diff_base(0) * t_row_base;

        t_m(0, 1) +=
            alpha * t_D(N2, N2, N1, N1) * t_col_diff_base(1) * t_row_base;

        t_m(0, 0) += alpha * (t_D(N2, N2, N2, N2) / r_cylinder) * t_col_base *
                     t_row_base;

        t_m(0, 0) +=
            alpha * t_D(N2, N2, N0, N1) * t_col_diff_base(1) * t_row_base;

        t_m(0, 1) +=
            alpha * t_D(N2, N2, N1, N0) * t_col_diff_base(0) * t_row_base;

        t_m(0, 0) +=
            alpha * t_D(N0, N1, N2, N2) * t_col_base * t_row_diff_base(1);

        t_m(1, 0) +=
            alpha * t_D(N1, N0, N2, N2) * t_col_base * t_row_diff_base(0);

        ++t_col_base;
        ++t_col_diff_base;
        ++t_m;
      }

      ++t_row_base;
      ++t_row_diff_base;
    }

    for (; rr < OpBase::nbRowBaseFunctions; ++rr) {
      ++t_row_base;
      ++t_row_diff_base;
    }

    ++t_coords;
    ++t_w; // move to another integration weight
    ++t_D;
  }

  MoFEMFunctionReturn(0);
}

template <ModelHypothesis MH, AssemblyType AT>
MoFEMErrorCode MFrontInterface<MH, AT>::opFactoryDomainRhs(
    boost::ptr_deque<ForcesAndSourcesCore::UserDataOperator> &pip,
    std::string field_name) {
  MoFEMFunctionBegin;

  auto jacobian = [&](const double r, const double, const double) {
    if (MH == AXISYMMETRICAL)
      return 2. * M_PI * r;
    else
      return 1.;
  };

  auto add_domain_ops_rhs = [&](auto &pipeline) {
    if (isFiniteKinematics)
      pipeline.push_back(
          new OpStressFiniteStrains<MH>(field_name, commonDataPtr, monitorPtr));
    else
      pipeline.push_back(
          new OpStressSmallStrains<MH>(field_name, commonDataPtr, monitorPtr));

    pipeline.push_back(
        new OpInternalForce(field_name, commonDataPtr->mStressPtr, jacobian));

    if (MH == AXISYMMETRICAL)
      pipeline.push_back(new OpAxisymmetricRhs<AT>(field_name, commonDataPtr));
  };

  auto add_domain_base_ops = [&](auto &pipeline) {
    pipeline.push_back(new OpCalculateVectorFieldValues<DIM>(
        field_name, commonDataPtr->mDispPtr));
    pipeline.push_back(new OpCalculateVectorFieldGradient<DIM, DIM>(
        field_name, commonDataPtr->mGradPtr));
  };

  add_domain_base_ops(pip);
  add_domain_ops_rhs(pip);

  MoFEMFunctionReturn(0);
}

template <ModelHypothesis MH, AssemblyType AT>
MoFEMErrorCode MFrontInterface<MH, AT>::opFactoryDomainLhs(
    boost::ptr_deque<ForcesAndSourcesCore::UserDataOperator> &pip,
    std::string field_name) {
  MoFEMFunctionBegin;

  auto jacobian = [&](const double r, const double, const double) {
    if (MH == AXISYMMETRICAL)
      return 2. * M_PI * r;
    else
      return 1.;
  };

  auto add_domain_ops_lhs = [&](auto &pipeline) {
    if (isFiniteKinematics) {
      pipeline.push_back(new OpTangentFiniteStrains<MH>(
          field_name, commonDataPtr, monitorPtr));
      pipeline.push_back(new OpAssembleLhsFiniteStrains(
          field_name, field_name, commonDataPtr->materialTangentPtr, jacobian));
    } else {
      pipeline.push_back(
          new OpTangentSmallStrains<MH>(field_name, commonDataPtr, monitorPtr));
      pipeline.push_back(new OpAssembleLhsSmallStrains(
          field_name, field_name, commonDataPtr->materialTangentPtr, nullptr,
          jacobian));
    }
    if (MH == AXISYMMETRICAL)
      pipeline.push_back(new OpAxisymmetricLhs<AT>(field_name, commonDataPtr));
  };

  auto add_domain_base_ops = [&](auto &pipeline) {
    pipeline.push_back(new OpCalculateVectorFieldValues<DIM>(
        field_name, commonDataPtr->mDispPtr));
    pipeline.push_back(new OpCalculateVectorFieldGradient<DIM, DIM>(
        field_name, commonDataPtr->mGradPtr));
  };

  add_domain_base_ops(pip);
  add_domain_ops_lhs(pip);

  MoFEMFunctionReturn(0);
}

template <ModelHypothesis MH, AssemblyType AT>
MoFEMErrorCode MFrontInterface<MH, AT>::setUpdateElementVariablesOperators(
    ForcesAndSourcesCore::RuleHookFun rule,
    std::string field_name) {
  MoFEMFunctionBegin;

  auto &moab_gauss = *moabGaussIntPtr;

  MoFEM::Interface &m_field = cOre;
  updateIntVariablesElePtr = boost::make_shared<DomainEle>(m_field);

  updateIntVariablesElePtr->getRuleHook = rule;

  CHKERR AddHOOps<DIM, DIM, DIM>::add(updateIntVariablesElePtr->getOpPtrVector(), {H1});

  updateIntVariablesElePtr->getOpPtrVector().push_back(
      new OpCalculateVectorFieldGradient<DIM, DIM>(field_name,
                                                   commonDataPtr->mGradPtr));
  updateIntVariablesElePtr->getOpPtrVector().push_back(
      new OpCalculateVectorFieldValues<DIM>(field_name,
                                            commonDataPtr->mDispPtr));
  if (isFiniteKinematics)
    updateIntVariablesElePtr->getOpPtrVector().push_back(
        new OpUpdateVariablesFiniteStrains<MH>(field_name, commonDataPtr,
                                              monitorPtr));
  else
    updateIntVariablesElePtr->getOpPtrVector().push_back(
        new OpUpdateVariablesSmallStrains<MH>(field_name, commonDataPtr,
                                             monitorPtr));
  if (saveGauss && (MH == TRIDIMENSIONAL)) {
    // FIXME: decide if needed for 2D
    // updateIntVariablesElePtr->getOpPtrVector().push_back(
    //     new OpSaveGaussPts<MH>(field_name, moab_gauss, commonDataPtr));
  }

  MoFEMFunctionReturn(0);
}

template <ModelHypothesis MH, AssemblyType AT>
MoFEMErrorCode
MFrontInterface<MH, AT>::updateElementVariables(SmartPetscObj<DM> dm,
                                                std::string fe_name) {
  MoFEMFunctionBegin;
  CHKERR DMoFEMLoopFiniteElements(dm, fe_name, updateIntVariablesElePtr);
  MoFEMFunctionReturn(0);
}

template struct MFrontInterface<TRIDIMENSIONAL, AssemblyType::PETSC>;
template struct MFrontInterface<AXISYMMETRICAL, AssemblyType::PETSC>;
template struct MFrontInterface<PLANESTRAIN, AssemblyType::PETSC>;

template struct MFrontInterface<TRIDIMENSIONAL, AssemblyType::BLOCK_SCHUR>;
template struct MFrontInterface<AXISYMMETRICAL, AssemblyType::BLOCK_SCHUR>;
template struct MFrontInterface<PLANESTRAIN, AssemblyType::BLOCK_SCHUR>;

} // namespace MoFEM

#endif // WITH_MGIS
