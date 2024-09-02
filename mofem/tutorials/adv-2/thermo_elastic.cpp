/**
 * \file thermo_elastic.cpp
 * \example thermo_elastic.cpp
 *
 * Thermo plasticity
 *
 */

#ifndef EXECUTABLE_DIMENSION
#define EXECUTABLE_DIMENSION 3
#endif

#ifndef FINITE_DEFORMATION_FLAG
#define FINITE_DEFORMATION_FLAG true
#endif

#include <MoFEM.hpp>
#include <MatrixFunction.hpp>

using namespace MoFEM;

constexpr int SPACE_DIM =
    EXECUTABLE_DIMENSION; //< Space dimension of problem, mesh

constexpr bool IS_LARGE_STRAINS =
    FINITE_DEFORMATION_FLAG; //< Flag to turn off/on geometric nonlinearities

using EntData = EntitiesFieldData::EntData;
using DomainEle = PipelineManager::ElementsAndOpsByDim<SPACE_DIM>::DomainEle;
using DomainEleOp = DomainEle::UserDataOperator;
using BoundaryEle =
    PipelineManager::ElementsAndOpsByDim<SPACE_DIM>::BoundaryEle;
using BoundaryEleOp = BoundaryEle::UserDataOperator;
using PostProcEle = PostProcBrokenMeshInMoab<DomainEle>;

using AssemblyDomainEleOp =
    FormsIntegrators<DomainEleOp>::Assembly<PETSC>::OpBase;

//! [Default input parameters]
constexpr AssemblyType AT = AssemblyType::PETSC; //< selected assembly type
constexpr IntegrationType IT =
    IntegrationType::GAUSS; //< selected integration type

double young_modulus = 1;
double poisson_ratio = 0.25;
double ref_temp = 0.0;

double coeff_expansion = 1;
double heat_conductivity =
    1; // Force / (time temperature )  or Power /
       // (length temperature) // Time unit is hour. force unit kN
double heat_capacity = 1; // length^2/(time^2 temperature) // length is
                          // millimeter time is hour

int order = 2; //< default approximation order
//! [Default input parameters]

//! [Linear elastic problem]
#include <HenckyOps.hpp> // Include Hencky operators
using namespace HenckyOps;
//! [Linear elastic problem]

// Include finite deformation operators
#include <FiniteThermalOps.hpp> // Include operators for finite strain diffusion problem
using namespace FiniteThermalOps;

//! [Thermal problem]
/**
 * @brief Integrate Lhs base of flux (1/k) base of flux (FLUX x FLUX)
 *
 */

// Add alias to allow same implementation for small and large strains
using OpHdivHdiv = OpHdivHdivImpl<SPACE_DIM, IS_LARGE_STRAINS>;

/**
 * @brief Integrate Lhs div of base of flux times base of temperature (FLUX x
 * T) and transpose of it, i.e. (T x FLUX)
 *
 */
using OpHdivT = FormsIntegrators<DomainEleOp>::Assembly<PETSC>::BiLinearForm<
    GAUSS>::OpMixDivTimesScalar<SPACE_DIM>;

/**
 * @brief Integrate Lhs of flux term coupled to displacement field
 *
 */
// Templated on IS_LARGE_STRAINS to allow different implementations with
// different numbers of arguments at compile time. An empty operator is called
// for IS_LARGE_STRAINS = false
using OpHdivU = OpCalculateQdotQLhs_dU<SPACE_DIM, GAUSS, AssemblyDomainEleOp,
                                       IS_LARGE_STRAINS>;

/**
 * @brief Integrate Lhs base of temperature times (heat capacity) times base of
 * temperature (T x T)
 *
 */
using OpCapacity = FormsIntegrators<DomainEleOp>::Assembly<PETSC>::BiLinearForm<
    GAUSS>::OpMass<1, 1>;

/**
 * @brief Integrating Rhs flux base (1/k) flux  (FLUX)
 */

// Add alias to allow same implementation for small and large strains
using OpHdivFlux = OpHdivFluxImpl<SPACE_DIM, IS_LARGE_STRAINS>;

/**
 * @brief  Integrate Rhs div flux base times temperature (T)
 *
 */
using OpHDivTemp = FormsIntegrators<DomainEleOp>::Assembly<PETSC>::LinearForm<
    GAUSS>::OpMixDivTimesU<3, 1, SPACE_DIM>;

/**
 * @brief Integrate Rhs base of temperature time heat capacity times heat rate
 * (T)
 *
 */
using OpBaseDotT = FormsIntegrators<DomainEleOp>::Assembly<PETSC>::LinearForm<
    GAUSS>::OpBaseTimesScalar<1>;

/**
 * @brief Integrate Rhs base of temperature times divergence of flux (T)
 *
 */
using OpBaseDivFlux = OpBaseDotT;

//! [Thermal problem]

//! [Body and heat source]
using DomainNaturalBCRhs =
    NaturalBC<DomainEleOp>::Assembly<PETSC>::LinearForm<GAUSS>;
using OpBodyForce =
    DomainNaturalBCRhs::OpFlux<NaturalMeshsetType<BLOCKSET>, 1, SPACE_DIM>;
using OpHeatSource =
    DomainNaturalBCRhs::OpFlux<NaturalMeshsetType<BLOCKSET>, 1, 1>;
using DomainNaturalBCLhs =
    NaturalBC<DomainEleOp>::Assembly<PETSC>::BiLinearForm<GAUSS>;
//! [Body and heat source]

//! [Natural boundary conditions]
using BoundaryNaturalBC =
    NaturalBC<BoundaryEleOp>::Assembly<PETSC>::LinearForm<GAUSS>;
using OpForce = BoundaryNaturalBC::OpFlux<NaturalForceMeshsets, 1, SPACE_DIM>;
using OpTemperatureBC =
    BoundaryNaturalBC::OpFlux<NaturalTemperatureMeshsets, 3, SPACE_DIM>;
//! [Natural boundary conditions]

//! [Essential boundary conditions (Least square approach)]
using OpEssentialFluxRhs =
    EssentialBC<BoundaryEleOp>::Assembly<PETSC>::LinearForm<
        GAUSS>::OpEssentialRhs<HeatFluxCubitBcData, 3, SPACE_DIM>;
using OpEssentialFluxLhs =
    EssentialBC<BoundaryEleOp>::Assembly<PETSC>::BiLinearForm<
        GAUSS>::OpEssentialLhs<HeatFluxCubitBcData, 3, SPACE_DIM>;
//! [Essential boundary conditions (Least square approach)]

double default_young_modulus = 1;
double default_poisson_ratio = 0.25;
double ref_temp = 0.0;
double init_temp = 0.0;

double default_coeff_expansion = 1;
double default_heat_conductivity =
    1; // Force / (time temperature )  or Power /
       // (length temperature) // Time unit is hour. force unit kN
double default_heat_capacity = 1; // length^2/(time^2 temperature) // length is
                                  // millimeter time is hour

int order = 2; //< default approximation order

#include <ThermoElasticOps.hpp>   //< additional coupling operators
using namespace ThermoElasticOps; //< name space of coupling operators

using OpSetTemperatureRhs =
    DomainNaturalBCRhs::OpFlux<SetTargetTemperature, 1, 1>;
using OpSetTemperatureLhs =
    DomainNaturalBCLhs::OpFlux<SetTargetTemperature, 1, 1>;

auto save_range = [](moab::Interface &moab, const std::string name,
                     const Range r) {
  MoFEMFunctionBegin;
  auto out_meshset = get_temp_meshset_ptr(moab);
  CHKERR moab.add_entities(*out_meshset, r);
  CHKERR moab.write_file(name.c_str(), "VTK", "", out_meshset->get_ptr(), 1);
  MoFEMFunctionReturn(0);
};

struct ThermoElasticProblem {

  ThermoElasticProblem(MoFEM::Interface &m_field) : mField(m_field) {}

  MoFEMErrorCode runProblem();

private:
  MoFEM::Interface &mField;

  PetscBool doEvalField;
  std::array<double, SPACE_DIM> fieldEvalCoords;
  boost::shared_ptr<FieldEvaluatorInterface::SetPtsData> fieldEvalData;

  boost::shared_ptr<VectorDouble> scalarFieldPtr;
  boost::shared_ptr<MatrixDouble> vectorFieldPtr;
  boost::shared_ptr<MatrixDouble> tensorFieldPtr;

  MoFEMErrorCode setupProblem();     ///< add fields
  MoFEMErrorCode createCommonData(); //< read global data from command line
  MoFEMErrorCode bC();               //< read boundary conditions
  MoFEMErrorCode OPs();              //< add operators to pipeline
  MoFEMErrorCode tsSolve();          //< time solver

  template <int DIM, AssemblyType A, IntegrationType I, typename DomainEleOp>
  MoFEMErrorCode opThermoElasticFactoryDomainRhs(
      MoFEM::Interface &m_field,
      boost::ptr_deque<ForcesAndSourcesCore::UserDataOperator> &pip,
      std::string field_name,
      boost::shared_ptr<HenckyOps::CommonData> elastic_common_ptr,
      boost::shared_ptr<BlockedParameters> thermal_common_ptr, Sev sev) {
    MoFEMFunctionBegin;

    using B = typename FormsIntegrators<DomainEleOp>::template Assembly<
        A>::template LinearForm<I>;
    using H = HenckyOps::HenkyIntegrators<DomainEleOp>;
    auto vec_temp_ptr = boost::make_shared<VectorDouble>();
    pip.push_back(new OpCalculateScalarFieldValues("T", vec_temp_ptr));
    auto coeff_expansion_ptr = thermal_common_ptr->getCoeffExpansionPtr();
    auto ref_temp_ptr = thermal_common_ptr->getRefTempPtr();
    pip.push_back(new
                  typename H::template OpCalculateHenckyThermalStress<DIM, I>(
                      "U", vec_temp_ptr, elastic_common_ptr,
                      coeff_expansion_ptr, ref_temp_ptr));
    pip.push_back(new typename H::template OpCalculatePiolaStress<DIM, I>(
        field_name, elastic_common_ptr));
    using OpInternalForcePiola =
        typename B::template OpGradTimesTensor<1, DIM, DIM>;
    pip.push_back(new OpInternalForcePiola(
        "U", elastic_common_ptr->getMatFirstPiolaStress()));

    MoFEMFunctionReturn(0);
  }

  template <int DIM, AssemblyType A, IntegrationType I, typename DomainEleOp>
  MoFEMErrorCode opThermoElasticFactoryDomainRhs(
      MoFEM::Interface &m_field,
      boost::ptr_deque<ForcesAndSourcesCore::UserDataOperator> &pip,
      std::string field_name, std::string elastic_block_name,
      std::string thermal_block_name, Sev sev, double scale = 1) {
    MoFEMFunctionBegin;

    auto elastic_common_ptr = commonDataFactory<DIM, I, DomainEleOp>(
        m_field, pip, field_name, elastic_block_name, sev, scale);
    auto thermal_common_ptr = boost::make_shared<BlockedParameters>();
    CHKERR addMatThermalBlockOps(m_field, pip, thermal_block_name,
                                 thermal_common_ptr, Sev::inform);
    CHKERR opThermoElasticFactoryDomainRhs<DIM, A, I, DomainEleOp>(
        m_field, pip, field_name, elastic_common_ptr, thermal_common_ptr, sev);

    MoFEMFunctionReturn(0);
  }

  template <int DIM, AssemblyType A, IntegrationType I, typename DomainEleOp>
  MoFEMErrorCode opThermoElasticFactoryDomainLhs(
      MoFEM::Interface &m_field,
      boost::ptr_deque<ForcesAndSourcesCore::UserDataOperator> &pip,
      std::string field_name, std::string coupled_field_name,
      boost::shared_ptr<HenckyOps::CommonData> elastic_common_ptr,
      boost::shared_ptr<BlockedParameters> thermal_common_ptr, Sev sev) {
    MoFEMFunctionBegin;

    using B = typename FormsIntegrators<DomainEleOp>::template Assembly<
        A>::template BiLinearForm<I>;
    using OpKPiola = typename B::template OpGradTensorGrad<1, DIM, DIM, 1>;

    using H = HenkyIntegrators<DomainEleOp>;
    auto vec_temp_ptr = boost::make_shared<VectorDouble>();
    pip.push_back(new OpCalculateScalarFieldValues("T", vec_temp_ptr));
    auto coeff_expansion_ptr = thermal_common_ptr->getCoeffExpansionPtr();
    auto ref_temp_ptr = thermal_common_ptr->getRefTempPtr();
    pip.push_back(new
                  typename H::template OpCalculateHenckyThermalStress<DIM, I>(
                      "U", vec_temp_ptr, elastic_common_ptr,
                      coeff_expansion_ptr, ref_temp_ptr));
    pip.push_back(new typename H::template OpCalculatePiolaStress<DIM, I>(
        field_name, elastic_common_ptr));
    pip.push_back(new typename H::template OpHenckyTangent<DIM, I>(
        field_name, elastic_common_ptr));
    pip.push_back(new OpKPiola(field_name, field_name,
                               elastic_common_ptr->getMatTangent()));
    pip.push_back(new typename HenckyOps::OpCalculateHenckyThermalStressdT<
                  DIM, I, AssemblyDomainEleOp>(field_name, coupled_field_name,
                                               elastic_common_ptr,
                                               coeff_expansion_ptr));

    MoFEMFunctionReturn(0);
  }

  template <int DIM, AssemblyType A, IntegrationType I, typename DomainEleOp>
  MoFEMErrorCode opThermoElasticFactoryDomainLhs(
      MoFEM::Interface &m_field,
      boost::ptr_deque<ForcesAndSourcesCore::UserDataOperator> &pip,
      std::string field_name, std::string coupled_field_name,
      std::string elastic_block_name, std::string thermal_block_name, Sev sev,
      double scale = 1) {
    MoFEMFunctionBegin;

    auto elastic_common_ptr = commonDataFactory<DIM, I, DomainEleOp>(
        m_field, pip, field_name, elastic_block_name, sev, scale);
    auto thermal_common_ptr = boost::make_shared<BlockedParameters>();
    CHKERR addMatThermalBlockOps(m_field, pip, thermal_block_name,
                                 thermal_common_ptr, Sev::inform);
    CHKERR opThermoElasticFactoryDomainLhs<DIM, A, I, DomainEleOp>(
        m_field, pip, field_name, coupled_field_name, elastic_common_ptr,
        thermal_common_ptr, sev);

    MoFEMFunctionReturn(0);
  }
};

//! [Run problem]
MoFEMErrorCode ThermoElasticProblem::runProblem() {
  MoFEMFunctionBegin;
  CHKERR setupProblem();
  CHKERR createCommonData();
  CHKERR bC();
  CHKERR OPs();
  CHKERR tsSolve();
  MoFEMFunctionReturn(0);
}
//! [Run problem]

//! [Set up problem]
MoFEMErrorCode ThermoElasticProblem::setupProblem() {
  MoFEMFunctionBegin;
  Simple *simple = mField.getInterface<Simple>();
  // Add field
  constexpr FieldApproximationBase base = AINSWORTH_LEGENDRE_BASE;
  // Mechanical fields
  CHKERR simple->addDomainField("U", H1, base, SPACE_DIM);
  CHKERR simple->addBoundaryField("U", H1, base, SPACE_DIM);
  // Temperature
  constexpr auto flux_space = (SPACE_DIM == 2) ? HCURL : HDIV;
  CHKERR simple->addDomainField("T", L2, AINSWORTH_LEGENDRE_BASE, 1);
  CHKERR simple->addDomainField("FLUX", flux_space, DEMKOWICZ_JACOBI_BASE, 1);
  CHKERR simple->addBoundaryField("FLUX", flux_space, DEMKOWICZ_JACOBI_BASE, 1);

  CHKERR PetscOptionsGetInt(PETSC_NULL, "", "-order", &order, PETSC_NULL);
  CHKERR simple->setFieldOrder("U", order + 1);
  CHKERR simple->setFieldOrder("FLUX", order + 1);
  CHKERR simple->setFieldOrder("T", order);
  CHKERR simple->setUp();

  int coords_dim = SPACE_DIM;
  CHKERR PetscOptionsGetRealArray(NULL, NULL, "-field_eval_coords",
                                  fieldEvalCoords.data(), &coords_dim,
                                  &doEvalField);

  scalarFieldPtr = boost::make_shared<VectorDouble>();
  vectorFieldPtr = boost::make_shared<MatrixDouble>();
  tensorFieldPtr = boost::make_shared<MatrixDouble>();

  if (doEvalField) {
    fieldEvalData =
        mField.getInterface<FieldEvaluatorInterface>()->getData<DomainEle>();

    if constexpr (SPACE_DIM == 3) {
      CHKERR mField.getInterface<FieldEvaluatorInterface>()->buildTree3D(
          fieldEvalData, simple->getDomainFEName());
    } else {
      CHKERR mField.getInterface<FieldEvaluatorInterface>()->buildTree2D(
          fieldEvalData, simple->getDomainFEName());
    }

    fieldEvalData->setEvalPoints(fieldEvalCoords.data(), 1);
    auto no_rule = [](int, int, int) { return -1; };

    auto field_eval_fe_ptr = fieldEvalData->feMethodPtr.lock();
    field_eval_fe_ptr->getRuleHook = no_rule;

    field_eval_fe_ptr->getOpPtrVector().push_back(
        new OpCalculateScalarFieldValues("T", scalarFieldPtr));
    field_eval_fe_ptr->getOpPtrVector().push_back(
        new OpCalculateVectorFieldValues<SPACE_DIM>("U", vectorFieldPtr));
    field_eval_fe_ptr->getOpPtrVector().push_back(
        new OpCalculateVectorFieldGradient<SPACE_DIM, SPACE_DIM>(
            "U", tensorFieldPtr));
  }

  MoFEMFunctionReturn(0);
}
//! [Set up problem]

//! [Create common data]
MoFEMErrorCode ThermoElasticProblem::createCommonData() {
  MoFEMFunctionBegin;

  auto get_command_line_parameters = [&]() {
    MoFEMFunctionBegin;
    CHKERR PetscOptionsGetScalar(PETSC_NULL, "", "-young_modulus",
                                 &young_modulus, PETSC_NULL);
    CHKERR PetscOptionsGetScalar(PETSC_NULL, "", "-poisson_ratio",
                                 &poisson_ratio, PETSC_NULL);
    CHKERR PetscOptionsGetScalar(PETSC_NULL, "", "-coeff_expansion",
                                 &coeff_expansion, PETSC_NULL);
    CHKERR PetscOptionsGetScalar(PETSC_NULL, "", "-ref_temp", &ref_temp,
                                 PETSC_NULL);
    CHKERR PetscOptionsGetScalar(PETSC_NULL, "", "-init_temp", &init_temp,
                                 PETSC_NULL);

    CHKERR PetscOptionsGetScalar(PETSC_NULL, "", "-capacity", &heat_capacity,
                                 PETSC_NULL);
    CHKERR PetscOptionsGetScalar(PETSC_NULL, "", "-conductivity",
                                 &heat_conductivity, PETSC_NULL);

    MOFEM_LOG("ThermoElastic", Sev::inform)
        << "Young modulus " << young_modulus;
    MOFEM_LOG("ThermoElastic", Sev::inform)
        << "Poisson ratio " << poisson_ratio;
    MOFEM_LOG("ThermoElastic", Sev::inform)
        << "Coeff of expansion " << coeff_expansion;
    MOFEM_LOG("ThermoElastic", Sev::inform) << "Capacity " << heat_capacity;
    MOFEM_LOG("ThermoElastic", Sev::inform)
        << "Heat conductivity " << heat_conductivity;

    MOFEM_LOG("ThermoElastic", Sev::inform)
        << "Reference temperature  " << ref_temp;
    MOFEM_LOG("ThermoElastic", Sev::inform)
        << "Initial temperature  " << init_temp;

    MoFEMFunctionReturn(0);
  };

  CHKERR get_command_line_parameters();

  MoFEMFunctionReturn(0);
}
//! [Create common data]

//! [Boundary condition]
MoFEMErrorCode ThermoElasticProblem::bC() {
  MoFEMFunctionBegin;

  MOFEM_LOG("SYNC", Sev::noisy) << "bC";
  MOFEM_LOG_SEVERITY_SYNC(mField.get_comm(), Sev::noisy);

  auto simple = mField.getInterface<Simple>();
  auto bc_mng = mField.getInterface<BcManager>();



  auto get_skin = [&]() {
    Range body_ents;
    CHKERR mField.get_moab().get_entities_by_dimension(0, SPACE_DIM, body_ents);
    Skinner skin(&mField.get_moab());
    Range skin_ents;
    CHKERR skin.find_skin(0, body_ents, false, skin_ents);
    return skin_ents;
  };

  auto filter_flux_blocks = [&](auto skin) {
    auto remove_cubit_blocks = [&](auto c) {
      MoFEMFunctionBegin;
      for (auto m :

           mField.getInterface<MeshsetsManager>()->getCubitMeshsetPtr(c)

      ) {
        Range ents;
        CHKERR mField.get_moab().get_entities_by_dimension(
            m->getMeshset(), SPACE_DIM - 1, ents, true);
        skin = subtract(skin, ents);
      }
      MoFEMFunctionReturn(0);
    };

    auto remove_named_blocks = [&](auto n) {
      MoFEMFunctionBegin;
      for (auto m : mField.getInterface<MeshsetsManager>()->getCubitMeshsetPtr(
               std::regex(

                   (boost::format("%s(.*)") % n).str()

                       ))

      ) {
        Range ents;
        CHKERR mField.get_moab().get_entities_by_dimension(
            m->getMeshset(), SPACE_DIM - 1, ents, true);
        skin = subtract(skin, ents);
      }
      MoFEMFunctionReturn(0);
    };

    CHK_THROW_MESSAGE(remove_cubit_blocks(NODESET | TEMPERATURESET),
                      "remove_cubit_blocks");
    CHK_THROW_MESSAGE(remove_cubit_blocks(SIDESET | HEATFLUXSET),
                      "remove_cubit_blocks");
    CHK_THROW_MESSAGE(remove_named_blocks("TEMPERATURE"),
                      "remove_named_blocks");
    CHK_THROW_MESSAGE(remove_named_blocks("HEATFLUX"), "remove_named_blocks");

    return skin;
  };

  auto filter_true_skin = [&](auto skin) {
    Range boundary_ents;
    ParallelComm *pcomm =
        ParallelComm::get_pcomm(&mField.get_moab(), MYPCOMM_INDEX);
    CHKERR pcomm->filter_pstatus(skin, PSTATUS_SHARED | PSTATUS_MULTISHARED,
                                 PSTATUS_NOT, -1, &boundary_ents);
    return boundary_ents;
  };

  auto remove_flux_ents = filter_true_skin(filter_flux_blocks(get_skin()));

  CHKERR mField.getInterface<CommInterface>()->synchroniseEntities(
      remove_flux_ents);

  MOFEM_LOG("SYNC", Sev::noisy) << remove_flux_ents << endl;
  MOFEM_LOG_SEVERITY_SYNC(mField.get_comm(), Sev::noisy);

#ifndef NDEBUG

  CHKERR save_range(
      mField.get_moab(),
      (boost::format("flux_remove_%d.vtk") % mField.get_comm_rank()).str(),
      remove_flux_ents);

#endif

  CHKERR mField.getInterface<ProblemsManager>()->removeDofsOnEntities(
      simple->getProblemName(), "FLUX", remove_flux_ents);

  auto set_init_temp = [](boost::shared_ptr<FieldEntity> field_entity_ptr) {
    field_entity_ptr->getEntFieldData()[0] = init_temp;
    return 0;
  };
  CHKERR mField.getInterface<FieldBlas>()->fieldLambdaOnEntities(set_init_temp,
                                                                 "T");

  CHKERR bc_mng->removeBlockDOFsOnEntities<DisplacementCubitBcData>(
      simple->getProblemName(), "U");
  CHKERR bc_mng->pushMarkDOFsOnEntities<HeatFluxCubitBcData>(
      simple->getProblemName(), "FLUX", false);

  MoFEMFunctionReturn(0);
}
//! [Boundary condition]

//! [Push operators to pipeline]
MoFEMErrorCode ThermoElasticProblem::OPs() {
  MoFEMFunctionBegin;

  MOFEM_LOG("SYNC", Sev::noisy) << "OPs";
  MOFEM_LOG_SEVERITY_SYNC(mField.get_comm(), Sev::noisy);

  auto pipeline_mng = mField.getInterface<PipelineManager>();
  auto simple = mField.getInterface<Simple>();
  auto bc_mng = mField.getInterface<BcManager>();

  auto boundary_marker =
      bc_mng->getMergedBlocksMarker(vector<string>{"HEATFLUX"});

  auto integration_rule = [](int, int, int approx_order) {
    return 2 * approx_order;
  };
  CHKERR pipeline_mng->setDomainRhsIntegrationRule(integration_rule);
  CHKERR pipeline_mng->setDomainLhsIntegrationRule(integration_rule);
  CHKERR pipeline_mng->setBoundaryRhsIntegrationRule(integration_rule);
  CHKERR pipeline_mng->setBoundaryLhsIntegrationRule(integration_rule);

  auto block_params = boost::make_shared<BlockedParameters>();
  auto coeff_expansion_ptr = block_params->getCoeffExpansionPtr();
  auto ref_temp_ptr = block_params->getRefTempPtr();
  auto heat_conductivity_ptr = block_params->getHeatConductivityPtr();
  auto heat_capacity_ptr = block_params->getHeatCapacityPtr();

  // Default time scaling of BCs and sources, from command line arguments
  auto time_scale = boost::make_shared<TimeScale>();

  // Files which define scaling for separate variables, if provided
  auto time_bodyforce_scaling =
      boost::make_shared<TimeScale>("bodyforce_scale.txt");
  auto time_heatsource_scaling =
      boost::make_shared<TimeScale>("heatsource_scale.txt");
  auto time_temperature_scaling =
      boost::make_shared<TimeScale>("temperature_bc_scale.txt");
  auto time_displacement_scaling =
      boost::make_shared<TimeScale>("displacement_bc_scale.txt");
  auto time_flux_scaling = boost::make_shared<TimeScale>("flux_bc_scale.txt");
  auto time_force_scaling = boost::make_shared<TimeScale>("force_bc_scale.txt");

  auto add_domain_rhs_ops = [&](auto &pipeline) {
    MoFEMFunctionBegin;
    CHKERR addMatThermalBlockOps(mField, pipeline, "MAT_THERMAL", block_params,
                                 Sev::inform);
    CHKERR AddHOOps<SPACE_DIM, SPACE_DIM, SPACE_DIM>::add(pipeline, {H1, HDIV});

    auto hencky_common_data_ptr =
        HenckyOps::commonDataFactory<SPACE_DIM, IT, DomainEleOp>(
            mField, pipeline, "U", "MAT_ELASTIC", Sev::inform);
    auto mat_D_ptr = hencky_common_data_ptr->matDPtr;
    auto mat_grad_ptr = hencky_common_data_ptr->matGradPtr;
    auto mat_strain_ptr = boost::make_shared<MatrixDouble>();
    auto mat_stress_ptr = boost::make_shared<MatrixDouble>();

    auto vec_temp_ptr = boost::make_shared<VectorDouble>();
    auto vec_temp_dot_ptr = boost::make_shared<VectorDouble>();
    auto mat_flux_ptr = boost::make_shared<MatrixDouble>();
    auto vec_temp_div_ptr = boost::make_shared<VectorDouble>();

    pipeline.push_back(new OpCalculateScalarFieldValues("T", vec_temp_ptr));
    pipeline.push_back(
        new OpCalculateScalarFieldValuesDot("T", vec_temp_dot_ptr));
    pipeline.push_back(new OpCalculateHdivVectorDivergence<3, SPACE_DIM>(
        "FLUX", vec_temp_div_ptr));
    pipeline.push_back(
        new OpCalculateHVecVectorField<3, SPACE_DIM>("FLUX", mat_flux_ptr));

    CHKERR opThermoElasticFactoryDomainRhs<SPACE_DIM, AT, IT, DomainEleOp>(
        mField, pipeline, "U", "MAT_ELASTIC", "MAT_THERMAL", Sev::inform);

    pipeline.push_back(new OpSetBc("FLUX", true, boundary_marker));

    auto resistance = [heat_conductivity_ptr](const double, const double,
                                              const double) {
      return (1. / (*heat_conductivity_ptr));
    };
    // negative value is consequence of symmetric system
    auto capacity = [heat_capacity_ptr](const double, const double,
                                        const double) {
      return -(*heat_capacity_ptr);
    };
    auto unity = [](const double, const double, const double) constexpr {
      return -1.;
    };
    pipeline.push_back(
        new OpHdivFlux("FLUX", mat_flux_ptr, resistance, mat_grad_ptr));
    pipeline.push_back(new OpHDivTemp("FLUX", vec_temp_ptr, unity));
    pipeline.push_back(new OpBaseDivFlux("T", vec_temp_div_ptr, unity));
    pipeline.push_back(new OpBaseDotT("T", vec_temp_dot_ptr, capacity));

    pipeline.push_back(new OpUnSetBc("FLUX"));

    // Set source terms
    CHKERR DomainNaturalBCRhs::AddFluxToPipeline<OpHeatSource>::add(
        pipeline, mField, "T", {time_scale, time_heatsource_scaling},
        "HEAT_SOURCE", Sev::inform);
    CHKERR DomainNaturalBCRhs::AddFluxToPipeline<OpBodyForce>::add(
        pipeline, mField, "U", {time_scale, time_bodyforce_scaling},
        "BODY_FORCE", Sev::inform);
    CHKERR DomainNaturalBCRhs::AddFluxToPipeline<OpSetTemperatureRhs>::add(
        pipeline, mField, "T", vec_temp_ptr, "SETTEMP", Sev::inform);

    MoFEMFunctionReturn(0);
  };

  auto add_domain_lhs_ops = [&](auto &pipeline) {
    MoFEMFunctionBegin;
    CHKERR addMatThermalBlockOps(mField, pipeline, "MAT_THERMAL", block_params,
                                 Sev::verbose);
    CHKERR AddHOOps<SPACE_DIM, SPACE_DIM, SPACE_DIM>::add(pipeline, {H1, HDIV});

    auto hencky_common_data_ptr =
        HenckyOps::commonDataFactory<SPACE_DIM, IT, DomainEleOp>(
            mField, pipeline, "U", "MAT_ELASTIC", Sev::inform, 1);
    auto mat_D_ptr = hencky_common_data_ptr->matDPtr;
    auto mat_grad_ptr = hencky_common_data_ptr->matGradPtr;

    pipeline.push_back(new OpSetBc("FLUX", true, boundary_marker));

    auto resistance = [heat_conductivity_ptr](const double, const double,
                                              const double) {
      return (1. / (*heat_conductivity_ptr));
    };
    auto capacity = [heat_capacity_ptr](const double, const double,
                                        const double) {
      return -(*heat_capacity_ptr);
    };

    pipeline.push_back(
        new OpHdivHdiv("FLUX", "FLUX", resistance, mat_grad_ptr));
    pipeline.push_back(
        new OpHdivT("FLUX", "T", []() constexpr { return -1; }, true));

    auto mat_flux_ptr = boost::make_shared<MatrixDouble>();
    pipeline.push_back(
        new OpCalculateHVecVectorField<3, SPACE_DIM>("FLUX", mat_flux_ptr));

    pipeline.push_back(
        new OpHdivU("FLUX", "U", mat_flux_ptr, resistance, mat_grad_ptr));

    CHKERR opThermoElasticFactoryDomainLhs<SPACE_DIM, AT, IT, DomainEleOp>(
        mField, pipeline, "U", "T", "MAT_ELASTIC", "MAT_THERMAL", Sev::inform);

    auto op_capacity = new OpCapacity("T", "T", capacity);
    op_capacity->feScalingFun = [](const FEMethod *fe_ptr) {
      return fe_ptr->ts_a;
    };
    pipeline.push_back(op_capacity);

    pipeline.push_back(new OpUnSetBc("FLUX"));

    auto vec_temp_ptr = boost::make_shared<VectorDouble>();
    pipeline.push_back(new OpCalculateScalarFieldValues("T", vec_temp_ptr));
    CHKERR DomainNaturalBCLhs::AddFluxToPipeline<OpSetTemperatureLhs>::add(
        pipeline, mField, "T", vec_temp_ptr, "SETTEMP", Sev::verbose);

    MoFEMFunctionReturn(0);
  };

  auto add_boundary_rhs_ops = [&](auto &pipeline) {
    MoFEMFunctionBegin;

    CHKERR AddHOOps<SPACE_DIM - 1, SPACE_DIM, SPACE_DIM>::add(pipeline, {HDIV});

    pipeline.push_back(new OpSetBc("FLUX", true, boundary_marker));

    // Set BCs using the least squares imposition
    CHKERR BoundaryNaturalBC::AddFluxToPipeline<OpForce>::add(
        pipeline_mng->getOpBoundaryRhsPipeline(), mField, "U",
        {time_scale, time_force_scaling}, "FORCE", Sev::inform);

    CHKERR BoundaryNaturalBC::AddFluxToPipeline<OpTemperatureBC>::add(
        pipeline_mng->getOpBoundaryRhsPipeline(), mField, "FLUX",
        {time_scale, time_temperature_scaling}, "TEMPERATURE", Sev::inform);

    pipeline.push_back(new OpUnSetBc("FLUX"));

    auto mat_flux_ptr = boost::make_shared<MatrixDouble>();
    pipeline.push_back(
        new OpCalculateHVecVectorField<3, SPACE_DIM>("FLUX", mat_flux_ptr));
    CHKERR EssentialBC<BoundaryEleOp>::Assembly<PETSC>::LinearForm<GAUSS>::
        AddEssentialToPipeline<OpEssentialFluxRhs>::add(
            mField, pipeline, simple->getProblemName(), "FLUX", mat_flux_ptr,
            {time_scale, time_flux_scaling});

    MoFEMFunctionReturn(0);
  };

  auto add_boundary_lhs_ops = [&](auto &pipeline) {
    MoFEMFunctionBegin;

    CHKERR AddHOOps<SPACE_DIM - 1, SPACE_DIM, SPACE_DIM>::add(pipeline, {HDIV});

    CHKERR EssentialBC<BoundaryEleOp>::Assembly<PETSC>::BiLinearForm<GAUSS>::
        AddEssentialToPipeline<OpEssentialFluxLhs>::add(
            mField, pipeline, simple->getProblemName(), "FLUX");

    MoFEMFunctionReturn(0);
  };

  // Set BCs by eliminating degrees of freedom
  auto get_bc_hook_rhs = [&]() {
    EssentialPreProc<DisplacementCubitBcData> hook(
        mField, pipeline_mng->getDomainRhsFE(),
        {time_scale, time_displacement_scaling});
    return hook;
  };
  auto get_bc_hook_lhs = [&]() {
    EssentialPreProc<DisplacementCubitBcData> hook(
        mField, pipeline_mng->getDomainLhsFE(),
        {time_scale, time_displacement_scaling});
    return hook;
  };

  pipeline_mng->getDomainRhsFE()->preProcessHook = get_bc_hook_rhs();
  pipeline_mng->getDomainLhsFE()->preProcessHook = get_bc_hook_lhs();

  CHKERR add_domain_rhs_ops(pipeline_mng->getOpDomainRhsPipeline());
  CHKERR add_domain_lhs_ops(pipeline_mng->getOpDomainLhsPipeline());
  CHKERR add_boundary_rhs_ops(pipeline_mng->getOpBoundaryRhsPipeline());
  CHKERR add_boundary_lhs_ops(pipeline_mng->getOpBoundaryLhsPipeline());

  MoFEMFunctionReturn(0);
}
//! [Push operators to pipeline]

//! [Solve]
MoFEMErrorCode ThermoElasticProblem::tsSolve() {
  MoFEMFunctionBegin;

  Simple *simple = mField.getInterface<Simple>();
  PipelineManager *pipeline_mng = mField.getInterface<PipelineManager>();
  ISManager *is_manager = mField.getInterface<ISManager>();

  auto dm = simple->getDM();
  auto solver = pipeline_mng->createTSIM();
  auto snes_ctx_ptr = getDMSnesCtx(dm);

  auto set_section_monitor = [&](auto solver) {
    MoFEMFunctionBegin;
    SNES snes;
    CHKERR TSGetSNES(solver, &snes);
    CHKERR SNESMonitorSet(snes,
                          (MoFEMErrorCode(*)(SNES, PetscInt, PetscReal,
                                             void *))MoFEMSNESMonitorFields,
                          (void *)(snes_ctx_ptr.get()), nullptr);
    MoFEMFunctionReturn(0);
  };

  auto create_post_process_element = [&]() {
    auto post_proc_fe = boost::make_shared<PostProcEle>(mField);

    auto block_params = boost::make_shared<BlockedParameters>();
    auto coeff_expansion_ptr = block_params->getCoeffExpansionPtr();
    auto ref_temp_ptr = block_params->getRefTempPtr();

    CHKERR addMatThermalBlockOps(mField, post_proc_fe->getOpPtrVector(),
                                 "MAT_THERMAL", block_params, Sev::verbose);

    CHKERR AddHOOps<SPACE_DIM, SPACE_DIM, SPACE_DIM>::add(
        post_proc_fe->getOpPtrVector(), {H1, HDIV});

    auto mat_grad_ptr = boost::make_shared<MatrixDouble>();
    auto mat_strain_ptr = boost::make_shared<MatrixDouble>();
    auto mat_stress_ptr = boost::make_shared<MatrixDouble>();

    auto vec_temp_ptr = boost::make_shared<VectorDouble>();
    auto mat_flux_ptr = boost::make_shared<MatrixDouble>();

    post_proc_fe->getOpPtrVector().push_back(
        new OpCalculateScalarFieldValues("T", vec_temp_ptr));
    post_proc_fe->getOpPtrVector().push_back(
        new OpCalculateHVecVectorField<3, SPACE_DIM>("FLUX", mat_flux_ptr));

    auto u_ptr = boost::make_shared<MatrixDouble>();
    post_proc_fe->getOpPtrVector().push_back(
        new OpCalculateVectorFieldValues<SPACE_DIM>("U", u_ptr));
    post_proc_fe->getOpPtrVector().push_back(
        new OpCalculateVectorFieldGradient<SPACE_DIM, SPACE_DIM>("U",
                                                                 mat_grad_ptr));
    auto elastic_common_ptr = commonDataFactory<SPACE_DIM, IT, DomainEleOp>(
        mField, post_proc_fe->getOpPtrVector(), "U", "MAT_ELASTIC",
        Sev::inform);
    using H = HenckyOps::HenkyIntegrators<DomainEleOp>;
    post_proc_fe->getOpPtrVector().push_back(
        new typename H::template OpCalculateHenckyThermalStress<SPACE_DIM, IT>(
            "U", vec_temp_ptr, elastic_common_ptr, coeff_expansion_ptr,
            ref_temp_ptr));
    post_proc_fe->getOpPtrVector().push_back(
        new typename H::template OpCalculatePiolaStress<SPACE_DIM, IT>(
            "U", elastic_common_ptr));

    using OpPPMap = OpPostProcMapInMoab<SPACE_DIM, SPACE_DIM>;

    post_proc_fe->getOpPtrVector().push_back(

        new OpPPMap(

            post_proc_fe->getPostProcMesh(), post_proc_fe->getMapGaussPts(),

            {{"T", vec_temp_ptr}},

            {{"U", u_ptr}, {"FLUX", mat_flux_ptr}},

            {{"PIOLA", elastic_common_ptr->getMatFirstPiolaStress()}},

            {{"HENCKY_STRAIN", elastic_common_ptr->getMatLogC()}}

            )

    );

    return post_proc_fe;
  };

  auto monitor_ptr = boost::make_shared<FEMethod>();
  auto post_proc_fe = create_post_process_element();

  auto set_time_monitor = [&](auto dm, auto solver) {
    MoFEMFunctionBegin;
    monitor_ptr->preProcessHook = [&]() {
      MoFEMFunctionBegin;

      CHKERR DMoFEMLoopFiniteElements(dm, simple->getDomainFEName(),
                                      post_proc_fe,
                                      monitor_ptr->getCacheWeakPtr());
      CHKERR post_proc_fe->writeFile(
          "out_" + boost::lexical_cast<std::string>(monitor_ptr->ts_step) +
          ".h5m");

      if (doEvalField) {
        if constexpr (SPACE_DIM == 3) {
          CHKERR mField.getInterface<FieldEvaluatorInterface>()
              ->evalFEAtThePoint3D(
                  fieldEvalCoords.data(), 1e-12, simple->getProblemName(),
                  simple->getDomainFEName(), fieldEvalData,
                  mField.get_comm_rank(), mField.get_comm_rank(), nullptr,
                  MF_EXIST, QUIET);
        } else {
          CHKERR mField.getInterface<FieldEvaluatorInterface>()
              ->evalFEAtThePoint2D(
                  fieldEvalCoords.data(), 1e-12, simple->getProblemName(),
                  simple->getDomainFEName(), fieldEvalData,
                  mField.get_comm_rank(), mField.get_comm_rank(), nullptr,
                  MF_EXIST, QUIET);
        }

        if (scalarFieldPtr->size()) {
          auto t_temp = getFTensor0FromVec(*scalarFieldPtr);
          MOFEM_LOG("ThermoElasticSync", Sev::inform)
              << "Eval point T: " << t_temp;
        }
        if (vectorFieldPtr->size1()) {
          FTensor::Index<'i', SPACE_DIM> i;
          auto t_disp = getFTensor1FromMat<SPACE_DIM>(*vectorFieldPtr);
          MOFEM_LOG("ThermoElasticSync", Sev::inform)
              << "Eval point U magnitude: " << sqrt(t_disp(i) * t_disp(i));
        }
        if (tensorFieldPtr->size1()) {
          FTensor::Index<'i', SPACE_DIM> i;
          auto t_disp_grad =
              getFTensor2FromMat<SPACE_DIM, SPACE_DIM>(*tensorFieldPtr);
          MOFEM_LOG("ThermoElasticSync", Sev::inform)
              << "Eval point U_GRAD trace: " << t_disp_grad(i, i);
        }

        MOFEM_LOG_SYNCHRONISE(mField.get_comm());
      }

      MoFEMFunctionReturn(0);
    };
    auto null = boost::shared_ptr<FEMethod>();
    CHKERR DMMoFEMTSSetMonitor(dm, solver, simple->getDomainFEName(), null,
                               monitor_ptr, null);
    MoFEMFunctionReturn(0);
  };

  auto set_fieldsplit_preconditioner = [&](auto solver) {
    MoFEMFunctionBeginHot;

    SNES snes;
    CHKERR TSGetSNES(solver, &snes);
    KSP ksp;
    CHKERR SNESGetKSP(snes, &ksp);
    PC pc;
    CHKERR KSPGetPC(ksp, &pc);
    PetscBool is_pcfs = PETSC_FALSE;
    PetscObjectTypeCompare((PetscObject)pc, PCFIELDSPLIT, &is_pcfs);

    // Setup fieldsplit (block) solver - optional: yes/no
    if (is_pcfs == PETSC_TRUE) {
      auto bc_mng = mField.getInterface<BcManager>();
      auto is_mng = mField.getInterface<ISManager>();
      auto name_prb = simple->getProblemName();

      SmartPetscObj<IS> is_u;
      CHKERR is_mng->isCreateProblemFieldAndRank(name_prb, ROW, "U", 0,
                                                 SPACE_DIM, is_u);
      SmartPetscObj<IS> is_flux;
      CHKERR is_mng->isCreateProblemFieldAndRank(name_prb, ROW, "FLUX", 0, 0,
                                                 is_flux);
      SmartPetscObj<IS> is_T;
      CHKERR is_mng->isCreateProblemFieldAndRank(name_prb, ROW, "T", 0, 0,
                                                 is_T);
      IS is_tmp;
      CHKERR ISExpand(is_T, is_flux, &is_tmp);
      auto is_TFlux = SmartPetscObj<IS>(is_tmp);

      auto is_all_bc = bc_mng->getBlockIS(name_prb, "HEATFLUX", "FLUX", 0, 0);
      int is_all_bc_size;
      CHKERR ISGetSize(is_all_bc, &is_all_bc_size);
      MOFEM_LOG("ThermoElastic", Sev::inform)
          << "Field split block size " << is_all_bc_size;
      if (is_all_bc_size) {
        IS is_tmp2;
        CHKERR ISDifference(is_TFlux, is_all_bc, &is_tmp2);
        is_TFlux = SmartPetscObj<IS>(is_tmp2);
        CHKERR PCFieldSplitSetIS(pc, PETSC_NULL,
                                 is_all_bc); // boundary block
      }

      CHKERR ISSort(is_u);
      CHKERR ISSort(is_TFlux);
      CHKERR PCFieldSplitSetIS(pc, PETSC_NULL, is_TFlux);
      CHKERR PCFieldSplitSetIS(pc, PETSC_NULL, is_u);
    }

    MoFEMFunctionReturnHot(0);
  };

  auto D = createDMVector(dm);
  CHKERR DMoFEMMeshToLocalVector(dm, D, INSERT_VALUES, SCATTER_FORWARD);
  CHKERR TSSetSolution(solver, D);
  CHKERR TSSetFromOptions(solver);

  CHKERR set_section_monitor(solver);
  CHKERR set_fieldsplit_preconditioner(solver);
  CHKERR set_time_monitor(dm, solver);

  CHKERR TSSetUp(solver);
  CHKERR TSSolve(solver, NULL);

  MoFEMFunctionReturn(0);
}
//! [Solve]

static char help[] = "...\n\n";

int main(int argc, char *argv[]) {

  // Initialisation of MoFEM/PETSc and MOAB data structures
  const char param_file[] = "param_file.petsc";
  MoFEM::Core::Initialize(&argc, &argv, param_file, help);

  // Add logging channel for example
  auto core_log = logging::core::get();
  core_log->add_sink(
      LogManager::createSink(LogManager::getStrmWorld(), "ThermoElastic"));
  LogManager::setLog("ThermoElastic");
  MOFEM_LOG_TAG("ThermoElastic", "ThermoElastic");

  core_log->add_sink(
      LogManager::createSink(LogManager::getStrmSync(), "ThermoElasticSync"));
  LogManager::setLog("ThermoElasticSync");
  MOFEM_LOG_TAG("ThermoElasticSync", "ThermoElasticSync");

  try {

    //! [Register MoFEM discrete manager in PETSc]
    DMType dm_name = "DMMOFEM";
    CHKERR DMRegister_MoFEM(dm_name);
    //! [Register MoFEM discrete manager in PETSc

    //! [Create MoAB]
    moab::Core mb_instance;              ///< mesh database
    moab::Interface &moab = mb_instance; ///< mesh database interface
    //! [Create MoAB]

    //! [Create MoFEM]
    MoFEM::Core core(moab);           ///< finite element database
    MoFEM::Interface &m_field = core; ///< finite element database interface
    //! [Create MoFEM]

    //! [Load mesh]
    Simple *simple = m_field.getInterface<Simple>();
    CHKERR simple->getOptions();
    CHKERR simple->loadFile();
    //! [Load mesh]

    //! [ThermoElasticProblem]
    ThermoElasticProblem ex(m_field);
    CHKERR ex.runProblem();
    //! [ThermoElasticProblem]
  }
  CATCH_ERRORS;

  CHKERR MoFEM::Core::Finalize();
}
