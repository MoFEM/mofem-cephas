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

#include <ThermalConvection.hpp>
#include <ThermalRadiation.hpp>

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
using SkinPostProcEle = PostProcBrokenMeshInMoab<BoundaryEle>;
using SideEle = PipelineManager::ElementsAndOpsByDim<SPACE_DIM>::FaceSideEle;

using AssemblyDomainEleOp =
    FormsIntegrators<DomainEleOp>::Assembly<PETSC>::OpBase;

constexpr AssemblyType AT = AssemblyType::PETSC; //< selected assembly type
constexpr IntegrationType IT =
    IntegrationType::GAUSS; //< selected integration type

//! [Linear elastic problem]
#include <HenckyOps.hpp> // Include Hencky operators
using namespace HenckyOps;
//! [Linear elastic problem]

// Include finite deformation operators
#include <FiniteThermalOps.hpp> // Include operators for finite strain diffusion problem
using namespace FiniteThermalOps;

//! [Thermal problem]
#include <ThermalOps.hpp> // Include operators for thermal problem which are agnostic to small/large strains
using namespace ThermalOps;
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

//! [Natural boundary conditions]

//! [Essential boundary conditions (Least square approach)]
using OpEssentialFluxRhs =
    EssentialBC<BoundaryEleOp>::Assembly<PETSC>::LinearForm<
        GAUSS>::OpEssentialRhs<HeatFluxCubitBcData, 3, SPACE_DIM>;
using OpEssentialFluxLhs =
    EssentialBC<BoundaryEleOp>::Assembly<PETSC>::BiLinearForm<
        GAUSS>::OpEssentialLhs<HeatFluxCubitBcData, 3, SPACE_DIM>;
//! [Essential boundary conditions (Least square approach)]

//! [Default input parameters]
double default_young_modulus = 1;
double default_poisson_ratio = 0.25;
double default_ref_temp = 0.0;
double default_init_temp = 0.0;

double default_coeff_expansion = 1;
double default_heat_conductivity =
    1; // Force / (time temperature )  or Power /
       // (length temperature) // Time unit is hour. force unit kN
double default_heat_capacity = 1; // length^2/(time^2 temperature) // length is
                                  // millimeter time is hour

int order_temp = 2; //< default approximation order for the temperature field
int order_flux = 3; //< default approximation order for heat flux field
int order_disp = 3; //< default approximation order for the displacement field

int atom_test = 0;

int save_every = 1; //< Save every n-th step
PetscBool do_output_domain;
PetscBool do_output_skin;
//! [Default input parameters]

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

  boost::shared_ptr<VectorDouble> tempFieldPtr;
  boost::shared_ptr<MatrixDouble> fluxFieldPtr;
  boost::shared_ptr<MatrixDouble> dispFieldPtr;
  boost::shared_ptr<MatrixDouble> dispGradPtr;
  boost::shared_ptr<MatrixDouble> strainFieldPtr;
  boost::shared_ptr<MatrixDouble> stressFieldPtr;

  MoFEMErrorCode setupProblem(); ///< add fields
  MoFEMErrorCode
  getCommandLineParameters(); //< read parameters from command line
  MoFEMErrorCode bC();        //< read boundary conditions
  MoFEMErrorCode OPs();       //< add operators to pipeline
  MoFEMErrorCode tsSolve();   //< time solver

  struct BlockedThermalParameters
      : public boost::enable_shared_from_this<BlockedThermalParameters> {
    double heatConductivity;
    double heatCapacity;

    inline auto getHeatConductivityPtr() {
      return boost::shared_ptr<double>(shared_from_this(), &heatConductivity);
    }

    inline auto getHeatCapacityPtr() {
      return boost::shared_ptr<double>(shared_from_this(), &heatCapacity);
    }
  };

  struct BlockedThermoElasticParameters
      : public boost::enable_shared_from_this<BlockedThermoElasticParameters> {
    VectorDouble coeffExpansion;
    double refTemp;

    inline auto getCoeffExpansionPtr() {
      return boost::shared_ptr<VectorDouble>(shared_from_this(),
                                             &coeffExpansion);
    }

    inline auto getRefTempPtr() {
      return boost::shared_ptr<double>(shared_from_this(), &refTemp);
    }
  };

  MoFEMErrorCode addMatThermalBlockOps(
      boost::ptr_deque<ForcesAndSourcesCore::UserDataOperator> &pipeline,
      std::string block_name,
      boost::shared_ptr<BlockedThermalParameters> blockedParamsPtr, Sev sev);

  MoFEMErrorCode addMatThermoElasticBlockOps(
      boost::ptr_deque<ForcesAndSourcesCore::UserDataOperator> &pipeline,
      std::string block_name,
      boost::shared_ptr<BlockedThermoElasticParameters> blockedParamsPtr,
      Sev sev);

  template <int DIM, AssemblyType A, IntegrationType I, typename DomainEleOp>
  MoFEMErrorCode opThermoElasticFactoryDomainRhs(
      MoFEM::Interface &m_field,
      boost::ptr_deque<ForcesAndSourcesCore::UserDataOperator> &pip,
      std::string field_name,
      boost::shared_ptr<HenckyOps::CommonData> elastic_common_ptr,
      boost::shared_ptr<ThermoElasticProblem::BlockedThermalParameters>
          thermal_common_ptr,
      boost::shared_ptr<ThermoElasticProblem::BlockedThermoElasticParameters>
          thermoelastic_common_ptr,
      Sev sev) {
    MoFEMFunctionBegin;

    using B = typename FormsIntegrators<DomainEleOp>::template Assembly<
        A>::template LinearForm<I>;
    using H = HenckyOps::HenckyIntegrators<DomainEleOp>;
    auto vec_temp_ptr = boost::make_shared<VectorDouble>();
    pip.push_back(new OpCalculateScalarFieldValues("T", vec_temp_ptr));
    auto coeff_expansion_ptr = thermoelastic_common_ptr->getCoeffExpansionPtr();
    auto ref_temp_ptr = thermoelastic_common_ptr->getRefTempPtr();
    pip.push_back(
        new typename H::template OpCalculateHenckyThermalStress<DIM, I, 0>(
            "U", vec_temp_ptr, elastic_common_ptr, coeff_expansion_ptr,
            ref_temp_ptr));
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
      std::string thermal_block_name, std::string thermoelastic_block_name,
      Sev sev, double scale = 1) {
    MoFEMFunctionBegin;

    auto elastic_common_ptr = commonDataFactory<DIM, I, DomainEleOp>(
        m_field, pip, field_name, elastic_block_name, sev, scale);
    auto thermal_common_ptr = boost::make_shared<BlockedThermalParameters>();
    CHKERR addMatThermalBlockOps(pip, thermal_block_name, thermal_common_ptr,
                                 Sev::inform);
    auto thermoelastic_common_ptr =
        boost::make_shared<BlockedThermoElasticParameters>();
    CHKERR addMatThermoElasticBlockOps(pip, thermoelastic_block_name,
                                       thermoelastic_common_ptr, Sev::inform);
    CHKERR opThermoElasticFactoryDomainRhs<DIM, A, I, DomainEleOp>(
        m_field, pip, field_name, elastic_common_ptr, thermal_common_ptr,
        thermoelastic_common_ptr, sev);

    MoFEMFunctionReturn(0);
  }

  template <int DIM, AssemblyType A, IntegrationType I, typename DomainEleOp>
  MoFEMErrorCode opThermoElasticFactoryDomainLhs(
      MoFEM::Interface &m_field,
      boost::ptr_deque<ForcesAndSourcesCore::UserDataOperator> &pip,
      std::string field_name, std::string coupled_field_name,
      boost::shared_ptr<HenckyOps::CommonData> elastic_common_ptr,
      boost::shared_ptr<ThermoElasticProblem::BlockedThermalParameters>
          thermal_common_ptr,
      boost::shared_ptr<ThermoElasticProblem::BlockedThermoElasticParameters>
          thermoelastic_common_ptr,
      Sev sev) {
    MoFEMFunctionBegin;

    using B = typename FormsIntegrators<DomainEleOp>::template Assembly<
        A>::template BiLinearForm<I>;
    using OpKPiola = typename B::template OpGradTensorGrad<1, DIM, DIM, 1>;

    using H = HenckyIntegrators<DomainEleOp>;
    auto vec_temp_ptr = boost::make_shared<VectorDouble>();
    pip.push_back(new OpCalculateScalarFieldValues("T", vec_temp_ptr));
    auto coeff_expansion_ptr = thermoelastic_common_ptr->getCoeffExpansionPtr();
    auto ref_temp_ptr = thermoelastic_common_ptr->getRefTempPtr();
    pip.push_back(
        new typename H::template OpCalculateHenckyThermalStress<DIM, I, 0>(
            "U", vec_temp_ptr, elastic_common_ptr, coeff_expansion_ptr,
            ref_temp_ptr));
    pip.push_back(new typename H::template OpHenckyTangent<DIM, I, 0>(
        field_name, elastic_common_ptr));
    pip.push_back(new OpKPiola(field_name, field_name,
                               elastic_common_ptr->getMatTangent()));
    pip.push_back(new typename H::template OpCalculateHenckyThermalStressdT<
                  DIM, I, AssemblyDomainEleOp, 0>(
        field_name, coupled_field_name, elastic_common_ptr,
        coeff_expansion_ptr));

    MoFEMFunctionReturn(0);
  }

  template <int DIM, AssemblyType A, IntegrationType I, typename DomainEleOp>
  MoFEMErrorCode opThermoElasticFactoryDomainLhs(
      MoFEM::Interface &m_field,
      boost::ptr_deque<ForcesAndSourcesCore::UserDataOperator> &pip,
      std::string field_name, std::string coupled_field_name,
      std::string elastic_block_name, std::string thermal_block_name,
      std::string thermoelastic_block_name, Sev sev, double scale = 1) {
    MoFEMFunctionBegin;

    auto elastic_common_ptr = commonDataFactory<DIM, I, DomainEleOp>(
        m_field, pip, field_name, elastic_block_name, sev, scale);
    auto thermal_common_ptr = boost::make_shared<BlockedThermalParameters>();
    CHKERR addMatThermalBlockOps(pip, thermal_block_name, thermal_common_ptr,
                                 Sev::inform);
    auto thermoelastic_common_ptr =
        boost::make_shared<BlockedThermoElasticParameters>();
    CHKERR addMatThermoElasticBlockOps(pip, thermoelastic_block_name,
                                       thermoelastic_common_ptr, Sev::inform);
    CHKERR opThermoElasticFactoryDomainLhs<DIM, A, I, DomainEleOp>(
        m_field, pip, field_name, coupled_field_name, elastic_common_ptr,
        thermal_common_ptr, thermoelastic_common_ptr, sev);

    MoFEMFunctionReturn(0);
  }
};

MoFEMErrorCode ThermoElasticProblem::addMatThermalBlockOps(
    boost::ptr_deque<ForcesAndSourcesCore::UserDataOperator> &pipeline,
    std::string block_name,
    boost::shared_ptr<BlockedThermalParameters> blockedParamsPtr, Sev sev) {
  MoFEMFunctionBegin;

  struct OpMatThermalBlocks : public DomainEleOp {
    OpMatThermalBlocks(boost::shared_ptr<double> conductivity_ptr,
                       boost::shared_ptr<double> capacity_ptr,
                       MoFEM::Interface &m_field, Sev sev,
                       std::vector<const CubitMeshSets *> meshset_vec_ptr)
        : DomainEleOp(NOSPACE, DomainEleOp::OPSPACE),
          conductivityPtr(conductivity_ptr), capacityPtr(capacity_ptr) {
      CHK_THROW_MESSAGE(extractThermalBlockData(m_field, meshset_vec_ptr, sev),
                        "Cannot get data from thermal block");
    }

    MoFEMErrorCode doWork(int side, EntityType type,
                          EntitiesFieldData::EntData &data) {
      MoFEMFunctionBegin;

      for (auto &b : blockData) {

        if (b.blockEnts.find(getFEEntityHandle()) != b.blockEnts.end()) {
          *conductivityPtr = b.conductivity;
          *capacityPtr = b.capacity;
          MoFEMFunctionReturnHot(0);
        }
      }

      *conductivityPtr = default_heat_conductivity;
      *capacityPtr = default_heat_capacity;

      MoFEMFunctionReturn(0);
    }

  private:
    struct BlockData {
      double conductivity;
      double capacity;
      Range blockEnts;
    };

    std::vector<BlockData> blockData;

    MoFEMErrorCode
    extractThermalBlockData(MoFEM::Interface &m_field,
                            std::vector<const CubitMeshSets *> meshset_vec_ptr,
                            Sev sev) {
      MoFEMFunctionBegin;

      for (auto m : meshset_vec_ptr) {
        MOFEM_TAG_AND_LOG("WORLD", sev, "Mat Thermal Block") << *m;
        std::vector<double> block_data;
        CHKERR m->getAttributes(block_data);
        if (block_data.size() < 2) {
          SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
                  "Expected that block has at least two attributes");
        }
        auto get_block_ents = [&]() {
          Range ents;
          CHKERR
          m_field.get_moab().get_entities_by_handle(m->meshset, ents, true);
          return ents;
        };

        MOFEM_TAG_AND_LOG("WORLD", sev, "Mat Thermal Block")
            << m->getName() << ": conductivity = " << block_data[0]
            << " capacity  = " << block_data[1];

        blockData.push_back({block_data[0], block_data[1], get_block_ents()});
      }
      MOFEM_LOG_CHANNEL("WORLD");
      MoFEMFunctionReturn(0);
    }

    boost::shared_ptr<double> conductivityPtr;
    boost::shared_ptr<double> capacityPtr;
  };

  pipeline.push_back(new OpMatThermalBlocks(
      blockedParamsPtr->getHeatConductivityPtr(),
      blockedParamsPtr->getHeatCapacityPtr(), mField, sev,

      // Get blockset using regular expression
      mField.getInterface<MeshsetsManager>()->getCubitMeshsetPtr(std::regex(

          (boost::format("%s(.*)") % block_name).str()

              ))

          ));

  MoFEMFunctionReturn(0);
}

MoFEMErrorCode ThermoElasticProblem::addMatThermoElasticBlockOps(
    boost::ptr_deque<ForcesAndSourcesCore::UserDataOperator> &pipeline,
    std::string block_name,
    boost::shared_ptr<BlockedThermoElasticParameters> blockedParamsPtr,
    Sev sev) {
  MoFEMFunctionBegin;

  struct OpMatThermoElasticBlocks : public DomainEleOp {
    OpMatThermoElasticBlocks(boost::shared_ptr<VectorDouble> expansion_ptr,
                             boost::shared_ptr<double> ref_temp_ptr,
                             MoFEM::Interface &m_field, Sev sev,
                             std::vector<const CubitMeshSets *> meshset_vec_ptr)
        : DomainEleOp(NOSPACE, DomainEleOp::OPSPACE),
          expansionPtr(expansion_ptr), refTempPtr(ref_temp_ptr) {
      CHK_THROW_MESSAGE(
          extractThermoElasticBlockData(m_field, meshset_vec_ptr, sev),
          "Cannot get data from thermoelastic block");
    }

    MoFEMErrorCode doWork(int side, EntityType type,
                          EntitiesFieldData::EntData &data) {
      MoFEMFunctionBegin;

      for (auto &b : blockData) {

        if (b.blockEnts.find(getFEEntityHandle()) != b.blockEnts.end()) {
          *expansionPtr = b.expansion;
          *refTempPtr = b.ref_temp;
          MoFEMFunctionReturnHot(0);
        }
      }

      *expansionPtr = VectorDouble(SPACE_DIM, default_coeff_expansion);
      *refTempPtr = default_ref_temp;

      MoFEMFunctionReturn(0);
    }

  private:
    struct BlockData {
      double ref_temp;
      VectorDouble expansion;
      Range blockEnts;
    };

    std::vector<BlockData> blockData;

    MoFEMErrorCode extractThermoElasticBlockData(
        MoFEM::Interface &m_field,
        std::vector<const CubitMeshSets *> meshset_vec_ptr, Sev sev) {
      MoFEMFunctionBegin;

      for (auto m : meshset_vec_ptr) {
        MOFEM_TAG_AND_LOG("WORLD", sev, "Mat Thermoelastic Block") << *m;
        std::vector<double> block_data;
        CHKERR m->getAttributes(block_data);
        if (block_data.size() < 2) {
          SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
                  "Expected that block has at least two attributes");
        }
        auto get_block_ents = [&]() {
          Range ents;
          CHKERR
          m_field.get_moab().get_entities_by_handle(m->meshset, ents, true);
          return ents;
        };

        auto get_expansion = [&]() {
          VectorDouble expansion(SPACE_DIM, block_data[1]);
          if (block_data.size() > 2) {
            expansion[1] = block_data[2];
          }
          if (SPACE_DIM == 3 && block_data.size() > 3) {
            expansion[2] = block_data[3];
          }
          return expansion;
        };

        auto coeff_exp_vec = get_expansion();

        MOFEM_TAG_AND_LOG("WORLD", sev, "Mat Thermoelastic Block")
            << " ref_temp = " << block_data[0]
            << " expansion = " << coeff_exp_vec;

        blockData.push_back({block_data[0], coeff_exp_vec, get_block_ents()});
      }
      MOFEM_LOG_CHANNEL("WORLD");
      MoFEMFunctionReturn(0);
    }

    boost::shared_ptr<VectorDouble> expansionPtr;
    boost::shared_ptr<double> refTempPtr;
  };

  pipeline.push_back(new OpMatThermoElasticBlocks(
      blockedParamsPtr->getCoeffExpansionPtr(),
      blockedParamsPtr->getRefTempPtr(), mField, sev,

      // Get blockset using regular expression
      mField.getInterface<MeshsetsManager>()->getCubitMeshsetPtr(std::regex(

          (boost::format("%s(.*)") % block_name).str()

              ))

          ));

  MoFEMFunctionReturn(0);
}

//! [Run problem]
MoFEMErrorCode ThermoElasticProblem::runProblem() {
  MoFEMFunctionBegin;
  CHKERR getCommandLineParameters();
  CHKERR setupProblem();
  CHKERR bC();
  CHKERR OPs();
  CHKERR tsSolve();
  MoFEMFunctionReturn(0);
}
//! [Run problem]

//! [Get command line parameters]
MoFEMErrorCode ThermoElasticProblem::getCommandLineParameters() {
  MoFEMFunctionBegin;

  auto get_command_line_parameters = [&]() {
    MoFEMFunctionBegin;

    CHKERR PetscOptionsGetInt(PETSC_NULL, "", "-order", &order_temp,
                              PETSC_NULL);
    CHKERR PetscOptionsGetInt(PETSC_NULL, "", "-order_temp", &order_temp,
                              PETSC_NULL);
    order_flux = order_temp + 1;
    CHKERR PetscOptionsGetInt(PETSC_NULL, "", "-order_flux", &order_flux,
                              PETSC_NULL);
    order_disp = order_temp + 1;
    CHKERR PetscOptionsGetInt(PETSC_NULL, "", "-order_disp", &order_disp,
                              PETSC_NULL);
    CHKERR PetscOptionsGetInt(PETSC_NULL, "", "-atom_test", &atom_test,
                              PETSC_NULL);
    CHKERR PetscOptionsGetInt(PETSC_NULL, "", "-save_every", &save_every,
                              PETSC_NULL);

    CHKERR PetscOptionsGetScalar(PETSC_NULL, "", "-young_modulus",
                                 &default_young_modulus, PETSC_NULL);
    CHKERR PetscOptionsGetScalar(PETSC_NULL, "", "-poisson_ratio",
                                 &default_poisson_ratio, PETSC_NULL);
    CHKERR PetscOptionsGetScalar(PETSC_NULL, "", "-coeff_expansion",
                                 &default_coeff_expansion, PETSC_NULL);
    CHKERR PetscOptionsGetScalar(PETSC_NULL, "", "-ref_temp", &default_ref_temp,
                                 PETSC_NULL);
    CHKERR PetscOptionsGetScalar(PETSC_NULL, "", "-init_temp",
                                 &default_init_temp, PETSC_NULL);

    CHKERR PetscOptionsGetScalar(PETSC_NULL, "", "-capacity",
                                 &default_heat_capacity, PETSC_NULL);
    CHKERR PetscOptionsGetScalar(PETSC_NULL, "", "-conductivity",
                                 &default_heat_conductivity, PETSC_NULL);

    if constexpr (SPACE_DIM == 2) {
      do_output_domain = PETSC_TRUE;
      do_output_skin = PETSC_FALSE;
    } else {
      do_output_domain = PETSC_FALSE;
      do_output_skin = PETSC_TRUE;
    }

    CHKERR PetscOptionsGetBool(PETSC_NULL, "", "-output_domain",
                               &do_output_domain, PETSC_NULL);
    CHKERR PetscOptionsGetBool(PETSC_NULL, "", "-output_skin", &do_output_skin,
                               PETSC_NULL);

    MOFEM_LOG("ThermoElastic", Sev::inform)
        << "Young's modulus " << default_young_modulus;
    MOFEM_LOG("ThermoElastic", Sev::inform)
        << "Poisson's ratio " << default_poisson_ratio;
    MOFEM_LOG("ThermoElastic", Sev::inform)
        << "Coeff of expansion " << default_coeff_expansion;
    MOFEM_LOG("ThermoElastic", Sev::inform)
        << "Default heat capacity " << default_heat_capacity;
    MOFEM_LOG("ThermoElastic", Sev::inform)
        << "Heat conductivity " << default_heat_conductivity;

    MOFEM_LOG("ThermoElastic", Sev::inform)
        << "Reference temperature  " << default_ref_temp;
    MOFEM_LOG("ThermoElastic", Sev::inform)
        << "Initial temperature  " << default_init_temp;

    MoFEMFunctionReturn(0);
  };

  CHKERR get_command_line_parameters();

  MoFEMFunctionReturn(0);
}
//! [Get command line parameters]

//! [Set up problem]
MoFEMErrorCode ThermoElasticProblem::setupProblem() {
  MoFEMFunctionBegin;
  Simple *simple = mField.getInterface<Simple>();
  // Add field
  constexpr FieldApproximationBase base = DEMKOWICZ_JACOBI_BASE;
  // Mechanical fields
  CHKERR simple->addDomainField("U", H1, AINSWORTH_LEGENDRE_BASE, SPACE_DIM);
  CHKERR simple->addBoundaryField("U", H1, AINSWORTH_LEGENDRE_BASE, SPACE_DIM);
  // Temperature
  constexpr auto flux_space = (SPACE_DIM == 2) ? HCURL : HDIV;
  CHKERR simple->addDomainField("T", L2, base, 1);
  CHKERR simple->addDomainField("FLUX", flux_space, base, 1);
  CHKERR simple->addBoundaryField("FLUX", flux_space, base, 1);
  CHKERR simple->addBoundaryField("TBC", L2, base, 1);

  CHKERR simple->setFieldOrder("U", order_disp);
  CHKERR simple->setFieldOrder("FLUX", order_flux);
  CHKERR simple->setFieldOrder("T", order_temp);
  CHKERR simple->setFieldOrder("TBC", order_temp);

  CHKERR simple->setUp();

  int coords_dim = SPACE_DIM;
  CHKERR PetscOptionsGetRealArray(NULL, NULL, "-field_eval_coords",
                                  fieldEvalCoords.data(), &coords_dim,
                                  &doEvalField);

  tempFieldPtr = boost::make_shared<VectorDouble>();
  fluxFieldPtr = boost::make_shared<MatrixDouble>();
  dispFieldPtr = boost::make_shared<MatrixDouble>();
  dispGradPtr = boost::make_shared<MatrixDouble>();
  strainFieldPtr = boost::make_shared<MatrixDouble>();
  stressFieldPtr = boost::make_shared<MatrixDouble>();

  if (doEvalField) {
    fieldEvalData =
        mField.getInterface<FieldEvaluatorInterface>()->getData<DomainEle>();

    CHKERR mField.getInterface<FieldEvaluatorInterface>()->buildTree<SPACE_DIM>(
        fieldEvalData, simple->getDomainFEName());

    fieldEvalData->setEvalPoints(fieldEvalCoords.data(), 1);
    auto no_rule = [](int, int, int) { return -1; };

    auto field_eval_fe_ptr = fieldEvalData->feMethodPtr.lock();
    field_eval_fe_ptr->getRuleHook = no_rule;

    auto block_thermal_params = boost::make_shared<BlockedThermalParameters>();
    auto block_thermoelastic_params =
        boost::make_shared<BlockedThermoElasticParameters>();
    auto coeff_expansion_ptr =
        block_thermoelastic_params->getCoeffExpansionPtr();
    auto ref_temp_ptr = block_thermoelastic_params->getRefTempPtr();

    CHKERR addMatThermalBlockOps(field_eval_fe_ptr->getOpPtrVector(),
                                 "MAT_THERMAL", block_thermal_params,
                                 Sev::verbose);
    CHKERR addMatThermoElasticBlockOps(
        field_eval_fe_ptr->getOpPtrVector(), "MAT_THERMOELASTIC",
        block_thermoelastic_params, Sev::verbose);

    CHKERR AddHOOps<SPACE_DIM, SPACE_DIM, SPACE_DIM>::add(
        field_eval_fe_ptr->getOpPtrVector(), {H1, HDIV});

    auto hencky_common_data_ptr =
        HenckyOps::commonDataFactory<SPACE_DIM, IT, DomainEleOp>(
            mField, field_eval_fe_ptr->getOpPtrVector(), "U", "MAT_ELASTIC",
            Sev::inform);
    auto mat_D_ptr = hencky_common_data_ptr->matDPtr;
    auto dispGradPtr = hencky_common_data_ptr->matGradPtr;
    auto mat_stress_ptr = boost::make_shared<MatrixDouble>();

    field_eval_fe_ptr->getOpPtrVector().push_back(
        new OpCalculateScalarFieldValues("T", tempFieldPtr));
    field_eval_fe_ptr->getOpPtrVector().push_back(
        new OpCalculateHVecVectorField<3, SPACE_DIM>("FLUX", fluxFieldPtr));
    field_eval_fe_ptr->getOpPtrVector().push_back(
        new OpCalculateVectorFieldValues<SPACE_DIM>("U", dispFieldPtr));
    field_eval_fe_ptr->getOpPtrVector().push_back(
        new OpCalculateVectorFieldGradient<SPACE_DIM, SPACE_DIM>("U",
                                                                 dispGradPtr));

    using H = HenckyOps::HenckyIntegrators<DomainEleOp>;

    field_eval_fe_ptr->getOpPtrVector().push_back(
        new
        typename H::template OpCalculateHenckyThermalStress<SPACE_DIM, IT, 0>(
            "U", tempFieldPtr, hencky_common_data_ptr, coeff_expansion_ptr,
            ref_temp_ptr));
    if (!IS_LARGE_STRAINS) {
      field_eval_fe_ptr->getOpPtrVector().push_back(
          new OpSymmetrizeTensor<SPACE_DIM>(
              hencky_common_data_ptr->getMatFirstPiolaStress(),
              stressFieldPtr));
      field_eval_fe_ptr->getOpPtrVector().push_back(
          new OpSymmetrizeTensor<SPACE_DIM>(dispGradPtr, strainFieldPtr));
    } else {
      field_eval_fe_ptr->getOpPtrVector().push_back(
          new typename H::template OpCalculateLogC<SPACE_DIM, IT>(
              "U", hencky_common_data_ptr));
      stressFieldPtr = hencky_common_data_ptr->getMatFirstPiolaStress();
      strainFieldPtr = hencky_common_data_ptr->getMatLogC();
    };
  }

  MoFEMFunctionReturn(0);
}
//! [Set up problem]

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

  auto filter_flux_blocks = [&](auto skin, bool temp_bc = false) {
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
    if (!temp_bc) {
      CHK_THROW_MESSAGE(remove_cubit_blocks(NODESET | TEMPERATURESET),
                        "remove_cubit_blocks");
      CHK_THROW_MESSAGE(remove_named_blocks("TEMPERATURE"),
                        "remove_named_blocks");
    }
    CHK_THROW_MESSAGE(remove_cubit_blocks(SIDESET | HEATFLUXSET),
                      "remove_cubit_blocks");
    CHK_THROW_MESSAGE(remove_named_blocks("HEATFLUX"), "remove_named_blocks");
    CHK_THROW_MESSAGE(remove_named_blocks("CONVECTION"), "remove_named_blocks");
    CHK_THROW_MESSAGE(remove_named_blocks("RADIATION"), "remove_named_blocks");
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
  auto remove_temp_bc_ents =
      filter_true_skin(filter_flux_blocks(get_skin(), true));

  CHKERR mField.getInterface<CommInterface>()->synchroniseEntities(
      remove_flux_ents);
  CHKERR mField.getInterface<CommInterface>()->synchroniseEntities(
      remove_temp_bc_ents);

  MOFEM_LOG("SYNC", Sev::noisy) << remove_flux_ents << endl;
  MOFEM_LOG_SEVERITY_SYNC(mField.get_comm(), Sev::noisy);

  MOFEM_LOG("SYNC", Sev::noisy) << remove_temp_bc_ents << endl;
  MOFEM_LOG_SEVERITY_SYNC(mField.get_comm(), Sev::noisy);

#ifndef NDEBUG

  CHKERR save_range(
      mField.get_moab(),
      (boost::format("flux_remove_%d.vtk") % mField.get_comm_rank()).str(),
      remove_flux_ents);

  CHKERR save_range(
      mField.get_moab(),
      (boost::format("temp_bc_remove_%d.vtk") % mField.get_comm_rank()).str(),
      remove_temp_bc_ents);

#endif

  CHKERR mField.getInterface<ProblemsManager>()->removeDofsOnEntities(
      simple->getProblemName(), "FLUX", remove_flux_ents);
  CHKERR mField.getInterface<ProblemsManager>()->removeDofsOnEntities(
      simple->getProblemName(), "TBC", remove_temp_bc_ents);

  auto set_init_temp = [](boost::shared_ptr<FieldEntity> field_entity_ptr) {
    field_entity_ptr->getEntFieldData()[0] = default_init_temp;
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

  auto block_thermal_params = boost::make_shared<BlockedThermalParameters>();
  auto heat_conductivity_ptr = block_thermal_params->getHeatConductivityPtr();
  auto heat_capacity_ptr = block_thermal_params->getHeatCapacityPtr();

  auto block_thermoelastic_params =
      boost::make_shared<BlockedThermoElasticParameters>();
  auto coeff_expansion_ptr = block_thermoelastic_params->getCoeffExpansionPtr();
  auto ref_temp_ptr = block_thermoelastic_params->getRefTempPtr();

  // Default time scaling of BCs and sources, from command line arguments
  auto time_scale =
      boost::make_shared<TimeScale>("", false, [](const double) { return 1; });
  auto def_time_scale = [time_scale](const double time) {
    return (!time_scale->argFileScale) ? time : 1;
  };
  auto def_file_name = [time_scale](const std::string &&name) {
    return (!time_scale->argFileScale) ? name : "";
  };

  // Files which define scaling for separate variables, if provided
  auto time_bodyforce_scaling = boost::make_shared<TimeScale>(
      def_file_name("bodyforce_scale.txt"), false, def_time_scale);
  auto time_heatsource_scaling = boost::make_shared<TimeScale>(
      def_file_name("heatsource_scale.txt"), false, def_time_scale);
  auto time_temperature_scaling = boost::make_shared<TimeScale>(
      def_file_name("temperature_bc_scale.txt"), false, def_time_scale);
  auto time_displacement_scaling = boost::make_shared<TimeScale>(
      def_file_name("displacement_bc_scale.txt"), false, def_time_scale);
  auto time_flux_scaling = boost::make_shared<TimeScale>(
      def_file_name("flux_bc_scale.txt"), false, def_time_scale);
  auto time_force_scaling = boost::make_shared<TimeScale>(
      def_file_name("force_bc_scale.txt"), false, def_time_scale);

  auto add_domain_rhs_ops = [&](auto &pipeline) {
    MoFEMFunctionBegin;
    CHKERR addMatThermalBlockOps(pipeline, "MAT_THERMAL", block_thermal_params,
                                 Sev::inform);
    CHKERR addMatThermoElasticBlockOps(pipeline, "MAT_THERMOELASTIC",
                                       block_thermoelastic_params, Sev::inform);
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
        mField, pipeline, "U", "MAT_ELASTIC", "MAT_THERMAL",
        "MAT_THERMOELASTIC", Sev::inform);

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
    CHKERR addMatThermalBlockOps(pipeline, "MAT_THERMAL", block_thermal_params,
                                 Sev::verbose);
    CHKERR addMatThermoElasticBlockOps(pipeline, "MAT_THERMOELASTIC",
                                       block_thermoelastic_params,
                                       Sev::verbose);
    CHKERR AddHOOps<SPACE_DIM, SPACE_DIM, SPACE_DIM>::add(pipeline, {H1, HDIV});

    auto hencky_common_data_ptr =
        HenckyOps::commonDataFactory<SPACE_DIM, IT, DomainEleOp>(
            mField, pipeline, "U", "MAT_ELASTIC", Sev::inform, 1);
    auto mat_D_ptr = hencky_common_data_ptr->matDPtr;
    auto mat_grad_ptr = hencky_common_data_ptr->matGradPtr;

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
        mField, pipeline, "U", "T", "MAT_ELASTIC", "MAT_THERMAL",
        "MAT_THERMOELASTIC", Sev::inform);

    auto op_capacity = new OpCapacity("T", "T", capacity);
    op_capacity->feScalingFun = [](const FEMethod *fe_ptr) {
      return fe_ptr->ts_a;
    };
    pipeline.push_back(op_capacity);

    auto vec_temp_ptr = boost::make_shared<VectorDouble>();
    pipeline.push_back(new OpCalculateScalarFieldValues("T", vec_temp_ptr));
    CHKERR DomainNaturalBCLhs::AddFluxToPipeline<OpSetTemperatureLhs>::add(
        pipeline, mField, "T", vec_temp_ptr, "SETTEMP", Sev::verbose);

    MoFEMFunctionReturn(0);
  };

  auto add_boundary_rhs_ops = [&](auto &pipeline) {
    MoFEMFunctionBegin;

    CHKERR AddHOOps<SPACE_DIM - 1, SPACE_DIM, SPACE_DIM>::add(pipeline, {HDIV});

    CHKERR BoundaryNaturalBC::AddFluxToPipeline<OpForce>::add(
        pipeline, mField, "U", {time_scale, time_force_scaling}, "FORCE",
        Sev::inform);

    // Temperature BC

    using OpTemperatureBC =
        BoundaryNaturalBC::OpFlux<NaturalTemperatureMeshsets, 3, SPACE_DIM>;
    CHKERR BoundaryNaturalBC::AddFluxToPipeline<OpTemperatureBC>::add(
        pipeline, mField, "FLUX", {time_scale, time_temperature_scaling},
        "TEMPERATURE", Sev::inform);

    // Note: fluxes for temperature should be aggregated. Here separate is
    // NaturalMeshsetType<HEATFLUXSET>, we should add version with BLOCKSET,
    // convection, see example tutorials/vec-0/src/NaturalBoundaryBC.hpp.
    // and vec-0/elastic.cpp

    using OpFluxBC =
        BoundaryNaturalBC::OpFlux<NaturalMeshsetType<HEATFLUXSET>, 1, 1>;
    CHKERR BoundaryNaturalBC::AddFluxToPipeline<OpFluxBC>::add(
        pipeline, mField, "TBC", {time_scale, time_flux_scaling}, "FLUX",
        Sev::inform);

    using T = NaturalBC<BoundaryEleOp>::Assembly<PETSC>::LinearForm<GAUSS>;
    using OpConvectionRhsBC =
        T::OpFlux<ThermoElasticOps::ConvectionBcType<BLOCKSET>, 1, 1>;
    using OpRadiationRhsBC =
        T::OpFlux<ThermoElasticOps::RadiationBcType<BLOCKSET>, 1, 1>;
    auto temp_bc_ptr = boost::make_shared<VectorDouble>();
    pipeline.push_back(new OpCalculateScalarFieldValues("TBC", temp_bc_ptr));
    T::AddFluxToPipeline<OpConvectionRhsBC>::add(
        pipeline, mField, "TBC", temp_bc_ptr, "CONVECTION", Sev::inform);
    T::AddFluxToPipeline<OpRadiationRhsBC>::add(
        pipeline, mField, "TBC", temp_bc_ptr, "RADIATION", Sev::inform);

    auto mat_flux_ptr = boost::make_shared<MatrixDouble>();
    pipeline.push_back(
        new OpCalculateHVecVectorField<3, SPACE_DIM>("FLUX", mat_flux_ptr));

    // This is temporary implementation. It be moved to LinearFormsIntegrators,
    // however for hybridised case, what is goal of this changes, such function
    // is implemented for fluxes in broken space. Thus ultimately this operator
    // would be not needed.

    struct OpTBCTimesNormalFlux
        : public FormsIntegrators<BoundaryEleOp>::Assembly<PETSC>::OpBase {

      using OP = FormsIntegrators<BoundaryEleOp>::Assembly<PETSC>::OpBase;

      OpTBCTimesNormalFlux(const std::string field_name,
                           boost::shared_ptr<MatrixDouble> vec,
                           boost::shared_ptr<Range> ents_ptr = nullptr)
          : OP(field_name, field_name, OP::OPROW, ents_ptr), sourceVec(vec) {}

      MoFEMErrorCode iNtegrate(EntitiesFieldData::EntData &row_data) {
        MoFEMFunctionBegin;
        FTENSOR_INDEX(SPACE_DIM, i);
        // get integration weights
        auto t_w = OP::getFTensor0IntegrationWeight();
        // get base function values on rows
        auto t_row_base = row_data.getFTensor0N();
        // get normal vector
        auto t_normal = OP::getFTensor1NormalsAtGaussPts();
        // get vector values
        auto t_vec = getFTensor1FromMat<SPACE_DIM, 1>(*sourceVec);
        // loop over integration points
        for (int gg = 0; gg != OP::nbIntegrationPts; gg++) {
          // take into account Jacobian
          const double alpha = t_w * (t_vec(i) * t_normal(i));
          // loop over rows base functions
          int rr = 0;
          for (; rr != OP::nbRows; ++rr) {
            OP::locF[rr] += alpha * t_row_base;
            ++t_row_base;
          }
          for (; rr < OP::nbRowBaseFunctions; ++rr)
            ++t_row_base;
          ++t_w; // move to another integration weight
          ++t_vec;
          ++t_normal;
        }
        EntityType fe_type = OP::getNumeredEntFiniteElementPtr()->getEntType();
        if (fe_type == MBTRI) {
          OP::locF /= 2;
        }
        MoFEMFunctionReturn(0);
      }

    protected:
      boost::shared_ptr<MatrixDouble> sourceVec;
    };
    pipeline.push_back(new OpTBCTimesNormalFlux("TBC", mat_flux_ptr));

    struct OpBaseNormalTimesTBC
        : public FormsIntegrators<BoundaryEleOp>::Assembly<PETSC>::OpBase {

      using OP = FormsIntegrators<BoundaryEleOp>::Assembly<PETSC>::OpBase;

      OpBaseNormalTimesTBC(const std::string field_name,
                           boost::shared_ptr<VectorDouble> vec,
                           boost::shared_ptr<Range> ents_ptr = nullptr)
          : OP(field_name, field_name, OP::OPROW, ents_ptr), sourceVec(vec) {}

      MoFEMErrorCode iNtegrate(EntitiesFieldData::EntData &row_data) {
        MoFEMFunctionBegin;
        FTENSOR_INDEX(SPACE_DIM, i);
        // get integration weights
        auto t_w = OP::getFTensor0IntegrationWeight();
        // get base function values on rows
        auto t_row_base = row_data.getFTensor1N<3>();
        // get normal vector
        auto t_normal = OP::getFTensor1NormalsAtGaussPts();
        // get vector values
        auto t_vec = getFTensor0FromVec(*sourceVec);
        // loop over integration points
        for (int gg = 0; gg != OP::nbIntegrationPts; gg++) {
          // take into account Jacobian
          const double alpha = t_w * t_vec;
          // loop over rows base functions
          int rr = 0;
          for (; rr != OP::nbRows; ++rr) {
            OP::locF[rr] += alpha * (t_row_base(i) * t_normal(i));
            ++t_row_base;
          }
          for (; rr < OP::nbRowBaseFunctions; ++rr)
            ++t_row_base;
          ++t_w; // move to another integration weight
          ++t_vec;
          ++t_normal;
        }
        EntityType fe_type = OP::getNumeredEntFiniteElementPtr()->getEntType();
        if (fe_type == MBTRI) {
          OP::locF /= 2;
        }
        MoFEMFunctionReturn(0);
      }

    protected:
      boost::shared_ptr<VectorDouble> sourceVec;
    };

    pipeline.push_back(new OpBaseNormalTimesTBC("FLUX", temp_bc_ptr));

    MoFEMFunctionReturn(0);
  };

  auto add_boundary_lhs_ops = [&](auto &pipeline) {
    MoFEMFunctionBegin;

    CHKERR AddHOOps<SPACE_DIM - 1, SPACE_DIM, SPACE_DIM>::add(pipeline, {HDIV});

    using T =
        NaturalBC<BoundaryEleOp>::template Assembly<PETSC>::BiLinearForm<GAUSS>;
    using OpConvectionLhsBC =
        T::OpFlux<ThermoElasticOps::ConvectionBcType<BLOCKSET>, 1, 1>;
    using OpRadiationLhsBC =
        T::OpFlux<ThermoElasticOps::RadiationBcType<BLOCKSET>, 1, 1>;
    auto temp_bc_ptr = boost::make_shared<VectorDouble>();
    pipeline.push_back(new OpCalculateScalarFieldValues("TBC", temp_bc_ptr));
    T::AddFluxToPipeline<OpConvectionLhsBC>::add(pipeline, mField, "TBC", "TBC",
                                                 "CONVECTION", Sev::verbose);
    T::AddFluxToPipeline<OpRadiationLhsBC>::add(
        pipeline, mField, "TBC", "TBC", temp_bc_ptr, "RADIATION", Sev::verbose);

    // This is temporary implementation. It be moved to LinearFormsIntegrators,
    // however for hybridised case, what is goal of this changes, such function
    // is implemented for fluxes in broken space. Thus ultimately this operator
    // would be not needed.

    struct OpTBCTimesNormalFlux
        : public FormsIntegrators<BoundaryEleOp>::Assembly<PETSC>::OpBase {

      using OP = FormsIntegrators<BoundaryEleOp>::Assembly<PETSC>::OpBase;

      OpTBCTimesNormalFlux(const std::string row_field_name,
                           const std::string col_field_name,
                           boost::shared_ptr<Range> ents_ptr = nullptr)
          : OP(row_field_name, col_field_name, OP::OPROWCOL, ents_ptr) {
        this->sYmm = false;
        this->assembleTranspose = true;
        this->onlyTranspose = false;
      }

      MoFEMErrorCode iNtegrate(EntitiesFieldData::EntData &row_data,
                               EntitiesFieldData::EntData &col_data) {
        MoFEMFunctionBegin;

        FTENSOR_INDEX(SPACE_DIM, i);

        // get integration weights
        auto t_w = OP::getFTensor0IntegrationWeight();
        // get base function values on rows
        auto t_row_base = row_data.getFTensor0N();
        // get normal vector
        auto t_normal = OP::getFTensor1NormalsAtGaussPts();
        // loop over integration points
        for (int gg = 0; gg != OP::nbIntegrationPts; gg++) {
          // loop over rows base functions
          auto a_mat_ptr = &*OP::locMat.data().begin();
          int rr = 0;
          for (; rr != OP::nbRows; rr++) {
            // take into account Jacobian
            const double alpha = t_w * t_row_base;
            // get column base functions values at gauss point gg
            auto t_col_base = col_data.getFTensor1N<3>(gg, 0);
            // loop over columns
            for (int cc = 0; cc != OP::nbCols; cc++) {
              // calculate element of local matrix
              // using L2 norm of normal vector for convenient area calculation
              // for quads, tris etc.
              *a_mat_ptr += alpha * (t_col_base(i) * t_normal(i));
              ++t_col_base;
              ++a_mat_ptr;
            }
            ++t_row_base;
          }
          for (; rr < OP::nbRowBaseFunctions; ++rr)
            ++t_row_base;
          ++t_normal;
          ++t_w; // move to another integration weight
        }
        EntityType fe_type = OP::getNumeredEntFiniteElementPtr()->getEntType();
        if (fe_type == MBTRI) {
          OP::locMat /= 2;
        }
        MoFEMFunctionReturn(0);
      }
    };

    pipeline.push_back(new OpTBCTimesNormalFlux("TBC", "FLUX"));

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

  auto create_post_process_elements = [&]() {
    auto block_thermal_params = boost::make_shared<BlockedThermalParameters>();

    auto block_thermoelastic_params =
        boost::make_shared<BlockedThermoElasticParameters>();
    auto coeff_expansion_ptr =
        block_thermoelastic_params->getCoeffExpansionPtr();
    auto ref_temp_ptr = block_thermoelastic_params->getRefTempPtr();

    auto u_ptr = boost::make_shared<MatrixDouble>();
    auto mat_grad_ptr = boost::make_shared<MatrixDouble>();
    auto mat_strain_ptr = boost::make_shared<MatrixDouble>();
    auto mat_stress_ptr = boost::make_shared<MatrixDouble>();
    auto vec_temp_ptr = boost::make_shared<VectorDouble>();
    auto mat_flux_ptr = boost::make_shared<MatrixDouble>();

    auto push_domain_ops = [&](auto &pp_fe) {
      MoFEMFunctionBegin;
      auto &pip = pp_fe->getOpPtrVector();

      CHKERR addMatThermalBlockOps(pip, "MAT_THERMAL", block_thermal_params,
                                   Sev::verbose);
      CHKERR addMatThermoElasticBlockOps(
          pip, "MAT_THERMOELASTIC", block_thermoelastic_params, Sev::verbose);

      CHKERR AddHOOps<SPACE_DIM, SPACE_DIM, SPACE_DIM>::add(pip, {H1, HDIV});

      pip.push_back(new OpCalculateScalarFieldValues("T", vec_temp_ptr));
      pip.push_back(
          new OpCalculateHVecVectorField<3, SPACE_DIM>("FLUX", mat_flux_ptr));

      pip.push_back(new OpCalculateVectorFieldValues<SPACE_DIM>("U", u_ptr));
      pip.push_back(new OpCalculateVectorFieldGradient<SPACE_DIM, SPACE_DIM>(
          "U", mat_grad_ptr));
      auto elastic_common_ptr = commonDataFactory<SPACE_DIM, IT, DomainEleOp>(
          mField, pip, "U", "MAT_ELASTIC", Sev::inform);
      using H = HenckyOps::HenckyIntegrators<DomainEleOp>;
      pip.push_back(
          new
          typename H::template OpCalculateHenckyThermalStress<SPACE_DIM, IT, 0>(
              "U", vec_temp_ptr, elastic_common_ptr, coeff_expansion_ptr,
              ref_temp_ptr));
      if (!IS_LARGE_STRAINS) {
        pip.push_back(new OpSymmetrizeTensor<SPACE_DIM>(
            elastic_common_ptr->getMatFirstPiolaStress(), mat_stress_ptr));
        pip.push_back(
            new OpSymmetrizeTensor<SPACE_DIM>(mat_grad_ptr, mat_strain_ptr));
      } else {
        mat_stress_ptr = elastic_common_ptr->getMatFirstPiolaStress();
        mat_strain_ptr = elastic_common_ptr->getMatLogC();
      }

      MoFEMFunctionReturn(0);
    };

    auto push_post_proc_ops = [&](auto &pp_fe) {
      MoFEMFunctionBegin;
      auto &pip = pp_fe->getOpPtrVector();
      using OpPPMap = OpPostProcMapInMoab<SPACE_DIM, SPACE_DIM>;

      if (!IS_LARGE_STRAINS) {
        pip.push_back(

            new OpPPMap(

                pp_fe->getPostProcMesh(), pp_fe->getMapGaussPts(),

                {{"T", vec_temp_ptr}},

                {{"U", u_ptr}, {"FLUX", mat_flux_ptr}},

                {},

                {{"CAUCHY", mat_stress_ptr}, {"STRAIN", mat_strain_ptr}}

                )

        );
      } else {
        pip.push_back(

            new OpPPMap(

                pp_fe->getPostProcMesh(), pp_fe->getMapGaussPts(),

                {{"T", vec_temp_ptr}},

                {{"U", u_ptr}, {"FLUX", mat_flux_ptr}},

                {{"PIOLA", mat_stress_ptr}},

                {{"HENCKY_STRAIN", mat_strain_ptr}}

                )

        );
      }

      MoFEMFunctionReturn(0);
    };

    auto domain_post_proc = [&]() {
      if (do_output_domain == PETSC_FALSE)
        return boost::shared_ptr<PostProcEle>();
      auto pp_fe = boost::make_shared<PostProcEle>(mField);
      CHK_MOAB_THROW(push_domain_ops(pp_fe),
                     "push domain ops to domain element");
      CHK_MOAB_THROW(push_post_proc_ops(pp_fe),
                     "push post proc ops to domain element");
      return pp_fe;
    };

    auto skin_post_proc = [&]() {
      if (do_output_skin == PETSC_FALSE)
        return boost::shared_ptr<SkinPostProcEle>();
      auto pp_fe = boost::make_shared<SkinPostProcEle>(mField);
      auto simple = mField.getInterface<Simple>();
      auto op_side = new OpLoopSide<SideEle>(mField, simple->getDomainFEName(),
                                             SPACE_DIM, Sev::verbose);
      CHK_MOAB_THROW(push_domain_ops(op_side),
                     "push domain ops to side element");
      pp_fe->getOpPtrVector().push_back(op_side);
      CHK_MOAB_THROW(push_post_proc_ops(pp_fe),
                     "push post proc ops to skin element");
      return pp_fe;
    };

    return std::make_pair(domain_post_proc(), skin_post_proc());
  };

  auto monitor_ptr = boost::make_shared<FEMethod>();

  auto set_time_monitor = [&](auto dm, auto solver, auto domain_post_proc_fe,
                              auto skin_post_proc_fe) {
    MoFEMFunctionBegin;
    monitor_ptr->preProcessHook = [&]() {
      MoFEMFunctionBegin;

      if (save_every && (monitor_ptr->ts_step % save_every == 0)) {
        if (do_output_domain) {
          CHKERR DMoFEMLoopFiniteElements(dm, simple->getDomainFEName(),
                                          domain_post_proc_fe,
                                          monitor_ptr->getCacheWeakPtr());
          CHKERR domain_post_proc_fe->writeFile(
              "out_" + boost::lexical_cast<std::string>(monitor_ptr->ts_step) +
              ".h5m");
        }
        if (do_output_skin) {
          CHKERR DMoFEMLoopFiniteElements(dm, simple->getBoundaryFEName(),
                                          skin_post_proc_fe,
                                          monitor_ptr->getCacheWeakPtr());
          CHKERR skin_post_proc_fe->writeFile(
              "out_skin_" +
              boost::lexical_cast<std::string>(monitor_ptr->ts_step) + ".h5m");
        }
      }

      if (doEvalField) {

        CHKERR mField.getInterface<FieldEvaluatorInterface>()
            ->evalFEAtThePoint<SPACE_DIM>(
                fieldEvalCoords.data(), 1e-12, simple->getProblemName(),
                simple->getDomainFEName(), fieldEvalData,
                mField.get_comm_rank(), mField.get_comm_rank(), nullptr,
                MF_EXIST, QUIET);

        if (atom_test) {
          auto eval_num_vec =
              createVectorMPI(mField.get_comm(), PETSC_DECIDE, 1);
          CHKERR VecZeroEntries(eval_num_vec);
          if (tempFieldPtr->size()) {
            CHKERR VecSetValue(eval_num_vec, 0, 1, ADD_VALUES);
          }
          CHKERR VecAssemblyBegin(eval_num_vec);
          CHKERR VecAssemblyEnd(eval_num_vec);

          double eval_num;
          CHKERR VecSum(eval_num_vec, &eval_num);
          if (!(int)eval_num) {
            SETERRQ1(PETSC_COMM_WORLD, MOFEM_ATOM_TEST_INVALID,
                     "atom test %d failed: did not find elements to evaluate "
                     "the field, check the coordinates",
                     atom_test);
          }
        }

        if (tempFieldPtr->size()) {
          auto t_temp = getFTensor0FromVec(*tempFieldPtr);
          MOFEM_LOG("ThermoElasticSync", Sev::inform)
              << "Eval point T: " << t_temp;
          if (atom_test && fabs(monitor_ptr->ts_t - 10) < 1e-12) {
            if (atom_test <= 3 && fabs(t_temp - 554.48) > 1e-2) {
              SETERRQ1(PETSC_COMM_WORLD, MOFEM_ATOM_TEST_INVALID,
                       "atom test %d failed: wrong temperature value",
                       atom_test);
            }
            if (atom_test == 4 && fabs(t_temp - 250) > 1e-2) {
              SETERRQ1(PETSC_COMM_WORLD, MOFEM_ATOM_TEST_INVALID,
                       "atom test %d failed: wrong temperature value",
                       atom_test);
            }
            if (atom_test == 5 && fabs(t_temp - 1) > 1e-2) {
              SETERRQ1(PETSC_COMM_WORLD, MOFEM_ATOM_TEST_INVALID,
                       "atom test %d failed: wrong temperature value",
                       atom_test);
            }
          }
        }
        if (fluxFieldPtr->size1()) {
          FTensor::Index<'i', SPACE_DIM> i;
          auto t_flux = getFTensor1FromMat<SPACE_DIM>(*fluxFieldPtr);
          auto flux_mag = sqrt(t_flux(i) * t_flux(i));
          MOFEM_LOG("ThermoElasticSync", Sev::inform)
              << "Eval point FLUX magnitude: " << flux_mag;
          if (atom_test && fabs(monitor_ptr->ts_t - 10) < 1e-12) {
            if (atom_test <= 3 && fabs(flux_mag - 27008.0) > 2e1) {
              SETERRQ1(PETSC_COMM_WORLD, MOFEM_ATOM_TEST_INVALID,
                       "atom test %d failed: wrong flux value", atom_test);
            }
            if (atom_test == 4 && fabs(flux_mag - 150e3) > 1e-1) {
              SETERRQ1(PETSC_COMM_WORLD, MOFEM_ATOM_TEST_INVALID,
                       "atom test %d failed: wrong flux value", atom_test);
            }
            if (atom_test == 5 && fabs(flux_mag) > 1e-6) {
              SETERRQ1(PETSC_COMM_WORLD, MOFEM_ATOM_TEST_INVALID,
                       "atom test %d failed: wrong flux value", atom_test);
            }
          }
        }
        if (dispFieldPtr->size1()) {
          FTensor::Index<'i', SPACE_DIM> i;
          auto t_disp = getFTensor1FromMat<SPACE_DIM>(*dispFieldPtr);
          auto disp_mag = sqrt(t_disp(i) * t_disp(i));
          MOFEM_LOG("ThermoElasticSync", Sev::inform)
              << "Eval point U magnitude: " << disp_mag;
          if (atom_test && fabs(monitor_ptr->ts_t - 10) < 1e-12) {
            if (atom_test == 1 && fabs(disp_mag - 0.00345) > 1e-5) {
              SETERRQ1(PETSC_COMM_WORLD, MOFEM_ATOM_TEST_INVALID,
                       "atom test %d failed: wrong displacement value",
                       atom_test);
            }
            if ((atom_test == 2 || atom_test == 3) &&
                fabs(disp_mag - 0.00265) > 1e-5) {
              SETERRQ1(PETSC_COMM_WORLD, MOFEM_ATOM_TEST_INVALID,
                       "atom test %d failed: wrong displacement value",
                       atom_test);
            }
            if ((atom_test == 5) &&
                fabs(t_disp(0) - (std::sqrt(std::exp(0.2)) - 1)) > 1e-5 &&
                fabs(t_disp(1) - (std::sqrt(std::exp(0.2)) - 1)) > 1e-5) {
              SETERRQ1(PETSC_COMM_WORLD, MOFEM_ATOM_TEST_INVALID,
                       "atom test %d failed: wrong displacement value",
                       atom_test);
            }
          }
        }
        if (strainFieldPtr->size1()) {
          FTensor::Index<'i', SPACE_DIM> i;
          auto t_strain =
              getFTensor2SymmetricFromMat<SPACE_DIM>(*strainFieldPtr);
          auto t_strain_trace = t_strain(i, i);
          if (atom_test && fabs(monitor_ptr->ts_t - 10) < 1e-12) {
            if (atom_test == 1 && fabs(t_strain_trace - 0.00679) > 1e-5) {
              SETERRQ1(PETSC_COMM_WORLD, MOFEM_ATOM_TEST_INVALID,
                       "atom test %d failed: wrong strain value", atom_test);
            }
            if ((atom_test == 2 || atom_test == 3) &&
                fabs(t_strain_trace - 0.00522) > 1e-5) {
              SETERRQ1(PETSC_COMM_WORLD, MOFEM_ATOM_TEST_INVALID,
                       "atom test %d failed: wrong strain value", atom_test);
            }
            if ((atom_test == 5) && fabs(t_strain_trace - 0.2) > 1e-5) {
              SETERRQ1(PETSC_COMM_WORLD, MOFEM_ATOM_TEST_INVALID,
                       "atom test %d failed: wrong strain value", atom_test);
            }
          }
        }
        if (stressFieldPtr->size1()) {
          double von_mises_stress;
          auto getVonMisesStress = [&](auto t_stress) {
            MoFEMFunctionBegin;
            von_mises_stress = std::sqrt(
                0.5 *
                ((t_stress(0, 0) - t_stress(1, 1)) *
                     (t_stress(0, 0) - t_stress(1, 1)) +
                 (SPACE_DIM == 3 ? (t_stress(1, 1) - t_stress(2, 2)) *
                                       (t_stress(1, 1) - t_stress(2, 2))
                                 : 0) +
                 (SPACE_DIM == 3 ? (t_stress(2, 2) - t_stress(0, 0)) *
                                       (t_stress(2, 2) - t_stress(0, 0))
                                 : 0) +
                 6 * (t_stress(0, 1) * t_stress(0, 1) +
                      (SPACE_DIM == 3 ? t_stress(1, 2) * t_stress(1, 2) : 0) +
                      (SPACE_DIM == 3 ? t_stress(2, 0) * t_stress(2, 0) : 0))));
            MoFEMFunctionReturn(0);
          };
          if (!IS_LARGE_STRAINS) {
            auto t_stress =
                getFTensor2SymmetricFromMat<SPACE_DIM>(*stressFieldPtr);
            CHKERR getVonMisesStress(t_stress);
          } else {
            auto t_stress =
                getFTensor2FromMat<SPACE_DIM, SPACE_DIM>(*stressFieldPtr);
            CHKERR getVonMisesStress(t_stress);
          }
          MOFEM_LOG("ThermoElasticSync", Sev::inform)
              << "Eval point von Mises Stress: " << von_mises_stress;
          if (atom_test && fabs(monitor_ptr->ts_t - 10) < 1e-12) {
            if (atom_test == 1 && fabs(von_mises_stress - 523.0) > 5e-1) {
              SETERRQ1(PETSC_COMM_WORLD, MOFEM_ATOM_TEST_INVALID,
                       "atom test %d failed: wrong von Mises stress value",
                       atom_test);
            }
            if (atom_test == 2 && fabs(von_mises_stress - 16.3) > 5e-2) {
              SETERRQ1(PETSC_COMM_WORLD, MOFEM_ATOM_TEST_INVALID,
                       "atom test %d failed: wrong von Mises stress value",
                       atom_test);
            }
            if (atom_test == 3 && fabs(von_mises_stress - 14.9) > 5e-2) {
              SETERRQ1(PETSC_COMM_WORLD, MOFEM_ATOM_TEST_INVALID,
                       "atom test %d failed: wrong von Mises stress value",
                       atom_test);
            }
            if (atom_test == 5 && fabs(von_mises_stress) > 5e-2) {
              SETERRQ1(PETSC_COMM_WORLD, MOFEM_ATOM_TEST_INVALID,
                       "atom test %d failed: wrong von Mises stress value",
                       atom_test);
            }
          }
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
      SmartPetscObj<IS> is_TBC;
      CHKERR is_mng->isCreateProblemFieldAndRank(name_prb, ROW, "TBC", 0, 0,
                                                 is_TBC);
      IS is_tmp, is_tmp2;
      CHKERR ISExpand(is_T, is_flux, &is_tmp);
      CHKERR ISExpand(is_TBC, is_tmp, &is_tmp2);
      CHKERR ISDestroy(&is_tmp);
      auto is_TFlux = SmartPetscObj<IS>(is_tmp2);

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

  auto [domain_post_proc_fe, skin_post_proc_fe] =
      create_post_process_elements();
  CHKERR set_time_monitor(dm, solver, domain_post_proc_fe, skin_post_proc_fe);

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
