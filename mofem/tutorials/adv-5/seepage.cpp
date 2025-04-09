/**
 * \file seepage.cpp
 * \example seepage.cpp
 *
 * Hydromechanical problem
 *
 */

/* This file is part of MoFEM.
 * MoFEM is free software: you can redistribute it and/or modify it under
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
 * License along with MoFEM. If not, see <http://www.gnu.org/licenses/>. */

#ifndef EXECUTABLE_DIMENSION
  #define EXECUTABLE_DIMENSION 2
#endif

#include <MoFEM.hpp>
#include <MatrixFunction.hpp>
#include <IntegrationRules.hpp>

using namespace MoFEM;

constexpr int SPACE_DIM =
    EXECUTABLE_DIMENSION; //< Space dimension of problem, mesh

using EntData = EntitiesFieldData::EntData;
using DomainEle = PipelineManager::ElementsAndOpsByDim<SPACE_DIM>::DomainEle;
using DomainEleOp = DomainEle::UserDataOperator;
using BoundaryEle =
    PipelineManager::ElementsAndOpsByDim<SPACE_DIM>::BoundaryEle;
using BoundaryEleOp = BoundaryEle::UserDataOperator;
using PostProcEle = PostProcBrokenMeshInMoab<DomainEle>;

using AssemblyDomainEleOp =
    FormsIntegrators<DomainEleOp>::Assembly<PETSC>::OpBase;

//! [Only used with Hooke equation (linear material model)]
using OpKCauchy = FormsIntegrators<DomainEleOp>::Assembly<PETSC>::BiLinearForm<
    GAUSS>::OpGradSymTensorGrad<1, SPACE_DIM, SPACE_DIM, 0>;
using OpInternalForceCauchy = FormsIntegrators<DomainEleOp>::Assembly<
    PETSC>::LinearForm<GAUSS>::OpGradTimesSymTensor<1, SPACE_DIM, SPACE_DIM>;
//! [Only used with Hooke equation (linear material model)]

//! [Only used for dynamics]
using OpMass = FormsIntegrators<DomainEleOp>::Assembly<PETSC>::BiLinearForm<
    GAUSS>::OpMass<1, SPACE_DIM>;
using OpInertiaForce = FormsIntegrators<DomainEleOp>::Assembly<
    PETSC>::LinearForm<GAUSS>::OpBaseTimesVector<1, SPACE_DIM, 1>;
//! [Only used for dynamics]

//! [Only used with Hencky/nonlinear material]
using OpKPiola = FormsIntegrators<DomainEleOp>::Assembly<PETSC>::BiLinearForm<
    GAUSS>::OpGradTensorGrad<1, SPACE_DIM, SPACE_DIM, 1>;
using OpInternalForcePiola = FormsIntegrators<DomainEleOp>::Assembly<
    PETSC>::LinearForm<GAUSS>::OpGradTimesTensor<1, SPACE_DIM, SPACE_DIM>;
//! [Only used with Hencky/nonlinear material]

//! [Essential boundary conditions]
using OpBoundaryMass = FormsIntegrators<BoundaryEleOp>::Assembly<
    PETSC>::BiLinearForm<GAUSS>::OpMass<1, SPACE_DIM>;
using OpBoundaryVec = FormsIntegrators<BoundaryEleOp>::Assembly<
    PETSC>::LinearForm<GAUSS>::OpBaseTimesVector<1, SPACE_DIM, 0>;
using OpBoundaryInternal = FormsIntegrators<BoundaryEleOp>::Assembly<
    PETSC>::LinearForm<GAUSS>::OpBaseTimesVector<1, SPACE_DIM, 1>;
//! [Essential boundary conditions]

using OpBaseDivU = FormsIntegrators<DomainEleOp>::Assembly<PETSC>::BiLinearForm<
    GAUSS>::OpMixScalarTimesDiv<SPACE_DIM>;

// Thermal operators
/**
 * @brief Integrate Lhs base of flux (1/k) base of flux (FLUX x FLUX)
 *
 */
using OpHdivHdiv = FormsIntegrators<DomainEleOp>::Assembly<PETSC>::BiLinearForm<
    GAUSS>::OpMass<3, 3>;

/**
 * @brief Integrate Lhs div of base of flux time base of temperature (FLUX x T)
 * and transpose of it, i.e. (T x FLAX)
 *
 */
using OpHdivQ = FormsIntegrators<DomainEleOp>::Assembly<PETSC>::BiLinearForm<
    GAUSS>::OpMixDivTimesScalar<SPACE_DIM>;

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
using OpHdivFlux = FormsIntegrators<DomainEleOp>::Assembly<PETSC>::LinearForm<
    GAUSS>::OpBaseTimesVector<3, 3, 1>;

/**
 * @brief  Integrate Rhs div flux base times temperature (T)
 *
 */
using OpHDivH = FormsIntegrators<DomainEleOp>::Assembly<PETSC>::LinearForm<
    GAUSS>::OpMixDivTimesU<3, 1, 2>;

/**
 * @brief Integrate Rhs base of temperature time heat capacity times heat rate
 * (T)
 *
 */
using OpBaseDotH = FormsIntegrators<DomainEleOp>::Assembly<PETSC>::LinearForm<
    GAUSS>::OpBaseTimesScalarField<1>;

/**
 * @brief Integrate Rhs base of temperature times divergent of flux (T)
 *
 */
using OpBaseDivFlux = OpBaseDotH;

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

double scale = 1.;

double default_young_modulus = 1;
double default_poisson_ratio = 0.25; // zero for verification
double default_conductivity = 1;
double fluid_density = 10;

// #include <OpPostProcElastic.hpp>
#include <SeepageOps.hpp>

auto save_range = [](moab::Interface &moab, const std::string name,
                     const Range r) {
  MoFEMFunctionBegin;
  auto out_meshset = get_temp_meshset_ptr(moab);
  CHKERR moab.add_entities(*out_meshset, r);
  CHKERR moab.write_file(name.c_str(), "VTK", "", out_meshset->get_ptr(), 1);
  MoFEMFunctionReturn(0);
};

struct Seepage {

  Seepage(MoFEM::Interface &m_field) : mField(m_field) {}

  MoFEMErrorCode runProblem();

private:
  MoFEM::Interface &mField;

  MoFEMErrorCode setupProblem();
  MoFEMErrorCode createCommonData();
  MoFEMErrorCode bC();
  MoFEMErrorCode OPs();
  MoFEMErrorCode tsSolve();

  PetscBool doEvalField;
  std::array<double, SPACE_DIM> fieldEvalCoords;

  std::tuple<SmartPetscObj<Vec>, SmartPetscObj<VecScatter>> uXScatter;
  std::tuple<SmartPetscObj<Vec>, SmartPetscObj<VecScatter>> uYScatter;
  std::tuple<SmartPetscObj<Vec>, SmartPetscObj<VecScatter>> uZScatter;

  struct BlockedParameters
      : public boost::enable_shared_from_this<BlockedParameters> {
    MatrixDouble mD;
    double fluidConductivity;

    inline auto getDPtr() {
      return boost::shared_ptr<MatrixDouble>(shared_from_this(), &mD);
    }

    inline auto getConductivityPtr() {
      return boost::shared_ptr<double>(shared_from_this(), &fluidConductivity);
    }
  };

  MoFEMErrorCode addMatBlockOps(
      boost::ptr_deque<ForcesAndSourcesCore::UserDataOperator> &pipeline,
      std::string block_elastic_name, std::string block_thermal_name,
      boost::shared_ptr<BlockedParameters> blockedParamsPtr, Sev sev);
};

MoFEMErrorCode Seepage::addMatBlockOps(
    boost::ptr_deque<ForcesAndSourcesCore::UserDataOperator> &pipeline,
    std::string block_elastic_name, std::string block_thermal_name,
    boost::shared_ptr<BlockedParameters> blockedParamsPtr, Sev sev) {
  MoFEMFunctionBegin;

  struct OpMatElasticBlocks : public DomainEleOp {
    OpMatElasticBlocks(boost::shared_ptr<MatrixDouble> m, double bulk_modulus_K,
                       double shear_modulus_G, MoFEM::Interface &m_field,
                       Sev sev,
                       std::vector<const CubitMeshSets *> meshset_vec_ptr)
        : DomainEleOp(NOSPACE, DomainEleOp::OPSPACE), matDPtr(m),
          bulkModulusKDefault(bulk_modulus_K),
          shearModulusGDefault(shear_modulus_G) {
      CHK_THROW_MESSAGE(extractElasticBlockData(m_field, meshset_vec_ptr, sev),
                        "Can not get data from block");
    }

    MoFEMErrorCode doWork(int side, EntityType type,
                          EntitiesFieldData::EntData &data) {
      MoFEMFunctionBegin;

      for (auto &b : blockData) {

        if (b.blockEnts.find(getFEEntityHandle()) != b.blockEnts.end()) {
          CHKERR getMatDPtr(matDPtr, b.bulkModulusK, b.shearModulusG);
          MoFEMFunctionReturnHot(0);
        }
      }

      CHKERR getMatDPtr(matDPtr, bulkModulusKDefault, shearModulusGDefault);
      MoFEMFunctionReturn(0);
    }

  private:
    boost::shared_ptr<MatrixDouble> matDPtr;

    struct BlockData {
      double bulkModulusK;
      double shearModulusG;
      Range blockEnts;
    };

    double bulkModulusKDefault;
    double shearModulusGDefault;
    std::vector<BlockData> blockData;

    MoFEMErrorCode
    extractElasticBlockData(MoFEM::Interface &m_field,
                            std::vector<const CubitMeshSets *> meshset_vec_ptr,
                            Sev sev) {
      MoFEMFunctionBegin;

      for (auto m : meshset_vec_ptr) {
        MOFEM_TAG_AND_LOG("WORLD", sev, "Mat Elastic Block") << *m;
        std::vector<double> block_data;
        CHKERR m->getAttributes(block_data);
        if (block_data.size() < 2) {
          SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
                  "Expected that block has two attributes");
        }
        auto get_block_ents = [&]() {
          Range ents;
          CHKERR
          m_field.get_moab().get_entities_by_handle(m->meshset, ents, true);
          return ents;
        };

        double young_modulus = block_data[0];
        double poisson_ratio = block_data[1];
        double bulk_modulus_K = young_modulus / (3 * (1 - 2 * poisson_ratio));
        double shear_modulus_G = young_modulus / (2 * (1 + poisson_ratio));

        MOFEM_TAG_AND_LOG("WORLD", sev, "Mat Elastic Block")
            << m->getName() << ": E = " << young_modulus
            << " nu = " << poisson_ratio;

        blockData.push_back(
            {bulk_modulus_K, shear_modulus_G, get_block_ents()});
      }
      MOFEM_LOG_CHANNEL("WORLD");
      MoFEMFunctionReturn(0);
    }

    MoFEMErrorCode getMatDPtr(boost::shared_ptr<MatrixDouble> mat_D_ptr,
                              double bulk_modulus_K, double shear_modulus_G) {
      MoFEMFunctionBegin;
      //! [Calculate elasticity tensor]
      auto set_material_stiffness = [&]() {
        FTensor::Index<'i', SPACE_DIM> i;
        FTensor::Index<'j', SPACE_DIM> j;
        FTensor::Index<'k', SPACE_DIM> k;
        FTensor::Index<'l', SPACE_DIM> l;
        constexpr auto t_kd = FTensor::Kronecker_Delta_symmetric<int>();
        double A = (SPACE_DIM == 2)
                       ? 2 * shear_modulus_G /
                             (bulk_modulus_K + (4. / 3.) * shear_modulus_G)
                       : 1;
        auto t_D = getFTensor4DdgFromMat<SPACE_DIM, SPACE_DIM, 0>(*mat_D_ptr);
        t_D(i, j, k, l) =
            2 * shear_modulus_G * ((t_kd(i, k) ^ t_kd(j, l)) / 4.) +
            A * (bulk_modulus_K - (2. / 3.) * shear_modulus_G) * t_kd(i, j) *
                t_kd(k, l);
      };
      //! [Calculate elasticity tensor]
      constexpr auto size_symm = (SPACE_DIM * (SPACE_DIM + 1)) / 2;
      mat_D_ptr->resize(size_symm * size_symm, 1);
      set_material_stiffness();
      MoFEMFunctionReturn(0);
    }
  };

  double default_bulk_modulus_K =
      default_young_modulus / (3 * (1 - 2 * default_poisson_ratio));
  double default_shear_modulus_G =
      default_young_modulus / (2 * (1 + default_poisson_ratio));

  pipeline.push_back(new OpMatElasticBlocks(
      blockedParamsPtr->getDPtr(), default_bulk_modulus_K,
      default_shear_modulus_G, mField, sev,

      // Get blockset using regular expression
      mField.getInterface<MeshsetsManager>()->getCubitMeshsetPtr(std::regex(

          (boost::format("%s(.*)") % block_elastic_name).str()

              ))

          ));

  struct OpMatFluidBlocks : public DomainEleOp {
    OpMatFluidBlocks(boost::shared_ptr<double> conductivity_ptr,
                     MoFEM::Interface &m_field, Sev sev,
                     std::vector<const CubitMeshSets *> meshset_vec_ptr)
        : DomainEleOp(NOSPACE, DomainEleOp::OPSPACE),
          conductivityPtr(conductivity_ptr) {
      CHK_THROW_MESSAGE(extractThermalBlockData(m_field, meshset_vec_ptr, sev),
                        "Can not get data from block");
    }

    MoFEMErrorCode doWork(int side, EntityType type,
                          EntitiesFieldData::EntData &data) {
      MoFEMFunctionBegin;

      for (auto &b : blockData) {

        if (b.blockEnts.find(getFEEntityHandle()) != b.blockEnts.end()) {
          *conductivityPtr = b.conductivity;
          MoFEMFunctionReturnHot(0);
        }
      }

      *conductivityPtr = default_conductivity;

      MoFEMFunctionReturn(0);
    }

  private:
    struct BlockData {
      double conductivity;
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
        if (block_data.size() < 1) {
          SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
                  "Expected that block has two attributes");
        }
        auto get_block_ents = [&]() {
          Range ents;
          CHKERR
          m_field.get_moab().get_entities_by_handle(m->meshset, ents, true);
          return ents;
        };

        MOFEM_TAG_AND_LOG("WORLD", sev, "Mat Thermal Block")
            << m->getName() << ": conductivity = " << block_data[0];

        blockData.push_back({block_data[0], get_block_ents()});
      }
      MOFEM_LOG_CHANNEL("WORLD");
      MoFEMFunctionReturn(0);
    }

    boost::shared_ptr<double> expansionPtr;
    boost::shared_ptr<double> conductivityPtr;
    boost::shared_ptr<double> capacityPtr;
  };

  pipeline.push_back(new OpMatFluidBlocks(
      blockedParamsPtr->getConductivityPtr(), mField, sev,

      // Get blockset using regular expression
      mField.getInterface<MeshsetsManager>()->getCubitMeshsetPtr(std::regex(

          (boost::format("%s(.*)") % block_thermal_name).str()

              ))

          ));

  MoFEMFunctionReturn(0);
}

//! [Run problem]
MoFEMErrorCode Seepage::runProblem() {
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
MoFEMErrorCode Seepage::setupProblem() {
  MoFEMFunctionBegin;
  Simple *simple = mField.getInterface<Simple>();
  // Add field
  constexpr FieldApproximationBase base = AINSWORTH_LEGENDRE_BASE;
  // Mechanical fields
  CHKERR simple->addDomainField("U", H1, base, SPACE_DIM);
  CHKERR simple->addBoundaryField("U", H1, base, SPACE_DIM);
  // Temerature
  const auto flux_space = (SPACE_DIM == 2) ? HCURL : HDIV;
  CHKERR simple->addDomainField("H", L2, AINSWORTH_LEGENDRE_BASE, 1);
  CHKERR simple->addDomainField("FLUX", flux_space, DEMKOWICZ_JACOBI_BASE, 1);
  CHKERR simple->addBoundaryField("FLUX", flux_space, DEMKOWICZ_JACOBI_BASE, 1);

  int order = 2.;
  CHKERR PetscOptionsGetInt(PETSC_NULL, "", "-order", &order, PETSC_NULL);
  CHKERR simple->setFieldOrder("U", order);
  CHKERR simple->setFieldOrder("H", order - 1);
  CHKERR simple->setFieldOrder("FLUX", order);

  CHKERR simple->setUp();

  MoFEMFunctionReturn(0);
}
//! [Set up problem]

//! [Create common data]
MoFEMErrorCode Seepage::createCommonData() {
  MoFEMFunctionBegin;

  auto get_command_line_parameters = [&]() {
    MoFEMFunctionBegin;
    CHKERR PetscOptionsGetScalar(PETSC_NULL, "", "-young_modulus",
                                 &default_young_modulus, PETSC_NULL);
    CHKERR PetscOptionsGetScalar(PETSC_NULL, "", "-poisson_ratio",
                                 &default_poisson_ratio, PETSC_NULL);
    CHKERR PetscOptionsGetScalar(PETSC_NULL, "", "-conductivity",
                                 &default_conductivity, PETSC_NULL);
    CHKERR PetscOptionsGetScalar(PETSC_NULL, "", "-fluid_density",
                                 &fluid_density, PETSC_NULL);

    MOFEM_LOG("Seepage", Sev::inform)
        << "Default Young modulus " << default_young_modulus;
    MOFEM_LOG("Seepage", Sev::inform)
        << "Default Poisson ratio " << default_poisson_ratio;
    MOFEM_LOG("Seepage", Sev::inform)
        << "Default conductivity " << default_conductivity;
    MOFEM_LOG("Seepage", Sev::inform) << "Fluid denisty " << fluid_density;

    int coords_dim = SPACE_DIM;
    CHKERR PetscOptionsGetRealArray(NULL, NULL, "-field_eval_coords",
                                    fieldEvalCoords.data(), &coords_dim,
                                    &doEvalField);

    MoFEMFunctionReturn(0);
  };

  CHKERR get_command_line_parameters();

  MoFEMFunctionReturn(0);
}
//! [Create common data]

//! [Boundary condition]
MoFEMErrorCode Seepage::bC() {
  MoFEMFunctionBegin;

  MOFEM_LOG("SYNC", Sev::noisy) << "bC";
  MOFEM_LOG_SEVERITY_SYNC(mField.get_comm(), Sev::noisy);

  auto simple = mField.getInterface<Simple>();
  auto bc_mng = mField.getInterface<BcManager>();

  CHKERR bc_mng->pushMarkDOFsOnEntities<DisplacementCubitBcData>(
      simple->getProblemName(), "U");
  CHKERR bc_mng->pushMarkDOFsOnEntities<HeatFluxCubitBcData>(
      simple->getProblemName(), "FLUX", false);

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
    CHK_THROW_MESSAGE(remove_named_blocks("PRESSURE"), "remove_named_blocks");
    CHK_THROW_MESSAGE(remove_named_blocks("FLUIDFLUX"), "remove_named_blocks");

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

  MoFEMFunctionReturn(0);
}
//! [Boundary condition]

//! [Push operators to pipeline]
MoFEMErrorCode Seepage::OPs() {
  MoFEMFunctionBegin;
  auto pipeline_mng = mField.getInterface<PipelineManager>();
  auto simple = mField.getInterface<Simple>();
  auto bc_mng = mField.getInterface<BcManager>();

  auto boundary_marker =
      bc_mng->getMergedBlocksMarker(vector<string>{"FLUIDFLUX"});

  auto u_grad_ptr = boost::make_shared<MatrixDouble>();
  auto dot_u_grad_ptr = boost::make_shared<MatrixDouble>();
  auto trace_dot_u_grad_ptr = boost::make_shared<VectorDouble>();
  auto h_ptr = boost::make_shared<VectorDouble>();
  auto dot_h_ptr = boost::make_shared<VectorDouble>();
  auto flux_ptr = boost::make_shared<MatrixDouble>();
  auto div_flux_ptr = boost::make_shared<VectorDouble>();
  auto strain_ptr = boost::make_shared<MatrixDouble>();
  auto stress_ptr = boost::make_shared<MatrixDouble>();

  auto time_scale = boost::make_shared<TimeScale>();

  auto block_params = boost::make_shared<BlockedParameters>();
  auto mDPtr = block_params->getDPtr();
  auto conductivity_ptr = block_params->getConductivityPtr();

  auto integration_rule = [](int, int, int approx_order) {
    return 2 * approx_order;
  };

  CHKERR pipeline_mng->setDomainRhsIntegrationRule(integration_rule);
  CHKERR pipeline_mng->setDomainLhsIntegrationRule(integration_rule);
  CHKERR pipeline_mng->setBoundaryLhsIntegrationRule(integration_rule);
  CHKERR pipeline_mng->setBoundaryRhsIntegrationRule(integration_rule);

  auto add_domain_base_ops = [&](auto &pip) {
    MoFEMFunctionBegin;

    CHKERR addMatBlockOps(pip, "MAT_ELASTIC", "MAT_FLUID", block_params,
                          Sev::inform);
    CHKERR AddHOOps<SPACE_DIM, SPACE_DIM, SPACE_DIM>::add(pip, {H1, HDIV});

    pip.push_back(new OpCalculateVectorFieldGradient<SPACE_DIM, SPACE_DIM>(
        "U", u_grad_ptr));
    pip.push_back(new OpSymmetrizeTensor<SPACE_DIM>(u_grad_ptr, strain_ptr));
    pip.push_back(new OpCalculateVectorFieldGradientDot<SPACE_DIM, SPACE_DIM>(
        "U", dot_u_grad_ptr));
    pip.push_back(new OpCalculateTraceFromMat<SPACE_DIM>(dot_u_grad_ptr,
                                                         trace_dot_u_grad_ptr));

    pip.push_back(new OpCalculateScalarFieldValues("H", h_ptr));
    pip.push_back(new OpCalculateScalarFieldValuesDot("H", dot_h_ptr));
    pip.push_back(new OpCalculateHVecVectorField<3>("FLUX", flux_ptr));
    pip.push_back(new OpCalculateHdivVectorDivergence<3, SPACE_DIM>(
        "FLUX", div_flux_ptr));

    MoFEMFunctionReturn(0);
  };

  auto add_domain_ops_lhs_mechanical = [&](auto &pip) {
    MoFEMFunctionBegin;
    pip.push_back(new OpKCauchy("U", "U", mDPtr));
    pip.push_back(new OpBaseDivU(
        "H", "U",
        [](const double, const double, const double) { return -9.81; }, true,
        true));
    MoFEMFunctionReturn(0);
  };

  auto add_domain_ops_rhs_mechanical = [&](auto &pip) {
    MoFEMFunctionBegin;

    CHKERR DomainNaturalBCRhs::AddFluxToPipeline<OpBodyForce>::add(
        pip, mField, "U", {time_scale}, "BODY_FORCE", Sev::inform);

    // Calculate internal forece
    pip.push_back(new OpTensorTimesSymmetricTensor<SPACE_DIM, SPACE_DIM>(
        strain_ptr, stress_ptr, mDPtr));
    pip.push_back(new OpInternalForceCauchy("U", stress_ptr));
    pip.push_back(
        new SeepageOps::OpDomainRhsHydrostaticStress<SPACE_DIM>("U", h_ptr));

    MoFEMFunctionReturn(0);
  };

  auto add_domain_ops_lhs_seepage = [&](auto &pip, auto &fe) {
    MoFEMFunctionBegin;
    auto resistance = [conductivity_ptr](const double, const double,
                                         const double) {
      return (1. / *(conductivity_ptr));
    };

    auto unity = []() constexpr { return -1; };
    pip.push_back(new OpHdivHdiv("FLUX", "FLUX", resistance));
    pip.push_back(new OpHdivQ("FLUX", "H", unity, true));
    auto op_base_div_u = new OpBaseDivU(
        "H", "U", [](double, double, double) constexpr { return -1; }, false,
        false);
    op_base_div_u->feScalingFun = [](const FEMethod *fe_ptr) {
      return fe_ptr->ts_a;
    };
    pip.push_back(op_base_div_u);

    MoFEMFunctionReturn(0);
  };

  auto add_domain_ops_rhs_seepage = [&](auto &pip) {
    MoFEMFunctionBegin;
    auto resistance = [conductivity_ptr](double, double, double) {
      return (1. / (*conductivity_ptr));
    };
    auto minus_one = [](double, double, double) constexpr { return -1; };

    pip.push_back(new OpHdivFlux("FLUX", flux_ptr, resistance));
    pip.push_back(new OpHDivH("FLUX", h_ptr, minus_one));
    pip.push_back(new OpBaseDotH("H", trace_dot_u_grad_ptr, minus_one));
    pip.push_back(new OpBaseDivFlux("H", div_flux_ptr, minus_one));

    MoFEMFunctionReturn(0);
  };

  auto add_boundary_rhs_ops = [&](auto &pip) {
    MoFEMFunctionBegin;

    CHKERR AddHOOps<SPACE_DIM - 1, SPACE_DIM, SPACE_DIM>::add(pip, {HDIV});

    pip.push_back(new OpSetBc("FLUX", true, boundary_marker));

    CHKERR BoundaryNaturalBC::AddFluxToPipeline<OpForce>::add(
        pipeline_mng->getOpBoundaryRhsPipeline(), mField, "U", {time_scale},
        "FORCE", Sev::inform);

    CHKERR BoundaryNaturalBC::AddFluxToPipeline<OpTemperatureBC>::add(
        pipeline_mng->getOpBoundaryRhsPipeline(), mField, "FLUX", {time_scale},
        "PRESSURE", Sev::inform);

    pip.push_back(new OpUnSetBc("FLUX"));

    auto mat_flux_ptr = boost::make_shared<MatrixDouble>();
    pip.push_back(
        new OpCalculateHVecVectorField<3, SPACE_DIM>("FLUX", mat_flux_ptr));
    CHKERR EssentialBC<BoundaryEleOp>::Assembly<PETSC>::LinearForm<GAUSS>::
        AddEssentialToPipeline<OpEssentialFluxRhs>::add(
            mField, pip, simple->getProblemName(), "FLUX", mat_flux_ptr,
            {time_scale});

    MoFEMFunctionReturn(0);
  };

  auto add_boundary_lhs_ops = [&](auto &pip) {
    MoFEMFunctionBegin;

    CHKERR AddHOOps<SPACE_DIM - 1, SPACE_DIM, SPACE_DIM>::add(pip, {HDIV});

    CHKERR EssentialBC<BoundaryEleOp>::Assembly<PETSC>::BiLinearForm<GAUSS>::
        AddEssentialToPipeline<OpEssentialFluxLhs>::add(
            mField, pip, simple->getProblemName(), "FLUX");

    MoFEMFunctionReturn(0);
  };

  // LHS
  CHKERR add_domain_base_ops(pipeline_mng->getOpDomainLhsPipeline());
  CHKERR add_domain_ops_lhs_mechanical(pipeline_mng->getOpDomainLhsPipeline());
  CHKERR add_domain_ops_lhs_seepage(pipeline_mng->getOpDomainLhsPipeline(),
                                    pipeline_mng->getDomainLhsFE());

  // RHS
  CHKERR add_domain_base_ops(pipeline_mng->getOpDomainRhsPipeline());
  CHKERR add_domain_ops_rhs_mechanical(pipeline_mng->getOpDomainRhsPipeline());
  CHKERR add_domain_ops_rhs_seepage(pipeline_mng->getOpDomainRhsPipeline());

  CHKERR add_boundary_rhs_ops(pipeline_mng->getOpBoundaryRhsPipeline());
  CHKERR add_boundary_lhs_ops(pipeline_mng->getOpBoundaryLhsPipeline());

  MoFEMFunctionReturn(0);
}
//! [Push operators to pipeline]

//! [Solve]
MoFEMErrorCode Seepage::tsSolve() {
  MoFEMFunctionBegin;
  auto simple = mField.getInterface<Simple>();
  PipelineManager *pipeline_mng = mField.getInterface<PipelineManager>();

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
    auto mDPtr = block_params->getDPtr();
    CHKERR addMatBlockOps(post_proc_fe->getOpPtrVector(), "MAT_ELASTIC",
                          "MAT_FLUID", block_params, Sev::verbose);
    CHKERR AddHOOps<SPACE_DIM, SPACE_DIM, SPACE_DIM>::add(
        post_proc_fe->getOpPtrVector(), {H1, HDIV});

    auto mat_grad_ptr = boost::make_shared<MatrixDouble>();
    auto mat_strain_ptr = boost::make_shared<MatrixDouble>();
    auto mat_stress_ptr = boost::make_shared<MatrixDouble>();

    auto h_ptr = boost::make_shared<VectorDouble>();
    auto mat_flux_ptr = boost::make_shared<MatrixDouble>();

    post_proc_fe->getOpPtrVector().push_back(
        new OpCalculateScalarFieldValues("H", h_ptr));
    post_proc_fe->getOpPtrVector().push_back(
        new OpCalculateHVecVectorField<3, SPACE_DIM>("FLUX", mat_flux_ptr));

    auto u_ptr = boost::make_shared<MatrixDouble>();
    post_proc_fe->getOpPtrVector().push_back(
        new OpCalculateVectorFieldValues<SPACE_DIM>("U", u_ptr));
    post_proc_fe->getOpPtrVector().push_back(
        new OpCalculateVectorFieldGradient<SPACE_DIM, SPACE_DIM>("U",
                                                                 mat_grad_ptr));
    post_proc_fe->getOpPtrVector().push_back(
        new OpSymmetrizeTensor<SPACE_DIM>(mat_grad_ptr, mat_strain_ptr));
    post_proc_fe->getOpPtrVector().push_back(
        new OpTensorTimesSymmetricTensor<SPACE_DIM, SPACE_DIM>(
            mat_strain_ptr, mat_stress_ptr, mDPtr));

    using OpPPMap = OpPostProcMapInMoab<SPACE_DIM, SPACE_DIM>;

    post_proc_fe->getOpPtrVector().push_back(

        new OpPPMap(

            post_proc_fe->getPostProcMesh(), post_proc_fe->getMapGaussPts(),

            {{"H", h_ptr}},

            {{"U", u_ptr}, {"FLUX", mat_flux_ptr}},

            {},

            {{"STRAIN", mat_strain_ptr}, {"STRESS", mat_stress_ptr}}

            )

    );

    return post_proc_fe;
  };

  auto create_creaction_fe = [&]() {
    auto fe_ptr = boost::make_shared<DomainEle>(mField);
    fe_ptr->getRuleHook = [](int, int, int o) { return 2 * o; };

    auto &pip = fe_ptr->getOpPtrVector();

    auto block_params = boost::make_shared<BlockedParameters>();
    auto mDPtr = block_params->getDPtr();
    CHKERR addMatBlockOps(pip, "MAT_ELASTIC", "MAT_FLUID", block_params,
                          Sev::verbose);
    CHKERR AddHOOps<SPACE_DIM, SPACE_DIM, SPACE_DIM>::add(pip, {H1, HDIV});

    auto u_grad_ptr = boost::make_shared<MatrixDouble>();
    auto strain_ptr = boost::make_shared<MatrixDouble>();
    auto stress_ptr = boost::make_shared<MatrixDouble>();

    pip.push_back(new OpCalculateVectorFieldGradient<SPACE_DIM, SPACE_DIM>(
        "U", u_grad_ptr));
    pip.push_back(new OpSymmetrizeTensor<SPACE_DIM>(u_grad_ptr, strain_ptr));

    // Calculate internal forece
    pip.push_back(new OpTensorTimesSymmetricTensor<SPACE_DIM, SPACE_DIM>(
        strain_ptr, stress_ptr, mDPtr));
    pip.push_back(new OpInternalForceCauchy("U", stress_ptr));

    fe_ptr->postProcessHook =
        EssentialPreProcReaction<DisplacementCubitBcData>(mField, fe_ptr);

    return fe_ptr;
  };

  auto monitor_ptr = boost::make_shared<FEMethod>();
  auto post_proc_fe = create_post_process_element();
  auto res = createDMVector(dm);
  auto rections_fe = create_creaction_fe();

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

      rections_fe->f = res;
      CHKERR DMoFEMLoopFiniteElements(dm, simple->getDomainFEName(),
                                      rections_fe,
                                      monitor_ptr->getCacheWeakPtr());

      double nrm;
      CHKERR VecNorm(res, NORM_2, &nrm);
      MOFEM_LOG("Seepage", Sev::verbose)
          << "Residual norm " << nrm << " at time step "
          << monitor_ptr->ts_step;

      if (doEvalField) {

        auto scalar_field_ptr = boost::make_shared<VectorDouble>();
        auto vector_field_ptr = boost::make_shared<MatrixDouble>();
        auto tensor_field_ptr = boost::make_shared<MatrixDouble>();

        auto field_eval_data = mField.getInterface<FieldEvaluatorInterface>()
                                   ->getData<DomainEle>();

        CHKERR mField.getInterface<FieldEvaluatorInterface>()
            ->buildTree<SPACE_DIM>(field_eval_data, simple->getDomainFEName());

        field_eval_data->setEvalPoints(fieldEvalCoords.data(), 1);
        auto no_rule = [](int, int, int) { return -1; };

        auto field_eval_ptr = field_eval_data->feMethodPtr.lock();
        field_eval_ptr->getRuleHook = no_rule;

        auto u_grad_ptr = boost::make_shared<MatrixDouble>();
        auto strain_ptr = boost::make_shared<MatrixDouble>();
        auto stress_ptr = boost::make_shared<MatrixDouble>();
        auto h_ptr = boost::make_shared<VectorDouble>();

        auto block_params = boost::make_shared<BlockedParameters>();
        auto mDPtr = block_params->getDPtr();
        CHKERR addMatBlockOps(field_eval_ptr->getOpPtrVector(), "MAT_ELASTIC",
                              "MAT_FLUID", block_params, Sev::noisy);
        CHKERR AddHOOps<SPACE_DIM, SPACE_DIM, SPACE_DIM>::add(
            field_eval_ptr->getOpPtrVector(), {H1, HDIV});
        field_eval_ptr->getOpPtrVector().push_back(
            new OpCalculateVectorFieldGradient<SPACE_DIM, SPACE_DIM>(
                "U", u_grad_ptr));
        field_eval_ptr->getOpPtrVector().push_back(
            new OpSymmetrizeTensor<SPACE_DIM>(u_grad_ptr, strain_ptr));
        field_eval_ptr->getOpPtrVector().push_back(
            new OpCalculateScalarFieldValues("H", h_ptr));
        field_eval_ptr->getOpPtrVector().push_back(
            new OpTensorTimesSymmetricTensor<SPACE_DIM, SPACE_DIM>(
                strain_ptr, stress_ptr, mDPtr));

        CHKERR mField.getInterface<FieldEvaluatorInterface>()
            ->evalFEAtThePoint<SPACE_DIM>(
                fieldEvalCoords.data(), 1e-12, simple->getProblemName(),
                simple->getDomainFEName(), field_eval_data,
                mField.get_comm_rank(), mField.get_comm_rank(), nullptr,
                MF_EXIST, QUIET);

        MOFEM_LOG("SeepageSync", Sev::inform)
            << "Eval point hydrostatic hight: " << *h_ptr;
        MOFEM_LOG("SeepageSync", Sev::inform)
            << "Eval point skeleton stress pressure: " << *stress_ptr;
        MOFEM_LOG_SEVERITY_SYNC(mField.get_comm(), Sev::inform);
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
      SmartPetscObj<IS> is_H;
      CHKERR is_mng->isCreateProblemFieldAndRank(name_prb, ROW, "H", 0, 0,
                                                 is_H);
      IS is_tmp;
      CHKERR ISExpand(is_H, is_flux, &is_tmp);
      auto is_Flux = SmartPetscObj<IS>(is_tmp);

      auto is_all_bc = bc_mng->getBlockIS(name_prb, "FLUIDFLUX", "FLUX", 0, 0);
      int is_all_bc_size;
      CHKERR ISGetSize(is_all_bc, &is_all_bc_size);
      MOFEM_LOG("ThermoElastic", Sev::inform)
          << "Field split block size " << is_all_bc_size;
      if (is_all_bc_size) {
        IS is_tmp2;
        CHKERR ISDifference(is_Flux, is_all_bc, &is_tmp2);
        is_Flux = SmartPetscObj<IS>(is_tmp2);
        CHKERR PCFieldSplitSetIS(pc, PETSC_NULL,
                                 is_all_bc); // boundary block
      }

      CHKERR ISSort(is_u);
      CHKERR ISSort(is_Flux);
      CHKERR PCFieldSplitSetIS(pc, PETSC_NULL, is_Flux);
      CHKERR PCFieldSplitSetIS(pc, PETSC_NULL, is_u);
    }

    MoFEMFunctionReturnHot(0);
  };

  auto pre_proc_ptr = boost::make_shared<FEMethod>();
  auto post_proc_rhs_ptr = boost::make_shared<FEMethod>();
  auto post_proc_lhs_ptr = boost::make_shared<FEMethod>();
  auto time_scale = boost::make_shared<TimeScale>();

  auto get_bc_hook_rhs = [this, pre_proc_ptr, time_scale]() {
    EssentialPreProc<DisplacementCubitBcData> hook(mField, pre_proc_ptr,
                                                   {time_scale}, false);
    return hook;
  };

  auto get_post_proc_hook_rhs = [this, post_proc_rhs_ptr]() {
    MoFEMFunctionBegin;
    CHKERR EssentialPreProcReaction<DisplacementCubitBcData>(
        mField, post_proc_rhs_ptr, nullptr, Sev::verbose)();
    CHKERR EssentialPostProcRhs<DisplacementCubitBcData>(
        mField, post_proc_rhs_ptr, 1.)();
    MoFEMFunctionReturn(0);
  };
  auto get_post_proc_hook_lhs = [this, post_proc_lhs_ptr]() {
    return EssentialPostProcLhs<DisplacementCubitBcData>(mField,
                                                         post_proc_lhs_ptr, 1.);
  };

  pre_proc_ptr->preProcessHook = get_bc_hook_rhs();
  post_proc_rhs_ptr->postProcessHook = get_post_proc_hook_rhs;
  post_proc_lhs_ptr->postProcessHook = get_post_proc_hook_lhs();

  auto ts_ctx_ptr = getDMTsCtx(dm);
  ts_ctx_ptr->getPreProcessIFunction().push_front(pre_proc_ptr);
  ts_ctx_ptr->getPreProcessIJacobian().push_front(pre_proc_ptr);
  ts_ctx_ptr->getPostProcessIFunction().push_back(post_proc_rhs_ptr);
  ts_ctx_ptr->getPostProcessIJacobian().push_back(post_proc_lhs_ptr);

  auto D = createDMVector(dm);
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
      LogManager::createSink(LogManager::getStrmWorld(), "Seepage"));
  LogManager::setLog("Seepage");
  MOFEM_LOG_TAG("Seepage", "seepage");
  core_log->add_sink(
      LogManager::createSink(LogManager::getStrmSync(), "SeepageSync"));
  LogManager::setLog("SeepageSync");
  MOFEM_LOG_TAG("SeepageSync", "seepage");

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
    MoFEM::Interface &m_field = core; ///< finite element database insterface
    //! [Create MoFEM]

    //! [Load mesh]
    Simple *simple = m_field.getInterface<Simple>();
    CHKERR simple->getOptions();
    CHKERR simple->loadFile();
    //! [Load mesh]

    //! [Seepage]
    Seepage ex(m_field);
    CHKERR ex.runProblem();
    //! [Seepage]
  }
  CATCH_ERRORS;

  CHKERR MoFEM::Core::Finalize();
}