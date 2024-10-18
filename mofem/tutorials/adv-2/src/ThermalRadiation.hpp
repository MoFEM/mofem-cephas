/**
 * @file ThermalRadiation.hpp
 * @author Andrei Shvarts (andrei.shvarts@glasgow.ac.uk)
 * @brief Implementation of thermal radiation bc
 * @version 0.14.0
 * @date 2024-08-23
 *
 * @copyright Copyright (c) 2024
 *
 */

namespace ThermoElasticOps {

template <CubitBC BC> struct RadiationBcType {};

template <CubitBC BC> struct GetRadiationParameters;

template <> struct GetRadiationParameters<BLOCKSET> {
  GetRadiationParameters() = delete;

  static MoFEMErrorCode
  getParameters(double &stefan_boltzmann_constant, double &emissivity_parameter,
                double &ambient_temperature, boost::shared_ptr<Range> &ents,
                MoFEM::Interface &m_field, int ms_id, Sev sev) {

    MoFEMFunctionBegin;

    auto cubit_meshset_ptr =
        m_field.getInterface<MeshsetsManager>()->getCubitMeshsetPtr(ms_id,
                                                                    BLOCKSET);
    std::vector<double> block_data;
    CHKERR cubit_meshset_ptr->getAttributes(block_data);
    if (block_data.size() < 3) {
      SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
              "Expected that radiation block has three attributes "
              "(Stefan–Boltzmann "
              "constant, emissivity parameter, and ambient temperature)");
    }
    stefan_boltzmann_constant = block_data[0];
    emissivity_parameter = block_data[1];
    ambient_temperature = block_data[2];

    MOFEM_LOG_CHANNEL("WORLD");
    MOFEM_TAG_AND_LOG("WORLD", Sev::inform, "RadiationBc")
        << "Stefan–Boltzmann constant " << stefan_boltzmann_constant;
    MOFEM_TAG_AND_LOG("WORLD", Sev::inform, "RadiationBc")
        << "Emissivity parameter " << emissivity_parameter;
    MOFEM_TAG_AND_LOG("WORLD", Sev::inform, "RadiationBc")
        << "Ambient temperature " << ambient_temperature;

    ents = boost::make_shared<Range>();
    CHKERR m_field.get_moab().get_entities_by_handle(cubit_meshset_ptr->meshset,
                                                     *(ents), true);

    MOFEM_LOG_CHANNEL("SYNC");
    MOFEM_TAG_AND_LOG("SYNC", Sev::noisy, "RadiationBc") << *ents;
    MOFEM_LOG_SEVERITY_SYNC(m_field.get_comm(), Sev::noisy);

    MoFEMFunctionReturn(0);
  }
};

} // namespace ThermoElasticOps

template <AssemblyType A, typename EleOp>
struct OpFluxRhsImpl<ThermoElasticOps::RadiationBcType<BLOCKSET>, 1, 1, A,
                     GAUSS, EleOp> : OpBaseImpl<A, EleOp> {

  using OpBase = OpBaseImpl<A, EleOp>;

  OpFluxRhsImpl(MoFEM::Interface &m_field, int ms_id, std::string field_name,
                boost::shared_ptr<VectorDouble> t_ptr);

protected:
  double stefanBoltzmanConstant;
  double emmisivityParameter;
  double ambientTemperature;
  boost::shared_ptr<VectorDouble> tPtr;
  MoFEMErrorCode iNtegrate(EntitiesFieldData::EntData &data);
};

template <AssemblyType A, typename EleOp>
struct OpFluxLhsImpl<ThermoElasticOps::RadiationBcType<BLOCKSET>, 1, 1, A,
                     GAUSS, EleOp> : OpBaseImpl<A, EleOp> {

  using OpBase = OpBaseImpl<A, EleOp>;

  OpFluxLhsImpl(MoFEM::Interface &m_field, int ms_id,
                std::string row_field_name, std::string col_field_name,
                boost::shared_ptr<VectorDouble> t_ptr);

protected:
  double stefanBoltzmanConstant;
  double emmisivityParameter;
  double ambientTemperature;
  boost::shared_ptr<VectorDouble>
      tPtr; // pointer for boundary temperature field
  MoFEMErrorCode iNtegrate(EntitiesFieldData::EntData &row_data,
                           EntitiesFieldData::EntData &col_data);
};

template <AssemblyType A, typename EleOp>
OpFluxRhsImpl<ThermoElasticOps::RadiationBcType<BLOCKSET>, 1, 1, A, GAUSS,
              EleOp>::OpFluxRhsImpl(MoFEM::Interface &m_field, int ms_id,
                                    std::string field_name,
                                    boost::shared_ptr<VectorDouble> t_ptr)
    : OpBase(field_name, field_name, OpBase::OPROW), tPtr(t_ptr) {
  CHK_THROW_MESSAGE(
      ThermoElasticOps::GetRadiationParameters<BLOCKSET>::getParameters(
          this->stefanBoltzmanConstant, this->emmisivityParameter,
          this->ambientTemperature, this->entsPtr, m_field, ms_id, Sev::inform),
      "Cannot read radiation data from blockset");
}

template <AssemblyType A, typename EleOp>
MoFEMErrorCode
MoFEM::OpFluxRhsImpl<ThermoElasticOps::RadiationBcType<BLOCKSET>, 1, 1, A,
                     GAUSS, EleOp>::iNtegrate(EntitiesFieldData::EntData
                                                  &row_data) {
  MoFEMFunctionBegin;

  // get integration weights
  auto t_w = OpBase::getFTensor0IntegrationWeight();
  // get base function values on rows
  auto t_row_base = row_data.getFTensor0N();

  // get temperature
  auto t_temp = getFTensor0FromVec(*tPtr);
  // loop over integration points
  for (int gg = 0; gg != OpBase::nbIntegrationPts; ++gg) {
    const auto alpha =
        -t_w * (std::pow(t_temp, 4) - std::pow(ambientTemperature, 4));

    int rr = 0;
    for (; rr != OpBase::nbRows; ++rr) {
      OpBase::locF[rr] += (alpha * t_row_base);
      ++t_row_base;
    }

    for (; rr < OpBase::nbRowBaseFunctions; ++rr)
      ++t_row_base;

    ++t_w;
    ++t_temp;
  }

  // get normal vector
  auto t_normal = OpBase::getFTensor1Normal();
  OpBase::locF *= t_normal.l2() * stefanBoltzmanConstant * emmisivityParameter;

  EntityType fe_type = OpBase::getNumeredEntFiniteElementPtr()->getEntType();
  if (fe_type == MBTRI) {
    OpBase::locF /= 2;
  }

  MoFEMFunctionReturn(0);
}

template <AssemblyType A, typename EleOp>
OpFluxLhsImpl<ThermoElasticOps::RadiationBcType<BLOCKSET>, 1, 1, A, GAUSS,
              EleOp>::OpFluxLhsImpl(MoFEM::Interface &m_field, int ms_id,
                                    std::string row_field_name,
                                    std::string col_field_name,
                                    boost::shared_ptr<VectorDouble> t_ptr)
    : OpBase(row_field_name, col_field_name, OpBase::OPROWCOL), tPtr(t_ptr) {

  this->sYmm = true;

  CHK_THROW_MESSAGE(
      ThermoElasticOps::GetRadiationParameters<BLOCKSET>::getParameters(
          this->stefanBoltzmanConstant, this->emmisivityParameter,
          this->ambientTemperature, this->entsPtr, m_field, ms_id,
          Sev::verbose),
      "Can not read radiation bc from blockset");
}

template <AssemblyType A, typename EleOp>
MoFEMErrorCode
MoFEM::OpFluxLhsImpl<ThermoElasticOps::RadiationBcType<BLOCKSET>, 1, 1, A,
                     GAUSS,
                     EleOp>::iNtegrate(EntitiesFieldData::EntData &row_data,
                                       EntitiesFieldData::EntData &col_data) {

  MoFEMFunctionBegin;
  // get element volume
  const double vol = OpBase::getMeasure();
  // get integration weights
  auto t_w = OpBase::getFTensor0IntegrationWeight();
  // get base function values on rows
  auto t_row_base = row_data.getFTensor0N();
  // get temperature
  auto t_temp = getFTensor0FromVec(*tPtr);
  // loop over integration points
  for (int gg = 0; gg != OpBase::nbIntegrationPts; ++gg) {
    const auto alpha = t_w * 4 * std::pow(t_temp, 3);
    ++t_w; // move to another integration weight
    ++t_temp;
    // loop over rows base functions
    int rr = 0;
    for (; rr != OpBase::nbRows; rr++) {
      const auto beta = alpha * t_row_base;
      // get column base functions values at gauss point gg
      auto t_col_base = col_data.getFTensor0N(gg, 0);
      // loop over columns
      for (int cc = 0; cc != OpBase::nbCols; cc++) {
        OpBase::locMat(rr, cc) += beta * t_col_base;
        ++t_col_base;
      }
      ++t_row_base;
    }
    for (; rr < OpBase::nbRowBaseFunctions; ++rr)
      ++t_row_base;
  }
  // get normal vector
  auto t_normal = OpBase::getFTensor1Normal();
  // using L2 norm of normal vector for convenient area calculation for quads,
  // tris etc.
  OpBase::locMat *=
      -t_normal.l2() * stefanBoltzmanConstant * emmisivityParameter;

  EntityType fe_type = OpBase::getNumeredEntFiniteElementPtr()->getEntType();
  if (fe_type == MBTRI) {
    OpBase::locMat /= 2;
  }

  MoFEMFunctionReturn(0);
}

template <CubitBC BCTYPE, int BASE_DIM, int FIELD_DIM, AssemblyType A,
          IntegrationType I, typename OpBase>
struct AddFluxToRhsPipelineImpl<

    OpFluxRhsImpl<ThermoElasticOps::RadiationBcType<BCTYPE>, BASE_DIM,
                  FIELD_DIM, A, I, OpBase>,
    A, I, OpBase

    > {

  AddFluxToRhsPipelineImpl() = delete;

  static MoFEMErrorCode add(

      boost::ptr_deque<ForcesAndSourcesCore::UserDataOperator> &pipeline,
      MoFEM::Interface &m_field, std::string field_name,
      boost::shared_ptr<VectorDouble> T_ptr, std::string block_name, Sev sev

  ) {
    MoFEMFunctionBegin;

    using OP = OpFluxRhsImpl<ThermoElasticOps::RadiationBcType<BLOCKSET>,
                             BASE_DIM, FIELD_DIM, PETSC, GAUSS, OpBase>;

    auto add_op = [&](auto &&meshset_vec_ptr) {
      for (auto m : meshset_vec_ptr) {
        MOFEM_TAG_AND_LOG("WORLD", sev, "OpRadiationRhs") << "Add " << *m;
        pipeline.push_back(
            new OP(m_field, m->getMeshsetId(), field_name, T_ptr));
      }
      MOFEM_LOG_CHANNEL("WORLD");
    };

    switch (BCTYPE) {
    case BLOCKSET:
      add_op(

          m_field.getInterface<MeshsetsManager>()->getCubitMeshsetPtr(
              std::regex(

                  (boost::format("%s(.*)") % block_name).str()

                      ))

      );

      break;
    default:
      SETERRQ(PETSC_COMM_SELF, MOFEM_NOT_IMPLEMENTED,
              "Handling of bc type not implemented");
      break;
    }
    MoFEMFunctionReturn(0);
  }
};

template <CubitBC BCTYPE, int BASE_DIM, int FIELD_DIM, AssemblyType A,
          IntegrationType I, typename OpBase>
struct AddFluxToLhsPipelineImpl<

    OpFluxLhsImpl<ThermoElasticOps::RadiationBcType<BCTYPE>, BASE_DIM,
                  FIELD_DIM, A, I, OpBase>,
    A, I, OpBase

    > {

  AddFluxToLhsPipelineImpl() = delete;

  static MoFEMErrorCode add(

      boost::ptr_deque<ForcesAndSourcesCore::UserDataOperator> &pipeline,
      MoFEM::Interface &m_field, std::string row_field_name,
      std::string col_field_name, boost::shared_ptr<VectorDouble> T_ptr,
      std::string block_name, Sev sev

  ) {
    MoFEMFunctionBegin;

    using OP = OpFluxLhsImpl<ThermoElasticOps::RadiationBcType<BLOCKSET>,
                             BASE_DIM, FIELD_DIM, PETSC, GAUSS, OpBase>;

    auto add_op = [&](auto &&meshset_vec_ptr) {
      for (auto m : meshset_vec_ptr) {
        MOFEM_TAG_AND_LOG("WORLD", sev, "OpRadiationLhs") << "Add " << *m;
        pipeline.push_back(new OP(m_field, m->getMeshsetId(), row_field_name,
                                  col_field_name, T_ptr));
      }
      MOFEM_LOG_CHANNEL("WORLD");
    };

    switch (BCTYPE) {
    case BLOCKSET:
      add_op(

          m_field.getInterface<MeshsetsManager>()->getCubitMeshsetPtr(
              std::regex(

                  (boost::format("%s(.*)") % block_name).str()

                      ))

      );

      break;
    default:
      SETERRQ(PETSC_COMM_SELF, MOFEM_NOT_IMPLEMENTED,
              "Handling of bc type not implemented");
      break;
    }
    MoFEMFunctionReturn(0);
  }
};
