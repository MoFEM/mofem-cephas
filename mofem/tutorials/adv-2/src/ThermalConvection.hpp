/**
 * @file ThermalConvection.hpp
 * @author Andrei Shvarts (andrei.shvarts@glasgow.ac.uk)
 * @brief Implementation of thermal convection bc
 * @version 0.14.0
 * @date 2024-08-23
 */

namespace ThermoElasticOps {

template <CubitBC BC> struct ConvectionBcType {};

template <CubitBC BC> struct GetConvectionParameters;

template <> struct GetConvectionParameters<BLOCKSET> {
  GetConvectionParameters() = delete;

  static MoFEMErrorCode getParameters(double &heat_transfer_coefficient,
                                      double &ambient_temperature,
                                      boost::shared_ptr<Range> &ents,
                                      MoFEM::Interface &m_field, int ms_id,
                                      Sev sev) {

    MoFEMFunctionBegin;

    auto cubit_meshset_ptr =
        m_field.getInterface<MeshsetsManager>()->getCubitMeshsetPtr(ms_id,
                                                                    BLOCKSET);
    std::vector<double> block_data;
    CHKERR cubit_meshset_ptr->getAttributes(block_data);
    if (block_data.size() < 2) {
      SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
              "Expected that convection block has two attributes (heat "
              "transfer coefficient and ambient temperature)");
    }
    heat_transfer_coefficient = block_data[0];
    ambient_temperature = block_data[1];

    MOFEM_LOG_CHANNEL("WORLD");
    MOFEM_TAG_AND_LOG("WORLD", Sev::inform, "ConvectionBc")
        << "Heat transfer coefficient " << heat_transfer_coefficient;
    MOFEM_TAG_AND_LOG("WORLD", Sev::inform, "ConvectionBc")
        << "Ambient temperature " << ambient_temperature;

    ents = boost::make_shared<Range>();
    CHKERR
    m_field.get_moab().get_entities_by_handle(cubit_meshset_ptr->meshset,
                                              *(ents), true);

    MOFEM_LOG_CHANNEL("SYNC");
    MOFEM_TAG_AND_LOG("SYNC", Sev::noisy, "ConvectionBc") << *ents;
    MOFEM_LOG_SEVERITY_SYNC(m_field.get_comm(), Sev::noisy);

    MoFEMFunctionReturn(0);
  }
};

} // namespace ThermoElasticOps

template <AssemblyType A, typename EleOp>
struct OpFluxRhsImpl<ThermoElasticOps::ConvectionBcType<BLOCKSET>, 1, 1, A,
                     GAUSS, EleOp> : OpBaseImpl<A, EleOp> {

  using OpBase = OpBaseImpl<A, EleOp>;

  OpFluxRhsImpl(MoFEM::Interface &m_field, int ms_id, std::string field_name,
                boost::shared_ptr<VectorDouble> t_ptr);

protected:
  double heatTransferCoefficient;
  double ambientTemperature;
  boost::shared_ptr<VectorDouble> tPtr;
  MoFEMErrorCode iNtegrate(EntitiesFieldData::EntData &data);
};

template <AssemblyType A, typename EleOp>
struct OpFluxLhsImpl<ThermoElasticOps::ConvectionBcType<BLOCKSET>, 1, 1, A,
                     GAUSS, EleOp> : OpBaseImpl<A, EleOp> {

  using OpBase = OpBaseImpl<A, EleOp>;

  OpFluxLhsImpl(MoFEM::Interface &m_field, int ms_id,
                std::string row_field_name, std::string col_field_name);

protected:
  double heatTransferCoefficient;
  double ambientTemperature;
  MoFEMErrorCode iNtegrate(EntitiesFieldData::EntData &row_data,
                           EntitiesFieldData::EntData &col_data);
};

template <AssemblyType A, typename EleOp>
OpFluxRhsImpl<ThermoElasticOps::ConvectionBcType<BLOCKSET>, 1, 1, A, GAUSS,
              EleOp>::OpFluxRhsImpl(MoFEM::Interface &m_field, int ms_id,
                                    std::string field_name,
                                    boost::shared_ptr<VectorDouble> t_ptr)
    : OpBase(field_name, field_name, OpBase::OPROW), tPtr(t_ptr) {
  CHK_THROW_MESSAGE(
      ThermoElasticOps::GetConvectionParameters<BLOCKSET>::getParameters(
          this->heatTransferCoefficient, this->ambientTemperature,
          this->entsPtr, m_field, ms_id, Sev::inform),
      "Cannot read convection data from blockset");
}

template <AssemblyType A, typename EleOp>
MoFEMErrorCode
MoFEM::OpFluxRhsImpl<ThermoElasticOps::ConvectionBcType<BLOCKSET>, 1, 1, A,
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
        -t_w * heatTransferCoefficient * (t_temp - ambientTemperature);
    int rr = 0;
    for (; rr != OpBase::nbRows; ++rr) {
      OpBase::locF[rr] += alpha * t_row_base;
      ++t_row_base;
    }

    for (; rr < OpBase::nbRowBaseFunctions; ++rr)
      ++t_row_base;

    ++t_w;
    ++t_temp;
  }

  // get normal vector
  auto t_normal = OpBase::getFTensor1Normal();
  OpBase::locF *= t_normal.l2();

  EntityType fe_type = OpBase::getNumeredEntFiniteElementPtr()->getEntType();
  if (fe_type == MBTRI) {
    OpBase::locF /= 2;
  }

  MoFEMFunctionReturn(0);
}

template <AssemblyType A, typename EleOp>
OpFluxLhsImpl<ThermoElasticOps::ConvectionBcType<BLOCKSET>, 1, 1, A, GAUSS,
              EleOp>::OpFluxLhsImpl(MoFEM::Interface &m_field, int ms_id,
                                    std::string row_field_name,
                                    std::string col_field_name)
    : OpBase(row_field_name, col_field_name, OpBase::OPROWCOL) {

  this->sYmm = true;

  CHK_THROW_MESSAGE(
      ThermoElasticOps::GetConvectionParameters<BLOCKSET>::getParameters(
          this->heatTransferCoefficient, this->ambientTemperature,
          this->entsPtr, m_field, ms_id, Sev::verbose),
      "Cannot read convection data from blockset");
}

template <AssemblyType A, typename EleOp>
MoFEMErrorCode
MoFEM::OpFluxLhsImpl<ThermoElasticOps::ConvectionBcType<BLOCKSET>, 1, 1, A,
                     GAUSS,
                     EleOp>::iNtegrate(EntitiesFieldData::EntData &row_data,
                                       EntitiesFieldData::EntData &col_data) {

  MoFEMFunctionBegin;
  // get integration weights
  auto t_w = OpBase::getFTensor0IntegrationWeight();
  // get base function values on rows
  auto t_row_base = row_data.getFTensor0N();
  // loop over integration points
  for (int gg = 0; gg != OpBase::nbIntegrationPts; ++gg) {
    // loop over rows base functions
    int rr = 0;
    for (; rr != OpBase::nbRows; rr++) {
      const auto alpha = t_w * t_row_base;
      // get column base functions values at gauss point gg
      auto t_col_base = col_data.getFTensor0N(gg, 0);
      // loop over columns
      for (int cc = 0; cc != OpBase::nbCols; cc++) {
        OpBase::locMat(rr, cc) += alpha * t_col_base;
        ++t_col_base;
      }
      ++t_row_base;
    }
    for (; rr < OpBase::nbRowBaseFunctions; ++rr)
      ++t_row_base;
    ++t_w; // move to another integration weight
  }
  // get normal vector
  auto t_normal = OpBase::getFTensor1Normal();
  // using L2 norm of normal vector for convenient area calculation for quads,
  // tris etc.
  OpBase::locMat *= -t_normal.l2() * heatTransferCoefficient;

  EntityType fe_type = OpBase::getNumeredEntFiniteElementPtr()->getEntType();
  if (fe_type == MBTRI) {
    OpBase::locMat /= 2;
  }

  MoFEMFunctionReturn(0);
}

template <CubitBC BCTYPE, int BASE_DIM, int FIELD_DIM, AssemblyType A,
          IntegrationType I, typename OpBase>
struct AddFluxToRhsPipelineImpl<

    OpFluxRhsImpl<ThermoElasticOps::ConvectionBcType<BCTYPE>, BASE_DIM,
                  FIELD_DIM, A, I, OpBase>,
    A, I, OpBase

    > {

  AddFluxToRhsPipelineImpl() = delete;

  static MoFEMErrorCode add(

      boost::ptr_deque<ForcesAndSourcesCore::UserDataOperator> &pipeline,
      MoFEM::Interface &m_field, std::string field_name,
      boost::shared_ptr<VectorDouble> u_ptr, std::string block_name, Sev sev

  ) {
    MoFEMFunctionBegin;

    using OP = OpFluxRhsImpl<ThermoElasticOps::ConvectionBcType<BLOCKSET>,
                             BASE_DIM, FIELD_DIM, PETSC, GAUSS, OpBase>;

    auto add_op = [&](auto &&meshset_vec_ptr) {
      for (auto m : meshset_vec_ptr) {
        MOFEM_TAG_AND_LOG("WORLD", sev, "OpConvectionRhs") << "Add " << *m;
        pipeline.push_back(
            new OP(m_field, m->getMeshsetId(), field_name, u_ptr));
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

    OpFluxLhsImpl<ThermoElasticOps::ConvectionBcType<BCTYPE>, BASE_DIM,
                  FIELD_DIM, A, I, OpBase>,
    A, I, OpBase

    > {

  AddFluxToLhsPipelineImpl() = delete;

  static MoFEMErrorCode add(

      boost::ptr_deque<ForcesAndSourcesCore::UserDataOperator> &pipeline,
      MoFEM::Interface &m_field, std::string row_field_name,
      std::string col_field_name, std::string block_name, Sev sev

  ) {
    MoFEMFunctionBegin;

    using OP = OpFluxLhsImpl<ThermoElasticOps::ConvectionBcType<BLOCKSET>,
                             BASE_DIM, FIELD_DIM, PETSC, GAUSS, OpBase>;

    auto add_op = [&](auto &&meshset_vec_ptr) {
      for (auto m : meshset_vec_ptr) {
        MOFEM_TAG_AND_LOG("WORLD", sev, "OpConvectionLhs") << "Add " << *m;
        pipeline.push_back(
            new OP(m_field, m->getMeshsetId(), row_field_name, col_field_name));
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
