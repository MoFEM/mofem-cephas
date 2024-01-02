/**
 * @file FluidLevel.hpp
 * @author your name (you@domain.com)
 * @brief
 * @version 0.1
 * @date 2024-01-02
 *
 * @copyright Copyright (c) 2024
 *
 */

namespace ElasticExample {

template <CubitBC BC> struct FluidLevelType {};

template <CubitBC BC> struct GetFluidLevel;

template <> struct GetFluidLevel<BLOCKSET> {
  GetFluidLevel() = delete;

  static MoFEMErrorCode getParams(double &density,
                                  std::array<double, 3> &gravity,
                                  std::array<double, 3> &zero_pressure,
                                  boost::shared_ptr<Range> &ents,
                                  MoFEM::Interface &m_field, int ms_id) {
    MoFEMFunctionBegin;

    auto cubit_meshset_ptr =
        m_field.getInterface<MeshsetsManager>()->getCubitMeshsetPtr(ms_id,
                                                                    BLOCKSET);
    std::vector<double> block_data;
    CHKERR cubit_meshset_ptr->getAttributes(block_data);
    if (block_data.size() != 1 + 3 + 3) {
      SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
              "Expected that block has two attribute");
    }
    density = block_data[0];
    gravity = {block_data[1], block_data[2], block_data[3]};
    zero_pressure = {block_data[4], block_data[5], block_data[6]};

    MOFEM_LOG_CHANNEL("WORLD");
    MOFEM_TAG_AND_LOG("WORLD", Sev::inform, "FluidPressure")
        << "Density " << density;
    MOFEM_TAG_AND_LOG("WORLD", Sev::inform, "FluidPressure")
        << "Gravity " << gravity[0] << " " << gravity[1] << " " << gravity[2];
    MOFEM_TAG_AND_LOG("WORLD", Sev::inform, "FluidPressure")
        << "Zero pressure " << zero_pressure[0] << " " << zero_pressure[1]
        << " " << zero_pressure[2];

    ents = boost::make_shared<Range>();
    CHKERR
    m_field.get_moab().get_entities_by_handle(cubit_meshset_ptr->meshset,
                                              *(ents), true);

    MOFEM_LOG_CHANNEL("SYNC");
    MOFEM_TAG_AND_LOG("SYNC", Sev::noisy, "FluidPressure") << *ents;
    MOFEM_LOG_SEVERITY_SYNC(m_field.get_comm(), Sev::noisy);

    MoFEMFunctionReturn(0);
  }
};

} // namespace ElasticExample

template <int FIELD_DIM, AssemblyType A, typename EleOp>
struct OpFluxRhsImpl<ElasticExample::FluidLevelType<BLOCKSET>, 1, FIELD_DIM, A,
                     GAUSS, EleOp> : OpBaseImpl<A, EleOp> {

  using OpBase = OpBaseImpl<A, EleOp>;

  OpFluxRhsImpl(MoFEM::Interface &m_field, int ms_id, std::string field_name,
                double scale);

protected:
  double fluidDensity;
  std::array<double, 3> gravityAcceleration;
  std::array<double, 3> zeroPressure;

  double rhsScale;
  boost::shared_ptr<MatrixDouble> uPtr;
  MoFEMErrorCode iNtegrate(EntitiesFieldData::EntData &data);
};

template <int FIELD_DIM, AssemblyType A, typename EleOp>
OpFluxRhsImpl<ElasticExample::FluidLevelType<BLOCKSET>, 1, FIELD_DIM, A, GAUSS,
              EleOp>::OpFluxRhsImpl(MoFEM::Interface &m_field, int ms_id,
                                    std::string field_name, double scale)
    : OpBase(field_name, field_name, OpBase::OPROW), rhsScale(scale) {
  CHK_THROW_MESSAGE(ElasticExample::GetFluidLevel<BLOCKSET>::getParams(
                        this->fluidDensity, this->gravityAcceleration,
                        this->zeroPressure, this->entsPtr, m_field, ms_id),
                    "Can not read spring stiffness data from blockset");
}

template <int FIELD_DIM, AssemblyType A, typename EleOp>
MoFEMErrorCode
MoFEM::OpFluxRhsImpl<ElasticExample::FluidLevelType<BLOCKSET>, 1, FIELD_DIM, A,
                     GAUSS, EleOp>::iNtegrate(EntitiesFieldData::EntData
                                                  &row_data) {
  FTensor::Index<'i', FIELD_DIM> i;
  FTensor::Index<'j', FIELD_DIM> j;
  MoFEMFunctionBegin;
  // get element volume
  const double vol = OpBase::getMeasure();
  // get integration weights
  auto t_w = OpBase::getFTensor0IntegrationWeight();
  // get base function gradient on rows
  auto t_row_base = row_data.getFTensor0N();

  // get coordinate at integration points
  auto t_normal = OpBase::getFTensor1NormalsAtGaussPts();
	// get coordinate at integration points
  auto t_coord = OpBase::getFTensor1CoordsAtGaussPts();

  auto t_gravity = getFTensor1FromPtr<FIELD_DIM>(gravityAcceleration.data());
  auto t_zero_pressure = getFTensor1FromPtr<FIELD_DIM>(zeroPressure.data());

  // loop over integration points
  for (int gg = 0; gg != OpBase::nbIntegrationPts; ++gg) {
    // take into account Jacobian
    const double alpha = t_w * rhsScale;

    FTensor::Tensor1<double, FIELD_DIM> t_level;
    t_level(i) = t_coord(i) - t_zero_pressure(i);
    auto dot = t_gravity(i) * t_level(i);

    FTensor::Tensor1<double, FIELD_DIM> t_force;
    t_force(i) = -fluidDensity * dot * t_normal(i);

    // loop over rows base functions
    auto t_nf = getFTensor1FromArray<FIELD_DIM, FIELD_DIM>(OpBase::locF);
    int rr = 0;
    for (; rr != OpBase::nbRows / FIELD_DIM; ++rr) {
      t_nf(i) += (alpha * t_row_base) * t_force(i);
      ++t_row_base;
      ++t_nf;
    }

    for (; rr < OpBase::nbRowBaseFunctions; ++rr)
      ++t_row_base;

    ++t_w;
    ++t_coord;
    ++t_normal;
  }
  MoFEMFunctionReturn(0);
}

template <CubitBC BCTYPE, int BASE_DIM, int FIELD_DIM, AssemblyType A,
          IntegrationType I, typename OpBase>
struct AddFluxToRhsPipelineImpl<

    OpFluxRhsImpl<ElasticExample::FluidLevelType<BCTYPE>, BASE_DIM, FIELD_DIM,
                  A, I, OpBase>,
    A, I, OpBase

    > {

  AddFluxToRhsPipelineImpl() = delete;

  static MoFEMErrorCode add(

      boost::ptr_deque<ForcesAndSourcesCore::UserDataOperator> &pipeline,
      MoFEM::Interface &m_field, std::string field_name, double scale,
      std::string block_name, Sev sev

  ) {
    MoFEMFunctionBegin;

    using OP = OpFluxRhsImpl<ElasticExample::FluidLevelType<BLOCKSET>, BASE_DIM,
                             FIELD_DIM, PETSC, GAUSS, OpBase>;

    auto add_op = [&](auto &&meshset_vec_ptr) {
      for (auto m : meshset_vec_ptr) {
        MOFEM_TAG_AND_LOG("WORLD", sev, "OFluidLevelRhs") << "Add " << *m;
        pipeline.push_back(
            new OP(m_field, m->getMeshsetId(), field_name, scale));
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
