/**
 * @file Natural.hpp
 * @brief Setting natural boundary conditions
 * @version 0.13.1
 * @date 2022-08-12
 *
 * @copyright Copyright (c) 2022a
 *
 */

#ifndef _NATURAL_BC_
#define _NATURAL_BC_

namespace MoFEM {

using VecOfTimeScalingMethods = std::vector<boost::shared_ptr<ScalingMethod>>;

template <CubitBC BCTYP> struct NaturalMeshsetType {};
struct NaturalForceMeshsets {};
struct NaturalTemperatureMeshsets {};

template <typename OP> struct FluxOpType {};

template <typename T, int BASE_DIM, int FIELD_DIM, AssemblyType A,
          IntegrationType I, typename OpBase>
struct OpFluxImpl;

template <typename T, AssemblyType A, IntegrationType I, typename OpBase>
struct AddFluxToPipelineImpl;

/**
 * @brief Natural boundary conditions
 * @ingroup mofem_forms
 *
 * @tparam EleOp
 */
template <typename EleOp> struct NaturalBC {

  using EntData = EntitiesFieldData::EntData;
  using OpType = typename EleOp::OpType;

  /**
   * @brief Assembly methods
   * @ingroup mofem_forms
   *
   * @tparam A
   */
  template <AssemblyType A> struct Assembly {

    template <IntegrationType I> struct LinearForm {

      template <typename T, int BASE_DIM, int FIELD_DIM>
      using OpFlux = OpFluxImpl<T, BASE_DIM, FIELD_DIM, A, I, EleOp>;

      template <typename T>
      using AddFluxToPipeline = AddFluxToPipelineImpl<T, A, I, EleOp>;

      template <typename T>
      static MoFEMErrorCode addFluxToPipeline(
          FluxOpType<T>,
          boost::ptr_vector<ForcesAndSourcesCore::UserDataOperator> &pipeline,
          MoFEM::Interface &m_field, const std::string field_name,
          std::vector<boost::shared_ptr<ScalingMethod>> smv,
          const std::string block_name = "", Sev sev = Sev::noisy);
    };

  }; // Assembly
};

// Declarations

template <int FIELD_DIM, AssemblyType A, IntegrationType I, typename OpBase>
struct OpFluxImpl<NaturalMeshsetType<UNKNOWNSET>, 1, FIELD_DIM, A, I, OpBase>
    : FormsIntegrators<OpBase>::template Assembly<A>::template LinearForm<
          I>::template OpSource<1, FIELD_DIM> {

  using OpSource = typename FormsIntegrators<OpBase>::template Assembly<
      A>::template LinearForm<I>::template OpSource<1, FIELD_DIM>;

  OpFluxImpl(const std::string field_name,
             FTensor::Tensor1<double, FIELD_DIM> t_force,
             boost::shared_ptr<Range> ents_ptr = nullptr);

  VecOfTimeScalingMethods &getVecOfTimeScalingMethods() {
    return vecOfTimeScalingMethods;
  }

protected:
  OpFluxImpl(const std::string field_name);
  FTensor::Tensor1<double, FIELD_DIM> tForce;
  VecOfTimeScalingMethods vecOfTimeScalingMethods;
};

template <int FIELD_DIM, AssemblyType A, IntegrationType I, typename OpBase>
struct OpFluxImpl<NaturalMeshsetType<FORCESET>, 1, FIELD_DIM, A, I, OpBase>
    : OpFluxImpl<NaturalMeshsetType<UNKNOWNSET>, 1, FIELD_DIM, A, I, OpBase> {

  OpFluxImpl(MoFEM::Interface &m_field, int ms_id,
             const std::string field_name);

protected:
  MoFEMErrorCode getMeshsetData(MoFEM::Interface &m_field, int ms_id);
};

template <int FIELD_DIM, AssemblyType A, IntegrationType I, typename OpBase>
struct OpFluxImpl<NaturalMeshsetType<BLOCKSET>, 1, FIELD_DIM, A, I, OpBase>
    : OpFluxImpl<NaturalMeshsetType<UNKNOWNSET>, 1, FIELD_DIM, A, I, OpBase> {

  OpFluxImpl(MoFEM::Interface &m_field, int ms_id,
             const std::string field_name);

protected:
  MoFEMErrorCode getMeshsetData(MoFEM::Interface &m_field, int ms_id);
};

template <int FIELD_DIM, AssemblyType A, IntegrationType I, typename OpBase>
struct OpFluxImpl<NaturalMeshsetType<PRESSURESET>, 1, FIELD_DIM, A, I, OpBase>
    : OpFluxImpl<NaturalMeshsetType<UNKNOWNSET>, 1, FIELD_DIM, A, I, OpBase> {

  using Parent =
      OpFluxImpl<NaturalMeshsetType<UNKNOWNSET>, 1, FIELD_DIM, A, I, OpBase>;
  using EntData = EntitiesFieldData::EntData;

  OpFluxImpl(MoFEM::Interface &m_field, int ms_id,
             const std::string field_name);

protected:
  MoFEMErrorCode iNtegrate(EntData &data);
  MoFEMErrorCode getMeshsetData(MoFEM::Interface &m_field, int ms_id);

  double surfacePressure;
};

template <int FIELD_DIM, AssemblyType A, IntegrationType I, typename OpBase>
struct OpFluxImpl<NaturalForceMeshsets, 1, FIELD_DIM, A, I, OpBase>;

template <int BASE_DIM, AssemblyType A, IntegrationType I, typename OpBase>
struct OpFluxImpl<NaturalTemperatureMeshsets, BASE_DIM, BASE_DIM, A, I, OpBase>;

template <AssemblyType A, IntegrationType I, typename OpBase>
struct OpFluxImpl<NaturalMeshsetType<UNKNOWNSET>, 1, 1, A, I, OpBase>
    : FormsIntegrators<OpBase>::template Assembly<A>::template LinearForm<
          I>::template OpSource<1, 1> {

  using OpSource = typename FormsIntegrators<OpBase>::template Assembly<
      A>::template LinearForm<I>::template OpSource<1, 1>;

  OpFluxImpl(const std::string field_name, double value,
             boost::shared_ptr<Range> ents_ptr = nullptr);

  VecOfTimeScalingMethods &getVecOfTimeScalingMethods() {
    return vecOfTimeScalingMethods;
  }

protected:
  OpFluxImpl(const std::string field_name);

  MoFEMErrorCode getMeshsetData(MoFEM::Interface &m_field, int ms_id);
  double scalarValue;

  VecOfTimeScalingMethods vecOfTimeScalingMethods;
};

template <AssemblyType A, IntegrationType I, typename OpBase>
struct OpFluxImpl<NaturalMeshsetType<BLOCKSET>, 1, 1, A, I, OpBase>
    : OpFluxImpl<NaturalMeshsetType<UNKNOWNSET>, 1, 1, A, I, OpBase> {

  using Parent = OpFluxImpl<NaturalMeshsetType<UNKNOWNSET>, 1, 1, A, I, OpBase>;

  OpFluxImpl(MoFEM::Interface &m_field, int ms_id,
             const std::string field_name);

protected:
  MoFEMErrorCode getMeshsetData(MoFEM::Interface &m_field, int ms_id);
};

template <AssemblyType A, IntegrationType I, typename OpBase>
struct OpFluxImpl<NaturalMeshsetType<HEATFLUXSET>, 1, 1, A, I, OpBase>
    : OpFluxImpl<NaturalMeshsetType<UNKNOWNSET>, 1, 1, A, I, OpBase> {

  using Parent = OpFluxImpl<NaturalMeshsetType<UNKNOWNSET>, 1, 1, A, I, OpBase>;

  OpFluxImpl(MoFEM::Interface &m_field, int ms_id,
             const std::string field_name);

protected:
  MoFEMErrorCode getMeshsetData(MoFEM::Interface &m_field, int ms_id);
};

template <int BASE_DIM, AssemblyType A, IntegrationType I, typename OpBase>
struct OpFluxImpl<NaturalMeshsetType<UNKNOWNSET>, BASE_DIM, BASE_DIM, A, I,
                  OpBase>
    : FormsIntegrators<OpBase>::template Assembly<A>::template LinearForm<
          I>::template OpNormalMixVecTimesScalar<BASE_DIM> {

  using OpSource = typename FormsIntegrators<OpBase>::template Assembly<
      A>::template LinearForm<I>::template OpNormalMixVecTimesScalar<BASE_DIM>;

  OpFluxImpl(const std::string field_name, const double value,
             boost::shared_ptr<Range> ents_ptr = nullptr);

  VecOfTimeScalingMethods &getVecOfTimeScalingMethods() {
    return vecOfTimeScalingMethods;
  }

protected:
  OpFluxImpl(const std::string field_name);
  double scalarValue;

  VecOfTimeScalingMethods vecOfTimeScalingMethods;
};

template <int BASE_DIM, AssemblyType A, IntegrationType I, typename OpBase>
struct OpFluxImpl<NaturalMeshsetType<BLOCKSET>, BASE_DIM, BASE_DIM, A, I,
                  OpBase> : OpFluxImpl<NaturalMeshsetType<UNKNOWNSET>, BASE_DIM,
                                       BASE_DIM, A, I, OpBase> {
  OpFluxImpl(MoFEM::Interface &m_field, int ms_id,
             const std::string field_name);

protected:
  MoFEMErrorCode getMeshsetData(MoFEM::Interface &m_field, int ms_id);
};

template <int BASE_DIM, AssemblyType A, IntegrationType I, typename OpBase>
struct OpFluxImpl<NaturalMeshsetType<TEMPERATURESET>, BASE_DIM, BASE_DIM, A, I,
                  OpBase> : OpFluxImpl<NaturalMeshsetType<UNKNOWNSET>, BASE_DIM,
                                       BASE_DIM, A, I, OpBase> {
  OpFluxImpl(MoFEM::Interface &m_field, int ms_id,
             const std::string field_name);

protected:
  MoFEMErrorCode getMeshsetData(MoFEM::Interface &m_field, int ms_id);
};

// Definitions

template <int FIELD_DIM, AssemblyType A, IntegrationType I, typename OpBase>
OpFluxImpl<NaturalMeshsetType<UNKNOWNSET>, 1, FIELD_DIM, A, I,
           OpBase>::OpFluxImpl(const std::string field_name,
                               FTensor::Tensor1<double, FIELD_DIM> t_force,
                               boost::shared_ptr<Range> ents_ptr)
    : OpFluxImpl(field_name) {
  FTensor::Index<'i', FIELD_DIM> i;
  this->tForce(i) = t_force(i);
  this->entsPtr = ents_ptr;
}

template <int FIELD_DIM, AssemblyType A, IntegrationType I, typename OpBase>
OpFluxImpl<NaturalMeshsetType<UNKNOWNSET>, 1, FIELD_DIM, A, I,
           OpBase>::OpFluxImpl(const std::string field_name)
    : OpSource(
          field_name,

          [this](const double, const double, const double) { return tForce; }

      ) {

  this->timeScalingFun = [this](const double t) {
    double s = 1;
    for (auto &o : vecOfTimeScalingMethods) {
      s *= o->getScale(t);
    }
    return s;
  };

  static_assert(FIELD_DIM > 1, "Is not implemented for scalar field");
}

template <int FIELD_DIM, AssemblyType A, IntegrationType I, typename OpBase>
OpFluxImpl<NaturalMeshsetType<FORCESET>, 1, FIELD_DIM, A, I,
           OpBase>::OpFluxImpl(MoFEM::Interface &m_field, int ms_id,
                               const std::string field_name)
    : OpFluxImpl<NaturalMeshsetType<UNKNOWNSET>, 1, FIELD_DIM, A, I, OpBase>(
          field_name) {
  CHK_THROW_MESSAGE(getMeshsetData(m_field, ms_id), "Get meshset data");
}

template <int FIELD_DIM, AssemblyType A, IntegrationType I, typename OpBase>
MoFEMErrorCode OpFluxImpl<NaturalMeshsetType<FORCESET>, 1, FIELD_DIM, A, I,
                          OpBase>::getMeshsetData(MoFEM::Interface &m_field,
                                                  int ms_id) {
  MoFEMFunctionBegin;

  auto cubit_meshset_ptr =
      m_field.getInterface<MeshsetsManager>()->getCubitMeshsetPtr(ms_id,
                                                                  NODESET);

  ForceCubitBcData bc_data;
  CHKERR cubit_meshset_ptr->getBcDataStructure(bc_data);

  this->tForce(0) = bc_data.data.value3;
  if constexpr (FIELD_DIM > 1)
    this->tForce(1) = bc_data.data.value4;
  if constexpr (FIELD_DIM > 2)
    this->tForce(2) = bc_data.data.value5;

  FTensor::Index<'i', FIELD_DIM> i;
  this->tForce(i) *= bc_data.data.value1;

  this->entsPtr = boost::make_shared<Range>();
  CHKERR m_field.get_moab().get_entities_by_handle(cubit_meshset_ptr->meshset,
                                                   *(this->entsPtr), true);

  MoFEMFunctionReturn(0);
}

template <int FIELD_DIM, AssemblyType A, IntegrationType I, typename OpBase>
OpFluxImpl<NaturalMeshsetType<BLOCKSET>, 1, FIELD_DIM, A, I,
           OpBase>::OpFluxImpl(MoFEM::Interface &m_field, int ms_id,
                               const std::string field_name)
    : OpFluxImpl<NaturalMeshsetType<UNKNOWNSET>, 1, FIELD_DIM, A, I, OpBase>(
          field_name) {
  CHK_THROW_MESSAGE(getMeshsetData(m_field, ms_id), "Get meshset data");
}

template <int FIELD_DIM, AssemblyType A, IntegrationType I, typename OpBase>
MoFEMErrorCode OpFluxImpl<NaturalMeshsetType<BLOCKSET>, 1, FIELD_DIM, A, I,
                          OpBase>::getMeshsetData(MoFEM::Interface &m_field,
                                                  int ms_id) {
  MoFEMFunctionBegin;

  auto cubit_meshset_ptr =
      m_field.getInterface<MeshsetsManager>()->getCubitMeshsetPtr(ms_id,
                                                                  BLOCKSET);

  std::vector<double> block_data;
  CHKERR cubit_meshset_ptr->getAttributes(block_data);
  if (block_data.size() != FIELD_DIM) {
    MOFEM_LOG("SELF", Sev::warning)
        << "BLOCKSET is expected to have " << FIELD_DIM
        << " attributes but has size " << block_data.size();
    if (block_data.size() < FIELD_DIM) {
      SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
              "Size of attribute in BLOCKSET is too small");
    }
  }

  MOFEM_LOG("WORLD", Sev::inform)
      << "Flux blockset " << cubit_meshset_ptr->getName();
  MOFEM_LOG("WORLD", Sev::inform)
      << "Number of attributes " << block_data.size();

  for (unsigned int ii = 0; ii != FIELD_DIM; ++ii) {
    this->tForce(ii) = block_data[ii];
  }

  this->entsPtr = boost::make_shared<Range>();
  CHKERR m_field.get_moab().get_entities_by_handle(cubit_meshset_ptr->meshset,
                                                   *(this->entsPtr), true);

  MoFEMFunctionReturn(0);
}

template <int FIELD_DIM, AssemblyType A, IntegrationType I, typename OpBase>
OpFluxImpl<NaturalMeshsetType<PRESSURESET>, 1, FIELD_DIM, A, I,
           OpBase>::OpFluxImpl(MoFEM::Interface &m_field, int ms_id,
                               const std::string field_name)
    : OpFluxImpl<NaturalMeshsetType<UNKNOWNSET>, 1, FIELD_DIM, A, I, OpBase>(
          field_name) {
  CHK_THROW_MESSAGE(getMeshsetData(m_field, ms_id), "Get meshset data");
}

template <int FIELD_DIM, AssemblyType A, IntegrationType I, typename OpBase>
MoFEMErrorCode OpFluxImpl<NaturalMeshsetType<PRESSURESET>, 1, FIELD_DIM, A, I,
                          OpBase>::iNtegrate(EntData &data) {
  MoFEMFunctionBegin;

  auto t_normal = this->getFTensor1Normal();
  FTensor::Index<'i', FIELD_DIM> i;

  this->sourceFun = [&](const double, const double, const double) {
    this->tForce(i) = t_normal(i);
    this->tForce.normalize();
    this->tForce(i) *= this->surfacePressure;
    return this->tForce;
  };

  CHKERR Parent::iNtegrate(data);
  MoFEMFunctionReturn(0);
}

template <int FIELD_DIM, AssemblyType A, IntegrationType I, typename OpBase>
MoFEMErrorCode OpFluxImpl<NaturalMeshsetType<PRESSURESET>, 1, FIELD_DIM, A, I,
                          OpBase>::getMeshsetData(MoFEM::Interface &m_field,
                                                  int ms_id) {
  MoFEMFunctionBegin;

  auto cubit_meshset_ptr =
      m_field.getInterface<MeshsetsManager>()->getCubitMeshsetPtr(ms_id,
                                                                  SIDESET);

  PressureCubitBcData pressure_data;
  CHKERR cubit_meshset_ptr->getBcDataStructure(pressure_data);
  this->surfacePressure = pressure_data.data.value1;

  this->entsPtr = boost::make_shared<Range>();
  CHKERR m_field.get_moab().get_entities_by_handle(cubit_meshset_ptr->meshset,
                                                   *(this->entsPtr), true);

  MoFEMFunctionReturn(0);
}

template <AssemblyType A, IntegrationType I, typename OpBase>
OpFluxImpl<NaturalMeshsetType<UNKNOWNSET>, 1, 1, A, I, OpBase>::OpFluxImpl(
    const std::string field_name, const double value,
    boost::shared_ptr<Range> ents_ptr)
    : OpFluxImpl(field_name) {
  this->scalarValue = value;
  this->entsPtr = ents_ptr;
}

template <AssemblyType A, IntegrationType I, typename OpBase>
OpFluxImpl<NaturalMeshsetType<UNKNOWNSET>, 1, 1, A, I, OpBase>::OpFluxImpl(
    const std::string field_name)
    : OpSource(field_name,

               [this](const double, const double, const double) {
                 return scalarValue;
               }

      ) {

  this->timeScalingFun = [this](const double t) {
    double s = 1;
    for (auto &o : vecOfTimeScalingMethods) {
      s *= o->getScale(t);
    }
    return s;
  };
}

template <AssemblyType A, IntegrationType I, typename OpBase>
OpFluxImpl<NaturalMeshsetType<BLOCKSET>, 1, 1, A, I, OpBase>::OpFluxImpl(
    MoFEM::Interface &m_field, int ms_id, const std::string field_name)
    : OpFluxImpl<NaturalMeshsetType<UNKNOWNSET>, 1, 1, A, I, OpBase>(
          field_name) {
  CHK_THROW_MESSAGE(getMeshsetData(m_field, ms_id), "Get meshset data");
}

template <AssemblyType A, IntegrationType I, typename OpBase>
MoFEMErrorCode
OpFluxImpl<NaturalMeshsetType<BLOCKSET>, 1, 1, A, I, OpBase>::getMeshsetData(
    MoFEM::Interface &m_field, int ms_id) {
  MoFEMFunctionBegin;

  auto cubit_meshset_ptr =
      m_field.getInterface<MeshsetsManager>()->getCubitMeshsetPtr(ms_id,
                                                                  BLOCKSET);
  std::vector<double> block_data;
  CHKERR cubit_meshset_ptr->getAttributes(block_data);
  if (block_data.size() != 1) {
    SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
            "Expected that block has one attribute, e.g. heat flux value");
  }
  this->scalarValue = block_data[0];

  this->entsPtr = boost::make_shared<Range>();
  CHKERR m_field.get_moab().get_entities_by_handle(cubit_meshset_ptr->meshset,
                                                   *(this->entsPtr), true);

  MoFEMFunctionReturn(0);
}

template <AssemblyType A, IntegrationType I, typename OpBase>
OpFluxImpl<NaturalMeshsetType<HEATFLUXSET>, 1, 1, A, I, OpBase>::OpFluxImpl(
    MoFEM::Interface &m_field, int ms_id, const std::string field_name)
    : OpFluxImpl<NaturalMeshsetType<UNKNOWNSET>, 1, 1, A, I, OpBase>(
          field_name) {
  CHK_THROW_MESSAGE(getMeshsetData(m_field, ms_id), "Get meshset data");
}

template <AssemblyType A, IntegrationType I, typename OpBase>
MoFEMErrorCode
OpFluxImpl<NaturalMeshsetType<HEATFLUXSET>, 1, 1, A, I, OpBase>::getMeshsetData(
    MoFEM::Interface &m_field, int ms_id) {
  MoFEMFunctionBegin;

  auto cubit_meshset_ptr =
      m_field.getInterface<MeshsetsManager>()->getCubitMeshsetPtr(ms_id,
                                                                  SIDESET);
  HeatFluxCubitBcData heat_flux;
  CHKERR cubit_meshset_ptr->getBcDataStructure(heat_flux);
  this->scalarValue = heat_flux.data.value1;

  this->entsPtr = boost::make_shared<Range>();
  CHKERR m_field.get_moab().get_entities_by_handle(cubit_meshset_ptr->meshset,
                                                   *(this->entsPtr), true);

  MoFEMFunctionReturn(0);
}

template <int BASE_DIM, AssemblyType A, IntegrationType I, typename OpBase>
OpFluxImpl<NaturalMeshsetType<UNKNOWNSET>, BASE_DIM, BASE_DIM, A, I,
           OpBase>::OpFluxImpl(const std::string field_name, const double value,
                               boost::shared_ptr<Range> ents_ptr)
    : OpFluxImpl(field_name) {
  this->scalarValue = value;
  this->entsPtr = ents_ptr;
}

template <int BASE_DIM, AssemblyType A, IntegrationType I, typename OpBase>
OpFluxImpl<NaturalMeshsetType<UNKNOWNSET>, BASE_DIM, BASE_DIM, A, I,
           OpBase>::OpFluxImpl(const std::string field_name)
    : OpSource(field_name,

               [this](const double, const double, const double) {
                 return scalarValue;
               }

      ) {

  this->timeScalingFun = [this](const double t) {
    double s = 1;
    for (auto &o : vecOfTimeScalingMethods) {
      s *= o->getScale(t);
    }
    return s;
  };
}

template <int BASE_DIM, AssemblyType A, IntegrationType I, typename OpBase>
OpFluxImpl<NaturalMeshsetType<BLOCKSET>, BASE_DIM, BASE_DIM, A, I,
           OpBase>::OpFluxImpl(MoFEM::Interface &m_field, int ms_id,
                               const std::string field_name)
    : OpFluxImpl<NaturalMeshsetType<UNKNOWNSET>, BASE_DIM, BASE_DIM, A, I,
                 OpBase>(field_name) {
  CHK_THROW_MESSAGE(getMeshsetData(m_field, ms_id), "Get meshset data");
}

template <int BASE_DIM, AssemblyType A, IntegrationType I, typename OpBase>
MoFEMErrorCode OpFluxImpl<NaturalMeshsetType<BLOCKSET>, BASE_DIM, BASE_DIM, A,
                          I, OpBase>::getMeshsetData(MoFEM::Interface &m_field,
                                                     int ms_id) {
  MoFEMFunctionBegin;

  auto cubit_meshset_ptr =
      m_field.getInterface<MeshsetsManager>()->getCubitMeshsetPtr(ms_id,
                                                                  BLOCKSET);
  std::vector<double> attr_vec;
  cubit_meshset_ptr->getAttributes(attr_vec);
  if (attr_vec.size() != 1)
    SETERRQ(PETSC_COMM_SELF, MOFEM_INVALID_DATA, "Should be one attribute");
  this->scalarValue = attr_vec[0];

  this->entsPtr = boost::make_shared<Range>();
  CHKERR m_field.get_moab().get_entities_by_handle(cubit_meshset_ptr->meshset,
                                                   *(this->entsPtr), true);

  MoFEMFunctionReturn(0);
}

template <int BASE_DIM, AssemblyType A, IntegrationType I, typename OpBase>
OpFluxImpl<NaturalMeshsetType<TEMPERATURESET>, BASE_DIM, BASE_DIM, A, I,
           OpBase>::OpFluxImpl(MoFEM::Interface &m_field, int ms_id,
                               const std::string field_name)
    : OpFluxImpl<NaturalMeshsetType<UNKNOWNSET>, BASE_DIM, BASE_DIM, A, I,
                 OpBase>(field_name) {
  CHK_THROW_MESSAGE(getMeshsetData(m_field, ms_id), "Get meshset data");
}

template <int BASE_DIM, AssemblyType A, IntegrationType I, typename OpBase>
MoFEMErrorCode OpFluxImpl<NaturalMeshsetType<TEMPERATURESET>, BASE_DIM, BASE_DIM, A,
                          I, OpBase>::getMeshsetData(MoFEM::Interface &m_field,
                                                     int ms_id) {
  MoFEMFunctionBegin;

  auto cubit_meshset_ptr =
      m_field.getInterface<MeshsetsManager>()->getCubitMeshsetPtr(ms_id,
                                                                  NODESET);

  TemperatureCubitBcData mydata;
  CHKERR cubit_meshset_ptr->getBcDataStructure(mydata);
  this->scalarValue = mydata.data.value1;

  this->entsPtr = boost::make_shared<Range>();
  CHKERR m_field.get_moab().get_entities_by_handle(cubit_meshset_ptr->meshset,
                                                   *(this->entsPtr), true);

  MoFEMFunctionReturn(0);
}

template <CubitBC BCTYPE, int BASE_DIM, int FIELD_DIM, AssemblyType A,
          IntegrationType I, typename OpBase>
struct AddFluxToPipelineImpl<

    OpFluxImpl<NaturalMeshsetType<BCTYPE>, BASE_DIM, FIELD_DIM, A, I, OpBase>,
    A, I, OpBase

    > {

  AddFluxToPipelineImpl() = delete;

  static MoFEMErrorCode add(

      FluxOpType<OpFluxImpl<NaturalMeshsetType<BCTYPE>, BASE_DIM, FIELD_DIM, A,
                            I, OpBase>>,

      boost::ptr_vector<ForcesAndSourcesCore::UserDataOperator> &pipeline,
      MoFEM::Interface &m_field, const std::string field_name,
      std::vector<boost::shared_ptr<ScalingMethod>> smv, std::string block_name,
      Sev sev

  ) {
    MoFEMFunctionBegin;

    using OP = typename NaturalBC<OpBase>::template Assembly<
        A>::template LinearForm<I>::template OpFlux<NaturalMeshsetType<BCTYPE>,
                                                    BASE_DIM, FIELD_DIM>;

    auto add_op = [&](auto &&meshset_vec_ptr) {
      for (auto m : meshset_vec_ptr) {
        MOFEM_TAG_AND_LOG("SELF", sev, "OpFlux") << "Add for meshset " << *m;
        auto op_ptr = new OP(m_field, m->getMeshsetId(), field_name);
        for (auto sm : smv)
          op_ptr->getVecOfTimeScalingMethods().push_back(sm);
        pipeline.push_back(op_ptr);
      }
      MOFEM_LOG_CHANNEL("SELF");
    };

    switch (BCTYPE) {
    case TEMPERATURESET:
    case FORCESET:
      add_op(m_field.getInterface<MeshsetsManager>()->getCubitMeshsetPtr(
          NODESET | BCTYPE));
    case PRESSURESET:
    case HEATFLUXSET:
      add_op(m_field.getInterface<MeshsetsManager>()->getCubitMeshsetPtr(
          SIDESET | BCTYPE));
      break;
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

template <int BASE_DIM, int FIELD_DIM, AssemblyType A, IntegrationType I,
          typename OpBase>
struct AddFluxToPipelineImpl<

    OpFluxImpl<NaturalForceMeshsets, BASE_DIM, FIELD_DIM, A, I, OpBase>, A, I,
    OpBase

    > {

  AddFluxToPipelineImpl() = delete;

  static MoFEMErrorCode add(

      FluxOpType<
          OpFluxImpl<NaturalForceMeshsets, BASE_DIM, FIELD_DIM, A, I, OpBase>>,

      boost::ptr_vector<ForcesAndSourcesCore::UserDataOperator> &pipeline,
      MoFEM::Interface &m_field, const std::string field_name,
      std::vector<boost::shared_ptr<ScalingMethod>> smv,
      const std::string block_name, Sev sev

  ) {
    MoFEMFunctionBegin;

    using OpFluxForceset =
        typename NaturalBC<OpBase>::template Assembly<A>::template LinearForm<
            I>::template OpFlux<NaturalMeshsetType<FORCESET>, BASE_DIM,
                                FIELD_DIM>;
    using OpFluxPressureset =
        typename NaturalBC<OpBase>::template Assembly<A>::template LinearForm<
            I>::template OpFlux<NaturalMeshsetType<PRESSURESET>, BASE_DIM,
                                FIELD_DIM>;
    using OpFluxBlockset =
        typename NaturalBC<OpBase>::template Assembly<A>::template LinearForm<
            I>::template OpFlux<NaturalMeshsetType<BLOCKSET>, BASE_DIM,
                                FIELD_DIM>;

    CHKERR
    NaturalBC<OpBase>::template Assembly<A>::template LinearForm<
        I>::addFluxToPipeline(FluxOpType<OpFluxForceset>(), pipeline, m_field,
                              field_name, smv, block_name, sev);
    CHKERR
    NaturalBC<OpBase>::template Assembly<A>::template LinearForm<
        I>::addFluxToPipeline(FluxOpType<OpFluxPressureset>(), pipeline,
                              m_field, field_name, smv, block_name, sev);
    CHKERR
    NaturalBC<OpBase>::template Assembly<A>::template LinearForm<
        I>::addFluxToPipeline(FluxOpType<OpFluxBlockset>(), pipeline, m_field,
                              field_name, smv, block_name, sev);

    MoFEMFunctionReturn(0);
  }
};

template <int BASE_DIM, int FIELD_DIM, AssemblyType A, IntegrationType I,
          typename OpBase>
struct AddFluxToPipelineImpl<

    OpFluxImpl<NaturalTemperatureMeshsets, BASE_DIM, FIELD_DIM, A, I, OpBase>,
    A, I, OpBase

    > {

  AddFluxToPipelineImpl() = delete;

  static MoFEMErrorCode add(

      FluxOpType<OpFluxImpl<NaturalTemperatureMeshsets, BASE_DIM, FIELD_DIM, A,
                            I, OpBase>>,

      boost::ptr_vector<ForcesAndSourcesCore::UserDataOperator> &pipeline,
      MoFEM::Interface &m_field, const std::string field_name,
      std::vector<boost::shared_ptr<ScalingMethod>> smv,
      const std::string block_name, Sev sev

  ) {
    MoFEMFunctionBegin;

    using OpFluxTempSet =
        typename NaturalBC<OpBase>::template Assembly<A>::template LinearForm<
            I>::template OpFlux<NaturalMeshsetType<TEMPERATURESET>, BASE_DIM,
                                FIELD_DIM>;
    using OpFluxBlockset =
        typename NaturalBC<OpBase>::template Assembly<A>::template LinearForm<
            I>::template OpFlux<NaturalMeshsetType<BLOCKSET>, BASE_DIM,
                                FIELD_DIM>;

    CHKERR
    NaturalBC<OpBase>::template Assembly<A>::template LinearForm<
        I>::addFluxToPipeline(FluxOpType<OpFluxTempSet>(), pipeline, m_field,
                              field_name, smv, block_name, sev);
    CHKERR
    NaturalBC<OpBase>::template Assembly<A>::template LinearForm<
        I>::addFluxToPipeline(FluxOpType<OpFluxBlockset>(), pipeline, m_field,
                              field_name, smv, block_name, sev);

    MoFEMFunctionReturn(0);
  }
};

template <typename OpBase>
template <AssemblyType A>
template <IntegrationType I>
template <typename T>
MoFEMErrorCode NaturalBC<OpBase>::Assembly<A>::LinearForm<I>::addFluxToPipeline(

    FluxOpType<T>,

    boost::ptr_vector<ForcesAndSourcesCore::UserDataOperator> &pipeline,
    MoFEM::Interface &m_field, const std::string field_name,
    std::vector<boost::shared_ptr<ScalingMethod>> smv,
    const std::string block_name, Sev sev

) {
  using ADD =
      NaturalBC<OpBase>::Assembly<A>::LinearForm<I>::AddFluxToPipeline<T>;
  return ADD::add(FluxOpType<T>(), pipeline, m_field, field_name, smv,
                  block_name, sev);
}

} // namespace MoFEM

#endif //_NATURAL_BC_