/**
 * @file NaturalMeshsetType.hpp
 * @brief Specialisation for NaturalMeshsetType
 * @version 0.13.2
 * @date 2022-09-18
 *
 * @copyright Copyright (c) 2022
 *
 */

#ifndef _NATURAL_MESHSET_TYPE_HPP_
#define _NATURAL_MESHSET_TYPE_HPP_

namespace MoFEM {

/**
 * @brief Type generating natural b.c. specialisations for meshsets
 *
 * @tparam BCTYP
 */
template <CubitBC BCTYP> struct NaturalMeshsetType {};

template <int FIELD_DIM, AssemblyType A, IntegrationType I, typename OpBase>
struct OpFluxRhsImpl<NaturalMeshsetType<FORCESET>, 1, FIELD_DIM, A, I, OpBase>
    : OpFluxRhsImpl<NaturalMeshsetType<UNKNOWNSET>, 1, FIELD_DIM, A, I,
                    OpBase> {

  OpFluxRhsImpl(MoFEM::Interface &m_field, int ms_id,
                const std::string field_name,
                std::vector<boost::shared_ptr<ScalingMethod>> smv);

protected:
  MoFEMErrorCode getMeshsetData(MoFEM::Interface &m_field, int ms_id);
};

template <int FIELD_DIM, AssemblyType A, IntegrationType I, typename OpBase>
struct OpFluxRhsImpl<NaturalMeshsetType<BLOCKSET>, 1, FIELD_DIM, A, I, OpBase>
    : OpFluxRhsImpl<NaturalMeshsetType<UNKNOWNSET>, 1, FIELD_DIM, A, I,
                    OpBase> {

  OpFluxRhsImpl(MoFEM::Interface &m_field, int ms_id,
                const std::string field_name,
                std::vector<boost::shared_ptr<ScalingMethod>> smv);

protected:
  MoFEMErrorCode getMeshsetData(MoFEM::Interface &m_field, int ms_id);
};

template <int FIELD_DIM, AssemblyType A, IntegrationType I, typename OpBase>
struct OpFluxRhsImpl<NaturalMeshsetType<PRESSURESET>, 1, FIELD_DIM, A, I,
                     OpBase> : OpFluxRhsImpl<NaturalMeshsetType<UNKNOWNSET>, 1,
                                             FIELD_DIM, A, I, OpBase> {

  using Parent =
      OpFluxRhsImpl<NaturalMeshsetType<UNKNOWNSET>, 1, FIELD_DIM, A, I, OpBase>;
  using EntData = EntitiesFieldData::EntData;

  OpFluxRhsImpl(MoFEM::Interface &m_field, int ms_id,
                const std::string field_name,
                std::vector<boost::shared_ptr<ScalingMethod>> smv);

protected:
  MoFEMErrorCode iNtegrate(EntData &data);
  MoFEMErrorCode getMeshsetData(MoFEM::Interface &m_field, int ms_id);

  double surfacePressure;
};

template <int FIELD_DIM, AssemblyType A, IntegrationType I, typename OpBase>
struct OpFluxRhsImpl<NaturalMeshsetType<UNKNOWNSET>, 1, FIELD_DIM, A, I, OpBase>
    : FormsIntegrators<OpBase>::template Assembly<A>::template LinearForm<
          I>::template OpSource<1, FIELD_DIM> {

  using OpSource = typename FormsIntegrators<OpBase>::template Assembly<
      A>::template LinearForm<I>::template OpSource<1, FIELD_DIM>;

  OpFluxRhsImpl(const std::string field_name,
                FTensor::Tensor1<double, FIELD_DIM> t_force,
                boost::shared_ptr<Range> ents_ptr,
                std::vector<boost::shared_ptr<ScalingMethod>> smv);

protected:
  OpFluxRhsImpl(const std::string field_name,
                std::vector<boost::shared_ptr<ScalingMethod>> smv);
  FTensor::Tensor1<double, FIELD_DIM> tForce;
  VecOfTimeScalingMethods vecOfTimeScalingMethods;
};

template <AssemblyType A, IntegrationType I, typename OpBase>
struct OpFluxRhsImpl<NaturalMeshsetType<UNKNOWNSET>, 1, 1, A, I, OpBase>
    : FormsIntegrators<OpBase>::template Assembly<A>::template LinearForm<
          I>::template OpSource<1, 1> {

  using OpSource = typename FormsIntegrators<OpBase>::template Assembly<
      A>::template LinearForm<I>::template OpSource<1, 1>;

  OpFluxRhsImpl(const std::string field_name, double value,
                boost::shared_ptr<Range> ents_ptr,
                std::vector<boost::shared_ptr<ScalingMethod>> smv);

protected:
  OpFluxRhsImpl(const std::string field_name,
                std::vector<boost::shared_ptr<ScalingMethod>> smv);

  MoFEMErrorCode getMeshsetData(MoFEM::Interface &m_field, int ms_id);
  double scalarValue;

  VecOfTimeScalingMethods vecOfTimeScalingMethods;
};

template <AssemblyType A, IntegrationType I, typename OpBase>
struct OpFluxRhsImpl<NaturalMeshsetType<BLOCKSET>, 1, 1, A, I, OpBase>
    : OpFluxRhsImpl<NaturalMeshsetType<UNKNOWNSET>, 1, 1, A, I, OpBase> {

  using Parent =
      OpFluxRhsImpl<NaturalMeshsetType<UNKNOWNSET>, 1, 1, A, I, OpBase>;

  OpFluxRhsImpl(MoFEM::Interface &m_field, int ms_id,
                const std::string field_name,
                std::vector<boost::shared_ptr<ScalingMethod>> smv);

protected:
  MoFEMErrorCode getMeshsetData(MoFEM::Interface &m_field, int ms_id);
};

template <AssemblyType A, IntegrationType I, typename OpBase>
struct OpFluxRhsImpl<NaturalMeshsetType<HEATFLUXSET>, 1, 1, A, I, OpBase>
    : OpFluxRhsImpl<NaturalMeshsetType<UNKNOWNSET>, 1, 1, A, I, OpBase> {

  using Parent =
      OpFluxRhsImpl<NaturalMeshsetType<UNKNOWNSET>, 1, 1, A, I, OpBase>;

  OpFluxRhsImpl(MoFEM::Interface &m_field, int ms_id,
                const std::string field_name,
                std::vector<boost::shared_ptr<ScalingMethod>> smv);

protected:
  MoFEMErrorCode getMeshsetData(MoFEM::Interface &m_field, int ms_id);
};

template <int FIELD_DIM, AssemblyType A, IntegrationType I, typename OpBase>
struct OpFluxRhsImpl<NaturalMeshsetType<UNKNOWNSET>, 3, FIELD_DIM, A, I, OpBase>
    : FormsIntegrators<OpBase>::template Assembly<A>::template LinearForm<
          I>::template OpNormalMixVecTimesScalar<FIELD_DIM> {

  using OpSource = typename FormsIntegrators<OpBase>::template Assembly<
      A>::template LinearForm<I>::template OpNormalMixVecTimesScalar<FIELD_DIM>;

  OpFluxRhsImpl(const std::string field_name, const double value,
                boost::shared_ptr<Range> ents_ptr,
                std::vector<boost::shared_ptr<ScalingMethod>> smv);

protected:
  OpFluxRhsImpl(const std::string field_name,
                std::vector<boost::shared_ptr<ScalingMethod>> smv);
  double scalarValue;

  VecOfTimeScalingMethods vecOfTimeScalingMethods;
};

template <int FIELD_DIM, AssemblyType A, IntegrationType I, typename OpBase>
struct OpFluxRhsImpl<NaturalMeshsetType<BLOCKSET>, 3, FIELD_DIM, A, I, OpBase>
    : OpFluxRhsImpl<NaturalMeshsetType<UNKNOWNSET>, 3, FIELD_DIM, A, I,
                    OpBase> {
  OpFluxRhsImpl(MoFEM::Interface &m_field, int ms_id,
                const std::string field_name,
                std::vector<boost::shared_ptr<ScalingMethod>> smv);

protected:
  MoFEMErrorCode getMeshsetData(MoFEM::Interface &m_field, int ms_id);
};

template <int FIELD_DIM, AssemblyType A, IntegrationType I, typename OpBase>
struct OpFluxRhsImpl<NaturalMeshsetType<TEMPERATURESET>, 3, FIELD_DIM, A, I,
                     OpBase> : OpFluxRhsImpl<NaturalMeshsetType<UNKNOWNSET>, 3,
                                             FIELD_DIM, A, I, OpBase> {
  OpFluxRhsImpl(MoFEM::Interface &m_field, int ms_id,
                const std::string field_name,
                std::vector<boost::shared_ptr<ScalingMethod>> smv);

protected:
  MoFEMErrorCode getMeshsetData(MoFEM::Interface &m_field, int ms_id);
};

template <CubitBC BCTYP> struct NaturalMeshsetTypeVectorScaling {};

template <int FIELD_DIM, AssemblyType A, IntegrationType I, typename OpBase>
struct OpFluxRhsImpl<NaturalMeshsetTypeVectorScaling<BLOCKSET>, 1, FIELD_DIM, A,
                     I, OpBase>
    : OpFluxRhsImpl<NaturalMeshsetTypeVectorScaling<UNKNOWNSET>, 1, FIELD_DIM,
                    A, I, OpBase> {

  OpFluxRhsImpl(MoFEM::Interface &m_field, int ms_id,
                const std::string field_name,
                std::vector<boost::shared_ptr<TimeScaleVector<FIELD_DIM>>> smv);

protected:
  MoFEMErrorCode getMeshsetData(MoFEM::Interface &m_field, int ms_id);
};

template <int FIELD_DIM, AssemblyType A, IntegrationType I, typename OpBase>
struct OpFluxRhsImpl<NaturalMeshsetTypeVectorScaling<UNKNOWNSET>, 1, FIELD_DIM,
                     A, I, OpBase>
    : FormsIntegrators<OpBase>::template Assembly<A>::template LinearForm<
          I>::template OpSource<1, FIELD_DIM> {

  using OpSource = typename FormsIntegrators<OpBase>::template Assembly<
      A>::template LinearForm<I>::template OpSource<1, FIELD_DIM>;

  OpFluxRhsImpl(const std::string field_name,
                FTensor::Tensor1<double, FIELD_DIM> t_force,
                boost::shared_ptr<Range> ents_ptr,
                std::vector<boost::shared_ptr<TimeScaleVector<FIELD_DIM>>> smv);

protected:
  OpFluxRhsImpl(const std::string field_name,
                std::vector<boost::shared_ptr<TimeScaleVector<FIELD_DIM>>> smv);

  FTensor::Tensor1<double, FIELD_DIM> tForce;
  VecOfTimeVectorScalingMethods<FIELD_DIM> vecOfTimeVectorScalingMethods;
};

// Definitions

template <int FIELD_DIM, AssemblyType A, IntegrationType I, typename OpBase>
OpFluxRhsImpl<NaturalMeshsetType<UNKNOWNSET>, 1, FIELD_DIM, A, I, OpBase>::
    OpFluxRhsImpl(const std::string field_name,
                  FTensor::Tensor1<double, FIELD_DIM> t_force,
                  boost::shared_ptr<Range> ents_ptr,
                  std::vector<boost::shared_ptr<ScalingMethod>> smv)
    : OpFluxRhsImpl(field_name, smv) {
  FTensor::Index<'i', FIELD_DIM> i;
  this->tForce(i) = t_force(i);
  this->entsPtr = ents_ptr;
}

template <int FIELD_DIM, AssemblyType A, IntegrationType I, typename OpBase>
OpFluxRhsImpl<NaturalMeshsetType<UNKNOWNSET>, 1, FIELD_DIM, A, I, OpBase>::
    OpFluxRhsImpl(const std::string field_name,
                  std::vector<boost::shared_ptr<ScalingMethod>> smv)
    : OpSource(
          field_name,

          [this](const double, const double, const double) { return tForce; }

          ),
      vecOfTimeScalingMethods(smv) {
  this->timeScalingFun = [this](const double t) {
    double s = 1;
    for (auto &o : vecOfTimeScalingMethods) {
      s *= o->getScale(t);
    }
    return s;
  };

  static_assert(FIELD_DIM > 1, "Not implemented for scalar field");
}

template <int FIELD_DIM, AssemblyType A, IntegrationType I, typename OpBase>
OpFluxRhsImpl<NaturalMeshsetType<FORCESET>, 1, FIELD_DIM, A, I, OpBase>::
    OpFluxRhsImpl(MoFEM::Interface &m_field, int ms_id,
                  const std::string field_name,
                  std::vector<boost::shared_ptr<ScalingMethod>> smv)
    : OpFluxRhsImpl<NaturalMeshsetType<UNKNOWNSET>, 1, FIELD_DIM, A, I, OpBase>(
          field_name, smv) {
  CHK_THROW_MESSAGE(getMeshsetData(m_field, ms_id), "Get meshset data");
}

template <int FIELD_DIM, AssemblyType A, IntegrationType I, typename OpBase>
MoFEMErrorCode OpFluxRhsImpl<NaturalMeshsetType<FORCESET>, 1, FIELD_DIM, A, I,
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
OpFluxRhsImpl<NaturalMeshsetType<BLOCKSET>, 1, FIELD_DIM, A, I, OpBase>::
    OpFluxRhsImpl(MoFEM::Interface &m_field, int ms_id,
                  const std::string field_name,
                  std::vector<boost::shared_ptr<ScalingMethod>> smv)
    : OpFluxRhsImpl<NaturalMeshsetType<UNKNOWNSET>, 1, FIELD_DIM, A, I, OpBase>(
          field_name, smv) {
  CHK_THROW_MESSAGE(getMeshsetData(m_field, ms_id), "Get meshset data");
}

template <int FIELD_DIM, AssemblyType A, IntegrationType I, typename OpBase>
MoFEMErrorCode OpFluxRhsImpl<NaturalMeshsetType<BLOCKSET>, 1, FIELD_DIM, A, I,
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
OpFluxRhsImpl<NaturalMeshsetType<PRESSURESET>, 1, FIELD_DIM, A, I, OpBase>::
    OpFluxRhsImpl(MoFEM::Interface &m_field, int ms_id,
                  const std::string field_name,
                  std::vector<boost::shared_ptr<ScalingMethod>> smv)
    : OpFluxRhsImpl<NaturalMeshsetType<UNKNOWNSET>, 1, FIELD_DIM, A, I, OpBase>(
          field_name, smv) {
  CHK_THROW_MESSAGE(getMeshsetData(m_field, ms_id), "Get meshset data");
}

template <int FIELD_DIM, AssemblyType A, IntegrationType I, typename OpBase>
MoFEMErrorCode OpFluxRhsImpl<NaturalMeshsetType<PRESSURESET>, 1, FIELD_DIM, A,
                             I, OpBase>::iNtegrate(EntData &data) {
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
MoFEMErrorCode
OpFluxRhsImpl<NaturalMeshsetType<PRESSURESET>, 1, FIELD_DIM, A, I,
              OpBase>::getMeshsetData(MoFEM::Interface &m_field, int ms_id) {
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
OpFluxRhsImpl<NaturalMeshsetType<UNKNOWNSET>, 1, 1, A, I, OpBase>::
    OpFluxRhsImpl(const std::string field_name, const double value,
                  boost::shared_ptr<Range> ents_ptr,
                  std::vector<boost::shared_ptr<ScalingMethod>> smv)
    : OpFluxRhsImpl(field_name, smv) {
  this->scalarValue = value;
  this->entsPtr = ents_ptr;
}

template <AssemblyType A, IntegrationType I, typename OpBase>
OpFluxRhsImpl<NaturalMeshsetType<UNKNOWNSET>, 1, 1, A, I, OpBase>::
    OpFluxRhsImpl(const std::string field_name,
                  std::vector<boost::shared_ptr<ScalingMethod>> smv)
    : OpSource(field_name,

               [this](const double, const double, const double) {
                 return scalarValue;
               }

               ),
      vecOfTimeScalingMethods(smv) {

  this->timeScalingFun = [this](const double t) {
    double s = 1;
    for (auto &o : vecOfTimeScalingMethods) {
      s *= o->getScale(t);
    }
    return s;
  };
}

template <AssemblyType A, IntegrationType I, typename OpBase>
OpFluxRhsImpl<NaturalMeshsetType<BLOCKSET>, 1, 1, A, I, OpBase>::OpFluxRhsImpl(
    MoFEM::Interface &m_field, int ms_id, const std::string field_name,
    std::vector<boost::shared_ptr<ScalingMethod>> smv)
    : OpFluxRhsImpl<NaturalMeshsetType<UNKNOWNSET>, 1, 1, A, I, OpBase>(
          field_name, smv) {
  CHK_THROW_MESSAGE(getMeshsetData(m_field, ms_id), "Get meshset data");
}

template <AssemblyType A, IntegrationType I, typename OpBase>
MoFEMErrorCode
OpFluxRhsImpl<NaturalMeshsetType<BLOCKSET>, 1, 1, A, I, OpBase>::getMeshsetData(
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

  MOFEM_LOG("WORLD", Sev::inform)
      << "Flux blockset " << cubit_meshset_ptr->getName();
  MOFEM_LOG("WORLD", Sev::inform)
      << "Scalar attribute value: " << this->scalarValue;

  MoFEMFunctionReturn(0);
}

template <AssemblyType A, IntegrationType I, typename OpBase>
OpFluxRhsImpl<NaturalMeshsetType<HEATFLUXSET>, 1, 1, A, I, OpBase>::
    OpFluxRhsImpl(MoFEM::Interface &m_field, int ms_id,
                  const std::string field_name,
                  std::vector<boost::shared_ptr<ScalingMethod>> smv)
    : OpFluxRhsImpl<NaturalMeshsetType<UNKNOWNSET>, 1, 1, A, I, OpBase>(
          field_name, smv) {
  CHK_THROW_MESSAGE(getMeshsetData(m_field, ms_id), "Get meshset data");
}

template <AssemblyType A, IntegrationType I, typename OpBase>
MoFEMErrorCode OpFluxRhsImpl<NaturalMeshsetType<HEATFLUXSET>, 1, 1, A, I,
                             OpBase>::getMeshsetData(MoFEM::Interface &m_field,
                                                     int ms_id) {
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

template <int FIELD_DIM, AssemblyType A, IntegrationType I, typename OpBase>
OpFluxRhsImpl<NaturalMeshsetType<UNKNOWNSET>, 3, FIELD_DIM, A, I, OpBase>::
    OpFluxRhsImpl(const std::string field_name, const double value,
                  boost::shared_ptr<Range> ents_ptr,
                  std::vector<boost::shared_ptr<ScalingMethod>> smv)
    : OpFluxRhsImpl(field_name, smv) {
  this->scalarValue = value;
  this->entsPtr = ents_ptr;
}

template <int FIELD_DIM, AssemblyType A, IntegrationType I, typename OpBase>
OpFluxRhsImpl<NaturalMeshsetType<UNKNOWNSET>, 3, FIELD_DIM, A, I, OpBase>::
    OpFluxRhsImpl(const std::string field_name,
                  std::vector<boost::shared_ptr<ScalingMethod>> smv)
    : OpSource(field_name,

               [this](const double, const double, const double) {
                 return scalarValue;
               }

               ),
      vecOfTimeScalingMethods(smv) {

  this->timeScalingFun = [this](const double t) {
    double s = 1;
    for (auto &o : vecOfTimeScalingMethods) {
      s *= o->getScale(t);
    }
    return s;
  };
}

template <int FIELD_DIM, AssemblyType A, IntegrationType I, typename OpBase>
OpFluxRhsImpl<NaturalMeshsetType<BLOCKSET>, 3, FIELD_DIM, A, I, OpBase>::
    OpFluxRhsImpl(MoFEM::Interface &m_field, int ms_id,
                  const std::string field_name,
                  std::vector<boost::shared_ptr<ScalingMethod>> smv)
    : OpFluxRhsImpl<NaturalMeshsetType<UNKNOWNSET>, 3, FIELD_DIM, A, I, OpBase>(
          field_name, smv) {
  CHK_THROW_MESSAGE(getMeshsetData(m_field, ms_id), "Get meshset data");
}

template <int FIELD_DIM, AssemblyType A, IntegrationType I, typename OpBase>
MoFEMErrorCode OpFluxRhsImpl<NaturalMeshsetType<BLOCKSET>, 3, FIELD_DIM, A, I,
                             OpBase>::getMeshsetData(MoFEM::Interface &m_field,
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

template <int FIELD_DIM, AssemblyType A, IntegrationType I, typename OpBase>
OpFluxRhsImpl<NaturalMeshsetType<TEMPERATURESET>, 3, FIELD_DIM, A, I, OpBase>::
    OpFluxRhsImpl(MoFEM::Interface &m_field, int ms_id,
                  const std::string field_name,
                  std::vector<boost::shared_ptr<ScalingMethod>> smv)
    : OpFluxRhsImpl<NaturalMeshsetType<UNKNOWNSET>, 3, FIELD_DIM, A, I, OpBase>(
          field_name, smv) {
  CHK_THROW_MESSAGE(getMeshsetData(m_field, ms_id), "Get meshset data");
}

template <int FIELD_DIM, AssemblyType A, IntegrationType I, typename OpBase>
MoFEMErrorCode
OpFluxRhsImpl<NaturalMeshsetType<TEMPERATURESET>, 3, FIELD_DIM, A, I,
              OpBase>::getMeshsetData(MoFEM::Interface &m_field, int ms_id) {
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

template <int FIELD_DIM, AssemblyType A, IntegrationType I, typename OpBase>
OpFluxRhsImpl<NaturalMeshsetTypeVectorScaling<UNKNOWNSET>, 1, FIELD_DIM, A, I,
              OpBase>::
    OpFluxRhsImpl(
        const std::string field_name,
        FTensor::Tensor1<double, FIELD_DIM> t_force,
        boost::shared_ptr<Range> ents_ptr,
        std::vector<boost::shared_ptr<TimeScaleVector<FIELD_DIM>>> smv)
    : OpFluxRhsImpl(field_name, smv) {
  FTensor::Index<'i', FIELD_DIM> i;
  this->tForce(i) = t_force(i);
  this->entsPtr = ents_ptr;
}

template <int FIELD_DIM, AssemblyType A, IntegrationType I, typename OpBase>
OpFluxRhsImpl<NaturalMeshsetTypeVectorScaling<UNKNOWNSET>, 1, FIELD_DIM, A, I,
              OpBase>::
    OpFluxRhsImpl(
        const std::string field_name,
        std::vector<boost::shared_ptr<TimeScaleVector<FIELD_DIM>>> smv)
    : OpSource(field_name,
               [this](double, double, double) {
                 FTensor::Index<'i', FIELD_DIM> i;
                 auto t = this->getFEMethod()->ts_t;
                 for (auto &o : vecOfTimeVectorScalingMethods) {
                   auto vec = o->getVector(t);
                   this->tForce(i) = vec(i);
                 }
                 return this->tForce;
               }),
      vecOfTimeVectorScalingMethods(smv) {

  static_assert(FIELD_DIM > 1, "Not implemented for scalar field");
}

template <int FIELD_DIM, AssemblyType A, IntegrationType I, typename OpBase>
OpFluxRhsImpl<NaturalMeshsetTypeVectorScaling<BLOCKSET>, 1, FIELD_DIM, A, I,
              OpBase>::
    OpFluxRhsImpl(
        MoFEM::Interface &m_field, int ms_id, const std::string field_name,
        std::vector<boost::shared_ptr<TimeScaleVector<FIELD_DIM>>> smv)
    : OpFluxRhsImpl<NaturalMeshsetTypeVectorScaling<UNKNOWNSET>, 1, FIELD_DIM,
                    A, I, OpBase>(field_name, smv) {
  CHK_THROW_MESSAGE(getMeshsetData(m_field, ms_id), "Get meshset data");
}

template <int FIELD_DIM, AssemblyType A, IntegrationType I, typename OpBase>
MoFEMErrorCode
OpFluxRhsImpl<NaturalMeshsetTypeVectorScaling<BLOCKSET>, 1, FIELD_DIM, A, I,
              OpBase>::getMeshsetData(MoFEM::Interface &m_field, int ms_id) {
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

  for (unsigned int ii = 0; ii != FIELD_DIM; ++ii) {
    this->tForce(ii) = block_data[ii];
  }

  MOFEM_LOG("WORLD", Sev::inform)
      << "Flux blockset " << cubit_meshset_ptr->getName();
  MOFEM_LOG("WORLD", Sev::inform)
      << "Number of attributes " << block_data.size();
  MOFEM_LOG("WORLD", Sev::inform)
      << "tForce vector initialised: " << this->tForce;

  this->entsPtr = boost::make_shared<Range>();
  CHKERR m_field.get_moab().get_entities_by_handle(cubit_meshset_ptr->meshset,
                                                   *(this->entsPtr), true);

  MoFEMFunctionReturn(0);
}

template <CubitBC BCTYPE, int BASE_DIM, int FIELD_DIM, AssemblyType A,
          IntegrationType I, typename OpBase>
struct AddFluxToRhsPipelineImpl<

    OpFluxRhsImpl<NaturalMeshsetType<BCTYPE>, BASE_DIM, FIELD_DIM, A, I,
                  OpBase>,
    A, I, OpBase

    > {

  AddFluxToRhsPipelineImpl() = delete;

  static MoFEMErrorCode add(

      boost::ptr_deque<ForcesAndSourcesCore::UserDataOperator> &pipeline,
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
        MOFEM_TAG_AND_LOG("WORLD", sev, "OpFlux") << "Add " << *m;
        pipeline.push_back(new OP(m_field, m->getMeshsetId(), field_name, smv));
      }
      MOFEM_LOG_CHANNEL("WORLD");
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

template <CubitBC BCTYPE, int BASE_DIM, int FIELD_DIM, AssemblyType A,
          IntegrationType I, typename OpBase>
struct AddFluxToRhsPipelineImpl<

    OpFluxRhsImpl<NaturalMeshsetTypeVectorScaling<BCTYPE>, BASE_DIM, FIELD_DIM,
                  A, I, OpBase>,
    A, I, OpBase

    > {

  AddFluxToRhsPipelineImpl() = delete;

  static MoFEMErrorCode add(

      boost::ptr_deque<ForcesAndSourcesCore::UserDataOperator> &pipeline,
      MoFEM::Interface &m_field, const std::string field_name,
      std::vector<boost::shared_ptr<ScalingMethod>> smv,
      std::vector<boost::shared_ptr<TimeScaleVector<FIELD_DIM>>> vsmv,
      std::string block_name, Sev sev

  ) {
    MoFEMFunctionBegin;

    using OPV =
        typename NaturalBC<OpBase>::template Assembly<A>::template LinearForm<
            I>::template OpFlux<NaturalMeshsetTypeVectorScaling<BCTYPE>,
                                BASE_DIM, FIELD_DIM>;

    auto add_opv = [&](auto &&meshset_vec_ptr) {
      for (auto m : meshset_vec_ptr) {
        MOFEM_TAG_AND_LOG("WORLD", sev, "OpFlux") << "Add " << *m;
        pipeline.push_back(
            new OPV(m_field, m->getMeshsetId(), field_name, vsmv));
      }
      MOFEM_LOG_CHANNEL("WORLD");
    };

    switch (BCTYPE) {
    case BLOCKSET:
      add_opv(

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

} // namespace MoFEM

#endif //_NATURAL_MESHSET_TYPE_HPP_