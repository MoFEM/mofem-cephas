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

/** Evaluate boundary integral on the right hand side*/
template <int FIELD_DIM, AssemblyType A, IntegrationType I, typename OpBase>
struct OpFluxRhsImpl<NaturalMeshsetType<FORCESET>, 1, FIELD_DIM, A, I, OpBase>
    : OpFluxRhsImpl<NaturalMeshsetType<UNKNOWNSET>, 1, FIELD_DIM, A, I,
                    OpBase> {

  OpFluxRhsImpl(
      MoFEM::Interface &m_field, int ms_id, const std::string field_name,
      std::vector<boost::shared_ptr<ScalingMethod>> smv,
      ScalarFun user_fun = [](double, double, double) { return 1; });

protected:
  /**
   * @brief Extract information stored on meshset (BLOCKSET)
   *
   * @param m_field
   * @param ms_id
   * @return MoFEMErrorCode
   */
  MoFEMErrorCode getMeshsetData(MoFEM::Interface &m_field, int ms_id);
};

template <int FIELD_DIM, AssemblyType A, IntegrationType I, typename OpBase>
struct OpFluxRhsImpl<NaturalMeshsetType<BLOCKSET>, 1, FIELD_DIM, A, I, OpBase>
    : OpFluxRhsImpl<NaturalMeshsetType<UNKNOWNSET>, 1, FIELD_DIM, A, I,
                    OpBase> {

  OpFluxRhsImpl(
      MoFEM::Interface &m_field, int ms_id, const std::string field_name,
      std::vector<boost::shared_ptr<ScalingMethod>> smv,
      ScalarFun user_fun = [](double, double, double) { return 1; });

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

  OpFluxRhsImpl(
      MoFEM::Interface &m_field, int ms_id, const std::string field_name,
      std::vector<boost::shared_ptr<ScalingMethod>> smv,
      ScalarFun user_fun = [](double, double, double) { return 1; });

protected:
  MoFEMErrorCode iNtegrate(EntData &data);
  MoFEMErrorCode getMeshsetData(MoFEM::Interface &m_field, int ms_id);
  double surfacePressure;
  ScalarFun userFun;
};

/**
 * @brief Base class for OpFluxRhsImpl<NaturalMeshsetType<T>, 1, FIELD_DIM, A,
 I, OpBase>
 *
 * This is only for scalar bases.
 *
 * @note It is derivitive from FormsIntegrators<OpBase>::template
 Assembly<A>::template LinearForm< I>::template OpSource<1, FIELD_DIM>
 *
 * @tparam FIELD_DIM field dimension
 * @tparam A
 * @tparam I
 * @tparam OpBase Base element operator for integration volume, face, edge, or
 vertex
 */
template <int FIELD_DIM, AssemblyType A, IntegrationType I, typename OpBase>
struct OpFluxRhsImpl<NaturalMeshsetType<UNKNOWNSET>, 1, FIELD_DIM, A, I, OpBase>
    : FormsIntegrators<OpBase>::template Assembly<A>::template LinearForm<
          I>::template OpSource<1, FIELD_DIM> {

  using OpSource = typename FormsIntegrators<OpBase>::template Assembly<
      A>::template LinearForm<I>::template OpSource<1, FIELD_DIM>;

  OpFluxRhsImpl(
      const std::string field_name, FTensor::Tensor1<double, FIELD_DIM> t_force,
      boost::shared_ptr<Range> ents_ptr,
      std::vector<boost::shared_ptr<ScalingMethod>> smv,
      ScalarFun user_fun = [](double, double, double) { return 1; });

protected:
  OpFluxRhsImpl(
      const std::string field_name,
      std::vector<boost::shared_ptr<ScalingMethod>> smv,
      ScalarFun user_fun = [](double, double, double) { return 1; });
  FTensor::Tensor1<double, FIELD_DIM> tForce;
  FTensor::Tensor1<double, FIELD_DIM> tScaledForce;
  VecOfTimeScalingMethods vecOfTimeScalingMethods;
  ScalarFun userFun;
};

template <AssemblyType A, IntegrationType I, typename OpBase>
struct OpFluxRhsImpl<NaturalMeshsetType<UNKNOWNSET>, 1, 1, A, I, OpBase>
    : FormsIntegrators<OpBase>::template Assembly<A>::template LinearForm<
          I>::template OpSource<1, 1> {

  using OpSource = typename FormsIntegrators<OpBase>::template Assembly<
      A>::template LinearForm<I>::template OpSource<1, 1>;

  OpFluxRhsImpl(
      const std::string field_name, double value,
      boost::shared_ptr<Range> ents_ptr,
      std::vector<boost::shared_ptr<ScalingMethod>> smv,
      ScalarFun user_fun = [](double, double, double) { return 1; });

protected:
  OpFluxRhsImpl(
      const std::string field_name,
      std::vector<boost::shared_ptr<ScalingMethod>> smv,
      ScalarFun user_fun = [](double, double, double) { return 1; });

  MoFEMErrorCode getMeshsetData(MoFEM::Interface &m_field, int ms_id);

  double scalarValue = 0;
  ScalarFun userFun;

  VecOfTimeScalingMethods vecOfTimeScalingMethods;
};

template <AssemblyType A, IntegrationType I, typename OpBase>
struct OpFluxRhsImpl<NaturalMeshsetType<BLOCKSET>, 1, 1, A, I, OpBase>
    : OpFluxRhsImpl<NaturalMeshsetType<UNKNOWNSET>, 1, 1, A, I, OpBase> {

  using Parent =
      OpFluxRhsImpl<NaturalMeshsetType<UNKNOWNSET>, 1, 1, A, I, OpBase>;

  OpFluxRhsImpl(
      MoFEM::Interface &m_field, int ms_id, const std::string field_name,
      std::vector<boost::shared_ptr<ScalingMethod>> smv,
      ScalarFun user_fun = [](double, double, double) { return 1; });

protected:
  MoFEMErrorCode getMeshsetData(MoFEM::Interface &m_field, int ms_id);
};

template <AssemblyType A, IntegrationType I, typename OpBase>
struct OpFluxRhsImpl<NaturalMeshsetType<HEATFLUXSET>, 1, 1, A, I, OpBase>
    : OpFluxRhsImpl<NaturalMeshsetType<UNKNOWNSET>, 1, 1, A, I, OpBase> {

  using Parent =
      OpFluxRhsImpl<NaturalMeshsetType<UNKNOWNSET>, 1, 1, A, I, OpBase>;

  OpFluxRhsImpl(
      MoFEM::Interface &m_field, int ms_id, const std::string field_name,
      std::vector<boost::shared_ptr<ScalingMethod>> smv,
      ScalarFun user_fun = [](double, double, double) { return 1; });

protected:
  MoFEMErrorCode getMeshsetData(MoFEM::Interface &m_field, int ms_id);
};

template <int FIELD_DIM, AssemblyType A, IntegrationType I, typename OpBase>
struct OpFluxRhsImpl<NaturalMeshsetType<UNKNOWNSET>, 3, FIELD_DIM, A, I, OpBase>
    : FormsIntegrators<OpBase>::template Assembly<A>::template LinearForm<
          I>::template OpNormalMixVecTimesScalar<FIELD_DIM> {

  using OpSource = typename FormsIntegrators<OpBase>::template Assembly<
      A>::template LinearForm<I>::template OpNormalMixVecTimesScalar<FIELD_DIM>;

  OpFluxRhsImpl(
      const std::string field_name, const double value,
      boost::shared_ptr<Range> ents_ptr,
      std::vector<boost::shared_ptr<ScalingMethod>> smv,
      ScalarFun user_fun = [](double, double, double) { return 1; });

protected:
  OpFluxRhsImpl(
      const std::string field_name,
      std::vector<boost::shared_ptr<ScalingMethod>> smv,
      ScalarFun user_fun = [](double, double, double) { return 1; });
  double scalarValue;
  VecOfTimeScalingMethods vecOfTimeScalingMethods;
  ScalarFun userFun;
};

template <int FIELD_DIM, AssemblyType A, IntegrationType I, typename OpBase>
struct OpFluxRhsImpl<NaturalMeshsetType<BLOCKSET>, 3, FIELD_DIM, A, I, OpBase>
    : OpFluxRhsImpl<NaturalMeshsetType<UNKNOWNSET>, 3, FIELD_DIM, A, I,
                    OpBase> {
  OpFluxRhsImpl(
      MoFEM::Interface &m_field, int ms_id, const std::string field_name,
      std::vector<boost::shared_ptr<ScalingMethod>> smv,
      ScalarFun user_fun = [](double, double, double) { return 1; });

protected:
  MoFEMErrorCode getMeshsetData(MoFEM::Interface &m_field, int ms_id);
};

template <int FIELD_DIM, AssemblyType A, IntegrationType I, typename OpBase>
struct OpFluxRhsImpl<NaturalMeshsetType<TEMPERATURESET>, 3, FIELD_DIM, A, I,
                     OpBase> : OpFluxRhsImpl<NaturalMeshsetType<UNKNOWNSET>, 3,
                                             FIELD_DIM, A, I, OpBase> {
  OpFluxRhsImpl(
      MoFEM::Interface &m_field, int ms_id, const std::string field_name,
      std::vector<boost::shared_ptr<ScalingMethod>> smv,
      ScalarFun user_fun = [](double, double, double) { return 1; });

protected:
  MoFEMErrorCode getMeshsetData(MoFEM::Interface &m_field, int ms_id);
};

template <CubitBC BCTYP> struct NaturalMeshsetTypeVectorScaling {};

template <int FIELD_DIM, AssemblyType A, IntegrationType I, typename OpBase>
struct OpFluxRhsImpl<NaturalMeshsetTypeVectorScaling<BLOCKSET>, 1, FIELD_DIM, A,
                     I, OpBase>
    : OpFluxRhsImpl<NaturalMeshsetTypeVectorScaling<UNKNOWNSET>, 1, FIELD_DIM,
                    A, I, OpBase> {

  OpFluxRhsImpl(
      MoFEM::Interface &m_field, int ms_id, const std::string field_name,
      std::vector<boost::shared_ptr<TimeScaleVector<FIELD_DIM>>> smv,
      ScalarFun user_fun = [](double, double, double) { return 1; });

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

  OpFluxRhsImpl(
      const std::string field_name, FTensor::Tensor1<double, FIELD_DIM> t_force,
      boost::shared_ptr<Range> ents_ptr,
      std::vector<boost::shared_ptr<TimeScaleVector<FIELD_DIM>>> smv,
      ScalarFun user_fun = [](double, double, double) { return 1; });

protected:
  OpFluxRhsImpl(
      const std::string field_name,
      std::vector<boost::shared_ptr<TimeScaleVector<FIELD_DIM>>> smv,
      ScalarFun user_fun = [](double, double, double) { return 1; });

  FTensor::Tensor1<double, FIELD_DIM> tForce;
  VecOfTimeVectorScalingMethods<FIELD_DIM> vecOfTimeVectorScalingMethods;
  ScalarFun userFun;
};

// Definitions

template <int FIELD_DIM, AssemblyType A, IntegrationType I, typename OpBase>
OpFluxRhsImpl<NaturalMeshsetType<UNKNOWNSET>, 1, FIELD_DIM, A, I, OpBase>::
    OpFluxRhsImpl(const std::string field_name,
                  FTensor::Tensor1<double, FIELD_DIM> t_force,
                  boost::shared_ptr<Range> ents_ptr,
                  std::vector<boost::shared_ptr<ScalingMethod>> smv,
                  ScalarFun user_fun)
    : OpFluxRhsImpl(field_name, smv, user_fun) {
  FTensor::Index<'i', FIELD_DIM> i;
  this->tForce(i) = t_force(i);
  this->entsPtr = ents_ptr;
}

template <int FIELD_DIM, AssemblyType A, IntegrationType I, typename OpBase>
OpFluxRhsImpl<NaturalMeshsetType<UNKNOWNSET>, 1, FIELD_DIM, A, I, OpBase>::
    OpFluxRhsImpl(const std::string field_name,
                  std::vector<boost::shared_ptr<ScalingMethod>> smv,
                  ScalarFun user_fun)
    : OpSource(field_name,

               [this](const double x, const double y, const double z) {
                 FTensor::Index<'i', FIELD_DIM> i;
                 tScaledForce(i) = tForce(i) * userFun(x, y, z);
                 return tScaledForce;
               }

               ),
      vecOfTimeScalingMethods(smv), userFun(user_fun) {
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
                  std::vector<boost::shared_ptr<ScalingMethod>> smv,
                  ScalarFun user_fun)
    : OpFluxRhsImpl<NaturalMeshsetType<UNKNOWNSET>, 1, FIELD_DIM, A, I, OpBase>(
          field_name, smv, user_fun) {
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
                  std::vector<boost::shared_ptr<ScalingMethod>> smv,
                  ScalarFun user_fun)
    : OpFluxRhsImpl<NaturalMeshsetType<UNKNOWNSET>, 1, FIELD_DIM, A, I, OpBase>(
          field_name, smv, user_fun) {
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

  MOFEM_LOG("WORLD", Sev::inform)
      << "tForce vector initialised: " << this->tForce;
  MOFEM_LOG("WORLD", Sev::inform)
      << "Number of elements " << this->entsPtr->size();

  MoFEMFunctionReturn(0);
}

template <int FIELD_DIM, AssemblyType A, IntegrationType I, typename OpBase>
OpFluxRhsImpl<NaturalMeshsetType<PRESSURESET>, 1, FIELD_DIM, A, I, OpBase>::
    OpFluxRhsImpl(MoFEM::Interface &m_field, int ms_id,
                  const std::string field_name,
                  std::vector<boost::shared_ptr<ScalingMethod>> smv,
                  ScalarFun user_fun)
    : OpFluxRhsImpl<NaturalMeshsetType<UNKNOWNSET>, 1, FIELD_DIM, A, I, OpBase>(
          field_name, smv, user_fun) {
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
                  std::vector<boost::shared_ptr<ScalingMethod>> smv,
                  ScalarFun user_fun)
    : OpFluxRhsImpl(field_name, smv, user_fun) {
  this->scalarValue = value;
  this->entsPtr = ents_ptr;
}

template <AssemblyType A, IntegrationType I, typename OpBase>
OpFluxRhsImpl<NaturalMeshsetType<UNKNOWNSET>, 1, 1, A, I, OpBase>::
    OpFluxRhsImpl(const std::string field_name,
                  std::vector<boost::shared_ptr<ScalingMethod>> smv,
                  ScalarFun user_fun)
    : OpSource(field_name,

               [this](const double x, const double y, const double z) {
                 return scalarValue * userFun(x, y, z);
               }

               ),
      vecOfTimeScalingMethods(smv), userFun(user_fun) {

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
    std::vector<boost::shared_ptr<ScalingMethod>> smv, ScalarFun user_fun)
    : OpFluxRhsImpl<NaturalMeshsetType<UNKNOWNSET>, 1, 1, A, I, OpBase>(
          field_name, smv, user_fun) {
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
                  std::vector<boost::shared_ptr<ScalingMethod>> smv,
                  ScalarFun user_fun)
    : OpFluxRhsImpl<NaturalMeshsetType<UNKNOWNSET>, 1, 1, A, I, OpBase>(
          field_name, smv, user_fun) {
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
                  std::vector<boost::shared_ptr<ScalingMethod>> smv,
                  ScalarFun user_fun)
    : OpFluxRhsImpl(field_name, smv) {
  this->scalarValue = value;
  this->entsPtr = ents_ptr;
}

template <int FIELD_DIM, AssemblyType A, IntegrationType I, typename OpBase>
OpFluxRhsImpl<NaturalMeshsetType<UNKNOWNSET>, 3, FIELD_DIM, A, I, OpBase>::
    OpFluxRhsImpl(const std::string field_name,
                  std::vector<boost::shared_ptr<ScalingMethod>> smv,
                  ScalarFun user_fun)
    : OpSource(field_name,

               [this](const double x, const double y, const double z) {
                 return scalarValue * userFun(x, y, z);
               }

               ),
      vecOfTimeScalingMethods(smv), userFun(user_fun) {

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
                  std::vector<boost::shared_ptr<ScalingMethod>> smv,
                  ScalarFun user_fun)
    : OpFluxRhsImpl<NaturalMeshsetType<UNKNOWNSET>, 3, FIELD_DIM, A, I, OpBase>(
          field_name, smv, user_fun) {
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
                  std::vector<boost::shared_ptr<ScalingMethod>> smv,
                  ScalarFun user_fun)
    : OpFluxRhsImpl<NaturalMeshsetType<UNKNOWNSET>, 3, FIELD_DIM, A, I, OpBase>(
          field_name, smv, user_fun) {
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

/**
 * @brief Function to push operators to rhs pipeline. This function user API.
 *
 * @tparam BCTYPE
 * @tparam BASE_DIM
 * @tparam FIELD_DIM
 * @tparam A
 * @tparam I
 * @tparam OpBase
 */
template <int FIELD_DIM, AssemblyType A, IntegrationType I, typename OpBase>
OpFluxRhsImpl<NaturalMeshsetTypeVectorScaling<UNKNOWNSET>, 1, FIELD_DIM, A, I,
              OpBase>::
    OpFluxRhsImpl(
        const std::string field_name,
        FTensor::Tensor1<double, FIELD_DIM> t_force,
        boost::shared_ptr<Range> ents_ptr,
        std::vector<boost::shared_ptr<TimeScaleVector<FIELD_DIM>>> smv,
        ScalarFun user_fun)
    : OpFluxRhsImpl(field_name, smv, user_fun) {
  FTensor::Index<'i', FIELD_DIM> i;
  this->tForce(i) = t_force(i);
  this->entsPtr = ents_ptr;
}

template <int FIELD_DIM, AssemblyType A, IntegrationType I, typename OpBase>
OpFluxRhsImpl<NaturalMeshsetTypeVectorScaling<UNKNOWNSET>, 1, FIELD_DIM, A, I,
              OpBase>::
    OpFluxRhsImpl(
        const std::string field_name,
        std::vector<boost::shared_ptr<TimeScaleVector<FIELD_DIM>>> smv,
        ScalarFun user_fun)
    : OpSource(field_name,
               [this](double x, double y, double z) {
                 FTensor::Index<'i', FIELD_DIM> i;
                 auto t = this->getFEMethod()->ts_t;
                 for (auto &o : vecOfTimeVectorScalingMethods) {
                   auto vec = o->getVector(t);
                   this->tForce(i) = vec(i) * userFun(x, y, z);
                 }
                 return this->tForce;
               }),
      vecOfTimeVectorScalingMethods(smv), userFun(user_fun) {

  static_assert(FIELD_DIM > 1, "Not implemented for scalar field");
}

template <int FIELD_DIM, AssemblyType A, IntegrationType I, typename OpBase>
OpFluxRhsImpl<NaturalMeshsetTypeVectorScaling<BLOCKSET>, 1, FIELD_DIM, A, I,
              OpBase>::
    OpFluxRhsImpl(
        MoFEM::Interface &m_field, int ms_id, const std::string field_name,
        std::vector<boost::shared_ptr<TimeScaleVector<FIELD_DIM>>> smv,
        ScalarFun user_fun)
    : OpFluxRhsImpl<NaturalMeshsetTypeVectorScaling<UNKNOWNSET>, 1, FIELD_DIM,
                    A, I, OpBase>(field_name, smv, user_fun) {
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
      ScalarFun user_fun, Sev sev

  ) {
    MoFEMFunctionBegin;

    using OP = typename NaturalBC<OpBase>::template Assembly<
        A>::template LinearForm<I>::template OpFlux<NaturalMeshsetType<BCTYPE>,
                                                    BASE_DIM, FIELD_DIM>;

    auto add_op = [&](auto &&meshset_vec_ptr) {
      for (auto m : meshset_vec_ptr) {
        MOFEM_TAG_AND_LOG("WORLD", sev, "OpFlux") << "Add " << *m;
        pipeline.push_back(
            new OP(m_field, m->getMeshsetId(), field_name, smv, user_fun));
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

  static MoFEMErrorCode
  add(boost::ptr_deque<ForcesAndSourcesCore::UserDataOperator> &pipeline,
      MoFEM::Interface &m_field, const std::string field_name,
      std::vector<boost::shared_ptr<ScalingMethod>> smv,
      std::string block_name, Sev sev) {
    return add(
        pipeline, m_field, field_name, smv, block_name,
        [](double, double, double) { return 1; }, sev);
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
      std::string block_name, ScalarFun user_fun, Sev sev

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
            new OPV(m_field, m->getMeshsetId(), field_name, vsmv, user_fun));
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

  static MoFEMErrorCode
  add(boost::ptr_deque<ForcesAndSourcesCore::UserDataOperator> &pipeline,
      MoFEM::Interface &m_field, const std::string field_name,
      std::vector<boost::shared_ptr<ScalingMethod>> smv,
      std::vector<boost::shared_ptr<TimeScaleVector<FIELD_DIM>>> vsmv,
      std::string block_name, Sev sev) {
    return add(
        pipeline, m_field, field_name, smv, vsmv, block_name,
        [](double, double, double) { return 1; }, sev);
  }
};

} // namespace MoFEM

#endif //_NATURAL_MESHSET_TYPE_HPP_