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

using VecOfTimeScalingMethods = boost::ptr_vector<ScalingMethod>;

template <int BCSET, int BASE_DIM, int FIELD_DIM, AssemblyType A,
          IntegrationType I, typename OpBase>
struct OpFluxImpl;

template <int BASE_DIM, int FIELD_DIM, AssemblyType A, IntegrationType I,
          typename OpBase>
struct OpFluxImpl<UNKNOWNSET, BASE_DIM, FIELD_DIM, A, I, OpBase>
    : FormsIntegrators<OpBase>::template Assembly<A>::template LinearForm<
          I>::template OpSource<BASE_DIM, FIELD_DIM> {

  using OpSource = typename FormsIntegrators<OpBase>::template Assembly<
      A>::template LinearForm<I>::template OpSource<BASE_DIM, FIELD_DIM>;

  OpFluxImpl(const std::string field_name,
             FTensor::Tensor1<double, FIELD_DIM> t_force);

  OpFluxImpl(const std::string field_name);

  VecOfTimeScalingMethods &getVecOfTimeScalingMethods() {
    return vecOfTimeScalingMethods;
  }

protected:
  FTensor::Tensor1<double, FIELD_DIM> tForce;
  VecOfTimeScalingMethods vecOfTimeScalingMethods;
};

template <int BASE_DIM, int FIELD_DIM, AssemblyType A, IntegrationType I,
          typename OpBase>
struct OpFluxImpl<FORCESET, BASE_DIM, FIELD_DIM, A, I, OpBase>
    : OpFluxImpl<UNKNOWNSET, BASE_DIM, FIELD_DIM, A, I, OpBase> {

  OpFluxImpl(MoFEM::Interface &m_field, int ms_id,
             const std::string field_name);

protected:
  MoFEMErrorCode getMeshsetData(MoFEM::Interface &m_field, int ms_id);
};

template <int BASE_DIM, int FIELD_DIM, AssemblyType A, IntegrationType I,
          typename OpBase>
struct OpFluxImpl<BLOCKSET, BASE_DIM, FIELD_DIM, A, I, OpBase>
    : OpFluxImpl<UNKNOWNSET, BASE_DIM, FIELD_DIM, A, I, OpBase> {

  OpFluxImpl(MoFEM::Interface &m_field, int ms_id,
             const std::string field_name);

protected:
  MoFEMErrorCode getMeshsetData(MoFEM::Interface &m_field, int ms_id);
};

template <int BASE_DIM, int FIELD_DIM, AssemblyType A, IntegrationType I,
          typename OpBase>
struct OpFluxImpl<PRESSURESET, BASE_DIM, FIELD_DIM, A, I, OpBase>
    : OpFluxImpl<UNKNOWNSET, BASE_DIM, FIELD_DIM, A, I, OpBase> {

  using Parent = OpFluxImpl<UNKNOWNSET, BASE_DIM, FIELD_DIM, A, I, OpBase>;
  using EntData = EntitiesFieldData::EntData;

  OpFluxImpl(MoFEM::Interface &m_field, int ms_id,
             const std::string field_name);

protected:
  MoFEMErrorCode iNtegrate(EntData &data);
  MoFEMErrorCode getMeshsetData(MoFEM::Interface &m_field, int ms_id);

  double surfacePressure;  

};

/**
 * @brief Integrator forms
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

      /**
       * @brief This operator apply flux, or body force (i.e. body force)
       * @ingroup mofem_forms
       *
       * @tparam I
       */
      template <int BCTYPE, int BASE_DIM, int FIELD_DIM>
      struct OpFlux
          : public OpFluxImpl<BCTYPE, BASE_DIM, FIELD_DIM, A, I, EleOp> {
        using OpFluxImpl<BCTYPE, BASE_DIM, FIELD_DIM, A, I, EleOp>::OpFluxImpl;
      };
    };

  }; // Assembly
};

template <int BASE_DIM, int FIELD_DIM, AssemblyType A, IntegrationType I,
          typename OpBase>
OpFluxImpl<UNKNOWNSET, BASE_DIM, FIELD_DIM, A, I, OpBase>::OpFluxImpl(
    const std::string field_name, FTensor::Tensor1<double, FIELD_DIM> t_force)
    : OpFluxImpl(field_name) {

  FTensor::Index<'i', FIELD_DIM> i;
  this->tForce(i) = t_force(i);
}

template <int BASE_DIM, int FIELD_DIM, AssemblyType A, IntegrationType I,
          typename OpBase>
OpFluxImpl<UNKNOWNSET, BASE_DIM, FIELD_DIM, A, I, OpBase>::OpFluxImpl(
    const std::string field_name)
    : OpSource(
          field_name,

          [this](const double, const double, const double) { return tForce; }

      ) {

  this->timeScalingFun = [this](const double t) {
    double s = 1;
    for (auto o = vecOfTimeScalingMethods.begin();
         o != vecOfTimeScalingMethods.end(); ++o) {
      s *= o->getScale(t);
    }
    return s;
  };

  static_assert(FIELD_DIM > 1, "Is not implemented for scalar field");
}

template <int BASE_DIM, int FIELD_DIM, AssemblyType A, IntegrationType I,
          typename OpBase>
OpFluxImpl<FORCESET, BASE_DIM, FIELD_DIM, A, I, OpBase>::OpFluxImpl(
    MoFEM::Interface &m_field, int ms_id, const std::string field_name)
    : OpFluxImpl<UNKNOWNSET, BASE_DIM, FIELD_DIM, A, I, OpBase>(field_name) {}

template <int BASE_DIM, int FIELD_DIM, AssemblyType A, IntegrationType I,
          typename OpBase>
MoFEMErrorCode
OpFluxImpl<FORCESET, BASE_DIM, FIELD_DIM, A, I, OpBase>::getMeshsetData(
    MoFEM::Interface &m_field, int ms_id) {
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

template <int BASE_DIM, int FIELD_DIM, AssemblyType A, IntegrationType I,
          typename OpBase>
OpFluxImpl<BLOCKSET, BASE_DIM, FIELD_DIM, A, I, OpBase>::OpFluxImpl(
    MoFEM::Interface &m_field, int ms_id, const std::string field_name)
    : OpFluxImpl<UNKNOWNSET, BASE_DIM, FIELD_DIM, A, I, OpBase>(field_name) {}

template <int BASE_DIM, int FIELD_DIM, AssemblyType A, IntegrationType I,
          typename OpBase>
MoFEMErrorCode
OpFluxImpl<BLOCKSET, BASE_DIM, FIELD_DIM, A, I, OpBase>::getMeshsetData(
    MoFEM::Interface &m_field, int ms_id) {
  MoFEMFunctionBegin;

  auto cubit_meshset_ptr =
      m_field.getInterface<MeshsetsManager>()->getCubitMeshsetPtr(
          ms_id, NODESET);

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

  for (unsigned int ii = 0; ii != block_data.size(); ++ii) {
    this->tForce(ii) = block_data[ii];
  }

  this->entsPtr = boost::make_shared<Range>();
  CHKERR m_field.get_moab().get_entities_by_handle(cubit_meshset_ptr->meshset,
                                                   *(this->entsPtr), true);

  MoFEMFunctionReturn(0);
}

template <int BASE_DIM, int FIELD_DIM, AssemblyType A, IntegrationType I,
          typename OpBase>
OpFluxImpl<PRESSURESET, BASE_DIM, FIELD_DIM, A, I, OpBase>::OpFluxImpl(
    MoFEM::Interface &m_field, int ms_id, const std::string field_name)
    : OpFluxImpl<UNKNOWNSET, BASE_DIM, FIELD_DIM, A, I, OpBase>(field_name) {}

template <int BASE_DIM, int FIELD_DIM, AssemblyType A, IntegrationType I,
          typename OpBase>
MoFEMErrorCode
OpFluxImpl<PRESSURESET, BASE_DIM, FIELD_DIM, A, I, OpBase>::iNtegrate(
    EntData &data) {
  MoFEMFunctionBegin;

  auto t_normal = this->getFTensor1NormalsAtGaussPts();
  FTensor::Index<'i', FIELD_DIM> i;

  this->sourceFun = [&](const double, const double, const double) {
    this->tFroce(i) = this->surfacePressure*t_normal(i);
    ++t_normal;
    return this->tFroce;
  };

  CHKERR Parent::iNtegrate(data);
  MoFEMFunctionReturn(0);
}

template <int BASE_DIM, int FIELD_DIM, AssemblyType A, IntegrationType I,
          typename OpBase>
MoFEMErrorCode
OpFluxImpl<PRESSURESET, BASE_DIM, FIELD_DIM, A, I, OpBase>::getMeshsetData(
    MoFEM::Interface &m_field, int ms_id) {
  MoFEMFunctionBegin;

  auto cubit_meshset_ptr =
      m_field.getInterface<MeshsetsManager>()->getCubitMeshsetPtr(ms_id,
                                                                  BLOCKSET);

  PressureCubitBcData pressure_data;
  CHKERR cubit_meshset_ptr->getBcDataStructure(pressure_data);
  this->surfacePressure = pressure_data.data.value1;

  this->entsPtr = boost::make_shared<Range>();
  CHKERR m_field.get_moab().get_entities_by_handle(cubit_meshset_ptr->meshset,
                                                   *(this->entsPtr), true);

  MoFEMFunctionReturn(0);
}

}; // namespace MoFEM

#endif //_NATURAL_BC_