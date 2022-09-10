/**
 * @file Essential.hpp
 * @brief Setting essential boundary conditions
 * @version 0.13.1
 * @date 2022-08-12
 *
 * @copyright Copyright (c) 2022
 *
 */

#ifndef _ESSENTIAL_BC_
#define _ESSENTIAL_BC_

namespace MoFEM {

template <CubitBC BC> struct EssentialMeshsetType {};

template <typename T> struct EssentialPreProc {};

template <> struct EssentialPreProc<DisplacementCubitBcData> {
  EssentialPreProc(MoFEM::Interface &m_field,
                   boost::shared_ptr<FEMethod> fe_ptr, bool get_coords = false);

  MoFEMErrorCode operator()();

  VecOfTimeScalingMethods &getVecOfTimeScalingMethods();

protected:
  MoFEM::Interface &mField;
  boost::weak_ptr<FEMethod> fePtr;
  VecOfTimeScalingMethods vecOfTimeScalingMethods;
  bool getCoords;
};

template <> struct EssentialPreProc<TemperatureCubitBcData> {
  EssentialPreProc(MoFEM::Interface &m_field,
                   boost::shared_ptr<FEMethod> fe_ptr);

  MoFEMErrorCode operator()();

  VecOfTimeScalingMethods &getVecOfTimeScalingMethods();

protected:
  MoFEM::Interface &mField;
  boost::weak_ptr<FEMethod> fePtr;
  VecOfTimeScalingMethods vecOfTimeScalingMethods;
};

template <typename T, int BASE_DIM, int FIELD_DIM, AssemblyType A,
          IntegrationType I, typename OpBase>
struct OpEssentialRhsImpl;

template <typename T, int BASE_DIM, int FIELD_DIM, AssemblyType A,
          IntegrationType I, typename OpBase>
struct OpEssentialLhsImpl;

template <typename T, AssemblyType A, IntegrationType I, typename OpBase>
struct AddEssentialToRhsPipelineImpl;

template <typename T, AssemblyType A, IntegrationType I, typename OpBase>
struct AddEssentialToLhsPipelineImpl;

/**
 * @brief Natural boundary conditions
 * @ingroup mofem_forms
 *
 * @tparam EleOp
 */
template <typename EleOp> struct EssentialBC {

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
      using OpEssentialRhs =
          OpEssentialRhsImpl<T, BASE_DIM, FIELD_DIM, A, I, EleOp>;
    };

    template <IntegrationType I> struct BiLinearForm {
      template <typename T, int BASE_DIM, int FIELD_DIM>
      using OpEssentialLhs =
          OpEssentialLhsImpl<T, BASE_DIM, FIELD_DIM, A, I, EleOp>;
    };
  };
};

template <int FIELD_DIM, AssemblyType A, IntegrationType I, typename OpBase>
struct OpEssentialRhsImpl<DisplacementCubitBcData, 1, FIELD_DIM, A, I, OpBase>
    : FormsIntegrators<OpBase>::template Assembly<A>::template LinearForm<
          I>::template OpSource<1, FIELD_DIM> {

  using OpSource = typename FormsIntegrators<OpBase>::template Assembly<
      A>::template LinearForm<I>::template OpSource<1, FIELD_DIM>;

  OpEssentialRhsImpl(const std::string field_name,
                     boost::shared_ptr<DisplacementCubitBcData> bc_data,
                     boost::shared_ptr<Range> ents_ptr);

private:
  FTensor::Tensor1<double, FIELD_DIM> tVal;
};

template <int FIELD_DIM, AssemblyType A, IntegrationType I, typename OpBase>
OpEssentialRhsImpl<DisplacementCubitBcData, 1, FIELD_DIM, A, I, OpBase>::
    OpEssentialRhsImpl(const std::string field_name,
                       boost::shared_ptr<DisplacementCubitBcData> bc_data,
                       boost::shared_ptr<Range> ents_ptr)
    : OpSource(
          field_name, [this](double, double, double) { return tVal; },
          ents_ptr) {
  static_assert(FIELD_DIM > 1, "Is not implemented for scalar field");

  FTensor::Index<'i', FIELD_DIM> i;
  tVal(i) = 0;
  if (bc_data->data.flag1 == 1)
    tVal(0) = bc_data->data.value1;
  if (bc_data->data.flag2 == 1 && FIELD_DIM > 1)
    tVal(1) = bc_data->data.value2;
  if (bc_data->data.flag3 == 1 && FIELD_DIM > 2)
    tVal(2) = bc_data->data.value3;
}

template <int BASE_DIM, int FIELD_DIM, AssemblyType A, IntegrationType I,
          typename OpBase>
struct OpEssentialLhsImpl<DisplacementCubitBcData, BASE_DIM, FIELD_DIM, A, I,
                          OpBase>
    : FormsIntegrators<OpBase>::template Assembly<A>::template BiLinearForm<
          I>::template OpMass<BASE_DIM, FIELD_DIM> {

  using OpMass = typename FormsIntegrators<OpBase>::template Assembly<
      A>::template BiLinearForm<I>::template OpMass<BASE_DIM, FIELD_DIM>;

  OpEssentialLhsImpl(const std::string field_name,
                     boost::shared_ptr<Range> ents_ptr);
};

template <int BASE_DIM, int FIELD_DIM, AssemblyType A, IntegrationType I,
          typename OpBase>
OpEssentialLhsImpl<DisplacementCubitBcData, BASE_DIM, FIELD_DIM, A, I,
                   OpBase>::OpEssentialLhsImpl(const std::string field_name,
                                               boost::shared_ptr<Range>
                                                   ents_ptr)
    : OpMass(
          field_name, field_name,

          [](double, double, double) constexpr { return 1; },

          ents_ptr) {}

// template <int BASE_DIM, int FIELD_DIM, AssemblyType A, IntegrationType I,
//           typename OpBase>
// struct AddEssentialToRhsPipelineImpl<DisplacementCubitBcData, BASE_DIM,
//                                      FIELD_DIM, A, I, OpBase> {

//   AddEssentialToRhsPipelineImpl() = delete;

//   static MoFEMErrorCodes add() {
//     MoFEMFunctionBegin;

//     using OP =
//         typename EssentialBC<OpBase>::template Assembly<A>::template LinearForm<
//             I>::template OpEssentialRhsImpl<DisplacementCubitBcData, BASE_DIM,
//                                             FIELD_DIM>;

//     auto add_op = [&](auto &&bcs) {
//       MoFEMFunctionBegin;
//       for (auto &m : bsc) {
//         if (auto bc = m.second->dispBcPtr) {
//           MOFEM_TAG_AND_LOG("SELF", sev, "OpEssentialRhs") << *bc;
//           pipeline.push_back(new OP(field_name, bc, m.second->getBcEdgesPtr()));
//         }
//       }
//       MOFEM_LOG_CHANNEL("SELF");
//       MoFEMFunctionReturn(0);
//     };

//     CHKERR add_op(m_field.getInterface<BcManager>()->getBcMapByBlockName());

//     MoFEMFunctionReturn(0);
//   };
// };

} // namespace MoFEM

#endif //_ESSENTIAL_BC_