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

/**
 * @brief Vector of time scaling methods
 *
 */
using VecOfTimeScalingMethods = std::vector<boost::shared_ptr<ScalingMethod>>;

/**
 * @brief Wrapper to generate natural b.c. specialisation based on operator type
 *
 * @tparam OP
 */
template <typename OP> struct FluxOpType {};

template <typename T, int BASE_DIM, int FIELD_DIM, AssemblyType A,
          IntegrationType I, typename OpBase>
struct OpFluxRhsImpl;

template <typename T, int BASE_DIM, int FIELD_DIM, AssemblyType A,
          IntegrationType I, typename OpBase>
struct OpFluxLhsImpl;

template <typename T, AssemblyType A, IntegrationType I, typename OpBase>
struct AddFluxToRhsPipelineImpl;

template <typename T, AssemblyType A, IntegrationType I, typename OpBase>
struct AddFluxToLhsPipelineImpl;

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
      using OpFlux = OpFluxRhsImpl<T, BASE_DIM, FIELD_DIM, A, I, EleOp>;

      template <typename T>
      using AddFluxToRhsPipeline = AddFluxToRhsPipelineImpl<T, A, I, EleOp>;

      template <typename T>
      static MoFEMErrorCode addFluxToRhsPipeline(
          FluxOpType<T>,
          boost::ptr_vector<ForcesAndSourcesCore::UserDataOperator> &pipeline,
          MoFEM::Interface &m_field, const std::string field_name,
          std::vector<boost::shared_ptr<ScalingMethod>> smv,
          const std::string block_name = "", Sev sev = Sev::verbose);
    };

    template <IntegrationType I> struct BiLinearForm {

      template <typename T, int BASE_DIM, int FIELD_DIM>
      using OpFlux = OpFluxLhsImpl<T, BASE_DIM, FIELD_DIM, A, I, EleOp>;

      template <typename T>
      using AddFluxToLhsPipeline = AddFluxToLhsPipelineImpl<T, A, I, EleOp>;

      template <typename T>
      static MoFEMErrorCode addFluxToLhsPipeline(
          FluxOpType<T>,
          boost::ptr_vector<ForcesAndSourcesCore::UserDataOperator> &pipeline,
          MoFEM::Interface &m_field, const std::string field_name,
          std::vector<boost::shared_ptr<ScalingMethod>> smv,
          const std::string block_name = "", Sev sev = Sev::verbose);
    };

  }; // Assembly
};

template <typename OpBase>
template <AssemblyType A>
template <IntegrationType I>
template <typename T>
MoFEMErrorCode
NaturalBC<OpBase>::Assembly<A>::LinearForm<I>::addFluxToRhsPipeline(

    FluxOpType<T>,

    boost::ptr_vector<ForcesAndSourcesCore::UserDataOperator> &pipeline,
    MoFEM::Interface &m_field, const std::string field_name,
    std::vector<boost::shared_ptr<ScalingMethod>> smv,
    const std::string block_name, Sev sev

) {
  using ADD =
      NaturalBC<OpBase>::Assembly<A>::LinearForm<I>::AddFluxToRhsPipeline<T>;
  return ADD::add(FluxOpType<T>(), pipeline, m_field, field_name, smv,
                  block_name, sev);
}

template <typename OpBase>
template <AssemblyType A>
template <IntegrationType I>
template <typename T>
MoFEMErrorCode
NaturalBC<OpBase>::Assembly<A>::BiLinearForm<I>::addFluxToLhsPipeline(

    FluxOpType<T>,

    boost::ptr_vector<ForcesAndSourcesCore::UserDataOperator> &pipeline,
    MoFEM::Interface &m_field, const std::string field_name,
    std::vector<boost::shared_ptr<ScalingMethod>> smv,
    const std::string block_name, Sev sev

) {
  using ADD =
      NaturalBC<OpBase>::Assembly<A>::BiLinearForm<I>::AddFluxToLhsPipeline<T>;
  return ADD::add(FluxOpType<T>(), pipeline, m_field, field_name, smv,
                  block_name, sev);
}

} // namespace MoFEM

#endif //_NATURAL_BC_