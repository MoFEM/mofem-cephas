/**
 * @file NaturalTemperatureMeshsets.hpp
 * @brief
 * @version 0.1
 * @date 2022-09-18
 *
 * @copyright Copyright (c) 2022
 *
 */

#ifndef _NATURAL_TEMPERATURE_MESHSETS_HPP_
#define _NATURAL_TEMPERATURE_MESHSETS_HPP_

namespace MoFEM {

/**
 * @brief Type generating specialisation for temperature meshsets
 *
 */
struct NaturalTemperatureMeshsets {};

template <int FIELD_DIM, AssemblyType A, IntegrationType I, typename OpBase>
struct OpFluxRhsImpl<NaturalTemperatureMeshsets, 3, FIELD_DIM, A, I, OpBase>;

template <int BASE_DIM, int FIELD_DIM, AssemblyType A, IntegrationType I,
          typename OpBase>
struct AddFluxToRhsPipelineImpl<

    OpFluxRhsImpl<NaturalTemperatureMeshsets, BASE_DIM, FIELD_DIM, A, I,
                  OpBase>,
    A, I, OpBase

    > {

  AddFluxToRhsPipelineImpl() = delete;

  static MoFEMErrorCode add(

      boost::ptr_deque<ForcesAndSourcesCore::UserDataOperator> &pipeline,
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
        I>::template AddFluxToPipeline<OpFluxTempSet>::add(pipeline,
                                                  m_field, field_name, smv, block_name, sev);
    CHKERR
    NaturalBC<OpBase>::template Assembly<A>::template LinearForm<
        I>::template AddFluxToPipeline<OpFluxBlockset>::add(pipeline,
                                                   m_field, field_name, smv, block_name, sev);

    MoFEMFunctionReturn(0);
  }
};

} // namespace MoFEM

#endif // _NATURAL_TEMPERATURE_MESHSETS_HPP_