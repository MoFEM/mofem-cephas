/**
 * @file NaturalForceMeshsets.hpp
 * @brief Specialization for NaturalForceMeshsets
 * @version 0.13.2
 * @date 2022-09-18
 *
 * @copyright Copyright (c) 2022
 *
 */

#ifndef _NATURAL_FORCE_MESHSETS_HPP_
#define _NATURAL_FORCE_MESHSETS_HPP_

namespace MoFEM {

/**
 * @brief Type generating specialisation for force meshsets
 *
 */
struct NaturalForceMeshsets {};

template <int FIELD_DIM, AssemblyType A, IntegrationType I, typename OpBase>
struct OpFluxRhsImpl<NaturalForceMeshsets, 1, FIELD_DIM, A, I, OpBase>;

template <int BASE_DIM, int FIELD_DIM, AssemblyType A, IntegrationType I,
          typename OpBase>
struct AddFluxToRhsPipelineImpl<

    OpFluxRhsImpl<NaturalForceMeshsets, BASE_DIM, FIELD_DIM, A, I, OpBase>, A,
    I, OpBase

    > {

  AddFluxToRhsPipelineImpl() = delete;

  static MoFEMErrorCode add(

      boost::ptr_deque<ForcesAndSourcesCore::UserDataOperator> &pipeline,
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
        I>::template AddFluxToPipeline<OpFluxForceset>::add(pipeline, m_field,
                                                            field_name, smv,
                                                            block_name, sev);
    CHKERR
    NaturalBC<OpBase>::template Assembly<A>::template LinearForm<
        I>::template AddFluxToPipeline<OpFluxPressureset>::add(pipeline,
                                                               m_field,
                                                               field_name, smv,
                                                               block_name, sev);
    CHKERR
    NaturalBC<OpBase>::template Assembly<A>::template LinearForm<
        I>::template AddFluxToPipeline<OpFluxBlockset>::add(pipeline, m_field,
                                                            field_name, smv,
                                                            block_name, sev);

    MoFEMFunctionReturn(0);
  }
};

/**
 * @brief Type generating specialisation for force meshsets
 *
 */
struct NaturalForceMeshsetsScalarAndVectorScaling {};

template <int FIELD_DIM, AssemblyType A, IntegrationType I, typename OpBase>
struct OpFluxRhsImpl<NaturalForceMeshsetsScalarAndVectorScaling, 1, FIELD_DIM,
                     A, I, OpBase>;

template <int BASE_DIM, int FIELD_DIM, AssemblyType A, IntegrationType I,
          typename OpBase>
struct AddFluxToRhsPipelineImpl<

    OpFluxRhsImpl<NaturalForceMeshsetsScalarAndVectorScaling, BASE_DIM,
                  FIELD_DIM, A, I, OpBase>,
    A, I, OpBase

    > {

  AddFluxToRhsPipelineImpl() = delete;

  static MoFEMErrorCode add(

      boost::ptr_deque<ForcesAndSourcesCore::UserDataOperator> &pipeline,
      MoFEM::Interface &m_field, const std::string field_name,
      std::vector<boost::shared_ptr<ScalingMethod>> smv,
      const std::string block_name, Sev sev

  ) {
    return add(pipeline, m_field, field_name, smv, block_name, sev);
  }

  static MoFEMErrorCode add(

      boost::ptr_deque<ForcesAndSourcesCore::UserDataOperator> &pipeline,
      MoFEM::Interface &m_field, const std::string field_name,
      std::vector<boost::shared_ptr<ScalingMethod>> smv,
      std::vector<boost::shared_ptr<TimeScaleVector<FIELD_DIM>>> vsmv,
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

    using OpFluxBlocksetV =
        typename NaturalBC<OpBase>::template Assembly<A>::template LinearForm<
            I>::template OpFlux<NaturalMeshsetTypeVectorScaling<BLOCKSET>,
                                BASE_DIM, FIELD_DIM>;

    CHKERR
    NaturalBC<OpBase>::template Assembly<A>::template LinearForm<
        I>::template AddFluxToPipeline<OpFluxForceset>::add(pipeline, m_field,
                                                            field_name, smv,
                                                            block_name, sev);
    CHKERR
    NaturalBC<OpBase>::template Assembly<A>::template LinearForm<
        I>::template AddFluxToPipeline<OpFluxPressureset>::add(pipeline,
                                                               m_field,
                                                               field_name, smv,
                                                               block_name, sev);
    CHKERR
    NaturalBC<OpBase>::template Assembly<A>::template LinearForm<
        I>::template AddFluxToPipeline<OpFluxBlockset>::add(pipeline, m_field,
                                                            field_name, smv,
                                                            block_name, sev);
    CHKERR
    NaturalBC<OpBase>::template Assembly<A>::template LinearForm<
        I>::template AddFluxToPipeline<OpFluxBlocksetV>::add(pipeline, m_field,
                                                             field_name, vsmv,
                                                             block_name, sev);

    MoFEMFunctionReturn(0);
  }
};

} // namespace MoFEM

#endif //_NATURAL_FORCE_MESHSETS_HPP_