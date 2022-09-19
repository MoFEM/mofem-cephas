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
      using AddFluxToPipeline = AddFluxToRhsPipelineImpl<T, A, I, EleOp>;
    };

    template <IntegrationType I> struct BiLinearForm {

      template <typename T, int BASE_DIM, int FIELD_DIM>
      using OpFlux = OpFluxLhsImpl<T, BASE_DIM, FIELD_DIM, A, I, EleOp>;
      template <typename T>
      using AddFluxToPipeline = AddFluxToLhsPipelineImpl<T, A, I, EleOp>;
    };

  }; // Assembly
};

} // namespace MoFEM

#endif //_NATURAL_BC_