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


/**
 * @brief Class (Function) to enforce essential constrains
 *
 * Class (Function) to enforce essential constrains for DOFs which were
 * removed from the system
 *
 * @tparam T
 */
template <typename T> struct EssentialPreProc;

/**
 * @brief Class (Function) to enforce essential constrains on the left hand side
 * diagonal
 *
 * @tparam T
 */
template <typename T> struct EssentialPostProcLhs;

/**
 * @brief Class (Function) to enforce essential constrains on the right hand
 * side diagonal
 *
 * @tparam T
 */
template <typename T> struct EssentialPostProcRhs;

/**
 * @brief Class (Function) to calculate residual
 * side diagonal
 *
 * @tparam T
 */
template <typename T> struct EssentialPreProcReaction;

/**
 * @brief Enforce essential constrains on rhs.
 *
 * This class is used when constrains are enforced by least square method.
 *
 * @tparam T
 * @tparam BASE_DIM
 * @tparam FIELD_DIM
 * @tparam A
 * @tparam I
 * @tparam OpBase
 */
template <typename T, int BASE_DIM, int FIELD_DIM, AssemblyType A,
          IntegrationType I, typename OpBase>
struct OpEssentialRhsImpl;

/**
 * @brief Enforce essential constrains on lhs
 *
 * This class is used when constrains are enforced by least square method.
 *
 * @tparam T
 * @tparam BASE_DIM
 * @tparam FIELD_DIM
 * @tparam A
 * @tparam I
 * @tparam OpBase
 */
template <typename T, int BASE_DIM, int FIELD_DIM, AssemblyType A,
          IntegrationType I, typename OpBase>
struct OpEssentialLhsImpl;

/**
 * @brief Function (factory) for setting operators for rhs pipeline
 *
 * @tparam T
 * @tparam A
 * @tparam I
 * @tparam OpBase
 */
template <typename T, AssemblyType A, IntegrationType I, typename OpBase>
struct AddEssentialToRhsPipelineImpl;

/**
 * @brief Function (factory) for setting operators for lhs pipeline
 *
 * @tparam T
 * @tparam A
 * @tparam I
 * @tparam OpBase
 */
template <typename T, AssemblyType A, IntegrationType I, typename OpBase>
struct AddEssentialToLhsPipelineImpl;

/**
 * @brief Essential boundary conditions
 * @ingroup mofem_forms
 *
 * @tparam EleOp
 */
template <typename EleOp> struct EssentialBC {

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
      template <typename T>
      using AddEssentialToPipeline =
          AddEssentialToRhsPipelineImpl<T, A, I, EleOp>;
    };

    template <IntegrationType I> struct BiLinearForm {
      template <typename T, int BASE_DIM, int FIELD_DIM>
      using OpEssentialLhs =
          OpEssentialLhsImpl<T, BASE_DIM, FIELD_DIM, A, I, EleOp>;
      template <typename T>
      using AddEssentialToPipeline =
          AddEssentialToLhsPipelineImpl<T, A, I, EleOp>;
    };
  };
};

} // namespace MoFEM

#endif //_ESSENTIAL_BC_