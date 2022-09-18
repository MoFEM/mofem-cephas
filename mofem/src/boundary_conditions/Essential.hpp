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
 * @brief Wrapper on user dtat (element) operator used to select specialization
 * of essential bc.
 *
 * @tparam T
 */
template <typename T> struct EssentialOpType {};

/**
 * @brief Specialisation for b.c. applied by different types of meshsets
 *
 * @tparam BC
 */
template <CubitBC BC> struct EssentialMeshsetType {};

/**
 * @brief Class (Function) to enforce essential constrains
 *
 * Class (Function) to enforce essential constrains for DOFs which were
 * removed from the system
 *
 * @tparam T
 */
template <typename T> struct EssentialPreProc {};

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

      template <typename T>
      static MoFEMErrorCode addEssentialToRhsPipeline(
          EssentialOpType<T>, MoFEM::Interface &m_field,
          boost::ptr_vector<ForcesAndSourcesCore::UserDataOperator> &pipeline,
          std::string problem_name, std::string field_name,
          boost::shared_ptr<MatrixDouble> field_mat_ptr,
          std::vector<boost::shared_ptr<ScalingMethod>> smv);
    };

    template <IntegrationType I> struct BiLinearForm {
      template <typename T, int BASE_DIM, int FIELD_DIM>
      using OpEssentialLhs =
          OpEssentialLhsImpl<T, BASE_DIM, FIELD_DIM, A, I, EleOp>;
      template <typename T>
      static MoFEMErrorCode addEssentialToLhsPipeline(
          EssentialOpType<T>, MoFEM::Interface &m_field,
          boost::ptr_vector<ForcesAndSourcesCore::UserDataOperator> &pipeline,
          std::string problem_name, std::string field_name);
    };
  };
};


template <typename OpBase>
template <AssemblyType A>
template <IntegrationType I>
template <typename T>
MoFEMErrorCode
EssentialBC<OpBase>::Assembly<A>::LinearForm<I>::addEssentialToRhsPipeline(
    EssentialOpType<T>, MoFEM::Interface &m_field,
    boost::ptr_vector<ForcesAndSourcesCore::UserDataOperator> &pipeline,
    std::string problem_name, std::string field_name,
    boost::shared_ptr<MatrixDouble> field_mat_ptr,
    std::vector<boost::shared_ptr<ScalingMethod>> smv) {
  return AddEssentialToRhsPipelineImpl<T, A, I, OpBase>::add(
      EssentialOpType<T>(), m_field, pipeline, problem_name, field_name,
      field_mat_ptr, smv);
}

template <typename OpBase>
template <AssemblyType A>
template <IntegrationType I>
template <typename T>
MoFEMErrorCode
EssentialBC<OpBase>::Assembly<A>::BiLinearForm<I>::addEssentialToLhsPipeline(
    EssentialOpType<T>, MoFEM::Interface &m_field,
    boost::ptr_vector<ForcesAndSourcesCore::UserDataOperator> &pipeline,
    std::string problem_name, std::string field_name) {
  return AddEssentialToLhsPipelineImpl<T, A, I, OpBase>::add(
      EssentialOpType<T>(), m_field, pipeline, problem_name, field_name);
}

} // namespace MoFEM

#endif //_ESSENTIAL_BC_