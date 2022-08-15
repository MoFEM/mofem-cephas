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

template <int BASE_DIM, int FIELD_DIM, IntegrationType I, typename OpBase>
struct OpTimeScaledImpl;

/**
 * @brief Integrate Neumann boundary condition or body forces scaled by time
 *
 */
template <int BASE_DIM, int FIELD_DIM, IntegrationType I, typename OpBase>
struct OpTimeScaledImpl<BASE_DIM, FIELD_DIM, GAUSS, OpBase>
    : public FormsIntegrators<OpBase>::Assembly<I>::LinearForm<I>::OpSource<
          BASE_DIM, FIELD_DIM> {
  /**
   * @brief Construct a new Op Source Impl object
   *
   * @param field_name
   * @param time_fun
   * @param source_fun
   * @param ents_ptr
   */
  OpTimeScaledImpl(const std::string field_name, TimeFun time_fun,
                   ScalarFun source_fun, std::string block_set_name);

protected:



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

    using OpBase = OpBaseImpl<A, EleOp>;

    /**
     * @brief Linear form
     * @ingroup mofem_forms
     *
     * @tparam I
     */
    template <IntegrationType I> struct LinearForm {

      template <int BASE_DIM, int FIELD_DIM>
      struct OpTimeScaledSource : public OpTimeScaledSourceImpl {
        using OpTimeScaledSourceImpl<BASE_DIM, FIELD_DIM>::TimeScaledSourceImpl;
      };
    };

  }; // Assembly
}

};

#endif //_NATURAL_BC_