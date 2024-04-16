/** \file BiLinearFormsIntegrators.hpp
  * \brief Bilinear forms integrators
  * \ingroup mofem_form

  \todo SSome operators could be optimised. To do that, we need to write tests
  and use Valgrind to profile code, shaking cache misses. For example, some
  operators should have iteration first over columns, then rows. ome operators.
  Since those operators are used in many problems, an implementation must be
  efficient.

*/

#ifndef __BILINEAR_FORMS_INTEGRATORS_HPP__
#define __BILINEAR_FORMS_INTEGRATORS_HPP__

#include <BiLinearFormsIntegratorsImpl.hpp>

namespace MoFEM {

/**
 * @brief Bilinear integrator form
 * @ingroup mofem_forms
 *
 * @tparam EleOp
 * @tparam A
 * @tparam I
 */
template <typename EleOp>
template <AssemblyType A>
template <IntegrationType I>
struct FormsIntegrators<EleOp>::Assembly<A>::BiLinearForm {

  /**
   * @brief Integrate \f$(v_{,i},\beta(\mathbf{x}) u_{,j}))_\Omega\f$
   * @ingroup mofem_forms
   *
   * @tparam SPACE_DIM
   */
  template <int BASE_DIM, int FIELD_DIM, int SPACE_DIM>
  using OpGradGrad = OpGradGradImpl<BASE_DIM, FIELD_DIM, SPACE_DIM, I, OpBase>;

  /**
   * @brief Integrate \f$(v_i,\beta(\mathbf{x}) u_j)_\Omega\f$
   * @ingroup mofem_forms
   *
   * @note That FIELD_DIM = 4 or 9 is assumed that OpMass is for tensorial field
   * 2x2 or 3x3
   *
   * @todo It should be considered another template parameter RANK which will
   * allow to distinguish between scalar, vectorial and tensorial fields
   *
   * @tparam BASE_DIM dimension of base
   * @tparam FIELD_DIM dimension of field
   */
  template <int BASE_DIM, int FIELD_DIM>
  using OpMass = OpMassImpl<BASE_DIM, FIELD_DIM, I, OpBase>;

  /**
   * @brief Integrate \f$(v_k,D_{ijkl} u_{,l})_\Omega\f$
   *
   * \note \f$D_{ijkl}\f$ is * tensor with minor & major symmetry
   *
   * @ingroup mofem_forms
   *
   * @tparam SPACE_DIM
   * @tparam S
   */
  template <int BASE_DIM, int FIELD_DIM, int SPACE_DIM, int S = 0>
  using OpGradSymTensorGrad =
      OpGradSymTensorGradImpl<BASE_DIM, FIELD_DIM, SPACE_DIM, S, I, OpBase>;

  /**
   * @brief Integrate \f$(v_{,ij},D_{ijkl} u_{,kl})_\Omega\f$
   *
   * \note \f$D_{ijkl}\f$ is * tensor with no symmetries
   *
   * @ingroup mofem_forms
   *
   * @tparam SPACE_DIM
   * @tparam S
   */
  template <int BASE_DIM, int FIELD_DIM, int SPACE_DIM, int S = 0>
  using OpGradGradSymTensorGradGrad =
      OpGradGradSymTensorGradGradImpl<BASE_DIM, FIELD_DIM, SPACE_DIM, S, I,
                                      OpBase>;

  /**
   * @brief Integrate \f$(v_{,j},D_{ijkl} u_{,l})_\Omega\f$
   *
   * \note \f$D_{ijkl}\f$ is * tensor with no symmetries
   *
   * @ingroup mofem_forms
   *
   * @tparam SPACE_DIM
   * @tparam S
   */
  template <int BASE_DIM, int FIELD_DIM, int SPACE_DIM, int S = 0>
  using OpGradTensorGrad =
      OpGradTensorGradImpl<BASE_DIM, FIELD_DIM, SPACE_DIM, S, I, OpBase>;

  /**
   * @brief Integrate \f$(\lambda_{i,i},u)_\Omega\f$
   *
   * @tparam SPACE_DIM
   */
  template <int SPACE_DIM>
  using OpMixDivTimesScalar = OpMixDivTimesScalarImpl<SPACE_DIM, I, OpBase>;

  /**
   * @brief Integrate \f$(\lambda_{ij,j},u_{i})_\Omega\f$
   *
   * @tparam SPACE_DIM
   */
  template <int SPACE_DIM, CoordinateTypes CoordSys = CARTESIAN>
  using OpMixDivTimesVec = OpMixDivTimesVecImpl<SPACE_DIM, I, OpBase, CoordSys>;

  /**
   * @brief Integrate \f$(\lambda,u_{i,i})_\Omega\f$
   *
   * @tparam SPACE_DIM
   */
  template <int SPACE_DIM, CoordinateTypes COORDINATE_SYSTEM = CARTESIAN>
  using OpMixScalarTimesDiv =
      OpMixScalarTimesDivImpl<SPACE_DIM, I, OpBase, COORDINATE_SYSTEM>;

  /**
   * @brief Integrate \f$(\lambda_{i},u_{,j})_\Omega\f$
   *
   * @tparam SPACE_DIM
   */
  template <int BASE_DIM, int FIELD_DIM, int SPACE_DIM>
  using OpMixVectorTimesGrad =
      OpMixVectorTimesGradImpl<BASE_DIM, FIELD_DIM, SPACE_DIM, I, OpBase>;

  /**
   * @brief Integrate \f$(\lambda_{ij},u_{i,j})_\Omega\f$
   *
   * @tparam SPACE_DIM
   */
  template <int SPACE_DIM>
  using OpMixTensorTimesGrad = OpMixTensorTimesGradImpl<SPACE_DIM, I, OpBase>;
};

} // namespace MoFEM

#endif //__BILINEAR_FORMS_INTEGRATORS_HPP__