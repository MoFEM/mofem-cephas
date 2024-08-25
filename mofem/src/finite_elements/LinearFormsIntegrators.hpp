/** \file LinearFormsIntegrators.hpp
  * \brief Linear forms integrators
  * \ingroup mofem_form

*/

#ifndef __LINEAR_FORMS_INTEGRATORS_HPP__
#define __LINEAR_FORMS_INTEGRATORS_HPP__



namespace MoFEM {

/**
 * @brief Linear integrator form
 * @ingroup mofem_forms
 *
 * @tparam EleOp
 * @tparam A
 * @tparam I
 */
template <typename EleOp>
template <AssemblyType A>
template <IntegrationType I>
struct FormsIntegrators<EleOp>::Assembly<A>::LinearForm {

  /**
   * @brief Integrate \f$(v,f(\mathbf{x}))_\Omega\f$, f is a scalar
   * @ingroup mofem_forms
   *
   * @note \f$f(\mathbf{x})\$ is scalar function of coordinates
   *
   * @tparam BASE_DIM
   * @tparam FIELD_DIM
   */
  template <int BASE_DIM, int FIELD_DIM,
            typename S = SourceFunctionSpecialization>
  using OpSource =
      OpSourceImpl<BASE_DIM, FIELD_DIM, I, typename S::template S<OpBase>>;

  /**
   * @brief Vector field integrator \f$(v_i,f)_\Omega\f$, f is a vector
   * @ingroup mofem_forms
   *
   * @tparam BASE_DIM
   */
  template <int BASE_DIM, int S = 1>
  using OpBaseTimesScalar = OpBaseTimesScalarImpl<BASE_DIM, S, I, OpBase>;

  /** @deprecated use instead OpBaseTimesScalar
   */
  template <int BASE_DIM, int S = 1>
  using OpBaseTimesScalarField = OpBaseTimesScalar<BASE_DIM, S>;

  /**
   * @brief Vector field integrator \f$(v,f_i)_\Omega\f$, f is a vector
   * @ingroup mofem_forms
   *
   * @tparam BASE_DIM
   * @tparam FIELD_DIM
   * @tparam 0
   */
  template <int BASE_DIM, int FIELD_DIM, int S>
  using OpBaseTimesVector =
      OpBaseTimesVectorImpl<BASE_DIM, FIELD_DIM, S, I, OpBase>;
  //! [Source operator]

  //! [Grad times tensor]

  /**
   * @brief Integrate \f$(v_{,i},f_{ij})_\Omega\f$, f is a vector
   * @ingroup mofem_forms
   *
   * @note \f$f_{ij}\$ is tensor at integration points
   *
   * @tparam BASE_DIM
   * @tparam FIELD_DIM
   * @tparam SPACE_DIM
   */
  template <int BASE_DIM, int FIELD_DIM, int SPACE_DIM, int S = 1>
  using OpGradTimesTensor =
      OpGradTimesTensorImpl<BASE_DIM, FIELD_DIM, SPACE_DIM, S, I, OpBase>;

  /**
   * @brief Integrate \f$(v_{,i},f_{ij})_\Omega\f$, f is symmetric tensor
   * @ingroup mofem_forms
   *
   * @note \f$f_{ij}\$ is tensor at integration points
   *
   * @tparam BASE_DIM
   * @tparam FIELD_DIM
   * @tparam SPACE_DIM
   */
  template <int BASE_DIM, int FIELD_DIM, int SPACE_DIM, int S = 1>
  using OpGradTimesSymTensor =
      OpGradTimesSymTensorImpl<BASE_DIM, FIELD_DIM, SPACE_DIM, S, I, OpBase>;

  /**
   * @brief Integrate \f$(\lambda_{ij,j},u_{i})_\Omega\f$
   *
   * @tparam BASE_DIM
   * @tparam FIELD_DIM
   * @tparam SPACE_DIM
   */
  template <int BASE_DIM, int FIELD_DIM, int SPACE_DIM,
            CoordinateTypes CoordSys = CARTESIAN>
  using OpMixDivTimesU =
      OpMixDivTimesUImpl<BASE_DIM, FIELD_DIM, SPACE_DIM, I, OpBase, CoordSys>;

  /**
   * @brief Integrate \f$(\lambda_{ij},u_{i,j})_\Omega\f$
   *
   * @tparam SPACE_DIM
   */
  template <int SPACE_DIM>
  using OpMixTensorTimesGradU = OpMixTensorTimesGradUImpl<SPACE_DIM, I, OpBase>;

  /**
   * @brief Integrate \f$(u_{i},\lambda_{ij,j})_\Omega\f$
   *
   * @tparam SPACE_DIM
   */
  template <int SPACE_DIM>
  using OpMixVecTimesDivLambda =
      OpMixVecTimesDivLambdaImpl<SPACE_DIM, I, OpBase>;

  /**
   * @brief Multiply vector times normal on the face times scalar function
   *
   * This operator typically will be used to evaluate natural boundary
   * conditions for mixed formulation.
   *
   * @tparam BASE_DIM
   * @tparam SPACE_DIM
   * @tparam OpBase
   */
  template <int SPACE_DIM>
  using OpNormalMixVecTimesScalar =
      OpNormalMixVecTimesScalarImpl<SPACE_DIM, I, OpBase>;

  /**
   * @brief Multiply vector times normal on the face times vector field
   *
   * This operator typically will be used to evaluate natural boundary
   * conditions for mixed formulation.
   *
   * @tparam BASE_DIM
   * @tparam SPACE_DIM
   * @tparam OpBase
   */
  template <int SPACE_DIM>
  using OpNormalMixVecTimesVectorField =
      OpNormalMixVecTimesVectorFieldImpl<SPACE_DIM, I, OpBase>;

  /**
   * @brief Convective term
   *
   * \f[
   * (v, u_i \mathbf{y}_{,i})
   * \f]
   * where \f$\mathbf{y}\$ can be scalar or vector.
   *
   * @tparam BASE_DIM
   * @tparam FIELD_DIM
   * @tparam SPACE_DIM
   */
  template <int BASE_DIM, int FIELD_DIM, int SPACE_DIM>
  using OpConvectiveTermRhs =
      OpConvectiveTermRhsImpl<BASE_DIM, FIELD_DIM, SPACE_DIM, I, OpBase>;

  /**
   * @brief Assemble Rhs for constraint matrix C while hybridisation.
   *
   * @tparam FIELD_DIM dimension of the field
   */
  template <int FIELD_DIM>
  using OpBrokenSpaceConstrainDHybrid =
      OpBrokenSpaceConstrainDHybridImpl<FIELD_DIM, I, OpBase>;

  /**
   * @brief Assemble Rhs for contraint matrix C^T whole hybridisation
   * 
   * @tparam FIELD_DIM 
   */
  template <int FIELD_DIM>
  using OpBrokenSpaceConstrainDFlux =
      OpBrokenSpaceConstrainDFluxImpl<FIELD_DIM, I, OpBrokenBase>;
};

} // namespace MoFEM

#endif // __LINEAR_FORMS_INTEGRATORS_HPP__