/**
 * @file ThermalOps.hpp
 * @author Ross Williams (ross.williams@glasgow.ac.uk)
 * @brief Thermal operators agnostic to small/large deformations
 * @version 0.14.0
 * @date 2025-03-31
 */

#ifndef __THERMAL_OPS_HPP__
#define __THERMAL_OPS_HPP__

namespace ThermalOps {

// Templated on IS_LARGE_STRAINS to allow different implementations with
// different numbers of arguments at compile time
/**
 * @brief Integrate Lhs base of flux (1/k) base of flux (FLUX x FLUX)
 *
 */
template <int SPACE_DIM, bool IS_LARGE_STRAINS> struct OpHdivHdivImpl;

template <int SPACE_DIM>
struct OpHdivHdivImpl<SPACE_DIM, false>
    : public FormsIntegrators<DomainEleOp>::Assembly<PETSC>::BiLinearForm<
          GAUSS>::OpMass<3, SPACE_DIM> {
  OpHdivHdivImpl(const std::string row_field_name,
                 const std::string col_field_name,
                 ScalarFun resistance_function,
                 boost::shared_ptr<MatrixDouble> mat_Grad_Ptr,
                 boost::shared_ptr<Range> ents_ptr = nullptr)
      : FormsIntegrators<DomainEleOp>::Assembly<PETSC>::BiLinearForm<
            GAUSS>::OpMass<3, SPACE_DIM>(row_field_name, col_field_name,
                                         resistance_function, ents_ptr) {}
};

template <int SPACE_DIM>
struct OpHdivHdivImpl<SPACE_DIM, true>
    : public OpCalculateQdotQLhs_dQ<SPACE_DIM, GAUSS, AssemblyDomainEleOp> {
  OpHdivHdivImpl(const std::string row_field_name,
                 const std::string col_field_name,
                 ScalarFun resistance_function,
                 boost::shared_ptr<MatrixDouble> mat_Grad_Ptr,
                 boost::shared_ptr<Range> ents_ptr = nullptr)
      : OpCalculateQdotQLhs_dQ<SPACE_DIM, GAUSS, AssemblyDomainEleOp>(
            row_field_name, col_field_name, resistance_function, mat_Grad_Ptr,
            ents_ptr) {}
};

// Add alias to allow same implementation for small and large strains
using OpHdivHdiv = OpHdivHdivImpl<SPACE_DIM, IS_LARGE_STRAINS>;

/**
 * @brief Integrate Lhs div of base of flux times base of temperature (FLUX x
 * T) and transpose of it, i.e. (T x FLUX)
 *
 */
using OpHdivT = FormsIntegrators<DomainEleOp>::Assembly<PETSC>::BiLinearForm<
    GAUSS>::OpMixDivTimesScalar<SPACE_DIM>;

// Templated on IS_LARGE_STRAINS to allow different implementations with
// different numbers of arguments at compile time. An empty operator is called
// for IS_LARGE_STRAINS = false
/**
 * @brief Integrate Lhs of flux term coupled to displacement field
 *
 */
using OpHdivU = OpCalculateQdotQLhs_dU<SPACE_DIM, GAUSS, AssemblyDomainEleOp,
                                       IS_LARGE_STRAINS>;

/**
 * @brief Integrate Lhs base of temperature times (heat capacity) times base of
 * temperature (T x T)
 *
 */
using OpCapacity = FormsIntegrators<DomainEleOp>::Assembly<PETSC>::BiLinearForm<
    GAUSS>::OpMass<1, 1>;

// Templated on IS_LARGE_STRAINS to allow different implementations with
// different numbers of arguments at compile time
/**
 * @brief Integrating Rhs flux base (1/k) flux  (FLUX)
 */
template <int SPACE_DIM, bool IS_LARGE_STRAINS> struct OpHdivFluxImpl;

template <int SPACE_DIM>
struct OpHdivFluxImpl<SPACE_DIM, false>
    : public FormsIntegrators<DomainEleOp>::Assembly<PETSC>::LinearForm<
          GAUSS>::OpBaseTimesVector<3, SPACE_DIM, 1> {
  OpHdivFluxImpl(const std::string field_name,
                 boost::shared_ptr<MatrixDouble> vec,
                 ScalarFun resistance_function,
                 boost::shared_ptr<MatrixDouble> mat_Grad_Ptr,
                 boost::shared_ptr<Range> ents_ptr = nullptr)
      : FormsIntegrators<DomainEleOp>::Assembly<PETSC>::LinearForm<
            GAUSS>::OpBaseTimesVector<3, SPACE_DIM, 1>(field_name, vec,
                                                       resistance_function,
                                                       ents_ptr) {}
};
template <int SPACE_DIM>
struct OpHdivFluxImpl<SPACE_DIM, true>
    : public OpCalculateQdotQRhs<SPACE_DIM, GAUSS, AssemblyDomainEleOp> {
  OpHdivFluxImpl(const std::string field_name,
                 boost::shared_ptr<MatrixDouble> vec,
                 ScalarFun resistance_function,
                 boost::shared_ptr<MatrixDouble> mat_Grad_Ptr,
                 boost::shared_ptr<Range> ents_ptr = nullptr)
      : OpCalculateQdotQRhs<SPACE_DIM, GAUSS, AssemblyDomainEleOp>(
            field_name, vec, resistance_function, mat_Grad_Ptr, ents_ptr) {}
};

// Add alias to allow same implementation for small and large strains
using OpHdivFlux = OpHdivFluxImpl<SPACE_DIM, IS_LARGE_STRAINS>;

/**
 * @brief  Integrate Rhs div flux base times temperature (T)
 *
 */
using OpHDivTemp = FormsIntegrators<DomainEleOp>::Assembly<PETSC>::LinearForm<
    GAUSS>::OpMixDivTimesU<3, 1, SPACE_DIM>;

/**
 * @brief Integrate Rhs base of temperature time heat capacity times heat rate
 * (T)
 *
 */
using OpBaseDotT = FormsIntegrators<DomainEleOp>::Assembly<PETSC>::LinearForm<
    GAUSS>::OpBaseTimesScalar<1>;

/**
 * @brief Integrate Rhs base of temperature times divergence of flux (T)
 *
 */
using OpBaseDivFlux = OpBaseDotT;

} // namespace ThermalOps

#endif // __THERMAL_OPS_HPP__