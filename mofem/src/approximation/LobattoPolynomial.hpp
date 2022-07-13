/** \file LobattoPolynomial.hpp
\brief Implementation of Lobatto polynomial

*/



#ifndef __LOBATTOPOLYNOMIALS_HPP__
#define __LOBATTOPOLYNOMIALS_HPP__

namespace MoFEM {

/**
 * \brief Class used to give arguments to Lobatto base functions
 * \ingroup mofem_base_functions
 */
struct LobattoPolynomialCtx : public LegendrePolynomialCtx {

  MoFEMErrorCode query_interface(boost::typeindex::type_index type_index,
                                 UnknownInterface **iface) const;

  LobattoPolynomialCtx(int p, double *diff_s, int dim,
                       boost::shared_ptr<MatrixDouble> base_fun_ptr,
                       boost::shared_ptr<MatrixDouble> base_diff_fun_ptr)
      : LegendrePolynomialCtx(p, diff_s, dim, base_fun_ptr, base_diff_fun_ptr) {
    basePolynomialsType0 = Lobatto_polynomials;
  }
  ~LobattoPolynomialCtx() {}
};

/**
 * \brief Calculating Lobatto base functions
 * \ingroup mofem_base_functions
 */
struct LobattoPolynomial : public LegendrePolynomial {

  MoFEMErrorCode query_interface(boost::typeindex::type_index type_index,
                                 UnknownInterface **iface) const;

  LobattoPolynomial() {}
  ~LobattoPolynomial() {}

  MoFEMErrorCode getValue(MatrixDouble &pts,
                          boost::shared_ptr<BaseFunctionCtx> ctx_ptr);
};

/**
 * \brief Class used to give arguments to Kernel Lobatto base functions
 * \ingroup mofem_base_functions
 */
struct KernelLobattoPolynomialCtx : public LegendrePolynomialCtx {

  MoFEMErrorCode query_interface(boost::typeindex::type_index type_index,
                                 UnknownInterface **iface) const;

  KernelLobattoPolynomialCtx(int p, double *diff_s, int dim,
                             boost::shared_ptr<MatrixDouble> base_fun_ptr,
                             boost::shared_ptr<MatrixDouble> base_diff_fun_ptr)
      : LegendrePolynomialCtx(p, diff_s, dim, base_fun_ptr, base_diff_fun_ptr) {
    basePolynomialsType0 = LobattoKernel_polynomials;
  }
  ~KernelLobattoPolynomialCtx() {}
};

/**
 * \brief Calculating Lobatto base functions
 * \ingroup mofem_base_functions
 */
struct KernelLobattoPolynomial : public LegendrePolynomial {

  MoFEMErrorCode query_interface(boost::typeindex::type_index type_index,
                                 UnknownInterface **iface) const;

  KernelLobattoPolynomial() {}
  ~KernelLobattoPolynomial() {}

  MoFEMErrorCode getValue(MatrixDouble &pts,
                          boost::shared_ptr<BaseFunctionCtx> ctx_ptr);
};

} // namespace MoFEM

#endif
