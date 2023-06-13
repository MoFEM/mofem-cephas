/** \file FunctionMatrix.hpp
 \brief Get function from matrix
 \ingroup ftensor

 For reference see \cite miehe2001algorithms

 Usage:

 To calculate exponent of matrix, first and second derivatives
 \code
 auto f = [](double v) { return exp(v); };
 auto d_f = [](double v) { return exp(v); };
 auto dd_f = [](double v) { return exp(v); };
 \endcode

 Calculate matrix here t_L are vector of eigen values, and t_N is matrix of
 eigen vectors.
 \code
 auto  t_A = EigenMatrix::getMat(t_L, t_N, f);
 \endcode
 Return t_A is symmetric tensor rank two.

 Calculate derivative
 \code
 auto t_P = EigenMatrix::getDiffMat(t_L, t_N, f, d_f ,nb);
 \endcode
 where return t_SL is 4th order tensor (symmetry on first two and
 second to indices, i.e. minor symmetrise)

 Calculate second derivative, L, such that S:L, for given S,
 \code
 FTensor::Tensor2<double, 3, 3> t_S{

    1., 0., 0.,

    0., 1., 0.,

    0., 0., 1.};

  auto t_SL = EigenMatrix::getDiffDiffMat( t_L, t_N, f, d_f, dd_f, t_S, nb)
 \endcode
 where return t_SL is 4th order tensor (symmetry on first two and
 second to indices, i.e. minor symmetrise)

 You can calculate eigen values using lapack.

 Eiegn values should be sorted such that unique values are first, and last eiegn
 value is reapitting. For example if eiegn values are \f[ \lambda = \{1,1,2} \f]
 should be sorted such that
 \f[
 \lambda = \{1,2,1}
 \f]
 Eigen vectors should be updated, such that order of eigen vectors follow order
 of eigen values.

 *
 */



#pragma once

namespace EigenMatrix {

template <typename T, int Dim> using Val = const FTensor::Tensor1<T, Dim>;
template <typename T, int Dim> using Vec = const FTensor::Tensor2<T, Dim, Dim>;
template <typename T> using Fun = boost::function<T(const T)>;

/**
 * @brief Get the Mat object
 *
 * \f[
 * \mathbf{B} = f(\mathbf{A})
 * \f]
 *
 * \f[
 * B_{ij} = \sum_{a}^d f(\lambda^a) n^a_i n^a_j
 * \f]
 * where \f$a\f$ is eigen value number.
 *
 * @param t_val eigen values
 * @param t_vec eigen vector
 * @param f function
 * @return FTensor::Tensor2_symmetric<double, 3>
 */
FTensor::Tensor2_symmetric<double, 3>
getMat(Val<double, 3> &t_val, Vec<double, 3> &t_vec, Fun<double> f);

/**
 * @copydoc EigenMatrix::getMat
 */
FTensor::Tensor2_symmetric<double, 3>
getMat(Val<FTensor::PackPtr<double *, 1>, 3> &t_val,
       Vec<FTensor::PackPtr<double *, 1>, 3> &t_vec, Fun<double> f);

/**
 * @brief Get the Diff Mat object
 *
 * \f[
 * P_{ijkl} = \frac{\partial B_{ij}}{\partial A_{kl}}
 * \f]
 *
 * \note Eiegn vetore are in rows.
 *
 * @param t_val eigen values
 * @param t_vec eigen vector
 * @param f function
 * @param d_f directive of function
 * @param nb number of nonequal eigen valuse
 * @return FTensor::Ddg<double, 3, 3>
 */
FTensor::Ddg<double, 3, 3> getDiffMat(Val<double, 3> &t_val,
                                      Vec<double, 3> &t_vec, Fun<double> f,
                                      Fun<double> d_f, const int nb);
/**
 * @copydoc EigenMatrix::getDiffMat
 */
FTensor::Ddg<double, 3, 3>
getDiffMat(Val<FTensor::PackPtr<double *, 1>, 3> &t_val,
           Vec<FTensor::PackPtr<double *, 1>, 3> &t_vec, Fun<double> f,
           Fun<double> d_f, const int nb);

/**
 * @brief Get the Diff Diff Mat object
 *
 * \f[
 * LS_{klmn} =
 * S_{ij} \frac{\partial^2 B_{ij}}{\partial A_{kl} \partial A_{mn} }
 * \f]
 *
 * \note Eiegn vetore are in rows.
 *
 * @param t_val eiegn values
 * @param t_vec eigen vectors
 * @param f function
 * @param d_f directive of function
 * @param dd_f second directive of function
 * @param t_S S tensor
 * @param nb number of nonzero eigen values
 * @return FTensor::Ddg<double, 3, 3>
 */
FTensor::Ddg<double, 3, 3>
getDiffDiffMat(Val<double, 3> &t_val, Vec<double, 3> &t_vec, Fun<double> f,
               Fun<double> d_f, Fun<double> dd_f,
               FTensor::Tensor2_symmetric<double, 3> &t_S, const int nb);

/**
 * @copydoc EigenMatrix::getDiffDiffMat
 */
FTensor::Ddg<double, 3, 3> getDiffDiffMat(Val<double, 3> &t_val,
                                          Vec<double, 3> &t_vec, Fun<double> f,
                                          Fun<double> d_f, Fun<double> dd_f,
                                          FTensor::Tensor2<double, 3, 3> &t_S,
                                          const int nb);
/**
 * @copydoc EigenMatrix::getDiffMat
 */
FTensor::Ddg<double, 3, 3>
getDiffDiffMat(Val<FTensor::PackPtr<double *, 1>, 3> &t_val,
               Vec<FTensor::PackPtr<double *, 1>, 3> &t_vec, Fun<double> f,
               Fun<double> d_f, Fun<double> dd_f,
               FTensor::Tensor2<double, 3, 3> &t_S, const int nb);

/**
 * @copydoc EigenMatrix::getMat
 */
FTensor::Tensor2_symmetric<double, 2>
getMat(Val<double, 2> &t_val, Vec<double, 2> &t_vec, Fun<double> f);

/**
 * @copydoc EigenMatrix::getDiffMat
 */
FTensor::Ddg<double, 2, 2> getDiffMat(Val<double, 2> &t_val,
                                      Vec<double, 2> &t_vec, Fun<double> f,
                                      Fun<double> d_f, const int nb);

/**
 * @copydoc EigenMatrix::getDiffDiffMat
 */
FTensor::Ddg<double, 2, 2> getDiffDiffMat(Val<double, 2> &t_val,
                                          Vec<double, 2> &t_vec, Fun<double> f,
                                          Fun<double> d_f, Fun<double> dd_f,
                                          FTensor::Tensor2<double, 2, 2> &t_S,
                                          const int nb);

/**
 * @copydoc EigenMatrix::getDiffDiffMat
 */
FTensor::Ddg<double, 2, 2>
getDiffDiffMat(Val<double, 2> &t_val, Vec<double, 2> &t_vec, Fun<double> f,
               Fun<double> d_f, Fun<double> dd_f,
               FTensor::Tensor2_symmetric<double, 2> &t_S, const int nb);

} // namespace EigenMatrix