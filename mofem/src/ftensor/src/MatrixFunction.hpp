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
 auto  t_A = EigenProjection<double, double, 3>::getMat(t_L, t_N, f);
 \endcode
 where <3> means that are three unique eigen values. Return t_A is symmetric
 tensor rank two.

 Calculate directive
 \code
 auto t_P = EigenProjection<double, double, 3>::getDiffMat(t_L, t_N, f, d_f);
 \endcode
 where return t_SL is 4th order tensor (symmetry on first two and
 second to indices, i.e. minor symmetrise)

 Calculate second derivative, L, such that S:L, for given S,
 \code
 FTensor::Tensor2<double, 3, 3> t_S{

    1., 0., 0.,

    0., 1., 0.,

    0., 0., 1.};

  auto t_SL = EigenProjection<double, double, 3>::getDiffDiffMat(
                  t_L, t_N, f, d_f, dd_f, t_S)
 \endcode
 where return t_SL is 4th order tensor (symmetry on first two and
 second to indices, i.e. minor symmetrise)

 You can calculate eigen values using lapack.

 *
 */

/* This file is part of MoFEM.
 * MoFEM is free software: you can redistribute it and/or modify it under
 * the terms of the GNU Lesser General Public License as published by the
 * Free Software Foundation, either version 3 of the License, or (at your
 * option) any later version.
 *
 * MoFEM is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
 * License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with MoFEM. If not, see <http://www.gnu.org/licenses/>. */

#pragma once

namespace EigenMatrix {

template <typename T> using Val = const FTensor::Tensor1<T, 3>;
template <typename T> using Vec = const FTensor::Tensor2<T, 3, 3>;
template <typename T> using Fun = boost::function<T(const T)>;

FTensor::Tensor2_symmetric<double, 3> getMat(Val<double> &t_val,
                                             Vec<double> &t_vec, Fun<double> f,
                                             const int nb);

FTensor::Ddg<double, 3, 3> getDiffMat(Val<double> &t_val, Vec<double> &t_vec,
                                      Fun<double> f, Fun<double> d_f,
                                      const int nb);

FTensor::Ddg<double, 3, 3> getDiffDiffMat(Val<double> &t_val,
                                          Vec<double> &t_vec, Fun<double> f,
                                          Fun<double> d_f, Fun<double> dd_f,
                                          FTensor::Tensor2<double, 3, 3> &t_S,
                                          const int nb);


} // namespace EigenMatrix