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
#include <type_traits>

template <typename E, typename C> struct d2MCoefficients {
  using Val = typename E::Val;
  using Vec = typename E::Vec;
  using Fun = typename E::Fun;

  template <int N> using Number = FTensor::Number<N>;

  d2MCoefficients() = delete;
  ~d2MCoefficients() = delete;

  template <int a, int b>
  static inline auto get(Val &t_val, const Number<a> &, const Number<b> &,
                         const Number<-1> &, const Number<-1> &,
                         const Number<3> &, Fun f, Fun d_f) {
    return f(E::L(t_val, Number<a>())) * E::F(t_val, Number<a>(), Number<b>());
  }

  template <int a, int b>
  static inline auto get(Val &t_val, const Number<a> &, const Number<b> &,
                         const Number<-1> &, const Number<-1> &,
                         const Number<2> &, Fun f, Fun d_f) {
    if (a == 1 || b == 1)
      return get(t_val, Number<a>(), Number<b>(), Number<-1>(), Number<-1>(),
                 Number<3>(), f, d_f);
    else
      return d_f(E::L(t_val, Number<a>())) / static_cast<C>(2);
  }

  template <int a, int b>
  static inline auto get(Val &t_val, const Number<a> &, const Number<b> &,
                         const Number<-1> &, const Number<-1> &,
                         const Number<1>, Fun f, Fun d_f) {
    return d_f(E::L(t_val, Number<a>())) / static_cast<C>(2);
  }
};

template <typename E, typename C> struct dd4MCoefficientsType1 {
  using Val = typename E::Val;
  using Vec = typename E::Vec;
  using Fun = typename E::Fun;

  template <int N> using Number = FTensor::Number<N>;

  dd4MCoefficientsType1() = delete;
  ~dd4MCoefficientsType1() = delete;

  template <int a, int b, int c, int d>
  static inline auto get(Val &t_val, const Number<a> &, const Number<b> &,
                         const Number<c> &, const Number<d> &,
                         const Number<3> &, Fun f, Fun dd_f) {
    return f(E::L(t_val, Number<c>())) * E::F(t_val, Number<c>(), Number<d>()) *
           E::F(t_val, Number<a>(), Number<b>());
  }

  template <int a, int b, int c, int d>
  static inline auto get(Val &t_val, const Number<a> &, const Number<b> &,
                         const Number<c> &, const Number<d> &,
                         const Number<2> &, Fun f, Fun dd_f) {
    if (a == 1 || b == 1)
      return get(t_val, Number<a>(), Number<b>(), Number<3>(), f, dd_f);
    else
      return get(t_val, Number<a>(), Number<b>(), Number<1>(), f, dd_f);
  }

  template <int a, int b, int c, int d>
  static inline auto get(Val &t_val, const Number<a> &, const Number<b> &,
                         const Number<c> &, const Number<d> &, const Number<1>,
                         Fun f, Fun dd_f) {
    return dd_f(E::L(t_val, Number<c>())) / static_cast<C>(2);
  }
};

template <typename E, typename C> struct dd4MCoefficientsType2 {
  using Val = typename E::Val;
  using Vec = typename E::Vec;
  using Fun = typename E::Fun;

  template <int N> using Number = FTensor::Number<N>;

  dd4MCoefficientsType2() = delete;
  ~dd4MCoefficientsType2() = delete;

  template <int a, int b, int m, int n>
  static inline auto get(Val &t_val, Vec &t_vec, const Number<a> &,
                         const Number<b> &, const Number<m> &,
                         const Number<n> &, const Number<3> &, Fun f,
                         Fun dd_f) {
    return f(E::L(t_val, Number<a>())) * E::dFdN(t_val, t_vec, Number<a>(),
                                                 Number<b>(), Number<m>(),
                                                 Number<n>());
  }

  template <int a, int b, int m, int n>
  static inline auto get(Val &t_val, Vec &t_vec, const Number<a> &,
                         const Number<b> &, const Number<m> &,
                         const Number<n> &, const Number<2> &, Fun f,
                         Fun dd_f) {
    if (a == 1 || b == 1)
      return get(t_val, Number<a>(), Number<b>(), Number<m>(), Number<n>(),
                 Number<3>(), f, dd_f);
    else
      return get(t_val, Number<a>(), Number<b>(), Number<m>(), Number<n>(),
                 Number<1>(), f, dd_f);
  }

  template <int a, int b, int m, int n>
  static inline auto get(Val &t_val, Vec &t_vec, const Number<a> &,
                         const Number<b> &, const Number<m> &,
                         const Number<n> &, const Number<1>, Fun f, Fun dd_f) {
    if (a != b) {
      return

          (

              E::M(t_vec, Number<a>(), Number<m>(), Number<n>()) *
                  dd_f(E::L(t_val, Number<a>())) +

              E::M(t_vec, Number<b>(), Number<m>(), Number<n>()) *
                  dd_f(E::L(t_val, Number<b>()))

                  ) /
          static_cast<C>(2);
    }
    return static_cast<C>(0);
  }
};

template <typename E, typename C, typename G, int a, int c, int d, int i, int j,
          int k, int l>
struct d2MImpl {
  using Val = typename E::Val;
  using Vec = typename E::Vec;
  using Fun = typename E::Fun;

  template <int N> using Number = FTensor::Number<N>;

  d2MImpl() = delete;
  ~d2MImpl() = delete;

  template <int b>
  static inline C term(Val &t_val, Vec &t_vec, Fun f, Fun d_f) {
    if (a != b) {
      return G::get(t_val, Number<a>(), Number<b>(), Number<c>(), Number<d>(),
                    typename E::NumberNb(), f, d_f) *
             E::S(t_vec, Number<a>(), Number<b>(), Number<i>(), Number<j>(),
                  Number<k>(), Number<l>());
    }
    return 0;
  }

  template <int nb>
  static inline C eval(Val &t_val, Vec &t_vec, Fun f, Fun d_f,
                       const Number<nb> &) {
    return term<nb - 1>(t_val, t_vec, f, d_f) +
           eval(t_val, t_vec, f, d_f, Number<nb - 1>());
  }

  static inline C eval(Val &t_val, Vec &t_vec, Fun f, Fun d_f,
                       const Number<1> &) {
    return term<0>(t_val, t_vec, f, d_f);
  }
};

template <typename E, typename C, typename G1, typename G2, int a, int i, int j,
          int k, int l, int m, int n>
struct dd4MImpl {
  using Val = typename E::Val;
  using Vec = typename E::Vec;
  using Fun = typename E::Fun;
  dd4MImpl() = delete;
  ~dd4MImpl() = delete;

  template <int N> using Number = FTensor::Number<N>;

  template <int b, int A, int I, int J, int K, int L>
  static inline auto d2M(Val &t_val, Vec &t_vec, Fun f, Fun dd_f) {
    using V =
        typename FTensor::promote<decltype(t_val(0)), decltype(t_vec(0, 0))>::V;
    return d2MImpl<E, V, G1, A, a, b, I, J, K, L>::eval(
        t_val, t_vec, f, dd_f, typename E::NumberDim());
  }

  template <int b, int A, int B, int I, int J, int K, int L, int M, int N>
  static inline auto d2G(Val &t_val, Vec &t_vec, Fun f, Fun dd_f) {
    return d2M<b, A, I, K, N, M>(t_val, t_vec, f, dd_f) *
               E::M(t_vec, Number<B>(), Number<J>(), Number<L>()) +
           E::M(t_vec, Number<A>(), Number<I>(), Number<K>()) *
               d2M<b, B, J, L, M, N>(t_val, t_vec, f, dd_f) +
           d2M<b, A, I, L, M, N>(t_val, t_vec, f, dd_f) *
               E::M(t_vec, Number<B>(), Number<J>(), Number<K>()) +
           E::M(t_vec, Number<A>(), Number<I>(), Number<L>()) *
               d2M<b, B, J, K, M, N>(t_val, t_vec, f, dd_f);
  }

  template <int A, int B, int I, int J, int K, int L, int M, int N>
  static inline auto d2S(Val &t_val, Vec &t_vec, Fun f, Fun dd_f) {
    return d2G<B, A, B, I, J, K, L, M, N>(t_val, t_vec, f, dd_f) +
           d2G<B, B, A, I, J, K, L, M, N>(t_val, t_vec, f, dd_f);
  }

  template <int b>
  static inline C term(Val &t_val, Vec &t_vec, Fun f, Fun d_f, Fun dd_f) {

    if (a != b) {
      return

          d2S<a, b, i, j, k, l, m, n>(t_val, t_vec, f, dd_f)

          +

          2 *

              G2::get(t_val, t_vec, Number<a>(), Number<b>(), Number<m>(),
                      Number<n>(), typename E::NumberNb(), f, dd_f) *
              E::S(t_vec, Number<a>(), Number<b>(), Number<i>(), Number<j>(),
                   Number<k>(), Number<l>());
    }

    return 0;
  }

  template <int nb>
  static inline C eval(Val &t_val, Vec &t_vec, Fun f, Fun d_f, Fun dd_f,
                       const Number<nb> &) {
    return term<nb - 1>(t_val, t_vec, f, d_f, dd_f) +
           eval(t_val, t_vec, f, d_f, dd_f, Number<nb - 1>());
  }

  static inline C eval(Val &t_val, Vec &t_vec, Fun f, Fun d_f, Fun dd_f,
                       const Number<1> &) {
    return term<0>(t_val, t_vec, f, d_f, dd_f);
  }
};

template <typename E, typename C, int i, int j> struct reconstructMatImpl {
  using Val = typename E::Val;
  using Vec = typename E::Vec;
  using Fun = typename E::Fun;

  template <int N> using Number = FTensor::Number<N>;

  reconstructMatImpl() = delete;
  ~reconstructMatImpl() = delete;

  template <int a> static inline C term(Val &t_val, Vec &t_vec, Fun f) {
    return E::M(t_vec, Number<a>(), Number<i>(), Number<j>()) *
           f(E::L(t_val, Number<a>()));
  }

  template <int nb>
  static inline C eval(Val &t_val, Vec &t_vec, Fun f, const Number<nb> &) {
    return term<nb - 1>(t_val, t_vec, f) +
           eval(t_val, t_vec, f, Number<nb - 1>());
  }

  static inline C eval(Val &t_val, Vec &t_vec, Fun f, const Number<1> &) {
    return term<0>(t_val, t_vec, f);
  }
};

template <typename E, typename C, int i, int j, int k, int l>
struct firstMatrixDirectiveImpl {
  using Val = typename E::Val;
  using Vec = typename E::Vec;
  using Fun = typename E::Fun;

  template <int N> using Number = FTensor::Number<N>;

  firstMatrixDirectiveImpl() = delete;
  ~firstMatrixDirectiveImpl() = delete;

  template <int a>
  static inline C term(Val &t_val, Vec &t_vec, Fun f, Fun d_f) {
    return

        E::M(t_vec, Number<a>(), Number<i>(), Number<j>()) *
            E::M(t_vec, Number<a>(), Number<k>(), Number<l>()) *
            d_f(E::L(t_val, Number<a>()))

        +

        d2MImpl<E, C, d2MCoefficients<E, C>, a, -1, -1, i, j, k, l>::eval(
            t_val, t_vec, f, d_f, Number<3>()) /
            static_cast<C>(2);
  }

  template <int nb>
  static inline C eval(Val &t_val, Vec &t_vec, Fun f, Fun d_f,
                       const Number<nb> &) {
    return term<nb - 1>(t_val, t_vec, f, d_f) +
           eval(t_val, t_vec, f, d_f, Number<nb - 1>());
  }

  static inline C eval(Val &t_val, Vec &t_vec, Fun f, Fun d_f,
                       const Number<1> &) {
    return term<0>(t_val, t_vec, f, d_f);
  }
};

template <typename E, typename C, int i, int j, int k, int l, int m, int n>
struct secondMatrixDirectiveImpl {
  using Val = typename E::Val;
  using Vec = typename E::Vec;
  using Fun = typename E::Fun;

  template <int N> using Number = FTensor::Number<N>;

  secondMatrixDirectiveImpl() = delete;
  ~secondMatrixDirectiveImpl() = delete;

  template <int a>
  static inline C term(Val &t_val, Vec &t_vec, Fun f, Fun d_f, Fun dd_f) {

    return

        (

            d2MImpl<E, C, d2MCoefficients<E, C>, a, -1, -1, i, j, m, n>::eval(
                t_val, t_vec, d_f, dd_f, Number<3>()) *
                E::M(t_vec, Number<a>(), Number<k>(), Number<l>())

            +

            E::M(t_vec, Number<a>(), Number<i>(), Number<j>()) *
                d2MImpl<E, C, d2MCoefficients<E, C>, a, -1, -1, k, l, m,
                        n>::eval(t_val, t_vec, d_f, dd_f, Number<3>())

                ) /
            static_cast<C>(2) +

        E::M(t_vec, Number<a>(), Number<i>(), Number<j>()) *
            E::M(t_vec, Number<a>(), Number<k>(), Number<l>()) *
            E::M(t_vec, Number<a>(), Number<m>(), Number<n>()) *
            dd_f(E::L(t_val, Number<a>())) +

        dd4MImpl<E, C, dd4MCoefficientsType1<E, C>,
                 dd4MCoefficientsType2<E, C>, a, i, j, k, l, m,
                 n>::eval(t_val, t_vec, f, d_f, dd_f, Number<3>()) /
            static_cast<C>(4) +

        d2MImpl<E, C, d2MCoefficients<E, C>, a, -1, -1, i, j, k, l>::eval(
            t_val, t_vec, d_f, dd_f, Number<3>()) *
            E::M(t_vec, Number<a>(), Number<m>(), Number<n>()) /
            static_cast<C>(2);
  }

  template <int nb>
  static inline C eval(Val &t_val, Vec &t_vec, Fun f, Fun d_f, Fun dd_f,
                       const Number<nb> &) {
    return term<nb - 1>(t_val, t_vec, f, d_f, dd_f)

           +

           secondMatrixDirectiveImpl<E, C, i, j, k, l, m, n>::eval(
               t_val, t_vec, f, d_f, dd_f, Number<nb - 1>());
  }

  static inline C eval(Val &t_val, Vec &t_vec, Fun f, Fun d_f, Fun dd_f,
                       const Number<1> &) {
    return term<0>(t_val, t_vec, f, d_f, dd_f);
  }
};

template <typename E, typename C> struct getMatImpl {
  using Val = typename E::Val;
  using Vec = typename E::Vec;
  using Fun = typename E::Fun;

  template <int N> using Number = FTensor::Number<N>;

  getMatImpl() = delete;
  ~getMatImpl() = delete;

  template <typename T, int I, int J>
  static inline void set(Val &t_val, Vec &t_vec, Fun f, T &t_a,
                         const Number<I> &, const Number<J> &) {
    set(t_val, t_vec, f, t_a, Number<I>(), Number<J - 1>());
    t_a(I - 1, J - 1) = reconstructMatImpl<E, C, I - 1, J - 1>::eval(
        t_val, t_vec, f, typename E::NumberNb());
  }

  template <typename T, int I>
  static inline void set(Val &t_val, Vec &t_vec, Fun f, T &t_a,
                         const Number<I> &, const Number<0> &) {
    set(t_val, t_vec, f, t_a, Number<I - 1>(), Number<I - 1>());
  }

  template <typename T>
  static inline void set(Val &t_val, Vec &t_vec, Fun f, T &t_a,
                         const Number<0> &, const Number<0> &) {}
};

template <typename E, typename C> struct getDiffMatImpl {
  using Val = typename E::Val;
  using Vec = typename E::Vec;
  using Fun = typename E::Fun;

  template <int N> using Number = FTensor::Number<N>;

  getDiffMatImpl() = delete;
  ~getDiffMatImpl() = delete;

  template <typename T, int I, int J, int K, int L>
  static inline void set(Val &t_val, Vec &t_vec, Fun f, Fun d_f, T &t_a,
                         const Number<I> &, const Number<J> &,
                         const Number<K> &, const Number<L> &) {
    set(t_val, t_vec, f, d_f, t_a, Number<I>(), Number<J>(), Number<K>(),
        Number<L - 1>());
    t_a(I - 1, J - 1, K - 1, L - 1) =
        firstMatrixDirectiveImpl<E, C, I - 1, J - 1, K - 1, L - 1>::eval(
            t_val, t_vec, f, d_f, typename E::NumberDim());
  }

  template <typename T, int I, int J, int K>
  static inline void set(Val &t_val, Vec &t_vec, Fun f, Fun d_f, T &t_a,
                         const Number<I> &, const Number<J> &,
                         const Number<K> &, const Number<0> &) {
    set(t_val, t_vec, f, d_f, t_a, Number<I>(), Number<J>(), Number<K - 1>(),
        Number<K - 1>());
  }

  template <typename T, int I, int J>
  static inline void set(Val &t_val, Vec &t_vec, Fun f, Fun d_f, T &t_a,
                         const Number<I> &, const Number<J> &,
                         const Number<0> &, const Number<0> &) {
    set(t_val, t_vec, f, d_f, t_a, Number<I>(), Number<J - 1>(),
        typename E::NumberDim(), typename E::NumberDim());
  }

  template <typename T, int I, int K, int L>
  static inline void set(Val &t_val, Vec &t_vec, Fun f, Fun d_f, T &t_a,
                         const Number<I> &, const Number<0> &,
                         const Number<K> &, const Number<L> &) {
    set(t_val, t_vec, f, d_f, t_a, Number<I - 1>(), Number<I - 1>(),
        Number<K>(), Number<L>());
  }

  template <typename T, int K, int L>
  static inline void set(Val &t_val, Vec &t_vec, Fun f, Fun d_f, T &t_a,
                         const Number<0> &, const Number<0> &,
                         const Number<K> &, const Number<L> &) {}
};

template <typename E, typename C> struct getDiffDiffMatImpl {
  using Val = typename E::Val;
  using Vec = typename E::Vec;
  using Fun = typename E::Fun;

  template <int N> using Number = FTensor::Number<N>;

  getDiffDiffMatImpl() = delete;
  ~getDiffDiffMatImpl() = delete;

  template <typename T1, typename T2, int I, int J, int K, int L, int M, int N>
  static inline auto add(Val &t_val, Vec &t_vec, Fun f, Fun d_f, Fun dd_f,
                         T1 &t_s, T2 &t_a, const Number<I> &, const Number<J> &,
                         const Number<K> &, const Number<L> &,
                         const Number<M> &, const Number<N> &) {
    return t_s(M - 1, N - 1) *
               secondMatrixDirectiveImpl<E, C, M - 1, N - 1, I - 1, J - 1,
                                         K - 1,
                                         L - 1>::eval(t_val, t_vec, f, d_f,
                                                      dd_f,
                                                      typename E::NumberNb())

           +

           add(t_val, t_vec, f, d_f, dd_f, t_s, t_a, Number<I>(), Number<J>(),
               Number<K>(), Number<L>(), Number<M>(), Number<N - 1>());
  }

  template <typename T1, typename T2, int I, int J, int K, int L, int M>
  static inline auto add(Val &t_val, Vec &t_vec, Fun f, Fun d_f, Fun dd_f,
                         T1 &t_s, T2 &t_a, const Number<I> &, const Number<J> &,
                         const Number<K> &, const Number<L> &,
                         const Number<M> &, const Number<1> &) {
    return t_s(M - 1, 0) *
               secondMatrixDirectiveImpl<E, C, M - 1, 0, I - 1, J - 1, K - 1,
                                         L - 1>::eval(t_val, t_vec, f, d_f,
                                                      dd_f,
                                                      typename E::NumberNb())

           +

           add(t_val, t_vec, f, d_f, dd_f, t_s, t_a, Number<I>(), Number<J>(),
               Number<K>(), Number<L>(), Number<M - 1>(),
               typename E::NumberDim());
  }

  template <typename T1, typename T2, int I, int J, int K, int L>
  static inline auto add(Val &t_val, Vec &t_vec, Fun f, Fun d_f, Fun dd_f,
                         T1 &t_s, T2 &t_a, const Number<I> &, const Number<J> &,
                         const Number<K> &, const Number<L> &,
                         const Number<1> &, const Number<1> &) {
    return t_s(0, 0) *
           secondMatrixDirectiveImpl<E, C, 0, 0, I - 1, J - 1, K - 1,
                                     L - 1>::eval(t_val, t_vec, f, d_f, dd_f,
                                                  typename E::NumberNb());
  }

  template <typename T1, typename T2, int I, int J, int K, int L>
  static inline void set(Val &t_val, Vec &t_vec, Fun f, Fun d_f, Fun dd_f,
                         T1 &t_s, T2 &t_a, const Number<I> &, const Number<J> &,
                         const Number<K> &, const Number<L> &) {
    set(t_val, t_vec, f, d_f, dd_f, t_s, t_a, Number<I>(), Number<J>(),
        Number<K>(), Number<L - 1>());
    t_a(I - 1, J - 1, K - 1, L - 1) =
        add(t_val, t_vec, f, d_f, dd_f, t_s, t_a, Number<I>(), Number<J>(),
            Number<K>(), Number<L>(), typename E::NumberDim(),
            typename E::NumberDim());
  }

  template <typename T1, typename T2, int I, int J, int K>
  static inline void set(Val &t_val, Vec &t_vec, Fun f, Fun d_f, Fun dd_f,
                         T1 &t_s, T2 &t_a, const Number<I> &, const Number<J> &,
                         const Number<K> &, const Number<0> &) {
    set(t_val, t_vec, f, d_f, dd_f, t_s, t_a,

        Number<I>(), Number<J>(), Number<K - 1>(), Number<K - 1>());
  }

  template <typename T1, typename T2, int I, int J>
  static inline void set(Val &t_val, Vec &t_vec, Fun f, Fun d_f, Fun dd_f,
                         T1 &t_s, T2 &t_a, const Number<I> &, const Number<J> &,
                         const Number<0> &, const Number<0> &) {
    set(t_val, t_vec, f, d_f, dd_f, t_s, t_a,

        Number<I>(), Number<J - 1>(), typename E::NumberDim(),
        typename E::NumberDim());
  }

  template <typename T1, typename T2, int I, int K, int L>
  static inline void set(Val &t_val, Vec &t_vec, Fun f, Fun d_f, Fun dd_f,
                         T1 &t_s, T2 &t_a, const Number<I> &, const Number<0> &,
                         const Number<K> &, const Number<L> &) {
    set(t_val, t_vec, f, d_f, dd_f, t_s, t_a,

        Number<I - 1>(), Number<I - 1>(), Number<K>(), Number<L>());
  }

  template <typename T1, typename T2, int K, int L>
  static inline void set(Val &t_val, Vec &t_vec, Fun f, Fun d_f, Fun dd_f,
                         T1 &t_s, T2 &t_a, const Number<0> &, const Number<0> &,
                         const Number<K> &, const Number<L> &) {}
};

template <typename T1, typename T2, int NB, int Dim = 3>
struct EigenProjection {

  using Val = const FTensor::Tensor1<T1, Dim>;
  using Vec = const FTensor::Tensor2<T2, Dim, Dim>;
  using Fun = boost::function<double(const double)>;

  template <int N> using Number = FTensor::Number<N>;
  template <char c> using I = typename FTensor::Index<c, Dim>;

  using NumberNb = Number<NB>;
  using NumberDim = Number<Dim>;

  /**
   * @brief Get matrix
   *
   * \f[
   * \mathbf{B} = f(\mathbf{A})
   * \f]
   *
   * \f[
   * B_{ij} = sum_{a}^3 f(\lambda^a) n^a_i n^a_j
   * \f]
   * where \f$a\f$ is eigen value number.
   *
   * @param t_val eiegn values vector
   * @param t_vec eigen vectors matrix
   * @param f function
   * @return auto function symmetric tensor rank two
   */
  static inline auto getMat(Val &t_val, Vec &t_vec, Fun f) {
    using V =
        typename FTensor::promote<decltype(t_val(0)), decltype(t_vec(0, 0))>::V;
    FTensor::Tensor2_symmetric<typename std::remove_const<V>::type, Dim> t_A;
    getMatImpl<EigenProjection<T1, T2, NB, Dim>, V>::set(
        t_val, t_vec, f, t_A, Number<Dim>(), Number<Dim>());
    return t_A;
  }

  /**
   * @brief Get derivative of matrix
   *
   * \f[
   * P_{ijkl} = \frac{\partial B_{ij}}{\partial A_{kl}}
   * \f]
   *
   * @param t_val eiegn values vector
   * @param t_vec eiegn vectors matrix
   * @param f function
   * @param d_f directive of function
   * @return auto direvatives, forth order tensor with minor simetries
   */
  static inline auto getDiffMat(Val &t_val, Vec &t_vec, Fun f, Fun d_f) {
    using V =
        typename FTensor::promote<decltype(t_val(0)), decltype(t_vec(0, 0))>::V;
    FTensor::Ddg<V, Dim, Dim> t_diff_A;
    getDiffMatImpl<EigenProjection<T1, T2, NB, Dim>, V>::set(
        t_val, t_vec, f, d_f, t_diff_A, Number<Dim>(), Number<Dim>(),
        Number<Dim>(), Number<Dim>());
    return t_diff_A;
  }

  /**
   * @brief Get second direvarive of matrix
   *
   * \f[
   * LS_{klmn} =
   * S_{ij} \frac{\partial^2 B_{ij}}{\partial A_{kl} \partial A_{mn} }
   * \f]
   *
   * @tparam T
   * @param t_val eiegn values vector
   * @param t_vec eiegn vectors matrix
   * @param f function
   * @param d_f derivative of function
   * @param dd_f second derivative of function
   * @param t_S second rank tensor S
   * @return auto second direvatives, forth order tensor with minor simetries
   */
  template <typename T>
  static inline auto getDiffDiffMat(Val &t_val, Vec &t_vec, Fun f, Fun d_f,
                                    Fun dd_f, T &t_S) {
    using V =
        typename FTensor::promote<decltype(t_val(0)), decltype(t_vec(0, 0))>::V;
    FTensor::Ddg<V, Dim, Dim> t_diff_A;
    getDiffDiffMatImpl<EigenProjection<T1, T2, NB, Dim>, V>::set(
        t_val, t_vec, f, d_f, dd_f, t_S, t_diff_A, Number<Dim>(), Number<Dim>(),
        Number<Dim>(), Number<Dim>());
    return t_diff_A;
  }

  // private:
  template <typename E, typename C> friend struct d2MCoefficients;
  template <typename E, typename C, typename G, int a, int c, int d, int i,
            int j, int k, int l>
  friend struct d2MImpl;
  template <typename E, typename C, typename G1, typename G2, int a, int i,
            int j, int k, int l, int m, int n>
  friend struct dd4MImpl;
  template <typename E, typename C, int i, int j>
  friend struct reconstructMatImpl;
  template <typename E, typename C, int i, int j, int k, int l>
  friend struct firstMatrixDirectiveImpl;
  template <typename E, typename C, int i, int j, int k, int l, int m, int n>
  friend struct secondMatrixDirectiveImpl;
  template <typename E, typename C> friend struct getDiffMatImpl;
  template <typename E, typename C> friend struct getDiffDiffMatImpl;

  template <int a, int i> static inline auto N(Vec &t_vec) {
    return t_vec(a, i);
  }

  template <int a> static inline auto L(Val &t_val, const Number<a> &) {
    return L<a>(t_val);
  }
  template <int a> static inline auto L(Val &t_val) { return t_val(a); }

  template <int a, int b>
  static inline auto F(Val &t_val, const Number<a> &, const Number<b> &) {
    return F<a, b>(t_val);
  }

  template <int a, int b> static inline auto F(Val &t_val) {
    return static_cast<decltype(t_val(0))>(1) / (L<a>(t_val) - L<b>(t_val));
  }

  template <int a, int i, int j>
  static inline auto M(Vec &t_vec, const Number<a> &, const Number<i> &,
                       const Number<j> &) {
    return M<a, i, j>(t_vec);
  }

  template <int a, int i, int j> static inline auto M(Vec &t_vec) {
    return N<a, i>(t_vec) * N<a, j>(t_vec);
  }

  template <int a, int b, int i, int j, int k, int l>
  static inline auto G(Vec &t_vec, const Number<a> &, const Number<b> &,
                       const Number<i> &, const Number<j> &, const Number<k> &,
                       const Number<l> &) {
    return M<a, i, k>(t_vec) * M<b, j, l>(t_vec) +
           M<a, i, l>(t_vec) * M<b, j, k>(t_vec);
  }

  template <int a, int b, int i, int j, int k, int l>
  static inline auto G(Vec &t_vec) {
    return M<a, i, k>(t_vec) * M<b, j, l>(t_vec) +
           M<a, i, l>(t_vec) * M<b, j, k>(t_vec);
  }

  template <int a, int b, int i, int j>
  static inline auto dFdN(Val &t_val, Vec &t_vec, const Number<a> &,
                          const Number<b> &, const Number<i> &,
                          const Number<j> &) {
    return dFdN<a, b, i, j>(t_val, t_vec);
  }

  template <int a, int b, int i, int j>
  static inline auto dFdN(Val &t_val, Vec &t_vec) {
    return dFdNa<a, b, i, j>(t_val, t_vec) + dFdNb<a, b, i, j>(t_val, t_vec);
  }

  template <int a, int b, int i, int j>
  static inline auto dFdNa(Val &t_val, Vec &t_vec) {
    return -M<a, i, j>(t_vec) /
           ((L<a>(t_val) - L<b>(t_val)) * (L<a>(t_val) - L<b>(t_val)));
  }

  template <int a, int b, int i, int j>
  static inline auto dFdNb(Val &t_val, Vec &t_vec) {
    return M<b, i, j>(t_vec) /
           ((L<a>(t_val) - L<b>(t_val)) * (L<a>(t_val) - L<b>(t_val)));
  }

  template <int a, int b, int i, int j, int k, int l>
  static inline auto S(Vec &t_vec, const Number<a> &, const Number<b> &,
                       const Number<i> &, const Number<j> &, const Number<k> &,
                       const Number<l> &) {
    return S<a, b, i, j, k, l>(t_vec);
  }

  template <int a, int b, int i, int j, int k, int l>
  static inline auto S(Vec &t_vec) {
    return G<a, b, i, j, k, l>(t_vec) + G<b, a, i, j, k, l>(t_vec);
  }
};
