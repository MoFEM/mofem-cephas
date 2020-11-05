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
                         const Number<3> &, Fun f, Fun d_f, Fun dd_f) {
    return f(E::L(t_val, Number<a>())) * E::F(t_val, Number<a>(), Number<b>());
  }

  template <int a, int b>
  static inline auto get(Val &t_val, const Number<a> &, const Number<b> &,
                         const Number<-1> &, const Number<-1> &,
                         const Number<2> &, Fun f, Fun d_f, Fun dd_f) {
    if (a == 1 || b == 1)
      return get(t_val, Number<a>(), Number<b>(), Number<-1>(), Number<-1>(),
                 Number<3>(), f, d_f, dd_f);
    else
      return get(t_val, Number<a>(), Number<b>(), Number<-1>(), Number<-1>(),
                 Number<1>(), f, d_f, dd_f);
  }

  template <int a, int b>
  static inline auto get(Val &t_val, const Number<a> &, const Number<b> &,
                         const Number<-1> &, const Number<-1> &,
                         const Number<1>, Fun f, Fun d_f, Fun dd_f) {
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
                         const Number<3> &, Fun f, Fun d_f, Fun dd_f) {
    return f(E::L(t_val, Number<c>())) * E::F(t_val, Number<c>(), Number<d>()) *
           E::F(t_val, Number<a>(), Number<b>());
  }

  template <int a, int b, int c, int d>
  static inline auto get(Val &t_val, const Number<a> &, const Number<b> &,
                         const Number<c> &, const Number<d> &,
                         const Number<2> &, Fun f, Fun d_f, Fun dd_f) {

    if ((c == 1 || d == 1) && (a == 1 || b == 1))
      return get(t_val, Number<a>(), Number<b>(), Number<c>(), Number<d>(),
                 Number<3>(), f, d_f, dd_f);

    if (c != 1 && d != 1 && a != 1 && b != 1)
      return get(t_val, Number<a>(), Number<b>(), Number<c>(), Number<d>(),
                 Number<1>(), f, d_f, dd_f);

    if ((c != 1 && d != 1) && (a == 1 || b == 1))
      return d_f(E::L(t_val, Number<c>())) *
             E::F(t_val, Number<a>(), Number<b>()) / static_cast<C>(2);

    if ((c == 1 || d == 1) && (a != 1 && b != 1)) {

      if ((c == 2 && d == 1) || (c == 2 && d == 1))
        return (

                   d_f(E::L(t_val, Number<c>()))

                   -
                   (f(E::L(t_val, Number<c>())) - f(E::L(t_val, Number<d>()))) *
                       E::F(t_val, Number<c>(), Number<d>())

                       ) *
               E::F(t_val, Number<c>(), Number<d>());
    }

    return static_cast<C>(0);
  }

  template <int a, int b, int c, int d>
  static inline auto get(Val &t_val, const Number<a> &, const Number<b> &,
                         const Number<c> &, const Number<d> &, const Number<1>,
                         Fun f, Fun d_f, Fun dd_f) {
    if ((a != b && b != d) && (a != d && b != c))
      return dd_f(E::L(t_val, Number<c>())) / static_cast<C>(4);
    else
      return static_cast<C>(0);
  }
};

template <typename E, typename C> struct dd4MCoefficientsType2 {
  using Val = typename E::Val;
  using Vec = typename E::Vec;
  using Fun = typename E::Fun;

  template <int N> using Number = FTensor::Number<N>;

  dd4MCoefficientsType2(E &e): e(e) {}
  E &e;

  template <int a, int b, int m, int n>
  inline auto get(Val &t_val, Vec &t_vec, const Number<a> &, const Number<b> &,
                  const Number<m> &, const Number<n> &, const Number<3> &,
                  Fun f, Fun d_f, Fun dd_f) {
    return e.fVal(a) * E::dFdN(t_val, t_vec, Number<a>(), Number<b>(),
                               Number<m>(), Number<n>());
  }

  template <int a, int b, int m, int n>
  inline auto get(Val &t_val, Vec &t_vec, const Number<a> &, const Number<b> &,
                  const Number<m> &, const Number<n> &, const Number<2> &,
                  Fun f, Fun d_f, Fun dd_f) {
    if (a == 1 || b == 1)
      return get(t_val, t_vec, Number<a>(), Number<b>(), Number<m>(),
                 Number<n>(), Number<3>(), f, d_f, dd_f);
    else
      return get(t_val, t_vec, Number<a>(), Number<b>(), Number<m>(),
                 Number<n>(), Number<1>(), f, d_f, dd_f);
  }

  template <int a, int b, int m, int n>
  inline auto get(Val &t_val, Vec &t_vec, const Number<a> &, const Number<b> &,
                  const Number<m> &, const Number<n> &, const Number<1>, Fun f,
                  Fun d_f, Fun dd_f) {
    return static_cast<C>(0);
  }
};

template <typename E, typename C, typename G>
struct d2MImpl {
  using Val = typename E::Val;
  using Vec = typename E::Vec;
  using Fun = typename E::Fun;

  template <int N> using Number = FTensor::Number<N>;
  d2MImpl(E &e): e(e) {}
  E &e;

  template <int b, int a, int c, int d, int i, int j, int k, int l>
  static inline C term(Val &t_val, Vec &t_vec, Fun f, Fun d_f, Fun dd_f) {
    if (a != b) {
      return G::get(t_val, Number<a>(), Number<b>(), Number<c>(), Number<d>(),
                    typename E::NumberNb(), f, d_f, dd_f) *
             E::S(t_vec, Number<a>(), Number<b>(), Number<i>(), Number<j>(),
                  Number<k>(), Number<l>());
    }
    return 0;
  }

  template <int nb, int a, int c, int d, int i, int j, int k, int l>
  static inline C eval(Val &t_val, Vec &t_vec, Fun f, Fun d_f, Fun dd_f,
                       const Number<nb> &, const Number<a> &, const Number<c> &,
                       const Number<d> &, const Number<i> &, const Number<j> &,
                       const Number<k> &, const Number<l> &) {
    return term<nb - 1, a, c, d, i, j, k, l>(t_val, t_vec, f, d_f, dd_f) +
           eval(t_val, t_vec, f, d_f, dd_f, Number<nb - 1>(), Number<a>(),
                Number<c>(), Number<d>(), Number<i>(), Number<j>(), Number<k>(),
                Number<l>());
  }

  template <int a, int c, int d, int i, int j, int k, int l>
  static inline C eval(Val &t_val, Vec &t_vec, Fun f, Fun d_f, Fun dd_f,
                       const Number<1> &, const Number<a> &, const Number<c> &,
                       const Number<d> &, const Number<i> &, const Number<j> &,
                       const Number<k> &, const Number<l> &) {
    return term<0, a, c, d, i, j, k, l>(t_val, t_vec, f, d_f, dd_f);
  }
};

template <typename E, typename C, typename G1, typename G2>
struct fdd4MImpl {
  using Val = typename E::Val;
  using Vec = typename E::Vec;
  using Fun = typename E::Fun;

  fdd4MImpl(E &e) : r(e), g2(e), e(e) {}
  d2MImpl<E, C, G1> r;
  G2 g2;
  E &e;

  template <int N> using Number = FTensor::Number<N>;

  template <int a, int b, int A, int I, int J, int K, int L>
  inline auto fd2M(Val &t_val, Vec &t_vec, Fun f, Fun d_f, Fun dd_f) {
    return r.eval(t_val, t_vec, f, d_f, dd_f, Number<3>(), Number<A>(),
                  Number<a>(), Number<b>(), Number<I>(), Number<J>(),
                  Number<K>(), Number<L>());
  }

  template <int a, int b, int A, int B, int I, int J, int K, int L, int M,
            int N>
  inline auto fd2G(Val &t_val, Vec &t_vec, Fun f, Fun d_f, Fun dd_f) {
    return fd2M<a, b, A, I, K, N, M>(t_val, t_vec, f, d_f, dd_f) *
               e.aM(B, J, L) +
           e.aM(A, I, K) *
               fd2M<a, b, B, J, L, M, N>(t_val, t_vec, f, d_f, dd_f) +
           fd2M<a, b, A, I, L, M, N>(t_val, t_vec, f, d_f, dd_f) *
               e.aM(B, J, K) +
           e.aM(A, I, L) *
               fd2M<a, b, B, J, K, M, N>(t_val, t_vec, f, d_f, dd_f);
  }

  template <int a, int A, int B, int I, int J, int K, int L, int M, int N>
  inline auto fd2S(Val &t_val, Vec &t_vec, Fun f, Fun d_f, Fun dd_f) {
    return fd2G<a, B, A, B, I, J, K, L, M, N>(t_val, t_vec, f, d_f, dd_f) +
           fd2G<a, B, B, A, I, J, K, L, M, N>(t_val, t_vec, f, d_f, dd_f);
  }

  template <int a, int b, int i, int j, int k, int l, int m, int n>
  inline C term(Val &t_val, Vec &t_vec, Fun f, Fun d_f, Fun dd_f) {

    if (a != b) {

      return

          fd2S<a, a, b, i, j, k, l, m, n>(t_val, t_vec, f, d_f, dd_f)

          +

          2 *

              g2.get(t_val, t_vec, Number<a>(), Number<b>(), Number<m>(),
                     Number<n>(), typename E::NumberNb(), f, d_f, dd_f) *
              E::S(t_vec, Number<a>(), Number<b>(), Number<i>(), Number<j>(),
                   Number<k>(), Number<l>());
    }

    return 0;
  }

  template <int nb, int a, int i, int j, int k, int l, int m, int n>
  inline C eval(Val &t_val, Vec &t_vec, Fun f, Fun d_f, Fun dd_f,
                const Number<nb> &, const Number<a> &, const Number<i> &,
                const Number<j> &, const Number<k> &, const Number<l> &,
                const Number<m> &, const Number<n> &) {
    return term<a, nb - 1, i, j, k, l, m, n>(t_val, t_vec, f, d_f, dd_f) +
           eval(t_val, t_vec, f, d_f, dd_f, Number<nb - 1>(), Number<a>(),
                Number<i>(), Number<j>(), Number<k>(), Number<l>(), Number<m>(),
                Number<n>());
  }

  template <int a, int i, int j, int k, int l, int m, int n>
  inline C eval(Val &t_val, Vec &t_vec, Fun f, Fun d_f, Fun dd_f,
                const Number<1> &, const Number<a> &, const Number<i> &,
                const Number<j> &, const Number<k> &, const Number<l> &,
                const Number<m> &, const Number<n> &) {
    return term<a, 0, i, j, k, l, m, n>(t_val, t_vec, f, d_f, dd_f);
  }
};

template <typename E, typename C> struct reconstructMatImpl {
  using Val = typename E::Val;
  using Vec = typename E::Vec;
  using Fun = typename E::Fun;

  template <int N> using Number = FTensor::Number<N>;

  reconstructMatImpl(E &e) : e(e) {}
  E &e;

  template <int a, int i, int j> inline C term() {
    return e.aM(a, i, j) * e.fVal(a);
  }

  template <int nb, int i, int j>
  inline C eval(const Number<nb> &, const Number<i> &, const Number<j> &) {
    return term<nb - 1, i, j>() +
           eval(Number<nb - 1>(), Number<i>(), Number<j>());
  }

  template <int i, int j>
  inline C eval(const Number<1> &, const Number<i> &, const Number<j> &) {
    return term<0, i, j>();
  }
};

template <typename E, typename C> struct firstMatrixDirectiveImpl {
  using Val = typename E::Val;
  using Vec = typename E::Vec;
  using Fun = typename E::Fun;

  template <int N> using Number = FTensor::Number<N>;

  firstMatrixDirectiveImpl(E &e) : r(e), e(e) {}
  d2MImpl<E, C, d2MCoefficients<E, C>> r;
  E &e;

  template <int a, int i, int j, int k, int l>
  inline C term(Val &t_val, Vec &t_vec, Fun f, Fun d_f) {
    return

        e.aM(a, i, j) * e.aM(a, k, l) * e.dfVal(a)

        +

        r.eval(t_val, t_vec, f, d_f, d_f, Number<3>(), Number<a>(),
               Number<-1>(), Number<-1>(), Number<i>(), Number<j>(),
               Number<k>(), Number<l>()) /
            static_cast<C>(2);
  }

  template <int nb, int i, int j, int k, int l>
  inline C eval(Val &t_val, Vec &t_vec, Fun f, Fun d_f, const Number<nb> &,
                const Number<i> &, const Number<j> &, const Number<k> &,
                const Number<l> &) {
    return term<nb - 1, i, j, k, l>(t_val, t_vec, f, d_f) +
           eval(t_val, t_vec, f, d_f, Number<nb - 1>(), Number<i>(),
                Number<j>(), Number<k>(), Number<l>());
  }

  template <int i, int j, int k, int l>
  inline C eval(Val &t_val, Vec &t_vec, Fun f, Fun d_f, const Number<1> &,
                const Number<i> &, const Number<j> &, const Number<k> &,
                const Number<l> &) {
    return term<0, i, j, k, l>(t_val, t_vec, f, d_f);
  }
};

template <typename E, typename C>
struct secondMatrixDirectiveImpl {
  using Val = typename E::Val;
  using Vec = typename E::Vec;
  using Fun = typename E::Fun;

  template <int N> using Number = FTensor::Number<N>;

  secondMatrixDirectiveImpl(E &e) : w(e), r(e), e(e) {}
  d2MImpl<E, C, d2MCoefficients<E, C>> w;
  fdd4MImpl<E, C, dd4MCoefficientsType1<E, C>, dd4MCoefficientsType2<E, C>> r;
  E &e;

  template <int a, int i, int j, int k, int l, int m, int n>
  inline C term(Fun f, Fun d_f, Fun dd_f) {

    return

        (

            w.eval(e.tVal, e.tVec, d_f, dd_f, dd_f, Number<3>(), Number<a>(),
                   Number<-1>(), Number<-1>(), Number<i>(), Number<j>(),
                   Number<m>(), Number<n>()) *
                e.aM(a, k, l)

            +

            e.aM(a, i, j) * w.eval(e.tVal, e.tVec, d_f, dd_f, dd_f, Number<3>(),
                                   Number<a>(), Number<-1>(), Number<-1>(),
                                   Number<k>(), Number<l>(), Number<m>(),
                                   Number<n>())) /
            static_cast<C>(2) +

        e.aM(a, i, j) * e.aM(a, k, l) * e.aM(a, m, n) * e.ddfVal(a)

        +

        r.eval(e.tVal, e.tVec, f, d_f, dd_f, Number<3>(), Number<a>(),
               Number<i>(), Number<j>(), Number<k>(), Number<l>(), Number<m>(),
               Number<n>()) /
            static_cast<C>(4) +

        w.eval(e.tVal, e.tVec, d_f, dd_f, dd_f, Number<3>(), Number<a>(),
               Number<-1>(), Number<-1>(), Number<i>(), Number<j>(),
               Number<k>(), Number<l>()) *
            e.aM(a, m, n) / static_cast<C>(2);
  }

  template <int nb, int i, int j, int k, int l, int m, int n>
  inline C eval(Fun f, Fun d_f, Fun dd_f, const Number<nb> &, const Number<i> &,
                const Number<j> &, const Number<k> &, const Number<l> &,
                const Number<m> &, const Number<n> &) {
    return term<nb - 1, i, j, k, l, m, n>(f, d_f, dd_f)

           +

           eval(f, d_f, dd_f, Number<nb - 1>(), Number<i>(), Number<j>(),
                Number<k>(), Number<l>(), Number<m>(), Number<n>());
  }

  template <int i, int j, int k, int l, int m, int n>
  inline C eval(Fun f, Fun d_f, Fun dd_f, const Number<1> &, const Number<i> &,
                const Number<j> &, const Number<k> &, const Number<l> &,
                const Number<m> &, const Number<n>) {
    return term<0, i, j, k, l, m, n>(f, d_f, dd_f);
  }
};

template <typename E, typename C, typename T> struct getMatImpl {
  using Val = typename E::Val;
  using Vec = typename E::Vec;
  using Fun = typename E::Fun;

  template <int N> using Number = FTensor::Number<N>;

  getMatImpl(E &e, T &t_a) : r(e), tA(t_a) {}
  reconstructMatImpl<E, C> r;
  T &tA;

  template <int I, int J>
  inline void set(const Number<I> &, const Number<J> &) {
    set(Number<I>(), Number<J - 1>());
    tA(I - 1, J - 1) = r.eval(Number<3>(), Number<I - 1>(), Number<J - 1>());
  }

  template <int I> inline void set(const Number<I> &, const Number<0> &) {
    set(Number<I - 1>(), Number<I - 1>());
  }

  inline void set(const Number<0> &, const Number<0> &) {}
};

template <typename E, typename C, typename T> struct getDiffMatImpl {
  using Val = typename E::Val;
  using Vec = typename E::Vec;
  using Fun = typename E::Fun;

  template <int N> using Number = FTensor::Number<N>;

  getDiffMatImpl(E &e, T &t_a) : r(e), tA(t_a) {}
  firstMatrixDirectiveImpl<E, C> r;
  T &tA;

  template <int I, int J, int K, int L>
  inline void set(Val &t_val, Vec &t_vec, Fun f, Fun d_f, T &t_a,
                  const Number<I> &, const Number<J> &, const Number<K> &,
                  const Number<L> &) {
    set(t_val, t_vec, f, d_f, t_a, Number<I>(), Number<J>(), Number<K>(),
        Number<L - 1>());
    t_a(I - 1, J - 1, K - 1, L - 1) =
        r.eval(t_val, t_vec, f, d_f, Number<3>(), Number<I - 1>(),
               Number<J - 1>(), Number<K - 1>(), Number<L - 1>());
  }

  template <int I, int J, int K>
  inline void set(Val &t_val, Vec &t_vec, Fun f, Fun d_f, T &t_a,
                  const Number<I> &, const Number<J> &, const Number<K> &,
                  const Number<0> &) {
    set(t_val, t_vec, f, d_f, t_a, Number<I>(), Number<J>(), Number<K - 1>(),
        Number<K - 1>());
  }

  template <int I, int J>
  inline void set(Val &t_val, Vec &t_vec, Fun f, Fun d_f, T &t_a,
                  const Number<I> &, const Number<J> &, const Number<0> &,
                  const Number<0> &) {
    set(t_val, t_vec, f, d_f, t_a, Number<I>(), Number<J - 1>(), Number<3>(),
        Number<3>());
  }

  template <int I, int K, int L>
  inline void set(Val &t_val, Vec &t_vec, Fun f, Fun d_f, T &t_a,
                  const Number<I> &, const Number<0> &, const Number<K> &,
                  const Number<L> &) {
    set(t_val, t_vec, f, d_f, t_a, Number<I - 1>(), Number<I - 1>(),
        Number<K>(), Number<L>());
  }

  template <int K, int L>
  inline void set(Val &t_val, Vec &t_vec, Fun f, Fun d_f, T &t_a,
                  const Number<0> &, const Number<0> &, const Number<K> &,
                  const Number<L> &) {}
};

template <typename E, typename C, typename T1, typename T2>
struct getDiffDiffMatImpl {
  using Val = typename E::Val;
  using Vec = typename E::Vec;
  using Fun = typename E::Fun;

  template <int N> using Number = FTensor::Number<N>;

  getDiffDiffMatImpl(E &e, T1 &t_a, T2 &t_S) : r(e), e(e), tA(t_a), tS(t_S) {}
  secondMatrixDirectiveImpl<E, C> r;
  E &e;
  T1 &tA;
  T2 &tS;

  template <int I, int J, int K, int L, int M, int N>
  inline auto add(Fun f, Fun d_f, Fun dd_f, const Number<I> &,
                  const Number<J> &, const Number<K> &, const Number<L> &,
                  const Number<M> &, const Number<N> &) {
    return tS(M - 1, N - 1) * r.eval(f, d_f, dd_f, Number<3>(), Number<M - 1>(),
                                     Number<N - 1>(), Number<I - 1>(),
                                     Number<J - 1>(), Number<K - 1>(),
                                     Number<L - 1>())

           +

           add(f, d_f, dd_f, Number<I>(), Number<J>(), Number<K>(), Number<L>(),
               Number<M>(), Number<N - 1>());
  }

  template <int I, int J, int K, int L, int M>
  inline auto add(Fun f, Fun d_f, Fun dd_f, const Number<I> &,
                  const Number<J> &, const Number<K> &, const Number<L> &,
                  const Number<M> &, const Number<1> &) {
    return tS(M - 1, 0) * r.eval(f, d_f, dd_f, Number<3>(), Number<M - 1>(),
                                 Number<0>(), Number<I - 1>(), Number<J - 1>(),
                                 Number<K - 1>(), Number<L - 1>())

           +

           add(f, d_f, dd_f, Number<I>(), Number<J>(), Number<K>(), Number<L>(),
               Number<M - 1>(), Number<3>());
  }

  template <int I, int J, int K, int L>
  inline auto add(Fun f, Fun d_f, Fun dd_f, const Number<I> &, const Number<J> &,
                  const Number<K> &, const Number<L> &, const Number<1> &,
                  const Number<1> &) {
    return tS(0, 0) * r.eval(f, d_f, dd_f, Number<3>(), Number<0>(),
                             Number<0>(), Number<I - 1>(), Number<J - 1>(),
                             Number<K - 1>(), Number<L - 1>());
  }

  template <int I, int J, int K, int L>
  inline void set(Fun f, Fun d_f, Fun dd_f, const Number<I> &,
                  const Number<J> &, const Number<K> &, const Number<L> &) {
    set(f, d_f, dd_f, Number<I>(), Number<J>(), Number<K>(), Number<L - 1>());
    tA(I - 1, J - 1, K - 1, L - 1) =
        add(f, d_f, dd_f, Number<I>(), Number<J>(), Number<K>(), Number<L>(),
            Number<3>(), Number<3>());
  }

  template <int I, int J, int K>
  inline void set(Fun f, Fun d_f, Fun dd_f, const Number<I> &,
                  const Number<J> &, const Number<K> &, const Number<0> &) {
    set(f, d_f, dd_f, Number<I>(), Number<J>(), Number<K - 1>(),
        Number<K - 1>());
  }

  template <int I, int J>
  inline void set(Fun f, Fun d_f, Fun dd_f, const Number<I> &,
                  const Number<J> &, const Number<0> &, const Number<0> &) {
    set(f, d_f, dd_f, Number<I>(), Number<J - 1>(), Number<3>(), Number<3>());
  }

  template <int I, int K, int L>
  inline void set(Fun f, Fun d_f, Fun dd_f, const Number<I> &,
                  const Number<0> &, const Number<K> &, const Number<L> &) {
    set(f, d_f, dd_f, Number<I - 1>(), Number<I - 1>(), Number<K>(),
        Number<L>());
  }

  template <int K, int L>
  inline void set(Fun f, Fun d_f, Fun dd_f, const Number<0> &,
                  const Number<0> &, const Number<K> &, const Number<L> &) {}
};

template <typename T1, typename T2, int NB> struct EigenProjection {

  static constexpr int Dim = 3;

  using Val = const FTensor::Tensor1<T1, Dim>;
  using Vec = const FTensor::Tensor2<T2, Dim, Dim>;
  using Fun = boost::function<double(const double)>;

  template <int N> using Number = FTensor::Number<N>;
  template <char c> using I = typename FTensor::Index<c, Dim>;

  using NumberNb = Number<NB>;
  using NumberDim = Number<Dim>;

  EigenProjection(Val &t_val, Vec &t_vec) : tVal(t_val), tVec(t_vec) {

    for (auto aa : {0, 1, 2})
      for (auto ii : {0, 1, 2})
        for (auto jj = 0; jj <= ii; ++jj)
          aM(aa, ii, jj) = tVec(aa, ii) * tVec(aa, jj);
  }

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
  inline auto getMat(Fun f) {

    for (auto aa : {0, 1, 2})
      fVal(aa) = f(tVal(aa));

    using V =
        typename FTensor::promote<decltype(tVal(0)), decltype(tVec(0, 0))>::V;
    using T3 =
        FTensor::Tensor2_symmetric<typename std::remove_const<V>::type, Dim>;
    T3 t_A;
    getMatImpl<EigenProjection<T1, T2, NB>, V, T3>(*this, t_A)
        .set(Number<Dim>(), Number<Dim>());
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
  inline auto getDiffMat(Val &t_val, Vec &t_vec, Fun f, Fun d_f) {

    for (auto aa : {0, 1, 2})
      fVal(aa) = f(tVal(aa));

    for (auto aa : {0, 1, 2})
      dfVal(aa) = d_f(tVal(aa));

    using V =
        typename FTensor::promote<decltype(t_val(0)), decltype(t_vec(0, 0))>::V;
    using T3 = FTensor::Ddg<V, Dim, Dim>;
    T3 t_diff_A;

    getDiffMatImpl<EigenProjection<T1, T2, NB>, V, T3>(*this, t_diff_A)
        .set(t_val, t_vec, f, d_f, t_diff_A, Number<Dim>(), Number<Dim>(),
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
  inline auto getDiffDiffMat(Val &t_val, Vec &t_vec, Fun f, Fun d_f, Fun dd_f,
                             T &t_S) {

    for (auto aa : {0, 1, 2})
      fVal(aa) = f(tVal(aa));

    for (auto aa : {0, 1, 2})
      dfVal(aa) = d_f(tVal(aa));

    for (auto aa : {0, 1, 2})
      ddfVal(aa) = dd_f(tVal(aa));

    using V =
        typename FTensor::promote<decltype(t_val(0)), decltype(t_vec(0, 0))>::V;
    using T3 = FTensor::Ddg<V, Dim, Dim>;
    T3 t_diff_A;
    getDiffDiffMatImpl<EigenProjection<T1, T2, NB>, V, T3, T>(*this, t_diff_A,
                                                              t_S)
        .set(f, d_f, dd_f, Number<Dim>(), Number<Dim>(), Number<Dim>(),
             Number<Dim>());
    return t_diff_A;
  }

private:
  Val &tVal;
  Vec &tVec;
  FTensor::Christof<T2, Dim, Dim> aM;
  FTensor::Tensor1<T1, Dim> fVal;
  FTensor::Tensor1<T1, Dim> dfVal;
  FTensor::Tensor1<T1, Dim> ddfVal;

  template <typename E, typename C> friend struct d2MCoefficients;
  template <typename E, typename C> friend struct dd4MCoefficientsType1;
  template <typename E, typename C> friend struct dd4MCoefficientsType2;
  template <typename E, typename C, typename G> friend struct d2MImpl;
  template <typename E, typename C, typename G1, typename G2>
  friend struct fdd4MImpl;
  template <typename E, typename C> friend struct reconstructMatImpl;
  template <typename E, typename C> friend struct firstMatrixDirectiveImpl;
  template <typename E, typename C> friend struct secondMatrixDirectiveImpl;
  template <typename E, typename C, typename T> friend struct getDiffMatImpl;
  template <typename E, typename C, typename T3, typename T4>
  friend struct getDiffDiffMatImpl;

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
