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
 auto  t_A = EigenMatrixImp<double, double, 3, 3>::getMat(t_L, t_N, f);
 \endcode
 where <3> means that are three unique eigen values. Return t_A is symmetric
 tensor rank two.

 Calculate directive
 \code
 auto t_P = EigenMatrixImp<double, double, 3, 3>::getDiffMat(t_L, t_N, f, d_f);
 \endcode
 where return t_SL is 4th order tensor (symmetry on first two and
 second to indices, i.e. minor symmetrise)

 Calculate second derivative, L, such that S:L, for given S,
 \code
 FTensor::Tensor2<double, 3, 3> t_S{

    1., 0., 0.,

    0., 1., 0.,

    0., 0., 1.};

  auto t_SL = EigenMatrixImp<double, double, 3, 3>::getDiffDiffMat(
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

namespace EigenMatrix {

template <typename E, typename C> struct d2MCoefficients {
  using Val = typename E::Val;
  using Vec = typename E::Vec;
  using Fun = typename E::Fun;

  template <int N> using Number = FTensor::Number<N>;

  d2MCoefficients(E &e) : e(e) {}
  E &e;

  template <int a, int b>
  inline auto get(const Number<a> &, const Number<b> &, const Number<-1> &,
                  const Number<-1> &, const Number<3> &) const {
    return e.fVal(a) * e.aF(a, b);
  }

  template <int a, int b>
  inline auto get(const Number<a> &, const Number<b> &, const Number<-1> &,
                  const Number<-1> &, const Number<2> &) const {
    if (a == 1 || b == 1)
      return get(Number<a>(), Number<b>(), Number<-1>(), Number<-1>(),
                 Number<3>());
    else
      return get(Number<a>(), Number<b>(), Number<-1>(), Number<-1>(),
                 Number<1>());
  }

  template <int a, int b>
  inline auto get(const Number<a> &, const Number<b> &, const Number<-1> &,
                  const Number<-1> &, const Number<1>) const {
    return e.dfVal(a) / static_cast<C>(2);
  }
};

template <typename E, typename C> struct d2MCoefficientsType0 {
  using Val = typename E::Val;
  using Vec = typename E::Vec;
  using Fun = typename E::Fun;

  template <int N> using Number = FTensor::Number<N>;

  d2MCoefficientsType0(E &e) : e(e) {}
  E &e;

  template <int a, int b>
  inline auto get(const Number<a> &, const Number<b> &, const Number<-1> &,
                  const Number<-1> &, const Number<3> &) const {
    return e.coefficientsType0(a, b);
  }

  template <int a, int b>
  inline auto get(const Number<a> &, const Number<b> &, const Number<-1> &,
                  const Number<-1> &, const Number<2> &) const {
    if (a == 1 || b == 1)
      return get(Number<a>(), Number<b>(), Number<-1>(), Number<-1>(),
                 Number<3>());
    else
      return get(Number<a>(), Number<b>(), Number<-1>(), Number<-1>(),
                 Number<1>());
  }

  template <int a, int b>
  inline auto get(const Number<a> &, const Number<b> &, const Number<-1> &,
                  const Number<-1> &, const Number<1> &) const {
    return e.ddfVal(a) / static_cast<C>(2);
  }
};

template <typename E, typename C> struct dd4MCoefficientsType1 {
  using Val = typename E::Val;
  using Vec = typename E::Vec;
  using Fun = typename E::Fun;

  template <int N> using Number = FTensor::Number<N>;

  dd4MCoefficientsType1(E &e) : e(e) {}
  E &e;

  template <int a, int b, int c, int d>
  inline auto get(const Number<a> &, const Number<b> &, const Number<c> &,
                  const Number<d> &, const Number<3> &) const {
    return e.coefficientsType1(a, b, c, d);
  }

  template <int a, int b, int c, int d>
  inline auto get(const Number<a> &, const Number<b> &, const Number<c> &,
                  const Number<d> &, const Number<2> &) const {
    return e.coefficientsType1(a, b, c, d);
  }

  template <int a, int b, int c, int d>
  inline auto get(const Number<a> &, const Number<b> &, const Number<c> &,
                  const Number<d> &, const Number<1>) const {
    return e.coefficientsType1(a, b, c, d);
  }
};

template <typename E, typename C> struct dd4MCoefficientsType2 {
  using Val = typename E::Val;
  using Vec = typename E::Vec;
  using Fun = typename E::Fun;

  template <int N> using Number = FTensor::Number<N>;

  dd4MCoefficientsType2(E &e) : e(e) {}
  E &e;

  template <int a, int b, int m, int n>
  inline auto get(const Number<a> &, const Number<b> &, const Number<m> &,
                  const Number<n> &, const Number<3> &) const {
    return e.fVal(a) *
           e.dFdN(Number<a>(), Number<b>(), Number<m>(), Number<n>());
  }

  template <int a, int b, int m, int n>
  inline auto get(const Number<a> &, const Number<b> &, const Number<m> &,
                  const Number<n> &, const Number<2> &) const {
    if (a == 1 || b == 1)
      return get(Number<a>(), Number<b>(), Number<m>(), Number<n>(),
                 Number<3>());
    else
      return get(Number<a>(), Number<b>(), Number<m>(), Number<n>(),
                 Number<1>());
  }

  template <int a, int b, int m, int n>
  inline auto get(const Number<a> &, const Number<b> &, const Number<m> &,
                  const Number<n> &, const Number<1>) const {
    return static_cast<C>(0);
  }
};

template <typename E, typename C, typename G> struct d2MImpl {
  using Val = typename E::Val;
  using Vec = typename E::Vec;
  using Fun = typename E::Fun;

  template <int N> using Number = FTensor::Number<N>;
  d2MImpl(E &e) : g(e), e(e) {}
  G g;
  E &e;

  template <int b, int a, int c, int d, int i, int j, int k, int l>
  inline C term() const {
    if (a != b) {
      return g.get(Number<a>(), Number<b>(), Number<c>(), Number<d>(),
                   typename E::NumberNb()) *
             e.S(Number<a>(), Number<b>(), Number<i>(), Number<j>(),
                 Number<k>(), Number<l>());
    }
    return 0;
  }

  template <int nb, int a, int c, int d, int i, int j, int k, int l>
  inline C eval(const Number<nb> &, const Number<a> &, const Number<c> &,
                const Number<d> &, const Number<i> &, const Number<j> &,
                const Number<k> &, const Number<l> &) const {
    return term<nb - 1, a, c, d, i, j, k, l>() +
           eval(Number<nb - 1>(), Number<a>(), Number<c>(), Number<d>(),
                Number<i>(), Number<j>(), Number<k>(), Number<l>());
  }

  template <int a, int c, int d, int i, int j, int k, int l>
  inline C eval(const Number<1> &, const Number<a> &, const Number<c> &,
                const Number<d> &, const Number<i> &, const Number<j> &,
                const Number<k> &, const Number<l> &) const {
    return term<0, a, c, d, i, j, k, l>();
  }
};

template <typename E, typename C, typename G1, typename G2> struct Fdd4MImpl {
  using Val = typename E::Val;
  using Vec = typename E::Vec;
  using Fun = typename E::Fun;

  Fdd4MImpl(E &e) : r(e), g2(e), e(e) {}
  d2MImpl<E, C, G1> r;
  G2 g2;
  E &e;

  template <int N> using Number = FTensor::Number<N>;

  template <int a, int b, int A, int I, int J, int K, int L>
  inline auto fd2M() const {
    return r.eval(typename E::NumberDim(), Number<A>(), Number<a>(),
                  Number<b>(), Number<I>(), Number<J>(), Number<K>(),
                  Number<L>());
  }

  template <int a, int b, int A, int B, int I, int J, int K, int L, int M,
            int N>
  inline auto fd2G() const {
    return fd2M<a, b, A, I, K, N, M>() * e.aM(B, J, L) +
           e.aM(A, I, K) * fd2M<a, b, B, J, L, M, N>() +
           fd2M<a, b, A, I, L, M, N>() * e.aM(B, J, K) +
           e.aM(A, I, L) * fd2M<a, b, B, J, K, M, N>();
  }

  template <int a, int B, int I, int J, int K, int L, int M, int N>
  inline auto fd2S() const {
    return fd2G<a, B, a, B, I, J, K, L, M, N>() +
           fd2G<a, B, B, a, I, J, K, L, M, N>();
  }

  template <int a, int b, int i, int j, int k, int l, int m, int n>
  inline C term() const {

    if (a != b) {

      return

          fd2S<a, b, i, j, k, l, m, n>()

          +

          2 *

              g2.get(Number<a>(), Number<b>(), Number<m>(), Number<n>(),
                     typename E::NumberNb()) * 
              e.S(Number<a>(), Number<b>(), Number<i>(), Number<j>(),
                  Number<k>(), Number<l>());
    }

    return 0;
  }

  template <int nb, int a, int i, int j, int k, int l, int m, int n>
  inline C eval(const Number<nb> &, const Number<a> &, const Number<i> &,
                const Number<j> &, const Number<k> &, const Number<l> &,
                const Number<m> &, const Number<n> &) const {
    return term<a, nb - 1, i, j, k, l, m, n>() +
           eval(Number<nb - 1>(), Number<a>(), Number<i>(), Number<j>(),
                Number<k>(), Number<l>(), Number<m>(), Number<n>());
  }

  template <int a, int i, int j, int k, int l, int m, int n>
  inline C eval(const Number<1> &, const Number<a> &, const Number<i> &,
                const Number<j> &, const Number<k> &, const Number<l> &,
                const Number<m> &, const Number<n> &) const {
    return term<a, 0, i, j, k, l, m, n>();
  }
};

template <typename E, typename C> struct ReconstructMatImpl {
  using Val = typename E::Val;
  using Vec = typename E::Vec;
  using Fun = typename E::Fun;

  template <int N> using Number = FTensor::Number<N>;

  ReconstructMatImpl(E &e) : e(e) {}
  E &e;

  template <int a, int i, int j> inline C term() const {
    return e.aM(a, i, j) * e.fVal(a);
  }

  template <int nb, int i, int j>
  inline C eval(const Number<nb> &, const Number<i> &,
                const Number<j> &) const {
    return term<nb - 1, i, j>() +
           eval(Number<nb - 1>(), Number<i>(), Number<j>());
  }

  template <int i, int j>
  inline C eval(const Number<1> &, const Number<i> &, const Number<j> &) const {
    return term<0, i, j>();
  }
};

template <typename E, typename C> struct FirstMatrixDirectiveImpl {
  using Val = typename E::Val;
  using Vec = typename E::Vec;
  using Fun = typename E::Fun;

  template <int N> using Number = FTensor::Number<N>;

  FirstMatrixDirectiveImpl(E &e) : r(e), e(e) {}
  d2MImpl<E, C, d2MCoefficients<E, C>> r;
  E &e;

  template <int a, int i, int j, int k, int l> inline C term() const {
    return

        e.aMM[a][a](i, j, k, l) * e.dfVal(a)

        +

        r.eval(typename E::NumberDim(), Number<a>(), Number<-1>(), Number<-1>(),
               Number<i>(), Number<j>(), Number<k>(), Number<l>()) /
            static_cast<C>(2);
  }

  template <int nb, int i, int j, int k, int l>
  inline C eval(const Number<nb> &, const Number<i> &, const Number<j> &,
                const Number<k> &, const Number<l> &) const {
    return term<nb - 1, i, j, k, l>() + eval(Number<nb - 1>(), Number<i>(),
                                             Number<j>(), Number<k>(),
                                             Number<l>());
  }

  template <int i, int j, int k, int l>
  inline C eval(const Number<1> &, const Number<i> &, const Number<j> &,
                const Number<k> &, const Number<l> &) const {
    return term<0, i, j, k, l>();
  }
};

template <typename E, typename C> struct SecondMatrixDirectiveImpl {
  using Val = typename E::Val;
  using Vec = typename E::Vec;
  using Fun = typename E::Fun;

  template <int N> using Number = FTensor::Number<N>;

  SecondMatrixDirectiveImpl(E &e) : w(e), r(e), e(e) {}
  d2MImpl<E, C, d2MCoefficientsType0<E, C>> w;
  Fdd4MImpl<E, C, dd4MCoefficientsType1<E, C>, dd4MCoefficientsType2<E, C>> r;
  E &e;

  template <int a, int i, int j, int k, int l, int m, int n>
  inline C term() const {

    return

        (

            w.eval(typename E::NumberDim(), Number<a>(), Number<-1>(),
                   Number<-1>(), Number<i>(), Number<j>(), Number<m>(),
                   Number<n>()) *
                e.aM(a, k, l)

            +

            e.aM(a, i, j) * w.eval(typename E::NumberDim(), Number<a>(),
                                   Number<-1>(), Number<-1>(), Number<k>(),
                                   Number<l>(), Number<m>(), Number<n>())) /
            static_cast<C>(2) +

        e.aMM[a][a](i, j, k, l) * e.aM(a, m, n) * e.ddfVal(a)

        +

        r.eval(typename E::NumberDim(), Number<a>(), Number<i>(), Number<j>(),
               Number<k>(), Number<l>(), Number<m>(), Number<n>()) /
            static_cast<C>(4) +

        w.eval(typename E::NumberDim(), Number<a>(), Number<-1>(), Number<-1>(),
               Number<i>(), Number<j>(), Number<k>(), Number<l>()) *
            e.aM(a, m, n) / static_cast<C>(2);
  }

  template <int nb, int i, int j, int k, int l, int m, int n>
  inline C eval(const Number<nb> &, const Number<i> &, const Number<j> &,
                const Number<k> &, const Number<l> &, const Number<m> &,
                const Number<n> &) const {
    return term<nb - 1, i, j, k, l, m, n>()

           +

           eval(Number<nb - 1>(), Number<i>(), Number<j>(), Number<k>(),
                Number<l>(), Number<m>(), Number<n>());
  }

  template <int i, int j, int k, int l, int m, int n>
  inline C eval(const Number<1> &, const Number<i> &, const Number<j> &,
                const Number<k> &, const Number<l> &, const Number<m> &,
                const Number<n>) const {
    return term<0, i, j, k, l, m, n>();
  }
};

template <typename E, typename C, typename T> struct GetMatImpl {
  using Val = typename E::Val;
  using Vec = typename E::Vec;
  using Fun = typename E::Fun;

  template <int N> using Number = FTensor::Number<N>;

  GetMatImpl(E &e, T &t_a) : r(e), tA(t_a) {}
  ReconstructMatImpl<E, C> r;
  T &tA;

  template <int I, int J>
  inline void set(const Number<I> &, const Number<J> &) {
    set(Number<I>(), Number<J - 1>());
    tA(I - 1, J - 1) =
        r.eval(typename E::NumberDim(), Number<I - 1>(), Number<J - 1>());
  }

  template <int I> inline void set(const Number<I> &, const Number<0> &) {
    set(Number<I - 1>(), Number<I - 1>());
  }

  inline void set(const Number<0> &, const Number<0> &) {}
};

template <typename E, typename C, typename T> struct GetDiffMatImpl {
  using Val = typename E::Val;
  using Vec = typename E::Vec;
  using Fun = typename E::Fun;

  template <int N> using Number = FTensor::Number<N>;

  GetDiffMatImpl(E &e, T &t_a) : r(e), tA(t_a) {}
  FirstMatrixDirectiveImpl<E, C> r;
  T &tA;

  template <int I, int J, int K, int L>
  inline void set(const Number<I> &, const Number<J> &, const Number<K> &,
                  const Number<L> &) {
    set(Number<I>(), Number<J>(), Number<K>(), Number<L - 1>());
    tA(I - 1, J - 1, K - 1, L - 1) =
        r.eval(typename E::NumberDim(), Number<I - 1>(), Number<J - 1>(),
               Number<K - 1>(), Number<L - 1>());
  }

  template <int I, int J, int K>
  inline void set(const Number<I> &, const Number<J> &, const Number<K> &,
                  const Number<0> &) {
    set(Number<I>(), Number<J>(), Number<K - 1>(), Number<K - 1>());
  }

  template <int I, int J>
  inline void set(const Number<I> &, const Number<J> &, const Number<0> &,
                  const Number<0> &) {
    set(Number<I>(), Number<J - 1>(), typename E::NumberDim(),
        typename E::NumberDim());
  }

  template <int I, int K, int L>
  inline void set(const Number<I> &, const Number<0> &, const Number<K> &,
                  const Number<L> &) {
    set(Number<I - 1>(), Number<I - 1>(), Number<K>(), Number<L>());
  }

  template <int K, int L>
  inline void set(const Number<0> &, const Number<0> &, const Number<K> &,
                  const Number<L> &) {}
};

template <typename E, typename C, typename T1, typename T2>
struct GetDiffDiffMatImpl {
  using Val = typename E::Val;
  using Vec = typename E::Vec;
  using Fun = typename E::Fun;

  template <int N> using Number = FTensor::Number<N>;

  GetDiffDiffMatImpl(E &e, T1 &t_a, T2 &t_S) : r(e), e(e), tA(t_a), tS(t_S) {}
  SecondMatrixDirectiveImpl<E, C> r;
  E &e;
  T1 &tA;
  T2 &tS;

  template <int I, int J, int K, int L, int M, int N>
  inline auto add(const Number<I> &, const Number<J> &, const Number<K> &,
                  const Number<L> &, const Number<M> &, const Number<N> &) {
    if (N != M)
      return (tS(M - 1, N - 1) + tS(N - 1, M - 1)) *
                 r.eval(typename E::NumberDim(), Number<M - 1>(),
                        Number<N - 1>(), Number<I - 1>(), Number<J - 1>(),
                        Number<K - 1>(), Number<L - 1>())

             +

             add(Number<I>(), Number<J>(), Number<K>(), Number<L>(),
                 Number<M>(), Number<N - 1>());
    else
      return tS(M - 1, N - 1) * r.eval(typename E::NumberDim(), Number<M - 1>(),
                                       Number<N - 1>(), Number<I - 1>(),
                                       Number<J - 1>(), Number<K - 1>(),
                                       Number<L - 1>())

             +

             add(Number<I>(), Number<J>(), Number<K>(), Number<L>(),
                 Number<M>(), Number<N - 1>());
  }

  template <int I, int J, int K, int L, int M>
  inline auto add(const Number<I> &, const Number<J> &, const Number<K> &,
                  const Number<L> &, const Number<M> &, const Number<1> &) {
    if (M != 1)
      return (tS(M - 1, 0) + tS(0, M - 1)) *
                 r.eval(typename E::NumberDim(), Number<M - 1>(), Number<0>(),
                        Number<I - 1>(), Number<J - 1>(), Number<K - 1>(),
                        Number<L - 1>())

             +

             add(Number<I>(), Number<J>(), Number<K>(), Number<L>(),
                 Number<M - 1>(), Number<M - 1>());

    else
      return tS(0, 0) * r.eval(typename E::NumberDim(), Number<M - 1>(),
                               Number<0>(), Number<I - 1>(), Number<J - 1>(),
                               Number<K - 1>(), Number<L - 1>())

             +

             add(Number<I>(), Number<J>(), Number<K>(), Number<L>(),
                 Number<M - 1>(), Number<M - 1>());
  }

  template <int I, int J, int K, int L>
  inline auto add(const Number<I> &, const Number<J> &, const Number<K> &,
                  const Number<L> &, const Number<1> &, const Number<1> &) {
    return tS(0, 0) * r.eval(typename E::NumberDim(), Number<0>(), Number<0>(),
                             Number<I - 1>(), Number<J - 1>(), Number<K - 1>(),
                             Number<L - 1>());
  }

  template <int I, int J, int K, int L>
  inline void set(const Number<I> &, const Number<J> &, const Number<K> &,
                  const Number<L> &) {
    set(Number<I>(), Number<J>(), Number<K>(), Number<L - 1>());
    tA(I - 1, J - 1, K - 1, L - 1) =
        add(Number<I>(), Number<J>(), Number<K>(), Number<L>(),
            typename E::NumberDim(), typename E::NumberDim());
    if (K != I || L != J)
      tA(K - 1, L - 1, I - 1, J - 1) = tA(I - 1, J - 1, K - 1, L - 1);
  }

  template <int I, int J, int K>
  inline void set(const Number<I> &, const Number<J> &, const Number<K> &,
                  const Number<0> &) {
    set(Number<I>(), Number<J>(), Number<K - 1>(), Number<K - 1>());
  }

  template <int I, int J>
  inline void set(const Number<I> &, const Number<J> &, const Number<0> &,
                  const Number<0> &) {
    set(Number<I>(), Number<J - 1>(), Number<I>(), Number<J>());
  }

  template <int I, int K, int L>
  inline void set(const Number<I> &, const Number<0> &, const Number<K> &,
                  const Number<L> &) {
    set(Number<I - 1>(), Number<I - 1>(), Number<K>(), Number<L>());
  }

  template <int K, int L>
  inline void set(const Number<0> &, const Number<0> &, const Number<K> &,
                  const Number<L> &) {}
};

template <typename T1, typename T2, int NB, int Dim> struct EigenMatrixImp {

  using Val = const FTensor::Tensor1<T1, Dim>;
  using Vec = const FTensor::Tensor2<T2, Dim, Dim>;
  using Fun = boost::function<double(const double)>;
  using V = typename FTensor::promote<T1, T2>::V;

  template <int N> using Number = FTensor::Number<N>;
  template <char c> using I = typename FTensor::Index<c, Dim>;

  using NumberNb = Number<NB>;
  using NumberDim = Number<Dim>;

  static constexpr int nB = NB;

  EigenMatrixImp(Val &t_val, Vec &t_vec) : tVal(t_val), tVec(t_vec) {

    for (auto aa = 0; aa != Dim; ++aa)
      for (auto ii = 0; ii != Dim; ++ii)
        for (auto jj = 0; jj <= ii; ++jj)
          aM(aa, ii, jj) = tVec(aa, ii) * tVec(aa, jj);

    for (auto aa = 0; aa != Dim; ++aa) {
      for (auto bb = 0; bb != Dim; ++bb) {
        auto &MM = aMM[aa][bb];
        for (int ii = 0; ii != Dim; ++ii) {
          for (int jj = ii; jj != Dim; ++jj) {
            for (int kk = 0; kk != Dim; ++kk) {
              for (int ll = kk; ll != Dim; ++ll) {
                MM(ii, jj, kk, ll) = aM(aa, ii, jj) * aM(bb, kk, ll);
              }
            }
          }
        }
      }
    }

    for (auto aa = 0; aa != Dim; ++aa) {
      for (auto bb = 0; bb != Dim; ++bb) {
        if (aa != bb) {
          auto &MM = aMM[aa][bb];
          auto &G = aG[aa][bb];
          for (int ii = 0; ii != Dim; ++ii) {
            for (int jj = ii; jj != Dim; ++jj) {
              for (int kk = 0; kk != Dim; ++kk) {
                for (int ll = kk; ll != Dim; ++ll) {
                  G(ii, jj, kk, ll) = MM(ii, kk, jj, ll) + MM(ii, ll, jj, kk);
                }
              }
            }
          }
        }
      }
    }

    for (auto aa = 0; aa != Dim; ++aa) {
      for (auto bb = 0; bb != Dim; ++bb) {
        if (aa != bb) {
          auto &Gab = aG[aa][bb];
          auto &Gba = aG[bb][aa];
          auto &S = aS[aa][bb];
          for (int ii = 0; ii != Dim; ++ii) {
            for (int jj = ii; jj != Dim; ++jj) {
              for (int kk = 0; kk != Dim; ++kk) {
                for (int ll = kk; ll != Dim; ++ll) {
                  S(ii, jj, kk, ll) = Gab(ii, jj, kk, ll) + Gba(ii, jj, kk, ll);
                }
              }
            }
          }
        }
      }
    }



  }

  /**
   * @brief Get matrix
   *
   * \f[
   * \mathbf{B} = f(\mathbf{A})
   * \f]
   *
   * \f[
   * B_{ij} = \sum_{a}^3 f(\lambda^a) n^a_i n^a_j
   * \f]
   * where \f$a\f$ is eigen value number.
   *
   * @param t_val eiegn values vector
   * @param t_vec eigen vectors matrix
   * @param f function
   * @return auto function symmetric tensor rank two
   */
  inline auto getMat(Fun f) {

    for (auto aa = 0; aa != Dim; ++aa)
      fVal(aa) = f(tVal(aa));

    using T3 =
        FTensor::Tensor2_symmetric<typename std::remove_const<V>::type, Dim>;
    T3 t_A;
    GetMatImpl<EigenMatrixImp<T1, T2, NB, Dim>, V, T3>(*this, t_A)
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
  inline auto getDiffMat(Fun f, Fun d_f) {

    for (auto aa = 0; aa != Dim; ++aa)
      fVal(aa) = f(tVal(aa));

    for (auto aa = 0; aa != Dim; ++aa)
      dfVal(aa) = d_f(tVal(aa));

    for (auto aa = 0; aa != Dim; ++aa)
      for (auto bb = 0; bb != aa; ++bb) {
        aF(aa, bb) = 1 / (tVal(aa) - tVal(bb));
        aF(bb, aa) = -aF(aa, bb);
        aF2(aa, bb) = aF2(bb, aa) = aF(aa, bb) * aF(aa, bb);
      }

    using T3 = FTensor::Ddg<V, Dim, Dim>;
    T3 t_diff_A;
    GetDiffMatImpl<EigenMatrixImp<T1, T2, NB, Dim>, V, T3>(*this, t_diff_A)
        .set(Number<Dim>(), Number<Dim>(), Number<Dim>(), Number<Dim>());
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
  inline auto getDiffDiffMat(Fun f, Fun d_f, Fun dd_f, T &t_S) {

    for (auto aa = 0; aa != Dim; ++aa)
      fVal(aa) = f(tVal(aa));

    for (auto aa = 0; aa != Dim; ++aa)
      dfVal(aa) = d_f(tVal(aa));

    for (auto aa = 0; aa != Dim; ++aa)
      ddfVal(aa) = dd_f(tVal(aa));

    for (auto aa = 0; aa != Dim; ++aa)
      for (auto bb = 0; bb != aa; ++bb) {
        aF(aa, bb) = 1 / (tVal(aa) - tVal(bb));
        aF(bb, aa) = -aF(aa, bb);
        aF2(aa, bb) = aF2(bb, aa) = aF(aa, bb) * aF(aa, bb);
      }


    for (auto aa = 0; aa != Dim; ++aa) {
      for (auto bb = 0; bb != Dim; ++bb) {
        if (aa != bb)
          coefficientsType0(aa, bb) = dfVal(aa) * aF(aa, bb);
      }
    }

    if (NB == 3)
      for (auto aa = 0; aa != Dim; ++aa) {
        for (auto bb = 0; bb != Dim; ++bb) {
          if (aa != bb)
            for (auto cc = 0; cc != Dim; ++cc) {
              for (auto dd = 0; dd != Dim; ++dd) {
                if (cc != dd)
                  coefficientsType1(aa, bb, cc, dd) =
                      fVal(cc) * aF(cc, dd) * aF(aa, bb);
              }
            }
        }
      }

    if (NB == 2)
      for (auto aa = 0; aa != Dim; ++aa) {
        for (auto bb = 0; bb != Dim; ++bb) {
          if (aa != bb)
            for (auto cc = 0; cc != Dim; ++cc) {
              for (auto dd = 0; dd != Dim; ++dd)
                if (cc != dd) {

                  auto &r = coefficientsType1(aa, bb, cc, dd);

                  if ((cc == 1 || dd == 1) && (aa == 1 || bb == 1))
                    r = fVal(cc) * aF(cc, dd) * aF(aa, bb);
                  else if (cc != 1 && dd != 1 && aa != 1 && bb != 1) {

                    if ((aa != bb && bb != dd) && (aa != dd && bb != cc))
                      r = ddfVal(cc) / 4;
                    else
                      r = 0;

                  } else if ((cc != 1 && dd != 1) && (aa == 1 || bb == 1))
                    r = dfVal(cc) * aF(aa, bb) / 2;
                  else if ((cc == 1 || dd == 1) && (aa != 1 && bb != 1)) {

                    if ((cc == 2 && dd == 1) || (cc == 2 && dd == 1))
                      r = (

                              dfVal(cc)

                              - (fVal(cc) - fVal(dd)) * aF(cc, dd)

                                  ) *
                          aF(cc, dd);

                    else
                      r = 0;

                  } else
                    r = 0;
                }
            }
        }
      }

    if (NB == 1)
      for (auto aa = 0; aa != Dim; ++aa) {
        for (auto bb = 0; bb != Dim; ++bb) {
          if (aa != bb) {
            for (auto cc = 0; cc != Dim; ++cc) {
              for (auto dd = 0; dd != Dim; ++dd) {
                if (cc != dd) {
                  auto &r = coefficientsType1(aa, bb, cc, dd);
                  if ((bb != dd) && (aa != dd && bb != cc))
                    r = ddfVal(cc) / 4;
                  else
                    r = 0;
                }
              }
            }
          }
        }
      }

    using THIS = EigenMatrixImp<T1, T2, NB, Dim>;
    using T3 = FTensor::Ddg<V, Dim, Dim>;
    using CT1 = dd4MCoefficientsType1<THIS, V>;

    T3 t_diff_A;
    GetDiffDiffMatImpl<THIS, V, T3, T>(*this, t_diff_A, t_S)
        .set(Number<Dim>(), Number<Dim>(), Number<Dim>(), Number<Dim>());
    return t_diff_A;
  }

private:
  Val &tVal;
  Vec &tVec;
  FTensor::Christof<T2, Dim, Dim> aM;
  std::array<std::array<FTensor::Ddg<V, Dim, Dim>, Dim>, Dim> aMM;
  std::array<std::array<FTensor::Ddg<V, Dim, Dim>, Dim>, Dim> aG;
  std::array<std::array<FTensor::Ddg<V, Dim, Dim>, Dim>, Dim> aS;
  FTensor::Tensor2<T1, Dim, Dim> aF;
  FTensor::Tensor2<T1, Dim, Dim> aF2;
  FTensor::Tensor1<T1, Dim> fVal;
  FTensor::Tensor1<T1, Dim> dfVal;
  FTensor::Tensor1<T1, Dim> ddfVal;
  FTensor::Tensor2<T1, Dim, Dim> coefficientsType0;
  FTensor::Tensor4<T1, Dim, Dim, Dim, Dim> coefficientsType1;

  template <typename E, typename C> friend struct d2MCoefficients;
  template <typename E, typename C> friend struct d2MCoefficientsType0;
  template <typename E, typename C> friend struct dd4MCoefficientsType1;
  template <typename E, typename C> friend struct dd4MCoefficientsType2;
  template <typename E, typename C, typename G> friend struct d2MImpl;
  template <typename E, typename C, typename G1, typename G2>
  friend struct Fdd4MImpl;
  template <typename E, typename C> friend struct ReconstructMatImpl;
  template <typename E, typename C> friend struct FirstMatrixDirectiveImpl;
  template <typename E, typename C> friend struct SecondMatrixDirectiveImpl;
  template <typename E, typename C, typename T> friend struct GetDiffMatImpl;
  template <typename E, typename C, typename T3, typename T4>
  friend struct GetDiffDiffMatImpl;

  template <int a, int b, int i, int j>
  inline auto dFdN(const Number<a> &, const Number<b> &, const Number<i> &,
                   const Number<j> &) {
    return dFdN<a, b, i, j>();
  }

  template <int a, int b, int i, int j> inline auto dFdN() {
    return (-aM(a, i, j) + aM(b, i, j)) * aF2(a, b);
  }

  template <int a, int b, int i, int j, int k, int l>
  inline auto G(const Number<a> &, const Number<b> &, const Number<i> &,
                const Number<j> &, const Number<k> &, const Number<l> &) {
    return G<a, b, i, j, k, l>();
  }

  template <int a, int b, int i, int j, int k, int l> inline auto G() {
    return aG[a][b](i, j, k, l);
  }

  template <int a, int b, int i, int j, int k, int l>
  inline auto S(const Number<a> &, const Number<b> &, const Number<i> &,
                const Number<j> &, const Number<k> &, const Number<l> &) {
    return S<a, b, i, j, k, l>();
  }

  template <int a, int b, int i, int j, int k, int l> inline auto S() {
    return aS[a][b](i, j, k, l);
  }
}; // namespace EigenMatrix
} // namespace EigenMatrix