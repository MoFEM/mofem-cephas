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
  inline auto get(const Number<a> &, const Number<b> &, const int i,
                  const int j, const int k, const int l,
                  const Number<3> &) const {
    return e.aS[a][b](i, j, k, l) * e.fVal(a) * e.aF(a, b);
  }

  template <int a, int b>
  inline auto get(const Number<a> &, const Number<b> &, const int i,
                  const int j, const int k, const int l,
                  const Number<2> &) const {
    if (a == 1 || b == 1)
      return get(Number<a>(), Number<b>(), i, j, k, l, Number<3>());
    else
      return get(Number<a>(), Number<b>(), i, j, k, l, Number<1>());
  }

  template <int a, int b>
  inline auto get(const Number<a> &, const Number<b> &, const int i,
                  const int j, const int k, const int l,
                  const Number<1>) const {
    return e.aS[a][b](i, j, k, l) * e.dfVal(a) / static_cast<C>(2);
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

  template <int b, int a>
  inline C term(const int i, const int j, const int k, const int l) const {
    if (a != b) {
      return g.get(Number<a>(), Number<b>(), i, j, k, l,
                   typename E::NumberNb());
    }
    return 0;
  }

  template <int nb, int a>
  inline C eval(const Number<nb> &, const Number<a> &, const int i, const int j,
                const int k, const int l) const {
    return term<nb - 1, a>(i, j, k, l) +
           eval(Number<nb - 1>(), Number<a>(), i, j, k, l);
  }

  template <int a>
  inline C eval(const Number<1> &, const Number<a> &, const int i,
                const int j, const int k, const int l) const {
    return term<0, a>(i, j, k, l);
  }
};

template <typename E, typename C> struct Fdd4MImpl {
  using Val = typename E::Val;
  using Vec = typename E::Vec;
  using Fun = typename E::Fun;

  Fdd4MImpl(E &e) : e(e) {}
  E &e;

  template <int N> using Number = FTensor::Number<N>;

  template <int A, int a, int b, int I, int J, int K, int L, int M, int N>
  inline auto fd2M() const {
    return e.d2MType1[A][a][b](I, J, K, L);
  }

  template <int a, int b, int A, int B, int I, int J, int K, int L, int M,
            int N>
  inline auto fd2G() const {
    return fd2M<A, a, b, I, K, N, M, J, L>() * e.aM[B](J, L) +
           e.aM[A](I, K) * fd2M<B, a, b, J, L, M, N, I, K>() +
           fd2M<A, a, b, I, L, M, N, J, K>() * e.aM[B](J, K) +
           e.aM[A](I, L) * fd2M<B, a, b, J, K, M, N, I, L>();
  }

  template <int a, int b, int I, int J, int K, int L, int M, int N>
  inline auto fd2S() const {
    return fd2G<a, b, a, b, I, J, K, L, M, N>() +
           fd2G<a, b, b, a, I, J, K, L, M, N>();
  }

  template <int a, int b, int i, int j, int k, int l, int m, int n>
  inline C term() const {

    if (a != b) {

      if (e.nB == 1) {
        return

            fd2S<a, b, i, j, k, l, m, n>();

      } else {
        return

            fd2S<a, b, i, j, k, l, m, n>()

            +

            e.d2MType2[a][b][m][n](i, j, k, l);
      }
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

  template <int a> inline C term(const int ii, const int jj) const {
    return e.aM[a](ii, jj) * e.fVal(a);
  }

  template <int nb>
  inline C eval(const Number<nb> &, const int ii, const int jj) const {
    return term<nb - 1>(ii, jj) + eval(Number<nb - 1>(), ii, jj);
  }

  inline C eval(const Number<1> &, const int ii, const int jj) const {
    return term<0>(ii, jj);
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

  template <int a>
  inline C term(const int i, const int j, const int k, const int l) const {
    return

        e.aMM[a][a](i, j, k, l) * e.dfVal(a)

        +

        r.eval(typename E::NumberDim(), Number<a>(), i, j, k, l) /
            static_cast<C>(2);
  }

  template <int nb>
  inline C eval(const Number<nb> &, const int i, const int j, const int k,
                const int l) const {
    return term<nb - 1>(i, j, k, l) + eval(Number<nb - 1>(), i, j, k, l);
  }

  inline C eval(const Number<1> &, const int i, const int j, const int k,
                const int l) const {
    return term<0>(i, j, k, l);
  }
};

template <typename E, typename C> struct SecondMatrixDirectiveImpl {
  using Val = typename E::Val;
  using Vec = typename E::Vec;
  using Fun = typename E::Fun;

  template <int N> using Number = FTensor::Number<N>;

  SecondMatrixDirectiveImpl(E &e) : r(e), e(e) {}
  Fdd4MImpl<E, C> r;
  E &e;

  template <int a, int i, int j, int k, int l, int m, int n>
  inline C term() const {

    return

        (

            e.d2MType0[a][k][l](i, j, m, n)

            +

            e.d2MType0[a][i][j](k, l, m, n)

            +

            e.d2MType0[a][m][n](i, j, k, l)

                ) /
            static_cast<C>(2) +

        e.aMM[a][a](i, j, k, l) * e.aM[a](m, n) * e.ddfVal(a)

        +

        r.eval(typename E::NumberDim(), Number<a>(), Number<i>(), Number<j>(),
               Number<k>(), Number<l>(), Number<m>(), Number<n>()) /
            static_cast<C>(4);
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

  template <int Dim> inline void set(const Number<Dim> &) {
    for (int ii = 0; ii != Dim; ++ii)
      for (int jj = ii; jj != Dim; ++jj)
        tA(ii, jj) = r.eval(typename E::NumberDim(), ii, jj);
  }
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
    return (tS(M - 1, 0) + tS(0, M - 1)) *
               r.eval(typename E::NumberDim(), Number<M - 1>(), Number<0>(),
                      Number<I - 1>(), Number<J - 1>(), Number<K - 1>(),
                      Number<L - 1>())

           +

           add(Number<I>(), Number<J>(), Number<K>(), Number<L>(),
               Number<M - 1>(), Number<M - 1>());

    ;
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
  using V = double; //typename FTensor::promote<T1, T2>::V;

  template <int N> using Number = FTensor::Number<N>;
  template <char c> using I = typename FTensor::Index<c, Dim>;

  using NumberNb = Number<NB>;
  using NumberDim = Number<Dim>;

  static constexpr int nB = NB;

  EigenMatrixImp(Val &t_val, Vec &t_vec) : tVal(t_val), tVec(t_vec) {

    for (auto aa = 0; aa != Dim; ++aa) {
      auto &M = aM[aa];
      for (auto ii = 0; ii != Dim; ++ii)
        for (auto jj = 0; jj <= ii; ++jj)
          M(ii, jj) = tVec(aa, ii) * tVec(aa, jj);
    }

    FTensor::Index<'i', Dim> i;
    FTensor::Index<'j', Dim> j;
    FTensor::Index<'k', Dim> k;
    FTensor::Index<'l', Dim> l;

    for (auto aa = 0; aa != Dim; ++aa) {
      for (auto bb = 0; bb != Dim; ++bb) {
        auto &Ma = aM[aa];
        auto &Mb = aM[bb];
        auto &MM = aMM[aa][bb];
        MM(i, j, k, l) = Ma(i, j) * Mb(k, l);
      }
    }

    for (auto aa = 0; aa != Dim; ++aa) {
      for (auto bb = 0; bb != Dim; ++bb) {
        if (aa != bb) {
          auto &MM = aMM[aa][bb];
          auto &G = aG[aa][bb];
          G(i, j, k, l) = MM(i, k, j, l) || MM(i, l, j, k);
        }
      }
    }

    for (auto aa = 0; aa != Dim; ++aa) {
      for (auto bb = 0; bb != Dim; ++bb) {
        if (aa != bb) {
          auto &Gab = aG[aa][bb];
          auto &Gba = aG[bb][aa];
          auto &S = aS[aa][bb];
          S(i, j, k, l) = Gab(i, j, k, l) + Gba(i, j, k, l);
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
        .set(Number<Dim>());
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

    FTensor::Index<'i', Dim> i;
    FTensor::Index<'j', Dim> j;
    FTensor::Index<'k', Dim> k;
    FTensor::Index<'l', Dim> l;

    for (auto aa = 0; aa != Dim; ++aa) {
      for (auto bb = 0; bb != Dim; ++bb) {
        if (aa != bb) {
          auto &S = aS[aa][bb];
          auto &M = aM[aa];
          for (auto mm = 0; mm != Dim; ++mm) {
            for (auto nn = mm; nn != Dim; ++nn) {
              auto &SMmn = aSM[aa][bb][mm][nn];
              auto &SMnm = aSM[aa][bb][nn][mm];
              SMmn(i, j, k, l) = aS[aa][bb](i, j, k, l) * M(mm, nn);
              SMnm(i, j, k, l) = SMmn(i, j, k, l);
            }
          }
        }
      }
    }

    for (auto aa = 0; aa != Dim; ++aa) {
      for (auto mm = 0; mm != Dim; ++mm) {
        for (auto nn = mm; nn != Dim; ++nn) {
          d2MType0[aa][nn][mm](i, j, k, l) = 0;
          d2MType0[aa][mm][nn](i, j, k, l) = 0;
        }
      }
    }

    if (NB == 3)
      for (auto aa = 0; aa != Dim; ++aa) {
        for (auto bb = 0; bb != Dim; ++bb) {
          if (aa != bb) {
            const V v = dfVal(aa) * aF(aa, bb);
            for (auto mm = 0; mm != Dim; ++mm) {
              for (auto nn = mm; nn != Dim; ++nn) {
                d2MType0[aa][mm][nn](i, j, k, l) +=
                    v * aSM[aa][bb][mm][nn](i, j, k, l);
                d2MType0[aa][nn][mm](i, j, k, l) =
                    d2MType0[aa][mm][nn](i, j, k, l);
              }
            }
          }
        }
      }

    if (NB == 2)
      for (auto aa = 0; aa != Dim; ++aa) {
        for (auto bb = 0; bb != Dim; ++bb) {
          if (aa != bb) {
            V v;
            if (aa == 1 || bb == 1)
              v = dfVal(aa) * aF(aa, bb);
            else
              v = ddfVal(aa) / 2;
            for (auto mm = 0; mm != Dim; ++mm) {
              for (auto nn = mm; nn != Dim; ++nn) {
                d2MType0[aa][mm][nn](i, j, k, l) +=
                    v * aSM[aa][bb][mm][nn](i, j, k, l);
                d2MType0[aa][nn][mm](i, j, k, l) =
                    d2MType0[aa][mm][nn](i, j, k, l);
              }
            }
          }
        }
      }

    if (NB == 1)
      for (auto aa = 0; aa != Dim; ++aa) {
        for (auto bb = 0; bb != Dim; ++bb) {
          if (aa != bb) {
            const V v = ddfVal(aa) / 2;
            for (auto mm = 0; mm != Dim; ++mm) {
              for (auto nn = mm; nn != Dim; ++nn) {
                d2MType0[aa][mm][nn](i, j, k, l) +=
                    v * aSM[aa][bb][mm][nn](i, j, k, l);
                d2MType0[aa][nn][mm](i, j, k, l) =
                    d2MType0[aa][mm][nn](i, j, k, l);
              }
            }
          }
        }
      }

    for (auto aa = 0; aa != Dim; ++aa) {
      for (auto mm = 0; mm != Dim; ++mm) {
        for (auto nn = mm; nn != Dim; ++nn) {
          d2MType1[aa][nn][mm](i, j, k, l) = 0;
          d2MType1[aa][mm][nn](i, j, k, l) = 0;
        }
      }
    }

    if (NB == 3)
      for (auto aa = 0; aa != Dim; ++aa) {
        for (auto bb = 0; bb != Dim; ++bb) {
          if (aa != bb) {
            const auto &S = aS[aa][bb];
            const auto v0 = aF(aa, bb);
            for (auto cc = 0; cc != Dim; ++cc) {
                for (auto dd = 0; dd != Dim; ++dd) {
                  if (cc != dd) {
                    const double v1 = fVal(cc) * aF(cc, dd);
                    d2MType1[aa][cc][dd](i, j, k, l) +=
                        (v1 * v0) * S(i, j, k, l);
                  }
                }
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

                  V r;

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

                  d2MType1[aa][cc][dd](i, j, k, l) +=
                      r * aS[aa][bb](i, j, k, l);
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
                  V r;
                  if ((bb != dd) && (aa != dd && bb != cc)) {
                    r = ddfVal(cc) / 4;
                    d2MType1[aa][cc][dd](i, j, k, l) +=
                        r * aS[aa][bb](i, j, k, l);
                  }
                }
              }
            }
          }
        }
      }

    if (NB == 3)
      for (auto aa = 0; aa != Dim; ++aa) {
        for (auto bb = 0; bb != Dim; ++bb) {
          if (aa != bb) {
            for (auto mm = 0; mm != Dim; ++mm) {
              for (auto nn = mm; nn != Dim; ++nn) {
                d2MType2[aa][bb][mm][nn](i, j, k, l) =
                    2 * (fVal(aa) * aF2(aa, bb)) *
                    (-aSM[aa][bb][mm][nn](i, j, k, l) +
                     aSM[bb][aa][mm][nn](i, j, k, l));
                d2MType2[aa][bb][nn][mm](i, j, k, l) =
                    d2MType2[aa][bb][mm][nn](i, j, k, l);
              }
            }
          }
        }
      }

    if (NB == 2)
      for (auto aa = 0; aa != Dim; ++aa) {
        for (auto bb = 0; bb != Dim; ++bb) {
          if (aa != bb) {
            for (auto mm = 0; mm != Dim; ++mm) {
              for (auto nn = mm; nn != Dim; ++nn) {
                if (aa == 1 || bb == 1) {
                  d2MType2[aa][bb][mm][nn](i, j, k, l) =
                      2 * (fVal(aa) * aF2(aa, bb)) *
                      (-aSM[aa][bb][mm][nn](i, j, k, l) +
                       aSM[bb][aa][mm][nn](i, j, k, l));
                  d2MType2[aa][bb][nn][mm](i, j, k, l) =
                      d2MType2[aa][bb][mm][nn](i, j, k, l);
                } else {
                  d2MType2[aa][bb][mm][nn](i, j, k, l) = 0;
                  d2MType2[aa][bb][nn][mm](i, j, k, l) = 0;
                }
              }
            }
          }
        }
      }

    using THIS = EigenMatrixImp<T1, T2, NB, Dim>;
    using T3 = FTensor::Ddg<V, Dim, Dim>;

    T3 t_diff_A;
    GetDiffDiffMatImpl<THIS, V, T3, T>(*this, t_diff_A, t_S)
        .set(Number<Dim>(), Number<Dim>(), Number<Dim>(), Number<Dim>());
    return t_diff_A;
  }

private:
  Val &tVal;
  Vec &tVec;
  FTensor::Tensor2_symmetric<V, Dim> aM[Dim];
  FTensor::Ddg<V, Dim, Dim> aMM[Dim][Dim];
  FTensor::Ddg<V, Dim, Dim> aG[Dim][Dim];
  FTensor::Ddg<V, Dim, Dim> aS[Dim][Dim];
  FTensor::Ddg<V, Dim, Dim> aSM[Dim][Dim][Dim][Dim];
  FTensor::Ddg<V, Dim, Dim> d2MType2[Dim][Dim][Dim][Dim];
  FTensor::Ddg<V, Dim, Dim> d2MType0[Dim][Dim][Dim];
  FTensor::Ddg<V, Dim, Dim> d2MType1[Dim][Dim][Dim];
  FTensor::Tensor2<V, Dim, Dim> aF;
  FTensor::Tensor2<V, Dim, Dim> aF2;
  FTensor::Tensor1<V, Dim> fVal;
  FTensor::Tensor1<V, Dim> dfVal;
  FTensor::Tensor1<V, Dim> ddfVal;

  template <typename E, typename C> friend struct d2MCoefficients;
  template <typename E, typename C, typename G> friend struct d2MImpl;
  template <typename E, typename C> friend struct Fdd4MImpl;
  template <typename E, typename C> friend struct ReconstructMatImpl;
  template <typename E, typename C> friend struct FirstMatrixDirectiveImpl;
  template <typename E, typename C> friend struct SecondMatrixDirectiveImpl;
  template <typename E, typename C, typename T> friend struct GetDiffMatImpl;
  template <typename E, typename C, typename T3, typename T4>
  friend struct GetDiffDiffMatImpl;

}; // namespace EigenMatrix
} // namespace EigenMatrix