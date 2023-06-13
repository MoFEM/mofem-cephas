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
 Return t_A is symmetric tensor rank two.

 Calculate derivarive
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



namespace EigenMatrix {

template <int N> using Number = FTensor::Number<N>;

template <int N1, int N2, int Dim>
inline auto get_sym_index(const Number<N1> &, const Number<N2> &,
                   const Number<Dim> &) {
  if constexpr (N1 > N2)
    return N1 + (N2 * (2 * Dim - N2 - 1)) / 2;
  else
    return N2 + (N1 * (2 * Dim - N1 - 1)) / 2;
}

inline auto get_sym_index(const int N1, const int N2, const int Dim) {
  if (N1 > N2)
    return N1 + (N2 * (2 * Dim - N2 - 1)) / 2;
  else
    return N2 + (N1 * (2 * Dim - N1 - 1)) / 2;
}

template <int N1, int N2, int Dim>
inline auto get_nodiag_index(const Number<N1> &, const Number<N2> &,
                   const Number<Dim> &) {
  static_assert(N1 != N2, "Bad index");
  if constexpr (N2 > N1)
    return (Dim - 1) * N1 + N2 - 1;
  else
    return (Dim - 1) * N1 + N2;
}

inline auto get_nodiag_index(const int N1, const int N2, int Dim) {
  if (N2 > N1)
    return (Dim - 1) * N1 + N2 - 1;
  else
    return (Dim - 1) * N1 + N2;
}

template <typename E, typename C> struct d2MCoefficients {
  using Val = typename E::Val;
  using Vec = typename E::Vec;
  using Fun = typename E::Fun;

  using NumberNb = typename E::NumberNb;
  using NumberDim = typename E::NumberDim;

  template <int N> using Number = FTensor::Number<N>;

  d2MCoefficients(E &e) : e(e) {}
  E &e;

  template <int a, int b>
  inline auto get(const Number<a> &, const Number<b> &, const int i,
                  const int j, const int k, const int l,
                  const Number<3> &) const {
    return e.aS[get_sym_index(Number<a>(), Number<b>(), NumberDim())](i, j, k,
                                                                      l) *
           e.fVal(a) * e.aF(a, b);
  }

  template <int a, int b>
  inline auto get(const Number<a> &, const Number<b> &, const int i,
                  const int j, const int k, const int l,
                  const Number<2> &) const {
    if constexpr (a == 1 || b == 1)
      return get(Number<a>(), Number<b>(), i, j, k, l, Number<3>());
    else
      return get(Number<a>(), Number<b>(), i, j, k, l, Number<1>());
  }

  template <int a, int b>
  inline auto get(const Number<a> &, const Number<b> &, const int i,
                  const int j, const int k, const int l,
                  const Number<1>) const {
    return e.aS[get_sym_index(Number<a>(), Number<b>(), NumberDim())](i, j, k,
                                                                      l) *
           e.dfVal(a) / static_cast<C>(2);
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
    if constexpr (a != b) {
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
  inline C eval(const Number<1> &, const Number<a> &, const int i, const int j,
                const int k, const int l) const {
    return term<0, a>(i, j, k, l);
  }
};

template <typename E, typename C> struct Fdd4MImpl {
  using Val = typename E::Val;
  using Vec = typename E::Vec;
  using Fun = typename E::Fun;

  Fdd4MImpl(E &e) : e(e) {}
  E &e;

  using NumberNb = typename E::NumberNb;
  using NumberDim = typename E::NumberDim;


  template <int A, int a, int b, int I, int J, int K, int L>
  inline auto fd2M() const {
    return e.d2MType1[b][get_sym_index(Number<A>(), Number<a>(), NumberDim())](
        Number<I>(), Number<J>(), Number<K>(), Number<L>());
  }

  template <int a, int b, int i, int j, int k, int l, int m, int n>
  inline auto term_fd2S(const Number<a> &, const Number<b> &, const Number<i> &,
                        const Number<j> &, const Number<k> &, const Number<l> &,
                        const Number<m> &, const Number<n> &) const {
    if constexpr (i == j && k == l) {
      return 4 *
             (

                 fd2M<a, a, b, i, k, m, n>() * e.aM[b](Number<j>(), Number<l>())

                 +

                 fd2M<b, a, b, i, k, m, n>() * e.aM[a](Number<j>(), Number<l>())

             );

    } else if constexpr (i == j)
      return 2 *
             (

                 fd2M<a, a, b, i, k, m, n>() *
                     e.aM[b](Number<j>(), Number<l>()) +
                 fd2M<b, a, b, j, l, m, n>() * e.aM[a](Number<i>(), Number<k>())

                 +

                 fd2M<b, a, b, i, k, m, n>() *
                     e.aM[a](Number<j>(), Number<l>()) +
                 fd2M<a, a, b, j, l, m, n>() * e.aM[b](Number<i>(), Number<k>())

             );
    else if constexpr (k == l)
      return 2 *
             (

                 fd2M<a, a, b, i, k, m, n>() *
                     e.aM[b](Number<j>(), Number<l>()) +
                 fd2M<b, a, b, j, l, m, n>() *
                     e.aM[a](Number<i>(), Number<k>()) +

                 +

                     fd2M<b, a, b, i, k, m, n>() *
                     e.aM[a](Number<j>(), Number<l>()) +
                 fd2M<a, a, b, j, l, m, n>() * e.aM[b](Number<i>(), Number<k>())

             );
    else
      return fd2M<a, a, b, i, k, m, n>() * e.aM[b](Number<j>(), Number<l>()) +
             fd2M<b, a, b, j, l, m, n>() * e.aM[a](Number<i>(), Number<k>()) +
             fd2M<a, a, b, i, l, m, n>() * e.aM[b](Number<j>(), Number<k>()) +
             fd2M<b, a, b, j, k, m, n>() * e.aM[a](Number<i>(), Number<l>())

             +

             fd2M<b, a, b, i, k, m, n>() * e.aM[a](Number<j>(), Number<l>()) +
             fd2M<a, a, b, j, l, m, n>() * e.aM[b](Number<i>(), Number<k>()) +
             fd2M<b, a, b, i, l, m, n>() * e.aM[a](Number<j>(), Number<k>()) +
             fd2M<a, a, b, j, k, m, n>() * e.aM[b](Number<i>(), Number<l>());
  }

  template <int NB, int a, int b, int i, int j, int k, int l, int m, int n>
  inline C term_SM(const Number<a> &, const Number<b> &, const Number<NB> &,
                   const Number<i> &, const Number<j> &, const Number<k> &,
                   const Number<l> &, const Number<m> &,
                   const Number<n> &) const {

    if constexpr (NB == 1) {
      return 0;

    } else if constexpr (NB == 2) {

      if constexpr (a == 1 || b == 1) {
        return

            e.aF2(Number<a>(), Number<b>()) *
            (e.aSM[get_nodiag_index(Number<b>(), Number<a>(), NumberDim())]
                  [get_sym_index(Number<m>(), Number<n>(), NumberDim())](
                      Number<i>(), Number<j>(), Number<k>(), Number<l>()) -
             e.aSM[get_nodiag_index(Number<a>(), Number<b>(), NumberDim())]
                  [get_sym_index(Number<m>(), Number<n>(), NumberDim())](
                      Number<i>(), Number<j>(), Number<k>(), Number<l>()));

      } else {
        return 0;
      }

    } else {

      return

          e.aF2(Number<a>(), Number<b>()) *
          (e.aSM[get_nodiag_index(Number<b>(), Number<a>(), NumberDim())]
                [get_sym_index(Number<m>(), Number<n>(), NumberDim())](
                    Number<i>(), Number<j>(), Number<k>(), Number<l>()) -
           e.aSM[get_nodiag_index(Number<a>(), Number<b>(), NumberDim())]
                [get_sym_index(Number<m>(), Number<n>(), NumberDim())](
                    Number<i>(), Number<j>(), Number<k>(), Number<l>()));
    }

    return 0;
  }

  template <int nb, int a, int i, int j, int k, int l, int m, int n>
  inline C eval_fdS2(const Number<nb> &, const Number<a> &, const Number<i> &,
                     const Number<j> &, const Number<k> &, const Number<l> &,
                     const Number<m> &, const Number<n> &) const {
    if constexpr (a != nb - 1)
      return term_fd2S(Number<a>(), Number<nb - 1>(), Number<i>(), Number<j>(),
                       Number<k>(), Number<l>(), Number<m>(), Number<n>()) +
             eval_fdS2(Number<nb - 1>(), Number<a>(), Number<i>(), Number<j>(),
                       Number<k>(), Number<l>(), Number<m>(), Number<n>());
    else
      return eval_fdS2(Number<nb - 1>(), Number<a>(), Number<i>(), Number<j>(),
                       Number<k>(), Number<l>(), Number<m>(), Number<n>());
  }

  template <int a, int i, int j, int k, int l, int m, int n>
  inline C eval_fdS2(const Number<1> &, const Number<a> &, const Number<i> &,
                     const Number<j> &, const Number<k> &, const Number<l> &,
                     const Number<m> &, const Number<n> &) const {
    if constexpr (a != 0)
      return term_fd2S(Number<a>(), Number<0>(), Number<i>(), Number<j>(),
                       Number<k>(), Number<l>(), Number<m>(), Number<n>());
    else
      return 0;
  }

  template <int nb, int a, int i, int j, int k, int l, int m, int n>
  inline C eval_SM(const Number<nb> &, const Number<a> &, const Number<i> &,
                   const Number<j> &, const Number<k> &, const Number<l> &,
                   const Number<m> &, const Number<n> &) const {
    if constexpr (a != nb - 1)
      return term_SM(Number<a>(), Number<nb - 1>(), NumberNb(), Number<i>(),
                     Number<j>(), Number<k>(), Number<l>(), Number<m>(),
                     Number<n>()) +
             eval_SM(Number<nb - 1>(), Number<a>(), Number<i>(), Number<j>(),
                     Number<k>(), Number<l>(), Number<m>(), Number<n>());

    else
      return eval_SM(Number<nb - 1>(), Number<a>(), Number<i>(), Number<j>(),
                     Number<k>(), Number<l>(), Number<m>(), Number<n>());
  }

  template <int a, int i, int j, int k, int l, int m, int n>
  inline C eval_SM(const Number<1> &, const Number<a> &, const Number<i> &,
                   const Number<j> &, const Number<k> &, const Number<l> &,
                   const Number<m> &, const Number<n> &) const {
    if constexpr (a != 0)
      return term_SM(Number<a>(), Number<0>(), NumberNb(), Number<i>(),
                     Number<j>(), Number<k>(), Number<l>(), Number<m>(),
                     Number<n>());
    else
      return 0;
  }

  template <int nb, int a, int i, int j, int k, int l, int m, int n>
  inline C eval(const Number<nb> &, const Number<a> &, const Number<i> &,
                const Number<j> &, const Number<k> &, const Number<l> &,
                const Number<m> &, const Number<n> &) const {
    return

        (2 * e.fVal(a)) * eval_SM(NumberDim(), Number<a>(), Number<i>(),
                                  Number<j>(), Number<k>(), Number<l>(),
                                  Number<m>(), Number<n>())

        +

        eval_fdS2(NumberDim(), Number<a>(), Number<i>(), Number<j>(),
                  Number<k>(), Number<l>(), Number<m>(), Number<n>());
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
    return e.aM[a](Number<i>(), Number<j>()) * e.fVal(a);
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

        r.eval(typename E::NumberDim(), Number<a>(), i, j, k, l) * 0.5;
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

  using NumberDim = typename E::NumberDim;

  template <int N> using Number = FTensor::Number<N>;

  SecondMatrixDirectiveImpl(E &e) : r(e), e(e) {}
  Fdd4MImpl<E, C> r;
  E &e;

  template <int a, int i, int j, int k, int l, int m, int n>
  inline C term1() const {
    return e.d2MType0[a][get_sym_index(Number<k>(), Number<l>(), NumberDim())](
        Number<i>(), Number<j>(), Number<m>(), Number<n>());
  };

  template <int a, int i, int j, int k, int l, int m, int n>
  inline C term2() const {
    return e.d2MType0[a][get_sym_index(Number<i>(), Number<j>(), NumberDim())](
        Number<k>(), Number<l>(), Number<m>(), Number<n>());
  }

  template <int a, int i, int j, int k, int l, int m, int n>
  inline C term3() const {
    return e.d2MType0[a][get_sym_index(Number<n>(), Number<m>(), NumberDim())](
        Number<i>(), Number<j>(), Number<k>(), Number<l>());
  }

  template <int a, int i, int j, int k, int l, int m, int n>
  inline C term() const {

    return

        (term1<a, i, j, k, l, m, n>() + term2<a, i, j, k, l, m, n>() +
         term3<a, i, j, k, l, m, n>()) *
            0.5

        +

        (e.aMM[a][a](Number<i>(), Number<j>(), Number<k>(), Number<l>()) *
         e.aM[a](Number<m>(), Number<n>())) *
            e.ddfVal(a)

        +

        r.eval(typename E::NumberDim(), Number<a>(), Number<i>(), Number<j>(),
               Number<k>(), Number<l>(), Number<m>(), Number<n>()) *
            0.25;
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
                const Number<n> &) const {
    return term<0, i, j, k, l, m, n>();
  }
};

template <typename E, typename C, typename T> struct GetMatImpl {
  using Val = typename E::Val;
  using Vec = typename E::Vec;
  using Fun = typename E::Fun;

  using NumberNb = typename E::NumberNb;
  using NumberDim = typename E::NumberDim;

  template <int N> using Number = FTensor::Number<N>;

  GetMatImpl(E &e, T &t_a) : r(e), tA(t_a) {}
  ReconstructMatImpl<E, C> r;
  T &tA;

  template <int i, int j>
  inline void set(const Number<i> &, const Number<j> &) {

    set(Number<i>(), Number<j - 1>());

    tA(Number<i - 1>(), Number<j - 1>()) =
        r.eval(NumberDim(), Number<i - 1>(), Number<j - 1>());
  }

  template <int i> inline void set(const Number<i> &, const Number<1> &) {

    set(Number<i - 1>(), Number<i - 1>());

    tA(Number<i - 1>(), Number<0>()) =
        r.eval(NumberDim(), Number<i - 1>(), Number<0>());
  }

  inline void set(const Number<1> &, const Number<1> &) {
    tA(Number<0>(), Number<0>()) =
        r.eval(NumberDim(), Number<0>(), Number<0>());
  }
};

template <typename E, typename C, typename T> struct GetDiffMatImpl {
  using Val = typename E::Val;
  using Vec = typename E::Vec;
  using Fun = typename E::Fun;

  using NumberNb = typename E::NumberNb;
  using NumberDim = typename E::NumberDim;

  template <int N> using Number = FTensor::Number<N>;

  GetDiffMatImpl(E &e, T &t_a) : r(e), tA(t_a) {}
  FirstMatrixDirectiveImpl<E, C> r;
  T &tA;

  template <int i, int j, int k, int l>
  inline void set(const Number<i> &, const Number<j> &, const Number<k> &,
                  const Number<l> &) {

    tA(Number<i - 1>(), Number<j - 1>(), Number<k - 1>(), Number<l - 1>()) =
        r.eval(NumberDim(), Number<i - 1>(), Number<j - 1>(), Number<k - 1>(),
               Number<l - 1>());

    set(Number<i - 1>(), Number<j>(), Number<k>(), Number<l>());
  }

  template <int j, int k, int l>
  inline void set(const Number<1> &, const Number<j> &, const Number<k> &,
           const Number<l> &) {

    tA(Number<0>(), Number<j - 1>(), Number<k - 1>(), Number<l - 1>()) =
        r.eval(NumberDim(), Number<0>(), Number<j - 1>(), Number<k - 1>(),
               Number<l - 1>());

    set(Number<j - 1>(), Number<j - 1>(), Number<k>(), Number<l>());
  }

  template <int k, int l>
  inline void set(const Number<1> &, const Number<1> &, const Number<k> &,
           const Number<l> &) {

    tA(Number<0>(), Number<0>(), Number<k - 1>(), Number<l - 1>()) =
        r.eval(NumberDim(), Number<0>(), Number<0>(), Number<k - 1>(),
               Number<l - 1>());

    set(NumberDim(), NumberDim(), Number<k - 1>(), Number<l>());
  }

  template <int l>
  inline void set(const Number<1> &, const Number<1> &, const Number<1> &,
                  const Number<l> &) {

    tA(Number<0>(), Number<0>(), Number<0>(), Number<l - 1>()) =
        r.eval(NumberDim(), Number<0>(), Number<0>(), Number<0>(),
               Number<l - 1>());

    set(NumberDim(), NumberDim(), Number<l - 1>(), Number<l - 1>());
  }

  inline void set(const Number<1> &, const Number<1> &, const Number<1> &,
                  const Number<1> &) {

    tA(Number<0>(), Number<0>(), Number<0>(), Number<0>()) =
        r.eval(NumberDim(), Number<0>(), Number<0>(), Number<0>(), Number<0>());
  }
};

template <typename E, typename C, typename T1, typename T2>
struct GetDiffDiffMatImpl {
  using Val = typename E::Val;
  using Vec = typename E::Vec;
  using Fun = typename E::Fun;

  using NumberNb = typename E::NumberNb;
  using NumberDim = typename E::NumberDim;

  template <int N> using Number = FTensor::Number<N>;

  GetDiffDiffMatImpl(E &e, T1 &t_a, T2 &t_S) : r(e), e(e), tA(t_a), tS(t_S) {}
  SecondMatrixDirectiveImpl<E, C> r;
  E &e;
  T1 &tA;
  T2 &tS;

  template <int I, int J, int K, int L, int M, int N>
  inline auto add(const Number<I> &, const Number<J> &, const Number<K> &,
                  const Number<L> &, const Number<M> &, const Number<N> &) {
    if constexpr (N != M)
      return (tS(M - 1, N - 1) + tS(N - 1, M - 1)) *
                 r.eval(NumberDim(), Number<M - 1>(), Number<N - 1>(),
                        Number<I - 1>(), Number<J - 1>(), Number<K - 1>(),
                        Number<L - 1>())

             +

             add(Number<I>(), Number<J>(), Number<K>(), Number<L>(),
                 Number<M>(), Number<N - 1>());
    else
      return tS(M - 1, N - 1) * r.eval(NumberDim(), Number<M - 1>(),
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
               r.eval(NumberDim(), Number<M - 1>(), Number<0>(),
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
    return tS(0, 0) * r.eval(NumberDim(), Number<0>(), Number<0>(),
                             Number<I - 1>(), Number<J - 1>(), Number<K - 1>(),
                             Number<L - 1>());
  }

  template <int I, int J, int K, int L>
  inline void set(const Number<I> &, const Number<J> &, const Number<K> &,
                  const Number<L> &) {
    set(Number<I>(), Number<J>(), Number<K>(), Number<L - 1>());
    tA(I - 1, J - 1, K - 1, L - 1) = add(Number<I>(), Number<J>(), Number<K>(),
                                         Number<L>(), NumberDim(), NumberDim());
    // Major symmetry
    if constexpr (K != I || L != J)
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
    set(Number<I>(), Number<J - 1>(), Number<I>(), Number<J - 1>());
  }

  template <int I, int K>
  inline void set(const Number<I> &, const Number<0> &, const Number<K> &,
                  const Number<0> &) {
    set(Number<I - 1>(), Number<I - 1>(), Number<I - 1>(), Number<I - 1>());
  }

  inline void set(const Number<0> &, const Number<0> &, const Number<0> &,
                  const Number<0> &) {}
};

template <typename E, typename C, typename T1, typename VT2, int DimT2>
struct GetDiffDiffMatImpl<E, C, T1, FTensor::Tensor2_symmetric<VT2, DimT2>> {
  using Val = typename E::Val;
  using Vec = typename E::Vec;
  using Fun = typename E::Fun;

  using NumberNb = typename E::NumberNb;
  using NumberDim = typename E::NumberDim;

  template <int N> using Number = FTensor::Number<N>;

  GetDiffDiffMatImpl(E &e, T1 &t_a, FTensor::Tensor2_symmetric<VT2, DimT2> &t_S)
      : r(e), e(e), tA(t_a), tS(t_S) {}
  SecondMatrixDirectiveImpl<E, C> r;
  E &e;
  T1 &tA;
  FTensor::Tensor2_symmetric<VT2, DimT2> &tS;

  template <int I, int J, int K, int L, int M, int N>
  inline auto add(const Number<I> &, const Number<J> &, const Number<K> &,
                  const Number<L> &, const Number<M> &, const Number<N> &) {

    if constexpr (N != M)
      return (2 * tS(Number<M - 1>(), Number<N - 1>())) *
                 r.eval(NumberDim(), Number<M - 1>(), Number<N - 1>(),
                        Number<I - 1>(), Number<J - 1>(), Number<K - 1>(),
                        Number<L - 1>())

             +

             add(Number<I>(), Number<J>(), Number<K>(), Number<L>(),
                 Number<M>(), Number<N - 1>());
    else
      return tS(Number<M - 1>(), Number<N - 1>()) *
                 r.eval(NumberDim(), Number<M - 1>(), Number<N - 1>(),
                        Number<I - 1>(), Number<J - 1>(), Number<K - 1>(),
                        Number<L - 1>())

             +

             add(Number<I>(), Number<J>(), Number<K>(), Number<L>(),
                 Number<M>(), Number<N - 1>());
  }

  template <int I, int J, int K, int L, int M>
  inline auto add(const Number<I> &, const Number<J> &, const Number<K> &,
                  const Number<L> &, const Number<M> &, const Number<1> &) {
    return (2 * tS(Number<M - 1>(), Number<0>())) *
               r.eval(NumberDim(), Number<M - 1>(), Number<0>(),
                      Number<I - 1>(), Number<J - 1>(), Number<K - 1>(),
                      Number<L - 1>())

           +

           add(Number<I>(), Number<J>(), Number<K>(), Number<L>(),
               Number<M - 1>(), Number<M - 1>());
  }

  template <int I, int J, int K, int L>
  inline auto add(const Number<I> &, const Number<J> &, const Number<K> &,
                  const Number<L> &, const Number<1> &, const Number<1> &) {
    return tS(Number<0>(), Number<0>()) *
           r.eval(NumberDim(), Number<0>(), Number<0>(), Number<I - 1>(),
                  Number<J - 1>(), Number<K - 1>(), Number<L - 1>());
  }

  template <int I, int J, int K, int L>
  inline void set(const Number<I> &, const Number<J> &, const Number<K> &,
                  const Number<L> &) {

    set(Number<I>(), Number<J>(), Number<K>(), Number<L - 1>());

    tA(Number<I - 1>(), Number<J - 1>(), Number<K - 1>(), Number<L - 1>()) =
        add(Number<I>(), Number<J>(), Number<K>(), Number<L>(), NumberDim(),
            NumberDim());

    // Major symmetry
    if constexpr (K != I || L != J)
      tA(Number<K - 1>(), Number<L - 1>(), Number<I - 1>(), Number<J - 1>()) =
          tA(Number<I - 1>(), Number<J - 1>(), Number<K - 1>(),
             Number<L - 1>());
  }

  template <int I, int J, int K>
  inline void set(const Number<I> &, const Number<J> &, const Number<K> &,
                  const Number<0> &) {
    set(Number<I>(), Number<J>(), Number<K - 1>(), Number<K - 1>());
  }

  template <int I, int J>
  inline void set(const Number<I> &, const Number<J> &, const Number<0> &,
                  const Number<0> &) {
    set(Number<I>(), Number<J - 1>(), Number<I>(), Number<J - 1>());
  }

  template <int I, int K>
  inline void set(const Number<I> &, const Number<0> &, const Number<K> &,
                  const Number<0> &) {
    set(Number<I - 1>(), Number<I - 1>(), Number<I - 1>(), Number<I - 1>());
  }

  inline void set(const Number<0> &, const Number<0> &, const Number<0> &,
                  const Number<0> &) {}
};

template <typename T1, typename T2, int NB, int Dim> struct EigenMatrixImp {

  using Val = const FTensor::Tensor1<T1, Dim>;
  using Vec = const FTensor::Tensor2<T2, Dim, Dim>;
  using Fun = boost::function<double(const double)>;
  using V = double; // typename FTensor::promote<T1, T2>::V;

  template <int N> using Number = FTensor::Number<N>;
  template <char c> using I = typename FTensor::Index<c, Dim>;

  using NumberNb = Number<NB>;
  using NumberDim = Number<Dim>;

  EigenMatrixImp(Val &t_val, Vec &t_vec) : tVal(t_val), tVec(t_vec) {

    FTensor::Index<'i', Dim> i;
    FTensor::Index<'j', Dim> j;
    FTensor::Index<'k', Dim> k;
    FTensor::Index<'l', Dim> l;

    for (auto aa = 0; aa != Dim; ++aa) {
      auto &M = aM[aa];
      for (auto ii = 0; ii != Dim; ++ii)
        for (auto jj = 0; jj <= ii; ++jj)
          M(ii, jj) = tVec(aa, ii) * tVec(aa, jj);
    }

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
          auto &S = aS[get_sym_index(aa, bb, Dim)];
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
        .set(NumberDim(), NumberDim());
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
   * @return auto derivatives, forth order tensor with minor simetries
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
        aF2(aa, bb) = aF(aa, bb) * aF(aa, bb);
      }

    using T3 = FTensor::Ddg<V, Dim, Dim>;
    T3 t_diff_A;
    GetDiffMatImpl<EigenMatrixImp<T1, T2, NB, Dim>, V, T3>(*this, t_diff_A)
        .set(NumberDim(), NumberDim(), NumberDim(), NumberDim());
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
   * @return auto second derivatives, forth order tensor with minor simetries
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
        aF2(aa, bb) = aF(aa, bb) * aF(aa, bb);
      }

    FTensor::Index<'i', Dim> i;
    FTensor::Index<'j', Dim> j;
    FTensor::Index<'k', Dim> k;
    FTensor::Index<'l', Dim> l;

    for (auto aa = 0; aa != Dim; ++aa) {
      for (auto bb = 0; bb != Dim; ++bb) {
        if (aa != bb) {
          const auto &S = aS[get_sym_index(aa, bb, Dim)];
          const auto &M = aM[aa];
          auto &SMmn = aSM[get_nodiag_index(aa, bb, Dim)];
          for (auto mm = 0; mm != Dim; ++mm) {
            for (auto nn = mm; nn != Dim; ++nn) {
              SMmn[get_sym_index(mm, nn, Dim)](i, j, k, l) =
                  S(i, j, k, l) * M(mm, nn);
            }
          }
        }
      }
    }

    for (auto aa = 0; aa != Dim; ++aa) {
      for (auto mm = 0; mm != Dim; ++mm) {
        for (auto nn = mm; nn != Dim; ++nn) {
          d2MType0[aa][get_sym_index(mm, nn, Dim)](i, j, k, l) = 0;
        }
      }
    }

    if constexpr (NB == 3)
      for (auto aa = 0; aa != Dim; ++aa) {
        for (auto bb = 0; bb != Dim; ++bb) {
          if (aa != bb) {
            const V v = dfVal(aa) * aF(aa, bb);
            for (auto mm = 0; mm != Dim; ++mm) {
              for (auto nn = mm; nn != Dim; ++nn) {
                d2MType0[aa][get_sym_index(mm, nn, Dim)](i, j, k, l) +=
                    v * aSM[get_nodiag_index(aa, bb, Dim)]
                           [get_sym_index(mm, nn, Dim)](i, j, k, l);
              }
            }
          }
        }
      }

    if constexpr (NB == 2)
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
                d2MType0[aa][get_sym_index(mm, nn, Dim)](i, j, k, l) +=
                    v * aSM[get_nodiag_index(aa, bb, Dim)]
                           [get_sym_index(mm, nn, Dim)](i, j, k, l);
              }
            }
          }
        }
      }

    if constexpr (NB == 1)
      for (auto aa = 0; aa != Dim; ++aa) {
        for (auto bb = 0; bb != Dim; ++bb) {
          if (aa != bb) {
            const V v = ddfVal(aa) / 2;
            for (auto mm = 0; mm != Dim; ++mm) {
              for (auto nn = mm; nn != Dim; ++nn) {
                d2MType0[aa][get_sym_index(mm, nn, Dim)](i, j, k, l) +=
                    v * aSM[get_nodiag_index(aa, bb, Dim)]
                           [get_sym_index(mm, nn, Dim)](i, j, k, l);
              }
            }
          }
        }
      }

    for (auto aa = 0; aa != Dim; ++aa) {
      for (auto mm = 0; mm != Dim; ++mm) {
        for (auto nn = 0; nn != Dim; ++nn) {
          if (nn != mm)
            d2MType1[mm][get_sym_index(aa, nn, Dim)](i, j, k, l) = 0;
        }
      }
    }

    if constexpr (NB == 3)
      for (auto aa = 0; aa != Dim; ++aa) {
        for (auto bb = 0; bb != Dim; ++bb) {
          if (aa != bb) {
            const auto &S = aS[get_sym_index(aa, bb, Dim)];
            const auto v0 = aF(aa, bb);
            for (auto cc = 0; cc != Dim; ++cc) {
              for (auto dd = 0; dd != Dim; ++dd) {
                if (cc != dd) {
                  const double v1 = fVal(cc) * aF(cc, dd);
                  d2MType1[dd][get_sym_index(aa, cc, Dim)](i, j, k, l) +=
                      (v1 * v0) * S(i, j, k, l);
                }
              }
            }
          }
        }
      }

    if constexpr (NB == 2)
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

                  if (r)
                    d2MType1[dd][get_sym_index(aa, cc, Dim)](i, j, k, l) +=
                        r * aS[get_sym_index(aa, bb, Dim)](i, j, k, l);
                }
            }
        }
      }

    if constexpr (NB == 1)
      for (auto aa = 0; aa != Dim; ++aa) {
        for (auto bb = 0; bb != Dim; ++bb) {
          if (aa != bb) {
            for (auto cc = 0; cc != Dim; ++cc) {
              for (auto dd = 0; dd != Dim; ++dd) {
                if (cc != dd) {
                  if ((bb != dd) && (aa != dd && bb != cc)) {
                    const double r = ddfVal(cc) / 4;
                    d2MType1[dd][get_sym_index(aa, cc, Dim)](i, j, k, l) +=
                        r * aS[get_sym_index(aa, bb, Dim)](i, j, k, l);
                  }
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
  FTensor::Ddg<V, Dim, Dim> aS[(Dim * (Dim + 1)) / 2];
  FTensor::Ddg<V, Dim, Dim> aSM[(Dim - 1) * Dim][(Dim * (Dim + 1)) / 2];
  FTensor::Ddg<V, Dim, Dim> d2MType0[Dim][(Dim * (Dim + 1)) / 2];
  FTensor::Ddg<V, Dim, Dim> d2MType1[Dim][(Dim * (Dim + 1)) / 2];
  FTensor::Tensor2_symmetric<V, Dim> aF2;
  FTensor::Tensor2<V, Dim, Dim> aF;
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