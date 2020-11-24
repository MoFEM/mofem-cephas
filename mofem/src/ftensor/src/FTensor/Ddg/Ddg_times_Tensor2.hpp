/* This file has all of the declarations for expressions like
   Ddg*Tensor2 and Tensor2*Ddg, yielding a
   Tensor2_symmetric. */

#pragma once

namespace FTensor
{
  /* A(i,j,k,l)*B(k,l) */

  template <class A, class B, class T, class U, int Dim01, int Dim23, char i,
            char j, char k, char l>
  class Ddg_times_Tensor2_23
  {
    Ddg_Expr<A, T, Dim01, Dim23, i, j, k, l> iterA;
    Tensor2_Expr<B, U, Dim23, Dim23, k, l> iterB;

    template <int Current_Dim0, int Current_Dim1>
    typename promote<T, U>::V
    eval(const int N1, const int N2, const Number<Current_Dim0> &,
         const Number<Current_Dim1> &) const
    {
      return iterA(N1, N2, Current_Dim0 - 1, Current_Dim1 - 1)
               * iterB(Current_Dim0 - 1, Current_Dim1 - 1)
             + eval(N1, N2, Number<Current_Dim0 - 1>(),
                    Number<Current_Dim1>());
    }
    template <int Current_Dim1>
    typename promote<T, U>::V
    eval(const int N1, const int N2, const Number<1> &,
         const Number<Current_Dim1> &) const
    {
      return iterA(N1, N2, 0, Current_Dim1 - 1) * iterB(0, Current_Dim1 - 1)
             + eval(N1, N2, Number<Dim23>(), Number<Current_Dim1 - 1>());
    }
    typename promote<T, U>::V eval(const int N1, const int N2,
                                   const Number<1> &, const Number<1> &) const
    {
      return iterA(N1, N2, 0, 0) * iterB(0, 0);
    }

  public:
    Ddg_times_Tensor2_23(const Ddg_Expr<A, T, Dim01, Dim23, i, j, k, l> &a,
                         const Tensor2_Expr<B, U, Dim23, Dim23, k, l> &b)
        : iterA(a), iterB(b)
    {}
    typename promote<T, U>::V operator()(const int N1, const int N2) const
    {
      return eval(N1, N2, Number<Dim23>(), Number<Dim23>());
    }
  };

  template <class A, class B, class T, class U, int Dim01, int Dim23, char i,
            char j, char k, char l>
  Tensor2_symmetric_Expr<
    Ddg_times_Tensor2_23<A, B, T, U, Dim01, Dim23, i, j, k, l>,
    typename promote<T, U>::V, Dim01, i, j>
  operator*(const Ddg_Expr<A, T, Dim01, Dim23, i, j, k, l> &a,
            const Tensor2_Expr<B, U, Dim23, Dim23, k, l> &b)
  {
    using TensorExpr
      = Ddg_times_Tensor2_23<A, B, T, U, Dim01, Dim23, i, j, k, l>;
    return Tensor2_symmetric_Expr<TensorExpr, typename promote<T, U>::V, Dim01,
                                  i, j>(TensorExpr(a, b));
  }

  /* B(k,l)*A(i,j,k,l) */

  template <class A, class B, class T, class U, int Dim01, int Dim23, char i,
            char j, char k, char l>
  Tensor2_symmetric_Expr<
    Ddg_times_Tensor2_23<A, B, T, U, Dim01, Dim23, i, j, k, l>,
    typename promote<T, U>::V, Dim01, i, j>
  operator*(const Tensor2_Expr<B, U, Dim23, Dim23, k, l> &b,
            const Ddg_Expr<A, T, Dim01, Dim23, i, j, k, l> &a)
  {
    using TensorExpr
      = Ddg_times_Tensor2_23<A, B, T, U, Dim01, Dim23, i, j, k, l>;
    return Tensor2_symmetric_Expr<TensorExpr, typename promote<T, U>::V, Dim01,
                                  i, j>(TensorExpr(a, b));
  }

  /* A(i,j,k,l)*B(l,k) */

  template <class A, class B, class T, class U, int Dim01, int Dim23, char i,
            char j, char k, char l>
  class Ddg_times_Tensor2_32
  {
    Ddg_Expr<A, T, Dim01, Dim23, i, j, k, l> iterA;
    Tensor2_Expr<B, U, Dim23, Dim23, l, k> iterB;

    template <int Current_Dim0, int Current_Dim1>
    typename promote<T, U>::V
    eval(const int N1, const int N2, const Number<Current_Dim0> &,
         const Number<Current_Dim1> &) const
    {
      return iterA(N1, N2, Current_Dim0 - 1, Current_Dim1 - 1)
               * iterB(Current_Dim0 - 1, Current_Dim1 - 1)
             + eval(N1, N2, Number<Current_Dim0 - 1>(),
                    Number<Current_Dim1>());
    }
    template <int Current_Dim1>
    typename promote<T, U>::V
    eval(const int N1, const int N2, const Number<1> &,
         const Number<Current_Dim1> &) const
    {
      return iterA(N1, N2, 0, Current_Dim1 - 1) * iterB(0, Current_Dim1 - 1)
             + eval(N1, N2, Number<Dim23>(), Number<Current_Dim1 - 1>());
    }
    typename promote<T, U>::V eval(const int N1, const int N2,
                                   const Number<1> &, const Number<1> &) const
    {
      return iterA(N1, N2, 0, 0) * iterB(0, 0);
    }

  public:
    Ddg_times_Tensor2_32(const Ddg_Expr<A, T, Dim01, Dim23, i, j, k, l> &a,
                         const Tensor2_Expr<B, U, Dim23, Dim23, l, k> &b)
        : iterA(a), iterB(b)
    {}
    typename promote<T, U>::V operator()(const int N1, const int N2) const
    {
      return eval(N1, N2, Number<Dim23>(), Number<Dim23>());
    }
  };

  template <class A, class B, class T, class U, int Dim01, int Dim23, char i,
            char j, char k, char l>
  Tensor2_symmetric_Expr<
    Ddg_times_Tensor2_32<A, B, T, U, Dim01, Dim23, i, j, k, l>,
    typename promote<T, U>::V, Dim01, i, j>
  operator*(const Ddg_Expr<A, T, Dim01, Dim23, i, j, k, l> &a,
            const Tensor2_Expr<B, U, Dim23, Dim23, l, k> &b)
  {
    using TensorExpr
      = Ddg_times_Tensor2_32<A, B, T, U, Dim01, Dim23, i, j, k, l>;
    return Tensor2_symmetric_Expr<TensorExpr, typename promote<T, U>::V, Dim01,
                                  i, j>(TensorExpr(a, b));
  }

  /* B(l,k)*A(i,j,k,l) */

  template <class A, class B, class T, class U, int Dim01, int Dim23, char i,
            char j, char k, char l>
  Tensor2_symmetric_Expr<
    Ddg_times_Tensor2_32<A, B, T, U, Dim01, Dim23, i, j, k, l>,
    typename promote<T, U>::V, Dim01, i, j>
  operator*(const Tensor2_Expr<B, U, Dim23, Dim23, l, k> &b,
            const Ddg_Expr<A, T, Dim01, Dim23, i, j, k, l> &a)
  {
    using TensorExpr
      = Ddg_times_Tensor2_32<A, B, T, U, Dim01, Dim23, i, j, k, l>;
    return Tensor2_symmetric_Expr<TensorExpr, typename promote<T, U>::V, Dim01,
                                  i, j>(TensorExpr(a, b));
  }

  /* A(i,j,k,l)*B(i,j) */

  template <class A, class B, class T, class U, int Dim01, int Dim23, char i,
            char j, char k, char l>
  class Ddg_times_Tensor2_01 {
    Ddg_Expr<A, T, Dim01, Dim23, i, j, k, l> iterA;
    Tensor2_Expr<B, U, Dim01, Dim01, i, j> iterB;

    template <int Current_Dim0, int Current_Dim1>
    typename promote<T, U>::V eval(const int N1, const int N2,
                                   const Number<Current_Dim0> &,
                                   const Number<Current_Dim1> &) const {
      return iterA(Current_Dim0 - 1, Current_Dim1 - 1, N1, N2) *
                 iterB(Current_Dim0 - 1, Current_Dim1 - 1) +
             eval(N1, N2, Number<Current_Dim0 - 1>(), Number<Current_Dim1>());
    }
    template <int Current_Dim1>
    typename promote<T, U>::V eval(const int N1, const int N2,
                                   const Number<1> &,
                                   const Number<Current_Dim1> &) const {
      return iterA(0, Current_Dim1 - 1, N1, N2) * iterB(0, Current_Dim1 - 1) +
             eval(N1, N2, Number<Dim01>(), Number<Current_Dim1 - 1>());
    }
    typename promote<T, U>::V eval(const int N1, const int N2,
                                   const Number<1> &, const Number<1> &) const {
      return iterA(0, 0, N1, N2) * iterB(0, 0);
    }

  public:
    Ddg_times_Tensor2_01(const Ddg_Expr<A, T, Dim01, Dim23, i, j, k, l> &a,
                         const Tensor2_Expr<B, U, Dim01, Dim01, i, j> &b)
        : iterA(a), iterB(b) {}
    typename promote<T, U>::V operator()(const int N1, const int N2) const {
      return eval(N1, N2, Number<Dim01>(), Number<Dim01>());
    }
  };

  template <class A, class B, class T, class U, int Dim01, int Dim23, char i,
            char j, char k, char l>
  Tensor2_symmetric_Expr<
      Ddg_times_Tensor2_01<A, B, T, U, Dim01, Dim23, i, j, k, l>,
      typename promote<T, U>::V, Dim23, k, l>
  operator*(const Ddg_Expr<A, T, Dim01, Dim23, i, j, k, l> &a,
            const Tensor2_Expr<B, U, Dim01, Dim01, i, j> &b) {
    using TensorExpr =
        Ddg_times_Tensor2_01<A, B, T, U, Dim01, Dim23, i, j, k, l>;
    return Tensor2_symmetric_Expr<TensorExpr, typename promote<T, U>::V, Dim23,
                                  k, l>(TensorExpr(a, b));
  }

  /* B(i,j)*A(i,j,k,l) */

  template <class A, class B, class T, class U, int Dim01, int Dim23, char i,
            char j, char k, char l>
  Tensor2_symmetric_Expr<
      Ddg_times_Tensor2_01<A, B, T, U, Dim01, Dim23, i, j, k, l>,
      typename promote<T, U>::V, Dim23, k, l>
  operator*(const Tensor2_Expr<B, U, Dim01, Dim01, i, j> &b,
            const Ddg_Expr<A, T, Dim01, Dim23, i, j, k, l> &a) {
    using TensorExpr =
        Ddg_times_Tensor2_01<A, B, T, U, Dim01, Dim23, i, j, k, l>;
    return Tensor2_symmetric_Expr<TensorExpr, typename promote<T, U>::V, Dim23,
                                  k, l>(TensorExpr(a, b));
  }

  /* A(j,i,k,l)*B(i,j) */

  template <class A, class B, class T, class U, int Dim01, int Dim23, char i,
            char j, char k, char l>
  Tensor2_symmetric_Expr<
      Ddg_times_Tensor2_01<A, B, T, U, Dim01, Dim23, i, j, k, l>,
      typename promote<T, U>::V, Dim23, k, l>
  operator*(const Ddg_Expr<A, T, Dim01, Dim23, j, i, k, l> &a,
            const Tensor2_Expr<B, U, Dim01, Dim01, i, j> &b) {
    using TensorExpr =
        Ddg_times_Tensor2_01<A, B, T, U, Dim01, Dim23, i, j, k, l>;
    return Tensor2_symmetric_Expr<TensorExpr, typename promote<T, U>::V, Dim23,
                                  k, l>(TensorExpr(a, b));
  }

  /* B(i,j)*A(j,i,k,l) */

  template <class A, class B, class T, class U, int Dim01, int Dim23, char i,
            char j, char k, char l>
  Tensor2_symmetric_Expr<
      Ddg_times_Tensor2_01<A, B, T, U, Dim01, Dim23, i, j, k, l>,
      typename promote<T, U>::V, Dim23, k, l>
  operator*(const Tensor2_Expr<B, U, Dim01, Dim01, i, j> &b,
            const Ddg_Expr<A, T, Dim01, Dim23, j, i, k, l> &a) {
    using TensorExpr =
        Ddg_times_Tensor2_01<A, B, T, U, Dim01, Dim23, i, j, k, l>;
    return Tensor2_symmetric_Expr<TensorExpr, typename promote<T, U>::V, Dim23,
                                  k, l>(TensorExpr(a, b));
  }

  /* This operatores will yield tensor 2 */

  /* A(i,j,k,l)*B(j,l) */

  template <class A, class B, class T, class U, int Dim01, int Dim23, char i,
            char j, char k, char l>
  class Ddg_times_Tensor2_13
  {
    Ddg_Expr<A, T, Dim01, Dim23, i, j, k, l> iterA;
    Tensor2_Expr<B, U, Dim01, Dim23, j, l> iterB;

    template <int Current_Dim0, int Current_Dim1>
    typename promote<T, U>::V
    eval(const int N1, const int N2, const Number<Current_Dim0> &,
         const Number<Current_Dim1> &) const
    {
      return iterA(N1, Current_Dim0 - 1, N2, Current_Dim1 - 1) *
                 iterB(Current_Dim0 - 1, Current_Dim1 - 1) +
             eval(N1, N2, Number<Current_Dim0 - 1>(), Number<Current_Dim1>());
    }
    template <int Current_Dim1>
    typename promote<T, U>::V
    eval(const int N1, const int N2, const Number<1> &,
         const Number<Current_Dim1> &) const
    {
      return iterA(N1, 0, N2, Current_Dim1 - 1) * iterB(0, Current_Dim1 - 1) +
             eval(N1, N2, Number<Dim01>(), Number<Current_Dim1 - 1>());
    }
    typename promote<T, U>::V eval(const int N1, const int N2,
                                   const Number<1> &, const Number<1> &) const
    {
      return iterA(N1, 0, N2, 0) * iterB(0, 0);
    }

  public:
    Ddg_times_Tensor2_13(const Ddg_Expr<A, T, Dim01, Dim23, i, j, k, l> &a,
                         const Tensor2_Expr<B, U, Dim01, Dim23, j, l> &b)
        : iterA(a), iterB(b)
    {}
    typename promote<T, U>::V operator()(const int N1, const int N2) const
    {
      return eval(N1, N2, Number<Dim01>(), Number<Dim23>());
    }
  };

  template <class A, class B, class T, class U, int Dim01, int Dim23, char i,
            char j, char k, char l>
  Tensor2_Expr<
    Ddg_times_Tensor2_13<A, B, T, U, Dim01, Dim23, i, j, k, l>,
    typename promote<T, U>::V, Dim01, Dim23, i, k>
  operator*(const Ddg_Expr<A, T, Dim01, Dim23, i, j, k, l> &a,
            const Tensor2_Expr<B, U, Dim01, Dim23, j, l> &b)
  {
    using TensorExpr
      = Ddg_times_Tensor2_13<A, B, T, U, Dim01, Dim23, i, j, k, l>;
    return Tensor2_Expr<TensorExpr, typename promote<T, U>::V, Dim01, Dim23, i,
                        k>(TensorExpr(a, b));
  }

  /* B(j,l)*A(i,j,k,l) */

  template <class A, class B, class T, class U, int Dim01, int Dim23, char i,
            char j, char k, char l>
  Tensor2_Expr<Ddg_times_Tensor2_13<A, B, T, U, Dim01, Dim23, i, j, k, l>,
               typename promote<T, U>::V, Dim01, Dim23, i, k>
  operator*(const Tensor2_Expr<B, U, Dim01, Dim23, j, l> &b,
            const Ddg_Expr<A, T, Dim01, Dim23, i, j, k, l> &a) {
    using TensorExpr
      = Ddg_times_Tensor2_13<A, B, T, U, Dim01, Dim23, i, j, k, l>;
    return Tensor2_Expr<TensorExpr, typename promote<T, U>::V, Dim01, Dim23, i,
                        k>(TensorExpr(a, b));
  }

  /* This file has all of the declarations for expressions like
     Ddg*Tensor2 and Tensor2*Ddg, yielding a
     Tensor4. */

  // FIXME: That should create T4 with first two indices symmetric two other
  // not. At this point will not create ideal code.
  // TODO: Not all possible permutations are included.
  // TODO: Check dimensions could be errors

  /* A(i,j,k,l)*B(l,m) */

  template <class A, class B, class T, class U, int Dim01, int Dim23, int Dim4,
            char i, char j, char k, char l, char m>
  class Ddg_times_Tensor2_3_1
  {
    const Ddg_Expr<A, T, Dim01, Dim23, i, j, k, l> iterA;
    const Tensor2_Expr<B, U, Dim23, Dim4, l, m> iterB;

    template <int Current_Dim0>
    typename promote<T, U>::V
    eval(const int N1, const int N2, const int N3, const int N4,
         const Number<Current_Dim0> &) const
    {
      return iterA(N1, N2, N3, Current_Dim0 - 1) * iterB(Current_Dim0 - 1, N4)
             + eval(N1, N2, N3, N4, Number<Current_Dim0 - 1>());
    }
    typename promote<T, U>::V eval(const int N1, const int N2, const int N3,
                                   const int N4, const Number<1> &) const
    {
      return iterA(N1, N2, N3, 0) * iterB(0, N4);
    }

  public:
    Ddg_times_Tensor2_3_1(const Ddg_Expr<A, T, Dim01, Dim23, i, j, k, l> &a,
                          const Tensor2_Expr<B, U, Dim23, Dim4, l, m> &b)
        : iterA(a), iterB(b)
    {}
    typename promote<T, U>::V
    operator()(const int N1, const int N2, const int N3, const int N4) const
    {
      return eval(N1, N2, N3, N4, Number<Dim23>());
    }
  };

  template <class A, class B, class T, class U, int Dim01, int Dim23, int Dim4,
            char i, char j, char k, char l, char m>
  inline const Tensor4_Expr<
    const Ddg_times_Tensor2_3_1<A, B, T, U, Dim01, Dim23, Dim4, i, j, k, l, m>,
    typename promote<T, U>::V, Dim01, Dim01, Dim23, Dim4, i, j, k, m>
  operator*(const Ddg_Expr<A, T, Dim01, Dim23, i, j, k, l> &a,
            const Tensor2_Expr<B, U, Dim23, Dim4, l, m> &b)
  {
    typedef const Ddg_times_Tensor2_3_1<A, B, T, U, Dim01, Dim23, Dim4, i, j,
                                        k, l, m>
      TensorExpr;
    return Tensor4_Expr<TensorExpr, typename promote<T, U>::V, Dim01, Dim01,
                        Dim23, Dim4, i, j, k, m>(TensorExpr(a, b));
  }

  /* B(l,m)*A(i,j,k,l) */

  template <class A, class B, class T, class U, int Dim01, int Dim23, int Dim4,
            char i, char j, char k, char l, char m>
  inline const Tensor4_Expr<
    const Ddg_times_Tensor2_3_1<A, B, T, U, Dim01, Dim23, Dim4, i, j, k, l, m>,
    typename promote<T, U>::V, Dim01, Dim01, Dim23, Dim4, i, j, k, m>
  operator*(const Tensor2_Expr<B, U, Dim23, Dim4, l, m> &b,
            const Ddg_Expr<A, T, Dim01, Dim23, i, j, k, l> &a)
  {
    typedef const Ddg_times_Tensor2_3_1<A, B, T, U, Dim01, Dim23, Dim4, i, j,
                                        k, l, m>
      TensorExpr;
    return Tensor4_Expr<TensorExpr, typename promote<T, U>::V, Dim01, Dim01,
                        Dim23, Dim4, i, j, k, m>(TensorExpr(a, b));
  }

  /* A(i,j,k,l)*B(m,l) */

  template <class A, class B, class T, class U, int Dim01, int Dim23, int Dim4,
            char i, char j, char k, char l, char m>
  class Ddg_times_Tensor2_3_0
  {
    const Ddg_Expr<A, T, Dim01, Dim23, i, j, k, l> iterA;
    const Tensor2_Expr<B, U, Dim4, Dim23, m, l> iterB;

    template <int Current_Dim0>
    typename promote<T, U>::V
    eval(const int N1, const int N2, const int N3, const int N4,
         const Number<Current_Dim0> &) const
    {
      return iterA(N1, N2, N3, Current_Dim0 - 1) * iterB(N4, Current_Dim0 - 1)
             + eval(N1, N2, N3, N4, Number<Current_Dim0 - 1>());
    }
    typename promote<T, U>::V eval(const int N1, const int N2, const int N3,
                                   const int N4, const Number<1> &) const
    {
      return iterA(N1, N2, N3, 0) * iterB(N4, 0);
    }

  public:
    Ddg_times_Tensor2_3_0(const Ddg_Expr<A, T, Dim01, Dim23, i, j, k, l> &a,
                          const Tensor2_Expr<B, U, Dim4, Dim23, m, l> &b)
        : iterA(a), iterB(b)
    {}
    typename promote<T, U>::V
    operator()(const int N1, const int N2, const int N3, const int N4) const
    {
      return eval(N1, N2, N3, N4, Number<Dim23>());
    }
  };

  template <class A, class B, class T, class U, int Dim01, int Dim23, int Dim4,
            char i, char j, char k, char l, char m>
  inline const Tensor4_Expr<
    const Ddg_times_Tensor2_3_0<A, B, T, U, Dim01, Dim23, Dim4, i, j, k, l, m>,
    typename promote<T, U>::V, Dim01, Dim01, Dim23, Dim4, i, j, k, m>
  operator*(const Ddg_Expr<A, T, Dim01, Dim23, i, j, k, l> &a,
            const Tensor2_Expr<B, U, Dim4, Dim23, m, l> &b)
  {
    typedef const Ddg_times_Tensor2_3_0<A, B, T, U, Dim01, Dim23, Dim4, i, j,
                                        k, l, m>
      TensorExpr;
    return Tensor4_Expr<TensorExpr, typename promote<T, U>::V, Dim01, Dim01,
                        Dim23, Dim4, i, j, k, m>(TensorExpr(a, b));
  }

  /* B(m,l)*A(i,j,k,l) */

  template <class A, class B, class T, class U, int Dim01, int Dim23, int Dim4,
            char i, char j, char k, char l, char m>
  inline const Tensor4_Expr<
    const Ddg_times_Tensor2_3_0<A, B, T, U, Dim01, Dim23, Dim4, i, j, k, l, m>,
    typename promote<T, U>::V, Dim01, Dim01, Dim23, Dim4, i, j, k, m>
  operator*(const Tensor2_Expr<B, U, Dim4, Dim23, m, l> &b,
            const Ddg_Expr<A, T, Dim01, Dim23, i, j, k, l> &a)
  {
    typedef const Ddg_times_Tensor2_3_0<A, B, T, U, Dim01, Dim23, Dim4, i, j,
                                        k, l, m>
      TensorExpr;
    return Tensor4_Expr<TensorExpr, typename promote<T, U>::V, Dim01, Dim01,
                        Dim23, Dim4, i, j, k, m>(TensorExpr(a, b));
  }

  /* A(i,j,k,l)*B(m,k) */

  template <class A, class B, class T, class U, int Dim01, int Dim23, int Dim4,
            char i, char j, char k, char l, char m>
  class Ddg_times_Tensor2_2_0
  {
    const Ddg_Expr<A, T, Dim01, Dim23, i, j, k, l> iterA;
    const Tensor2_Expr<B, U, Dim4, Dim23, m, k> iterB;

    template <int Current_Dim0>
    typename promote<T, U>::V
    eval(const int N1, const int N2, const int N3, const int N4,
         const Number<Current_Dim0> &) const
    {
      return iterA(N1, N2, Current_Dim0 - 1, N4) * iterB(N3, Current_Dim0 - 1)
             + eval(N1, N2, N3, N4, Number<Current_Dim0 - 1>());
    }
    typename promote<T, U>::V eval(const int N1, const int N2, const int N3,
                                   const int N4, const Number<1> &) const
    {
      return iterA(N1, N2, 0, N4) * iterB(N3, 0);
    }

  public:
    Ddg_times_Tensor2_2_0(const Ddg_Expr<A, T, Dim01, Dim23, i, j, k, l> &a,
                          const Tensor2_Expr<B, U, Dim4, Dim23, m, k> &b)
        : iterA(a), iterB(b)
    {}
    typename promote<T, U>::V
    operator()(const int N1, const int N2, const int N3, const int N4) const
    {
      return eval(N1, N2, N3, N4, Number<Dim23>());
    }
  };

  template <class A, class B, class T, class U, int Dim01, int Dim23, int Dim4,
            char i, char j, char k, char l, char m>
  inline const Tensor4_Expr<
    const Ddg_times_Tensor2_2_0<A, B, T, U, Dim01, Dim23, Dim4, i, j, k, l, m>,
    typename promote<T, U>::V, Dim01, Dim01, Dim23, Dim4, i, j, m, l>
  operator*(const Ddg_Expr<A, T, Dim01, Dim23, i, j, k, l> &a,
            const Tensor2_Expr<B, U, Dim4, Dim23, m, k> &b)
  {
    typedef const Ddg_times_Tensor2_2_0<A, B, T, U, Dim01, Dim23, Dim4, i, j,
                                        k, l, m>
      TensorExpr;
    return Tensor4_Expr<TensorExpr, typename promote<T, U>::V, Dim01, Dim01,
                        Dim23, Dim4, i, j, m, l>(TensorExpr(a, b));
  }

  /* B(m,k)*A(i,j,k,l) */

  template <class A, class B, class T, class U, int Dim01, int Dim23, int Dim4,
            char i, char j, char k, char l, char m>
  inline const Tensor4_Expr<
    const Ddg_times_Tensor2_2_0<A, B, T, U, Dim01, Dim23, Dim4, i, j, k, l, m>,
    typename promote<T, U>::V, Dim01, Dim01, Dim23, Dim4, i, j, m, l>
  operator*(const Tensor2_Expr<B, U, Dim4, Dim23, m, k> &b,
            const Ddg_Expr<A, T, Dim01, Dim23, i, j, k, l> &a)
  {
    typedef const Ddg_times_Tensor2_2_0<A, B, T, U, Dim01, Dim23, Dim4, i, j,
                                        k, l, m>
      TensorExpr;
    return Tensor4_Expr<TensorExpr, typename promote<T, U>::V, Dim01, Dim01,
                        Dim23, Dim4, i, j, m, l>(TensorExpr(a, b));
  }

  /* A(i,j,k,l)*B(k,m) */

  template <class A, class B, class T, class U, int Dim01, int Dim23, int Dim4,
            char i, char j, char k, char l, char m>
  class Ddg_times_Tensor2_2_1
  {
    const Ddg_Expr<A, T, Dim01, Dim23, i, j, k, l> iterA;
    const Tensor2_Expr<B, U, Dim4, Dim23, k, m> iterB;

    template <int Current_Dim0>
    typename promote<T, U>::V
    eval(const int N1, const int N2, const int N3, const int N4,
         const Number<Current_Dim0> &) const
    {
      return iterA(N1, N2, Current_Dim0 - 1, N4) * iterB(Current_Dim0 - 1, N3)
             + eval(N1, N2, N3, N4, Number<Current_Dim0 - 1>());
    }
    typename promote<T, U>::V eval(const int N1, const int N2, const int N3,
                                   const int N4, const Number<1> &) const
    {
      return iterA(N1, N2, 0, N4) * iterB(0, N3);
    }

  public:
    Ddg_times_Tensor2_2_1(const Ddg_Expr<A, T, Dim01, Dim23, i, j, k, l> &a,
                          const Tensor2_Expr<B, U, Dim4, Dim23, k, m> &b)
        : iterA(a), iterB(b)
    {}
    typename promote<T, U>::V
    operator()(const int N1, const int N2, const int N3, const int N4) const
    {
      return eval(N1, N2, N3, N4, Number<Dim23>());
    }
  };

  template <class A, class B, class T, class U, int Dim01, int Dim23, int Dim4,
            char i, char j, char k, char l, char m>
  inline const Tensor4_Expr<
    const Ddg_times_Tensor2_2_1<A, B, T, U, Dim01, Dim23, Dim4, i, j, k, l, m>,
    typename promote<T, U>::V, Dim01, Dim01, Dim23, Dim4, i, j, m, l>
  operator*(const Ddg_Expr<A, T, Dim01, Dim23, i, j, k, l> &a,
            const Tensor2_Expr<B, U, Dim4, Dim23, k, m> &b)
  {
    typedef const Ddg_times_Tensor2_2_1<A, B, T, U, Dim01, Dim23, Dim4, i, j,
                                        k, l, m>
      TensorExpr;
    return Tensor4_Expr<TensorExpr, typename promote<T, U>::V, Dim01, Dim01,
                        Dim23, Dim4, i, j, m, l>(TensorExpr(a, b));
  }

  /* B(k,m)*A(i,j,k,l) */

  template <class A, class B, class T, class U, int Dim01, int Dim23, int Dim4,
            char i, char j, char k, char l, char m>
  inline const Tensor4_Expr<
    const Ddg_times_Tensor2_2_1<A, B, T, U, Dim01, Dim23, Dim4, i, j, k, l, m>,
    typename promote<T, U>::V, Dim01, Dim01, Dim23, Dim4, i, j, m, l>
  operator*(const Tensor2_Expr<B, U, Dim4, Dim23, k, m> &b,
            const Ddg_Expr<A, T, Dim01, Dim23, i, j, k, l> &a)
  {
    typedef const Ddg_times_Tensor2_2_1<A, B, T, U, Dim01, Dim23, Dim4, i, j,
                                        k, l, m>
      TensorExpr;
    return Tensor4_Expr<TensorExpr, typename promote<T, U>::V, Dim01, Dim01,
                        Dim23, Dim4, i, j, m, l>(TensorExpr(a, b));
  }

  /* A(i,j,k,l)*B(j,m) */

  template <class A, class B, class T, class U, int Dim01, int Dim23, int Dim2,
            char i, char j, char k, char l, char m>
  class Ddg_times_Tensor2_1_0
  {
    const Ddg_Expr<A, T, Dim01, Dim23, i, j, k, l> iterA;
    const Tensor2_Expr<B, U, Dim01, Dim2, j, m> iterB;

    template <int Current_Dim0>
    typename promote<T, U>::V
    eval(const int N1, const int N2, const int N3, const int N4,
         const Number<Current_Dim0> &) const
    {
      return iterA(N1, Current_Dim0 - 1, N3, N4) * iterB(Current_Dim0 - 1, N2)
             + eval(N1, N2, N3, N4, Number<Current_Dim0 - 1>());
    }
    typename promote<T, U>::V eval(const int N1, const int N2, const int N3,
                                   const int N4, const Number<1> &) const
    {
      return iterA(N1, 0, N3, N4) * iterB(0, N2);
    }

  public:
    Ddg_times_Tensor2_1_0(const Ddg_Expr<A, T, Dim01, Dim23, i, j, k, l> &a,
                          const Tensor2_Expr<B, U, Dim01, Dim2, j, m> &b)
        : iterA(a), iterB(b)
    {}
    typename promote<T, U>::V
    operator()(const int N1, const int N2, const int N3, const int N4) const
    {
      return eval(N1, N2, N3, N4, Number<Dim01>());
    }
  };

  template <class A, class B, class T, class U, int Dim01, int Dim23, int Dim2,
            char i, char j, char k, char l, char m>
  inline const Tensor4_Expr<
    const Ddg_times_Tensor2_1_0<A, B, T, U, Dim01, Dim23, Dim2, i, j, k, l, m>,
    typename promote<T, U>::V, Dim01, Dim2, Dim23, Dim23, i, m, k, l>
  operator*(const Ddg_Expr<A, T, Dim01, Dim23, i, j, k, l> &a,
            const Tensor2_Expr<B, U, Dim01, Dim2, j, m> &b)
  {
    typedef const Ddg_times_Tensor2_1_0<A, B, T, U, Dim01, Dim23, Dim2, i, j,
                                        k, l, m>
      TensorExpr;
    return Tensor4_Expr<TensorExpr, typename promote<T, U>::V, Dim01, Dim01,
                        Dim23, Dim2, i, m, k, l>(TensorExpr(a, b));
  }

  /* B(j,m)*A(i,j,k,l) */

  template <class A, class B, class T, class U, int Dim01, int Dim23, int Dim2,
            char i, char j, char k, char l, char m>
  inline const Tensor4_Expr<
    const Ddg_times_Tensor2_1_0<A, B, T, U, Dim01, Dim23, Dim2, i, j, k, l, m>,
    typename promote<T, U>::V, Dim01, Dim2, Dim23, Dim23, i, m, k, l>
  operator*(const Tensor2_Expr<B, U, Dim01, Dim2, j, m> &b,
            const Ddg_Expr<A, T, Dim01, Dim23, i, j, k, l> &a)
  {
    typedef const Ddg_times_Tensor2_1_0<A, B, T, U, Dim01, Dim23, Dim2, i, j,
                                        k, l, m>
      TensorExpr;
    return Tensor4_Expr<TensorExpr, typename promote<T, U>::V, Dim01, Dim01,
                        Dim23, Dim2, i, m, k, l>(TensorExpr(a, b));
  }

  /* A(i,j,k,l)*B(m,j) */

  template <class A, class B, class T, class U, int Dim01, int Dim23, int Dim2,
            char i, char j, char k, char l, char m>
  class Ddg_times_Tensor2_1_1
  {
    const Ddg_Expr<A, T, Dim01, Dim23, i, j, k, l> iterA;
    const Tensor2_Expr<B, U, Dim01, Dim2, m, j> iterB;

    template <int Current_Dim0>
    typename promote<T, U>::V
    eval(const int N1, const int N2, const int N3, const int N4,
         const Number<Current_Dim0> &) const
    {
      return iterA(N1, Current_Dim0 - 1, N3, N4) * iterB(N2, Current_Dim0 - 1)
             + eval(N1, N2, N3, N4, Number<Current_Dim0 - 1>());
    }
    typename promote<T, U>::V eval(const int N1, const int N2, const int N3,
                                   const int N4, const Number<1> &) const
    {
      return iterA(N1, 0, N3, N4) * iterB(N2, 0);
    }

  public:
    Ddg_times_Tensor2_1_1(const Ddg_Expr<A, T, Dim01, Dim23, i, j, k, l> &a,
                          const Tensor2_Expr<B, U, Dim01, Dim2, m, j> &b)
        : iterA(a), iterB(b)
    {}
    typename promote<T, U>::V
    operator()(const int N1, const int N2, const int N3, const int N4) const
    {
      return eval(N1, N2, N3, N4, Number<Dim01>());
    }
  };

  template <class A, class B, class T, class U, int Dim01, int Dim23, int Dim2,
            char i, char j, char k, char l, char m>
  inline const Tensor4_Expr<
    const Ddg_times_Tensor2_1_1<A, B, T, U, Dim01, Dim23, Dim2, i, j, k, l, m>,
    typename promote<T, U>::V, Dim01, Dim2, Dim23, Dim23, i, m, k, l>
  operator*(const Ddg_Expr<A, T, Dim01, Dim23, i, j, k, l> &a,
            const Tensor2_Expr<B, U, Dim01, Dim2, m, j> &b)
  {
    typedef const Ddg_times_Tensor2_1_1<A, B, T, U, Dim01, Dim23, Dim2, i, j,
                                        k, l, m>
      TensorExpr;
    return Tensor4_Expr<TensorExpr, typename promote<T, U>::V, Dim01, Dim2,
                        Dim23, Dim23, i, m, k, l>(TensorExpr(a, b));
  }

  /* B(m,j)*A(i,j,k,l) */

  template <class A, class B, class T, class U, int Dim01, int Dim23, int Dim2,
            char i, char j, char k, char l, char m>
  inline const Tensor4_Expr<
    const Ddg_times_Tensor2_1_1<A, B, T, U, Dim01, Dim23, Dim2, i, j, k, l, m>,
    typename promote<T, U>::V, Dim01, Dim2, Dim23, Dim23, i, m, k, l>
  operator*(const Tensor2_Expr<B, U, Dim01, Dim2, m, j> &b,
            const Ddg_Expr<A, T, Dim01, Dim23, i, j, k, l> &a)
  {
    typedef const Ddg_times_Tensor2_1_1<A, B, T, U, Dim01, Dim23, Dim2, i, j,
                                        k, l, m>
      TensorExpr;
    return Tensor4_Expr<TensorExpr, typename promote<T, U>::V, Dim01, Dim01,
                        Dim23, Dim2, i, m, k, l>(TensorExpr(a, b));
  }

  /* A(i,j,k,l)*B(i,m) */

  template <class A, class B, class T, class U, int Dim01, int Dim23, int Dim1,
            char i, char j, char k, char l, char m>
  class Ddg_times_Tensor2_0_0
  {
    const Ddg_Expr<A, T, Dim01, Dim23, i, j, k, l> iterA;
    const Tensor2_Expr<B, U, Dim01, Dim1, i, m> iterB;

    template <int Current_Dim0>
    typename promote<T, U>::V
    eval(const int N1, const int N2, const int N3, const int N4,
         const Number<Current_Dim0> &) const
    {
      return iterA(Current_Dim0 - 1, N2, N3, N4) * iterB(Current_Dim0 - 1, N1)
             + eval(N1, N2, N3, N4, Number<Current_Dim0 - 1>());
    }
    typename promote<T, U>::V eval(const int N1, const int N2, const int N3,
                                   const int N4, const Number<1> &) const
    {
      return iterA(0, N2, N3, N4) * iterB(0, N1);
    }

  public:
    Ddg_times_Tensor2_0_0(const Ddg_Expr<A, T, Dim01, Dim23, i, j, k, l> &a,
                          const Tensor2_Expr<B, U, Dim01, Dim1, i, m> &b)
        : iterA(a), iterB(b)
    {}
    typename promote<T, U>::V
    operator()(const int N1, const int N2, const int N3, const int N4) const
    {
      return eval(N1, N2, N3, N4, Number<Dim01>());
    }
  };

  template <class A, class B, class T, class U, int Dim01, int Dim23, int Dim1,
            char i, char j, char k, char l, char m>
  inline const Tensor4_Expr<
    const Ddg_times_Tensor2_0_0<A, B, T, U, Dim01, Dim23, Dim1, i, j, k, l, m>,
    typename promote<T, U>::V, Dim1, Dim01, Dim23, Dim23, m, j, k, l>
  operator*(const Ddg_Expr<A, T, Dim01, Dim23, i, j, k, l> &a,
            const Tensor2_Expr<B, U, Dim01, Dim1, i, m> &b)
  {
    typedef const Ddg_times_Tensor2_0_0<A, B, T, U, Dim01, Dim23, Dim1, i, j,
                                        k, l, m>
      TensorExpr;
    return Tensor4_Expr<TensorExpr, typename promote<T, U>::V, Dim1, Dim01,
                        Dim23, Dim23, m, j, k, l>(TensorExpr(a, b));
  }

  /* B(i,m)*A(i,j,k,l) */

  template <class A, class B, class T, class U, int Dim01, int Dim23, int Dim1,
            char i, char j, char k, char l, char m>
  inline const Tensor4_Expr<
    const Ddg_times_Tensor2_0_0<A, B, T, U, Dim01, Dim23, Dim1, i, j, k, l, m>,
    typename promote<T, U>::V, Dim01, Dim01, Dim23, Dim1, m, j, k, l>
  operator*(const Tensor2_Expr<B, U, Dim01, Dim1, i, m> &b,
            const Ddg_Expr<A, T, Dim01, Dim23, i, j, k, l> &a)
  {
    typedef const Ddg_times_Tensor2_0_0<A, B, T, U, Dim01, Dim23, Dim1, i, j,
                                        k, l, m>
      TensorExpr;
    return Tensor4_Expr<TensorExpr, typename promote<T, U>::V, Dim1, Dim01,
                        Dim23, Dim23, m, j, k, l>(TensorExpr(a, b));
  }

  /* A(i,j,k,l)*B(m,i) */

  template <class A, class B, class T, class U, int Dim01, int Dim23, int Dim1,
            char i, char j, char k, char l, char m>
  class Ddg_times_Tensor2_0_1
  {
    const Ddg_Expr<A, T, Dim01, Dim23, i, j, k, l> iterA;
    const Tensor2_Expr<B, U, Dim01, Dim1, m, i> iterB;

    template <int Current_Dim0>
    typename promote<T, U>::V
    eval(const int N1, const int N2, const int N3, const int N4,
         const Number<Current_Dim0> &) const
    {
      return iterA(Current_Dim0 - 1, N2, N3, N4) * iterB(N1, Current_Dim0 - 1)
             + eval(N1, N2, N3, N4, Number<Current_Dim0 - 1>());
    }
    typename promote<T, U>::V eval(const int N1, const int N2, const int N3,
                                   const int N4, const Number<1> &) const
    {
      return iterA(0, N2, N3, N4) * iterB(N1, 0);
    }

  public:
    Ddg_times_Tensor2_0_1(const Ddg_Expr<A, T, Dim01, Dim23, i, j, k, l> &a,
                          const Tensor2_Expr<B, U, Dim01, Dim1, m, i> &b)
        : iterA(a), iterB(b)
    {}
    typename promote<T, U>::V
    operator()(const int N1, const int N2, const int N3, const int N4) const
    {
      return eval(N1, N2, N3, N4, Number<Dim01>());
    }
  };

  template <class A, class B, class T, class U, int Dim01, int Dim23, int Dim1,
            char i, char j, char k, char l, char m>
  inline const Tensor4_Expr<
    const Ddg_times_Tensor2_0_1<A, B, T, U, Dim01, Dim23, Dim1, i, j, k, l, m>,
    typename promote<T, U>::V, Dim1, Dim01, Dim23, Dim23, m, j, k, l>
  operator*(const Ddg_Expr<A, T, Dim01, Dim23, i, j, k, l> &a,
            const Tensor2_Expr<B, U, Dim01, Dim1, m, i> &b)
  {
    typedef const Ddg_times_Tensor2_0_1<A, B, T, U, Dim01, Dim23, Dim1, i, j,
                                        k, l, m>
      TensorExpr;
    return Tensor4_Expr<TensorExpr, typename promote<T, U>::V, Dim01, Dim01,
                        Dim23, Dim1, m, j, k, l>(TensorExpr(a, b));
  }

  /* B(i,m)*A(i,j,k,l) */

  template <class A, class B, class T, class U, int Dim01, int Dim23, int Dim1,
            char i, char j, char k, char l, char m>
  inline const Tensor4_Expr<
    const Ddg_times_Tensor2_0_1<A, B, T, U, Dim01, Dim23, Dim1, i, j, k, l, m>,
    typename promote<T, U>::V, Dim1, Dim01, Dim23, Dim23, m, j, k, l>
  operator*(const Tensor2_Expr<B, U, Dim01, Dim1, m, i> &b,
            const Ddg_Expr<A, T, Dim01, Dim23, i, j, k, l> &a)
  {
    typedef const Ddg_times_Tensor2_0_1<A, B, T, U, Dim01, Dim23, Dim1, i, j,
                                        k, l, m>
      TensorExpr;
    return Tensor4_Expr<TensorExpr, typename promote<T, U>::V, Dim1, Dim01,
                        Dim23, Dim23, m, j, k, l>(TensorExpr(a, b));
  }
}
