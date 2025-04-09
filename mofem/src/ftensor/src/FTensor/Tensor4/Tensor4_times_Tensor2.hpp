/* This file has all of the declarations for expressions like
   Tensor4*Tensor2 and Tensor2*Tensor4, yielding a
   Tensor2. */

#pragma once

namespace FTensor
{
  /* A(i,j,k,l)*B(k,l) -> Tensor2 */

  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2,
            int Dim3, char i, char j, char k, char l>
  class Tensor4_times_Tensor2_23
  {
    const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> iterA;
    const Tensor2_Expr<B, U, Dim2, Dim3, k, l> iterB;

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
             + eval(N1, N2, Number<Dim2>(), Number<Current_Dim1 - 1>());
    }
    typename promote<T, U>::V eval(const int N1, const int N2,
                                   const Number<1> &, const Number<1> &) const
    {
      return iterA(N1, N2, 0, 0) * iterB(0, 0);
    }

  public:
    typename promote<T, U>::V operator()(const int N1, const int N2) const
    {
      return eval(N1, N2, Number<Dim2>(), Number<Dim3>());
    }

    Tensor4_times_Tensor2_23(
      const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a,
      const Tensor2_Expr<B, U, Dim2, Dim3, k, l> &b)
        : iterA(a), iterB(b)
    {}
  };

  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2,
            int Dim3, char i, char j, char k, char l>
  inline Tensor2_Expr<const Tensor4_times_Tensor2_23<A, B, T, U, Dim0, Dim1,
                                                     Dim2, Dim3, i, j, k, l>,
                      typename promote<T, U>::V, Dim0, Dim1, i, j>
  operator*(const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a,
            const Tensor2_Expr<B, U, Dim2, Dim3, k, l> &b)
  {
    typedef const Tensor4_times_Tensor2_23<A, B, T, U, Dim0, Dim1, Dim2, Dim3,
                                           i, j, k, l>
      TensorExpr;
    return Tensor2_Expr<TensorExpr, typename promote<T, U>::V, Dim0, Dim1, i,
                        j>(TensorExpr(a, b));
  }

  /* B(k,l)*A(i,j,k,l) -> Tensor2 */

  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2,
            int Dim3, char i, char j, char k, char l>
  inline Tensor2_Expr<const Tensor4_times_Tensor2_23<A, B, T, U, Dim0, Dim1,
                                                     Dim2, Dim3, i, j, k, l>,
                      typename promote<T, U>::V, Dim0, Dim1, i, j>
  operator*(const Tensor2_Expr<B, U, Dim2, Dim3, k, l> &b,
            const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a)
  {
    typedef const Tensor4_times_Tensor2_23<A, B, T, U, Dim0, Dim1, Dim2, Dim3,
                                           i, j, k, l>
      TensorExpr;
    return Tensor2_Expr<TensorExpr, typename promote<T, U>::V, Dim0, Dim1, i,
                        j>(TensorExpr(a, b));
  }

  /* A(i,j,k,l)*B(l,k) -> Tensor2 */

  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2,
            int Dim3, char i, char j, char k, char l>
  class Tensor4_times_Tensor2_32
  {
    const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> iterA;
    const Tensor2_Expr<B, U, Dim3, Dim2, l, k> iterB;

    template <int Current_Dim0, int Current_Dim1>
    typename promote<T, U>::V
    eval(const int N1, const int N2, const Number<Current_Dim0> &,
         const Number<Current_Dim1> &) const
    {
      return iterA(N1, N2, Current_Dim0 - 1, Current_Dim1 - 1)
               * iterB(Current_Dim1 - 1, Current_Dim0 - 1)
             + eval(N1, N2, Number<Current_Dim0 - 1>(),
                    Number<Current_Dim1>());
    }
    template <int Current_Dim1>
    typename promote<T, U>::V
    eval(const int N1, const int N2, const Number<1> &,
         const Number<Current_Dim1> &) const
    {
      return iterA(N1, N2, 0, Current_Dim1 - 1) * iterB(Current_Dim1 - 1, 0)
             + eval(N1, N2, Number<Dim2>(), Number<Current_Dim1 - 1>());
    }
    typename promote<T, U>::V eval(const int N1, const int N2,
                                   const Number<1> &, const Number<1> &) const
    {
      return iterA(N1, N2, 0, 0) * iterB(0, 0);
    }

  public:
    typename promote<T, U>::V operator()(const int N1, const int N2) const
    {
      return eval(N1, N2, Number<Dim2>(), Number<Dim3>());
    }

    Tensor4_times_Tensor2_32(
      const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a,
      const Tensor2_Expr<B, U, Dim3, Dim2, l, k> &b)
        : iterA(a), iterB(b)
    {}
  };

  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2,
            int Dim3, char i, char j, char k, char l>
  inline Tensor2_Expr<const Tensor4_times_Tensor2_32<A, B, T, U, Dim0, Dim1,
                                                     Dim2, Dim3, i, j, k, l>,
                      typename promote<T, U>::V, Dim0, Dim1, i, j>
  operator*(const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a,
            const Tensor2_Expr<B, U, Dim3, Dim2, l, k> &b)
  {
    typedef const Tensor4_times_Tensor2_32<A, B, T, U, Dim0, Dim1, Dim2, Dim3,
                                           i, j, k, l>
      TensorExpr;
    return Tensor2_Expr<TensorExpr, typename promote<T, U>::V, Dim0, Dim1, i,
                        j>(TensorExpr(a, b));
  }

  /* B(l,k)*A(i,j,k,l) -> Tensor2 */

  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2,
            int Dim3, char i, char j, char k, char l>
  inline Tensor2_Expr<const Tensor4_times_Tensor2_32<A, B, T, U, Dim0, Dim1,
                                                     Dim2, Dim3, i, j, k, l>,
                      typename promote<T, U>::V, Dim0, Dim1, i, j>
  operator*(const Tensor2_Expr<B, U, Dim3, Dim2, l, k> &b,
            const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a)
  {
    typedef const Tensor4_times_Tensor2_32<A, B, T, U, Dim0, Dim1, Dim2, Dim3,
                                           i, j, k, l>
      TensorExpr;
    return Tensor2_Expr<TensorExpr, typename promote<T, U>::V, Dim0, Dim1, i,
                        j>(TensorExpr(a, b));
  }

  /* A(i,j,k,l)*B(i,l) -> Tensor2 */

  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2,
            int Dim3, char i, char j, char k, char l>
  class Tensor4_times_Tensor2_03
  {
    const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> iterA;
    const Tensor2_Expr<B, U, Dim0, Dim3, i, l> iterB;

    template <int Current_Dim0, int Current_Dim1>
    typename promote<T, U>::V
    eval(const int N1, const int N2, const Number<Current_Dim0> &,
         const Number<Current_Dim1> &) const
    {
      return iterA(Current_Dim0 - 1, N1, N2, Current_Dim1 - 1)
               * iterB(Current_Dim0 - 1, Current_Dim1 - 1)
             + eval(N1, N2, Number<Current_Dim0 - 1>(),
                    Number<Current_Dim1>());
    }
    template <int Current_Dim1>
    typename promote<T, U>::V
    eval(const int N1, const int N2, const Number<1> &,
         const Number<Current_Dim1> &) const
    {
      return iterA(0, N1, N2, Current_Dim1 - 1) * iterB(0, Current_Dim1 - 1)
             + eval(N1, N2, Number<Dim0>(), Number<Current_Dim1 - 1>());
    }
    typename promote<T, U>::V eval(const int N1, const int N2,
                                   const Number<1> &, const Number<1> &) const
    {
      return iterA(0, N1, N2, 0) * iterB(0, 0);
    }

  public:
    typename promote<T, U>::V operator()(const int N1, const int N2) const
    {
      return eval(N1, N2, Number<Dim0>(), Number<Dim3>());
    }

    Tensor4_times_Tensor2_03(
      const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a,
      const Tensor2_Expr<B, U, Dim0, Dim3, i, l> &b)
        : iterA(a), iterB(b)
    {}
  };

  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2,
            int Dim3, char i, char j, char k, char l>
  inline Tensor2_Expr<const Tensor4_times_Tensor2_03<A, B, T, U, Dim0, Dim1,
                                                     Dim2, Dim3, i, j, k, l>,
                      typename promote<T, U>::V, Dim1, Dim2, j, k>
  operator*(const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a,
            const Tensor2_Expr<B, U, Dim0, Dim3, i, l> &b)
  {
    typedef const Tensor4_times_Tensor2_03<A, B, T, U, Dim0, Dim1, Dim2, Dim3,
                                           i, j, k, l>
      TensorExpr;
    return Tensor2_Expr<TensorExpr, typename promote<T, U>::V, Dim1, Dim2, j,
                        k>(TensorExpr(a, b));
  }

  /* B(i,l)*A(i,j,k,l) -> Tensor2 */

  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2,
            int Dim3, char i, char j, char k, char l>
  inline Tensor2_Expr<const Tensor4_times_Tensor2_03<A, B, T, U, Dim0, Dim1,
                                                     Dim2, Dim3, i, j, k, l>,
                      typename promote<T, U>::V, Dim1, Dim2, j, k>
  operator*(const Tensor2_Expr<B, U, Dim0, Dim3, i, l> &b,
            const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a)
  {
    typedef const Tensor4_times_Tensor2_03<A, B, T, U, Dim0, Dim1, Dim2, Dim3,
                                           i, j, k, l>
      TensorExpr;
    return Tensor2_Expr<TensorExpr, typename promote<T, U>::V, Dim1, Dim2, j,
                        k>(TensorExpr(a, b));
  }

  /* A(i,j,k,l)*B(l,i) -> Tensor2 */

  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2,
            int Dim3, char i, char j, char k, char l>
  class Tensor4_times_Tensor2_30
  {
    const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> iterA;
    const Tensor2_Expr<B, U, Dim3, Dim0, l, i> iterB;

    template <int Current_Dim0, int Current_Dim1>
    typename promote<T, U>::V
    eval(const int N1, const int N2, const Number<Current_Dim0> &,
         const Number<Current_Dim1> &) const
    {
      return iterA(Current_Dim0 - 1, N1, N2, Current_Dim1 - 1)
               * iterB(Current_Dim1 - 1, Current_Dim0 - 1)
             + eval(N1, N2, Number<Current_Dim0 - 1>(),
                    Number<Current_Dim1>());
    }
    template <int Current_Dim1>
    typename promote<T, U>::V
    eval(const int N1, const int N2, const Number<1> &,
         const Number<Current_Dim1> &) const
    {
      return iterA(0, N1, N2, Current_Dim1 - 1) * iterB(Current_Dim1 - 1, 0)
             + eval(N1, N2, Number<Dim0>(), Number<Current_Dim1 - 1>());
    }
    typename promote<T, U>::V eval(const int N1, const int N2,
                                   const Number<1> &, const Number<1> &) const
    {
      return iterA(0, N1, N2, 0) * iterB(0, 0);
    }

  public:
    typename promote<T, U>::V operator()(const int N1, const int N2) const
    {
      return eval(N1, N2, Number<Dim0>(), Number<Dim3>());
    }

    Tensor4_times_Tensor2_30(
      const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a,
      const Tensor2_Expr<B, U, Dim3, Dim0, l, i> &b)
        : iterA(a), iterB(b)
    {}
  };

  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2,
            int Dim3, char i, char j, char k, char l>
  inline Tensor2_Expr<const Tensor4_times_Tensor2_30<A, B, T, U, Dim0, Dim1,
                                                     Dim2, Dim3, i, j, k, l>,
                      typename promote<T, U>::V, Dim1, Dim2, j, k>
  operator*(const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a,
            const Tensor2_Expr<B, U, Dim3, Dim0, l, i> &b)
  {
    typedef const Tensor4_times_Tensor2_30<A, B, T, U, Dim0, Dim1, Dim2, Dim3,
                                           i, j, k, l>
      TensorExpr;
    return Tensor2_Expr<TensorExpr, typename promote<T, U>::V, Dim1, Dim2, j,
                        k>(TensorExpr(a, b));
  }

  /* B(l,i)*A(i,j,k,l) -> Tensor2 */

  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2,
            int Dim3, char i, char j, char k, char l>
  inline Tensor2_Expr<const Tensor4_times_Tensor2_30<A, B, T, U, Dim0, Dim1,
                                                     Dim2, Dim3, i, j, k, l>,
                      typename promote<T, U>::V, Dim1, Dim2, j, k>
  operator*(const Tensor2_Expr<B, U, Dim3, Dim0, l, i> &b,
            const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a)
  {
    typedef const Tensor4_times_Tensor2_30<A, B, T, U, Dim0, Dim1, Dim2, Dim3,
                                           i, j, k, l>
      TensorExpr;
    return Tensor2_Expr<TensorExpr, typename promote<T, U>::V, Dim1, Dim2, j,
                        k>(TensorExpr(a, b));
  }

  /* A(i,j,k,l)*B(j,l) -> Tensor2 */

  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2,
            int Dim3, char i, char j, char k, char l>
  class Tensor4_times_Tensor2_13
  {
    const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> iterA;
    const Tensor2_Expr<B, U, Dim2, Dim3, j, l> iterB;

    template <int Current_Dim0, int Current_Dim1>
    typename promote<T, U>::V
    eval(const int N1, const int N2, const Number<Current_Dim0> &,
         const Number<Current_Dim1> &) const
    {
      return iterA(N1, Current_Dim0 - 1, N2, Current_Dim1 - 1)
               * iterB(Current_Dim0 - 1, Current_Dim1 - 1)
             + eval(N1, N2, Number<Current_Dim0 - 1>(),
                    Number<Current_Dim1>());
    }
    template <int Current_Dim1>
    typename promote<T, U>::V
    eval(const int N1, const int N2, const Number<1> &,
         const Number<Current_Dim1> &) const
    {
      return iterA(N1, 0, N2, Current_Dim1 - 1) * iterB(0, Current_Dim1 - 1)
             + eval(N1, N2, Number<Dim0>(), Number<Current_Dim1 - 1>());
    }
    typename promote<T, U>::V eval(const int N1, const int N2,
                                   const Number<1> &, const Number<1> &) const
    {
      return iterA(N1, 0, N2, 0) * iterB(0, 0);
    }

  public:
    typename promote<T, U>::V operator()(const int N1, const int N2) const
    {
      return eval(N1, N2, Number<Dim0>(), Number<Dim3>());
    }

    Tensor4_times_Tensor2_13(
      const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a,
      const Tensor2_Expr<B, U, Dim2, Dim3, j, l> &b)
        : iterA(a), iterB(b)
    {}
  };

  /* B(j,l)*A(i,j,k,l) -> Tensor2 */

  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2,
            int Dim3, char i, char j, char k, char l>
  inline Tensor2_Expr<const Tensor4_times_Tensor2_13<A, B, T, U, Dim0, Dim1,
                                                     Dim2, Dim3, i, j, k, l>,
                      typename promote<T, U>::V, Dim0, Dim2, i, k>
  operator*(const Tensor2_Expr<B, U, Dim1, Dim3, j, l> &b,
            const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a)
  {
    typedef const Tensor4_times_Tensor2_13<A, B, T, U, Dim0, Dim1, Dim2, Dim3,
                                           i, j, k, l>
      TensorExpr;
    return Tensor2_Expr<TensorExpr, typename promote<T, U>::V, Dim1, Dim2, i,
                        k>(TensorExpr(a, b));
  }

  /* This file has all of the declarations for expressions like
     Tensor4*Tensor2 and Tensor2*Tensor4, yielding a
     Tensor4. */

  // TODO: Check dimensions could be errors

  /* A(i,j,k,l)*B(l,m) -> Tensor4 */

  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2,
            int Dim3, int Dim4, char i, char j, char k, char l, char m>
  class Tensor4_times_Tensor2_3_1
  {
    const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> iterA;
    const Tensor2_Expr<B, U, Dim3, Dim4, l, m> iterB;

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
    Tensor4_times_Tensor2_3_1(
      const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a,
      const Tensor2_Expr<B, U, Dim3, Dim4, l, m> &b)
        : iterA(a), iterB(b)
    {}
    typename promote<T, U>::V
    operator()(const int N1, const int N2, const int N3, const int N4) const
    {
      return eval(N1, N2, N3, N4, Number<Dim3>());
    }
  };

  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2,
            int Dim3, int Dim4, char i, char j, char k, char l, char m>
  inline const Tensor4_Expr<
    const Tensor4_times_Tensor2_3_1<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4,
                                    i, j, k, l, m>,
    typename promote<T, U>::V, Dim0, Dim1, Dim2, Dim4, i, j, k, m>
  operator*(const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a,
            const Tensor2_Expr<B, U, Dim3, Dim4, l, m> &b)
  {
    typedef const Tensor4_times_Tensor2_3_1<A, B, T, U, Dim0, Dim1, Dim2, Dim3,
                                            Dim4, i, j, k, l, m>
      TensorExpr;
    return Tensor4_Expr<TensorExpr, typename promote<T, U>::V, Dim0, Dim1,
                        Dim2, Dim4, i, j, k, m>(TensorExpr(a, b));
  }

  /* B(l,m)*A(i,j,k,l) -> Tensor4 */

  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2,
            int Dim3, int Dim4, char i, char j, char k, char l, char m>
  inline const Tensor4_Expr<
    const Tensor4_times_Tensor2_3_1<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4,
                                    i, j, k, l, m>,
    typename promote<T, U>::V, Dim0, Dim1, Dim2, Dim4, i, j, k, m>
  operator*(const Tensor2_Expr<B, U, Dim3, Dim4, l, m> &b,
            const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a)
  {
    typedef const Tensor4_times_Tensor2_3_1<A, B, T, U, Dim0, Dim1, Dim2, Dim3,
                                            Dim4, i, j, k, l, m>
      TensorExpr;
    return Tensor4_Expr<TensorExpr, typename promote<T, U>::V, Dim0, Dim1,
                        Dim2, Dim4, i, j, k, m>(TensorExpr(a, b));
  }

  /* A(i,j,k,l)*B(m,l) -> Tensor4 */

  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2,
            int Dim3, int Dim4, char i, char j, char k, char l, char m>
  class Tensor4_times_Tensor2_3_0
  {
    const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> iterA;
    const Tensor2_Expr<B, U, Dim4, Dim3, m, l> iterB;

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
    Tensor4_times_Tensor2_3_0(
      const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a,
      const Tensor2_Expr<B, U, Dim4, Dim3, m, l> &b)
        : iterA(a), iterB(b)
    {}
    typename promote<T, U>::V
    operator()(const int N1, const int N2, const int N3, const int N4) const
    {
      return eval(N1, N2, N3, N4, Number<Dim3>());
    }
  };

  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2,
            int Dim3, int Dim4, char i, char j, char k, char l, char m>
  inline const Tensor4_Expr<
    const Tensor4_times_Tensor2_3_0<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4,
                                    i, j, k, l, m>,
    typename promote<T, U>::V, Dim0, Dim1, Dim2, Dim4, i, j, k, m>
  operator*(const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a,
            const Tensor2_Expr<B, U, Dim4, Dim3, m, l> &b)
  {
    typedef const Tensor4_times_Tensor2_3_0<A, B, T, U, Dim0, Dim1, Dim2, Dim3,
                                            Dim4, i, j, k, l, m>
      TensorExpr;
    return Tensor4_Expr<TensorExpr, typename promote<T, U>::V, Dim0, Dim1,
                        Dim2, Dim4, i, j, k, m>(TensorExpr(a, b));
  }

  /* B(m,l)*A(i,j,k,l) -> Tensor4 */

  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2,
            int Dim3, int Dim4, char i, char j, char k, char l, char m>
  inline const Tensor4_Expr<
    const Tensor4_times_Tensor2_3_0<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4,
                                    i, j, k, l, m>,
    typename promote<T, U>::V, Dim0, Dim1, Dim2, Dim4, i, j, k, m>
  operator*(const Tensor2_Expr<B, U, Dim4, Dim3, m, l> &b,
            const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a)
  {
    typedef const Tensor4_times_Tensor2_3_0<A, B, T, U, Dim0, Dim1, Dim2, Dim3,
                                            Dim4, i, j, k, l, m>
      TensorExpr;
    return Tensor4_Expr<TensorExpr, typename promote<T, U>::V, Dim0, Dim1,
                        Dim2, Dim4, i, j, k, m>(TensorExpr(a, b));
  }

  /* A(i,j,k,l)*B(j,m) -> Tensor4 */

  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2,
            int Dim3, int Dim4, char i, char j, char k, char l, char m>
  class Tensor4_times_Tensor2_1_0
  {
    const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> iterA;
    const Tensor2_Expr<B, U, Dim1, Dim4, j, m> iterB;

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
    Tensor4_times_Tensor2_1_0(
      const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a,
      const Tensor2_Expr<B, U, Dim1, Dim4, j, m> &b)
        : iterA(a), iterB(b)
    {}
    typename promote<T, U>::V
    operator()(const int N1, const int N2, const int N3, const int N4) const
    {
      return eval(N1, N2, N3, N4, Number<Dim1>());
    }
  };

  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2,
            int Dim3, int Dim4, char i, char j, char k, char l, char m>
  inline const Tensor4_Expr<
    const Tensor4_times_Tensor2_1_0<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4,
                                    i, j, k, l, m>,
    typename promote<T, U>::V, Dim0, Dim4, Dim2, Dim3, i, m, k, l>
  operator*(const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a,
            const Tensor2_Expr<B, U, Dim1, Dim4, j, m> &b)
  {
    typedef const Tensor4_times_Tensor2_1_0<A, B, T, U, Dim0, Dim1, Dim2, Dim3,
                                            Dim4, i, j, k, l, m>
      TensorExpr;
    return Tensor4_Expr<TensorExpr, typename promote<T, U>::V, Dim0, Dim1,
                        Dim4, Dim3, i, m, k, l>(TensorExpr(a, b));
  }

  /* B(j,m)*A(i,j,k,l) -> Tensor4 */

  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2,
            int Dim3, int Dim4, char i, char j, char k, char l, char m>
  inline const Tensor4_Expr<
    const Tensor4_times_Tensor2_1_0<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4,
                                    i, j, k, l, m>,
    typename promote<T, U>::V, Dim0, Dim4, Dim2, Dim3, i, m, k, l>
  operator*(const Tensor2_Expr<B, U, Dim1, Dim4, j, m> &b,
            const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a)
  {
    typedef const Tensor4_times_Tensor2_1_0<A, B, T, U, Dim0, Dim1, Dim2, Dim3,
                                            Dim4, i, j, k, l, m>
      TensorExpr;
    return Tensor4_Expr<TensorExpr, typename promote<T, U>::V, Dim0, Dim1,
                        Dim4, Dim3, i, m, k, l>(TensorExpr(a, b));
  }

  /* A(i,j,k,l)*B(m,j) -> Tensor4 */

  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2,
            int Dim3, int Dim4, char i, char j, char k, char l, char m>
  class Tensor4_times_Tensor2_1_1
  {
    const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> iterA;
    const Tensor2_Expr<B, U, Dim4, Dim1, m, j> iterB;

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
    Tensor4_times_Tensor2_1_1(
      const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a,
      const Tensor2_Expr<B, U, Dim4, Dim1, m, j> &b)
        : iterA(a), iterB(b)
    {}
    typename promote<T, U>::V
    operator()(const int N1, const int N2, const int N3, const int N4) const
    {
      return eval(N1, N2, N3, N4, Number<Dim1>());
    }
  };

  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2,
            int Dim3, int Dim4, char i, char j, char k, char l, char m>
  inline const Tensor4_Expr<
    const Tensor4_times_Tensor2_1_1<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4,
                                    i, j, k, l, m>,
    typename promote<T, U>::V, Dim0, Dim4, Dim2, Dim3, i, m, k, l>
  operator*(const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a,
            const Tensor2_Expr<B, U, Dim4, Dim1, m, j> &b)
  {
    typedef const Tensor4_times_Tensor2_1_1<A, B, T, U, Dim0, Dim1, Dim2, Dim3,
                                            Dim4, i, j, k, l, m>
      TensorExpr;
    return Tensor4_Expr<TensorExpr, typename promote<T, U>::V, Dim0, Dim4,
                        Dim2, Dim3, i, m, k, l>(TensorExpr(a, b));
  }

  /* B(m,j)*A(i,j,k,l) -> Tensor4 */

  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2,
            int Dim3, int Dim4, char i, char j, char k, char l, char m>
  inline const Tensor4_Expr<
    const Tensor4_times_Tensor2_1_1<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4,
                                    i, j, k, l, m>,
    typename promote<T, U>::V, Dim0, Dim4, Dim2, Dim3, i, m, k, l>
  operator*(const Tensor2_Expr<B, U, Dim4, Dim1, m, j> &b,
            const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a)
  {
    typedef const Tensor4_times_Tensor2_1_1<A, B, T, U, Dim0, Dim1, Dim2, Dim3,
                                            Dim4, i, j, k, l, m>
      TensorExpr;
    return Tensor4_Expr<TensorExpr, typename promote<T, U>::V, Dim0, Dim4,
                        Dim2, Dim3, i, m, k, l>(TensorExpr(a, b));
  }

  /* A(i,j,k,l)*B(i,m) -> Tensor4 */

  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2,
            int Dim3, int Dim4, char i, char j, char k, char l, char m>
  class Tensor4_times_Tensor2_0_0
  {
    const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> iterA;
    const Tensor2_Expr<B, U, Dim0, Dim4, i, m> iterB;

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
    Tensor4_times_Tensor2_0_0(
      const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a,
      const Tensor2_Expr<B, U, Dim0, Dim4, i, m> &b)
        : iterA(a), iterB(b)
    {}
    typename promote<T, U>::V
    operator()(const int N1, const int N2, const int N3, const int N4) const
    {
      return eval(N1, N2, N3, N4, Number<Dim0>());
    }
  };

  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2,
            int Dim3, int Dim4, char i, char j, char k, char l, char m>
  inline const Tensor4_Expr<
    const Tensor4_times_Tensor2_0_0<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4,
                                    i, j, k, l, m>,
    typename promote<T, U>::V, Dim1, Dim1, Dim2, Dim3, m, j, k, l>
  operator*(const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a,
            const Tensor2_Expr<B, U, Dim0, Dim4, i, m> &b)
  {
    typedef const Tensor4_times_Tensor2_0_0<A, B, T, U, Dim0, Dim1, Dim2, Dim3,
                                            Dim4, i, j, k, l, m>
      TensorExpr;
    return Tensor4_Expr<TensorExpr, typename promote<T, U>::V, Dim4, Dim1,
                        Dim2, Dim3, m, j, k, l>(TensorExpr(a, b));
  }

  /* B(i,m)*A(i,j,k,l) -> Tensor4 */

  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2,
            int Dim3, int Dim4, char i, char j, char k, char l, char m>
  inline const Tensor4_Expr<
    const Tensor4_times_Tensor2_0_0<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4,
                                    i, j, k, l, m>,
    typename promote<T, U>::V, Dim1, Dim1, Dim2, Dim3, m, j, k, l>
  operator*(const Tensor2_Expr<B, U, Dim0, Dim4, i, m> &b,
            const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a)
  {
    typedef const Tensor4_times_Tensor2_0_0<A, B, T, U, Dim0, Dim1, Dim2, Dim3,
                                            Dim4, i, j, k, l, m>
      TensorExpr;
    return Tensor4_Expr<TensorExpr, typename promote<T, U>::V, Dim4, Dim1,
                        Dim2, Dim3, m, j, k, l>(TensorExpr(a, b));
  }

  /* A(i,j,k,l)*B(m,i) -> Tensor4 */

  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2,
            int Dim3, int Dim4, char i, char j, char k, char l, char m>
  class Tensor4_times_Tensor2_0_1
  {
    const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> iterA;
    const Tensor2_Expr<B, U, Dim4, Dim0, m, i> iterB;

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
    Tensor4_times_Tensor2_0_1(
      const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a,
      const Tensor2_Expr<B, U, Dim4, Dim1, m, i> &b)
        : iterA(a), iterB(b)
    {}
    typename promote<T, U>::V
    operator()(const int N1, const int N2, const int N3, const int N4) const
    {
      return eval(N1, N2, N3, N4, Number<Dim0>());
    }
  };

  // FIXME: Template is ambiguous
  // template<
  // class A, class B, class T, class U, int Dim0,int Dim1, int Dim2,int Dim3,
  // int Dim4,char i, char j, char k, char l,char m> inline const Tensor4_Expr
  // <const
  // Tensor4_times_Tensor2_0_1<A,B,T,U,Dim0,Dim1,Dim2,Dim3,Dim4,i,j,k,l,m>,typename
  // promote<T,U>::V,Dim4,Dim1,Dim2,Dim3,m,j,k,l> operator*(
  //   const Tensor4_Expr<A,T,Dim0,Dim1,Dim2,Dim3,i,j,k,l> &a,
  //   const Tensor2_Expr<B,U,Dim4,Dim1,m,i> &b
  // ) {
  //   typedef const
  //   Tensor4_times_Tensor2_0_1<A,B,T,U,Dim0,Dim1,Dim2,Dim3,Dim4,i,j,k,l,m>
  //   TensorExpr; return Tensor4_Expr<TensorExpr,typename
  //   promote<T,U>::V,Dim4,Dim1,Dim2,Dim3,m,j,k,l> (TensorExpr(a,b));
  // }

  /* B(i,m)*A(i,j,k,l) -> Tensor4 */

  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2,
            int Dim3, int Dim4, char i, char j, char k, char l, char m>
  inline const Tensor4_Expr<
    const Tensor4_times_Tensor2_0_1<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4,
                                    i, j, k, l, m>,
    typename promote<T, U>::V, Dim4, Dim1, Dim2, Dim3, m, j, k, l>
  operator*(const Tensor2_Expr<B, U, Dim4, Dim1, m, i> &b,
            const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a)
  {
    typedef const Tensor4_times_Tensor2_0_1<A, B, T, U, Dim0, Dim1, Dim2, Dim3,
                                            Dim4, i, j, k, l, m>
      TensorExpr;
    return Tensor4_Expr<TensorExpr, typename promote<T, U>::V, Dim4, Dim1,
                        Dim2, Dim3, m, j, k, l>(TensorExpr(a, b));
  }

  /* A(i,j,k,l)*B(k,m) -> Tensor4 */

  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2,
            int Dim3, int Dim4, char i, char j, char k, char l, char m>
  class Tensor4_times_Tensor2_2_0
  {
    const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> iterA;
    const Tensor2_Expr<B, U, Dim2, Dim4, k, m> iterB;

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
    Tensor4_times_Tensor2_2_0(
      const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a,
      const Tensor2_Expr<B, U, Dim2, Dim4, k, m> &b)
        : iterA(a), iterB(b)
    {}
    typename promote<T, U>::V
    operator()(const int N1, const int N2, const int N3, const int N4) const
    {
      return eval(N1, N2, N3, N4, Number<Dim3>());
    }
  };

  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2,
            int Dim3, int Dim4, char i, char j, char k, char l, char m>
  inline const Tensor4_Expr<
    const Tensor4_times_Tensor2_2_0<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4,
                                    i, j, k, l, m>,
    typename promote<T, U>::V, Dim0, Dim1, Dim4, Dim3, i, j, m, l>
  operator*(const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a,
            const Tensor2_Expr<B, U, Dim2, Dim4, k, m> &b)
  {
    typedef const Tensor4_times_Tensor2_2_0<A, B, T, U, Dim0, Dim1, Dim2, Dim3,
                                            Dim4, i, j, k, l, m>
      TensorExpr;
    return Tensor4_Expr<TensorExpr, typename promote<T, U>::V, Dim0, Dim1,
                        Dim4, Dim3, i, j, m, l>(TensorExpr(a, b));
  }

  /* B(k,m)*A(i,j,k,l) -> Tensor4 */

  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2,
            int Dim3, int Dim4, char i, char j, char k, char l, char m>
  inline const Tensor4_Expr<
    const Tensor4_times_Tensor2_2_0<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4,
                                    i, j, k, l, m>,
    typename promote<T, U>::V, Dim0, Dim1, Dim4, Dim3, i, j, m, l>
  operator*(const Tensor2_Expr<B, U, Dim2, Dim4, k, m> &b,
            const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a)
  {
    typedef const Tensor4_times_Tensor2_2_0<A, B, T, U, Dim0, Dim1, Dim2, Dim3,
                                            Dim4, i, j, k, l, m>
      TensorExpr;
    return Tensor4_Expr<TensorExpr, typename promote<T, U>::V, Dim0, Dim1,
                        Dim4, Dim3, i, j, m, l>(TensorExpr(a, b));
  }

  /* A(i,j,k,l)*B(m,k) -> Tensor4 */

  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2,
            int Dim3, int Dim4, char i, char j, char k, char l, char m>
  class Tensor4_times_Tensor2_2_1
  {
    const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> iterA;
    const Tensor2_Expr<B, U, Dim4, Dim2, m, k> iterB;

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
    Tensor4_times_Tensor2_2_1(
      const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a,
      const Tensor2_Expr<B, U, Dim4, Dim2, m, k> &b)
        : iterA(a), iterB(b)
    {}
    typename promote<T, U>::V
    operator()(const int N1, const int N2, const int N3, const int N4) const
    {
      return eval(N1, N2, N3, N4, Number<Dim3>());
    }
  };

  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2,
            int Dim3, int Dim4, char i, char j, char k, char l, char m>
  inline const Tensor4_Expr<
    const Tensor4_times_Tensor2_2_1<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4,
                                    i, j, k, l, m>,
    typename promote<T, U>::V, Dim0, Dim1, Dim4, Dim3, i, j, m, l>
  operator*(const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a,
            const Tensor2_Expr<B, U, Dim4, Dim2, m, k> &b)
  {
    typedef const Tensor4_times_Tensor2_2_1<A, B, T, U, Dim0, Dim1, Dim2, Dim3,
                                            Dim4, i, j, k, l, m>
      TensorExpr;
    return Tensor4_Expr<TensorExpr, typename promote<T, U>::V, Dim0, Dim1,
                        Dim4, Dim3, i, j, m, l>(TensorExpr(a, b));
  }

  /* B(m,k)*A(i,j,k,l) -> Tensor4 */

  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2,
            int Dim3, int Dim4, char i, char j, char k, char l, char m>
  inline const Tensor4_Expr<
    const Tensor4_times_Tensor2_2_1<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4,
                                    i, j, k, l, m>,
    typename promote<T, U>::V, Dim0, Dim1, Dim4, Dim3, i, j, m, l>
  operator*(const Tensor2_Expr<B, U, Dim4, Dim2, m, k> &b,
            const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a)
  {
    typedef const Tensor4_times_Tensor2_2_1<A, B, T, U, Dim0, Dim1, Dim2, Dim3,
                                            Dim4, i, j, k, l, m>
      TensorExpr;
    return Tensor4_Expr<TensorExpr, typename promote<T, U>::V, Dim0, Dim1,
                        Dim4, Dim3, i, j, m, l>(TensorExpr(a, b));
  }

} // namespace FTensor