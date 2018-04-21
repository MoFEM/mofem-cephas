/* This file has all of the declarations for expressions like
   Tensor3*Tensor2 and Tensor2*Tensor3, yielding a
   Tensor3 or Tensor1. */

#pragma once

namespace FTensor
{
  /* A(i,j,k)*B(j,k)->Tensor1 */

  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2,
            char i, char j, char k>
  class Tensor3_times_Tensor2_12
  {
    const Tensor3_Expr<A, T, Dim0, Dim1, Dim2, i, j, k> iterA;
    const Tensor2_Expr<B, U, Dim1, Dim2, j, k> iterB;

    template <int Current_Dim0, int Current_Dim1>
    typename promote<T, U>::V eval(const int N1, const Number<Current_Dim0> &,
                                   const Number<Current_Dim1> &) const
    {
      return iterA(N1, Current_Dim0 - 1, Current_Dim1 - 1)
               * iterB(Current_Dim0 - 1, Current_Dim1 - 1)
             + eval(N1, Number<Current_Dim0 - 1>(), Number<Current_Dim1>());
    }
    template <int Current_Dim1>
    typename promote<T, U>::V
    eval(const int N1, const Number<1> &, const Number<Current_Dim1> &) const
    {
      return iterA(N1, 0, Current_Dim1 - 1) * iterB(0, Current_Dim1 - 1)
             + eval(N1, Number<Dim1>(), Number<Current_Dim1 - 1>());
    }
    typename promote<T, U>::V
    eval(const int N1, const Number<1> &, const Number<1> &) const
    {
      return iterA(N1, 0, 0) * iterB(0, 0);
    }

  public:
    Tensor3_times_Tensor2_12(
      const Tensor3_Expr<A, T, Dim0, Dim1, Dim2, i, j, k> &a,
      const Tensor2_Expr<B, U, Dim1, Dim2, j, k> &b)
        : iterA(a), iterB(b)
    {}
    typename promote<T, U>::V operator()(const int N1) const
    {
      return eval(N1, Number<Dim1>(), Number<Dim2>());
    }
  };

  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2,
            char i, char j, char k>
  inline const Tensor1_Expr<
    const Tensor3_times_Tensor2_12<A, B, T, U, Dim0, Dim1, Dim2, i, j, k>,
    typename promote<T, U>::V, Dim0, i>
  operator*(const Tensor3_Expr<A, T, Dim0, Dim1, Dim2, i, j, k> &a,
            const Tensor2_Expr<B, U, Dim1, Dim2, j, k> &b)
  {
    typedef const Tensor3_times_Tensor2_12<A, B, T, U, Dim0, Dim1, Dim2, i, j,
                                           k>
      TensorExpr;
    return Tensor1_Expr<TensorExpr, typename promote<T, U>::V, Dim0, i>(
      TensorExpr(a, b));
  }

  /* B(j,k)*A(i,j,k)->Tensor1 */

  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2,
            char i, char j, char k>
  inline const Tensor1_Expr<
    const Tensor3_times_Tensor2_12<A, B, T, U, Dim0, Dim1, Dim2, i, j, k>,
    typename promote<T, U>::V, Dim0, i>
  operator*(const Tensor2_Expr<B, U, Dim1, Dim2, j, k> &b,
            const Tensor3_Expr<A, T, Dim0, Dim1, Dim2, i, j, k> &a)
  {
    typedef const Tensor3_times_Tensor2_12<A, B, T, U, Dim0, Dim1, Dim2, i, j,
                                           k>
      TensorExpr;
    return Tensor1_Expr<TensorExpr, typename promote<T, U>::V, Dim0, i>(
      TensorExpr(a, b));
  }

  /* A(j,i,k)*B(j,k)->Tensor1 */

  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2,
            char i, char j, char k>
  class Tensor3_times_Tensor2_02
  {
    const Tensor3_Expr<A, T, Dim0, Dim1, Dim2, j, i, k> iterA;
    const Tensor2_Expr<B, U, Dim0, Dim2, j, k> iterB;

    template <int Current_Dim0, int Current_Dim1>
    typename promote<T, U>::V eval(const int N1, const Number<Current_Dim0> &,
                                   const Number<Current_Dim1> &) const
    {
      return iterA(Current_Dim0 - 1, N1, Current_Dim1 - 1)
               * iterB(Current_Dim0 - 1, Current_Dim1 - 1)
             + eval(N1, Number<Current_Dim0 - 1>(), Number<Current_Dim1>());
    }
    template <int Current_Dim1>
    typename promote<T, U>::V
    eval(const int N1, const Number<1> &, const Number<Current_Dim1> &) const
    {
      return iterA(0, N1, Current_Dim1 - 1) * iterB(0, Current_Dim1 - 1)
             + eval(N1, Number<Dim0>(), Number<Current_Dim1 - 1>());
    }
    typename promote<T, U>::V
    eval(const int N1, const Number<1> &, const Number<1> &) const
    {
      return iterA(0, N1, 0) * iterB(0, 0);
    }

  public:
    Tensor3_times_Tensor2_02(
      const Tensor3_Expr<A, T, Dim0, Dim1, Dim2, j, i, k> &a,
      const Tensor2_Expr<B, U, Dim0, Dim2, j, k> &b)
        : iterA(a), iterB(b)
    {}
    typename promote<T, U>::V operator()(const int N1) const
    {
      return eval(N1, Number<Dim0>(), Number<Dim2>());
    }
  };

  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2,
            char i, char j, char k>
  inline const Tensor1_Expr<
    const Tensor3_times_Tensor2_02<A, B, T, U, Dim0, Dim1, Dim2, i, j, k>,
    typename promote<T, U>::V, Dim1, i>
  operator*(const Tensor3_Expr<A, T, Dim0, Dim1, Dim2, j, i, k> &a,
            const Tensor2_Expr<B, U, Dim0, Dim2, j, k> &b)
  {
    typedef const Tensor3_times_Tensor2_02<A, B, T, U, Dim0, Dim1, Dim2, i, j,
                                           k>
      TensorExpr;
    return Tensor1_Expr<TensorExpr, typename promote<T, U>::V, Dim1, i>(
      TensorExpr(a, b));
  }

  /* B(j,k)*A(j,i,k)->Tensor1 */

  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2,
            char i, char j, char k>
  inline const Tensor1_Expr<
    const Tensor3_times_Tensor2_02<A, B, T, U, Dim0, Dim1, Dim2, i, j, k>,
    typename promote<T, U>::V, Dim1, i>
  operator*(const Tensor2_Expr<B, U, Dim0, Dim2, j, k> &b,
            const Tensor3_Expr<A, T, Dim0, Dim1, Dim2, j, i, k> &a)
  {
    typedef const Tensor3_times_Tensor2_02<A, B, T, U, Dim0, Dim1, Dim2, i, j,
                                           k>
      TensorExpr;
    return Tensor1_Expr<TensorExpr, typename promote<T, U>::V, Dim1, i>(
      TensorExpr(a, b));
  }

  /* A(j,k,i)*B(j,k)->Tensor1 */

  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2,
            char i, char j, char k>
  class Tensor3_times_Tensor2_01
  {
    const Tensor3_Expr<A, T, Dim0, Dim1, Dim2, j, k, i> iterA;
    const Tensor2_Expr<B, U, Dim0, Dim1, j, k> iterB;

    template <int Current_Dim0, int Current_Dim1>
    typename promote<T, U>::V eval(const int N1, const Number<Current_Dim0> &,
                                   const Number<Current_Dim1> &) const
    {
      return iterA(Current_Dim0 - 1, Current_Dim1 - 1, N1)
               * iterB(Current_Dim0 - 1, Current_Dim1 - 1)
             + eval(N1, Number<Current_Dim0 - 1>(), Number<Current_Dim1>());
    }
    template <int Current_Dim1>
    typename promote<T, U>::V
    eval(const int N1, const Number<1> &, const Number<Current_Dim1> &) const
    {
      return iterA(0, Current_Dim1 - 1, N1) * iterB(0, Current_Dim1 - 1)
             + eval(N1, Number<Dim0>(), Number<Current_Dim1 - 1>());
    }
    typename promote<T, U>::V
    eval(const int N1, const Number<1> &, const Number<1> &) const
    {
      return iterA(0, 0, N1) * iterB(0, 0);
    }

  public:
    Tensor3_times_Tensor2_01(
      const Tensor3_Expr<A, T, Dim0, Dim1, Dim2, j, k, i> &a,
      const Tensor2_Expr<B, U, Dim0, Dim1, j, k> &b)
        : iterA(a), iterB(b)
    {}
    typename promote<T, U>::V operator()(const int N1) const
    {
      return eval(N1, Number<Dim0>(), Number<Dim1>());
    }
  };

  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2,
            char i, char j, char k>
  inline const Tensor1_Expr<
    const Tensor3_times_Tensor2_01<A, B, T, U, Dim0, Dim1, Dim2, i, j, k>,
    typename promote<T, U>::V, Dim2, i>
  operator*(const Tensor3_Expr<A, T, Dim0, Dim1, Dim2, j, k, i> &a,
            const Tensor2_Expr<B, U, Dim0, Dim1, j, k> &b)
  {
    typedef const Tensor3_times_Tensor2_01<A, B, T, U, Dim0, Dim1, Dim2, i, j,
                                           k>
      TensorExpr;
    return Tensor1_Expr<TensorExpr, typename promote<T, U>::V, Dim2, i>(
      TensorExpr(a, b));
  }

  /* B(j,k)*A(j,k,i)->Tensor1 */

  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2,
            char i, char j, char k>
  inline const Tensor1_Expr<
    const Tensor3_times_Tensor2_01<A, B, T, U, Dim0, Dim1, Dim2, i, j, k>,
    typename promote<T, U>::V, Dim2, i>
  operator*(const Tensor2_Expr<B, U, Dim0, Dim1, j, k> &b,
            const Tensor3_Expr<A, T, Dim0, Dim1, Dim2, j, k, i> &a)
  {
    typedef const Tensor3_times_Tensor2_01<A, B, T, U, Dim0, Dim1, Dim2, i, j,
                                           k>
      TensorExpr;
    return Tensor1_Expr<TensorExpr, typename promote<T, U>::V, Dim2, i>(
      TensorExpr(a, b));
  }

  /* A(i,j,k)*B(k,l)->Tensor3 */

  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2,
            int Dim3, char i, char j, char k, char l>
  class Tensor3_times_Tensor2_2_01
  {
    const Tensor3_Expr<A, T, Dim0, Dim1, Dim2, i, j, k> iterA;
    const Tensor2_Expr<B, U, Dim2, Dim3, k, l> iterB;

    template <int Current_Dim>
    typename promote<T, U>::V eval(const int N1, const int N2, const int N3,
                                   const Number<Current_Dim> &) const
    {
      return iterA(N1, N2, Current_Dim - 1) * iterB(Current_Dim - 1, N3)
             + eval(N1, N2, N3, Number<Current_Dim - 1>());
    }
    typename promote<T, U>::V
    eval(const int N1, const int N2, const int N3, const Number<1> &) const
    {
      return iterA(N1, N2, 0) * iterB(0, N3);
    }

  public:
    Tensor3_times_Tensor2_2_01(
      const Tensor3_Expr<A, T, Dim0, Dim1, Dim2, i, j, k> &a,
      const Tensor2_Expr<B, U, Dim2, Dim3, k, l> &b)
        : iterA(a), iterB(b)
    {}
    typename promote<T, U>::V
    operator()(const int N1, const int N2, const int N3) const
    {
      return eval(N1, N2, N3, Number<Dim2>());
    }
  };

  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2,
            int Dim3, char i, char j, char k, char l>
  inline const Tensor3_Expr<const Tensor3_times_Tensor2_2_01<
                              A, B, T, U, Dim0, Dim1, Dim2, Dim3, i, j, k, l>,
                            typename promote<T, U>::V, Dim0, Dim1, Dim3, i, j,
                            l>
  operator*(const Tensor3_Expr<A, T, Dim0, Dim1, Dim2, i, j, k> &a,
            const Tensor2_Expr<B, U, Dim2, Dim3, k, l> &b)
  {
    typedef const Tensor3_times_Tensor2_2_01<A, B, T, U, Dim0, Dim1, Dim2,
                                             Dim3, i, j, k, l>
      TensorExpr;
    return Tensor3_Expr<TensorExpr, typename promote<T, U>::V, Dim0, Dim1,
                        Dim3, i, j, l>(TensorExpr(a, b));
  }

  /* B(k,l)*A(i,j,k)->Tensor3 */

  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2,
            int Dim3, char i, char j, char k, char l>
  inline const Tensor3_Expr<const Tensor3_times_Tensor2_2_01<
                              A, B, T, U, Dim0, Dim1, Dim2, Dim3, i, j, k, l>,
                            typename promote<T, U>::V, Dim0, Dim1, Dim3, i, j,
                            l>
  operator*(const Tensor2_Expr<B, U, Dim2, Dim3, k, l> &b,
            const Tensor3_Expr<A, T, Dim0, Dim1, Dim2, i, j, k> &a)
  {
    typedef const Tensor3_times_Tensor2_2_01<A, B, T, U, Dim0, Dim1, Dim2,
                                             Dim3, i, j, k, l>
      TensorExpr;
    return Tensor3_Expr<TensorExpr, typename promote<T, U>::V, Dim0, Dim1,
                        Dim3, i, j, l>(TensorExpr(a, b));
  }

  /* A(i,j,k)*B(l,k)->Tensor3 */

  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2,
            int Dim3, char i, char j, char k, char l>
  class Tensor3_times_Tensor2_2_10
  {
    const Tensor3_Expr<A, T, Dim0, Dim1, Dim2, i, j, k> iterA;
    const Tensor2_Expr<B, U, Dim3, Dim2, l, k> iterB;

    template <int Current_Dim>
    typename promote<T, U>::V eval(const int N1, const int N2, const int N3,
                                   const Number<Current_Dim> &) const
    {
      return iterA(N1, N2, Current_Dim - 1) * iterB(N3, Current_Dim - 1)
             + eval(N1, N2, N3, Number<Current_Dim - 1>());
    }
    typename promote<T, U>::V
    eval(const int N1, const int N2, const int N3, const Number<1> &) const
    {
      return iterA(N1, N2, 0) * iterB(N3, 0);
    }

  public:
    Tensor3_times_Tensor2_2_10(
      const Tensor3_Expr<A, T, Dim0, Dim1, Dim2, i, j, k> &a,
      const Tensor2_Expr<B, U, Dim2, Dim3, l, k> &b)
        : iterA(a), iterB(b)
    {}
    typename promote<T, U>::V
    operator()(const int N1, const int N2, const int N3) const
    {
      return eval(N1, N2, N3, Number<Dim2>());
    }
  };

  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2,
            int Dim3, char i, char j, char k, char l>
  inline const Tensor3_Expr<const Tensor3_times_Tensor2_2_10<
                              A, B, T, U, Dim0, Dim1, Dim2, Dim3, i, j, k, l>,
                            typename promote<T, U>::V, Dim0, Dim1, Dim3, i, j,
                            l>
  operator*(const Tensor3_Expr<A, T, Dim0, Dim1, Dim2, i, j, k> &a,
            const Tensor2_Expr<B, U, Dim3, Dim2, l, k> &b)
  {
    typedef const Tensor3_times_Tensor2_2_10<A, B, T, U, Dim0, Dim1, Dim2,
                                             Dim3, i, j, k, l>
      TensorExpr;
    return Tensor3_Expr<TensorExpr, typename promote<T, U>::V, Dim0, Dim1,
                        Dim3, i, j, l>(TensorExpr(a, b));
  }

  /* B(l,k)*A(i,j,k)->Tensor3 */

  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2,
            int Dim3, char i, char j, char k, char l>
  inline const Tensor3_Expr<const Tensor3_times_Tensor2_2_10<
                              A, B, T, U, Dim0, Dim1, Dim2, Dim3, i, j, k, l>,
                            typename promote<T, U>::V, Dim0, Dim1, Dim3, i, j,
                            l>
  operator*(const Tensor2_Expr<B, U, Dim3, Dim2, l, k> &b,
            const Tensor3_Expr<A, T, Dim0, Dim1, Dim2, i, j, k> &a)
  {
    typedef const Tensor3_times_Tensor2_2_10<A, B, T, U, Dim0, Dim1, Dim2,
                                             Dim3, i, j, k, l>
      TensorExpr;
    return Tensor3_Expr<TensorExpr, typename promote<T, U>::V, Dim0, Dim1,
                        Dim3, i, j, l>(TensorExpr(a, b));
  }

  /* A(i,j,l)*B(j,k)->Tensor3 */

  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2,
            int Dim3, char i, char j, char k, char l>
  class Tensor3_times_Tensor2_1_01
  {
    const Tensor3_Expr<A, T, Dim0, Dim1, Dim3, i, j, l> iterA;
    const Tensor2_Expr<B, U, Dim1, Dim2, j, k> iterB;

    template <int Current_Dim>
    typename promote<T, U>::V eval(const int N1, const int N2, const int N3,
                                   const Number<Current_Dim> &) const
    {
      return iterA(N1, Current_Dim - 1, N3) * iterB(Current_Dim - 1, N2)
             + eval(N1, N2, N3, Number<Current_Dim - 1>());
    }
    typename promote<T, U>::V
    eval(const int N1, const int N2, const int N3, const Number<1> &) const
    {
      return iterA(N1, 0, N3) * iterB(0, N2);
    }

  public:
    Tensor3_times_Tensor2_1_01(
      const Tensor3_Expr<A, T, Dim0, Dim1, Dim3, i, j, l> &a,
      const Tensor2_Expr<B, U, Dim1, Dim2, j, k> &b)
        : iterA(a), iterB(b)
    {}
    typename promote<T, U>::V
    operator()(const int N1, const int N2, const int N3) const
    {
      return eval(N1, N2, N3, Number<Dim1>());
    }
  };

  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2,
            int Dim3, char i, char j, char k, char l>
  inline const Tensor3_Expr<const Tensor3_times_Tensor2_1_01<
                              A, B, T, U, Dim0, Dim1, Dim2, Dim3, i, j, k, l>,
                            typename promote<T, U>::V, Dim0, Dim2, Dim3, i, k,
                            l>
  operator*(const Tensor3_Expr<A, T, Dim0, Dim1, Dim3, i, j, l> &a,
            const Tensor2_Expr<B, U, Dim1, Dim2, j, k> &b)
  {
    typedef const Tensor3_times_Tensor2_1_01<A, B, T, U, Dim0, Dim1, Dim2,
                                             Dim3, i, j, k, l>
      TensorExpr;
    return Tensor3_Expr<TensorExpr, typename promote<T, U>::V, Dim0, Dim2,
                        Dim3, i, k, l>(TensorExpr(a, b));
  }

  /* B(j,k)*A(i,j,l)->Tensor3 */

  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2,
            int Dim3, char i, char j, char k, char l>
  inline const Tensor3_Expr<const Tensor3_times_Tensor2_1_01<
                              A, B, T, U, Dim0, Dim1, Dim2, Dim3, i, j, k, l>,
                            typename promote<T, U>::V, Dim0, Dim2, Dim3, i, k,
                            l>
  operator*(const Tensor2_Expr<B, U, Dim1, Dim2, j, k> &b,
            const Tensor3_Expr<A, T, Dim0, Dim1, Dim3, i, j, l> &a)
  {
    typedef const Tensor3_times_Tensor2_1_01<A, B, T, U, Dim0, Dim1, Dim2,
                                             Dim3, i, j, k, l>
      TensorExpr;
    return Tensor3_Expr<TensorExpr, typename promote<T, U>::V, Dim0, Dim2,
                        Dim3, i, k, l>(TensorExpr(a, b));
  }

  /* A(i,j,l)*B(k,j)->Tensor3 */

  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2,
            int Dim3, char i, char j, char k, char l>
  class Tensor3_times_Tensor2_1_10
  {
    const Tensor3_Expr<A, T, Dim0, Dim1, Dim3, i, j, l> iterA;
    const Tensor2_Expr<B, U, Dim2, Dim1, k, j> iterB;

    template <int Current_Dim>
    typename promote<T, U>::V eval(const int N1, const int N2, const int N3,
                                   const Number<Current_Dim> &) const
    {
      return iterA(N1, Current_Dim - 1, N3) * iterB(N2, Current_Dim - 1)
             + eval(N1, N2, N3, Number<Current_Dim - 1>());
    }
    typename promote<T, U>::V
    eval(const int N1, const int N2, const int N3, const Number<1> &) const
    {
      return iterA(N1, 0, N3) * iterB(N2, 0);
    }

  public:
    Tensor3_times_Tensor2_1_10(
      const Tensor3_Expr<A, T, Dim0, Dim1, Dim3, i, j, l> &a,
      const Tensor2_Expr<B, U, Dim2, Dim1, k, j> &b)
        : iterA(a), iterB(b)
    {}
    typename promote<T, U>::V
    operator()(const int N1, const int N2, const int N3) const
    {
      return eval(N1, N2, N3, Number<Dim1>());
    }
  };

  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2,
            int Dim3, char i, char j, char k, char l>
  inline const Tensor3_Expr<const Tensor3_times_Tensor2_1_10<
                              A, B, T, U, Dim0, Dim1, Dim2, Dim3, i, j, k, l>,
                            typename promote<T, U>::V, Dim0, Dim2, Dim3, i, k,
                            l>
  operator*(const Tensor3_Expr<A, T, Dim0, Dim1, Dim3, i, j, l> &a,
            const Tensor2_Expr<B, U, Dim2, Dim1, k, j> &b)
  {
    typedef const Tensor3_times_Tensor2_1_10<A, B, T, U, Dim0, Dim1, Dim2,
                                             Dim3, i, j, k, l>
      TensorExpr;
    return Tensor3_Expr<TensorExpr, typename promote<T, U>::V, Dim0, Dim2,
                        Dim3, i, k, l>(TensorExpr(a, b));
  }

  /* B(k,j)*A(i,j,l)->Tensor3 */

  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2,
            int Dim3, char i, char j, char k, char l>
  inline const Tensor3_Expr<const Tensor3_times_Tensor2_1_10<
                              A, B, T, U, Dim0, Dim1, Dim2, Dim3, i, j, k, l>,
                            typename promote<T, U>::V, Dim0, Dim2, Dim3, i, k,
                            l>
  operator*(const Tensor2_Expr<B, U, Dim2, Dim1, k, j> &b,
            const Tensor3_Expr<A, T, Dim0, Dim1, Dim3, i, j, l> &a)
  {
    typedef const Tensor3_times_Tensor2_1_10<A, B, T, U, Dim0, Dim1, Dim2,
                                             Dim3, i, j, k, l>
      TensorExpr;
    return Tensor3_Expr<TensorExpr, typename promote<T, U>::V, Dim0, Dim2,
                        Dim3, i, k, l>(TensorExpr(a, b));
  }

  /* A(i,k,l)*B(i,j)->Tensor3 */

  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2,
            int Dim3, char i, char j, char k, char l>
  class Tensor3_times_Tensor2_0_01
  {
    const Tensor3_Expr<A, T, Dim0, Dim2, Dim3, i, k, l> iterA;
    const Tensor2_Expr<B, U, Dim0, Dim1, i, j> iterB;

    template <int Current_Dim>
    typename promote<T, U>::V eval(const int N1, const int N2, const int N3,
                                   const Number<Current_Dim> &) const
    {
      return iterA(Current_Dim - 1, N2, N3) * iterB(Current_Dim - 1, N1)
             + eval(N1, N2, N3, Number<Current_Dim - 1>());
    }
    typename promote<T, U>::V
    eval(const int N1, const int N2, const int N3, const Number<1> &) const
    {
      return iterA(0, N2, N3) * iterB(0, N1);
    }

  public:
    Tensor3_times_Tensor2_0_01(
      const Tensor3_Expr<A, T, Dim0, Dim2, Dim3, i, k, l> &a,
      const Tensor2_Expr<B, U, Dim0, Dim1, i, j> &b)
        : iterA(a), iterB(b)
    {}
    typename promote<T, U>::V
    operator()(const int N1, const int N2, const int N3) const
    {
      return eval(N1, N2, N3, Number<Dim0>());
    }
  };

  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2,
            int Dim3, char i, char j, char k, char l>
  inline const Tensor3_Expr<const Tensor3_times_Tensor2_0_01<
                              A, B, T, U, Dim0, Dim1, Dim2, Dim3, i, j, k, l>,
                            typename promote<T, U>::V, Dim1, Dim2, Dim3, j, k,
                            l>
  operator*(const Tensor3_Expr<A, T, Dim0, Dim2, Dim3, i, k, l> &a,
            const Tensor2_Expr<B, U, Dim0, Dim1, i, j> &b)
  {
    typedef const Tensor3_times_Tensor2_0_01<A, B, T, U, Dim0, Dim1, Dim2,
                                             Dim3, i, j, k, l>
      TensorExpr;
    return Tensor3_Expr<TensorExpr, typename promote<T, U>::V, Dim1, Dim2,
                        Dim3, j, k, l>(TensorExpr(a, b));
  }

  /* B(i,j)*A(i,k,l)->Tensor3 */

  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2,
            int Dim3, char i, char j, char k, char l>
  inline const Tensor3_Expr<const Tensor3_times_Tensor2_0_01<
                              A, B, T, U, Dim0, Dim1, Dim2, Dim3, i, j, k, l>,
                            typename promote<T, U>::V, Dim1, Dim2, Dim3, j, k,
                            l>
  operator*(const Tensor2_Expr<B, U, Dim0, Dim1, i, j> &b,
            const Tensor3_Expr<A, T, Dim0, Dim2, Dim3, i, k, l> &a)
  {
    typedef const Tensor3_times_Tensor2_0_01<A, B, T, U, Dim0, Dim1, Dim2,
                                             Dim3, i, j, k, l>
      TensorExpr;
    return Tensor3_Expr<TensorExpr, typename promote<T, U>::V, Dim1, Dim2,
                        Dim3, j, k, l>(TensorExpr(a, b));
  }

  /* A(i,k,l)*B(j,i)->Tensor3 */

  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2,
            int Dim3, char i, char j, char k, char l>
  class Tensor3_times_Tensor2_0_10
  {
    const Tensor3_Expr<A, T, Dim0, Dim2, Dim3, i, k, l> iterA;
    const Tensor2_Expr<B, U, Dim1, Dim0, j, i> iterB;

    template <int Current_Dim>
    typename promote<T, U>::V eval(const int N1, const int N2, const int N3,
                                   const Number<Current_Dim> &) const
    {
      return iterA(Current_Dim - 1, N2, N3) * iterB(N1, Current_Dim - 1)
             + eval(N1, N2, N3, Number<Current_Dim - 1>());
    }
    typename promote<T, U>::V
    eval(const int N1, const int N2, const int N3, const Number<1> &) const
    {
      return iterA(0, N2, N3) * iterB(N1, 0);
    }

  public:
    Tensor3_times_Tensor2_0_10(
      const Tensor3_Expr<A, T, Dim0, Dim2, Dim3, i, k, l> &a,
      const Tensor2_Expr<B, U, Dim1, Dim0, j, i> &b)
        : iterA(a), iterB(b)
    {}
    typename promote<T, U>::V
    operator()(const int N1, const int N2, const int N3) const
    {
      return eval(N1, N2, N3, Number<Dim0>());
    }
  };

  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2,
            int Dim3, char i, char j, char k, char l>
  inline const Tensor3_Expr<const Tensor3_times_Tensor2_0_10<
                              A, B, T, U, Dim0, Dim1, Dim2, Dim3, i, j, k, l>,
                            typename promote<T, U>::V, Dim1, Dim2, Dim3, j, k,
                            l>
  operator*(const Tensor3_Expr<A, T, Dim0, Dim2, Dim3, i, k, l> &a,
            const Tensor2_Expr<B, U, Dim1, Dim0, j, i> &b)
  {
    typedef const Tensor3_times_Tensor2_0_10<A, B, T, U, Dim0, Dim1, Dim2,
                                             Dim3, i, j, k, l>
      TensorExpr;
    return Tensor3_Expr<TensorExpr, typename promote<T, U>::V, Dim1, Dim2,
                        Dim3, j, k, l>(TensorExpr(a, b));
  }

  /* B(j,i)*A(i,k,l)->Tensor3 */

  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2,
            int Dim3, char i, char j, char k, char l>
  inline const Tensor3_Expr<const Tensor3_times_Tensor2_0_10<
                              A, B, T, U, Dim0, Dim1, Dim2, Dim3, i, j, k, l>,
                            typename promote<T, U>::V, Dim1, Dim2, Dim3, j, k,
                            l>
  operator*(const Tensor2_Expr<B, U, Dim1, Dim0, j, i> &b,
            const Tensor3_Expr<A, T, Dim0, Dim2, Dim3, i, k, l> &a)
  {
    typedef const Tensor3_times_Tensor2_0_10<A, B, T, U, Dim0, Dim1, Dim2,
                                             Dim3, i, j, k, l>
      TensorExpr;
    return Tensor3_Expr<TensorExpr, typename promote<T, U>::V, Dim1, Dim2,
                        Dim3, j, k, l>(TensorExpr(a, b));
  }
}
