/* Fully contracts a Tensor3 with a Tensor3, yielding a typename
   promote<T,U>::V. */

#pragma once

namespace FTensor
{
  /* A(i,j,k)*B(i,j,k) */

  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2,
            char i, char j, char k, int Current_Dim0, int Current_Dim1,
            int Current_Dim2>
  typename promote<T, U>::V
  T3_times_T3_012(const Tensor3_Expr<A, T, Dim0, Dim1, Dim2, i, j, k> &a,
                  const Tensor3_Expr<B, U, Dim0, Dim1, Dim2, i, j, k> &b,
                  const Number<Current_Dim0> &, const Number<Current_Dim1> &,
                  const Number<Current_Dim2> &)
  {
    return a(Current_Dim0 - 1, Current_Dim1 - 1, Current_Dim2 - 1)
             * b(Current_Dim0 - 1, Current_Dim1 - 1, Current_Dim2 - 1)
           + T3_times_T3_012(a, b, Number<Current_Dim0>(),
                             Number<Current_Dim1 - 1>(),
                             Number<Current_Dim2>());
  }

  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2,
            char i, char j, char k, int Current_Dim0, int Current_Dim2>
  typename promote<T, U>::V
  T3_times_T3_012(const Tensor3_Expr<A, T, Dim0, Dim1, Dim2, i, j, k> &a,
                  const Tensor3_Expr<B, U, Dim0, Dim1, Dim2, i, j, k> &b,
                  const Number<Current_Dim0> &, const Number<1> &,
                  const Number<Current_Dim2> &)
  {
    return a(Current_Dim0 - 1, 0, Current_Dim2 - 1)
             * b(Current_Dim0 - 1, 0, Current_Dim2 - 1)
           + T3_times_T3_012(a, b, Number<Current_Dim0>(), Number<Dim1>(),
                             Number<Current_Dim2 - 1>());
  }

  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2,
            char i, char j, char k, int Current_Dim0>
  typename promote<T, U>::V
  T3_times_T3_012(const Tensor3_Expr<A, T, Dim0, Dim1, Dim2, i, j, k> &a,
                  const Tensor3_Expr<B, U, Dim0, Dim1, Dim2, i, j, k> &b,
                  const Number<Current_Dim0> &, const Number<1> &,
                  const Number<1> &)
  {
    return a(Current_Dim0 - 1, 0, 0) * b(Current_Dim0 - 1, 0, 0)
           + T3_times_T3_012(a, b, Number<Current_Dim0 - 1>(), Number<Dim1>(),
                             Number<Dim2>());
  }

  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2,
            char i, char j, char k>
  typename promote<T, U>::V
  T3_times_T3_012(const Tensor3_Expr<A, T, Dim0, Dim1, Dim2, i, j, k> &a,
                  const Tensor3_Expr<B, U, Dim0, Dim1, Dim2, i, j, k> &b,
                  const Number<1> &, const Number<1> &, const Number<1> &)
  {
    return a(0, 0, 0) * b(0, 0, 0);
  }

  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2,
            char i, char j, char k>
  typename promote<T, U>::V
  operator*(const Tensor3_Expr<A, T, Dim0, Dim1, Dim2, i, j, k> &a,
            const Tensor3_Expr<B, U, Dim0, Dim1, Dim2, i, j, k> &b)
  {
    return T3_times_T3_012(a, b, Number<Dim0>(), Number<Dim1>(),
                           Number<Dim2>());
  }

  /* A(i,j,k)*B(k,i,j) */

  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2,
            char i, char j, char k, int Current_Dim0, int Current_Dim1,
            int Current_Dim2>
  typename promote<T, U>::V
  T3_times_T3_201(const Tensor3_Expr<A, T, Dim0, Dim1, Dim2, i, j, k> &a,
                  const Tensor3_Expr<B, U, Dim2, Dim0, Dim1, k, i, j> &b,
                  const Number<Current_Dim0> &, const Number<Current_Dim1> &,
                  const Number<Current_Dim2> &)
  {
    return a(Current_Dim0 - 1, Current_Dim1 - 1, Current_Dim2 - 1)
             * b(Current_Dim2 - 1, Current_Dim0 - 1, Current_Dim1 - 1)
           + T3_times_T3_201(a, b, Number<Current_Dim0>(),
                             Number<Current_Dim1 - 1>(),
                             Number<Current_Dim2>());
  }

  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2,
            char i, char j, char k, int Current_Dim0, int Current_Dim2>
  typename promote<T, U>::V
  T3_times_T3_201(const Tensor3_Expr<A, T, Dim0, Dim1, Dim2, i, j, k> &a,
                  const Tensor3_Expr<B, U, Dim2, Dim0, Dim1, k, i, j> &b,
                  const Number<Current_Dim0> &, const Number<1> &,
                  const Number<Current_Dim2> &)
  {
    return a(Current_Dim0 - 1, 0, Current_Dim2 - 1)
             * b(Current_Dim2 - 1, Current_Dim0 - 1, 0)
           + T3_times_T3_201(a, b, Number<Current_Dim0>(), Number<Dim1>(),
                             Number<Current_Dim2 - 1>());
  }

  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2,
            char i, char j, char k, int Current_Dim0>
  typename promote<T, U>::V
  T3_times_T3_201(const Tensor3_Expr<A, T, Dim0, Dim1, Dim2, i, j, k> &a,
                  const Tensor3_Expr<B, U, Dim2, Dim0, Dim1, k, i, j> &b,
                  const Number<Current_Dim0> &, const Number<1> &,
                  const Number<1> &)
  {
    return a(Current_Dim0 - 1, 0, 0) * b(0, Current_Dim0 - 1, 0)
           + T3_times_T3_201(a, b, Number<Current_Dim0 - 1>(), Number<Dim1>(),
                             Number<Dim2>());
  }

  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2,
            char i, char j, char k>
  typename promote<T, U>::V
  T3_times_T3_201(const Tensor3_Expr<A, T, Dim0, Dim1, Dim2, i, j, k> &a,
                  const Tensor3_Expr<B, U, Dim2, Dim0, Dim1, k, i, j> &b,
                  const Number<1> &, const Number<1> &, const Number<1> &)
  {
    return a(0, 0, 0) * b(0, 0, 0);
  }

  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2,
            char i, char j, char k>
  typename promote<T, U>::V
  operator*(const Tensor3_Expr<A, T, Dim0, Dim1, Dim2, i, j, k> &a,
            const Tensor3_Expr<B, U, Dim2, Dim0, Dim1, k, i, j> &b)
  {
    return T3_times_T3_201(a, b, Number<Dim0>(), Number<Dim1>(),
                           Number<Dim2>());
  }

  /* A(i,j,k)*B(j,k,i) */

  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2,
            char i, char j, char k, int Current_Dim0, int Current_Dim1,
            int Current_Dim2>
  typename promote<T, U>::V
  T3_times_T3_120(const Tensor3_Expr<A, T, Dim0, Dim1, Dim2, i, j, k> &a,
                  const Tensor3_Expr<B, U, Dim1, Dim2, Dim0, j, k, i> &b,
                  const Number<Current_Dim0> &, const Number<Current_Dim1> &,
                  const Number<Current_Dim2> &)
  {
    return a(Current_Dim0 - 1, Current_Dim1 - 1, Current_Dim2 - 1)
             * b(Current_Dim1 - 1, Current_Dim2 - 1, Current_Dim0 - 1)
           + T3_times_T3_120(a, b, Number<Current_Dim0>(),
                             Number<Current_Dim1 - 1>(),
                             Number<Current_Dim2>());
  }

  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2,
            char i, char j, char k, int Current_Dim0, int Current_Dim2>
  typename promote<T, U>::V
  T3_times_T3_120(const Tensor3_Expr<A, T, Dim0, Dim1, Dim2, i, j, k> &a,
                  const Tensor3_Expr<B, U, Dim1, Dim2, Dim0, j, k, i> &b,
                  const Number<Current_Dim0> &, const Number<1> &,
                  const Number<Current_Dim2> &)
  {
    return a(Current_Dim0 - 1, 0, Current_Dim2 - 1)
             * b(0, Current_Dim2 - 1, Current_Dim0 - 1)
           + T3_times_T3_120(a, b, Number<Current_Dim0>(), Number<Dim1>(),
                             Number<Current_Dim2 - 1>());
  }

  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2,
            char i, char j, char k, int Current_Dim0>
  typename promote<T, U>::V
  T3_times_T3_120(const Tensor3_Expr<A, T, Dim0, Dim1, Dim2, i, j, k> &a,
                  const Tensor3_Expr<B, U, Dim1, Dim2, Dim0, j, k, i> &b,
                  const Number<Current_Dim0> &, const Number<1> &,
                  const Number<1> &)
  {
    return a(Current_Dim0 - 1, 0, 0) * b(0, 0, Current_Dim0 - 1)
           + T3_times_T3_120(a, b, Number<Current_Dim0 - 1>(), Number<Dim1>(),
                             Number<Dim2>());
  }

  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2,
            char i, char j, char k>
  typename promote<T, U>::V
  T3_times_T3_120(const Tensor3_Expr<A, T, Dim0, Dim1, Dim2, i, j, k> &a,
                  const Tensor3_Expr<B, U, Dim1, Dim2, Dim0, j, k, i> &b,
                  const Number<1> &, const Number<1> &, const Number<1> &)
  {
    return a(0, 0, 0) * b(0, 0, 0);
  }

  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2,
            char i, char j, char k>
  typename promote<T, U>::V
  operator*(const Tensor3_Expr<A, T, Dim0, Dim1, Dim2, i, j, k> &a,
            const Tensor3_Expr<B, U, Dim1, Dim2, Dim0, j, k, i> &b)
  {
    return T3_times_T3_120(a, b, Number<Dim0>(), Number<Dim1>(),
                           Number<Dim2>());
  }

  /* A(i,j,k)*B(j,i,k) */

  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2,
            char i, char j, char k, int Current_Dim0, int Current_Dim1,
            int Current_Dim2>
  typename promote<T, U>::V
  T3_times_T3_102(const Tensor3_Expr<A, T, Dim0, Dim1, Dim2, i, j, k> &a,
                  const Tensor3_Expr<B, U, Dim1, Dim0, Dim2, j, i, k> &b,
                  const Number<Current_Dim0> &, const Number<Current_Dim1> &,
                  const Number<Current_Dim2> &)
  {
    return a(Current_Dim0 - 1, Current_Dim1 - 1, Current_Dim2 - 1)
             * b(Current_Dim1 - 1, Current_Dim0 - 1, Current_Dim2 - 1)
           + T3_times_T3_102(a, b, Number<Current_Dim0>(),
                             Number<Current_Dim1 - 1>(),
                             Number<Current_Dim2>());
  }

  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2,
            char i, char j, char k, int Current_Dim0, int Current_Dim2>
  typename promote<T, U>::V
  T3_times_T3_102(const Tensor3_Expr<A, T, Dim0, Dim1, Dim2, i, j, k> &a,
                  const Tensor3_Expr<B, U, Dim1, Dim0, Dim2, j, i, k> &b,
                  const Number<Current_Dim0> &, const Number<1> &,
                  const Number<Current_Dim2> &)
  {
    return a(Current_Dim0 - 1, 0, Current_Dim2 - 1)
             * b(0, Current_Dim0 - 1, Current_Dim2 - 1)
           + T3_times_T3_102(a, b, Number<Current_Dim0>(), Number<Dim1>(),
                             Number<Current_Dim2 - 1>());
  }

  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2,
            char i, char j, char k, int Current_Dim0>
  typename promote<T, U>::V
  T3_times_T3_102(const Tensor3_Expr<A, T, Dim0, Dim1, Dim2, i, j, k> &a,
                  const Tensor3_Expr<B, U, Dim1, Dim0, Dim2, j, i, k> &b,
                  const Number<Current_Dim0> &, const Number<1> &,
                  const Number<1> &)
  {
    return a(Current_Dim0 - 1, 0, 0) * b(0, Current_Dim0 - 1, 0)
           + T3_times_T3_102(a, b, Number<Current_Dim0 - 1>(), Number<Dim1>(),
                             Number<Dim2>());
  }

  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2,
            char i, char j, char k>
  typename promote<T, U>::V
  T3_times_T3_102(const Tensor3_Expr<A, T, Dim0, Dim1, Dim2, i, j, k> &a,
                  const Tensor3_Expr<B, U, Dim1, Dim0, Dim2, j, i, k> &b,
                  const Number<1> &, const Number<1> &, const Number<1> &)
  {
    return a(0, 0, 0) * b(0, 0, 0);
  }

  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2,
            char i, char j, char k>
  typename promote<T, U>::V
  operator*(const Tensor3_Expr<A, T, Dim0, Dim1, Dim2, i, j, k> &a,
            const Tensor3_Expr<B, U, Dim1, Dim0, Dim2, j, i, k> &b)
  {
    return T3_times_T3_102(a, b, Number<Dim0>(), Number<Dim1>(),
                           Number<Dim2>());
  }

  /* A(i,j,k)*B(k,j,i) */

  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2,
            char i, char j, char k, int Current_Dim0, int Current_Dim1,
            int Current_Dim2>
  typename promote<T, U>::V
  T3_times_T3_210(const Tensor3_Expr<A, T, Dim0, Dim1, Dim2, i, j, k> &a,
                  const Tensor3_Expr<B, U, Dim2, Dim1, Dim0, k, j, i> &b,
                  const Number<Current_Dim0> &, const Number<Current_Dim1> &,
                  const Number<Current_Dim2> &)
  {
    return a(Current_Dim0 - 1, Current_Dim1 - 1, Current_Dim2 - 1)
             * b(Current_Dim2 - 1, Current_Dim1 - 1, Current_Dim0 - 1)
           + T3_times_T3_210(a, b, Number<Current_Dim0>(),
                             Number<Current_Dim1 - 1>(),
                             Number<Current_Dim2>());
  }

  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2,
            char i, char j, char k, int Current_Dim0, int Current_Dim2>
  typename promote<T, U>::V
  T3_times_T3_210(const Tensor3_Expr<A, T, Dim0, Dim1, Dim2, i, j, k> &a,
                  const Tensor3_Expr<B, U, Dim2, Dim1, Dim0, k, j, i> &b,
                  const Number<Current_Dim0> &, const Number<1> &,
                  const Number<Current_Dim2> &)
  {
    return a(Current_Dim0 - 1, 0, Current_Dim2 - 1)
             * b(Current_Dim2 - 1, 0, Current_Dim0 - 1)
           + T3_times_T3_210(a, b, Number<Current_Dim0>(), Number<Dim1>(),
                             Number<Current_Dim2 - 1>());
  }

  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2,
            char i, char j, char k, int Current_Dim0>
  typename promote<T, U>::V
  T3_times_T3_210(const Tensor3_Expr<A, T, Dim0, Dim1, Dim2, i, j, k> &a,
                  const Tensor3_Expr<B, U, Dim2, Dim1, Dim0, k, j, i> &b,
                  const Number<Current_Dim0> &, const Number<1> &,
                  const Number<1> &)
  {
    return a(Current_Dim0 - 1, 0, 0) * b(0, 0, Current_Dim0 - 1)
           + T3_times_T3_210(a, b, Number<Current_Dim0 - 1>(), Number<Dim1>(),
                             Number<Dim2>());
  }

  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2,
            char i, char j, char k>
  typename promote<T, U>::V
  T3_times_T3_210(const Tensor3_Expr<A, T, Dim0, Dim1, Dim2, i, j, k> &a,
                  const Tensor3_Expr<B, U, Dim2, Dim1, Dim0, k, j, i> &b,
                  const Number<1> &, const Number<1> &, const Number<1> &)
  {
    return a(0, 0, 0) * b(0, 0, 0);
  }

  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2,
            char i, char j, char k>
  typename promote<T, U>::V
  operator*(const Tensor3_Expr<A, T, Dim0, Dim1, Dim2, i, j, k> &a,
            const Tensor3_Expr<B, U, Dim2, Dim1, Dim0, k, j, i> &b)
  {
    return T3_times_T3_210(a, b, Number<Dim0>(), Number<Dim1>(),
                           Number<Dim2>());
  }

  /* A(i,j,k)*B(i,k,j) */

  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2,
            char i, char j, char k, int Current_Dim0, int Current_Dim1,
            int Current_Dim2>
  typename promote<T, U>::V
  T3_times_T3_021(const Tensor3_Expr<A, T, Dim0, Dim1, Dim2, i, j, k> &a,
                  const Tensor3_Expr<B, U, Dim0, Dim2, Dim1, i, k, j> &b,
                  const Number<Current_Dim0> &, const Number<Current_Dim1> &,
                  const Number<Current_Dim2> &)
  {
    return a(Current_Dim0 - 1, Current_Dim1 - 1, Current_Dim2 - 1)
             * b(Current_Dim0 - 1, Current_Dim2 - 1, Current_Dim1 - 1)
           + T3_times_T3_021(a, b, Number<Current_Dim0>(),
                             Number<Current_Dim1 - 1>(),
                             Number<Current_Dim2>());
  }

  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2,
            char i, char j, char k, int Current_Dim0, int Current_Dim2>
  typename promote<T, U>::V
  T3_times_T3_021(const Tensor3_Expr<A, T, Dim0, Dim1, Dim2, i, j, k> &a,
                  const Tensor3_Expr<B, U, Dim0, Dim2, Dim1, i, k, j> &b,
                  const Number<Current_Dim0> &, const Number<1> &,
                  const Number<Current_Dim2> &)
  {
    return a(Current_Dim0 - 1, 0, Current_Dim2 - 1)
             * b(Current_Dim0 - 1, Current_Dim2 - 1, 0)
           + T3_times_T3_021(a, b, Number<Current_Dim0>(), Number<Dim1>(),
                             Number<Current_Dim2 - 1>());
  }

  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2,
            char i, char j, char k, int Current_Dim0>
  typename promote<T, U>::V
  T3_times_T3_021(const Tensor3_Expr<A, T, Dim0, Dim1, Dim2, i, j, k> &a,
                  const Tensor3_Expr<B, U, Dim0, Dim2, Dim1, i, k, j> &b,
                  const Number<Current_Dim0> &, const Number<1> &,
                  const Number<1> &)
  {
    return a(Current_Dim0 - 1, 0, 0) * b(Current_Dim0 - 1, 0, 0)
           + T3_times_T3_021(a, b, Number<Current_Dim0 - 1>(), Number<Dim1>(),
                             Number<Dim2>());
  }

  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2,
            char i, char j, char k>
  typename promote<T, U>::V
  T3_times_T3_021(const Tensor3_Expr<A, T, Dim0, Dim1, Dim2, i, j, k> &a,
                  const Tensor3_Expr<B, U, Dim0, Dim2, Dim1, i, k, j> &b,
                  const Number<1> &, const Number<1> &, const Number<1> &)
  {
    return a(0, 0, 0) * b(0, 0, 0);
  }

  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2,
            char i, char j, char k>
  typename promote<T, U>::V
  operator*(const Tensor3_Expr<A, T, Dim0, Dim1, Dim2, i, j, k> &a,
            const Tensor3_Expr<B, U, Dim0, Dim2, Dim1, i, k, j> &b)
  {
    return T3_times_T3_021(a, b, Number<Dim0>(), Number<Dim1>(),
                           Number<Dim2>());
  }

  /* A(i,j,k)*B(j,i,l) -> Tensor2 */

  template <class A, class B, class T, class U, int Dim04, int Dim13, int Dim2,
            int Dim5, char i, char j, char k, char l>
  class Tensor3_times_Tensor3_12_21
  {
    Tensor3_Expr<A, T, Dim04, Dim13, Dim2, i, j, k> iterA;
    Tensor3_Expr<B, U, Dim13, Dim04, Dim5, j, i, l> iterB;

    template <int CurrentDim0, int CurrentDim1>
    typename promote<T, U>::V eval(const int N1, const int N2,
                                   const Number<CurrentDim0> &,
                                   const Number<CurrentDim1> &) const {
      return iterA(CurrentDim0 - 1, CurrentDim1 - 1, N1) *
                 iterB(CurrentDim1 - 1, CurrentDim0 - 1, N2) +
             eval(N1, N2, Number<CurrentDim0>(), Number<CurrentDim1 - 1>());
    }
    template <int CurrentDim0>
    typename promote<T, U>::V eval(const int N1, const int N2,
                                   const Number<CurrentDim0> &,
                                   const Number<1> &) const {
      return iterA(CurrentDim0 - 1, 0, N1) * iterB(0, CurrentDim0 - 1, N2) +
             eval(N1, N2, Number<CurrentDim0 - 1>(), Number<Dim13>());
    }
    typename promote<T, U>::V eval(const int N1, const int N2,
                                   const Number<1> &, const Number<1> &) const {
      return iterA(0, 0, N1) * iterB(0, 0, N2);
    }

  public:
    Tensor3_times_Tensor3_12_21(
      const Tensor3_Expr<A, T, Dim04, Dim13, Dim2, i, j, k> &a,
      const Tensor3_Expr<B, U, Dim13, Dim04, Dim5, j, i, l> &b)
        : iterA(a), iterB(b)
    {}
    typename promote<T, U>::V operator()(const int &N1, const int &N2) const
    {
      return eval(N1, N2, Number<Dim04>(), Number<Dim13>());
    }
  };

  template <class A, class B, class T, class U, int Dim04, int Dim13, int Dim2,
            int Dim5, char i, char j, char k, char l>
  Tensor2_Expr<Tensor3_times_Tensor3_12_21<A, B, T, U, Dim04, Dim13, Dim2, Dim5,
                                           i, j, k, l>,
               typename promote<T, U>::V, Dim2, Dim5, k, l>
  operator*(const Tensor3_Expr<A, T, Dim04, Dim13, Dim2, i, j, k> &a,
            const Tensor3_Expr<B, U, Dim13, Dim04, Dim5, j, i, l> &b) {
    using TensorExpr = Tensor3_times_Tensor3_12_21<A, B, T, U, Dim04, Dim13,
                                                   Dim2, Dim5, i, j, k, l>;
    return Tensor2_Expr<TensorExpr, typename promote<T, U>::V, Dim2, Dim5, k,
                        l>(TensorExpr(a, b));
  };

    /* A(j,l,k)*B(i,k,l) -> Tensor2 */

  template <class A, class B, class T, class U, int Dim15, int Dim24, int Dim0,
            int Dim3, char i, char j, char k, char l>
  class Tensor3_times_Tensor3_23_32 {
    Tensor3_Expr<A, T, Dim0, Dim15, Dim24, i, k, l> iterA;
    Tensor3_Expr<B, U, Dim3, Dim24, Dim15, j, l, k> iterB;

    template <int CurrentDim0, int CurrentDim1>
    typename promote<T, U>::V eval(const int N1, const int N2,
                                   const Number<CurrentDim0> &,
                                   const Number<CurrentDim1> &) const {
      return iterA(N1, CurrentDim0 - 1, CurrentDim1 - 1) *
                 iterB(N2, CurrentDim1 - 1, CurrentDim0 - 1) +
             eval(N1, N2, Number<CurrentDim0>(), Number<CurrentDim1 - 1>());
    }
    template <int CurrentDim0>
    typename promote<T, U>::V eval(const int N1, const int N2,
                                   const Number<CurrentDim0> &,
                                   const Number<1> &) const {
      return iterA(N1, CurrentDim0 - 1, 0) * iterB(N2, 0, CurrentDim0 - 1) +
             eval(N1, N2, Number<CurrentDim0 - 1>(), Number<Dim24>());
    }
    typename promote<T, U>::V eval(const int N1, const int N2,
                                   const Number<1> &, const Number<1> &) const {
      return iterA(N1, 0, 0) * iterB(N2, 0, 0);
    }

  public:
    Tensor3_times_Tensor3_23_32(
      const Tensor3_Expr<A, T, Dim0, Dim15, Dim24, i, k, l> &a,
      const Tensor3_Expr<B, U, Dim3, Dim24, Dim15, j, l, k> &b)
        : iterA(a), iterB(b)
    {}
    typename promote<T, U>::V operator()(const int &N1, const int &N2) const
    {
      return eval(N1, N2, Number<Dim15>(), Number<Dim24>());
    }
  };

  template <class A, class B, class T, class U, int Dim15, int Dim24, int Dim0,
            int Dim3, char i, char j, char k, char l>
  Tensor2_Expr<Tensor3_times_Tensor3_23_32<A, B, T, U, Dim15, Dim24, Dim0, Dim3,
                                           i, j, k, l>,
               typename promote<T, U>::V, Dim0, Dim3, i, j>
  operator*(const Tensor3_Expr<A, T, Dim0, Dim15, Dim24, i, k, l> &a,
            const Tensor3_Expr<B, U, Dim3, Dim24, Dim15, j, l, k> &b) {
    using TensorExpr = Tensor3_times_Tensor3_23_32<A, B, T, U, Dim15, Dim24,
                                                   Dim0, Dim3, i, j, k, l>;
    return Tensor2_Expr<TensorExpr, typename promote<T, U>::V, Dim0, Dim3, i,
                        j>(TensorExpr(a, b));
  };

  /* A(l,i,j)*B(k,i,j) -> Tensor2 */

  template <class A, class B, class T, class U, int Dim14, int Dim25, int Dim0,
            int Dim3, char i, char j, char k, char l>
  class Tensor3_times_Tensor3_23_23 {
    Tensor3_Expr<A, T, Dim0, Dim14, Dim25, k, i, j> iterA;
    Tensor3_Expr<B, U, Dim3, Dim14, Dim25, l, i, j> iterB;

    template <int CurrentDim0, int CurrentDim1>
    typename promote<T, U>::V eval(const int N1, const int N2,
                                   const Number<CurrentDim0> &,
                                   const Number<CurrentDim1> &) const {
      return iterA(N1, CurrentDim0 - 1, CurrentDim1 - 1) *
                 iterB(N2, CurrentDim0 - 1, CurrentDim1 - 1) +
             eval(N1, N2, Number<CurrentDim0>(), Number<CurrentDim1 - 1>());
    }
    template <int CurrentDim0>
    typename promote<T, U>::V eval(const int N1, const int N2,
                                   const Number<CurrentDim0> &,
                                   const Number<1> &) const {
      return iterA(N1, CurrentDim0 - 1, 0) * iterB(N2, CurrentDim0 - 1, 0) +
             eval(N1, N2, Number<CurrentDim0 - 1>(), Number<Dim25>());
    }
    typename promote<T, U>::V eval(const int N1, const int N2,
                                   const Number<1> &, const Number<1> &) const {
      return iterA(N1, 0, 0) * iterB(N2, 0, 0);
    }

  public:
    Tensor3_times_Tensor3_23_23(
      const Tensor3_Expr<A, T, Dim0, Dim14, Dim25, k, i, j> &a,
      const Tensor3_Expr<B, U, Dim3, Dim14, Dim25, l, i, j> &b)
        : iterA(a), iterB(b)
    {}
    typename promote<T, U>::V operator()(const int &N1, const int &N2) const
    {
      return eval(N1, N2, Number<Dim14>(), Number<Dim25>());
    }
  };

  template <class A, class B, class T, class U, int Dim14, int Dim25, int Dim0,
            int Dim3, char i, char j, char k, char l>
  Tensor2_Expr<Tensor3_times_Tensor3_23_23<A, B, T, U, Dim14, Dim25, Dim0, Dim3,
                                           i, j, k, l>,
               typename promote<T, U>::V, Dim0, Dim3, k, l>
  operator*(const Tensor3_Expr<A, T, Dim0, Dim14, Dim25, k, i, j> &a,
            const Tensor3_Expr<B, U, Dim3, Dim14, Dim25, l, i, j> &b) {
    using TensorExpr = Tensor3_times_Tensor3_23_23<A, B, T, U, Dim14, Dim25,
                                                   Dim0, Dim3, i, j, k, l>;
    return Tensor2_Expr<TensorExpr, typename promote<T, U>::V, Dim0, Dim3, k,
                        l>(TensorExpr(a, b));
  };

    /* A(i,j,k)*B(i,j,l) -> Tensor2 */

  template <class A, class B, class T, class U, int Dim03, int Dim14, int Dim2,
            int Dim5, char i, char j, char k, char l>
  class Tensor3_times_Tensor3_12_12 {
    Tensor3_Expr<A, T, Dim03, Dim14, Dim2, i, j, k> iterA;
    Tensor3_Expr<B, U, Dim03, Dim14, Dim5, i, j, l> iterB;

    template <int CurrentDim0, int CurrentDim1>
    typename promote<T, U>::V eval(const int N1, const int N2,
                                   const Number<CurrentDim0> &,
                                   const Number<CurrentDim1> &) const {
      return iterA(CurrentDim0 - 1, CurrentDim1 - 1, N1) *
                 iterB(CurrentDim0 - 1, CurrentDim1 - 1, N2) +
             eval(N1, N2, Number<CurrentDim0>(), Number<CurrentDim1 - 1>());
    }
    template <int CurrentDim0>
    typename promote<T, U>::V eval(const int N1, const int N2,
                                   const Number<CurrentDim0> &,
                                   const Number<1> &) const {
      return iterA(CurrentDim0 - 1, 0, N1) * iterB(CurrentDim0 - 1, 0, N2) +
             eval(N1, N2, Number<CurrentDim0 - 1>(), Number<Dim14>());
    }
    typename promote<T, U>::V eval(const int N1, const int N2,
                                   const Number<1> &, const Number<1> &) const {
      return iterA(0, 0, N1) * iterB(0, 0, N2);
    }

  public:
    Tensor3_times_Tensor3_12_12(
      const Tensor3_Expr<A, T, Dim03, Dim14, Dim2, i, j, k> &a,
      const Tensor3_Expr<B, U, Dim03, Dim14, Dim5, i, j, l> &b)
        : iterA(a), iterB(b)
    {}
    typename promote<T, U>::V operator()(const int &N1, const int &N2) const
    {
      return eval(N1, N2, Number<Dim03>(), Number<Dim14>());
    }
  };

  template <class A, class B, class T, class U, int Dim03, int Dim14, int Dim2,
            int Dim5, char i, char j, char k, char l>
  Tensor2_Expr<Tensor3_times_Tensor3_12_12<A, B, T, U, Dim03, Dim14, Dim2, Dim5,
                                           i, j, k, l>,
               typename promote<T, U>::V, Dim2, Dim5, k, l>
  operator*(const Tensor3_Expr<A, T, Dim03, Dim14, Dim2, i, j, k> &a,
            const Tensor3_Expr<B, U, Dim03, Dim14, Dim5, i, j, l> &b) {
    using TensorExpr = Tensor3_times_Tensor3_12_12<A, B, T, U, Dim03, Dim14,
                                                   Dim2, Dim5, i, j, k, l>;
    return Tensor2_Expr<TensorExpr, typename promote<T, U>::V, Dim2, Dim5, k,
                        l>(TensorExpr(a, b));
  };

  /* A(i,j,k)*B(k,l,m) -> Tensor4 */

  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim23,
            int Dim4, int Dim5, char i, char j, char k, char l, char m>
  class Tensor3_times_Tensor3_21
  {
    Tensor3_Expr<A, T, Dim0, Dim1, Dim23, i, j, k> iterA;
    Tensor3_Expr<B, U, Dim23, Dim4, Dim5, k, l, m> iterB;

    template <int CurrentDim>
    typename promote<T, U>::V
    eval(const int N1, const int N2, const int N3, const int N4,
         const Number<CurrentDim> &) const
    {
      return iterA(N1, N2, CurrentDim - 1) * iterB(CurrentDim - 1, N3, N4)
             + eval(N1, N2, N3, N4, Number<CurrentDim - 1>());
    }
    typename promote<T, U>::V eval(const int N1, const int N2, const int N3,
                                   const int N4, const Number<1> &) const
    {
      return iterA(N1, N2, 0) * iterB(0, N3, N4);
    }

  public:
    Tensor3_times_Tensor3_21(
      const Tensor3_Expr<A, T, Dim0, Dim1, Dim23, i, j, k> &a,
      const Tensor3_Expr<B, U, Dim23, Dim4, Dim5, k, l, m> &b)
        : iterA(a), iterB(b)
    {}
    typename promote<T, U>::V operator()(const int &N1, const int &N2,
                                         const int &N3, const int &N4) const
    {
      return eval(N1, N2, N3, N4, Number<Dim23>());
    }
  };

  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim23,
            int Dim4, int Dim5, char i, char j, char k, char l, char m>
  Tensor4_Expr<Tensor3_times_Tensor3_21<A, B, T, U, Dim0, Dim1, Dim23, Dim4,
                                        Dim5, i, j, k, l, m>,
               typename promote<T, U>::V, Dim0, Dim1, Dim4, Dim5, i, j, l, m>
  operator*(const Tensor3_Expr<A, T, Dim0, Dim1, Dim23, i, j, k> &a,
            const Tensor3_Expr<B, U, Dim23, Dim4, Dim5, k, l, m> &b)
  {
    using TensorExpr = Tensor3_times_Tensor3_21<A, B, T, U, Dim0, Dim1, Dim23,
                                                Dim4, Dim5, i, j, k, l, m>;
    return Tensor4_Expr<TensorExpr, typename promote<T, U>::V, Dim0, Dim1,
                        Dim4, Dim5, i, j, l, m>(TensorExpr(a, b));
  };
}
