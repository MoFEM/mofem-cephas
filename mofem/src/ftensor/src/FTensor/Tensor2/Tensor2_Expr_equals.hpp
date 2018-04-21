/* Various assignment operators. */

/* T2=T2 */

#pragma once

namespace FTensor
{
  template <class A, class B, class U, int Dim0, int Dim1, char i, char j,
            int Current_Dim0, int Current_Dim1>
  inline void
  T2_equals_T2(A &iter, const Tensor2_Expr<B, U, Dim0, Dim1, i, j> &result,
               const Number<Current_Dim0> &, const Number<Current_Dim1> &)
  {
    iter(Current_Dim0 - 1, Current_Dim1 - 1)
      = result(Current_Dim0 - 1, Current_Dim1 - 1);
    T2_equals_T2(iter, result, Number<Current_Dim0 - 1>(),
                 Number<Current_Dim1>());
  }

  template <class A, class B, class U, int Dim0, int Dim1, char i, char j,
            int Current_Dim1>
  inline void
  T2_equals_T2(A &iter, const Tensor2_Expr<B, U, Dim0, Dim1, i, j> &result,
               const Number<1> &, const Number<Current_Dim1> &)
  {
    iter(0, Current_Dim1 - 1) = result(0, Current_Dim1 - 1);
    T2_equals_T2(iter, result, Number<Dim0>(), Number<Current_Dim1 - 1>());
  }

  template <class A, class B, class U, int Dim0, int Dim1, char i, char j>
  inline void
  T2_equals_T2(A &iter, const Tensor2_Expr<B, U, Dim0, Dim1, i, j> &result,
               const Number<1> &, const Number<1> &)
  {
    iter(0, 0) = result(0, 0);
  }

  template <class A, class T, int Tensor_Dim0, int Tensor_Dim1, int Dim0,
            int Dim1, char i, char j, Layout layout>
  template <class B, class U>
  inline Tensor2_Expr<Tensor2<A, Tensor_Dim0, Tensor_Dim1, layout>, T, Dim0,
                      Dim1, i, j> &
  Tensor2_Expr<Tensor2<A, Tensor_Dim0, Tensor_Dim1, layout>, T, Dim0, Dim1, i,
               j>::
  operator=(const Tensor2_Expr<B, U, Dim0, Dim1, i, j> &result)
  {
    T2_equals_T2(iter, result, Number<Dim0>(), Number<Dim1>());
    return *this;
  }

  /* ADOL-C assignment T2<<=T2 (left) */

  template <class A, class B, class U, int Dim0, int Dim1, char i, char j,
            int Current_Dim0, int Current_Dim1>
  inline void
  T2_equals_adolc_left_T2(A &iter,
                          const Tensor2_Expr<B, U, Dim0, Dim1, i, j> &result,
                          const Number<Current_Dim0> &,
                          const Number<Current_Dim1> &)
  {
    iter(Current_Dim0 - 1, Current_Dim1 - 1)
      <<= result(Current_Dim0 - 1, Current_Dim1 - 1);
    T2_equals_adolc_left_T2(iter, result, Number<Current_Dim0 - 1>(),
                            Number<Current_Dim1>());
  }

  template <class A, class B, class U, int Dim0, int Dim1, char i, char j,
            int Current_Dim1>
  inline void
  T2_equals_adolc_left_T2(A &iter,
                          const Tensor2_Expr<B, U, Dim0, Dim1, i, j> &result,
                          const Number<1> &, const Number<Current_Dim1> &)
  {
    iter(0, Current_Dim1 - 1) <<= result(0, Current_Dim1 - 1);
    T2_equals_adolc_left_T2(iter, result, Number<Dim0>(),
                            Number<Current_Dim1 - 1>());
  }

  template <class A, class B, class U, int Dim0, int Dim1, char i, char j>
  inline void
  T2_equals_adolc_left_T2(A &iter,
                          const Tensor2_Expr<B, U, Dim0, Dim1, i, j> &result,
                          const Number<1> &, const Number<1> &)
  {
    iter(0, 0) <<= result(0, 0);
  }

  template <class A, class T, int Tensor_Dim0, int Tensor_Dim1, int Dim0,
            int Dim1, char i, char j, Layout layout>
  template <class B, class U>
  inline Tensor2_Expr<Tensor2<A, Tensor_Dim0, Tensor_Dim1, layout>, T, Dim0,
                      Dim1, i, j> &
  Tensor2_Expr<Tensor2<A, Tensor_Dim0, Tensor_Dim1, layout>, T, Dim0, Dim1, i,
               j>::
  operator<<=(const Tensor2_Expr<B, U, Dim0, Dim1, i, j> &result)
  {
    T2_equals_adolc_left_T2(iter, result, Number<Dim0>(), Number<Dim1>());
    return *this;
  }

  /* ADOL-C assignment T2>>=T2 (right) */

  template <class A, class B, class U, int Dim0, int Dim1, char i, char j,
            int Current_Dim0, int Current_Dim1>
  inline void
  T2_equals_adolc_rigth_T2(A &iter,
                           const Tensor2_Expr<B, U, Dim0, Dim1, i, j> &result,
                           const Number<Current_Dim0> &,
                           const Number<Current_Dim1> &)
  {
    iter(Current_Dim0 - 1, Current_Dim1 - 1)
      <<= result(Current_Dim0 - 1, Current_Dim1 - 1);
    T2_equals_adolc_rigth_T2(iter, result, Number<Current_Dim0 - 1>(),
                             Number<Current_Dim1>());
  }

  template <class A, class B, class U, int Dim0, int Dim1, char i, char j,
            int Current_Dim1>
  inline void
  T2_equals_adolc_rigth_T2(A &iter,
                           const Tensor2_Expr<B, U, Dim0, Dim1, i, j> &result,
                           const Number<1> &, const Number<Current_Dim1> &)
  {
    iter(0, Current_Dim1 - 1) <<= result(0, Current_Dim1 - 1);
    T2_equals_adolc_rigth_T2(iter, result, Number<Dim0>(),
                             Number<Current_Dim1 - 1>());
  }

  template <class A, class B, class U, int Dim0, int Dim1, char i, char j>
  inline void
  T2_equals_adolc_rigth_T2(A &iter,
                           const Tensor2_Expr<B, U, Dim0, Dim1, i, j> &result,
                           const Number<1> &, const Number<1> &)
  {
    iter(0, 0) <<= result(0, 0);
  }

  template <class A, class T, int Tensor_Dim0, int Tensor_Dim1, int Dim0,
            int Dim1, char i, char j, Layout layout>
  template <class B, class U>
  inline Tensor2_Expr<Tensor2<A, Tensor_Dim0, Tensor_Dim1, layout>, T, Dim0,
                      Dim1, i, j> &
  Tensor2_Expr<Tensor2<A, Tensor_Dim0, Tensor_Dim1, layout>, T, Dim0, Dim1, i,
               j>::
  operator>>=(const Tensor2_Expr<B, U, Dim0, Dim1, i, j> &result)
  {
    T2_equals_adolc_rigth_T2(iter, result, Number<Dim0>(), Number<Dim1>());
    return *this;
  }

  /* T2=T2_Expr(T2) ADOL-C assignment */

  template <class A, class T, int Tensor_Dim0, int Tensor_Dim1, int Dim0,
            int Dim1, char i, char j, Layout layout>
  inline Tensor2_Expr<Tensor2<A, Tensor_Dim0, Tensor_Dim1, layout>, T, Dim0,
                      Dim1, i, j> &
  Tensor2_Expr<Tensor2<A, Tensor_Dim0, Tensor_Dim1, layout>, T, Dim0, Dim1, i,
               j>::
  operator<<=(const Tensor2_Expr<Tensor2<A, Tensor_Dim0, Tensor_Dim1, layout>,
                                 T, Dim0, Dim1, i, j> &result)
  {
    return operator<<=<Tensor2<A, Tensor_Dim0, Tensor_Dim1, layout>, T>(
      result);
  }

  template <class A, class T, int Tensor_Dim0, int Tensor_Dim1, int Dim0,
            int Dim1, char i, char j, Layout layout>
  inline Tensor2_Expr<Tensor2<A, Tensor_Dim0, Tensor_Dim1, layout>, T, Dim0,
                      Dim1, i, j> &
  Tensor2_Expr<Tensor2<A, Tensor_Dim0, Tensor_Dim1, layout>, T, Dim0, Dim1, i,
               j>::
  operator>>=(const Tensor2_Expr<Tensor2<A, Tensor_Dim0, Tensor_Dim1, layout>,
                                 T, Dim0, Dim1, i, j> &result)
  {
    return operator>>=<Tensor2<A, Tensor_Dim0, Tensor_Dim1, layout>, T>(
      result);
  }

  /* T2=T2_Expr(T2) */

  template <class A, class T, int Tensor_Dim0, int Tensor_Dim1, int Dim0,
            int Dim1, char i, char j, Layout layout>
  inline Tensor2_Expr<Tensor2<A, Tensor_Dim0, Tensor_Dim1, layout>, T, Dim0,
                      Dim1, i, j> &
  Tensor2_Expr<Tensor2<A, Tensor_Dim0, Tensor_Dim1, layout>, T, Dim0, Dim1, i,
               j>::
  operator=(const Tensor2_Expr<Tensor2<A, Tensor_Dim0, Tensor_Dim1, layout>, T,
                               Dim0, Dim1, i, j> &result)
  {
    return operator=<Tensor2<A, Tensor_Dim0, Tensor_Dim1, layout>, T>(result);
  }

  /* T2+=T2 */

  template <class A, class B, class U, int Dim0, int Dim1, char i, char j,
            int Current_Dim0, int Current_Dim1>
  inline void
  T2_plus_equals_T2(A &iter,
                    const Tensor2_Expr<B, U, Dim0, Dim1, i, j> &result,
                    const Number<Current_Dim0> &, const Number<Current_Dim1> &)
  {
    iter(Current_Dim0 - 1, Current_Dim1 - 1)
      += result(Current_Dim0 - 1, Current_Dim1 - 1);
    T2_plus_equals_T2(iter, result, Number<Current_Dim0 - 1>(),
                      Number<Current_Dim1>());
  }

  template <class A, class B, class U, int Dim0, int Dim1, char i, char j,
            int Current_Dim1>
  inline void
  T2_plus_equals_T2(A &iter,
                    const Tensor2_Expr<B, U, Dim0, Dim1, i, j> &result,
                    const Number<1> &, const Number<Current_Dim1> &)
  {
    iter(0, Current_Dim1 - 1) += result(0, Current_Dim1 - 1);
    T2_plus_equals_T2(iter, result, Number<Dim0>(),
                      Number<Current_Dim1 - 1>());
  }

  template <class A, class B, class U, int Dim0, int Dim1, char i, char j>
  inline void
  T2_plus_equals_T2(A &iter,
                    const Tensor2_Expr<B, U, Dim0, Dim1, i, j> &result,
                    const Number<1> &, const Number<1> &)
  {
    iter(0, 0) += result(0, 0);
  }

  template <class A, class T, int Tensor_Dim0, int Tensor_Dim1, int Dim0,
            int Dim1, char i, char j, Layout layout>
  template <class B, class U>
  inline Tensor2_Expr<Tensor2<A, Tensor_Dim0, Tensor_Dim1, layout>, T, Dim0,
                      Dim1, i, j> &
  Tensor2_Expr<Tensor2<A, Tensor_Dim0, Tensor_Dim1, layout>, T, Dim0, Dim1, i,
               j>::
  operator+=(const Tensor2_Expr<B, U, Dim0, Dim1, i, j> &result)
  {
    T2_plus_equals_T2(iter, result, Number<Dim0>(), Number<Dim1>());
    return *this;
  }

  /* T2-=T2 */

  template <class A, class B, class U, int Dim0, int Dim1, char i, char j,
            int Current_Dim0, int Current_Dim1>
  inline void
  T2_minus_equals_T2(A &iter,
                     const Tensor2_Expr<B, U, Dim0, Dim1, i, j> &result,
                     const Number<Current_Dim0> &,
                     const Number<Current_Dim1> &)
  {
    iter(Current_Dim0 - 1, Current_Dim1 - 1)
      -= result(Current_Dim0 - 1, Current_Dim1 - 1);
    T2_minus_equals_T2(iter, result, Number<Current_Dim0 - 1>(),
                       Number<Current_Dim1>());
  }

  template <class A, class B, class U, int Dim0, int Dim1, char i, char j,
            int Current_Dim1>
  inline void
  T2_minus_equals_T2(A &iter,
                     const Tensor2_Expr<B, U, Dim0, Dim1, i, j> &result,
                     const Number<1> &, const Number<Current_Dim1> &)
  {
    iter(0, Current_Dim1 - 1) -= result(0, Current_Dim1 - 1);
    T2_minus_equals_T2(iter, result, Number<Dim0>(),
                       Number<Current_Dim1 - 1>());
  }

  template <class A, class B, class U, int Dim0, int Dim1, char i, char j>
  inline void
  T2_minus_equals_T2(A &iter,
                     const Tensor2_Expr<B, U, Dim0, Dim1, i, j> &result,
                     const Number<1> &, const Number<1> &)
  {
    iter(0, 0) -= result(0, 0);
  }

  template <class A, class T, int Tensor_Dim0, int Tensor_Dim1, int Dim0,
            int Dim1, char i, char j, Layout layout>
  template <class B, class U>
  inline Tensor2_Expr<Tensor2<A, Tensor_Dim0, Tensor_Dim1, layout>, T, Dim0,
                      Dim1, i, j> &
  Tensor2_Expr<Tensor2<A, Tensor_Dim0, Tensor_Dim1, layout>, T, Dim0, Dim1, i,
               j>::
  operator-=(const Tensor2_Expr<B, U, Dim0, Dim1, i, j> &result)
  {
    T2_minus_equals_T2(iter, result, Number<Dim0>(), Number<Dim1>());
    return *this;
  }

  /* This is for when the indices are switched (i,j) -> (j,i). */

  /* T2(i,j)=T2(j,i) */

  template <class A, class B, class U, int Dim0, int Dim1, char i, char j,
            int Current_Dim0, int Current_Dim1>
  inline void
  T2_equals_switched_T2(A &iter,
                        const Tensor2_Expr<B, U, Dim0, Dim1, j, i> &result,
                        const Number<Current_Dim0> &,
                        const Number<Current_Dim1> &)
  {
    iter(Current_Dim0 - 1, Current_Dim1 - 1)
      = result(Current_Dim1 - 1, Current_Dim0 - 1);
    T2_equals_switched_T2(iter, result, Number<Current_Dim0 - 1>(),
                          Number<Current_Dim1>());
  }

  template <class A, class B, class U, int Dim0, int Dim1, char i, char j,
            int Current_Dim1>
  inline void
  T2_equals_switched_T2(A &iter,
                        const Tensor2_Expr<B, U, Dim0, Dim1, j, i> &result,
                        const Number<1> &, const Number<Current_Dim1> &)
  {
    iter(0, Current_Dim1 - 1) = result(Current_Dim1 - 1, 0);
    T2_equals_switched_T2(iter, result, Number<Dim0>(),
                          Number<Current_Dim1 - 1>());
  }

  template <class A, class B, class U, int Dim0, int Dim1, char i, char j>
  inline void
  T2_equals_switched_T2(A &iter,
                        const Tensor2_Expr<B, U, Dim0, Dim1, j, i> &result,
                        const Number<1> &, const Number<1> &)
  {
    iter(0, 0) = result(0, 0);
  }

  template <class A, class T, int Tensor_Dim0, int Tensor_Dim1, int Dim0,
            int Dim1, char i, char j, Layout layout>
  template <class B, class U>
  inline Tensor2_Expr<Tensor2<A, Tensor_Dim0, Tensor_Dim1, layout>, T, Dim0,
                      Dim1, i, j> &
  Tensor2_Expr<Tensor2<A, Tensor_Dim0, Tensor_Dim1, layout>, T, Dim0, Dim1, i,
               j>::
  operator=(const Tensor2_Expr<B, U, Dim0, Dim1, j, i> &result)
  {
    T2_equals_switched_T2(iter, result, Number<Dim0>(), Number<Dim1>());
    return *this;
  }

  /* T2(i,j)+=T2(j,i) */

  template <class A, class B, class U, int Dim0, int Dim1, char i, char j,
            int Current_Dim0, int Current_Dim1>
  inline void T2_plus_equals_switched_T2(
    A &iter, const Tensor2_Expr<B, U, Dim0, Dim1, j, i> &result,
    const Number<Current_Dim0> &, const Number<Current_Dim1> &)
  {
    iter(Current_Dim0 - 1, Current_Dim1 - 1)
      += result(Current_Dim1 - 1, Current_Dim0 - 1);
    T2_plus_equals_switched_T2(iter, result, Number<Current_Dim0 - 1>(),
                               Number<Current_Dim1>());
  }

  template <class A, class B, class U, int Dim0, int Dim1, char i, char j,
            int Current_Dim1>
  inline void T2_plus_equals_switched_T2(
    A &iter, const Tensor2_Expr<B, U, Dim0, Dim1, j, i> &result,
    const Number<1> &, const Number<Current_Dim1> &)
  {
    iter(0, Current_Dim1 - 1) += result(Current_Dim1 - 1, 0);
    T2_plus_equals_switched_T2(iter, result, Number<Dim0>(),
                               Number<Current_Dim1 - 1>());
  }

  template <class A, class B, class U, int Dim0, int Dim1, char i, char j>
  inline void T2_plus_equals_switched_T2(
    A &iter, const Tensor2_Expr<B, U, Dim0, Dim1, j, i> &result,
    const Number<1> &, const Number<1> &)
  {
    iter(0, 0) += result(0, 0);
  }

  template <class A, class T, int Tensor_Dim0, int Tensor_Dim1, int Dim0,
            int Dim1, char i, char j, Layout layout>
  template <class B, class U>
  inline Tensor2_Expr<Tensor2<A, Tensor_Dim0, Tensor_Dim1, layout>, T, Dim0,
                      Dim1, i, j> &
  Tensor2_Expr<Tensor2<A, Tensor_Dim0, Tensor_Dim1, layout>, T, Dim0, Dim1, i,
               j>::
  operator+=(const Tensor2_Expr<B, U, Dim0, Dim1, j, i> &result)
  {
    T2_plus_equals_switched_T2(iter, result, Number<Dim0>(), Number<Dim1>());
    return *this;
  }

  /* T2(i,j)-=T2(j,i) */

  template <class A, class B, class U, int Dim0, int Dim1, char i, char j,
            int Current_Dim0, int Current_Dim1>
  inline void T2_minus_equals_switched_T2(
    A &iter, const Tensor2_Expr<B, U, Dim0, Dim1, j, i> &result,
    const Number<Current_Dim0> &, const Number<Current_Dim1> &)
  {
    iter(Current_Dim0 - 1, Current_Dim1 - 1)
      -= result(Current_Dim1 - 1, Current_Dim0 - 1);
    T2_minus_equals_switched_T2(iter, result, Number<Current_Dim0 - 1>(),
                                Number<Current_Dim1>());
  }

  template <class A, class B, class U, int Dim0, int Dim1, char i, char j,
            int Current_Dim1>
  inline void T2_minus_equals_switched_T2(
    A &iter, const Tensor2_Expr<B, U, Dim0, Dim1, j, i> &result,
    const Number<1> &, const Number<Current_Dim1> &)
  {
    iter(0, Current_Dim1 - 1) -= result(Current_Dim1 - 1, 0);
    T2_minus_equals_switched_T2(iter, result, Number<Dim0>(),
                                Number<Current_Dim1 - 1>());
  }

  template <class A, class B, class U, int Dim0, int Dim1, char i, char j>
  inline void T2_minus_equals_switched_T2(
    A &iter, const Tensor2_Expr<B, U, Dim0, Dim1, j, i> &result,
    const Number<1> &, const Number<1> &)
  {
    iter(0, 0) -= result(0, 0);
  }

  template <class A, class T, int Tensor_Dim0, int Tensor_Dim1, int Dim0,
            int Dim1, char i, char j, Layout layout>
  template <class B, class U>
  inline Tensor2_Expr<Tensor2<A, Tensor_Dim0, Tensor_Dim1, layout>, T, Dim0,
                      Dim1, i, j> &
  Tensor2_Expr<Tensor2<A, Tensor_Dim0, Tensor_Dim1, layout>, T, Dim0, Dim1, i,
               j>::
  operator-=(const Tensor2_Expr<B, U, Dim0, Dim1, j, i> &result)
  {
    T2_minus_equals_switched_T2(iter, result, Number<Dim0>(), Number<Dim1>());
    return *this;
  }

  /* For int's, double's, etc. */

  /* T2=U */

  template <class A, class U, int Dim0, int Dim1, int Current_Dim0,
            int Current_Dim1, Layout layout>
  inline void
  T2_equals_generic(Tensor2<A, Dim0, Dim1, layout> &iter, const U &u,
                    const Number<Current_Dim0> &, const Number<Current_Dim1> &)
  {
    iter(Current_Dim0 - 1, Current_Dim1 - 1) = u;
    T2_equals_generic(iter, u, Number<Current_Dim0 - 1>(),
                      Number<Current_Dim1>());
  }

  template <class A, class U, int Dim0, int Dim1, int Current_Dim1,
            Layout layout>
  inline void
  T2_equals_generic(Tensor2<A, Dim0, Dim1, layout> &iter, const U &u,
                    const Number<1> &, const Number<Current_Dim1> &)
  {
    iter(0, Current_Dim1 - 1) = u;
    T2_equals_generic(iter, u, Number<Dim0>(), Number<Current_Dim1 - 1>());
  }

  template <class A, class U, int Dim0, int Dim1, Layout layout>
  inline void
  T2_equals_generic(Tensor2<A, Dim0, Dim1, layout> &iter, const U &u,
                    const Number<1> &, const Number<1> &)
  {
    iter(0, 0) = u;
  }

  template <class A, class T, int Tensor_Dim0, int Tensor_Dim1, int Dim0,
            int Dim1, char i, char j, Layout layout>
  template <class U>
  inline Tensor2_Expr<Tensor2<A, Tensor_Dim0, Tensor_Dim1, layout>, T, Dim0,
                      Dim1, i, j> &
  Tensor2_Expr<Tensor2<A, Tensor_Dim0, Tensor_Dim1, layout>, T, Dim0, Dim1, i,
               j>::operator=(const U &u)
  {
    T2_equals_generic(iter, u, Number<Dim0>(), Number<Dim1>());
    return *this;
  }

  /* T2+=U */

  template <class A, class U, int Dim0, int Dim1, int Current_Dim0,
            int Current_Dim1, Layout layout>
  inline void T2_plus_equals_generic(Tensor2<A, Dim0, Dim1, layout> &iter,
                                     const U &u, const Number<Current_Dim0> &,
                                     const Number<Current_Dim1> &)
  {
    iter(Current_Dim0 - 1, Current_Dim1 - 1) += u;
    T2_plus_equals_generic(iter, u, Number<Current_Dim0 - 1>(),
                           Number<Current_Dim1>());
  }

  template <class A, class U, int Dim0, int Dim1, int Current_Dim1,
            Layout layout>
  inline void
  T2_plus_equals_generic(Tensor2<A, Dim0, Dim1, layout> &iter, const U &u,
                         const Number<1> &, const Number<Current_Dim1> &)
  {
    iter(0, Current_Dim1 - 1) += u;
    T2_plus_equals_generic(iter, u, Number<Dim0>(),
                           Number<Current_Dim1 - 1>());
  }

  template <class A, class U, int Dim0, int Dim1, Layout layout>
  inline void
  T2_plus_equals_generic(Tensor2<A, Dim0, Dim1, layout> &iter, const U &u,
                         const Number<1> &, const Number<1> &)
  {
    iter(0, 0) += u;
  }

  template <class A, class T, int Tensor_Dim0, int Tensor_Dim1, int Dim0,
            int Dim1, char i, char j, Layout layout>
  template <class U>
  inline Tensor2_Expr<Tensor2<A, Tensor_Dim0, Tensor_Dim1, layout>, T, Dim0,
                      Dim1, i, j> &
  Tensor2_Expr<Tensor2<A, Tensor_Dim0, Tensor_Dim1, layout>, T, Dim0, Dim1, i,
               j>::operator+=(const U &u)
  {
    T2_plus_equals_generic(iter, u, Number<Dim0>(), Number<Dim1>());
    return *this;
  }

  /* T2-=U */

  template <class A, class U, int Dim0, int Dim1, int Current_Dim0,
            int Current_Dim1, Layout layout>
  inline void T2_minus_equals_generic(Tensor2<A, Dim0, Dim1, layout> &iter,
                                      const U &u, const Number<Current_Dim0> &,
                                      const Number<Current_Dim1> &)
  {
    iter(Current_Dim0 - 1, Current_Dim1 - 1) -= u;
    T2_minus_equals_generic(iter, u, Number<Current_Dim0 - 1>(),
                            Number<Current_Dim1>());
  }

  template <class A, class U, int Dim0, int Dim1, int Current_Dim1,
            Layout layout>
  inline void
  T2_minus_equals_generic(Tensor2<A, Dim0, Dim1, layout> &iter, const U &u,
                          const Number<1> &, const Number<Current_Dim1> &)
  {
    iter(0, Current_Dim1 - 1) -= u;
    T2_minus_equals_generic(iter, u, Number<Dim0>(),
                            Number<Current_Dim1 - 1>());
  }

  template <class A, class U, int Dim0, int Dim1, Layout layout>
  inline void
  T2_minus_equals_generic(Tensor2<A, Dim0, Dim1, layout> &iter, const U &u,
                          const Number<1> &, const Number<1> &)
  {
    iter(0, 0) -= u;
  }

  template <class A, class T, int Tensor_Dim0, int Tensor_Dim1, int Dim0,
            int Dim1, char i, char j, Layout layout>
  template <class U>
  inline Tensor2_Expr<Tensor2<A, Tensor_Dim0, Tensor_Dim1, layout>, T, Dim0,
                      Dim1, i, j> &
  Tensor2_Expr<Tensor2<A, Tensor_Dim0, Tensor_Dim1, layout>, T, Dim0, Dim1, i,
               j>::operator-=(const U &u)
  {
    T2_minus_equals_generic(iter, u, Number<Dim0>(), Number<Dim1>());
    return *this;
  }

  /* T2*=U */

  template <class A, class U, int Dim0, int Dim1, int Current_Dim0,
            int Current_Dim1, Layout layout>
  inline void T2_times_equals_generic(Tensor2<A, Dim0, Dim1, layout> &iter,
                                      const U &u, const Number<Current_Dim0> &,
                                      const Number<Current_Dim1> &)
  {
    iter(Current_Dim0 - 1, Current_Dim1 - 1) *= u;
    T2_times_equals_generic(iter, u, Number<Current_Dim0 - 1>(),
                            Number<Current_Dim1>());
  }

  template <class A, class U, int Dim0, int Dim1, int Current_Dim1,
            Layout layout>
  inline void
  T2_times_equals_generic(Tensor2<A, Dim0, Dim1, layout> &iter, const U &u,
                          const Number<1> &, const Number<Current_Dim1> &)
  {
    iter(0, Current_Dim1 - 1) *= u;
    T2_times_equals_generic(iter, u, Number<Dim0>(),
                            Number<Current_Dim1 - 1>());
  }

  template <class A, class U, int Dim0, int Dim1, Layout layout>
  inline void
  T2_times_equals_generic(Tensor2<A, Dim0, Dim1, layout> &iter, const U &u,
                          const Number<1> &, const Number<1> &)
  {
    iter(0, 0) *= u;
  }

  template <class A, class T, int Tensor_Dim0, int Tensor_Dim1, int Dim0,
            int Dim1, char i, char j, Layout layout>
  template <class U>
  inline Tensor2_Expr<Tensor2<A, Tensor_Dim0, Tensor_Dim1, layout>, T, Dim0,
                      Dim1, i, j> &
  Tensor2_Expr<Tensor2<A, Tensor_Dim0, Tensor_Dim1, layout>, T, Dim0, Dim1, i,
               j>::operator*=(const U &u)
  {
    T2_times_equals_generic(iter, u, Number<Dim0>(), Number<Dim1>());
    return *this;
  }

  /* T2/=U */

  template <class A, class U, int Dim0, int Dim1, int Current_Dim0,
            int Current_Dim1, Layout layout>
  inline void
  T2_divide_equals_generic(Tensor2<A, Dim0, Dim1, layout> &iter, const U &u,
                           const Number<Current_Dim0> &,
                           const Number<Current_Dim1> &)
  {
    iter(Current_Dim0 - 1, Current_Dim1 - 1) /= u;
    T2_divide_equals_generic(iter, u, Number<Current_Dim0 - 1>(),
                             Number<Current_Dim1>());
  }

  template <class A, class U, int Dim0, int Dim1, int Current_Dim1,
            Layout layout>
  inline void
  T2_divide_equals_generic(Tensor2<A, Dim0, Dim1, layout> &iter, const U &u,
                           const Number<1> &, const Number<Current_Dim1> &)
  {
    iter(0, Current_Dim1 - 1) /= u;
    T2_divide_equals_generic(iter, u, Number<Dim0>(),
                             Number<Current_Dim1 - 1>());
  }

  template <class A, class U, int Dim0, int Dim1, Layout layout>
  inline void
  T2_divide_equals_generic(Tensor2<A, Dim0, Dim1, layout> &iter, const U &u,
                           const Number<1> &, const Number<1> &)
  {
    iter(0, 0) /= u;
  }

  template <class A, class T, int Tensor_Dim0, int Tensor_Dim1, int Dim0,
            int Dim1, char i, char j, Layout layout>
  template <class U>
  inline Tensor2_Expr<Tensor2<A, Tensor_Dim0, Tensor_Dim1, layout>, T, Dim0,
                      Dim1, i, j> &
  Tensor2_Expr<Tensor2<A, Tensor_Dim0, Tensor_Dim1, layout>, T, Dim0, Dim1, i,
               j>::operator/=(const U &u)
  {
    T2_divide_equals_generic(iter, u, Number<Dim0>(), Number<Dim1>());
    return *this;
  }

  /* Various assignment operators for Dg_number_rhs_0.  I have
     to explicitly declare the second operator= because otherwise the
     compiler will generate its own and not use the template code. */

  template <class A, class B, class U, int Dim0, int Dim1, char i, char j,
            int N, int Current_Dim0, int Current_Dim1>
  inline void
  T3dgrhs0_equals_T2(A &iter,
                     const Tensor2_Expr<B, U, Dim0, Dim1, i, j> &result,
                     const Number<Current_Dim0> &,
                     const Number<Current_Dim1> &, const Number<N> &)
  {
    iter(N, Current_Dim0 - 1, Current_Dim1 - 1)
      = result(Current_Dim0 - 1, Current_Dim1 - 1);
    T3dgrhs0_equals_T2(iter, result, Number<Current_Dim0 - 1>(),
                       Number<Current_Dim1>(), Number<N>());
  }

  template <class A, class B, class U, int Dim0, int Dim1, char i, char j,
            int N, int Current_Dim1>
  inline void
  T3dgrhs0_equals_T2(A &iter,
                     const Tensor2_Expr<B, U, Dim0, Dim1, i, j> &result,
                     const Number<1> &, const Number<Current_Dim1> &,
                     const Number<N> &)
  {
    iter(N, 0, Current_Dim1 - 1) = result(0, Current_Dim1 - 1);
    T3dgrhs0_equals_T2(iter, result, Number<Dim0>(),
                       Number<Current_Dim1 - 1>(), Number<N>());
  }

  template <class A, class B, class U, int Dim0, int Dim1, char i, char j,
            int N>
  inline void
  T3dgrhs0_equals_T2(A &iter,
                     const Tensor2_Expr<B, U, Dim0, Dim1, i, j> &result,
                     const Number<1> &, const Number<1> &, const Number<N> &)
  {
    iter(N, 0, 0) = result(0, 0);
  }

  template <class A, class T, int Dim0, int Dim1, char i, char j, int N>
  template <class B, class U>
  inline Tensor2_Expr<Dg_number_rhs_0<A, T, N>, T, Dim0, Dim1, i, j> &
  Tensor2_Expr<Dg_number_rhs_0<A, T, N>, T, Dim0, Dim1, i, j>::
  operator=(const Tensor2_Expr<B, U, Dim0, Dim1, i, j> &result)
  {
    T3dgrhs0_equals_T2(iter, result, Number<Dim0>(), Number<Dim1>(),
                       Number<N>());
    return *this;
  }

  template <class A, class T, int Dim0, int Dim1, char i, char j, int N>
  inline Tensor2_Expr<Dg_number_rhs_0<A, T, N>, T, Dim0, Dim1, i, j> &
  Tensor2_Expr<Dg_number_rhs_0<A, T, N>, T, Dim0, Dim1, i, j>::operator=(
    const Tensor2_Expr<Dg_number_rhs_0<A, T, N>, T, Dim0, Dim1, i, j> &result)
  {
    return operator=<Dg_number_rhs_0<A, T, N>, T>(result);
  }
}
