/* Declare a wrapper class for generic rank 4 Tensor expressions. */

#pragma once

#include "Tensor4_plus_Tensor4.hpp"
#include "Tensor4_minus_Tensor4.hpp"
#include "Tensor4_times_Tensor2_symmetric.hpp"
#include "Tensor4_times_Tensor2.hpp"

#include "Tensor4_times_Tensor1.hpp"

namespace FTensor
{
  template <class A, class T, int Dim0, int Dim1, int Dim2, int Dim3, char i,
            char j, char k, char l>
  class Tensor4_Expr
  {
    A iter;

  public:
    Tensor4_Expr(A &a) : iter(a) {}
    T operator()(const int N1, const int N2, const int N3, const int N4) const
    {
      return iter(N1, N2, N3, N4);
    }
  };

  template <class A, class T, int Dim0, int Dim1, int Dim2, int Dim3, char i,
            char j, char k, char l>
  class Tensor4_Expr<Tensor4<A, Dim0, Dim1, Dim2, Dim3>, T, Dim0, Dim1, Dim2,
                     Dim3, i, j, k, l>
  {
    Tensor4<A, Dim0, Dim1, Dim2, Dim3> &iter;

  public:
    Tensor4_Expr(Tensor4<A, Dim0, Dim1, Dim2, Dim3> &a) : iter(a) {}
    T &operator()(const int N1, const int N2, const int N3, const int N4)
    {
      return iter(N1, N2, N3, N4);
    }
    T operator()(const int N1, const int N2, const int N3, const int N4) const
    {
      return iter(N1, N2, N3, N4);
    }

    /* Various assignment operators.  I have to explicitly declare the
       second operator= because otherwise the compiler will generate its
       own and not use the template code. */

    template <class B, class U>
    const Tensor4_Expr<Tensor4<A, Dim0, Dim1, Dim2, Dim3>, T, Dim0, Dim1, Dim2,
                       Dim3, i, j, k, l> &
    operator=(
      const Tensor4_Expr<B, U, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &result);
    const Tensor4_Expr<Tensor4<A, Dim0, Dim1, Dim2, Dim3>, T, Dim0, Dim1, Dim2,
                       Dim3, i, j, k, l> &
    operator=(const Tensor4_Expr<Tensor4<A, Dim0, Dim1, Dim2, Dim3>, T, Dim0,
                                 Dim1, Dim2, Dim3, i, j, k, l> &result);

    template <class B, class U>
    const Tensor4_Expr<Tensor4<A, Dim0, Dim1, Dim2, Dim3>, T, Dim0, Dim1, Dim2,
                       Dim3, i, j, k, l> &
    operator+=(
      const Tensor4_Expr<B, U, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &result);

    template <class B, class U>
    const Tensor4_Expr<Tensor4<A, Dim0, Dim1, Dim2, Dim3>, T, Dim0, Dim1, Dim2,
                       Dim3, i, j, k, l> &
    operator-=(
      const Tensor4_Expr<B, U, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &result);

    template <class B, class U>
    const Tensor4_Expr<Tensor4<A, Dim0, Dim1, Dim2, Dim3>, T, Dim0, Dim1, Dim2,
                       Dim3, i, j, k, l> &
    operator=(
      const Tensor4_Expr<B, U, Dim0, Dim3, Dim1, Dim2, i, l, j, k> &result);

    /* permutations of assignment for tensor4 */

#define EQUAL(D0, D1, D2, D3, I, J, K, L)                                     \
  template <class B, class U>                                                 \
  inline const Tensor4_Expr<Tensor4<A, Dim0, Dim1, Dim2, Dim3>, T, Dim0,      \
                            Dim1, Dim2, Dim3, i, j, k, l>                     \
    &operator=(const Tensor4_Expr<B, U, D0, D1, D2, D3, I, J, K, L> &result);

    EQUAL(Dim2, Dim0, Dim1, Dim3, k, i, j, l); // kijl

#undef EQUAL
  };
}

#include "Tensor4_Expr_equals.hpp"
