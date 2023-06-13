/**
 * @file Tensor4_times_Tensor4_double.hpp
 * @brief Tensor4 times Tensor4 yields Tensor4
 * @date 2022-10-18
 *
 * @copyright Copyright (c) 2022
 *
 */

#pragma once

namespace FTensor {

/* Double contraction. */

/* A(i,j,k,l)*B(k,l,m,n) */

template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2,
          int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m,
          char n>
class Tensor4_times_Tensor4_2345 {
  const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> iterA;
  const Tensor4_Expr<B, U, Dim2, Dim3, Dim4, Dim5, k, l, m, n> iterB;

  template <int Current_Dim0, int Current_Dim1>
  typename promote<T, U>::V eval(const int N1, const int N2, const int N3,
                                 const int N4, const Number<Current_Dim0> &,
                                 const Number<Current_Dim1> &) const {
    return iterA(N1, N2, Current_Dim0 - 1, Current_Dim1 - 1) *
               iterB(Current_Dim0 - 1, Current_Dim1 - 1, N3, N4) +
           eval(N1, N2, N3, N4, Number<Current_Dim0 - 1>(),
                Number<Current_Dim1>());
  }
  template <int Current_Dim1>
  typename promote<T, U>::V eval(const int N1, const int N2, const int N3,
                                 const int N4, const Number<1> &,
                                 const Number<Current_Dim1> &) const {
    return iterA(N1, N2, 0, Current_Dim1 - 1) *
               iterB(0, Current_Dim1 - 1, N3, N4) +
           eval(N1, N2, N3, N4, Number<Dim2>(), Number<Current_Dim1 - 1>());
  }
  typename promote<T, U>::V eval(const int N1, const int N2, const int N3,
                                 const int N4, const Number<1> &,
                                 const Number<1> &) const {
    return iterA(N1, N2, 0, 0) * iterB(0, 0, N3, N4);
  }

public:
  typename promote<T, U>::V operator()(const int N1, const int N2, const int N3,
                                       const int N4) const {
    return eval(N1, N2, N3, N4, Number<Dim2>(), Number<Dim3>());
  }

  Tensor4_times_Tensor4_2345(
      const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a,
      const Tensor4_Expr<B, U, Dim2, Dim3, Dim4, Dim5, k, l, m, n> &b)
      : iterA(a), iterB(b) {}
};

template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2,
          int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m,
          char n>
inline auto
operator*(const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a,
          const Tensor4_Expr<B, U, Dim2, Dim3, Dim4, Dim5, k, l, m, n> &b) {
  using TensorExpr =
      Tensor4_times_Tensor4_2345<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5,
                                 i, j, k, l, m, n>;
  return Tensor4_Expr<TensorExpr, typename promote<T, U>::V, Dim0, Dim1, Dim4,
                      Dim5, i, j, m, n>(TensorExpr(a, b));
}
}