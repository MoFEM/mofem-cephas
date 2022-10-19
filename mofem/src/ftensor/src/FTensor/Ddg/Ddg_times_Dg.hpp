/**
 * @file Ddg_times_Dg.hpp
 * @brief Ddg times Dg
 * @date 2022-10-19
 *
 * @copyright Copyright (c) 2022
 *
 */

#pragma once

namespace FTensor {
/* A(i,j,k,l)*B(k,l,m) */

template <class A, class B, class T, class U, int Dim01, int Dim23, int Dim4,
          char i, char j, char k, char l, char m>
class Ddg_times_Dg_23 {
  const Ddg_Expr<A, T, Dim01, Dim23, i, j, k, l> iterA;
  const Dg_Expr<B, U, Dim23, Dim4, k, l, m> iterB;

  template <int Current_Dim0, int Current_Dim1>
  typename promote<T, U>::V eval(const int N1, const int N2, const int N3,
                                 const Number<Current_Dim0> &,
                                 const Number<Current_Dim1> &) const {
    return iterA(N1, N2, Current_Dim0 - 1, Current_Dim1 - 1) *
               iterB(Current_Dim0 - 1, Current_Dim1 - 1, N3) +
           eval(N1, N2, N3, Number<Current_Dim0 - 1>(), Number<Current_Dim1>());
  }
  template <int Current_Dim1>
  typename promote<T, U>::V eval(const int N1, const int N2, const int N3,
                                 const Number<1> &,
                                 const Number<Current_Dim1> &) const {
    return iterA(N1, N2, 0, Current_Dim1 - 1) * iterB(0, Current_Dim1 - 1, N3) +
           eval(N1, N2, N3, Number<Dim23>(), Number<Current_Dim1 - 1>());
  }
  typename promote<T, U>::V eval(const int N1, const int N2, const int N3,
                                 const Number<1> &, const Number<1> &) const {
    return iterA(N1, N2, 0, 0) * iterB(0, 0, N3);
  }

public:
  typename promote<T, U>::V operator()(const int N1, const int N2,
                                       const int N3) const {
    return eval(N1, N2, N3, Number<Dim23>(), Number<Dim23>());
  }

  Ddg_times_Dg_23(const Ddg_Expr<A, T, Dim01, Dim23, i, j, k, l> &a,
                  const Dg_Expr<B, U, Dim23, Dim4, k, l, m> &b)
      : iterA(a), iterB(b) {}
};

template <class A, class B, class T, class U, int Dim01, int Dim23, int Dim4,
          char i, char j, char k, char l, char m>
inline auto operator*(const Ddg_Expr<A, T, Dim01, Dim23, i, j, k, l> &a,
                      const Dg_Expr<B, U, Dim23, Dim4, k, l, m> &b) {
  using TensorExpr =
      Ddg_times_Dg_23<A, B, T, U, Dim01, Dim23, Dim4, i, j, k, l, m>;
  return Dg_Expr<TensorExpr, typename promote<T, U>::V, Dim01, Dim4, i, j, m>(
      TensorExpr(a, b));
}

/* B(k,l, m)*A(i,j,k,l) */

template <class A, class B, class T, class U, int Dim01, int Dim23, int Dim4,
          char i, char j, char k, char l, char m>
inline auto operator*(const Dg_Expr<B, U, Dim23, Dim4, k, l, m> &b,
                      const Ddg_Expr<A, T, Dim01, Dim23, i, j, k, l> &a) {
  using TensorExpr =
      Ddg_times_Dg_23<A, B, T, U, Dim01, Dim23, Dim4, i, j, k, l, m>;
  return Dg_Expr<TensorExpr, typename promote<T, U>::V, Dim01, Dim4, i, j, m>(
      TensorExpr(a, b));
}

} // namespace FTensor