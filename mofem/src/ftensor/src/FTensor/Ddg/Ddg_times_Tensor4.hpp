/* This file has all of the declarations for expressions like
   Ddg*Tensor4 and Tensor4*Ddg, yielding a
   Tensor4. */

#pragma once

namespace FTensor {

/* S(i,j,m,n)*A(m,n,k,l) */

template <class A, class B, class T, class U, int Dim01, int Dim23, int Dim2,
          int Dim3, char i, char j, char k, char l, char m, char n>
class Ddg_times_Tensor4_2301 {
  Ddg_Expr<A, T, Dim01, Dim23, i, j, m, n> iterA;
  Tensor4_Expr<B, U, Dim23, Dim23, Dim2, Dim3, m, n, k, l> iterB;

  template <int Current_Dim0, int Current_Dim1>
  inline typename promote<T, U>::V
  eval(const int N1, const int N2, const int N3, const int N4,
       const Number<Current_Dim0> &, const Number<Current_Dim1> &) const {
    return iterA(N1, N2, Current_Dim0 - 1, Current_Dim1 - 1) *
               iterB(Current_Dim0 - 1, Current_Dim1 - 1, N3, N4) +
           eval(N1, N2, N3, N4, Number<Current_Dim0 - 1>(),
                Number<Current_Dim1>());
  }
  template <int Current_Dim1>
  inline typename promote<T, U>::V
  eval(const int N1, const int N2, const int N3, const int N4,
       const Number<1> &, const Number<Current_Dim1> &) const {
    return iterA(N1, N2, 0, Current_Dim1 - 1) *
               iterB(0, Current_Dim1 - 1, N3, N4) +
           eval(N1, N2, N3, N4, Number<Dim23>(), Number<Current_Dim1 - 1>());
  }
  inline typename promote<T, U>::V eval(const int N1, const int N2,
                                        const int N3, const int N4,
                                        const Number<1> &,
                                        const Number<1> &) const {
    return iterA(N1, N2, 0, 0) * iterB(0, 0, N3, N4);
  }

public:
  Ddg_times_Tensor4_2301(
      const Ddg_Expr<A, T, Dim01, Dim23, i, j, m, n> &a,
      const Tensor4_Expr<B, U, Dim23, Dim23, Dim2, Dim3, m, n, k, l> &b)
      : iterA(a), iterB(b) {}
  typename promote<T, U>::V operator()(const int N1, const int N2, const int N3,
                                       const int N4) const {
    return eval(N1, N2, N3, N4, Number<Dim23>(), Number<Dim23>());
  }
};

template <class A, class B, class T, class U, int Dim01, int Dim23, int Dim2,
          int Dim3, char i, char j, char k, char l, char m, char n>
Tensor4_Expr<Ddg_times_Tensor4_2301<A, B, T, U, Dim01, Dim23, Dim2, Dim3, i, j,
                                    k, l, m, n>,
             typename promote<T, U>::V, Dim01, Dim01, Dim2, Dim3, i, j, k, l>
operator*(const Ddg_Expr<A, T, Dim01, Dim23, i, j, m, n> &a,
          const Tensor4_Expr<B, U, Dim23, Dim23, Dim2, Dim3, m, n, k, l> &b) {
  using TensorExpr = Ddg_times_Tensor4_2301<A, B, T, U, Dim01, Dim23, Dim2,
                                            Dim3, i, j, k, l, m, n>;
  return Tensor4_Expr<TensorExpr, typename promote<T, U>::V, Dim01, Dim01, Dim2,
                      Dim3, i, j, k, l>(TensorExpr(a, b));
}

} // namespace FTensor