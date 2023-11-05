/* This file has all of the declarations for expressions like
   Ddg*Tensor4 and Tensor4*Ddg, yielding a
   Tensor4. */

#pragma once

namespace FTensor {

/* S(i,j,m,n)*A(m,n,k,l) */

template <class A, class B, class T, class U, int Dim01, int Dim23, int Dim2,
          int Dim3, char i, char j, char k, char l, char m, char n>
class Ddg_times_Tensor4_2301_ijkl {
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
  Ddg_times_Tensor4_2301_ijkl(
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
Tensor4_Expr<Ddg_times_Tensor4_2301_ijkl<A, B, T, U, Dim01, Dim23, Dim2, Dim3, i, j,
                                    k, l, m, n>,
             typename promote<T, U>::V, Dim01, Dim01, Dim2, Dim3, i, j, k, l>
operator*(const Ddg_Expr<A, T, Dim01, Dim23, i, j, m, n> &a,
          const Tensor4_Expr<B, U, Dim23, Dim23, Dim2, Dim3, m, n, k, l> &b) {
  using TensorExpr = Ddg_times_Tensor4_2301_ijkl<A, B, T, U, Dim01, Dim23, Dim2,
                                            Dim3, i, j, k, l, m, n>;
  return Tensor4_Expr<TensorExpr, typename promote<T, U>::V, Dim01, Dim01, Dim2,
                      Dim3, i, j, k, l>(TensorExpr(a, b));
}

/* A(m,n,k,l)*S(i,j,m,n) */

template <class A, class B, class T, class U, int Dim01, int Dim23, int Dim2,
          int Dim3, char i, char j, char k, char l, char m, char n>
Tensor4_Expr<Ddg_times_Tensor4_2301_ijkl<A, B, T, U, Dim01, Dim23, Dim2, Dim3, i, j,
                                    k, l, m, n>,
             typename promote<T, U>::V, Dim01, Dim01, Dim2, Dim3, i, j, k, l>
operator*(const Tensor4_Expr<B, U, Dim23, Dim23, Dim2, Dim3, m, n, k, l> &b,
          const Ddg_Expr<A, T, Dim01, Dim23, i, j, m, n> &a) {
  using TensorExpr = Ddg_times_Tensor4_2301_ijkl<A, B, T, U, Dim01, Dim23, Dim2,
                                            Dim3, i, j, k, l, m, n>;
  return Tensor4_Expr<TensorExpr, typename promote<T, U>::V, Dim01, Dim01, Dim2,
                      Dim3, i, j, k, l>(TensorExpr(a, b));
}

/* S(m,n,k,l)*A(i,j,m,n) */
/* A(i,j,m,n)*S(m,n,k,l) */

template <class A, class B, class T, class U, int Dim01, int Dim23, int Dim0,
          int Dim1, char i, char j, char k, char l, char m, char n>
class Ddg_times_Tensor4_2323_klij {
  Ddg_Expr<A, T, Dim01, Dim23, m, n, k, l> iterA;
  Tensor4_Expr<B, U, Dim0, Dim1, Dim23, Dim23, i, j, m, n> iterB;

  template <int Current_Dim0, int Current_Dim1>
  inline typename promote<T, U>::V
  eval(const int N1, const int N2, const int N3, const int N4,
       const Number<Current_Dim0> &, const Number<Current_Dim1> &) const {
    return iterA(Current_Dim0 - 1, Current_Dim1 - 1, N3, N4) *
               iterB(N1, N2, Current_Dim0 - 1, Current_Dim1 - 1) +
           eval(N1, N2, N3, N4, Number<Current_Dim0 - 1>(),
                Number<Current_Dim1>());
  }
  template <int Current_Dim1>
  inline typename promote<T, U>::V
  eval(const int N1, const int N2, const int N3, const int N4,
       const Number<1> &, const Number<Current_Dim1> &) const {
    return iterA(0, Current_Dim1 - 1, N3, N4) *
               iterB(N1, N2, 0, Current_Dim1 - 1) +
           eval(N1, N2, N3, N4, Number<Dim23>(), Number<Current_Dim1 - 1>());
  }
  inline typename promote<T, U>::V eval(const int N1, const int N2,
                                        const int N3, const int N4,
                                        const Number<1> &,
                                        const Number<1> &) const {
    return iterA(0, 0, N3, N4) * iterB(N1, N2, 0, 0);
  }

public:
  Ddg_times_Tensor4_2323_klij(
      const Ddg_Expr<A, T, Dim01, Dim23, m, n, k, l> &a,
      const Tensor4_Expr<B, U, Dim0, Dim1, Dim23, Dim23, i, j, m, n> &b)
      : iterA(a), iterB(b) {}
  typename promote<T, U>::V operator()(const int N1, const int N2, const int N3,
                                       const int N4) const {
    return eval(N1, N2, N3, N4, Number<Dim01>(), Number<Dim01>());
  }
};

/* S(k,l,m,n)*A(i,j,m,n) */

template <class A, class B, class T, class U, int Dim01, int Dim23, int Dim0,
          int Dim1, char i, char j, char k, char l, char m, char n>
Tensor4_Expr<Ddg_times_Tensor4_2323_klij<A, B, T, U, Dim01, Dim23, Dim0, Dim1,
                                         i, j, k, l, m, n>,
             typename promote<T, U>::V, Dim0, Dim1, Dim23, Dim23, i, j, k, l>
operator*(const Ddg_Expr<A, T, Dim01, Dim23, m, n, k, l> &a,
          const Tensor4_Expr<B, U, Dim0, Dim1, Dim23, Dim23, i, j, m, n> &b) {
  using TensorExpr = Ddg_times_Tensor4_2323_klij<A, B, T, U, Dim01, Dim23, Dim0,
                                                 Dim1, i, j, k, l, m, n>;
  return Tensor4_Expr<TensorExpr, typename promote<T, U>::V, Dim0, Dim1, Dim23,
                      Dim23, i, j, k, l>(TensorExpr(a, b));
}

/* A(i,j,m,n)*S(k,l,m,n) */

template <class A, class B, class T, class U, int Dim01, int Dim23, int Dim0,
          int Dim1, char i, char j, char k, char l, char m, char n>
Tensor4_Expr<Ddg_times_Tensor4_2323_klij<A, B, T, U, Dim01, Dim23, Dim0, Dim1,
                                         i, j, k, l, m, n>,
             typename promote<T, U>::V, Dim0, Dim1, Dim23, Dim23, i, j, k, l>
operator*(const Tensor4_Expr<B, U, Dim0, Dim1, Dim23, Dim23, i, j, m, n> &b,
          const Ddg_Expr<A, T, Dim01, Dim23, m, n, k, l> &a) {
  using TensorExpr = Ddg_times_Tensor4_2323_klij<A, B, T, U, Dim01, Dim23, Dim0,
                                                 Dim1, i, j, k, l, m, n>;
  return Tensor4_Expr<TensorExpr, typename promote<T, U>::V, Dim0, Dim1, Dim23,
                      Dim23, i, j, k, l>(TensorExpr(a, b));
}

} // namespace FTensor