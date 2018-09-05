/* This file has all of the declarations for expressions like
   Tensor4*Tensor3 and Tensor3*Tensor4, yielding a
   Tensor3. */

#pragma once

namespace FTensor
{
  /* A(i,j,k,l)*B(k,l,m) */

  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2,
            int Dim3, int Dim4, char i, char j, char k, char l, char m>
  class Tensor4_times_Tensor3_23
  {
    const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> iterA;
    const Tensor3_Expr<B, U, Dim2, Dim3, Dim4, k, l, m> iterB;

    template <int Current_Dim0, int Current_Dim1>
    typename promote<T, U>::V
    eval(const int N1, const int N2, const int N3, const Number<Current_Dim0> &,
         const Number<Current_Dim1> &) const
    {
      return iterA(N1, N2, Current_Dim0 - 1, Current_Dim1 - 1)
               * iterB(Current_Dim0 - 1, Current_Dim1 - 1, N3)
             + eval(N1, N2, N3, Number<Current_Dim0 - 1>(),
                    Number<Current_Dim1>());
    }
    template <int Current_Dim1>
    typename promote<T, U>::V
    eval(const int N1, const int N2, const int N3, const Number<1> &,
         const Number<Current_Dim1> &) const
    {
      return iterA(N1, N2, 0, Current_Dim1 - 1) * iterB(0, Current_Dim1 - 1, N3)
             + eval(N1, N2, N3, Number<Dim2>(), Number<Current_Dim1 - 1>());
    }
    typename promote<T, U>::V eval(const int N1, const int N2, const int N3,
                                   const Number<1> &, const Number<1> &) const
    {
      return iterA(N1, N2, 0, 0) * iterB(0, 0, N3);
    }

  public:
    typename promote<T, U>::V operator()(const int N1, const int N2,
                                         const int N3) const {
      return eval(N1, N2, N3, Number<Dim2>(), Number<Dim3>());
    }

    Tensor4_times_Tensor3_23(
      const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a,
      const Tensor3_Expr<B, U, Dim2, Dim3, Dim4, k, l, m> &b)
        : iterA(a), iterB(b)
    {}
  };

  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2,
            int Dim3, int Dim4, char i, char j, char k, char l, char m>
  inline Tensor3_Expr<
      const Tensor4_times_Tensor3_23<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4,
                                     i, j, k, l, m>,
      typename promote<T, U>::V, Dim0, Dim1, Dim4, i, j, m>
  operator*(const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a,
            const Tensor3_Expr<B, U, Dim2, Dim3, Dim4, k, l, m> &b) {
    typedef const Tensor4_times_Tensor3_23<A, B, T, U, Dim0, Dim1, Dim2, Dim3,
                                           Dim4, i, j, k, l, m>
        TensorExpr;
    return Tensor3_Expr<TensorExpr, typename promote<T, U>::V, Dim0, Dim1, Dim4,
                        i, j, m>(TensorExpr(a, b));
  }

  /* B(k,l, m)*A(i,j,k,l) */

  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2,
            int Dim3, int Dim4, char i, char j, char k, char l, char m>
  inline Tensor3_Expr<
      const Tensor4_times_Tensor3_23<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4,
                                     i, j, k, l, m>,
      typename promote<T, U>::V, Dim0, Dim1, Dim4, i, j, m>
  operator*(const Tensor3_Expr<B, U, Dim2, Dim3, Dim4, k, l, m> &b,
            const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a) {
    typedef const Tensor4_times_Tensor3_23<A, B, T, U, Dim0, Dim1, Dim2, Dim3,
                                           Dim4, i, j, k, l, m>
        TensorExpr;
    return Tensor3_Expr<TensorExpr, typename promote<T, U>::V, Dim0, Dim1, Dim4,
                        i, j, m>(TensorExpr(a, b));
  }

} // namespace FTensor