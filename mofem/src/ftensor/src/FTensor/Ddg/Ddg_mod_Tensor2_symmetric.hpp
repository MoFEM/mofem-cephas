/* Divide a Ddg by a Tensor2_symmetric without contracting, yielding a
   Ddg. */

#pragma once

namespace FTensor
{
  /* A(i,j,k,l) % B(i,j) -> Ddg */

  template <class A, class B, class T, class U, int Dim01, int Dim23, char i,
            char j, char k, char l>
  class Ddg_mod_Tensor2_symmetric_01
  {
    const Ddg_Expr<A, T, Dim01, Dim23, i, j, k, l> iterA;
    const Tensor2_symmetric_Expr<B, U, Dim01, i, j> iterB;

  public:
    typename promote<T, U>::V
    operator()(const int N1, const int N2, const int N3, const int N4) const
    {
      return iterA(N1, N2, N3, N4) / iterB(N1, N2);
    }

    Ddg_mod_Tensor2_symmetric_01(
      const Ddg_Expr<A, T, Dim01, Dim23, i, j, k, l> &a,
      const Tensor2_symmetric_Expr<B, U, Dim01, i, j> &b)
        : iterA(a), iterB(b)
    {}
  };

  template <class A, class B, class T, class U, int Dim01, int Dim23, char i,
            char j, char k, char l>
  inline const Ddg_Expr<
    const Ddg_mod_Tensor2_symmetric_01<A, B, T, U, Dim01, Dim23, i, j, k, l>,
    typename promote<T, U>::V, Dim01, Dim23, i, j, k, l>
  operator%(const Ddg_Expr<A, T, Dim01, Dim23, i, j, k, l> &a,
            const Tensor2_symmetric_Expr<B, U, Dim01, i, j> &b)
  {
    typedef const Ddg_mod_Tensor2_symmetric_01<A, B, T, U, Dim01, Dim23, i, j,
                                               k, l>
      TensorExpr;
    return Ddg_Expr<TensorExpr, typename promote<T, U>::V, Dim01, Dim23, i, j,
                    k, l>(TensorExpr(a, b));
  }

  /* B(i,j) % A(i,j,k,l) -> Ddg */

  template <class A, class B, class T, class U, int Dim01, int Dim23, char i,
            char j, char k, char l>
  inline const Ddg_Expr<
    const Ddg_mod_Tensor2_symmetric_01<A, B, T, U, Dim01, Dim23, i, j, k, l>,
    typename promote<T, U>::V, Dim01, Dim23, i, j, k, l>
  operator%(const Tensor2_symmetric_Expr<B, U, Dim01, i, j> &b,
            const Ddg_Expr<A, T, Dim01, Dim23, i, j, k, l> &a)
  {
    typedef const Ddg_mod_Tensor2_symmetric_01<A, B, T, U, Dim01, Dim23, i, j,
                                               k, l>
      TensorExpr;
    return Ddg_Expr<TensorExpr, typename promote<T, U>::V, Dim01, Dim23, i, j,
                    k, l>(TensorExpr(a, b));
  }

  /* A(i,j,k,l) % B(k,l) -> Ddg */

  template <class A, class B, class T, class U, int Dim01, int Dim23, char i,
            char j, char k, char l>
  class Ddg_mod_Tensor2_symmetric_23
  {
    const Ddg_Expr<A, T, Dim01, Dim23, i, j, k, l> iterA;
    const Tensor2_symmetric_Expr<B, U, Dim01, k, l> iterB;

  public:
    typename promote<T, U>::V
    operator()(const int N1, const int N2, const int N3, const int N4) const
    {
      return iterA(N1, N2, N3, N4) / iterB(N3, N4);
    }

    Ddg_mod_Tensor2_symmetric_23(
      const Ddg_Expr<A, T, Dim01, Dim23, i, j, k, l> &a,
      const Tensor2_symmetric_Expr<B, U, Dim01, k, l> &b)
        : iterA(a), iterB(b)
    {}
  };

  template <class A, class B, class T, class U, int Dim01, int Dim23, char i,
            char j, char k, char l>
  inline const Ddg_Expr<
    const Ddg_mod_Tensor2_symmetric_23<A, B, T, U, Dim01, Dim23, i, j, k, l>,
    typename promote<T, U>::V, Dim01, Dim23, i, j, k, l>
  operator%(const Ddg_Expr<A, T, Dim01, Dim23, i, j, k, l> &a,
            const Tensor2_symmetric_Expr<B, U, Dim01, k, l> &b)
  {
    typedef const Ddg_mod_Tensor2_symmetric_23<A, B, T, U, Dim01, Dim23, i, j,
                                               k, l>
      TensorExpr;
    return Ddg_Expr<TensorExpr, typename promote<T, U>::V, Dim01, Dim23, i, j,
                    k, l>(TensorExpr(a, b));
  }

  /* B(k,l) % A(i,j,k,l) -> Ddg */

  template <class A, class B, class T, class U, int Dim01, int Dim23, char i,
            char j, char k, char l>
  inline const Ddg_Expr<
    const Ddg_mod_Tensor2_symmetric_23<A, B, T, U, Dim01, Dim23, i, j, k, l>,
    typename promote<T, U>::V, Dim01, Dim23, i, j, k, l>
  operator%(const Tensor2_symmetric_Expr<B, U, Dim01, k, l> &b,
            const Ddg_Expr<A, T, Dim01, Dim23, i, j, k, l> &a)
  {
    typedef const Ddg_mod_Tensor2_symmetric_23<A, B, T, U, Dim01, Dim23, i, j,
                                               k, l>
      TensorExpr;
    return Ddg_Expr<TensorExpr, typename promote<T, U>::V, Dim01, Dim23, i, j,
                    k, l>(TensorExpr(a, b));
  }
}
