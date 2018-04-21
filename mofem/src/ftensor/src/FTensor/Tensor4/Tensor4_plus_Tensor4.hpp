/* Adds Tensor4+Tensor4 -> Tensor4 */

#pragma once

namespace FTensor
{
  /* A(i,j,k,l)+B(i,j,k,l) */

  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2,
            int Dim3, char i, char j, char k, char l>
  class Tensor4_plus_Tensor4
  {
    const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> iterA;
    const Tensor4_Expr<B, U, Dim0, Dim1, Dim2, Dim3, i, j, k, l> iterB;

  public:
    typename promote<T, U>::V
    operator()(const int N1, const int N2, const int N3, const int N4) const
    {
      return iterA(N1, N2, N3, N4) + iterB(N1, N2, N3, N4);
    }

    Tensor4_plus_Tensor4(
      const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a,
      const Tensor4_Expr<B, U, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &b)
        : iterA(a), iterB(b)
    {}
  };

  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2,
            int Dim3, char i, char j, char k, char l>
  inline Tensor4_Expr<
    const Tensor4_plus_Tensor4<A, B, T, U, Dim0, Dim1, Dim2, Dim3, i, j, k, l>,
    typename promote<T, U>::V, Dim0, Dim1, Dim2, Dim3, i, j, k, l>
  operator+(const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a,
            const Tensor4_Expr<B, U, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &b)
  {
    typedef const Tensor4_plus_Tensor4<A, B, T, U, Dim0, Dim1, Dim2, Dim3, i,
                                       j, k, l>
      TensorExpr;
    return Tensor4_Expr<TensorExpr, typename promote<T, U>::V, Dim0, Dim1,
                        Dim2, Dim3, i, j, k, l>(TensorExpr(a, b));
  }

  /* A(i,j,k,l)+B(i,l,j,k) */

  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2,
            int Dim3, char i, char j, char k, char l>
  class Tensor4_plus_Tensor4_0312
  {
    const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> iterA;
    const Tensor4_Expr<B, U, Dim0, Dim3, Dim1, Dim2, i, l, j, k> iterB;

  public:
    typename promote<T, U>::V
    operator()(const int N1, const int N2, const int N3, const int N4) const
    {
      return iterA(N1, N2, N3, N4) + iterB(N1, N4, N2, N3);
    }

    Tensor4_plus_Tensor4_0312(
      const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a,
      const Tensor4_Expr<B, U, Dim0, Dim3, Dim1, Dim2, i, l, j, k> &b)
        : iterA(a), iterB(b)
    {}
  };

  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2,
            int Dim3, char i, char j, char k, char l>
  inline Tensor4_Expr<const Tensor4_plus_Tensor4_0312<A, B, T, U, Dim0, Dim1,
                                                      Dim2, Dim3, i, j, k, l>,
                      typename promote<T, U>::V, Dim0, Dim1, Dim2, Dim3, i, j,
                      k, l>
  operator+(const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a,
            const Tensor4_Expr<B, U, Dim0, Dim3, Dim1, Dim2, i, l, j, k> &b)
  {
    typedef const Tensor4_plus_Tensor4_0312<A, B, T, U, Dim0, Dim1, Dim2, Dim3,
                                            i, j, k, l>
      TensorExpr;
    return Tensor4_Expr<TensorExpr, typename promote<T, U>::V, Dim0, Dim1,
                        Dim2, Dim3, i, j, k, l>(TensorExpr(a, b));
  }

  /* A(i,j,k,l)+B(i,l,k,j) */

  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2,
            int Dim3, char i, char j, char k, char l>
  class Tensor4_plus_Tensor4_0321
  {
    const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> iterA;
    const Tensor4_Expr<B, U, Dim0, Dim3, Dim2, Dim1, i, l, k, j> iterB;

  public:
    typename promote<T, U>::V
    operator()(const int N1, const int N2, const int N3, const int N4) const
    {
      return iterA(N1, N2, N3, N4) + iterB(N1, N4, N3, N2);
    }

    Tensor4_plus_Tensor4_0321(
      const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a,
      const Tensor4_Expr<B, U, Dim0, Dim3, Dim2, Dim1, i, l, k, j> &b)
        : iterA(a), iterB(b)
    {}
  };

  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2,
            int Dim3, char i, char j, char k, char l>
  inline Tensor4_Expr<const Tensor4_plus_Tensor4_0321<A, B, T, U, Dim0, Dim1,
                                                      Dim2, Dim3, i, j, k, l>,
                      typename promote<T, U>::V, Dim0, Dim1, Dim2, Dim3, i, j,
                      k, l>
  operator+(const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a,
            const Tensor4_Expr<B, U, Dim0, Dim3, Dim2, Dim1, i, l, k, j> &b)
  {
    typedef const Tensor4_plus_Tensor4_0321<A, B, T, U, Dim0, Dim1, Dim2, Dim3,
                                            i, j, k, l>
      TensorExpr;
    return Tensor4_Expr<TensorExpr, typename promote<T, U>::V, Dim0, Dim1,
                        Dim2, Dim3, i, j, k, l>(TensorExpr(a, b));
  }

  /* A(i,j,k,l)+B(i,k,j,l) */

  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2,
            int Dim3, char i, char j, char k, char l>
  class Tensor4_plus_Tensor4_0213
  {
    const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> iterA;
    const Tensor4_Expr<B, U, Dim0, Dim2, Dim1, Dim3, i, k, j, l> iterB;

  public:
    typename promote<T, U>::V
    operator()(const int N1, const int N2, const int N3, const int N4) const
    {
      return iterA(N1, N2, N3, N4) + iterB(N1, N3, N2, N4);
    }

    Tensor4_plus_Tensor4_0213(
      const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a,
      const Tensor4_Expr<B, U, Dim0, Dim2, Dim1, Dim3, i, k, j, l> &b)
        : iterA(a), iterB(b)
    {}
  };

  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2,
            int Dim3, char i, char j, char k, char l>
  inline Tensor4_Expr<const Tensor4_plus_Tensor4_0213<A, B, T, U, Dim0, Dim1,
                                                      Dim2, Dim3, i, j, k, l>,
                      typename promote<T, U>::V, Dim0, Dim1, Dim2, Dim3, i, j,
                      k, l>
  operator+(const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a,
            const Tensor4_Expr<B, U, Dim0, Dim2, Dim1, Dim3, i, k, j, l> &b)
  {
    typedef const Tensor4_plus_Tensor4_0213<A, B, T, U, Dim0, Dim1, Dim2, Dim3,
                                            i, j, k, l>
      TensorExpr;
    return Tensor4_Expr<TensorExpr, typename promote<T, U>::V, Dim0, Dim1,
                        Dim2, Dim3, i, j, k, l>(TensorExpr(a, b));
  }

  /* A(i,j,k,l)+B(l,k,i,j) */

  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2,
            int Dim3, char i, char j, char k, char l>
  class Tensor4_plus_Tensor4_3201
  {
    const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> iterA;
    const Tensor4_Expr<B, U, Dim3, Dim2, Dim0, Dim1, l, k, i, j> iterB;

  public:
    typename promote<T, U>::V
    operator()(const int N1, const int N2, const int N3, const int N4) const
    {
      return iterA(N1, N2, N3, N4) + iterB(N4, N3, N1, N2);
    }

    Tensor4_plus_Tensor4_3201(
      const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a,
      const Tensor4_Expr<B, U, Dim3, Dim2, Dim0, Dim1, l, k, i, j> &b)
        : iterA(a), iterB(b)
    {}
  };

  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2,
            int Dim3, char i, char j, char k, char l>
  inline Tensor4_Expr<const Tensor4_plus_Tensor4_3201<A, B, T, U, Dim0, Dim1,
                                                      Dim2, Dim3, i, j, k, l>,
                      typename promote<T, U>::V, Dim0, Dim1, Dim2, Dim3, i, j,
                      k, l>
  operator+(const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a,
            const Tensor4_Expr<B, U, Dim3, Dim2, Dim0, Dim1, l, k, i, j> &b)
  {
    typedef const Tensor4_plus_Tensor4_3201<A, B, T, U, Dim0, Dim1, Dim2, Dim3,
                                            i, j, k, l>
      TensorExpr;
    return Tensor4_Expr<TensorExpr, typename promote<T, U>::V, Dim0, Dim1,
                        Dim2, Dim3, i, j, k, l>(TensorExpr(a, b));
  }

  /* A(i,j,k,l)+B(i,j,l,k) */

  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2,
            int Dim3, char i, char j, char k, char l>
  class Tensor4_plus_Tensor4_0132
  {
    const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> iterA;
    const Tensor4_Expr<B, U, Dim0, Dim1, Dim3, Dim2, i, j, l, k> iterB;

  public:
    typename promote<T, U>::V
    operator()(const int N1, const int N2, const int N3, const int N4) const
    {
      return iterA(N1, N2, N3, N4) + iterB(N1, N2, N4, N3);
    }

    Tensor4_plus_Tensor4_0132(
      const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a,
      const Tensor4_Expr<B, U, Dim0, Dim1, Dim3, Dim2, i, j, l, k> &b)
        : iterA(a), iterB(b)
    {}
  };

  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2,
            int Dim3, char i, char j, char k, char l>
  inline Tensor4_Expr<const Tensor4_plus_Tensor4_0132<A, B, T, U, Dim0, Dim1,
                                                      Dim2, Dim3, i, j, k, l>,
                      typename promote<T, U>::V, Dim0, Dim1, Dim2, Dim3, i, j,
                      k, l>
  operator+(const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a,
            const Tensor4_Expr<B, U, Dim0, Dim1, Dim3, Dim2, i, j, l, k> &b)
  {
    typedef const Tensor4_plus_Tensor4_0132<A, B, T, U, Dim0, Dim1, Dim2, Dim3,
                                            i, j, k, l>
      TensorExpr;
    return Tensor4_Expr<TensorExpr, typename promote<T, U>::V, Dim0, Dim1,
                        Dim2, Dim3, i, j, k, l>(TensorExpr(a, b));
  }
}
