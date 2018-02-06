/* Adds Ddg+Ddg -> Ddg */

#pragma once

namespace FTensor
{
  /* A(i,j,k,l)+B(i,j,k,l) -> Ddg */

  template<class A, class B, class T, class U, int Dim01, int Dim23,
           char i, char j, char k, char l>
  class Ddg_plus_Ddg
  {
    const Ddg_Expr<A,T,Dim01,Dim23,i,j,k,l> iterA;
    const Ddg_Expr<B,U,Dim01,Dim23,i,j,k,l> iterB;
  public:
    typename promote<T,U>::V operator()(const int N1, const int N2, const int N3,
                                        const int N4) const
    {
      return iterA(N1,N2,N3,N4)+iterB(N1,N2,N3,N4);
    }

    Ddg_plus_Ddg(const Ddg_Expr<A,T,Dim01,Dim23,i,j,k,l> &a, const Ddg_Expr<B,U,Dim01,Dim23,i,j,k,l> &b):
      iterA(a), iterB(b) {}
  };

  template<class A, class B, class T, class U, int Dim01, int Dim23,
           char i, char j, char k, char l>
  inline const Ddg_Expr
  <const Ddg_plus_Ddg<A,B,T,U,Dim01,Dim23,i,j,k,l>,
   typename promote<T,U>::V,Dim01,Dim23,i,j,k,l>
  operator+(const Ddg_Expr<A,T,Dim01,Dim23,i,j,k,l> &a,
            const Ddg_Expr<B,U,Dim01,Dim23,i,j,k,l> &b)
  {
    typedef const Ddg_plus_Ddg<A,B,T,U,Dim01,Dim23,i,j,k,l>
      TensorExpr;
    return Ddg_Expr<TensorExpr,typename promote<T,U>::V,Dim01,Dim23,i,j,k,l>
      (TensorExpr(a,b));
  }

  /* A(i,j,k,l)+B(k,l,i,j) -> Ddg */

  template<class A, class B, class T, class U, int Dim01, int Dim23,
           char i, char j, char k, char l>
  class Ddg_plus_Ddg_2301
  {
    const Ddg_Expr<A,T,Dim01,Dim23,i,j,k,l> iterA;
    const Ddg_Expr<B,U,Dim23,Dim01,k,l,i,j> iterB;
  public:
    typename promote<T,U>::V operator()(const int N1, const int N2, const int N3,
                                        const int N4) const
    {
      return iterA(N1,N2,N3,N4)+iterB(N3,N4,N1,N2);
    }

    Ddg_plus_Ddg_2301
    (const Ddg_Expr<A,T,Dim01,Dim23,i,j,k,l> &a,
     const Ddg_Expr<B,U,Dim23,Dim01,k,l,i,j> &b): iterA(a), iterB(b) {}
  };

  template<class A, class B, class T, class U, int Dim01, int Dim23,
           char i, char j, char k, char l>
  inline const Ddg_Expr
  <const Ddg_plus_Ddg_2301<A,B,T,U,Dim01,Dim23,i,j,k,l>,
   typename promote<T,U>::V,Dim01,Dim23,i,j,k,l>
  operator+(const Ddg_Expr<A,T,Dim01,Dim23,i,j,k,l> &a,
            const Ddg_Expr<B,U,Dim23,Dim01,k,l,i,j> &b)
  {
    typedef const Ddg_plus_Ddg_2301<A,B,T,U,Dim01,Dim23,i,j,k,l>
      TensorExpr;
    return Ddg_Expr<TensorExpr,typename promote<T,U>::V,Dim01,Dim23,i,j,k,l>
      (TensorExpr(a,b));
  }
}
