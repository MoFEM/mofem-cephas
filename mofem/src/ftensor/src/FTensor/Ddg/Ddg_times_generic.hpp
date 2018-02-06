/* Multiplies a Ddg with a generic, yielding a
   Ddg. */

#pragma once

namespace FTensor
{
  template<class A, class T, class U, int Dim01, int Dim23,
           char i, char j, char k, char l>
  class Ddg_times_generic
  {
    const Ddg_Expr<A,T,Dim01,Dim23,i,j,k,l> iterA;
    const U d;
  public:
    typename promote<T,U>::V operator()(const int N1, const int N2, const int N3,
                                        const int N4) const
    {
      return iterA(N1,N2,N3,N4)*d;
    }

    Ddg_times_generic(const Ddg_Expr<A,T,Dim01,Dim23,i,j,k,l> &a,
                              const U &d0): iterA(a), d(d0) {}
  };

  template<class A, class T, class U, int Dim01, int Dim23,
           char i, char j, char k, char l>
  inline const Ddg_Expr
  <const Ddg_times_generic<A,T,U,Dim01,Dim23,i,j,k,l>,typename promote<T,U>::V,Dim01,Dim23,i,j,k,l>
  operator*(const Ddg_Expr<A,T,Dim01,Dim23,i,j,k,l> &a, const U &d0)
  {
    typedef const Ddg_times_generic<A,T,U,Dim01,Dim23,i,j,k,l>
      TensorExpr;
    return Ddg_Expr<TensorExpr,typename promote<T,U>::V,Dim01,Dim23,i,j,k,l>
      (TensorExpr(a,d0));
  }

  template<class A, class T, class U, int Dim01, int Dim23,
           char i, char j, char k, char l>
  inline const Ddg_Expr
  <const Ddg_times_generic<A,T,U,Dim01,Dim23,i,j,k,l>,typename promote<T,U>::V,Dim01,Dim23,i,j,k,l>
  operator*(const U &d0, const Ddg_Expr<A,T,Dim01,Dim23,i,j,k,l> &a)
  {
    typedef const Ddg_times_generic<A,T,U,Dim01,Dim23,i,j,k,l>
      TensorExpr;
    return Ddg_Expr<TensorExpr,typename promote<T,U>::V,Dim01,Dim23,i,j,k,l>
      (TensorExpr(a,d0));
  }
}
