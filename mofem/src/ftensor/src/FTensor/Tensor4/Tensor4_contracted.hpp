/* A helper class that automatically returns the internal contraction
   of a rank 4 tensor. */

#pragma once

namespace FTensor
{
  template<class A,class T,int Dim>
  class Tensor4_contracted_02
  {
    const A iterA;

    template<int Current_Dim>
    T eval(const Number<Current_Dim> &, const int N0,
           const Number<Current_Dim> &, const int N1) const
    {
      return iterA(Current_Dim-1,N0,Current_Dim-1,N1) +
             eval(Number<Current_Dim-1>(),N0,Number<Current_Dim-1>(),N1);
    }
    T eval(const Number<1> &, const int N0,
           const Number<1> &, const int N1) const
    {
      return iterA(0,N0,0,N1);
    }
  public:
    T operator()(const int N0, const int N1) const
    {
      return eval(Number<Dim>(), N0, Number<Dim>(), N1);
    }
    Tensor4_contracted_02(const A &a): iterA(a) {}
  };
}
