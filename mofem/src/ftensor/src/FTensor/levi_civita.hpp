/* Special Fully Antisimetric
 * Levi-Civita Tensor for Vector Product */

#pragma once

#include "LeviCivita.hpp"

namespace FTensor
{
  template<class T=int,char i,char j,int Dim0,int Dim1>
  inline constexpr
  typename std::enable_if<(Dim0 <= 2 && Dim1 <= 2),
          Tensor2_Expr<Levi_Civita<T>,T,Dim0,Dim1,i,j>>::type
  levi_civita(const Index<i,Dim0>&, const Index<j,Dim1>&)
  {
    return Levi_Civita<T>()(Index<i,Dim0>(),Index<j,Dim1>());
  }

  template<class T=int,char i,char j,char k,int Dim0,int Dim1,int Dim2>
  inline constexpr
  typename std::enable_if<(Dim0 <= 3 && Dim1 <= 3 && Dim2 <= 3),
          Tensor3_Expr<Levi_Civita<T>,T,Dim0,Dim1,Dim2,i,j,k>>::type
  levi_civita(const Index<i,Dim0>&, const Index<j,Dim1>&, const Index<k,Dim2>&)
  {
    return Levi_Civita<T>()(Index<i,Dim0>(),Index<j,Dim1>(),Index<k,Dim2>());
  }

  template<class T=int,char i,char j,char k,char l,int Dim0,int Dim1,int Dim2,int Dim3>
  inline constexpr
  typename std::enable_if<(Dim0 <= 4 && Dim1 <= 4 && Dim2 <= 4 && Dim3 <= 4),
          Tensor4_Expr<Levi_Civita<T>,T,Dim0,Dim1,Dim2,Dim3,i,j,k,l>>::type
  levi_civita(const Index <i, Dim0> &, const Index <j, Dim1> &, const Index <k, Dim2> &, const Index <l, Dim3> &)
  {
    return Levi_Civita<T>()(Index<i,Dim0>(),Index<j,Dim1>(),Index<k,Dim2>(),Index<l,Dim3>());
  }
}

