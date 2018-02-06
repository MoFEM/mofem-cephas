/* Declare a wrapper class for generic rank 3 Tensor expressions.
   There isn't a Tensor3 class yet.  I only use Tensor3_Expr as an
   intermediate expression which immediately get contracted with
   something to make a Tensor2 or Tensor1. */

#pragma once

#include "Tensor3_times_generic.hpp"
#include "Tensor3_times_Tensor1.hpp"
#include "Tensor3_times_Tensor2_symmetric.hpp"
#include "Tensor3_times_Tensor2.hpp"
#include "Tensor3_times_Tensor3.hpp"
#include "Tensor3_times_Dg.hpp"
#include "Tensor3_plus_Tensor3.hpp"
#include "Tensor3_or_Tensor3.hpp"
#include "Tensor3_minus_Dg.hpp"

namespace FTensor
{
  template<class A, class T, int Dim0, int Dim1, int Dim2,
           char i, char j, char k>
  class Tensor3_Expr
  {
    A iter;
  public:
    Tensor3_Expr(A &a): iter(a) {}
    T operator()(const int N1, const int N2, const int N3) const
    {
      return iter(N1,N2,N3);
    }
  };

  template<class A, class T, int Tensor_Dim0, int Tensor_Dim1, int Tensor_Dim2,
           int Dim0, int Dim1, int Dim2, char i, char j, char k>
  class Tensor3_Expr<Tensor3<A,Tensor_Dim0,Tensor_Dim1,Tensor_Dim2>,T,
                     Dim0,Dim1,Dim2,i,j,k>
  {
    Tensor3<A,Tensor_Dim0,Tensor_Dim1,Tensor_Dim2> &iter;
  public:
    Tensor3_Expr(Tensor3<A,Tensor_Dim0,Tensor_Dim1,Tensor_Dim2> &a): iter(a) {}
    T & operator()(const int N1, const int N2, const int N3)
    {
      return iter(N1,N2,N3);
    }
    T operator()(const int N1, const int N2, const int N3) const
    {
      return iter(N1,N2,N3);
    }

    /* Various assignment operators.  I have to explicitly declare the
       second operator= because otherwise the compiler will generate its
       own and not use the template code. */

    template<class B, class U> inline
    const Tensor3_Expr<Tensor3<A,Tensor_Dim0,Tensor_Dim1,Tensor_Dim2>,T,
                       Dim0,Dim1,Dim2,i,j,k> &
    operator=(const Tensor3_Expr<B,U,Dim0,Dim1,Dim2,i,j,k> &result);

    const Tensor3_Expr<Tensor3<A,Tensor_Dim0,Tensor_Dim1,Tensor_Dim2>,T,
                       Dim0,Dim1,Dim2,i,j,k> &
    operator=(const Tensor3_Expr<Tensor3<A,Tensor_Dim0,Tensor_Dim1,Tensor_Dim2>,T,
              Dim0,Dim1,Dim2,i,j,k> &result);

    template<class B, class U> inline
    const Tensor3_Expr<Tensor3<A,Tensor_Dim0,Tensor_Dim1,Tensor_Dim2>,T,
                       Dim0,Dim1,Dim2,i,j,k> &
    operator=(const Tensor3_Expr<B,U,Dim0,Dim2,Dim1,i,k,j> &result);

    template<class B, class U> inline
    const Tensor3_Expr<Tensor3<A,Tensor_Dim0,Tensor_Dim1,Tensor_Dim2>,T,
                       Dim0,Dim1,Dim2,i,j,k> &
    operator=(const Tensor3_Expr<B,U,Dim1,Dim0,Dim2,j,i,k> &result);

    template<class B, class U> inline
    const Tensor3_Expr<Tensor3<A,Tensor_Dim0,Tensor_Dim1,Tensor_Dim2>,T,
                       Dim0,Dim1,Dim2,i,j,k> &
    operator=(const Tensor3_Expr<B,U,Dim1,Dim2,Dim0,j,k,i> &result);

    template<class B, class U> inline
    const Tensor3_Expr<Tensor3<A,Tensor_Dim0,Tensor_Dim1,Tensor_Dim2>,T,
                       Dim0,Dim1,Dim2,i,j,k> &
    operator=(const Tensor3_Expr<B,U,Dim2,Dim0,Dim1,k,i,j> &result);

    template<class B, class U> inline
    const Tensor3_Expr<Tensor3<A,Tensor_Dim0,Tensor_Dim1,Tensor_Dim2>,T,
                       Dim0,Dim1,Dim2,i,j,k> &
    operator=(const Tensor3_Expr<B,U,Dim2,Dim1,Dim0,k,j,i> &result);

    template <class U> inline
    const Tensor3_Expr<Tensor3<A,Tensor_Dim0,Tensor_Dim1,Tensor_Dim2>,T,
                       Dim0,Dim1,Dim2,i,j,k> &
    operator=(const U &d);
  };
}

#include "Tensor3_Expr_equals.hpp"
