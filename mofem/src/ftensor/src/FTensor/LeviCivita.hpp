/* Special Fully Antisimetric
 * Levi-Civita Tensor for Vector Product */

#pragma once

namespace FTensor
{
  template <class T=int>
  class Levi_Civita
  {
  public:
    /* A Tensor2 Implementation */
    constexpr
    T operator()(const int N1, const int N2) const
    {
      return (N1==N2) ? T(0) : ((N1==0) ? T(1) : T(-1));
    }

    template<char i,char j,int Dim0,int Dim1>
    typename std::enable_if<(Dim0 <= 2 && Dim1 <= 2 ),
            Tensor2_Expr<Levi_Civita<T>,T,Dim0,Dim1,i,j>>::type
    operator()(const Index<i,Dim0>&, const Index<j,Dim1>&)
    {
      return Tensor2_Expr<Levi_Civita<T>,T,Dim0,Dim1,i,j> (*this);
    };

    /* A Tensor3 Implementation */
    constexpr
    T operator()(const int N1, const int N2, const int N3) const
    {
      return (N1==N2 || N1==N3 || N2==N3) ? T(0) : (((N1+1)%3==N2) ? T(1) : T(-1));
    }

    template<char i,char j,char k,int Dim0,int Dim1,int Dim2>
    typename std::enable_if<(Dim0 <= 3 && Dim1 <= 3 && Dim2 <= 3),
            Tensor3_Expr<Levi_Civita<T>,T,Dim0,Dim1,Dim2,i,j,k>>::type
    operator()(const Index<i,Dim0>&, const Index<j,Dim1>&, const Index<k,Dim2>&)
    {
      return Tensor3_Expr<Levi_Civita<T>,T,Dim0,Dim1,Dim2,i,j,k> (*this);
    };

    /* A Tensor4 Implementation */
    constexpr
    T operator()(const int N1, const int N2, const int N3, const int N4) const
    {
      return (N1==N2 || N1==N3 || N1==N4 || N2==N3 || N2==N4 || N3==N4) ? T(0) :
        ((N1+N2==1 || N1+N2==5) ? (((N2+N3)%4==3) ? T(1):T(-1)) :
         ((N1+N2==2 || N1+N2==4) ? (((N2+N3)%4==1) ? T(1):T(-1)) :
          (N1+N2==3 ? (((N2+N3)%4!=1) ? T(1):T(-1)) :
           T(0))));
    }

    template<char i,char j,char k,char l,int Dim0,int Dim1,int Dim2,int Dim3>
    typename std::enable_if<(Dim0 <= 4 && Dim1 <= 4 && Dim2 <= 4 && Dim3 <= 4),
            Tensor4_Expr<Levi_Civita<T>,T,Dim0,Dim1,Dim2,Dim3,i,j,k,l>>::type
    operator()(const Index<i,Dim0>&, const Index<j,Dim1>&, const Index<k,Dim2>&, const Index<l,Dim3>&)
    {
      return Tensor4_Expr<Levi_Civita<T>,T,Dim0,Dim1,Dim2,Dim3,i,j,k,l> (*this);
    };
  };
}

