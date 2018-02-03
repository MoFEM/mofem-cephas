/* A general version, not for pointers. */

#include <ostream>

#ifdef FTENSOR_DEBUG
#include <sstream>
#include <stdexcept>
#endif

#pragma once

namespace FTensor
{
  template<class T,int Tensor_Dim0,int Tensor_Dim1,int Tensor_Dim2,int Tensor_Dim3>
  class Tensor4
  {
    T data[Tensor_Dim0][Tensor_Dim1][Tensor_Dim2][Tensor_Dim3];
  public:
    /* Initialization operators
     * TODO revisit this operators*/
    Tensor4(T d0000, T d0001, T d0002, T d0010, T d0011, T d0012, T d0020, T d0021, T d0022,
            T d0100, T d0101, T d0102, T d0110, T d0111, T d0112, T d0120, T d0121, T d0122,
            T d0200, T d0201, T d0202, T d0210, T d0211, T d0212, T d0220, T d0221, T d0222,
            T d1000, T d1001, T d1002, T d1010, T d1011, T d1012, T d1020, T d1021, T d1022,
            T d1100, T d1101, T d1102, T d1110, T d1111, T d1112, T d1120, T d1121, T d1122,
            T d1200, T d1201, T d1202, T d1210, T d1211, T d1212, T d1220, T d1221, T d1222,
            T d2000, T d2001, T d2002, T d2010, T d2011, T d2012, T d2020, T d2021, T d2022,
            T d2100, T d2101, T d2102, T d2110, T d2111, T d2112, T d2120, T d2121, T d2122,
            T d2200, T d2201, T d2202, T d2210, T d2211, T d2212, T d2220, T d2221, T d2222)
    {
      data[0][0][0][0]=d0000;data[0][0][0][1]=d0001;data[0][0][0][2]=d0002;data[0][0][1][0]=d0010;data[0][0][1][1]=d0011;data[0][0][1][2]=d0012;data[0][0][2][0]=d0020;data[0][0][2][1]=d0021;data[0][0][2][2]=d0022;
      data[0][1][0][0]=d0100;data[0][1][0][1]=d0101;data[0][1][0][2]=d0102;data[0][1][1][0]=d0110;data[0][1][1][1]=d0111;data[0][1][1][2]=d0112;data[0][1][2][0]=d0120;data[0][1][2][1]=d0121;data[0][1][2][2]=d0122;
      data[0][2][0][0]=d0200;data[0][2][0][1]=d0201;data[0][2][0][2]=d0202;data[0][2][1][0]=d0210;data[0][2][1][1]=d0211;data[0][2][1][2]=d0212;data[0][2][2][0]=d0220;data[0][2][2][1]=d0221;data[0][2][2][2]=d0222;
      data[1][0][0][0]=d1000;data[1][0][0][1]=d1001;data[1][0][0][2]=d1002;data[1][0][1][0]=d1010;data[1][0][1][1]=d1011;data[1][0][1][2]=d1012;data[1][0][2][0]=d1020;data[1][0][2][1]=d1021;data[1][0][2][2]=d1022;
      data[1][1][0][0]=d1100;data[1][1][0][1]=d1101;data[1][1][0][2]=d1102;data[1][1][1][0]=d1110;data[1][1][1][1]=d1111;data[1][1][1][2]=d1112;data[1][1][2][0]=d1120;data[1][1][2][1]=d1121;data[1][1][2][2]=d1122;
      data[1][2][0][0]=d1200;data[1][2][0][1]=d1201;data[1][2][0][2]=d1202;data[1][2][1][0]=d1210;data[1][2][1][1]=d1211;data[1][2][1][2]=d1212;data[1][2][2][0]=d1220;data[1][2][2][1]=d1221;data[1][2][2][2]=d1222;
      data[2][0][0][0]=d2000;data[2][0][0][1]=d2001;data[2][0][0][2]=d2002;data[2][0][1][0]=d2010;data[2][0][1][1]=d2011;data[2][0][1][2]=d2012;data[2][0][2][0]=d2020;data[2][0][2][1]=d2021;data[2][0][2][2]=d2022;
      data[2][1][0][0]=d2100;data[2][1][0][1]=d2101;data[2][1][0][2]=d2102;data[2][1][1][0]=d2110;data[2][1][1][1]=d2111;data[2][1][1][2]=d2112;data[2][1][2][0]=d2120;data[2][1][2][1]=d2121;data[2][1][2][2]=d2122;
      data[2][2][0][0]=d2200;data[2][2][0][1]=d2201;data[2][2][0][2]=d2202;data[2][2][1][0]=d2210;data[2][2][1][1]=d2211;data[2][2][1][2]=d2212;data[2][2][2][0]=d2220;data[2][2][2][1]=d2221;data[2][2][2][2]=d2222;
    }
    Tensor4() {}
    /* There are two operator(int,int,int,int)'s, one for non-consts
       that lets you change the value, and one for consts that
       doesn't. */
    T & operator()(const int N1,const int N2,const int N3,const int N4)
    {
#ifdef FTENSOR_DEBUG
      if(N1>=Tensor_Dim0 || N1<0 || N2>=Tensor_Dim1 || N2<0
         || N3>=Tensor_Dim2 || N3<0 || N4>=Tensor_Dim3 || N4<0)
        {
          std::stringstream s;
          s << "Bad index in Tensor4<T,"
            << Tensor_Dim0 << "," << Tensor_Dim1 << "," << Tensor_Dim2 << "," << Tensor_Dim3
            << ">.operator("
            << N1 << "," << N2 << "," << N3 << "," << N4
            << ") const"
            << std::endl;
          throw std::runtime_error(s.str());
        }
#endif
      return data[N1][N2][N3][N4];
    }

    T operator()(const int N1,const int N2,const int N3,const int N4)
            const
    {
#ifdef FTENSOR_DEBUG
      if(N1>=Tensor_Dim0 || N1<0 || N2>=Tensor_Dim1 || N2<0
         || N3>=Tensor_Dim2 || N3<0 || N4>=Tensor_Dim3 || N4<0)
        {
          std::stringstream s;
          s << "Bad index in Tensor4<T,"
            << Tensor_Dim0 << "," << Tensor_Dim1 << "," << Tensor_Dim2 << "," << Tensor_Dim3
            << ">.operator("
            << N1 << "," << N2 << "," << N3 << "," << N4 << ")"
            << std::endl;
          throw std::runtime_error(s.str());
        }
#endif
      return data[N1][N2][N3][N4];
    }

    /* These operator()'s are the first part in constructing template
       expressions.  They can be used to slice off lower dimensional
       parts. */

    template<char i,char j,char k,char l,int Dim0,int Dim1,int Dim2,int Dim3>
    typename std::enable_if<(Tensor_Dim0 >= Dim0 && Tensor_Dim1 >= Dim1 && Tensor_Dim2 >= Dim2 && Tensor_Dim3 >= Dim3),
            Tensor4_Expr<Tensor4<T,Tensor_Dim0,Tensor_Dim1,Tensor_Dim2,Tensor_Dim3>,
                   T,Dim0,Dim1,Dim2,Dim3,i,j,k,l> >::type
    operator()(const Index<i,Dim0> , const Index<j,Dim1> ,
               const Index<k,Dim2> , const Index<l,Dim3>)
    {
      return Tensor4_Expr<Tensor4<T,Tensor_Dim0,Tensor_Dim1,Tensor_Dim2,Tensor_Dim3>,
                          T,Dim0,Dim1,Dim2,Dim3,i,j,k,l> (*this);
    };

    template<char i,char j,char k,char l,int Dim0,int Dim1,int Dim2,int Dim3>
    typename std::enable_if<(Tensor_Dim0 >= Dim0 && Tensor_Dim1 >= Dim1 && Tensor_Dim2 >= Dim2 && Tensor_Dim3 >= Dim3),
            Tensor4_Expr<const Tensor4<T,Tensor_Dim0,Tensor_Dim1,Tensor_Dim2,Tensor_Dim3>,
            T,Dim0,Dim1,Dim2,Dim3,i,j,k,l> >::type
    operator()(const Index<i,Dim0> , const Index<j,Dim1> ,
               const Index<k,Dim2> , const Index<l,Dim3>) const
    {
      return Tensor4_Expr<const Tensor4<T,Tensor_Dim0,Tensor_Dim1,Tensor_Dim2,Tensor_Dim3>,
                          T,Dim0,Dim1,Dim2,Dim3,i,j,k,l> (*this);
    };

    /* These operators are for internal contractions, resulting in a Tensor2.
     * For example something like A(k,i,k,j) */


    template<char i,char j,char k,int Dim02,int Dim1,int Dim3>
    inline
    typename std::enable_if<(Tensor_Dim0 >= Dim02 && Tensor_Dim1 >= Dim1 && Tensor_Dim2 >= Dim02 && Tensor_Dim3 >= Dim3),
            Tensor2_Expr<const Tensor4_contracted_02<Tensor4<T,Tensor_Dim0,Tensor_Dim1,Tensor_Dim2,Tensor_Dim3>,
            T,Dim02>,T,Dim1,Dim3,i,j> >::type
    operator()(const Index<k,Dim02> , const Index<i,Dim1> ,
               const Index<k,Dim02> , const Index<j,Dim3>) const
    {
      typedef const Tensor4_contracted_02<Tensor4<T,Tensor_Dim0,Tensor_Dim1,Tensor_Dim2,Tensor_Dim3>,T,Dim02> TensorExpr;
      return Tensor2_Expr<TensorExpr,T,Dim1,Dim3,i,j>(TensorExpr(*this));
    };


    /* This is for expressions where a number is used for one slot, and
       an index for another, yielding a Tensor3_Expr.  The non-const
       versions don't actually create a Tensor3_number_rhs_[0123] object.
       They create a Tensor3_Expr directly, which provides the
       appropriate indexing operators.  The const versions do create a
       Tensor3_number_[0123]. */
    // TODO

  };
}

/// JSON compatible output

namespace FTensor
{
  template<class T,int Tensor_Dim0,int Tensor_Dim1,int Tensor_Dim2,int Tensor_Dim3>
  std::ostream & Tensor4_0001(std::ostream &os,
                              const Tensor4<T,Tensor_Dim0,Tensor_Dim1,Tensor_Dim2,Tensor_Dim3> &t,
                              const int &iterator0,
                              const int &iterator1,
                              const int &iterator2)
  {
    os << '[';
    for (int i = 0; i < Tensor_Dim3-1; ++i) {
      os << t(iterator0,iterator1,iterator2,i);
      os << ',';
    }
    os << t(iterator0,iterator1,iterator2,Tensor_Dim3-1);
    os << ']';

    return os;
  }

  template<class T,int Tensor_Dim0,int Tensor_Dim1,int Tensor_Dim2,int Tensor_Dim3>
  std::ostream & Tensor4_0010(std::ostream &os,
                              const Tensor4<T,Tensor_Dim0,Tensor_Dim1,Tensor_Dim2,Tensor_Dim3> &t,
                              const int &iterator0,
                              const int &iterator1)
  {
    os << '[';
    for (int i = 0; i < Tensor_Dim2-1; ++i) {
      FTensor::Tensor4_0001(os, t, iterator0, iterator1, i);
      os << ',';
    }
    FTensor::Tensor4_0001(os, t, iterator0, iterator1, Tensor_Dim2-1);
    os << ']';

    return os;
  }

  template<class T,int Tensor_Dim0,int Tensor_Dim1,int Tensor_Dim2,int Tensor_Dim3>
  std::ostream & Tensor4_0100(std::ostream &os,
                              const Tensor4<T,Tensor_Dim0,Tensor_Dim1,Tensor_Dim2,Tensor_Dim3> &t,
                              const int &iterator0)
  {
    os << '[';
    for (int i = 0; i < Tensor_Dim1-1; ++i) {
      FTensor::Tensor4_0010(os, t, iterator0, i);
      os << ',';
    }
    FTensor::Tensor4_0010(os, t, iterator0, Tensor_Dim1-1);
    os << ']';

    return os;
  }

}

template<class T,int Tensor_Dim0,int Tensor_Dim1,int Tensor_Dim2,int Tensor_Dim3>
std::ostream & operator<<(std::ostream &os,
                          const FTensor::Tensor4<T,Tensor_Dim0,Tensor_Dim1,Tensor_Dim2,Tensor_Dim3> &t)
{
  os << '[';
  for (int i = 0; i < Tensor_Dim0-1; ++i) {
    FTensor::Tensor4_0100(os, t, i);
    os << ',';
  }
  FTensor::Tensor4_0100(os, t, Tensor_Dim0-1);
  os << ']';

  return os;
}
