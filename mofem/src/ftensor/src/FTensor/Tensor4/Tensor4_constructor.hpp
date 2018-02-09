/* A helper class that allows simple initialization of the Tensor4,
   but only if it has the correct number of elements. */

#pragma once

namespace FTensor {

template<class T, int Tensor_Dim0, int Tensor_Dim1, int Tensor_Dim2, int Tensor_Dim3>
class Tensor4_constructor;

template<class T>
class Tensor4_constructor<T,2,2,2,2>
{
public:
  Tensor4_constructor(T data[2][2][2][2],
    T d0000,T d0001,T d0010,T d0011,
    T d0100,T d0101,T d0110,T d0111,
    T d1000,T d1001,T d1010,T d1011,
    T d1100,T d1101,T d1110,T d1111
  )
  {
    data[0][0][0][0] = d0000; data[0][0][0][1] = d0001; data[0][0][1][0] = d0010; data[0][0][1][1] = d0011;
    data[0][1][0][0] = d0100; data[0][1][0][1] = d0101; data[0][1][1][0] = d0110; data[0][1][1][1] = d0111;
    data[1][0][0][0] = d1000; data[1][0][0][1] = d1001; data[1][0][1][0] = d1010; data[1][0][1][1] = d1011;
    data[1][1][0][0] = d1100; data[1][1][0][1] = d1101; data[1][1][1][0] = d1110; data[1][1][1][1] = d1111;
  }
};

template<class T>
class Tensor4_constructor<T,3,3,3,3>
{
public:
  Tensor4_constructor(T data[3][3][3][3],
    T d0000, T d0001, T d0002, T d0010, T d0011, T d0012, T d0020, T d0021, T d0022,
    T d0100, T d0101, T d0102, T d0110, T d0111, T d0112, T d0120, T d0121, T d0122,
    T d0200, T d0201, T d0202, T d0210, T d0211, T d0212, T d0220, T d0221, T d0222,
    T d1000, T d1001, T d1002, T d1010, T d1011, T d1012, T d1020, T d1021, T d1022,
    T d1100, T d1101, T d1102, T d1110, T d1111, T d1112, T d1120, T d1121, T d1122,
    T d1200, T d1201, T d1202, T d1210, T d1211, T d1212, T d1220, T d1221, T d1222,
    T d2000, T d2001, T d2002, T d2010, T d2011, T d2012, T d2020, T d2021, T d2022,
    T d2100, T d2101, T d2102, T d2110, T d2111, T d2112, T d2120, T d2121, T d2122,
    T d2200, T d2201, T d2202, T d2210, T d2211, T d2212, T d2220, T d2221, T d2222
  )
  {
    data[0][0][0][0] = d0000; data[0][0][0][1] = d0001; data[0][0][0][2] = d0002; data[0][0][1][0] = d0010; data[0][0][1][1] = d0011; data[0][0][1][2] = d0012; data[0][0][2][0] = d0020; data[0][0][2][1] = d0021; data[0][0][2][2] = d0022;
    data[0][1][0][0] = d0100; data[0][1][0][1] = d0101; data[0][1][0][2] = d0102; data[0][1][1][0] = d0110; data[0][1][1][1] = d0111; data[0][1][1][2] = d0112; data[0][1][2][0] = d0120; data[0][1][2][1] = d0121; data[0][1][2][2] = d0122;
    data[0][2][0][0] = d0200; data[0][2][0][1] = d0201; data[0][2][0][2] = d0202; data[0][2][1][0] = d0210; data[0][2][1][1] = d0211; data[0][2][1][2] = d0212; data[0][2][2][0] = d0220; data[0][2][2][1] = d0221; data[0][2][2][2] = d0222;
    data[1][0][0][0] = d1000; data[1][0][0][1] = d1001; data[1][0][0][2] = d1002; data[1][0][1][0] = d1010; data[1][0][1][1] = d1011; data[1][0][1][2] = d1012; data[1][0][2][0] = d1020; data[1][0][2][1] = d1021; data[1][0][2][2] = d1022;
    data[1][1][0][0] = d1100; data[1][1][0][1] = d1101; data[1][1][0][2] = d1102; data[1][1][1][0] = d1110; data[1][1][1][1] = d1111; data[1][1][1][2] = d1112; data[1][1][2][0] = d1120; data[1][1][2][1] = d1121; data[1][1][2][2] = d1122;
    data[1][2][0][0] = d1200; data[1][2][0][1] = d1201; data[1][2][0][2] = d1202; data[1][2][1][0] = d1210; data[1][2][1][1] = d1211; data[1][2][1][2] = d1212; data[1][2][2][0] = d1220; data[1][2][2][1] = d1221; data[1][2][2][2] = d1222;
    data[2][0][0][0] = d2000; data[2][0][0][1] = d2001; data[2][0][0][2] = d2002; data[2][0][1][0] = d2010; data[2][0][1][1] = d2011; data[2][0][1][2] = d2012; data[2][0][2][0] = d2020; data[2][0][2][1] = d2021; data[2][0][2][2] = d2022;
    data[2][1][0][0] = d2100; data[2][1][0][1] = d2101; data[2][1][0][2] = d2102; data[2][1][1][0] = d2110; data[2][1][1][1] = d2111; data[2][1][1][2] = d2112; data[2][1][2][0] = d2120; data[2][1][2][1] = d2121; data[2][1][2][2] = d2122;
    data[2][2][0][0] = d2200; data[2][2][0][1] = d2201; data[2][2][0][2] = d2202; data[2][2][1][0] = d2210; data[2][2][1][1] = d2211; data[2][2][1][2] = d2212; data[2][2][2][0] = d2220; data[2][2][2][1] = d2221; data[2][2][2][2] = d2222;
  }
};

}