/* Declare a wrapper class for generic rank 4 Tensor expressions. */

#include "Tensor4_plus_Tensor4.hpp"
#include "Tensor4_times_Tensor2_symmetric.hpp"
#include "Tensor4_times_Tensor1.hpp"
#include "Tensor4_times_Tensor2.hpp"

template<class A, class T, int Dim0, int Dim1, int Dim2, int Dim3, char i, char j, char k, char l>
class Tensor4_Expr
{
  A iter;
public:
  Tensor4_Expr(A &a): iter(a) {}
  T operator()(const int N1, const int N2, const int N3, const int N4) const
  {
    return iter(N1,N2,N3,N4);
  }
};

template<class A, class T, int Dim0, int Dim1, int Dim2, int Dim3,
  char i, char j, char k, char l>
class Tensor4_Expr<Tensor4<A,Dim0,Dim1,Dim2,Dim3>,T,Dim0,Dim1,Dim2,Dim3,i,j,k,l>
{
  Tensor4<A,Dim0,Dim1,Dim2,Dim3> &iter;
public:
  Tensor4_Expr(Tensor4<A,Dim0,Dim1,Dim2,Dim3> &a): iter(a) {}
  T & operator()(const int N1, const int N2, const int N3, const int N4)
  {
    return iter(N1,N2,N3,N4);
  }
  T operator()(const int N1, const int N2, const int N3, const int N4) const
  {
    return iter(N1,N2,N3,N4);
  }

  /* Various assignment operators.  I have to explicitly declare the
     second operator= because otherwise the compiler will generate its
     own and not use the template code. */

  template<class B, class U>
  const Tensor4_Expr<Tensor4<A,Dim0,Dim1,Dim2,Dim3>,T,Dim0,Dim1,Dim2,Dim3,i,j,k,l> &
  operator=(const Tensor4_Expr<B,U,Dim0,Dim1,Dim2,Dim3,i,j,k,l> &result);

  const Tensor4_Expr<Tensor4<A,Dim0,Dim1,Dim2,Dim3>,T,Dim0,Dim1,Dim2,Dim3,i,j,k,l> &
  operator=(const Tensor4_Expr<Tensor4<A,Dim0,Dim1,Dim2,Dim3>,T,Dim0,Dim1,Dim2,Dim3,i,j,k,l> &result);

  /* 24 permutations of assignment for tensor4 */

  #define EQUAL(D0,D1,D2,D3,I,J,K,L) \
  template<class B, class U> inline \
  const Tensor4_Expr<Tensor4<A,Dim0,Dim1,Dim2,Dim3>,T, \
                     Dim0,Dim1,Dim2,Dim3,i,j,k,l> & \
  operator=(const Tensor4_Expr<B,U,D0,D1,D2,D3,I,J,K,L> &result);

  // ij jk kl
  // ik jl lk
  // il kj
  // ji lj
  // ki
  // li
  //
  // ijik ijlk ikjl iklj iljk ilkj jikl jilk kijl kilj lijk likj (12)
  // jkil jkli jlik jlki kjil kjli ljik ljki (8)
  // klij klji lkij lkji

  EQUAL(Dim1,Dim0,Dim2,Dim3,j,i,k,l); // jikl
  EQUAL(Dim1,Dim2,Dim0,Dim3,j,k,i,l); // jkil
  EQUAL(Dim1,Dim2,Dim3,Dim0,j,k,l,i); // jkli
  EQUAL(Dim2,Dim1,Dim3,Dim0,k,j,l,i); // kjli
  EQUAL(Dim2,Dim3,Dim1,Dim0,k,l,j,i); // klji
  EQUAL(Dim2,Dim3,Dim0,Dim1,k,l,i,j); // klij
  EQUAL(Dim3,Dim2,Dim0,Dim1,l,k,i,j); // lkij
  EQUAL(Dim3,Dim0,Dim2,Dim1,l,i,k,j); // likj
  EQUAL(Dim3,Dim0,Dim1,Dim2,l,i,j,k); // lijk
  EQUAL(Dim0,Dim3,Dim1,Dim2,i,l,j,k); // iljk
  EQUAL(Dim0,Dim1,Dim3,Dim2,i,j,l,k); // ijlk
  EQUAL(Dim3,Dim2,Dim1,Dim0,l,k,j,i); // lkji
  EQUAL(Dim0,Dim2,Dim1,Dim3,i,k,j,l); // ikjl
  EQUAL(Dim0,Dim2,Dim3,Dim1,i,k,l,j); // iklj
  EQUAL(Dim0,Dim3,Dim2,Dim1,i,l,k,j); // ilkj
  EQUAL(Dim1,Dim0,Dim3,Dim2,j,i,l,k); // jilk
  EQUAL(Dim2,Dim0,Dim1,Dim3,k,i,j,l); // kijl
  EQUAL(Dim2,Dim0,Dim3,Dim1,k,i,l,j); // kilj
  EQUAL(Dim1,Dim3,Dim0,Dim2,j,l,i,k); // jlik
  EQUAL(Dim1,Dim3,Dim2,Dim0,j,l,k,i); // jlki
  EQUAL(Dim2,Dim1,Dim0,Dim3,k,j,i,l); // kjil
  EQUAL(Dim3,Dim1,Dim0,Dim2,l,j,i,k); // ljik
  EQUAL(Dim3,Dim1,Dim2,Dim0,l,j,k,i); // ljki

  #undef EQUAL

  /* This is for int's, double's, etc. */

  template <class U> inline
  const Tensor4_Expr<Tensor4<A,Dim0,Dim1,Dim2,Dim3>,T,Dim0,Dim1,Dim2,Dim3,i,j,k,l> &
  operator=(const U &d);

  template <class U> inline
  const Tensor4_Expr<Tensor4<A,Dim0,Dim1,Dim2,Dim3>,T,Dim0,Dim1,Dim2,Dim3,i,j,k,l> &
  operator*=(const U &d);


};

#include "Tensor4_Expr_equals.hpp"
