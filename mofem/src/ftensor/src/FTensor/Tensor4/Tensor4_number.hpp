/* This is for expressions where a number is used for one, two or three
   slots, and an index for the others, yielding a Tensor1_Expr, Tenors2_Expr or
   Tensor3_Expr. */

#pragma once

namespace FTensor
{

  /* Third slot. */

  template <class A, class T, int N> class Tensor4_number_2
  {
    A iterA;

  public:
    T operator()(const int N0, const int N1, const int N3) const {
      return iterA(N0, N1, N, N3);
    }
    Tensor4_number_2(const A &a) : iterA(a) {}
  };

  template <class A, class T, int N> class Tensor4_number_rhs_2
  {};

  /* Forth slot. */

  template <class A, class T, int N> class Tensor4_number_3
  {
    A iterA;

  public:
    T operator()(const int N0, const int N1, const int N2) const {
      return iterA(N0, N1, N2, N);
    }
    Tensor4_number_3(const A &a) : iterA(a) {}
  };

  template <class A, class T, int N> class Tensor4_number_rhs_3
  {};

}
