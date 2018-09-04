/* This is for expressions where a number is used for one. two or three
   slots, and an index for the others, yielding a Tensor1_Expr, Tensor2_Expr or
   Tensor3_Expr. */

#pragma once

namespace FTensor
{

  /* Third slot. */

  template <class A, class T> class Tensor4_numeral_2
  {
    A iterA;
    int N;

  public:
    T operator()(const int N0, const int N1, const int N2) const {
      return iterA(N0, N1, N, N2);
    }
    Tensor4_numeral_2(const A &a, const int NN) : iterA(a), N(NN) {}
  };

  /* Forth slot. */

  template <class A, class T> class Tensor4_numeral_3
  {
    A iterA;
    int N;

  public:
    T operator()(const int N0, const int N1, const int N2) const {
      return iterA(N0, N1, N2, N);
    }
    Tensor4_numeral_3(const A &a, const int NN) : iterA(a), N(NN) {}
  };


}
