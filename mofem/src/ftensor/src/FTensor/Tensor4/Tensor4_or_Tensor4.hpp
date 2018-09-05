/* Adds a Tensor4 to a Tensor4, yielding a Ddg. */

#pragma once

namespace FTensor
{
  /* Base Template */
template <class A, class B, class T, class U, int Dim0_0, int Dim1_0,
          int Dim2_0, int Dim3_0, int Dim0_1, int Dim1_1, int Dim2_1,
          int Dim3_1, char i0, char j0, char k0, char l0, char i1, char j1,
          char k1, char l1>
class Tensor4_or_Tensor4 {};

/* A(i,j,k,l)+B(j,i,l,k)->Ddg */

template <class A, class B, class T, class U, int Dim01, int Dim23, char i,
          char j, char k, char l>
class Tensor4_or_Tensor4<A, B, T, U, Dim01, Dim01, Dim23, Dim23, Dim01, Dim01,
                         Dim23, Dim23, i, j, k, l, j, i, l, k> {
  Tensor4_Expr<A, T, Dim01, Dim01, Dim23, Dim23, i, j, k, l> iterA;
  Tensor4_Expr<B, U, Dim01, Dim01, Dim23, Dim23, j, i, l, k> iterB;

public:
  typename promote<T, U>::V operator()(const int N1, const int N2, const int N3,
                                       const int N4) const {
    return iterA(N1, N2, N3, N4) + iterB(N2, N1, N4, N3);
  }

  Tensor4_or_Tensor4(
      const Tensor4_Expr<A, T, Dim01, Dim01, Dim23, Dim23, i, j, k, l> &a,
      const Tensor4_Expr<B, U, Dim01, Dim01, Dim23, Dim23, j, i, l, k> &b)
      : iterA(a), iterB(b) {}
  };

  template <class A, class B, class T, class U, int Dim0_0, int Dim1_0,
            int Dim2_0, int Dim3_0, int Dim0_1, int Dim1_1, int Dim2_1,
            int Dim3_1, char i0, char j0, char k0, char l0, char i1, char j1,
            char k1, char l1>
  auto operator||(const Tensor4_Expr<A, T, Dim0_0, Dim1_0, Dim2_0, Dim3_0, i0,
                                     j0, k0, l0> &a,
                  const Tensor4_Expr<B, U, Dim0_1, Dim1_1, Dim2_1, Dim3_1, i1,
                                     j1, k1, l1> &b) {
    using TensorExpr =
        Tensor4_or_Tensor4<A, B, T, U, Dim0_0, Dim1_0, Dim2_0, Dim3_0, Dim0_1,
                           Dim1_1, Dim2_1, Dim3_1, i0, j0, k0, l0, i1, j1, k1,
                           l1>;
    static_assert(
        !std::is_empty<TensorExpr>::value,
        "Indexes or Dimensions are not compatible with the || operator");

    return Ddg_Expr<TensorExpr, typename promote<T, U>::V, Dim0_0, Dim2_1, i0,
                   j0, k0, l0>(TensorExpr(a, b));
  }
}
