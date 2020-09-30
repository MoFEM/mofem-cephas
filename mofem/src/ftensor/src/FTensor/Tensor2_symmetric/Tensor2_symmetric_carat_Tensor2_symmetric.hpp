/* Multiply a Tensor2_symmetric cart a Tensor2_symmetric together and yieald DDg
 * tensor. */

/* A(i,j) ^ B(i,j) -> Ddg */

#pragma once

namespace FTensor {

template <class A, class B, class T, class U, int Dim, char i, char j, char k,
          char l>
class Tensor2_symmetric_carat_Tensor2_symmetric {
  const Tensor2_symmetric_Expr<A, T, Dim, i, k> iterA;
  const Tensor2_symmetric_Expr<B, T, Dim, j, l> iterB;

  inline typename promote<T, U>::V eval(const int N1, const int N2,
                                        const int N3, const int N4) const {
    auto small_eval = [&](auto n1, auto n2, auto n3, auto n4) {
      return iterA(n1, n3) * iterB(n2, n4);
    };
    return small_eval(N1, N2, N3, N4) + small_eval(N2, N1, N3, N4) +
           small_eval(N1, N2, N4, N3) + small_eval(N2, N1, N4, N3);
  }

public:
  Tensor2_symmetric_carat_Tensor2_symmetric(
      const Tensor2_symmetric_Expr<A, T, Dim, i, k> &a,
      const Tensor2_symmetric_Expr<B, U, Dim, j, l> &b)
      : iterA(a), iterB(b) {}
  typename promote<T, U>::V operator()(const int N1, const int N2, const int N3,
                                       const int N4) const {
    return eval(N1, N2, N3, N4);
  }
};

template <class A, class B, class T, class U, int Dim, char i, char j, char k,
          char l>
Ddg_Expr<Tensor2_symmetric_carat_Tensor2_symmetric<A, B, T, U, Dim, i, j, k, l>,
         typename promote<T, U>::V, Dim, Dim, i, j, k, l>
operator^(const Tensor2_symmetric_Expr<A, T, Dim, i, k> &a,
          const Tensor2_symmetric_Expr<B, U, Dim, j, l> &b) {
  using TensorExpr =
      Tensor2_symmetric_carat_Tensor2_symmetric<A, B, T, U, Dim, i, j, k, l>;
  return Ddg_Expr<TensorExpr, typename promote<T, U>::V, Dim, Dim, i, j, k, l>(
      TensorExpr(a, b));
}

} // namespace FTensor