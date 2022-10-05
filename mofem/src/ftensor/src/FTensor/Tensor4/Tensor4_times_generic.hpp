namespace FTensor {
/* A(i,j,k,l)*generic */

template <class A, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3,
          char i, char j, char k, char l>
auto operator*(const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a,
               const U &d0) {
  auto TensorExpr = [&a, &d0](const int N1, const int N2, const int N3,
                              const int N4) {
    return a.operator()(N1, N2, N3, N4) * d0;
  };
  return Tensor4_Expr<decltype(TensorExpr), typename promote<T, U>::V, Dim0,
                      Dim1, Dim2, Dim3, i, j, k, l>(TensorExpr);
}

/* generic*A(i,j,k,l) */

template <class A, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3,
          char i, char j, char k, char l>
auto operator*(
    const U &d0,
    const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a) {
  auto TensorExpr = [&a, &d0](const int N1, const int N2, const int N3,
                              const int N4) {
    return d0 * a.operator()(N1, N2, N3, N4);
  };
  return Tensor4_Expr<decltype(TensorExpr), typename promote<T, U>::V, Dim0,
                      Dim1, Dim2, Dim3, i, j, k, l>(TensorExpr);
}

/* A(i,j,k,l)/generic */

template <class A, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3,
          char i, char j, char k, char l>
auto operator/(const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a,
               const U &d0) {
  auto TensorExpr = [&a, &d0](const int N1, const int N2, const int N3,
                              const int N4) {
    return a.operator()(N1, N2, N3) / d0;
  };
  return Tensor4_Expr<decltype(TensorExpr), typename promote<T, U>::V, Dim0,
                      Dim1, Dim2, Dim3, i, j, k, l>(TensorExpr);
}

/* generic/A(i,j,k,l) */

template <class A, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3,
          char i, char j, char k, char l>
auto operator/(
    const U &d0,
    const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a) {
  auto TensorExpr = [&a, &d0](const int N1, const int N2, const int N3,
                              const int N4) {
    return d0 / a.operator()(N1, N2, N3, N4);
  };
  return Tensor4_Expr<decltype(TensorExpr), typename promote<T, U>::V, Dim0,
                      Dim1, Dim2, Dim3, i, j, k, l>(TensorExpr);
}

} // namespace FTensor