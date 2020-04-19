/**
 * @file Kronecker_Delta.hpp
 *
 */

#pragma once

namespace FTensor {

/**
 * @brief Kronecker Delta class
 *
 * @tparam int
 */
template <class T = int> class Kronecker_Delta {
public:
  constexpr T operator()(const int N1, const int N2) const {
    return (N1 == N2) ? T(1) : T(0);
  }

  template <char i, char j, int Dim0, int Dim1>
  Tensor2_Expr<Kronecker_Delta<T>, T, Dim0, Dim1, i, j>
  operator()(const Index<i, Dim0> &, const Index<j, Dim1> &) const {
    return Tensor2_Expr<Kronecker_Delta<T>, T, Dim0, Dim1, i, j>(*this);
  };

  template <char i, int Dim0>
  constexpr auto operator()(const Index<i, Dim0> &, const int &N1) const {
    auto TensorExpr = [this, N1](const int &N0) {
      return this->operator()(N0, N1);
    };
    return Tensor1_Expr<decltype(TensorExpr), T, Dim0, i>(TensorExpr);
  };

  template <char j, int Dim1>
  constexpr auto operator()(const int &N0, const Index<j, Dim1> &) const {
    auto TensorExpr = [this, N0](const int &N1) {
      return this->operator()(N0, N1);
    };
    return Tensor1_Expr<decltype(TensorExpr), T, Dim1, j>{TensorExpr};
  };
};

/**
 * @brief Kronecker Delta class symmetric
 *
 * @tparam int
 */
template <class T = int> class Kronecker_Delta_symmetric {
public:
  constexpr T operator()(const int N1, const int N2) const {
    return (N1 == N2) ? T(1) : T(0);
  }

  template <char i, char j, int Dim>
  Tensor2_symmetric_Expr<Kronecker_Delta_symmetric<T>, T, Dim, i, j>
  operator()(const Index<i, Dim> &, const Index<j, Dim> &) const {
    return Tensor2_symmetric_Expr<Kronecker_Delta_symmetric<T>, T, Dim, i, j>(*this);
  };

  template <char i, int Dim0>
  constexpr auto operator()(const Index<i, Dim0> &, const int &N1) const {
    auto TensorExpr = [this, N1](const int &N0) {
      return this->operator()(N0, N1);
    };
    return Tensor1_Expr<decltype(TensorExpr), T, Dim0, i>(TensorExpr);
  };

  template <char j, int Dim1>
  constexpr auto operator()(const int &N0, const Index<j, Dim1> &) const {
    auto TensorExpr = [this, N0](const int &N1) {
      return this->operator()(N0, N1);
    };
    return Tensor1_Expr<decltype(TensorExpr), T, Dim1, j>{TensorExpr};
  };
};

/// Rank 2
template <class T = int, char i, char j, int Dim0, int Dim1>
Tensor2_Expr<Kronecker_Delta<T>, T, Dim0, Dim1, i, j>
kronecker_delta(const Index<i, Dim0> &, const Index<j, Dim1> &) {
  return Kronecker_Delta<T>()(Index<i, Dim0>(), Index<j, Dim1>());
}

template <class T = int, char i, int Dim0>
constexpr auto kronecker_delta(const Index<i, Dim0> &, const int &N1) {
  return Kronecker_Delta<T>()(Index<i, Dim0>(), N1);
}

template <class T = int, char j, int Dim1>
constexpr auto kronecker_delta(const int &N0, const Index<j, Dim1> &) {
  return Kronecker_Delta<T>()(N0, Index<j, Dim1>());
}

template <class T = int, char i, char j, int Dim>
Tensor2_symmetric_Expr<Kronecker_Delta_symmetric<T>, T, Dim, i, j>
kronecker_delta_symmetric(const Index<i, Dim> &, const Index<j, Dim> &) {
  return Kronecker_Delta_symmetric<T>()(Index<i, Dim>(), Index<j, Dim>());
}

template <class T = int, char i, int Dim0>
constexpr auto kronecker_delta_symmetric(const Index<i, Dim0> &,
                                         const int &N1) {
  return Kronecker_Delta_symmetric<T>()(Index<i, Dim0>(), N1);
}

template <class T = int, char j, int Dim1>
constexpr auto kronecker_delta_symmetric(const int &N0,
                                         const Index<j, Dim1> &) {
  return Kronecker_Delta_symmetric<T>()(N0, Index<j, Dim1>());
}

} // namespace FTensor