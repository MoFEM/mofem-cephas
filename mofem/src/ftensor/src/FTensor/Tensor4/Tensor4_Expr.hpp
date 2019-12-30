/* Declare a wrapper class for generic rank 4 Tensor expressions. */

#pragma once

#include "Tensor4_minus_Tensor4.hpp"
#include "Tensor4_plus_Tensor4.hpp"
#include "Tensor4_or_Tensor4.hpp"
#include "Tensor4_times_Tensor2_single.hpp"
#include "Tensor4_times_Tensor2_double.hpp"
#include "Tensor4_times_Tensor2_symmetric.hpp"
#include "Tensor4_times_Tensor3_triple.hpp"
#include "Tensor4_times_Tensor3_double.hpp"

#include "../permute.hpp"

#include "Tensor4_times_Tensor1.hpp"

namespace FTensor {
template <class A, class T, int Dim0, int Dim1, int Dim2, int Dim3, char i,
          char j, char k, char l>
class Tensor4_Expr {
  A iter;

public:
  Tensor4_Expr(const A &a) : iter(a) {}
  T operator()(const int N1, const int N2, const int N3, const int N4) const {
    return iter(N1, N2, N3, N4);
  }
};

template <class A, class T, int Dim0, int Dim1, int Dim2, int Dim3, char i,
          char j, char k, char l>
class Tensor4_Expr<Tensor4<A, Dim0, Dim1, Dim2, Dim3>, T, Dim0, Dim1, Dim2,
                   Dim3, i, j, k, l> {
  Tensor4<A, Dim0, Dim1, Dim2, Dim3> &iter;

public:
  Tensor4_Expr(Tensor4<A, Dim0, Dim1, Dim2, Dim3> &a) : iter(a) {}
  T &operator()(const int N1, const int N2, const int N3, const int N4) {
    return iter(N1, N2, N3, N4);
  }
  T operator()(const int N1, const int N2, const int N3, const int N4) const {
    return iter(N1, N2, N3, N4);
  }

  /* Various assignment operators.  I have to explicitly declare the
     second operator= because otherwise the compiler will generate its
     own and not use the template code. */

  template <class B, class U, int Dim1_0, int Dim1_1, int Dim1_2, int Dim1_3,
            char i_1, char j_1, char k_1, char l_1>
  auto &equals(const Tensor4_Expr<B, U, Dim1_0, Dim1_1, Dim1_2, Dim1_3, i_1,
                                  j_1, k_1, l_1> &rhs) {
    auto index_sequence0 = std::make_index_sequence<Dim0>();
    auto index_sequence1 = std::make_index_sequence<Dim1>();
    auto index_sequence2 = std::make_index_sequence<Dim2>();
    auto index_sequence3 = std::make_index_sequence<Dim3>();

    auto l0 = [&](auto N0) {
      auto l1 = [&](auto N1) {
        auto l2 = [&](auto N2) {
          auto l3 = [&](auto N3) {
            iter(N0, N1, N2, N3) = permute(*this, rhs, N0, N1, N2, N3);
          };
          boost::hana::for_each(index_sequence3, l3);
        };
        boost::hana::for_each(index_sequence2, l2);
      };
      boost::hana::for_each(index_sequence1, l1);
    };
    boost::hana::for_each(index_sequence0, l0);

    return *this;
  }

  template <class B, class U, int Dim1_0, int Dim1_1, int Dim1_2, int Dim1_3,
            char i_1, char j_1, char k_1, char l_1>
  auto &operator=(const Tensor4_Expr<B, U, Dim1_0, Dim1_1, Dim1_2, Dim1_3, i_1,
                                     j_1, k_1, l_1> &rhs) {
    return equals(rhs);
  }

  auto &operator=(const Tensor4_Expr<Tensor4<A, Dim0, Dim1, Dim2, Dim3>, T,
                                     Dim0, Dim1, Dim2, Dim3, i, j, k, l> &rhs) {
    return equals(rhs);
  }

  template <class B, class U, int Dim1_0, int Dim1_1, int Dim1_2, int Dim1_3,
            char i_1, char j_1, char k_1, char l_1>
  auto &operator+=(const Tensor4_Expr<B, U, Dim1_0, Dim1_1, Dim1_2, Dim1_3, i_1,
                                      j_1, k_1, l_1> &rhs) {

    auto index_sequence0 = std::make_index_sequence<Dim0>();
    auto index_sequence1 = std::make_index_sequence<Dim1>();
    auto index_sequence2 = std::make_index_sequence<Dim2>();
    auto index_sequence3 = std::make_index_sequence<Dim3>();

    auto l0 = [&](auto N0) {
      auto l1 = [&](auto N1) {
        auto l2 = [&](auto N2) {
          auto l3 = [&](auto N3) {
            iter(N0, N1, N2, N3) += permute(*this, rhs, N0, N1, N2, N3);
          };
          boost::hana::for_each(index_sequence3, l3);
        };
        boost::hana::for_each(index_sequence2, l2);
      };
      boost::hana::for_each(index_sequence1, l1);
    };
    boost::hana::for_each(index_sequence0, l0);

    return *this;
  }

  template <class B, class U, int Dim1_0, int Dim1_1, int Dim1_2, int Dim1_3,
            char i_1, char j_1, char k_1, char l_1>
  auto &operator-=(const Tensor4_Expr<B, U, Dim1_0, Dim1_1, Dim1_2, Dim1_3, i_1,
                                      j_1, k_1, l_1> &rhs) {

    auto index_sequence0 = std::make_index_sequence<Dim0>();
    auto index_sequence1 = std::make_index_sequence<Dim1>();
    auto index_sequence2 = std::make_index_sequence<Dim2>();
    auto index_sequence3 = std::make_index_sequence<Dim3>();

    auto l0 = [&](auto N0) {
      auto l1 = [&](auto N1) {
        auto l2 = [&](auto N2) {
          auto l3 = [&](auto N3) {
            iter(N0, N1, N2, N3) -= permute(*this, rhs, N0, N1, N2, N3);
          };
          boost::hana::for_each(index_sequence3, l3);
        };
        boost::hana::for_each(index_sequence2, l2);
      };
      boost::hana::for_each(index_sequence1, l1);
    };
    boost::hana::for_each(index_sequence0, l0);

    return *this;
  }

  template <class U> auto &operator=(const U &d) {

    auto index_sequence0 = std::make_index_sequence<Dim0>();
    auto index_sequence1 = std::make_index_sequence<Dim1>();
    auto index_sequence2 = std::make_index_sequence<Dim2>();
    auto index_sequence3 = std::make_index_sequence<Dim3>();

    auto l0 = [&](auto N0) {
      auto l1 = [&](auto N1) {
        auto l2 = [&](auto N2) {
          auto l3 = [&](auto N3) { iter(N0, N1, N2, N3) = d; };
          boost::hana::for_each(index_sequence3, l3);
        };
        boost::hana::for_each(index_sequence2, l2);
      };
      boost::hana::for_each(index_sequence1, l1);
    };
    boost::hana::for_each(index_sequence0, l0);

    return *this;
  }

  template <class U> auto &operator*=(const U &d) {

    auto index_sequence0 = std::make_index_sequence<Dim0>();
    auto index_sequence1 = std::make_index_sequence<Dim1>();
    auto index_sequence2 = std::make_index_sequence<Dim2>();
    auto index_sequence3 = std::make_index_sequence<Dim3>();

    auto l0 = [&](auto N0) {
      auto l1 = [&](auto N1) {
        auto l2 = [&](auto N2) {
          auto l3 = [&](auto N3) { iter(N0, N1, N2, N3) *= d; };
          boost::hana::for_each(index_sequence3, l3);
        };
        boost::hana::for_each(index_sequence2, l2);
      };
      boost::hana::for_each(index_sequence1, l1);
    };
    boost::hana::for_each(index_sequence0, l0);

    return *this;
  }
};
} // namespace FTensor
