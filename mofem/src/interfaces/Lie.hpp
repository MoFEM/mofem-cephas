/**
 * @file Lie.hpp
 * @brief Lie algebra implementation
 * @version 0.1
 * @date 2024-12-22
 *
 * @copyright Copyright (c) 2024
 *
 */

#pragma once

namespace LieGroups {

template <class T> struct TensorTypeExtractor {
  typedef typename std::remove_pointer<T>::type Type;
};
template <class T, int I> struct TensorTypeExtractor<FTensor::PackPtr<T, I>> {
  typedef typename std::remove_pointer<T>::type Type;
};

struct SO3 {

  SO3() = delete;
  ~SO3() = delete;

  template <typename T> inline static auto getVee(T &&w1, T &&w2, T &&w3) {
    return FTensor::Tensor1<T, dim>(

        std::forward<T>(w1), std::forward<T>(w2), std::forward<T>(w3)

    );
  }

  template <typename T> inline static auto getHat(T &&w1, T &&w2, T &&w3) {
    return getHatImpl(std::forward<A>(getVee(w1, w2, w3)));
  }

  template <typename A> inline static auto getVee(A &&t_w_hat) {
    return getVeeImpl(std::forward<A>(t_w_hat));
  }

  template <typename A> inline static auto getHat(A &&t_w_vee) {
    return getHatImpl(std::forward<A>(t_w_vee));
  }

  template <typename A, typename B>
  inline static auto exp(A &&t_w_vee, B &&theta) {
    return expImpl(std::forward<A>(t_w_vee), std::forward<B>(theta));
  }

  template <typename A, typename B>
  inline static auto Jl(A &&t_w_vee, B &&theta) {
    return JlImpl(std::forward<A>(t_w_vee), std::forward<B>(theta));
  }

  template <typename A, typename B>
  inline static auto Jr(A &&t_w_vee, B &&theta) {
    return JrImpl(std::forward<A>(t_w_vee), std::forward<B>(theta));
  }

  template <typename A, typename B, typename C>
  inline static auto action(A &&t_w_vee, B &&theta, C &&t_A) {
    return actionImpl(std::forward<A>(t_w_vee), std::forward<B>(theta),
                      std::forward<C>(t_A));
  }

  template <typename A, typename B>
  inline static auto dActionJl(A &&t_w_vee, B &&t_A) {
    return dActionJlImpl(std::forward<A>(t_w_vee), std::forward<B>(t_A));
  }

  template <typename A, typename B>
  inline static auto dActionJr(A &&t_w_vee, B &&t_A) {
    return dActionJrImpl(std::forward<A>(t_w_vee), std::forward<B>(t_A));
  }

  template <typename A, typename B>
  inline static auto diffExp(A &&t_w_vee, B &&theta) {
    return diffExpImpl(std::forward<A>(t_w_vee), std::forward<B>(theta));
  }

private:
  inline static constexpr int dim = 3;
  inline static FTENSOR_INDEX(dim, i);
  inline static FTENSOR_INDEX(dim, j);
  inline static FTENSOR_INDEX(dim, k);
  inline static FTENSOR_INDEX(dim, l);
  inline static FTENSOR_INDEX(dim, m);
  inline static FTENSOR_INDEX(dim, n);

  template <typename T>
  inline static auto getVeeImpl(FTensor::Tensor2<T, dim, dim> &t_w_hat) {
    using D = typename TensorTypeExtractor<T>::Type;
    FTensor::Tensor1<D, dim> t_w_vee;
    t_w_vee(k) = (levi_civita(i, j, k) * t_w_hat(i, j)) / 2;
    return t_w_vee;
  }

  template <typename T>
  inline static auto getHatImpl(FTensor::Tensor1<T, dim> &t_w_vee) {
    using D = typename TensorTypeExtractor<T>::Type;
    FTensor::Tensor2<D, dim, dim> t_w_hat;
    t_w_hat(i, j) = levi_civita(i, j, k) * t_w_vee(k);
    return t_w_hat;
  }

  template <typename T1, typename T2, typename T3>
  inline static auto genericFormImpl(FTensor::Tensor1<T1, dim> &t_w_vee,
                                     const T2 alpha, const T3 beta) {
    using D = typename TensorTypeExtractor<T1>::Type;
    FTensor::Tensor2<D, dim, dim> t_X;
    auto t_hat = getHat(t_w_vee);
    t_X(i, j) = FTensor::Kronecker_Delta<int>()(i, j) + alpha * t_hat(i, j) +
                beta * (t_hat(i, k) * t_hat(k, j));
    return t_X;
  }

  template <typename T1, typename T2>
  inline static auto expImpl(FTensor::Tensor1<T1, dim> &t_w_vee,
                             const T2 theta) {
    if (std::fabs(theta) < std::numeric_limits<T2>::epsilon()) {
      return genericFormImpl(t_w_vee, 1, 0.5);
    }
    const auto s = sin(theta);
    const auto s_half = sin(theta / 2);
    const auto a = s / theta;
    const auto b = 2 * (s_half / theta) * (s_half / theta);
    return genericFormImpl(t_w_vee, a, b);
  }

  template <typename T1, typename T2>
  inline static auto diffExpImpl(FTensor::Tensor1<T1, dim> &t_w_vee,
                                 const T2 theta) {

    auto get_tensor = [&t_w_vee](auto a, auto diff_a, auto b, auto diff_b) {
      FTENSOR_INDEX(3, i);
      FTENSOR_INDEX(3, j);
      FTENSOR_INDEX(3, k);
      FTENSOR_INDEX(3, l);

      using D = typename TensorTypeExtractor<T1>::Type;
      FTensor::Tensor3<D, 3, 3, 3> t_diff_exp;
      auto t_hat = getHat(t_w_vee);
      t_diff_exp(i, j, k) =

          a * FTensor::levi_civita<int>(i, j, k)

          +

          diff_a * t_hat(i, j) * t_w_vee(k)

          +

          b * (t_hat(i, l) * FTensor::levi_civita<int>(l, j, k) +
               FTensor::levi_civita<int>(i, l, k) * t_hat(l, j))

          +

          diff_b * t_hat(i, l) * t_hat(l, j) * t_w_vee(k);

      return t_diff_exp;
    };

    if(std::fabs(theta) < std::numeric_limits<T2>::epsilon()){
      return get_tensor(1., -1. / 3., 0., 1. / 1.2);
    }

    const auto ss = sin(theta);
    const auto a = ss / theta;

    const auto theta2 = theta * theta;
    const auto cc = cos(theta);
    const auto diff_a = (theta * cc - ss) / (theta2 * theta);

    const auto ss_2 = sin(theta / 2.);
    const auto cc_2 = cos(theta / 2.);
    const auto b = 2. * ss_2 * ss_2 / theta2;
    const auto diff_b =
        (2. * theta * ss_2 * cc_2 - 4. * ss_2 * ss_2) / (theta2 * theta2);

    return get_tensor(a, diff_a, b, diff_b);
  }

  template <typename T1, typename T2>
  inline static auto JlImpl(FTensor::Tensor1<T1, dim> &t_w_vee,
                            const T2 &theta) {
    if (std::fabs(theta) < std::numeric_limits<T2>::epsilon()) {
      return genericFormImpl(t_w_vee, 0.5, 1. / 6.);
    }
    const auto s = sin(theta);
    const auto s_half = sin(theta / 2);
    const auto a = 2 * (s_half / theta) * (s_half / theta);
    const auto b = ((theta - s) / theta) / theta / theta;
    return genericFormImpl(t_w_vee, a, b);
  }

  template <typename T1, typename T2>
  inline static auto JrImpl(FTensor::Tensor1<T1, dim> &t_w_vee,
                            const T2 theta) {
    if (std::fabs(theta) < std::numeric_limits<T2>::epsilon()) {
      return genericFormImpl(t_w_vee, -0.5, 1. / 6.);
    }
    const auto s = sin(theta);
    const auto s_half = sin(theta / 2);
    const auto a = 2 * (s_half / theta) * (s_half / theta);
    const auto b = ((theta - s) / theta) / theta / theta;
    return genericFormImpl(t_w_vee, -a, b);
  }

  template <typename T1, typename T2, typename T3>
  inline static auto actionImpl(FTensor::Tensor1<T1, dim> &t_w_vee,
                                const T2 theta,
                                FTensor::Tensor2_symmetric<T3, dim> &t_A) {
    using D = typename TensorTypeExtractor<T3>::Type;
    FTensor::Tensor2<D, dim, dim> t_B;
    t_B(i, j) = exp(t_w_vee, theta)(i,k) * t_A(k, j);
    return t_B;
  }


}; // namespace SO3

}; // namespace LieGroups