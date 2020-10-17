/** \file FunctionMatrix.hpp
 * \brief Get function from matrix
 * \ingroup ftensor
 *
 */

/* This file is part of MoFEM.
 * MoFEM is free software: you can redistribute it and/or modify it under
 * the terms of the GNU Lesser General Public License as published by the
 * Free Software Foundation, either version 3 of the License, or (at your
 * option) any later version.
 *
 * MoFEM is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
 * License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with MoFEM. If not, see <http://www.gnu.org/licenses/>. */

#pragma once
#include <type_traits>

template <typename E, typename C, int a, int i, int j, int k, int l>
struct d2MImpl {
  using Val = typename E::Val;
  using Vec = typename E::Vec;
 
  template <int N> using Number = FTensor::Number<N>;

  d2MImpl() = delete;
  ~d2MImpl() = delete;

  template <int b> static inline C term(Val &t_val, Vec &t_vec) {
    if (a != b)
      return E::F(t_val, Number<a>(), Number<b>()) *
             E::S(t_vec, Number<a>(), Number<b>(), Number<i>(), Number<j>(),
                  Number<k>(), Number<l>());
    else
      return 0;
  }

  template <int nb>
  static inline C eval(Val &t_val, Vec &t_vec, const Number<nb> &) {
    return term<nb - 1>(t_val, t_vec) + eval(t_val, t_vec, Number<nb - 1>());
  }

  static inline C eval(Val &t_val, Vec &t_vec, const Number<1> &) {
    return term<0>(t_val, t_vec);
  }
};

// template <typename E, typename C, int NB, int a, int i, int j, int k, int l,
//           int m, int n>
// struct dd4MImpl : public dd4MImpl<E, C, 1, a, i, j, k, l, m, n> {
//   using I = dd4MImpl<E, C, 1, a, i, j, k, l, m, n>;
//   using Val = typename E::Val;
//   using Vec = typename E::Vec;
//   dd4MImpl() = delete;
//   ~dd4MImpl() = delete;
//   static inline C eval(Val &t_val, Vec &t_vec) {
//     return I::term<NB - 1>(t_val, t_vec) +
//            dd4MImpl<E, C, NB - 1, a, i, j, k, l, m, n>::eval(t_val, t_vec);
//   }
// };

template <typename E, typename C, int i, int j> struct reconstructMatImpl {
  using Val = typename E::Val;
  using Vec = typename E::Vec;
  using Fun = typename E::Fun;

  template <int N> using Number = FTensor::Number<N>;

  reconstructMatImpl() = delete;
  ~reconstructMatImpl() = delete;

  template <int a> static inline C term(Val &t_val, Vec &t_vec, Fun f) {
    return E::M(t_vec, Number<a>(), Number<a>(), Number<i>(), Number<j>()) *
           f(E::L(t_val, Number<a>()));
  }

  template <int nb>
  static inline C eval(Val &t_val, Vec &t_vec, Fun f, const Number<nb> &) {
    return term<nb - 1>(t_val, t_vec, f) +
           eval(t_val, t_vec, f, Number<nb - 1>());
  }

  static inline C eval(Val &t_val, Vec &t_vec, Fun f, const Number<1> &) {
    return term<0>(t_val, t_vec, f);
  }
};

template <typename E, typename C, int NB, int i, int j, int k, int l>
struct firstMatrixDirectiveImpl {
  using Val = typename E::Val;
  using Vec = typename E::Vec;
  using Fun = typename E::Fun;

  template <int N> using Number = FTensor::Number<N>;

  firstMatrixDirectiveImpl() = delete;
  ~firstMatrixDirectiveImpl() = delete;

  template <int a>
  static inline C term(Val &t_val, Vec &t_vec, Fun f, Fun d_f) {
    return

        E::M(t_vec, Number<a>(), Number<a>(), Number<i>(), Number<j>()) *
            E::M(t_vec, Number<a>(), Number<a>(), Number<k>(), Number<l>()) *
            d_f(E::L(t_val, Number<a>()))

        +

        E::d2M(t_val, t_vec, Number<NB>(), Number<a>(), Number<i>(),
               Number<j>(), Number<k>(), Number<l>()) *
            f(E::L(t_val, Number<a>())) / static_cast<C>(2);
  }

  template <int nb>
  static inline C eval(Val &t_val, Vec &t_vec, Fun f, Fun d_f,
                       const Number<nb> &) {
    return term<nb - 1>(t_val, t_vec, f, d_f) +
           eval(t_val, t_vec, f, d_f, Number<nb - 1>());
  }

  static inline C eval(Val &t_val, Vec &t_vec, Fun f, Fun d_f,
                       const Number<1> &) {
    return term<0>(t_val, t_vec, f, d_f);
  }
};

// template <typename E, typename C, int NB, int i, int j, int k, int l, int m,
//           int n>
// struct secondMatrixDirectiveImpl
//     : public secondMatrixDirectiveImpl<E, C, 1, i, j, k, l, m, n> {
//   using I = secondMatrixDirectiveImpl<E, C, 1, i, j, k, l, m, n>;
//   using Val = typename E::Val;
//   using Vec = typename E::Vec;
//   using Fun = typename E::Fun;
//   secondMatrixDirectiveImpl() = delete;
//   ~secondMatrixDirectiveImpl() = delete;
//   static inline C eval(Val &t_val, Vec &t_vec, Fun f, Fun d_f, Fun dd_f) {
//     return I::term<NB - 1>(t_val, t_vec, f, d_f, dd_f) +
//            secondMatrixDirectiveImpl<E, C, NB - 1, i, j, k, l, m, n>::eval(
//                t_val, t_vec, f, d_f, dd_f);
//   }
// };

template <typename E, typename C, int NB> struct getMatImpl {
  using Val = typename E::Val;
  using Vec = typename E::Vec;
  using Fun = typename E::Fun;

  template <int N> using Number = FTensor::Number<N>;

  getMatImpl() = delete;
  ~getMatImpl() = delete;

  template <typename T, int I, int J>
  static inline void set(Val &t_val, Vec &t_vec, Fun f, T &t_a,
                         const Number<I> &, const Number<J> &) {
    set(t_val, t_vec, f, t_a, Number<I>(), Number<J - 1>());
    t_a(I - 1, J - 1) = reconstructMatImpl<E, C, I - 1, J - 1>::eval(
        t_val, t_vec, f, Number<NB>());
  };

  template <typename T, int I>
  static inline void set(Val &t_val, Vec &t_vec, Fun f, T &t_a,
                         const Number<I> &, const Number<1> &) {
    set(t_val, t_vec, f, t_a, Number<I - 1>(), Number<I - 1>());
    t_a(I - 1, 0) =
        reconstructMatImpl<E, C, I - 1, 0>::eval(t_val, t_vec, f, Number<NB>());
  };

  template <typename T>
  static inline void set(Val &t_val, Vec &t_vec, Fun f, T &t_a,
                         const Number<1> &, const Number<1> &) {
    t_a(0, 0) =
        reconstructMatImpl<E, C, 0, 0>::eval(t_val, t_vec, f, Number<NB>());
  };
};

template <typename E, typename C, int NB, int Dim> struct getDiffMatImpl {
  using Val = typename E::Val;
  using Vec = typename E::Vec;
  using Fun = typename E::Fun;

  template <int N> using Number = FTensor::Number<N>;

  getDiffMatImpl() = delete;
  ~getDiffMatImpl() = delete;

  template <typename T, int I, int J, int K, int L>
  static inline void set(Val &t_val, Vec &t_vec, Fun f, Fun d_f, T &t_a,
                         const Number<I> &, const Number<J> &,
                         const Number<K> &, const Number<L> &) {
    set(t_val, t_vec, f, d_f, t_a, Number<I>(), Number<J>(), Number<K>(),
        Number<L - 1>());
    t_a(I - 1, J - 1, K - 1, L - 1) =
        firstMatrixDirectiveImpl<E, C, NB, I - 1, J - 1, K - 1, L - 1>::eval(
            t_val, t_vec, f, d_f, Number<NB>());
  };

  template <typename T, int I, int J, int K>
  static inline void set(Val &t_val, Vec &t_vec, Fun f, Fun d_f, T &t_a,
                         const Number<I> &, const Number<J> &,
                         const Number<K> &, const Number<1> &) {
    // set(t_val, t_vec, f, d_f, t_a, Number<I>(), Number<J>(), Number<K - 1>(),
    //     Number<K - 1>());
    t_a(I - 1, J - 1, K - 1, 0) =
        firstMatrixDirectiveImpl<E, C, NB, I - 1, J - 1, K - 1, 0>::eval(
            t_val, t_vec, f, d_f, Number<NB>());
  };

  template <typename T, int I, int J>
  static inline void set(Val &t_val, Vec &t_vec, Fun f, Fun d_f, T &t_a,
                         const Number<I> &, const Number<I> &,
                         const Number<1> &, const Number<1> &) {
    // set(t_val, t_vec, f, d_f, t_a, Number<I>(), Number<J - 1>(), Number<Dim>(),
        // Number<Dim>());
    t_a(I - 1, J - 1, 0, 0) =
        firstMatrixDirectiveImpl<E, C, NB, I - 1, J - 1, 0, 0>::eval(
            t_val, t_vec, f, d_f, Number<NB>());
  };

  template <typename T, int I, int K, int L>
  static inline void set(Val &t_val, Vec &t_vec, Fun f, Fun d_f, T &t_a,
                         const Number<I> &, const Number<1> &,
                         const Number<K> &, const Number<L> &) {
    // set(t_val, t_vec, f, d_f, t_a, Number<I - 1>(), Number<I - 1>(),
    //     Number<K>(), Number<L>());
    t_a(I - 1, 0, K - 1, L - 1) =
        firstMatrixDirectiveImpl<E, C, NB, I - 1, 0, K - 1, L - 1>::eval(
            t_val, t_vec, f, d_f, Number<NB>());
  };

  template <typename T>
  static inline void set(Val &t_val, Vec &t_vec, Fun f, Fun d_f, T &t_a,
                         const Number<1> &, const Number<1> &,
                         const Number<1> &, const Number<1> &) {
    t_a(0, 0, 0, 0) = firstMatrixDirectiveImpl<E, C, NB, 0, 0, 0, 0>::eval(
        t_val, t_vec, f, d_f, Number<NB>());
  };
};

template <typename T1, typename T2, int Dim = 3>
struct EigenProjection {

  using Val = const FTensor::Tensor1<T1, Dim>;
  using Vec = const FTensor::Tensor2<T2, Dim, Dim>;
  using Fun = boost::function<double(const double)>;

  template <int N> using Number = FTensor::Number<N>;

  template <int a, int i> static inline auto N(Vec &t_vec) {
    return t_vec(a, i);
  }

  template <int a> static inline auto L(Val &t_val, const  Number<a> &) {
    return L<a>(t_val);
  }
  template <int a> static inline auto L(Val &t_val) { return t_val(a); }

  template <int a, int b>
  static inline auto F(Val &t_val, const Number<a> &, const Number<b> &) {
    return F<a, b>(t_val);
  }

  template <int a, int b> static inline auto F(Val &t_val) {
    return static_cast<decltype(t_val(0))>(1) / (L<a>(t_val) - L<b>(t_val));
  }

  template <int a, int b, int i, int j>
  static inline auto M(Vec &t_vec, const Number<a> &, const Number<b> &,
                       const Number<i> &, const Number<j> &) {
    return M<a, b, i, j>(t_vec);
  }
  template <int a, int b, int i, int j> static inline auto M(Vec &t_vec) {
    return N<a, i>(t_vec) * N<b, j>(t_vec);
  }

  template <int a, int b, int i, int j>
  static inline auto dFdN(Val &t_val, Vec &t_vec) {
    return dFdNa<a, b, i, j>(t_val, t_vec) + dFdNb<a, b, i, j>(t_val, t_vec);
  }

  template <int a, int b, int i, int j, int k, int l>
  static inline auto G(Vec &t_vec) {
    return M<a, a, i, k>(t_vec) * M<b, b, j, l>(t_vec) +
           M<a, a, i, l>(t_vec) * M<b, b, j, k>(t_vec);
  }

  template <int nb, int a, int i, int j, int k, int l>
  static inline auto d2M(Val &t_val, Vec &t_vec, const Number<nb> &,
                         const Number<a> &, const Number<i> &,
                         const Number<j> &, const Number<k> &,
                         const Number<l> &) {
    return d2M<nb, a, i, j, k, l>(t_val, t_vec);
  }

  template <int nb, int a, int i, int j, int k, int l>
  static inline auto d2M(Val &t_val, Vec &t_vec) {
    using V =
        typename FTensor::promote<decltype(t_val(0)), decltype(t_vec(0, 0))>::V;
    return d2MImpl<EigenProjection<T1, T2, Dim>, V, a, i, j, k, l>::eval(
        t_val, t_vec, Number<nb>());
  }

  // template <int NB, int a, int i, int j, int k, int l, int m, int n>
  // static inline auto dd4M(Val &t_val, Vec &t_vec) {
  //   using V =
  //       typename FTensor::promote<decltype(t_val(0)), decltype(t_vec(0, 0))>::V;
  //   return dd4MImpl<EigenProjection<T1, T2, Dim>, V, NB, a, i, j, k, l, m,
  //                   n>::eval(t_val, t_vec);
  // }

  // template <int NB, int i, int j, int k, int l, int m, int n>
  // static inline auto secondMatrixDirective(Val &t_val, Vec &t_vec, Fun f,
  //                                          Fun d_f, Fun dd_f) {
  //   using V =
  //       typename FTensor::promote<decltype(t_val(0)), decltype(t_vec(0, 0))>::V;
  //   return secondMatrixDirectiveImpl<EigenProjection<T1, T2, Dim>, V, NB, i, j,
  //                                    k, l, m, n>::eval(t_val, t_vec, f, d_f,
  //                                                      dd_f);
  // }

  template<int nb>
  static inline auto getMat(Val &t_val, Vec &t_vec, Fun f) {
    using V =
        typename FTensor::promote<decltype(t_val(0)), decltype(t_vec(0, 0))>::V;
    FTensor::Tensor2_symmetric<typename std::remove_const<V>::type, Dim> t_A;
    getMatImpl<EigenProjection<T1, T2, Dim>, V, nb>::set(
        t_val, t_vec, f, t_A, Number<Dim>(), Number<Dim>());
    return t_A;
  }

  template <int nb>
  static inline auto getDiffMat(Val &t_val, Vec &t_vec, Fun f, Fun d_f) {
    using V =
        typename FTensor::promote<decltype(t_val(0)), decltype(t_vec(0, 0))>::V;
    FTensor::Ddg<V, Dim, Dim> t_diff_A;
    getDiffMatImpl<EigenProjection<T1, T2, Dim>, V, nb, Dim>::set(
        t_val, t_vec, f, d_f, t_diff_A, Number<Dim>(), Number<Dim>(),
        Number<Dim>(), Number<Dim>());
    return t_diff_A;
  }

  // template <int NB, typename A>
  // static inline auto getDiffDiffMat(FTensor::Tensor2_symmetric<A, Dim> &t_s,
  //                                   Val &t_val, Vec &t_vec, Fun f, Fun d_f,
  //                                   Fun dd_f, FTensor::Number<NB>) {

  //   using V =
  //       typename FTensor::promote<decltype(t_val(0)), decltype(t_vec(0, 0))>::V;

  //   FTensor::Ddg<V, Dim, Dim> t_diff_diff_a;
  //   {
  //     FTensor::Index<'i', Dim> i;
  //     FTensor::Index<'j', Dim> j;
  //     FTensor::Index<'k', Dim> k;
  //     FTensor::Index<'l', Dim> l;
  //     t_diff_A(i, j, k, l) = 0;
  //   }

  //   // boost::hana::for_each(

  //   //     boost::hana::make_range(boost::hana::int_c<0>,
  //   //                             boost::hana::int_c<Dim * Dim * Dim * Dim>),

  //   //     [&](auto s) {
  //   //       constexpr int k = (s - s % (3 * 3 * 3)) / (3 * 3 * 3);
  //   //       constexpr int l = (s - k * 3 * 3 * 3 - s % (3 * 3)) / (3 * 3);
  //   //       constexpr int m = (s - k * 3 * 3 * 3 - l * 3 * 3 - s % 3) / 3;
  //   //       constexpr int n = s - k * 3 * 3 * 3 - l * 3 * 3 - m * 3;

  //   //       if (k >= l && m >= n) {

  //   //         boost::hana::for_each(

  //   //             boost::hana::make_range(boost::hana::int_c<0>,
  //   //                                     boost::hana::int_c<Dim * Dim>),

  //   //             [&](auto r) {
  //   //               constexpr int i = (s - s % Dim) / Dim;
  //   //               constexpr int j = s - 3 * j;

  //   //               double a;
  //   //               if (i >= j) {
  //   //                 a = secondMatrixDirective<NB, i, j, k, l, m, n>(
  //   //                     t_val, t_vec, f, d_f, dd_f);
  //   //                 if (i != j)
  //   //                   t_diff_diff_a(k, l, m, n) += 2 * a;
  //   //                 else
  //   //                   t_diff_diff_a(k, l, m, n) += a;
  //   //               }
  //   //             }

  //   //         );
  //   //       }
  //   //     }

  //   // );

  //   return t_diff_diff_a;
  // }

  template <int a, int b, int i, int j>
  static inline auto dFdNa(Val &t_val, Vec &t_vec) {
    return -M<a, a, i, j>(t_vec) /
           ((L<a>(t_val) - L<b>(t_val)) * (L<a>(t_val) - L<b>(t_val)));
  }

  template <int a, int b, int i, int j>
  static inline auto dFdNb(Val &t_val, Vec &t_vec) {
    return M<b, b, i, j>(t_vec) /
           ((L<a>(t_val) - L<b>(t_val)) * (L<a>(t_val) - L<b>(t_val)));
  }

  template <int a, int b, int i, int j, int k, int l>
  static inline auto S(Vec &t_vec, const Number<a> &, const Number<b> &,
                       const Number<i> &, const Number<j> &, const Number<k> &,
                       const Number<l> &) {
    return S<a, b, i, j, k, l>(t_vec);
  }

  template <int a, int b, int i, int j, int k, int l>
  static inline auto S(Vec &t_vec) {
    return G<a, b, i, j, k, l>(t_vec) + G<b, a, i, j, k, l>(t_vec);
  }

  // template <int a, int b, int i, int j, int k, int l, int m, int n>
  // static inline auto d2G(Val &t_val, Vec &t_vec) {
  //   return d2M<a, i, k, n, m>(t_val, t_vec) * M<b, j, l>(t_vec) +
  //          M<a, i, k>(t_vec) * d2M<b, j, l, m, n>(t_val, t_vec) +
  //          d2M<a, i, l, m, n>(t_val, t_vec) * M<b, j, k>(t_val) +
  //          M<a, i, l>(t_vec) * d2M<b, j, k, m, n>(t_val, t_vec);
  // }

  // template <int a, int b, int i, int j, int k, int l, int m, int n>
  // static inline auto d2S(Val &t_val, Vec &t_vec) {
  //   return d2G<a, b, i, j, k, l, m, n>(t_val, t_vec) +
  //          d2G<b, a, i, j, k, l, m, n>(t_val, t_vec);
  // }
};

//   static inline C eval(Val &t_val, Vec &t_vec) { return term<0>(t_val, t_vec); }
// };

// template <typename E, typename C, int a, int i, int j, int k, int l, int m,
//           int n>
// struct dd4MImpl<E, C, 1, a, i, j, k, l, m, n> {
//   using Val = typename E::Val;
//   using Vec = typename E::Vec;
//   dd4MImpl() = delete;
//   ~dd4MImpl() = delete;
//   template <int b> static inline C term(Val &t_val, Vec &t_vec) {
//     if (a != b)
//       return E::F<a, b>(t_val) * E::d2S<a, b, i, j, k, l, m, n>(t_val, t_vec) +
//              2 * E::dFdN<a, b, m, n>(t_val, t_vec) *
//                  E::S<a, b, i, j, k, l>(t_val, t_vec);
//     else
//       return 0;
//   }
//   static inline C eval(Val &t_val, Vec &t_vec) { return term<0>(t_val, t_vec); }
// };


// template <typename E, typename C, int i, int j, int k, int l, int m, int n>
// struct secondMatrixDirectiveImpl<E, C, 1, i, j, k, l, m, n> {
//   using Val = typename E::Val;
//   using Vec = typename E::Vec;
//   using Fun = typename E::Fun;
//   secondMatrixDirectiveImpl() = delete;
//   ~secondMatrixDirectiveImpl() = delete;
//   template <int a>
//   static inline C term(Val &t_val, Vec &t_vec, Fun f, Fun d_f, Fun dd_f) {
//     return

//         E::d2M<a, i, j, m, n>(t_vec) * E::M<a, k, l>(t_vec) *
//             d_f(E::L<a>(t_val)) / static_cast<C>(2) +

//         E::M<a, i, j>(t_vec) * E::d2M<a, k, l, m, n>(t_vec) *
//             d_f(E::L<a>(t_val)) / static_cast<C>(2) +

//         E::M<a, i, j>(t_vec) * E::M<a, k, l>(t_vec) * dd_f(E::L<a>(t_val)) +

//         E::dd4M<a, i, j, k, l, m, n>(t_val, t_vec) * f(E::L<a>(t_val)) /
//             static_cast<C>(4) +

//         E::d2M<a, i, j, k, l>(t_val, t_vec) * d_f(E::L<a>(t_val)) /
//             static_cast<C>(2);
//   }

//   static inline C eval(Val &t_val, Vec &t_vec, Fun f, Fun d_f, Fun dd_f) {
//     return term<0>(t_val, t_vec, f, d_f, dd_f);
//   }
// };