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

template <typename T, int Dim, int Shift> struct EigenProjection {

  using Val = FTensor::Tensor1<FTensor::PackPtr<T *, Shift>, Dim>;
  using Vec = FTensor::Tensor2<FTensor::PackPtr<T *, Shift>, Dim, Dim>;
  using Fun = boost::function<double(const double)>;

  template <int a, int i> static inline T N(Vec &t_vec) { return t_vec(a, i); }

  template <int a> static inline T L(Val &&t_val) { return t_val(a); }

  template <int a, int b> static inline T F(Val &t_val) {
    return static_cast<T>(1) / (F<a>(t_val) - F<b>(t_val));
  }

  template <int a, int b, int i, int j> static T M(Vec &t_vec) {
    return N<a, i>(t_vec) * N<b, j>(t_vec);
  }

  template <int a, int b, int i, int j>
  static inline T dFdN(Val &t_val, Vec &t_vec) {
    return dFdNa<a, b, i, j>(t_val, t_vec) + dFdNb<a, b, i, j>(t_val, t_vec);
  }

  template <int a, int b, int i, int j, int k, int l>
  static inline T G(Vec &t_vec) {
    return M<a, i, k>(t_vec) * M<b, j, l>(t_vec) +
           M<a, i, l>(t_vec) * M<b, j, k>(t_vec);
  }

  template <int a, int i, int j, int k, int l>
  static inline T d2M(Val &t_val, Vec &t_vec) {
    T ret = 0;

    boost::hana::for_each(

        boost::hana::make_range(boost::hana::int_c<0>, boost::hana::int_c<Dim>),

        [&](auto b) {
          if (a != b) {
            ret += F<a, b>(t_val) * S<a, b, i, j, k, l>(t_val);
          }
        }

    );

    return ret;
  }

  template <int a, int i, int j, int k, int l, int m, int n>
  static inline T dd4M(Val &t_val, Vec &t_vec) {
    T ret = 0;

    boost::hana::for_each(

        boost::hana::make_range(boost::hana::int_c<0>, boost::hana::int_c<Dim>),

        [&](auto b) {
          if (a != b) {
            ret += F<a, b>(t_val) * d2S<a, b, i, j, k, l, m, n>(t_val, t_vec) +
                   2 * dFdN<a, b, m, n>(t_val, t_vec) *
                       S<a, b, i, j, k, l>(t_val, t_vec);
          }
        }

    );

    return ret;
  }

  template <int NB, int i, int j>
  static inline T reconstructMatrix(Val &t_val, Vec &t_vec, Fun f) {
    T ret = 0;

    boost::hana::for_each(

        boost::hana::make_range(boost::hana::int_c<0>, boost::hana::int_c<NB>),

        [&](auto a) { ret += M<a, i, j>(t_vec) * f(L<a>(t_val)); }

    );

    return ret;
  }

  template <int NB, int i, int j, int k, int l>
  static inline T firstMatrixDirective(Val &t_val, Vec &t_vec, Fun f, Fun d_f) {
    T ret = 0;

    boost::hana::for_each(

        boost::hana::make_range(boost::hana::int_c<0>, boost::hana::int_c<NB>),

        [&](auto a) {
          ret += M<a, i, j>(t_vec) * M<a, k, l>(t_vec) * d_f(L<a>(t_val)) +
                 d2M<a, i, j, k, l>(t_val, t_vec) * f(L<a>(t_val)) /
                     static_cast<T>(2);
        }

    );

    return ret;
  }

  template <int NB, int i, int j, int k, int l, int m, int n>
  static inline T secondMatrixDirective(Val &t_val, Vec &t_vec, Fun f, Fun d_f,
                                        Fun dd_f) {
    T ret = 0;

    boost::hana::for_each(

        boost::hana::make_range(boost::hana::int_c<0>, boost::hana::int_c<NB>),

        [&](auto a) {
          ret +=

              d2M<a, i, j, m, n>(t_vec) * M<a, k, l>(t_vec) * d_f(L<a>(t_val)) /
                  static_cast<T>(2) +

              M<a, i, j>(t_vec) * d2M<a, k, l, m, n>(t_vec) * d_f(L<a>(t_val)) /
                  static_cast<T>(2) +

              M<a, i, j>(t_vec) * M<a, k, l>(t_vec) * dd_f(L<a>(t_val)) +

              dd4M<a, i, j, k, l, m, n>(t_val, t_vec) * f(L<a>(t_val)) /
                  static_cast<T>(4) +

              d2M<a, i, j, k, l>(t_val, t_vec) * d_f(L<a>(t_val)) /
                  static_cast<T>(2);
        }

    );

    return ret;
  }

  template <int NB> static inline auto getMat(Val &t_val, Vec &t_vec, Fun f) {
    FTensor::Tensor2<T, Dim, Dim> t_A;

    boost::hana::for_each(

        boost::hana::make_range(boost::hana::int_c<0>,
                                boost::hana::int_c<Dim * Dim>),

        [&](constexpr auto s) {
          constexpr int i = (s - s % Dim) / Dim;
          constexpr int j = s - 3 * j;
          t_A(i, j) = reconstructMatrix<NB, i, j>(t_val, t_vec, f);
        }

    );

    return t_A;
  }

  template <int NB>
  static inline auto getDiffMat(Val &t_val, Vec &t_vec, Fun f, Fun d_f,
                                FTensor::Number<NB>) {
    FTensor::Ddg<T, Dim, Dim> t_diff_A;

    boost::hana::for_each(

        boost::hana::make_range(boost::hana::int_c<0>,
                                boost::hana::int_c<Dim * Dim * Dim * Dim>),

        [&](auto s) {
          constexpr int i = (s - s % (3 * 3 * 3)) / (3 * 3 * 3);
          constexpr int j = (s - i * 3 * 3 * 3 - s % (3 * 3)) / (3 * 3);
          constexpr int k = (s - i * 3 * 3 * 3 - j * 3 * 3 - s % 3) / 3;
          constexpr int l = s - i * 3 * 3 * 3 - j * 3 * 3 - k * 3;

          if (i >= j && k >= l) {
            t_diff_A(i, j) =
                firstMatrixDirective<NB, i, j, k, l>(t_val, t_vec, f, d_f);
          }
        }

    );

    return t_diff_A;
  }

  template <int NB, typename A>
  static inline auto getDiffDiffMat(FTensor::Tensor2_symmetric<A, Dim> &t_s,
                                    Val &t_val, Vec &t_vec, Fun f, Fun d_f,
                                    Fun dd_f, FTensor::Number<NB>) {
    FTensor::Ddg<T, Dim, Dim> t_diff_diff_a;
    {
      FTensor::Index<'i', Dim> i;
      FTensor::Index<'j', Dim> j;
      FTensor::Index<'k', Dim> k;
      FTensor::Index<'l', Dim> l;
      t_diff_A(i, j, k, l) = 0;
    }

    boost::hana::for_each(

        boost::hana::make_range(boost::hana::int_c<0>,
                                boost::hana::int_c<Dim * Dim * Dim * Dim>),

        [&](auto s) {
          constexpr int k = (s - s % (3 * 3 * 3)) / (3 * 3 * 3);
          constexpr int l = (s - k * 3 * 3 * 3 - s % (3 * 3)) / (3 * 3);
          constexpr int m = (s - k * 3 * 3 * 3 - l * 3 * 3 - s % 3) / 3;
          constexpr int n = s - k * 3 * 3 * 3 - l * 3 * 3 - m * 3;

          if (k >= l && m >= n) {

            boost::hana::for_each(

                boost::hana::make_range(boost::hana::int_c<0>,
                                        boost::hana::int_c<Dim * Dim>),

                [&](auto r) {
                  constexpr int i = (s - s % Dim) / Dim;
                  constexpr int j = s - 3 * j;

                  double a;
                  if (i >= j) {
                    a = secondMatrixDirective<NB, i, j, k, l, m, n>(
                        t_val, t_vec, f, d_f, dd_f);
                    if (i != j)
                      t_diff_diff_a(k, l, m, n) += 2 * a;
                    else
                      t_diff_diff_a(k, l, m, n) += a;
                  }
                }

            );
          }
        }

    );

    return t_diff_diff_a;
  }

private:
  template <int a, int b, int i, int j>
  static inline T dFdNa(Val &t_val, Vec &t_vec) {
    return -M<a, a, i, j>(t_vec) /
           ((L<a>(t_val) - L<b>(t_val)) * (L<a>(t_val) - L<b>(t_val)));
  }

  template <int a, int b, int i, int j>
  static inline T dFdNb(Val &t_val, Vec &t_vec) {
    return M<b, b, i, j>(t_vec) /
           ((L<a>(t_val) - L<b>(t_val)) * (L<a>(t_val) - L<b>(t_val)));
  }

  template <int a, int b, int i, int j, int k, int l>
  static inline T S(Vec &t_vec) {
    return G<a, b, i, j, k, l>(t_vec) + G<b, a, i, j, k, l>(t_vec);
  }

  template <int a, int b, int i, int j, int k, int l, int m, int n>
  static inline T d2G(Val &t_val, Vec &t_vec) {
    return d2M<a, i, k, n, m>(t_val, t_vec) * M<b, j, l>(t_vec) +
           M<a, i, k>(t_vec) * d2M<b, j, l, m, n>(t_val, t_vec) +
           d2M<a, i, l, m, n>(t_val, t_vec) * M<b, j, k>(t_val) +
           M<a, i, l>(t_vec) * d2M<b, j, k, m, n>(t_val, t_vec);
  }

  template <int a, int b, int i, int j, int k, int l, int m, int n>
  static inline T d2S(Val &t_val, Vec &t_vec) {
    return d2G<a, b, i, j, k, l, m, n>(t_val, t_vec) +
           d2G<b, a, i, j, k, l, m, n>(t_val, t_vec);
  }
};

