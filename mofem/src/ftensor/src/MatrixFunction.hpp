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

namespace MatrixFunction {

template <typename T, int DIM, int SHIFT> struct EigenProjection {

  using Val = FTensor::Tensor1<FTensor::PackPtr<T *, SHIFT>, DIM>;
  using Vec = FTensor::Tensor2<FTensor::PackPtr<T *, SHIFT>, DIM, DIM>;
  using Fun = boost::function<double(const double)>;

  template <int a, int i> static inline T N(Vec &vec) { return vec(a, i); }

  template <int a> static inline T L(Val &&val) { return val(a); }

  template <int a, int b> static inline T F(Val &val) {
    return static_cast<T>(1) / (F<a>(val) - F<b>(val));
  }

  template <int a, int b, int i, int j> static T M(Vec &vec) {
    return N<a, i>(vec) * N<b, j>(vec);
  }

  template <int a, int b, int i, int j>
  static inline T dFdN(Val &val, Vec &vec) {
    return dFdNa<a, b, i, j>(val, vec) + dFdNb<a, b, i, j>(val, vec);
  }

  template <int a, int b, int i, int j, int k, int l>
  static inline T G(Vec &vec) {
    return M<a, i, k>(vec) * M<b, j, l>(vec) +
           M<a, i, l>(vec) * M<b, j, k>(vec);
  }

  template <int a, int i, int j, int k, int l>
  static inline T d2M(Val &val, Vec &vec) {
    T ret = 0;

    boost::hana::for_each(

        boost::hana::make_range(boost::hana::int_c<0>, boost::hana::int_c<DIM>),

        [&](auto b) {
          if (a != b) {
            ret += F<a, b>(val) * S<a, b, i, j, k, l>(val);
          }
        }

    );

    return ret;
  }

  template <int a, int i, int j, int k, int l, int m, int n>
  static inline T dd4M(Val &val, Vec &vec) {
    T ret = 0;

    boost::hana::for_each(

        boost::hana::make_range(boost::hana::int_c<0>, boost::hana::int_c<DIM>),

        [&](auto b) {
          if (a != b) {
            ret +=
                F<a, b>(val) * d2S<a, b, i, j, k, l, m, n>(val, vec) +
                2 * dFdN<a, b, m, n>(val, vec) * S<a, b, i, j, k, l>(val, vec);
          }
        }

    );

    return ret;
  }

  template <int NB, int i, int j>
  static inline T reconstructMatrix(Val &val, Vec &vec, Fun f) {
    T ret = 0;

    boost::hana::for_each(

        boost::hana::make_range(boost::hana::int_c<0>, boost::hana::int_c<NB>),

        [&](auto a) { ret += M<a, i, j>(vec) * f(L<a>(val)); }

    );

    return ret;
  }

  template <int NB, int i, int j, int k, int l>
  static inline T firstMatrixDirective(Val &val, Vec &vec, Fun f, Fun d_f) {
    T ret = 0;

    boost::hana::for_each(

        boost::hana::make_range(boost::hana::int_c<0>, boost::hana::int_c<NB>),

        [&](auto a) {
          ret +=
              M<a, i, j>(vec) * M<a, k, l>(vec) * d_f(L<a>(val)) +
              d2M<a, i, j, k, l>(val, vec) * f(L<a>(val)) / static_cast<T>(2);
        }

    );

    return ret;
  }

  template <int NB, int i, int j, int k, int l, int m, int n>
  static inline T secondMatrixDirective(Val &val, Vec &vec, Fun f, Fun d_f,
                                        Fun dd_f) {
    T ret = 0;

    boost::hana::for_each(

        boost::hana::make_range(boost::hana::int_c<0>, boost::hana::int_c<NB>),

        [&](auto a) {
          ret +=

              d2M<a, i, j, m, n>(vec) * M<a, k, l>(vec) * d_f(L<a>(val)) /
                  static_cast<T>(2) +

              M<a, i, j>(vec) * d2M<a, k, l, m, n>(vec) * d_f(L<a>(val)) /
                  static_cast<T>(2) +

              M<a, i, j>(vec) * M<a, k, l>(vec) * dd_f(L<a>(val)) +

              dd4M<a, i, j, k, l, m, n>(val, vec) * f(L<a>(val)) /
                  static_cast<T>(4) +

              d2M<a, i, j, k, l>(val, vec) * d_f(L<a>(val)) / static_cast<T>(2);
        }

    );

    return ret;
  }

private:
  template <int a, int b, int i, int j>
  static inline T dFdNa(Val &val, Vec &vec) {
    return -M<a, a, i, j>(vec) /
           ((L<a>(val) - L<b>(val)) * (L<a>(val) - L<b>(val)));
  }

  template <int a, int b, int i, int j>
  static inline T dFdNb(Val &val, Vec &vec) {
    return M<b, b, i, j>(vec) /
           ((L<a>(val) - L<b>(val)) * (L<a>(val) - L<b>(val)));
  }

  template <int a, int b, int i, int j, int k, int l>
  static inline T S(Vec &vec) {
    return G<a, b, i, j, k, l>(vec) + G<b, a, i, j, k, l>(vec);
  }

  template <int a, int b, int i, int j, int k, int l, int m, int n>
  static inline T d2G(Val &val, Vec &vec) {
    return d2M<a, i, k, n, m>(val, vec) * M<b, j, l>(vec) +
           M<a, i, k>(vec) * d2M<b, j, l, m, n>(val, vec) +
           d2M<a, i, l, m, n>(val, vec) * M<b, j, k>(val) +
           M<a, i, l>(vec) * d2M<b, j, k, m, n>(val, vec);
  }

  template <int a, int b, int i, int j, int k, int l, int m, int n>
  static inline T d2S(Val &val, Vec &vec) {
    return d2G<a, b, i, j, k, l, m, n>(val, vec) +
           d2G<b, a, i, j, k, l, m, n>(val, vec);
  }
};

} // namespace MatrixFunction
