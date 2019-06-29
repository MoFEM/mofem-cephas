/** \file Templates.hpp
 * \brief Templates declarations
 *
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
 * License along with MoFEM. If not, see <http://www.gnu.org/licenses/>
 */

#ifndef __TEMPLATES_HPP__
#define __TEMPLATES_HPP__

namespace MoFEM {

/**
 * @brief Get Vector adaptor
 *
 * \code
 *
 * double *a;
 * CHKERR VecGetArray(v,&a);
 *
 * for(int n = 0; n != nodes; ++n) {
 *
 *   auto a = getVectorAdaptor(&a[3*n], 3);
 *   double dot = inner_prod(a, a);
 *
 * }
 *
 * CHKERR VecRetsoreArray(v,&a);
 * \endcode
 *
 */
template <typename T1>
inline auto getVectorAdaptor(T1 ptr, const size_t n) {
  typedef typename std::remove_pointer<T1>::type T;
  return VectorShallowArrayAdaptor<T>(n,
                                      ublas::shallow_array_adaptor<T>(n, ptr));
};

/**
 * @brief Get Matrix adaptor
 *
 * \code
 *
 * double *a;
 * CHKERR VecGetArray(v,&a);
 *
 * for(int n = 0; n != nodes; ++n) {
 *
 *   auto F = getMatrixAdaptor(&a[3*3*n], 3, 3);
 *   MatrixDouble C = prod(F, trans(F));
 *
 * }
 *
 * CHKERR VecRetsoreArray(v,&a);
 * \endcode
 *
 */
template <typename T1>
inline auto getMatrixAdaptor(T1 ptr, const size_t n, const size_t m) {
  typedef typename std::remove_pointer<T1>::type T;
  return MatrixShallowArrayAdaptor<T>(
      n, m, ublas::shallow_array_adaptor<T>(n * m, ptr));
};

/**
 * This small utility that cascades two key extractors will be
 * used throughout the boost example
 * <a
 * href=http://www.boost.org/doc/libs/1_53_0/libs/multi_index/example/complex_structs.cpp>
 * http://www.boost.org/doc/libs/1_53_0/libs/multi_index/example/complex_structs.cpp
 * </a>
 */
template <class KeyExtractor1, class KeyExtractor2> struct KeyFromKey {
public:
  typedef typename KeyExtractor1::result_type result_type;

  KeyFromKey(const KeyExtractor1 &key1_ = KeyExtractor1(),
             const KeyExtractor2 &key2_ = KeyExtractor2())
      : key1(key1_), key2(key2_) {}

  template <typename Arg> result_type operator()(Arg &arg) const {
    return key1(key2(arg));
  }

private:
  KeyExtractor1 key1;
  KeyExtractor2 key2;
};

template <typename id_type> struct LtBit {
  inline bool operator()(const id_type &valueA, const id_type &valueB) const {
    return valueA.to_ulong() < valueB.to_ulong();
  }
};

template <typename id_type> struct EqBit {
  inline bool operator()(const id_type &valueA, const id_type &valueB) const {
    return valueA.to_ulong() == valueB.to_ulong();
  }
};

template <typename id_type> struct HashBit {
  inline unsigned int operator()(const id_type &value) const {
    return value.to_ulong();
  }
};

template <class X> inline std::string toString(X x) {
  std::ostringstream buffer;
  buffer << x;
  return buffer.str();
}

/**
* \brief Get tensor rank 0 (scalar) form data vector
* \ingroup mofem_forces_and_sources_user_data_operators

Example how to use it.
\code
VectorDouble vec;
vec.resize(nb_gauss_pts,false);
vec.clear();
auto t0 = getFTensor0FromData(data);
for(int gg = 0;gg!=nb_gauss_pts;gg++) {

  ++t0;
}
\endcode

*/
template <class T, class A>
static inline FTensor::Tensor0<FTensor::PackPtr<double *, 1>>
getFTensor0FromVec(ublas::vector<T, A> &data) {
  static_assert(!std::is_same<T, T>::value, "not implemented");
}

template <>
inline FTensor::Tensor0<FTensor::PackPtr<double *, 1>>
getFTensor0FromVec<double, DoubleAllocator>(
    ublas::vector<double, DoubleAllocator> &data) {
  return FTensor::Tensor0<FTensor::PackPtr<double *, 1>>(&*data.data().begin());
}

/**
 * \brief Get tensor rank 1 (vector) form data matrix
 * \ingroup mofem_forces_and_sources_user_data_operators
 */
template <int Tensor_Dim, class T, class L, class A>
static inline FTensor::Tensor1<FTensor::PackPtr<T *, 1>, Tensor_Dim>
getFTensor1FromMat(ublas::matrix<T, L, A> &data) {
  static_assert(!std::is_same<T, T>::value, "not implemented");
}

/**
 * \brief Get tensor rank 1 (vector) form data matrix (specialization)
 * \ingroup mofem_forces_and_sources_user_data_operators
 */
template <int Tensor_Dim>
static inline FTensor::Tensor1<FTensor::PackPtr<double *, 1>, Tensor_Dim>
getFTensor1FromMat(MatrixDouble &data) {
  return getFTensor1FromMat<Tensor_Dim, double, ublas::row_major,
                            DoubleAllocator>(data);
}

template <>
inline FTensor::Tensor1<FTensor::PackPtr<double *, 1>, 3>
getFTensor1FromMat<3, double, ublas::row_major, DoubleAllocator>(
    MatrixDouble &data) {
  if (data.size1() != 3) {
    THROW_MESSAGE("Wrong size of data matrix");
  }
  return FTensor::Tensor1<FTensor::PackPtr<double *, 1>, 3>(
      &data(0, 0), &data(1, 0), &data(2, 0));
}

template <>
inline FTensor::Tensor1<FTensor::PackPtr<double *, 1>, 2>
getFTensor1FromMat<2, double, ublas::row_major, DoubleAllocator>(
    MatrixDouble &data) {
  if (data.size1() != 2) {
    THROW_MESSAGE("Wrong size of data matrix");
  }
  return FTensor::Tensor1<FTensor::PackPtr<double *, 1>, 2>(&data(0, 0),
                                                            &data(1, 0));
}

// /**
//  * @deprecated Name change to getFTensor1FromMat
//  */
// template <int Tensor_Dim>
// DEPRECATED static inline FTensor::Tensor1<FTensor::PackPtr<double *, 1>,
//                                           Tensor_Dim>
// getTensor1FormData(MatrixDouble &data) {
//   return getFTensor1FromMat<Tensor_Dim>(data);
// }

/**
 * \brief Get tensor rank 2 (matrix) form data matrix
 * \ingroup mofem_forces_and_sources_user_data_operators
 */
template <int Tensor_Dim0, int Tensor_Dim1, class T, class L, class A>
static inline FTensor::Tensor2<FTensor::PackPtr<T *, 1>, Tensor_Dim0,
                               Tensor_Dim1>
getFTensor2FromMat(ublas::matrix<T, L, A> &data) {
  static_assert(!std::is_same<T, T>::value, "not implemented");
}

/**
 * Template specialization for getFTensor2FromMat
 * \ingroup mofem_forces_and_sources_user_data_operators
 *
 */
template <>
inline FTensor::Tensor2<FTensor::PackPtr<double *, 1>, 3, 3>
getFTensor2FromMat(MatrixDouble &data) {
  if (data.size1() != 9) {
    THROW_MESSAGE("Wrong size of data matrix; numer of rows is " +
                  boost::lexical_cast<std::string>(data.size1()));
  }
  return FTensor::Tensor2<FTensor::PackPtr<double *, 1>, 3, 3>(
      &data(0, 0), &data(1, 0), &data(2, 0), &data(3, 0), &data(4, 0),
      &data(5, 0), &data(6, 0), &data(7, 0), &data(8, 0));
}

/**
 * Template specialization for getFTensor2FromMat
 */
template <>
inline FTensor::Tensor2<FTensor::PackPtr<double *, 1>, 3, 2>
getFTensor2FromMat(MatrixDouble &data) {
  if (data.size1() != 6) {
    THROW_MESSAGE("Wrong size of data matrix");
  }
  return FTensor::Tensor2<FTensor::PackPtr<double *, 1>, 3, 2>(
      &data(0, 0), &data(1, 0), &data(2, 0), &data(3, 0), &data(4, 0),
      &data(5, 0));
}

/**
 * \brief Get tensor rank 2 (matrix) form data matrix (specialization)
 * \ingroup mofem_forces_and_sources_user_data_operators
 */
template <int Tensor_Dim0, int Tensor_Dim1>
static inline FTensor::Tensor2<FTensor::PackPtr<double *, 1>, Tensor_Dim0,
                               Tensor_Dim1>
getFTensor2FromMat(MatrixDouble &data) {
  return getFTensor2FromMat<Tensor_Dim0, Tensor_Dim1, double, ublas::row_major,
                            DoubleAllocator>(data);
}

// /**
//  * @deprecated Name change to getFTensor1FromMat
//  */
// template <int Tensor_Dim0, int Tensor_Dim1>
// static inline DEPRECATED
//     FTensor::Tensor2<FTensor::PackPtr<double *, 1>, Tensor_Dim0, Tensor_Dim1>
//     getTensor2FormData(MatrixDouble &data) {
//   return getFTensor2FromMat<Tensor_Dim0, Tensor_Dim1>(data);
// }

/**
 * \brief Get symmetric tensor rank 2 (matrix) form data matrix
 * \ingroup mofem_forces_and_sources_user_data_operators
 */
template <int Tensor_Dim, class T, class L, class A>
static inline FTensor::Tensor2_symmetric<FTensor::PackPtr<T *, 1>, Tensor_Dim>
getFTensor2SymmetricFromMat(ublas::matrix<T, L, A> &data) {
  static_assert(!std::is_same<T, T>::value, "not implemented");
}

/**
 * @brief Get symmetric tensor rank 2 form matrix of for dimension 3
 *
 * Specialisation for symmetric tensor 2
 *
 * @tparam
 * @param data
 * @return FTensor::Tensor2_symmetric<FTensor::PackPtr<double *, 1>, 3>
 */
template <>
inline FTensor::Tensor2_symmetric<FTensor::PackPtr<double *, 1>, 3>
getFTensor2SymmetricFromMat(MatrixDouble &data) {
  if (data.size1() != 6) {
    THROW_MESSAGE("Wrong size of data matrix");
  }
  return FTensor::Tensor2_symmetric<FTensor::PackPtr<double *, 1>, 3>(
      &data(0, 0), &data(1, 0), &data(2, 0), &data(3, 0), &data(4, 0),
      &data(5, 0));
}

/**
 * @brief Get symmetric tensor rank 2 form matrix
 *
 * Specialisation for symmetric tensor 2
 *
 * @tparam Tensor_Dim
 * @param data
 * @return FTensor::Tensor2_symmetric<FTensor::PackPtr<double *, 1>, Tensor_Dim>
 */
template <int Tensor_Dim>
static inline FTensor::Tensor2_symmetric<FTensor::PackPtr<double *, 1>,
                                         Tensor_Dim>
getFTensor2SymmetricFromMat(MatrixDouble &data) {
  return getFTensor2SymmetricFromMat<Tensor_Dim, double, ublas::row_major,
                                     DoubleAllocator>(data);
}

/**
 * @brief Calculate the determinant of a 3x3 matrix or a tensor of rank 2
 *
 * @tparam T
 * @param t
 * @return double
 */
template <class T> static inline double dEterminant(T &t) {
  return t(0, 0) * t(1, 1) * t(2, 2) + t(1, 0) * t(2, 1) * t(0, 2) +
         t(2, 0) * t(0, 1) * t(1, 2) - t(0, 0) * t(2, 1) * t(1, 2) -
         t(2, 0) * t(1, 1) * t(0, 2) - t(1, 0) * t(0, 1) * t(2, 2);
}

/**
 * \brief Calculate inverse of tensor rank 2 at integration points

 * \ingroup mofem_forces_and_sources
 */
template <int Tensor_Dim, class T, class L, class A>
inline MoFEMErrorCode invertTensor3by3(ublas::matrix<T, L, A> &jac_data,
                                       ublas::vector<T, A> &det_data,
                                       ublas::matrix<T, L, A> &inv_jac_data) {
  MoFEMFunctionBeginHot;
  SETERRQ(PETSC_COMM_SELF, MOFEM_NOT_IMPLEMENTED,
          "Specialization for this template not yet implemented");
  MoFEMFunctionReturnHot(0);
}

template <>
inline MoFEMErrorCode
invertTensor3by3<3, double, ublas::row_major, DoubleAllocator>(
    MatrixDouble &jac_data, VectorDouble &det_data, MatrixDouble &inv_jac_data);

/**
 * \brief Calculate determinant

 * \ingroup mofem_forces_and_sources
 */
template <class T1, class T2>
inline MoFEMErrorCode determinantTensor3by3(T1 &t, T2 &det) {
  MoFEMFunctionBeginHot;
  det = +t(0, 0) * t(1, 1) * t(2, 2) + t(1, 0) * t(2, 1) * t(0, 2) +
        t(2, 0) * t(0, 1) * t(1, 2) - t(0, 0) * t(2, 1) * t(1, 2) -
        t(2, 0) * t(1, 1) * t(0, 2) - t(1, 0) * t(0, 1) * t(2, 2);
  MoFEMFunctionReturnHot(0);
}

/**
 * \brief Calculate matrix inverse

 * \ingroup mofem_forces_and_sources
 */
template <class T1, class T2, class T3>
inline MoFEMErrorCode invertTensor3by3(T1 &t, T2 &det, T3 &inv_t) {
  MoFEMFunctionBeginHot;
  const auto inv_det = 1. / det;
  inv_t(0, 0) = (t(1, 1) * t(2, 2) - t(1, 2) * t(2, 1)) * inv_det;
  inv_t(0, 1) = (t(0, 2) * t(2, 1) - t(0, 1) * t(2, 2)) * inv_det;
  inv_t(0, 2) = (t(0, 1) * t(1, 2) - t(0, 2) * t(1, 1)) * inv_det;
  inv_t(1, 0) = (t(1, 2) * t(2, 0) - t(1, 0) * t(2, 2)) * inv_det;
  inv_t(1, 1) = (t(0, 0) * t(2, 2) - t(0, 2) * t(2, 0)) * inv_det;
  inv_t(1, 2) = (t(0, 2) * t(1, 0) - t(0, 0) * t(1, 2)) * inv_det;
  inv_t(2, 0) = (t(1, 0) * t(2, 1) - t(1, 1) * t(2, 0)) * inv_det;
  inv_t(2, 1) = (t(0, 1) * t(2, 0) - t(0, 0) * t(2, 1)) * inv_det;
  inv_t(2, 2) = (t(0, 0) * t(1, 1) - t(0, 1) * t(1, 0)) * inv_det;
  MoFEMFunctionReturnHot(0);
}

#ifdef WITH_ADOL_C

/**
 * \brief Calculate matrix inverse, specialization for adouble tensor

 * \ingroup mofem_forces_and_sources
 */
template <>
inline MoFEMErrorCode invertTensor3by3<FTensor::Tensor2<adouble, 3, 3>, adouble,
                                       FTensor::Tensor2<adouble, 3, 3>>(
    FTensor::Tensor2<adouble, 3, 3> &t, adouble &det,
    FTensor::Tensor2<adouble, 3, 3> &inv_t) {
  MoFEMFunctionBeginHot;
  inv_t(0, 0) = (t(1, 1) * t(2, 2) - t(1, 2) * t(2, 1)) / det;
  inv_t(0, 1) = (t(0, 2) * t(2, 1) - t(0, 1) * t(2, 2)) / det;
  inv_t(0, 2) = (t(0, 1) * t(1, 2) - t(0, 2) * t(1, 1)) / det;
  inv_t(1, 0) = (t(1, 2) * t(2, 0) - t(1, 0) * t(2, 2)) / det;
  inv_t(1, 1) = (t(0, 0) * t(2, 2) - t(0, 2) * t(2, 0)) / det;
  inv_t(1, 2) = (t(0, 2) * t(1, 0) - t(0, 0) * t(1, 2)) / det;
  inv_t(2, 0) = (t(1, 0) * t(2, 1) - t(1, 1) * t(2, 0)) / det;
  inv_t(2, 1) = (t(0, 1) * t(2, 0) - t(0, 0) * t(2, 1)) / det;
  inv_t(2, 2) = (t(0, 0) * t(1, 1) - t(0, 1) * t(1, 0)) / det;
  MoFEMFunctionReturnHot(0);
}

#endif

/**
 * \brief Calculate matrix inverse, specialization for symmetric tensor

 * \ingroup mofem_forces_and_sources
 */
template <>
inline MoFEMErrorCode
invertTensor3by3<FTensor::Tensor2_symmetric<double, 3>, double,
                 FTensor::Tensor2_symmetric<double, 3>>(
    FTensor::Tensor2_symmetric<double, 3> &t, double &det,
    FTensor::Tensor2_symmetric<double, 3> &inv_t) {
  MoFEMFunctionBeginHot;
  const auto inv_det = 1. / det;
  inv_t(0, 0) = (t(1, 1) * t(2, 2) - t(1, 2) * t(2, 1)) * inv_det;
  inv_t(0, 1) = (t(0, 2) * t(2, 1) - t(0, 1) * t(2, 2)) * inv_det;
  inv_t(0, 2) = (t(0, 1) * t(1, 2) - t(0, 2) * t(1, 1)) * inv_det;
  inv_t(1, 1) = (t(0, 0) * t(2, 2) - t(0, 2) * t(2, 0)) * inv_det;
  inv_t(1, 2) = (t(0, 2) * t(1, 0) - t(0, 0) * t(1, 2)) * inv_det;
  inv_t(2, 2) = (t(0, 0) * t(1, 1) - t(0, 1) * t(1, 0)) * inv_det;
  MoFEMFunctionReturnHot(0);
}

#ifdef WITH_ADOL_C

/**
 * \brief Calculate matrix inverse, specialization for adouble symmetric tensor

 * \ingroup mofem_forces_and_sources
 */
template <>
inline MoFEMErrorCode
invertTensor3by3<FTensor::Tensor2_symmetric<adouble, 3>, adouble,
                 FTensor::Tensor2_symmetric<adouble, 3>>(
    FTensor::Tensor2_symmetric<adouble, 3> &t, adouble &det,
    FTensor::Tensor2_symmetric<adouble, 3> &inv_t) {
  MoFEMFunctionBeginHot;
  inv_t(0, 0) = (t(1, 1) * t(2, 2) - t(1, 2) * t(2, 1)) / det;
  inv_t(0, 1) = (t(0, 2) * t(2, 1) - t(0, 1) * t(2, 2)) / det;
  inv_t(0, 2) = (t(0, 1) * t(1, 2) - t(0, 2) * t(1, 1)) / det;
  inv_t(1, 1) = (t(0, 0) * t(2, 2) - t(0, 2) * t(2, 0)) / det;
  inv_t(1, 2) = (t(0, 2) * t(1, 0) - t(0, 0) * t(1, 2)) / det;
  inv_t(2, 2) = (t(0, 0) * t(1, 1) - t(0, 1) * t(1, 0)) / det;
  MoFEMFunctionReturnHot(0);
}

#endif

/**
 * \brief Calculate matrix inverse, specialization for symmetric (pointer)
 tensor

 * \ingroup mofem_forces_and_sources
 */
template <>
inline MoFEMErrorCode
invertTensor3by3<FTensor::Tensor2_symmetric<double, 3>, double,
                 FTensor::Tensor2_symmetric<double *, 3>>(
    FTensor::Tensor2_symmetric<double, 3> &t, double &det,
    FTensor::Tensor2_symmetric<double *, 3> &inv_t) {
  MoFEMFunctionBeginHot;
  const auto inv_det = 1. / det;
  inv_t(0, 0) = (t(1, 1) * t(2, 2) - t(1, 2) * t(2, 1)) * inv_det;
  inv_t(0, 1) = (t(0, 2) * t(2, 1) - t(0, 1) * t(2, 2)) * inv_det;
  inv_t(0, 2) = (t(0, 1) * t(1, 2) - t(0, 2) * t(1, 1)) * inv_det;
  inv_t(1, 1) = (t(0, 0) * t(2, 2) - t(0, 2) * t(2, 0)) * inv_det;
  inv_t(1, 2) = (t(0, 2) * t(1, 0) - t(0, 0) * t(1, 2)) * inv_det;
  inv_t(2, 2) = (t(0, 0) * t(1, 1) - t(0, 1) * t(1, 0)) * inv_det;
  MoFEMFunctionReturnHot(0);
}

} // namespace MoFEM

#endif //__TEMPLATES_HPP__