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
template <typename T1> inline auto getVectorAdaptor(T1 ptr, const size_t n) {
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
  return FTensor::Tensor0<FTensor::PackPtr<double *, 1>>(nullptr);
}

template <>
inline FTensor::Tensor0<FTensor::PackPtr<double *, 1>>
getFTensor0FromVec<double, DoubleAllocator>(
    ublas::vector<double, DoubleAllocator> &data) {
  return FTensor::Tensor0<FTensor::PackPtr<double *, 1>>(&*data.data().begin());
}

/**
 * \brief Get tensor rank 1 (vector) form data matrix
 */
template <int Tensor_Dim, class T, class L, class A>
static inline FTensor::Tensor1<FTensor::PackPtr<T *, 1>, Tensor_Dim>
getFTensor1FromMat(ublas::matrix<T, L, A> &data) {
  static_assert(!std::is_same<T, T>::value, "not implemented");
}

/**
 * \brief Get tensor rank 1 (vector) form data matrix (specialization)
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
  if (data.size1() != 3)
    THROW_MESSAGE("getFTensor1FromMat<3>: wrong size of data matrix, number of "
                  "rows should be 3 but is %d" +
                  boost::lexical_cast<std::string>(data.size1()));

  return FTensor::Tensor1<FTensor::PackPtr<double *, 1>, 3>(
      &data(0, 0), &data(1, 0), &data(2, 0));
}

template <>
inline FTensor::Tensor1<FTensor::PackPtr<double *, 1>, 2>
getFTensor1FromMat<2, double, ublas::row_major, DoubleAllocator>(
    MatrixDouble &data) {
  if (data.size1() != 2)
    THROW_MESSAGE("getFTensor1FromMat<2>: wrong size of data matrix, number of "
                  "rows should be 2 but is %d" +
                  boost::lexical_cast<std::string>(data.size1()));

  return FTensor::Tensor1<FTensor::PackPtr<double *, 1>, 2>(&data(0, 0),
                                                            &data(1, 0));
}

/**
 * \brief Get tensor rank 2 (matrix) form data matrix
 */
template <int Tensor_Dim0, int Tensor_Dim1, class T, class L, class A>
static inline FTensor::Tensor2<FTensor::PackPtr<T *, 1>, Tensor_Dim0,
                               Tensor_Dim1>
getFTensor2FromMat(ublas::matrix<T, L, A> &data) {
  static_assert(!std::is_same<T, T>::value,
                "Such getFTensor2FromMat specialisation is not implemented");
}

/**
 * Template specialization for getFTensor2FromMat
 *
 */
template <>
inline FTensor::Tensor2<FTensor::PackPtr<double *, 1>, 3, 3>
getFTensor2FromMat(MatrixDouble &data) {
  if (data.size1() != 9)
    THROW_MESSAGE("getFTensor2FromMat<3,3>: wrong size of data matrix; numer "
                  "of rows should be 9 but is " +
                  boost::lexical_cast<std::string>(data.size1()));

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
  if (data.size1() != 6)
    THROW_MESSAGE("getFTensor2FromMat<3,3>: wrong size of data matrix, numer "
                  "of rows should be 6 but is " +
                  boost::lexical_cast<std::string>(data.size1()));

  return FTensor::Tensor2<FTensor::PackPtr<double *, 1>, 3, 2>(
      &data(0, 0), &data(1, 0), &data(2, 0), &data(3, 0), &data(4, 0),
      &data(5, 0));
}

/**
 * Template specialization for getFTensor2FromMat
 */
template <>
inline FTensor::Tensor2<FTensor::PackPtr<double *, 1>, 2, 2>
getFTensor2FromMat(MatrixDouble &data) {
  if (data.size1() != 4)
    THROW_MESSAGE("getFTensor2FromMat<2,2>: wrong size of data matrix, numer "
                  "of rows should be 4 but is " +
                  boost::lexical_cast<std::string>(data.size1()));

  return FTensor::Tensor2<FTensor::PackPtr<double *, 1>, 2, 2>(
      &data(0, 0), &data(1, 0), &data(2, 0), &data(3, 0));
}

/**
 * \brief Get tensor rank 2 (matrix) form data matrix (specialization)
 */
template <int Tensor_Dim0, int Tensor_Dim1>
static inline FTensor::Tensor2<FTensor::PackPtr<double *, 1>, Tensor_Dim0,
                               Tensor_Dim1>
getFTensor2FromMat(MatrixDouble &data) {
  return getFTensor2FromMat<Tensor_Dim0, Tensor_Dim1, double, ublas::row_major,
                            DoubleAllocator>(data);
}

/**
 * \brief Get symmetric tensor rank 2 (matrix) form data matrix
 */
template <int Tensor_Dim, class T, class L, class A>
static inline FTensor::Tensor2_symmetric<FTensor::PackPtr<T *, 1>, Tensor_Dim>
getFTensor2SymmetricFromMat(ublas::matrix<T, L, A> &data) {
  static_assert(
      !std::is_same<T, T>::value,
      "Such getFTensor2SymmetricFromMat specialisation is not implemented");
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
  if (data.size1() != 6)
    THROW_MESSAGE(
        "getFTensor2SymmetricFromMat<3>: wrong size of data matrix, numer "
        "of rows should be 6 but is " +
        boost::lexical_cast<std::string>(data.size1()));

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
 * @brief Get symmetric tensor rank 4  on first two and last indices from
 * form data matrix
 *
 * @tparam Tensor_Dim01 dimension of frirst two indicies
 * @tparam Tensor_Dim23 dimension of second two indicies
 * @tparam T the type of object stored
 * @tparam L the storage organization
 * @tparam A 	the type of Storage array
 * @param data data container
 * @return FTensor::Ddg<FTensor::PackPtr<T *, 1>, Tensor_Dim01, TensorDim23>
 */
template <int Tensor_Dim01, int Tensor_Dim23, class T, class L, class A>
static inline FTensor::Ddg<FTensor::PackPtr<T *, 1>, Tensor_Dim01, Tensor_Dim23>
getFTensor4DdgFromMat(ublas::matrix<T, L, A> &data) {
  static_assert(
      !std::is_same<T, T>::value,
      "Such getFTensor4DdgFromMat specialisation is not implemented");
}

/**
 * @brief Get symmetric tensor rank 4  on first two and last indices from
 * form data matrix
 *
 * @param data matrix container which has 36 rows
 * @return FTensor::Ddg<FTensor::PackPtr<T *, 1>, Tensor_Dim01, TensorDim23>
 */
template <>
inline FTensor::Ddg<FTensor::PackPtr<double *, 1>, 3, 3>
getFTensor4DdgFromMat(MatrixDouble &data) {
  if (data.size1() != 36)
    THROW_MESSAGE(
        "getFTensor4DdgFromMat<3, 3>: wrong size of data matrix, number "
        "of rows should be 36 but is " +
        boost::lexical_cast<std::string>(data.size1()));

  return FTensor::Ddg<FTensor::PackPtr<double *, 1>, 3, 3>(
      &data(0, 0), &data(1, 0), &data(2, 0), &data(3, 0), &data(4, 0),
      &data(5, 0), &data(6, 0), &data(7, 0), &data(8, 0), &data(9, 0),
      &data(10, 0), &data(11, 0), &data(12, 0), &data(13, 0), &data(14, 0),
      &data(15, 0), &data(16, 0), &data(17, 0), &data(18, 0), &data(19, 0),
      &data(20, 0), &data(21, 0), &data(22, 0), &data(23, 0), &data(24, 0),
      &data(25, 0), &data(26, 0), &data(27, 0), &data(28, 0), &data(29, 0),
      &data(30, 0), &data(31, 0), &data(32, 0), &data(33, 0), &data(34, 0),
      &data(35, 0));
}

/**
 * @brief Make Tensor1 from pointer
 *
 * @tparam DIM
 * @param ptr
 * @return FTensor::Tensor2<FTensor::PackPtr<double *, 3 * DIM>, 3, DIM>
 */
template <int DIM>
inline FTensor::Tensor1<FTensor::PackPtr<double *, DIM>, DIM>
getFTensor1FromPtr(double *ptr) {
  static_assert(DIM != 3,
                "Such getFTensor1FromPtr specialization is not implemented");
};

template <>
inline FTensor::Tensor1<FTensor::PackPtr<double *, 3>, 3>
getFTensor1FromPtr<3>(double *ptr) {
  return FTensor::Tensor1<FTensor::PackPtr<double *, 3>, 3>(
      &ptr[HVEC0], &ptr[HVEC1], &ptr[HVEC2]);
};

/**
 * @brief Make Tensor2 from pointer
 *
 * @tparam DIM
 * @param ptr
 * @return FTensor::Tensor2<FTensor::PackPtr<double *, DIM1 * DIM2>, DIM1, DIM2>
 */
template <int DIM1, int DIM2>
inline FTensor::Tensor2<FTensor::PackPtr<double *, DIM1 * DIM2>, DIM1, DIM2>
getFTensor2FromPtr(double *ptr) {
  static_assert(DIM1 != 3, "Such getFTensor2FromPtr is not implemented");
  static_assert(DIM2 >= 2 && DIM2 <= 3,
                "Such getFTensor2FromPtr is not implemented");
};

template <>
FTensor::Tensor2<FTensor::PackPtr<double *, 6>, 3,
                 2> inline getFTensor2FromPtr<3, 2>(double *ptr) {
  return FTensor::Tensor2<FTensor::PackPtr<double *, 6>, 3, 2>(
      &ptr[HVEC0_0], &ptr[HVEC0_1],

      &ptr[HVEC1_0], &ptr[HVEC1_1],

      &ptr[HVEC2_0], &ptr[HVEC2_1]);
};

template <>
FTensor::Tensor2<FTensor::PackPtr<double *, 9>, 3,
                 3> inline getFTensor2FromPtr<3, 3>(double *ptr) {
  return FTensor::Tensor2<FTensor::PackPtr<double *, 9>, 3, 3>(
      &ptr[HVEC0_0], &ptr[HVEC0_1], &ptr[HVEC0_2],

      &ptr[HVEC1_0], &ptr[HVEC1_1], &ptr[HVEC1_2],

      &ptr[HVEC2_0], &ptr[HVEC2_1], &ptr[HVEC2_2]);
};

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
 * \brief Calculate determinant 3 by 3

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
 * \brief Calculate determinant 2 by 2

 */
template <class T1, class T2>
inline MoFEMErrorCode determinantTensor2by2(T1 &t, T2 &det) {
  MoFEMFunctionBeginHot;
  det = t(0, 0) * t(1, 1) - t(0, 1) * t(1, 0);
  MoFEMFunctionReturnHot(0);
}

/**
 * \brief Calculate matrix inverse 3 by 3

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

/**
 * \brief Calculate matrix inverse 2 by 2

 */
template <class T1, class T2, class T3>
inline MoFEMErrorCode invertTensor2by2(T1 &t, T2 &det, T3 &inv_t) {
  MoFEMFunctionBeginHot;
  const auto inv_det = 1. / det;
  inv_t(0, 0) = t(1, 1) * inv_det;
  inv_t(0, 1) = -t(0, 1) * inv_det;
  inv_t(1, 0) = -t(1, 0) * inv_det;
  inv_t(1, 1) = t(0, 0) * inv_det;
  MoFEMFunctionReturnHot(0);
}

#ifdef WITH_ADOL_C

/**
 * \brief Calculate matrix inverse, specialization for adouble tensor

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

/**
 * @brief Extract entity handle form multi-index container
 *
 */
struct RefEntExtractor {
  template <typename Iterator>
  static inline EntityHandle extract(const Iterator &it) {
    return (*it)->getRefEnt();
  }
};

/**
 * @brief Insert ordered mofem multi-index into range
 *
 * \code
 * auto hi_rit = refEntsPtr->upper_bound(start);
 * auto hi_rit = refEntsPtr->upper_bound(end);
 * Range to_erase;
 * insertOrdered(to_erase, RefEntExtractor(), rit, hi_rit);
 * \endcode
 *
 * @tparam Iterator
 * @param r
 * @param begin_iter
 * @param end_iter
 * @return moab::Range::iterator
 */
template <typename Extractor, typename Iterator>
moab::Range::iterator insertOrdered(Range &r, Extractor, Iterator begin_iter,
                                    Iterator end_iter) {
  moab::Range::iterator hint = r.begin();
  while (begin_iter != end_iter) {
    size_t j = 0;
    auto bi = Extractor::extract(begin_iter);
    Iterator pj = begin_iter;
    while (pj != end_iter && (bi + j) == Extractor::extract(pj)) {
      ++pj;
      ++j;
    }
    hint = r.insert(hint, bi, bi + (j - 1));
    begin_iter = pj;
  }
  return hint;
};

} // namespace MoFEM

#endif //__TEMPLATES_HPP__