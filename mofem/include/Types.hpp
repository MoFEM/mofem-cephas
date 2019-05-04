/** \file Types.hpp
 * \brief Types 
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

#ifndef __TYPES_HPP__
#define __TYPES_HPP__

namespace MoFEM {

/**
 * @brief Types 
 * 
 */
namespace Types {

typedef int DofIdx;                  ///< Index of DOF
typedef int MoFEMDofIdx;             ///< Index of DOF using mofem native index
typedef int PetscLocalDofIdx;        ///< Index of DOF using local petsc index
typedef int PetscGlobalDofIdx;       ///< Index of DOF using global pets index
typedef int FEIdx;                   ///< Index of the element
typedef int EntIdx;                  ///< Index of DOF on the entity
typedef int EntPart;                 ///< Partition owning entity
typedef double FieldData;            ///< Field data type
typedef int ApproximationOrder;      ///< Approximation on the entity
typedef int FieldCoefficientsNumber; ///< Number of field coefficients

// typedef checked_uint128_t UId;
typedef uint128_t UId; ///< Unique Id
typedef int ShortId;   ///< Unique Id in the field


typedef std::bitset<BITREFEDGES_SIZE> BitRefEdges;

/**
 * \brief Bit structure attached to each entity identifying to what mesh entity
 * is attached.
 */
typedef std::bitset<BITREFLEVEL_SIZE> BitRefLevel;

typedef std::bitset<BITFIELDID_SIZE> BitFieldId;     ///< Field Id
typedef std::bitset<BITFEID_SIZE> BitFEId;           ///< Finite element Id
typedef std::bitset<BITPROBLEMID_SIZE> BitProblemId; ///< Problem Id
typedef std::bitset<BITINTERFACEUID_SIZE> BitIntefaceId;

/**
 * \typedef CubitBCType
 * bc & material meshsets
 *
 */
typedef std::bitset<32> CubitBCType;

// array with std allocators (i.e. concept of capacity is useful here)
// typedef ublas::unbounded_array<int,std::allocator<int> > IntAllocator;
// typedef ublas::unbounded_array<double,std::allocator<double> >
// DoubleAllocator;
typedef std::vector<int, std::allocator<int>> IntAllocator;
typedef std::vector<double, std::allocator<double>> DoubleAllocator;
// DEPRECATED typedef IntAllocator
//     IntAllacator; ///< \deprecated Do not use spelling mistake
// DEPRECATED typedef DoubleAllocator
//     DoubleAllacator; ///< \deprecated Do not use spelling mistake
typedef ublas::vector<int, IntAllocator> VectorInt;
typedef ublas::vector<double, DoubleAllocator> VectorDouble;
typedef ublas::matrix<double, ublas::row_major, DoubleAllocator> MatrixDouble;

// bounded vector & matrices
template <typename T, size_t N>
using VectorBoundedArray = ublas::vector<T, ublas::bounded_array<T, N>>;

typedef VectorBoundedArray<int, 3> VectorInt3;
typedef VectorBoundedArray<int, 4> VectorInt4;
typedef VectorBoundedArray<int, 5> VectorInt5;
typedef VectorBoundedArray<int, 6> VectorInt6;
typedef VectorBoundedArray<int, 9> VectorInt9;
typedef VectorBoundedArray<double, 3> VectorDouble3;
typedef VectorBoundedArray<double, 4> VectorDouble4;
typedef VectorBoundedArray<double, 5> VectorDouble5;
typedef VectorBoundedArray<double, 6> VectorDouble6;
typedef VectorBoundedArray<double, 9> VectorDouble9;
typedef VectorBoundedArray<double, 12> VectorDouble12;

template <typename T, size_t N>
using MatrixBoundedArray =
    ublas::matrix<T, ublas::row_major, ublas::bounded_array<T, N>>;
typedef MatrixBoundedArray<double, 9> MatrixDouble3by3;

// shallow adaptor classes
template <typename T>
using VectorShallowArrayAdaptor =
    ublas::vector<T, ublas::shallow_array_adaptor<T>>;
typedef VectorShallowArrayAdaptor<double> VectorAdaptor;
typedef VectorShallowArrayAdaptor<int> VectorIntAdaptor;

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
auto getVectorAdaptor = [](auto ptr, const int n) {
  typedef typename std::remove_pointer<decltype(ptr)>::type T;
  return VectorShallowArrayAdaptor<T>(n,
                                      ublas::shallow_array_adaptor<T>(n, ptr));
};

template <typename T>
using MatrixShallowArrayAdaptor =
    ublas::matrix<double, ublas::row_major,
                  ublas::shallow_array_adaptor<double>>;

/**
 * @brief Matrix adaptor.
 *
 * \code
 * MatrixAdaptor mat = MatrixAdaptor(3, 3,
 *    ublas::shallow_array_adaptor<double>(9, ptr));
 * \endcode
 *
 */
typedef MatrixShallowArrayAdaptor<double> MatrixAdaptor;

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

// /**
//  * @deprecated Name change to getFTensor0FromVec
//  */
// template <class T, class A>
// DEPRECATED static inline FTensor::Tensor0<FTensor::PackPtr<double *, 1>>
// getTensor0FormData(ublas::vector<T, A> &data) {
//   return getFTensor0FromVec(data);
// }

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


} // namespace Types

using namespace Types;

} // namespace MoFEM

#endif // __TYPES_HPP__