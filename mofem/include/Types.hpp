/** \file Types.hpp
 * \brief Types
 */

/* MIT License
 *
 * Copyright (c) 2022
 * 
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 * 
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 * 
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
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
typedef char FieldBitNumber;         ///< Field bit number

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
template <typename T>
using VecAllocator = std::vector<T, std::allocator<T>>;

using IntAllocator = VecAllocator<int>;
using DoubleAllocator = VecAllocator<double>;
using ComplexDoubleAllocator = VecAllocator<std::complex<double>>;

template <typename T>
using UBlasVector = ublas::vector<T, VecAllocator<T>>;
using VectorInt = UBlasVector<int>;
using VectorDouble = UBlasVector<double>;
using VectorComplexDouble = UBlasVector<std::complex<double>>;
#ifdef WITH_ADOL_C
using VectorADouble = UBlasVector<adouble>;
#endif

template <typename T>
using UBlasMatrix = ublas::matrix<T, ublas::row_major, VecAllocator<T>>;
using MatrixInt = UBlasMatrix<int>;
using MatrixDouble = UBlasMatrix<double>;
using MatrixComplexDouble = UBlasMatrix<std::complex<double>>;
#ifdef WITH_ADOL_C
using MatrixADouble = UBlasMatrix<adouble>;
#endif

// bounded vector & matrices
template <typename T, size_t N>
using VectorBoundedArray = ublas::vector<T, ublas::bounded_array<T, N>>;

using VectorInt3 = VectorBoundedArray<int, 3>;
using VectorInt4 = VectorBoundedArray<int, 4>;
using VectorInt5 = VectorBoundedArray<int, 5>;
using VectorInt6 = VectorBoundedArray<int, 6>;
using VectorInt9 = VectorBoundedArray<int, 9>;
using VectorDouble3 = VectorBoundedArray<double, 3>;
using VectorDouble4 = VectorBoundedArray<double, 4>;
using VectorDouble5 = VectorBoundedArray<double, 5>;
using VectorDouble6 = VectorBoundedArray<double, 6>;
using VectorDouble9 = VectorBoundedArray<double, 9>;
using VectorDouble12 = VectorBoundedArray<double, 12>;
#ifdef WITH_ADOL_C
using VectorADouble9 = VectorBoundedArray<adouble, 9>;
#endif

template <typename T, size_t N>
using MatrixBoundedArray =
    ublas::matrix<T, ublas::row_major, ublas::bounded_array<T, N>>;
using MatrixDouble3by3 = MatrixBoundedArray<double, 9>;
using MatrixComplexDouble3by3 = MatrixBoundedArray<std::complex<double>, 9>;
#ifdef WITH_ADOL_C
using MatrixADouble3by3 = MatrixBoundedArray<adouble, 9>;
#endif

// shallow adaptor classes
template <typename T>
using VectorShallowArrayAdaptor =
    ublas::vector<T, ublas::shallow_array_adaptor<T>>;
using VectorAdaptor = VectorShallowArrayAdaptor<double>;
using VectorIntAdaptor = VectorShallowArrayAdaptor<int>;

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

} // namespace Types

using namespace Types;

} // namespace MoFEM

#endif // __TYPES_HPP__