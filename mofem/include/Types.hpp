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
typedef std::vector<int, std::allocator<int>> IntAllocator;
typedef std::vector<double, std::allocator<double>> DoubleAllocator;
typedef std::vector<std::complex<double>, std::allocator<std::complex<double>>>
    ComplexDoubleAllocator;
typedef ublas::vector<int, IntAllocator> VectorInt;
typedef ublas::vector<double, DoubleAllocator> VectorDouble;
typedef ublas::matrix<int, ublas::row_major, IntAllocator> MatrixInt;
typedef ublas::matrix<double, ublas::row_major, DoubleAllocator> MatrixDouble;
typedef ublas::vector<std::complex<double>, ComplexDoubleAllocator>
    VectorComplexDouble;
typedef ublas::matrix<std::complex<double>, ublas::row_major,
                      ComplexDoubleAllocator>
    MatrixComplexDouble;

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
typedef MatrixBoundedArray<std::complex<double>, 9> MatrixComplexDouble3by3;

// shallow adaptor classes
template <typename T>
using VectorShallowArrayAdaptor =
    ublas::vector<T, ublas::shallow_array_adaptor<T>>;
typedef VectorShallowArrayAdaptor<double> VectorAdaptor;
typedef VectorShallowArrayAdaptor<int> VectorIntAdaptor;

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