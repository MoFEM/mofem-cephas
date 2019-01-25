/** \file Common.hpp
 * \brief Basic structures and data
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

#ifndef __COMMON_HPP__
#define __COMMON_HPP__

namespace MoFEM {

/**
 * @brief No entity handle is indicated by zero handle, i.e. root meshset
 *
 */
const EntityHandle no_handle = 0;

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
DEPRECATED typedef IntAllocator
    IntAllacator; ///< \deprecated Do not use spelling mistake
DEPRECATED typedef DoubleAllocator
    DoubleAllacator; ///< \deprecated Do not use spelling mistake
typedef ublas::vector<int, IntAllocator> VectorInt;
typedef ublas::vector<double, DoubleAllocator> VectorDouble;
typedef ublas::matrix<double, ublas::row_major, DoubleAllocator> MatrixDouble;

// bounded vector & matrices
template <typename T, size_t N>
using VectorBoundedArray = ublas::vector<T, ublas::bounded_array<T, N>>;

typedef VectorBoundedArray<int, 3> VectorInt3;
typedef VectorBoundedArray<int, 9> VectorInt9;
typedef VectorBoundedArray<double, 3> VectorDouble3;
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

template <class X> inline std::string toString(X x) {
  std::ostringstream buffer;
  buffer << x;
  return buffer.str();
}
} // namespace MoFEM

#endif //__COMMON_HPP__

/***************************************************************************/ /**
                                                                               * \defgroup mofem MoFEM
                                                                               ******************************************************************************/
