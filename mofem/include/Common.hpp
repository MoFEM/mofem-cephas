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
 * \brief Exception to catch
 */
struct MoFEMException : public std::exception {
  const int errorCode;
  char errorMessage[255];
  MoFEMException(const MoFEMErrorCodes error_code) : errorCode(error_code) {
    strcpy(errorMessage, "Huston we have a problem, something is wrong");
  }
  MoFEMException(const MoFEMErrorCodes error_code, const char error_message[])
      : errorCode(error_code) {
    strcpy(errorMessage, error_message);
  }
  const char *what() const throw() { return errorMessage; }

protected:
  MoFEMException(const int error_code) : errorCode(error_code) {
    strcpy(errorMessage, "Huston we have a problem, something is wrong");
  }
};

struct MoFEMExceptionRepeat : public MoFEMException {
  const int lINE;
  MoFEMExceptionRepeat(const int error_code, const int line)
      : MoFEMException(error_code), lINE(line) {
    strcpy(errorMessage, " ");
  }
};

struct MoFEMExceptionInitial : public MoFEMExceptionRepeat {
  MoFEMExceptionInitial(const int error_code, const char error_message[],
                        const int line)
      : MoFEMExceptionRepeat(error_code, line) {
    strcpy(errorMessage, error_message);
  }
};

typedef moab::ErrorCode MoABErrorCode; ///< MoAB error code
typedef PetscErrorCode MoFEMErrorCode; ///< MoFEM/PETSc error code

template <typename TYPE> struct MoFEMErrorCodeGeneric {
  MoFEMErrorCodeGeneric(const TYPE) {}
};

template <> struct MoFEMErrorCodeGeneric<PetscErrorCode> {
  PetscErrorCode iERR;
  MoFEMErrorCodeGeneric(const PetscErrorCode ierr) : iERR(ierr) {}
  inline operator PetscErrorCode() const { return iERR; }
};

template <> struct MoFEMErrorCodeGeneric<moab::ErrorCode> {
  moab::ErrorCode rVAL;
  MoFEMErrorCodeGeneric(const moab::ErrorCode rval) : rVAL(rval) {}
  inline operator moab::ErrorCode() const { return rVAL; }
};

static MoFEMErrorCodeGeneric<moab::ErrorCode> rval =
    MoFEMErrorCodeGeneric<moab::ErrorCode>(MB_SUCCESS);
static MoFEMErrorCodeGeneric<PetscErrorCode> ierr =
    MoFEMErrorCodeGeneric<PetscErrorCode>(0);

/**
 * \brief Error check for inline function check.
 *
 * This class is not used directly, it is called in CHKERR. In case of the error
 * pass line number and that is catch at the end of the function. Information is
 * enriched by function name and file name. Then error is pushed to PETSc error
 * stack.
 *
 * \note This class has no variables and line number is set at compilation.
 * Adding variables to this function will reduce efficiency of the code. Do
 * not do that.
 *
 */
template <int LINE> struct ErrorCheckerCode {

  /**
   * @brief Operator for handling PetscErrorCode and MoFEMErrorCode
   *
   */
  inline void operator<<(const MoFEMErrorCode err) {
    if (PetscUnlikely(err)) {
      throw MoFEMExceptionRepeat(err, LINE);
    }
    return;
  }

  /**
   * @brief Operator for handling moab::ErrorCode
   *
   */
  inline void operator<<(const moab::ErrorCode err) {
    if (PetscLikely(MB_SUCCESS != err)) {
      std::string error_str = (unsigned)err <= (unsigned)MB_FAILURE
                                  ? moab::ErrorCodeStr[err]
                                  : "INVALID ERROR CODE";
      std::string str("MOAB error (" + boost::lexical_cast<std::string>(err) +
                      ") " + error_str);
      throw MoFEMExceptionInitial(MOFEM_MOAB_ERROR, str.c_str(), LINE);
    }
    return;
  }
};

typedef int DofIdx;                  ///< Index of DOF
typedef int MoFEMDofIdx;             ///< Index of DOF using mofem native index
typedef int PetscLocalDodIdx;        ///< Index of DOF using local petsc index
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

#define UID_DOF_MAK 0x1FF

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

// AUX STRUCTURES

/* This small utility that cascades two key extractors will be
 * used throughout the boost example
 * http://www.boost.org/doc/libs/1_53_0/libs/multi_index/example/complex_structs.cpp
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
