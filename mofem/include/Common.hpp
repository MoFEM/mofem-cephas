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
 * \brief Exception to catch
 */
struct MoFEMException : public std::exception {
  MoFEMErrorCodes errorCode;
  char errorMessage[255];
  MoFEMException(const MoFEMErrorCodes error_code) : errorCode(error_code) {
    strcpy(errorMessage, "Huston we have a problem, something is wrong");
  }
  MoFEMException(const MoFEMErrorCodes error_code, const char error_message[])
      : errorCode(error_code) {
    strcpy(errorMessage, error_message);
  }
  const char *what() const throw() { return errorMessage; }
};

struct MoFEMExceptionRepeat : public MoFEMException {
  const int lINE;
  const char* fILE;
  const char* fUN;
  MoFEMExceptionRepeat(const MoFEMErrorCodes error_code, const int line,
                       const char *file, const char *fun)
      : MoFEMException(error_code), lINE(line), fILE(file), fUN(fun) {
    strcpy(errorMessage, " ");
  }
  const char *what() const throw() { return errorMessage; }
};

struct MoFEMExceptionNoRepeat : public MoFEMExceptionRepeat {
  MoFEMExceptionNoRepeat(const MoFEMErrorCodes error_code, const int line,
                         const char *file, const char *fun)
      : MoFEMExceptionRepeat(error_code, line, file, fun) {
    strcpy(errorMessage, " ");
  }
  const char *what() const throw() { return errorMessage; }
};

typedef moab::ErrorCode MoABErrorCode;  ///< MoAB error code
typedef PetscErrorCode MoFEMErrorCode;  ///< MoFEM/PETSc error code

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
 * \brief Error check form inline erro check
 *
 * \note This class is not used directly, only for testing. See
 * ERROR_CHECKER_CODE, ERROR_RUNNER_CODE, ERROR_CHECKER_AND_RUNNER_CODE and
 * CHKERR for details.
 *
 * A sequence of overloaded operators << is called, starting from line number
 * file name, function name and error code. Use this using definition CHKERR,
 * for example \code CHKERR fun_moab(); CHKERR fun_petsc(); CHKERR fun_mofem();
 * \endcode
 *
 */
template <int LINE> struct ErrorCheckerCode {

/**
 * @brief Definition for error checker class to create checker operator for
 * PetscErrorCode and MoFEMErrorCode
 *
 */
#define OP_ERR_MOFEM_ERROR_CODE(LINE, FILE, FUNC)                              \
  inline void operator<<(const MoFEMErrorCode err) {                           \
    if (PetscUnlikely(err)) {                                                  \
      switch (err) {                                                           \
      case MOFEM_DATA_INCONSISTENCY:                                           \
        throw MoFEMExceptionRepeat(MOFEM_DATA_INCONSISTENCY, LINE, FILE,       \
                                   FUNC);                                      \
      case MOFEM_NOT_IMPLEMENTED:                                              \
        throw MoFEMExceptionRepeat(MOFEM_NOT_IMPLEMENTED, LINE, FILE, FUNC);   \
      case MOFEM_NOT_FOUND:                                                    \
        throw MoFEMExceptionRepeat(MOFEM_NOT_FOUND, LINE, FILE, FUNC);         \
      case MOFEM_OPERATION_UNSUCCESSFUL:                                       \
        throw MoFEMExceptionRepeat(MOFEM_OPERATION_UNSUCCESSFUL, LINE, FILE,   \
                                   FUNC);                                      \
      case MOFEM_IMPOSIBLE_CASE:                                               \
        throw MoFEMExceptionRepeat(MOFEM_IMPOSIBLE_CASE, LINE, FILE, FUNC);    \
      case MOFEM_INVALID_DATA:                                                 \
        throw MoFEMExceptionRepeat(MOFEM_INVALID_DATA, LINE, FILE, FUNC);      \
      default:                                                                 \
        throw MoFEMExceptionRepeat(MOFEM_MOFEMEXCEPTION_THROW, LINE, FILE,     \
                                   FUNC);                                      \
      }                                                                        \
    }                                                                          \
    return;                                                                    \
  }

/**
 * @brief Definition for error checker class to create checker operator for
 * moab::ErrorCode
 *
 */
#define OP_ERR_MOAB_ERROR_CODE(LINE, FILE, FUNC)                               \
  inline void operator<<(const moab::ErrorCode err) {                          \
    if (PetscLikely(MB_SUCCESS != err)) {                                      \
      std::string str("MOAB error " + boost::lexical_cast<std::string>(err));  \
      PetscError(PETSC_COMM_SELF, LINE, FUNC, FILE, MOFEM_MOAB_ERROR,          \
                 PETSC_ERROR_INITIAL, str.c_str());                            \
      throw MoFEMExceptionNoRepeat(MOFEM_MOAB_ERROR, LINE, FILE, FUNC);        \
    }                                                                          \
    return;                                                                    \
  }

  OP_ERR_MOFEM_ERROR_CODE(LINE,fILE,fUNC);
  OP_ERR_MOAB_ERROR_CODE(LINE,fILE,fUNC);

  inline ErrorCheckerCode(const char *func, const char *file)
      : fUNC(func), fILE(file) {}
  const char *PETSC_CXX_RESTRICT fUNC;
  const char *PETSC_CXX_RESTRICT fILE;
};

typedef int DofIdx;                  ///< Index of DOF
typedef int FEIdx;                   ///< Index of the element
typedef int EntIdx;                  ///< Index of DOF on the entity
typedef int EntPart;                 ///< Partition owning entity
typedef double FieldData;            ///< Field data type
typedef int ApproximationOrder;      ///< Approximation on the entity
typedef int FieldCoefficientsNumber; ///< Number of field coefficients
const EntityHandle no_handle = 0;

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

typedef std::bitset<BITFIELDID_SIZE> BitFieldId;
typedef std::bitset<BITFEID_SIZE> BitFEId;
typedef std::bitset<BITPROBLEMID_SIZE> BitProblemId;
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
typedef std::vector<int, std::allocator<int> > IntAllocator;
typedef std::vector<double, std::allocator<double> > DoubleAllocator;
DEPRECATED typedef IntAllocator IntAllacator; ///< \deprecated Do not use spelling mistake
DEPRECATED typedef DoubleAllocator DoubleAllacator; ///< \deprecated Do not use spelling mistake
typedef ublas::vector<int, IntAllocator> VectorInt;
typedef ublas::vector<double, DoubleAllocator> VectorDouble;
typedef ublas::matrix<double, ublas::row_major, DoubleAllocator> MatrixDouble;

// bounded vector & matrices
typedef ublas::vector<int, ublas::bounded_array<int, 3> > VectorInt3;
typedef ublas::vector<int, ublas::bounded_array<int, 9> > VectorInt9;
typedef ublas::matrix<double, ublas::row_major, ublas::bounded_array<double, 9> >
    MatrixDouble3by3;
typedef ublas::vector<double, ublas::bounded_array<double, 3> > VectorDouble3;
typedef ublas::vector<double, ublas::bounded_array<double, 9> > VectorDouble9;

// shallow adaptor classes
typedef ublas::vector<double, ublas::shallow_array_adaptor<double> >
    VectorAdaptor;
typedef ublas::matrix<double, ublas::row_major,
                      ublas::shallow_array_adaptor<double> >
    MatrixAdaptor;
typedef ublas::vector<int, ublas::shallow_array_adaptor<int> > VectorIntAdaptor;

typedef std::vector<boost::shared_ptr<MatrixDouble> > ShapeFunctionBasesVector;

template <class X> inline std::string toString(X x) {
  std::ostringstream buffer;
  buffer << x;
  return buffer.str();
}
}

#endif //__COMMON_HPP__

/***************************************************************************//**
 * \defgroup mofem MoFEM
 ******************************************************************************/
