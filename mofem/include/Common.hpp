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

struct MoFEMException : public std::exception {
  MoFEMErrorCodes errorCode;
  char errorMessage[255];
  MoFEMException(MoFEMErrorCodes error_code) : errorCode(error_code) {
    strcpy(errorMessage, "Huston we have a problem, something is wrong");
  }
  MoFEMException(MoFEMErrorCodes error_code, const char error_message[])
      : errorCode(error_code) {
    strcpy(errorMessage, error_message);
  }
  const char *what() const throw() { return errorMessage; }
};

typedef moab::ErrorCode MoABErrorCode;

template <typename TYPE> struct MoFEMErrorCodeGeneric {
  MoFEMErrorCodeGeneric(const TYPE) {}
};

template <> struct MoFEMErrorCodeGeneric<PetscErrorCode> {
  PetscErrorCode iERR;
  MoFEMErrorCodeGeneric(const PetscErrorCode ierr) : iERR(ierr) {}
  inline operator PetscErrorCode() const { return iERR; }
  inline const MoFEMErrorTypes errorType() const { return PETSC_ERROR; }
  inline bool getSuccess() const {
    if (PetscUnlikely(iERR))
      return false;
    return true;
  }
};

template <> struct MoFEMErrorCodeGeneric<moab::ErrorCode> {
  moab::ErrorCode rVAL;
  MoFEMErrorCodeGeneric(const moab::ErrorCode rval) : rVAL(rval) {}
  inline operator moab::ErrorCode() const { return rVAL; }
  inline const MoFEMErrorTypes errorType() const { return MOAB_ERROR; }
  inline bool getSuccess() const {
    if (MB_SUCCESS == rVAL)
      return true;
    return false;
  }
};

typedef MoFEMErrorCodeGeneric<PetscErrorCode> MoFEMErrorCode;

static MoFEMErrorCodeGeneric<moab::ErrorCode> rval =
    MoFEMErrorCodeGeneric<moab::ErrorCode>(MB_SUCCESS);
static MoFEMErrorCode ierr = MoFEMErrorCode(0);


struct ErrorChecker {
  ErrorChecker(const int line): lINE(line) {}
  inline void 
  operator<<(const MoFEMErrorCodeGeneric<PetscErrorCode> err) {
    if (PetscUnlikely(!err.getSuccess())) {
      std::string str("MoFEM-PETSc error " +
                      boost::lexical_cast<std::string>(err) + " at line " +
                      boost::lexical_cast<std::string>(lINE));
      throw MoFEMException(MOFEM_MOFEMEXCEPTION_THROW, str.c_str());
    }
    return;
  }
  inline void
  operator<<(const MoFEMErrorCodeGeneric<moab::ErrorCode> err) {
    if (PetscUnlikely(!err.getSuccess())) {
      std::string str("MOFEM-MOAB error " +
                      boost::lexical_cast<std::string>(err) + " at line " +
                      boost::lexical_cast<std::string>(lINE));
      throw MoFEMException(MOFEM_MOFEMEXCEPTION_THROW, str.c_str());
    }
    return;
  }
  inline void operator<<(const PetscErrorCode err) {
    if (PetscUnlikely(err)) {
      std::string str("PETSc error " + boost::lexical_cast<std::string>(err) +
                      " at line " + boost::lexical_cast<std::string>(lINE));
      throw MoFEMException(MOFEM_MOFEMEXCEPTION_THROW, str.c_str());
    }
    return;
  }
  inline void
  operator<<(const moab::ErrorCode err) {
    if (PetscLikely(MB_SUCCESS != err)) {
      std::string str("MOAB error " + boost::lexical_cast<std::string>(err) +
                      " at line " + boost::lexical_cast<std::string>(lINE));
      throw MoFEMException(MOFEM_MOFEMEXCEPTION_THROW, str.c_str());
    }
    return;
  }
private:
  int lINE;
};

struct ErrorCheckerLine {
  inline ErrorChecker operator<<(int line) {
    return ErrorChecker(line);
  }
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
