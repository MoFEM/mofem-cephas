/** \file common.hpp
 * \brief Mylti-index containers, data structures and other low-level functions
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

  typedef ErrorCode MoABErrorCode;

  //typedefs
  typedef int DofIdx; ///< Index of DOF
  typedef int FEIdx;  ///< Index of the element
  typedef int EntIdx;   ///< Index of DOF on the entity
  typedef int EntPart;  ///< Partition owning entity
  typedef double FieldData;  ///< Field data type
  typedef int ApproximationOrder; ///< Approximation on the entity
  typedef int FieldCoefficientsNumber; ///< Number of field coefficients

  //typedefs
  const EntityHandle no_handle = 0;
  //typedef checked_uint128_t UId;
  typedef uint128_t UId;
  typedef UId LocalUId; ///< Local unique id
  typedef UId GlobalUId; ///< Global unique id
  typedef int ShortId;

  #define UID_DOF_MAK 0x1FF

  typedef std::bitset<BITREFEDGES_SIZE> BitRefEdges;
  typedef std::bitset<BITREFLEVEL_SIZE> BitRefLevel;
  typedef std::bitset<BITFIELDID_SIZE> BitFieldId;
  typedef std::bitset<BITFEID_SIZE> BitFEId;
  typedef std::bitset<BITPROBLEMID_SIZE> BitProblemId;
  typedef std::bitset<BITINTERFACEUID_SIZE> BitIntefaceId;

  //AUX STRUCTURES

  /* This small utility that cascades two key extractors will be
  * used throughout the boost example
  * http://www.boost.org/doc/libs/1_53_0/libs/multi_index/example/complex_structs.cpp
  */
  template<class KeyExtractor1,class KeyExtractor2>
  struct KeyFromKey {
  public:
    typedef typename KeyExtractor1::result_type result_type;

    KeyFromKey(
      const KeyExtractor1& key1_=KeyExtractor1(),
      const KeyExtractor2& key2_=KeyExtractor2()
    ):
    key1(key1_),
    key2(key2_) {}

    template<typename Arg>
    result_type operator()(Arg& arg) const {
      return key1(key2(arg));
    }

  private:
    KeyExtractor1 key1;
    KeyExtractor2 key2;
  };

  template <typename id_type>
  struct LtBit {
    inline bool operator()(const id_type& valueA,const id_type& valueB) const {
      return valueA.to_ulong()<valueB.to_ulong();
    }
  };

  template <typename id_type>
  struct EqBit {
    inline bool operator()(const id_type& valueA,const id_type& valueB) const {
      return valueA.to_ulong() == valueB.to_ulong();
    }
  };

  template <typename id_type>
  struct HashBit {
    inline unsigned int operator()(const id_type& value) const {
      return value.to_ulong();
    }
  };

  struct MoFEMException: public std::exception {
    MoFEMErrorCodes errorCode;
    char errorMessage[255];
    MoFEMException(MoFEMErrorCodes error_code):
    errorCode(error_code) {
      strcpy(errorMessage,"Huston we have a problem, something is wrong");
    }
    MoFEMException(MoFEMErrorCodes error_code,const char error_message[]):
    errorCode(error_code) {
      strcpy(errorMessage,error_message);
    }
    const char* what() const throw() {
      return errorMessage;
    }
  };

  /**
  * \typedef CubitBCType
  * bc & material meshsets
  *
  */
  typedef std::bitset<32> CubitBCType;

  // array with std allocators (i.e. concept of capaity is useful here)
  typedef ublas::unbounded_array<int,std::allocator<int> > IntAllacator;
  typedef ublas::unbounded_array<double,std::allocator<double> > DoubleAllacator;
  typedef ublas::unbounded_array<double,std::allocator<double> > DoubleMatrixAllacator;

  // bounded vector
  typedef ublas::vector<int,IntAllacator > VectorInt;
  typedef ublas::vector<double,DoubleAllacator > VectorDouble;
  typedef ublas::matrix<double,ublas::row_major, DoubleMatrixAllacator > MatrixDouble;
  typedef ublas::matrix<double,ublas::row_major,ublas::bounded_array<double,9> > MatrixDouble3by3;
  typedef ublas::vector<double,ublas::bounded_array<double,3> > VectorDouble3;

  // shallow adaptor classes
  typedef ublas::vector<double,ublas::shallow_array_adaptor<double> > VectorAdaptor;
  typedef ublas::matrix<double,ublas::row_major,ublas::shallow_array_adaptor<double> > MatrixAdaptor;
  typedef ublas::vector<int,ublas::shallow_array_adaptor<int> > VectorIntAdaptor;

  typedef std::vector<boost::shared_ptr<MatrixDouble> > ShapeFunctionBasesVector;

  template<class X>
  inline std::string toString(X x) {
    std::ostringstream buffer;
    buffer << x;
    return buffer.str();
  }

  #if PETSC_VERSION_GE(3,7,0)

  /**
  \deprected Funtion is deprected use <http://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/Sys/PetscOptionsGetInt.html>
  */
  DEPRECATED inline PetscErrorCode  PetscOptionsGetInt(const char pre[],const char name[],PetscInt *ivalue,PetscBool  *set) {
    PetscErrorCode ierr;
    PetscFunctionBegin;
    ierr = ::PetscOptionsGetInt(PETSC_NULL,pre,name,ivalue,set); CHKERRQ(ierr);
    PetscFunctionReturn(0);
  }

  /**
  \deprected Funtion is deprected use <http://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/Sys/PetscOptionsGetReal.html>
  */
  DEPRECATED inline PetscErrorCode PetscOptionsGetReal(const char pre[],const char name[],PetscReal *dval,PetscBool *set) {
    PetscErrorCode ierr;
    PetscFunctionBegin;
    ierr  = ::PetscOptionsGetReal(PETSC_NULL,pre,name,dval,set); CHKERRQ(ierr);
    PetscFunctionReturn(0);
  }

  /**
  \deprected Funtion is deprected use <http://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/Sys/PetscOptionsGetScalar.html>
  */
  DEPRECATED inline PetscErrorCode PetscOptionsGetScalar(const char pre[],const char name[],PetscScalar *dval,PetscBool *set) {
    PetscErrorCode ierr;
    PetscFunctionBegin;
    ierr  = ::PetscOptionsGetScalar(PETSC_NULL,pre,name,dval,set); CHKERRQ(ierr);
    PetscFunctionReturn(0);
  }

  /**
  \deprected Funtion is deprected use <http://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/Sys/PetscOptionsGetString.html>
  */
  DEPRECATED inline PetscErrorCode PetscOptionsGetString(const char pre[],const char name[],char str[],size_t size,PetscBool *set) {
    PetscErrorCode ierr;
    PetscFunctionBegin;
    ierr  = ::PetscOptionsGetString(PETSC_NULL,pre,name,str,size,set); CHKERRQ(ierr);
    PetscFunctionReturn(0);
  }

  /**
  \deprected Funtion is deprected use <http://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/Sys/PetscOptionsGetBool.html>
  */
  DEPRECATED inline PetscErrorCode PetscOptionsGetBool(const char pre[],const char name[],PetscBool  *bval,PetscBool *set) {
    PetscErrorCode ierr;
    PetscFunctionBegin;
    ierr  = ::PetscOptionsGetBool(PETSC_NULL,pre,name,bval,set); CHKERRQ(ierr);
    PetscFunctionReturn(0);
  }

  /**
  \deprected Funtion is deprected use <http://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/Sys/PetscOptionsGetRealArray.html>
  */
  DEPRECATED inline PetscErrorCode PetscOptionsGetRealArray(const char pre[],const char name[],PetscReal dval[],PetscInt *nmax,PetscBool *set) {
    PetscErrorCode ierr;
    PetscFunctionBegin;
    ierr  = ::PetscOptionsGetRealArray(PETSC_NULL,pre,name,dval,nmax,set); CHKERRQ(ierr);
    PetscFunctionReturn(0);
  }

  /**
  \deprected Funtion is deprected use <http://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/Sys/PetscOptionsGetEList.html>
  */
  DEPRECATED inline PetscErrorCode PetscOptionsGetEList(
    const char pre[],const char name[],const char*const* list,PetscInt next,PetscInt *value,PetscBool *set
  ) {
    PetscErrorCode ierr;
    PetscFunctionBegin;
    ierr  = ::PetscOptionsGetEList(PETSC_NULL,pre,name,list,next,value,set); CHKERRQ(ierr);
    PetscFunctionReturn(0);
  }

  /**
  \deprected Funtion is deprected use <http://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/Sys/PetscOptionsGetIntArray.html>
  */
  DEPRECATED inline PetscErrorCode PetscOptionsGetIntArray(
    const char pre[],const char name[],PetscInt dvalue[],PetscInt *nmax,PetscBool  *set
  ) {
    PetscErrorCode ierr;
    PetscFunctionBegin;
    ierr = ::PetscOptionsGetIntArray(PETSC_NULL,pre,name,dvalue,nmax,set); CHKERRQ(ierr);
    PetscFunctionReturn(0);
  }

  /**
  \deprected Funtion is deprected use <http://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/Sys/PetscOptionsGetScalarArray.html>
  */
  DEPRECATED inline PetscErrorCode PetscOptionsGetScalarArray(
    const char pre[],const char name[],PetscScalar dvalue[],PetscInt *nmax,PetscBool  *set
  ) {
    PetscErrorCode ierr;
    PetscFunctionBegin;
    ierr = ::PetscOptionsGetScalarArray(PETSC_NULL,pre,name,dvalue,nmax,set); CHKERRQ(ierr);
    PetscFunctionReturn(0);
  }

  #else

  inline PetscErrorCode  PetscOptionsGetInt(PetscOptions *,const char pre[],const char name[],PetscInt *ivalue,PetscBool  *set) {
    PetscErrorCode ierr;
    PetscFunctionBegin;
    ierr = ::PetscOptionsGetInt(pre,name,ivalue,set); CHKERRQ(ierr);
    PetscFunctionReturn(0);
  }

  inline PetscErrorCode PetscOptionsGetReal(PetscOptions *,const char pre[],const char name[],PetscReal *dval,PetscBool *set) {
    PetscErrorCode ierr;
    PetscFunctionBegin;
    ierr  = ::PetscOptionsGetReal(pre,name,dval,set); CHKERRQ(ierr);
    PetscFunctionReturn(0);
  }

  inline PetscErrorCode PetscOptionsGetScalar(PetscOptions *,const char pre[],const char name[],PetscScalar *dval,PetscBool *set) {
    PetscErrorCode ierr;
    PetscFunctionBegin;
    ierr  = ::PetscOptionsGetScalar(pre,name,dval,set); CHKERRQ(ierr);
    PetscFunctionReturn(0);
  }

  inline PetscErrorCode PetscOptionsGetString(PetscOptions *,const char pre[],const char name[],char str[],size_t size,PetscBool *set) {
    PetscErrorCode ierr;
    PetscFunctionBegin;
    ierr  = ::PetscOptionsGetString(pre,name,str,size,set); CHKERRQ(ierr);
    PetscFunctionReturn(0);
  }

  inline PetscErrorCode PetscOptionsGetBool(PetscOptions *,const char pre[],const char name[],PetscBool  *bval,PetscBool *set) {
    PetscErrorCode ierr;
    PetscFunctionBegin;
    ierr  = ::PetscOptionsGetBool(pre,name,bval,set); CHKERRQ(ierr);
    PetscFunctionReturn(0);
  }

  inline PetscErrorCode PetscOptionsGetRealArray(PetscOptions *,const char pre[],const char name[],PetscReal dval[],PetscInt *nmax,PetscBool *set) {
    PetscErrorCode ierr;
    PetscFunctionBegin;
    ierr  = ::PetscOptionsGetRealArray(pre,name,dval,nmax,set); CHKERRQ(ierr);
    PetscFunctionReturn(0);
  }

  inline PetscErrorCode PetscOptionsGetEList(
    PetscOptions *,const char pre[],const char name[],const char*const* list,PetscInt next,PetscInt *value,PetscBool *set
  ) {
    PetscErrorCode ierr;
    PetscFunctionBegin;
    ierr  = ::PetscOptionsGetEList(pre,name,list,next,value,set); CHKERRQ(ierr);
    PetscFunctionReturn(0);
  }

  inline PetscErrorCode PetscOptionsGetIntArray(
    PetscOptions options,const char pre[],const char name[],PetscInt dvalue[],PetscInt *nmax,PetscBool  *set
  ) {
    PetscErrorCode ierr;
    PetscFunctionBegin;
    ierr = ::PetscOptionsGetIntArray(pre,name,dvalue,nmax,set); CHKERRQ(ierr);
    PetscFunctionReturn(0);
  }

  inline PetscErrorCode PetscOptionsGetScalarArray(
    PetscOptions options,const char pre[],const char name[],PetscScalar dvalue[],PetscInt *nmax,PetscBool  *set
  ) {
    PetscErrorCode ierr;
    PetscFunctionBegin;
    ierr = ::PetscOptionsGetScalarArray(pre,name,dvalue,nmax,set); CHKERRQ(ierr);
    PetscFunctionReturn(0);

  }

  #endif

}

#endif //__COMMON_HPP__

/***************************************************************************//**
 * \defgroup mofem MoFEM
 ******************************************************************************/
