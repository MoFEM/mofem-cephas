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
typedef PetscInt DofIdx;
typedef int FEIdx;
typedef int EntIdx;
typedef int EntPart;
typedef PetscScalar FieldData;
typedef int ApproximationOrder;
typedef int ApproximationRank;

//consts
const int max_ApproximationOrder = 10;

//typedefs
const EntityHandle no_handle = 0;
typedef uint128_t UId;
//typedef checked_uint128_tUId;
typedef int ShortId;

/** \brief local unique id
  *
  * It is based on local entity handle
  */
struct LocalUId: public UId {
  LocalUId(): UId() {}
  LocalUId(const UId &u): UId(u) {}
  friend bool operator< (const LocalUId& lhs, const LocalUId& rhs);
  friend bool operator> (const LocalUId& lhs, const LocalUId& rhs);
  friend bool operator==(const LocalUId& lhs, const LocalUId& rhs);
  friend bool operator!=(const LocalUId& lhs, const LocalUId& rhs);
};

inline bool operator< (const LocalUId& lhs, const LocalUId& rhs){ return (UId)rhs > (UId)lhs; }
inline bool operator> (const LocalUId& lhs, const LocalUId& rhs){ return rhs < lhs; }
inline bool operator==(const LocalUId& lhs, const LocalUId& rhs) { return (UId)lhs == (UId)rhs; }
inline bool operator!=(const LocalUId& lhs, const LocalUId& rhs) { return !(lhs == rhs); }


/** \brief local unique id
  *
  * It is based on owner entity handle. Each entity is own by some
  * proc/partition, which own set entity handles, which are not unique across
  * mesh on different processors.
  *
  */
struct GlobalUId: public UId {
  GlobalUId(): UId() {}
  GlobalUId(const UId &u): UId(u) {}
  friend bool operator< (const GlobalUId& lhs, const GlobalUId& rhs);
  friend bool operator> (const GlobalUId& lhs, const GlobalUId& rhs);
  friend bool operator==(const GlobalUId& lhs, const GlobalUId& rhs);
  friend bool operator!=(const GlobalUId& lhs, const GlobalUId& rhs);
};

inline bool operator< (const GlobalUId& lhs, const GlobalUId& rhs){ return (UId)rhs > (UId)lhs; }
inline bool operator> (const GlobalUId& lhs, const GlobalUId& rhs){ return rhs < lhs; }
inline bool operator==(const GlobalUId& lhs, const GlobalUId& rhs) { return (UId)lhs == (UId)rhs; }
inline bool operator!=(const GlobalUId& lhs, const GlobalUId& rhs) { return !(lhs == rhs); }

typedef bitset<BITREFEDGES_SIZE> BitRefEdges;
typedef bitset<BITREFLEVEL_SIZE> BitRefLevel;
typedef bitset<BITFIELDID_SIZE> BitFieldId;
typedef bitset<BITFEID_SIZE> BitFEId;
typedef bitset<BITPROBLEMID_SIZE> BitProblemId;
typedef bitset<BITINTERFACEUID_SIZE> BitIntefaceId;
typedef bitset<BITCOORDSYS_SIZE> BitCoordSysId;

//AUX STRUCTURES

/* This small utility that cascades two key extractors will be
 * used throughout the boost example
 * http://www.boost.org/doc/libs/1_53_0/libs/multi_index/example/complex_structs.cpp
 */
template<class KeyExtractor1,class KeyExtractor2>
struct KeyFromKey
{
public:
  typedef typename KeyExtractor1::result_type result_type;

  KeyFromKey(
    const KeyExtractor1& key1_=KeyExtractor1(),
    const KeyExtractor2& key2_=KeyExtractor2()):
    key1(key1_),key2(key2_)
  {}

  template<typename Arg>
  result_type operator()(Arg& arg)const
  {
    return key1(key2(arg));
  }

private:
  KeyExtractor1 key1;
  KeyExtractor2 key2;
};

template <typename id_type>
struct LtBit
{ inline bool operator()(const id_type& valueA,const id_type& valueB) const {
  return valueA.to_ulong()<valueB.to_ulong(); } };

template <typename id_type>
struct EqBit {
  inline bool operator()(const id_type& valueA,const id_type& valueB) const {
    return valueA.to_ulong() == valueB.to_ulong();
  }
};

template <typename id_type>
struct HashBit
{ inline unsigned int operator()(const id_type& value) const {
  return value.to_ulong(); } };

struct MoFEMException: public std::exception {
  MoFEMErrorCode errorCode;
  char errorMessage[255];
  MoFEMException(MoFEMErrorCode error_code):
    errorCode(error_code) {
    strcpy(errorMessage,"Huston we have a problem, something is wrong");
  }
  MoFEMException(MoFEMErrorCode error_code,const char error_message[]):
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
typedef bitset<32> CubitBCType;


}

#endif //__COMMON_HPP__

/***************************************************************************//**
 * \defgroup mofem MoFEM
 ******************************************************************************/
