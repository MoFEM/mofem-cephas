/** \file common.hpp
 * \brief Myltindex containes, data structures and other low-level functions 
 * 
 * Copyright (C) 2013, Lukasz Kaczmarczyk (likask AT wp.pl) <br>
 *
 * The MoFEM package is copyrighted by Lukasz Kaczmarczyk. 
 * It can be freely used for educational and research purposes 
 * by other institutions. If you use this softwre pleas cite my work. 
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

#include "config.h"

//STL
#include<string>
#include<ostream>
#include<sstream>
#include<algorithm>
#include<set>
#include<float.h>
#include<limits.h>
#include<bitset>
#include<exception>

//BOOST
#include <boost/multi_index_container.hpp>
#include <boost/multi_index/ordered_index.hpp>
#include <boost/multi_index/hashed_index.hpp>
#include <boost/multi_index/sequenced_index.hpp>
#include <boost/multi_index/identity.hpp>
#include <boost/multi_index/member.hpp>
#include <boost/multi_index/mem_fun.hpp>
#include <boost/multi_index/global_fun.hpp>
#include <boost/multi_index/composite_key.hpp>
#include <boost/iterator/transform_iterator.hpp>
#include <boost/multiprecision/cpp_int.hpp>
#include <boost/utility/string_ref.hpp>

//MOAB
#include<moab_mpi.h>
#include<moab/ParallelComm.hpp>
#include<MBParallelConventions.h>
#include<moab/Core.hpp>
#include<moab/Interface.hpp>
#include<moab/Skinner.hpp>
#include<moab/GeomUtil.hpp>
#include<moab/Range.hpp>
#include<moab/MeshTopoUtil.hpp>
#include<moab/MergeMesh.hpp>
#include<moab/AdaptiveKDTree.hpp>
#include<MBTagConventions.hpp>
#include<io/Tqdcfr.hpp>

//PETSC
#include<petscmat.h>
#include<petscao.h>
#include<petscbt.h>
#include<petscmat.h>
#include<petscao.h>
#include<petscbt.h>
#include<petsclog.h>
#include<petscsnes.h>
#include<petscts.h>

//BLAS
#ifdef __APPLE__
  #include <Accelerate/Accelerate.h>
  #include<lapack_wrap.h>
#else 
  #include<cblas.h>
  #include<lapack_wrap.h>
#endif

//MOFEM
#include<FEM.h>
#include<H1HdivHcurlL2.h>

//DEFINES
#define MYPCOMM_INDEX 0
//This Is form MOAB
#define MB_TYPE_WIDTH 4
#define MB_ID_WIDTH (8*sizeof(EntityHandle)-MB_TYPE_WIDTH)
#define MB_TYPE_MASK ((EntityHandle)0xF << MB_ID_WIDTH)
//             2^MB_TYPE_WIDTH-1 ------^

#define MB_START_ID ((EntityID)1)        //!< All entity id's currently start at 1
#define MB_END_ID ((EntityID)MB_ID_MASK) //!< Last id is the complement of the MASK
#define MB_ID_MASK (~MB_TYPE_MASK)

#define NOT_USED(x) ( (void)(x) )

/** \brief set barrier start
 *
 * Run code in sequence, starting from process 0, and ends on last process.
 */
#define BARRIER_RANK_START(PCMB) \
  { for(unsigned int i = 0; \
  i<PCMB->proc_config().proc_rank(); i++) MPI_Barrier(PCMB->proc_config().proc_comm()); };
/// set barrier end
#define BARRIER_RANK_END(PCMB) \
  { for(unsigned int i = PCMB->proc_config().proc_rank(); \
  i<PCMB->proc_config().proc_size(); i++) MPI_Barrier(PCMB->proc_config().proc_comm()); };


//ERROR
/// check moab error
#define CHKERR(a) do { \
  ErrorCode val = (a); \
  if (MB_SUCCESS != val) { \
    std::cerr << "Error code  " << val << " at " << __FILE__ << ":" << __LINE__ << std::endl; \
    assert(1); \
  } \
} while (false) 

/// check moab error and communicate it using petsc interface
#define CHKERR_PETSC(a) do { \
  ErrorCode val = (a); \
  if (MB_SUCCESS != val) { \
    std::ostringstream ss; \
    ss << "Error code  " << val << " at " << __FILE__ << ":" << __LINE__ << std::endl; \
    std::string str(ss.str()); \
    SETERRQ(PETSC_COMM_SELF,1,str.c_str()); \
  } \
} while (false)

#define CHKERR_THROW(a) do { \
  ErrorCode val = (a); \
  if (MB_SUCCESS != val) { \
    std::ostringstream ss; \
    ss << "Error code  " << val << " at " << __FILE__ << ":" << __LINE__ << std::endl; \
    std::string str(ss.str()); \
    throw str.c_str(); \
  } \
} while (false)

#define THROW_AT_LINE(a) { \
  std::ostringstream ss; \
  ss << a << " " << " at " << __FILE__ << ":" << __LINE__ << std::endl; \
  std::string str(ss.str()); \
  throw str.c_str(); \
}

//set that with care, it turns off check for ublas
//#define BOOST_UBLAS_NDEBUG

using namespace moab;
using namespace std;
using boost::multi_index_container;
using namespace boost::multi_index;
using namespace boost::multiprecision;

namespace MoFEM {

//CONSTS

const int max_ApproximationOrder = 5;
const EntityHandle no_handle = (EntityHandle)-1;

//TYPEDEFS
typedef PetscInt DofIdx;
typedef int FEIdx;
typedef int EntIdx;
typedef int EntPart;
//typedef uint128_t UId;
typedef checked_uint128_t UId;
typedef bitset<6> BitRefEdges;
typedef bitset<32/*max number of refinments*/> BitRefLevel;
typedef bitset<32/*max number of fields*/> BitFieldId;
typedef PetscScalar FieldData;
typedef bitset<32/*max number of finite elements*/> BitFEId;
typedef bitset<32/*max number of problems*/> BitProblemId;
typedef short ApproximationOrder;
typedef short ApproximationRank;

//ENUMS
/// approximation space 
enum FieldSpace { 
  NoField = 1, 	///< signel scalar or vector of scalars describe state
  H1, 		///< continuous field
  Hdiv,		///< field with continuous normal traction
  Hcurl,	///< field with continuous tangents
  L2,		///< field with C-1 continuity
  LastSpace 	///< FieldSpace in [ 0, LastSpace )
}; 

enum MoFEMTypes {
  MF_ZERO = 0,
  MF_EXCL = 1<<0
};

/// \brief RowColData
enum RowColData {
  Row,Col,Data,LastRowColData
};

enum by_what { 
  by_row = 1<<0, by_col = 1<<1, by_data = 1<<2,
  by_row_data = 1<<0|1<<2, by_col_data = 1<<1|1<<2, by_row_col = 1<<0|1<<1,
  by_all = 1<<0|1<<1|1<<2 
};

//AUX STRUCTURES

/* This small utility that cascades two key extractors will be
 * used throughout the boost example 
 * http://www.boost.org/doc/libs/1_53_0/libs/multi_index/example/complex_structs.cpp
 */
template<class KeyExtractor1,class KeyExtractor2>
struct key_from_key
{
public:
  typedef typename KeyExtractor1::result_type result_type;

  key_from_key(
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
struct ltbit 
{ inline bool operator()(const id_type& valueA,const id_type& valueB) const {
  return valueA.to_ulong()<valueB.to_ulong(); } };

template <typename id_type>
struct eqbit { 
  inline bool operator()(const id_type& valueA,const id_type& valueB) const {
    return valueA.to_ulong() == valueB.to_ulong();
  }
};

template <typename id_type> 
struct hashbit 
{ inline bool operator()(const id_type& value) const {
  return value.to_ulong(); } };

}

//MULTIINDICES
#include "TagMultiIndices.hpp"
#include "FieldMultiIndices.hpp"
#include "EntsMultiIndices.hpp"
#include "DofsMultiIndices.hpp"
#include "FEMMultiIndices.hpp"
#include "ProblemsMultiIndices.hpp"
#include "AdjacencyMultiIndices.hpp"
#include "BCMultiIndices.hpp"

#endif //__COMMON_HPP__
