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

#include<string>
#include<ostream>
#include<sstream>
#include<algorithm>
#include<set>
#include<float.h>
#include<limits.h>
#include<bitset>
#include<exception>

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

#define NOT_USED(x) ( (void)(x) )

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

//This Is form MOAB
#define MB_TYPE_WIDTH 4
#define MB_ID_WIDTH (8*sizeof(EntityHandle)-MB_TYPE_WIDTH)
#define MB_TYPE_MASK ((EntityHandle)0xF << MB_ID_WIDTH)
//             2^MB_TYPE_WIDTH-1 ------^

#define MB_START_ID ((EntityID)1)        //!< All entity id's currently start at 1
#define MB_END_ID ((EntityID)MB_ID_MASK) //!< Last id is the complement of the MASK
#define MB_ID_MASK (~MB_TYPE_MASK)


#include<petscmat.h>
#include<petscao.h>
#include<petscbt.h>

#include<petscmat.h>
#include<petscao.h>
#include<petscbt.h>
#include<petsclog.h>
#include<petscsnes.h>

#include<petscts.h>

#ifdef __APPLE__
  #include <Accelerate/Accelerate.h>
  #include<lapack_wrap.h>
#else 
  #include<cblas.h>
  #include<lapack_wrap.h>
#endif

#include<FEM.h>
#include<H1HdivHcurlL2.h>

#define MYPCOMM_INDEX 0

using namespace moab;
using namespace std;
using boost::multi_index_container;
using namespace boost::multi_index;
using namespace boost::multiprecision;


//set that with care, it turns off check for ublas
//#define BOOST_UBLAS_NDEBUG

namespace MoFEM {

const int max_ApproximationOrder = 5;

typedef PetscInt DofIdx;
typedef int FEIdx;
typedef int EntIdx;
typedef int EntPart;
//typedef uint128_t UId;
typedef checked_uint128_t UId;

typedef bitset<6> BitRefEdges;
typedef bitset<16/*max number of refinments*/> BitRefLevel;
typedef bitset<16/*max number of fields*/> BitFieldId;

/// approximation space 
enum FieldSpace { 
  NoField = 1, 	///< signel scalar or vector of scalars describe state
  H1, 		///< continuous field
  Hdiv,		///< field with continuous normal traction
  Hcurl,	///< field with continuous tangents
  L2,		///< field with C-1 continuity
  LastSpace 	///< FieldSpace in [ 0, LastSpace )
}; 
typedef PetscScalar FieldData;

//
typedef bitset<16/*max number of finite elements*/> BitFEId;
typedef bitset<16/*max number of problems*/> BitProblemId;
//
typedef short ApproximationOrder;
typedef short ApproximationRank;

const EntityHandle no_handle = (EntityHandle)-1;

/// \brief RowColData
enum RowColData {
  Row,Col,Data,LastRowColData
};

enum by_what { 
  by_row = 1<<0, by_col = 1<<1, by_data = 1<<2,
  by_row_data = 1<<0|1<<2, by_col_data = 1<<1|1<<2, by_row_col = 1<<0|1<<1,
  by_all = 1<<0|1<<1|1<<2 
};

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

// tags and unaryFunction functions

/// MultiIndex Tag for field id 
struct CubitMeshSets_mi_tag {};
struct CubitMeshSets_mask_meshset_mi_tag {};
struct CubitMeshSets_bc_data_mi_tag {};
struct BitFieldId_mi_tag {};
struct Unique_mi_tag {};
struct MoABEnt_mi_tag {};
struct EntType_mi_tag {};
struct Composite_unique_mi_tag {};
struct Composite_mi_tag {};
struct Composite_mi_tag2 {};
struct Composite_mi_tag3 {};
struct MoFEMFiniteElement_Meshset_mi_tag {};
struct BitFEId_mi_tag {};
struct MoFEMFiniteElement_name_mi_tag {};
struct SideNumber_mi_tag{};
struct MoABEnt_MoABEnt_mi_tag {};
struct MoABEnt_mi_tag2 {};
/// MultiIndex Tag for dof MoFEM index
struct Idx_mi_tag { 
  static const bool IamNotPartitioned;
  /// extract dof index from iterator 
  template<class IT>
  static DofIdx get_index(const IT &it) { return it->dof_idx; }
};
struct PetscGlobalIdx_mi_tag {};
struct PetscLocalIdx_mi_tag {};
/// MultiIndex Tag for partition number
struct Part_mi_tag {
  static const bool IamNotPartitioned;
  /// extract global dof index from iterator 
  template<class IT>
  static DofIdx get_index(const IT &it) { return it->petsc_gloabl_dof_idx; }
};
struct Unique_MoABEnt_mi_tag {};
struct Unique_MoFEMFiniteElement_mi_tag {};
struct MoABEnt_MoFEMFiniteElement_mi_tag {};
struct Meshset_mi_tag {};
/// MultiIndex Tag for field name
struct FieldName_mi_tag {};
struct BitFieldId_space_mi_tag {};
struct MoFEMFiniteElement_Part_mi_tag {};
struct BitProblemId_mi_tag {};
struct MoFEMProblem_mi_tag {};

/**
 * \brief keeps information about side number for the finite element
 */
struct SideNumber {
  EntityHandle ent;
  int side_number;
  int sense;
  int offset;
  int brother_side_number;
  inline EntityType get_ent_type() const { return (EntityType)((ent&MB_TYPE_MASK)>>MB_ID_WIDTH); }
  SideNumber(EntityHandle _ent,int _side_number,int _sense,int _offset):
    ent(_ent),side_number(_side_number),sense(_sense),offset(_offset),brother_side_number(-1) {};
};


/** 
 * @relates multi_index_container
 * \brief SideNumber_multiIndex for SideNumber
 *
 *  \param hashed_unique<
 *     member<SideNumber,EntityHandle,&SideNumber::ent> >,
 *  \param ordered_non_unique<
 *     composite_key<
 *	SideNumber,
 *	const_mem_fun<SideNumber,EntityType,&SideNumber::get_ent_type>,
 *	member<SideNumber,int,&SideNumber::side_number> > >,
 *  \param ordered_non_unique<
 *     const_mem_fun<SideNumber,EntityType,&SideNumber::get_ent_type> >
 */
typedef multi_index_container<
  SideNumber,
  indexed_by<
    hashed_unique<
      member<SideNumber,EntityHandle,&SideNumber::ent> >,
    ordered_non_unique<
      composite_key<
	SideNumber,
	const_mem_fun<SideNumber,EntityType,&SideNumber::get_ent_type>,
	member<SideNumber,int,&SideNumber::side_number> > >,
    ordered_non_unique<
      const_mem_fun<SideNumber,EntityType,&SideNumber::get_ent_type> >
  > > SideNumber_multiIndex; 

/** 
 * \brief this struct keeps basic methods for moab entity
 */
struct BasicMoFEMEntity {
  EntityHandle ent;
  /// \param ent handle to moab entity
  BasicMoFEMEntity(const EntityHandle _ent);
  /// get entity type
  inline EntityType get_ent_type() const { return (EntityType)((ent&MB_TYPE_MASK)>>MB_ID_WIDTH); }
  /// get entity id
  inline EntityID get_ent_id() const { return (EntityID)(ent&MB_ID_MASK); };
};

/** 
 * \brief struct keeps data about selected prism adjacencies, and potentially other entities
 */
struct AdjacencyMapForBasicMoFEMEntity: public BasicMoFEMEntity {
  BasicMoFEMEntity Adj; ///< adjacent entity to this AdjacencyMapForBasicMoFEMEntity
  AdjacencyMapForBasicMoFEMEntity(const EntityHandle _ent,const EntityHandle adj):
    BasicMoFEMEntity(_ent), Adj(adj) {};
  inline EntityHandle get_adj() const { return Adj.ent; };
  inline EntityType get_adj_type() const { return Adj.get_ent_type(); };
};

/** 
 * \brief struct keeps handle to refined handle.
 */
struct RefMoFEMEntity: public BasicMoFEMEntity {
  const EntityHandle *tag_parent_ent;
  BitRefLevel *tag_BitRefLevel;
  RefMoFEMEntity(Interface &moab,const EntityHandle _ent);
  /// get entity
  inline EntityHandle get_ref_ent() const { return ent; }
  /// get patent entity
  inline EntityType get_parent_ent_type() const { return (EntityType)((*tag_parent_ent&MB_TYPE_MASK)>>MB_ID_WIDTH); }
  /// get entity ref bit refinment signature
  inline const BitRefLevel& get_BitRefLevel() const { return *tag_BitRefLevel; }
  /// get parent entity, i.e. entity form one refinment level up
  inline EntityHandle get_parent_ent() const { return *tag_parent_ent; }
  const RefMoFEMEntity* get_RefMoFEMEntity_ptr() { return this; }
  friend ostream& operator<<(ostream& os,const RefMoFEMEntity& e);
};

/** 
 * \brief interface to RefMoFEMEntity
 */
template <typename T>
struct interface_RefMoFEMEntity {
  const T *ref_ptr;
  interface_RefMoFEMEntity(const T *_ref_ptr): ref_ptr(_ref_ptr) {}
  inline EntityHandle get_ref_ent() const { return ref_ptr->get_ref_ent(); }
  inline EntityHandle get_parent_ent() const { return ref_ptr->get_parent_ent(); }
  inline const BitRefLevel& get_BitRefLevel() const { return ref_ptr->get_BitRefLevel(); }
  inline EntityType get_ent_type() const { return ref_ptr->get_ent_type(); };
  inline EntityType get_parent_ent_type() const { return ref_ptr->get_parent_ent_type(); };
  inline EntityID get_ent_id() const { return ref_ptr->get_ent_id(); };
  inline const RefMoFEMEntity* get_RefMoFEMEntity_ptr() { return ref_ptr->get_RefMoFEMEntity_ptr(); }
};

/**
 * \brief keeps data about abstract refined finite element
 */
struct RefMoFEMElement: public interface_RefMoFEMEntity<RefMoFEMEntity> {
  typedef interface_RefMoFEMEntity<RefMoFEMEntity> interface_type_RefMoFEMEntity;
  BitRefEdges *tag_BitRefEdges;
  SideNumber_multiIndex side_number_table;
  RefMoFEMElement(Interface &moab,const RefMoFEMEntity *_RefMoFEMEntity_ptr);
  inline const BitRefEdges& get_BitRefEdges() const { return *tag_BitRefEdges; }
  int get_BitRefEdges_ulong() const { return get_BitRefEdges().to_ulong(); }
  SideNumber_multiIndex &get_side_number_table() const { return const_cast<SideNumber_multiIndex&>(side_number_table); };
  virtual SideNumber* get_side_number_ptr(Interface &moab,EntityHandle ent) const {
    NOT_USED(moab);
    NOT_USED(ent);
    return NULL; 
  };
  const RefMoFEMElement* get_RefMoFEMElement() const { return this; }
  friend ostream& operator<<(ostream& os,const RefMoFEMElement& e);
};

/**
 * \brief keeps data about abstract MESHSET finite element
 */
struct RefMoFEMElement_MESHSET: public RefMoFEMElement {
  RefMoFEMElement_MESHSET(Interface &moab,const RefMoFEMEntity *_RefMoFEMEntity_ptr);
  const RefMoFEMElement* get_RefMoFEMElement() const { return this; }
  SideNumber* get_side_number_ptr(Interface &moab,EntityHandle ent) const;
};

/**
 * \brief keeps data about abstract PRISM finite element
 */
struct RefMoFEMElement_PRISM: public RefMoFEMElement {
  RefMoFEMElement_PRISM(Interface &moab,const RefMoFEMEntity *_RefMoFEMEntity_ptr);
  const RefMoFEMElement* get_RefMoFEMElement() const { return this; }
  SideNumber* get_side_number_ptr(Interface &moab,EntityHandle ent) const;
};

/**
 * \brief keeps data about abstract TET finite element
 */
struct RefMoFEMElement_TET: public RefMoFEMElement {
  const int* tag_type_data;
  RefMoFEMElement_TET(Interface &moab,const RefMoFEMEntity *_RefMoFEMEntity_ptr);
  const RefMoFEMElement* get_RefMoFEMElement() const { return this; }
  SideNumber* get_side_number_ptr(Interface &moab,EntityHandle ent) const;
  inline int get_ref_type() const { return tag_type_data[0]; }
  inline int get_ref_sub_type() const { return tag_type_data[1]; }
  friend ostream& operator<<(ostream& os,const RefMoFEMElement_TET& e);
};

/**
 * \brief intrface to RefMoFEMElement
 */
template<typename T>
struct interface_RefMoFEMElement: interface_RefMoFEMEntity<T> {
  interface_RefMoFEMElement(const T *_ref_ptr): interface_RefMoFEMEntity<T>(_ref_ptr) {}
  int get_BitRefEdges_ulong() const { return interface_RefMoFEMEntity<T>::ref_ptr->get_BitRefEdges_ulong(); }
  SideNumber_multiIndex &get_side_number_table() const { return interface_RefMoFEMEntity<T>::ref_ptr->get_side_number_table(); }
  SideNumber* get_side_number_ptr(Interface &moab,EntityHandle ent) const { return interface_RefMoFEMEntity<T>::ref_ptr->get_side_number_ptr(moab,ent); }
  inline const RefMoFEMElement* get_RefMoFEMElement() const { return interface_RefMoFEMEntity<T>::ref_ptr->get_RefMoFEMElement(); }
};

/// \brief keeps data about field
struct MoFEMField {
  EntityHandle meshset; 		///< keeps entities for this meshset
  Tag th_FieldData,th_AppOrder;
  Tag th_AppDofOrder,th_DofRank;
  BitFieldId* tag_id_data; 		///< tag keeps field id
  FieldSpace* tag_space_data;		///< tag keeps field space
  ApproximationRank* tag_rank_data; 	///< tag keeps field rank (dimension, f.e. temerature field has rank 1, displacemenst field in 3d has rank 3)
  const void* tag_name_data; 		///< tag keeps name of the field
  int tag_name_size; 			///< number of bits necessery to keep field name
  const void* tag_name_prefix_data; 	///< tag keeps name prefix of the field
  int tag_name_prefix_size; 		///< number of bits necessery to keep field name prefix
  int (*forder_entityset)(int); 	///< nb. dofs on meshset for given space
  int (*forder_vertex)(int); 		///< nb. dofs on node for given space
  int (*forder_edge)(int); 		///< nb. dofs on edge for given space
  int (*forder_face)(int); 		///< nb. dofs on face for given space
  int (*forder_elem)(int); 		///< nb. dofs on elem for given space
  /**
    * \brief constructor for moab field
    *
    * \param _meshset meshset which keeps entities for this field
    */
  MoFEMField(Interface &moab,const EntityHandle _meshset);					
  inline EntityHandle get_meshset() const { return meshset; };
  inline const BitFieldId& get_id() const { return *((BitFieldId*)tag_id_data); }; 			
  inline boost::string_ref get_name_ref() const { return boost::string_ref((char *)tag_name_data,tag_name_size); };	
  inline string get_name() const { return string((char *)tag_name_data,tag_name_size); };	
  inline FieldSpace get_space() const { return *tag_space_data; };
  inline ApproximationRank get_max_rank() const { return *tag_rank_data; };
  inline unsigned int get_bit_number() const { return ffsl(((BitFieldId*)tag_id_data)->to_ulong()); }
  const MoFEMField* get_MoFEMField_ptr() const { return this; };
  friend ostream& operator<<(ostream& os,const MoFEMField& e);
};

/**
 * \brief interface for MoFEMField
 */
template <typename T> 
struct interface_MoFEMField {
  const T *field_ptr;
  interface_MoFEMField(const T *_field_ptr): field_ptr(_field_ptr) {};
  inline EntityHandle get_meshset() const { return field_ptr->get_meshset(); };
  inline const BitFieldId& get_id() const { return field_ptr->get_id(); };
  inline unsigned int get_bit_number() const { return field_ptr->get_bit_number(); }
  inline boost::string_ref get_name_ref() const { return field_ptr->get_name_ref(); };
  inline string get_name() const { return field_ptr->get_name(); };
  inline FieldSpace get_space() const { return field_ptr->get_space(); };
  inline ApproximationRank get_max_rank() const { return field_ptr->get_max_rank(); };
  inline int forder_entityset(int p) const { return field_ptr->forder_entityset(p); };
  inline int forder_vertex(int p) const { return field_ptr->forder_vertex(p); };
  inline int forder_edge(int p) const { return field_ptr->forder_edge(p); };
  inline int forder_face(int p) const { return field_ptr->forder_face(p); };
  inline int forder_elem(int p) const { return field_ptr->forder_elem(p); };
  inline const MoFEMField* get_MoFEMField_ptr() const { return field_ptr->get_MoFEMField_ptr(); };
};

/**
 * \brief struct keeps handle to entity in the field.
 */
struct MoFEMEntity: public interface_MoFEMField<MoFEMField>, interface_RefMoFEMEntity<RefMoFEMEntity> {
  typedef interface_MoFEMField<MoFEMField> interface_type_MoFEMField;
  const RefMoFEMEntity *ref_mab_ent_ptr;
  const ApproximationOrder* tag_order_data;
  const FieldData* tag_FieldData;
  int tag_FieldData_size;
  const ApproximationOrder* tag_dof_order_data;
  const ApproximationRank* tag_dof_rank_data;
  int (*forder)(int);
  UId uid;
  MoFEMEntity(Interface &moab,const MoFEMField *_FieldData,const RefMoFEMEntity *_ref_mab_ent_ptr);
  inline EntityHandle get_ent() const { return get_ref_ent(); }
  inline FieldData* get_FieldData() const { return const_cast<FieldData*>(tag_FieldData); }
  inline int get_order_nb_dofs(int order) const { return forder(order); }
  inline int get_order_nb_dofs_diff(int order) const { return forder(order)-forder(order-1); }
  inline ApproximationOrder get_max_order() const { return *((ApproximationOrder*)tag_order_data); }
  inline const RefMoFEMEntity* get_RefMoFEMEntity_ptr() const { return ref_mab_ent_ptr; }
  const UId& get_unique_id() const { return uid; }
  UId get_unique_id_calculate() const {
    char bit_number = get_bit_number();
    assert(bit_number<=32);
    UId _uid_ = (ref_ptr->ent)|(((UId)bit_number)<<(8*sizeof(EntityHandle)));
    return _uid_;
  }
  const MoFEMEntity* get_MoFEMEntity_ptr() const { return this; };
  friend ostream& operator<<(ostream& os,const MoFEMEntity& e);
};

/**
 * \brief interface to MoFEMEntity
 *
 * interface to MoFEMEntity
 */
template <typename T>
struct interface_MoFEMEntity: public interface_MoFEMField<T>,interface_RefMoFEMEntity<RefMoFEMEntity> {
  interface_MoFEMEntity(const T *_ptr): interface_MoFEMField<T>(_ptr),interface_RefMoFEMEntity<RefMoFEMEntity>(_ptr->get_RefMoFEMEntity_ptr()) {};
  inline EntityHandle get_ent() const { return interface_MoFEMField<T>::get_ent(); }
  inline FieldData* get_FieldData() const { return interface_MoFEMField<T>::field_ptr->get_FieldData(); }
  inline int get_order_nb_dofs(int order) const { return interface_MoFEMField<T>::field_ptr->get_order_nb_dofs(order); }
  inline int get_order_nb_dofs_diff(int order) const { return interface_MoFEMField<T>::field_ptr->get_order_nb_dofs_diff(order); }
  inline ApproximationOrder get_max_order() const { return interface_MoFEMField<T>::field_ptr->get_max_order(); }
  inline const UId& get_unique_id() const { return interface_MoFEMField<T>::field_ptr->get_unique_id(); }
  inline const MoFEMEntity* get_MoFEMEntity_ptr() const { return interface_MoFEMField<T>::field_ptr->get_MoFEMEntity_ptr(); };
  inline const RefMoFEMEntity* get_RefMoFEMEntity_ptr() const { return interface_MoFEMField<T>::field_ptr->get_RefMoFEMEntity_ptr(); }
};

/**
 * \brief structure to chane MoFEMEntity order
 */
struct MoFEMEntity_change_order {
  Interface& moab;
  ApproximationOrder order;
  vector<FieldData> data;
  vector<ApproximationOrder> data_dof_order;
  vector<ApproximationRank> data_dof_rank;
  MoFEMEntity_change_order(Interface& _moab,ApproximationOrder _order): moab(_moab),order(_order) {};
  void operator()(MoFEMEntity &e);
};

/**
 * \brief keeps information about indexed dofs
 */
struct DofMoFEMEntity: public interface_MoFEMEntity<MoFEMEntity> {
  typedef interface_MoFEMField<MoFEMEntity> interface_type_MoFEMField;
  typedef interface_MoFEMEntity<MoFEMEntity> interface_type_MoFEMEntity;
  typedef interface_RefMoFEMEntity<RefMoFEMEntity> interface_type_RefMoFEMEntity;
  static UId get_unique_id_calculate(const DofIdx _dof_,const MoFEMEntity *_ent_ptr_) {
    if(_dof_>=256) THROW_AT_LINE("_dof>=256");
    UId _uid_ = ((UId)_dof_)|((_ent_ptr_->get_unique_id())<<8);
    return _uid_;
  }
  //
  DofIdx dof;
  bool active;
  UId uid;
  DofMoFEMEntity(const MoFEMEntity *_MoFEMEntity_ptr,const ApproximationOrder _dof_order,const ApproximationRank _dof_rank,const DofIdx _dof);
  inline DofIdx get_EntDofIdx() const { return dof; }
  inline FieldData& get_FieldData() const { return const_cast<FieldData&>(field_ptr->tag_FieldData[dof]); }
  UId get_unique_id() const { return uid; };
  UId get_unique_id_calculate() const { return get_unique_id_calculate(dof,get_MoFEMEntity_ptr()); }
  inline EntityHandle get_ent() const { return field_ptr->get_ent(); };
  inline ApproximationOrder get_dof_order() const { return ((ApproximationOrder*)field_ptr->tag_dof_order_data)[dof]; };
  inline ApproximationRank get_dof_rank() const { return ((ApproximationRank*)field_ptr->tag_dof_rank_data)[dof]; };
  inline int get_active() const { return active ? 1 : 0; }
  inline const DofMoFEMEntity* get_DofMoFEMEntity_ptr() const { return this; };
  friend ostream& operator<<(ostream& os,const DofMoFEMEntity& e);
};

/**
 * \brief interface to DofMoFEMEntitys
 */
template <typename T>
struct interface_DofMoFEMEntity: public interface_MoFEMEntity<T> {
  interface_DofMoFEMEntity(const T *_ptr): interface_MoFEMEntity<T>(_ptr) {};
  UId get_unique_id() const { return interface_MoFEMEntity<T>::field_ptr->get_unique_id(); }
  inline DofIdx get_EntDofIdx() const { return interface_MoFEMEntity<T>::field_ptr->get_EntDofIdx(); }
  inline FieldData& get_FieldData() const { return interface_MoFEMEntity<T>::field_ptr->get_FieldData(); }
  inline EntityHandle get_ent() const { return interface_MoFEMEntity<T>::field_ptr->get_ent(); };
  inline ApproximationOrder get_dof_order() const { return interface_MoFEMEntity<T>::field_ptr->get_dof_order(); };
  inline ApproximationRank get_dof_rank() const { return interface_MoFEMEntity<T>::field_ptr->get_dof_rank(); };
  inline int get_active() const { return interface_MoFEMEntity<T>::field_ptr->get_active(); }
  inline const DofMoFEMEntity* get_DofMoFEMEntity_ptr() const { return interface_MoFEMEntity<T>::field_ptr->get_DofMoFEMEntity_ptr(); };
};

/**
 * \brief keeps information about indexed dofs for the problem
 */
struct NumeredDofMoFEMEntity: public interface_DofMoFEMEntity<DofMoFEMEntity> {
  typedef interface_MoFEMField<DofMoFEMEntity> interface_type_MoFEMField;
  typedef interface_MoFEMEntity<DofMoFEMEntity> interface_type_MoFEMEntity;
  typedef interface_DofMoFEMEntity<DofMoFEMEntity> interface_type_DofMoFEMEntity;
  DofIdx dof_idx;
  DofIdx petsc_gloabl_dof_idx;
  DofIdx petsc_local_dof_idx;
  unsigned int part;
  inline DofIdx get_petsc_gloabl_dof_idx() const { return petsc_gloabl_dof_idx;  }
  inline DofIdx get_petsc_local_dof_idx() const { return petsc_local_dof_idx; }
  inline unsigned int get_part() const { return part;  }
  NumeredDofMoFEMEntity(const DofIdx idx,const DofMoFEMEntity* _DofMoFEMEntity_ptr);
  inline const NumeredDofMoFEMEntity* get_NumeredDofMoFEMEntity_ptr() const { return this; };
  inline bool operator<(const NumeredDofMoFEMEntity& _dof) const { return get_unique_id()<_dof.get_unique_id(); }
  friend ostream& operator<<(ostream& os,const NumeredDofMoFEMEntity& e);
};

/**
 * \brief interface to NumeredDofMoFEMEntity
 */
template <typename T>
struct interface_NumeredDofMoFEMEntity: public interface_DofMoFEMEntity<T> {
  interface_NumeredDofMoFEMEntity(const T *_ptr): interface_DofMoFEMEntity<T>(_ptr) {};
  inline DofIdx get_dof_idx() const { return interface_DofMoFEMEntity<T>::field_ptr->dof_idx; }
  inline DofIdx get_petsc_gloabl_dof_idx() const { return interface_DofMoFEMEntity<T>::field_ptr->get_petsc_gloabl_dof_idx();  }
  inline DofIdx get_petsc_local_dof_idx() const { return interface_DofMoFEMEntity<T>::field_ptr->get_petsc_local_dof_idx(); }
  inline unsigned int get_part() const { return interface_DofMoFEMEntity<T>::field_ptr->get_part();  }
  inline const NumeredDofMoFEMEntity* get_NumeredDofMoFEMEntity_ptr() const { return interface_DofMoFEMEntity<T>::field_ptr->get_NumeredDofMoFEMEntity_ptr(); };
};

/**
 * \brief keeps basic information about indexed dofs for the finite element
 */
struct BaseFEDofMoFEMEntity {
  BaseFEDofMoFEMEntity(SideNumber *_side_number_ptr): side_number_ptr(_side_number_ptr) {};
  SideNumber *side_number_ptr;
};

/**
 * \brief keeps information about indexed dofs for the finite element
 */
struct FEDofMoFEMEntity: public BaseFEDofMoFEMEntity,interface_DofMoFEMEntity<DofMoFEMEntity> {
  typedef interface_MoFEMField<DofMoFEMEntity> interface_type_MoFEMField;
  typedef interface_DofMoFEMEntity<DofMoFEMEntity> interface_type_DofMoFEMEntity;
  typedef interface_RefMoFEMEntity<RefMoFEMEntity> interface_type_RefMoFEMEntity;
  FEDofMoFEMEntity(
    SideNumber *_side_number_ptr,
    const DofMoFEMEntity *_DofMoFEMEntity_ptr);
  friend ostream& operator<<(ostream& os,const FEDofMoFEMEntity& e);
};

/**
 * \brief keeps information about indexed dofs for the finite element
 */
struct FENumeredDofMoFEMEntity: public BaseFEDofMoFEMEntity,interface_NumeredDofMoFEMEntity<NumeredDofMoFEMEntity> {
  typedef interface_MoFEMField<NumeredDofMoFEMEntity> interface_type_MoFEMField;
  typedef interface_DofMoFEMEntity<NumeredDofMoFEMEntity> interface_type_DofMoFEMEntity;
  typedef interface_RefMoFEMEntity<RefMoFEMEntity> interface_type_RefMoFEMEntity;
  typedef interface_NumeredDofMoFEMEntity<NumeredDofMoFEMEntity> interface_type_NumeredDofMoFEMEntity;
  FENumeredDofMoFEMEntity(
    SideNumber *_side_number_ptr,
    const NumeredDofMoFEMEntity *_NumeredDofMoFEMEntity_ptr);
  friend ostream& operator<<(ostream& os,const FENumeredDofMoFEMEntity& e);
};

// multi_index_containers

/**
 * @relates multi_index_container
 * \brief MoFEMField_multiIndex for MoFEMField
 *
 * \param hashed_unique<
 *     tag<BitFieldId_mi_tag>, const_mem_fun<MoFEMField,const BitFieldId&,&MoFEMField::get_id>, hashbit<BitFieldId>, eqbit<BitFieldId> >,
 * \param   ordered_unique<
 *     tag<Meshset_mi_tag>, member<MoFEMField,EntityHandle,&MoFEMField::meshset> >,
 * \param hashed_unique<
 *     tag<FieldName_mi_tag>, const_mem_fun<MoFEMField,boost::string_ref,&MoFEMField::get_name_ref> >,
 * \param ordered_non_unique<
 *     tag<BitFieldId_space_mi_tag>, const_mem_fun<MoFEMField,FieldSpace,&MoFEMField::get_space> >
 */
typedef multi_index_container<
  MoFEMField,
  indexed_by<
    hashed_unique<
      tag<BitFieldId_mi_tag>, const_mem_fun<MoFEMField,const BitFieldId&,&MoFEMField::get_id>, hashbit<BitFieldId>, eqbit<BitFieldId> >,
    ordered_unique<
      tag<Meshset_mi_tag>, member<MoFEMField,EntityHandle,&MoFEMField::meshset> >,
    ordered_unique<
      tag<FieldName_mi_tag>, const_mem_fun<MoFEMField,boost::string_ref,&MoFEMField::get_name_ref> >,
    ordered_non_unique<
      tag<BitFieldId_space_mi_tag>, const_mem_fun<MoFEMField,FieldSpace,&MoFEMField::get_space> >
  > > MoFEMField_multiIndex;

typedef multi_index_container<
  const MoFEMField*,
  indexed_by<
    ordered_unique<
      tag<BitFieldId_mi_tag>, const_mem_fun<MoFEMField,const BitFieldId&,&MoFEMField::get_id>, ltbit<BitFieldId> >
   > > MoFEMField_multiIndex_view;

/** 
 * @relates multi_index_container
 * \brief MultiIndex container keeps MoFEMEntity
 *
 * \param ordered_unique<
 *    tag<Unique_mi_tag>, member<MoFEMEntity,UId,&MoFEMEntity::uid> >,
 * \param ordered_non_unique<
 *    tag<BitFieldId_mi_tag>, const_mem_fun<MoFEMEntity::interface_type_MoFEMField,const BitFieldId&,&MoFEMEntity::get_id>, ltbit<BitFieldId> >,
 * \param ordered_non_unique<
 *    tag<FieldName_mi_tag>, const_mem_fun<MoFEMEntity::interface_type_MoFEMField,boost::string_ref,&MoFEMEntity::get_name_ref> >,
 * \param hashed_non_unique<
 *    tag<MoABEnt_mi_tag>, const_mem_fun<MoFEMEntity,EntityHandle,&MoFEMEntity::get_ent> >
 */
typedef multi_index_container<
  MoFEMEntity,
  indexed_by<
    ordered_unique<
      tag<Unique_mi_tag>, member<MoFEMEntity,UId,&MoFEMEntity::uid> >,
    ordered_non_unique<
      tag<BitFieldId_mi_tag>, const_mem_fun<MoFEMEntity::interface_type_MoFEMField,const BitFieldId&,&MoFEMEntity::get_id>, ltbit<BitFieldId> >,
    ordered_non_unique<
      tag<FieldName_mi_tag>, const_mem_fun<MoFEMEntity::interface_type_MoFEMField,boost::string_ref,&MoFEMEntity::get_name_ref> >,
    hashed_non_unique<
      tag<MoABEnt_mi_tag>, const_mem_fun<MoFEMEntity,EntityHandle,&MoFEMEntity::get_ent> >
  > > MoFEMEntity_multiIndex;

typedef multi_index_container<
  DofMoFEMEntity,
  indexed_by<
    ordered_unique< 
      tag<Unique_mi_tag>, member<DofMoFEMEntity,UId,&DofMoFEMEntity::uid> >,
    ordered_non_unique<
      tag<FieldName_mi_tag>, const_mem_fun<DofMoFEMEntity::interface_type_MoFEMField,boost::string_ref,&DofMoFEMEntity::get_name_ref> >,
    ordered_non_unique<
      tag<MoABEnt_mi_tag>, const_mem_fun<DofMoFEMEntity,EntityHandle,&DofMoFEMEntity::get_ent> >,
    ordered_non_unique<
      tag<BitFieldId_mi_tag>, const_mem_fun<DofMoFEMEntity::interface_type_MoFEMField,const BitFieldId&,&DofMoFEMEntity::get_id>, ltbit<BitFieldId> >,
    ordered_non_unique<
      tag<Composite_mi_tag>, 
      composite_key<
	DofMoFEMEntity,
	const_mem_fun<DofMoFEMEntity::interface_type_MoFEMField,boost::string_ref,&DofMoFEMEntity::get_name_ref>,
	const_mem_fun<DofMoFEMEntity,EntityHandle,&DofMoFEMEntity::get_ent>,
	const_mem_fun<DofMoFEMEntity,DofIdx,&DofMoFEMEntity::get_EntDofIdx> 
      > >,
    ordered_non_unique<
      tag<Composite_mi_tag2>, 
      composite_key<
	DofMoFEMEntity,
	const_mem_fun<DofMoFEMEntity::interface_type_MoFEMField,boost::string_ref,&DofMoFEMEntity::get_name_ref>,
	const_mem_fun<DofMoFEMEntity,EntityHandle,&DofMoFEMEntity::get_ent>
      > >
  > > DofMoFEMEntity_multiIndex;

typedef multi_index_container<
  const DofMoFEMEntity*,
  indexed_by<
    ordered_unique< 
      member<DofMoFEMEntity,const UId,&DofMoFEMEntity::uid> >
  > > DofMoFEMEntity_multiIndex_uid_view;

/** 
 * @relates multi_index_container
 * \brief MultiIndex container keeps FEDofMoFEMEntity
 *
 * \param ordered_unique< 
 *     tag<Unique_mi_tag>, const_mem_fun<FEDofMoFEMEntity::interface_type_DofMoFEMEntity,UId,&FEDofMoFEMEntity::get_unique_id> >,
 * \param ordered_non_unique<
 *    tag<MoABEnt_mi_tag>, const_mem_fun<FEDofMoFEMEntity::interface_type_DofMoFEMEntity,EntityHandle,&FEDofMoFEMEntity::get_ent> >,
 * \param ordered_non_unique<
 *    tag<FieldName_mi_tag>, const_mem_fun<FEDofMoFEMEntity::interface_type_MoFEMField,string,&FEDofMoFEMEntity::get_name> >,
 * \param ordered_non_unique<
 *    tag<Composite_mi_tag>, <br>
 *     composite_key<  
 *	FEDofMoFEMEntity,  <br>
 *	  const_mem_fun<FEDofMoFEMEntity::interface_type_MoFEMField,string,&FEDofMoFEMEntity::get_name>,  <br>
 *	  const_mem_fun<FEDofMoFEMEntity::interface_type_RefMoFEMEntity,EntityType,&FEDofMoFEMEntity::get_ent_type>,  <br>
 *	  key_from_key< <br>
 *	    member<SideNumber,int,&SideNumber::side_number>,  <br>
 *	    member<FEDofMoFEMEntity::BaseFEDofMoFEMEntity,SideNumber *,&FEDofMoFEMEntity::side_number_ptr>
 *	  >
 *     > >,
 * \param ordered_non_unique<
 *     tag<Composite_mi_tag2>,  <br>
 *     composite_key<
 *	FEDofMoFEMEntity, <br> 
 *	  const_mem_fun<FEDofMoFEMEntity::interface_type_MoFEMField,boost::string_ref,&FEDofMoFEMEntity::get_name>,  <br>
 *	  const_mem_fun<FEDofMoFEMEntity::interface_type_RefMoFEMEntity,EntityType,&FEDofMoFEMEntity::get_ent_type>  <br>
 *	> >,
 * \param ordered_non_unique<
 *     tag<Composite_mi_tag3>,  <br>
 *     composite_key<
 *	FEDofMoFEMEntity,  <br>
 *	  const_mem_fun<FEDofMoFEMEntity::interface_type_MoFEMField,boost::string_ref,&FEDofMoFEMEntity::get_name>,  <br>
 *	  const_mem_fun<FEDofMoFEMEntity::interface_type_DofMoFEMEntity,EntityHandle,&FEDofMoFEMEntity::get_ent>
 *	> >
 */
typedef multi_index_container<
  FEDofMoFEMEntity,
  indexed_by<
    ordered_unique< 
      tag<Unique_mi_tag>, const_mem_fun<FEDofMoFEMEntity::interface_type_DofMoFEMEntity,UId,&FEDofMoFEMEntity::get_unique_id> >,
    ordered_non_unique<
      tag<MoABEnt_mi_tag>, const_mem_fun<FEDofMoFEMEntity::interface_type_DofMoFEMEntity,EntityHandle,&FEDofMoFEMEntity::get_ent> >,
    ordered_non_unique<
      tag<FieldName_mi_tag>, const_mem_fun<FEDofMoFEMEntity::interface_type_MoFEMField,boost::string_ref,&FEDofMoFEMEntity::get_name_ref> >,
    ordered_non_unique<
      tag<Composite_mi_tag>, 
      composite_key<
	FEDofMoFEMEntity,
	  const_mem_fun<FEDofMoFEMEntity::interface_type_MoFEMField,boost::string_ref,&FEDofMoFEMEntity::get_name_ref>,
	  const_mem_fun<FEDofMoFEMEntity::interface_type_RefMoFEMEntity,EntityType,&FEDofMoFEMEntity::get_ent_type>,
	  key_from_key<
	    member<SideNumber,int,&SideNumber::side_number>,
	    member<FEDofMoFEMEntity::BaseFEDofMoFEMEntity,SideNumber *,&FEDofMoFEMEntity::side_number_ptr>
	  >
      > >,
    ordered_non_unique<
      tag<Composite_mi_tag2>, 
      composite_key<
	FEDofMoFEMEntity,
	  const_mem_fun<FEDofMoFEMEntity::interface_type_MoFEMField,boost::string_ref,&FEDofMoFEMEntity::get_name_ref>,
	  const_mem_fun<FEDofMoFEMEntity::interface_type_RefMoFEMEntity,EntityType,&FEDofMoFEMEntity::get_ent_type>
	> >,
    ordered_non_unique<
      tag<Composite_mi_tag3>, 
      composite_key<
	FEDofMoFEMEntity,
	  const_mem_fun<FEDofMoFEMEntity::interface_type_MoFEMField,boost::string_ref,&FEDofMoFEMEntity::get_name_ref>,
	  const_mem_fun<FEDofMoFEMEntity::interface_type_DofMoFEMEntity,EntityHandle,&FEDofMoFEMEntity::get_ent>
	> >
  > > FEDofMoFEMEntity_multiIndex;

/** 
 * @relates multi_index_container
 * \brief MultiIndex container keeps FENumeredDofMoFEMEntity
 *
 * \param ordered_unique< 
 *     tag<Unique_mi_tag>, const_mem_fun<FEDofMoFEMEntity::interface_type_DofMoFEMEntity,UId,&FEDofMoFEMEntity::get_unique_id> >,
 * \param ordered_non_unique<
 *    tag<MoABEnt_mi_tag>, const_mem_fun<FEDofMoFEMEntity::interface_type_DofMoFEMEntity,EntityHandle,&FEDofMoFEMEntity::get_ent> >,
 * \param ordered_non_unique<
 *    tag<FieldName_mi_tag>, const_mem_fun<FEDofMoFEMEntity::interface_type_MoFEMField,string,&FEDofMoFEMEntity::get_name> >,
 * \param ordered_non_unique<
 *    tag<Composite_mi_tag>, <br>
 *     composite_key<  
 *	FEDofMoFEMEntity,  <br>
 *	  const_mem_fun<FEDofMoFEMEntity::interface_type_MoFEMField,boost::string_ref,&FEDofMoFEMEntity::get_name>,  <br>
 *	  const_mem_fun<FEDofMoFEMEntity::interface_type_RefMoFEMEntity,EntityType,&FEDofMoFEMEntity::get_ent_type>,  <br>
 *	  key_from_key< <br>
 *	    member<SideNumber,int,&SideNumber::side_number>,  <br>
 *	    member<FEDofMoFEMEntity::BaseFEDofMoFEMEntity,SideNumber *,&FEDofMoFEMEntity::side_number_ptr>
 *	  >
 *     > >,
 * \param ordered_non_unique<
 *     tag<Composite_mi_tag2>,  <br>
 *     composite_key<
 *	FEDofMoFEMEntity, <br> 
 *	  const_mem_fun<FEDofMoFEMEntity::interface_type_MoFEMField,boost::string_ref,&FEDofMoFEMEntity::get_name>,  <br>
 *	  const_mem_fun<FEDofMoFEMEntity::interface_type_RefMoFEMEntity,EntityType,&FEDofMoFEMEntity::get_ent_type>  <br>
 *	> >,
 * \param ordered_non_unique<
 *     tag<Composite_mi_tag3>,  <br>
 *     composite_key<
 *	FEDofMoFEMEntity,  <br>
 *	  const_mem_fun<FEDofMoFEMEntity::interface_type_MoFEMField,boost::string_ref,&FEDofMoFEMEntity::get_name>,  <br>
 *	  const_mem_fun<FEDofMoFEMEntity::interface_type_DofMoFEMEntity,EntityHandle,&FEDofMoFEMEntity::get_ent>
 *	> >
 */
typedef multi_index_container<
  FENumeredDofMoFEMEntity,
  indexed_by<
    ordered_unique< 
      tag<Unique_mi_tag>, const_mem_fun<FENumeredDofMoFEMEntity::interface_type_DofMoFEMEntity,UId,&FENumeredDofMoFEMEntity::get_unique_id> >,
    ordered_non_unique<
      tag<MoABEnt_mi_tag>, const_mem_fun<FENumeredDofMoFEMEntity::interface_type_DofMoFEMEntity,EntityHandle,&FENumeredDofMoFEMEntity::get_ent> >,
    ordered_non_unique<
      tag<FieldName_mi_tag>, const_mem_fun<FENumeredDofMoFEMEntity::interface_type_MoFEMField,boost::string_ref,&FENumeredDofMoFEMEntity::get_name_ref> >,
    ordered_non_unique<
      tag<Composite_mi_tag>, 
      composite_key<
	FENumeredDofMoFEMEntity,
	  const_mem_fun<FENumeredDofMoFEMEntity::interface_type_MoFEMField,boost::string_ref,&FENumeredDofMoFEMEntity::get_name_ref>,
	  const_mem_fun<FENumeredDofMoFEMEntity::interface_type_RefMoFEMEntity,EntityType,&FENumeredDofMoFEMEntity::get_ent_type>,
	  key_from_key<
	    member<SideNumber,int,&SideNumber::side_number>,
	    member<FENumeredDofMoFEMEntity::BaseFEDofMoFEMEntity,SideNumber*,&FENumeredDofMoFEMEntity::side_number_ptr>
	  >
      > >,
    ordered_non_unique<
      tag<Composite_mi_tag2>, 
      composite_key<
	FENumeredDofMoFEMEntity,
	  const_mem_fun<FENumeredDofMoFEMEntity::interface_type_MoFEMField,boost::string_ref,&FENumeredDofMoFEMEntity::get_name_ref>,
	  const_mem_fun<FENumeredDofMoFEMEntity::interface_type_RefMoFEMEntity,EntityType,&FENumeredDofMoFEMEntity::get_ent_type>
	> >,
    ordered_non_unique<
      tag<Composite_mi_tag3>, 
      composite_key<
	FENumeredDofMoFEMEntity,
	  const_mem_fun<FENumeredDofMoFEMEntity::interface_type_MoFEMField,boost::string_ref,&FENumeredDofMoFEMEntity::get_name_ref>,
	  const_mem_fun<FENumeredDofMoFEMEntity::interface_type_DofMoFEMEntity,EntityHandle,&FENumeredDofMoFEMEntity::get_ent>
	> >
  > > FENumeredDofMoFEMEntity_multiIndex;

/** 
 * @relates multi_index_container
 * \brief MultiIndex container keeps NumeredDofMoFEMEntity
 *
 * \param    ordered_unique< 
      tag<Unique_mi_tag>, const_mem_fun<NumeredDofMoFEMEntity::interface_type_DofMoFEMEntity,UId,&NumeredDofMoFEMEntity::get_unique_id> >,
 * \param    ordered_unique< 
      tag<Idx_mi_tag>, member<NumeredDofMoFEMEntity,DofIdx,&NumeredDofMoFEMEntity::dof_idx> >,
 * \param    ordered_non_unique<
      tag<FieldName_mi_tag>, const_mem_fun<NumeredDofMoFEMEntity::interface_type_MoFEMField,boost::string_ref,&NumeredDofMoFEMEntity::get_name> >,
 * \param    ordered_non_unique< 
      tag<PetscGlobalIdx_mi_tag>, member<NumeredDofMoFEMEntity,DofIdx,&NumeredDofMoFEMEntity::petsc_gloabl_dof_idx> >,
 * \param    ordered_non_unique< 
      tag<PetscLocalIdx_mi_tag>, member<NumeredDofMoFEMEntity,DofIdx,&NumeredDofMoFEMEntity::petsc_local_dof_idx> >,
 * \param    ordered_non_unique< 
      tag<Part_mi_tag>, member<NumeredDofMoFEMEntity,unsigned int,&NumeredDofMoFEMEntity::part> >,
 * \param    ordered_non_unique<
      tag<MoABEnt_mi_tag>, const_mem_fun<NumeredDofMoFEMEntity::interface_type_DofMoFEMEntity,EntityHandle,&NumeredDofMoFEMEntity::get_ent> >
 *
 */
typedef multi_index_container<
  NumeredDofMoFEMEntity,
  indexed_by<
    ordered_unique< 
      tag<Unique_mi_tag>, const_mem_fun<NumeredDofMoFEMEntity::interface_type_DofMoFEMEntity,UId,&NumeredDofMoFEMEntity::get_unique_id> >,
    ordered_unique< 
      tag<Idx_mi_tag>, member<NumeredDofMoFEMEntity,DofIdx,&NumeredDofMoFEMEntity::dof_idx> >,
    ordered_non_unique<
      tag<FieldName_mi_tag>, const_mem_fun<NumeredDofMoFEMEntity::interface_type_MoFEMField,boost::string_ref,&NumeredDofMoFEMEntity::get_name_ref> >,
    ordered_non_unique< 
      tag<PetscGlobalIdx_mi_tag>, member<NumeredDofMoFEMEntity,DofIdx,&NumeredDofMoFEMEntity::petsc_gloabl_dof_idx> >,
    ordered_non_unique< 
      tag<PetscLocalIdx_mi_tag>, member<NumeredDofMoFEMEntity,DofIdx,&NumeredDofMoFEMEntity::petsc_local_dof_idx> >,
    ordered_non_unique< 
      tag<Part_mi_tag>, member<NumeredDofMoFEMEntity,unsigned int,&NumeredDofMoFEMEntity::part> >,
    ordered_non_unique<
      tag<MoABEnt_mi_tag>, const_mem_fun<NumeredDofMoFEMEntity::interface_type_DofMoFEMEntity,EntityHandle,&NumeredDofMoFEMEntity::get_ent> >
  > > NumeredDofMoFEMEntity_multiIndex;

typedef multi_index_container<
  const NumeredDofMoFEMEntity*,
  indexed_by<
    ordered_unique< 
      const_mem_fun<NumeredDofMoFEMEntity::interface_type_DofMoFEMEntity,UId,&NumeredDofMoFEMEntity::get_unique_id> >
  > > NumeredDofMoFEMEntity_multiIndex_uid_view;

/** 
 * \brief Finite element definition
 */
struct MoFEMFiniteElement {
  EntityHandle meshset; ///< meshset stores FE ents 
  BitFEId* tag_id_data; ///< ptr to tag storing FE id
  void* tag_name_data; ///< ptr to tag storing FE name
  int tag_name_size; ///< numer of characters in FE name
  BitFieldId* tag_BitFieldId_col_data; ///< tag stores col id_id for fields
  BitFieldId* tag_BitFieldId_row_data;  ///< tag stores row id_id for fields
  BitFieldId* tag_BitFieldId_data; ///< tag stores data id_id for fields
  Tag th_FEMatData,th_FEVecData;
  Tag th_DofUidRow,th_DofUidCol,th_DofUidData;
  MoFEMFiniteElement(Interface &moab,const EntityHandle _meshset);
  inline BitFEId get_id() const { return *tag_id_data; };
  /// get meshset
  inline EntityHandle get_meshset() const { return meshset; }
  /// get FE name
  inline boost::string_ref get_name_ref() const { return boost::string_ref((char *)tag_name_data,tag_name_size); }
  inline string get_name() const { return string((char *)tag_name_data,tag_name_size); }
  /// get BitFieldId col
  inline BitFieldId get_BitFieldId_col() const { return *((BitFieldId*)tag_BitFieldId_col_data); }
  /// get BitFieldId row
  inline BitFieldId get_BitFieldId_row() const { return *((BitFieldId*)tag_BitFieldId_row_data); }
  /// get BitFieldId data
  inline BitFieldId get_BitFieldId_data() const { return *((BitFieldId*)tag_BitFieldId_data); }
  friend ostream& operator<<(ostream& os,const MoFEMFiniteElement& e);
};

/**
 * \brief Inetface for FE
 */
template <typename T>
struct interface_MoFEMFiniteElement {
  const T *fe_ptr;
  interface_MoFEMFiniteElement(const T *_ptr): fe_ptr(_ptr) {};
  inline BitFEId get_id() const { return fe_ptr->get_id(); }
  inline EntityHandle get_meshset() const { return fe_ptr->get_meshset(); }
  inline boost::string_ref get_name_ref() const { return fe_ptr->get_name_ref(); }
  inline string get_name() const { return fe_ptr->get_name(); }
  inline BitFieldId get_BitFieldId_col() const { return fe_ptr->get_BitFieldId_col(); }
  inline BitFieldId get_BitFieldId_row() const { return fe_ptr->get_BitFieldId_row(); }
  inline BitFieldId get_BitFieldId_data() const { return fe_ptr->get_BitFieldId_data(); }
};

/**
 * \brief Finite element data for entitiy
 */
struct EntMoFEMFiniteElement: public interface_MoFEMFiniteElement<MoFEMFiniteElement>,interface_RefMoFEMElement<RefMoFEMElement> {
  typedef interface_RefMoFEMEntity<RefMoFEMElement> interface_type_RefMoFEMEntity;
  typedef interface_RefMoFEMElement<RefMoFEMElement> interface_type_RefMoFEMElement;
  typedef interface_MoFEMFiniteElement<MoFEMFiniteElement> interface_type_MoFEMFiniteElement;
  const UId* tag_row_uids_data;
  int tag_row_uids_size;
  const UId* tag_col_uids_data;
  int tag_col_uids_size;
  const UId* tag_data_uids_data;
  int tag_data_uids_size;
  FEDofMoFEMEntity_multiIndex data_dofs;
  EntMoFEMFiniteElement(Interface &moab,const RefMoFEMElement *_ref_MoFEMFiniteElement,const MoFEMFiniteElement *_MoFEMFiniteElement_ptr);
  inline EntityHandle get_ent() const { return get_ref_ent(); }
  inline DofIdx get_nb_dofs_row() const { return tag_row_uids_size/sizeof(UId); }
  inline DofIdx get_nb_dofs_col() const { return tag_col_uids_size/sizeof(UId); }
  inline DofIdx get_nb_dofs_data() const { return tag_data_uids_size/sizeof(UId); }
  friend ostream& operator<<(ostream& os,const EntMoFEMFiniteElement& e);
  PetscErrorCode get_MoFEMFiniteElement_row_dof_uid_view(
    const DofMoFEMEntity_multiIndex &dofs,DofMoFEMEntity_multiIndex_uid_view &dofs_view,
    const int operation_type = Interface::UNION) const;
  PetscErrorCode get_MoFEMFiniteElement_col_dof_uid_view(
    const DofMoFEMEntity_multiIndex &dofs,DofMoFEMEntity_multiIndex_uid_view &dofs_view,
    const int operation_type = Interface::UNION) const;
  PetscErrorCode get_MoFEMFiniteElement_row_dof_uid_view(
    const NumeredDofMoFEMEntity_multiIndex &dofs,NumeredDofMoFEMEntity_multiIndex_uid_view &dofs_view,
    const int operation_type = Interface::UNION) const;
  PetscErrorCode get_MoFEMFiniteElement_col_dof_uid_view(
    const NumeredDofMoFEMEntity_multiIndex &dofs,NumeredDofMoFEMEntity_multiIndex_uid_view &dofs_view,
    const int operation_type = Interface::UNION) const;
  //
  PetscErrorCode get_uid_side_number(
    Interface &moab,const UId uid,
    const DofMoFEMEntity_multiIndex &dofs_moabfield,
    int &side_number, int &sense, int &offset) const;
};

/**
 * \brief interface to EntMoFEMFiniteElement
 */
template <typename T>
struct interface_EntMoFEMFiniteElement:public interface_MoFEMFiniteElement<T>,interface_RefMoFEMElement<T> {
  interface_EntMoFEMFiniteElement(const T *_ptr): interface_MoFEMFiniteElement<T>(_ptr),interface_RefMoFEMElement<T>(_ptr) {};
  inline EntityID get_ent_id() const { return interface_MoFEMFiniteElement<T>::fe_ptr->get_ent_id(); }
  inline EntityType get_ent_type() const { return interface_MoFEMFiniteElement<T>::fe_ptr->get_ent_type(); }
  //
  inline DofIdx get_nb_dofs_row() const { return interface_MoFEMFiniteElement<T>::fe_ptr->get_nb_dofs_row(); }
  inline DofIdx get_nb_dofs_col() const { return interface_MoFEMFiniteElement<T>::fe_ptr->get_nb_dofs_col(); }
  inline DofIdx get_nb_dofs_data() const { return interface_MoFEMFiniteElement<T>::fe_ptr->get_nb_dofs_data(); }
  inline EntityHandle get_ent() const { return interface_MoFEMFiniteElement<T>::fe_ptr->get_ref_ent(); };
  //
  SideNumber_multiIndex &get_side_number_table() const { return interface_MoFEMFiniteElement<T>::fe_ptr->get_side_number_table(); }
  SideNumber* get_side_number_ptr(Interface &moab,EntityHandle ent) const { return interface_MoFEMFiniteElement<T>::fe_ptr->get_side_number_ptr(moab,ent); }
};

/// Partitioned Finite Element in Problem
struct NumeredMoFEMFiniteElement: public interface_EntMoFEMFiniteElement<EntMoFEMFiniteElement> {
  typedef interface_MoFEMFiniteElement<EntMoFEMFiniteElement> interface_type_MoFEMFiniteElement;
  typedef interface_EntMoFEMFiniteElement<EntMoFEMFiniteElement> interface_type_EntMoFEMFiniteElement;
  unsigned int part;
  FENumeredDofMoFEMEntity_multiIndex rows_dofs;
  FENumeredDofMoFEMEntity_multiIndex cols_dofs;
  NumeredMoFEMFiniteElement(const EntMoFEMFiniteElement *EntMoFEMFiniteElement_ptr): interface_EntMoFEMFiniteElement<EntMoFEMFiniteElement>(EntMoFEMFiniteElement_ptr),part(-1) {};
  inline unsigned int get_part() const { return part; };
  friend ostream& operator<<(ostream& os,const NumeredMoFEMFiniteElement& e) {
    os << "part " << e.part << " " << *(e.fe_ptr);
    return os;
  }
};

/// interface for NumeredMoFEMFiniteElement
template <typename T>
struct interface_NumeredMoFEMFiniteElement: public interface_EntMoFEMFiniteElement<T> {
  const T *ptr;
  interface_NumeredMoFEMFiniteElement(const T *_ptr): interface_EntMoFEMFiniteElement<T>(_ptr) {};
  inline unsigned int get_part() const { return ptr->get_part(); }
};

/**
 * @relates multi_index_container
 * \brief MultiIndex container for finite element enetities
 *
 * \param  hashed_unique<
      tag<Composite_unique_mi_tag>,  
      composite_key< 
	EntMoFEMFiniteElement, <br>
	const_mem_fun<EntMoFEMFiniteElement::interface_type_MoFEMFiniteElement,EntityHandle,&EntMoFEMFiniteElement::get_meshset>,
	const_mem_fun<EntMoFEMFiniteElement,EntityHandle,&EntMoFEMFiniteElement::get_ent> > >,
 * \param    ordered_non_unique<
      tag<MoABEnt_mi_tag>, <br> const_mem_fun<EntMoFEMFiniteElement,EntityHandle,&EntMoFEMFiniteElement::get_ent> >,
 * \param    ordered_non_unique<
      tag<MoFEMFiniteElement_name_mi_tag>, <br> const_mem_fun<EntMoFEMFiniteElement::interface_type_MoFEMFiniteElement,boost::string_ref,&EntMoFEMFiniteElement::get_name_ref> >,
 * \param    ordered_non_unique<
      tag<BitFEId_mi_tag>, <br> const_mem_fun<EntMoFEMFiniteElement::interface_type_MoFEMFiniteElement,BitFEId,&EntMoFEMFiniteElement::get_id>, ltbit<BitFEId> >,
 * \param    ordered_non_unique<
      tag<EntType_mi_tag>, <br> const_mem_fun<EntMoFEMFiniteElement::interface_type_RefMoFEMEntity,EntityType,&EntMoFEMFiniteElement::get_ent_type> >,
 * \param    ordered_non_unique<
      tag<Composite_mi_tag>, 
      composite_key<
	EntMoFEMFiniteElement, <br>
	const_mem_fun<EntMoFEMFiniteElement,EntityHandle,&EntMoFEMFiniteElement::get_ent>,
	const_mem_fun<EntMoFEMFiniteElement::interface_type_MoFEMFiniteElement,boost::string_ref,&EntMoFEMFiniteElement::get_name_ref> > >
 */
typedef multi_index_container<
  EntMoFEMFiniteElement,
  indexed_by<
    hashed_unique<
      tag<Composite_unique_mi_tag>,       
      composite_key<
	EntMoFEMFiniteElement,
	const_mem_fun<EntMoFEMFiniteElement::interface_type_MoFEMFiniteElement,EntityHandle,&EntMoFEMFiniteElement::get_meshset>,
	const_mem_fun<EntMoFEMFiniteElement,EntityHandle,&EntMoFEMFiniteElement::get_ent> > >,
    ordered_non_unique<
      tag<MoABEnt_mi_tag>, const_mem_fun<EntMoFEMFiniteElement,EntityHandle,&EntMoFEMFiniteElement::get_ent> >,
    ordered_non_unique<
      tag<MoFEMFiniteElement_name_mi_tag>, const_mem_fun<EntMoFEMFiniteElement::interface_type_MoFEMFiniteElement,boost::string_ref,&EntMoFEMFiniteElement::get_name_ref> >,
    ordered_non_unique<
      tag<BitFEId_mi_tag>, const_mem_fun<EntMoFEMFiniteElement::interface_type_MoFEMFiniteElement,BitFEId,&EntMoFEMFiniteElement::get_id>, ltbit<BitFEId> >,
    ordered_non_unique<
      tag<EntType_mi_tag>, const_mem_fun<EntMoFEMFiniteElement::interface_type_RefMoFEMEntity,EntityType,&EntMoFEMFiniteElement::get_ent_type> >,
    ordered_non_unique<
      tag<Composite_mi_tag>, 
      composite_key<
	EntMoFEMFiniteElement,
	const_mem_fun<EntMoFEMFiniteElement,EntityHandle,&EntMoFEMFiniteElement::get_ent>,
	const_mem_fun<EntMoFEMFiniteElement::interface_type_MoFEMFiniteElement,boost::string_ref,&EntMoFEMFiniteElement::get_name_ref> > >
  > > EntMoFEMFiniteElement_multiIndex;

/** 
 * @relates multi_index_container
 * \brief MultiIndex for entities uin the problem.
 *
 * \param   hashed_unique<
      tag<Composite_unique_mi_tag>,       
      composite_key< <br>
	NumeredMoFEMFiniteElement,
	const_mem_fun<NumeredMoFEMFiniteElement::interface_type_MoFEMFiniteElement,EntityHandle,&NumeredMoFEMFiniteElement::get_meshset>,
	const_mem_fun<NumeredMoFEMFiniteElement::interface_type_EntMoFEMFiniteElement,EntityHandle,&NumeredMoFEMFiniteElement::get_ent> > >,
 * \param    ordered_non_unique<
      tag<MoFEMFiniteElement_name_mi_tag>, <br> const_mem_fun<NumeredMoFEMFiniteElement::interface_type_MoFEMFiniteElement,string,&NumeredMoFEMFiniteElement::get_name> >,
 * \param    ordered_non_unique<
      tag<MoFEMFiniteElement_Part_mi_tag>, <br> member<NumeredMoFEMFiniteElement,unsigned int,&NumeredMoFEMFiniteElement::part> >,
 * \param    ordered_non_unique<
      tag<Composite_mi_tag>,       
      composite_key< <br>
	NumeredMoFEMFiniteElement,
	const_mem_fun<NumeredMoFEMFiniteElement::interface_type_MoFEMFiniteElement,boost::string_ref,&NumeredMoFEMFiniteElement::get_name_ref>,
	member<NumeredMoFEMFiniteElement,unsigned int,&NumeredMoFEMFiniteElement::part> > >
 */
typedef multi_index_container<
  NumeredMoFEMFiniteElement,
  indexed_by<
    hashed_unique<
      tag<Composite_unique_mi_tag>,       
      composite_key<
	NumeredMoFEMFiniteElement,
	const_mem_fun<NumeredMoFEMFiniteElement::interface_type_MoFEMFiniteElement,EntityHandle,&NumeredMoFEMFiniteElement::get_meshset>,
	const_mem_fun<NumeredMoFEMFiniteElement::interface_type_EntMoFEMFiniteElement,EntityHandle,&NumeredMoFEMFiniteElement::get_ent> > >,
    ordered_non_unique<
      tag<MoFEMFiniteElement_name_mi_tag>, const_mem_fun<NumeredMoFEMFiniteElement::interface_type_MoFEMFiniteElement,boost::string_ref,&NumeredMoFEMFiniteElement::get_name_ref> >,
    ordered_non_unique<
      tag<MoFEMFiniteElement_Part_mi_tag>, member<NumeredMoFEMFiniteElement,unsigned int,&NumeredMoFEMFiniteElement::part> >,
    ordered_non_unique<
      tag<Composite_mi_tag>,       
      composite_key<
	NumeredMoFEMFiniteElement,
	const_mem_fun<NumeredMoFEMFiniteElement::interface_type_MoFEMFiniteElement,boost::string_ref,&NumeredMoFEMFiniteElement::get_name_ref>,
	member<NumeredMoFEMFiniteElement,unsigned int,&NumeredMoFEMFiniteElement::part> > >
  > > NumeredMoFEMFiniteElement_multiIndex;

/// \brief keeps data about problem
struct MoFEMProblem {
  EntityHandle meshset;
  BitProblemId* tag_id_data;
  const void* tag_name_data;
  int tag_name_size;
  DofIdx* tag_nbdof_data_row;
  DofIdx* tag_nbdof_data_col;
  DofIdx* tag_local_nbdof_data_row;
  DofIdx* tag_local_nbdof_data_col;
  DofIdx* tag_ghost_nbdof_data_row;
  DofIdx* tag_ghost_nbdof_data_col;
  BitFEId* tag_BitFEId_data;
  BitRefLevel* tag_BitRefLevel;
  NumeredDofMoFEMEntity_multiIndex numered_dofs_rows;
  NumeredDofMoFEMEntity_multiIndex numered_dofs_cols;
  NumeredMoFEMFiniteElement_multiIndex numered_finite_elements;

  /**
    * use with loops to iterate problem fes 
    *
    * for(_IT_NUMEREDFEMOFEMENTITY_FOR_LOOP_(MOFEMPROBLEM,NAME,IT)) {
    *   ...
    * }
    *
    */
  #define _IT_NUMEREDFEMOFEMENTITY_BY_NAME_FOR_LOOP_(MOFEMPROBLEM,NAME,IT) \
    NumeredMoFEMFiniteElement_multiIndex::index<MoFEMFiniteElement_name_mi_tag>::type::iterator IT = MOFEMPROBLEM->get_numered_fes_begin(NAME); \
    IT!=MOFEMPROBLEM->get_numered_fes_begin(NAME); IT++

  NumeredMoFEMFiniteElement_multiIndex::index<MoFEMFiniteElement_name_mi_tag>::type::iterator get_numered_fes_begin(string fe_name) const { 
    return numered_finite_elements.get<MoFEMFiniteElement_name_mi_tag>().lower_bound(fe_name);
  }

  NumeredMoFEMFiniteElement_multiIndex::index<MoFEMFiniteElement_name_mi_tag>::type::iterator get_numered_fes_end(string fe_name) const { 
    return numered_finite_elements.get<MoFEMFiniteElement_name_mi_tag>().upper_bound(fe_name);
  }

  /**
    * use with loops to iterate problem fes 
    *
    * for(_IT_NUMEREDFEMOFEMENTITY_BY_NAME_AND_PART_FOR_LOOP_(MOFEMPROBLEM,NAME,PART,IT)) {
    *   ...
    * }
    *
    */
  #define _IT_NUMEREDFEMOFEMENTITY_BY_NAME_AND_PART_FOR_LOOP_(MOFEMPROBLEM,NAME,PART,IT) \
    NumeredMoFEMFiniteElement_multiIndex::index<Composite_mi_tag>::type::iterator IT = MOFEMPROBLEM->get_numered_fes_begin(NAME,PART); \
    IT!=MOFEMPROBLEM->get_numered_fes_begin(NAME,PART); IT++

  NumeredMoFEMFiniteElement_multiIndex::index<Composite_mi_tag>::type::iterator get_numered_fes_begin(string fe_name,int part) const { 
    return numered_finite_elements.get<Composite_mi_tag>().lower_bound(boost::make_tuple(fe_name,part));
  }

  NumeredMoFEMFiniteElement_multiIndex::index<Composite_mi_tag>::type::iterator get_numered_fes_end(string fe_name,int part) const { 
    return numered_finite_elements.get<Composite_mi_tag>().upper_bound(boost::make_tuple(fe_name,part));
  }

  /**
    * use with loops to iterate row dofs 
    *
    * for(_IT_NUMEREDDOFMOFEMENTITY_ROW_FOR_LOOP_(MOFEMPROBLEM,IT)) {
    *   ...
    * }
    *
    */
  #define _IT_NUMEREDDOFMOFEMENTITY_ROW_FOR_LOOP_(MOFEMPROBLEM,IT) \
    NumeredDofMoFEMEntity_multiIndex::iterator IT = MOFEMPROBLEM->get_numered_dofs_rows_begin(); \
    IT!=MOFEMPROBLEM->get_numered_dofs_rows_end(); IT++

  /**
    * use with loops to iterate row dofs 
    *
    * for(_IT_NUMEREDDOFMOFEMENTITY_COL_FOR_LOOP_(MOFEMPROBLEM,IT)) {
    *   ...
    * }
    *
    */
  #define _IT_NUMEREDDOFMOFEMENTITY_COL_FOR_LOOP_(MOFEMPROBLEM,IT) \
    NumeredDofMoFEMEntity_multiIndex::iterator IT = MOFEMPROBLEM->get_numered_dofs_cols_begin(); \
    IT!=MOFEMPROBLEM->get_numered_dofs_cols_end(); IT++

  /// get begin iterator for numered_dofs_rows (insted you can use #_IT_NUMEREDDOFMOFEMENTITY_ROW_FOR_LOOP_ for loops)
  NumeredDofMoFEMEntity_multiIndex::iterator get_numered_dofs_rows_begin() const { return numered_dofs_rows.begin(); }

  /// get end iterator for numered_dofs_rows (insted you can use #_IT_NUMEREDDOFMOFEMENTITY_ROW_FOR_LOOP_ for loops)
  NumeredDofMoFEMEntity_multiIndex::iterator get_numered_dofs_rows_end() const { return numered_dofs_rows.end(); }

  /// get begin iterator for numered_dofs_cols (insted you can use #_IT_NUMEREDDOFMOFEMENTITY_COL_FOR_LOOP_ for loops)
  NumeredDofMoFEMEntity_multiIndex::iterator get_numered_dofs_cols_begin() const { return numered_dofs_cols.begin(); }

  /// get end iterator for numered_dofs_cols (insted you can use #_IT_NUMEREDDOFMOFEMENTITY_COL_FOR_LOOP_ for loops)
  NumeredDofMoFEMEntity_multiIndex::iterator get_numered_dofs_cols_end() const { return numered_dofs_cols.end(); }

  /**
    * get iterator of dof in row by uid
    */
  #define _IT_NUMEREDDOFMOFEMENTITY_ROW_BY_UID_(MOFEMPROBLEM,UID,IT) \
    NumeredDofMoFEMEntity_multiIndex::index<Unique_mi_tag>::type::iterator IT = MOFEMPROBLEM->get_row_dof_by_uid(UID);

  /**
    * get iterator of dof in col by uid
    */
  #define _IT_NUMEREDDOFMOFEMENTITY_COL_BY_UID_(MOFEMPROBLEM,UID,IT) \
    NumeredDofMoFEMEntity_multiIndex::index<Unique_mi_tag>::type::iterator IT = MOFEMPROBLEM->get_col_dof_by_uid(UID);

  /// get iterator of dof in row by uid (instead you can use #_IT_NUMEREDDOFMOFEMENTITY_ROW_BY_UID_FOR_LOOP_)
  NumeredDofMoFEMEntity_multiIndex::index<Unique_mi_tag>::type::iterator get_row_dof_by_uid(UId uid) const { return numered_dofs_rows.get<Unique_mi_tag>().find(uid); };

  /// get iterator of dof in column by uid (instead you can use #_IT_NUMEREDDOFMOFEMENTITY_COL_BY_UID_FOR_LOOP_)
  NumeredDofMoFEMEntity_multiIndex::index<Unique_mi_tag>::type::iterator get_col_dof_by_uid(UId uid) const { return numered_dofs_cols.get<Unique_mi_tag>().find(uid); };

  /**
    * use with loops to iterate row dofs 
    *
    * for(_IT_NUMEREDDOFMOFEMENTITY_BY_LOCIDX_ROW_FOR_LOOP_(MOFEMPROBLEM,IT)) {
    *   ...
    * }
    *
    */
  #define _IT_NUMEREDDOFMOFEMENTITY_ROW_BY_LOCIDX_FOR_LOOP_(MOFEMPROBLEM,IT) \
    NumeredDofMoFEMEntity_multiIndex::index<PetscLocalIdx_mi_tag>::type::iterator IT = MOFEMPROBLEM->get_numered_dofs_rows_by_locidx_begin(0); \
    IT!=MOFEMPROBLEM->get_numered_dofs_rows_by_locidx_end(MOFEMPROBLEM->get_nb_local_dofs_row()); IT++

  /**
    * use with loops to iterate row dofs 
    *
    * for(_IT_NUMEREDDOFMOFEMENTITY_COL_BY_LOCIDX_FOR_LOOP_(MOFEMPROBLEM,IT)) {
    *   ...
    * }
    *
    */
  #define _IT_NUMEREDDOFMOFEMENTITY_COL_BY_LOCIDX_FOR_LOOP_(MOFEMPROBLEM,IT) \
    NumeredDofMoFEMEntity_multiIndex::index<PetscLocalIdx_mi_tag>::type::iterator IT = MOFEMPROBLEM->get_numered_dofs_cols_by_locidx_begin(0); \
    IT!=MOFEMPROBLEM->get_numered_dofs_cols_by_locidx_end(MOFEMPROBLEM->get_nb_local_dofs_row()); IT++

  /// get begin iterator for numered_dofs_rows (insted you can use #_IT_NUMEREDDOFMOFEMENTITY_ROW_FOR_LOOP_ for loops)
  NumeredDofMoFEMEntity_multiIndex::index<PetscLocalIdx_mi_tag>::type::iterator get_numered_dofs_rows_by_locidx_begin(const DofIdx locidx) const 
    { return numered_dofs_rows.get<PetscLocalIdx_mi_tag>().lower_bound(locidx); }

  /// get end iterator for numered_dofs_rows (insted you can use #_IT_NUMEREDDOFMOFEMENTITY_ROW_FOR_LOOP_ for loops)
  NumeredDofMoFEMEntity_multiIndex::index<PetscLocalIdx_mi_tag>::type::iterator get_numered_dofs_rows_by_locidx_end(const DofIdx locidx) const 
    { return numered_dofs_rows.get<PetscLocalIdx_mi_tag>().upper_bound(locidx); }

  /// get begin iterator for numered_dofs_cols (insted you can use #_IT_NUMEREDDOFMOFEMENTITY_COL_FOR_LOOP_ for loops)
  NumeredDofMoFEMEntity_multiIndex::index<PetscLocalIdx_mi_tag>::type::iterator get_numered_dofs_cols_by_locidx_begin(const DofIdx locidx) const 
    { return numered_dofs_cols.get<PetscLocalIdx_mi_tag>().lower_bound(locidx); }

  /// get end iterator for numered_dofs_cols (insted you can use #_IT_NUMEREDDOFMOFEMENTITY_COL_FOR_LOOP_ for loops)
  NumeredDofMoFEMEntity_multiIndex::index<PetscLocalIdx_mi_tag>::type::iterator get_numered_dofs_cols_by_locidx_end(const DofIdx locidx) const 
    { return numered_dofs_cols.get<PetscLocalIdx_mi_tag>().upper_bound(locidx); }

  /**
    * use with loops to iterate row dofs 
    *
    * for(_IT_NUMEREDDOFMOFEMENTITY_BY_ENT_ROW_FOR_LOOP_(MOFEMPROBLEM,IT)) {
    *   ...
    * }
    *
    */
  #define _IT_NUMEREDDOFMOFEMENTITY_ROW_BY_ENT_FOR_LOOP_(MOFEMPROBLEM,ENT,IT) \
    NumeredDofMoFEMEntity_multiIndex::index<MoABEnt_mi_tag>::type::iterator IT = MOFEMPROBLEM->get_numered_dofs_rows_by_ent_begin(ENT); \
    IT!=MOFEMPROBLEM->get_numered_dofs_rows_by_ent_end(ENT); IT++

  /**
    * use with loops to iterate row dofs 
    *
    * for(_IT_NUMEREDDOFMOFEMENTITY_COL_BY_ENT_FOR_LOOP_(MOFEMPROBLEM,IT)) {
    *   ...
    * }
    *
    */
  #define _IT_NUMEREDDOFMOFEMENTITY_COL_BY_ENT_FOR_LOOP_(MOFEMPROBLEM,ENT,IT) \
    NumeredDofMoFEMEntity_multiIndex::index<MoABEnt_mi_tag>::type::iterator IT = MOFEMPROBLEM->get_numered_dofs_cols_by_ent_begin(ENT); \
    IT!=MOFEMPROBLEM->get_numered_dofs_cols_by_ent_end(ENT); IT++

  /// get begin iterator for numered_dofs_rows (insted you can use #_IT_NUMEREDDOFMOFEMENTITY_ROW_BY_ENT_FOR_LOOP_ for loops)
  NumeredDofMoFEMEntity_multiIndex::index<MoABEnt_mi_tag>::type::iterator get_numered_dofs_rows_by_ent_begin(const EntityHandle ent) const 
    { return numered_dofs_rows.get<MoABEnt_mi_tag>().lower_bound(ent); }

  /// get end iterator for numered_dofs_rows (insted you can use #_IT_NUMEREDDOFMOFEMENTITY_ROW_BY_ENT_FOR_LOOP_ for loops)
  NumeredDofMoFEMEntity_multiIndex::index<MoABEnt_mi_tag>::type::iterator get_numered_dofs_rows_by_ent_end(const EntityHandle ent) const 
    { return numered_dofs_rows.get<MoABEnt_mi_tag>().upper_bound(ent); }

  /// get begin iterator for numered_dofs_cols (insted you can use #_IT_NUMEREDDOFMOFEMENTITY_COL_BY_ENT_FOR_LOOP_ for loops)
  NumeredDofMoFEMEntity_multiIndex::index<MoABEnt_mi_tag>::type::iterator get_numered_dofs_cols_by_ent_begin(const EntityHandle ent) const 
    { return numered_dofs_cols.get<MoABEnt_mi_tag>().lower_bound(ent); }

  /// get end iterator for numered_dofs_cols (insted you can use #_IT_NUMEREDDOFMOFEMENTITY_COL_BY_ENT_FOR_LOOP_ for loops)
  NumeredDofMoFEMEntity_multiIndex::index<MoABEnt_mi_tag>::type::iterator get_numered_dofs_cols_by_ent_end(const EntityHandle ent) const 
    { return numered_dofs_cols.get<MoABEnt_mi_tag>().upper_bound(ent); }


  MoFEMProblem(Interface &moab,const EntityHandle _meshset);
  inline BitProblemId get_id() const { return *((BitProblemId*)tag_id_data); }
  inline string get_name() const { return string((char *)tag_name_data,tag_name_size); }
  inline DofIdx get_nb_dofs_row() const { return *((DofIdx*)tag_nbdof_data_row); }
  inline DofIdx get_nb_dofs_col() const { return *((DofIdx*)tag_nbdof_data_col); }
  inline DofIdx get_nb_local_dofs_row() const { return *((DofIdx*)tag_local_nbdof_data_row); }
  inline DofIdx get_nb_local_dofs_col() const { return *((DofIdx*)tag_local_nbdof_data_col); }
  inline DofIdx get_nb_ghost_dofs_row() const { return *((DofIdx*)tag_ghost_nbdof_data_row); }
  inline DofIdx get_nb_ghost_dofs_col() const { return *((DofIdx*)tag_ghost_nbdof_data_col); }
  inline BitRefLevel get_BitRefLevel() const { return *tag_BitRefLevel; }
  BitFEId get_BitFEId() const;
  friend ostream& operator<<(ostream& os,const MoFEMProblem& e);
};

// indexes
typedef multi_index_container<
  MoFEMFiniteElement,
  indexed_by<
    hashed_unique<
      tag<MoFEMFiniteElement_Meshset_mi_tag>, member<MoFEMFiniteElement,EntityHandle,&MoFEMFiniteElement::meshset> >,
    hashed_unique<
      tag<BitFEId_mi_tag>, const_mem_fun<MoFEMFiniteElement,BitFEId,&MoFEMFiniteElement::get_id>, hashbit<BitFEId>, eqbit<BitFEId> >,
    ordered_unique<
      tag<MoFEMFiniteElement_name_mi_tag>, const_mem_fun<MoFEMFiniteElement,boost::string_ref,&MoFEMFiniteElement::get_name_ref> >
  > > MoFEMFiniteElement_multiIndex;


/** 
 * @relates multi_index_container
 * \brief MultiIndex container keeps AdjacencyMapForBasicMoFEMEntity
 *
 * \param   hashed_non_unique<
      tag<MoABEnt_mi_tag>, 
      member<AdjacencyMapForBasicMoFEMEntity::BasicMoFEMEntity,EntityHandle,&AdjacencyMapForBasicMoFEMEntity::ent> >,
 * \param    hashed_non_unique<
      tag<MoABEnt_mi_tag2>, 
      const_mem_fun<AdjacencyMapForBasicMoFEMEntity,EntityHandle,&AdjacencyMapForBasicMoFEMEntity::get_adj> >,
 * \param    ordered_non_unique<
      tag<EntType_mi_tag>, 
      const_mem_fun<AdjacencyMapForBasicMoFEMEntity,EntityType,&AdjacencyMapForBasicMoFEMEntity::get_adj_type> >,
 * \param    hashed_unique< tag<Composite_mi_tag>, <br>
      composite_key< <br>
	AdjacencyMapForBasicMoFEMEntity,
      	member<AdjacencyMapForBasicMoFEMEntity::BasicMoFEMEntity,EntityHandle,&AdjacencyMapForBasicMoFEMEntity::ent>,
	const_mem_fun<AdjacencyMapForBasicMoFEMEntity,EntityHandle,&AdjacencyMapForBasicMoFEMEntity::get_adj> > >
 *
 */
typedef multi_index_container<
  AdjacencyMapForBasicMoFEMEntity,
  indexed_by<
    hashed_non_unique<
      tag<MoABEnt_mi_tag>, 
      member<AdjacencyMapForBasicMoFEMEntity::BasicMoFEMEntity,EntityHandle,&AdjacencyMapForBasicMoFEMEntity::ent> >,
    hashed_non_unique<
      tag<MoABEnt_mi_tag2>, 
      const_mem_fun<AdjacencyMapForBasicMoFEMEntity,EntityHandle,&AdjacencyMapForBasicMoFEMEntity::get_adj> >,
    ordered_non_unique<
      tag<EntType_mi_tag>, 
      const_mem_fun<AdjacencyMapForBasicMoFEMEntity,EntityType,&AdjacencyMapForBasicMoFEMEntity::get_adj_type> >,
    hashed_unique<
      tag<Composite_mi_tag>, 
      composite_key<
	AdjacencyMapForBasicMoFEMEntity,
      	member<AdjacencyMapForBasicMoFEMEntity::BasicMoFEMEntity,EntityHandle,&AdjacencyMapForBasicMoFEMEntity::ent>,
	const_mem_fun<AdjacencyMapForBasicMoFEMEntity,EntityHandle,&AdjacencyMapForBasicMoFEMEntity::get_adj> > >
  > > AdjacencyMapForBasicMoFEMEntity_multiIndex;

// 

/**
  * \brief MoFEMAdjacencies of mofem finite element and entities
  *
  */
struct MoFEMAdjacencies {
  unsigned int by_other;
  const MoFEMEntity *MoFEMEntity_ptr; ///< field entity
  const EntMoFEMFiniteElement *EntMoFEMFiniteElement_ptr; ///< finite element entity
  MoFEMAdjacencies(const MoFEMEntity *_MoFEMEntity_ptr,const EntMoFEMFiniteElement *_EntMoFEMFiniteElement_ptr);
  inline EntityHandle get_MoFEMFiniteElement_meshset() const { return EntMoFEMFiniteElement_ptr->get_meshset(); }
  inline EntityHandle get_MoFEMFiniteElement_entity_handle() const { return EntMoFEMFiniteElement_ptr->get_ent(); }
  inline EntityHandle get_ent_meshset() const { return MoFEMEntity_ptr->get_meshset(); };
  inline EntityHandle get_ent_entity_handle() const { return MoFEMEntity_ptr->get_ent(); };
  BitFieldId get_ent_id() const { return MoFEMEntity_ptr->get_id(); }
  BitFEId get_BitFEId() const { return EntMoFEMFiniteElement_ptr->get_id(); }
  friend ostream& operator<<(ostream& os,const MoFEMAdjacencies &e);
};

/** 
 * @relates multi_index_container
 * \brief MultiIndex container keeps Adjacencies Element and dof entities adjacencies and vice veras.
 *
 * \param    hashed_unique< tag<Composite_unique_mi_tag>, <br>
      composite_key<
	MoFEMAdjacencies, <br>
	const_mem_fun<MoFEMAdjacencies,EntityHandle,&MoFEMAdjacencies::get_ent_meshset>, <br>
	const_mem_fun<MoFEMAdjacencies,EntityHandle,&MoFEMAdjacencies::get_ent_entity_handle>, <br>
	const_mem_fun<MoFEMAdjacencies,EntityHandle,&MoFEMAdjacencies::get_MoFEMFiniteElement_meshset>, <br>
	const_mem_fun<MoFEMAdjacencies,EntityHandle,&MoFEMAdjacencies::get_MoFEMFiniteElement_entity_handle> > >, <br>
 * \param    ordered_non_unique<
      tag<Composite_mi_tag>, <br>
       composite_key<
	MoFEMAdjacencies, <br>
	const_mem_fun<MoFEMAdjacencies,EntityHandle,&MoFEMAdjacencies::get_ent_meshset>, <br>
	const_mem_fun<MoFEMAdjacencies,EntityHandle,&MoFEMAdjacencies::get_ent_entity_handle> > > <br>
 *
 */
typedef multi_index_container<
  MoFEMAdjacencies,
  indexed_by<
    hashed_unique<
      tag<Composite_unique_mi_tag>,       
      composite_key<
	MoFEMAdjacencies,
	const_mem_fun<MoFEMAdjacencies,EntityHandle,&MoFEMAdjacencies::get_ent_meshset>,
	const_mem_fun<MoFEMAdjacencies,EntityHandle,&MoFEMAdjacencies::get_ent_entity_handle>,
	const_mem_fun<MoFEMAdjacencies,EntityHandle,&MoFEMAdjacencies::get_MoFEMFiniteElement_meshset>,
	const_mem_fun<MoFEMAdjacencies,EntityHandle,&MoFEMAdjacencies::get_MoFEMFiniteElement_entity_handle> > >,
    ordered_non_unique<
      tag<Composite_mi_tag>,
       composite_key<
	MoFEMAdjacencies,
	const_mem_fun<MoFEMAdjacencies,EntityHandle,&MoFEMAdjacencies::get_ent_meshset>,
	const_mem_fun<MoFEMAdjacencies,EntityHandle,&MoFEMAdjacencies::get_ent_entity_handle> > >
  > > MoFEMAdjacencies_multiIndex;

typedef multi_index_container<
  MoFEMProblem,
  indexed_by<
    ordered_unique<
      tag<Meshset_mi_tag>, member<MoFEMProblem,EntityHandle,&MoFEMProblem::meshset> >,
    hashed_unique<
      tag<BitProblemId_mi_tag>, const_mem_fun<MoFEMProblem,BitProblemId,&MoFEMProblem::get_id>, hashbit<BitProblemId>, eqbit<BitProblemId> >,
    hashed_unique<
      tag<MoFEMProblem_mi_tag>, const_mem_fun<MoFEMProblem,string,&MoFEMProblem::get_name> >
  > > MoFEMProblem_multiIndex;

//CUBIT BC DATA

/** 
 * \typedef Cubit_BC_bitset
 * bc & material meshsets
 *
 */
typedef bitset<16> Cubit_BC_bitset;
enum Cubit_BC {
  UnknownSet = 0,
  NodeSet = 1<<0,
  SideSet = 1<<1,
  BlockSet = 1<<2,
  MaterialSet = 1<<3,
  DisplacementSet = 1<<4,
  ForceSet = 1<<5,
  PressureSet = 1<<6,
  VelocitySet = 1<<7,
  AccelerationSet = 1<<8,
  TemperatureSet = 1<<9,
  HeatfluxSet = 1<<10,
  InterfaceSet = 1<<11,
  LastSet
};

/*! \struct generic_cubit_bc_data
 *  \brief Generic bc data structure
 */
struct generic_cubit_bc_data {
    PetscErrorCode ierr;
    
    virtual PetscErrorCode fill_data(const vector<char>& bc_data) {
        PetscFunctionBegin;
        SETERRQ(PETSC_COMM_SELF,1,"It makes no sense for the generic bc type");
        PetscFunctionReturn(0);
    }
    
};

/*! \struct displacement_cubit_bc_data
 *  \brief Definition of the displacement bc data structure
 */
struct displacement_cubit_bc_data: public generic_cubit_bc_data {
    struct __attribute__ ((packed)) _data_{
    char name[12]; // 12 characters for "Displacement"
    char pre1; // Always zero
    char pre2; // pre-processing flags for modification of displacement bcs. They should not affect analysis, i.e. safe to ignore; 1: smallest combine, 2: average, 3: largest combine, 4: overwrite or no combination defined (default)
    char flag1; // Flag for X-Translation (0: N/A, 1: specified)
    char flag2; // Flag for Y-Translation (0: N/A, 1: specified)
    char flag3; // Flag for Z-Translation (0: N/A, 1: specified)
    char flag4; // Flag for X-Rotation (0: N/A, 1: specified)
    char flag5; // Flag for Y-Rotation (0: N/A, 1: specified)
    char flag6; // Flag for Z-Rotation (0: N/A, 1: specified)
    double value1; // Value for X-Translation
    double value2; // Value for Y-Translation
    double value3; // Value for Z-Translation
    double value4; // Value for X-Rotation
    double value5; // Value for Y-Rotation
    double value6; // Value for Z-Rotation
    };
    
    _data_ data;

    const Cubit_BC_bitset type;
    displacement_cubit_bc_data(): type(DisplacementSet) {};
    
    virtual PetscErrorCode fill_data(const vector<char>& bc_data) {
        PetscFunctionBegin;
        //Fill data
	if(bc_data.size()!=sizeof(data)) SETERRQ(PETSC_COMM_SELF,1,"data inconsistency");
        memcpy(&data, &bc_data[0], sizeof(data));
        PetscFunctionReturn(0);
    }
 
    /*! \brief Print displacement bc data
     */
    friend ostream& operator<<(ostream& os,const displacement_cubit_bc_data& e);
    
};

/*! \struct force_cubit_bc_data
 *  \brief Definition of the force bc data structure
 */
struct force_cubit_bc_data: public generic_cubit_bc_data {
    struct __attribute__ ((packed)) _data_{
    char name[5]; // 5 characters for "Force"
    char zero[3]; // 3 zeros
    double value1; // Force magnitude
    double value2; // Moment magnitude
    double value3; // X-component of force direction vector
    double value4; // Y-component of force direction vector
    double value5; // Z-component of force direction vector
    double value6; // X-component of moment direction vector
    double value7; // Y-component of moment direction vector
    double value8; // Z-component of moment direction vector
    char zero2; // 0
    };
    
    _data_ data;
    const Cubit_BC_bitset type;
    force_cubit_bc_data(): type(ForceSet) {};

    virtual PetscErrorCode fill_data(const vector<char>& bc_data) {
        PetscFunctionBegin;
        //Fill data
	if(bc_data.size()!=sizeof(data)) SETERRQ(PETSC_COMM_SELF,1,"data inconsistency");
        memcpy(&data, &bc_data[0], sizeof(data));
        PetscFunctionReturn(0);
    }
 
    /*! \brief Print force bc data
    */
    friend ostream& operator<<(ostream& os,const force_cubit_bc_data& e);
    
};

/*! \struct velocity_cubit_bc_data
 *  \brief Definition of the velocity bc data structure
 */
struct velocity_cubit_bc_data: public generic_cubit_bc_data {
    struct __attribute__ ((packed)) _data_{
    char name[8]; // 8 characters for "Velocity"
    char pre1; // Always zero
    char pre2; // pre-processing flags for modification of displacement bcs. They should not affect analysis, i.e. safe to ignore; 1: smallest combine, 2: average, 3: largest combine, 4: overwrite or no combination defined (default)
    char flag1; // Flag for X-Translation (0: N/A, 1: specified)
    char flag2; // Flag for Y-Translation (0: N/A, 1: specified)
    char flag3; // Flag for Z-Translation (0: N/A, 1: specified)
    char flag4; // Flag for X-Rotation (0: N/A, 1: specified)
    char flag5; // Flag for Y-Rotation (0: N/A, 1: specified)
    char flag6; // Flag for Z-Rotation (0: N/A, 1: specified)
    double value1; // Value for X-Translation
    double value2; // Value for Y-Translation
    double value3; // Value for Z-Translation
    double value4; // Value for X-Rotation
    double value5; // Value for Y-Rotation
    double value6; // Value for Z-Rotation
    };
    
    _data_ data;
    const Cubit_BC_bitset type;
    velocity_cubit_bc_data(): type(VelocitySet) {};
   
    virtual PetscErrorCode fill_data(const vector<char>& bc_data) {
        PetscFunctionBegin;
        //Fill data
	if(bc_data.size()!=sizeof(data)) SETERRQ(PETSC_COMM_SELF,1,"data inconsistency");
        memcpy(&data, &bc_data[0], sizeof(data));
        PetscFunctionReturn(0);
    }
 
    /*! \brief Print velocity bc data
    */
    friend ostream& operator<<(ostream& os,const velocity_cubit_bc_data& e);
    
};  

/*! \struct acceleration_cubit_bc_data
 *  \brief Definition of the acceleration bc data structure
 */    
struct acceleration_cubit_bc_data: public generic_cubit_bc_data {
    struct __attribute__ ((packed)) _data_{
    char name[12]; // 12 characters for "Acceleration"
    char pre1; // Always zero
    char pre2; // pre-processing flags for modification of displacement bcs. They should not affect analysis, i.e. safe to ignore; 1: smallest combine, 2: average, 3: largest combine, 4: overwrite or no combination defined (default)
    char flag1; // Flag for X-Translation (0: N/A, 1: specified)
    char flag2; // Flag for Y-Translation (0: N/A, 1: specified)
    char flag3; // Flag for Z-Translation (0: N/A, 1: specified)
    char flag4; // Flag for X-Rotation (0: N/A, 1: specified)
    char flag5; // Flag for Y-Rotation (0: N/A, 1: specified)
    char flag6; // Flag for Z-Rotation (0: N/A, 1: specified)
    double value1; // Value for X-Translation
    double value2; // Value for Y-Translation
    double value3; // Value for Z-Translation
    double value4; // Value for X-Rotation
    double value5; // Value for Y-Rotation
    double value6; // Value for Z-Rotation
    };
    
    _data_ data;
    const Cubit_BC_bitset type;
    acceleration_cubit_bc_data(): type(AccelerationSet) {};

    virtual PetscErrorCode fill_data(const vector<char>& bc_data) {
        PetscFunctionBegin;
        //Fill data
	if(bc_data.size()!=sizeof(data)) SETERRQ(PETSC_COMM_SELF,1,"data inconsistency");
        memcpy(&data, &bc_data[0], sizeof(data));
        PetscFunctionReturn(0);
    }
    
    /*! \brief Print acceleration bc data
    */
    friend ostream& operator<<(ostream& os,const acceleration_cubit_bc_data& e);
    
};

/*! \struct temperature_cubit_bc_data
 *  \brief Definition of the temperature bc data structure
 */
struct temperature_cubit_bc_data: public generic_cubit_bc_data {
    struct __attribute__ ((packed)) _data_{
    char name[11]; // 11 characters for "Temperature"
    char pre1; // This is always zero
    char pre2; // 0: temperature is not applied on thin shells (default); 1: temperature is applied on thin shells
    char flag1; // 0: N/A, 1: temperature value applied (not on thin shells)
    char flag2; // 0: N/A, 1: temperature applied on thin shell middle
    char flag3; // 0: N/A, 1: thin shell temperature gradient specified
    char flag4; // 0: N/A, 1: top thin shell temperature
    char flag5; // 0: N/A, 1: bottom thin shell temperature
    char flag6; // This is always zero
    double value1; // Temperature (default case - no thin shells)
    double value2; // Temperature for middle of thin shells
    double value3; // Temperature gradient for thin shells
    double value4; // Temperature for top of thin shells
    double value5; // Temperature for bottom of thin shells
    double value6; // This is always zero, i.e. ignore
    };
    
    _data_ data;
    const Cubit_BC_bitset type;
    temperature_cubit_bc_data(): type(TemperatureSet) {};

    virtual PetscErrorCode fill_data(const vector<char>& bc_data) {
        PetscFunctionBegin;
        //Fill data
	if(bc_data.size()!=sizeof(data)) SETERRQ(PETSC_COMM_SELF,1,"data inconsistency");
        memcpy(&data, &bc_data[0], sizeof(data));
        PetscFunctionReturn(0);
    }
    
    /*! \brief Print temperature bc data
    */
    friend ostream& operator<<(ostream& os,const temperature_cubit_bc_data& e);
};

/*! \struct pressure_cubit_bc_data
 *  \brief Definition of the pressure bc data structure
 */
struct pressure_cubit_bc_data: public generic_cubit_bc_data {
    struct __attribute__ ((packed)) _data_{
    char name[8]; // 8 characters for "Pressure"
    char flag1; // This is always zero
    char flag2; // 0: Pressure is interpeted as pure pressure 1: pressure is interpreted as total force
    double value1; // Pressure value
    char zero; // This is always zero
    };
    
    _data_ data;
    const Cubit_BC_bitset type;
    pressure_cubit_bc_data(): type(PressureSet) {};
   
    virtual PetscErrorCode fill_data(const vector<char>& bc_data) {
        PetscFunctionBegin;
        //Fill data
	if(bc_data.size()!=sizeof(data)) SETERRQ(PETSC_COMM_SELF,1,"data inconsistency");
        memcpy(&data, &bc_data[0], sizeof(data));
        PetscFunctionReturn(0);
    }
 
    /*! \brief Print pressure bc data
    */
    friend ostream& operator<<(ostream& os,const pressure_cubit_bc_data& e);
    
};

/*! \struct heatflux_cubit_bc_data
 *  \brief Definition of the heat flux bc data structure
 */
struct heatflux_cubit_bc_data: public generic_cubit_bc_data {
    struct __attribute__ ((packed)) _data_{
    char name[8]; // 8 characters for "HeatFlux" (no space)
    char pre1; // This is always zero
    char pre2; // 0: heat flux is not applied on thin shells (default); 1: heat flux is applied on thin shells
    char flag1; // 0: N/A, 1: normal heat flux case (i.e. single value, case without thin shells)
    char flag2; // 0: N/A, 1: Thin shell top heat flux specified
    char flag3; // 0: N/A, 1: Thin shell bottom heat flux specidied
    double value1; // Heat flux value for default case (no thin shells)
    double value2; // Heat flux (thin shell top)
    double value3; // Heat flux (thin shell bottom)
    };
    
    _data_ data;
    const Cubit_BC_bitset type;
    heatflux_cubit_bc_data(): type(HeatfluxSet) {};

    virtual PetscErrorCode fill_data(const vector<char>& bc_data) {
        PetscFunctionBegin;
        //Fill data
	if(bc_data.size()!=sizeof(data)) SETERRQ(PETSC_COMM_SELF,1,"data inconsistency");
        memcpy(&data, &bc_data[0], sizeof(data));
        PetscFunctionReturn(0);
    }
 
    /*! \brief Print heat flux bc data
    */
    friend ostream& operator<<(ostream& os,const heatflux_cubit_bc_data& e);
    
};

    /*! \struct interface_cubit_bc_data
     *  \brief Definition of the interface (cfd_bc) data structure
     */
    struct interface_cubit_bc_data: public generic_cubit_bc_data {
        struct __attribute__ ((packed)) _data_{
            char name[6]; // 6 characters for "cfd_bc"
            char pre1; // This is always zero
            char pre2; // 6
        };
        
        _data_ data;
        const Cubit_BC_bitset type;
        interface_cubit_bc_data(): type(InterfaceSet) {};
        
        virtual PetscErrorCode fill_data(const vector<char>& bc_data) {
            PetscFunctionBegin;
            //Fill data
            if(bc_data.size()!=sizeof(data)) SETERRQ(PETSC_COMM_SELF,1,"data inconsistency");
            memcpy(&data, &bc_data[0], sizeof(data));
            PetscFunctionReturn(0);
        }
        
        /*! \brief Print interface bc data
         */
        friend ostream& operator<<(ostream& os,const interface_cubit_bc_data& e);
        
};
    
/** 
 * \brief this struct keeps basic methods for moab meshset about material and boundary conditions
 */
struct CubitMeshSets {
  EntityHandle meshset;
  Cubit_BC_bitset CubitBCType;
  vector<Tag> tag_handles;
  int *msId;
  char* tag_bc_data;
  int tag_bc_size;
  unsigned int *tag_block_header_data;
  double* tag_block_attributes;
  int tag_block_attributes_size;
  const Cubit_BC_bitset meshsets_mask;
  CubitMeshSets(Interface &moab,const EntityHandle _meshset);
  inline int get_msId() const { return *msId; }
  inline Cubit_BC_bitset get_CubitBCType() const { return CubitBCType; }

  inline unsigned long int get_CubitBCType_ulong() const { return CubitBCType.to_ulong(); }
  inline unsigned long int get_CubitBCType_mask_meshset_types_ulong() const { return (CubitBCType&meshsets_mask).to_ulong(); }
  inline unsigned long int get_CubitBCType_bc_data_types_ulong() const { return (CubitBCType&(~meshsets_mask)).to_ulong(); }

  PetscErrorCode get_Cubit_msId_entities_by_dimension(Interface &moab,const int dimension,Range &entities,const bool recursive = false) const;
  PetscErrorCode get_Cubit_msId_entities_by_dimension(Interface &moab,Range &entities,const bool recursive = false)  const;

  /** 
   *  \brief Function that returns the Cubit_BC_bitset type of the contents of bc_data
   */
  PetscErrorCode get_type_from_bc_data(const vector<char> &bc_data,Cubit_BC_bitset &type) const;

  /** 
   *  \brief Function that returns the Cubit_BC_bitset type of the contents of bc_data
  */
  PetscErrorCode get_type_from_bc_data(Cubit_BC_bitset &type) const;
    
  /**
   * \brief get bc_data vector from MoFEM database
   * 
   * \param bc_data is the in/out vector were bc_data will be stored
   */
  PetscErrorCode get_Cubit_bc_data(vector<char>& bc_data) const;
    
  /**
   * \brief print bc_data int stream given by os
   *
   * f.e. it->print_Cubit_bc_data(cout), i.e. printing to standard output
   * f.e. it->print_Cubit_bc_data(cerr), i.e. printing to standard error output
   */
  PetscErrorCode print_Cubit_bc_data(ostream& os) const;

  template<class _CUBIT_BC_DATA_TYPE_>
  PetscErrorCode get_cubit_bc_data_structure(_CUBIT_BC_DATA_TYPE_& data) const {
    PetscFunctionBegin;
    PetscErrorCode ierr;
    if((CubitBCType&data.type).none()) {
      SETERRQ(PETSC_COMM_SELF,1,"bc_data are not for _CUBIT_BC_DATA_TYPE_ structure");  
    }
    vector<char> bc_data;
    get_Cubit_bc_data(bc_data);
    ierr = data.fill_data(bc_data); CHKERRQ(ierr);
    PetscFunctionReturn(0);
  }
  
  /**
   * \brief get Cubit block attributes
   *
   * \param attributes is the vector where the block attribute data will be stored
   */
  PetscErrorCode get_Cubit_attributes(vector<double> &attributes) const;

  /**
   * \brief print the attributes vector
   *
   * f.e. it->print_Cubit_attributes(cout), i.e. printing to standard output
   * f.e. it->print_Cubit_attributes(cerr), i.e. printing to standard error output
   */
  PetscErrorCode print_Cubit_attributes(ostream& os) const;
    
  friend ostream& operator<<(ostream& os,const CubitMeshSets& e);
    
};
    
/**
 * @relates multi_index_container
 * \brief moabCubitMeshSet_multiIndex
 *
 * \param hashed_unique<
      tag<Meshset_mi_tag>, member<CubitMeshSets,EntityHandle,&CubitMeshSets::meshset> >,
 * \param ordered_non_unique<
      tag<CubitMeshSets_mi_tag>, const_mem_fun<CubitMeshSets,unsigned long int,&CubitMeshSets::get_CubitBCType_ulong> >,
 * \param ordered_non_unique<
      tag<CubitMeshSets_mask_meshset_mi_tag>, const_mem_fun<CubitMeshSets,unsigned long int,&CubitMeshSets::get_CubitBCType_mask_meshset_types_ulong> >,
 * \param ordered_non_unique<
      tag<CubitMeshSets_bc_data_mi_tag>, const_mem_fun<CubitMeshSets,unsigned long int,&CubitMeshSets::get_CubitBCType_bc_data_types_ulong> >,
 *
 * \param    hashed_unique<
      tag<Composite_mi_tag>,       
      composite_key<
	CubitMeshSets, <br>
	  const_mem_fun<CubitMeshSets,int,&CubitMeshSets::get_msId>,
	  const_mem_fun<CubitMeshSets,unsigned long int,&CubitMeshSets::get_CubitBCType_ulong> > >
 *
 */
typedef multi_index_container<
  CubitMeshSets,
  indexed_by<
    hashed_unique<
      tag<Meshset_mi_tag>, member<CubitMeshSets,EntityHandle,&CubitMeshSets::meshset> >,
    ordered_non_unique<
      tag<CubitMeshSets_mi_tag>, const_mem_fun<CubitMeshSets,unsigned long int,&CubitMeshSets::get_CubitBCType_ulong> >,
    ordered_non_unique<
      tag<CubitMeshSets_mask_meshset_mi_tag>, const_mem_fun<CubitMeshSets,unsigned long int,&CubitMeshSets::get_CubitBCType_mask_meshset_types_ulong> >,
    ordered_non_unique<
      tag<CubitMeshSets_bc_data_mi_tag>, const_mem_fun<CubitMeshSets,unsigned long int,&CubitMeshSets::get_CubitBCType_bc_data_types_ulong> >,
    hashed_unique<
      tag<Composite_mi_tag>,       
      composite_key<
	CubitMeshSets,
	  const_mem_fun<CubitMeshSets,int,&CubitMeshSets::get_msId>,
	  const_mem_fun<CubitMeshSets,unsigned long int,&CubitMeshSets::get_CubitBCType_mask_meshset_types_ulong> > >
  > > moabCubitMeshSet_multiIndex;

}

#endif //__COMMON_HPP__
