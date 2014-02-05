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

#ifndef __FIELDMULTIINDICES_HPP__
#define __FIELDMULTIINDICES_HPP__

namespace MoFEM {

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


}

#endif // __FIELDMULTIINDICES_HPP__
