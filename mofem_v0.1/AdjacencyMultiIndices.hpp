/** \file AdjacencyMultiIndices.hpp
 * \brief Myltindex containes, data structures for mofem adjacencies and other low-level functions 
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

#ifndef __ADJACENCYMULTIINDICES_HPP__
#define __ADJACENCYMULTIINDICES_HPP__


namespace MoFEM {

/** 
 * \brief struct keeps data about selected prism adjacencies, and potentially other entities
 */
struct BasicMoFEMEntityAdjacenctMap: public BasicMoFEMEntity {
  BasicMoFEMEntity Adj; ///< adjacent entity to this BasicMoFEMEntityAdjacenctMap
  BasicMoFEMEntityAdjacenctMap(const EntityHandle _ent,const EntityHandle adj):
    BasicMoFEMEntity(_ent), Adj(adj) {};
  inline EntityHandle get_adj() const { return Adj.ent; };
  inline EntityType get_adj_type() const { return Adj.get_ent_type(); };
};

/** 
 * @relates multi_index_container
 * \brief MultiIndex container keeps BasicMoFEMEntityAdjacenctMap
 *
 * \param   hashed_non_unique<
      tag<MoABEnt_mi_tag>, 
      member<BasicMoFEMEntityAdjacenctMap::BasicMoFEMEntity,EntityHandle,&BasicMoFEMEntityAdjacenctMap::ent> >,
 * \param    hashed_non_unique<
      tag<MoABEnt_mi_tag2>, 
      const_mem_fun<BasicMoFEMEntityAdjacenctMap,EntityHandle,&BasicMoFEMEntityAdjacenctMap::get_adj> >,
 * \param    ordered_non_unique<
      tag<EntType_mi_tag>, 
      const_mem_fun<BasicMoFEMEntityAdjacenctMap,EntityType,&BasicMoFEMEntityAdjacenctMap::get_adj_type> >,
 * \param    hashed_unique< tag<Composite_mi_tag>, <br>
      composite_key< <br>
	BasicMoFEMEntityAdjacenctMap,
      	member<BasicMoFEMEntityAdjacenctMap::BasicMoFEMEntity,EntityHandle,&BasicMoFEMEntityAdjacenctMap::ent>,
	const_mem_fun<BasicMoFEMEntityAdjacenctMap,EntityHandle,&BasicMoFEMEntityAdjacenctMap::get_adj> > >
 *
 */
typedef multi_index_container<
  BasicMoFEMEntityAdjacenctMap,
  indexed_by<
    hashed_non_unique<
      tag<MoABEnt_mi_tag>, 
      member<BasicMoFEMEntityAdjacenctMap::BasicMoFEMEntity,EntityHandle,&BasicMoFEMEntityAdjacenctMap::ent> >,
    hashed_non_unique<
      tag<MoABEnt_mi_tag2>, 
      const_mem_fun<BasicMoFEMEntityAdjacenctMap,EntityHandle,&BasicMoFEMEntityAdjacenctMap::get_adj> >,
    ordered_non_unique<
      tag<EntType_mi_tag>, 
      const_mem_fun<BasicMoFEMEntityAdjacenctMap,EntityType,&BasicMoFEMEntityAdjacenctMap::get_adj_type> >,
    hashed_unique<
      tag<Composite_mi_tag>, 
      composite_key<
	BasicMoFEMEntityAdjacenctMap,
      	member<BasicMoFEMEntityAdjacenctMap::BasicMoFEMEntity,EntityHandle,&BasicMoFEMEntityAdjacenctMap::ent>,
	const_mem_fun<BasicMoFEMEntityAdjacenctMap,EntityHandle,&BasicMoFEMEntityAdjacenctMap::get_adj> > >
  > > BasicMoFEMEntityAdjacenctMap_multiIndex;

/**
  * \brief MoFEMEntityEntMoFEMFiniteElementAdjacencyMap of mofem finite element and entities
  *
  */
struct MoFEMEntityEntMoFEMFiniteElementAdjacencyMap {
  unsigned int by_other;
  const MoFEMEntity *MoFEMEntity_ptr; ///< field entity
  const EntMoFEMFiniteElement *EntMoFEMFiniteElement_ptr; ///< finite element entity
  MoFEMEntityEntMoFEMFiniteElementAdjacencyMap(const MoFEMEntity *_MoFEMEntity_ptr,const EntMoFEMFiniteElement *_EntMoFEMFiniteElement_ptr);
  inline UId get_MoFEMFiniteElement_unique_id() const { return EntMoFEMFiniteElement_ptr->get_unique_id(); }
  inline EntityHandle get_MoFEMFiniteElement_meshset() const { return EntMoFEMFiniteElement_ptr->get_meshset(); }
  inline EntityHandle get_MoFEMFiniteElement_entity_handle() const { return EntMoFEMFiniteElement_ptr->get_ent(); }
  inline UId get_ent_unique_id() const { return MoFEMEntity_ptr->get_unique_id(); };
  inline EntityHandle get_ent_meshset() const { return MoFEMEntity_ptr->get_meshset(); };
  inline EntityHandle get_ent_entity_handle() const { return MoFEMEntity_ptr->get_ent(); };
  BitFieldId get_ent_id() const { return MoFEMEntity_ptr->get_id(); }
  BitFEId get_BitFEId() const { return EntMoFEMFiniteElement_ptr->get_id(); }
  friend ostream& operator<<(ostream& os,const MoFEMEntityEntMoFEMFiniteElementAdjacencyMap &e);
};

/** 
 * @relates multi_index_container
 * \brief MultiIndex container keeps Adjacencies Element and dof entities adjacencies and vice veras.
 *
 * \param    hashed_unique< tag<Composite_unique_mi_tag>, <br>
      composite_key<
	MoFEMEntityEntMoFEMFiniteElementAdjacencyMap, <br>
	const_mem_fun<MoFEMEntityEntMoFEMFiniteElementAdjacencyMap,UId,&MoFEMEntityEntMoFEMFiniteElementAdjacencyMap::get_ent_unique_id>, <br>
	const_mem_fun<MoFEMEntityEntMoFEMFiniteElementAdjacencyMap,UId,&MoFEMEntityEntMoFEMFiniteElementAdjacencyMap::get_MoFEMFiniteElement_unique_id> > >, <br>
 * \param    ordered_non_unique<
      tag<Unique_mi_tag>, const_mem_fun<MoFEMEntityEntMoFEMFiniteElementAdjacencyMap,UId,&MoFEMEntityEntMoFEMFiniteElementAdjacencyMap::get_ent_unique_id> >
 *
 */
typedef multi_index_container<
  MoFEMEntityEntMoFEMFiniteElementAdjacencyMap,
  indexed_by<
    ordered_non_unique<
      tag<Composite_unique_mi_tag>,       
      composite_key<
	MoFEMEntityEntMoFEMFiniteElementAdjacencyMap,
	const_mem_fun<MoFEMEntityEntMoFEMFiniteElementAdjacencyMap,UId,&MoFEMEntityEntMoFEMFiniteElementAdjacencyMap::get_ent_unique_id>,
	const_mem_fun<MoFEMEntityEntMoFEMFiniteElementAdjacencyMap,UId,&MoFEMEntityEntMoFEMFiniteElementAdjacencyMap::get_MoFEMFiniteElement_unique_id> > >,
    ordered_non_unique<
      tag<Unique_mi_tag>, const_mem_fun<MoFEMEntityEntMoFEMFiniteElementAdjacencyMap,UId,&MoFEMEntityEntMoFEMFiniteElementAdjacencyMap::get_ent_unique_id> >,
    ordered_non_unique<
      tag<MoABFEEnt_mi_tag>, const_mem_fun<MoFEMEntityEntMoFEMFiniteElementAdjacencyMap,EntityHandle,&MoFEMEntityEntMoFEMFiniteElementAdjacencyMap::get_MoFEMFiniteElement_entity_handle> >,
    ordered_non_unique<
      tag<MoABEnt_mi_tag>, const_mem_fun<MoFEMEntityEntMoFEMFiniteElementAdjacencyMap,EntityHandle,&MoFEMEntityEntMoFEMFiniteElementAdjacencyMap::get_ent_entity_handle> >
  > > MoFEMEntityEntMoFEMFiniteElementAdjacencyMap_multiIndex;

}

#endif // __ADJACENCYMULTIINDICES_HPP__
