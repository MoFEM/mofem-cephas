/** \file AdjacencyMultiIndices.hpp
 * \brief Myltindex containes, data structures for mofem adjacencies and other low-level functions
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
  * \brief MoFEMEntityEntFiniteElementAdjacencyMap of mofem finite element and entities
  *
  */
struct MoFEMEntityEntFiniteElementAdjacencyMap {
  unsigned int by_other;
  const MoFEMEntity *MoFEMEntity_ptr; ///< field entity
  const EntFiniteElement *EntFiniteElement_ptr; ///< finite element entity
  MoFEMEntityEntFiniteElementAdjacencyMap(const MoFEMEntity *_MoFEMEntity_ptr,const EntFiniteElement *_EntFiniteElement_ptr);
  inline GlobalUId get_MoFEMFiniteElement_unique_id() const { return EntFiniteElement_ptr->get_global_unique_id(); }
  inline EntityHandle get_MoFEMFiniteElement_meshset() const { return EntFiniteElement_ptr->get_meshset(); }
  inline EntityHandle get_MoFEMFiniteElement_entity_handle() const { return EntFiniteElement_ptr->get_ent(); }
  inline GlobalUId get_ent_unique_id() const { return MoFEMEntity_ptr->get_global_unique_id(); };
  inline EntityHandle get_ent_meshset() const { return MoFEMEntity_ptr->get_meshset(); };
  inline EntityHandle get_ent_entity_handle() const { return MoFEMEntity_ptr->get_ent(); };
  BitFieldId get_ent_id() const { return MoFEMEntity_ptr->get_id(); }
  BitFEId get_BitFEId() const { return EntFiniteElement_ptr->get_id(); }
  friend ostream& operator<<(ostream& os,const MoFEMEntityEntFiniteElementAdjacencyMap &e);
};

/**
 * @relates multi_index_container
 * \brief MultiIndex container keeps Adjacencies Element and dof entities adjacencies and vice veras.

 */
typedef multi_index_container<
  MoFEMEntityEntFiniteElementAdjacencyMap,
  indexed_by<
    ordered_unique<
      tag<Composite_Unique_mi_tag>,
      composite_key<
	MoFEMEntityEntFiniteElementAdjacencyMap,
	const_mem_fun<MoFEMEntityEntFiniteElementAdjacencyMap,GlobalUId,&MoFEMEntityEntFiniteElementAdjacencyMap::get_ent_unique_id>,
	const_mem_fun<MoFEMEntityEntFiniteElementAdjacencyMap,GlobalUId,&MoFEMEntityEntFiniteElementAdjacencyMap::get_MoFEMFiniteElement_unique_id> > >,
    ordered_non_unique<
      tag<Unique_mi_tag>, const_mem_fun<MoFEMEntityEntFiniteElementAdjacencyMap,GlobalUId,&MoFEMEntityEntFiniteElementAdjacencyMap::get_ent_unique_id> >,
    ordered_non_unique<
      tag<FEEnt_mi_tag>, const_mem_fun<MoFEMEntityEntFiniteElementAdjacencyMap,EntityHandle,&MoFEMEntityEntFiniteElementAdjacencyMap::get_MoFEMFiniteElement_entity_handle> >,
    ordered_non_unique<
      tag<Ent_mi_tag>, const_mem_fun<MoFEMEntityEntFiniteElementAdjacencyMap,EntityHandle,&MoFEMEntityEntFiniteElementAdjacencyMap::get_ent_entity_handle> >
  > > MoFEMEntityEntFiniteElementAdjacencyMap_multiIndex;

  struct MoFEMEntityEntFiniteElementAdjacencyMap_change_ByWhat {
    ByWhat by;
    MoFEMEntityEntFiniteElementAdjacencyMap_change_ByWhat(const ByWhat _by): by(_by) {}
    void operator()(MoFEMEntityEntFiniteElementAdjacencyMap &e) { e.by_other |= by; }
  };

}

#endif // __ADJACENCYMULTIINDICES_HPP__
