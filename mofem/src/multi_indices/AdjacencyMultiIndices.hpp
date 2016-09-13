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
  const boost::shared_ptr<MoFEMEntity> mofemEntPtr; ///< field entity
  const boost::shared_ptr<EntFiniteElement> entFePtr; ///< finite element entity
  MoFEMEntityEntFiniteElementAdjacencyMap(
    const boost::shared_ptr<MoFEMEntity> mofem_ent_ptr,
    const boost::shared_ptr<EntFiniteElement> ent_fe_ptr
  );
  inline GlobalUId get_MoFEMFiniteElement_unique_id() const { return entFePtr->getGlobalUniqueId(); }
  inline EntityHandle get_MoFEMFiniteElement_meshset() const { return entFePtr->getMeshset(); }
  inline EntityHandle get_MoFEMFiniteElement_entity_handle() const { return entFePtr->getEnt(); }
  inline GlobalUId get_ent_unique_id() const { return mofemEntPtr->getGlobalUniqueId(); };
  inline EntityHandle get_ent_meshset() const { return mofemEntPtr->getMeshset(); };
  inline EntityHandle get_ent_entity_handle() const { return mofemEntPtr->getEnt(); };
  BitFieldId getEntId() const { return mofemEntPtr->getId(); }
  BitFEId get_BitFEId() const { return entFePtr->getId(); }
  friend std::ostream& operator<<(std::ostream& os,const MoFEMEntityEntFiniteElementAdjacencyMap &e);
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
    int bY;
    MoFEMEntityEntFiniteElementAdjacencyMap_change_ByWhat(const int by): bY(by) {}
    void operator()(MoFEMEntityEntFiniteElementAdjacencyMap &e) {
      e.by_other |= bY;
    }
  };

}

#endif // __ADJACENCYMULTIINDICES_HPP__
