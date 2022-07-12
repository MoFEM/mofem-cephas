/** \file AdjacencyMultiIndices.hpp
 * \brief Myltindex containes, data structures for mofem adjacencies and other
 * low-level functions
 */

/* MIT License
 *
 * Copyright (c) 2022
 * 
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 * 
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 * 
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */

#ifndef __ADJACENCYMULTIINDICES_HPP__
#define __ADJACENCYMULTIINDICES_HPP__

namespace MoFEM {

/**
 * \brief FieldEntityEntFiniteElementAdjacencyMap of mofem finite element and
 * entities
 *
 */
struct FieldEntityEntFiniteElementAdjacencyMap {
  unsigned int byWhat;                              ///< see options \ref ByWhat
  const boost::shared_ptr<FieldEntity> entFieldPtr; ///< field entity
  const boost::shared_ptr<EntFiniteElement> entFePtr; ///< finite element entity
  FieldEntityEntFiniteElementAdjacencyMap(
      const boost::shared_ptr<FieldEntity> &ent_field_ptr,
      const boost::shared_ptr<EntFiniteElement> &ent_fe_ptr);

  virtual ~FieldEntityEntFiniteElementAdjacencyMap() = default;

  /**
   * \brief get unique iD of finite element entity
   */
  inline UId getFeUniqueId() const { return entFePtr->getLocalUniqueId(); }

  /**
   * \brief get meshset of finite element
   */
  inline EntityHandle getFeMeshset() const { return entFePtr->getMeshset(); }

  /**
   * \brief get finite element handle
   */
  inline EntityHandle getFeHandle() const { return entFePtr->getEnt(); }

  /**
   * \brief get unique iD of entity on field
   */
  inline UId getEntUniqueId() const { return entFieldPtr->getLocalUniqueId(); }

  /**
   * \brief get entity meshset carrying its field
   */
  inline EntityHandle getEntMeshset() const {
    return entFieldPtr->getMeshset();
  }

  /**
   * \brief get entity handle
   */
  inline EntityHandle getEntHandle() const { return entFieldPtr->getEnt(); }

  /**
   * \brief get field iD
   */
  BitFieldId getEntId() const { return entFieldPtr->getId(); }

  /**
   * \brief get finite element iD
   */
  BitFEId getBitFEId() const { return entFePtr->getId(); }

  friend std::ostream &
  operator<<(std::ostream &os,
             const FieldEntityEntFiniteElementAdjacencyMap &e);
};

/**
 * @relates multi_index_container
 * \brief MultiIndex container keeps Adjacencies Element and dof entities
 adjacencies and vice versa.

 */
typedef multi_index_container<
    FieldEntityEntFiniteElementAdjacencyMap,
    indexed_by<
        ordered_unique<
            tag<Composite_Unique_mi_tag>,
            composite_key<
                FieldEntityEntFiniteElementAdjacencyMap,
                const_mem_fun<
                    FieldEntityEntFiniteElementAdjacencyMap, UId,
                    &FieldEntityEntFiniteElementAdjacencyMap::getEntUniqueId>,
                const_mem_fun<
                    FieldEntityEntFiniteElementAdjacencyMap, UId,
                    &FieldEntityEntFiniteElementAdjacencyMap::getFeUniqueId>>>,
        ordered_non_unique<
            tag<Unique_mi_tag>,
            const_mem_fun<
                FieldEntityEntFiniteElementAdjacencyMap, UId,
                &FieldEntityEntFiniteElementAdjacencyMap::getEntUniqueId>>,
        ordered_non_unique<
            tag<FE_Unique_mi_tag>,
            const_mem_fun<
                FieldEntityEntFiniteElementAdjacencyMap, UId,
                &FieldEntityEntFiniteElementAdjacencyMap::getFeUniqueId>>,
        ordered_non_unique<
            tag<FEEnt_mi_tag>,
            const_mem_fun<
                FieldEntityEntFiniteElementAdjacencyMap, EntityHandle,
                &FieldEntityEntFiniteElementAdjacencyMap::getFeHandle>>,
        ordered_non_unique<
            tag<Ent_mi_tag>,
            const_mem_fun<
                FieldEntityEntFiniteElementAdjacencyMap, EntityHandle,
                &FieldEntityEntFiniteElementAdjacencyMap::getEntHandle>>>>
    FieldEntityEntFiniteElementAdjacencyMap_multiIndex;

struct FieldEntityEntFiniteElementAdjacencyMap_change_ByWhat {
  int bY;
  FieldEntityEntFiniteElementAdjacencyMap_change_ByWhat(const int by)
      : bY(by) {}
  void operator()(FieldEntityEntFiniteElementAdjacencyMap &e) {
    e.byWhat |= bY;
  }
};

} // namespace MoFEM

#endif // __ADJACENCYMULTIINDICES_HPP__
