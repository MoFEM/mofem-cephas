/** \file RefElementMultiIndices.hpp
 * \brief Multi-index contains, data structures for mofem finite elements and
 * other low-level functions
 */

/* MoFEM is free software: you can redistribute it and/or modify it under
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

#ifndef __REF_ELEMENT_MULTIINDICES_HPP__
#define __REF_ELEMENT_MULTIINDICES_HPP__

namespace MoFEM {

/**
 * \brief keeps data about abstract refined finite element
 * \ingroup fe_multi_indices
 */
struct RefElement : public interface_RefEntity<RefEntity> {

  typedef interface_RefEntity<RefEntity> interface_type_RefEntity;

  static BitRefEdges DummyBitRefEdges;

  SideNumber_multiIndex side_number_table;
  RefElement(const boost::shared_ptr<RefEntity> &ref_ents_ptr);
  virtual ~RefElement() = default;

  virtual const BitRefEdges &getBitRefEdges() const { return DummyBitRefEdges; }

  virtual int getBitRefEdgesUlong() const { return 0; }

  SideNumber_multiIndex &getSideNumberTable() const {
    return const_cast<SideNumber_multiIndex &>(side_number_table);
  }

  static const boost::shared_ptr<SideNumber> nullSideNumber;

  virtual const boost::shared_ptr<SideNumber> &
  getSideNumberPtr(const EntityHandle ent) const {
    NOT_USED(ent);
    return nullSideNumber;
  };

  /**
   * \brief Get pointer to RefEntity
   */
  inline boost::shared_ptr<RefEntity> &getRefEntityPtr() const {
    return this->sPtr;
  }

  friend std::ostream &operator<<(std::ostream &os, const RefElement &e);
};

/**
 * \brief keeps data about abstract MESHSET finite element
 * \ingroup fe_multi_indices
 */
struct RefElement_MESHSET : public RefElement {
  RefElement_MESHSET(const boost::shared_ptr<RefEntity> &ref_ents_ptr);
  virtual ~RefElement_MESHSET() = default;
  const boost::shared_ptr<SideNumber> &
  getSideNumberPtr(const EntityHandle ent) const;
};
/**
 * \brief keeps data about abstract PRISM finite element
 * \ingroup fe_multi_indices
 */
struct RefElement_PRISM : public RefElement {
  BitRefEdges *tag_BitRefEdges;
  RefElement_PRISM(const boost::shared_ptr<RefEntity> &ref_ents_ptr);
  virtual ~RefElement_PRISM() = default;

  const boost::shared_ptr<SideNumber> &
  getSideNumberPtr(const EntityHandle ent) const;
  const BitRefEdges &getBitRefEdges() const { return *tag_BitRefEdges; }
  int getBitRefEdgesUlong() const { return getBitRefEdges().to_ulong(); }
};

/**
 * \brief keeps data about abstract TET finite element
 * \ingroup fe_multi_indices
 */
struct RefElementVolume : public RefElement {
  BitRefEdges *tag_BitRefEdges;
  RefElementVolume(const boost::shared_ptr<RefEntity> &ref_ents_ptr);
  virtual ~RefElementVolume() = default;

  const boost::shared_ptr<SideNumber> &
  getSideNumberPtr(const EntityHandle ent) const;
  SideNumber_multiIndex &getSideNumberTable() const {
    return const_cast<SideNumber_multiIndex &>(side_number_table);
  };
  const BitRefEdges &getBitRefEdges() const { return *tag_BitRefEdges; }
  int getBitRefEdgesUlong() const { return getBitRefEdges().to_ulong(); }
};

/**
 * \brief keeps data about abstract TRI finite element
 * \ingroup fe_multi_indices
 */
struct RefElementFace : public RefElement {
  RefElementFace(const boost::shared_ptr<RefEntity> &ref_ents_ptr);
  virtual ~RefElementFace() = default;
  const boost::shared_ptr<SideNumber> &
  getSideNumberPtr(const EntityHandle ent) const;
  friend std::ostream &operator<<(std::ostream &os, const RefElementFace &e);
};

/**
 * \brief keeps data about abstract EDGE finite element
 * \ingroup fe_multi_indices
 */
struct RefElement_EDGE : public RefElement {
  RefElement_EDGE(const boost::shared_ptr<RefEntity> &ref_ents_ptr);
  virtual ~RefElement_EDGE() = default;
  const boost::shared_ptr<SideNumber> &
  getSideNumberPtr(const EntityHandle ent) const;
  friend std::ostream &operator<<(std::ostream &os, const RefElement_EDGE &e);
};

/**
 * \brief keeps data about abstract VERTEX finite element
 * \ingroup fe_multi_indices
 */
struct RefElement_VERTEX : public RefElement {
  RefElement_VERTEX(const boost::shared_ptr<RefEntity> &ref_ents_ptr);
  virtual ~RefElement_VERTEX() = default;
  const boost::shared_ptr<SideNumber> &
  getSideNumberPtr(const EntityHandle ent) const;
  friend std::ostream &operator<<(std::ostream &os, const RefElement_VERTEX &e);
};

/**
 * \brief intrface to RefElement
 * \ingroup fe_multi_indices
 */
template <typename T> struct interface_RefElement : interface_RefEntity<T> {

  typedef interface_RefEntity<T> interface_type_RefEntity;
  typedef interface_RefElement<T> interface_type_RefElement;

  interface_RefElement(const boost::shared_ptr<T> &sptr)
      : interface_RefEntity<T>(sptr) {}
  virtual ~interface_RefElement() = default;

  inline int getBitRefEdgesUlong() const {
    return this->sPtr->getBitRefEdgesUlong();
  }

  inline SideNumber_multiIndex &getSideNumberTable() const {
    return this->sPtr->getSideNumberTable();
  }

  inline const boost::shared_ptr<SideNumber> &
  getSideNumberPtr(const EntityHandle ent) const {
    return this->sPtr->getSideNumberPtr(ent);
  }

  inline boost::shared_ptr<RefEntity> &getRefEntityPtr() const {
    return this->sPtr->getRefEntityPtr();
  }

  inline boost::shared_ptr<T> &getRefElement() const { return this->sPtr; }
};

/**
 * \typedef RefElement_multiIndex
 * type multiIndex container for RefElement
 * \ingroup fe_multi_indices
 *
 * \param hashed_unique Ent_mi_tag
 * \param ordered_non_unique Meshset_mi_tag
 * \param ordered_non_unique Ent_Ent_mi_tag
 * \param ordered_non_unique Composite_ParentEnt_And_BitsOfRefinedEdges_mi_tag
 */
typedef multi_index_container<
    boost::shared_ptr<RefElement>,
    // ptrWrapperRefElement,
    indexed_by<ordered_unique<
        tag<Ent_mi_tag>, const_mem_fun<RefElement::interface_type_RefEntity,
                                       EntityHandle, &RefElement::getEnt>>>>
    RefElement_multiIndex;



} // namespace MoFEM

#endif // __REF_ELEMENT_MULTIINDICES_HPP__

/**
 * \defgroup fe_multi_indices Finite elements structures and multi-indices
 * \ingroup mofem
 **/