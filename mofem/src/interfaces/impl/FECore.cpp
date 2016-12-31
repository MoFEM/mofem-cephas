/** \file DeleteCore.cpp
 * \brief Core interface methods for managing deletions and insertion dofs
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

#include <version.h>
#include <Includes.hpp>
#include <definitions.h>
#include <Common.hpp>

#include <h1_hdiv_hcurl_l2.h>

#include <UnknownInterface.hpp>

#include <MaterialBlocks.hpp>
#include <BCData.hpp>
#include <TagMultiIndices.hpp>
#include <CoordSysMultiIndices.hpp>
#include <FieldMultiIndices.hpp>
#include <EntsMultiIndices.hpp>
#include <DofsMultiIndices.hpp>
#include <FEMultiIndices.hpp>
#include <ProblemsMultiIndices.hpp>
#include <AdjacencyMultiIndices.hpp>
#include <BCMultiIndices.hpp>
#include <CoreDataStructures.hpp>
#include <SeriesMultiIndices.hpp>

#include <LoopMethods.hpp>
#include <Interface.hpp>
#include <MeshRefinement.hpp>
#include <PrismInterface.hpp>
#include <SeriesRecorder.hpp>
#include <Core.hpp>

namespace MoFEM {

  bool Core::check_finite_element(const std::string &name) const {
    typedef FiniteElement_multiIndex::index<FiniteElement_name_mi_tag>::type FeSetByName;
    const FeSetByName &set = finiteElements.get<FiniteElement_name_mi_tag>();
    FeSetByName::iterator miit = set.find(name);
    if(miit==set.end()) return false;
    return true;
  }

  PetscErrorCode Core::add_finite_element(const std::string &fe_name,enum MoFEMTypes bh) {
    PetscFunctionBegin;
    *buildMoFEM &= 1<<0;
    typedef FiniteElement_multiIndex::index<FiniteElement_name_mi_tag>::type FiniteElements_by_name;
    FiniteElements_by_name &finite_element_name_set = finiteElements.get<FiniteElement_name_mi_tag>();
    FiniteElements_by_name::iterator it_fe = finite_element_name_set.find(fe_name);
    if(bh == MF_EXCL) {
      if(it_fe!=finite_element_name_set.end()) {
        SETERRQ1(comm,MOFEM_NOT_FOUND,"this < %s > is there",fe_name.c_str());
      }
    } else {
      if(it_fe!=finite_element_name_set.end()) PetscFunctionReturn(0);
    }
    EntityHandle meshset;
    rval = moab.create_meshset(MESHSET_SET|MESHSET_TRACK_OWNER,meshset); CHKERRQ_MOAB(rval);
    //id
    BitFEId id = getFEShift();
    rval = moab.tag_set_data(th_FEId,&meshset,1,&id); CHKERRQ_MOAB(rval);
    //id name
    void const* tag_data[] = { fe_name.c_str() };
    int tag_sizes[1]; tag_sizes[0] = fe_name.size();
    rval = moab.tag_set_by_ptr(th_FEName,&meshset,1,tag_data,tag_sizes); CHKERRQ_MOAB(rval);
    //add FiniteElement
    std::pair<FiniteElement_multiIndex::iterator,bool> p = finiteElements.insert(
      boost::shared_ptr<FiniteElement>(new FiniteElement(moab,meshset))
    );
    if(!p.second) SETERRQ(comm,MOFEM_OPERATION_UNSUCCESSFUL,"FiniteElement not inserted");
    if(verbose>0) {
      std::ostringstream ss;
      ss << "add finite element: " << fe_name << std::endl;
      PetscPrintf(comm,ss.str().c_str());
      //list_finiteElements();
    }
    PetscFunctionReturn(0);
  }

  PetscErrorCode Core::modify_finite_element_adjacency_table(const std::string &fe_name,const EntityType type,ElementAdjacencyFunct function) {
    PetscFunctionBegin;
    *buildMoFEM &= 1<<0;
    typedef FiniteElement_multiIndex::index<FiniteElement_name_mi_tag>::type FiniteElements_by_name;
    FiniteElements_by_name &finite_element_name_set = finiteElements.get<FiniteElement_name_mi_tag>();
    FiniteElements_by_name::iterator it_fe = finite_element_name_set.find(fe_name);
    if(it_fe==finite_element_name_set.end()) {
      SETERRQ(comm,MOFEM_NOT_FOUND,"this FiniteElement is there");
    }
    boost::shared_ptr<FiniteElement> fe;
    fe = *it_fe;
    fe->element_adjacency_table[type] = function;
    PetscFunctionReturn(0);
  }

  PetscErrorCode Core::modify_finite_element_add_field_data(const std::string &fe_name,const std::string &name_data) {
    PetscFunctionBegin;
    *buildMoFEM &= 1<<0;
    typedef FiniteElement_multiIndex::index<FiniteElement_name_mi_tag>::type FiniteElements_by_name;
    FiniteElements_by_name &finite_element_name_set = finiteElements.get<FiniteElement_name_mi_tag>();
    FiniteElements_by_name::iterator it_fe = finite_element_name_set.find(fe_name);
    if(it_fe==finite_element_name_set.end()) SETERRQ(comm,MOFEM_NOT_FOUND,"this FiniteElement is there");
    try {
      bool success = finite_element_name_set.modify(it_fe,MoFEMFiniteElement_change_bit_add(get_BitFieldId(name_data)));
      if(!success) SETERRQ(comm,MOFEM_OPERATION_UNSUCCESSFUL,"modification unsuccessful");
    } catch (MoFEMException const &e) {
      SETERRQ(comm,e.errorCode,e.errorMessage);
    }
    PetscFunctionReturn(0);
  }

  PetscErrorCode Core::modify_finite_element_add_field_row(const std::string &fe_name,const std::string &name_row) {
    PetscFunctionBegin;
    *buildMoFEM &= 1<<0;
    typedef FiniteElement_multiIndex::index<FiniteElement_name_mi_tag>::type FiniteElements_by_name;
    FiniteElements_by_name &finite_element_name_set = finiteElements.get<FiniteElement_name_mi_tag>();
    FiniteElements_by_name::iterator it_fe = finite_element_name_set.find(fe_name);
    if(it_fe==finite_element_name_set.end()) SETERRQ1(comm,MOFEM_NOT_FOUND,"this < %s > is not there",fe_name.c_str());
    try {
      bool success = finite_element_name_set.modify(it_fe,MoFEMFiniteElement_row_change_bit_add(get_BitFieldId(name_row)));
      if(!success) SETERRQ(comm,MOFEM_OPERATION_UNSUCCESSFUL,"modification unsuccessful");
    } catch (MoFEMException const &e) {
      SETERRQ(comm,e.errorCode,e.errorMessage);
    }
    PetscFunctionReturn(0);
  }

  PetscErrorCode Core::modify_finite_element_add_field_col(const std::string &fe_name,const std::string &name_col) {
    PetscFunctionBegin;
    *buildMoFEM &= 1<<0;
    typedef FiniteElement_multiIndex::index<FiniteElement_name_mi_tag>::type FiniteElements_by_name;
    FiniteElements_by_name &finite_element_name_set = finiteElements.get<FiniteElement_name_mi_tag>();
    FiniteElements_by_name::iterator it_fe = finite_element_name_set.find(fe_name);
    if(it_fe==finite_element_name_set.end()) SETERRQ(comm,MOFEM_OPERATION_UNSUCCESSFUL,"this FiniteElement is there");
    try {
      bool success = finite_element_name_set.modify(it_fe,MoFEMFiniteElement_col_change_bit_add(get_BitFieldId(name_col)));
      if(!success) SETERRQ(comm,MOFEM_OPERATION_UNSUCCESSFUL,"modification unsuccessful");
    } catch (MoFEMException const &e) {
      SETERRQ(comm,e.errorCode,e.errorMessage);
    }
    PetscFunctionReturn(0);
  }

  PetscErrorCode Core::modify_finite_element_off_field_data(const std::string &fe_name,const std::string &name_data) {
    PetscFunctionBegin;
    *buildMoFEM &= 1<<0;
    typedef FiniteElement_multiIndex::index<FiniteElement_name_mi_tag>::type FiniteElements_by_name;
    FiniteElements_by_name &finite_element_name_set = finiteElements.get<FiniteElement_name_mi_tag>();
    FiniteElements_by_name::iterator it_fe = finite_element_name_set.find(fe_name);
    if(it_fe==finite_element_name_set.end()) SETERRQ(comm,MOFEM_NOT_FOUND,"this FiniteElement is there");
    try {
      bool success = finite_element_name_set.modify(it_fe,MoFEMFiniteElement_change_bit_off(get_BitFieldId(name_data)));
      if(!success) SETERRQ(comm,MOFEM_OPERATION_UNSUCCESSFUL,"modification unsuccessful");
    } catch (MoFEMException const &e) {
      SETERRQ(comm,e.errorCode,e.errorMessage);
    }
    PetscFunctionReturn(0);
  }

  PetscErrorCode Core::modify_finite_element_off_field_row(const std::string &fe_name,const std::string &name_row) {
    PetscFunctionBegin;
    *buildMoFEM &= 1<<0;
    typedef FiniteElement_multiIndex::index<FiniteElement_name_mi_tag>::type FiniteElements_by_name;
    FiniteElements_by_name &finite_element_name_set = finiteElements.get<FiniteElement_name_mi_tag>();
    FiniteElements_by_name::iterator it_fe = finite_element_name_set.find(fe_name);
    if(it_fe==finite_element_name_set.end()) SETERRQ1(comm,MOFEM_NOT_FOUND,"this < %s > is not there",fe_name.c_str());
    try {
      bool success = finite_element_name_set.modify(it_fe,MoFEMFiniteElement_row_change_bit_off(get_BitFieldId(name_row)));
      if(!success) SETERRQ(comm,MOFEM_OPERATION_UNSUCCESSFUL,"modification unsuccessful");
    } catch (MoFEMException const &e) {
      SETERRQ(comm,e.errorCode,e.errorMessage);
    }
    PetscFunctionReturn(0);
  }

  PetscErrorCode Core::modify_finite_element_off_field_col(const std::string &fe_name,const std::string &name_col) {
    PetscFunctionBegin;
    *buildMoFEM &= 1<<0;
    typedef FiniteElement_multiIndex::index<FiniteElement_name_mi_tag>::type FiniteElements_by_name;
    FiniteElements_by_name &finite_element_name_set = finiteElements.get<FiniteElement_name_mi_tag>();
    FiniteElements_by_name::iterator it_fe = finite_element_name_set.find(fe_name);
    if(it_fe==finite_element_name_set.end()) SETERRQ(comm,MOFEM_NOT_FOUND,"this FiniteElement is there");
    try {
      bool success = finite_element_name_set.modify(it_fe,MoFEMFiniteElement_col_change_bit_off(get_BitFieldId(name_col)));
      if(!success) SETERRQ(comm,MOFEM_OPERATION_UNSUCCESSFUL,"modification unsuccessful");
    } catch (MoFEMException const &e) {
      SETERRQ(comm,e.errorCode,e.errorMessage);
    }
    PetscFunctionReturn(0);
  }

  BitFEId Core::getBitFEId(const std::string& name) const {
    typedef FiniteElement_multiIndex::index<FiniteElement_name_mi_tag>::type FiniteElements_by_name;
    const FiniteElements_by_name& set = finiteElements.get<FiniteElement_name_mi_tag>();
    FiniteElements_by_name::iterator miit = set.find(name);
    if(miit==set.end()) THROW_MESSAGE(("finite element < "+name+" > not found (top tip: check spelling)").c_str());
    return (*miit)->getId();
  }

  std::string Core::getBitFEId_name(const BitFEId id) const {
    typedef FiniteElement_multiIndex::index<BitFEId_mi_tag>::type finiteElements_by_id;
    const finiteElements_by_id& set = finiteElements.get<BitFEId_mi_tag>();
    finiteElements_by_id::iterator miit = set.find(id);
    assert(miit!=set.end());
    return (*miit)->getName();
  }

  EntityHandle Core::get_finite_element_meshset(const BitFEId id) const {
    typedef FiniteElement_multiIndex::index<BitFEId_mi_tag>::type finiteElements_by_id;
    const finiteElements_by_id& set = finiteElements.get<BitFEId_mi_tag>();
    finiteElements_by_id::iterator miit = set.find(id);
    if(miit==set.end()) THROW_MESSAGE("finite element not found");
    return (*miit)->meshset;
  }

  EntityHandle Core::get_finite_element_meshset(const std::string& name) const {
    return get_finite_element_meshset(getBitFEId(name));
  }

  PetscErrorCode Core::get_finite_element_entities_by_dimension(const std::string name,int dim,Range &ents) const {
    MoABErrorCode rval;
    PetscFunctionBegin;
    try {
      EntityHandle meshset = get_finite_element_meshset(name);
      rval = moab.get_entities_by_dimension(meshset,dim,ents,true); CHKERRQ_MOAB(rval);
    } catch (MoFEMException const &e) {
      SETERRQ(PETSC_COMM_SELF,e.errorCode,e.errorMessage);
    }
    PetscFunctionReturn(0);
  }

  PetscErrorCode Core::get_finite_element_entities_by_type(const std::string name,EntityType type,Range &ents) const {
    MoABErrorCode rval;
    PetscFunctionBegin;
    try {
      EntityHandle meshset = get_finite_element_meshset(name);
      rval = moab.get_entities_by_type(meshset,type,ents,true); CHKERRQ_MOAB(rval);
    } catch (MoFEMException const &e) {
      SETERRQ(PETSC_COMM_SELF,e.errorCode,e.errorMessage);
    }
    PetscFunctionReturn(0);
  }

  PetscErrorCode Core::get_finite_element_entities_by_handle(const std::string name,Range &ents) const {
    MoABErrorCode rval;
    PetscFunctionBegin;
    try {
      EntityHandle meshset = get_finite_element_meshset(name);
      rval = moab.get_entities_by_handle(meshset,ents,true); CHKERRQ_MOAB(rval);
    } catch (MoFEMException const &e) {
      SETERRQ(PETSC_COMM_SELF,e.errorCode,e.errorMessage);
    }
    PetscFunctionReturn(0);
  }

  PetscErrorCode Core::list_finite_elements() const {
    PetscFunctionBegin;
    typedef FiniteElement_multiIndex::index<BitFEId_mi_tag>::type finiteElements_by_id;
    const finiteElements_by_id &BitFEId_set = finiteElements.get<BitFEId_mi_tag>();
    finiteElements_by_id::iterator miit = BitFEId_set.begin();
    for(;miit!=BitFEId_set.end();miit++) {
      std::ostringstream ss;
      ss << *miit << std::endl;
      PetscSynchronizedPrintf(comm,ss.str().c_str());
    }
    PetscSynchronizedFlush(comm,PETSC_STDOUT);
    PetscFunctionReturn(0);
  }

  PetscErrorCode Core::add_ents_to_finite_element_by_EDGEs(const Range& edges,const BitFEId id) {
    PetscFunctionBegin;
    *buildMoFEM &= 1<<0;
    const EntityHandle idm = get_finite_element_meshset(id);
    rval = moab.add_entities(idm,edges.subset_by_type(MBEDGE)); CHKERRQ_MOAB(rval);
    PetscFunctionReturn(0);
  }
  PetscErrorCode Core::add_ents_to_finite_element_by_EDGEs(const Range& edges,const std::string &name) {
    PetscFunctionBegin;
    *buildMoFEM &= 1<<0;
    try {
      ierr = seed_finite_elements(edges.subset_by_type(MBEDGE)); CHKERRQ(ierr);
      ierr = add_ents_to_finite_element_by_EDGEs(edges,getBitFEId(name));  CHKERRQ(ierr);
    } catch (MoFEMException const &e) {
      SETERRQ(comm,e.errorCode,e.errorMessage);
    }
    PetscFunctionReturn(0);
  }
  PetscErrorCode Core::add_ents_to_finite_element_by_VERTICEs(const Range& vert,const BitFEId id) {
    PetscFunctionBegin;
    *buildMoFEM &= 1<<0;
    const EntityHandle idm = get_finite_element_meshset(id);
    ierr = seed_finite_elements(vert.subset_by_type(MBVERTEX)); CHKERRQ(ierr);
    rval = moab.add_entities(idm,vert.subset_by_type(MBVERTEX)); CHKERRQ_MOAB(rval);
    PetscFunctionReturn(0);
  }
  PetscErrorCode Core::add_ents_to_finite_element_by_VERTICEs(const Range& vert,const std::string &name) {
    PetscFunctionBegin;
    *buildMoFEM &= 1<<0;
    try {
      ierr = add_ents_to_finite_element_by_VERTICEs(vert,getBitFEId(name));  CHKERRQ(ierr);
    } catch (MoFEMException const &e) {
      SETERRQ(comm,e.errorCode,e.errorMessage);
    }
    PetscFunctionReturn(0);
  }
  PetscErrorCode Core::add_ents_to_finite_element_by_TRIs(const Range& tris,const BitFEId id) {
    PetscFunctionBegin;
    *buildMoFEM &= 1<<0;
    const EntityHandle idm = get_finite_element_meshset(id);
    ierr = seed_finite_elements(tris.subset_by_type(MBTRI)); CHKERRQ(ierr);
    rval = moab.add_entities(idm,tris.subset_by_type(MBTRI)); CHKERRQ_MOAB(rval);
    PetscFunctionReturn(0);
  }
  PetscErrorCode Core::add_ents_to_finite_element_by_TRIs(const Range& tris,const std::string &name) {
    PetscFunctionBegin;
    *buildMoFEM &= 1<<0;
    try {
      ierr = add_ents_to_finite_element_by_TRIs(tris,getBitFEId(name));  CHKERRQ(ierr);
    } catch (MoFEMException const &e) {
      SETERRQ(comm,e.errorCode,e.errorMessage);
    }
    PetscFunctionReturn(0);
  }
  PetscErrorCode Core::add_ents_to_finite_element_by_TRIs(const EntityHandle meshset,const std::string &name,const bool recursive) {
    PetscFunctionBegin;
    *buildMoFEM &= 1<<0;
    EntityHandle idm = no_handle;
    try {
      idm = get_finite_element_meshset(getBitFEId(name));
    } catch (MoFEMException const &e) {
      SETERRQ(comm,e.errorCode,e.errorMessage);
    }
    Range tris;
    rval = moab.get_entities_by_type(meshset,MBTRI,tris,recursive); CHKERRQ_MOAB(rval);
    rval = moab.add_entities(idm,tris); CHKERRQ_MOAB(rval);
    PetscFunctionReturn(0);
  }
  PetscErrorCode Core::add_ents_to_finite_element_by_TETs(const EntityHandle meshset,const BitFEId id,const bool recursive) {
    PetscFunctionBegin;
    *buildMoFEM &= 1<<0;
    EntityHandle idm = no_handle;
    try {
      idm = get_finite_element_meshset(id);
    } catch (MoFEMException const &e) {
      SETERRQ(comm,e.errorCode,e.errorMessage);
    }
    Range tets;
    rval = moab.get_entities_by_type(meshset,MBTET,tets,recursive); CHKERRQ_MOAB(rval);
    rval = moab.add_entities(idm,tets); CHKERRQ_MOAB(rval);
    PetscFunctionReturn(0);
  }
  PetscErrorCode Core::add_ents_to_finite_element_by_TETs(const Range& tets,const BitFEId id) {
    PetscFunctionBegin;
    *buildMoFEM &= 1<<0;
    const EntityHandle idm = get_finite_element_meshset(id);
    rval = moab.add_entities(idm,tets.subset_by_type(MBTET)); CHKERRQ_MOAB(rval);
    PetscFunctionReturn(0);
  }
  PetscErrorCode Core::add_ents_to_finite_element_by_TETs(const Range& tets,const std::string &name) {
    PetscFunctionBegin;
    *buildMoFEM &= 1<<0;
    try {
      ierr = add_ents_to_finite_element_by_TETs(tets,getBitFEId(name));  CHKERRQ(ierr);
    } catch (MoFEMException const &e) {
      SETERRQ(comm,e.errorCode,e.errorMessage);
    }
    PetscFunctionReturn(0);
  }
  PetscErrorCode Core::add_ents_to_finite_element_by_TETs(const EntityHandle meshset,const std::string &name,const bool recursive) {
    PetscFunctionBegin;
    *buildMoFEM &= 1<<0;
    try {
      ierr = add_ents_to_finite_element_by_TETs(meshset,getBitFEId(name),recursive);  CHKERRQ(ierr);
    } catch (MoFEMException const &e) {
      SETERRQ(comm,e.errorCode,e.errorMessage);
    }
    PetscFunctionReturn(0);
  }

  PetscErrorCode Core::add_ents_to_finite_element_by_PRISMs(const EntityHandle meshset,const BitFEId id,const bool recursive) {
    PetscFunctionBegin;
    *buildMoFEM &= 1<<0;
    EntityHandle idm = no_handle;
    try {
      idm = get_finite_element_meshset(id);
    } catch (MoFEMException const &e) {
      SETERRQ(comm,e.errorCode,e.errorMessage);
    }
    Range prisms;
    rval = moab.get_entities_by_type(meshset,MBPRISM,prisms,recursive); CHKERRQ_MOAB(rval);
    rval = moab.add_entities(idm,prisms); CHKERRQ_MOAB(rval);
    PetscFunctionReturn(0);
  }

  PetscErrorCode Core::add_ents_to_finite_element_by_PRISMs(const Range& tets,const BitFEId id) {
    PetscFunctionBegin;
    *buildMoFEM &= 1<<0;
    const EntityHandle idm = get_finite_element_meshset(id);
    rval = moab.add_entities(idm,tets.subset_by_type(MBPRISM)); CHKERRQ_MOAB(rval);
    PetscFunctionReturn(0);
  }

  PetscErrorCode Core::add_ents_to_finite_element_by_PRISMs(const Range& prims,const std::string &name) {
    PetscFunctionBegin;
    *buildMoFEM &= 1<<0;
    try {
      ierr = add_ents_to_finite_element_by_PRISMs(prims,getBitFEId(name));  CHKERRQ(ierr);
    } catch (MoFEMException const &e) {
      SETERRQ(comm,e.errorCode,e.errorMessage);
    }
    PetscFunctionReturn(0);
  }

  PetscErrorCode Core::add_ents_to_finite_element_by_PRISMs(const EntityHandle meshset,const std::string &name,const bool recursive) {
    PetscFunctionBegin;
    *buildMoFEM &= 1<<0;
    try {
      ierr = add_ents_to_finite_element_by_PRISMs(meshset,getBitFEId(name),recursive);  CHKERRQ(ierr);
    } catch (MoFEMException const &e) {
      SETERRQ(comm,e.errorCode,e.errorMessage);
    }
    PetscFunctionReturn(0);
  }

  PetscErrorCode Core::add_ents_to_finite_element_EntType_by_bit_ref(const BitRefLevel &bit,const std::string &name,EntityType type,int verb) {
    PetscFunctionBegin;
    ierr = add_ents_to_finite_element_EntType_by_bit_ref(bit,BitRefLevel().set(),name,type,verb); CHKERRQ(ierr);
    PetscFunctionReturn(0);
  }

  PetscErrorCode Core::add_ents_to_finite_element_EntType_by_bit_ref(const BitRefLevel &bit,const BitRefLevel &mask,const std::string &name,EntityType type,int verb) {
    PetscFunctionBegin;
    try {
      if(verb==-1) verb = verbose;
      *buildMoFEM &= 1<<0;
      const BitFEId id = getBitFEId(name);
      const EntityHandle idm = get_finite_element_meshset(id);
      typedef RefElement_multiIndex::index<EntType_mi_tag>::type refMoabFE_by_type;
      refMoabFE_by_type &ref_MoFEMFiniteElement = refinedFiniteElements.get<EntType_mi_tag>();
      refMoabFE_by_type::iterator miit = ref_MoFEMFiniteElement.lower_bound(type);
      refMoabFE_by_type::iterator hi_miit = ref_MoFEMFiniteElement.upper_bound(type);
      if(verb > 1) {
        PetscSynchronizedPrintf(comm,"nb. ref elements in database %d\n",distance(miit,hi_miit));
      }
      int nb_add_FEs = 0;
      for(;miit!=hi_miit;miit++) {
        BitRefLevel bit2 = miit->getBitRefLevel();
        //check if all bits in mask are ib fe bit2
        //if((miit->getBitRefLevel()&bit)!=bit) continue;
        if((bit2&mask) != bit2) continue;
        if((bit2&bit).any()) {
          EntityHandle ent = miit->getRefEnt();
          rval = moab.add_entities(idm,&ent,1); CHKERRQ_MOAB(rval);
          nb_add_FEs++;
        }
      }
      if(verb > 0) {
        std::ostringstream ss;
        ss << "Add Nb. FEs " << nb_add_FEs << " form BitRef " << bit << std::endl;
        PetscSynchronizedPrintf(comm,"%s",ss.str().c_str());
        PetscSynchronizedFlush(comm,PETSC_STDOUT);
      }
    } catch (MoFEMException const &e) {
      SETERRQ(comm,e.errorCode,e.errorMessage);
    }
    PetscFunctionReturn(0);
  }

  PetscErrorCode Core::add_ents_to_finite_element_by_MESHSET(const EntityHandle meshset,const std::string& name,const bool recursive) {
    PetscFunctionBegin;
    *buildMoFEM &= 1<<0;
    const BitFEId id = getBitFEId(name);
    const EntityHandle idm = get_finite_element_meshset(id);
    if(recursive==false){
      rval = moab.add_entities(idm,&meshset,1); CHKERRQ_MOAB(rval);
    } else {
      Range meshsets;
      rval = moab.get_entities_by_type(meshset,MBENTITYSET,meshsets,false); CHKERRQ_MOAB(rval);
      rval = moab.add_entities(idm,meshsets); CHKERRQ_MOAB(rval);
    }
    PetscFunctionReturn(0);
  }

  PetscErrorCode Core::build_finite_elements(
    const boost::shared_ptr<FiniteElement> fe,const Range *ents_ptr,int verb
  ) {
    PetscFunctionBegin;
    if(verb==-1) verb = verbose;

    typedef RefElement_multiIndex::index<Ent_mi_tag>::type RefFiniteElementByEnt;
    typedef Field_multiIndex::index<BitFieldId_mi_tag>::type FieldById;
    FieldById &fields_by_id = fIelds.get<BitFieldId_mi_tag>();

    //get id of mofem fields for row, col and data
    enum IntLoop { ROW = 0,COL,DATA,LAST };
    BitFieldId fe_fields[LAST] = {
      fe.get()->getBitFieldIdRow(),
      fe.get()->getBitFieldIdCol(),
      fe.get()->getBitFieldIdData()
    };

    //get finite element meshset
    EntityHandle meshset = get_finite_element_meshset(fe.get()->getId());
    // get entities from finite element meshset // if meshset
    Range fe_ents;
    rval = moab.get_entities_by_handle(meshset,fe_ents,false); CHKERRQ_MOAB(rval);
    if(ents_ptr) fe_ents = intersect(fe_ents,*ents_ptr);

    // map entity uid to pointers
    std::map<UId,std::vector<boost::weak_ptr<EntFiniteElement> > > map_uid_fe;
    std::map<EntityHandle,int> data_dofs_size;

    //loop meshset Ents and add finite elements
    for(Range::iterator eit = fe_ents.begin();eit!=fe_ents.end();eit++) {

      // note: iterator is a wrapper
      // check if is in refinedFiniteElements database
      RefFiniteElementByEnt::iterator ref_fe_miit;
      ref_fe_miit = refinedFiniteElements.get<Ent_mi_tag>().find(*eit);
      if(ref_fe_miit == refinedFiniteElements.get<Ent_mi_tag>().end()) {
        std::ostringstream ss;
        ss << "ref FiniteElement not in database ent = " << *eit;
        ss << " type " << moab.type_from_handle(*eit);
        ss << " " << *fe;
        SETERRQ(comm,MOFEM_DATA_INCONSISTENCY,ss.str().c_str());
      }
      std::pair<EntFiniteElement_multiIndex::iterator,bool> p = entsFiniteElements.insert(
        boost::make_shared<EntFiniteElement>(ref_fe_miit->getRefElement(),fe)
      );

      if(fe_fields[ROW]==fe_fields[COL]) {
        p.first->get()->col_dof_view = p.first->get()->row_dof_view;
      }

      if(
        fe_fields[ROW]!=fe_fields[COL] &&
        p.first->get()->col_dof_view == p.first->get()->row_dof_view
      ) {
        p.first->get()->col_dof_view
        = boost::shared_ptr<DofEntity_multiIndex_uid_view>(new DofEntity_multiIndex_uid_view());
      }

      p.first->get()->row_dof_view->clear();
      p.first->get()->col_dof_view->clear();
      p.first->get()->data_dofs.clear();

      for(unsigned int ii = 0;ii<BitFieldId().size();ii++) {
        // Common field id for ROW, COL and DATA
        BitFieldId id_common = 0;
        // Check if the field (ii) is added to finite element
        for(int ss = 0;ss<LAST;ss++) {
          id_common |= fe_fields[ss]&BitFieldId().set(ii);
        }
        if( id_common.none() ) continue;
        // Find in database data associated with the field (ii)
        const BitFieldId field_id = BitFieldId().set(ii);
        FieldById::iterator miit = fields_by_id.find(field_id);
        if(miit==fields_by_id.end()) {
          SETERRQ(comm,MOFEM_DATA_INCONSISTENCY,"Data inconsistency");
        }

        // Entities adjacent to entities
        Range adj_ents;

        // Resolve entities on element, those entities are used to build tag with dof
        // uids on finite element tag
        ierr = p.first->get()->getElementAdjacency(*miit,adj_ents); CHKERRQ(ierr);

        // Loop over adjacencies of element and find field entities on those
        // adjacencies, that create hash map map_uid_fe which is used later
        std::string field_name = miit->get()->getName();
        for(Range::iterator eit = adj_ents.begin();eit!=adj_ents.end();eit++) {
          MoFEMEntity_multiIndex::index<Composite_Name_And_Ent_mi_tag>::type::iterator meit;
          meit = entsFields.get<Composite_Name_And_Ent_mi_tag>().find(boost::make_tuple(field_name,*eit));
          if(meit!=entsFields.get<Composite_Name_And_Ent_mi_tag>().end()) {
            UId uid = meit->get()->getGlobalUniqueId();
            map_uid_fe[uid].push_back(*p.first);
            if((field_id&p.first->get()->getBitFieldIdData()).any()) {
              data_dofs_size[p.first->get()->getEnt()] += meit->get()->getNbDofsOnEnt();
            }
          }
        }
      }

    }

    std::map<
      EntityHandle,
      boost::shared_ptr<std::vector<FEDofEntity> >
    > data_dofs_array;
    for(
      std::map<EntityHandle,int>::iterator
      mit = data_dofs_size.begin();mit!=data_dofs_size.end();mit++
    ) {
      if(mit->second>0) {
        data_dofs_array[mit->first]=boost::make_shared<std::vector<FEDofEntity> >();
        data_dofs_array[mit->first]->reserve(mit->second);
      }
    }
    std::vector<boost::shared_ptr<FEDofEntity> > data_dofs_shared_array;

    typedef DofEntity_multiIndex::index<Unique_Ent_mi_tag>::type DofsByEntUId;
    DofsByEntUId& dofs_by_ent_uid = dofsField.get<Unique_Ent_mi_tag>();

    // Loop over hash map, which has all entities on given elemnts
    boost::shared_ptr<SideNumber> side_number_ptr;
    for(
      std::map<UId,std::vector<boost::weak_ptr<EntFiniteElement> > >::iterator
      mit = map_uid_fe.begin();mit!=map_uid_fe.end();mit++
    ) {
      DofsByEntUId::iterator dit,hi_dit;
      dit = dofs_by_ent_uid.lower_bound(mit->first);
      hi_dit = dofs_by_ent_uid.upper_bound(mit->first);
      for(;dit!=hi_dit;dit++) {
        // cerr << mit->first << endl;
        // cerr << **dit << endl;
        const BitFieldId field_id = dit->get()->getId();
        const EntityHandle dof_ent = dit->get()->getEnt();
        std::vector<boost::weak_ptr<EntFiniteElement> >::iterator fe_it,hi_fe_it;
        fe_it = mit->second.begin();
        hi_fe_it = mit->second.end();
        for(;fe_it!=hi_fe_it;fe_it++) {

          // if rows and columns of finite element are the same, then
          // we exploit that case
          if((field_id&fe_it->lock().get()->getBitFieldIdRow()).any()) {
            fe_it->lock().get()->row_dof_view->insert(fe_it->lock().get()->row_dof_view->end(),*dit);
          }
          if(fe_it->lock().get()->col_dof_view!=fe_it->lock().get()->row_dof_view) {
            if((field_id&fe_it->lock().get()->getBitFieldIdCol()).any()) {
              fe_it->lock().get()->col_dof_view->insert(fe_it->lock().get()->col_dof_view->end(),*dit);
            }
          }

          // Add FEDofEntity, first create dofs, one by one, note that memory
          // is already reserved. Then create shared pointers and finally add th_FEName
          // to element multi-index
          const EntityHandle fe_ent = fe_it->lock().get()->getEnt();
          if(
            data_dofs_size[fe_ent]!=0 &&
            (field_id&fe_it->lock().get()->getBitFieldIdData()).any()
          ) {
            // There are data dofs on this element

            side_number_ptr = fe_it->lock().get()->getSideNumberPtr(dof_ent);
            data_dofs_array[fe_ent]->push_back(FEDofEntity(side_number_ptr,*dit));
            if(
              data_dofs_array[fe_ent]->size()==data_dofs_size[fe_ent]
            ) {
              // That means that FEDofEntity vector is full, and can be added to
              // multi-index

              // Create shared pointers vector
              data_dofs_shared_array.clear();
              data_dofs_shared_array.reserve(data_dofs_size[fe_ent]);
              for(
                std::vector<FEDofEntity>::iterator
                vit=data_dofs_array[fe_ent]->begin();
                vit!=data_dofs_array[fe_ent]->end();
                vit++
              ) {
                data_dofs_shared_array.push_back(
                  boost::shared_ptr<FEDofEntity>(data_dofs_array[fe_ent],&(*vit))
                );
              }
              fe_it->lock().get()->data_dofs.get<Unique_mi_tag>().insert(
                data_dofs_shared_array.begin(),data_dofs_shared_array.end()
              );
              fe_it->lock().get()->getDofsSeqence()=data_dofs_array[fe_ent];
            }
            // // add dofs to finite element multi_index database
            // fe_it->lock().get()->data_dofs.get<Unique_mi_tag>().insert(
            //   boost::shared_ptr<FEDofEntity>(new FEDofEntity(side_number_ptr,*dit))
            // );
          }
        }
      }

    }

    PetscFunctionReturn(0);
  }

  PetscErrorCode Core::build_finite_elements(int verb) {
    PetscFunctionBegin;
    if(verb==-1) verb = verbose;

    FiniteElement_multiIndex::iterator fe_miit = finiteElements.begin();

    // loop Finite Elements
    for(;fe_miit!=finiteElements.end();fe_miit++) {
      if(verb>0) PetscPrintf(comm,"Build Finite Elements %s\n",(*fe_miit)->getName().c_str());
      ierr = build_finite_elements(*fe_miit,NULL,verb); CHKERRQ(ierr);
    }

    if(verb>0) {
      PetscSynchronizedPrintf(comm,"Nb. FEs %u\n",entsFiniteElements.size());
      PetscSynchronizedFlush(comm,PETSC_STDOUT);
      typedef EntFiniteElement_multiIndex::index<BitFEId_mi_tag>::type FiniteElementById;
      FiniteElementById &finite_elements_by_id = entsFiniteElements.get<BitFEId_mi_tag>();
      FiniteElement_multiIndex::iterator fe_id_it = finiteElements.begin();
      for(;fe_id_it!=finiteElements.end();fe_id_it++) {
        FiniteElementById::iterator miit = finite_elements_by_id.lower_bound((*fe_id_it)->getId());
        FiniteElementById::iterator hi_miit = finite_elements_by_id.upper_bound((*fe_id_it)->getId());
        int count = distance(miit,hi_miit);
        std::ostringstream ss;
        ss << *(*fe_id_it) << " Nb. FEs " << count << std::endl;
        PetscSynchronizedPrintf(comm,ss.str().c_str());
        PetscSynchronizedFlush(comm,PETSC_STDOUT);
      }
    }

    *buildMoFEM |= 1<<1;
    PetscFunctionReturn(0);
  }

  PetscErrorCode Core::build_finite_elements(const BitRefLevel &bit,int verb) {
    PetscFunctionBegin;
    SETERRQ(comm,MOFEM_NOT_IMPLEMENTED,"Not yet implemented");
    PetscFunctionReturn(0);
  }

  PetscErrorCode Core::build_finite_elements(const string fe_name,const Range *ents_ptr,int verb) {
    PetscFunctionBegin;
    if(verb==-1) verb = verbose;

    FiniteElement_multiIndex::index<FiniteElement_name_mi_tag>::type::iterator
    fe_miit = finiteElements.get<FiniteElement_name_mi_tag>().find(fe_name);
    if(fe_miit==finiteElements.get<FiniteElement_name_mi_tag>().end()) {
      SETERRQ1(
        comm,MOFEM_NOT_FOUND,
        "Finite element <%s> not found",fe_name.c_str()
      );
    }

    if(verb>0) PetscPrintf(comm,"Build Finite Elements %s\n",fe_name.c_str());
    ierr = build_finite_elements(*fe_miit,ents_ptr,verb); CHKERRQ(ierr);
    if(verb>0) {
      typedef EntFiniteElement_multiIndex::index<BitFEId_mi_tag>::type FiniteElementById;
      FiniteElementById &finite_elements_by_id = entsFiniteElements.get<BitFEId_mi_tag>();
      FiniteElementById::iterator miit = finite_elements_by_id.lower_bound((*fe_miit)->getId());
      FiniteElementById::iterator hi_miit = finite_elements_by_id.upper_bound((*fe_miit)->getId());
      int count = distance(miit,hi_miit);
      std::ostringstream ss;
      ss << *(*fe_miit) << " Nb. FEs " << count << std::endl;
      PetscSynchronizedPrintf(comm,ss.str().c_str());
      PetscSynchronizedFlush(comm,PETSC_STDOUT);
    }

    *buildMoFEM |= 1<<1;
    PetscFunctionReturn(0);
  }

  PetscErrorCode Core::build_adjacencies(const Range &ents,int verb) {
    PetscFunctionBegin;
    if(verb==-1) verb = verbose;
    if(!((*buildMoFEM)&(1<<0))) SETERRQ(comm,MOFEM_NOT_FOUND,"field not build");
    if(!((*buildMoFEM)&(1<<1))) SETERRQ(comm,MOFEM_NOT_FOUND,"fe not build");
    EntFiniteElement_multiIndex::iterator fit = entsFiniteElements.begin();
    for(;fit!=entsFiniteElements.end();fit++) {
      if(
        (*fit)->getBitFieldIdRow().none() &&
        (*fit)->getBitFieldIdCol().none() &&
        (*fit)->getBitFieldIdData().none()
      ) continue;
      if(!ents.empty()) {
        if(ents.find((*fit)->getEnt())==ents.end()) continue;
      }
      int by = BYROW;
      if((*fit)->getBitFieldIdRow()!=(*fit)->getBitFieldIdCol()) by |= BYCOL;
      if((*fit)->getBitFieldIdRow()!=(*fit)->getBitFieldIdData()) by |= BYDATA;
      MoFEMEntityEntFiniteElementAdjacencyMap_change_ByWhat modify_row(by);
      GlobalUId ent_uid = UId(0);
      DofEntity_multiIndex_uid_view::iterator rvit;
      rvit = (*fit)->row_dof_view->begin();
      for(;rvit!=(*fit)->row_dof_view->end();rvit++) {
        if(ent_uid == (*rvit)->getMoFEMEntityPtr()->getGlobalUniqueId()) continue;
        ent_uid = (*rvit)->getMoFEMEntityPtr()->getGlobalUniqueId();
        std::pair<MoFEMEntityEntFiniteElementAdjacencyMap_multiIndex::iterator,bool> p;
        p = entFEAdjacencies.insert(
          MoFEMEntityEntFiniteElementAdjacencyMap((*rvit)->getMoFEMEntityPtr(),*fit)
        );
        bool success = entFEAdjacencies.modify(p.first,modify_row);
        if(!success) SETERRQ(comm,MOFEM_OPERATION_UNSUCCESSFUL,"modification unsuccessful");
      }
      if((*fit)->getBitFieldIdRow()!=(*fit)->getBitFieldIdCol()) {
        int by = BYCOL;
        if((*fit)->getBitFieldIdCol()!=(*fit)->getBitFieldIdData()) by |= BYDATA;
        MoFEMEntityEntFiniteElementAdjacencyMap_change_ByWhat modify_col(by);
        ent_uid = UId(0);
        DofEntity_multiIndex_uid_view::iterator cvit;
        cvit = (*fit)->col_dof_view->begin();
        for(;cvit!=(*fit)->col_dof_view->end();cvit++) {
          if( ent_uid == (*cvit)->getMoFEMEntityPtr()->getGlobalUniqueId()) continue;
          ent_uid = (*cvit)->getMoFEMEntityPtr()->getGlobalUniqueId();
          std::pair<MoFEMEntityEntFiniteElementAdjacencyMap_multiIndex::iterator,bool> p;
          p = entFEAdjacencies.insert(
            MoFEMEntityEntFiniteElementAdjacencyMap((*cvit)->getMoFEMEntityPtr(),*fit)
          );
          bool success = entFEAdjacencies.modify(p.first,modify_col);
          if(!success) SETERRQ(comm,MOFEM_OPERATION_UNSUCCESSFUL,"modification unsuccessful");
        }
      }
      if(
        (*fit)->getBitFieldIdRow()!=(*fit)->getBitFieldIdData()||
        (*fit)->getBitFieldIdCol()!=(*fit)->getBitFieldIdData()
      ) {
        MoFEMEntityEntFiniteElementAdjacencyMap_change_ByWhat modify_data(BYDATA);
        ent_uid = UId(0);
        FEDofEntity_multiIndex::iterator dvit;
        dvit = (*fit)->data_dofs.begin();
        for(;dvit!=(*fit)->data_dofs.end();dvit++) {
          if(ent_uid == (*dvit)->getMoFEMEntityPtr()->getGlobalUniqueId()) continue;
          ent_uid = (*dvit)->getMoFEMEntityPtr()->getGlobalUniqueId();
          std::pair<MoFEMEntityEntFiniteElementAdjacencyMap_multiIndex::iterator,bool> p;
          p = entFEAdjacencies.insert(
            MoFEMEntityEntFiniteElementAdjacencyMap((*dvit)->getMoFEMEntityPtr(),*fit)
          );
          bool success = entFEAdjacencies.modify(p.first,modify_data);
          if(!success) SETERRQ(comm,MOFEM_OPERATION_UNSUCCESSFUL,"modification unsuccessful");
        }
      }
    }
    if(verbose>1) {
      list_adjacencies();
    }
    if(verbose>0) {
      PetscSynchronizedPrintf(comm,"Nb. entFEAdjacencies %u\n",entFEAdjacencies.size());
      PetscSynchronizedFlush(comm,PETSC_STDOUT);
    }
    *buildMoFEM |= 1<<2;
    PetscFunctionReturn(0);
  }

  PetscErrorCode Core::build_adjacencies(const BitRefLevel &bit,const BitRefLevel &mask,int verb) {
    PetscFunctionBegin;
    if(verb==-1) verb = verbose;
    Range ents;
    ierr = get_entities_by_ref_level(bit,mask,ents); CHKERRQ(ierr);
    ierr = build_adjacencies(ents,verb); CHKERRQ(ierr);
    PetscFunctionReturn(0);
  }
  PetscErrorCode Core::build_adjacencies(const BitRefLevel &bit,int verb) {
    PetscFunctionBegin;
    if(verb==-1) verb = verbose;
    ierr = build_adjacencies(bit,BitRefLevel().set(),verb); CHKERRQ(ierr);
    PetscFunctionReturn(0);
  }

  PetscErrorCode Core::seed_finite_elements(const Range &entities,int verb) {
    PetscFunctionBegin;
    for(Range::iterator eit = entities.begin();eit!=entities.end();eit++) {
      RefEntity_multiIndex::index<Ent_mi_tag>::type::iterator
        eiit = refinedEntities.get<Ent_mi_tag>().find(*eit);
      if(eiit == refinedEntities.get<Ent_mi_tag>().end())  {
        SETERRQ(comm,MOFEM_NOT_FOUND,"entity is not in database");
      }
      if((*eiit)->getBitRefLevel().none()) continue;
      std::pair<RefElement_multiIndex::iterator,bool> p_MoFEMFiniteElement;
      switch ((*eiit)->getEntType()) {
        case MBVERTEX:
        p_MoFEMFiniteElement = refinedFiniteElements.insert(ptrWrapperRefElement(
          boost::shared_ptr<RefElement>(new RefElement_VERTEX(*eiit)))
        );
        break;
        case MBEDGE:
        p_MoFEMFiniteElement = refinedFiniteElements.insert(ptrWrapperRefElement(
          boost::shared_ptr<RefElement>(new RefElement_EDGE(*eiit)))
        );
        break;
        case MBTRI:
        p_MoFEMFiniteElement = refinedFiniteElements.insert(ptrWrapperRefElement(
          boost::shared_ptr<RefElement>(new RefElement_TRI(*eiit)))
        );
        break;
        case MBTET:
        p_MoFEMFiniteElement = refinedFiniteElements.insert(ptrWrapperRefElement(
          boost::shared_ptr<RefElement>(new RefElement_TET(*eiit)))
        );
        break;
        case MBPRISM:
        p_MoFEMFiniteElement = refinedFiniteElements.insert(ptrWrapperRefElement(
          boost::shared_ptr<RefElement>(new RefElement_PRISM(*eiit)))
        );
        break;
        default:
        SETERRQ(comm,MOFEM_DATA_INCONSISTENCY,"not implemented");
      }
    }
    PetscFunctionReturn(0);
  }

  PetscErrorCode Core::seed_finite_elements(const EntityHandle meshset,int verb) {
    PetscFunctionBegin;
    Range entities;
    ierr = moab.get_entities_by_handle(meshset,entities,true); CHKERRQ(ierr);
    ierr = seed_finite_elements(entities,verb); CHKERRQ(ierr);
    PetscFunctionReturn(0);
  }

  PetscErrorCode Core::get_finite_elements(const FiniteElement_multiIndex **finiteElements_ptr) const {
    PetscFunctionBegin;
    *finiteElements_ptr = &finiteElements;
    PetscFunctionReturn(0);
  }
  EntFiniteElementbyName::iterator
  Core::get_fe_by_name_begin(const std::string &fe_name) const {
    return entsFiniteElements.get<FiniteElement_name_mi_tag>().lower_bound(fe_name);
  }
  EntFiniteElementbyName::iterator
  Core::get_fe_by_name_end(const std::string &fe_name) const {
    return entsFiniteElements.get<FiniteElement_name_mi_tag>().upper_bound(fe_name);
  }

  PetscErrorCode Core::check_number_of_ents_in_ents_finite_element(const std::string& name) const {
    PetscFunctionBegin;
    FiniteElement_multiIndex::index<FiniteElement_name_mi_tag>::type::iterator it;
    it = finiteElements.get<FiniteElement_name_mi_tag>().find(name);
    if(it == finiteElements.get<FiniteElement_name_mi_tag>().end()) {
      SETERRQ1(comm,1,"finite element not found < %s >",name.c_str());
    }
    EntityHandle meshset = (*it)->getMeshset();
    MoABErrorCode rval;
    int num_entities;
    rval = moab.get_number_entities_by_handle(meshset,num_entities); CHKERRQ_MOAB(rval);
    if(
      entsFiniteElements.get<FiniteElement_name_mi_tag>().count((*it)->getName().c_str())
      != (unsigned int)num_entities
    ) {
      SETERRQ1(comm,1,"not equal number of entities in meshset and finite elements multiindex < %s >",(*it)->getName().c_str());
    }
    PetscFunctionReturn(0);
  }

  PetscErrorCode Core::check_number_of_ents_in_ents_finite_element() const {
    PetscFunctionBegin;
    FiniteElement_multiIndex::index<FiniteElement_name_mi_tag>::type::iterator it;
    it = finiteElements.get<FiniteElement_name_mi_tag>().begin();
    for(;it!=finiteElements.get<FiniteElement_name_mi_tag>().end();it++) {
      EntityHandle meshset = (*it)->getMeshset();
      MoABErrorCode rval;
      int num_entities;
      rval = moab.get_number_entities_by_handle(meshset,num_entities); CHKERRQ_MOAB(rval);
      if(entsFiniteElements.get<FiniteElement_name_mi_tag>().count((*it)->getName().c_str())
        != (unsigned int)num_entities) {
        SETERRQ1(comm,1,"not equal number of entities in meshset and finite elements multiindex < %s >",(*it)->getName().c_str());
      }
    }
    PetscFunctionReturn(0);
  }

}
