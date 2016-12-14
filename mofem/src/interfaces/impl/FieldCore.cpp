/** \file FieldCore.cpp
 * \brief Core interface methods for managing fields.
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
#include <FEMMultiIndices.hpp>
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

BitFieldId Core::get_BitFieldId(const std::string& name) const {
  typedef Field_multiIndex::index<FieldName_mi_tag>::type FieldSetByName;
  const FieldSetByName &set = fIelds.get<FieldName_mi_tag>();
  FieldSetByName::iterator miit = set.find(name);
  if(miit==set.end()) {
    THROW_MESSAGE("field < "+name+" > not in database (top tip: check spelling)");
  }
  return (*miit)->getId();
}
std::string Core::get_BitFieldId_name(const BitFieldId id) const {
  typedef Field_multiIndex::index<BitFieldId_mi_tag>::type FieldSetById;
  const FieldSetById &set = fIelds.get<BitFieldId_mi_tag>();
  FieldSetById::iterator miit = set.find(id);
  return (*miit)->getName();
}
EntityHandle Core::get_field_meshset(const BitFieldId id) const {
  typedef Field_multiIndex::index<BitFieldId_mi_tag>::type FieldSetById;
  const FieldSetById &set = fIelds.get<BitFieldId_mi_tag>();
  FieldSetById::iterator miit = set.find(id);
  if(miit==set.end()) THROW_MESSAGE("field not in database (top tip: check spelling)");
  return (*miit)->meshSet;
}
EntityHandle Core::get_field_meshset(const std::string& name) const {
  return get_field_meshset(get_BitFieldId(name));
}

bool Core::check_field(const std::string &name) const {
  typedef Field_multiIndex::index<FieldName_mi_tag>::type FieldSetByName;
  const FieldSetByName &set = fIelds.get<FieldName_mi_tag>();
  FieldSetByName::iterator miit = set.find(name);
  if(miit==set.end()) return false;
  return true;
}

const Field* Core::get_field_structure(const std::string& name) {
  typedef Field_multiIndex::index<FieldName_mi_tag>::type FieldSetByName;
  const FieldSetByName &set = fIelds.get<FieldName_mi_tag>();
  FieldSetByName::iterator miit = set.find(name);
  if(miit==set.end()) {
    throw MoFEMException(
      MOFEM_NOT_FOUND,
      std::string("field < "+name+" > not in databse (top tip: check spelling)").c_str()
    );
  }
  return miit->get();
}

PetscErrorCode Core::get_field_entities_by_dimension(const std::string name,int dim,Range &ents) const {
  MoABErrorCode rval;
  PetscFunctionBegin;
  try {
    EntityHandle meshset = get_field_meshset(name);
    rval = moab.get_entities_by_dimension(meshset,dim,ents,true); CHKERRQ_MOAB(rval);
  } catch (MoFEMException const &e) {
    SETERRQ(PETSC_COMM_SELF,e.errorCode,e.errorMessage);
  }
  PetscFunctionReturn(0);
}

PetscErrorCode Core::get_field_entities_by_type(const std::string name,EntityType type,Range &ents) const {
  MoABErrorCode rval;
  PetscFunctionBegin;
  try {
    EntityHandle meshset = get_field_meshset(name);
    rval = moab.get_entities_by_type(meshset,type,ents,true); CHKERRQ_MOAB(rval);
  } catch (MoFEMException const &e) {
    SETERRQ(PETSC_COMM_SELF,e.errorCode,e.errorMessage);
  }
  PetscFunctionReturn(0);
}

PetscErrorCode Core::get_field_entities_by_handle(const std::string name,Range &ents) const {
  MoABErrorCode rval;
  PetscFunctionBegin;
  try {
    EntityHandle meshset = get_field_meshset(name);
    rval = moab.get_entities_by_handle(meshset,ents,true); CHKERRQ_MOAB(rval);
  } catch (MoFEMException const &e) {
    SETERRQ(PETSC_COMM_SELF,e.errorCode,e.errorMessage);
  }
  PetscFunctionReturn(0);
}

PetscErrorCode Core::add_field(
  const std::string& name,
  const FieldSpace space,
  const FieldApproximationBase base,
  const FieldCoefficientsNumber nb_of_coefficients,
  const TagType tag_type,
  const enum MoFEMTypes bh,
  int verb
) {
  PetscFunctionBegin;
  if(verb==-1) verb = verbose;
  *buildMoFEM = 0;
  Field_multiIndex::index<FieldName_mi_tag>::type::iterator fit;
  fit = fIelds.get<FieldName_mi_tag>().find(name);
  if(fit != fIelds.get<FieldName_mi_tag>().end() ) {
    if(bh == MF_EXCL) {
      SETERRQ1(comm,MOFEM_OPERATION_UNSUCCESSFUL,"field is <%s> in database",name.c_str());
    }
  } else {
    EntityHandle meshset;
    rval = moab.create_meshset(MESHSET_SET|MESHSET_TRACK_OWNER,meshset); CHKERRQ_MOAB(rval);
    //id
    BitFieldId id = getFieldShift();
    rval = moab.tag_set_data(th_FieldId,&meshset,1,&id); CHKERRQ_MOAB(rval);
    //space
    rval = moab.tag_set_data(th_FieldSpace,&meshset,1,&space); CHKERRQ_MOAB(rval);
    //base
    rval = moab.tag_set_data(th_FieldBase,&meshset,1,&base); CHKERRQ_MOAB(rval);
    //name
    void const* tag_data[] = { name.c_str() };
    int tag_sizes[1]; tag_sizes[0] = name.size();
    rval = moab.tag_set_by_ptr(th_FieldName,&meshset,1,tag_data,tag_sizes); CHKERRQ_MOAB(rval);
    //name data prefix
    std::string name_data_prefix("_App_Data");
    void const* tag_prefix_data[] = { name_data_prefix.c_str() };
    int tag_prefix_sizes[1]; tag_prefix_sizes[0] = name_data_prefix.size();
    rval = moab.tag_set_by_ptr(th_FieldName_DataNamePrefix,&meshset,1,tag_prefix_data,tag_prefix_sizes); CHKERRQ_MOAB(rval);
    Tag th_AppOrder,th_FieldData,th_Rank,th_AppDofOrder,th_DofRank;
    //data
    std::string Tag_data_name = name_data_prefix+name;
    const int def_len = 0;
    rval = moab.tag_get_handle(
      Tag_data_name.c_str(),
      def_len,
      MB_TYPE_OPAQUE,
      th_FieldData,
      MB_TAG_CREAT|MB_TAG_BYTES|MB_TAG_VARLEN|MB_TAG_SPARSE,
      NULL
    ); CHKERRQ_MOAB(rval);
    //order
    ApproximationOrder def_ApproximationOrder = -1;
    std::string Tag_ApproximationOrder_name = "_App_Order_"+name;
    rval = moab.tag_get_handle(
      Tag_ApproximationOrder_name.c_str(),
      sizeof(ApproximationOrder),
      MB_TYPE_OPAQUE,
      th_AppOrder,MB_TAG_CREAT|MB_TAG_BYTES|tag_type,
      &def_ApproximationOrder
    ); CHKERRQ_MOAB(rval);
    //dof order
    std::string Tag_dof_ApproximationOrder_name = "_App_Dof_Order"+name;
    rval = moab.tag_get_handle(
      Tag_dof_ApproximationOrder_name.c_str(),
      def_len,MB_TYPE_OPAQUE,
      th_AppDofOrder,
      MB_TAG_CREAT|MB_TAG_BYTES|MB_TAG_VARLEN|MB_TAG_SPARSE,
      NULL
    ); CHKERRQ_MOAB(rval);
    //rank
    int def_rank = 1;
    std::string Tag_rank_name = "_Field_Rank_"+name;
    rval = moab.tag_get_handle(
      Tag_rank_name.c_str(),
      sizeof(FieldCoefficientsNumber),
      MB_TYPE_OPAQUE,
      th_Rank,
      MB_TAG_CREAT|MB_TAG_BYTES|tag_type,
      &def_rank
    ); CHKERRQ_MOAB(rval);
    rval = moab.tag_set_data(th_Rank,&meshset,1,&nb_of_coefficients); CHKERRQ_MOAB(rval);
    //dof rank
    std::string Tag_dof_rank_name = "_Field_Dof_Rank_"+name;
    rval = moab.tag_get_handle(
      Tag_dof_rank_name.c_str(),
      def_len,MB_TYPE_OPAQUE,
      th_DofRank,
      MB_TAG_CREAT|MB_TAG_BYTES|MB_TAG_VARLEN|MB_TAG_SPARSE,
      NULL
    ); CHKERRQ_MOAB(rval);
    //add meshset
    std::pair<Field_multiIndex::iterator,bool> p;
    try {
      CoordSystemsManager *cs_manger_ptr;
      ierr = query_interface(cs_manger_ptr); CHKERRQ(ierr);
      boost::shared_ptr<CoordSys > undefined_cs_ptr;
      ierr = cs_manger_ptr->getCoordSysPtr("UNDEFINED",undefined_cs_ptr); CHKERRQ(ierr);
      int sys_name_size[1];
      sys_name_size[0] = undefined_cs_ptr->getName().size();
      void const* sys_name[] = { &*undefined_cs_ptr->getNameRef().begin() };
      rval = moab.tag_set_by_ptr(
        cs_manger_ptr->get_th_CoordSysName(),&meshset,1,sys_name,sys_name_size
      ); CHKERRQ_MOAB(rval);
      EntityHandle coord_sys_id = undefined_cs_ptr->getMeshset();
      rval = moab.add_entities(coord_sys_id,&meshset,1); CHKERR_MOAB(rval);
      p = fIelds.insert(boost::shared_ptr<Field>(new Field(moab,meshset,undefined_cs_ptr)));
      if(bh == MF_EXCL) {
        if(!p.second) SETERRQ1(
          comm,MOFEM_NOT_FOUND,
          "field not inserted %s (top tip, it could be already there)",
          Field(moab,meshset,undefined_cs_ptr).getName().c_str()
        );
      }
    } catch (MoFEMException const &e) {
      SETERRQ(comm,e.errorCode,e.errorMessage);
    }
    if(verbose > 0) {
      std::ostringstream ss;
      ss << "add: " << **p.first << std::endl;
      PetscPrintf(comm,ss.str().c_str());
    }
  }
  //unt
  PetscFunctionReturn(0);
}

PetscErrorCode Core::add_ents_to_field_by_EDGEs(const Range &edges,const BitFieldId id,int verb) {
  PetscFunctionBegin;
  if(verb==-1) verb = verbose;
  *buildMoFEM = 0;
  EntityHandle idm = no_handle;
  try {
    idm = get_field_meshset(id);
  } catch (MoFEMException const &e) {
    SETERRQ(comm,e.errorCode,e.errorMessage);
  }
  FieldSpace space;
  rval = moab.tag_get_data(th_FieldSpace,&idm,1,&space); CHKERRQ_MOAB(rval);
  Range nodes;
  switch (space) {
    case L2:
      rval = moab.add_entities(idm,edges); CHKERRQ_MOAB(rval);
      if(verb>1) {
        std::ostringstream ss;
        ss << "add entities to field " << get_BitFieldId_name(id);
        ss << " nb. add edges " << edges.size();
        ss << std::endl;
        PetscPrintf(comm,ss.str().c_str());
      }
      break;
    case H1:
      rval = moab.add_entities(idm,edges); CHKERRQ_MOAB(rval);
      //rval = moab.get_connectivity(edges,nodes,true); CHKERRQ_MOAB(rval);
      //use get adjacencies, this will allow take in account adjacencies set user
      rval = moab.get_adjacencies(edges,0,false,nodes,moab::Interface::UNION); CHKERRQ_MOAB(rval);
      {
        Range topo_nodes;
        rval = moab.get_connectivity(edges,topo_nodes,true); CHKERRQ_MOAB(rval);
        Range mid_nodes;
        rval = moab.get_connectivity(edges,mid_nodes,false); CHKERRQ_MOAB(rval);
        mid_nodes = subtract(mid_nodes,topo_nodes);
        nodes = subtract(nodes,mid_nodes);
      }
      rval = moab.add_entities(idm,nodes); CHKERRQ_MOAB(rval);
      if(verb>1) {
        std::ostringstream ss;
        ss << "add entities to field " << get_BitFieldId_name(id);
        ss << " nb. add edges " << edges.size();
        ss << " nb. add nodes " << nodes.size();
        ss << std::endl;
        PetscPrintf(comm,ss.str().c_str());
      }
      break;
    case HCURL:

      // You would add hcurl space on edge here if you would have 1d h-curl space.
      // Ass far I know that really makes no sense. H-curl space is spanning on
      // edges, but you would add it by adding space through entities of higher
      // dimension.

      SETERRQ(
        comm,MOFEM_NOT_IMPLEMENTED,
        "sorry, not implemented, HCURL not implemented for edge"
      );

      break;
    case HDIV:

      // Look to comments above about h-curl spce. That apply as well to h-div space.

      SETERRQ(
        comm,MOFEM_NOT_IMPLEMENTED,
        "sorry, not implemented, HDIV not implemented for edge"
      );
      break;
    default:
      SETERRQ(comm,MOFEM_DATA_INCONSISTENCY,"sorry, adds unknown field to edge");
  }
  PetscFunctionReturn(0);
}
PetscErrorCode Core::add_ents_to_field_by_EDGEs(const Range &edges,const std::string& name,int verb) {
  PetscFunctionBegin;
  if(verb==-1) verb = verbose;
  *buildMoFEM = 0;
  try {
    ierr = add_ents_to_field_by_EDGEs(edges,get_BitFieldId(name),verb);  CHKERRQ(ierr);
  } catch (MoFEMException const &e) {
    SETERRQ(comm,e.errorCode,e.errorMessage);
  }
  PetscFunctionReturn(0);
}
PetscErrorCode Core::add_ents_to_field_by_EDGEs(const EntityHandle meshset,const BitFieldId id,int verb) {
  PetscFunctionBegin;
  if(verb==-1) verb = verbose;
  *buildMoFEM = 0;
  Range edges;
  rval = moab.get_entities_by_type(meshset,MBEDGE,edges,true); CHKERRQ_MOAB(rval);
  try {
    ierr = add_ents_to_field_by_EDGEs(edges,id,verb); CHKERRQ(ierr);
  } catch (MoFEMException const &e) {
    SETERRQ(comm,e.errorCode,e.errorMessage);
  }
  PetscFunctionReturn(0);
}
PetscErrorCode Core::add_ents_to_field_by_EDGEs(const EntityHandle meshset,const std::string& name,int verb) {
  PetscFunctionBegin;
  if(verb==-1) verb = verbose;
  *buildMoFEM = 0;
  try {
    ierr = add_ents_to_field_by_EDGEs(meshset,get_BitFieldId(name),verb);  CHKERRQ(ierr);
  } catch (MoFEMException const &e) {
    SETERRQ(comm,e.errorCode,e.errorMessage);
  }
  PetscFunctionReturn(0);
}
PetscErrorCode Core::add_ents_to_field_by_TRIs(const EntityHandle meshset,const BitFieldId id,int verb) {
  PetscFunctionBegin;
  Range tris;
  rval = moab.get_entities_by_type(meshset,MBTRI,tris,true); CHKERRQ_MOAB(rval);
  ierr = add_ents_to_field_by_TRIs(tris,id,verb); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
PetscErrorCode Core::add_ents_to_field_by_TRIs(const Range &tris,const BitFieldId id,int verb) {
  PetscFunctionBegin;
  if(verb==-1) verb = verbose;
  *buildMoFEM = 0;
  EntityHandle idm = no_handle;
  try {
    idm = get_field_meshset(id);
  } catch (MoFEMException const &e) {
    SETERRQ(comm,e.errorCode,e.errorMessage);
  }
  FieldSpace space;
  rval = moab.tag_get_data(th_FieldSpace,&idm,1,&space); CHKERRQ_MOAB(rval);
  Range nodes,edges;
  switch(space) {
    case L2:
    rval = moab.add_entities(idm,tris); CHKERRQ_MOAB(rval);
    if(verb>1) {
      std::ostringstream ss;
      ss << "add entities to field " << get_BitFieldId_name(id);
      ss << " nb. add tris " << tris.size();
      ss << std::endl;
      PetscPrintf(comm,ss.str().c_str());
    }
    break;
    case H1:
    rval = moab.add_entities(idm,tris); CHKERRQ_MOAB(rval);
    //rval = moab.get_connectivity(tris,nodes,true); CHKERRQ_MOAB(rval);
    //use get adjacencies, this will allow take in account adjacencies set user
    rval = moab.get_adjacencies(tris,0,false,nodes,moab::Interface::UNION); CHKERRQ_MOAB(rval);
    {
      Range topo_nodes;
      rval = moab.get_connectivity(tris,topo_nodes,true); CHKERRQ_MOAB(rval);
      Range mid_nodes;
      rval = moab.get_connectivity(tris,mid_nodes,false); CHKERRQ_MOAB(rval);
      mid_nodes = subtract(mid_nodes,topo_nodes);
      nodes = subtract(nodes,mid_nodes);
    }
    rval = moab.add_entities(idm,nodes); CHKERRQ_MOAB(rval);
    rval = moab.get_adjacencies(tris,1,false,edges,moab::Interface::UNION); CHKERRQ_MOAB(rval);
    rval = moab.add_entities(idm,edges); CHKERRQ_MOAB(rval);
    if(verb>1) {
      std::ostringstream ss;
      ss << "add entities to field " << get_BitFieldId_name(id);
      ss << " nb. add tris " << tris.size();
      ss << " nb. add edges " << edges.size();
      ss << " nb. add nodes " << nodes.size();
      ss << std::endl;
      PetscPrintf(comm,ss.str().c_str());
    }
    break;
    case HCURL:
    rval = moab.add_entities(idm,tris); CHKERRQ_MOAB(rval);
    rval = moab.get_adjacencies(tris,1,false,edges,moab::Interface::UNION); CHKERRQ_MOAB(rval);
    rval = moab.add_entities(idm,edges); CHKERRQ_MOAB(rval);
    if(verb>1) {
      std::ostringstream ss;
      ss << "add entities to field " << get_BitFieldId_name(id);
      ss << " nb. add tris " << tris.size();
      ss << " nb. add edges " << edges.size();
      ss << std::endl;
      PetscPrintf(comm,ss.str().c_str());
    }
    break;
    case HDIV:
    rval = moab.add_entities(idm,tris); CHKERRQ_MOAB(rval);
    if(verb>1) {
      std::ostringstream ss;
      ss << "add entities to field " << get_BitFieldId_name(id);
      ss << " nb. add tris " << tris.size();
      ss << std::endl;
      PetscPrintf(comm,ss.str().c_str());
    }
    break;
    default:
    SETERRQ(comm,MOFEM_DATA_INCONSISTENCY,"sorry, unknown field is applied to triangle entity");
  }
  PetscFunctionReturn(0);
}
PetscErrorCode Core::add_ents_to_field_by_TRIs(const EntityHandle meshset,const std::string& name,int verb) {
  PetscFunctionBegin;
  if(verb==-1) verb = verbose;
  *buildMoFEM = 0;
  try {
    ierr = add_ents_to_field_by_TRIs(meshset,get_BitFieldId(name),verb);  CHKERRQ(ierr);
  } catch (MoFEMException const &e) {
    SETERRQ(comm,e.errorCode,e.errorMessage);
  }
  PetscFunctionReturn(0);
}
PetscErrorCode Core::add_ents_to_field_by_TRIs(const Range &tris,const std::string& name,int verb) {
  PetscFunctionBegin;
  if(verb==-1) verb = verbose;
  *buildMoFEM = 0;
  try {
    ierr = add_ents_to_field_by_TRIs(tris,get_BitFieldId(name),verb);  CHKERRQ(ierr);
  } catch (MoFEMException const &e) {
    SETERRQ(comm,e.errorCode,e.errorMessage);
  }
  PetscFunctionReturn(0);
}
PetscErrorCode Core::add_ents_to_field_by_VERTICEs(const Range &nodes,const BitFieldId id,int verb) {
  PetscFunctionBegin;
  if(verb==-1) verb = verbose;
  *buildMoFEM = 0;
  EntityHandle idm = no_handle;
  try {
    idm = get_field_meshset(id);
  } catch (MoFEMException const &e) {
    SETERRQ(comm,e.errorCode,e.errorMessage);
  }
  FieldSpace space;
  rval = moab.tag_get_data(th_FieldSpace,&idm,1,&space); CHKERRQ_MOAB(rval);
  switch (space) {
    case L2:
    case H1:
    rval = moab.add_entities(idm,nodes); CHKERRQ_MOAB(rval);
    if(verb>1) {
      std::ostringstream ss;
      ss << "add entities to field " << get_BitFieldId_name(id);
      ss << " nb. add nodes " << nodes.size();
      ss << std::endl;
      PetscPrintf(comm,ss.str().c_str());
    }
    break;
    default:
    SETERRQ(
      comm,MOFEM_DATA_INCONSISTENCY,
      "sorry, if it try to span space on vertex entity which makes no sense");
  }
  PetscFunctionReturn(0);
}
PetscErrorCode Core::add_ents_to_field_by_VERTICEs(const EntityHandle meshset,const BitFieldId id,int verb) {
  PetscFunctionBegin;
  Range nodes;
  rval = moab.get_entities_by_type(meshset,MBVERTEX,nodes,true); CHKERRQ_MOAB(rval);
  ierr = add_ents_to_field_by_VERTICEs(nodes,id,verb); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
PetscErrorCode Core::add_ents_to_field_by_VERTICEs(const Range &nodes,const std::string& name,int verb) {
  PetscFunctionBegin;
  if(verb==-1) verb = verbose;
  *buildMoFEM = 0;
  try {
    ierr = add_ents_to_field_by_VERTICEs(nodes,get_BitFieldId(name),verb);  CHKERRQ(ierr);
  } catch (MoFEMException const &e) {
    SETERRQ(comm,e.errorCode,e.errorMessage);
  }
  PetscFunctionReturn(0);
}
PetscErrorCode Core::add_ents_to_field_by_VERTICEs(const EntityHandle meshset,const std::string& name,int verb) {
  PetscFunctionBegin;
  if(verb==-1) verb = verbose;
  *buildMoFEM = 0;
  try {
    ierr = add_ents_to_field_by_VERTICEs(meshset,get_BitFieldId(name),verb);  CHKERRQ(ierr);
  } catch (MoFEMException const &e) {
    SETERRQ(comm,e.errorCode,e.errorMessage);
  }
  PetscFunctionReturn(0);
}

PetscErrorCode Core::add_ents_to_field_by_TETs(const Range &tets,const BitFieldId id,int verb) {
  PetscFunctionBegin;
  if(verb==-1) verb = verbose;
  *buildMoFEM = 0;
  EntityHandle idm = no_handle;
  try {
    idm = get_field_meshset(id);
  } catch (MoFEMException const &e) {
    SETERRQ(comm,e.errorCode,e.errorMessage);
  }
  FieldSpace space;
  rval = moab.tag_get_data(th_FieldSpace,&idm,1,&space); CHKERRQ_MOAB(rval);
  Range nodes,tris,edges;
  switch(space) {
    case L2:
    rval = moab.add_entities(idm,tets); CHKERRQ_MOAB(rval);
    if(verb>1) {
      std::ostringstream ss;
      ss << "add entities to field " << get_BitFieldId_name(id);
      ss << " nb. add tets " << tets.size();
      ss << std::endl;
      PetscSynchronizedPrintf(comm,ss.str().c_str());
    }
    break;
    case H1:
    rval = moab.add_entities(idm,tets); CHKERRQ_MOAB(rval);
    //rval = moab.get_connectivity(tets,nodes,true); CHKERRQ_MOAB(rval);
    //use get adjacencies, this will allow take in account adjacencies set user
    rval = moab.get_adjacencies(tets,0,false,nodes,moab::Interface::UNION); CHKERRQ_MOAB(rval);
    {
      Range topo_nodes;
      rval = moab.get_connectivity(tets,topo_nodes,true); CHKERRQ_MOAB(rval);
      Range mid_nodes;
      rval = moab.get_connectivity(tets,mid_nodes,false); CHKERRQ_MOAB(rval);
      mid_nodes = subtract(mid_nodes,topo_nodes);
      nodes = subtract(nodes,mid_nodes);
    }
    rval = moab.add_entities(idm,nodes); CHKERRQ_MOAB(rval);
    rval = moab.get_adjacencies(tets,2,false,tris,moab::Interface::UNION); CHKERRQ_MOAB(rval);
    rval = moab.add_entities(idm,tris); CHKERRQ_MOAB(rval);
    rval = moab.get_adjacencies(tets,1,false,edges,moab::Interface::UNION); CHKERRQ_MOAB(rval);
    rval = moab.add_entities(idm,edges); CHKERRQ_MOAB(rval);
    if(verb>1) {
      std::ostringstream ss;
      ss << "add entities to field " << get_BitFieldId_name(id);
      ss << " nb. add tets " << tets.size();
      ss << " nb. add tris " << tris.size();
      ss << " nb. add edges " << edges.size();
      ss << " nb. add nodes " << nodes.size();
      ss << std::endl;
      PetscSynchronizedPrintf(comm,ss.str().c_str());
    }
    break;
    case HCURL:
    rval = moab.add_entities(idm,tets); CHKERRQ_MOAB(rval);
    rval = moab.get_adjacencies(tets,2,false,tris,moab::Interface::UNION); CHKERRQ_MOAB(rval);
    rval = moab.add_entities(idm,tris); CHKERRQ_MOAB(rval);
    rval = moab.get_adjacencies(tets,1,false,edges,moab::Interface::UNION); CHKERRQ_MOAB(rval);
    rval = moab.add_entities(idm,edges); CHKERRQ_MOAB(rval);
    if(verb>1) {
      std::ostringstream ss;
      ss << "add entities to field " << get_BitFieldId_name(id);
      ss << " nb. add tets " << tets.size();
      ss << " nb. add tris " << tris.size();
      ss << " nb. add edges " << edges.size();
      ss << std::endl;
      PetscSynchronizedPrintf(comm,ss.str().c_str());
    }
    break;
    case HDIV:
    rval = moab.add_entities(idm,tets); CHKERRQ_MOAB(rval);
    rval = moab.get_adjacencies(tets,2,false,tris,moab::Interface::UNION); CHKERRQ_MOAB(rval);
    rval = moab.add_entities(idm,tris); CHKERRQ_MOAB(rval);
    if(verb>1) {
      std::ostringstream ss;
      ss << "add entities to field " << get_BitFieldId_name(id);
      ss << " nb. add tets " << tets.size();
      ss << " nb. add tris " << tris.size();
      ss << std::endl;
      PetscSynchronizedPrintf(comm,ss.str().c_str());
    }
    break;
    default:
    SETERRQ(comm,MOFEM_DATA_INCONSISTENCY,"sorry, unknown space added to entity");
  }
  if(verb>1) {
    PetscSynchronizedFlush(comm,PETSC_STDOUT);
  }
  PetscFunctionReturn(0);
}
PetscErrorCode Core::add_ents_to_field_by_TETs(const EntityHandle meshset,const BitFieldId id,int verb) {
  PetscFunctionBegin;
  if(verb==-1) verb = verbose;
  Range tets;
  rval = moab.get_entities_by_type(meshset,MBTET,tets,true); CHKERRQ_MOAB(rval);
  if(verb>3) {
    PetscSynchronizedPrintf(comm,"nb. of tets %d\n",tets.size());
  }
  ierr = add_ents_to_field_by_TETs(tets,id,verb); CHKERRQ(ierr);
  if(verb>3) {
    PetscSynchronizedFlush(comm,PETSC_STDOUT);
  }
  PetscFunctionReturn(0);
}
PetscErrorCode Core::add_ents_to_field_by_TETs(const Range &tets,const std::string& name,int verb) {
  PetscFunctionBegin;
  if(verb==-1) verb = verbose;
  *buildMoFEM = 0;
  try {
    ierr = add_ents_to_field_by_TETs(tets,get_BitFieldId(name),verb);  CHKERRQ(ierr);
  } catch (MoFEMException const &e) {
    SETERRQ(comm,e.errorCode,e.errorMessage);
  }
  PetscFunctionReturn(0);
}
PetscErrorCode Core::add_ents_to_field_by_TETs(const EntityHandle meshset,const std::string& name,int verb) {
  PetscFunctionBegin;
  if(verb==-1) verb = verbose;
  *buildMoFEM = 0;
  try {
    ierr = add_ents_to_field_by_TETs(meshset,get_BitFieldId(name),verb);  CHKERRQ(ierr);
  } catch (MoFEMException const &e) {
    SETERRQ(comm,e.errorCode,e.errorMessage);
  }
  PetscFunctionReturn(0);
}
PetscErrorCode Core::add_ents_to_field_by_QUADs(const Range &quads,const BitFieldId id,int verb) {
  PetscFunctionBegin;
  if(verb==-1) verb = verbose;
  *buildMoFEM = 0;
  EntityHandle idm = no_handle;
  try {
    idm = get_field_meshset(id);
  } catch (MoFEMException const &e) {
    SETERRQ(comm,e.errorCode,e.errorMessage);
  }
  FieldSpace space;
  rval = moab.tag_get_data(th_FieldSpace,&idm,1,&space); CHKERRQ_MOAB(rval);
  Range nodes,faces,edges;
  switch(space) {
    case L2:
    rval = moab.add_entities(idm,quads); CHKERRQ_MOAB(rval);
    if(verb>1) {
      std::ostringstream ss;
      ss << "add entities to field " << get_BitFieldId_name(id);
      ss << " nb. add quads " << quads.size();
      ss << std::endl;
      PetscSynchronizedPrintf(comm,ss.str().c_str());
    }
    break;
    case H1:
    rval = moab.add_entities(idm,quads); CHKERRQ_MOAB(rval);
    //rval = moab.get_connectivity(quads,nodes,true); CHKERRQ_MOAB(rval);
    //use get adjacencies, this will allow take in account adjacencies set user
    rval = moab.get_adjacencies(quads,0,false,nodes,moab::Interface::UNION); CHKERRQ_MOAB(rval);
    {
      Range topo_nodes;
      rval = moab.get_connectivity(quads,topo_nodes,true); CHKERRQ_MOAB(rval);
      Range mid_nodes;
      rval = moab.get_connectivity(quads,mid_nodes,false); CHKERRQ_MOAB(rval);
      mid_nodes = subtract(mid_nodes,topo_nodes);
      nodes = subtract(nodes,mid_nodes);
    }
    rval = moab.add_entities(idm,nodes); CHKERRQ_MOAB(rval);
    rval = moab.get_adjacencies(quads,1,true,edges,moab::Interface::UNION); CHKERRQ_MOAB(rval);
    rval = moab.add_entities(idm,edges); CHKERRQ_MOAB(rval);
    if(verb>1) {
      std::ostringstream ss;
      ss << "add entities to field " << get_BitFieldId_name(id);
      ss << " nb. add quads " << quads.size();
      ss << " nb. add edges " << edges.size();
      ss << " nb. add nodes " << nodes.size();
      ss << std::endl;
      PetscSynchronizedPrintf(comm,ss.str().c_str());
    }
    break;
    case HCURL:
    SETERRQ(comm,MOFEM_NOT_IMPLEMENTED,"not implemented");
    break;
    case HDIV:
    SETERRQ(comm,MOFEM_NOT_IMPLEMENTED,"not implemented");
    break;
    default:
    SETERRQ(comm,MOFEM_DATA_INCONSISTENCY,"add_ents_to_field_by_TETs this field not work for TETs");
  }
  if(verb>1) {
    PetscSynchronizedFlush(comm,PETSC_STDOUT);
  }
  PetscFunctionReturn(0);
}
PetscErrorCode Core::add_ents_to_field_by_QUADs(const Range &quads,const std::string& name,int verb) {
  PetscFunctionBegin;
  if(verb==-1) verb = verbose;
  *buildMoFEM = 0;
  try {
    ierr = add_ents_to_field_by_QUADs(quads,get_BitFieldId(name),verb);  CHKERRQ(ierr);
  } catch (MoFEMException const &e) {
    SETERRQ(comm,e.errorCode,e.errorMessage);
  }
  PetscFunctionReturn(0);
}
PetscErrorCode Core::add_ents_to_field_by_QUADs(EntityHandle meshset,const std::string& name,int verb) {
  PetscFunctionBegin;
  if(verb==-1) verb = verbose;
  *buildMoFEM = 0;
  try {
    Range quads;
    rval = moab.get_entities_by_type(meshset,MBQUAD,quads,true); CHKERRQ_MOAB(rval);
    if(verb>3) {
      PetscSynchronizedPrintf(comm,"nb. of quads %d\n",quads.size());
    }
    ierr = add_ents_to_field_by_QUADs(quads,name,verb);  CHKERRQ(ierr);
    if(verb>3) {
      PetscSynchronizedFlush(comm,PETSC_STDOUT);
    }
  } catch (MoFEMException const &e) {
    SETERRQ(comm,e.errorCode,e.errorMessage);
  }
  PetscFunctionReturn(0);
}
PetscErrorCode Core::add_ents_to_field_by_PRISMs(const Range &prisms,const BitFieldId id,int verb) {
  PetscFunctionBegin;
  if(verb==-1) verb = verbose;
  *buildMoFEM = 0;
  EntityHandle idm = no_handle;
  try {
    idm = get_field_meshset(id);
  } catch (MoFEMException const &e) {
    SETERRQ(comm,e.errorCode,e.errorMessage);
  }
  FieldSpace space;
  rval = moab.tag_get_data(th_FieldSpace,&idm,1,&space); CHKERRQ_MOAB(rval);
  Range nodes,faces,edges;
  switch(space) {
    case L2:
    rval = moab.add_entities(idm,prisms); CHKERRQ_MOAB(rval);
    if(verb>1) {
      std::ostringstream ss;
      ss << "add entities to field " << get_BitFieldId_name(id);
      ss << " nb. add prisms " << prisms.size();
      ss << std::endl;
      PetscSynchronizedPrintf(comm,ss.str().c_str());
    }
    break;
    case H1:
    rval = moab.add_entities(idm,prisms); CHKERRQ_MOAB(rval);
    //rval = moab.get_connectivity(prisms,nodes,true); CHKERRQ_MOAB(rval);
    //use get adjacencies, this will allow take in account adjacencies set user
    rval = moab.get_adjacencies(prisms,0,false,nodes,moab::Interface::UNION); CHKERRQ_MOAB(rval);
    {
      Range topo_nodes;
      rval = moab.get_connectivity(prisms,topo_nodes,true); CHKERRQ_MOAB(rval);
      Range mid_nodes;
      rval = moab.get_connectivity(prisms,mid_nodes,false); CHKERRQ_MOAB(rval);
      mid_nodes = subtract(mid_nodes,topo_nodes);
      nodes = subtract(nodes,mid_nodes);
    }
    rval = moab.add_entities(idm,nodes); CHKERRQ_MOAB(rval);
    rval = moab.get_adjacencies(prisms,2,true,faces,moab::Interface::UNION); CHKERRQ_MOAB(rval);
    rval = moab.add_entities(idm,faces); CHKERRQ_MOAB(rval);
    rval = moab.get_adjacencies(prisms,1,true,edges,moab::Interface::UNION); CHKERRQ_MOAB(rval);
    rval = moab.add_entities(idm,edges); CHKERRQ_MOAB(rval);
    if(verb>1) {
      std::ostringstream ss;
      ss << "add entities to field " << get_BitFieldId_name(id);
      ss << " nb. add prisms " << prisms.size();
      ss << " nb. add faces " << faces.size();
      ss << " nb. add edges " << edges.size();
      ss << " nb. add nodes " << nodes.size();
      ss << std::endl;
      PetscSynchronizedPrintf(comm,ss.str().c_str());
    }
    break;
    case HCURL:
    SETERRQ(comm,MOFEM_NOT_IMPLEMENTED,"not implemented");
    break;
    case HDIV:
    SETERRQ(comm,MOFEM_NOT_IMPLEMENTED,"not implemented");
    break;
    default:
    SETERRQ(comm,MOFEM_DATA_INCONSISTENCY,"add_ents_to_field_by_TETs this field not work for TETs");
  }
  if(verb>1) {
    PetscSynchronizedFlush(comm,PETSC_STDOUT);
  }
  PetscFunctionReturn(0);
}
PetscErrorCode Core::add_ents_to_field_by_PRISMs(const Range &prisms,const std::string& name,int verb) {
  PetscFunctionBegin;
  if(verb==-1) verb = verbose;
  *buildMoFEM = 0;
  try {
    ierr = add_ents_to_field_by_PRISMs(prisms,get_BitFieldId(name),verb);  CHKERRQ(ierr);
  } catch (MoFEMException const &e) {
    SETERRQ(comm,e.errorCode,e.errorMessage);
  }
  PetscFunctionReturn(0);
}
PetscErrorCode Core::add_ents_to_field_by_PRISMs(EntityHandle meshset,const std::string& name,int verb) {
  PetscFunctionBegin;
  if(verb==-1) verb = verbose;
  *buildMoFEM = 0;
  try {
    Range prisms;
    rval = moab.get_entities_by_type(meshset,MBPRISM,prisms,true); CHKERRQ_MOAB(rval);
    if(verb>3) {
      PetscSynchronizedPrintf(comm,"nb. of prisms %d\n",prisms.size());
    }
    ierr = add_ents_to_field_by_PRISMs(prisms,name,verb);  CHKERRQ(ierr);
    if(verb>3) {
      PetscSynchronizedFlush(comm,PETSC_STDOUT);
    }
  } catch (MoFEMException const &e) {
    SETERRQ(comm,e.errorCode,e.errorMessage);
  }
  PetscFunctionReturn(0);
}
PetscErrorCode Core::set_field_order(const Range &ents,const BitFieldId id,const ApproximationOrder order,int verb) {
  PetscFunctionBegin;
  if(verb==-1) verb = verbose;
  *buildMoFEM = 0;

  //check field & meshset
  typedef Field_multiIndex::index<BitFieldId_mi_tag>::type FieldSetById;
  const FieldSetById &set_id = fIelds.get<BitFieldId_mi_tag>();
  FieldSetById::iterator miit = set_id.find(id);
  if(miit==set_id.end()) SETERRQ(comm,MOFEM_NOT_FOUND,"no filed found");
  EntityHandle idm;
  try {
   idm = get_field_meshset(id);
  } catch (MoFEMException const &e) {
    SETERRQ(comm,e.errorCode,e.errorMessage);
  }

  //intersection with field meshset
  Range ents_of_id_meshset;
  rval = moab.get_entities_by_handle(idm,ents_of_id_meshset,false); CHKERRQ_MOAB(rval);
  Range ents_ = intersect(ents,ents_of_id_meshset);
  if(verb>1) {
    PetscSynchronizedPrintf(
      comm,"nb. of ents for order change in the field <%s> %d\n",
      miit->get()->getName().c_str(),ents_.size()
  );
  }

  //ent view by field id (in set all MoabEnts has the same FieldId)
  typedef MoFEMEntity_multiIndex::index<BitFieldId_mi_tag>::type ent_set_by_id;
  ent_set_by_id& set = entsFields.get<BitFieldId_mi_tag>();
  ent_set_by_id::iterator eiit = set.lower_bound(id);
  MoFEMEntity_multiIndex_ent_view ents_id_view;
  if(eiit != set.end()) {
    ent_set_by_id::iterator hi_eiit = set.upper_bound(id);
    for(;eiit!=hi_eiit;eiit++) {
      ents_id_view.insert(*eiit);
    }
  }
  if(verb>1) {
    PetscSynchronizedPrintf(
      comm,"nb. of ents in the multi index field <%s> %d\n",
      miit->get()->getName().c_str(),ents_id_view.size()
    );
  }


  // get tags on entities
  std::vector<ApproximationOrder*> tag_data_order(ents_.size());
  rval = moab.tag_get_by_ptr(
    (*miit)->th_AppOrder,ents_,(const void **)&tag_data_order[0]
  ); CHKERRQ_MOAB(rval);

  //loop over ents
  int nb_ents_set_order_up = 0;
  int nb_ents_set_order_down = 0;
  int nb_ents_set_order_new = 0;
  Range::iterator eit = ents_.begin();
  for(unsigned int ee = 0;ee<ents_.size();ee++,eit++) {

    //sanity check
    switch((*miit)->getSpace()) {
      case H1:
      if(moab.type_from_handle(*eit)==MBVERTEX) {
        if(order!=1) {
          SETERRQ(comm,MOFEM_DATA_INCONSISTENCY,"approximation order for H1 space and vertex different than 1 makes not sense");
        }
      }
      break;
      case HCURL:
      if(moab.type_from_handle(*eit)==MBVERTEX) {
        SETERRQ(comm,MOFEM_DATA_INCONSISTENCY,"HDIV space on vertices makes no sense");
      }
      break;
      case HDIV:
      if(moab.type_from_handle(*eit)==MBVERTEX) {
        SETERRQ(comm,MOFEM_DATA_INCONSISTENCY,"HDIV space on vertices makes no sense");
      }
      if(moab.type_from_handle(*eit)==MBEDGE) {
        SETERRQ(comm,MOFEM_DATA_INCONSISTENCY,"HDIV space on edges makes no sense");
      }
      break;
      default:
      break;
    }

    MoFEMEntity_multiIndex_ent_view::iterator vit = ents_id_view.find(*eit);
    if(vit!=ents_id_view.end()) {

      //entity is in database and order is changed or reset
      const ApproximationOrder old_approximation_order = (*vit)->getMaxOrder();
      if(old_approximation_order==order) continue;
      MoFEMEntity_multiIndex::iterator miit = entsFields.get<Unique_mi_tag>().
      find((*vit)->getGlobalUniqueId());

      if((*miit)->getMaxOrder()<order) nb_ents_set_order_up++;
      if((*miit)->getMaxOrder()>order) nb_ents_set_order_down++;

      {

      	//set dofs inactive if order is reduced, and set new order to entity if
      	//order is increased (note that dofs are not build if order is
      	//increased)

        typedef DofEntityByNameAndEnt dof_set_type;
        dof_set_type& set_set = dofsField.get<Composite_Name_And_Ent_mi_tag>();
        dof_set_type::iterator dit = set_set.lower_bound(
          boost::make_tuple((*miit)->getNameRef(),(*miit)->getEnt())
        );
        dof_set_type::iterator hi_dit = set_set.upper_bound(
          boost::make_tuple((*miit)->getNameRef(),(*miit)->getEnt())
        );

        for(;dit!=hi_dit;dit++) {
          if((*dit)->getDofOrder()<=order) continue;
          bool success = dofsField.modify(dofsField.project<0>(dit),DofEntity_active_change(false));
          if(!success) SETERRQ(comm,MOFEM_OPERATION_UNSUCCESSFUL,"modification unsuccessful");
        }

        bool success = entsFields.modify(entsFields.project<0>(miit),MoFEMEntity_change_order(order));
        if(!success) SETERRQ(comm,MOFEM_OPERATION_UNSUCCESSFUL,"modification unsuccessful");

      }

    } else {

      //entity is not in database and order is changed or reset
      *tag_data_order[ee] = order;
      RefEntity_multiIndex::index<Ent_mi_tag>::type::iterator miit_ref_ent =
      refinedEntities.get<Ent_mi_tag>().find(*eit);
      if(miit_ref_ent==refinedEntities.get<Ent_mi_tag>().end()) {

        RefEntity ref_ent(basicEntityDataPtr,*eit);
        // FIXME: need some consistent policy in that case
        if(ref_ent.getBitRefLevel().none()) continue; // not on any mesh and not in database
        std::cerr << ref_ent << std::endl;
        std::cerr << "bit level " << ref_ent.getBitRefLevel() << std::endl;
        SETERRQ(comm,MOFEM_DATA_INCONSISTENCY,"Try to add entities which are not seeded or added to database");

      }

      try {

        // increase order
        boost::shared_ptr<MoFEMEntity> moabent(new MoFEMEntity(*miit,*miit_ref_ent));
        std::pair<MoFEMEntity_multiIndex::iterator,bool> e_miit = entsFields.insert(moabent);
        bool success = entsFields.modify(e_miit.first,MoFEMEntity_change_order(order));
        if(!success) SETERRQ(comm,MOFEM_OPERATION_UNSUCCESSFUL,"modification unsuccessful");
        nb_ents_set_order_new++;

      } catch (MoFEMException const &e) {
        SETERRQ(comm,e.errorCode,e.errorMessage);
      } catch (const std::exception& ex) {
        std::ostringstream ss;
        ss << ex.what() << std::endl;
        SETERRQ(comm,MOFEM_STD_EXCEPTION_THROW,ss.str().c_str());
      }

    }
  }

  if(verb>1) {
    PetscSynchronizedPrintf(
      comm,"nb. of entities in field <%s> for which order was increased %d (order %d)\n",
      miit->get()->getName().c_str(),nb_ents_set_order_up,order
    );
    PetscSynchronizedPrintf(
      comm,"nb. of entities in field <%s> for which order was reduced %d (order %d)\n",
      miit->get()->getName().c_str(),nb_ents_set_order_down,order
    );
    PetscSynchronizedPrintf(
      comm,"nb. of entities in field <%s> for which order set %d (order %d)\n",
      miit->get()->getName().c_str(),nb_ents_set_order_new,order
    );
    PetscSynchronizedFlush(comm,PETSC_STDOUT);
  }

  PetscFunctionReturn(0);
}
PetscErrorCode Core::set_field_order(const EntityHandle meshset,const EntityType type,const BitFieldId id,const ApproximationOrder order,int verb) {
  PetscFunctionBegin;
  if(verb==-1) verb = verbose;
  *buildMoFEM = 0;
  Range ents;
  rval = moab.get_entities_by_type(meshset,type,ents); CHKERRQ_MOAB(rval);
  if(verb>1) {
    PetscSynchronizedPrintf(comm,"nb. of ents for order change %d\n",ents.size());
  }
  try{
    ierr = set_field_order(ents,id,order,verb); CHKERRQ(ierr);
  } catch (MoFEMException const &e) {
    SETERRQ(comm,e.errorCode,e.errorMessage);
  }
  if(verb>1) {
    PetscSynchronizedFlush(comm,PETSC_STDOUT);
  }
  PetscFunctionReturn(0);
}
PetscErrorCode Core::set_field_order(const EntityHandle meshset,const EntityType type,const std::string& name,const ApproximationOrder order,int verb) {
  PetscFunctionBegin;
  if(verb==-1) verb = verbose;
  *buildMoFEM = 0;
  try{
    ierr = set_field_order(meshset,type,get_BitFieldId(name),order,verb); CHKERRQ(ierr);
  } catch (MoFEMException const &e) {
    SETERRQ(comm,e.errorCode,e.errorMessage);
  }
  PetscFunctionReturn(0);
}
PetscErrorCode Core::set_field_order(const Range &ents,const std::string& name,const ApproximationOrder order,int verb) {
  PetscFunctionBegin;
  if(verb==-1) verb = verbose;
  *buildMoFEM = 0;
  try{
    ierr = set_field_order(ents,get_BitFieldId(name),order,verb); CHKERRQ(ierr);
  } catch (MoFEMException const &e) {
    SETERRQ(comm,e.errorCode,e.errorMessage);
  }
  PetscFunctionReturn(0);
}
PetscErrorCode Core::set_field_order_by_entity_type_and_bit_ref(
  const BitRefLevel &bit,const BitRefLevel &mask,
  const EntityType type,const BitFieldId id,const ApproximationOrder order,int verb
) {
  PetscFunctionBegin;
  if(verb==-1) verb = verbose;
  *buildMoFEM = 0;
  Range ents;
  ierr = get_entities_by_type_and_ref_level(bit,mask,type,ents,verb); CHKERRQ(ierr);
  try{
    ierr = set_field_order(ents,id,order,verb); CHKERRQ(ierr);
  } catch (MoFEMException const &e) {
    SETERRQ(comm,e.errorCode,e.errorMessage);
  }
  PetscFunctionReturn(0);
}
PetscErrorCode Core::set_field_order_by_entity_type_and_bit_ref(
  const BitRefLevel &bit,const BitRefLevel &mask,
  const EntityType type,const std::string& name,const ApproximationOrder order,int verb
) {
  PetscFunctionBegin;
  if(verb==-1) verb = verbose;
  *buildMoFEM = 0;
  Range ents;
  ierr = get_entities_by_type_and_ref_level(bit,mask,type,ents,verb); CHKERRQ(ierr);
  try{
    ierr = set_field_order(ents,get_BitFieldId(name),order,verb); CHKERRQ(ierr);
  } catch (MoFEMException const &e) {
    SETERRQ(comm,e.errorCode,e.errorMessage);
  }
  PetscFunctionReturn(0);
}
PetscErrorCode Core::BuildFieldForNoField(
  const BitFieldId id,std::map<EntityType,int> &dof_counter,int verb
) {
  PetscFunctionBegin;
  if(verb==-1) verb = verbose;
  //field it
  typedef Field_multiIndex::index<BitFieldId_mi_tag>::type FieldSetById;
  const FieldSetById &set_id = fIelds.get<BitFieldId_mi_tag>();
  //find fiels
  FieldSetById::iterator miit = set_id.find(id);
  if(miit == set_id.end()) {
    SETERRQ(comm,MOFEM_NOT_FOUND,"field not found");
  }

  //ents in the field meshset
  Range ents_of_id_meshset;
  rval = moab.get_entities_by_handle((*miit)->meshSet,ents_of_id_meshset,false); CHKERRQ_MOAB(rval);
  if(verb>5) {
    PetscSynchronizedPrintf(
      comm,"ents in field %s meshset %d\n",(*miit)->getName().c_str(),ents_of_id_meshset.size()
    );
  }
  for(
    Range::iterator eit = ents_of_id_meshset.begin();
    eit != ents_of_id_meshset.end();eit++
  ) {
    //serch if field meshset is in database
    RefEntity_multiIndex::index<Ent_mi_tag>::type::iterator miit_ref_ent;
    miit_ref_ent = refinedEntities.get<Ent_mi_tag>().find(*eit);
    if(miit_ref_ent==refinedEntities.get<Ent_mi_tag>().end()) {
      SETERRQ(comm,MOFEM_DATA_INCONSISTENCY,"database inconsistency");
    }
    std::pair<MoFEMEntity_multiIndex::iterator,bool> e_miit;
    try {
      //create database entity
      boost::shared_ptr<MoFEMEntity> moabent(new MoFEMEntity(*miit,*miit_ref_ent));
      e_miit = entsFields.insert(moabent);
    } catch (MoFEMException const &e) {
      SETERRQ(comm,e.errorCode,e.errorMessage);
    } catch (const std::exception& ex) {
      std::ostringstream ss;
      ss << ex.what() << std::endl;
      SETERRQ(comm,MOFEM_STD_EXCEPTION_THROW,ss.str().c_str());
    }
    //this is nor real field in space (set order to zero)
    bool success = entsFields.modify(e_miit.first,MoFEMEntity_change_order(0));
    if(!success) SETERRQ(comm,MOFEM_OPERATION_UNSUCCESSFUL,"modification unsuccessful");
    FieldCoefficientsNumber rank = 0;
    //create dofs on this entity (nb. of dofs is equal to rank)
    for(;rank<(*e_miit.first)->getNbOfCoeffs();rank++) {
      std::pair<DofEntity_multiIndex::iterator,bool> d_miit;
      //check if dof is in darabase
      d_miit.first = dofsField.project<0>(
        dofsField.get<Unique_mi_tag>().find(DofEntity::getGlobalUniqueIdCalculate(rank,*(e_miit.first)))
      );
      //if dof is not in databse
      if(d_miit.first==dofsField.end()) {
        //insert dof
        d_miit = dofsField.insert(boost::shared_ptr<DofEntity>(new DofEntity(*(e_miit.first),0,rank,rank)));
        if(d_miit.second) {
          dof_counter[MBENTITYSET]++; // Count entities in the meshset
        }
        bool success = dofsField.modify(d_miit.first,DofEntity_active_change(true));
        if(!success) SETERRQ(comm,MOFEM_OPERATION_UNSUCCESSFUL,"modification unsuccessful");
      }
      //check consistency
      assert((*d_miit.first)->getEntType()==(*e_miit.first)->getEntType());
      assert((*d_miit.first)->getId()==(*e_miit.first)->getId());
      assert((*d_miit.first)->getMaxOrder()==0);
    }
  }
  if(verb>2) {
    typedef DofEntity_multiIndex::index<BitFieldId_mi_tag>::type dof_set_by_id;
    dof_set_by_id &set = dofsField.get<BitFieldId_mi_tag>();
    dof_set_by_id::iterator miit2 = set.lower_bound(id);
    dof_set_by_id::iterator hi_miit2 = set.upper_bound(id);
    assert(miit2!=hi_miit2);
    for(;miit2!=hi_miit2;miit2++) {
      std::ostringstream ss;
      ss << *miit2 << std::endl;;
      PetscSynchronizedPrintf(comm,ss.str().c_str());
    }
    PetscSynchronizedFlush(comm,PETSC_STDOUT);
  }
  PetscFunctionReturn(0);
}
PetscErrorCode Core::BuildFieldForL2H1HcurlHdiv(
  const BitFieldId id,
  std::map<EntityType,int> &dof_counter,
  std::map<EntityType,int> &inactive_dof_counter,
  int verb
) {
  PetscFunctionBegin;
  if(verb==-1) verb = verbose;
  //field it
  typedef Field_multiIndex::index<BitFieldId_mi_tag>::type FieldSetById;
  //find field
  const FieldSetById &set_id = fIelds.get<BitFieldId_mi_tag>();
  FieldSetById::iterator miit = set_id.find(id);
  if(miit == set_id.end()) {
    SETERRQ(comm,MOFEM_NOT_FOUND,"field not found");
  }
  //ents in the field meshset
  Range ents_of_id_meshset;
  rval = moab.get_entities_by_handle((*miit)->meshSet,ents_of_id_meshset,false); CHKERRQ_MOAB(rval);
  if(verb>5) {
    PetscSynchronizedPrintf(
      comm,"ents in field %s meshset %d\n",(*miit)->getName().c_str(),ents_of_id_meshset.size()
    );
  }

  MoFEMEntity_multiIndex::iterator eit_insert_hint = entsFields.end();
  DofEntity_multiIndex::iterator dit_insert_hint = dofsField.end();

  //create dofsField
  Range::iterator eit = ents_of_id_meshset.begin();
  for(;eit!=ents_of_id_meshset.end();eit++) {

    // check if ent is in ref meshset
    RefEntity_multiIndex::index<Ent_mi_tag>::type::iterator miit_ref_ent;
    miit_ref_ent = refinedEntities.get<Ent_mi_tag>().find(*eit);
    if(miit_ref_ent==refinedEntities.get<Ent_mi_tag>().end()) {
      RefEntity ref_ent(basicEntityDataPtr,*eit);
      if(ref_ent.getBitRefLevel().none()) {
        continue; // not on any mesh and not in database
      }
      std::cerr << ref_ent << std::endl;
      std::cerr << "bit level " << ref_ent.getBitRefLevel() << std::endl;
      SETERRQ(comm,MOFEM_DATA_INCONSISTENCY,"database inconsistency");
    }

    boost::shared_ptr<MoFEMEntity> moabent(new MoFEMEntity(*miit,*miit_ref_ent));
    // create mofem entity linked to ref ent
    MoFEMEntity_multiIndex::iterator e_miit;
    try {
      e_miit = entsFields.find(moabent->getGlobalUniqueId());
    } catch (MoFEMException const &e) {
      SETERRQ(comm,e.errorCode,e.errorMessage);
    } catch (const std::exception& ex) {
      std::ostringstream ss;
      ss << ex.what() << std::endl;
      SETERRQ(comm,MOFEM_OPERATION_UNSUCCESSFUL,ss.str().c_str());
    }
    // add field entity if not exist
    if(e_miit == entsFields.end()) {
      ApproximationOrder order = -1;
      rval = moab.tag_set_data((*miit)->th_AppOrder,&*eit,1,&order); CHKERRQ_MOAB(rval);
      try {
        // boost::shared_ptr<MoFEMEntity> moabent(new MoFEMEntity(*miit,*miit_ref_ent));
        e_miit = entsFields.insert(eit_insert_hint,moabent);
      } catch (MoFEMException const &e) {
        SETERRQ(comm,e.errorCode,e.errorMessage);
      } catch (const std::exception& ex) {
        std::ostringstream ss;
        ss << ex.what() << std::endl;
        SETERRQ(comm,MOFEM_OPERATION_UNSUCCESSFUL,ss.str().c_str());
      }
      bool success = entsFields.modify(e_miit,MoFEMEntity_change_order(-1));
      if(!success) SETERRQ(comm,MOFEM_OPERATION_UNSUCCESSFUL,"modification unsuccessful");
    }
    eit_insert_hint = e_miit;
    eit_insert_hint++;

    // insert dofmoabent into mofem databse
    int nb_dofs_on_ent = (*e_miit)->getNbDofsOnEnt();
    int nb_active_dosf_on_ent = (*e_miit)->getNbOfCoeffs()*(*e_miit)->getOrderNbDofs((*e_miit)->getMaxOrder());
    int DD = 0;
    int oo = 0;
    // loop orders (loop until max entity order is set)
    for(;oo<=(*e_miit)->getMaxOrder()||DD<nb_dofs_on_ent;oo++) {
      //loop nb. dofs at order oo
      for(int dd = 0;dd<(*e_miit)->getOrderNbDofsDiff(oo);dd++) {
        //loop rank
        for(int rr = 0;rr<(*e_miit)->getNbOfCoeffs();rr++,DD++) {
          try {
            boost::shared_ptr<DofEntity> mdof(new DofEntity(*(e_miit),oo,rr,DD));
            DofEntity_multiIndex::iterator d_miit;
            d_miit = dofsField.insert(dit_insert_hint,mdof);
            dit_insert_hint = d_miit; // hint for next insertion
            dit_insert_hint++;
            bool is_active;
            if(DD<nb_active_dosf_on_ent) {
              is_active = true;
              dof_counter[(*d_miit)->getEntType()]++;
            } else {
              is_active = false;
              inactive_dof_counter[(*d_miit)->getEntType()]++;
            }
            bool success = dofsField.modify(d_miit,DofEntity_active_change(is_active));
            if(!success) SETERRQ(comm,MOFEM_OPERATION_UNSUCCESSFUL,"modification unsuccessful");
            //check ent
            if((*d_miit)->getEnt()!=(*e_miit)->getEnt()) {
              SETERRQ(comm,MOFEM_DATA_INCONSISTENCY,"data inconsistency");
            }
            if((*d_miit)->getEntType()!=(*e_miit)->getEntType()) {
              SETERRQ(comm,MOFEM_DATA_INCONSISTENCY,"data inconsistency");
            }
            if((*d_miit)->getId()!=(*e_miit)->getId()) {
              SETERRQ(comm,MOFEM_DATA_INCONSISTENCY,"data inconsistency");
            }
            //check dof
            if((*d_miit)->getDofOrder()!=oo) {
              std::ostringstream ss;
              ss << "data inconsistency!" << std::endl;
              ss << "should be " << mdof << std::endl;
              ss << "but is " << *(*d_miit) << std::endl;
              SETERRQ(comm,MOFEM_DATA_INCONSISTENCY,ss.str().c_str());
            }
            if((*d_miit)->getMaxOrder()!=(*e_miit)->getMaxOrder()) {
              SETERRQ(comm,MOFEM_DATA_INCONSISTENCY,"data inconsistency");
            }
          } catch (MoFEMException const &e) {
            SETERRQ(comm,e.errorCode,e.errorMessage);
          }
        }
      }
    }
    if(DD != (*e_miit)->getNbDofsOnEnt()) {
      SETERRQ(comm,MOFEM_DATA_INCONSISTENCY,"data inconsistency");
    }
  }
  PetscFunctionReturn(0);
}
PetscErrorCode Core::build_fields(int verb) {
  PetscFunctionBegin;
  if(verb==-1) verb = verbose;
  typedef Field_multiIndex::index<BitFieldId_mi_tag>::type FieldSetById;
  FieldSetById &set_id = fIelds.get<BitFieldId_mi_tag>();
  FieldSetById::iterator miit = set_id.begin();
  for(;miit!=set_id.end();miit++) {
    std::map<EntityType,int> dof_counter;
    std::map<EntityType,int> inactive_dof_counter;
    if (verb > 0) {
      PetscSynchronizedPrintf(
        comm,"Build Field %s (rank %d)\n",(*miit)->getName().c_str(),rAnk
      );
    }
    switch ((*miit)->getSpace()) {
      case NOFIELD:
      ierr = BuildFieldForNoField((*miit)->getId(),dof_counter,verb); CHKERRQ(ierr);
      break;
      case L2:
      case H1:
      case HCURL:
      case HDIV:
      ierr = BuildFieldForL2H1HcurlHdiv((*miit)->getId(),dof_counter,inactive_dof_counter,verb); CHKERRQ(ierr);
      break;
      default:
      SETERRQ(comm,MOFEM_NOT_IMPLEMENTED,"not implemented");
    }
    if (verb > 0) {
      int nb_added_dofs = 0;
      int nb_inactive_added_dofs = 0;
      for(std::map<EntityType,int>::iterator it = dof_counter.begin();it!=dof_counter.end();it++) {
        switch (it->first) {
          case MBVERTEX:
          PetscSynchronizedPrintf(
            comm,
            "nb added dofs (vertices) %d (inactive %d)\n",
            it->second,inactive_dof_counter[it->first]
          );
          break;
          case MBEDGE:
          PetscSynchronizedPrintf(
            comm,
            "nb added dofs (edges) %d (inactive %d)\n",
            it->second,inactive_dof_counter[it->first]
          );
          break;
          case MBTRI:
          PetscSynchronizedPrintf(
            comm,
            "nb added dofs (triangles) %d (inactive %d)\n",
            it->second,inactive_dof_counter[it->first]
          );
          break;
          case MBQUAD:
          PetscSynchronizedPrintf(
            comm,
            "nb added dofs (quads) %d (inactive %d)\n",
            it->second,inactive_dof_counter[it->first]
          );
          break;
          case MBTET:
          PetscSynchronizedPrintf(
            comm,
            "nb added dofs (tets) %d (inactive %d)\n",
            it->second,inactive_dof_counter[it->first]
          );
          break;
          case MBPRISM:
          PetscSynchronizedPrintf(
            comm,
            "nb added dofs (prisms) %d (inactive %d)\n",
            it->second,inactive_dof_counter[it->first]
          );
          break;
          case MBENTITYSET:
          PetscSynchronizedPrintf(
            comm,
            "nb added dofs (meshsets) %d (inactive %d)\n",
            it->second,inactive_dof_counter[it->first]
          );
          break;
          default:
          SETERRQ(comm,MOFEM_NOT_IMPLEMENTED,"not implemented");
        }
        nb_added_dofs += it->second;
        nb_inactive_added_dofs += inactive_dof_counter[it->first];
      }
      if(verbose>0) {
        PetscSynchronizedPrintf(
          comm,
          "nb added dofs %d (number of inactive dofs %d)\n",
          nb_added_dofs,nb_inactive_added_dofs
        );
      }
    }
  }
  *buildMoFEM = 1<<0;
  if(verb>0) {
    PetscSynchronizedPrintf(comm,"Nb. dofs %u\n",dofsField.size());
  }
  PetscSynchronizedFlush(comm,PETSC_STDOUT);
  PetscFunctionReturn(0);
  //return 0;
}
PetscErrorCode Core::list_dofs_by_field_name(const std::string &field_name) const {
  PetscFunctionBegin;
  DofEntityByFieldName::iterator dit,hi_dit;
  dit = dofsField.get<FieldName_mi_tag>().lower_bound(field_name);
  hi_dit = dofsField.get<FieldName_mi_tag>().upper_bound(field_name);
  for(;dit!=hi_dit;dit++) {
    std::ostringstream ss;
    ss << "rank " << rAnk << " ";
    ss << *dit << std::endl;
    PetscSynchronizedPrintf(comm,ss.str().c_str());
  }
  PetscSynchronizedFlush(comm,PETSC_STDOUT);
  PetscFunctionReturn(0);
}
PetscErrorCode Core::list_fields() const {
  PetscFunctionBegin;
  typedef Field_multiIndex::index<BitFieldId_mi_tag>::type FieldSetById;
  const FieldSetById &set_id = fIelds.get<BitFieldId_mi_tag>();
  FieldSetById::iterator miit = set_id.begin();
  for(;miit!=set_id.end();miit++) {
    std::ostringstream ss;
    ss << *miit << std::endl;
    PetscSynchronizedPrintf(comm,ss.str().c_str());
  }
  PetscSynchronizedFlush(comm,PETSC_STDOUT);
  PetscFunctionReturn(0);
}


PetscErrorCode Core::list_adjacencies() const {
  PetscFunctionBegin;
  MoFEMEntityEntFiniteElementAdjacencyMap_multiIndex::iterator miit = entFEAdjacencies.begin();
  for(;miit!=entFEAdjacencies.end();miit++) {
    std::ostringstream ss;
    ss << *miit << std::endl;
    PetscSynchronizedPrintf(comm,ss.str().c_str());
  }
  PetscSynchronizedFlush(comm,PETSC_STDOUT);
  PetscFunctionReturn(0);
}

PetscErrorCode Core::shift_left_bit_ref(const int shift,int verb) {
  PetscFunctionBegin;
  SETERRQ(comm,MOFEM_NOT_IMPLEMENTED,"not implemented");
  PetscFunctionReturn(0);
}
PetscErrorCode Core::shift_right_bit_ref(const int shift,int verb) {
  PetscFunctionBegin;
  if(verb==-1) verb = verbose;
  BitRefLevel delete_bits;
  for(int ii = 0;ii<shift;ii++) {
    delete_bits.set(0);
    ierr = delete_ents_by_bit_ref(delete_bits,delete_bits,verb); CHKERRQ(ierr);
  }
  RefEntity_multiIndex::iterator ent_it = refinedEntities.begin();
  for(;ent_it!=refinedEntities.end();ent_it++) {
    if(verb>5) {
      std::cout << (*ent_it)->getBitRefLevel() << " : ";
    }
    bool success = refinedEntities.modify(ent_it,RefEntity_change_right_shift(shift));
    if(!success) SETERRQ(comm,MOFEM_DATA_INCONSISTENCY,"inconsistency in data");
    if(verb>5) {
      std::cout << (*ent_it)->getBitRefLevel() << std::endl;
    }
  }
  PetscFunctionReturn(0);
}
PetscErrorCode Core::get_entities_by_type_and_ref_level(const BitRefLevel &bit,const BitRefLevel &mask,const EntityType type,const EntityHandle meshset,int verb) {
  PetscFunctionBegin;
  Range ents;
  ierr = get_entities_by_type_and_ref_level(bit,mask,type,ents,verb); CHKERRQ(ierr);
  rval = moab.add_entities(meshset,ents); CHKERRQ_MOAB(rval);
  PetscFunctionReturn(0);
}
PetscErrorCode Core::get_entities_by_type_and_ref_level(const BitRefLevel &bit,const BitRefLevel &mask,const EntityType type,Range &ents,int verb) {
  PetscFunctionBegin;
  if(verb==-1) verb = verbose;
  ierr = moab.get_entities_by_type(0,type,ents,false); CHKERRQ(ierr);
  Range::iterator eit = ents.begin();
  for(;eit!=ents.end();) {
    BitRefLevel bit2;
    rval = moab.tag_get_data(th_RefBitLevel,&*eit,1,&bit2); CHKERRQ_MOAB(rval);
    if(mask.any()&&bit2.none()) {
      eit = ents.erase(eit);
      continue;
    }
    // Not masked
    if((bit2&mask) != bit2) {
      eit = ents.erase(eit);
      continue;
    }
    // Not in bit
    if((bit2&bit).none()) {
      eit = ents.erase(eit);
      continue;
    }
    eit++;
  }
  PetscFunctionReturn(0);
}
PetscErrorCode Core::get_entities_by_ref_level(const BitRefLevel &bit,const BitRefLevel &mask,const EntityHandle meshset) {
  PetscFunctionBegin;
  Range ents;
  ierr = get_entities_by_ref_level(bit,mask,ents); CHKERRQ(ierr);
  rval = moab.add_entities(meshset,ents); CHKERRQ_MOAB(rval);
  PetscFunctionReturn(0);
}
PetscErrorCode Core::get_entities_by_ref_level(const BitRefLevel &bit,const BitRefLevel &mask,Range &ents) {
  PetscFunctionBegin;
  Range meshset_ents;
  ierr = moab.get_entities_by_type(0,MBENTITYSET,meshset_ents,false); CHKERRQ(ierr);
  ierr = moab.get_entities_by_handle(0,ents,false); CHKERRQ(ierr);
  ents.merge(meshset_ents);
  Range::iterator eit = ents.begin();
  for(;eit!=ents.end();) {
    switch (moab.type_from_handle(*eit)) {
      case MBVERTEX:
      case MBEDGE:
      case MBTRI:
      case MBQUAD:
      case MBTET:
      case MBPRISM:
      break;
      case MBENTITYSET:
      break;
      default:
      eit = ents.erase(eit);
      continue;
    }
    BitRefLevel bit2;
    rval = moab.tag_get_data(th_RefBitLevel,&*eit,1,&bit2); CHKERRQ_MOAB(rval);
    if(mask.any()&&bit2.none()) {
      eit = ents.erase(eit);
      continue;
    }
    if((bit2&mask) != bit2) {
      eit = ents.erase(eit);
      continue;
    }
    if((bit2&bit).none()) {
      eit = ents.erase(eit);
      continue;
    }
    eit++;
  }
  PetscFunctionReturn(0);
}
PetscErrorCode Core::add_ref_level_to_entities(const BitRefLevel &bit,Range &ents) {
  PetscFunctionBegin;
  Range::iterator eit = ents.begin();
  for(;eit!=ents.end();eit++) {
    BitRefLevel bit2;
    rval = moab.tag_get_data(th_RefBitLevel,&*eit,1,&bit2); CHKERRQ_MOAB(rval);
    bit2 |= bit;
    rval = moab.tag_set_data(th_RefBitLevel,&*eit,1,&bit2); CHKERRQ_MOAB(rval);
  }
  PetscFunctionReturn(0);
}
PetscErrorCode Core::set_ref_level_to_entities(const BitRefLevel &bit,Range &ents) {
  PetscFunctionBegin;
  Range::iterator eit = ents.begin();
  for(;eit!=ents.end();eit++) {
    BitRefLevel bit2;
    rval = moab.tag_get_data(th_RefBitLevel,&*eit,1,&bit2); CHKERRQ_MOAB(rval);
    bit2 = bit;
    rval = moab.tag_set_data(th_RefBitLevel,&*eit,1,&bit2); CHKERRQ_MOAB(rval);
  }
  PetscFunctionReturn(0);
}
PetscErrorCode Core::update_meshset_by_entities_children(
    const EntityHandle parent, const BitRefLevel &child_bit,const EntityHandle child, EntityType child_type,const bool recursive,int verb) {
  PetscFunctionBegin;
  if(verb==-1) verb = verbose;

  Range ents;
  rval = moab.get_entities_by_handle(parent,ents,recursive);
  if(rval != MB_SUCCESS) {
    std::cerr << parent << std::endl;
    std::cerr << moab.type_from_handle(parent) <<  " " << MBENTITYSET << std::endl;
  } CHKERRQ_MOAB(rval);

  typedef RefEntity_multiIndex::index<Composite_ParentEnt_And_EntType_mi_tag>::type RefEntsByComposite;
  RefEntsByComposite &ref_ents = refinedEntities.get<Composite_ParentEnt_And_EntType_mi_tag>();
  Range::iterator eit = ents.begin();
  for(;eit!=ents.end();eit++) {
    if(verb>2) {
      std::ostringstream ss;
      ss << "ent " << *eit << std::endl;;
      PetscPrintf(comm,ss.str().c_str());
    }
    RefEntsByComposite::iterator miit = ref_ents.lower_bound(boost::make_tuple(*eit,child_type));
    RefEntsByComposite::iterator hi_miit = ref_ents.upper_bound(boost::make_tuple(*eit,child_type));
    for(;miit!=hi_miit;miit++) {
      if(verb>2) {
        std::ostringstream ss;
        ss << "any bit " << *miit << std::endl;;
        PetscPrintf(comm,ss.str().c_str());
      }
      if(((*miit)->getBitRefLevel()&child_bit).any()) {
        EntityHandle ref_ent = (*miit)->getRefEnt();
        if(ref_ent == *eit) continue;
        if(ref_ent == 0) {
          SETERRQ(comm,MOFEM_IMPOSIBLE_CASE,"this should not happen");
        }
        if(moab.type_from_handle(*eit)==MBENTITYSET) {
          SETERRQ(comm,MOFEM_IMPOSIBLE_CASE,"this should not happen");
        }
        rval = moab.add_entities(child,&ref_ent,1); CHKERRQ_MOAB(rval);
        if(verb>1) {
          std::ostringstream ss;
          ss << "good bit " << *miit << std::endl;
          PetscPrintf(comm,ss.str().c_str());
        }
      }
    }
  }
  PetscFunctionReturn(0);
}
PetscErrorCode Core::update_field_meshset_by_entities_children(const BitRefLevel &child_bit,int verb) {
  PetscFunctionBegin;
  Field_multiIndex::iterator fit = fIelds.begin();
  for(;fit!=fIelds.end();fit++) {
    EntityHandle meshset = (*fit)->getMeshset();
    ierr = update_meshset_by_entities_children(meshset,child_bit,meshset,MBTET,false,verb);  CHKERRQ(ierr);
    ierr = update_meshset_by_entities_children(meshset,child_bit,meshset,MBTRI,false,verb);  CHKERRQ(ierr);
    ierr = update_meshset_by_entities_children(meshset,child_bit,meshset,MBEDGE,false,verb);  CHKERRQ(ierr);
    ierr = update_meshset_by_entities_children(meshset,child_bit,meshset,MBVERTEX,false,verb);  CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}
PetscErrorCode Core::update_field_meshset_by_entities_children(const std::string name,const BitRefLevel &child_bit,int verb) {
  PetscFunctionBegin;
  Field_multiIndex::index<FieldName_mi_tag>::type::iterator fit;
  fit = fIelds.get<FieldName_mi_tag>().find(name);
  EntityHandle meshset = (*fit)->getMeshset();
  ierr = update_meshset_by_entities_children(meshset,child_bit,meshset,MBTET,false,verb);  CHKERRQ(ierr);
  ierr = update_meshset_by_entities_children(meshset,child_bit,meshset,MBTRI,false,verb);  CHKERRQ(ierr);
  ierr = update_meshset_by_entities_children(meshset,child_bit,meshset,MBEDGE,false,verb);  CHKERRQ(ierr);
  ierr = update_meshset_by_entities_children(meshset,child_bit,meshset,MBVERTEX,false,verb);  CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
PetscErrorCode Core::update_finite_element_meshset_by_entities_children(const std::string name,const BitRefLevel &child_bit,const EntityType fe_ent_type,int verb) {
  PetscFunctionBegin;
  typedef FiniteElement_multiIndex::index<FiniteElement_name_mi_tag>::type FiniteElementsByName;
  const FiniteElementsByName& set = finiteElements.get<FiniteElement_name_mi_tag>();
  FiniteElementsByName::iterator miit = set.find(name);
  if(miit==set.end()) THROW_MESSAGE(("finite element < "+name+" > not found (top tip: check spelling)").c_str());
  EntityHandle meshset = (*miit)->getMeshset();
  ierr = update_meshset_by_entities_children(meshset,child_bit,meshset,fe_ent_type,false,verb);  CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
PetscErrorCode Core::get_problem_finite_elements_entities(const std::string &problem_name,const std::string &fe_name,const EntityHandle meshset) {
  PetscFunctionBegin;
  typedef MoFEMProblem_multiIndex::index<Problem_mi_tag>::type ProblemsByName;
  ProblemsByName &pRoblems_set = pRoblems.get<Problem_mi_tag>();
  ProblemsByName::iterator p_miit = pRoblems_set.find(problem_name);
  if(p_miit == pRoblems_set.end()) SETERRQ1(comm,1,"no such problem like < %s >",problem_name.c_str());
  NumeredEntFiniteElement_multiIndex &numeredFiniteElements = const_cast<NumeredEntFiniteElement_multiIndex&>(p_miit->numeredFiniteElements);
  NumeredEntFiniteElementbyName::iterator miit = numeredFiniteElements.get<FiniteElement_name_mi_tag>().lower_bound(fe_name);
  for(;miit!=numeredFiniteElements.get<FiniteElement_name_mi_tag>().upper_bound(fe_name);miit++) {
    EntityHandle ent = (*miit)->getEnt();
    rval = moab.add_entities(meshset,&ent,1); CHKERRQ_MOAB(rval);
    int part = (*miit)->getPart();
    rval = moab.tag_set_data(th_Part,&ent,1,&part); CHKERRQ_MOAB(rval);
  }
  PetscFunctionReturn(0);
}

MoFEMEntityByFieldName::iterator Core::get_ent_moabfield_by_name_begin(const std::string &field_name) const {
  return entsFields.get<FieldName_mi_tag>().lower_bound(field_name);
}
MoFEMEntityByFieldName::iterator Core::get_ent_moabfield_by_name_end(const std::string &field_name) const {
  return entsFields.get<FieldName_mi_tag>().upper_bound(field_name);
}
DofEntityByFieldName::iterator Core::get_dofs_by_name_begin(const std::string &field_name) const {
  return dofsField.get<FieldName_mi_tag>().lower_bound(field_name);
}
DofEntityByFieldName::iterator Core::get_dofs_by_name_end(const std::string &field_name) const {
  return dofsField.get<FieldName_mi_tag>().upper_bound(field_name);
}
DofEntityByNameAndEnt::iterator
Core::get_dofs_by_name_and_ent_begin(const std::string &field_name,const EntityHandle ent) const {
  return dofsField.get<Composite_Name_And_Ent_mi_tag>().lower_bound(boost::make_tuple(field_name,ent));
}
DofEntityByNameAndEnt::iterator
Core::get_dofs_by_name_and_ent_end(const std::string &field_name,const EntityHandle ent) const {
  return dofsField.get<Composite_Name_And_Ent_mi_tag>().upper_bound(boost::make_tuple(field_name,ent));
}
DofEntityByNameAndType::iterator
Core::get_dofs_by_name_and_type_begin(const std::string &field_name,const EntityType type) const {
  return dofsField.get<Composite_Name_And_Type_mi_tag>().lower_bound(boost::make_tuple(field_name,type));
}
DofEntityByNameAndType::iterator
Core::get_dofs_by_name_and_type_end(const std::string &field_name,const EntityType type) const {
  return dofsField.get<Composite_Name_And_Type_mi_tag>().upper_bound(boost::make_tuple(field_name,type));
}
PetscErrorCode Core::check_number_of_ents_in_ents_field(const std::string& name) const {
  PetscFunctionBegin;
  Field_multiIndex::index<FieldName_mi_tag>::type::iterator it = fIelds.get<FieldName_mi_tag>().find(name);
  if(it == fIelds.get<FieldName_mi_tag>().end()) {
    SETERRQ1(comm,1,"field not found < %s >",name.c_str());
  }
  EntityHandle meshset = (*it)->getMeshset();

  int num_entities;
  MoABErrorCode rval;
  rval = moab.get_number_entities_by_handle(meshset,num_entities); CHKERRQ_MOAB(rval);
  if(entsFields.get<FieldName_mi_tag>().count((*it)->getName())
    != (unsigned int)num_entities) {
    SETERRQ1(comm,1,"not equal number of entities in meshset and field multiindex < %s >",name.c_str());
  }
  PetscFunctionReturn(0);
}
PetscErrorCode Core::check_number_of_ents_in_ents_field() const {
  PetscFunctionBegin;
  Field_multiIndex::index<FieldName_mi_tag>::type::iterator it = fIelds.get<FieldName_mi_tag>().begin();
  for(;it!=fIelds.get<FieldName_mi_tag>().end();it++) {
    if((*it)->getSpace() == NOFIELD) continue; //FIXME: should be treated properly, not test is just skipped for this NOFIELD space
    EntityHandle meshset = (*it)->getMeshset();
    MoABErrorCode rval;
    int num_entities;
    rval = moab.get_number_entities_by_handle(meshset,num_entities); CHKERRQ_MOAB(rval);
    if(entsFields.get<FieldName_mi_tag>().count((*it)->getName()) != (unsigned int)num_entities) {
      SETERRQ1(comm,1,"not equal number of entities in meshset and field multiindex < %s >",(*it)->getName().c_str());
    }
  }
  PetscFunctionReturn(0);
}


}
