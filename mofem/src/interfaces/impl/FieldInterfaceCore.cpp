/** \file FieldInterfaceCore.cpp
 * \brief Mylti-index containers, data structures and other low-level functions
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
#include <CubitBCData.hpp>
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
#include <FieldInterface.hpp>
#include <MeshRefinment.hpp>
#include <PrismInterface.hpp>
#include <SeriesRecorder.hpp>
#include <Core.hpp>

namespace MoFEM {

const static int debug = 1;

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
      SETERRQ1(PETSC_COMM_SELF,MOFEM_OPERATION_UNSUCCESSFUL,"field is <%s> in database",name.c_str());
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
      CoordSys_multiIndex::index<CoordSysName_mi_tag>::type::iterator undefined_cs_it;
      undefined_cs_it = coordinateSystems.get<CoordSysName_mi_tag>().find("UNDEFINED");
      if(undefined_cs_it==coordinateSystems.get<CoordSysName_mi_tag>().end()) {
        SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"Undefined system not found");
      }
      EntityHandle coord_sys_id = (*undefined_cs_it)->getMeshSet();
      rval = moab.tag_set_data(th_CoordSysMeshSet,&meshset,1,&coord_sys_id); CHKERRQ_MOAB(rval);
      p = fIelds.insert(boost::shared_ptr<Field>(new Field(moab,meshset,*undefined_cs_it)));
      if(bh == MF_EXCL) {
        if(!p.second) SETERRQ1(
          PETSC_COMM_SELF,1,
          "field not inserted %s (top tip, it could be already there)",
          Field(moab,meshset,*undefined_cs_it).getName().c_str()
        );
      }
    } catch (MoFEMException const &e) {
      SETERRQ(PETSC_COMM_SELF,e.errorCode,e.errorMessage);
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
    SETERRQ(PETSC_COMM_SELF,e.errorCode,e.errorMessage);
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
      rval = moab.get_adjacencies(edges,0,false,nodes,Interface::UNION); CHKERRQ_MOAB(rval);
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
      SETERRQ(PETSC_COMM_SELF,MOFEM_NOT_IMPLEMENTED,"sorry, not implemented");
      break;
    case HDIV:
      SETERRQ(PETSC_COMM_SELF,MOFEM_NOT_IMPLEMENTED,"sorry, not implemented, HDIV not implemented for EDGEs");
      break;
    default:
      SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"add_ents_to_field_by_EDGEs this field not work for EDGEs");
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
    SETERRQ(PETSC_COMM_SELF,e.errorCode,e.errorMessage);
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
    SETERRQ(PETSC_COMM_SELF,e.errorCode,e.errorMessage);
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
    SETERRQ(PETSC_COMM_SELF,e.errorCode,e.errorMessage);
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
    SETERRQ(PETSC_COMM_SELF,e.errorCode,e.errorMessage);
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
    rval = moab.get_adjacencies(tris,0,false,nodes,Interface::UNION); CHKERRQ_MOAB(rval);
    {
      Range topo_nodes;
      rval = moab.get_connectivity(tris,topo_nodes,true); CHKERRQ_MOAB(rval);
      Range mid_nodes;
      rval = moab.get_connectivity(tris,mid_nodes,false); CHKERRQ_MOAB(rval);
      mid_nodes = subtract(mid_nodes,topo_nodes);
      nodes = subtract(nodes,mid_nodes);
    }
    rval = moab.add_entities(idm,nodes); CHKERRQ_MOAB(rval);
    rval = moab.get_adjacencies(tris,1,false,edges,Interface::UNION); CHKERRQ_MOAB(rval);
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
    SETERRQ(PETSC_COMM_SELF,MOFEM_NOT_IMPLEMENTED,"sorry, not implemented, HCURL not implemented for TRI");
    break;
    case HDIV:
    SETERRQ(PETSC_COMM_SELF,MOFEM_NOT_IMPLEMENTED,"sorry, not implemented, HDIV not implemented for TRI");
    break;
    default:
    SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"add_ents_to_field_by_TRIs this field not work for TRIs");
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
    SETERRQ(PETSC_COMM_SELF,e.errorCode,e.errorMessage);
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
    SETERRQ(PETSC_COMM_SELF,e.errorCode,e.errorMessage);
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
    SETERRQ(PETSC_COMM_SELF,e.errorCode,e.errorMessage);
  }
  FieldSpace space;
  rval = moab.tag_get_data(th_FieldSpace,&idm,1,&space); CHKERRQ_MOAB(rval);
  switch (space) {
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
    SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"add_ents_to_field_by_TRIs this field not work for TRIs");
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
    SETERRQ(PETSC_COMM_SELF,e.errorCode,e.errorMessage);
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
    SETERRQ(PETSC_COMM_SELF,e.errorCode,e.errorMessage);
  }
  PetscFunctionReturn(0);
}
PetscErrorCode Core::synchronise_field_entities(const BitFieldId id,int verb) {
  PetscFunctionBegin;
  if(verb==-1) verb = verbose;
  EntityHandle idm;
  try {
    idm = get_field_meshset(id);
  } catch (MoFEMException const &e) {
    SETERRQ(PETSC_COMM_SELF,e.errorCode,e.errorMessage);
  }
  Range ents;
  ierr = moab.get_entities_by_handle(idm,ents,false); CHKERRQ(ierr);
  ierr = synchronise_entities(ents,verb); CHKERRQ(ierr);
  rval = moab.add_entities(idm,ents); CHKERRQ_MOAB(rval);
  PetscFunctionReturn(0);
}
PetscErrorCode Core::synchronise_field_entities(const std::string& name,int verb) {
  PetscFunctionBegin;
  if(verb==-1) verb = verbose;
  *buildMoFEM = 0;
  try {
    ierr = synchronise_field_entities(get_BitFieldId(name),verb);  CHKERRQ(ierr);
  } catch (MoFEMException const &e) {
    SETERRQ(PETSC_COMM_SELF,e.errorCode,e.errorMessage);
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
    SETERRQ(PETSC_COMM_SELF,e.errorCode,e.errorMessage);
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
    rval = moab.get_adjacencies(tets,0,false,nodes,Interface::UNION); CHKERRQ_MOAB(rval);
    {
      Range topo_nodes;
      rval = moab.get_connectivity(tets,topo_nodes,true); CHKERRQ_MOAB(rval);
      Range mid_nodes;
      rval = moab.get_connectivity(tets,mid_nodes,false); CHKERRQ_MOAB(rval);
      mid_nodes = subtract(mid_nodes,topo_nodes);
      nodes = subtract(nodes,mid_nodes);
    }
    rval = moab.add_entities(idm,nodes); CHKERRQ_MOAB(rval);
    rval = moab.get_adjacencies(tets,2,false,tris,Interface::UNION); CHKERRQ_MOAB(rval);
    rval = moab.add_entities(idm,tris); CHKERRQ_MOAB(rval);
    rval = moab.get_adjacencies(tets,1,false,edges,Interface::UNION); CHKERRQ_MOAB(rval);
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
    rval = moab.get_adjacencies(tets,2,false,tris,Interface::UNION); CHKERRQ_MOAB(rval);
    rval = moab.add_entities(idm,tris); CHKERRQ_MOAB(rval);
    rval = moab.get_adjacencies(tets,1,false,edges,Interface::UNION); CHKERRQ_MOAB(rval);
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
    rval = moab.get_adjacencies(tets,2,false,tris,Interface::UNION); CHKERRQ_MOAB(rval);
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
    SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"add_ents_to_field_by_TETs this field not work for TETs");
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
    SETERRQ(PETSC_COMM_SELF,e.errorCode,e.errorMessage);
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
    SETERRQ(PETSC_COMM_SELF,e.errorCode,e.errorMessage);
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
    SETERRQ(PETSC_COMM_SELF,e.errorCode,e.errorMessage);
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
    rval = moab.get_adjacencies(quads,0,false,nodes,Interface::UNION); CHKERRQ_MOAB(rval);
    {
      Range topo_nodes;
      rval = moab.get_connectivity(quads,topo_nodes,true); CHKERRQ_MOAB(rval);
      Range mid_nodes;
      rval = moab.get_connectivity(quads,mid_nodes,false); CHKERRQ_MOAB(rval);
      mid_nodes = subtract(mid_nodes,topo_nodes);
      nodes = subtract(nodes,mid_nodes);
    }
    rval = moab.add_entities(idm,nodes); CHKERRQ_MOAB(rval);
    rval = moab.get_adjacencies(quads,1,true,edges,Interface::UNION); CHKERRQ_MOAB(rval);
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
    SETERRQ(PETSC_COMM_SELF,MOFEM_NOT_IMPLEMENTED,"not implemented");
    break;
    case HDIV:
    SETERRQ(PETSC_COMM_SELF,MOFEM_NOT_IMPLEMENTED,"not implemented");
    break;
    default:
    SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"add_ents_to_field_by_TETs this field not work for TETs");
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
    SETERRQ(PETSC_COMM_SELF,e.errorCode,e.errorMessage);
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
    SETERRQ(PETSC_COMM_SELF,e.errorCode,e.errorMessage);
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
    SETERRQ(PETSC_COMM_SELF,e.errorCode,e.errorMessage);
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
    rval = moab.get_adjacencies(prisms,0,false,nodes,Interface::UNION); CHKERRQ_MOAB(rval);
    {
      Range topo_nodes;
      rval = moab.get_connectivity(prisms,topo_nodes,true); CHKERRQ_MOAB(rval);
      Range mid_nodes;
      rval = moab.get_connectivity(prisms,mid_nodes,false); CHKERRQ_MOAB(rval);
      mid_nodes = subtract(mid_nodes,topo_nodes);
      nodes = subtract(nodes,mid_nodes);
    }
    rval = moab.add_entities(idm,nodes); CHKERRQ_MOAB(rval);
    rval = moab.get_adjacencies(prisms,2,true,faces,Interface::UNION); CHKERRQ_MOAB(rval);
    rval = moab.add_entities(idm,faces); CHKERRQ_MOAB(rval);
    rval = moab.get_adjacencies(prisms,1,true,edges,Interface::UNION); CHKERRQ_MOAB(rval);
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
    SETERRQ(PETSC_COMM_SELF,MOFEM_NOT_IMPLEMENTED,"not implemented");
    break;
    case HDIV:
    SETERRQ(PETSC_COMM_SELF,MOFEM_NOT_IMPLEMENTED,"not implemented");
    break;
    default:
    SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"add_ents_to_field_by_TETs this field not work for TETs");
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
    SETERRQ(PETSC_COMM_SELF,e.errorCode,e.errorMessage);
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
    SETERRQ(PETSC_COMM_SELF,e.errorCode,e.errorMessage);
  }
  PetscFunctionReturn(0);
}
PetscErrorCode Core::set_field_order(const Range &ents,const BitFieldId id,const ApproximationOrder order,int verb) {
  PetscFunctionBegin;
  if(verb==-1) verb = verbose;
  *buildMoFEM = 0;

  //check field & meshset
  typedef Field_multiIndex::index<BitFieldId_mi_tag>::type field_set_by_id;
  const field_set_by_id &set_id = fIelds.get<BitFieldId_mi_tag>();
  field_set_by_id::iterator miit = set_id.find(id);
  if(miit==set_id.end()) SETERRQ(PETSC_COMM_SELF,MOFEM_NOT_FOUND,"no filed found");
  EntityHandle idm;
  try {
   idm = get_field_meshset(id);
  } catch (MoFEMException const &e) {
    SETERRQ(PETSC_COMM_SELF,e.errorCode,e.errorMessage);
  }

  //intersection with field meshset
  Range ents_of_id_meshset;
  rval = moab.get_entities_by_handle(idm,ents_of_id_meshset,false); CHKERRQ_MOAB(rval);
  Range ents_ = intersect(ents,ents_of_id_meshset);
  if(verb>1) {
    PetscSynchronizedPrintf(comm,"nb. of ents for order change in the field %d\n",ents_.size());
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
    PetscSynchronizedPrintf(comm,"nb. of ents in the multi index field %d\n",ents_id_view.size());
  }


  // get tags on entities
  std::vector<ApproximationOrder*> tag_data_order(ents_.size());
  rval = moab.tag_get_by_ptr((*miit)->th_AppOrder,ents_,(const void **)&tag_data_order[0]); CHKERRQ_MOAB(rval);

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
          SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"approximation order for H1 space and vertex different than 1 makes not sense");
        }
      }
      break;
      case HDIV:
      if(moab.type_from_handle(*eit)==MBVERTEX) {
        SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"HDIV space on vertices makes no sense");
      }
      if(moab.type_from_handle(*eit)==MBEDGE) {
        SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"HDIV space on edges makes no sense");
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
      MoFEMEntity_multiIndex::iterator miit = entsFields.get<Unique_mi_tag>().find((*vit)->getGlobalUniqueId());
      if((*miit)->getMaxOrder()<order) nb_ents_set_order_up++;
      if((*miit)->getMaxOrder()>order) nb_ents_set_order_down++;

      {

      	//set dofs inactive if order is reduced, and set new order to entity if
      	//order is increased (note that dofs are not build if order is
      	//increased)

        typedef DofEntity_multiIndex::index<Composite_Name_And_Ent_mi_tag>::type dof_set_type;
        dof_set_type& set_set = dofsField.get<Composite_Name_And_Ent_mi_tag>();
        dof_set_type::iterator dit = set_set.lower_bound(boost::make_tuple((*miit)->getNameRef(),(*miit)->getEnt()));
        dof_set_type::iterator hi_dit = set_set.upper_bound(boost::make_tuple((*miit)->getNameRef(),(*miit)->getEnt()));

        for(;dit!=hi_dit;dit++) {
          if((*dit)->get_dof_order()<=order) continue;
          bool success = dofsField.modify(dofsField.project<0>(dit),DofEntity_active_change(false));
          if(!success) SETERRQ(PETSC_COMM_SELF,MOFEM_OPERATION_UNSUCCESSFUL,"modification unsuccessful");
        }

        bool success = entsFields.modify(entsFields.project<0>(miit),MoFEMEntity_change_order(order));
        if(!success) SETERRQ(PETSC_COMM_SELF,MOFEM_OPERATION_UNSUCCESSFUL,"modification unsuccessful");

      }

    } else {

      //entity is not in database and order is changed or reset
      *tag_data_order[ee] = order;
      RefEntity_multiIndex::index<Ent_mi_tag>::type::iterator miit_ref_ent = refinedEntities.get<Ent_mi_tag>().find(*eit);
      if(miit_ref_ent==refinedEntities.get<Ent_mi_tag>().end()) {

        RefEntity ref_ent(basicEntityDataPtr,*eit);
        // FIXME: need some consistent policy in that case
        if(ref_ent.getBitRefLevel().none()) continue; // not on any mesh and not in database
        std::cerr << ref_ent << std::endl;
        std::cerr << "bit level " << ref_ent.getBitRefLevel() << std::endl;
        SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"Try to add entities which are not seeded or added to database");

      }

      try {

        // increase order
        boost::shared_ptr<MoFEMEntity> moabent(new MoFEMEntity(*miit,*miit_ref_ent));
        std::pair<MoFEMEntity_multiIndex::iterator,bool> e_miit = entsFields.insert(moabent);
        bool success = entsFields.modify(e_miit.first,MoFEMEntity_change_order(order));
        if(!success) SETERRQ(PETSC_COMM_SELF,MOFEM_OPERATION_UNSUCCESSFUL,"modification unsuccessful");
        nb_ents_set_order_new++;

      } catch (MoFEMException const &e) {
        SETERRQ(PETSC_COMM_SELF,e.errorCode,e.errorMessage);
      } catch (const std::exception& ex) {
        std::ostringstream ss;
        ss << ex.what() << std::endl;
        SETERRQ(PETSC_COMM_SELF,MOFEM_STD_EXCEPTION_THROW,ss.str().c_str());
      }

    }
  }

  if(verb>1) {
    PetscSynchronizedPrintf(comm,"nb. of entities for which order was increased %d (order %d)\n",nb_ents_set_order_up,order);
    PetscSynchronizedPrintf(comm,"nb. of entities for which order was reduced %d (order %d)\n",nb_ents_set_order_down,order);
    PetscSynchronizedPrintf(comm,"nb. of entities for which order set %d (order %d)\n",nb_ents_set_order_new,order);
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
    SETERRQ(PETSC_COMM_SELF,e.errorCode,e.errorMessage);
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
    SETERRQ(PETSC_COMM_SELF,e.errorCode,e.errorMessage);
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
    SETERRQ(PETSC_COMM_SELF,e.errorCode,e.errorMessage);
  }
  PetscFunctionReturn(0);
}
PetscErrorCode Core::set_field_order_by_entity_type_and_bit_ref(const BitRefLevel &bit,const BitRefLevel &mask,const EntityType type,const BitFieldId id,const ApproximationOrder order,int verb) {
  PetscFunctionBegin;
  if(verb==-1) verb = verbose;
  *buildMoFEM = 0;
  Range ents;
  ierr = get_entities_by_type_and_ref_level(bit,mask,type,ents,verb); CHKERRQ(ierr);
  try{
    ierr = set_field_order(ents,id,order,verb); CHKERRQ(ierr);
  } catch (MoFEMException const &e) {
    SETERRQ(PETSC_COMM_SELF,e.errorCode,e.errorMessage);
  }
  PetscFunctionReturn(0);
}
PetscErrorCode Core::set_field_order_by_entity_type_and_bit_ref(const BitRefLevel &bit,const BitRefLevel &mask,const EntityType type,const std::string& name,const ApproximationOrder order,int verb) {
  PetscFunctionBegin;
  if(verb==-1) verb = verbose;
  *buildMoFEM = 0;
  Range ents;
  ierr = get_entities_by_type_and_ref_level(bit,mask,type,ents,verb); CHKERRQ(ierr);
  try{
    ierr = set_field_order(ents,get_BitFieldId(name),order,verb); CHKERRQ(ierr);
  } catch (MoFEMException const &e) {
    SETERRQ(PETSC_COMM_SELF,e.errorCode,e.errorMessage);
  }
  PetscFunctionReturn(0);
}
PetscErrorCode Core::dofs_NoField(const BitFieldId id,std::map<EntityType,int> &dof_counter,int verb) {
  PetscFunctionBegin;
  if(verb==-1) verb = verbose;
  //field it
  typedef Field_multiIndex::index<BitFieldId_mi_tag>::type field_set_by_id;
  const field_set_by_id &set_id = fIelds.get<BitFieldId_mi_tag>();
  //find fiels
  field_set_by_id::iterator miit = set_id.find(id);
  if(miit == set_id.end()) {
    SETERRQ(PETSC_COMM_SELF,MOFEM_NOT_FOUND,"field not found");
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
      SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"database inconsistency");
    }
    std::pair<MoFEMEntity_multiIndex::iterator,bool> e_miit;
    try {
      //create database entity
      boost::shared_ptr<MoFEMEntity> moabent(new MoFEMEntity(*miit,*miit_ref_ent));
      e_miit = entsFields.insert(moabent);
    } catch (MoFEMException const &e) {
      SETERRQ(PETSC_COMM_SELF,e.errorCode,e.errorMessage);
    } catch (const std::exception& ex) {
      std::ostringstream ss;
      ss << ex.what() << std::endl;
      SETERRQ(PETSC_COMM_SELF,MOFEM_STD_EXCEPTION_THROW,ss.str().c_str());
    }
    //this is nor real field in space (set order to zero)
    bool success = entsFields.modify(e_miit.first,MoFEMEntity_change_order(0));
    if(!success) SETERRQ(PETSC_COMM_SELF,MOFEM_OPERATION_UNSUCCESSFUL,"modification unsuccessful");
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
        if(!success) SETERRQ(PETSC_COMM_SELF,MOFEM_OPERATION_UNSUCCESSFUL,"modification unsuccessful");
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
PetscErrorCode Core::dofs_L2H1HcurlHdiv(
  const BitFieldId id,
  std::map<EntityType,int> &dof_counter,
  std::map<EntityType,int> &inactive_dof_counter,
  int verb
) {
  PetscFunctionBegin;
  if(verb==-1) verb = verbose;
  //field it
  typedef Field_multiIndex::index<BitFieldId_mi_tag>::type field_set_by_id;
  //find field
  const field_set_by_id &set_id = fIelds.get<BitFieldId_mi_tag>();
  field_set_by_id::iterator miit = set_id.find(id);
  if(miit == set_id.end()) {
    SETERRQ(PETSC_COMM_SELF,MOFEM_NOT_FOUND,"field not found");
  }
  //ents in the field meshset
  Range ents_of_id_meshset;
  rval = moab.get_entities_by_handle((*miit)->meshSet,ents_of_id_meshset,false); CHKERRQ_MOAB(rval);
  if(verb>5) {
    PetscSynchronizedPrintf(
      comm,"ents in field %s meshset %d\n",(*miit)->getName().c_str(),ents_of_id_meshset.size()
    );
  }
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
      SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"database inconsistency");
    }
    // create mofem entity linked to ref ent
    MoFEMEntity_multiIndex::iterator e_miit;
    try {
      e_miit = entsFields.find(MoFEMEntity(*miit,*miit_ref_ent).getGlobalUniqueId());
    } catch (MoFEMException const &e) {
      SETERRQ(PETSC_COMM_SELF,e.errorCode,e.errorMessage);
    } catch (const std::exception& ex) {
      std::ostringstream ss;
      ss << ex.what() << std::endl;
      SETERRQ(PETSC_COMM_SELF,MOFEM_OPERATION_UNSUCCESSFUL,ss.str().c_str());
    }
    if(e_miit == entsFields.end()) {
      ApproximationOrder order = -1;
      rval = moab.tag_set_data((*miit)->th_AppOrder,&*eit,1,&order); CHKERRQ_MOAB(rval);
      std::pair<MoFEMEntity_multiIndex::iterator,bool> p_e_miit;
      try {
        boost::shared_ptr<MoFEMEntity> moabent(new MoFEMEntity(*miit,*miit_ref_ent));
        p_e_miit = entsFields.insert(moabent);
      } catch (MoFEMException const &e) {
        SETERRQ(PETSC_COMM_SELF,e.errorCode,e.errorMessage);
      } catch (const std::exception& ex) {
        std::ostringstream ss;
        ss << ex.what() << std::endl;
        SETERRQ(PETSC_COMM_SELF,MOFEM_OPERATION_UNSUCCESSFUL,ss.str().c_str());
      }
      if(!p_e_miit.second) SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"data inconsistency");
      bool success = entsFields.modify(p_e_miit.first,MoFEMEntity_change_order(-1));
      if(!success) SETERRQ(PETSC_COMM_SELF,MOFEM_OPERATION_UNSUCCESSFUL,"modification unsuccessful");
      e_miit = p_e_miit.first;
    }
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
          std::pair<DofEntity_multiIndex::iterator,bool> d_miit;
          try {
            boost::shared_ptr<DofEntity> mdof(new DofEntity(*(e_miit),oo,rr,DD));
            d_miit = dofsField.insert(mdof);
            bool is_active;
            if(d_miit.second) {
              if(DD<nb_active_dosf_on_ent) {
                is_active = true;
                dof_counter[(*d_miit.first)->getEntType()]++;
              } else {
                is_active = false;
                inactive_dof_counter[(*d_miit.first)->getEntType()]++;
              }
              bool success = dofsField.modify(d_miit.first,DofEntity_active_change(is_active));
              if(!success) SETERRQ(PETSC_COMM_SELF,MOFEM_OPERATION_UNSUCCESSFUL,"modification unsuccessful");
            } else {
              if(DD<nb_active_dosf_on_ent) {
              } else {
                if((*d_miit.first)->get_active()) {
                  is_active = false;
                  inactive_dof_counter[(*d_miit.first)->getEntType()]++;
                  bool success = dofsField.modify(d_miit.first,DofEntity_active_change(is_active));
                  if(!success) SETERRQ(PETSC_COMM_SELF,MOFEM_OPERATION_UNSUCCESSFUL,"modification unsuccessful");
                }
              }
            }
            //check ent
            if((*d_miit.first)->getEnt()!=(*e_miit)->getEnt()) {
              SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"data inconsistency");
            }
            if((*d_miit.first)->getEntType()!=(*e_miit)->getEntType()) {
              SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"data inconsistency");
            }
            if((*d_miit.first)->getId()!=(*e_miit)->getId()) {
              SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"data inconsistency");
            }
            //check dof
            if((*d_miit.first)->get_dof_order()!=oo) {
              std::ostringstream ss;
              ss << "data inconsistency!" << std::endl;
              ss << "should be " << mdof << std::endl;
              ss << "but is " << *(*d_miit.first) << std::endl;
              SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,ss.str().c_str());
            }
            if((*d_miit.first)->getMaxOrder()!=(*e_miit)->getMaxOrder()) {
              SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"data inconsistency");
            }
          } catch (MoFEMException const &e) {
            SETERRQ(PETSC_COMM_SELF,e.errorCode,e.errorMessage);
          }
        }
      }
    }
    if(DD != (*e_miit)->getNbDofsOnEnt()) {
      SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"data inconsistency");
    }
  }
  PetscFunctionReturn(0);
}
PetscErrorCode Core::build_fields(int verb) {
  PetscFunctionBegin;
  if(verb==-1) verb = verbose;
  typedef Field_multiIndex::index<BitFieldId_mi_tag>::type field_set_by_id;
  field_set_by_id &set_id = fIelds.get<BitFieldId_mi_tag>();
  field_set_by_id::iterator miit = set_id.begin();
  for(;miit!=set_id.end();miit++) {
    std::map<EntityType,int> dof_counter;
    std::map<EntityType,int> inactive_dof_counter;
    if (verb > 0) {
      PetscSynchronizedPrintf(comm,"Build Field %s (rank %d)\n",(*miit)->getName().c_str(),rAnk);
    }
    switch ((*miit)->getSpace()) {
      case NOFIELD:
      ierr = dofs_NoField((*miit)->getId(),dof_counter,verb); CHKERRQ(ierr);
      break;
      case L2:
      case H1:
      case HCURL:
      case HDIV:
      ierr = dofs_L2H1HcurlHdiv((*miit)->getId(),dof_counter,inactive_dof_counter,verb); CHKERRQ(ierr);
      break;
      default:
      SETERRQ(PETSC_COMM_SELF,MOFEM_NOT_IMPLEMENTED,"not implemented");
    }
    if (verb > 0) {
      int nb_added_dofs = 0;
      int nb_inactive_added_dofs = 0;
      for(std::map<EntityType,int>::iterator it = dof_counter.begin();it!=dof_counter.end();it++) {
        switch (it->first) {
          case MBVERTEX:
          PetscSynchronizedPrintf(comm,"nb added dofs (vertices) %d (inactive %d)\n",it->second,inactive_dof_counter[it->first]);
          break;
          case MBEDGE:
          PetscSynchronizedPrintf(comm,"nb added dofs (edges) %d (inactive %d)\n",it->second,inactive_dof_counter[it->first]);
          break;
          case MBTRI:
          PetscSynchronizedPrintf(comm,"nb added dofs (triangles) %d (inactive %d)\n",it->second,inactive_dof_counter[it->first]);
          break;
          case MBQUAD:
          PetscSynchronizedPrintf(comm,"nb added dofs (quads) %d (inactive %d)\n",it->second,inactive_dof_counter[it->first]);
          break;
          case MBTET:
          PetscSynchronizedPrintf(comm,"nb added dofs (tets) %d (inactive %d)\n",it->second,inactive_dof_counter[it->first]);
          break;
          case MBPRISM:
          PetscSynchronizedPrintf(comm,"nb added dofs (prisms) %d (inactive %d)\n",it->second,inactive_dof_counter[it->first]);
          break;
          case MBENTITYSET:
          PetscSynchronizedPrintf(comm,"nb added dofs (meshsets) %d (inactive %d)\n",it->second,inactive_dof_counter[it->first]);
          break;
          default:
          SETERRQ(PETSC_COMM_SELF,MOFEM_NOT_IMPLEMENTED,"not implemented");
        }
        nb_added_dofs += it->second;
        nb_inactive_added_dofs += inactive_dof_counter[it->first];
      }
      PetscSynchronizedPrintf(comm,"nb added dofs %d (number of inactive dofs %d)\n",nb_added_dofs,nb_inactive_added_dofs);
    }
  }
  *buildMoFEM = 1<<0;
  if(verbose>0) {
    PetscSynchronizedPrintf(comm,"Nb. dofs %u\n",dofsField.size());
    PetscSynchronizedFlush(comm,PETSC_STDOUT);
  }
  PetscFunctionReturn(0);
  //return 0;
}
PetscErrorCode Core::list_dofs_by_field_name(const std::string &field_name) const {
  PetscFunctionBegin;
  DofEntity_multiIndex::index<FieldName_mi_tag>::type::iterator dit,hi_dit;
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
  typedef Field_multiIndex::index<BitFieldId_mi_tag>::type field_set_by_id;
  const field_set_by_id &set_id = fIelds.get<BitFieldId_mi_tag>();
  field_set_by_id::iterator miit = set_id.begin();
  for(;miit!=set_id.end();miit++) {
    std::ostringstream ss;
    ss << *miit << std::endl;
    PetscSynchronizedPrintf(comm,ss.str().c_str());
  }
  PetscSynchronizedFlush(comm,PETSC_STDOUT);
  PetscFunctionReturn(0);
}
PetscErrorCode Core::add_finite_element(const std::string &fe_name,enum MoFEMTypes bh) {
  PetscFunctionBegin;
  *buildMoFEM &= 1<<0;
  typedef FiniteElement_multiIndex::index<FiniteElement_name_mi_tag>::type FiniteElements_by_name;
  FiniteElements_by_name &finite_element_name_set = finiteElements.get<FiniteElement_name_mi_tag>();
  FiniteElements_by_name::iterator it_fe = finite_element_name_set.find(fe_name);
  if(bh == MF_EXCL) {
    if(it_fe!=finite_element_name_set.end()) {
      SETERRQ1(PETSC_COMM_SELF,MOFEM_NOT_FOUND,"this < %s > is there",fe_name.c_str());
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
  if(!p.second) SETERRQ(PETSC_COMM_SELF,MOFEM_OPERATION_UNSUCCESSFUL,"FiniteElement not inserted");
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
    SETERRQ(PETSC_COMM_SELF,MOFEM_NOT_FOUND,"this FiniteElement is there");
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
  if(it_fe==finite_element_name_set.end()) SETERRQ(PETSC_COMM_SELF,MOFEM_NOT_FOUND,"this FiniteElement is there");
  try {
    bool success = finite_element_name_set.modify(it_fe,MoFEMFiniteElement_change_bit_add(get_BitFieldId(name_data)));
    if(!success) SETERRQ(PETSC_COMM_SELF,MOFEM_OPERATION_UNSUCCESSFUL,"modification unsuccessful");
  } catch (MoFEMException const &e) {
    SETERRQ(PETSC_COMM_SELF,e.errorCode,e.errorMessage);
  }
  PetscFunctionReturn(0);
}
PetscErrorCode Core::modify_finite_element_add_field_row(const std::string &fe_name,const std::string &name_row) {
  PetscFunctionBegin;
  *buildMoFEM &= 1<<0;
  typedef FiniteElement_multiIndex::index<FiniteElement_name_mi_tag>::type FiniteElements_by_name;
  FiniteElements_by_name &finite_element_name_set = finiteElements.get<FiniteElement_name_mi_tag>();
  FiniteElements_by_name::iterator it_fe = finite_element_name_set.find(fe_name);
  if(it_fe==finite_element_name_set.end()) SETERRQ1(PETSC_COMM_SELF,MOFEM_NOT_FOUND,"this < %s > is not there",fe_name.c_str());
  try {
    bool success = finite_element_name_set.modify(it_fe,MoFEMFiniteElement_row_change_bit_add(get_BitFieldId(name_row)));
    if(!success) SETERRQ(PETSC_COMM_SELF,MOFEM_OPERATION_UNSUCCESSFUL,"modification unsuccessful");
  } catch (MoFEMException const &e) {
    SETERRQ(PETSC_COMM_SELF,e.errorCode,e.errorMessage);
  }
  PetscFunctionReturn(0);
}
PetscErrorCode Core::modify_finite_element_add_field_col(const std::string &fe_name,const std::string &name_col) {
  PetscFunctionBegin;
  *buildMoFEM &= 1<<0;
  typedef FiniteElement_multiIndex::index<FiniteElement_name_mi_tag>::type FiniteElements_by_name;
  FiniteElements_by_name &finite_element_name_set = finiteElements.get<FiniteElement_name_mi_tag>();
  FiniteElements_by_name::iterator it_fe = finite_element_name_set.find(fe_name);
  if(it_fe==finite_element_name_set.end()) SETERRQ(PETSC_COMM_SELF,MOFEM_OPERATION_UNSUCCESSFUL,"this FiniteElement is there");
  try {
    bool success = finite_element_name_set.modify(it_fe,MoFEMFiniteElement_col_change_bit_add(get_BitFieldId(name_col)));
    if(!success) SETERRQ(PETSC_COMM_SELF,MOFEM_OPERATION_UNSUCCESSFUL,"modification unsuccessful");
  } catch (MoFEMException const &e) {
    SETERRQ(PETSC_COMM_SELF,e.errorCode,e.errorMessage);
  }
  PetscFunctionReturn(0);
}
PetscErrorCode Core::modify_finite_element_off_field_data(const std::string &fe_name,const std::string &name_data) {
  PetscFunctionBegin;
  *buildMoFEM &= 1<<0;
  typedef FiniteElement_multiIndex::index<FiniteElement_name_mi_tag>::type FiniteElements_by_name;
  FiniteElements_by_name &finite_element_name_set = finiteElements.get<FiniteElement_name_mi_tag>();
  FiniteElements_by_name::iterator it_fe = finite_element_name_set.find(fe_name);
  if(it_fe==finite_element_name_set.end()) SETERRQ(PETSC_COMM_SELF,MOFEM_NOT_FOUND,"this FiniteElement is there");
  try {
    bool success = finite_element_name_set.modify(it_fe,MoFEMFiniteElement_change_bit_off(get_BitFieldId(name_data)));
    if(!success) SETERRQ(PETSC_COMM_SELF,MOFEM_OPERATION_UNSUCCESSFUL,"modification unsuccessful");
  } catch (MoFEMException const &e) {
    SETERRQ(PETSC_COMM_SELF,e.errorCode,e.errorMessage);
  }
  PetscFunctionReturn(0);
}
PetscErrorCode Core::modify_finite_element_off_field_row(const std::string &fe_name,const std::string &name_row) {
  PetscFunctionBegin;
  *buildMoFEM &= 1<<0;
  typedef FiniteElement_multiIndex::index<FiniteElement_name_mi_tag>::type FiniteElements_by_name;
  FiniteElements_by_name &finite_element_name_set = finiteElements.get<FiniteElement_name_mi_tag>();
  FiniteElements_by_name::iterator it_fe = finite_element_name_set.find(fe_name);
  if(it_fe==finite_element_name_set.end()) SETERRQ1(PETSC_COMM_SELF,MOFEM_NOT_FOUND,"this < %s > is not there",fe_name.c_str());
  try {
    bool success = finite_element_name_set.modify(it_fe,MoFEMFiniteElement_row_change_bit_off(get_BitFieldId(name_row)));
    if(!success) SETERRQ(PETSC_COMM_SELF,MOFEM_OPERATION_UNSUCCESSFUL,"modification unsuccessful");
  } catch (MoFEMException const &e) {
    SETERRQ(PETSC_COMM_SELF,e.errorCode,e.errorMessage);
  }
  PetscFunctionReturn(0);
}
PetscErrorCode Core::modify_finite_element_off_field_col(const std::string &fe_name,const std::string &name_col) {
  PetscFunctionBegin;
  *buildMoFEM &= 1<<0;
  typedef FiniteElement_multiIndex::index<FiniteElement_name_mi_tag>::type FiniteElements_by_name;
  FiniteElements_by_name &finite_element_name_set = finiteElements.get<FiniteElement_name_mi_tag>();
  FiniteElements_by_name::iterator it_fe = finite_element_name_set.find(fe_name);
  if(it_fe==finite_element_name_set.end()) SETERRQ(PETSC_COMM_SELF,MOFEM_NOT_FOUND,"this FiniteElement is there");
  try {
    bool success = finite_element_name_set.modify(it_fe,MoFEMFiniteElement_col_change_bit_off(get_BitFieldId(name_col)));
    if(!success) SETERRQ(PETSC_COMM_SELF,MOFEM_OPERATION_UNSUCCESSFUL,"modification unsuccessful");
  } catch (MoFEMException const &e) {
    SETERRQ(PETSC_COMM_SELF,e.errorCode,e.errorMessage);
  }
  PetscFunctionReturn(0);
}
BitFEId Core::get_BitFEId(const std::string& name) const {
  typedef FiniteElement_multiIndex::index<FiniteElement_name_mi_tag>::type FiniteElements_by_name;
  const FiniteElements_by_name& set = finiteElements.get<FiniteElement_name_mi_tag>();
  FiniteElements_by_name::iterator miit = set.find(name);
  if(miit==set.end()) THROW_MESSAGE(("finite element < "+name+" > not found (top tip: check spelling)").c_str());
  return (*miit)->getId();
}
std::string Core::get_BitFEId_name(const BitFEId id) const {
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
  return get_finite_element_meshset(get_BitFEId(name));
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
PetscErrorCode Core::add_problem(const BitProblemId id,const std::string& name) {
  PetscFunctionBegin;
  EntityHandle meshset;
  rval = moab.create_meshset(MESHSET_SET|MESHSET_TRACK_OWNER,meshset); CHKERRQ_MOAB(rval);
  rval = moab.tag_set_data(th_ProblemId,&meshset,1,&id); CHKERRQ_MOAB(rval);
  void const* tag_data[] = { name.c_str() };
  int tag_sizes[1]; tag_sizes[0] = name.size();
  rval = moab.tag_set_by_ptr(th_ProblemName,&meshset,1,tag_data,tag_sizes); CHKERRQ_MOAB(rval);
  //create entry
  std::pair<MoFEMProblem_multiIndex::iterator,bool> p = pRoblems.insert(MoFEMProblem(moab,meshset));
  NOT_USED(p);
  assert(p.second);
  if(verbose>0) {
    std::ostringstream ss;
    ss << "add problem: " << name << std::endl;
    PetscPrintf(comm,ss.str().c_str());
  }
  PetscFunctionReturn(0);
}
PetscErrorCode Core::add_problem(const std::string& name,enum MoFEMTypes bh,int verb) {
  PetscFunctionBegin;
  if(verb==-1) verb = verbose;
  typedef MoFEMProblem_multiIndex::index<Problem_mi_tag>::type mofem_problems_by_name;
  const mofem_problems_by_name& set = pRoblems.get<Problem_mi_tag>();
  mofem_problems_by_name::iterator miit = set.find(name);
  if(miit==set.end()) {
    BitProblemId id = getProblemShift();
    ierr = add_problem(id,name); CHKERRQ(ierr);
  } else if(bh == MF_EXCL) {
    SETERRQ1(PETSC_COMM_SELF,MOFEM_NOT_FOUND,"problem is in database %s",name.c_str());
  }
  PetscFunctionReturn(0);
}
PetscErrorCode Core::delete_problem(const std::string name) {
  PetscFunctionBegin;
  typedef MoFEMProblem_multiIndex::index<Problem_mi_tag>::type mofem_problems_by_name;
  mofem_problems_by_name &mofem_problems_set = pRoblems.get<Problem_mi_tag>();
  mofem_problems_by_name::iterator p_miit = mofem_problems_set.find(name);
  if(p_miit == mofem_problems_set.end()) {
    SETERRQ1(PETSC_COMM_SELF,1,"no such problem like < %s >",name.c_str());
  }
  EntityHandle meshset = p_miit->meshset;
  mofem_problems_set.erase(p_miit);
  rval = moab.delete_entities(&meshset,1); CHKERRQ_MOAB(rval);
  PetscFunctionReturn(0);
}
BitProblemId Core::get_BitProblemId(const std::string& name) const {
  typedef MoFEMProblem_multiIndex::index<Problem_mi_tag>::type mofem_problems_by_name;
  const mofem_problems_by_name& set = pRoblems.get<Problem_mi_tag>();
  mofem_problems_by_name::iterator miit = set.find(name);
  return miit->getId();
}
PetscErrorCode Core::list_problem() const {
  PetscFunctionBegin;
  typedef MoFEMProblem_multiIndex::index<BitProblemId_mi_tag>::type problem_set_by_id;
  const problem_set_by_id &set_id = pRoblems.get<BitProblemId_mi_tag>();
  problem_set_by_id::iterator miit = set_id.begin();
  for(;miit!=set_id.end();miit++) {
    std::ostringstream ss;
    ss << *miit << std::endl;
    PetscPrintf(comm,ss.str().c_str());
  }
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
    ierr = add_ents_to_finite_element_by_EDGEs(edges,get_BitFEId(name));  CHKERRQ(ierr);
  } catch (MoFEMException const &e) {
    SETERRQ(PETSC_COMM_SELF,e.errorCode,e.errorMessage);
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
    ierr = add_ents_to_finite_element_by_VERTICEs(vert,get_BitFEId(name));  CHKERRQ(ierr);
  } catch (MoFEMException const &e) {
    SETERRQ(PETSC_COMM_SELF,e.errorCode,e.errorMessage);
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
    ierr = add_ents_to_finite_element_by_TRIs(tris,get_BitFEId(name));  CHKERRQ(ierr);
  } catch (MoFEMException const &e) {
    SETERRQ(PETSC_COMM_SELF,e.errorCode,e.errorMessage);
  }
  PetscFunctionReturn(0);
}
PetscErrorCode Core::add_ents_to_finite_element_by_TRIs(const EntityHandle meshset,const std::string &name,const bool recursive) {
  PetscFunctionBegin;
  *buildMoFEM &= 1<<0;
  EntityHandle idm = no_handle;
  try {
    idm = get_finite_element_meshset(get_BitFEId(name));
  } catch (MoFEMException const &e) {
    SETERRQ(PETSC_COMM_SELF,e.errorCode,e.errorMessage);
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
    SETERRQ(PETSC_COMM_SELF,e.errorCode,e.errorMessage);
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
    ierr = add_ents_to_finite_element_by_TETs(tets,get_BitFEId(name));  CHKERRQ(ierr);
  } catch (MoFEMException const &e) {
    SETERRQ(PETSC_COMM_SELF,e.errorCode,e.errorMessage);
  }
  PetscFunctionReturn(0);
}
PetscErrorCode Core::add_ents_to_finite_element_by_TETs(const EntityHandle meshset,const std::string &name,const bool recursive) {
  PetscFunctionBegin;
  *buildMoFEM &= 1<<0;
  try {
    ierr = add_ents_to_finite_element_by_TETs(meshset,get_BitFEId(name),recursive);  CHKERRQ(ierr);
  } catch (MoFEMException const &e) {
    SETERRQ(PETSC_COMM_SELF,e.errorCode,e.errorMessage);
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
    SETERRQ(PETSC_COMM_SELF,e.errorCode,e.errorMessage);
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
    ierr = add_ents_to_finite_element_by_PRISMs(prims,get_BitFEId(name));  CHKERRQ(ierr);
  } catch (MoFEMException const &e) {
    SETERRQ(PETSC_COMM_SELF,e.errorCode,e.errorMessage);
  }
  PetscFunctionReturn(0);
}
PetscErrorCode Core::add_ents_to_finite_element_by_PRISMs(const EntityHandle meshset,const std::string &name,const bool recursive) {
  PetscFunctionBegin;
  *buildMoFEM &= 1<<0;
  try {
    ierr = add_ents_to_finite_element_by_PRISMs(meshset,get_BitFEId(name),recursive);  CHKERRQ(ierr);
  } catch (MoFEMException const &e) {
    SETERRQ(PETSC_COMM_SELF,e.errorCode,e.errorMessage);
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
    const BitFEId id = get_BitFEId(name);
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
    SETERRQ(PETSC_COMM_SELF,e.errorCode,e.errorMessage);
  }
  PetscFunctionReturn(0);
}
PetscErrorCode Core::add_ents_to_finite_element_by_MESHSET(const EntityHandle meshset,const std::string& name,const bool recursive) {
  PetscFunctionBegin;
  *buildMoFEM &= 1<<0;
  const BitFEId id = get_BitFEId(name);
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
PetscErrorCode Core::modify_problem_add_finite_element(const std::string &name_problem,const std::string &fe_name) {
  PetscFunctionBegin;
  try {
    typedef MoFEMProblem_multiIndex::index<Problem_mi_tag>::type mofem_problems_by_name;
    mofem_problems_by_name& set = pRoblems.get<Problem_mi_tag>();
    mofem_problems_by_name::iterator miit = set.find(name_problem);
    if(miit==set.end()) {
      SETERRQ1(PETSC_COMM_SELF,MOFEM_NOT_FOUND,"this problem <%s> is not there",name_problem.c_str());
    }
    BitFEId f_id = get_BitFEId(fe_name);
    bool success = set.modify(miit,ProblemFiniteElementChangeBitAdd(f_id));
    if(!success) SETERRQ(PETSC_COMM_SELF,MOFEM_OPERATION_UNSUCCESSFUL,"modification unsuccessful");
  } catch (MoFEMException const &e) {
    SETERRQ(PETSC_COMM_SELF,e.errorCode,e.errorMessage);
  }
  PetscFunctionReturn(0);
}
PetscErrorCode Core::modify_problem_unset_finite_element(const std::string &name_problem,const std::string &fe_name) {
  PetscFunctionBegin;
  try {
    typedef MoFEMProblem_multiIndex::index<Problem_mi_tag>::type mofem_problems_by_name;
    mofem_problems_by_name& set = pRoblems.get<Problem_mi_tag>();
    mofem_problems_by_name::iterator miit = set.find(name_problem);
    if(miit==set.end()) {
      SETERRQ1(PETSC_COMM_SELF,MOFEM_NOT_FOUND,"this problem <%s> is not there",name_problem.c_str());
    }
    BitFEId f_id = get_BitFEId(fe_name);
    bool success = set.modify(miit,ProblemFiniteElementChangeBitUnSet(f_id));
    if(!success) SETERRQ(PETSC_COMM_SELF,MOFEM_OPERATION_UNSUCCESSFUL,"modification unsuccessful");
  } catch (MoFEMException const &e) {
    SETERRQ(PETSC_COMM_SELF,e.errorCode,e.errorMessage);
  }
  PetscFunctionReturn(0);
}
PetscErrorCode Core::modify_problem_ref_level_add_bit(const std::string &name_problem,const BitRefLevel &bit) {
  PetscFunctionBegin;
  typedef MoFEMProblem_multiIndex::index<Problem_mi_tag>::type mofem_problems_by_name;
  mofem_problems_by_name& set = pRoblems.get<Problem_mi_tag>();
  mofem_problems_by_name::iterator miit = set.find(name_problem);
  std::ostringstream ss;
  ss << name_problem;
  if(miit==set.end()) SETERRQ1(PETSC_COMM_SELF,MOFEM_NOT_FOUND,"this problem <%s> is there",ss.str().c_str());
  bool success = set.modify(miit,ProblemChangeRefLevelBitAdd(bit));
  if(!success) SETERRQ(PETSC_COMM_SELF,1,"modification unsuccessful");
  PetscFunctionReturn(0);
}
PetscErrorCode Core::modify_problem_ref_level_set_bit(const std::string &name_problem,const BitRefLevel &bit) {
  PetscFunctionBegin;
  typedef MoFEMProblem_multiIndex::index<Problem_mi_tag>::type mofem_problems_by_name;
  mofem_problems_by_name& set = pRoblems.get<Problem_mi_tag>();
  mofem_problems_by_name::iterator miit = set.find(name_problem);
  std::ostringstream ss;
  ss << name_problem;
  if(miit==set.end()) SETERRQ1(PETSC_COMM_SELF,MOFEM_NOT_FOUND,"this problem <%s> is there",ss.str().c_str());
  bool success = set.modify(miit,ProblemChangeRefLevelBitSet(bit));
  if(!success) SETERRQ(PETSC_COMM_SELF,MOFEM_OPERATION_UNSUCCESSFUL,"modification unsuccessful");
  PetscFunctionReturn(0);
}
PetscErrorCode Core::modify_problem_dof_mask_ref_level_set_bit(const std::string &name_problem,const BitRefLevel &bit) {
  PetscFunctionBegin;
  typedef MoFEMProblem_multiIndex::index<Problem_mi_tag>::type mofem_problems_by_name;
  mofem_problems_by_name& set = pRoblems.get<Problem_mi_tag>();
  mofem_problems_by_name::iterator miit = set.find(name_problem);
  std::ostringstream ss;
  ss << name_problem;
  if(miit==set.end()) SETERRQ1(PETSC_COMM_SELF,MOFEM_NOT_FOUND,"this problem <%s> is there",ss.str().c_str());
  bool success = set.modify(miit,ProblemChangeRefLevelBitDofMaskSet(bit));
  if(!success) SETERRQ(PETSC_COMM_SELF,MOFEM_OPERATION_UNSUCCESSFUL,"modification unsuccessful");
  PetscFunctionReturn(0);
}
PetscErrorCode Core::build_finite_elements(const boost::shared_ptr<FiniteElement> fe,const Range *ents_ptr,int verb) {
  PetscFunctionBegin;
  if(verb==-1) verb = verbose;

  typedef RefElement_multiIndex::index<Ent_mi_tag>::type RefFiniteElementByEnt;
  typedef Field_multiIndex::index<BitFieldId_mi_tag>::type FieldById;
  FieldById &fields_by_id = fIelds.get<BitFieldId_mi_tag>();

  //get id of mofem fields for row, col and data
  enum IntLoop { ROW = 0,COL,DATA,LAST };
  BitFieldId fe_fields[LAST] = {
    fe.get()->get_BitFieldId_row(),
    fe.get()->get_BitFieldId_col(),
    fe.get()->get_BitFieldId_data()
  };

  //get finite element meshset
  EntityHandle meshset = get_finite_element_meshset(fe.get()->getId());
  // get entities from finite element meshset // if meshset
  Range fe_ents;
  rval = moab.get_entities_by_handle(meshset,fe_ents,false); CHKERRQ_MOAB(rval);
  if(ents_ptr) fe_ents = intersect(fe_ents,*ents_ptr);

  // map entity uid to pointers
  std::map<UId,std::vector<boost::shared_ptr<EntFiniteElement> > > map_uid_fe;

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
      SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,ss.str().c_str());
    }
    boost::shared_ptr<EntFiniteElement> ent_fe = boost::shared_ptr<EntFiniteElement>(
      new EntFiniteElement(moab,ref_fe_miit->get_RefElement(),fe)
    );
    std::pair<EntFiniteElement_multiIndex::iterator,bool> p = entsFiniteElements.insert(ent_fe);

    if(fe_fields[ROW]==fe_fields[COL]) {
       p.first->get()->col_dof_view = p.first->get()->row_dof_view;
    }
    if(fe_fields[ROW]==fe_fields[DATA]) {
      p.first->get()->data_dof_view = p.first->get()->row_dof_view;
    } else if(fe_fields[COL]==fe_fields[DATA]) {
      p.first->get()->data_dof_view = p.first->get()->col_dof_view;
    }

    if(fe_fields[ROW]!=fe_fields[COL] && p.first->get()->col_dof_view == p.first->get()->row_dof_view) {
      p.first->get()->col_dof_view
      = boost::shared_ptr<DofEntity_multiIndex_uid_view>(new DofEntity_multiIndex_uid_view());
    }
    if(fe_fields[ROW]!=fe_fields[DATA] && p.first->get()->data_dof_view == p.first->get()->row_dof_view) {
      p.first->get()->data_dof_view
      = boost::shared_ptr<DofEntity_multiIndex_uid_view>(new DofEntity_multiIndex_uid_view());
    } else if(fe_fields[COL]!=fe_fields[DATA] && p.first->get()->data_dof_view == p.first->get()->col_dof_view) {
      p.first->get()->data_dof_view
      = boost::shared_ptr<DofEntity_multiIndex_uid_view>(new DofEntity_multiIndex_uid_view());
    }

    p.first->get()->row_dof_view->clear();
    p.first->get()->col_dof_view->clear();
    p.first->get()->data_dof_view->clear();
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
      FieldById::iterator miit = fields_by_id.find(BitFieldId().set(ii));
      if(miit==fields_by_id.end()) {
        SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"Data inconsistency");
      }

      // Entities adjacent to entities
      Range adj_ents;

      // Resolve entities on element, those entities are used to build tag with dof
      // uids on finite element tag
      ierr = p.first->get()->get_element_adjacency(moab,*miit,adj_ents); CHKERRQ(ierr);

      std::string field_name = miit->get()->getName();
      for(Range::iterator eit = adj_ents.begin();eit!=adj_ents.end();eit++) {
        MoFEMEntity_multiIndex::index<Composite_Name_And_Ent_mi_tag>::type::iterator meit;
        meit = entsFields.get<Composite_Name_And_Ent_mi_tag>().find(boost::make_tuple(field_name,*eit));
        if(meit!=entsFields.get<Composite_Name_And_Ent_mi_tag>().end()) {
          UId uid = meit->get()->getGlobalUniqueId();
          map_uid_fe[uid].push_back(*p.first);
        }
      }
    }

  }

  typedef DofEntity_multiIndex::index<Unique_Ent_mi_tag>::type DofsByEntUId;
  DofsByEntUId& dofs_by_ent_uid = dofsField.get<Unique_Ent_mi_tag>();

  boost::shared_ptr<SideNumber> side_number_ptr;
  for(
    std::map<UId,std::vector<boost::shared_ptr<EntFiniteElement> > >::iterator
    mit = map_uid_fe.begin();mit!=map_uid_fe.end();mit++
  ) {
    DofsByEntUId::iterator dit,hi_dit;
    dit = dofs_by_ent_uid.lower_bound(mit->first);
    hi_dit = dofs_by_ent_uid.upper_bound(mit->first);
    for(;dit!=hi_dit;dit++) {
      // cerr << mit->first << endl;
      // cerr << **dit << endl;
      BitFieldId field_id = dit->get()->getId();
      std::vector<boost::shared_ptr<EntFiniteElement> >::iterator fe_it,hi_fe_it;
      fe_it = mit->second.begin();
      hi_fe_it = mit->second.end();
      for(;fe_it!=hi_fe_it;fe_it++) {
        if((field_id&fe_it->get()->get_BitFieldId_row()).any()) {
          fe_it->get()->row_dof_view->insert(*dit);
        }
        if(fe_it->get()->col_dof_view!=fe_it->get()->row_dof_view) {
          if((field_id&fe_it->get()->get_BitFieldId_col()).any()) {
            fe_it->get()->col_dof_view->insert(*dit);
          }
        }
        if(fe_it->get()->data_dof_view!=fe_it->get()->row_dof_view && fe_it->get()->data_dof_view!=fe_it->get()->col_dof_view) {
          if((field_id&fe_it->get()->get_BitFieldId_data()).any()) {
            fe_it->get()->data_dof_view->insert(*dit);
          }
        }
        if((field_id&fe_it->get()->get_BitFieldId_data()).any()) {
          side_number_ptr = fe_it->get()->get_side_number_ptr(moab,(*dit)->getEnt());
          //add dofs to finite element multi_index database
          fe_it->get()->data_dofs.get<Unique_mi_tag>().insert(
            boost::shared_ptr<FEDofEntity>(new FEDofEntity(side_number_ptr,*dit))
          );
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
PetscErrorCode Core::build_finite_elements(const string fe_name,const Range *ents_ptr,int verb) {
  PetscFunctionBegin;
  if(verb==-1) verb = verbose;

  FiniteElement_multiIndex::index<FiniteElement_name_mi_tag>::type::iterator
  fe_miit = finiteElements.get<FiniteElement_name_mi_tag>().find(fe_name);
  if(fe_miit==finiteElements.get<FiniteElement_name_mi_tag>().end()) {
    SETERRQ1(
      PETSC_COMM_SELF,MOFEM_NOT_FOUND,
      "Finite element <%s> not found",fe_name.c_str()
    );
  }

  if(verb>0) PetscPrintf(comm,"Build Finite Elements %s\n",fe_name.c_str());
  ierr = build_finite_elements(*fe_miit,ents_ptr,verb); CHKERRQ(ierr);

  *buildMoFEM |= 1<<1;
  PetscFunctionReturn(0);
}

PetscErrorCode Core::build_adjacencies(const Range &ents,int verb) {
  PetscFunctionBegin;
  if(verb==-1) verb = verbose;
  if(!((*buildMoFEM)&(1<<0))) SETERRQ(PETSC_COMM_SELF,MOFEM_NOT_FOUND,"field not build");
  if(!((*buildMoFEM)&(1<<1))) SETERRQ(PETSC_COMM_SELF,MOFEM_NOT_FOUND,"fe not build");
  //typedef MoFEMEntity_multiIndex::index<Unique_mi_tag>::type ents_by_uid;
  EntFiniteElement_multiIndex::iterator fit = entsFiniteElements.begin();
  for(;fit!=entsFiniteElements.end();fit++) {
    if(!ents.empty()) {
      if(ents.find((*fit)->getEnt())==ents.end()) {
        continue;
      }
    }
    int by = BYROW;
    if((*fit)->get_BitFieldId_row()!=(*fit)->get_BitFieldId_col()) by |= BYCOL;
    if((*fit)->get_BitFieldId_row()!=(*fit)->get_BitFieldId_data()) by |= BYDATA;
    MoFEMEntityEntFiniteElementAdjacencyMap_change_ByWhat modify_row(by);
    GlobalUId ent_uid = UId(0);
    DofEntity_multiIndex_uid_view::iterator rvit;
    rvit = (*fit)->row_dof_view->begin();
    for(;rvit!=(*fit)->row_dof_view->end();rvit++) {
      if( ent_uid == (*rvit)->getMoFEMEntityPtr()->getGlobalUniqueId()) continue;
      ent_uid = (*rvit)->getMoFEMEntityPtr()->getGlobalUniqueId();
      std::pair<MoFEMEntityEntFiniteElementAdjacencyMap_multiIndex::iterator,bool> p;
      p = entFEAdjacencies.insert(
        MoFEMEntityEntFiniteElementAdjacencyMap((*rvit)->getMoFEMEntityPtr(),*fit)
      );
      bool success = entFEAdjacencies.modify(p.first,modify_row);
      if(!success) SETERRQ(PETSC_COMM_SELF,MOFEM_OPERATION_UNSUCCESSFUL,"modification unsuccessful");
    }
    if((*fit)->get_BitFieldId_row()!=(*fit)->get_BitFieldId_col()) {
      int by = BYCOL;
      if((*fit)->get_BitFieldId_col()!=(*fit)->get_BitFieldId_data()) by |= BYDATA;
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
        if(!success) SETERRQ(PETSC_COMM_SELF,MOFEM_OPERATION_UNSUCCESSFUL,"modification unsuccessful");
      }
    }
    if(
      (*fit)->get_BitFieldId_row()!=(*fit)->get_BitFieldId_data()||
      (*fit)->get_BitFieldId_col()!=(*fit)->get_BitFieldId_data()
    ) {
      MoFEMEntityEntFiniteElementAdjacencyMap_change_ByWhat modify_data(BYDATA);
      ent_uid = UId(0);
      DofEntity_multiIndex_uid_view::iterator dvit;
      dvit = (*fit)->data_dof_view->begin();
      for(;dvit!=(*fit)->data_dof_view->end();dvit++) {
        if( ent_uid == (*dvit)->getMoFEMEntityPtr()->getGlobalUniqueId()) continue;
        ent_uid = (*dvit)->getMoFEMEntityPtr()->getGlobalUniqueId();
        std::pair<MoFEMEntityEntFiniteElementAdjacencyMap_multiIndex::iterator,bool> p;
        p = entFEAdjacencies.insert(
          MoFEMEntityEntFiniteElementAdjacencyMap((*dvit)->getMoFEMEntityPtr(),*fit)
        );
        bool success = entFEAdjacencies.modify(p.first,modify_data);
        if(!success) SETERRQ(PETSC_COMM_SELF,MOFEM_OPERATION_UNSUCCESSFUL,"modification unsuccessful");
      }
    }
  }
  if(verbose>1) {
    list_adjacencies();
  }
  if(verbose>0) {
    PetscSynchronizedPrintf(comm,"Nb. entFEAdjacencies %u\n",entFEAdjacencies.size());
    PetscSynchronizedFlush(comm,PETSC_STDOUT);
    //PetscSynchronizedPrintf(comm,"Nb. entFEAdjacencies %u\n",entFEAdjacencies.size());
    //PetscSynchronizedFlush(comm,PETSC_STDOUT);
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
PetscErrorCode Core::partition_finite_elements(
  const std::string &name,bool part_from_moab,int low_proc,int hi_proc,int verb
) {
  PetscFunctionBegin;
  if(verb==-1) verb = verbose;

  if(!(*buildMoFEM&BUILD_FIELD)) SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"fields not build");
  if(!(*buildMoFEM&BUILD_FE)) SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"FEs not build");
  if(!(*buildMoFEM&BUILD_ADJ)) SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"adjacencies not build");
  if(!(*buildMoFEM&BUILD_PROBLEM)) SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"problem not build");
  if(!(*buildMoFEM&PARTITION_PROBLEM)) SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"problem not partitioned");

  if(low_proc == -1) low_proc = rAnk;
  if(hi_proc == -1) hi_proc = rAnk;
  typedef MoFEMProblem_multiIndex::index<Problem_mi_tag>::type mofem_problems_by_name;
  //find p_miit
  mofem_problems_by_name &pRoblems_set = pRoblems.get<Problem_mi_tag>();
  mofem_problems_by_name::iterator p_miit = pRoblems_set.find(name);
  if(p_miit == pRoblems_set.end()) {
    SETERRQ1(
      PETSC_COMM_SELF,1,"problem < %s > not found (top tip: check spelling)",name.c_str()
    );
  }
  NumeredEntFiniteElement_multiIndex& numeredFiniteElements
  = const_cast<NumeredEntFiniteElement_multiIndex&>(p_miit->numeredFiniteElements);
  numeredFiniteElements.clear();

  bool do_cols_prob = true;
  if(p_miit->numered_dofs_rows == p_miit->numered_dofs_cols) {
    do_cols_prob = false;
  }

  //FiniteElement set
  EntFiniteElement_multiIndex::iterator miit2 = entsFiniteElements.begin();
  EntFiniteElement_multiIndex::iterator hi_miit2 = entsFiniteElements.end();
  for(;miit2!=hi_miit2;miit2++) {
    if(((*miit2)->getId()&p_miit->get_BitFEId()).none()) continue; // if element is not part of problem
    if(((*miit2)->getBitRefLevel()&p_miit->getBitRefLevel())!=p_miit->getBitRefLevel()) continue; // if entity is not problem refinement level
    boost::shared_ptr<NumeredEntFiniteElement> numered_fe(new NumeredEntFiniteElement(*miit2));
    bool do_cols_fe = true;
    if(numered_fe->sPtr->row_dof_view == numered_fe->sPtr->col_dof_view && do_cols_prob) {
      do_cols_fe = false;
      numered_fe->cols_dofs = numered_fe->rows_dofs;
    } else {
      numered_fe->cols_dofs = boost::shared_ptr<FENumeredDofEntity_multiIndex>(new FENumeredDofEntity_multiIndex());
    }
    boost::shared_ptr<FENumeredDofEntity_multiIndex> rows_dofs = numered_fe->rows_dofs;
    boost::shared_ptr<FENumeredDofEntity_multiIndex> cols_dofs = numered_fe->cols_dofs;
    rows_dofs->clear();
    if(do_cols_fe) {
      cols_dofs->clear();
    }
    {
      NumeredDofEntity_multiIndex_uid_view_ordered rows_view;
      NumeredDofEntity_multiIndex_uid_view_ordered::iterator viit_rows;
      if(part_from_moab) {
        int proc = (*miit2)->getOwnerProc();
        NumeredEntFiniteElement_change_part(proc).operator()(numered_fe);
      } else {
        //rows_view
        ierr = (*miit2)->get_MoFEMFiniteElement_row_dof_view(
          *(p_miit->numered_dofs_rows),rows_view,Interface::UNION
        ); CHKERRQ(ierr);
        std::vector<int> parts(sIze,0);
        viit_rows = rows_view.begin();
        for(;viit_rows!=rows_view.end();viit_rows++) {
          parts[(*viit_rows)->part]++;
        }
        std::vector<int>::iterator pos = max_element(parts.begin(),parts.end());
        unsigned int max_part = distance(parts.begin(),pos);
        NumeredEntFiniteElement_change_part(max_part).operator()(numered_fe);
      }
      if(
        (numered_fe->get_part()>=(unsigned int)low_proc)&&
        (numered_fe->get_part()<=(unsigned int)hi_proc)
      ) {
        if(part_from_moab) {
          //rows_view
          ierr = (*miit2)->get_MoFEMFiniteElement_row_dof_view(
            *(p_miit->numered_dofs_rows),rows_view,Interface::UNION
          ); CHKERRQ(ierr);
        }
        //rows element dof multi-indices
        viit_rows = rows_view.begin();
        for(;viit_rows!=rows_view.end();viit_rows++) {
          try {
            boost::shared_ptr<SideNumber> side_number_ptr = (*miit2)->get_side_number_ptr(moab,(*viit_rows)->getEnt());
            rows_dofs->insert(boost::shared_ptr<FENumeredDofEntity>(new FENumeredDofEntity(side_number_ptr,*viit_rows)));
          } catch (MoFEMException const &e) {
            SETERRQ(PETSC_COMM_SELF,e.errorCode,e.errorMessage);
          }
        }
        if(do_cols_fe) {
          //cols_views
          NumeredDofEntity_multiIndex_uid_view_ordered cols_view;
          ierr = (*miit2)->get_MoFEMFiniteElement_col_dof_view(
            *(p_miit->numered_dofs_cols),cols_view,Interface::UNION
          ); CHKERRQ(ierr);
          //cols element dof multi-indices
          NumeredDofEntity_multiIndex_uid_view_ordered::iterator viit_cols;;
          viit_cols = cols_view.begin();
          for(;viit_cols!=cols_view.end();viit_cols++) {
            try {
              boost::shared_ptr<SideNumber> side_number_ptr = (*miit2)->get_side_number_ptr(moab,(*viit_cols)->getEnt());
              cols_dofs->insert(boost::shared_ptr<FENumeredDofEntity>(new FENumeredDofEntity(side_number_ptr,*viit_cols)));
            } catch (MoFEMException const &e) {
              SETERRQ(PETSC_COMM_SELF,e.errorCode,e.errorMessage);
            }
          }
        }
      }
      std::pair<NumeredEntFiniteElement_multiIndex::iterator,bool> p;
      p = numeredFiniteElements.insert(numered_fe);
      if(!p.second) {
        SETERRQ(PETSC_COMM_SELF,MOFEM_NOT_FOUND,"element is there");
      }
      if(verb>1) {
        std::ostringstream ss;
        ss << *p_miit << std::endl;
        ss << *p.first << std::endl;
        typedef FENumeredDofEntity_multiIndex::index<Unique_mi_tag>::type FENumeredDofEntity_multiIndex_by_Unique_mi_tag;
        FENumeredDofEntity_multiIndex_by_Unique_mi_tag::iterator miit = (*p.first)->rows_dofs->get<Unique_mi_tag>().begin();
        for(;miit!= (*p.first)->rows_dofs->get<Unique_mi_tag>().end();miit++) ss << "rows: " << *(*miit) << std::endl;
        miit = (*p.first)->cols_dofs->get<Unique_mi_tag>().begin();
        for(;miit!=(*p.first)->cols_dofs->get<Unique_mi_tag>().end();miit++) ss << "cols: " << *(*miit) << std::endl;
        PetscSynchronizedPrintf(comm,ss.str().c_str());
      }
    }
  }
  if(verb>0) {
    typedef NumeredEntFiniteElement_multiIndex::index<FiniteElement_Part_mi_tag>::type NumeredEntFiniteElement_multiIndex_by_part;
    NumeredEntFiniteElement_multiIndex_by_part::iterator MoFEMFiniteElement_miit = numeredFiniteElements.get<FiniteElement_Part_mi_tag>().lower_bound(rAnk);
    NumeredEntFiniteElement_multiIndex_by_part::iterator hi_MoMoFEMFiniteElement_miitFEMFE_miit = numeredFiniteElements.get<FiniteElement_Part_mi_tag>().upper_bound(rAnk);
    int count = distance(MoFEMFiniteElement_miit,hi_MoMoFEMFiniteElement_miitFEMFE_miit);
    std::ostringstream ss;
    ss << *p_miit;
    ss << " Nb. elems " << count << " on proc " << rAnk << std::endl;
    PetscSynchronizedPrintf(comm,ss.str().c_str());
    PetscSynchronizedFlush(comm,PETSC_STDOUT);
  }
  *buildMoFEM |= PARTITION_FE;
  PetscFunctionReturn(0);
}
PetscErrorCode Core::partition_ghost_dofs(const std::string &name,int verb) {
  PetscFunctionBegin;
  if(verb==-1) verb = verbose;
  if(!(*buildMoFEM&BUILD_FIELD)) SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"fields not build");
  if(!(*buildMoFEM&BUILD_FE)) SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"FEs not build");
  if(!(*buildMoFEM&BUILD_ADJ)) SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"adjacencies not build");
  if(!(*buildMoFEM&BUILD_PROBLEM)) SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"problem not build");
  if(!(*buildMoFEM&PARTITION_PROBLEM)) SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"partition of problem not build");
  if(!(*buildMoFEM&PARTITION_FE)) SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"partitions finite elements not build");
  typedef MoFEMProblem_multiIndex::index<Problem_mi_tag>::type mofem_problems_by_name;
  //find p_miit
  mofem_problems_by_name &pRoblems_set = pRoblems.get<Problem_mi_tag>();
  mofem_problems_by_name::iterator p_miit = pRoblems_set.find(name);
  //
  DofIdx &nb_row_ghost_dofs = *((DofIdx*)p_miit->tag_ghost_nbdof_data_row);
  DofIdx &nb_col_ghost_dofs = *((DofIdx*)p_miit->tag_ghost_nbdof_data_col);
  //
  nb_row_ghost_dofs = 0;
  nb_col_ghost_dofs = 0;
  if(sIze>1) {
    NumeredDofEntity_multiIndex_uid_view_ordered ghost_idx_col_view,ghost_idx_row_view;
    NumeredEntFiniteElement_multiIndex::index<FiniteElement_Part_mi_tag>::type::iterator fe_it,hi_fe_it;
    fe_it = p_miit->numeredFiniteElements.get<FiniteElement_Part_mi_tag>().lower_bound(rAnk);
    hi_fe_it = p_miit->numeredFiniteElements.get<FiniteElement_Part_mi_tag>().upper_bound(rAnk);
    for(;fe_it!=hi_fe_it;fe_it++) {
      typedef FENumeredDofEntity_multiIndex::iterator dof_it;
      if((*fe_it)->rows_dofs->size()>0) {
        dof_it rowdofit,hi_rowdofit;
        rowdofit = (*fe_it)->rows_dofs->begin();
        hi_rowdofit = (*fe_it)->rows_dofs->end();
        for(;rowdofit!=hi_rowdofit;rowdofit++) {
          if((*rowdofit)->get_part()==(unsigned int)rAnk) continue;
          ghost_idx_row_view.insert((*rowdofit)->get_NumeredDofEntity_ptr());
        }
      }
      if((*fe_it)->cols_dofs->size()>0) {
        dof_it coldofit,hi_coldofit;
        coldofit = (*fe_it)->cols_dofs->begin();
        hi_coldofit = (*fe_it)->cols_dofs->end();
        for(;coldofit!=hi_coldofit;coldofit++) {
          if((*coldofit)->get_part()==(unsigned int)rAnk) continue;
          ghost_idx_col_view.insert((*coldofit)->get_NumeredDofEntity_ptr());
        }
      }
    }
    DofIdx *nb_ghost_dofs[2] = { &nb_col_ghost_dofs, &nb_row_ghost_dofs };
    DofIdx nb_local_dofs[2] = { *((DofIdx*)p_miit->tag_local_nbdof_data_col), *((DofIdx*)p_miit->tag_local_nbdof_data_row) };
    NumeredDofEntity_multiIndex_uid_view_ordered *ghost_idx_view[2] = { &ghost_idx_col_view, &ghost_idx_row_view };
    typedef NumeredDofEntity_multiIndex::index<Unique_mi_tag>::type NumeredDofEntitys_by_unique_id;
    NumeredDofEntitys_by_unique_id *dof_by_uid_no_const[2] = {
      const_cast<NumeredDofEntitys_by_unique_id*>(&p_miit->numered_dofs_cols->get<Unique_mi_tag>()),
      const_cast<NumeredDofEntitys_by_unique_id*>(&p_miit->numered_dofs_rows->get<Unique_mi_tag>())
    };
    int loop_size = 2;
    if(p_miit->numered_dofs_cols==p_miit->numered_dofs_rows) {
      loop_size = 1;
    }
    for(int ss = 0;ss<loop_size;ss++) {
      NumeredDofEntity_multiIndex_uid_view_ordered::iterator ghost_idx_miit = ghost_idx_view[ss]->begin();
      for(;ghost_idx_miit!=ghost_idx_view[ss]->end();ghost_idx_miit++) {
        NumeredDofEntitys_by_unique_id::iterator diit = dof_by_uid_no_const[ss]->find((*ghost_idx_miit)->getGlobalUniqueId());
        if((*diit)->petsc_local_dof_idx!=(DofIdx)-1) {
          SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"inconsistent data, ghost dof already set");
        }
        bool success = dof_by_uid_no_const[ss]->modify(diit,NumeredDofEntity_local_idx_change(nb_local_dofs[ss]++));
        if(!success) SETERRQ(PETSC_COMM_SELF,MOFEM_OPERATION_UNSUCCESSFUL,"modification unsuccessful");
        (*nb_ghost_dofs[ss])++;
      }
    }
    if(loop_size==1) {
      (*nb_ghost_dofs[1]) = (*nb_ghost_dofs[0]);
    }
  }
  if(verb>0) {
    std::ostringstream ss;
    ss << "partition_ghost_col_dofs: rank = " << rAnk
    << " FEs col ghost dofs "<< *p_miit
    << " Nb. col ghost dof " << p_miit->get_nb_ghost_dofs_col()
    << " Nb. local dof " << p_miit->get_nb_local_dofs_col() << std::endl;
    ss << "partition_ghost_row_dofs: rank = " << rAnk
    << " FEs row ghost dofs "<< *p_miit
    << " Nb. row ghost dof " << p_miit->get_nb_ghost_dofs_row()
    << " Nb. local dof " << p_miit->get_nb_local_dofs_row() << std::endl;
    if(verb>1) {
      NumeredDofEntity_multiIndex::iterator miit_dd_col = p_miit->numered_dofs_cols->begin();
      for(;miit_dd_col!=p_miit->numered_dofs_cols->end();miit_dd_col++) {
        if((*miit_dd_col)->part==(unsigned int)rAnk) continue;
        if((*miit_dd_col)->petsc_local_dof_idx==(DofIdx)-1) continue;
        ss<<*(*miit_dd_col)<<std::endl;
      }
    }
    PetscSynchronizedPrintf(comm,ss.str().c_str());
    PetscSynchronizedFlush(comm,PETSC_STDOUT);
  }
  *buildMoFEM |= PARTITION_GHOST_DOFS;
  PetscFunctionReturn(0);
}
PetscErrorCode Core::seed_finite_elements(const EntityHandle meshset,int verb) {
  PetscFunctionBegin;
  Range entities;
  ierr = moab.get_entities_by_handle(meshset,entities,true); CHKERRQ(ierr);
  ierr = seed_finite_elements(entities,verb); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
PetscErrorCode Core::seed_finite_elements(const Range &entities,int verb) {
  PetscFunctionBegin;
  for(Range::iterator eit = entities.begin();eit!=entities.end();eit++) {
    RefEntity_multiIndex::index<Ent_mi_tag>::type::iterator
      eiit = refinedEntities.get<Ent_mi_tag>().find(*eit);
    if(eiit == refinedEntities.get<Ent_mi_tag>().end())  {
      SETERRQ(PETSC_COMM_SELF,MOFEM_NOT_FOUND,"entity is not in database");
    }
    if((*eiit)->getBitRefLevel().none()) continue;
    std::pair<RefElement_multiIndex::iterator,bool> p_MoFEMFiniteElement;
    switch ((*eiit)->getEntType()) {
      case MBVERTEX:
      p_MoFEMFiniteElement = refinedFiniteElements.insert(ptrWrapperRefElement(
        boost::shared_ptr<RefElement>(new RefElement_VERTEX(moab,*eiit)))
      );
      break;
      case MBEDGE:
      p_MoFEMFiniteElement = refinedFiniteElements.insert(ptrWrapperRefElement(
        boost::shared_ptr<RefElement>(new RefElement_EDGE(moab,*eiit)))
      );
      break;
      case MBTRI:
      p_MoFEMFiniteElement = refinedFiniteElements.insert(ptrWrapperRefElement(
        boost::shared_ptr<RefElement>(new RefElement_TRI(moab,*eiit)))
      );
      break;
      case MBTET:
      p_MoFEMFiniteElement = refinedFiniteElements.insert(ptrWrapperRefElement(
        boost::shared_ptr<RefElement>(new RefElement_TET(moab,*eiit)))
      );
      break;
      case MBPRISM:
      p_MoFEMFiniteElement = refinedFiniteElements.insert(ptrWrapperRefElement(
        boost::shared_ptr<RefElement>(new RefElement_PRISM(moab,*eiit)))
      );
      break;
      default:
      SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"not implemented");
    }
  }
  PetscFunctionReturn(0);
}
PetscErrorCode Core::seed_ref_level_2D(const EntityHandle meshset,const BitRefLevel &bit,int verb) {
  PetscFunctionBegin;
  Range ents2d;
  rval = moab.get_entities_by_dimension(meshset,2,ents2d,false); CHKERRQ_MOAB(rval);
  ierr = seed_ref_level(ents2d,bit,verb); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
PetscErrorCode Core::seed_ref_level(const Range &ents,const BitRefLevel &bit,int verb) {
  PetscFunctionBegin;
  if(verb==-1) verb = verbose;
  Range seeded_ents;
  try {
    if(verb > 1) {
      PetscSynchronizedPrintf(comm,"nb. entities for seed %d\n",ents.size());
    }
    Range::iterator tit = ents.begin();
    for(;tit!=ents.end();tit++) {
      boost::shared_ptr<RefEntity> ref_ent(new RefEntity(basicEntityDataPtr,*tit));
      std::bitset<8> ent_pstat(ref_ent->getPStatus());
      ent_pstat.flip(0);
      std::pair<RefEntity_multiIndex::iterator,bool> p_ent = refinedEntities.insert(ref_ent);
      if(debug > 0) {
        ierr = test_moab(moab,*tit); CHKERRQ(ierr);
      }
      bool success = refinedEntities.modify(p_ent.first,RefEntity_change_add_bit(bit));
      if(!success) SETERRQ(PETSC_COMM_SELF,MOFEM_OPERATION_UNSUCCESSFUL,"modification unsuccessful");
      if(verb>2) {
        std::ostringstream ss;
        ss << **p_ent.first;
        PetscSynchronizedPrintf(comm,"%s\n",ss.str().c_str());
      }
      std::pair<RefElement_multiIndex::iterator,bool> p_MoFEMFiniteElement;
      switch((*p_ent.first)->getEntType()) {
        case MBVERTEX:
        p_MoFEMFiniteElement = refinedFiniteElements.insert(ptrWrapperRefElement(
          boost::shared_ptr<RefElement>(new RefElement_VERTEX(moab,*p_ent.first)))
        );
        seeded_ents.insert(*tit);
        break;
        case MBEDGE:
        p_MoFEMFiniteElement = refinedFiniteElements.insert(ptrWrapperRefElement(
          boost::shared_ptr<RefElement>(new RefElement_EDGE(moab,*p_ent.first)))
        );
        seeded_ents.insert(*tit);
        break;
        case MBTRI:
        p_MoFEMFiniteElement = refinedFiniteElements.insert(ptrWrapperRefElement(
          boost::shared_ptr<RefElement>(new RefElement_TRI(moab,*p_ent.first)))
        );
        seeded_ents.insert(*tit);
        break;
        case MBTET:
        p_MoFEMFiniteElement = refinedFiniteElements.insert(ptrWrapperRefElement(
          boost::shared_ptr<RefElement>(new RefElement_TET(moab,*p_ent.first)))
        );
        seeded_ents.insert(*tit);
        break;
        case MBPRISM:
        p_MoFEMFiniteElement = refinedFiniteElements.insert(ptrWrapperRefElement(
          boost::shared_ptr<RefElement>(new RefElement_PRISM(moab,*p_ent.first)))
        );
        break;
        case MBENTITYSET:
        p_MoFEMFiniteElement = refinedFiniteElements.insert(ptrWrapperRefElement(
          boost::shared_ptr<RefElement>(new RefElement_MESHSET(moab,*p_ent.first)))
        );
        break;
        default:
        SETERRQ(PETSC_COMM_SELF,MOFEM_NOT_IMPLEMENTED,"not implemented");
      }
      if(verb>3) {
        std::ostringstream ss;
        ss << *(p_MoFEMFiniteElement.first->get_RefElement());
        PetscSynchronizedPrintf(comm,"%s\n",ss.str().c_str());
      }
    }
    if(!seeded_ents.empty()) {
      int dim = moab.dimension_from_handle(seeded_ents[0]);
      for(int dd = 0;dd<dim;dd++) {
        Range ents;
        rval = moab.get_adjacencies(seeded_ents,dd,true,ents,Interface::UNION); CHKERRQ_MOAB(rval);
        if(dd == 2) {
          // currently only works with triangles
          ents = ents.subset_by_type(MBTRI);
        }
        Range::iterator eit = ents.begin();
        for(;eit!=ents.end();eit++) {
          std::pair<RefEntity_multiIndex::iterator,bool> p_ent = refinedEntities.insert(
            boost::shared_ptr<RefEntity>(new RefEntity(basicEntityDataPtr,*eit))
          );
          bool success = refinedEntities.modify(p_ent.first,RefEntity_change_add_bit(bit));
          if(!success) SETERRQ(PETSC_COMM_SELF,MOFEM_OPERATION_UNSUCCESSFUL,"modification unsuccessful");
          if(verb>2) {
            std::ostringstream ss;
            ss << *(p_ent.first);
            PetscSynchronizedPrintf(comm,"%s\n",ss.str().c_str());
          }
        }
      }
    }
    if(verb>2) {
      PetscSynchronizedPrintf(comm,"\n");
      PetscSynchronizedFlush(comm,PETSC_STDOUT);
    }
  } catch (MoFEMException const &e) {
    SETERRQ(PETSC_COMM_SELF,e.errorCode,e.errorMessage);
  }
  PetscFunctionReturn(0);
}
PetscErrorCode Core::seed_ref_level_3D(const EntityHandle meshset,const BitRefLevel &bit,int verb) {
  PetscFunctionBegin;
  Range ents3d;
  rval = moab.get_entities_by_dimension(meshset,3,ents3d,false); CHKERRQ_MOAB(rval);
  ierr = seed_ref_level(ents3d,bit,verb); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
PetscErrorCode Core::seed_ref_level_MESHSET(const EntityHandle meshset,const BitRefLevel &bit,int verb) {
  PetscFunctionBegin;
  if(verb==-1) verb = verbose;
  std::pair<RefEntity_multiIndex::iterator,bool> p_ent = refinedEntities.insert(
    boost::shared_ptr<RefEntity>(new RefEntity(basicEntityDataPtr,meshset))
  );
  refinedEntities.modify(p_ent.first,RefEntity_change_add_bit(bit));
  ptrWrapperRefElement pack_fe(
    boost::shared_ptr<RefElement>(new RefElement_MESHSET(moab,*p_ent.first))
  );
  std::pair<RefElement_multiIndex::iterator,bool> p_MoFEMFiniteElement = refinedFiniteElements.insert(pack_fe);
  if(verbose > 0) {
    std::ostringstream ss;
    ss << "add meshset as ref_ent " << *(p_MoFEMFiniteElement.first->get_RefElement()) << std::endl;
    PetscPrintf(comm,ss.str().c_str());
  }
  PetscFunctionReturn(0);
}
PetscErrorCode Core::shift_left_bit_ref(const int shift,int verb) {
  PetscFunctionBegin;
  SETERRQ(PETSC_COMM_SELF,MOFEM_NOT_IMPLEMENTED,"not implemented");
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
    if(!success) SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"inconsistency in data");
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

  typedef RefEntity_multiIndex::index<Composite_Ent_And_ParentEntType_mi_tag>::type ref_ents_by_composite;
  ref_ents_by_composite &ref_ents = refinedEntities.get<Composite_Ent_And_ParentEntType_mi_tag>();
  Range::iterator eit = ents.begin();
  for(;eit!=ents.end();eit++) {
    if(verb>2) {
      std::ostringstream ss;
      ss << "ent " << *eit << std::endl;;
      PetscPrintf(comm,ss.str().c_str());
    }
    ref_ents_by_composite::iterator miit = ref_ents.lower_bound(boost::make_tuple(*eit,child_type));
    ref_ents_by_composite::iterator hi_miit = ref_ents.upper_bound(boost::make_tuple(*eit,child_type));
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
          SETERRQ(PETSC_COMM_SELF,MOFEM_IMPOSIBLE_CASE,"this should not happen");
        }
        if(moab.type_from_handle(*eit)==MBENTITYSET) {
          SETERRQ(PETSC_COMM_SELF,MOFEM_IMPOSIBLE_CASE,"this should not happen");
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
    EntityHandle meshset = (*fit)->getMeshSet();
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
  EntityHandle meshset = (*fit)->getMeshSet();
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
  EntityHandle meshset = (*miit)->getMeshSet();
  ierr = update_meshset_by_entities_children(meshset,child_bit,meshset,fe_ent_type,false,verb);  CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
PetscErrorCode Core::get_problem_finite_elements_entities(const std::string &problem_name,const std::string &fe_name,const EntityHandle meshset) {
  PetscFunctionBegin;
  typedef MoFEMProblem_multiIndex::index<Problem_mi_tag>::type mofem_problems_by_name;
  mofem_problems_by_name &pRoblems_set = pRoblems.get<Problem_mi_tag>();
  mofem_problems_by_name::iterator p_miit = pRoblems_set.find(problem_name);
  if(p_miit == pRoblems_set.end()) SETERRQ1(PETSC_COMM_SELF,1,"no such problem like < %s >",problem_name.c_str());
  NumeredEntFiniteElement_multiIndex &numeredFiniteElements = const_cast<NumeredEntFiniteElement_multiIndex&>(p_miit->numeredFiniteElements);
  NumeredEntFiniteElement_multiIndex::index<FiniteElement_name_mi_tag>::type::iterator miit = numeredFiniteElements.get<FiniteElement_name_mi_tag>().lower_bound(fe_name);
  for(;miit!=numeredFiniteElements.get<FiniteElement_name_mi_tag>().upper_bound(fe_name);miit++) {
    EntityHandle ent = (*miit)->getEnt();
    rval = moab.add_entities(meshset,&ent,1); CHKERRQ_MOAB(rval);
    int part = (*miit)->get_part();
    rval = moab.tag_set_data(th_Part,&ent,1,&part); CHKERRQ_MOAB(rval);
  }
  PetscFunctionReturn(0);
}
bool Core::check_msId_meshset(const int ms_id,const CubitBCType cubit_bc_type) {
  CubitMeshSet_multiIndex::index<Composite_Cubit_msId_And_MeshSetType_mi_tag>::type::iterator
    miit = cubitMeshsets.get<Composite_Cubit_msId_And_MeshSetType_mi_tag>().find(boost::make_tuple(ms_id,cubit_bc_type.to_ulong()));
  if(miit!=cubitMeshsets.get<Composite_Cubit_msId_And_MeshSetType_mi_tag>().end()) {
    return true;
  }
  return false;
}
PetscErrorCode Core::add_cubit_msId(const CubitBCType cubit_bc_type,const int ms_id,const std::string name) {
  PetscFunctionBegin;
  if(check_msId_meshset(ms_id,cubit_bc_type)) {
    SETERRQ1(PETSC_COMM_SELF,1,"such cubit meshset is already there",ms_id);
  }
  try {
    CubitMeshSets cmeshset(moab,cubit_bc_type,ms_id);
    if((cmeshset.cubitBcType&CubitBCType(NODESET|SIDESET|BLOCKSET)).any()) {
      std::pair<CubitMeshSet_multiIndex::iterator,bool> p = cubitMeshsets.insert(cmeshset);
      if(!p.second) {
        SETERRQ(PETSC_COMM_SELF,MOFEM_OPERATION_UNSUCCESSFUL,"meshset not inserted");
      }
      if(name.size()>0) {
        bool success  = cubitMeshsets.modify(p.first,CubitMeshSets_change_name(moab,name));
        if(!success) {
          SETERRQ(PETSC_COMM_SELF,MOFEM_OPERATION_UNSUCCESSFUL,"name to cubit meshset can not be set");
        }
      }
    }
  } catch (MoFEMException const &e) {
    SETERRQ(PETSC_COMM_SELF,e.errorCode,e.errorMessage);
  }
  PetscFunctionReturn(0);
}
PetscErrorCode Core::delete_cubit_msId(const CubitBCType cubit_bc_type,const int ms_id) {
  PetscFunctionBegin;
  CubitMeshSet_multiIndex::index<Composite_Cubit_msId_And_MeshSetType_mi_tag>::type::iterator
    miit = cubitMeshsets.get<Composite_Cubit_msId_And_MeshSetType_mi_tag>().find(boost::make_tuple(ms_id,cubit_bc_type.to_ulong()));
  if(miit==cubitMeshsets.get<Composite_Cubit_msId_And_MeshSetType_mi_tag>().end()) {
    SETERRQ1(PETSC_COMM_SELF,1,"such cubit meshset is already there",ms_id);
  }
  EntityHandle meshset = miit->getMeshSet();
  cubitMeshsets.get<Composite_Cubit_msId_And_MeshSetType_mi_tag>().erase(miit);
  rval = moab.delete_entities(&meshset,1); CHKERRQ_MOAB(rval);
  PetscFunctionReturn(0);
}
PetscErrorCode Core::get_cubit_msId(const int ms_id,const CubitBCType cubit_bc_type,const CubitMeshSets **cubit_meshset_ptr) {
  PetscFunctionBegin;
  CubitMeshSet_multiIndex::index<Composite_Cubit_msId_And_MeshSetType_mi_tag>::type::iterator
    miit = cubitMeshsets.get<Composite_Cubit_msId_And_MeshSetType_mi_tag>().find(boost::make_tuple(ms_id,cubit_bc_type.to_ulong()));
  if(miit!=cubitMeshsets.get<Composite_Cubit_msId_And_MeshSetType_mi_tag>().end()) {
    *cubit_meshset_ptr = &*miit;
  } else {
    SETERRQ1(PETSC_COMM_SELF,1,"msId = %d is not there",ms_id);
  }
  PetscFunctionReturn(0);
}
PetscErrorCode Core::get_cubit_msId_entities_by_dimension(const int msId,const CubitBCType cubit_bc_type,
  const int dimension,Range &entities,const bool recursive) {
  PetscFunctionBegin;
  CubitMeshSet_multiIndex::index<Composite_Cubit_msId_And_MeshSetType_mi_tag>::type::iterator
    miit = cubitMeshsets.get<Composite_Cubit_msId_And_MeshSetType_mi_tag>().find(boost::make_tuple(msId,cubit_bc_type.to_ulong()));
  if(miit!=cubitMeshsets.get<Composite_Cubit_msId_And_MeshSetType_mi_tag>().end()) {
    ierr = miit->get_cubit_msId_entities_by_dimension(moab,dimension,entities,recursive); CHKERRQ(ierr);
  } else {
    SETERRQ1(PETSC_COMM_SELF,1,"msId = %d is not there",msId);
  }
  PetscFunctionReturn(0);
}
PetscErrorCode Core::get_cubit_msId_entities_by_dimension(const int msId,const CubitBCType cubit_bc_type,Range &entities,const bool recursive) {
  PetscFunctionBegin;
  CubitMeshSet_multiIndex::index<Composite_Cubit_msId_And_MeshSetType_mi_tag>::type::iterator
    miit = cubitMeshsets.get<Composite_Cubit_msId_And_MeshSetType_mi_tag>().find(boost::make_tuple(msId,cubit_bc_type.to_ulong()));
  if(miit!=cubitMeshsets.get<Composite_Cubit_msId_And_MeshSetType_mi_tag>().end()) {
    ierr = miit->get_cubit_msId_entities_by_dimension(moab,entities,recursive); CHKERRQ(ierr);
  } else {
    SETERRQ1(PETSC_COMM_SELF,1,"msId = %d is not there",msId);
  }
  PetscFunctionReturn(0);
}
PetscErrorCode Core::get_cubit_msId_entities_by_dimension(const int ms_id,const unsigned int cubit_bc_type,
  const int dimension,Range &entities,const bool recursive) {
  PetscFunctionBegin;
  ierr = get_cubit_msId_entities_by_dimension(ms_id,CubitBCType(cubit_bc_type),dimension,entities,recursive); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
PetscErrorCode Core::get_cubit_msId_entities_by_dimension(const int ms_id,const unsigned int cubit_bc_type,
  Range &entities,const bool recursive) {
  PetscFunctionBegin;
  ierr = get_cubit_msId_entities_by_dimension(ms_id,CubitBCType(cubit_bc_type),entities,recursive); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
PetscErrorCode Core::get_cubit_msId_meshset(const int ms_id,const unsigned int cubit_bc_type,EntityHandle &meshset) {
  PetscFunctionBegin;
  CubitMeshSet_multiIndex::index<Composite_Cubit_msId_And_MeshSetType_mi_tag>::type::iterator
    miit = cubitMeshsets.get<Composite_Cubit_msId_And_MeshSetType_mi_tag>().find(boost::make_tuple(ms_id,cubit_bc_type));
  if(miit!=cubitMeshsets.get<Composite_Cubit_msId_And_MeshSetType_mi_tag>().end()) {
    meshset = miit->meshset;
  } else {
    SETERRQ1(PETSC_COMM_SELF,1,"msId = %d is not there",ms_id);
  }
  PetscFunctionReturn(0);
}
PetscErrorCode Core::get_cubit_meshsets(const unsigned int cubit_bc_type,Range &meshsets) {
  PetscFunctionBegin;
  CubitMeshSet_multiIndex::index<CubitMeshSets_mi_tag>::type::iterator
    miit = cubitMeshsets.get<CubitMeshSets_mi_tag>().lower_bound(cubit_bc_type);
  CubitMeshSet_multiIndex::index<CubitMeshSets_mi_tag>::type::iterator
    hi_miit = cubitMeshsets.get<CubitMeshSets_mi_tag>().upper_bound(cubit_bc_type);
  for(;miit!=hi_miit;miit++) {
    meshsets.insert(miit->meshset);
  }
  PetscFunctionReturn(0);
}

#define SET_BASIC_METHOD(METHOD,PROBLEM_PTR) \
  { \
    METHOD.rAnk = rAnk; \
    METHOD.sIze = sIze; \
    METHOD.problemPtr = PROBLEM_PTR; \
    METHOD.fieldsPtr = &fIelds; \
    METHOD.refinedEntitiesPtr = &refinedEntities; \
    METHOD.entitiesPtr = &entsFields; \
    METHOD.dofsPtr = &dofsField; \
    METHOD.refinedFiniteElementsPtr = &refinedFiniteElements; \
    METHOD.finiteElementsPtr = &finiteElements; \
    METHOD.finiteElementsEntitiesPtr = &entsFiniteElements; \
    METHOD.adjacenciesPtr = &entFEAdjacencies; \
  }

PetscErrorCode Core::problem_basic_method_preProcess(const MoFEMProblem *problem_ptr,BasicMethod &method,int verb) {
  PetscFunctionBegin;
  if(verb==-1) verb = verbose;
  // finite element
  SET_BASIC_METHOD(method,problem_ptr)
  PetscLogEventBegin(USER_EVENT_preProcess,0,0,0,0);
  ierr = method.preProcess(); CHKERRQ(ierr);
  PetscLogEventEnd(USER_EVENT_preProcess,0,0,0,0);
  PetscFunctionReturn(0);
}

PetscErrorCode Core::problem_basic_method_preProcess(const std::string &problem_name,BasicMethod &method,int verb) {
  PetscFunctionBegin;
  if(verb==-1) verb = verbose;
  typedef MoFEMProblem_multiIndex::index<Problem_mi_tag>::type mofem_problems_by_name;
  // find p_miit
  mofem_problems_by_name &pRoblems_set = pRoblems.get<Problem_mi_tag>();
  mofem_problems_by_name::iterator p_miit = pRoblems_set.find(problem_name);
  if(p_miit == pRoblems_set.end()) SETERRQ1(PETSC_COMM_SELF,1,"problem is not in database %s",problem_name.c_str());
  ierr = problem_basic_method_preProcess(&*p_miit,method,verb); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
PetscErrorCode Core::problem_basic_method_postProcess(const MoFEMProblem *problem_ptr,BasicMethod &method,int verb) {
  PetscFunctionBegin;
  SET_BASIC_METHOD(method,problem_ptr)

  PetscLogEventBegin(USER_EVENT_postProcess,0,0,0,0);
  ierr = method.postProcess(); CHKERRQ(ierr);
  PetscLogEventEnd(USER_EVENT_postProcess,0,0,0,0);

  PetscFunctionReturn(0);
}
PetscErrorCode Core::problem_basic_method_postProcess(const std::string &problem_name,BasicMethod &method,int verb) {
  PetscFunctionBegin;
  if(verb==-1) verb = verbose;
  typedef MoFEMProblem_multiIndex::index<Problem_mi_tag>::type mofem_problems_by_name;

  // find p_miit
  mofem_problems_by_name &pRoblems_set = pRoblems.get<Problem_mi_tag>();
  mofem_problems_by_name::iterator p_miit = pRoblems_set.find(problem_name);
  if(p_miit == pRoblems_set.end()) SETERRQ1(PETSC_COMM_SELF,1,"problem is not in database %s",problem_name.c_str());

  ierr = problem_basic_method_postProcess(&*p_miit,method,verb); CHKERRQ(ierr);

  PetscFunctionReturn(0);
}
PetscErrorCode Core::loop_finite_elements(const std::string &problem_name,const std::string &fe_name,FEMethod &method,int verb) {
  PetscFunctionBegin;
  if(verb==-1) verb = verbose;

  ierr = loop_finite_elements(problem_name,fe_name,method,rAnk,rAnk,verb); CHKERRQ(ierr);

  PetscFunctionReturn(0);
}
PetscErrorCode Core::loop_finite_elements(
  const MoFEMProblem *problem_ptr,const std::string &fe_name,FEMethod &method,int lower_rank,int upper_rank,int verb) {
  PetscFunctionBegin;
  if(verb==-1) verb = verbose;
  // finite element

  method.feName = fe_name;
  SET_BASIC_METHOD(method,&*problem_ptr)
  PetscLogEventBegin(USER_EVENT_preProcess,0,0,0,0);
  ierr = method.preProcess(); CHKERRQ(ierr);
  PetscLogEventEnd(USER_EVENT_preProcess,0,0,0,0);

  typedef NumeredEntFiniteElement_multiIndex::index<Composite_Name_And_Part_mi_tag>::type FEByComposite;

  FEByComposite &numeredFiniteElements =
    (const_cast<NumeredEntFiniteElement_multiIndex&>(problem_ptr->numeredFiniteElements)).get<Composite_Name_And_Part_mi_tag>();
  FEByComposite::iterator miit = numeredFiniteElements.lower_bound(boost::make_tuple(fe_name,lower_rank));
  FEByComposite::iterator hi_miit = numeredFiniteElements.upper_bound(boost::make_tuple(fe_name,upper_rank));

  // FEByComposite::iterator back_miit = hi_miit;
  method.loopSize = distance(miit,hi_miit);
  for(int nn = 0;miit!=hi_miit;miit++,nn++) {

    // back_miit--;

    method.nInTheLoop = nn;
    method.numeredEntFiniteElementPtr = &*(*miit);
    method.dataPtr = &((*miit)->sPtr->data_dofs);
    method.rowPtr = (*miit)->rows_dofs.get();
    method.colPtr = (*miit)->cols_dofs.get();

    try {
      PetscLogEventBegin(USER_EVENT_operator,0,0,0,0);
      ierr = method(); CHKERRQ(ierr);
      PetscLogEventEnd(USER_EVENT_operator,0,0,0,0);
    } catch (const std::exception& ex) {
      std::ostringstream ss;
      ss << "throw in method: " << ex.what() << " at line " << __LINE__ << " in file " << __FILE__ << std::endl;
      SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,ss.str().c_str());
    }

  }

  PetscLogEventBegin(USER_EVENT_postProcess,0,0,0,0);
  ierr = method.postProcess(); CHKERRQ(ierr);
  PetscLogEventEnd(USER_EVENT_postProcess,0,0,0,0);

  PetscFunctionReturn(0);
}
PetscErrorCode Core::loop_finite_elements(
  const std::string &problem_name,const std::string &fe_name,FEMethod &method,int lower_rank,int upper_rank,int verb) {
  PetscFunctionBegin;
  if(verb==-1) verb = verbose;
  typedef MoFEMProblem_multiIndex::index<Problem_mi_tag>::type mofem_problems_by_name;
  // find p_miit
  mofem_problems_by_name &pRoblems_set = pRoblems.get<Problem_mi_tag>();
  mofem_problems_by_name::iterator p_miit = pRoblems_set.find(problem_name);
  if(p_miit == pRoblems_set.end()) SETERRQ1(PETSC_COMM_SELF,1,"problem is not in database %s",problem_name.c_str());

  ierr = loop_finite_elements(&*p_miit,fe_name,method,lower_rank,upper_rank,verb); CHKERRQ(ierr);

  PetscFunctionReturn(0);
}
PetscErrorCode Core::loop_dofs(
  const MoFEMProblem *problem_ptr,const std::string &field_name,RowColData rc,EntMethod &method,int lower_rank,int upper_rank,int verb
) {
  PetscFunctionBegin;
  SET_BASIC_METHOD(method,&*problem_ptr);
  typedef NumeredDofEntity_multiIndex::index<Composite_Name_And_Part_mi_tag>::type numerd_dofs;
  numerd_dofs *dofs;
  switch (rc) {
    case ROW:
      dofs = const_cast<numerd_dofs*>(&problem_ptr->numered_dofs_rows->get<Composite_Name_And_Part_mi_tag>());
      break;
    case COL:
      dofs = const_cast<numerd_dofs*>(&problem_ptr->numered_dofs_cols->get<Composite_Name_And_Part_mi_tag>());
      break;
    default:
     SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"not implemented");
  }
  numerd_dofs::iterator miit = dofs->lower_bound(boost::make_tuple(field_name,lower_rank));
  numerd_dofs::iterator hi_miit = dofs->upper_bound(boost::make_tuple(field_name,upper_rank));
  ierr = method.preProcess(); CHKERRQ(ierr);
  for(;miit!=hi_miit;miit++) {
    method.dofPtr = &(*(*miit)->get_DofEntity_ptr());
    method.dofNumeredPtr = &*(*miit);
    try {
      ierr = method(); CHKERRQ(ierr);
    } catch (const std::exception& ex) {
      std::ostringstream ss;
      ss << "throw in method: " << ex.what() << std::endl;
      SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,ss.str().c_str());
    }
  }
  ierr = method.postProcess(); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
PetscErrorCode Core::loop_dofs(
  const std::string &problem_name,const std::string &field_name,RowColData rc,EntMethod &method,int lower_rank,int upper_rank,int verb
) {
  PetscFunctionBegin;
  if(verb==-1) verb = verbose;
  typedef MoFEMProblem_multiIndex::index<Problem_mi_tag>::type mofem_problems_by_name;
  // find p_miit
  mofem_problems_by_name &pRoblems_set = pRoblems.get<Problem_mi_tag>();
  mofem_problems_by_name::iterator p_miit = pRoblems_set.find(problem_name);
  if(p_miit == pRoblems_set.end()) SETERRQ1(PETSC_COMM_SELF,1,"problem not in database %s",problem_name.c_str());
  ierr = loop_dofs(&*p_miit,field_name,rc,method,lower_rank,upper_rank,verb); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
PetscErrorCode Core::loop_dofs(const std::string &problem_name,const std::string &field_name,RowColData rc,EntMethod &method,int verb) {
  PetscFunctionBegin;
  if(verb==-1) verb = verbose;
  ierr = loop_dofs(problem_name,field_name,rc,method,0,sIze,verb); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
PetscErrorCode Core::loop_dofs(const std::string &field_name,EntMethod &method,int verb) {
  PetscFunctionBegin;
  if(verb==-1) verb = verbose;
  SET_BASIC_METHOD(method,NULL);
  DofEntity_multiIndex::index<FieldName_mi_tag>::type::iterator miit,hi_miit;
  miit = dofsField.get<FieldName_mi_tag>().lower_bound(field_name);
  hi_miit = dofsField.get<FieldName_mi_tag>().upper_bound(field_name);
  method.loopSize = distance(miit,hi_miit);
  ierr = method.preProcess(); CHKERRQ(ierr);
  for(int nn = 0; miit!=hi_miit; miit++, nn++) {
    method.nInTheLoop = nn;
    method.dofPtr = &*(*miit);
    method.dofNumeredPtr = NULL;
    try {
      ierr = method(); CHKERRQ(ierr);
    } catch (const std::exception& ex) {
      std::ostringstream ss;
      ss << "throw in method: " << ex.what() << std::endl;
      SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,ss.str().c_str());
    }
  }
  ierr = method.postProcess(); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
PetscErrorCode Core::get_ref_ents(const RefEntity_multiIndex **refined_entities_ptr) const {
  PetscFunctionBegin;
  *refined_entities_ptr = &refinedEntities;
  PetscFunctionReturn(0);
}
PetscErrorCode Core::get_ref_finite_elements(const RefElement_multiIndex **refined_finite_elements_ptr) const {
  PetscFunctionBegin;
  *refined_finite_elements_ptr = &refinedFiniteElements;
  PetscFunctionReturn(0);
}

PetscErrorCode Core::get_problem(const std::string &problem_name,const MoFEMProblem **problem_ptr) const {
  PetscFunctionBegin;
  typedef MoFEMProblem_multiIndex::index<Problem_mi_tag>::type ProblemsByName;
  const ProblemsByName &problems = pRoblems.get<Problem_mi_tag>();
  ProblemsByName::iterator p_miit = problems.find(problem_name);
  if(p_miit == problems.end()) {
    SETERRQ1(
      PETSC_COMM_SELF,MOFEM_OPERATION_UNSUCCESSFUL,
      "problem < %s > not found, (top tip: check spelling)",problem_name.c_str()
    );
  }
  *problem_ptr = &*p_miit;
  PetscFunctionReturn(0);
}
PetscErrorCode Core::get_field_ents(const MoFEMEntity_multiIndex **field_ents) const {
  PetscFunctionBegin;
  *field_ents = &entsFields;
  PetscFunctionReturn(0);
}
PetscErrorCode Core::get_dofs(const DofEntity_multiIndex **dofs_ptr) const {
  PetscFunctionBegin;
  *dofs_ptr = &dofsField;
  PetscFunctionReturn(0);
}
MoFEMEntity_multiIndex::index<FieldName_mi_tag>::type::iterator Core::get_ent_moabfield_by_name_begin(const std::string &field_name) const {
  return entsFields.get<FieldName_mi_tag>().lower_bound(field_name);
}
MoFEMEntity_multiIndex::index<FieldName_mi_tag>::type::iterator Core::get_ent_moabfield_by_name_end(const std::string &field_name) const {
  return entsFields.get<FieldName_mi_tag>().upper_bound(field_name);
}
DofEntity_multiIndex::index<FieldName_mi_tag>::type::iterator Core::get_dofs_by_name_begin(const std::string &field_name) const {
  return dofsField.get<FieldName_mi_tag>().lower_bound(field_name);
}
DofEntity_multiIndex::index<FieldName_mi_tag>::type::iterator Core::get_dofs_by_name_end(const std::string &field_name) const {
  return dofsField.get<FieldName_mi_tag>().upper_bound(field_name);
}
DofEntity_multiIndex::index<Composite_Name_And_Ent_mi_tag>::type::iterator
Core::get_dofs_by_name_and_ent_begin(const std::string &field_name,const EntityHandle ent) const {
  return dofsField.get<Composite_Name_And_Ent_mi_tag>().lower_bound(boost::make_tuple(field_name,ent));
}
DofEntity_multiIndex::index<Composite_Name_And_Ent_mi_tag>::type::iterator
Core::get_dofs_by_name_and_ent_end(const std::string &field_name,const EntityHandle ent) const {
  return dofsField.get<Composite_Name_And_Ent_mi_tag>().upper_bound(boost::make_tuple(field_name,ent));
}
DofEntity_multiIndex::index<Composite_Name_And_Type_mi_tag>::type::iterator
Core::get_dofs_by_name_and_type_begin(const std::string &field_name,const EntityType type) const {
  return dofsField.get<Composite_Name_And_Type_mi_tag>().lower_bound(boost::make_tuple(field_name,type));
}
DofEntity_multiIndex::index<Composite_Name_And_Type_mi_tag>::type::iterator
Core::get_dofs_by_name_and_type_end(const std::string &field_name,const EntityType type) const {
  return dofsField.get<Composite_Name_And_Type_mi_tag>().upper_bound(boost::make_tuple(field_name,type));
}
PetscErrorCode Core::get_finite_elements(const FiniteElement_multiIndex **finiteElements_ptr) const {
  PetscFunctionBegin;
  *finiteElements_ptr = &finiteElements;
  PetscFunctionReturn(0);
}
EntFiniteElement_multiIndex::index<FiniteElement_name_mi_tag>::type::iterator
Core::get_fe_by_name_begin(const std::string &fe_name) const {
  return entsFiniteElements.get<FiniteElement_name_mi_tag>().lower_bound(fe_name);
}
EntFiniteElement_multiIndex::index<FiniteElement_name_mi_tag>::type::iterator
Core::get_fe_by_name_end(const std::string &fe_name) const {
  return entsFiniteElements.get<FiniteElement_name_mi_tag>().upper_bound(fe_name);
}
PetscErrorCode Core::check_number_of_ents_in_ents_field(const std::string& name) const {
  PetscFunctionBegin;
  Field_multiIndex::index<FieldName_mi_tag>::type::iterator it = fIelds.get<FieldName_mi_tag>().find(name);
  if(it == fIelds.get<FieldName_mi_tag>().end()) {
    SETERRQ1(PETSC_COMM_SELF,1,"field not found < %s >",name.c_str());
  }
  EntityHandle meshset = (*it)->getMeshSet();

  int num_entities;
  MoABErrorCode rval;
  rval = moab.get_number_entities_by_handle(meshset,num_entities); CHKERRQ_MOAB(rval);
  if(entsFields.get<FieldName_mi_tag>().count((*it)->getName())
    != (unsigned int)num_entities) {
    SETERRQ1(PETSC_COMM_SELF,1,"not equal number of entities in meshset and field multiindex < %s >",name.c_str());
  }
  PetscFunctionReturn(0);
}
PetscErrorCode Core::check_number_of_ents_in_ents_field() const {
  PetscFunctionBegin;
  Field_multiIndex::index<FieldName_mi_tag>::type::iterator it = fIelds.get<FieldName_mi_tag>().begin();
  for(;it!=fIelds.get<FieldName_mi_tag>().end();it++) {
    if((*it)->getSpace() == NOFIELD) continue; //FIXME: should be treated properly, not test is just skipped for this NOFIELD space
    EntityHandle meshset = (*it)->getMeshSet();
    MoABErrorCode rval;
    int num_entities;
    rval = moab.get_number_entities_by_handle(meshset,num_entities); CHKERRQ_MOAB(rval);
    if(entsFields.get<FieldName_mi_tag>().count((*it)->getName()) != (unsigned int)num_entities) {
      SETERRQ1(PETSC_COMM_SELF,1,"not equal number of entities in meshset and field multiindex < %s >",(*it)->getName().c_str());
    }
  }
  PetscFunctionReturn(0);
}
PetscErrorCode Core::check_number_of_ents_in_ents_finite_element(const std::string& name) const {
  PetscFunctionBegin;
  FiniteElement_multiIndex::index<FiniteElement_name_mi_tag>::type::iterator it;
  it = finiteElements.get<FiniteElement_name_mi_tag>().find(name);
  if(it == finiteElements.get<FiniteElement_name_mi_tag>().end()) {
    SETERRQ1(PETSC_COMM_SELF,1,"finite element not found < %s >",name.c_str());
  }
  EntityHandle meshset = (*it)->getMeshSet();
  MoABErrorCode rval;
  int num_entities;
  rval = moab.get_number_entities_by_handle(meshset,num_entities); CHKERRQ_MOAB(rval);
  if(
    entsFiniteElements.get<FiniteElement_name_mi_tag>().count((*it)->getName().c_str())
    != (unsigned int)num_entities
  ) {
    SETERRQ1(PETSC_COMM_SELF,1,"not equal number of entities in meshset and finite elements multiindex < %s >",(*it)->getName().c_str());
  }
  PetscFunctionReturn(0);
}
PetscErrorCode Core::check_number_of_ents_in_ents_finite_element() const {
  PetscFunctionBegin;
  FiniteElement_multiIndex::index<FiniteElement_name_mi_tag>::type::iterator it;
  it = finiteElements.get<FiniteElement_name_mi_tag>().begin();
  for(;it!=finiteElements.get<FiniteElement_name_mi_tag>().end();it++) {
    EntityHandle meshset = (*it)->getMeshSet();
    MoABErrorCode rval;
    int num_entities;
    rval = moab.get_number_entities_by_handle(meshset,num_entities); CHKERRQ_MOAB(rval);
    if(entsFiniteElements.get<FiniteElement_name_mi_tag>().count((*it)->getName().c_str())
      != (unsigned int)num_entities) {
      SETERRQ1(PETSC_COMM_SELF,1,"not equal number of entities in meshset and finite elements multiindex < %s >",(*it)->getName().c_str());
    }
  }
  PetscFunctionReturn(0);
}

//clear,remove and delete

PetscErrorCode Core::clear_inactive_dofs(int verb) {
  PetscFunctionBegin;
  if(verb==-1) verb = verbose;
  DofEntity_multiIndex::iterator dit;
  dit = dofsField.begin();
  for(;dit!=dofsField.end();dit++) {
    if(!(*dit)->get_active()) {
      MoFEMEntityEntFiniteElementAdjacencyMap_multiIndex::index<Unique_mi_tag>::type::iterator ait,hi_ait;
      ait = entFEAdjacencies.get<Unique_mi_tag>().lower_bound((*dit)->getMoFEMEntityPtr()->getGlobalUniqueId());
      hi_ait = entFEAdjacencies.get<Unique_mi_tag>().upper_bound((*dit)->getMoFEMEntityPtr()->getGlobalUniqueId());
      for(;ait!=hi_ait;ait++) {
        boost::shared_ptr<EntFiniteElement> ent_fe_ptr;
        ent_fe_ptr = ait->entFePtr;
        ent_fe_ptr->row_dof_view->erase((*dit)->getGlobalUniqueId());
        if(ent_fe_ptr->row_dof_view!=ent_fe_ptr->col_dof_view) {
          ent_fe_ptr->col_dof_view->erase((*dit)->getGlobalUniqueId());
        }
        if(
          ent_fe_ptr->row_dof_view!=ent_fe_ptr->data_dof_view||
          ent_fe_ptr->col_dof_view!=ent_fe_ptr->data_dof_view
        ) {
          ent_fe_ptr->data_dof_view->erase((*dit)->getGlobalUniqueId());
        }
        ent_fe_ptr->data_dofs.get<Unique_mi_tag>().erase((*dit)->getGlobalUniqueId());
      }
      dit = dofsField.erase(dit);
      if(dit==dofsField.end()) break;
    }
  }
  PetscFunctionReturn(0);
}

PetscErrorCode Core::clear_dofs_fields(const BitRefLevel &bit,const BitRefLevel &mask,int verb) {
  PetscFunctionBegin;
  if(verb==-1) verb = verbose;
  {
    DofEntity_multiIndex::iterator dit;
    dit = dofsField.begin();
    for(;dit!=dofsField.end();) {
      BitRefLevel bit2 = (*dit)->getBitRefLevel();
      if((*dit)->getEntType()==MBENTITYSET) {
        dit++;
        continue;
      }
      if((bit2&mask)!=bit2) {
        dit++;
        continue;
      }
      if((bit2&bit).none()) {
        dit++;
        continue;
      }
      MoFEMEntityEntFiniteElementAdjacencyMap_multiIndex::index<Unique_mi_tag>::type::iterator ait,hi_ait;
      ait = entFEAdjacencies.get<Unique_mi_tag>().lower_bound((*dit)->getMoFEMEntityPtr()->getGlobalUniqueId());
      hi_ait = entFEAdjacencies.get<Unique_mi_tag>().upper_bound((*dit)->getMoFEMEntityPtr()->getGlobalUniqueId());
      for(;ait!=hi_ait;ait++) {
        boost::shared_ptr<EntFiniteElement> ent_fe_ptr;
        ent_fe_ptr = ait->entFePtr;
        ent_fe_ptr->row_dof_view->erase((*dit)->getGlobalUniqueId());
        if(ent_fe_ptr->row_dof_view!=ent_fe_ptr->col_dof_view) {
          ent_fe_ptr->col_dof_view->erase((*dit)->getGlobalUniqueId());
        }
        if(
          ent_fe_ptr->row_dof_view!=ent_fe_ptr->data_dof_view||
          ent_fe_ptr->col_dof_view!=ent_fe_ptr->data_dof_view
        ) {
          ent_fe_ptr->data_dof_view->erase((*dit)->getGlobalUniqueId());
        }
        ent_fe_ptr->data_dofs.get<Unique_mi_tag>().erase((*dit)->getGlobalUniqueId());
      }
      dit = dofsField.erase(dit);
    }
  }
  PetscFunctionReturn(0);
}
PetscErrorCode Core::clear_dofs_fields(const std::string &name,const Range ents,int verb) {
  PetscFunctionBegin;
  if(verb==-1) verb = verbose;
  Range::iterator eit = ents.begin();
  for(;eit!=ents.end();eit++) {
    DofEntity_multiIndex::index<Composite_Name_And_Ent_mi_tag>::type::iterator dit,hi_dit;
    dit = dofsField.get<Composite_Name_And_Ent_mi_tag>().lower_bound(boost::make_tuple(name,*eit));
    hi_dit = dofsField.get<Composite_Name_And_Ent_mi_tag>().upper_bound(boost::make_tuple(name,*eit));
    for(;dit!=hi_dit;) {
      MoFEMEntityEntFiniteElementAdjacencyMap_multiIndex::index<Unique_mi_tag>::type::iterator ait,hi_ait;
      ait = entFEAdjacencies.get<Unique_mi_tag>().lower_bound((*dit)->getMoFEMEntityPtr()->getGlobalUniqueId());
      hi_ait = entFEAdjacencies.get<Unique_mi_tag>().upper_bound((*dit)->getMoFEMEntityPtr()->getGlobalUniqueId());
      for(;ait!=hi_ait;ait++) {
        boost::shared_ptr<EntFiniteElement> ent_fe_ptr;
        ent_fe_ptr = ait->entFePtr;
        ent_fe_ptr->row_dof_view->erase((*dit)->getGlobalUniqueId());
        if(ent_fe_ptr->row_dof_view!=ent_fe_ptr->col_dof_view) {
          ent_fe_ptr->col_dof_view->erase((*dit)->getGlobalUniqueId());
        }
        if(
          ent_fe_ptr->row_dof_view!=ent_fe_ptr->data_dof_view||
          ent_fe_ptr->col_dof_view!=ent_fe_ptr->data_dof_view
        ) {
          ent_fe_ptr->data_dof_view->erase((*dit)->getGlobalUniqueId());
        }
        ent_fe_ptr->data_dofs.get<Unique_mi_tag>().erase((*dit)->getGlobalUniqueId());
      }
      dit = dofsField.get<Composite_Name_And_Ent_mi_tag>().erase(dit);
    }
  }
  PetscFunctionReturn(0);
}
PetscErrorCode Core::clear_ents_fields(const BitRefLevel &bit,const BitRefLevel &mask,int verb) {
  PetscFunctionBegin;
  if(verb==-1) verb = verbose;
  ierr = clear_dofs_fields(bit,mask,verb); CHKERRQ(ierr);
  ierr = clear_adjacencies_entities(bit,mask,verb); CHKERRQ(ierr);
  MoFEMEntity_multiIndex::iterator eit;
  eit = entsFields.begin();
  for(;eit!=entsFields.end();) {
    if((*eit)->getEntType()==MBENTITYSET) {
      eit++;
      continue;
    }
    BitRefLevel bit2 = (*eit)->getBitRefLevel();
    if((bit2&mask)!=bit2) {
      eit++;
      continue;
    }
    if((bit2&bit).none()) {
      eit++;
      continue;
    }
    EntityHandle ent = (*eit)->getEnt();
    rval = moab.tag_delete_data((*eit)->sFieldPtr->th_AppOrder,&ent,1); CHKERRQ_MOAB(rval);
    if((*eit)->tag_FieldData_size>0) {
      rval = moab.tag_delete_data((*eit)->sFieldPtr->th_FieldData,&ent,1); CHKERRQ_MOAB(rval);
    }
    eit = entsFields.erase(eit);
  }
  PetscFunctionReturn(0);
}
PetscErrorCode Core::clear_ents_fields(const std::string &name,const Range ents,int verb) {
  PetscFunctionBegin;
  if(verb==-1) verb = verbose;
  ierr = clear_dofs_fields(name,ents,verb); CHKERRQ(ierr);
  ierr = clear_adjacencies_entities(name,ents,verb); CHKERRQ(ierr);
  Range::iterator eit = ents.begin();
  for(;eit!=ents.end();eit++) {
    MoFEMEntity_multiIndex::index<Composite_Name_And_Ent_mi_tag>::type::iterator dit,hi_dit;
    dit = entsFields.get<Composite_Name_And_Ent_mi_tag>().lower_bound(boost::make_tuple(name,*eit));
    hi_dit = entsFields.get<Composite_Name_And_Ent_mi_tag>().upper_bound(boost::make_tuple(name,*eit));
    for(;dit!=hi_dit;) {
      dit = entsFields.get<Composite_Name_And_Ent_mi_tag>().erase(dit);
    }
  }
  PetscFunctionReturn(0);
}
PetscErrorCode Core::clear_finite_elements(const BitRefLevel &bit,const BitRefLevel &mask,int verb) {
  PetscFunctionBegin;
  if(verb==-1) verb = verbose;
  ierr = clear_adjacencies_finite_elements(bit,mask,verb); CHKERRQ(ierr);
  EntFiniteElement_multiIndex::iterator fe_it = entsFiniteElements.begin();
  for(;fe_it!=entsFiniteElements.end();) {
    BitRefLevel bit2 = (*fe_it)->getBitRefLevel();
    if((*fe_it)->getEntType()==MBENTITYSET) {
      fe_it++;
      continue;
    }
    if((bit2&mask)!=bit2) {
      fe_it++;
      continue;
    }
    if((bit2&bit).none()) {
      fe_it++;
      continue;
    }
    fe_it = entsFiniteElements.erase(fe_it);
  }
  PetscFunctionReturn(0);
}
PetscErrorCode Core::clear_finite_elements(const std::string &name,const Range &ents,int verb) {
  PetscFunctionBegin;
  if(verb==-1) verb = verbose;
  ierr = clear_adjacencies_finite_elements(name,ents,verb); CHKERRQ(ierr);
  Range::iterator eit = ents.begin();
  for(;eit!=ents.end();eit++) {
    EntFiniteElement_multiIndex::index<Composite_Name_And_Ent_mi_tag>::type::iterator fit,hi_fit;
    fit = entsFiniteElements.get<Composite_Name_And_Ent_mi_tag>().lower_bound(boost::make_tuple(name,*eit));
    hi_fit = entsFiniteElements.get<Composite_Name_And_Ent_mi_tag>().upper_bound(boost::make_tuple(name,*eit));
    for(;fit!=hi_fit;) {
      fit = entsFiniteElements.get<Composite_Name_And_Ent_mi_tag>().erase(fit);
    }
  }
  PetscFunctionReturn(0);
}
PetscErrorCode Core::clear_adjacencies_finite_elements(const BitRefLevel &bit,const BitRefLevel &mask,int verb) {
  PetscFunctionBegin;
  MoFEMEntityEntFiniteElementAdjacencyMap_multiIndex::iterator ait;
  ait = entFEAdjacencies.begin();
  for(;ait!=entFEAdjacencies.end();) {
    BitRefLevel bit2 = ait->entFePtr->getBitRefLevel();
    if(ait->entFePtr->getEntType()==MBENTITYSET) {
      ait++;
      continue;
    }
    if((bit2&mask)!=bit2) {
      ait++;
      continue;
    }
    if((bit2&bit).none()) {
      ait++;
      continue;
    }
    ait = entFEAdjacencies.erase(ait);
  }
  PetscFunctionReturn(0);
}
PetscErrorCode Core::clear_adjacencies_finite_elements(const std::string &name,const Range &ents,int verb) {
  PetscFunctionBegin;
  if(verb==-1) verb = verbose;
  Range::iterator eit = ents.begin();
  for(;eit!=ents.end();eit++) {
    MoFEMEntityEntFiniteElementAdjacencyMap_multiIndex::index<FEEnt_mi_tag>::type::iterator ait,hi_ait;
    ait = entFEAdjacencies.get<FEEnt_mi_tag>().lower_bound(*eit);
    hi_ait = entFEAdjacencies.get<FEEnt_mi_tag>().upper_bound(*eit);
    for(;ait!=hi_ait;) {
      if(ait->entFePtr->getName() == name) {
        ait = entFEAdjacencies.get<FEEnt_mi_tag>().erase(ait);
      } else {
        ait++;
      }
    }
  }
  PetscFunctionReturn(0);
}
PetscErrorCode Core::clear_adjacencies_entities(const BitRefLevel &bit,const BitRefLevel &mask,int verb ) {
  PetscFunctionBegin;
  if(verb==-1) verb = verbose;
  MoFEMEntityEntFiniteElementAdjacencyMap_multiIndex::iterator ait;
  ait = entFEAdjacencies.begin();
  for(;ait!=entFEAdjacencies.end();) {
    BitRefLevel bit2 = ait->mofemEntPtr->getBitRefLevel();
    if(ait->mofemEntPtr->getEntType()==MBENTITYSET) {
      ait++;
      continue;
    }
    if((bit2&mask)!=bit2) {
      ait++;
      continue;
    }
    if((bit2&bit).none()) {
      ait++;
      continue;
    }
    ait = entFEAdjacencies.erase(ait);
  }
  PetscFunctionReturn(0);
}
PetscErrorCode Core::clear_adjacencies_entities(const std::string &name,const Range &ents,int verb) {
  PetscFunctionBegin;
  if(verb==-1) verb = verbose;
  Range::iterator eit = ents.begin();
  for(;eit!=ents.end();eit++) {
    MoFEMEntityEntFiniteElementAdjacencyMap_multiIndex::index<Ent_mi_tag>::type::iterator ait,hi_ait;
    ait = entFEAdjacencies.get<Ent_mi_tag>().lower_bound(*eit);
    hi_ait = entFEAdjacencies.get<Ent_mi_tag>().upper_bound(*eit);
    for(;ait!=hi_ait;) {
      if(ait->mofemEntPtr->getName() == name) {
        ait = entFEAdjacencies.get<Ent_mi_tag>().erase(ait);
      } else {
        ait++;
      }
    }
  }
  PetscFunctionReturn(0);
}
PetscErrorCode Core::remove_ents_from_field_by_bit_ref(const BitRefLevel &bit,const BitRefLevel &mask,int verb) {
  PetscFunctionBegin;
  if(verb==-1) verb = verbose;
  ierr = clear_ents_fields(bit,mask,verb); CHKERRQ(ierr);
  Field_multiIndex::iterator f_it = fIelds.begin();
  for(;f_it!=fIelds.end();f_it++) {
    EntityHandle meshset = (*f_it)->getMeshSet();
    Range ents_to_remove;
    rval = moab.get_entities_by_handle(
      meshset,ents_to_remove,false); CHKERRQ_MOAB(rval);
    Range::iterator eit = ents_to_remove.begin();
    for(;eit!=ents_to_remove.end();) {
      if(moab.type_from_handle(*eit)==MBENTITYSET) {
        eit = ents_to_remove.erase(eit);
        continue;
      }
      BitRefLevel bit2;
      rval = moab.tag_get_data(th_RefBitLevel,&*eit,1,&bit2); CHKERRQ_MOAB(rval);
      if((bit2&mask)!=bit2) {
        eit = ents_to_remove.erase(eit);
        continue;
      }
      if((bit2&bit).none()) {
        eit = ents_to_remove.erase(eit);
        continue;
      }
      MoFEMEntity_multiIndex::index<Composite_Name_And_Ent_mi_tag>::type::iterator iit;
      iit = entsFields.get<Composite_Name_And_Ent_mi_tag>().find(boost::make_tuple((*f_it)->getName(),*eit));
      if(iit != entsFields.get<Composite_Name_And_Ent_mi_tag>().end()) {
        SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"data inconsistency");
      }
      eit++;
    }
    rval = moab.remove_entities(meshset,ents_to_remove); CHKERRQ_MOAB(rval);
    if(verb>0) {
      PetscPrintf(comm,
        "number of removed entities = %u from field %s\n",
        ents_to_remove.size(),
        (*f_it)->getName().c_str()
      );
      if(verb>1) {
        int num_entities;
        rval = moab.get_number_entities_by_handle(meshset,num_entities); CHKERRQ_MOAB(rval);
        PetscPrintf(comm,"\tnumber of entities in database = %u and meshset = %u\n",
        entsFields.get<BitFieldId_mi_tag>().count((*f_it)->getId()),num_entities);
      }
    }
  }
  PetscFunctionReturn(0);
}
PetscErrorCode Core::remove_ents_from_field(const std::string& name,const EntityHandle meshset,const EntityType type,int verb) {
  PetscFunctionBegin;
  if(verb==-1) verb = verbose;
  Range ents;
  rval = moab.get_entities_by_type(meshset,type,ents); CHKERRQ_MOAB(rval);
  ierr = remove_ents_from_field(name,ents,verb); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
PetscErrorCode Core::remove_ents_from_field(const std::string& name,const Range &ents,int verb) {
  PetscFunctionBegin;
  if(verb==-1) verb = verbose;
  EntityHandle meshset;
  try {
    meshset = get_field_meshset(name);
  } catch (MoFEMException const &e) {
    SETERRQ(PETSC_COMM_SELF,e.errorCode,e.errorMessage);
  }
  rval = moab.remove_entities(meshset,ents); CHKERRQ_MOAB(rval);
  ierr = clear_ents_fields(name,ents,verb); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
PetscErrorCode Core::remove_ents_from_finite_element_by_bit_ref(const BitRefLevel &bit,const BitRefLevel &mask,int verb) {
  PetscFunctionBegin;
  if(verb==-1) verb = verbose;
  ierr = clear_finite_elements(bit,mask,verb); CHKERRQ(ierr);
  FiniteElement_multiIndex::iterator fe_it = finiteElements.begin();
  for(;fe_it!=finiteElements.end();fe_it++) {
    EntityHandle meshset = (*fe_it)->getMeshSet();
    Range ents_to_remove;
    rval = moab.get_entities_by_handle(
      meshset,ents_to_remove,false
    ); CHKERRQ_MOAB(rval);
    Range::iterator eit = ents_to_remove.begin();
    for(;eit!=ents_to_remove.end();) {
      if(moab.type_from_handle(*eit)==MBENTITYSET) {
        eit = ents_to_remove.erase(eit);
        continue;
      }
      BitRefLevel bit2;
      rval = moab.tag_get_data(th_RefBitLevel,&*eit,1,&bit2); CHKERRQ_MOAB(rval);
      if((bit2&mask)!=bit2) {
        eit = ents_to_remove.erase(eit);
        continue;
      }
      if((bit2&bit).none()) {
        eit = ents_to_remove.erase(eit);
        continue;
      }
      EntFiniteElement_multiIndex::index<Composite_Name_And_Ent_mi_tag>::type::iterator iit;
      iit = entsFiniteElements.get<Composite_Name_And_Ent_mi_tag>().find(
        boost::make_tuple((*fe_it)->getName(),*eit)
      );
      if(iit != entsFiniteElements.get<Composite_Name_And_Ent_mi_tag>().end()) {
        SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"data inconsistency");
      }
      eit++;
    }
    rval = moab.remove_entities(meshset,ents_to_remove); CHKERRQ_MOAB(rval);
    if(verb>0) {
      PetscPrintf(comm,
        "number of removed entities = %u from finite element %s\n",
        ents_to_remove.size(),(*fe_it)->getName().c_str()
      );
    }
  }
  PetscFunctionReturn(0);
}
PetscErrorCode Core::remove_ents_from_finite_element(const std::string &name,const EntityHandle meshset,const EntityType type,int verb) {
  PetscFunctionBegin;
  if(verb==-1) verb = verbose;
  Range ents;
  rval = moab.get_entities_by_type(meshset,type,ents,false); CHKERRQ_MOAB(rval);
  ierr = remove_ents_from_finite_element(name,ents,verb); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
PetscErrorCode Core::remove_ents_from_finite_element(const std::string &name,const Range &ents,int verb) {
  PetscFunctionBegin;
  if(verb==-1) verb = verbose;
  ierr = clear_finite_elements(name,ents,verb); CHKERRQ(ierr);
  const BitFEId id = get_BitFEId(name);
  const EntityHandle idm = get_finite_element_meshset(id);
  rval = moab.remove_entities(idm,ents); CHKERRQ_MOAB(rval);
  PetscFunctionReturn(0);
}
PetscErrorCode Core::remove_ents_by_bit_ref(const BitRefLevel &bit,const BitRefLevel &mask,int verb) {
  PetscFunctionBegin;
  if(verb==-1) verb = verbose;
  ierr = delete_finite_elements_by_bit_ref(bit,mask,verb); CHKERRQ(ierr);
  ierr = remove_ents_from_field_by_bit_ref(bit,mask,verb); CHKERRQ(ierr);
  RefEntity_multiIndex::iterator ent_it = refinedEntities.begin();
  for(;ent_it!=refinedEntities.end();) {
    BitRefLevel bit2 = (*ent_it)->getBitRefLevel();
    if((*ent_it)->getEntType()==MBENTITYSET) {
      ent_it++;
      continue;
    }
    if((bit2&mask)!=bit2) {
      ent_it++;
      continue;
    }
    if((bit2&bit).none()) {
      ent_it++;
      continue;
    }
    ent_it = refinedEntities.erase(ent_it);
  }
  PetscFunctionReturn(0);
}
PetscErrorCode Core::delete_ents_by_bit_ref(
  const BitRefLevel &bit,const BitRefLevel &mask,const bool remove_parent,int verb
) {
  PetscFunctionBegin;
  Range ents_to_delete;
  rval = moab.get_entities_by_handle(0,ents_to_delete,false); CHKERRQ_MOAB(rval);
  {
    Range::iterator eit = ents_to_delete.begin();
    for(;eit!=ents_to_delete.end();) {
      if(moab.type_from_handle(*eit)==MBENTITYSET) {
        eit = ents_to_delete.erase(eit);
        continue;
      }
      BitRefLevel bit2;
      rval = moab.tag_get_data(th_RefBitLevel,&*eit,1,&bit2); CHKERRQ_MOAB(rval);
      if((bit2&mask)!=bit2) {
        eit = ents_to_delete.erase(eit);
        continue;
      }
      if((bit2&bit).none()) {
        eit = ents_to_delete.erase(eit);
        continue;
      }
      eit++;
    }
  }
  if(remove_parent) { //remove parent
    Range::iterator eit = ents_to_delete.begin();
    for(;eit != ents_to_delete.end();eit++) {
      RefEntity_multiIndex::index<Ent_Ent_mi_tag>::type::iterator pit,hi_pit;
      pit = refinedEntities.get<Ent_Ent_mi_tag>().lower_bound(*eit);
      hi_pit = refinedEntities.get<Ent_Ent_mi_tag>().upper_bound(*eit);
      for(;pit!=hi_pit;pit++) {
        EntityHandle ent = (*pit)->getRefEnt();
        if(ents_to_delete.find(ent) != ents_to_delete.end()) {
          continue;
        }
        /*if(rAnk==0) {
        EntityHandle out_meshset;
        rval = moab.create_meshset(MESHSET_SET,out_meshset); CHKERRQ_MOAB(rval);
        rval = moab.add_entities(out_meshset,&ent,1); CHKERRQ_MOAB(rval);
        rval = moab.add_entities(out_meshset,&*eit,1); CHKERRQ_MOAB(rval);
        rval = moab.write_file("error.vtk","VTK","",&out_meshset,1); CHKERRQ_MOAB(rval);
        rval = moab.delete_entities(&out_meshset,1); CHKERRQ_MOAB(rval);
        }
        std::ostringstream ss;
        ss << "child:\n" << *pit << std::endl;
        ss << "parent:\n" << RefEntity(moab,*eit) << std::endl;
        SETERRQ1(PETSC_COMM_SELF,1,
        "entity can not be removed, it is parent for some other entity\n%s",ss.str().c_str());*/
        bool success = refinedEntities.modify(
          refinedEntities.project<0>(pit),RefEntity_change_remove_parent()
        );
        if(!success) {
          SETERRQ(PETSC_COMM_SELF,MOFEM_OPERATION_UNSUCCESSFUL,"mofification unsucessfull");
        }
      }
    }
  }
  { //remove deleted entities form cubit meshsets
    CubitMeshSet_multiIndex::iterator cubit_it;
    cubit_it = cubitMeshsets.begin();
    for(;cubit_it!=cubitMeshsets.end();cubit_it++) {
      EntityHandle cubit_meshset = cubit_it->meshset;
      rval = moab.remove_entities(cubit_meshset,ents_to_delete); CHKERRQ_MOAB(rval);
      Range meshsets;
      rval = moab.get_entities_by_type(cubit_meshset,MBENTITYSET,meshsets);  CHKERRQ_MOAB(rval);
      for(Range::iterator mit = meshsets.begin();mit!=meshsets.end();mit++) {
        rval = moab.remove_entities(*mit,ents_to_delete); CHKERRQ_MOAB(rval);
      }
    }
  }
  ierr = remove_ents_by_bit_ref(bit,mask,verb); CHKERRQ(ierr);
  if(verb>0) {
    PetscPrintf(comm,"number of deleted entities = %u\n",ents_to_delete.size());
  }
  if(verb>2) {
    EntityHandle out_meshset;
    rval = moab.create_meshset(MESHSET_SET,out_meshset); CHKERRQ_MOAB(rval);
    rval = moab.add_entities(out_meshset,ents_to_delete.subset_by_type(MBTET)); CHKERRQ_MOAB(rval);
    rval = moab.write_file("debug_ents_to_delete.vtk","VTK","",&out_meshset,1); CHKERRQ_MOAB(rval);
    rval = moab.delete_entities(&out_meshset,1); CHKERRQ_MOAB(rval);
  }
  //delete entities form moab
  for(int dd = 3;dd>=0;dd--) {
    rval = moab.delete_entities(ents_to_delete.subset_by_dimension(dd)); //CHKERRQ_MOAB(rval);
  }
  PetscFunctionReturn(0);
}
PetscErrorCode Core::delete_finite_elements_by_bit_ref(
  const BitRefLevel &bit,const BitRefLevel &mask,int verb
) {
  PetscFunctionBegin;
  if(verb==-1) verb = verbose;
  ierr = remove_ents_from_finite_element_by_bit_ref(bit,mask,verb); CHKERRQ(ierr);
  RefElement_multiIndex::iterator fe_it = refinedFiniteElements.begin();
  for(;fe_it!=refinedFiniteElements.end();) {
    BitRefLevel bit2 = fe_it->getBitRefLevel();
    if(fe_it->getEntType()==MBENTITYSET) {
      fe_it++;
      continue;
    }
    if((bit2&mask)!=bit2) {
      fe_it++;
      continue;
    }
    if((bit2&bit).none()) {
      fe_it++;
      continue;
    }
    fe_it = refinedFiniteElements.erase(fe_it);
  }
  PetscFunctionReturn(0);
}
PetscErrorCode Core::delete_finite_element(const std::string name,int verb) {
  PetscFunctionBegin;
  typedef FiniteElement_multiIndex::index<FiniteElement_name_mi_tag>::type FiniteElements_by_name;
  FiniteElements_by_name& fe = finiteElements.get<FiniteElement_name_mi_tag>();
  FiniteElements_by_name::iterator miit = fe.find(name);
  if(miit==fe.end()) {
    SETERRQ1(
      PETSC_COMM_SELF,MOFEM_OPERATION_UNSUCCESSFUL,
      "finite element <%s> not found",name.c_str()
    );
  }
  EntityHandle meshset = (*miit)->getMeshSet();
  Range ents;
  rval = moab.get_entities_by_handle(meshset,ents,false); CHKERRQ_MOAB(rval);
  ierr = remove_ents_from_finite_element(name,ents,verb); CHKERRQ(ierr);
  fe.erase(miit);
  rval = moab.delete_entities(&meshset,1); CHKERRQ_MOAB(rval);
  PetscFunctionReturn(0);
}


}
