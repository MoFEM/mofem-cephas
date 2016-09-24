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
#include <MeshRefinment.hpp>
#include <PrismInterface.hpp>
#include <SeriesRecorder.hpp>
#include <Core.hpp>

namespace MoFEM {


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
      EntityHandle coord_sys_id = undefined_cs_ptr->getMeshset();
      rval = moab.tag_set_data(
        cs_manger_ptr->get_th_CoordSysMeshset(),&meshset,1,&coord_sys_id
      ); CHKERRQ_MOAB(rval);
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

    // You would add h-div space on edge here if you would have 2d hdiv space. At the
    // moment attention is focussed on 3d problems.

    SETERRQ(comm,MOFEM_NOT_IMPLEMENTED,"sorry, not implemented, HCURL not implemented for triangle");

    break;
    case HDIV:

    // You would add h-div space on edge here if you would have 2d hdiv space. At the MOFEM_NOT_IMPLEMENTED
    // attention is focussed on 3d problems.

    SETERRQ(comm,MOFEM_NOT_IMPLEMENTED,"sorry, not implemented, HDIV not implemented for triangle");
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
  typedef Field_multiIndex::index<BitFieldId_mi_tag>::type field_set_by_id;
  const field_set_by_id &set_id = fIelds.get<BitFieldId_mi_tag>();
  field_set_by_id::iterator miit = set_id.find(id);
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
          SETERRQ(comm,MOFEM_DATA_INCONSISTENCY,"approximation order for H1 space and vertex different than 1 makes not sense");
        }
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
      RefEntity_multiIndex::index<Ent_mi_tag>::type::iterator miit_ref_ent = refinedEntities.get<Ent_mi_tag>().find(*eit);
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
PetscErrorCode Core::set_field_order_by_entity_type_and_bit_ref(const BitRefLevel &bit,const BitRefLevel &mask,const EntityType type,const BitFieldId id,const ApproximationOrder order,int verb) {
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
PetscErrorCode Core::set_field_order_by_entity_type_and_bit_ref(const BitRefLevel &bit,const BitRefLevel &mask,const EntityType type,const std::string& name,const ApproximationOrder order,int verb) {
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
PetscErrorCode Core::dofs_NoField(const BitFieldId id,std::map<EntityType,int> &dof_counter,int verb) {
  PetscFunctionBegin;
  if(verb==-1) verb = verbose;
  //field it
  typedef Field_multiIndex::index<BitFieldId_mi_tag>::type field_set_by_id;
  const field_set_by_id &set_id = fIelds.get<BitFieldId_mi_tag>();
  //find fiels
  field_set_by_id::iterator miit = set_id.find(id);
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
    // create mofem entity linked to ref ent
    MoFEMEntity_multiIndex::iterator e_miit;
    try {
      e_miit = entsFields.find(MoFEMEntity(*miit,*miit_ref_ent).getGlobalUniqueId());
    } catch (MoFEMException const &e) {
      SETERRQ(comm,e.errorCode,e.errorMessage);
    } catch (const std::exception& ex) {
      std::ostringstream ss;
      ss << ex.what() << std::endl;
      SETERRQ(comm,MOFEM_OPERATION_UNSUCCESSFUL,ss.str().c_str());
    }
    if(e_miit == entsFields.end()) {
      ApproximationOrder order = -1;
      rval = moab.tag_set_data((*miit)->th_AppOrder,&*eit,1,&order); CHKERRQ_MOAB(rval);
      std::pair<MoFEMEntity_multiIndex::iterator,bool> p_e_miit;
      try {
        boost::shared_ptr<MoFEMEntity> moabent(new MoFEMEntity(*miit,*miit_ref_ent));
        p_e_miit = entsFields.insert(moabent);
      } catch (MoFEMException const &e) {
        SETERRQ(comm,e.errorCode,e.errorMessage);
      } catch (const std::exception& ex) {
        std::ostringstream ss;
        ss << ex.what() << std::endl;
        SETERRQ(comm,MOFEM_OPERATION_UNSUCCESSFUL,ss.str().c_str());
      }
      if(!p_e_miit.second) SETERRQ(comm,MOFEM_DATA_INCONSISTENCY,"data inconsistency");
      bool success = entsFields.modify(p_e_miit.first,MoFEMEntity_change_order(-1));
      if(!success) SETERRQ(comm,MOFEM_OPERATION_UNSUCCESSFUL,"modification unsuccessful");
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
              if(!success) SETERRQ(comm,MOFEM_OPERATION_UNSUCCESSFUL,"modification unsuccessful");
            } else {
              if(DD<nb_active_dosf_on_ent) {
              } else {
                if((*d_miit.first)->getActive()) {
                  is_active = false;
                  inactive_dof_counter[(*d_miit.first)->getEntType()]++;
                  bool success = dofsField.modify(d_miit.first,DofEntity_active_change(is_active));
                  if(!success) SETERRQ(comm,MOFEM_OPERATION_UNSUCCESSFUL,"modification unsuccessful");
                }
              }
            }
            //check ent
            if((*d_miit.first)->getEnt()!=(*e_miit)->getEnt()) {
              SETERRQ(comm,MOFEM_DATA_INCONSISTENCY,"data inconsistency");
            }
            if((*d_miit.first)->getEntType()!=(*e_miit)->getEntType()) {
              SETERRQ(comm,MOFEM_DATA_INCONSISTENCY,"data inconsistency");
            }
            if((*d_miit.first)->getId()!=(*e_miit)->getId()) {
              SETERRQ(comm,MOFEM_DATA_INCONSISTENCY,"data inconsistency");
            }
            //check dof
            if((*d_miit.first)->getDofOrder()!=oo) {
              std::ostringstream ss;
              ss << "data inconsistency!" << std::endl;
              ss << "should be " << mdof << std::endl;
              ss << "but is " << *(*d_miit.first) << std::endl;
              SETERRQ(comm,MOFEM_DATA_INCONSISTENCY,ss.str().c_str());
            }
            if((*d_miit.first)->getMaxOrder()!=(*e_miit)->getMaxOrder()) {
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
      SETERRQ(comm,MOFEM_NOT_IMPLEMENTED,"not implemented");
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
          SETERRQ(comm,MOFEM_NOT_IMPLEMENTED,"not implemented");
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

  if(!(*buildMoFEM&BUILD_FIELD)) SETERRQ(comm,MOFEM_DATA_INCONSISTENCY,"fields not build");
  if(!(*buildMoFEM&BUILD_FE)) SETERRQ(comm,MOFEM_DATA_INCONSISTENCY,"FEs not build");
  if(!(*buildMoFEM&BUILD_ADJ)) SETERRQ(comm,MOFEM_DATA_INCONSISTENCY,"adjacencies not build");
  if(!(*buildMoFEM&BUILD_PROBLEM)) SETERRQ(comm,MOFEM_DATA_INCONSISTENCY,"problem not build");
  if(!(*buildMoFEM&PARTITION_PROBLEM)) SETERRQ(comm,MOFEM_DATA_INCONSISTENCY,"problem not partitioned");

  if(low_proc == -1) low_proc = rAnk;
  if(hi_proc == -1) hi_proc = rAnk;
  typedef MoFEMProblem_multiIndex::index<Problem_mi_tag>::type mofem_problems_by_name;
  //find p_miit
  mofem_problems_by_name &pRoblems_set = pRoblems.get<Problem_mi_tag>();
  mofem_problems_by_name::iterator p_miit = pRoblems_set.find(name);
  if(p_miit == pRoblems_set.end()) {
    SETERRQ1(
      comm,1,"problem < %s > not found (top tip: check spelling)",name.c_str()
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
    if(((*miit2)->getId()&p_miit->getBitFEId()).none()) continue; // if element is not part of problem
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
        ierr = (*miit2)->getRowDofView(
            *(p_miit->numered_dofs_rows),rows_view,moab::Interface::UNION
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
        (numered_fe->getPart()>=(unsigned int)low_proc)&&
        (numered_fe->getPart()<=(unsigned int)hi_proc)
      ) {
        if(part_from_moab) {
          //rows_view
          ierr = (*miit2)->getRowDofView(
            *(p_miit->numered_dofs_rows),rows_view,moab::Interface::UNION
          ); CHKERRQ(ierr);
        }
        //rows element dof multi-indices
        viit_rows = rows_view.begin();
        for(;viit_rows!=rows_view.end();viit_rows++) {
          try {
            boost::shared_ptr<SideNumber> side_number_ptr = (*miit2)->getSideNumberPtr(moab,(*viit_rows)->getEnt());
            rows_dofs->insert(boost::shared_ptr<FENumeredDofEntity>(new FENumeredDofEntity(side_number_ptr,*viit_rows)));
          } catch (MoFEMException const &e) {
            SETERRQ(comm,e.errorCode,e.errorMessage);
          }
        }
        if(do_cols_fe) {
          //cols_views
          NumeredDofEntity_multiIndex_uid_view_ordered cols_view;
          ierr = (*miit2)->getColDofView(
            *(p_miit->numered_dofs_cols),cols_view,moab::Interface::UNION
          ); CHKERRQ(ierr);
          //cols element dof multi-indices
          NumeredDofEntity_multiIndex_uid_view_ordered::iterator viit_cols;;
          viit_cols = cols_view.begin();
          for(;viit_cols!=cols_view.end();viit_cols++) {
            try {
              boost::shared_ptr<SideNumber> side_number_ptr = (*miit2)->getSideNumberPtr(moab,(*viit_cols)->getEnt());
              cols_dofs->insert(boost::shared_ptr<FENumeredDofEntity>(new FENumeredDofEntity(side_number_ptr,*viit_cols)));
            } catch (MoFEMException const &e) {
              SETERRQ(comm,e.errorCode,e.errorMessage);
            }
          }
        }
      }
      std::pair<NumeredEntFiniteElement_multiIndex::iterator,bool> p;
      p = numeredFiniteElements.insert(numered_fe);
      if(!p.second) {
        SETERRQ(comm,MOFEM_NOT_FOUND,"element is there");
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
  if(!(*buildMoFEM&BUILD_FIELD)) SETERRQ(comm,MOFEM_DATA_INCONSISTENCY,"fields not build");
  if(!(*buildMoFEM&BUILD_FE)) SETERRQ(comm,MOFEM_DATA_INCONSISTENCY,"FEs not build");
  if(!(*buildMoFEM&BUILD_ADJ)) SETERRQ(comm,MOFEM_DATA_INCONSISTENCY,"adjacencies not build");
  if(!(*buildMoFEM&BUILD_PROBLEM)) SETERRQ(comm,MOFEM_DATA_INCONSISTENCY,"problem not build");
  if(!(*buildMoFEM&PARTITION_PROBLEM)) SETERRQ(comm,MOFEM_DATA_INCONSISTENCY,"partition of problem not build");
  if(!(*buildMoFEM&PARTITION_FE)) SETERRQ(comm,MOFEM_DATA_INCONSISTENCY,"partitions finite elements not build");
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
          if((*rowdofit)->getPart()==(unsigned int)rAnk) continue;
          ghost_idx_row_view.insert((*rowdofit)->getNumeredDofEntityPtr());
        }
      }
      if((*fe_it)->cols_dofs->size()>0) {
        dof_it coldofit,hi_coldofit;
        coldofit = (*fe_it)->cols_dofs->begin();
        hi_coldofit = (*fe_it)->cols_dofs->end();
        for(;coldofit!=hi_coldofit;coldofit++) {
          if((*coldofit)->getPart()==(unsigned int)rAnk) continue;
          ghost_idx_col_view.insert((*coldofit)->getNumeredDofEntityPtr());
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
          SETERRQ(comm,MOFEM_DATA_INCONSISTENCY,"inconsistent data, ghost dof already set");
        }
        bool success = dof_by_uid_no_const[ss]->modify(diit,NumeredDofEntity_local_idx_change(nb_local_dofs[ss]++));
        if(!success) SETERRQ(comm,MOFEM_OPERATION_UNSUCCESSFUL,"modification unsuccessful");
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
    << " Nb. col ghost dof " << p_miit->getNbGhostDofsCol()
    << " Nb. local dof " << p_miit->getNbLocalDofsCol() << std::endl;
    ss << "partition_ghost_row_dofs: rank = " << rAnk
    << " FEs row ghost dofs "<< *p_miit
    << " Nb. row ghost dof " << p_miit->getNbGhostDofsRow()
    << " Nb. local dof " << p_miit->getNbLocalDofsRow() << std::endl;
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
      if(!success) {
        SETERRQ(
          comm,MOFEM_OPERATION_UNSUCCESSFUL,"modification unsuccessful"
        );
      }
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
        SETERRQ(comm,MOFEM_NOT_IMPLEMENTED,"not implemented");
      }
      if(verb>3) {
        std::ostringstream ss;
        ss << *(p_MoFEMFiniteElement.first->getRefElement());
        PetscSynchronizedPrintf(comm,"%s\n",ss.str().c_str());
      }
    }
    if(!seeded_ents.empty()) {
      int dim = moab.dimension_from_handle(seeded_ents[0]);
      for(int dd = 0;dd<dim;dd++) {
        Range ents;
        rval = moab.get_adjacencies(seeded_ents,dd,true,ents,moab::Interface::UNION); CHKERRQ_MOAB(rval);
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
          if(!success) SETERRQ(comm,MOFEM_OPERATION_UNSUCCESSFUL,"modification unsuccessful");
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
    SETERRQ(comm,e.errorCode,e.errorMessage);
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
    ss << "add meshset as ref_ent " << *(p_MoFEMFiniteElement.first->getRefElement()) << std::endl;
    PetscPrintf(comm,ss.str().c_str());
  }
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

  typedef RefEntity_multiIndex::index<Composite_Ent_And_ParentEntType_mi_tag>::type RefEntsByComposite;
  RefEntsByComposite &ref_ents = refinedEntities.get<Composite_Ent_And_ParentEntType_mi_tag>();
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
  typedef MoFEMProblem_multiIndex::index<Problem_mi_tag>::type mofem_problems_by_name;
  mofem_problems_by_name &pRoblems_set = pRoblems.get<Problem_mi_tag>();
  mofem_problems_by_name::iterator p_miit = pRoblems_set.find(problem_name);
  if(p_miit == pRoblems_set.end()) SETERRQ1(comm,1,"no such problem like < %s >",problem_name.c_str());
  NumeredEntFiniteElement_multiIndex &numeredFiniteElements = const_cast<NumeredEntFiniteElement_multiIndex&>(p_miit->numeredFiniteElements);
  NumeredEntFiniteElement_multiIndex::index<FiniteElement_name_mi_tag>::type::iterator miit = numeredFiniteElements.get<FiniteElement_name_mi_tag>().lower_bound(fe_name);
  for(;miit!=numeredFiniteElements.get<FiniteElement_name_mi_tag>().upper_bound(fe_name);miit++) {
    EntityHandle ent = (*miit)->getEnt();
    rval = moab.add_entities(meshset,&ent,1); CHKERRQ_MOAB(rval);
    int part = (*miit)->getPart();
    rval = moab.tag_set_data(th_Part,&ent,1,&part); CHKERRQ_MOAB(rval);
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
  if(p_miit == pRoblems_set.end()) SETERRQ1(comm,1,"problem is not in database %s",problem_name.c_str());
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
  if(p_miit == pRoblems_set.end()) SETERRQ1(comm,1,"problem is not in database %s",problem_name.c_str());

  ierr = problem_basic_method_postProcess(&*p_miit,method,verb); CHKERRQ(ierr);

  PetscFunctionReturn(0);
}
PetscErrorCode Core::loop_finite_elements(const std::string &problem_name,const std::string &fe_name,FEMethod &method,MoFEMTypes bh,int verb) {
  PetscFunctionBegin;
  if(verb==-1) verb = verbose;

  ierr = loop_finite_elements(problem_name,fe_name,method,rAnk,rAnk,bh,verb); CHKERRQ(ierr);

  PetscFunctionReturn(0);
}
PetscErrorCode Core::loop_finite_elements(
  const MoFEMProblem *problem_ptr,const std::string &fe_name,FEMethod &method,int lower_rank,int upper_rank,MoFEMTypes bh,int verb) {
  PetscFunctionBegin;
  if(verb==-1) verb = verbose;
  // finite element

  method.feName = fe_name;
  SET_BASIC_METHOD(method,&*problem_ptr)
  PetscLogEventBegin(USER_EVENT_preProcess,0,0,0,0);
  ierr = method.preProcess(); CHKERRQ(ierr);
  PetscLogEventEnd(USER_EVENT_preProcess,0,0,0,0);

  typedef NumeredEntFiniteElement_multiIndex::index<Composite_Name_And_Part_mi_tag>::type FEByComposite;

  FEByComposite &numered_fe =
    (const_cast<NumeredEntFiniteElement_multiIndex&>(problem_ptr->numeredFiniteElements)).get<Composite_Name_And_Part_mi_tag>();
  FEByComposite::iterator miit = numered_fe.lower_bound(boost::make_tuple(fe_name,lower_rank));
  FEByComposite::iterator hi_miit = numered_fe.upper_bound(boost::make_tuple(fe_name,upper_rank));

  if(miit==hi_miit && bh&MF_EXIST) {
    if(!check_finite_element(fe_name)) {
      SETERRQ1(comm,MOFEM_NOT_FOUND,"finite element < %s > not found",fe_name.c_str());
    }
  }

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
      SETERRQ(comm,MOFEM_DATA_INCONSISTENCY,ss.str().c_str());
    }

  }

  PetscLogEventBegin(USER_EVENT_postProcess,0,0,0,0);
  ierr = method.postProcess(); CHKERRQ(ierr);
  PetscLogEventEnd(USER_EVENT_postProcess,0,0,0,0);

  PetscFunctionReturn(0);
}
PetscErrorCode Core::loop_finite_elements(
  const std::string &problem_name,const std::string &fe_name,FEMethod &method,int lower_rank,int upper_rank,MoFEMTypes bh,int verb) {
  PetscFunctionBegin;
  if(verb==-1) verb = verbose;
  typedef MoFEMProblem_multiIndex::index<Problem_mi_tag>::type mofem_problems_by_name;
  // find p_miit
  mofem_problems_by_name &pRoblems_set = pRoblems.get<Problem_mi_tag>();
  mofem_problems_by_name::iterator p_miit = pRoblems_set.find(problem_name);
  if(p_miit == pRoblems_set.end()) SETERRQ1(comm,1,"problem is not in database %s",problem_name.c_str());

  ierr = loop_finite_elements(&*p_miit,fe_name,method,lower_rank,upper_rank,bh,verb); CHKERRQ(ierr);

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
     SETERRQ(comm,MOFEM_DATA_INCONSISTENCY,"not implemented");
  }
  numerd_dofs::iterator miit = dofs->lower_bound(boost::make_tuple(field_name,lower_rank));
  numerd_dofs::iterator hi_miit = dofs->upper_bound(boost::make_tuple(field_name,upper_rank));
  ierr = method.preProcess(); CHKERRQ(ierr);
  for(;miit!=hi_miit;miit++) {
    method.dofPtr = &(*(*miit)->getDofEntityPtr());
    method.dofNumeredPtr = &*(*miit);
    try {
      ierr = method(); CHKERRQ(ierr);
    } catch (const std::exception& ex) {
      std::ostringstream ss;
      ss << "throw in method: " << ex.what() << std::endl;
      SETERRQ(comm,MOFEM_DATA_INCONSISTENCY,ss.str().c_str());
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
  if(p_miit == pRoblems_set.end()) SETERRQ1(comm,1,"problem not in database %s",problem_name.c_str());
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
    } catch (MoFEMException const &e) {
      SETERRQ(comm,e.errorCode,e.errorMessage);
    } catch (const std::exception& ex) {
      std::ostringstream ss;
      ss << "throw in method: " << ex.what() << std::endl;
      SETERRQ(comm,MOFEM_DATA_INCONSISTENCY,ss.str().c_str());
    }
  }
  ierr = method.postProcess(); CHKERRQ(ierr);
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
