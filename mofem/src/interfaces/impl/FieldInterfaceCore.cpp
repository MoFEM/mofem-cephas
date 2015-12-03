/** \file FieldInterfaceCore.cpp
 * \brief Mylti-index containers, data structures and other low-level functions
 */

/* The MoFEM package is copyrighted by Lukasz Kaczmarczyk.
 * It can be freely used for educational and research purposes
 * by other institutions. If you use this software pleas cite my work.
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

#include <Includes.hpp>
#include <definitions.h>
#include <Common.hpp>

#include <h1_hdiv_hcurl_l2.h>

#include <MaterialBlocks.hpp>
#include <CubitBCData.hpp>
#include <TagMultiIndices.hpp>
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

PetscErrorCode Core::add_field(const string& name,const FieldSpace space,const ApproximationRank rank,enum MoFEMTypes bh,int verb) {
  PetscFunctionBegin;
  if(verb==-1) verb = verbose;
  *build_MoFEM = 0;
  MoFEMField_multiIndex::index<FieldName_mi_tag>::type::iterator fit;
  fit = moabFields.get<FieldName_mi_tag>().find(name);
  if(fit != moabFields.get<FieldName_mi_tag>().end() ) {
    if(bh == MF_EXCL) {
      SETERRQ1(PETSC_COMM_SELF,1,"field is <%s> in database",name.c_str());
    }
  } else {
    EntityHandle meshset;
    rval = moab.create_meshset(MESHSET_SET|MESHSET_TRACK_OWNER,meshset); CHKERR_PETSC(rval);
    //id
    BitFieldId id = getFieldShift();
    rval = moab.tag_set_data(th_FieldId,&meshset,1,&id); CHKERR_PETSC(rval);
    //space
    rval = moab.tag_set_data(th_FieldSpace,&meshset,1,&space); CHKERR_PETSC(rval);
    //add meshset to ref_ents // meshset dof on all level sets
    if(space == NOFIELD) {
      pair<RefMoFEMEntity_multiIndex::iterator,bool> p_ref_ent = refinedEntities.insert(RefMoFEMEntity(moab,meshset));
      bool success = refinedEntities.modify(p_ref_ent.first,RefMoFEMEntity_change_add_bit(BitRefLevel().set()));
      if(!success) SETERRQ(PETSC_COMM_SELF,MOFEM_OPERATION_UNSUCCESSFUL,"modification unsucceeded");
    }
    //name
    void const* tag_data[] = { name.c_str() };
    int tag_sizes[1]; tag_sizes[0] = name.size();
    rval = moab.tag_set_by_ptr(th_FieldName,&meshset,1,tag_data,tag_sizes); CHKERR_PETSC(rval);
    //name data prefix
    string name_data_prefix("_App_Data");
    void const* tag_prefix_data[] = { name_data_prefix.c_str() };
    int tag_prefix_sizes[1]; tag_prefix_sizes[0] = name_data_prefix.size();
    rval = moab.tag_set_by_ptr(th_FieldName_DataNamePrefix,&meshset,1,tag_prefix_data,tag_prefix_sizes); CHKERR_PETSC(rval);
    Tag th_AppOrder,th_FieldData,th_Rank,th_AppDofOrder,th_DofRank;
    //data
    string Tag_data_name = name_data_prefix+name;
    const int def_len = 0;
    rval = moab.tag_get_handle(Tag_data_name.c_str(),def_len,MB_TYPE_OPAQUE,
      th_FieldData,MB_TAG_CREAT|MB_TAG_SPARSE|MB_TAG_BYTES|MB_TAG_VARLEN,NULL); CHKERR_PETSC(rval);
    //order
    ApproximationOrder def_ApproximationOrder = -1;
    string Tag_ApproximationOrder_name = "_App_Order_"+name;
    rval = moab.tag_get_handle(Tag_ApproximationOrder_name.c_str(),sizeof(ApproximationOrder),MB_TYPE_OPAQUE,
      th_AppOrder,MB_TAG_CREAT|MB_TAG_SPARSE|MB_TAG_BYTES,&def_ApproximationOrder); CHKERR_PETSC(rval);
    //dof order
    string Tag_dof_ApproximationOrder_name = "_App_Dof_Order"+name;
    rval = moab.tag_get_handle(Tag_dof_ApproximationOrder_name.c_str(),def_len,MB_TYPE_OPAQUE,
      th_AppDofOrder,MB_TAG_CREAT|MB_TAG_SPARSE|MB_TAG_BYTES|MB_TAG_VARLEN,NULL); CHKERR_PETSC(rval);
    //rank
    int def_rank = 1;
    string Tag_rank_name = "_Field_Rank_"+name;
    rval = moab.tag_get_handle(Tag_rank_name.c_str(),sizeof(ApproximationRank),MB_TYPE_OPAQUE,
      th_Rank,MB_TAG_CREAT|MB_TAG_SPARSE|MB_TAG_BYTES,&def_rank); CHKERR_PETSC(rval);
    rval = moab.tag_set_data(th_Rank,&meshset,1,&rank); CHKERR_PETSC(rval);
    //dof rank
    string Tag_dof_rank_name = "_Field_Dof_Rank_"+name;
    rval = moab.tag_get_handle(Tag_dof_rank_name.c_str(),def_len,MB_TYPE_OPAQUE,
      th_DofRank,MB_TAG_CREAT|MB_TAG_SPARSE|MB_TAG_BYTES|MB_TAG_VARLEN,NULL); CHKERR_PETSC(rval);
    //add meshset
    pair<MoFEMField_multiIndex::iterator,bool> p;
    try {
      p = moabFields.insert(MoFEMField(moab,meshset));
      if(bh == MF_EXCL) {
        if(!p.second) SETERRQ1(PETSC_COMM_SELF,1,
  	"field not inserted %s (top tip, it could be already there)",
  	MoFEMField(moab,meshset).get_name().c_str());
      }
    } catch (MoFEMException const &e) {
      SETERRQ(PETSC_COMM_SELF,e.errorCode,e.errorMessage);
    }
    if(verbose > 0) {
      ostringstream ss;
      ss << "add: " << *p.first << endl;
      PetscPrintf(comm,ss.str().c_str());
    }
  }
  //
  PetscFunctionReturn(0);
}
PetscErrorCode Core::rebuild_database(int verb) {
  PetscFunctionBegin;
  if(verb==-1) verb = verbose;
  ierr = clearMap(); CHKERRQ(ierr);
  ierr = initialiseDatabseInformationFromMesh(verb); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
PetscErrorCode Core::initialiseDatabseInformationFromMesh(int verb) {
  PetscFunctionBegin;
  if(verb==-1) verb = verbose;
  Range meshsets;
  rval = moab.get_entities_by_type(0,MBENTITYSET,meshsets,true);  CHKERR_PETSC(rval);
  //loop all meshsets in moab database
  Range::iterator mit = meshsets.begin();
  for(;mit!=meshsets.end();mit++) {
    try {
      //check if meshset is cubit meshset
      CubitMeshSets base_meshset(moab,*mit);
      if((base_meshset.cubit_bc_type&CubitBCType(NODESET|SIDESET|BLOCKSET)).any()) {
        pair<CubitMeshSet_multiIndex::iterator,bool> p = cubitMeshsets.insert(base_meshset);
        if(!p.second) {
          SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"meshset not inserted");
        }
        if(verb > 0) {
          ostringstream ss;
          ss << "read cubit" << base_meshset << endl;
          //PetscSynchronizedPrintf(comm,ss.str().c_str());
          PetscPrintf(comm,ss.str().c_str());
        }
      }
    } catch (MoFEMException const &e) {
      SETERRQ(PETSC_COMM_SELF,e.errorCode,e.errorMessage);
    }
  }
  //PetscSynchronizedFlush(comm,PETSC_STDOUT);
  mit = meshsets.begin();
  for(;mit!=meshsets.end();mit++) {
    BitFieldId field_id;
    //get bit id form field tag
    rval = moab.tag_get_data(th_FieldId,&*mit,1,&field_id); CHKERR_PETSC(rval);
    //check if meshset if field meshset
    if(field_id!=0) {
      pair<MoFEMField_multiIndex::iterator,bool> p;
      try {
        p = moabFields.insert(MoFEMField(moab,*mit));
        if(verb > 0) {
          ostringstream ss;
          ss << "read field " << *p.first << endl;;
          PetscPrintf(comm,ss.str().c_str());
        }
      } catch (MoFEMException const &e) {
        SETERRQ(PETSC_COMM_SELF,e.errorCode,e.errorMessage);
      }
      if(p.first->get_space()==NOFIELD) {
        assert(p.first->meshset == *mit);
        //add field to ref ents
        pair<RefMoFEMEntity_multiIndex::iterator,bool> p_ref_ent = refinedEntities.insert(RefMoFEMEntity(moab,*mit));
        NOT_USED(p_ref_ent);
      } else {
        Range ents;
        rval = moab.get_entities_by_handle(*mit,ents,false); CHKERR_PETSC(rval);
        if(verb > 1) {
          ostringstream ss;
          ss << "read field ents " << ents.size() << endl;;
          PetscPrintf(comm,ss.str().c_str());
        }
        Range::iterator eit = ents.begin();
        for(;eit!=ents.end();eit++) {
          pair<RefMoFEMEntity_multiIndex::iterator,bool> p_ref_ent = refinedEntities.insert(RefMoFEMEntity(moab,*eit));
          try {
            MoFEMEntity moabent(moab,&*p.first,&*p_ref_ent.first);
            pair<MoFEMEntity_multiIndex::iterator,bool> p_ent = entsMoabField.insert(moabent);
            NOT_USED(p_ent);
          } catch (const std::exception& ex) {
            ostringstream ss;
            ss << ex.what() << endl;
            SETERRQ(PETSC_COMM_SELF,MOFEM_STD_EXCEPTION_THROW,ss.str().c_str());
          }
        }
      }
    }
    BitFieldId fe_id;
    //get bit id from fe tag
    rval = moab.tag_get_data(th_FEId,&*mit,1,&fe_id); CHKERR_PETSC(rval);
    //check if meshset is finite element meshset
    if(fe_id!=0) {
      pair<MoFEMFiniteElement_multiIndex::iterator,bool> p = finiteElements.insert(MoFEMFiniteElement(moab,*mit));
      if(verb > 0) {
        ostringstream ss;
        ss << "read finite element " << *p.first << endl;;
        PetscPrintf(comm,ss.str().c_str());
      }
      NOT_USED(p);
      assert(p.first->meshset == *mit);
      Range ents;
      rval = moab.get_entities_by_type(*mit,MBENTITYSET,ents,false); CHKERR_PETSC(rval);
      rval = moab.get_entities_by_handle(*mit,ents,true); CHKERR_PETSC(rval);
      Range::iterator eit = ents.begin();
      for(;eit!=ents.end();eit++) {
        pair<RefMoFEMEntity_multiIndex::iterator,bool> p_ref_ent = refinedEntities.insert(RefMoFEMEntity(moab,*eit));
        pair<RefMoFEMElement_multiIndex::iterator,bool> p_MoFEMFiniteElement;
        try {
          switch (moab.type_from_handle(*eit)) {
            case MBVERTEX:
            p_MoFEMFiniteElement = refinedFiniteElements.insert(ptrWrapperRefMoFEMElement(new RefMoFEMElement_VERTEX(moab,&*p_ref_ent.first)));
            break;
            case MBEDGE:
            p_MoFEMFiniteElement = refinedFiniteElements.insert(ptrWrapperRefMoFEMElement(new RefMoFEMElement_EDGE(moab,&*p_ref_ent.first)));
            break;
            case MBTRI:
            p_MoFEMFiniteElement = refinedFiniteElements.insert(ptrWrapperRefMoFEMElement(new RefMoFEMElement_TRI(moab,&*p_ref_ent.first)));
            break;
            case MBTET:
            p_MoFEMFiniteElement = refinedFiniteElements.insert(ptrWrapperRefMoFEMElement(new RefMoFEMElement_TET(moab,&*p_ref_ent.first)));
            break;
            case MBPRISM:
            ierr = addPrismToDatabase(*eit,verb); CHKERRQ(ierr);
            p_MoFEMFiniteElement = refinedFiniteElements.insert(ptrWrapperRefMoFEMElement(new RefMoFEMElement_PRISM(moab,&*p_ref_ent.first)));
            break;
            case MBENTITYSET:
            p_MoFEMFiniteElement = refinedFiniteElements.insert(ptrWrapperRefMoFEMElement(new RefMoFEMElement_MESHSET(moab,&*p_ref_ent.first)));
            break;
            default:
            SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"Only finite elements of type MBTET, MBPRISM and MBENTITYSET are implemented");
          }
          if(p_MoFEMFiniteElement.second) {
            //PetscPrintf(comm,"Warrning: this entity should be already in refined finite elements database");
            //SETERRQ(PETSC_COMM_SELF,1,"data inconsistency, this entity should be already in refined finite elements database");
          }
        } catch (MoFEMException const &e) {
          SETERRQ(PETSC_COMM_SELF,e.errorCode,e.errorMessage);
        }
      }
    }
    BitProblemId problem_id;
    //get bit id form problem tag
    rval = moab.tag_get_data(th_ProblemId,&*mit,1,&problem_id); CHKERR_PETSC(rval);
    //check if meshset if problem meshset
    if(problem_id!=0) {
      pair<MoFEMProblem_multiIndex::iterator,bool> p = moFEMProblems.insert(MoFEMProblem(moab,*mit));
      if(verb > 0) {
        ostringstream ss;
        ss << "read problem " << *p.first << endl;;
        PetscPrintf(comm,ss.str().c_str());
      }
    }
    //check if meshset is Series meshset
    {
      const void* tag_name_data;
      int tag_name_size;
      rval = moab.tag_get_by_ptr(th_SeriesName,&*mit,1,(const void **)&tag_name_data,&tag_name_size);
      if(rval == MB_SUCCESS) {
        pair<Series_multiIndex::iterator,bool> p = sEries.insert(MoFEMSeries(moab,*mit));
        if(verb > 0) {
          ostringstream ss;
          ss << "read series " << *p.first << endl;
          PetscPrintf(comm,ss.str().c_str());
        }
      }
    }
  }
  //build ref entities meshset
  for(int dd = 0;dd<=3;dd++) {
    Range ents;
    rval = moab.get_entities_by_dimension(0,dd,ents,false); CHKERR_PETSC(rval);
    Range::iterator eit = ents.begin();
    for(;eit!=ents.end();eit++) {
      switch (moab.type_from_handle(*eit)) {
        case MBVERTEX:
        case MBEDGE:
        case MBTRI:
        case MBTET:
        case MBPRISM:
        break;
        default:
        continue;
      }
      RefMoFEMEntity mofem_ent(moab,*eit);
      BitRefLevel bit = mofem_ent.get_BitRefLevel();
      if(bit.none()) {
        continue;
      }
      pair<RefMoFEMEntity_multiIndex::iterator,bool> p;
      p = refinedEntities.insert(mofem_ent);
    }
  }
  //build series steps
  for(Series_multiIndex::iterator sit = sEries.begin();sit!=sEries.end();sit++) {
    int nb_steps;
    ierr = sit->get_nb_steps(moab,nb_steps); CHKERRQ(ierr);
    int ss = 0;
    for(;ss<nb_steps;ss++) {
      pair<SeriesStep_multiIndex::iterator,bool> p = seriesSteps.insert(MoFEMSeriesStep(moab,&*sit,ss));
      if(verb > 0) {
        ostringstream ss;
        ss << "add series step " << *p.first << endl;
        PetscPrintf(comm,ss.str().c_str());
      }
    }
  }
  if(verb > 2) {
    list_fields();
    list_finite_elements();
    list_problem();
  }
  PetscFunctionReturn(0);
}
PetscErrorCode Core::add_ents_to_field_by_EDGEs(const EntityHandle meshset,const BitFieldId id,int verb) {
  PetscFunctionBegin;
  if(verb==-1) verb = verbose;
  *build_MoFEM = 0;
  EntityHandle idm = no_handle;
  try {
    idm = get_field_meshset(id);
  } catch (MoFEMException const &e) {
    SETERRQ(PETSC_COMM_SELF,e.errorCode,e.errorMessage);
  }
  FieldSpace space;
  rval = moab.tag_get_data(th_FieldSpace,&idm,1,&space); CHKERR_PETSC(rval);
  Range nodes,edges;
  rval = moab.get_entities_by_type(meshset,MBEDGE,edges,true); CHKERR_PETSC(rval);
  switch (space) {
    case L2:
      rval = moab.add_entities(idm,edges); CHKERR_PETSC(rval);
      if(verb>1) {
        ostringstream ss;
        ss << "add entities to field " << get_BitFieldId_name(id);
        ss << " nb. add edges " << edges.size();
        ss << endl;
        PetscPrintf(comm,ss.str().c_str());
      }
      break;
    case H1:
      rval = moab.add_entities(idm,edges); CHKERR_PETSC(rval);
      //rval = moab.get_connectivity(edges,nodes,true); CHKERR_PETSC(rval);
      //use get adjacencies, this will allow take in account adjacencies set user
      rval = moab.get_adjacencies(edges,0,false,nodes,Interface::UNION); CHKERR_PETSC(rval);
      {
        Range topo_nodes;
        rval = moab.get_connectivity(edges,topo_nodes,true); CHKERR_PETSC(rval);
        Range mid_nodes;
        rval = moab.get_connectivity(edges,mid_nodes,false); CHKERR_PETSC(rval);
        mid_nodes = subtract(mid_nodes,topo_nodes);
        nodes = subtract(nodes,mid_nodes);
      }
      rval = moab.add_entities(idm,nodes); CHKERR_PETSC(rval);
      if(verb>1) {
        ostringstream ss;
        ss << "add entities to field " << get_BitFieldId_name(id);
        ss << " nb. add edges " << edges.size();
        ss << " nb. add nodes " << nodes.size();
        ss << endl;
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
PetscErrorCode Core::add_ents_to_field_by_EDGEs(const EntityHandle meshset,const string& name,int verb) {
  PetscFunctionBegin;
  if(verb==-1) verb = verbose;
  *build_MoFEM = 0;
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
  rval = moab.get_entities_by_type(meshset,MBTRI,tris,true); CHKERR_PETSC(rval);
  ierr = add_ents_to_field_by_TRIs(tris,id,verb); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
PetscErrorCode Core::add_ents_to_field_by_TRIs(const Range &tris,const BitFieldId id,int verb) {
  PetscFunctionBegin;
  if(verb==-1) verb = verbose;
  *build_MoFEM = 0;
  EntityHandle idm = no_handle;
  try {
    idm = get_field_meshset(id);
  } catch (MoFEMException const &e) {
    SETERRQ(PETSC_COMM_SELF,e.errorCode,e.errorMessage);
  }
  FieldSpace space;
  rval = moab.tag_get_data(th_FieldSpace,&idm,1,&space); CHKERR_PETSC(rval);
  Range nodes,edges;
  switch(space) {
    case L2:
    rval = moab.add_entities(idm,tris); CHKERR_PETSC(rval);
    if(verb>1) {
      ostringstream ss;
      ss << "add entities to field " << get_BitFieldId_name(id);
      ss << " nb. add tris " << tris.size();
      ss << endl;
      PetscPrintf(comm,ss.str().c_str());
    }
    break;
    case H1:
    rval = moab.add_entities(idm,tris); CHKERR_PETSC(rval);
    //rval = moab.get_connectivity(tris,nodes,true); CHKERR_PETSC(rval);
    //use get adjacencies, this will allow take in account adjacencies set user
    rval = moab.get_adjacencies(tris,0,false,nodes,Interface::UNION); CHKERR_PETSC(rval);
    {
      Range topo_nodes;
      rval = moab.get_connectivity(tris,topo_nodes,true); CHKERR_PETSC(rval);
      Range mid_nodes;
      rval = moab.get_connectivity(tris,mid_nodes,false); CHKERR_PETSC(rval);
      mid_nodes = subtract(mid_nodes,topo_nodes);
      nodes = subtract(nodes,mid_nodes);
    }
    rval = moab.add_entities(idm,nodes); CHKERR_PETSC(rval);
    rval = moab.get_adjacencies(tris,1,false,edges,Interface::UNION); CHKERR_PETSC(rval);
    rval = moab.add_entities(idm,edges); CHKERR_PETSC(rval);
    if(verb>1) {
      ostringstream ss;
      ss << "add entities to field " << get_BitFieldId_name(id);
      ss << " nb. add tris " << tris.size();
      ss << " nb. add edges " << edges.size();
      ss << " nb. add nodes " << nodes.size();
      ss << endl;
      PetscPrintf(comm,ss.str().c_str());
    }
    break;
    case HCURL:
    SETERRQ(PETSC_COMM_SELF,MOFEM_NOT_IMPLEMENTED,"sorry, not implemented, HCURL not implented for TRI");
    break;
    case HDIV:
    SETERRQ(PETSC_COMM_SELF,MOFEM_NOT_IMPLEMENTED,"sorry, not implemented, HDIV not implemented for TRI");
    break;
    default:
    SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"add_ents_to_field_by_TRIs this field not work for TRIs");
  }
  PetscFunctionReturn(0);
}
PetscErrorCode Core::add_ents_to_field_by_TRIs(const EntityHandle meshset,const string& name,int verb) {
  PetscFunctionBegin;
  if(verb==-1) verb = verbose;
  *build_MoFEM = 0;
  try {
    ierr = add_ents_to_field_by_TRIs(meshset,get_BitFieldId(name),verb);  CHKERRQ(ierr);
  } catch (MoFEMException const &e) {
    SETERRQ(PETSC_COMM_SELF,e.errorCode,e.errorMessage);
  }
  PetscFunctionReturn(0);
}
PetscErrorCode Core::add_ents_to_field_by_TRIs(const Range &tris,const string& name,int verb) {
  PetscFunctionBegin;
  if(verb==-1) verb = verbose;
  *build_MoFEM = 0;
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
  *build_MoFEM = 0;
  EntityHandle idm = no_handle;
  try {
    idm = get_field_meshset(id);
  } catch (MoFEMException const &e) {
    SETERRQ(PETSC_COMM_SELF,e.errorCode,e.errorMessage);
  }
  FieldSpace space;
  rval = moab.tag_get_data(th_FieldSpace,&idm,1,&space); CHKERR_PETSC(rval);
  switch (space) {
    case H1:
    rval = moab.add_entities(idm,nodes); CHKERR_PETSC(rval);
    if(verb>1) {
      ostringstream ss;
      ss << "add entities to field " << get_BitFieldId_name(id);
      ss << " nb. add nodes " << nodes.size();
      ss << endl;
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
  rval = moab.get_entities_by_type(meshset,MBVERTEX,nodes,true); CHKERR_PETSC(rval);
  ierr = add_ents_to_field_by_VERTICEs(nodes,id,verb); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
PetscErrorCode Core::add_ents_to_field_by_VERTICEs(const Range &nodes,const string& name,int verb) {
  PetscFunctionBegin;
  if(verb==-1) verb = verbose;
  *build_MoFEM = 0;
  try {
    ierr = add_ents_to_field_by_VERTICEs(nodes,get_BitFieldId(name),verb);  CHKERRQ(ierr);
  } catch (MoFEMException const &e) {
    SETERRQ(PETSC_COMM_SELF,e.errorCode,e.errorMessage);
  }
  PetscFunctionReturn(0);
}
PetscErrorCode Core::add_ents_to_field_by_VERTICEs(const EntityHandle meshset,const string& name,int verb) {
  PetscFunctionBegin;
  if(verb==-1) verb = verbose;
  *build_MoFEM = 0;
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
  rval = moab.add_entities(idm,ents); CHKERR_PETSC(rval);
  PetscFunctionReturn(0);
}
PetscErrorCode Core::synchronise_field_entities(const string& name,int verb) {
  PetscFunctionBegin;
  if(verb==-1) verb = verbose;
  *build_MoFEM = 0;
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
  *build_MoFEM = 0;
  EntityHandle idm = no_handle;
  try {
    idm = get_field_meshset(id);
  } catch (MoFEMException const &e) {
    SETERRQ(PETSC_COMM_SELF,e.errorCode,e.errorMessage);
  }
  FieldSpace space;
  rval = moab.tag_get_data(th_FieldSpace,&idm,1,&space); CHKERR_PETSC(rval);
  Range nodes,tris,edges;
  switch(space) {
    case L2:
    rval = moab.add_entities(idm,tets); CHKERR_PETSC(rval);
    if(verb>1) {
      ostringstream ss;
      ss << "add entities to field " << get_BitFieldId_name(id);
      ss << " nb. add tets " << tets.size();
      ss << endl;
      PetscSynchronizedPrintf(comm,ss.str().c_str());
    }
    break;
    case H1:
    rval = moab.add_entities(idm,tets); CHKERR_PETSC(rval);
    //rval = moab.get_connectivity(tets,nodes,true); CHKERR_PETSC(rval);
    //use get adjacencies, this will allow take in account adjacencies set user
    rval = moab.get_adjacencies(tets,0,false,nodes,Interface::UNION); CHKERR_PETSC(rval);
    {
      Range topo_nodes;
      rval = moab.get_connectivity(tets,topo_nodes,true); CHKERR_PETSC(rval);
      Range mid_nodes;
      rval = moab.get_connectivity(tets,mid_nodes,false); CHKERR_PETSC(rval);
      mid_nodes = subtract(mid_nodes,topo_nodes);
      nodes = subtract(nodes,mid_nodes);
    }
    rval = moab.add_entities(idm,nodes); CHKERR_PETSC(rval);
    rval = moab.get_adjacencies(tets,2,false,tris,Interface::UNION); CHKERR_PETSC(rval);
    rval = moab.add_entities(idm,tris); CHKERR_PETSC(rval);
    rval = moab.get_adjacencies(tets,1,false,edges,Interface::UNION); CHKERR_PETSC(rval);
    rval = moab.add_entities(idm,edges); CHKERR_PETSC(rval);
    if(verb>1) {
      ostringstream ss;
      ss << "add entities to field " << get_BitFieldId_name(id);
      ss << " nb. add tets " << tets.size();
      ss << " nb. add tris " << tris.size();
      ss << " nb. add edges " << edges.size();
      ss << " nb. add nodes " << nodes.size();
      ss << endl;
      PetscSynchronizedPrintf(comm,ss.str().c_str());
    }
    break;
    case HCURL:
    rval = moab.add_entities(idm,tets); CHKERR_PETSC(rval);
    rval = moab.get_adjacencies(tets,2,false,tris,Interface::UNION); CHKERR_PETSC(rval);
    rval = moab.add_entities(idm,tris); CHKERR_PETSC(rval);
    rval = moab.get_adjacencies(tets,1,false,edges,Interface::UNION); CHKERR_PETSC(rval);
    rval = moab.add_entities(idm,edges); CHKERR_PETSC(rval);
    if(verb>1) {
      ostringstream ss;
      ss << "add entities to field " << get_BitFieldId_name(id);
      ss << " nb. add tets " << tets.size();
      ss << " nb. add tris " << tris.size();
      ss << " nb. add edges " << edges.size();
      ss << endl;
      PetscSynchronizedPrintf(comm,ss.str().c_str());
    }
    break;
    case HDIV:
    rval = moab.add_entities(idm,tets); CHKERR_PETSC(rval);
    rval = moab.get_adjacencies(tets,2,false,tris,Interface::UNION); CHKERR_PETSC(rval);
    rval = moab.add_entities(idm,tris); CHKERR_PETSC(rval);
    if(verb>1) {
      ostringstream ss;
      ss << "add entities to field " << get_BitFieldId_name(id);
      ss << " nb. add tets " << tets.size();
      ss << " nb. add tris " << tris.size();
      ss << endl;
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
  rval = moab.get_entities_by_type(meshset,MBTET,tets,true); CHKERR_PETSC(rval);
  if(verb>3) {
    PetscSynchronizedPrintf(comm,"nb. of tets %d\n",tets.size());
  }
  ierr = add_ents_to_field_by_TETs(tets,id,verb); CHKERRQ(ierr);
  if(verb>3) {
    PetscSynchronizedFlush(comm,PETSC_STDOUT);
  }
  PetscFunctionReturn(0);
}
PetscErrorCode Core::add_ents_to_field_by_TETs(const Range &tets,const string& name,int verb) {
  PetscFunctionBegin;
  if(verb==-1) verb = verbose;
  *build_MoFEM = 0;
  try {
    ierr = add_ents_to_field_by_TETs(tets,get_BitFieldId(name),verb);  CHKERRQ(ierr);
  } catch (MoFEMException const &e) {
    SETERRQ(PETSC_COMM_SELF,e.errorCode,e.errorMessage);
  }
  PetscFunctionReturn(0);
}
PetscErrorCode Core::add_ents_to_field_by_TETs(const EntityHandle meshset,const string& name,int verb) {
  PetscFunctionBegin;
  if(verb==-1) verb = verbose;
  *build_MoFEM = 0;
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
  *build_MoFEM = 0;
  EntityHandle idm = no_handle;
  try {
    idm = get_field_meshset(id);
  } catch (MoFEMException const &e) {
    SETERRQ(PETSC_COMM_SELF,e.errorCode,e.errorMessage);
  }
  FieldSpace space;
  rval = moab.tag_get_data(th_FieldSpace,&idm,1,&space); CHKERR_PETSC(rval);
  Range nodes,faces,edges;
  switch(space) {
    case L2:
    rval = moab.add_entities(idm,quads); CHKERR_PETSC(rval);
    if(verb>1) {
      ostringstream ss;
      ss << "add entities to field " << get_BitFieldId_name(id);
      ss << " nb. add quads " << quads.size();
      ss << endl;
      PetscSynchronizedPrintf(comm,ss.str().c_str());
    }
    break;
    case H1:
    rval = moab.add_entities(idm,quads); CHKERR_PETSC(rval);
    //rval = moab.get_connectivity(quads,nodes,true); CHKERR_PETSC(rval);
    //use get adjacencies, this will allow take in account adjacencies set user
    rval = moab.get_adjacencies(quads,0,false,nodes,Interface::UNION); CHKERR_PETSC(rval);
    {
      Range topo_nodes;
      rval = moab.get_connectivity(quads,topo_nodes,true); CHKERR_PETSC(rval);
      Range mid_nodes;
      rval = moab.get_connectivity(quads,mid_nodes,false); CHKERR_PETSC(rval);
      mid_nodes = subtract(mid_nodes,topo_nodes);
      nodes = subtract(nodes,mid_nodes);
    }
    rval = moab.add_entities(idm,nodes); CHKERR_PETSC(rval);
    rval = moab.get_adjacencies(quads,1,true,edges,Interface::UNION); CHKERR_PETSC(rval);
    rval = moab.add_entities(idm,edges); CHKERR_PETSC(rval);
    if(verb>1) {
      ostringstream ss;
      ss << "add entities to field " << get_BitFieldId_name(id);
      ss << " nb. add quads " << quads.size();
      ss << " nb. add edges " << edges.size();
      ss << " nb. add nodes " << nodes.size();
      ss << endl;
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
PetscErrorCode Core::add_ents_to_field_by_QUADs(const Range &quads,const string& name,int verb) {
  PetscFunctionBegin;
  if(verb==-1) verb = verbose;
  *build_MoFEM = 0;
  try {
    ierr = add_ents_to_field_by_QUADs(quads,get_BitFieldId(name),verb);  CHKERRQ(ierr);
  } catch (MoFEMException const &e) {
    SETERRQ(PETSC_COMM_SELF,e.errorCode,e.errorMessage);
  }
  PetscFunctionReturn(0);
}
PetscErrorCode Core::add_ents_to_field_by_QUADs(EntityHandle meshset,const string& name,int verb) {
  PetscFunctionBegin;
  if(verb==-1) verb = verbose;
  *build_MoFEM = 0;
  try {
    Range quads;
    rval = moab.get_entities_by_type(meshset,MBQUAD,quads,true); CHKERR_PETSC(rval);
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
  *build_MoFEM = 0;
  EntityHandle idm = no_handle;
  try {
    idm = get_field_meshset(id);
  } catch (MoFEMException const &e) {
    SETERRQ(PETSC_COMM_SELF,e.errorCode,e.errorMessage);
  }
  FieldSpace space;
  rval = moab.tag_get_data(th_FieldSpace,&idm,1,&space); CHKERR_PETSC(rval);
  Range nodes,faces,edges;
  switch(space) {
    case L2:
    rval = moab.add_entities(idm,prisms); CHKERR_PETSC(rval);
    if(verb>1) {
      ostringstream ss;
      ss << "add entities to field " << get_BitFieldId_name(id);
      ss << " nb. add prisms " << prisms.size();
      ss << endl;
      PetscSynchronizedPrintf(comm,ss.str().c_str());
    }
    break;
    case H1:
    rval = moab.add_entities(idm,prisms); CHKERR_PETSC(rval);
    //rval = moab.get_connectivity(prisms,nodes,true); CHKERR_PETSC(rval);
    //use get adjacencies, this will allow take in account adjacencies set user
    rval = moab.get_adjacencies(prisms,0,false,nodes,Interface::UNION); CHKERR_PETSC(rval);
    {
      Range topo_nodes;
      rval = moab.get_connectivity(prisms,topo_nodes,true); CHKERR_PETSC(rval);
      Range mid_nodes;
      rval = moab.get_connectivity(prisms,mid_nodes,false); CHKERR_PETSC(rval);
      mid_nodes = subtract(mid_nodes,topo_nodes);
      nodes = subtract(nodes,mid_nodes);
    }
    rval = moab.add_entities(idm,nodes); CHKERR_PETSC(rval);
    rval = moab.get_adjacencies(prisms,2,true,faces,Interface::UNION); CHKERR_PETSC(rval);
    rval = moab.add_entities(idm,faces); CHKERR_PETSC(rval);
    rval = moab.get_adjacencies(prisms,1,true,edges,Interface::UNION); CHKERR_PETSC(rval);
    rval = moab.add_entities(idm,edges); CHKERR_PETSC(rval);
    if(verb>1) {
      ostringstream ss;
      ss << "add entities to field " << get_BitFieldId_name(id);
      ss << " nb. add prisms " << prisms.size();
      ss << " nb. add faces " << faces.size();
      ss << " nb. add edges " << edges.size();
      ss << " nb. add nodes " << nodes.size();
      ss << endl;
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
PetscErrorCode Core::add_ents_to_field_by_PRISMs(const Range &prisms,const string& name,int verb) {
  PetscFunctionBegin;
  if(verb==-1) verb = verbose;
  *build_MoFEM = 0;
  try {
    ierr = add_ents_to_field_by_PRISMs(prisms,get_BitFieldId(name),verb);  CHKERRQ(ierr);
  } catch (MoFEMException const &e) {
    SETERRQ(PETSC_COMM_SELF,e.errorCode,e.errorMessage);
  }
  PetscFunctionReturn(0);
}
PetscErrorCode Core::add_ents_to_field_by_PRISMs(EntityHandle meshset,const string& name,int verb) {
  PetscFunctionBegin;
  if(verb==-1) verb = verbose;
  *build_MoFEM = 0;
  try {
    Range prisms;
    rval = moab.get_entities_by_type(meshset,MBPRISM,prisms,true); CHKERR_PETSC(rval);
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
  *build_MoFEM = 0;

  //check field & meshset
  typedef MoFEMField_multiIndex::index<BitFieldId_mi_tag>::type field_set_by_id;
  const field_set_by_id &set_id = moabFields.get<BitFieldId_mi_tag>();
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
  rval = moab.get_entities_by_handle(idm,ents_of_id_meshset,false); CHKERR_PETSC(rval);
  Range ents_ = intersect(ents,ents_of_id_meshset);
  if(verb>1) {
    PetscSynchronizedPrintf(comm,"nb. of ents for order change in the field %d\n",ents_.size());
  }

  //ent view by field id (in set all MoabEnts has the same FieldId)
  typedef MoFEMEntity_multiIndex::index<BitFieldId_mi_tag>::type ent_set_by_id;
  ent_set_by_id& set = entsMoabField.get<BitFieldId_mi_tag>();
  ent_set_by_id::iterator eiit = set.lower_bound(id);
  MoFEMEntity_multiIndex_ent_view ents_id_view;
  if(eiit != set.end()) {
    ent_set_by_id::iterator hi_eiit = set.upper_bound(id);
    for(;eiit!=hi_eiit;eiit++) {
      ents_id_view.insert(&*eiit);
    }
  }
  if(verb>1) {
    PetscSynchronizedPrintf(comm,"nb. of ents in the multi index field %d\n",ents_id_view.size());
  }


  // get tags on entities
  vector<ApproximationOrder*> tag_data_order(ents_.size());
  rval = moab.tag_get_by_ptr(miit->th_AppOrder,ents_,(const void **)&tag_data_order[0]); CHKERR_PETSC(rval);

  //loop over ents
  int nb_ents_set_order_up = 0;
  int nb_ents_set_order_down = 0;
  int nb_ents_set_order_new = 0;
  Range::iterator eit = ents_.begin();
  for(unsigned int ee = 0;ee<ents_.size();ee++,eit++) {

    //sanity check
    switch(miit->get_space()) {
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
      const ApproximationOrder old_approximation_order = (*vit)->get_max_order();
      if(old_approximation_order==order) continue;
      MoFEMEntity_multiIndex::iterator miit = entsMoabField.get<Unique_mi_tag>().find((*vit)->get_global_unique_id());
      if(miit->get_max_order()<order) nb_ents_set_order_up++;
      if(miit->get_max_order()>order) nb_ents_set_order_down++;

      {

      	//set dofs inactive if order is reduced, and set new order to entity if
      	//order is increased (note that dofs are not build if order is
      	//increased)

        typedef DofMoFEMEntity_multiIndex::index<Composite_Name_And_Ent_mi_tag>::type dof_set_type;
        dof_set_type& set_set = dofsMoabField.get<Composite_Name_And_Ent_mi_tag>();
        dof_set_type::iterator dit = set_set.lower_bound(boost::make_tuple(miit->get_name_ref(),miit->get_ent()));
        dof_set_type::iterator hi_dit = set_set.upper_bound(boost::make_tuple(miit->get_name_ref(),miit->get_ent()));

        for(;dit!=hi_dit;dit++) {
          if(dit->get_dof_order()<=order) continue;
          bool success = dofsMoabField.modify(dofsMoabField.project<0>(dit),DofMoFEMEntity_active_change(false));
          if(!success) SETERRQ(PETSC_COMM_SELF,MOFEM_OPERATION_UNSUCCESSFUL,"modification unsuccessful");
        }

        bool success = entsMoabField.modify(entsMoabField.project<0>(miit),MoFEMEntity_change_order(moab,order));
        if(!success) SETERRQ(PETSC_COMM_SELF,MOFEM_OPERATION_UNSUCCESSFUL,"modification unsuccessful");

      }

    } else {

      //entity is not in database and order is changed or reset
      *tag_data_order[ee] = order;
      RefMoFEMEntity_multiIndex::index<Ent_mi_tag>::type::iterator miit_ref_ent = refinedEntities.get<Ent_mi_tag>().find(*eit);
      if(miit_ref_ent==refinedEntities.get<Ent_mi_tag>().end()) {

        RefMoFEMEntity ref_ent(moab,*eit);
        // FIXME: need some consistent policy in that case
        if(ref_ent.get_BitRefLevel().none()) continue; // not on any mesh and not in database
        cerr << ref_ent << endl;
        cerr << "bit level " << ref_ent.get_BitRefLevel() << endl;
        SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"database inconsistency");

      }

      try {

        // increase order
        MoFEMEntity moabent(moab,&*miit,&*miit_ref_ent);
        pair<MoFEMEntity_multiIndex::iterator,bool> e_miit = entsMoabField.insert(moabent);
        bool success = entsMoabField.modify(e_miit.first,MoFEMEntity_change_order(moab,order));
        if(!success) SETERRQ(PETSC_COMM_SELF,MOFEM_OPERATION_UNSUCCESSFUL,"modification unsuccessful");
        nb_ents_set_order_new++;

      } catch (MoFEMException const &e) {
        SETERRQ(PETSC_COMM_SELF,e.errorCode,e.errorMessage);
      } catch (const std::exception& ex) {
        ostringstream ss;
        ss << ex.what() << endl;
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
  *build_MoFEM = 0;
  Range ents;
  rval = moab.get_entities_by_type(meshset,type,ents); CHKERR_PETSC(rval);
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
PetscErrorCode Core::set_field_order(const EntityHandle meshset,const EntityType type,const string& name,const ApproximationOrder order,int verb) {
  PetscFunctionBegin;
  if(verb==-1) verb = verbose;
  *build_MoFEM = 0;
  try{
    ierr = set_field_order(meshset,type,get_BitFieldId(name),order,verb); CHKERRQ(ierr);
  } catch (MoFEMException const &e) {
    SETERRQ(PETSC_COMM_SELF,e.errorCode,e.errorMessage);
  }
  PetscFunctionReturn(0);
}
PetscErrorCode Core::set_field_order(const Range &ents,const string& name,const ApproximationOrder order,int verb) {
  PetscFunctionBegin;
  if(verb==-1) verb = verbose;
  *build_MoFEM = 0;
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
  *build_MoFEM = 0;
  Range ents;
  ierr = get_entities_by_type_and_ref_level(bit,mask,type,ents,verb); CHKERRQ(ierr);
  try{
    ierr = set_field_order(ents,id,order,verb); CHKERRQ(ierr);
  } catch (MoFEMException const &e) {
    SETERRQ(PETSC_COMM_SELF,e.errorCode,e.errorMessage);
  }
  PetscFunctionReturn(0);
}
PetscErrorCode Core::set_field_order_by_entity_type_and_bit_ref(const BitRefLevel &bit,const BitRefLevel &mask,const EntityType type,const string& name,const ApproximationOrder order,int verb) {
  PetscFunctionBegin;
  if(verb==-1) verb = verbose;
  *build_MoFEM = 0;
  Range ents;
  ierr = get_entities_by_type_and_ref_level(bit,mask,type,ents,verb); CHKERRQ(ierr);
  try{
    ierr = set_field_order(ents,get_BitFieldId(name),order,verb); CHKERRQ(ierr);
  } catch (MoFEMException const &e) {
    SETERRQ(PETSC_COMM_SELF,e.errorCode,e.errorMessage);
  }
  PetscFunctionReturn(0);
}
PetscErrorCode Core::dofs_NoField(const BitFieldId id,map<EntityType,int> &dof_counter,int verb) {
  PetscFunctionBegin;
  if(verb==-1) verb = verbose;
  //field it
  typedef MoFEMField_multiIndex::index<BitFieldId_mi_tag>::type field_set_by_id;
  const field_set_by_id &set_id = moabFields.get<BitFieldId_mi_tag>();
  //find fiels
  field_set_by_id::iterator miit = set_id.find(id);
  if(miit == set_id.end()) SETERRQ(PETSC_COMM_SELF,MOFEM_NOT_FOUND,"field not found");
  //serch if field meshset is in database
  RefMoFEMEntity_multiIndex::index<Ent_mi_tag>::type::iterator miit_ref_ent = refinedEntities.get<Ent_mi_tag>().find(miit->meshset);
  if(miit_ref_ent==refinedEntities.get<Ent_mi_tag>().end()) SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"database inconsistency");
  pair<MoFEMEntity_multiIndex::iterator,bool> e_miit;
  try {
    //create database entity
    MoFEMEntity moabent(moab,&*miit,&*miit_ref_ent);
    e_miit = entsMoabField.insert(moabent);
    //this is nor real field in space (set order to zero)
    bool success = entsMoabField.modify(e_miit.first,MoFEMEntity_change_order(moab,0));
    if(!success) SETERRQ(PETSC_COMM_SELF,MOFEM_OPERATION_UNSUCCESSFUL,"modification unsuccessful");
  } catch (MoFEMException const &e) {
    SETERRQ(PETSC_COMM_SELF,e.errorCode,e.errorMessage);
  } catch (const std::exception& ex) {
    ostringstream ss;
    ss << ex.what() << endl;
    SETERRQ(PETSC_COMM_SELF,MOFEM_STD_EXCEPTION_THROW,ss.str().c_str());
  }
  assert(e_miit.first->get_ent()==miit->meshset);
  ApproximationRank rank = 0;
  //create dofs on this entity (nb. of dofs is equal to rank)
  for(;rank<e_miit.first->get_max_rank();rank++) {
    pair<DofMoFEMEntity_multiIndex::iterator,bool> d_miit;
    //check if dof is in darabase
    d_miit.first = dofsMoabField.project<0>(
      dofsMoabField.get<Unique_mi_tag>().find(DofMoFEMEntity::get_global_unique_id_calculate(rank,&*(e_miit.first)))
    );
    //if dof is not in databse
    if(d_miit.first==dofsMoabField.end()) {
      //insert dof
      d_miit = dofsMoabField.insert(DofMoFEMEntity(&*(e_miit.first),0,rank,rank));
      if(d_miit.second) {
        dof_counter[d_miit.first->get_ent_type()]++;
      }
      bool success = dofsMoabField.modify(d_miit.first,DofMoFEMEntity_active_change(true));
      if(!success) SETERRQ(PETSC_COMM_SELF,MOFEM_OPERATION_UNSUCCESSFUL,"modification unsuccessful");
    }
    //check consistency
    assert(d_miit.first->get_ent()==e_miit.first->get_MoFEMField_ptr()->meshset);
    assert(d_miit.first->get_ent_type()==e_miit.first->get_ent_type());
    assert(d_miit.first->get_id()==e_miit.first->get_id());
  }
  if(verb>2) {
    typedef DofMoFEMEntity_multiIndex::index<BitFieldId_mi_tag>::type dof_set_by_id;
    dof_set_by_id &set = dofsMoabField.get<BitFieldId_mi_tag>();
    dof_set_by_id::iterator miit2 = set.lower_bound(id);
    dof_set_by_id::iterator hi_miit2 = set.upper_bound(id);
    assert(miit2!=hi_miit2);
    for(;miit2!=hi_miit2;miit2++) {
      ostringstream ss;
      ss << *miit2 << endl;;
      PetscSynchronizedPrintf(comm,ss.str().c_str());
    }
    PetscSynchronizedFlush(comm,PETSC_STDOUT);
  }
  PetscFunctionReturn(0);
}
PetscErrorCode Core::dofs_L2H1HcurlHdiv(const BitFieldId id,map<EntityType,int> &dof_counter,int verb) {
  PetscFunctionBegin;
  if(verb==-1) verb = verbose;
  //field it
  typedef MoFEMField_multiIndex::index<BitFieldId_mi_tag>::type field_set_by_id;
  //find field
  const field_set_by_id &set_id = moabFields.get<BitFieldId_mi_tag>();
  field_set_by_id::iterator miit = set_id.find(id);
  if(miit == set_id.end()) {
    SETERRQ(PETSC_COMM_SELF,MOFEM_NOT_FOUND,"field not found");
  }
  //ents in the field meshset
  Range ents_of_id_meshset;
  rval = moab.get_entities_by_handle(miit->meshset,ents_of_id_meshset,false); CHKERR_PETSC(rval);
  if(verb>5) {
    PetscSynchronizedPrintf(
      comm,"ents in field %s meshset %d\n",miit->get_name().c_str(),ents_of_id_meshset.size()
    );
  }
  //create dofsMoabField
  Range::iterator eit = ents_of_id_meshset.begin();
  for(;eit!=ents_of_id_meshset.end();eit++) {
    // check if ent is in ref meshset
    RefMoFEMEntity_multiIndex::index<Ent_mi_tag>::type::iterator miit_ref_ent;
    miit_ref_ent = refinedEntities.get<Ent_mi_tag>().find(*eit);
    if(miit_ref_ent==refinedEntities.get<Ent_mi_tag>().end()) {
      RefMoFEMEntity ref_ent(moab,*eit);
      if(ref_ent.get_BitRefLevel().none()) {
        continue; // not on any mesh and not in database
      }
      cerr << ref_ent << endl;
      cerr << "bit level " << ref_ent.get_BitRefLevel() << endl;
      SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"database inconsistency");
    }
    // create mofem entity linked to ref ent
    MoFEMEntity_multiIndex::iterator e_miit;
    try {
      e_miit = entsMoabField.find(MoFEMEntity(moab,&*miit,&*miit_ref_ent).get_global_unique_id());
    } catch (MoFEMException const &e) {
      SETERRQ(PETSC_COMM_SELF,e.errorCode,e.errorMessage);
    } catch (const std::exception& ex) {
      ostringstream ss;
      ss << ex.what() << endl;
      SETERRQ(PETSC_COMM_SELF,MOFEM_OPERATION_UNSUCCESSFUL,ss.str().c_str());
    }
    if(e_miit == entsMoabField.end()) {
      ApproximationOrder order = -1;
      rval = moab.tag_set_data(miit->th_AppOrder,&*eit,1,&order); CHKERR_PETSC(rval);
      pair<MoFEMEntity_multiIndex::iterator,bool> p_e_miit;
      try {
        MoFEMEntity moabent(moab,&*miit,&*miit_ref_ent);
        p_e_miit = entsMoabField.insert(moabent);
      } catch (MoFEMException const &e) {
        SETERRQ(PETSC_COMM_SELF,e.errorCode,e.errorMessage);
      } catch (const std::exception& ex) {
        ostringstream ss;
        ss << ex.what() << endl;
        SETERRQ(PETSC_COMM_SELF,MOFEM_OPERATION_UNSUCCESSFUL,ss.str().c_str());
      }
      if(!p_e_miit.second) SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"data inconsistency");
      bool success = entsMoabField.modify(p_e_miit.first,MoFEMEntity_change_order(moab,-1));
      if(!success) SETERRQ(PETSC_COMM_SELF,MOFEM_OPERATION_UNSUCCESSFUL,"modification unsuccessful");
      e_miit = p_e_miit.first;
    }
    // insert dofmoabent into mofem databse
    int DD = 0;
    int oo = 0;
    // loop orders
    for(;oo<=e_miit->get_max_order();oo++) {
      //loop nb. dofs at order oo
      for(int dd = 0;dd<e_miit->get_order_nb_dofs_diff(oo);dd++) {
        //loop rank
        for(int rr = 0;rr<e_miit->get_max_rank();rr++,DD++) {
          pair<DofMoFEMEntity_multiIndex::iterator,bool> d_miit;
          try {
            DofMoFEMEntity mdof(&*(e_miit),oo,rr,DD);
            d_miit = dofsMoabField.insert(mdof);
            if(d_miit.second) {
              dof_counter[d_miit.first->get_ent_type()]++;
              bool success = dofsMoabField.modify(d_miit.first,DofMoFEMEntity_active_change(true));
              if(!success) SETERRQ(PETSC_COMM_SELF,MOFEM_OPERATION_UNSUCCESSFUL,"modification unsuccessful");
            }
            //check ent
            if(d_miit.first->get_ent()!=e_miit->get_ent()) {
              SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"data inconsistency");
            }
            if(d_miit.first->get_ent_type()!=e_miit->get_ent_type()) {
              SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"data inconsistency");
            }
            if(d_miit.first->get_id()!=e_miit->get_id()) {
              SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"data inconsistency");
            }
            //check dof
            if(d_miit.first->get_dof_order()!=oo) {
              ostringstream ss;
              ss << "data inconsistency!" << endl;
              ss << "should be " << mdof << endl;
              ss << "but is " << *d_miit.first << endl;
              SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,ss.str().c_str());
            }
            if(d_miit.first->get_max_order()!=e_miit->get_max_order()) {
              SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"data inconsistency");
            }
          } catch (MoFEMException const &e) {
            SETERRQ(PETSC_COMM_SELF,e.errorCode,e.errorMessage);
          }
        }
      }
    }
    if(DD != e_miit->get_max_rank()*e_miit->get_order_nb_dofs(e_miit->get_max_order())) {
      SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"data inconsistency");
    }
  }
  PetscFunctionReturn(0);
}
PetscErrorCode Core::build_fields(int verb) {
  PetscFunctionBegin;
  if(verb==-1) verb = verbose;
  typedef MoFEMField_multiIndex::index<BitFieldId_mi_tag>::type field_set_by_id;
  field_set_by_id &set_id = moabFields.get<BitFieldId_mi_tag>();
  field_set_by_id::iterator miit = set_id.begin();
  for(;miit!=set_id.end();miit++) {
    map<EntityType,int> dof_counter;
    if (verb > 0) {
      PetscSynchronizedPrintf(comm,"Build Field %s (rank %d)\n",miit->get_name().c_str(),rAnk);
    }
    switch (miit->get_space()) {
      case NOFIELD:
      ierr = dofs_NoField(miit->get_id(),dof_counter,verb); CHKERRQ(ierr);
      break;
      case L2:
      case H1:
      case HCURL:
      case HDIV:
      ierr = dofs_L2H1HcurlHdiv(miit->get_id(),dof_counter,verb); CHKERRQ(ierr);
      break;
      default:
      SETERRQ(PETSC_COMM_SELF,MOFEM_NOT_IMPLEMENTED,"not implemented");
    }
    if (verb > 0) {
      int _dof_counter_ = 0;
      for (map<EntityType,int>::iterator it = dof_counter.begin();it!=dof_counter.end();it++) {
        switch (it->first) {
          case MBVERTEX:
          PetscSynchronizedPrintf(comm,"nb added dofs (vertices) %d\n",it->second);
          break;
          case MBEDGE:
          PetscSynchronizedPrintf(comm,"nb added dofs (edges) %d\n",it->second);
          break;
          case MBTRI:
          PetscSynchronizedPrintf(comm,"nb added dofs (triangles) %d\n",it->second);
          break;
          case MBQUAD:
          PetscSynchronizedPrintf(comm,"nb added dofs (quads) %d\n",it->second);
          break;
          case MBTET:
          PetscSynchronizedPrintf(comm,"nb added dofs (tets) %d\n",it->second);
          break;
          case MBPRISM:
          PetscSynchronizedPrintf(comm,"nb added dofs (prisms) %d\n",it->second);
          break;
          case MBENTITYSET:
          PetscSynchronizedPrintf(comm,"nb added dofs (meshsets) %d\n",it->second);
          break;
          default:
          SETERRQ(PETSC_COMM_SELF,MOFEM_NOT_IMPLEMENTED,"not implemented");
        }
        _dof_counter_ += it->second;
      }
      PetscSynchronizedPrintf(comm,"nb added dofs %d\n",_dof_counter_);
    }
  }
  *build_MoFEM = 1<<0;
  if(verbose>0) {
    PetscSynchronizedPrintf(comm,"Nb. dofs %u\n",dofsMoabField.size());
    PetscSynchronizedFlush(comm,PETSC_STDOUT);
  }
  PetscFunctionReturn(0);
  //return 0;
}
PetscErrorCode Core::list_dofs_by_field_name(const string &field_name) const {
  PetscFunctionBegin;
  DofMoFEMEntity_multiIndex::index<FieldName_mi_tag>::type::iterator dit,hi_dit;
  dit = dofsMoabField.get<FieldName_mi_tag>().lower_bound(field_name);
  hi_dit = dofsMoabField.get<FieldName_mi_tag>().upper_bound(field_name);
  for(;dit!=hi_dit;dit++) {
    ostringstream ss;
    ss << "rank " << rAnk << " ";
    ss << *dit << endl;
    PetscSynchronizedPrintf(comm,ss.str().c_str());
  }
  PetscSynchronizedFlush(comm,PETSC_STDOUT);
  PetscFunctionReturn(0);
}
PetscErrorCode Core::list_fields() const {
  PetscFunctionBegin;
  typedef MoFEMField_multiIndex::index<BitFieldId_mi_tag>::type field_set_by_id;
  const field_set_by_id &set_id = moabFields.get<BitFieldId_mi_tag>();
  field_set_by_id::iterator miit = set_id.begin();
  for(;miit!=set_id.end();miit++) {
    ostringstream ss;
    ss << *miit << endl;
    PetscSynchronizedPrintf(comm,ss.str().c_str());
  }
  PetscSynchronizedFlush(comm,PETSC_STDOUT);
  PetscFunctionReturn(0);
}
PetscErrorCode Core::add_finite_element(const string &MoFEMFiniteElement_name,enum MoFEMTypes bh) {
  PetscFunctionBegin;
  *build_MoFEM &= 1<<0;
  typedef MoFEMFiniteElement_multiIndex::index<FiniteElement_name_mi_tag>::type finiteElements_by_name;
  finiteElements_by_name &MoFEMFiniteElement_name_set = finiteElements.get<FiniteElement_name_mi_tag>();
  finiteElements_by_name::iterator it_MoFEMFiniteElement = MoFEMFiniteElement_name_set.find(MoFEMFiniteElement_name);
  if(bh == MF_EXCL) {
    if(it_MoFEMFiniteElement!=MoFEMFiniteElement_name_set.end()) {
      SETERRQ1(PETSC_COMM_SELF,MOFEM_NOT_FOUND,"this < %s > is there",MoFEMFiniteElement_name.c_str());
    }
  } else {
    if(it_MoFEMFiniteElement!=MoFEMFiniteElement_name_set.end()) PetscFunctionReturn(0);
  }
  EntityHandle meshset;
  rval = moab.create_meshset(MESHSET_SET|MESHSET_TRACK_OWNER,meshset); CHKERR_PETSC(rval);
  //id
  BitFEId id = getFEShift();
  rval = moab.tag_set_data(th_FEId,&meshset,1,&id); CHKERR_PETSC(rval);
  //id name
  void const* tag_data[] = { MoFEMFiniteElement_name.c_str() };
  int tag_sizes[1]; tag_sizes[0] = MoFEMFiniteElement_name.size();
  rval = moab.tag_set_by_ptr(th_FEName,&meshset,1,tag_data,tag_sizes); CHKERR_PETSC(rval);
  //add MoFEMFiniteElement
  pair<MoFEMFiniteElement_multiIndex::iterator,bool> p = finiteElements.insert(MoFEMFiniteElement(moab,meshset));
  if(!p.second) SETERRQ(PETSC_COMM_SELF,MOFEM_OPERATION_UNSUCCESSFUL,"MoFEMFiniteElement not inserted");
  if(verbose>0) {
    ostringstream ss;
    ss << "add finite element: " << MoFEMFiniteElement_name << endl;
    PetscPrintf(comm,ss.str().c_str());
    //list_finiteElements();
  }
  PetscFunctionReturn(0);
}
PetscErrorCode Core::modify_finite_element_adjacency_table(const string &MoFEMFiniteElement_name,const EntityType type,ElementAdjacencyFunct function) {
  PetscFunctionBegin;
  *build_MoFEM &= 1<<0;
  typedef MoFEMFiniteElement_multiIndex::index<FiniteElement_name_mi_tag>::type finiteElements_by_name;
  finiteElements_by_name &MoFEMFiniteElement_name_set = finiteElements.get<FiniteElement_name_mi_tag>();
  finiteElements_by_name::iterator it_MoFEMFiniteElement = MoFEMFiniteElement_name_set.find(MoFEMFiniteElement_name);
  if(it_MoFEMFiniteElement==MoFEMFiniteElement_name_set.end()) SETERRQ(PETSC_COMM_SELF,MOFEM_NOT_FOUND,"this MoFEMFiniteElement is there");
  const_cast<MoFEMFiniteElement*>(&*it_MoFEMFiniteElement)->element_adjacency_table[type] = function;
  PetscFunctionReturn(0);
}
PetscErrorCode Core::modify_finite_element_add_field_data(const string &MoFEMFiniteElement_name,const string &name_data) {
  PetscFunctionBegin;
  *build_MoFEM &= 1<<0;
  typedef MoFEMFiniteElement_multiIndex::index<FiniteElement_name_mi_tag>::type finiteElements_by_name;
  finiteElements_by_name &MoFEMFiniteElement_name_set = finiteElements.get<FiniteElement_name_mi_tag>();
  finiteElements_by_name::iterator it_MoFEMFiniteElement = MoFEMFiniteElement_name_set.find(MoFEMFiniteElement_name);
  if(it_MoFEMFiniteElement==MoFEMFiniteElement_name_set.end()) SETERRQ(PETSC_COMM_SELF,MOFEM_NOT_FOUND,"this MoFEMFiniteElement is there");
  try {
    bool success = MoFEMFiniteElement_name_set.modify(it_MoFEMFiniteElement,EntMoFEMFiniteElement_change_bit_add(get_BitFieldId(name_data)));
    if(!success) SETERRQ(PETSC_COMM_SELF,MOFEM_OPERATION_UNSUCCESSFUL,"modification unsuccessful");
  } catch (MoFEMException const &e) {
    SETERRQ(PETSC_COMM_SELF,e.errorCode,e.errorMessage);
  }
  PetscFunctionReturn(0);
}
PetscErrorCode Core::modify_finite_element_add_field_row(const string &MoFEMFiniteElement_name,const string &name_row) {
  PetscFunctionBegin;
  *build_MoFEM &= 1<<0;
  typedef MoFEMFiniteElement_multiIndex::index<FiniteElement_name_mi_tag>::type finiteElements_by_name;
  finiteElements_by_name &MoFEMFiniteElement_name_set = finiteElements.get<FiniteElement_name_mi_tag>();
  finiteElements_by_name::iterator it_MoFEMFiniteElement = MoFEMFiniteElement_name_set.find(MoFEMFiniteElement_name);
  if(it_MoFEMFiniteElement==MoFEMFiniteElement_name_set.end()) SETERRQ1(PETSC_COMM_SELF,MOFEM_NOT_FOUND,"this < %s > is not there",MoFEMFiniteElement_name.c_str());
  try {
    bool success = MoFEMFiniteElement_name_set.modify(it_MoFEMFiniteElement,MoFEMFiniteElement_row_change_bit_add(get_BitFieldId(name_row)));
    if(!success) SETERRQ(PETSC_COMM_SELF,MOFEM_OPERATION_UNSUCCESSFUL,"modification unsuccessful");
  } catch (MoFEMException const &e) {
    SETERRQ(PETSC_COMM_SELF,e.errorCode,e.errorMessage);
  }
  PetscFunctionReturn(0);
}
PetscErrorCode Core::modify_finite_element_add_field_col(const string &MoFEMFiniteElement_name,const string &name_col) {
  PetscFunctionBegin;
  *build_MoFEM &= 1<<0;
  typedef MoFEMFiniteElement_multiIndex::index<FiniteElement_name_mi_tag>::type finiteElements_by_name;
  finiteElements_by_name &MoFEMFiniteElement_name_set = finiteElements.get<FiniteElement_name_mi_tag>();
  finiteElements_by_name::iterator it_MoFEMFiniteElement = MoFEMFiniteElement_name_set.find(MoFEMFiniteElement_name);
  if(it_MoFEMFiniteElement==MoFEMFiniteElement_name_set.end()) SETERRQ(PETSC_COMM_SELF,MOFEM_OPERATION_UNSUCCESSFUL,"this MoFEMFiniteElement is there");
  try {
    bool success = MoFEMFiniteElement_name_set.modify(it_MoFEMFiniteElement,MoFEMFiniteElement_col_change_bit_add(get_BitFieldId(name_col)));
    if(!success) SETERRQ(PETSC_COMM_SELF,MOFEM_OPERATION_UNSUCCESSFUL,"modification unsuccessful");
  } catch (MoFEMException const &e) {
    SETERRQ(PETSC_COMM_SELF,e.errorCode,e.errorMessage);
  }
  PetscFunctionReturn(0);
}
PetscErrorCode Core::modify_finite_element_off_field_data(const string &MoFEMFiniteElement_name,const string &name_data) {
  PetscFunctionBegin;
  *build_MoFEM &= 1<<0;
  typedef MoFEMFiniteElement_multiIndex::index<FiniteElement_name_mi_tag>::type finiteElements_by_name;
  finiteElements_by_name &MoFEMFiniteElement_name_set = finiteElements.get<FiniteElement_name_mi_tag>();
  finiteElements_by_name::iterator it_MoFEMFiniteElement = MoFEMFiniteElement_name_set.find(MoFEMFiniteElement_name);
  if(it_MoFEMFiniteElement==MoFEMFiniteElement_name_set.end()) SETERRQ(PETSC_COMM_SELF,MOFEM_NOT_FOUND,"this MoFEMFiniteElement is there");
  try {
    bool success = MoFEMFiniteElement_name_set.modify(it_MoFEMFiniteElement,EntMoFEMFiniteElement_change_bit_off(get_BitFieldId(name_data)));
    if(!success) SETERRQ(PETSC_COMM_SELF,MOFEM_OPERATION_UNSUCCESSFUL,"modification unsuccessful");
  } catch (MoFEMException const &e) {
    SETERRQ(PETSC_COMM_SELF,e.errorCode,e.errorMessage);
  }
  PetscFunctionReturn(0);
}
PetscErrorCode Core::modify_finite_element_off_field_row(const string &MoFEMFiniteElement_name,const string &name_row) {
  PetscFunctionBegin;
  *build_MoFEM &= 1<<0;
  typedef MoFEMFiniteElement_multiIndex::index<FiniteElement_name_mi_tag>::type finiteElements_by_name;
  finiteElements_by_name &MoFEMFiniteElement_name_set = finiteElements.get<FiniteElement_name_mi_tag>();
  finiteElements_by_name::iterator it_MoFEMFiniteElement = MoFEMFiniteElement_name_set.find(MoFEMFiniteElement_name);
  if(it_MoFEMFiniteElement==MoFEMFiniteElement_name_set.end()) SETERRQ1(PETSC_COMM_SELF,MOFEM_NOT_FOUND,"this < %s > is not there",MoFEMFiniteElement_name.c_str());
  try {
    bool success = MoFEMFiniteElement_name_set.modify(it_MoFEMFiniteElement,MoFEMFiniteElement_row_change_bit_off(get_BitFieldId(name_row)));
    if(!success) SETERRQ(PETSC_COMM_SELF,MOFEM_OPERATION_UNSUCCESSFUL,"modification unsuccessful");
  } catch (MoFEMException const &e) {
    SETERRQ(PETSC_COMM_SELF,e.errorCode,e.errorMessage);
  }
  PetscFunctionReturn(0);
}
PetscErrorCode Core::modify_finite_element_off_field_col(const string &MoFEMFiniteElement_name,const string &name_col) {
  PetscFunctionBegin;
  *build_MoFEM &= 1<<0;
  typedef MoFEMFiniteElement_multiIndex::index<FiniteElement_name_mi_tag>::type finiteElements_by_name;
  finiteElements_by_name &MoFEMFiniteElement_name_set = finiteElements.get<FiniteElement_name_mi_tag>();
  finiteElements_by_name::iterator it_MoFEMFiniteElement = MoFEMFiniteElement_name_set.find(MoFEMFiniteElement_name);
  if(it_MoFEMFiniteElement==MoFEMFiniteElement_name_set.end()) SETERRQ(PETSC_COMM_SELF,MOFEM_NOT_FOUND,"this MoFEMFiniteElement is there");
  try {
    bool success = MoFEMFiniteElement_name_set.modify(it_MoFEMFiniteElement,MoFEMFiniteElement_col_change_bit_off(get_BitFieldId(name_col)));
    if(!success) SETERRQ(PETSC_COMM_SELF,MOFEM_OPERATION_UNSUCCESSFUL,"modification unsuccessful");
  } catch (MoFEMException const &e) {
    SETERRQ(PETSC_COMM_SELF,e.errorCode,e.errorMessage);
  }
  PetscFunctionReturn(0);
}
BitFEId Core::get_BitFEId(const string& name) const {
  typedef MoFEMFiniteElement_multiIndex::index<FiniteElement_name_mi_tag>::type finiteElements_by_name;
  const finiteElements_by_name& set = finiteElements.get<FiniteElement_name_mi_tag>();
  finiteElements_by_name::iterator miit = set.find(name);
  if(miit==set.end()) THROW_AT_LINE(("finite element < "+name+" > not found (top tip: check spelling)").c_str());
  return miit->get_id();
}
string Core::get_BitFEId_name(const BitFEId id) const {
  typedef MoFEMFiniteElement_multiIndex::index<BitFEId_mi_tag>::type finiteElements_by_id;
  const finiteElements_by_id& set = finiteElements.get<BitFEId_mi_tag>();
  finiteElements_by_id::iterator miit = set.find(id);
  assert(miit!=set.end());
  return miit->get_name();
}
EntityHandle Core::get_finite_element_meshset(const BitFEId id) const {
  typedef MoFEMFiniteElement_multiIndex::index<BitFEId_mi_tag>::type finiteElements_by_id;
  const finiteElements_by_id& set = finiteElements.get<BitFEId_mi_tag>();
  finiteElements_by_id::iterator miit = set.find(id);
  if(miit==set.end()) THROW_AT_LINE("finite element not found");
  return miit->meshset;
}
EntityHandle Core::get_finite_element_meshset(const string& name) const {
  return get_finite_element_meshset(get_BitFEId(name));
}
PetscErrorCode Core::list_finite_elements() const {
  PetscFunctionBegin;
  typedef MoFEMFiniteElement_multiIndex::index<BitFEId_mi_tag>::type finiteElements_by_id;
  const finiteElements_by_id &BitFEId_set = finiteElements.get<BitFEId_mi_tag>();
  finiteElements_by_id::iterator miit = BitFEId_set.begin();
  for(;miit!=BitFEId_set.end();miit++) {
    ostringstream ss;
    ss << *miit << endl;
    PetscSynchronizedPrintf(comm,ss.str().c_str());
  }
  PetscSynchronizedFlush(comm,PETSC_STDOUT);
  PetscFunctionReturn(0);
}
PetscErrorCode Core::add_problem(const BitProblemId id,const string& name) {
  PetscFunctionBegin;
  EntityHandle meshset;
  rval = moab.create_meshset(MESHSET_SET|MESHSET_TRACK_OWNER,meshset); CHKERR_PETSC(rval);
  rval = moab.tag_set_data(th_ProblemId,&meshset,1,&id); CHKERR_PETSC(rval);
  void const* tag_data[] = { name.c_str() };
  int tag_sizes[1]; tag_sizes[0] = name.size();
  rval = moab.tag_set_by_ptr(th_ProblemName,&meshset,1,tag_data,tag_sizes); CHKERR_PETSC(rval);
  //create entry
  pair<MoFEMProblem_multiIndex::iterator,bool> p = moFEMProblems.insert(MoFEMProblem(moab,meshset));
  NOT_USED(p);
  assert(p.second);
  if(verbose>0) {
    ostringstream ss;
    ss << "add problem: " << name << endl;
    PetscPrintf(comm,ss.str().c_str());
  }
  PetscFunctionReturn(0);
}
PetscErrorCode Core::add_problem(const string& name,enum MoFEMTypes bh,int verb) {
  PetscFunctionBegin;
  if(verb==-1) verb = verbose;
  typedef MoFEMProblem_multiIndex::index<Problem_mi_tag>::type mofem_problems_by_name;
  const mofem_problems_by_name& set = moFEMProblems.get<Problem_mi_tag>();
  mofem_problems_by_name::iterator miit = set.find(name);
  if(miit==set.end()) {
    BitProblemId id = getProblemShift();
    ierr = add_problem(id,name); CHKERRQ(ierr);
  } else if(bh == MF_EXCL) {
    SETERRQ1(PETSC_COMM_SELF,MOFEM_NOT_FOUND,"problem is in database %s",name.c_str());
  }
  PetscFunctionReturn(0);
}
PetscErrorCode Core::delete_problem(const string& name) {
  PetscFunctionBegin;
  typedef MoFEMProblem_multiIndex::index<Problem_mi_tag>::type mofem_problems_by_name;
  mofem_problems_by_name &mofem_problems_set = moFEMProblems.get<Problem_mi_tag>();
  mofem_problems_by_name::iterator p_miit = mofem_problems_set.find(name);
  if(p_miit == mofem_problems_set.end()) {
    SETERRQ1(PETSC_COMM_SELF,1,"no such problem like < %s >",name.c_str());
  }
  EntityHandle meshset = p_miit->meshset;
  mofem_problems_set.erase(p_miit);
  rval = moab.delete_entities(&meshset,1); CHKERR_PETSC(rval);
  PetscFunctionReturn(0);
}
BitProblemId Core::get_BitProblemId(const string& name) const {
  typedef MoFEMProblem_multiIndex::index<Problem_mi_tag>::type mofem_problems_by_name;
  const mofem_problems_by_name& set = moFEMProblems.get<Problem_mi_tag>();
  mofem_problems_by_name::iterator miit = set.find(name);
  return miit->get_id();
}
PetscErrorCode Core::list_problem() const {
  PetscFunctionBegin;
  typedef MoFEMProblem_multiIndex::index<BitProblemId_mi_tag>::type problem_set_by_id;
  const problem_set_by_id &set_id = moFEMProblems.get<BitProblemId_mi_tag>();
  problem_set_by_id::iterator miit = set_id.begin();
  for(;miit!=set_id.end();miit++) {
    ostringstream ss;
    ss << *miit << endl;
    PetscPrintf(comm,ss.str().c_str());
  }
  PetscFunctionReturn(0);
}
PetscErrorCode Core::add_ents_to_finite_element_by_EDGEs(const Range& edges,const BitFEId id) {
  PetscFunctionBegin;
  *build_MoFEM &= 1<<0;
  const EntityHandle idm = get_finite_element_meshset(id);
  rval = moab.add_entities(idm,edges.subset_by_type(MBEDGE)); CHKERR_PETSC(rval);
  PetscFunctionReturn(0);
}
PetscErrorCode Core::add_ents_to_finite_element_by_EDGEs(const Range& edges,const string &name) {
  PetscFunctionBegin;
  *build_MoFEM &= 1<<0;
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
  *build_MoFEM &= 1<<0;
  const EntityHandle idm = get_finite_element_meshset(id);
  ierr = seed_finite_elements(vert.subset_by_type(MBVERTEX)); CHKERRQ(ierr);
  rval = moab.add_entities(idm,vert.subset_by_type(MBVERTEX)); CHKERR_PETSC(rval);
  PetscFunctionReturn(0);
}
PetscErrorCode Core::add_ents_to_finite_element_by_VERTICEs(const Range& vert,const string &name) {
  PetscFunctionBegin;
  *build_MoFEM &= 1<<0;
  try {
    ierr = add_ents_to_finite_element_by_VERTICEs(vert,get_BitFEId(name));  CHKERRQ(ierr);
  } catch (MoFEMException const &e) {
    SETERRQ(PETSC_COMM_SELF,e.errorCode,e.errorMessage);
  }
  PetscFunctionReturn(0);
}
PetscErrorCode Core::add_ents_to_finite_element_by_TRIs(const Range& tris,const BitFEId id) {
  PetscFunctionBegin;
  *build_MoFEM &= 1<<0;
  const EntityHandle idm = get_finite_element_meshset(id);
  ierr = seed_finite_elements(tris.subset_by_type(MBTRI)); CHKERRQ(ierr);
  rval = moab.add_entities(idm,tris.subset_by_type(MBTRI)); CHKERR_PETSC(rval);
  PetscFunctionReturn(0);
}
PetscErrorCode Core::add_ents_to_finite_element_by_TRIs(const Range& tris,const string &name) {
  PetscFunctionBegin;
  *build_MoFEM &= 1<<0;
  try {
    ierr = add_ents_to_finite_element_by_TRIs(tris,get_BitFEId(name));  CHKERRQ(ierr);
  } catch (MoFEMException const &e) {
    SETERRQ(PETSC_COMM_SELF,e.errorCode,e.errorMessage);
  }
  PetscFunctionReturn(0);
}
PetscErrorCode Core::add_ents_to_finite_element_by_TETs(const EntityHandle meshset,const BitFEId id,const bool recursive) {
  PetscFunctionBegin;
  *build_MoFEM &= 1<<0;
  EntityHandle idm = no_handle;
  try {
    idm = get_finite_element_meshset(id);
  } catch (MoFEMException const &e) {
    SETERRQ(PETSC_COMM_SELF,e.errorCode,e.errorMessage);
  }
  Range tets;
  rval = moab.get_entities_by_type(meshset,MBTET,tets,recursive); CHKERR_PETSC(rval);
  rval = moab.add_entities(idm,tets); CHKERR_PETSC(rval);
  PetscFunctionReturn(0);
}
PetscErrorCode Core::add_ents_to_finite_element_by_TETs(const Range& tets,const BitFEId id) {
  PetscFunctionBegin;
  *build_MoFEM &= 1<<0;
  const EntityHandle idm = get_finite_element_meshset(id);
  rval = moab.add_entities(idm,tets.subset_by_type(MBTET)); CHKERR_PETSC(rval);
  PetscFunctionReturn(0);
}
PetscErrorCode Core::add_ents_to_finite_element_by_TETs(const Range& tets,const string &name) {
  PetscFunctionBegin;
  *build_MoFEM &= 1<<0;
  try {
    ierr = add_ents_to_finite_element_by_TETs(tets,get_BitFEId(name));  CHKERRQ(ierr);
  } catch (MoFEMException const &e) {
    SETERRQ(PETSC_COMM_SELF,e.errorCode,e.errorMessage);
  }
  PetscFunctionReturn(0);
}
PetscErrorCode Core::add_ents_to_finite_element_by_TETs(const EntityHandle meshset,const string &name,const bool recursive) {
  PetscFunctionBegin;
  *build_MoFEM &= 1<<0;
  try {
    ierr = add_ents_to_finite_element_by_TETs(meshset,get_BitFEId(name),recursive);  CHKERRQ(ierr);
  } catch (MoFEMException const &e) {
    SETERRQ(PETSC_COMM_SELF,e.errorCode,e.errorMessage);
  }
  PetscFunctionReturn(0);
}

PetscErrorCode Core::add_ents_to_finite_element_by_PRISMs(const EntityHandle meshset,const BitFEId id,const bool recursive) {
  PetscFunctionBegin;
  *build_MoFEM &= 1<<0;
  EntityHandle idm = no_handle;
  try {
    idm = get_finite_element_meshset(id);
  } catch (MoFEMException const &e) {
    SETERRQ(PETSC_COMM_SELF,e.errorCode,e.errorMessage);
  }
  Range prisms;
  rval = moab.get_entities_by_type(meshset,MBPRISM,prisms,recursive); CHKERR_PETSC(rval);
  rval = moab.add_entities(idm,prisms); CHKERR_PETSC(rval);
  PetscFunctionReturn(0);
}
PetscErrorCode Core::add_ents_to_finite_element_by_PRISMs(const Range& tets,const BitFEId id) {
  PetscFunctionBegin;
  *build_MoFEM &= 1<<0;
  const EntityHandle idm = get_finite_element_meshset(id);
  rval = moab.add_entities(idm,tets.subset_by_type(MBPRISM)); CHKERR_PETSC(rval);
  PetscFunctionReturn(0);
}
PetscErrorCode Core::add_ents_to_finite_element_by_PRISMs(const Range& prims,const string &name) {
  PetscFunctionBegin;
  *build_MoFEM &= 1<<0;
  try {
    ierr = add_ents_to_finite_element_by_PRISMs(prims,get_BitFEId(name));  CHKERRQ(ierr);
  } catch (MoFEMException const &e) {
    SETERRQ(PETSC_COMM_SELF,e.errorCode,e.errorMessage);
  }
  PetscFunctionReturn(0);
}
PetscErrorCode Core::add_ents_to_finite_element_by_PRISMs(const EntityHandle meshset,const string &name,const bool recursive) {
  PetscFunctionBegin;
  *build_MoFEM &= 1<<0;
  try {
    ierr = add_ents_to_finite_element_by_PRISMs(meshset,get_BitFEId(name),recursive);  CHKERRQ(ierr);
  } catch (MoFEMException const &e) {
    SETERRQ(PETSC_COMM_SELF,e.errorCode,e.errorMessage);
  }
  PetscFunctionReturn(0);
}
PetscErrorCode Core::add_ents_to_finite_element_EntType_by_bit_ref(const BitRefLevel &bit,const string &name,EntityType type,int verb) {
  PetscFunctionBegin;
  ierr = add_ents_to_finite_element_EntType_by_bit_ref(bit,BitRefLevel().set(),name,type,verb); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
PetscErrorCode Core::add_ents_to_finite_element_EntType_by_bit_ref(const BitRefLevel &bit,const BitRefLevel &mask,const string &name,EntityType type,int verb) {
  PetscFunctionBegin;
  try {
    if(verb==-1) verb = verbose;
    *build_MoFEM &= 1<<0;
    const BitFEId id = get_BitFEId(name);
    const EntityHandle idm = get_finite_element_meshset(id);
    typedef RefMoFEMElement_multiIndex::index<EntType_mi_tag>::type refMoabFE_by_type;
    refMoabFE_by_type &ref_MoFEMFiniteElement = refinedFiniteElements.get<EntType_mi_tag>();
    refMoabFE_by_type::iterator miit = ref_MoFEMFiniteElement.lower_bound(type);
    refMoabFE_by_type::iterator hi_miit = ref_MoFEMFiniteElement.upper_bound(type);
    if(verb > 1) {
      PetscSynchronizedPrintf(comm,"nb. ref elements in database %d\n",distance(miit,hi_miit));
    }
    int nb_add_FEs = 0;
    for(;miit!=hi_miit;miit++) {
      BitRefLevel bit2 = miit->get_BitRefLevel();
      //check if all bits in mask are ib fe bit2
      //if((miit->get_BitRefLevel()&bit)!=bit) continue;
      if((bit2&mask) != bit2) continue;
      if((bit2&bit).any()) {
        EntityHandle ent = miit->get_ref_ent();
        rval = moab.add_entities(idm,&ent,1); CHKERR_PETSC(rval);
        nb_add_FEs++;
      }
    }
    if(verb > 0) {
      ostringstream ss;
      ss << "Add Nb. FEs " << nb_add_FEs << " form BitRef " << bit << endl;
      PetscSynchronizedPrintf(comm,"%s",ss.str().c_str());
      PetscSynchronizedFlush(comm,PETSC_STDOUT);
    }
  } catch (MoFEMException const &e) {
    SETERRQ(PETSC_COMM_SELF,e.errorCode,e.errorMessage);
  }
  PetscFunctionReturn(0);
}
PetscErrorCode Core::add_ents_to_finite_element_by_MESHSET(const EntityHandle meshset,const string& name,const bool recursive) {
  PetscFunctionBegin;
  *build_MoFEM &= 1<<0;
  const BitFEId id = get_BitFEId(name);
  const EntityHandle idm = get_finite_element_meshset(id);
  if(recursive==false){
    rval = moab.add_entities(idm,&meshset,1); CHKERR_PETSC(rval);
  } else {
    Range meshsets;
    rval = moab.get_entities_by_type(meshset,MBENTITYSET,meshsets,false); CHKERR_PETSC(rval);
    rval = moab.add_entities(idm,meshsets); CHKERR_PETSC(rval);
  }
  PetscFunctionReturn(0);
}
PetscErrorCode Core::modify_problem_add_finite_element(const string &name_problem,const string &MoFEMFiniteElement_name) {
  PetscFunctionBegin;
  try {
    typedef MoFEMProblem_multiIndex::index<Problem_mi_tag>::type mofem_problems_by_name;
    mofem_problems_by_name& set = moFEMProblems.get<Problem_mi_tag>();
    mofem_problems_by_name::iterator miit = set.find(name_problem);
    ostringstream ss;
    ss << name_problem;
    if(miit==set.end()) SETERRQ1(PETSC_COMM_SELF,MOFEM_NOT_FOUND,"this problem <%s> is not there",ss.str().c_str());
    BitFEId f_id = get_BitFEId(MoFEMFiniteElement_name);
    bool success = set.modify(miit,problem_MoFEMFiniteElement_change_bit_add(f_id));
    if(!success) SETERRQ(PETSC_COMM_SELF,MOFEM_OPERATION_UNSUCCESSFUL,"modification unsuccessful");
  } catch (MoFEMException const &e) {
    SETERRQ(PETSC_COMM_SELF,e.errorCode,e.errorMessage);
  }
  PetscFunctionReturn(0);
}
PetscErrorCode Core::modify_problem_ref_level_add_bit(const string &name_problem,const BitRefLevel &bit) {
  PetscFunctionBegin;
  typedef MoFEMProblem_multiIndex::index<Problem_mi_tag>::type mofem_problems_by_name;
  mofem_problems_by_name& set = moFEMProblems.get<Problem_mi_tag>();
  mofem_problems_by_name::iterator miit = set.find(name_problem);
  ostringstream ss;
  ss << name_problem;
  if(miit==set.end()) SETERRQ1(PETSC_COMM_SELF,MOFEM_NOT_FOUND,"this problem <%s> is there",ss.str().c_str());
  bool success = set.modify(miit,problem_change_ref_level_bit_add(bit));
  if(!success) SETERRQ(PETSC_COMM_SELF,1,"modification unsuccessful");
  PetscFunctionReturn(0);
}
PetscErrorCode Core::modify_problem_ref_level_set_bit(const string &name_problem,const BitRefLevel &bit) {
  PetscFunctionBegin;
  typedef MoFEMProblem_multiIndex::index<Problem_mi_tag>::type mofem_problems_by_name;
  mofem_problems_by_name& set = moFEMProblems.get<Problem_mi_tag>();
  mofem_problems_by_name::iterator miit = set.find(name_problem);
  ostringstream ss;
  ss << name_problem;
  if(miit==set.end()) SETERRQ1(PETSC_COMM_SELF,MOFEM_NOT_FOUND,"this problem <%s> is there",ss.str().c_str());
  bool success = set.modify(miit,problem_change_ref_level_bit_set(bit));
  if(!success) SETERRQ(PETSC_COMM_SELF,MOFEM_OPERATION_UNSUCCESSFUL,"modification unsuccessful");
  PetscFunctionReturn(0);
}
PetscErrorCode Core::modify_problem_dof_mask_ref_level_set_bit(const string &name_problem,const BitRefLevel &bit) {
  PetscFunctionBegin;
  typedef MoFEMProblem_multiIndex::index<Problem_mi_tag>::type mofem_problems_by_name;
  mofem_problems_by_name& set = moFEMProblems.get<Problem_mi_tag>();
  mofem_problems_by_name::iterator miit = set.find(name_problem);
  ostringstream ss;
  ss << name_problem;
  if(miit==set.end()) SETERRQ1(PETSC_COMM_SELF,MOFEM_NOT_FOUND,"this problem <%s> is there",ss.str().c_str());
  bool success = set.modify(miit,problem_change_ref_level_bit_dof_mask_set(bit));
  if(!success) SETERRQ(PETSC_COMM_SELF,MOFEM_OPERATION_UNSUCCESSFUL,"modification unsuccessful");
  PetscFunctionReturn(0);
}
PetscErrorCode Core::build_finite_element_data_dofs(EntMoFEMFiniteElement &ent_fe,int verb) {
  PetscFunctionBegin;

  if(!(*build_MoFEM)&(1<<0)) SETERRQ(PETSC_COMM_SELF,MOFEM_NOT_FOUND,"fields not build");
  FEDofMoFEMEntity_multiIndex &data_dofs = const_cast<FEDofMoFEMEntity_multiIndex&>(ent_fe.data_dofs);
  data_dofs.clear(); //clear data dofs multi-index //FIXME should be cleaned when dofs are cleaned form datasets

  DofMoFEMEntity_multiIndex_active_view data_view;
  ierr = ent_fe.get_MoFEMFiniteElement_data_dof_view(dofsMoabField,data_view,Interface::UNION); CHKERRQ(ierr);
  DofMoFEMEntity_multiIndex_active_view::nth_index<1>::type::iterator viit_data,hi_viit_data;

  //loops over active dofs only
  viit_data = data_view.get<1>().lower_bound(1);
  hi_viit_data = data_view.get<1>().upper_bound(1);
  for(;viit_data!=hi_viit_data;viit_data++) {
    try {
      switch((*viit_data)->get_space()) {
        case H1:
        case HDIV:
        case HCURL:
        case L2:
        case NOFIELD:
        {
          SideNumber *side_number_ptr = ent_fe.get_side_number_ptr(moab,(*viit_data)->get_ent());
          //add dofs to finite element multi_index database
          data_dofs.get<Unique_mi_tag>().insert(data_dofs.end(),boost::make_tuple(side_number_ptr,&**viit_data));
        }
        break;
        default:
        SETERRQ(PETSC_COMM_SELF,MOFEM_NOT_IMPLEMENTED,"not implemented");
      }
    } catch (MoFEMException const &e) {
      SETERRQ(PETSC_COMM_SELF,e.errorCode,e.errorMessage);
    }
  }
  viit_data = data_view.get<1>().lower_bound(1);
  if(data_dofs.size()!=(unsigned int)distance(viit_data,hi_viit_data)) {
    SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"data inconsistency");
  }
  PetscFunctionReturn(0);
}
PetscErrorCode Core::build_finite_element_uids_view(EntMoFEMFiniteElement &ent_fe,int verb) {
  PetscFunctionBegin;
  if(!(*build_MoFEM)&(1<<0)) SETERRQ(PETSC_COMM_SELF,MOFEM_NOT_FOUND,"fields not build");
  typedef MoFEMField_multiIndex::index<BitFieldId_mi_tag>::type field_by_id;
  typedef RefMoFEMEntity_multiIndex::index<Ent_mi_tag>::type ref_ent_by_ent;
  typedef DofMoFEMEntity_multiIndex::index<Composite_Name_And_Ent_mi_tag>::type dof_set_type;
  field_by_id &moabFields_by_id = moabFields.get<BitFieldId_mi_tag>();
  dof_set_type& dof_set = dofsMoabField.get<Composite_Name_And_Ent_mi_tag>();
  //get id of mofem fields for row, col and data
  enum IntLoop { Row = 0,Col,Data,Last };
  BitFieldId FEAdj_fields[Last] = {
    ent_fe.get_BitFieldId_row(),
    ent_fe.get_BitFieldId_col(),
    ent_fe.get_BitFieldId_data()
  };
  //get refinment level
  const BitRefLevel& bit_ref_MoFEMFiniteElement = ent_fe.get_BitRefLevel();
  Range tets,faces,edges,nodes,meshsets,adj_ents,ent_ents;
  Range::iterator eit_eit;
  DofMoFEMEntity_multiIndex_uid_view* MoFEMFiniteElement_dof_uid_view[Last] = {
    &ent_fe.row_dof_view,
    &ent_fe.col_dof_view,
    &ent_fe.data_dof_view
  };
  unsigned int nb_view_dofs[Last];
  for(int ss = 0;ss<Last;ss++) {
    MoFEMFiniteElement_dof_uid_view[ss]->clear();
    nb_view_dofs[ss] = 0;
  }
  //lopp over all fields in database
  for(unsigned int ii = 0;ii<BitFieldId().size();ii++) {
    // common field id for Row, Col and Data
    BitFieldId id_common = 0;
    //check if the field (ii) is added to finite element
    for(int ss = 0;ss<Last;ss++) {
      id_common |= FEAdj_fields[ss]&BitFieldId().set(ii);
    }
    if( id_common.none() ) continue;
    //find in database data associated with the field (ii)
    field_by_id::iterator miit = moabFields_by_id.find(BitFieldId().set(ii));
    if(miit==moabFields_by_id.end()) {
      SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"data inconsistency");
    }
    //resolve entities on element, those entities are used to build tag with dof
    //uids on finite element tag
    ierr = ent_fe.get_element_adjacency(moab,&*miit,adj_ents); CHKERRQ(ierr);
    //loop over adjacent to finite entities, and find dofs on those entities
    //this part is to build MoFEMFiniteElement_dof_uid_view
    Range::iterator eit2 = adj_ents.begin();
    for(;eit2!=adj_ents.end();eit2++) {
      ref_ent_by_ent::iterator ref_ent_miit = refinedEntities.get<Ent_mi_tag>().find(*eit2);
      if(ref_ent_miit==refinedEntities.get<Ent_mi_tag>().end()) {
        RefMoFEMEntity ref_ent(moab,*eit2);
        if(!ref_ent.get_BitRefLevel().any()) {
          continue;
        }
        cerr << adj_ents << endl;
        cerr << ent_fe << endl;
        cerr << "bit level " << ent_fe.get_BitRefLevel() << endl;
        cerr << ref_ent << endl;
        cerr << "bit level " << ref_ent.get_BitRefLevel() << endl;
        SETERRQ(PETSC_COMM_SELF,MOFEM_NOT_FOUND,"ref ent not in database");
      }
      const BitRefLevel& bit_ref_ent = ref_ent_miit->get_BitRefLevel();
      if(!(bit_ref_MoFEMFiniteElement&bit_ref_ent).any()) {
        ostringstream ss;
        ss << "top tip: check if you seed mesh with the elements for bit ref level1" << endl;
        ss << "inconsitency in database entity" << " type "
        << moab.type_from_handle(*eit2) << " bits ENT " << bit_ref_ent << endl;
        ss << "inconsitency in database entity" << " type "
        << moab.type_from_handle(ent_fe.get_ent()) << " bits FE  " << bit_ref_MoFEMFiniteElement << endl;
        SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,ss.str().c_str());
      }
      dof_set_type::iterator ents_miit2 = dof_set.lower_bound(boost::make_tuple(miit->get_name_ref(),ref_ent_miit->get_ref_ent()));
      dof_set_type::iterator ents_hi_miit2 = dof_set.upper_bound(boost::make_tuple(miit->get_name_ref(),ref_ent_miit->get_ref_ent()));
      for(int ss = 0;ss<Last;ss++) {
        if( !(FEAdj_fields[ss].test(ii)) ) continue;
        dof_set_type::iterator ents_miit3 = ents_miit2;
        for(;ents_miit3!=ents_hi_miit2;ents_miit3++) {
          MoFEMFiniteElement_dof_uid_view[ss]->insert(&*ents_miit3);
          nb_view_dofs[ss]++;
        }
      }
    }
  }

  for(int ss = 0;ss<Last;ss++) {
    if(nb_view_dofs[ss] != MoFEMFiniteElement_dof_uid_view[ss]->size()) {
      SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"data insonsistency");
    }
  }

  if(verb>2) {
    ostringstream ss;
    ss << "add: FE data"  << endl << ent_fe << endl;
    //rows
    DofMoFEMEntity_multiIndex_uid_view MoFEMFiniteElement_row_dof_uid_view;
    ierr = ent_fe.get_MoFEMFiniteElement_row_dof_view(dofsMoabField,MoFEMFiniteElement_row_dof_uid_view); CHKERRQ(ierr);
    DofMoFEMEntity_multiIndex_uid_view::iterator miit_row = MoFEMFiniteElement_row_dof_uid_view.begin();
    ss << "rows dofs" << endl;
    for(;miit_row!=MoFEMFiniteElement_row_dof_uid_view.end();miit_row++) ss << **miit_row << endl;
    //cols
    DofMoFEMEntity_multiIndex_uid_view MoFEMFiniteElement_col_dof_uid_view;
    ierr = ent_fe.get_MoFEMFiniteElement_col_dof_view(dofsMoabField,MoFEMFiniteElement_col_dof_uid_view); CHKERRQ(ierr);
    DofMoFEMEntity_multiIndex_uid_view::iterator miit_col = MoFEMFiniteElement_col_dof_uid_view.begin();
    ss << "cols dofs" << endl;
    for(;miit_col!=MoFEMFiniteElement_col_dof_uid_view.end();miit_col++) ss << **miit_col << endl;
    PetscPrintf(comm,ss.str().c_str());
  }

  PetscFunctionReturn(0);
}
PetscErrorCode Core::build_finite_elements(int verb) {
  PetscFunctionBegin;
  if(verb==-1) verb = verbose;

  typedef RefMoFEMElement_multiIndex::index<Ent_mi_tag>::type ref_MoFEMFiniteElement_by_ent;
  MoFEMFiniteElement_multiIndex::iterator MoFEMFiniteElement_miit = finiteElements.begin();
  // loop Finite Elements
  for(;MoFEMFiniteElement_miit!=finiteElements.end();MoFEMFiniteElement_miit++) {

    if(verb>0) PetscPrintf(comm,"Build Finite Elements %s\n",MoFEMFiniteElement_miit->get_name().c_str());
    //get finite element meshset
    EntityHandle meshset = get_finite_element_meshset(MoFEMFiniteElement_miit->get_id());
    // get entities from finite element meshset // if meshset
    Range MoFEMFiniteElement_ents;
    rval = moab.get_entities_by_handle(meshset,MoFEMFiniteElement_ents,false); CHKERR_PETSC(rval);
    //loop meshset Ents and add finite elements
    Range::iterator eit = MoFEMFiniteElement_ents.begin();
    for(;eit!=MoFEMFiniteElement_ents.end();eit++) {
      // note: iterator is a wrapper
      // check if is in refinedFiniteElements database
      ref_MoFEMFiniteElement_by_ent::iterator ref_MoFEMFiniteElement_miit;
      ref_MoFEMFiniteElement_miit = refinedFiniteElements.get<Ent_mi_tag>().find(*eit);
      if(ref_MoFEMFiniteElement_miit == refinedFiniteElements.get<Ent_mi_tag>().end()) {
        ostringstream ss;
        ss << "ref MoFEMFiniteElement not in database ent = " << *eit;
        ss << " type " << moab.type_from_handle(*eit);
        ss << " " << *MoFEMFiniteElement_miit;
        SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,ss.str().c_str());
      }
      EntMoFEMFiniteElement ent_fe(moab,ref_MoFEMFiniteElement_miit->get_RefMoFEMElement(),&*MoFEMFiniteElement_miit);
      pair<EntMoFEMFiniteElement_multiIndex::iterator,bool> p = finiteElementsMoFEMEnts.insert(ent_fe);
      ierr = build_finite_element_uids_view(const_cast<EntMoFEMFiniteElement&>(*p.first),verb); CHKERRQ(ierr);
      ierr = build_finite_element_data_dofs(const_cast<EntMoFEMFiniteElement&>(*p.first),verb); CHKERRQ(ierr);
    }

  }

  if(verb>0) {
    PetscSynchronizedPrintf(comm,"Nb. FEs %u\n",finiteElementsMoFEMEnts.size());
    PetscSynchronizedFlush(comm,PETSC_STDOUT);
    typedef EntMoFEMFiniteElement_multiIndex::index<BitFEId_mi_tag>::type MoFEMFiniteElement_by_id;
    MoFEMFiniteElement_by_id &MoFEMFiniteElements = finiteElementsMoFEMEnts.get<BitFEId_mi_tag>();
    MoFEMFiniteElement_multiIndex::iterator id_MoFEMFiniteElement = finiteElements.begin();
    for(;id_MoFEMFiniteElement!=finiteElements.end();id_MoFEMFiniteElement++) {
      MoFEMFiniteElement_by_id::iterator miit = MoFEMFiniteElements.lower_bound(id_MoFEMFiniteElement->get_id());
      MoFEMFiniteElement_by_id::iterator hi_miit = MoFEMFiniteElements.upper_bound(id_MoFEMFiniteElement->get_id());
      int count = distance(miit,hi_miit);
      ostringstream ss;
      ss << *id_MoFEMFiniteElement << " Nb. FEs " << count << endl;
      PetscSynchronizedPrintf(comm,ss.str().c_str());
      PetscSynchronizedFlush(comm,PETSC_STDOUT);
    }
  }

  *build_MoFEM |= 1<<1;
  PetscFunctionReturn(0);
}
PetscErrorCode Core::build_adjacencies(const Range &ents,int verb) {
  PetscFunctionBegin;
  if(verb==-1) verb = verbose;
  if(!(*build_MoFEM)&(1<<0)) SETERRQ(PETSC_COMM_SELF,MOFEM_NOT_FOUND,"field not build");
  if(!(*build_MoFEM)&(1<<1)) SETERRQ(PETSC_COMM_SELF,MOFEM_NOT_FOUND,"fe not build");
  //typedef MoFEMEntity_multiIndex::index<Unique_mi_tag>::type ents_by_uid;
  EntMoFEMFiniteElement_multiIndex::iterator fit = finiteElementsMoFEMEnts.begin();
  for(;fit!=finiteElementsMoFEMEnts.end();fit++) {
    if(!ents.empty()) {
      if(ents.find(fit->get_ent())==ents.end()) {
        continue;
      }
    }
    GlobalUId ent_uid = UId(0);
    DofMoFEMEntity_multiIndex_uid_view::iterator rvit;
    rvit = fit->row_dof_view.begin();
    for(;rvit!=fit->row_dof_view.end();rvit++) {
      if( ent_uid == (*rvit)->get_MoFEMEntity_ptr()->get_global_unique_id()) continue;
      ent_uid = (*rvit)->get_MoFEMEntity_ptr()->get_global_unique_id();
      pair<MoFEMEntityEntMoFEMFiniteElementAdjacencyMap_multiIndex::iterator,bool> p;
      p = entFEAdjacencies.insert(MoFEMEntityEntMoFEMFiniteElementAdjacencyMap((*rvit)->get_MoFEMEntity_ptr(),&*fit));
      bool success = entFEAdjacencies.modify(p.first,MoFEMEntityEntMoFEMFiniteElementAdjacencyMap_change_ByWhat(BYROW));
      if(!success) SETERRQ(PETSC_COMM_SELF,MOFEM_OPERATION_UNSUCCESSFUL,"modification unsuccessful");
    }
    ent_uid = UId(0);
    DofMoFEMEntity_multiIndex_uid_view::iterator cvit;
    cvit = fit->col_dof_view.begin();
    for(;cvit!=fit->col_dof_view.end();cvit++) {
      if( ent_uid == (*cvit)->get_MoFEMEntity_ptr()->get_global_unique_id()) continue;
      ent_uid = (*cvit)->get_MoFEMEntity_ptr()->get_global_unique_id();
      pair<MoFEMEntityEntMoFEMFiniteElementAdjacencyMap_multiIndex::iterator,bool> p;
      p = entFEAdjacencies.insert(MoFEMEntityEntMoFEMFiniteElementAdjacencyMap((*cvit)->get_MoFEMEntity_ptr(),&*fit));
      bool success = entFEAdjacencies.modify(p.first,MoFEMEntityEntMoFEMFiniteElementAdjacencyMap_change_ByWhat(BYCOL));
      if(!success) SETERRQ(PETSC_COMM_SELF,MOFEM_OPERATION_UNSUCCESSFUL,"modification unsuccessful");
    }
    ent_uid = UId(0);
    DofMoFEMEntity_multiIndex_uid_view::iterator dvit;
    dvit = fit->data_dof_view.begin();
    for(;dvit!=fit->data_dof_view.end();dvit++) {
      if( ent_uid == (*dvit)->get_MoFEMEntity_ptr()->get_global_unique_id()) continue;
      ent_uid = (*dvit)->get_MoFEMEntity_ptr()->get_global_unique_id();
      pair<MoFEMEntityEntMoFEMFiniteElementAdjacencyMap_multiIndex::iterator,bool> p;
      p = entFEAdjacencies.insert(MoFEMEntityEntMoFEMFiniteElementAdjacencyMap((*dvit)->get_MoFEMEntity_ptr(),&*fit));
      bool success = entFEAdjacencies.modify(p.first,MoFEMEntityEntMoFEMFiniteElementAdjacencyMap_change_ByWhat(BYDATA));
      if(!success) SETERRQ(PETSC_COMM_SELF,MOFEM_OPERATION_UNSUCCESSFUL,"modification unsuccessful");
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
  *build_MoFEM |= 1<<2;
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
  MoFEMEntityEntMoFEMFiniteElementAdjacencyMap_multiIndex::iterator miit = entFEAdjacencies.begin();
  for(;miit!=entFEAdjacencies.end();miit++) {
    ostringstream ss;
    ss << *miit << endl;
    PetscSynchronizedPrintf(comm,ss.str().c_str());
  }
  PetscSynchronizedFlush(comm,PETSC_STDOUT);
  PetscFunctionReturn(0);
}
PetscErrorCode Core::partition_finite_elements(
  const string &name,bool part_from_moab,int low_proc,int hi_proc,int verb) {
  PetscFunctionBegin;
  if(verb==-1) verb = verbose;
  if(!(*build_MoFEM&(1<<0))) SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"fields not build");
  if(!(*build_MoFEM&(1<<1))) SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"FEs not build");
  if(!(*build_MoFEM&(1<<2))) SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"entFEAdjacencies not build");
  if(!(*build_MoFEM&(1<<3))) SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"partitions not build");
  if(!(*build_MoFEM&(1<<4))) SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"partitions moFEMProblems not build");
  if(low_proc == -1) low_proc = rAnk;
  if(hi_proc == -1) hi_proc = rAnk;
  typedef MoFEMProblem_multiIndex::index<Problem_mi_tag>::type mofem_problems_by_name;
  //find p_miit
  mofem_problems_by_name &moFEMProblems_set = moFEMProblems.get<Problem_mi_tag>();
  mofem_problems_by_name::iterator p_miit = moFEMProblems_set.find(name);
  if(p_miit == moFEMProblems_set.end()) {
    SETERRQ1(
      PETSC_COMM_SELF,1,"problem < %s > not found (top tip: check spelling)",name.c_str()
    );
  }
  NumeredMoFEMFiniteElement_multiIndex& numeredFiniteElements
  = const_cast<NumeredMoFEMFiniteElement_multiIndex&>(p_miit->numeredFiniteElements);
  numeredFiniteElements.clear();
  //MoFEMFiniteElement set
  EntMoFEMFiniteElement_multiIndex::iterator miit2 = finiteElementsMoFEMEnts.begin();
  EntMoFEMFiniteElement_multiIndex::iterator hi_miit2 = finiteElementsMoFEMEnts.end();
  for(;miit2!=hi_miit2;miit2++) {
    if((miit2->get_id()&p_miit->get_BitFEId()).none()) continue; // if element is not part of prblem
    if((miit2->get_BitRefLevel()&p_miit->get_BitRefLevel())!=p_miit->get_BitRefLevel()) continue; // if entity is not problem refinment level
    NumeredMoFEMFiniteElement numered_fe(&*miit2);
    FENumeredDofMoFEMEntity_multiIndex &rows_dofs = numered_fe.rows_dofs;
    FENumeredDofMoFEMEntity_multiIndex &cols_dofs = numered_fe.cols_dofs;
    rows_dofs.clear();
    cols_dofs.clear();
    {
      NumeredDofMoFEMEntity_multiIndex_uid_view_ordered rows_view;
      //rows_view
      ierr = miit2->get_MoFEMFiniteElement_row_dof_view(p_miit->numered_dofs_rows,rows_view,Interface::UNION); CHKERRQ(ierr);
      NumeredDofMoFEMEntity_multiIndex_uid_view_ordered::iterator viit_rows;
      if(part_from_moab) {
        int proc = miit2->get_owner_proc();
        NumeredMoFEMFiniteElement_change_part(proc).operator()(numered_fe);
      } else {
        vector<int> parts(sIze,0);
        viit_rows = rows_view.begin();
        for(;viit_rows!=rows_view.end();viit_rows++) {
          parts[(*viit_rows)->part]++;
        }
        vector<int>::iterator pos = max_element(parts.begin(),parts.end());
        unsigned int max_part = distance(parts.begin(),pos);
        NumeredMoFEMFiniteElement_change_part(max_part).operator()(numered_fe);
      }
      if(
        (numered_fe.get_part()>=(unsigned int)low_proc)&&
        (numered_fe.get_part()<=(unsigned int)hi_proc)
      ) {
        //rows element dof multi-indices
        viit_rows = rows_view.begin();
        for(;viit_rows!=rows_view.end();viit_rows++) {
          try {
            SideNumber *side_number_ptr = miit2->get_side_number_ptr(moab,(*viit_rows)->get_ent());
            rows_dofs.insert(boost::make_tuple(side_number_ptr,&**viit_rows));
          } catch (MoFEMException const &e) {
            SETERRQ(PETSC_COMM_SELF,e.errorCode,e.errorMessage);
          }
        }
        //cols_views
        NumeredDofMoFEMEntity_multiIndex_uid_view_ordered cols_view;
        ierr = miit2->get_MoFEMFiniteElement_col_dof_view(
          p_miit->numered_dofs_cols,cols_view,Interface::UNION
        ); CHKERRQ(ierr);
        //cols element dof multi-indices
        NumeredDofMoFEMEntity_multiIndex_uid_view_ordered::iterator viit_cols;;
        viit_cols = cols_view.begin();
        for(;viit_cols!=cols_view.end();viit_cols++) {
          try {
            SideNumber *side_number_ptr = miit2->get_side_number_ptr(moab,(*viit_cols)->get_ent());
            cols_dofs.insert(boost::make_tuple(side_number_ptr,&**viit_cols));
          } catch (MoFEMException const &e) {
            SETERRQ(PETSC_COMM_SELF,e.errorCode,e.errorMessage);
          }
        }
      }
      pair<NumeredMoFEMFiniteElement_multiIndex::iterator,bool> p;
      p = numeredFiniteElements.insert(numered_fe);
      if(!p.second) {
        SETERRQ(PETSC_COMM_SELF,MOFEM_NOT_FOUND,"element is there");
      }
      if(verb>1) {
        ostringstream ss;
        ss << *p_miit << endl;
        ss << *p.first << endl;
        typedef FENumeredDofMoFEMEntity_multiIndex::index<Unique_mi_tag>::type FENumeredDofMoFEMEntity_multiIndex_by_Unique_mi_tag;
        FENumeredDofMoFEMEntity_multiIndex_by_Unique_mi_tag::iterator miit = p.first->rows_dofs.get<Unique_mi_tag>().begin();
        for(;miit!= p.first->rows_dofs.get<Unique_mi_tag>().end();miit++) ss << "rows: " << *miit << endl;
        miit = p.first->cols_dofs.get<Unique_mi_tag>().begin();
        for(;miit!=p.first->cols_dofs.get<Unique_mi_tag>().end();miit++) ss << "cols: " << *miit << endl;
        PetscSynchronizedPrintf(comm,ss.str().c_str());
      }
    }
  }
  if(verb>0) {
    typedef NumeredMoFEMFiniteElement_multiIndex::index<FiniteElement_Part_mi_tag>::type NumeredMoFEMFiniteElement_multiIndex_by_part;
    NumeredMoFEMFiniteElement_multiIndex_by_part::iterator MoFEMFiniteElement_miit = numeredFiniteElements.get<FiniteElement_Part_mi_tag>().lower_bound(rAnk);
    NumeredMoFEMFiniteElement_multiIndex_by_part::iterator hi_MoMoFEMFiniteElement_miitFEMFE_miit = numeredFiniteElements.get<FiniteElement_Part_mi_tag>().upper_bound(rAnk);
    int count = distance(MoFEMFiniteElement_miit,hi_MoMoFEMFiniteElement_miitFEMFE_miit);
    ostringstream ss;
    ss << *p_miit;
    ss << " Nb. elems " << count << " on proc " << rAnk << endl;
    PetscSynchronizedPrintf(comm,ss.str().c_str());
    PetscSynchronizedFlush(comm,PETSC_STDOUT);
  }
  *build_MoFEM |= 1<<5;
  PetscFunctionReturn(0);
}
PetscErrorCode Core::partition_ghost_dofs(const string &name,int verb) {
  PetscFunctionBegin;
  if(verb==-1) verb = verbose;
  if(!(*build_MoFEM&(1<<0))) SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"fields not build");
  if(!(*build_MoFEM&(1<<1))) SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"FEs not build");
  if(!(*build_MoFEM&(1<<2))) SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"entFEAdjacencies not build");
  if(!(*build_MoFEM&(1<<3))) SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"moFEMProblems not build");
  if(!(*build_MoFEM&(1<<4))) SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"partitions moFEMProblems not build");
  if(!(*build_MoFEM&(1<<5))) SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"partitions finite elements not build");
  typedef MoFEMProblem_multiIndex::index<Problem_mi_tag>::type mofem_problems_by_name;
  //find p_miit
  mofem_problems_by_name &moFEMProblems_set = moFEMProblems.get<Problem_mi_tag>();
  mofem_problems_by_name::iterator p_miit = moFEMProblems_set.find(name);
  //
  DofIdx &nb_row_ghost_dofs = *((DofIdx*)p_miit->tag_ghost_nbdof_data_row);
  DofIdx &nb_col_ghost_dofs = *((DofIdx*)p_miit->tag_ghost_nbdof_data_col);
  //
  nb_row_ghost_dofs = 0;
  nb_col_ghost_dofs = 0;
  if(sIze>1) {
    NumeredDofMoFEMEntity_multiIndex_uid_view_ordered ghost_idx_col_view,ghost_idx_row_view;
    NumeredMoFEMFiniteElement_multiIndex::index<FiniteElement_Part_mi_tag>::type::iterator fe_it,hi_fe_it;
    fe_it = p_miit->numeredFiniteElements.get<FiniteElement_Part_mi_tag>().lower_bound(rAnk);
    hi_fe_it = p_miit->numeredFiniteElements.get<FiniteElement_Part_mi_tag>().upper_bound(rAnk);
    for(;fe_it!=hi_fe_it;fe_it++) {
      typedef FENumeredDofMoFEMEntity_multiIndex::iterator dof_it;
      if(fe_it->rows_dofs.size()>0) {
        dof_it rowdofit,hi_rowdofit;
        rowdofit = fe_it->rows_dofs.begin();
        hi_rowdofit = fe_it->rows_dofs.end();
        for(;rowdofit!=hi_rowdofit;rowdofit++) {
          if(rowdofit->get_part()==(unsigned int)rAnk) continue;
          ghost_idx_row_view.insert(rowdofit->get_NumeredDofMoFEMEntity_ptr());
        }
      }
      if(fe_it->cols_dofs.size()>0) {
        dof_it coldofit,hi_coldofit;
        coldofit = fe_it->cols_dofs.begin();
        hi_coldofit = fe_it->cols_dofs.end();
        for(;coldofit!=hi_coldofit;coldofit++) {
          if(coldofit->get_part()==(unsigned int)rAnk) continue;
          ghost_idx_col_view.insert(coldofit->get_NumeredDofMoFEMEntity_ptr());
        }
      }
    }
    DofIdx *nb_ghost_dofs[2] = { &nb_col_ghost_dofs, &nb_row_ghost_dofs };
    DofIdx nb_local_dofs[2] = { *((DofIdx*)p_miit->tag_local_nbdof_data_col), *((DofIdx*)p_miit->tag_local_nbdof_data_row) };
    NumeredDofMoFEMEntity_multiIndex_uid_view_ordered *ghost_idx_view[2] = { &ghost_idx_col_view, &ghost_idx_row_view };
    typedef NumeredDofMoFEMEntity_multiIndex::index<Unique_mi_tag>::type NumeredDofMoFEMEntitys_by_unique_id;
    NumeredDofMoFEMEntitys_by_unique_id *dof_by_uid_no_const[2] = {
      const_cast<NumeredDofMoFEMEntitys_by_unique_id*>(&p_miit->numered_dofs_cols.get<Unique_mi_tag>()),
      const_cast<NumeredDofMoFEMEntitys_by_unique_id*>(&p_miit->numered_dofs_rows.get<Unique_mi_tag>())
    };
    for(int ss = 0;ss<2;ss++) {
      NumeredDofMoFEMEntity_multiIndex_uid_view_ordered::iterator ghost_idx_miit = ghost_idx_view[ss]->begin();
      for(;ghost_idx_miit!=ghost_idx_view[ss]->end();ghost_idx_miit++) {
        NumeredDofMoFEMEntitys_by_unique_id::iterator diit = dof_by_uid_no_const[ss]->find((*ghost_idx_miit)->get_global_unique_id());
        if(diit->petsc_local_dof_idx!=(DofIdx)-1) {
          SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"inconsistent data, ghost dof already set");
        }
        bool success = dof_by_uid_no_const[ss]->modify(diit,NumeredDofMoFEMEntity_local_idx_change(nb_local_dofs[ss]++));
        if(!success) SETERRQ(PETSC_COMM_SELF,MOFEM_OPERATION_UNSUCCESSFUL,"modification unsuccessful");
        (*nb_ghost_dofs[ss])++;
      }
    }
  }
  if(verb>0) {
    ostringstream ss;
    ss << "partition_ghost_col_dofs: rank = " << rAnk
    << " FEs col ghost dofs "<< *p_miit
    << " Nb. col ghost dof " << p_miit->get_nb_ghost_dofs_col()
    << " Nb. local dof " << p_miit->get_nb_local_dofs_col() << endl;
    ss << "partition_ghost_row_dofs: rank = " << rAnk
    << " FEs row ghost dofs "<< *p_miit
    << " Nb. row ghost dof " << p_miit->get_nb_ghost_dofs_row()
    << " Nb. local dof " << p_miit->get_nb_local_dofs_row() << endl;
    if(verb>1) {
      NumeredDofMoFEMEntity_multiIndex::iterator miit_dd_col = p_miit->numered_dofs_cols.begin();
      for(;miit_dd_col!=p_miit->numered_dofs_cols.end();miit_dd_col++) {
        if(miit_dd_col->part==(unsigned int)rAnk) continue;
        if(miit_dd_col->petsc_local_dof_idx==(DofIdx)-1) continue;
        ss<<*miit_dd_col<<endl;
      }
    }
    PetscSynchronizedPrintf(comm,ss.str().c_str());
    PetscSynchronizedFlush(comm,PETSC_STDOUT);
  }
  *build_MoFEM |= 1<<6;
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
    RefMoFEMEntity_multiIndex::index<Ent_mi_tag>::type::iterator
      eiit = refinedEntities.get<Ent_mi_tag>().find(*eit);
    if(eiit == refinedEntities.get<Ent_mi_tag>().end())  {
      SETERRQ(PETSC_COMM_SELF,MOFEM_NOT_FOUND,"entity is not in database");
    }
    if(eiit->get_BitRefLevel().none()) continue;
    pair<RefMoFEMElement_multiIndex::iterator,bool> p_MoFEMFiniteElement;
    switch (eiit->get_ent_type()) {
      case MBVERTEX:
      p_MoFEMFiniteElement = refinedFiniteElements.insert(ptrWrapperRefMoFEMElement(new RefMoFEMElement_VERTEX(moab,&*eiit)));
      break;
      case MBEDGE:
      p_MoFEMFiniteElement = refinedFiniteElements.insert(ptrWrapperRefMoFEMElement(new RefMoFEMElement_EDGE(moab,&*eiit)));
      break;
      case MBTRI:
      p_MoFEMFiniteElement = refinedFiniteElements.insert(ptrWrapperRefMoFEMElement(new RefMoFEMElement_TRI(moab,&*eiit)));
      break;
      case MBTET:
      p_MoFEMFiniteElement = refinedFiniteElements.insert(ptrWrapperRefMoFEMElement(new RefMoFEMElement_TET(moab,&*eiit)));
      break;
      case MBPRISM:
      p_MoFEMFiniteElement = refinedFiniteElements.insert(ptrWrapperRefMoFEMElement(new RefMoFEMElement_PRISM(moab,&*eiit)));
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
  rval = moab.get_entities_by_dimension(meshset,2,ents2d,false); CHKERR_PETSC(rval);
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
      RefMoFEMEntity ref_ent(moab,*tit);
      bitset<8> ent_pstat(ref_ent.get_pstatus());
      ent_pstat.flip(0);
      pair<RefMoFEMEntity_multiIndex::iterator,bool> p_ent = refinedEntities.insert(ref_ent);
      if(debug > 0) {
        ierr = test_moab(moab,*tit); CHKERRQ(ierr);
      }
      bool success = refinedEntities.modify(p_ent.first,RefMoFEMEntity_change_add_bit(bit));
      if(!success) SETERRQ(PETSC_COMM_SELF,MOFEM_OPERATION_UNSUCCESSFUL,"modification unsuccessful");
      if(verb>2) {
        ostringstream ss;
        ss << *(p_ent.first);
        PetscSynchronizedPrintf(comm,"%s\n",ss.str().c_str());
      }
      pair<RefMoFEMElement_multiIndex::iterator,bool> p_MoFEMFiniteElement;
      switch (p_ent.first->get_ent_type()) {
        case MBVERTEX:
        p_MoFEMFiniteElement = refinedFiniteElements.insert(ptrWrapperRefMoFEMElement(new RefMoFEMElement_VERTEX(moab,&*p_ent.first)));
        seeded_ents.insert(*tit);
        break;
        case MBEDGE:
        p_MoFEMFiniteElement = refinedFiniteElements.insert(ptrWrapperRefMoFEMElement(new RefMoFEMElement_EDGE(moab,&*p_ent.first)));
        seeded_ents.insert(*tit);
        break;
        case MBTRI:
        p_MoFEMFiniteElement = refinedFiniteElements.insert(ptrWrapperRefMoFEMElement(new RefMoFEMElement_TRI(moab,&*p_ent.first)));
        seeded_ents.insert(*tit);
        break;
        case MBTET:
        p_MoFEMFiniteElement = refinedFiniteElements.insert(ptrWrapperRefMoFEMElement(new RefMoFEMElement_TET(moab,&*p_ent.first)));
        seeded_ents.insert(*tit);
        break;
        case MBPRISM:
        p_MoFEMFiniteElement = refinedFiniteElements.insert(ptrWrapperRefMoFEMElement(new RefMoFEMElement_PRISM(moab,&*p_ent.first)));
        break;
        case MBENTITYSET:
        p_MoFEMFiniteElement = refinedFiniteElements.insert(ptrWrapperRefMoFEMElement(new RefMoFEMElement_MESHSET(moab,&*p_ent.first)));
        break;
        default:
        SETERRQ(PETSC_COMM_SELF,MOFEM_NOT_IMPLEMENTED,"not implemented");
      }
      if(verb>3) {
        ostringstream ss;
        ss << *(p_MoFEMFiniteElement.first->get_RefMoFEMElement());
        PetscSynchronizedPrintf(comm,"%s\n",ss.str().c_str());
      }
    }
    if(!seeded_ents.empty()) {
      int dim = moab.dimension_from_handle(seeded_ents[0]);
      for(int dd = 0;dd<dim;dd++) {
        Range ents;
        rval = moab.get_adjacencies(seeded_ents,dd,true,ents,Interface::UNION); CHKERR_PETSC(rval);
        if(dd == 2) {
          // currently only works with triangles
          ents = ents.subset_by_type(MBTRI);
        }
        Range::iterator eit = ents.begin();
        for(;eit!=ents.end();eit++) {
          pair<RefMoFEMEntity_multiIndex::iterator,bool> p_ent = refinedEntities.insert(RefMoFEMEntity(moab,*eit));
          bool success = refinedEntities.modify(p_ent.first,RefMoFEMEntity_change_add_bit(bit));
          if(!success) SETERRQ(PETSC_COMM_SELF,MOFEM_OPERATION_UNSUCCESSFUL,"modification unsuccessful");
          if(verb>2) {
            ostringstream ss;
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
  rval = moab.get_entities_by_dimension(meshset,3,ents3d,false); CHKERR_PETSC(rval);
  ierr = seed_ref_level(ents3d,bit,verb); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
PetscErrorCode Core::seed_ref_level_MESHSET(const EntityHandle meshset,const BitRefLevel &bit,int verb) {
  PetscFunctionBegin;
  if(verb==-1) verb = verbose;
  pair<RefMoFEMEntity_multiIndex::iterator,bool> p_ent = refinedEntities.insert(RefMoFEMEntity(moab,meshset));
  refinedEntities.modify(p_ent.first,RefMoFEMEntity_change_add_bit(bit));
  ptrWrapperRefMoFEMElement pack_fe(new RefMoFEMElement_MESHSET(moab,&*p_ent.first));
  pair<RefMoFEMElement_multiIndex::iterator,bool> p_MoFEMFiniteElement = refinedFiniteElements.insert(pack_fe);
  if(verbose > 0) {
    ostringstream ss;
    ss << "add meshset as ref_ent " << *(p_MoFEMFiniteElement.first->get_RefMoFEMElement()) << endl;
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
  RefMoFEMEntity_multiIndex::iterator ent_it = refinedEntities.begin();
  for(;ent_it!=refinedEntities.end();ent_it++) {
    if(verb>5) {
      cout << ent_it->get_BitRefLevel() << " : ";
    }
    bool success = refinedEntities.modify(ent_it,RefMoFEMEntity_change_right_shift(shift));
    if(!success) SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"inconsistency in data");
    if(verb>5) {
      cout << ent_it->get_BitRefLevel() << endl;
    }
  }
  PetscFunctionReturn(0);
}
PetscErrorCode Core::get_entities_by_type_and_ref_level(const BitRefLevel &bit,const BitRefLevel &mask,const EntityType type,const EntityHandle meshset,int verb) {
  PetscFunctionBegin;
  Range ents;
  ierr = get_entities_by_type_and_ref_level(bit,mask,type,ents,verb); CHKERRQ(ierr);
  rval = moab.add_entities(meshset,ents); CHKERR_PETSC(rval);
  PetscFunctionReturn(0);
}
PetscErrorCode Core::get_entities_by_type_and_ref_level(const BitRefLevel &bit,const BitRefLevel &mask,const EntityType type,Range &ents,int verb) {
  PetscFunctionBegin;
  if(verb==-1) verb = verbose;
  ierr = moab.get_entities_by_type(0,type,ents,false); CHKERRQ(ierr);
  Range::iterator eit = ents.begin();
  for(;eit!=ents.end();) {
    BitRefLevel bit2;
    rval = moab.tag_get_data(th_RefBitLevel,&*eit,1,&bit2); CHKERR_PETSC(rval);
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
PetscErrorCode Core::get_entities_by_ref_level(const BitRefLevel &bit,const BitRefLevel &mask,const EntityHandle meshset) {
  PetscFunctionBegin;
  Range ents;
  ierr = get_entities_by_ref_level(bit,mask,ents); CHKERRQ(ierr);
  rval = moab.add_entities(meshset,ents); CHKERR_PETSC(rval);
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
    rval = moab.tag_get_data(th_RefBitLevel,&*eit,1,&bit2); CHKERR_PETSC(rval);
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
    rval = moab.tag_get_data(th_RefBitLevel,&*eit,1,&bit2); CHKERR_PETSC(rval);
    bit2 |= bit;
    rval = moab.tag_set_data(th_RefBitLevel,&*eit,1,&bit2); CHKERR_PETSC(rval);
  }
  PetscFunctionReturn(0);
}
PetscErrorCode Core::set_ref_level_to_entities(const BitRefLevel &bit,Range &ents) {
  PetscFunctionBegin;
  Range::iterator eit = ents.begin();
  for(;eit!=ents.end();eit++) {
    BitRefLevel bit2;
    rval = moab.tag_get_data(th_RefBitLevel,&*eit,1,&bit2); CHKERR_PETSC(rval);
    bit2 = bit;
    rval = moab.tag_set_data(th_RefBitLevel,&*eit,1,&bit2); CHKERR_PETSC(rval);
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
    cerr << parent << endl;
    cerr << moab.type_from_handle(parent) <<  " " << MBENTITYSET << endl;
  } CHKERR_PETSC(rval);

  typedef RefMoFEMEntity_multiIndex::index<Composite_Ent_And_ParentEntType_mi_tag>::type ref_ents_by_composite;
  ref_ents_by_composite &ref_ents = refinedEntities.get<Composite_Ent_And_ParentEntType_mi_tag>();
  Range::iterator eit = ents.begin();
  for(;eit!=ents.end();eit++) {
    if(verb>2) {
      ostringstream ss;
      ss << "ent " << *eit << endl;;
      PetscPrintf(comm,ss.str().c_str());
    }
    ref_ents_by_composite::iterator miit = ref_ents.lower_bound(boost::make_tuple(*eit,child_type));
    ref_ents_by_composite::iterator hi_miit = ref_ents.upper_bound(boost::make_tuple(*eit,child_type));
    for(;miit!=hi_miit;miit++) {
      if(verb>2) {
	ostringstream ss;
	ss << "any bit " << *miit << endl;;
	PetscPrintf(comm,ss.str().c_str());
      }
      if((miit->get_BitRefLevel()&child_bit).any()) {
	EntityHandle ref_ent = miit->get_ref_ent();
	if(ref_ent == *eit) continue;
	if(ref_ent == 0) {
	  SETERRQ(PETSC_COMM_SELF,MOFEM_IMPOSIBLE_CASE,"this should not happen");
	}
	if(moab.type_from_handle(*eit)==MBENTITYSET) {
	  SETERRQ(PETSC_COMM_SELF,MOFEM_IMPOSIBLE_CASE,"this should not happen");
	}
	rval = moab.add_entities(child,&ref_ent,1); CHKERR_PETSC(rval);
	if(verb>1) {
	  ostringstream ss;
	  ss << "good bit " << *miit << endl;
	  PetscPrintf(comm,ss.str().c_str());
	}
      }
    }
  }
  PetscFunctionReturn(0);
}
PetscErrorCode Core::update_field_meshset_by_entities_children(const BitRefLevel &child_bit,int verb) {
  PetscFunctionBegin;
  MoFEMField_multiIndex::iterator fit = moabFields.begin();
  for(;fit!=moabFields.end();fit++) {
    EntityHandle meshset = fit->get_meshset();
    ierr = update_meshset_by_entities_children(meshset,child_bit,meshset,MBTET,false,verb);  CHKERRQ(ierr);
    ierr = update_meshset_by_entities_children(meshset,child_bit,meshset,MBTRI,false,verb);  CHKERRQ(ierr);
    ierr = update_meshset_by_entities_children(meshset,child_bit,meshset,MBEDGE,false,verb);  CHKERRQ(ierr);
    ierr = update_meshset_by_entities_children(meshset,child_bit,meshset,MBVERTEX,false,verb);  CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}
PetscErrorCode Core::update_field_meshset_by_entities_children(const string name,const BitRefLevel &child_bit,int verb) {
  PetscFunctionBegin;
  MoFEMField_multiIndex::index<FieldName_mi_tag>::type::iterator miit;
  miit = moabFields.get<FieldName_mi_tag>().find(name);
  EntityHandle meshset = miit->get_meshset();
  ierr = update_meshset_by_entities_children(meshset,child_bit,meshset,MBTET,false,verb);  CHKERRQ(ierr);
  ierr = update_meshset_by_entities_children(meshset,child_bit,meshset,MBTRI,false,verb);  CHKERRQ(ierr);
  ierr = update_meshset_by_entities_children(meshset,child_bit,meshset,MBEDGE,false,verb);  CHKERRQ(ierr);
  ierr = update_meshset_by_entities_children(meshset,child_bit,meshset,MBVERTEX,false,verb);  CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
PetscErrorCode Core::update_finite_element_meshset_by_entities_children(const string name,const BitRefLevel &child_bit,const EntityType fe_ent_type,int verb) {
  PetscFunctionBegin;
  typedef MoFEMFiniteElement_multiIndex::index<FiniteElement_name_mi_tag>::type finiteElements_by_name;
  const finiteElements_by_name& set = finiteElements.get<FiniteElement_name_mi_tag>();
  finiteElements_by_name::iterator miit = set.find(name);
  if(miit==set.end()) THROW_AT_LINE(("finite element < "+name+" > not found (top tip: check spelling)").c_str());
  EntityHandle meshset = miit->get_meshset();
  ierr = update_meshset_by_entities_children(meshset,child_bit,meshset,fe_ent_type,false,verb);  CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
PetscErrorCode Core::get_problem_finite_elements_entities(const string &problem_name,const string &fe_name,const EntityHandle meshset) {
  PetscFunctionBegin;
  typedef MoFEMProblem_multiIndex::index<Problem_mi_tag>::type mofem_problems_by_name;
  mofem_problems_by_name &moFEMProblems_set = moFEMProblems.get<Problem_mi_tag>();
  mofem_problems_by_name::iterator p_miit = moFEMProblems_set.find(problem_name);
  if(p_miit == moFEMProblems_set.end()) SETERRQ1(PETSC_COMM_SELF,1,"no such problem like < %s >",problem_name.c_str());
  NumeredMoFEMFiniteElement_multiIndex &numeredFiniteElements = const_cast<NumeredMoFEMFiniteElement_multiIndex&>(p_miit->numeredFiniteElements);
  NumeredMoFEMFiniteElement_multiIndex::index<FiniteElement_name_mi_tag>::type::iterator miit = numeredFiniteElements.get<FiniteElement_name_mi_tag>().lower_bound(fe_name);
  for(;miit!=numeredFiniteElements.get<FiniteElement_name_mi_tag>().upper_bound(fe_name);miit++) {
    EntityHandle ent = miit->get_ent();
    rval = moab.add_entities(meshset,&ent,1); CHKERR_PETSC(rval);
    int part = miit->get_part();
    rval = moab.tag_set_data(th_Part,&ent,1,&part); CHKERR_PETSC(rval);
  }
  PetscFunctionReturn(0);
}
bool Core::check_msId_meshset(const int msId,const CubitBCType cubit_bc_type) {
  CubitMeshSet_multiIndex::index<Composite_Cubit_msId_And_MeshSetType_mi_tag>::type::iterator
    miit = cubitMeshsets.get<Composite_Cubit_msId_And_MeshSetType_mi_tag>().find(boost::make_tuple(msId,cubit_bc_type.to_ulong()));
  if(miit!=cubitMeshsets.get<Composite_Cubit_msId_And_MeshSetType_mi_tag>().end()) {
    return true;
  }
  return false;
}
PetscErrorCode Core::add_cubit_msId(const CubitBCType cubit_bc_type,const int msId) {
  PetscFunctionBegin;
  if(check_msId_meshset(msId,cubit_bc_type)) {
    SETERRQ1(PETSC_COMM_SELF,1,"such cubit meshset is already there",msId);
  }
  try {
    CubitMeshSets cmeshset(moab,cubit_bc_type,msId);
    if((cmeshset.cubit_bc_type&CubitBCType(NODESET|SIDESET|BLOCKSET)).any()) {
      pair<CubitMeshSet_multiIndex::iterator,bool> p = cubitMeshsets.insert(cmeshset);
      if(!p.second) SETERRQ(PETSC_COMM_SELF,MOFEM_OPERATION_UNSUCCESSFUL,"meshset not inserted");
    }
  } catch (MoFEMException const &e) {
    SETERRQ(PETSC_COMM_SELF,e.errorCode,e.errorMessage);
  }
  PetscFunctionReturn(0);
}
PetscErrorCode Core::delete_cubit_msId(const CubitBCType cubit_bc_type,const int msId) {
  PetscFunctionBegin;
  CubitMeshSet_multiIndex::index<Composite_Cubit_msId_And_MeshSetType_mi_tag>::type::iterator
    miit = cubitMeshsets.get<Composite_Cubit_msId_And_MeshSetType_mi_tag>().find(boost::make_tuple(msId,cubit_bc_type.to_ulong()));
  if(miit==cubitMeshsets.get<Composite_Cubit_msId_And_MeshSetType_mi_tag>().end()) {
    SETERRQ1(PETSC_COMM_SELF,1,"such cubit meshset is already there",msId);
  }
  EntityHandle meshset = miit->get_meshset();
  cubitMeshsets.get<Composite_Cubit_msId_And_MeshSetType_mi_tag>().erase(miit);
  rval = moab.delete_entities(&meshset,1); CHKERR_PETSC(rval);
  PetscFunctionReturn(0);
}
PetscErrorCode Core::get_cubit_msId(const int msId,const CubitBCType cubit_bc_type,const CubitMeshSets **cubit_meshset_ptr) {
  PetscFunctionBegin;
  CubitMeshSet_multiIndex::index<Composite_Cubit_msId_And_MeshSetType_mi_tag>::type::iterator
    miit = cubitMeshsets.get<Composite_Cubit_msId_And_MeshSetType_mi_tag>().find(boost::make_tuple(msId,cubit_bc_type.to_ulong()));
  if(miit!=cubitMeshsets.get<Composite_Cubit_msId_And_MeshSetType_mi_tag>().end()) {
    *cubit_meshset_ptr = &*miit;
  } else {
    SETERRQ1(PETSC_COMM_SELF,1,"msId = %d is not there",msId);
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
PetscErrorCode Core::get_cubit_msId_entities_by_dimension(const int msId,const unsigned int cubit_bc_type,
  const int dimension,Range &entities,const bool recursive) {
  PetscFunctionBegin;
  ierr = get_cubit_msId_entities_by_dimension(msId,CubitBCType(cubit_bc_type),dimension,entities,recursive); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
PetscErrorCode Core::get_cubit_msId_entities_by_dimension(const int msId,const unsigned int cubit_bc_type,
  Range &entities,const bool recursive) {
  PetscFunctionBegin;
  ierr = get_cubit_msId_entities_by_dimension(msId,CubitBCType(cubit_bc_type),entities,recursive); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
PetscErrorCode Core::get_cubit_msId_meshset(const int msId,const unsigned int cubit_bc_type,EntityHandle &meshset) {
  PetscFunctionBegin;
  CubitMeshSet_multiIndex::index<Composite_Cubit_msId_And_MeshSetType_mi_tag>::type::iterator
    miit = cubitMeshsets.get<Composite_Cubit_msId_And_MeshSetType_mi_tag>().find(boost::make_tuple(msId,cubit_bc_type));
  if(miit!=cubitMeshsets.get<Composite_Cubit_msId_And_MeshSetType_mi_tag>().end()) {
    meshset = miit->meshset;
  } else {
    SETERRQ1(PETSC_COMM_SELF,1,"msId = %d is not there",msId);
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
    METHOD.fieldsPtr = &moabFields; \
    METHOD.refinedEntitiesPtr = &refinedEntities; \
    METHOD.entitiesPtr = &entsMoabField; \
    METHOD.dofsPtr = &dofsMoabField; \
    METHOD.refinedFiniteElementsPtr = &refinedFiniteElements; \
    METHOD.finiteElementsPtr = &finiteElements; \
    METHOD.finiteElementsEntitiesPtr = &finiteElementsMoFEMEnts; \
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

PetscErrorCode Core::problem_basic_method_preProcess(const string &problem_name,BasicMethod &method,int verb) {
  PetscFunctionBegin;
  if(verb==-1) verb = verbose;
  typedef MoFEMProblem_multiIndex::index<Problem_mi_tag>::type mofem_problems_by_name;
  // find p_miit
  mofem_problems_by_name &moFEMProblems_set = moFEMProblems.get<Problem_mi_tag>();
  mofem_problems_by_name::iterator p_miit = moFEMProblems_set.find(problem_name);
  if(p_miit == moFEMProblems_set.end()) SETERRQ1(PETSC_COMM_SELF,1,"problem is not in database %s",problem_name.c_str());
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
PetscErrorCode Core::problem_basic_method_postProcess(const string &problem_name,BasicMethod &method,int verb) {
  PetscFunctionBegin;
  if(verb==-1) verb = verbose;
  typedef MoFEMProblem_multiIndex::index<Problem_mi_tag>::type mofem_problems_by_name;

  // find p_miit
  mofem_problems_by_name &moFEMProblems_set = moFEMProblems.get<Problem_mi_tag>();
  mofem_problems_by_name::iterator p_miit = moFEMProblems_set.find(problem_name);
  if(p_miit == moFEMProblems_set.end()) SETERRQ1(PETSC_COMM_SELF,1,"problem is not in database %s",problem_name.c_str());

  ierr = problem_basic_method_postProcess(&*p_miit,method,verb); CHKERRQ(ierr);

  PetscFunctionReturn(0);
}
PetscErrorCode Core::loop_finite_elements(const string &problem_name,const string &fe_name,FEMethod &method,int verb) {
  PetscFunctionBegin;
  if(verb==-1) verb = verbose;

  ierr = loop_finite_elements(problem_name,fe_name,method,rAnk,rAnk,verb); CHKERRQ(ierr);

  PetscFunctionReturn(0);
}
PetscErrorCode Core::loop_finite_elements(
  const MoFEMProblem *problem_ptr,const string &fe_name,FEMethod &method,int lower_rank,int upper_rank,int verb) {
  PetscFunctionBegin;
  if(verb==-1) verb = verbose;
  // finite element
  typedef NumeredMoFEMFiniteElement_multiIndex::index<Composite_Name_And_Part_mi_tag>::type FEs_by_composite;

  method.feName = fe_name;
  SET_BASIC_METHOD(method,&*problem_ptr)
  PetscLogEventBegin(USER_EVENT_preProcess,0,0,0,0);
  ierr = method.preProcess(); CHKERRQ(ierr);
  PetscLogEventEnd(USER_EVENT_preProcess,0,0,0,0);

  FEs_by_composite &numeredFiniteElements =
    (const_cast<NumeredMoFEMFiniteElement_multiIndex&>(problem_ptr->numeredFiniteElements)).get<Composite_Name_And_Part_mi_tag>();
  FEs_by_composite::iterator miit = numeredFiniteElements.lower_bound(boost::make_tuple(fe_name,lower_rank));
  FEs_by_composite::iterator hi_miit = numeredFiniteElements.upper_bound(boost::make_tuple(fe_name,upper_rank));

  method.loopSize = distance(miit,hi_miit);
  for(int nn = 0;miit!=hi_miit;miit++,nn++) {

    method.nInTheLoop = nn;
    method.fePtr = &*miit;
    method.dataPtr = const_cast<FEDofMoFEMEntity_multiIndex*>(&(miit->fe_ptr->data_dofs));
    method.rowPtr = const_cast<FENumeredDofMoFEMEntity_multiIndex*>(&(miit->rows_dofs));
    method.colPtr = const_cast<FENumeredDofMoFEMEntity_multiIndex*>(&(miit->cols_dofs));

    try {
      PetscLogEventBegin(USER_EVENT_operator,0,0,0,0);
      ierr = method(); CHKERRQ(ierr);
      PetscLogEventEnd(USER_EVENT_operator,0,0,0,0);
    } catch (const std::exception& ex) {
      ostringstream ss;
      ss << "throw in method: " << ex.what() << " at line " << __LINE__ << " in file " << __FILE__ << endl;
      SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,ss.str().c_str());
    }

  }

  PetscLogEventBegin(USER_EVENT_postProcess,0,0,0,0);
  ierr = method.postProcess(); CHKERRQ(ierr);
  PetscLogEventEnd(USER_EVENT_postProcess,0,0,0,0);

  PetscFunctionReturn(0);
}
PetscErrorCode Core::loop_finite_elements(
  const string &problem_name,const string &fe_name,FEMethod &method,int lower_rank,int upper_rank,int verb) {
  PetscFunctionBegin;
  if(verb==-1) verb = verbose;
  typedef MoFEMProblem_multiIndex::index<Problem_mi_tag>::type mofem_problems_by_name;
  // find p_miit
  mofem_problems_by_name &moFEMProblems_set = moFEMProblems.get<Problem_mi_tag>();
  mofem_problems_by_name::iterator p_miit = moFEMProblems_set.find(problem_name);
  if(p_miit == moFEMProblems_set.end()) SETERRQ1(PETSC_COMM_SELF,1,"problem is not in database %s",problem_name.c_str());

  ierr = loop_finite_elements(&*p_miit,fe_name,method,lower_rank,upper_rank,verb); CHKERRQ(ierr);

  PetscFunctionReturn(0);
}
PetscErrorCode Core::loop_dofs(
  const MoFEMProblem *problem_ptr,const string &field_name,RowColData rc,EntMethod &method,int lower_rank,int upper_rank,int verb
) {
  PetscFunctionBegin;
  SET_BASIC_METHOD(method,&*problem_ptr);
  typedef NumeredDofMoFEMEntity_multiIndex::index<Composite_Name_And_Part_mi_tag>::type numerd_dofs;
  numerd_dofs *dofs;
  switch (rc) {
    case ROW:
      dofs = const_cast<numerd_dofs*>(&problem_ptr->numered_dofs_rows.get<Composite_Name_And_Part_mi_tag>());
      break;
    case COL:
      dofs = const_cast<numerd_dofs*>(&problem_ptr->numered_dofs_cols.get<Composite_Name_And_Part_mi_tag>());
      break;
    default:
     SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"not implemented");
  }
  numerd_dofs::iterator miit = dofs->lower_bound(boost::make_tuple(field_name,lower_rank));
  numerd_dofs::iterator hi_miit = dofs->upper_bound(boost::make_tuple(field_name,upper_rank));
  ierr = method.preProcess(); CHKERRQ(ierr);
  for(;miit!=hi_miit;miit++) {
    method.dofPtr = miit->get_DofMoFEMEntity_ptr();
    method.dofNumeredPtr = &*miit;
    try {
      ierr = method(); CHKERRQ(ierr);
    } catch (const std::exception& ex) {
      ostringstream ss;
      ss << "throw in method: " << ex.what() << endl;
      SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,ss.str().c_str());
    }
  }
  ierr = method.postProcess(); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
PetscErrorCode Core::loop_dofs(
  const string &problem_name,const string &field_name,RowColData rc,EntMethod &method,int lower_rank,int upper_rank,int verb
) {
  PetscFunctionBegin;
  if(verb==-1) verb = verbose;
  typedef MoFEMProblem_multiIndex::index<Problem_mi_tag>::type mofem_problems_by_name;
  // find p_miit
  mofem_problems_by_name &moFEMProblems_set = moFEMProblems.get<Problem_mi_tag>();
  mofem_problems_by_name::iterator p_miit = moFEMProblems_set.find(problem_name);
  if(p_miit == moFEMProblems_set.end()) SETERRQ1(PETSC_COMM_SELF,1,"problem not in database %s",problem_name.c_str());
  ierr = loop_dofs(&*p_miit,field_name,rc,method,lower_rank,upper_rank,verb); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
PetscErrorCode Core::loop_dofs(const string &problem_name,const string &field_name,RowColData rc,EntMethod &method,int verb) {
  PetscFunctionBegin;
  if(verb==-1) verb = verbose;
  ierr = loop_dofs(problem_name,field_name,rc,method,0,sIze,verb); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
PetscErrorCode Core::loop_dofs(const string &field_name,EntMethod &method,int verb) {
  PetscFunctionBegin;
  if(verb==-1) verb = verbose;
  SET_BASIC_METHOD(method,NULL);
  DofMoFEMEntity_multiIndex::index<FieldName_mi_tag>::type::iterator miit,hi_miit;
  miit = dofsMoabField.get<FieldName_mi_tag>().lower_bound(field_name);
  hi_miit = dofsMoabField.get<FieldName_mi_tag>().upper_bound(field_name);
  method.loopSize = distance(miit,hi_miit);
  ierr = method.preProcess(); CHKERRQ(ierr);
  for(int nn = 0; miit!=hi_miit; miit++, nn++) {
    method.nInTheLoop = nn;
    method.dofPtr = &*miit;
    method.dofNumeredPtr = NULL;
    try {
      ierr = method(); CHKERRQ(ierr);
    } catch (const std::exception& ex) {
      ostringstream ss;
      ss << "throw in method: " << ex.what() << endl;
      SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,ss.str().c_str());
    }
  }
  ierr = method.postProcess(); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
PetscErrorCode Core::get_ref_ents(const RefMoFEMEntity_multiIndex **refined_entities_ptr) {
  PetscFunctionBegin;
  *refined_entities_ptr = &refinedEntities;
  PetscFunctionReturn(0);
}
PetscErrorCode Core::get_ref_finite_elements(const RefMoFEMElement_multiIndex **refined_finite_elements_ptr) {
  PetscFunctionBegin;
  *refined_finite_elements_ptr = &refinedFiniteElements;
  PetscFunctionReturn(0);
}

PetscErrorCode Core::get_problem(const string &problem_name,const MoFEMProblem **problem_ptr) {
  PetscFunctionBegin;
  typedef MoFEMProblem_multiIndex::index<Problem_mi_tag>::type mofem_problems_by_name;
  mofem_problems_by_name &moFEMProblems_set = moFEMProblems.get<Problem_mi_tag>();
  mofem_problems_by_name::iterator p_miit = moFEMProblems_set.find(problem_name);
  if(p_miit == moFEMProblems_set.end()) SETERRQ1(PETSC_COMM_SELF,1,"problem < %s > not found, (top tip: check spelling)",problem_name.c_str());
  *problem_ptr = &*p_miit;
  PetscFunctionReturn(0);
}
PetscErrorCode Core::get_dofs(const DofMoFEMEntity_multiIndex **dofs_ptr) {
  PetscFunctionBegin;
  *dofs_ptr = &dofsMoabField;
  PetscFunctionReturn(0);
}
MoFEMEntity_multiIndex::index<FieldName_mi_tag>::type::iterator Core::get_ent_moabfield_by_name_begin(const string &field_name) {
  return entsMoabField.get<FieldName_mi_tag>().lower_bound(field_name);
}
MoFEMEntity_multiIndex::index<FieldName_mi_tag>::type::iterator Core::get_ent_moabfield_by_name_end(const string &field_name) {
  return entsMoabField.get<FieldName_mi_tag>().upper_bound(field_name);
}
DofMoFEMEntity_multiIndex::index<FieldName_mi_tag>::type::iterator Core::get_dofs_by_name_begin(const string &field_name) const {
  return dofsMoabField.get<FieldName_mi_tag>().lower_bound(field_name);
}
DofMoFEMEntity_multiIndex::index<FieldName_mi_tag>::type::iterator Core::get_dofs_by_name_end(const string &field_name) const {
  return dofsMoabField.get<FieldName_mi_tag>().upper_bound(field_name);
}
DofMoFEMEntity_multiIndex::index<Composite_Name_And_Ent_mi_tag>::type::iterator Core::get_dofs_by_name_and_ent_begin(const string &field_name,const EntityHandle ent) {
  return dofsMoabField.get<Composite_Name_And_Ent_mi_tag>().lower_bound(boost::make_tuple(field_name,ent));
}
DofMoFEMEntity_multiIndex::index<Composite_Name_And_Ent_mi_tag>::type::iterator Core::get_dofs_by_name_and_ent_end(const string &field_name,const EntityHandle ent) {
  return dofsMoabField.get<Composite_Name_And_Ent_mi_tag>().upper_bound(boost::make_tuple(field_name,ent));
}
DofMoFEMEntity_multiIndex::index<Composite_Name_And_Type_mi_tag>::type::iterator Core::get_dofs_by_name_and_type_begin(const string &field_name,const EntityType type) {
  return dofsMoabField.get<Composite_Name_And_Type_mi_tag>().lower_bound(boost::make_tuple(field_name,type));
}
DofMoFEMEntity_multiIndex::index<Composite_Name_And_Type_mi_tag>::type::iterator Core::get_dofs_by_name_and_type_end(const string &field_name,const EntityType type) {
  return dofsMoabField.get<Composite_Name_And_Type_mi_tag>().upper_bound(boost::make_tuple(field_name,type));
}
PetscErrorCode Core::get_finite_elements(const MoFEMFiniteElement_multiIndex **finiteElements_ptr) {
  PetscFunctionBegin;
  *finiteElements_ptr = &finiteElements;
  PetscFunctionReturn(0);
}
EntMoFEMFiniteElement_multiIndex::index<FiniteElement_name_mi_tag>::type::iterator Core::get_fe_by_name_begin(const string &fe_name) {
  return finiteElementsMoFEMEnts.get<FiniteElement_name_mi_tag>().lower_bound(fe_name);
}
EntMoFEMFiniteElement_multiIndex::index<FiniteElement_name_mi_tag>::type::iterator Core::get_fe_by_name_end(const string &fe_name) {
  return finiteElementsMoFEMEnts.get<FiniteElement_name_mi_tag>().upper_bound(fe_name);
}
PetscErrorCode Core::check_number_of_ents_in_ents_field(const string& name) {
  PetscFunctionBegin;
  MoFEMField_multiIndex::index<FieldName_mi_tag>::type::iterator it = moabFields.get<FieldName_mi_tag>().find(name);
  if(it == moabFields.get<FieldName_mi_tag>().end()) {
    SETERRQ1(PETSC_COMM_SELF,1,"field not found < %s >",name.c_str());
  }
  EntityHandle meshset = it->get_meshset();
  int num_entities;
  rval = moab.get_number_entities_by_handle(meshset,num_entities); CHKERR_PETSC(rval);
  if(entsMoabField.get<FieldName_mi_tag>().count(it->get_name())
    != (unsigned int)num_entities) {
    SETERRQ1(PETSC_COMM_SELF,1,"not equal number of entities in meshset and field multiindex < %s >",name.c_str());
  }
  PetscFunctionReturn(0);
}
PetscErrorCode Core::check_number_of_ents_in_ents_field() {
  PetscFunctionBegin;
  MoFEMField_multiIndex::index<FieldName_mi_tag>::type::iterator it = moabFields.get<FieldName_mi_tag>().begin();
  for(;it!=moabFields.get<FieldName_mi_tag>().end();it++) {
    if(it->get_space() == NOFIELD) continue; //FIXME: should be treated proprly, not test is just skiped for this NOFIELD space
    EntityHandle meshset = it->get_meshset();
    int num_entities;
    rval = moab.get_number_entities_by_handle(meshset,num_entities); CHKERR_PETSC(rval);
    if(entsMoabField.get<FieldName_mi_tag>().count(it->get_name()) != (unsigned int)num_entities) {
      SETERRQ1(PETSC_COMM_SELF,1,"not equal number of entities in meshset and field multiindex < %s >",it->get_name().c_str());
    }
  }
  PetscFunctionReturn(0);
}
PetscErrorCode Core::check_number_of_ents_in_ents_finite_element(const string& name) {
  PetscFunctionBegin;
  MoFEMFiniteElement_multiIndex::index<FiniteElement_name_mi_tag>::type::iterator it;
  it = finiteElements.get<FiniteElement_name_mi_tag>().find(name);
  if(it == finiteElements.get<FiniteElement_name_mi_tag>().end()) {
    SETERRQ1(PETSC_COMM_SELF,1,"finite element not found < %s >",name.c_str());
  }
  EntityHandle meshset = it->get_meshset();
  int num_entities;
  rval = moab.get_number_entities_by_handle(meshset,num_entities); CHKERR_PETSC(rval);
  if(finiteElementsMoFEMEnts.get<FiniteElement_name_mi_tag>().count(it->get_name().c_str())
    != (unsigned int)num_entities) {
    SETERRQ1(PETSC_COMM_SELF,1,"not equal number of entities in meshset and finite elements multiindex < %s >",it->get_name().c_str());
  }
  PetscFunctionReturn(0);
}
PetscErrorCode Core::check_number_of_ents_in_ents_finite_element() {
  PetscFunctionBegin;
  MoFEMFiniteElement_multiIndex::index<FiniteElement_name_mi_tag>::type::iterator it;
  it = finiteElements.get<FiniteElement_name_mi_tag>().begin();
  for(;it!=finiteElements.get<FiniteElement_name_mi_tag>().end();it++) {
    EntityHandle meshset = it->get_meshset();
    int num_entities;
    rval = moab.get_number_entities_by_handle(meshset,num_entities); CHKERR_PETSC(rval);
    if(finiteElementsMoFEMEnts.get<FiniteElement_name_mi_tag>().count(it->get_name().c_str())
      != (unsigned int)num_entities) {
      SETERRQ1(PETSC_COMM_SELF,1,"not equal number of entities in meshset and finite elements multiindex < %s >",it->get_name().c_str());
    }
  }
  PetscFunctionReturn(0);
}

//clear,remove and delete

PetscErrorCode Core::clear_dofs_fields(const BitRefLevel &bit,const BitRefLevel &mask,int verb) {
  PetscFunctionBegin;
  if(verb==-1) verb = verbose;
  {
    DofMoFEMEntity_multiIndex::iterator dit;
    dit = dofsMoabField.begin();
    for(;dit!=dofsMoabField.end();) {
      BitRefLevel bit2 = dit->get_BitRefLevel();
      if(dit->get_ent_type()==MBENTITYSET) {
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
      MoFEMEntityEntMoFEMFiniteElementAdjacencyMap_multiIndex::index<Unique_mi_tag>::type::iterator ait,hi_ait;
      ait = entFEAdjacencies.get<Unique_mi_tag>().lower_bound(dit->get_MoFEMEntity_ptr()->get_global_unique_id());
      hi_ait = entFEAdjacencies.get<Unique_mi_tag>().upper_bound(dit->get_MoFEMEntity_ptr()->get_global_unique_id());
      for(;ait!=hi_ait;ait++) {
        EntMoFEMFiniteElement *EntMoFEMFiniteElement_ptr;
        EntMoFEMFiniteElement_ptr = const_cast<EntMoFEMFiniteElement *>(ait->EntMoFEMFiniteElement_ptr);
        EntMoFEMFiniteElement_ptr->row_dof_view.erase(dit->get_global_unique_id());
        EntMoFEMFiniteElement_ptr->col_dof_view.erase(dit->get_global_unique_id());
        EntMoFEMFiniteElement_ptr->data_dof_view.erase(dit->get_global_unique_id());
        EntMoFEMFiniteElement_ptr->data_dofs.get<Unique_mi_tag>().erase(dit->get_global_unique_id());
      }
      dit = dofsMoabField.erase(dit);
    }
  }
  PetscFunctionReturn(0);
}
PetscErrorCode Core::clear_dofs_fields(const string &name,const Range ents,int verb) {
  PetscFunctionBegin;
  if(verb==-1) verb = verbose;
  Range::iterator eit = ents.begin();
  for(;eit!=ents.end();eit++) {
    DofMoFEMEntity_multiIndex::index<Composite_Name_And_Ent_mi_tag>::type::iterator dit,hi_dit;
    dit = dofsMoabField.get<Composite_Name_And_Ent_mi_tag>().lower_bound(boost::make_tuple(name,*eit));
    hi_dit = dofsMoabField.get<Composite_Name_And_Ent_mi_tag>().upper_bound(boost::make_tuple(name,*eit));
    for(;dit!=hi_dit;) {
      MoFEMEntityEntMoFEMFiniteElementAdjacencyMap_multiIndex::index<Unique_mi_tag>::type::iterator ait,hi_ait;
      ait = entFEAdjacencies.get<Unique_mi_tag>().lower_bound(dit->get_MoFEMEntity_ptr()->get_global_unique_id());
      hi_ait = entFEAdjacencies.get<Unique_mi_tag>().upper_bound(dit->get_MoFEMEntity_ptr()->get_global_unique_id());
      for(;ait!=hi_ait;ait++) {
  	EntMoFEMFiniteElement *EntMoFEMFiniteElement_ptr;
	EntMoFEMFiniteElement_ptr = const_cast<EntMoFEMFiniteElement *>(ait->EntMoFEMFiniteElement_ptr);
	EntMoFEMFiniteElement_ptr->row_dof_view.erase(dit->get_global_unique_id());
	EntMoFEMFiniteElement_ptr->col_dof_view.erase(dit->get_global_unique_id());
	EntMoFEMFiniteElement_ptr->data_dof_view.erase(dit->get_global_unique_id());
	EntMoFEMFiniteElement_ptr->data_dofs.get<Unique_mi_tag>().erase(dit->get_global_unique_id());
      }
      dit = dofsMoabField.get<Composite_Name_And_Ent_mi_tag>().erase(dit);
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
  eit = entsMoabField.begin();
  for(;eit!=entsMoabField.end();) {
    if(eit->get_ent_type()==MBENTITYSET) {
      eit++;
      continue;
    }
    BitRefLevel bit2 = eit->get_BitRefLevel();
    if((bit2&mask)!=bit2) {
      eit++;
      continue;
    }
    if((bit2&bit).none()) {
      eit++;
      continue;
    }
    EntityHandle ent = eit->get_ent();
    rval = moab.tag_delete_data(eit->field_ptr->th_AppOrder,&ent,1); CHKERR_PETSC(rval);
    if(eit->tag_FieldData_size>0) {
      rval = moab.tag_delete_data(eit->field_ptr->th_FieldData,&ent,1); CHKERR_PETSC(rval);
    }
    eit = entsMoabField.erase(eit);
  }
  PetscFunctionReturn(0);
}
PetscErrorCode Core::clear_ents_fields(const string &name,const Range ents,int verb) {
  PetscFunctionBegin;
  if(verb==-1) verb = verbose;
  ierr = clear_dofs_fields(name,ents,verb); CHKERRQ(ierr);
  ierr = clear_adjacencies_entities(name,ents,verb); CHKERRQ(ierr);
  Range::iterator eit = ents.begin();
  for(;eit!=ents.end();eit++) {
    MoFEMEntity_multiIndex::index<Composite_Name_And_Ent_mi_tag>::type::iterator dit,hi_dit;
    dit = entsMoabField.get<Composite_Name_And_Ent_mi_tag>().lower_bound(boost::make_tuple(name,*eit));
    hi_dit = entsMoabField.get<Composite_Name_And_Ent_mi_tag>().upper_bound(boost::make_tuple(name,*eit));
    for(;dit!=hi_dit;) {
      dit = entsMoabField.get<Composite_Name_And_Ent_mi_tag>().erase(dit);
    }
  }
  PetscFunctionReturn(0);
}
PetscErrorCode Core::clear_finite_elements(const BitRefLevel &bit,const BitRefLevel &mask,int verb) {
  PetscFunctionBegin;
  if(verb==-1) verb = verbose;
  ierr = clear_adjacencies_finite_elements(bit,mask,verb); CHKERRQ(ierr);
  EntMoFEMFiniteElement_multiIndex::iterator fe_it = finiteElementsMoFEMEnts.begin();
  for(;fe_it!=finiteElementsMoFEMEnts.end();) {
    BitRefLevel bit2 = fe_it->get_BitRefLevel();
    if(fe_it->get_ent_type()==MBENTITYSET) {
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
    fe_it = finiteElementsMoFEMEnts.erase(fe_it);
  }
  PetscFunctionReturn(0);
}
PetscErrorCode Core::clear_finite_elements(const string &name,const Range &ents,int verb) {
  PetscFunctionBegin;
  if(verb==-1) verb = verbose;
  ierr = clear_adjacencies_finite_elements(name,ents,verb); CHKERRQ(ierr);
  Range::iterator eit = ents.begin();
  for(;eit!=ents.end();eit++) {
    EntMoFEMFiniteElement_multiIndex::index<Composite_Name_And_Ent_mi_tag>::type::iterator fit,hi_fit;
    fit = finiteElementsMoFEMEnts.get<Composite_Name_And_Ent_mi_tag>().lower_bound(boost::make_tuple(name,*eit));
    hi_fit = finiteElementsMoFEMEnts.get<Composite_Name_And_Ent_mi_tag>().upper_bound(boost::make_tuple(name,*eit));
    for(;fit!=hi_fit;) {
      fit = finiteElementsMoFEMEnts.get<Composite_Name_And_Ent_mi_tag>().erase(fit);
    }
  }
  PetscFunctionReturn(0);
}
PetscErrorCode Core::clear_adjacencies_finite_elements(const BitRefLevel &bit,const BitRefLevel &mask,int verb) {
  PetscFunctionBegin;
  MoFEMEntityEntMoFEMFiniteElementAdjacencyMap_multiIndex::iterator ait;
  ait = entFEAdjacencies.begin();
  for(;ait!=entFEAdjacencies.end();) {
    BitRefLevel bit2 = ait->EntMoFEMFiniteElement_ptr->get_BitRefLevel();
    if(ait->EntMoFEMFiniteElement_ptr->get_ent_type()==MBENTITYSET) {
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
PetscErrorCode Core::clear_adjacencies_finite_elements(const string &name,const Range &ents,int verb) {
  PetscFunctionBegin;
  if(verb==-1) verb = verbose;
  Range::iterator eit = ents.begin();
  for(;eit!=ents.end();eit++) {
    MoFEMEntityEntMoFEMFiniteElementAdjacencyMap_multiIndex::index<FEEnt_mi_tag>::type::iterator ait,hi_ait;
    ait = entFEAdjacencies.get<FEEnt_mi_tag>().lower_bound(*eit);
    hi_ait = entFEAdjacencies.get<FEEnt_mi_tag>().upper_bound(*eit);
    for(;ait!=hi_ait;) {
      if(ait->EntMoFEMFiniteElement_ptr->get_name() == name) {
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
  MoFEMEntityEntMoFEMFiniteElementAdjacencyMap_multiIndex::iterator ait;
  ait = entFEAdjacencies.begin();
  for(;ait!=entFEAdjacencies.end();) {
    BitRefLevel bit2 = ait->MoFEMEntity_ptr->get_BitRefLevel();
    if(ait->MoFEMEntity_ptr->get_ent_type()==MBENTITYSET) {
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
PetscErrorCode Core::clear_adjacencies_entities(const string &name,const Range &ents,int verb) {
  PetscFunctionBegin;
  if(verb==-1) verb = verbose;
  Range::iterator eit = ents.begin();
  for(;eit!=ents.end();eit++) {
    MoFEMEntityEntMoFEMFiniteElementAdjacencyMap_multiIndex::index<Ent_mi_tag>::type::iterator ait,hi_ait;
    ait = entFEAdjacencies.get<Ent_mi_tag>().lower_bound(*eit);
    hi_ait = entFEAdjacencies.get<Ent_mi_tag>().upper_bound(*eit);
    for(;ait!=hi_ait;) {
      if(ait->MoFEMEntity_ptr->get_name() == name) {
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
  MoFEMField_multiIndex::iterator f_it = moabFields.begin();
  for(;f_it!=moabFields.end();f_it++) {
    EntityHandle meshset = f_it->get_meshset();
    Range ents_to_remove;
    rval = moab.get_entities_by_handle(
      meshset,ents_to_remove,false); CHKERR_PETSC(rval);
    Range::iterator eit = ents_to_remove.begin();
    for(;eit!=ents_to_remove.end();) {
      if(moab.type_from_handle(*eit)==MBENTITYSET) {
        eit = ents_to_remove.erase(eit);
        continue;
      }
      BitRefLevel bit2;
      rval = moab.tag_get_data(th_RefBitLevel,&*eit,1,&bit2); CHKERR_PETSC(rval);
      if((bit2&mask)!=bit2) {
        eit = ents_to_remove.erase(eit);
        continue;
      }
      if((bit2&bit).none()) {
        eit = ents_to_remove.erase(eit);
        continue;
      }
      MoFEMEntity_multiIndex::index<Composite_Name_And_Ent_mi_tag>::type::iterator iit;
      iit = entsMoabField.get<Composite_Name_And_Ent_mi_tag>().find(boost::make_tuple(f_it->get_name(),*eit));
      if(iit != entsMoabField.get<Composite_Name_And_Ent_mi_tag>().end()) {
        SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"data inconsistency");
      }
      eit++;
    }
    rval = moab.remove_entities(meshset,ents_to_remove); CHKERR_PETSC(rval);
    if(verb>0) {
      PetscPrintf(comm,
        "number of removed entities = %u from field %s\n",ents_to_remove.size(),f_it->get_name().c_str());
        if(verb>1) {
          int num_entities;
          rval = moab.get_number_entities_by_handle(meshset,num_entities); CHKERR_PETSC(rval);
          PetscPrintf(comm,"\tnumber of entities in database = %u and meshset = %u\n",
          entsMoabField.get<BitFieldId_mi_tag>().count(f_it->get_id()),num_entities);
        }
      }
    }
  PetscFunctionReturn(0);
}
PetscErrorCode Core::remove_ents_from_field(const string& name,const EntityHandle meshset,const EntityType type,int verb) {
  PetscFunctionBegin;
  if(verb==-1) verb = verbose;
  Range ents;
  rval = moab.get_entities_by_type(meshset,type,ents); CHKERR_PETSC(rval);
  ierr = remove_ents_from_field(name,ents,verb); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
PetscErrorCode Core::remove_ents_from_field(const string& name,const Range &ents,int verb) {
  PetscFunctionBegin;
  if(verb==-1) verb = verbose;
  EntityHandle meshset;
  try {
    meshset = get_field_meshset(name);
  } catch (MoFEMException const &e) {
    SETERRQ(PETSC_COMM_SELF,e.errorCode,e.errorMessage);
  }
  rval = moab.remove_entities(meshset,ents); CHKERR_PETSC(rval);
  ierr = clear_ents_fields(name,ents,verb); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
PetscErrorCode Core::remove_ents_from_finite_element_by_bit_ref(const BitRefLevel &bit,const BitRefLevel &mask,int verb) {
  PetscFunctionBegin;
  if(verb==-1) verb = verbose;
  ierr = clear_finite_elements(bit,mask,verb); CHKERRQ(ierr);
  MoFEMFiniteElement_multiIndex::iterator fe_it = finiteElements.begin();
  for(;fe_it!=finiteElements.end();fe_it++) {
    EntityHandle meshset = fe_it->get_meshset();
    Range ents_to_remove;
    rval = moab.get_entities_by_handle(
      meshset,ents_to_remove,false
    ); CHKERR_PETSC(rval);
    Range::iterator eit = ents_to_remove.begin();
    for(;eit!=ents_to_remove.end();) {
      if(moab.type_from_handle(*eit)==MBENTITYSET) {
        eit = ents_to_remove.erase(eit);
        continue;
      }
      BitRefLevel bit2;
      rval = moab.tag_get_data(th_RefBitLevel,&*eit,1,&bit2); CHKERR_PETSC(rval);
      if((bit2&mask)!=bit2) {
        eit = ents_to_remove.erase(eit);
        continue;
      }
      if((bit2&bit).none()) {
        eit = ents_to_remove.erase(eit);
        continue;
      }
      EntMoFEMFiniteElement_multiIndex::index<Composite_Name_And_Ent_mi_tag>::type::iterator iit;
      iit = finiteElementsMoFEMEnts.get<Composite_Name_And_Ent_mi_tag>().find(
        boost::make_tuple(fe_it->get_name(),*eit)
      );
      if(iit != finiteElementsMoFEMEnts.get<Composite_Name_And_Ent_mi_tag>().end()) {
        SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"data inconsistency");
      }
      eit++;
    }
    rval = moab.remove_entities(meshset,ents_to_remove); CHKERR_PETSC(rval);
    if(verb>0) {
      PetscPrintf(comm,
        "number of removed entities = %u from finite element %s\n",
        ents_to_remove.size(),fe_it->get_name().c_str()
      );
    }
  }
  PetscFunctionReturn(0);
}
PetscErrorCode Core::remove_ents_from_finite_element(const string &name,const EntityHandle meshset,const EntityType type,int verb) {
  PetscFunctionBegin;
  if(verb==-1) verb = verbose;
  Range ents;
  rval = moab.get_entities_by_type(meshset,type,ents,false); CHKERR_PETSC(rval);
  ierr = remove_ents_from_finite_element(name,ents,verb); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
PetscErrorCode Core::remove_ents_from_finite_element(const string &name,const Range &ents,int verb) {
  PetscFunctionBegin;
  if(verb==-1) verb = verbose;
  ierr = clear_finite_elements(name,ents,verb); CHKERRQ(ierr);
  const BitFEId id = get_BitFEId(name);
  const EntityHandle idm = get_finite_element_meshset(id);
  rval = moab.remove_entities(idm,ents); CHKERR_PETSC(rval);
  PetscFunctionReturn(0);
}
PetscErrorCode Core::remove_ents_by_bit_ref(const BitRefLevel &bit,const BitRefLevel &mask,int verb) {
  PetscFunctionBegin;
  if(verb==-1) verb = verbose;
  ierr = delete_finite_elements_by_bit_ref(bit,mask,verb); CHKERRQ(ierr);
  ierr = remove_ents_from_field_by_bit_ref(bit,mask,verb); CHKERRQ(ierr);
  RefMoFEMEntity_multiIndex::iterator ent_it = refinedEntities.begin();
  for(;ent_it!=refinedEntities.end();) {
    BitRefLevel bit2 = ent_it->get_BitRefLevel();
    if(ent_it->get_ent_type()==MBENTITYSET) {
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
PetscErrorCode Core::delete_ents_by_bit_ref(const BitRefLevel &bit,const BitRefLevel &mask,const bool remove_parent,int verb) {
  PetscFunctionBegin;
  Range ents_to_delete;
  rval = moab.get_entities_by_handle(0,ents_to_delete,false); CHKERR_PETSC(rval);
  {
    Range::iterator eit = ents_to_delete.begin();
    for(;eit!=ents_to_delete.end();) {
      if(moab.type_from_handle(*eit)==MBENTITYSET) {
        eit = ents_to_delete.erase(eit);
        continue;
      }
      BitRefLevel bit2;
      rval = moab.tag_get_data(th_RefBitLevel,&*eit,1,&bit2); CHKERR_PETSC(rval);
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
      RefMoFEMEntity_multiIndex::index<Ent_Ent_mi_tag>::type::iterator pit,hi_pit;
      pit = refinedEntities.get<Ent_Ent_mi_tag>().lower_bound(*eit);
      hi_pit = refinedEntities.get<Ent_Ent_mi_tag>().upper_bound(*eit);
      for(;pit!=hi_pit;pit++) {
        EntityHandle ent = pit->get_ref_ent();
        if(ents_to_delete.find(ent) != ents_to_delete.end()) {
          continue;
        }
        /*if(rAnk==0) {
        EntityHandle out_meshset;
        rval = moab.create_meshset(MESHSET_SET,out_meshset); CHKERR_PETSC(rval);
        rval = moab.add_entities(out_meshset,&ent,1); CHKERR_PETSC(rval);
        rval = moab.add_entities(out_meshset,&*eit,1); CHKERR_PETSC(rval);
        rval = moab.write_file("error.vtk","VTK","",&out_meshset,1); CHKERR_PETSC(rval);
        rval = moab.delete_entities(&out_meshset,1); CHKERR_PETSC(rval);
        }
        ostringstream ss;
        ss << "child:\n" << *pit << endl;
        ss << "parent:\n" << RefMoFEMEntity(moab,*eit) << endl;
        SETERRQ1(PETSC_COMM_SELF,1,
        "entity can not be removed, it is parent for some other entity\n%s",ss.str().c_str());*/
        bool success = refinedEntities.modify(
          refinedEntities.project<0>(pit),RefMoFEMEntity_change_remove_parent(moab)
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
      rval = moab.remove_entities(cubit_meshset,ents_to_delete); CHKERR_PETSC(rval);
      Range meshsets;
      rval = moab.get_entities_by_type(cubit_meshset,MBENTITYSET,meshsets);  CHKERR_PETSC(rval);
      for(Range::iterator mit = meshsets.begin();mit!=meshsets.end();mit++) {
        rval = moab.remove_entities(*mit,ents_to_delete); CHKERR_PETSC(rval);
      }
    }
  }
  ierr = remove_ents_by_bit_ref(bit,mask,verb); CHKERRQ(ierr);
  if(verb>0) {
    PetscPrintf(comm,"number of deleted entities = %u\n",ents_to_delete.size());
  }
  if(verb>2) {
    EntityHandle out_meshset;
    rval = moab.create_meshset(MESHSET_SET,out_meshset); CHKERR_PETSC(rval);
    rval = moab.add_entities(out_meshset,ents_to_delete.subset_by_type(MBTET)); CHKERR_PETSC(rval);
    rval = moab.write_file("debug_ents_to_delete.vtk","VTK","",&out_meshset,1); CHKERR_PETSC(rval);
    rval = moab.delete_entities(&out_meshset,1); CHKERR_PETSC(rval);
  }
  //delete entities form moab
  for(int dd = 3;dd>=0;dd--) {
    rval = moab.delete_entities(ents_to_delete.subset_by_dimension(dd)); //CHKERR_PETSC(rval);
  }
  PetscFunctionReturn(0);
}
PetscErrorCode Core::delete_finite_elements_by_bit_ref(const BitRefLevel &bit,const BitRefLevel &mask,int verb) {
  PetscFunctionBegin;
  if(verb==-1) verb = verbose;
  ierr = remove_ents_from_finite_element_by_bit_ref(bit,mask,verb); CHKERRQ(ierr);
  RefMoFEMElement_multiIndex::iterator fe_it = refinedFiniteElements.begin();
  for(;fe_it!=refinedFiniteElements.end();) {
    BitRefLevel bit2 = fe_it->get_BitRefLevel();
    if(fe_it->get_ent_type()==MBENTITYSET) {
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


}
