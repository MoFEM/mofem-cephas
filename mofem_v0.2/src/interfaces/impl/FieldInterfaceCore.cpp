/** \file FieldInterfaceCore.cpp
 * \brief Myltindex containes, data structures and other low-level functions 
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

#include <moab/ParallelComm.hpp>

#include <petscsys.h>
#include <petscvec.h> 
#include <petscmat.h> 
#include <petscsnes.h> 
#include <petscts.h> 

#include <definitions.h>
#include <h1_hdiv_hcurl_l2.h>

#include <Common.hpp>
#include <LoopMethods.hpp>

#include <boost/ptr_container/ptr_map.hpp>
#include <Core.hpp>

#include <CoreDataStructures.hpp>

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
    BitFieldId id = get_field_shift();
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
    } catch (const char* msg) {
      SETERRQ(PETSC_COMM_SELF,MOFEM_CHAR_THROW,msg);
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
  ierr = clear_map(); CHKERRQ(ierr);
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
      if((base_meshset.CubitBCType&CubitBC_BitSet(NODESET|SIDESET|BLOCKSET)).any()) {
	pair<CubitMeshSet_multiIndex::iterator,bool> p = cubit_meshsets.insert(base_meshset);
	if(!p.second) SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INSONSISTENCY,"meshset not inserted");
	if(verb > 0) {
	  ostringstream ss;
	  ss << "read cubit" << base_meshset << endl;
	  //PetscSynchronizedPrintf(comm,ss.str().c_str());
	  PetscPrintf(comm,ss.str().c_str());
	}
      }
    } catch (const char* msg) {
      SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INSONSISTENCY,msg);
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
      } catch (const char* msg) {
	SETERRQ(PETSC_COMM_SELF,MOFEM_CHAR_THROW,msg);
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
	    ierr = add_prism_to_mofem_database(*eit,verb); CHKERRQ(ierr);
	    p_MoFEMFiniteElement = refinedFiniteElements.insert(ptrWrapperRefMoFEMElement(new RefMoFEMElement_PRISM(moab,&*p_ref_ent.first)));
	    break;
	  case MBENTITYSET:
  	    p_MoFEMFiniteElement = refinedFiniteElements.insert(ptrWrapperRefMoFEMElement(new RefMoFEMElement_MESHSET(moab,&*p_ref_ent.first)));
	    break;
	  default:
	    SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INSONSISTENCY,"Only finite elements of type MBTET, MBPRISM and MBENTITYSET are implemented");
	}
	if(p_MoFEMFiniteElement.second) {
	  //PetscPrintf(comm,"Warrning: this entity should be already in refined finite elements database");
	  //SETERRQ(PETSC_COMM_SELF,1,"data inconsistency, this entity should be already in refined finite elements database");
	}
	} catch (const char* msg) {
	  SETERRQ(PETSC_COMM_SELF,MOFEM_CHAR_THROW,msg);
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
	pair<Series_multiIndex::iterator,bool> p = series.insert(MoFEMSeries(moab,*mit));
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
  for(Series_multiIndex::iterator sit = series.begin();sit!=series.end();sit++) {
    int nb_steps;
    ierr = sit->get_nb_steps(moab,nb_steps); CHKERRQ(ierr);
    int ss = 0;
    for(;ss<nb_steps;ss++) {
      pair<SeriesStep_multiIndex::iterator,bool> p = series_steps.insert(MoFEMSeriesStep(moab,&*sit,ss));
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
  } catch (const char* msg) {
    SETERRQ(PETSC_COMM_SELF,MOFEM_CHAR_THROW,msg);
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
      SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INSONSISTENCY,"add_ents_to_field_by_EDGEs this field not work for EDGEs");
  }
  PetscFunctionReturn(0);
}
PetscErrorCode Core::add_ents_to_field_by_EDGEs(const EntityHandle meshset,const string& name,int verb) {
  PetscFunctionBegin;
  if(verb==-1) verb = verbose;
  *build_MoFEM = 0;
  try {
    ierr = add_ents_to_field_by_EDGEs(meshset,get_BitFieldId(name),verb);  CHKERRQ(ierr);
  } catch  (const char* msg) {
    SETERRQ(PETSC_COMM_SELF,MOFEM_CHAR_THROW ,msg);
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
  } catch (const char* msg) {
    SETERRQ(PETSC_COMM_SELF,MOFEM_CHAR_THROW ,msg);
  }
  FieldSpace space;
  rval = moab.tag_get_data(th_FieldSpace,&idm,1,&space); CHKERR_PETSC(rval);
  Range nodes,edges;
  switch (space) {
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
      SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INSONSISTENCY,"add_ents_to_field_by_TRIs this field not work for TRIs");
  }
  PetscFunctionReturn(0);
}
PetscErrorCode Core::add_ents_to_field_by_TRIs(const EntityHandle meshset,const string& name,int verb) {
  PetscFunctionBegin;
  if(verb==-1) verb = verbose;
  *build_MoFEM = 0;
  try {
    ierr = add_ents_to_field_by_TRIs(meshset,get_BitFieldId(name),verb);  CHKERRQ(ierr);
  } catch  (const char* msg) {
    SETERRQ(PETSC_COMM_SELF,1,msg);
  }
  PetscFunctionReturn(0);
}
PetscErrorCode Core::add_ents_to_field_by_TRIs(const Range &tris,const string& name,int verb) {
  PetscFunctionBegin;
  if(verb==-1) verb = verbose;
  *build_MoFEM = 0;
  try {
    ierr = add_ents_to_field_by_TRIs(tris,get_BitFieldId(name),verb);  CHKERRQ(ierr);
  } catch  (const char* msg) {
    SETERRQ(PETSC_COMM_SELF,1,msg);
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
  } catch (const char* msg) {
    SETERRQ(PETSC_COMM_SELF,1,msg);
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
      SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INSONSISTENCY,"add_ents_to_field_by_TRIs this field not work for TRIs");
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
  } catch  (const char* msg) {
    SETERRQ(PETSC_COMM_SELF,MOFEM_CHAR_THROW,msg);
  }
  PetscFunctionReturn(0);
}
PetscErrorCode Core::add_ents_to_field_by_VERTICEs(const EntityHandle meshset,const string& name,int verb) {
  PetscFunctionBegin;
  if(verb==-1) verb = verbose;
  *build_MoFEM = 0;
  try {
    ierr = add_ents_to_field_by_VERTICEs(meshset,get_BitFieldId(name),verb);  CHKERRQ(ierr);
  } catch  (const char* msg) {
    SETERRQ(PETSC_COMM_SELF,MOFEM_CHAR_THROW,msg);
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
  } catch (const char* msg) {
    SETERRQ(PETSC_COMM_SELF,MOFEM_CHAR_THROW,msg);
  }
  FieldSpace space;
  rval = moab.tag_get_data(th_FieldSpace,&idm,1,&space); CHKERR_PETSC(rval);
  Range nodes,tris,edges;
  switch (space) {
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
      SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INSONSISTENCY,"add_ents_to_field_by_TETs this field not work for TETs");
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
  } catch  (const char* msg) {
    SETERRQ(PETSC_COMM_SELF,MOFEM_CHAR_THROW,msg);
  }
  PetscFunctionReturn(0);
}
PetscErrorCode Core::add_ents_to_field_by_TETs(const EntityHandle meshset,const string& name,int verb) {
  PetscFunctionBegin;
  if(verb==-1) verb = verbose;
  *build_MoFEM = 0;
  try {
    ierr = add_ents_to_field_by_TETs(meshset,get_BitFieldId(name),verb);  CHKERRQ(ierr);
  } catch  (const char* msg) {
    SETERRQ(PETSC_COMM_SELF,MOFEM_CHAR_THROW ,msg);
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
  if(miit==set_id.end()) SETERRQ(PETSC_COMM_SELF,MOFEM_NOT_FOUND,"no id found"); 
  EntityHandle idm = no_handle;
  try {
   idm = get_field_meshset(id);
  } catch (const char* msg) {
    SETERRQ(PETSC_COMM_SELF,MOFEM_CHAR_THROW,msg);
  }
  //itersection with field meshset
  Range ents_of_id_meshset;
  rval = moab.get_entities_by_handle(idm,ents_of_id_meshset,false); CHKERR_PETSC(rval);
  Range ents_ = intersect(ents,ents_of_id_meshset);
  if(verb>1) {
    PetscSynchronizedPrintf(comm,"nb. of ents for order change in the field %d\n",ents_.size());
  }
  vector<const void*> tag_data_order(ents.size());
  rval = moab.tag_get_by_ptr(miit->th_AppOrder,ents,&tag_data_order[0]); CHKERR_PETSC(rval);
  //ent view by field id (in set all MoabEnts has the same FieldId)
  typedef MoFEMEntity_multiIndex::index<BitFieldId_mi_tag>::type ent_set_by_id;
  ent_set_by_id& set = entsMoabField.get<BitFieldId_mi_tag>();
  ent_set_by_id::iterator miit2 = set.lower_bound(id);
  MoFEMEntity_multiIndex_ent_view ents_id_view;
  if(miit2 != set.end()) {
    ent_set_by_id::iterator hi_miit2 = set.upper_bound(id);
    for(;miit2!=hi_miit2;miit2++) {
      ents_id_view.insert(&*miit2);
    }
  }
  if(verb>1) {
    PetscSynchronizedPrintf(comm,"nb. of ents in the multi index field %d\n",ents_id_view.size());
  }
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
	    SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INSONSISTENCY,"approximation order for H1 sapce and vertex diffrent than 1 makes not sense"); 
	  }
	}
	break;
      case HDIV:
	if(moab.type_from_handle(*eit)==MBVERTEX) {
	  SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INSONSISTENCY,"HDIV space on vertices makes no sense"); 
	} 
	if(moab.type_from_handle(*eit)==MBEDGE) {
	  SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INSONSISTENCY,"HDIV space on edges makes no sense"); 
	} 
	break;
      default:
	break;
    }
    //
    MoFEMEntity_multiIndex_ent_view::iterator miit3 = ents_id_view.find(*eit);
    if(miit3!=ents_id_view.end()) {
      const ApproximationOrder old_ApproximationOrder = (*miit3)->get_max_order();
      if(old_ApproximationOrder==order) continue;
      MoFEMEntity_multiIndex::iterator miit4 = entsMoabField.get<Unique_mi_tag>().find((*miit3)->get_global_unique_id());
      assert(miit4!=entsMoabField.end());
      if(miit4->get_max_order()<order) nb_ents_set_order_up++;
      if(miit4->get_max_order()>order) nb_ents_set_order_down++;
      typedef DofMoFEMEntity_multiIndex::index<Composite_Name_And_Ent_mi_tag>::type dof_set_type;
      dof_set_type& set_set = dofsMoabField.get<Composite_Name_And_Ent_mi_tag>();
      dof_set_type::iterator miit5 = set_set.lower_bound(boost::make_tuple(miit4->get_name_ref(),miit4->get_ent()));
      dof_set_type::iterator hi_miit6 = set_set.upper_bound(boost::make_tuple(miit4->get_name_ref(),miit4->get_ent()));
      for(;miit5!=hi_miit6;miit5++) {
	if(miit5->get_dof_order()<=order) continue;
	bool success = dofsMoabField.modify(dofsMoabField.project<0>(miit5),DofMoFEMEntity_active_change(false));
	if(!success) SETERRQ(PETSC_COMM_SELF,MOFEM_OPERATION_UNSUCCESSFUL,"modification unsuccessful");
      }
      bool success = entsMoabField.modify(entsMoabField.project<0>(miit4),MoFEMEntity_change_order(moab,order));
      if(!success) SETERRQ(PETSC_COMM_SELF,MOFEM_OPERATION_UNSUCCESSFUL,"modification unsuccessful");
    } else {
      *(ApproximationOrder*)tag_data_order[ee] = order;
      RefMoFEMEntity_multiIndex::index<Ent_mi_tag>::type::iterator miit_ref_ent = refinedEntities.get<Ent_mi_tag>().find(*eit);
      if(miit_ref_ent==refinedEntities.get<Ent_mi_tag>().end()) {
	RefMoFEMEntity ref_ent(moab,*eit);
	if(ref_ent.get_BitRefLevel().none()) continue; // not on any mesh and not in database 
	cerr << ref_ent << endl;
	cerr << "bit level " << ref_ent.get_BitRefLevel() << endl;
	SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INSONSISTENCY,"database inconsistency");
      }
      try { 
	MoFEMEntity moabent(moab,&*miit,&*miit_ref_ent);
	//if(moabent.get_order_nb_dofs(moabent.get_max_order())==0) continue; 
	pair<MoFEMEntity_multiIndex::iterator,bool> e_miit = entsMoabField.insert(moabent);
	bool success = entsMoabField.modify(e_miit.first,MoFEMEntity_change_order(moab,order));
	if(!success) SETERRQ(PETSC_COMM_SELF,MOFEM_OPERATION_UNSUCCESSFUL,"modification unsuccessful");
	nb_ents_set_order_new++;
      } catch (const char* msg) {
	SETERRQ(PETSC_COMM_SELF,MOFEM_CHAR_THROW ,msg);
      } catch (const std::exception& ex) {
	ostringstream ss;
	ss << ex.what() << endl;
	SETERRQ(PETSC_COMM_SELF,MOFEM_STD_EXCEPTION_THROW,ss.str().c_str());
      }
    }
  }
  if(verb>1) {
    PetscSynchronizedPrintf(comm,"nb. of ents for which order was increased %d (order %d)\n",nb_ents_set_order_up,order);
    PetscSynchronizedPrintf(comm,"nb. of ents for which order was reduced %d (order %d)\n",nb_ents_set_order_down,order);
    PetscSynchronizedPrintf(comm,"nb. of ents for which order set %d (order %d)\n",nb_ents_set_order_new,order);
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
  } catch (const char* msg) {
    SETERRQ(PETSC_COMM_SELF,MOFEM_CHAR_THROW,msg);
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
  } catch (const char* msg) {
    SETERRQ(PETSC_COMM_SELF,MOFEM_CHAR_THROW,msg);
  } 
  PetscFunctionReturn(0);
}
PetscErrorCode Core::set_field_order(const Range &ents,const string& name,const ApproximationOrder order,int verb) {
  PetscFunctionBegin;
  if(verb==-1) verb = verbose;
  *build_MoFEM = 0;
  try{
    ierr = set_field_order(ents,get_BitFieldId(name),order,verb); CHKERRQ(ierr);
  } catch (const char* msg) {
    SETERRQ(PETSC_COMM_SELF,MOFEM_CHAR_THROW,msg);
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
  } catch (const char* msg) {
    SETERRQ(PETSC_COMM_SELF,MOFEM_CHAR_THROW,msg);
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
  } catch (const char* msg) {
    SETERRQ(PETSC_COMM_SELF,MOFEM_CHAR_THROW,msg);
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
  if(miit_ref_ent==refinedEntities.get<Ent_mi_tag>().end()) SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INSONSISTENCY,"database inconsistency");
  pair<MoFEMEntity_multiIndex::iterator,bool> e_miit;
  try {
    //create database entity
    MoFEMEntity moabent(moab,&*miit,&*miit_ref_ent);
    e_miit = entsMoabField.insert(moabent);
    //this is nor real field in space (set order to zero)
    bool success = entsMoabField.modify(e_miit.first,MoFEMEntity_change_order(moab,0));
    if(!success) SETERRQ(PETSC_COMM_SELF,MOFEM_OPERATION_UNSUCCESSFUL,"modification unsuccessful");
  } catch (const char* msg) {
    SETERRQ(PETSC_COMM_SELF,MOFEM_CHAR_THROW,msg);
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
  //typedef RefMoFEMEntity_multiIndex::index<Ent_mi_tag>::type ref_ents_by_ents;
  //find field
  const field_set_by_id &set_id = moabFields.get<BitFieldId_mi_tag>();
  field_set_by_id::iterator miit = set_id.find(id);
  if(miit == set_id.end()) SETERRQ(PETSC_COMM_SELF,MOFEM_NOT_FOUND,"field not found");
  //ents in the field meshset
  Range ents_of_id_meshset;
  rval = moab.get_entities_by_handle(miit->meshset,ents_of_id_meshset,false); CHKERR_PETSC(rval);
  //create dofsMoabField
  Range::iterator eit = ents_of_id_meshset.begin();
  for(;eit!=ents_of_id_meshset.end();eit++) {
    // check if ent is in ref meshset
    RefMoFEMEntity_multiIndex::index<Ent_mi_tag>::type::iterator miit_ref_ent = refinedEntities.get<Ent_mi_tag>().find(*eit);
    if(miit_ref_ent==refinedEntities.get<Ent_mi_tag>().end()) {
      RefMoFEMEntity ref_ent(moab,*eit);
      if(ref_ent.get_BitRefLevel().none()) continue; // not on any mesh and not in database 
      cerr << ref_ent << endl;
      cerr << "bit level " << ref_ent.get_BitRefLevel() << endl;
      SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INSONSISTENCY,"database inconsistency");
    }
    // create mofem entity linked to ref ent
    MoFEMEntity_multiIndex::iterator e_miit;
    try {
      e_miit = entsMoabField.find(MoFEMEntity(moab,&*miit,&*miit_ref_ent).get_global_unique_id());
    } catch (const char* msg) {
	SETERRQ(PETSC_COMM_SELF,MOFEM_CHAR_THROW ,msg);
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
      } catch (const char* msg) {
	SETERRQ(PETSC_COMM_SELF,MOFEM_CHAR_THROW,msg);
      } catch (const std::exception& ex) {
	ostringstream ss;
	ss << ex.what() << endl;
	SETERRQ(PETSC_COMM_SELF,MOFEM_OPERATION_UNSUCCESSFUL,ss.str().c_str());
      }
      if(!p_e_miit.second) SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INSONSISTENCY,"data inconsistency");
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
	      SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INSONSISTENCY,"data inconsistency");
	    }
	    if(d_miit.first->get_ent_type()!=e_miit->get_ent_type()) {
	      SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INSONSISTENCY,"data inconsistency");
	    }
	    if(d_miit.first->get_id()!=e_miit->get_id()) {
	      SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INSONSISTENCY,"data inconsistency");
	    }
	    //check dof
	    if(d_miit.first->get_dof_order()!=oo) {
	      ostringstream ss;
	      ss << "data inconsistency!" << endl;
	      ss << "should be " << mdof << endl;
	      ss << "but is " << *d_miit.first << endl;
	      SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INSONSISTENCY,ss.str().c_str());
	    }
	    if(d_miit.first->get_max_order()!=e_miit->get_max_order()) {
	      SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INSONSISTENCY,"data inconsistency");
	    }
	  } catch (const char* msg) {
	    SETERRQ(PETSC_COMM_SELF,MOFEM_CHAR_THROW ,msg);
	  }
	}
      }
    }
    if(DD != e_miit->get_max_rank()*e_miit->get_order_nb_dofs(e_miit->get_max_order())) {
      SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INSONSISTENCY,"data inconsistency");
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
    if(verbose>0) {
      PetscPrintf(comm,"Build Field %s\n",miit->get_name().c_str());
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
    if(verbose>0) {
      int _dof_counter_ = 0;
      for(map<EntityType,int>::iterator it = dof_counter.begin();
	it!=dof_counter.end();it++) {
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
	  case MBTET:
	    PetscSynchronizedPrintf(comm,"nb added dofs (tets) %d\n",it->second);
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
  PetscSynchronizedPrintf(comm,"Nb. dofs %u\n",dofsMoabField.size());
  *build_MoFEM = 1<<0;
  if(verbose>0) {
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
    PetscSynchronizedFlush(comm,PETSC_STDOUT); 
  }
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
  BitFEId id = get_BitFEId();
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
  } catch (const char* msg) {
    SETERRQ(PETSC_COMM_SELF,MOFEM_CHAR_THROW,msg);
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
  } catch (const char* msg) {
    SETERRQ(PETSC_COMM_SELF,MOFEM_CHAR_THROW,msg);
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
  } catch (const char* msg) {
    SETERRQ(PETSC_COMM_SELF,MOFEM_CHAR_THROW,msg);
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
  } catch (const char* msg) {
    SETERRQ(PETSC_COMM_SELF,MOFEM_CHAR_THROW,msg);
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
  } catch (const char* msg) {
    SETERRQ(PETSC_COMM_SELF,MOFEM_CHAR_THROW,msg);
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
  } catch (const char* msg) {
    SETERRQ(PETSC_COMM_SELF,MOFEM_CHAR_THROW,msg);
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
  typedef MoFEMProblem_multiIndex::index<Problem_mi_tag>::type moFEMProblems_by_name;
  const moFEMProblems_by_name& set = moFEMProblems.get<Problem_mi_tag>();
  moFEMProblems_by_name::iterator miit = set.find(name);
  if(miit==set.end()) {
    BitProblemId id = get_problem_shift();
    ierr = add_problem(id,name); CHKERRQ(ierr);
  } else if(bh == MF_EXCL) {
    SETERRQ1(PETSC_COMM_SELF,MOFEM_NOT_FOUND,"problem is in database %s",name.c_str());
  }
  PetscFunctionReturn(0);
}
BitProblemId Core::get_BitProblemId(const string& name) const {
  typedef MoFEMProblem_multiIndex::index<Problem_mi_tag>::type moFEMProblems_by_name;
  const moFEMProblems_by_name& set = moFEMProblems.get<Problem_mi_tag>();
  moFEMProblems_by_name::iterator miit = set.find(name);
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
  } catch (const char* msg) {
    SETERRQ(PETSC_COMM_SELF,MOFEM_CHAR_THROW,msg);
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
  } catch (const char* msg) {
    SETERRQ(PETSC_COMM_SELF,MOFEM_CHAR_THROW,msg);
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
  } catch (const char* msg) {
    SETERRQ(PETSC_COMM_SELF,MOFEM_CHAR_THROW,msg);
  }
  PetscFunctionReturn(0);
}
PetscErrorCode Core::add_ents_to_finite_element_by_TETs(const EntityHandle meshset,const BitFEId id,const bool recursive) {
  PetscFunctionBegin;
  *build_MoFEM &= 1<<0;
  EntityHandle idm = no_handle;
  try {
    idm = get_finite_element_meshset(id);
  } catch (const char* msg) {
    SETERRQ(PETSC_COMM_SELF,MOFEM_CHAR_THROW,msg);
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
  } catch (const char* msg) {
    SETERRQ(PETSC_COMM_SELF,MOFEM_CHAR_THROW,msg);
  }
  PetscFunctionReturn(0);
}
PetscErrorCode Core::add_ents_to_finite_element_by_TETs(const EntityHandle meshset,const string &name,const bool recursive) {
  PetscFunctionBegin;
  *build_MoFEM &= 1<<0;
  try {
    ierr = add_ents_to_finite_element_by_TETs(meshset,get_BitFEId(name),recursive);  CHKERRQ(ierr);
  } catch (const char* msg) {
    SETERRQ(PETSC_COMM_SELF,MOFEM_CHAR_THROW,msg);
  }
  PetscFunctionReturn(0);
}

PetscErrorCode Core::add_ents_to_finite_element_by_PRISMs(const EntityHandle meshset,const BitFEId id,const bool recursive) {
  PetscFunctionBegin;
  *build_MoFEM &= 1<<0;
  EntityHandle idm = no_handle;
  try {
    idm = get_finite_element_meshset(id);
  } catch (const char* msg) {
    SETERRQ(PETSC_COMM_SELF,MOFEM_CHAR_THROW,msg);
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
  } catch (const char* msg) {
    SETERRQ(PETSC_COMM_SELF,MOFEM_CHAR_THROW,msg);
  }
  PetscFunctionReturn(0);
}
PetscErrorCode Core::add_ents_to_finite_element_by_PRISMs(const EntityHandle meshset,const string &name,const bool recursive) {
  PetscFunctionBegin;
  *build_MoFEM &= 1<<0;
  try {
    ierr = add_ents_to_finite_element_by_PRISMs(meshset,get_BitFEId(name),recursive);  CHKERRQ(ierr);
  } catch (const char* msg) {
    SETERRQ(PETSC_COMM_SELF,MOFEM_CHAR_THROW,msg);
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
  } catch (const char* msg) {
    SETERRQ(PETSC_COMM_SELF,MOFEM_CHAR_THROW,msg);
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
    typedef MoFEMProblem_multiIndex::index<Problem_mi_tag>::type moFEMProblems_by_name;
    moFEMProblems_by_name& set = moFEMProblems.get<Problem_mi_tag>();
    moFEMProblems_by_name::iterator miit = set.find(name_problem);
    ostringstream ss;
    ss << name_problem;
    if(miit==set.end()) SETERRQ1(PETSC_COMM_SELF,MOFEM_NOT_FOUND,"this problem <%s> is not there",ss.str().c_str());
    BitFEId f_id = get_BitFEId(MoFEMFiniteElement_name);
    bool success = set.modify(miit,problem_MoFEMFiniteElement_change_bit_add(f_id));
    if(!success) SETERRQ(PETSC_COMM_SELF,MOFEM_OPERATION_UNSUCCESSFUL,"modification unsuccessful");
  } catch (const char* msg) {
    SETERRQ(PETSC_COMM_SELF,MOFEM_CHAR_THROW,msg);
  }
  PetscFunctionReturn(0);
}
PetscErrorCode Core::modify_problem_ref_level_add_bit(const string &name_problem,const BitRefLevel &bit) {
  PetscFunctionBegin;
  typedef MoFEMProblem_multiIndex::index<Problem_mi_tag>::type moFEMProblems_by_name;
  moFEMProblems_by_name& set = moFEMProblems.get<Problem_mi_tag>();
  moFEMProblems_by_name::iterator miit = set.find(name_problem);
  ostringstream ss;
  ss << name_problem;
  if(miit==set.end()) SETERRQ1(PETSC_COMM_SELF,MOFEM_NOT_FOUND,"this problem <%s> is there",ss.str().c_str());
  bool success = set.modify(miit,problem_change_ref_level_bit_add(bit));
  if(!success) SETERRQ(PETSC_COMM_SELF,1,"modification unsuccessful");
  PetscFunctionReturn(0);
}
PetscErrorCode Core::modify_problem_ref_level_set_bit(const string &name_problem,const BitRefLevel &bit) {
  PetscFunctionBegin;
  typedef MoFEMProblem_multiIndex::index<Problem_mi_tag>::type moFEMProblems_by_name;
  moFEMProblems_by_name& set = moFEMProblems.get<Problem_mi_tag>();
  moFEMProblems_by_name::iterator miit = set.find(name_problem);
  ostringstream ss;
  ss << name_problem;
  if(miit==set.end()) SETERRQ1(PETSC_COMM_SELF,MOFEM_NOT_FOUND,"this problem <%s> is there",ss.str().c_str());
  bool success = set.modify(miit,problem_change_ref_level_bit_set(bit));
  if(!success) SETERRQ(PETSC_COMM_SELF,MOFEM_OPERATION_UNSUCCESSFUL,"modification unsuccessful");
  PetscFunctionReturn(0);
}
PetscErrorCode Core::modify_problem_dof_mask_ref_level_set_bit(const string &name_problem,const BitRefLevel &bit) {
  PetscFunctionBegin;
  typedef MoFEMProblem_multiIndex::index<Problem_mi_tag>::type moFEMProblems_by_name;
  moFEMProblems_by_name& set = moFEMProblems.get<Problem_mi_tag>();
  moFEMProblems_by_name::iterator miit = set.find(name_problem);
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
  data_dofs.clear(); //clear data dofs multiindex //FIXME should be cleaned when dofs are cleaned form datasets
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
    } catch (const char* msg) {
      SETERRQ(PETSC_COMM_SELF,MOFEM_CHAR_THROW,msg);
    }
  }
  viit_data = data_view.get<1>().lower_bound(1); 
  if(data_dofs.size()!=(unsigned int)distance(viit_data,hi_viit_data)) {
    SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INSONSISTENCY,"data inconsistency");
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
    for(int ss = 0;ss<Last;ss++) id_common |= FEAdj_fields[ss]&BitFieldId().set(ii);
    if( id_common.none() ) continue;
    //find in database data associated with the field (ii)
    field_by_id::iterator miit = moabFields_by_id.find(BitFieldId().set(ii));
    if(miit==moabFields_by_id.end()) {
      SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INSONSISTENCY,"data incosistency");
    }
    //resolve entities on element, those entities are used to build tag with dof uids on finite element tag
    ierr = ent_fe.get_element_adjacency(moab,&*miit,adj_ents); CHKERRQ(ierr);
    //loop over adjacent to finite entities, and find dofs on those entities
    //this part is to build MoFEMFiniteElement_dof_uid_view
    Range::iterator eit2 = adj_ents.begin();
    for(;eit2!=adj_ents.end();eit2++) {
      ref_ent_by_ent::iterator ref_ent_miit = refinedEntities.get<Ent_mi_tag>().find(*eit2);
      if(ref_ent_miit==refinedEntities.get<Ent_mi_tag>().end()) {
	cerr << adj_ents << endl;
	cerr << ent_fe << endl;
	cerr << "bit level " << ent_fe.get_BitRefLevel() << endl;
	RefMoFEMEntity ref_ent(moab,*eit2);
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
	SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INSONSISTENCY,ss.str().c_str());
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
      SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INSONSISTENCY,"data insonsistency"); 
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
    if(verbose>0) PetscPrintf(comm,"Build Finite Elements %s\n",MoFEMFiniteElement_miit->get_name().c_str());
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
	SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INSONSISTENCY,ss.str().c_str());
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
      bool success = entFEAdjacencies.modify(
	p.first,MoFEMEntityEntMoFEMFiniteElementAdjacencyMap_change_ByWhat(BYROW));
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
      bool success = entFEAdjacencies.modify(
	p.first,MoFEMEntityEntMoFEMFiniteElementAdjacencyMap_change_ByWhat(BYCOL));
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
      bool success = entFEAdjacencies.modify(
	p.first,MoFEMEntityEntMoFEMFiniteElementAdjacencyMap_change_ByWhat(BYDATA));
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
PetscErrorCode Core::build_partitioned_problems(int verb) {
  PetscFunctionBegin;
  if(verb==-1) verb = verbose;
  if(!(*build_MoFEM&(1<<0))) SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INSONSISTENCY,"fields not build");
  if(!(*build_MoFEM&(1<<1))) SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INSONSISTENCY,"FEs not build");
  if(!(*build_MoFEM&(1<<2))) SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INSONSISTENCY,"entFEAdjacencies not build");
  DofMoFEMEntity_multiIndex_active_view dofs_rows,dofs_cols;
  MoFEMProblem_multiIndex::iterator p_miit = moFEMProblems.begin();
  for(;p_miit!=moFEMProblems.end();p_miit++) {
    if(p_miit->get_BitRefLevel().none()) {
      SETERRQ1(PETSC_COMM_SELF,1,"problem <%s> refinement level not set",p_miit->get_name().c_str());
    }
    //zero finite elements
    bool success = moFEMProblems.modify(p_miit,problem_clear_numered_finiteElementsPtr_change());
    if(!success) SETERRQ(PETSC_COMM_SELF,MOFEM_OPERATION_UNSUCCESSFUL,"modification unsuccessful");
    //miit2 iterator for finite elements
    EntMoFEMFiniteElement_multiIndex::iterator miit2 = finiteElementsMoFEMEnts.begin();
    EntMoFEMFiniteElement_multiIndex::iterator hi_miit2 = finiteElementsMoFEMEnts.end();
    //DofMoFEMEntity_multiIndex_active_view dofs_rows,dofs_cols;
    dofs_rows.clear();
    dofs_cols.clear();
    EntMoFEMFiniteElement_multiIndex::iterator miit3 = miit2;
    //iterate all finite elemen entities in database
    for(;miit3!=hi_miit2;miit3++) {
      //if element is in problem
      if((miit3->get_id()&p_miit->get_BitFEId()).any()) {
	//if finite element bit level has all refined bits sets
	if((miit3->get_BitRefLevel()&p_miit->get_BitRefLevel())==p_miit->get_BitRefLevel()) {
	  //get dof uids for rows and columns
	  ierr = miit3->get_MoFEMFiniteElement_row_dof_view(dofsMoabField,dofs_rows); CHKERRQ(ierr);
	  ierr = miit3->get_MoFEMFiniteElement_col_dof_view(dofsMoabField,dofs_cols); CHKERRQ(ierr);
	}
      }
    }
    //zero rows
    MoFEMProblem& problem = const_cast<MoFEMProblem&>(*p_miit);
    (*(DofIdx*)problem.tag_nbdof_data_row) = 0;
    (*(DofIdx*)problem.tag_local_nbdof_data_row) = 0;
    (*(DofIdx*)problem.tag_ghost_nbdof_data_row) = 0;
    problem.numered_dofs_rows.clear();
    //zero cols
    (*(DofIdx*)problem.tag_nbdof_data_col) = 0;
    (*(DofIdx*)problem.tag_local_nbdof_data_col) = 0;
    (*(DofIdx*)problem.tag_ghost_nbdof_data_col) = 0;
    problem.numered_dofs_cols.clear();
    //add dofs for rows
    DofMoFEMEntity_multiIndex_active_view::nth_index<1>::type::iterator miit4,hi_miit4;
    miit4 = dofs_rows.get<1>().lower_bound(1);
    hi_miit4 = dofs_rows.get<1>().upper_bound(1);
    for(;miit4!=hi_miit4;miit4++) {
      if(((*miit4)->get_BitRefLevel()&p_miit->get_DofMask_BitRefLevel())!=(*miit4)->get_BitRefLevel()) {
	continue;
      }
      NumeredDofMoFEMEntity dof(&**miit4);
      pair<NumeredDofMoFEMEntity_multiIndex::iterator,bool> p = problem.numered_dofs_rows.insert(dof); 
      int owner_proc = p.first->get_owner_proc();
      if(owner_proc<0) {
	SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INSONSISTENCY,"data inconsistency");
      }
      bool success = problem.numered_dofs_rows.modify(p.first,NumeredDofMoFEMEntity_part_change(owner_proc,-1));
      if(!success) SETERRQ(PETSC_COMM_SELF,MOFEM_OPERATION_UNSUCCESSFUL,"modification unsuccessful");
    }
    //add dofs for cols
    DofMoFEMEntity_multiIndex_active_view::nth_index<1>::type::iterator miit5,hi_miit5;
    miit5 = dofs_cols.get<1>().lower_bound(1);
    hi_miit5 = dofs_cols.get<1>().upper_bound(1);
    for(;miit5!=hi_miit5;miit5++) {
      if(((*miit5)->get_BitRefLevel()&p_miit->get_DofMask_BitRefLevel())!=(*miit5)->get_BitRefLevel()) {
	continue;
      }
      NumeredDofMoFEMEntity dof(&**miit5);
      pair<NumeredDofMoFEMEntity_multiIndex::iterator,bool> p = problem.numered_dofs_cols.insert(dof); 
      int owner_proc = p.first->get_owner_proc();
      if(owner_proc<0) {
	SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INSONSISTENCY,"data inconsistency");
      }
      bool success = problem.numered_dofs_cols.modify(p.first,NumeredDofMoFEMEntity_part_change(owner_proc,-1));
      if(!success) SETERRQ(PETSC_COMM_SELF,MOFEM_OPERATION_UNSUCCESSFUL,"modification unsuccessful");
    }
    //number dofs
    NumeredDofMoFEMEntity_multiIndex *numered_dofs[] = { 
      &(problem.numered_dofs_rows), &(problem.numered_dofs_cols) };
    DofIdx *nb_global_dofs[] = {
      (DofIdx*)problem.tag_nbdof_data_row, (DofIdx*)problem.tag_nbdof_data_col };
    DofIdx *nb_local_dofs[] = {
      (DofIdx*)p_miit->tag_local_nbdof_data_row, (DofIdx*)p_miit->tag_local_nbdof_data_col };
    for(int ss = 0;ss<2;ss++) {
      NumeredDofMoFEMEntity_multiIndex::index<Part_mi_tag>::type::iterator diit,hi_diit;
      diit = numered_dofs[ss]->get<Part_mi_tag>().lower_bound(rAnk);
      hi_diit = numered_dofs[ss]->get<Part_mi_tag>().upper_bound(rAnk);
      //cerr << "OOOOOOOOOOOOOOOOOOO " << distance(diit,hi_diit) << " r " << rank << endl;
      vector<GlobalUId> uids;
      uids.resize(distance(diit,hi_diit));
      int &local_idx = *(nb_local_dofs[ss]);
      local_idx = 0;
      for(int ii = 0;diit!=hi_diit;diit++,ii++) {
        bool success;
        success = numered_dofs[ss]->modify(
	  numered_dofs[ss]->project<0>(diit),NumeredDofMoFEMEntity_local_idx_change(local_idx++));
        if(!success) SETERRQ(PETSC_COMM_SELF,MOFEM_OPERATION_UNSUCCESSFUL,"modification unsuccessful");
        uids[ii] = diit->get_global_unique_id();
      }
      {
        IS is,isout;
        int *petscint_ptr = (int*)&*uids.begin();
        ierr = ISCreateGeneral(comm,
	  uids.size()*sizeof(GlobalUId)/sizeof(int),
	 petscint_ptr,PETSC_USE_POINTER,&is); CHKERRQ(ierr);
        ierr = ISAllGather(is,&isout); CHKERRQ(ierr);
        int isout_size;
        ierr = ISGetSize(isout,&isout_size); CHKERRQ(ierr);
        isout_size = isout_size*sizeof(int)/sizeof(GlobalUId);
	*(nb_global_dofs[ss]) = isout_size;
        {	
	  const int *ptr;
	  ierr = ISGetIndices(isout,&ptr); CHKERRQ(ierr);
	  GlobalUId *uid_ptr = const_cast<GlobalUId*>((const GlobalUId*)ptr);
	  NumeredDofMoFEMEntity_multiIndex::iterator diit,hi_diit;
	  diit = numered_dofs[ss]->begin();
	  hi_diit = numered_dofs[ss]->end();
	  for(;diit!=hi_diit;diit++) {
	    GlobalUId uid = diit->get_global_unique_id();
	    GlobalUId *ptr;
	    ptr = find(uid_ptr,uid_ptr+isout_size,uid);
	    if(ptr == uid_ptr+isout_size) {
	      SETERRQ(PETSC_COMM_SELF,MOFEM_OPERATION_UNSUCCESSFUL,"could not find uid");
	    }
	    DofIdx idx = distance(uid_ptr,ptr);
	    int part = diit->get_part();
	    bool success;
	    success = numered_dofs[ss]->modify(
	      diit,NumeredDofMoFEMEntity_mofem_index_change(idx));
	    if(!success) SETERRQ(PETSC_COMM_SELF,MOFEM_OPERATION_UNSUCCESSFUL,"modification unsuccessful");
	    success = numered_dofs[ss]->modify(
	      diit,NumeredDofMoFEMEntity_part_change(part,idx));
	    if(!success) SETERRQ(PETSC_COMM_SELF,MOFEM_OPERATION_UNSUCCESSFUL,"modification unsuccessful");
	  }
	  ierr = ISRestoreIndices(isout,&ptr); CHKERRQ(ierr);
        }
        ierr = ISDestroy(&is); CHKERRQ(ierr);
        ierr = ISDestroy(&isout); CHKERRQ(ierr);
      }
    }
    if(verb>2) {
      NumeredDofMoFEMEntity_multiIndex::iterator it,hi_it;
      for(int ss = 0;ss<2;ss++) {
	it = numered_dofs[ss]->begin();
	hi_it = numered_dofs[ss]->end();
	for(;it!=hi_it;it++) {
	  ostringstream zz;
	  zz << ss << " " << "rank " << rAnk << " ";
	  zz << *it << endl;
	  PetscSynchronizedPrintf(comm,zz.str().c_str());
	}
      }
    }
    if(verb>0) {
      PetscSynchronizedPrintf(comm,"Problem %s Nb. rows %u  Nb. cols %u\n",
	p_miit->get_name().c_str(),
	p_miit->get_nb_dofs_row(),p_miit->get_nb_dofs_col());
      PetscSynchronizedFlush(comm,PETSC_STDOUT); 
    }
  }
  *build_MoFEM |= 1<<3;
  *build_MoFEM |= 1<<4;
  PetscFunctionReturn(0);
}
PetscErrorCode Core::build_problems(int verb) {
  PetscFunctionBegin;
  if(verb==-1) verb = verbose;
  if(!(*build_MoFEM&(1<<0))) SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INSONSISTENCY,"fields not build");
  if(!(*build_MoFEM&(1<<1))) SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INSONSISTENCY,"FEs not build");
  if(!(*build_MoFEM&(1<<2))) SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INSONSISTENCY,"entFEAdjacencies not build");
  //iterate problems
  DofMoFEMEntity_multiIndex_active_view dofs_rows,dofs_cols;
  MoFEMProblem_multiIndex::iterator p_miit = moFEMProblems.begin();
  for(;p_miit!=moFEMProblems.end();p_miit++) {
    if(p_miit->get_BitRefLevel().none()) {
      SETERRQ1(PETSC_COMM_SELF,1,"problem <%s> refinement level not set",p_miit->get_name().c_str());
    }
    //zero finite elements
    bool success = moFEMProblems.modify(p_miit,problem_clear_numered_finiteElementsPtr_change());
    if(!success) SETERRQ(PETSC_COMM_SELF,MOFEM_OPERATION_UNSUCCESSFUL,"modification unsuccessful");
    //miit2 iterator for finite elements
    EntMoFEMFiniteElement_multiIndex::iterator miit2 = finiteElementsMoFEMEnts.begin();
    EntMoFEMFiniteElement_multiIndex::iterator hi_miit2 = finiteElementsMoFEMEnts.end();
    //DofMoFEMEntity_multiIndex_active_view dofs_rows,dofs_cols;
    dofs_rows.clear();
    dofs_cols.clear();
    EntMoFEMFiniteElement_multiIndex::iterator miit3 = miit2;
    //iterate all finite elemen entities in database
    for(;miit3!=hi_miit2;miit3++) {
      //if element is in problem
      if((miit3->get_id()&p_miit->get_BitFEId()).any()) {
	//if finite element bit level has all refined bits sets
	if((miit3->get_BitRefLevel()&p_miit->get_BitRefLevel())==p_miit->get_BitRefLevel()) {
	  //get dof uids for rows and columns
	  ierr = miit3->get_MoFEMFiniteElement_row_dof_view(dofsMoabField,dofs_rows); CHKERRQ(ierr);
	  ierr = miit3->get_MoFEMFiniteElement_col_dof_view(dofsMoabField,dofs_cols); CHKERRQ(ierr);
	}
      }
    }
    //zero rows
    success = moFEMProblems.modify(p_miit,problem_zero_nb_rows_change());
    if(!success) SETERRQ(PETSC_COMM_SELF,MOFEM_OPERATION_UNSUCCESSFUL,"modification unsuccessful");
    //zero cols
    success = moFEMProblems.modify(p_miit,problem_zero_nb_cols_change());
    if(!success) SETERRQ(PETSC_COMM_SELF,MOFEM_OPERATION_UNSUCCESSFUL,"modification unsuccessful");
    //add dofs for rows
    DofMoFEMEntity_multiIndex_active_view::nth_index<1>::type::iterator miit4,hi_miit4;
    miit4 = dofs_rows.get<1>().lower_bound(1);
    hi_miit4 = dofs_rows.get<1>().upper_bound(1);
    for(;miit4!=hi_miit4;miit4++) {
      if(((*miit4)->get_BitRefLevel()&p_miit->get_DofMask_BitRefLevel())!=(*miit4)->get_BitRefLevel()) {
	continue;
      }
      success = moFEMProblems.modify(p_miit,problem_row_change(&**miit4));
      if(!success) SETERRQ(PETSC_COMM_SELF,MOFEM_OPERATION_UNSUCCESSFUL,"modification unsuccessful");
    }
    //add dofs for cols
    DofMoFEMEntity_multiIndex_active_view::nth_index<1>::type::iterator miit5,hi_miit5;
    miit5 = dofs_cols.get<1>().lower_bound(1);
    hi_miit5 = dofs_cols.get<1>().upper_bound(1);
    for(;miit5!=hi_miit5;miit5++) {
      if(((*miit5)->get_BitRefLevel()&p_miit->get_DofMask_BitRefLevel())!=(*miit5)->get_BitRefLevel()) {
	continue;
      }
      success = moFEMProblems.modify(p_miit,problem_col_change(&**miit5));
      if(!success) SETERRQ(PETSC_COMM_SELF,MOFEM_OPERATION_UNSUCCESSFUL,"modification unsuccessful");
    }
    //number dofs on rows and collumns
    success = moFEMProblems.modify(p_miit,problem_row_number_change());
    if(!success) SETERRQ(PETSC_COMM_SELF,MOFEM_OPERATION_UNSUCCESSFUL,"modification unsuccessful");
    success = moFEMProblems.modify(p_miit,problem_col_number_change());
    if(!success) SETERRQ(PETSC_COMM_SELF,MOFEM_OPERATION_UNSUCCESSFUL,"modification unsuccessful");
    //job done, some debugging and postprocessing
    if(verbose>0) {
      PetscSynchronizedPrintf(comm,"Problem %s Nb. rows %u Nb. cols %u\n",
	p_miit->get_name().c_str(),
	p_miit->numered_dofs_rows.size(),p_miit->numered_dofs_cols.size());
    }
    if(verb>2) {
      EntMoFEMFiniteElement_multiIndex::iterator miit_ss = miit2;
      ostringstream ss;
      ss << "rank " << rAnk << " ";
      ss << "FEs data for problem " << *p_miit << endl;
      for(;miit_ss!=hi_miit2;miit_ss++) {
	ss << "rank " << rAnk << " ";
	ss << *miit_ss << endl;
      }
      ss << "rank " << rAnk << " ";
      ss << "FEs row dofs "<< *p_miit << " Nb. row dof " << p_miit->get_nb_dofs_row() << endl;
      NumeredDofMoFEMEntity_multiIndex::iterator miit_dd_row = p_miit->numered_dofs_rows.begin();
      for(;miit_dd_row!=p_miit->numered_dofs_rows.end();miit_dd_row++) {
	ss << "rank " << rAnk << " ";
	ss<<*miit_dd_row<<endl;
      }
      ss << "rank " << rAnk << " ";
      ss << "FEs col dofs "<< *p_miit << " Nb. col dof " << p_miit->get_nb_dofs_col() << endl;
      NumeredDofMoFEMEntity_multiIndex::iterator miit_dd_col = p_miit->numered_dofs_cols.begin();
      for(;miit_dd_col!=p_miit->numered_dofs_cols.end();miit_dd_col++) {
	ss << "rank " << rAnk << " ";
	ss<<*miit_dd_col<<endl;
      }
      PetscSynchronizedPrintf(comm,ss.str().c_str());
    }
  }
  if(verb>0) {
    PetscSynchronizedFlush(comm,PETSC_STDOUT); 
  }
  *build_MoFEM |= 1<<3;
  PetscFunctionReturn(0);
}
PetscErrorCode Core::clear_problems(int verb) {
  PetscFunctionBegin;
  if(verb==-1) verb = verbose;
  MoFEMProblem_multiIndex::iterator p_miit = moFEMProblems.begin();
  //iterate problems
  for(;p_miit!=moFEMProblems.end();p_miit++) {
    //zero rows
    bool success = moFEMProblems.modify(p_miit,problem_zero_nb_rows_change());
    if(!success) SETERRQ(PETSC_COMM_SELF,MOFEM_OPERATION_UNSUCCESSFUL,"modification unsuccessful");
    //zero cols
    success = moFEMProblems.modify(p_miit,problem_zero_nb_cols_change());
    if(!success) SETERRQ(PETSC_COMM_SELF,MOFEM_OPERATION_UNSUCCESSFUL,"modification unsuccessful");
    //clear finite elements
    success = moFEMProblems.modify(p_miit,problem_clear_numered_finiteElementsPtr_change());
    if(!success) SETERRQ(PETSC_COMM_SELF,MOFEM_OPERATION_UNSUCCESSFUL,"modification unsuccessful");
  }
  PetscFunctionReturn(0);
}
PetscErrorCode Core::simple_partition_problem(const string &name,int verb) {
  PetscFunctionBegin;
  if(verb==-1) verb = verbose;
  if(!(*build_MoFEM&(1<<0))) SETERRQ(PETSC_COMM_SELF,MOFEM_NOT_FOUND,"fields not build");
  if(!(*build_MoFEM&(1<<1))) SETERRQ(PETSC_COMM_SELF,MOFEM_NOT_FOUND,"FEs not build");
  if(!(*build_MoFEM&(1<<2))) SETERRQ(PETSC_COMM_SELF,MOFEM_NOT_FOUND,"entFEAdjacencies not build");
  if(!(*build_MoFEM&(1<<3))) SETERRQ(PETSC_COMM_SELF,MOFEM_NOT_FOUND,"moFEMProblems not build");
  if(verb>0) {
    PetscPrintf(comm,"Simple partition problem %s\n",name.c_str());
  }
  // find p_miit
  typedef MoFEMProblem_multiIndex::index<Problem_mi_tag>::type moFEMProblems_by_name;
  moFEMProblems_by_name &moFEMProblems_set = moFEMProblems.get<Problem_mi_tag>();
  moFEMProblems_by_name::iterator p_miit = moFEMProblems_set.find(name);
  if(p_miit==moFEMProblems_set.end()) SETERRQ1(PETSC_COMM_SELF,1,"problem < %s > is not found (top tip: check spelling)",name.c_str());
  typedef boost::multi_index::index<NumeredDofMoFEMEntity_multiIndex,Idx_mi_tag>::type NumeredDofMoFEMEntitys_by_idx;
  NumeredDofMoFEMEntitys_by_idx &dofs_row_by_idx = const_cast<NumeredDofMoFEMEntitys_by_idx&>(p_miit->numered_dofs_rows.get<Idx_mi_tag>());
  NumeredDofMoFEMEntitys_by_idx &dofs_col_by_idx = const_cast<NumeredDofMoFEMEntitys_by_idx&>(p_miit->numered_dofs_cols.get<Idx_mi_tag>());
  boost::multi_index::index<NumeredDofMoFEMEntity_multiIndex,Idx_mi_tag>::type::iterator miit_row,hi_miit_row;
  boost::multi_index::index<NumeredDofMoFEMEntity_multiIndex,Idx_mi_tag>::type::iterator miit_col,hi_miit_col;
  DofIdx &nb_row_local_dofs = *((DofIdx*)p_miit->tag_local_nbdof_data_row);
  DofIdx &nb_row_ghost_dofs = *((DofIdx*)p_miit->tag_ghost_nbdof_data_row);
  nb_row_local_dofs = 0;
  nb_row_ghost_dofs = 0;
  DofIdx &nb_col_local_dofs = *((DofIdx*)p_miit->tag_local_nbdof_data_col);
  DofIdx &nb_col_ghost_dofs = *((DofIdx*)p_miit->tag_ghost_nbdof_data_col);
  nb_col_local_dofs = 0;
  nb_col_ghost_dofs = 0;
  //get row range of local indices
  DofIdx nb_dofs_row = dofs_row_by_idx.size();
  PetscLayout layout_row;
  ierr = PetscLayoutCreate(comm,&layout_row); CHKERRQ(ierr);
  ierr = PetscLayoutSetBlockSize(layout_row,1); CHKERRQ(ierr);
  ierr = PetscLayoutSetSize(layout_row,nb_dofs_row); CHKERRQ(ierr);
  ierr = PetscLayoutSetUp(layout_row); CHKERRQ(ierr);
  const int *ranges_row;
  ierr = PetscLayoutGetRanges(layout_row,&ranges_row); CHKERRQ(ierr);
  //get col range of local indices
  DofIdx nb_dofs_col = dofs_col_by_idx.size();
  PetscLayout layout_col;
  ierr = PetscLayoutCreate(comm,&layout_col); CHKERRQ(ierr);
  ierr = PetscLayoutSetBlockSize(layout_col,1); CHKERRQ(ierr);
  ierr = PetscLayoutSetSize(layout_col,nb_dofs_col); CHKERRQ(ierr);
  ierr = PetscLayoutSetUp(layout_col); CHKERRQ(ierr);
  const int *ranges_col;
  ierr = PetscLayoutGetRanges(layout_col,&ranges_col); CHKERRQ(ierr);
  int size;
  MPI_Comm_size(comm,&size);
  int rank;
  MPI_Comm_rank(comm,&rank);
  for(unsigned int part = 0;part<size;part++) {
    miit_row = dofs_row_by_idx.lower_bound(ranges_row[part]);
    hi_miit_row = dofs_row_by_idx.lower_bound(ranges_row[part+1]);
    if(distance(miit_row,hi_miit_row) != ranges_row[part+1]-ranges_row[part]) {
      SETERRQ4(PETSC_COMM_SELF,PETSC_ERR_ARG_SIZ,
	  "data inconsistency, distance(miit_row,hi_miit_row) != rend - rstart (%d != %d - %d = %d) ",
	  distance(miit_row,hi_miit_row),ranges_row[part+1],ranges_row[part],ranges_row[part+1]-ranges_row[part]);
    }
    // loop rows
    for(;miit_row!=hi_miit_row;miit_row++) {
      bool success = dofs_row_by_idx.modify(miit_row,NumeredDofMoFEMEntity_part_change(part,miit_row->dof_idx));
      if(!success) SETERRQ(PETSC_COMM_SELF,MOFEM_OPERATION_UNSUCCESSFUL,"modification unsuccessful");
      if(part == rank) {
	success = dofs_row_by_idx.modify(miit_row,NumeredDofMoFEMEntity_local_idx_change(nb_row_local_dofs++));
	if(!success) SETERRQ(PETSC_COMM_SELF,MOFEM_OPERATION_UNSUCCESSFUL,"modification unsuccessful");
      }
    }
    miit_col = dofs_col_by_idx.lower_bound(ranges_col[part]);
    hi_miit_col = dofs_col_by_idx.lower_bound(ranges_col[part+1]);
    if(distance(miit_col,hi_miit_col) != ranges_col[part+1]-ranges_col[part]) {
      SETERRQ4(PETSC_COMM_SELF,PETSC_ERR_ARG_SIZ,
	  "data inconsistency, distance(miit_col,hi_miit_col) != rend - rstart (%d != %d - %d = %d) ",
	  distance(miit_col,hi_miit_col),ranges_col[part+1],ranges_col[part],ranges_col[part+1]-ranges_col[part]);
    }
    // loop cols
    for(;miit_col!=hi_miit_col;miit_col++) {
      bool success = dofs_col_by_idx.modify(miit_col,NumeredDofMoFEMEntity_part_change(part,miit_col->dof_idx));
      if(!success) SETERRQ(PETSC_COMM_SELF,MOFEM_OPERATION_UNSUCCESSFUL,"modification unsuccessful");
      if(part == rank) {
	success = dofs_col_by_idx.modify(miit_col,NumeredDofMoFEMEntity_local_idx_change(nb_col_local_dofs++));
	if(!success) SETERRQ(PETSC_COMM_SELF,MOFEM_OPERATION_UNSUCCESSFUL,"modification unsuccessful");
      }
    }
  }
  ierr = PetscLayoutDestroy(&layout_row); CHKERRQ(ierr);
  ierr = PetscLayoutDestroy(&layout_col); CHKERRQ(ierr);
  if(verbose>0) {
    ostringstream ss;
    ss << "simple_partition_problem: rank = " << rank << " FEs row ghost dofs "<< *p_miit 
      << " Nb. local dof " << p_miit->get_nb_local_dofs_row() << " nb global row dofs " << p_miit->get_nb_dofs_row() << endl;
    ss << "simple_partition_problem: rank = " << rank << " FEs col ghost dofs " << *p_miit 
      << " Nb. local dof " << p_miit->get_nb_local_dofs_col() << " nb global col dofs " << p_miit->get_nb_dofs_col() << endl;
    PetscSynchronizedPrintf(comm,ss.str().c_str());
    PetscSynchronizedFlush(comm,PETSC_STDOUT); 
  }
  *build_MoFEM |= 1<<4;
  PetscFunctionReturn(0);
}
PetscErrorCode Core::compose_problem(const string &name,const string &problem_for_rows,const string &problem_for_cols,int verb) {
  PetscFunctionBegin;
  ierr = compose_problem(name,problem_for_rows,false,problem_for_cols,false,verb); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
PetscErrorCode Core::compose_problem(const string &name,const string &problem_for_rows,bool copy_rows,const string &problem_for_cols,bool copy_cols,int verb) {
  PetscFunctionBegin;
  if(verb==-1) verb = verbose;
  if(!(*build_MoFEM&(1<<0))) SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INSONSISTENCY,"fields not build");
  if(!(*build_MoFEM&(1<<1))) SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INSONSISTENCY,"FEs not build");
  if(!(*build_MoFEM&(1<<2))) SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INSONSISTENCY,"entFEAdjacencies not build");
  if(!(*build_MoFEM&(1<<3))) SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INSONSISTENCY,"moFEMProblems not build");
  typedef MoFEMProblem_multiIndex::index<Problem_mi_tag>::type moFEMProblems_by_name;
  typedef NumeredDofMoFEMEntity_multiIndex::index<Unique_mi_tag>::type NumeredDofMoFEMEntitys_by_uid;
  //find p_miit
  moFEMProblems_by_name &moFEMProblems_set = moFEMProblems.get<Problem_mi_tag>();
  moFEMProblems_by_name::iterator p_miit = moFEMProblems_set.find(name);
  if(p_miit==moFEMProblems_set.end()) SETERRQ1(PETSC_COMM_SELF,1,"problem with name %s not defined (top tip check spelling)",name.c_str());
  if(verb>0) {
    PetscPrintf(comm,"Compose problem %s from rows of %s and columns of %s\n",
      p_miit->get_name().c_str(),problem_for_rows.c_str(),problem_for_cols.c_str());
  }
  //find p_miit_row
  moFEMProblems_by_name::iterator p_miit_row = moFEMProblems_set.find(problem_for_rows);
  if(p_miit_row==moFEMProblems_set.end()) SETERRQ1(PETSC_COMM_SELF,1,"problem with name %s not defined (top tip check spelling)",problem_for_rows.c_str());
  moFEMProblems_by_name::iterator p_miit_col = moFEMProblems_set.find(problem_for_cols);
  //find p_mit_col
  if(p_miit_col==moFEMProblems_set.end()) SETERRQ1(PETSC_COMM_SELF,1,"problem with name %s not defined (top tip check spelling)",problem_for_cols.c_str());
  const NumeredDofMoFEMEntity_multiIndex &dofs_col = p_miit_col->numered_dofs_cols;
  //do rows
  map<DofIdx,const NumeredDofMoFEMEntity*> rows_problem_map;
  const NumeredDofMoFEMEntity_multiIndex &dofs_row = p_miit_row->numered_dofs_rows;
  if(!copy_rows) {
    MoFEMEntity *MoFEMEntity_ptr = NULL;
    NumeredDofMoFEMEntity_multiIndex::iterator miit_row = dofs_row.begin();
    NumeredDofMoFEMEntity_multiIndex::iterator hi_miit_row = dofs_row.end();
    NumeredDofMoFEMEntity_multiIndex_uid_view_ordered row_dof_view;
    for(;miit_row!=hi_miit_row;miit_row++) {
      if( (MoFEMEntity_ptr == NULL) ? 1 : (MoFEMEntity_ptr->get_global_unique_id() != miit_row->field_ptr->field_ptr->get_global_unique_id()) ) {
	MoFEMEntity_ptr = const_cast<MoFEMEntity*>(miit_row->field_ptr->field_ptr);
	typedef MoFEMEntityEntMoFEMFiniteElementAdjacencyMap_multiIndex::index<Unique_mi_tag>::type adj_by_ent;
	adj_by_ent::iterator adj_miit = entFEAdjacencies.get<Unique_mi_tag>().lower_bound(MoFEMEntity_ptr->get_global_unique_id());
	adj_by_ent::iterator hi_adj_miit = entFEAdjacencies.get<Unique_mi_tag>().upper_bound(MoFEMEntity_ptr->get_global_unique_id());
	for(;adj_miit!=hi_adj_miit;adj_miit++) {
	  if(!(adj_miit->by_other&BYROW)) {
	    // if it is not row if element
	    continue; 
	  }
	  if((adj_miit->EntMoFEMFiniteElement_ptr->get_id()&p_miit->get_BitFEId()).none()) {
	    // if element is not part of prblem
	    continue; 
	  }
	  if((adj_miit->EntMoFEMFiniteElement_ptr->get_BitRefLevel()&miit_row->get_BitRefLevel()).none()) {
	    // if entity is not problem refinment level
	    continue; 
	  }
	  //NumeredDofMoFEMEntity_multiIndex_uid_view_ordered row_dof_view;
	  row_dof_view.clear();
	  ierr = adj_miit->EntMoFEMFiniteElement_ptr->get_MoFEMFiniteElement_row_dof_view( 
	    dofs_row,row_dof_view,Interface::UNION); CHKERRQ(ierr);
	  NumeredDofMoFEMEntity_multiIndex_uid_view_ordered::iterator rdvit;
	  rdvit = row_dof_view.begin();
	  for(;rdvit!=row_dof_view.end();rdvit++) {
	    DofIdx petsc_global_idx = (*rdvit)->get_petsc_gloabl_dof_idx();
	    rows_problem_map[petsc_global_idx] = (*rdvit);
	  }
	}
      }
    }
    if(rows_problem_map.size() != (unsigned int)p_miit->get_nb_dofs_row()) {
      SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INSONSISTENCY,"data inconsistency");
    }
  } 
  //do cols
  map<DofIdx,const NumeredDofMoFEMEntity*> cols_problem_map;
  if(!copy_cols) {
    MoFEMEntity *MoFEMEntity_ptr = NULL;
    NumeredDofMoFEMEntity_multiIndex::iterator miit_col = dofs_col.begin();
    NumeredDofMoFEMEntity_multiIndex::iterator hi_miit_col = dofs_col.end();
    NumeredDofMoFEMEntity_multiIndex_uid_view_ordered col_dof_view;
    for(;miit_col!=hi_miit_col;miit_col++) {
      if( (MoFEMEntity_ptr == NULL) ? 1 : (MoFEMEntity_ptr->get_global_unique_id() != miit_col->field_ptr->field_ptr->get_global_unique_id()) ) {
	MoFEMEntity_ptr = const_cast<MoFEMEntity*>(miit_col->field_ptr->field_ptr);
	typedef MoFEMEntityEntMoFEMFiniteElementAdjacencyMap_multiIndex::index<Unique_mi_tag>::type adj_by_ent;
	adj_by_ent::iterator adj_miit = entFEAdjacencies.get<Unique_mi_tag>().lower_bound(MoFEMEntity_ptr->get_global_unique_id());
	adj_by_ent::iterator hi_adj_miit = entFEAdjacencies.get<Unique_mi_tag>().upper_bound(MoFEMEntity_ptr->get_global_unique_id());
	for(;adj_miit!=hi_adj_miit;adj_miit++) {
	  if(!(adj_miit->by_other&BYCOL)) {
	    // if it is not row if element
	    continue; 
	  }
	  if((adj_miit->EntMoFEMFiniteElement_ptr->get_id()&p_miit->get_BitFEId()).none()) {
	    // if element is not part of prblem
	    continue; 
	  }
	  if((adj_miit->EntMoFEMFiniteElement_ptr->get_BitRefLevel()&miit_col->get_BitRefLevel()).none()) {
	    // if entity is not problem refinment level
	    continue; 
	  }
	  //NumeredDofMoFEMEntity_multiIndex_uid_view_ordered col_dof_view;
	  col_dof_view.clear();
	  ierr = adj_miit->EntMoFEMFiniteElement_ptr->get_MoFEMFiniteElement_col_dof_view( 
	    dofs_col,col_dof_view,Interface::UNION); CHKERRQ(ierr);
	  NumeredDofMoFEMEntity_multiIndex_uid_view_ordered::iterator cdvit;
	  cdvit = col_dof_view.begin();
	  for(;cdvit!=col_dof_view.end();cdvit++) {
	    DofIdx petsc_global_idx = (*cdvit)->get_petsc_gloabl_dof_idx();
	    cols_problem_map[petsc_global_idx] = (*cdvit);
	  }
	}
      }
    }
    if(cols_problem_map.size() != (unsigned int)p_miit->get_nb_dofs_col()) {
      SETERRQ2(PETSC_COMM_SELF,1,"data inconsistency, %u != %u",cols_problem_map.size(),(unsigned int)p_miit->get_nb_dofs_col());
    }
  }
  // build indices 
  NumeredDofMoFEMEntitys_by_uid &dofs_row_by_uid = const_cast<NumeredDofMoFEMEntitys_by_uid&>(p_miit->numered_dofs_rows.get<Unique_mi_tag>());
  DofIdx &nb_row_local_dofs = *((DofIdx*)p_miit->tag_local_nbdof_data_row);
  nb_row_local_dofs = 0;
  if(!copy_rows) {
  //rows
  map<DofIdx,const NumeredDofMoFEMEntity*>::iterator miit_map_row = rows_problem_map.begin();
  for(;miit_map_row!=rows_problem_map.end();miit_map_row++) {
    NumeredDofMoFEMEntitys_by_uid::iterator pr_dof = dofs_row_by_uid.find(miit_map_row->second->get_global_unique_id());
    if(pr_dof == dofs_row_by_uid.end()) SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INSONSISTENCY,"data inconsistency");
    int part_number = miit_map_row->second->get_part();
    int petsc_global_dof = distance(rows_problem_map.begin(),miit_map_row);
    if(verb>1) {
      PetscPrintf(comm,"Row Problem Glob Idx %d Problem Glob Idx %d\n",miit_map_row->second->get_petsc_gloabl_dof_idx(),petsc_global_dof);
    }
    bool success = dofs_row_by_uid.modify(pr_dof,NumeredDofMoFEMEntity_part_change(part_number,petsc_global_dof));
    if(!success) SETERRQ(PETSC_COMM_SELF,MOFEM_OPERATION_UNSUCCESSFUL,"modification unsuccessful");
    if(pr_dof->get_part() == rAnk) {
      success = dofs_row_by_uid.modify(pr_dof,NumeredDofMoFEMEntity_local_idx_change(nb_row_local_dofs++));
      if(!success) SETERRQ(PETSC_COMM_SELF,MOFEM_OPERATION_UNSUCCESSFUL,"modification unsuccessful");
    }
  }
  } else {
    for(_IT_NUMEREDDOFMOFEMENTITY_ROW_FOR_LOOP_(p_miit_row,diit)) {
      bool success = moFEMProblems.modify(moFEMProblems.project<0>(p_miit),problem_row_change(diit->get_DofMoFEMEntity_ptr()));
      if(!success) SETERRQ(PETSC_COMM_SELF,MOFEM_OPERATION_UNSUCCESSFUL,"modification unsuccessful");
      int part_number = diit->get_part();
      int petsc_global_dof = diit->get_petsc_gloabl_dof_idx();
      NumeredDofMoFEMEntitys_by_uid::iterator pr_dof = dofs_row_by_uid.find(diit->get_global_unique_id());
      if(pr_dof == dofs_row_by_uid.end()) SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INSONSISTENCY,"data inconsistency");
      success = dofs_row_by_uid.modify(pr_dof,NumeredDofMoFEMEntity_part_change(part_number,petsc_global_dof));
      if(!success) SETERRQ(PETSC_COMM_SELF,MOFEM_OPERATION_UNSUCCESSFUL,"modification unsuccessful");
      if(pr_dof->get_part() == rAnk) {
	success = dofs_row_by_uid.modify(pr_dof,NumeredDofMoFEMEntity_local_idx_change(nb_row_local_dofs++));
	if(!success) SETERRQ(PETSC_COMM_SELF,MOFEM_OPERATION_UNSUCCESSFUL,"modification unsuccessful");
      }
    }
  }
  //cols
  NumeredDofMoFEMEntitys_by_uid &dofs_col_by_uid = const_cast<NumeredDofMoFEMEntitys_by_uid&>(p_miit->numered_dofs_cols.get<Unique_mi_tag>());
  DofIdx &nb_col_local_dofs = *((DofIdx*)p_miit->tag_local_nbdof_data_col);
  nb_col_local_dofs = 0;
  if(!copy_cols) {
  map<DofIdx,const NumeredDofMoFEMEntity*>::iterator miit_map_col = cols_problem_map.begin();
  for(;miit_map_col!=cols_problem_map.end();miit_map_col++) {
    NumeredDofMoFEMEntitys_by_uid::iterator pr_dof = dofs_col_by_uid.find(miit_map_col->second->get_global_unique_id());
    if(pr_dof == dofs_col_by_uid.end()) SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INSONSISTENCY,"data inconsistency");
    int part_number = miit_map_col->second->get_part();
    int petsc_global_dof = distance(cols_problem_map.begin(),miit_map_col);
    if(verb>1) {
      PetscPrintf(comm,"Col Problem Glob Idx %d Problem Glob Idx %d\n",miit_map_col->second->get_petsc_gloabl_dof_idx(),petsc_global_dof);
    }
    bool success = dofs_col_by_uid.modify(pr_dof,NumeredDofMoFEMEntity_part_change(part_number,petsc_global_dof));
    if(!success) SETERRQ(PETSC_COMM_SELF,MOFEM_OPERATION_UNSUCCESSFUL,"modification unsuccessful");
    if(pr_dof->get_part() == rAnk) {
      success = dofs_col_by_uid.modify(pr_dof,NumeredDofMoFEMEntity_local_idx_change(nb_col_local_dofs++));
      if(!success) SETERRQ(PETSC_COMM_SELF,MOFEM_OPERATION_UNSUCCESSFUL,"modification unsuccessful");
    }
  }
  } else {
    for(_IT_NUMEREDDOFMOFEMENTITY_COL_FOR_LOOP_(p_miit_col,diit)) {
      bool success = moFEMProblems.modify(moFEMProblems.project<0>(p_miit),problem_col_change(diit->get_DofMoFEMEntity_ptr()));
      if(!success) SETERRQ(PETSC_COMM_SELF,MOFEM_OPERATION_UNSUCCESSFUL,"modification unsuccessful");
      int part_number = diit->get_part();
      int petsc_global_dof = diit->get_petsc_gloabl_dof_idx();
      NumeredDofMoFEMEntitys_by_uid::iterator pr_dof = dofs_col_by_uid.find(diit->get_global_unique_id());
      success = dofs_col_by_uid.modify(pr_dof,NumeredDofMoFEMEntity_part_change(part_number,petsc_global_dof));
      if(!success) SETERRQ(PETSC_COMM_SELF,MOFEM_OPERATION_UNSUCCESSFUL,"modification unsuccessful");
      if(pr_dof->get_part() == rAnk) {
	success = dofs_col_by_uid.modify(pr_dof,NumeredDofMoFEMEntity_local_idx_change(nb_col_local_dofs++));
	if(!success) SETERRQ(PETSC_COMM_SELF,MOFEM_OPERATION_UNSUCCESSFUL,"modification unsuccessful");
      }
    }
  }
  if(verbose>0) {
    ostringstream ss;
    ss << "partition_problem: rank = " << rAnk << " FEs row ghost dofs "<< *p_miit 
      << " Nb. local dof " << p_miit->get_nb_local_dofs_row() << " nb global row dofs " << p_miit->get_nb_dofs_row() << endl;
    ss << "partition_problem: rank = " << rAnk << " FEs col ghost dofs " << *p_miit 
      << " Nb. local dof " << p_miit->get_nb_local_dofs_col() << " nb global col dofs " << p_miit->get_nb_dofs_col() << endl;
    PetscSynchronizedPrintf(comm,ss.str().c_str());
    PetscSynchronizedFlush(comm,PETSC_STDOUT); 
  }
  if(verb>2) {
    ostringstream ss;
    ss << "rank = " << rAnk << " FEs row dofs "<< *p_miit << " Nb. row dof " << p_miit->get_nb_dofs_row() 
	<< " Nb. local dof " << p_miit->get_nb_local_dofs_row() << endl;
    NumeredDofMoFEMEntity_multiIndex::iterator miit_dd_row = p_miit->numered_dofs_rows.begin();
    for(;miit_dd_row!=p_miit->numered_dofs_rows.end();miit_dd_row++) {
	ss<<*miit_dd_row<<endl;
    }
    ss << "rank = " << rAnk << " FEs col dofs "<< *p_miit << " Nb. col dof " << p_miit->get_nb_dofs_col() 
	<< " Nb. local dof " << p_miit->get_nb_local_dofs_col() << endl;
    NumeredDofMoFEMEntity_multiIndex::iterator miit_dd_col = p_miit->numered_dofs_cols.begin();
    for(;miit_dd_col!=p_miit->numered_dofs_cols.end();miit_dd_col++) {
	ss<<*miit_dd_col<<endl;
    }
    PetscSynchronizedPrintf(comm,ss.str().c_str());
    PetscSynchronizedFlush(comm,PETSC_STDOUT); 
  }
  if(debug>0) {
    typedef NumeredDofMoFEMEntity_multiIndex::index<Idx_mi_tag>::type NumeredDofMoFEMEntitys_by_idx;
    NumeredDofMoFEMEntitys_by_idx::iterator dit,hi_dit;
    dit = p_miit->numered_dofs_rows.get<Idx_mi_tag>().begin();
    hi_dit = p_miit->numered_dofs_rows.get<Idx_mi_tag>().end();
    for(;dit!=hi_dit;dit++) {
      if(dit->get_part()==rAnk) {
	if(dit->get_petsc_local_dof_idx()<0) {
	  ostringstream ss;
	  ss << "rank " << rAnk << " " << *dit;
	  SETERRQ1(PETSC_COMM_SELF,1,"local dof index for row not set\n %s",ss.str().c_str());
	}
      }
    }
    dit = p_miit->numered_dofs_cols.get<Idx_mi_tag>().begin();
    hi_dit = p_miit->numered_dofs_cols.get<Idx_mi_tag>().end();
    for(;dit!=hi_dit;dit++) {
      if(dit->get_part()==rAnk) {
	if(dit->get_petsc_local_dof_idx()<0) {
	  ostringstream ss;
	  ss << "rank " << rAnk << " " << *dit;
	  SETERRQ1(PETSC_COMM_SELF,1,"local dof index for col not set\n %s",ss.str().c_str());
	}
      }
    }
  }
  PetscFunctionReturn(0);
}
PetscErrorCode Core::partition_finite_elements(
  const string &name,bool part_from_moab,int low_proc,int hi_proc,int verb) {
  PetscFunctionBegin;
  if(verb==-1) verb = verbose;
  if(!(*build_MoFEM&(1<<0))) SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INSONSISTENCY,"fields not build");
  if(!(*build_MoFEM&(1<<1))) SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INSONSISTENCY,"FEs not build");
  if(!(*build_MoFEM&(1<<2))) SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INSONSISTENCY,"entFEAdjacencies not build");
  if(!(*build_MoFEM&(1<<3))) SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INSONSISTENCY,"partitions not build");
  if(!(*build_MoFEM&(1<<4))) SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INSONSISTENCY,"partitions moFEMProblems not build");
  if(low_proc == -1) low_proc = rAnk;
  if(hi_proc == -1) hi_proc = rAnk;
  typedef MoFEMProblem_multiIndex::index<Problem_mi_tag>::type moFEMProblems_by_name;
  //find p_miit
  moFEMProblems_by_name &moFEMProblems_set = moFEMProblems.get<Problem_mi_tag>();
  moFEMProblems_by_name::iterator p_miit = moFEMProblems_set.find(name);
  if(p_miit == moFEMProblems_set.end()) SETERRQ1(PETSC_COMM_SELF,1,"problem < %s > not found (top tip: check spelling)",name.c_str());
  NumeredMoFEMFiniteElement_multiIndex& numeredFiniteElements = const_cast<NumeredMoFEMFiniteElement_multiIndex&>(p_miit->numeredFiniteElements);
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
      if( (numered_fe.get_part()>=(unsigned int)low_proc)&&(numered_fe.get_part()<=(unsigned int)hi_proc) ) {
	//rows element dof multiindices
	viit_rows = rows_view.begin();
	for(;viit_rows!=rows_view.end();viit_rows++) {
	  try {
	    SideNumber *side_number_ptr = miit2->get_side_number_ptr(moab,(*viit_rows)->get_ent());
	    rows_dofs.insert(boost::make_tuple(side_number_ptr,&**viit_rows));
	  } catch (const char* msg) {
	    SETERRQ(PETSC_COMM_SELF,MOFEM_CHAR_THROW,msg);
	  }
	}
	//cols_views
	NumeredDofMoFEMEntity_multiIndex_uid_view_ordered cols_view;
	ierr = miit2->get_MoFEMFiniteElement_col_dof_view(p_miit->numered_dofs_cols,cols_view,Interface::UNION); CHKERRQ(ierr);
	//cols elmeny dof multiindices
	NumeredDofMoFEMEntity_multiIndex_uid_view_ordered::iterator viit_cols;;
	viit_cols = cols_view.begin();
	for(;viit_cols!=cols_view.end();viit_cols++) {
	  try {
	    SideNumber *side_number_ptr = miit2->get_side_number_ptr(moab,(*viit_cols)->get_ent());
	    cols_dofs.insert(boost::make_tuple(side_number_ptr,&**viit_cols));
	  } catch (const char* msg) {
	    SETERRQ(PETSC_COMM_SELF,MOFEM_CHAR_THROW,msg);
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
  if(!(*build_MoFEM&(1<<0))) SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INSONSISTENCY,"fields not build");
  if(!(*build_MoFEM&(1<<1))) SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INSONSISTENCY,"FEs not build");
  if(!(*build_MoFEM&(1<<2))) SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INSONSISTENCY,"entFEAdjacencies not build");
  if(!(*build_MoFEM&(1<<3))) SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INSONSISTENCY,"moFEMProblems not build");
  if(!(*build_MoFEM&(1<<4))) SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INSONSISTENCY,"partitions moFEMProblems not build");
  if(!(*build_MoFEM&(1<<5))) SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INSONSISTENCY,"partitions finite elements not build");
  typedef MoFEMProblem_multiIndex::index<Problem_mi_tag>::type moFEMProblems_by_name;
  //find p_miit
  moFEMProblems_by_name &moFEMProblems_set = moFEMProblems.get<Problem_mi_tag>();
  moFEMProblems_by_name::iterator p_miit = moFEMProblems_set.find(name);
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
	  if(rowdofit->get_part()==rAnk) continue; 
	  ghost_idx_row_view.insert(rowdofit->get_NumeredDofMoFEMEntity_ptr());
	}
      }
      if(fe_it->cols_dofs.size()>0) {
	dof_it coldofit,hi_coldofit;
	coldofit = fe_it->cols_dofs.begin();
	hi_coldofit = fe_it->cols_dofs.end();
	for(;coldofit!=hi_coldofit;coldofit++) {
	  if(coldofit->get_part()==rAnk) continue;
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
	  SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INSONSISTENCY,"inconsistent data, ghost dof already set");
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
	if(miit_dd_col->part==rAnk) continue;
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
	SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INSONSISTENCY,"not implemented");
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
  } catch (const char* msg) {
    SETERRQ(PETSC_COMM_SELF,MOFEM_CHAR_THROW,msg);
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
  if(!((p_ent.first->get_BitRefLevel()&bit)==bit)) {
    refinedEntities.modify(p_ent.first,RefMoFEMEntity_change_add_bit(bit));
  }
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
    if(!success) SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INSONSISTENCY,"inconsistency in data");
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
  ierr = moab.get_entities_by_handle(0,ents,false); CHKERRQ(ierr);
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

  typedef RefMoFEMEntity_multiIndex::index<Composite_EntityHandle_And_ParentEntityType_mi_tag>::type ref_ents_by_composite;
  ref_ents_by_composite &ref_ents = refinedEntities.get<Composite_EntityHandle_And_ParentEntityType_mi_tag>();
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
PetscErrorCode Core::problem_get_FE(const string &problem_name,const string &fe_name,const EntityHandle meshset) {
  PetscFunctionBegin;
  typedef MoFEMProblem_multiIndex::index<Problem_mi_tag>::type moFEMProblems_by_name;
  moFEMProblems_by_name &moFEMProblems_set = moFEMProblems.get<Problem_mi_tag>();
  moFEMProblems_by_name::iterator p_miit = moFEMProblems_set.find(problem_name);
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
bool Core::check_msId_meshset(const int msId,const CubitBC_BitSet CubitBCType) {
  CubitMeshSet_multiIndex::index<Composite_Cubit_msId_and_MeshSetType_mi_tag>::type::iterator 
    miit = cubit_meshsets.get<Composite_Cubit_msId_and_MeshSetType_mi_tag>().find(boost::make_tuple(msId,CubitBCType.to_ulong()));
  if(miit!=cubit_meshsets.get<Composite_Cubit_msId_and_MeshSetType_mi_tag>().end()) {
    return true;
  } 
  return false;
}
PetscErrorCode Core::add_Cubit_msId(const CubitBC_BitSet CubitBCType,const int msId) {
  PetscFunctionBegin;
  if(check_msId_meshset(msId,CubitBCType)) {
    SETERRQ1(PETSC_COMM_SELF,1,"such cubit meshset is already there",msId);
  } 
  try {
    CubitMeshSets cmeshset(moab,CubitBCType,msId);
    if((cmeshset.CubitBCType&CubitBC_BitSet(NODESET|SIDESET|BLOCKSET)).any()) {
      pair<CubitMeshSet_multiIndex::iterator,bool> p = cubit_meshsets.insert(cmeshset);
      if(!p.second) SETERRQ(PETSC_COMM_SELF,MOFEM_OPERATION_UNSUCCESSFUL,"meshset not inserted");
    }
  } catch (const char* msg) {
    SETERRQ(PETSC_COMM_SELF,MOFEM_CHAR_THROW,msg);
  }
  PetscFunctionReturn(0);
}
PetscErrorCode Core::delete_Cubit_msId(const CubitBC_BitSet CubitBCType,const int msId) {
  PetscFunctionBegin;
  CubitMeshSet_multiIndex::index<Composite_Cubit_msId_and_MeshSetType_mi_tag>::type::iterator 
    miit = cubit_meshsets.get<Composite_Cubit_msId_and_MeshSetType_mi_tag>().find(boost::make_tuple(msId,CubitBCType.to_ulong()));
  if(miit==cubit_meshsets.get<Composite_Cubit_msId_and_MeshSetType_mi_tag>().end()) {
    SETERRQ1(PETSC_COMM_SELF,1,"such cubit meshset is already there",msId);
  }
  EntityHandle meshset = miit->get_meshset();
  cubit_meshsets.get<Composite_Cubit_msId_and_MeshSetType_mi_tag>().erase(miit);
  rval = moab.delete_entities(&meshset,1); CHKERR_PETSC(rval);
  PetscFunctionReturn(0);
}
PetscErrorCode Core::get_Cubit_msId(const int msId,const CubitBC_BitSet CubitBCType,const CubitMeshSets **cubit_meshset_ptr) {
  PetscFunctionBegin;
  CubitMeshSet_multiIndex::index<Composite_Cubit_msId_and_MeshSetType_mi_tag>::type::iterator 
    miit = cubit_meshsets.get<Composite_Cubit_msId_and_MeshSetType_mi_tag>().find(boost::make_tuple(msId,CubitBCType.to_ulong()));
  if(miit!=cubit_meshsets.get<Composite_Cubit_msId_and_MeshSetType_mi_tag>().end()) {
    *cubit_meshset_ptr = &*miit;
  } else {
    SETERRQ1(PETSC_COMM_SELF,1,"msId = %d is not there",msId);
  }
  PetscFunctionReturn(0);
}
PetscErrorCode Core::get_Cubit_msId_entities_by_dimension(const int msId,const CubitBC_BitSet CubitBCType,
  const int dimension,Range &entities,const bool recursive) {
  PetscFunctionBegin;
  CubitMeshSet_multiIndex::index<Composite_Cubit_msId_and_MeshSetType_mi_tag>::type::iterator 
    miit = cubit_meshsets.get<Composite_Cubit_msId_and_MeshSetType_mi_tag>().find(boost::make_tuple(msId,CubitBCType.to_ulong()));
  if(miit!=cubit_meshsets.get<Composite_Cubit_msId_and_MeshSetType_mi_tag>().end()) {
    ierr = miit->get_Cubit_msId_entities_by_dimension(moab,dimension,entities,recursive); CHKERRQ(ierr);
  } else {
    SETERRQ1(PETSC_COMM_SELF,1,"msId = %d is not there",msId);
  }
  PetscFunctionReturn(0);
}
PetscErrorCode Core::get_Cubit_msId_entities_by_dimension(const int msId,const CubitBC_BitSet CubitBCType,Range &entities,const bool recursive) {
  PetscFunctionBegin;
  CubitMeshSet_multiIndex::index<Composite_Cubit_msId_and_MeshSetType_mi_tag>::type::iterator 
    miit = cubit_meshsets.get<Composite_Cubit_msId_and_MeshSetType_mi_tag>().find(boost::make_tuple(msId,CubitBCType.to_ulong()));
  if(miit!=cubit_meshsets.get<Composite_Cubit_msId_and_MeshSetType_mi_tag>().end()) {
    ierr = miit->get_Cubit_msId_entities_by_dimension(moab,entities,recursive); CHKERRQ(ierr);
  } else {
    SETERRQ1(PETSC_COMM_SELF,1,"msId = %d is not there",msId);
  }
  PetscFunctionReturn(0);
}
PetscErrorCode Core::get_Cubit_msId_entities_by_dimension(const int msId,const unsigned int CubitBCType,
  const int dimension,Range &entities,const bool recursive) {
  PetscFunctionBegin;
  ierr = get_Cubit_msId_entities_by_dimension(msId,CubitBC_BitSet(CubitBCType),dimension,entities,recursive); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
PetscErrorCode Core::get_Cubit_msId_entities_by_dimension(const int msId,const unsigned int CubitBCType,
  Range &entities,const bool recursive) {
  PetscFunctionBegin;
  ierr = get_Cubit_msId_entities_by_dimension(msId,CubitBC_BitSet(CubitBCType),entities,recursive); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
PetscErrorCode Core::get_Cubit_msId_meshset(const int msId,const unsigned int CubitBCType,EntityHandle &meshset) {
  PetscFunctionBegin;
  CubitMeshSet_multiIndex::index<Composite_Cubit_msId_and_MeshSetType_mi_tag>::type::iterator 
    miit = cubit_meshsets.get<Composite_Cubit_msId_and_MeshSetType_mi_tag>().find(boost::make_tuple(msId,CubitBCType));
  if(miit!=cubit_meshsets.get<Composite_Cubit_msId_and_MeshSetType_mi_tag>().end()) {
    meshset = miit->meshset;
  } else {
    SETERRQ1(PETSC_COMM_SELF,1,"msId = %d is not there",msId);
  }
  PetscFunctionReturn(0);
}
PetscErrorCode Core::get_Cubit_meshsets(const unsigned int CubitBCType,Range &meshsets) {
  PetscFunctionBegin;
  CubitMeshSet_multiIndex::index<CubitMeshSets_mi_tag>::type::iterator 
    miit = cubit_meshsets.get<CubitMeshSets_mi_tag>().lower_bound(CubitBCType);
  CubitMeshSet_multiIndex::index<CubitMeshSets_mi_tag>::type::iterator 
    hi_miit = cubit_meshsets.get<CubitMeshSets_mi_tag>().upper_bound(CubitBCType);
  for(;miit!=hi_miit;miit++) {
    meshsets.insert(miit->meshset);
  }
  PetscFunctionReturn(0);
}

#define SET_BASIC_METHOD(PROBLEM_PTR) \
  { \
    method.problemPtr = PROBLEM_PTR; \
    method.fieldsPtr = &moabFields; \
    method.refinedEntitiesPtr = &refinedEntities; \
    method.entitiesPtr = &entsMoabField; \
    method.dofsPtr = &dofsMoabField; \
    method.refinedFiniteElementsPtr = &refinedFiniteElements; \
    method.finiteElementsPtr = &finiteElements; \
    method.finiteElementsEntitiesPtr = &finiteElementsMoFEMEnts; \
    method.adjacenciesPtr = &entFEAdjacencies; \
  }

PetscErrorCode Core::problem_basic_method_preProcess(const MoFEMProblem *problem_ptr,BasicMethod &method,int verb) {
  PetscFunctionBegin;
  if(verb==-1) verb = verbose;
  // finite element
  SET_BASIC_METHOD(problem_ptr)
  PetscLogEventBegin(USER_EVENT_preProcess,0,0,0,0);
  ierr = method.preProcess(); CHKERRQ(ierr);
  PetscLogEventEnd(USER_EVENT_preProcess,0,0,0,0);
  PetscFunctionReturn(0);
}

PetscErrorCode Core::problem_basic_method_preProcess(const string &problem_name,BasicMethod &method,int verb) {
  PetscFunctionBegin;
  if(verb==-1) verb = verbose;
  typedef MoFEMProblem_multiIndex::index<Problem_mi_tag>::type moFEMProblems_by_name;
  // find p_miit
  moFEMProblems_by_name &moFEMProblems_set = moFEMProblems.get<Problem_mi_tag>();
  moFEMProblems_by_name::iterator p_miit = moFEMProblems_set.find(problem_name);
  if(p_miit == moFEMProblems_set.end()) SETERRQ1(PETSC_COMM_SELF,1,"problem is not in database %s",problem_name.c_str());
  ierr = problem_basic_method_preProcess(&*p_miit,method,verb); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
PetscErrorCode Core::problem_basic_method_postProcess(const MoFEMProblem *problem_ptr,BasicMethod &method,int verb) {
  PetscFunctionBegin;
  SET_BASIC_METHOD(problem_ptr)
  PetscLogEventBegin(USER_EVENT_postProcess,0,0,0,0);
  ierr = method.postProcess(); CHKERRQ(ierr);
  PetscLogEventEnd(USER_EVENT_postProcess,0,0,0,0);
  PetscFunctionReturn(0);
}
PetscErrorCode Core::problem_basic_method_postProcess(const string &problem_name,BasicMethod &method,int verb) {
  PetscFunctionBegin;
  if(verb==-1) verb = verbose;
  typedef MoFEMProblem_multiIndex::index<Problem_mi_tag>::type moFEMProblems_by_name;
  // find p_miit
  moFEMProblems_by_name &moFEMProblems_set = moFEMProblems.get<Problem_mi_tag>();
  moFEMProblems_by_name::iterator p_miit = moFEMProblems_set.find(problem_name);
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
  typedef NumeredMoFEMFiniteElement_multiIndex::index<Composite_mi_tag>::type FEs_by_composite;
  method.feName = fe_name;
  SET_BASIC_METHOD(&*problem_ptr)
  PetscLogEventBegin(USER_EVENT_preProcess,0,0,0,0);
  ierr = method.preProcess(); CHKERRQ(ierr);
  PetscLogEventEnd(USER_EVENT_preProcess,0,0,0,0);
  FEs_by_composite &numeredFiniteElements = 
    (const_cast<NumeredMoFEMFiniteElement_multiIndex&>(problem_ptr->numeredFiniteElements)).get<Composite_mi_tag>();
  FEs_by_composite::iterator miit = numeredFiniteElements.lower_bound(boost::make_tuple(fe_name,lower_rank));
  FEs_by_composite::iterator hi_miit = numeredFiniteElements.upper_bound(boost::make_tuple(fe_name,upper_rank));
  for(;miit!=hi_miit;miit++) {
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
      SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INSONSISTENCY,ss.str().c_str());
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
  typedef MoFEMProblem_multiIndex::index<Problem_mi_tag>::type moFEMProblems_by_name;
  // find p_miit
  moFEMProblems_by_name &moFEMProblems_set = moFEMProblems.get<Problem_mi_tag>();
  moFEMProblems_by_name::iterator p_miit = moFEMProblems_set.find(problem_name);
  if(p_miit == moFEMProblems_set.end()) SETERRQ1(PETSC_COMM_SELF,1,"problem is not in database %s",problem_name.c_str());
  ierr = loop_finite_elements(&*p_miit,fe_name,method,lower_rank,upper_rank,verb); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
PetscErrorCode Core::loop_dofs(const MoFEMProblem *problem_ptr,const string &field_name,RowColData rc,EntMethod &method,int lower_rank,int upper_rank,int verb) {
  PetscFunctionBegin;
  SET_BASIC_METHOD(&*problem_ptr);
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
     SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INSONSISTENCY,"not implemented");
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
      SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INSONSISTENCY,ss.str().c_str());
    }
  }
  ierr = method.postProcess(); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
PetscErrorCode Core::loop_dofs(const string &problem_name,const string &field_name,RowColData rc,EntMethod &method,int lower_rank,int upper_rank,int verb) {
  PetscFunctionBegin;
  if(verb==-1) verb = verbose;
  typedef MoFEMProblem_multiIndex::index<Problem_mi_tag>::type moFEMProblems_by_name;
  // find p_miit
  moFEMProblems_by_name &moFEMProblems_set = moFEMProblems.get<Problem_mi_tag>();
  moFEMProblems_by_name::iterator p_miit = moFEMProblems_set.find(problem_name);
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
  SET_BASIC_METHOD(NULL);
  DofMoFEMEntity_multiIndex::index<FieldName_mi_tag>::type::iterator miit,hi_miit;
  miit = dofsMoabField.get<FieldName_mi_tag>().lower_bound(field_name);
  hi_miit = dofsMoabField.get<FieldName_mi_tag>().upper_bound(field_name);
  ierr = method.preProcess(); CHKERRQ(ierr);
  for(;miit!=hi_miit;miit++) {
    method.dofPtr = &*miit;
    method.dofNumeredPtr = NULL;
    try {
      ierr = method(); CHKERRQ(ierr);
    } catch (const std::exception& ex) {
      ostringstream ss;
      ss << "throw in method: " << ex.what() << endl;
      SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INSONSISTENCY,ss.str().c_str());
    }
  }
  ierr = method.postProcess(); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
PetscErrorCode Core::get_ref_ents(const RefMoFEMEntity_multiIndex **refinedEntities_ptr) {
  PetscFunctionBegin;
  *refinedEntities_ptr = &refinedEntities;
  PetscFunctionReturn(0);
}
PetscErrorCode Core::get_problem(const string &problem_name,const MoFEMProblem **problem_ptr) {
  PetscFunctionBegin;
  typedef MoFEMProblem_multiIndex::index<Problem_mi_tag>::type moFEMProblems_by_name;
  moFEMProblems_by_name &moFEMProblems_set = moFEMProblems.get<Problem_mi_tag>();
  moFEMProblems_by_name::iterator p_miit = moFEMProblems_set.find(problem_name);
  if(p_miit == moFEMProblems_set.end()) SETERRQ1(PETSC_COMM_SELF,1,"problem < %s > not found, (top tip: check spelling)",problem_name.c_str());
  *problem_ptr = &*p_miit;
  PetscFunctionReturn(0);
}
PetscErrorCode Core::get_dofs(const DofMoFEMEntity_multiIndex **dofsMoabField_ptr) {
  PetscFunctionBegin;
  *dofsMoabField_ptr = &dofsMoabField;
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
EntMoFEMFiniteElement_multiIndex::index<FiniteElement_name_mi_tag>::type::iterator Core::get_fes_moabfield_by_name_begin(const string &fe_name) {
  return finiteElementsMoFEMEnts.get<FiniteElement_name_mi_tag>().lower_bound(fe_name);
}
EntMoFEMFiniteElement_multiIndex::index<FiniteElement_name_mi_tag>::type::iterator Core::get_fes_moabfield_by_name_end(const string &fe_name) {
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
	SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INSONSISTENCY,"data inconsistency");
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
  } catch (const char* msg) {
    SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INSONSISTENCY,msg);
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
      EntMoFEMFiniteElement_multiIndex::index<Composite_Name_And_Ent_mi_tag>::type::iterator iit;
      iit = finiteElementsMoFEMEnts.get<Composite_Name_And_Ent_mi_tag>().find(
	boost::make_tuple(fe_it->get_name(),*eit));
      if(iit != finiteElementsMoFEMEnts.get<Composite_Name_And_Ent_mi_tag>().end()) {
	SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INSONSISTENCY,"data inconsistency");
      }
      eit++;
    }
    rval = moab.remove_entities(meshset,ents_to_remove); CHKERR_PETSC(rval);
    if(verb>0) {
      PetscPrintf(comm,
	"number of removed entities = %u from finite element %s\n",
	ents_to_remove.size(),fe_it->get_name().c_str());
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
	  refinedEntities.project<0>(pit),RefMoFEMEntity_change_remove_parent(moab));
	if(!success) {
	  SETERRQ(PETSC_COMM_SELF,MOFEM_OPERATION_UNSUCCESSFUL,"mofification unsucessfull");
	}
      }
    }
  }
  { //remove deleted entities form cubit meshsets
    CubitMeshSet_multiIndex::iterator cubit_it;
    cubit_it = cubit_meshsets.begin();
    for(;cubit_it!=cubit_meshsets.end();cubit_it++) {
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
