/** \file FieldCore.cpp
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

#include<FieldCore.hpp>
#include<FEM.h>
#include<version.h>

namespace MoFEM {

PetscErrorCode FieldCore::add_verices_in_the_middel_of_edges(const EntityHandle meshset,const BitRefLevel &bit,const bool recursive,int verb) {
  PetscFunctionBegin;
  if(verb==-1) verb = verbose;
  Range edges;
  rval = moab.get_entities_by_type(meshset,MBEDGE,edges,recursive);  CHKERR_PETSC(rval);
  if(edges.empty()) {
    Range tets;
    rval = moab.get_entities_by_type(meshset,MBTET,tets,recursive); CHKERR_PETSC(rval);
    rval = moab.get_adjacencies(tets,1,true,edges,Interface::UNION); CHKERR_PETSC(rval);
    if(tets.empty()) {
      Range prisms;
      rval = moab.get_entities_by_type(meshset,MBPRISM,prisms,recursive); CHKERR_PETSC(rval);
      for(Range::iterator pit = prisms.begin();pit!=prisms.end();pit++) {
        const EntityHandle* conn; 
        int num_nodes; 
        rval = moab.get_connectivity(*pit,conn,num_nodes,true);  CHKERR_PETSC(rval);
        assert(num_nodes==6);
        //
        Range edge;
        rval = moab.get_adjacencies(&conn[0],2,1,true,edge); CHKERR_PETSC(rval);
        assert(edge.size()==1);
        edges.insert(edge[0]);
        edge.clear();
        rval = moab.get_adjacencies(&conn[1],2,1,true,edge); CHKERR_PETSC(rval);
        assert(edge.size()==1);
        edges.insert(edge[0]);
        EntityHandle conn_edge2[] = { conn[2], conn[0] };
        edge.clear();
        rval = moab.get_adjacencies(conn_edge2,2,1,true,edge); CHKERR_PETSC(rval);
        assert(edge.size()==1);
        edges.insert(edge[0]);
        //
        edge.clear();
        rval = moab.get_adjacencies(&conn[3],2,1,true,edge); CHKERR_PETSC(rval);
        assert(edge.size()==1);
        edges.insert(edge[0]);
        edge.clear();
        rval = moab.get_adjacencies(&conn[4],2,1,true,edge); CHKERR_PETSC(rval);
        assert(edge.size()==1);
        edges.insert(edge[0]);
        EntityHandle conn_edge8[] = { conn[5], conn[3] };
        edge.clear();
        rval = moab.get_adjacencies(conn_edge8,2,1,true,edge); CHKERR_PETSC(rval);
        assert(edge.size()==1);
        edges.insert(edge[0]);
      }
    }
  }
  ierr = add_verices_in_the_middel_of_edges(edges,bit,verb); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
PetscErrorCode FieldCore::add_verices_in_the_middel_of_edges(const Range &_edges,const BitRefLevel &bit,int verb) {
  PetscFunctionBegin;
  Range edges = _edges;
  typedef RefMoFEMEntity_multiIndex::index<Composite_EntityType_And_ParentEntityType_mi_tag>::type ref_ents_by_composite;
  ref_ents_by_composite &ref_ents = refinedMoFemEntities.get<Composite_EntityType_And_ParentEntityType_mi_tag>();
  ref_ents_by_composite::iterator miit = ref_ents.lower_bound(boost::make_tuple(MBVERTEX,MBEDGE));
  ref_ents_by_composite::iterator hi_miit = ref_ents.upper_bound(boost::make_tuple(MBVERTEX,MBEDGE));
  RefMoFEMEntity_multiIndex_view_by_parent_entity ref_parent_ents_view;
  for(;miit!=hi_miit;miit++) {
    pair<RefMoFEMEntity_multiIndex_view_by_parent_entity::iterator,bool> p_ref_ent_view;
    p_ref_ent_view = ref_parent_ents_view.insert(&*miit);
    if(!p_ref_ent_view.second) {
      SETERRQ(PETSC_COMM_SELF,1,"non uniqe insertion");
    }
  }
  if(verb > 0) {
    ostringstream ss;
    ss << "ref level " << bit << " nb. edges to refine " << edges.size() << endl;
    PetscPrintf(PETSC_COMM_WORLD,ss.str().c_str());
  }
  Range::iterator eit = edges.begin();
  for(;eit!=edges.end();eit++) {
    RefMoFEMEntity_multiIndex_view_by_parent_entity::iterator miit_view = ref_parent_ents_view.find(*eit);
    const EntityHandle* conn; 
    int num_nodes; 
    rval = moab.get_connectivity(*eit,conn,num_nodes,true);  CHKERR_PETSC(rval);
    assert(num_nodes==2);
    if(miit_view == ref_parent_ents_view.end()||ref_parent_ents_view.empty()) {
      double coords[num_nodes*3]; 
      rval = moab.get_coords(conn,num_nodes,coords); CHKERR_PETSC(rval);
      cblas_daxpy(3,1.,&coords[3],1,coords,1);
      cblas_dscal(3,0.5,coords,1);
      EntityHandle node;
      rval = moab.create_vertex(coords,node); CHKERR_PETSC(rval);
      rval = moab.tag_set_data(th_RefParentHandle,&node,1,&*eit); CHKERR_PETSC(rval);
      rval = moab.tag_set_data(th_RefBitLevel,&node,1,&bit); CHKERR_PETSC(rval);
      pair<RefMoFEMEntity_multiIndex::iterator,bool> p_ent = refinedMoFemEntities.insert(RefMoFEMEntity(moab,node));
      if(!p_ent.second) SETERRQ(PETSC_COMM_SELF,1,"this entity is there");
      if(verbose>2) {
	ostringstream ss;
	ss << *(p_ent.first) << endl;
	PetscPrintf(PETSC_COMM_WORLD,ss.str().c_str());
      }
    } else {
      const EntityHandle node = (*miit_view)->get_ref_ent();
      bool success = refinedMoFemEntities.modify(refinedMoFemEntities.get<MoABEnt_mi_tag>().find(node),RefMoFEMEntity_change_add_bit(bit));
      if(!success) SETERRQ(PETSC_COMM_SELF,1,"inconsistency in data");
    }
  }
  PetscFunctionReturn(0);
}
PetscErrorCode FieldCore::refine_TET(const EntityHandle meshset,const BitRefLevel &bit,const bool respect_interface) {
  PetscFunctionBegin;
  Range tets;
  rval = moab.get_entities_by_type(meshset,MBTET,tets,false); CHKERR_PETSC(rval);
  ierr = refine_TET(tets,bit,respect_interface); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
PetscErrorCode FieldCore::refine_TET(const Range &_tets,const BitRefLevel &bit,const bool respect_interface) {
  PetscFunctionBegin;
  //FIXME: refinment is based on entity handlers, should work on global ids of nodes, this will allow parallelize agortihm in the future
  PetscFunctionBegin;
  typedef RefMoFEMEntity_multiIndex::index<MoABEnt_mi_tag>::type ref_ents_by_ent;
  ref_ents_by_ent &ref_ents_ent = refinedMoFemEntities.get<MoABEnt_mi_tag>();
  // find all verices which parent is edge
  typedef RefMoFEMEntity_multiIndex::index<Composite_EntityType_And_ParentEntityType_mi_tag>::type ref_ents_by_composite;
  ref_ents_by_composite &ref_ents = refinedMoFemEntities.get<Composite_EntityType_And_ParentEntityType_mi_tag>();
  ref_ents_by_composite::iterator miit = ref_ents.lower_bound(boost::make_tuple(MBVERTEX,MBEDGE));
  ref_ents_by_composite::iterator hi_miit = ref_ents.upper_bound(boost::make_tuple(MBVERTEX,MBEDGE));
  RefMoFEMEntity_multiIndex_view_by_parent_entity ref_parent_ents_view;
  for(;miit!=hi_miit;miit++) {
    pair<RefMoFEMEntity_multiIndex_view_by_parent_entity::iterator,bool> p_ref_ent_view;
    p_ref_ent_view = ref_parent_ents_view.insert(&*miit);
    if(!p_ref_ent_view.second) {
      SETERRQ(PETSC_COMM_SELF,1,"non uniqe insertion");
    }
  }
  typedef RefMoFEMElement_multiIndex::index<MoABEnt_mi_tag>::type ref_MoFEMFiniteElement_by_ent;
  ref_MoFEMFiniteElement_by_ent &ref_MoFEMFiniteElement = refinedMoFemElements.get<MoABEnt_mi_tag>();
  typedef RefMoFEMElement_multiIndex::index<Composite_EntType_mi_tag_and_ParentEntType_mi_tag>::type ref_ent_by_parent;
  ref_ent_by_parent &by_parent = refinedMoFemElements.get<Composite_EntType_mi_tag_and_ParentEntType_mi_tag>();
  typedef RefMoFEMElement_multiIndex::index<Composite_of_ParentEnt_And_BitsOfRefinedEdges_mi_tag>::type ref_ent_by_composite;
  ref_ent_by_composite &by_composite = refinedMoFemElements.get<Composite_of_ParentEnt_And_BitsOfRefinedEdges_mi_tag>();
  //
  if(respect_interface) {
    SETERRQ(PETSC_COMM_SELF,1,"not implemented, set last parameter in refine_TET to false");
  }
  //
  Range tets = _tets.subset_by_type(MBTET);
  Range::iterator tit = tets.begin();
  for(;tit!=tets.end();tit++) {
    ref_MoFEMFiniteElement_by_ent::iterator miit2 = ref_MoFEMFiniteElement.find(*tit);
    if(miit2==ref_MoFEMFiniteElement.end()) SETERRQ(PETSC_COMM_SELF,1,"this tet is not in refinedMoFemElements");
    //connectivity
    const EntityHandle* conn; 
    int num_nodes; 
    moab.get_connectivity(*tit,conn,num_nodes,true); 
    assert(num_nodes==4);
    for(int nn = 0;nn<num_nodes;nn++) {
      bool success = refinedMoFemEntities.modify(refinedMoFemEntities.get<MoABEnt_mi_tag>().find(conn[nn]),RefMoFEMEntity_change_add_bit(bit));
      if(!success) SETERRQ(PETSC_COMM_SELF,1,"can not set refinement bit level to tet node");
    }
    //get edges
    BitRefEdges parent_edges_bit(0);
    EntityHandle edge_new_nodes[6];
    fill(&edge_new_nodes[0],&edge_new_nodes[6],no_handle); 
    int split_edges[6];  
    fill(&split_edges[0],&split_edges[6],-1); 
    //hash map of nodes (RefMoFEMEntity) by edges (EntityHandle)
    map<EntityHandle /*edge*/,const RefMoFEMEntity* /*node*/> map_ref_nodes_by_edges; 
    for(int ee = 0;ee<6;ee++) { 
      EntityHandle edge = no_handle;
      rval = moab.side_element(*tit,1,ee,edge);  CHKERR_PETSC(rval);
      RefMoFEMEntity_multiIndex_view_by_parent_entity::iterator miit_view;
      miit_view = ref_parent_ents_view.find(edge);
      if(miit_view != ref_parent_ents_view.end()) {
	if(((*miit_view)->get_BitRefLevel()&bit).any()) {
	  edge_new_nodes[ee] = (*miit_view)->get_ref_ent(); 
	  map_ref_nodes_by_edges[(*miit_view)->get_parent_ent()] = &**miit_view;
	  split_edges[parent_edges_bit.count()] = ee;
	  parent_edges_bit.set(ee,1);
	}
      }
    }
    // swap nodes forward
    EntityHandle _conn_[4];
    copy(&conn[0],&conn[4],&_conn_[0]);
    // build connectivity for rf tets
    EntityHandle new_tets_conns[8*4];
    fill(&new_tets_conns[0],&new_tets_conns[8*4],no_handle);
    int sub_type = -1,nb_new_tets = 0;
    switch (parent_edges_bit.count()) {
      case 0: {
	  ref_ents_by_ent::iterator tit_miit;
	  tit_miit = ref_ents_ent.find(*tit);
	  if(tit_miit==ref_ents_ent.end()) SETERRQ(PETSC_COMM_SELF,1,"data inconsistency");
	  bool success = refinedMoFemEntities.modify(tit_miit,RefMoFEMEntity_change_add_bit(bit));
	  if(!success) SETERRQ(PETSC_COMM_SELF,1,"impossible tet");
	  Range tit_conn;
	  rval = moab.get_connectivity(&*tit,1,tit_conn,true); CHKERR_PETSC(rval);
	  for(Range::iterator nit = tit_conn.begin();nit!=tit_conn.end();nit++) {
	    ref_ents_by_ent::iterator nit_miit = ref_ents_ent.find(*nit);
	    if(nit_miit==ref_ents_ent.end()) SETERRQ(PETSC_COMM_SELF,1,"can not find face in refinedMoFemEntities");
	    bool success = refinedMoFemEntities.modify(nit_miit,RefMoFEMEntity_change_add_bit(bit));
	    if(!success) SETERRQ(PETSC_COMM_SELF,1,"impossible node");
	  }
	  Range tit_edges;
	  rval = moab.get_adjacencies(&*tit,1,1,false,tit_edges); CHKERR_PETSC(rval);
	  for(Range::iterator eit = tit_edges.begin();eit!=tit_edges.end();eit++) {
	    ref_ents_by_ent::iterator eit_miit = ref_ents_ent.find(*eit);
	    if(eit_miit==ref_ents_ent.end()) SETERRQ(PETSC_COMM_SELF,1,"can not find face in refinedMoFemEntities");
	    bool success = refinedMoFemEntities.modify(eit_miit,RefMoFEMEntity_change_add_bit(bit));
	    if(!success) SETERRQ(PETSC_COMM_SELF,1,"impossible edge");
	  }
	  Range tit_faces;
	  rval = moab.get_adjacencies(&*tit,1,2,false,tit_faces); CHKERR_PETSC(rval);
	  if(tit_faces.size()!=4) SETERRQ(PETSC_COMM_SELF,1,"existing tet in mofem database should have 4 adjacent edges");
	  for(Range::iterator fit = tit_faces.begin();fit!=tit_faces.end();fit++) {
	    ref_ents_by_ent::iterator fit_miit = ref_ents_ent.find(*fit);
	    if(fit_miit==ref_ents_ent.end()) SETERRQ(PETSC_COMM_SELF,1,"can not find face in refinedMoFemEntities");
	    bool success = refinedMoFemEntities.modify(fit_miit,RefMoFEMEntity_change_add_bit(bit));
	    if(!success) SETERRQ(PETSC_COMM_SELF,1,"impossible face");
	  }
	  continue;
	}
	break;
      case 1:
	sub_type = 0;
	tet_type_1(_conn_,split_edges[0],edge_new_nodes[split_edges[0]],new_tets_conns);
	nb_new_tets = 2;
	break;
      case 2:
	sub_type = tet_type_2(_conn_,split_edges,edge_new_nodes,new_tets_conns);
	if(sub_type&(4|8|16)) {
	  nb_new_tets = 3;
	  break;
	} else if(sub_type == 1) {
	  nb_new_tets = 4;
	  break; 
	};
	assert(0);
	break;
      case 3:
	sub_type = tet_type_3(_conn_,split_edges,edge_new_nodes,new_tets_conns);
	if(sub_type <= 4 ) {
	  nb_new_tets = 4;
	  break;
	} else if(sub_type <= 7 ) {
	  nb_new_tets = 5;
	  break;
	}
	assert(0);
      case 4:
	sub_type = tet_type_4(_conn_,split_edges,edge_new_nodes,new_tets_conns);
	if(sub_type == 0) {
	  nb_new_tets = 5;
	  break;
	} else if(sub_type <= 7) {
	  nb_new_tets = 6;
	  break;
	}  
	assert(0);
      case 5:
	sub_type = tet_type_5(moab,_conn_,edge_new_nodes,new_tets_conns);
	nb_new_tets = 7;
	break;
      case 6:
	sub_type = 0;
	tet_type_6(moab,_conn_,edge_new_nodes,new_tets_conns);
	nb_new_tets = 8;
	break;
      default:
	assert(0);
    }
    // find that tets
    bitset<8> ref_tets_bit(0);
    ref_ent_by_composite::iterator miit_composite = by_composite.lower_bound(boost::make_tuple(*tit,parent_edges_bit.to_ulong()));
    ref_ent_by_composite::iterator hi_miit_composite = by_composite.upper_bound(boost::make_tuple(*tit,parent_edges_bit.to_ulong()));
    ref_ent_by_composite::iterator miit_composite2 = miit_composite;
    for(int tt = 0;miit_composite2!=hi_miit_composite;miit_composite2++,tt++) {
      //add this tet if exist to this ref level
      EntityHandle tet = miit_composite2->get_ref_ent();
      refinedMoFemEntities.modify(refinedMoFemEntities.find(tet),RefMoFEMEntity_change_add_bit(bit));
      //set bit that this element is in databse - no need to create it
      ref_tets_bit.set(tt,1);
      if(verbose>2) {
	ostringstream ss;
	ss << miit_composite2->get_RefMoFEMElement() << endl;
	PetscPrintf(PETSC_COMM_WORLD,ss.str().c_str());
      }
    }
    if(miit_composite!=hi_miit_composite) {
      //if that tet has the same pattern of splitted edges it has to have the same number of refined 
      //children elements - if not thorw an error
      if(ref_tets_bit.count()!=(unsigned int)nb_new_tets) SETERRQ(PETSC_COMM_SELF,1,"data inconsistency");
    } else {
      //if this element was not refined or was refined with diffrent patterns of splitted edges create new elements
      EntityHandle ref_tets[8];
      for(int tt = 0;tt<nb_new_tets;tt++) {
	if(!ref_tets_bit.test(tt)) {
	  rval = moab.create_element(MBTET,&new_tets_conns[4*tt],4,ref_tets[tt]); CHKERR_PETSC(rval);
	  /*double coords[12];
	  ierr = moab.get_coords(&new_tets_conns[4*tt],4,coords); CHKERRQ(ierr);
	  double V = Shape_intVolumeMBTET(diffN_TET,coords); 
	  if(V<=0) {
	    ParallelComm* pcomm = ParallelComm::get_pcomm(&moab,MYPCOMM_INDEX);
	    if(pcomm->rank()==0) {
	      EntityHandle meshset_error_out;
	      rval = moab.create_meshset(MESHSET_SET,meshset_error_out); CHKERR_PETSC(rval);
	      rval = moab.add_entities(meshset_error_out,&*tit,1); CHKERR_PETSC(rval);
	      ierr = moab.write_file("error_out.vtk","VTK","",&meshset_error_out,1); CHKERRQ(ierr);
	    }
	    ierr = PetscBarrier(PETSC_NULL); CHKERRQ(ierr);
	    ostringstream ss;
	    ss << "tit " << new_tets_conns[4*tt] << "\n";
	    ss << coords[0] << " " << coords[1] << " " << coords[2] << "\n";
	    ss << coords[3] << " " << coords[4] << " " << coords[5] << "\n"; 
	    ss << coords[6] << " " << coords[7] << " " << coords[8] << "\n"; 
	    ss << coords[9] << " " << coords[10] << " " << coords[11] << "\n";
	    ss << "error tet saved to error_out.vtk"  << "\n";
	    SETERRQ(PETSC_COMM_SELF,1,ss.str().c_str());
	    assert(V>0); 
	  }*/
	  int ref_type[2];
	  ref_type[0] = parent_edges_bit.count();
	  ref_type[1] = sub_type; 
	  rval = moab.tag_set_data(th_RefType,&ref_tets[tt],1,ref_type); CHKERR_PETSC(rval);
	  rval = moab.tag_set_data(th_RefParentHandle,&ref_tets[tt],1,&*tit); CHKERR_PETSC(rval);
	  rval = moab.tag_set_data(th_RefBitLevel,&ref_tets[tt],1,&bit); CHKERR_PETSC(rval);
	  rval = moab.tag_set_data(th_RefBitEdge,&ref_tets[tt],1,&parent_edges_bit); CHKERR_PETSC(rval);
	  //add refined entity
	  pair<RefMoFEMEntity_multiIndex::iterator,bool> p_MoFEMEntity = refinedMoFemEntities.insert(RefMoFEMEntity(moab,ref_tets[tt]));
	  //add refined element
	  pair<RefMoFEMElement_multiIndex::iterator,bool> p_MoFEMFiniteElement;
	  try {
	    p_MoFEMFiniteElement = refinedMoFemElements.insert(ptrWrapperRefMoFEMElement(new RefMoFEMElement_TET(moab,&*p_MoFEMEntity.first)));
	  } catch (const char* msg) {
	    SETERRQ(PETSC_COMM_SELF,1,msg);
	  }
	  //set bit that this element is now in databse
	  ref_tets_bit.set(tt);
	  if(verbose>2) {
	    ostringstream ss;
	    ss << "add tet: " << *(p_MoFEMFiniteElement.first->get_RefMoFEMElement()) << endl;
	    PetscPrintf(PETSC_COMM_WORLD,ss.str().c_str());
	  }
	}
      }
      //find parents for new edges and faces
      //get tet edges and faces
      Range tit_edges,tit_faces;
      rval = moab.get_adjacencies(&*tit,1,1,false,tit_edges); CHKERR_PETSC(rval);
      rval = moab.get_adjacencies(&*tit,1,2,false,tit_faces); CHKERR_PETSC(rval);
      Range edges_nodes[6],faces_nodes[4];
      //for edges - add ref nodes
      //edges_nodes[ee] - contains all nodes on edge ee inluding mid nodes if exist
      Range::iterator eit = tit_edges.begin();
      for(int ee = 0;eit!=tit_edges.end();eit++,ee++) {
	rval = moab.get_connectivity(&*eit,1,edges_nodes[ee],true); CHKERR_PETSC(rval);
	map<EntityHandle,const RefMoFEMEntity*>::iterator map_miit = map_ref_nodes_by_edges.find(*eit);
	if(map_miit!=map_ref_nodes_by_edges.end()) {
	  edges_nodes[ee].insert(map_miit->second->get_ref_ent());
	}
      }
      //for faces - add ref nodes
      //faces_nodes[ff] - contains all nodes on face ff inluding mid nodes if exist
      Range::iterator fit=tit_faces.begin();
      for(int ff = 0;fit!=tit_faces.end();fit++,ff++) {
	rval = moab.get_connectivity(&*fit,1,faces_nodes[ff],true); CHKERR_PETSC(rval);
	Range fit_edges;
	rval = moab.get_adjacencies(&*fit,1,1,false,fit_edges); CHKERR_PETSC(rval);
	for(Range::iterator eit2 =  fit_edges.begin();eit2 != fit_edges.end();eit2++) {
	  map<EntityHandle,const RefMoFEMEntity*>::iterator map_miit = map_ref_nodes_by_edges.find(*eit2);
	  if(map_miit!=map_ref_nodes_by_edges.end()) {
	    faces_nodes[ff].insert(map_miit->second->get_ref_ent());
	  }
	}
      }
      //add ref nodes to tet
      //tet_nodes contains all nodes on tet inluding mid edge nodes
      Range tet_nodes;
      rval = moab.get_connectivity(&*tit,1,tet_nodes,true); CHKERR_PETSC(rval);
      for(map<EntityHandle,const RefMoFEMEntity*>::iterator map_miit = map_ref_nodes_by_edges.begin();
	map_miit != map_ref_nodes_by_edges.end();map_miit++) {
	tet_nodes.insert(map_miit->second->get_ref_ent());
      }
      Range ref_edges;
      //get all all edges of refined tets
      rval = moab.get_adjacencies(ref_tets,nb_new_tets,1,true,ref_edges,Interface::UNION); CHKERR_PETSC(rval);
      //check for all ref edge and set parents
      for(Range::iterator reit = ref_edges.begin();reit!=ref_edges.end();reit++) {
	Range ref_edges_nodes;
	rval = moab.get_connectivity(&*reit,1,ref_edges_nodes,true); CHKERR_PETSC(rval);
	//check if ref edge is an coarse edge
	int ee = 0;
	for(;ee<6;ee++) {
	  //two nodes are common (node[0],node[1],ref_node (if exist))
	  //this tests if given edge is contained by edge of refined tetrahedral
	  if(intersect(edges_nodes[ee],ref_edges_nodes).size()==2) {
	    EntityHandle edge = tit_edges[ee];
	    rval = moab.tag_set_data(th_RefParentHandle,&*reit,1,&edge); CHKERR_PETSC(rval);
	    pair<RefMoFEMEntity_multiIndex::iterator,bool> p_ent = refinedMoFemEntities.insert(RefMoFEMEntity(moab,*reit));
	    bool success = refinedMoFemEntities.modify(p_ent.first,RefMoFEMEntity_change_add_bit(bit));
	    if(!success) SETERRQ(PETSC_COMM_SELF,1,"impossible to set edge pranet");
	    if(p_ent.second) {
	      if(verbose>2) {
		ostringstream ss;
		ss << "edge parent: " << *(p_ent.first) << endl;
		PetscPrintf(PETSC_COMM_WORLD,ss.str().c_str());
	    }}
	    break;
	  }  
	}
	if(ee<6) continue; //this refined edge is contined by edge of tetrahedral
	//check if ref edge is in coarse face
	int ff = 0;
	for(;ff<4;ff++) {
	  //two nodes are common (node[0],node[1],ref_node (if exist))
	  //thi tests if givem edge is contained by face of  tetrahedral
	  if(intersect(faces_nodes[ff],ref_edges_nodes).size()==2) {
	    EntityHandle face = tit_faces[ff];
	    rval = moab.tag_set_data(th_RefParentHandle,&*reit,1,&face); CHKERR_PETSC(rval);
	    //add edge to refinedMoFemEntities
	    pair<RefMoFEMEntity_multiIndex::iterator,bool> p_ent = refinedMoFemEntities.insert(RefMoFEMEntity(moab,*reit));
	    bool success = refinedMoFemEntities.modify(p_ent.first,RefMoFEMEntity_change_add_bit(bit));
	    if(!success) SETERRQ(PETSC_COMM_SELF,1,"impossible to set edge parent");
	    if(p_ent.second) {
	      if(verbose>2) {
		ostringstream ss;
		ss << "face parent: " << *(p_ent.first) << endl;
		PetscPrintf(PETSC_COMM_WORLD,ss.str().c_str());
	    }}
	    break;
	  }
	}
	if(ff<4) continue; //this refined egde is contained by face of tetrahedral
	// check if ref edge is in coarse tetrahedral (i.e. that is internal edge of refined tetrahedral)
	if(intersect(tet_nodes,ref_edges_nodes).size()==2) {
	  rval = moab.tag_set_data(th_RefParentHandle,&*reit,1,&*tit); CHKERR_PETSC(rval);
	  //add edge to refinedMoFemEntities
	  pair<RefMoFEMEntity_multiIndex::iterator,bool> p_ent = refinedMoFemEntities.insert(RefMoFEMEntity(moab,*reit));
	  bool success = refinedMoFemEntities.modify(p_ent.first,RefMoFEMEntity_change_add_bit(bit));
	  if(!success) SETERRQ(PETSC_COMM_SELF,1,"impossible to set edge parent");
	  if(p_ent.second) {
	    if(verbose>2) {
	      ostringstream ss;
	      ss << "tet parent: " << *(p_ent.first) << endl;
	      PetscPrintf(PETSC_COMM_WORLD,ss.str().c_str());
	  }}
	  continue;
	}
	//refined edge is not child of any edge, face or tetrahedral, this is imposible edge
	SETERRQ(PETSC_COMM_SELF,1,"impossible refined edge");
      }
      Range ref_faces;
      rval = moab.get_adjacencies(ref_tets,nb_new_tets,2,true,ref_faces,Interface::UNION); CHKERR_PETSC(rval);
      Tag th_interface_side;
      const int def_side[] = {0};
      rval = moab.tag_get_handle("INTERFACE_SIDE",1,MB_TYPE_INTEGER,
	th_interface_side,MB_TAG_CREAT|MB_TAG_SPARSE,def_side); CHKERR_PETSC(rval);
      // check for all ref faces
      for(Range::iterator rfit = ref_faces.begin();rfit!=ref_faces.end();rfit++) {
	Range ref_faces_nodes;
	rval = moab.get_connectivity(&*rfit,1,ref_faces_nodes,true); CHKERR_PETSC(rval);
	// check if ref face is in coarse face
	int ff = 0;
	for(;ff<4;ff++) {
	  //check if refined edge is contained by face of tetrahedral
	  if(intersect(faces_nodes[ff],ref_faces_nodes).size()==3) {
	    EntityHandle face = tit_faces[ff];
	    rval = moab.tag_set_data(th_RefParentHandle,&*rfit,1,&face); CHKERR_PETSC(rval);
	    int side = 0;
	    //set face side if it is on inteface
	    rval = moab.tag_get_data(th_interface_side,&face,1,&side); CHKERR_PETSC(rval);
	    rval = moab.tag_set_data(th_interface_side,&*rfit,1,&side); CHKERR_PETSC(rval);
	    //add face to refinedMoFemEntities
	    pair<RefMoFEMEntity_multiIndex::iterator,bool> p_ent = refinedMoFemEntities.insert(RefMoFEMEntity(moab,*rfit));
	    bool success = refinedMoFemEntities.modify(p_ent.first,RefMoFEMEntity_change_add_bit(bit));
	    if(!success) SETERRQ(PETSC_COMM_SELF,1,"impossible to set face parent");
	    if(p_ent.second) {
	      if(verbose>2) {
		ostringstream ss;
		ss << "face parent: " << *(p_ent.first) << endl;
		PetscPrintf(PETSC_COMM_WORLD,ss.str().c_str());
	    }}
	    break;
	  }
	}
	if(ff<4) continue; //this face is contained by one of tetrahedrals 
	//check if ref face is in coarse tetrahedral
	//this is ref face which is contained by tetrahedral volume
	if(intersect(tet_nodes,ref_faces_nodes).size()==3) {
	  rval = moab.tag_set_data(th_RefParentHandle,&*rfit,1,&*tit); CHKERR_PETSC(rval);
	  pair<RefMoFEMEntity_multiIndex::iterator,bool> p_ent = refinedMoFemEntities.insert(RefMoFEMEntity(moab,*rfit));
	  //add face to refinedMoFemEntities
	  bool success = refinedMoFemEntities.modify(p_ent.first,RefMoFEMEntity_change_add_bit(bit));
	  if(!success) SETERRQ(PETSC_COMM_SELF,1,"impossible to set face parent");
	  if(p_ent.second) {
	    if(verbose>2) {
	      ostringstream ss;
	      ss << "tet parent: " << *(p_ent.first) << endl;
	      PetscPrintf(PETSC_COMM_WORLD,ss.str().c_str());
	  }}
	  continue;
	}
	SETERRQ(PETSC_COMM_SELF,1,"impossible refined face");
      }
    }
  }
  PetscFunctionReturn(0);
}
PetscErrorCode FieldCore::refine_PRISM(const EntityHandle meshset,const BitRefLevel &bit,int verb) {
  //FIXME: refinment is based on entity handlers, should work on global ids of nodes, this will allow parallelize agortihm in the future
  PetscFunctionBegin;
  if(verb==-1) verb = verbose;
  typedef RefMoFEMEntity_multiIndex::index<MoABEnt_mi_tag>::type ref_ENTs_by_ent;
  typedef RefMoFEMElement_multiIndex::index<Composite_of_ParentEnt_And_BitsOfRefinedEdges_mi_tag>::type ref_fe_by_composite;
  ref_fe_by_composite &ref_fe_by_comp = refinedMoFemElements.get<Composite_of_ParentEnt_And_BitsOfRefinedEdges_mi_tag>();
  //find all verices which parent is edge
  typedef RefMoFEMEntity_multiIndex::index<Composite_EntityType_And_ParentEntityType_mi_tag>::type ref_ents_by_composite;
  ref_ents_by_composite &ref_ents_by_comp = refinedMoFemEntities.get<Composite_EntityType_And_ParentEntityType_mi_tag>();
  ref_ents_by_composite::iterator miit = ref_ents_by_comp.lower_bound(boost::make_tuple(MBVERTEX,MBEDGE));
  ref_ents_by_composite::iterator hi_miit = ref_ents_by_comp.upper_bound(boost::make_tuple(MBVERTEX,MBEDGE));
  RefMoFEMEntity_multiIndex_view_by_parent_entity ref_parent_ents_view;
  for(;miit!=hi_miit;miit++) {
    pair<RefMoFEMEntity_multiIndex_view_by_parent_entity::iterator,bool> p_ref_ent_view;
    p_ref_ent_view = ref_parent_ents_view.insert(&*miit);
    if(!p_ref_ent_view.second) {
      SETERRQ(PETSC_COMM_SELF,1,"non uniqe insertion");
    }
  }
  Range prisms;
  rval = moab.get_entities_by_type(meshset,MBPRISM,prisms,false); CHKERR_PETSC(rval);
  Range::iterator pit = prisms.begin();
  for(;pit!=prisms.end();pit++) {
    ref_ENTs_by_ent::iterator miit_prism = refinedMoFemEntities.get<MoABEnt_mi_tag>().find(*pit);   
    if(miit_prism==refinedMoFemEntities.end()) SETERRQ(PETSC_COMM_SELF,1,"this prism is not in ref database");
    if(verb>3) {
      ostringstream ss;
      ss << "ref prism " << *miit << endl;
      PetscPrintf(PETSC_COMM_WORLD,ss.str().c_str());
    }
    //prism connectivity
    int num_nodes; 
    const EntityHandle* conn;
    rval = moab.get_connectivity(*pit,conn,num_nodes,true); CHKERR_PETSC(rval);
    assert(num_nodes==6);
    // edges connectivity
    EntityHandle edges[6];
    for(int ee = 0;ee<3; ee++) {
      rval = moab.side_element(*pit,1,ee,edges[ee]); CHKERR_PETSC(rval);
    }
    for(int ee = 6;ee<9; ee++) {
      rval = moab.side_element(*pit,1,ee,edges[ee-3]); CHKERR_PETSC(rval);
    }
    // detetct split edges
    BitRefEdges split_edges(0);
    EntityHandle edge_nodes[6];
    fill(&edge_nodes[0],&edge_nodes[6],no_handle);
    for(int ee = 0;ee<6;ee++) {
      RefMoFEMEntity_multiIndex_view_by_parent_entity::iterator miit_view = ref_parent_ents_view.find(edges[ee]);
      if(miit_view != ref_parent_ents_view.end()) {
	if(((*miit_view)->get_BitRefLevel()&bit).any()) {
	  edge_nodes[ee] = (*miit_view)->get_ref_ent(); 
	  split_edges.set(ee);
	}
      }
    }
    if(split_edges.count()==0) {
      refinedMoFemEntities.modify(miit_prism,RefMoFEMEntity_change_add_bit(bit));
      if(verb>6) PetscPrintf(PETSC_COMM_WORLD,"no refinement");
      continue;
    } 
    //check consitency
    if(verb>3) {
      ostringstream ss;
      ss << "prism split edges " << split_edges << " count " << split_edges.count() << endl;
      PetscPrintf(PETSC_COMM_WORLD,ss.str().c_str());
    }
    // prism ref
    EntityHandle new_prism_conn[4*6];
    fill(&new_prism_conn[0],&new_prism_conn[4*6],no_handle);
    int nb_new_prisms = 0;
    switch (split_edges.count()) {
      case 0:
	break;
      case 2:
	ierr = prism_type_1(conn,split_edges,edge_nodes,new_prism_conn); CHKERRQ(ierr);
	nb_new_prisms = 2;
	break;
      case 4:
	ierr = prism_type_2(conn,split_edges,edge_nodes,new_prism_conn); CHKERRQ(ierr);
	nb_new_prisms = 3;
	break;
      case 6:
	ierr = prism_type_3(conn,split_edges,edge_nodes,new_prism_conn); CHKERRQ(ierr);
	nb_new_prisms = 4;
	break;
      default:
	ostringstream ss;
	ss << split_edges << " : [ " 
	  << conn[0] << " "
	  << conn[1] << " "
	  << conn[2] << " "
	  << conn[3] << " "
	  << conn[4] << " "
	  << conn[5] << " ]";
	SETERRQ(PETSC_COMM_SELF,1,ss.str().c_str());
    }
    // find that prism
    bitset<4> ref_prism_bit(0);
    ref_fe_by_composite::iterator miit_composite = ref_fe_by_comp.lower_bound(boost::make_tuple(*pit,split_edges.to_ulong()));
    ref_fe_by_composite::iterator hi_miit_composite = ref_fe_by_comp.upper_bound(boost::make_tuple(*pit,split_edges.to_ulong()));
    ref_fe_by_composite::iterator miit_composite2 = miit_composite;
    for(int pp = 0;miit_composite2!=hi_miit_composite;miit_composite2++,pp++) {
      //add this tet to this ref
      refinedMoFemEntities.modify(refinedMoFemEntities.find(miit_composite2->get_ref_ent()),RefMoFEMEntity_change_add_bit(bit));
      ref_prism_bit.set(pp,1);
      if(verb>2) {
	ostringstream ss;
	ss << "is refined " << *(miit_composite2->get_RefMoFEMElement()) << endl;
	PetscPrintf(PETSC_COMM_WORLD,ss.str().c_str());
      }
    }
    if(miit_composite!=hi_miit_composite) {
      if(ref_prism_bit.count()!=(unsigned int)nb_new_prisms) SETERRQ(PETSC_COMM_SELF,1,"data inconsistency");
    } else {
      EntityHandle ref_prisms[4];
      // create prism
      for(int pp = 0;pp<nb_new_prisms;pp++) {
	if(verb>3) {
	  ostringstream ss;
	  ss << "ref prism " << ref_prism_bit << endl;
	  PetscPrintf(PETSC_COMM_WORLD,ss.str().c_str());
	}
	if(!ref_prism_bit.test(pp)) {
	  rval = moab.create_element(MBPRISM,&new_prism_conn[6*pp],6,ref_prisms[pp]); CHKERR_PETSC(rval);
	  rval = moab.tag_set_data(th_RefParentHandle,&ref_prisms[pp],1,&*pit); CHKERR_PETSC(rval);
	  rval = moab.tag_set_data(th_RefBitLevel,&ref_prisms[pp],1,&bit); CHKERR_PETSC(rval);
	  rval = moab.tag_set_data(th_RefBitEdge,&ref_prisms[pp],1,&split_edges); CHKERR_PETSC(rval);
	  pair<RefMoFEMEntity_multiIndex::iterator,bool> p_ent = refinedMoFemEntities.insert(RefMoFEMEntity(moab,ref_prisms[pp]));
	  pair<RefMoFEMElement_multiIndex::iterator,bool> p_MoFEMFiniteElement;
	  try {
	    p_MoFEMFiniteElement = refinedMoFemElements.insert(ptrWrapperRefMoFEMElement(new RefMoFEMElement_PRISM(moab,&*p_ent.first)));
	  } catch (const char* msg) {
	    SETERRQ(PETSC_COMM_SELF,1,msg);
	  }
	  ref_prism_bit.set(pp);
	  ierr = add_prism_to_mofem_database(ref_prisms[pp]); CHKERRQ(ierr);
	  if(verb>2) {
	    ostringstream ss;
	    ss << "add prism: " << *(p_MoFEMFiniteElement.first->get_RefMoFEMElement()) << endl;
	    if(verb>7) {
	      for(int nn = 0;nn<6;nn++) {
		ss << new_prism_conn[nn] << " ";
	      }
	      ss << endl;
	    }
	    PetscPrintf(PETSC_COMM_WORLD,ss.str().c_str());
	  }
	}
      }
    }
  }
  PetscFunctionReturn(0);
}
PetscErrorCode FieldCore::refine_MESHSET(const EntityHandle meshset,const BitRefLevel &bit,const bool recursive,int verb) {
  PetscFunctionBegin;
  if(verb==-1) verb = verbose;
  typedef RefMoFEMEntity_multiIndex::index<MoABEnt_mi_tag>::type ref_ENTs_by_ent;
  ref_ENTs_by_ent::iterator miit = refinedMoFemEntities.find(meshset);
  if(miit==refinedMoFemEntities.end()) SETERRQ(PETSC_COMM_SELF,1,"this meshset is not in ref database");
  ierr = update_meshset_by_entities_children(meshset,bit,meshset,MBEDGE,recursive,verb); CHKERRQ(ierr);
  ierr = update_meshset_by_entities_children(meshset,bit,meshset,MBTRI,recursive,verb); CHKERRQ(ierr);
  ierr = update_meshset_by_entities_children(meshset,bit,meshset,MBTET,recursive,verb); CHKERRQ(ierr);
  refinedMoFemEntities.modify(miit,RefMoFEMEntity_change_add_bit(bit));
  PetscFunctionReturn(0);
}
PetscErrorCode FieldCore::get_msId_3dENTS_sides(const int msId,const Cubit_BC_bitset CubitBCType,const BitRefLevel mesh_bit_level,const bool recursive,int verb) {
  PetscFunctionBegin;
  if(verb==-1) verb = verbose;
  moabCubitMeshSet_multiIndex::index<Composite_Cubit_msId_and_MeshSetType_mi_tag>::type::iterator 
    miit = cubit_meshsets.get<Composite_Cubit_msId_and_MeshSetType_mi_tag>().find(boost::make_tuple(msId,CubitBCType.to_ulong()));
  if(miit!=cubit_meshsets.get<Composite_Cubit_msId_and_MeshSetType_mi_tag>().end()) {
    ierr = FieldCore::get_msId_3dENTS_sides(miit->meshset,mesh_bit_level,recursive,verb); CHKERRQ(ierr);
  } else {
    SETERRQ(PETSC_COMM_SELF,1,"msId is not there");
  }
  PetscFunctionReturn(0);
}
PetscErrorCode FieldCore::get_msId_3dENTS_sides(const EntityHandle SideSet,const BitRefLevel mesh_bit_level,const bool recursive,int verb) {
  PetscFunctionBegin;
  if(verb==-1) verb = verbose;
  Range mesh_level_ents3d;
  Range mesh_level_tris;
  Range mesh_level_edges;
  Range mesh_level_nodes;
  if(mesh_bit_level.any()) {
    ierr = get_entities_by_type_and_ref_level(mesh_bit_level,BitRefLevel().set(),MBTET,mesh_level_ents3d); CHKERRQ(ierr);
    ierr = get_entities_by_type_and_ref_level(mesh_bit_level,BitRefLevel().set(),MBTRI,mesh_level_tris); CHKERRQ(ierr);
    ierr = get_entities_by_type_and_ref_level(mesh_bit_level,BitRefLevel().set(),MBEDGE,mesh_level_edges); CHKERRQ(ierr);
    ierr = get_entities_by_type_and_ref_level(mesh_bit_level,BitRefLevel().set(),MBVERTEX,mesh_level_nodes); CHKERRQ(ierr);
  }
  Range mesh_level_prisms;
  if(mesh_bit_level.any()) {
    ierr = get_entities_by_type_and_ref_level(mesh_bit_level,BitRefLevel().set(),MBPRISM,mesh_level_prisms); CHKERRQ(ierr);
    mesh_level_ents3d.merge(mesh_level_prisms);
  }
  Skinner skin(&moab);
  //get interface triangles from side set
  Range triangles;
  rval = moab.get_entities_by_type(SideSet,MBTRI,triangles,recursive);  CHKERR_PETSC(rval);
  if(mesh_bit_level.any()) {
    triangles = intersect(triangles,mesh_level_tris);
  }
  if(verb>1) {
    PetscPrintf(PETSC_COMM_WORLD,"Nb. of triangles in set %u\n",triangles.size());
  }
  //get nodes, edges and 3d ents (i.e. tets and prisms)
  Range nodes; // nodes from triangles
  rval = moab.get_connectivity(triangles,nodes,true); CHKERR_PETSC(rval);
  Range ents3d,ents3d_with_prisms; // 3d ents form nodes
  rval = moab.get_adjacencies(nodes,3,false,ents3d_with_prisms,Interface::UNION); CHKERR_PETSC(rval);
  if(mesh_bit_level.any()) {
    ents3d_with_prisms = intersect(ents3d_with_prisms,mesh_level_ents3d);
  }
  ents3d = ents3d_with_prisms.subset_by_type(MBTET); // take only tets, add prism later
  //take skin faces
  Range skin_faces; // skin faces from 3d ents
  rval = skin.find_skin(ents3d,false,skin_faces); CHKERR(rval);
  //take skin edges (boundary of surface if there is any)
  Range skin_edges_boundary; //skin edges from triangles
  rval = skin.find_skin(triangles,false,skin_edges_boundary); CHKERR(rval);
  if(verb>3) PetscPrintf(PETSC_COMM_WORLD,"skin_edges_boundary %u\n",skin_edges_boundary.size());
  //take all edges on skin faces (i.e. skin surface)
  Range skin_faces_edges; //edges from skin faces of 3d ents
  rval = moab.get_adjacencies(skin_faces,1,false,skin_faces_edges,Interface::UNION); CHKERR_PETSC(rval);
  if(mesh_bit_level.any()) {
    skin_faces_edges = intersect(skin_faces_edges,mesh_level_edges);
  }
  if(verb>3) PetscPrintf(PETSC_COMM_WORLD,"skin_faces_edges %u\n",skin_faces_edges.size());
  //note: that skin faces edges do not contain internal boundary
  //note: that prisms are not included in ents3d, so if ents3d have border with other inteface is like external boundary 
  //skin edges bondart are internal edge <- skin_faces_edges contains edges which are on the body boundary <- that is the trick
  skin_edges_boundary = subtract(skin_edges_boundary,skin_faces_edges); // from skin edges subtract edges from skin faces of 3d ents (only internal edges)
  if(verb>3) {
    EntityHandle out_meshset;
    rval = moab.create_meshset(MESHSET_SET,out_meshset); CHKERR_PETSC(rval);
    rval = moab.add_entities(out_meshset,triangles); CHKERR_PETSC(rval);
    rval = moab.write_file("triangles.vtk","VTK","",&out_meshset,1); CHKERR_PETSC(rval);
    rval = moab.delete_entities(&out_meshset,1); CHKERR_PETSC(rval);
    rval = moab.create_meshset(MESHSET_SET,out_meshset); CHKERR_PETSC(rval);
    rval = moab.add_entities(out_meshset,ents3d); CHKERR_PETSC(rval);
    rval = moab.write_file("ents3d.vtk","VTK","",&out_meshset,1); CHKERR_PETSC(rval);
    rval = moab.delete_entities(&out_meshset,1); CHKERR_PETSC(rval);
    rval = moab.create_meshset(MESHSET_SET,out_meshset); CHKERR_PETSC(rval);
    rval = moab.add_entities(out_meshset,skin_edges_boundary); CHKERR_PETSC(rval);
    rval = moab.write_file("skin_edges_boundary.vtk","VTK","",&out_meshset,1); CHKERR_PETSC(rval);
    rval = moab.delete_entities(&out_meshset,1); CHKERR_PETSC(rval);
  }
  if(verb>3) PetscPrintf(PETSC_COMM_WORLD,"subtract skin_edges_boundary %u\n",skin_edges_boundary.size());
  //Get nodes on boundary edge
  Range skin_nodes_boundary;
  rval = moab.get_connectivity(skin_edges_boundary,skin_nodes_boundary,true); CHKERR_PETSC(rval);
  //Remove noded which are bondary with other existing interface
  Range prisms_nodes;
  rval = moab.get_connectivity(ents3d_with_prisms.subset_by_type(MBPRISM),prisms_nodes,true); CHKERR_PETSC(rval);
  skin_nodes_boundary = subtract(skin_nodes_boundary,prisms_nodes);
  if(verb>3) PetscPrintf(PETSC_COMM_WORLD,"subtract skin_nodes_boundary %u\n",skin_nodes_boundary.size());
  //use nodes on body boundary and interface (without internal boundary) to find adjacent tets
  Range nodes_without_front = subtract(nodes,skin_nodes_boundary); // nodes_without_front adjacent to all splitted face edges except those on internal edge
  if(verb>3) PetscPrintf(PETSC_COMM_WORLD,"adj. node if ents3d but not on the internal edge %u\n",nodes_without_front.size());
  //ents3 that are adjacent to front nodes on splitted faces but not those which are on the front nodes on internal edgea
  ents3d.clear();
  ents3d_with_prisms.clear();
  rval = moab.get_adjacencies(nodes_without_front,3,false,ents3d_with_prisms,Interface::UNION); CHKERR_PETSC(rval);
  if(mesh_bit_level.any()) {
    ents3d_with_prisms = intersect(ents3d_with_prisms,mesh_level_ents3d);
  }
  ents3d = ents3d_with_prisms.subset_by_type(MBTET);
  if(verb>3) PetscPrintf(PETSC_COMM_WORLD,"adj. ents3d to fornt nodes %u\n",ents3d.size());
  Range side_ents3d;
  unsigned int nb_side_ents3d = side_ents3d.size();
  side_ents3d.insert(*ents3d.begin());
  do {
    Range adj_tris,adj_ents3d;
    nb_side_ents3d = side_ents3d.size();
    if(verb>2) PetscPrintf(PETSC_COMM_WORLD,"nb_side_ents3d %u\n",nb_side_ents3d);
    //get faces
    rval = moab.get_adjacencies(side_ents3d.subset_by_type(MBTET),2,false,adj_tris,Interface::UNION); CHKERR_PETSC(rval);
    if(mesh_bit_level.any()) {
      adj_tris = intersect(adj_tris,mesh_level_tris);
    }
    //subtrace from faces interface
    adj_tris = subtract(adj_tris,triangles);
    if(verb>2) PetscPrintf(PETSC_COMM_WORLD,"adj_tris %u\n",adj_tris.size());
    //get tets adjacent to faces
    rval = moab.get_adjacencies(adj_tris,3,true,adj_ents3d,Interface::UNION); CHKERR_PETSC(rval);
    //intersect tets with tets adjacent to inetface
    adj_ents3d = intersect(adj_ents3d,ents3d_with_prisms);
    if(verb>2) PetscPrintf(PETSC_COMM_WORLD,"adj_ents3d %u\n",adj_ents3d.size());
    //add tets to side
    side_ents3d.insert(adj_ents3d.begin(),adj_ents3d.end());
  } while (nb_side_ents3d != side_ents3d.size());
  //other side ents
  Range other_side = subtract(ents3d_with_prisms,side_ents3d);
  //side nodes
  Range side_nodes;
  rval = moab.get_connectivity(side_ents3d.subset_by_type(MBTET),side_nodes,true); CHKERR_PETSC(rval);
  //nodes on crack surface without front
  nodes_without_front = intersect(nodes_without_front,side_nodes);
  Range side_edges;
  rval = moab.get_adjacencies(side_ents3d.subset_by_type(MBTET),1,false,side_edges,Interface::UNION); CHKERR_PETSC(rval);
  skin_edges_boundary = intersect(skin_edges_boundary,side_edges);
  //make child meshsets
  vector<EntityHandle> children;
  rval = moab.get_child_meshsets(SideSet,children);  CHKERR_PETSC(rval);
  if(children.empty()) {
    children.resize(3);
    rval = moab.create_meshset(MESHSET_SET,children[0]); CHKERR_PETSC(rval);
    rval = moab.create_meshset(MESHSET_SET,children[1]); CHKERR_PETSC(rval);
    rval = moab.create_meshset(MESHSET_SET,children[2]); CHKERR_PETSC(rval);
    rval = moab.add_child_meshset(SideSet,children[0]); CHKERR_PETSC(rval);
    rval = moab.add_child_meshset(SideSet,children[1]); CHKERR_PETSC(rval);
    rval = moab.add_child_meshset(SideSet,children[2]); CHKERR_PETSC(rval);
  } else { 
    if(children.size()!=3) {
      SETERRQ(PETSC_COMM_SELF,1,"this meshset shuld have 3 children meshsets");
    }
    children.resize(3);
    ierr = moab.clear_meshset(&children[0],3); CHKERRQ(ierr);
  }
  EntityHandle &child_side = children[0];
  EntityHandle &child_other_side = children[1];
  EntityHandle &child_nodes_and_skin_edges = children[2];
  rval = moab.add_entities(child_side,side_ents3d); CHKERR_PETSC(rval);
  rval = moab.add_entities(child_other_side,other_side); CHKERR_PETSC(rval);
  rval = moab.add_entities(child_nodes_and_skin_edges,nodes_without_front); CHKERR_PETSC(rval);
  rval = moab.add_entities(child_nodes_and_skin_edges,skin_edges_boundary); CHKERR_PETSC(rval);
  if(verb>1) {
    PetscPrintf(PETSC_COMM_WORLD,"Nb. of side ents3d in set %u\n",side_ents3d.size());
    PetscPrintf(PETSC_COMM_WORLD,"Nb. of other side ents3d in set %u\n",other_side.size());
    PetscPrintf(PETSC_COMM_WORLD,"Nb. of boudary edges %u\n",skin_edges_boundary.size());
  }
  if(verb>3) {
    ierr = moab.write_file("side.vtk","VTK","",&children[0],1); CHKERRQ(ierr);
    ierr = moab.write_file("other_side.vtk","VTK","",&children[1],1); CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}
PetscErrorCode FieldCore::get_msId_3dENTS_split_sides(
  const EntityHandle meshset,const BitRefLevel &bit,
  const int msId,const Cubit_BC_bitset CubitBCType,const bool add_iterfece_entities,const bool recursive,int verb) {
  PetscFunctionBegin;
  if(verb==-1) verb = verbose;
  moabCubitMeshSet_multiIndex::index<Composite_Cubit_msId_and_MeshSetType_mi_tag>::type::iterator 
    miit = cubit_meshsets.get<Composite_Cubit_msId_and_MeshSetType_mi_tag>().find(boost::make_tuple(msId,CubitBCType.to_ulong()));
  if(miit!=cubit_meshsets.get<Composite_Cubit_msId_and_MeshSetType_mi_tag>().end()) {
    ierr = FieldCore::get_msId_3dENTS_split_sides(
      meshset,bit,miit->meshset,add_iterfece_entities,recursive,verb); CHKERRQ(ierr);
  } else {
    SETERRQ(PETSC_COMM_SELF,1,"msId is not there");
  }
  PetscFunctionReturn(0);
}
PetscErrorCode FieldCore::get_msId_3dENTS_split_sides(
  const EntityHandle meshset,const BitRefLevel &bit,
  const EntityHandle SideSet,const bool add_iterfece_entities,const bool recursive,
  int verb) {
  PetscFunctionBegin;
  ierr = get_msId_3dENTS_split_sides(meshset,bit,BitRefLevel(),SideSet,add_iterfece_entities,recursive,verb); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
PetscErrorCode FieldCore::get_msId_3dENTS_split_sides(
  const EntityHandle meshset,const BitRefLevel &bit,const BitRefLevel &inheret_from_bit_level,
  const EntityHandle SideSet,const bool add_iterfece_entities,const bool recursive,
  int verb) {
  PetscFunctionBegin;
  if(verb==-1) verb = verbose;
  vector<EntityHandle> children;
  //get children meshsets
  rval = moab.get_child_meshsets(SideSet,children);  CHKERR_PETSC(rval);
  if(children.size()!=3) {
    SETERRQ(PETSC_COMM_SELF,1,"should be 3 child meshsets, each of them contains tets on two sides of interface");
  }
  //faces of interface
  Range triangles;
  rval = moab.get_entities_by_type(SideSet,MBTRI,triangles,recursive);  CHKERR_PETSC(rval);
  //3d ents on "father" side
  Range side_ents3d;
  rval = moab.get_entities_by_handle(children[0],side_ents3d,false);  CHKERR_PETSC(rval);
  //3d ents on "mather" side
  Range other_ents3d;
  rval = moab.get_entities_by_handle(children[1],other_ents3d,false);  CHKERR_PETSC(rval);
  //nodes on interface but not on crack front (those should not be splitted)
  Range nodes;
  rval = moab.get_entities_by_type(children[2],MBVERTEX,nodes,false);  CHKERR_PETSC(rval);
  Range meshset_3d_ents,meshset_2d_ents;
  rval = moab.get_entities_by_dimension(meshset,3,meshset_3d_ents,true); CHKERR_PETSC(rval);
  Range meshset_tets = meshset_3d_ents.subset_by_type(MBTET);
  rval = moab.get_adjacencies(meshset_tets,2,false,meshset_2d_ents,moab::Interface::UNION); CHKERR_PETSC(rval);
  side_ents3d = intersect(meshset_3d_ents,side_ents3d);
  other_ents3d = intersect(meshset_3d_ents,other_ents3d); 
  triangles = intersect(meshset_2d_ents,triangles);
  if(verb>3) {
    PetscPrintf(PETSC_COMM_WORLD,"triangles %u\n",triangles.size());
    PetscPrintf(PETSC_COMM_WORLD,"side_ents3d %u\n",side_ents3d.size());
    PetscPrintf(PETSC_COMM_WORLD,"nodes %u\n",nodes.size());
  }
  typedef RefMoFEMEntity_multiIndex::index<MoABEnt_mi_tag>::type ref_ents_by_ent_type;
  typedef RefMoFEMEntity_multiIndex::index<Composite_EntityType_And_ParentEntityType_mi_tag>::type ref_ents_by_composite;
  ref_ents_by_ent_type &ref_ents_by_ent = refinedMoFemEntities.get<MoABEnt_mi_tag>();
  RefMoFEMEntity_multiIndex_view_by_parent_entity ref_parent_ents_view;
  //create view index by parent entity
  {
    ref_ents_by_composite &ref_ents = refinedMoFemEntities.get<Composite_EntityType_And_ParentEntityType_mi_tag>();
    ref_ents_by_composite::iterator miit;
    ref_ents_by_composite::iterator hi_miit;
    //view by parent type (VERTEX)
    if(inheret_from_bit_level.any()) {
      miit = ref_ents.lower_bound(boost::make_tuple(MBVERTEX,MBVERTEX));
      hi_miit = ref_ents.upper_bound(boost::make_tuple(MBVERTEX,MBVERTEX));
      for(;miit!=hi_miit;miit++) {
	if((miit->get_BitRefLevel()&inheret_from_bit_level) == miit->get_BitRefLevel()) {
	  pair<RefMoFEMEntity_multiIndex_view_by_parent_entity::iterator,bool> p_ref_ent_view;
	  p_ref_ent_view = ref_parent_ents_view.insert(&*miit);
	  if(!p_ref_ent_view.second) {
	    SETERRQ(PETSC_COMM_SELF,1,"non uniqe insertion");
	  }
	}
      }
    }
    //view by parent type (TET and PRISM) 
    {
      miit = ref_ents.lower_bound(boost::make_tuple(MBTET,MBTET));
      hi_miit = ref_ents.upper_bound(boost::make_tuple(MBTET,MBTET));
      for(;miit!=hi_miit;miit++) {
	if((miit->get_BitRefLevel()&inheret_from_bit_level) == miit->get_BitRefLevel()) {
	  pair<RefMoFEMEntity_multiIndex_view_by_parent_entity::iterator,bool> p_ref_ent_view;
	  p_ref_ent_view = ref_parent_ents_view.insert(&*miit);
	  if(!p_ref_ent_view.second) {
	    SETERRQ(PETSC_COMM_SELF,1,"non uniqe insertion");
	  }
	}
      }
      miit = ref_ents.lower_bound(boost::make_tuple(MBPRISM,MBPRISM));
      hi_miit = ref_ents.upper_bound(boost::make_tuple(MBPRISM,MBPRISM));
      for(;miit!=hi_miit;miit++) {
	if((miit->get_BitRefLevel()&inheret_from_bit_level) == miit->get_BitRefLevel()) {
	  pair<RefMoFEMEntity_multiIndex_view_by_parent_entity::iterator,bool> p_ref_ent_view;
	  p_ref_ent_view = ref_parent_ents_view.insert(&*miit);
	  if(!p_ref_ent_view.second) {
	    SETERRQ(PETSC_COMM_SELF,1,"non uniqe insertion");
	  }
	}
      }
    }
  }
  //maps nodes on "father" and "mather" side
  map<
    EntityHandle, /*node on "mather" side*/
    EntityHandle /*node on "father" side*/
    > map_nodes;
  //add new nodes on interface and create map
  Range::iterator nit = nodes.begin();
  double coord[3];
  for(;nit!=nodes.end();nit++) {
    //find ref enet
    ref_ents_by_ent_type::iterator miit_ref_ent = ref_ents_by_ent.find(*nit);
    if(miit_ref_ent == ref_ents_by_ent.end()) {
      SETERRQ(PETSC_COMM_SELF,1,"can not find node in MoFEM database");
    }
    EntityHandle child_entity = 0;
    RefMoFEMEntity_multiIndex::iterator child_it;
    if(inheret_from_bit_level.any()) {
      RefMoFEMEntity_multiIndex_view_by_parent_entity::iterator child_iit;
      child_iit = ref_parent_ents_view.find(*nit);
      if(child_iit != ref_parent_ents_view.end()) {
	child_it = refinedMoFemEntities.find((*child_iit)->get_ref_ent());
	BitRefLevel bit_child = child_it->get_BitRefLevel();
	if( (inheret_from_bit_level&bit_child).any() ) {
	  child_entity = child_it->get_ref_ent();
	}
      }
    }
    //
    bool success;
    if(child_entity == 0) {
      rval = moab.get_coords(&*nit,1,coord); CHKERR_PETSC(rval);	
      EntityHandle new_node;
      rval = moab.create_vertex(coord,new_node); CHKERR(rval);
      map_nodes[*nit] = new_node;
      //create new node on "father" side
      //parent is node on "mather" side
      rval = moab.tag_set_data(th_RefParentHandle,&new_node,1,&*nit); CHKERR_PETSC(rval);
      pair<RefMoFEMEntity_multiIndex::iterator,bool> p_ref_ent = refinedMoFemEntities.insert(RefMoFEMEntity(moab,new_node));
      //set ref bit level to node on "father" side
      success = refinedMoFemEntities.modify(p_ref_ent.first,RefMoFEMEntity_change_add_bit(bit));
      if(!success) SETERRQ(PETSC_COMM_SELF,1,"modification unsucceeded");
    } else {
      map_nodes[*nit] = child_entity;
      //set ref bit level to node on "father" side
      success = refinedMoFemEntities.modify(child_it,RefMoFEMEntity_change_add_bit(bit));
      if(!success) SETERRQ(PETSC_COMM_SELF,1,"modification unsucceeded");
    }
    //set ref bit level to node on "mather" side
    success = refinedMoFemEntities.modify(miit_ref_ent,RefMoFEMEntity_change_add_bit(bit));
    if(!success) SETERRQ(PETSC_COMM_SELF,1,"modification unsucceeded");
  }
  //crete meshset for new mesh bit level
  EntityHandle meshset_for_bit_level;
  rval = moab.create_meshset(MESHSET_SET,meshset_for_bit_level); CHKERR_PETSC(rval);
  //subtract those elements which will be refined, i.e. disconetcted form other side elements, and connected to new prisms, if they area created
  meshset_3d_ents = subtract(meshset_3d_ents,side_ents3d);
  rval = moab.add_entities(meshset_for_bit_level,meshset_3d_ents); CHKERR_PETSC(rval);
  for(int dd = 0;dd<3;dd++) {
    Range ents_dd;
    rval = moab.get_adjacencies(meshset_3d_ents,dd,false,ents_dd,moab::Interface::UNION); CHKERR_PETSC(rval);
    rval = moab.add_entities(meshset_for_bit_level,ents_dd); CHKERR_PETSC(rval);
  }
  //
  //typedef RefMoFEMEntity_multiIndex::index<Composite_EntityHandle_And_ParentEntityType_mi_tag>::type ref_ent_by_composite;
  //ref_ent_by_composite &by_composite = refinedMoFemEntities.get<Composite_EntityHandle_And_ParentEntityType_mi_tag>();
  //create new 3d ents on "father" side
  Range new_3d_ents;
  Range::iterator eit3d = side_ents3d.begin();
  for(;eit3d!=side_ents3d.end();eit3d++) {
    ref_ents_by_ent_type::iterator miit_ref_ent = ref_ents_by_ent.find(*eit3d);
    if(miit_ref_ent==ref_ents_by_ent.end()) SETERRQ(PETSC_COMM_SELF,1,"tet not in database");
    int num_nodes; 
    const EntityHandle* conn;
    rval = moab.get_connectivity(*eit3d,conn,num_nodes,true); CHKERR_PETSC(rval);
    EntityHandle new_conn[num_nodes];
    int nb_new_conn = 0;
    int ii = 0;
    for(; ii<num_nodes; ii++) {
      map<EntityHandle,EntityHandle>::iterator mit = map_nodes.find(conn[ii]);
      if(mit != map_nodes.end()) {
	new_conn[ii] = mit->second;
	nb_new_conn++;
	if(verb>6) {
	  PetscPrintf(PETSC_COMM_WORLD,"nodes %u -> %d\n",conn[ii],new_conn[ii]);
	}
      } else {
	new_conn[ii] = conn[ii];
      }
    }
    if(nb_new_conn==0) {
      if(verb>3) {
	EntityHandle meshset_error_out;
	rval = moab.create_meshset(MESHSET_SET,meshset_error_out); CHKERR_PETSC(rval);
	rval = moab.add_entities(meshset_error_out,&*eit3d,1); CHKERR_PETSC(rval);
	ierr = moab.write_file("error_out.vtk","VTK","",&meshset_error_out,1); CHKERRQ(ierr);
      }
      SETERRQ1(PETSC_COMM_SELF,1,"database inconsistency, in side_ent3 is a tet which has no common node with interface, num_nodes = %d",num_nodes);
    }
    //here is created new or prism is on inteface
    EntityHandle existing_ent = 0;
    /* check if tet element whith new connectivity is in database*/ 
    RefMoFEMEntity_multiIndex_view_by_parent_entity::iterator child_iit;
    EntityHandle parent_ent = *eit3d;
    do {
      child_iit = ref_parent_ents_view.find(parent_ent);
      if(child_iit != ref_parent_ents_view.end()) {
	const EntityHandle* conn_ref_tet;
	rval = moab.get_connectivity((*child_iit)->get_ref_ent(),conn_ref_tet,num_nodes,true); CHKERR_PETSC(rval);
	int nn = 0;
	for(;nn<num_nodes;nn++) {
	  if(conn_ref_tet[nn]!=new_conn[nn]) {
	    break;
	  }
	}
	if(nn == num_nodes) {
	  existing_ent = (*child_iit)->get_ref_ent();
	  break;
	} else {
	  parent_ent = (*child_iit)->get_ref_ent();
	}
      }
    } while (child_iit != ref_parent_ents_view.end());
    switch (moab.type_from_handle(*eit3d)) {
      case MBTET: {
	ref_ents_by_ent_type::iterator child_it;
	EntityHandle tet;
	if(existing_ent == 0) {
	  Range new_conn_tet;
	  rval = moab.get_adjacencies(new_conn,4,3,false,new_conn_tet); CHKERR(rval);
	  if(new_conn_tet.empty()) {
	    rval = moab.create_element(MBTET,new_conn,4,tet); CHKERR_PETSC(rval);
	    rval = moab.tag_set_data(th_RefParentHandle,&tet,1,&*eit3d); CHKERR_PETSC(rval);
	  } else {
	    RefMoFEMElement_multiIndex::index<MoABEnt_mi_tag>::type::iterator rit,new_rit;
	    rit = refinedMoFemElements.get<MoABEnt_mi_tag>().find(*eit3d);
	    if(rit==refinedMoFemElements.get<MoABEnt_mi_tag>().end()) {
	      SETERRQ(PETSC_COMM_SELF,1,"can't find this in database");
	    }
	    new_rit  = refinedMoFemElements.get<MoABEnt_mi_tag>().find(*new_conn_tet.begin());
	    if(new_rit==refinedMoFemElements.get<MoABEnt_mi_tag>().end()) {
	      SETERRQ(PETSC_COMM_SELF,1,"can't find this in database");
	    }    
	    tet = *new_conn_tet.begin();
	    /*ostringstream ss;
	    ss << "nb new conns: " << nb_new_conn << endl;
	    ss << "new_conn_tets.size() " << new_conn_tet.size() << endl;
	    ss << "data inconsistency\n";
	    ss << "this ent:\n";
	    ss << *rit->ref_ptr << endl;
	    ss << "found this ent:\n";
	    ss << *new_rit->ref_ptr << endl;
	    SETERRQ(PETSC_COMM_SELF,1,ss.str().c_str());*/
	  }
	} else {
	  tet = existing_ent;
	}
	rval = moab.add_entities(meshset_for_bit_level,&tet,1); CHKERR_PETSC(rval);
	rval = moab.add_entities(meshset_for_bit_level,new_conn,4); CHKERR_PETSC(rval);
	new_3d_ents.insert(tet);
      } break;
      case MBPRISM: {
	EntityHandle prism;
	if(verb>3) {
	  PetscPrintf(PETSC_COMM_WORLD,"prims nb_new_nodes %d\n",nb_new_conn);
	}
	if(existing_ent == 0) {
	  Range new_conn_prism;
	  rval = moab.get_adjacencies(new_conn,6,3,false,new_conn_prism); CHKERR(rval);
	  if(new_conn_prism.empty()) {
	    rval = moab.create_element(MBPRISM,new_conn,6,prism); CHKERR_PETSC(rval);
	    rval = moab.tag_set_data(th_RefParentHandle,&prism,1,&*eit3d); CHKERR_PETSC(rval);
	  } else {
	    SETERRQ(PETSC_COMM_SELF,1,"data inconsistency");
	  }
	} else {
	  prism = existing_ent;
	}
	rval = moab.add_entities(meshset_for_bit_level,&prism,1); CHKERR_PETSC(rval);
	rval = moab.add_entities(meshset_for_bit_level,new_conn,4); CHKERR_PETSC(rval);
	new_3d_ents.insert(prism);
      } break;
      default: 
	SETERRQ(PETSC_COMM_SELF,1,"not implemented");
    }
  }
  Range new_ents; 
  //create new entities by adjecies form new tets
  rval = moab.get_adjacencies(new_3d_ents.subset_by_type(MBTET),2,true,new_ents,Interface::UNION); CHKERR_PETSC(rval);
  rval = moab.get_adjacencies(new_3d_ents.subset_by_type(MBTET),1,true,new_ents,Interface::UNION); CHKERR_PETSC(rval);
  //Tags for setting side
  Tag th_interface_side;
  const int def_side[] = {0};
  rval = moab.tag_get_handle("INTERFACE_SIDE",1,MB_TYPE_INTEGER,
      th_interface_side,MB_TAG_CREAT|MB_TAG_SPARSE,def_side); CHKERR_PETSC(rval);
  //add new edges and triangles to mofem database
  Range ents; 
  rval = moab.get_adjacencies(triangles,1,false,ents,Interface::UNION); CHKERR_PETSC(rval);
  ents.insert(triangles.begin(),triangles.end());
  Range new_ents_in_database; //this range contains all new entities
  Range::iterator eit = ents.begin();
  for(;eit!=ents.end();eit++) {
    int num_nodes; 
    const EntityHandle* conn;
    rval = moab.get_connectivity(*eit,conn,num_nodes,true); CHKERR_PETSC(rval);
    EntityHandle new_conn[num_nodes];
    int nb_new_conn = 0;
    int ii = 0;
    for(;ii<num_nodes; ii++) {
      map<EntityHandle,EntityHandle>::iterator mit = map_nodes.find(conn[ii]);
      if(mit != map_nodes.end()) {
	new_conn[ii] = mit->second;
	nb_new_conn++;
	if(verb>6) {
	  PetscPrintf(PETSC_COMM_WORLD,"nodes %u -> %d\n",conn[ii],new_conn[ii]);
	}
      } else {
	new_conn[ii] = conn[ii];
      }
    }
    if(nb_new_conn==0) continue;
    ref_ents_by_ent_type::iterator miit_ref_ent = ref_ents_by_ent.find(*eit);
    if(miit_ref_ent == ref_ents_by_ent.end()) {
      SETERRQ(PETSC_COMM_SELF,1,"this entity (edge or tri) should be already in database");
    }
    Range new_ent; //contains all entities (edges or triangles) added to mofem database
    switch (moab.type_from_handle(*eit)) {
      case MBTRI: {
	  //get entity based on its connectivity
	  rval = moab.get_adjacencies(new_conn,3,2,false,new_ent); CHKERR_PETSC(rval);
	  if(new_ent.size() != 1) SETERRQ(PETSC_COMM_SELF,1,"this tri should be in moab database"); 
	  int new_side = 1;
	  rval = moab.tag_set_data(th_interface_side,&*new_ent.begin(),1,&new_side); CHKERR_PETSC(rval);
	  if(verb>3) PetscPrintf(PETSC_COMM_WORLD,"new_ent %u\n",new_ent.size());
	  //add prism element
	  if(add_iterfece_entities) {
	    if(inheret_from_bit_level.any()) {
	      SETERRQ(PETSC_COMM_SELF,1,"not implemented for inheret_from_bit_level");
	    }
	    //set prism connectivity
	    EntityHandle prism_conn[6] = { 
	      conn[0],conn[1],conn[2],
	      new_conn[0],new_conn[1],new_conn[2] 
	    };
	    Range new_conn_prism;
	    rval = moab.get_adjacencies(new_conn,6,3,false,new_conn_prism); CHKERR(rval);
	    EntityHandle prism;
	    if(new_conn_prism.empty()) {
	      rval = moab.create_element(MBPRISM,prism_conn,6,prism); CHKERR_PETSC(rval);
	    } else {
	      SETERRQ(PETSC_COMM_SELF,1,"data inconsistency");
	    }
	    ierr = add_prism_to_mofem_database(prism,verb); CHKERRQ(ierr);
	    rval = moab.add_entities(meshset_for_bit_level,&prism,1); CHKERR_PETSC(rval);
	  }
	} break;
      case MBEDGE: {
	  rval = moab.get_adjacencies(new_conn,2,1,false,new_ent); CHKERR_PETSC(rval);
	  if(new_ent.size()!=1) {
	    ParallelComm* pcomm = ParallelComm::get_pcomm(&moab,MYPCOMM_INDEX);
	    if(pcomm->rank()==0) {
	      EntityHandle out_meshset;
	      rval = moab.create_meshset(MESHSET_SET,out_meshset); CHKERR_PETSC(rval);
	      rval = moab.add_entities(out_meshset,&*eit,1); CHKERR_PETSC(rval);
	      rval = moab.write_file("debug_get_msId_3dENTS_split_sides.vtk","VTK","",&out_meshset,1); CHKERR_PETSC(rval);
	      rval = moab.delete_entities(&out_meshset,1); CHKERR_PETSC(rval);
	      rval = moab.create_meshset(MESHSET_SET,out_meshset); CHKERR_PETSC(rval);
	      rval = moab.add_entities(out_meshset,side_ents3d); CHKERR_PETSC(rval);
	      rval = moab.write_file("debug_get_msId_3dENTS_split_sides_side_ents3d.vtk","VTK","",&out_meshset,1); CHKERR_PETSC(rval);
	      rval = moab.delete_entities(&out_meshset,1); CHKERR_PETSC(rval);
	      rval = moab.create_meshset(MESHSET_SET,out_meshset); CHKERR_PETSC(rval);
	      rval = moab.add_entities(out_meshset,other_ents3d); CHKERR_PETSC(rval);
	      rval = moab.write_file("debug_get_msId_3dENTS_split_sides_other_ents3d.vtk","VTK","",&out_meshset,1); CHKERR_PETSC(rval);
	      rval = moab.delete_entities(&out_meshset,1); CHKERR_PETSC(rval);
	      rval = moab.create_meshset(MESHSET_SET,out_meshset); CHKERR_PETSC(rval);
	      rval = moab.add_entities(out_meshset,triangles); CHKERR_PETSC(rval);
	      rval = moab.write_file("debug_get_msId_3dENTS_split_triangles.vtk","VTK","",&out_meshset,1); CHKERR_PETSC(rval);
	      rval = moab.delete_entities(&out_meshset,1); CHKERR_PETSC(rval);
	    }
	    SETERRQ2(PETSC_COMM_SELF,1,
	      "this edge should be in moab database new_ent.size() = %u nb_new_conn = %d",
	      new_ent.size(),nb_new_conn);
	  }
	} break;
      default:
	SETERRQ(PETSC_COMM_SELF,1,"houston we have a problem !!!");
    }
    if(new_ent.size()!=1) {
      SETERRQ1(PETSC_COMM_SELF,1,"new_ent.size() = %u, size always should be 1",new_ent.size());
    }
    //set parent 
    rval = moab.tag_set_data(th_RefParentHandle,&*new_ent.begin(),1,&*eit); CHKERR_PETSC(rval);
    //add to database
    pair<RefMoFEMEntity_multiIndex::iterator,bool> p_ref_ent = refinedMoFemEntities.insert(RefMoFEMEntity(moab,new_ent[0]));
    refinedMoFemEntities.modify(p_ref_ent.first,RefMoFEMEntity_change_add_bit(bit));
    new_ents_in_database.insert(new_ent.begin(),new_ent.end());
  }
  //all other entities, some ents like triangles and faces on the side of tets
  Range side_adj_faces_and_edges;
  rval = moab.get_adjacencies(side_ents3d.subset_by_type(MBTET),1,true,side_adj_faces_and_edges,Interface::UNION); CHKERR_PETSC(rval);
  rval = moab.get_adjacencies(side_ents3d.subset_by_type(MBTET),2,true,side_adj_faces_and_edges,Interface::UNION); CHKERR_PETSC(rval);
  //subtract entities already added to mofem database
  side_adj_faces_and_edges = subtract(side_adj_faces_and_edges,new_ents_in_database);
  eit = side_adj_faces_and_edges.begin();
  for(;eit!=side_adj_faces_and_edges.end();eit++) {
    int num_nodes; 
    const EntityHandle* conn;
    rval = moab.get_connectivity(*eit,conn,num_nodes,true); CHKERR_PETSC(rval);
    EntityHandle new_conn[num_nodes];
    int nb_new_conn = 0;
    int ii = 0;
    for(;ii<num_nodes; ii++) {
      map<EntityHandle,EntityHandle>::iterator mit = map_nodes.find(conn[ii]);
      if(mit != map_nodes.end()) {
	new_conn[ii] = mit->second;
	nb_new_conn++;
	if(verb>6) {
	  PetscPrintf(PETSC_COMM_WORLD,"nodes %u -> %d\n",conn[ii],new_conn[ii]);
	}
      } else {
	new_conn[ii] = conn[ii];
      }
    }
    if(nb_new_conn==0) continue;
    ref_ents_by_ent_type::iterator miit_ref_ent = ref_ents_by_ent.find(*eit);
    if(miit_ref_ent == ref_ents_by_ent.end()) {
      SETERRQ1(PETSC_COMM_SELF,1,"entity should be in MoFem database, num_nodes = %d",num_nodes);
    }
    Range new_ent;
    switch (moab.type_from_handle(*eit)) {
      case MBTRI: {
	  rval = moab.get_adjacencies(new_conn,3,2,false,new_ent); CHKERR_PETSC(rval);
	}
	break;
      case MBEDGE: {
	  rval = moab.get_adjacencies(new_conn,2,1,false,new_ent); CHKERR_PETSC(rval);
	}
	break;
      default:
	SETERRQ(PETSC_COMM_SELF,1,"houston we have a problem");
    }
    if(new_ent.size()!=1) {
      SETERRQ1(PETSC_COMM_SELF,1,"database inconsistency, new_ent.size() = %u",new_ent.size());
    }
    //add entity to mofem database
    rval = moab.tag_set_data(th_RefParentHandle,&*new_ent.begin(),1,&*eit); CHKERR_PETSC(rval);
    pair<RefMoFEMEntity_multiIndex::iterator,bool> p_ref_ent 
      = refinedMoFemEntities.insert(RefMoFEMEntity(moab,new_ent[0]));
    refinedMoFemEntities.modify(p_ref_ent.first,RefMoFEMEntity_change_add_bit(bit));
    if(verb>3) PetscPrintf(PETSC_COMM_WORLD,"new_ent %u\n",new_ent.size());
    new_ents_in_database.insert(new_ent.begin(),new_ent.end());
  }
  //add new prisms which parents are part of other intefaces
  Range new_3d_prims = new_3d_ents.subset_by_type(MBPRISM);
  for(Range::iterator pit = new_3d_prims.begin();pit!=new_3d_prims.end();pit++) {
    ierr = add_prism_to_mofem_database(*pit,verb); CHKERRQ(ierr);
    //get parent entity
    EntityHandle parent_prism;
    rval = moab.tag_get_data(th_RefParentHandle,&*pit,1,&parent_prism); CHKERR_PETSC(rval);
    const EntityHandle root_meshset = moab.get_root_set();
    if(parent_prism == root_meshset)  {
      SETERRQ(PETSC_COMM_SELF,1,"this prism should have parent");
    }
    if(moab.type_from_handle(parent_prism)!=MBPRISM) {
      SETERRQ(PETSC_COMM_SELF,1,"this prism should have parent which is prism as well");
    }
    int num_nodes;
    //parent prism
    const EntityHandle* conn_parent;
    rval = moab.get_connectivity(parent_prism,conn_parent,num_nodes,true); CHKERR_THROW(rval);
    Range face_side3_parent,face_side4_parent;
    rval = moab.get_adjacencies(conn_parent,3,2,false,face_side3_parent); CHKERR_PETSC(rval);
    rval = moab.get_adjacencies(&conn_parent[3],3,2,false,face_side4_parent); CHKERR_PETSC(rval);
    if(face_side3_parent.size()!=1) {
      SETERRQ1(PETSC_COMM_SELF,1,"parent face3.size() = %u",face_side3_parent.size());
    }
    if(face_side4_parent.size()!=1) {
      SETERRQ1(PETSC_COMM_SELF,1,"parent face4.size() = %u",face_side4_parent.size());
    }
    //new prism
    const EntityHandle* conn;
    rval = moab.get_connectivity(*pit,conn,num_nodes,true); CHKERR_THROW(rval);
    Range face_side3,face_side4;
    rval = moab.get_adjacencies(conn,3,2,false,face_side3); CHKERR_PETSC(rval);
    rval = moab.get_adjacencies(&conn[3],3,2,false,face_side4); CHKERR_PETSC(rval);
    if(face_side3.size()!=1) {
      SETERRQ(PETSC_COMM_SELF,1,"face3 is missing");
    }
    if(face_side4.size()!=1) {
      SETERRQ(PETSC_COMM_SELF,1,"face4 is missing");
    }
    //
    vector<EntityHandle> face(2),parent_face(2);
    face[0] = *face_side3.begin();
    face[1] = *face_side4.begin();
    parent_face[0] = *face_side3_parent.begin();
    parent_face[1] = *face_side4_parent.begin();
    for(int ff = 0;ff<2;ff++) {
      if(parent_face[ff] == face[ff]) continue;
      int interface_side;
      rval = moab.tag_get_data(th_interface_side,&parent_face[ff],1,&interface_side); CHKERR_PETSC(rval);
      rval = moab.tag_set_data(th_interface_side,&face[ff],1,&interface_side); CHKERR_PETSC(rval);
      EntityHandle parent_tri;
      rval = moab.tag_get_data(th_RefParentHandle,&face[ff],1,&parent_tri); CHKERR_PETSC(rval);
      if(parent_tri != parent_face[ff]) {
	SETERRQ1(PETSC_COMM_SELF,1,"wrong parent %lu",parent_tri);
      }
      if(new_ents_in_database.find(face[ff])==new_ents_in_database.end()) {
	RefMoFEMEntity_multiIndex::index<MoABEnt_mi_tag>::type::iterator miit_ref_ent 
	    = refinedMoFemEntities.get<MoABEnt_mi_tag>().find(face[ff]);
	if(miit_ref_ent==refinedMoFemEntities.get<MoABEnt_mi_tag>().end()) {
	  SETERRQ(PETSC_COMM_SELF,1,"this is not in database, but should not be");
	}
      }
    }
  }
  //finalise by adding new tets and prism ti bitlelvel
  ierr = seed_ref_level_3D(meshset_for_bit_level,bit); CHKERRQ(ierr);
  rval = moab.delete_entities(&meshset_for_bit_level,1); CHKERR_PETSC(rval);
  ierr = moab.clear_meshset(&children[0],3); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
PetscErrorCode FieldCore::add_prism_to_mofem_database(const EntityHandle prism,int verb) {
  PetscFunctionBegin;
  if(verb==-1) verb = verbose;
  try {
    pair<RefMoFEMEntity_multiIndex::iterator,bool> p_ent = refinedMoFemEntities.insert(RefMoFEMEntity(moab,prism));
    if(p_ent.second) {
      pair<RefMoFEMElement_multiIndex::iterator,bool> p_MoFEMFiniteElement;
      p_MoFEMFiniteElement = refinedMoFemElements.insert(
	ptrWrapperRefMoFEMElement(new RefMoFEMElement_PRISM(moab,&*p_ent.first)));
      int num_nodes;
      const EntityHandle* conn;
      rval = moab.get_connectivity(prism,conn,num_nodes,true); CHKERR_THROW(rval);
      Range face_side3,face_side4;
      rval = moab.get_adjacencies(conn,3,2,false,face_side3); CHKERR_PETSC(rval);
      rval = moab.get_adjacencies(&conn[3],3,2,false,face_side4); CHKERR_PETSC(rval);
      if(face_side3.size()!=1) SETERRQ(PETSC_COMM_SELF,1,"prism don't have side face 3");
      if(face_side4.size()!=1) SETERRQ(PETSC_COMM_SELF,1,"prims don't have side face 4");
      p_MoFEMFiniteElement.first->get_side_number_ptr(moab,*face_side3.begin());
      p_MoFEMFiniteElement.first->get_side_number_ptr(moab,*face_side4.begin());
    } 
  } catch (const char* msg) {
    SETERRQ(PETSC_COMM_SELF,1,msg);
  }
  PetscFunctionReturn(0);
}

}
