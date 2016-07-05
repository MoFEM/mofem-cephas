/** \file MeshRefinmentCore.cpp
 * \brief FIXME this is not so good implementation
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

#include <EntityRefine.hpp>

namespace MoFEM {

MeshRefinment::MeshRefinment(moab::Interface &moab) {
  MoABErrorCode rval;
  const int def_type[] = {0,0};
  rval = moab.tag_get_handle(
    "_RefType",
    2,
    MB_TYPE_INTEGER,
    th_RefType,
    MB_TAG_CREAT|MB_TAG_SPARSE,
    def_type
  ); MOAB_THROW(rval);
}


PetscErrorCode Core::add_verices_in_the_middel_of_edges(const EntityHandle meshset,const BitRefLevel &bit,const bool recursive,int verb) {
  PetscFunctionBegin;
  if(verb==-1) verb = verbose;
  Range edges;
  rval = moab.get_entities_by_type(meshset,MBEDGE,edges,recursive);  CHKERRQ_MOAB(rval);
  if(edges.empty()) {
    Range tets;
    rval = moab.get_entities_by_type(meshset,MBTET,tets,recursive); CHKERRQ_MOAB(rval);
    rval = moab.get_adjacencies(tets,1,true,edges,Interface::UNION); CHKERRQ_MOAB(rval);
    if(tets.empty()) {
      Range prisms;
      rval = moab.get_entities_by_type(meshset,MBPRISM,prisms,recursive); CHKERRQ_MOAB(rval);
      for(Range::iterator pit = prisms.begin();pit!=prisms.end();pit++) {
        const EntityHandle* conn;
        int num_nodes;
        rval = moab.get_connectivity(*pit,conn,num_nodes,true);  CHKERRQ_MOAB(rval);
        assert(num_nodes==6);
        //
        Range edge;
        rval = moab.get_adjacencies(&conn[0],2,1,true,edge); CHKERRQ_MOAB(rval);
        assert(edge.size()==1);
        edges.insert(edge[0]);
        edge.clear();
        rval = moab.get_adjacencies(&conn[1],2,1,true,edge); CHKERRQ_MOAB(rval);
        assert(edge.size()==1);
        edges.insert(edge[0]);
        EntityHandle conn_edge2[] = { conn[2], conn[0] };
        edge.clear();
        rval = moab.get_adjacencies(conn_edge2,2,1,true,edge); CHKERRQ_MOAB(rval);
        assert(edge.size()==1);
        edges.insert(edge[0]);
        //
        edge.clear();
        rval = moab.get_adjacencies(&conn[3],2,1,true,edge); CHKERRQ_MOAB(rval);
        assert(edge.size()==1);
        edges.insert(edge[0]);
        edge.clear();
        rval = moab.get_adjacencies(&conn[4],2,1,true,edge); CHKERRQ_MOAB(rval);
        assert(edge.size()==1);
        edges.insert(edge[0]);
        EntityHandle conn_edge8[] = { conn[5], conn[3] };
        edge.clear();
        rval = moab.get_adjacencies(conn_edge8,2,1,true,edge); CHKERRQ_MOAB(rval);
        assert(edge.size()==1);
        edges.insert(edge[0]);
      }
    }
  }
  ierr = add_verices_in_the_middel_of_edges(edges,bit,verb); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
PetscErrorCode Core::add_verices_in_the_middel_of_edges(const Range &_edges,const BitRefLevel &bit,int verb) {
  PetscFunctionBegin;
  Range edges = _edges.subset_by_type(MBEDGE);
  typedef RefEntity_multiIndex::index<Composite_EntType_and_ParentEntType_mi_tag>::type ref_ents_by_composite;
  ref_ents_by_composite &ref_ents = refinedEntities.get<Composite_EntType_and_ParentEntType_mi_tag>();
  ref_ents_by_composite::iterator miit = ref_ents.lower_bound(boost::make_tuple(MBVERTEX,MBEDGE));
  ref_ents_by_composite::iterator hi_miit = ref_ents.upper_bound(boost::make_tuple(MBVERTEX,MBEDGE));
  RefEntity_multiIndex_view_by_parent_entity ref_parent_ents_view;
  for(;miit!=hi_miit;miit++) {
    std::pair<RefEntity_multiIndex_view_by_parent_entity::iterator,bool> p_ref_ent_view;
    p_ref_ent_view = ref_parent_ents_view.insert(*miit);
    if(!p_ref_ent_view.second) {
      SETERRQ(PETSC_COMM_SELF,1,"non uniqe insertion");
    }
  }
  if(verb > 0) {
    std::ostringstream ss;
    ss << "ref level " << bit << " nb. edges to refine " << edges.size() << std::endl;
    PetscPrintf(comm,ss.str().c_str());
  }
  Range::iterator eit = edges.begin();
  for(;eit!=edges.end();eit++) {
    //bool add_vertex = false;
    RefEntity_multiIndex_view_by_parent_entity::iterator miit_view;
    if(ref_parent_ents_view.empty()) {
      //add_vertex = true;
    } else {
      miit_view = ref_parent_ents_view.find(*eit);
    }
    if(ref_parent_ents_view.empty()||miit_view == ref_parent_ents_view.end()) {
      const EntityHandle* conn;
      int num_nodes;
      rval = moab.get_connectivity(*eit,conn,num_nodes,true);  CHKERRQ_MOAB(rval);
      if(num_nodes !=2 ) {
        SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"edge should have 2 edges");
      }
      double coords[num_nodes*3];
      rval = moab.get_coords(conn,num_nodes,coords); CHKERRQ_MOAB(rval);
      cblas_daxpy(3,1.,&coords[3],1,coords,1);
      cblas_dscal(3,0.5,coords,1);
      EntityHandle node;
      rval = moab.create_vertex(coords,node); CHKERRQ_MOAB(rval);
      rval = moab.tag_set_data(th_RefParentHandle,&node,1,&*eit); CHKERRQ_MOAB(rval);
      rval = moab.tag_set_data(th_RefBitLevel,&node,1,&bit); CHKERRQ_MOAB(rval);
      std::pair<RefEntity_multiIndex::iterator,bool> p_ent = refinedEntities.insert(
        boost::shared_ptr<RefEntity>(new RefEntity(basicEntityDataPtr,node))
      );
      if(!p_ent.second) SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"this entity is there");
      if(verbose>2) {
        std::ostringstream ss;
        ss << *(p_ent.first) << std::endl;
        PetscPrintf(comm,ss.str().c_str());
      }
    } else {
      const EntityHandle node = (*miit_view)->getRefEnt();
      if((*miit_view)->getEntType() != MBVERTEX) {
        SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"child of edge should be vertex");
      }
      bool success = refinedEntities.modify(refinedEntities.get<Ent_mi_tag>().find(node),RefEntity_change_add_bit(bit));
      if(!success) SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"inconsistency in data");
    }
  }
  PetscFunctionReturn(0);
}
PetscErrorCode Core::refine_TET(const EntityHandle meshset,const BitRefLevel &bit,const bool respect_interface) {
  PetscFunctionBegin;
  Range tets;
  rval = moab.get_entities_by_type(meshset,MBTET,tets,false); CHKERRQ_MOAB(rval);
  ierr = refine_TET(tets,bit,respect_interface); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
PetscErrorCode Core::refine_TET(const Range &_tets,const BitRefLevel &bit,const bool respect_interface) {
  PetscFunctionBegin;
  //FIXME: refinement is based on entity handlers, should work on global ids of nodes, this will allow parallelize agortihm in the future
  PetscFunctionBegin;
  typedef RefEntity_multiIndex::index<Ent_mi_tag>::type ref_ents_by_ent;
  ref_ents_by_ent &ref_ents_ent = refinedEntities.get<Ent_mi_tag>();
  // find all vertices which parent is edge
  typedef RefEntity_multiIndex::index<Composite_EntType_and_ParentEntType_mi_tag>::type ref_ents_by_composite;
  ref_ents_by_composite &ref_ents = refinedEntities.get<Composite_EntType_and_ParentEntType_mi_tag>();
  ref_ents_by_composite::iterator miit = ref_ents.lower_bound(boost::make_tuple(MBVERTEX,MBEDGE));
  ref_ents_by_composite::iterator hi_miit = ref_ents.upper_bound(boost::make_tuple(MBVERTEX,MBEDGE));
  RefEntity_multiIndex_view_by_parent_entity ref_parent_ents_view;
  for(;miit!=hi_miit;miit++) {
    std::pair<RefEntity_multiIndex_view_by_parent_entity::iterator,bool> p_ref_ent_view;
    p_ref_ent_view = ref_parent_ents_view.insert(*miit);
    if(!p_ref_ent_view.second) {
      SETERRQ(PETSC_COMM_SELF,1,"non uniqe insertion");
    }
  }
  typedef RefElement_multiIndex::index<Ent_mi_tag>::type ref_MoFEMFiniteElement_by_ent;
  ref_MoFEMFiniteElement_by_ent &ref_MoFEMFiniteElement = refinedFiniteElements.get<Ent_mi_tag>();
  typedef RefElement_multiIndex::index<Composite_ParentEnt_And_BitsOfRefinedEdges_mi_tag>::type ref_ent_by_composite;
  ref_ent_by_composite &by_composite = refinedFiniteElements.get<Composite_ParentEnt_And_BitsOfRefinedEdges_mi_tag>();
  //
  if(respect_interface) {
    SETERRQ(PETSC_COMM_SELF,1,"not implemented, set last parameter in refine_TET to false");
  }
  //
  Range tets = _tets.subset_by_type(MBTET);
  Range::iterator tit = tets.begin();
  for(;tit!=tets.end();tit++) {
    ref_MoFEMFiniteElement_by_ent::iterator miit2 = ref_MoFEMFiniteElement.find(*tit);
    if(miit2==ref_MoFEMFiniteElement.end()) SETERRQ(PETSC_COMM_SELF,1,"this tet is not in refinedFiniteElements");
    //connectivity
    const EntityHandle* conn;
    int num_nodes;
    moab.get_connectivity(*tit,conn,num_nodes,true);
    //Range _tit_edges;
    //rval = moab.get_adjacencies(&*tit,1,1,true,_tit_edges); CHKERRQ_MOAB(rval);
    for(int nn = 0;nn<num_nodes;nn++) {
      bool success = refinedEntities.modify(refinedEntities.get<Ent_mi_tag>().find(conn[nn]),RefEntity_change_add_bit(bit));
      if(!success) SETERRQ(PETSC_COMM_SELF,1,"can not set refinement bit level to tet node");
    }
    //get edges
    BitRefEdges parent_edges_bit(0);
    EntityHandle edge_new_nodes[6];
    std::fill(&edge_new_nodes[0],&edge_new_nodes[6],no_handle);
    int split_edges[6];
    std::fill(&split_edges[0],&split_edges[6],-1);
    //hash map of nodes (RefEntity) by edges (EntityHandle)
    std::map<EntityHandle /*edge*/,const RefEntity* /*node*/> map_ref_nodes_by_edges;
    for(int ee = 0;ee<6;ee++) {
      EntityHandle edge;
      rval = moab.side_element(*tit,1,ee,edge);  CHKERRQ_MOAB(rval);
      RefEntity_multiIndex_view_by_parent_entity::iterator miit_view;
      miit_view = ref_parent_ents_view.find(edge);
      if(miit_view != ref_parent_ents_view.end()) {
        if(((*miit_view)->getBitRefLevel()&bit).any()) {
          edge_new_nodes[ee] = (*miit_view)->getRefEnt();
          map_ref_nodes_by_edges[(*miit_view)->getParentEnt()] = &**miit_view;
          {
            const EntityHandle* conn_edge;
            int num_nodes;
            moab.get_connectivity(edge,conn_edge,num_nodes,true);
            if(conn_edge[0] == edge_new_nodes[ee]) {
              SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"data inconsistency");
            }
            if(conn_edge[1] == edge_new_nodes[ee]) {
              SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"data inconsistency");
            }
          }
          split_edges[parent_edges_bit.count()] = ee;
          parent_edges_bit.set(ee,1);
        }
      }
    }
    for(int ee = 0;ee<6;ee++) {
      if(edge_new_nodes[ee] == no_handle) continue;
      for(int nn = 0;nn<4;nn++) {
        if(conn[nn] == edge_new_nodes[ee]) {
          std::cerr << conn[0] << " "
          << conn[1] << " "
          << conn[2] << " "
          << conn[3] << " : "
          << edge_new_nodes[ee] << std::endl;
          for(int eee = 0;eee<6;eee++) {
            EntityHandle edge;
            rval = moab.side_element(*tit,1,eee,edge);  CHKERRQ_MOAB(rval);
            const EntityHandle* conn_edge;
            int num_nodes;
            moab.get_connectivity(edge,conn_edge,num_nodes,true);
            std::cerr << conn_edge[0] << " " << conn_edge[1] << std::endl;
          }
          SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"data inconsistency");
        }
      }
    }

    // swap nodes forward
    EntityHandle _conn_[4];
    std::copy(&conn[0],&conn[4],&_conn_[0]);
    // build connectivity for rf tets
    EntityHandle new_tets_conns[8*4];
    std::fill(&new_tets_conns[0],&new_tets_conns[8*4],no_handle);
    int sub_type = -1,nb_new_tets = 0;
    switch (parent_edges_bit.count()) {
      case 0: {
        ref_ents_by_ent::iterator tit_miit;
        tit_miit = ref_ents_ent.find(*tit);
        if(tit_miit==ref_ents_ent.end()) SETERRQ(PETSC_COMM_SELF,1,"data inconsistency");
        bool success = refinedEntities.modify(tit_miit,RefEntity_change_add_bit(bit));
        if(!success) SETERRQ(PETSC_COMM_SELF,1,"impossible tet");
        Range tit_conn;
        rval = moab.get_connectivity(&*tit,1,tit_conn,true); CHKERRQ_MOAB(rval);
        for(Range::iterator nit = tit_conn.begin();nit!=tit_conn.end();nit++) {
          ref_ents_by_ent::iterator nit_miit = ref_ents_ent.find(*nit);
          if(nit_miit==ref_ents_ent.end()) SETERRQ(PETSC_COMM_SELF,1,"can not find face in refinedEntities");
          bool success = refinedEntities.modify(nit_miit,RefEntity_change_add_bit(bit));
          if(!success) SETERRQ(PETSC_COMM_SELF,1,"impossible node");
        }
        Range tit_edges;
        rval = moab.get_adjacencies(&*tit,1,1,false,tit_edges); CHKERRQ_MOAB(rval);
        for(Range::iterator eit = tit_edges.begin();eit!=tit_edges.end();eit++) {
          ref_ents_by_ent::iterator eit_miit = ref_ents_ent.find(*eit);
          if(eit_miit==ref_ents_ent.end()) SETERRQ(PETSC_COMM_SELF,1,"can not find face in refinedEntities");
          bool success = refinedEntities.modify(eit_miit,RefEntity_change_add_bit(bit));
          if(!success) SETERRQ(PETSC_COMM_SELF,1,"impossible edge");
        }
        Range tit_faces;
        rval = moab.get_adjacencies(&*tit,1,2,false,tit_faces); CHKERRQ_MOAB(rval);
        if(tit_faces.size()!=4) SETERRQ(PETSC_COMM_SELF,1,"existing tet in mofem database should have 4 adjacent edges");
        for(Range::iterator fit = tit_faces.begin();fit!=tit_faces.end();fit++) {
          ref_ents_by_ent::iterator fit_miit = ref_ents_ent.find(*fit);
          if(fit_miit==ref_ents_ent.end()) SETERRQ(PETSC_COMM_SELF,1,"can not find face in refinedEntities");
          bool success = refinedEntities.modify(fit_miit,RefEntity_change_add_bit(bit));
          if(!success) SETERRQ(PETSC_COMM_SELF,1,"impossible face");
        }
        continue;
      }
      break;
      case 1:
      sub_type = 0;
      // for(int nn = 0;nn<4;nn++) {
      //   if(_conn_[nn] == edge_new_nodes[split_edges[0]]) {
      //     std::cerr << _conn_[0] << " "
      //     << _conn_[1] << " "
      //     << _conn_[2] << " "
      //     << _conn_[3] << " : "
      //     << edge_new_nodes[split_edges[0]] << std::endl;
      //     SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"data inconsistency");
      //   }
      // }
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
    EntityHandle ref_tets[8];
    std::bitset<8> ref_tets_bit(0);
    ref_ent_by_composite::iterator miit_composite = by_composite.lower_bound(boost::make_tuple(*tit,parent_edges_bit.to_ulong()));
    ref_ent_by_composite::iterator hi_miit_composite = by_composite.upper_bound(boost::make_tuple(*tit,parent_edges_bit.to_ulong()));
    // if(miit_composite!=hi_miit_composite) {
    //   if(distance(miit_composite,hi_miit_composite)!=(unsigned int)nb_new_tets) {
    //     SETERRQ2(PETSC_COMM_SELF,1,"data inconsistency %u != %u",
    //     distance(miit_composite,hi_miit_composite),(unsigned int)nb_new_tets);
    //   }
    // }
    if(distance(miit_composite,hi_miit_composite)==(unsigned int)nb_new_tets) {
      for(int tt = 0;miit_composite!=hi_miit_composite;miit_composite++,tt++) {
        EntityHandle tet = miit_composite->getRefEnt();
        //set ref tets entities
        ref_tets[tt] = tet;
        ref_tets_bit.set(tt,1);
        //add this tet if exist to this ref level
        RefEntity_multiIndex::iterator ref_tet_it;
        ref_tet_it = refinedEntities.find(tet);
        if(ref_tet_it == refinedEntities.end()) {
          SETERRQ(PETSC_COMM_SELF,1,"data inconsistency");
        }
        bool success = refinedEntities.modify(
          ref_tet_it,RefEntity_change_add_bit(bit)
        );
        if(!success) {
          SETERRQ(PETSC_COMM_SELF,1,"modification unsuccessful");
        }
        //verbose
        if(verbose>2) {
          std::ostringstream ss;
          ss << miit_composite->getRefElement() << std::endl;
          PetscPrintf(comm,ss.str().c_str());
        }
      }
    } else {
      //if this element was not refined or was refined with different patterns of split edges create new elements
      for(int tt = 0;tt<nb_new_tets;tt++) {
        if(!ref_tets_bit.test(tt)) {
          if(miit_composite!=hi_miit_composite) {
            Range new_tets_conns_tet;
            rval = moab.get_adjacencies(&new_tets_conns[4*tt],4,2,false,new_tets_conns_tet); CHKERRQ_MOAB(rval);
            if(new_tets_conns_tet.empty()) {
              rval = moab.create_element(MBTET,&new_tets_conns[4*tt],4,ref_tets[tt]); CHKERRQ_MOAB(rval);
            } else {
              if(new_tets_conns_tet.size()!=1) {
                SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"data inconsistency");
              }
              ref_tets[tt] = new_tets_conns_tet[0];
            }
          } else {
            rval = moab.create_element(MBTET,&new_tets_conns[4*tt],4,ref_tets[tt]); CHKERRQ_MOAB(rval);
          }
          Range ref_tet_nodes;
          rval = moab.get_connectivity(&ref_tets[tt],1,ref_tet_nodes,true); CHKERRQ_MOAB(rval);
          if(ref_tet_nodes.size()!=4) {
            std::cerr << tt << " " << nb_new_tets << " " << sub_type << " " << parent_edges_bit << std::endl;
            std::cerr << _conn_[0] << " "
            << _conn_[1] << " "
            << _conn_[2] << " "
            << _conn_[3] << std::endl;

            std::cerr << edge_new_nodes[0] << " "
            << edge_new_nodes[1] << " "
            << edge_new_nodes[2] << " "
            << edge_new_nodes[3] << " "
            << edge_new_nodes[4] << " "
            << edge_new_nodes[5] << std::endl;

            std::cerr << split_edges[0] << " "
            << split_edges[1] << " "
            << split_edges[2] << " "
            << split_edges[3] << " "
            << split_edges[4] << " "
            << split_edges[5] << std::endl;

            std::cerr << new_tets_conns[4*tt+0] << " "
            << new_tets_conns[4*tt+1] << " "
            << new_tets_conns[4*tt+2] << " "
            << new_tets_conns[4*tt+3] << std::endl << std::endl;
            SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"tetrahedral should have 4 nodes");
          }
          int ref_type[2];
          ref_type[0] = parent_edges_bit.count();
          ref_type[1] = sub_type;
          rval = moab.tag_set_data(th_RefType,&ref_tets[tt],1,ref_type); CHKERRQ_MOAB(rval);
          rval = moab.tag_set_data(th_RefParentHandle,&ref_tets[tt],1,&*tit); CHKERRQ_MOAB(rval);
          rval = moab.tag_set_data(th_RefBitLevel,&ref_tets[tt],1,&bit); CHKERRQ_MOAB(rval);
          rval = moab.tag_set_data(th_RefBitEdge,&ref_tets[tt],1,&parent_edges_bit); CHKERRQ_MOAB(rval);
          //add refined entity
          std::pair<RefEntity_multiIndex::iterator,bool> p_MoFEMEntity = refinedEntities.insert(
            boost::shared_ptr<RefEntity>(new RefEntity(basicEntityDataPtr,ref_tets[tt]))
          );
          //add refined element
          std::pair<RefElement_multiIndex::iterator,bool> p_MoFEMFiniteElement;
          try {
            p_MoFEMFiniteElement = refinedFiniteElements.insert(ptrWrapperRefElement(
              boost::shared_ptr<RefElement>(new RefElement_TET(moab,*p_MoFEMEntity.first)))
            );
          } catch (MoFEMException const &e) {
            SETERRQ(PETSC_COMM_SELF,e.errorCode,e.errorMessage);
          }
          //set bit that this element is now in databse
          ref_tets_bit.set(tt);
          if(verbose>2) {
            std::ostringstream ss;
            ss << "add tet: " << *(p_MoFEMFiniteElement.first->getRefElement()) << std::endl;
            PetscPrintf(comm,ss.str().c_str());
          }
        }
      }
    }
    // //debug
    // miit_composite = by_composite.lower_bound(boost::make_tuple(*tit,parent_edges_bit.to_ulong()));
    // hi_miit_composite = by_composite.upper_bound(boost::make_tuple(*tit,parent_edges_bit.to_ulong()));
    // if(miit_composite!=hi_miit_composite) {
    //   if(distance(miit_composite,hi_miit_composite)!=(unsigned int)nb_new_tets) {
    //     SETERRQ2(PETSC_COMM_SELF,1,"data inconsistency %u != %u",
    //     distance(miit_composite,hi_miit_composite),(unsigned int)nb_new_tets);
    //   }
    // }
    //find parents for new edges and faces
    //get tet edges and faces
    Range tit_edges,tit_faces;
    rval = moab.get_adjacencies(&*tit,1,1,false,tit_edges); CHKERRQ_MOAB(rval);
    rval = moab.get_adjacencies(&*tit,1,2,false,tit_faces); CHKERRQ_MOAB(rval);
    Range edges_nodes[6],faces_nodes[4];
    //for edges - add ref nodes
    //edges_nodes[ee] - contains all nodes on edge ee including mid nodes if exist
    Range::iterator eit = tit_edges.begin();
    for(int ee = 0;eit!=tit_edges.end();eit++,ee++) {
      rval = moab.get_connectivity(&*eit,1,edges_nodes[ee],true); CHKERRQ_MOAB(rval);
      std::map<EntityHandle,const RefEntity*>::iterator map_miit = map_ref_nodes_by_edges.find(*eit);
      if(map_miit!=map_ref_nodes_by_edges.end()) {
        edges_nodes[ee].insert(map_miit->second->getRefEnt());
      }
    }
    //for faces - add ref nodes
    //faces_nodes[ff] - contains all nodes on face ff including mid nodes if exist
    Range::iterator fit=tit_faces.begin();
    for(int ff = 0;fit!=tit_faces.end();fit++,ff++) {
      rval = moab.get_connectivity(&*fit,1,faces_nodes[ff],true); CHKERRQ_MOAB(rval);
      Range fit_edges;
      rval = moab.get_adjacencies(&*fit,1,1,false,fit_edges); CHKERRQ_MOAB(rval);
      for(Range::iterator eit2 =  fit_edges.begin();eit2 != fit_edges.end();eit2++) {
        std::map<EntityHandle,const RefEntity*>::iterator map_miit = map_ref_nodes_by_edges.find(*eit2);
        if(map_miit!=map_ref_nodes_by_edges.end()) {
          faces_nodes[ff].insert(map_miit->second->getRefEnt());
        }
      }
    }
    //add ref nodes to tet
    //tet_nodes contains all nodes on tet inluding mid edge nodes
    Range tet_nodes;
    rval = moab.get_connectivity(&*tit,1,tet_nodes,true); CHKERRQ_MOAB(rval);
    for(std::map<EntityHandle,const RefEntity*>::iterator map_miit = map_ref_nodes_by_edges.begin();
    map_miit != map_ref_nodes_by_edges.end();map_miit++) {
      tet_nodes.insert(map_miit->second->getRefEnt());
    }
    Range ref_edges;
    //get all all edges of refined tets
    rval = moab.get_adjacencies(ref_tets,nb_new_tets,1,true,ref_edges,Interface::UNION); CHKERRQ_MOAB(rval);
    //check for all ref edge and set parents
    for(Range::iterator reit = ref_edges.begin();reit!=ref_edges.end();reit++) {
      Range ref_edges_nodes;
      rval = moab.get_connectivity(&*reit,1,ref_edges_nodes,true); CHKERRQ_MOAB(rval);
      if(ref_edges_nodes.size()!=2) {
        SETERRQ(PETSC_COMM_SELF,1,"data inconsistency, edge should have 2 nodes");
      }
      //check if ref edge is an coarse edge
      int ee = 0;
      for(;ee<6;ee++) {
        //two nodes are common (node[0],node[1],ref_node (if exist))
        //this tests if given edge is contained by edge of refined tetrahedral
        if(intersect(edges_nodes[ee],ref_edges_nodes).size()==2) {
          EntityHandle edge = tit_edges[ee];
          rval = moab.tag_set_data(th_RefParentHandle,&*reit,1,&edge); CHKERRQ_MOAB(rval);
          std::pair<RefEntity_multiIndex::iterator,bool> p_ent = refinedEntities.insert(
            boost::shared_ptr<RefEntity>(new RefEntity(basicEntityDataPtr,*reit))
          );
          bool success = refinedEntities.modify(p_ent.first,RefEntity_change_add_bit(bit));
          if(!success) SETERRQ(PETSC_COMM_SELF,1,"impossible to set edge pranet");
          if(p_ent.second) {
            if(verbose>2) {
              std::ostringstream ss;
              ss << "edge parent: " << *(p_ent.first) << std::endl;
              PetscPrintf(comm,ss.str().c_str());
            }
          }
          break;
        }
      }
      if(ee<6) continue; //this refined edge is contained by edge of tetrahedral
      //check if ref edge is in coarse face
      int ff = 0;
      for(;ff<4;ff++) {
        //two nodes are common (node[0],node[1],ref_node (if exist))
        //thi tests if givem edge is contained by face of  tetrahedral
        if(intersect(faces_nodes[ff],ref_edges_nodes).size()==2) {
          EntityHandle face = tit_faces[ff];
          rval = moab.tag_set_data(th_RefParentHandle,&*reit,1,&face); CHKERRQ_MOAB(rval);
          //add edge to refinedEntities
          std::pair<RefEntity_multiIndex::iterator,bool> p_ent = refinedEntities.insert(
            boost::shared_ptr<RefEntity>(new RefEntity(basicEntityDataPtr,*reit))
          );
          bool success = refinedEntities.modify(p_ent.first,RefEntity_change_add_bit(bit));
          if(!success) SETERRQ(PETSC_COMM_SELF,1,"impossible to set edge parent");
          if(p_ent.second) {
            if(verbose>2) {
              std::ostringstream ss;
              ss << "face parent: " << *(p_ent.first) << std::endl;
              PetscPrintf(comm,ss.str().c_str());
            }
          }
          break;
        }
      }
      if(ff<4) continue; //this refined edge is contained by face of tetrahedral
      // check if ref edge is in coarse tetrahedral (i.e. that is internal edge of refined tetrahedral)
      if(intersect(tet_nodes,ref_edges_nodes).size()==2) {
        rval = moab.tag_set_data(th_RefParentHandle,&*reit,1,&*tit); CHKERRQ_MOAB(rval);
        //add edge to refinedEntities
        std::pair<RefEntity_multiIndex::iterator,bool> p_ent = refinedEntities.insert(
          boost::shared_ptr<RefEntity>(new RefEntity(basicEntityDataPtr,*reit))
        );
        bool success = refinedEntities.modify(p_ent.first,RefEntity_change_add_bit(bit));
        if(!success) SETERRQ(PETSC_COMM_SELF,1,"impossible to set edge parent");
        if(p_ent.second) {
          if(verbose>2) {
            std::ostringstream ss;
            ss << "tet parent: " << *(p_ent.first) << std::endl;
            PetscPrintf(comm,ss.str().c_str());
          }
        }
        continue;
      }

      //this will help to debug error
      Range ref_edges_nodes_tets;
      rval = moab.get_adjacencies(ref_edges_nodes,3,true,ref_edges_nodes_tets,Interface::UNION); CHKERRQ_MOAB(rval);
      ref_edges_nodes_tets = intersect(ref_edges_nodes_tets,tets);
      ref_edges_nodes_tets.insert(*reit);

      EntityHandle meshset_out;
      rval = moab.create_meshset(MESHSET_SET,meshset_out); CHKERRQ_MOAB(rval);
      rval = moab.add_entities(meshset_out,ref_edges_nodes_tets); CHKERRQ_MOAB(rval);
      rval = moab.write_file("debug_error.vtk","VTK","",&meshset_out,1); CHKERRQ_MOAB(rval);
      rval = moab.delete_entities(&meshset_out,1); CHKERRQ_MOAB(rval);

      //refined edge is not child of any edge, face or tetrahedral, this is imposible edge
      SETERRQ(PETSC_COMM_SELF,1,"impossible refined edge");
    }
    Range ref_faces;
    rval = moab.get_adjacencies(ref_tets,nb_new_tets,2,true,ref_faces,Interface::UNION); CHKERRQ_MOAB(rval);
    Tag th_interface_side;
    const int def_side[] = {0};
    rval = moab.tag_get_handle("INTERFACE_SIDE",1,MB_TYPE_INTEGER,
    th_interface_side,MB_TAG_CREAT|MB_TAG_SPARSE,def_side); CHKERRQ_MOAB(rval);
    // check for all ref faces
    for(Range::iterator rfit = ref_faces.begin();rfit!=ref_faces.end();rfit++) {
      Range ref_faces_nodes;
      rval = moab.get_connectivity(&*rfit,1,ref_faces_nodes,true); CHKERRQ_MOAB(rval);
      // check if ref face is in coarse face
      int ff = 0;
      for(;ff<4;ff++) {
        //check if refined edge is contained by face of tetrahedral
        if(intersect(faces_nodes[ff],ref_faces_nodes).size()==3) {
          EntityHandle face = tit_faces[ff];
          rval = moab.tag_set_data(th_RefParentHandle,&*rfit,1,&face); CHKERRQ_MOAB(rval);
          int side = 0;
          //set face side if it is on interface
          rval = moab.tag_get_data(th_interface_side,&face,1,&side); CHKERRQ_MOAB(rval);
          rval = moab.tag_set_data(th_interface_side,&*rfit,1,&side); CHKERRQ_MOAB(rval);
          //add face to refinedEntities
          std::pair<RefEntity_multiIndex::iterator,bool> p_ent = refinedEntities.insert(
            boost::shared_ptr<RefEntity>(new RefEntity(basicEntityDataPtr,*rfit))
          );
          bool success = refinedEntities.modify(p_ent.first,RefEntity_change_add_bit(bit));
          if(!success) SETERRQ(PETSC_COMM_SELF,1,"impossible to set face parent");
          if(p_ent.second) {
            if(verbose>2) {
              std::ostringstream ss;
              ss << "face parent: " << *(p_ent.first) << std::endl;
              PetscPrintf(comm,ss.str().c_str());
            }
          }
          break;
        }
      }
      if(ff<4) continue; //this face is contained by one of tetrahedrons
      //check if ref face is in coarse tetrahedral
      //this is ref face which is contained by tetrahedral volume
      if(intersect(tet_nodes,ref_faces_nodes).size()==3) {
        rval = moab.tag_set_data(th_RefParentHandle,&*rfit,1,&*tit); CHKERRQ_MOAB(rval);
        std::pair<RefEntity_multiIndex::iterator,bool> p_ent = refinedEntities.insert(
          boost::shared_ptr<RefEntity>(new RefEntity(basicEntityDataPtr,*rfit))
        );
        //add face to refinedEntities
        bool success = refinedEntities.modify(p_ent.first,RefEntity_change_add_bit(bit));
        if(!success) SETERRQ(PETSC_COMM_SELF,1,"impossible to set face parent");
        if(p_ent.second) {
          if(verbose>2) {
            std::ostringstream ss;
            ss << "tet parent: " << *(p_ent.first) << std::endl;
            PetscPrintf(comm,ss.str().c_str());
          }
        }
        continue;
      }
      SETERRQ(PETSC_COMM_SELF,1,"impossible refined face");
    }
  }
  PetscFunctionReturn(0);
}
PetscErrorCode Core::refine_PRISM(const EntityHandle meshset,const BitRefLevel &bit,int verb) {
  //FIXME: refinement is based on entity handlers, should work on global ids of nodes, this will allow parallelize algorithm in the future
  PetscFunctionBegin;
  if(verb==-1) verb = verbose;
  typedef RefEntity_multiIndex::index<Ent_mi_tag>::type ref_ENTs_by_ent;
  typedef RefElement_multiIndex::index<Composite_ParentEnt_And_BitsOfRefinedEdges_mi_tag>::type ref_fe_by_composite;
  ref_fe_by_composite &ref_fe_by_comp = refinedFiniteElements.get<Composite_ParentEnt_And_BitsOfRefinedEdges_mi_tag>();
  //find all vertices which parent is edge
  typedef RefEntity_multiIndex::index<Composite_EntType_and_ParentEntType_mi_tag>::type ref_ents_by_composite;
  ref_ents_by_composite &ref_ents_by_comp = refinedEntities.get<Composite_EntType_and_ParentEntType_mi_tag>();
  ref_ents_by_composite::iterator miit = ref_ents_by_comp.lower_bound(boost::make_tuple(MBVERTEX,MBEDGE));
  ref_ents_by_composite::iterator hi_miit = ref_ents_by_comp.upper_bound(boost::make_tuple(MBVERTEX,MBEDGE));
  RefEntity_multiIndex_view_by_parent_entity ref_parent_ents_view;
  for(;miit!=hi_miit;miit++) {
    std::pair<RefEntity_multiIndex_view_by_parent_entity::iterator,bool> p_ref_ent_view;
    p_ref_ent_view = ref_parent_ents_view.insert(*miit);
    if(!p_ref_ent_view.second) {
      SETERRQ(PETSC_COMM_SELF,1,"non uniqe insertion");
    }
  }
  Range prisms;
  rval = moab.get_entities_by_type(meshset,MBPRISM,prisms,false); CHKERRQ_MOAB(rval);
  Range::iterator pit = prisms.begin();
  for(;pit!=prisms.end();pit++) {
    ref_ENTs_by_ent::iterator miit_prism = refinedEntities.get<Ent_mi_tag>().find(*pit);
    if(miit_prism==refinedEntities.end()) SETERRQ(PETSC_COMM_SELF,1,"this prism is not in ref database");
    if(verb>3) {
      std::ostringstream ss;
      ss << "ref prism " << *miit << std::endl;
      PetscPrintf(comm,ss.str().c_str());
    }
    //prism connectivity
    int num_nodes;
    const EntityHandle* conn;
    rval = moab.get_connectivity(*pit,conn,num_nodes,true); CHKERRQ_MOAB(rval);
    assert(num_nodes==6);
    // edges connectivity
    EntityHandle edges[6];
    for(int ee = 0;ee<3; ee++) {
      rval = moab.side_element(*pit,1,ee,edges[ee]); CHKERRQ_MOAB(rval);
    }
    for(int ee = 6;ee<9; ee++) {
      rval = moab.side_element(*pit,1,ee,edges[ee-3]); CHKERRQ_MOAB(rval);
    }
    // detect split edges
    BitRefEdges split_edges(0);
    EntityHandle edge_nodes[6];
    std::fill(&edge_nodes[0],&edge_nodes[6],no_handle);
    for(int ee = 0;ee<6;ee++) {
      RefEntity_multiIndex_view_by_parent_entity::iterator miit_view = ref_parent_ents_view.find(edges[ee]);
      if(miit_view != ref_parent_ents_view.end()) {
	if(((*miit_view)->getBitRefLevel()&bit).any()) {
	  edge_nodes[ee] = (*miit_view)->getRefEnt();
	  split_edges.set(ee);
	}
      }
    }
    if(split_edges.count()==0) {
      refinedEntities.modify(miit_prism,RefEntity_change_add_bit(bit));
      if(verb>6) PetscPrintf(comm,"no refinement");
      continue;
    }
    //check consistency
    if(verb>3) {
      std::ostringstream ss;
      ss << "prism split edges " << split_edges << " count " << split_edges.count() << std::endl;
      PetscPrintf(comm,ss.str().c_str());
    }
    // prism ref
    EntityHandle new_prism_conn[4*6];
    std::fill(&new_prism_conn[0],&new_prism_conn[4*6],no_handle);
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
	std::ostringstream ss;
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
    std::bitset<4> ref_prism_bit(0);
    ref_fe_by_composite::iterator miit_composite = ref_fe_by_comp.lower_bound(boost::make_tuple(*pit,split_edges.to_ulong()));
    ref_fe_by_composite::iterator hi_miit_composite = ref_fe_by_comp.upper_bound(boost::make_tuple(*pit,split_edges.to_ulong()));
    ref_fe_by_composite::iterator miit_composite2 = miit_composite;
    for(int pp = 0;miit_composite2!=hi_miit_composite;miit_composite2++,pp++) {
      //add this tet to this ref
      refinedEntities.modify(refinedEntities.find(miit_composite2->getRefEnt()),RefEntity_change_add_bit(bit));
      ref_prism_bit.set(pp,1);
      if(verb>2) {
	std::ostringstream ss;
	ss << "is refined " << *(miit_composite2->getRefElement()) << std::endl;
	PetscPrintf(comm,ss.str().c_str());
      }
    }
    if(miit_composite!=hi_miit_composite) {
      if(ref_prism_bit.count()!=(unsigned int)nb_new_prisms) SETERRQ(PETSC_COMM_SELF,1,"data inconsistency");
    } else {
      EntityHandle ref_prisms[4];
      // create prism
      for(int pp = 0;pp<nb_new_prisms;pp++) {
	if(verb>3) {
	  std::ostringstream ss;
	  ss << "ref prism " << ref_prism_bit << std::endl;
	  PetscPrintf(comm,ss.str().c_str());
	}
	if(!ref_prism_bit.test(pp)) {
	  rval = moab.create_element(MBPRISM,&new_prism_conn[6*pp],6,ref_prisms[pp]); CHKERRQ_MOAB(rval);
	  rval = moab.tag_set_data(th_RefParentHandle,&ref_prisms[pp],1,&*pit); CHKERRQ_MOAB(rval);
	  rval = moab.tag_set_data(th_RefBitLevel,&ref_prisms[pp],1,&bit); CHKERRQ_MOAB(rval);
	  rval = moab.tag_set_data(th_RefBitEdge,&ref_prisms[pp],1,&split_edges); CHKERRQ_MOAB(rval);
	  std::pair<RefEntity_multiIndex::iterator,bool> p_ent = refinedEntities.insert(
      boost::shared_ptr<RefEntity>(new RefEntity(basicEntityDataPtr,ref_prisms[pp]))
    );
	  std::pair<RefElement_multiIndex::iterator,bool> p_MoFEMFiniteElement;
	  try {
	    p_MoFEMFiniteElement = refinedFiniteElements.insert(ptrWrapperRefElement(
        boost::shared_ptr<RefElement>(new RefElement_PRISM(moab,*p_ent.first)))
      );
	  } catch (MoFEMException const &e) {
	    SETERRQ(PETSC_COMM_SELF,e.errorCode,e.errorMessage);
	  }
	  ref_prism_bit.set(pp);
	  ierr = addPrismToDatabase(ref_prisms[pp]); CHKERRQ(ierr);
	  if(verb>2) {
	    std::ostringstream ss;
	    ss << "add prism: " << *(p_MoFEMFiniteElement.first->getRefElement()) << std::endl;
      if(verb>7) {
        for(int nn = 0;nn<6;nn++) {
          ss << new_prism_conn[nn] << " ";
        }
        ss << std::endl;
      }
	    PetscPrintf(comm,ss.str().c_str());
	  }
	}
      }
    }
  }
  PetscFunctionReturn(0);
}
PetscErrorCode Core::refine_MESHSET(const EntityHandle meshset,const BitRefLevel &bit,const bool recursive,int verb) {
  PetscFunctionBegin;
  if(verb==-1) verb = verbose;
  typedef RefEntity_multiIndex::index<Ent_mi_tag>::type ref_ENTs_by_ent;
  ref_ENTs_by_ent::iterator miit = refinedEntities.find(meshset);
  if(miit==refinedEntities.end()) SETERRQ(PETSC_COMM_SELF,1,"this meshset is not in ref database");
  ierr = update_meshset_by_entities_children(meshset,bit,meshset,MBEDGE,recursive,verb); CHKERRQ(ierr);
  ierr = update_meshset_by_entities_children(meshset,bit,meshset,MBTRI,recursive,verb); CHKERRQ(ierr);
  ierr = update_meshset_by_entities_children(meshset,bit,meshset,MBTET,recursive,verb); CHKERRQ(ierr);
  refinedEntities.modify(miit,RefEntity_change_add_bit(bit));
  PetscFunctionReturn(0);
}

}
