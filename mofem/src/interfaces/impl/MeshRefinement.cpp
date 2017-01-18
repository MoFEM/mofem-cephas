/** \file MeshRefinementCore.cpp
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
#include <PrismInterface.hpp>
#include <SeriesRecorder.hpp>
#include <Core.hpp>

#include <MeshRefinement.hpp>
#include <EntityRefine.hpp>

namespace MoFEM {

PetscErrorCode MeshRefinement::queryInterface(const MOFEMuuid& uuid, UnknownInterface** iface) {
  PetscFunctionBegin;
  *iface = NULL;
  if(uuid == IDD_MOFEMMeshRefine) {
    *iface = dynamic_cast<MeshRefinement*>(this);
    PetscFunctionReturn(0);
  }
  if(uuid == IDD_MOFEMUnknown) {
    *iface = dynamic_cast<UnknownInterface*>(this);
    PetscFunctionReturn(0);
  }
  SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"unknown interface");
  PetscFunctionReturn(0);
}

MeshRefinement::MeshRefinement(const MoFEM::Core &core):
cOre(const_cast<MoFEM::Core&>(core)) {
}

PetscErrorCode MeshRefinement::add_verices_in_the_middel_of_edges(const EntityHandle meshset,const BitRefLevel &bit,const bool recursive,int verb) {
  PetscErrorCode ierr;
  MoABErrorCode rval;
  MoFEM::Interface &m_field = cOre;
  moab::Interface &moab = m_field.get_moab();
  PetscFunctionBegin;
  Range edges;
  rval = moab.get_entities_by_type(meshset,MBEDGE,edges,recursive);  CHKERRQ_MOAB(rval);
  if(edges.empty()) {
    Range tets;
    rval = moab.get_entities_by_type(meshset,MBTET,tets,recursive); CHKERRQ_MOAB(rval);
    rval = moab.get_adjacencies(tets,1,true,edges,moab::Interface::UNION); CHKERRQ_MOAB(rval);
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
PetscErrorCode MeshRefinement::add_verices_in_the_middel_of_edges(const Range &_edges,const BitRefLevel &bit,int verb) {
  PetscErrorCode ierr;
  MoABErrorCode rval;
  MoFEM::Interface &m_field = cOre;
  moab::Interface &moab = m_field.get_moab();
  const RefEntity_multiIndex *refined_ents_ptr;
  PetscFunctionBegin;
  ierr = m_field.get_ref_ents(&refined_ents_ptr); CHKERRQ(ierr);
  Range edges = _edges.subset_by_type(MBEDGE);
  typedef const RefEntity_multiIndex::index<Composite_EntType_and_ParentEntType_mi_tag>::type RefEntsByComposite;
  RefEntsByComposite &ref_ents = refined_ents_ptr->get<Composite_EntType_and_ParentEntType_mi_tag>();
  RefEntsByComposite::iterator miit = ref_ents.lower_bound(boost::make_tuple(MBVERTEX,MBEDGE));
  RefEntsByComposite::iterator hi_miit = ref_ents.upper_bound(boost::make_tuple(MBVERTEX,MBEDGE));
  RefEntity_multiIndex_view_by_parent_entity ref_parent_ents_view;
  for(;miit!=hi_miit;miit++) {
    std::pair<RefEntity_multiIndex_view_by_parent_entity::iterator,bool> p_ref_ent_view;
    p_ref_ent_view = ref_parent_ents_view.insert(*miit);
    if(!p_ref_ent_view.second) {
      SETERRQ(m_field.get_comm(),MOFEM_DATA_INCONSISTENCY,"non unique insertion");
    }
  }
  if(verb > 0) {
    std::ostringstream ss;
    ss << "ref level " << bit << " nb. edges to refine " << edges.size() << std::endl;
    PetscPrintf(m_field.get_comm(),ss.str().c_str());
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
        SETERRQ(m_field.get_comm(),MOFEM_DATA_INCONSISTENCY,"edge should have 2 edges");
      }
      double coords[num_nodes*3];
      rval = moab.get_coords(conn,num_nodes,coords); CHKERRQ_MOAB(rval);
      cblas_daxpy(3,1.,&coords[3],1,coords,1);
      cblas_dscal(3,0.5,coords,1);
      EntityHandle node;
      rval = moab.create_vertex(coords,node); CHKERRQ_MOAB(rval);
      rval = moab.tag_set_data(cOre.get_th_RefParentHandle(),&node,1,&*eit); CHKERRQ_MOAB(rval);
      rval = moab.tag_set_data(cOre.get_th_RefBitLevel(),&node,1,&bit); CHKERRQ_MOAB(rval);
      std::pair<RefEntity_multiIndex::iterator,bool> p_ent =
      const_cast<RefEntity_multiIndex*>(refined_ents_ptr)->insert(
        boost::shared_ptr<RefEntity>(new RefEntity(m_field.get_basic_entity_data_ptr(),node))
      );
      if(!p_ent.second) SETERRQ(m_field.get_comm(),MOFEM_DATA_INCONSISTENCY,"this entity is there");
      if(verb>2) {
        std::ostringstream ss;
        ss << *(p_ent.first) << std::endl;
        PetscPrintf(m_field.get_comm(),ss.str().c_str());
      }
    } else {
      const EntityHandle node = (*miit_view)->getRefEnt();
      if((*miit_view)->getEntType() != MBVERTEX) {
        SETERRQ(m_field.get_comm(),MOFEM_DATA_INCONSISTENCY,"child of edge should be vertex");
      }
      bool success = const_cast<RefEntity_multiIndex*>(refined_ents_ptr)->modify(
        refined_ents_ptr->get<Ent_mi_tag>().find(node),RefEntity_change_add_bit(bit)
      );
      if(!success) SETERRQ(m_field.get_comm(),MOFEM_DATA_INCONSISTENCY,"inconsistency in data");
    }
  }
  PetscFunctionReturn(0);
}
PetscErrorCode MeshRefinement::refine_TET(
  const EntityHandle meshset,const BitRefLevel &bit,const bool respect_interface,int verb
) {
  PetscErrorCode ierr;
  MoABErrorCode rval;
  MoFEM::Interface &m_field = cOre;
  moab::Interface &moab = m_field.get_moab();
  PetscFunctionBegin;
  Range tets;
  rval = moab.get_entities_by_type(meshset,MBTET,tets,false); CHKERRQ_MOAB(rval);
  ierr = refine_TET(tets,bit,respect_interface); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
PetscErrorCode MeshRefinement::refine_TET(
  const Range &_tets,const BitRefLevel &bit,const bool respect_interface,int verb
) {
  PetscErrorCode ierr;
  MoABErrorCode rval;
  MoFEM::Interface &m_field = cOre;
  moab::Interface &moab = m_field.get_moab();
  const RefEntity_multiIndex *refined_ents_ptr;
  const RefElement_multiIndex *refined_finite_elements_ptr;
  PetscFunctionBegin;

  ierr = m_field.get_ref_ents(&refined_ents_ptr); CHKERRQ(ierr);
  ierr = m_field.get_ref_finite_elements(&refined_finite_elements_ptr);

  //FIXME: refinement is based on entity handlers, should work on global ids of
  //nodes, this will allow parallelize algorithm in the future

  typedef const RefEntity_multiIndex::index<Ent_mi_tag>::type RefEntsByEnt;
  RefEntsByEnt &ref_ents_ent = refined_ents_ptr->get<Ent_mi_tag>();

  // Find all vertices which parent is edge
  typedef const RefEntity_multiIndex::index<Composite_EntType_and_ParentEntType_mi_tag>::type RefEntsByComposite;
  RefEntsByComposite &ref_ents = refined_ents_ptr->get<Composite_EntType_and_ParentEntType_mi_tag>();
  RefEntsByComposite::iterator miit = ref_ents.lower_bound(boost::make_tuple(MBVERTEX,MBEDGE));
  RefEntsByComposite::iterator hi_miit = ref_ents.upper_bound(boost::make_tuple(MBVERTEX,MBEDGE));
  RefEntity_multiIndex_view_by_parent_entity ref_parent_ents_view;
  for(;miit!=hi_miit;miit++) {
    std::pair<RefEntity_multiIndex_view_by_parent_entity::iterator,bool> p_ref_ent_view;
    p_ref_ent_view = ref_parent_ents_view.insert(*miit);
    if(!p_ref_ent_view.second) {
      SETERRQ(m_field.get_comm(),MOFEM_DATA_INCONSISTENCY,"non uniqe insertion");
    }
  }
  typedef const RefElement_multiIndex::index<Ent_mi_tag>::type RefElementByEnt;
  RefElementByEnt &ref_finite_element = refined_finite_elements_ptr->get<Ent_mi_tag>();
  typedef const RefElement_multiIndex::index<Composite_ParentEnt_And_BitsOfRefinedEdges_mi_tag>::type RefEntByComposite;
  RefEntByComposite &by_composite = refined_finite_elements_ptr->get<Composite_ParentEnt_And_BitsOfRefinedEdges_mi_tag>();

  if(respect_interface) {
    SETERRQ(m_field.get_comm(),MOFEM_DATA_INCONSISTENCY,"not implemented, set last parameter in refine_TET to false");
  }

  // make loop over all tets which going to be refined
  Range tets = _tets.subset_by_type(MBTET);
  Range::iterator tit = tets.begin();
  for(;tit!=tets.end();tit++) {

    RefElementByEnt::iterator miit2 = ref_finite_element.find(*tit);
    if(miit2==ref_finite_element.end()) {
      SETERRQ(m_field.get_comm(),MOFEM_DATA_INCONSISTENCY,"this tet is not in refinedFiniteElements");
    }

    // get tet connectivity
    const EntityHandle* conn;
    int num_nodes;
    moab.get_connectivity(*tit,conn,num_nodes,true);

    for(int nn = 0;nn<num_nodes;nn++) {
      bool success = const_cast<RefEntity_multiIndex*>(refined_ents_ptr)->modify(
        refined_ents_ptr->get<Ent_mi_tag>().find(conn[nn]),RefEntity_change_add_bit(bit)
      );
      if(!success) {
        SETERRQ(
          m_field.get_comm(),
          MOFEM_DATA_INCONSISTENCY,
          "can not set refinement bit level to tet node"
        );
      }
    }

    // get edges
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
              SETERRQ(m_field.get_comm(),MOFEM_DATA_INCONSISTENCY,"data inconsistency");
            }
            if(conn_edge[1] == edge_new_nodes[ee]) {
              SETERRQ(m_field.get_comm(),MOFEM_DATA_INCONSISTENCY,"data inconsistency");
            }
          }
          split_edges[parent_edges_bit.count()] = ee;
          parent_edges_bit.set(ee,1);
        }
      }
    }
    // test if nodes used to refine are not part of tet
    for(int ee = 0;ee<6;ee++) {
      if(edge_new_nodes[ee] == no_handle) continue;
      for(int nn = 0;nn<4;nn++) {
        if(conn[nn] == edge_new_nodes[ee]) {
          // std::cerr << conn[0] << " "
          // << conn[1] << " "
          // << conn[2] << " "
          // << conn[3] << " : "
          // << edge_new_nodes[ee] << std::endl;
          for(int eee = 0;eee<6;eee++) {
            EntityHandle edge;
            rval = moab.side_element(*tit,1,eee,edge);  CHKERRQ_MOAB(rval);
            const EntityHandle* conn_edge;
            int num_nodes;
            moab.get_connectivity(edge,conn_edge,num_nodes,true);
            std::cerr << conn_edge[0] << " " << conn_edge[1] << std::endl;
          }
          SETERRQ(m_field.get_comm(),MOFEM_DATA_INCONSISTENCY,"data inconsistency");
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
        RefEntsByEnt::iterator tit_miit;
        tit_miit = ref_ents_ent.find(*tit);
        if(tit_miit==ref_ents_ent.end()) SETERRQ(PETSC_COMM_SELF,1,"data inconsistency");
        bool success = const_cast<RefEntity_multiIndex*>(refined_ents_ptr)->modify(
          tit_miit,RefEntity_change_add_bit(bit)
        );
        if(!success) SETERRQ(m_field.get_comm(),MOFEM_DATA_INCONSISTENCY,"impossible tet");
        Range tit_conn;
        rval = moab.get_connectivity(&*tit,1,tit_conn,true); CHKERRQ_MOAB(rval);
        for(Range::iterator nit = tit_conn.begin();nit!=tit_conn.end();nit++) {
          RefEntsByEnt::iterator nit_miit = ref_ents_ent.find(*nit);
          if(nit_miit==ref_ents_ent.end()) {
            SETERRQ(m_field.get_comm(),
            MOFEM_DATA_INCONSISTENCY,"can not find face in refinedEntities");
          }
          bool success = const_cast<RefEntity_multiIndex*>(refined_ents_ptr)->modify(
            nit_miit,RefEntity_change_add_bit(bit)
          );
          if(!success) {
            SETERRQ(m_field.get_comm(),MOFEM_DATA_INCONSISTENCY,"impossible node");
          }
        }
        Range tit_edges;
        rval = moab.get_adjacencies(&*tit,1,1,false,tit_edges); CHKERRQ_MOAB(rval);
        for(Range::iterator eit = tit_edges.begin();eit!=tit_edges.end();eit++) {
          RefEntsByEnt::iterator eit_miit = ref_ents_ent.find(*eit);
          if(eit_miit==ref_ents_ent.end()) {
            SETERRQ(
              m_field.get_comm(),
              MOFEM_DATA_INCONSISTENCY,"can not find face in refinedEntities"
            );
          }
          bool success = const_cast<RefEntity_multiIndex*>(refined_ents_ptr)->modify(
            eit_miit,RefEntity_change_add_bit(bit)
          );
          if(!success) SETERRQ(m_field.get_comm(),MOFEM_DATA_INCONSISTENCY,"impossible edge");
        }
        Range tit_faces;
        rval = moab.get_adjacencies(&*tit,1,2,false,tit_faces); CHKERRQ_MOAB(rval);
        if(tit_faces.size()!=4) {
          SETERRQ(
            m_field.get_comm(),
            MOFEM_DATA_INCONSISTENCY,
            "existing tet in mofem database should have 4 adjacent edges"
          );
        }
        for(Range::iterator fit = tit_faces.begin();fit!=tit_faces.end();fit++) {
          RefEntsByEnt::iterator fit_miit = ref_ents_ent.find(*fit);
          if(fit_miit==ref_ents_ent.end()) {
            SETERRQ(m_field.get_comm(),MOFEM_DATA_INCONSISTENCY,"can not find face in refinedEntities");
          }
          bool success = const_cast<RefEntity_multiIndex*>(refined_ents_ptr)->modify(
            fit_miit,RefEntity_change_add_bit(bit)
          );
          if(!success) SETERRQ(m_field.get_comm(),MOFEM_DATA_INCONSISTENCY,"impossible face");
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
    RefEntByComposite::iterator miit_composite = by_composite.lower_bound(
      boost::make_tuple(*tit,parent_edges_bit.to_ulong())
    );
    RefEntByComposite::iterator hi_miit_composite = by_composite.upper_bound(
      boost::make_tuple(*tit,parent_edges_bit.to_ulong())
    );
    // check if tet with this refinement shame already exits
    if(distance(miit_composite,hi_miit_composite)==(unsigned int)nb_new_tets) {
      for(int tt = 0;miit_composite!=hi_miit_composite;miit_composite++,tt++) {
        EntityHandle tet = miit_composite->getRefEnt();
        //set ref tets entities
        ref_tets[tt] = tet;
        ref_tets_bit.set(tt,1);
        //add this tet if exist to this ref level
        RefEntity_multiIndex::iterator ref_tet_it;
        ref_tet_it = refined_ents_ptr->find(tet);
        if(ref_tet_it == refined_ents_ptr->end()) {
          SETERRQ(m_field.get_comm(),MOFEM_DATA_INCONSISTENCY,"data inconsistency");
        }
        bool success = const_cast<RefEntity_multiIndex*>(refined_ents_ptr)->modify(
          ref_tet_it,RefEntity_change_add_bit(bit)
        );
        if(!success) {
          SETERRQ(m_field.get_comm(),MOFEM_DATA_INCONSISTENCY,"modification unsuccessful");
        }
        // change bit for adjacent entities of existing tet
        Range adj_ents;
        rval = moab.get_adjacencies(&tet,1,1,false,adj_ents); CHKERRQ_MOAB(rval);
        rval = moab.get_adjacencies(&tet,1,2,false,adj_ents); CHKERRQ_MOAB(rval);
        for(Range::iterator ait = adj_ents.begin();ait!=adj_ents.end();ait++) {
          RefEntity_multiIndex::iterator ref_it;
          ref_it = refined_ents_ptr->find(*ait);
          if(ref_it == refined_ents_ptr->end()) {
            SETERRQ(m_field.get_comm(),MOFEM_DATA_INCONSISTENCY,"data inconsistency");
          }
          bool success = const_cast<RefEntity_multiIndex*>(refined_ents_ptr)->modify(
            ref_it,RefEntity_change_add_bit(bit)
          );
          if(!success) {
            SETERRQ(m_field.get_comm(),MOFEM_DATA_INCONSISTENCY,"modification unsuccessful");
          }
        }
        //verbose
        if(verb>2) {
          std::ostringstream ss;
          ss << miit_composite->getRefElement() << std::endl;
          PetscPrintf(m_field.get_comm(),ss.str().c_str());
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
                SETERRQ(m_field.get_comm(),MOFEM_DATA_INCONSISTENCY,"data inconsistency");
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
            SETERRQ(m_field.get_comm(),MOFEM_DATA_INCONSISTENCY,"tetrahedral should have 4 nodes");
          }
          int ref_type[2];
          ref_type[0] = parent_edges_bit.count();
          ref_type[1] = sub_type;
          rval = moab.tag_set_data(cOre.get_th_RefType(),&ref_tets[tt],1,ref_type); CHKERRQ_MOAB(rval);
          rval = moab.tag_set_data(cOre.get_th_RefParentHandle(),&ref_tets[tt],1,&*tit); CHKERRQ_MOAB(rval);
          rval = moab.tag_set_data(cOre.get_th_RefBitLevel(),&ref_tets[tt],1,&bit); CHKERRQ_MOAB(rval);
          rval = moab.tag_set_data(cOre.get_th_RefBitEdge(),&ref_tets[tt],1,&parent_edges_bit); CHKERRQ_MOAB(rval);
          //add refined entity
          std::pair<RefEntity_multiIndex::iterator,bool> p_MoFEMEntity =
          const_cast<RefEntity_multiIndex*>(refined_ents_ptr)->insert(
            boost::shared_ptr<RefEntity>(new RefEntity(
              m_field.get_basic_entity_data_ptr(),ref_tets[tt]
            ))
          );
          //add refined element
          std::pair<RefElement_multiIndex::iterator,bool> p_MoFEMFiniteElement;
          try {
            p_MoFEMFiniteElement =  const_cast<RefElement_multiIndex*>(refined_finite_elements_ptr)->
            insert(ptrWrapperRefElement(
              boost::shared_ptr<RefElement>(new RefElement_TET(*p_MoFEMEntity.first)))
            );
          } catch (MoFEMException const &e) {
            SETERRQ(m_field.get_comm(),e.errorCode,e.errorMessage);
          }
          //set bit that this element is now in databse
          ref_tets_bit.set(tt);
          if(verb>2) {
            std::ostringstream ss;
            ss << "add tet: " << *(p_MoFEMFiniteElement.first->getRefElement()) << std::endl;
            PetscPrintf(m_field.get_comm(),ss.str().c_str());
          }
        }
      }
    }
    // //debug
    // miit_composite = by_composite.lower_bound(boost::make_tuple(*tit,parent_edges_bit.to_ulong()));
    // hi_miit_composite = by_composite.upper_bound(boost::make_tuple(*tit,parent_edges_bit.to_ulong()));
    // if(miit_composite!=hi_miit_composite) {
    //   if(distance(miit_composite,hi_miit_composite)!=(unsigned int)nb_new_tets) {
    //     SETERRQ2(m_field.get_comm(),1,"data inconsistency %u != %u",
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
      // Get edges on face and loop over those edges to add mid-nodes to range
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
    //tet_nodes contains all nodes on tet including mid edge nodes
    Range tet_nodes;
    rval = moab.get_connectivity(&*tit,1,tet_nodes,true); CHKERRQ_MOAB(rval);
    for(
      std::map<EntityHandle,const RefEntity*>::iterator
      map_miit = map_ref_nodes_by_edges.begin();
      map_miit != map_ref_nodes_by_edges.end();map_miit++
    ) {
      tet_nodes.insert(map_miit->second->getRefEnt());
    }
    Range ref_edges;
    // Get all all edges of refined tets
    rval = moab.get_adjacencies(
      ref_tets,nb_new_tets,1,true,ref_edges,moab::Interface::UNION
    ); CHKERRQ_MOAB(rval);
    // Check for all ref edge and set parents
    for(Range::iterator reit = ref_edges.begin();reit!=ref_edges.end();reit++) {
      Range ref_edges_nodes;
      rval = moab.get_connectivity(&*reit,1,ref_edges_nodes,true); CHKERRQ_MOAB(rval);
      if(ref_edges_nodes.size()!=2) {
        SETERRQ(
          m_field.get_comm(),MOFEM_DATA_INCONSISTENCY,
          "data inconsistency, edge should have 2 nodes"
        );
      }
      // Check if ref edge is an coarse edge (loop over coarse tet edges)
      int ee = 0;
      for(;ee<6;ee++) {
        // Two nodes are common (node[0],node[1],ref_node (if exist))
        // this tests if given edge is contained by edge of refined tetrahedral
        if(intersect(edges_nodes[ee],ref_edges_nodes).size()==2) {
          EntityHandle edge = tit_edges[ee];
          std::pair<RefEntity_multiIndex::iterator,bool> p_ent =
          const_cast<RefEntity_multiIndex*>(refined_ents_ptr)->insert(
            boost::shared_ptr<RefEntity>(new RefEntity(
              m_field.get_basic_entity_data_ptr(),*reit
            ))
          );
          if(p_ent.second) {
            bool success = const_cast<RefEntity_multiIndex*>(refined_ents_ptr)->
            modify(p_ent.first,RefEntity_change_parent(edge));
            if(!success) {
              SETERRQ(
                m_field.get_comm(),
                MOFEM_DATA_INCONSISTENCY,"impossible to set edge parent"
              );
            }
          }
          bool success = const_cast<RefEntity_multiIndex*>(refined_ents_ptr)->
          modify(p_ent.first,RefEntity_change_add_bit(bit));
          if(!success) {
            SETERRQ(
              m_field.get_comm(),
              MOFEM_DATA_INCONSISTENCY,"impossible to set edge bit"
            );
          }
          if(p_ent.second) {
            if(verb>2) {
              std::ostringstream ss;
              ss << "edge parent: " << *(p_ent.first) << std::endl;
              PetscPrintf(m_field.get_comm(),ss.str().c_str());
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
        //this tests if given edge is contained by face of  tetrahedral
        if(intersect(faces_nodes[ff],ref_edges_nodes).size()==2) {
          EntityHandle face = tit_faces[ff];
          //add edge to refinedEntities
          std::pair<RefEntity_multiIndex::iterator,bool> p_ent =
          const_cast<RefEntity_multiIndex*>(refined_ents_ptr)->insert(
            boost::shared_ptr<RefEntity>(new RefEntity(m_field.get_basic_entity_data_ptr(),*reit))
          );
          if(p_ent.second) {
            bool success = const_cast<RefEntity_multiIndex*>(refined_ents_ptr)->
            modify(p_ent.first,RefEntity_change_parent(face));
            if(!success) {
              SETERRQ(
                m_field.get_comm(),
                MOFEM_DATA_INCONSISTENCY,"impossible to set edge parent"
              );
            }
          }
          bool success = const_cast<RefEntity_multiIndex*>(refined_ents_ptr)->
          modify(p_ent.first,RefEntity_change_add_bit(bit));
          if(!success) {
            SETERRQ(m_field.get_comm(),MOFEM_DATA_INCONSISTENCY,"impossible to set edge bit");
          }
          if(p_ent.second) {
            if(verb>2) {
              std::ostringstream ss;
              ss << "face parent: " << *(p_ent.first) << std::endl;
              PetscPrintf(m_field.get_comm(),ss.str().c_str());
            }
          }
          break;
        }
      }
      if(ff<4) continue; //this refined edge is contained by face of tetrahedral
      // check if ref edge is in coarse tetrahedral (i.e. that is internal edge of refined tetrahedral)
      if(intersect(tet_nodes,ref_edges_nodes).size()==2) {
        rval = moab.tag_set_data(cOre.get_th_RefParentHandle(),&*reit,1,&*tit); CHKERRQ_MOAB(rval);
        //add edge to refinedEntities
        std::pair<RefEntity_multiIndex::iterator,bool> p_ent =
        const_cast<RefEntity_multiIndex*>(refined_ents_ptr)->insert(
          boost::shared_ptr<RefEntity>(new RefEntity(m_field.get_basic_entity_data_ptr(),*reit))
        );
        bool success = const_cast<RefEntity_multiIndex*>(refined_ents_ptr)->
        modify(p_ent.first,RefEntity_change_add_bit(bit));
        if(!success) {
          SETERRQ(m_field.get_comm(),MOFEM_DATA_INCONSISTENCY,"impossible to set edge parent");
        }
        if(p_ent.second) {
          if(verb>2) {
            std::ostringstream ss;
            ss << "tet parent: " << *(p_ent.first) << std::endl;
            PetscPrintf(m_field.get_comm(),ss.str().c_str());
          }
        }
        continue;
      }

      // This will help to debug error
      Range ref_edges_nodes_tets;
      rval = moab.get_adjacencies(ref_edges_nodes,3,true,ref_edges_nodes_tets,moab::Interface::UNION); CHKERRQ_MOAB(rval);
      ref_edges_nodes_tets = intersect(ref_edges_nodes_tets,tets);
      ref_edges_nodes_tets.insert(*reit);

      EntityHandle meshset_out;
      rval = moab.create_meshset(MESHSET_SET|MESHSET_TRACK_OWNER,meshset_out); CHKERRQ_MOAB(rval);
      rval = moab.add_entities(meshset_out,ref_edges_nodes_tets); CHKERRQ_MOAB(rval);
      rval = moab.write_file("debug_error.vtk","VTK","",&meshset_out,1); CHKERRQ_MOAB(rval);
      rval = moab.delete_entities(&meshset_out,1); CHKERRQ_MOAB(rval);

      // Refined edge is not child of any edge, face or tetrahedral, this is imposible edge
      SETERRQ(m_field.get_comm(),MOFEM_DATA_INCONSISTENCY,"impossible refined edge");
    }

    Range ref_faces;
    rval = moab.get_adjacencies(
      ref_tets,nb_new_tets,2,true,ref_faces,moab::Interface::UNION
    ); CHKERRQ_MOAB(rval);
    Tag th_interface_side;
    const int def_side[] = {0};
    rval = moab.tag_get_handle(
      "INTERFACE_SIDE",1,MB_TYPE_INTEGER,
      th_interface_side,MB_TAG_CREAT|MB_TAG_SPARSE,def_side
    ); CHKERRQ_MOAB(rval);
    // Check for all ref faces
    for(Range::iterator rfit = ref_faces.begin();rfit!=ref_faces.end();rfit++) {
      Range ref_faces_nodes;
      rval = moab.get_connectivity(&*rfit,1,ref_faces_nodes,true); CHKERRQ_MOAB(rval);
      // Check if ref face is in coarse face
      int ff = 0;
      for(;ff<4;ff++) {
        // Check if refined triangle is contained by face of tetrahedral
        if(intersect(faces_nodes[ff],ref_faces_nodes).size()==3) {
          EntityHandle face = tit_faces[ff];
          rval = moab.tag_set_data(cOre.get_th_RefParentHandle(),&*rfit,1,&face); CHKERRQ_MOAB(rval);
          int side = 0;
          // Set face side if it is on interface
          rval = moab.tag_get_data(th_interface_side,&face,1,&side); CHKERRQ_MOAB(rval);
          rval = moab.tag_set_data(th_interface_side,&*rfit,1,&side); CHKERRQ_MOAB(rval);
          // Add face to refinedEntities
          std::pair<RefEntity_multiIndex::iterator,bool> p_ent =
          const_cast<RefEntity_multiIndex*>(refined_ents_ptr)->insert(
            boost::shared_ptr<RefEntity>(new RefEntity(m_field.get_basic_entity_data_ptr(),*rfit))
          );
          bool success = const_cast<RefEntity_multiIndex*>(refined_ents_ptr)->
          modify(p_ent.first,RefEntity_change_add_bit(bit));
          if(!success) {
            SETERRQ(m_field.get_comm(),MOFEM_DATA_INCONSISTENCY,"impossible to set face bit level");
          }
          if(p_ent.second) {
            if(verb>2) {
              std::ostringstream ss;
              ss << "face: " << *(p_ent.first) << std::endl;
              PetscPrintf(m_field.get_comm(),ss.str().c_str());
            }
          }
          break;
        }
      }
      if(ff<4) continue; //this face is contained by one of tetrahedrons
      //check if ref face is in coarse tetrahedral
      //this is ref face which is contained by tetrahedral volume
      if(intersect(tet_nodes,ref_faces_nodes).size()==3) {
        rval = moab.tag_set_data(cOre.get_th_RefParentHandle(),&*rfit,1,&*tit); CHKERRQ_MOAB(rval);
        std::pair<RefEntity_multiIndex::iterator,bool> p_ent =
        const_cast<RefEntity_multiIndex*>(refined_ents_ptr)->insert(
          boost::shared_ptr<RefEntity>(new RefEntity(m_field.get_basic_entity_data_ptr(),*rfit))
        );
        //add face to refinedEntities
        bool success = const_cast<RefEntity_multiIndex*>(refined_ents_ptr)->
        modify(p_ent.first,RefEntity_change_add_bit(bit));
        if(!success) SETERRQ(m_field.get_comm(),MOFEM_DATA_INCONSISTENCY,"impossible to set face parent");
        if(p_ent.second) {
          if(verb>2) {
            std::ostringstream ss;
            ss << "tet parent: " << *(p_ent.first) << std::endl;
            PetscPrintf(m_field.get_comm(),ss.str().c_str());
          }
        }
        continue;
      }
      SETERRQ(m_field.get_comm(),MOFEM_DATA_INCONSISTENCY,"impossible refined face");
    }
  }
  PetscFunctionReturn(0);
}
PetscErrorCode MeshRefinement::refine_PRISM(const EntityHandle meshset,const BitRefLevel &bit,int verb) {
  PetscErrorCode ierr;
  MoABErrorCode rval;
  MoFEM::Interface &m_field = cOre;
  moab::Interface &moab = m_field.get_moab();
  const RefEntity_multiIndex *refined_ents_ptr;
  const RefElement_multiIndex *refined_finite_elements_ptr;

  //FIXME: refinement is based on entity handlers, should work on global ids of nodes, this will allow parallelize algorithm in the future

  PetscFunctionBegin;

  ierr = m_field.get_ref_ents(&refined_ents_ptr); CHKERRQ(ierr);
  ierr = m_field.get_ref_finite_elements(&refined_finite_elements_ptr);

  typedef const RefEntity_multiIndex::index<Ent_mi_tag>::type RefEntsByEnt;
  typedef const RefElement_multiIndex::index<Composite_ParentEnt_And_BitsOfRefinedEdges_mi_tag>::type RefFeByComposite;
  RefFeByComposite &ref_fe_by_comp = refined_finite_elements_ptr->get<Composite_ParentEnt_And_BitsOfRefinedEdges_mi_tag>();
  //find all vertices which parent is edge
  typedef const RefEntity_multiIndex::index<Composite_EntType_and_ParentEntType_mi_tag>::type RefEntsByComposite;
  RefEntsByComposite &ref_ents_by_comp = refined_ents_ptr->get<Composite_EntType_and_ParentEntType_mi_tag>();
  RefEntsByComposite::iterator miit = ref_ents_by_comp.lower_bound(boost::make_tuple(MBVERTEX,MBEDGE));
  RefEntsByComposite::iterator hi_miit = ref_ents_by_comp.upper_bound(boost::make_tuple(MBVERTEX,MBEDGE));
  RefEntity_multiIndex_view_by_parent_entity ref_parent_ents_view;
  for(;miit!=hi_miit;miit++) {
    std::pair<RefEntity_multiIndex_view_by_parent_entity::iterator,bool> p_ref_ent_view;
    p_ref_ent_view = ref_parent_ents_view.insert(*miit);
    if(!p_ref_ent_view.second) {
      SETERRQ(m_field.get_comm(),MOFEM_DATA_INCONSISTENCY,"non uniqe insertion");
    }
  }
  Range prisms;
  rval = moab.get_entities_by_type(meshset,MBPRISM,prisms,false); CHKERRQ_MOAB(rval);
  Range::iterator pit = prisms.begin();
  for(;pit!=prisms.end();pit++) {
    RefEntsByEnt::iterator miit_prism = refined_ents_ptr->get<Ent_mi_tag>().find(*pit);
    if(miit_prism==refined_ents_ptr->end()) {
      SETERRQ(m_field.get_comm(),MOFEM_DATA_INCONSISTENCY,"this prism is not in ref database");
    }
    if(verb>3) {
      std::ostringstream ss;
      ss << "ref prism " << *miit << std::endl;
      PetscPrintf(m_field.get_comm(),ss.str().c_str());
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
      const_cast<RefEntity_multiIndex*>(refined_ents_ptr)->
      modify(miit_prism,RefEntity_change_add_bit(bit));
      if(verb>6) PetscPrintf(m_field.get_comm(),"no refinement");
      continue;
    }
    //check consistency
    if(verb>3) {
      std::ostringstream ss;
      ss << "prism split edges " << split_edges << " count " << split_edges.count() << std::endl;
      PetscPrintf(m_field.get_comm(),ss.str().c_str());
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
      SETERRQ(m_field.get_comm(),MOFEM_DATA_INCONSISTENCY,ss.str().c_str());
    }
    // find that prism
    std::bitset<4> ref_prism_bit(0);
    RefFeByComposite::iterator miit_composite = ref_fe_by_comp.lower_bound(boost::make_tuple(*pit,split_edges.to_ulong()));
    RefFeByComposite::iterator hi_miit_composite = ref_fe_by_comp.upper_bound(boost::make_tuple(*pit,split_edges.to_ulong()));
    RefFeByComposite::iterator miit_composite2 = miit_composite;
    for(int pp = 0;miit_composite2!=hi_miit_composite;miit_composite2++,pp++) {
      //add this tet to this ref
      const_cast<RefEntity_multiIndex*>(refined_ents_ptr)->
      modify(refined_ents_ptr->find(miit_composite2->getRefEnt()),RefEntity_change_add_bit(bit));
      ref_prism_bit.set(pp,1);
      if(verb>2) {
        std::ostringstream ss;
        ss << "is refined " << *(miit_composite2->getRefElement()) << std::endl;
        PetscPrintf(m_field.get_comm(),ss.str().c_str());
      }
    }
    if(miit_composite!=hi_miit_composite) {
      if(ref_prism_bit.count()!=(unsigned int)nb_new_prisms) SETERRQ(m_field.get_comm(),MOFEM_DATA_INCONSISTENCY,"data inconsistency");
    } else {
      EntityHandle ref_prisms[4];
      // create prism
      for(int pp = 0;pp<nb_new_prisms;pp++) {
        if(verb>3) {
          std::ostringstream ss;
          ss << "ref prism " << ref_prism_bit << std::endl;
          PetscPrintf(m_field.get_comm(),ss.str().c_str());
        }
        if(!ref_prism_bit.test(pp)) {
          rval = moab.create_element(MBPRISM,&new_prism_conn[6*pp],6,ref_prisms[pp]); CHKERRQ_MOAB(rval);
          rval = moab.tag_set_data(cOre.get_th_RefParentHandle(),&ref_prisms[pp],1,&*pit); CHKERRQ_MOAB(rval);
          rval = moab.tag_set_data(cOre.get_th_RefBitLevel(),&ref_prisms[pp],1,&bit); CHKERRQ_MOAB(rval);
          rval = moab.tag_set_data(cOre.get_th_RefBitEdge(),&ref_prisms[pp],1,&split_edges); CHKERRQ_MOAB(rval);
          std::pair<RefEntity_multiIndex::iterator,bool> p_ent =
          const_cast<RefEntity_multiIndex*>(refined_ents_ptr)->insert(
            boost::shared_ptr<RefEntity>(new RefEntity(m_field.get_basic_entity_data_ptr(),ref_prisms[pp]))
          );
          std::pair<RefElement_multiIndex::iterator,bool> p_fe;
          try {
            p_fe = const_cast<RefElement_multiIndex*>(refined_finite_elements_ptr)->
            insert(ptrWrapperRefElement(
              boost::shared_ptr<RefElement>(new RefElement_PRISM(*p_ent.first)))
            );
          } catch (MoFEMException const &e) {
            SETERRQ(PETSC_COMM_SELF,e.errorCode,e.errorMessage);
          }
          ref_prism_bit.set(pp);
          ierr = cOre.addPrismToDatabase(ref_prisms[pp]); CHKERRQ(ierr);
          if(verb>2) {
            std::ostringstream ss;
            ss << "add prism: " << *(p_fe.first->getRefElement()) << std::endl;
            if(verb>7) {
              for(int nn = 0;nn<6;nn++) {
                ss << new_prism_conn[nn] << " ";
              }
              ss << std::endl;
            }
            PetscPrintf(m_field.get_comm(),ss.str().c_str());
          }
        }
      }
    }
  }
  PetscFunctionReturn(0);
}
PetscErrorCode MeshRefinement::refine_MESHSET(
  const EntityHandle meshset,const BitRefLevel &bit,const bool recursive,int verb
) {
  PetscErrorCode ierr;
  MoABErrorCode rval;
  MoFEM::Interface &m_field = cOre;
  moab::Interface &moab = m_field.get_moab();
  const RefEntity_multiIndex *refined_ents_ptr;
  PetscFunctionBegin;
  ierr = m_field.get_ref_ents(&refined_ents_ptr); CHKERRQ(ierr);
  typedef const RefEntity_multiIndex::index<Ent_mi_tag>::type RefEntsByEnt;
  RefEntsByEnt::iterator miit = refined_ents_ptr->find(meshset);
  if(miit==refined_ents_ptr->end()) {
    SETERRQ(m_field.get_comm(),MOFEM_DATA_INCONSISTENCY,"this meshset is not in ref database");
  }
  ierr = m_field.update_meshset_by_entities_children(meshset,bit,meshset,MBEDGE,recursive,verb); CHKERRQ(ierr);
  ierr = m_field.update_meshset_by_entities_children(meshset,bit,meshset,MBTRI,recursive,verb); CHKERRQ(ierr);
  ierr = m_field.update_meshset_by_entities_children(meshset,bit,meshset,MBTET,recursive,verb); CHKERRQ(ierr);
  const_cast<RefEntity_multiIndex*>(refined_ents_ptr)->modify(miit,RefEntity_change_add_bit(bit));
  PetscFunctionReturn(0);
}

}
