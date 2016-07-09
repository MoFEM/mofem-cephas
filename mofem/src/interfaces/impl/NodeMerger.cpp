/** \file NodeMerger.cpp
 * \brief Interface for merging nodes
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
#include <version.h>
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
#include <Interface.hpp>
#include <MeshRefinment.hpp>
#include <PrismInterface.hpp>
#include <SeriesRecorder.hpp>
#include <Core.hpp>

#include <NodeMerger.hpp>

namespace MoFEM {

PetscErrorCode NodeMergerInterface::queryInterface(const MOFEMuuid& uuid, UnknownInterface** iface) {
  PetscFunctionBegin;
  *iface = NULL;
  if(uuid == IDD_MOFEMNodeMerger) {
    *iface = dynamic_cast<NodeMergerInterface*>(this);
    PetscFunctionReturn(0);
  }
  if(uuid == IDD_MOFEMUnknown) {
    *iface = dynamic_cast<UnknownInterface*>(this);
    PetscFunctionReturn(0);
  }
  SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"unknown interface");
  PetscFunctionReturn(0);
}

PetscErrorCode NodeMergerInterface::mergeNodes(EntityHandle father,EntityHandle mother,BitRefLevel bit,Range *tets_ptr) {
  PetscFunctionBegin;

  MoFEM::Interface& m_field = cOre;
  PetscErrorCode ierr;
  ErrorCode rval;

  Range father_edges;
  rval = m_field.get_moab().get_adjacencies(&father,1,1,false,father_edges); CHKERRQ_MOAB(rval);
  Range mother_edges;
  rval = m_field.get_moab().get_adjacencies(&mother,1,1,false,mother_edges); CHKERRQ_MOAB(rval);
  Range common_edge;
  common_edge = intersect(father_edges,mother_edges);
  if(tets_ptr != NULL) {
    Range tets_edges;
    rval = m_field.get_moab().get_adjacencies(
      *tets_ptr,1,false,tets_edges,moab::Interface::UNION
    ); CHKERRQ_MOAB(rval);
    common_edge = intersect(common_edge,tets_edges);
    father_edges = intersect(father_edges,tets_edges);
    mother_edges = intersect(mother_edges,tets_edges);
  }
  if(common_edge.empty()) {
    SETERRQ(PETSC_COMM_SELF,1,"no common edge between nodes");
  }

  Range father_tets;
  rval = m_field.get_moab().get_adjacencies(&father,1,3,false,father_tets); CHKERRQ_MOAB(rval);
  Range mother_tets;
  rval = m_field.get_moab().get_adjacencies(&mother,1,3,false,mother_tets); CHKERRQ_MOAB(rval);
  Range edge_tets;
  rval = m_field.get_moab().get_adjacencies(common_edge,3,true,edge_tets); CHKERRQ_MOAB(rval);
  mother_tets = subtract(mother_tets,edge_tets);
  if(tets_ptr!=NULL) {
    father_tets = intersect(father_tets,*tets_ptr);
    mother_tets = intersect(mother_tets,*tets_ptr);
    edge_tets = intersect(edge_tets,*tets_ptr);
  }

  Range created_tets;
  for(
    Range::iterator tit = mother_tets.begin();
    tit!=mother_tets.end();tit++
  ) {
    const EntityHandle* conn;
    int num_nodes;
    rval = m_field.get_moab().get_connectivity(*tit,conn,num_nodes,true); CHKERRQ_MOAB(rval);
    EntityHandle new_conn[4];
    for(int nn = 0;nn<4;nn++) {
      if(conn[nn] == mother) {
        new_conn[nn] = father;
      } else {
        new_conn[nn] = conn[nn];
      }
    }
    EntityHandle tet;
    rval = m_field.get_moab().create_element(MBTET,new_conn,4,tet); CHKERRQ_MOAB(rval);
    rval = m_field.get_moab().tag_set_data(cOre.get_th_RefParentHandle(),&tet,1,&*tit); CHKERRQ_MOAB(rval);
    created_tets.insert(tet);
  }

  Range adj_ents;
  rval = m_field.get_moab().get_adjacencies(father_tets,1,false,adj_ents,moab::Interface::UNION); CHKERRQ_MOAB(rval);
  rval = m_field.get_moab().get_adjacencies(father_tets,2,false,adj_ents,moab::Interface::UNION); CHKERRQ_MOAB(rval);
  for(Range::iterator eit = adj_ents.begin();eit!=adj_ents.end();eit++) {
    const EntityHandle* conn;
    int num_nodes;
    rval = m_field.get_moab().get_connectivity(*eit,conn,num_nodes,true); CHKERRQ_MOAB(rval);
    EntityHandle new_conn[num_nodes];
    int nb_new_node = 0;
    int nn = 0;
    for(;nn<num_nodes;nn++) {
      if(conn[nn] == mother) {
        nb_new_node = 0;
        break;
      } else if(conn[nn] == father) {
        new_conn[nn] = mother;
        nb_new_node++;
      } else {
        new_conn[nn] = conn[nn];
      }
    }
    if(nb_new_node > 0) {
      int dim = m_field.get_moab().dimension_from_handle(*eit);
      Range new_ent;
      rval = m_field.get_moab().get_adjacencies(new_conn,num_nodes,dim,true,new_ent); CHKERRQ_MOAB(rval);
      if(new_ent.empty()) continue;
      if(new_ent.size()!=1) {
        SETERRQ1(PETSC_COMM_SELF,1,"data inconsistency %u",new_ent.size());
      }
      rval = m_field.get_moab().tag_set_data(cOre.get_th_RefParentHandle(),&*eit,1,&*new_ent.begin()); CHKERRQ_MOAB(rval);
    }
  }
  adj_ents.clear();
  rval = m_field.get_moab().get_adjacencies(mother_tets,1,false,adj_ents,moab::Interface::UNION); CHKERRQ_MOAB(rval);
  rval = m_field.get_moab().get_adjacencies(mother_tets,2,false,adj_ents,moab::Interface::UNION); CHKERRQ_MOAB(rval);
  rval = m_field.get_moab().get_adjacencies(edge_tets,1,false,adj_ents,moab::Interface::UNION); CHKERRQ_MOAB(rval);
  rval = m_field.get_moab().get_adjacencies(edge_tets,2,false,adj_ents,moab::Interface::UNION); CHKERRQ_MOAB(rval);
  for(Range::iterator eit = adj_ents.begin();eit!=adj_ents.end();eit++) {
    const EntityHandle* conn;
    int num_nodes;
    rval = m_field.get_moab().get_connectivity(*eit,conn,num_nodes,true); CHKERRQ_MOAB(rval);
    EntityHandle new_conn[num_nodes];
    int nb_new_node = 0;
    int nn = 0;
    for(;nn<num_nodes;nn++) {
      if(conn[nn] == father) {
        nb_new_node = 0;
        break;
      } else if(conn[nn] == mother) {
        new_conn[nn] = father;
        nb_new_node++;
      } else {
        new_conn[nn] = conn[nn];
      }
    }
    if(nb_new_node > 0) {
      int dim = m_field.get_moab().dimension_from_handle(*eit);
      Range new_ent;
      rval = m_field.get_moab().get_adjacencies(new_conn,num_nodes,dim,true,new_ent); CHKERRQ_MOAB(rval);
      if(new_ent.empty()) continue;
      if(new_ent.size()!=1) {
        SETERRQ1(PETSC_COMM_SELF,1,"data inconsistency %u",new_ent.size());
      }
      rval = m_field.get_moab().tag_set_data(cOre.get_th_RefParentHandle(),&*new_ent.begin(),1,&*eit); CHKERRQ_MOAB(rval);
    }
  }

  Range seed_tets;
  if(tets_ptr!=NULL) {
    seed_tets.merge(*tets_ptr);
  }
  seed_tets = subtract(seed_tets,mother_tets);
  seed_tets = subtract(seed_tets,edge_tets);
  seed_tets.merge(created_tets);

  ierr = m_field.seed_ref_level(seed_tets,bit); CHKERRQ(ierr);

  PetscFunctionReturn(0);
}
PetscErrorCode NodeMergerInterface::mergeNodes(EntityHandle father,EntityHandle mother,BitRefLevel bit,BitRefLevel tets_from_bit_ref_level) {
  PetscFunctionBegin;
  MoFEM::Interface& m_field = cOre;
  PetscErrorCode ierr;

  Range level_tets;
  ierr = m_field.get_entities_by_type_and_ref_level(tets_from_bit_ref_level,BitRefLevel().set(),MBTET,level_tets); CHKERRQ(ierr);
  ierr = mergeNodes(father,mother,bit,&level_tets); CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

}
