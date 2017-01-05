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

#include <NodeMerger.hpp>

#include <FTensor.hpp>
#include <fem_tools.h>

namespace MoFEM {

template<class T>
static inline double determinant(T &t) {
  return
  +t(0,0)*t(1,1)*t(2,2) + t(1,0)*t(2,1)*t(0,2)
  +t(2,0)*t(0,1)*t(1,2) - t(0,0)*t(2,1)*t(1,2)
  -t(2,0)*t(1,1)*t(0,2) - t(1,0)*t(0,1)*t(2,2);
}

double volume_length_quality(double *coords) {
  PetscFunctionBegin;
  double lrms = 0;
  for(int dd = 0;dd!=3;dd++) {
    lrms +=
    pow(coords[0*3+dd]-coords[1*3+dd],2)+
    pow(coords[0*3+dd]-coords[2*3+dd],2)+
    pow(coords[0*3+dd]-coords[3*3+dd],2)+
    pow(coords[1*3+dd]-coords[2*3+dd],2)+
    pow(coords[1*3+dd]-coords[3*3+dd],2)+
    pow(coords[2*3+dd]-coords[3*3+dd],2);
  }
  lrms = sqrt((1./6.)*lrms);
  double diff_n[12];
  ShapeDiffMBTET(diff_n);
  FTensor::Tensor1<double*,3> t_diff_n(&diff_n[0],&diff_n[1],&diff_n[2],3);
  FTensor::Tensor1<double*,3> t_coords(&coords[0],&coords[1],&coords[2],3);
  FTensor::Tensor2<double,3,3> jac;
  FTensor::Index<'i',3> i;
  FTensor::Index<'j',3> j;
  jac(i,j) = 0;
  for(int nn = 0;nn!=4;nn++) {
    jac(i,j) += t_coords(i)*t_diff_n(j);
    ++t_coords;
    ++t_diff_n;
  }
  double volume = determinant(jac)*G_TET_W1[0]/6.;
  return 6.*sqrt(2.)*volume/pow(lrms,3);
}

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

PetscErrorCode NodeMergerInterface::mergeNodes(
  EntityHandle father,
  EntityHandle mother,
  BitRefLevel bit,
  Range *tets_ptr,
  const bool only_if_improve_quality,
  const double move
) {
  PetscFunctionBegin;

  MoFEM::Interface& m_field = cOre;
  PetscErrorCode ierr;
  ErrorCode rval;

  // Get adges adjacent to father and mother, i.e. mother is merged to father.
  Range father_edges;
  rval = m_field.get_moab().get_adjacencies(&father,1,1,false,father_edges); CHKERRQ_MOAB(rval);
  Range mother_edges;
  rval = m_field.get_moab().get_adjacencies(&mother,1,1,false,mother_edges); CHKERRQ_MOAB(rval);

  // Find common edge
  Range common_edge;
  common_edge = intersect(father_edges,mother_edges);
  if(tets_ptr!=NULL) {
    Range tets_edges;
    rval = m_field.get_moab().get_adjacencies(
      *tets_ptr,1,false,tets_edges,moab::Interface::UNION
    ); CHKERRQ_MOAB(rval);
    common_edge = intersect(common_edge,tets_edges);
    father_edges = intersect(father_edges,tets_edges);
    mother_edges = intersect(mother_edges,tets_edges);
  }

  // No common edge, merge no possible
  if(errorIfNoCommonEdge && common_edge.empty()) {
    SETERRQ(PETSC_COMM_SELF,1,"no common edge between nodes");
  } else if(common_edge.empty()) {
    Range seed_tets;
    if(tets_ptr!=NULL) {
      seed_tets.merge(*tets_ptr);
    }
    ierr = m_field.seed_ref_level(seed_tets,bit); CHKERRQ(ierr);
    successMerge = false;
    PetscFunctionReturn(0);
  }

  if(move) {
    EntityHandle conn[] = {father,mother};
    double coords[6];
    rval = m_field.get_moab().get_coords(conn,2,coords); CHKERRQ_MOAB(rval);
    coords[0] += move*(coords[3]-coords[0]);
    coords[1] += move*(coords[4]-coords[1]);
    coords[2] += move*(coords[5]-coords[2]);
    rval = m_field.get_moab().set_coords(conn,1,coords); CHKERRQ_MOAB(rval);
  }

  // Get tets adjacent to mother and father
  Range father_tets;
  rval = m_field.get_moab().get_adjacencies(&father,1,3,false,father_tets); CHKERRQ_MOAB(rval);
  Range mother_tets;
  rval = m_field.get_moab().get_adjacencies(&mother,1,3,false,mother_tets); CHKERRQ_MOAB(rval);

  // Common edge tets, that tests will be squashed
  Range edge_tets;
  rval = m_field.get_moab().get_adjacencies(common_edge,3,true,edge_tets); CHKERRQ_MOAB(rval);
  // Mother tets, has only one mother vertex and no father vertex.
  mother_tets = subtract(mother_tets,edge_tets);

  // Intersect with ptr_tets (usually assciated with some bit level)
  if(tets_ptr!=NULL) {
    father_tets = intersect(father_tets,*tets_ptr);
    mother_tets = intersect(mother_tets,*tets_ptr);
    edge_tets = intersect(edge_tets,*tets_ptr);
  }

  if(only_if_improve_quality) {
    double min_quality0 = 1;
    double coords[12];
    for(Range::iterator tit = mother_tets.begin(); tit!=mother_tets.end();tit++) {
      const EntityHandle* conn;
      int num_nodes;
      rval = m_field.get_moab().get_connectivity(*tit,conn,num_nodes,true); CHKERRQ_MOAB(rval);
      rval = m_field.get_moab().get_coords(conn,num_nodes,coords); CHKERRQ_MOAB(rval);
      double quality = volume_length_quality(coords);
      min_quality0 = (min_quality0>quality) ? quality : min_quality0;
    }
    for(Range::iterator tit = edge_tets.begin();tit!=edge_tets.end();tit++) {
      const EntityHandle* conn;
      int num_nodes;
      rval = m_field.get_moab().get_connectivity(*tit,conn,num_nodes,true); CHKERRQ_MOAB(rval);
      rval = m_field.get_moab().get_coords(conn,num_nodes,coords); CHKERRQ_MOAB(rval);
      double quality = volume_length_quality(coords);
      min_quality0 = (min_quality0>quality) ? quality : min_quality0;
      // cerr << "min_quality0 " << min_quality0 << endl;
    }
    double min_quality = 1;
    for(Range::iterator tit = mother_tets.begin();tit!=mother_tets.end();tit++) {
      const EntityHandle* conn;
      int num_nodes;
      rval = m_field.get_moab().get_connectivity(*tit,conn,num_nodes,true); CHKERRQ_MOAB(rval);
      EntityHandle new_conn[4];
      // Replace mother vertices by father vertices
      int nb_mother_verts = 0;
      for(int nn = 0;nn<4;nn++) {
        if(conn[nn] == father) {
          SETERRQ(
            PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,
            "Tet has father vertex, impossible but here it is"
          );
        }
        if(conn[nn] == mother) {
          new_conn[nn] = father;
          nb_mother_verts++;
        } else {
          new_conn[nn] = conn[nn];
        }
      }
      if(nb_mother_verts!=1) {
        SETERRQ1(
          PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,
          "Tet should have only one vertex but have %d",
          nb_mother_verts
        );
      }
      rval = m_field.get_moab().get_coords(new_conn,num_nodes,coords); CHKERRQ_MOAB(rval);
      double quality = volume_length_quality(coords);
      min_quality = (min_quality>quality) ? quality : min_quality;
      // cerr << "min_quality " << min_quality << endl;
    }
    if(min_quality<min_quality0) {
      Range seed_tets;
      if(tets_ptr!=NULL) {
        seed_tets.merge(*tets_ptr);
      }
      ierr = m_field.seed_ref_level(seed_tets,bit); CHKERRQ(ierr);
      successMerge = false;
      PetscFunctionReturn(0);
    }
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
    // Replace mother vertices by father vertices
    int nb_mother_verts = 0;
    for(int nn = 0;nn<4;nn++) {
      if(conn[nn] == father) {
        SETERRQ(
          PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,
          "Tet has father vertex, impossible but here it is"
        );
      }
      if(conn[nn] == mother) {
        new_conn[nn] = father;
        nb_mother_verts++;
      } else {
        new_conn[nn] = conn[nn];
      }
    }
    if(nb_mother_verts!=1) {
      SETERRQ1(
        PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,
        "Tet should have only one vertex but have %d",
        nb_mother_verts
      );
    }
    // Create tet with new connectivity
    EntityHandle tet;
    rval = m_field.get_moab().create_element(MBTET,new_conn,4,tet); CHKERRQ_MOAB(rval);
    rval = m_field.get_moab().tag_set_data(cOre.get_th_RefParentHandle(),&tet,1,&*tit); CHKERRQ_MOAB(rval);
    created_tets.insert(tet);
  }

  // Loop over father adjacent entities to use them as parents
  Range adj_ents;
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
        nb_new_node=0;
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
      rval = m_field.get_moab().get_adjacencies(
        new_conn,num_nodes,dim,true,new_ent
      ); CHKERRQ_MOAB(rval);
      // if(new_ent.empty()) continue;
      if(new_ent.size()!=1) {
        SETERRQ1(PETSC_COMM_SELF,1,"data inconsistency %u",new_ent.size());
      }
      rval = m_field.get_moab().tag_set_data(
        cOre.get_th_RefParentHandle(),&*new_ent.begin(),1,&*eit
      ); CHKERRQ_MOAB(rval);
    }
  }

  // Seed tets to given bit level
  Range seed_tets;
  if(tets_ptr!=NULL) {
    seed_tets.merge(*tets_ptr);
  }
  seed_tets = subtract(seed_tets,mother_tets);
  seed_tets = subtract(seed_tets,edge_tets);
  seed_tets.merge(created_tets);

  ierr = m_field.seed_ref_level(seed_tets,bit); CHKERRQ(ierr);

  successMerge = true;

  PetscFunctionReturn(0);
}
PetscErrorCode NodeMergerInterface::mergeNodes(
  EntityHandle father,
  EntityHandle mother,
  BitRefLevel bit,
  BitRefLevel tets_from_bit_ref_level,
  const bool only_if_improve_quality,const double move
) {
  PetscFunctionBegin;
  MoFEM::Interface& m_field = cOre;
  PetscErrorCode ierr;

  Range level_tets;
  ierr = m_field.get_entities_by_type_and_ref_level(
    tets_from_bit_ref_level,BitRefLevel().set(),MBTET,level_tets
  ); CHKERRQ(ierr);
  ierr = mergeNodes(father,mother,bit,&level_tets,only_if_improve_quality,move); CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

}
