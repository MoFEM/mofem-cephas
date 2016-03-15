/** \file PrismsFromSurfaceInterface.cpp
 * \brief Interface for creating prisms from surface elements
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
#include <FieldInterface.hpp>
#include <MeshRefinment.hpp>
#include <PrismInterface.hpp>
#include <SeriesRecorder.hpp>
#include <Core.hpp>

#include <PrismsFromSurfaceInterface.hpp>

namespace MoFEM {

PetscErrorCode PrismsFromSurfaceInterface::queryInterface(const MOFEMuuid& uuid, UnknownInterface** iface) {
  PetscFunctionBegin;
  *iface = NULL;
  if(uuid == IDD_MOFEMPrismsFromSurface) {
    *iface = dynamic_cast<PrismsFromSurfaceInterface*>(this);
    PetscFunctionReturn(0);
  }
  if(uuid == IDD_MOFEMUnknown) {
    *iface = dynamic_cast<UnknownInterface*>(this);
    PetscFunctionReturn(0);
  }
  SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"unknown interface");
  PetscFunctionReturn(0);
}

PetscErrorCode PrismsFromSurfaceInterface::createPrisms(const Range &ents,Range &prisms,int verb) {
  PetscFunctionBegin;
  MoABErrorCode rval;
  FieldInterface& m_field = cOre;
  Range tris = ents.subset_by_type(MBTRI);
  for(Range::iterator tit = tris.begin();tit!=tris.end();tit++) {
    const EntityHandle* conn;
    int number_nodes = 0;
    rval = m_field.get_moab().get_connectivity(*tit,conn,number_nodes,false); CHKERRQ_MOAB(rval);
    double coords[3*number_nodes];
    rval = m_field.get_moab().get_coords(conn,number_nodes,coords); CHKERRQ_MOAB(rval);
    EntityHandle prism_nodes[6];
    for(int nn = 0;nn<3;nn++) {
      prism_nodes[nn] = conn[nn];
      if(createdVertices.find(conn[nn])!=createdVertices.end()) {
        prism_nodes[3+nn] = createdVertices[prism_nodes[nn]];
      } else {
        rval = m_field.get_moab().create_vertex(&coords[3*nn],prism_nodes[3+nn]); CHKERRQ_MOAB(rval);
        createdVertices[conn[nn]] = prism_nodes[3+nn];
        rval = m_field.get_moab().tag_set_data(
          cOre.get_th_RefParentHandle(),&prism_nodes[3+nn],1,&prism_nodes[nn]
        ); CHKERRQ_MOAB(rval);
      }
    }
    EntityHandle prism;
    rval = m_field.get_moab().create_element(MBPRISM,prism_nodes,6,prism); CHKERRQ(rval);
    Range edges;
    rval = m_field.get_moab().get_adjacencies(&prism,1,1,true,edges,Interface::UNION); CHKERRQ_MOAB(rval);
    Range faces;
    rval = m_field.get_moab().get_adjacencies(&prism,1,2,true,faces,Interface::UNION); CHKERRQ_MOAB(rval);
    prisms.insert(prism);
    for(int ee = 0;ee<=2;ee++) {
      EntityHandle e1;
      rval = m_field.get_moab().side_element(prism,1,ee,e1);
      EntityHandle e2;
      rval = m_field.get_moab().side_element(prism,1,ee+6,e2); CHKERRQ_MOAB(rval);
      rval = m_field.get_moab().tag_set_data(cOre.get_th_RefParentHandle(),&e2,1,&e1); CHKERRQ_MOAB(rval);
    }
    EntityHandle f3,f4;
    {
      rval = m_field.get_moab().side_element(prism,2,3,f3);
      rval = m_field.get_moab().side_element(prism,2,4,f4);
      rval = m_field.get_moab().tag_set_data(cOre.get_th_RefParentHandle(),&f4,1,&f3); CHKERRQ_MOAB(rval);
    }
    if(number_nodes>3) {
      EntityHandle meshset;
      rval = m_field.get_moab().create_meshset(MESHSET_SET,meshset); CHKERRQ_MOAB(rval);
      rval = m_field.get_moab().add_entities(meshset,&f4,1); CHKERRQ_MOAB(rval);
      for(int ee = 0;ee<=2;ee++) {
        EntityHandle e2;
        rval = m_field.get_moab().side_element(prism,1,ee+6,e2); CHKERRQ_MOAB(rval);
        rval = m_field.get_moab().add_entities(meshset,&e2,1); CHKERRQ_MOAB(rval);
      }
      rval = m_field.get_moab().convert_entities(meshset,true,false,false); CHKERRQ_MOAB(rval);
      rval = m_field.get_moab().delete_entities(&meshset,1); CHKERRQ_MOAB(rval);
      const EntityHandle* conn_f4;
      int number_nodes_f4 = 0;
      rval = m_field.get_moab().get_connectivity(f4,conn_f4,number_nodes_f4,false); CHKERRQ_MOAB(rval);
      if(number_nodes_f4 != number_nodes) {
        SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"data inconsistency");
      }
      rval = m_field.get_moab().set_coords(&conn_f4[3],3,&coords[9]); CHKERRQ_MOAB(rval);
    }
  }
  PetscFunctionReturn(0);
}

PetscErrorCode PrismsFromSurfaceInterface::seedPrismsEntities(Range &prisms,const BitRefLevel &bit,int verb) {
  PetscFunctionBegin;
  PetscErrorCode ierr;
  MoABErrorCode rval;
  FieldInterface& m_field = cOre;
  const RefMoFEMEntity_multiIndex *const_refined_entities_ptr;
  ierr = m_field.get_ref_ents(&const_refined_entities_ptr); CHKERRQ(ierr);
  MPI_Comm comm = m_field.get_comm();
  RefMoFEMEntity_multiIndex *refined_entities_ptr;
  refined_entities_ptr = const_cast<RefMoFEMEntity_multiIndex *>(const_refined_entities_ptr);
  if(!prisms.empty()) {
    int dim = m_field.get_moab().dimension_from_handle(prisms[0]);
    for(int dd = 0;dd<=dim;dd++) {
      Range ents;
      rval = m_field.get_moab().get_adjacencies(prisms,dd,true,ents,Interface::UNION); CHKERRQ_MOAB(rval);
      Range::iterator eit = ents.begin();
      for(;eit!=ents.end();eit++) {
        pair<RefMoFEMEntity_multiIndex::iterator,bool> p_ent = refined_entities_ptr->insert(
          RefMoFEMEntity(m_field.get_moab(),*eit)
        );
        bool success = refined_entities_ptr->modify(p_ent.first,RefMoFEMEntity_change_add_bit(bit));
        if(!success) {
          SETERRQ(PETSC_COMM_SELF,MOFEM_OPERATION_UNSUCCESSFUL,"modification unsuccessful");
        }
        if(verb>2) {
          ostringstream ss;
          ss << *(p_ent.first);
          PetscSynchronizedPrintf(comm,"%s\n",ss.str().c_str());
        }
      }
    }
  }
  PetscFunctionReturn(0);
}


}
