/** \file BitLevelCoupler.cpp
 * \brief BitLevelCoupler interface implementation

 *
 * MoFEM is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
 * License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with MoFEM. If not, see <http://www.gnu.org/licenses/>
*/

#include <Common.hpp>
#include <LoopMethods.hpp>
#include <Core.hpp>
#include <FieldInterface.hpp>

//Boundary Volume Tree
#include <BVHTree.hpp>

#include <fem_tools.h>

PetscErrorCode ierr;
ErrorCode rval;

namespace MoFEM {

PetscErrorCode BitLevelCouplerInterface::queryInterface(const MOFEMuuid& uuid, FieldUnknownInterface** iface) {
  PetscFunctionBegin;
  *iface = NULL;
  if(uuid == IDD_MOFENBitLevelCoupler) {
    *iface = dynamic_cast<BitLevelCouplerInterface*>(this);
    PetscFunctionReturn(0);
  }
  if(uuid == IDD_MOFEMUnknown) {
    *iface = dynamic_cast<FieldUnknownInterface*>(this);
    PetscFunctionReturn(0);
  }
  SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCT,"unknown interface");
  PetscFunctionReturn(0);
}

PetscErrorCode BitLevelCouplerInterface::buidlAdjacenciesVerticesOnTets(const BitRefLevel &parent_level,Range &children,
    const double iter_tol,const double inside_tol,bool vertex_elements,int verb) {

  PetscFunctionBegin;
  FieldInterface& m_field = cOre;
  //build BVHTree
  bool init_tree = false;
  BVHTree *tree_ptr;
  
  //find parents of all nodes, if node has no parent then tetrahedral containing that node is searched
  //node on tetrahedra my by part of face or edge on that tetrahedral, this need to be verified

  FieldInterface& m_field = cOre;
  const RefMoFEMEntity_multiIndex *refined_ptr;
  ierr = m_field.get_ref_ents(&refined_ptr); CHKERRQ(ierr);
  RefMoFEMEntity_multiIndex::index<EntType_mi_tag>::type::iterator it = refined_ptr->get()<EntType_mi_tag>.lower_bound(MBVERTEX);
  RefMoFEMEntity_multiIndex::index<EntType_mi_tag>::type::iterator hi_it = refined_ptr->get()<EntType_mi_tag>.upper_bound(MBVERTEX);
  for(;it!=get_refined_ptr.end();it++) {

    //that vertex is on parent bit level, no need to process
    if((it->get_BitRefLevel()&parent_level).any()) continue;

    //check if vertex has a parent and parent is on parent bit level
    EntityHandle parent_ent;
    parent_ent = it->get_parent_ent();
    RefMoFEMEntity_multiIndex *pit = refined_ptr.get()<Ent_mi_tag>().find(parent_ent);
    if((pit->get_BitRefLevel()&parent_level).any()) {
      continue;
    }

    //check if vertex is on child entities set
    EntityHandle node = it->get_ent();
    if(children.find(node)==refined_ptr.end()) {
      continue;
    }

    //build a boundary volume Tree
    if(!init_tree) {
      tree_ptr = new BVHTree(m_field.get_moab());
      Range tets;
      rval = m_field.get_moab().get_entities_by_type_and_ref_level(
	parent_level,BitRefLevel().set(),MBTET,tets); CHKERR(rval);
      rval = tree_ptr->build_tree(tets); CHKERR(rval);
      if(verb > 0) {
	rval = tree_ptr->print(); CHKERR(rval);
      }
    }

    //find a leaf
    double coords[3];
    rval = m_field.get_moab().get_coords(&node,1,coords); CHKERR(rval);
    EntityHandle leaf_out;
    rval = tree_ptr->point_search(coords,leaf_out,iter_tol,inside_tol); CHKERR(rval);
    if(verb>0) {
      cout << "leaf_out " << leaf_out << endl;
    }
    if(vertex_elements) {
      RefMoFEMElement_multiIndex *refined_finite_elements_ptr;
      ierr = m_field.get_ref_finite_elements(&refined_finite_elements_ptr); CHKERRQ(ierr);
      RefMoFEMElement_multiIndex::index<Ent_mi_tag>::type::ietrator eit;
      eit = refined_finite_elements_ptr.get<Ent_mi_tag>().find(node);
      if(eit!=refined_finite_elements_ptr.get<Ent_mi_tag>().end()) {
	bool success = refined_finite_elements_ptr->modify(refined_finite_elements_ptr->project<0>(it),
	  RefMoFEMElement_change_parent(m_field.get_moab(),refined_ptr,refined_ptr->project<0>(eit),leaf_out);
	if(!success) {
	  SETERRQ(PETSC_COMM_SELF,MOFEM_OPERATION_UNSUCCESSFUL,"unsuccessful operation");
	}
      }
    } else {
      bool success = RefMoFEMElement_change_parent(m_field.get_moab(),refined_ptr,refined_ptr->project<0>(eit),leaf_out);
      if(!success) {
	SETERRQ(PETSC_COMM_SELF,MOFEM_OPERATION_UNSUCCESSFUL,"unsuccessful operation");
      }
    }
  }
  if(init_tree) {
    delete tree_ptr;
  }
  PetscFunctionReturn(0);
}

PetscErrorCode buidlAdjacenciesVerticesOnFacesEdgesVolumes(
  const BitRefLevel &parent_level,Range &children,bool vertex_elements,const double inside_tol,int verb = 0) {
  PetscFunctionBegin;

  //calculate shape functions and derivatives
  double diffN[9];
  ierr = ShapeDiffMBTET(diffN); CHKERRQ(ierr);
  double N[4];
  double loc[] = { 0,0,0 };
  ierr = ShapeMBTET(N,&loc[0],&loc[1],&loc[2],1);; CHKERRQ(ierr);

  //access to ref dofs multi-index
  const RefMoFEMEntity_multiIndex *refined_ptr;
  ierr = m_field.get_ref_ents(&refined_ptr); CHKERRQ(ierr);

  FieldInterface& m_field = cOre;
  RefMoFEMEntity_multiIndex::index<Composite_EntType_and_ParentEntType_mi_tag>::type::iterator it,hi_it;
  it = refined_ptr->get<Composite_EntType_and_ParentEntType_mi_tag>().lower_bound(boost::make_tuple(MBVERTEX,MBTET));
  hi_it = refined_ptr->get<Composite_EntType_and_ParentEntType_mi_tag>().upper_bound(boost::make_tuple(MBVERTEX,MBTET));
  for(;it!=hi_it;it++) {

    //that vertex is on parent bit level, no need to process
    if((it->get_BitRefLevel()&parent_level).any()) continue;

    //check if vertex has a parent and parent is on parent bit level
    EntityHandle parent_ent;
    parent_ent = it->get_parent_ent();
    RefMoFEMEntity_multiIndex *pit = refined_ptr.get()<Ent_mi_tag>().find(parent_ent);
    if((pit->get_BitRefLevel()&parent_level).none()) {
      continue;
    }

    //check if vertex is on child entities set
    EntityHandle node = it->get_ent();
    if(children.find(node)==refined_ptr.end()) {
      continue;
    }

    double node_coords[3];
    rval = m_field.get_moab().get_coords(&node,1,node_coords); CHKERR(rval);


    Range parent_faces;
    rval = moab.get_adjacencies(&parent_ent,1,2,false,parent_faces); CHKERR_PETSC(rval);
    Range::iterator fit,hi_fit;
    fit = parent_faces.begin();
    hi_fit = parent_faces.end();
    while(fit!=hi_fit) {

      EntityType type;
      rval = m_field.get_moab().get_entity_type(*fit,type); CHKERR(rval);
      if(type != MBTRI) SETERRQ(PETSC_COMM_SELF,MOFEM_NOT_IMPLEMENTED,"not implemented");

      int count;
      double *x,*y,*z;
      rval = m_field.get_moab().coords_iterate(fit,hi_fit,x,y,z,count); CHKERR(rval);

      double coords[9] = { 
	0, 0, 0, x[1]-x[0], y[1]-z[0], z[1]-z[0], x[2]-x[0], y[2]-y[0], z[2]-z[0],
	node_coords[0]-x[0], node_coords[1]-z[0], node_coords[2]-y[0] };
      ierr = ShapeMBTET_inverse(N,diffN,elem_coords,glob_coords,loc); CHKERRQ(ierr);

      if(fabs(loc[2]) < inside_tol) continue;

      EntityHandle parent = 0;

      //check if is on face
      if(
	fabs(loc[0]) < inside_tol || 
	fabs(loc[1]) < inside_tol ||
	fabs(loc[0]-1) < inside_tol ||
	fabs(loc[1]-1) < inside_tol {

	const EntityHandle *conn;
	int num_nodes;
	ierr = m_field.get_moab().get_connectivity(&*fit,conn,num_nodes); CHKERRQ(ierr);

	//check vertices
	if(fabs(loc[0]) < inside_tol && fabs(loc[1]) < inside_tol) parent = conn[0]
	else
	if(fabs(loc[1]-1) < inside_tol && fabs(loc[0]) < inside_tol) parent = conn[1]
	else
	if(fabs(loc[0]-1) < inside_tol && fabs(loc[1]) < inside_tol) perent = conn[2];

	if(parent !=0) break;

	//check edges
	if(fabs(loc[1]) < inside_tol) {
	  ierr = m_field.get_moab().side_element(*fit,1,0,parent); CHKERRQ(ierr);
	} else 
	if(fabs(loc[0]) < inside_tol) {
	  ierr = m_field.get_moab().side_element(*fit,1,2,parent); CHKERRQ(ierr);
	} else 
	if(fabs(loc[0]+loc[1]-1) < inside_tol) {
	  ierr = m_field.get_moab().side_element(*fit,1,1,parent); CHKERRQ(ierr);
	}
	
	if(parent !=0) break;
	
	if(loc[0] > inside_tol &&
	  loc[1] > inside_tol &&
	  loc[0]+loc[1] < 1-inside_tol) {
	  parent = *fit;
	}

      }

      if(parent != 0) {
	if(vertex_elements) {
	  RefMoFEMElement_multiIndex *refined_finite_elements_ptr;
	  ierr = m_field.get_ref_finite_elements(&refined_finite_elements_ptr); CHKERRQ(ierr);
	  RefMoFEMElement_multiIndex::index<Ent_mi_tag>::type::ietrator eit;
	  eit = refined_finite_elements_ptr.get<Ent_mi_tag>().find(node);
	  if(eit!=refined_finite_elements_ptr.get<Ent_mi_tag>().end()) {
	    bool success = refined_finite_elements_ptr->modify(refined_finite_elements_ptr->project<0>(it),
	    RefMoFEMElement_change_parent(m_field.get_moab(),refined_ptr,refined_ptr->project<0>(eit),*fit);
	    if(!success) {
	      SETERRQ(PETSC_COMM_SELF,MOFEM_OPERATION_UNSUCCESSFUL,"unsuccessful operation");
	    }
	  }
	} else {
	  bool success = RefMoFEMElement_change_parent(m_field.get_moab(),refined_ptr,refined_ptr->project<0>(eit),*fit);
	  if(!success) {
	    SETERRQ(PETSC_COMM_SELF,MOFEM_OPERATION_UNSUCCESSFUL,"unsuccessful operation");
	  }
	}
	break;
      }
      fit += count;
    }

  }
  PetscFunctionReturn(0);
}

PetscErroCode buidlAdjacenciesEdgesFacesVolumes(
  const BitRefLevel &parent_level,Range &children,bool elements,int verb) {
  PetscFunctionBegin;

  //access to ref dofs multi-index
  const RefMoFEMEntity_multiIndex *refined_ptr;
  ierr = m_field.get_ref_ents(&refined_ptr); CHKERRQ(ierr);

  vector<EntityHandle> conn_parents;

  Range::iterator eit,hi_eit;
  eit = children.begin();
  hi_eit = children.end();
  for(;it!=hi_it;it++) {

    //check entity type
    EntityType type;
    rval = m_field.get_moab().get_entity_type(*eit,type); CHKERR(rval);
    switch (type) {
      case MBEDGE:
      case MBTRI:
      case MBTET:
      break;
      default:
      SETERRQ(PETSC_COMM_SELF,MOFEM_NOT_IMPLEMENTED,"not implemented");
    }

    //get ref entity iterator
    RefMoFEMEntity_multiIndex *it;
    it = refined_ptr.get()<Ent_mi_tag>().find(*eit);
    //that entity is on parent bit level, no need to process
    if((it->get_BitRefLevel()&parent_level).any()) continue;

    //check if vertex has a parent and parent is on parent bit level
    EntityHandle parent_ent;
    parent_ent = it->get_parent_ent();
    RefMoFEMEntity_multiIndex *pit;
    pit  = refined_ptr.get()<Ent_mi_tag>().find(parent_ent);
    if((pit->get_BitRefLevel()&parent_level).none()) {
      continue;
    }

    //connectivity
    int max_dim = m_field.get_moab().dimension_from_handle(*fit)
    const EntityHandle *conn;
    int num_nodes;
    ierr = m_field.get_moab().get_connectivity(&*fit,conn,num_nodes); CHKERRQ(ierr);
    conn_parents.resize(num_nodes);
    for(int nn = 0;nn<num_nodes;nn++) {
      RefMoFEMEntity_multiIndex *cit;
      cit = refined_ptr.get()<Ent_mi_tag>().find(conn[nn]);
      if(cit->get_BitRefLevel()&parent_level.none()) {
	SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCT,"parent of vertex is not on parent bit level");
      }
      conn_parents[nn] = cit->get_parent_ent();
      int ent_dim = m_field.get_moab().dimension_from_handle(conn_parents[nn]);
      max_dim = ent_dim > max_dim ? ent_dim : max_dim;
    }

    if(max_dim > 1) {
      for(;max_dim<=3;max_dim++) {
      
	Range parent_ents;
	rval = moab.get_adjacencies(&*conn_parent_ents.begin(),num_nodes,max_dim,false,parent_ents); CHKERR_PETSC(rval);
	if(!parent_ents.empty()) {
	  if(elements) {
	    RefMoFEMElement_multiIndex *refined_finite_elements_ptr;
	    ierr = m_field.get_ref_finite_elements(&refined_finite_elements_ptr); CHKERRQ(ierr);
	    RefMoFEMElement_multiIndex::index<Ent_mi_tag>::type::ietrator eit;
	    eit = refined_finite_elements_ptr.get<Ent_mi_tag>().find(*fit);
	    if(eit!=refined_finite_elements_ptr.get<Ent_mi_tag>().end()) {
	      bool success = refined_finite_elements_ptr->modify(refined_finite_elements_ptr->project<0>(eit),
	      RefMoFEMElement_change_parent(m_field.get_moab(),refined_ptr,refined_ptr->project<0>(eit),&*parent_ents.begin());
	      if(!success) {
		SETERRQ(PETSC_COMM_SELF,MOFEM_OPERATION_UNSUCCESSFUL,"unsuccessful operation");
	      }
	   }
	  } else {
	    bool success = RefMoFEMElement_change_parent(m_field.get_moab(),refined_ptr,refined_ptr->project<0>(eit),&*parent_ents.begin());
	    if(!success) {
	      SETERRQ(PETSC_COMM_SELF,MOFEM_OPERATION_UNSUCCESSFUL,"unsuccessful operation");
	    }
	  }
	  break;
	}

      }
    }

  }

  PetscFunctionReturn(0);
}

}
