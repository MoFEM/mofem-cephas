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

#include <petscsys.h>
#include <petscvec.h> 
#include <petscmat.h> 
#include <petscsnes.h> 
#include <petscts.h> 

#include <moab/ParallelComm.hpp>

#include <definitions.h>
#include <h1_hdiv_hcurl_l2.h>
#include <fem_tools.h>

#include <Common.hpp>

#include <LoopMethods.hpp>
#include <FieldInterface.hpp>

#include <boost/ptr_container/ptr_vector.hpp>
#include <boost/ptr_container/ptr_map.hpp>
#include <Core.hpp>

#include <BitLevelCoupler.hpp>

//Boundary Volume Tree
#include <moab/BVHTree.hpp>

#include <fem_tools.h>

//static bool debug = true;

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

  const RefMoFEMEntity_multiIndex *refined_ptr;
  ierr = m_field.get_ref_ents(&refined_ptr); CHKERRQ(ierr);
  RefMoFEMEntity_multiIndex::index<EntType_mi_tag>::type::iterator it = refined_ptr->get<EntType_mi_tag>().lower_bound(MBVERTEX);
  RefMoFEMEntity_multiIndex::index<EntType_mi_tag>::type::iterator hi_it = refined_ptr->get<EntType_mi_tag>().upper_bound(MBVERTEX);
  for(;it!=hi_it;it++) {

    //that vertex is on parent bit level, no need to process
    if((it->get_BitRefLevel()&parent_level).any()) continue;

    //check if vertex has a parent and parent is on parent bit level
    EntityHandle parent_ent;
    parent_ent = it->get_parent_ent();
    RefMoFEMEntity_multiIndex::index<Ent_mi_tag>::type::iterator pit;
    pit = refined_ptr->get<Ent_mi_tag>().find(parent_ent);
    if((pit->get_BitRefLevel()&parent_level).any()) {
      continue;
    }

    //check if vertex is on child entities set
    EntityHandle node = it->get_ref_ent();
    if(children.find(node)==children.end()) {
      continue;
    }

    //build a boundary volume Tree
    if(!init_tree) {
      tree_ptr = new BVHTree(&m_field.get_moab());
      Range tets;
      ierr = m_field.get_entities_by_type_and_ref_level(
	parent_level,BitRefLevel().set(),MBTET,tets); CHKERRQ(ierr);
      rval = tree_ptr->build_tree(tets); CHKERR_PETSC(rval);
      if(verb > 0) {
	rval = tree_ptr->print(); CHKERR_PETSC(rval);
      }
    }

    //find a leaf
    double coords[3];
    rval = m_field.get_moab().get_coords(&node,1,coords); CHKERR_PETSC(rval);
    EntityHandle leaf_out;
    rval = tree_ptr->point_search(coords,leaf_out,iter_tol,inside_tol); CHKERR_PETSC(rval);
    if(verb>0) {
      cout << "leaf_out " << leaf_out << endl;
    }
    if(vertex_elements) {
      const RefMoFEMElement_multiIndex *refined_finite_elements_ptr;
      ierr = m_field.get_ref_finite_elements(&refined_finite_elements_ptr); CHKERRQ(ierr);
      RefMoFEMElement_multiIndex::index<Ent_mi_tag>::type::iterator eit;
      eit = refined_finite_elements_ptr->get<Ent_mi_tag>().find(node);
      if(eit!=refined_finite_elements_ptr->get<Ent_mi_tag>().end()) {
	RefMoFEMElement_change_parent modifier(m_field.get_moab(),refined_ptr,refined_ptr->project<0>(it),leaf_out);
	bool success = const_cast<RefMoFEMElement_multiIndex*>(refined_finite_elements_ptr)->modify(refined_finite_elements_ptr->project<0>(eit),modifier);
	if(!success) {
	  SETERRQ(PETSC_COMM_SELF,MOFEM_OPERATION_UNSUCCESSFUL,"unsuccessful operation");
	}
      }
    } else {
      RefMoFEMEntity_change_parent modifier(m_field.get_moab(),leaf_out);
      bool success = const_cast<RefMoFEMEntity_multiIndex*>(refined_ptr)->modify(refined_ptr->project<0>(it),modifier);
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

PetscErrorCode BitLevelCouplerInterface::buidlAdjacenciesVerticesOnFacesEdgesVolumes(
  const BitRefLevel &parent_level,Range &children,bool vertex_elements,const double inside_tol,bool throw_error,int verb) {
  PetscFunctionBegin;

  FieldInterface& m_field = cOre;

  //calculate shape functions and derivatives
  double diffN[9];
  ierr = ShapeDiffMBTET(diffN); CHKERRQ(ierr);
  double N[4];
  
  //access to ref dofs multi-index
  const RefMoFEMEntity_multiIndex *refined_ptr;
  ierr = m_field.get_ref_ents(&refined_ptr); CHKERRQ(ierr);

  RefMoFEMEntity_multiIndex::index<Composite_EntType_and_ParentEntType_mi_tag>::type::iterator it,hi_it;
  it = refined_ptr->get<Composite_EntType_and_ParentEntType_mi_tag>().lower_bound(boost::make_tuple(MBVERTEX,MBTET));
  hi_it = refined_ptr->get<Composite_EntType_and_ParentEntType_mi_tag>().upper_bound(boost::make_tuple(MBVERTEX,MBTET));
  for(;it!=hi_it;it++) {

    //that vertex is on parent bit level, no need to process
    if((it->get_BitRefLevel()&parent_level).any()) continue;

    //check if vertex has a parent and parent is on parent bit level
    EntityHandle parent_ent;
    parent_ent = it->get_parent_ent();
    RefMoFEMEntity_multiIndex::index<Ent_mi_tag>::type::iterator pit;
    pit = refined_ptr->get<Ent_mi_tag>().find(parent_ent);
    if((pit->get_BitRefLevel()&parent_level).none()) {
      continue;
    }

    //check if vertex is on child entities set
    EntityHandle node = it->get_ref_ent();
    if(children.find(node)==children.end()) {
      continue;
    }

    int num_nodes;
    const EntityHandle *conn;
    rval = m_field.get_moab().get_connectivity(parent_ent,conn,num_nodes,true); CHKERR_PETSC(rval);
    double coords[12+3];
    rval = m_field.get_moab().get_coords(conn,num_nodes,coords); CHKERR_PETSC(rval);
    rval = m_field.get_moab().get_coords(&node,1,&coords[12]); CHKERR_PETSC(rval);
    for(int nn = 1;nn<4;nn++) {
      cblas_daxpy(3,-1,coords,1,&coords[nn*3],1);
    }
    for(int dd =0;dd<3;dd++) {
      coords[dd] = 0;
    }

    double loc[3];
    ierr = ShapeMBTET_inverse(N,diffN,coords,&coords[12],loc); CHKERRQ(ierr);
    ierr = ShapeMBTET(N,&loc[0],&loc[1],&loc[2],1);; CHKERRQ(ierr);

    EntityHandle parent = 0;

    if(N[0] > inside_tol && N[1] < inside_tol && N[2] < inside_tol && N[3] < inside_tol) {
      //vertex 0
      rval = m_field.get_moab().side_element(parent_ent,0,0,parent); CHKERR_PETSC(rval);
    } else
    if(N[0] < inside_tol && N[1] > inside_tol && N[2] < inside_tol && N[3] < inside_tol) {
      //vertex 1
      rval = m_field.get_moab().side_element(parent_ent,0,1,parent); CHKERR_PETSC(rval);
    } else
    if(N[0] < inside_tol && N[1] < inside_tol && N[2] > inside_tol && N[3] < inside_tol) {
      //vertex 2
      rval = m_field.get_moab().side_element(parent_ent,0,2,parent); CHKERR_PETSC(rval);
    } else
    if(N[0] < inside_tol && N[1] < inside_tol && N[2] < inside_tol && N[3] > inside_tol) {
      //vertex 3
      rval = m_field.get_moab().side_element(parent_ent,0,3,parent); CHKERR_PETSC(rval);
    } else
    if(N[0] > inside_tol && N[1] > inside_tol && N[2] < inside_tol && N[3] < inside_tol) {
      //edge 0
      rval = m_field.get_moab().side_element(parent_ent,1,0,parent); CHKERR_PETSC(rval);
    } else
    if(N[0] < inside_tol && N[1] > inside_tol && N[2] > inside_tol && N[3] < inside_tol) {
      //edge 1
      rval = m_field.get_moab().side_element(parent_ent,1,1,parent); CHKERR_PETSC(rval);
    } else
    if(N[0] > inside_tol && N[1] < inside_tol && N[2] > inside_tol && N[3] < inside_tol) {
      //edge 2
      rval = m_field.get_moab().side_element(parent_ent,1,2,parent); CHKERR_PETSC(rval);
    } else
    if(N[0] > inside_tol && N[1] < inside_tol && N[2] < inside_tol && N[3] > inside_tol) {
      //edge 3
      rval = m_field.get_moab().side_element(parent_ent,1,3,parent); CHKERR_PETSC(rval);
    } else
    if(N[0] < inside_tol && N[1] > inside_tol && N[2] < inside_tol && N[3] < inside_tol) {
      //edge 4
      rval = m_field.get_moab().side_element(parent_ent,1,4,parent); CHKERR_PETSC(rval);
    } else
    if(N[0] > inside_tol && N[1] < inside_tol && N[2] > inside_tol && N[3] < inside_tol) {
      //edge 5
      rval = m_field.get_moab().side_element(parent_ent,1,5,parent); CHKERR_PETSC(rval);
    } else
    if(N[0] > inside_tol && N[1] > inside_tol && N[2] < inside_tol && N[3] > inside_tol) {
      //face 0
      rval = m_field.get_moab().side_element(parent_ent,2,0,parent); CHKERR_PETSC(rval);
    } else
    if(N[0] < inside_tol && N[1] > inside_tol && N[2] > inside_tol && N[3] < inside_tol) {
      //face 1
      rval = m_field.get_moab().side_element(parent_ent,2,1,parent); CHKERR_PETSC(rval);
    } else
    if(N[0] > inside_tol && N[1] < inside_tol && N[2] > inside_tol && N[3] > inside_tol) {
      //face 2
      rval = m_field.get_moab().side_element(parent_ent,2,2,parent); CHKERR_PETSC(rval);
    } else
    if(N[0] > inside_tol && N[1] > inside_tol && N[2] > inside_tol && N[3] < inside_tol) {
      //face 3
      rval = m_field.get_moab().side_element(parent_ent,2,3,parent); CHKERR_PETSC(rval);
    } else
    if(N[0] > inside_tol && N[1] > inside_tol && N[2] > inside_tol && N[3] < inside_tol) {
      //volume
      continue;
    } else {
      if(throw_error) {
	SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCT,
	  "vertex is not on any entity of volume incling volume itself");
      }
    }

    if(!parent) continue;

    if(vertex_elements) {
      const RefMoFEMElement_multiIndex *refined_finite_elements_ptr;
      ierr = m_field.get_ref_finite_elements(&refined_finite_elements_ptr); CHKERRQ(ierr);
      RefMoFEMElement_multiIndex::index<Ent_mi_tag>::type::iterator eit;
      eit = refined_finite_elements_ptr->get<Ent_mi_tag>().find(node);
      if(eit!=refined_finite_elements_ptr->get<Ent_mi_tag>().end()) {
	RefMoFEMElement_change_parent modifier(m_field.get_moab(),refined_ptr,refined_ptr->project<0>(it),parent);
	bool success = const_cast<RefMoFEMElement_multiIndex*>(refined_finite_elements_ptr)->modify(refined_finite_elements_ptr->project<0>(eit),modifier);
	if(!success) {
	  SETERRQ(PETSC_COMM_SELF,MOFEM_OPERATION_UNSUCCESSFUL,"unsuccessful operation");
	}
      }
    } else {
      RefMoFEMEntity_change_parent modifier(m_field.get_moab(),parent);
      bool success = const_cast<RefMoFEMEntity_multiIndex*>(refined_ptr)->modify(refined_ptr->project<0>(it),modifier);
      if(!success) {
	SETERRQ(PETSC_COMM_SELF,MOFEM_OPERATION_UNSUCCESSFUL,"unsuccessful operation");
      }
    }

  }
  PetscFunctionReturn(0);
}

PetscErrorCode BitLevelCouplerInterface::buidlAdjacenciesEdgesFacesVolumes(
  const BitRefLevel &parent_level,Range &children,bool elements,int verb) {
  PetscFunctionBegin;

  FieldInterface& m_field = cOre;

  //access to ref dofs multi-index
  const RefMoFEMEntity_multiIndex *refined_ptr;
  ierr = m_field.get_ref_ents(&refined_ptr); CHKERRQ(ierr);

  vector<EntityHandle> conn_parents;

  Range::iterator eit,hi_eit;
  eit = children.begin();
  hi_eit = children.end();
  for(;eit!=hi_eit;eit++) {

    //check entity type
    EntityType type;
    type = m_field.get_moab().type_from_handle(*eit);
    switch (type) {
      case MBEDGE:
      case MBTRI:
      case MBTET:
      break;
      default:
      SETERRQ(PETSC_COMM_SELF,MOFEM_NOT_IMPLEMENTED,"not implemented");
    }

    //get ref entity iterator
    RefMoFEMEntity_multiIndex::index<Ent_mi_tag>::type::iterator it;
    it = refined_ptr->get<Ent_mi_tag>().find(*eit);
    //that entity is on parent bit level, no need to process
    if((it->get_BitRefLevel()&parent_level).any()) continue;

    //check if vertex has a parent and parent is on parent bit level
    EntityHandle parent_ent;
    parent_ent = it->get_parent_ent();
    RefMoFEMEntity_multiIndex::index<Ent_mi_tag>::type::iterator pit;
    pit = refined_ptr->get<Ent_mi_tag>().find(parent_ent);
    if((pit->get_BitRefLevel()&parent_level).none()) {
      continue;
    }

    //connectivity
    int max_dim = m_field.get_moab().dimension_from_handle(*eit);
    const EntityHandle *conn;
    int num_nodes;
    ierr = m_field.get_moab().get_connectivity(*eit,conn,num_nodes); CHKERRQ(ierr);
    conn_parents.resize(num_nodes);
    for(int nn = 0;nn<num_nodes;nn++) {
      RefMoFEMEntity_multiIndex::index<Ent_mi_tag>::type::iterator cit;
      cit = refined_ptr->get<Ent_mi_tag>().find(conn[nn]);
      if((cit->get_BitRefLevel()&parent_level).none()) {
	SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCT,"parent of vertex is not on parent bit level");
      }
      conn_parents[nn] = cit->get_parent_ent();
      int ent_dim = m_field.get_moab().dimension_from_handle(conn_parents[nn]);
      max_dim = ent_dim > max_dim ? ent_dim : max_dim;
    }

    if(max_dim > 1) {
      for(;max_dim<=3;max_dim++) {
      
	Range parent_ents;
	rval = m_field.get_moab().get_adjacencies(&*conn_parents.begin(),num_nodes,max_dim,false,parent_ents); CHKERR_PETSC(rval);
	if(!parent_ents.empty()) {
	  if(elements) {

	  } else {

	  }
	  break;
	}

      }
    }

  }

  PetscFunctionReturn(0);
}

}
