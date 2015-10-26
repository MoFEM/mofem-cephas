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

//Tree
//#include <moab/BVHTree.hpp>
#include <moab/AdaptiveKDTree.hpp>
#include <boost/scoped_ptr.hpp>
#include <BitLevelCoupler.hpp>

#include <fem_tools.h>

//static bool debug = true;

static PetscErrorCode ierr;
static ErrorCode rval;

namespace MoFEM {

PetscErrorCode BitLevelCouplerInterface::queryInterface(const MOFEMuuid& uuid, FieldUnknownInterface** iface) {
  PetscFunctionBegin;
  *iface = NULL;
  if(uuid == IDD_MOFEMBitLevelCoupler) {
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

PetscErrorCode BitLevelCouplerInterface::buildTree(const BitRefLevel &parent_level,int verb) {
  PetscFunctionBegin;
  FieldInterface& m_field = cOre;
  treePtr.reset(new AdaptiveKDTree(&m_field.get_moab()));
  Range tets;
  ierr = m_field.get_entities_by_type_and_ref_level(
    parent_level,BitRefLevel().set(),MBTET,tets); CHKERRQ(ierr);
  rval = treePtr->build_tree(tets); CHKERR_PETSC(rval);
  if(verb > 2) {
    rval = treePtr->print(); CHKERR_PETSC(rval);
  }
  PetscFunctionReturn(0);
}

PetscErrorCode BitLevelCouplerInterface::resetTree(const BitRefLevel &parent_level,int verb) {
  PetscFunctionBegin;
  treePtr->reset_tree();
  treePtr.reset();
  PetscFunctionReturn(0);
}

PetscErrorCode BitLevelCouplerInterface::getParent(const double *coords,EntityHandle &parent,
  bool tet_only,const double iter_tol,const double inside_tol,int verb) {
  PetscFunctionBegin;
  FieldInterface& m_field = cOre;
  EntityHandle leaf_out;
  rval = treePtr->point_search(coords,leaf_out,iter_tol,inside_tol); CHKERR_PETSC(rval);
  bool is_in;
  Range tets;
  ierr = m_field.get_moab().get_entities_by_type(leaf_out,MBTET,tets); CHKERRQ(ierr);
  Range::iterator tit = tets.begin();
  for(;tit!=tets.end();tit++) {
    ierr = getLocCoordsOnTet(*tit,coords,verb); CHKERRQ(ierr);
    is_in = true;
    for(int nn = 0;nn<4;nn++) {
      if(N[nn] < -inside_tol || N[nn] > 1+inside_tol)   {
	is_in = false;
	break;
      }
      if(!is_in) break;
    }
    parent = 0;
    if(is_in) {
      if(!tet_only) {
	//vertices
	if(fabs(N[0]-1) < inside_tol && fabs(N[1])<inside_tol && fabs(N[2])<inside_tol && fabs(N[3])<inside_tol) {
	  if(verb>1) cout << "node 0 found " << endl;
	  rval = m_field.get_moab().side_element(*tit,0,0,parent); CHKERR_PETSC(rval);
	  PetscFunctionReturn(0);
	}
	if(fabs(N[0]) < inside_tol && fabs(N[1]-1)<inside_tol && fabs(N[2])<inside_tol && fabs(N[3])<inside_tol) {
	  if(verb>1) cout << "node 1 found " << endl;
	  rval = m_field.get_moab().side_element(*tit,0,1,parent); CHKERR_PETSC(rval);
	  PetscFunctionReturn(0);
	}
	if(fabs(N[0]) < inside_tol && fabs(N[1])<inside_tol && fabs(N[2]-2)<inside_tol && fabs(N[3])<inside_tol) {
	  if(verb>1) cout << "node 2 found " << endl;
	  rval = m_field.get_moab().side_element(*tit,0,2,parent); CHKERR_PETSC(rval);
	  PetscFunctionReturn(0);
	}
	if(fabs(N[0]) < inside_tol && fabs(N[1])<inside_tol && fabs(N[2])<inside_tol && fabs(N[3]-1)<inside_tol) {
	  if(verb>1) cout << "node 3 found " << endl;
	  rval = m_field.get_moab().side_element(*tit,0,3,parent); CHKERR_PETSC(rval);
	  PetscFunctionReturn(0);
	}
	//edges
	if(fabs(N[0])>inside_tol && fabs(N[1])>inside_tol && fabs(N[2])<inside_tol && fabs(N[3])<inside_tol) {
	  if(verb>1) cout << "edge 0 found " << endl;
	  rval = m_field.get_moab().side_element(*tit,1,0,parent); CHKERR_PETSC(rval);
	  PetscFunctionReturn(0);
	}
	if(fabs(N[0])<inside_tol && fabs(N[1])>inside_tol && fabs(N[2])>inside_tol && fabs(N[3])<inside_tol) {
	  if(verb>1) cout << "edge 1 found " << endl;
	  rval = m_field.get_moab().side_element(*tit,1,1,parent); CHKERR_PETSC(rval);
	  PetscFunctionReturn(0);
	}
	if(fabs(N[0])>inside_tol && fabs(N[1])<inside_tol && fabs(N[2])>inside_tol && fabs(N[3])<inside_tol) {
	  if(verb>1) cout << "edge 2 found " << endl;
	  rval = m_field.get_moab().side_element(*tit,1,2,parent); CHKERR_PETSC(rval);
	  PetscFunctionReturn(0);
	}
	if(fabs(N[0])>inside_tol && fabs(N[1])<inside_tol && fabs(N[2])<inside_tol && fabs(N[3])>inside_tol) {
	  if(verb>1) cout << "edge 3 found " << endl;
	  rval = m_field.get_moab().side_element(*tit,1,3,parent); CHKERR_PETSC(rval);
	  PetscFunctionReturn(0);
	}
	if(fabs(N[0])<inside_tol && fabs(N[1])>inside_tol && fabs(N[2])<inside_tol && fabs(N[3])>inside_tol) {
	  if(verb>1) cout << "edge 4 found " << endl;
	  rval = m_field.get_moab().side_element(*tit,1,4,parent); CHKERR_PETSC(rval);
	  PetscFunctionReturn(0);
	}
	if(fabs(N[0])<inside_tol && fabs(N[1])<inside_tol && fabs(N[2])>inside_tol && fabs(N[3])>inside_tol) {
	  if(verb>1) cout << "edge 5 found " << endl;
	  rval = m_field.get_moab().side_element(*tit,1,5,parent); CHKERR_PETSC(rval);
	  PetscFunctionReturn(0);
	}
	//faces
	if(fabs(N[0])>inside_tol && fabs(N[1])>inside_tol && fabs(N[2])<inside_tol && fabs(N[3])>inside_tol) {
	  if(verb>1) cout << "face 0 found " << endl;
	  rval = m_field.get_moab().side_element(*tit,2,0,parent); CHKERR_PETSC(rval);
	  PetscFunctionReturn(0);
	}
	if(fabs(N[0])<inside_tol && fabs(N[1])>inside_tol && fabs(N[2])>inside_tol && fabs(N[3])>inside_tol) {
	  if(verb>1) cout << "face 1 found " << endl;
	  rval = m_field.get_moab().side_element(*tit,2,1,parent); CHKERR_PETSC(rval);
	  PetscFunctionReturn(0);
	}
	if(fabs(N[0])>inside_tol && fabs(N[1])<inside_tol && fabs(N[2])>inside_tol && fabs(N[3])>inside_tol) {
	  if(verb>1) cout << "face 2 found " << endl;
	  rval = m_field.get_moab().side_element(*tit,2,2,parent); CHKERR_PETSC(rval);
	  PetscFunctionReturn(0);
	}
	if(fabs(N[0])>inside_tol && fabs(N[1])>inside_tol && fabs(N[2])>inside_tol && fabs(N[3])<inside_tol) {
	  if(verb>1) cout << "face 3 found " << endl;
	  rval = m_field.get_moab().side_element(*tit,2,2,parent); CHKERR_PETSC(rval);
	  PetscFunctionReturn(0);
	}
	//set parent
	if(parent!=0) {
	  break;
	}
      }
      if(verb>1) cout << "tet found " << endl;
      parent = *tit;
      break;
    }
  }
  PetscFunctionReturn(0);
}

PetscErrorCode BitLevelCouplerInterface::buidlAdjacenciesVerticesOnTets(const BitRefLevel &parent_level,Range &children,
    bool vertex_elements,const double iter_tol,const double inside_tol,bool throw_error,int verb) {
  PetscFunctionBegin;
  FieldInterface& m_field = cOre;
  //build Tree
  bool init_tree = false;
  
  //find parents of all nodes, if node has no parent then tetrahedral containing that node is searched
  //node on tetrahedra my by part of face or edge on that tetrahedral, this need to be verified
  const RefMoFEMEntity_multiIndex *refined_ptr;
  ierr = m_field.get_ref_ents(&refined_ptr); CHKERRQ(ierr);
  RefMoFEMEntity_multiIndex::index<EntType_mi_tag>::type::iterator it,hi_it;
  it = refined_ptr->get<EntType_mi_tag>().lower_bound(MBVERTEX);
  hi_it = refined_ptr->get<EntType_mi_tag>().upper_bound(MBVERTEX);

  for(;it!=hi_it;it++) {

    //entity on parent level, can be parent to yourself
    if((it->get_BitRefLevel()&parent_level).any()) continue;

    if(verb > 1) {
      cout << *it << " " << it->get_BitRefLevel() << endl;
    }

    //that vertex is on parent bit level, no need to process
    //check if vertex has a parent and parent is on parent bit level
    EntityHandle parent_ent;
    parent_ent = it->get_parent_ent();
    const RefMoFEMEntity ref_parent_ent(m_field.get_moab(),parent_ent);
    if((ref_parent_ent.get_BitRefLevel()&parent_level).any()) {
      continue;
    }

    //check if vertex is on child entities set
    EntityHandle node = it->get_ref_ent();
    if(children.find(node)==children.end()) {
      continue;
    }

    //build a boundary volume Tree
    //if(!treePtr) {
      ierr = buildTree(parent_level,verb); CHKERRQ(ierr);
      init_tree = true;
    //}

    double coords[3];
    rval = m_field.get_moab().get_coords(&node,1,coords); CHKERR_PETSC(rval);
    EntityHandle parent = 0;
    ierr = getParent(coords,parent,false,iter_tol,inside_tol,verb); CHKERRQ(ierr);
    ierr = chanegParent(refined_ptr->project<0>(it),parent,vertex_elements); CHKERRQ(ierr);
    if(throw_error && parent == 0) {
      SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCT,
  	  "tets or any other entity for node not found");
    }

  }

  if(init_tree) {
    treePtr->reset_tree();
  }
  PetscFunctionReturn(0);
}

PetscErrorCode BitLevelCouplerInterface::buidlAdjacenciesEdgesFacesVolumes(
  const BitRefLevel &parent_level,Range &children,bool elements,int verb) {
  PetscFunctionBegin;

  if(verb>2) cout << children << endl;

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
      continue;
    }

    //get ref entity iterator
    RefMoFEMEntity_multiIndex::index<Ent_mi_tag>::type::iterator it;
    it = refined_ptr->get<Ent_mi_tag>().find(*eit);
    //that entity is on parent bit level, no need to process
    if((it->get_BitRefLevel()&parent_level).any()) {
      continue;
    }

    if(verb>1) {
      cout << "before: " << *it << endl;
    }

    //check if entity has a parent and parent is on parent bit level
    EntityHandle parent_ent;
    parent_ent = it->get_parent_ent();
    const RefMoFEMEntity ref_parent_ent(m_field.get_moab(),parent_ent);
    if((ref_parent_ent.get_BitRefLevel()&parent_level).any()) {
      if(!vErify) continue;
    }

    //connectivity
    int max_dim = m_field.get_moab().dimension_from_handle(*eit);
    const EntityHandle *conn;
    int num_nodes;
    ierr = m_field.get_moab().get_connectivity(*eit,conn,num_nodes); CHKERRQ(ierr);
    conn_parents.resize(num_nodes);
    for(int nn = 0;nn<num_nodes;nn++) {
      const RefMoFEMEntity ent(m_field.get_moab(),conn[nn]);
      conn_parents[nn] = ent.get_parent_ent();
      RefMoFEMEntity_multiIndex::index<Ent_mi_tag>::type::iterator cit;
      cit = refined_ptr->get<Ent_mi_tag>().find(conn_parents[nn]);
      if(cit == refined_ptr->end()) {
	conn_parents[nn] = conn[nn];
	cit = refined_ptr->get<Ent_mi_tag>().find(conn_parents[nn]);
      }
      if((cit->get_BitRefLevel()&parent_level).none()) {
	SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCT,"parent of vertex is not on parent bit level");
      }
      int ent_dim = m_field.get_moab().dimension_from_handle(conn_parents[nn]);
      max_dim = ent_dim > max_dim ? ent_dim : max_dim;
    }

    if(verb>1) cout << "max_dim " << max_dim << endl;

    if(max_dim > 0) {

      for(;max_dim<=3;max_dim++) {
	Range parent_ents;
	rval = m_field.get_moab().get_adjacencies(
	  &*conn_parents.begin(),num_nodes,max_dim,false,parent_ents); CHKERR_PETSC(rval);
	parent_ents.erase(it->get_ref_ent());
	if(!parent_ents.empty()) {
	  ierr = chanegParent(refined_ptr->project<0>(it),*parent_ents.begin(),elements); CHKERRQ(ierr);
	  if(verb > 1) {
	    cout << "after " << *it << endl << endl;
	  }
	  break;
	}
      }

      if(!vErify && max_dim>3) {
	ierr = chanegParent(refined_ptr->project<0>(it),0,elements); CHKERRQ(ierr);
	if(verb > 1) {
	  cout << "parent not found\n";
	}	
      }

    }

  }

  PetscFunctionReturn(0);
}

PetscErrorCode BitLevelCouplerInterface::chanegParent(RefMoFEMEntity_multiIndex::iterator it,EntityHandle parent,bool element) {
  PetscFunctionBegin;

  FieldInterface& m_field = cOre;
  const RefMoFEMEntity_multiIndex *refined_ptr;
  ierr = m_field.get_ref_ents(&refined_ptr); CHKERRQ(ierr);

  if(vErify) {
    ierr = verifyParent(it,parent); CHKERRQ(ierr);
  }

  bool parent_is_set = false;
  if(element) {
    EntityHandle ent;
    ent = it->get_ref_ent();
    const RefMoFEMElement_multiIndex *refined_finite_elements_ptr;
    ierr = m_field.get_ref_finite_elements(&refined_finite_elements_ptr); CHKERRQ(ierr);
    RefMoFEMElement_multiIndex::index<Ent_mi_tag>::type::iterator eit;
    eit = refined_finite_elements_ptr->get<Ent_mi_tag>().find(ent);
    if(eit!=refined_finite_elements_ptr->get<Ent_mi_tag>().end()) {
      RefMoFEMElement_change_parent modifier(m_field.get_moab(),refined_ptr,it,parent);
      bool success;
      success = const_cast<RefMoFEMElement_multiIndex*>(refined_finite_elements_ptr)->modify(refined_finite_elements_ptr->project<0>(eit),modifier);
      if(!success) {
	SETERRQ(PETSC_COMM_SELF,MOFEM_OPERATION_UNSUCCESSFUL,"unsuccessful operation");
      }
      parent_is_set = true;
    }
  } 

  if(!parent_is_set) {
    RefMoFEMEntity_change_parent modifier(m_field.get_moab(),parent);
    bool success = const_cast<RefMoFEMEntity_multiIndex*>(refined_ptr)->modify(it,modifier);
    if(!success) {
      SETERRQ(PETSC_COMM_SELF,MOFEM_OPERATION_UNSUCCESSFUL,"unsuccessful operation");
    }
  }


  PetscFunctionReturn(0);
}

PetscErrorCode BitLevelCouplerInterface::resetParents(Range &children,bool elements,int verb) {
  PetscFunctionBegin;

  FieldInterface& m_field = cOre;

  //access to ref dofs multi-index
  const RefMoFEMEntity_multiIndex *refined_ptr;
  ierr = m_field.get_ref_ents(&refined_ptr); CHKERRQ(ierr);

  Range::iterator eit,hi_eit;
  eit = children.begin();
  hi_eit = children.end();
  for(;eit!=hi_eit;eit++) {

    //get ref entity iterator
    RefMoFEMEntity_multiIndex::index<Ent_mi_tag>::type::iterator it;
    it = refined_ptr->get<Ent_mi_tag>().find(*eit);

    //resent entity parent
    ierr = chanegParent(refined_ptr->project<0>(it),0,elements); CHKERRQ(ierr);

  }

  PetscFunctionReturn(0);
}


PetscErrorCode BitLevelCouplerInterface::verifyParent(RefMoFEMEntity_multiIndex::iterator it,EntityHandle parent) {
  PetscFunctionBegin;

  if(parent != it->get_parent_ent()) {
    SETERRQ3(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCT,"data inconsistency %lu != %lu for ent %lu",
      parent,it->get_parent_ent(),it->get_ref_ent());
  }

  PetscFunctionReturn(0);
}

PetscErrorCode BitLevelCouplerInterface::getLocCoordsOnTet(EntityHandle tet,const double *glob_coords,int verb) {
  PetscFunctionBegin;

  FieldInterface& m_field = cOre;

  int num_nodes;
  rval = m_field.get_moab().get_connectivity(tet,cOnn,num_nodes,true); CHKERR_PETSC(rval);
  rval = m_field.get_moab().get_coords(cOnn,num_nodes,cOords); CHKERR_PETSC(rval);
  double shifted_glob_coors[3];
  cblas_dcopy(3,glob_coords,1,shifted_glob_coors,1);
  for(int nn = 1;nn<4;nn++) {
    cblas_daxpy(3,-1,cOords,1,&cOords[nn*3],1);
  }
  cblas_daxpy(3,-1,cOords,1,shifted_glob_coors,1);
  for(int dd =0;dd<3;dd++) {
    cOords[dd] = 0;
  }
  ierr = ShapeDiffMBTET(diffN); CHKERRQ(ierr);
  locCoords[0] = locCoords[1] = locCoords[2] = 0;
  ierr = ShapeMBTET(N,&locCoords[0],&locCoords[1],&locCoords[2],1);; CHKERRQ(ierr);
  ierr = ShapeMBTET_inverse(N,diffN,cOords,shifted_glob_coors,locCoords); CHKERRQ(ierr);
  ierr = ShapeMBTET(N,&locCoords[0],&locCoords[1],&locCoords[2],1);; CHKERRQ(ierr);

  if(verb>1) {
    cout << "N " << N[0] << " " << N[1] << " " << N[2] << " " << N[3] << endl;
  }

  PetscFunctionReturn(0);
}


}
