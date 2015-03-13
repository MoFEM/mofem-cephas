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

PetscErrorCode BitLevelCouplerInterface::buidlDirectAdjacencies(const BitRefLevel &parent_level,Range &children) {
  PetscFunctionBegin;
  FieldInterface& m_field = cOre;
  const RefMoFEMEntity_multiIndex *refined_ptr;
  ierr = m_field.get_ref_ents(&refined_ptr); CHKERRQ(ierr);
  RefMoFEMEntity_multiIndex::iterator it = refined_ptr.end();
  for(;it!=get_refined_ptr.end();it++) {
    if(children.find(it->get_ent())==refined_ptr.end()) {
      continue;
    }
    EntityHandle parent_ent;
    parent_ent = it->get_parent_ent();
    RefMoFEMEntity_multiIndex *pit = refined_ptr.get()<Ent_mi_tag>().find(parent_ent);
    if((pit->get_BitRefLevel()&parent_level).any()) {
      continue;
    }
    const EntityHandle* conn;
    int number_nodes = 0;
    rval = m_field.get_moab().get_connectivity(it->get_ent(),conn,number_nodes);  CHKERR(rval);
    switch(it->get_ent_type()) {
      case MBEDGE: {
	Range face;
	rval = m_field.get_moab().get_adjacencies(conn,number_nodes,2,face,false,Interface::INTERSECT); CHKERR(rval);
	if(face.empty()) {

	}
	
      }
      break;
  

    }
    if(test.empty()) {

    }

  }
  PetscFunctionReturn(0);
}

PetscErrorCode BitLevelCouplerInterface::buidlDirectAdjacencies(const BitRefLevel &parent_level,const BitRefLevel &children_level) {
  PetscFunctionBegin;
  FieldInterface& m_field = cOre;
  Range &children_ents;
  ierr = m_field->get_entities_by_ref_level(children_level,BitRefLevel().set(),children_ents); CHKERRQ(ierr);
  ierr = buidlDirectAdjacencies(parent_level,children); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

}
