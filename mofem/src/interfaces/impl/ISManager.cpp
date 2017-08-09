/** \file ISManager.cpp
 * \brief Managing complexities for problem
 * \ingroup mofem_section_manager
 */

/* MoFEM is free software: you can redistribute it and/or modify it under
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
#include <fem_tools.h>

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
#include <Core.hpp>

#include <ISManager.hpp>

#include <moab/MeshTopoUtil.hpp>

namespace MoFEM {

  PetscErrorCode ISManager::queryInterface(const MOFEMuuid& uuid, UnknownInterface** iface) {
    PetscFunctionBegin;
    *iface = NULL;
    if(uuid == IDD_MOFEMISManager) {
      *iface = dynamic_cast<ISManager*>(this);
      PetscFunctionReturn(0);
    }
    if(uuid == IDD_MOFEMUnknown) {
      *iface = dynamic_cast<UnknownInterface*>(this);
      PetscFunctionReturn(0);
    }
    SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"unknown interface");
    PetscFunctionReturn(0);
  }

  ISManager::ISManager(const MoFEM::Core& core):
  cOre(const_cast<MoFEM::Core&>(core)),
  dEbug(false) {
  }
  ISManager::~ISManager() {
  }

  PetscErrorCode ISManager::sectionCreate(
    const std::string &problem_name,
    PetscSection *s,
    const RowColData row_col
  ) const {
    const MoFEM::Interface &m_field = cOre;
    const Problem *problem_ptr;
    const Field_multiIndex *fields_ptr;
    const FiniteElement_multiIndex *fe_ptr;
    PetscFunctionBegin;
    ierr = m_field.get_problem(problem_name,&problem_ptr); CHKERRQ(ierr);
    ierr = m_field.get_fields(&fields_ptr); CHKERRQ(ierr);
    ierr = m_field.get_finite_elements(&fe_ptr); CHKERRQ(ierr);
    boost::shared_ptr<NumeredDofEntity_multiIndex> dofs;
    BitFieldId fields_ids;
    switch (row_col) {
      case ROW:
      dofs = problem_ptr->numeredDofsRows;
      for(FiniteElement_multiIndex::iterator fit = fe_ptr->begin();fit!=fe_ptr->end();fit++) {
        if((fit->get()->getId()&problem_ptr->getBitFEId()).any()) {
          fields_ids |= fit->get()->getBitFieldIdRow();
        }
      }
      break;
      case COL:
      dofs = problem_ptr->numeredDofsCols;
      for(FiniteElement_multiIndex::iterator fit = fe_ptr->begin();fit!=fe_ptr->end();fit++) {
        if((fit->get()->getId()&problem_ptr->getBitFEId()).any()) {
          fields_ids |= fit->get()->getBitFieldIdCol();
        }
      }
      break;
      default:
      SETERRQ(PETSC_COMM_WORLD,MOFEM_DATA_INCONSISTENCY,"Has to be ROW or COLUMN");
    }
    // get fields names on the problem
    map<std::string,std::pair<int,int> > fields_map;
    {
      int field = 0;
      for(Field_multiIndex::iterator fit = fields_ptr->begin();fit!=fields_ptr->end();fit++) {
        if((fit->get()->getId()&fields_ids).any()) {
          fields_map[fit->get()->getName()].first = field++;
          fields_map[fit->get()->getName()].second = fit->get()->getNbOfCoeffs();
        }
      }
    }
    const int proc = m_field.get_comm_rank();
    ierr = PetscSectionCreate(PETSC_COMM_WORLD,s); CHKERRQ(ierr);
    ierr = PetscSectionSetNumFields(*s,fields_map.size()); CHKERRQ(ierr);
    for(map<std::string,std::pair<int,int> >::iterator mit = fields_map.begin();mit!=fields_map.end();mit++) {
      ierr = PetscSectionSetFieldName(*s,mit->second.first,mit->first.c_str()); CHKERRQ(ierr);
      ierr = PetscSectionSetFieldComponents(*s,mit->second.first,mit->second.second); CHKERRQ(ierr);
    }
    // determine number of points
    int nb_charts = 0;
    {
      NumeredDofEntity_multiIndex::iterator dit,hi_dit;
      dit = dofs->begin();
      hi_dit = dofs->end();
      for(;dit!=hi_dit;) {
        EntityHandle ent = dit->get()->getEnt();
        if(dit->get()->getPart()==proc&&dit->get()->getEntDofIdx()==0) {
          while(dit!=hi_dit&&ent==dit->get()->getEnt()) {
            const int nb_of_dofs_on_ent = dit->get()->getNbDofsOnEnt();
            for(int dd = 0;dd!=nb_of_dofs_on_ent;dd++,dit++) {}
          }
          ++nb_charts;
        } else {
          ++dit;
        }
      }
    }
    // get layout, i.e. chart
    PetscLayout layout;
    ierr = PetscLayoutCreate(PETSC_COMM_WORLD,&layout); CHKERRQ(ierr);
    ierr = PetscLayoutSetBlockSize(layout,1); CHKERRQ(ierr);
    ierr = PetscLayoutSetLocalSize(layout,nb_charts); CHKERRQ(ierr);
    ierr = PetscLayoutSetUp(layout); CHKERRQ(ierr);
    int rstart,rend;
    ierr = PetscLayoutGetRange(layout,&rstart,&rend); CHKERRQ(ierr);
    ierr = PetscLayoutDestroy(&layout); CHKERRQ(ierr);
    ierr = PetscSectionSetChart(*s,rstart,rend); CHKERRQ(ierr);
    // cerr << rstart << " " << rend << " " << proc << endl;
    // loop of all dofs
    {
      NumeredDofEntity_multiIndex::iterator dit,hi_dit;
      dit = dofs->begin();
      hi_dit = dofs->end();
      int point = rstart;
      for(;dit!=hi_dit;) {
        EntityHandle ent = dit->get()->getEnt();
        if(dit->get()->getPart()==proc&&dit->get()->getEntDofIdx()==0) {
          // exploit that does are continuously stored on entity
          // that includes fields
          while(dit!=hi_dit&&ent==dit->get()->getEnt()) {
            const int nb_of_dofs_on_ent = dit->get()->getNbDofsOnEnt();
            ierr = PetscSectionAddDof(*s,point,nb_of_dofs_on_ent); CHKERRQ(ierr);
            int field = fields_map.at(dit->get()->getName()).first;
            ierr = PetscSectionSetFieldDof(*s,point,field,nb_of_dofs_on_ent); CHKERRQ(ierr);
            for(int dd = 0;dd!=nb_of_dofs_on_ent;dd++,dit++) {}
            // cerr << point << endl;
          }
          ++point;
        } else {
          ++dit;
        }
      }
    }
    // cerr << "done " << proc << endl;
    ierr = PetscSectionSetUp(*s); CHKERRQ(ierr);
    // cerr << "end " << proc << endl;
    PetscFunctionReturn(0);
  }

  PetscErrorCode ISManager::isCreateProblemOrder(
    const std::string &problem,RowColData rc,int min_order,int max_order,IS *is
  ) const {
    const MoFEM::Interface &m_field = cOre;
    const Problem *problem_ptr;
    PetscFunctionBegin;
    ierr = m_field.get_problem(problem,&problem_ptr); CHKERRQ(ierr);
    typedef NumeredDofEntity_multiIndex::index<Composite_Part_And_Order_mi_tag>::type dofs_order;
    int rank = m_field.get_comm_rank();
    dofs_order::iterator it,hi_it;
    switch(rc) {
      case ROW:
      it = problem_ptr->numeredDofsRows->get<Composite_Part_And_Order_mi_tag>().lower_bound(boost::make_tuple(rank,min_order));
      hi_it = problem_ptr->numeredDofsRows->get<Composite_Part_And_Order_mi_tag>().upper_bound(boost::make_tuple(rank,max_order));
      break;
      case COL:
      it = problem_ptr->numeredDofsCols->get<Composite_Part_And_Order_mi_tag>().lower_bound(boost::make_tuple(rank,min_order));
      hi_it = problem_ptr->numeredDofsCols->get<Composite_Part_And_Order_mi_tag>().upper_bound(boost::make_tuple(rank,max_order));
      break;
      default:
      SETERRQ(PETSC_COMM_SELF,MOFEM_NOT_IMPLEMENTED,"not implemented");
    }
    NumeredDofEntity_multiIndex_petsc_local_dof_view_ordered_non_unique dof_loc_idx_view;
    for(;it!=hi_it;it++) {
      std::pair<NumeredDofEntity_multiIndex_petsc_local_dof_view_ordered_non_unique::iterator,bool> p;
      if((*it)->getPart()!=(unsigned int)rank) continue;
      p = dof_loc_idx_view.insert(*it);
    }
    NumeredDofEntity_multiIndex_petsc_local_dof_view_ordered_non_unique::iterator vit,hi_vit;
    vit = dof_loc_idx_view.begin();
    hi_vit = dof_loc_idx_view.end();
    int size = distance(vit,hi_vit);
    int *id;
    ierr = PetscMalloc(size*sizeof(int),&id); CHKERRQ(ierr);
    for(int ii = 0;vit!=hi_vit;vit++) {
      id[ii++] = (*vit)->getPetscGlobalDofIdx();
    }
    ierr = ISCreateGeneral(PETSC_COMM_WORLD,size,id,PETSC_OWN_POINTER,is); CHKERRQ(ierr);
    PetscFunctionReturn(0);
  }

  PetscErrorCode ISManager::isCreateProblemFieldAndRank(
    const std::string &problem,
    RowColData rc,
    const std::string &field,
    int min_coeff_idx,
    int max_coeff_idx,
    IS *is
  ) const {
    const MoFEM::Interface &m_field = cOre;
    const Problem *problem_ptr;
    PetscFunctionBegin;
    ierr = m_field.get_problem(problem,&problem_ptr); CHKERRQ(ierr);
    typedef NumeredDofEntity_multiIndex::index<Composite_Name_Part_And_CoeffIdx_mi_tag>::type DofsByNamePartAndCoeffIdx;
    int rank = m_field.get_comm_rank();
    DofsByNamePartAndCoeffIdx::iterator it,hi_it;
    switch(rc) {
      case ROW:
      it = problem_ptr->numeredDofsRows->get<Composite_Name_Part_And_CoeffIdx_mi_tag>().lower_bound(
        boost::make_tuple(field,rank,min_coeff_idx)
      );
      hi_it = problem_ptr->numeredDofsRows->get<Composite_Name_Part_And_CoeffIdx_mi_tag>().upper_bound(
        boost::make_tuple(field,rank,max_coeff_idx)
      );
      break;
      case COL:
      it = problem_ptr->numeredDofsCols->get<Composite_Name_Part_And_CoeffIdx_mi_tag>().lower_bound(
        boost::make_tuple(field,rank,min_coeff_idx)
      );
      hi_it = problem_ptr->numeredDofsCols->get<Composite_Name_Part_And_CoeffIdx_mi_tag>().upper_bound(
        boost::make_tuple(field,rank,max_coeff_idx)
      );
      break;
      default:
      SETERRQ(PETSC_COMM_SELF,MOFEM_NOT_IMPLEMENTED,"not implemented");
    }

    // Sort by local index
    NumeredDofEntity_multiIndex_petsc_local_dof_view_ordered_non_unique dof_loc_idx_view;
    dof_loc_idx_view.insert(it,hi_it);
    NumeredDofEntity_multiIndex_petsc_local_dof_view_ordered_non_unique::iterator vit,hi_vit;
    vit = dof_loc_idx_view.begin();
    hi_vit = dof_loc_idx_view.end();

    // create IS
    int size = distance(it,hi_it);
    int *id;
    ierr = PetscMalloc(size*sizeof(int),&id); CHKERRQ(ierr);
    for(int ii = 0;vit!=hi_vit;vit++) {
      id[ii++] = (*vit)->getPetscGlobalDofIdx();
    }
    ierr = ISCreateGeneral(PETSC_COMM_WORLD,size,id,PETSC_OWN_POINTER,is); CHKERRQ(ierr);
    PetscFunctionReturn(0);
  }

  PetscErrorCode ISManager::isCreateFromProblemFieldToOtherProblemField(
    const std::string &x_problem,
    const std::string &x_field_name,RowColData x_rc,
    const std::string &y_problem,
    const std::string &y_field_name,RowColData y_rc,
    std::vector<int> &idx,std::vector<int> &idy
  ) const {
    const MoFEM::Interface &m_field = cOre;
    const Problem *px_ptr;
    const Problem *py_ptr;
    PetscFunctionBegin;
    ierr = m_field.get_problem(x_problem,&px_ptr); CHKERRQ(ierr);
    ierr = m_field.get_problem(y_problem,&py_ptr); CHKERRQ(ierr);
    NumeredDofEntityByLocalIdx::iterator y_dit,hi_y_dit;
    switch (y_rc) {
      case ROW:
        y_dit = py_ptr->numeredDofsRows->get<PetscLocalIdx_mi_tag>().lower_bound(0);
        hi_y_dit = py_ptr->numeredDofsRows->get<PetscLocalIdx_mi_tag>().
        upper_bound(py_ptr->getNbLocalDofsRow()-1);
        break;
      case COL:
        y_dit = py_ptr->numeredDofsCols->get<PetscLocalIdx_mi_tag>().lower_bound(0);
        hi_y_dit = py_ptr->numeredDofsCols->get<PetscLocalIdx_mi_tag>().
        upper_bound(py_ptr->getNbLocalDofsCol()-1);
        break;
      default:
       SETERRQ(PETSC_COMM_SELF,MOFEM_NOT_IMPLEMENTED,"only makes sense for ROWS and COLS");
    }
    typedef NumeredDofEntity_multiIndex::index<Composite_Name_And_Ent_And_EntDofIdx_mi_tag>::type DofsByNameAndEntDofIdx;
    const DofsByNameAndEntDofIdx* x_numered_dofs_by_ent_name_dof;
    switch (x_rc) {
      case ROW:
        x_numered_dofs_by_ent_name_dof =
        &(px_ptr->numeredDofsRows->get<Composite_Name_And_Ent_And_EntDofIdx_mi_tag>());
        break;
      case COL:
        x_numered_dofs_by_ent_name_dof =
        &(px_ptr->numeredDofsCols->get<Composite_Name_And_Ent_And_EntDofIdx_mi_tag>());
        break;
      default:
       SETERRQ(PETSC_COMM_SELF,MOFEM_NOT_IMPLEMENTED,"only makes sense for ROWS and COLS");
    }
    std::map<int,int> global_dofs_map;
    for(;y_dit!=hi_y_dit;y_dit++) {
      if((*y_dit)->getPart()!=(unsigned int)m_field.get_comm_rank()) {
        SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"data inconsistency");
      }
      if((*y_dit)->getName()!=y_field_name) continue;
      DofsByNameAndEntDofIdx::iterator x_dit;
      x_dit = x_numered_dofs_by_ent_name_dof->find(
        boost::make_tuple(x_field_name,(*y_dit)->getEnt(),(*y_dit)->getEntDofIdx())
      );
      if(x_dit==x_numered_dofs_by_ent_name_dof->end()) continue;
      global_dofs_map[(*x_dit)->getPetscGlobalDofIdx()] = (*y_dit)->getPetscGlobalDofIdx();
    }
    idx.resize(global_dofs_map.size());
    idy.resize(global_dofs_map.size());
    {
      std::vector<int>::iterator ix,iy;
      ix = idx.begin();
      iy = idy.begin();
      map<int,int>::iterator mit = global_dofs_map.begin();
      for(;mit!=global_dofs_map.end();mit++,ix++,iy++) {
        *ix = mit->first;
        *iy = mit->second;
      }
    }
    PetscFunctionReturn(0);
  }

  PetscErrorCode ISManager::isCreateFromProblemFieldToOtherProblemField(
    const std::string &x_problem,const std::string &x_field_name,RowColData x_rc,
    const std::string &y_problem,const std::string &y_field_name,RowColData y_rc,
    IS *ix,IS *iy
  ) const {
    PetscFunctionBegin;
    std::vector<int> idx(0),idy(0);
    ierr = isCreateFromProblemFieldToOtherProblemField(
      x_problem,x_field_name,x_rc,
      y_problem,y_field_name,y_rc,idx,idy
    ); CHKERRQ(ierr);
    if(ix!=PETSC_NULL) {
      ierr = ISCreateGeneral(
        PETSC_COMM_WORLD,idx.size(),&idx[0],PETSC_COPY_VALUES,ix
      ); CHKERRQ(ierr);
    }
    ierr = ISCreateGeneral(
      PETSC_COMM_WORLD,idy.size(),&idy[0],PETSC_COPY_VALUES,iy
    ); CHKERRQ(ierr);
    if(dEbug) {
      ISView(*ix,PETSC_VIEWER_STDOUT_WORLD);
      ISView(*iy,PETSC_VIEWER_STDOUT_WORLD);
    }
    PetscFunctionReturn(0);
  }

  PetscErrorCode ISManager::isCreateFromProblemToOtherProblem(
    const std::string &x_problem,RowColData x_rc,
    const std::string &y_problem,RowColData y_rc,
    std::vector<int> &idx,std::vector<int> &idy
  ) const {
    const MoFEM::Interface &m_field = cOre;
    const Problem *px_ptr;
    const Problem *py_ptr;
    PetscFunctionBegin;
    ierr = m_field.get_problem(x_problem,&px_ptr); CHKERRQ(ierr);
    ierr = m_field.get_problem(y_problem,&py_ptr); CHKERRQ(ierr);
    NumeredDofEntityByLocalIdx::iterator y_dit,hi_y_dit;
    switch (y_rc) {
      case ROW:
        y_dit = py_ptr->numeredDofsRows->get<PetscLocalIdx_mi_tag>().lower_bound(0);
        hi_y_dit = py_ptr->numeredDofsRows->get<PetscLocalIdx_mi_tag>().lower_bound(py_ptr->getNbLocalDofsRow()); // should be lower
        break;
      case COL:
        y_dit = py_ptr->numeredDofsCols->get<PetscLocalIdx_mi_tag>().lower_bound(0);
        hi_y_dit = py_ptr->numeredDofsCols->get<PetscLocalIdx_mi_tag>().lower_bound(py_ptr->getNbLocalDofsCol()); // should be lower
        break;
      default:
       SETERRQ(PETSC_COMM_SELF,MOFEM_NOT_IMPLEMENTED,"not implemented");
    }
    const NumeredDofEntityByUId* x_numered_dofs_by_uid;
    switch (x_rc) {
      case ROW:
        x_numered_dofs_by_uid = &(px_ptr->numeredDofsRows->get<Unique_mi_tag>());
        break;
      case COL:
        x_numered_dofs_by_uid = &(px_ptr->numeredDofsCols->get<Unique_mi_tag>());
        break;
      default:
       SETERRQ(PETSC_COMM_SELF,MOFEM_NOT_IMPLEMENTED,"not implemented");
    }
    for(;y_dit!=hi_y_dit;y_dit++) {
      if((*y_dit)->getPart()!=(unsigned int)m_field.get_comm_rank()) {
        SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"data inconsistency");
      }
      NumeredDofEntityByUId::iterator x_dit;
      x_dit = x_numered_dofs_by_uid->find((*y_dit)->getGlobalUniqueId());
      if(x_dit==x_numered_dofs_by_uid->end()) continue;
      idx.push_back((*x_dit)->getPetscGlobalDofIdx());
      idy.push_back((*y_dit)->getPetscGlobalDofIdx());
    }
    PetscFunctionReturn(0);
  }

  PetscErrorCode ISManager::isCreateFromProblemToOtherProblem(
    const std::string &x_problem,RowColData x_rc,
    const std::string &y_problem,RowColData y_rc,
    IS *ix,IS *iy
  ) const {
    PetscFunctionBegin;
    std::vector<int> idx(0),idy(0);
    ierr = isCreateFromProblemToOtherProblem(x_problem,x_rc,y_problem,y_rc,idx,idy); CHKERRQ(ierr);
    ierr = ISCreateGeneral(PETSC_COMM_WORLD,idx.size(),&idx[0],PETSC_COPY_VALUES,ix); CHKERRQ(ierr);
    ierr = ISCreateGeneral(PETSC_COMM_WORLD,idy.size(),&idy[0],PETSC_COPY_VALUES,iy); CHKERRQ(ierr);
    if(dEbug) {
      ISView(*ix,PETSC_VIEWER_STDOUT_WORLD);
      ISView(*iy,PETSC_VIEWER_STDOUT_WORLD);
    }
    PetscFunctionReturn(0);
  }

}
