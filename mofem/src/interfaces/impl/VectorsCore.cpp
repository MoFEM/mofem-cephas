/** \file Vectors.cpp
 * \brief Managing Vec, IS and Scatter
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

namespace MoFEM {

// const static int debug = 1;

PetscErrorCode Core::VecCreateSeq(const std::string &name,RowColData rc,Vec *V) const {
  PetscErrorCode ierr;
  PetscFunctionBegin;
  typedef MoFEMProblem_multiIndex::index<Problem_mi_tag>::type ProblemsByName;
  const ProblemsByName &problems_set = pRoblems.get<Problem_mi_tag>();
  ProblemsByName::iterator p_miit = problems_set.find(name);
  if(p_miit==problems_set.end()) {
    SETERRQ1(
      PETSC_COMM_SELF,
      MOFEM_DATA_INCONSISTENCY,
      "No such problem %s (top tip check spelling)",
      name.c_str()
    );
  }
  DofIdx nb_local_dofs,nb_ghost_dofs;
  switch (rc) {
    case ROW:
      nb_local_dofs = p_miit->getNbLocalDofsRow();
      nb_ghost_dofs = p_miit->getNbGhostDofsRow();
      break;
    case COL:
      nb_local_dofs = p_miit->getNbLocalDofsCol();
      nb_ghost_dofs = p_miit->getNbGhostDofsCol();
      break;
    default:
     SETERRQ(PETSC_COMM_SELF,MOFEM_NOT_IMPLEMENTED,"Not implemented");
  }
  ierr = ::VecCreateSeq(PETSC_COMM_SELF,nb_local_dofs+nb_ghost_dofs,V); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
PetscErrorCode Core::VecCreateGhost(const std::string &name,RowColData rc,Vec *V) const {
  PetscErrorCode ierr;
  PetscFunctionBegin;
  typedef MoFEMProblem_multiIndex::index<Problem_mi_tag>::type ProblemsByName;
  const ProblemsByName &problems_set = pRoblems.get<Problem_mi_tag>();
  ProblemsByName::iterator p_miit = problems_set.find(name);
  if(p_miit==problems_set.end()) {
    SETERRQ1(
      PETSC_COMM_SELF,
      MOFEM_DATA_INCONSISTENCY,
      "No such problem %s (top tip check spelling)",
      name.c_str()
    );
  }
  DofIdx nb_dofs,nb_local_dofs,nb_ghost_dofs;
  NumeredDofEntityByLocalIdx *dofs;
  switch (rc) {
    case ROW:
      nb_dofs = p_miit->getNbDofsRow();
      nb_local_dofs = p_miit->getNbLocalDofsRow();
      nb_ghost_dofs = p_miit->getNbGhostDofsRow();
      dofs = const_cast<NumeredDofEntityByLocalIdx*>(&p_miit->numered_dofs_rows->get<PetscLocalIdx_mi_tag>());
      break;
    case COL:
      nb_dofs = p_miit->getNbDofsCol();
      nb_local_dofs = p_miit->getNbLocalDofsCol();
      nb_ghost_dofs = p_miit->getNbGhostDofsCol();
      dofs = const_cast<NumeredDofEntityByLocalIdx*>(&p_miit->numered_dofs_cols->get<PetscLocalIdx_mi_tag>());
      break;
    default:
     SETERRQ(PETSC_COMM_SELF,MOFEM_NOT_IMPLEMENTED,"not implemented");
  }
  NumeredDofEntityByLocalIdx::iterator miit = dofs->lower_bound(nb_local_dofs);
  NumeredDofEntityByLocalIdx::iterator hi_miit = dofs->upper_bound(nb_local_dofs+nb_ghost_dofs);
  int count = distance(miit,hi_miit);
  if(count != nb_ghost_dofs) {
    SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"data inconsistency");
  }
  std::vector<DofIdx> ghost_idx(count);
  std::vector<DofIdx>::iterator vit = ghost_idx.begin();
  for(;miit!=hi_miit;miit++,vit++) {
    *vit = (*miit)->petscGloablDofIdx;
  }
  ierr = ::VecCreateGhost(comm,nb_local_dofs,nb_dofs,nb_ghost_dofs,&ghost_idx[0],V); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
PetscErrorCode Core::ISCreateProblemOrder(
  const std::string &problem,RowColData rc,int min_order,int max_order,IS *is,int verb
) const {
  PetscErrorCode ierr;
  PetscFunctionBegin;
  if(verb==-1) verb = verbose;
  typedef MoFEMProblem_multiIndex::index<Problem_mi_tag>::type ProblemsByName;
  const ProblemsByName &problems_set = pRoblems.get<Problem_mi_tag>();
  ProblemsByName::iterator p = problems_set.find(problem);
  if(p==problems_set.end()) SETERRQ1(PETSC_COMM_SELF,1,"no such problem %s (top tip check spelling)",problem.c_str());
  typedef NumeredDofEntity_multiIndex::index<Composite_Part_And_Order_mi_tag>::type dofs_order;
  dofs_order::iterator it,hi_it;
  switch(rc) {
    case ROW:
    it = p->numered_dofs_rows->get<Composite_Part_And_Order_mi_tag>().lower_bound(boost::make_tuple(rAnk,min_order));
    hi_it = p->numered_dofs_rows->get<Composite_Part_And_Order_mi_tag>().upper_bound(boost::make_tuple(rAnk,max_order));
    break;
    case COL:
    it = p->numered_dofs_cols->get<Composite_Part_And_Order_mi_tag>().lower_bound(boost::make_tuple(rAnk,min_order));
    hi_it = p->numered_dofs_cols->get<Composite_Part_And_Order_mi_tag>().upper_bound(boost::make_tuple(rAnk,max_order));
    break;
    default:
     SETERRQ(PETSC_COMM_SELF,MOFEM_NOT_IMPLEMENTED,"not implemented");
  }
  NumeredDofEntity_multiIndex_petsc_local_dof_view_ordered_non_unique dof_loc_idx_view;
  for(;it!=hi_it;it++) {
    std::pair<NumeredDofEntity_multiIndex_petsc_local_dof_view_ordered_non_unique::iterator,bool> p;
    if((*it)->getPart()!=(unsigned int)rAnk) continue;
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
  ierr = ISCreateGeneral(comm,size,id,PETSC_OWN_POINTER,is); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
PetscErrorCode Core::ISCreateProblemFieldAndRank(
  const std::string &problem,
  RowColData rc,
  const std::string &field,
  int min_coeff_idx,
  int max_coeff_idx,
  IS *is,
  int verb
) const {
  PetscErrorCode ierr;
  PetscFunctionBegin;
  if(verb==-1) verb = verbose;
  typedef MoFEMProblem_multiIndex::index<Problem_mi_tag>::type ProblemsByName;
  const ProblemsByName &problems_set = pRoblems.get<Problem_mi_tag>();
  ProblemsByName::iterator p = problems_set.find(problem);
  if(p==problems_set.end()) SETERRQ1(PETSC_COMM_SELF,1,"no such problem %s (top tip check spelling)",problem.c_str());
  typedef NumeredDofEntity_multiIndex::index<Composite_Name_Part_And_CoeffIdx_mi_tag>::type dofs_by_name_and_rank;
  dofs_by_name_and_rank::iterator it,hi_it;
  switch(rc) {
    case ROW:
    it = p->numered_dofs_rows->get<Composite_Name_Part_And_CoeffIdx_mi_tag>().lower_bound(boost::make_tuple(field,rAnk,min_coeff_idx));
    hi_it = p->numered_dofs_rows->get<Composite_Name_Part_And_CoeffIdx_mi_tag>().upper_bound(boost::make_tuple(field,rAnk,max_coeff_idx));
    break;
    case COL:
    it = p->numered_dofs_cols->get<Composite_Name_Part_And_CoeffIdx_mi_tag>().lower_bound(boost::make_tuple(field,rAnk,min_coeff_idx));
    hi_it = p->numered_dofs_cols->get<Composite_Name_Part_And_CoeffIdx_mi_tag>().upper_bound(boost::make_tuple(field,rAnk,max_coeff_idx));
    break;
    default:
     SETERRQ(PETSC_COMM_SELF,MOFEM_NOT_IMPLEMENTED,"not implemented");
  }

  NumeredDofEntity_multiIndex_petsc_local_dof_view_ordered_non_unique dof_loc_idx_view;
  for(;it!=hi_it;it++) {
    std::pair<NumeredDofEntity_multiIndex_petsc_local_dof_view_ordered_non_unique::iterator,bool> p;
    if((*it)->getPart()!=(unsigned int)rAnk) continue;
    if((*it)->getNameRef() != field) continue;
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

  ierr = ISCreateGeneral(comm,size,id,PETSC_OWN_POINTER,is); CHKERRQ(ierr);

  PetscFunctionReturn(0);
}
PetscErrorCode Core::ISCreateFromProblemFieldToOtherProblemField(
  const std::string &x_problem,const std::string &x_field_name,RowColData x_rc,
  const std::string &y_problem,const std::string &y_field_name,RowColData y_rc,
  std::vector<int> &idx,std::vector<int> &idy,int verb
) const {
  //PetscErrorCode ierr;
  PetscFunctionBegin;
  if(verb==-1) verb = verbose;
  typedef MoFEMProblem_multiIndex::index<Problem_mi_tag>::type ProblemsByName;
  const ProblemsByName &problems_set = pRoblems.get<Problem_mi_tag>();
  ProblemsByName::iterator p_x = problems_set.find(x_problem);
  if(p_x==problems_set.end()) {
    SETERRQ1(
      PETSC_COMM_SELF,
      MOFEM_DATA_INCONSISTENCY,
      "no such problem %s (top tip check spelling)",
      x_problem.c_str());
  }
  ProblemsByName::iterator p_y = problems_set.find(y_problem);
  if(p_y==problems_set.end()) {
    SETERRQ1(
      PETSC_COMM_SELF,
      MOFEM_DATA_INCONSISTENCY,
      "no such problem %s (top tip check spelling)",
      y_problem.c_str()
    );
  }
  NumeredDofEntityByLocalIdx::iterator y_dit,hi_y_dit;
  switch (y_rc) {
    case ROW:
      y_dit = p_y->numered_dofs_rows->get<PetscLocalIdx_mi_tag>().lower_bound(0);
      hi_y_dit = p_y->numered_dofs_rows->get<PetscLocalIdx_mi_tag>().
      upper_bound(p_y->getNbLocalDofsRow()-1);
      break;
    case COL:
      y_dit = p_y->numered_dofs_cols->get<PetscLocalIdx_mi_tag>().lower_bound(0);
      hi_y_dit = p_y->numered_dofs_cols->get<PetscLocalIdx_mi_tag>().
      upper_bound(p_y->getNbLocalDofsCol()-1);
      break;
    default:
     SETERRQ(PETSC_COMM_SELF,MOFEM_NOT_IMPLEMENTED,"only makes sense for ROWS and COLS");
  }
  typedef NumeredDofEntity_multiIndex::index<Composite_Name_And_Ent_And_EntDofIdx_mi_tag>::type DofsByNameAndEntDofIdx;
  const DofsByNameAndEntDofIdx* x_numered_dofs_by_ent_name_dof;
  switch (x_rc) {
    case ROW:
      x_numered_dofs_by_ent_name_dof =
      &(p_x->numered_dofs_rows->get<Composite_Name_And_Ent_And_EntDofIdx_mi_tag>());
      break;
    case COL:
      x_numered_dofs_by_ent_name_dof =
      &(p_x->numered_dofs_cols->get<Composite_Name_And_Ent_And_EntDofIdx_mi_tag>());
      break;
    default:
     SETERRQ(PETSC_COMM_SELF,MOFEM_NOT_IMPLEMENTED,"only makes sense for ROWS and COLS");
  }
  std::map<int,int> global_dofs_map;
  for(;y_dit!=hi_y_dit;y_dit++) {
    if((*y_dit)->getPart()!=(unsigned int)rAnk) {
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
PetscErrorCode Core::ISCreateFromProblemFieldToOtherProblemField(
  const std::string &x_problem,const std::string &x_field_name,RowColData x_rc,
  const std::string &y_problem,const std::string &y_field_name,RowColData y_rc,
  IS *ix,IS *iy,int verb
) const {
  PetscErrorCode ierr;
  PetscFunctionBegin;

  std::vector<int> idx(0),idy(0);
  ierr = ISCreateFromProblemFieldToOtherProblemField(
    x_problem,x_field_name,x_rc,
    y_problem,y_field_name,y_rc,idx,idy,
    verb
  ); CHKERRQ(ierr);

  if(ix!=PETSC_NULL) {
    ierr = ISCreateGeneral(comm,idx.size(),&idx[0],PETSC_COPY_VALUES,ix); CHKERRQ(ierr);
  }
  ierr = ISCreateGeneral(comm,idy.size(),&idy[0],PETSC_COPY_VALUES,iy); CHKERRQ(ierr);
  if(verb>2) {
    ISView(*ix,PETSC_VIEWER_STDOUT_WORLD);
    ISView(*iy,PETSC_VIEWER_STDOUT_WORLD);
  }
  PetscFunctionReturn(0);
}
PetscErrorCode Core::VecScatterCreate(
  Vec xin,const std::string &x_problem,const std::string &x_field_name,RowColData x_rc,
  Vec yin,const std::string &y_problem,const std::string &y_field_name,RowColData y_rc,
  VecScatter *newctx,int verb
) const {
  PetscErrorCode ierr;
  PetscFunctionBegin;
  if(verb==-1) verb = verbose;
  std::vector<int> idx(0),idy(0);
  ierr = ISCreateFromProblemFieldToOtherProblemField(
    x_problem,x_field_name,x_rc,y_problem,y_field_name,y_rc,
    idx,idy,verb
  ); CHKERRQ(ierr);
  IS ix,iy;
  ierr = ISCreateGeneral(comm,idx.size(),&idx[0],PETSC_USE_POINTER,&ix); CHKERRQ(ierr);
  ierr = ISCreateGeneral(comm,idy.size(),&idy[0],PETSC_USE_POINTER,&iy); CHKERRQ(ierr);
  ierr = ::VecScatterCreate(xin,ix,yin,iy,newctx); CHKERRQ(ierr);
  ierr = ISDestroy(&ix); CHKERRQ(ierr);
  ierr = ISDestroy(&iy); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
PetscErrorCode Core::ISCreateFromProblemToOtherProblem(
  const std::string &x_problem,
  RowColData x_rc,
  const std::string &y_problem,
  RowColData y_rc,
  std::vector<int> &idx,
  std::vector<int> &idy,
  int verb
) const {
  //PetscErrorCode ierr;
  PetscFunctionBegin;
  if(verb==-1) verb = verbose;
  typedef MoFEMProblem_multiIndex::index<Problem_mi_tag>::type ProblemsByName;
  const ProblemsByName &problems_set = pRoblems.get<Problem_mi_tag>();
  ProblemsByName::iterator p_x = problems_set.find(x_problem);
  if(p_x==problems_set.end()) SETERRQ1(PETSC_COMM_SELF,1,"no such problem %s (top tip check spelling)",x_problem.c_str());
  ProblemsByName::iterator p_y = problems_set.find(y_problem);
  if(p_y==problems_set.end()) SETERRQ1(PETSC_COMM_SELF,1,"no such problem %s (top tip check spelling)",y_problem.c_str());
  NumeredDofEntityByLocalIdx::iterator y_dit,hi_y_dit;
  switch (y_rc) {
    case ROW:
      y_dit = p_y->numered_dofs_rows->get<PetscLocalIdx_mi_tag>().lower_bound(0);
      hi_y_dit = p_y->numered_dofs_rows->get<PetscLocalIdx_mi_tag>().lower_bound(p_y->getNbLocalDofsRow()); // should be lower
      break;
    case COL:
      y_dit = p_y->numered_dofs_cols->get<PetscLocalIdx_mi_tag>().lower_bound(0);
      hi_y_dit = p_y->numered_dofs_cols->get<PetscLocalIdx_mi_tag>().lower_bound(p_y->getNbLocalDofsCol()); // should be lower
      break;
    default:
     SETERRQ(PETSC_COMM_SELF,MOFEM_NOT_IMPLEMENTED,"not implemented");
  }
  const NumeredDofEntityByUId* x_numered_dofs_by_uid;
  switch (x_rc) {
    case ROW:
      x_numered_dofs_by_uid = &(p_x->numered_dofs_rows->get<Unique_mi_tag>());
      break;
    case COL:
      x_numered_dofs_by_uid = &(p_x->numered_dofs_cols->get<Unique_mi_tag>());
      break;
    default:
     SETERRQ(PETSC_COMM_SELF,MOFEM_NOT_IMPLEMENTED,"not implemented");
  }
  for(;y_dit!=hi_y_dit;y_dit++) {
    if((*y_dit)->getPart()!=(unsigned int)rAnk) SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"data inconsistency");
    NumeredDofEntityByUId::iterator x_dit;
    x_dit = x_numered_dofs_by_uid->find((*y_dit)->getGlobalUniqueId());
    if(x_dit==x_numered_dofs_by_uid->end()) continue;
    idx.push_back((*x_dit)->getPetscGlobalDofIdx());
    idy.push_back((*y_dit)->getPetscGlobalDofIdx());
  }
  PetscFunctionReturn(0);
}
PetscErrorCode Core::ISCreateFromProblemToOtherProblem(
  const std::string &x_problem,
  RowColData x_rc,
  const std::string &y_problem,
  RowColData y_rc,
  IS *ix,
  IS *iy,
  int verb
) const {
  PetscErrorCode ierr;
  PetscFunctionBegin;
  std::vector<int> idx(0),idy(0);
  ierr = ISCreateFromProblemToOtherProblem(x_problem,x_rc,y_problem,y_rc,idx,idy,verb); CHKERRQ(ierr);
  ierr = ISCreateGeneral(comm,idx.size(),&idx[0],PETSC_COPY_VALUES,ix); CHKERRQ(ierr);
  ierr = ISCreateGeneral(comm,idy.size(),&idy[0],PETSC_COPY_VALUES,iy); CHKERRQ(ierr);
  if(verb>2) {
    ISView(*ix,PETSC_VIEWER_STDOUT_WORLD);
    ISView(*iy,PETSC_VIEWER_STDOUT_WORLD);
  }
  PetscFunctionReturn(0);
}
PetscErrorCode Core::VecScatterCreate(
  Vec xin,
  const std::string &x_problem,
  RowColData x_rc,
  Vec yin,
  const std::string &y_problem,
  RowColData y_rc,
  VecScatter *newctx,
  int verb
) const {
  PetscErrorCode ierr;
  PetscFunctionBegin;
  if(verb==-1) verb = verbose;
  std::vector<int> idx(0),idy(0);
  ierr = ISCreateFromProblemToOtherProblem(x_problem,x_rc,y_problem,y_rc,idx,idy,verb); CHKERRQ(ierr);
  IS ix,iy;
  ierr = ISCreateGeneral(comm,idx.size(),&idx[0],PETSC_USE_POINTER,&ix); CHKERRQ(ierr);
  ierr = ISCreateGeneral(comm,idy.size(),&idy[0],PETSC_USE_POINTER,&iy); CHKERRQ(ierr);
  if(verb>2) {
    ISView(ix,PETSC_VIEWER_STDOUT_WORLD);
    ISView(iy,PETSC_VIEWER_STDOUT_WORLD);
  }
  ierr = ::VecScatterCreate(xin,ix,yin,iy,newctx); CHKERRQ(ierr);
  ierr = ISDestroy(&ix); CHKERRQ(ierr);
  ierr = ISDestroy(&iy); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
PetscErrorCode Core::set_local_ghost_vector(
  const MoFEMProblem *problem_ptr,RowColData rc,Vec V,InsertMode mode,ScatterMode scatter_mode
) const {
  PetscErrorCode ierr;
  PetscFunctionBegin;
  NumeredDofEntityByLocalIdx *dofs;
  DofIdx nb_local_dofs,nb_ghost_dofs;
  switch (rc) {
    case ROW:
      nb_local_dofs = problem_ptr->getNbLocalDofsRow();
      nb_ghost_dofs = problem_ptr->getNbGhostDofsRow();
      dofs = const_cast<NumeredDofEntityByLocalIdx*>(&problem_ptr->numered_dofs_rows->get<PetscLocalIdx_mi_tag>());
      break;
    case COL:
      nb_local_dofs = problem_ptr->getNbLocalDofsCol();
      nb_ghost_dofs = problem_ptr->getNbGhostDofsCol();
      dofs = const_cast<NumeredDofEntityByLocalIdx*>(&problem_ptr->numered_dofs_cols->get<PetscLocalIdx_mi_tag>());
      break;
    default:
     SETERRQ(PETSC_COMM_SELF,MOFEM_NOT_IMPLEMENTED,"not implemented");
  }
  Vec Vlocal;
  ierr = VecGhostGetLocalForm(V,&Vlocal); CHKERRQ(ierr);
  int size;
  ierr = VecGetLocalSize(V,&size); CHKERRQ(ierr);
  if(size!=nb_local_dofs) {
    SETERRQ(
      PETSC_COMM_SELF,
      MOFEM_DATA_INCONSISTENCY,
      "data inconsistency: check ghost vector, problem with nb. of local nodes"
    );
  }
  ierr = VecGetLocalSize(Vlocal,&size); CHKERRQ(ierr);
  if(size!=nb_local_dofs+nb_ghost_dofs) {
    SETERRQ(
      PETSC_COMM_SELF,
      MOFEM_DATA_INCONSISTENCY,
      "data inconsistency: check ghost vector, problem with nb. of ghost nodes"
    );
  }
  NumeredDofEntityByLocalIdx::iterator miit = dofs->lower_bound(0);
  NumeredDofEntityByLocalIdx::iterator hi_miit = dofs->upper_bound(nb_local_dofs+nb_ghost_dofs);
  PetscScalar *array;
  VecGetArray(Vlocal,&array);
  DofIdx ii = 0;
  switch (scatter_mode) {
    case SCATTER_FORWARD:
    switch (mode) {
      case INSERT_VALUES:
      for(;miit!=hi_miit;miit++,ii++) array[ii] = (*miit)->getFieldData();
      break;
      case ADD_VALUES:
      for(;miit!=hi_miit;miit++,ii++) array[ii] += (*miit)->getFieldData();
      break;
      default:
      SETERRQ(PETSC_COMM_SELF,MOFEM_NOT_IMPLEMENTED,"not implemented");
    }
    break;
    case SCATTER_REVERSE:
    switch (mode) {
      case INSERT_VALUES:
      for(;miit!=hi_miit;miit++,ii++) {
        //std::cerr << *miit << std::endl;
        //std::cerr << array[ii] << std::endl;
        (*miit)->getFieldData() = array[ii];
      }
      break;
      case ADD_VALUES:
      for(;miit!=hi_miit;miit++,ii++) (*miit)->getFieldData() += array[ii];
      break;
      default:
      SETERRQ(PETSC_COMM_SELF,MOFEM_NOT_IMPLEMENTED,"not implemented");
    }
    break;
    default:
    SETERRQ(PETSC_COMM_SELF,MOFEM_NOT_IMPLEMENTED,"not implemented");
  }
  VecRestoreArray(Vlocal,&array);
  VecDestroy(&Vlocal);
  PetscFunctionReturn(0);
}
PetscErrorCode Core::set_local_ghost_vector(
  const std::string &name,RowColData rc,Vec V,InsertMode mode,ScatterMode scatter_mode
) const {
  PetscErrorCode ierr;
  PetscFunctionBegin;
  typedef MoFEMProblem_multiIndex::index<Problem_mi_tag>::type ProblemsByName;
  const ProblemsByName &problems_set = pRoblems.get<Problem_mi_tag>();
  ProblemsByName::iterator p_miit = problems_set.find(name);
  if(p_miit==problems_set.end()) {
    SETERRQ1(
      PETSC_COMM_SELF,
      MOFEM_DATA_INCONSISTENCY,
      "problem < %s > not found (top tip: check spelling)",
      name.c_str()
    );
  }
  ierr = set_local_ghost_vector(&*p_miit,rc,V,mode,scatter_mode); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
PetscErrorCode Core::set_global_ghost_vector(
  const MoFEMProblem *problem_ptr,RowColData rc,Vec V,InsertMode mode,ScatterMode scatter_mode
) const {
  PetscErrorCode ierr;
  PetscFunctionBegin;
  typedef NumeredDofEntity_multiIndex::index<PetscGlobalIdx_mi_tag>::type DofsByGlobalIdx;
  DofsByGlobalIdx *dofs;
  DofIdx nb_dofs;
  switch (rc) {
    case ROW:
      nb_dofs = problem_ptr->getNbDofsRow();
      dofs = const_cast<DofsByGlobalIdx*>(&problem_ptr->numered_dofs_rows->get<PetscGlobalIdx_mi_tag>());
      break;
    case COL:
      nb_dofs = problem_ptr->getNbDofsCol();
      dofs = const_cast<DofsByGlobalIdx*>(&problem_ptr->numered_dofs_cols->get<PetscGlobalIdx_mi_tag>());
      break;
    default:
     SETERRQ(PETSC_COMM_SELF,MOFEM_NOT_IMPLEMENTED,"not implemented");
  }
  DofsByGlobalIdx::iterator miit = dofs->lower_bound(0);
  DofsByGlobalIdx::iterator hi_miit = dofs->upper_bound(nb_dofs);
  switch (scatter_mode) {
    case SCATTER_REVERSE: {
      VecScatter ctx;
      Vec V_glob;
      ierr = VecScatterCreateToAll(V,&ctx,&V_glob); CHKERRQ(ierr);
      ierr = VecScatterBegin(ctx,V,V_glob,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
      ierr = VecScatterEnd(ctx,V,V_glob,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
      int size;
      ierr = VecGetSize(V_glob,&size); CHKERRQ(ierr);
      if(size!=nb_dofs) {
        SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"data inconsistency: nb. of dofs and declared nb. dofs in database");
      }
      PetscScalar *array;
      ierr = VecGetArray(V_glob,&array); CHKERRQ(ierr);
      switch (mode) {
        case INSERT_VALUES:
        for(;miit!=hi_miit;miit++) (*miit)->getFieldData() = array[(*miit)->getPetscGlobalDofIdx()];
        break;
        case ADD_VALUES:
        for(;miit!=hi_miit;miit++) (*miit)->getFieldData() += array[(*miit)->getPetscGlobalDofIdx()];
        break;
        default:
        SETERRQ(PETSC_COMM_SELF,MOFEM_NOT_IMPLEMENTED,"not implemented");
      }
      ierr = VecRestoreArray(V_glob,&array); CHKERRQ(ierr);
      ierr = VecScatterDestroy(&ctx); CHKERRQ(ierr);
      ierr = VecDestroy(&V_glob); CHKERRQ(ierr);
      break;
    }
    default:
    SETERRQ(PETSC_COMM_SELF,MOFEM_NOT_IMPLEMENTED,"not implemented");
  }
  PetscFunctionReturn(0);
}
PetscErrorCode Core::set_global_ghost_vector(
  const std::string &name,RowColData rc,Vec V,InsertMode mode,ScatterMode scatter_mode
) const {
  PetscErrorCode ierr;
  PetscFunctionBegin;
  typedef MoFEMProblem_multiIndex::index<Problem_mi_tag>::type ProblemsByName;
  const ProblemsByName &problems_set = pRoblems.get<Problem_mi_tag>();
  ProblemsByName::iterator p_miit = problems_set.find(name);
  if(p_miit==problems_set.end()) SETERRQ1(PETSC_COMM_SELF,1,"problem < %s > not found (top tip: check spelling)",name.c_str());
  ierr = set_global_ghost_vector(&*p_miit,rc,V,mode,scatter_mode); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
PetscErrorCode Core::set_other_local_ghost_vector(
  const MoFEMProblem *problem_ptr,const std::string& field_name,const std::string& cpy_field_name,RowColData rc,Vec V,InsertMode mode,ScatterMode scatter_mode,int verb) {
  PetscFunctionBegin;
  typedef NumeredDofEntity_multiIndex::index<Composite_Name_And_HasLocalIdx_mi_tag>::type DofsByNameAndLocalIdx;
  DofsByNameAndLocalIdx *dofs;
  switch (rc) {
    case ROW:
      dofs = const_cast<DofsByNameAndLocalIdx*>(&problem_ptr->numered_dofs_rows->get<Composite_Name_And_HasLocalIdx_mi_tag>());
      break;
    case COL:
      dofs = const_cast<DofsByNameAndLocalIdx*>(&problem_ptr->numered_dofs_cols->get<Composite_Name_And_HasLocalIdx_mi_tag>());
      break;
    default:
     SETERRQ(PETSC_COMM_SELF,MOFEM_NOT_IMPLEMENTED,"not implemented");
  }
  Field_multiIndex::index<FieldName_mi_tag>::type::iterator cpy_fit = fIelds.get<FieldName_mi_tag>().find(cpy_field_name);
  if(cpy_fit==fIelds.get<FieldName_mi_tag>().end()) {
    SETERRQ1(PETSC_COMM_SELF,MOFEM_NOT_FOUND,"cpy field < %s > not found, (top tip: check spelling)",cpy_field_name.c_str());
  }
  DofsByNameAndLocalIdx::iterator miit = dofs->lower_bound(boost::make_tuple(field_name,1));
  if(miit==dofs->end()) {
    SETERRQ1(PETSC_COMM_SELF,MOFEM_NOT_FOUND,"cpy field < %s > not found, (top tip: check spelling)",field_name.c_str());
  }
  DofsByNameAndLocalIdx::iterator hi_miit = dofs->upper_bound(boost::make_tuple(field_name,1));
  if((*miit)->getSpace() != (*cpy_fit)->getSpace()) {
    SETERRQ4(
      PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,
      "fields have to have same space (%s) %s != (%s) %s",
      (*miit)->getName().c_str(),
      FieldSpaceNames[(*miit)->getSpace()],
      cpy_field_name.c_str(),
      FieldSpaceNames[(*cpy_fit)->getSpace()]
    );
  }
  if((*miit)->getNbOfCoeffs() != (*cpy_fit)->getNbOfCoeffs()) {
    SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"fields have to have same rank");
  }
  switch(scatter_mode) {
    case SCATTER_REVERSE: {
      bool alpha = true;
      switch (mode) {
        case INSERT_VALUES:
        break;
        case ADD_VALUES:
        alpha = false;
        break;
        default:
        SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"not implemented");
      }
      PetscScalar *array;
      VecGetArray(V,&array);
      for(;miit!=hi_miit;miit++) {
        //if(miit->getNameRef()!=field_name) continue;
        DofEntity_multiIndex::index<Composite_Name_And_Ent_And_EntDofIdx_mi_tag>::type::iterator diiiit;
        diiiit = dofsField.get<Composite_Name_And_Ent_And_EntDofIdx_mi_tag>().find(
          boost::make_tuple(cpy_field_name,(*miit)->getEnt(),(*miit)->getEntDofIdx())
        );
        if(diiiit==dofsField.get<Composite_Name_And_Ent_And_EntDofIdx_mi_tag>().end()) {
          SETERRQ(
            PETSC_COMM_SELF,MOFEM_NOT_FOUND,
            "equivalent dof does not exist, use set_other_global_ghost_vector to create dofs entries"
          );
        }
        if(alpha) {
          (*diiiit)->getFieldData() = array[(*miit)->getPetscLocalDofIdx()];
        } else {
          (*diiiit)->getFieldData() += array[(*miit)->getPetscLocalDofIdx()];
        }
      }
      ierr = VecRestoreArray(V,&array); CHKERRQ(ierr);
    }
    break;
    case SCATTER_FORWARD: {
      for(;miit!=hi_miit;miit++) {
        //if(miit->getNameRef()!=field_name) continue;
        DofEntity_multiIndex::index<Composite_Name_And_Ent_And_EntDofIdx_mi_tag>::type::iterator diiiit;
        diiiit = dofsField.get<Composite_Name_And_Ent_And_EntDofIdx_mi_tag>().find(
          boost::make_tuple(cpy_field_name,(*miit)->getEnt(),(*miit)->getEntDofIdx())
        );
        if(diiiit==dofsField.get<Composite_Name_And_Ent_And_EntDofIdx_mi_tag>().end()) {
          SETERRQ(
            PETSC_COMM_SELF,
            MOFEM_DATA_INCONSISTENCY,
            "no data to fill the vector (top tip: you want scatter forward or scatter reverse?)"
          );
        }
        ierr = VecSetValue(V,(*miit)->getPetscGlobalDofIdx(),(*diiiit)->getFieldData(),mode); CHKERRQ(ierr);
      }
      ierr = VecAssemblyBegin(V); CHKERRQ(ierr);
      ierr = VecAssemblyEnd(V); CHKERRQ(ierr);
    }
    break;
    default:
    SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"not implemented");
  }
  PetscFunctionReturn(0);
}
PetscErrorCode Core::set_other_local_ghost_vector(
  const std::string &name,const std::string& field_name,const std::string& cpy_field_name,RowColData rc,Vec V,InsertMode mode,ScatterMode scatter_mode,int verb
) {
  PetscFunctionBegin;
  if(verb==-1) verb = verbose;
  typedef MoFEMProblem_multiIndex::index<Problem_mi_tag>::type pRoblems_by_name;
  pRoblems_by_name &pRoblems_set = pRoblems.get<Problem_mi_tag>();
  pRoblems_by_name::iterator p_miit = pRoblems_set.find(name);
  if(p_miit==pRoblems_set.end()) SETERRQ1(PETSC_COMM_SELF,1,"problem < %s > not found",name.c_str());
  ierr = set_other_local_ghost_vector(&*p_miit,field_name,cpy_field_name,rc,V,mode,scatter_mode,verb); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
PetscErrorCode Core::set_other_global_ghost_vector(
  const MoFEMProblem *problem_ptr,
  const std::string& field_name,
  const std::string& cpy_field_name,
  RowColData rc,
  Vec V,
  InsertMode mode,
  ScatterMode scatter_mode,
  int verb
) {
  MoABErrorCode rval;
  PetscErrorCode ierr;
  PetscFunctionBegin;
  if(verb==-1) verb = verbose;
  typedef NumeredDofEntityByFieldName DofsByName;
  DofsByName *dofs;
  DofIdx nb_dofs;
  switch (rc) {
    case ROW:
      nb_dofs = problem_ptr->getNbDofsRow();
      dofs = const_cast<DofsByName*>(&problem_ptr->numered_dofs_rows->get<FieldName_mi_tag>());
      break;
    case COL:
      nb_dofs = problem_ptr->getNbDofsCol();
      dofs = const_cast<DofsByName*>(&problem_ptr->numered_dofs_cols->get<FieldName_mi_tag>());
      break;
    default:
     SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"not implemented");
  }
  Field_multiIndex::index<FieldName_mi_tag>::type::iterator cpy_fit = fIelds.get<FieldName_mi_tag>().find(cpy_field_name);
  if(cpy_fit==fIelds.get<FieldName_mi_tag>().end()) {
    SETERRQ1(PETSC_COMM_SELF,MOFEM_NOT_FOUND,"cpy field < %s > not found, (top tip: check spelling)",cpy_field_name.c_str());
  }
  DofsByName::iterator miit = dofs->lower_bound(field_name);
  if(miit==dofs->end()) {
    SETERRQ1(PETSC_COMM_SELF,MOFEM_NOT_FOUND,"problem field < %s > not found, (top tip: check spelling)",field_name.c_str());
  }
  DofsByName::iterator hi_miit = dofs->upper_bound(field_name);
  if((*miit)->getSpace() != (*cpy_fit)->getSpace()) {
    SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"fields have to have same space");
  }
  if((*miit)->getNbOfCoeffs() != (*cpy_fit)->getNbOfCoeffs()) {
    SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"fields have to have same rank");
  }
  switch (scatter_mode) {
    case SCATTER_REVERSE: {
      Vec V_glob;
      VecScatter ctx;
      ierr = VecScatterCreateToAll(V,&ctx,&V_glob); CHKERRQ(ierr);
      ierr = VecScatterBegin(ctx,V,V_glob,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
      ierr = VecScatterEnd(ctx,V,V_glob,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
      int size;
      ierr = VecGetSize(V_glob,&size); CHKERRQ(ierr);
      if(size!=nb_dofs) SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"data inconsistency: nb. of dofs and declared nb. dofs in database");
      PetscScalar *array;
      VecGetArray(V_glob,&array);
      bool alpha = true;
      switch (mode) {
        case INSERT_VALUES:
        break;
        case ADD_VALUES:
        alpha = false;
        break;
        default:
        SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"not implemented");
      }
      for(;miit!=hi_miit;miit++) {
        if((*miit)->getPetscGlobalDofIdx()>=size) {
          SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"data inconsistency: nb. of dofs and declared nb. dofs in database");
        }
        DofEntity_multiIndex::index<Composite_Name_And_Ent_And_EntDofIdx_mi_tag>::type::iterator diiiit;
        diiiit = dofsField.get<Composite_Name_And_Ent_And_EntDofIdx_mi_tag>().find(boost::make_tuple(cpy_field_name,(*miit)->getEnt(),(*miit)->getEntDofIdx()));
        if(diiiit==dofsField.get<Composite_Name_And_Ent_And_EntDofIdx_mi_tag>().end()) {
          EntityHandle ent = (*miit)->getEnt();
          rval = moab.add_entities((*cpy_fit)->getMeshset(),&ent,1); CHKERRQ_MOAB(rval);
          //create field moabent
          ApproximationOrder order = (*miit)->getMaxOrder();
          std::pair<MoFEMEntity_multiIndex::iterator,bool> p_e_miit;
          try {
            boost::shared_ptr<MoFEMEntity> moabent(
              new MoFEMEntity(*cpy_fit,(*miit)->getRefEntityPtr())
            );
            p_e_miit = entsFields.insert(moabent);
          } catch (const std::exception& ex) {
            std::ostringstream ss;
            ss << "throw in method: " << ex.what() << " at line " << __LINE__ << " in file " << __FILE__ << std::endl;
            SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,ss.str().c_str());
          }
          if((*p_e_miit.first)->getMaxOrder()<order) {
            bool success = entsFields.modify(p_e_miit.first,MoFEMEntity_change_order(order));
            if(!success) SETERRQ(PETSC_COMM_SELF,MOFEM_OPERATION_UNSUCCESSFUL,"modification unsuccessful");
          }
          //create field moabdof
          DofEntityByNameAndEnt::iterator hi_diit,diit;
          diit = dofsField.get<Composite_Name_And_Ent_mi_tag>().lower_bound(boost::make_tuple(field_name,(*miit)->getEnt()));
          hi_diit = dofsField.get<Composite_Name_And_Ent_mi_tag>().upper_bound(boost::make_tuple(field_name,(*miit)->getEnt()));
          for(;diit!=hi_diit;diit++) {
            boost::shared_ptr<DofEntity> mdof =
            boost::shared_ptr<DofEntity>(
              new DofEntity(
                *(p_e_miit.first),
                (*diit)->getDofOrder(),
                (*diit)->getDofCoeffIdx(),
                (*diit)->getEntDofIdx()
              )
            );
            std::pair<DofEntity_multiIndex::iterator,bool> cpy_p_diit;
            cpy_p_diit = dofsField.insert(mdof);
            if(cpy_p_diit.second) {
              bool success = dofsField.modify(cpy_p_diit.first,DofEntity_active_change(true));
              if(!success) SETERRQ(PETSC_COMM_SELF,MOFEM_OPERATION_UNSUCCESSFUL,"modification unsuccessful");
            }
          }
          diiiit = dofsField.get<Composite_Name_And_Ent_And_EntDofIdx_mi_tag>().find(boost::make_tuple(cpy_field_name,(*miit)->getEnt(),(*miit)->getEntDofIdx()));
          if(diiiit==dofsField.get<Composite_Name_And_Ent_And_EntDofIdx_mi_tag>().end()) SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"data inconsistency");
        }
        if(alpha) (*diiiit)->getFieldData() = 0;
        (*diiiit)->getFieldData() += array[(*miit)->getPetscGlobalDofIdx()];
        if(verb > 1) {
          std::ostringstream ss;
          ss << *(*diiiit) << "set " << array[(*miit)->getPetscGlobalDofIdx()] << std::endl;
          PetscPrintf(comm,ss.str().c_str());
        }
      }
      ierr = VecRestoreArray(V_glob,&array); CHKERRQ(ierr);
      ierr = VecDestroy(&V_glob); CHKERRQ(ierr);
      ierr = VecScatterDestroy(&ctx); CHKERRQ(ierr);
    }
    break;
    case SCATTER_FORWARD: {
      for(;miit!=hi_miit;miit++) {
        DofEntity_multiIndex::index<Composite_Name_And_Ent_And_EntDofIdx_mi_tag>::type::iterator diiiit;
        diiiit = dofsField.get<Composite_Name_And_Ent_And_EntDofIdx_mi_tag>().find(boost::make_tuple(cpy_field_name,(*miit)->getEnt(),(*miit)->getEntDofIdx()));
        if(diiiit==dofsField.get<Composite_Name_And_Ent_And_EntDofIdx_mi_tag>().end()) {
          SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"no data to fill the vector (top tip: you want scatter forward or scatter reverse?)");
        }
        ierr = VecSetValue(V,(*miit)->getPetscGlobalDofIdx(),(*diiiit)->getFieldData(),mode); CHKERRQ(ierr);
      }
      ierr = VecAssemblyBegin(V); CHKERRQ(ierr);
      ierr = VecAssemblyEnd(V); CHKERRQ(ierr);
    }
    break;
    default:
    SETERRQ(PETSC_COMM_SELF,MOFEM_NOT_IMPLEMENTED,"not implemented");

  }
  PetscFunctionReturn(0);
}
PetscErrorCode Core::set_other_global_ghost_vector(
  const std::string &name,
  const std::string& field_name,
  const std::string& cpy_field_name,
  RowColData rc,
  Vec V,
  InsertMode mode,
  ScatterMode scatter_mode,
  int verb
) {
  PetscErrorCode ierr;
  PetscFunctionBegin;
  if(verb==-1) verb = verbose;
  typedef MoFEMProblem_multiIndex::index<Problem_mi_tag>::type ProblemsByName;
  ProblemsByName &problems_set = pRoblems.get<Problem_mi_tag>();
  ProblemsByName::iterator p_miit = problems_set.find(name);
  if(p_miit==problems_set.end()) SETERRQ1(PETSC_COMM_SELF,1,"problem < %s > not found",name.c_str());
  ierr = set_other_global_ghost_vector(&*p_miit,field_name,cpy_field_name,rc,V,mode,scatter_mode,verb); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}


}
