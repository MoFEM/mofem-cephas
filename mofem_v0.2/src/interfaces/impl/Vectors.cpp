/** \file Vectors.cpp
 * \brief Myltindex containes, data structures and other low-level functions
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

#include <moab/ParallelComm.hpp>

#include <petscsys.h>
#include <petscvec.h>
#include <petscmat.h>
#include <petscsnes.h>
#include <petscts.h>

#include <definitions.h>
#include <h1_hdiv_hcurl_l2.h>

#include <Common.hpp>
#include <LoopMethods.hpp>

#include <boost/ptr_container/ptr_map.hpp>
#include <Core.hpp>

#include <CoreDataStructures.hpp>

namespace MoFEM {

const static int debug = 1;

PetscErrorCode Core::VecCreateSeq(const string &name,RowColData rc,Vec *V) {
  PetscFunctionBegin;
  typedef MoFEMProblem_multiIndex::index<Problem_mi_tag>::type moFEMProblems_by_name;
  //typedef NumeredDofMoFEMEntity_multiIndex::index<PetscLocalIdx_mi_tag>::type dofs_by_local_idx;
  moFEMProblems_by_name &moFEMProblems_set = moFEMProblems.get<Problem_mi_tag>();
  moFEMProblems_by_name::iterator p_miit = moFEMProblems_set.find(name);
  if(p_miit==moFEMProblems_set.end()) SETERRQ1(PETSC_COMM_SELF,1,"no such problem %s (top tip check spelling)",name.c_str());
  DofIdx nb_local_dofs,nb_ghost_dofs;
  switch (rc) {
    case ROW:
      nb_local_dofs = p_miit->get_nb_local_dofs_row();
      nb_ghost_dofs = p_miit->get_nb_ghost_dofs_row();
      break;
    case COL:
      nb_local_dofs = p_miit->get_nb_local_dofs_col();
      nb_ghost_dofs = p_miit->get_nb_ghost_dofs_col();
      break;
    default:
     SETERRQ(PETSC_COMM_SELF,MOFEM_NOT_IMPLEMENTED,"not implemented");
  }
  ierr = ::VecCreateSeq(PETSC_COMM_SELF,nb_local_dofs+nb_ghost_dofs,V); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
PetscErrorCode Core::VecCreateGhost(const string &name,RowColData rc,Vec *V) {
  PetscFunctionBegin;
  typedef MoFEMProblem_multiIndex::index<Problem_mi_tag>::type moFEMProblems_by_name;
  typedef NumeredDofMoFEMEntity_multiIndex::index<PetscLocalIdx_mi_tag>::type dofs_by_local_idx;
  moFEMProblems_by_name &moFEMProblems_set = moFEMProblems.get<Problem_mi_tag>();
  moFEMProblems_by_name::iterator p_miit = moFEMProblems_set.find(name);
  if(p_miit==moFEMProblems_set.end()) SETERRQ1(PETSC_COMM_SELF,1,"no such problem %s (top tip check spelling)",name.c_str());
  DofIdx nb_dofs,nb_local_dofs,nb_ghost_dofs;
  dofs_by_local_idx *dofs;
  switch (rc) {
    case ROW:
      nb_dofs = p_miit->get_nb_dofs_row();
      nb_local_dofs = p_miit->get_nb_local_dofs_row();
      nb_ghost_dofs = p_miit->get_nb_ghost_dofs_row();
      dofs = const_cast<dofs_by_local_idx*>(&p_miit->numered_dofs_rows.get<PetscLocalIdx_mi_tag>());
      break;
    case COL:
      nb_dofs = p_miit->get_nb_dofs_col();
      nb_local_dofs = p_miit->get_nb_local_dofs_col();
      nb_ghost_dofs = p_miit->get_nb_ghost_dofs_col();
      dofs = const_cast<dofs_by_local_idx*>(&p_miit->numered_dofs_cols.get<PetscLocalIdx_mi_tag>());
      break;
    default:
     SETERRQ(PETSC_COMM_SELF,MOFEM_NOT_IMPLEMENTED,"not implemented");
  }
  dofs_by_local_idx::iterator miit = dofs->lower_bound(nb_local_dofs);
  dofs_by_local_idx::iterator hi_miit = dofs->upper_bound(nb_local_dofs+nb_ghost_dofs);
  int count = distance(miit,hi_miit);
  if(count != nb_ghost_dofs) SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCT,"data inconsistency");
  vector<DofIdx> ghost_idx(count);
  vector<DofIdx>::iterator vit = ghost_idx.begin();
  for(;miit!=hi_miit;miit++,vit++) *vit = miit->petsc_gloabl_dof_idx;
  ierr = ::VecCreateGhost(comm,nb_local_dofs,nb_dofs,nb_ghost_dofs,&ghost_idx[0],V); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
PetscErrorCode Core::ISCreateProblemOrder(const string &problem,RowColData rc,int min_order,int max_order,IS *is,int verb) {
  PetscFunctionBegin;
  if(verb==-1) verb = verbose;
  typedef MoFEMProblem_multiIndex::index<Problem_mi_tag>::type moFEMProblems_by_name;
  moFEMProblems_by_name &moFEMProblems_set = moFEMProblems.get<Problem_mi_tag>();
  moFEMProblems_by_name::iterator p = moFEMProblems_set.find(problem);
  if(p==moFEMProblems_set.end()) SETERRQ1(PETSC_COMM_SELF,1,"no such problem %s (top tip check spelling)",problem.c_str());
  typedef NumeredDofMoFEMEntity_multiIndex::index<Composite_Part_And_Oder_mi_tag>::type dofs_order;
  dofs_order::iterator it,hi_it;
  switch(rc) {
    case ROW:
    it = p->numered_dofs_rows.get<Composite_Part_And_Oder_mi_tag>().lower_bound(boost::make_tuple(rAnk,min_order));
    hi_it = p->numered_dofs_rows.get<Composite_Part_And_Oder_mi_tag>().upper_bound(boost::make_tuple(rAnk,max_order));
    break;
    case COL:
    it = p->numered_dofs_cols.get<Composite_Part_And_Oder_mi_tag>().lower_bound(boost::make_tuple(rAnk,min_order));
    hi_it = p->numered_dofs_cols.get<Composite_Part_And_Oder_mi_tag>().upper_bound(boost::make_tuple(rAnk,max_order));
    break;
    default:
     SETERRQ(PETSC_COMM_SELF,MOFEM_NOT_IMPLEMENTED,"not implemented");
  }

  NumeredDofMoFEMEntity_multiIndex_petsc_local_dof_view_ordered_non_unique dof_loc_idx_view;
  for(;it!=hi_it;it++) {
    pair<NumeredDofMoFEMEntity_multiIndex_petsc_local_dof_view_ordered_non_unique::iterator,bool> p;
    if(it->get_part()!=(unsigned int)rAnk) continue;
    p = dof_loc_idx_view.insert(&*it);
  }
  NumeredDofMoFEMEntity_multiIndex_petsc_local_dof_view_ordered_non_unique::iterator vit,hi_vit;
  vit = dof_loc_idx_view.begin();
  hi_vit = dof_loc_idx_view.end();

  int size = distance(vit,hi_vit);
  int *id;
  ierr = PetscMalloc(size*sizeof(int),&id); CHKERRQ(ierr);
  for(int ii = 0;vit!=hi_vit;vit++) {
    id[ii++] = (*vit)->get_petsc_gloabl_dof_idx();
  }

  ierr = ISCreateGeneral(comm,size,id,PETSC_OWN_POINTER,is); CHKERRQ(ierr);

  PetscFunctionReturn(0);
}
PetscErrorCode Core::ISCreateProblemFieldAndRank(const string &problem,RowColData rc,const string &field,int min_rank,int max_rank,IS *is,int verb) {
  PetscFunctionBegin;

  if(verb==-1) verb = verbose;
  typedef MoFEMProblem_multiIndex::index<Problem_mi_tag>::type moFEMProblems_by_name;
  moFEMProblems_by_name &moFEMProblems_set = moFEMProblems.get<Problem_mi_tag>();
  moFEMProblems_by_name::iterator p = moFEMProblems_set.find(problem);
  if(p==moFEMProblems_set.end()) SETERRQ1(PETSC_COMM_SELF,1,"no such problem %s (top tip check spelling)",problem.c_str());
  typedef NumeredDofMoFEMEntity_multiIndex::index<Composite_Name_Part_And_Rank_mi_tag>::type dofs_by_name_and_rank;
  dofs_by_name_and_rank::iterator it,hi_it;
  switch(rc) {
    case ROW:
    it = p->numered_dofs_rows.get<Composite_Name_Part_And_Rank_mi_tag>().lower_bound(boost::make_tuple(field,rAnk,min_rank));
    hi_it = p->numered_dofs_rows.get<Composite_Name_Part_And_Rank_mi_tag>().upper_bound(boost::make_tuple(field,rAnk,max_rank));
    break;
    case COL:
    it = p->numered_dofs_cols.get<Composite_Name_Part_And_Rank_mi_tag>().lower_bound(boost::make_tuple(field,rAnk,min_rank));
    hi_it = p->numered_dofs_cols.get<Composite_Name_Part_And_Rank_mi_tag>().upper_bound(boost::make_tuple(field,rAnk,max_rank));
    break;
    default:
     SETERRQ(PETSC_COMM_SELF,MOFEM_NOT_IMPLEMENTED,"not implemented");
  }


  NumeredDofMoFEMEntity_multiIndex_petsc_local_dof_view_ordered_non_unique dof_loc_idx_view;
  for(;it!=hi_it;it++) {
    pair<NumeredDofMoFEMEntity_multiIndex_petsc_local_dof_view_ordered_non_unique::iterator,bool> p;
    if(it->get_part()!=(unsigned int)rAnk) continue;
    if(it->get_name_ref() != field) continue;
    p = dof_loc_idx_view.insert(&*it);
  }
  NumeredDofMoFEMEntity_multiIndex_petsc_local_dof_view_ordered_non_unique::iterator vit,hi_vit;
  vit = dof_loc_idx_view.begin();
  hi_vit = dof_loc_idx_view.end();

  int size = distance(vit,hi_vit);
  int *id;
  ierr = PetscMalloc(size*sizeof(int),&id); CHKERRQ(ierr);
  for(int ii = 0;vit!=hi_vit;vit++) {
    id[ii++] = (*vit)->get_petsc_gloabl_dof_idx();
  }

  ierr = ISCreateGeneral(comm,size,id,PETSC_OWN_POINTER,is); CHKERRQ(ierr);

  PetscFunctionReturn(0);
}
PetscErrorCode Core::ISCreateFromProblemFieldToOtherProblemField(
    const string &x_problem,const string &x_field_name,RowColData x_rc,
    const string &y_problem,const string &y_field_name,RowColData y_rc,
    vector<int> &idx,vector<int> &idy,int verb) {
  PetscFunctionBegin;
  if(verb==-1) verb = verbose;
  typedef MoFEMProblem_multiIndex::index<Problem_mi_tag>::type moFEMProblems_by_name;
  moFEMProblems_by_name &moFEMProblems_set = moFEMProblems.get<Problem_mi_tag>();
  moFEMProblems_by_name::iterator p_x = moFEMProblems_set.find(x_problem);
  if(p_x==moFEMProblems_set.end()) SETERRQ1(PETSC_COMM_SELF,1,"no such problem %s (top tip check spelling)",x_problem.c_str());
  moFEMProblems_by_name::iterator p_y = moFEMProblems_set.find(y_problem);
  if(p_y==moFEMProblems_set.end()) SETERRQ1(PETSC_COMM_SELF,1,"no such problem %s (top tip check spelling)",y_problem.c_str());
  typedef NumeredDofMoFEMEntity_multiIndex::index<PetscLocalIdx_mi_tag>::type dofs_by_glob_idx;
  dofs_by_glob_idx::iterator y_dit,hi_y_dit;
  switch (y_rc) {
    case ROW:
      y_dit = p_y->numered_dofs_rows.get<PetscLocalIdx_mi_tag>().lower_bound(0);
      hi_y_dit = p_y->numered_dofs_rows.get<PetscLocalIdx_mi_tag>().upper_bound(p_y->get_nb_local_dofs_row()-1);
      break;
    case COL:
      y_dit = p_y->numered_dofs_cols.get<PetscLocalIdx_mi_tag>().lower_bound(0);
      hi_y_dit = p_y->numered_dofs_cols.get<PetscLocalIdx_mi_tag>().upper_bound(p_y->get_nb_local_dofs_col()-1);
      break;
    default:
     SETERRQ(PETSC_COMM_SELF,MOFEM_NOT_IMPLEMENTED,"not implemented");
  }
  typedef NumeredDofMoFEMEntity_multiIndex::index<Composite_Name_And_Ent_And_EndDofIdx_mi_tag>::type dofs_by_name_ent_dof;
  const dofs_by_name_ent_dof* x_numered_dofs_by_ent_name_dof;
  switch (x_rc) {
    case ROW:
      x_numered_dofs_by_ent_name_dof = &(p_x->numered_dofs_rows.get<Composite_Name_And_Ent_And_EndDofIdx_mi_tag>());
      break;
    case COL:
      x_numered_dofs_by_ent_name_dof = &(p_x->numered_dofs_cols.get<Composite_Name_And_Ent_And_EndDofIdx_mi_tag>());
      break;
    default:
     SETERRQ(PETSC_COMM_SELF,MOFEM_NOT_IMPLEMENTED,"not implemented");
  }
  for(;y_dit!=hi_y_dit;y_dit++) {
    if(y_dit->get_part()!=(unsigned int)rAnk) SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCT,"data inconsistency");
    if(y_dit->get_name()!=y_field_name) continue;
    dofs_by_name_ent_dof::iterator x_dit;
    x_dit = x_numered_dofs_by_ent_name_dof->find(boost::make_tuple(x_field_name,y_dit->get_ent(),y_dit->get_EntDofIdx()));
    if(x_dit==x_numered_dofs_by_ent_name_dof->end()) continue;
    idx.push_back(x_dit->get_petsc_gloabl_dof_idx());
    idy.push_back(y_dit->get_petsc_gloabl_dof_idx());
  }
  PetscFunctionReturn(0);
}
PetscErrorCode Core::ISCreateFromProblemFieldToOtherProblemField(
  const string &x_problem,const string &x_field_name,RowColData x_rc,
  const string &y_problem,const string &y_field_name,RowColData y_rc,
  IS *ix,IS *iy,int verb) {
  PetscFunctionBegin;
  vector<int> idx(0),idy(0);
  ierr = ISCreateFromProblemFieldToOtherProblemField(
    x_problem,x_field_name,x_rc,y_problem,y_field_name,y_rc,
    idx,idy,verb); CHKERRQ(ierr);
  ierr = ISCreateGeneral(comm,idx.size(),&idx[0],PETSC_COPY_VALUES,ix); CHKERRQ(ierr);
  ierr = ISCreateGeneral(comm,idy.size(),&idy[0],PETSC_COPY_VALUES,iy); CHKERRQ(ierr);
  if(verb>2) {
    ISView(*ix,PETSC_VIEWER_STDOUT_WORLD);
    ISView(*iy,PETSC_VIEWER_STDOUT_WORLD);
  }
  PetscFunctionReturn(0);
}
PetscErrorCode Core::VecScatterCreate(Vec xin,const string &x_problem,const string &x_field_name,RowColData x_rc,
  Vec yin,const string &y_problem,const string &y_field_name,RowColData y_rc,VecScatter *newctx,int verb) {
  PetscFunctionBegin;
  if(verb==-1) verb = verbose;
  vector<int> idx(0),idy(0);
  ierr = ISCreateFromProblemFieldToOtherProblemField(
    x_problem,x_field_name,x_rc,y_problem,y_field_name,y_rc,
    idx,idy,verb); CHKERRQ(ierr);
  IS ix,iy;
  ierr = ISCreateGeneral(comm,idx.size(),&idx[0],PETSC_USE_POINTER,&ix); CHKERRQ(ierr);
  ierr = ISCreateGeneral(comm,idy.size(),&idy[0],PETSC_USE_POINTER,&iy); CHKERRQ(ierr);
  ierr = ::VecScatterCreate(xin,ix,yin,iy,newctx); CHKERRQ(ierr);
  ierr = ISDestroy(&ix); CHKERRQ(ierr);
  ierr = ISDestroy(&iy); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
PetscErrorCode Core::ISCreateFromProblemToOtherProblem(
  const string &x_problem,RowColData x_rc,const string &y_problem,RowColData y_rc,vector<int> &idx,vector<int> &idy,int verb) {
  PetscFunctionBegin;
  if(verb==-1) verb = verbose;
  typedef MoFEMProblem_multiIndex::index<Problem_mi_tag>::type moFEMProblems_by_name;
  moFEMProblems_by_name &moFEMProblems_set = moFEMProblems.get<Problem_mi_tag>();
  moFEMProblems_by_name::iterator p_x = moFEMProblems_set.find(x_problem);
  if(p_x==moFEMProblems_set.end()) SETERRQ1(PETSC_COMM_SELF,1,"no such problem %s (top tip check spelling)",x_problem.c_str());
  moFEMProblems_by_name::iterator p_y = moFEMProblems_set.find(y_problem);
  if(p_y==moFEMProblems_set.end()) SETERRQ1(PETSC_COMM_SELF,1,"no such problem %s (top tip check spelling)",y_problem.c_str());
  typedef NumeredDofMoFEMEntity_multiIndex::index<PetscLocalIdx_mi_tag>::type dofs_by_glob_idx;
  dofs_by_glob_idx::iterator y_dit,hi_y_dit;
  switch (y_rc) {
    case ROW:
      y_dit = p_y->numered_dofs_rows.get<PetscLocalIdx_mi_tag>().lower_bound(0);
      hi_y_dit = p_y->numered_dofs_rows.get<PetscLocalIdx_mi_tag>().lower_bound(p_y->get_nb_local_dofs_row()); // should be lower
      break;
    case COL:
      y_dit = p_y->numered_dofs_cols.get<PetscLocalIdx_mi_tag>().lower_bound(0);
      hi_y_dit = p_y->numered_dofs_cols.get<PetscLocalIdx_mi_tag>().lower_bound(p_y->get_nb_local_dofs_col()); // should be lower
      break;
    default:
     SETERRQ(PETSC_COMM_SELF,MOFEM_NOT_IMPLEMENTED,"not implemented");
  }
  typedef NumeredDofMoFEMEntity_multiIndex::index<Unique_mi_tag>::type dofs_by_uid;
  const dofs_by_uid* x_numered_dofs_by_uid;
  switch (x_rc) {
    case ROW:
      x_numered_dofs_by_uid = &(p_x->numered_dofs_rows.get<Unique_mi_tag>());
      break;
    case COL:
      x_numered_dofs_by_uid = &(p_x->numered_dofs_cols.get<Unique_mi_tag>());
      break;
    default:
     SETERRQ(PETSC_COMM_SELF,MOFEM_NOT_IMPLEMENTED,"not implemented");
  }
  for(;y_dit!=hi_y_dit;y_dit++) {
    if(y_dit->get_part()!=(unsigned int)rAnk) SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCT,"data inconsistency");
    dofs_by_uid::iterator x_dit;
    x_dit = x_numered_dofs_by_uid->find(y_dit->get_global_unique_id());
    if(x_dit==x_numered_dofs_by_uid->end()) continue;
    idx.push_back(x_dit->get_petsc_gloabl_dof_idx());
    idy.push_back(y_dit->get_petsc_gloabl_dof_idx());
  }
  PetscFunctionReturn(0);
}
PetscErrorCode Core::ISCreateFromProblemToOtherProblem(
  const string &x_problem,RowColData x_rc,const string &y_problem,RowColData y_rc,IS *ix,IS *iy,int verb) {
  PetscFunctionBegin;
  vector<int> idx(0),idy(0);
  ierr = ISCreateFromProblemToOtherProblem(x_problem,x_rc,y_problem,y_rc,idx,idy,verb); CHKERRQ(ierr);
  ierr = ISCreateGeneral(comm,idx.size(),&idx[0],PETSC_COPY_VALUES,ix); CHKERRQ(ierr);
  ierr = ISCreateGeneral(comm,idy.size(),&idy[0],PETSC_COPY_VALUES,iy); CHKERRQ(ierr);
  if(verb>2) {
    ISView(*ix,PETSC_VIEWER_STDOUT_WORLD);
    ISView(*iy,PETSC_VIEWER_STDOUT_WORLD);
  }
  PetscFunctionReturn(0);
}
PetscErrorCode Core::VecScatterCreate(Vec xin,const string &x_problem,RowColData x_rc,Vec yin,const string &y_problem,RowColData y_rc,VecScatter *newctx,int verb) {
  PetscFunctionBegin;
  if(verb==-1) verb = verbose;
  vector<int> idx(0),idy(0);
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
PetscErrorCode Core::set_local_ghost_vector(const MoFEMProblem *problem_ptr,RowColData rc,Vec V,InsertMode mode,ScatterMode scatter_mode) {
  PetscFunctionBegin;
  typedef NumeredDofMoFEMEntity_multiIndex::index<PetscLocalIdx_mi_tag>::type dofs_by_local_idx;
  dofs_by_local_idx *dofs;
  DofIdx nb_local_dofs,nb_ghost_dofs;
  switch (rc) {
    case ROW:
      nb_local_dofs = problem_ptr->get_nb_local_dofs_row();
      nb_ghost_dofs = problem_ptr->get_nb_ghost_dofs_row();
      dofs = const_cast<dofs_by_local_idx*>(&problem_ptr->numered_dofs_rows.get<PetscLocalIdx_mi_tag>());
      break;
    case COL:
      nb_local_dofs = problem_ptr->get_nb_local_dofs_col();
      nb_ghost_dofs = problem_ptr->get_nb_ghost_dofs_col();
      dofs = const_cast<dofs_by_local_idx*>(&problem_ptr->numered_dofs_cols.get<PetscLocalIdx_mi_tag>());
      break;
    default:
     SETERRQ(PETSC_COMM_SELF,MOFEM_NOT_IMPLEMENTED,"not implemented");
  }
  Vec Vlocal;
  ierr = VecGhostGetLocalForm(V,&Vlocal); CHKERRQ(ierr);
  int size;
  ierr = VecGetLocalSize(V,&size); CHKERRQ(ierr);
  if(size!=nb_local_dofs) SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCT,"data inconsistency: check ghost vector, problem with nb. of local nodes");
  ierr = VecGetLocalSize(Vlocal,&size); CHKERRQ(ierr);
  if(size!=nb_local_dofs+nb_ghost_dofs) SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCT,"data inconsistency: check ghost vector, problem with nb. of ghost nodes");
  dofs_by_local_idx::iterator miit = dofs->lower_bound(0);
  dofs_by_local_idx::iterator hi_miit = dofs->upper_bound(nb_local_dofs+nb_ghost_dofs);
  PetscScalar *array;
  VecGetArray(Vlocal,&array);
  DofIdx ii = 0;
  switch (scatter_mode) {
    case SCATTER_FORWARD:
      switch (mode) {
	case INSERT_VALUES:
	  for(;miit!=hi_miit;miit++,ii++) array[ii] = miit->get_FieldData();
	  break;
	case ADD_VALUES:
	  for(;miit!=hi_miit;miit++,ii++) array[ii] += miit->get_FieldData();
	  break;
	default:
	  SETERRQ(PETSC_COMM_SELF,MOFEM_NOT_IMPLEMENTED,"not implemented");
      }
    break;
    case SCATTER_REVERSE:
      switch (mode) {
	case INSERT_VALUES:
	  for(;miit!=hi_miit;miit++,ii++) {
	    //cerr << *miit << endl;
	    //cerr << array[ii] << endl;
	    miit->get_FieldData() = array[ii];
	  }
	  break;
	case ADD_VALUES:
	  for(;miit!=hi_miit;miit++,ii++) miit->get_FieldData() += array[ii];
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
PetscErrorCode Core::set_local_ghost_vector(const string &name,RowColData rc,Vec V,InsertMode mode,ScatterMode scatter_mode) {
  PetscFunctionBegin;
  typedef MoFEMProblem_multiIndex::index<Problem_mi_tag>::type moFEMProblems_by_name;
  moFEMProblems_by_name &moFEMProblems_set = moFEMProblems.get<Problem_mi_tag>();
  moFEMProblems_by_name::iterator p_miit = moFEMProblems_set.find(name);
  if(p_miit==moFEMProblems_set.end()) SETERRQ1(PETSC_COMM_SELF,1,"problem < %s > not found (top tip: check spelling)",name.c_str());
  ierr = set_local_ghost_vector(&*p_miit,rc,V,mode,scatter_mode); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
PetscErrorCode Core::set_global_ghost_vector(const MoFEMProblem *problem_ptr,RowColData rc,Vec V,InsertMode mode,ScatterMode scatter_mode) {
  PetscFunctionBegin;
  typedef NumeredDofMoFEMEntity_multiIndex::index<PetscGlobalIdx_mi_tag>::type dofs_by_global_idx;
  dofs_by_global_idx *dofs;
  DofIdx nb_dofs;
  switch (rc) {
    case ROW:
      nb_dofs = problem_ptr->get_nb_dofs_row();
      dofs = const_cast<dofs_by_global_idx*>(&problem_ptr->numered_dofs_rows.get<PetscGlobalIdx_mi_tag>());
      break;
    case COL:
      nb_dofs = problem_ptr->get_nb_dofs_col();
      dofs = const_cast<dofs_by_global_idx*>(&problem_ptr->numered_dofs_cols.get<PetscGlobalIdx_mi_tag>());
      break;
    default:
     SETERRQ(PETSC_COMM_SELF,MOFEM_NOT_IMPLEMENTED,"not implemented");
  }
  dofs_by_global_idx::iterator miit = dofs->lower_bound(0);
  dofs_by_global_idx::iterator hi_miit = dofs->upper_bound(nb_dofs);
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
        SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCT,"data inconsistency: nb. of dofs and declared nb. dofs in database");
      }
      PetscScalar *array;
      ierr = VecGetArray(V_glob,&array); CHKERRQ(ierr);
      switch (mode) {
        case INSERT_VALUES:
        for(;miit!=hi_miit;miit++) miit->get_FieldData() = array[miit->get_petsc_gloabl_dof_idx()];
        break;
        case ADD_VALUES:
        for(;miit!=hi_miit;miit++) miit->get_FieldData() += array[miit->get_petsc_gloabl_dof_idx()];
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
PetscErrorCode Core::set_global_ghost_vector(const string &name,RowColData rc,Vec V,InsertMode mode,ScatterMode scatter_mode) {
  PetscFunctionBegin;
  typedef MoFEMProblem_multiIndex::index<Problem_mi_tag>::type moFEMProblems_by_name;
  moFEMProblems_by_name &moFEMProblems_set = moFEMProblems.get<Problem_mi_tag>();
  moFEMProblems_by_name::iterator p_miit = moFEMProblems_set.find(name);
  if(p_miit==moFEMProblems_set.end()) SETERRQ1(PETSC_COMM_SELF,1,"problem < %s > not found (top tip: check spelling)",name.c_str());
  ierr = set_global_ghost_vector(&*p_miit,rc,V,mode,scatter_mode); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
PetscErrorCode Core::set_other_local_ghost_vector(
  const MoFEMProblem *problem_ptr,const string& field_name,const string& cpy_field_name,RowColData rc,Vec V,InsertMode mode,ScatterMode scatter_mode,int verb) {
  PetscFunctionBegin;
  typedef NumeredDofMoFEMEntity_multiIndex::index<Composite_Name_And_HasLocalIdx_mi_tag>::type DofsByNameAndLocalIdx;
  DofsByNameAndLocalIdx *dofs;
  switch (rc) {
    case ROW:
      dofs = const_cast<DofsByNameAndLocalIdx*>(&problem_ptr->numered_dofs_rows.get<Composite_Name_And_HasLocalIdx_mi_tag>());
      break;
    case COL:
      dofs = const_cast<DofsByNameAndLocalIdx*>(&problem_ptr->numered_dofs_cols.get<Composite_Name_And_HasLocalIdx_mi_tag>());
      break;
    default:
     SETERRQ(PETSC_COMM_SELF,MOFEM_NOT_IMPLEMENTED,"not implemented");
  }
  MoFEMField_multiIndex::index<FieldName_mi_tag>::type::iterator cpy_fit = moabFields.get<FieldName_mi_tag>().find(cpy_field_name);
  if(cpy_fit==moabFields.get<FieldName_mi_tag>().end()) {
    SETERRQ1(PETSC_COMM_SELF,MOFEM_NOT_FOUND,"cpy field < %s > not found, (top tip: check spelling)",cpy_field_name.c_str());
  }
  DofsByNameAndLocalIdx::iterator miit = dofs->lower_bound(boost::make_tuple(field_name,1));
  if(miit==dofs->end()) {
    SETERRQ1(PETSC_COMM_SELF,MOFEM_NOT_FOUND,"cpy field < %s > not found, (top tip: check spelling)",field_name.c_str());
  }
  DofsByNameAndLocalIdx::iterator hi_miit = dofs->upper_bound(boost::make_tuple(field_name,1));
  if(miit->get_space() != cpy_fit->get_space()) {
    SETERRQ4(
      PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCT,
      "fields have to have same space (%s) %s != (%s) %s",
      miit->get_name().c_str(),
      FieldSpaceNames[miit->get_space()],
      cpy_field_name.c_str(),
      FieldSpaceNames[cpy_fit->get_space()]
    );
  }
  if(miit->get_max_rank() != cpy_fit->get_max_rank()) {
    SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCT,"fields have to have same rank");
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
        SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCT,"not implemented");
      }
      PetscScalar *array;
      VecGetArray(V,&array);
      for(;miit!=hi_miit;miit++) {
        //if(miit->get_name_ref()!=field_name) continue;
        DofMoFEMEntity_multiIndex::index<Composite_Name_And_Ent_And_EndDofIdx_mi_tag>::type::iterator diiiit;
        diiiit = dofsMoabField.get<Composite_Name_And_Ent_And_EndDofIdx_mi_tag>().find(
          boost::make_tuple(cpy_field_name,miit->get_ent(),miit->get_EntDofIdx())
        );
        if(diiiit==dofsMoabField.get<Composite_Name_And_Ent_And_EndDofIdx_mi_tag>().end()) {
          SETERRQ(
            PETSC_COMM_SELF,MOFEM_NOT_FOUND,
            "equivalent dof does not exist, use set_other_global_ghost_vector to create dofs entries"
          );
        }
        if(alpha) {
          diiiit->get_FieldData() = array[miit->get_petsc_local_dof_idx()];
        } else {
          diiiit->get_FieldData() += array[miit->get_petsc_local_dof_idx()];
        }
      }
      ierr = VecRestoreArray(V,&array); CHKERRQ(ierr);
    }
    break;
    case SCATTER_FORWARD: {
      for(;miit!=hi_miit;miit++) {
        //if(miit->get_name_ref()!=field_name) continue;
        DofMoFEMEntity_multiIndex::index<Composite_Name_And_Ent_And_EndDofIdx_mi_tag>::type::iterator diiiit;
        diiiit = dofsMoabField.get<Composite_Name_And_Ent_And_EndDofIdx_mi_tag>().find(
          boost::make_tuple(cpy_field_name,miit->get_ent(),miit->get_EntDofIdx())
        );
        if(diiiit==dofsMoabField.get<Composite_Name_And_Ent_And_EndDofIdx_mi_tag>().end()) {
          SETERRQ(
            PETSC_COMM_SELF,
            MOFEM_DATA_INCONSISTENCT,
            "no data to fill the vector (top tip: you want scatter forward or scatter reverse?)"
          );
        }
        ierr = VecSetValue(V,miit->get_petsc_gloabl_dof_idx(),diiiit->get_FieldData(),mode); CHKERRQ(ierr);
      }
      ierr = VecAssemblyBegin(V); CHKERRQ(ierr);
      ierr = VecAssemblyEnd(V); CHKERRQ(ierr);
    }
    break;
    default:
    SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCT,"not implemented");
  }
  PetscFunctionReturn(0);
}
PetscErrorCode Core::set_other_local_ghost_vector(
  const string &name,const string& field_name,const string& cpy_field_name,RowColData rc,Vec V,InsertMode mode,ScatterMode scatter_mode,int verb
) {
  PetscFunctionBegin;
  if(verb==-1) verb = verbose;
  typedef MoFEMProblem_multiIndex::index<Problem_mi_tag>::type moFEMProblems_by_name;
  moFEMProblems_by_name &moFEMProblems_set = moFEMProblems.get<Problem_mi_tag>();
  moFEMProblems_by_name::iterator p_miit = moFEMProblems_set.find(name);
  if(p_miit==moFEMProblems_set.end()) SETERRQ1(PETSC_COMM_SELF,1,"problem < %s > not found",name.c_str());
  ierr = set_other_local_ghost_vector(&*p_miit,field_name,cpy_field_name,rc,V,mode,scatter_mode,verb); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
PetscErrorCode Core::set_other_global_ghost_vector(
  const MoFEMProblem *problem_ptr,const string& field_name,const string& cpy_field_name,RowColData rc,Vec V,InsertMode mode,ScatterMode scatter_mode,
  int verb) {
  PetscFunctionBegin;
  if(verb==-1) verb = verbose;
  typedef NumeredDofMoFEMEntity_multiIndex::index<FieldName_mi_tag>::type dofs_by_name;
  dofs_by_name *dofs;
  DofIdx nb_dofs;
  switch (rc) {
    case ROW:
      nb_dofs = problem_ptr->get_nb_dofs_row();
      dofs = const_cast<dofs_by_name*>(&problem_ptr->numered_dofs_rows.get<FieldName_mi_tag>());
      break;
    case COL:
      nb_dofs = problem_ptr->get_nb_dofs_col();
      dofs = const_cast<dofs_by_name*>(&problem_ptr->numered_dofs_cols.get<FieldName_mi_tag>());
      break;
    default:
     SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCT,"not implemented");
  }
  MoFEMField_multiIndex::index<FieldName_mi_tag>::type::iterator cpy_fit = moabFields.get<FieldName_mi_tag>().find(cpy_field_name);
  if(cpy_fit==moabFields.get<FieldName_mi_tag>().end()) {
    SETERRQ1(PETSC_COMM_SELF,MOFEM_NOT_FOUND,"cpy field < %s > not found, (top tip: check spelling)",cpy_field_name.c_str());
  }
  dofs_by_name::iterator miit = dofs->lower_bound(field_name);
  if(miit==dofs->end()) {
    SETERRQ1(PETSC_COMM_SELF,MOFEM_NOT_FOUND,"problem field < %s > not found, (top tip: check spelling)",field_name.c_str());
  }
  dofs_by_name::iterator hi_miit = dofs->upper_bound(field_name);
  if(miit->get_space() != cpy_fit->get_space()) {
    SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCT,"fields have to have same space");
  }
  if(miit->get_max_rank() != cpy_fit->get_max_rank()) {
    SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCT,"fields have to have same rank");
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
      if(size!=nb_dofs) SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCT,"data inconsistency: nb. of dofs and declared nb. dofs in database");
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
	  SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCT,"not implemented");
      }
      for(;miit!=hi_miit;miit++) {
        if(miit->get_petsc_gloabl_dof_idx()>=size) {
          SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCT,"data inconsistency: nb. of dofs and declared nb. dofs in database");
        }
        DofMoFEMEntity_multiIndex::index<Composite_Name_And_Ent_And_EndDofIdx_mi_tag>::type::iterator diiiit;
        diiiit = dofsMoabField.get<Composite_Name_And_Ent_And_EndDofIdx_mi_tag>().find(boost::make_tuple(cpy_field_name,miit->get_ent(),miit->get_EntDofIdx()));
        if(diiiit==dofsMoabField.get<Composite_Name_And_Ent_And_EndDofIdx_mi_tag>().end()) {
          EntityHandle ent = miit->get_ent();
          rval = moab.add_entities(cpy_fit->get_meshset(),&ent,1); CHKERR_PETSC(rval);
          //create field moabent
          ApproximationOrder order = miit->get_max_order();
          pair<MoFEMEntity_multiIndex::iterator,bool> p_e_miit;
          try {
            MoFEMEntity moabent(moab,cpy_fit->get_MoFEMField_ptr(),miit->get_RefMoFEMEntity_ptr());
            p_e_miit = entsMoabField.insert(moabent);
          } catch (const std::exception& ex) {
            ostringstream ss;
            ss << "throw in method: " << ex.what() << " at line " << __LINE__ << " in file " << __FILE__ << endl;
            SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCT,ss.str().c_str());
          }
          if(p_e_miit.first->get_max_order()<order) {
            bool success = entsMoabField.modify(p_e_miit.first,MoFEMEntity_change_order(moab,order));
            if(!success) SETERRQ(PETSC_COMM_SELF,MOFEM_OPERATION_UNSUCCESSFUL,"modification unsuccessful");
          }
          //create field moabdof
          DofMoFEMEntity_multiIndex::index<Composite_Name_And_Ent_mi_tag>::type::iterator hi_diit,diit;
          diit = dofsMoabField.get<Composite_Name_And_Ent_mi_tag>().lower_bound(boost::make_tuple(field_name,miit->get_ent()));
          hi_diit = dofsMoabField.get<Composite_Name_And_Ent_mi_tag>().upper_bound(boost::make_tuple(field_name,miit->get_ent()));
          for(;diit!=hi_diit;diit++) {
            DofMoFEMEntity mdof(&*(p_e_miit.first),diit->get_dof_order(),diit->get_dof_rank(),diit->get_EntDofIdx());
            pair<DofMoFEMEntity_multiIndex::iterator,bool> cpy_p_diit;
            cpy_p_diit = dofsMoabField.insert(mdof);
            if(cpy_p_diit.second) {
              bool success = dofsMoabField.modify(cpy_p_diit.first,DofMoFEMEntity_active_change(true));
              if(!success) SETERRQ(PETSC_COMM_SELF,MOFEM_OPERATION_UNSUCCESSFUL,"modification unsuccessful");
            }
          }
          diiiit = dofsMoabField.get<Composite_Name_And_Ent_And_EndDofIdx_mi_tag>().find(boost::make_tuple(cpy_field_name,miit->get_ent(),miit->get_EntDofIdx()));
          if(diiiit==dofsMoabField.get<Composite_Name_And_Ent_And_EndDofIdx_mi_tag>().end()) SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCT,"data inconsistency");
        }
        if(alpha) diiiit->get_FieldData() = 0;
        diiiit->get_FieldData() += array[miit->get_petsc_gloabl_dof_idx()];
        if(verb > 1) {
          ostringstream ss;
          ss << *diiiit << "set " << array[miit->get_petsc_gloabl_dof_idx()] << endl;
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
        DofMoFEMEntity_multiIndex::index<Composite_Name_And_Ent_And_EndDofIdx_mi_tag>::type::iterator diiiit;
        diiiit = dofsMoabField.get<Composite_Name_And_Ent_And_EndDofIdx_mi_tag>().find(boost::make_tuple(cpy_field_name,miit->get_ent(),miit->get_EntDofIdx()));
        if(diiiit==dofsMoabField.get<Composite_Name_And_Ent_And_EndDofIdx_mi_tag>().end()) {
          SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCT,"no data to fill the vector (top tip: you want scatter forward or scatter reverse?)");
        }
        ierr = VecSetValue(V,miit->get_petsc_gloabl_dof_idx(),diiiit->get_FieldData(),mode); CHKERRQ(ierr);
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
  const string &name,const string& field_name,const string& cpy_field_name,RowColData rc,Vec V,InsertMode mode,ScatterMode scatter_mode,
  int verb) {
  PetscFunctionBegin;
  if(verb==-1) verb = verbose;
  typedef MoFEMProblem_multiIndex::index<Problem_mi_tag>::type moFEMProblems_by_name;
  moFEMProblems_by_name &moFEMProblems_set = moFEMProblems.get<Problem_mi_tag>();
  moFEMProblems_by_name::iterator p_miit = moFEMProblems_set.find(name);
  if(p_miit==moFEMProblems_set.end()) SETERRQ1(PETSC_COMM_SELF,1,"problem < %s > not found",name.c_str());
  ierr = set_other_global_ghost_vector(&*p_miit,field_name,cpy_field_name,rc,V,mode,scatter_mode,verb); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}


}
