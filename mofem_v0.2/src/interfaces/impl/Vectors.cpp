/** \file Vectors.cpp
 * \brief Myltindex containes, data structures and other low-level functions 
 * 
 * Copyright (C) 2013, Lukasz Kaczmarczyk (likask AT wp.pl) <br>
 *
 * The MoFEM package is copyrighted by Lukasz Kaczmarczyk. 
 * It can be freely used for educational and research purposes 
 * by other institutions. If you use this softwre pleas cite my work. 
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
  typedef NumeredDofMoFEMEntity_multiIndex::index<PetscLocalIdx_mi_tag>::type dofs_by_local_idx;
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
  if(count != nb_ghost_dofs) SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INSONSISTENCY,"data inconsistency");
  vector<DofIdx> ghost_idx(count);
  vector<DofIdx>::iterator vit = ghost_idx.begin();
  for(;miit!=hi_miit;miit++,vit++) *vit = miit->petsc_gloabl_dof_idx;
  ierr = ::VecCreateGhost(PETSC_COMM_WORLD,nb_local_dofs,nb_dofs,nb_ghost_dofs,&ghost_idx[0],V); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
PetscErrorCode Core::VecScatterCreate(Vec xin,string &x_problem,RowColData x_rc,Vec yin,string &y_problem,RowColData y_rc,VecScatter *newctx,int verb) {
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
  vector<int> idx(0),idy(0);
  ParallelComm* pcomm = ParallelComm::get_pcomm(&moab,MYPCOMM_INDEX);
  for(;y_dit!=hi_y_dit;y_dit++) {
    if(y_dit->get_part()!=pcomm->rank()) SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INSONSISTENCY,"data inconsistency");
    dofs_by_uid::iterator x_dit;
    x_dit = x_numered_dofs_by_uid->find(y_dit->get_global_unique_id());
    if(x_dit==x_numered_dofs_by_uid->end()) continue;
    idx.push_back(x_dit->get_petsc_gloabl_dof_idx());
    idy.push_back(y_dit->get_petsc_gloabl_dof_idx());
  }
  IS ix,iy;
  ierr = ISCreateGeneral(PETSC_COMM_WORLD,idx.size(),&idx[0],PETSC_USE_POINTER,&ix); CHKERRQ(ierr);
  ierr = ISCreateGeneral(PETSC_COMM_WORLD,idy.size(),&idy[0],PETSC_USE_POINTER,&iy); CHKERRQ(ierr);
  if(verb>3) {
    ISView(ix,PETSC_VIEWER_STDOUT_WORLD);
    ISView(iy,PETSC_VIEWER_STDOUT_WORLD);
  }
  ierr = ::VecScatterCreate(xin,ix,yin,iy,newctx); CHKERRQ(ierr);
  ierr = ISDestroy(&ix); CHKERRQ(ierr);
  ierr = ISDestroy(&iy); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
PetscErrorCode Core::set_local_VecCreateGhost(const MoFEMProblem *problem_ptr,RowColData rc,Vec V,InsertMode mode,ScatterMode scatter_mode) {
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
  if(size!=nb_local_dofs) SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INSONSISTENCY,"data inconsistency: check ghost vector, problem with nb. of local nodes");
  ierr = VecGetLocalSize(Vlocal,&size); CHKERRQ(ierr);
  if(size!=nb_local_dofs+nb_ghost_dofs) SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INSONSISTENCY,"data inconsistency: check ghost vector, problem with nb. of ghost nodes");
  dofs_by_local_idx::iterator miit = dofs->lower_bound(0);
  dofs_by_local_idx::iterator hi_miit = dofs->upper_bound(nb_local_dofs+nb_ghost_dofs);
  /*cerr << "AAAAAAAAAAAAA " << distance(miit,hi_miit) << " " << problem_ptr->numered_dofs_rows.size() << " " << pcomm->rank() << endl;
  {
    NumeredDofMoFEMEntity_multiIndex::iterator it = problem_ptr->numered_dofs_rows.begin();
    NumeredDofMoFEMEntity_multiIndex::iterator hi_it = problem_ptr->numered_dofs_rows.end();
    for(;it!=problem_ptr->numered_dofs_rows.end();it++) {
      cerr << "BBB " << *it << " " << problem_ptr->numered_dofs_rows.size() << " " << pcomm->rank() << endl;
    }
  }*/
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
PetscErrorCode Core::set_local_VecCreateGhost(const string &name,RowColData rc,Vec V,InsertMode mode,ScatterMode scatter_mode) {
  PetscFunctionBegin;
  //ParallelComm* pcomm = ParallelComm::get_pcomm(&moab,MYPCOMM_INDEX);
  typedef MoFEMProblem_multiIndex::index<Problem_mi_tag>::type moFEMProblems_by_name;
  moFEMProblems_by_name &moFEMProblems_set = moFEMProblems.get<Problem_mi_tag>();
  moFEMProblems_by_name::iterator p_miit = moFEMProblems_set.find(name);
  if(p_miit==moFEMProblems_set.end()) SETERRQ1(PETSC_COMM_SELF,1,"problem < %s > not found (top tip: check spelling)",name.c_str()); 
  ierr = set_local_VecCreateGhost(&*p_miit,rc,V,mode,scatter_mode); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
PetscErrorCode Core::set_global_VecCreateGhost(const MoFEMProblem *problem_ptr,RowColData rc,Vec V,InsertMode mode,ScatterMode scatter_mode) {
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
	SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INSONSISTENCY,"data inconsistency: nb. of dofs and declared nb. dofs in database");
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
PetscErrorCode Core::set_global_VecCreateGhost(const string &name,RowColData rc,Vec V,InsertMode mode,ScatterMode scatter_mode) {
  PetscFunctionBegin;
  typedef MoFEMProblem_multiIndex::index<Problem_mi_tag>::type moFEMProblems_by_name;
  moFEMProblems_by_name &moFEMProblems_set = moFEMProblems.get<Problem_mi_tag>();
  moFEMProblems_by_name::iterator p_miit = moFEMProblems_set.find(name);
  if(p_miit==moFEMProblems_set.end()) SETERRQ1(PETSC_COMM_SELF,1,"problem < %s > not found (top tip: check spelling)",name.c_str());
  ierr = set_global_VecCreateGhost(&*p_miit,rc,V,mode,scatter_mode); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
PetscErrorCode Core::set_other_local_VecCreateGhost(
  const MoFEMProblem *problem_ptr,const string& field_name,const string& cpy_field_name,RowColData rc,Vec V,InsertMode mode,ScatterMode scatter_mode,int verb) {
  PetscFunctionBegin;
  typedef NumeredDofMoFEMEntity_multiIndex::index<PetscLocalIdx_mi_tag>::type dofs_by_loc_petsc_index;
  dofs_by_loc_petsc_index *dofs;
  int nb_local_dofs,nb_ghost_dofs;
  switch (rc) {
    case ROW:
      nb_local_dofs = problem_ptr->get_nb_local_dofs_row();
      nb_ghost_dofs = problem_ptr->get_nb_ghost_dofs_row();
      dofs = const_cast<dofs_by_loc_petsc_index*>(&problem_ptr->numered_dofs_rows.get<PetscLocalIdx_mi_tag>());
      break;
    case COL:
      nb_local_dofs = problem_ptr->get_nb_local_dofs_col();
      nb_ghost_dofs = problem_ptr->get_nb_ghost_dofs_col();
      dofs = const_cast<dofs_by_loc_petsc_index*>(&problem_ptr->numered_dofs_cols.get<PetscLocalIdx_mi_tag>());
      break;
    default:
     SETERRQ(PETSC_COMM_SELF,MOFEM_NOT_IMPLEMENTED,"not implemented");
  }
  MoFEMField_multiIndex::index<FieldName_mi_tag>::type::iterator cpy_fit = moabFields.get<FieldName_mi_tag>().find(cpy_field_name);
  if(cpy_fit==moabFields.get<FieldName_mi_tag>().end()) {
    SETERRQ1(PETSC_COMM_SELF,MOFEM_NOT_FOUND,"cpy field < %s > not found, (top tip: check spelling)",cpy_field_name.c_str());
  }
  dofs_by_loc_petsc_index::iterator miit = dofs->lower_bound(0);
  if(miit==dofs->end()) {
    SETERRQ1(PETSC_COMM_SELF,MOFEM_NOT_FOUND,"cpy field < %s > not found, (top tip: check spelling)",field_name.c_str());
  }
  dofs_by_loc_petsc_index::iterator hi_miit = dofs->upper_bound(nb_local_dofs+nb_ghost_dofs);
  if(miit->get_space() != cpy_fit->get_space()) {
    SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INSONSISTENCY,"fiedls has to have same space");
  }
  if(miit->get_max_rank() != cpy_fit->get_max_rank()) {
    SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INSONSISTENCY,"fiedls has to have same rank");
  }
  switch (scatter_mode) {
    case SCATTER_REVERSE: {
      PetscScalar *array;
      VecGetArray(V,&array);
      switch (mode) {
	case INSERT_VALUES:
	  for(;miit!=hi_miit;miit++) {
	    if(miit->get_name()!=field_name) continue;
	    DofMoFEMEntity_multiIndex::index<Composite_Name_And_Ent_And_EndDofIdx_mi_tag>::type::iterator diiiit;
	    diiiit = dofsMoabField.get<Composite_Name_And_Ent_And_EndDofIdx_mi_tag>().find(boost::make_tuple(cpy_field_name,miit->get_ent(),miit->get_EntDofIdx()));
	    if(diiiit==dofsMoabField.get<Composite_Name_And_Ent_And_EndDofIdx_mi_tag>().end()) {
	      SETERRQ(PETSC_COMM_SELF,MOFEM_NOT_FOUND,"equivalalent dof does not exist, dof has to be creates, use set_other_global_VecCreateGhost to create dofs entries");
	    }
	    diiiit->get_FieldData() = array[miit->get_petsc_local_dof_idx()];
	    if(verb > 1) {
	      ostringstream ss;
	      ss << *diiiit << "set " << array[miit->get_petsc_local_dof_idx()] << endl;
	      PetscPrintf(PETSC_COMM_WORLD,ss.str().c_str());
	    }
	  }
	    //if(verb > 0) {
	      //cerr << "AAAAAAAAAAAA\n";
	      //ierr = check_number_of_ents_in_ents_field(cpy_field_name); CHKERRQ(ierr);
	    //}
	  break;
	default:
	  SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INSONSISTENCY,"not implemented");
      }
      ierr = VecRestoreArray(V,&array); CHKERRQ(ierr);
    }
    break;
    case SCATTER_FORWARD: {
	for(;miit!=hi_miit;miit++) {
	  if(miit->get_name()!=field_name) continue;
	  DofMoFEMEntity_multiIndex::index<Composite_Name_And_Ent_And_EndDofIdx_mi_tag>::type::iterator diiiit;
	  diiiit = dofsMoabField.get<Composite_Name_And_Ent_And_EndDofIdx_mi_tag>().find(boost::make_tuple(cpy_field_name,miit->get_ent(),miit->get_EntDofIdx()));
	  if(diiiit==dofsMoabField.get<Composite_Name_And_Ent_And_EndDofIdx_mi_tag>().end()) {
	    SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INSONSISTENCY,"no data to fill the vector (top tip: you want scatter forward or scatter reverse?)");
	  }
	  ierr = VecSetValue(V,miit->get_petsc_gloabl_dof_idx(),diiiit->get_FieldData(),mode); CHKERRQ(ierr);
	}
	ierr = VecAssemblyBegin(V); CHKERRQ(ierr);
	ierr = VecAssemblyEnd(V); CHKERRQ(ierr);
      } 
      break;  
    default:
     SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INSONSISTENCY,"not implemented");
  }


  PetscFunctionReturn(0);
}
PetscErrorCode Core::set_other_local_VecCreateGhost(
  const string &name,const string& field_name,const string& cpy_field_name,RowColData rc,Vec V,InsertMode mode,ScatterMode scatter_mode,int verb) {
  PetscFunctionBegin;
  if(verb==-1) verb = verbose;
  typedef MoFEMProblem_multiIndex::index<Problem_mi_tag>::type moFEMProblems_by_name;
  moFEMProblems_by_name &moFEMProblems_set = moFEMProblems.get<Problem_mi_tag>();
  moFEMProblems_by_name::iterator p_miit = moFEMProblems_set.find(name);
  if(p_miit==moFEMProblems_set.end()) SETERRQ1(PETSC_COMM_SELF,1,"problem < %s > not found",name.c_str());
  ierr = set_other_local_VecCreateGhost(&*p_miit,field_name,cpy_field_name,rc,V,mode,scatter_mode,verb); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
PetscErrorCode Core::set_other_global_VecCreateGhost(
  const string &name,const string& field_name,const string& cpy_field_name,RowColData rc,Vec V,InsertMode mode,ScatterMode scatter_mode,
  int verb) {
  PetscFunctionBegin;
  if(verb==-1) verb = verbose;
  ParallelComm* pcomm = ParallelComm::get_pcomm(&moab,MYPCOMM_INDEX);
  typedef MoFEMProblem_multiIndex::index<Problem_mi_tag>::type moFEMProblems_by_name;
  typedef NumeredDofMoFEMEntity_multiIndex::index<FieldName_mi_tag>::type dofs_by_name;
  moFEMProblems_by_name &moFEMProblems_set = moFEMProblems.get<Problem_mi_tag>();
  moFEMProblems_by_name::iterator p_miit = moFEMProblems_set.find(name);
  if(p_miit==moFEMProblems_set.end()) SETERRQ1(PETSC_COMM_SELF,1,"problem < %s > not found",name.c_str());
  dofs_by_name *dofs;
  DofIdx nb_dofs;
  switch (rc) {
    case ROW:
      nb_dofs = p_miit->get_nb_dofs_row();
      dofs = const_cast<dofs_by_name*>(&p_miit->numered_dofs_rows.get<FieldName_mi_tag>());
      break;
    case COL:
      nb_dofs = p_miit->get_nb_dofs_col();
      dofs = const_cast<dofs_by_name*>(&p_miit->numered_dofs_cols.get<FieldName_mi_tag>());
      break;
    default:
     SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INSONSISTENCY,"not implemented");
  }
  MoFEMField_multiIndex::index<FieldName_mi_tag>::type::iterator cpy_fit = moabFields.get<FieldName_mi_tag>().find(cpy_field_name);
  if(cpy_fit==moabFields.get<FieldName_mi_tag>().end()) {
    SETERRQ1(PETSC_COMM_SELF,MOFEM_NOT_FOUND,"cpy field < %s > not found, (top tip: check spelling)",cpy_field_name.c_str());
  }
  dofs_by_name::iterator miit = dofs->lower_bound(field_name);
  if(miit==dofs->end()) {
    SETERRQ1(PETSC_COMM_SELF,MOFEM_NOT_FOUND,"cpy field < %s > not found, (top tip: check spelling)",field_name.c_str());
  }
  dofs_by_name::iterator hi_miit = dofs->upper_bound(field_name);
  if(miit->get_space() != cpy_fit->get_space()) {
    SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INSONSISTENCY,"fiedls has to have same space");
  }
  if(miit->get_max_rank() != cpy_fit->get_max_rank()) {
    SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INSONSISTENCY,"fiedls has to have same rank");
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
      if(size!=nb_dofs) SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INSONSISTENCY,"data inconsistency: nb. of dofs and declared nb. dofs in database");
      PetscScalar *array;
      VecGetArray(V_glob,&array);
      switch (mode) {
	case INSERT_VALUES:
	  for(;miit!=hi_miit;miit++) {
	    if(miit->get_petsc_gloabl_dof_idx()>=size) {
	      SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INSONSISTENCY,"data inconsistency: nb. of dofs and declared nb. dofs in database");
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
		SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INSONSISTENCY,ss.str().c_str());
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
	      if(diiiit==dofsMoabField.get<Composite_Name_And_Ent_And_EndDofIdx_mi_tag>().end()) SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INSONSISTENCY,"data inconsistency");
	    }
	    diiiit->get_FieldData() = array[miit->get_petsc_gloabl_dof_idx()];
	    if(verb > 1) {
	      ostringstream ss;
	      ss << *diiiit << "set " << array[miit->get_petsc_gloabl_dof_idx()] << endl;
	      PetscPrintf(PETSC_COMM_WORLD,ss.str().c_str());
	    }
	  }
	    //if(verb > 0) {
	      //cerr << "AAAAAAAAAAAA\n";
	      //ierr = check_number_of_ents_in_ents_field(cpy_field_name); CHKERRQ(ierr);
	    //}
	  break;
	default:
	  SETERRQ(PETSC_COMM_SELF,MOFEM_NOT_IMPLEMENTED,"not implemented");
      }
      ierr = VecRestoreArray(V_glob,&array); CHKERRQ(ierr);
      ierr = VecDestroy(&V_glob); CHKERRQ(ierr);
      ierr = VecScatterDestroy(&ctx); CHKERRQ(ierr);
    }
    break;
    case SCATTER_FORWARD: {
	for(;miit!=hi_miit;miit++) {
	  if(pcomm->rank()!=miit->get_part()) continue;
	  DofMoFEMEntity_multiIndex::index<Composite_Name_And_Ent_And_EndDofIdx_mi_tag>::type::iterator diiiit;
	  diiiit = dofsMoabField.get<Composite_Name_And_Ent_And_EndDofIdx_mi_tag>().find(boost::make_tuple(cpy_field_name,miit->get_ent(),miit->get_EntDofIdx()));
	  if(diiiit==dofsMoabField.get<Composite_Name_And_Ent_And_EndDofIdx_mi_tag>().end()) {
	    SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INSONSISTENCY,"no data to fill the vector (top tip: you want scatter forward or scatter reverse?)");
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

}

