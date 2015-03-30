/** \file FieldInterfaceCore.cpp
 * \brief Myltindex containes, data structures and other low-level functions 
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

PetscErrorCode Core::build_partitioned_problem(const string &name,int verb) {
  PetscFunctionBegin;
  if(verb==-1) verb = verbose;
  if(!(*build_MoFEM&(1<<0))) SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCT,"fields not build");
  if(!(*build_MoFEM&(1<<1))) SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCT,"FEs not build");
  if(!(*build_MoFEM&(1<<2))) SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCT,"entFEAdjacencies not build");
  const MoFEMProblem *problem_ptr;
  ierr = get_problem(name,&problem_ptr); CHKERRQ(ierr);
  ierr = build_partitioned_problem(const_cast<MoFEMProblem*>(problem_ptr),verb); CHKERRQ(ierr);
  *build_MoFEM |= 1<<3;
  *build_MoFEM |= 1<<4; // It is assumed that user who uses this function knows what he is doing
  PetscFunctionReturn(0);
}
PetscErrorCode Core::build_partitioned_problem(MoFEMProblem *problem_ptr,int verb) {
  PetscFunctionBegin;
  if(verb==-1) verb = verbose;
  if(problem_ptr->get_BitRefLevel().none()) {
    SETERRQ1(PETSC_COMM_SELF,1,"problem <%s> refinement level not set",problem_ptr->get_name().c_str());
  }
  //zero finite elements
  problem_clear_numered_finiteElementsPtr_change().operator()(*problem_ptr);
  //bool success = moFEMProblems.modify(p_miit,problem_clear_numered_finiteElementsPtr_change());
  //if(!success) SETERRQ(PETSC_COMM_SELF,MOFEM_OPERATION_UNSUCCESSFUL,"modification unsuccessful");
  //miit2 iterator for finite elements
  EntMoFEMFiniteElement_multiIndex::iterator miit2 = finiteElementsMoFEMEnts.begin();
  EntMoFEMFiniteElement_multiIndex::iterator hi_miit2 = finiteElementsMoFEMEnts.end();
  DofMoFEMEntity_multiIndex_active_view dofs_rows,dofs_cols;
  dofs_rows.clear();
  dofs_cols.clear();
  EntMoFEMFiniteElement_multiIndex::iterator miit3 = miit2;
  //iterate all finite elemen entities in database
  for(;miit3!=hi_miit2;miit3++) {
    //if element is in problem
    if((miit3->get_id()&problem_ptr->get_BitFEId()).any()) {
      //if finite element bit level has all refined bits sets
      if((miit3->get_BitRefLevel()&problem_ptr->get_BitRefLevel())==problem_ptr->get_BitRefLevel()) {
	//get dof uids for rows and columns
	ierr = miit3->get_MoFEMFiniteElement_row_dof_view(dofsMoabField,dofs_rows); CHKERRQ(ierr);
	ierr = miit3->get_MoFEMFiniteElement_col_dof_view(dofsMoabField,dofs_cols); CHKERRQ(ierr);
      }
    }
  }
  //zero rows
  *(problem_ptr->tag_nbdof_data_row) = 0;
  *(problem_ptr->tag_local_nbdof_data_row) = 0;
  *(problem_ptr->tag_ghost_nbdof_data_row) = 0;
  problem_ptr->numered_dofs_rows.clear();
  //zero cols
  (*problem_ptr->tag_nbdof_data_col) = 0;
  (*problem_ptr->tag_local_nbdof_data_col) = 0;
  (*problem_ptr->tag_ghost_nbdof_data_col) = 0;
  problem_ptr->numered_dofs_cols.clear();
  //add dofs for rows
  DofMoFEMEntity_multiIndex_active_view::nth_index<1>::type::iterator miit4,hi_miit4;
  miit4 = dofs_rows.get<1>().lower_bound(1);
  hi_miit4 = dofs_rows.get<1>().upper_bound(1);
  for(;miit4!=hi_miit4;miit4++) {
    if(((*miit4)->get_BitRefLevel()&problem_ptr->get_DofMask_BitRefLevel())!=(*miit4)->get_BitRefLevel()) {
      continue;
    }
    NumeredDofMoFEMEntity dof(&**miit4);
    pair<NumeredDofMoFEMEntity_multiIndex::iterator,bool> p = problem_ptr->numered_dofs_rows.insert(dof); 
    int owner_proc = p.first->get_owner_proc();
    if(owner_proc<0) {
      SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCT,"data inconsistency");
    }
    bool success = problem_ptr->numered_dofs_rows.modify(p.first,NumeredDofMoFEMEntity_part_change(owner_proc,-1));
     if(!success) SETERRQ(PETSC_COMM_SELF,MOFEM_OPERATION_UNSUCCESSFUL,"modification unsuccessful");
  }
  //add dofs for cols
  DofMoFEMEntity_multiIndex_active_view::nth_index<1>::type::iterator miit5,hi_miit5;
  miit5 = dofs_cols.get<1>().lower_bound(1);
  hi_miit5 = dofs_cols.get<1>().upper_bound(1);
  for(;miit5!=hi_miit5;miit5++) {
    if(((*miit5)->get_BitRefLevel()&problem_ptr->get_DofMask_BitRefLevel())!=(*miit5)->get_BitRefLevel()) {
      continue;
    }
    NumeredDofMoFEMEntity dof(&**miit5);
    pair<NumeredDofMoFEMEntity_multiIndex::iterator,bool> p = problem_ptr->numered_dofs_cols.insert(dof); 
    int owner_proc = p.first->get_owner_proc();
    if(owner_proc<0) {
      SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCT,"data inconsistency");
    }
    bool success = problem_ptr->numered_dofs_cols.modify(p.first,NumeredDofMoFEMEntity_part_change(owner_proc,-1));
    if(!success) SETERRQ(PETSC_COMM_SELF,MOFEM_OPERATION_UNSUCCESSFUL,"modification unsuccessful");
  }
  //number dofs
  NumeredDofMoFEMEntity_multiIndex *numered_dofs[] = { 
    &(problem_ptr->numered_dofs_rows), &(problem_ptr->numered_dofs_cols) };
  DofIdx *nb_global_dofs[] = {
    problem_ptr->tag_nbdof_data_row, problem_ptr->tag_nbdof_data_col };
  DofIdx *nb_local_dofs[] = {
    problem_ptr->tag_local_nbdof_data_row, problem_ptr->tag_local_nbdof_data_col };
  vector<GlobalUId> uids;
  for(int ss = 0;ss<2;ss++) {
    NumeredDofMoFEMEntity_multiIndex::index<Part_mi_tag>::type::iterator diit,hi_diit;
    diit = numered_dofs[ss]->get<Part_mi_tag>().lower_bound(rAnk);
    hi_diit = numered_dofs[ss]->get<Part_mi_tag>().upper_bound(rAnk);
    uids.resize(distance(diit,hi_diit));
    int &local_idx = *(nb_local_dofs[ss]);
    local_idx = 0;
    for(int ii = 0;diit!=hi_diit;diit++,ii++) {
      bool success;
      success = numered_dofs[ss]->modify(
	numered_dofs[ss]->project<0>(diit),NumeredDofMoFEMEntity_local_idx_change(local_idx++));
      if(!success) SETERRQ(PETSC_COMM_SELF,MOFEM_OPERATION_UNSUCCESSFUL,"modification unsuccessful");
      uids[ii] = diit->get_global_unique_id();
    }
    /// FIXME: this should be done with collective communication,
    {
      IS is,isout;
      int *petscint_ptr = (int*)&*uids.begin();
      ierr = ISCreateGeneral(comm,
	uids.size()*sizeof(GlobalUId)/sizeof(int),
	petscint_ptr,PETSC_USE_POINTER,&is); CHKERRQ(ierr);
      ierr = ISAllGather(is,&isout); CHKERRQ(ierr);
      int isout_size;
      ierr = ISGetSize(isout,&isout_size); CHKERRQ(ierr);
      isout_size = isout_size*sizeof(int)/sizeof(GlobalUId);
      *(nb_global_dofs[ss]) = isout_size;
      {	
	const int *ptr;
	ierr = ISGetIndices(isout,&ptr); CHKERRQ(ierr);
	GlobalUId *uid_ptr = const_cast<GlobalUId*>((const GlobalUId*)ptr);
	NumeredDofMoFEMEntity_multiIndex::iterator diit,hi_diit;
	diit = numered_dofs[ss]->begin();
	hi_diit = numered_dofs[ss]->end();
	for(;diit!=hi_diit;diit++) {
	  GlobalUId uid = diit->get_global_unique_id();
	  GlobalUId *ptr;
	  ptr = find(uid_ptr,uid_ptr+isout_size,uid);
	  if(ptr == uid_ptr+isout_size) {
	    SETERRQ(PETSC_COMM_SELF,MOFEM_OPERATION_UNSUCCESSFUL,"could not find uid");
	  }
	  DofIdx idx = distance(uid_ptr,ptr);
	  int part = diit->get_part();
	  bool success;
	  success = numered_dofs[ss]->modify(
	    diit,NumeredDofMoFEMEntity_mofem_index_change(idx));
	  if(!success) SETERRQ(PETSC_COMM_SELF,MOFEM_OPERATION_UNSUCCESSFUL,"modification unsuccessful");
	  success = numered_dofs[ss]->modify(
	    diit,NumeredDofMoFEMEntity_part_change(part,idx));
	  if(!success) SETERRQ(PETSC_COMM_SELF,MOFEM_OPERATION_UNSUCCESSFUL,"modification unsuccessful");
	}
	ierr = ISRestoreIndices(isout,&ptr); CHKERRQ(ierr);
      }
      ierr = ISDestroy(&is); CHKERRQ(ierr);
      ierr = ISDestroy(&isout); CHKERRQ(ierr);
    }
  }
  ierr = print_partitioned_problem(problem_ptr,verb); CHKERRQ(ierr);
  if(verb>2) {
    NumeredDofMoFEMEntity_multiIndex::iterator it,hi_it;
    for(int ss = 0;ss<2;ss++) {
      it = numered_dofs[ss]->begin();
      hi_it = numered_dofs[ss]->end();
      for(;it!=hi_it;it++) {
	ostringstream zz;
	zz << ss << " " << "rank " << rAnk << " ";
	zz << *it << endl;
	PetscSynchronizedPrintf(comm,zz.str().c_str());
      }
    }
  }
  *build_MoFEM |= 1<<3;
  *build_MoFEM |= 1<<4; // It is assumed that user who uses this function knows what he is doing
  PetscFunctionReturn(0);
}
PetscErrorCode Core::build_partitioned_problems(int verb) {
  PetscFunctionBegin;
  if(verb==-1) verb = verbose;
  DofMoFEMEntity_multiIndex_active_view dofs_rows,dofs_cols;
  MoFEMProblem_multiIndex::iterator p_miit = moFEMProblems.begin();
  for(;p_miit!=moFEMProblems.end();p_miit++) {
    ierr = build_partitioned_problem(const_cast<MoFEMProblem*>(&*p_miit),verb); CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}
PetscErrorCode Core::build_problem(MoFEMProblem *problem_ptr,int verb) {
  PetscFunctionBegin;
  // Note: Only allowe changes on problem_ptr structure which not influence multindex
  // indexing are allowd.
  if(verb==-1) verb = verbose;
  if(problem_ptr->get_BitRefLevel().none()) {
    SETERRQ1(PETSC_COMM_SELF,1,"problem <%s> refinement level not set",problem_ptr->get_name().c_str());
  }
  if(problem_ptr->get_BitRefLevel().none()) {
    SETERRQ1(PETSC_COMM_SELF,1,"problem <%s> refinement level not set",problem_ptr->get_name().c_str());
  }
  //zero finite elements
  problem_ptr->numeredFiniteElements.clear();
  //miit2 iterator for finite elements
  EntMoFEMFiniteElement_multiIndex::iterator miit2 = finiteElementsMoFEMEnts.begin();
  EntMoFEMFiniteElement_multiIndex::iterator hi_miit2 = finiteElementsMoFEMEnts.end();
  DofMoFEMEntity_multiIndex_active_view dofs_rows,dofs_cols;
  EntMoFEMFiniteElement_multiIndex::iterator miit3 = miit2;
  //iterate all finite elemen entities in database
  for(;miit3!=hi_miit2;miit3++) {
    //if element is in problem
    if((miit3->get_id()&problem_ptr->get_BitFEId()).any()) {
	//if finite element bit level has all refined bits sets
	if((miit3->get_BitRefLevel()&problem_ptr->get_BitRefLevel())==problem_ptr->get_BitRefLevel()) {
	  //get dof uids for rows and columns
	  ierr = miit3->get_MoFEMFiniteElement_row_dof_view(dofsMoabField,dofs_rows); CHKERRQ(ierr);
	  ierr = miit3->get_MoFEMFiniteElement_col_dof_view(dofsMoabField,dofs_cols); CHKERRQ(ierr);
	}
    }
  }
  //zero rows
  *problem_ptr->tag_nbdof_data_row = 0;
  *problem_ptr->tag_local_nbdof_data_row = 0;
  *problem_ptr->tag_ghost_nbdof_data_row = 0;
  problem_ptr->numered_dofs_rows.clear();
  //zero cols
  *problem_ptr->tag_nbdof_data_col = 0;
  *problem_ptr->tag_local_nbdof_data_col = 0;
  *problem_ptr->tag_ghost_nbdof_data_col = 0;
  problem_ptr->numered_dofs_cols.clear();
  //add dofs for rows
  DofMoFEMEntity_multiIndex_active_view::nth_index<1>::type::iterator miit4,hi_miit4;
  miit4 = dofs_rows.get<1>().lower_bound(1);
  hi_miit4 = dofs_rows.get<1>().upper_bound(1);
  for(;miit4!=hi_miit4;miit4++) {
    if(((*miit4)->get_BitRefLevel()&problem_ptr->get_DofMask_BitRefLevel())!=(*miit4)->get_BitRefLevel()) {
      continue;
    }
    problem_add_row_dof(&**miit4).operator()(*problem_ptr);
    //bool success = moFEMProblems.modify(problem_ptr,problem_add_row_dof(&**miit4));
    //if(!success) SETERRQ(PETSC_COMM_SELF,MOFEM_OPERATION_UNSUCCESSFUL,"modification unsuccessful");
  }
  //add dofs for cols
  DofMoFEMEntity_multiIndex_active_view::nth_index<1>::type::iterator miit5,hi_miit5;
  miit5 = dofs_cols.get<1>().lower_bound(1);
  hi_miit5 = dofs_cols.get<1>().upper_bound(1);
  for(;miit5!=hi_miit5;miit5++) {
    if(((*miit5)->get_BitRefLevel()&problem_ptr->get_DofMask_BitRefLevel())!=(*miit5)->get_BitRefLevel()) {
      continue;
    }
    problem_add_col_dof(&**miit5).operator()(*problem_ptr);
    //success = moFEMProblems.modify(problem_ptr,problem_add_col_dof(&**miit5));
    //if(!success) SETERRQ(PETSC_COMM_SELF,MOFEM_OPERATION_UNSUCCESSFUL,"modification unsuccessful");
  }
  //number dofs on rows and columns
  problem_row_number_change().operator()(*problem_ptr);
  problem_col_number_change().operator()(*problem_ptr);
  //success = moFEMProblems.modify(problem_ptr,problem_row_number_change());
  //if(!success) SETERRQ(PETSC_COMM_SELF,MOFEM_OPERATION_UNSUCCESSFUL,"modification unsuccessful");
  //success = moFEMProblems.modify(problem_ptr,problem_col_number_change());
  //if(!success) SETERRQ(PETSC_COMM_SELF,MOFEM_OPERATION_UNSUCCESSFUL,"modification unsuccessful");
  //job done, some debugging and postprocessing
  if(verbose>0) {
    PetscSynchronizedPrintf(comm,"Problem %s Nb. rows %u Nb. cols %u\n",
	problem_ptr->get_name().c_str(),
	problem_ptr->numered_dofs_rows.size(),problem_ptr->numered_dofs_cols.size());
  }
  if(verb>1) {
    EntMoFEMFiniteElement_multiIndex::iterator miit_ss = miit2;
    ostringstream ss;
    ss << "rank " << rAnk << " ";
    ss << "FEs data for problem " << *problem_ptr << endl;
    for(;miit_ss!=hi_miit2;miit_ss++) {
	ss << "rank " << rAnk << " ";
	ss << *miit_ss << endl;
    }
    ss << "rank " << rAnk << " ";
    ss << "FEs row dofs "<< *problem_ptr << " Nb. row dof " << problem_ptr->get_nb_dofs_row() << endl;
    NumeredDofMoFEMEntity_multiIndex::iterator miit_dd_row = problem_ptr->numered_dofs_rows.begin();
    for(;miit_dd_row!=problem_ptr->numered_dofs_rows.end();miit_dd_row++) {
	ss << "rank " << rAnk << " ";
	ss<<*miit_dd_row<<endl;
    }
    ss << "rank " << rAnk << " ";
    ss << "FEs col dofs "<< *problem_ptr << " Nb. col dof " << problem_ptr->get_nb_dofs_col() << endl;
    NumeredDofMoFEMEntity_multiIndex::iterator miit_dd_col = problem_ptr->numered_dofs_cols.begin();
    for(;miit_dd_col!=problem_ptr->numered_dofs_cols.end();miit_dd_col++) {
	ss << "rank " << rAnk << " ";
	ss<<*miit_dd_col<<endl;
    }
    PetscSynchronizedPrintf(comm,ss.str().c_str());
  }
  if(verb>0) {
    PetscSynchronizedFlush(comm,PETSC_STDOUT); 
  }
  *build_MoFEM |= 1<<3; // It is assumed that user who uses this function knows what he is doing
  PetscFunctionReturn(0);
}
PetscErrorCode Core::build_problem(const string &problem_name,int verb) {
  PetscFunctionBegin;
  if(verb==-1) verb = verbose;
  if(!(*build_MoFEM&(1<<0))) SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCT,"fields not build");
  if(!(*build_MoFEM&(1<<1))) SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCT,"FEs not build");
  if(!(*build_MoFEM&(1<<2))) SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCT,"entFEAdjacencies not build");
  const MoFEMProblem *problem_ptr;
  ierr = get_problem(problem_name,&problem_ptr); CHKERRQ(ierr);
  ierr = build_problem(const_cast<MoFEMProblem*>(problem_ptr),verb); CHKERRQ(ierr);
  *build_MoFEM |= 1<<3; // It is assumed that user who uses this function knows what he is doing
  PetscFunctionReturn(0);
}
PetscErrorCode Core::build_problems(int verb) {
  PetscFunctionBegin;
  if(verb==-1) verb = verbose;
  if(!(*build_MoFEM&(1<<0))) SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCT,"fields not build");
  if(!(*build_MoFEM&(1<<1))) SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCT,"FEs not build");
  if(!(*build_MoFEM&(1<<2))) SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCT,"entFEAdjacencies not build");
  //iterate problems
  MoFEMProblem_multiIndex::iterator p_miit = moFEMProblems.begin();
  for(;p_miit!=moFEMProblems.end();p_miit++) {
    MoFEMProblem *problem_ptr =  const_cast<MoFEMProblem*>(&*p_miit);
    ierr = build_problem(problem_ptr,verb); CHKERRQ(ierr);
  }
  *build_MoFEM |= 1<<3;
  PetscFunctionReturn(0);
}
PetscErrorCode Core::clear_problems(int verb) {
  PetscFunctionBegin;
  if(verb==-1) verb = verbose;
  MoFEMProblem_multiIndex::iterator p_miit = moFEMProblems.begin();
  //iterate problems
  for(;p_miit!=moFEMProblems.end();p_miit++) {
    //zero rows
    bool success = moFEMProblems.modify(p_miit,problem_zero_nb_rows_change());
    if(!success) SETERRQ(PETSC_COMM_SELF,MOFEM_OPERATION_UNSUCCESSFUL,"modification unsuccessful");
    //zero cols
    success = moFEMProblems.modify(p_miit,problem_zero_nb_cols_change());
    if(!success) SETERRQ(PETSC_COMM_SELF,MOFEM_OPERATION_UNSUCCESSFUL,"modification unsuccessful");
    //clear finite elements
    success = moFEMProblems.modify(p_miit,problem_clear_numered_finiteElementsPtr_change());
    if(!success) SETERRQ(PETSC_COMM_SELF,MOFEM_OPERATION_UNSUCCESSFUL,"modification unsuccessful");
  }
  PetscFunctionReturn(0);
}
PetscErrorCode Core::simple_partition_problem(const string &name,int verb) {
  PetscFunctionBegin;
  if(verb==-1) verb = verbose;
  if(!(*build_MoFEM&(1<<0))) SETERRQ(PETSC_COMM_SELF,MOFEM_NOT_FOUND,"fields not build");
  if(!(*build_MoFEM&(1<<1))) SETERRQ(PETSC_COMM_SELF,MOFEM_NOT_FOUND,"FEs not build");
  if(!(*build_MoFEM&(1<<2))) SETERRQ(PETSC_COMM_SELF,MOFEM_NOT_FOUND,"entFEAdjacencies not build");
  if(!(*build_MoFEM&(1<<3))) SETERRQ(PETSC_COMM_SELF,MOFEM_NOT_FOUND,"moFEMProblems not build");
  if(verb>0) {
    PetscPrintf(comm,"Simple partition problem %s\n",name.c_str());
  }
  // find p_miit
  typedef MoFEMProblem_multiIndex::index<Problem_mi_tag>::type MoFEMProblem_multiIndex_by_name;
  MoFEMProblem_multiIndex_by_name &moFEMProblems_set = moFEMProblems.get<Problem_mi_tag>();
  MoFEMProblem_multiIndex_by_name::iterator p_miit = moFEMProblems_set.find(name);
  if(p_miit==moFEMProblems_set.end()) SETERRQ1(PETSC_COMM_SELF,1,"problem < %s > is not found (top tip: check spelling)",name.c_str());
  typedef boost::multi_index::index<NumeredDofMoFEMEntity_multiIndex,Idx_mi_tag>::type NumeredDofMoFEMEntitys_by_idx;
  NumeredDofMoFEMEntitys_by_idx &dofs_row_by_idx = const_cast<NumeredDofMoFEMEntitys_by_idx&>(p_miit->numered_dofs_rows.get<Idx_mi_tag>());
  NumeredDofMoFEMEntitys_by_idx &dofs_col_by_idx = const_cast<NumeredDofMoFEMEntitys_by_idx&>(p_miit->numered_dofs_cols.get<Idx_mi_tag>());
  boost::multi_index::index<NumeredDofMoFEMEntity_multiIndex,Idx_mi_tag>::type::iterator miit_row,hi_miit_row;
  boost::multi_index::index<NumeredDofMoFEMEntity_multiIndex,Idx_mi_tag>::type::iterator miit_col,hi_miit_col;
  DofIdx &nb_row_local_dofs = *((DofIdx*)p_miit->tag_local_nbdof_data_row);
  DofIdx &nb_row_ghost_dofs = *((DofIdx*)p_miit->tag_ghost_nbdof_data_row);
  nb_row_local_dofs = 0;
  nb_row_ghost_dofs = 0;
  DofIdx &nb_col_local_dofs = *((DofIdx*)p_miit->tag_local_nbdof_data_col);
  DofIdx &nb_col_ghost_dofs = *((DofIdx*)p_miit->tag_ghost_nbdof_data_col);
  nb_col_local_dofs = 0;
  nb_col_ghost_dofs = 0;
  //get row range of local indices
  DofIdx nb_dofs_row = dofs_row_by_idx.size();
  PetscLayout layout_row;
  ierr = PetscLayoutCreate(comm,&layout_row); CHKERRQ(ierr);
  ierr = PetscLayoutSetBlockSize(layout_row,1); CHKERRQ(ierr);
  ierr = PetscLayoutSetSize(layout_row,nb_dofs_row); CHKERRQ(ierr);
  ierr = PetscLayoutSetUp(layout_row); CHKERRQ(ierr);
  const int *ranges_row;
  ierr = PetscLayoutGetRanges(layout_row,&ranges_row); CHKERRQ(ierr);
  //get col range of local indices
  DofIdx nb_dofs_col = dofs_col_by_idx.size();
  PetscLayout layout_col;
  ierr = PetscLayoutCreate(comm,&layout_col); CHKERRQ(ierr);
  ierr = PetscLayoutSetBlockSize(layout_col,1); CHKERRQ(ierr);
  ierr = PetscLayoutSetSize(layout_col,nb_dofs_col); CHKERRQ(ierr);
  ierr = PetscLayoutSetUp(layout_col); CHKERRQ(ierr);
  const int *ranges_col;
  ierr = PetscLayoutGetRanges(layout_col,&ranges_col); CHKERRQ(ierr);
  for(unsigned int part = 0;part<(unsigned int)sIze;part++) {
    miit_row = dofs_row_by_idx.lower_bound(ranges_row[part]);
    hi_miit_row = dofs_row_by_idx.lower_bound(ranges_row[part+1]);
    if(distance(miit_row,hi_miit_row) != ranges_row[part+1]-ranges_row[part]) {
      SETERRQ4(PETSC_COMM_SELF,PETSC_ERR_ARG_SIZ,
	  "data inconsistency, distance(miit_row,hi_miit_row) != rend - rstart (%d != %d - %d = %d) ",
	  distance(miit_row,hi_miit_row),ranges_row[part+1],ranges_row[part],ranges_row[part+1]-ranges_row[part]);
    }
    // loop rows
    for(;miit_row!=hi_miit_row;miit_row++) {
      bool success = dofs_row_by_idx.modify(miit_row,NumeredDofMoFEMEntity_part_change(part,miit_row->dof_idx));
      if(!success) SETERRQ(PETSC_COMM_SELF,MOFEM_OPERATION_UNSUCCESSFUL,"modification unsuccessful");
      if(part == (unsigned int)rAnk) {
	success = dofs_row_by_idx.modify(miit_row,NumeredDofMoFEMEntity_local_idx_change(nb_row_local_dofs++));
	if(!success) SETERRQ(PETSC_COMM_SELF,MOFEM_OPERATION_UNSUCCESSFUL,"modification unsuccessful");
      }
    }
    miit_col = dofs_col_by_idx.lower_bound(ranges_col[part]);
    hi_miit_col = dofs_col_by_idx.lower_bound(ranges_col[part+1]);
    if(distance(miit_col,hi_miit_col) != ranges_col[part+1]-ranges_col[part]) {
      SETERRQ4(PETSC_COMM_SELF,PETSC_ERR_ARG_SIZ,
	  "data inconsistency, distance(miit_col,hi_miit_col) != rend - rstart (%d != %d - %d = %d) ",
	  distance(miit_col,hi_miit_col),ranges_col[part+1],ranges_col[part],ranges_col[part+1]-ranges_col[part]);
    }
    // loop cols
    for(;miit_col!=hi_miit_col;miit_col++) {
      bool success = dofs_col_by_idx.modify(miit_col,NumeredDofMoFEMEntity_part_change(part,miit_col->dof_idx));
      if(!success) SETERRQ(PETSC_COMM_SELF,MOFEM_OPERATION_UNSUCCESSFUL,"modification unsuccessful");
      if(part == (unsigned int)rAnk) {
	success = dofs_col_by_idx.modify(miit_col,NumeredDofMoFEMEntity_local_idx_change(nb_col_local_dofs++));
	if(!success) SETERRQ(PETSC_COMM_SELF,MOFEM_OPERATION_UNSUCCESSFUL,"modification unsuccessful");
      }
    }
  }
  ierr = PetscLayoutDestroy(&layout_row); CHKERRQ(ierr);
  ierr = PetscLayoutDestroy(&layout_col); CHKERRQ(ierr);
  ierr = print_partitioned_problem(&*p_miit,verb); CHKERRQ(ierr);
  *build_MoFEM |= 1<<4;
  PetscFunctionReturn(0);
}
PetscErrorCode Core::compose_problem(const string &name,const string &problem_for_rows,bool copy_rows,const string &problem_for_cols,bool copy_cols,int verb) {
  PetscFunctionBegin;
  if(verb==-1) verb = verbose;
  if(!(*build_MoFEM&(1<<0))) SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCT,"fields not build");
  if(!(*build_MoFEM&(1<<1))) SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCT,"FEs not build");
  if(!(*build_MoFEM&(1<<2))) SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCT,"entFEAdjacencies not build");
  if(!(*build_MoFEM&(1<<3))) SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCT,"moFEMProblems not build");
  typedef MoFEMProblem_multiIndex::index<Problem_mi_tag>::type MoFEMProblem_multiIndex_by_name;
  typedef NumeredDofMoFEMEntity_multiIndex::index<Unique_mi_tag>::type NumeredDofMoFEMEntitys_by_uid;
  //find p_miit
  MoFEMProblem_multiIndex_by_name &moFEMProblems_set = moFEMProblems.get<Problem_mi_tag>();
  MoFEMProblem_multiIndex_by_name::iterator p_miit = moFEMProblems_set.find(name);
  if(p_miit==moFEMProblems_set.end()) SETERRQ1(PETSC_COMM_SELF,1,"problem with name < %s > not defined (top tip check spelling)",name.c_str());
  if(verb>0) {
    PetscPrintf(comm,"Compose problem %s from rows of %s and columns of %s\n",
      p_miit->get_name().c_str(),problem_for_rows.c_str(),problem_for_cols.c_str());
  }
  //find p_miit_row
  MoFEMProblem_multiIndex_by_name::iterator p_miit_row = moFEMProblems_set.find(problem_for_rows);
  if(p_miit_row==moFEMProblems_set.end()) SETERRQ1(PETSC_COMM_SELF,1,"problem with name < %s > not defined (top tip check spelling)",problem_for_rows.c_str());
  const NumeredDofMoFEMEntity_multiIndex &dofs_row = p_miit_row->numered_dofs_rows;
  //find p_mit_col
  MoFEMProblem_multiIndex_by_name::iterator p_miit_col = moFEMProblems_set.find(problem_for_cols);
  if(p_miit_col==moFEMProblems_set.end()) SETERRQ1(PETSC_COMM_SELF,1,"problem with name < %s > not defined (top tip check spelling)",problem_for_cols.c_str());
  const NumeredDofMoFEMEntity_multiIndex &dofs_col = p_miit_col->numered_dofs_cols;
  bool copy[] = { copy_rows, copy_cols };
  NumeredDofMoFEMEntity_multiIndex* composed_dofs[] = {
    const_cast<NumeredDofMoFEMEntity_multiIndex*>(&p_miit->numered_dofs_rows), 
    const_cast<NumeredDofMoFEMEntity_multiIndex*>(&p_miit->numered_dofs_cols) };
  int* nb_local_dofs[] = { p_miit->tag_local_nbdof_data_row, p_miit->tag_local_nbdof_data_col };
  int* nb_dofs[] = { p_miit->tag_nbdof_data_row, p_miit->tag_nbdof_data_col };
  const NumeredDofMoFEMEntity_multiIndex* copied_dofs[] = { &dofs_row, &dofs_col };
  for(int ss = 0; ss<2;ss++) {
    // build indices 
    *nb_local_dofs[ss] = 0;
    if(!copy[ss]) {
      NumeredDofMoFEMEntitys_by_uid &dofs_by_uid 
	= const_cast<NumeredDofMoFEMEntitys_by_uid&>(copied_dofs[ss]->get<Unique_mi_tag>());
      for(NumeredDofMoFEMEntity_multiIndex::iterator
	dit = composed_dofs[ss]->begin();dit!=composed_dofs[ss]->end();dit++) {
	NumeredDofMoFEMEntitys_by_uid::iterator diit = dofs_by_uid.find(dit->get_global_unique_id());
	if(diit==dofs_by_uid.end()) {
	  SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCT,"data inconsistency, could not find dof in composite problem");
	}
	int part_number = diit->get_part(); // get part number
	int petsc_global_dof = dit->get_dof_idx();
	bool success;
	success = composed_dofs[ss]->modify(dit,NumeredDofMoFEMEntity_part_change(part_number,petsc_global_dof));
	if(!success) SETERRQ(PETSC_COMM_SELF,MOFEM_OPERATION_UNSUCCESSFUL,"modification unsuccessful");
	if(dit->get_part() == (unsigned int)rAnk) {
	  success =composed_dofs[ss]->modify(dit,NumeredDofMoFEMEntity_local_idx_change((*nb_local_dofs[ss])++));
	  if(!success) SETERRQ(PETSC_COMM_SELF,MOFEM_OPERATION_UNSUCCESSFUL,"modification unsuccessful");
	}
      }
    } else {
      for(NumeredDofMoFEMEntity_multiIndex::iterator
	dit = copied_dofs[ss]->begin();dit!=copied_dofs[ss]->end();dit++) {
	pair<NumeredDofMoFEMEntity_multiIndex::iterator,bool> p;
	p = composed_dofs[ss]->insert(NumeredDofMoFEMEntity(dit->get_DofMoFEMEntity_ptr())); 
	if(p.second) {
	  (*nb_dofs[ss])++;
	}
	int dof_idx = dit->get_dof_idx();
	int part_number = dit->get_part(); // get part number
	int petsc_global_dof = dit->get_petsc_gloabl_dof_idx();
	bool success;
	success = composed_dofs[ss]->modify(p.first,NumeredDofMoFEMEntity_mofem_index_change(dof_idx));
	if(!success) SETERRQ(PETSC_COMM_SELF,MOFEM_OPERATION_UNSUCCESSFUL,"modification unsuccessful");
	success = composed_dofs[ss]->modify(p.first,NumeredDofMoFEMEntity_part_change(part_number,petsc_global_dof));
	if(!success) SETERRQ(PETSC_COMM_SELF,MOFEM_OPERATION_UNSUCCESSFUL,"modification unsuccessful");
	if(p.first->get_part() == (unsigned int)rAnk) {
	  success =composed_dofs[ss]->modify(p.first,NumeredDofMoFEMEntity_local_idx_change((*nb_local_dofs[ss])++));
	  if(!success) SETERRQ(PETSC_COMM_SELF,MOFEM_OPERATION_UNSUCCESSFUL,"modification unsuccessful");
	}
      }
    }
  }
  ierr = print_partitioned_problem(&*p_miit,verb); CHKERRQ(ierr);
  ierr = debug_partitioned_problem(&*p_miit,verb); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
PetscErrorCode Core::block_problem(const string &name,const vector<string> block_problems,int verb) {
  PetscFunctionBegin;
  if(verb==-1) verb = verbose;
  if(!(*build_MoFEM&(1<<0))) SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCT,"fields not build");
  if(!(*build_MoFEM&(1<<1))) SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCT,"FEs not build");
  if(!(*build_MoFEM&(1<<2))) SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCT,"entFEAdjacencies not build");
  if(!(*build_MoFEM&(1<<3))) SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCT,"moFEMProblems not build");
  typedef MoFEMProblem_multiIndex::index<Problem_mi_tag>::type MoFEMProblem_multiIndex_by_name;
  //find pit 
  MoFEMProblem_multiIndex_by_name &moFEMProblems_set = moFEMProblems.get<Problem_mi_tag>();
  MoFEMProblem_multiIndex_by_name::iterator pit = moFEMProblems_set.find(name);
  if(pit==moFEMProblems_set.end()) {
    SETERRQ1(PETSC_COMM_SELF,MOFEM_NOT_FOUND,"problem with name < %s > not defined (top tip check spelling)",name.c_str());
  }
  Vec shift_rows_and_cols_global;
  int ghost[] = {0,1};
  if(rAnk == 0) {
    ierr = ::VecCreateGhost(comm,2,2,2,&ghost[0],&shift_rows_and_cols_global); CHKERRQ(ierr);
  } else {
    ierr = ::VecCreateGhost(comm,0,2,2,&ghost[0],&shift_rows_and_cols_global); CHKERRQ(ierr);
  }
  NumeredDofMoFEMEntity_multiIndex *pit_rows_and_cols[] = {
    const_cast<NumeredDofMoFEMEntity_multiIndex*>(&pit->numered_dofs_rows), 
    const_cast<NumeredDofMoFEMEntity_multiIndex*>(&pit->numered_dofs_cols) };
  int* pit_nb_dofs[] = { pit->tag_nbdof_data_row, pit->tag_nbdof_data_col };
  int* pit_nb_local_dofs[] = { pit->tag_local_nbdof_data_row, pit->tag_local_nbdof_data_col };
  for(int ss = 0;ss<2;ss++) {
    (*pit_nb_dofs[ss]) = 0;
    (*pit_nb_local_dofs[ss]) = 0;
  }
  for(vector<string>::const_iterator vit = block_problems.begin();vit!=block_problems.end();vit++) {
    MoFEMProblem_multiIndex_by_name::iterator piit = moFEMProblems_set.find(*vit);
    if(piit==moFEMProblems_set.end()) {
      SETERRQ1(PETSC_COMM_SELF,MOFEM_NOT_FOUND,"blocked problem with name < %s > not defined (top tip check spelling)",vit->c_str());
    }
    const NumeredDofMoFEMEntity_multiIndex *piit_rows_and_cols[] = {
      &piit->numered_dofs_rows, &piit->numered_dofs_cols };
    ierr = VecZeroEntries(shift_rows_and_cols_global); CHKERRQ(ierr);
    ierr = VecAssemblyBegin(shift_rows_and_cols_global); CHKERRQ(ierr);
    ierr = VecAssemblyEnd(shift_rows_and_cols_global); CHKERRQ(ierr);
    ierr = VecGhostUpdateBegin(shift_rows_and_cols_global,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
    ierr = VecGhostUpdateEnd(shift_rows_and_cols_global,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
    //cerr << piit->get_nb_dofs_row() << " WW " << piit->get_nb_dofs_col() << endl;
    //double *array;
    //ierr = VecGetArray(shift_rows_and_cols_global,&array); CHKERRQ(ierr);
    for(int ss = 0;ss<2;ss++) {
      NumeredDofMoFEMEntity_multiIndex::iterator dit,hi_dit;
      dit = piit_rows_and_cols[ss]->begin();
      hi_dit = piit_rows_and_cols[ss]->end();
      int shift = (*pit_nb_dofs[ss]);
      for(;dit!=hi_dit;dit++) {
	int part_number = dit->get_part();
	if(part_number == rAnk) {
	  ierr = VecSetValue(shift_rows_and_cols_global,ss,1,ADD_VALUES); CHKERRQ(ierr);
	}
	pair<NumeredDofMoFEMEntity_multiIndex::iterator,bool> p;
	p = pit_rows_and_cols[ss]->insert(NumeredDofMoFEMEntity(dit->get_DofMoFEMEntity_ptr())); 
	bool success;
	int dof_idx = shift+dit->get_dof_idx();
	int petsc_global_dof = shift+dit->get_petsc_gloabl_dof_idx();
	success = pit_rows_and_cols[ss]->modify(p.first,NumeredDofMoFEMEntity_mofem_index_change(dof_idx));
	if(!success) SETERRQ(PETSC_COMM_SELF,MOFEM_OPERATION_UNSUCCESSFUL,"modification unsuccessful");
	success = pit_rows_and_cols[ss]->modify(p.first,NumeredDofMoFEMEntity_part_change(part_number,petsc_global_dof));
	if(!success) SETERRQ(PETSC_COMM_SELF,MOFEM_OPERATION_UNSUCCESSFUL,"modification unsuccessful");
    	if(part_number == rAnk) {
	  success = pit_rows_and_cols[ss]->modify(p.first,NumeredDofMoFEMEntity_local_idx_change((*pit_nb_local_dofs[ss])++));
	  if(!success) SETERRQ(PETSC_COMM_SELF,MOFEM_OPERATION_UNSUCCESSFUL,"modification unsuccessful");
	}
      }
    }
    //ierr = VecRestoreArray(shift_rows_and_cols_global,&array); CHKERRQ(ierr);
    ierr = VecAssemblyBegin(shift_rows_and_cols_global); CHKERRQ(ierr);
    ierr = VecAssemblyEnd(shift_rows_and_cols_global); CHKERRQ(ierr);
    ierr = VecGhostUpdateBegin(shift_rows_and_cols_global,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
    ierr = VecGhostUpdateEnd(shift_rows_and_cols_global,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
    ierr = VecGhostUpdateBegin(shift_rows_and_cols_global,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
    ierr = VecGhostUpdateEnd(shift_rows_and_cols_global,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
    double *array;
    ierr = VecGetArray(shift_rows_and_cols_global,&array); CHKERRQ(ierr);
    for(int ss = 0;ss<2;ss++) {
      (*pit_nb_dofs[ss]) += array[ss];
    }
    ierr = VecRestoreArray(shift_rows_and_cols_global,&array); CHKERRQ(ierr);
  }
  ierr = VecDestroy(&shift_rows_and_cols_global); CHKERRQ(ierr);
  ierr = print_partitioned_problem(&*pit,verb); CHKERRQ(ierr);
  ierr = debug_partitioned_problem(&*pit,verb); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
PetscErrorCode Core::print_partitioned_problem(const MoFEMProblem *problem_ptr,int verb) {
  PetscFunctionBegin;
  if(verb==-1) verb = verbose;
  if(verbose>0) {
    ostringstream ss;
    ss << "partition_problem: rank = " << rAnk << " FEs row ghost dofs "<< *problem_ptr 
      << " Nb. local dof " << problem_ptr->get_nb_local_dofs_row() << " nb global row dofs " << problem_ptr->get_nb_dofs_row() << endl;
    ss << "partition_problem: rank = " << rAnk << " FEs col ghost dofs " << *problem_ptr 
      << " Nb. local dof " << problem_ptr->get_nb_local_dofs_col() << " nb global col dofs " << problem_ptr->get_nb_dofs_col() << endl;
    PetscSynchronizedPrintf(comm,ss.str().c_str());
    PetscSynchronizedFlush(comm,PETSC_STDOUT); 
  }
  if(verb>1) {
    ostringstream ss;
    ss << "rank = " << rAnk << " FEs row dofs "<< *problem_ptr << " Nb. row dof " << problem_ptr->get_nb_dofs_row() 
	<< " Nb. local dof " << problem_ptr->get_nb_local_dofs_row() << endl;
    NumeredDofMoFEMEntity_multiIndex::iterator miit_dd_row = problem_ptr->numered_dofs_rows.begin();
    for(;miit_dd_row!=problem_ptr->numered_dofs_rows.end();miit_dd_row++) {
	ss<<*miit_dd_row<<endl;
    }
    ss << "rank = " << rAnk << " FEs col dofs "<< *problem_ptr << " Nb. col dof " << problem_ptr->get_nb_dofs_col() 
	<< " Nb. local dof " << problem_ptr->get_nb_local_dofs_col() << endl;
    NumeredDofMoFEMEntity_multiIndex::iterator miit_dd_col = problem_ptr->numered_dofs_cols.begin();
    for(;miit_dd_col!=problem_ptr->numered_dofs_cols.end();miit_dd_col++) {
	ss<<*miit_dd_col<<endl;
    }
    PetscSynchronizedPrintf(comm,ss.str().c_str());
    PetscSynchronizedFlush(comm,PETSC_STDOUT); 
  }
  PetscFunctionReturn(0);
}
PetscErrorCode Core::debug_partitioned_problem(const MoFEMProblem *problem_ptr,int verb) {
  PetscFunctionBegin;
  if(verb==-1) verb = verbose;
  if(debug>0) {
    typedef NumeredDofMoFEMEntity_multiIndex::index<Idx_mi_tag>::type NumeredDofMoFEMEntitys_by_idx;
    NumeredDofMoFEMEntitys_by_idx::iterator dit,hi_dit;
    dit = problem_ptr->numered_dofs_rows.get<Idx_mi_tag>().begin();
    hi_dit = problem_ptr->numered_dofs_rows.get<Idx_mi_tag>().end();
    for(;dit!=hi_dit;dit++) {
      if(dit->get_part()==(unsigned int)rAnk) {
	if(dit->get_petsc_local_dof_idx()<0) {
	  ostringstream ss;
	  ss << "rank " << rAnk << " " << *dit;
	  SETERRQ1(PETSC_COMM_SELF,1,"local dof index for row not set\n %s",ss.str().c_str());
	}
      }
    }
    dit = problem_ptr->numered_dofs_cols.get<Idx_mi_tag>().begin();
    hi_dit = problem_ptr->numered_dofs_cols.get<Idx_mi_tag>().end();
    for(;dit!=hi_dit;dit++) {
      if(dit->get_part()==(unsigned int)rAnk) {
	if(dit->get_petsc_local_dof_idx()<0) {
	  ostringstream ss;
	  ss << "rank " << rAnk << " " << *dit;
	  SETERRQ1(PETSC_COMM_SELF,1,"local dof index for col not set\n %s",ss.str().c_str());
	}
      }
    }
  }
  PetscFunctionReturn(0);
}

}
