/** 
 * \brief MoFEM interface 
 * 
 * Low level data structures not used directly by user
 *
 * Copyright (C) 2013, Lukasz Kaczmarczyk (likask AT wp.pl) <br>
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

struct CoreTemplates: public Core {

  CoreTemplates(Interface& moab,MPI_Comm _comm = PETSC_COMM_WORLD,int _verbose = 1): 
    Core(moab,_comm,_verbose) {};

  template<typename Tag> 
  PetscErrorCode create_Mat(
    const string &name,Mat *M,const MatType type,PetscInt **_i,PetscInt **_j,PetscScalar **_v,const bool no_diagonals = true,int verb = -1);

};

template<typename Tag> 
PetscErrorCode CoreTemplates::create_Mat(
  const string &name,Mat *M,const MatType type,PetscInt **_i,PetscInt **_j,PetscScalar **_v,
  const bool no_diagonals,int verb) {
  PetscFunctionBegin;
  PetscLogEventBegin(USER_EVENT_createMat,0,0,0,0);
  if(verb==-1) verb = verbose;
  typedef typename boost::multi_index::index<NumeredDofMoFEMEntity_multiIndex,Tag>::type NumeredDofMoFEMEntitys_by_idx;
  typedef MoFEMEntityEntMoFEMFiniteElementAdjacencyMap_multiIndex::index<Unique_mi_tag>::type adj_by_ent;
  //find p_miit
  typedef MoFEMProblem_multiIndex::index<Problem_mi_tag>::type moFEMProblems_by_name;
  moFEMProblems_by_name &moFEMProblems_set = moFEMProblems.get<Problem_mi_tag>();
  moFEMProblems_by_name::iterator p_miit = moFEMProblems_set.find(name);
  if(p_miit==moFEMProblems_set.end()) SETERRQ1(PETSC_COMM_SELF,MOFEM_NOT_FOUND,"problem < %s > is not found (top tip: check spelling)",name.c_str());
  //
  const NumeredDofMoFEMEntitys_by_idx &dofs_row_by_idx = p_miit->numered_dofs_rows.get<Tag>();
  const NumeredDofMoFEMEntitys_by_idx &dofs_col_by_idx = p_miit->numered_dofs_cols.get<Tag>();
  DofIdx nb_dofs_row = p_miit->get_nb_dofs_row();
  if(nb_dofs_row == 0) {
    SETERRQ1(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCT,"problem <%s> has zero rows",name.c_str());
  }
  map<int,vector<int> > adjacent_dofs_on_other_parts;
  typename boost::multi_index::index<NumeredDofMoFEMEntity_multiIndex,Tag>::type::iterator miit_row,hi_miit_row;
  if(Tag::IamNotPartitioned) {
    //get range of local indices
    PetscLayout layout;
    ierr = PetscLayoutCreate(comm,&layout); CHKERRQ(ierr);
    ierr = PetscLayoutSetBlockSize(layout,1); CHKERRQ(ierr);
    ierr = PetscLayoutSetSize(layout,nb_dofs_row); CHKERRQ(ierr);
    ierr = PetscLayoutSetUp(layout); CHKERRQ(ierr);
    PetscInt rstart,rend;
    ierr = PetscLayoutGetRange(layout,&rstart,&rend); CHKERRQ(ierr);
    ierr = PetscLayoutDestroy(&layout); CHKERRQ(ierr);
    if(verb > 0) {
	PetscSynchronizedPrintf(comm,"\tcreate_Mat: row lower %d row upper %d\n",rstart,rend);
	//PetscSynchronizedFlush(comm,PETSC_STDOUT); 
    }
    miit_row = dofs_row_by_idx.lower_bound(rstart);
    hi_miit_row = dofs_row_by_idx.lower_bound(rend);
    if(distance(miit_row,hi_miit_row) != rend-rstart) {
	SETERRQ4(PETSC_COMM_SELF,PETSC_ERR_ARG_SIZ,
	  "data inconsistency, distance(miit_row,hi_miit_row) != rend - rstart (%d != %d - %d = %d) ",
	  distance(miit_row,hi_miit_row),rend,rstart,rend-rstart);
    }
  } else {
    miit_row = dofs_row_by_idx.lower_bound(rAnk);
    hi_miit_row = dofs_row_by_idx.upper_bound(rAnk);
    //get adjacent nodes on other partitions
    vector<int> dofs_vec;
    MoFEMEntity *mofem_ent_ptr = NULL;
    NumeredDofMoFEMEntity_multiIndex_uid_view_hashed dofs_col_view;
    typename boost::multi_index::index<NumeredDofMoFEMEntity_multiIndex,Tag>::type::iterator mit_row,hi_mit_row;
    mit_row = dofs_row_by_idx.begin();
    hi_mit_row = dofs_row_by_idx.end();
    for(;mit_row!=hi_mit_row;mit_row++) {
	bitset<8> pstatus(mit_row->get_pstatus());
	if(pstatus.test(0)) {
	  if(pstatus.test(1)||pstatus.test(2)) {
	    if( (mofem_ent_ptr == NULL) ? 1 : (mofem_ent_ptr->get_global_unique_id() != mit_row->get_MoFEMEntity_ptr()->get_global_unique_id()) ) {
	      mofem_ent_ptr = const_cast<MoFEMEntity*>(mit_row->get_MoFEMEntity_ptr());
	      adj_by_ent::iterator adj_miit = entFEAdjacencies.get<Unique_mi_tag>().lower_bound(mofem_ent_ptr->get_global_unique_id());
	      adj_by_ent::iterator hi_adj_miit = entFEAdjacencies.get<Unique_mi_tag>().upper_bound(mofem_ent_ptr->get_global_unique_id());
	      dofs_col_view.clear();
	      for(;adj_miit!=hi_adj_miit;adj_miit++) {
		if(adj_miit->by_other&BYROW) {
		  if((adj_miit->EntMoFEMFiniteElement_ptr->get_id()&p_miit->get_BitFEId()).none()) {
		    // if element is not part of problem
		    continue; 
		  }
		  if((adj_miit->EntMoFEMFiniteElement_ptr->get_BitRefLevel()&mit_row->get_BitRefLevel()).none()) {
		    // if entity is not problem refinment level
		    continue; 
		  }
		  ierr = adj_miit->EntMoFEMFiniteElement_ptr->get_MoFEMFiniteElement_col_dof_view( 
		    p_miit->numered_dofs_cols,dofs_col_view,Interface::UNION); CHKERRQ(ierr);
		}
	      }
	      dofs_vec.push_back(Tag::get_index(mit_row));
	      dofs_vec.push_back(dofs_col_view.size());
	      NumeredDofMoFEMEntity_multiIndex_uid_view_hashed::iterator cvit;
	      cvit = dofs_col_view.begin();
	      for(;cvit!=dofs_col_view.end();cvit++) {
		int idx = Tag::get_index(*cvit);
		dofs_vec.push_back(idx);
		if(idx<0) {
		  SETERRQ(PETSC_COMM_SELF,PETSC_ERR_ARG_SIZ,"data inconsistency");
		}
		if(idx>=p_miit->get_nb_dofs_col()) {
		  SETERRQ(PETSC_COMM_SELF,PETSC_ERR_ARG_SIZ,"data inconsistency");
		}
	      }
	    }
	  }
	}
    }
    //gather information from other processors
    IS is,isout;
    ierr = ISCreateGeneral(comm,
	dofs_vec.size(),&*dofs_vec.begin(),
	PETSC_USE_POINTER,&is); CHKERRQ(ierr);
    ierr = ISAllGather(is,&isout); CHKERRQ(ierr);
    int isout_size;
    ierr = ISGetSize(isout,&isout_size); CHKERRQ(ierr);
    const int *ptr;
    ierr = ISGetIndices(isout,&ptr); CHKERRQ(ierr);
    for(int ii = 0;ii<isout_size;) {
	int row_idx = ptr[ii++];
	int nb_adj_dofs = ptr[ii++];
	for(int jj = 0;jj<nb_adj_dofs;jj++) {
	  adjacent_dofs_on_other_parts[row_idx].push_back(ptr[ii++]);
	}
    }
    ierr = ISRestoreIndices(isout,&ptr); CHKERRQ(ierr);
    ierr = ISDestroy(&is); CHKERRQ(ierr);
    ierr = ISDestroy(&isout); CHKERRQ(ierr);
  }
  int nb_loc_row_from_iterators = distance(miit_row,hi_miit_row);
  MoFEMEntity *mofem_ent_ptr = NULL;
  int row_last_evaluated_idx = -1;
  vector<PetscInt> i,j;
  vector<DofIdx> dofs_vec;
  NumeredDofMoFEMEntity_multiIndex_uid_view_hashed dofs_col_view;
  // loop local rows
  unsigned int rows_to_fill = distance(miit_row,hi_miit_row);
  i.reserve( rows_to_fill+1 );
  for(;miit_row!=hi_miit_row;miit_row++) {
    i.push_back(j.size());
    if(strcmp(type,MATMPIADJ)==0) {
	DofIdx idx = Tag::get_index(miit_row);
	if(dofs_col_by_idx.find(idx)->get_global_unique_id()!=miit_row->get_global_unique_id()) {
	  SETERRQ(PETSC_COMM_SELF,PETSC_ERR_ARG_SIZ,"data insonsistency");
	}
    }
    if( (mofem_ent_ptr == NULL) ? 1 : (mofem_ent_ptr->get_global_unique_id() != miit_row->get_MoFEMEntity_ptr()->get_global_unique_id()) ) {
	// get field ptr
	mofem_ent_ptr = const_cast<MoFEMEntity*>(miit_row->get_MoFEMEntity_ptr());
	row_last_evaluated_idx = Tag::get_index(miit_row);
	adj_by_ent::iterator adj_miit = entFEAdjacencies.get<Unique_mi_tag>().lower_bound(mofem_ent_ptr->get_global_unique_id());
	adj_by_ent::iterator hi_adj_miit = entFEAdjacencies.get<Unique_mi_tag>().upper_bound(mofem_ent_ptr->get_global_unique_id());
	dofs_col_view.clear();
	for(;adj_miit!=hi_adj_miit;adj_miit++) {
	  if(adj_miit->by_other&BYROW) {
	    if((adj_miit->EntMoFEMFiniteElement_ptr->get_id()&p_miit->get_BitFEId()).none()) {
	      // if element is not part of problem
	      continue; 
	    }
	    if((adj_miit->EntMoFEMFiniteElement_ptr->get_BitRefLevel()&miit_row->get_BitRefLevel()).none()) {
	      // if entity is not problem refinment level
	      continue; 
	    }
	    ierr = adj_miit->EntMoFEMFiniteElement_ptr->get_MoFEMFiniteElement_col_dof_view( 
	      p_miit->numered_dofs_cols,dofs_col_view,Interface::UNION); CHKERRQ(ierr);
	  }
	}
	dofs_vec.resize(0);
	NumeredDofMoFEMEntity_multiIndex_uid_view_hashed::iterator cvit;
	cvit = dofs_col_view.begin();
	for(;cvit!=dofs_col_view.end();cvit++) {
	  int idx = Tag::get_index(*cvit);
	  dofs_vec.push_back(idx);
	  if(idx<0) {
	    SETERRQ(PETSC_COMM_SELF,PETSC_ERR_ARG_SIZ,"data inconsistency");
	  }
	  if(idx>=p_miit->get_nb_dofs_col()) {
	    SETERRQ(PETSC_COMM_SELF,PETSC_ERR_ARG_SIZ,"data inconsistency");
	  }
	}
	bitset<8> pstatus(miit_row->get_pstatus());
	if(!pstatus.test(0)) {
	  if(pstatus.test(1)||pstatus.test(2)) {
	    map<int,vector<int> >::iterator mit;
	    mit = adjacent_dofs_on_other_parts.find(row_last_evaluated_idx);
	    if(mit == adjacent_dofs_on_other_parts.end()) {
	      SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCT,"data inconsistency");
	    }
	    dofs_vec.insert(dofs_vec.end(),mit->second.begin(),mit->second.end());
	  }
	}
	sort(dofs_vec.begin(),dofs_vec.end());
	if(!pstatus.test(0)) {
	  if(pstatus.test(1)||pstatus.test(2)) {
	    unique(dofs_vec.begin(),dofs_vec.end());
	  }
	}
    }
    //if(dofs_vec.size()==0) {
	//SETERRQ1(PETSC_COMM_SELF,PETSC_ERR_ARG_SIZ,"zero dofs at row %d",Tag::get_index(miit_row));
    //}
    if( j.capacity() < j.size() + dofs_vec.size() ) {
	unsigned int nb_nonzero = j.size() + dofs_vec.size();
	unsigned int average_row_fill = nb_nonzero/i.size() + nb_nonzero % i.size();
	if( j.capacity() < rows_to_fill*average_row_fill ) {
	  j.reserve( rows_to_fill*average_row_fill );
	}
    }
    vector<DofIdx>::iterator diit,hi_diit;
    diit = dofs_vec.begin();
    hi_diit = dofs_vec.end();
    for(;diit!=hi_diit;diit++) {
	if(no_diagonals) {
	  if(*diit == Tag::get_index(miit_row)) {
	    continue;
	  }
	}
	j.push_back(*diit);
    }
  }
  //build adj matrix
  i.push_back(j.size());
  ierr = PetscMalloc(i.size()*sizeof(PetscInt),_i); CHKERRQ(ierr);
  ierr = PetscMalloc(j.size()*sizeof(PetscInt),_j); CHKERRQ(ierr);
  copy(i.begin(),i.end(),*_i);
  copy(j.begin(),j.end(),*_j);
  PetscInt nb_row_dofs = p_miit->get_nb_dofs_row();
  PetscInt nb_col_dofs = p_miit->get_nb_dofs_col();
  if(strcmp(type,MATMPIADJ)==0) { 
    if(i.size()-1 != (unsigned int)nb_loc_row_from_iterators) {
	SETERRQ(PETSC_COMM_SELF,PETSC_ERR_ARG_SIZ,"data inconsistency");
    }
    ierr = MatCreateMPIAdj(comm,i.size()-1,nb_col_dofs,*_i,*_j,PETSC_NULL,M); CHKERRQ(ierr);
    ierr = MatSetOption(*M,MAT_STRUCTURALLY_SYMMETRIC,PETSC_TRUE); CHKERRQ(ierr);
  } else if(strcmp(type,MATMPIAIJ)==0) {
    if(i.size()-1 != (unsigned int)nb_loc_row_from_iterators) {
	SETERRQ(PETSC_COMM_SELF,PETSC_ERR_ARG_SIZ,"data inconsistency");
    }
    PetscInt nb_local_dofs_row = p_miit->get_nb_local_dofs_row();
    if((unsigned int)nb_local_dofs_row!=i.size()-1) {
	SETERRQ(PETSC_COMM_SELF,PETSC_ERR_ARG_SIZ,"data inconsistency");
    }
    PetscInt nb_local_dofs_col = p_miit->get_nb_local_dofs_col();
    ierr = ::MatCreateMPIAIJWithArrays(comm,nb_local_dofs_row,nb_local_dofs_col,nb_row_dofs,nb_col_dofs,*_i,*_j,PETSC_NULL,M); CHKERRQ(ierr);
  } else if(strcmp(type,MATAIJ)==0) {
    if(i.size()-1 != (unsigned int)nb_loc_row_from_iterators) {
	SETERRQ(PETSC_COMM_SELF,PETSC_ERR_ARG_SIZ,"data inconsistency");
    }
    PetscInt nb_local_dofs_row = p_miit->get_nb_local_dofs_row();
    if((unsigned int)nb_local_dofs_row!=i.size()-1) {
	SETERRQ(PETSC_COMM_SELF,PETSC_ERR_ARG_SIZ,"data inconsistency");
    }
    PetscInt nb_local_dofs_col = p_miit->get_nb_local_dofs_col();
    ierr = ::MatCreateSeqAIJWithArrays(PETSC_COMM_SELF,nb_local_dofs_row,nb_local_dofs_col,*_i,*_j,PETSC_NULL,M); CHKERRQ(ierr);
  } else {
    SETERRQ(PETSC_COMM_SELF,PETSC_ERR_ARG_NULL,"not implemented");
  }

  //MatView(*M,PETSC_VIEWER_STDOUT_WORLD);


  PetscLogEventEnd(USER_EVENT_createMat,0,0,0,0);
  PetscFunctionReturn(0);
}

PetscErrorCode Core::MatCreateMPIAIJWithArrays(const string &name,Mat *Aij,int verb) {
  PetscFunctionBegin;
  if(verb==-1) verb = verbose;
  int *_i,*_j;
  CoreTemplates *core_ptr = static_cast<CoreTemplates*>(const_cast<Core*>(this));
  ierr = core_ptr->create_Mat<Part_mi_tag>(name,Aij,MATMPIAIJ,&_i,&_j,PETSC_NULL,false,verb); CHKERRQ(ierr);
  ierr = PetscFree(_i); CHKERRQ(ierr);
  ierr = PetscFree(_j); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
PetscErrorCode Core::MatCreateSeqAIJWithArrays(const string &name,Mat *Aij,PetscInt **i,PetscInt **j,PetscScalar **v,int verb) {
  PetscFunctionBegin;
  if(verb==-1) verb = verbose;
  CoreTemplates *core_ptr = static_cast<CoreTemplates*>(const_cast<Core*>(this));
  ierr = core_ptr->create_Mat<PetscLocalIdx_mi_tag>(name,Aij,MATAIJ,i,j,v,false,verb); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode Core::partition_problem(const string &name,int verb) {
  PetscFunctionBegin;
  if(verb==-1) verb = verbose;
  if(!(*build_MoFEM&(1<<0))) SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCT,"fields not build");
  if(!(*build_MoFEM&(1<<1))) SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCT,"FEs not build");
  if(!(*build_MoFEM&(1<<2))) SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCT,"entFEAdjacencies not build");
  if(!(*build_MoFEM&(1<<3))) SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCT,"moFEMProblems not build");
  if(verb>0) {
    PetscPrintf(comm,"Partition problem %s\n",name.c_str());
  }
  typedef NumeredDofMoFEMEntity_multiIndex::index<Idx_mi_tag>::type NumeredDofMoFEMEntitys_by_idx;
  typedef MoFEMProblem_multiIndex::index<Problem_mi_tag>::type moFEMProblems_by_name;
  //find p_miit
  moFEMProblems_by_name &moFEMProblems_set = moFEMProblems.get<Problem_mi_tag>();
  moFEMProblems_by_name::iterator p_miit = moFEMProblems_set.find(name);
  if(p_miit==moFEMProblems_set.end()) SETERRQ1(PETSC_COMM_SELF,1,"problem with name %s not defined (top tip check spelling)",name.c_str());
  DofIdx nb_dofs_row = p_miit->get_nb_dofs_row();
  int *i,*j;
  Mat Adj;
  if(verb>1) {
    PetscPrintf(comm,"\tcreate Adj matrix\n");
  }
  try {
    CoreTemplates *core_ptr = static_cast<CoreTemplates*>(const_cast<Core*>(this));
    ierr = core_ptr->create_Mat<Idx_mi_tag>(name,&Adj,MATMPIADJ,&i,&j,PETSC_NULL,true,verb); CHKERRQ(ierr);
  } catch (const char* msg) {
    SETERRQ(PETSC_COMM_SELF,MOFEM_CHAR_THROW,msg);
  } catch (const std::exception& ex) {
    ostringstream ss;
    ss << "throw in method: " << ex.what() << " at line " << __LINE__ << " in file " << __FILE__ << endl;
    SETERRQ(PETSC_COMM_SELF,MOFEM_STD_EXCEPTION_THROW,ss.str().c_str());
  }
  if(verb>1) {
    PetscPrintf(comm,"\t<- done\n");
  }
  int m,n;
  ierr = MatGetSize(Adj,&m,&n); CHKERRQ(ierr);
  if(verb>2) {
    MatView(Adj,PETSC_VIEWER_STDOUT_WORLD);
  }
  //partitioning
  MatPartitioning part;
  IS is;
  ierr = MatPartitioningCreate(comm,&part); CHKERRQ(ierr);
  ierr = MatPartitioningSetAdjacency(part,Adj); CHKERRQ(ierr);
  ierr = MatPartitioningSetFromOptions(part); CHKERRQ(ierr);
  ierr = MatPartitioningSetNParts(part,sIze); CHKERRQ(ierr);
  ierr = MatPartitioningApply(part,&is); CHKERRQ(ierr);
  if(verb>2) {
    ISView(is,PETSC_VIEWER_STDOUT_WORLD);
  }
  //gather
  IS is_gather,is_num,is_gather_num;
  ierr = ISAllGather(is,&is_gather); CHKERRQ(ierr);
  ierr = ISPartitioningToNumbering(is,&is_num); CHKERRQ(ierr);
  ierr = ISAllGather(is_num,&is_gather_num); CHKERRQ(ierr);
  const int *part_number,*petsc_idx;
  ierr = ISGetIndices(is_gather,&part_number);  CHKERRQ(ierr);
  ierr = ISGetIndices(is_gather_num,&petsc_idx);  CHKERRQ(ierr);
  int size_is_num,size_is_gather;
  ISGetSize(is_gather,&size_is_gather);
  if(size_is_gather != (int)nb_dofs_row) {
    SETERRQ2(PETSC_COMM_SELF,PETSC_ERR_ARG_SIZ,"data inconsistency %d != %d",size_is_gather,nb_dofs_row);
  }
  ISGetSize(is_num,&size_is_num);
  if(size_is_num != (int)nb_dofs_row) {
    SETERRQ2(PETSC_COMM_SELF,PETSC_ERR_ARG_SIZ,"data inconsistency %d != %d",size_is_num,nb_dofs_row);
  }
  //set petsc global indicies
  NumeredDofMoFEMEntitys_by_idx &dofs_row_by_idx_no_const = const_cast<NumeredDofMoFEMEntitys_by_idx&>(p_miit->numered_dofs_rows.get<Idx_mi_tag>());
  NumeredDofMoFEMEntitys_by_idx &dofs_col_by_idx_no_const = const_cast<NumeredDofMoFEMEntitys_by_idx&>(p_miit->numered_dofs_cols.get<Idx_mi_tag>());
  DofIdx &nb_row_local_dofs = *((DofIdx*)p_miit->tag_local_nbdof_data_row);
  DofIdx &nb_col_local_dofs = *((DofIdx*)p_miit->tag_local_nbdof_data_col);
  nb_row_local_dofs = 0;
  nb_col_local_dofs = 0;
  DofIdx &nb_row_ghost_dofs = *((DofIdx*)p_miit->tag_ghost_nbdof_data_row);
  DofIdx &nb_col_ghost_dofs = *((DofIdx*)p_miit->tag_ghost_nbdof_data_col);
  nb_row_ghost_dofs = 0;
  nb_col_ghost_dofs = 0;
  NumeredDofMoFEMEntitys_by_idx::iterator miit_dofs_row = dofs_row_by_idx_no_const.begin();
  NumeredDofMoFEMEntitys_by_idx::iterator miit_dofs_col = dofs_col_by_idx_no_const.begin();
  if(verb>1) {
    PetscPrintf(comm,"\tloop problem dofs");
  }
  try {
  for(;miit_dofs_row!=dofs_row_by_idx_no_const.end();miit_dofs_row++,miit_dofs_col++) {
    if(miit_dofs_col==dofs_col_by_idx_no_const.end()) {
      SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCT,"check finite element definition, nb. of rows is not equal to number for columns");
    }
    if(miit_dofs_row->get_global_unique_id()!=miit_dofs_col->get_global_unique_id()) {
      SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCT,"check finite element definition, nb. of rows is not equal to columns");
    }
    if(miit_dofs_row->dof_idx!=miit_dofs_col->dof_idx) {
      SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCT,"check finite element definition, nb. of rows is not equal to columns");
    }
    assert(petsc_idx[miit_dofs_row->dof_idx]>=0);
    assert(petsc_idx[miit_dofs_row->dof_idx]<(int)p_miit->get_nb_dofs_row());
    bool success = dofs_row_by_idx_no_const.modify(miit_dofs_row,NumeredDofMoFEMEntity_part_change(part_number[miit_dofs_row->dof_idx],petsc_idx[miit_dofs_row->dof_idx]));
    if(!success) {
      SETERRQ(PETSC_COMM_SELF,MOFEM_OPERATION_UNSUCCESSFUL,"modification unsuccessful");
    }
    success = dofs_col_by_idx_no_const.modify(miit_dofs_col,NumeredDofMoFEMEntity_part_change(part_number[miit_dofs_col->dof_idx],petsc_idx[miit_dofs_col->dof_idx]));
    if(!success) {
      SETERRQ(PETSC_COMM_SELF,MOFEM_OPERATION_UNSUCCESSFUL,"modification unsuccessful");
    }
    if(miit_dofs_row->part == (unsigned int)rAnk) {
      assert(miit_dofs_row->part==miit_dofs_col->part);
      assert(miit_dofs_row->petsc_gloabl_dof_idx==miit_dofs_col->petsc_gloabl_dof_idx);
      success = dofs_row_by_idx_no_const.modify(miit_dofs_row,NumeredDofMoFEMEntity_local_idx_change(nb_row_local_dofs++));
      if(!success) {
	SETERRQ(PETSC_COMM_SELF,MOFEM_OPERATION_UNSUCCESSFUL,"modification unsuccessful");
      }
      success = dofs_col_by_idx_no_const.modify(miit_dofs_col,NumeredDofMoFEMEntity_local_idx_change(nb_col_local_dofs++));
      if(!success) {
	SETERRQ(PETSC_COMM_SELF,MOFEM_OPERATION_UNSUCCESSFUL,"modification unsuccessful");
      }
    }
  }
  } catch (const char* msg) {
    SETERRQ(PETSC_COMM_SELF,MOFEM_CHAR_THROW,msg);
  } catch (const std::exception& ex) {
    ostringstream ss;
    ss << "throw in method: " << ex.what() << " at line " << __LINE__ << " in file " << __FILE__ << endl;
    SETERRQ(PETSC_COMM_SELF,MOFEM_STD_EXCEPTION_THROW,ss.str().c_str());
  }
  if(verb>1) {
    PetscPrintf(comm," <- done\n");
  }
  if(verbose>0) {
    ostringstream ss;
    ss << "partition_problem: rank = " << rAnk << " FEs row ghost dofs "<< *p_miit 
	<< " Nb. local dof " << p_miit->get_nb_local_dofs_row() << " nb global row dofs " << p_miit->get_nb_dofs_row() << endl;
    ss << "partition_problem: rank = " << rAnk << " FEs col ghost dofs " << *p_miit 
	<< " Nb. local dof " << p_miit->get_nb_local_dofs_col() << " nb global col dofs " << p_miit->get_nb_dofs_col() << endl;
    PetscSynchronizedPrintf(comm,ss.str().c_str());
    PetscSynchronizedFlush(comm,PETSC_STDOUT); 
  }
  if(verb>2) {
    ostringstream ss;
    ss << "rank = " << rAnk << " FEs row dofs "<< *p_miit << " Nb. row dof " << p_miit->get_nb_dofs_row() 
	<< " Nb. local dof " << p_miit->get_nb_local_dofs_row() << endl;
    NumeredDofMoFEMEntity_multiIndex::iterator miit_dd_row = p_miit->numered_dofs_rows.begin();
    for(;miit_dd_row!=p_miit->numered_dofs_rows.end();miit_dd_row++) {
	ss<<*miit_dd_row<<endl;
    }
    ss << "rank = " << rAnk << " FEs col dofs "<< *p_miit << " Nb. col dof " << p_miit->get_nb_dofs_col() 
	<< " Nb. local dof " << p_miit->get_nb_local_dofs_col() << endl;
    NumeredDofMoFEMEntity_multiIndex::iterator miit_dd_col = p_miit->numered_dofs_cols.begin();
    for(;miit_dd_col!=p_miit->numered_dofs_cols.end();miit_dd_col++) {
	ss<<*miit_dd_col<<endl;
    }
    PetscSynchronizedPrintf(comm,ss.str().c_str());
    PetscSynchronizedFlush(comm,PETSC_STDOUT); 
  }
  ierr = ISRestoreIndices(is_gather,&part_number);  CHKERRQ(ierr);
  ierr = ISRestoreIndices(is_gather_num,&petsc_idx);  CHKERRQ(ierr);
  ierr = ISDestroy(&is_num); CHKERRQ(ierr);
  ierr = ISDestroy(&is_gather_num); CHKERRQ(ierr);
  ierr = ISDestroy(&is_gather); CHKERRQ(ierr);
  ierr = ISDestroy(&is); CHKERRQ(ierr);
  ierr = MatPartitioningDestroy(&part); CHKERRQ(ierr);
  ierr = MatDestroy(&Adj); CHKERRQ(ierr);
  if(debug>0) {
    try {
      NumeredDofMoFEMEntitys_by_idx::iterator dit,hi_dit;
      dit = p_miit->numered_dofs_rows.get<Idx_mi_tag>().begin();
      hi_dit = p_miit->numered_dofs_rows.get<Idx_mi_tag>().end();
      for(;dit!=hi_dit;dit++) {
	if(dit->get_part()==(unsigned int)rAnk) {
	  if(dit->get_petsc_local_dof_idx()<0) {
	    ostringstream ss;
	    ss << "rank " << rAnk << " " << *dit;
	    SETERRQ1(PETSC_COMM_SELF,1,"local dof index for row not set\n %s",ss.str().c_str());
	  }
	}
      }
      dit = p_miit->numered_dofs_cols.get<Idx_mi_tag>().begin();
      hi_dit = p_miit->numered_dofs_cols.get<Idx_mi_tag>().end();
      for(;dit!=hi_dit;dit++) {
	if(dit->get_part()==(unsigned int)rAnk) {
	  if(dit->get_petsc_local_dof_idx()<0) {
	    ostringstream ss;
	    ss << "rank " << rAnk << " " << *dit;
	    SETERRQ1(PETSC_COMM_SELF,1,"local dof index for col not set\n %s",ss.str().c_str());
	  }
	}
      }
    } catch (const char* msg) {
      SETERRQ(PETSC_COMM_SELF,MOFEM_CHAR_THROW,msg);
    } catch (const std::exception& ex) {
      ostringstream ss;
      ss << "throw in method: " << ex.what() << " at line " << __LINE__ << " in file " << __FILE__ << endl;
      SETERRQ(PETSC_COMM_SELF,MOFEM_STD_EXCEPTION_THROW,ss.str().c_str());
    }
  }
  *build_MoFEM |= 1<<4;
  PetscFunctionReturn(0);
}

PetscErrorCode Core::partition_check_matrix_fill_in(const string &problem_name,int verb) {
  PetscFunctionBegin;
  if(verb==-1) verb = verbose;

  struct TestMatrixFillIn:public FEMethod {
    FieldInterface *mFieldPtr;

    Mat A;
    PetscErrorCode ierr;
    ErrorCode rval;

    TestMatrixFillIn(FieldInterface *m_field_ptr,Mat _A): mFieldPtr(m_field_ptr),A(_A) {};
    PetscErrorCode preProcess() {
      PetscFunctionBegin;
      PetscFunctionReturn(0);
    }

    PetscErrorCode operator()() {
      PetscFunctionBegin;
      if(refinedFiniteElementsPtr->find(fePtr->get_ent())==refinedFiniteElementsPtr->end()) {
	SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCT,"data inconsistency");
      }
      FENumeredDofMoFEMEntity_multiIndex::iterator rit = rowPtr->begin();
      for(;rit!=rowPtr->end();rit++) {
	if(refinedEntitiesPtr->find(rit->get_ent())==refinedEntitiesPtr->end()) {
	  SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCT,"data inconsistency");
	}
	if(!rit->get_active()) {
	  SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCT,"data inconsistency");
	}
	MoFEMEntityEntMoFEMFiniteElementAdjacencyMap_multiIndex::index<Composite_unique_mi_tag>::type::iterator ait;
	ait = adjacenciesPtr->get<Composite_unique_mi_tag>().find(boost::make_tuple(
	      rit->get_MoFEMEntity_ptr()->get_global_unique_id(),fePtr->get_global_unique_id()));
	if(ait==adjacenciesPtr->end()) {
	    ostringstream ss;
	    ss << *rit << endl;
	    ss << *fePtr << endl;
	    ss << "dof: " << rit->get_BitRefLevel() << endl;
	    ss << "fe: " << fePtr->get_BitRefLevel() << endl;
	    ss << "problem: " << problemPtr->get_BitRefLevel() << endl;
	    PetscPrintf(mFieldPtr->get_comm(),"%s",ss.str().c_str());
	    SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCT,"adjacencies data inconsistency");
	} else {
	  LocalUId uid = ait->get_ent_unique_id();
	  if(entitiesPtr->find(uid) == entitiesPtr->end()) {
	    SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCT,"data inconsistency");
	  } 
	  if(dofsPtr->find(rit->get_global_unique_id())==dofsPtr->end()) {
	    SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCT,"data inconsistency");
	  }
	}
	int row = rit->get_petsc_gloabl_dof_idx();
	FENumeredDofMoFEMEntity_multiIndex::iterator cit = colPtr->begin();
	for(;cit!=colPtr->end();cit++) {
	  if(refinedEntitiesPtr->find(cit->get_ent())==refinedEntitiesPtr->end()) {
	    SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCT,"data inconsistency");
	  }
	  if(!cit->get_active()) {
	    SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCT,"data inconsistency");
	  }
	  int col = cit->get_petsc_gloabl_dof_idx();
	  ait = adjacenciesPtr->get<Composite_unique_mi_tag>().find(boost::make_tuple(
	      cit->get_MoFEMEntity_ptr()->get_global_unique_id(),fePtr->get_global_unique_id()));
	  if(ait==adjacenciesPtr->end()) {
	    SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCT,"adjacencies data inconsistency");
	  } else {
	    LocalUId uid = ait->get_ent_unique_id();
	    if(entitiesPtr->find(uid) == entitiesPtr->end()) {
	      SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCT,"data inconsistency");
	    } 
	    if(dofsPtr->find(cit->get_global_unique_id())==dofsPtr->end()) {
	      SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCT,"data inconsistency");
	    }
	  }
  	  ierr = MatSetValue(A,row,col,1,INSERT_VALUES);

	  if(ierr!=0) {
	  //if(row == 87 && col == 909) {
	  
	    EntityHandle ent = fePtr->get_ent();
      
	    ostringstream ss;
	    ss << "fe:\n" << *fePtr << endl;
	    ss << "row:\n" << *rit << endl;
	    ss << "col:\n" << *cit << endl;

	    ss << "fe:\n" << fePtr->get_BitRefLevel() << endl;
	    ss << "row:\n" << rit->get_BitRefLevel() << endl;
	    ss << "col:\n" << cit->get_BitRefLevel() << endl;

	    ss << "edges:\n";
	    for(int ee = 0;ee<6;ee++) {
	      EntityHandle edge;
	      rval = mFieldPtr->get_moab().side_element(ent,1,ee,edge); CHKERR_THROW(rval);
	      ss << edge << " ";
	    }
	    ss << endl;

	    PetscPrintf(mFieldPtr->get_comm(),"%s\n",ss.str().c_str());
	  //}
	  }
	}
      }
      PetscFunctionReturn(0);
    }

    PetscErrorCode postProcess() {
      PetscFunctionBegin;
      PetscErrorCode ierr;

      ierr = MatAssemblyBegin(A,MAT_FLUSH_ASSEMBLY); CHKERRQ(ierr);
      ierr = MatAssemblyEnd(A,MAT_FLUSH_ASSEMBLY); CHKERRQ(ierr);

      PetscFunctionReturn(0);
    }


  };

  Mat A;
  ierr = MatCreateMPIAIJWithArrays(problem_name,&A); CHKERRQ(ierr);
  ierr = MatSetOption(A,MAT_NEW_NONZERO_ALLOCATION_ERR,PETSC_TRUE);  CHKERRQ(ierr);
  TestMatrixFillIn method(this,A);

  typedef MoFEMProblem_multiIndex::index<Problem_mi_tag>::type moFEMProblems_by_name;
  //find p_miit
  moFEMProblems_by_name &moFEMProblems_set = moFEMProblems.get<Problem_mi_tag>();
  moFEMProblems_by_name::iterator p_miit = moFEMProblems_set.find(problem_name);
  if(p_miit == moFEMProblems_set.end()) {
    SETERRQ1(PETSC_COMM_SELF,1,"problem < %s > not found (top tip: check spelling)",problem_name.c_str());
  }
  if(verb>0) {
    PetscPrintf(comm,"check problem < %s >\n",problem_name.c_str());
  }
  //MoFEMFiniteElement set
  MoFEMFiniteElement_multiIndex::iterator fe = finiteElements.begin();
  MoFEMFiniteElement_multiIndex::iterator hi_fe = finiteElements.end();
  for(;fe!=hi_fe;fe++) {
    if(verb>0) {
      PetscPrintf(comm,"\tcheck element %s\n",fe->get_name().c_str());
    }

    ierr = loop_finite_elements(problem_name,fe->get_name(),method,0,sIze,verb);  CHKERRQ(ierr);

  }

  ierr = MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
  ierr = MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);

  ierr = MatDestroy(&A); CHKERRQ(ierr);
 
  PetscFunctionReturn(0);
}


}
