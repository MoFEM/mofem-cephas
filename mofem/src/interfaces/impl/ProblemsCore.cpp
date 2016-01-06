/** \file FieldInterfaceCore.cpp
 * \brief Mylti-index containers, data structures and other low-level functions
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

#include <MaterialBlocks.hpp>
#include <CubitBCData.hpp>
#include <TagMultiIndices.hpp>
#include <CoordSysMultiIndices.hpp>
#include <FieldMultiIndices.hpp>
#include <EntsMultiIndices.hpp>
#include <DofsMultiIndices.hpp>
#include <FEMMultiIndices.hpp>
#include <ProblemsMultiIndices.hpp>
#include <AdjacencyMultiIndices.hpp>
#include <BCMultiIndices.hpp>
#include <CoreDataStructures.hpp>
#include <SeriesMultiIndices.hpp>

#include <LoopMethods.hpp>
#include <FieldInterface.hpp>
#include <MeshRefinment.hpp>
#include <PrismInterface.hpp>
#include <SeriesRecorder.hpp>
#include <Core.hpp>

namespace MoFEM {

#ifndef __INTEL_COMPILER
  #define SIZEOFUID 16
#else
  #define SIZEOFUID 24
#endif

struct __attribute__ ((__packed__)) IdxDataType {
  int globalDof;
  char uId[SIZEOFUID];
  IdxDataType() {}
  IdxDataType(const GlobalUId &uid,int global_dof):
    globalDof(global_dof) {
    bcopy(&uid,uId,SIZEOFUID);
  }
};


const static int debug = 1;

PetscErrorCode Core::build_partitioned_problem(const string &name,bool square_matrix,int verb) {
  PetscFunctionBegin;
  if(verb==-1) verb = verbose;
  if(!(*build_MoFEM&(1<<0))) SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"fields not build");
  if(!(*build_MoFEM&(1<<1))) SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"FEs not build");
  if(!(*build_MoFEM&(1<<2))) SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"entFEAdjacencies not build");
  const MoFEMProblem *problem_ptr;
  ierr = get_problem(name,&problem_ptr); CHKERRQ(ierr);
  ierr = build_partitioned_problem(const_cast<MoFEMProblem*>(problem_ptr),square_matrix,verb); CHKERRQ(ierr);
  *build_MoFEM |= 1<<3;
  *build_MoFEM |= 1<<4; // It is assumed that user who uses this function knows what he is doing
  PetscFunctionReturn(0);
}
PetscErrorCode Core::build_partitioned_problem(MoFEMProblem *problem_ptr,bool square_matrix,int verb) {
  PetscFunctionBegin;
  if(verb==-1) verb = verbose;
  if(problem_ptr->get_BitRefLevel().none()) {
    SETERRQ1(PETSC_COMM_SELF,1,"problem <%s> refinement level not set",problem_ptr->get_name().c_str());
  }
  //zero finite elements
  ProblemClearNumeredFiniteElementsChange().operator()(*problem_ptr);
  //bool success = pRoblems.modify(p_miit,ProblemClearNumeredFiniteElementsChange());
  //if(!success) SETERRQ(PETSC_COMM_SELF,MOFEM_OPERATION_UNSUCCESSFUL,"modification unsuccessful");

  //get rows and cols dofs view
  DofMoFEMEntity_multiIndex_active_view dofs_rows,dofs_cols;
  {
    //fe_miit iterator for finite elements
    EntMoFEMFiniteElement_multiIndex::iterator fe_miit = entsFiniteElements.begin();
    EntMoFEMFiniteElement_multiIndex::iterator hi_fe_miit = entsFiniteElements.end();
    //iterate all finite elemen entities in database
    for(;fe_miit!=hi_fe_miit;fe_miit++) {
      //if element is in problem
      if((fe_miit->get_id()&problem_ptr->get_BitFEId()).any()) {
        //if finite element bit level has all refined bits sets
        if((fe_miit->get_BitRefLevel()&problem_ptr->get_BitRefLevel())==problem_ptr->get_BitRefLevel()) {
          //get dof uids for rows and columns
          ierr = fe_miit->get_MoFEMFiniteElement_row_dof_view(dofsField,dofs_rows); CHKERRQ(ierr);
          ierr = fe_miit->get_MoFEMFiniteElement_col_dof_view(dofsField,dofs_cols); CHKERRQ(ierr);
        }
      }
    }
  }

  //get problem bit level
  const BitRefLevel &problem_bit_level = problem_ptr->get_DofMask_BitRefLevel();

  //add dofs for rows and cols and set ownership
  DofMoFEMEntity_multiIndex_active_view* dofs_ptr[] = { &dofs_rows, &dofs_cols };
  NumeredDofMoFEMEntity_multiIndex *numered_dofs_ptr[] = {
    &(problem_ptr->numered_dofs_rows), &(problem_ptr->numered_dofs_cols) };
  int* nbdof_ptr[] = {
    problem_ptr->tag_nbdof_data_row, problem_ptr->tag_nbdof_data_col };
  int* local_nbdof_ptr[] = {
    problem_ptr->tag_local_nbdof_data_row, problem_ptr->tag_local_nbdof_data_col };
  int* ghost_nbdof_ptr[] = {
    problem_ptr->tag_ghost_nbdof_data_row, problem_ptr->tag_ghost_nbdof_data_col };
  int nb_local_dofs[] = { 0,0 };
  for(int ss = 0;ss<2;ss++) {
    *(nbdof_ptr[ss]) = 0;
    *(local_nbdof_ptr[ss]) = 0;
    *(ghost_nbdof_ptr[ss]) = 0;
    DofMoFEMEntity_multiIndex_active_view::nth_index<1>::type::iterator miit,hi_miit;
    miit = dofs_ptr[ss]->get<1>().lower_bound(1);
    hi_miit = dofs_ptr[ss]->get<1>().upper_bound(1);
    for(;miit!=hi_miit;miit++) {
      const BitRefLevel &dof_bit_level = (*miit)->get_BitRefLevel();
      if((dof_bit_level&problem_bit_level)!=dof_bit_level) {
        continue;
      }
      NumeredDofMoFEMEntity dof(&**miit);
      pair<NumeredDofMoFEMEntity_multiIndex::iterator,bool> p = numered_dofs_ptr[ss]->insert(dof);
      int owner_proc = p.first->get_owner_proc();
      if(owner_proc<0) {
        SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"data inconsistency");
      }
      if(owner_proc == rAnk) {
        nb_local_dofs[ss]++;
      }
      bool success;
      success = numered_dofs_ptr[ss]->modify(p.first,NumeredDofMoFEMEntity_part_change(owner_proc,-1));
      if(!success) SETERRQ(PETSC_COMM_SELF,MOFEM_OPERATION_UNSUCCESSFUL,"modification unsuccessful");
      success = numered_dofs_ptr[ss]->modify(p.first,NumeredDofMoFEMEntity_mofem_index_change(-1));
      if(!success) SETERRQ(PETSC_COMM_SELF,MOFEM_OPERATION_UNSUCCESSFUL,"modification unsuccessful");
      success = numered_dofs_ptr[ss]->modify(p.first,NumeredDofMoFEMEntity_local_idx_change(-1));
      if(!success) SETERRQ(PETSC_COMM_SELF,MOFEM_OPERATION_UNSUCCESSFUL,"modification unsuccessful");
    }
  }

  //get layout
  int start_ranges[2],end_ranges[2];
  for(int ss = 0;ss<2;ss++) {
    PetscLayout layout;
    ierr = PetscLayoutCreate(comm,&layout); CHKERRQ(ierr);
    ierr = PetscLayoutSetBlockSize(layout,1); CHKERRQ(ierr);
    ierr = PetscLayoutSetLocalSize(layout,nb_local_dofs[ss]); CHKERRQ(ierr);
    ierr = PetscLayoutSetUp(layout); CHKERRQ(ierr);
    ierr = PetscLayoutGetSize(layout,&*nbdof_ptr[ss]); CHKERRQ(ierr); // get global size
    ierr = PetscLayoutGetRange(layout,&start_ranges[ss],&end_ranges[ss]); CHKERRQ(ierr); //get ranges
    ierr = PetscLayoutDestroy(&layout); CHKERRQ(ierr);
  }

  if(sizeof(UId) != SIZEOFUID) {
    SETERRQ1(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"check size of UId, size of UId is %u",sizeof(UId));
  }

  //set local and global indices on own dofs
  const size_t idx_data_type_size = sizeof(IdxDataType);
  const size_t data_block_size = idx_data_type_size/sizeof(int);
  if(sizeof(IdxDataType) % sizeof(int)) {
    SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"data inconsistency");
  }
  vector<vector<IdxDataType> > ids_data_packed_rows(sIze),ids_data_packed_cols(sIze);

  //ParallelComm* pcomm = ParallelComm::get_pcomm(&moab,MYPCOMM_INDEX);
  for(int ss = 0;ss<2;ss++) {

    NumeredDofMoFEMEntity_multiIndex::index<Part_mi_tag>::type::iterator mit,hi_mit;
    mit = numered_dofs_ptr[ss]->get<Part_mi_tag>().lower_bound(rAnk);
    hi_mit = numered_dofs_ptr[ss]->get<Part_mi_tag>().upper_bound(rAnk);
    int &local_idx = *local_nbdof_ptr[ss];
    for(;mit!=hi_mit;mit++) {

      bool success;
      success = numered_dofs_ptr[ss]->modify( numered_dofs_ptr[ss]->project<0>(mit),NumeredDofMoFEMEntity_local_idx_change(local_idx));
      if(!success) SETERRQ(PETSC_COMM_SELF,MOFEM_OPERATION_UNSUCCESSFUL,"modification unsuccessful");
      int glob_idx = start_ranges[ss]+local_idx;
      success = numered_dofs_ptr[ss]->modify( numered_dofs_ptr[ss]->project<0>(mit),NumeredDofMoFEMEntity_mofem_index_change(glob_idx));
      if(!success) SETERRQ(PETSC_COMM_SELF,MOFEM_OPERATION_UNSUCCESSFUL,"modification unsuccessful");
      success = numered_dofs_ptr[ss]->modify( numered_dofs_ptr[ss]->project<0>(mit),NumeredDofMoFEMEntity_part_change(mit->get_part(),glob_idx));
      if(!success) SETERRQ(PETSC_COMM_SELF,MOFEM_OPERATION_UNSUCCESSFUL,"modification unsuccessful");
      local_idx++;
      unsigned char pstatus = mit->get_pstatus();

      //check id dof is shared
      if(pstatus>0) {

        for(int proc = 0; proc<MAX_SHARING_PROCS && -1 != mit->get_sharing_procs_ptr()[proc]; proc++) {
          if(ss == 0) {
            ids_data_packed_rows[mit->get_sharing_procs_ptr()[proc]].push_back(IdxDataType(mit->get_global_unique_id(),glob_idx));
          } else {
            ids_data_packed_cols[mit->get_sharing_procs_ptr()[proc]].push_back(IdxDataType(mit->get_global_unique_id(),glob_idx));
          }
          if(!(pstatus&PSTATUS_MULTISHARED)) {
            break;
          }
        }

      }
    }

  }

  int nsends_rows = 0,nsends_cols = 0;
  // Non zero lengths[i] represent a message to i of length lengths[i].
  vector<int> lengths_rows(sIze),lengths_cols(sIze);
  lengths_rows.clear();
  lengths_cols.clear();
  for(int proc = 0;proc<sIze;proc++) {
    lengths_rows[proc] = ids_data_packed_rows[proc].size()*data_block_size;
    lengths_cols[proc] = ids_data_packed_cols[proc].size()*data_block_size;
    if(!ids_data_packed_rows[proc].empty()) nsends_rows++;
    if(!ids_data_packed_cols[proc].empty()) nsends_cols++;
  }

  MPI_Status *status;
  ierr = PetscMalloc1(sIze,&status);CHKERRQ(ierr);

  // make sure it is a PETSc comm
  ierr = PetscCommDuplicate(comm,&comm,NULL); CHKERRQ(ierr);

  //rows

  // Computes the number of messages a node expects to receive
  int nrecvs_rows;	// number of messages received
  ierr = PetscGatherNumberOfMessages(comm,NULL,&lengths_rows[0],&nrecvs_rows); CHKERRQ(ierr);
  //cerr << nrecvs_rows << endl;

  // Computes info about messages that a MPI-node will receive, including (from-id,length) pairs for each message.
  int *onodes_rows;	// list of node-ids from which messages are expected
  int *olengths_rows;	// corresponding message lengths
  ierr = PetscGatherMessageLengths(comm,
    nsends_rows,nrecvs_rows,
    &lengths_rows[0],&onodes_rows,&olengths_rows);  CHKERRQ(ierr);

  // Gets a unique new tag from a PETSc communicator. All processors that share
  // the communicator MUST call this routine EXACTLY the same number of times.
  // This tag should only be used with the current objects communicator; do NOT
  // use it with any other MPI communicator.
  int tag_row;
  ierr = PetscCommGetNewTag(comm,&tag_row); CHKERRQ(ierr);

  // Allocate a buffer sufficient to hold messages of size specified in
  // olengths. And post Irecvs on these buffers using node info from onodes
  int **rbuf_row;	// must bee freed by user
  MPI_Request *r_waits_row; // must bee freed by user
  // rbuf has a pointers to messeges. It has size of of nrecvs (number of
  // messages) +1. In the first index a block is allocated,
  // such that rbuf[i] = rbuf[i-1]+olengths[i-1].
  ierr = PetscPostIrecvInt(comm,tag_row,nrecvs_rows,onodes_rows,olengths_rows,&rbuf_row,&r_waits_row); CHKERRQ(ierr);
  ierr = PetscFree(onodes_rows); CHKERRQ(ierr);


  MPI_Request *s_waits_row; // status of sens messages
  ierr = PetscMalloc1(nsends_rows,&s_waits_row);CHKERRQ(ierr);

  // Send messeges
  for(int proc=0,kk=0; proc<sIze; proc++) {
    if(!lengths_rows[proc]) continue; 	// no message to send to this proc
    ierr = MPI_Isend(
      &(ids_data_packed_rows[proc])[0], // buffer to send
      lengths_rows[proc], 		// message length
      MPIU_INT,proc, 			// to proc
      tag_row,comm,s_waits_row+kk); CHKERRQ(ierr);
    kk++;
  }

  if(nrecvs_rows) {
    ierr = MPI_Waitall(nrecvs_rows,r_waits_row,status);CHKERRQ(ierr);
  }
  if(nsends_rows) {
    ierr = MPI_Waitall(nsends_rows,s_waits_row,status);CHKERRQ(ierr);
  }

  ierr = PetscFree(r_waits_row); CHKERRQ(ierr);
  ierr = PetscFree(s_waits_row); CHKERRQ(ierr);

  //cols
  int nrecvs_cols = nrecvs_rows;
  int *olengths_cols = olengths_rows;
  PetscInt **rbuf_col = rbuf_row;
  if(!square_matrix) {

    // Computes the number of messages a node expects to receive
    ierr = PetscGatherNumberOfMessages(comm,NULL,&lengths_cols[0],&nrecvs_cols); CHKERRQ(ierr);

    // Computes info about messages that a MPI-node will receive, including (from-id,length) pairs for each message.
    int *onodes_cols;
    ierr = PetscGatherMessageLengths(comm,
      nsends_cols,nrecvs_cols,
      &lengths_cols[0],&onodes_cols,&olengths_cols);  CHKERRQ(ierr);

    // Gets a unique new tag from a PETSc communicator.
    int tag_col;
    ierr = PetscCommGetNewTag(comm,&tag_col); CHKERRQ(ierr);

    // Allocate a buffer sufficient to hold messages of size specified in
    // olengths. And post Irecvs on these buffers using node info from onodes
    MPI_Request *r_waits_col; // must bee freed by user
    ierr = PetscPostIrecvInt(comm,tag_col,nrecvs_cols,onodes_cols,olengths_cols,&rbuf_col,&r_waits_col); CHKERRQ(ierr);
    ierr = PetscFree(onodes_cols); CHKERRQ(ierr);

    MPI_Request *s_waits_col; // status of sens messages
    ierr = PetscMalloc1(nsends_cols,&s_waits_col);CHKERRQ(ierr);

    // Send messeges
    for(int proc=0,kk=0; proc<sIze; proc++) {
      if(!lengths_cols[proc]) continue; 	// no message to send to this proc
      ierr = MPI_Isend(
        &(ids_data_packed_cols[proc])[0],	// buffer to send
        lengths_cols[proc], 			// message length
        MPIU_INT,proc, 				// to proc
        tag_col,comm,s_waits_col+kk); CHKERRQ(ierr);
        kk++;
      }

      if(nrecvs_cols) {
      ierr = MPI_Waitall(nrecvs_cols,r_waits_col,status);CHKERRQ(ierr);
    }
    if(nsends_cols) {
      ierr = MPI_Waitall(nsends_cols,s_waits_col,status);CHKERRQ(ierr);
    }

    ierr = PetscFree(r_waits_col); CHKERRQ(ierr);
    ierr = PetscFree(s_waits_col); CHKERRQ(ierr);

  }

  // set values
  for(int ss = 0;ss<2;ss++) {

    int nrecvs;
    int *olengths;
    int **data_procs;
    if(ss == 0) {
      nrecvs = nrecvs_rows;
      olengths = olengths_rows;
      data_procs = rbuf_row;
    } else {
      nrecvs = nrecvs_cols;
      olengths = olengths_cols;
      data_procs = rbuf_col;
    }

    IdxDataType *idx_data;
    GlobalUId uid;

    NumeredDofMoFEMEntity_multiIndex::iterator dit;
    for(int kk=0; kk<nrecvs; kk++) {

      int len = olengths[kk];
      int *data_from_proc = data_procs[kk];
      for(int dd = 0;dd<len;dd+=data_block_size) {

        idx_data = (IdxDataType*)(&data_from_proc[dd]);
        bcopy(idx_data->uId,&uid,sizeof(UId));
        dit = numered_dofs_ptr[ss]->find(uid);
        if(dit == numered_dofs_ptr[ss]->end()) {
          DofMoFEMEntity_multiIndex::iterator ddit = dofsField.find(uid);
          if(ddit!=dofsField.end()) {
            cerr << *ddit << endl;
          } else {
            ostringstream zz;
            zz << uid << endl;
            SETERRQ1(PETSC_COMM_SELF,MOFEM_OPERATION_UNSUCCESSFUL,"no such dof %s in mofem database",zz.str().c_str());
          }
          ostringstream zz;
          zz << uid << endl;
          SETERRQ1(PETSC_COMM_SELF,MOFEM_OPERATION_UNSUCCESSFUL,"dof %s not found",zz.str().c_str());
        }
        int global_idx = idx_data->globalDof;
        bool success;
        success = numered_dofs_ptr[ss]->modify(dit,NumeredDofMoFEMEntity_mofem_index_change(global_idx));
        if(!success) SETERRQ(PETSC_COMM_SELF,MOFEM_OPERATION_UNSUCCESSFUL,"modification unsuccessful");
        success = numered_dofs_ptr[ss]->modify(dit,NumeredDofMoFEMEntity_part_change(dit->get_part(),global_idx));
        if(!success) SETERRQ(PETSC_COMM_SELF,MOFEM_OPERATION_UNSUCCESSFUL,"modification unsuccessful");

      }

    }

  }

  ierr = PetscFree(olengths_rows); CHKERRQ(ierr);
  ierr = PetscFree(rbuf_row[0]); CHKERRQ(ierr);
  ierr = PetscFree(rbuf_row); CHKERRQ(ierr);
  if(!square_matrix) {
    ierr = PetscFree(olengths_cols); CHKERRQ(ierr);
    ierr = PetscFree(rbuf_col[0]); CHKERRQ(ierr);
    ierr = PetscFree(rbuf_col); CHKERRQ(ierr);
  }

  ierr = PetscFree(status); CHKERRQ(ierr);

  ierr = printPartitionedProblem(problem_ptr,verb); CHKERRQ(ierr);
  ierr = debugPartitionedProblem(problem_ptr,verb); CHKERRQ(ierr);

  *build_MoFEM |= 1<<3;
  *build_MoFEM |= 1<<4; // It is assumed that user who uses this function knows what he is doing
  PetscFunctionReturn(0);
}
PetscErrorCode Core::build_partitioned_problems(int verb) {
  PetscFunctionBegin;
  if(verb==-1) verb = verbose;
  DofMoFEMEntity_multiIndex_active_view dofs_rows,dofs_cols;
  MoFEMProblem_multiIndex::iterator p_miit = pRoblems.begin();
  for(;p_miit!=pRoblems.end();p_miit++) {
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
  EntMoFEMFiniteElement_multiIndex::iterator miit2 = entsFiniteElements.begin();
  EntMoFEMFiniteElement_multiIndex::iterator hi_miit2 = entsFiniteElements.end();
  DofMoFEMEntity_multiIndex_active_view dofs_rows,dofs_cols;
  EntMoFEMFiniteElement_multiIndex::iterator miit3 = miit2;
  //iterate all finite elemen entities in database
  for(;miit3!=hi_miit2;miit3++) {
    //if element is in problem
    if((miit3->get_id()&problem_ptr->get_BitFEId()).any()) {
      //if finite element bit level has all refined bits sets
      if((miit3->get_BitRefLevel()&problem_ptr->get_BitRefLevel())==problem_ptr->get_BitRefLevel()) {
        //get dof uids for rows and columns
        ierr = miit3->get_MoFEMFiniteElement_row_dof_view(dofsField,dofs_rows); CHKERRQ(ierr);
        ierr = miit3->get_MoFEMFiniteElement_col_dof_view(dofsField,dofs_cols); CHKERRQ(ierr);
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
    ProblemAddRowDof(&**miit4).operator()(*problem_ptr);
  }

  //add dofs for cols
  DofMoFEMEntity_multiIndex_active_view::nth_index<1>::type::iterator miit5,hi_miit5;
  miit5 = dofs_cols.get<1>().lower_bound(1);
  hi_miit5 = dofs_cols.get<1>().upper_bound(1);
  for(;miit5!=hi_miit5;miit5++) {
    if(((*miit5)->get_BitRefLevel()&problem_ptr->get_DofMask_BitRefLevel())!=(*miit5)->get_BitRefLevel()) {
      continue;
    }
    ProblemAddColDof(&**miit5).operator()(*problem_ptr);
  }

  //number dofs on rows and columns
  ProblemRowNumberChange().operator()(*problem_ptr);
  ProblemColNumberChange().operator()(*problem_ptr);

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
  if(!(*build_MoFEM&(1<<0))) SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"fields not build");
  if(!(*build_MoFEM&(1<<1))) SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"FEs not build");
  if(!(*build_MoFEM&(1<<2))) SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"entFEAdjacencies not build");
  const MoFEMProblem *problem_ptr;
  ierr = get_problem(problem_name,&problem_ptr); CHKERRQ(ierr);
  ierr = build_problem(const_cast<MoFEMProblem*>(problem_ptr),verb); CHKERRQ(ierr);
  *build_MoFEM |= 1<<3; // It is assumed that user who uses this function knows what he is doing
  PetscFunctionReturn(0);
}
PetscErrorCode Core::clear_problem(const string &problem_name,int verb) {
  PetscFunctionBegin;
  if(verb==-1) verb = verbose;
  typedef MoFEMProblem_multiIndex::index<Problem_mi_tag>::type pRoblemsByName;
  pRoblemsByName &prob_by_name = pRoblems.get<Problem_mi_tag>();
  pRoblemsByName::iterator p_miit = prob_by_name.find(problem_name);
  if(p_miit == prob_by_name.end()) {
    SETERRQ1(
      PETSC_COMM_SELF,MOFEM_OPERATION_UNSUCCESSFUL,
      "problem < %s > not found, (top tip: check spelling)",problem_name.c_str()
    );
  }
  //zero rows
  bool success = prob_by_name.modify(p_miit,ProblemZeroNbRowsChange());
  if(!success) SETERRQ(PETSC_COMM_SELF,MOFEM_OPERATION_UNSUCCESSFUL,"modification unsuccessful");
  //zero cols
  success = prob_by_name.modify(p_miit,ProblemZeroNbColsChange());
  if(!success) SETERRQ(PETSC_COMM_SELF,MOFEM_OPERATION_UNSUCCESSFUL,"modification unsuccessful");
  //clear finite elements
  success = prob_by_name.modify(p_miit,ProblemClearNumeredFiniteElementsChange());
  if(!success) SETERRQ(PETSC_COMM_SELF,MOFEM_OPERATION_UNSUCCESSFUL,"modification unsuccessful");
  PetscFunctionReturn(0);
}
PetscErrorCode Core::build_problems(int verb) {
  PetscFunctionBegin;
  if(verb==-1) verb = verbose;
  if(!(*build_MoFEM&(1<<0))) SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"fields not build");
  if(!(*build_MoFEM&(1<<1))) SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"FEs not build");
  if(!(*build_MoFEM&(1<<2))) SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"entFEAdjacencies not build");
  //iterate problems
  MoFEMProblem_multiIndex::iterator p_miit = pRoblems.begin();
  for(;p_miit!=pRoblems.end();p_miit++) {
    MoFEMProblem *problem_ptr =  const_cast<MoFEMProblem*>(&*p_miit);
    ierr = build_problem(problem_ptr,verb); CHKERRQ(ierr);
  }
  *build_MoFEM |= 1<<3;
  PetscFunctionReturn(0);
}
PetscErrorCode Core::clear_problems(int verb) {
  PetscFunctionBegin;
  if(verb==-1) verb = verbose;
  MoFEMProblem_multiIndex::iterator p_miit = pRoblems.begin();
  //iterate problems
  for(;p_miit!=pRoblems.end();p_miit++) {
    //zero rows
    bool success = pRoblems.modify(p_miit,ProblemZeroNbRowsChange());
    if(!success) SETERRQ(PETSC_COMM_SELF,MOFEM_OPERATION_UNSUCCESSFUL,"modification unsuccessful");
    //zero cols
    success = pRoblems.modify(p_miit,ProblemZeroNbColsChange());
    if(!success) SETERRQ(PETSC_COMM_SELF,MOFEM_OPERATION_UNSUCCESSFUL,"modification unsuccessful");
    //clear finite elements
    success = pRoblems.modify(p_miit,ProblemClearNumeredFiniteElementsChange());
    if(!success) SETERRQ(PETSC_COMM_SELF,MOFEM_OPERATION_UNSUCCESSFUL,"modification unsuccessful");
  }
  PetscFunctionReturn(0);
}
PetscErrorCode Core::partition_simple_problem(const string &name,int verb) {
  PetscFunctionBegin;
  if(verb==-1) verb = verbose;
  if(!(*build_MoFEM&(1<<0))) SETERRQ(PETSC_COMM_SELF,MOFEM_NOT_FOUND,"fields not build");
  if(!(*build_MoFEM&(1<<1))) SETERRQ(PETSC_COMM_SELF,MOFEM_NOT_FOUND,"FEs not build");
  if(!(*build_MoFEM&(1<<2))) SETERRQ(PETSC_COMM_SELF,MOFEM_NOT_FOUND,"entFEAdjacencies not build");
  if(!(*build_MoFEM&(1<<3))) SETERRQ(PETSC_COMM_SELF,MOFEM_NOT_FOUND,"pRoblems not build");
  if(verb>0) {
    PetscPrintf(comm,"Simple partition problem %s\n",name.c_str());
  }
  // find p_miit
  typedef MoFEMProblem_multiIndex::index<Problem_mi_tag>::type MoFEMProblem_multiIndex_by_name;
  MoFEMProblem_multiIndex_by_name &pRoblems_set = pRoblems.get<Problem_mi_tag>();
  MoFEMProblem_multiIndex_by_name::iterator p_miit = pRoblems_set.find(name);
  if(p_miit==pRoblems_set.end()) SETERRQ1(PETSC_COMM_SELF,1,"problem < %s > is not found (top tip: check spelling)",name.c_str());
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
          distance(miit_col,hi_miit_col),ranges_col[part+1],ranges_col[part],ranges_col[part+1]-ranges_col[part]
        );
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
    ierr = printPartitionedProblem(&*p_miit,verb); CHKERRQ(ierr);
    *build_MoFEM |= 1<<4;
    PetscFunctionReturn(0);
}
PetscErrorCode Core::partition_compose_problem(const string &name,const string &problem_for_rows,bool copy_rows,const string &problem_for_cols,bool copy_cols,int verb) {
  PetscFunctionBegin;
  if(verb==-1) verb = verbose;
  if(!(*build_MoFEM&(1<<0))) SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"fields not build");
  if(!(*build_MoFEM&(1<<1))) SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"FEs not build");
  if(!(*build_MoFEM&(1<<2))) SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"entFEAdjacencies not build");
  if(!(*build_MoFEM&(1<<3))) SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"pRoblems not build");

  typedef MoFEMProblem_multiIndex::index<Problem_mi_tag>::type MoFEMProblem_multiIndex_by_name;
  typedef NumeredDofMoFEMEntity_multiIndex::index<Unique_mi_tag>::type NumeredDofMoFEMEntitys_by_uid;

  //find p_miit
  MoFEMProblem_multiIndex_by_name &pRoblems_set = pRoblems.get<Problem_mi_tag>();
  MoFEMProblem_multiIndex_by_name::iterator p_miit = pRoblems_set.find(name);
  if(p_miit==pRoblems_set.end()) SETERRQ1(PETSC_COMM_SELF,1,"problem with name < %s > not defined (top tip check spelling)",name.c_str());
  if(verb>0) {
    PetscPrintf(
      comm,"Compose problem %s from rows of %s and columns of %s\n",
      p_miit->get_name().c_str(),problem_for_rows.c_str(),problem_for_cols.c_str()
    );
  }

  //find p_miit_row
  MoFEMProblem_multiIndex_by_name::iterator p_miit_row = pRoblems_set.find(problem_for_rows);
  if(p_miit_row==pRoblems_set.end()) SETERRQ1(PETSC_COMM_SELF,1,"problem with name < %s > not defined (top tip check spelling)",problem_for_rows.c_str());
  const NumeredDofMoFEMEntity_multiIndex &dofs_row = p_miit_row->numered_dofs_rows;
  //find p_mit_col
  MoFEMProblem_multiIndex_by_name::iterator p_miit_col = pRoblems_set.find(problem_for_cols);
  if(p_miit_col==pRoblems_set.end()) SETERRQ1(PETSC_COMM_SELF,1,"problem with name < %s > not defined (top tip check spelling)",problem_for_cols.c_str());
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

      // only copy indices which are belong to some elements if this porblem
      vector<int> is_local,is_new;

      NumeredDofMoFEMEntitys_by_uid &dofs_by_uid = const_cast<NumeredDofMoFEMEntitys_by_uid&>(copied_dofs[ss]->get<Unique_mi_tag>());
      for(NumeredDofMoFEMEntity_multiIndex::iterator dit = composed_dofs[ss]->begin();dit!=composed_dofs[ss]->end();dit++) {

        NumeredDofMoFEMEntitys_by_uid::iterator diit = dofs_by_uid.find(dit->get_global_unique_id());
        if(diit==dofs_by_uid.end()) {
          SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"data inconsistency, could not find dof in composite problem");
        }
        int part_number = diit->get_part(); // get part number
        int petsc_global_dof = diit->get_petsc_gloabl_dof_idx();
        bool success;
        success = composed_dofs[ss]->modify(dit,NumeredDofMoFEMEntity_part_change(part_number,petsc_global_dof));
        if(!success) SETERRQ(PETSC_COMM_SELF,MOFEM_OPERATION_UNSUCCESSFUL,"modification unsuccessful");
        if(dit->get_part() == (unsigned int)rAnk) {
          success =composed_dofs[ss]->modify(dit,NumeredDofMoFEMEntity_local_idx_change((*nb_local_dofs[ss])++));
          if(!success) SETERRQ(PETSC_COMM_SELF,MOFEM_OPERATION_UNSUCCESSFUL,"modification unsuccessful");
          is_local.push_back(petsc_global_dof);
        }

      }

      AO ao;
      ierr = AOCreateMapping(comm,is_local.size(),&is_local[0],NULL,&ao); CHKERRQ(ierr);

      // apply local to global mapping
      is_local.resize(0);
      for(NumeredDofMoFEMEntity_multiIndex::iterator dit = composed_dofs[ss]->begin();dit!=composed_dofs[ss]->end();dit++) {

        is_local.push_back(dit->get_petsc_gloabl_dof_idx());

      }
      ierr = AOPetscToApplication(ao,is_local.size(),&is_local[0]); CHKERRQ(ierr);
      int idx2 = 0;
      for(NumeredDofMoFEMEntity_multiIndex::iterator dit = composed_dofs[ss]->begin();dit!=composed_dofs[ss]->end();dit++) {

        int part_number = dit->get_part(); // get part number
        int petsc_global_dof = is_local[idx2++];
        bool success;
        success = composed_dofs[ss]->modify(dit,NumeredDofMoFEMEntity_part_change(part_number,petsc_global_dof));
        if(!success) SETERRQ(PETSC_COMM_SELF,MOFEM_OPERATION_UNSUCCESSFUL,"modification unsuccessful");

      }

      ierr = AODestroy(&ao); CHKERRQ(ierr);

    } else {

      for(NumeredDofMoFEMEntity_multiIndex::iterator dit = copied_dofs[ss]->begin();dit!=copied_dofs[ss]->end();dit++) {
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

  ierr = printPartitionedProblem(&*p_miit,verb); CHKERRQ(ierr);
  ierr = debugPartitionedProblem(&*p_miit,verb); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
PetscErrorCode Core::printPartitionedProblem(const MoFEMProblem *problem_ptr,int verb) {
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
PetscErrorCode Core::debugPartitionedProblem(const MoFEMProblem *problem_ptr,int verb) {
  PetscFunctionBegin;
  if(verb==-1) verb = verbose;
  if(debug>0) {

    typedef NumeredDofMoFEMEntity_multiIndex::index<Idx_mi_tag>::type NumeredDofMoFEMEntitys_by_idx;
    NumeredDofMoFEMEntitys_by_idx::iterator dit,hi_dit;
    const NumeredDofMoFEMEntitys_by_idx* numered_dofs_ptr[] = {
      &(problem_ptr->numered_dofs_rows.get<Idx_mi_tag>()), &(problem_ptr->numered_dofs_rows.get<Idx_mi_tag>())
    };

    int* nbdof_ptr[] = {
      problem_ptr->tag_nbdof_data_row, problem_ptr->tag_nbdof_data_col
    };
    int* local_nbdof_ptr[] = {
      problem_ptr->tag_local_nbdof_data_row, problem_ptr->tag_local_nbdof_data_col
    };

    for(int ss = 0;ss<2;ss++) {

      dit = numered_dofs_ptr[ss]->begin();
      hi_dit = problem_ptr->numered_dofs_rows.get<Idx_mi_tag>().end();
      for(;dit!=hi_dit;dit++) {
        if(dit->get_part()==(unsigned int)rAnk) {
          if(dit->get_petsc_local_dof_idx()<0) {
            ostringstream zz;
            zz << "rank " << rAnk << " " << *dit;
            SETERRQ2(PETSC_COMM_SELF,MOFEM_IMPOSIBLE_CASE,"local dof index for %d (0-row, 1-col) not set, i.e. has negative value\n %s",ss,zz.str().c_str());
          }
          if(dit->get_petsc_local_dof_idx()>=*local_nbdof_ptr[ss]) {
            ostringstream zz;
            zz << "rank " << rAnk << " " << *dit;
            SETERRQ2(PETSC_COMM_SELF,MOFEM_IMPOSIBLE_CASE,"local dofs for %d (0-row, 1-col) out of range\n %s",ss,zz.str().c_str());
          }
        } else {
          if(dit->get_petsc_gloabl_dof_idx()<0) {
            ostringstream zz;
            zz << "rank " << rAnk << " " << *dit;
            SETERRQ2(PETSC_COMM_SELF,MOFEM_IMPOSIBLE_CASE,"global dof index for %d (0-row, 1-col) row not set, i.e. has negative value\n %s",ss,zz.str().c_str());
          }
          if(dit->get_petsc_gloabl_dof_idx()>=*nbdof_ptr[ss]) {
            ostringstream zz;
            zz << "rank " << rAnk << " " << *dit;
            SETERRQ2(PETSC_COMM_SELF,MOFEM_IMPOSIBLE_CASE,"global dofs for %d (0-row, 1-col) out of range\n %s",ss,zz.str().c_str());
          }
        }
      }
    }

  }
  PetscFunctionReturn(0);
}

}
