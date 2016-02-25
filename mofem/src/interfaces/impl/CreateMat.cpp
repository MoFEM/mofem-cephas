/**
 * \brief Create adjacent matrices using different indices

  Low level data structures not used directly by user
  TODO: Make it work to create local sequential matrices using local indices.

 * MoFEM is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
 * License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with MoFEM. If not, see <http://www.gnu.org/licenses/>

*/

#include <Includes.hpp>
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

const static bool debug = false;

/** \brief Create compressed matrix

  NOTE: Only function class members are allowed in this class. NO VARIABLES.

  */
struct CreateRowComressedADJMatrix: public Core {

  CreateRowComressedADJMatrix(Interface& moab,MPI_Comm _comm = PETSC_COMM_WORLD,int _verbose = 1):
    Core(moab,_comm,_verbose) {};

  typedef MoFEMEntityEntMoFEMFiniteElementAdjacencyMap_multiIndex::index<Unique_mi_tag>::type AdjByEnt;
  typedef MoFEMProblem_multiIndex::index<Problem_mi_tag>::type ProblemsByName;
  typedef NumeredDofMoFEMEntity_multiIndex::index<PetscGlobalIdx_mi_tag>::type DofByGlobalPetscIndex;

  /** \brief Create matrix adjacencies

    Depending on TAG type, which some structure used to number dofs, matrix is
    partitioned using part number stored in multi-index, or is partitioned on parts
    based only on index number.

    See: Idx_mi_tag  PetscGlobalIdx_mi_tag and PetscLocalIdx_mi_tag

    */
  template<typename TAG>
  PetscErrorCode createMat(
    const string &name,Mat *M,const MatType type,
    PetscInt **_i,PetscInt **_j,PetscScalar **_v,
    const bool no_diagonals = true,int verb = -1);

  /** \brief Get element adjacencies
    */
  template<typename TAG>
  PetscErrorCode getEntityAdjacenies(
    ProblemsByName::iterator p_miit,
    typename boost::multi_index::index<NumeredDofMoFEMEntity_multiIndex,TAG>::type::iterator mit_row,
    MoFEMEntity* &mofem_ent_ptr,NumeredDofMoFEMEntity_multiIndex_uid_view_hashed &dofs_col_view,int verb);


};

template<typename TAG>
PetscErrorCode CreateRowComressedADJMatrix::getEntityAdjacenies(
  ProblemsByName::iterator p_miit,
  typename boost::multi_index::index<NumeredDofMoFEMEntity_multiIndex,TAG>::type::iterator mit_row,
  MoFEMEntity* &mofem_ent_ptr,NumeredDofMoFEMEntity_multiIndex_uid_view_hashed &dofs_col_view,
  int verb
) {
  PetscFunctionBegin;

  // get adjeacent element
  mofem_ent_ptr = const_cast<MoFEMEntity*>(mit_row->get_MoFEMEntity_ptr());

  AdjByEnt::iterator adj_miit,hi_adj_miit;
  adj_miit = entFEAdjacencies.get<Unique_mi_tag>().lower_bound(mofem_ent_ptr->get_global_unique_id());
  hi_adj_miit = entFEAdjacencies.get<Unique_mi_tag>().upper_bound(mofem_ent_ptr->get_global_unique_id());

  dofs_col_view.clear();
  for (; adj_miit != hi_adj_miit; adj_miit++) {
    if (adj_miit->by_other&BYROW) {
      if ((adj_miit->EntMoFEMFiniteElement_ptr->get_id()&p_miit->get_BitFEId()).none()) {
        // if element is not part of problem
        continue;
      }
      if ((adj_miit->EntMoFEMFiniteElement_ptr->get_BitRefLevel()&mit_row->get_BitRefLevel()).none()) {
        // if entity is not problem refinement level
        continue;
      }

      if (verb > 2) {
          stringstream ss;

          ss << "rank " << rAnk << ":  numered_dofs_cols" << endl;
          DofMoFEMEntity_multiIndex_uid_view::iterator dit, hi_dit;
          dit = adj_miit->EntMoFEMFiniteElement_ptr->col_dof_view.begin();
          hi_dit = adj_miit->EntMoFEMFiniteElement_ptr->col_dof_view.end();

          for (; dit != hi_dit; dit++) {
            ss << "\t" << **dit << endl;
          }
          PetscSynchronizedPrintf(comm, "%s", ss.str().c_str());
      }

      ierr = adj_miit->EntMoFEMFiniteElement_ptr->
        get_MoFEMFiniteElement_col_dof_view(
          p_miit->numered_dofs_cols,
          dofs_col_view,
          Interface::UNION); CHKERRQ(ierr);

    }

  }

  PetscFunctionReturn(0);
}

template<typename TAG>
PetscErrorCode CreateRowComressedADJMatrix::createMat(
  const string &name,Mat *M,const MatType type,PetscInt **_i,PetscInt **_j,PetscScalar **_v,
  const bool no_diagonals,int verb
) {
  PetscFunctionBegin;
  PetscLogEventBegin(USER_EVENT_createMat,0,0,0,0);
  if(verb==-1) verb = verbose;

  typedef typename boost::multi_index::index<NumeredDofMoFEMEntity_multiIndex,TAG>::type NumeredDofMoFEMEntitysByIdx;

  // Find problem by name FIXME: this should be outsourced to other function,
  // where to this function only problem pointer is passed.
  ProblemsByName &pRoblems_set = pRoblems.get<Problem_mi_tag>();
  ProblemsByName::iterator p_miit = pRoblems_set.find(name);
  if(p_miit==pRoblems_set.end()) SETERRQ1(PETSC_COMM_SELF,MOFEM_NOT_FOUND,"problem < %s > is not found (top tip: check spelling)",name.c_str());

  // Get multi-indices for rows and columns
  const NumeredDofMoFEMEntitysByIdx &dofs_row_by_idx = p_miit->numered_dofs_rows.get<TAG>();
  const NumeredDofMoFEMEntitysByIdx &dofs_col_by_idx = p_miit->numered_dofs_cols.get<TAG>();
  DofIdx nb_dofs_row = p_miit->get_nb_dofs_row();
  if(nb_dofs_row == 0) {
    SETERRQ1(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"problem <%s> has zero rows",name.c_str());
  }

  // Get adjacencies form other processors
  map<int, vector<int> > adjacent_dofs_on_other_parts;

  // If not partitioned set petsc layout for matrix. If partitioned need to get
  // adjacencies form other parts. Note if algebra is only partitioned no need
  // to collect adjacencies form other entities. Those are already on mesh
  // which is assumed that is on each processor the same.
  typename boost::multi_index::index<NumeredDofMoFEMEntity_multiIndex,TAG>::type::iterator miit_row,hi_miit_row;
  if(TAG::IamNotPartitioned) {

    // Get range of local indices
    PetscLayout layout;
    ierr = PetscLayoutCreate(comm, &layout); CHKERRQ(ierr);
    ierr = PetscLayoutSetBlockSize(layout, 1); CHKERRQ(ierr);
    ierr = PetscLayoutSetSize(layout, nb_dofs_row); CHKERRQ(ierr);
    ierr = PetscLayoutSetUp(layout); CHKERRQ(ierr);
    PetscInt rstart, rend;
    ierr = PetscLayoutGetRange(layout, &rstart, &rend); CHKERRQ(ierr);
    ierr = PetscLayoutDestroy(&layout); CHKERRQ(ierr);
    if (verb > 0) {
      PetscSynchronizedPrintf(comm, "\tcreate_Mat: row lower %d row upper %d\n", rstart, rend);
      PetscSynchronizedFlush(comm, PETSC_STDOUT);
    }
    miit_row = dofs_row_by_idx.lower_bound(rstart);
    hi_miit_row = dofs_row_by_idx.lower_bound(rend);
    if (distance(miit_row, hi_miit_row) != rend-rstart) {
      SETERRQ4(PETSC_COMM_SELF,PETSC_ERR_ARG_SIZ,
        "data inconsistency, distance(miit_row,hi_miit_row) != rend - rstart (%d != %d - %d = %d) ",
        distance(miit_row, hi_miit_row), rend, rstart, rend-rstart
      );
    }
  } else {
    //get adjacent nodes on other partitions
    vector<vector<int> > dofs_vec(sIze);

    MoFEMEntity *mofem_ent_ptr = NULL;
    NumeredDofMoFEMEntity_multiIndex_uid_view_hashed dofs_col_view;

    typename boost::multi_index::index<NumeredDofMoFEMEntity_multiIndex,TAG>::type::iterator mit_row,hi_mit_row;
    mit_row = dofs_row_by_idx.begin();
    hi_mit_row = dofs_row_by_idx.end();
    for(;mit_row!=hi_mit_row;mit_row++) {

      // Shared or multishared row and not owned. Those data should be send to other side.

      // Get entity adjacencies, no need to repeat that operation for dofs when
      // are on the same entity. For simplicity is assumed that those sheered the
      // same adjacencies.
      unsigned char pstatus = mit_row->get_pstatus();
      if((pstatus & PSTATUS_NOT_OWNED) && (pstatus&(PSTATUS_SHARED|PSTATUS_MULTISHARED))) {

        if( (mofem_ent_ptr == NULL) ? 1 : mofem_ent_ptr->get_global_unique_id() != mit_row->get_MoFEMEntity_ptr()->get_global_unique_id() ) {
          // get entity adjacencies
          ierr = getEntityAdjacenies<TAG>(p_miit,mit_row,mofem_ent_ptr,dofs_col_view,verb); CHKERRQ(ierr);
          // Add that row. Patterns is that first index is row index, second is
          // size of adjacencies after that follows column adjacencies.
          int owner = mit_row->get_owner_proc();
          dofs_vec[owner].push_back(TAG::get_index(mit_row)); 	// row index
          dofs_vec[owner].push_back(dofs_col_view.size()); 	// nb. of column adjacencies
          // add adjacent cools
          NumeredDofMoFEMEntity_multiIndex_uid_view_hashed::iterator cvit;
          cvit = dofs_col_view.begin();
          for(; cvit!=dofs_col_view.end(); cvit++) {

            int col_idx = TAG::get_index(*cvit);
            if(col_idx<0) {
              ostringstream zz;
              zz << "rank " << rAnk << " ";
              zz << *(*cvit) << endl;
              SETERRQ(PETSC_COMM_SELF,PETSC_ERR_ARG_SIZ,zz.str().c_str());
            }
            if(col_idx>=p_miit->get_nb_dofs_col()) {
              ostringstream zz;
              zz << "rank " << rAnk << " ";
              zz << *(*cvit) << endl;
              SETERRQ(PETSC_COMM_SELF,PETSC_ERR_ARG_SIZ,zz.str().c_str());
            }
            dofs_vec[owner].push_back(col_idx);
          }
        }
      }
    }

    // Make sure it is a PETSc comm
    ierr = PetscCommDuplicate(comm,&comm,NULL); CHKERRQ(ierr);

    int nsends = 0; 			// number of messages to send
    vector<int> dofs_vec_length(sIze);	// length of the message to proc
    for(int proc = 0;proc<sIze;proc++) {

      if(!dofs_vec[proc].empty()) {

        dofs_vec_length[proc] = dofs_vec[proc].size();
        nsends++;

      } else {

        dofs_vec_length[proc] = 0;

      }

    }

    vector<MPI_Status> status(sIze);

    // Computes the number of messages a node expects to receive
    int nrecvs;	// number of messages received
    ierr = PetscGatherNumberOfMessages(comm,NULL,&dofs_vec_length[0],&nrecvs); CHKERRQ(ierr);

    // Computes info about messages that a MPI-node will receive, including (from-id,length) pairs for each message.
    int *onodes;	// list of node-ids from which messages are expected
    int *olengths;	// corresponding message lengths
    ierr = PetscGatherMessageLengths(comm,nsends,nrecvs,&dofs_vec_length[0],&onodes,&olengths);  CHKERRQ(ierr);

    // Gets a unique new tag from a PETSc communicator.
    int tag;
    ierr = PetscCommGetNewTag(comm,&tag); CHKERRQ(ierr);

    // Allocate a buffer sufficient to hold messages of size specified in
    // olengths. And post Irecvs on these buffers using node info from onodes
    int **rbuf;		// must bee freed by user
    MPI_Request *r_waits; // must bee freed by user

    // rbuf has a pointers to messages. It has size of of nrecvs (number of
    // messages) +1. In the first index a block is allocated,
    // such that rbuf[i] = rbuf[i-1]+olengths[i-1].
    ierr = PetscPostIrecvInt(comm,tag,nrecvs,onodes,olengths,&rbuf,&r_waits); CHKERRQ(ierr);

    MPI_Request *s_waits; // status of sens messages
    ierr = PetscMalloc1(nsends,&s_waits);CHKERRQ(ierr);

    // Send messages
    for(int proc=0,kk=0; proc<sIze; proc++) {
      if(!dofs_vec_length[proc]) continue; // no message to send to this proc
      ierr = MPI_Isend(
        &(dofs_vec[proc])[0], 	// buffer to send
        dofs_vec_length[proc], 	// message length
        MPIU_INT,proc,       	// to proc
        tag,comm,s_waits+kk
      ); CHKERRQ(ierr);
      kk++;
    }

    // Wait for received
    if(nrecvs) {
      ierr = MPI_Waitall(nrecvs,r_waits,&status[0]);CHKERRQ(ierr);
    }
    // Wait for send messages
    if(nsends) {
      ierr = MPI_Waitall(nsends,s_waits,&status[0]);CHKERRQ(ierr);
    }

    for(int kk = 0;kk<nrecvs;kk++) {

      int len = olengths[kk];
      int *data_from_proc = rbuf[kk];

      for(int ii = 0;ii<len;) {

        int row_idx = data_from_proc[ii++];	// get row number
        int nb_adj_dofs = data_from_proc[ii++];	// get nb. of adjacent dofs

        if(debug) {

          DofByGlobalPetscIndex::iterator dit;
          dit = p_miit->numered_dofs_rows.get<PetscGlobalIdx_mi_tag>().find(row_idx);
          if(dit==p_miit->numered_dofs_rows.get<PetscGlobalIdx_mi_tag>().end()) {
            SETERRQ1(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"dof %d can not be found in problem",row_idx);
          }

        }

        for(int jj = 0;jj<nb_adj_dofs;jj++) {
          adjacent_dofs_on_other_parts[row_idx].push_back(data_from_proc[ii++]);
        }

      }

    }

    // Cleaning
    ierr = PetscFree(s_waits); CHKERRQ(ierr);
    ierr = PetscFree(rbuf[0]); CHKERRQ(ierr);
    ierr = PetscFree(rbuf); CHKERRQ(ierr);
    ierr = PetscFree(r_waits); CHKERRQ(ierr);
    ierr = PetscFree(onodes); CHKERRQ(ierr);
    ierr = PetscFree(olengths); CHKERRQ(ierr);

    miit_row = dofs_row_by_idx.lower_bound(rAnk);
    hi_miit_row = dofs_row_by_idx.upper_bound(rAnk);

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

    // add next row to compressed matrix
    i.push_back(j.size());
    if(strcmp(type,MATMPIADJ)==0) {
      DofIdx idx = TAG::get_index(miit_row);
      if(dofs_col_by_idx.find(idx)->get_global_unique_id()!=miit_row->get_global_unique_id()) {
        SETERRQ(PETSC_COMM_SELF,PETSC_ERR_ARG_SIZ,"data insonsistency");
      }
    }

    // Get entity adjacencies, no need to repeat that operation for dofs when
    // are on the same entity. For simplicity is assumed that those share the
    // same adjacencies.
    if( (mofem_ent_ptr == NULL) ? 1 : (mofem_ent_ptr->get_global_unique_id() != miit_row->get_MoFEMEntity_ptr()->get_global_unique_id()) ) {

      // get entity adjacencies
      ierr = getEntityAdjacenies<TAG>(p_miit,miit_row,mofem_ent_ptr,dofs_col_view,verb); CHKERRQ(ierr);
      row_last_evaluated_idx = TAG::get_index(miit_row);

      dofs_vec.resize(0);
      NumeredDofMoFEMEntity_multiIndex_uid_view_hashed::iterator cvit;

      cvit = dofs_col_view.begin();
      for(;cvit!=dofs_col_view.end();cvit++) {

        int idx = TAG::get_index(*cvit);
        dofs_vec.push_back(idx);

        if(idx<0) {
          SETERRQ1(
            PETSC_COMM_SELF,PETSC_ERR_ARG_SIZ,"data inconsistency, dof index is smaller than 0, problem name < %s >",name.c_str()
          );
        }
        if(idx>=p_miit->get_nb_dofs_col()) {

          ostringstream ss;
          ss << "Notes: " << endl;
          ss << *(*cvit) << endl;
          PetscPrintf(comm,"%s\n",ss.str().c_str());
          SETERRQ1(
            PETSC_COMM_SELF,PETSC_ERR_ARG_SIZ,"data inconsistency, dof index is bigger than size of problem, problem name < %s >",name.c_str()
          );

        }

      }

      unsigned char pstatus = miit_row->get_pstatus();
      if( pstatus>0 ) {
        map<int,vector<int> >::iterator mit;
        mit = adjacent_dofs_on_other_parts.find(row_last_evaluated_idx);
        if(mit == adjacent_dofs_on_other_parts.end()) {
          cerr << *miit_row << endl;
          SETERRQ1(
            PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,
            "data inconsistency row_last_evaluated_idx = %d",
            row_last_evaluated_idx
          );
        } else {
          dofs_vec.insert(dofs_vec.end(),mit->second.begin(),mit->second.end());
        }
      }

      sort(dofs_vec.begin(),dofs_vec.end());
      vector<DofIdx>::iterator new_end = unique(dofs_vec.begin(),dofs_vec.end());
      int new_size = distance(dofs_vec.begin(),new_end);
      dofs_vec.resize(new_size);
      if(verb>2) {
        stringstream ss;
        ss << "rank " << rAnk << ": dofs_vec for " << *mofem_ent_ptr << endl;
        PetscSynchronizedPrintf(comm,"%s",ss.str().c_str());
      }

    }

    // Try to be smart reserving memory
    if( j.capacity() < j.size() + dofs_vec.size() ) {

      unsigned int nb_nonzero = j.size() + dofs_vec.size();
      unsigned int average_row_fill = nb_nonzero/i.size() + nb_nonzero % i.size();
      if( j.capacity() < rows_to_fill*average_row_fill ) {
        j.reserve( rows_to_fill*average_row_fill );
      }

    }

    // add indices to compressed matrix
    if(verb>1) {
      PetscSynchronizedPrintf(comm,"rank %d: ",rAnk);
    }
    vector<DofIdx>::iterator diit,hi_diit;
    diit = dofs_vec.begin();
    hi_diit = dofs_vec.end();
    for(;diit!=hi_diit;diit++) {

      if(no_diagonals) {
        if(*diit == TAG::get_index(miit_row)) {
          continue;
        }
      }
      j.push_back(*diit);

      if(verb>1) {
        PetscSynchronizedPrintf(comm,"%d ",*diit);
      }

    }
    if(verb>1) {
      PetscSynchronizedPrintf(comm,"\n",*diit);
    }

  }

  if(verb>1) {
    PetscSynchronizedFlush(comm,PETSC_STDOUT);
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

    // Adjacency matrix used to partition problems, f.e. METIS
    if(i.size()-1 != (unsigned int)nb_loc_row_from_iterators) {
      SETERRQ(PETSC_COMM_SELF,PETSC_ERR_ARG_SIZ,"data inconsistency");
    }
    ierr = MatCreateMPIAdj(comm,i.size()-1,nb_col_dofs,*_i,*_j,PETSC_NULL,M); CHKERRQ(ierr);
    ierr = MatSetOption(*M,MAT_STRUCTURALLY_SYMMETRIC,PETSC_TRUE); CHKERRQ(ierr);

  } else if(strcmp(type,MATMPIAIJ)==0) {

    // Compressed MPIADJ matrix
    if(i.size()-1 != (unsigned int)nb_loc_row_from_iterators) {
      SETERRQ(PETSC_COMM_SELF,PETSC_ERR_ARG_SIZ,"data inconsistency");
    }
    PetscInt nb_local_dofs_row = p_miit->get_nb_local_dofs_row();
    if((unsigned int)nb_local_dofs_row!=i.size()-1) {
      SETERRQ(PETSC_COMM_SELF,PETSC_ERR_ARG_SIZ,"data inconsistency");
    }
    PetscInt nb_local_dofs_col = p_miit->get_nb_local_dofs_col();
    ierr = ::MatCreateMPIAIJWithArrays(
      comm,nb_local_dofs_row,nb_local_dofs_col,nb_row_dofs,nb_col_dofs,*_i,*_j,PETSC_NULL,M
    ); CHKERRQ(ierr);

  } else if(strcmp(type,MATAIJ)==0) {

    // Sequential compressed ADJ matrix
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
  CreateRowComressedADJMatrix *core_ptr = static_cast<CreateRowComressedADJMatrix*>(const_cast<Core*>(this));
  ierr = core_ptr->createMat<Part_mi_tag>(name,Aij,MATMPIAIJ,&_i,&_j,PETSC_NULL,false,verb); CHKERRQ(ierr);
  ierr = PetscFree(_i); CHKERRQ(ierr);
  ierr = PetscFree(_j); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
PetscErrorCode Core::MatCreateSeqAIJWithArrays(const string &name,Mat *Aij,PetscInt **i,PetscInt **j,PetscScalar **v,int verb) {
  PetscFunctionBegin;
  if(verb==-1) verb = verbose;
  CreateRowComressedADJMatrix *core_ptr = static_cast<CreateRowComressedADJMatrix*>(const_cast<Core*>(this));
  ierr = core_ptr->createMat<PetscLocalIdx_mi_tag>(name,Aij,MATAIJ,i,j,v,false,verb); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode Core::partition_problem(const string &name,int verb) {
  PetscFunctionBegin;
  if(verb==-1) verb = verbose;
  if(!(*build_MoFEM&(1<<0))) SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"fields not build");
  if(!(*build_MoFEM&(1<<1))) SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"FEs not build");
  if(!(*build_MoFEM&(1<<2))) SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"entFEAdjacencies not build");
  if(!(*build_MoFEM&(1<<3))) SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"pRoblems not build");
  if(verb>0) {
    PetscPrintf(comm,"Partition problem %s\n",name.c_str());
  }

  typedef NumeredDofMoFEMEntity_multiIndex::index<Idx_mi_tag>::type NumeredDofMoFEMEntitysByIdx;
  typedef MoFEMProblem_multiIndex::index<Problem_mi_tag>::type ProblemsByName;

  // Find problem pointer by name
  ProblemsByName &pRoblems_set = pRoblems.get<Problem_mi_tag>();
  ProblemsByName::iterator p_miit = pRoblems_set.find(name);
  if(p_miit==pRoblems_set.end()) SETERRQ1(PETSC_COMM_SELF,1,"problem with name %s not defined (top tip check spelling)",name.c_str());
  DofIdx nb_dofs_row = p_miit->get_nb_dofs_row();

  int *i,*j;
  Mat Adj;
  if(verb>1) {
    PetscPrintf(comm,"\tcreate Adj matrix\n");
  }

  try {
    CreateRowComressedADJMatrix *core_ptr = static_cast<CreateRowComressedADJMatrix*>(const_cast<Core*>(this));
    ierr = core_ptr->createMat<Idx_mi_tag>(name,&Adj,MATMPIADJ,&i,&j,PETSC_NULL,true,verb); CHKERRQ(ierr);
  } catch (MoFEMException const &e) {
    SETERRQ(PETSC_COMM_SELF,e.errorCode,e.errorMessage);
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
  //#ifdef __APPLE__
  //ierr = PetscBarrier((PetscObject)Adj); CHKERRQ(ierr);
  //#endif
  ierr = MatPartitioningSetAdjacency(part,Adj); CHKERRQ(ierr);
  ierr = MatPartitioningSetFromOptions(part); CHKERRQ(ierr);
  ierr = MatPartitioningSetNParts(part,sIze); CHKERRQ(ierr);
  //#ifdef __APPLE__
  //ierr = PetscBarrier((PetscObject)part); CHKERRQ(ierr);
  //#endif
  ierr = MatPartitioningApply(part,&is); CHKERRQ(ierr);
  if(verb>2) {
    ISView(is,PETSC_VIEWER_STDOUT_WORLD);
  }
  #ifdef __APPLE__
  ierr = PetscBarrier((PetscObject)is); CHKERRQ(ierr);
  #endif

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
  NumeredDofMoFEMEntitysByIdx &dofs_row_by_idx_no_const = const_cast<NumeredDofMoFEMEntitysByIdx&>(p_miit->numered_dofs_rows.get<Idx_mi_tag>());
  NumeredDofMoFEMEntitysByIdx &dofs_col_by_idx_no_const = const_cast<NumeredDofMoFEMEntitysByIdx&>(p_miit->numered_dofs_cols.get<Idx_mi_tag>());
  DofIdx &nb_row_local_dofs = *((DofIdx*)p_miit->tag_local_nbdof_data_row);
  DofIdx &nb_col_local_dofs = *((DofIdx*)p_miit->tag_local_nbdof_data_col);
  nb_row_local_dofs = 0;
  nb_col_local_dofs = 0;
  DofIdx &nb_row_ghost_dofs = *((DofIdx*)p_miit->tag_ghost_nbdof_data_row);
  DofIdx &nb_col_ghost_dofs = *((DofIdx*)p_miit->tag_ghost_nbdof_data_col);
  nb_row_ghost_dofs = 0;
  nb_col_ghost_dofs = 0;
  NumeredDofMoFEMEntitysByIdx::iterator miit_dofs_row = dofs_row_by_idx_no_const.begin();
  NumeredDofMoFEMEntitysByIdx::iterator miit_dofs_col = dofs_col_by_idx_no_const.begin();
  if(verb>1) {
    PetscPrintf(comm,"\tloop problem dofs");
  }

  try {

    for(;miit_dofs_row!=dofs_row_by_idx_no_const.end();miit_dofs_row++,miit_dofs_col++) {
      if(miit_dofs_col==dofs_col_by_idx_no_const.end()) {
        SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"check finite element definition, nb. of rows is not equal to number for columns");
      }
      if(miit_dofs_row->get_global_unique_id()!=miit_dofs_col->get_global_unique_id()) {
        SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"check finite element definition, nb. of rows is not equal to columns");
      }
      if(miit_dofs_row->dof_idx!=miit_dofs_col->dof_idx) {
        SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"check finite element definition, nb. of rows is not equal to columns");
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

  } catch (MoFEMException const &e) {
    SETERRQ(PETSC_COMM_SELF,e.errorCode,e.errorMessage);
  } catch (const std::exception& ex) {
    ostringstream ss;
    ss << "throw in method: " << ex.what() << " at line " << __LINE__ << " in file " << __FILE__ << endl;
    SETERRQ(PETSC_COMM_SELF,MOFEM_STD_EXCEPTION_THROW,ss.str().c_str());
  }

  if(verb>1) {
    PetscPrintf(comm," <- done\n");
  }

  ierr = ISRestoreIndices(is_gather,&part_number);  CHKERRQ(ierr);
  ierr = ISRestoreIndices(is_gather_num,&petsc_idx);  CHKERRQ(ierr);
  ierr = ISDestroy(&is_num); CHKERRQ(ierr);
  ierr = ISDestroy(&is_gather_num); CHKERRQ(ierr);
  ierr = ISDestroy(&is_gather); CHKERRQ(ierr);
  ierr = ISDestroy(&is); CHKERRQ(ierr);
  ierr = MatPartitioningDestroy(&part); CHKERRQ(ierr);
  ierr = MatDestroy(&Adj); CHKERRQ(ierr);
  ierr = printPartitionedProblem(&*p_miit,verb); CHKERRQ(ierr);
  *build_MoFEM |= 1<<4;
  PetscFunctionReturn(0);
}
PetscErrorCode Core::partition_check_matrix_fill_in(const string &problem_name,int row_print,int col_print,int verb) {
  PetscFunctionBegin;
  if(verb==-1) verb = verbose;

  struct TestMatrixFillIn: public FEMethod {
    FieldInterface *mFieldPtr;

    Mat A;
    PetscErrorCode ierr;
    ErrorCode rval;

    int rowPrint,colPrint;

    TestMatrixFillIn(FieldInterface *m_field_ptr,Mat a,int row_print,int col_print):
      mFieldPtr(m_field_ptr),A(a),
      rowPrint(row_print),colPrint(col_print) {};

    PetscErrorCode preProcess() {
      PetscFunctionBegin;
      PetscFunctionReturn(0);
    }

    PetscErrorCode operator()() {
      PetscFunctionBegin;

      if(refinedFiniteElementsPtr->find(fePtr->get_ent())==refinedFiniteElementsPtr->end()) {
        SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"data inconsistency");
      }

      FENumeredDofMoFEMEntity_multiIndex::iterator rit = rowPtr->begin();
      for(;rit!=rowPtr->end();rit++) {

        if(refinedEntitiesPtr->find(rit->get_ent())==refinedEntitiesPtr->end()) {
          SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"data inconsistency");
        }
        if(!rit->get_active()) {
          SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"data inconsistency");
        }

        MoFEMEntityEntMoFEMFiniteElementAdjacencyMap_multiIndex::index<Composite_Unique_mi_tag>::type::iterator ait;
        ait = adjacenciesPtr->get<Composite_Unique_mi_tag>().find(boost::make_tuple(rit->get_MoFEMEntity_ptr()->get_global_unique_id(),fePtr->get_global_unique_id()));
        if(ait==adjacenciesPtr->end()) {
          ostringstream ss;
          ss << *rit << endl;
          ss << *fePtr << endl;
          ss << "dof: " << rit->get_BitRefLevel() << endl;
          ss << "fe: " << fePtr->get_BitRefLevel() << endl;
          ss << "problem: " << problemPtr->get_BitRefLevel() << endl;
          PetscPrintf(mFieldPtr->get_comm(),"%s",ss.str().c_str());
          SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"adjacencies data inconsistency");
        } else {
          LocalUId uid = ait->get_ent_unique_id();
          if(entitiesPtr->find(uid) == entitiesPtr->end()) {
            SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"data inconsistency");
          }
          if(dofsPtr->find(rit->get_global_unique_id())==dofsPtr->end()) {
            SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"data inconsistency");
          }
        }
        int row = rit->get_petsc_gloabl_dof_idx();

        FENumeredDofMoFEMEntity_multiIndex::iterator cit = colPtr->begin();
        for(;cit!=colPtr->end();cit++) {

          if(refinedEntitiesPtr->find(cit->get_ent())==refinedEntitiesPtr->end()) {
            SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"data inconsistency");
          }
          if(!cit->get_active()) {
            SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"data inconsistency");
          }
          int col = cit->get_petsc_gloabl_dof_idx();
          ait = adjacenciesPtr->get<Composite_Unique_mi_tag>().find(boost::make_tuple(cit->get_MoFEMEntity_ptr()->get_global_unique_id(),fePtr->get_global_unique_id()));
          if(ait==adjacenciesPtr->end()) {
            SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"adjacencies data inconsistency");
          } else {
            LocalUId uid = ait->get_ent_unique_id();
            if(entitiesPtr->find(uid) == entitiesPtr->end()) {
              SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"data inconsistency");
            }
            if(dofsPtr->find(cit->get_global_unique_id())==dofsPtr->end()) {
              SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"data inconsistency");
            }
          }

          if(row == rowPrint && col == colPrint) {

            ostringstream ss;
            ss << "fe:\n" << *fePtr << endl;
            ss << "row:\n" << *rit << endl;
            ss << "col:\n" << *cit << endl;

            ss << "fe:\n" << fePtr->get_BitRefLevel() << endl;
            ss << "row:\n" << rit->get_BitRefLevel() << endl;
            ss << "col:\n" << cit->get_BitRefLevel() << endl;

            cerr << ss.str() << endl;

            //PetscPrintf(mFieldPtr->get_comm(),"%s\n",ss.str().c_str());

          }

          ierr = MatSetValue(A,row,col,1,INSERT_VALUES);

          if(cit->get_ent_type()!=MBVERTEX) {


            FENumeredDofMoFEMEntity_multiIndex::index<Composite_Name_Type_And_Side_Number_mi_tag>::type::iterator dit,hi_dit;
            dit = colPtr->get<Composite_Name_Type_And_Side_Number_mi_tag>().lower_bound(boost::make_tuple(cit->get_name(),cit->get_ent_type(),cit->side_number_ptr->side_number));
            hi_dit = colPtr->get<Composite_Name_Type_And_Side_Number_mi_tag>().upper_bound(boost::make_tuple(cit->get_name(),cit->get_ent_type(),cit->side_number_ptr->side_number));
            int nb_dofs_on_ent = distance(dit,hi_dit);

            int max_order = cit->get_max_order();
            if(cit->get_nb_of_coeffs()*cit->get_order_nb_dofs(max_order)!=nb_dofs_on_ent) {
              cerr << "Warning: Number of Dofs in Col diffrent than number of dofs for given entity order "
              << cit->get_nb_of_coeffs()*cit->get_order_nb_dofs(max_order) << " " << nb_dofs_on_ent  << endl;
            }

          }

        }

        if(rit->get_ent_type()!=MBVERTEX) {

          FENumeredDofMoFEMEntity_multiIndex::index<Composite_Name_Type_And_Side_Number_mi_tag>::type::iterator dit,hi_dit;
          dit = rowPtr->get<Composite_Name_Type_And_Side_Number_mi_tag>().lower_bound(boost::make_tuple(rit->get_name(),rit->get_ent_type(),rit->side_number_ptr->side_number));
          hi_dit = rowPtr->get<Composite_Name_Type_And_Side_Number_mi_tag>().upper_bound(boost::make_tuple(rit->get_name(),rit->get_ent_type(),rit->side_number_ptr->side_number));
          int nb_dofs_on_ent = distance(dit,hi_dit);

          int max_order = rit->get_max_order();
          if(rit->get_nb_of_coeffs()*rit->get_order_nb_dofs(max_order) != nb_dofs_on_ent) {
            cerr << "Warning: Number of Dofs in Row diffrent than number of dofs for given entity order "
            << rit->get_nb_of_coeffs()*rit->get_order_nb_dofs(max_order) << " " << nb_dofs_on_ent << endl;
          }

        }

      }

      if(fePtr->fe_ptr->row_dof_view.size()!=fePtr->rows_dofs.size()) {
        cerr << "Warning: FEDof Row size != NumeredFEDof RowSize" << endl;
      }

      if(fePtr->fe_ptr->col_dof_view.size()!=fePtr->cols_dofs.size()) {
        cerr << "Warning: FEDof Row size != NumeredFEDof RowSize" << endl;
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

  if(verb>1) {
    MatView(A,PETSC_VIEWER_STDOUT_WORLD);
  }

  if(verb>2) {
    MatView(A,PETSC_VIEWER_DRAW_WORLD);
    std::string wait;
    std::cin >> wait;
  }

  TestMatrixFillIn method(this,A,row_print,col_print);

  typedef MoFEMProblem_multiIndex::index<Problem_mi_tag>::type ProblemsByName;
  //find p_miit
  ProblemsByName &pRoblems_set = pRoblems.get<Problem_mi_tag>();
  ProblemsByName::iterator p_miit = pRoblems_set.find(problem_name);
  if(p_miit == pRoblems_set.end()) {
    SETERRQ1(PETSC_COMM_SELF,1,"problem < %s > not found (top tip: check spelling)",problem_name.c_str());
  }
  if(verb>0) {
    PetscPrintf(comm,"check problem < %s >\n",problem_name.c_str());
  }

  //Loop all elements in problem and check if assemble is without error
  NumeredMoFEMFiniteElement_multiIndex::iterator fe = p_miit->numeredFiniteElements.begin();
  NumeredMoFEMFiniteElement_multiIndex::iterator hi_fe = p_miit->numeredFiniteElements.end();
  for(;fe!=hi_fe;fe++) {

    if(verb>0) {
      PetscPrintf(comm,"\tcheck element %s\n",fe->get_name().c_str());
    }

    ierr = loop_finite_elements(problem_name,fe->get_name(),method,verb);  CHKERRQ(ierr);

  }

  ierr = MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
  ierr = MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);

  ierr = MatDestroy(&A); CHKERRQ(ierr);

  PetscFunctionReturn(0);
}


}
