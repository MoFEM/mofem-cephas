/** \file PCMGSetUpViaApproxOrders.cpp
 * \brief useful compiler directives and definitions
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

#include <MoFEM.hpp>
using namespace MoFEM;
#include <PCMGSetUpViaApproxOrders.hpp>

static PetscErrorCode ierr;
static ErrorCode rval;

struct PCMGSetUpViaApproxOrdersCtx {

  FieldInterface *mFieldPtr;		///< MoFEM interface
  string problemName;			///< Problem name

  PCMGSetUpViaApproxOrdersCtx(FieldInterface *mfield_ptr,string &problem_name): 
    mFieldPtr(mfield_ptr),problemName(problem_name) {
  }

  int nbLevels;				///< number of multi-grid levels

  PetscErrorCode getOptions() {
    PetscFunctionBegin;
    ierr = PetscOptionsBegin(PETSC_COMM_WORLD,"","MOFEM Multi-Grid (Orders) pre-conditioner","none"); CHKERRQ(ierr);

    nbLevels = 2;
    ierr = PetscOptionsInt("-mg_levels",
      "nb levels of multi-grid solver","",
      2,&nbLevels,PETSC_NULL); CHKERRQ(ierr);

    ierr = PetscOptionsEnd(); CHKERRQ(ierr);
    PetscFunctionReturn(0);
  }

  struct ShellMatCtx {

    VecScatter isLevel;
    VecScatter isLevelP1;

    Vec vLevel;
    Vec vLevelP1;
  
  }

  vector<ShellMatCtx> shellMatCtxVec;

  PetscErroCode buildProlongationOperator(PC pc) {
    PetscFunctionBegin;

    vector<IS> is_vec(nbLevels);
    vectro<int> is_glob_size(nbLevels),is_loc_size(nbLevels);
    for(int kk = 1;kk<=nbLevels; kk++) {

      int next_level = kk;
      if(kk == nbLevels) {
	next_level = 100;
      }

      ierr = ISCreateProblemOrder(problemName,ROW,0,next_level,&is_vec[kk-1]); CHKERRQ(ierr);
      ierr = ISGetSize(is[kk],&is_glob_size[kk-1]); CHKERRQ(ierr);
      ierr = ISGetLocalSize(is[kk],&is_loc_size[kk-1]); CHKERRQ(ierr);

      if(is_glob_size[kk-1]==0) {
	ierr = ISDestroy(&is_vec[kk-1]); CHKERRQ(ierr);
	is_vec.resize(kk);
	break;
      }

    }

    shellMatCtxVec.resize(kk);

    for(int kk = 0;kk<is_vec.size()-1;kk++) {

      Mat R;	///< prolongation matrix, transpose is used as restriction
      ierr = MatCreateShell(mFieldPtr->get_comm(),
	is_loc_size[kk],is_loc_size[kk+1],
	is_glob_size[kk],is_glob_size[kk+1],
	(void*)&shellMatCtxVec[kk],&R); CHKERRQ(ierr);

      #if PETSC_VERSION_GE(3,5,3) 
	ierr = MatCreateVecs(R,&shellMatCtxVec[kk].vecLevel,&shellMatCtxVec[kk].vecLevelP1); CHKERRQ(ierr);
      #else 
	ierr = MatGetVecs(R,&shellMatCtxVec[kk].vecLevel,&shellMatCtxVec[kk].vecLevelP1); CHKERRQ(ierr);
      #endif

      ierr = VecScatterCreate(shellMatCtxVec[kk].vecLevel,PETSC_NULL,is_vec[kk],&shellMatCtxVec[kk].isLevel); CHKERRQ(ierr);
      ierr = VecScatterCreate(shellMatCtxVec[kk].vecLevelP1,PETSC_NULL,is_vec[kk],&shellMatCtxVec[kk].isLevelP1); CHKERRQ(ierr);

      ierr = PCMGSetInterpolation(pc,k+1,R); CHKERRQ(ierr);
      ierr = MatDestroy(&R); CHKERRQ(ierr); // matrix exist, i.e. petsc version of smart pointer

    }

    PetscFunctionReturn(0);
  }

  friend PetscErrorCode MultOpR(Mat A,Vec x,Vec f);

};

static PetscErrorCode MultOpR(Mat A,Vec x,Vec f) {
  PetscFunctionBegin;
  PetscFunctionReturn(0);
}



PetcErrorCode PCMGSetUpViaApproxOrders(PC pc,FieldInterface *mfild_ptr,const char problem_name[]) {
  PetscErrorCode ierr;
  PetscValidHeaderSpecific(pc,PC_CLASSID,1);
  PetscFunctionBegin;

  MPI_Comm comm;
  ierr = PetscObjectGetComm((PetscObject)pc,&comm); CHKERRQ(ierr);
  int result = 0;
  MPI_Comm_compare(comm,mfield_ptr->get_comm(),&result);
  if(result > MPI_CONGRUENT) {
    SETERRQ(PETSC_COMM_SELF,MOFEM_NOT_IMPLEMENTED,"MoFEM and PC have to use the same communicator");
  }
  
  rval = PCMGSetUpViaApproxOrdersCtx ctx(m_field_ptr,problem_name); CHKERR(rval);
  ierr = ctx.getOptions(); CHKERRQ(ierr);





  
  PetscFunctionReturn(0);
}
