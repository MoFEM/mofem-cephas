/** \file PCMGSetUpViaApproxOrders.cpp
 * \brief implementation of multi-grid solver for p- adaptivity 
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

#include <petsc-private/petscimpl.h> 

static PetscErrorCode ierr;
//static ErrorCode rval;

struct PCMGSetUpViaApproxOrdersCtx {

  FieldInterface *mFieldPtr;		///< MoFEM interface
  string problemName;			///< Problem name

  PCMGSetUpViaApproxOrdersCtx(FieldInterface *mfield_ptr,string problem_name): 
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

    bool iNitialized;
    ShellMatCtx(): 
      iNitialized(false)  {};

    IS isLevel;
    VecScatter scatterLevel; 			///< vector scatter at level

    PetscErrorCode setUp(IS is) {
      PetscValidHeaderSpecific(is,VEC_CLASSID,1);
      PetscFunctionBegin;

      ierr = PetscObjectReference((PetscObject)is); CHKERRQ(ierr);		
      isLevel = is;

      PetscFunctionReturn(0);
    }

    static PetscErrorCode opMultR(Mat R,Vec x,Vec f) {
      PetscFunctionBegin;
      void *void_ctx;
      ierr = MatShellGetContext(R,&void_ctx); CHKERRQ(ierr);
      ShellMatCtx *ctx = (ShellMatCtx*)void_ctx;
      if(!ctx->iNitialized) {
	ierr = VecScatterCreate(x,PETSC_NULL,f,ctx->isLevel,&ctx->scatterLevel); CHKERRQ(ierr);
	ierr = ISDestroy(&ctx->isLevel); CHKERRQ(ierr);
	ctx->iNitialized = true;
      }
      ierr = VecScatterBegin(ctx->scatterLevel,x,f,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
      ierr = VecScatterEnd(ctx->scatterLevel,x,f,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
      PetscFunctionReturn(0);
    }

    static PetscErrorCode opMultTransR(Mat R,Vec x,Vec f) {
      PetscFunctionBegin;
      void *void_ctx;
      ierr = MatShellGetContext(R,&void_ctx); CHKERRQ(ierr);
      ShellMatCtx *ctx = (ShellMatCtx*)void_ctx;
      if(!ctx->iNitialized) {
	ierr = VecScatterCreate(f,PETSC_NULL,x,ctx->isLevel,&ctx->scatterLevel); CHKERRQ(ierr);
	ierr = ISDestroy(&ctx->isLevel); CHKERRQ(ierr);
	ctx->iNitialized = true;
      }
      ierr = VecScatterBegin(ctx->scatterLevel,x,f,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
      ierr = VecScatterEnd(ctx->scatterLevel,x,f,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
      PetscFunctionReturn(0);
    }

    ///< destroy shell matrix context data
    static PetscErrorCode opDestroyR(Mat R) {
      PetscFunctionBegin;
      void *void_ctx;
      ierr = MatShellGetContext(R,&void_ctx); CHKERRQ(ierr);
      ShellMatCtx *ctx = (ShellMatCtx*)void_ctx;
      if(ctx->iNitialized) {
	ierr = VecScatterDestroy(&ctx->scatterLevel); CHKERRQ(ierr);
      }
      PetscFunctionReturn(0);
    }

  };

  vector<ShellMatCtx> shellMatCtxVec;

  PetscErrorCode buildProlongationOperator(PC pc) {
    PetscFunctionBegin;

    vector<IS> is_vec(nbLevels);
    vector<int> is_glob_size(nbLevels),is_loc_size(nbLevels);

    for(int kk = 0;kk<nbLevels; kk++) {

      //if is last level, take all remaining orders dofs, if any left
      int next_level = kk;
      if(kk == nbLevels) {
	next_level = 100;
      }

      //get indices up to up to give approximation order
      ierr = mFieldPtr->ISCreateProblemOrder(problemName,ROW,0,next_level,&is_vec[kk]); CHKERRQ(ierr);
      ierr = ISGetSize(is_vec[kk],&is_glob_size[kk]); CHKERRQ(ierr);
      ierr = ISGetLocalSize(is_vec[kk],&is_loc_size[kk]); CHKERRQ(ierr);

      //if no dofs on level kk finish here
      if(is_glob_size[kk]==0) {
	ierr = ISDestroy(&is_vec[kk]); CHKERRQ(ierr);
	is_vec.resize(kk);
	break;
      }

    }

    

    // prolongation matrices 
    shellMatCtxVec.resize(is_vec.size());
    vector<int> indices_from_level_to_next_levevl;

    for(unsigned int kk = 1;kk<is_vec.size();kk++) {

      Mat R; ///< prolongation matrix, transpose is used as restriction
      // I know, this makes first index not used, small loss, code simplicity more important
      ierr = MatCreateShell(mFieldPtr->get_comm(),
	is_loc_size[kk],is_loc_size[kk],
	is_glob_size[kk],is_glob_size[kk],
	(void*)&shellMatCtxVec[kk],&R); CHKERRQ(ierr);

      // set operators
      ierr = MatShellSetOperation(R,MATOP_DESTROY,(void(*)(void))ShellMatCtx::opDestroyR); CHKERRQ(ierr);
      ierr = MatShellSetOperation(R,MATOP_MULT,(void(*)(void))ShellMatCtx::opMultR); CHKERRQ(ierr);
      ierr = MatShellSetOperation(R,MATOP_MULT_TRANSPOSE,(void(*)(void))ShellMatCtx::opMultTransR); CHKERRQ(ierr);

      const int *next_indices_ptr,*indices_ptr;
      ierr = ISGetIndices(is_vec[kk-1],&indices_ptr); CHKERRQ(ierr);
      ierr = ISGetIndices(is_vec[kk],&next_indices_ptr); CHKERRQ(ierr);

      int loc_size = is_loc_size[kk-1];
      indices_from_level_to_next_levevl.resize(loc_size);

      int ii = 0,jj = 0;
      for(;ii<loc_size;ii++) {

	// it can be done like that since indices are sorted
	if(indices_ptr[ii] == next_indices_ptr[jj]) {
	  indices_from_level_to_next_levevl[ii] = jj;
	} else {
	  jj++;
	}

      }
	  
      ierr = ISRestoreIndices(is_vec[kk],&indices_ptr); CHKERRQ(ierr);
      ierr = ISRestoreIndices(is_vec[kk],&indices_ptr); CHKERRQ(ierr);

      IS is;
      ISCreateGeneral(mFieldPtr->get_comm(),loc_size,&indices_from_level_to_next_levevl[0],PETSC_USE_POINTER,&is);
      ierr = shellMatCtxVec[kk].setUp(is); CHKERRQ(ierr);
      ierr = ISDestroy(&is); CHKERRQ(ierr);

      ierr = PCMGSetInterpolation(pc,kk,R); CHKERRQ(ierr);
      ierr = MatDestroy(&R); CHKERRQ(ierr);

    }

    PetscFunctionReturn(0);
  }

};

PetscErrorCode PCMGSetUpViaApproxOrders(PC pc,FieldInterface *mfield_ptr,const char problem_name[]) {
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
  
  PCMGSetUpViaApproxOrdersCtx ctx(mfield_ptr,problem_name); 
  ierr = ctx.getOptions(); CHKERRQ(ierr);
  ierr = ctx.buildProlongationOperator(pc); CHKERRQ(ierr);

  PetscFunctionReturn(0);
}
