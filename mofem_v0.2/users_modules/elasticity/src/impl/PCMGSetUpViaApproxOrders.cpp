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
  int nbLevels;				///< number of multi-grid levels
  int verboseLevel;			

  PCMGSetUpViaApproxOrdersCtx(FieldInterface *mfield_ptr,string problem_name): 
    mFieldPtr(mfield_ptr),problemName(problem_name),nbLevels(2),verboseLevel(0) {
  }


  PetscErrorCode getOptions() {
    PetscFunctionBegin;
    ierr = PetscOptionsBegin(PETSC_COMM_WORLD,"",
      "MOFEM Multi-Grid (Orders) pre-conditioner","none"); CHKERRQ(ierr);

    ierr = PetscOptionsInt("-mofem_mg_levels",
      "nb levels of multi-grid solver","",
      2,&nbLevels,PETSC_NULL); CHKERRQ(ierr);
    ierr = PetscOptionsInt("-mofem_mg_verbose",
      "nb levels of multi-grid solver","",
      0,&verboseLevel,PETSC_NULL); CHKERRQ(ierr);

    ierr = PetscOptionsEnd(); CHKERRQ(ierr);
    PetscFunctionReturn(0);
  }

  PetscErrorCode buildProlongationOperator(PC pc,int verb = 0) {
    PetscFunctionBegin;
    verb = verb > verboseLevel ? verb : verboseLevel;

    int sIze,rAnk;
    MPI_Comm_size(mFieldPtr->get_comm(),&sIze);
    MPI_Comm_rank(mFieldPtr->get_comm(),&rAnk);

    vector<IS> is_vec(nbLevels);
    vector<int> is_glob_size(nbLevels),is_loc_size(nbLevels),is_loc_shift(nbLevels);

    for(int kk = 0;kk<nbLevels; kk++) {

      //if is last level, take all remaining orders dofs, if any left
      int order_at_next_level = kk+1;
      if(kk == nbLevels-1) {
	order_at_next_level = 100;
      }

      //get indices up to up to give approximation order
      ierr = mFieldPtr->ISCreateProblemOrder(problemName,ROW,0,order_at_next_level,&is_vec[kk]); CHKERRQ(ierr);
      ierr = ISGetSize(is_vec[kk],&is_glob_size[kk]); CHKERRQ(ierr);
      ierr = ISGetLocalSize(is_vec[kk],&is_loc_size[kk]); CHKERRQ(ierr);

      if(verb>0) {
	PetscSynchronizedPrintf(mFieldPtr->get_comm(),"Nb. dofs at level [ %d ] global %u local %d\n",kk,is_glob_size[kk],is_loc_size[kk]);
      }

      //if no dofs on level kk finish here
      if(is_glob_size[kk]==0) {
	ierr = ISDestroy(&is_vec[kk]); CHKERRQ(ierr);
	is_vec.resize(kk);
	break;
      }

    }

    if(verb>0) {
      PetscPrintf(mFieldPtr->get_comm(),"set MG levels %u\n",is_vec.size());
    }

    ierr = PCMGSetLevels(pc,is_vec.size(),NULL);  CHKERRQ(ierr);
    // prolongation and restriction uses the same matrices
    ierr = PCMGSetGalerkin(pc,PETSC_TRUE); CHKERRQ(ierr); 

    for(unsigned int kk = 1;kk<is_vec.size();kk++) {

      int row_loc_size = is_loc_size[kk];
      int col_loc_size = is_loc_size[kk-1];
      int row_glob_size = is_glob_size[kk];
      int col_glob_size = is_glob_size[kk-1];

      Mat R;
      ierr = MatCreate(mFieldPtr->get_comm(),&R); CHKERRQ(ierr);
      ierr = MatSetSizes(R,row_loc_size,col_loc_size,row_glob_size,col_glob_size); CHKERRQ(ierr);
      ierr = MatSetType(R,MATMPIAIJ); CHKERRQ(ierr);
      ierr = MatMPIAIJSetPreallocation(R,1,PETSC_NULL,0,PETSC_NULL); CHKERRQ(ierr);
    
      //get matrix layout
      PetscLayout rmap,cmap;
      ierr = MatGetLayouts(R,&rmap,&cmap); CHKERRQ(ierr);
      int rstart,rend,cstart,cend;
      ierr = PetscLayoutGetRange(rmap,&rstart,&rend); CHKERRQ(ierr);
      ierr = PetscLayoutGetRange(cmap,&cstart,&cend); CHKERRQ(ierr);

      if(verb>0) {
	PetscSynchronizedPrintf(mFieldPtr->get_comm(),"level %d row start %d row end %d\n",kk,rstart,rend);
	PetscSynchronizedPrintf(mFieldPtr->get_comm(),"level %d col start %d col end %d\n",kk,cstart,cend);
      }

      const int *next_indices_ptr,*indices_ptr;
      ierr = ISGetIndices(is_vec[kk-1],&indices_ptr); CHKERRQ(ierr);
      ierr = ISGetIndices(is_vec[kk],&next_indices_ptr); CHKERRQ(ierr);

      int ii = 0;
      for(;ii<col_loc_size;ii++) {
	int row = rstart+ii;
	int col = cstart+ii;
	ierr = MatSetValue(R,row,col,1,INSERT_VALUES); CHKERRQ(ierr);
      }
      ierr = ISRestoreIndices(is_vec[kk-1],&indices_ptr); CHKERRQ(ierr);
      ierr = ISRestoreIndices(is_vec[kk],&next_indices_ptr); CHKERRQ(ierr);

      ierr = MatAssemblyBegin(R,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
      ierr = MatAssemblyEnd(R,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);

      if(verb>1) {
	//MatView(R,PETSC_VIEWER_STDOUT_WORLD);
	MatView(R,PETSC_VIEWER_DRAW_WORLD);
	std::string wait;
	std::cin >> wait;
      }

      ierr = PCMGSetInterpolation(pc,kk,R); CHKERRQ(ierr);

      { 
	// FIXME: If restriction is not set, MatPtAP generate error. This is rather PETSc than MoFEM problem.
	// Petsc Development GIT revision: v3.5.3-1524-gee900cc  GIT Date: 2015-01-31 17:44:15 -0600

	// [0]PETSC ERROR: [0] MatPtAPSymbolic_MPIAIJ_MPIAIJ line 124 /opt/petsc/src/mat/impls/aij/mpi/mpiptap.c
	// [0]PETSC ERROR: [0] MatPtAP_MPIAIJ_MPIAIJ line 80 /opt/petsc/src/mat/impls/aij/mpi/mpiptap.c
	// [0]PETSC ERROR: [0] MatPtAP line 8458 /opt/petsc/src/mat/interface/matrix.c
	// [0]PETSC ERROR: [0] PCSetUp_MG line 552 /opt/petsc/src/ksp/pc/impls/mg/mg.c
	// [0]PETSC ERROR: [0] KSPSetUp line 220 /opt/petsc/src/ksp/ksp/interface/itfunc.c

	// ==4284== Invalid read of size 8
	// ==4284==    at 0x5CAC873: MatPtAPSymbolic_MPIAIJ_MPIAIJ (mpiptap.c:154)
	// ==4284==    by 0x5CABA74: MatPtAP_MPIAIJ_MPIAIJ (mpiptap.c:83)
	// ==4284==    by 0x5D566D1: MatPtAP (matrix.c:8537)
	// ==4284==    by 0x61F27D6: PCSetUp_MG (mg.c:642)
	// ==4284==    by 0x612BD9C: PCSetUp (precon.c:909)
	// ==4284==    by 0x62C1A51: KSPSetUp (itfunc.c:306)
	// ==4284==    by 0xB98326: main (elasticity.cpp:403)
	// ==4284==  Address 0x208 is not stack'd, malloc'd or (recently) free'd

	Mat RT;
	ierr = MatTranspose(R,MAT_INITIAL_MATRIX,&RT); CHKERRQ(ierr);
	ierr = PCMGSetRestriction(pc,kk,RT); CHKERRQ(ierr);

	if(verb>1) {
	  //MatView(R,PETSC_VIEWER_STDOUT_WORLD);
	  MatView(RT,PETSC_VIEWER_DRAW_WORLD);
	  std::string wait;
	  std::cin >> wait;
	}

	ierr = MatDestroy(&RT); CHKERRQ(ierr);
      }

      ierr = MatDestroy(&R); CHKERRQ(ierr);

    }

    for(unsigned int kk = 0;kk<is_vec.size();kk++) {
      ierr = ISDestroy(&is_vec[kk]); CHKERRQ(ierr);
    }

    if(verb>0) {
      PetscSynchronizedFlush(mFieldPtr->get_comm(),PETSC_STDOUT); 
    }

    PetscFunctionReturn(0);
  }

  PetscErrorCode setCoarseKsp(PC pc,int verb = 0) {
    PetscFunctionBegin;
    verb = verb > verboseLevel ? verb : verboseLevel;
  
    KSP coarse_ksp;
    ierr = PCMGGetCoarseSolve(pc,&coarse_ksp); CHKERRQ(ierr);
    PC coarse_pc;
    ierr = KSPGetPC(coarse_ksp,&coarse_pc); CHKERRQ(ierr);
    ierr = PCFactorSetMatSolverPackage(coarse_pc,MATSOLVERMUMPS); CHKERRQ(ierr);

    PetscFunctionReturn(0);
  }

};

PetscErrorCode PCMGSetUpViaApproxOrders(PC pc,FieldInterface *mfield_ptr,const char problem_name[],int verb) {
  PetscErrorCode ierr;
  PetscValidHeaderSpecific(pc,PC_CLASSID,1);
  PetscFunctionBegin;

  if(verb>0) {
    PetscPrintf(mfield_ptr->get_comm(),"Start PCMGSetUpViaApproxOrders\n");
  }

  MPI_Comm comm;
  ierr = PetscObjectGetComm((PetscObject)pc,&comm); CHKERRQ(ierr);
  int result = 0;
  MPI_Comm_compare(comm,mfield_ptr->get_comm(),&result);
  if(result > MPI_CONGRUENT) {
    SETERRQ(PETSC_COMM_SELF,MOFEM_NOT_IMPLEMENTED,"MoFEM and PC have to use the same communicator");
  }
  

  PCMGSetUpViaApproxOrdersCtx ctx(mfield_ptr,problem_name); 
  ierr = ctx.getOptions(); CHKERRQ(ierr);
  ierr = ctx.buildProlongationOperator(pc,verb); CHKERRQ(ierr);
  ierr = ctx.setCoarseKsp(pc,verb); CHKERRQ(ierr);

  if(verb>0) {
    PetscPrintf(mfield_ptr->get_comm(),"End PCMGSetUpViaApproxOrders\n");
  }

  PetscFunctionReturn(0);
}
