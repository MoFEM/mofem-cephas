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
#include <UnknownInterface.hpp>
using namespace MoFEM;
#include <PCMGSetUpViaApproxOrders.hpp>

#undef PETSC_VERSION_RELEASE
#define PETSC_VERSION_RELEASE 1

#if PETSC_VERSION_GE(3,6,0)
  #include <petsc/private/petscimpl.h>
#else
  #include <petsc-private/petscimpl.h>
#endif

#if PETSC_VERSION_GE(3,6,0)
  #include <petsc/private/dmimpl.h> /*I  "petscdm.h"   I*/
  // #include <petsc/private/vecimpl.h> /*I  "petscdm.h"   I*/
#else
  #include <petsc-private/dmimpl.h> /*I  "petscdm.h"   I*/
  #include <petsc-private/vecimpl.h> /*I  "petscdm.h"   I*/
#endif

DMMGViaApproxOrdersCtx::DMMGViaApproxOrdersCtx(): MoFEM::DMCtx() {
}
DMMGViaApproxOrdersCtx::~DMMGViaApproxOrdersCtx() {
  PetscErrorCode ierr;
  for(unsigned int ii = 0;ii<kspOperators.size();ii++) {
    ierr = MatDestroy(&kspOperators[ii]); CHKERRABORT(PETSC_COMM_WORLD,ierr);
  }
}

PetscErrorCode DMMGViaApproxOrdersCtx::queryInterface(const MOFEMuuid& uuid,MoFEM::UnknownInterface** iface) {
  PetscFunctionBegin;
  *iface = NULL;
  if(uuid == IDD_DMMGVIAAPPROXORDERSCTX) {
    *iface = dynamic_cast<DMMGViaApproxOrdersCtx*>(this);
    PetscFunctionReturn(0);
  } else {
    SETERRQ(PETSC_COMM_WORLD,MOFEM_DATA_INCONSISTENCY,"wrong interference");
  }
  PetscErrorCode ierr;
  ierr = DMCtx::queryInterface(uuid,iface); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#define GET_DM_FIELD(DM) \
  MoFEM::UnknownInterface *iface; \
  ierr = ((DMCtx*)DM->data)->queryInterface(IDD_DMMGVIAAPPROXORDERSCTX,&iface); CHKERRQ(ierr); \
  DMMGViaApproxOrdersCtx *dm_field = reinterpret_cast<DMMGViaApproxOrdersCtx*>(iface)

PetscErrorCode DMMGViaApproxOrdersPushBackCoarseningIS(DM dm,IS is,Mat A,Mat *subA) {
  PetscErrorCode ierr;
  PetscValidHeaderSpecific(dm,DM_CLASSID,1);
  PetscFunctionBegin;
  GET_DM_FIELD(dm);
  dm_field->coarseningIS.push_back(is);
  if(is) {
    ierr = PetscObjectReference((PetscObject)is); CHKERRQ(ierr);
  }
  // FIXME: If is not the coarse level it would be better to have shell matrix.
  // It would save memory.
  if(is) {
    ierr = MatGetSubMatrix(A,is,is,MAT_INITIAL_MATRIX,subA); CHKERRQ(ierr);
    dm_field->kspOperators.push_back(*subA);
  } else {
    SETERRQ(PETSC_COMM_WORLD,MOFEM_DATA_INCONSISTENCY,"data inconsistency");
  }
  ierr = PetscObjectReference((PetscObject)dm_field->kspOperators.back()); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode DMMGViaApproxOrdersPopBackCoarseningIS(DM dm) {
  PetscErrorCode ierr;
  PetscValidHeaderSpecific(dm,DM_CLASSID,1);
  PetscFunctionBegin;
  GET_DM_FIELD(dm);
  if(dm_field->coarseningIS.back()) {
    ierr = ISDestroy(&dm_field->coarseningIS.back()); CHKERRQ(ierr);
    dm_field->coarseningIS.pop_back();
  }
  if(dm_field->kspOperators.back()) {
    ierr = MatDestroy(&dm_field->kspOperators.back()); CHKERRQ(ierr);
  }
  dm_field->kspOperators.pop_back();
  PetscFunctionReturn(0);
}

PetscErrorCode DMRegister_MGViaApproxOrders(const char sname[]) {
  PetscErrorCode ierr;
  PetscFunctionBegin;
  ierr = DMRegister(sname,DMCreate_MGViaApproxOrders); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

static PetscErrorCode ksp_set_operators(KSP ksp,Mat A,Mat B,void *ctx) {
  PetscErrorCode ierr;
  PetscFunctionBegin;
  // do nothing
  int M,N;
  ierr = MatGetSize(B,&N,&M); CHKERRQ(ierr);
  cerr << "Mat Op " << M << " " << " " << N << endl;
  PetscFunctionReturn(0);
}

PetscErrorCode DMCreate_MGViaApproxOrders(DM dm) {
  PetscErrorCode ierr;
  PetscValidHeaderSpecific(dm,DM_CLASSID,1);
  PetscFunctionBegin;
  if(!dm->data) {
    dm->data = new DMMGViaApproxOrdersCtx();
    cerr << "create dm" << endl;
  } else {
    ((DMCtx*)(dm->data))->referenceNumber++;
    cerr << "ref number " << ((DMCtx*)(dm->data))->referenceNumber << endl;
  }
  ierr = DMSetOperators_MoFEM(dm); CHKERRQ(ierr);
  dm->ops->creatematrix = DMCreateMatrix_MGViaApproxOrders;
  dm->ops->coarsen = DMCoarsen_MGViaApproxOrders;
  dm->ops->createinterpolation = DMCreateInterpolation_MGViaApproxOrders;
  ierr = DMKSPSetComputeOperators(dm,ksp_set_operators,NULL); CHKERRQ(ierr);
  ierr = DMCoarsenHookAdd(dm,NULL,DMRestrict_MGViaApproxOrders,NULL); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode DMCreateMatrix_MGViaApproxOrders(DM dm,Mat *M) {
  PetscErrorCode ierr;
  PetscValidHeaderSpecific(dm,DM_CLASSID,1);
  PetscFunctionBegin;
  GET_DM_FIELD(dm);

  MPI_Comm comm;
  ierr = PetscObjectGetComm((PetscObject)dm,&comm); CHKERRQ(ierr);
  if(dm_field->kspOperators.empty()) {
    SETERRQ(comm,MOFEM_DATA_INCONSISTENCY,"data inconsistency, operator can not be set");
  }
  int leveldown = dm->leveldown;
  if(dm_field->kspOperators.size()<leveldown) {
    SETERRQ(comm,MOFEM_DATA_INCONSISTENCY,"data inconsistency, no IS for that level");
  }
  *M = dm_field->kspOperators[leveldown];
  ierr = PetscObjectReference((PetscObject)*M); CHKERRQ(ierr);

  int m,n;
  ierr = MatGetSize(*M,&m,&n); CHKERRQ(ierr);

  cerr << "Create matrix " << dm->leveldown << " " << m << " " << n << endl;

  PetscFunctionReturn(0);
}

PetscErrorCode DMCoarsen_MGViaApproxOrders(DM dm, MPI_Comm comm, DM *dmc) {
  PetscErrorCode ierr;
  PetscValidHeaderSpecific(dm,DM_CLASSID,1);
  PetscFunctionBegin;
  GET_DM_FIELD(dm);

  ierr = PetscObjectGetComm((PetscObject)dm,&comm); CHKERRQ(ierr);
  ierr = DMCreate(comm,dmc);CHKERRQ(ierr);
  (*dmc)->data = dm->data;
  ierr = DMSetType(*dmc,(dm_field->problemName).c_str()); CHKERRQ(ierr);
  cerr << "Coarsen " << dm->leveldown << " " << dm->levelup << " " << endl;

  PetscFunctionReturn(0);
}

struct MGShellProjectionMatrix {

  int levelUp,levelDown;
  IS isUp,isDown;
  VecScatter sCatter;
  MGShellProjectionMatrix():
  sCatter(PETSC_NULL) {
  }
  virtual ~MGShellProjectionMatrix() {
  }

};

static PetscErrorCode inerpolation_matrix_destroy(Mat mat) {
  PetscErrorCode ierr;
  PetscFunctionBegin;
  void *void_ctx;
  ierr = MatShellGetContext(mat,&void_ctx); CHKERRQ(ierr);
  MGShellProjectionMatrix *ctx = (MGShellProjectionMatrix*)void_ctx;
  if(ctx->sCatter) {
    ierr = VecScatterDestroy(&ctx->sCatter); CHKERRQ(ierr);
  }
  delete ctx;
  PetscFunctionReturn(0);
}

static PetscErrorCode inerpolation_matrix_mult_generic(Mat mat,Vec x,Vec f,InsertMode addv,ScatterMode mode) {
  PetscErrorCode ierr;
  PetscFunctionBegin;

  void *void_ctx;
  ierr = MatShellGetContext(mat,&void_ctx); CHKERRQ(ierr);
  MGShellProjectionMatrix *ctx = (MGShellProjectionMatrix*)void_ctx;

  int M,N;
  ierr = ISGetSize(ctx->isUp,&M); CHKERRQ(ierr);
  ierr = ISGetSize(ctx->isDown,&N); CHKERRQ(ierr);
  int K,L;
  ierr = VecGetSize(x,&K); CHKERRQ(ierr);
  ierr = VecGetSize(f,&L); CHKERRQ(ierr);

  if(!ctx->sCatter) {

    cerr << "inerpolation_matrix_mult_generic ";
    cerr << ctx->levelDown << " " << ctx->levelUp << " : ";
    cerr << M << " " << N << " : ";
    cerr << K << " " << L << " : ";
    cerr << mode << endl;

    IS is_to_map,is_to_scatter;
    int size;
    if(N>M) {
      size = M;
      is_to_map = ctx->isDown;
      is_to_scatter = ctx->isUp;
    } else {
      size = N;
      is_to_map = ctx->isUp;
      is_to_scatter = ctx->isDown;
    }

    ISLocalToGlobalMapping map;
    ierr = ISLocalToGlobalMappingCreateIS(is_to_map,&map); CHKERRQ(ierr);

    vector<int> coarse_idx(size);
    vector<int> coarse_loc_idx(size);
    const int *idx;
    ierr = ISGetIndices(is_to_scatter,&idx); CHKERRQ(ierr);
    copy(&idx[0],&idx[size],coarse_idx.begin());
    ierr = ISRestoreIndices(is_to_scatter,&idx); CHKERRQ(ierr);
    int nout;
    ierr = ISGlobalToLocalMappingApply(
      map,IS_GTOLM_MASK,size,&*coarse_idx.begin(),&nout,&*coarse_loc_idx.begin()
    ); CHKERRQ(ierr);

    ierr = ISLocalToGlobalMappingDestroy(&map); CHKERRQ(ierr);

    MPI_Comm comm;
    ierr = PetscObjectGetComm((PetscObject)mat,&comm); CHKERRQ(ierr);
    IS is;
    ierr = ISCreateGeneral(comm,size,&*coarse_loc_idx.begin(),PETSC_USE_POINTER,&is); CHKERRQ(ierr);

    if(mode == SCATTER_FORWARD) {
      if(K>L) {
        ierr = VecScatterCreate(x,is,f,PETSC_NULL,&ctx->sCatter); CHKERRQ(ierr);
      } else {
        ierr = VecScatterCreate(x,PETSC_NULL,f,is,&ctx->sCatter); CHKERRQ(ierr);
      }
    } else {
      if(K>L) {
        ierr = VecScatterCreate(f,PETSC_NULL,x,is,&ctx->sCatter); CHKERRQ(ierr);
      } else {
        ierr = VecScatterCreate(f,is,x,PETSC_NULL,&ctx->sCatter); CHKERRQ(ierr);
      }
    }

    ierr = ISDestroy(&is); CHKERRQ(ierr);

  }

  ierr = VecScatterBegin(ctx->sCatter,x,f,addv,mode); CHKERRQ(ierr);
  ierr = VecScatterEnd(ctx->sCatter,x,f,addv,mode); CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

static PetscErrorCode inerpolation_matrix_mult(Mat mat,Vec x,Vec f) {
  PetscErrorCode ierr;
  PetscFunctionBegin;
  cerr << "mult ";
  int M,N;
  ierr = MatGetSize(mat,&N,&M); CHKERRQ(ierr);
  cerr << M << " " << N << " : ";
  ierr = VecGetSize(x,&N); CHKERRQ(ierr);
  ierr = VecGetSize(f,&M); CHKERRQ(ierr);
  cerr << M << " " << N << endl;
  ierr = inerpolation_matrix_mult_generic(mat,x,f,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

static PetscErrorCode inerpolation_matrix_mult_transpose(Mat mat,Vec x,Vec f) {
  PetscErrorCode ierr;
  PetscFunctionBegin;
  cerr << "trans mult" << endl;
  ierr = inerpolation_matrix_mult_generic(mat,x,f,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

static PetscErrorCode inerpolation_matrix_mult_add(Mat mat,Vec x,Vec f) {
  PetscErrorCode ierr;
  PetscFunctionBegin;
  ierr = inerpolation_matrix_mult_generic(mat,x,f,ADD_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

static PetscErrorCode inerpolation_matrix_mult_transpose_add(Mat mat,Vec x,Vec f) {
  PetscErrorCode ierr;
  PetscFunctionBegin;
  ierr = inerpolation_matrix_mult_generic(mat,x,f,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode DMCreateInterpolation_MGViaApproxOrders(DM dm1,DM dm2,Mat *mat,Vec *vec) {
  PetscErrorCode ierr;
  PetscValidHeaderSpecific(dm1,DM_CLASSID,1);
  PetscValidHeaderSpecific(dm2,DM_CLASSID,1);
  PetscFunctionBegin;

  int dm1_leveldown = dm1->leveldown;
  int dm2_leveldown = dm2->leveldown;

  cerr << "Create interpolation matrix ";
  cerr << dm1_leveldown << " " << dm2_leveldown << " ";

  MPI_Comm comm;
  ierr = PetscObjectGetComm((PetscObject)dm1,&comm); CHKERRQ(ierr);

  int m,n,M,N;

  MGShellProjectionMatrix *mat_ctx = new MGShellProjectionMatrix();
  {
    // Coarser mesh
    GET_DM_FIELD(dm1);
    if(dm_field->coarseningIS.size()<dm1_leveldown) {
      SETERRQ(PETSC_COMM_WORLD,MOFEM_DATA_INCONSISTENCY,"data inconsistency");
    }
    mat_ctx->isDown = dm_field->coarseningIS[dm1_leveldown];
    mat_ctx->levelDown = dm1_leveldown;
    ierr = ISGetSize(mat_ctx->isDown,&M); CHKERRQ(ierr);
    ierr = ISGetLocalSize(mat_ctx->isDown,&m); CHKERRQ(ierr);
  }
  {
    // Finer mesh
    GET_DM_FIELD(dm2);
    if(dm_field->coarseningIS.size()<dm2_leveldown) {
      SETERRQ(PETSC_COMM_WORLD,MOFEM_DATA_INCONSISTENCY,"data inconsistency");
    }
    mat_ctx->isUp = dm_field->coarseningIS[dm2_leveldown];
    mat_ctx->levelUp = dm2_leveldown;
    ierr = ISGetSize(mat_ctx->isUp,&N); CHKERRQ(ierr);
    ierr = ISGetLocalSize(mat_ctx->isUp,&n); CHKERRQ(ierr);
  }

  // if(N<M) {
  //   SETERRQ(comm,MOFEM_DATA_INCONSISTENCY,"data inconsistency");
  // }

  cerr << N << " " << M << endl;

  ierr = MatCreateShell(comm,m,n,M,N,(void*)mat_ctx,mat); CHKERRQ(ierr);

  Vec right,left;
  ierr = MatCreateVecs(*mat,&right,&left); CHKERRQ(ierr);
  ierr = inerpolation_matrix_mult_generic(*mat,right,left,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = VecDestroy(&right); CHKERRQ(ierr);
  ierr = VecDestroy(&left); CHKERRQ(ierr);

  ierr = MatShellSetOperation(*mat,MATOP_DESTROY,(void(*)(void))inerpolation_matrix_destroy); CHKERRQ(ierr);
  ierr = MatShellSetOperation(*mat,MATOP_MULT,(void(*)(void))inerpolation_matrix_mult); CHKERRQ(ierr);
  ierr = MatShellSetOperation(*mat,MATOP_MULT_TRANSPOSE,(void(*)(void))inerpolation_matrix_mult_transpose); CHKERRQ(ierr);
  // ierr = MatShellSetOperation(*mat,MATOP_MULT_ADD,(void(*)(void))inerpolation_matrix_mult_add); CHKERRQ(ierr);
  // ierr = MatShellSetOperation(*mat,MATOP_MULT_TRANSPOSE_ADD,(void(*)(void))inerpolation_matrix_mult_transpose_add); CHKERRQ(ierr);

  *vec = PETSC_NULL;

  PetscFunctionReturn(0);
}

PetscErrorCode DMRestrict_MGViaApproxOrders(DM fine,Mat mat,Vec vec,Mat mat2,DM coarse,void *ctx) {
  PetscErrorCode ierr;
  PetscFunctionBegin;
  cerr << "Restrict ";
  int  M,N;
  ierr = MatGetSize(mat,&M,&N); CHKERRQ(ierr);
  cerr << M << " " << " " << N << endl;
  PetscFunctionReturn(0);
}

PetscErrorCode PCMGSetUpViaApproxOrdersCtx::getOptions() {
  PetscFunctionBegin;
  ierr = PetscOptionsBegin(
    PETSC_COMM_WORLD,"",
    "MOFEM Multi-Grid (Orders) pre-conditioner","none"
  ); CHKERRQ(ierr);

  ierr = PetscOptionsInt("-mofem_mg_levels",
    "nb levels of multi-grid solver","",
    2,&nbLevels,PETSC_NULL
  ); CHKERRQ(ierr);

  ierr = PetscOptionsInt(
    "-mofem_mg_coarse_order",
    "approximation order of coarse level","",
    2,&coarseOrder,PETSC_NULL
  ); CHKERRQ(ierr);

  ierr = PetscOptionsInt(
    "-mofem_mg_order_at_last_level",
    "order at last level","",
    100,&orderAtLastLevel,PETSC_NULL
  ); CHKERRQ(ierr);


  ierr = PetscOptionsInt(
    "-mofem_mg_verbose",
    "nb levels of multi-grid solver","",
    0,&verboseLevel,PETSC_NULL
  ); CHKERRQ(ierr);

  ierr = PetscOptionsEnd(); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode PCMGSetUpViaApproxOrdersCtx::createIsAtLevel(int kk,IS *is) {
  PetscFunctionBegin;
  //if is last level, take all remaining orders dofs, if any left
  int order_at_next_level = kk+coarseOrder;
  if(kk == nbLevels-1) {
    order_at_next_level = orderAtLastLevel;
  }
  const MoFEM::FieldInterface *m_field_ptr;
  ierr = DMoFEMGetFieldInterfacePtr(dM,&m_field_ptr); CHKERRQ(ierr);
  const MoFEM::MoFEMProblem *problem_ptr;
  ierr = DMMoFEMGetProblemPtr(dM,&problem_ptr); CHKERRQ(ierr);
  string problem_name = problem_ptr->get_name();
  ierr = const_cast<MoFEM::FieldInterface *>(m_field_ptr)
  ->ISCreateProblemOrder(problem_name,ROW,0,order_at_next_level,is); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode PCMGSetUpViaApproxOrdersCtx::destroyIsAtLevel(int kk,IS *is) {
  PetscFunctionBegin;
  ierr = ISDestroy(is);
  PetscFunctionReturn(0);
}

PetscErrorCode PCMGSetUpViaApproxOrdersCtx::buildProlongationOperator(PC pc,int verb) {
  PetscFunctionBegin;
  verb = verb > verboseLevel ? verb : verboseLevel;

  MPI_Comm comm;
  ierr = PetscObjectGetComm((PetscObject)dM,&comm); CHKERRQ(ierr);

  if(verb>0) {
    PetscPrintf(comm,"set MG levels %u\n",nbLevels);
  }

  int sIze,rAnk;
  MPI_Comm_size(comm,&sIze);
  MPI_Comm_rank(comm,&rAnk);

  vector<IS> is_vec(nbLevels+1);
  vector<int> is_glob_size(nbLevels+1),is_loc_size(nbLevels+1);

  for(int kk = 0;kk<nbLevels;kk++) {

    //get indices up to up to give approximation order
    ierr = createIsAtLevel(kk,&is_vec[kk]); CHKERRQ(ierr);
    ierr = ISGetSize(is_vec[kk],&is_glob_size[kk]); CHKERRQ(ierr);
    ierr = ISGetLocalSize(is_vec[kk],&is_loc_size[kk]); CHKERRQ(ierr);

    if(verb>0) {
      PetscSynchronizedPrintf(comm,
      "Nb. dofs at level [ %d ] global %u local %d\n",
      kk,is_glob_size[kk],is_loc_size[kk]);
    }

    //if no dofs on level kk finish here
    if(is_glob_size[kk]==0) {
      SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"no dofs at level");
    }

  }

  ierr = PCMGSetGalerkin(pc,PETSC_FALSE); CHKERRQ(ierr);
  ierr = PCMGSetLevels(pc,nbLevels,NULL);  CHKERRQ(ierr);

  for(int kk = nbLevels-1;kk!=-1;kk--) {
  // for(int kk = 0;kk!=nbLevels;kk++) {
    Mat subA;
    ierr = DMMGViaApproxOrdersPushBackCoarseningIS(dM,is_vec[kk],A,&subA); CHKERRQ(ierr);
    cerr << "kk " << kk << " " << is_vec[kk] << endl;
    if(subA) {
      ierr = MatDestroy(&subA); CHKERRQ(ierr);
    }
  }
  cerr << "end\n\n";

  /*// prolongation and restriction uses the same matrices
  ierr = PCMGSetGalerkin(pc,PETSC_TRUE); CHKERRQ(ierr);

  map<int,int> idx_map;
  int kk = 1;
  while(kk<nbLevels) {

    if(verb>0) {
      PetscPrintf(comm,"level %d\n",kk);
    }

    int row_loc_size,row_glob_size,col_loc_size,col_glob_size;

    row_loc_size = is_loc_size[kk];		//bigger
    row_glob_size = is_glob_size[kk];
    col_loc_size = is_loc_size[kk-1];		//smaller
    col_glob_size = is_glob_size[kk-1];

    // FIXME: Use MatCreateMPIAIJWithArrays
    Mat R;  //small to big
    ierr = MatCreate(comm,&R); CHKERRQ(ierr);
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
      PetscSynchronizedPrintf(comm,"level %d row start %d row end %d\n",kk,rstart,rend);
      PetscSynchronizedPrintf(comm,"level %d col start %d col end %d\n",kk,cstart,cend);
    }

    const int *row_indices_ptr,*col_indices_ptr;
    ierr = ISGetIndices(is_vec[kk],&row_indices_ptr); CHKERRQ(ierr);
    ierr = ISGetIndices(is_vec[kk-1],&col_indices_ptr); CHKERRQ(ierr);

    idx_map.clear();
    for(int ii = 0;ii<row_loc_size;ii++) {
      idx_map[row_indices_ptr[ii]] = rstart+ii;
    }

    // FIXME: Use MatCreateMPIAIJWithArrays and set array directly
    for(int jj = 0;jj<col_loc_size;jj++) {
      map<int,int>::iterator mit = idx_map.find(col_indices_ptr[jj]);
      if(mit != idx_map.end()) {
        ierr = MatSetValue(R,mit->second,cstart+jj,dIag[kk],INSERT_VALUES); CHKERRQ(ierr);
      }
    }

    ierr = ISRestoreIndices(is_vec[kk],&row_indices_ptr); CHKERRQ(ierr);
    ierr = ISRestoreIndices(is_vec[kk-1],&col_indices_ptr); CHKERRQ(ierr);

    ierr = MatAssemblyBegin(R,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
    ierr = MatAssemblyEnd(R,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);

    if(verb>1) {
      //MatView(R,PETSC_VIEWER_STDOUT_WORLD);
      MatView(R,PETSC_VIEWER_DRAW_WORLD);
      std::string wait;
      std::cin >> wait;
    }

    ierr = PCMGSetInterpolation(pc,kk,R); CHKERRQ(ierr);


    #if PETSC_VERSION_LE(3,5,3)
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
    #endif

    ierr = MatDestroy(&R); CHKERRQ(ierr);

    kk++;

  }*/

  for(unsigned int kk = 0;kk<is_vec.size();kk++) {
    ierr = destroyIsAtLevel(kk,&is_vec[kk]); CHKERRQ(ierr);
  }

  if(verb>0) {
    PetscSynchronizedFlush(comm,PETSC_STDOUT);
  }

  PetscFunctionReturn(0);
}

PetscErrorCode PCMGSetUpViaApproxOrders(PC pc,PCMGSetUpViaApproxOrdersCtx *ctx,int verb) {
  PetscErrorCode ierr;
  PetscValidHeaderSpecific(pc,PC_CLASSID,1);
  PetscFunctionBegin;

  MPI_Comm comm;
  ierr = PetscObjectGetComm((PetscObject)pc,&comm); CHKERRQ(ierr);
  if(verb>0) {
    PetscPrintf(comm,"Start PCMGSetUpViaApproxOrders\n");
  }
  ierr = ctx->getOptions(); CHKERRQ(ierr);
  ierr = ctx->buildProlongationOperator(pc,verb); CHKERRQ(ierr);
  if(verb>0) {
    PetscPrintf(comm,"End PCMGSetUpViaApproxOrders\n");
  }

  PetscFunctionReturn(0);
}
