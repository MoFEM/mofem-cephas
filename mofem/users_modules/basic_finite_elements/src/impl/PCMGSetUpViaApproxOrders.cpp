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

DMMGViaApproxOrdersCtx::DMMGViaApproxOrdersCtx():
  MoFEM::DMCtx(),
  aO(PETSC_NULL) {
    // std::cerr << "create dm\n";
}
DMMGViaApproxOrdersCtx::~DMMGViaApproxOrdersCtx() {
  PetscErrorCode ierr;
  for(unsigned int ii = 0;ii<coarseningIS.size();ii++) {
    ierr = ISDestroy(&coarseningIS[ii]); CHKERRABORT(PETSC_COMM_WORLD,ierr);
  }
  for(unsigned int ii = 0;ii<kspOperators.size();ii++) {
    ierr = MatDestroy(&kspOperators[ii]); CHKERRABORT(PETSC_COMM_WORLD,ierr);
  }
  if(aO) {
    // std::cerr << "destroy ao\n";
    ierr = AODestroy(&aO); CHKERRABORT(PETSC_COMM_WORLD,ierr);
  }
  // std::cerr << "destroy dm\n";
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

PetscErrorCode DMMGViaApproxOrdersSetAO(DM dm,AO ao) {
  PetscErrorCode ierr;
  PetscValidHeaderSpecific(dm,DM_CLASSID,1);
  PetscFunctionBegin;
  GET_DM_FIELD(dm);
  if(dm_field->aO) {
    //std::cerr << dm_field->aO << std::endl;
    ierr = AODestroy(&dm_field->aO); CHKERRQ(ierr);
    // std::cerr << "destroy ao when adding\n";
  }
  dm_field->aO = ao;
  ierr = PetscObjectReference((PetscObject)ao); CHKERRQ(ierr);
  // std::cerr << "add ao\n";
  PetscFunctionReturn(0);
}

PetscErrorCode DMMGViaApproxOrdersGetCoarseningISSize(DM dm,int *size) {
  PetscErrorCode ierr;
  PetscValidHeaderSpecific(dm,DM_CLASSID,1);
  PetscFunctionBegin;
  GET_DM_FIELD(dm);
  *size = dm_field->coarseningIS.size();
  PetscFunctionReturn(0);
}

PetscErrorCode DMMGViaApproxOrdersPushBackCoarseningIS(DM dm,IS is,Mat A,Mat *subA,bool create_sub_matrix) {
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
    IS is2 = is;
    if(dm_field->aO) {
      ierr = ISDuplicate(is,&is2); CHKERRQ(ierr);
      ierr = ISCopy(is,is2); CHKERRQ(ierr);
      ierr = AOApplicationToPetscIS(dm_field->aO,is2); CHKERRQ(ierr);
    }
    if(create_sub_matrix) {
      ierr = MatGetSubMatrix(A,is2,is2,MAT_INITIAL_MATRIX,subA); CHKERRQ(ierr);
    }
    if(dm_field->aO) {
      ierr = ISDestroy(&is2); CHKERRQ(ierr);
    }
    dm_field->kspOperators.push_back(*subA);
    ierr = PetscObjectReference((PetscObject)(*subA)); CHKERRQ(ierr);
  } else {
    SETERRQ(PETSC_COMM_WORLD,MOFEM_DATA_INCONSISTENCY,"data inconsistency");
  }
  PetscInfo(dm,"Push back IS to DMMGViaApproxOrders\n");
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
  PetscInfo(dm,"Pop back IS to DMMGViaApproxOrders\n");
  PetscFunctionReturn(0);
}

PetscErrorCode DMMGViaApproxOrdersReplaceCoarseningIS(DM dm,IS *is_vec,int nb_elems,Mat A,int verb) {
  PetscErrorCode ierr;
  PetscValidHeaderSpecific(dm,DM_CLASSID,1);
  PetscFunctionBegin;
  GET_DM_FIELD(dm);
  int nb_no_changed = 0;
  int nb_replaced = 0;
  int nb_deleted = 0;
  int nb_added = 0;
  std::vector<IS>::iterator it;
  it = dm_field->coarseningIS.begin();
  int ii = 0;
  for(;it!=dm_field->coarseningIS.end();it++,ii++) {
    if(ii<nb_elems) {
      PetscBool  flg;
      ierr = ISEqual(*it,is_vec[ii],&flg); CHKERRQ(ierr);
      if(!flg) {
        ierr = ISDestroy(&*it); CHKERRQ(ierr);
        ierr = MatDestroy(&dm_field->kspOperators[ii]); CHKERRQ(ierr);
        *it = is_vec[ii];
        ierr = PetscObjectReference((PetscObject)is_vec[ii]); CHKERRQ(ierr);
        if(ii<nb_elems-1) {
          IS is = is_vec[ii];
          if(dm_field->aO) {
            ierr = ISDuplicate(is_vec[ii],&is); CHKERRQ(ierr);
            ierr = ISCopy(is_vec[ii],is); CHKERRQ(ierr);
            ierr = AOApplicationToPetscIS(dm_field->aO,is); CHKERRQ(ierr);
          }
          Mat subA;
          ierr = MatGetSubMatrix(A,is,is,MAT_INITIAL_MATRIX,&subA); CHKERRQ(ierr);
          ierr = PetscObjectReference((PetscObject)subA); CHKERRQ(ierr);
          dm_field->kspOperators[ii] = subA;
          ierr = MatDestroy(&subA); CHKERRQ(ierr);
          if(dm_field->aO) {
            ierr = ISDestroy(&is); CHKERRQ(ierr);
          }
        } else {
          ierr = PetscObjectReference((PetscObject)A); CHKERRQ(ierr);
          dm_field->kspOperators[ii] = A;
        }
        nb_replaced++;
      }
    } else {
      nb_no_changed++;
      continue;
    }
  }
  if(dm_field->coarseningIS.size()<nb_elems) {
    for(;ii<nb_elems-1;ii++) {
      Mat subA;
      ierr = DMMGViaApproxOrdersPushBackCoarseningIS(dm,is_vec[ii],A,&subA,true); CHKERRQ(ierr);
      ierr = MatDestroy(&subA); CHKERRQ(ierr);
      nb_added++;
    }
    ierr = DMMGViaApproxOrdersPushBackCoarseningIS(dm,is_vec[ii],A,&A,false); CHKERRQ(ierr);
    nb_added++;
  } else {
    for(;ii<dm_field->coarseningIS.size();ii++) {
      ierr = DMMGViaApproxOrdersPopBackCoarseningIS(dm); CHKERRQ(ierr);
      nb_deleted++;
    }
  }
  MPI_Comm comm;
  ierr = PetscObjectGetComm((PetscObject)dm,&comm); CHKERRQ(ierr);
  if(verb>0) {
    PetscPrintf(
      comm,"DMMGViaApproxOrders nb_no_changed = %d, nb_replaced = %d, nb_added = %d, nb_deleted = %d, size = %d\n",
      nb_no_changed,nb_replaced,nb_added,nb_deleted,dm_field->coarseningIS.size()
    );
  }
  PetscInfo(dm,"Replace IS to DMMGViaApproxOrders\n");
  PetscFunctionReturn(0);
}

PetscErrorCode DMRegister_MGViaApproxOrders(const char sname[]) {
  PetscErrorCode ierr;
  PetscFunctionBegin;
  ierr = DMRegister(sname,DMCreate_MGViaApproxOrders); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

static PetscErrorCode ksp_set_operators(KSP ksp,Mat A,Mat B,void *ctx) {
  PetscFunctionBegin;
  PetscFunctionReturn(0);
}

PetscErrorCode DMCreate_MGViaApproxOrders(DM dm) {
  PetscErrorCode ierr;
  PetscValidHeaderSpecific(dm,DM_CLASSID,1);
  PetscFunctionBegin;
  if(!dm->data) {
    dm->data = new DMMGViaApproxOrdersCtx();
  } else {
    ((DMCtx*)(dm->data))->referenceNumber++;
  }
  ierr = DMSetOperators_MoFEM(dm); CHKERRQ(ierr);
  dm->ops->creatematrix = DMCreateMatrix_MGViaApproxOrders;
  dm->ops->createglobalvector = DMCreateGlobalVector_MGViaApproxOrders;
  dm->ops->coarsen = DMCoarsen_MGViaApproxOrders;
  dm->ops->createinterpolation = DMCreateInterpolation_MGViaApproxOrders;
  ierr = DMKSPSetComputeOperators(dm,ksp_set_operators,NULL); CHKERRQ(ierr);
  PetscInfo1(dm,"Create DMMGViaApproxOrders reference = %d\n",((DMCtx*)(dm->data))->referenceNumber);
  PetscFunctionReturn(0);
}

PetscErrorCode DMCreateMatrix_MGViaApproxOrders(DM dm,Mat *M) {
  PetscErrorCode ierr;
  PetscValidHeaderSpecific(dm,DM_CLASSID,1);
  PetscFunctionBegin;
  GET_DM_FIELD(dm);

  int leveldown = dm->leveldown;

  if(dm_field->kspOperators.empty()) {
    ierr = DMCreateMatrix_MoFEM(dm,M); CHKERRQ(ierr);
  } else {
    MPI_Comm comm;
    ierr = PetscObjectGetComm((PetscObject)dm,&comm); CHKERRQ(ierr);
    if(dm_field->kspOperators.empty()) {
      SETERRQ(comm,MOFEM_DATA_INCONSISTENCY,"data inconsistency, operator can not be set");
    }
    if(dm_field->kspOperators.size()<leveldown) {
      SETERRQ(comm,MOFEM_DATA_INCONSISTENCY,"data inconsistency, no IS for that level");
    }
    *M = dm_field->kspOperators[dm_field->kspOperators.size()-1-leveldown];
    ierr = PetscObjectReference((PetscObject)*M); CHKERRQ(ierr);
  }

  PetscInfo1(dm,"Create Matrix DMMGViaApproxOrders leveldown = %d\n",leveldown);

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
  ierr = PetscObjectReference((PetscObject)(*dmc)); CHKERRQ(ierr);
  PetscInfo1(dm,"Coarsen DMMGViaApproxOrders leveldown = %d\n",dm->leveldown);

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
    if(sCatter) {
      PetscErrorCode ierr;
      ierr = VecScatterDestroy(&sCatter); CHKERRABORT(PETSC_COMM_WORLD,ierr);
    }
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
    ctx->sCatter = PETSC_NULL;
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

  if(!ctx->sCatter) {

    int M,N;
    ierr = ISGetSize(ctx->isUp,&M); CHKERRQ(ierr);
    ierr = ISGetSize(ctx->isDown,&N); CHKERRQ(ierr);
    int K,L;
    ierr = VecGetSize(x,&K); CHKERRQ(ierr);
    ierr = VecGetSize(f,&L); CHKERRQ(ierr);

    int m,n;
    ierr = ISGetLocalSize(ctx->isUp,&m); CHKERRQ(ierr);
    ierr = ISGetLocalSize(ctx->isDown,&n); CHKERRQ(ierr);
    int k,l;
    ierr = VecGetLocalSize(x,&k); CHKERRQ(ierr);
    ierr = VecGetLocalSize(f,&l); CHKERRQ(ierr);

    IS is_to_map,is_to_scatter;
    int size;
    if(N>M) {
      size = m;
      is_to_map = ctx->isDown;
      is_to_scatter = ctx->isUp;
    } else {
      size = n;
      is_to_map = ctx->isUp;
      is_to_scatter = ctx->isDown;
    }

    IS is;
    ierr = ISDuplicate(is_to_scatter,&is); CHKERRQ(ierr);
    ierr = ISCopy(is_to_scatter,is); CHKERRQ(ierr);

    AO ao;
    ierr = AOCreateMappingIS(is_to_map,PETSC_NULL,&ao); CHKERRQ(ierr);
    ierr = AOApplicationToPetscIS(ao,is); CHKERRQ(ierr);
    ierr = AODestroy(&ao); CHKERRQ(ierr);

    MPI_Comm comm;
    ierr = PetscObjectGetComm((PetscObject)mat,&comm); CHKERRQ(ierr);
    // IS is;
    // ierr = ISCreateGeneral(comm,size,&*coarse_loc_idx.begin(),PETSC_USE_POINTER,&is); CHKERRQ(ierr);

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
  int M,N;
  ierr = MatGetSize(mat,&N,&M); CHKERRQ(ierr);
  ierr = VecGetSize(x,&N); CHKERRQ(ierr);
  ierr = VecGetSize(f,&M); CHKERRQ(ierr);
  ierr = inerpolation_matrix_mult_generic(mat,x,f,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

static PetscErrorCode inerpolation_matrix_mult_transpose(Mat mat,Vec x,Vec f) {
  PetscErrorCode ierr;
  PetscFunctionBegin;
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

  MPI_Comm comm;
  ierr = PetscObjectGetComm((PetscObject)dm1,&comm); CHKERRQ(ierr);

  int m,n,M,N;

  DM dm_down = dm1;
  DM dm_up = dm2;

  int dm_down_leveldown = dm_down->leveldown;
  int dm_up_leveldown = dm_up->leveldown;

  PetscInfo2(
    dm1,"Create interpolation DMMGViaApproxOrders dm1_leveldown = %d dm2_leveldown = %d\n",
    dm_down_leveldown,dm_up_leveldown
  );

  MGShellProjectionMatrix *mat_ctx = new MGShellProjectionMatrix();
  {
    // Coarser mesh
    GET_DM_FIELD(dm_down);
    if(dm_field->coarseningIS.size()<dm_down_leveldown) {
      SETERRQ(PETSC_COMM_WORLD,MOFEM_DATA_INCONSISTENCY,"data inconsistency");
    }
    mat_ctx->isDown = dm_field->coarseningIS[dm_field->coarseningIS.size()-1-dm_down_leveldown];
    mat_ctx->levelDown = dm_down_leveldown;
    ierr = ISGetSize(mat_ctx->isDown,&M); CHKERRQ(ierr);
    ierr = ISGetLocalSize(mat_ctx->isDown,&m); CHKERRQ(ierr);
  }
  {
    // Finer mesh
    GET_DM_FIELD(dm_up);
    if(dm_field->coarseningIS.size()<dm_up_leveldown) {
      SETERRQ(PETSC_COMM_WORLD,MOFEM_DATA_INCONSISTENCY,"data inconsistency");
    }
    mat_ctx->isUp = dm_field->coarseningIS[dm_field->coarseningIS.size()-1-dm_up_leveldown];
    mat_ctx->levelUp = dm_up_leveldown;
    ierr = ISGetSize(mat_ctx->isUp,&N); CHKERRQ(ierr);
    ierr = ISGetLocalSize(mat_ctx->isUp,&n); CHKERRQ(ierr);
  }

  ierr = MatCreateShell(comm,m,n,M,N,(void*)mat_ctx,mat); CHKERRQ(ierr);

  ierr = MatShellSetOperation(*mat,MATOP_DESTROY,(void(*)(void))inerpolation_matrix_destroy); CHKERRQ(ierr);
  ierr = MatShellSetOperation(*mat,MATOP_MULT,(void(*)(void))inerpolation_matrix_mult); CHKERRQ(ierr);
  ierr = MatShellSetOperation(*mat,MATOP_MULT_TRANSPOSE,(void(*)(void))inerpolation_matrix_mult_transpose); CHKERRQ(ierr);
  ierr = MatShellSetOperation(*mat,MATOP_MULT_ADD,(void(*)(void))inerpolation_matrix_mult_add); CHKERRQ(ierr);
  ierr = MatShellSetOperation(*mat,MATOP_MULT_TRANSPOSE_ADD,(void(*)(void))inerpolation_matrix_mult_transpose_add); CHKERRQ(ierr);

  *vec = PETSC_NULL;

  PetscFunctionReturn(0);
}

PetscErrorCode DMCreateGlobalVector_MGViaApproxOrders(DM dm,Vec *g) {
  PetscErrorCode ierr;
  PetscValidHeaderSpecific(dm,DM_CLASSID,1);
  PetscFunctionBegin;
  int leveldown = dm->leveldown;
  GET_DM_FIELD(dm);
  if(dm_field->kspOperators.empty()) {
    ierr = DMCreateGlobalVector_MoFEM(dm,g); CHKERRQ(ierr);
  } else {
    #if PETSC_VERSION_GE(3,5,3)
    ierr = MatCreateVecs(dm_field->kspOperators[dm_field->kspOperators.size()-1-leveldown],g,NULL); CHKERRQ(ierr);
    #else
    ierr = MatGetVecs(dm_field->kspOperators[dm_field->kspOperators.size()-1-leveldown],g,NULL); CHKERRQ(ierr);
    #endif
  }
  PetscInfo1(
    dm,"Create global vector DMMGViaApproxOrders leveldown = %d\n",dm->leveldown
  );
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
  const MoFEM::Interface *m_field_ptr;
  ierr = DMoFEMGetInterfacePtr(dM,&m_field_ptr); CHKERRQ(ierr);
  const MoFEMProblem *problem_ptr;
  ierr = DMMoFEMGetProblemPtr(dM,&problem_ptr); CHKERRQ(ierr);
  string problem_name = problem_ptr->getName();
  ierr = const_cast<MoFEM::Interface *>(m_field_ptr)
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

  std::vector<IS> is_vec(nbLevels+1);
  std::vector<int> is_glob_size(nbLevels+1),is_loc_size(nbLevels+1);

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

  for(int kk = 0;kk!=nbLevels;kk++) {
    Mat subA;
    ierr = DMMGViaApproxOrdersPushBackCoarseningIS(dM,is_vec[kk],A,&subA,true); CHKERRQ(ierr);
    if(subA) {
      ierr = MatDestroy(&subA); CHKERRQ(ierr);
    }
  }

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
