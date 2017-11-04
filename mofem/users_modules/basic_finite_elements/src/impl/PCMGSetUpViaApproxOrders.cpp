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

PCMGSubMatrixCtx::PCMGSubMatrixCtx(Mat a,IS is):
  A(a),
  iS(is) {
  // Increase reference of petsc object (works like shared_ptr but unique for PETSc)
  ierr = PetscObjectReference((PetscObject)A); CHKERRABORT(PETSC_COMM_WORLD,ierr);
  ierr = PetscObjectReference((PetscObject)iS); CHKERRABORT(PETSC_COMM_WORLD,ierr);
}

PCMGSubMatrixCtx::~PCMGSubMatrixCtx() {

  ierr = MatDestroy(&A); CHKERRABORT(PETSC_COMM_WORLD,ierr);
  ierr = ISDestroy(&iS); CHKERRABORT(PETSC_COMM_WORLD,ierr);
}

struct PCMGSubMatrixCtx_private: public PCMGSubMatrixCtx {
  PCMGSubMatrixCtx_private(Mat a,IS is):
  PCMGSubMatrixCtx(a,is),
  isInitisalised(false) {
    PetscLogEventRegister("PCMGSubMatrixCtx_mult",0,&MOFEM_EVENT_mult);
    PetscLogEventRegister("PCMGSubMatrixCtx_sor",0,&MOFEM_EVENT_sor);
  }
  ~PCMGSubMatrixCtx_private() {
    if(isInitisalised) {
      ierr = VecScatterDestroy(&sCat); CHKERRABORT(PETSC_COMM_WORLD,ierr);
      ierr = VecDestroy(&X); CHKERRABORT(PETSC_COMM_WORLD,ierr);
      ierr = VecDestroy(&F); CHKERRABORT(PETSC_COMM_WORLD,ierr);
    }
  }
  template<InsertMode MODE>
  friend MoFEMErrorCode sub_mat_mult_generic(Mat a,Vec x,Vec f);
  friend MoFEMErrorCode sub_mat_sor(
    Mat mat,Vec b,PetscReal omega,MatSORType flag,PetscReal shift,PetscInt its,PetscInt lits,Vec x
  );
public:
  MoFEMErrorCode initData(Vec x) {

    MoFEMFunctionBeginHot;
    if(!isInitisalised) {
      ierr = MatCreateVecs(A,&X,&F); CHKERRG(ierr);
      ierr = VecScatterCreate(X,iS,x,PETSC_NULL,&sCat); CHKERRG(ierr);
      isInitisalised = true;
    }
    MoFEMFunctionReturnHot(0);
  }
  PetscLogEvent MOFEM_EVENT_mult;
  PetscLogEvent MOFEM_EVENT_sor;
  bool isInitisalised;
};

template<InsertMode MODE>
MoFEMErrorCode sub_mat_mult_generic(Mat a,Vec x,Vec f) {
  void *void_ctx;

  MoFEMFunctionBeginHot;
  ierr = MatShellGetContext(a,&void_ctx); CHKERRG(ierr);
  PCMGSubMatrixCtx_private *ctx = (PCMGSubMatrixCtx_private*)void_ctx;
  if(!ctx->isInitisalised) {
    ierr = ctx->initData(x); CHKERRG(ierr);
  }
  PetscLogEventBegin(ctx->MOFEM_EVENT_mult,0,0,0,0);
  ierr = VecScatterBegin(ctx->sCat,x,ctx->X,INSERT_VALUES,SCATTER_REVERSE); CHKERRG(ierr);
  ierr = VecScatterEnd(ctx->sCat,x,ctx->X,INSERT_VALUES,SCATTER_REVERSE); CHKERRG(ierr);
  ierr = MatMult(ctx->A,ctx->X,ctx->F); CHKERRG(ierr);
  ierr = VecScatterBegin(ctx->sCat,ctx->F,f,MODE,SCATTER_FORWARD); CHKERRG(ierr);
  ierr = VecScatterEnd(ctx->sCat,ctx->F,f,MODE,SCATTER_FORWARD); CHKERRG(ierr);
  PetscLogEventEnd(ctx->MOFEM_EVENT_mult,0,0,0,0);
  MoFEMFunctionReturnHot(0);
}

MoFEMErrorCode sub_mat_mult(Mat a,Vec x,Vec f) {
  return sub_mat_mult_generic<INSERT_VALUES>(a,x,f);
}

MoFEMErrorCode sub_mat_mult_add(Mat a,Vec x,Vec f) {
  return sub_mat_mult_generic<ADD_VALUES>(a,x,f);
}

MoFEMErrorCode sub_mat_sor(
  Mat mat,Vec b,PetscReal omega,MatSORType flag,PetscReal shift,PetscInt its,PetscInt lits,Vec x
) {
  void *void_ctx;

  MoFEMFunctionBeginHot;
  ierr = MatShellGetContext(mat,&void_ctx); CHKERRG(ierr);
  PCMGSubMatrixCtx_private *ctx = (PCMGSubMatrixCtx_private*)void_ctx;
  if(!ctx->isInitisalised) {
    ierr = ctx->initData(x); CHKERRG(ierr);
  }
  PetscLogEventBegin(ctx->MOFEM_EVENT_sor,0,0,0,0);
  ierr = VecScatterBegin(ctx->sCat,b,ctx->X,INSERT_VALUES,SCATTER_REVERSE); CHKERRG(ierr);
  ierr = VecScatterEnd(ctx->sCat,b,ctx->X,INSERT_VALUES,SCATTER_REVERSE); CHKERRG(ierr);
  ierr = MatSOR(ctx->A,ctx->X,omega,flag,shift,its,lits,ctx->F); CHKERRG(ierr);
  ierr = VecScatterBegin(ctx->sCat,ctx->F,x,INSERT_VALUES,SCATTER_FORWARD); CHKERRG(ierr);
  ierr = VecScatterEnd(ctx->sCat,ctx->F,x,INSERT_VALUES,SCATTER_FORWARD); CHKERRG(ierr);
  PetscLogEventEnd(ctx->MOFEM_EVENT_sor,0,0,0,0);
  MoFEMFunctionReturnHot(0);
}

DMMGViaApproxOrdersCtx::DMMGViaApproxOrdersCtx():
  MoFEM::DMCtx(),
  aO(PETSC_NULL) {
    // std::cerr << "create dm\n";
}
DMMGViaApproxOrdersCtx::~DMMGViaApproxOrdersCtx() {
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
  // std::cerr << "destroy dm data\n";
}

MoFEMErrorCode DMMGViaApproxOrdersCtx::query_interface(const MOFEMuuid& uuid,MoFEM::UnknownInterface** iface) const {
  MoFEMFunctionBeginHot;
  *iface = NULL;
  if(uuid == IDD_DMMGVIAAPPROXORDERSCTX) {
    *iface = static_cast<DMMGViaApproxOrdersCtx*>(const_cast<DMMGViaApproxOrdersCtx*>(this));
    MoFEMFunctionReturnHot(0);
  } else {
    SETERRQ(PETSC_COMM_WORLD,MOFEM_DATA_INCONSISTENCY,"wrong interference");
  }

  ierr = DMCtx::query_interface(uuid,iface); CHKERRG(ierr);
  MoFEMFunctionReturnHot(0);
}

#define GET_DM_FIELD(DM) \
  MoFEM::UnknownInterface *iface; \
  ierr = ((DMCtx*)DM->data)->query_interface(IDD_DMMGVIAAPPROXORDERSCTX,&iface); CHKERRG(ierr); \
  DMMGViaApproxOrdersCtx *dm_field = static_cast<DMMGViaApproxOrdersCtx*>(iface)


MoFEMErrorCode DMMGViaApproxOrdersGetCtx(DM dm,DMMGViaApproxOrdersCtx **ctx) {
  PetscValidHeaderSpecific(dm,DM_CLASSID,1);
  MoFEMFunctionBeginHot;
  GET_DM_FIELD(dm);
  *ctx = dm_field;
  MoFEMFunctionReturnHot(0);
}

MoFEMErrorCode DMMGViaApproxOrdersSetAO(DM dm,AO ao) {
  PetscValidHeaderSpecific(dm,DM_CLASSID,1);
  MoFEMFunctionBeginHot;
  GET_DM_FIELD(dm);
  if(dm_field->aO) {
    //std::cerr << dm_field->aO << std::endl;
    ierr = AODestroy(&dm_field->aO); CHKERRG(ierr);
    // std::cerr << "destroy ao when adding\n";
  }
  dm_field->aO = ao;
  ierr = PetscObjectReference((PetscObject)ao); CHKERRG(ierr);
  // std::cerr << "add ao\n";
  MoFEMFunctionReturnHot(0);
}

MoFEMErrorCode DMMGViaApproxOrdersGetCoarseningISSize(DM dm,int *size) {
  PetscValidHeaderSpecific(dm,DM_CLASSID,1);
  MoFEMFunctionBeginHot;
  GET_DM_FIELD(dm);
  *size = dm_field->coarseningIS.size();
  MoFEMFunctionReturnHot(0);
}

MoFEMErrorCode DMMGViaApproxOrdersPushBackCoarseningIS(
  DM dm,IS is,Mat A,Mat *subA,bool create_sub_matrix,bool shell_sub_a
) {
  PetscValidHeaderSpecific(dm,DM_CLASSID,1);
  MoFEMFunctionBeginHot;
  GET_DM_FIELD(dm);
  dm_field->coarseningIS.push_back(is);
  dm_field->shellMatrixCtxPtr.push_back(new PCMGSubMatrixCtx_private(A,is));
  if(is) {
    ierr = PetscObjectReference((PetscObject)is); CHKERRG(ierr);
  }
  if(is) {
    IS is2 = is;
    if(dm_field->aO) {
      ierr = ISDuplicate(is,&is2); CHKERRG(ierr);
      ierr = ISCopy(is,is2); CHKERRG(ierr);
      ierr = AOApplicationToPetscIS(dm_field->aO,is2); CHKERRG(ierr);
    }
    if(create_sub_matrix) {
      if(shell_sub_a) {
        int n,N;
        ierr = ISGetSize(is,&N); CHKERRG(ierr);
        ierr = ISGetLocalSize(is,&n); CHKERRG(ierr);
        MPI_Comm comm;
        ierr = PetscObjectGetComm((PetscObject)A,&comm); CHKERRG(ierr);
        ierr = MatCreateShell(
          comm,n,n,N,N,&(dm_field->shellMatrixCtxPtr.back()),subA
        ); CHKERRG(ierr);
        ierr = MatShellSetOperation(
          *subA,MATOP_MULT,(void(*)(void))sub_mat_mult
        ); CHKERRG(ierr);
        ierr = MatShellSetOperation(
          *subA,MATOP_MULT_ADD,(void(*)(void))sub_mat_mult_add
        ); CHKERRG(ierr);
        ierr = MatShellSetOperation(
          *subA,MATOP_SOR,(void(*)(void))sub_mat_sor
        ); CHKERRG(ierr);
      } else {
        ierr = MatGetSubMatrix(A,is2,is2,MAT_INITIAL_MATRIX,subA); CHKERRG(ierr);
      }
    }
    if(dm_field->aO) {
      ierr = ISDestroy(&is2); CHKERRG(ierr);
    }
    dm_field->kspOperators.push_back(*subA);
    ierr = PetscObjectReference((PetscObject)(*subA)); CHKERRG(ierr);
  } else {
    SETERRQ(PETSC_COMM_WORLD,MOFEM_DATA_INCONSISTENCY,"data inconsistency");
  }
  PetscInfo(dm,"Push back IS to DMMGViaApproxOrders\n");
  MoFEMFunctionReturnHot(0);
}

MoFEMErrorCode DMMGViaApproxOrdersPopBackCoarseningIS(DM dm) {
  PetscValidHeaderSpecific(dm,DM_CLASSID,1);
  MoFEMFunctionBeginHot;
  GET_DM_FIELD(dm);
  if(dm_field->coarseningIS.back()) {
    ierr = ISDestroy(&dm_field->coarseningIS.back()); CHKERRG(ierr);
    dm_field->coarseningIS.pop_back();
  }
  if(dm_field->kspOperators.back()) {
    ierr = MatDestroy(&dm_field->kspOperators.back()); CHKERRG(ierr);
  }
  dm_field->kspOperators.pop_back();
  PetscInfo(dm,"Pop back IS to DMMGViaApproxOrders\n");
  MoFEMFunctionReturnHot(0);
}

MoFEMErrorCode DMMGViaApproxOrdersReplaceCoarseningIS(DM dm,IS *is_vec,int nb_elems,Mat A,int verb) {
  PetscValidHeaderSpecific(dm,DM_CLASSID,1);
  MoFEMFunctionBeginHot;
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
      ierr = ISEqual(*it,is_vec[ii],&flg); CHKERRG(ierr);
      if(!flg) {
        ierr = ISDestroy(&*it); CHKERRG(ierr);
        ierr = MatDestroy(&dm_field->kspOperators[ii]); CHKERRG(ierr);
        *it = is_vec[ii];
        ierr = PetscObjectReference((PetscObject)is_vec[ii]); CHKERRG(ierr);
        if(ii<nb_elems-1) {
          IS is = is_vec[ii];
          if(dm_field->aO) {
            ierr = ISDuplicate(is_vec[ii],&is); CHKERRG(ierr);
            ierr = ISCopy(is_vec[ii],is); CHKERRG(ierr);
            ierr = AOApplicationToPetscIS(dm_field->aO,is); CHKERRG(ierr);
          }
          Mat subA;
          ierr = MatGetSubMatrix(A,is,is,MAT_INITIAL_MATRIX,&subA); CHKERRG(ierr);
          ierr = PetscObjectReference((PetscObject)subA); CHKERRG(ierr);
          dm_field->kspOperators[ii] = subA;
          ierr = MatDestroy(&subA); CHKERRG(ierr);
          if(dm_field->aO) {
            ierr = ISDestroy(&is); CHKERRG(ierr);
          }
        } else {
          ierr = PetscObjectReference((PetscObject)A); CHKERRG(ierr);
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
      ierr = DMMGViaApproxOrdersPushBackCoarseningIS(dm,is_vec[ii],A,&subA,true,false); CHKERRG(ierr);
      ierr = MatDestroy(&subA); CHKERRG(ierr);
      nb_added++;
    }
    ierr = DMMGViaApproxOrdersPushBackCoarseningIS(dm,is_vec[ii],A,&A,false,false); CHKERRG(ierr);
    nb_added++;
  } else {
    for(;ii<dm_field->coarseningIS.size();ii++) {
      ierr = DMMGViaApproxOrdersPopBackCoarseningIS(dm); CHKERRG(ierr);
      nb_deleted++;
    }
  }
  MPI_Comm comm;
  ierr = PetscObjectGetComm((PetscObject)dm,&comm); CHKERRG(ierr);
  if(verb>0) {
    PetscPrintf(
      comm,"DMMGViaApproxOrders nb_no_changed = %d, nb_replaced = %d, nb_added = %d, nb_deleted = %d, size = %d\n",
      nb_no_changed,nb_replaced,nb_added,nb_deleted,dm_field->coarseningIS.size()
    );
  }
  PetscInfo(dm,"Replace IS to DMMGViaApproxOrders\n");
  MoFEMFunctionReturnHot(0);
}

MoFEMErrorCode DMMGViaApproxOrdersGetCtx(DM dm,const DMMGViaApproxOrdersCtx **ctx) {
  PetscValidHeaderSpecific(dm,DM_CLASSID,1);
  MoFEMFunctionBeginHot;
  GET_DM_FIELD(dm);
  *ctx = dm_field;
  MoFEMFunctionReturnHot(0);
}

MoFEMErrorCode DMRegister_MGViaApproxOrders(const char sname[]) {
  MoFEMFunctionBeginHot;
  ierr = DMRegister(sname,DMCreate_MGViaApproxOrders); CHKERRG(ierr);
  MoFEMFunctionReturnHot(0);
}

// MoFEMErrorCode DMDestroy_MGViaApproxOrders(DM dm) {
//   PetscValidHeaderSpecific(dm,DM_CLASSID,1);
//   MoFEMFunctionBeginHot;
//   if(!((DMMGViaApproxOrdersCtx*)dm->data)->referenceNumber) {
//     DMMGViaApproxOrdersCtx *dm_field = (DMMGViaApproxOrdersCtx*)dm->data;
//     if(dm_field->destroyProblem) {
//       if(dm_field->mField_ptr->check_problem(dm_field->problemName)) {
//         dm_field->mField_ptr->delete_problem(dm_field->problemName);
//       } // else problem has to be deleted by the user
//     }
//     cerr << "Destroy " << dm_field->problemName << endl;
//     delete (DMMGViaApproxOrdersCtx*)dm->data;
//   } else {
//     DMMGViaApproxOrdersCtx *dm_field = (DMMGViaApproxOrdersCtx*)dm->data;
//
//     cerr << "Derefrence " << dm_field->problemName << " "  << ((DMCtx*)dm->data)->referenceNumber << endl;
//     (((DMMGViaApproxOrdersCtx*)dm->data)->referenceNumber)--;
//   }
//   MoFEMFunctionReturnHot(0);
// }

static MoFEMErrorCode ksp_set_operators(KSP ksp,Mat A,Mat B,void *ctx) {
  MoFEMFunctionBeginHot;
  MoFEMFunctionReturnHot(0);
}

MoFEMErrorCode DMCreate_MGViaApproxOrders(DM dm) {
  PetscValidHeaderSpecific(dm,DM_CLASSID,1);
  MoFEMFunctionBeginHot;
  if(!dm->data) {
    dm->data = new DMMGViaApproxOrdersCtx();
  } else {
    ((DMCtx*)(dm->data))->referenceNumber++;
  }
  // cerr << "Create " << ((DMCtx*)(dm->data))->referenceNumber << endl;
  ierr = DMSetOperators_MoFEM(dm); CHKERRG(ierr);
  dm->ops->creatematrix = DMCreateMatrix_MGViaApproxOrders;
  dm->ops->createglobalvector = DMCreateGlobalVector_MGViaApproxOrders;
  dm->ops->coarsen = DMCoarsen_MGViaApproxOrders;
  // dm->ops->destroy = DMDestroy_MGViaApproxOrders;
  dm->ops->createinterpolation = DMCreateInterpolation_MGViaApproxOrders;
  ierr = DMKSPSetComputeOperators(dm,ksp_set_operators,NULL); CHKERRG(ierr);
  PetscInfo1(dm,"Create DMMGViaApproxOrders reference = %d\n",((DMCtx*)(dm->data))->referenceNumber);
  MoFEMFunctionReturnHot(0);
}

MoFEMErrorCode DMCreateMatrix_MGViaApproxOrders(DM dm,Mat *M) {

  PetscValidHeaderSpecific(dm,DM_CLASSID,1);
  MoFEMFunctionBeginHot;
  GET_DM_FIELD(dm);

  int leveldown = dm->leveldown;

  if(dm_field->kspOperators.empty()) {
    ierr = DMCreateMatrix_MoFEM(dm,M); CHKERRG(ierr);
  } else {
    MPI_Comm comm;
    ierr = PetscObjectGetComm((PetscObject)dm,&comm); CHKERRG(ierr);
    if(dm_field->kspOperators.empty()) {
      SETERRQ(comm,MOFEM_DATA_INCONSISTENCY,"data inconsistency, operator can not be set");
    }
    if(dm_field->kspOperators.size()<leveldown) {
      SETERRQ(comm,MOFEM_DATA_INCONSISTENCY,"data inconsistency, no IS for that level");
    }
    *M = dm_field->kspOperators[dm_field->kspOperators.size()-1-leveldown];
    ierr = PetscObjectReference((PetscObject)*M); CHKERRG(ierr);
  }

  PetscInfo1(dm,"Create Matrix DMMGViaApproxOrders leveldown = %d\n",leveldown);

  MoFEMFunctionReturnHot(0);
}

MoFEMErrorCode DMCoarsen_MGViaApproxOrders(DM dm, MPI_Comm comm, DM *dmc) {
  PetscValidHeaderSpecific(dm,DM_CLASSID,1);
  MoFEMFunctionBeginHot;
  GET_DM_FIELD(dm);
  ierr = PetscObjectGetComm((PetscObject)dm,&comm); CHKERRG(ierr);
  ierr = DMCreate(comm,dmc);CHKERRG(ierr);
  (*dmc)->data = dm->data;
  DMType type;
  ierr = DMGetType(dm,&type); CHKERRG(ierr);
  ierr = DMSetType(*dmc,type); CHKERRG(ierr);
  ierr = PetscObjectReference((PetscObject)(*dmc));CHKERRG(ierr);
  PetscInfo1(dm,"Coarsen DMMGViaApproxOrders leveldown = %d\n",dm->leveldown);
  MoFEMFunctionReturnHot(0);
}

MoFEMErrorCode DMCreateInterpolation_MGViaApproxOrders(DM dm1,DM dm2,Mat *mat,Vec *vec) {
  PetscValidHeaderSpecific(dm1,DM_CLASSID,1);
  PetscValidHeaderSpecific(dm2,DM_CLASSID,1);
  MoFEMFunctionBeginHot;

  MPI_Comm comm;
  ierr = PetscObjectGetComm((PetscObject)dm1,&comm); CHKERRG(ierr);

  int m,n,M,N;

  DM dm_down = dm1;
  DM dm_up = dm2;

  int dm_down_leveldown = dm_down->leveldown;
  int dm_up_leveldown = dm_up->leveldown;

  PetscInfo2(
    dm1,"Create interpolation DMMGViaApproxOrders dm1_leveldown = %d dm2_leveldown = %d\n",
    dm_down_leveldown,dm_up_leveldown
  );

  IS is_down,is_up;
  {
    // Coarser mesh
    GET_DM_FIELD(dm_down);
    if(dm_field->coarseningIS.size()<dm_down_leveldown) {
      SETERRQ(PETSC_COMM_WORLD,MOFEM_DATA_INCONSISTENCY,"data inconsistency");
    }
    is_down = dm_field->coarseningIS[dm_field->coarseningIS.size()-1-dm_down_leveldown];
    ierr = ISGetSize(is_down,&M); CHKERRG(ierr);
    ierr = ISGetLocalSize(is_down,&m); CHKERRG(ierr);
  }
  {
    // Finer mesh
    GET_DM_FIELD(dm_up);
    if(dm_field->coarseningIS.size()<dm_up_leveldown) {
      SETERRQ(PETSC_COMM_WORLD,MOFEM_DATA_INCONSISTENCY,"data inconsistency");
    }
    is_up = dm_field->coarseningIS[dm_field->coarseningIS.size()-1-dm_up_leveldown];
    ierr = ISGetSize(is_up,&N); CHKERRG(ierr);
    ierr = ISGetLocalSize(is_up,&n); CHKERRG(ierr);
  }

  // is_dow rows
  // is_up columns

  ierr = MatCreate(comm,mat); CHKERRG(ierr);
  ierr = MatSetSizes(*mat,m,n,M,N); CHKERRG(ierr);
  ierr = MatSetType(*mat,MATMPIAIJ); CHKERRG(ierr);
  ierr = MatMPIAIJSetPreallocation(*mat,1,PETSC_NULL,0,PETSC_NULL); CHKERRG(ierr);

  //get matrix layout
  PetscLayout rmap,cmap;
  ierr = MatGetLayouts(*mat,&rmap,&cmap); CHKERRG(ierr);
  int rstart,rend,cstart,cend;
  ierr = PetscLayoutGetRange(rmap,&rstart,&rend); CHKERRG(ierr);
  ierr = PetscLayoutGetRange(cmap,&cstart,&cend); CHKERRG(ierr);

  // if(verb>0) {
  //   PetscSynchronizedPrintf(comm,"level %d row start %d row end %d\n",kk,rstart,rend);
  //   PetscSynchronizedPrintf(comm,"level %d col start %d col end %d\n",kk,cstart,cend);
  // }

  const int *row_indices_ptr,*col_indices_ptr;
  ierr = ISGetIndices(is_down,&row_indices_ptr); CHKERRG(ierr);
  ierr = ISGetIndices(is_up,&col_indices_ptr); CHKERRG(ierr);

  map<int,int> idx_map;
  for(int ii = 0;ii<m;ii++) {
    idx_map[row_indices_ptr[ii]] = rstart+ii;
  }

  ierr = MatZeroEntries(*mat); CHKERRG(ierr);
  // FIXME: Use MatCreateMPIAIJWithArrays and set array directly
  for(int jj = 0;jj<n;jj++) {
    map<int,int>::iterator mit = idx_map.find(col_indices_ptr[jj]);
    if(mit != idx_map.end()) {
      ierr = MatSetValue(*mat,mit->second,cstart+jj,1,INSERT_VALUES); CHKERRG(ierr);
    }
  }

  ierr = ISRestoreIndices(is_down,&row_indices_ptr); CHKERRG(ierr);
  ierr = ISRestoreIndices(is_up,&col_indices_ptr); CHKERRG(ierr);

  ierr = MatAssemblyBegin(*mat,MAT_FINAL_ASSEMBLY); CHKERRG(ierr);
  ierr = MatAssemblyEnd(*mat,MAT_FINAL_ASSEMBLY); CHKERRG(ierr);

  if(vec!=NULL) {
    *vec = PETSC_NULL;
  }

  MoFEMFunctionReturnHot(0);
}

MoFEMErrorCode DMCreateGlobalVector_MGViaApproxOrders(DM dm,Vec *g) {

  PetscValidHeaderSpecific(dm,DM_CLASSID,1);
  MoFEMFunctionBeginHot;
  int leveldown = dm->leveldown;
  GET_DM_FIELD(dm);
  if(dm_field->kspOperators.empty()) {
    ierr = DMCreateGlobalVector_MoFEM(dm,g); CHKERRG(ierr);
  } else {
    #if PETSC_VERSION_GE(3,5,3)
    ierr = MatCreateVecs(dm_field->kspOperators[dm_field->kspOperators.size()-1-leveldown],g,NULL); CHKERRG(ierr);
    #else
    ierr = MatGetVecs(dm_field->kspOperators[dm_field->kspOperators.size()-1-leveldown],g,NULL); CHKERRG(ierr);
    #endif
  }
  PetscInfo1(
    dm,"Create global vector DMMGViaApproxOrders leveldown = %d\n",dm->leveldown
  );
  MoFEMFunctionReturnHot(0);
}

MoFEMErrorCode PCMGSetUpViaApproxOrdersCtx::getOptions() {
  MoFEMFunctionBeginHot;
  ierr = PetscOptionsBegin(
    PETSC_COMM_WORLD,"",
    "MOFEM Multi-Grid (Orders) pre-conditioner","none"
  ); CHKERRG(ierr);

  ierr = PetscOptionsInt("-mofem_mg_levels",
    "nb levels of multi-grid solver","",
    2,&nbLevels,PETSC_NULL
  ); CHKERRG(ierr);

  ierr = PetscOptionsInt(
    "-mofem_mg_coarse_order",
    "approximation order of coarse level","",
    2,&coarseOrder,PETSC_NULL
  ); CHKERRG(ierr);

  ierr = PetscOptionsInt(
    "-mofem_mg_order_at_last_level",
    "order at last level","",
    100,&orderAtLastLevel,PETSC_NULL
  ); CHKERRG(ierr);


  ierr = PetscOptionsInt(
    "-mofem_mg_verbose",
    "nb levels of multi-grid solver","",
    0,&verboseLevel,PETSC_NULL
  ); CHKERRG(ierr);

  PetscBool shell_sub_a = shellSubA ? PETSC_TRUE : PETSC_FALSE;
  ierr = PetscOptionsBool(
    "-mofem_mg_shell_a",
    "use shell matrix as sub matrix","",
    shell_sub_a,&shell_sub_a,NULL
  ); CHKERRG(ierr);
  shellSubA = (shellSubA == PETSC_TRUE);


  ierr = PetscOptionsEnd(); CHKERRG(ierr);
  MoFEMFunctionReturnHot(0);
}

MoFEMErrorCode PCMGSetUpViaApproxOrdersCtx::createIsAtLevel(int kk,IS *is) {
  MoFEM::Interface *m_field_ptr;
  MoFEM::ISManager *is_manager_ptr;
  MoFEMFunctionBeginHot;
  //if is last level, take all remaining orders dofs, if any left
  ierr = DMoFEMGetInterfacePtr(dM,&m_field_ptr); CHKERRG(ierr);
  ierr = m_field_ptr->getInterface(is_manager_ptr); CHKERRG(ierr);
  const Problem *problem_ptr;
  ierr = DMMoFEMGetProblemPtr(dM,&problem_ptr); CHKERRG(ierr);
  int order_at_next_level = kk+coarseOrder;
  if(kk == nbLevels-1) {
    int first = problem_ptr->getNumeredDofsRows()->
    get<PetscLocalIdx_mi_tag>().find(0)->get()->getPetscGlobalDofIdx();
    ierr = ISCreateStride(
      PETSC_COMM_WORLD,problem_ptr->getNbLocalDofsRow(),first,1,is
    ); CHKERRG(ierr);
    MoFEMFunctionReturnHot(0);
    // order_at_next_level = orderAtLastLevel;
  }
  string problem_name = problem_ptr->getName();
  ierr = is_manager_ptr->isCreateProblemOrder(
    problem_name,ROW,0,order_at_next_level,is
  ); CHKERRG(ierr);
  MoFEMFunctionReturnHot(0);
}

MoFEMErrorCode PCMGSetUpViaApproxOrdersCtx::destroyIsAtLevel(int kk,IS *is) {
  MoFEMFunctionBeginHot;
  ierr = ISDestroy(is);
  MoFEMFunctionReturnHot(0);
}

MoFEMErrorCode PCMGSetUpViaApproxOrdersCtx::buildProlongationOperator(
  bool use_mat_a,int verb
) {
  MoFEMFunctionBeginHot;
  verb = verb > verboseLevel ? verb : verboseLevel;

  MPI_Comm comm;
  ierr = PetscObjectGetComm((PetscObject)dM,&comm); CHKERRG(ierr);

  if(verb>0) {
    PetscPrintf(comm,"set MG levels %u\n",nbLevels);
  }

  std::vector<IS> is_vec(nbLevels+1);
  std::vector<int> is_glob_size(nbLevels+1),is_loc_size(nbLevels+1);

  for(int kk = 0;kk<nbLevels;kk++) {

    //get indices up to up to give approximation order
    ierr = createIsAtLevel(kk,&is_vec[kk]); CHKERRG(ierr);
    ierr = ISGetSize(is_vec[kk],&is_glob_size[kk]); CHKERRG(ierr);
    ierr = ISGetLocalSize(is_vec[kk],&is_loc_size[kk]); CHKERRG(ierr);

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

  for(int kk = 0;kk!=nbLevels;kk++) {
    Mat subA;
    if(kk==nbLevels-1&&use_mat_a) {
      subA = A;
      ierr = DMMGViaApproxOrdersPushBackCoarseningIS(dM,is_vec[kk],A,&subA,false,false); CHKERRG(ierr);
    } else {
      if(kk>0) {
        // Not coarse level
        ierr = DMMGViaApproxOrdersPushBackCoarseningIS(dM,is_vec[kk],A,&subA,true,shellSubA); CHKERRG(ierr);
      } else {
        // Coarse lave is compressed matrix allowing for factorization when needed
        ierr = DMMGViaApproxOrdersPushBackCoarseningIS(dM,is_vec[kk],A,&subA,true,false); CHKERRG(ierr);
      }
      if(subA) {
        ierr = MatDestroy(&subA); CHKERRG(ierr);
      }
    }
  }

  for(unsigned int kk = 0;kk<is_vec.size();kk++) {
    ierr = destroyIsAtLevel(kk,&is_vec[kk]); CHKERRG(ierr);
  }

  if(verb>0) {
    PetscSynchronizedFlush(comm,PETSC_STDOUT);
  }

  MoFEMFunctionReturnHot(0);
}

MoFEMErrorCode PCMGSetUpViaApproxOrders(PC pc, PCMGSetUpViaApproxOrdersCtx *ctx,
                                        int verb) {

  PetscValidHeaderSpecific(pc, PC_CLASSID, 1);
  MoFEMFunctionBegin;

  MPI_Comm comm;
  CHKERR PetscObjectGetComm((PetscObject)pc, &comm);
  if (verb > 0) {
    PetscPrintf(comm, "Start PCMGSetUpViaApproxOrders\n");
  }

  CHKERR ctx->getOptions();
  CHKERR ctx->buildProlongationOperator(true, verb);

#if PETSC_VERSION_GE(3, 8, 0)
  CHKERR PCMGSetGalerkin(pc, PC_MG_GALERKIN_BOTH);
#else
  CHKERR PCMGSetGalerkin(pc, PETSC_FALSE);
#endif

  CHKERR PCMGSetLevels(pc, ctx->nbLevels, NULL);

  if (verb > 0) {
    PetscPrintf(comm, "End PCMGSetUpViaApproxOrders\n");
  }

  MoFEMFunctionReturn(0);
}
