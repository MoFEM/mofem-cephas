/** \file PCMGSetUpViaApproxOrders.cpp
 * \brief implementation of multi-grid solver for p- adaptivity
 */


#undef PETSC_VERSION_RELEASE
#define PETSC_VERSION_RELEASE 1

#if PETSC_VERSION_GE(3, 6, 0)
#include <petsc/private/petscimpl.h>
#else
#include <petsc-private/petscimpl.h>
#endif

#if PETSC_VERSION_GE(3, 6, 0)
#include <petsc/private/dmimpl.h> /*I  "petscdm.h"   I*/
// #include <petsc/private/vecimpl.h> /*I  "petscdm.h"   I*/
#else
#include <petsc-private/dmimpl.h>  /*I  "petscdm.h"   I*/
#include <petsc-private/vecimpl.h> /*I  "petscdm.h"   I*/
#endif

#include <DMMoFEM.hpp>
#include <DMCtxImpl.hpp>
#include <PCMGSetUpViaApproxOrders.hpp>

namespace MoFEM {

struct PCMGSubMatrixCtx {

  PCMGSubMatrixCtx(Mat a, IS is);
  virtual ~PCMGSubMatrixCtx() = default;

  template <InsertMode MODE>
  friend MoFEMErrorCode sub_mat_mult_generic(Mat a, Vec x, Vec f);
  friend MoFEMErrorCode sub_mat_sor(Mat mat, Vec b, PetscReal omega,
                                    MatSORType flag, PetscReal shift,
                                    PetscInt its, PetscInt lits, Vec x);

  MoFEMErrorCode initData(Vec x);

  SmartPetscObj<Mat> A;
  SmartPetscObj<IS> iS;
  SmartPetscObj<VecScatter> sCat;
  SmartPetscObj<Vec> X;
  SmartPetscObj<Vec> F;

  PetscLogEvent MOFEM_EVENT_mult;
  PetscLogEvent MOFEM_EVENT_sor;
};

PCMGSubMatrixCtx::PCMGSubMatrixCtx(Mat a, IS is) : A(a, true), iS(is, true) {
  PetscLogEventRegister("PCMGSubMatrixCtx_mult", 0, &MOFEM_EVENT_mult);
  PetscLogEventRegister("PCMGSubMatrixCtx_sor", 0, &MOFEM_EVENT_sor);
}

MoFEMErrorCode PCMGSubMatrixCtx::initData(Vec x) {
  MoFEMFunctionBegin;
  if (sCat.use_count() == 0) {
    auto [x_tmp, f_tmp] = matCreateVecs(A);
    X = x_tmp;
    F = f_tmp;
    sCat = createVecScatter(X, iS, x, PETSC_NULL);
  }
  MoFEMFunctionReturn(0);
}

/**
 * \brief Structure for DM for multi-grid via approximation orders
 * \ingroup dm
 */
struct DMMGViaApproxOrdersCtx : public MoFEM::DMCtxImpl {

  MoFEMErrorCode query_interface(boost::typeindex::type_index type_index,
                                 MoFEM::UnknownInterface **iface) const;

  DMMGViaApproxOrdersCtx();
  virtual ~DMMGViaApproxOrdersCtx();

  MoFEMErrorCode destroyCoarseningIS();

  SmartPetscObj<DM> coarsenDM;
  SmartPetscObj<AO> aO;
  std::vector<SmartPetscObj<IS>> coarseningIS;  ///< Coarsening IS
  std::vector<SmartPetscObj<Mat>> kspOperators; ///< Get KSP operators
  boost::ptr_vector<PCMGSubMatrixCtx>
      shellMatrixCtxPtr; ///< Shell sub-matrix context
};

/* \brief Set data structures of MG pre-conditioner via approximation orders
 */
struct PCMGSetUpViaApproxOrdersCtx {

  PCMGSetUpViaApproxOrdersCtx(DM dm, Mat a, bool shell_sub_a)
      : dM(dm), A(a), nbLevels(2), coarseOrder(2), orderAtLastLevel(1000),
        shellSubA(shell_sub_a), verboseLevel(0) {}

  virtual ~PCMGSetUpViaApproxOrdersCtx() = default;

  /**
   * \brief get options from line command
   * @return error code
   */
  virtual MoFEMErrorCode getOptions();

  /**
   * \brief Set IS for levels
   * @param  kk level
   * @param  is pointer to IS
   * @return    error code
   */
  virtual MoFEMErrorCode createIsAtLevel(int kk, IS *is);

  /**
   * \brief Destroy IS if internally created
   * @param  kk level
   * @param  is pointer to is
   * @return    error code
   */
  virtual MoFEMErrorCode destroyIsAtLevel(int kk, IS *is);

  /**
   * \brief Set up data structures for MG
   * @param  pc   MG pre-conditioner
   * <http://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/PC/PCMG.html>
   * @param  verb verbosity level
   * @return      error code
   */
  virtual MoFEMErrorCode buildProlongationOperator(bool use_mat_a,
                                                   int verb = 0);

  DM dM; ///< Distributed mesh manager
  Mat A; ///< Matrix at fine level

  int nbLevels;         ///< number of multi-grid levels
  int coarseOrder;      ///< approximation order of coarse level
  int orderAtLastLevel; ///< set maximal evaluated order

  bool shellSubA;
  int verboseLevel;
};

template <InsertMode MODE>
MoFEMErrorCode sub_mat_mult_generic(Mat a, Vec x, Vec f) {
  void *void_ctx;
  MoFEMFunctionBegin;
  CHKERR MatShellGetContext(a, &void_ctx);
  PCMGSubMatrixCtx *ctx = (PCMGSubMatrixCtx *)void_ctx;
  PetscLogEventBegin(ctx->MOFEM_EVENT_mult, 0, 0, 0, 0);
  CHKERR ctx->initData(x);
  CHKERR VecScatterBegin(ctx->sCat, x, ctx->X, INSERT_VALUES, SCATTER_REVERSE);
  CHKERR VecScatterEnd(ctx->sCat, x, ctx->X, INSERT_VALUES, SCATTER_REVERSE);
  CHKERR MatMult(ctx->A, ctx->X, ctx->F);
  CHKERR VecScatterBegin(ctx->sCat, ctx->F, f, MODE, SCATTER_FORWARD);
  CHKERR VecScatterEnd(ctx->sCat, ctx->F, f, MODE, SCATTER_FORWARD);
  PetscLogEventEnd(ctx->MOFEM_EVENT_mult, 0, 0, 0, 0);
  MoFEMFunctionReturn(0);
}

MoFEMErrorCode sub_mat_mult(Mat a, Vec x, Vec f) {
  return sub_mat_mult_generic<INSERT_VALUES>(a, x, f);
}

MoFEMErrorCode sub_mat_mult_add(Mat a, Vec x, Vec f) {
  return sub_mat_mult_generic<ADD_VALUES>(a, x, f);
}

MoFEMErrorCode sub_mat_sor(Mat mat, Vec b, PetscReal omega, MatSORType flag,
                           PetscReal shift, PetscInt its, PetscInt lits,
                           Vec x) {

  //FIXME: that is crap implementation of SOR

  void *void_ctx;
  MoFEMFunctionBegin;
  CHKERR MatShellGetContext(mat, &void_ctx);
  PCMGSubMatrixCtx *ctx = (PCMGSubMatrixCtx *)void_ctx;
  PetscLogEventBegin(ctx->MOFEM_EVENT_sor, 0, 0, 0, 0);
  CHKERR ctx->initData(x);
  CHKERR VecScatterBegin(ctx->sCat, b, ctx->X, INSERT_VALUES, SCATTER_REVERSE);
  CHKERR VecScatterEnd(ctx->sCat, b, ctx->X, INSERT_VALUES, SCATTER_REVERSE);
  CHKERR MatSOR(ctx->A, ctx->X, omega, flag, shift, its, lits, ctx->F);
  CHKERR VecScatterBegin(ctx->sCat, ctx->F, x, INSERT_VALUES, SCATTER_FORWARD);
  CHKERR VecScatterEnd(ctx->sCat, ctx->F, x, INSERT_VALUES, SCATTER_FORWARD);
  PetscLogEventEnd(ctx->MOFEM_EVENT_sor, 0, 0, 0, 0);
  MoFEMFunctionReturn(0);
}

DMMGViaApproxOrdersCtx::DMMGViaApproxOrdersCtx() : MoFEM::DMCtxImpl() {}

DMMGViaApproxOrdersCtx::~DMMGViaApproxOrdersCtx() {
  ierr = destroyCoarseningIS();
  CHKERRABORT(PETSC_COMM_WORLD, ierr);
}

MoFEMErrorCode DMMGViaApproxOrdersCtx::destroyCoarseningIS() {
  MoFEMFunctionBegin;
  coarsenDM.reset();
  aO.reset();
  coarseningIS.clear();
  kspOperators.clear();
  shellMatrixCtxPtr.clear();
  MoFEMFunctionReturn(0);
}

MoFEMErrorCode
DMMGViaApproxOrdersCtx::query_interface(boost::typeindex::type_index type_index,
                                        MoFEM::UnknownInterface **iface) const {
  *iface = static_cast<DMMGViaApproxOrdersCtx *>(
      const_cast<DMMGViaApproxOrdersCtx *>(this));
  return 0;
}

#define GET_DM_FIELD(DM)                                                       \
  auto dm_field =                                                              \
      static_cast<DMCtx *>(DM->data)->getInterface<DMMGViaApproxOrdersCtx>();  \
  NOT_USED(dm_field)

MoFEMErrorCode DMMGViaApproxOrdersSetAO(DM dm, AO ao) {
  PetscValidHeaderSpecific(dm, DM_CLASSID, 1);
  MoFEMFunctionBegin;
  GET_DM_FIELD(dm);
  dm_field->aO = SmartPetscObj<AO>(ao, true);
  MoFEMFunctionReturn(0);
}

MoFEMErrorCode DMMGViaApproxOrdersPushBackCoarseningIS(DM dm, IS is, Mat A,
                                                       bool create_sub_matrix,
                                                       bool shell_sub_a) {
  PetscValidHeaderSpecific(dm, DM_CLASSID, 1);
  MoFEMFunctionBegin;
  GET_DM_FIELD(dm);
  dm_field->coarseningIS.emplace_back(SmartPetscObj<IS>(is, true));
  dm_field->shellMatrixCtxPtr.push_back(new PCMGSubMatrixCtx(A, is));
  if (is) {
    auto is2 = SmartPetscObj<IS>(is, true);
    if (dm_field->aO) {
      is2 = isDuplicate(is);
      CHKERR ISCopy(is, is2);
      CHKERR AOApplicationToPetscIS(dm_field->aO, is2);
    }
    if (create_sub_matrix) {
      if (shell_sub_a) {
        int n, N;
        CHKERR ISGetSize(is, &N);
        CHKERR ISGetLocalSize(is, &n);
        MPI_Comm comm;
        CHKERR PetscObjectGetComm((PetscObject)A, &comm);
        Mat sub_a_raw;
        CHKERR MatCreateShell(comm, n, n, N, N,
                              &(dm_field->shellMatrixCtxPtr.back()),
                              &sub_a_raw);
        CHKERR MatShellSetOperation(sub_a_raw, MATOP_MULT,
                                    (void (*)(void))sub_mat_mult);
        CHKERR MatShellSetOperation(sub_a_raw, MATOP_MULT_ADD,
                                    (void (*)(void))sub_mat_mult_add);
        CHKERR MatShellSetOperation(sub_a_raw, MATOP_SOR,
                                    (void (*)(void))sub_mat_sor);
        dm_field->kspOperators.emplace_back(
            SmartPetscObj<Mat>(sub_a_raw, false));
      } else {
        Mat sub_a_raw;
#if PETSC_VERSION_GE(3, 8, 0)
        CHKERR MatCreateSubMatrix(A, is2, is2, MAT_INITIAL_MATRIX, &sub_a_raw);
#else
        CHKERR MatGetSubMatrix(A, is2, is2, MAT_INITIAL_MATRIX, &sub_a_raw);
#endif
        dm_field->kspOperators.emplace_back(
            SmartPetscObj<Mat>(sub_a_raw, false));
      }
    } else {
      dm_field->kspOperators.emplace_back(SmartPetscObj<Mat>(A, true));
    }
  } else {
    SETERRQ(PETSC_COMM_WORLD, MOFEM_DATA_INCONSISTENCY, "data inconsistency");
  }
  PetscInfo(dm, "Push back IS to DMMGViaApproxOrders\n");
  MoFEMFunctionReturn(0);
}

MoFEMErrorCode DMMGViaApproxOrdersPopBackCoarseningIS(DM dm) {
  PetscValidHeaderSpecific(dm, DM_CLASSID, 1);
  MoFEMFunctionBegin;
  GET_DM_FIELD(dm);
  if (dm_field->coarseningIS.back()) {
    dm_field->coarseningIS.pop_back();
  }
  dm_field->kspOperators.pop_back();
  PetscInfo(dm, "Pop back IS to DMMGViaApproxOrders\n");
  MoFEMFunctionReturn(0);
}

MoFEMErrorCode DMMGViaApproxOrdersClearCoarseningIS(DM dm) {
  PetscValidHeaderSpecific(dm, DM_CLASSID, 1);
  MoFEMFunctionBegin;
  GET_DM_FIELD(dm);
  CHKERR dm_field->destroyCoarseningIS();
  PetscInfo(dm, "Clear DMs data structures\n");
  MoFEMFunctionReturn(0);
}

MoFEMErrorCode DMRegister_MGViaApproxOrders(const char sname[]) {
  MoFEMFunctionBegin;
  CHKERR DMRegister(sname, DMCreate_MGViaApproxOrders);
  MoFEMFunctionReturn(0);
}

static MoFEMErrorCode ksp_set_operators(KSP ksp, Mat A, Mat B, void *ctx) {
  MoFEMFunctionBeginHot;
  MoFEMFunctionReturnHot(0);
}

MoFEMErrorCode DMCreate_MGViaApproxOrders(DM dm) {
  PetscValidHeaderSpecific(dm, DM_CLASSID, 1);
  MoFEMFunctionBegin;
  if (!dm->data) {
    dm->data = new DMMGViaApproxOrdersCtx();
  } else {
    ((DMCtxImpl *)(dm->data))->incrementReference();
  }
  CHKERR DMSetOperators_MoFEM(dm);

  dm->ops->creatematrix = DMCreateMatrix_MGViaApproxOrders;
  dm->ops->createglobalvector = DMCreateGlobalVector_MGViaApproxOrders;
  dm->ops->coarsen = DMCoarsen_MGViaApproxOrders;
  dm->ops->createinterpolation = DMCreateInterpolation_MGViaApproxOrders;
  dm->ops->destroy = DMDestroy_MGViaApproxOrders;

  CHKERR DMKSPSetComputeOperators(dm, ksp_set_operators, NULL);
  PetscInfo1(dm, "Create DMMGViaApproxOrders reference = %d\n",
             ((DMCtxImpl *)(dm->data))->useCount());
  MoFEMFunctionReturn(0);
}

PetscErrorCode DMDestroy_MGViaApproxOrders(DM dm) {
  PetscValidHeaderSpecific(dm, DM_CLASSID, 1);
  GET_DM_FIELD(dm);
  MoFEMFunctionBeginHot;
  CHKERR dm_field->destroyCoarseningIS();
  CHKERR DMDestroy_MoFEM(dm);
  MoFEMFunctionReturnHot(0);
}

MoFEMErrorCode DMCreateMatrix_MGViaApproxOrders(DM dm, Mat *M) {

  PetscValidHeaderSpecific(dm, DM_CLASSID, 1);
  MoFEMFunctionBegin;
  GET_DM_FIELD(dm);

  int leveldown = dm->leveldown;

  if (dm_field->kspOperators.empty()) {
    CHKERR DMCreateMatrix_MoFEM(dm, M);
  } else {
    MPI_Comm comm;
    CHKERR PetscObjectGetComm((PetscObject)dm, &comm);
    if (dm_field->kspOperators.empty()) {
      SETERRQ(comm, MOFEM_DATA_INCONSISTENCY,
              "data inconsistency, operator can not be set");
    }
    if (static_cast<int>(dm_field->kspOperators.size()) < leveldown) {
      SETERRQ(comm, MOFEM_DATA_INCONSISTENCY,
              "data inconsistency, no IS for that level");
    }
    *M = dm_field->kspOperators[dm_field->kspOperators.size() - 1 - leveldown];
    CHKERR PetscObjectReference((PetscObject)*M);
    CHKERR MatSetDM(*M, dm);
  }

  PetscInfo1(dm, "Create Matrix DMMGViaApproxOrders leveldown = %d\n",
             leveldown);

  MoFEMFunctionReturn(0);
}

MoFEMErrorCode DMCoarsen_MGViaApproxOrders(DM dm, MPI_Comm comm, DM *dmc) {
  PetscValidHeaderSpecific(dm, DM_CLASSID, 1);
  MoFEMFunctionBegin;
  GET_DM_FIELD(dm);
  CHKERR PetscObjectGetComm((PetscObject)dm, &comm);
  CHKERR DMCreate(comm, dmc);
  (*dmc)->data = dm->data;
  DMType type;
  CHKERR DMGetType(dm, &type);
  CHKERR DMSetType(*dmc, type);
  dm_field->coarsenDM = SmartPetscObj<DM>(*dmc, false);
  PetscInfo1(dm, "Coarsen DMMGViaApproxOrders leveldown = %d\n", dm->leveldown);
  MoFEMFunctionReturn(0);
}

MoFEMErrorCode DMCreateInterpolation_MGViaApproxOrders(DM dm1, DM dm2, Mat *mat,
                                                       Vec *vec) {
  PetscValidHeaderSpecific(dm1, DM_CLASSID, 1);
  PetscValidHeaderSpecific(dm2, DM_CLASSID, 1);
  MoFEMFunctionBegin;

  MPI_Comm comm;
  CHKERR PetscObjectGetComm((PetscObject)dm1, &comm);

  int m, n, M, N;

  DM dm_down = dm1;
  DM dm_up = dm2;

  int dm_down_leveldown = dm_down->leveldown;
  int dm_up_leveldown = dm_up->leveldown;

  PetscInfo2(dm1,
             "Create interpolation DMMGViaApproxOrders dm1_leveldown = %d "
             "dm2_leveldown = %d\n",
             dm_down_leveldown, dm_up_leveldown);

  IS is_down, is_up;
  {
    // Coarser mesh
    GET_DM_FIELD(dm_down);
    if (static_cast<int>(dm_field->coarseningIS.size()) < dm_down_leveldown) {
      SETERRQ(PETSC_COMM_WORLD, MOFEM_DATA_INCONSISTENCY, "data inconsistency");
    }
    is_down = dm_field->coarseningIS[dm_field->coarseningIS.size() - 1 -
                                     dm_down_leveldown];
    CHKERR ISGetSize(is_down, &M);
    CHKERR ISGetLocalSize(is_down, &m);
  }
  {
    // Finer mesh
    GET_DM_FIELD(dm_up);
    if (static_cast<int>(dm_field->coarseningIS.size()) < dm_up_leveldown) {
      SETERRQ(PETSC_COMM_WORLD, MOFEM_DATA_INCONSISTENCY, "data inconsistency");
    }
    is_up =
        dm_field
            ->coarseningIS[dm_field->coarseningIS.size() - 1 - dm_up_leveldown];
    CHKERR ISGetSize(is_up, &N);
    CHKERR ISGetLocalSize(is_up, &n);
  }

  // is_dow rows
  // is_up columns

  CHKERR MatCreate(comm, mat);
  CHKERR MatSetSizes(*mat, m, n, M, N);
  CHKERR MatSetType(*mat, MATMPIAIJ);
  CHKERR MatMPIAIJSetPreallocation(*mat, 1, PETSC_NULL, 0, PETSC_NULL);

  // get matrix layout
  PetscLayout rmap, cmap;
  CHKERR MatGetLayouts(*mat, &rmap, &cmap);
  int rstart, rend, cstart, cend;
  CHKERR PetscLayoutGetRange(rmap, &rstart, &rend);
  CHKERR PetscLayoutGetRange(cmap, &cstart, &cend);

  const int *row_indices_ptr, *col_indices_ptr;
  CHKERR ISGetIndices(is_down, &row_indices_ptr);
  CHKERR ISGetIndices(is_up, &col_indices_ptr);

  map<int, int> idx_map;
  for (int ii = 0; ii < m; ii++) {
    idx_map[row_indices_ptr[ii]] = rstart + ii;
  }

  CHKERR MatZeroEntries(*mat);
  // FIXME: Use MatCreateMPIAIJWithArrays and set array directly
  for (int jj = 0; jj < n; jj++) {
    map<int, int>::iterator mit = idx_map.find(col_indices_ptr[jj]);
    if (mit != idx_map.end()) {
      CHKERR MatSetValue(*mat, mit->second, cstart + jj, 1, INSERT_VALUES);
    }
  }

  CHKERR ISRestoreIndices(is_down, &row_indices_ptr);
  CHKERR ISRestoreIndices(is_up, &col_indices_ptr);

  CHKERR MatAssemblyBegin(*mat, MAT_FINAL_ASSEMBLY);
  CHKERR MatAssemblyEnd(*mat, MAT_FINAL_ASSEMBLY);

  if (vec != NULL) {
    *vec = PETSC_NULL;
  }

  MoFEMFunctionReturn(0);
}

MoFEMErrorCode DMCreateGlobalVector_MGViaApproxOrders(DM dm, Vec *g) {
  PetscValidHeaderSpecific(dm, DM_CLASSID, 1);
  MoFEMFunctionBegin;
  int leveldown = dm->leveldown;
  GET_DM_FIELD(dm);
  if (dm_field->kspOperators.empty()) {
    CHKERR DMCreateGlobalVector_MoFEM(dm, g);
  } else {
#if PETSC_VERSION_GE(3, 5, 3)
    CHKERR MatCreateVecs(
        dm_field->kspOperators[dm_field->kspOperators.size() - 1 - leveldown],
        g, NULL);
#else
    CHKERR MatGetVecs(
        dm_field->kspOperators[dm_field->kspOperators.size() - 1 - leveldown],
        g, NULL);
#endif
    CHKERR VecSetDM(*g, dm);
  }
  PetscInfo1(dm, "Create global vector DMMGViaApproxOrders leveldown = %d\n",
             dm->leveldown);
  MoFEMFunctionReturn(0);
}

MoFEMErrorCode PCMGSetUpViaApproxOrdersCtx::getOptions() {
  MoFEMFunctionBeginHot;
  ierr = PetscOptionsBegin(PETSC_COMM_WORLD, "",
                           "MOFEM Multi-Grid (Orders) pre-conditioner", "none");
  CHKERRG(ierr);

  CHKERR PetscOptionsInt("-mofem_mg_levels", "nb levels of multi-grid solver",
                         "", 2, &nbLevels, PETSC_NULL);
  CHKERR PetscOptionsInt("-mofem_mg_coarse_order",
                         "approximation order of coarse level", "", 2,
                         &coarseOrder, PETSC_NULL);
  CHKERR PetscOptionsInt("-mofem_mg_order_at_last_level", "order at last level",
                         "", 100, &orderAtLastLevel, PETSC_NULL);
  CHKERR PetscOptionsInt("-mofem_mg_verbose", "nb levels of multi-grid solver",
                         "", 0, &verboseLevel, PETSC_NULL);
  PetscBool shell_sub_a = shellSubA ? PETSC_TRUE : PETSC_FALSE;
  CHKERR PetscOptionsBool("-mofem_mg_shell_a", "use shell matrix as sub matrix",
                          "", shell_sub_a, &shell_sub_a, NULL);
  shellSubA = (shellSubA == PETSC_TRUE);

  ierr = PetscOptionsEnd();
  CHKERRG(ierr);
  MoFEMFunctionReturnHot(0);
}

MoFEMErrorCode PCMGSetUpViaApproxOrdersCtx::createIsAtLevel(int kk, IS *is) {
  MoFEM::Interface *m_field_ptr;
  MoFEM::ISManager *is_manager_ptr;
  MoFEMFunctionBegin;
  // if is last level, take all remaining orders dofs, if any left
  CHKERR DMoFEMGetInterfacePtr(dM, &m_field_ptr);
  CHKERR m_field_ptr->getInterface(is_manager_ptr);
  const Problem *problem_ptr;
  CHKERR DMMoFEMGetProblemPtr(dM, &problem_ptr);
  int order_at_next_level = kk + coarseOrder;
  if (kk == nbLevels - 1) {
    int first = problem_ptr->getNumeredRowDofsPtr()
                    ->get<PetscGlobalIdx_mi_tag>()
                    .lower_bound(0)
                    ->get()
                    ->getPetscGlobalDofIdx();
    CHKERR ISCreateStride(PETSC_COMM_WORLD, problem_ptr->getNbLocalDofsRow(),
                          first, 1, is);
    MoFEMFunctionReturnHot(0);
    // order_at_next_level = orderAtLastLevel;
  }
  string problem_name = problem_ptr->getName();
  CHKERR is_manager_ptr->isCreateProblemOrder(problem_name, ROW, 0,
                                              order_at_next_level, is);
  MoFEMFunctionReturn(0);
}

MoFEMErrorCode PCMGSetUpViaApproxOrdersCtx::destroyIsAtLevel(int kk, IS *is) {
  MoFEMFunctionBegin;
  CHKERR ISDestroy(is);
  MoFEMFunctionReturn(0);
}

MoFEMErrorCode
PCMGSetUpViaApproxOrdersCtx::buildProlongationOperator(bool use_mat_a,
                                                       int verb) {
  MoFEMFunctionBegin;
  verb = verb > verboseLevel ? verb : verboseLevel;

  MPI_Comm comm;
  CHKERR PetscObjectGetComm((PetscObject)dM, &comm);

  MOFEM_LOG_C("PCMGViaApproximationOrders", Sev::inform, "set MG levels %u\n",
              nbLevels);

  MOFEM_LOG_CHANNEL("SYNC");
  MOFEM_LOG_TAG("SYNC", "PCMGViaApproximationOrders")

  std::vector<IS> is_vec(nbLevels + 1);
  std::vector<int> is_glob_size(nbLevels + 1), is_loc_size(nbLevels + 1);

  for (int kk = 0; kk < nbLevels; kk++) {

    // get indices up to up to give approximation order
    CHKERR createIsAtLevel(kk, &is_vec[kk]);
    CHKERR ISGetSize(is_vec[kk], &is_glob_size[kk]);
    CHKERR ISGetLocalSize(is_vec[kk], &is_loc_size[kk]);

    MOFEM_LOG_C("SYNC", Sev::inform,
                "Nb. dofs at level [ %d ] global %u local %d\n", kk,
                is_glob_size[kk], is_loc_size[kk]);

    // if no dofs on level kk finish here
    if (is_glob_size[kk] == 0) {
      SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY, "no dofs at level");
    }
  }

  for (int kk = 0; kk != nbLevels; kk++) {
    if (kk == nbLevels - 1 && use_mat_a) {
      CHKERR DMMGViaApproxOrdersPushBackCoarseningIS(dM, is_vec[kk], A, false,
                                                     false);
    } else {
      if (kk > 0) {
        // Not coarse level
        CHKERR DMMGViaApproxOrdersPushBackCoarseningIS(dM, is_vec[kk], A, true,
                                                       shellSubA);
      } else {
        // Coarse lave is compressed matrix allowing for factorization when
        // needed
        CHKERR DMMGViaApproxOrdersPushBackCoarseningIS(dM, is_vec[kk], A, true,
                                                       false);
      }
    }
  }

  for (unsigned int kk = 0; kk < is_vec.size(); kk++) {
    CHKERR destroyIsAtLevel(kk, &is_vec[kk]);
  }

  MOFEM_LOG_SEVERITY_SYNC(comm, Sev::inform);
  MOFEM_LOG_CHANNEL("SYNC");

  MoFEMFunctionReturn(0);
}

boost::shared_ptr<PCMGSetUpViaApproxOrdersCtx>
createPCMGSetUpViaApproxOrdersCtx(DM dm, Mat A, bool use_shell_mat) {
  return boost::make_shared<PCMGSetUpViaApproxOrdersCtx>(dm, A, use_shell_mat);
}

MoFEMErrorCode PCMGSetUpViaApproxOrders(
    PC pc, boost::shared_ptr<PCMGSetUpViaApproxOrdersCtx> ctx, int verb) {

  PetscValidHeaderSpecific(pc, PC_CLASSID, 1);
  MoFEMFunctionBegin;

  auto core_log = logging::core::get();
  core_log->add_sink(LogManager::createSink(LogManager::getStrmWorld(),
                                            "PCMGViaApproximationOrders"));
  LogManager::setLog("PCMGViaApproximationOrders");
  MOFEM_LOG_TAG("PCMGViaApproximationOrders", "PCMGViaApproximationOrders");

  MOFEM_LOG("PCMGViaApproximationOrders", Sev::verbose)
      << "Setup PCMGSetUpViaApproxOrders";

  CHKERR ctx->getOptions();
  CHKERR ctx->buildProlongationOperator(true, verb);

#if PETSC_VERSION_GE(3, 8, 0)
  CHKERR PCMGSetGalerkin(pc, PC_MG_GALERKIN_NONE);
#else
  CHKERR PCMGSetGalerkin(pc, PETSC_FALSE);
#endif

  CHKERR PCMGSetLevels(pc, ctx->nbLevels, NULL);

  MOFEM_LOG("PCMGViaApproximationOrders", Sev::noisy)
      << "Setup PCMGSetUpViaApproxOrders <- end";

  MoFEMFunctionReturn(0);
}

} // namespace MoFEM
