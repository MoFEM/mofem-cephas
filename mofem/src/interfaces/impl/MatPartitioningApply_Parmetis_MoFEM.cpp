#ifdef PARMETIS

#include <definitions.h>
#include <Includes.hpp>

#if PETSC_VERSION_GE(3, 6, 0)
#include <petsc/private/matimpl.h>
#else
#include <petsc-private/matimpl.h>
#endif

#include <parmetis.h>

/*
      The first 5 elements of this structure are the input control array to
   Metis
*/
typedef struct {
  PetscInt cuts; /* number of cuts made (output) */
  PetscInt foldfactor;
  PetscInt parallel; /* use parallel partitioner for coarse problem */
  PetscInt indexing; /* 0 indicates C indexing, 1 Fortran */
  PetscInt printout; /* indicates if one wishes Metis to print info */
  PetscBool repartition;
} MatPartitioning_Parmetis;

/*
  MATMPIAdj format - Compressed row storage for storing adjacency lists, and
  possibly weights This is for grid reorderings (to reduce bandwidth) grid
  partitionings, etc. This is NOT currently a dynamic data-structure.

*/
typedef struct {
  PetscInt nz;
  PetscInt *diag;      /* pointers to diagonal elements, if they exist */
  PetscInt *i;         /* pointer to beginning of each row */
  PetscInt *j;         /* column values: j + i[k] is start of row k */
  PetscInt *values;    /* numerical values */
  PetscBool symmetric; /* user indicates the nonzero structure is symmetric */
  PetscBool freeaij;   /* free a, i,j at destroy */
  PetscBool freeaijwithfree; /* use free() to free i,j instead of PetscFree() */
  PetscScalar *rowvalues;    /* scalar work space for MatGetRow() */
  PetscInt rowvalues_alloc;
} Mat_MPIAdj;

#define CHKERRQPARMETIS(n, func)                                               \
  if (n == METIS_ERROR_INPUT)                                                  \
    SETERRQ(PETSC_COMM_SELF, PETSC_ERR_LIB,                                   \
             "ParMETIS error due to wrong inputs and/or options for %s",       \
             func);                                                            \
  else if (n == METIS_ERROR_MEMORY)                                            \
    SETERRQ(PETSC_COMM_SELF, PETSC_ERR_LIB,                                   \
             "ParMETIS error due to insufficient memory in %s", func);         \
  else if (n == METIS_ERROR)                                                   \
    SETERRQ(PETSC_COMM_SELF, PETSC_ERR_LIB, "ParMETIS general error in %s",   \
             func);

#define PetscStackCallParmetis(func, args)                                     \
  do {                                                                         \
    PetscStackPush(#func);                                                     \
    status = func args;                                                        \
    PetscStackPop;                                                             \
    CHKERRQPARMETIS(status, #func);                                            \
  } while (0)

namespace MoFEM {

#undef __FUNCT__
#define __FUNCT__ "MatPartitioningApply_Parmetis_MoFEM"
PetscErrorCode MatPartitioningApply_Parmetis_MoFEM(MatPartitioning part,
                                                   IS *partitioning) {
  MatPartitioning_Parmetis *pmetis = (MatPartitioning_Parmetis *)part->data;
  PetscErrorCode ierr;
  PetscInt *locals = NULL;
  Mat mat = part->adj, amat, pmat;
  PetscBool flg;
  PetscInt bs = 1;

  MoFEMFunctionBeginHot;
  ierr = PetscObjectTypeCompare((PetscObject)mat, MATMPIADJ, &flg);
  CHKERRQ(ierr);
  if (flg) {
    amat = mat;
    PetscObjectReference((PetscObject)amat);
    CHKERRQ(ierr);
  } else {
    /* bs indicates if the converted matrix is "reduced" from the original and
       hence the resulting partition results need to be stretched to match the
       original matrix */
    ierr = MatConvert(mat, MATMPIADJ, MAT_INITIAL_MATRIX, &amat);
    CHKERRQ(ierr);
    if (amat->rmap->n > 0)
      bs = mat->rmap->n / amat->rmap->n;
  }
  ierr = MatMPIAdjCreateNonemptySubcommMat(amat, &pmat);
  CHKERRQ(ierr);
  ierr = MPI_Barrier(PetscObjectComm((PetscObject)part));
  CHKERRQ(ierr);

  if (pmat) {
    MPI_Comm pcomm, comm;
    Mat_MPIAdj *adj = (Mat_MPIAdj *)pmat->data;
    PetscInt *vtxdist = pmat->rmap->range;
    PetscInt *xadj = adj->i;
    PetscInt *adjncy = adj->j;
    PetscInt itmp = 0, wgtflag = 0, numflag = 0, ncon = 1, nparts = part->n,
             options[24], i, j;
    real_t *tpwgts, *ubvec, itr = 0.1;
    int status;

    ierr = PetscObjectGetComm((PetscObject)pmat, &pcomm);
    CHKERRQ(ierr);
#if defined(PETSC_USE_DEBUG)
    /* check that matrix has no diagonal entries */
    {
      PetscInt rstart;
      ierr = MatGetOwnershipRange(pmat, &rstart, NULL);
      CHKERRQ(ierr);
      for (i = 0; i < pmat->rmap->n; i++) {
        for (j = xadj[i]; j < xadj[i + 1]; j++) {
          if (adjncy[j] == i + rstart)
            SETERRQ(
                PETSC_COMM_SELF, PETSC_ERR_ARG_WRONG,
                "Row %d has diagonal entry; Parmetis forbids diagonal entry",
                i + rstart);
        }
      }
    }
#endif

    // use vertex_weights
    if (part->vertex_weights) {
      wgtflag = 2;
    }

    ierr = PetscMalloc1(amat->rmap->n, &locals);
    CHKERRQ(ierr);

    if (PetscLogPrintInfo) {
      itmp = pmetis->printout;
      pmetis->printout = 127;
    }
    ierr = PetscMalloc1(ncon * nparts, &tpwgts);
    CHKERRQ(ierr);
    for (i = 0; i < ncon; i++) {
      for (j = 0; j < nparts; j++) {
        if (part->part_weights) {
          tpwgts[i * nparts + j] = part->part_weights[i * nparts + j];
        } else {
          tpwgts[i * nparts + j] = 1. / nparts;
        }
      }
    }
    ierr = PetscMalloc1(ncon, &ubvec);
    CHKERRQ(ierr);
    for (i = 0; i < ncon; i++) {
      ubvec[i] = 1.05;
    }
    /* This sets the defaults */
    options[0] = 0;
    for (i = 1; i < 24; i++) {
      options[i] = -1;
    }
    /* Duplicate the communicator to be sure that ParMETIS attribute caching
     * does not interfere with PETSc. */
    ierr = MPI_Comm_dup(pcomm, &comm);
    CHKERRQ(ierr);
    if (pmetis->repartition) {
      PetscStackCallParmetis(
          ParMETIS_V3_AdaptiveRepart,
          ((idx_t *)vtxdist, (idx_t *)xadj, (idx_t *)adjncy,
           (idx_t *)part->vertex_weights, (idx_t *)part->vertex_weights,
           (idx_t *)adj->values, (idx_t *)&wgtflag, (idx_t *)&numflag,
           (idx_t *)&ncon, (idx_t *)&nparts, tpwgts, ubvec, &itr,
           (idx_t *)options, (idx_t *)&pmetis->cuts, (idx_t *)locals, &comm));
    } else {
      PetscStackCallParmetis(ParMETIS_V3_PartKway,
                             ((idx_t *)vtxdist, (idx_t *)xadj, (idx_t *)adjncy,
                              (idx_t *)part->vertex_weights,
                              (idx_t *)adj->values, (idx_t *)&wgtflag,
                              (idx_t *)&numflag, (idx_t *)&ncon,
                              (idx_t *)&nparts, tpwgts, ubvec, (idx_t *)options,
                              (idx_t *)&pmetis->cuts, (idx_t *)locals, &comm));
    }
    ierr = MPI_Comm_free(&comm);
    CHKERRQ(ierr);

    ierr = PetscFree(tpwgts);
    CHKERRQ(ierr);
    ierr = PetscFree(ubvec);
    CHKERRQ(ierr);
    if (PetscLogPrintInfo)
      pmetis->printout = itmp;
  }

  if (bs > 1) {
    PetscInt i, j, *newlocals;
    ierr = PetscMalloc1(bs * amat->rmap->n, &newlocals);
    CHKERRQ(ierr);
    for (i = 0; i < amat->rmap->n; i++) {
      for (j = 0; j < bs; j++) {
        newlocals[bs * i + j] = locals[i];
      }
    }
    ierr = PetscFree(locals);
    CHKERRQ(ierr);
    ierr =
        ISCreateGeneral(PetscObjectComm((PetscObject)part), bs * amat->rmap->n,
                        newlocals, PETSC_OWN_POINTER, partitioning);
    CHKERRQ(ierr);
  } else {
    ierr = ISCreateGeneral(PetscObjectComm((PetscObject)part), amat->rmap->n,
                           locals, PETSC_OWN_POINTER, partitioning);
    CHKERRQ(ierr);
  }

  ierr = MatDestroy(&pmat);
  CHKERRQ(ierr);
  ierr = MatDestroy(&amat);
  CHKERRQ(ierr);
  MoFEMFunctionReturnHot(0);
}

} // namespace MoFEM

#endif // PARMETIS
