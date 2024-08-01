/** \file PCMGSetUpViaApproxOrders.hpp
 * \brief header of multi-grid solver for p- adaptivity
 */

#ifndef __PCMGSETUP_VIA_APPROX_ORDERS_HPP__
#define __PCMGSETUP_VIA_APPROX_ORDERS_HPP__

namespace MoFEM {

/**
 * \brief Set DM ordering
 *
 * AO will transform indices from
 * coarseningIS ordering to ordering used to construct fine matrix.
 *
 * @param  dm [description]
 * @param  ao [description]
 * @return    [description]
 */
MoFEMErrorCode DMMGViaApproxOrdersSetAO(DM dm, AO ao);

/**
 * \brief Push back coarsening level to MG via approximation orders
 *
 * @param  DM discrete manager
 * @param  is Push back IS used for coarsening
 * @param  A  Get sub-matrix of A using is (that sets operators for coarsening
 * levels)
 * @param  subA  Returning pointer to created sub matrix
 * @param  subA  If true create sub matrix, otherwise in subA has to be valid
 * pointer to subA
 * @return Error code
 *
 * \ingroup dm
 */
MoFEMErrorCode DMMGViaApproxOrdersPushBackCoarseningIS(DM, IS is, Mat A,
                                                       bool create_sub_matrix,
                                                       bool shell_sub_a);

/**
 * \brief Pop IS form MG via approximation orders
 * @param  DM dm
 * @return    error code
 *
 * \ingroup dm
 */
MoFEMErrorCode DMMGViaApproxOrdersPopBackCoarseningIS(DM);

/**
 * \brief Clear approximation orders
 * @param  DM dm
 * @return Error code
 *
 * \ingroup dm
 */
MoFEMErrorCode DMMGViaApproxOrdersClearCoarseningIS(DM);

/**
 * \brief Register DM for Multi-Grid via approximation orders
 * @param  sname problem/dm registered name
 * @return       error code
 * \ingroup dm
 */
MoFEMErrorCode DMRegister_MGViaApproxOrders(const char sname[]);

/**
 * \brief Create DM data structure for Multi-Grid via approximation orders
 *
 * It set data structure and operators needed
 *
 * @param  dm Discrete manager
 * @return    Error code
 */
MoFEMErrorCode DMCreate_MGViaApproxOrders(DM dm);

/**
 * @brief  Destroy DM
 * 
 * @param dm 
 * @return * PetscErrorCode 
 */
PetscErrorCode DMDestroy_MGViaApproxOrders(DM dm);

/**
 * \brief Create matrix for Multi-Grid via approximation orders
 *
 * Not used directly by user
 *
 * @param  dm  Distributed mesh data structure
 * @param  M  Matrix
 * @return    Error code
 * \ingroup dm
 */
MoFEMErrorCode DMCreateMatrix_MGViaApproxOrders(DM dm, Mat *M);

/**
 * \brief Coarsen DM
 *
 * Not used directly by user
 *
 * @param  dm   Distributed mesh data structure
 * @param  comm Communicator
 * @param  dmc  Coarse distributed mesh data structure
 * @return      Error code
 *
 * \ingroup dm
 */
MoFEMErrorCode DMCoarsen_MGViaApproxOrders(DM dm, MPI_Comm comm, DM *dmc);

/**
 * \brief Create interpolation matrix between data managers dm1 and dm2
 * @param  dm1 Distributed mesh data structure
 * @param  dm2 Distributed mesh data structure
 * @param  mat Pointer to returned interpolation matrix
 * @param  vec Pointer to scaling vector here returned NULL
 * @return     Error code
 */
MoFEMErrorCode DMCreateInterpolation_MGViaApproxOrders(DM dm1, DM dm2, Mat *mat,
                                                       Vec *vec);

/**
 * \brief Create global vector for DMGViaApproxOrders
 * @param  dm Distributed mesh data structure
 * @param  g  returned pointer to vector
 * @return    Error code
 */
MoFEMErrorCode DMCreateGlobalVector_MGViaApproxOrders(DM dm, Vec *g);

struct PCMGSetUpViaApproxOrdersCtx;

/**
 * @brief createPCMGSetUpViaApproxOrdersCtx
 * 
 * @param dm 
 * @param A 
 * @param use_shell_mat 
 * @return PCMGSetUpViaApproxOrdersCtx 
 */
boost::shared_ptr<PCMGSetUpViaApproxOrdersCtx>
createPCMGSetUpViaApproxOrdersCtx(DM dm, Mat A, bool use_shell_mat);

/**
 * \brief Function build MG structure
 * @param  pc   MG pre-conditioner
 * <http://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/PC/PCMG.html>
 * @param  ctx  MoFEM data structure for MG
 * @param  verb verbosity level
 * @return      error code
 */
MoFEMErrorCode PCMGSetUpViaApproxOrders(
    PC pc, boost::shared_ptr<PCMGSetUpViaApproxOrdersCtx> ctx, int verb = 0);

} // namespace MoFEM

#endif //__PCMGSETUP_VIA_APPROX_ORDERS_HPP__
