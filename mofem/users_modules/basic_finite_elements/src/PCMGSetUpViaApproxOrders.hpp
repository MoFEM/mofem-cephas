/** \file PCMGSetUpViaApproxOrders.hpp
 * \brief header of multi-grid solver for p- adaptivity
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

#ifndef __PCMGSETUP_VIA_APPROX_ORDERS_HPP__
#define __PCMGSETUP_VIA_APPROX_ORDERS_HPP__

static const int DMMGVIAAPPROXORDERSCTX_INTERFACE = 1<<1;
static const MOFEMuuid IDD_DMMGVIAAPPROXORDERSCTX = MOFEMuuid(BitIntefaceId(DMMGVIAAPPROXORDERSCTX_INTERFACE));

/**
 * \brief Structure for DM for multi-grid via approximation orders
 * \ingroup dm
 */
struct DMMGViaApproxOrdersCtx: public MoFEM::DMCtx {

  PetscErrorCode queryInterface(const MOFEMuuid& uuid,MoFEM::UnknownInterface** iface);

  DMMGViaApproxOrdersCtx();
  ~DMMGViaApproxOrdersCtx();

  AO aO;
  vector<IS> coarseningIS;   ///< Coarsening IS
  vector<Mat> kspOperators;  ///< Get KSP operators

};

/**
 * \brief Set DM ordering
 *
 * IS can be given is some other ordering, AO will transform indices from coarseningIS ordering
 * to ordering used to construct fine matrix.
 *
 * @param  dm [description]
 * @param  ao [description]
 * @return    [description]
 */
PetscErrorCode DMMGViaApproxOrdersSetAO(DM dm,AO ao);

/**
 * \brief Gets size of coarseningIS in internal data struture DMMGViaApproxOrdersCtx
 * @param  dm   DM
 * @param  size size of coarseningIS
 * @return      Error code
 */
PetscErrorCode DMMGViaApproxOrdersGetCoarseningISSize(DM dm,int *size);

/**
 * \brief Push back coarsening level to MG via approximation orders
 *
 * @param  DM discrete manager
 * @param  is Push back IS used for coarsening
 * @param  A  Get sub-matrix of A using is (that sets operators for coarsening levels)
 * @param  subA  Returning pointer to created sub matrix
 * @param  subA  If true create sub matrix, otherwise in subA has to be valid pointer to subA
 * @return Error code
 *
 * \ingroup dm
 */
PetscErrorCode DMMGViaApproxOrdersPushBackCoarseningIS(DM,IS is,Mat A,Mat *subA,bool create_sub_matrix);

/**
 * \brief Pop is form MG via approximation orders
 * @param  DM dm
 * @param  is pop back IS
 * @return    error code
 *
 * \ingroup dm
 */
PetscErrorCode DMMGViaApproxOrdersPopBackCoarseningIS(DM);

/**
 * \brief Replace coarsening IS in DM via approximation orders
 * @param  dm       dm
 * @param  is_vec   Pointer to vector of is
 * @param  nb_elems Number of elements
 * @param  A        Fine matrix
 * @return          Error code
 */
PetscErrorCode DMMGViaApproxOrdersReplaceCoarseningIS(DM dm,IS *is_vec,int nb_elems,Mat A,int verb = 0);

/**
 * \brief Register DM for Multi-Grid via approximation orders
 * @param  sname problem/dm registered name
 * @return       error code
 * \ingroup dm
 */
PetscErrorCode DMRegister_MGViaApproxOrders(const char sname[]);

/**
 * \brief Create DM data structure for Multi-Grid via approximation orders
 *
 * It set data structure and operators needed
 *
 * @param  dm Discrete manager
 * @return    Error code
 */
PetscErrorCode DMCreate_MGViaApproxOrders(DM dm);

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
PetscErrorCode DMCreateMatrix_MGViaApproxOrders(DM dm,Mat *M);

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
PetscErrorCode DMCoarsen_MGViaApproxOrders(DM dm, MPI_Comm comm, DM *dmc);

/**
 * \brief Create interpolation matrix between data managers dm1 and dm2
 * @param  dm1 Distributed mesh data structure
 * @param  dm2 Distributed mesh data structure
 * @param  mat Pointer to returned interpolation matrix
 * @param  vec Pointer to scaling vector here returned NULL
 * @return     Error code
 */
PetscErrorCode DMCreateInterpolation_MGViaApproxOrders(DM dm1,DM dm2,Mat *mat,Vec *vec);

/**
 * \brief Create global vector for DMGViaApproxOrders
 * @param  dm Distributed mesh data structure
 * @param  g  returned pointer to vector
 * @return    Error code
 */
PetscErrorCode DMCreateGlobalVector_MGViaApproxOrders(DM dm,Vec *g);

/**
 * \brief Set data structures of MG pre-conditioner via approximation orders
 */
struct PCMGSetUpViaApproxOrdersCtx {

  // FieldInterface *mFieldPtr;		///< MoFEM interface
  // string problemName;			      ///< Problem name

  DM dM;  ///< Distributed mesh manager
  Mat A;  ///< Matrix at fine level

  PCMGSetUpViaApproxOrdersCtx(
    DM dm,Mat a
  ):
  // mFieldPtr(mfield_ptr),
  // problemName(problem_name),
  dM(dm),
  A(a),
  nbLevels(2),
  coarseOrder(2),
  orderAtLastLevel(1000),
  verboseLevel(0) {
  }

  ~PCMGSetUpViaApproxOrdersCtx() {
  }

  int nbLevels;				///< number of multi-grid levels
  int coarseOrder;			///< approximation order of coarse level
  int orderAtLastLevel;  ///< set maximal evaluated order

  int verboseLevel;

  PetscErrorCode ierr;

  /**
   * \brief get options from line command
   * @return error code
   */
  virtual PetscErrorCode getOptions();

  /**
   * \brief Set IS for levels
   * @param  kk level
   * @param  is pointer to IS
   * @return    error code
   */
  virtual PetscErrorCode createIsAtLevel(int kk,IS *is);

  /**
   * \brief Destroy IS if internally created
   * @param  kk level
   * @param  is pointer to is
   * @return    error code
   */
  virtual PetscErrorCode destroyIsAtLevel(int kk,IS *is);

  /**
   * \brief Set up data structures for MG
   * @param  pc   MG pre-conditioner <http://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/PC/PCMG.html>
   * @param  verb verbosity level
   * @return      error code
   */
  virtual PetscErrorCode buildProlongationOperator(PC pc,int verb = 0);

};

/**
 * \brief Function build MG structure
 * @param  pc   MG pre-conditioner <http://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/PC/PCMG.html>
 * @param  ctx  MoFEM data structure for MG
 * @param  verb verbosity level
 * @return      error code
 */
PetscErrorCode PCMGSetUpViaApproxOrders(
  PC pc,PCMGSetUpViaApproxOrdersCtx *ctx,int verb = 0
);

#endif //__PCMGSETUP_VIA_APPROX_ORDERS_HPP__
