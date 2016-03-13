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
 */
struct DMMGViaApproxOrdersCtx: public MoFEM::DMCtx {

  PetscErrorCode queryInterface(const MOFEMuuid& uuid,MoFEM::UnknownInterface** iface);

  DMMGViaApproxOrdersCtx();
  ~DMMGViaApproxOrdersCtx();

  Mat fineMatrix;            ///< Assembled matrix at fine level
  vector<IS> coarseningIS;   ///< Coarsening IS
  vector<Mat> kspOperators;  ///< Get KSP operators

};


/**
 * \brief Push back coarsening level
 *
 * @param  DM discrete manager
 * @param  is Push back IS used for coarsening
 * @param  A  Get sub-matrix of A using is (that sets operators for coarsening levels)
 * @return error code
 *
 * \ingroup dm
 */
PetscErrorCode DMMGViaApproxOrdersPushBackCoarseningIS(DM,IS is,Mat A,Mat *subA);

/**
 * [DMMoFEMPushBackCoarseningIS description]
 * @param  DM dm
 * @param  is pop back IS
 * @return    error code
 */
PetscErrorCode DMMGViaApproxOrdersPopBackCoarseningIS(DM);


PetscErrorCode DMRegister_MGViaApproxOrders(const char sname[]);

PetscErrorCode DMCreate_MGViaApproxOrders(DM dm);

PetscErrorCode DMCreateMatrix_MGViaApproxOrders(DM dm,Mat *M);

PetscErrorCode DMCoarsen_MGViaApproxOrders(DM dm, MPI_Comm comm, DM *dmc);

PetscErrorCode DMCreateInterpolation_MGViaApproxOrders(DM dm1,DM dm2,Mat *mat,Vec *vec);

PetscErrorCode DMRestrict_MGViaApproxOrders(DM fine,Mat mat,Vec vec,Mat mat2,DM coarse,void *ctx);


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
