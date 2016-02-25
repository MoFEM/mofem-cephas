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

struct PCMGSetUpViaApproxOrdersCtx {

  FieldInterface *mFieldPtr;		///< MoFEM interface
  string problemName;			///< Problem name

  PCMGSetUpViaApproxOrdersCtx(FieldInterface *mfield_ptr,string problem_name):
  mFieldPtr(mfield_ptr),
  problemName(problem_name),
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

  virtual PetscErrorCode getOptions();
  virtual PetscErrorCode createIsAtLevel(int kk,IS *is);
  virtual PetscErrorCode destroyIsAtLevel(int kk,IS *is);
  virtual PetscErrorCode buildProlongationOperator(PC pc,int verb = 0);

};

PetscErrorCode PCMGSetUpViaApproxOrders(
  PC pc,PCMGSetUpViaApproxOrdersCtx *ctx,int verb = 0
);

#endif //__PCMGSETUP_VIA_APPROX_ORDERS_HPP__
