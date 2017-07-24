/* \fiele MethodForForceScaling.hpp

*/

/* This file is part of MoFEM.
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
 * License along with MoFEM. If not, see <http://www.gnu.org/licenses/>. */

#ifndef __METHOD_FOR_FORCE_SCALING_HPP__
#define __METHOD_FOR_FORCE_SCALING_HPP__

/// Class used to scale loads, f.e. in arc-length control
struct MethodForForceScaling {

  virtual PetscErrorCode scaleNf(const FEMethod *fe,VectorDouble& Nf) = 0;
  virtual PetscErrorCode getForceScale(const double ts_t,double& scale) {
    PetscFunctionBegin;
    SETERRQ(PETSC_COMM_SELF,MOFEM_NOT_IMPLEMENTED,"not implemented");
    PetscFunctionReturn(0);
  }

  static PetscErrorCode applyScale(
    const FEMethod *fe,
    boost::ptr_vector<MethodForForceScaling> &methodsOp,VectorDouble &Nf) {
      
      PetscFunctionBegin;
      boost::ptr_vector<MethodForForceScaling>::iterator vit = methodsOp.begin();
      for(;vit!=methodsOp.end();vit++) {
        ierr = vit->scaleNf(fe,Nf); CHKERRQ(ierr);
      }
      PetscFunctionReturn(0);
    }

    virtual ~MethodForForceScaling() {}

  };

#endif //__METHOD_FOR_FORCE_SCALING_HPP__
