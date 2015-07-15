/** \file Gels.hpp
  \brief Implementation of Gel finite element
  \ingroup gel
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


#ifndef __UESRGELMODEL_HPP__
#define __UESRGELMODEL_HPP__

#ifndef WITH_ADOL_C
  #error "MoFEM need to be compiled with ADOL-C"
#endif

template<typename TYPE>
struct UserGelConstitutiveEquation: public Gel::ConstitutiveEquation<TYPE>  {

  UserGelConstitutiveEquation(Gel::BlockMaterialData &data):
  Gel::ConstitutiveEquation<TYPE>(data) {
  }

  /*PetscErrorCode calculateStressAlpha() {
    PetscFunctionBegin;
    PetscFunctionReturn(0);
  }*/

  /*PetscErrorCode calculateStressBeta() {
    PetscFunctionBegin;
    PetscFunctionReturn(0);
  }*/

  /*PetscErrorCode calculateStrainHatFlux() {
    PetscFunctionBegin;
    PetscFunctionReturn(0);
  }*/

  /*PetscErrorCode calculateStressBetaHat() {
    PetscFunctionBegin;
    PetscFunctionReturn(0);
  }*/

  /*PetscErrorCode calculateSolventFlux() {
    PetscFunctionBegin;
    PetscFunctionReturn(0);
  }*/

  /*PetscErrorCode calculateVolumeDot() {
    PetscFunctionBegin;
    volumeDot = traceStrainTotalDot;
    PetscFunctionReturn(0);
  }*/

};


#endif // __UESRGELMODEL_HPP__
