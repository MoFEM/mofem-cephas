/** \file NeoHookean.hpp 
 * \brief Implementation of Neo-Hookean elastic material 
 *
 * This file is part of MoFEM.
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

#ifndef __NEOHOOKEAN_HPP__
#define __NEOHOOKEAN_HPP__

template<typename TYPE>
struct NeoHookean: public NonlinearElasticElement::FunctionsToCalulatePiolaKirchhoffI<TYPE> {

    NeoHookean(): NonlinearElasticElement::FunctionsToCalulatePiolaKirchhoffI<TYPE>() {}

    TYPE detC;
    ublas::matrix<TYPE> invC;
    
    PetscErrorCode NeoHooke_PiolaKirchhoffII() {
      PetscFunctionBegin;
      PetscErrorCode ierr;
      invC.resize(3,3);
      this->S.resize(3,3);
      ierr = this->dEterminatnt(this->C,detC); CHKERRQ(ierr);
      ierr = this->iNvert(detC,this->C,invC); CHKERRQ(ierr);
      ierr = this->dEterminatnt(this->F,this->J); CHKERRQ(ierr);
      for(int i = 0;i<3;i++) {
	for(int j = 0;j<3;j++) {
	  this->S(i,j) = this->mu*( ((i==j) ? 1 : 0) - invC(i,j) ) + this->lambda*log(this->J)*invC(i,j);
	}
      }
      PetscFunctionReturn(0);
    }

    virtual PetscErrorCode CalualteP_PiolaKirchhoffI(
      const NonlinearElasticElement::BlockData block_data,
      const NumeredMoFEMFiniteElement *fe_ptr) {
      PetscFunctionBegin;
      PetscErrorCode ierr;
      this->lambda = LAMBDA(block_data.E,block_data.PoissonRatio);
      this->mu = MU(block_data.E,block_data.PoissonRatio);
      ierr = this->CalulateC_CauchyDefromationTensor(); CHKERRQ(ierr);
      ierr = this->NeoHooke_PiolaKirchhoffII(); CHKERRQ(ierr);
      this->P.resize(3,3);
      noalias(this->P) = prod(this->F,this->S);
      //cerr << "P: " << P << endl;
      PetscFunctionReturn(0);
    }

};



#endif
