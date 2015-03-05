/** \file Hooke.hpp 
 * \brief Implementation of linear elastic material 
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

#ifndef __HOOKE_HPP__
#define __HOOKE_HPP__

#ifndef WITH_ADOL_C
  #error "MoFEM need to be compiled with ADOL-C"
#endif 

template<typename TYPE>
struct Hooke: public NonlinearElasticElement::FunctionsToCalulatePiolaKirchhoffI<TYPE> {

    Hooke(): NonlinearElasticElement::FunctionsToCalulatePiolaKirchhoffI<TYPE>() {}

    ublas::matrix<TYPE> Eps;
    TYPE tr;
    
    virtual PetscErrorCode CalualteP_PiolaKirchhoffI(
      const NonlinearElasticElement::BlockData block_data,
      const NumeredMoFEMFiniteElement *fe_ptr) {
      PetscFunctionBegin;
      //PetscErrorCode ierr;
      this->lambda = LAMBDA(block_data.E,block_data.PoissonRatio);
      this->mu = MU(block_data.E,block_data.PoissonRatio);
      Eps.resize(3,3);
      noalias(Eps) = 0.5*(this->F + trans(this->F));
      this->P.resize(3,3);
      noalias(this->P) = 2*this->mu*Eps;
      tr = 0;
      for(int dd = 0;dd<3;dd++) {
	tr += this->lambda*Eps(dd,dd);
      }
      for(int dd =0;dd<3;dd++) {
	this->P(dd,dd) += tr;
      } 
      PetscFunctionReturn(0);
    }

};

#endif //__HOOKE_HPP__
