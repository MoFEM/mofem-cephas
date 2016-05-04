/** \file NeoHookean.hpp
 * \ingroup nonlinear_elastic_elem
 * \brief Implementation of Neo-Hookean elastic material
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

#ifndef __NEOHOOKEAN_HPP__
#define __NEOHOOKEAN_HPP__

/** \brief NeoHookan equation
  * \ingroup nonlinear_elastic_elem
  */
template<typename TYPE>
struct NeoHookean: public NonlinearElasticElement::FunctionsToCalculatePiolaKirchhoffI<TYPE> {

    NeoHookean(): NonlinearElasticElement::FunctionsToCalculatePiolaKirchhoffI<TYPE>() {}

    TYPE detC;
    ublas::matrix<TYPE> invC;

    /** \brief calculate second Piola Kirchoff
      *
      * \f$\mathbf{S} = \mu(\mathbf{I}-\mathbf{C}^{-1})+\lambda(\ln{J})\mathbf{C}^{-1}\f$

      For details look to: <br>
      NONLINEAR CONTINUUM MECHANICS FOR FINITE ELEMENT ANALYSIS, Javier Bonet,
      Richard D. Wood

      */
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

    virtual PetscErrorCode calculateP_PiolaKirchhoffI(
      const NonlinearElasticElement::BlockData block_data,
      const NumeredEntFiniteElement *fe_ptr) {
      PetscFunctionBegin;
      PetscErrorCode ierr;
      this->lambda = LAMBDA(block_data.E,block_data.PoissonRatio);
      this->mu = MU(block_data.E,block_data.PoissonRatio);
      ierr = this->calculateC_CauchyDefromationTensor(); CHKERRQ(ierr);
      ierr = this->NeoHooke_PiolaKirchhoffII(); CHKERRQ(ierr);
      this->P.resize(3,3);
      noalias(this->P) = prod(this->F,this->S);
      //std::cerr << "P: " << P << std::endl;
      PetscFunctionReturn(0);
    }

    TYPE logJ;

   /** \brief calculate elastic energy density
    *

    For details look to: <br>
    NONLINEAR CONTINUUM MECHANICS FOR FINITE ELEMENT ANALYSIS, Javier Bonet,
    Richard D. Wood

*/
    PetscErrorCode NeoHookean_ElasticEnergy(){
        PetscFunctionBegin;
        this->eNergy = 0;
        for(int ii = 0;ii<3;ii++) {
            this->eNergy += this->C(ii,ii);
        }
        this->eNergy = 0.5*this->mu*(this->eNergy-3);
        logJ = log(this->J);
        this->eNergy += -this->mu*logJ + 0.5*this->lambda*pow(logJ,2);
        PetscFunctionReturn(0);
    }

    virtual PetscErrorCode calculateElasticEnergy(const NonlinearElasticElement::BlockData block_data,
      const NumeredEntFiniteElement *fe_ptr) {
      PetscFunctionBegin;
      PetscErrorCode ierr;
      this->lambda = LAMBDA(block_data.E,block_data.PoissonRatio);
      this->mu = MU(block_data.E,block_data.PoissonRatio);
      ierr = this->calculateC_CauchyDefromationTensor(); CHKERRQ(ierr);
      ierr = this->dEterminatnt(this->F,this->J); CHKERRQ(ierr);
      ierr = this->NeoHookean_ElasticEnergy(); CHKERRQ(ierr);
      PetscFunctionReturn(0);
    }

};



#endif
