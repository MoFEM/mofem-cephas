/** \file MooneyRivlin.hpp
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

#ifndef __MOONEYRIVLIN_HPP__
#define __MOONEYRIVLIN_HPP__

/** \brief Mooney-Rivlin equation
  * \ingroup nonlinear_elastic_elem
  */
template<typename TYPE>
struct MooneyRivlin: public NonlinearElasticElement::FunctionsToCalculatePiolaKirchhoffI<TYPE> {

    MooneyRivlin(): NonlinearElasticElement::FunctionsToCalculatePiolaKirchhoffI<TYPE>() {}


    ublas::matrix<TYPE> SiGm;
    ublas::matrix<TYPE> invF;
    /** \brief calculate second Piola Kirchoff
      *
      Applied_mechanics_of_solids
      Allan_F_Bower
      */
    PetscErrorCode MooneyRivlin_PiolaKirchhoffII() {
      PetscFunctionBegin;
      PetscErrorCode ierr;
      this->S.resize(3,3);
      SiGm.resize(3,3);
      invF.resize(3,3);
      // double C1 = this->mu * (3 * this->lambda + 2 * this-> mu) /(this->lambda + this->mu);
      // double C2 = this->lambda / (2*(this->lambda + this->mu));
      double C1 = this->mu /4;
      double C2 = this->mu /4;
      double D = 2 / (this->lambda + 2 * this->mu / 3);
      ierr = this->dEterminatnt(this->F,this->J); CHKERRQ(ierr);
      ierr = this->iNvert(this->J,this->F,invF); CHKERRQ(ierr);
      for(int i = 0;i<3;i++) {
        for(int j = 0;j<3;j++) {
          int Kr = ((i==j) ? 1 : 0);
          for(int k = 0;k<3;k++) {
            for(int l = 0;l<3;l++) {
              SiGm(i,j) += C1 / pow(this->J,5/3) * (this->B(i,j) - (1/3)*this->B(k,k)*Kr)
              + C2 / pow(this->J,7/3) * ( this->B(k,k)*this->B(i,j) - (1/3)*this->B(k,k)*this->B(k,k)*Kr
              - this->B(i,k)*this->B(k,j) + (1/3)*this->B(k,l)*this->B(l,k)*Kr)
              + 2*D*(this->J -1)*Kr;
            }
          }
          for(int k = 0;k<3;k++) {
            for(int l = 0;l<k;l++) {
              this->S(k,l) += (this->J) * invF(k,i) * SiGm(i,j) * invF(l,j);
              this->S(l,k) = this->S(k,l);
            }
          }
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
      ierr = this->calculateB_CauchyDefromationTensor(); CHKERRQ(ierr);
      ierr = this->calculateC_CauchyDefromationTensor(); CHKERRQ(ierr);
      ierr = this->MooneyRivlin_PiolaKirchhoffII(); CHKERRQ(ierr);
      this->P.resize(3,3);
      noalias(this->P) = prod(this->F,this->S);
      //std::cerr << "P: " << P << std::endl;
      PetscFunctionReturn(0);
    }

    TYPE logJ;

   /** \brief calculate elastic energy density
    *

    For details look to: <br>
    Applied_mechanics_of_solids
    Allan_F_Bower

*/
    PetscErrorCode MooneyRivlin_ElasticEnergy(){
        PetscFunctionBegin;
        this->eNergy = 0;
        TYPE trB = 0;
        // double C1 = this->mu * (3 * this->lambda + 2 * this-> mu) /(this->lambda + this->mu);
        // double C2 = this->lambda / (2*(this->lambda + this->mu));
        double C1 = this->mu /4;
        double C2 = this->mu /4;
        double D = 2 / (this->lambda + 2 * this->mu / 3);
        for(int ii = 0;ii<3;ii++) {
            trB += this->B(ii,ii);
        }
        TYPE I2 = 0;
        for(int i = 0;i<3;i++) {
          for(int j = 0;j<3;j++) {
            I2 += 0.5*(trB*trB - this->B(i,j)*this->B(j,i));
          }
        }
        this->eNergy = C1*(trB - 3) + C2*(I2 - 3) + 1/D * (this->J -1) * (this->J -1);
        PetscFunctionReturn(0);
    }

    virtual PetscErrorCode calculateElasticEnergy(const NonlinearElasticElement::BlockData block_data,
      const NumeredEntFiniteElement *fe_ptr) {
      PetscFunctionBegin;
      PetscErrorCode ierr;
      this->lambda = LAMBDA(block_data.E,block_data.PoissonRatio);
      this->mu = MU(block_data.E,block_data.PoissonRatio);
      ierr = this->calculateB_CauchyDefromationTensor(); CHKERRQ(ierr);
      ierr = this->dEterminatnt(this->F,this->J); CHKERRQ(ierr);
      ierr = this->MooneyRivlin_ElasticEnergy(); CHKERRQ(ierr);
      PetscFunctionReturn(0);
    }

};



#endif
