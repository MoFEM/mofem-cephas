/** \file VolumeLengthQuality.hpp
 * \ingroup nonlinear_elastic_elem
 * \brief Implementation of Volume-Lebgth-Quality measure with barrier
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

#ifndef __VOLUME_LENGTH_QUALITY_HPP__
#define __VOLUME_LENGTH_QUALITY_HPP__

#ifndef WITH_ADOL_C
  #error "MoFEM need to be compiled with ADOL-C"
#endif

/** \brief Volume Length Quality
  * \ingroup nonlinear_elastic_elem
  */
template<typename TYPE>
struct VolumeLengthQuality: public NonlinearElasticElement::FunctionsToCalulatePiolaKirchhoffI<TYPE> {

    enum Type {
      QUALITY,
      BARRIER_AND_QUALITY,
      BARRIER_AND_CHANGE_QUALITY,
      BARRIER_AND_CHANGE_QUALITY_SCALED_BY_VOLUME
    };

    Type tYpe;

    VolumeLengthQuality():
      NonlinearElasticElement::FunctionsToCalulatePiolaKirchhoffI<TYPE>(),
      tYpe(BARRIER_AND_QUALITY) {}

    double gAmma;
    ublas::vector<double> coordsEdges;

    ublas::vector<TYPE> deltaChi;
    ublas::vector<TYPE> deltaX;
    ublas::matrix<TYPE> Q,dXdChiT;
    TYPE lrmsSquared,q,b,detF;

    /** Get coordinates of edges using cannonical element numeration
     */
    PetscErrorCode getEdgesFromElemCoords() {
      PetscFunctionBegin;
      if(coordsEdges.empty()) {
        coordsEdges.resize(6*2*3,false);
      }
      cblas_dcopy(3,&this->opPtr->getCoords()[0*3],1,&coordsEdges[0* 3*2+0],1);
      cblas_dcopy(3,&this->opPtr->getCoords()[1*3],1,&coordsEdges[0* 3*2+3],1);
      cblas_dcopy(3,&this->opPtr->getCoords()[0*3],1,&coordsEdges[1* 3*2+0],1);
      cblas_dcopy(3,&this->opPtr->getCoords()[2*3],1,&coordsEdges[1* 3*2+3],1);
      cblas_dcopy(3,&this->opPtr->getCoords()[0*3],1,&coordsEdges[2* 3*2+0],1);
      cblas_dcopy(3,&this->opPtr->getCoords()[3*3],1,&coordsEdges[2* 3*2+3],1);
      cblas_dcopy(3,&this->opPtr->getCoords()[1*3],1,&coordsEdges[3* 3*2+0],1);
      cblas_dcopy(3,&this->opPtr->getCoords()[2*3],1,&coordsEdges[3* 3*2+3],1);
      cblas_dcopy(3,&this->opPtr->getCoords()[1*3],1,&coordsEdges[4* 3*2+0],1);
      cblas_dcopy(3,&this->opPtr->getCoords()[3*3],1,&coordsEdges[4* 3*2+3],1);
      cblas_dcopy(3,&this->opPtr->getCoords()[2*3],1,&coordsEdges[5* 3*2+0],1);
      cblas_dcopy(3,&this->opPtr->getCoords()[3*3],1,&coordsEdges[5* 3*2+3],1);
      PetscFunctionReturn(0);
    }


    /** \brief Calculate mean element edge length

      \f[
      \Delta \boldsymbol\chi = \boldsymbol\chi^1 - \boldsymbol\chi^2
      \f]

      \f[
      \Delta X = \mathbf{F} \Delta \boldsymbol\chi
      \f]

      \f[
      l_\textrm{rms} = \sqrt{\frac{1}{6} \sum_{i=0}^6 l_i^2 } = L_\textrm{rms}dl_\textrm{rms}
      \f]

     */
     PetscErrorCode calculateLrms() {
       PetscFunctionBegin;
       if(deltaChi.empty()) {
         deltaChi.resize(3);
         deltaX.resize(3);
         dXdChiT.resize(3,3);
       }
       lrmsSquared = 0;
       dXdChiT.clear();
       for(int ee = 0;ee<6;ee++) {
         for(int dd = 0;dd<3;dd++) {
           deltaChi[dd] = coordsEdges[6*ee+dd] - coordsEdges[6*dd+3+dd];
         }
         noalias(deltaX) = prod(this->F,deltaChi);
         noalias(dXdChiT) += outer_prod(deltaX,deltaChi);
         for(int dd = 0;dd<3;dd++) {
           lrmsSquared += deltaX[3*ee+dd]*deltaX[3*ee+dd];
         }
       }
       lrmsSquared = (1./6.)*lrmsSquared;
       PetscFunctionReturn(0);
     }

     /** \brief Calculate Q

     \f[
     \tilde{\mathbf{Q}} =
      \mathbf{F}^{-\mathsf{T}}
      -
      \frac{1}{2}
      \frac{1}{l^2_\textrm{rms}}
      \sum_i^6
        \Delta\mathbf{X}_i
        \Delta\boldsymbol\chi_i^\mathsf{T}
     \f]

     */
     PetscErrorCode calculateQ() {
       PetscFunctionBegin;
       if(Q.empty()) {
         Q.resize(3,3,false);
       }
       noalias(Q) = trans(this->invH) + (0.5/lrmsSquared)*dXdChiT;
       PetscFunctionReturn(0);
     }

    /** \brief Volume Length Quality

      Based on:
      Three‐dimensional brittle fracture: configurational‐force‐driven crack propagation
      International Journal for Numerical Methods in Engineering 97 (7), 531-550

      \f[
      \mathcal{B}(a)=\frac{a}{2(1-\gamma)}-\ln{(a-\gamma)}
      \f]

      \f[
      q = q_0 b,
      \quad q_0 = 6\sqrt{2}\frac{V_0}{L^3_\textrm{rms}},
      \quad b = \frac{\textrm{det}(\mathbf{F})}{\textrm{d}l^3_\textrm{rms}}
      \f]

      \f[
      \mathbf{P} = \mathcal{B}(a)\mathbf{Q},
      \f]
      where \f$a\f$ depending on problem could be \f$q\f$ or \f$b\f$.

      */
    virtual PetscErrorCode calculateP_PiolaKirchhoffI(
      const NonlinearElasticElement::BlockData block_data,
      const NumeredMoFEMFiniteElement *fe_ptr
    ) {
      PetscFunctionBegin;

      PetscErrorCode ierr;
      ierr = dEterminatnt(this->F,detF); CHKERRQ(ierr);
      ierr = getEdgesFromElemCoords(); CHKERRQ(ierr);
      ierr = calculateLrms(); CHKERRQ(ierr);
      ierr = calculateQ(); CHKERRQ(ierr);

      b = this->detF/(lrmsSquared*sqrt(lrmsSquared));
      q = 6.*sqrt(2.)*this->opPtr->getVolume()*b;

      switch(tYpe) {
        case QUALITY:
        noalias(this->P) = q*Q;
        break;
        case BARRIER_AND_QUALITY:
        noalias(this->P) = q/(2*(1-gAmma))-log(q-gAmma)*Q;
        break;
        case BARRIER_AND_CHANGE_QUALITY:
        noalias(this->P) = b/(2*(1-gAmma))-log(b-gAmma)*Q;
        break;
        case BARRIER_AND_CHANGE_QUALITY_SCALED_BY_VOLUME:
        noalias(this->P) = this->opPtr->getVolume()*this->detH*b/(2*(1-gAmma))-log(b-gAmma)*Q;
        break;
      }

      PetscFunctionReturn(0);
    }


};

#endif //__VOLUME_LENGTH_QUALITY_HPP__
