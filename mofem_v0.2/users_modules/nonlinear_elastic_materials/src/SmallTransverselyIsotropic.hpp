/** \file Hooke.hpp
 * \ingroup nonlinear_elastic_elem
 * \brief Implementation of linear elastic material
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

#ifndef __HOOKE_HPP__
#define __HOOKE_HPP__

#ifndef WITH_ADOL_C
  #error "MoFEM need to be compiled with ADOL-C"
#endif

/** \brief Hook equation
  * \ingroup nonlinear_elastic_elem
  */
template<typename TYPE>
struct SmallStrainTranverslyIsotropic: public NonlinearElasticElement::FunctionsToCalculatePiolaKirchhoffI<TYPE> {

    SmallStrainTranverslyIsotropic(): NonlinearElasticElement::FunctionsToCalculatePiolaKirchhoffI<TYPE>() {}

    ublas::matrix<TYPE> sTrain;
    ublas::vector<TYPE> voightStrain;
    TYPE tR;

    PetscErrorCode calculateStrain()
      PetscFunctionBegin;
      sTrain.resize(3,3,false);
      noalias(sTrain) = this->F;
      for(int dd = 0;dd<3;dd++) {
        sTrain(dd,dd) -= 1;
      }
      sTrain += trans(sTrain);
      sTrain *= 0.5;
      voightStrain.resize(6,false);
      voightStrain[0] = Strain(0,0);
      voightStrain[1] = Strain(1,1);
      voightStrain[2] = Strain(2,2);
      voightStrain[3] = 2*Strain(0,1);
      voightStrain[4] = 2*Strain(1,2);
      voightStrain[5] = 2*Strain(2,0);
      PetscFunctionReturn(0);
    )

    double nu_p, nu_pz, E_p, E_z, G_zp;
    ublas::symmetric_matrix<TYPE,ublas::upper> localStiffnessMatrix;
    PetscErrorCode calculateLocalStiffnesMatrix() {
      PetscFunctionBegin;
      double nu_zp=(nu_pz*E_z)/E_p;
      double delta=((1+nu_p)*(1-nu_p-(2*nu_pz*nu_zp)))/(E_p*E_p*E_z);

      localStiffnessMatrix.resize(6);
      localStiffnessMatrix.clear();
      localStiffnessMatrix(0,0)=localStiffnessMatrix(1,1)=(1-nu_pz*nu_zp)/(E_p*E_z*delta);
      localStiffnessMatrix(2,2)=(1-nu_p*nu_p)/(E_p*E_p*delta);

      localStiffnessMatrix(0,1)=localStiffnessMatrix(1,0)=(nu_p+nu_zp*nu_pz)/(E_p*E_z*delta);
      localStiffnessMatrix(0,2)=localStiffnessMatrix(2,0)=localStiffnessMatrix(1,2)=localStiffnessMatrix(2,1)=(nu_zp+nu_p*nu_zp)/(E_p*E_z*delta);

      localStiffnessMatrix(3,3)=E_p/(2*(1+nu_p));
      localStiffnessMatrix(4,4)=localStiffnessMatrix(5,5)=G_zp;

      PetscFunctionReturn(0);
    }

    ublas::matrix<TYPE> aARotMat;
    ublas::vector<TYPE> axVector;
    TYPE axAngle;

    /**
    * \brief Function to Calculate the Rotation Matrix at a given axis and angle of rotation

    * This function computes the rotational matrix for a given axis of rotation
    * and angle of rotation about that angle <br>

    *\param axVector A vector representing the axis of rotation
    *\param axAngle Angle of rotation along the axis (in radians)
    */
    PetscErrorCode calculateAxisAngleRotationalMatrix(ublas::matrix<TYPE> &axVector,TYPE axAngle){
      PetscFunctionBegin;

      aARotMat.resize(3,3,false);
      aARotMat.clear();

      aARotMat(0,0) = 1-((1-cos(axAngle))*(pow(axVector[1],2)+pow(axVector[2],2))/pow(norm_axVector,2));
      aARotMat(1,1) = 1-((1-cos(axAngle))*(pow(axVector[0],2)+pow(axVector[2],2))/pow(norm_axVector,2));
      aARotMat(2,2) = 1-((1-cos(axAngle))*(pow(axVector[0],2)+pow(axVector[1],2))/pow(norm_axVector,2));

      aARotMat(0,1) = ((1-cos(axAngle))*axVector[0]*axVector[1]-norm_axVector*axVector[2]*sin(axAngle))/pow(norm_axVector,2);
      aARotMat(1,0) = ((1-cos(axAngle))*axVector[0]*axVector[1]+norm_axVector*axVector[2]*sin(axAngle))/pow(norm_axVector,2);

      aARotMat(0,2) = ((1-cos(axAngle))*axVector[0]*axVector[2]+norm_axVector*axVector[1]*sin(axAngle))/pow(norm_axVector,2);
      aARotMat(2,0) = ((1-cos(axAngle))*axVector[0]*axVector[2]-norm_axVector*axVector[1]*sin(axAngle))/pow(norm_axVector,2);

      aARotMat(1,2) = ((1-cos(axAngle))*axVector[1]*axVector[2]-norm_axVector*axVector[0]*sin(axAngle))/pow(norm_axVector,2);
      aARotMat(2,1) = ((1-cos(axAngle))*axVector[1]*axVector[2]+norm_axVector*axVector[0]*sin(axAngle))/pow(norm_axVector,2);

      PetscFunctionReturn(0);
    }


    ublas::matrix<TYPE> stressRotMat;

    /**
    * \brief Function to Calculate Stress Transformation Matrix
    * This function computes the stress transformation Matrix at a give axis and angle of rotation <br>
    * One can also output the axis/angle rotational Matrix
    */
    PetscFunctionBegin stressTransformation() {
      PetscFunctionBegin;

      stressRotMat.resize(6,6,false);
      stressRotMat.clear();

      stressRotMat(0, 0) =       aARotMat(0,0) * aARotMat(0,0);
      stressRotMat(0, 1) =       aARotMat(1,0) * aARotMat(1,0);
      stressRotMat(0, 2) =       aARotMat(2,0) * aARotMat(2,0);
      stressRotMat(0, 3) = 2.0 * aARotMat(1,0) * aARotMat(0,0);
      stressRotMat(0, 4) = 2.0 * aARotMat(2,0) * aARotMat(1,0);
      stressRotMat(0, 5) = 2.0 * aARotMat(0,0) * aARotMat(2,0);

      stressRotMat(1, 0) =       aARotMat(0,1) * aARotMat(0,1);
      stressRotMat(1, 1) =       aARotMat(1,1) * aARotMat(1,1);
      stressRotMat(1, 2) =       aARotMat(2,1) * aARotMat(2,1);
      stressRotMat(1, 3) = 2.0 * aARotMat(1,1) * aARotMat(0,1);
      stressRotMat(1, 4) = 2.0 * aARotMat(2,1) * aARotMat(1,1);
      stressRotMat(1, 5) = 2.0 * aARotMat(0,1) * aARotMat(2,1);

      stressRotMat(2, 0) =       aARotMat(0,2) * aARotMat(0,2);
      stressRotMat(2, 1) =       aARotMat(1,2) * aARotMat(1,2);
      stressRotMat(2, 2) =       aARotMat(2,2) * aARotMat(2,2);
      stressRotMat(2, 3) = 2.0 * aARotMat(1,2) * aARotMat(0,2);
      stressRotMat(2, 4) = 2.0 * aARotMat(2,2) * aARotMat(1,2);
      stressRotMat(2, 5) = 2.0 * aARotMat(0,2) * aARotMat(2,2);

      stressRotMat(3, 0) =   aARotMat(0,1) * aARotMat(0,0);
      stressRotMat(3, 1) =   aARotMat(1,1) * aARotMat(1,0);
      stressRotMat(3, 2) =   aARotMat(2,1) * aARotMat(2,0);
      stressRotMat(3, 3) = ( aARotMat(1,1) * aARotMat(0,0) + aARotMat(0,1) * aARotMat(1,0) );
      stressRotMat(3, 4) = ( aARotMat(2,1) * aARotMat(1,0) + aARotMat(1,1) * aARotMat(2,0) );
      stressRotMat(3, 5) = ( aARotMat(0,1) * aARotMat(2,0) + aARotMat(2,1) * aARotMat(0,0) );

      stressRotMat(4, 0) =   aARotMat(0,2) * aARotMat(0,1);
      stressRotMat(4, 1) =   aARotMat(1,2) * aARotMat(1,1);
      stressRotMat(4, 2) =   aARotMat(2,2) * aARotMat(2,1);
      stressRotMat(4, 3) = ( aARotMat(1,2) * aARotMat(0,1) + aARotMat(0,2) * aARotMat(1,1) );
      stressRotMat(4, 4) = ( aARotMat(2,2) * aARotMat(1,1) + aARotMat(1,2) * aARotMat(2,1) );
      stressRotMat(4, 5) = ( aARotMat(0,2) * aARotMat(2,1) + aARotMat(2,2) * aARotMat(0,1) );

      stressRotMat(5, 0) =   aARotMat(0,0) * aARotMat(0,2);
      stressRotMat(5, 1) =   aARotMat(1,0) * aARotMat(1,2);
      stressRotMat(5, 2) =   aARotMat(2,0) * aARotMat(2,2);
      stressRotMat(5, 3) = ( aARotMat(1,0) * aARotMat(0,2) + aARotMat(0,0) * aARotMat(1,2) );
      stressRotMat(5, 4) = ( aARotMat(2,0) * aARotMat(1,2) + aARotMat(1,0) * aARotMat(2,2) );
      stressRotMat(5, 5) = ( aARotMat(0,0) * aARotMat(2,2) + aARotMat(2,0) * aARotMat(0,2) );

      PetscFunctionReturn(0);
    }

    ublas::matrix<TYPE> strainRotMat;

    /**
    * \brief Function to Calculate Strain Transformation Matrix<br>
    * This function computes the strain transformation Matrix at a give axis and angle of rotation <br>
    * One can also output the axis/angle rotational Matrix
    */
    PetscFunctionBegin strainTransformation() {
      PetscFunctionBegin;

      strainRotMat.resize(6,6,false);
      strainRotMat.clear();

      strainRotMat(0, 0) = aARotMat(0,0) * aARotMat(0,0);
      strainRotMat(0, 1) = aARotMat(1,0) * aARotMat(1,0);
      strainRotMat(0, 2) = aARotMat(2,0) * aARotMat(2,0);
      strainRotMat(0, 3) = aARotMat(1,0) * aARotMat(0,0);
      strainRotMat(0, 4) = aARotMat(2,0) * aARotMat(1,0);
      strainRotMat(0, 5) = aARotMat(0,0) * aARotMat(2,0);

      strainRotMat(1, 0) = aARotMat(0,1) * aARotMat(0,1);
      strainRotMat(1, 1) = aARotMat(1,1) * aARotMat(1,1);
      strainRotMat(1, 2) = aARotMat(2,1) * aARotMat(2,1);
      strainRotMat(1, 3) = aARotMat(1,1) * aARotMat(0,1);
      strainRotMat(1, 4) = aARotMat(2,1) * aARotMat(1,1);
      strainRotMat(1, 5) = aARotMat(0,1) * aARotMat(2,1);

      strainRotMat(2, 0) = aARotMat(0,2) * aARotMat(0,2);
      strainRotMat(2, 1) = aARotMat(1,2) * aARotMat(1,2);
      strainRotMat(2, 2) = aARotMat(2,2) * aARotMat(2,2);
      strainRotMat(2, 3) = aARotMat(1,2) * aARotMat(0,2);
      strainRotMat(2, 4) = aARotMat(2,2) * aARotMat(1,2);
      strainRotMat(2, 5) = aARotMat(0,2) * aARotMat(2,2);

      strainRotMat(3, 0) = 2.0 * aARotMat(0,1) * aARotMat(0,0);
      strainRotMat(3, 1) = 2.0 * aARotMat(1,1) * aARotMat(1,0);
      strainRotMat(3, 2) = 2.0 * aARotMat(2,1) * aARotMat(2,0);
      strainRotMat(3, 3) =     ( aARotMat(1,1) * aARotMat(0,0) + aARotMat(0,1) * aARotMat(1,0) );
      strainRotMat(3, 4) =     ( aARotMat(2,1) * aARotMat(1,0) + aARotMat(1,1) * aARotMat(2,0) );
      strainRotMat(3, 5) =     ( aARotMat(0,1) * aARotMat(2,0) + aARotMat(2,1) * aARotMat(0,0) );

      strainRotMat(4, 0) = 2.0 * aARotMat(0,2) * aARotMat(0,1);
      strainRotMat(4, 1) = 2.0 * aARotMat(1,2) * aARotMat(1,1);
      strainRotMat(4, 2) = 2.0 * aARotMat(2,2) * aARotMat(2,1);
      strainRotMat(4, 3) =     ( aARotMat(1,2) * aARotMat(0,1) + aARotMat(0,2) * aARotMat(1,1) );
      strainRotMat(4, 4) =     ( aARotMat(2,2) * aARotMat(1,1) + aARotMat(1,2) * aARotMat(2,1) );
      strainRotMat(4, 5) =     ( aARotMat(0,2) * aARotMat(2,1) + aARotMat(2,2) * aARotMat(0,1) );

      strainRotMat(5, 0) = 2.0 * aARotMat(0,0) * aARotMat(0,2);
      strainRotMat(5, 1) = 2.0 * aARotMat(1,0) * aARotMat(1,2);
      strainRotMat(5, 2) = 2.0 * aARotMat(2,0) * aARotMat(2,2);
      strainRotMat(5, 3) =     ( aARotMat(1,0) * aARotMat(0,2) + aARotMat(0,0) * aARotMat(1,2) );
      strainRotMat(5, 4) =     ( aARotMat(2,0) * aARotMat(1,2) + aARotMat(1,0) * aARotMat(2,2) );
      strainRotMat(5, 5) =     ( aARotMat(0,0) * aARotMat(2,2) + aARotMat(2,0) * aARotMat(0,2) );

      PetscFunctionReturn(0);
    }

  };


  ublas::matrix<TYPE> globalStiffnessMatrix;

  PetscErrorCode calculateGlobalStiffnesMatrix() {
    PetscFunctionBegin;

    globalStiffnessMatrix.resize(6,6,false);
    noalias(globalStiffnessMatrix) = prod(stressRotMat,prod(localStiffnessMatrix,strainRotMat));

    PetscFunctionReturn(0);
  }


  ublas::vector<TYPE> voigtStress;

  /** \brief Calculate global stress
  *
  */
  virtual PetscErrorCode calculateP_PiolaKirchhoffI(
    const NonlinearElasticElement::BlockData block_data,
    const NumeredMoFEMFiniteElement *fe_ptr
  ) {
    PetscFunctionBegin;
    ierr = calculateStrain(); CHKERRQ(ierr);
    ierr = calculateLocalStiffnesMatrix(); CHKERRQ(ierr);
    ierr = calculateAxisAngleRotationalMatrix(); CHKERRQ(ierr);
    ierr = stressTransformation(); CHKERRQ(ierr);
    ierr = strainTransformation(); CHKERRQ(ierr);
    ierr = calculateGlobalStiffnesMatrix(); CHKERRQ(ierr);

    voigtStress.resize(6,false);
    noalias(voigtStress) = prod(globalStiffnessMatrix,voightStrain);
    this->P.resize(3,3,false);
    P(0,0) = voigtStress[0];
    P(1,1) = voigtStress[1];
    P(2,2) = voigtStress[2];
    P(0,1) = voigtStress[3];
    P(1,2) = voigtStress[4];
    P(0,2) = voigtStress[5];
    P(1,0) = P(0,1);
    P(2,1) = P(1,2);
    P(2,0) = P(0,2);

    PetscFunctionReturn(0);
  }

  /** \brief calculate density of strain energy
  *
  */
  virtual PetscErrorCode calculateElasticEnergy(
    const NonlinearElasticElement::BlockData block_data,
    const NumeredMoFEMFiniteElement *fe_ptr
  ) {
    PetscFunctionBegin;
    ierr = calculateStrain(); CHKERRQ(ierr);
    ierr = calculateLocalStiffnesMatrix(); CHKERRQ(ierr);
    ierr = calculateAxisAngleRotationalMatrix(); CHKERRQ(ierr);
    ierr = stressTransformation(); CHKERRQ(ierr);
    ierr = strainTransformation(); CHKERRQ(ierr);
    ierr = calculateGlobalStiffnesMatrix(); CHKERRQ(ierr);

    voigtStress.resize(6,false);
    noalias(voigtStress) = prod(globalStiffnessMatrix,voightStrain);
    this->eNergy += 0.5*(voigtStress,voightStrain);
    PetscFunctionReturn(0);
  }

};

#endif //__HOOKE_HPP__
