/* Copyright (C) 2013, Michel Cortis <mikecortis at gmail.com>
 * --------------------------------------------------------------
 * FIXME: DESCRIPTION
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

#ifndef __ELASTICFEMETHODTRANSISO_HPP__
#define __ELASTICFEMETHODTRANSISO_HPP__

#include <boost/numeric/ublas/symmetric.hpp>
#include "ElasticFEMethod.hpp"
#include "PostProcDisplacementAndStrainOnRefindedMesh.hpp"

#include "SnesCtx.hpp"
#include "ArcLengthTools.hpp"

namespace ObosleteUsersModules {

  /**
   * \brief Function to Calculate Transverse Isotropic Stiffness Matrix
   * this is similiar to Orthotropic Stiffness Matrix but material parameters in x and y are identical
   * hence it is used for modelling of fibre and wood
   *
   *\param E_p Young's Modulus in x and y direction
   *\param nu_p Poisson's Ratio in x-y plane
   *\param E_z Young's Modulus in z-direction (direction of fibre)
   *\param n_pz Poisson's Ratio in z-direction
   *\param G_zp Shear Modulus in z-direction
   */
  struct TransverseIsotropicStiffnessMatrix {

    double nu_p, nu_pz, E_p, E_z, G_zp;

    ublas::symmetric_matrix<FieldData,ublas::upper> StiffnessMatrix;

    TransverseIsotropicStiffnessMatrix(double nu_p, double nu_pz, double E_p, double E_z, double G_zp){

      double nu_zp=(nu_pz*E_z)/E_p;
      //double nu_zp = nu_pz;
      //nu_pz=(nu_zp*E_p)/E_z;
      double delta=((1+nu_p)*(1-nu_p-(2*nu_pz*nu_zp)))/(E_p*E_p*E_z);

      StiffnessMatrix.resize(6);
      StiffnessMatrix.clear();
      StiffnessMatrix(0,0)=StiffnessMatrix(1,1)=(1-nu_pz*nu_zp)/(E_p*E_z*delta);
      StiffnessMatrix(2,2)=(1-nu_p*nu_p)/(E_p*E_p*delta);

      StiffnessMatrix(0,1)=StiffnessMatrix(1,0)=(nu_p+nu_zp*nu_pz)/(E_p*E_z*delta);
      StiffnessMatrix(0,2)=StiffnessMatrix(2,0)=StiffnessMatrix(1,2)=StiffnessMatrix(2,1)=(nu_zp+nu_p*nu_zp)/(E_p*E_z*delta);

      StiffnessMatrix(3,3)=E_p/(2*(1+nu_p));
      StiffnessMatrix(4,4)=StiffnessMatrix(5,5)=G_zp;
    }

  };

  struct TransverseIsotropicComplianceMatrix {

    double nu_p, nu_pz, E_p, E_z, G_zp;

    ublas::symmetric_matrix<FieldData,ublas::upper> ComplianceMatrix;

    TransverseIsotropicComplianceMatrix(double nu_p, double nu_pz, double E_p, double E_z, double G_zp){

      //double nu_zp = (nu_pz*E_z)/E_p;
      double nu_zp = nu_pz;
      nu_pz=(nu_zp*E_p)/E_z;

      double Gp = E_p/(2*(1+nu_p));

      ComplianceMatrix.resize(6);
      ComplianceMatrix.clear();
      ComplianceMatrix(0,0)=ComplianceMatrix(1,1)=1/E_p;
      ComplianceMatrix(2,2)=1/E_z;

      ComplianceMatrix(0,1)=ComplianceMatrix(1,0)=-nu_p/E_p;
      ComplianceMatrix(0,2)=ComplianceMatrix(2,0)=-nu_zp/E_z;
      ComplianceMatrix(1,2)=ComplianceMatrix(2,1)=-nu_pz/E_p;

      ComplianceMatrix(3,3)=1/Gp;
      ComplianceMatrix(4,4)=ComplianceMatrix(5,5)=1/G_zp;
    }

  };

  /**
   * \brief Function to Calculate Isotropic Stiffness Matrix
   *
   *\param lambda is the Lame's first parameter, computed using LAMBDA(Young's Modulus,Poisson's Ratio)
   *\param mu is computed using MU(Young's Modulus,Poisson's Ratio)
   */
  struct IsotropicStiffnessMatrix {

    double lambda, mu;

    ublas::symmetric_matrix<FieldData,ublas::upper> StiffnessMatrix;

    IsotropicStiffnessMatrix(double lambda, double mu){

      ublas::symmetric_matrix<FieldData,ublas::upper> D_lambda,D_mu;

      D_lambda.resize(6);
      D_lambda.clear();
      for(int rr = 0;rr<3;rr++) {
        for(int cc = 0;cc<3;cc++) {
          D_lambda(rr,cc) = 1;
        }
      }
      D_mu.resize(6);
      D_mu.clear();
      for(int rr = 0;rr<6;rr++) {
        D_mu(rr,rr) = rr<3 ? 2 : 1;
      }
      StiffnessMatrix.resize(6);
      StiffnessMatrix.clear();
      StiffnessMatrix = ublas::zero_matrix<FieldData>(6,6);
      StiffnessMatrix = lambda*D_lambda + mu*D_mu;
    }

  };

  /**
   * \brief Function to Calculate the Rotation Matrix at a given axis and angle of rotation
   * This function computes the rotational matrix for a given axis of rotation and angle of rotation about that angle <br>

   *\param AxVector A vector representing the axis of rotation
   *\param AxAngle Angle of rotation along the axis (in radians)
   */
  struct AxisAngleRotationalMatrix {

    ublas::matrix<double> AARotMat;

    AxisAngleRotationalMatrix(double *AxVector, double AxAngle){

      double norm_AxVector = sqrt(pow(AxVector[0],2) + pow(AxVector[1],2) + pow(AxVector[2],2));

      AARotMat = ublas::zero_matrix<FieldData>(3,3);

      AARotMat(0,0) = 1-((1-cos(AxAngle))*(pow(AxVector[1],2)+pow(AxVector[2],2))/pow(norm_AxVector,2));
      AARotMat(1,1) = 1-((1-cos(AxAngle))*(pow(AxVector[0],2)+pow(AxVector[2],2))/pow(norm_AxVector,2));
      AARotMat(2,2) = 1-((1-cos(AxAngle))*(pow(AxVector[0],2)+pow(AxVector[1],2))/pow(norm_AxVector,2));

      AARotMat(0,1) = ((1-cos(AxAngle))*AxVector[0]*AxVector[1]-norm_AxVector*AxVector[2]*sin(AxAngle))/pow(norm_AxVector,2);
      AARotMat(1,0) = ((1-cos(AxAngle))*AxVector[0]*AxVector[1]+norm_AxVector*AxVector[2]*sin(AxAngle))/pow(norm_AxVector,2);

      AARotMat(0,2) = ((1-cos(AxAngle))*AxVector[0]*AxVector[2]+norm_AxVector*AxVector[1]*sin(AxAngle))/pow(norm_AxVector,2);
      AARotMat(2,0) = ((1-cos(AxAngle))*AxVector[0]*AxVector[2]-norm_AxVector*AxVector[1]*sin(AxAngle))/pow(norm_AxVector,2);

      AARotMat(1,2) = ((1-cos(AxAngle))*AxVector[1]*AxVector[2]-norm_AxVector*AxVector[0]*sin(AxAngle))/pow(norm_AxVector,2);
      AARotMat(2,1) = ((1-cos(AxAngle))*AxVector[1]*AxVector[2]+norm_AxVector*AxVector[0]*sin(AxAngle))/pow(norm_AxVector,2);

    }

  };

  /**
   * \brief Function to Calculate Stress Transformation Matrix
   * This function computes the stress transformation Matrix at a give axis and angle of rotation <br>
   * One can also output the axis/angle rotational Matrix
   *
   *\param AxVector A vector representing the axis of rotation
   *\param AxAngle Angle of rotation along the axis (in radians)
   */

  struct StressTransformation {

    ublas::matrix<double> StressRotMat;
    ublas::matrix<double> AARotMat;

    StressTransformation(double *AxVector, double AxAngle){

      AARotMat = ublas::zero_matrix<FieldData>(3,3);
      AxisAngleRotationalMatrix RotMat(&AxVector[0], AxAngle);
      AARotMat=RotMat.AARotMat;

      StressRotMat = ublas::zero_matrix<FieldData>(6,6);
      /*
      StressRotMat(0, 0) =       AARotMat(0,0) * AARotMat(0,0);
      StressRotMat(0, 1) =       AARotMat(0,1) * AARotMat(0,1);
      StressRotMat(0, 2) =       AARotMat(0,2) * AARotMat(0,2);
      StressRotMat(0, 3) = 2.0 * AARotMat(0,1) * AARotMat(0,2);
      StressRotMat(0, 4) = 2.0 * AARotMat(0,0) * AARotMat(0,2);
      StressRotMat(0, 5) = 2.0 * AARotMat(0,0) * AARotMat(0,1);

      StressRotMat(1, 0) =       AARotMat(1,0) * AARotMat(1,0);
      StressRotMat(1, 1) =       AARotMat(1,1) * AARotMat(1,1);
      StressRotMat(1, 2) =       AARotMat(1,2) * AARotMat(1,2);
      StressRotMat(1, 3) = 2.0 * AARotMat(1,1) * AARotMat(1,2);
      StressRotMat(1, 4) = 2.0 * AARotMat(1,0) * AARotMat(1,2);
      StressRotMat(1, 5) = 2.0 * AARotMat(1,0) * AARotMat(1,1);

      StressRotMat(2, 0) =       AARotMat(2,0) * AARotMat(2,0);
      StressRotMat(2, 1) =       AARotMat(2,1) * AARotMat(2,1);
      StressRotMat(2, 2) =       AARotMat(2,2) * AARotMat(2,2);
      StressRotMat(2, 3) = 2.0 * AARotMat(2,1) * AARotMat(2,2);
      StressRotMat(2, 4) = 2.0 * AARotMat(2,0) * AARotMat(2,2);
      StressRotMat(2, 5) = 2.0 * AARotMat(2,0) * AARotMat(2,1);

      StressRotMat(3, 0) =   AARotMat(1,0) * AARotMat(2,0);
      StressRotMat(3, 1) =   AARotMat(1,1) * AARotMat(2,1);
      StressRotMat(3, 2) =   AARotMat(1,2) * AARotMat(2,2);
      StressRotMat(3, 3) = ( AARotMat(1,1) * AARotMat(2,2) + AARotMat(1,2) * AARotMat(2,1) );
      StressRotMat(3, 4) = ( AARotMat(1,0) * AARotMat(2,2) + AARotMat(1,2) * AARotMat(2,0) );
      StressRotMat(3, 5) = ( AARotMat(1,0) * AARotMat(2,1) + AARotMat(1,1) * AARotMat(2,0) );

      StressRotMat(4, 0) =   AARotMat(0,0) * AARotMat(2,0);
      StressRotMat(4, 1) =   AARotMat(0,1) * AARotMat(2,1);
      StressRotMat(4, 2) =   AARotMat(0,2) * AARotMat(2,2);
      StressRotMat(4, 3) = ( AARotMat(0,1) * AARotMat(2,2) + AARotMat(0,2) * AARotMat(2,1) );
      StressRotMat(4, 4) = ( AARotMat(0,0) * AARotMat(2,2) + AARotMat(0,2) * AARotMat(2,0) );
      StressRotMat(4, 5) = ( AARotMat(0,0) * AARotMat(2,1) + AARotMat(0,1) * AARotMat(2,0) );

      StressRotMat(5, 0) =   AARotMat(0,0) * AARotMat(2,0);
      StressRotMat(5, 1) =   AARotMat(0,1) * AARotMat(1,1);
      StressRotMat(5, 2) =   AARotMat(0,2) * AARotMat(1,2);
      StressRotMat(5, 3) = ( AARotMat(0,1) * AARotMat(1,2) + AARotMat(0,2) * AARotMat(1,1) );
      StressRotMat(5, 4) = ( AARotMat(0,0) * AARotMat(1,2) + AARotMat(0,2) * AARotMat(1,0) );
      StressRotMat(5, 5) = ( AARotMat(0,0) * AARotMat(1,1) + AARotMat(0,1) * AARotMat(1,0) );
      */
      // -------------
      StressRotMat(0, 0) =       AARotMat(0,0) * AARotMat(0,0);
      StressRotMat(0, 1) =       AARotMat(1,0) * AARotMat(1,0);
      StressRotMat(0, 2) =       AARotMat(2,0) * AARotMat(2,0);
      StressRotMat(0, 3) = 2.0 * AARotMat(1,0) * AARotMat(0,0);
      StressRotMat(0, 4) = 2.0 * AARotMat(2,0) * AARotMat(1,0);
      StressRotMat(0, 5) = 2.0 * AARotMat(0,0) * AARotMat(2,0);

      StressRotMat(1, 0) =       AARotMat(0,1) * AARotMat(0,1);
      StressRotMat(1, 1) =       AARotMat(1,1) * AARotMat(1,1);
      StressRotMat(1, 2) =       AARotMat(2,1) * AARotMat(2,1);
      StressRotMat(1, 3) = 2.0 * AARotMat(1,1) * AARotMat(0,1);
      StressRotMat(1, 4) = 2.0 * AARotMat(2,1) * AARotMat(1,1);
      StressRotMat(1, 5) = 2.0 * AARotMat(0,1) * AARotMat(2,1);

      StressRotMat(2, 0) =       AARotMat(0,2) * AARotMat(0,2);
      StressRotMat(2, 1) =       AARotMat(1,2) * AARotMat(1,2);
      StressRotMat(2, 2) =       AARotMat(2,2) * AARotMat(2,2);
      StressRotMat(2, 3) = 2.0 * AARotMat(1,2) * AARotMat(0,2);
      StressRotMat(2, 4) = 2.0 * AARotMat(2,2) * AARotMat(1,2);
      StressRotMat(2, 5) = 2.0 * AARotMat(0,2) * AARotMat(2,2);

      StressRotMat(3, 0) =   AARotMat(0,1) * AARotMat(0,0);
      StressRotMat(3, 1) =   AARotMat(1,1) * AARotMat(1,0);
      StressRotMat(3, 2) =   AARotMat(2,1) * AARotMat(2,0);
      StressRotMat(3, 3) = ( AARotMat(1,1) * AARotMat(0,0) + AARotMat(0,1) * AARotMat(1,0) );
      StressRotMat(3, 4) = ( AARotMat(2,1) * AARotMat(1,0) + AARotMat(1,1) * AARotMat(2,0) );
      StressRotMat(3, 5) = ( AARotMat(0,1) * AARotMat(2,0) + AARotMat(2,1) * AARotMat(0,0) );

      StressRotMat(4, 0) =   AARotMat(0,2) * AARotMat(0,1);
      StressRotMat(4, 1) =   AARotMat(1,2) * AARotMat(1,1);
      StressRotMat(4, 2) =   AARotMat(2,2) * AARotMat(2,1);
      StressRotMat(4, 3) = ( AARotMat(1,2) * AARotMat(0,1) + AARotMat(0,2) * AARotMat(1,1) );
      StressRotMat(4, 4) = ( AARotMat(2,2) * AARotMat(1,1) + AARotMat(1,2) * AARotMat(2,1) );
      StressRotMat(4, 5) = ( AARotMat(0,2) * AARotMat(2,1) + AARotMat(2,2) * AARotMat(0,1) );

      StressRotMat(5, 0) =   AARotMat(0,0) * AARotMat(0,2);
      StressRotMat(5, 1) =   AARotMat(1,0) * AARotMat(1,2);
      StressRotMat(5, 2) =   AARotMat(2,0) * AARotMat(2,2);
      StressRotMat(5, 3) = ( AARotMat(1,0) * AARotMat(0,2) + AARotMat(0,0) * AARotMat(1,2) );
      StressRotMat(5, 4) = ( AARotMat(2,0) * AARotMat(1,2) + AARotMat(1,0) * AARotMat(2,2) );
      StressRotMat(5, 5) = ( AARotMat(0,0) * AARotMat(2,2) + AARotMat(2,0) * AARotMat(0,2) );

    }

  };

  /**
   * \brief Function to Calculate Strain Transformation Matrix<br>
   * This function computes the strain transformation Matrix at a give axis and angle of rotation <br>
   * One can also output the axis/angle rotational Matrix
   *
   *\param AxVector A vector representing the axis of rotation
   *\param AxAngle Angle of rotation along the axis (in radians)
   */
  struct StrainTransformation {

    ublas::matrix<double> StrainRotMat;
    ublas::matrix<double> AARotMat;

    StrainTransformation(double *AxVector, double AxAngle){

      AARotMat = ublas::zero_matrix<FieldData>(3,3);
      AxisAngleRotationalMatrix RotMat(&AxVector[0], AxAngle);
      AARotMat=RotMat.AARotMat;

      StrainRotMat = ublas::zero_matrix<FieldData>(6,6);

//      StrainRotMat(0, 0) = AARotMat(0,0) * AARotMat(0,0);
//      StrainRotMat(0, 1) = AARotMat(0,1) * AARotMat(0,1);
//      StrainRotMat(0, 2) = AARotMat(0,2) * AARotMat(0,2);
//      StrainRotMat(0, 3) = AARotMat(0,1) * AARotMat(0,2);
//      StrainRotMat(0, 4) = AARotMat(0,0) * AARotMat(0,2);
//      StrainRotMat(0, 5) = AARotMat(0,0) * AARotMat(0,1);
//
//      StrainRotMat(1, 0) = AARotMat(1,0) * AARotMat(1,0);
//      StrainRotMat(1, 1) = AARotMat(1,1) * AARotMat(1,1);
//      StrainRotMat(1, 2) = AARotMat(1,2) * AARotMat(1,2);
//      StrainRotMat(1, 3) = AARotMat(1,1) * AARotMat(1,2);
//      StrainRotMat(1, 4) = AARotMat(1,0) * AARotMat(1,2);
//      StrainRotMat(1, 5) = AARotMat(1,0) * AARotMat(1,1);
//
//      StrainRotMat(2, 0) = AARotMat(2,0) * AARotMat(2,0);
//      StrainRotMat(2, 1) = AARotMat(2,1) * AARotMat(2,1);
//      StrainRotMat(2, 2) = AARotMat(2,2) * AARotMat(2,2);
//      StrainRotMat(2, 3) = AARotMat(2,1) * AARotMat(2,2);
//      StrainRotMat(2, 4) = AARotMat(2,0) * AARotMat(2,2);
//      StrainRotMat(2, 5) = AARotMat(2,0) * AARotMat(2,1);
//
//      StrainRotMat(3, 0) = 2.0 * AARotMat(1,0) * AARotMat(2,0);
//      StrainRotMat(3, 1) = 2.0 * AARotMat(1,1) * AARotMat(2,1);
//      StrainRotMat(3, 2) = 2.0 * AARotMat(1,2) * AARotMat(2,2);
//      StrainRotMat(3, 3) =     ( AARotMat(1,1) * AARotMat(2,2) + AARotMat(1,2) * AARotMat(2,1) );
//      StrainRotMat(3, 4) =     ( AARotMat(1,0) * AARotMat(2,2) + AARotMat(1,2) * AARotMat(2,0) );
//      StrainRotMat(3, 5) =     ( AARotMat(1,0) * AARotMat(2,1) + AARotMat(1,1) * AARotMat(2,0) );
//
//      StrainRotMat(4, 0) = 2.0 * AARotMat(0,0) * AARotMat(2,0);
//      StrainRotMat(4, 1) = 2.0 * AARotMat(0,1) * AARotMat(2,1);
//      StrainRotMat(4, 2) = 2.0 * AARotMat(0,2) * AARotMat(2,2);
//      StrainRotMat(4, 3) =     ( AARotMat(0,1) * AARotMat(2,2) + AARotMat(0,2) * AARotMat(2,1) );
//      StrainRotMat(4, 4) =     ( AARotMat(0,0) * AARotMat(2,2) + AARotMat(0,2) * AARotMat(2,0) );
//      StrainRotMat(4, 5) =     ( AARotMat(0,0) * AARotMat(2,1) + AARotMat(0,1) * AARotMat(2,0) );
//
//      StrainRotMat(5, 0) = 2.0 * AARotMat(0,0) * AARotMat(1,0);
//      StrainRotMat(5, 1) = 2.0 * AARotMat(0,1) * AARotMat(1,1);
//      StrainRotMat(5, 2) = 2.0 * AARotMat(0,2) * AARotMat(1,2);
//      StrainRotMat(5, 3) =     ( AARotMat(0,1) * AARotMat(1,2) + AARotMat(0,2) * AARotMat(1,1) );
//      StrainRotMat(5, 4) =     ( AARotMat(0,0) * AARotMat(1,2) + AARotMat(0,2) * AARotMat(1,0) );
//      StrainRotMat(5, 5) =     ( AARotMat(0,0) * AARotMat(1,1) + AARotMat(0,1) * AARotMat(1,0) );

      // ------
      StrainRotMat(0, 0) = AARotMat(0,0) * AARotMat(0,0);
      StrainRotMat(0, 1) = AARotMat(1,0) * AARotMat(1,0);
      StrainRotMat(0, 2) = AARotMat(2,0) * AARotMat(2,0);
      StrainRotMat(0, 3) = AARotMat(1,0) * AARotMat(0,0);
      StrainRotMat(0, 4) = AARotMat(2,0) * AARotMat(1,0);
      StrainRotMat(0, 5) = AARotMat(0,0) * AARotMat(2,0);

      StrainRotMat(1, 0) = AARotMat(0,1) * AARotMat(0,1);
      StrainRotMat(1, 1) = AARotMat(1,1) * AARotMat(1,1);
      StrainRotMat(1, 2) = AARotMat(2,1) * AARotMat(2,1);
      StrainRotMat(1, 3) = AARotMat(1,1) * AARotMat(0,1);
      StrainRotMat(1, 4) = AARotMat(2,1) * AARotMat(1,1);
      StrainRotMat(1, 5) = AARotMat(0,1) * AARotMat(2,1);

      StrainRotMat(2, 0) = AARotMat(0,2) * AARotMat(0,2);
      StrainRotMat(2, 1) = AARotMat(1,2) * AARotMat(1,2);
      StrainRotMat(2, 2) = AARotMat(2,2) * AARotMat(2,2);
      StrainRotMat(2, 3) = AARotMat(1,2) * AARotMat(0,2);
      StrainRotMat(2, 4) = AARotMat(2,2) * AARotMat(1,2);
      StrainRotMat(2, 5) = AARotMat(0,2) * AARotMat(2,2);

      StrainRotMat(3, 0) = 2.0 * AARotMat(0,1) * AARotMat(0,0);
      StrainRotMat(3, 1) = 2.0 * AARotMat(1,1) * AARotMat(1,0);
      StrainRotMat(3, 2) = 2.0 * AARotMat(2,1) * AARotMat(2,0);
      StrainRotMat(3, 3) =     ( AARotMat(1,1) * AARotMat(0,0) + AARotMat(0,1) * AARotMat(1,0) );
      StrainRotMat(3, 4) =     ( AARotMat(2,1) * AARotMat(1,0) + AARotMat(1,1) * AARotMat(2,0) );
      StrainRotMat(3, 5) =     ( AARotMat(0,1) * AARotMat(2,0) + AARotMat(2,1) * AARotMat(0,0) );

      StrainRotMat(4, 0) = 2.0 * AARotMat(0,2) * AARotMat(0,1);
      StrainRotMat(4, 1) = 2.0 * AARotMat(1,2) * AARotMat(1,1);
      StrainRotMat(4, 2) = 2.0 * AARotMat(2,2) * AARotMat(2,1);
      StrainRotMat(4, 3) =     ( AARotMat(1,2) * AARotMat(0,1) + AARotMat(0,2) * AARotMat(1,1) );
      StrainRotMat(4, 4) =     ( AARotMat(2,2) * AARotMat(1,1) + AARotMat(1,2) * AARotMat(2,1) );
      StrainRotMat(4, 5) =     ( AARotMat(0,2) * AARotMat(2,1) + AARotMat(2,2) * AARotMat(0,1) );

      StrainRotMat(5, 0) = 2.0 * AARotMat(0,0) * AARotMat(0,2);
      StrainRotMat(5, 1) = 2.0 * AARotMat(1,0) * AARotMat(1,2);
      StrainRotMat(5, 2) = 2.0 * AARotMat(2,0) * AARotMat(2,2);
      StrainRotMat(5, 3) =     ( AARotMat(1,0) * AARotMat(0,2) + AARotMat(0,0) * AARotMat(1,2) );
      StrainRotMat(5, 4) =     ( AARotMat(2,0) * AARotMat(1,2) + AARotMat(1,0) * AARotMat(2,2) );
      StrainRotMat(5, 5) =     ( AARotMat(0,0) * AARotMat(2,2) + AARotMat(2,0) * AARotMat(0,2) );


    }

  };

  /**
   * \brief Function to build up the Transverse Isotropic Stiffness Matrix and rotated when necessary
   *
   * \param E_p Young's Modulus xy plane
   * \param E_z Young's Modulus z-axis
   * \param nu_p Poisson's Ratio xy plane
   * \param nu_pz Poisson's Ration z-axis
   * \param G_zp Shear Modulus z-direction
   * \param noAA Number of Axis-Angle Rotations (Could be 0 for no rotation)
   * \param AxVector An array containing all the axes of rotation ( ex: AxVector[6] = {(1st Axis) 1,0,0 , (2nd Axis) 0,1,0} )
   * \param AxAngle An array containg all the angles of rotation ( ex: AxAngle[2] = {(1st Angle) 0.5*M_PI, (2nd Angle) 0.25*M_PI} )
   */
  struct TranIsotropicAxisAngleRotElasticFEMethod: public ElasticFEMethod {

    int noAA;
    double *AxVector, *AxAngle;
    bool propeties_from_BLOCKSET_MAT_ELASTICSET;

    TranIsotropicAxisAngleRotElasticFEMethod(FieldInterface& _mField,Mat &_Aij,Vec _D,Vec _F,
                                             int _noAA, double *_AxVector, double *_AxAngle):
    ElasticFEMethod(_mField,_Aij,_D,_F,0,0), noAA(_noAA), AxVector(_AxVector), AxAngle(_AxAngle)  {
      propeties_from_BLOCKSET_MAT_ELASTICSET = false;
      for(_IT_CUBITMESHSETS_BY_BCDATA_TYPE_FOR_LOOP_(mField,BLOCKSET|MAT_ELASTICSET,it)) {
        propeties_from_BLOCKSET_MAT_ELASTICSET = true;
      }
    }

//--------------------------------------------------------------------------------------------------------------------------------------------------//
    PetscErrorCode GetMatParameters(double *_E_p, double *_E_z, double *_nu_p, double *_nu_pz, double *_G_zp) {
      PetscFunctionBegin;

      if(propeties_from_BLOCKSET_MAT_ELASTICSET) {
        EntityHandle ent = fePtr->get_ent();
        for(_IT_CUBITMESHSETS_BY_BCDATA_TYPE_FOR_LOOP_(mField,BLOCKSET|MAT_ELASTICSET,it)) {

          if(it->get_name().compare(0,20,"MAT_ELASTIC_TRANSISO") == 0) {

            Mat_Elastic_TransIso mydata;
            ierr = it->get_attribute_data_structure(mydata); CHKERRQ(ierr);

            Range meshsets;
            rval = moab.get_entities_by_type(it->meshset,MBENTITYSET,meshsets,true); CHKERR_PETSC(rval);
            meshsets.insert(it->meshset);
            for(Range::iterator mit = meshsets.begin();mit != meshsets.end(); mit++) {
              if( moab.contains_entities(*mit,&ent,1) ) {
                *_E_p = mydata.data.Youngp;
                *_E_z = mydata.data.Youngz;
                *_nu_p = mydata.data.Poissonp;
                *_nu_pz = mydata.data.Poissonpz;
                *_G_zp = mydata.data.Shearzp;
                PetscFunctionReturn(0);
              }
            }
          }
        }

        SETERRQ(PETSC_COMM_SELF,1,
                "Element is not in transervely isotropic block, however you run linear transversely isotropic analysis with that element\n"
                "top tip: check if you update block sets after mesh refinments or interface insertion");

      }

      PetscFunctionReturn(0);
    }
//--------------------------------------------------------------------------------------------------------------------------------------------------//
    PetscErrorCode calculateD(double _E_p,double _E_z, double _nu_p,double _nu_pz, double _G_zp) {
      PetscFunctionBegin;

      ///Get Stiffness Matrix
      ublas::symmetric_matrix<FieldData,ublas::upper> StiffnessMatrix;
      StiffnessMatrix.resize(6);
      StiffnessMatrix.clear();
      TransverseIsotropicStiffnessMatrix TranIsoMat(_nu_p, _nu_pz, _E_p, _E_z, _G_zp);
      StiffnessMatrix=TranIsoMat.StiffnessMatrix;

      ///Rotating the Stiffness matrix according a set of axes of rotations and their respective angle

      int noOfRotations = noAA; //Number of Rotations
      double negAxAngle[noOfRotations];
      for (int aa=0; aa<noOfRotations; aa++) negAxAngle[aa]=-AxAngle[aa];

      ublas::matrix<double> DummyMatrix,DummyMatrix2;
      DummyMatrix = ublas::zero_matrix<FieldData>(6,6);
      DummyMatrix = StiffnessMatrix;

      ///Rotating Stiffness over a number of axis/angle rotations
      for (int aa=0; aa<noOfRotations; aa++) {

        StressTransformation StressRotMat(&AxVector[3*aa], AxAngle[aa]);
        StrainTransformation invStrainRotMat(&AxVector[3*aa], negAxAngle[aa]);

        ublas::matrix<double> TrpMatrixStress;
        TrpMatrixStress = ublas::zero_matrix<FieldData>(6,6);
        TrpMatrixStress=StressRotMat.StressRotMat;

        ublas::matrix<double> TrpMatrixInvStrain;
        TrpMatrixInvStrain = ublas::zero_matrix<FieldData>(6,6);
        TrpMatrixInvStrain=invStrainRotMat.StrainRotMat;

        DummyMatrix2 = ublas::zero_matrix<FieldData>(6,6);
        ublas::matrix< FieldData > dummyA = prod( DummyMatrix , TrpMatrixInvStrain );
        DummyMatrix2 = prod(TrpMatrixStress,dummyA);
        DummyMatrix = ublas::zero_matrix<FieldData>(6,6);
        DummyMatrix = DummyMatrix2;

      }

      D.resize(6,6);
      D.clear();
      D = DummyMatrix;

      PetscFunctionReturn(0);
    }
//--------------------------------------------------------------------------------------------------------------------------------------------------//
    PetscErrorCode Fint() {
      PetscFunctionBegin;

      try {

        //Higher order approximation of geometry
        ierr = GetHierarchicalGeometryApproximation(invH,detH); CHKERRQ(ierr);

        double _E_p, _E_z, _nu_p, _nu_pz, _G_zp;
        ierr = GetMatParameters(&_E_p,&_E_z,&_nu_p,&_nu_pz,&_G_zp); CHKERRQ(ierr);
        ierr = calculateD(_E_p,_E_z,_nu_p,_nu_pz,_G_zp); CHKERRQ(ierr);

        //Gradient at Gauss points;
        vector< ublas::matrix< FieldData > > GradU_at_GaussPt;
        ierr = GetGaussDiffDataVector("DISPLACEMENT",GradU_at_GaussPt); CHKERRQ(ierr);
        unsigned int g_dim = g_NTET.size()/4;
        assert(GradU_at_GaussPt.size() == g_dim);
        NOT_USED(g_dim);
        vector< ublas::matrix< FieldData > >::iterator viit = GradU_at_GaussPt.begin();
        int gg = 0;
        for(;viit!=GradU_at_GaussPt.end();viit++,gg++) {
          try {
            ublas::matrix< FieldData > GradU = *viit;
            if(!invH.empty()) {
              //GradU =
              //[ dU/dChi1 dU/dChi2 dU/dChi3 ]
              //[ dV/dChi1 dV/dChi2 dU/dChi3 ]
              //[ dW/dChi1 dW/dChi2 dW/dChi3 ]
              //H =
              //[ dX1/dChi1 dX1/dChi2 dX1/dChi3 ]
              //[ dX2/dChi1 dX2/dChi2 dX2/dChi3 ]
              //[ dX3/dChi1 dX3/dChi2 dX3/dChi3 ]
              //invH =
              //[ dChi1/dX1 dChi1/dX2 dChi1/dX3 ]
              //[ dChi2/dX1 dChi2/dX2 dChi2/dX3 ]
              //[ dChi3/dX1 dChi3/dX2 dChi3/dX3 ]
              //GradU =
              //[ dU/dX1 dU/dX2 dU/dX3 ]
              //[ dV/dX1 dV/dX2 dV/dX3 ] = GradU * invH
              //[ dW/dX1 dW/dX2 dW/dX3 ]
              GradU = prod( GradU, invH[gg] );
            }
            ublas::matrix< FieldData > Strain = 0.5*( GradU + trans(GradU) );
            ublas::vector< FieldData > VoightStrain(6);
            VoightStrain[0] = Strain(0,0);
            VoightStrain[1] = Strain(1,1);
            VoightStrain[2] = Strain(2,2);
            VoightStrain[3] = 2*Strain(0,1);
            VoightStrain[4] = 2*Strain(1,2);
            VoightStrain[5] = 2*Strain(2,0);
            double w = V*G_TET_W[gg];
            ublas::vector<FieldData> VoightStress = prod(w*D,VoightStrain);
            //BT * VoigtStress
            f_int.resize(row_mat);
            for(int rr = 0;rr<row_mat;rr++) {
              if(RowGlob[rr].size()==0) continue;
              ublas::matrix<FieldData> &B = (rowBMatrices[rr])[gg];
              if(gg == 0) {
                f_int[rr] = prod( trans(B), VoightStress );
              } else {
                f_int[rr] += prod( trans(B), VoightStress );
              }
            }
          } catch (const std::exception& ex) {
            ostringstream ss;
            ss << "thorw in method: " << ex.what() << " at line " << __LINE__ << " in file " << __FILE__ << endl;
            SETERRQ(PETSC_COMM_SELF,1,ss.str().c_str());
          }
        }

      } catch (const std::exception& ex) {
        ostringstream ss;
        ss << "thorw in method: " << ex.what() << " at line " << __LINE__ << " in file " << __FILE__ << endl;
        SETERRQ(PETSC_COMM_SELF,1,ss.str().c_str());
      }

      PetscFunctionReturn(0);
    }
//--------------------------------------------------------------------------------------------------------------------------------------------------//
    PetscErrorCode Stiffness() {
      PetscFunctionBegin;

      double _E_p, _E_z, _nu_p, _nu_pz, _G_zp;
      ierr = GetMatParameters(&_E_p,&_E_z,&_nu_p,&_nu_pz,&_G_zp); CHKERRQ(ierr);
      ierr = calculateD(_E_p,_E_z,_nu_p,_nu_pz,_G_zp); CHKERRQ(ierr);

      K.resize(row_mat,col_mat);
      int g_dim = g_NTET.size()/4;
      for(int rr = 0;rr<row_mat;rr++) {
        if(RowGlob[rr].size()==0) continue;
        for(int gg = 0;gg<g_dim;gg++) {
          ublas::matrix<FieldData> &row_Mat = (rowBMatrices[rr])[gg];
          double w = V*G_TET_W[gg];
          if(detH.size()>0) {
            w *= detH[gg];
          }
          BD.resize(6,row_Mat.size2());
          //ublas::noalias(BD) = prod( w*D,row_Mat );
          cblas_dsymm(CblasRowMajor,CblasLeft,CblasUpper,
                      BD.size1(),BD.size2(),
                      w,&*D.data().begin(),D.size2(),
                      &*row_Mat.data().begin(),row_Mat.size2(),
                      0.,&*BD.data().begin(),BD.size2());
          for(int cc = rr;cc<col_mat;cc++) {
            if(ColGlob[cc].size()==0) continue;
            ublas::matrix<FieldData> &col_Mat = (colBMatrices[cc])[gg];
            if(gg == 0) {
              K(rr,cc).resize(BD.size2(),col_Mat.size2());
              //ublas::noalias(K(rr,cc)) = prod(trans(BD) , col_Mat ); // int BT*D*B
              cblas_dgemm(CblasRowMajor,CblasTrans,CblasNoTrans,
                          BD.size2(),col_Mat.size2(),BD.size1(),
                          1.,&*BD.data().begin(),BD.size2(),
                          &*col_Mat.data().begin(),col_Mat.size2(),
                          0.,&*K(rr,cc).data().begin(),K(rr,cc).size2());
            } else {
              //ublas::noalias(K(rr,cc)) += prod(trans(BTD) , col_Mat ); // int BT*D*B
              cblas_dgemm(CblasRowMajor,CblasTrans,CblasNoTrans,
                          BD.size2(),col_Mat.size2(),BD.size1(),
                          1.,&*BD.data().begin(),BD.size2(),
                          &*col_Mat.data().begin(),col_Mat.size2(),
                          1.,&*K(rr,cc).data().begin(),K(rr,cc).size2());
            }
          }
        }
      }
      PetscFunctionReturn(0);
    }

  };

  //==================================================================================================================================================//

  /**
   * \brief Function to build up the Transverse Isotropic Stiffness Matrix and rotate it according to the fibre direction.<br>
   * This class, same as \ref<TranIsotropicAxisAngleRotElasticFEMethod>TranIsotropicAxisAngleRotElasticFEMethod
   * but will rotate the stiffness matrix using fibre direction which are computed from phi (that are computed for simplest potential flow problem). <br>
   * A fibre direction would be computed for every gauss point and hence, stiffness matrix is rotated at every gauss point. <br>
   *
   * \param E_p Young's Modulus xy plane
   * \param E_z Young's Modulus z-axis
   * \param nu_p Poisson's Ratio xy plane
   * \param nu_pz Poisson's Ration z-axis
   * \param G_zp Shear Modulus z-direction
   */
  struct TranIsotropicFibreDirRotElasticFEMethod: public ElasticFEMethod {

    Tag th_fibre_dir;
    bool propeties_from_BLOCKSET_MAT_ELASTICSET;

    TranIsotropicFibreDirRotElasticFEMethod(FieldInterface& _mField,Mat &_Aij,Vec _D,Vec _F, string _field_name = "DISPLACEMENT"):
      ElasticFEMethod(_mField,_Aij,_D,_F,0,0,_field_name) {

      propeties_from_BLOCKSET_MAT_ELASTICSET = false;
      for(_IT_CUBITMESHSETS_BY_BCDATA_TYPE_FOR_LOOP_(mField,BLOCKSET|MAT_ELASTICSET,it)) {
        propeties_from_BLOCKSET_MAT_ELASTICSET = true;
      }

      double def_VAL2[3] = {0,0,0};
      rval = moab.tag_get_handle( "POT_FLOW_FIBRE_DIR",3,MB_TYPE_DOUBLE,th_fibre_dir,MB_TAG_CREAT|MB_TAG_SPARSE,&def_VAL2); CHKERR_THROW(rval);
      //cout<<"The fibre direction"<<endl;

      };

    //--------------------------------------------------------------------------------------------------------------------------------------------------//
    PetscErrorCode GetMatParameters(double *_E_p, double *_E_z, double *_nu_p, double *_nu_pz, double *_G_zp) {
      PetscFunctionBegin;

      if(propeties_from_BLOCKSET_MAT_ELASTICSET) {
        EntityHandle ent = fePtr->get_ent();
        for(_IT_CUBITMESHSETS_BY_BCDATA_TYPE_FOR_LOOP_(mField,BLOCKSET|MAT_ELASTICSET,it)) {

          if(it->get_name().compare(0,20,"MAT_ELASTIC_TRANSISO") == 0) {

            Mat_Elastic_TransIso mydata;
            ierr = it->get_attribute_data_structure(mydata); CHKERRQ(ierr);

            Range meshsets;
            rval = moab.get_entities_by_type(it->meshset,MBENTITYSET,meshsets,true); CHKERR_PETSC(rval);
            meshsets.insert(it->meshset);
            for(Range::iterator mit = meshsets.begin();mit != meshsets.end(); mit++) {
              if( moab.contains_entities(*mit,&ent,1) ) {
                *_E_p = mydata.data.Youngp;
                *_E_z = mydata.data.Youngz;
                *_nu_p = mydata.data.Poissonp;
                *_nu_pz = mydata.data.Poissonpz;
                *_G_zp = mydata.data.Shearzp;
                PetscFunctionReturn(0);
              }
            }
          }
        }

        SETERRQ(PETSC_COMM_SELF,1,
                "Element is not in transervely isotropic block, however you run linear transversely isotropic analysis with that element\n"
                "top tip: check if you update block sets after mesh refinments or interface insertion");

      }

      PetscFunctionReturn(0);
    }
    //--------------------------------------------------------------------------------------------------------------------------------------------------//
    vector< ublas::matrix<FieldData> > D_At_GaussPoint;
    PetscErrorCode calculateD(double _E_p,double _E_z, double _nu_p,double _nu_pz, double _G_zp) {
      PetscFunctionBegin;

      ///Get Stiffness Matrix
      ublas::symmetric_matrix<FieldData,ublas::upper> StiffnessMatrix;
      StiffnessMatrix.resize(6);
      StiffnessMatrix.clear();
      TransverseIsotropicStiffnessMatrix TranIsoMat(_nu_p,_nu_pz,_E_p,_E_z,_G_zp);
      StiffnessMatrix=TranIsoMat.StiffnessMatrix;

      ///Rotating the Stiffness matrix according a set of axes of rotations and their respective angle

      D_At_GaussPoint.resize(coords_at_Gauss_nodes.size());

      vector< ublas::matrix< double > > normalized_phi;
      normalized_phi.resize(coords_at_Gauss_nodes.size());
      ierr = ComputeFibreDirection(normalized_phi); CHKERRQ(ierr);

      for(unsigned int gg=0;gg<coords_at_Gauss_nodes.size();gg++){

        int noOfRotations = 1; //Number of Rotations

        double zVec[3]={0.0,0.0,1.0};
        double AxVector[3]={normalized_phi[gg](0,1)*zVec[2]-normalized_phi[gg](0,2)*zVec[1] , normalized_phi[gg](0,2)*zVec[0]-normalized_phi[gg](0,0)*zVec[2] , normalized_phi[gg](0,0)*zVec[1]-normalized_phi[gg](0,1)*zVec[0]};
        double AxAngle[1]= {asin((sqrt(pow(AxVector[0],2)+pow(AxVector[1],2)+pow(AxVector[2],2)))/(sqrt(pow(normalized_phi[gg](0,0),2)+pow(normalized_phi[gg](0,1),2)+pow(normalized_phi[gg](0,2),2)))*(sqrt(pow(zVec[0],2)+pow(zVec[1],2)+pow(zVec[2],2))))};

        double negAxAngle[noOfRotations];
        for (int aa=0; aa<noOfRotations; aa++) negAxAngle[aa]=-AxAngle[aa];

        ublas::matrix<double> DummyMatrix,DummyMatrix2;
        DummyMatrix = ublas::zero_matrix<FieldData>(6,6);
        DummyMatrix = StiffnessMatrix;

        ///Rotating Stiffness over a number of axis/angle rotations
        for (int aa=0; aa<noOfRotations; aa++) {

          StressTransformation StressRotMat(&AxVector[3*aa], AxAngle[aa]);
          StrainTransformation invStrainRotMat(&AxVector[3*aa], negAxAngle[aa]);

          ublas::matrix<double> TrpMatrixStress;
          TrpMatrixStress = ublas::zero_matrix<FieldData>(6,6);
          TrpMatrixStress=StressRotMat.StressRotMat;

          ublas::matrix<double> TrpMatrixInvStrain;
          TrpMatrixInvStrain = ublas::zero_matrix<FieldData>(6,6);
          TrpMatrixInvStrain=invStrainRotMat.StrainRotMat;

          DummyMatrix2 = ublas::zero_matrix<FieldData>(6,6);
          ublas::matrix< FieldData > dummyA = prod( DummyMatrix , TrpMatrixInvStrain );
          DummyMatrix2 = prod(TrpMatrixStress,dummyA);
          DummyMatrix = ublas::zero_matrix<FieldData>(6,6);
          DummyMatrix = DummyMatrix2;
        }

        D_At_GaussPoint[gg].resize(6,6);
        D_At_GaussPoint[gg].clear();
        D_At_GaussPoint[gg] = DummyMatrix;
      }

      PetscFunctionReturn(0);
    }
    //--------------------------------------------------------------------------------------------------------------------------------------------------//
    PetscErrorCode Fint() {
      PetscFunctionBegin;

      try {

        //Higher order approximation of geometry
        ierr = GetHierarchicalGeometryApproximation(invH,detH); CHKERRQ(ierr);

        double _E_p, _E_z, _nu_p, _nu_pz, _G_zp;
        ierr = GetMatParameters(&_E_p,&_E_z,&_nu_p,&_nu_pz,&_G_zp); CHKERRQ(ierr);
        ierr = calculateD(_E_p,_E_z,_nu_p,_nu_pz,_G_zp); CHKERRQ(ierr);

        //Gradient at Gauss points;
        vector< ublas::matrix< FieldData > > GradU_at_GaussPt;
        ierr = GetGaussDiffDataVector(fieldName,GradU_at_GaussPt); CHKERRQ(ierr);
        unsigned int g_dim = g_NTET.size()/4;
        assert(GradU_at_GaussPt.size() == g_dim);
        NOT_USED(g_dim);
        vector< ublas::matrix< FieldData > >::iterator viit = GradU_at_GaussPt.begin();
        int gg = 0;
        for(;viit!=GradU_at_GaussPt.end();viit++,gg++) {
          try {
            ublas::matrix< FieldData > GradU = *viit;
            if(!invH.empty()) {
              //GradU =
              //[ dU/dChi1 dU/dChi2 dU/dChi3 ]
              //[ dV/dChi1 dV/dChi2 dU/dChi3 ]
              //[ dW/dChi1 dW/dChi2 dW/dChi3 ]
              //H =
              //[ dX1/dChi1 dX1/dChi2 dX1/dChi3 ]
              //[ dX2/dChi1 dX2/dChi2 dX2/dChi3 ]
              //[ dX3/dChi1 dX3/dChi2 dX3/dChi3 ]
              //invH =
              //[ dChi1/dX1 dChi1/dX2 dChi1/dX3 ]
              //[ dChi2/dX1 dChi2/dX2 dChi2/dX3 ]
              //[ dChi3/dX1 dChi3/dX2 dChi3/dX3 ]
              //GradU =
              //[ dU/dX1 dU/dX2 dU/dX3 ]
              //[ dV/dX1 dV/dX2 dV/dX3 ] = GradU * invH
              //[ dW/dX1 dW/dX2 dW/dX3 ]
              GradU = prod( GradU, invH[gg] );
            }
            ublas::matrix< FieldData > Strain = 0.5*( GradU + trans(GradU) );
            ublas::vector< FieldData > VoightStrain(6);
            VoightStrain[0] = Strain(0,0);
            VoightStrain[1] = Strain(1,1);
            VoightStrain[2] = Strain(2,2);
            VoightStrain[3] = 2*Strain(0,1);
            VoightStrain[4] = 2*Strain(1,2);
            VoightStrain[5] = 2*Strain(2,0);
            double w = V*G_TET_W[gg];
            ublas::vector<FieldData> VoightStress = prod(w*D_At_GaussPoint[gg],VoightStrain);
            //BT * VoigtStress
            f_int.resize(row_mat);
            for(int rr = 0;rr<row_mat;rr++) {
              if(RowGlob[rr].size()==0) continue;
              ublas::matrix<FieldData> &B = (rowBMatrices[rr])[gg];
              if(gg == 0) {
                f_int[rr] = prod( trans(B), VoightStress );
              } else {
                f_int[rr] += prod( trans(B), VoightStress );
              }
            }
          } catch (const std::exception& ex) {
            ostringstream ss;
            ss << "thorw in method: " << ex.what() << " at line " << __LINE__ << " in file " << __FILE__ << endl;
            SETERRQ(PETSC_COMM_SELF,1,ss.str().c_str());
          }
        }

      } catch (const std::exception& ex) {
        ostringstream ss;
        ss << "thorw in method: " << ex.what() << " at line " << __LINE__ << " in file " << __FILE__ << endl;
        SETERRQ(PETSC_COMM_SELF,1,ss.str().c_str());
      }

      PetscFunctionReturn(0);
    }
    //--------------------------------------------------------------------------------------------------------------------------------------------------//
    PetscErrorCode Stiffness() {
      PetscFunctionBegin;

      double _E_p, _E_z, _nu_p, _nu_pz, _G_zp;
      ierr = GetMatParameters(&_E_p,&_E_z,&_nu_p,&_nu_pz,&_G_zp); CHKERRQ(ierr);
      ierr = calculateD(_E_p,_E_z,_nu_p,_nu_pz,_G_zp); CHKERRQ(ierr);

      K.resize(row_mat,col_mat);
      int g_dim = g_NTET.size()/4;
      for(int rr = 0;rr<row_mat;rr++) {
        if(RowGlob[rr].size()==0) continue;
        for(int gg = 0;gg<g_dim;gg++) {
          ublas::matrix<FieldData> &row_Mat = (rowBMatrices[rr])[gg];
          double w = V*G_TET_W[gg];
          if(detH.size()>0) {
            w *= detH[gg];
          }
          BD.resize(6,row_Mat.size2());
          //ublas::noalias(BD) = prod( w*D,row_Mat );
          cblas_dsymm(CblasRowMajor,CblasLeft,CblasUpper,
                      BD.size1(),BD.size2(),
                      w,&*D_At_GaussPoint[gg].data().begin(),D_At_GaussPoint[gg].size2(),
                      &*row_Mat.data().begin(),row_Mat.size2(),
                      0.,&*BD.data().begin(),BD.size2());
          for(int cc = rr;cc<col_mat;cc++) {
            if(ColGlob[cc].size()==0) continue;
            ublas::matrix<FieldData> &col_Mat = (colBMatrices[cc])[gg];
            if(gg == 0) {
              K(rr,cc).resize(BD.size2(),col_Mat.size2());
              //ublas::noalias(K(rr,cc)) = prod(trans(BD) , col_Mat ); // int BT*D*B
              cblas_dgemm(CblasRowMajor,CblasTrans,CblasNoTrans,
                          BD.size2(),col_Mat.size2(),BD.size1(),
                          1.,&*BD.data().begin(),BD.size2(),
                          &*col_Mat.data().begin(),col_Mat.size2(),
                          0.,&*K(rr,cc).data().begin(),K(rr,cc).size2());
            } else {
              //ublas::noalias(K(rr,cc)) += prod(trans(BTD) , col_Mat ); // int BT*D*B
              cblas_dgemm(CblasRowMajor,CblasTrans,CblasNoTrans,
                          BD.size2(),col_Mat.size2(),BD.size1(),
                          1.,&*BD.data().begin(),BD.size2(),
                          &*col_Mat.data().begin(),col_Mat.size2(),
                          1.,&*K(rr,cc).data().begin(),K(rr,cc).size2());
            }
          }
        }
      }
      PetscFunctionReturn(0);

    }
    //--------------------------------------------------------------------------------------------------------------------------------------------------//
    PetscErrorCode ComputeFibreDirection(vector<ublas::matrix<double> > &normalized_phi) {
      PetscFunctionBegin;
      ierr = OpStudentStart_TET(g_NTET); CHKERRQ(ierr);

      EntityHandle fe_handle = fePtr->get_ent();

      Range tetNodes;
      rval = moab.get_connectivity(&fe_handle,1,tetNodes); CHKERR_THROW(rval);

      vector< ublas::matrix< FieldData > > phi;
      ierr = GetGaussDiffDataVector("POTENTIAL_FIELD",phi); CHKERRQ(ierr);
      double fibreVector[3];

      for (unsigned int gg=0; gg<phi.size(); gg++) {
        normalized_phi[gg].resize(1,3);
        for (int ii=0; ii<3; ii++) {
          normalized_phi[gg](0,ii) = -phi[gg](0,ii)/sqrt(pow(phi[gg](0,0),2)+pow(phi[gg](0,1),2)+pow(phi[gg](0,2),2));
          fibreVector[ii] = normalized_phi[0](0,ii);
        }
      }

      for(Range::iterator niit1 = tetNodes.begin();niit1!=tetNodes.end();niit1++){
        rval = moab.tag_set_data(th_fibre_dir,&*niit1,1,&fibreVector[0]); CHKERR_PETSC(rval);
      }

      ierr = OpStudentEnd(); CHKERRQ(ierr);
      PetscFunctionReturn(0);
    }

  };

  //==================================================================================================================================================//

  /**
   * \brief PostProcessing Functions which computes the Strains and Elastic Linear Stresses for Transverse Isotropic Material and rotated when necessary
   *
   * Look at \ref<TranIsotropicAxisAngleRotElasticFEMethod>TranIsotropicAxisAngleRotElasticFEMethod Function for parameters
   */
  struct TranIso_PostProc_AxisAngle_OnRefMesh: public PostProcDisplacemenysAndStarinAndElasticLinearStressOnRefMesh {

    double E_p, E_z, nu_p, nu_pz, G_zp;
    int noAA;
    double *AxVector, *AxAngle;
    Tag th_fibre_orientation;

    TranIso_PostProc_AxisAngle_OnRefMesh( FieldInterface& _mField,double _lambda,double _mu, double _E_p,double _E_z, double _nu_p, double _nu_pz, double _G_zp, int _noAA, double *_AxVector, double *_AxAngle):
    PostProcDisplacemenysAndStarinAndElasticLinearStressOnRefMesh(_mField,"DISPLACEMENT",_lambda,_mu),E_p(_E_p),E_z(_E_z),nu_p(_nu_p),nu_pz(_nu_pz),G_zp(_G_zp), noAA(_noAA), AxVector(_AxVector), AxAngle(_AxAngle) {

      double def_VAL2[3] = {0,0,0};
      rval = moab_post_proc.tag_get_handle("FIBRE_DIRECTION",3,MB_TYPE_DOUBLE,th_fibre_orientation,MB_TAG_CREAT|MB_TAG_SPARSE,&def_VAL2); CHKERR_THROW(rval);

    };

    PetscErrorCode operator()() {
      PetscFunctionBegin;

      ierr = do_operator(); CHKERRQ(ierr);
      ierr = OpStudentStart_TET(g_NTET); CHKERRQ(ierr);

      ///Get Stiffness Matrix
      ublas::symmetric_matrix<FieldData,ublas::upper> StiffnessMatrix;
      StiffnessMatrix.resize(6);
      StiffnessMatrix.clear();
      TransverseIsotropicStiffnessMatrix TranIsoMat(nu_p, nu_pz, E_p, E_z, G_zp);
      StiffnessMatrix=TranIsoMat.StiffnessMatrix;
      //        IsotropicStiffnessMatrix IsoMat(lambda, mu);
      //        StiffnessMatrix=IsoMat.StiffnessMatrix;

      ///Rotating the Stiffness matrix according a set of axes of rotations and their respective angle

      int noOfRotations = noAA; //Number of Rotations
      double negAxAngle[noOfRotations];
      for (int aa=0; aa<noOfRotations; aa++) negAxAngle[aa]=-AxAngle[aa];

      ublas::matrix<double> DummyMatrix,DummyMatrix2;
      DummyMatrix = ublas::zero_matrix<FieldData>(6,6);
      DummyMatrix = StiffnessMatrix;

      ///Rotating Stiffness over a number of axis/angle rotations
      for (int aa=0; aa<noOfRotations; aa++) {

        StressTransformation StressRotMat(&AxVector[3*aa], AxAngle[aa]);
        StrainTransformation invStrainRotMat(&AxVector[3*aa], negAxAngle[aa]);

        ublas::matrix<double> TrpMatrixStress;
        TrpMatrixStress = ublas::zero_matrix<FieldData>(6,6);
        TrpMatrixStress=StressRotMat.StressRotMat;

        ublas::matrix<double> TrpMatrixInvStrain;
        TrpMatrixInvStrain = ublas::zero_matrix<FieldData>(6,6);
        TrpMatrixInvStrain=invStrainRotMat.StrainRotMat;

        DummyMatrix2 = ublas::zero_matrix<FieldData>(6,6);
        ublas::matrix< FieldData > dummyA = prod( DummyMatrix , TrpMatrixInvStrain );
        DummyMatrix2 = prod(TrpMatrixStress,dummyA);
        DummyMatrix = ublas::zero_matrix<FieldData>(6,6);
        DummyMatrix = DummyMatrix2;
      }

      D.resize(6,6);
      D.clear();
      D = DummyMatrix;

      int gg=0;
      vector< ublas::matrix< FieldData > > GradU_at_GaussPt;
      ierr = GetGaussDiffDataVector(field_name,GradU_at_GaussPt); CHKERRQ(ierr);
      vector< ublas::matrix< FieldData > >::iterator viit = GradU_at_GaussPt.begin();
      map<EntityHandle,EntityHandle>::iterator mit = node_map.begin();

      for(;viit!=GradU_at_GaussPt.end();viit++,mit++,gg++) {

        ///Compute Strains and save them on TAG
        ublas::matrix< FieldData > GradU = *viit;
        ublas::matrix< FieldData > Strain = 0.5*( GradU + trans(GradU) );
        rval = moab_post_proc.tag_set_data(th_strain,&mit->second,1,&(Strain.data()[0])); CHKERR_PETSC(rval);

        ublas::matrix<double> AARotMatrix;
        AARotMatrix = ublas::identity_matrix<FieldData>(3);

        for (int aa=0; aa<noOfRotations; aa++) {

          AxisAngleRotationalMatrix RotMatrix(&AxVector[3*aa], AxAngle[aa]);

          ublas::matrix<double> rotationalMat;
          rotationalMat = ublas::zero_matrix<FieldData>(3,3);
          rotationalMat=RotMatrix.AARotMat;

          ublas::matrix<double> AARotMatrix1;
          AARotMatrix1 = prod(rotationalMat,AARotMatrix);
          AARotMatrix = ublas::zero_matrix<FieldData>(3,3);
          AARotMatrix = AARotMatrix1;
        }

        ///Rotate AxisYVector[0,1,0] to the direction of the fibre and save in TAG
        ublas::vector<FieldData> AxisYVector(3);
        AxisYVector[0]=0; AxisYVector[1]=0;AxisYVector[2]=1;
        ublas::vector<FieldData> Fibre = prod(AARotMatrix,AxisYVector);

        rval = moab_post_proc.tag_set_data(th_fibre_orientation,&mit->second,1,&Fibre[0]); CHKERR_PETSC(rval);

        ///calculate stress and save it into tag
        ublas::vector<FieldData> Strain_VectorNotation(6);
        Strain_VectorNotation[0] = Strain(0,0);
        Strain_VectorNotation[1] = Strain(1,1);
        Strain_VectorNotation[2] = Strain(2,2);
        Strain_VectorNotation[3] = 2*Strain(0,1);
        Strain_VectorNotation[4] = 2*Strain(1,2);
        Strain_VectorNotation[5] = 2*Strain(2,0);
        ublas::vector< FieldData > Stress_VectorNotation = prod( D, Strain_VectorNotation );
        ublas::matrix< FieldData > Stress = ublas::zero_matrix<FieldData>(3,3);
        Stress(0,0) = Stress_VectorNotation[0];
        Stress(1,1) = Stress_VectorNotation[1];
        Stress(2,2) = Stress_VectorNotation[2];
        Stress(0,1) = Stress(1,0) = Stress_VectorNotation[3];
        Stress(1,2) = Stress(2,1) = Stress_VectorNotation[4];
        Stress(2,0) = Stress(0,2) = Stress_VectorNotation[5];

        rval = moab_post_proc.tag_set_data(th_stress,&mit->second,1,&(Stress.data()[0])); CHKERR_PETSC(rval);
      }

      ierr = OpStudentEnd(); CHKERRQ(ierr);
      PetscFunctionReturn(0);
    }

  };

  /**
   * \brief PostProcessing Functions which computes the Strains and Elastic Linear Stresses for rotated Transverse Isotropic Material <br>
   * Using direction computed in potential flow problem
   *
   * Look at \ref<TranIsotropicAxisAngleRotElasticFEMethod>TranIsotropicAxisAngleRotElasticFEMethod Function for parameters
   */
  struct TranIso_PostProc_FibreDirRot_OnRefMesh: public PostProcDisplacemenysAndStarinAndElasticLinearStressOnRefMesh {

    double E_p, E_z, nu_p, nu_pz, G_zp;
    Tag th_fibre_orientation;

    TranIso_PostProc_FibreDirRot_OnRefMesh( FieldInterface& _mField,double _lambda,double _mu, double _E_p,double _E_z, double _nu_p, double _nu_pz, double _G_zp):
    PostProcDisplacemenysAndStarinAndElasticLinearStressOnRefMesh(_mField,"DISPLACEMENT",_lambda,_mu),E_p(_E_p),E_z(_E_z),nu_p(_nu_p),nu_pz(_nu_pz),G_zp(_G_zp) {

      double def_VAL2[3] = {0,0,0};
      rval = moab_post_proc.tag_get_handle("FIBRE_DIRECTION",3,MB_TYPE_DOUBLE,th_fibre_orientation,MB_TAG_CREAT|MB_TAG_SPARSE,&def_VAL2); CHKERR_THROW(rval);

    };

    PetscErrorCode ComputeGradient(vector<ublas::matrix<double> > &normalized_phi) {
      PetscFunctionBegin;
      ierr = OpStudentStart_TET(g_NTET); CHKERRQ(ierr);

      EntityHandle fe_handle = fePtr->get_ent();

      Range tetNodes;
      rval = moab.get_connectivity(&fe_handle,1,tetNodes); CHKERR_THROW(rval);

      vector< ublas::matrix< FieldData > > phi;
      ierr = GetGaussDiffDataVector("POTENTIAL_FIELD",phi); CHKERRQ(ierr);

      for (unsigned int gg=0; gg<phi.size(); gg++) {
        normalized_phi[gg].resize(1,3);
        for (int ii=0; ii<3; ii++) {
          normalized_phi[gg](0,ii) = -phi[gg](0,ii)/sqrt(pow(phi[gg](0,0),2)+pow(phi[gg](0,1),2)+pow(phi[gg](0,2),2));
        }
      }

      ierr = OpStudentEnd(); CHKERRQ(ierr);
      PetscFunctionReturn(0);
    }


    PetscErrorCode operator()() {
      PetscFunctionBegin;

      ierr = do_operator(); CHKERRQ(ierr);
      ierr = OpStudentStart_TET(g_NTET); CHKERRQ(ierr);

      ///Get Stiffness Matrix
      ublas::symmetric_matrix<FieldData,ublas::upper> StiffnessMatrix;
      StiffnessMatrix.resize(6);
      StiffnessMatrix.clear();
      TransverseIsotropicStiffnessMatrix TranIsoMat(nu_p, nu_pz, E_p, E_z, G_zp);
      StiffnessMatrix=TranIsoMat.StiffnessMatrix;
      //        IsotropicStiffnessMatrix IsoMat(lambda, mu);
      //        StiffnessMatrix=IsoMat.StiffnessMatrix;

      int gg=0;
      vector< ublas::matrix< FieldData > > GradU_at_GaussPt;
      ierr = GetGaussDiffDataVector(field_name,GradU_at_GaussPt); CHKERRQ(ierr);
      vector< ublas::matrix< FieldData > >::iterator viit = GradU_at_GaussPt.begin();
      map<EntityHandle,EntityHandle>::iterator mit = node_map.begin();

      vector< ublas::matrix< double > > normalized_phi;
      normalized_phi.resize(GradU_at_GaussPt.size());
      ierr = ComputeGradient(normalized_phi); CHKERRQ(ierr);

      cout<<GradU_at_GaussPt.size()<<endl;

      for(;viit!=GradU_at_GaussPt.end();viit++,mit++,gg++) {

        double zVec[3]={0.0,0.0,1.0};

        double AxVector[3]={normalized_phi[gg](0,1)*zVec[2]-normalized_phi[gg](0,2)*zVec[1] , normalized_phi[gg](0,2)*zVec[0]-normalized_phi[gg](0,0)*zVec[2] , normalized_phi[gg](0,0)*zVec[1]-normalized_phi[gg](0,1)*zVec[0]};
        double AxAngle[1]= {asin((sqrt(pow(AxVector[0],2)+pow(AxVector[1],2)+pow(AxVector[2],2)))/(sqrt(pow(normalized_phi[gg](0,0),2)+pow(normalized_phi[gg](0,1),2)+pow(normalized_phi[gg](0,2),2)))*(sqrt(pow(zVec[0],2)+pow(zVec[1],2)+pow(zVec[2],2))))};

        ///Rotating the Stiffness matrix according a set of axes of rotations and their respective angle

        int noOfRotations = 1; //Number of Rotations
        double negAxAngle[noOfRotations];
        for (int aa=0; aa<noOfRotations; aa++) negAxAngle[aa]=-AxAngle[aa];

        ublas::matrix<double> DummyMatrix,DummyMatrix2;
        DummyMatrix = ublas::zero_matrix<FieldData>(6,6);
        DummyMatrix = StiffnessMatrix;

        ///Rotating Stiffness over a number of axis/angle rotations
        for (int aa=0; aa<noOfRotations; aa++) {

          StressTransformation StressRotMat(&AxVector[3*aa], AxAngle[aa]);
          StrainTransformation invStrainRotMat(&AxVector[3*aa], negAxAngle[aa]);

          ublas::matrix<double> TrpMatrixStress;
          TrpMatrixStress = ublas::zero_matrix<FieldData>(6,6);
          TrpMatrixStress=StressRotMat.StressRotMat;

          ublas::matrix<double> TrpMatrixInvStrain;
          TrpMatrixInvStrain = ublas::zero_matrix<FieldData>(6,6);
          TrpMatrixInvStrain=invStrainRotMat.StrainRotMat;

          DummyMatrix2 = ublas::zero_matrix<FieldData>(6,6);
          ublas::matrix< FieldData > dummyA = prod( DummyMatrix , TrpMatrixInvStrain );
          DummyMatrix2 = prod(TrpMatrixStress,dummyA);
          DummyMatrix = ublas::zero_matrix<FieldData>(6,6);
          DummyMatrix = DummyMatrix2;
        }

        D.resize(6,6);
        D.clear();
        D = DummyMatrix;

        ///Compute Strains and save them on TAG
        ublas::matrix< FieldData > GradU = *viit;
        ublas::matrix< FieldData > Strain = 0.5*( GradU + trans(GradU) );
        rval = moab_post_proc.tag_set_data(th_strain,&mit->second,1,&(Strain.data()[0])); CHKERR_PETSC(rval);

        ublas::matrix<double> AARotMatrix;
        AARotMatrix = ublas::identity_matrix<FieldData>(3);

        for (int aa=0; aa<noOfRotations; aa++) {

          AxisAngleRotationalMatrix RotMatrix(&AxVector[3*aa], AxAngle[aa]);

          ublas::matrix<double> rotationalMat;
          rotationalMat = ublas::zero_matrix<FieldData>(3,3);
          rotationalMat=RotMatrix.AARotMat;

          ublas::matrix<double> AARotMatrix1;
          AARotMatrix1 = prod(rotationalMat,AARotMatrix);
          AARotMatrix = ublas::zero_matrix<FieldData>(3,3);
          AARotMatrix = AARotMatrix1;
        }

        ///Rotate AxisYVector[0,1,0] to the direction of the fibre and save in TAG
        ublas::vector<FieldData> AxisYVector(3);
        AxisYVector[0]=0; AxisYVector[1]=0;AxisYVector[2]=1;
        ublas::vector<FieldData> Fibre = prod(AARotMatrix,AxisYVector);

        rval = moab_post_proc.tag_set_data(th_fibre_orientation,&mit->second,1,&Fibre[0]); CHKERR_PETSC(rval);

        ///calculate stress and save it into tag
        ublas::vector<FieldData> Strain_VectorNotation(6);
        Strain_VectorNotation[0] = Strain(0,0);
        Strain_VectorNotation[1] = Strain(1,1);
        Strain_VectorNotation[2] = Strain(2,2);
        Strain_VectorNotation[3] = 2*Strain(0,1);
        Strain_VectorNotation[4] = 2*Strain(1,2);
        Strain_VectorNotation[5] = 2*Strain(2,0);
        ublas::vector< FieldData > Stress_VectorNotation = prod( D, Strain_VectorNotation );
        ublas::matrix< FieldData > Stress = ublas::zero_matrix<FieldData>(3,3);
        Stress(0,0) = Stress_VectorNotation[0];
        Stress(1,1) = Stress_VectorNotation[1];
        Stress(2,2) = Stress_VectorNotation[2];
        Stress(0,1) = Stress(1,0) = Stress_VectorNotation[3];
        Stress(1,2) = Stress(2,1) = Stress_VectorNotation[4];
        Stress(2,0) = Stress(0,2) = Stress_VectorNotation[5];

        rval = moab_post_proc.tag_set_data(th_stress,&mit->second,1,&(Stress.data()[0])); CHKERR_PETSC(rval);
      }

      ierr = OpStudentEnd(); CHKERRQ(ierr);
      PetscFunctionReturn(0);
    }

  };

  //============================================================================
  //
  // Solving elastic FE problem for textile RVE consisting of two material
  // phases matrix and fibre. To use meso-scale FE, yarn comprising of fibre
  // bundle and immersed matrix is considered as a single phase material, and
  // its material or mechanical properties are predicted by homogenization
  // method. Here we use the Mori-Tanaka theory to approximate/predict the
  // effective elastic properties.
  //
  // Reference:
  //   [1]
  //
  //============================================================================

  /*
   * Function to calculate stiffness matrix of yarn consisting of fibre and
   * matrix using Mori-Tanaka method (asymptotic homogenization method) to
   * predict the effective elastic properties
   *
   * parameters:
   *
   *         (1) transversely isotropic material
   *       nu_p: transverse Poisson's ratio
   *      nu_pz: axial Poisson's ratio
   *        E_p: transverse elastic modulus
   *        E_z: axial/longitudinal elastic modulus
   *       G_zp: shear modulus
   *
   *         (2) isotropic material
   *     lambda: lame parameter
   *         mu: shear modulus
   *
   *         (3) geomtrical parameters
   *         vf: volume fraction of fibre
   *         vm: volume fraction of matrix
   */
  struct YarnStiffnessMatrix {
    double nu_p, nu_pz, E_p, E_z, G_zp; // engineering parameters for trans-iso
    double lambda, mu;                 // lame parameters for isotropic material
    double vf;                                          // fibre volume fraction
    ublas::symmetric_matrix<FieldData,ublas::upper> StiffnessMatrix;

    YarnStiffnessMatrix(double nu_p, double nu_pz, double E_p, double E_z, double G_zp,
                        double lambda, double mu,
                        double vf){
      StiffnessMatrix.resize(6);
      StiffnessMatrix.clear();
      double k_t, l_t, n_t, m_t, p_t; // Hill's moduli for transversely isotropic material
      double k_i, l_i, n_i, m_i, p_i; // Hill's moduli for isotropic material

      // Calculate Hill's moduli for tran-iso material using its engineering parameters
      double nu_zp=(nu_pz*E_z)/E_p;
      double G_p = E_p/(2*(1+nu_p));

      k_t = 1/(2*(1-nu_p)/E_p-4*nu_zp*nu_zp/E_z);
      l_t = 2*k_t*nu_zp;
      n_t = E_z+l_t*l_t/k_t;
      m_t = G_p;
      p_t = G_zp; // need to be checked G_21 = G_12 [?]

      // Calculate Hill's moduli for isotropic material using its engineering parameters
      double K_value = lambda+2*mu/3;
      double G_value = mu;

      k_i = K_value + G_value/3;
      l_i = K_value - 2*G_value/3;
      n_i = K_value + 4*G_value/3;
      m_i = G_value;
      p_i = G_value;

      /*
       * Calculate Hill's modulus for two-phase composite with fibre of
       * transversely isotropic material and matrix of isotropic material using
       * Mori-Tanaka asymptotic method based on Dvorak derived simple formula
       */

      double k_c, l_c, n_c, m_c, p_c;
      double k_f, l_f, n_f, m_f, p_f;            // Hill's moduli for fibre yarn
      double k_m, l_m, n_m, m_m, p_m;                // Hill's moduli for matrix
      double vm;                                       // matrix volume fraction
      k_f=k_t; l_f=l_t; n_f=n_t; m_f=m_t; p_f=p_t;
      k_m=k_i; l_m=l_i; n_m=n_i; m_m=m_i; p_m=p_i;

      //vf=0.6;
      vm = 1-vf;
      //cout<<"Volume fraction: \t"<<vf<<endl;

      k_c = (vf*k_f*(k_m+m_m)+vm*k_m*(k_f+m_m))/(vf*(k_m+m_m)+vm*(k_f+m_m));
      m_c = (m_f*m_m*(k_m+2*m_m)+k_m*m_m*(vf*m_f+vm*m_m))/(k_m*m_m+(k_m+2*m_m)*(vf*m_m+vm*m_f));
      p_c = (2*vf*p_f*p_m+vm*(p_f*p_m+p_m*p_m))/(2*vf*p_m+vm*(p_f+p_m));
      l_c = (vf*l_f*(k_m+m_m)+vm*l_m*(k_f+m_m))/(vf*(k_m+m_m)+vm*(p_f+p_m));
      n_c = vf*n_f+vm*n_m+(l_c-vf*l_f-vm*l_m)*(l_f-l_m)/(k_f-k_m);

	  // case 1: fibre direction in x-axis
      /*StiffnessMatrix(0,0) = n_c;
      StiffnessMatrix(0,1) = StiffnessMatrix(0,2) = l_c;
      StiffnessMatrix(1,1) = StiffnessMatrix(2,2) = k_c + m_c;
      StiffnessMatrix(1,2) = k_c - m_c;
      StiffnessMatrix(3,3) = m_c;
      StiffnessMatrix(4,4) = StiffnessMatrix(5,5) = p_c;*/
	  // case 2: fibre direction in z-axis
	  StiffnessMatrix(0,0) = StiffnessMatrix(1,1) = k_c + m_c;
	  StiffnessMatrix(0,1) = k_c - m_c;
	  StiffnessMatrix(0,2) = StiffnessMatrix(1,2) = l_c;
	  StiffnessMatrix(2,2) = n_c;
	  StiffnessMatrix(3,3) = m_c;
	  StiffnessMatrix(4,4) = StiffnessMatrix(5,5) = p_c;
      // cout<<"Yarn C matrix \t"<<StiffnessMatrix<<endl;
    }

  };

  struct YarnComplianceMatrix {
    double nu_p, nu_pz, E_p, E_z, G_zp; // engineering parameters for trans-iso
    double lambda, mu;                 // lame parameters for isotropic material
    double vf;                                          // fibre volume fraction
    ublas::symmetric_matrix<FieldData,ublas::upper> ComplianceMatrix;

    YarnComplianceMatrix(double nu_p, double nu_pz, double E_p, double E_z, double G_zp,
                        double lambda, double mu,
                        double vf){
      ComplianceMatrix.resize(6);
      ComplianceMatrix.clear();
      double k_t, l_t, n_t, m_t, p_t; // Hill's moduli for transversely isotropic material
      double k_i, l_i, n_i, m_i, p_i; // Hill's moduli for isotropic material

      // Calculate Hill's moduli for tran-iso material using its engineering parameters
      double nu_zp=(nu_pz*E_z)/E_p;
      double G_p = E_p/(2*(1+nu_p));

      k_t = 1/(2*(1-nu_p)/E_p-4*nu_zp*nu_zp/E_z);
      l_t = 2*k_t*nu_zp;
      n_t = E_z+l_t*l_t/k_t;
      m_t = G_p;
      p_t = G_zp; // need to be checked G_21 = G_12 [?]

      // Calculate Hill's moduli for isotropic material using its engineering parameters
      double K_value = lambda+2*mu/3;
      double G_value = mu;

      k_i = K_value + G_value/3;
      l_i = K_value - 2*G_value/3;
      n_i = K_value + 4*G_value/3;
      m_i = G_value;
      p_i = G_value;

      /*
       * Calculate Hill's modulus for two-phase composite with fibre of
       * transversely isotropic material and matrix of isotropic material using
       * Mori-Tanaka asymptotic method based on Dvorak derived simple formula
       */

      double k_c, l_c, n_c, m_c, p_c;
      double k_f, l_f, n_f, m_f, p_f;            // Hill's moduli for fibre yarn
      double k_m, l_m, n_m, m_m, p_m;                // Hill's moduli for matrix
      double vm;                                       // matrix volume fraction
      k_f=k_t; l_f=l_t; n_f=n_t; m_f=m_t; p_f=p_t;
      k_m=k_i; l_m=l_i; n_m=n_i; m_m=m_i; p_m=p_i;

      vm = 1-vf;

      k_c = (vf*k_f*(k_m+m_m)+vm*k_m*(k_f+m_m))/(vf*(k_m+m_m)+vm*(k_f+m_m));
      m_c = (m_f*m_m*(k_m+2*m_m)+k_m*m_m*(vf*m_f+vm*m_m))/(k_m*m_m+(k_m+2*m_m)*(vf*m_m+vm*m_f));
      p_c = (2*vf*p_f*p_m+vm*(p_f*p_m+p_m*p_m))/(2*vf*p_m+vm*(p_f+p_m));
      l_c = (vf*l_f*(k_m+m_m)+vm*l_m*(k_f+m_m))/(vf*(k_m+m_m)+vm*(p_f+p_m));
      n_c = vf*n_f+vm*n_m+(l_c-vf*l_f-vm*l_m)*(l_f-l_m)/(k_f-k_m);

      double E1, E2;// E3;
      double NU12, NU23; //NU13, NU23;
      double G12; //G13, G23;

      E1 = n_c - l_c*l_c/k_c;
      E2 = 4*m_c*(k_c*n_c-l_c*l_c)/(n_c*(k_c+m_c)-l_c*l_c);
      //E3 = E2;
      NU12 = l_c/k_c/2;
      //NU13 = NU12;
      NU23 = (n_c*(k_c-m_c)-l_c*l_c)/(n_c*(k_c+m_c)-l_c*l_c);
      G12 = p_c;
      //G13 = G12;
      //G23 = m_c;

	  // case 1: fibre direction in x-axis
//      ComplianceMatrix(0,0) = 1/E1;
//      ComplianceMatrix(0,1) = ComplianceMatrix(1,0) = -NU12/E1;
//      ComplianceMatrix(0,2) = ComplianceMatrix(2,0) = -NU13/E1;
//      ComplianceMatrix(1,1) = 1/E2;
//      ComplianceMatrix(1,2) = ComplianceMatrix(2,1) = -NU23/E2;
//      ComplianceMatrix(2,2) = 1/E3;
//      ComplianceMatrix(3,3) = 2*(1+NU23)/E2;
//      ComplianceMatrix(4,4) = ComplianceMatrix(5,5) = 1/G12;

	  // case 2: fibre direction in z-axis
	  ComplianceMatrix(0,0) = ComplianceMatrix(1,1) = 1/E2;
	  ComplianceMatrix(0,1) = -NU23/E2;
	  ComplianceMatrix(0,2) = ComplianceMatrix(1,2) = -NU12/E1;
	  ComplianceMatrix(2,2) = 1/E1;
      ComplianceMatrix(3,3) = 2*(1+NU23)/E2;
      ComplianceMatrix(4,4) = ComplianceMatrix(5,5) = 1/G12;
      // cout<<"Yarn S matrix \t"<<ComplianceMatrix<<endl;
    }
  };

  /*
   * Function to build up the stiffness matrix (using Mori-Tanaka theory)
   *   rotate it with respect to
   *          (1) misalignment angle of fibre
   *          (2) fibre yarn direction
   */
  struct TranIso_YarnMisalignmentElasticFEMethod: public ElasticFEMethod {

    Tag th_fibre_dir;
    bool propeties_from_BLOCKSET_MAT_ELASTICSET;

    TranIso_YarnMisalignmentElasticFEMethod(FieldInterface& _mField, Mat &_Aij, Vec _D, Vec _F, string _field_name = "DISPLACEMENT"):ElasticFEMethod(_mField, _Aij, _D, _F, 0, 0, _field_name) {


      propeties_from_BLOCKSET_MAT_ELASTICSET = false;
      for(_IT_CUBITMESHSETS_BY_BCDATA_TYPE_FOR_LOOP_(mField,BLOCKSET|MAT_ELASTICSET,it)) {
        propeties_from_BLOCKSET_MAT_ELASTICSET = true;
      }

      double def_VAL2[3] = {0,0,0};
      rval = moab.tag_get_handle("POT_FLOW_FIBRE_DIR",3,MB_TYPE_DOUBLE,
                                 th_fibre_dir,MB_TAG_CREAT|MB_TAG_SPARSE,
                                 &def_VAL2); CHKERR_THROW(rval);
      cout<<"From misalignment"<<endl;
    }

    /* =========================================================================
     *
     * Get material paramters from FE model
     *
     * =======================================================================*/
    PetscErrorCode GetMatParameters(double *_E_p, double *_E_z, double *_nu_p,
                                    double *_nu_pz, double *_G_zp,
                                    double *_lambda, double *_mu,
                                    double *_vf,
                                    double *_theta_f, double *_WavinessFactor) {
      PetscFunctionBegin;
      double YoungModulus_1, PoissonRatio_1;       // For 1st isotropic material
      // double YoungModulus_2, PoissonRatio_2;       // For 2nd isotropic material
      //double FibreVolFrac;                              // Fibre volume fraction

      if (propeties_from_BLOCKSET_MAT_ELASTICSET) {
        EntityHandle ent = fePtr->get_ent();
        for (_IT_CUBITMESHSETS_BY_BCDATA_TYPE_FOR_LOOP_(mField,BLOCKSET|MAT_ELASTICSET,it)) {
          // isotropic material 1 - matrix
          //cout<<it->get_name()<<endl;
//          if (it->get_name().compare(0,13,"MAT_ELASTIC_1") == 0) {cout<<"Isotropic"<<endl;
//            Mat_Elastic mydata;
//            ierr = it->get_attribute_data_structure(mydata); CHKERRQ(ierr);
//
//            Range meshsets;
//            rval = moab.get_entities_by_type(it->meshset,MBENTITYSET,meshsets,true); CHKERR_PETSC(rval);
//            meshsets.insert(it->meshset);
//            for (Range::iterator mit = meshsets.begin(); mit != meshsets.end(); mit++) {
//              if (moab.contains_entities(*mit,&ent,1)) {cout<<"Something"<<endl;
//                *_lambda = LAMBDA(mydata.data.Young,mydata.data.Poisson);
//                *_mu = MU(mydata.data.Young,mydata.data.Poisson);
//                //*_vf = mydata.data.User1;
//                cout<<"Fibre volume fraction \t"<<mydata.data.Poisson<<endl;
//                PetscFunctionReturn(0);
//              }
//            }
//          }
          // transversly isotropic material
         if (it->get_name().compare(0,22,"MAT_ELASTIC_TRANSISO_1") == 0) {

           //Mat_Elastic_TransIso mydata;
           Mat_Elastic mydata;
           ierr = it->get_attribute_data_structure(mydata); CHKERRQ(ierr);

           Range meshsets;
           rval = moab.get_entities_by_type(it->meshset,MBENTITYSET,meshsets,true); CHKERR_PETSC(rval);
           meshsets.insert(it->meshset);
           for (Range::iterator mit = meshsets.begin(); mit != meshsets.end(); mit++) {
             if (moab.contains_entities(*mit,&ent,1)) {
               //*_E_p   = mydata.data.Youngp;
               //*_E_z   = mydata.data.Youngz;
               //*_nu_p  = mydata.data.Poissonp;
               //*_nu_pz = mydata.data.Poissonpz;
               //*_G_zp  = mydata.data.Shearzp;

               YoungModulus_1 = mydata.data.Young;
               PoissonRatio_1 = mydata.data.Poisson;

               *_lambda = LAMBDA(YoungModulus_1,PoissonRatio_1);
               *_mu     = MU(YoungModulus_1,PoissonRatio_1);

               *_WavinessFactor = mydata.data.ThermalExpansion; // waviness factor

               *_E_p     = mydata.data.User1;
               *_E_z     = mydata.data.User2;
               *_nu_p    = mydata.data.User3;
               *_nu_pz   = mydata.data.User4;
               *_G_zp    = mydata.data.User5;

               *_vf      = mydata.data.User6; // Fibre volume fraction

               *_theta_f = mydata.data.User7; // Misalignment angle

//               cout<<"Young modulus \t = \t"<<YoungModulus_1<<endl;
//               cout<<"Poisson's ratio \t = \t"<<PoissonRatio_1<<endl;
//               cout<<"Axial modulus \t = \t"<<*_E_p<<endl;
//               cout<<"Transverse modulus \t = \t"<<*_E_z<<endl;
//               cout<<"Transverse Poisson \t = \t"<<*_nu_p<<endl;
//               cout<<"Axial Poisson \t = \t"<<*_nu_pz<<endl;
//               cout<<"Fibre volue fraction \t = \t"<<*_vf<<endl;
//               cout<<"Matrix modulus \t = \t"<<YoungModulus_1<<endl;
//               cout<<"Matrix Poisson \t = \t"<<PoissonRatio_1<<endl;

              }
            }
          }
        }
      }
      //cout<<"Get material parameters"<<endl;
      PetscFunctionReturn(0);
    }

    /* =========================================================================
     *
     * Calculate constitutent matrix D
     *
     * =======================================================================*/
    vector< ublas::matrix<FieldData> > D_At_GaussPoint;
    PetscErrorCode calculateD(double _E_p, double _E_z, double _nu_p,
                              double _nu_pz, double _G_zp,
                              double _lambda, double _mu,
                              double _vf,
                              double _theta_f, double _WavinessFactor) {
      PetscFunctionBegin;

      // ------
     // 1. Get stiffness matrix
      ublas::symmetric_matrix<FieldData,ublas::upper> StiffnessMatrix;
      StiffnessMatrix.resize(6);
      StiffnessMatrix.clear();
      YarnStiffnessMatrix TranIsoMat(_nu_p, _nu_pz, _E_p, _E_z, _G_zp,
                                     _lambda, _mu, _vf);
      StiffnessMatrix = TranIsoMat.StiffnessMatrix;

      //cout<<"Fibre volume fraction = \t"<<_vf<<endl;
      // ------
	  // 2. Rotating stiffness matrix according to misalignment angle
      //
      int noOfRotations_Mis = 1; // Number of rotation times
      double AxVector_Mis[3] = {1.0,0.0,0.0};// direction: the axis of rotation
      double AxAngle_Mis[1] = {_theta_f};    // magnitude: the angle of rotation
      double negAxAngle_Mis[noOfRotations_Mis];
      for (int aa=0; aa<noOfRotations_Mis; aa++) negAxAngle_Mis[aa] = -AxAngle_Mis[aa];

      ublas::matrix<double> DummyMatrix0, DummyMatrix3;
      DummyMatrix0 = ublas::zero_matrix<FieldData>(6,6);
      DummyMatrix0 = StiffnessMatrix;

      ublas::matrix<double> AARotMat_Mis;
      AARotMat_Mis = ublas::zero_matrix<FieldData>(3,3);
      AxisAngleRotationalMatrix RotMat_Mis(&AxVector_Mis[0], AxAngle_Mis[0]);
      AARotMat_Mis=RotMat_Mis.AARotMat;

      //cout<<AARotMat_Mis<<endl;

      for (int aa=0; aa<noOfRotations_Mis; aa++) {
        StressTransformation StressRotMat(&AxVector_Mis[3*aa], AxAngle_Mis[aa]);
        StrainTransformation invStrainRotMat(&AxVector_Mis[3*aa], negAxAngle_Mis[aa]);

        ublas::matrix<double> TrpMatrixStress;
        TrpMatrixStress = ublas::zero_matrix<FieldData>(6,6);
        TrpMatrixStress = StressRotMat.StressRotMat;

        ublas::matrix<double> TrpMatrixInvStrain;
        TrpMatrixInvStrain = ublas::zero_matrix<FieldData>(6,6);
        TrpMatrixInvStrain = invStrainRotMat.StrainRotMat;

        DummyMatrix3 = ublas::zero_matrix<FieldData>(6,6);
        ublas::matrix< FieldData > dummyAA = prod( DummyMatrix0, TrpMatrixInvStrain);
        DummyMatrix3 = prod(TrpMatrixStress, dummyAA);
        DummyMatrix0  = ublas::zero_matrix<FieldData>(6,6);
        DummyMatrix0  = DummyMatrix3;
      }
      //cout<<"Original stiffness matrix"<<endl;
      //cout<<StiffnessMatrix<<endl;

      StiffnessMatrix.resize(6);
      StiffnessMatrix.clear();
      StiffnessMatrix = DummyMatrix0;


      //cout<<"Modified stiffness matrix"<<endl;
      //cout<<StiffnessMatrix<<endl;

      // 3. Rotating the stiffness matrix according to misalignment information
      //    includeing a set of axes of rotations and their respective angle

      D_At_GaussPoint.resize(coords_at_Gauss_nodes.size());

      vector< ublas::matrix< double > > normalized_phi;
      normalized_phi.resize(coords_at_Gauss_nodes.size());
      ierr = ComputeFibreDirection(normalized_phi); CHKERRQ(ierr);

      for (unsigned int gg=0; gg<coords_at_Gauss_nodes.size(); gg++) {

        int noOfRotations = 1; // Number of rotations ???

        double zVec[3] = {0.0, 0.0, 1.0};
        double AxVector[3]={normalized_phi[gg](0,1)*zVec[2]-normalized_phi[gg](0,2)*zVec[1] , normalized_phi[gg](0,2)*zVec[0]-normalized_phi[gg](0,0)*zVec[2] , normalized_phi[gg](0,0)*zVec[1]-normalized_phi[gg](0,1)*zVec[0]};
        double AxAngle[1]= {asin((sqrt(pow(AxVector[0],2)+pow(AxVector[1],2)+pow(AxVector[2],2)))/(sqrt(pow(normalized_phi[gg](0,0),2)+pow(normalized_phi[gg](0,1),2)+pow(normalized_phi[gg](0,2),2)))*(sqrt(pow(zVec[0],2)+pow(zVec[1],2)+pow(zVec[2],2))))};
        //cout<<"Rotation vector = \t"<<AxVector[0]<<"\t"<<AxVector[1]<<"\t"<<AxVector[2]<<"\t"<<AxAngle[0]<<endl;
        //cout<<"Rontation angle = \t"<<AxAngle[0]<<endl;

        double negAxAngle[noOfRotations];
        for (int aa=0; aa<noOfRotations; aa++) negAxAngle[aa]=-AxAngle[aa];

        ublas::matrix<double> DummyMatrix, DummyMatrix2;
        DummyMatrix = ublas::zero_matrix<FieldData>(6,6);
        DummyMatrix = StiffnessMatrix;

        // Rotating stiffness matrix over a number of axis/angle rotations
        for (int aa=0; aa<noOfRotations; aa++) {

          StressTransformation StressRotMat(&AxVector[3*aa], AxAngle[aa]);
          StrainTransformation invStrainRotMat(&AxVector[3*aa], negAxAngle[aa]);

          ublas::matrix<double> TrpMatrixStress;
          TrpMatrixStress = ublas::zero_matrix<FieldData>(6,6);
          TrpMatrixStress = StressRotMat.StressRotMat;

          ublas::matrix<double> TrpMatrixInvStrain;
          TrpMatrixInvStrain = ublas::zero_matrix<FieldData>(6,6);
          TrpMatrixInvStrain = invStrainRotMat.StrainRotMat;

          DummyMatrix2 = ublas::zero_matrix<FieldData>(6,6);
          ublas::matrix< FieldData > dummyA = prod( DummyMatrix, TrpMatrixInvStrain);
          DummyMatrix2 = prod(TrpMatrixStress, dummyA);
          DummyMatrix  = ublas::zero_matrix<FieldData>(6,6);
          DummyMatrix  = DummyMatrix2;
        }

        D_At_GaussPoint[gg].resize(6,6);
        D_At_GaussPoint[gg].clear();
        D_At_GaussPoint[gg] = DummyMatrix;
      }

      PetscFunctionReturn(0);
    }

    /* =========================================================================
     *
     *
     *
     * =======================================================================*/

    PetscErrorCode Fint() {
      PetscFunctionBegin;

      try {

        //Higher order approximation of geometry
        ierr = GetHierarchicalGeometryApproximation(invH,detH); CHKERRQ(ierr);

        double _E_p, _E_z, _nu_p, _nu_pz, _G_zp;
        double _lambda, _mu;
        double _vf;
        double _theta_f;
        double _WavinessFactor;
        //cout<<"Fint"<<endl;
        ierr = GetMatParameters(&_E_p, &_E_z, &_nu_p, &_nu_pz, &_G_zp,
                                &_lambda, &_mu, &_vf, &_theta_f, &_WavinessFactor); CHKERRQ(ierr);
        ierr = calculateD(_E_p, _E_z, _nu_p, _nu_pz, _G_zp,
                          _lambda, _mu, _vf, _theta_f, _WavinessFactor); CHKERRQ(ierr);

        //Gradient at Gauss points;
        vector< ublas::matrix< FieldData > > GradU_at_GaussPt;
        ierr = GetGaussDiffDataVector(fieldName,GradU_at_GaussPt); CHKERRQ(ierr);
        unsigned int g_dim = g_NTET.size()/4;
        assert(GradU_at_GaussPt.size() == g_dim);
        NOT_USED(g_dim);
        vector< ublas::matrix< FieldData > >::iterator viit = GradU_at_GaussPt.begin();
        int gg = 0;
        for(;viit!=GradU_at_GaussPt.end();viit++,gg++) {
          try {
            ublas::matrix< FieldData > GradU = *viit;
            if(!invH.empty()) {
              //GradU =
              //[ dU/dChi1 dU/dChi2 dU/dChi3 ]
              //[ dV/dChi1 dV/dChi2 dU/dChi3 ]
              //[ dW/dChi1 dW/dChi2 dW/dChi3 ]
              //H =
              //[ dX1/dChi1 dX1/dChi2 dX1/dChi3 ]
              //[ dX2/dChi1 dX2/dChi2 dX2/dChi3 ]
              //[ dX3/dChi1 dX3/dChi2 dX3/dChi3 ]
              //invH =
              //[ dChi1/dX1 dChi1/dX2 dChi1/dX3 ]
              //[ dChi2/dX1 dChi2/dX2 dChi2/dX3 ]
              //[ dChi3/dX1 dChi3/dX2 dChi3/dX3 ]
              //GradU =
              //[ dU/dX1 dU/dX2 dU/dX3 ]
              //[ dV/dX1 dV/dX2 dV/dX3 ] = GradU * invH
              //[ dW/dX1 dW/dX2 dW/dX3 ]
              GradU = prod( GradU, invH[gg] );
            }
            ublas::matrix< FieldData > Strain = 0.5*( GradU + trans(GradU) );
            ublas::vector< FieldData > VoightStrain(6);
            VoightStrain[0] = Strain(0,0);
            VoightStrain[1] = Strain(1,1);
            VoightStrain[2] = Strain(2,2);
            VoightStrain[3] = 2*Strain(0,1);
            VoightStrain[4] = 2*Strain(1,2);
            VoightStrain[5] = 2*Strain(2,0);
            double w = V*G_TET_W[gg];
            ublas::vector<FieldData> VoightStress = prod(w*D_At_GaussPoint[gg],VoightStrain);
            //BT * VoigtStress
            f_int.resize(row_mat);
            for(int rr = 0;rr<row_mat;rr++) {
              if(RowGlob[rr].size()==0) continue;
              ublas::matrix<FieldData> &B = (rowBMatrices[rr])[gg];
              if(gg == 0) {
                f_int[rr] = prod( trans(B), VoightStress );
              } else {
                f_int[rr] += prod( trans(B), VoightStress );
              }
            }
          } catch (const std::exception& ex) {
            ostringstream ss;
            ss << "thorw in method: " << ex.what() << " at line " << __LINE__ << " in file " << __FILE__ << endl;
            SETERRQ(PETSC_COMM_SELF,1,ss.str().c_str());
          }
        }

      } catch (const std::exception& ex) {
        ostringstream ss;
        ss << "thorw in method: " << ex.what() << " at line " << __LINE__ << " in file " << __FILE__ << endl;
        SETERRQ(PETSC_COMM_SELF,1,ss.str().c_str());
      }

      PetscFunctionReturn(0);
    }

    /* =========================================================================
     *
     * Element stiffness matrix
     *
     * =======================================================================*/
    PetscErrorCode Stiffness() {
      PetscFunctionBegin;

      double _E_p, _E_z, _nu_p, _nu_pz, _G_zp;
      double _lambda, _mu;
      double _vf;
      double _theta_f;
      double _WavinessFactor;
      //cout<<"Stiffness"<<endl;
      ierr = GetMatParameters(&_E_p, &_E_z, &_nu_p, &_nu_pz, &_G_zp,
                              &_lambda, &_mu, &_vf, &_theta_f, &_WavinessFactor); CHKERRQ(ierr);
      ierr = calculateD(_E_p, _E_z, _nu_p, _nu_pz, _G_zp,
                        _lambda, _mu, _vf, _theta_f, _WavinessFactor); CHKERRQ(ierr);
      //cout<<"K calculation"<<endl;
      K.resize(row_mat,col_mat);
      int g_dim = g_NTET.size()/4;
      for (int rr = 0; rr<row_mat; rr++) {
        if (RowGlob[rr].size()==0) continue;
        for (int gg=0; gg<g_dim; gg++) {

          //cout<<"gg = \t"<<gg<<endl;
          //cout<<"D matrix"<<endl;cout<<D_At_GaussPoint[gg]<<endl;
          ublas::matrix<FieldData> &row_Mat = (rowBMatrices[rr])[gg];
          double w = V*G_TET_W[gg];
          if (detH.size()>0) {
            w *= detH[gg];
          }
          BD.resize(6,row_Mat.size2());
          //cout<<"BD_Resize"<<endl;
          cblas_dsymm(CblasRowMajor,CblasLeft,CblasUpper,
                      BD.size1(),BD.size2(),
                      w,&*D_At_GaussPoint[gg].data().begin(),D_At_GaussPoint[gg].size2(),
                      &*row_Mat.data().begin(),row_Mat.size2(),
                      0.,&*BD.data().begin(),BD.size2());
          for(int cc = rr;cc<col_mat;cc++) {
            if(ColGlob[cc].size()==0) continue;
            ublas::matrix<FieldData> &col_Mat = (colBMatrices[cc])[gg];
            if(gg == 0) {
              K(rr,cc).resize(BD.size2(),col_Mat.size2());
              //ublas::noalias(K(rr,cc)) = prod(trans(BD) , col_Mat ); // int BT*D*B
              cblas_dgemm(CblasRowMajor,CblasTrans,CblasNoTrans,
                          BD.size2(),col_Mat.size2(),BD.size1(),
                          1.,&*BD.data().begin(),BD.size2(),
                          &*col_Mat.data().begin(),col_Mat.size2(),
                          0.,&*K(rr,cc).data().begin(),K(rr,cc).size2());
            } else {
              //ublas::noalias(K(rr,cc)) += prod(trans(BTD) , col_Mat ); // int BT*D*B
              cblas_dgemm(CblasRowMajor,CblasTrans,CblasNoTrans,
                          BD.size2(),col_Mat.size2(),BD.size1(),
                          1.,&*BD.data().begin(),BD.size2(),
                          &*col_Mat.data().begin(),col_Mat.size2(),
                          1.,&*K(rr,cc).data().begin(),K(rr,cc).size2());
            }
          }
        }
      }

      //cout<<"K calculation finish"<<endl;

      PetscFunctionReturn(0);
    }

    /* =========================================================================
     *
     * Compute fibre direction
     *
     * =======================================================================*/

    PetscErrorCode ComputeFibreDirection(vector<ublas::matrix<double> > &normalized_phi) {
      PetscFunctionBegin;
      ierr = OpStudentStart_TET(g_NTET); CHKERRQ(ierr);

      EntityHandle fe_handle = fePtr->get_ent();

      Range tetNodes;
      rval = moab.get_connectivity(&fe_handle,1,tetNodes); CHKERR_THROW(rval);

      vector< ublas::matrix< FieldData > > phi;
      ierr = GetGaussDiffDataVector("POTENTIAL_FIELD",phi); CHKERRQ(ierr);
      double fibreVector[3];

      for (unsigned int gg=0; gg<phi.size(); gg++) {
        normalized_phi[gg].resize(1,3);
        for (int ii=0; ii<3; ii++) {
          normalized_phi[gg](0,ii) = -phi[gg](0,ii)/sqrt(pow(phi[gg](0,0),2)+pow(phi[gg](0,1),2)+pow(phi[gg](0,2),2));
          fibreVector[ii] = normalized_phi[0](0,ii);
        }
      }

      for(Range::iterator niit1 = tetNodes.begin();niit1!=tetNodes.end();niit1++){
        rval = moab.tag_set_data(th_fibre_dir,&*niit1,1,&fibreVector[0]); CHKERR_PETSC(rval);
      }

      ierr = OpStudentEnd(); CHKERRQ(ierr);
      PetscFunctionReturn(0);
    }

   };

  struct TranIso_PostProc_Misalignment_OnRefMesh: public PostProcDisplacemenysAndStarinAndElasticLinearStressOnRefMesh {
    double E_p, E_z, nu_p, nu_pz, G_zp;
    double lambda, mu;
    double vf;
    Tag th_fibre_orientation;

    TranIso_PostProc_Misalignment_OnRefMesh (FieldInterface& _mField,
                                             double _lambda, double _mu,
                                             double _E_p, double _E_z,
                                             double _nu_p, double _nu_pz,
                                             double _G_zp):
    PostProcDisplacemenysAndStarinAndElasticLinearStressOnRefMesh (_mField, "DISPLACEMENT",
                                                                   _lambda, _mu),
    E_p(_E_p),E_z(_E_z),nu_p(_nu_p),nu_pz(_nu_pz),G_zp(_G_zp) {

      double def_VAL2[3] = {0,0,0};
      rval = moab_post_proc.tag_get_handle("FIBRE_DIRECTION",3,MB_TYPE_DOUBLE,
                                           th_fibre_orientation,MB_TAG_CREAT|MB_TAG_SPARSE,
                                           &def_VAL2); CHKERR_THROW(rval);
    }

    PetscErrorCode ComputeGradient(vector<ublas::matrix<double> > &normalized_phi) {

      PetscFunctionBegin;
      ierr = OpStudentStart_TET(g_NTET); CHKERRQ(ierr);

      EntityHandle fe_handle = fePtr->get_ent();

      Range tetNodes;
      rval = moab.get_connectivity(&fe_handle,1,tetNodes); CHKERR_THROW(rval);

      vector< ublas::matrix< FieldData> > phi;
      ierr = GetGaussDiffDataVector("POTENTIAL_FIELD",phi); CHKERRQ(ierr);

      for (unsigned int gg=0; gg<phi.size(); gg++) {
        normalized_phi[gg].resize(1,3);
        for (int ii=0; ii<3; ii++) {
          normalized_phi[gg](0,ii) = -phi[gg](0,ii)/sqrt(pow(phi[gg](0,0),2)+
                                                         pow(phi[gg](0,1),2)+
                                                         pow(phi[gg](0,2),2));
        }
      }

      ierr = OpStudentEnd(); CHKERRQ(ierr);

      PetscFunctionReturn(0);
    }

    PetscErrorCode operator()() {
      PetscFunctionBegin;

      ierr = do_operator(); CHKERRQ(ierr);
      ierr = OpStudentStart_TET(g_NTET); CHKERRQ(ierr);

      // Get stifffness matrix
      ublas::symmetric_matrix<FieldData,ublas::upper> StiffnessMatrix;
      StiffnessMatrix.resize(6);
      StiffnessMatrix.clear();
      YarnStiffnessMatrix TranIsoMat(nu_p, nu_pz, E_p, E_z, G_zp, lambda, mu, vf);
      StiffnessMatrix = TranIsoMat.StiffnessMatrix;

      int gg=0;
      vector< ublas::matrix< FieldData > > GradU_at_GaussPt;
      ierr = GetGaussDiffDataVector(field_name,GradU_at_GaussPt); CHKERRQ(ierr);
      vector< ublas::matrix< FieldData > >::iterator viit = GradU_at_GaussPt.begin();
      map<EntityHandle, EntityHandle>::iterator mit = node_map.begin();

      vector< ublas::matrix< double > > normalized_phi;
      normalized_phi.resize(GradU_at_GaussPt.size());
      ierr = ComputeGradient(normalized_phi); CHKERRQ(ierr);

      cout<<GradU_at_GaussPt.size()<<endl;

      for (; viit!=GradU_at_GaussPt.end(); viit++, mit++, gg++) {

        double zVec[3]={0.0,0.0,1.0};

        double AxVector[3]={normalized_phi[gg](0,1)*zVec[2]-normalized_phi[gg](0,2)*zVec[1] , normalized_phi[gg](0,2)*zVec[0]-normalized_phi[gg](0,0)*zVec[2] , normalized_phi[gg](0,0)*zVec[1]-normalized_phi[gg](0,1)*zVec[0]};
        double AxAngle[1]= {asin((sqrt(pow(AxVector[0],2)+pow(AxVector[1],2)+pow(AxVector[2],2)))/(sqrt(pow(normalized_phi[gg](0,0),2)+pow(normalized_phi[gg](0,1),2)+pow(normalized_phi[gg](0,2),2)))*(sqrt(pow(zVec[0],2)+pow(zVec[1],2)+pow(zVec[2],2))))};

        ///Rotating the Stiffness matrix according a set of axes of rotations and their respective angle

        int noOfRotations = 1; //Number of Rotations
        double negAxAngle[noOfRotations];
        for (int aa=0; aa<noOfRotations; aa++) negAxAngle[aa]=-AxAngle[aa];

        ublas::matrix<double> DummyMatrix,DummyMatrix2;
        DummyMatrix = ublas::zero_matrix<FieldData>(6,6);
        DummyMatrix = StiffnessMatrix;

        ///Rotating Stiffness over a number of axis/angle rotations
        for (int aa=0; aa<noOfRotations; aa++) {

          StressTransformation StressRotMat(&AxVector[3*aa], AxAngle[aa]);
          StrainTransformation invStrainRotMat(&AxVector[3*aa], negAxAngle[aa]);

          ublas::matrix<double> TrpMatrixStress;
          TrpMatrixStress = ublas::zero_matrix<FieldData>(6,6);
          TrpMatrixStress=StressRotMat.StressRotMat;

          ublas::matrix<double> TrpMatrixInvStrain;
          TrpMatrixInvStrain = ublas::zero_matrix<FieldData>(6,6);
          TrpMatrixInvStrain=invStrainRotMat.StrainRotMat;

          DummyMatrix2 = ublas::zero_matrix<FieldData>(6,6);
          ublas::matrix< FieldData > dummyA = prod( DummyMatrix , TrpMatrixInvStrain );
          DummyMatrix2 = prod(TrpMatrixStress,dummyA);
          DummyMatrix = ublas::zero_matrix<FieldData>(6,6);
          DummyMatrix = DummyMatrix2;
        }

        D.resize(6,6);
        D.clear();
        D = DummyMatrix;

        ///Compute Strains and save them on TAG
        ublas::matrix< FieldData > GradU = *viit;
        ublas::matrix< FieldData > Strain = 0.5*( GradU + trans(GradU) );
        rval = moab_post_proc.tag_set_data(th_strain,&mit->second,1,&(Strain.data()[0])); CHKERR_PETSC(rval);

        ublas::matrix<double> AARotMatrix;
        AARotMatrix = ublas::identity_matrix<FieldData>(3);

        for (int aa=0; aa<noOfRotations; aa++) {

          AxisAngleRotationalMatrix RotMatrix(&AxVector[3*aa], AxAngle[aa]);

          ublas::matrix<double> rotationalMat;
          rotationalMat = ublas::zero_matrix<FieldData>(3,3);
          rotationalMat=RotMatrix.AARotMat;

          ublas::matrix<double> AARotMatrix1;
          AARotMatrix1 = prod(rotationalMat,AARotMatrix);
          AARotMatrix = ublas::zero_matrix<FieldData>(3,3);
          AARotMatrix = AARotMatrix1;
        }

        ///Rotate AxisYVector[0,1,0] to the direction of the fibre and save in TAG
        ublas::vector<FieldData> AxisYVector(3);
        AxisYVector[0]=0; AxisYVector[1]=0;AxisYVector[2]=1;
        ublas::vector<FieldData> Fibre = prod(AARotMatrix,AxisYVector);

        rval = moab_post_proc.tag_set_data(th_fibre_orientation,&mit->second,1,&Fibre[0]); CHKERR_PETSC(rval);

        ///calculate stress and save it into tag
        ublas::vector<FieldData> Strain_VectorNotation(6);
        Strain_VectorNotation[0] = Strain(0,0);
        Strain_VectorNotation[1] = Strain(1,1);
        Strain_VectorNotation[2] = Strain(2,2);
        Strain_VectorNotation[3] = 2*Strain(0,1);
        Strain_VectorNotation[4] = 2*Strain(1,2);
        Strain_VectorNotation[5] = 2*Strain(2,0);
        ublas::vector< FieldData > Stress_VectorNotation = prod( D, Strain_VectorNotation );
        ublas::matrix< FieldData > Stress = ublas::zero_matrix<FieldData>(3,3);
        Stress(0,0) = Stress_VectorNotation[0];
        Stress(1,1) = Stress_VectorNotation[1];
        Stress(2,2) = Stress_VectorNotation[2];
        Stress(0,1) = Stress(1,0) = Stress_VectorNotation[3];
        Stress(1,2) = Stress(2,1) = Stress_VectorNotation[4];
        Stress(2,0) = Stress(0,2) = Stress_VectorNotation[5];

        rval = moab_post_proc.tag_set_data(th_stress,&mit->second,1,&(Stress.data()[0])); CHKERR_PETSC(rval);
      }

      ierr = OpStudentEnd(); CHKERRQ(ierr);
      PetscFunctionReturn(0);
    }
  };

  /***************************************************************************
   *
   * Function to build up the stiffness matrix (using Mori-Tanaka theory) with
   *   consideration of fibre waviness in yarn
   *
   **************************************************************************/
  struct TranIso_FibreWavinessElasticFEMethod: public ElasticFEMethod {

    Tag th_fibre_dir;
    bool propeties_from_BLOCKSET_MAT_ELASTICSET;

    TranIso_FibreWavinessElasticFEMethod(FieldInterface& _mField, Mat &_Aij, Vec _D,
                                         Vec _F, string _field_name = "DISPLACEMENT"):
										 ElasticFEMethod(_mField, _Aij, _D, _F, 0, 0, _field_name) {


      propeties_from_BLOCKSET_MAT_ELASTICSET = false;
      for(_IT_CUBITMESHSETS_BY_BCDATA_TYPE_FOR_LOOP_(mField,BLOCKSET|MAT_ELASTICSET,it)) {
        propeties_from_BLOCKSET_MAT_ELASTICSET = true;
      }

      double def_VAL2[3] = {0,0,0};
      rval = moab.tag_get_handle("POT_FLOW_FIBRE_DIR", 3, MB_TYPE_DOUBLE,
                                 th_fibre_dir, MB_TAG_CREAT|MB_TAG_SPARSE,
                                 &def_VAL2); CHKERR_THROW(rval);
      cout<<"From fibre waviness"<<endl;
    }

    /* =========================================================================
     *
     * Get material paramters from FE model
     *
     * =======================================================================*/
    PetscErrorCode GetMatParameters(double *_E_p, double *_E_z, double *_nu_p,
                                    double *_nu_pz, double *_G_zp,
                                    double *_lambda, double *_mu,
                                    double *_vf,
                                    double *_theta_f, double *_WavinessFactor) {
      PetscFunctionBegin;
      double YoungModulus_1, PoissonRatio_1;       // For 1st isotropic material
      // double YoungModulus_2, PoissonRatio_2;       // For 2nd isotropic material
      //double FibreVolFrac;                              // Fibre volume fraction

      if (propeties_from_BLOCKSET_MAT_ELASTICSET) {
        EntityHandle ent = fePtr->get_ent();
        for (_IT_CUBITMESHSETS_BY_BCDATA_TYPE_FOR_LOOP_(mField,BLOCKSET|MAT_ELASTICSET,it)) {
		  string name = it->get_name();
		   // if (name.compare(0,13,"MAT_ELASTIC_1") == 0) {
		   //  Mat_Elastic mydata;
		   // ierr = it->get_attribute_data_structure(mydata); CHKERRQ(ierr);
		   // *_lambda = LAMBDA(mydata.data.Young,mydata.data.Poisson);
		   // *_mu = MU(mydata.data.Young,mydata.data.Poisson);
		   // }
          // isotropic material 1 - matrix
          //cout<<it->get_name()<<endl;
          //          if (it->get_name().compare(0,13,"MAT_ELASTIC_1") == 0) {cout<<"Isotropic"<<endl;
          //            Mat_Elastic mydata;
          //            ierr = it->get_attribute_data_structure(mydata); CHKERRQ(ierr);
          //
          //            Range meshsets;
          //            rval = moab.get_entities_by_type(it->meshset,MBENTITYSET,meshsets,true); CHKERR_PETSC(rval);
          //            meshsets.insert(it->meshset);
          //            for (Range::iterator mit = meshsets.begin(); mit != meshsets.end(); mit++) {
          //              if (moab.contains_entities(*mit,&ent,1)) {cout<<"Something"<<endl;
          //                *_lambda = LAMBDA(mydata.data.Young,mydata.data.Poisson);
          //                *_mu = MU(mydata.data.Young,mydata.data.Poisson);
          //                //*_vf = mydata.data.User1;
          //                cout<<"Fibre volume fraction \t"<<mydata.data.Poisson<<endl;
          //                PetscFunctionReturn(0);
          //              }
          //            }
          //          }
          // transversly isotropic material
          if (name.compare(0,22,"MAT_ELASTIC_TRANSISO_1") == 0) {

            //Mat_Elastic_TransIso mydata;
            Mat_Elastic mydata;
            ierr = it->get_attribute_data_structure(mydata); CHKERRQ(ierr);

            Range meshsets;
            rval = moab.get_entities_by_type(it->meshset,MBENTITYSET,meshsets,true); CHKERR_PETSC(rval);
            meshsets.insert(it->meshset);
            for (Range::iterator mit = meshsets.begin(); mit != meshsets.end(); mit++) {
              if (moab.contains_entities(*mit,&ent,1)) {
                //*_E_p   = mydata.data.Youngp;
                //*_E_z   = mydata.data.Youngz;
                //*_nu_p  = mydata.data.Poissonp;
                //*_nu_pz = mydata.data.Poissonpz;
                //*_G_zp  = mydata.data.Shearzp;

                YoungModulus_1 = mydata.data.Young;
                PoissonRatio_1 = mydata.data.Poisson;
                *_lambda = LAMBDA(YoungModulus_1,PoissonRatio_1);
                *_mu     = MU(YoungModulus_1,PoissonRatio_1);

                *_WavinessFactor = mydata.data.ThermalExpansion; // waviness factor

                *_E_p     = mydata.data.User1;
                *_E_z     = mydata.data.User2;
                *_nu_p    = mydata.data.User3;
                *_nu_pz   = mydata.data.User4;
                *_G_zp    = mydata.data.User5;

                *_vf      = mydata.data.User6; // Fibre volume fraction

                *_theta_f = mydata.data.User7; // Misalignment angle
              }
            }
          }
        }
      }
      //cout<<"Get material parameters"<<endl;
      PetscFunctionReturn(0);
    }

    /* =========================================================================
     *
     * Calculate constitutent matrix D
     *
     * =======================================================================*/
    vector< ublas::matrix<FieldData> > D_At_GaussPoint;
    PetscErrorCode calculateD(double _E_p, double _E_z, double _nu_p,
                              double _nu_pz, double _G_zp,
                              double _lambda, double _mu,
                              double _vf,
                              double _theta_f, double _WavinessFactor) {
      PetscFunctionBegin;

			// --------
      // 1. Get stiffness matrix & compliance matrix
      ublas::symmetric_matrix<FieldData,ublas::upper> StiffnessMatrix;
      StiffnessMatrix.resize(6);
      StiffnessMatrix.clear();

      //TransverseIsotropicStiffnessMatrix TranIsoMat(_nu_p,_nu_pz,_E_p,_E_z,_G_zp);

      YarnStiffnessMatrix TranIsoMat(_nu_p, _nu_pz, _E_p, _E_z, _G_zp,
                                     _lambda, _mu, _vf);
      StiffnessMatrix = TranIsoMat.StiffnessMatrix;

      //cout<<"Original stiffness matrix"<<endl;
      // cout<<"Yarn C matrix \t"<<StiffnessMatrix<<endl;

      // Get compliance matrix
      ublas::symmetric_matrix<FieldData,ublas::upper> SMat;
      SMat.resize(6);
      SMat.clear();

      //TransverseIsotropicComplianceMatrix TranIsoMat_S(_nu_p,_nu_pz,_E_p,_E_z,_G_zp);
      YarnComplianceMatrix TranIsoMat_S(_nu_p, _nu_pz, _E_p, _E_z, _G_zp,
                                        _lambda, _mu, _vf);
      SMat = TranIsoMat_S.ComplianceMatrix;
      //cout<<"Original compliance matrix"<<endl;
      // cout<<"Yarn S matrix \t"<<SMat<<endl;

			// --------
      // 2. Get waviness parameter
      ublas::symmetric_matrix<FieldData,ublas::upper> ComplianceMatrix;
      ComplianceMatrix.resize(6);
      ComplianceMatrix.clear();

      double I1, I3, I5, I6, I8;
      double W_mgn, L_wave, alpha_w;
      // W_mgn  =  _amplitude;//1.19;   // unit: mm
      L_wave = 27.9 ;   // unit: mm
      W_mgn  = L_wave*_WavinessFactor;

      alpha_w = 2*M_PI*W_mgn/L_wave;

      I1 = (1+alpha_w*alpha_w/2)/pow((1+alpha_w*alpha_w),1.5);        // m^44
      I3 = (alpha_w*alpha_w/2)/pow((1+alpha_w*alpha_w),1.5);          // m^2 n^2
      I5 = 1-(1+3*alpha_w*alpha_w/2)/pow((1+alpha_w*alpha_w),1.5);    // n^4
      I6 = 1/sqrt(1+alpha_w*alpha_w);                                 // m^2
      I8 = 1-1/sqrt(1+alpha_w*alpha_w);                               // n^2
      // cout<<"alpha ="<<alpha_w<<"\t I1 = "<<I1<<"\t I3 = "<<I3<<"\t I5 = "<<I5<<"\t I6 = "<<I6<<"\t I8 = "<<I8<<endl;

      // --------
      // 3. Calculate transformed compliance matrix
	  // rotate about y-axis
      ComplianceMatrix(0,0) = SMat(0,0)*I1 + (2*SMat(0,2)+SMat(5,5))*I3 + SMat(2,2)*I5;
      ComplianceMatrix(0,1) = SMat(0,1)*I6 + SMat(0,2)*I8;
      ComplianceMatrix(0,2) = SMat(0,2)*(I1+I5) + (SMat(0,0)+SMat(2,2)-SMat(5,5))*I3;
      ComplianceMatrix(1,1) = SMat(0,0);
      ComplianceMatrix(1,2) = SMat(0,2)*I6 + SMat(0,1)*I8;
      ComplianceMatrix(2,2) = SMat(2,2)*I1 + SMat(0,0)*I5 + (2*SMat(0,2)+SMat(5,5))*I3;
      ComplianceMatrix(3,3) = 2*(SMat(0,0)-SMat(0,1))*I6 + SMat(5,5)*I8;
      ComplianceMatrix(4,4) = 2*(2*SMat(0,0)-4*SMat(0,2)+2*SMat(2,2)-SMat(5,5))*I3 + SMat(5,5)*(I1+I5);
      ComplianceMatrix(5,5) = SMat(5,5)*I6 + 2*(SMat(0,0)-SMat(0,1))*I8;

	  // rotate about x-axis
//      ComplianceMatrix(0,0) = SMat(0,0);
//      ComplianceMatrix(0,1) = SMat(0,1)*I6 + SMat(0,2)*I8;
//      ComplianceMatrix(0,2) = SMat(0,2)*I6 + SMat(0,1)*I8;
//      ComplianceMatrix(1,1) = SMat(0,0)*I1 + 2*(SMat(0,0)-SMat(0,1)+SMat(0,2))*I3 + SMat(2,2)*I5;
//      ComplianceMatrix(1,2) = SMat(0,2)*I1 + (SMat(2,2)-SMat(0,0)+2*SMat(0,1))*I3 + SMat(0,2)*I5;
//      ComplianceMatrix(2,2) = SMat(2,2)*I1 + SMat(0,0)*I5 + 2*(SMat(0,0)-SMat(0,1)+SMat(0,2))*I3;
//      ComplianceMatrix(3,3) = 2*(SMat(0,0)-SMat(0,1))*(I1-2*I3+I5)+4*(SMat(0,0)-2*SMat(0,2)+SMat(2,2))*I3;
//      ComplianceMatrix(4,4) = SMat(5,5);
//      ComplianceMatrix(5,5) = SMat(5,5);

      // cout<<"Transformed compliance matrix"<<endl;
      // cout<<ComplianceMatrix<<endl;

      // --------
      // 4. Update the stiffness matrix of yarn with the consideration of waviness
      double SVal;
      double S11, S12, S13, S21, S22, S23, S31, S32, S33;
      S11 = ComplianceMatrix(0,0);
      S12 = ComplianceMatrix(0,1);
      S13 = ComplianceMatrix(0,2);

      S21 = ComplianceMatrix(1,0);
      S22 = ComplianceMatrix(1,1);
      S23 = ComplianceMatrix(1,2);

      S31 = ComplianceMatrix(2,0);
      S32 = ComplianceMatrix(2,1);
      S33 = ComplianceMatrix(2,2);

      SVal = S11*(S22 * S33 - S23 * S32) + S12*(S23 * S31 - S21 * S33)
              + S13*(S21 * S32 - S22 * S31);

      StiffnessMatrix.resize(6);
      StiffnessMatrix.clear();

      StiffnessMatrix(0,0) = (ComplianceMatrix(1,1)*ComplianceMatrix(2,2) - pow(ComplianceMatrix(1,2),2))/SVal;
      StiffnessMatrix(0,1) = (ComplianceMatrix(0,2)*ComplianceMatrix(1,2) - ComplianceMatrix(0,1)*ComplianceMatrix(2,2))/SVal;
      StiffnessMatrix(0,2) = (ComplianceMatrix(0,1)*ComplianceMatrix(1,2) - ComplianceMatrix(0,2)*ComplianceMatrix(1,1))/SVal;

      StiffnessMatrix(1,1) = (ComplianceMatrix(0,0)*ComplianceMatrix(2,2) - pow(ComplianceMatrix(0,2),2))/SVal;
      StiffnessMatrix(1,2) = (ComplianceMatrix(0,1)*ComplianceMatrix(0,2) - ComplianceMatrix(1,2)*ComplianceMatrix(0,0))/SVal;

      StiffnessMatrix(2,2) = (ComplianceMatrix(0,0)*ComplianceMatrix(1,1) - pow(ComplianceMatrix(0,1),2))/SVal;

      StiffnessMatrix(3,3) = 1/ComplianceMatrix(3,3);

      StiffnessMatrix(4,4) = 1/ComplianceMatrix(4,4);

      StiffnessMatrix(5,5) = 1/ComplianceMatrix(5,5);

      //cout<<"Transformed C matrix"<<endl;
	  //cout<<StiffnessMatrix<<endl;

			// ----------
      // 5. Rotating the stiffness matrix according to misalignment information
      //   includeing a set of axes of rotations and their respective angle

      D_At_GaussPoint.resize(coords_at_Gauss_nodes.size());

      vector< ublas::matrix< double > > normalized_phi;
      normalized_phi.resize(coords_at_Gauss_nodes.size());
      ierr = ComputeFibreDirection(normalized_phi); CHKERRQ(ierr);

      for (unsigned int gg=0; gg<coords_at_Gauss_nodes.size(); gg++) {

        int noOfRotations = 1; // Number of rotations ???

        double zVec[3] = {0.0, 0.0, 1.0};
        double AxVector[3]={normalized_phi[gg](0,1)*zVec[2]-normalized_phi[gg](0,2)*zVec[1] , normalized_phi[gg](0,2)*zVec[0]-normalized_phi[gg](0,0)*zVec[2] , normalized_phi[gg](0,0)*zVec[1]-normalized_phi[gg](0,1)*zVec[0]};
        double AxAngle[1]= {asin((sqrt(pow(AxVector[0],2)+pow(AxVector[1],2)+pow(AxVector[2],2)))/(sqrt(pow(normalized_phi[gg](0,0),2)+pow(normalized_phi[gg](0,1),2)+pow(normalized_phi[gg](0,2),2)))*(sqrt(pow(zVec[0],2)+pow(zVec[1],2)+pow(zVec[2],2))))};
        //cout<<"Rotation vector = \t"<<AxVector[0]<<"\t"<<AxVector[1]<<"\t"<<AxVector[2]<<"\t"<<AxAngle[0]<<'\t'<<coords_at_Gauss_nodes.size()<<endl;
		//cout<<"Fibre direction = \t"<<normalized_phi[gg]<<endl;
        //cout<<"Rontation angle = \t"<<AxAngle[0]<<endl;

        double negAxAngle[noOfRotations];
        for (int aa=0; aa<noOfRotations; aa++) negAxAngle[aa]=-AxAngle[aa];

        ublas::matrix<double> DummyMatrix, DummyMatrix2;
        DummyMatrix = ublas::zero_matrix<FieldData>(6,6);
        DummyMatrix = StiffnessMatrix;

        // Rotating stiffness matrix over a number of axis/angle rotations
        for (int aa=0; aa<noOfRotations; aa++) {

          StressTransformation StressRotMat(&AxVector[3*aa], AxAngle[aa]);
          StrainTransformation invStrainRotMat(&AxVector[3*aa], negAxAngle[aa]);

          ublas::matrix<double> TrpMatrixStress;
          TrpMatrixStress = ublas::zero_matrix<FieldData>(6,6);
          TrpMatrixStress = StressRotMat.StressRotMat;

          ublas::matrix<double> TrpMatrixInvStrain;
          TrpMatrixInvStrain = ublas::zero_matrix<FieldData>(6,6);
          TrpMatrixInvStrain = invStrainRotMat.StrainRotMat;

          DummyMatrix2 = ublas::zero_matrix<FieldData>(6,6);
          ublas::matrix< FieldData > dummyA = prod( DummyMatrix, TrpMatrixInvStrain);
          DummyMatrix2 = prod(TrpMatrixStress, dummyA);
          DummyMatrix  = ublas::zero_matrix<FieldData>(6,6);
          DummyMatrix  = DummyMatrix2;
        }

        D_At_GaussPoint[gg].resize(6,6);
        D_At_GaussPoint[gg].clear();
        D_At_GaussPoint[gg] = DummyMatrix;
      }

      PetscFunctionReturn(0);
    }

    /* =========================================================================
     *
     *
     *
     * =======================================================================*/
    PetscErrorCode Fint() {
      PetscFunctionBegin;

      try {

        //Higher order approximation of geometry
        ierr = GetHierarchicalGeometryApproximation(invH,detH); CHKERRQ(ierr);

        double _E_p, _E_z, _nu_p, _nu_pz, _G_zp;
        double _lambda, _mu;
        double _vf;
        double _theta_f;
        double _WavinessFactor;
        //cout<<"Fint"<<endl;
        ierr = GetMatParameters(&_E_p, &_E_z, &_nu_p, &_nu_pz, &_G_zp,
                                &_lambda, &_mu, &_vf, &_theta_f, &_WavinessFactor); CHKERRQ(ierr);
        ierr = calculateD(_E_p, _E_z, _nu_p, _nu_pz, _G_zp,
                          _lambda, _mu, _vf, _theta_f, _WavinessFactor); CHKERRQ(ierr);

        //Gradient at Gauss points;
        vector< ublas::matrix< FieldData > > GradU_at_GaussPt;
        ierr = GetGaussDiffDataVector(fieldName,GradU_at_GaussPt); CHKERRQ(ierr);
        unsigned int g_dim = g_NTET.size()/4;
        assert(GradU_at_GaussPt.size() == g_dim);
        NOT_USED(g_dim);
        vector< ublas::matrix< FieldData > >::iterator viit = GradU_at_GaussPt.begin();
        int gg = 0;
        for(;viit!=GradU_at_GaussPt.end();viit++,gg++) {
          try {
            ublas::matrix< FieldData > GradU = *viit;
            if(!invH.empty()) {
              //GradU =
              //[ dU/dChi1 dU/dChi2 dU/dChi3 ]
              //[ dV/dChi1 dV/dChi2 dU/dChi3 ]
              //[ dW/dChi1 dW/dChi2 dW/dChi3 ]
              //H =
              //[ dX1/dChi1 dX1/dChi2 dX1/dChi3 ]
              //[ dX2/dChi1 dX2/dChi2 dX2/dChi3 ]
              //[ dX3/dChi1 dX3/dChi2 dX3/dChi3 ]
              //invH =
              //[ dChi1/dX1 dChi1/dX2 dChi1/dX3 ]
              //[ dChi2/dX1 dChi2/dX2 dChi2/dX3 ]
              //[ dChi3/dX1 dChi3/dX2 dChi3/dX3 ]
              //GradU =
              //[ dU/dX1 dU/dX2 dU/dX3 ]
              //[ dV/dX1 dV/dX2 dV/dX3 ] = GradU * invH
              //[ dW/dX1 dW/dX2 dW/dX3 ]
              GradU = prod( GradU, invH[gg] );
            }
            ublas::matrix< FieldData > Strain = 0.5*( GradU + trans(GradU) );
            ublas::vector< FieldData > VoightStrain(6);
            VoightStrain[0] = Strain(0,0);
            VoightStrain[1] = Strain(1,1);
            VoightStrain[2] = Strain(2,2);
            VoightStrain[3] = 2*Strain(0,1);
            VoightStrain[4] = 2*Strain(1,2);
            VoightStrain[5] = 2*Strain(2,0);
            double w = V*G_TET_W[gg];
            ublas::vector<FieldData> VoightStress = prod(w*D_At_GaussPoint[gg],VoightStrain);
            //BT * VoigtStress
            f_int.resize(row_mat);
            for(int rr = 0;rr<row_mat;rr++) {
              if(RowGlob[rr].size()==0) continue;
              ublas::matrix<FieldData> &B = (rowBMatrices[rr])[gg];
              if(gg == 0) {
                f_int[rr] = prod( trans(B), VoightStress );
              } else {
                f_int[rr] += prod( trans(B), VoightStress );
              }
            }
          } catch (const std::exception& ex) {
            ostringstream ss;
            ss << "thorw in method: " << ex.what() << " at line " << __LINE__ << " in file " << __FILE__ << endl;
            SETERRQ(PETSC_COMM_SELF,1,ss.str().c_str());
          }
        }

      } catch (const std::exception& ex) {
        ostringstream ss;
        ss << "thorw in method: " << ex.what() << " at line " << __LINE__ << " in file " << __FILE__ << endl;
        SETERRQ(PETSC_COMM_SELF,1,ss.str().c_str());
      }

      PetscFunctionReturn(0);
    }

    /* =========================================================================
     *
     * Element stiffness matrix
     *
     * =======================================================================*/
    PetscErrorCode Stiffness() {
      PetscFunctionBegin;

      double _E_p, _E_z, _nu_p, _nu_pz, _G_zp;
      double _lambda, _mu;
      double _vf;
      double _theta_f;
      double _WavinessFactor;
      // -----------
      // 1. Get material parameters
      //cout<<"Stiffness"<<endl;
      ierr = GetMatParameters(&_E_p, &_E_z, &_nu_p, &_nu_pz, &_G_zp,
                              &_lambda, &_mu, &_vf, &_theta_f, &_WavinessFactor); CHKERRQ(ierr);
      // -----------
      // 2. Calculate constitutive matrix
      ierr = calculateD(_E_p, _E_z, _nu_p, _nu_pz, _G_zp,
                        _lambda, _mu, _vf, _theta_f, _WavinessFactor); CHKERRQ(ierr);

      // -----------
      // 3. Calculate stiffness matrix
      //cout<<"K calculation"<<endl;
      K.resize(row_mat,col_mat);
      int g_dim = g_NTET.size()/4;
      for (int rr = 0; rr<row_mat; rr++) {
        if (RowGlob[rr].size()==0) continue;
        for (int gg=0; gg<g_dim; gg++) {

          //cout<<"gg = \t"<<gg<<endl;
          //cout<<"D matrix"<<endl;cout<<D_At_GaussPoint[gg]<<endl;
          ublas::matrix<FieldData> &row_Mat = (rowBMatrices[rr])[gg];
          double w = V*G_TET_W[gg];
          if (detH.size()>0) {
            w *= detH[gg];
          }
          BD.resize(6,row_Mat.size2());
          //cout<<"BD_Resize"<<endl;
          cblas_dsymm(CblasRowMajor,CblasLeft,CblasUpper,
                      BD.size1(),BD.size2(),
                      w,&*D_At_GaussPoint[gg].data().begin(),D_At_GaussPoint[gg].size2(),
                      &*row_Mat.data().begin(),row_Mat.size2(),
                      0.,&*BD.data().begin(),BD.size2());
          for(int cc = rr;cc<col_mat;cc++) {
            if(ColGlob[cc].size()==0) continue;
            ublas::matrix<FieldData> &col_Mat = (colBMatrices[cc])[gg];
            if(gg == 0) {
              K(rr,cc).resize(BD.size2(),col_Mat.size2());
              //ublas::noalias(K(rr,cc)) = prod(trans(BD) , col_Mat ); // int BT*D*B
              cblas_dgemm(CblasRowMajor,CblasTrans,CblasNoTrans,
                          BD.size2(),col_Mat.size2(),BD.size1(),
                          1.,&*BD.data().begin(),BD.size2(),
                          &*col_Mat.data().begin(),col_Mat.size2(),
                          0.,&*K(rr,cc).data().begin(),K(rr,cc).size2());
            } else {
              //ublas::noalias(K(rr,cc)) += prod(trans(BTD) , col_Mat ); // int BT*D*B
              cblas_dgemm(CblasRowMajor,CblasTrans,CblasNoTrans,
                          BD.size2(),col_Mat.size2(),BD.size1(),
                          1.,&*BD.data().begin(),BD.size2(),
                          &*col_Mat.data().begin(),col_Mat.size2(),
                          1.,&*K(rr,cc).data().begin(),K(rr,cc).size2());
            }
          }
        }
      }

      //cout<<"K calculation finish"<<endl;

      PetscFunctionReturn(0);
    }

    /* =========================================================================
     *
     * Compute fibre direction
     *
     * =======================================================================*/
    PetscErrorCode ComputeFibreDirection(vector<ublas::matrix<double> > &normalized_phi) {
      PetscFunctionBegin;
      ierr = OpStudentStart_TET(g_NTET); CHKERRQ(ierr);

      EntityHandle fe_handle = fePtr->get_ent();

      Range tetNodes;
      rval = moab.get_connectivity(&fe_handle,1,tetNodes); CHKERR_THROW(rval);

      vector< ublas::matrix< FieldData > > phi;
      ierr = GetGaussDiffDataVector("POTENTIAL_FIELD",phi); CHKERRQ(ierr);
      double fibreVector[3];

      for (unsigned int gg=0; gg<phi.size(); gg++) {
        normalized_phi[gg].resize(1,3);
        for (int ii=0; ii<3; ii++) {
          normalized_phi[gg](0,ii) = -phi[gg](0,ii)/sqrt(pow(phi[gg](0,0),2)+pow(phi[gg](0,1),2)+pow(phi[gg](0,2),2));
          fibreVector[ii] = normalized_phi[0](0,ii);
        }
      }

      for(Range::iterator niit1 = tetNodes.begin();niit1!=tetNodes.end();niit1++){
        rval = moab.tag_set_data(th_fibre_dir,&*niit1,1,&fibreVector[0]); CHKERR_PETSC(rval);
      }

      ierr = OpStudentEnd(); CHKERRQ(ierr);
      PetscFunctionReturn(0);
    }

  };

  struct TranIso_PostProc_Waviness_OnRefMesh: public PostProcDisplacemenysAndStarinAndElasticLinearStressOnRefMesh {
    double E_p, E_z, nu_p, nu_pz, G_zp;
    double lambda, mu;
    double vf;
    Tag th_fibre_orientation;

    TranIso_PostProc_Waviness_OnRefMesh (FieldInterface& _mField,
                                             double _lambda, double _mu,
                                             double _E_p, double _E_z,
                                             double _nu_p, double _nu_pz,
                                             double _G_zp):
    PostProcDisplacemenysAndStarinAndElasticLinearStressOnRefMesh (_mField, "DISPLACEMENT",
                                                                   _lambda, _mu),
    E_p(_E_p),E_z(_E_z),nu_p(_nu_p),nu_pz(_nu_pz),G_zp(_G_zp) {

      double def_VAL2[3] = {0,0,0};
      rval = moab_post_proc.tag_get_handle("FIBRE_DIRECTION",3,MB_TYPE_DOUBLE,
                                           th_fibre_orientation,MB_TAG_CREAT|MB_TAG_SPARSE,
                                           &def_VAL2); CHKERR_THROW(rval);
    }

    PetscErrorCode ComputeGradient(vector<ublas::matrix<double> > &normalized_phi) {

      PetscFunctionBegin;
      ierr = OpStudentStart_TET(g_NTET); CHKERRQ(ierr);

      EntityHandle fe_handle = fePtr->get_ent();

      Range tetNodes;
      rval = moab.get_connectivity(&fe_handle,1,tetNodes); CHKERR_THROW(rval);

      vector< ublas::matrix< FieldData> > phi;
      ierr = GetGaussDiffDataVector("POTENTIAL_FIELD",phi); CHKERRQ(ierr);

      for (unsigned int gg=0; gg<phi.size(); gg++) {
        normalized_phi[gg].resize(1,3);
        for (int ii=0; ii<3; ii++) {
          normalized_phi[gg](0,ii) = -phi[gg](0,ii)/sqrt(pow(phi[gg](0,0),2)+
                                                         pow(phi[gg](0,1),2)+
                                                         pow(phi[gg](0,2),2));
        }
      }

      ierr = OpStudentEnd(); CHKERRQ(ierr);

      PetscFunctionReturn(0);
    }

    PetscErrorCode operator()() {
      PetscFunctionBegin;

      ierr = do_operator(); CHKERRQ(ierr);
      ierr = OpStudentStart_TET(g_NTET); CHKERRQ(ierr);

      // Get stifffness matrix
      ublas::symmetric_matrix<FieldData,ublas::upper> StiffnessMatrix;
      StiffnessMatrix.resize(6);
      StiffnessMatrix.clear();
      YarnStiffnessMatrix TranIsoMat(nu_p, nu_pz, E_p, E_z, G_zp, lambda, mu, vf);
      StiffnessMatrix = TranIsoMat.StiffnessMatrix;

      int gg=0;
      vector< ublas::matrix< FieldData > > GradU_at_GaussPt;
      ierr = GetGaussDiffDataVector(field_name,GradU_at_GaussPt); CHKERRQ(ierr);
      vector< ublas::matrix< FieldData > >::iterator viit = GradU_at_GaussPt.begin();
      map<EntityHandle, EntityHandle>::iterator mit = node_map.begin();

      vector< ublas::matrix< double > > normalized_phi;
      normalized_phi.resize(GradU_at_GaussPt.size());
      ierr = ComputeGradient(normalized_phi); CHKERRQ(ierr);

      cout<<GradU_at_GaussPt.size()<<endl;

      for (; viit!=GradU_at_GaussPt.end(); viit++, mit++, gg++) {

        double zVec[3]={0.0,0.0,1.0};

        double AxVector[3]={normalized_phi[gg](0,1)*zVec[2]-normalized_phi[gg](0,2)*zVec[1] , normalized_phi[gg](0,2)*zVec[0]-normalized_phi[gg](0,0)*zVec[2] , normalized_phi[gg](0,0)*zVec[1]-normalized_phi[gg](0,1)*zVec[0]};
        double AxAngle[1]= {asin((sqrt(pow(AxVector[0],2)+pow(AxVector[1],2)+pow(AxVector[2],2)))/(sqrt(pow(normalized_phi[gg](0,0),2)+pow(normalized_phi[gg](0,1),2)+pow(normalized_phi[gg](0,2),2)))*(sqrt(pow(zVec[0],2)+pow(zVec[1],2)+pow(zVec[2],2))))};

        ///Rotating the Stiffness matrix according a set of axes of rotations and their respective angle

        int noOfRotations = 1; //Number of Rotations
        double negAxAngle[noOfRotations];
        for (int aa=0; aa<noOfRotations; aa++) negAxAngle[aa]=-AxAngle[aa];

        ublas::matrix<double> DummyMatrix,DummyMatrix2;
        DummyMatrix = ublas::zero_matrix<FieldData>(6,6);
        DummyMatrix = StiffnessMatrix;

        ///Rotating Stiffness over a number of axis/angle rotations
        for (int aa=0; aa<noOfRotations; aa++) {

          StressTransformation StressRotMat(&AxVector[3*aa], AxAngle[aa]);
          StrainTransformation invStrainRotMat(&AxVector[3*aa], negAxAngle[aa]);

          ublas::matrix<double> TrpMatrixStress;
          TrpMatrixStress = ublas::zero_matrix<FieldData>(6,6);
          TrpMatrixStress = StressRotMat.StressRotMat;

          ublas::matrix<double> TrpMatrixInvStrain;
          TrpMatrixInvStrain = ublas::zero_matrix<FieldData>(6,6);
          TrpMatrixInvStrain = invStrainRotMat.StrainRotMat;

          DummyMatrix2 = ublas::zero_matrix<FieldData>(6,6);
          ublas::matrix< FieldData > dummyA = prod( DummyMatrix , TrpMatrixInvStrain );
          DummyMatrix2 = prod(TrpMatrixStress,dummyA);
          DummyMatrix = ublas::zero_matrix<FieldData>(6,6);
          DummyMatrix = DummyMatrix2;
        }

        D.resize(6,6);
        D.clear();
        D = DummyMatrix;

        ///Compute Strains and save them on TAG
        ublas::matrix< FieldData > GradU = *viit;
        ublas::matrix< FieldData > Strain = 0.5*( GradU + trans(GradU) );
        rval = moab_post_proc.tag_set_data(th_strain,&mit->second,1,&(Strain.data()[0])); CHKERR_PETSC(rval);

        ublas::matrix<double> AARotMatrix;
        AARotMatrix = ublas::identity_matrix<FieldData>(3);

        for (int aa=0; aa<noOfRotations; aa++) {

          AxisAngleRotationalMatrix RotMatrix(&AxVector[3*aa], AxAngle[aa]);

          ublas::matrix<double> rotationalMat;
          rotationalMat = ublas::zero_matrix<FieldData>(3,3);
          rotationalMat=RotMatrix.AARotMat;

          ublas::matrix<double> AARotMatrix1;
          AARotMatrix1 = prod(rotationalMat,AARotMatrix);
          AARotMatrix = ublas::zero_matrix<FieldData>(3,3);
          AARotMatrix = AARotMatrix1;
        }

        ///Rotate AxisYVector[0,1,0] to the direction of the fibre and save in TAG
        ublas::vector<FieldData> AxisYVector(3);
        AxisYVector[0]=0; AxisYVector[1]=0;AxisYVector[2]=1;
        ublas::vector<FieldData> Fibre = prod(AARotMatrix,AxisYVector);

        rval = moab_post_proc.tag_set_data(th_fibre_orientation,&mit->second,1,&Fibre[0]); CHKERR_PETSC(rval);

        ///calculate stress and save it into tag
        ublas::vector<FieldData> Strain_VectorNotation(6);
        Strain_VectorNotation[0] = Strain(0,0);
        Strain_VectorNotation[1] = Strain(1,1);
        Strain_VectorNotation[2] = Strain(2,2);
        Strain_VectorNotation[3] = 2*Strain(0,1);
        Strain_VectorNotation[4] = 2*Strain(1,2);
        Strain_VectorNotation[5] = 2*Strain(2,0);
        ublas::vector< FieldData > Stress_VectorNotation = prod( D, Strain_VectorNotation );
        ublas::matrix< FieldData > Stress = ublas::zero_matrix<FieldData>(3,3);
        Stress(0,0) = Stress_VectorNotation[0];
        Stress(1,1) = Stress_VectorNotation[1];
        Stress(2,2) = Stress_VectorNotation[2];
        Stress(0,1) = Stress(1,0) = Stress_VectorNotation[3];
        Stress(1,2) = Stress(2,1) = Stress_VectorNotation[4];
        Stress(2,0) = Stress(0,2) = Stress_VectorNotation[5];

        rval = moab_post_proc.tag_set_data(th_stress,&mit->second,1,&(Stress.data()[0])); CHKERR_PETSC(rval);
      }

      ierr = OpStudentEnd(); CHKERRQ(ierr);
      PetscFunctionReturn(0);
    }
  };

}
#endif //__ELASTICFEMETHODTRANSISO_HPP__
