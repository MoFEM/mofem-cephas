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

namespace MoFEM {

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
        double delta=((1+nu_p)*(1-nu_p-(2*nu_pz*nu_zp)))/(E_p*E_p*E_z);
        
        StiffnessMatrix.resize(6);
        StiffnessMatrix.clear();
        StiffnessMatrix(0,0)=StiffnessMatrix(1,1)=(1-nu_pz*nu_zp)/(E_p*E_z*delta);
        StiffnessMatrix(2,2)=(1-nu_p*nu_p)/(E_p*E_p*delta);
        
        StiffnessMatrix(0,1)=StiffnessMatrix(1,0)=(nu_p+nu_zp*nu_pz)/(E_p*E_z*delta);
        StiffnessMatrix(0,2)=StiffnessMatrix(2,0)=StiffnessMatrix(1,2)=StiffnessMatrix(2,1)=(nu_zp+nu_p*nu_zp)/(E_p*E_z*delta);
        
        StiffnessMatrix(3,3)=E_p/(2*(1+nu_p));
        StiffnessMatrix(4,4)=StiffnessMatrix(5,5)=G_zp;
    };
};

struct TransverseIsotropicComplianceMatrix {
    
    double nu_p, nu_pz, E_p, E_z, G_zp;
    
    ublas::symmetric_matrix<FieldData,ublas::upper> ComplianceMatrix;
    
    TransverseIsotropicComplianceMatrix(double nu_p, double nu_pz, double E_p, double E_z, double G_zp){
        
        double nu_zp = (nu_pz*E_z)/E_p;
        double Gp = E_p/(2*(1+nu_p));
        
        ComplianceMatrix.resize(6);
        ComplianceMatrix.clear();
        ComplianceMatrix(0,0)=ComplianceMatrix(1,1)=1/E_p;
        ComplianceMatrix(2,2)=1/E_z;
        
        ComplianceMatrix(0,1)=ComplianceMatrix(1,0)=-nu_p/E_p;
        ComplianceMatrix(0,2)=ComplianceMatrix(2,0)=-nu_zp/E_z;
        ComplianceMatrix(1,2)=ComplianceMatrix(2,1)=-nu_pz/E_z;
        
        ComplianceMatrix(3,3)=1/Gp;
        ComplianceMatrix(4,4)=ComplianceMatrix(5,5)=1/G_zp;
    };
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
    };
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
        
    };
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
        
    };  
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
        
        
    };  
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
          
    double E_p, E_z, nu_p, nu_pz, G_zp;
    int noAA;
    double *AxVector, *AxAngle;
    
    TranIsotropicAxisAngleRotElasticFEMethod(
                                     FieldInterface& _mField,BaseDirihletBC *_dirihlet_ptr,Mat &_Aij,Vec& _D,Vec& _F,
                                     double _lambda,double _mu,double _E_p,double _E_z, double _nu_p,double _nu_pz, double _G_zp, int _noAA, double *_AxVector, double *_AxAngle): 
    ElasticFEMethod(_mField,_dirihlet_ptr,_Aij,_D,_F,_lambda,_mu), E_p(_E_p), E_z(_E_z), nu_p(_nu_p), nu_pz(_nu_pz), G_zp(_G_zp), noAA(_noAA), AxVector(_AxVector), AxAngle(_AxAngle)  {
            
    };
        
    PetscErrorCode preProcess() {
        PetscFunctionBegin;
        ierr = ElasticFEMethod::preProcess();  CHKERRQ(ierr);
        PetscSynchronizedFlush(PETSC_COMM_WORLD);
        ierr = MatAssemblyBegin(Aij,MAT_FLUSH_ASSEMBLY); CHKERRQ(ierr);
        ierr = MatAssemblyEnd(Aij,MAT_FLUSH_ASSEMBLY); CHKERRQ(ierr);
        
        PetscFunctionReturn(0);
        }
        
    PetscErrorCode Fint(Vec F_int) {
        PetscFunctionBegin;
        ierr = ElasticFEMethod::Fint(); CHKERRQ(ierr);
        for(int rr = 0;rr<row_mat;rr++) {
            if(RowGlob[rr].size()!=f_int[rr].size()) SETERRQ(PETSC_COMM_SELF,1,"data inconsistency");
            if(RowGlob[rr].size()==0) continue;
            f_int[rr] *= -1; //This is not SNES we solve K*D = -RES
            ierr = VecSetValues(F_int,RowGlob[rr].size(),&(RowGlob[rr])[0],&(f_int[rr].data()[0]),ADD_VALUES); CHKERRQ(ierr);
        }
        PetscFunctionReturn(0);
    }
        
    PetscErrorCode calulateD(double _lambda,double _mu) {
        PetscFunctionBegin;
        PetscFunctionReturn(0);
    }
        
    PetscErrorCode operator()() {
        PetscFunctionBegin;
        ierr = OpStudentStart_TET(g_NTET); CHKERRQ(ierr);
        ierr = GetMatrices(); CHKERRQ(ierr);
            
        //Dirihlet Boundary Condition
        ierr = dirihlet_bc_method_ptr->SetDirihletBC_to_ElementIndicies(this,RowGlob,ColGlob,DirihletBC); CHKERRQ(ierr);
            
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
            
        //Assembly Aij and F
        ierr = RhsAndLhs(); CHKERRQ(ierr);
        
        //Neumann Boundary Conditions
        ierr = NeumannBC(F); CHKERRQ(ierr);
            
        ierr = OpStudentEnd(); CHKERRQ(ierr);
            
        PetscFunctionReturn(0); }
        
};

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
        
    double E_p, E_z, nu_p, nu_pz, G_zp;
    Tag th_fibre_dir;
    
    TranIsotropicFibreDirRotElasticFEMethod(
                                            FieldInterface& _mField,BaseDirihletBC *_dirihlet_ptr,Mat &_Aij,Vec& _D,Vec& _F,
                                            double _lambda,double _mu,double _E_p,double _E_z, double _nu_p,double _nu_pz, double _G_zp): 
    ElasticFEMethod(_mField,_dirihlet_ptr,_Aij,_D,_F,_lambda,_mu), E_p(_E_p), E_z(_E_z), nu_p(_nu_p), nu_pz(_nu_pz), G_zp(_G_zp)  {
            
        double def_VAL2[3] = {0,0,0};
        // create a new tag th_fibre_dir
        rval = moab.tag_get_handle( "POT_FLOW_FIBRE_DIR",3,MB_TYPE_DOUBLE,th_fibre_dir,MB_TAG_CREAT|MB_TAG_SPARSE,&def_VAL2); CHKERR_THROW(rval);
    };
        
    vector< ublas::matrix<FieldData> > D_At_GaussPoint;
    
    PetscErrorCode preProcess() {
        PetscFunctionBegin;
        ierr = ElasticFEMethod::preProcess();  CHKERRQ(ierr);
        PetscSynchronizedFlush(PETSC_COMM_WORLD);
        ierr = MatAssemblyBegin(Aij,MAT_FLUSH_ASSEMBLY); CHKERRQ(ierr);
        ierr = MatAssemblyEnd(Aij,MAT_FLUSH_ASSEMBLY); CHKERRQ(ierr);
            
        PetscFunctionReturn(0);
    }
    
    PetscErrorCode calulateD(double _lambda,double _mu) {
        PetscFunctionBegin;
        PetscFunctionReturn(0);
    }
 
    PetscErrorCode GetMatParameters(double *_lambda,double *_mu) {
      PetscFunctionBegin;

      *_lambda = lambda;
      *_mu = mu;

      //FIXME
      //Here You should Get Parameters from atributes of Mat_TransIsoSet.
      //All parameters which you geting in constructor can be retrived for each element,
      //	this will allow to have problem, whith heterogenous material properties,
      //	f.e. each rope can have diffrent properties
      //Look to ElasticFEMethid how is done.

      PetscFunctionReturn(0);
    }
   
    PetscErrorCode Fint() {
        PetscFunctionBegin;
        
        double _lambda,_mu;
        ierr = GetMatParameters(&_lambda,&_mu); CHKERRQ(ierr);
        ierr = calulateD(_lambda,_mu); CHKERRQ(ierr); 

        //Gradient at Gauss points; 
        vector< ublas::matrix< FieldData > > GradU_at_GaussPt;
        ierr = GetGaussDiffDataVector("DISPLACEMENT",GradU_at_GaussPt); CHKERRQ(ierr);
        unsigned int g_dim = g_NTET.size()/4;
        assert(GradU_at_GaussPt.size() == g_dim);
        NOT_USED(g_dim);
        vector< ublas::matrix< FieldData > >::iterator viit = GradU_at_GaussPt.begin();
        int gg = 0;
        for(;viit!=GradU_at_GaussPt.end();viit++,gg++) {
            ublas::matrix< FieldData > GradU = *viit;
            ublas::matrix< FieldData > Strain = 0.5*( GradU + trans(GradU) );
            ublas::vector< FieldData > VoightStrain(6);
            VoightStrain[0] = Strain(0,0);
            VoightStrain[1] = Strain(1,1);
            VoightStrain[2] = Strain(2,2);
            VoightStrain[3] = 2*Strain(0,1);
            VoightStrain[4] = 2*Strain(1,2);
            VoightStrain[5] = 2*Strain(2,0);
            double w = V*G_W_TET[gg];
            ublas::vector<FieldData> VoightStress = prod(w*D_At_GaussPoint[gg],VoightStrain);
            //BT * VoigtStress
            for(int rr = 0;rr<row_mat;rr++) {
                f_int.resize(row_mat);
                ublas::matrix<FieldData> &B = (rowBMatrices[rr])[gg];
                if(gg == 0) {
                    f_int[rr] = prod( trans(B), VoightStress );
                } else {
                    f_int[rr] += prod( trans(B), VoightStress );
                }
            }
        }
        
        PetscFunctionReturn(0);
    }
        
    PetscErrorCode Fint(Vec F_int) {
        PetscFunctionBegin;
        ierr = TranIsotropicFibreDirRotElasticFEMethod::Fint(); CHKERRQ(ierr);
        for(int rr = 0;rr<row_mat;rr++) {
            if(RowGlob[rr].size()!=f_int[rr].size()) SETERRQ(PETSC_COMM_SELF,1,"data inconsistency");
            if(RowGlob[rr].size()==0) continue;
            f_int[rr] *= -1; //This is not SNES we solve K*D = -RES
            ierr = VecSetValues(F_int,RowGlob[rr].size(),&(RowGlob[rr])[0],&(f_int[rr].data()[0]),ADD_VALUES); CHKERRQ(ierr);

        }
        PetscFunctionReturn(0);
    }

    PetscErrorCode Stiffness() {
        PetscFunctionBegin;
        
        double _lambda,_mu;
        ierr = GetMatParameters(&_lambda,&_mu); CHKERRQ(ierr);
        ierr = calulateD(_lambda,_mu); CHKERRQ(ierr);
        
        K.resize(row_mat,col_mat);
        int g_dim = g_NTET.size()/4;
        for(int rr = 0;rr<row_mat;rr++) {
            for(int gg = 0;gg<g_dim;gg++) {
                ublas::matrix<FieldData> &row_Mat = (rowBMatrices[rr])[gg];
                double w = V*G_W_TET[gg];
                BD.resize(6,row_Mat.size2());
                //ublas::noalias(BD) = prod( w*D,row_Mat );
                cblas_dsymm(CblasRowMajor,CblasLeft,CblasUpper,
                            BD.size1(),BD.size2(),
                            w,&*D_At_GaussPoint[gg].data().begin(),D_At_GaussPoint[gg].size2(),
                            &*row_Mat.data().begin(),row_Mat.size2(),
                            0.,&*BD.data().begin(),BD.size2());
                for(int cc = rr;cc<col_mat;cc++) {
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

    PetscErrorCode ComputeFibreDirection(double *fibreVector) {
        PetscFunctionBegin;
        ierr = OpStudentStart_TET(g_NTET); CHKERRQ(ierr);
        
        EntityHandle fe_handle = fe_ptr->get_ent();
        
        Range tetNodes;
        rval = moab.get_connectivity(&fe_handle,1,tetNodes); CHKERR_THROW(rval);
        
        vector< ublas::matrix< FieldData > > phi;  // size is 45<3>
        
        ierr = GetGaussDiffDataVector("POTENTIAL_FIELD",phi); CHKERRQ(ierr);  // calculate the derivatives for potential field for all the gauss points of one element
        //cout<<"\n\n";
        //cout<<"phi.size() "<<phi.size();
        //cout<<phi[0] <<"\n";
        //cout<<"\n\n";
        
        for (int ii=0; ii<3; ii++) fibreVector[ii] = -phi[0](0,ii)/sqrt(pow(phi[0](0,0),2)+pow(phi[0](0,1),2)+pow(phi[0](0,2),2)); //normalizing the vector (to find a unit vector)
        
        for(Range::iterator niit1 = tetNodes.begin();niit1!=tetNodes.end();niit1++){   // this is to just save the fiber dirction data on the elements
            rval = moab.tag_set_data(th_fibre_dir,&*niit1,1,&fibreVector[0]); CHKERR_PETSC(rval);  
        }
        
        ierr = OpStudentEnd(); CHKERRQ(ierr);
        PetscFunctionReturn(0); 
    }  
    
    PetscErrorCode operator()() {
        PetscFunctionBegin;
        ierr = OpStudentStart_TET(g_NTET); CHKERRQ(ierr);
        ierr = GetMatrices(); CHKERRQ(ierr);

        //Dirihlet Boundary Condition
        ierr = dirihlet_bc_method_ptr->SetDirihletBC_to_ElementIndicies(this,RowGlob,ColGlob,DirihletBC); CHKERRQ(ierr);
            
        ///Get Stiffness Matrix
        ublas::symmetric_matrix<FieldData,ublas::upper> StiffnessMatrix;
        StiffnessMatrix.resize(6);
        StiffnessMatrix.clear();
        TransverseIsotropicStiffnessMatrix TranIsoMat(nu_p, nu_pz, E_p, E_z, G_zp);
        StiffnessMatrix=TranIsoMat.StiffnessMatrix;
//        IsotropicStiffnessMatrix IsoMat(lambda, mu);
//        StiffnessMatrix=IsoMat.StiffnessMatrix;
            
        ///Rotating the Stiffness matrix according a set of axes of rotations and their respective angle
                   
        D_At_GaussPoint.resize(coords_at_Gauss_nodes.size());
        
        for(unsigned int gg=0;gg<coords_at_Gauss_nodes.size();gg++){
            
            double fVec[3];
        ierr = ComputeFibreDirection(&fVec[0]); CHKERRQ(ierr);
            
            int noOfRotations = 1; //Number of Rotations

            double zVec[3]={0.0,0.0,1.0};
            double AxVector[3]={fVec[1]*zVec[2]-fVec[2]*zVec[1] , fVec[2]*zVec[0]-fVec[0]*zVec[2] , fVec[0]*zVec[1]-fVec[1]*zVec[0]};
            double AxAngle[1]= {asin((sqrt(pow(AxVector[0],2)+pow(AxVector[1],2)+pow(AxVector[2],2)))/(sqrt(pow(fVec[0],2)+pow(fVec[1],2)+pow(fVec[2],2)))*(sqrt(pow(zVec[0],2)+pow(zVec[1],2)+pow(zVec[2],2))))};
                        
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
        
        //Assembly Aij and F
        ierr = RhsAndLhs(); CHKERRQ(ierr);
            
        //Neumann Boundary Conditions
        ierr = NeumannBC(F); CHKERRQ(ierr);
        
        ierr = OpStudentEnd(); CHKERRQ(ierr);
        
        PetscFunctionReturn(0); }
        
};

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
        
    TranIso_PostProc_AxisAngle_OnRefMesh( Interface& _moab,double _lambda,double _mu, double _E_p,double _E_z, double _nu_p, double _nu_pz, double _G_zp, int _noAA, double *_AxVector, double *_AxAngle):
    PostProcDisplacemenysAndStarinAndElasticLinearStressOnRefMesh(_moab,_lambda,_mu),E_p(_E_p),E_z(_E_z),nu_p(_nu_p),nu_pz(_nu_pz),G_zp(_G_zp), noAA(_noAA), AxVector(_AxVector), AxAngle(_AxAngle) {
            
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
        
    TranIso_PostProc_FibreDirRot_OnRefMesh( Interface& _moab,double _lambda,double _mu, double _E_p,double _E_z, double _nu_p, double _nu_pz, double _G_zp):
    PostProcDisplacemenysAndStarinAndElasticLinearStressOnRefMesh(_moab,_lambda,_mu),E_p(_E_p),E_z(_E_z),nu_p(_nu_p),nu_pz(_nu_pz),G_zp(_G_zp) {
            
        double def_VAL2[3] = {0,0,0};
        rval = moab_post_proc.tag_get_handle("FIBRE_DIRECTION",3,MB_TYPE_DOUBLE,th_fibre_orientation,MB_TAG_CREAT|MB_TAG_SPARSE,&def_VAL2); CHKERR_THROW(rval);
            
    };
       
    PetscErrorCode ComputeGradient(double *fibreVector) {
        PetscFunctionBegin;
        ierr = OpStudentStart_TET(g_NTET); CHKERRQ(ierr);
                
        vector< ublas::matrix< FieldData > > phi;
        ierr = GetGaussDiffDataVector("POTENTIAL_FIELD",phi); CHKERRQ(ierr);  
        for (int ii=0; ii<3; ii++) fibreVector[ii] = phi[0](0,ii);

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
            
        for(;viit!=GradU_at_GaussPt.end();viit++,mit++,gg++) {
            
            double fVec[3];
            ierr = ComputeGradient(&fVec[0]); CHKERRQ(ierr);
            
            double zVec[3]={0.0,0.0,1.0};
            double AxVector[3]={fVec[1]*zVec[2]-fVec[2]*zVec[1] , fVec[2]*zVec[0]-fVec[0]*zVec[2] , fVec[0]*zVec[1]-fVec[1]*zVec[0]};
            double AxAngle[1]= {asin((sqrt(pow(AxVector[0],2)+pow(AxVector[1],2)+pow(AxVector[2],2)))/(sqrt(pow(fVec[0],2)+pow(fVec[1],2)+pow(fVec[2],2)))*(sqrt(pow(zVec[0],2)+pow(zVec[1],2)+pow(zVec[2],2))))};
            
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

}
#endif //__ELASTICFEMETHODTRANSISO_HPP__

