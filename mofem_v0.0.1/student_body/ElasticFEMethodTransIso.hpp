/* Copyright (C) 2013, Michel Cortis
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
 * \brief Function to Calculate Strain Transformation Matrix
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

struct TranIsotropicElasticFEMethod: public ElasticFEMethod {
    
    double E_p, E_z, nu_p, nu_pz, G_zp;
    
    TranIsotropicElasticFEMethod(
                                 FieldInterface& _mField,BaseDirihletBC *_dirihlet_ptr,Mat &_Aij,Vec& _D,Vec& _F,
                                 double _lambda,double _mu,double _E_p,double _E_z, double _nu_p,double _nu_pz, double _G_zp): 
    ElasticFEMethod(_mField,_dirihlet_ptr,_Aij,_D,_F,_lambda,_mu), E_p(_E_p), E_z(_E_z), nu_p(_nu_p), nu_pz(_nu_pz), G_zp(_G_zp)  {
        
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
        
        ///Insert Axes and Angles of Rotation
        int noOfRotations = 1; //Number of Rotations
        double AxVector[3*noOfRotations]; //First Rotational Axis
        AxVector[0]=1; AxVector[1]=0; AxVector[2]=0; //First Rotational Axis 
        //        AxVector[3]=1; AxVector[4]=0; AxVector[5]=0; //Second Rotational Axis
        double AxAngle[noOfRotations];
        AxAngle[0] = -0.25*M_PI; //First Rotational Angle
        //        AxAngle[1] = -0.25*M_PI; //Second Rotational Angle
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
}
#endif //__ELASTICFEMETHODTRANSISO_HPP__

