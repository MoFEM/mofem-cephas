/* Copyright (C) 2013, Lukasz Kaczmarczyk (likask AT wp.pl)
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

#include "FieldInterface.hpp"
#include "FieldCore.hpp"
#include "moabFEMethod_UpLevelStudent.hpp"
#include "cholesky.hpp"
#include <petscksp.h>

#include "ElasticFEMethod.hpp"
#include "PostProcVertexMethod.hpp"
#include "PostProcDisplacementAndStrainOnRefindedMesh.hpp"

using namespace MoFEM;

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

#define RND_EPS 1e-6
double roundn(double n) {
    
    //break n into fractional part (fract) and integral part (intp)
    double fract, intp;
    fract = modf(n,&intp);
    
    // case where n approximates zero, set n to "positive" zero
    if (abs(intp)==0) {
        if(abs(fract)<=RND_EPS) {
            n=0.000;
        }
    }
    
    return n;
}

struct TranIsotropic_Fibre_PostProcDisplacemenysAndStarinAndElasticLinearStressOnRefMesh: public PostProcDisplacemenysAndStarinAndElasticLinearStressOnRefMesh {
    
    double E_p, E_z, nu_p, nu_pz, G_zp;
    Tag th_fibre_orientation;
    ofstream myfile;

    TranIsotropic_Fibre_PostProcDisplacemenysAndStarinAndElasticLinearStressOnRefMesh( Interface& _moab,double _lambda,double _mu, double _E_p,double _E_z, double _nu_p, double _nu_pz, double _G_zp):
    PostProcDisplacemenysAndStarinAndElasticLinearStressOnRefMesh(_moab,_lambda,_mu),E_p(_E_p),E_z(_E_z),nu_p(_nu_p),nu_pz(_nu_pz),G_zp(_G_zp) {

        double def_VAL2[3] = {0,0,0};
        rval = moab_post_proc.tag_get_handle("FIBRE_DIRECTION",3,MB_TYPE_DOUBLE,th_fibre_orientation,MB_TAG_CREAT|MB_TAG_SPARSE,&def_VAL2); CHKERR_THROW(rval);
        
        myfile.open("transIso_Mat_Blessed_File.txt");
        myfile<<"Blessed File for Testing Transversely Isotropic Material Implimentation on moFEM "<<endl;
        myfile<<"-------------------------------------------------------------------------------- "<<endl<<endl;
        myfile<<"Analytical Solution was computed using and equivalent rotated Compliance Stiffness Matrix multiplied by the applied stress tensor (0,0,1,0,0,0)"<<endl<<endl;
        myfile<<"    \t\t\tXX \tYY \tZZ \tXY \tXZ \tYZ"<<endl;
        
    };
    
    ~TranIsotropic_Fibre_PostProcDisplacemenysAndStarinAndElasticLinearStressOnRefMesh() {
        myfile.close();
    }
    
    PetscErrorCode operator()() {
        PetscFunctionBegin;
        
        ierr = do_operator(); CHKERRQ(ierr);
        ierr = OpStudentStart_TET(g_NTET); CHKERRQ(ierr);
        
        EntityHandle fe_handle = fe_ptr->get_ent();
        
        ///Get Stiffness Matrix
        ublas::symmetric_matrix<FieldData,ublas::upper> StiffnessMatrix;
        StiffnessMatrix.resize(6);
        StiffnessMatrix.clear();
        TransverseIsotropicStiffnessMatrix TranIsoMat(nu_p, nu_pz, E_p, E_z, G_zp);
        StiffnessMatrix=TranIsoMat.StiffnessMatrix;
//        IsotropicStiffnessMatrix IsoMat(lambda, mu);
//        StiffnessMatrix=IsoMat.StiffnessMatrix;
        
        ///Get Compliance Matrix
        ublas::symmetric_matrix<FieldData,ublas::upper> ComplianceMatrix;
        ComplianceMatrix.resize(6);
        ComplianceMatrix.clear();
        TransverseIsotropicComplianceMatrix TranIsoComMat(nu_p, nu_pz, E_p, E_z, G_zp);
        ComplianceMatrix=TranIsoComMat.ComplianceMatrix;
        
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
        ublas::matrix<double> DummyMatrixCom,DummyMatrixCom2; 
        DummyMatrix = ublas::zero_matrix<FieldData>(6,6); 
        DummyMatrixCom = ublas::zero_matrix<FieldData>(6,6); 
        DummyMatrix = StiffnessMatrix;
        DummyMatrixCom = ComplianceMatrix;
        
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
            
            StressTransformation invStressRotMat(&AxVector[3*aa], negAxAngle[aa]);
            StrainTransformation StrainRotMat(&AxVector[3*aa], AxAngle[aa]);
            
            ublas::matrix<double> TrpMatrixInvStress; 
            TrpMatrixInvStress = ublas::zero_matrix<FieldData>(6,6);
            TrpMatrixInvStress=invStressRotMat.StressRotMat;
            
            ublas::matrix<double> TrpMatrixStrain; 
            TrpMatrixStrain = ublas::zero_matrix<FieldData>(6,6); 
            TrpMatrixStrain=StrainRotMat.StrainRotMat;
            
            DummyMatrixCom2 = ublas::zero_matrix<FieldData>(6,6); 
            ublas::matrix< FieldData > dummyComA = prod( DummyMatrixCom , TrpMatrixInvStress );
            DummyMatrixCom2 = prod(TrpMatrixStrain,dummyComA);
            DummyMatrixCom = ublas::zero_matrix<FieldData>(6,6); 
            DummyMatrixCom = DummyMatrixCom2;
            
        }
        
        D.resize(6,6);
        D.clear();
        D = DummyMatrix;
        
        //Analytical Solution
        ublas::vector<FieldData> AppliedStress(6);
        AppliedStress.clear();
        AppliedStress[2]=1.0;
//        cout<<D<<endl;
//        cout<<"Applied Stress: "<<AppliedStress<<endl;
        ublas::vector<FieldData> ResultStrain=prod(DummyMatrixCom,AppliedStress);
//        cout<<"Result Strain: "<<ResultStrain<<endl<<endl;
        myfile<<"Analytical Strain: \t"<<boost::format("%.3lf") % roundn( ResultStrain[0] )<<"\t"<<boost::format("%.3lf") % roundn( ResultStrain[1] )<<"\t"<<boost::format("%.3lf") % roundn( ResultStrain[2] )<<"\t"<<boost::format("%.3lf") % roundn( ResultStrain[3] )<<"\t"<<boost::format("%.3lf") % roundn( ResultStrain[4] )<<"\t"<<boost::format("%.3lf") % roundn( ResultStrain[5] )<<endl;
        myfile<<"Analytical Stress: \t"<<boost::format("%.3lf") % roundn( AppliedStress[0] )<<"\t"<<boost::format("%.3lf") % roundn( AppliedStress[1] )<<"\t"<<boost::format("%.3lf") % roundn( AppliedStress[2] )<<"\t"<<boost::format("%.3lf") % roundn( AppliedStress[3] )<<"\t"<<boost::format("%.3lf") % roundn( AppliedStress[4] )<<"\t"<<boost::format("%.3lf") % roundn( AppliedStress[5] )<<endl;
        myfile<<endl;
        
        int gg=0;
        vector< ublas::matrix< FieldData > > GradU_at_GaussPt;
        ierr = GetGaussDiffDataVector(field_name,GradU_at_GaussPt); CHKERRQ(ierr);
        vector< ublas::matrix< FieldData > >::iterator viit = GradU_at_GaussPt.begin();
        map<EntityHandle,EntityHandle>::iterator mit = node_map.begin();
        
        double tetStrain[6], tetStress[6];
        
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
            
            tetStrain[0]=Strain(0,0); tetStrain[1]=Strain(1,1); tetStrain[2]=Strain(2,2);
            tetStrain[3]=Strain(0,1); tetStrain[4]=Strain(0,2); tetStrain[5]=Strain(1,2);
            tetStress[0]=Stress(0,0); tetStress[1]=Stress(1,1); tetStress[2]=Stress(2,2);
            tetStress[3]=Stress(0,1); tetStress[4]=Stress(0,2); tetStress[5]=Stress(1,2);
        }
        
        myfile<<"Numerical Strain: "<<"\t"<<boost::format("%.3lf") % roundn( tetStrain[0] )<<"\t"<<boost::format("%.3lf") % roundn( tetStrain[1] )<<"\t"<<boost::format("%.3lf") % roundn( tetStrain[2] )<<"\t"<<boost::format("%.3lf") % roundn( tetStrain[3] )<<"\t"<<boost::format("%.3lf") % roundn( tetStrain[4] )<<"\t"<<boost::format("%.3lf") % roundn( tetStrain[5] )<<endl;
        myfile<<"Numerical Stress: "<<"\t"<<boost::format("%.3lf") % roundn( tetStress[0] )<<"\t"<<boost::format("%.3lf") % roundn( tetStress[1] )<<"\t"<<boost::format("%.3lf") % roundn( tetStress[2] )<<"\t"<<boost::format("%.3lf") % roundn( tetStress[3] )<<"\t"<<boost::format("%.3lf") % roundn( tetStress[4] )<<"\t"<<boost::format("%.3lf") % roundn( tetStress[5] )<<endl;
        myfile<<endl;
        
        ierr = OpStudentEnd(); CHKERRQ(ierr);
        PetscFunctionReturn(0); 
    }

};

ErrorCode rval;
PetscErrorCode ierr;

static char help[] = "...\n\n";

int main(int argc, char *argv[]) {

  PetscInitialize(&argc,&argv,(char *)0,help);

  Core mb_instance;
  Interface& moab = mb_instance;
  int rank;
  MPI_Comm_rank(PETSC_COMM_WORLD,&rank);

  //Reade parameters from line command
  PetscBool flg = PETSC_TRUE;
  char mesh_file_name[255];
  ierr = PetscOptionsGetString(PETSC_NULL,"-my_file",mesh_file_name,255,&flg); CHKERRQ(ierr);
  if(flg != PETSC_TRUE) {
    SETERRQ(PETSC_COMM_SELF,1,"*** ERROR -my_file (MESH FILE NEEDED)");
  }

  //Read mesh to MOAB
  const char *option;
  option = "";//"PARALLEL=BCAST;";//;DEBUG_IO";
  rval = moab.load_file(mesh_file_name, 0, option); CHKERR_PETSC(rval); 
  ParallelComm* pcomm = ParallelComm::get_pcomm(&moab,MYPCOMM_INDEX);
  if(pcomm == NULL) pcomm =  new ParallelComm(&moab,PETSC_COMM_WORLD);

  //We need that for code profiling
  PetscLogDouble t1,t2;
  PetscLogDouble v1,v2;
  ierr = PetscGetTime(&v1); CHKERRQ(ierr);
  ierr = PetscGetCPUTime(&t1); CHKERRQ(ierr);

  //Create MoFEM (Joseph) database
  FieldCore core(moab);
  FieldInterface& mField = core;
    
    //Open mesh_file_name.txt for writing
//    ofstream myfile;
//    myfile.open ((string(mesh_file_name)+".txt").c_str());

  //ref meshset ref level 0
  ierr = mField.seed_ref_level_3D(0,0); CHKERRQ(ierr);

  // stl::bitset see for more details
  BitRefLevel bit_level0;
  bit_level0.set(0);
  EntityHandle meshset_level0;
  rval = moab.create_meshset(MESHSET_SET,meshset_level0); CHKERR_PETSC(rval);
  ierr = mField.seed_ref_level_3D(0,bit_level0); CHKERRQ(ierr);
  ierr = mField.refine_get_ents(bit_level0,BitRefLevel().set(),meshset_level0); CHKERRQ(ierr);

    
  for(_IT_CUBITMESHSETS_FOR_LOOP_(mField,cubit_it)) {
      EntityHandle cubit_meshset = cubit_it->meshset; 
      ierr = mField.refine_get_childern(cubit_meshset,bit_level0,cubit_meshset,MBTRI,true); CHKERRQ(ierr);
      ierr = mField.refine_get_childern(cubit_meshset,bit_level0,cubit_meshset,MBVERTEX,true); CHKERRQ(ierr);
      ierr = mField.refine_get_childern(cubit_meshset,bit_level0,cubit_meshset,MBEDGE,true); CHKERRQ(ierr);
  }

  /***/
  //Define problem

  //Fields
  ierr = mField.add_field("DISPLACEMENT",H1,3); CHKERRQ(ierr);

  //FE
  ierr = mField.add_finite_element("TRAN_ISOTROPIC_ELASTIC"); CHKERRQ(ierr);

  //Define rows/cols and element data
  ierr = mField.modify_finite_element_add_field_row("TRAN_ISOTROPIC_ELASTIC","DISPLACEMENT"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_col("TRAN_ISOTROPIC_ELASTIC","DISPLACEMENT"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_data("TRAN_ISOTROPIC_ELASTIC","DISPLACEMENT"); CHKERRQ(ierr);
    
  //define problems
  ierr = mField.add_problem("ELASTIC_MECHANICS"); CHKERRQ(ierr);

  //set finite elements for problem
  ierr = mField.modify_problem_add_finite_element("ELASTIC_MECHANICS","TRAN_ISOTROPIC_ELASTIC"); CHKERRQ(ierr);

  //set refinment level for problem
  ierr = mField.modify_problem_ref_level_add_bit("ELASTIC_MECHANICS",bit_level0); CHKERRQ(ierr);

  /***/
  //Declare problem

  //add entitities (by tets) to the field
  ierr = mField.add_ents_to_field_by_TETs(0,"DISPLACEMENT"); CHKERRQ(ierr);

  //add finite elements entities
    ierr = mField.add_ents_to_finite_element_EntType_by_bit_ref(bit_level0,"TRAN_ISOTROPIC_ELASTIC",MBTET); CHKERRQ(ierr);

  //set app. order
  //see Hierarchic Finite Element Bases on Unstructured Tetrahedral Meshes (Mark Ainsworth & Joe Coyle)
  ierr = mField.set_field_order(0,MBTET,"DISPLACEMENT",1); CHKERRQ(ierr);
  ierr = mField.set_field_order(0,MBTRI,"DISPLACEMENT",1); CHKERRQ(ierr);
  ierr = mField.set_field_order(0,MBEDGE,"DISPLACEMENT",1); CHKERRQ(ierr);
  ierr = mField.set_field_order(0,MBVERTEX,"DISPLACEMENT",1); CHKERRQ(ierr);
    
  /****/
  //build database

  //build field
  ierr = mField.build_fields(); CHKERRQ(ierr);

  //build finite elemnts
  ierr = mField.build_finite_elements(); CHKERRQ(ierr);

  //build adjacencies
  ierr = mField.build_adjacencies(bit_level0); CHKERRQ(ierr);

  //build problem
  ierr = mField.build_problems(); CHKERRQ(ierr);

  /****/
  //mesh partitioning 

  //partition
  ierr = mField.partition_problem("ELASTIC_MECHANICS"); CHKERRQ(ierr);
  ierr = mField.partition_finite_elements("ELASTIC_MECHANICS"); CHKERRQ(ierr);
  //what are ghost nodes, see Petsc Manual
  ierr = mField.partition_ghost_dofs("ELASTIC_MECHANICS"); CHKERRQ(ierr);

  //create matrices
    Vec D,F;
    ierr = mField.VecCreateGhost("ELASTIC_MECHANICS",Row,&F); CHKERRQ(ierr);
    ierr = mField.VecCreateGhost("ELASTIC_MECHANICS",Col,&D); CHKERRQ(ierr);
    Mat Aij;
    ierr = mField.MatCreateMPIAIJWithArrays("ELASTIC_MECHANICS",&Aij); CHKERRQ(ierr);
    
    CubitDisplacementDirihletBC myDirihletBC(mField,"ELASTIC_MECHANICS","DISPLACEMENT");
    ierr = myDirihletBC.Init(); CHKERRQ(ierr);
    
//    vector<double> attributes;
    double YoungModulusP;
    double PoissonRatioP;
    double YoungModulusZ;
    double PoissonRatioPZ;
    double ShearModulusZP;

    for(_IT_CUBITMESHSETS_BY_SET_TYPE_FOR_LOOP_(mField,BlockSet,it))
    {
        cout << endl << *it << endl;
        
        //Get block name
        string name = it->get_Cubit_name();
        
        if (name.compare(0,12,"MAT_TRANSISO") == 0)
        {
            Mat_TransIso mydata;
            ierr = it->get_attribute_data_structure(mydata); CHKERRQ(ierr);
            cout << mydata;
            YoungModulusP=mydata.data.Youngp;
            YoungModulusZ=mydata.data.Youngz;
            PoissonRatioP=mydata.data.Poissonp;
            PoissonRatioPZ=mydata.data.Poissonpz;
            ShearModulusZP=mydata.data.Shearzp;
        }
        
        else SETERRQ(PETSC_COMM_SELF,1,"Error: Unrecognizable Material type");
    }
    
  //Assemble F and Aij
    
  TranIsotropicElasticFEMethod MyTIsotFE(mField,&myDirihletBC,Aij,D,F,LAMBDA(YoungModulusP,PoissonRatioP),MU(YoungModulusP,PoissonRatioP),YoungModulusP,YoungModulusZ,PoissonRatioP,PoissonRatioPZ,ShearModulusZP);

  ierr = VecZeroEntries(F); CHKERRQ(ierr);
  ierr = VecGhostUpdateBegin(F,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = VecGhostUpdateEnd(F,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = MatZeroEntries(Aij); CHKERRQ(ierr);

    ierr = mField.loop_finite_elements("ELASTIC_MECHANICS","TRAN_ISOTROPIC_ELASTIC",MyTIsotFE);  CHKERRQ(ierr);

  ierr = VecGhostUpdateBegin(F,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
  ierr = VecGhostUpdateEnd(F,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
  ierr = VecAssemblyBegin(F); CHKERRQ(ierr);
  ierr = VecAssemblyEnd(F); CHKERRQ(ierr);
  ierr = MatAssemblyBegin(Aij,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
  ierr = MatAssemblyEnd(Aij,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);

  PetscSynchronizedFlush(PETSC_COMM_WORLD);

  //Matrix View
  //MatView(Aij,PETSC_VIEWER_DRAW_WORLD);//PETSC_VIEWER_STDOUT_WORLD);
  //std::string wait;
  //std::cin >> wait;

  //Solver
  KSP solver;
  ierr = KSPCreate(PETSC_COMM_WORLD,&solver); CHKERRQ(ierr);
  ierr = KSPSetOperators(solver,Aij,Aij,SAME_NONZERO_PATTERN); CHKERRQ(ierr);
  ierr = KSPSetFromOptions(solver); CHKERRQ(ierr);
  ierr = KSPSetUp(solver); CHKERRQ(ierr);

  ierr = KSPSolve(solver,F,D); CHKERRQ(ierr);
  ierr = VecGhostUpdateBegin(D,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = VecGhostUpdateEnd(D,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);

  //Save data on mesh
  ierr = mField.set_global_VecCreateGhost("ELASTIC_MECHANICS",Row,D,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
  //ierr = VecView(F,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);

  PostProcVertexMethod ent_method(moab);
  ierr = mField.loop_dofs("ELASTIC_MECHANICS","DISPLACEMENT",Row,ent_method); CHKERRQ(ierr);

  if(pcomm->rank()==0) {
    EntityHandle out_meshset;
    rval = moab.create_meshset(MESHSET_SET,out_meshset); CHKERR_PETSC(rval);
    ierr = mField.problem_get_FE("ELASTIC_MECHANICS","TRAN_ISOTROPIC_ELASTIC",out_meshset); CHKERRQ(ierr);
    rval = moab.write_file("out.vtk","VTK","",&out_meshset,1); CHKERR_PETSC(rval);
    rval = moab.delete_entities(&out_meshset,1); CHKERR_PETSC(rval);
  }

    TranIsotropic_Fibre_PostProcDisplacemenysAndStarinAndElasticLinearStressOnRefMesh fe_fibre_post_proc_method( moab, LAMBDA(YoungModulusP,PoissonRatioP),MU(YoungModulusP,PoissonRatioP),YoungModulusP,YoungModulusZ,PoissonRatioP,PoissonRatioPZ,ShearModulusZP);
    
  ierr = mField.loop_finite_elements("ELASTIC_MECHANICS","TRAN_ISOTROPIC_ELASTIC",fe_fibre_post_proc_method);  CHKERRQ(ierr);

  PetscSynchronizedFlush(PETSC_COMM_WORLD);
  if(pcomm->rank()==0) {
      rval = fe_fibre_post_proc_method.moab_post_proc.write_file("out_post_proc.vtk","VTK",""); CHKERR_PETSC(rval);
  }

  //detroy matrices
  ierr = VecDestroy(&F); CHKERRQ(ierr);
  ierr = VecDestroy(&D); CHKERRQ(ierr);
  ierr = MatDestroy(&Aij); CHKERRQ(ierr);
  ierr = KSPDestroy(&solver); CHKERRQ(ierr);


  ierr = PetscGetTime(&v2);CHKERRQ(ierr);
  ierr = PetscGetCPUTime(&t2);CHKERRQ(ierr);

  PetscSynchronizedPrintf(PETSC_COMM_WORLD,"Total Rank %d Time = %f CPU Time = %f\n",pcomm->rank(),v2-v1,t2-t1);
  PetscSynchronizedFlush(PETSC_COMM_WORLD);

  PetscFinalize();

}

