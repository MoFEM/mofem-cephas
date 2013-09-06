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

#include "moabField.hpp"
#include "moabField_Core.hpp"
#include "moabFEMethod_UpLevelStudent.hpp"
#include "cholesky.hpp"
#include <petscksp.h>

//#include <boost/numeric/ublas/matrix.hpp>
//#include <boost/numeric/ublas/lu.hpp>
//#include <boost/numeric/ublas/io.hpp>

#include "ElasticFEMethod.hpp"
#include "PostProcVertexMethod.hpp"
#include "PostProcDisplacementAndStrainOnRefindedMesh.hpp"

using namespace MoFEM;

//using namespace std;

///* Matrix inversion routine.
// Uses lu_factorize and lu_substitute in uBLAS to invert a matrix */
//template<class T>
//bool InvertMatrix(const ublas::matrix<T>& input, ublas::matrix<T>& inverse)
//{
//	typedef ublas::permutation_matrix<std::size_t> pmatrix;
//    
//	// create a working copy of the input
//    ublas::matrix<T> A(input);
//    
//	// create a permutation matrix for the LU-factorization
//	pmatrix pm(A.size1());
//    
//	// perform LU-factorization
//	int res = lu_factorize(A, pm);
//	if (res != 0)
//		return false;
//    
//	// create identity matrix of "inverse"
//    inverse.assign(ublas::identity_matrix<T> (A.size1()));
//    
//	// backsubstitute to get the inverse
//	lu_substitute(A, pm, inverse);
//    
//	return true;
//}

struct TranIsotropicElasticFEMethod: public ElasticFEMethod {
    
    double delta, nu_p, nu_pz, E_p, E_z, G_zp, AngleX, AngleY, AngleZ;
            
    TranIsotropicElasticFEMethod(
                             Interface& _moab,BaseDirihletBC *_dirihlet_ptr,Mat &_Aij,Vec& _F,
                             double _lambda,double _mu,double _E_p,double _E_z,double _delta,double _nu_p,double _nu_pz, double _G_zp, double _AngleX, double _AngleY, double _AngleZ, Range &_SideSet1,Range &_SideSet2): 
    ElasticFEMethod(_moab,_dirihlet_ptr,_Aij,_F,_lambda,_mu,_SideSet1,_SideSet2), E_p(_E_p), E_z(_E_z), delta(_delta), nu_p(_nu_p), nu_pz(_nu_pz), G_zp(_G_zp), AngleX(_AngleX),AngleY(_AngleY),AngleZ(_AngleZ) {};
    
    PetscErrorCode NeumannBC() {
        PetscFunctionBegin;
        ublas::vector<FieldData,ublas::bounded_array<double,3> > traction(3);
        //Set Direction of Traction On SideSet2
        traction[0] = 0; //X
        traction[1] = 1; //Y 
        traction[2] = 0; //Z
        //ElasticFEMethod::NeumannBC(...) function calulating external forces (see file ElasticFEMethod.hpp)
        ierr = ElasticFEMethod::NeumannBC(F,traction,SideSet2); CHKERRQ(ierr);
        PetscFunctionReturn(0);
    }
    
    PetscErrorCode preProcess() {
        PetscFunctionBegin;
        
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"Start Assembly\n",pcomm->rank(),v2-v1,t2-t1);
        
        ierr = PetscGetTime(&v1); CHKERRQ(ierr);
        ierr = PetscGetCPUTime(&t1); CHKERRQ(ierr);
        g_NTET.resize(4*45);
        ShapeMBTET(&g_NTET[0],G_TET_X45,G_TET_Y45,G_TET_Z45,45);
        g_NTRI.resize(3*13);
        ShapeMBTRI(&g_NTRI[0],G_TRI_X13,G_TRI_Y13,13); 
                        
//        double c1 = cos (0.5*AngleX);
//        double c2 = cos (0.5*AngleY);
//        double c3 = cos (0.5*AngleZ);
//        double s1 = sin (0.5*AngleX);
//        double s2 = sin (0.5*AngleY);
//        double s3 = sin (0.5*AngleZ);
//        
//        double AxAngle = 2*acos(c1*c2*c3 - s1*s2*s3);
//        double norm_AxVector;
////        double AxVector[3]={s1*s2*c3 + c1*c2*s3, s1*c2*c3 + c1*s2*s3, c1*s2*c3 + s1*c2*s3};
//        double AxVector[3]={s1*c2*c3 + c1*s2*s3 , c1*s2*c3 - s1*c2*s3 , s1*s2*c3 + c1*c2*s3};
//        norm_AxVector = sqrt(pow(AxVector[0],2) + pow(AxVector[1],2) + pow(AxVector[2],2));
//        
////        printf("AxVector: %f \t %f \t %f\n",AxVector[0]/norm_AxVector,AxVector[1]/norm_AxVector,AxVector[2]/norm_AxVector);
////        printf("AxAngle: %f\n",180*AxAngle/M_PI);
//        
////        AxVector[0]=1;AxVector[1]=0;AxVector[2]=0;
//        norm_AxVector = sqrt(pow(AxVector[0],2)+pow(AxVector[1],2)+pow(AxVector[2],2));
////        printf("norm_AxVector: %f\n",norm_AxVector);
////        AxAngle=0.5*M_PI;
//
//        ublas::matrix<double> Rotational_Matrix;
//        Rotational_Matrix = ublas::zero_matrix<FieldData>(3,3);
////        Rotational_Matrix.resize(3,3);
////        Rotational_Matrix.clear();
//        
//        if (norm_AxVector==0) {
//            for(int rr = 0;rr<3;rr++) {
//                Rotational_Matrix(rr,rr) = 1;
//            }
//        }
//        else{
//        //Axis-Angle Rotational Matrix
//        Rotational_Matrix(0,0) = 1-((1-cos(AxAngle))*(pow(AxVector[1],2)+pow(AxVector[2],2))/pow(norm_AxVector,2)); 
//        Rotational_Matrix(1,1) = 1-((1-cos(AxAngle))*(pow(AxVector[0],2)+pow(AxVector[2],2))/pow(norm_AxVector,2)); 
//        Rotational_Matrix(2,2) = 1-((1-cos(AxAngle))*(pow(AxVector[0],2)+pow(AxVector[1],2))/pow(norm_AxVector,2)); 
//        
//        Rotational_Matrix(0,1) = ((1-cos(AxAngle))*AxVector[0]*AxVector[1]-norm_AxVector*AxVector[2]*sin(AxAngle))/pow(norm_AxVector,2);
//        Rotational_Matrix(1,0) = ((1-cos(AxAngle))*AxVector[0]*AxVector[1]+norm_AxVector*AxVector[2]*sin(AxAngle))/pow(norm_AxVector,2);
//
//        Rotational_Matrix(0,2) = ((1-cos(AxAngle))*AxVector[0]*AxVector[2]+norm_AxVector*AxVector[1]*sin(AxAngle))/pow(norm_AxVector,2);
//        Rotational_Matrix(2,0) = ((1-cos(AxAngle))*AxVector[0]*AxVector[2]-norm_AxVector*AxVector[1]*sin(AxAngle))/pow(norm_AxVector,2);
//
//        Rotational_Matrix(1,2) = ((1-cos(AxAngle))*AxVector[1]*AxVector[2]-norm_AxVector*AxVector[0]*sin(AxAngle))/pow(norm_AxVector,2);
//        Rotational_Matrix(2,1) = ((1-cos(AxAngle))*AxVector[1]*AxVector[2]+norm_AxVector*AxVector[0]*sin(AxAngle))/pow(norm_AxVector,2);
//        
//        //Eulerian Rotational Matrix//
//        //Rotation about x-axis
////        Rotational_Matrix(0,0)=1;
////        Rotational_Matrix(1,1)=Rotational_Matrix(2,2)=cos(AxAngle);
////        Rotational_Matrix(1,2)=-sin(AxAngle);
////        Rotational_Matrix(2,1)=sin(AxAngle);
//        
//        //Rotation about y-axis
////        Rotational_Matrix(0,0)=Rotational_Matrix(2,2)=cos(AxAngle);
////        Rotational_Matrix(1,1)=1;
////        Rotational_Matrix(0,2)=sin(AxAngle);
////        Rotational_Matrix(2,0)=-sin(AxAngle);
//        
//        //Rotation about z-axis
////        Rotational_Matrix(0,0)=Rotational_Matrix(1,1)=cos(AxAngle);
////        Rotational_Matrix(2,2)=1;
////        Rotational_Matrix(0,1)=-sin(AxAngle);
////        Rotational_Matrix(1,0)=sin(AxAngle);
//        }
//        
//        //Check Rotational Matrix
////        ublas::matrix< FieldData > rTr = prod( trans(Rotational_Matrix) , Rotational_Matrix );
////
////        for (int rr=0; rr<3; rr++) {
////            for(int cc= 0;cc<3;cc++) {
////                printf("%f \t",rTr(rr,cc));
////            }printf("\n");}
////        printf("\n");
//        
//        ublas::matrix<double> TrpMatrix;
//        TrpMatrix = ublas::zero_matrix<FieldData>(6,6);
//                    
//        TrpMatrix(0, 0) = Rotational_Matrix(0,0) * Rotational_Matrix(0,0);
//        TrpMatrix(0, 1) = Rotational_Matrix(1,0) * Rotational_Matrix(1,0);
//        TrpMatrix(0, 2) = Rotational_Matrix(2,0) * Rotational_Matrix(2,0);
//        TrpMatrix(0, 3) = Rotational_Matrix(1,0) * Rotational_Matrix(2,0);
//        TrpMatrix(0, 4) = Rotational_Matrix(0,0) * Rotational_Matrix(2,0);
//        TrpMatrix(0, 5) = Rotational_Matrix(0,0) * Rotational_Matrix(1,0);
//         
//        TrpMatrix(1, 0) = Rotational_Matrix(0,1) * Rotational_Matrix(0,1);
//        TrpMatrix(1, 1) = Rotational_Matrix(1,1) * Rotational_Matrix(1,1);
//        TrpMatrix(1, 2) = Rotational_Matrix(2,1) * Rotational_Matrix(2,1);
//        TrpMatrix(1, 3) = Rotational_Matrix(1,1) * Rotational_Matrix(2,1);
//        TrpMatrix(1, 4) = Rotational_Matrix(0,1) * Rotational_Matrix(2,1);
//        TrpMatrix(1, 5) = Rotational_Matrix(0,1) * Rotational_Matrix(1,1);
//        
//        TrpMatrix(2, 0) = Rotational_Matrix(0,2) * Rotational_Matrix(0,2);
//        TrpMatrix(2, 1) = Rotational_Matrix(1,2) * Rotational_Matrix(1,2);
//        TrpMatrix(2, 2) = Rotational_Matrix(2,2) * Rotational_Matrix(2,2);
//        TrpMatrix(2, 3) = Rotational_Matrix(1,2) * Rotational_Matrix(2,2);
//        TrpMatrix(2, 4) = Rotational_Matrix(0,2) * Rotational_Matrix(2,2);
//        TrpMatrix(2, 5) = Rotational_Matrix(0,2) * Rotational_Matrix(1,2);
//        
//        TrpMatrix(3, 0) = 2.0 * Rotational_Matrix(0,1) * Rotational_Matrix(0,2);
//        TrpMatrix(3, 1) = 2.0 * Rotational_Matrix(1,1) * Rotational_Matrix(1,2);
//        TrpMatrix(3, 2) = 2.0 * Rotational_Matrix(2,1) * Rotational_Matrix(2,2);
//        TrpMatrix(3, 3) = ( Rotational_Matrix(1,1) * Rotational_Matrix(2,2) + Rotational_Matrix(2,1) * Rotational_Matrix(1,2) );
//        TrpMatrix(3, 4) = ( Rotational_Matrix(0,1) * Rotational_Matrix(2,2) + Rotational_Matrix(2,1) * Rotational_Matrix(0,2) );
//        TrpMatrix(3, 5) = ( Rotational_Matrix(0,1) * Rotational_Matrix(1,2) + Rotational_Matrix(1,1) * Rotational_Matrix(0,2) );
//    
//        TrpMatrix(4, 0) = 2.0 * Rotational_Matrix(0,0) * Rotational_Matrix(0,2);
//        TrpMatrix(4, 1) = 2.0 * Rotational_Matrix(1,0) * Rotational_Matrix(1,2);
//        TrpMatrix(4, 2) = 2.0 * Rotational_Matrix(2,0) * Rotational_Matrix(2,2);
//        TrpMatrix(4, 3) = ( Rotational_Matrix(1,0) * Rotational_Matrix(2,2) + Rotational_Matrix(2,0) * Rotational_Matrix(1,2) );
//        TrpMatrix(4, 4) = ( Rotational_Matrix(0,0) * Rotational_Matrix(2,2) + Rotational_Matrix(2,0) * Rotational_Matrix(0,2) );
//        TrpMatrix(4, 5) = ( Rotational_Matrix(0,0) * Rotational_Matrix(1,2) + Rotational_Matrix(1,0) * Rotational_Matrix(0,2) );
//        
//        TrpMatrix(5, 0) = 2.0 * Rotational_Matrix(0,0) * Rotational_Matrix(0,1);
//        TrpMatrix(5, 1) = 2.0 * Rotational_Matrix(1,0) * Rotational_Matrix(1,1);
//        TrpMatrix(5, 2) = 2.0 * Rotational_Matrix(2,0) * Rotational_Matrix(2,1);
//        TrpMatrix(5, 3) = ( Rotational_Matrix(1,0) * Rotational_Matrix(2,1) + Rotational_Matrix(2,0) * Rotational_Matrix(1,1) );
//        TrpMatrix(5, 4) = ( Rotational_Matrix(0,0) * Rotational_Matrix(2,1) + Rotational_Matrix(2,0) * Rotational_Matrix(0,1) );
//        TrpMatrix(5, 5) = ( Rotational_Matrix(0,0) * Rotational_Matrix(1,1) + Rotational_Matrix(1,0) * Rotational_Matrix(0,1) );
        
        //Transverse-Isotropic Stiffness Matrix (about y-axis)
//        ublas::matrix<double> StiffnessMatrix;
//        StiffnessMatrix = ublas::zero_matrix<FieldData>(6,6);
//        double nu_zp=(nu_pz*E_z)/E_p;
//
//        StiffnessMatrix(0,0)=StiffnessMatrix(1,1)=(1-nu_pz*nu_zp)/(E_p*E_z*delta);
//        StiffnessMatrix(2,2)=(1-nu_p*nu_p)/(E_p*E_p*delta);
//        
//        StiffnessMatrix(0,1)=StiffnessMatrix(1,0)=(nu_p+nu_zp*nu_pz)/(E_p*E_z*delta);
//        StiffnessMatrix(0,2)=StiffnessMatrix(2,0)=StiffnessMatrix(1,2)=StiffnessMatrix(2,1)=(nu_zp+nu_p*nu_zp)/(E_p*E_z*delta);
//
//        StiffnessMatrix(3,3)=E_p/(1+nu_p);
//        StiffnessMatrix(5,5)=StiffnessMatrix(4,4)=2*G_zp;
        
        //Isotropic Stiffness Matrix
//        D_lambda.resize(6);
//        D_lambda.clear();
//        for(int rr = 0;rr<3;rr++) {
//            for(int cc = 0;cc<3;cc++) {
//                D_lambda(rr,cc) = 1;
//            }
//        }
//        D_mu.resize(6);
//        D_mu.clear();
//        for(int rr = 0;rr<6;rr++) {
//            D_mu(rr,rr) = rr<3 ? 2 : 1;
//        }
//        StiffnessMatrix = lambda*D_lambda + mu*D_mu;
        
//        D.resize(6);
//        D.clear();
//        ublas::matrix< FieldData > dummy2 = prod( StiffnessMatrix , TrpMatrix );
//        D=prod( trans(TrpMatrix) , dummy2 );
        
//        for (int rr=0; rr<6; rr++) {
//            for(int cc= 0;cc<6;cc++) {
//                printf("%f \t",D(rr,cc));
//            }printf("\n");}
//        printf("\n");

        ierr = VecDuplicate(F,&Diagonal); CHKERRQ(ierr);
        PetscFunctionReturn(0);        
    }
    
    PetscErrorCode operator()() {
        PetscFunctionBegin;
        ierr = OpStudentStart_TET(g_NTET); CHKERRQ(ierr);
        ierr = GetMatrices(); CHKERRQ(ierr);
        
        double c1 = cos (0.5*AngleX);
        double c2 = cos (0.5*AngleY);
        double c3 = cos (0.5*AngleZ);
        double s1 = sin (0.5*AngleX);
        double s2 = sin (0.5*AngleY);
        double s3 = sin (0.5*AngleZ);
        
        double AxAngle = 2*acos(c1*c2*c3 - s1*s2*s3);
        double norm_AxVector;
        //        double AxVector[3]={s1*s2*c3 + c1*c2*s3, s1*c2*c3 + c1*s2*s3, c1*s2*c3 + s1*c2*s3};
        double AxVector[3]={s1*c2*c3 + c1*s2*s3 , c1*s2*c3 - s1*c2*s3 , s1*s2*c3 + c1*c2*s3};
        norm_AxVector = sqrt(pow(AxVector[0],2) + pow(AxVector[1],2) + pow(AxVector[2],2));
        
        ublas::matrix<double> Rotational_Matrix;
        Rotational_Matrix = ublas::zero_matrix<FieldData>(3,3);
        
        if (norm_AxVector==0) {
            for(int rr = 0;rr<3;rr++) {
                Rotational_Matrix(rr,rr) = 1;
            }
        }
        else{
            //Axis-Angle Rotational Matrix
            Rotational_Matrix(0,0) = 1-((1-cos(AxAngle))*(pow(AxVector[1],2)+pow(AxVector[2],2))/pow(norm_AxVector,2)); 
            Rotational_Matrix(1,1) = 1-((1-cos(AxAngle))*(pow(AxVector[0],2)+pow(AxVector[2],2))/pow(norm_AxVector,2)); 
            Rotational_Matrix(2,2) = 1-((1-cos(AxAngle))*(pow(AxVector[0],2)+pow(AxVector[1],2))/pow(norm_AxVector,2)); 
            
            Rotational_Matrix(0,1) = ((1-cos(AxAngle))*AxVector[0]*AxVector[1]-norm_AxVector*AxVector[2]*sin(AxAngle))/pow(norm_AxVector,2);
            Rotational_Matrix(1,0) = ((1-cos(AxAngle))*AxVector[0]*AxVector[1]+norm_AxVector*AxVector[2]*sin(AxAngle))/pow(norm_AxVector,2);
            
            Rotational_Matrix(0,2) = ((1-cos(AxAngle))*AxVector[0]*AxVector[2]+norm_AxVector*AxVector[1]*sin(AxAngle))/pow(norm_AxVector,2);
            Rotational_Matrix(2,0) = ((1-cos(AxAngle))*AxVector[0]*AxVector[2]-norm_AxVector*AxVector[1]*sin(AxAngle))/pow(norm_AxVector,2);
            
            Rotational_Matrix(1,2) = ((1-cos(AxAngle))*AxVector[1]*AxVector[2]-norm_AxVector*AxVector[0]*sin(AxAngle))/pow(norm_AxVector,2);
            Rotational_Matrix(2,1) = ((1-cos(AxAngle))*AxVector[1]*AxVector[2]+norm_AxVector*AxVector[0]*sin(AxAngle))/pow(norm_AxVector,2);            
            
            //Eulerian Rotational Matrix//
            //Rotation about x-axis
//              Rotational_Matrix(0,0)=1;
//              Rotational_Matrix(1,1)=Rotational_Matrix(2,2)=cos(AxAngle);
//              Rotational_Matrix(1,2)=-sin(AxAngle);
//              Rotational_Matrix(2,1)=sin(AxAngle);
            
            //Rotation about y-axis
//              Rotational_Matrix(0,0)=Rotational_Matrix(2,2)=cos(AxAngle);
//              Rotational_Matrix(1,1)=1;
//              Rotational_Matrix(0,2)=sin(AxAngle);
//              Rotational_Matrix(2,0)=-sin(AxAngle);
            
            //Rotation about z-axis
//              Rotational_Matrix(0,0)=Rotational_Matrix(1,1)=cos(AxAngle);
//              Rotational_Matrix(2,2)=1;
//              Rotational_Matrix(0,1)=-sin(AxAngle);
//              Rotational_Matrix(1,0)=sin(AxAngle);
        }
        
        //Check Rotational Matrix
        //        ublas::matrix< FieldData > rTr = prod( trans(Rotational_Matrix) , Rotational_Matrix );
        //
        //        for (int rr=0; rr<3; rr++) {
        //            for(int cc= 0;cc<3;cc++) {
        //                printf("%f \t",rTr(rr,cc));
        //            }printf("\n");}
        //        printf("\n");
        
        ublas::matrix<double> TrpMatrix;
        TrpMatrix = ublas::zero_matrix<FieldData>(6,6);
        
        TrpMatrix(0, 0) = Rotational_Matrix(0,0) * Rotational_Matrix(0,0);
        TrpMatrix(0, 1) = Rotational_Matrix(1,0) * Rotational_Matrix(1,0);
        TrpMatrix(0, 2) = Rotational_Matrix(2,0) * Rotational_Matrix(2,0);
        TrpMatrix(0, 3) = Rotational_Matrix(1,0) * Rotational_Matrix(2,0);
        TrpMatrix(0, 4) = Rotational_Matrix(0,0) * Rotational_Matrix(2,0);
        TrpMatrix(0, 5) = Rotational_Matrix(0,0) * Rotational_Matrix(1,0);
        
        TrpMatrix(1, 0) = Rotational_Matrix(0,1) * Rotational_Matrix(0,1);
        TrpMatrix(1, 1) = Rotational_Matrix(1,1) * Rotational_Matrix(1,1);
        TrpMatrix(1, 2) = Rotational_Matrix(2,1) * Rotational_Matrix(2,1);
        TrpMatrix(1, 3) = Rotational_Matrix(1,1) * Rotational_Matrix(2,1);
        TrpMatrix(1, 4) = Rotational_Matrix(0,1) * Rotational_Matrix(2,1);
        TrpMatrix(1, 5) = Rotational_Matrix(0,1) * Rotational_Matrix(1,1);
        
        TrpMatrix(2, 0) = Rotational_Matrix(0,2) * Rotational_Matrix(0,2);
        TrpMatrix(2, 1) = Rotational_Matrix(1,2) * Rotational_Matrix(1,2);
        TrpMatrix(2, 2) = Rotational_Matrix(2,2) * Rotational_Matrix(2,2);
        TrpMatrix(2, 3) = Rotational_Matrix(1,2) * Rotational_Matrix(2,2);
        TrpMatrix(2, 4) = Rotational_Matrix(0,2) * Rotational_Matrix(2,2);
        TrpMatrix(2, 5) = Rotational_Matrix(0,2) * Rotational_Matrix(1,2);
        
        TrpMatrix(3, 0) = 2.0 * Rotational_Matrix(0,1) * Rotational_Matrix(0,2);
        TrpMatrix(3, 1) = 2.0 * Rotational_Matrix(1,1) * Rotational_Matrix(1,2);
        TrpMatrix(3, 2) = 2.0 * Rotational_Matrix(2,1) * Rotational_Matrix(2,2);
        TrpMatrix(3, 3) = ( Rotational_Matrix(1,1) * Rotational_Matrix(2,2) + Rotational_Matrix(2,1) * Rotational_Matrix(1,2) );
        TrpMatrix(3, 4) = ( Rotational_Matrix(0,1) * Rotational_Matrix(2,2) + Rotational_Matrix(2,1) * Rotational_Matrix(0,2) );
        TrpMatrix(3, 5) = ( Rotational_Matrix(0,1) * Rotational_Matrix(1,2) + Rotational_Matrix(1,1) * Rotational_Matrix(0,2) );
        
        TrpMatrix(4, 0) = 2.0 * Rotational_Matrix(0,0) * Rotational_Matrix(0,2);
        TrpMatrix(4, 1) = 2.0 * Rotational_Matrix(1,0) * Rotational_Matrix(1,2);
        TrpMatrix(4, 2) = 2.0 * Rotational_Matrix(2,0) * Rotational_Matrix(2,2);
        TrpMatrix(4, 3) = ( Rotational_Matrix(1,0) * Rotational_Matrix(2,2) + Rotational_Matrix(2,0) * Rotational_Matrix(1,2) );
        TrpMatrix(4, 4) = ( Rotational_Matrix(0,0) * Rotational_Matrix(2,2) + Rotational_Matrix(2,0) * Rotational_Matrix(0,2) );
        TrpMatrix(4, 5) = ( Rotational_Matrix(0,0) * Rotational_Matrix(1,2) + Rotational_Matrix(1,0) * Rotational_Matrix(0,2) );
        
        TrpMatrix(5, 0) = 2.0 * Rotational_Matrix(0,0) * Rotational_Matrix(0,1);
        TrpMatrix(5, 1) = 2.0 * Rotational_Matrix(1,0) * Rotational_Matrix(1,1);
        TrpMatrix(5, 2) = 2.0 * Rotational_Matrix(2,0) * Rotational_Matrix(2,1);
        TrpMatrix(5, 3) = ( Rotational_Matrix(1,0) * Rotational_Matrix(2,1) + Rotational_Matrix(2,0) * Rotational_Matrix(1,1) );
        TrpMatrix(5, 4) = ( Rotational_Matrix(0,0) * Rotational_Matrix(2,1) + Rotational_Matrix(2,0) * Rotational_Matrix(0,1) );
        TrpMatrix(5, 5) = ( Rotational_Matrix(0,0) * Rotational_Matrix(1,1) + Rotational_Matrix(1,0) * Rotational_Matrix(0,1) );
        
        ublas::matrix<double> StiffnessMatrix;
        StiffnessMatrix = ublas::zero_matrix<FieldData>(6,6);
        double nu_zp=(nu_pz*E_z)/E_p;
        
        StiffnessMatrix(0,0)=StiffnessMatrix(1,1)=(1-nu_pz*nu_zp)/(E_p*E_z*delta);
        StiffnessMatrix(2,2)=(1-nu_p*nu_p)/(E_p*E_p*delta);
        
        StiffnessMatrix(0,1)=StiffnessMatrix(1,0)=(nu_p+nu_zp*nu_pz)/(E_p*E_z*delta);
        StiffnessMatrix(0,2)=StiffnessMatrix(2,0)=StiffnessMatrix(1,2)=StiffnessMatrix(2,1)=(nu_zp+nu_p*nu_zp)/(E_p*E_z*delta);
        
        StiffnessMatrix(3,3)=E_p/(1+nu_p);
        StiffnessMatrix(5,5)=StiffnessMatrix(4,4)=2*G_zp;
        
        D.resize(6);
        D.clear();
        ublas::matrix< FieldData > dummy2 = prod( StiffnessMatrix , TrpMatrix );
        D=prod( trans(TrpMatrix) , dummy2 );
        
        //Dirihlet Boundary Condition
        ierr = SetDirihletBC_to_ElementIndicies(); CHKERRQ(ierr);
        if(Diagonal!=PETSC_NULL) {
            if(DirihletBC.size()>0) {
                DirihletBCDiagVal.resize(DirihletBC.size());
                fill(DirihletBCDiagVal.begin(),DirihletBCDiagVal.end(),1);
                ierr = VecSetValues(Diagonal,DirihletBC.size(),&(DirihletBC[0]),&DirihletBCDiagVal[0],INSERT_VALUES); CHKERRQ(ierr);
            }
        }
        //Assembly Aij and F
        ierr = RhsAndLhs(); CHKERRQ(ierr);
        
        //Neumann Boundary Conditions
        ierr = NeumannBC(); CHKERRQ(ierr);
        
        ierr = OpStudentEnd(); CHKERRQ(ierr);
        PetscFunctionReturn(0); }
            
};

/*struct AxisOrthotropicElasticFEMethod: public ElasticFEMethod {
    
    double E1,E2,nu1,nu2,G12;
    
    TranIsotropicElasticFEMethod(
                                 Interface& _moab,BaseDirihletBC *_dirihlet_ptr,Mat &_Aij,Vec& _F,
                                 double _lambda,double _mu,double _E1,double _E2,double _nu1,double _nu2, double _G12, Range &_SideSet1,Range &_SideSet2): 
    ElasticFEMethod(_moab,_dirihlet_ptr,_Aij,_F,_lambda,_mu,_SideSet1,_SideSet2), E1(_E1), E2(_E2), nu1(_nu1), nu2(_nu2), G12(_G12) {};
    
    PetscErrorCode NeumannBC() {
        PetscFunctionBegin;
        
        ublas::vector<FieldData,ublas::bounded_array<double,3> > traction2(3);
        traction2[0] = 0;
        traction2[1] = +1;
        traction2[2] = 0;
        ierr = ElasticFEMethod::NeumannBC(F,traction2,SideSet2); CHKERRQ(ierr);
        
        PetscFunctionReturn(0);
    }
    
    PetscErrorCode preProcess() {
        PetscFunctionBegin;
        
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"Start Assembly\n",pcomm->rank(),v2-v1,t2-t1);
        
        ierr = PetscGetTime(&v1); CHKERRQ(ierr);
        ierr = PetscGetCPUTime(&t1); CHKERRQ(ierr);
        g_NTET.resize(4*45);
        ShapeMBTET(&g_NTET[0],G_TET_X45,G_TET_Y45,G_TET_Z45,45);
        g_NTRI.resize(3*13);
        ShapeMBTRI(&g_NTRI[0],G_TRI_X13,G_TRI_Y13,13); 
        
        double denominator = (1/(E1*E2*E2*G12))-((2*nu1*nu1)/(E1*E1*E2*G12))-((2*nu1*nu1*nu2)/(E1*E1*E2*G12))-((nu2*nu2)/(E1*E2*E2*G12));
        
        
        //Polar Coord Axisymmetric Orthotropic Material
//        double nu_zp=(nu_pz*E_z)/E_p;
//        D.resize(6);
//        D.clear();
//        D(0,0)=(1-nu_pz*nu_zp)/(E_p*E_z*delta);
//        D(1,1)=(1-nu_zp*nu_pz)/(E_z*E_p*delta);
//        D(2,2)=(1-nu_p*nu_p)/(E_p*E_p*delta);
//        
//        D(0,1)=D(1,0)=(nu_p+nu_zp*nu_pz)/(E_p*E_z*delta);
//        D(0,2)=D(2,0)=(nu_zp+nu_p*nu_zp)/(E_p*E_z*delta);
//        D(1,2)=D(2,1)=(nu_zp+nu_zp*nu_p)/(E_z*E_p*delta);
//        
//        D(3,3)=E_p/(1+nu_p);
//        D(5,5)=D(4,4)=2*G_zp;
        
        
        
        ierr = VecDuplicate(F,&Diagonal); CHKERRQ(ierr);
        
        PetscFunctionReturn(0);
    }
};
*/
struct TranIsotropicPostProcDisplacemenysAndStarinAndElasticLinearStressOnRefMesh: public PostProcDisplacemenysAndStarinAndElasticLinearStressOnRefMesh {
    
    double delta, nu_p, nu_pz, E_p, E_z, G_zp, AngleX, AngleY, AngleZ;
        
    TranIsotropicPostProcDisplacemenysAndStarinAndElasticLinearStressOnRefMesh( Interface& _moab,double _lambda,double _mu, double _E_p,double _E_z, double _delta, double _nu_p, double _nu_pz, double _G_zp,double _AngleX, double _AngleY, double _AngleZ):
    PostProcDisplacemenysAndStarinAndElasticLinearStressOnRefMesh(_moab,_lambda,_mu),E_p(_E_p),E_z(_E_z),delta(_delta),nu_p(_nu_p),nu_pz(_nu_pz),G_zp(_G_zp),AngleX(_AngleX),AngleY(_AngleY),AngleZ(_AngleZ) {
        
        
        double c1 = cos (0.5*AngleX);
        double c2 = cos (0.5*AngleY);
        double c3 = cos (0.5*AngleZ);
        double s1 = sin (0.5*AngleX);
        double s2 = sin (0.5*AngleY);
        double s3 = sin (0.5*AngleZ);
        
        double AxAngle = 2*acos(c1*c2*c3 - s1*s2*s3);
        double norm_AxVector;
        double AxVector[3]={s1*c2*c3 + c1*s2*s3 , c1*s2*c3 - s1*c2*s3 , s1*s2*c3 + c1*c2*s3};
        norm_AxVector = sqrt(pow(AxVector[0],2) + pow(AxVector[1],2) + pow(AxVector[2],2));
        
        ublas::matrix<double> Rotational_Matrix;
        Rotational_Matrix = ublas::zero_matrix<FieldData>(3,3);
            
        if (norm_AxVector==0) {
            for(int rr = 0;rr<3;rr++) {
                Rotational_Matrix(rr,rr) = 1;
            }
        }
        else{
            //Axis-Angle Rotational Matrix
            Rotational_Matrix(0,0) = 1-((1-cos(AxAngle))*(pow(AxVector[1],2)+pow(AxVector[2],2))/pow(norm_AxVector,2)); 
            Rotational_Matrix(1,1) = 1-((1-cos(AxAngle))*(pow(AxVector[0],2)+pow(AxVector[2],2))/pow(norm_AxVector,2)); 
            Rotational_Matrix(2,2) = 1-((1-cos(AxAngle))*(pow(AxVector[0],2)+pow(AxVector[1],2))/pow(norm_AxVector,2)); 
            
            Rotational_Matrix(0,1) = ((1-cos(AxAngle))*AxVector[0]*AxVector[1]-norm_AxVector*AxVector[2]*sin(AxAngle))/pow(norm_AxVector,2);
            Rotational_Matrix(1,0) = ((1-cos(AxAngle))*AxVector[0]*AxVector[1]+norm_AxVector*AxVector[2]*sin(AxAngle))/pow(norm_AxVector,2);
            
            Rotational_Matrix(0,2) = ((1-cos(AxAngle))*AxVector[0]*AxVector[2]+norm_AxVector*AxVector[1]*sin(AxAngle))/pow(norm_AxVector,2);
            Rotational_Matrix(2,0) = ((1-cos(AxAngle))*AxVector[0]*AxVector[2]-norm_AxVector*AxVector[1]*sin(AxAngle))/pow(norm_AxVector,2);
            
            Rotational_Matrix(1,2) = ((1-cos(AxAngle))*AxVector[1]*AxVector[2]-norm_AxVector*AxVector[0]*sin(AxAngle))/pow(norm_AxVector,2);
            Rotational_Matrix(2,1) = ((1-cos(AxAngle))*AxVector[1]*AxVector[2]+norm_AxVector*AxVector[0]*sin(AxAngle))/pow(norm_AxVector,2);
        
        }
        
        //Check Rotational Matrix
        //        ublas::matrix< FieldData > rTr = prod( trans(Rotational_Matrix) , Rotational_Matrix );
        //
        //        for (int rr=0; rr<3; rr++) {
        //            for(int cc= 0;cc<3;cc++) {
        //                printf("%f \t",rTr(rr,cc));
        //            }printf("\n");}
        //        printf("\n");
        
        ublas::matrix<double> TrpMatrix;
        TrpMatrix = ublas::zero_matrix<FieldData>(6,6);
        
        TrpMatrix(0, 0) = Rotational_Matrix(0,0) * Rotational_Matrix(0,0);
        TrpMatrix(0, 1) = Rotational_Matrix(1,0) * Rotational_Matrix(1,0);
        TrpMatrix(0, 2) = Rotational_Matrix(2,0) * Rotational_Matrix(2,0);
        TrpMatrix(0, 3) = Rotational_Matrix(1,0) * Rotational_Matrix(2,0);
        TrpMatrix(0, 4) = Rotational_Matrix(0,0) * Rotational_Matrix(2,0);
        TrpMatrix(0, 5) = Rotational_Matrix(0,0) * Rotational_Matrix(1,0);
        
        TrpMatrix(1, 0) = Rotational_Matrix(0,1) * Rotational_Matrix(0,1);
        TrpMatrix(1, 1) = Rotational_Matrix(1,1) * Rotational_Matrix(1,1);
        TrpMatrix(1, 2) = Rotational_Matrix(2,1) * Rotational_Matrix(2,1);
        TrpMatrix(1, 3) = Rotational_Matrix(1,1) * Rotational_Matrix(2,1);
        TrpMatrix(1, 4) = Rotational_Matrix(0,1) * Rotational_Matrix(2,1);
        TrpMatrix(1, 5) = Rotational_Matrix(0,1) * Rotational_Matrix(1,1);
        
        TrpMatrix(2, 0) = Rotational_Matrix(0,2) * Rotational_Matrix(0,2);
        TrpMatrix(2, 1) = Rotational_Matrix(1,2) * Rotational_Matrix(1,2);
        TrpMatrix(2, 2) = Rotational_Matrix(2,2) * Rotational_Matrix(2,2);
        TrpMatrix(2, 3) = Rotational_Matrix(1,2) * Rotational_Matrix(2,2);
        TrpMatrix(2, 4) = Rotational_Matrix(0,2) * Rotational_Matrix(2,2);
        TrpMatrix(2, 5) = Rotational_Matrix(0,2) * Rotational_Matrix(1,2);
        
        TrpMatrix(3, 0) = 2.0 * Rotational_Matrix(0,1) * Rotational_Matrix(0,2);
        TrpMatrix(3, 1) = 2.0 * Rotational_Matrix(1,1) * Rotational_Matrix(1,2);
        TrpMatrix(3, 2) = 2.0 * Rotational_Matrix(2,1) * Rotational_Matrix(2,2);
        TrpMatrix(3, 3) = ( Rotational_Matrix(1,1) * Rotational_Matrix(2,2) + Rotational_Matrix(2,1) * Rotational_Matrix(1,2) );
        TrpMatrix(3, 4) = ( Rotational_Matrix(0,1) * Rotational_Matrix(2,2) + Rotational_Matrix(2,1) * Rotational_Matrix(0,2) );
        TrpMatrix(3, 5) = ( Rotational_Matrix(0,1) * Rotational_Matrix(1,2) + Rotational_Matrix(1,1) * Rotational_Matrix(0,2) );
        
        TrpMatrix(4, 0) = 2.0 * Rotational_Matrix(0,0) * Rotational_Matrix(0,2);
        TrpMatrix(4, 1) = 2.0 * Rotational_Matrix(1,0) * Rotational_Matrix(1,2);
        TrpMatrix(4, 2) = 2.0 * Rotational_Matrix(2,0) * Rotational_Matrix(2,2);
        TrpMatrix(4, 3) = ( Rotational_Matrix(1,0) * Rotational_Matrix(2,2) + Rotational_Matrix(2,0) * Rotational_Matrix(1,2) );
        TrpMatrix(4, 4) = ( Rotational_Matrix(0,0) * Rotational_Matrix(2,2) + Rotational_Matrix(2,0) * Rotational_Matrix(0,2) );
        TrpMatrix(4, 5) = ( Rotational_Matrix(0,0) * Rotational_Matrix(1,2) + Rotational_Matrix(1,0) * Rotational_Matrix(0,2) );
        
        TrpMatrix(5, 0) = 2.0 * Rotational_Matrix(0,0) * Rotational_Matrix(0,1);
        TrpMatrix(5, 1) = 2.0 * Rotational_Matrix(1,0) * Rotational_Matrix(1,1);
        TrpMatrix(5, 2) = 2.0 * Rotational_Matrix(2,0) * Rotational_Matrix(2,1);
        TrpMatrix(5, 3) = ( Rotational_Matrix(1,0) * Rotational_Matrix(2,1) + Rotational_Matrix(2,0) * Rotational_Matrix(1,1) );
        TrpMatrix(5, 4) = ( Rotational_Matrix(0,0) * Rotational_Matrix(2,1) + Rotational_Matrix(2,0) * Rotational_Matrix(0,1) );
        TrpMatrix(5, 5) = ( Rotational_Matrix(0,0) * Rotational_Matrix(1,1) + Rotational_Matrix(1,0) * Rotational_Matrix(0,1) );
        
        //Transverse-Isotropic Stiffness Matrix (about y-axis)
        ublas::matrix<double> StiffnessMatrix;
        StiffnessMatrix = ublas::zero_matrix<FieldData>(6,6);
        double nu_zp=(nu_pz*E_z)/E_p;
        
        StiffnessMatrix(0,0)=StiffnessMatrix(1,1)=(1-nu_pz*nu_zp)/(E_p*E_z*delta);
        StiffnessMatrix(2,2)=(1-nu_p*nu_p)/(E_p*E_p*delta);
        
        StiffnessMatrix(0,1)=StiffnessMatrix(1,0)=(nu_p+nu_zp*nu_pz)/(E_p*E_z*delta);
        StiffnessMatrix(0,2)=StiffnessMatrix(2,0)=StiffnessMatrix(1,2)=StiffnessMatrix(2,1)=(nu_zp+nu_p*nu_zp)/(E_p*E_z*delta);
        
        StiffnessMatrix(3,3)=E_p/(1+nu_p);
        StiffnessMatrix(5,5)=StiffnessMatrix(4,4)=2*G_zp;
        
        //Isotropic Stiffness Matrix
        //        D_lambda.resize(6);
        //        D_lambda.clear();
        //        for(int rr = 0;rr<3;rr++) {
        //            for(int cc = 0;cc<3;cc++) {
        //                D_lambda(rr,cc) = 1;
        //            }
        //        }
        //        D_mu.resize(6);
        //        D_mu.clear();
        //        for(int rr = 0;rr<6;rr++) {
        //            D_mu(rr,rr) = rr<3 ? 2 : 1;
        //        }
        //        StiffnessMatrix = lambda*D_lambda + mu*D_mu;

        D.resize(6,6);
        D.clear();
        ublas::matrix< FieldData > dummy2 = prod( StiffnessMatrix , TrpMatrix );
        D=prod( trans(TrpMatrix) , dummy2 );
        
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
  moabField_Core core(moab);
  moabField& mField = core;

  //ref meshset ref level 0
  ierr = mField.seed_ref_level_3D(0,0); CHKERRQ(ierr);

  // stl::bitset see for more details
  BitRefLevel bit_level0;
  bit_level0.set(0);
  EntityHandle meshset_level0;
  rval = moab.create_meshset(MESHSET_SET,meshset_level0); CHKERR_PETSC(rval);
  ierr = mField.seed_ref_level_3D(0,bit_level0); CHKERRQ(ierr);
  ierr = mField.refine_get_ents(bit_level0,meshset_level0); CHKERRQ(ierr);

  /***/
  //Define problem

  //Fields
  ierr = mField.add_field("DISPLACEMENT",H1,3); CHKERRQ(ierr);

  //FE
  ierr = mField.add_finite_element("ELASTIC"); CHKERRQ(ierr);
  ierr = mField.add_finite_element("TRAN_ISOTROPIC_ELASTIC"); CHKERRQ(ierr);

  //Define rows/cols and element data
  ierr = mField.modify_finite_element_add_field_row("ELASTIC","DISPLACEMENT"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_col("ELASTIC","DISPLACEMENT"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_data("ELASTIC","DISPLACEMENT"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_row("TRAN_ISOTROPIC_ELASTIC","DISPLACEMENT"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_col("TRAN_ISOTROPIC_ELASTIC","DISPLACEMENT"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_data("TRAN_ISOTROPIC_ELASTIC","DISPLACEMENT"); CHKERRQ(ierr);
    
  //define problems
  ierr = mField.add_problem("ELASTIC_MECHANICS"); CHKERRQ(ierr);

  //set finite elements for problem
  ierr = mField.modify_problem_add_finite_element("ELASTIC_MECHANICS","ELASTIC"); CHKERRQ(ierr);
  ierr = mField.modify_problem_add_finite_element("ELASTIC_MECHANICS","TRAN_ISOTROPIC_ELASTIC"); CHKERRQ(ierr);

  //set refinment level for problem
  ierr = mField.modify_problem_ref_level_add_bit("ELASTIC_MECHANICS",bit_level0); CHKERRQ(ierr);

  /***/
  //Declare problem

  //add entitities (by tets) to the field
  ierr = mField.add_ents_to_field_by_TETs(0,"DISPLACEMENT"); CHKERRQ(ierr);

  //add finite elements entities
//  ierr = mField.add_ents_to_finite_element_EntType_by_bit_ref(bit_level0,"ELASTIC",MBTET); CHKERRQ(ierr);
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
  Vec F;
  ierr = mField.VecCreateGhost("ELASTIC_MECHANICS",Row,&F); CHKERRQ(ierr);
  Mat Aij;
  ierr = mField.MatCreateMPIAIJWithArrays("ELASTIC_MECHANICS",&Aij); CHKERRQ(ierr);

  //Get SideSet 1 and SideSet 2 defined in CUBIT
  Range SideSet1,SideSet2,SideSet3,SideSet4,NodeSet1;
  ierr = mField.get_Cubit_msId_entities_by_dimension(1,SideSet,2,SideSet1,true); CHKERRQ(ierr);
  ierr = mField.get_Cubit_msId_entities_by_dimension(2,SideSet,2,SideSet2,true); CHKERRQ(ierr);
  ierr = mField.get_Cubit_msId_entities_by_dimension(3,SideSet,1,SideSet3,true); CHKERRQ(ierr);
  ierr = mField.get_Cubit_msId_entities_by_dimension(4,SideSet,1,SideSet4,true); CHKERRQ(ierr);
  ierr = mField.get_Cubit_msId_entities_by_dimension(1,NodeSet,0,NodeSet1,true); CHKERRQ(ierr);

  PetscPrintf(PETSC_COMM_WORLD,"Nb. faces in SideSet 1 : %u\n",SideSet1.size());
  PetscPrintf(PETSC_COMM_WORLD,"Nb. faces in SideSet 2 : %u\n",SideSet2.size());
  PetscPrintf(PETSC_COMM_WORLD,"Nb. edges in SideSet 3 : %u\n",SideSet3.size());
  PetscPrintf(PETSC_COMM_WORLD,"Nb. edges in SideSet 4 : %u\n",SideSet4.size());
  PetscPrintf(PETSC_COMM_WORLD,"Nb. nodes in SideSet 1 : %u\n",NodeSet1.size());

/*  struct MyElasticFEMethod: public ElasticFEMethod {
    MyElasticFEMethod(Interface& _moab,BaseDirihletBC *_dirihlet_ptr,
      Mat &_Aij,Vec& _F,
      double _lambda,double _mu,Range &_SideSet1,Range &_SideSet2): 
      ElasticFEMethod(_moab,_dirihlet_ptr,_Aij,_F,_lambda,_mu,
      _SideSet1,_SideSet2) {};

    /// Set Neumann Boundary Conditions on SideSet2
    PetscErrorCode NeumannBC() {
      PetscFunctionBegin;
      ublas::vector<FieldData,ublas::bounded_array<double,3> > traction(3);
      //Set Direction of Traction On SideSet2
      traction[0] = 0; //X
      traction[1] = 1; //Y 
      traction[2] = 0; //Z
      //ElasticFEMethod::NeumannBC(...) function calulating external forces (see file ElasticFEMethod.hpp)
      ierr = ElasticFEMethod::NeumannBC(F,traction,SideSet2); CHKERRQ(ierr);
      PetscFunctionReturn(0);
    }

    PetscErrorCode operator()() {
      PetscFunctionBegin;
      ierr = OpStudentStart_TET(g_NTET); CHKERRQ(ierr);
      ierr = GetMatrices(); CHKERRQ(ierr);

      //Dirihlet Boundary Condition
      ierr = SetDirihletBC_to_ElementIndicies(); CHKERRQ(ierr);
      if(Diagonal!=PETSC_NULL) {
	if(DirihletBC.size()>0) {
	  DirihletBCDiagVal.resize(DirihletBC.size());
	  fill(DirihletBCDiagVal.begin(),DirihletBCDiagVal.end(),1);
	  ierr = VecSetValues(Diagonal,DirihletBC.size(),&(DirihletBC[0]),&DirihletBCDiagVal[0],INSERT_VALUES); CHKERRQ(ierr);
	}
      }

      //Assembly Aij and F
      ierr = RhsAndLhs(); CHKERRQ(ierr);

      //Neumann Boundary Conditions
      ierr = NeumannBC(); CHKERRQ(ierr);

      ierr = OpStudentEnd(); CHKERRQ(ierr);
      PetscFunctionReturn(0);
    }

  };*/

  //Assemble F and Aij
  const double YoungModulusP = 1.0;
  const double PoissonRatioP = 0.3;
  const double YoungModulusZ = 1.0;
  const double PoissonRatioPZ = 0.3;
  const double ShearModulusZP = 0.5;
    
  //Note: if angles X,Y & Z are all zero, the fibre direction is in the y-axis
  const double AngleX = (45.0 * M_PI)/180; //theta
  const double AngleY = (0.0 * M_PI)/180; //alpha
  const double AngleZ = (0.0 * M_PI)/180; //psi

    struct MyDirihletBC: public BaseDirihletBC {
        Range& SideSet1,SideSet3,SideSet4,NodeSet1;
        Range SideSet1_,SideSet3_,SideSet4_;
        
        // Constructor
        MyDirihletBC(Interface &moab,Range& _SideSet1,Range& _SideSet3,Range& _SideSet4, Range& _NodeSet1): 
        BaseDirihletBC(),SideSet1(_SideSet1), SideSet3(_SideSet3), SideSet4(_SideSet4),NodeSet1(_NodeSet1){
            
            //Add to SideSet1_ nodes,edges, and faces, where dirihilet boundary conditions are applied.
            //Note that SideSet1 consist only faces in this particular example.
            ErrorCode rval;
            Range SideSet1Edges,SideSet1Nodes;
            Range SideSet3Edges,SideSet3Nodes;
            Range SideSet4Edges,SideSet4Nodes;
            
            rval = moab.get_adjacencies(SideSet1,1,false,SideSet1Edges,Interface::UNION); CHKERR_THROW(rval);
            rval = moab.get_connectivity(SideSet1,SideSet1Nodes,true); CHKERR_THROW(rval);
            SideSet1_.insert(SideSet1.begin(),SideSet1.end());
            SideSet1_.insert(SideSet1Edges.begin(),SideSet1Edges.end());
            SideSet1_.insert(SideSet1Nodes.begin(),SideSet1Nodes.end());
            
            rval = moab.get_adjacencies(SideSet3,1,false,SideSet3Edges,Interface::UNION); CHKERR_THROW(rval);
            rval = moab.get_connectivity(SideSet3,SideSet3Nodes,true); CHKERR_THROW(rval);
            SideSet3_.insert(SideSet3.begin(),SideSet3.end());
            SideSet3_.insert(SideSet3Edges.begin(),SideSet3Edges.end());
            SideSet3_.insert(SideSet3Nodes.begin(),SideSet3Nodes.end());
            
            rval = moab.get_adjacencies(SideSet4,1,false,SideSet4Edges,Interface::UNION); CHKERR_THROW(rval);
            rval = moab.get_connectivity(SideSet4,SideSet4Nodes,true); CHKERR_THROW(rval);
            SideSet4_.insert(SideSet4.begin(),SideSet4.end());
            SideSet4_.insert(SideSet4Edges.begin(),SideSet4Edges.end());
            SideSet4_.insert(SideSet4Nodes.begin(),SideSet4Nodes.end()); 
        }
        
        //This method is called insiaid finite element loop to apply boundary conditions
        PetscErrorCode SetDirihletBC_to_ElementIndicies(
                                                        moabField::FEMethod *fe_method_ptr,string field_name,
                                                        vector<vector<DofIdx> > &RowGlob,vector<vector<DofIdx> > &ColGlob,vector<DofIdx>& DirihletBC) {
            PetscFunctionBegin;
            
            ierr = InternalClassBCSet(fe_method_ptr,RowGlob,ColGlob,DirihletBC,field_name,SideSet1_,fixed_y,false); CHKERRQ(ierr);
            ierr = InternalClassBCSet(fe_method_ptr,RowGlob,ColGlob,DirihletBC,field_name,SideSet3_,fixed_x,false); CHKERRQ(ierr);
            ierr = InternalClassBCSet(fe_method_ptr,RowGlob,ColGlob,DirihletBC,field_name,SideSet4_,fixed_z,false); CHKERRQ(ierr);
            ierr = InternalClassBCSet(fe_method_ptr,RowGlob,ColGlob,DirihletBC,field_name,NodeSet1,fixed_x|fixed_y|fixed_z,false); CHKERRQ(ierr);
            
            // You can call this functitions more times for other Rnages (SideSets)
            //
            // For example ... fixing only x-direction
            //
            // ierr = __SetDirihletBC_to_ElementIndicies__("DISPLACEMENT",RowGlob,ColGlob,DirihletBC,field_name,SideSet2_,fixed_x,false); CHKERRQ(ierr);
            //
            // or fixing x- and y- direction
            //
            // ierr = __SetDirihletBC_to_ElementIndicies__("DISPLACEMENT",RowGlob,ColGlob,DirihletBC,field_name,SideSet3_,fixed_x|fixed_y,false); CHKERRQ(ierr);
            //
            // Note that for supsequent calss of __SetDirihletBC_to_ElementIndicies__ last paramater have to be set to false,
            // indicating that you adding boundary conditions.
            
            PetscFunctionReturn(0);
            
            
        }
        
    private:
        //Only to use in this auxiliary function
        enum bc_type { fixed_x = 1,fixed_y = 1<<1, fixed_z = 1<<2 };
        PetscErrorCode InternalClassBCSet(
                                          moabField::FEMethod *fe_method_ptr,vector<vector<DofIdx> > &RowGlob,vector<vector<DofIdx> > &ColGlob,vector<DofIdx>& DirihletBC,
                                          string field_name,Range &SideSet,unsigned int bc = fixed_x|fixed_y|fixed_z,bool zero_bc = true) {
            PetscFunctionBegin;
            //Dirihlet form SideSet1
            if(zero_bc) DirihletBC.resize(0);
            Range::iterator siit1 = SideSet.begin();
            for(;siit1!=SideSet.end();siit1++) {
                FENumeredDofMoFEMEntity_multiIndex::index<MoABEnt_mi_tag>::type::iterator riit = fe_method_ptr->row_multiIndex->get<MoABEnt_mi_tag>().lower_bound(*siit1);
                FENumeredDofMoFEMEntity_multiIndex::index<MoABEnt_mi_tag>::type::iterator hi_riit = fe_method_ptr->row_multiIndex->get<MoABEnt_mi_tag>().upper_bound(*siit1);
                for(;riit!=hi_riit;riit++) {
                    if(riit->get_name()!=field_name) continue;
                    unsigned int my_bc = 0;
                    switch(riit->get_dof_rank()) {
                        case 0: my_bc = fixed_x; break;
                        case 1: my_bc = fixed_y; break;
                        case 2: my_bc = fixed_z; break;
                        default:
                            SETERRQ(PETSC_COMM_SELF,1,"not implemented");
                    }
                    if((bc&my_bc) == 0) continue;
                    // all fixed
                    // if some ranks are selected then we could apply BC in particular direction
                    DirihletBC.push_back(riit->get_petsc_gloabl_dof_idx());
                    for(unsigned int cc = 0;cc<ColGlob.size();cc++) {
                        vector<DofIdx>::iterator it = find(ColGlob[cc].begin(),ColGlob[cc].end(),riit->get_petsc_gloabl_dof_idx());
                        if( it!=ColGlob[cc].end() ) *it = -1; // of idx is set -1 column is not assembled
                    }
                    for(unsigned int rr = 0;rr<RowGlob.size();rr++) {
                        vector<DofIdx>::iterator it = find(RowGlob[rr].begin(),RowGlob[rr].end(),riit->get_petsc_gloabl_dof_idx());
                        if( it!=RowGlob[rr].end() ) *it = -1; // of idx is set -1 row is not assembled
                    }
                }
            }
            PetscFunctionReturn(0);
        }
        
    };
    
    MyDirihletBC myDirihletBC(moab,SideSet1,SideSet3,SideSet4,NodeSet1);

//    ExampleDiriheltBC myDirihletBC(moab,SideSet1);
//  MyElasticFEMethod MyFE(moab,&myDirihletBC,Aij,F,LAMBDA(YoungModulus,PoissonRatio),MU(YoungModulus,PoissonRatio),SideSet1,SideSet2);
//  MyElasticFEMethod MyTIsotFE(moab,&myDirihletBC,Aij,F,LAMBDA(YoungModulusP,PoissonRatioP),MU(YoungModulusP,PoissonRatioP),SideSet1,SideSet2);
    TranIsotropicElasticFEMethod MyTIsotFE(moab,&myDirihletBC,Aij,F,LAMBDA(YoungModulusP,PoissonRatioP),MU(YoungModulusP,PoissonRatioP),YoungModulusP,YoungModulusZ,DELTA(PoissonRatioP,PoissonRatioPZ,YoungModulusP,YoungModulusZ),PoissonRatioP,PoissonRatioPZ,ShearModulusZP,AngleX,AngleY,AngleZ,SideSet1,SideSet2);

  ierr = VecZeroEntries(F); CHKERRQ(ierr);
  ierr = VecGhostUpdateBegin(F,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = VecGhostUpdateEnd(F,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = MatZeroEntries(Aij); CHKERRQ(ierr);

//  ierr = mField.loop_finite_elements("ELASTIC_MECHANICS","ELASTIC",MyFE);  CHKERRQ(ierr);
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

  Vec D;
  ierr = VecDuplicate(F,&D); CHKERRQ(ierr);
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
    ierr = mField.problem_get_FE("ELASTIC_MECHANICS","ELASTIC",out_meshset); CHKERRQ(ierr);
    ierr = mField.problem_get_FE("ELASTIC_MECHANICS","TRAN_ISOTROPIC_ELASTIC",out_meshset); CHKERRQ(ierr);
    rval = moab.write_file("out.vtk","VTK","",&out_meshset,1); CHKERR_PETSC(rval);
    rval = moab.delete_entities(&out_meshset,1); CHKERR_PETSC(rval);
  }

  //PostProcDisplacemenysAndStarinOnRefMesh fe_post_proc_method(moab);
//  PostProcDisplacemenysAndStarinAndElasticLinearStressOnRefMesh fe_post_proc_method(moab,LAMBDA(YoungModulusP,PoissonRatioP),MU(YoungModulusP,PoissonRatioP));
    TranIsotropicPostProcDisplacemenysAndStarinAndElasticLinearStressOnRefMesh fe_post_proc_method( moab, LAMBDA(YoungModulusP,PoissonRatioP),MU(YoungModulusP,PoissonRatioP), YoungModulusP,YoungModulusZ,DELTA(PoissonRatioP,PoissonRatioPZ,YoungModulusP,YoungModulusZ), PoissonRatioP,PoissonRatioPZ,ShearModulusZP,AngleX,AngleY,AngleZ);
    
  ierr = mField.loop_finite_elements("ELASTIC_MECHANICS","ELASTIC",fe_post_proc_method);  CHKERRQ(ierr);
  ierr = mField.loop_finite_elements("ELASTIC_MECHANICS","TRAN_ISOTROPIC_ELASTIC",fe_post_proc_method);  CHKERRQ(ierr);

  PetscSynchronizedFlush(PETSC_COMM_WORLD);
  if(pcomm->rank()==0) {
    rval = fe_post_proc_method.moab_post_proc.write_file("out_post_proc.vtk","VTK",""); CHKERR_PETSC(rval);
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

