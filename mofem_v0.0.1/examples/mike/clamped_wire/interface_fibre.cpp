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

#include "ElasticFEMethod.hpp"
#include "PostProcVertexMethod.hpp"
#include "PostProcDisplacementAndStrainOnRefindedMesh.hpp"

#include "ElasticFEMethodForInterface.hpp"

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
        
        const double nu_zp=(nu_pz*E_z)/E_p;
        const double delta=((1+nu_p)*(1-nu_p-2*nu_pz*(nu_pz*E_z/E_p)))/(E_p*E_p*E_z);
        
        StiffnessMatrix.resize(6);
        StiffnessMatrix.clear();
        StiffnessMatrix(0,0)=StiffnessMatrix(1,1)=(1-nu_pz*nu_zp)/(E_p*E_z*delta);
        StiffnessMatrix(2,2)=(1-nu_p*nu_p)/(E_p*E_p*delta);
        
        StiffnessMatrix(0,1)=StiffnessMatrix(1,0)=(nu_p+nu_zp*nu_pz)/(E_p*E_z*delta);
        StiffnessMatrix(0,2)=StiffnessMatrix(2,0)=StiffnessMatrix(1,2)=StiffnessMatrix(2,1)=(nu_zp+nu_p*nu_zp)/(E_p*E_z*delta);
        
        StiffnessMatrix(4,4)=E_p/(2*(1+nu_p));
        StiffnessMatrix(5,5)=StiffnessMatrix(3,3)=G_zp;
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
 * \brief Function to Calculate the Rotation Matrix at every gauss point 
 * to rotate the Transverse Isotroptic Stiffness Matrix in the direction of the fibre
 * z-direction of the stiffness matrix would be orientated along the fibre
 * This transformation depends on the geometry and hence, in the case
 * the pitch size of a 3 strand wire was specified, while the wire was oriented along the z-direction
 
 *\param coordinateX is the x-coordinate of the gauss point at which the rotational matrix is constructed
 *\param coordinateZ is the y-coordinate of the gauss point at which the rotational matrix is constructed
 */
struct RotationMatrixForTransverseIsotropy {
    
    double coordinatesX, coordinatesY;
    
    ublas::matrix<double> TrpMatrix,Rotational_Matrix; 
    
    RotationMatrixForTransverseIsotropy(double coordinatesX, double coordinatesY){
                
        //setting vectors of rotation and angle (axis-angle rotation)
        double AxVector1[3]={0,1,0};
        double norm_AxVector1 = 1.0;
        double AxAngle1 = 0.5*M_PI;
        
        double AxVector2[3]={coordinatesX,coordinatesY,0};
        double norm_AxVector2 = sqrt(pow(AxVector2[0],2) + pow(AxVector2[1],2) + pow(AxVector2[2],2)); 
        const double pitch = 150;
        double radius = sqrt(pow(coordinatesX,2)+pow(coordinatesY,2));
        double AxAngle2 = -atan(pitch/(2*M_PI*radius)); //Angle of helix
        
        double AxVector3[3]={0,0,1};
        double norm_AxVector3 = sqrt(pow(AxVector3[0],2) + pow(AxVector3[1],2) + pow(AxVector3[2],2));
        double AxAngle3 = -acos(coordinatesY/sqrt(pow(coordinatesX,2)+pow(coordinatesY,2)));
        if (coordinatesX<0) {
            AxAngle3=-AxAngle3;
        }else{}
        
        //Rotational Matrices
        ublas::matrix<double> RotationC; //Rotate about z-axis by 90 degrees
        RotationC = ublas::zero_matrix<FieldData>(3,3);
        ublas::matrix<double> RotationD; //Rotate about z-axis by the angle between [0,1,0] and [coordX,coordY,0]
        RotationD = ublas::zero_matrix<FieldData>(3,3);
        ublas::matrix<double> RotationE; //Rotate about [coordX,coordY,0] by angle of helix
        RotationE = ublas::zero_matrix<FieldData>(3,3);
        
        RotationC(0,0) = 1-((1-cos(AxAngle1))*(pow(AxVector1[1],2)+pow(AxVector1[2],2))/pow(norm_AxVector1,2)); 
        RotationC(1,1) = 1-((1-cos(AxAngle1))*(pow(AxVector1[0],2)+pow(AxVector1[2],2))/pow(norm_AxVector1,2)); 
        RotationC(2,2) = 1-((1-cos(AxAngle1))*(pow(AxVector1[0],2)+pow(AxVector1[1],2))/pow(norm_AxVector1,2)); 
        
        RotationC(0,1) = ((1-cos(AxAngle1))*AxVector1[0]*AxVector1[1]-norm_AxVector1*AxVector1[2]*sin(AxAngle1))/pow(norm_AxVector1,2);
        RotationC(1,0) = ((1-cos(AxAngle1))*AxVector1[0]*AxVector1[1]+norm_AxVector1*AxVector1[2]*sin(AxAngle1))/pow(norm_AxVector1,2);
        
        RotationC(0,2) = ((1-cos(AxAngle1))*AxVector1[0]*AxVector1[2]+norm_AxVector1*AxVector1[1]*sin(AxAngle1))/pow(norm_AxVector1,2);
        RotationC(2,0) = ((1-cos(AxAngle1))*AxVector1[0]*AxVector1[2]-norm_AxVector1*AxVector1[1]*sin(AxAngle1))/pow(norm_AxVector1,2);
        
        RotationC(1,2) = ((1-cos(AxAngle1))*AxVector1[1]*AxVector1[2]-norm_AxVector1*AxVector1[0]*sin(AxAngle1))/pow(norm_AxVector1,2);
        RotationC(2,1) = ((1-cos(AxAngle1))*AxVector1[1]*AxVector1[2]+norm_AxVector1*AxVector1[0]*sin(AxAngle1))/pow(norm_AxVector1,2); 
        
        RotationD(0,0) = 1-((1-cos(AxAngle2))*(pow(AxVector2[1],2)+pow(AxVector2[2],2))/pow(norm_AxVector2,2)); 
        RotationD(1,1) = 1-((1-cos(AxAngle2))*(pow(AxVector2[0],2)+pow(AxVector2[2],2))/pow(norm_AxVector2,2)); 
        RotationD(2,2) = 1-((1-cos(AxAngle2))*(pow(AxVector2[0],2)+pow(AxVector2[1],2))/pow(norm_AxVector2,2)); 
        
        RotationD(0,1) = ((1-cos(AxAngle2))*AxVector2[0]*AxVector2[1]-norm_AxVector2*AxVector2[2]*sin(AxAngle2))/pow(norm_AxVector2,2);
        RotationD(1,0) = ((1-cos(AxAngle2))*AxVector2[0]*AxVector2[1]+norm_AxVector2*AxVector2[2]*sin(AxAngle2))/pow(norm_AxVector2,2);
        
        RotationD(0,2) = ((1-cos(AxAngle2))*AxVector2[0]*AxVector2[2]+norm_AxVector2*AxVector2[1]*sin(AxAngle2))/pow(norm_AxVector2,2);
        RotationD(2,0) = ((1-cos(AxAngle2))*AxVector2[0]*AxVector2[2]-norm_AxVector2*AxVector2[1]*sin(AxAngle2))/pow(norm_AxVector2,2);
        
        RotationD(1,2) = ((1-cos(AxAngle2))*AxVector2[1]*AxVector2[2]-norm_AxVector2*AxVector2[0]*sin(AxAngle2))/pow(norm_AxVector2,2);
        RotationD(2,1) = ((1-cos(AxAngle2))*AxVector2[1]*AxVector2[2]+norm_AxVector2*AxVector2[0]*sin(AxAngle2))/pow(norm_AxVector2,2); 
        
        RotationE(0,0) = 1-((1-cos(AxAngle3))*(pow(AxVector3[1],2)+pow(AxVector3[2],2))/pow(norm_AxVector3,2)); 
        RotationE(1,1) = 1-((1-cos(AxAngle3))*(pow(AxVector3[0],2)+pow(AxVector3[2],2))/pow(norm_AxVector3,2)); 
        RotationE(2,2) = 1-((1-cos(AxAngle3))*(pow(AxVector3[0],2)+pow(AxVector3[1],2))/pow(norm_AxVector3,2)); 
        
        RotationE(0,1) = ((1-cos(AxAngle3))*AxVector3[0]*AxVector3[1]-norm_AxVector3*AxVector3[2]*sin(AxAngle3))/pow(norm_AxVector3,2);
        RotationE(1,0) = ((1-cos(AxAngle3))*AxVector3[0]*AxVector3[1]+norm_AxVector3*AxVector3[2]*sin(AxAngle3))/pow(norm_AxVector3,2);
        
        RotationE(0,2) = ((1-cos(AxAngle3))*AxVector3[0]*AxVector3[2]+norm_AxVector3*AxVector3[1]*sin(AxAngle3))/pow(norm_AxVector3,2);
        RotationE(2,0) = ((1-cos(AxAngle3))*AxVector3[0]*AxVector3[2]-norm_AxVector3*AxVector3[1]*sin(AxAngle3))/pow(norm_AxVector3,2);
        
        RotationE(1,2) = ((1-cos(AxAngle3))*AxVector3[1]*AxVector3[2]-norm_AxVector3*AxVector3[0]*sin(AxAngle3))/pow(norm_AxVector3,2);
        RotationE(2,1) = ((1-cos(AxAngle3))*AxVector3[1]*AxVector3[2]+norm_AxVector3*AxVector3[0]*sin(AxAngle3))/pow(norm_AxVector3,2); 
        
        //Combine Rotational Matrices to rotate the stiffness matrix
        ublas::matrix<double> Rotational_Matrix1 = prod(RotationE,RotationC);
        Rotational_Matrix = ublas::zero_matrix<FieldData>(3,3);
        Rotational_Matrix = prod(RotationD,Rotational_Matrix1);            
        
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
        
    };
};

struct TranIsotropicElasticFEMethod: public ElasticFEMethod {
    
    double nu_p, nu_pz, E_p, E_z, G_zp;
        
    TranIsotropicElasticFEMethod(
                                 Interface& _moab,BaseDirihletBC *_dirihlet_ptr,Mat &_Aij,Vec& _F,
                                 double _lambda,double _mu,double _E_p,double _E_z, double _nu_p,double _nu_pz, double _G_zp, Range &_SideSet1,Range &_SideSet2 ): 
    ElasticFEMethod(_moab,_dirihlet_ptr,_Aij,_F,_lambda,_mu,_SideSet1,_SideSet2), E_p(_E_p), E_z(_E_z), nu_p(_nu_p), nu_pz(_nu_pz), G_zp(_G_zp) {};
        
    vector< ublas::symmetric_matrix<FieldData,ublas::upper> > D_At_GaussPoint;

    PetscErrorCode NeumannBC() {
        PetscFunctionBegin;
        ublas::vector<FieldData,ublas::bounded_array<double,3> > traction2(3);
        //Set Direction of Traction On SideSet2
        traction2[0] = 0; //X
        traction2[1] = 0; //Y 
        traction2[2] = +100; //Z
        //ElasticFEMethod::NeumannBC(...) function calulating external forces (see file ElasticFEMethod.hpp)
        ierr = ElasticFEMethod::NeumannBC(F,traction2,SideSet2); CHKERRQ(ierr);
        
        PetscFunctionReturn(0);
    }
    
    PetscErrorCode Stiffness() {
        PetscFunctionBegin;
        K.resize(row_mat,col_mat);
        int g_dim = g_NTET.size()/4;
        ublas::matrix<FieldData> BTD;
        for(int rr = 0;rr<row_mat;rr++) {
            for(int gg = 0;gg<g_dim;gg++) {
                ublas::matrix<FieldData> &row_Mat = (rowBMatrices[rr])[gg];
                double w = V*G_W_TET[gg];
                BTD.resize(row_Mat.size2(),6);
                ublas::noalias(BTD) = prod( trans(row_Mat), w*D_At_GaussPoint[gg] );
                for(int cc = rr;cc<col_mat;cc++) {
                    ublas::matrix<FieldData> &col_Mat = (colBMatrices[cc])[gg];
                    if(gg == 0) {
                        K(rr,cc).resize(BTD.size1(),col_Mat.size2());
                        ublas::noalias(K(rr,cc)) = prod(BTD , col_Mat ); // int BT*D*B
                    } else {
                        ublas::noalias(K(rr,cc)) += prod(BTD , col_Mat ); // int BT*D*B
                    }
                }
            }
        }
        PetscFunctionReturn(0);
    }
    
    PetscErrorCode operator()() {
        PetscFunctionBegin;
        ierr = OpStudentStart_TET(g_NTET); CHKERRQ(ierr);
        ierr = GetMatrices(); CHKERRQ(ierr);
        
        D_At_GaussPoint.resize(coords_at_Gauss_nodes.size());

        for(int gg=0;gg<coords_at_Gauss_nodes.size();gg++){
            
            ///Get the Axis and Angles according to the position of gauss point
            double coordinates[3]={(coords_at_Gauss_nodes[gg]).data()[0],(coords_at_Gauss_nodes[gg]).data()[1],(coords_at_Gauss_nodes[gg]).data()[2]};
            
            ///Get Rotation matrix according to coordinate of Gauss Point
            ublas::matrix<double> TrpMatrix; 
            TrpMatrix = ublas::zero_matrix<FieldData>(6,6);
            RotationMatrixForTransverseIsotropy RotMat(coordinates[0],coordinates[1]);
            TrpMatrix=RotMat.TrpMatrix;

            ///Get Stiffness Matrix
            ublas::symmetric_matrix<FieldData,ublas::upper> StiffnessMatrix;
            StiffnessMatrix.resize(6);
            StiffnessMatrix.clear();
            TransverseIsotropicStiffnessMatrix TranIsoMat(nu_p, nu_pz, E_p, E_z, G_zp);
            StiffnessMatrix=TranIsoMat.StiffnessMatrix;
//            IsotropicStiffnessMatrix IsoMat(lambda, mu);
//            StiffnessMatrix=IsoMat.StiffnessMatrix;
            
            ///Rotating the Stiffness matrix according to the fibre direction
            D_At_GaussPoint[gg].resize(6);
            D_At_GaussPoint[gg].clear();
            ublas::matrix< FieldData > dummy2 = prod( StiffnessMatrix , TrpMatrix );
            D_At_GaussPoint[gg] = prod( trans(TrpMatrix) , dummy2 );
            
        }
        
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

struct TranIsotropicPostProcDisplacemenysAndStarinAndElasticLinearStressOnRefMesh: public PostProcDisplacemenysAndStarinAndElasticLinearStressOnRefMesh {
    
    double nu_p, nu_pz, E_p, E_z, G_zp;
    Tag th_fibre_orientation;
    
    TranIsotropicPostProcDisplacemenysAndStarinAndElasticLinearStressOnRefMesh( Interface& _moab,double _lambda,double _mu, double _E_p,double _E_z, double _nu_p, double _nu_pz, double _G_zp):
    PostProcDisplacemenysAndStarinAndElasticLinearStressOnRefMesh(_moab,_lambda,_mu),E_p(_E_p),E_z(_E_z),nu_p(_nu_p),nu_pz(_nu_pz),G_zp(_G_zp) {
        
        double def_VAL2[3] = {0,0,0};
        rval = moab_post_proc.tag_get_handle("FIBRE_DIRECTION",3,MB_TYPE_DOUBLE,th_fibre_orientation,MB_TAG_CREAT|MB_TAG_SPARSE,&def_VAL2); CHKERR_THROW(rval);
    }
    
    PetscErrorCode operator()() {
        PetscFunctionBegin;
        
        ierr = do_operator(); CHKERRQ(ierr);
        ierr = OpStudentStart_TET(g_NTET); CHKERRQ(ierr);
        
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
            
            ///Get the Axis and Angles according to the position of gauss point
            double coordinates[3]={(coords_at_Gauss_nodes[gg]).data()[0],(coords_at_Gauss_nodes[gg]).data()[1],(coords_at_Gauss_nodes[gg]).data()[2]};
            
            ///Get Rotation matrix according to coordinate of Gauss Point
            ublas::matrix<double> TrpMatrix,Rotational_Matrix;
            TrpMatrix = ublas::zero_matrix<FieldData>(6,6);
            Rotational_Matrix = ublas::zero_matrix<FieldData>(3,3);
            RotationMatrixForTransverseIsotropy RotMat(coordinates[0],coordinates[1]);
            
            TrpMatrix=RotMat.TrpMatrix;
            Rotational_Matrix=RotMat.Rotational_Matrix;
            
            ///Rotate AxisYVector[0,1,0] to the direction of the fibre and save in TAG
            ublas::vector<FieldData> AxisYVector(3);
            AxisYVector[0]=0; AxisYVector[1]=0;AxisYVector[2]=1;
            ublas::vector<FieldData> Fibre = prod(Rotational_Matrix,AxisYVector);
            
            rval = moab_post_proc.tag_set_data(th_fibre_orientation,&mit->second,1,&Fibre[0]); CHKERR_PETSC(rval);
            
            ///Get Stiffness Matrix
            ublas::symmetric_matrix<FieldData,ublas::upper> StiffnessMatrix;
            StiffnessMatrix.resize(6);
            StiffnessMatrix.clear();
            TransverseIsotropicStiffnessMatrix TranIsoMat(nu_p, nu_pz, E_p, E_z, G_zp);
            StiffnessMatrix=TranIsoMat.StiffnessMatrix;            
//            IsotropicStiffnessMatrix IsoMat(lambda, mu);
//            StiffnessMatrix=IsoMat.StiffnessMatrix;
            
            ///Rotating the Stiffness matrix according to the fibre direction
            D.resize(6,6);
            D.clear();
            ublas::matrix< FieldData > dummy2 = prod( StiffnessMatrix , TrpMatrix );
            D=prod( trans(TrpMatrix) , dummy2 );
            
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
     PetscInt order;
  ierr = PetscOptionsGetInt(PETSC_NULL,"-my_order",&order,&flg); CHKERRQ(ierr);
  if(flg != PETSC_TRUE) {
    order = 3;
  } 
  
   char outName[PETSC_MAX_PATH_LEN]="out.vtk";
     ierr = PetscOptionsGetString(PETSC_NULL,"-my_out",outName,sizeof(outName),&flg); CHKERRQ(ierr);
     
       char outName2[PETSC_MAX_PATH_LEN]="out_post_proc.vtk";
     ierr = PetscOptionsGetString(PETSC_NULL,"-my_post_out",outName2,sizeof(outName2),&flg); CHKERRQ(ierr);

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

    EntityHandle meshset_SideSet1; //Dirihlet BC is there
    ierr = mField.get_msId_meshset(1,SideSet,meshset_SideSet1); CHKERRQ(ierr);
    ierr = mField.refine_get_childern(meshset_SideSet1,bit_level0,meshset_SideSet1,MBTRI,true); CHKERRQ(ierr);
    EntityHandle meshset_SideSet2;
    ierr = mField.get_msId_meshset(2,SideSet,meshset_SideSet2); CHKERRQ(ierr);
    ierr = mField.refine_get_childern(meshset_SideSet2,bit_level0,meshset_SideSet2,MBTRI,true); CHKERRQ(ierr);
//    EntityHandle meshset_SideSet5;
//    ierr = mField.get_msId_meshset(5,SideSet,meshset_SideSet5); CHKERRQ(ierr);
//    ierr = mField.refine_get_childern(meshset_SideSet5,bit_level0,meshset_SideSet5,MBTRI,true); CHKERRQ(ierr);
//    EntityHandle meshset_SideSet6;
//    ierr = mField.get_msId_meshset(6,SideSet,meshset_SideSet6); CHKERRQ(ierr);
//    ierr = mField.refine_get_childern(meshset_SideSet6,bit_level0,meshset_SideSet6,MBTRI,true); CHKERRQ(ierr);
//    EntityHandle meshset_NodeSet1; 
//    ierr = mField.get_msId_meshset(1,NodeSet,meshset_NodeSet1); CHKERRQ(ierr);
//    ierr = mField.refine_get_childern(meshset_NodeSet1,bit_level0,meshset_NodeSet1,MBVERTEX,true); CHKERRQ(ierr);
    EntityHandle meshset_BlockSet1; 
    ierr = mField.get_msId_meshset(1,BlockSet,meshset_BlockSet1); CHKERRQ(ierr);
    ierr = mField.refine_get_childern(meshset_BlockSet1,bit_level0,meshset_BlockSet1,MBTET,true); CHKERRQ(ierr);
    EntityHandle meshset_BlockSet2; 
    ierr = mField.get_msId_meshset(2,BlockSet,meshset_BlockSet2); CHKERRQ(ierr);
    ierr = mField.refine_get_childern(meshset_BlockSet2,bit_level0,meshset_BlockSet2,MBTET,true); CHKERRQ(ierr);

    
    EntityHandle meshset_level0;
    rval = moab.create_meshset(MESHSET_SET,meshset_level0); CHKERR_PETSC(rval);
    ierr = mField.seed_ref_level_3D(0,bit_level0); CHKERRQ(ierr);
    ierr = mField.refine_get_ents(bit_level0,meshset_level0); CHKERRQ(ierr);
    
    BitRefLevel bit_level1;
    bit_level1.set(1);    

    ierr = mField.add_verices_in_the_middel_of_edges(meshset_level0,bit_level1); CHKERRQ(ierr);
    ierr = mField.refine_TET(meshset_level0,bit_level1); CHKERRQ(ierr);
    ierr = mField.refine_get_childern(meshset_SideSet1, bit_level1,meshset_SideSet1,MBTRI,true); CHKERRQ(ierr);
    ierr = mField.refine_get_childern(meshset_SideSet2, bit_level1,meshset_SideSet2,MBTRI,true); CHKERRQ(ierr);
//    ierr = mField.refine_get_childern(meshset_SideSet5, bit_level1,meshset_SideSet5,MBEDGE,true); CHKERRQ(ierr);
//    ierr = mField.refine_get_childern(meshset_SideSet6, bit_level1,meshset_SideSet6,MBEDGE,true); CHKERRQ(ierr);
//    ierr = mField.refine_get_childern(meshset_NodeSet1, bit_level1,meshset_NodeSet1,MBVERTEX,true); CHKERRQ(ierr);
    ierr = mField.refine_get_childern(meshset_BlockSet1,bit_level1,meshset_BlockSet1,MBTET,true); CHKERRQ(ierr);
    ierr = mField.refine_get_childern(meshset_BlockSet2,bit_level1,meshset_BlockSet2,MBTET,true); CHKERRQ(ierr);
    
    EntityHandle meshset_level1;
    rval = moab.create_meshset(MESHSET_SET,meshset_level1); CHKERR_PETSC(rval);
    ierr = mField.seed_ref_level_3D(0,bit_level1); CHKERRQ(ierr);
    ierr = mField.refine_get_ents(bit_level1,meshset_level1); CHKERRQ(ierr);
    
    BitRefLevel bit_level2;
    bit_level2.set(2);    
    
    ierr = mField.add_verices_in_the_middel_of_edges(meshset_level1,bit_level2); CHKERRQ(ierr);
    ierr = mField.refine_TET(meshset_level1,bit_level2); CHKERRQ(ierr);
    ierr = mField.refine_get_childern(meshset_SideSet1, bit_level2,meshset_SideSet1,MBTRI,true); CHKERRQ(ierr);
    ierr = mField.refine_get_childern(meshset_SideSet2, bit_level2,meshset_SideSet2,MBTRI,true); CHKERRQ(ierr);
//    ierr = mField.refine_get_childern(meshset_SideSet5, bit_level2,meshset_SideSet5,MBEDGE,true); CHKERRQ(ierr);
//    ierr = mField.refine_get_childern(meshset_SideSet6, bit_level2,meshset_SideSet6,MBEDGE,true); CHKERRQ(ierr);
//    ierr = mField.refine_get_childern(meshset_NodeSet1, bit_level2,meshset_NodeSet1,MBVERTEX,true); CHKERRQ(ierr);
    ierr = mField.refine_get_childern(meshset_BlockSet1,bit_level2,meshset_BlockSet1,MBTET,true); CHKERRQ(ierr);
    ierr = mField.refine_get_childern(meshset_BlockSet2,bit_level2,meshset_BlockSet2,MBTET,true); CHKERRQ(ierr);

    BitRefLevel problem_bit_level = bit_level2;
    
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
    //FE Transverse Isotropic
    ierr = mField.modify_finite_element_add_field_row("TRAN_ISOTROPIC_ELASTIC","DISPLACEMENT"); CHKERRQ(ierr);
    ierr = mField.modify_finite_element_add_field_col("TRAN_ISOTROPIC_ELASTIC","DISPLACEMENT"); CHKERRQ(ierr);
    ierr = mField.modify_finite_element_add_field_data("TRAN_ISOTROPIC_ELASTIC","DISPLACEMENT"); CHKERRQ(ierr);
    
    //define problems
    ierr = mField.add_problem("ELASTIC_MECHANICS"); CHKERRQ(ierr);
    
    //set finite elements for problem
    ierr = mField.modify_problem_add_finite_element("ELASTIC_MECHANICS","ELASTIC"); CHKERRQ(ierr);
    ierr = mField.modify_problem_add_finite_element("ELASTIC_MECHANICS","TRAN_ISOTROPIC_ELASTIC"); CHKERRQ(ierr);
    
    //set refinment level for problem
    ierr = mField.modify_problem_ref_level_add_bit("ELASTIC_MECHANICS",problem_bit_level); CHKERRQ(ierr);
    
    /***/
    //Declare problem
    
    //add entitities (by tets) to the field
    ierr = mField.add_ents_to_field_by_TETs(0,"DISPLACEMENT"); CHKERRQ(ierr);
    
    //add finite elements entities
    ierr = mField.add_ents_to_finite_element_by_TETs(meshset_BlockSet1,"TRAN_ISOTROPIC_ELASTIC",true); CHKERRQ(ierr);
    ierr = mField.add_ents_to_finite_element_by_TETs(meshset_BlockSet2,"ELASTIC",true); CHKERRQ(ierr);
    
    //set app. order
    //see Hierarchic Finite Element Bases on Unstructured Tetrahedral Meshes (Mark Ainsworth & Joe Coyle)
    ierr = mField.set_field_order(0,MBTET,"DISPLACEMENT",order); CHKERRQ(ierr);
    ierr = mField.set_field_order(0,MBTRI,"DISPLACEMENT",order); CHKERRQ(ierr);
    ierr = mField.set_field_order(0,MBEDGE,"DISPLACEMENT",order); CHKERRQ(ierr);
    ierr = mField.set_field_order(0,MBVERTEX,"DISPLACEMENT",1); CHKERRQ(ierr);
    
    /****/
    //build database
    
    //build field
    ierr = mField.build_fields(); CHKERRQ(ierr);
    
    //build finite elemnts
    ierr = mField.build_finite_elements(); CHKERRQ(ierr);
    
    //build adjacencies
    ierr = mField.build_adjacencies(problem_bit_level); CHKERRQ(ierr);
    
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
    Range SideSet1,SideSet2,SideSet5,SideSet6,NodeSet1;
    ierr = mField.get_Cubit_msId_entities_by_dimension(1,SideSet,2,SideSet1,true); CHKERRQ(ierr);
    ierr = mField.get_Cubit_msId_entities_by_dimension(2,SideSet,2,SideSet2,true); CHKERRQ(ierr);
//    ierr = mField.get_Cubit_msId_entities_by_dimension(5,SideSet,1,SideSet5,true); CHKERRQ(ierr);
//    ierr = mField.get_Cubit_msId_entities_by_dimension(6,SideSet,1,SideSet6,true); CHKERRQ(ierr);
//    ierr = mField.get_Cubit_msId_entities_by_dimension(1,NodeSet,0,NodeSet1,true); CHKERRQ(ierr);
    
    PetscPrintf(PETSC_COMM_WORLD,"Nb. faces in SideSet 1 : %u\n",SideSet1.size());
    PetscPrintf(PETSC_COMM_WORLD,"Nb. faces in SideSet 2 : %u\n",SideSet2.size());
//    PetscPrintf(PETSC_COMM_WORLD,"Nb. faces in SideSet 5 : %u\n",SideSet5.size());
//    PetscPrintf(PETSC_COMM_WORLD,"Nb. faces in SideSet 6 : %u\n",SideSet6.size());
//    PetscPrintf(PETSC_COMM_WORLD,"Nb. faces in NodeSet 1 : %u\n",NodeSet1.size());

    
    //Assemble F and Aij
    const double YoungModulusP = 135000;
    const double PoissonRatioP = 0.77;
    const double YoungModulusZ = 135000;
    const double PoissonRatioPZ = 0.2;
    const double ShearModulusZP = 5000;
    const double YoungModulus = 200000;
    const double PoissonRatio = 0.3;
    
    struct MyElasticFEMethod: public ElasticFEMethod {
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
            traction[1] = 0; //Y 
            traction[2] = +100; //Z
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
        
    };

    struct MyDirihletBC: public BaseDirihletBC {
        Range& SideSet1;
        Range SideSet1_;
        
        // Constructor
        MyDirihletBC(Interface &moab,Range& _SideSet1): 
        BaseDirihletBC(),SideSet1(_SideSet1){
            
            //Add to SideSet1_ nodes,edges, and faces, where dirihilet boundary conditions are applied.
            //Note that SideSet1 consist only faces in this particular example.
            ErrorCode rval;
            Range SideSet1Edges,SideSet1Nodes;
            
            rval = moab.get_adjacencies(SideSet1,1,false,SideSet1Edges,Interface::UNION); CHKERR_THROW(rval);
            rval = moab.get_connectivity(SideSet1,SideSet1Nodes,true); CHKERR_THROW(rval);
            SideSet1_.insert(SideSet1.begin(),SideSet1.end());
            SideSet1_.insert(SideSet1Edges.begin(),SideSet1Edges.end());
            SideSet1_.insert(SideSet1Nodes.begin(),SideSet1Nodes.end());
            
        }
        
        //This method is called insiaid finite element loop to apply boundary conditions
        PetscErrorCode SetDirihletBC_to_ElementIndicies(
                                                        moabField::FEMethod *fe_method_ptr,string field_name,
                                                        vector<vector<DofIdx> > &RowGlob,vector<vector<DofIdx> > &ColGlob,vector<DofIdx>& DirihletBC) {
            PetscFunctionBegin;
            
            ierr = InternalClassBCSet(fe_method_ptr,RowGlob,ColGlob,DirihletBC,field_name,SideSet1_,fixed_x|fixed_y|fixed_z,false); CHKERRQ(ierr);
            
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
    MyDirihletBC myDirihletBC(moab,SideSet1);
    
    /*struct MyDirihletBC: public BaseDirihletBC {
        Range& SideSet1,SideSet5,SideSet6,NodeSet1;
        Range SideSet1_,SideSet5_,SideSet6_;
        
        // Constructor
        MyDirihletBC(Interface &moab,Range& _SideSet1,Range& _SideSet5,Range& _SideSet6, Range& _NodeSet1): 
        BaseDirihletBC(),SideSet1(_SideSet1), SideSet5(_SideSet5), SideSet6(_SideSet6),NodeSet1(_NodeSet1){
            
            //Add to SideSet1_ nodes,edges, and faces, where dirihilet boundary conditions are applied.
            //Note that SideSet1 consist only faces in this particular example.
            ErrorCode rval;
            Range SideSet1Edges,SideSet1Nodes;
            Range SideSet5Edges,SideSet5Nodes;
            Range SideSet6Edges,SideSet6Nodes;
            
            rval = moab.get_adjacencies(SideSet1,1,false,SideSet1Edges,Interface::UNION); CHKERR_THROW(rval);
            rval = moab.get_connectivity(SideSet1,SideSet1Nodes,true); CHKERR_THROW(rval);
            SideSet1_.insert(SideSet1.begin(),SideSet1.end());
            SideSet1_.insert(SideSet1Edges.begin(),SideSet1Edges.end());
            SideSet1_.insert(SideSet1Nodes.begin(),SideSet1Nodes.end());
            
            rval = moab.get_adjacencies(SideSet5,1,false,SideSet5Edges,Interface::UNION); CHKERR_THROW(rval);
            rval = moab.get_connectivity(SideSet5,SideSet5Nodes,true); CHKERR_THROW(rval);
            SideSet5_.insert(SideSet5.begin(),SideSet5.end());
            SideSet5_.insert(SideSet5Edges.begin(),SideSet5Edges.end());
            SideSet5_.insert(SideSet5Nodes.begin(),SideSet5Nodes.end());
            
            rval = moab.get_adjacencies(SideSet6,1,false,SideSet6Edges,Interface::UNION); CHKERR_THROW(rval);
            rval = moab.get_connectivity(SideSet6,SideSet6Nodes,true); CHKERR_THROW(rval);
            SideSet6_.insert(SideSet6.begin(),SideSet6.end());
            SideSet6_.insert(SideSet6Edges.begin(),SideSet6Edges.end());
            SideSet6_.insert(SideSet6Nodes.begin(),SideSet6Nodes.end()); 
        }
        
        //This method is called insiaid finite element loop to apply boundary conditions
        PetscErrorCode SetDirihletBC_to_ElementIndicies(
                                                        moabField::FEMethod *fe_method_ptr,string field_name,
                                                        vector<vector<DofIdx> > &RowGlob,vector<vector<DofIdx> > &ColGlob,vector<DofIdx>& DirihletBC) {
            PetscFunctionBegin;
            
            ierr = InternalClassBCSet(fe_method_ptr,RowGlob,ColGlob,DirihletBC,field_name,SideSet1_,fixed_x|fixed_y|fixed_z,false); CHKERRQ(ierr);
//            ierr = InternalClassBCSet(fe_method_ptr,RowGlob,ColGlob,DirihletBC,field_name,SideSet5_,fixed_y,false); CHKERRQ(ierr);
//            ierr = InternalClassBCSet(fe_method_ptr,RowGlob,ColGlob,DirihletBC,field_name,SideSet6_,fixed_z,false); CHKERRQ(ierr);
//            ierr = InternalClassBCSet(fe_method_ptr,RowGlob,ColGlob,DirihletBC,field_name,NodeSet1,fixed_x|fixed_y|fixed_z,false); CHKERRQ(ierr);
            
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
    MyDirihletBC myDirihletBC(moab,SideSet1,SideSet5,SideSet6,NodeSet1);*/
    
    TranIsotropicElasticFEMethod MyTIsotFE(moab,&myDirihletBC,Aij,F,LAMBDA(YoungModulusP,PoissonRatioP),MU(YoungModulusP,PoissonRatioP),YoungModulusP,YoungModulusZ,PoissonRatioP,PoissonRatioPZ,ShearModulusZP,SideSet1,SideSet2);
    MyElasticFEMethod MyFE(moab,&myDirihletBC,Aij,F,LAMBDA(YoungModulus,PoissonRatio),MU(YoungModulus,PoissonRatio),SideSet1,SideSet2);
    
    ierr = VecZeroEntries(F); CHKERRQ(ierr);
    ierr = VecGhostUpdateBegin(F,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
    ierr = VecGhostUpdateEnd(F,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
    ierr = MatZeroEntries(Aij); CHKERRQ(ierr);
    
    ierr = mField.loop_finite_elements("ELASTIC_MECHANICS","ELASTIC",MyFE);  CHKERRQ(ierr);
    ierr = mField.loop_finite_elements("ELASTIC_MECHANICS","TRAN_ISOTROPIC_ELASTIC",MyTIsotFE);  CHKERRQ(ierr);
    PetscSynchronizedFlush(PETSC_COMM_WORLD);
    
    ierr = VecGhostUpdateBegin(F,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
    ierr = VecGhostUpdateEnd(F,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
    ierr = VecAssemblyBegin(F); CHKERRQ(ierr);
    ierr = VecAssemblyEnd(F); CHKERRQ(ierr);
    ierr = MatAssemblyBegin(Aij,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
    ierr = MatAssemblyEnd(Aij,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
    
    
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
        rval = moab.write_file(outName,"VTK","",&out_meshset,1); CHKERR_PETSC(rval);
        rval = moab.delete_entities(&out_meshset,1); CHKERR_PETSC(rval);
    }
    
    //  PostProcDisplacemenysAndStarinOnRefMesh fe_post_proc_method(moab);
    TranIsotropicPostProcDisplacemenysAndStarinAndElasticLinearStressOnRefMesh fe_post_proc_method_tran_iso( moab, LAMBDA(YoungModulusP,PoissonRatioP),MU(YoungModulusP,PoissonRatioP), YoungModulusP,YoungModulusZ,PoissonRatioP,PoissonRatioPZ,ShearModulusZP);
    PostProcDisplacemenysAndStarinAndElasticLinearStressOnRefMesh fe_post_proc_method_elastic(moab,LAMBDA(YoungModulus,PoissonRatio),MU(YoungModulus,PoissonRatio));
    ierr = mField.loop_finite_elements("ELASTIC_MECHANICS","ELASTIC",fe_post_proc_method_elastic);  CHKERRQ(ierr);
    ierr = mField.loop_finite_elements("ELASTIC_MECHANICS","TRAN_ISOTROPIC_ELASTIC",fe_post_proc_method_tran_iso);  CHKERRQ(ierr);
    
    PetscSynchronizedFlush(PETSC_COMM_WORLD);
    if(pcomm->rank()==0) {
        rval = fe_post_proc_method_elastic.fe_post_proc_method_tran_iso.moab_post_proc.write_file(outName2,"VTK",""); CHKERR_PETSC(rval);
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
