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

struct TranIsotropicElasticFEMethod: public InterfaceElasticFEMethod {
    
    double nu_p, nu_pz, E_p, E_z, G_zp;
        
    TranIsotropicElasticFEMethod(
                                 Interface& _moab,BaseDirihletBC *_dirihlet_ptr,Mat &_Aij,Vec& _F,
                                 double _lambda,double _mu,double _E_p,double _E_z, double _nu_p,double _nu_pz, double _G_zp, Range &_SideSet1,Range &_SideSet2): 
    InterfaceElasticFEMethod(_moab,_dirihlet_ptr,_Aij,_F,_lambda,_mu,_SideSet1,_SideSet2,dummy), E_p(_E_p), E_z(_E_z), nu_p(_nu_p), nu_pz(_nu_pz), G_zp(_G_zp) {};
        
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

struct TranIsotropic_Fibre_PostProcDisplacemenysAndStarinAndElasticLinearStressOnRefMesh: public PostProcDisplacemenysAndStarinAndElasticLinearStressOnRefMesh {
    
    double nu_p, nu_pz, E_p, E_z, G_zp;
    Tag th_fibre_orientation;
    Tag th_gradient;
    
    TranIsotropic_Fibre_PostProcDisplacemenysAndStarinAndElasticLinearStressOnRefMesh( Interface& _moab,double _lambda,double _mu, double _E_p,double _E_z, double _nu_p, double _nu_pz, double _G_zp):
    PostProcDisplacemenysAndStarinAndElasticLinearStressOnRefMesh(_moab,_lambda,_mu),E_p(_E_p),E_z(_E_z),nu_p(_nu_p),nu_pz(_nu_pz),G_zp(_G_zp) {
        
        double def_VAL2[3] = {0,0,0};
        rval = moab_post_proc.tag_get_handle("FIBRE_DIRECTION",3,MB_TYPE_DOUBLE,th_fibre_orientation,MB_TAG_CREAT|MB_TAG_SPARSE,&def_VAL2); CHKERR_THROW(rval);
        rval = moab_post_proc.tag_get_handle("GRADIENT",3,MB_TYPE_DOUBLE,th_gradient,MB_TAG_CREAT|MB_TAG_SPARSE,&def_VAL2); CHKERR_THROW(rval);
        
    }
    
    PetscErrorCode operator()() {
        PetscFunctionBegin;
        
        ierr = do_operator(); CHKERRQ(ierr);
        ierr = OpStudentStart_TET(g_NTET); CHKERRQ(ierr);
        
        ///Getting distance of nodes of tethra from the nearest tri fibre surface
        EntityHandle fe_handle = fe_ptr->get_ent();
        Range tetNodes;
        rval = moab.get_connectivity(&fe_handle,1,tetNodes); CHKERR_THROW(rval);
        
        ublas::vector<FieldData> surf_disp_node(4);
        int iii=0;
        Range::iterator niit1 = tetNodes.begin();
        for(;niit1!=tetNodes.end();niit1++,iii++) {
            FEDofMoFEMEntity_multiIndex::index<MoABEnt_mi_tag>::type::iterator it = data_multiIndex->get<MoABEnt_mi_tag>().lower_bound(*niit1);
            FEDofMoFEMEntity_multiIndex::index<MoABEnt_mi_tag>::type::iterator hi_it = data_multiIndex->get<MoABEnt_mi_tag>().upper_bound(*niit1);
            for(;it!=hi_it;it++) {
                
                //            for(FEDofMoFEMEntity_multiIndex::index<MoABEnt_mi_tag>::type::iterator it = data_multiIndex->get<MoABEnt_mi_tag>().lower_bound(*niit1); it != data_multiIndex->get<MoABEnt_mi_tag>().upper_bound(*niit1);++it){ 
                
                if(it->get_name()!="SURFACE_DISTANCE") continue;
                double &distVal = it->get_FieldData();
                //                cout<<it->get_name()<<"   "<<distVal<<endl;
                surf_disp_node[iii]=distVal;
            }
        }
        
        ublas::matrix<double> JacMat; //Rotate about z-axis by 90 degrees
        JacMat = ublas::zero_matrix<FieldData>(3,3);
        ublas::matrix<double> diffNMat; //Rotate about z-axis by 90 degrees
        diffNMat = ublas::zero_matrix<FieldData>(3,4);
        
        double diffN[12];
        ShapeDiffMBTET(diffN);
        for (int rr=0; rr<3; rr++) {
            for(int cc=0; cc<4; cc++) {
                diffNMat(rr,cc)=diffN[cc+4*rr];}}
        
        double Jac[9];
        double coords[12];
        rval = moab.get_coords(tetNodes,coords); CHKERR_PETSC(rval);
        
        ShapeJacMBTET(diffN,coords,Jac);
        
        for (int rr=0; rr<3; rr++) {
            for(int cc=0; cc<3; cc++){
                JacMat(rr,cc)=Jac[cc+3*rr];}}
        
        //        cout<<JacMat(0,0)<<"  "<<JacMat(0,1)<<"  "<<JacMat(0,2)<<endl;
        //        cout<<JacMat(1,0)<<"  "<<JacMat(1,1)<<"  "<<JacMat(1,2)<<endl;
        //        cout<<JacMat(2,0)<<"  "<<JacMat(2,1)<<"  "<<JacMat(2,2)<<endl;
        //        cout<<endl;
        
        ublas::matrix<double> diffN_Jac = prod(trans(diffNMat),JacMat);
        ublas::vector<FieldData> Gradient = prod(diffN_Jac,surf_disp_node);
        
        //        cout<<Gradient[0]<<"  "<<Gradient[1]<<"  "<<Gradient[2]<<endl;
//        double prevGradient[3];
//        Range::iterator niit2 = tetNodes.begin();
//        map<EntityHandle,EntityHandle>::iterator mit2 = node_map.begin();
//        for(;niit2!=tetNodes.end();niit2++,mit2++) {
//            rval = moab_post_proc.tag_get_data(th_gradient,&mit2->second,1,prevGradient); CHKERR_PETSC(rval);  
//            double tempGradient[3];
//            for (int ii=0; ii<3; ii++) {
//                tempGradient[ii] = 0.5*(prevGradient[ii]+Gradient[ii]);
//            }
//            rval = moab_post_proc.tag_set_data(th_gradient,&mit2->second,1,&tempGradient[0]); CHKERR_PETSC(rval);  
//        }
        
        int gg3=0;
        vector< ublas::matrix< FieldData > > surfDistGrad;
        ierr = GetGaussDiffDataVector("SURFACE_DISTANCE",surfDistGrad); CHKERRQ(ierr);        
        vector< ublas::matrix< FieldData > >::iterator viit3 = surfDistGrad.begin();
        map<EntityHandle,EntityHandle>::iterator mit3 = node_map.begin();
        for(;viit3!=surfDistGrad.end();viit3++,mit3++,gg3++) {
            ublas::matrix< FieldData > GradU = *viit3;
            cout << GradU << endl;
//            ublas::vector<FieldData> GradU2 = prod(GradU,surf_disp_node);
            double prevGradient[3];
//            rval = moab_post_proc.tag_get_data(th_gradient,&mit3->second,1,prevGradient); CHKERR_PETSC(rval);  
            double tempGradient[3];
            for (int ii=0; ii<3; ii++) {
//                tempGradient[ii] = 0.5*(prevGradient[ii]+GradU[ii]);
                tempGradient[ii] = GradU(0,ii);
            }
            cout<<tempGradient[0]<<"  "<<tempGradient[1]<<"  "<<tempGradient[2]<<endl;
            rval = moab_post_proc.tag_set_data(th_gradient,&mit3->second,1,&tempGradient[0]); CHKERR_PETSC(rval);  
        }
        cout<<endl;
        
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
            //            StiffnessMatrix.resize(6);
            //            StiffnessMatrix.clear();
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
            
            ///Get Gradients for fibre_direction
            
            //            double surf_disp_gauss_point = g_NTET[4*gg+0]*surf_disp_node[0] + g_NTET[4*gg+1]*surf_disp_node[1] + g_NTET[4*gg+2]*surf_disp_node[2] + g_NTET[4*gg+3]*surf_disp_node[3];
            
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
        
    //Interface meshset 4
    EntityHandle meshset_interface;
    ierr = mField.get_msId_meshset(4,SideSet,meshset_interface); CHKERRQ(ierr);
    ierr = mField.get_msId_3dENTS_sides(meshset_interface,true); CHKERRQ(ierr);

    // stl::bitset see for more details
    BitRefLevel bit_level_interface;
    bit_level_interface.set(0);
    ierr = mField.get_msId_3dENTS_split_sides(0,bit_level_interface,meshset_interface,true,true); CHKERRQ(ierr);
    EntityHandle meshset_level_interface;
    rval = moab.create_meshset(MESHSET_SET,meshset_level_interface); CHKERR_PETSC(rval);
    ierr = mField.refine_get_ents(bit_level_interface,meshset_level_interface); CHKERRQ(ierr);
    
    //update BC for refined (with interface) mesh
    EntityHandle meshset_SideSet1; //Dirihlet BC is there
    ierr = mField.get_msId_meshset(1,SideSet,meshset_SideSet1); CHKERRQ(ierr);
    ierr = mField.refine_get_childern(meshset_SideSet1,bit_level_interface,meshset_SideSet1,MBTRI,true); CHKERRQ(ierr);
    EntityHandle meshset_SideSet2; //Dirihlet BC is there
    ierr = mField.get_msId_meshset(2,SideSet,meshset_SideSet2); CHKERRQ(ierr);
    ierr = mField.refine_get_childern(meshset_SideSet2,bit_level_interface,meshset_SideSet2,MBTRI,true); CHKERRQ(ierr);
    EntityHandle meshset_BlockSet1; //Dirihlet BC is there
    ierr = mField.get_msId_meshset(1,BlockSet,meshset_BlockSet1); CHKERRQ(ierr);
    ierr = mField.refine_get_childern(meshset_BlockSet1,bit_level_interface,meshset_BlockSet1,MBTET,true); CHKERRQ(ierr);
    EntityHandle meshset_BlockSet2; //Dirihlet BC is there
    ierr = mField.get_msId_meshset(2,BlockSet,meshset_BlockSet2); CHKERRQ(ierr);
    ierr = mField.refine_get_childern(meshset_BlockSet2,bit_level_interface,meshset_BlockSet2,MBTET,true); CHKERRQ(ierr);
//    EntityHandle meshset_BlockSet3; //Dirihlet BC is there
//    ierr = mField.get_msId_meshset(3,BlockSet,meshset_BlockSet3); CHKERRQ(ierr);
//    ierr = mField.refine_get_childern(meshset_BlockSet3,bit_level_interface,meshset_BlockSet3,MBTET,true); CHKERRQ(ierr);
    
    //Interface meshset 5
    EntityHandle meshset_interface1;
    ierr = mField.get_msId_meshset(5,SideSet,meshset_interface1); CHKERRQ(ierr);
    ierr = mField.get_msId_3dENTS_sides(meshset_interface1,true); CHKERRQ(ierr);

    BitRefLevel bit_level_interface1;
    bit_level_interface1.set(1);
    ierr = mField.get_msId_3dENTS_split_sides(0,bit_level_interface1,meshset_interface1,true,true); CHKERRQ(ierr);
    EntityHandle meshset_level_interface1;
    rval = moab.create_meshset(MESHSET_SET,meshset_level_interface1); CHKERR_PETSC(rval);
    ierr = mField.refine_get_ents(bit_level_interface1,meshset_level_interface1); CHKERRQ(ierr);
    
    ierr = mField.refine_get_childern(meshset_SideSet1, bit_level_interface1,meshset_SideSet1,MBTRI,true); CHKERRQ(ierr);
    ierr = mField.refine_get_childern(meshset_SideSet2, bit_level_interface1,meshset_SideSet2,MBTRI,true); CHKERRQ(ierr);
    ierr = mField.refine_get_childern(meshset_BlockSet1,bit_level_interface1,meshset_BlockSet1,MBTET,true); CHKERRQ(ierr);
    ierr = mField.refine_get_childern(meshset_BlockSet2,bit_level_interface1,meshset_BlockSet2,MBTET,true); CHKERRQ(ierr);
//    ierr = mField.refine_get_childern(meshset_BlockSet3,bit_level_interface1,meshset_BlockSet3,MBTET,true); CHKERRQ(ierr);

    //Interface meshset 6
    EntityHandle meshset_interface2;
    ierr = mField.get_msId_meshset(6,SideSet,meshset_interface2); CHKERRQ(ierr);
    ierr = mField.get_msId_3dENTS_sides(meshset_interface2,true); CHKERRQ(ierr);
    
    BitRefLevel bit_level_interface2;
    bit_level_interface1.set(2);
    ierr = mField.get_msId_3dENTS_split_sides(0,bit_level_interface2,meshset_interface2,true,true); CHKERRQ(ierr);
    EntityHandle meshset_level_interface2;
    rval = moab.create_meshset(MESHSET_SET,meshset_level_interface2); CHKERR_PETSC(rval);
    ierr = mField.refine_get_ents(bit_level_interface2,meshset_level_interface2); CHKERRQ(ierr);

    ierr = mField.refine_get_childern(meshset_SideSet1, bit_level_interface2,meshset_SideSet1,MBTRI,true); CHKERRQ(ierr);
    ierr = mField.refine_get_childern(meshset_SideSet2, bit_level_interface2,meshset_SideSet2,MBTRI,true); CHKERRQ(ierr);
    ierr = mField.refine_get_childern(meshset_BlockSet1,bit_level_interface2,meshset_BlockSet1,MBTET,true); CHKERRQ(ierr);
    ierr = mField.refine_get_childern(meshset_BlockSet2,bit_level_interface2,meshset_BlockSet2,MBTET,true); CHKERRQ(ierr);
//    ierr = mField.refine_get_childern(meshset_BlockSet3,bit_level_interface2,meshset_BlockSet3,MBTET,true); CHKERRQ(ierr);
    
//    rval = moab.write_file("AAA.vtk","VTK","",&meshset_level_interface,1); CHKERR_PETSC(rval);
    
    // stl::bitset see for more details
    BitRefLevel bit_level0;
    bit_level0.set(3);
    EntityHandle meshset_level0;
    rval = moab.create_meshset(MESHSET_SET,meshset_level0); CHKERR_PETSC(rval);
    ierr = mField.seed_ref_level_3D(meshset_level_interface2,bit_level0); CHKERRQ(ierr);
    ierr = mField.refine_get_ents(bit_level0,meshset_level0); CHKERRQ(ierr);
    
    Range TETSFormLevel0;
    rval = moab.get_entities_by_type(meshset_level0,MBTET,TETSFormLevel0,true); CHKERR_PETSC(rval);
    Range TETsFromBlockSet1onmeshset_level0;
    rval = moab.get_entities_by_type(meshset_BlockSet1,MBTET,TETsFromBlockSet1onmeshset_level0,true); CHKERR_PETSC(rval);
    TETsFromBlockSet1onmeshset_level0 = intersect(TETsFromBlockSet1onmeshset_level0,TETSFormLevel0);
    EntityHandle meshset_BlockSet1OnLevel0;
    rval = moab.create_meshset(MESHSET_SET,meshset_BlockSet1OnLevel0); CHKERR_PETSC(rval); CHKERR_PETSC(rval);
    rval = moab.add_entities(meshset_BlockSet1OnLevel0,TETsFromBlockSet1onmeshset_level0); CHKERR_PETSC(rval);
    
    Range TETsFromBlockSet2onmeshset_level0;
    rval = moab.get_entities_by_type(meshset_BlockSet2,MBTET,TETsFromBlockSet2onmeshset_level0,true); CHKERR_PETSC(rval);
    EntityHandle meshset_BlockSet2OnLevel0;
    rval = moab.create_meshset(MESHSET_SET,meshset_BlockSet2OnLevel0); CHKERR_PETSC(rval); CHKERR_PETSC(rval);
    TETsFromBlockSet2onmeshset_level0 = intersect(TETsFromBlockSet2onmeshset_level0,TETSFormLevel0);
    rval = moab.add_entities(meshset_BlockSet2OnLevel0,TETsFromBlockSet2onmeshset_level0); CHKERR_PETSC(rval);
    
//    Range TETsFromBlockSet3onmeshset_level0;
//    rval = moab.get_entities_by_type(meshset_BlockSet3,MBTET,TETsFromBlockSet3onmeshset_level0,true); CHKERR_PETSC(rval);
//    EntityHandle meshset_BlockSet3OnLevel0;
//    rval = moab.create_meshset(MESHSET_SET,meshset_BlockSet3OnLevel0); CHKERR_PETSC(rval); CHKERR_PETSC(rval);
//    TETsFromBlockSet3onmeshset_level0 = intersect(TETsFromBlockSet3onmeshset_level0,TETSFormLevel0);
//    rval = moab.add_entities(meshset_BlockSet3OnLevel0,TETsFromBlockSet3onmeshset_level0); CHKERR_PETSC(rval);

    //rval = moab.write_file("AAA.vtk","VTK","",&meshset_BlockSet1,1); CHKERR_PETSC(rval);
    
    BitRefLevel problem_bit_level = bit_level0;
    
    /***/
    //Define problem
    
    //Fields
    ierr = mField.add_field("DISPLACEMENT",H1,3); CHKERRQ(ierr);
    ierr = mField.add_field("SURFACE_DISTANCE",H1,1); CHKERRQ(ierr);

    //FE
    ierr = mField.add_finite_element("ELASTIC"); CHKERRQ(ierr);
    ierr = mField.add_finite_element("TRAN_ISOTROPIC_ELASTIC"); CHKERRQ(ierr);
    ierr = mField.add_finite_element("INTERFACE"); CHKERRQ(ierr);

    //Define rows/cols and element data
    ierr = mField.modify_finite_element_add_field_row("ELASTIC","DISPLACEMENT"); CHKERRQ(ierr);
    ierr = mField.modify_finite_element_add_field_col("ELASTIC","DISPLACEMENT"); CHKERRQ(ierr);
    ierr = mField.modify_finite_element_add_field_data("ELASTIC","DISPLACEMENT"); CHKERRQ(ierr);
    //FE Transverse Isotropic
    ierr = mField.modify_finite_element_add_field_row("TRAN_ISOTROPIC_ELASTIC","DISPLACEMENT"); CHKERRQ(ierr);
    ierr = mField.modify_finite_element_add_field_col("TRAN_ISOTROPIC_ELASTIC","DISPLACEMENT"); CHKERRQ(ierr);
    ierr = mField.modify_finite_element_add_field_data("TRAN_ISOTROPIC_ELASTIC","DISPLACEMENT"); CHKERRQ(ierr);
    ierr = mField.modify_finite_element_add_field_data("TRAN_ISOTROPIC_ELASTIC","SURFACE_DISTANCE"); CHKERRQ(ierr);
    //FE Interface
    ierr = mField.modify_finite_element_add_field_row("INTERFACE","DISPLACEMENT"); CHKERRQ(ierr);
    ierr = mField.modify_finite_element_add_field_col("INTERFACE","DISPLACEMENT"); CHKERRQ(ierr);
    ierr = mField.modify_finite_element_add_field_data("INTERFACE","DISPLACEMENT"); CHKERRQ(ierr);
    
    //define problems
    ierr = mField.add_problem("ELASTIC_MECHANICS"); CHKERRQ(ierr);
    
    //set finite elements for problem
    ierr = mField.modify_problem_add_finite_element("ELASTIC_MECHANICS","ELASTIC"); CHKERRQ(ierr);
    ierr = mField.modify_problem_add_finite_element("ELASTIC_MECHANICS","TRAN_ISOTROPIC_ELASTIC"); CHKERRQ(ierr);
    ierr = mField.modify_problem_add_finite_element("ELASTIC_MECHANICS","INTERFACE"); CHKERRQ(ierr);

    
    //set refinment level for problem
    ierr = mField.modify_problem_ref_level_add_bit("ELASTIC_MECHANICS",problem_bit_level); CHKERRQ(ierr);
    
    /***/
    //Declare problem
    
    //add entitities (by tets) to the field
    ierr = mField.add_ents_to_field_by_TETs(0,"DISPLACEMENT"); CHKERRQ(ierr);
    ierr = mField.add_ents_to_field_by_TETs(meshset_BlockSet2OnLevel0,"SURFACE_DISTANCE"); CHKERRQ(ierr);

    //add finite elements entities
    ierr = mField.add_ents_to_finite_element_by_TETs(meshset_BlockSet2OnLevel0,"TRAN_ISOTROPIC_ELASTIC",true); CHKERRQ(ierr);
    ierr = mField.add_ents_to_finite_element_by_TETs(meshset_BlockSet1OnLevel0,"ELASTIC",true); CHKERRQ(ierr);
//    ierr = mField.add_ents_to_finite_element_by_TETs(meshset_BlockSet3OnLevel0,"TRAN_ISOTROPIC_ELASTIC",true); CHKERRQ(ierr);
    //    ierr = mField.add_ents_to_finite_element_EntType_by_bit_ref(problem_bit_level,"ELASTIC",MBTET); CHKERRQ(ierr);
    ierr = mField.add_ents_to_finite_element_EntType_by_bit_ref(problem_bit_level,"INTERFACE",MBPRISM); CHKERRQ(ierr);
    
    //set app. order
    //see Hierarchic Finite Element Bases on Unstructured Tetrahedral Meshes (Mark Ainsworth & Joe Coyle)
    ierr = mField.set_field_order(0,MBTET,"DISPLACEMENT",order); CHKERRQ(ierr);
    ierr = mField.set_field_order(0,MBTRI,"DISPLACEMENT",order); CHKERRQ(ierr);
    ierr = mField.set_field_order(0,MBEDGE,"DISPLACEMENT",order); CHKERRQ(ierr);
    ierr = mField.set_field_order(0,MBVERTEX,"DISPLACEMENT",1); CHKERRQ(ierr);
    
    ierr = mField.set_field_order(0,MBTET,"SURFACE_DISTANCE",1); CHKERRQ(ierr);
    ierr = mField.set_field_order(0,MBTRI,"SURFACE_DISTANCE",1); CHKERRQ(ierr);
    ierr = mField.set_field_order(0,MBEDGE,"SURFACE_DISTANCE",1); CHKERRQ(ierr);
    ierr = mField.set_field_order(0,MBVERTEX,"SURFACE_DISTANCE",1); CHKERRQ(ierr);
    
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
    Range SideSet1,SideSet2,SideSet3,SideSet4,SideSet5,SideSet6,SurfaceMesh;
    ierr = mField.get_Cubit_msId_entities_by_dimension(1,SideSet,2,SideSet1,true); CHKERRQ(ierr);
    ierr = mField.get_Cubit_msId_entities_by_dimension(2,SideSet,2,SideSet2,true); CHKERRQ(ierr);
    ierr = mField.get_Cubit_msId_entities_by_dimension(3,SideSet,2,SideSet3,true); CHKERRQ(ierr);
    ierr = mField.get_Cubit_msId_entities_by_dimension(4,SideSet,2,SideSet4,true); CHKERRQ(ierr);
    ierr = mField.get_Cubit_msId_entities_by_dimension(5,SideSet,2,SideSet5,true); CHKERRQ(ierr);
    ierr = mField.get_Cubit_msId_entities_by_dimension(6,SideSet,2,SideSet6,true); CHKERRQ(ierr);
    ierr = mField.get_Cubit_msId_entities_by_dimension(13,SideSet,2,SurfaceMesh,true); CHKERRQ(ierr);

    PetscPrintf(PETSC_COMM_WORLD,"Nb. faces in SideSet 1 : %u\n",SideSet1.size());
    PetscPrintf(PETSC_COMM_WORLD,"Nb. faces in SideSet 2 : %u\n",SideSet2.size());
    PetscPrintf(PETSC_COMM_WORLD,"Nb. faces in SideSet 3 : %u\n",SideSet3.size());
    PetscPrintf(PETSC_COMM_WORLD,"Nb. faces in SideSet 4 : %u\n",SideSet4.size());
    PetscPrintf(PETSC_COMM_WORLD,"Nb. faces in SideSet 5 : %u\n",SideSet5.size());
    PetscPrintf(PETSC_COMM_WORLD,"Nb. faces in SideSet 6 : %u\n",SideSet6.size());
    
    //Assemble F and Aij
    const double YoungModulusP = 135000;
    const double PoissonRatioP = 0.77;
    const double YoungModulusZ = 135000;
    const double PoissonRatioPZ = 0.2;
    const double ShearModulusZP = 5000;
    const double YoungModulus = 200000;
    const double PoissonRatio = 0.3;
    const double alpha = 0.05;
    
    //*********************** Saving Shortest distance from fibre outer surface to SURFACE_DISTANCE field ************************//
    
    Tag th_dist_from_surface;
    double def_VAL3=0.0;
    rval = moab.tag_get_handle("SURFACE_DISTANCE",1,MB_TYPE_DOUBLE,th_dist_from_surface,MB_TAG_CREAT|MB_TAG_SPARSE,&def_VAL3); CHKERR_THROW(rval);
    
    EntityHandle tree_root;
    AdaptiveKDTree tool( &moab );
    rval = tool.build_tree( SurfaceMesh, tree_root ); CHKERR(rval);
    
    EntityHandle node = 0;
    double coords[3];
    for(_IT_GET_DOFS_MOABFIELD_BY_NAME_FOR_LOOP_(mField,"SURFACE_DISTANCE",dof_ptr)) {
        if(dof_ptr->get_ent_type()!=MBVERTEX) continue;
        EntityHandle ent = dof_ptr->get_ent();
        double &fval = dof_ptr->get_FieldData();
        double distance2;
//        int dof_rank = dof_ptr->get_dof_rank();
        if(node!=ent) {
            rval = moab.get_coords(&ent,1,coords); CHKERR_PETSC(rval);
            node = ent;
            double closestPoint1[3];
            EntityHandle closestTri1;
            rval = tool.closest_triangle(tree_root,coords,closestPoint1,closestTri1);CHKERR(rval);
            distance2 = sqrt(pow(closestPoint1[0]-coords[0],2)+pow(closestPoint1[1]-coords[1],2)+pow(closestPoint1[2]-coords[2],2));
//            cout<<distance2<<endl;
            rval = moab.tag_set_data(th_dist_from_surface,&ent,1,&distance2); CHKERR_PETSC(rval);  
        }
        fval = distance2;
    }
    
    
    //**************************** Function to Set Boundary Conditions ******************************//  
    
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
    TranIsotropicElasticFEMethod MyTIsotFE(moab,&myDirihletBC,Aij,F,LAMBDA(YoungModulusP,PoissonRatioP),MU(YoungModulusP,PoissonRatioP),YoungModulusP,YoungModulusZ,PoissonRatioP,PoissonRatioPZ,ShearModulusZP,SideSet1,SideSet2);
    MyElasticFEMethod MyFE(moab,&myDirihletBC,Aij,F,LAMBDA(YoungModulus,PoissonRatio),MU(YoungModulus,PoissonRatio),SideSet1,SideSet2);
    Range dummy;
    InterfaceFEMethod IntMyFE(moab,&myDirihletBC,Aij,F,YoungModulus*alpha,SideSet1,SideSet2,dummy);

    ierr = VecZeroEntries(F); CHKERRQ(ierr);
    ierr = VecGhostUpdateBegin(F,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
    ierr = VecGhostUpdateEnd(F,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
    ierr = MatZeroEntries(Aij); CHKERRQ(ierr);
    
    ierr = mField.loop_finite_elements("ELASTIC_MECHANICS","INTERFACE",IntMyFE);  CHKERRQ(ierr);
    PetscSynchronizedFlush(PETSC_COMM_WORLD);
    
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
        ierr = mField.problem_get_FE("ELASTIC_MECHANICS","INTERFACE",out_meshset); CHKERRQ(ierr);
        rval = moab.write_file(outName,"VTK","",&out_meshset,1); CHKERR_PETSC(rval);
        rval = moab.delete_entities(&out_meshset,1); CHKERR_PETSC(rval);
    }
    
    //  PostProcDisplacemenysAndStarinOnRefMesh fe_post_proc_method(moab);
    PostProcDisplacemenysAndStarinAndElasticLinearStressOnRefMesh fe_post_proc_method(moab,LAMBDA(YoungModulus,PoissonRatio),MU(YoungModulus,PoissonRatio));
    TranIsotropic_Fibre_PostProcDisplacemenysAndStarinAndElasticLinearStressOnRefMesh fe_fibre_post_proc_method( moab, LAMBDA(YoungModulusP,PoissonRatioP),MU(YoungModulusP,PoissonRatioP), YoungModulusP,YoungModulusZ,PoissonRatioP,PoissonRatioPZ,ShearModulusZP);
    
    ierr = mField.loop_finite_elements("ELASTIC_MECHANICS","ELASTIC",fe_post_proc_method);  CHKERRQ(ierr);
    ierr = mField.loop_finite_elements("ELASTIC_MECHANICS","TRAN_ISOTROPIC_ELASTIC",fe_fibre_post_proc_method);  CHKERRQ(ierr);
    
    PetscSynchronizedFlush(PETSC_COMM_WORLD);
    if(pcomm->rank()==0) {
        rval = fe_fibre_post_proc_method.moab_post_proc.write_file(outName2,"VTK",""); CHKERR_PETSC(rval);
    }

    PostProcCohesiveForces fe_post_proc_prisms(moab,YoungModulus*alpha);
    ierr = mField.loop_finite_elements("ELASTIC_MECHANICS","INTERFACE",fe_post_proc_prisms);  CHKERRQ(ierr);
    PetscSynchronizedFlush(PETSC_COMM_WORLD);
    if(pcomm->rank()==0) {
        rval = fe_post_proc_prisms.moab_post_proc.write_file("out_post_proc_prisms.vtk","VTK",""); CHKERR_PETSC(rval);
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
