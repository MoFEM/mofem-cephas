/* Copyright (C) 2013, Zahur Ullah (Zahur.Ullah AT glasgow.ac.uk)
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

#ifndef __ElasticFE_RVELagrange_Homogenized_Stress_Disp_HPP__
#define __ElasticFE_RVELagrange_Homogenized_Stress_Disp__

#include <boost/numeric/ublas/symmetric.hpp>
#include "ElasticFE_RVELagrange_Disp.hpp"

namespace MoFEM {

struct ElasticFE_RVELagrange_Homogenized_Stress_Disp: public ElasticFE_RVELagrange_Disp {

    Vec DVec;
    double *RVE_volume;

    ElasticFE_RVELagrange_Homogenized_Stress_Disp(FieldInterface& _mField,BaseDirihletBC *_dirihlet_ptr,Mat &_Aij,Vec &_D,Vec& _F,double *_RVE_volume):
    ElasticFE_RVELagrange_Disp(_mField, _dirihlet_ptr,_Aij, _D, _F), DVec(_D),RVE_volume(_RVE_volume){};
    
    
    
    PetscErrorCode preProcess() {
        PetscFunctionBegin;
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"Start Assembly\n");
        PetscSynchronizedFlush(PETSC_COMM_WORLD);
        ierr = PetscTime(&v1); CHKERRQ(ierr);
        ierr = PetscGetCPUTime(&t1); CHKERRQ(ierr);
        
        Stress_Homo.resize(6);   Stress_Homo.clear();
        
        PetscFunctionReturn(0);
    }

    
    PetscErrorCode postProcess() {
        PetscFunctionBegin;
        // Note MAT_FLUSH_ASSEMBLY
        ierr = MatAssemblyBegin(Aij,MAT_FLUSH_ASSEMBLY); CHKERRQ(ierr);
        ierr = MatAssemblyEnd(Aij,MAT_FLUSH_ASSEMBLY); CHKERRQ(ierr);
        ierr = VecAssemblyBegin(F); CHKERRQ(ierr);
        ierr = VecAssemblyEnd(F); CHKERRQ(ierr);
        ierr = VecAssemblyBegin(F); CHKERRQ(ierr);
        ierr = VecAssemblyEnd(F); CHKERRQ(ierr);
        ierr = PetscTime(&v2); CHKERRQ(ierr);
        ierr = PetscGetCPUTime(&t2); CHKERRQ(ierr);
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"End Assembly: Rank %d Time = %f CPU Time = %f\n",pcomm->rank(),v2-v1,t2-t1);
        PetscSynchronizedFlush(PETSC_COMM_WORLD);
        
        cout<<"RVE_volume in postProcess = "<<*RVE_volume<<endl;
        Stress_Homo=(1.0/(*RVE_volume))*Stress_Homo;
        
        cout<< " Stress_Homo =  "<<endl;
        for(int ii=0; ii<6; ii++) cout<<Stress_Homo(ii)<<endl; 
        
        
        PetscFunctionReturn(0);
    }
 
    
    
    ublas::matrix<FieldData> X_mat, nodes_coord, gauss_coord;
    ublas::vector<ublas::matrix<FieldData> > D_mat;
    ublas::vector<FieldData>  Stress_Homo;
    ublas::vector<ublas::vector<FieldData> > Lamda;
    
    
    virtual PetscErrorCode Calculate_Homo_Stress() {
        PetscFunctionBegin;
        X_mat.resize(3,6);    X_mat.clear();
        nodes_coord.resize(3,3);
        gauss_coord.resize(3,g_TRI_dim);
        D_mat.resize(row_mat);
        Lamda.resize(row_mat);

        //used to calculate the coordinates of a Gauss points
        nodes_coord(0,0)=coords_face[0]; nodes_coord(0,1)=coords_face[3]; nodes_coord(0,2)=coords_face[6];
        nodes_coord(1,0)=coords_face[1]; nodes_coord(1,1)=coords_face[4]; nodes_coord(1,2)=coords_face[7];
        nodes_coord(2,0)=coords_face[2]; nodes_coord(2,1)=coords_face[5]; nodes_coord(2,2)=coords_face[8];
        
        //coordinates for all gauss points
        gauss_coord=prod(nodes_coord, g_NTRI_mat);
        
//        cout<<"g_NTRI_mat "<<g_NTRI_mat<<endl<<endl;
//        cout<<"nodes_coord "<<nodes_coord<<endl<<endl;
//        cout<<"gauss_coord "<<gauss_coord<<endl<<endl;
//        std::string wait;
//        std::cin >> wait;

        for(int rr=0; rr<row_mat; rr++){
            for(int gg = 0;gg<g_TRI_dim;gg++) {
                double w = area*G_W_TRI[gg];
                X_mat(0,0)=2.0*gauss_coord(0,gg);  X_mat(0,3)=gauss_coord(1,gg);  X_mat(0,4)=gauss_coord(2,gg);
                X_mat(1,1)=2.0*gauss_coord(1,gg);  X_mat(1,3)=gauss_coord(0,gg);  X_mat(1,5)=gauss_coord(2,gg);
                X_mat(2,2)=2.0*gauss_coord(2,gg);  X_mat(2,4)=gauss_coord(0,gg);  X_mat(2,5)=gauss_coord(1,gg);
                X_mat=0.5*X_mat;
                
                ublas::matrix<FieldData> &row_Mat = (rowNMatrices[rr])[gg];
                ublas::matrix<FieldData> &col_Mat = X_mat;
                
                ublas::matrix<FieldData>  D_mat1;    //Dmat1=NT*X_mat
                D_mat1.resize(row_Mat.size2(),col_Mat.size2());
                
                //Integrate D_mat
                if(gg == 0) {
                    D_mat[rr].resize(H_mat[rr].size1(),D_mat1.size2());
//                    cout<<"\n row_Mat "<<row_Mat<<endl;
//                    cout<<"\n col_Mat "<<rr<<col_Mat;
//                    cout<<"\n w "<<w<<col_Mat;
                    D_mat1=prod(w*trans(row_Mat), col_Mat);
                    D_mat[rr]=prod(H_mat[rr], D_mat1);
                } else {
                    //calculate (D_mat1= w * NT * X_mat)
                    D_mat1=prod(w*trans(row_Mat), col_Mat);
                    //calculate (D_mat = H_mat * D_mat1)
                    D_mat[rr]+=prod(H_mat[rr], D_mat1);
                }
            }
//            cout<<"row_mat  =  "<<row_mat<<endl;
//            cout<< " D_mat[rr] =  "<<D_mat[rr]<<endl;
            Lamda[rr].resize(RowGlob[rr].size());
            ierr = VecGetValues(DVec,RowGlob[rr].size(),&(RowGlob[rr])[0],&(Lamda[rr].data())[0]); CHKERRQ(ierr);
//            cout<< " Lamda[rr] =  "<<Lamda[rr]<<endl;
            
//            cout<< " Stress_Homo before =  "<<Stress_Homo<<endl;
//            cout<< "  New part          =  "<<prod(trans(D_mat[rr]), Lamda[rr])<<endl;
            Stress_Homo+=prod(trans(D_mat[rr]), -1*Lamda[rr]);   //Lamda is reaction force (so multiply for -1 to get the force)
//            cout<< " Stress_Homo after  =   "<<Stress_Homo<<endl;
        }
        
        PetscFunctionReturn(0);
    }

    
    PetscErrorCode operator()() {
        PetscFunctionBegin;
//        cout<<"Hi from class ElasticFE_RVELagrange_Homogenized_Stress"<<endl;
        ierr = GetN_and_Indices(); CHKERRQ(ierr);
        ierr = Get_H_mat();   //It will be used from the Class ElasticFE_RVELagrange_Disp
        ierr = Calculate_Homo_Stress(); CHKERRQ(ierr);
        PetscFunctionReturn(0);
    }

    
};

    
}

#endif //__ElasticFE_RVELagrange_Periodic_RigidBodyMotion
