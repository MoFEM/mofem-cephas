/* Copyright (C) 2014, Zahur Ullah (Zahur.Ullah AT glasgow.ac.uk)
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


#include <MoFEM.hpp>
using namespace MoFEM;
#include <FEMethod_LowLevelStudent.hpp>
#include <FEMethod_UpLevelStudent.hpp>
#include <ElasticFEMethod.hpp>

#include <ElasticFE_RVELagrange_Disp.hpp>
#include <ElasticFE_RVELagrange_Disp_Multi_Rhs.hpp>

using namespace ObosleteUsersModules;

namespace MoFEM {
  
  //********************************************************************************
  //Post-process Force vectors and stiffness matrix 
  PetscErrorCode ElasticFE_RVELagrange_Disp_Multi_Rhs::postProcess() {
    PetscFunctionBegin;
    
//    cout<<"Hi from  ElasticFE_RVELagrange_Disp_Multi_Rhs::postProcess "<<endl;
    switch(snes_ctx) {
      case CTX_SNESNONE: {
        ierr = MatAssemblyBegin(snes_B,MAT_FLUSH_ASSEMBLY); CHKERRQ(ierr);
        ierr = MatAssemblyEnd(snes_B,MAT_FLUSH_ASSEMBLY); CHKERRQ(ierr);
        ierr = VecAssemblyBegin(F1); CHKERRQ(ierr);
        ierr = VecAssemblyEnd(F1); CHKERRQ(ierr);
        
        ierr = VecAssemblyBegin(F2); CHKERRQ(ierr);
        ierr = VecAssemblyEnd(F2); CHKERRQ(ierr);
        
        ierr = VecAssemblyBegin(F3); CHKERRQ(ierr);
        ierr = VecAssemblyEnd(F3); CHKERRQ(ierr);

        if(rank_field==3){ //This poriton will execute for mechanical problem (rank=3) only
          ierr = VecAssemblyBegin(F4); CHKERRQ(ierr);
          ierr = VecAssemblyEnd(F4); CHKERRQ(ierr);
          
          ierr = VecAssemblyBegin(F5); CHKERRQ(ierr);
          ierr = VecAssemblyEnd(F5); CHKERRQ(ierr);
          
          ierr = VecAssemblyBegin(F6); CHKERRQ(ierr);
          ierr = VecAssemblyEnd(F6); CHKERRQ(ierr);
        }

      }
        break;
      case CTX_SNESSETFUNCTION: {
        ierr = VecAssemblyBegin(snes_f); CHKERRQ(ierr);
        ierr = VecAssemblyEnd(snes_f); CHKERRQ(ierr);
      }
        break;
      case CTX_SNESSETJACOBIAN: {
        // Note MAT_FLUSH_ASSEMBLY
        ierr = MatAssemblyBegin(snes_B,MAT_FLUSH_ASSEMBLY); CHKERRQ(ierr);
        ierr = MatAssemblyEnd(snes_B,MAT_FLUSH_ASSEMBLY); CHKERRQ(ierr);
      }
        break;
      default:
        SETERRQ(PETSC_COMM_SELF,1,"not implemented");
    }
//    cout<<"End of  ElasticFE_RVELagrange_Disp_Multi_Rhs::postProcess "<<endl;

    PetscFunctionReturn(0);
  }
  
  //********************************************************************************
  //Calculate the right hand side vector, i.e. f=D_max * applied_strain and assemble it into the global force vector F
  PetscErrorCode ElasticFE_RVELagrange_Disp_Multi_Rhs::Rhs() {
//    cout<<"Hi from  ElasticFE_RVELagrange_Disp_Multi_Rhs::Rhs "<<endl;

    PetscFunctionBegin;
    X_mat.resize(rank_field,1.5*rank_field+1.5);    X_mat.clear();
    nodes_coord.resize(3,3);
    gauss_coord.resize(3,g_TRI_dim);
    D_mat.resize(row_mat);
    
    //used to calculate the coordinates of a Gauss points
    nodes_coord(0,0)=coords_face[0]; nodes_coord(0,1)=coords_face[3]; nodes_coord(0,2)=coords_face[6];
    nodes_coord(1,0)=coords_face[1]; nodes_coord(1,1)=coords_face[4]; nodes_coord(1,2)=coords_face[7];
    nodes_coord(2,0)=coords_face[2]; nodes_coord(2,1)=coords_face[5]; nodes_coord(2,2)=coords_face[8];
    
    //coordinates for all gauss points
    gauss_coord=prod(nodes_coord, g_NTRI_mat);
    
    //        cout<<"g_NTRI_mat "<<g_NTRI_mat<<endl<<endl;
    //        cout<<"nodes_coord "<<nodes_coord<<endl<<endl;
    
    //cout<<"area "<<area << endl;
    for(int rr=0; rr<row_mat; rr++){
      for(int gg = 0;gg<g_TRI_dim;gg++) {
        double w = area*G_W_TRI[gg];
        
        switch(rank_field) {
          case 3:
            X_mat(0,0)=2.0*gauss_coord(0,gg);  X_mat(0,3)=gauss_coord(1,gg);  X_mat(0,4)=gauss_coord(2,gg);
            X_mat(1,1)=2.0*gauss_coord(1,gg);  X_mat(1,3)=gauss_coord(0,gg);  X_mat(1,5)=gauss_coord(2,gg);
            X_mat(2,2)=2.0*gauss_coord(2,gg);  X_mat(2,4)=gauss_coord(0,gg);  X_mat(2,5)=gauss_coord(1,gg);
            X_mat=0.5*X_mat;
            break;
          case 1:
            X_mat(0,0)=gauss_coord(0,gg);  X_mat(0,1)=gauss_coord(1,gg);  X_mat(0,2)=gauss_coord(2,gg);
            break;
          default:
            SETERRQ(PETSC_COMM_SELF,1,"not implemented");
        }
        
        
        ublas::matrix<FieldData> &row_Mat = (rowNMatrices[rr])[gg];
        ublas::matrix<FieldData> &col_Mat = X_mat;
        
        ublas::matrix<FieldData>  D_mat1;    //Dmat1=NT*X_mat
        D_mat1.resize(row_Mat.size2(),col_Mat.size2());
        
        //Integrate D_mat
        if(gg == 0) {
          D_mat[rr].resize(row_Mat.size2(),col_Mat.size2());
          
          //calculate (D_mat1= w * NT * X_mat)
          cblas_dgemm(CblasRowMajor,CblasTrans,CblasNoTrans,
                      row_Mat.size2(),col_Mat.size2(),row_Mat.size1(),
                      w,&*row_Mat.data().begin(),row_Mat.size2(),
                      &*col_Mat.data().begin(),col_Mat.size2(),
                      0.,&*D_mat1.data().begin(),D_mat1.size2());
          
          //calculate (D_mat = H_mat * D_mat1)
          cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,
                      H_mat[rr].size2(),D_mat1.size2(),H_mat[rr].size1(),
                      1.0,&*H_mat[rr].data().begin(),H_mat[rr].size2(),
                      &*D_mat1.data().begin(),D_mat1.size2(),
                      0.,&*D_mat[rr].data().begin(),D_mat[rr].size2());
          
        } else {
          
          //calculate (D_mat1= w * NT * X_mat)
          cblas_dgemm(CblasRowMajor,CblasTrans,CblasNoTrans,
                      row_Mat.size2(),col_Mat.size2(),row_Mat.size1(),
                      w,&*row_Mat.data().begin(),row_Mat.size2(),
                      &*col_Mat.data().begin(),col_Mat.size2(),
                      0.,&*D_mat1.data().begin(),D_mat1.size2());
          
          //calculate (D_mat = H_mat * D_mat1)
          cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,
                      H_mat[rr].size2(),D_mat1.size2(),H_mat[rr].size1(),
                      1.0,&*H_mat[rr].data().begin(),H_mat[rr].size2(),
                      &*D_mat1.data().begin(),D_mat1.size2(),
                      1.,&*D_mat[rr].data().begin(),D_mat[rr].size2());
        }
      }
      //cout<< " D_mat[rr] =  "<<D_mat[rr]<<endl<<endl;
      
      applied_strain.clear();
      applied_strain(0)=1.0;
//      cout<<"applied_strain = "<<applied_strain<<endl;
      f=prod(D_mat[rr], applied_strain);
      ierr = VecSetValues(F1,RowGlob[rr].size(),&(RowGlob[rr])[0],&(f.data())[0],ADD_VALUES); CHKERRQ(ierr);
//      cout<<"Hello "<<endl;
//      std::string wait;
//      std::cin >> wait;
      applied_strain.clear();
      applied_strain(1)=1.0;
//      cout<<"applied_strain = "<<applied_strain<<endl;
      f=prod(D_mat[rr], applied_strain);
//      cout<<"f = "<<f<<endl;
      ierr = VecSetValues(F2,RowGlob[rr].size(),&(RowGlob[rr])[0],&(f.data())[0],ADD_VALUES); CHKERRQ(ierr);

      applied_strain.clear();
      applied_strain(2)=1.0;
//      cout<<"applied_strain = "<<applied_strain<<endl;
      f=prod(D_mat[rr], applied_strain);
//      cout<<"f = "<<f<<endl;
      ierr = VecSetValues(F3,RowGlob[rr].size(),&(RowGlob[rr])[0],&(f.data())[0],ADD_VALUES); CHKERRQ(ierr);

      if(rank_field==3){ //This poriton will execute for mechanical problem (rank=3) only
        applied_strain.clear();
        applied_strain(3)=1.0;
        //      cout<<"applied_strain = "<<applied_strain<<endl;
        f=prod(D_mat[rr], applied_strain);
        //      cout<<"f = "<<f<<endl;
        ierr = VecSetValues(F4,RowGlob[rr].size(),&(RowGlob[rr])[0],&(f.data())[0],ADD_VALUES); CHKERRQ(ierr);
        
        applied_strain.clear();
        applied_strain(4)=1.0;
        //      cout<<"applied_strain = "<<applied_strain<<endl;
        f=prod(D_mat[rr], applied_strain);
        //      cout<<"f = "<<f<<endl;
        ierr = VecSetValues(F5,RowGlob[rr].size(),&(RowGlob[rr])[0],&(f.data())[0],ADD_VALUES); CHKERRQ(ierr);
        
        applied_strain.clear();
        applied_strain(5)=1.0;
        //      cout<<"applied_strain = "<<applied_strain<<endl;
        f=prod(D_mat[rr], applied_strain);
        //      cout<<"f = "<<f<<endl;
        ierr = VecSetValues(F6,RowGlob[rr].size(),&(RowGlob[rr])[0],&(f.data())[0],ADD_VALUES); CHKERRQ(ierr);
      }

    }
//    cout<<"end of  ElasticFE_RVELagrange_Disp_Multi_Rhs::Rhs "<<endl;
    PetscFunctionReturn(0);
  }
  
}

