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

#ifndef __FE2_RHS_R_PSFEM_HPP__
#define __FE2_RHS_R_PSFEM_HPP__

using namespace ObosleteUsersModules;

namespace MoFEM {
  
  struct FE2_Rhs_r_PSFEM: public FE2_ElasticFEMethod {
    
    ublas::matrix<FieldData> Dmat_r;
    Vec dF;
    string zeroth_field;
    FE2_Rhs_r_PSFEM(FieldInterface& _mField,
                    Mat &_A,
                    Vec _D,
                    Vec _F,
                    ublas::matrix<FieldData> _Dmat,
                    string _field_name):
    FE2_ElasticFEMethod(_mField,
                      _A,
                      _D,
                      _F,
                      _Dmat,
                      _field_name),dF(_F),Dmat_r(_Dmat),zeroth_field(_field_name){};
    
	ublas::matrix<ublas::matrix<FieldData> > K_r;
	virtual PetscErrorCode Stiffness_K_r() {
      PetscFunctionBegin;
    K_r.resize(row_mat,col_mat);//cout<<row_mat<<col_mat<<endl;
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
          cblas_dsymm(CblasRowMajor,CblasLeft,CblasUpper,
                      BD.size1(),BD.size2(),
                      w,&*Dmat_r.data().begin(),Dmat_r.size2(),
                      &*row_Mat.data().begin(),row_Mat.size2(),
                      0.,&*BD.data().begin(),BD.size2());
          for(int cc = rr;cc<col_mat;cc++) {
            if(ColGlob[cc].size()==0) continue;
            ublas::matrix<FieldData> &col_Mat = (colBMatrices[cc])[gg];
            if(gg == 0) {
              K_r(rr,cc).resize(BD.size2(),col_Mat.size2());
              cblas_dgemm(CblasRowMajor,CblasTrans,CblasNoTrans,
                          BD.size2(),col_Mat.size2(),BD.size1(),
                          1.,&*BD.data().begin(),BD.size2(),
                          &*col_Mat.data().begin(),col_Mat.size2(),
                          0.,&*K_r(rr,cc).data().begin(),K_r(rr,cc).size2());
            } else {
              cblas_dgemm(CblasRowMajor,CblasTrans,CblasNoTrans,
                          BD.size2(),col_Mat.size2(),BD.size1(),
                          1.,&*BD.data().begin(),BD.size2(),
                          &*col_Mat.data().begin(),col_Mat.size2(),
                          1.,&*K_r(rr,cc).data().begin(),K_r(rr,cc).size2());
            }
          }
        }
      }
      PetscFunctionReturn(0);
    }

  
    vector<ublas::vector<FieldData> > f_re;
    virtual PetscErrorCode Rhs() {
      PetscFunctionBegin;
      ierr = Stiffness_K_r(); CHKERRQ(ierr);
      
      vector<ublas::vector<FieldData> > D_elm;
      
      D_elm.resize(col_mat);
      int col_mat1 = 0;  //nodes
      ierr = GetDataVector(zeroth_field,D_elm[col_mat1]); CHKERRQ(ierr);
      col_mat1++;
      
      
      for(int ee=0; ee<6; ee++) { //edges
        if(ColGlob[col_mat1].size()!=0) {
          ierr = GetDataVector(zeroth_field,MBEDGE,D_elm[col_mat1],ee); CHKERRQ(ierr);
          col_mat1++;
        }
      }
      
      for(int ff=0; ff<4; ff++) { //faces
        if(ColGlob[col_mat1].size()!=0) {
          ierr = GetDataVector(zeroth_field,MBTRI,D_elm[col_mat1],ff); CHKERRQ(ierr);
          col_mat1++;
        }
      }
      
      if(ColGlob[col_mat1].size()!=0) {
        ierr = GetDataVector(zeroth_field,MBTET,D_elm[col_mat1]); CHKERRQ(ierr);
      }
      
      f_re.resize(row_mat);
      for(int rr = 0;rr<row_mat;rr++) {
        if(RowGlob[rr].size()==0) continue;
        
        int rr_start=0;
        for(int cc = 0;cc<col_mat;cc++) {
          if(ColGlob[cc].size()==0) continue;
          if(rr_start == 0) {
            f_re[rr] =  -prod(K_r(rr,cc),D_elm[cc]);
            rr_start++;
          } else {
            f_re[rr] -= prod(K_r(rr,cc),D_elm[cc]);
          }
        }
      }
      
      for(int rr = 0;rr<row_mat;rr++) {
        if(RowGlob[rr].size()==0) continue;
        if(RowGlob[rr].size()!=f_re[rr].size()) SETERRQ(PETSC_COMM_SELF,1,"data inconsistency");
        ierr = VecSetValues(dF,RowGlob[rr].size(),&(RowGlob[rr])[0],&(f_re[rr].data()[0]),ADD_VALUES); CHKERRQ(ierr);
      }
      
      PetscFunctionReturn(0);
    }
    
    
    PetscErrorCode postProcess() {
      PetscFunctionBegin;
      PetscFunctionReturn(0);
    }

    PetscErrorCode operator()() {
      PetscFunctionBegin;
      
      ierr = Get_g_NTET(); CHKERRQ(ierr);
      ierr = OpStudentStart_TET(g_NTET); CHKERRQ(ierr);
      // cout<<"Hi from K_rPoissonP_ElasticFEMethodTransIso "<<endl;
      ierr = GetMatrices(); CHKERRQ(ierr);
      ierr = Rhs(); CHKERRQ(ierr);
//      std::string wait;
//      std::cin >> wait;
      ierr = OpStudentEnd(); CHKERRQ(ierr);
      PetscFunctionReturn(0);
    }

  };
  
  
}

#endif //__FE2_RHS_R_PSFEM_HPP__
