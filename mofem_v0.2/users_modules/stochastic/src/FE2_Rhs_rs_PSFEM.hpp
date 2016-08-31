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

#ifndef __FE2_RHS_RS_PSFEM_HPP__
#define __FE2_RHS_RS_PSFEM_HPP__

using namespace ObosleteUsersModules;

namespace MoFEM {
  
  struct FE2_Rhs_rr_PSFEM: public FE2_ElasticFEMethod {
    
    ublas::matrix<FieldData> Dmat_r; // 1st-order derivative of D matrix w.r.t. the r-th variable
    ublas::matrix<FieldData> Dmat_s; // 1st-order derivative of D matrix w.r.t. the s-th variable
    ublas::matrix<FieldData> Dmat_rs; // 2nd-order derivative of D matrix w.r.t. the r- and s-th variables
    Vec ddF;
    string zeroth_field;
    string first_field_r;
    string first_field_s;
    
    FE2_Rhs_rr_PSFEM(FieldInterface& _mField,
                     Mat &_Aij,
                     Vec _X,
                     Vec _F,
                     ublas::matrix<FieldData> _Dmat_r,
                     ublas::matrix<FieldData> _Dmat_s,
                     string _zeroth_field,
                     ublas::matrix<FieldData> _Dmat_rs,
                     string _first_field_r,
                     string _first_field_s):
    FE2_ElasticFEMethod(_mField,
                      _Aij,
                      _X,
                      _F,
                      _Dmat_r,
                      _zeroth_field),
                      ddF(_F),
                      Dmat_r(_Dmat_r),
                      Dmat_s(_Dmat_s),
                      zeroth_field(_zeroth_field),
                      Dmat_rs(_Dmat_rs),
                      first_field_r(_first_field_r),
                      first_field_s(_first_field_s){};
    
    
    ublas::matrix<ublas::matrix<FieldData> > K_r;
    virtual PetscErrorCode StiffnessK_r() {
      PetscFunctionBegin;
      K_r.resize(row_mat,col_mat);
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
    
    
    ublas::matrix<ublas::matrix<FieldData> > K_s;
    virtual PetscErrorCode StiffnessK_s() {
      PetscFunctionBegin;
      K_s.resize(row_mat,col_mat);
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
                      w,&*Dmat_s.data().begin(),Dmat_s.size2(),
                      &*row_Mat.data().begin(),row_Mat.size2(),
                      0.,&*BD.data().begin(),BD.size2());
          for(int cc = rr;cc<col_mat;cc++) {
            if(ColGlob[cc].size()==0) continue;
            ublas::matrix<FieldData> &col_Mat = (colBMatrices[cc])[gg];
            if(gg == 0) {
              K_s(rr,cc).resize(BD.size2(),col_Mat.size2());
              cblas_dgemm(CblasRowMajor,CblasTrans,CblasNoTrans,
                          BD.size2(),col_Mat.size2(),BD.size1(),
                          1.,&*BD.data().begin(),BD.size2(),
                          &*col_Mat.data().begin(),col_Mat.size2(),
                          0.,&*K_s(rr,cc).data().begin(),K_s(rr,cc).size2());
            } else {
              cblas_dgemm(CblasRowMajor,CblasTrans,CblasNoTrans,
                          BD.size2(),col_Mat.size2(),BD.size1(),
                          1.,&*BD.data().begin(),BD.size2(),
                          &*col_Mat.data().begin(),col_Mat.size2(),
                          1.,&*K_s(rr,cc).data().begin(),K_s(rr,cc).size2());
            }
          }
        }
      }
      PetscFunctionReturn(0);
    }
    
    
    
    ublas::matrix<ublas::matrix<FieldData> > K_rs;
    virtual PetscErrorCode StiffnessK_rs() {
      PetscFunctionBegin;
      K_rs.resize(row_mat,col_mat);
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
                      w,&*Dmat_rs.data().begin(),Dmat_rs.size2(),
                      &*row_Mat.data().begin(),row_Mat.size2(),
                      0.,&*BD.data().begin(),BD.size2());
          for(int cc = rr;cc<col_mat;cc++) {
            if(ColGlob[cc].size()==0) continue;
            ublas::matrix<FieldData> &col_Mat = (colBMatrices[cc])[gg];
            if(gg == 0) {
              K_rs(rr,cc).resize(BD.size2(),col_Mat.size2());
              cblas_dgemm(CblasRowMajor,CblasTrans,CblasNoTrans,
                          BD.size2(),col_Mat.size2(),BD.size1(),
                          1.,&*BD.data().begin(),BD.size2(),
                          &*col_Mat.data().begin(),col_Mat.size2(),
                          0.,&*K_rs(rr,cc).data().begin(),K_rs(rr,cc).size2());
            } else {
              cblas_dgemm(CblasRowMajor,CblasTrans,CblasNoTrans,
                          BD.size2(),col_Mat.size2(),BD.size1(),
                          1.,&*BD.data().begin(),BD.size2(),
                          &*col_Mat.data().begin(),col_Mat.size2(),
                          1.,&*K_rs(rr,cc).data().begin(),K_rs(rr,cc).size2());
            }
          }
        }
      }
      PetscFunctionReturn(0);
    }
    
    
    vector<ublas::vector<FieldData> > f_el_rs; // element force
    virtual PetscErrorCode Rhs() {
      PetscFunctionBegin;
      //cout<<" Rhs() "<<endl;
      ierr = StiffnessK_r(); CHKERRQ(ierr);    // get K_r
      ierr = StiffnessK_s(); CHKERRQ(ierr);
      ierr = StiffnessK_rs(); CHKERRQ(ierr);   // get K_rs
      
      // displacements for nodes in each element and,
      // first-order derivative of displacements for nodes in each element
      vector<ublas::vector<FieldData> > D_elm;
      vector<ublas::vector<FieldData> > D_elm_r;
      vector<ublas::vector<FieldData> > D_elm_s;
      //     cout<<"col_mat = "<< col_mat << endl;
      D_elm.resize(col_mat);
      D_elm_r.resize(col_mat);
      D_elm_s.resize(col_mat);
      
      int col_mat1 = 0;  //only nodes (1st order)
      ierr = GetDataVector(zeroth_field,D_elm[col_mat1]); CHKERRQ(ierr);
      ierr = GetDataVector(first_field_r,D_elm_r[col_mat1]); CHKERRQ(ierr);
      ierr = GetDataVector(first_field_s,D_elm_s[col_mat1]); CHKERRQ(ierr);
      //    cout<<"D_elm[col_mat] = "<< D_elm[col_mat1] << endl;
      col_mat1++;
      
      for(int ee=0; ee<6; ee++) { //edges
        if(ColGlob[col_mat1].size()!=0) {
          ierr = GetDataVector(zeroth_field,MBEDGE,D_elm[col_mat1],ee); CHKERRQ(ierr);
          ierr = GetDataVector(first_field_r,MBEDGE,D_elm_r[col_mat1],ee); CHKERRQ(ierr);
          ierr = GetDataVector(first_field_s,MBEDGE,D_elm_s[col_mat1],ee); CHKERRQ(ierr);
          //          cout<<"Edges D_elm[col_mat1] = "<< D_elm[col_mat1] << endl;
          col_mat1++;
        }
      }
      
      for(int ff=0; ff<4; ff++) { //faces
        if(ColGlob[col_mat1].size()!=0) {
          ierr = GetDataVector(zeroth_field,MBTRI,D_elm[col_mat1],ff); CHKERRQ(ierr);
          ierr = GetDataVector(first_field_r,MBTRI,D_elm_r[col_mat1],ff); CHKERRQ(ierr);
          ierr = GetDataVector(first_field_s,MBTRI,D_elm_s[col_mat1],ff); CHKERRQ(ierr);
          //          cout<<"Faces D_elm[col_mat1] = "<< D_elm[col_mat1] << endl;
          col_mat1++;
        }
      }
      
      if(ColGlob[col_mat1].size()!=0) { // volumes
        ierr = GetDataVector(zeroth_field,MBTET,D_elm[col_mat1]); CHKERRQ(ierr);
        ierr = GetDataVector(first_field_r,MBTET,D_elm_r[col_mat1]); CHKERRQ(ierr);
        ierr = GetDataVector(first_field_s,MBTET,D_elm_s[col_mat1]); CHKERRQ(ierr);
        //    cout<<"Faces D_elm[col_mat] = "<< D_elm[col_mat1] << endl;
      }
      
      // calculate element nodal forces, f_el_rs
      f_el_rs.resize(row_mat);
      for(int rr = 0;rr<row_mat;rr++) {
        if(RowGlob[rr].size()==0) continue;
        int rr_start=0;
        for(int cc = 0;cc<col_mat;cc++) {
          if(ColGlob[cc].size()==0) continue;
          //          cout<<"rr "<<rr<<endl;
          //          cout<<"cc "<<cc<<endl;
          if(rr_start == 0) {
            //            cout<<"K_r(rr,cc) "<<K_r(rr,cc)<<endl;
            //            cout<<"D_elm_r[cc] "<<D_elm_r[cc]<<endl;
            f_el_rs[rr] = - prod( K_rs(rr,cc), D_elm[cc] )
            - prod(K_r(rr,cc), D_elm_s[cc]) - prod(K_s(rr,cc), D_elm_r[cc]);
            rr_start++;
          } else {
            f_el_rs[rr] -= prod( K_rs(rr,cc), D_elm[cc] )
            + prod(K_r(rr,cc), D_elm_s[cc]) + prod(K_s(rr,cc), D_elm_r[cc]);
          }
        }
        //      cout<<"f_el_rs[rr] "<<f_el_rs[rr]<<endl;
      }
      
      // assemble the obtained element nodal forces, fe_rs, into the
      // global nodal force vector, ddF,
      for(int rr = 0;rr<row_mat;rr++) {
        if(RowGlob[rr].size()==0) continue;
        if(RowGlob[rr].size()!=f_el_rs[rr].size()) SETERRQ(PETSC_COMM_SELF,1,"data inconsistency");
        ierr = VecSetValues(ddF,RowGlob[rr].size(),&(RowGlob[rr])[0],&(f_el_rs[rr].data()[0]),ADD_VALUES); CHKERRQ(ierr);
      }
      PetscFunctionReturn(0);
    }


    virtual PetscErrorCode Fint() {
      PetscFunctionBegin;
//      cout<<"Fint "<<endl;
      try {cout<<fieldName<<endl;
        
        //Higher order approximation of geometry
        ierr = GetHierarchicalGeometryApproximation(invH,detH); CHKERRQ(ierr);
        
        //Gradient at Gauss points;
        vector< ublas::matrix< FieldData > > GradU_at_GaussPt;
        ierr = GetGaussDiffDataVector(fieldName,GradU_at_GaussPt); CHKERRQ(ierr);//fieldname
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
            ublas::vector<FieldData> VoightStress = prod(w*Dmat_rs,VoightStrain);
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

    PetscErrorCode postProcess() {
      PetscFunctionBegin;
      PetscFunctionReturn(0);
    }

    PetscErrorCode operator()() {
      PetscFunctionBegin;
      
      ierr = Get_g_NTET(); CHKERRQ(ierr);
      ierr = OpStudentStart_TET(g_NTET); CHKERRQ(ierr);
//      cout<<"Hi from K_rPoissonP_ElasticFEMethodTransIso "<<endl;
      ierr = GetMatrices(); CHKERRQ(ierr);
      ierr = Rhs(); CHKERRQ(ierr);
//      std::string wait;
//      std::cin >> wait;
      ierr = OpStudentEnd(); CHKERRQ(ierr);
      PetscFunctionReturn(0);
    }

  };
  
  
  struct FE2_Rhs_rs_PSFEM: public FE2_ElasticFEMethod {
    
    ublas::matrix<FieldData> Dmat_r; // 1st-order derivative of D matrix w.r.t. the r-th variable
    ublas::matrix<FieldData> Dmat_rs; // 2nd-order derivative of D matrix w.r.t. the r- and s-th variables
    Vec ddF;
    string zeroth_field;
    string first_field_r;
    
    FE2_Rhs_rs_PSFEM(FieldInterface& _mField,
                     Mat &_Aij,
                     Vec _X,
                     Vec _F,
                     ublas::matrix<FieldData> _Dmat_r,
                     string _zeroth_field,
                     ublas::matrix<FieldData> _Dmat_rs,
                     string _first_field_r):
    FE2_ElasticFEMethod(_mField,
                        _Aij,
                        _X,
                        _F,
                        _Dmat_r,
                        _zeroth_field),
    ddF(_F),
    Dmat_r(_Dmat_r),
    zeroth_field(_zeroth_field),
    Dmat_rs(_Dmat_rs),
    first_field_r(_first_field_r){};
    
    
    ublas::matrix<ublas::matrix<FieldData> > K_r;
    virtual PetscErrorCode StiffnessK_r() {
      PetscFunctionBegin;
      K_r.resize(row_mat,col_mat);
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
    
    ublas::matrix<ublas::matrix<FieldData> > K_rs;
    virtual PetscErrorCode StiffnessK_rs() {
      PetscFunctionBegin;
      K_rs.resize(row_mat,col_mat);
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
                      w,&*Dmat_rs.data().begin(),Dmat_rs.size2(),
                      &*row_Mat.data().begin(),row_Mat.size2(),
                      0.,&*BD.data().begin(),BD.size2());
          for(int cc = rr;cc<col_mat;cc++) {
            if(ColGlob[cc].size()==0) continue;
            ublas::matrix<FieldData> &col_Mat = (colBMatrices[cc])[gg];
            if(gg == 0) {
              K_rs(rr,cc).resize(BD.size2(),col_Mat.size2());
              cblas_dgemm(CblasRowMajor,CblasTrans,CblasNoTrans,
                          BD.size2(),col_Mat.size2(),BD.size1(),
                          1.,&*BD.data().begin(),BD.size2(),
                          &*col_Mat.data().begin(),col_Mat.size2(),
                          0.,&*K_rs(rr,cc).data().begin(),K_rs(rr,cc).size2());
            } else {
              cblas_dgemm(CblasRowMajor,CblasTrans,CblasNoTrans,
                          BD.size2(),col_Mat.size2(),BD.size1(),
                          1.,&*BD.data().begin(),BD.size2(),
                          &*col_Mat.data().begin(),col_Mat.size2(),
                          1.,&*K_rs(rr,cc).data().begin(),K_rs(rr,cc).size2());
            }
          }
        }
      }
      PetscFunctionReturn(0);
    }
    
    
    vector<ublas::vector<FieldData> > f_el_rs; // element force
    virtual PetscErrorCode Rhs() {
      PetscFunctionBegin;
      //cout<<" Rhs() "<<endl;
      ierr = StiffnessK_r(); CHKERRQ(ierr);    // get K_r
      ierr = StiffnessK_rs(); CHKERRQ(ierr);   // get K_rs
      
      // displacements for nodes in each element and,
      // first-order derivative of displacements for nodes in each element
      vector<ublas::vector<FieldData> > D_elm;
      vector<ublas::vector<FieldData> > D_elm_r;

      //     cout<<"col_mat = "<< col_mat << endl;
      D_elm.resize(col_mat);
      D_elm_r.resize(col_mat);
      
      int col_mat1 = 0;  //only nodes (1st order)
      ierr = GetDataVector(zeroth_field,D_elm[col_mat1]); CHKERRQ(ierr);
      ierr = GetDataVector(first_field_r,D_elm_r[col_mat1]); CHKERRQ(ierr);
      //    cout<<"D_elm[col_mat] = "<< D_elm[col_mat1] << endl;
      col_mat1++;
      
      for(int ee=0; ee<6; ee++) { //edges
        if(ColGlob[col_mat1].size()!=0) {
          ierr = GetDataVector(zeroth_field,MBEDGE,D_elm[col_mat1],ee); CHKERRQ(ierr);
          ierr = GetDataVector(first_field_r,MBEDGE,D_elm_r[col_mat1],ee); CHKERRQ(ierr);
          //          cout<<"Edges D_elm[col_mat1] = "<< D_elm[col_mat1] << endl;
          col_mat1++;
        }
      }
      
      for(int ff=0; ff<4; ff++) { //faces
        if(ColGlob[col_mat1].size()!=0) {
          ierr = GetDataVector(zeroth_field,MBTRI,D_elm[col_mat1],ff); CHKERRQ(ierr);
          ierr = GetDataVector(first_field_r,MBTRI,D_elm_r[col_mat1],ff); CHKERRQ(ierr);
          //          cout<<"Faces D_elm[col_mat1] = "<< D_elm[col_mat1] << endl;
          col_mat1++;
        }
      }
      
      if(ColGlob[col_mat1].size()!=0) { // volumes
        ierr = GetDataVector(zeroth_field,MBTET,D_elm[col_mat1]); CHKERRQ(ierr);
        ierr = GetDataVector(first_field_r,MBTET,D_elm_r[col_mat1]); CHKERRQ(ierr);
        //    cout<<"Faces D_elm[col_mat] = "<< D_elm[col_mat1] << endl;
      }
      
      // calculate element nodal forces, f_el_rs
      f_el_rs.resize(row_mat);
      for(int rr = 0;rr<row_mat;rr++) {
        if(RowGlob[rr].size()==0) continue;
        int rr_start=0;
        for(int cc = 0;cc<col_mat;cc++) {
          if(ColGlob[cc].size()==0) continue;
          //          cout<<"rr "<<rr<<endl;
          //          cout<<"cc "<<cc<<endl;
          if(rr_start == 0) {
            //            cout<<"K_r(rr,cc) "<<K_r(rr,cc)<<endl;
            //            cout<<"D_elm_r[cc] "<<D_elm_r[cc]<<endl;
            f_el_rs[rr] = - prod( K_rs(rr,cc), D_elm[cc] )
            - 2*prod(K_r(rr,cc), D_elm_r[cc]);
            rr_start++;
          } else {
            f_el_rs[rr] -= prod( K_rs(rr,cc), D_elm[cc] )
            + 2*prod(K_r(rr,cc), D_elm_r[cc]);
          }
        }
        //      cout<<"f_el_rs[rr] "<<f_el_rs[rr]<<endl;
      }
      
      // assemble the obtained element nodal forces, fe_rs, into the
      // global nodal force vector, ddF,
      for(int rr = 0;rr<row_mat;rr++) {
        if(RowGlob[rr].size()==0) continue;
        if(RowGlob[rr].size()!=f_el_rs[rr].size()) SETERRQ(PETSC_COMM_SELF,1,"data inconsistency");
        ierr = VecSetValues(ddF,RowGlob[rr].size(),&(RowGlob[rr])[0],&(f_el_rs[rr].data()[0]),ADD_VALUES); CHKERRQ(ierr);
      }
      PetscFunctionReturn(0);
    }
    
    
    virtual PetscErrorCode Fint() {
      PetscFunctionBegin;
      //      cout<<"Fint "<<endl;
      try {cout<<fieldName<<endl;
        
        //Higher order approximation of geometry
        ierr = GetHierarchicalGeometryApproximation(invH,detH); CHKERRQ(ierr);
        
        //Gradient at Gauss points;
        vector< ublas::matrix< FieldData > > GradU_at_GaussPt;
        ierr = GetGaussDiffDataVector(fieldName,GradU_at_GaussPt); CHKERRQ(ierr);//fieldname
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
            ublas::vector<FieldData> VoightStress = prod(w*Dmat_rs,VoightStrain);
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
    
    PetscErrorCode postProcess() {
      PetscFunctionBegin;
      PetscFunctionReturn(0);
    }
    
    PetscErrorCode operator()() {
      PetscFunctionBegin;
      
      ierr = Get_g_NTET(); CHKERRQ(ierr);
      ierr = OpStudentStart_TET(g_NTET); CHKERRQ(ierr);
      //      cout<<"Hi from K_rPoissonP_ElasticFEMethodTransIso "<<endl;
      ierr = GetMatrices(); CHKERRQ(ierr);
      ierr = Rhs(); CHKERRQ(ierr);
      //      std::string wait;
      //      std::cin >> wait;
      ierr = OpStudentEnd(); CHKERRQ(ierr);
      PetscFunctionReturn(0);
    }
    
  };
  
}

#endif //__FE2_RHS_RS_PSFEM_HPP__
