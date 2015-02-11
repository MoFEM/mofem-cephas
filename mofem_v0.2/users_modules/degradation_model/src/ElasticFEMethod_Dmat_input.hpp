/* Copyright (C) 2015, Zahur Ullah (Zahur.Ullah AT glasgow.ac.uk)
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

#ifndef __ELASTICFEMETHOD_DMAT_INPUT_HPP__
#define __ELASTICFEMETHOD_DMAT_INPUT_HPP__

namespace ObosleteUsersModules {
  
  
  
  struct ElasticFEMethod_Dmat_input: public ElasticFEMethod {
    
    map<EntityHandle, ublas::vector<ublas::matrix<double> > > &Dmat_RVE;
    ElasticFEMethod_Dmat_input(FieldInterface& _mField,Mat _Aij,Vec _X,Vec _F,double _lambda,double _mu,
      map<EntityHandle, ublas::vector<ublas::matrix<double> > > &_Dmat_RVE, string _field_name):
      ElasticFEMethod(_mField,_Aij,_X,_F,_lambda,_mu,_field_name), Dmat_RVE(_Dmat_RVE){}
    
    
    virtual PetscErrorCode calculateD(int gg) {
      PetscFunctionBegin;
      
      EntityHandle fe_ent;  fe_ent=fePtr->get_ent();
//      cout << "fe_ent =  "<< fe_ent<<endl;
//      cout << "gg =  "<< gg<<endl;
      if(Dmat_RVE.find(fe_ent)==Dmat_RVE.end()) {
        SETERRQ1(PETSC_COMM_SELF,1,"RVE material matrix not calulated for entity %ld",fe_ent);
      }
//      cerr << "Dmat_RVE[fe_ent].size() = "<<Dmat_RVE[fe_ent].size() << endl;
//      cerr << "Dmat_RVE[fe_ent][gg].size1() = "<< Dmat_RVE[fe_ent][gg].size1()<<endl;
//      cerr << "Dmat_RVE[fe_ent][gg].size2() = "<< Dmat_RVE[fe_ent][gg].size2()<<endl;
//      cout << "Dmat_RVE[fe_ent](gg) =  "<< Dmat_RVE[fe_ent](gg)<<endl;
      //accessing values from STL map Dmat_RVE (given Element and Gauss point)
      D=Dmat_RVE[fe_ent](gg);
      D=1e6*D; //conversion from MPa to Pa
//      cout << "D =  "<< D<<endl;
      PetscFunctionReturn(0);
    }
    
    
    
    
    
    virtual PetscErrorCode Stiffness() {
      PetscFunctionBegin;
      K.resize(row_mat,col_mat);
      int g_dim = g_NTET.size()/4;
//      cout<<"g_dim   =  "<<g_dim<<endl;

      for(int rr = 0;rr<row_mat;rr++) {
        if(RowGlob[rr].size()==0) continue;
        for(int gg = 0;gg<g_dim;gg++) {
          ierr = calculateD(gg); CHKERRQ(ierr);
          ublas::matrix<FieldData> &row_Mat = (rowBMatrices[rr])[gg];
          double w = V*G_TET_W[gg];
          if(detH.size()>0) {
            w *= detH[gg];
          }
          BD.resize(6,row_Mat.size2());
          //ublas::noalias(BD) = prod( w*D,row_Mat );
          cblas_dsymm(CblasRowMajor,CblasLeft,CblasUpper,
                      BD.size1(),BD.size2(),
                      w,&*D.data().begin(),D.size2(),
                      &*row_Mat.data().begin(),row_Mat.size2(),
                      0.,&*BD.data().begin(),BD.size2());
          for(int cc = rr;cc<col_mat;cc++) {
            if(ColGlob[cc].size()==0) continue;
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
    
    
    virtual PetscErrorCode Fint() {
      PetscFunctionBegin;

      try {

        //Higher order approximation of geometry
        ierr = GetHierarchicalGeometryApproximation(invH,detH); CHKERRQ(ierr);
        //Gradient at Gauss points;
        vector< ublas::matrix< FieldData > > GradU_at_GaussPt;
        ierr = GetGaussDiffDataVector(fieldName,GradU_at_GaussPt); CHKERRQ(ierr);
        unsigned int g_dim = g_NTET.size()/4;
        assert(GradU_at_GaussPt.size() == g_dim);
        NOT_USED(g_dim);
        vector< ublas::matrix< FieldData > >::iterator viit = GradU_at_GaussPt.begin();
        
//        cout<<"g_dim   =  "<<g_dim<<endl;
        
        int gg = 0;
        for(;viit!=GradU_at_GaussPt.end();viit++,gg++) {
          try {
            ierr = calculateD(gg); CHKERRQ(ierr);
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
            ublas::vector<FieldData> VoightStress = prod(w*D,VoightStrain);
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
  };
}

#endif //__ELASTICFEMETHOD_DMAT_INPUT_HPP__
