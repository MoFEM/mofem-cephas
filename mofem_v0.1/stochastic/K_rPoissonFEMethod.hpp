/* Copyright (C) 2014, Zahur Ullah (Zahur.Ullah@glasgow.ac.uk)
 * --------------------------------------------------------------
 * Implemnetation of the stochastic finite element componenet K_r
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

#ifndef __K_RFEMETHOD_HPP__
#define __K_RFEMETHOD_HPP__

#include <boost/numeric/ublas/symmetric.hpp>
extern "C" {
#include <gm_rule.h>
}

namespace MoFEM {
  
  struct K_rPoissonFEMethod: public ElasticFEMethod {
    
    double young;
    double pois;
    Vec dF;
    
    K_rPoissonFEMethod( FieldInterface& _mField,Mat &_Aij,Vec _X,Vec _F,double _young,double _pois):
    ElasticFEMethod( _mField,_Aij,_X,_F,0,0,"DISPLACEMENT"),young(_young),pois(_pois),dF(_F) {
      
    }
    
    PetscErrorCode calculateD(double young, double nu) {
      PetscFunctionBegin;
      D.resize(6,6);
      D.clear();
      double D00,D01,D33;
      D00=(2*young*(nu*nu + 2*nu - 2)) /pow((nu + 1),2);
      D01=-(young*(2*nu*nu + 4*nu - 1))/pow((nu + 1),2);
      D33=(young*(4*nu*nu + 8*nu - 5))/(2*pow((nu + 1),2));
      
      D(0,0)=D00;  D(0,1)=D01;  D(0,2)=D01;
      D(1,0)=D01;  D(1,1)=D00;  D(1,2)=nu;
      D(2,0)=D01;  D(2,1)=D01;   D(2,2)=D00;
      D(3,3)=D33;
      D(4,4)=D33;
      D(5,5)=D33;
      //      cout<<"D = "<<D;
      PetscFunctionReturn(0);
    }
    
    virtual PetscErrorCode GetMatParameters(double *_young,double *_pois) {
      PetscFunctionBegin;
      
      *_young = young;
      *_pois = pois;
      
      
      if(propeties_from_BLOCKSET_MAT_ELASTICSET) {
        EntityHandle ent = fePtr->get_ent();
        for(_IT_CUBITMESHSETS_BY_BCDATA_TYPE_FOR_LOOP_(mField,BLOCKSET|MAT_ELASTICSET,it)) {
          
          Mat_Elastic mydata;
          ierr = it->get_attribute_data_structure(mydata); CHKERRQ(ierr);
          
          Range meshsets;
          rval = moab.get_entities_by_type(it->meshset,MBENTITYSET,meshsets,true); CHKERR_PETSC(rval);
          meshsets.insert(it->meshset);
          for(Range::iterator mit = meshsets.begin();mit != meshsets.end(); mit++) {
            if( moab.contains_entities(*mit,&ent,1) ) {
              *_young = mydata.data.Young;
              *_pois  = mydata.data.Poisson;
              PetscFunctionReturn(0);
            }
          }
          
        }
        
        SETERRQ(PETSC_COMM_SELF,1,
                "Element is not in elestic block, however you run linear elastic analysis with that element\n"
                "top tip: check if you update block sets after mesh refinments or interface insertion");
        
      }
      
      PetscFunctionReturn(0);
    }
    
    
    ublas::matrix<ublas::matrix<FieldData> > K;
    ublas::matrix<FieldData> BD;
    virtual PetscErrorCode Stiffness() {
      PetscFunctionBegin;
      
      double _young,_pois;
      ierr = GetMatParameters(&_young,&_pois); CHKERRQ(ierr);
      ierr = calculateD(_young,_pois); CHKERRQ(ierr);
      //      cout<<" D "<<D<<endl;
      //      cout<<" row_mat "<<row_mat<<endl;
      //      cout<<" col_mat "<<col_mat<<endl;
      K.resize(row_mat,col_mat);
      int g_dim = g_NTET.size()/4;
      //      cout<<" g_dim "<<g_dim<<endl;
      for(int rr = 0;rr<row_mat;rr++) {
        if(RowGlob[rr].size()==0) continue;
        for(int gg = 0;gg<g_dim;gg++) {
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
    
    
    vector<ublas::vector<FieldData> > f_re;
    virtual PetscErrorCode Rhs() {
      PetscFunctionBegin;
      //      cout<<" Rhs() "<<endl;
      ierr = Stiffness(); CHKERRQ(ierr);
      
      vector<ublas::vector<FieldData> > D_elm;
      //       cout<<"col_mat = "<< col_mat << endl;
      
      D_elm.resize(col_mat);
      
      int col_mat1 = 0;  //only nodes (1st order)
      ierr = GetDataVector("DISPLACEMENT",D_elm[col_mat1]); CHKERRQ(ierr);
      //       cout<<"D_elm[col_mat] = "<< D_elm[col_mat1] << endl;
      
      f_re.resize(row_mat);
      for(int rr = 0;rr<row_mat;rr++) {
        if(RowGlob[rr].size()==0) continue;
        for(int cc = rr;cc<col_mat;cc++) {
          if(ColGlob[cc].size()==0) continue;
          //                  cout<<""<<prod( K(rr,cc), D_elm[cc])<<endl;
          f_re[rr] = -prod( K(rr,cc), D_elm[cc] );
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
      ierr = OpStudentStart_TET(g_NTET); CHKERRQ(ierr);
//      cout<<"Hi from K_rPoissonFEMethod "<<endl;
      //      cout<<"fieldName "<<fieldName<<endl;
      ierr = GetMatrices(); CHKERRQ(ierr);
      
      ierr = Rhs(); CHKERRQ(ierr);
      
      ierr = OpStudentEnd(); CHKERRQ(ierr);
      PetscFunctionReturn(0);
    }
    
    
  };
  
  
}

#endif //__K_RFEMETHOD_HPP__
