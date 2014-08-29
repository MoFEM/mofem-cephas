/* Copyright (C) 2014, Xiao-Yi Zhou (xiaoyi.zhou AT newcastle.ac.uk)
 * --------------------------------------------------------------
 * This routine serves as a general template to construct the 1st-order partial
 * derivative of material constitutive matrix, D_r, the corresponding 1st-order
 * derivative of stiffness matrix, K_r, and the right-hand side of the 1st-order
 * finite-element equilibrium equation, Rhs = [K_r][U] for implementing 
 * perturbation-based stochastic finite element (PSFE).
 *
 * HISTORY
 *
 * 2014.08.29 (first version)
 *
 * REFERENCES
 * 1. Kleiber M. and Hien T. D. (1992) The stochastic finite element method - 
 *      Basic perturbation technique and computer implementation. John Wiley & 
 *      Sons.
 *
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

#ifndef __K_R_PSFEMETHOD_HPP__
#define __K_R_PSFEMETHOD_HPP__

#include <boost/numeric/ublas/symmetric.hpp>
extern "C" {
#include <gm_rule.h>
}

namespace MoFEM {
  
  struct K_r_PSFEMethod: public ElasticFEMethod {
    
    double young;
    double pois;
    Vec dF;
    K_r_PSFEMethod( FieldInterface& _mField,Mat &_Aij,Vec _X,Vec _F,double _young,double _pois):
    ElasticFEMethod( _mField,_Aij,_X,_F,0,0,"DISPLACEMENT"),young(_young),pois(_pois),dF(_F) {
      
    }

    //F: calculate the first-order derivative of constitutive matrix
    virtual PetscErrorCode calculateD_r(double young, double nu) {
      PetscFunctionBegin;
      D.resize(6,6);
      D.clear(); // Initiate D_r to zero
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
    
    //F: calculate the first-order partial derivative of stiffness matrix
    //   with respect to Poisson's ratio
    //   K_r = ϝ [Β]'[D_r][B]dV
    ublas::matrix<ublas::matrix<FieldData> > K_r;
    ublas::matrix<FieldData> BD; // get value for [B]' * [D_r]
    virtual PetscErrorCode StiffnessK_r() {
      PetscFunctionBegin;
      
      double _young,_pois;
      ierr = GetMatParameters(&_young,&_pois); CHKERRQ(ierr);
      ierr = calculateD_r(_young,_pois); CHKERRQ(ierr);
      //      cout<<" D "<<D<<endl;
      //      cout<<" row_mat "<<row_mat<<endl;
      //      cout<<" col_mat "<<col_mat<<endl;
      K_r.resize(row_mat,col_mat);
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
          for(int cc = 0;cc<col_mat;cc++) {
            if(ColGlob[cc].size()==0) continue;
            ublas::matrix<FieldData> &col_Mat = (colBMatrices[cc])[gg];
            if(gg == 0) {
              K_r(rr,cc).resize(BD.size2(),col_Mat.size2());
              //ublas::noalias(K(rr,cc)) = prod(trans(BD) , col_Mat ); // int BT*D*B
              cblas_dgemm(CblasRowMajor,CblasTrans,CblasNoTrans,
                          BD.size2(),col_Mat.size2(),BD.size1(),
                          1.,&*BD.data().begin(),BD.size2(),
                          &*col_Mat.data().begin(),col_Mat.size2(),
                          0.,&*K_r(rr,cc).data().begin(),K_r(rr,cc).size2());
            } else {
              //ublas::noalias(K(rr,cc)) += prod(trans(BTD) , col_Mat ); // int BT*D*B
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
    
    //F: calculate the right-hand side of the first-order finite element 
    //   equilibrium [K][U_r] = [F_r]
    //   [F_r] = - [K_r][U]
    //   Note: the applied force is assumed to be deterministic.  
    vector<ublas::vector<FieldData> > f_el_r;
    virtual PetscErrorCode Rhs() {
      PetscFunctionBegin;
      //      cout<<" Rhs() "<<endl;
      ierr = StiffnessK_r(); CHKERRQ(ierr);
      
      vector<ublas::vector<FieldData> > D_elm;
      //       cout<<"col_mat = "<< col_mat << endl;
      
      D_elm.resize(col_mat);
      
      int col_mat1 = 0;  //nodes
      ierr = GetDataVector("DISPLACEMENT",D_elm[col_mat1]); CHKERRQ(ierr);
      //       cout<<"D_elm[col_mat1] = "<< D_elm[col_mat1] << endl;
      col_mat1++;
      
      for(int ee=0; ee<6; ee++) { //edges
        if(ColGlob[col_mat1].size()!=0) {
          ierr = GetDataVector("DISPLACEMENT",MBEDGE,D_elm[col_mat1],ee); CHKERRQ(ierr);
          //          cout<<"Edges D_elm[col_mat1] = "<< D_elm[col_mat1] << endl;
          col_mat1++;
        }
      }
      
      for(int ff=0; ff<4; ff++) { //faces
        if(ColGlob[col_mat1].size()!=0) {
          ierr = GetDataVector("DISPLACEMENT",MBTRI,D_elm[col_mat1],ff); CHKERRQ(ierr);
          //          cout<<"Faces D_elm[col_mat1] = "<< D_elm[col_mat1] << endl;
          col_mat1++;
        }
      }
      
      if(ColGlob[col_mat1].size()!=0) {
        ierr = GetDataVector("DISPLACEMENT",MBTET,D_elm[col_mat1]); CHKERRQ(ierr);
        //        cout<<"Faces D_elm[col_mat] = "<< D_elm[col_mat1] << endl;
      }
      
      f_el_r.resize(row_mat);
      for(int rr = 0;rr<row_mat;rr++) {
        if(RowGlob[rr].size()==0) continue;
        
        int rr_start=0;
        for(int cc = 0;cc<col_mat;cc++) {
          if(ColGlob[cc].size()==0) continue;
//          cout<<"rr "<<rr<<endl;
//          cout<<"cc "<<cc<<endl;
          if(rr_start == 0) {
//            cout<<"K(rr,cc) "<<K(rr,cc)<<endl;
//            cout<<"D_elm[cc] "<<D_elm[cc]<<endl;
            f_el_r[rr] =  -prod(K_r(rr,cc),D_elm[cc]);
            rr_start++;
          } else {
            f_el_r[rr] -= prod(K_r(rr,cc),D_elm[cc]);
          }
         }
//        cout<<"f_el_r[rr] "<<f_el_r[rr]<<endl;
      }

      
      for(int rr = 0;rr<row_mat;rr++) {
        if(RowGlob[rr].size()==0) continue;
        if(RowGlob[rr].size()!=f_el_r[rr].size()) SETERRQ(PETSC_COMM_SELF,1,"data inconsistency");
        ierr = VecSetValues(dF,RowGlob[rr].size(),&(RowGlob[rr])[0],&(f_el_r[rr].data()[0]),ADD_VALUES); CHKERRQ(ierr);
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

#endif //__K_R_PSFEMETHOD_HPP__
