/* Copyright (C) 2014, 
 Zahur Ullah (Zahur.Ullah@glasgow.ac.uk)
 Xiao-Yi Zhou (xiaoyi.zhou@ncl.ac.uk)
 * --------------------------------------------------------------
 * This routine calculates the second-order derivative of stiffness matrix,
 * K_rs, with respect to Poisson's ratio for isotropic material for the 
 * stochastic finite element
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

#ifndef __ISO_RHS_RS_PSFEM_HPP__
#define __ISO_RHS_RS_PSFEM_HPP__

#include <boost/numeric/ublas/symmetric.hpp>
extern "C" {
#include <gm_rule.h>
}

namespace MoFEM {

struct K_rsPoissonFEMethod: public ElasticFEMethod {

    double young;
    double pois;
    Vec ddF; // 
    const string second_field;

    K_rsPoissonFEMethod( FieldInterface& _mField,Mat &_Aij,Vec _X,Vec _F,double _young,double _pois,const string& _second_field):
      ElasticFEMethod( _mField,_Aij,_X,_F,0,0,"DISPLACEMENT"),young(_young),pois(_pois),ddF(_F),second_field(_second_field){
        
    }
    //F: calculate the first-order derivative of constitutive matrix  
    //   with respect to Poisson's ratio
    virtual PetscErrorCode calculateD_r(double young, double nu) {
      PetscFunctionBegin;
      D.resize(6,6);
      D.clear();
      double D_r00,D_r01,D_r33;
      D_r00 = -(2*young*nu*(nu - 2))/pow((2*nu*nu + nu - 1),2);
      D_r01 = (young*(2*nu*nu + 1))/pow((2*nu*nu + nu - 1),2);
      D_r33 = -young/(2*pow((nu + 1),2));
      
      D(0,0)=D_r00;  D(0,1)=D_r01;   D(0,2)=D_r01;
      D(1,0)=D_r01;  D(1,1)=D_r00;   D(1,2)=D_r01;
      D(2,0)=D_r01;  D(2,1)=D_r01;   D(2,2)=D_r00;
      D(3,3)=D_r33;
      D(4,4)=D_r33;
      D(5,5)=D_r33;
      // cout<<"D_r = "<<D;
      PetscFunctionReturn(0);
    }

    //F: calculate the second-order derivative of constitutive matrix 
    //   with respect to Poisson's ratio
    virtual PetscErrorCode calculateD_rs(double young, double nu) {
      PetscFunctionBegin;
      D.resize(6,6);
      D.clear(); // Initiate D_rs to zero

      double D_rs00, D_rs01, D_rs33;

      D_rs00 = -(4*young*(- 2*nu*nu*nu + 6*nu*nu + 1))/pow((2*nu*nu + nu - 1),3);
      D_rs01 = -(2*young*(4*nu*nu*nu + 6*nu + 1))/pow((2*nu*nu + nu - 1),3);
      D_rs33 = young/pow((nu + 1),3);
      
      
       // Set nonzero values
      D(0,0) = D_rs00; D(0,1) = D_rs01;  D(0,2) = D_rs01; 
      D(1,0) = D_rs01; D(1,1) = D_rs00;  D(1,2) = D_rs01; 
      D(2,0) = D_rs01; D(2,1) = D_rs01;  D(2,2) = D_rs00;
      D(3,3) = D_rs33;
      D(4,4) = D_rs33;
      D(5,5) = D_rs33;

      // cout<<"D_rs = "<<D;

      PetscFunctionReturn(0);
    }
    
    //F: Retrieve material properties
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

    //F: calculate the second-order partial derivative of stiffness matrix
    //   with respect to Poisson's ratio
    //   K_rs = ϝ [Β]'[D_rs][B]dV
    ublas::matrix<ublas::matrix<FieldData> > K_rs;
    // ublas::matrix<FieldData> BD; // get value for [B]' * [D_r]
    virtual PetscErrorCode StiffnessK_rs() {
      PetscFunctionBegin;
      
      double _young,_pois;
      ierr = GetMatParameters(&_young,&_pois); CHKERRQ(ierr);
      ierr = calculateD_rs(_young,_pois); CHKERRQ(ierr);
//      cout<<" D_rs "<<D<<endl;
      //      cout<<" row_mat "<<row_mat<<endl;
      //      cout<<" col_mat "<<col_mat<<endl;
      K_rs.resize(row_mat,col_mat);
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
              K_rs(rr,cc).resize(BD.size2(),col_Mat.size2());
              //ublas::noalias(K(rr,cc)) = prod(trans(BD) , col_Mat ); // int BT*D*B
              cblas_dgemm(CblasRowMajor,CblasTrans,CblasNoTrans,
                          BD.size2(),col_Mat.size2(),BD.size1(),
                          1.,&*BD.data().begin(),BD.size2(),
                          &*col_Mat.data().begin(),col_Mat.size2(),
                          0.,&*K_rs(rr,cc).data().begin(),K_rs(rr,cc).size2());
            } else {
              //ublas::noalias(K(rr,cc)) += prod(trans(BTD) , col_Mat ); // int BT*D*B
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
    //F: calculate the second-order partial derivative of "external force" which
    //   is referred as right-hand side in the algebriac equation [K][U_rs] = [F_rs]
    //   [F_rs] = - [K_rs][U] - 2[K_r][U_s]  
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
    //       cout<<"col_mat = "<< col_mat << endl;
    D_elm.resize(col_mat);
    D_elm_r.resize(col_mat);
    
    int col_mat1 = 0;  //only nodes (1st order)
    ierr = GetDataVector("DISPLACEMENT",D_elm[col_mat1]); CHKERRQ(ierr);
    ierr = GetDataVector(second_field,D_elm_r[col_mat1]); CHKERRQ(ierr);
    col_mat1++;
    //       cout<<"D_elm[col_mat] = "<< D_elm[col_mat1] << endl;
    
    for(int ee=0; ee<6; ee++) { //edges
      if(ColGlob[col_mat1].size()!=0) {
        ierr = GetDataVector("DISPLACEMENT",MBEDGE,D_elm[col_mat1],ee); CHKERRQ(ierr);
        ierr = GetDataVector(second_field,MBEDGE,D_elm_r[col_mat1],ee); CHKERRQ(ierr);
        //          cout<<"Edges D_elm[col_mat1] = "<< D_elm[col_mat1] << endl;
        col_mat1++;
      }
    }
    
    for(int ff=0; ff<4; ff++) { //faces
      if(ColGlob[col_mat1].size()!=0) {
        ierr = GetDataVector("DISPLACEMENT",MBTRI,D_elm[col_mat1],ff); CHKERRQ(ierr);
        ierr = GetDataVector(second_field,MBTRI,D_elm_r[col_mat1],ff); CHKERRQ(ierr);
        //          cout<<"Faces D_elm[col_mat1] = "<< D_elm[col_mat1] << endl;
        col_mat1++;
      }
    }
    
    if(ColGlob[col_mat1].size()!=0) { // volumes
      ierr = GetDataVector("DISPLACEMENT",MBTET,D_elm[col_mat1]); CHKERRQ(ierr);
      ierr = GetDataVector(second_field,MBTET,D_elm_r[col_mat1]); CHKERRQ(ierr);
      //cout<<"Faces D_elm[col_mat] = "<< D_elm[col_mat1] << endl;
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
  

    PetscErrorCode postProcess() {
      PetscFunctionBegin;
      PetscFunctionReturn(0);
    }

  
     PetscErrorCode operator()() {
      PetscFunctionBegin;
      ierr = OpStudentStart_TET(g_NTET); CHKERRQ(ierr);
      //cout<<"Hi from K_rPoissonFEMethod "<<endl;
      //      cout<<"fieldName "<<fieldName<<endl;
      ierr = GetMatrices(); CHKERRQ(ierr);
      ierr = Rhs(); CHKERRQ(ierr);
      ierr = OpStudentEnd(); CHKERRQ(ierr);
      PetscFunctionReturn(0);
    }
};

  
  
  
struct K_rsYoungFEMethod: public K_rsPoissonFEMethod {
    K_rsYoungFEMethod( FieldInterface& _mField,Mat &_Aij,Vec _X,Vec _F,double _young,double _pois,const string& _second_field):
    K_rsPoissonFEMethod( _mField,_Aij,_X,_F,_young,_pois,_second_field){
      
    }
    //   with respect to young modulus
    virtual PetscErrorCode calculateD_r(double young, double nu) {
      PetscFunctionBegin;
      D.resize(6,6);
      D.clear();

      double D_r00,D_r01,D_r33,constt;
      constt=1/((1+nu)*(1-2*nu));
      
      D_r00 =constt*(1-nu);
      D_r01 = constt*nu;
      D_r33 = constt*(1-2*nu)/2;
      
      D(0,0)=D_r00;  D(0,1)=D_r01;   D(0,2)=D_r01;
      D(1,0)=D_r01;  D(1,1)=D_r00;   D(1,2)=D_r01;
      D(2,0)=D_r01;  D(2,1)=D_r01;   D(2,2)=D_r00;
      D(3,3)=D_r33;
      D(4,4)=D_r33;
      D(5,5)=D_r33;
      // cout<<"D_r = "<<D;
      PetscFunctionReturn(0);
    }
    
    //F: calculate the second-order derivative of constitutive matrix
    //   with respect to young modulus 2 times
    virtual PetscErrorCode calculateD_rs(double young, double nu) {
      PetscFunctionBegin;
      D.resize(6,6);
      D.clear(); // Initiate D_rs to zero
      PetscFunctionReturn(0);
    }
    
  };

  
  
  
  struct K_rYoungPoissonFEMethod: public K_rsPoissonFEMethod {
    K_rYoungPoissonFEMethod( FieldInterface& _mField,Mat &_Aij,Vec _X,Vec _F,double _young,double _pois,const string& _second_field):
    K_rsPoissonFEMethod( _mField,_Aij,_X,_F,_young,_pois,_second_field){
      
    }
    //F: calculate the first-order derivative of constitutive matrix
    //   with respect to Poisson's ratio
    virtual PetscErrorCode calculateD_r(double young, double nu) {
      PetscFunctionBegin;
      D.resize(6,6);
      D.clear();
      double D_r00,D_r01,D_r33,constt;
      constt=1/((1+nu)*(1-2*nu));
      
      D_r00 = constt*(1-nu);
      D_r01 = constt*nu;
      D_r33 = constt*(1-2*nu)/2;
      
      D(0,0)=D_r00;  D(0,1)=D_r01;   D(0,2)=D_r01;
      D(1,0)=D_r01;  D(1,1)=D_r00;   D(1,2)=D_r01;
      D(2,0)=D_r01;  D(2,1)=D_r01;   D(2,2)=D_r00;
      D(3,3)=D_r33;
      D(4,4)=D_r33;
      D(5,5)=D_r33;
      // cout<<"D_r = "<<D;
      PetscFunctionReturn(0);
    }
    
    //F: calculate the second-order derivative of constitutive matrix
    //   with respect to Poisson's ratio
    virtual PetscErrorCode calculateD_rs(double young, double nu) {
      PetscFunctionBegin;
      D.resize(6,6);
      D.clear(); // Initiate D_rs to zero
      
      double D_rs00, D_rs01, D_rs33;
      
      D_rs00 = -(2*nu*(nu - 2))/pow((2*nu*nu + nu - 1),2);
      D_rs01 = (2*nu*nu + 1)/pow((2*nu*nu + nu - 1),2);
      D_rs33 = -1/(2*pow((nu + 1),2));
      
      
      // Set nonzero values
      D(0,0) = D_rs00; D(0,1) = D_rs01;  D(0,2) = D_rs01;
      D(1,0) = D_rs01; D(1,1) = D_rs00;  D(1,2) = D_rs01;
      D(2,0) = D_rs01; D(2,1) = D_rs01;  D(2,2) = D_rs00;
      D(3,3) = D_rs33;
      D(4,4) = D_rs33;
      D(5,5) = D_rs33;
      // cout<<"D_rs = "<<D;
      PetscFunctionReturn(0);
    }
  };

  
  
  struct K_rs_EmEPf_FEMethod: public K_rsPoissonFEMethod {
    K_rs_EmEPf_FEMethod( FieldInterface& _mField,Mat &_Aij,Vec _X,Vec _F,double _young,double _pois,const string& _second_field):
    K_rsPoissonFEMethod( _mField,_Aij,_X,_F,_young,_pois,_second_field){
      
    }
    //F: calculate the first-order derivative of constitutive matrix
    //   with respect to Young's modulus
    virtual PetscErrorCode calculateD_r(double young, double nu) {
      PetscFunctionBegin;
      D.resize(6,6);
      D.clear();
      double D_r00,D_r01,D_r33;
      D_r00 = (1-nu)/((1+nu)*(1-2*nu));
      D_r01 = (nu)/((1+nu)*(1-2*nu));
      D_r33 = (1-2*nu)/(2*(1+nu)*(1-2*nu));
      D(0,0)=D_r00;  D(0,1)=D_r01;   D(0,2)=D_r01;
      D(1,0)=D_r01;  D(1,1)=D_r00;   D(1,2)=D_r01;
      D(2,0)=D_r01;  D(2,1)=D_r01;   D(2,2)=D_r00;
      D(3,3)=D_r33;
      D(4,4)=D_r33;
      D(5,5)=D_r33;
      // cout<<"D_r = "<<D;
      PetscFunctionReturn(0);
    }
    //F: calculate the second-order derivative of constitutive matrix
    //   with respect to Young's modulus of matrix and Young's modulus of fibre
    virtual PetscErrorCode calculateD_rs(double young, double nu) {
      PetscFunctionBegin;
      D.resize(6,6);
      D.clear(); // Initiate D_rs to zero
      // cout<<"D_rs = "<<D;
      PetscFunctionReturn(0);
    }
  };

  
  struct K_rs_PmEPf_FEMethod: public K_rsPoissonFEMethod {
    
    K_rs_PmEPf_FEMethod( FieldInterface& _mField,Mat &_Aij,Vec _X,Vec _F,double _young,double _pois,const string& _second_field):
    K_rsPoissonFEMethod( _mField,_Aij,_X,_F,_young,_pois,_second_field){
      
    }
    //F: calculate the first-order derivative of constitutive matrix
    //   with respect to Poisson's ratio of matrix
    virtual PetscErrorCode calculateD_r(double young, double nu) {
      PetscFunctionBegin;
      D.resize(6,6);
      D.clear();
      
      double D_r00,D_r01,D_r33;
      D_r00 = -(2*young*nu*(nu - 2))/pow((2*nu*nu + nu - 1),2);
      D_r01 = (young*(2*nu*nu + 1))/pow((2*nu*nu + nu - 1),2);
      D_r33 = -young/(2*pow((nu + 1),2));
      
      D(0,0)=D_r00;  D(0,1)=D_r01;   D(0,2)=D_r01;
      D(1,0)=D_r01;  D(1,1)=D_r00;   D(1,2)=D_r01;
      D(2,0)=D_r01;  D(2,1)=D_r01;   D(2,2)=D_r00;
      D(3,3)=D_r33;
      D(4,4)=D_r33;
      D(5,5)=D_r33;
      // cout<<"D_r = "<<D;
      PetscFunctionReturn(0);
    }
    
    //F: calculate the second-order derivative of constitutive matrix
    //   with respect to Poisson's ratio of matrix and Young's modulus of fibre
    virtual PetscErrorCode calculateD_rs(double young, double nu) {
      PetscFunctionBegin;
      D.resize(6,6);
      D.clear(); // Initiate D_rs to zero
      // cout<<"D_rs = "<<D;
      PetscFunctionReturn(0);
    }
  };
  

}

#endif //__ISO_RHS_RS_PSFEM_HPP__
