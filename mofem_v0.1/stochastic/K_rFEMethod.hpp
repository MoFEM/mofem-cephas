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

struct K_rFEMethod: public ElasticFEMethod {

    double young;
    double pois;
  
    K_rFEMethod( FieldInterface& _mField,Mat &_Aij,Vec _X,Vec _F,double _young,double _pois):
      ElasticFEMethod( _mField,_Aij,_X,_F,0,0,"DISP_r"),young(_young),pois(_pois) {
        
    }

    PetscErrorCode calculateD(double young, double nu) {
      PetscFunctionBegin;

      D.resize(6,6);
      
      double E1=1/((1+nu)*(1-2*nu));
      D(0,0)=1-nu; D(0,1)=nu;   D(0,2)=nu;    D(0,3)=0;          D(0,4)=0;          D(0,5)=0;
      D(1,0)=nu;   D(1,1)=1-nu; D(1,2)=nu;    D(1,3)=0;          D(1,4)=0;          D(1,5)=0;
      D(2,0)=nu;   D(2,1)=nu;   D(2,2)=1-nu;  D(2,3)=0;          D(2,4)=0;          D(2,5)=0;
      D(3,0)=0;    D(3,1)=0;    D(3,2)=0;     D(3,3)=(1-2*nu)/2; D(3,4)=0;          D(3,5)=0;
      D(4,0)=0;    D(4,1)=0;    D(4,2)=0;     D(4,3)=0;          D(4,4)=(1-2*nu)/2; D(4,5)=0;
      D(5,0)=0;    D(5,1)=0;    D(5,2)=0;     D(5,3)=0;          D(5,4)=0;          D(5,5)=(1-2*nu)/2;
      D=E1*D;
      //cerr << D_lambda << endl;
      //cerr << D_mu << endl;
      //cerr << D << endl;

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
//      ierr = calculateD(_young,_pois); CHKERRQ(ierr);

//      K.resize(row_mat,col_mat);
//      int g_dim = g_NTET.size()/4;
//      for(int rr = 0;rr<row_mat;rr++) {
//	if(RowGlob[rr].size()==0) continue;
//	for(int gg = 0;gg<g_dim;gg++) {
//	  ublas::matrix<FieldData> &row_Mat = (rowBMatrices[rr])[gg];
//	  double w = V*G_TET_W[gg];
//	  if(detH.size()>0) {
//	    w *= detH[gg];
//	  }
//	  BD.resize(6,row_Mat.size2());
//	  //ublas::noalias(BD) = prod( w*D,row_Mat );
//	  cblas_dsymm(CblasRowMajor,CblasLeft,CblasUpper,
//	    BD.size1(),BD.size2(),
//	    w,&*D.data().begin(),D.size2(),
//	    &*row_Mat.data().begin(),row_Mat.size2(),
//	    0.,&*BD.data().begin(),BD.size2());
//	  for(int cc = rr;cc<col_mat;cc++) {
//	    if(ColGlob[cc].size()==0) continue;
//	    ublas::matrix<FieldData> &col_Mat = (colBMatrices[cc])[gg];
//	    if(gg == 0) {
//	      K(rr,cc).resize(BD.size2(),col_Mat.size2());
//	      //ublas::noalias(K(rr,cc)) = prod(trans(BD) , col_Mat ); // int BT*D*B
//	      cblas_dgemm(CblasRowMajor,CblasTrans,CblasNoTrans,
//		BD.size2(),col_Mat.size2(),BD.size1(),
//		1.,&*BD.data().begin(),BD.size2(),
//		&*col_Mat.data().begin(),col_Mat.size2(),
//		0.,&*K(rr,cc).data().begin(),K(rr,cc).size2());
//	    } else {
//	      //ublas::noalias(K(rr,cc)) += prod(trans(BTD) , col_Mat ); // int BT*D*B
//	      cblas_dgemm(CblasRowMajor,CblasTrans,CblasNoTrans,
//		BD.size2(),col_Mat.size2(),BD.size1(),
//		1.,&*BD.data().begin(),BD.size2(),
//		&*col_Mat.data().begin(),col_Mat.size2(),
//		1.,&*K(rr,cc).data().begin(),K(rr,cc).size2());
//	    }
//	  }
//	}
//      }
      PetscFunctionReturn(0);
    }

    virtual PetscErrorCode Lhs() {
      PetscFunctionBegin;
      ierr = Stiffness(); CHKERRQ(ierr);
//      for(int rr = 0;rr<row_mat;rr++) {
//	if(RowGlob[rr].size()==0) continue;
//	for(int cc = rr;cc<col_mat;cc++) {
//	  if(ColGlob[cc].size()==0) continue;
//	  if(RowGlob[rr].size()!=K(rr,cc).size1()) SETERRQ(PETSC_COMM_SELF,1,"data inconsistency");
//	  if(ColGlob[cc].size()!=K(rr,cc).size2()) SETERRQ(PETSC_COMM_SELF,1,"data inconsistency");
//	  ierr = MatSetValues(*snes_B,RowGlob[rr].size(),&(RowGlob[rr])[0],ColGlob[cc].size(),&(ColGlob[cc])[0],&(K(rr,cc).data())[0],ADD_VALUES); CHKERRQ(ierr);
//	  if(rr!=cc) {
//	    K(cc,rr) = trans(K(rr,cc));
//	    ierr = MatSetValues(*snes_B,ColGlob[cc].size(),&(ColGlob[cc])[0],RowGlob[rr].size(),&(RowGlob[rr])[0],&(K(cc,rr).data())[0],ADD_VALUES); CHKERRQ(ierr);
//	  }
//	}
//      }
      PetscFunctionReturn(0);
    }

//    virtual PetscErrorCode Rhs() {
//      PetscFunctionBegin;
//      ierr = Fint(snes_f); CHKERRQ(ierr);
//      PetscFunctionReturn(0);
//    }

//    vector<ublas::vector<FieldData> > f_int;
//    virtual PetscErrorCode Fint() {
//      PetscFunctionBegin;
//
//      try {
//
//      //Higher order approximation of geometry
//      ierr = GetHierarchicalGeometryApproximation(invH,detH); CHKERRQ(ierr);
//
//      double _lambda,_mu;
//      ierr = GetMatParameters(&_lambda,&_mu); CHKERRQ(ierr);
//      ierr = calculateD(_lambda,_mu); CHKERRQ(ierr);
//
//      //Gradient at Gauss points; 
//      vector< ublas::matrix< FieldData > > GradU_at_GaussPt;
//      ierr = GetGaussDiffDataVector("DISP_r",GradU_at_GaussPt); CHKERRQ(ierr);
//      unsigned int g_dim = g_NTET.size()/4;
//      assert(GradU_at_GaussPt.size() == g_dim);
//      NOT_USED(g_dim);
//      vector< ublas::matrix< FieldData > >::iterator viit = GradU_at_GaussPt.begin();
//      int gg = 0;
//      for(;viit!=GradU_at_GaussPt.end();viit++,gg++) {
//	try {
//	  ublas::matrix< FieldData > GradU = *viit;
//	  if(!invH.empty()) {
//	    //GradU = 
//	      //[ dU/dChi1 dU/dChi2 dU/dChi3 ]
//	      //[ dV/dChi1 dV/dChi2 dU/dChi3 ]
//	      //[ dW/dChi1 dW/dChi2 dW/dChi3 ]
//	    //H = 
//	      //[ dX1/dChi1 dX1/dChi2 dX1/dChi3 ]
//	      //[ dX2/dChi1 dX2/dChi2 dX2/dChi3 ]
//	      //[ dX3/dChi1 dX3/dChi2 dX3/dChi3 ]    
//	    //invH = 
//	      //[ dChi1/dX1 dChi1/dX2 dChi1/dX3 ]
//	      //[ dChi2/dX1 dChi2/dX2 dChi2/dX3 ]
//	      //[ dChi3/dX1 dChi3/dX2 dChi3/dX3 ]
//	    //GradU = 
//	      //[ dU/dX1 dU/dX2 dU/dX3 ]
//	      //[ dV/dX1 dV/dX2 dV/dX3 ] = GradU * invH
//	      //[ dW/dX1 dW/dX2 dW/dX3 ] 
//	    GradU = prod( GradU, invH[gg] ); 
//	  }
//	  ublas::matrix< FieldData > Strain = 0.5*( GradU + trans(GradU) );
//	  ublas::vector< FieldData > VoightStrain(6);
//	  VoightStrain[0] = Strain(0,0);
//	  VoightStrain[1] = Strain(1,1);
//	  VoightStrain[2] = Strain(2,2);
//	  VoightStrain[3] = 2*Strain(0,1);
//	  VoightStrain[4] = 2*Strain(1,2);
//	  VoightStrain[5] = 2*Strain(2,0);
//	  double w = V*G_TET_W[gg];
//	  ublas::vector<FieldData> VoightStress = prod(w*D,VoightStrain);
//	  //BT * VoigtStress
//	  f_int.resize(row_mat);
//	  for(int rr = 0;rr<row_mat;rr++) {
//	    if(RowGlob[rr].size()==0) continue;
//	    ublas::matrix<FieldData> &B = (rowBMatrices[rr])[gg];
//	    if(gg == 0) {
//	      f_int[rr] = prod( trans(B), VoightStress );
//	    } else {
//	      f_int[rr] += prod( trans(B), VoightStress );
//	    }
//	  }
//	} catch (const std::exception& ex) {
//	  ostringstream ss;
//	  ss << "thorw in method: " << ex.what() << " at line " << __LINE__ << " in file " << __FILE__ << endl;
//	  SETERRQ(PETSC_COMM_SELF,1,ss.str().c_str());
//	} 
//      }
//
//      } catch (const std::exception& ex) {
//	ostringstream ss;
//	ss << "thorw in method: " << ex.what() << " at line " << __LINE__ << " in file " << __FILE__ << endl;
//	SETERRQ(PETSC_COMM_SELF,1,ss.str().c_str());
//      }
//
//      PetscFunctionReturn(0);
//    }

//    virtual PetscErrorCode Fint(Vec F_int) {
//      PetscFunctionBegin;
//      try {
//	ierr = Fint(); CHKERRQ(ierr);
//	for(int rr = 0;rr<row_mat;rr++) {
//	  if(RowGlob[rr].size()==0) continue;
//	  if(RowGlob[rr].size()!=f_int[rr].size()) SETERRQ(PETSC_COMM_SELF,1,"data inconsistency");
//	  ierr = VecSetValues(F_int,RowGlob[rr].size(),&(RowGlob[rr])[0],&(f_int[rr].data()[0]),ADD_VALUES); CHKERRQ(ierr);
//	}
//      } catch (const std::exception& ex) {
//	ostringstream ss;
//	ss << "thorw in method: " << ex.what() << " at line " << __LINE__ << " in file " << __FILE__ << endl;
//	SETERRQ(PETSC_COMM_SELF,1,ss.str().c_str());
//      }
//      PetscFunctionReturn(0);
//    }


//    virtual PetscErrorCode RhsAndLhs() {
//      PetscFunctionBegin;
//
//      ierr = Rhs(); CHKERRQ(ierr);
//      ierr = Lhs(); CHKERRQ(ierr);
//
//      PetscFunctionReturn(0);
//    }

//    ublas::matrix<double> gaussPts;
//    virtual PetscErrorCode Get_g_NTET() {
//      PetscFunctionBegin;
//
//      int order = 1;
//      for(_IT_GET_FEDATA_BY_NAME_DOFS_FOR_LOOP_(this,"DISP_r",dof)) {
//	order = max(order,dof->get_max_order());
//      }
//
//      int rule = max(0,order-1);
//      if( 2*rule + 1 < 2*(order-1) ) {
//	SETERRQ2(PETSC_COMM_SELF,1,"wrong rule %d %d",order,rule);
//      }
//      int nb_gauss_pts = gm_rule_size(rule,3);
//      if(gaussPts.size2() == (unsigned int)nb_gauss_pts) {
//	PetscFunctionReturn(0);
//      }
//      gaussPts.resize(4,nb_gauss_pts);
//      ierr = Grundmann_Moeller_integration_points_3D_TET(
//	rule,&gaussPts(0,0),&gaussPts(1,0),&gaussPts(2,0),&gaussPts(3,0)); CHKERRQ(ierr);
//
//      g_NTET.resize(4*nb_gauss_pts);
//      ierr = ShapeMBTET(&g_NTET[0],&gaussPts(0,0),&gaussPts(1,0),&gaussPts(2,0),nb_gauss_pts); CHKERRQ(ierr);
//      G_TET_W = &gaussPts(3,0);
//
//      PetscFunctionReturn(0);
//    }

//    PetscErrorCode preProcess() {
//      PetscFunctionBegin;
//
//      g_NTET.resize(4*45);
//      ierr = ShapeMBTET(&g_NTET[0],G_TET_X45,G_TET_Y45,G_TET_Z45,45); CHKERRQ(ierr);
//      G_TET_W = G_TET_W45;
//			
//      // See FEAP - - A Finite Element Analysis Program
//      D_lambda.resize(6,6);
//      D_lambda.clear();
//      for(int rr = 0;rr<3;rr++) {
//	for(int cc = 0;cc<3;cc++) {
//	  D_lambda(rr,cc) = 1;
//	}
//      }
//      D_mu.resize(6,6);
//      D_mu.clear();
//      for(int rr = 0;rr<6;rr++) {
//	D_mu(rr,rr) = rr<3 ? 2 : 1;
//      }
//
//      PetscFunctionReturn(0);
//    }




};

    
}

#endif //__K_RFEMETHOD_HPP__
