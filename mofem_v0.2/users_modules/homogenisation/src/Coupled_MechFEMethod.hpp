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

#ifndef __COUPLED_MECHFEMETHOD_HPP__
#define __COUPLED_MECHFEMETHOD_HPP__

#include "ElasticFEMethod.hpp"
//#include <boost/numeric/ublas/symmetric.hpp>
//extern "C" {
//#include <gm_rule.h>
//}
using namespace ObosleteUsersModules;

namespace MoFEM {
  
  struct Coupled_MechFEMethod: public ElasticFEMethod {
    double young;
    double pois;

    Coupled_MechFEMethod( FieldInterface& _mField,Mat &_Aij,Vec _X,Vec _F,double _lambda,double _mu):
    ElasticFEMethod(_mField,_Aij,_X,_F,_lambda,_mu) {};
    

    virtual PetscErrorCode calculateD_mech(double young, double nu, double conc_gauss) {
      PetscFunctionBegin;
      D.resize(6,6);
      D.clear();
      
      //Assuming linear relation between moisture and E
      young=young-3000*conc_gauss;
//      cout<<"conc_gauss "<<conc_gauss<<endl;
//      cout<<"young "<<young<<endl<<endl;

      double D00,D01,D33,constt;
      constt=young/((1+nu)*(1-2*nu));

      D00=constt*(1-nu);
      D01=constt*nu;
      D33=constt*(1-2*nu)/2;
      
      D(0,0)=D00;  D(0,1)=D01;  D(0,2)=D01;
      D(1,0)=D01;  D(1,1)=D00;  D(1,2)=D01;
      D(2,0)=D01;  D(2,1)=D01;  D(2,2)=D00;
      D(3,3)=D33;
      D(4,4)=D33;
      D(5,5)=D33;
      //      cout<<"D = "<<D;
      PetscFunctionReturn(0);
    }

    virtual PetscErrorCode GetMatParameters_mech(double *_young,double *_pois) {
      PetscFunctionBegin;
      *_young = young;
      *_pois = pois;
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
      PetscFunctionReturn(0);
    }

    
    virtual PetscErrorCode Stiffness() {
      PetscFunctionBegin;
      double _young,_pois;
      ierr = GetMatParameters_mech(&_young,&_pois); CHKERRQ(ierr);
//      cout<<"_young "<<_young<<endl;
//      cout<<"_pois "<<_pois<<endl;

      //moisture concentration at gasus point
      vector< ublas::vector< FieldData > > conc;
      ierr = GetGaussDataVector("CONC",conc); CHKERRQ(ierr);
//      cout<<"conc[0] "<<conc[0]<<endl;

      double conc_gauss;
      K.resize(row_mat,col_mat);
      int g_dim = g_NTET.size()/4;
//      cout<<"g_dim  = " <<g_dim<<endl;
      for(int rr = 0;rr<row_mat;rr++) {
        if(RowGlob[rr].size()==0) continue;
        for(int gg = 0;gg<g_dim;gg++) {
          conc_gauss=(conc[gg])[0];
//          cout<<"conc_gauss "<<conc_gauss<<endl;

          ierr = calculateD_mech(_young,_pois,conc_gauss); CHKERRQ(ierr);
//          cout<<"D "<<D<<endl;

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
//      cout<<"K(rr,cc) "<<K(0,0)<<endl;
//      cout<<"end of kstiffness"<<endl;
      PetscFunctionReturn(0);
    }
    
    
    virtual PetscErrorCode Fint() {
      PetscFunctionBegin;
      try {
        //Higher order approximation of geometry
        ierr = GetHierarchicalGeometryApproximation(invH,detH); CHKERRQ(ierr);
        
        double _young,_pois;
        ierr = GetMatParameters_mech(&_young,&_pois); CHKERRQ(ierr);
        
        //moisture concentration at gasus point
        vector< ublas::vector< FieldData > > conc;
        ierr = GetGaussDataVector("CONC",conc); CHKERRQ(ierr);
        double conc_gauss;
//      cout<<"conc[0] "<<conc[0]<<endl;
        
        //Gradient at Gauss points;
        vector< ublas::matrix< FieldData > > GradU_at_GaussPt;
        ierr = GetGaussDiffDataVector("DISPLACEMENT",GradU_at_GaussPt); CHKERRQ(ierr);
        unsigned int g_dim = g_NTET.size()/4;
        assert(GradU_at_GaussPt.size() == g_dim);
        NOT_USED(g_dim);
        vector< ublas::matrix< FieldData > >::iterator viit = GradU_at_GaussPt.begin();
        int gg = 0;
        for(;viit!=GradU_at_GaussPt.end();viit++,gg++) {
          try {
            ublas::matrix< FieldData > GradU = *viit;
//            cout<<"GradU "<<GradU<<endl;
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
            conc_gauss=(conc[gg])[0];
            ierr = calculateD_mech(_young,_pois,conc_gauss); CHKERRQ(ierr);
            //          cout<<"D "<<D<<endl;
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

#endif //__COUPLED_MECHFEMETHOD_HPP__
