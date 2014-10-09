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

#ifndef __COUPLED_MECHMOISTUREFEMETHOD_HPP__
#define __COUPLED_MECHMOISTUREFEMETHOD_HPP__

#include <boost/numeric/ublas/symmetric.hpp>

namespace MoFEM {
  
  struct Coupled_MechMoistureFEMethod: public FEMethod_UpLevelStudent {
    
    FieldInterface& mField;
    Mat Aij;
    Vec Data,F;
    double young;
    double pois;
    
    Coupled_MechMoistureFEMethod(FieldInterface& _mField,Mat &_Aij,Vec &_D,Vec& _F):FEMethod_UpLevelStudent(_mField.get_moab(),1), mField(_mField),Aij(_Aij),Data(_D),F(_F){
      
      snes_B = &Aij;
      snes_x = Data;
      snes_f = F;

      pcomm = ParallelComm::get_pcomm(&moab,MYPCOMM_INDEX);
      
      RowGlob.resize(1+6+4+1);
      RowLocal.resize(1+6+4+1);
      rowNMatrices.resize(1+6+4+1);
      rowDiffNMatrices.resize(1+6+4+1);
      rowBMatrices.resize(1+6+4+1);
      
      ColGlob.resize(1+6+4+1);
      ColLocal.resize(1+6+4+1);
      colNMatrices.resize(1+6+4+1);
      colDiffNMatrices.resize(1+6+4+1);
      
      g_NTET.resize(4*45);
      ShapeMBTET(&g_NTET[0],G_TET_X45,G_TET_Y45,G_TET_Z45,45);
      //      cout<<"g_NTET[0] "<<g_NTET[0]<<endl;
      G_W_TET = G_TET_W45;
      
      dE_dc=-3000;
      beta=0.0001;
//      dE_dc=0;
//      beta=0;

    };
    
    ErrorCode rval;
    ParallelComm* pcomm;
    PetscLogDouble t1,t2;
    PetscLogDouble v1,v2;
    
    int row_mat,col_mat;
    vector<vector<DofIdx> > RowGlob;
    vector<vector<DofIdx> > RowLocal;
    vector<vector<ublas::matrix<FieldData> > > rowNMatrices;
    vector<vector<ublas::matrix<FieldData> > > rowDiffNMatrices;
    vector<vector<ublas::matrix<FieldData> > > rowBMatrices;
    
    vector<vector<DofIdx> > ColGlob;
    vector<vector<DofIdx> > ColLocal;
    vector<vector<ublas::matrix<FieldData> > > colNMatrices;
    vector<vector<ublas::matrix<FieldData> > > colDiffNMatrices;
    
    vector<double> g_NTET,g_NTRI;
    const double* G_W_TET;
    const double* G_W_TRI;
    double dE_dc;
    double beta;
    
    
    //******************************************************************************************
    PetscErrorCode preProcess() {
      PetscFunctionBegin;
//      PetscSynchronizedPrintf(PETSC_COMM_WORLD,"Start Assembly\n");
//      ierr = PetscTime(&v1); CHKERRQ(ierr);
//      ierr = PetscGetCPUTime(&t1); CHKERRQ(ierr);
      PetscFunctionReturn(0);
    }
    
    //******************************************************************************************
    PetscErrorCode postProcess() {
      PetscFunctionBegin;
      
      switch(snes_ctx) {
        case CTX_SNESNONE: {
          // Note MAT_FLUSH_ASSEMBLY
          ierr = MatAssemblyBegin(*snes_B,MAT_FLUSH_ASSEMBLY); CHKERRQ(ierr);
          ierr = MatAssemblyEnd(*snes_B,MAT_FLUSH_ASSEMBLY); CHKERRQ(ierr);
          ierr = VecAssemblyBegin(snes_f); CHKERRQ(ierr);
          ierr = VecAssemblyEnd(snes_f); CHKERRQ(ierr);
        }
          break;
        case CTX_SNESSETFUNCTION: {
          ierr = VecAssemblyBegin(snes_f); CHKERRQ(ierr);
          ierr = VecAssemblyEnd(snes_f); CHKERRQ(ierr);
        }
          break;
        case CTX_SNESSETJACOBIAN: {
          // Note MAT_FLUSH_ASSEMBLY
          ierr = MatAssemblyBegin(*snes_B,MAT_FLUSH_ASSEMBLY); CHKERRQ(ierr);
          ierr = MatAssemblyEnd(*snes_B,MAT_FLUSH_ASSEMBLY); CHKERRQ(ierr);
        }
          break;
        default:
          SETERRQ(PETSC_COMM_SELF,1,"not implemented");
      }
      
      PetscFunctionReturn(0);
    }

    
    
    //******************************************************************************************
    ublas::matrix<FieldData> D;
    virtual PetscErrorCode calculateD_mech(double young, double nu, double conc_gauss) {
      PetscFunctionBegin;
      D.resize(6,6);
      D.clear();
      
      //Assuming linear relation between moisture and E
      young=young+dE_dc*conc_gauss;
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
    
    //******************************************************************************************
    ublas::matrix<FieldData> dD_dE;
    virtual PetscErrorCode calculatedD_dE(double young, double nu) {
      PetscFunctionBegin;
      dD_dE.resize(6,6);
      dD_dE.clear();
      
      double D00,D01,D33,constt;
      constt=1/((1+nu)*(1-2*nu));
      
      D00=constt*(1-nu);
      D01=constt*nu;
      D33=constt*(1-2*nu)/2;
      
      dD_dE(0,0)=D00;  dD_dE(0,1)=D01;  dD_dE(0,2)=D01;
      dD_dE(1,0)=D01;  dD_dE(1,1)=D00;  dD_dE(1,2)=D01;
      dD_dE(2,0)=D01;  dD_dE(2,1)=D01;  dD_dE(2,2)=D00;
      dD_dE(3,3)=D33;
      dD_dE(4,4)=D33;
      dD_dE(5,5)=D33;
      //      cout<<"D = "<<D;
      PetscFunctionReturn(0);
    }
    
    //******************************************************************************************
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
    
    //******************************************************************************************
    virtual PetscErrorCode GetMatricesRows() {
      PetscFunctionBegin;
      //indicies ROWS
      //cout<<"GetMatricesRows " <<endl;
      row_mat = 0;
      ierr = GetRowGlobalIndices("DISPLACEMENT",RowGlob[row_mat]); CHKERRQ(ierr);
      ierr = GetRowLocalIndices("DISPLACEMENT",RowLocal[row_mat]); CHKERRQ(ierr);
      ierr = GetGaussRowNMatrix("DISPLACEMENT",rowNMatrices[row_mat]); CHKERRQ(ierr);
      ierr = GetGaussRowDiffNMatrix("DISPLACEMENT",rowDiffNMatrices[row_mat]); CHKERRQ(ierr);
      ierr = MakeBMatrix3D("DISPLACEMENT",rowDiffNMatrices[row_mat],rowBMatrices[row_mat]);  CHKERRQ(ierr);
      //      cout<<"RowLocal "<< (RowLocal[0])[0]<<"\t"<<(RowLocal[0])[1]<<"\t"<<(RowLocal[0])[2]<<"\t"<<(RowLocal[0])[3]<<endl;
      //      cout<<"RowGlob  "<< (RowGlob[0]) [0]<<"\t"<<(RowGlob[0]) [1]<<"\t"<<(RowGlob[0]) [2]<<"\t"<<(RowGlob[0]) [3]<<endl;
      
      row_mat++;
      for(int ee = 0;ee<6;ee++) { //edges matrices
        ierr = GetRowGlobalIndices("DISPLACEMENT",MBEDGE,RowGlob[row_mat],ee); CHKERRQ(ierr);
        ierr = GetRowLocalIndices("DISPLACEMENT",MBEDGE,RowLocal[row_mat],ee); CHKERRQ(ierr);
        if(RowGlob[row_mat].size()!=0) {
          ierr = GetGaussRowNMatrix("DISPLACEMENT",MBEDGE,rowNMatrices[row_mat],ee); CHKERRQ(ierr);
          ierr = GetGaussRowDiffNMatrix("DISPLACEMENT",MBEDGE,rowDiffNMatrices[row_mat],ee); CHKERRQ(ierr);
          ierr = MakeBMatrix3D("DISPLACEMENT",rowDiffNMatrices[row_mat],rowBMatrices[row_mat]);  CHKERRQ(ierr);
          row_mat++;
        }
      }
      for(int ff = 0;ff<4;ff++) { //faces matrices
        ierr = GetRowGlobalIndices("DISPLACEMENT",MBTRI,RowGlob[row_mat],ff); CHKERRQ(ierr);
        ierr = GetRowLocalIndices("DISPLACEMENT",MBTRI,RowLocal[row_mat],ff); CHKERRQ(ierr);
        if(RowGlob[row_mat].size()!=0) {
          ierr = GetGaussRowNMatrix("DISPLACEMENT",MBTRI,rowNMatrices[row_mat],ff); CHKERRQ(ierr);
          ierr = GetGaussRowDiffNMatrix("DISPLACEMENT",MBTRI,rowDiffNMatrices[row_mat],ff); CHKERRQ(ierr);
          ierr = MakeBMatrix3D("DISPLACEMENT",rowDiffNMatrices[row_mat],rowBMatrices[row_mat]);  CHKERRQ(ierr);
          row_mat++;
        }
      }
      ierr = GetRowGlobalIndices("DISPLACEMENT",MBTET,RowGlob[row_mat]); CHKERRQ(ierr);
      ierr = GetRowLocalIndices("DISPLACEMENT",MBTET,RowLocal[row_mat]); CHKERRQ(ierr);
      if(RowGlob[row_mat].size() != 0) { //volume matrices
        ierr = GetGaussRowNMatrix("DISPLACEMENT",MBTET,rowNMatrices[row_mat]); CHKERRQ(ierr);
        ierr = GetGaussRowDiffNMatrix("DISPLACEMENT",MBTET,rowDiffNMatrices[row_mat]); CHKERRQ(ierr);
        ierr = MakeBMatrix3D("DISPLACEMENT",rowDiffNMatrices[row_mat],rowBMatrices[row_mat]);  CHKERRQ(ierr);
        row_mat++;
      }
      PetscFunctionReturn(0);
    }
    
    //******************************************************************************************
    virtual PetscErrorCode GetMatricesCols() {
      PetscFunctionBegin;
      //indicies COLS
      col_mat = 0;
      ierr = GetColGlobalIndices("CONC",ColGlob[col_mat]); CHKERRQ(ierr);
      ierr = GetColLocalIndices("CONC",ColLocal[col_mat]); CHKERRQ(ierr);
      ierr = GetGaussColNMatrix("CONC",colNMatrices[col_mat]); CHKERRQ(ierr);
      ierr = GetGaussColDiffNMatrix("CONC",colDiffNMatrices[col_mat]); CHKERRQ(ierr);
      //      cout<<"colNMatrices row "<< (colNMatrices[col_mat])[0] <<endl;
      col_mat++;
      for(int ee = 0;ee<6;ee++) { //edges matrices
        ierr = GetColGlobalIndices("CONC",MBEDGE,ColGlob[col_mat],ee); CHKERRQ(ierr);
        ierr = GetColLocalIndices("CONC",MBEDGE,ColLocal[col_mat],ee); CHKERRQ(ierr);
        if(ColGlob[col_mat].size()!=0) {
          ierr = GetGaussColNMatrix("CONC",MBEDGE,colNMatrices[col_mat],ee); CHKERRQ(ierr);
          ierr = GetGaussColDiffNMatrix("CONC",MBEDGE,colDiffNMatrices[col_mat],ee); CHKERRQ(ierr);
          //          cout<<"colNMatrices row "<< (colNMatrices[col_mat])[0] <<endl;
          col_mat++;
        }
      }
      for(int ff = 0;ff<4;ff++) { //faces matrices
        ierr = GetColGlobalIndices("CONC",MBTRI,ColGlob[col_mat],ff); CHKERRQ(ierr);
        ierr = GetColLocalIndices("CONC",MBTRI,ColLocal[col_mat],ff); CHKERRQ(ierr);
        if(ColGlob[col_mat].size()!=0) {
          ierr = GetGaussColNMatrix("CONC",MBTRI,colNMatrices[col_mat],ff); CHKERRQ(ierr);
          ierr = GetGaussColDiffNMatrix("CONC",MBTRI,colDiffNMatrices[col_mat],ff); CHKERRQ(ierr);
          col_mat++;
        }
      }
      ierr = GetColGlobalIndices("CONC",MBTET,ColGlob[col_mat]); CHKERRQ(ierr);
      ierr = GetColLocalIndices("CONC",MBTET,ColLocal[col_mat]); CHKERRQ(ierr);
      if(ColGlob[col_mat].size() != 0) { //volume matrices
        ierr = GetGaussColNMatrix("CONC",MBTET,colNMatrices[col_mat]); CHKERRQ(ierr);
        ierr = GetGaussColDiffNMatrix("CONC",MBTET,colDiffNMatrices[col_mat]); CHKERRQ(ierr);
        col_mat++;
      }
      PetscFunctionReturn(0);
    }
    
    virtual PetscErrorCode GetMatrices() {
      PetscFunctionBegin;
      ierr = GetMatricesRows(); CHKERRQ(ierr);
      ierr = GetMatricesCols(); CHKERRQ(ierr);
      PetscFunctionReturn(0);
    }
    
    /******************************************************************************************
     Kuc calculation
     ******************************************************************************************/
    ublas::matrix< FieldData > I;
    ublas::matrix<ublas::matrix<FieldData> > K;
    ublas::matrix<FieldData> row_Mat;
    virtual PetscErrorCode Stiffness() {
      PetscFunctionBegin;
      I.resize(6,1); I.clear();
      I(0,0) = 1.0;  I(0,1) = 1.0;  I(0,2) = 1.0;
      
      double _young,_pois;
      ierr = GetMatParameters_mech(&_young,&_pois); CHKERRQ(ierr);
      //      cout<< "Young's modulus "<<_young<<endl;
      //      cout<< "Poisson's ratio "<<_pois<<endl;
      
      ierr = calculatedD_dE(_young, _pois); CHKERRQ(ierr);
      //      cout<< "dD_dE "<<dD_dE<<endl;
      
      //moisture concentration at gasus point
      vector< ublas::vector< FieldData > > conc;
      ierr = GetGaussDataVector("CONC",conc); CHKERRQ(ierr);
//      cout<<"conc[0] "<<conc[0]<<endl;
      
      //calculate strain
      //Gradient at Gauss points;
      vector< ublas::matrix< FieldData > > GradU_at_GaussPt;
      ierr = GetGaussDiffDataVector("DISPLACEMENT",GradU_at_GaussPt); CHKERRQ(ierr);
      vector< ublas::matrix< FieldData > >::iterator viit = GradU_at_GaussPt.begin();
      //      cout<<"row_mat "<<row_mat<<endl;
      //      cout<<"col_mat "<<col_mat<<endl;
      K.resize(row_mat,col_mat);
      int g_dim = g_NTET.size()/4;
      for(int gg = 0;gg<g_dim;gg++,viit++) {
        //strain claculation
        ublas::matrix< FieldData > GradU = *viit;
//        cout<<"GradU "<<GradU<<endl;

        ublas::matrix< FieldData > Strain = 0.5*( GradU + trans(GradU) );
        ublas::matrix< FieldData > VoightStrain(6,1);
        VoightStrain(0,0) = Strain(0,0);
        VoightStrain(0,1) = Strain(1,1);
        VoightStrain(0,2) = Strain(2,2);
        VoightStrain(0,3) = 2*Strain(0,1);
        VoightStrain(0,4) = 2*Strain(1,2);
        VoightStrain(0,5) = 2*Strain(2,0);
        
        double conc_gauss=(conc[gg])[0];
        ierr = calculateD_mech(_young,_pois,conc_gauss); CHKERRQ(ierr);
        //        cout<<"D "<<D<<endl;
        //          cout<<"VoightStrain "<<VoightStrain<<endl;
        for(int rr = 0;rr<row_mat;rr++) {
          ublas::matrix<FieldData> &Bmat = (rowBMatrices[rr])[gg];
          //          cout<<"(rowBMatrices[rr])[gg]"<<(rowBMatrices[rr])[gg]<<endl;
          
//          cout<<"prod(D,beta*I) "<<prod(D,beta*I)<<endl;
//          cout<<"prod(dD_dE*dE_dc,(VoightStrain-beta*conc_gauss*I)) "<<prod(dD_dE*dE_dc,(VoightStrain-beta*conc_gauss*I))<<endl;
          ublas::matrix<FieldData> AA=prod(dD_dE*dE_dc,(VoightStrain-beta*conc_gauss*I))-prod(D,beta*I);
          double w = V*G_W_TET[gg];
          row_Mat = prod(w*trans(Bmat),AA);
          
          for(int cc = 0;cc<col_mat;cc++) {
            ublas::matrix<FieldData> &col_Mat=(colNMatrices[cc])[gg];
            if(gg == 0) {
              K(rr,cc).resize(row_Mat.size1(),col_Mat.size2());
              K(rr,cc).clear();
            }
            K(rr,cc)+=prod(row_Mat,col_Mat);
//            cout<<"K(rr,cc) "<<K(rr,cc)<<endl;
          }
        }
      }
      PetscFunctionReturn(0);
    }
  //******************************************************************************************
    virtual PetscErrorCode Lhs() {
      PetscFunctionBegin;
      ierr = Stiffness(); CHKERRQ(ierr);
      for(int rr = 0;rr<row_mat;rr++) {
        if(RowGlob[rr].size()==0) continue;
        for(int cc = 0;cc<col_mat;cc++) {
          if(ColGlob[cc].size()==0) continue;
          if(RowGlob[rr].size()!=K(rr,cc).size1()) SETERRQ(PETSC_COMM_SELF,1,"data inconsistency");
          if(ColGlob[cc].size()!=K(rr,cc).size2()) SETERRQ(PETSC_COMM_SELF,1,"data inconsistency");
          ierr = MatSetValues(*snes_B,RowGlob[rr].size(),&(RowGlob[rr])[0],ColGlob[cc].size(),&(ColGlob[cc])[0],&(K(rr,cc).data())[0],ADD_VALUES); CHKERRQ(ierr);
        }
      }
      PetscFunctionReturn(0);
    }
  //******************************************************************************************
    virtual PetscErrorCode Rhs() {
      PetscFunctionBegin;
      ierr = Fint(snes_f); CHKERRQ(ierr);
      PetscFunctionReturn(0);
    }

     //******************************************************************************************
    vector<ublas::vector<FieldData> > f_int;
    virtual PetscErrorCode Fint() {
      PetscFunctionBegin;
      try {
        ublas::vector< FieldData >  Ivec;  Ivec.resize(6);  Ivec.clear();
        Ivec[0] = 1;  Ivec[1] = 1;  Ivec[2] = 1;
//        cout<<"Ivec "<<Ivec<<endl;
        double _young,_pois;
        ierr = GetMatParameters_mech(&_young,&_pois); CHKERRQ(ierr);
        
        //moisture concentration at gasus point
        vector< ublas::vector< FieldData > > conc;
        ierr = GetGaussDataVector("CONC",conc); CHKERRQ(ierr);
        double conc_gauss;
//        cout<<"conc[0] "<<conc[0]<<endl;

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
            ublas::matrix< FieldData > Strain = 0.5*( GradU + trans(GradU) );
//            cout<<"GradU "<<GradU<<endl;

            ublas::vector< FieldData > VoightStrain(6);
            VoightStrain[0] = Strain(0,0);
            VoightStrain[1] = Strain(1,1);
            VoightStrain[2] = Strain(2,2);
            VoightStrain[3] = 2*Strain(0,1);
            VoightStrain[4] = 2*Strain(1,2);
            VoightStrain[5] = 2*Strain(2,0);
            
            double conc_gauss=(conc[gg])[0];
            ierr = calculateD_mech(_young,_pois,conc_gauss); CHKERRQ(ierr);
            //        cout<<"D "<<D<<endl;
            
            ublas::vector< FieldData > StrianSwell=beta*Ivec*conc_gauss;
//            cout<<"VoightStrain "<<VoightStrain<<endl;
//            cout<<"StrianSwell "<<StrianSwell<<endl;
            
            VoightStrain=VoightStrain-StrianSwell;
            
            double w = V*G_W_TET[gg];
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
   
    
     //******************************************************************************************
    virtual PetscErrorCode Fint(Vec F_int) {
      PetscFunctionBegin;
      try {
        ierr = Fint(); CHKERRQ(ierr);
        for(int rr = 0;rr<row_mat;rr++) {
          if(RowGlob[rr].size()==0) continue;
//          cout<<"f_int[rr] "<<f_int[rr]<<endl;
          
          f_int[rr]*=1;
          
          if(RowGlob[rr].size()!=f_int[rr].size()) SETERRQ(PETSC_COMM_SELF,1,"data inconsistency");
          ierr = VecSetValues(F_int,RowGlob[rr].size(),&(RowGlob[rr])[0],&(f_int[rr].data()[0]),ADD_VALUES); CHKERRQ(ierr);
        }
      } catch (const std::exception& ex) {
        ostringstream ss;
        ss << "thorw in method: " << ex.what() << " at line " << __LINE__ << " in file " << __FILE__ << endl;
        SETERRQ(PETSC_COMM_SELF,1,ss.str().c_str());
      }
      PetscFunctionReturn(0);
    }
    //******************************************************************************************
    PetscErrorCode operator()() {
      PetscFunctionBegin;
      ierr = OpStudentStart_TET(g_NTET); CHKERRQ(ierr);
//      cout<<"Hi from Coupled_MechMoistureFEMethod Class "<<endl;
      ierr = GetMatrices(); CHKERRQ(ierr);
      
      switch(snes_ctx) {
        case CTX_SNESNONE: {
//          cout<<"CTX_SNESNONE "<<endl;
          ierr = Rhs(); CHKERRQ(ierr);
          ierr = Lhs(); CHKERRQ(ierr);
        }
        break;
        case CTX_SNESSETFUNCTION: {
//          cout<<"CTX_SNESSETFUNCTION "<<endl;
          ierr = Rhs(); CHKERRQ(ierr);
        }
        break;
        case CTX_SNESSETJACOBIAN: {
//          cout<<"CTX_SNESSETJACOBIAN "<<endl;
          ierr = Lhs(); CHKERRQ(ierr);
        }
        break;
        default:
          SETERRQ(PETSC_COMM_SELF,1,"not implemented");
      }
      //      std::string wait;
      //      std::cin >> wait;
      ierr = OpStudentEnd(); CHKERRQ(ierr);
      PetscFunctionReturn(0);
    }
//******************************************************************************************

  };
  
  
}

#endif //__COUPLED_MECHMOISTUREFEMETHOD_HPP__
