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

#ifndef __BCS_RVE_HOMOGENIZED_STRESS_TRAC_HPP
#define __BCS_RVE_HOMOGENIZED_STRESS_TRAC_HPP

namespace MoFEM {
  struct BCs_RVE_Homogenized_Stress_Trac: public BCs_RVELagrange_Trac {
    FieldInterface &mField;
    BCs_RVE_Homogenized_Stress_Trac(FieldInterface &m_field):BCs_RVELagrange_Trac(m_field),
    mField(m_field){}

    /// \biref operator to calculate the RVE homogenised stress
    struct OpRVEHomoStress_Cal:public FaceElementForcesAndSourcesCore::UserDataOperator {
      RVEBC_Data &dAta;
      bool ho_geometry;
      Vec Stress_Homo;
      CommonData &commonData;
      CommonFunctions common_functions;

      OpRVEHomoStress_Cal(const string field_name, RVEBC_Data &data, CommonData &_commonData, CommonFunctions &_common_functions, bool _ho_geometry = false):
      FaceElementForcesAndSourcesCore::UserDataOperator(field_name, UserDataOperator::OPCOL), dAta(data), commonData(_commonData), common_functions(_common_functions), ho_geometry(_ho_geometry){}

      /** brief calculate Convection condition on the right hand side
       *  R=int_S N^T*alpha*N_f  dS **/
      PetscErrorCode doWork(int side,EntityType type,DataForcesAndSurcesCore::EntData &data) {
        PetscFunctionBegin;
        if(data.getIndices().size()==0) PetscFunctionReturn(0);
        if(type!=MBVERTEX) PetscFunctionReturn(0);
        if(dAta.tRis.find(getMoFEMFEPtr()->get_ent())==dAta.tRis.end()) PetscFunctionReturn(0);
//        cout<<"OpRVEBCsRhs "<<endl;

        PetscErrorCode ierr;
        const FENumeredDofMoFEMEntity *dof_ptr;
        ierr = getMoFEMFEPtr()->get_row_dofs_by_petsc_gloabl_dof_idx(data.getIndices()[0],&dof_ptr); CHKERRQ(ierr);
        int rank = dof_ptr->get_max_rank();
        commonData.rank=rank;
        ublas::matrix<FieldData> X_mat, N_mat, NTX, H_mat;
        X_mat.resize(rank,1.5*rank+1.5);    X_mat.clear();

        double x,y,z;
        for(unsigned int gg = 0;gg<data.getN().size1();gg++) {
          double area;
          if(ho_geometry) {
            area = norm_2(getNormals_at_GaussPt(gg))*0.5;
          } else {
            area = getArea();
          }
          double val = getGaussPts()(2,gg)*area;

          x = getCoordsAtGaussPts()(gg,0);
          y = getCoordsAtGaussPts()(gg,1);
          z = getCoordsAtGaussPts()(gg,2);
          //          cout<<"Gauss point coordinates "<< "  x=    "<< x << "   y=   "<< y << "  z=    "<<z<<endl;
          
          switch(rank) {
            case 3: //mech problem
              X_mat(0,0)=2.0*x;  X_mat(0,3)=y;  X_mat(0,5)=z;
              X_mat(1,1)=2.0*y;  X_mat(1,3)=x;  X_mat(1,4)=z;
              X_mat(2,2)=2.0*z;  X_mat(2,4)=y;  X_mat(2,5)=x;
              X_mat=0.5*X_mat;
              break;
            case 1:  //moisture transport or thermal problem
              X_mat(0,0)=x;  X_mat(0,1)=y;  X_mat(0,2)=z;
              break;
            default:
              SETERRQ(PETSC_COMM_SELF,1,"not implemented");
          }
          ierr = common_functions.shapeMat(rank, gg,  data, N_mat); CHKERRQ(ierr);
//          cout<<"N_mat "<<N_mat<<endl;
          ierr = common_functions.HMat(getNormals_at_GaussPt(gg), rank, N_mat, H_mat); CHKERRQ(ierr);
//          cout<<"H_mat "<<H_mat<<endl;
          if(gg==0){
            NTX=val*prod(trans(N_mat),X_mat);
            commonData.D_mat=prod(H_mat, NTX);
          }else{
            NTX=val*prod(trans(N_mat),X_mat);
            commonData.D_mat+=prod(H_mat, NTX);
          }
        }
//        cout<<"commonData.D_mat "<<commonData.D_mat<<endl;

        PetscFunctionReturn(0);
      }
    };
    
    /// \biref operator to calculate the RVE homogenised stress
    struct OpRVEHomoStress_Assemble:public FaceElementForcesAndSourcesCore::UserDataOperator {
      RVEBC_Data &dAta;
      bool ho_geometry;
      Vec Stress_Homo;
      CommonData &commonData;
      
      OpRVEHomoStress_Assemble(const string lagrang_field_name, Vec _Stress_Homo, RVEBC_Data &data, CommonData &_commonData, bool _ho_geometry = false):
      FaceElementForcesAndSourcesCore::UserDataOperator(lagrang_field_name, UserDataOperator::OPROW), dAta(data), commonData(_commonData), Stress_Homo(_Stress_Homo), ho_geometry(_ho_geometry){}
      
      /** brief calculate Convection condition on the right hand side
       *  R=int_S N^T*alpha*N_f  dS **/
      PetscErrorCode doWork(int side,EntityType type,DataForcesAndSurcesCore::EntData &data) {
        PetscFunctionBegin;
        if(data.getIndices().size()==0) PetscFunctionReturn(0);
        if(dAta.tRis.find(getMoFEMFEPtr()->get_ent())==dAta.tRis.end()) PetscFunctionReturn(0);
//        cout<<"OpRVEBCsRhs "<<endl;
        ublas::vector<FieldData> Lamda;
        ublas::vector<FieldData>  Stress_Homo_elem;
        
        int rank=commonData.rank;
        Stress_Homo_elem.resize(1.5*rank+1.5);   Stress_Homo_elem.clear();   //homogenised stress for one element (triangle)
        
//        cout<<"commonData.D_mat "<<commonData.D_mat<<endl;
//        cout<<"data.getFieldData() "<<data.getFieldData()<<endl;

        Stress_Homo_elem=prod(trans(commonData.D_mat), -1*data.getFieldData());   //Lamda=data.getFieldData() is reaction force (so multiply for -1 to get the force)
        commonData.D_mat.clear(); 
//        cout<<"Stress_Homo_elem "<<Stress_Homo_elem<<endl;
        int Indices6[6]={0, 1, 2, 3, 4, 5};
        int Indices3[3]={0, 1, 2};
        PetscErrorCode ierr;

        switch(rank) {
          case 3:
            ierr = VecSetValues(Stress_Homo,6,Indices6,&(Stress_Homo_elem.data())[0],ADD_VALUES); CHKERRQ(ierr);
            break;
          case 1:
            ierr = VecSetValues(Stress_Homo,3,Indices3,&(Stress_Homo_elem.data())[0],ADD_VALUES); CHKERRQ(ierr);
            break;
          default:
            SETERRQ(PETSC_COMM_SELF,1,"not implemented");
        }

      PetscFunctionReturn(0);
      }
  };

    
    
    PetscErrorCode setRVEBCsHomoStressOperators(string field_name,string lagrang_field_name, Vec Stress_Homo, map<int,RVEBC_Data> &setOfRVEBC, const string mesh_nodals_positions) {
      PetscFunctionBegin;
      
      bool ho_geometry = false;
      if(mField.check_field(mesh_nodals_positions)) {
        ho_geometry = true;
      }
      
      //      cout<<"Hi from setRVEBCsHomoStressOperators "<<endl;
      map<int,RVEBC_Data>::iterator sit = setOfRVEBC.begin();
      for(;sit!=setOfRVEBC.end();sit++) {
        //        cout<<"Hi from setOfRVEBC "<<endl;
//        commonData.D_mat.resize(6,6);  commonData.D_mat.clear();
        feRVEBCRhs.getOpPtrVector().push_back(new OpRVEHomoStress_Cal(field_name, sit->second, commonData, common_functions, ho_geometry));
        feRVEBCRhs.getOpPtrVector().push_back(new OpRVEHomoStress_Assemble(lagrang_field_name, Stress_Homo, sit->second, commonData, ho_geometry));

      }
      PetscFunctionReturn(0);
    }
    
    
  };
}

#endif
