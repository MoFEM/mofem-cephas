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

#ifndef __BCS_RVE_HOMOGENIZED_STRESS_DISP_HPP__
#define __BCS_RVE_HOMOGENIZED_STRESS_DISP_HPP__

namespace MoFEM {
  struct BCs_RVE_Homogenized_Stress_Disp: public BCs_RVELagrange_Disp {
    
    
    FieldInterface &mField;
    BCs_RVE_Homogenized_Stress_Disp(FieldInterface &m_field):BCs_RVELagrange_Disp(m_field),
    mField(m_field){}
    
    
    /// \biref operator to calculate the RVE homogenised stress
    struct OpRVEHomoStress:public FaceElementForcesAndSourcesCore::UserDataOperator {
      RVEBC_Data &dAta;
      bool ho_geometry;
      Vec Stress_Homo;
      
      OpRVEHomoStress(const string field_name, const string lagrang_field_name, Vec _Stress_Homo, RVEBC_Data &data, bool _ho_geometry = false):
      FaceElementForcesAndSourcesCore::UserDataOperator(lagrang_field_name, UserDataOperator::OPROW),Stress_Homo(_Stress_Homo), dAta(data), ho_geometry(_ho_geometry){}
      
      /** brief calculate Convection condition on the right hand side
       *  R=int_S N^T*alpha*N_f  dS **/
      PetscErrorCode doWork(int side,EntityType type,DataForcesAndSurcesCore::EntData &data) {
        PetscFunctionBegin;
        if(data.getIndices().size()==0) PetscFunctionReturn(0);
        if(dAta.tRis.find(getMoFEMFEPtr()->get_ent())==dAta.tRis.end()) PetscFunctionReturn(0);
        //        cout<<"OpRVEBCsRhs "<<endl;
        
        PetscErrorCode ierr;
        const FENumeredDofMoFEMEntity *dof_ptr;
        ierr = getMoFEMFEPtr()->get_row_dofs_by_petsc_gloabl_dof_idx(data.getIndices()[0],&dof_ptr); CHKERRQ(ierr);
        int rank = dof_ptr->get_max_rank();
        int nb_row_dofs = data.getIndices().size()/rank;
        
        ublas::matrix<FieldData> X_mat, N_mat;
        X_mat.resize(rank,1.5*rank+1.5);    X_mat.clear();
        ublas::matrix<FieldData>  D_mat;
        ublas::vector<FieldData> Lamda;
        ublas::vector<FieldData>  Stress_Homo_elem;
        Stress_Homo_elem.resize(1.5*rank+1.5);   Stress_Homo_elem.clear();   //homogenised stress for one element (triangle)
        
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
          
          int shape_size=data.getN().size2();
          N_mat.resize(rank,shape_size*rank);   N_mat.clear();
          int gg1=0;       int kk=0;
          for(int ii=0; ii<shape_size; ii++){ //number of shape funcitons
            for(int jj=0; jj<rank; jj++){
              //              cout<<"ii "<<ii<<endl;
              N_mat(jj,kk)=data.getN()(gg,gg1); kk++;
            }
            gg1++;
          }
          //          cout<<"N_mat "<<N_mat<<endl;
          if(gg==0){
            D_mat=val*prod(trans(N_mat),X_mat);
          }else{
            D_mat+=val*prod(trans(N_mat),X_mat);
          }
        }
        
        //        cout<<"data.getFieldData() "<<data.getFieldData()<<endl;
        Stress_Homo_elem=prod(trans(D_mat), -1*data.getFieldData());   //Lamda=data.getFieldData() is reaction force (so multiply for -1 to get the force)
        //        cout<<"Stress_Homo_elem "<<Stress_Homo_elem<<endl;
        int Indices6[6]={0, 1, 2, 3, 4, 5};
        int Indices3[3]={0, 1, 2};
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
        feRVEBCRhs.getOpPtrVector().push_back(new OpRVEHomoStress(field_name, lagrang_field_name, Stress_Homo, sit->second, ho_geometry));
        
      }
      PetscFunctionReturn(0);
    }
    
    


    
  };
}

#endif
