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

#ifndef __BCS_RVELAGRANGE_DISP_HPP
#define __BCS_RVELAGRANGE_DISP_HPP

namespace MoFEM {
  
  struct BCs_RVELagrange_Disp {
    
    FieldInterface &mField;
    
    struct MyTriFE: public FaceElementForcesAndSourcesCore {
      MyTriFE(FieldInterface &_mField): FaceElementForcesAndSourcesCore(_mField) {}
      int getRule(int order) { return order; };
    };

    MyTriFE feRVEBCRhs; //To calculate the Rhs or RVE BCs
    MyTriFE feRVEBCLhs; //To calculate the Lhs or RVE BCs
    
    MyTriFE& getLoopFeRVEBCRhs() { return feRVEBCRhs; }
    MyTriFE& getLoopFeRVEBCLhs() { return feRVEBCLhs; }


    BCs_RVELagrange_Disp(FieldInterface &m_field):
    feRVEBCRhs(m_field),feRVEBCLhs(m_field),
    mField(m_field) {}

    struct RVEBC_Data {
      Range tRis; // All boundary surfaces
    };
    map<int,RVEBC_Data> setOfRVEBC; ///< maps side set id with appropriate FluxData
    
    

    
  
    
    PetscErrorCode addLagrangiangElement(const string element_name, const string field_name, const string lagrang_field_name, const string mesh_nodals_positions) {
      PetscFunctionBegin;
      
      PetscErrorCode ierr;
      ErrorCode rval;
      
      ierr = mField.add_finite_element(element_name,MF_ZERO); CHKERRQ(ierr);
      //============================================================================================================
      //C row as Lagrange_mul_disp and col as DISPLACEMENT
      ierr = mField.modify_finite_element_add_field_row(element_name,lagrang_field_name); CHKERRQ(ierr);
      ierr = mField.modify_finite_element_add_field_col(element_name,field_name); CHKERRQ(ierr);
      //CT col as Lagrange_mul_disp and row as DISPLACEMENT
      ierr = mField.modify_finite_element_add_field_col(element_name,lagrang_field_name); CHKERRQ(ierr);
      ierr = mField.modify_finite_element_add_field_row(element_name,field_name); CHKERRQ(ierr);
      //data
      ierr = mField.modify_finite_element_add_field_data(element_name,lagrang_field_name); CHKERRQ(ierr);
      ierr = mField.modify_finite_element_add_field_data(element_name,field_name); CHKERRQ(ierr);
      //============================================================================================================

      if(mField.check_field(mesh_nodals_positions)) { //for high order geometry
        ierr = mField.modify_finite_element_add_field_data(element_name,mesh_nodals_positions); CHKERRQ(ierr);
      }
      //this is alternative method for setting boundary conditions, to bypass bu in cubit file reader.
      //not elegant, but good enough
      for(_IT_CUBITMESHSETS_BY_SET_TYPE_FOR_LOOP_(mField,SIDESET,it)) {
        if(it->get_name().compare(0,12,"AllBoundSurf") == 0) {
//          cout<<"Hi from addLagrangiangElement "<<endl;
          rval = mField.get_moab().get_entities_by_type(it->meshset,MBTRI,setOfRVEBC[it->get_msId()].tRis,true); CHKERR_PETSC(rval);
          ierr = mField.add_ents_to_finite_element_by_TRIs(setOfRVEBC[it->get_msId()].tRis,element_name); CHKERRQ(ierr);
          
        }
      }
       PetscFunctionReturn(0);
    }

    
    
    
    /// \biref operator to calculate convection therms on body surface and assemble to lhs of equations
    struct OpRVEBCsLhs:public FaceElementForcesAndSourcesCore::UserDataOperator {
      
      RVEBC_Data &dAta;
      bool ho_geometry;
      Mat Aij;

      OpRVEBCsLhs(const string field_name, const string lagrang_field_name, Mat _Aij,
                      RVEBC_Data &data,bool _ho_geometry = false):
      FaceElementForcesAndSourcesCore::UserDataOperator(field_name,lagrang_field_name),Aij(_Aij),
      dAta(data),ho_geometry(_ho_geometry){}
      
      ublas::matrix<double> K,transK;
      /** \brief calculate thermal convection term in the lhs of equations
       *
       * K = intS N^T alpha N dS
       */
      PetscErrorCode doWork(
                            int row_side,int col_side,
                            EntityType row_type,EntityType col_type,
                            DataForcesAndSurcesCore::EntData &row_data,
                            DataForcesAndSurcesCore::EntData &col_data) {
        PetscFunctionBegin;

        try {
//          cout<<"OpRVEBCsLhs "<<endl;
          if(row_data.getIndices().size()==0) PetscFunctionReturn(0);
          if(col_data.getIndices().size()==0) PetscFunctionReturn(0);
          
          int nb_row = row_data.getN().size2();
          int nb_col = col_data.getN().size2();
          K.resize(nb_row,nb_col);
          K.clear();
//          cout<<"no of Gauss points per element = "<<row_data.getN().size1()<<endl;
          
          for(unsigned int gg = 0;gg<row_data.getN().size1();gg++) {
            double area;
            if(ho_geometry) {
              area = norm_2(getNormals_at_GaussPt(gg))*0.5;
//              cout<<"area with ho_geometry = "<<area<<endl;
            }   else {
              area = getArea();
//              cout<<"area without ho_geometry = "<<area<<endl;
            }
            double val = getGaussPts()(2,gg)*area;
            noalias(K) += val*outer_prod( row_data.getN(gg,nb_row),col_data.getN(gg,nb_col) );
            
          }
          
          
//          cout<<"&row_data.getIndices()[0] "<< endl;
//            for (int ii = nb_row,ii<nb_row; ii++) row_data.getIndices()[ii];
//          cout<<"&col_data.getIndices()[0] "<<&col_data.getIndices()[0]<<endl;

          PetscErrorCode ierr;
          ierr = MatSetValues(
                              Aij,
                              nb_row,&row_data.getIndices()[0],
                              nb_col,&col_data.getIndices()[0],
                              &K(0,0),ADD_VALUES); CHKERRQ(ierr);
          if(row_side != col_side || row_type != col_type) {
            transK.resize(nb_col,nb_row);
            noalias(transK) = trans( K );
            ierr = MatSetValues(
                                Aij,
                                nb_col,&col_data.getIndices()[0],
                                nb_row,&row_data.getIndices()[0],
                                &transK(0,0),ADD_VALUES); CHKERRQ(ierr);
          }
        } catch (const std::exception& ex) {
          ostringstream ss;
          ss << "throw in method: " << ex.what() << endl;
          SETERRQ(PETSC_COMM_SELF,1,ss.str().c_str());
        }
        
        PetscFunctionReturn(0);
      }
    };

    
    
    
    
    
    
    
    
    
    
    PetscErrorCode setRVEBCsOperators(string field_name,string lagrang_field_name,Mat _Aij, const string mesh_nodals_positions) {
      PetscFunctionBegin;
      
      bool ho_geometry = false;
      if(mField.check_field(mesh_nodals_positions)) {
        ho_geometry = true;
      }
//      cout<<"Hi 1 from setRVEBCsOperators "<<endl;
      map<int,RVEBC_Data>::iterator sit = setOfRVEBC.begin();
      for(;sit!=setOfRVEBC.end();sit++) {
//        cout<<"Hi from setRVEBCsOperators "<<endl;
        feRVEBCLhs.getRowColOpPtrVector().push_back(new OpRVEBCsLhs(field_name,lagrang_field_name, _Aij, sit->second,ho_geometry));
      }
      PetscFunctionReturn(0);
    }
    
    
    
    
    
    
  };

}

#endif
