/* Copyright (C) 2015, Zahur Ullah (Zahur.Ullah AT glasgow.ac.uk)
 * --------------------------------------------------------------
 *
 * Description: Implementation of thermal stress, i.e. right hand side as result of thermal stresses
 *
 * This is not exactly procedure for linear elatic dynamics, since jacobian is
 * evaluated at every time step and snes procedure is involved. However it is
 * implemented like that, to test methodology for general nonlinear problem.
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

#ifndef __ELASTIC_ELEMENT_HPP
#define __ELASTIC_ELEMENT_HPP

namespace MoFEM {
  struct ElasticElement {
    
    struct MyVolumeFE: public TetElementForcesAndSourcesCore {
      MyVolumeFE(FieldInterface &_mField): TetElementForcesAndSourcesCore(_mField) {}
      int getRule(int order) { return order-1; };
    };
    
    
    MyVolumeFE feRhs; ///< cauclate right hand side for tetrahedral elements
    MyVolumeFE& getLoopFeRhs() { return feRhs; } ///< get rhs volume element
    MyVolumeFE feLhs; //< calculate left hand side for tetrahedral elements
    MyVolumeFE& getLoopFeLhs() { return feLhs; } ///< get lhs volume element

    struct MyTriFE: public TriElementForcesAndSurcesCore {
      MyTriFE(FieldInterface &_mField): TriElementForcesAndSurcesCore(_mField) {}
      int getRule(int order) { return order; };
    };
    MyTriFE feFlux; //< heat flux element
    MyTriFE& getLoopFeFlux() { return feFlux; } //< get heat flux element

    FieldInterface &mField;
    ElasticElement(FieldInterface &m_field):
    feRhs(m_field),feLhs(m_field),
    feFlux(m_field),
    mField(m_field) {}

    
    struct BlockData {
      double yOung_modulus;
      double pOisson_ratio;
      Range tEts; ///< constatins elements in block set
    };
    map<int,BlockData> setOfBlocks; ///< maps block set id with appropiate BlockData

    
    struct CommonData {
      ublas::vector<double> wtAtGaussPts;  //Here we will get Wt which will used to Dmat= Wt*Dmat0
      ublas::matrix<double> gradAtGaussPts;

    };
    CommonData commonData;

    
    PetscErrorCode addElasticElements(const string field_name,const string mesh_nodals_positions = "MESH_NODE_POSITIONS") {
      PetscFunctionBegin;
      
      PetscErrorCode ierr;
      ErrorCode rval;
      
      ierr = mField.add_finite_element("ELASTIC_FE_MACRO",MF_ZERO); CHKERRQ(ierr);
      ierr = mField.modify_finite_element_add_field_row("ELASTIC_FE_MACRO",field_name); CHKERRQ(ierr);
      ierr = mField.modify_finite_element_add_field_col("ELASTIC_FE_MACRO",field_name); CHKERRQ(ierr);
      ierr = mField.modify_finite_element_add_field_data("ELASTIC_FE_MACRO",field_name); CHKERRQ(ierr);
      if(mField.check_field(mesh_nodals_positions)) {
        ierr = mField.modify_finite_element_add_field_data("ELASTIC_FE_MACRO",mesh_nodals_positions); CHKERRQ(ierr);
      }
      
      // loop over all blocksets and get data which name is MAT_ELASTICSET
      for(_IT_CUBITMESHSETS_BY_BCDATA_TYPE_FOR_LOOP_(mField,BLOCKSET|MAT_ELASTICSET,it)) {
        cout<<"MAT_ELASTICSET  "<<endl;
        Mat_Elastic elastic_data;
        ierr = it->get_attribute_data_structure(elastic_data); CHKERRQ(ierr);
        setOfBlocks[it->get_msId()].yOung_modulus = elastic_data.data.Young;
        setOfBlocks[it->get_msId()].pOisson_ratio = elastic_data.data.Poisson;
        rval = mField.get_moab().get_entities_by_type(it->meshset,MBTET,setOfBlocks[it->get_msId()].tEts,true); CHKERR_PETSC(rval);
        ierr = mField.add_ents_to_finite_element_by_TETs(setOfBlocks[it->get_msId()].tEts,"ELASTIC_FE_MACRO"); CHKERRQ(ierr);
      }
      
      PetscFunctionReturn(0);
    }

    
    
    
    template<typename OP>
    struct OpGetFieldAtGaussPts: public OP::UserDataOperator {
      
      ublas::vector<double> &fieldAtGaussPts;
      OpGetFieldAtGaussPts(const string field_name,ublas::vector<double> &field_at_gauss_pts):
      OP::UserDataOperator(field_name),
      fieldAtGaussPts(field_at_gauss_pts) {}
      
      /** \brief operator calculating temperature and rate of temperature
       *
       * temperature temperature or rate of temperature is calculated multiplyingshape functions by degrees of freedom
       */
      PetscErrorCode doWork(
                            int side,EntityType type,DataForcesAndSurcesCore::EntData &data) {
        PetscFunctionBegin;
        try {
//          cout<<"form OpGetFieldAtGaussPts "<<endl;
          
          if(data.getFieldData().size()==0) PetscFunctionReturn(0);
          int nb_dofs = data.getFieldData().size();
          int nb_gauss_pts = data.getN().size1();
          
          //initialize
          fieldAtGaussPts.resize(nb_gauss_pts);
          if(type == MBVERTEX) {
            //loop over shape functions on entities allways start from
            //vertices, so if nodal shape functions are processed, vector of
            //field values is zeroad at initialization
            fill(fieldAtGaussPts.begin(),fieldAtGaussPts.end(),0);
          }
          
          for(int gg = 0;gg<nb_gauss_pts;gg++) {
            fieldAtGaussPts[gg] += inner_prod(data.getN(gg,nb_dofs),data.getFieldData());
            
//            cout<<"fieldAtGaussPts[gg] "<<fieldAtGaussPts[gg] <<endl;
            
          }
          
        } catch (const std::exception& ex) {
          ostringstream ss;
          ss << "throw in method: " << ex.what() << endl;
          SETERRQ(PETSC_COMM_SELF,MOFEM_STD_EXCEPTION_THROW,ss.str().c_str());
        }
        
        PetscFunctionReturn(0);
      }
      
    };

    
    struct OpGetWtAtGaussPts: public OpGetFieldAtGaussPts<TetElementForcesAndSourcesCore> {
      OpGetWtAtGaussPts(const string wt_field_name,CommonData &common_data):
      OpGetFieldAtGaussPts<TetElementForcesAndSourcesCore>(wt_field_name,common_data.wtAtGaussPts) {}
    };

    
    
    
    struct OpGetGradAtGaussPts: public TetElementForcesAndSourcesCore::UserDataOperator {
      
      CommonData &commonData;
      OpGetGradAtGaussPts(const string field_name,CommonData &common_data):
      TetElementForcesAndSourcesCore::UserDataOperator(field_name),
      commonData(common_data) {}
      
      /** \brief operator calculating temperature gradients
       *
       * temerature gradient is calculated multiplying derivatives of shape functions by degrees of freedom.
       */
      PetscErrorCode doWork(
                            int side,EntityType type,DataForcesAndSurcesCore::EntData &data) {
        PetscFunctionBegin;
        try {

          if(data.getIndices().size()==0) PetscFunctionReturn(0);
          
          PetscErrorCode ierr;
          const FENumeredDofMoFEMEntity *dof_ptr;
          ierr = getMoFEMFEPtr()->get_row_dofs_by_petsc_gloabl_dof_idx(data.getIndices()[0],&dof_ptr); CHKERRQ(ierr);

          int rank = dof_ptr->get_max_rank();
          int nb_dofs = data.getFieldData().size()/rank;
          int nb_gauss_pts = data.getN().size1();
          
//          cout<<" rank "<< rank<<endl;
//          cout<<" nb_dofs "<< nb_dofs<<endl;
//          cout<<" nb_gauss_pts "<< nb_gauss_pts<<endl;

          //for each gauss point
          //[dux/dx  duy/dx  dyz/dx]
          //[dux/dy  duy/dy  dyz/dy]
          //[dux/dz  duy/dz  dyz/dz]
          
          //initialize
          commonData.gradAtGaussPts.resize(3*nb_gauss_pts,3);
          if(type == MBVERTEX) {
            fill(commonData.gradAtGaussPts.data().begin(),commonData.gradAtGaussPts.data().end(),0);
          }


          for(int gg = 0;gg<nb_gauss_pts;gg++) {
            cout<<"Hi from gauss points  "<<endl;
            cout<<"data.getDiffN(gg,nb_dofs)  "<<data.getDiffN(gg,nb_dofs)<<endl;
            cout<<" data.getFieldData()  "<< data.getFieldData()<<endl;
            //for each gauss point
            //[u1x  u1y  u1z]
            //[u2x  u2y  u2z]
            //[u3x  u3y  u3z]
            //[u4x  u4y  u4z]
            ublas::vector<double> uvec=data.getFieldData();
            ublas::matrix<double> umat;
            umat.resize(4,3);
            umat(0,0)=uvec(0); umat(0,1)=uvec(1); umat(0,2)=uvec(2);
            umat(1,0)=uvec(3); umat(1,1)=uvec(4); umat(1,2)=uvec(5);
            umat(2,0)=uvec(6); umat(2,1)=uvec(7); umat(2,2)=uvec(8);
            umat(3,0)=uvec(9); umat(3,1)=uvec(10); umat(3,2)=uvec(11);
            
            cout<<"umat  "<<umat<<endl;
            cout<<"prod(trans(data.getDiffN(gg,nb_dofs)), data.getFieldData() )  "<<prod(trans(data.getDiffN(gg,nb_dofs)), data.getFieldData() )<<endl;

            

            
////            ublas::noalias(commonData.getGradAtGaussPts(gg)) += prod( trans(data.getDiffN(gg,nb_dofs)), data.getFieldData() );
            string wait;
            cin >>wait;

          }
          
        } catch (const std::exception& ex) {
          ostringstream ss;
          ss << "throw in method: " << ex.what() << endl;
          SETERRQ(PETSC_COMM_SELF,1,ss.str().c_str());
        }
        
        PetscFunctionReturn(0);
      }
      
    };

    
    
    
    
    
    
    PetscErrorCode setElasticFiniteElementRhsOperators(string field_name,string wt_field_name, Vec &F) {
      PetscFunctionBegin;
      feRhs.get_op_to_do_Rhs().push_back(new OpGetWtAtGaussPts(wt_field_name,commonData));
      feRhs.get_op_to_do_Rhs().push_back(new OpGetGradAtGaussPts(field_name,commonData));

//      map<int,BlockData>::iterator sit = setOfBlocks.begin();
//      feRhs.get_op_to_do_Rhs().push_back(new OpGetGradAtGaussPts(field_name,commonData));
//      for(;sit!=setOfBlocks.end();sit++) {
//        //add finite element
//        feRhs.get_op_to_do_Rhs().push_back(new OpThermalRhs(field_name,F,sit->second,commonData));
//      }
      PetscFunctionReturn(0);
    }

    
    
    
    PetscErrorCode setElasticFiniteElementLhsOperators(string field_name,string wt_field_name, Mat A) {
      PetscFunctionBegin;
      
//      feLhs.get_op_to_do_Lhs().push_back(new OpGetWtAtGaussPts(wt_field_name,commonData));

//      map<int,BlockData>::iterator sit = setOfBlocks.begin();
//      for(;sit!=setOfBlocks.end();sit++) {
        //add finite elemen
//        feLhs.get_op_to_do_Lhs().push_back(new OpThermalLhs(field_name,A,sit->second,commonData));
//      }
      PetscFunctionReturn(0);
    }

    
    
    
  };
  
}

#endif //__ELASTIC_ELEMENT_HPP



