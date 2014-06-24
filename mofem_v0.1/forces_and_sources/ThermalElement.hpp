/* Copyright (C) 2013, Lukasz Kaczmarczyk (likask AT wp.pl)
 * --------------------------------------------------------------
 *
 * Test for linar elastic dynamics.
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


#ifndef __BODY_FORCE_HPP
#define __BODY_FORCE_HPP

#include "ForcesAndSurcesCore.hpp"

namespace MoFEM {

struct ThermalElement {

  struct MyVolumeFE: public TetElementForcesAndSurcesCore {
    MyVolumeFE(FieldInterface &_mField): TetElementForcesAndSurcesCore(_mField) {}
    int getRule(int order) { return order-1; };
  };

  MyVolumeFE feRhs;
  MyVolumeFE& getLoopFeRhs() { return feRhs; }
  MyVolumeFE feLhs;
  MyVolumeFE& getLoopFeLhs() { return feLhs; }

  FieldInterface &mField;
  ThermalElement(
    FieldInterface &m_field):
    feRhs(m_field),feLhs(m_field),mField(m_field) {}


  struct BlockData {
    double cOnductivity;
    Range tEts;
  };
  map<int,BlockData> setOfBlocks;

  struct CommonData {
    ublas::vector<double> temperatureAtGaussPts;
    ublas::matrix<double> gradAtGaussPts;
  };
  CommonData commonData;

  struct OpGetGradAtGaussPts: public TetElementForcesAndSurcesCore::UserDataOperator {

    CommonData &commonData;
    OpGetGradAtGaussPts(const string field_name,CommonData &common_data):
      TetElementForcesAndSurcesCore::UserDataOperator(field_name),
      commonData(common_data) {}

    PetscErrorCode doWork(
      int side,EntityType type,DataForcesAndSurcesCore::EntData &data) {
      PetscFunctionBegin;
      try {

      if(data.getIndices().size()==0) PetscFunctionReturn(0);

      commonData.gradAtGaussPts.resize(data.getN().size1(),3);

      switch(type) {
	case MBVERTEX:
	for(int dd = 0;dd<3;dd++) {
	  commonData.gradAtGaussPts(0,dd) = cblas_ddot(4,&data.getDiffN()(0,dd),3,&data.getFieldData()[0],1);
	  for(unsigned int gg = 1;gg<data.getN().size1();gg++) {
	    commonData.gradAtGaussPts(gg,dd) = commonData.gradAtGaussPts(gg,dd); 
	  }
	}
	break;
	default:
	for(unsigned int gg = 0;gg<data.getN().size1();gg++) {
	  for(int dd = 0;dd<3;dd++) {
	   commonData.gradAtGaussPts(gg,dd) += cblas_ddot(data.getN().size2(),&data.getDiffN()(gg,dd),3,&data.getFieldData()[0],1);
	  }
	}
      }

      } catch (const std::exception& ex) {
	ostringstream ss;
	ss << "throw in method: " << ex.what() << endl;
	SETERRQ(PETSC_COMM_SELF,1,ss.str().c_str());
      }

      PetscFunctionReturn(0);
    }

  };

  struct OpGetTemperatureAtGaussPts:  public TetElementForcesAndSurcesCore::UserDataOperator {

    CommonData &commonData;
    OpGetTemperatureAtGaussPts(const string field_name,CommonData &common_data):
      TetElementForcesAndSurcesCore::UserDataOperator(field_name),
      commonData(common_data) {}

    PetscErrorCode doWork(
      int side,EntityType type,DataForcesAndSurcesCore::EntData &data) {
      PetscFunctionBegin;

      try {

      if(data.getIndices().size()==0) PetscFunctionReturn(0);

      commonData.temperatureAtGaussPts.resize(data.getN().size1());
      for(unsigned int gg = 0;gg<data.getN().size1();gg++) {
	commonData.temperatureAtGaussPts(gg) += cblas_ddot(data.getN().size2(),&data.getN()(gg,0),1,&data.getFieldData()[0],1);
      }

      } catch (const std::exception& ex) {
	ostringstream ss;
	ss << "throw in method: " << ex.what() << endl;
	SETERRQ(PETSC_COMM_SELF,1,ss.str().c_str());
      }

      PetscFunctionReturn(0);
    }

  };

  struct OpThernalRhs: public TetElementForcesAndSurcesCore::UserDataOperator {

    Vec &F;
    BlockData &dAta;
    CommonData &commonData;
    OpThernalRhs(const string field_name,Vec &_F,BlockData &data,CommonData &common_data):
      TetElementForcesAndSurcesCore::UserDataOperator(field_name),
      F(_F),dAta(data),commonData(common_data) {}

    ublas::vector<double> Nf;
    PetscErrorCode doWork(
      int side,EntityType type,DataForcesAndSurcesCore::EntData &data) {
      PetscFunctionBegin;

      try {

      if(data.getIndices().size()==0) PetscFunctionReturn(0);
      if(dAta.tEts.find(getMoFEMFEPtr()->get_ent())==dAta.tEts.end()) PetscFunctionReturn(0);

      PetscErrorCode ierr;

      int nb_row_dofs = data.getIndices().size();
      Nf.resize(data.getIndices().size());
      bzero(&*Nf.data().begin(),data.getIndices().size()*sizeof(FieldData));

      for(unsigned int gg = 0;gg<data.getN().size1();gg++) {

	double val = -dAta.cOnductivity*getVolume()*getGaussPts()(3,gg)*getVolume();
	switch(type) {
	  case MBVERTEX:
	    cblas_dgemv(CblasRowMajor,CblasNoTrans,nb_row_dofs,3,val,
	      &data.getDiffN()(0,0),3,&commonData.gradAtGaussPts(gg,0),1,
	      1.,&Nf[0],1);
	  break;
	  default:
	    cblas_dgemv(CblasRowMajor,CblasNoTrans,nb_row_dofs,3,val,
	      &data.getDiffN()(gg,0),3,&commonData.gradAtGaussPts(gg,0),1,
	      1.,&Nf[0],1);
	}
      }

      ierr = VecSetValues(F,data.getIndices().size(),
	&data.getIndices()[0],&Nf[0],ADD_VALUES); CHKERRQ(ierr);

      } catch (const std::exception& ex) {
	ostringstream ss;
	ss << "throw in method: " << ex.what() << endl;
	SETERRQ(PETSC_COMM_SELF,1,ss.str().c_str());
      }

      PetscFunctionReturn(0);
    }

  };

  struct OpThernalLhs: public TetElementForcesAndSurcesCore::UserDataOperator {

    Mat &A;
    BlockData &dAta;
    CommonData &commonData;
    OpThernalLhs(const string field_name,Mat &_A,BlockData &data,CommonData &common_data):
      TetElementForcesAndSurcesCore::UserDataOperator(field_name),
      A(_A),dAta(data),commonData(common_data) {}

    ublas::matrix<double> K,transK;
    PetscErrorCode doWork(
      int row_side,int col_side,
      EntityType row_type,EntityType col_type,
      DataForcesAndSurcesCore::EntData &row_data,
      DataForcesAndSurcesCore::EntData &col_data) {
      PetscFunctionBegin;

      try {
  
      if(row_data.getIndices().size()==0) PetscFunctionReturn(0);
      if(col_data.getIndices().size()==0) PetscFunctionReturn(0);
  
      int nb_row = row_data.getN().size2();
      int nb_col = col_data.getN().size2();
      K.resize(nb_row,nb_col);
      bzero(&*K.data().begin(),nb_row*nb_col*sizeof(double));
  
      for(unsigned int gg = 0;gg<row_data.getN().size1();gg++) {

	double *diff_N_row,*diff_N_col;

	if(row_type==MBVERTEX) {
	  diff_N_row = &row_data.getDiffN()(0,0);
	} else {
	  diff_N_row = &row_data.getDiffN()(gg,0);
	}

	if(col_type==MBVERTEX) {
	  diff_N_col = &col_data.getDiffN()(0,0);
	} else {
	  diff_N_col = &col_data.getDiffN()(gg,0);
	}
  

        double val = -dAta.cOnductivity*getVolume()*getGaussPts()(3,gg)*getVolume();
        cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasTrans,
	  nb_row,nb_col,3,
	  val,diff_N_row,3,diff_N_col,3,1.,&K(0,0),nb_col);
  
      }
  
      PetscErrorCode ierr;
      ierr = MatSetValues(
        A,
        nb_row,&row_data.getIndices()[0],
        nb_col,&col_data.getIndices()[0],
        &K(0,0),ADD_VALUES); CHKERRQ(ierr);
      if(row_side != col_side || row_type != col_type) {
	transK.resize(nb_col,nb_row);
	noalias(transK) = trans( K );
	ierr = MatSetValues(
	  A,
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

  PetscErrorCode addThermalElements(
    const string problem_name,const string field_name,const string mesh_nodals_positions = "MESH_NODE_POSITIONS") {
    PetscFunctionBegin;

    PetscErrorCode ierr;
    ErrorCode rval;

    ierr = mField.add_finite_element("THERMAL_FE",MF_ZERO); CHKERRQ(ierr);
    ierr = mField.modify_finite_element_add_field_row("THERMAL_FE",field_name); CHKERRQ(ierr);
    ierr = mField.modify_finite_element_add_field_col("THERMAL_FE",field_name); CHKERRQ(ierr);
    ierr = mField.modify_finite_element_add_field_data("THERMAL_FE",field_name); CHKERRQ(ierr);
    if(mField.check_field(mesh_nodals_positions)) {
      ierr = mField.modify_finite_element_add_field_data("THERMAL_FE",mesh_nodals_positions); CHKERRQ(ierr);
    }
    ierr = mField.modify_problem_add_finite_element(problem_name,"THERMAL_FE"); CHKERRQ(ierr);

    //takes skin of block of entities
    Skinner skin(&mField.get_moab());
    // loop over all blocksets and get data which name is FluidPressure
    for(_IT_CUBITMESHSETS_BY_BCDATA_TYPE_FOR_LOOP_(mField,BlockSet|Mat_ThermalSet,it)) {

      Mat_Thermal temp_data;
      ierr = it->get_attribute_data_structure(temp_data); CHKERRQ(ierr);  
      setOfBlocks[it->get_msId()].cOnductivity = temp_data.data.Conductivity;
      rval = mField.get_moab().get_entities_by_type(it->meshset,MBTET,setOfBlocks[it->get_msId()].tEts,true); CHKERR_PETSC(rval);
      ierr = mField.add_ents_to_finite_element_by_TETs(setOfBlocks[it->get_msId()].tEts,"THERMAL_FE"); CHKERRQ(ierr);

    }

    PetscFunctionReturn(0);
  }

  PetscErrorCode setThermalFiniteElementRhsOperators(string field_name,Vec &F) {
    PetscFunctionBegin;
    map<int,BlockData>::iterator sit = setOfBlocks.begin();
    for(;sit!=setOfBlocks.end();sit++) {
      //add finite element
      feRhs.get_op_to_do_Rhs().push_back(new OpGetGradAtGaussPts(field_name,commonData));
      feRhs.get_op_to_do_Rhs().push_back(new OpThernalRhs(field_name,F,sit->second,commonData));
    }
    PetscFunctionReturn(0);
  }


  PetscErrorCode setThermalFiniteElementLhsOperators(string field_name,Mat &A) {
    PetscFunctionBegin;
    map<int,BlockData>::iterator sit = setOfBlocks.begin();
    for(;sit!=setOfBlocks.end();sit++) {
      //add finite element
      feLhs.get_op_to_do_Lhs().push_back(new OpThernalLhs(field_name,A,sit->second,commonData));
    }
    PetscFunctionReturn(0);
  }

};

}

#endif //__BODY_FORCE_HPP

