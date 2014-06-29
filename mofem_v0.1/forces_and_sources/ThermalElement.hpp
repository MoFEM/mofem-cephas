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
#include "TsCtx.hpp"

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

  struct MyTriFE: public TriElementForcesAndSurcesCore {
    MyTriFE(FieldInterface &_mField): TriElementForcesAndSurcesCore(_mField) {}
    int getRule(int order) { return ceil(order/2); };
  };
 
  MyTriFE feFlux;
  MyTriFE& getLoopFeFlux() { return feFlux; }

  FieldInterface &mField;
  ThermalElement(
    FieldInterface &m_field):
    feRhs(m_field),feLhs(m_field),feFlux(m_field),mField(m_field) {}


  struct BlockData {
    double cOnductivity;
    double cApacity;
    Range tEts;
  };
  map<int,BlockData> setOfBlocks;

  struct FluxData {
    heatflux_cubit_bc_data dAta;
    Range tRis;
  };
  map<int,FluxData> setOfFluxes;

  struct CommonData {
    ublas::vector<double> temperatureAtGaussPts;
    ublas::vector<double> temperatureRateAtGaussPts;
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
	int nb_dof = data.getFieldData().size();
  
        switch(type) {
	  case MBVERTEX:
	  for(int dd = 0;dd<3;dd++) {
	    commonData.gradAtGaussPts(0,dd) = cblas_ddot(4,&data.getDiffN()(0,dd),3,&data.getFieldData()[0],1);
	    for(unsigned int gg = 1;gg<data.getN().size1();gg++) {
	      commonData.gradAtGaussPts(gg,dd) = commonData.gradAtGaussPts(0,dd); 
	    }
	  }
	  break;
	  default:
	  for(unsigned int gg = 0;gg<data.getN().size1();gg++) {
	    for(int dd = 0;dd<3;dd++) {
	      commonData.gradAtGaussPts(gg,dd) += cblas_ddot(nb_dof,&data.getDiffN()(gg,dd),3,&data.getFieldData()[0],1);
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

  struct OpGetRateAtGaussPts: public TetElementForcesAndSurcesCore::UserDataOperator {

    CommonData &commonData;
    OpGetRateAtGaussPts(const string field_name,CommonData &common_data):
      TetElementForcesAndSurcesCore::UserDataOperator(field_name),
      commonData(common_data) {}

    PetscErrorCode doWork(
      int side,EntityType type,DataForcesAndSurcesCore::EntData &data) {
      PetscFunctionBegin;
      try {

        if(data.getFieldData().size()==0) PetscFunctionReturn(0);
        commonData.temperatureRateAtGaussPts.resize(data.getN().size1());
	int nb_dof = data.getFieldData().size();

        switch(type) {
	  case MBVERTEX:
	  for(unsigned int gg = 0;gg<data.getN().size1();gg++) {
	    commonData.temperatureRateAtGaussPts[gg] 
	      = cblas_ddot(nb_dof,&data.getN()(gg,0),1,&data.getFieldData()[0],1);
	  }
	  break;
	  default:
	  for(unsigned int gg = 0;gg<data.getN().size1();gg++) {
	    commonData.temperatureRateAtGaussPts[gg] 
	      += cblas_ddot(nb_dof,&data.getN()(gg,0),1,&data.getFieldData()[0],1);
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


  struct OpThermalRhs: public TetElementForcesAndSurcesCore::UserDataOperator {

    Vec F;
    BlockData &dAta;
    CommonData &commonData;
    bool useTsF;
    OpThermalRhs(const string field_name,Vec _F,BlockData &data,CommonData &common_data):
      TetElementForcesAndSurcesCore::UserDataOperator(field_name),
      F(_F),dAta(data),commonData(common_data),useTsF(false) { }
    OpThermalRhs(const string field_name,BlockData &data,CommonData &common_data):
      TetElementForcesAndSurcesCore::UserDataOperator(field_name),
      F(PETSC_NULL),dAta(data),commonData(common_data),useTsF(true) { }

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

      //cerr << data.getIndices() << endl;
      //cerr << data.getDiffN() << endl;

      for(unsigned int gg = 0;gg<data.getN().size1();gg++) {

	double val = dAta.cOnductivity*getVolume()*getGaussPts()(3,gg);

	if(getHoGaussPtsDetJac().size()>0) {
	  val *= getHoGaussPtsDetJac()[gg]; ///< higher order geometry
	}

	//cerr << val << endl;
	//cerr << data.getDiffN() << endl;
	//cerr << data.getIndices() << endl;
	//cerr << commonData.gradAtGaussPts << endl;
	cblas_dgemv(CblasRowMajor,CblasNoTrans,nb_row_dofs,3,val,
	  &data.getDiffN()(gg,0),3,&commonData.gradAtGaussPts(gg,0),1,
	  1.,&Nf[0],1);

      }
      
      //cerr << Nf << endl;
      if(useTsF) {
	ierr = VecSetValues(getFEMethod()->ts_F,data.getIndices().size(),
	  &data.getIndices()[0],&Nf[0],ADD_VALUES); CHKERRQ(ierr);
      } else {
	ierr = VecSetValues(F,data.getIndices().size(),
	  &data.getIndices()[0],&Nf[0],ADD_VALUES); CHKERRQ(ierr);

      }

      } catch (const std::exception& ex) {
	ostringstream ss;
	ss << "throw in method: " << ex.what() << endl;
	SETERRQ(PETSC_COMM_SELF,1,ss.str().c_str());
      }

      PetscFunctionReturn(0);
    }

  };

  struct OpThermalLhs: public TetElementForcesAndSurcesCore::UserDataOperator {

    Mat *A;
    BlockData &dAta;
    CommonData &commonData;
    bool useTSB;
    OpThermalLhs(const string field_name,Mat *_A,BlockData &data,CommonData &common_data):
      TetElementForcesAndSurcesCore::UserDataOperator(field_name),
      A(_A),dAta(data),commonData(common_data),useTSB(false) { }
    OpThermalLhs(const string field_name,BlockData &data,CommonData &common_data):
      TetElementForcesAndSurcesCore::UserDataOperator(field_name),
      A(PETSC_NULL),dAta(data),commonData(common_data),useTSB(true) { }

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
	  diff_N_row = &row_data.getDiffN()(gg,0);
	  diff_N_col = &col_data.getDiffN()(gg,0);

	  double val = dAta.cOnductivity*getVolume()*getGaussPts()(3,gg);
	  if(getHoGaussPtsDetJac().size()>0) {
	    val *= getHoGaussPtsDetJac()[gg]; ///< higher order geometry
	  }
	  cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasTrans,
	    nb_row,nb_col,3,
	    val,diff_N_row,3,diff_N_col,3,1.,&K(0,0),nb_col);
  
	}
  
	Mat M;
	if(useTSB) {
	  M = *(getFEMethod()->ts_B);
	} else {
	  M = *A;
	}

	PetscErrorCode ierr;
	ierr = MatSetValues(
	  M,
	  nb_row,&row_data.getIndices()[0],
	  nb_col,&col_data.getIndices()[0],
	  &K(0,0),ADD_VALUES); CHKERRQ(ierr);
	if(row_side != col_side || row_type != col_type) {
	  transK.resize(nb_col,nb_row);
	  noalias(transK) = trans( K );
	  ierr = MatSetValues(
	    M,
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


  struct OpHeatCapacityRhs: public TetElementForcesAndSurcesCore::UserDataOperator {

    BlockData &dAta;
    CommonData &commonData;
    OpHeatCapacityRhs(const string field_name,BlockData &data,CommonData &common_data):
      TetElementForcesAndSurcesCore::UserDataOperator(field_name),
      dAta(data),commonData(common_data) {}

    ublas::vector<double> Nf;
    PetscErrorCode doWork(
      int side,EntityType type,DataForcesAndSurcesCore::EntData &data) {
      PetscFunctionBegin;

      try {
  
	if(data.getIndices().size()==0) PetscFunctionReturn(0);
	int nb_row = data.getN().size2();
	Nf.resize(nb_row);
	bzero(&Nf[0],nb_row*sizeof(double));
	for(unsigned int gg = 0;gg<data.getN().size1();gg++) {
	  double val = dAta.cApacity*getVolume()*getGaussPts()(3,gg);
	  if(getHoGaussPtsDetJac().size()>0) {
	    val *= getHoGaussPtsDetJac()[gg]; ///< higher order geometry
	  }
	  val *= commonData.temperatureRateAtGaussPts[gg];
	  cblas_daxpy(nb_row,val,&data.getN()(gg,0),1,&*Nf.data().begin(),1);
	}
	PetscErrorCode ierr;
	ierr = VecSetValues(getFEMethod()->ts_F,data.getIndices().size(),
	  &data.getIndices()[0],&Nf[0],ADD_VALUES); CHKERRQ(ierr);
      

      } catch (const std::exception& ex) {
	ostringstream ss;
	ss << "throw in method: " << ex.what() << endl;
	SETERRQ(PETSC_COMM_SELF,1,ss.str().c_str());
      }


      PetscFunctionReturn(0);
    }

  };

  struct OpHeatCapacityLsh: public TetElementForcesAndSurcesCore::UserDataOperator {

    BlockData &dAta;
    CommonData &commonData;
    OpHeatCapacityLsh(const string field_name,BlockData &data,CommonData &common_data):
      TetElementForcesAndSurcesCore::UserDataOperator(field_name),
      dAta(data),commonData(common_data) {}

    ublas::matrix<double> M,transM;
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
	M.resize(nb_row,nb_col);
	bzero(&*M.data().begin(),nb_row*nb_col*sizeof(double));
  
	for(unsigned int gg = 0;gg<row_data.getN().size1();gg++) {

	  double *N_row,*N_col;
	  N_row = &row_data.getN()(gg,0);
	  N_col = &col_data.getN()(gg,0);

	  double val = dAta.cApacity*getVolume()*getGaussPts()(3,gg);
	  if(getHoGaussPtsDetJac().size()>0) {
	    val *= getHoGaussPtsDetJac()[gg]; ///< higher order geometry
	  }
	  val *= getFEMethod()->ts_a;
	  cblas_dger(CblasRowMajor,
	    nb_row,nb_col,val,N_row,1,N_col,1,&M(0,0),nb_col);

	}
  
	PetscErrorCode ierr;
	ierr = MatSetValues(
	  *(getFEMethod()->ts_B),
	  nb_row,&row_data.getIndices()[0],
	  nb_col,&col_data.getIndices()[0],
	  &M(0,0),ADD_VALUES); CHKERRQ(ierr);
	if(row_side != col_side || row_type != col_type) {
	  transM.resize(nb_col,nb_row);
	  noalias(transM) = trans(M);
	  ierr = MatSetValues(
	    *(getFEMethod()->ts_B),
	    nb_col,&col_data.getIndices()[0],
	    nb_row,&row_data.getIndices()[0],
	    &transM(0,0),ADD_VALUES); CHKERRQ(ierr);
	}
      

      } catch (const std::exception& ex) {
	ostringstream ss;
	ss << "throw in method: " << ex.what() << endl;
	SETERRQ(PETSC_COMM_SELF,1,ss.str().c_str());
      }


      PetscFunctionReturn(0);
    }

  };

  struct OpHeatFlux:public TriElementForcesAndSurcesCore::UserDataOperator {

    Vec F;
    FluxData &dAta;
    bool ho_geometry;
    bool useTsF;
    OpHeatFlux(const string field_name,Vec _F,
      FluxData &data,bool _ho_geometry = false):
      TriElementForcesAndSurcesCore::UserDataOperator(field_name),
      F(_F),dAta(data),ho_geometry(_ho_geometry),useTsF(false) { }
    OpHeatFlux(const string field_name,FluxData &data,bool _ho_geometry = false):
      TriElementForcesAndSurcesCore::UserDataOperator(field_name),
      F(PETSC_NULL),dAta(data),ho_geometry(_ho_geometry),useTsF(true) { }

    ublas::vector<FieldData> Nf;

    PetscErrorCode doWork(
      int side,EntityType type,DataForcesAndSurcesCore::EntData &data) {
      PetscFunctionBegin;

      if(data.getIndices().size()==0) PetscFunctionReturn(0);
      if(dAta.tRis.find(getMoFEMFEPtr()->get_ent())==dAta.tRis.end()) PetscFunctionReturn(0);

      PetscErrorCode ierr;

      const FENumeredDofMoFEMEntity *dof_ptr;
      ierr = getMoFEMFEPtr()->get_row_dofs_by_petsc_gloabl_dof_idx(data.getIndices()[0],&dof_ptr); CHKERRQ(ierr);
      int rank = dof_ptr->get_max_rank();

      int nb_row_dofs = data.getIndices().size()/rank;
      
      Nf.resize(data.getIndices().size());
      bzero(&*Nf.data().begin(),data.getIndices().size()*sizeof(FieldData));

      //cerr << getNormal() << endl;
      //cerr << getNormals_at_GaussPt() << endl;

      for(unsigned int gg = 0;gg<data.getN().size1();gg++) {

	double val = getGaussPts()(2,gg);
	double flux;
	if(ho_geometry) {
	  double area = cblas_dnrm2(3,&getNormals_at_GaussPt()(gg,0),1);
	  flux = dAta.dAta.data.value1*area;
	} else {
	  flux = dAta.dAta.data.value1*getArea();
	}
	cblas_daxpy(nb_row_dofs,val*flux,&data.getN()(gg,0),1,&*Nf.data().begin(),1);

      }
    

    
      //cerr << "VecSetValues\n";
      //cerr << Nf << endl;
      //cerr << data.getIndices() << endl;

      if(useTsF) {
	ierr = VecSetValues(getFEMethod()->ts_F,data.getIndices().size(),
	  &data.getIndices()[0],&Nf[0],ADD_VALUES); CHKERRQ(ierr);
      } else {
	ierr = VecSetValues(F,data.getIndices().size(),
	  &data.getIndices()[0],&Nf[0],ADD_VALUES); CHKERRQ(ierr);
      }

      PetscFunctionReturn(0);
    }

  };

  struct UpdateAndControl: public FieldInterface::FEMethod {

    FieldInterface& mField;
    TS tS;
    const string tempName;
    const string rateName;
    int jacobianLag;
    UpdateAndControl(FieldInterface& _mField,TS _ts,
      const string temp_name,const string rate_name): mField(_mField),tS(_ts),
      tempName(temp_name),rateName(rate_name),jacobianLag(-1) {}

    PetscErrorCode preProcess() {
      PetscFunctionBegin;
      PetscErrorCode ierr;
      ierr = mField.set_other_local_VecCreateGhost(
	problem_ptr,tempName,rateName,Row,ts_u_t,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
      PetscFunctionReturn(0);
    }

    PetscErrorCode postProcess() {
      PetscFunctionBegin;
      PetscErrorCode ierr;
      SNES snes;
      ierr = TSGetSNES(tS,&snes); CHKERRQ(ierr);
      ierr = SNESSetLagJacobian(snes,jacobianLag); CHKERRQ(ierr);
      PetscFunctionReturn(0);
    }

  };

  struct TimeSeriesMonitor: public FieldInterface::FEMethod {

    FieldInterface &mField;
    const string seriesName;
    const string tempName;
    BitRefLevel mask;

    TimeSeriesMonitor(FieldInterface &m_field,const string series_name,const string temp_name):
      mField(m_field),seriesName(series_name),tempName(temp_name) {
      mask.set();
    }

    PetscErrorCode postProcess() {
      PetscFunctionBegin;
      PetscErrorCode ierr;

      ierr = mField.set_global_VecCreateGhost(
	problem_ptr,Row,ts_u,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);

      BitRefLevel proble_bit_level = problem_ptr->get_BitRefLevel();
      ierr = mField.record_begin(seriesName); CHKERRQ(ierr);
      ierr = mField.record_field(seriesName,tempName,proble_bit_level,mask); CHKERRQ(ierr);
      ierr = mField.record_end(seriesName); CHKERRQ(ierr);

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
      setOfBlocks[it->get_msId()].cApacity = temp_data.data.HeatCapacity;
      rval = mField.get_moab().get_entities_by_type(it->meshset,MBTET,setOfBlocks[it->get_msId()].tEts,true); CHKERR_PETSC(rval);
      ierr = mField.add_ents_to_finite_element_by_TETs(setOfBlocks[it->get_msId()].tEts,"THERMAL_FE"); CHKERRQ(ierr);

    }

    PetscFunctionReturn(0);
  }

  PetscErrorCode addThermalFluxElement(
    const string problem_name,const string field_name,const string mesh_nodals_positions = "MESH_NODE_POSITIONS") {
    PetscFunctionBegin;

    PetscErrorCode ierr;
    ErrorCode rval;

    ierr = mField.add_finite_element("THERMAL_FLUX_FE",MF_ZERO); CHKERRQ(ierr);
    ierr = mField.modify_finite_element_add_field_row("THERMAL_FLUX_FE",field_name); CHKERRQ(ierr);
    ierr = mField.modify_finite_element_add_field_col("THERMAL_FLUX_FE",field_name); CHKERRQ(ierr);
    ierr = mField.modify_finite_element_add_field_data("THERMAL_FLUX_FE",field_name); CHKERRQ(ierr);
    if(mField.check_field(mesh_nodals_positions)) {
      ierr = mField.modify_finite_element_add_field_data("THERMAL_FLUX_FE",mesh_nodals_positions); CHKERRQ(ierr);
    }
    ierr = mField.modify_problem_add_finite_element(problem_name,"THERMAL_FLUX_FE"); CHKERRQ(ierr);

    for(_IT_CUBITMESHSETS_BY_BCDATA_TYPE_FOR_LOOP_(mField,SideSet|HeatfluxSet,it)) {
      ierr = it->get_cubit_bc_data_structure(setOfFluxes[it->get_msId()].dAta); CHKERRQ(ierr);
      rval = mField.get_moab().get_entities_by_type(it->meshset,MBTRI,setOfFluxes[it->get_msId()].tRis,true); CHKERR_PETSC(rval);
      ierr = mField.add_ents_to_finite_element_by_TRIs(setOfFluxes[it->get_msId()].tRis,"THERMAL_FLUX_FE"); CHKERRQ(ierr);
    }

    PetscFunctionReturn(0);
  }

  PetscErrorCode setThermalFiniteElementRhsOperators(string field_name,Vec &F) {
    PetscFunctionBegin;
    map<int,BlockData>::iterator sit = setOfBlocks.begin();
    for(;sit!=setOfBlocks.end();sit++) {
      //add finite element
      feRhs.get_op_to_do_Rhs().push_back(new OpGetGradAtGaussPts(field_name,commonData));
      feRhs.get_op_to_do_Rhs().push_back(new OpThermalRhs(field_name,F,sit->second,commonData));
    }
    PetscFunctionReturn(0);
  }

  PetscErrorCode setThermalFiniteElementLhsOperators(string field_name,Mat *A) {
    PetscFunctionBegin;
    map<int,BlockData>::iterator sit = setOfBlocks.begin();
    for(;sit!=setOfBlocks.end();sit++) {
      //add finite elemen
      feLhs.get_op_to_do_Lhs().push_back(new OpThermalLhs(field_name,A,sit->second,commonData));
    }
    PetscFunctionReturn(0);
  }

  PetscErrorCode setThermalFluxFiniteElementLhsOperators(string field_name,Vec &F,const string mesh_nodals_positions = "MESH_NODE_POSITIONS") {
    PetscFunctionBegin;
    bool ho_geometry = false;
    if(mField.check_field(mesh_nodals_positions)) {
      ho_geometry = true;
    }
    map<int,FluxData>::iterator sit = setOfFluxes.begin();
    for(;sit!=setOfFluxes.end();sit++) {
      //add finite element
      feFlux.get_op_to_do_Rhs().push_back(new OpHeatFlux(field_name,F,sit->second,ho_geometry));
    }
    PetscFunctionReturn(0);
  }

  PetscErrorCode setTimeSteppingProblem(TsCtx &ts_ctx,string field_name,string rate_name,const string mesh_nodals_positions = "MESH_NODE_POSITIONS") {
    PetscFunctionBegin;

    {
      map<int,BlockData>::iterator sit = setOfBlocks.begin();
      for(;sit!=setOfBlocks.end();sit++) {
	//add finite element
	feLhs.get_op_to_do_Lhs().push_back(new OpThermalLhs(field_name,sit->second,commonData));
	feLhs.get_op_to_do_Lhs().push_back(new OpHeatCapacityLsh(field_name,sit->second,commonData));
	feRhs.get_op_to_do_Rhs().push_back(new OpGetRateAtGaussPts(rate_name,commonData));
	feRhs.get_op_to_do_Rhs().push_back(new OpGetGradAtGaussPts(field_name,commonData));
	feRhs.get_op_to_do_Rhs().push_back(new OpThermalRhs(field_name,sit->second,commonData));
	feRhs.get_op_to_do_Rhs().push_back(new OpHeatCapacityRhs(field_name,sit->second,commonData));
      }
    }
    {
      bool ho_geometry = false;
      if(mField.check_field(mesh_nodals_positions)) {
	ho_geometry = true;
      }
      map<int,FluxData>::iterator sit = setOfFluxes.begin();
      for(;sit!=setOfFluxes.end();sit++) {
	//add finite element
	feFlux.get_op_to_do_Rhs().push_back(new OpHeatFlux(field_name,sit->second,ho_geometry));
      }
    }

    //rhs
    TsCtx::loops_to_do_type& loops_to_do_Rhs = ts_ctx.get_loops_to_do_IFunction();
    loops_to_do_Rhs.push_back(TsCtx::loop_pair_type("THERMAL_FE",&feRhs));
    loops_to_do_Rhs.push_back(TsCtx::loop_pair_type("THERMAL_FLUX_FE",&feFlux));

    //lhs
    TsCtx::loops_to_do_type& loops_to_do_Mat = ts_ctx.get_loops_to_do_IJacobian();
    loops_to_do_Mat.push_back(TsCtx::loop_pair_type("THERMAL_FE",&feLhs));

    //monitor
    //TsCtx::loops_to_do_type& loops_to_do_Monitor = ts_ctx.get_loops_to_do_Monitor();

    PetscFunctionReturn(0);
  }

};

}

#endif //__BODY_FORCE_HPP

