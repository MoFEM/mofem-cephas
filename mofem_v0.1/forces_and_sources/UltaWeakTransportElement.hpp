/** \file UltraWeakTransportElement.hpp
 * \brief Ultra weak implementation of transport element
 *
 */

/* Copyright (C) 2013, Lukasz Kaczmarczyk (likask AT wp.pl)
 * --------------------------------------------------------------
 *
 * This file is part of MoFEM.
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

#ifndef __ULTRAWEAK_TARNSPORT_ELEMENT_HPP
#define __ULTRAWEAK_TARNSPORT_ELEMENT_HPP

namespace MoFEM {

struct UltraWeakTransportElement {

  FieldInterface &mField;

  /// \brief  definition of volume element
  struct MyVolumeFE: public TetElementForcesAndSourcesCore {
    MyVolumeFE(FieldInterface &m_field): TetElementForcesAndSourcesCore(m_field) {}
    int getRule(int order) { return order+1; };
  };

  MyVolumeFE feVolFluxFlux;   
  MyVolumeFE feVolValueFlux;  

  /** \brief define surface element
    *
    */
  struct MyTriFE: public TriElementForcesAndSurcesCore {
    MyTriFE(FieldInterface &m_field): TriElementForcesAndSurcesCore(m_field) {}
    int getRule(int order) { return order+1; };
  };

  UltraWeakTransportElement(FieldInterface &m_field): mField(m_field),
    feVolFluxFlux(m_field),feVolValueFlux(m_field) {};

  virtual PetscErrorCode getFlux(const double x,const double y,const double z,double &flux) {
    PetscFunctionBegin;
    flux  = 0;
    PetscFunctionReturn(0);
  } 

  /** \brief data for calulation het conductivity and heat capacity elements
    * \infroup mofem_thermal_elem 
    */
  struct BlockData {
    double cOnductivity;
    double cApacity;
    Range tEts; ///< constatins elements in block set
  }; 
  map<int,BlockData> setOfBlocks; ///< maps block set id with appropiate BlockData

  /// \brief add finite elements
  PetscErrorCode addFiniteElements(
    const string &fluxes_name,const string &values_name,const string mesh_nodals_positions = "MESH_NODE_POSITIONS") {
    PetscFunctionBegin;

    PetscErrorCode ierr;
    ErrorCode rval;

    ierr = mField.add_finite_element("ULTRAWEAK_FLUXFLUX",MF_ZERO); CHKERRQ(ierr);
    ierr = mField.modify_finite_element_add_field_row("ULTRAWEAK_FLUXFLUX",fluxes_name); CHKERRQ(ierr);
    ierr = mField.modify_finite_element_add_field_col("ULTRAWEAK_FLUXFLUX",fluxes_name); CHKERRQ(ierr);
    ierr = mField.modify_finite_element_add_field_data("ULTRAWEAK_FLUXFLUX",fluxes_name); CHKERRQ(ierr);
    if(mField.check_field(mesh_nodals_positions)) {
      ierr = mField.modify_finite_element_add_field_data("ULTRAWEAK_FLUXFLUX",fluxes_name); CHKERRQ(ierr);
    }

    ierr = mField.add_finite_element("ULTRAWEAK_FLUXVALUE",MF_ZERO); CHKERRQ(ierr);
    ierr = mField.modify_finite_element_add_field_row("ULTRAWEAK_FLUXVALUE",fluxes_name); CHKERRQ(ierr);
    ierr = mField.modify_finite_element_add_field_col("ULTRAWEAK_FLUXVALUE",values_name); CHKERRQ(ierr);
    ierr = mField.modify_finite_element_add_field_data("ULTRAWEAK_FLUXVALUE",fluxes_name); CHKERRQ(ierr);
    ierr = mField.modify_finite_element_add_field_data("ULTRAWEAK_FLUXVALUE",values_name); CHKERRQ(ierr);
    if(mField.check_field(mesh_nodals_positions)) {
      ierr = mField.modify_finite_element_add_field_data("ULTRAWEAK_FLUXVALUE",fluxes_name); CHKERRQ(ierr);
    }

    ierr = mField.add_finite_element("ULTRAWEAK_VALUEFLUX",MF_ZERO); CHKERRQ(ierr);
    ierr = mField.modify_finite_element_add_field_row("ULTRAWEAK_VALUEFLUX",values_name); CHKERRQ(ierr);
    ierr = mField.modify_finite_element_add_field_col("ULTRAWEAK_VALUEFLUX",fluxes_name); CHKERRQ(ierr);
    ierr = mField.modify_finite_element_add_field_data("ULTRAWEAK_VALUEFLUX",fluxes_name); CHKERRQ(ierr);
    ierr = mField.modify_finite_element_add_field_data("ULTRAWEAK_VALUEFLUX",values_name); CHKERRQ(ierr);
    if(mField.check_field(mesh_nodals_positions)) {
      ierr = mField.modify_finite_element_add_field_data("ULTRAWEAK_VALUEFLUX",fluxes_name); CHKERRQ(ierr);
    }

    for(_IT_CUBITMESHSETS_BY_BCDATA_TYPE_FOR_LOOP_(mField,BLOCKSET|MAT_THERMALSET,it)) {

      Mat_Thermal temp_data;
      ierr = it->get_attribute_data_structure(temp_data); CHKERRQ(ierr);  
      setOfBlocks[it->get_msId()].cOnductivity = temp_data.data.Conductivity;
      setOfBlocks[it->get_msId()].cApacity = temp_data.data.HeatCapacity;
      rval = mField.get_moab().get_entities_by_type(it->meshset,MBTET,setOfBlocks[it->get_msId()].tEts,true); CHKERR_PETSC(rval);
      ierr = mField.add_ents_to_finite_element_by_TETs(setOfBlocks[it->get_msId()].tEts,"ULTRAWEAK_FLUXFLUX"); CHKERRQ(ierr);
      ierr = mField.add_ents_to_finite_element_by_TETs(setOfBlocks[it->get_msId()].tEts,"ULTRAWEAK_FLUXVALUE"); CHKERRQ(ierr);
      ierr = mField.add_ents_to_finite_element_by_TETs(setOfBlocks[it->get_msId()].tEts,"ULTRAWEAK_VALUEFLUX"); CHKERRQ(ierr);

    }

    PetscFunctionReturn(0);
  }

  /** \brief tau,sigma in Hdiv, calulates Aij = Asemble int sigma_dot_tau dTet 
    */
  struct OpLhsTauDotSigma_HdivHdiv: public TetElementForcesAndSourcesCore::UserDataOperator {

    UltraWeakTransportElement &cTx;
    Mat Aij;

    OpLhsTauDotSigma_HdivHdiv(
      UltraWeakTransportElement &ctx,
      const string field_name,Mat _Aij):
      TetElementForcesAndSourcesCore::UserDataOperator(field_name),
      cTx(ctx),Aij(_Aij) {}

    ublas::matrix<FieldData> NN,transNN;
    PetscErrorCode doWork(
      int row_side,int col_side,
      EntityType row_type,EntityType col_type,
      DataForcesAndSurcesCore::EntData &row_data,
      DataForcesAndSurcesCore::EntData &col_data) {
      PetscFunctionBegin;

      PetscErrorCode ierr;

      try {

	if(row_data.getFieldData().size()==0) PetscFunctionReturn(0);
	if(col_data.getFieldData().size()==0) PetscFunctionReturn(0);

	int nb_row = row_data.getFieldData().size();
	int nb_col = col_data.getFieldData().size();
	NN.resize(nb_row,nb_col);
	bzero(&NN(0,0),NN.data().size()*sizeof(FieldData));

	int nb_gauss_pts = row_data.getHdivN().size1();
	int gg = 0;
	for(;gg<nb_gauss_pts;gg++) {

	  double w = getGaussPts()(3,gg)*getVolume();
	  if(getHoGaussPtsDetJac().size()>0) {
	    w *= getHoGaussPtsDetJac()(gg);
	  }
      
	  //FIXME this multiplication should be done in blas or ublas
	  /*for(int rr = 0;rr<nb_row;rr++) {
	    for(int cc = 0;cc<nb_col;cc++) {
	      NN(rr,cc) += w*(
		row_data.getHdivN(gg)(rr,0)*col_data.getHdivN(gg)(cc,0)+
		row_data.getHdivN(gg)(rr,1)*col_data.getHdivN(gg)(cc,1)+
		row_data.getHdivN(gg)(rr,2)*col_data.getHdivN(gg)(cc,2) );
	    }
	  }*/
	  noalias(NN) += w*prod(row_data.getHdivN(gg),trans(col_data.getHdivN(gg))); 


	}

	ierr = MatSetValues(
	  Aij,
	  nb_row,&row_data.getIndices()[0],
	  nb_col,&col_data.getIndices()[0],
	  &NN(0,0),ADD_VALUES); CHKERRQ(ierr);

	if(row_side != col_side || row_type != col_type) {
	  transNN.resize(nb_col,nb_row);
	  noalias(transNN) = trans(NN);
	  ierr = MatSetValues(
	    Aij,
	    nb_col,&col_data.getIndices()[0],
	    nb_row,&row_data.getIndices()[0],
	    &transNN(0,0),ADD_VALUES); CHKERRQ(ierr);
	}


      } catch (const std::exception& ex) {
	ostringstream ss;
	ss << "throw in method: " << ex.what() << endl;
	SETERRQ(PETSC_COMM_SELF,1,ss.str().c_str());
      }

      PetscFunctionReturn(0);
    }

  };


  /** \brief u in L2 and tau in Hdiv, calulates Aij = Asemble int u * div(tau) dTet 
    */
  struct OpLhsVDotDivSigma_L2Hdiv: public TetElementForcesAndSourcesCore::UserDataOperator {

    UltraWeakTransportElement &cTx;
    Mat Aij;

    OpLhsVDotDivSigma_L2Hdiv(
      UltraWeakTransportElement &ctx,
      const string field_name_row,string field_name_col,Mat _Aij):
      TetElementForcesAndSourcesCore::UserDataOperator(field_name_row,field_name_col),
      cTx(ctx),Aij(_Aij) {
  
      //this operator is not symmetric settig this varible makes element
      //operator to loop over element entities (subsimplicies) without
      //assumption that off-diagonal matrices are symmetric.
      symm = false;

    }

    ublas::matrix<FieldData> NN,transNN;
    ublas::vector<FieldData> div_vec;

    PetscErrorCode doWork(
      int row_side,int col_side,
      EntityType row_type,EntityType col_type,
      DataForcesAndSurcesCore::EntData &row_data,
      DataForcesAndSurcesCore::EntData &col_data) {
      PetscFunctionBegin;

      PetscErrorCode ierr;

      try {

	if(row_data.getFieldData().size()==0) PetscFunctionReturn(0);
	if(col_data.getFieldData().size()==0) PetscFunctionReturn(0);

	int nb_row = row_data.getFieldData().size();
	int nb_col = col_data.getFieldData().size();
	NN.resize(nb_row,nb_col);
	bzero(&NN(0,0),NN.data().size()*sizeof(FieldData));

	div_vec.resize(nb_col,0);

	int nb_gauss_pts = row_data.getHdivN().size1();
	int gg = 0;
	for(;gg<nb_gauss_pts;gg++) {

	  double w = getGaussPts()(3,gg)*getVolume();
	  if(getHoGaussPtsDetJac().size()>0) {
	    w *= getHoGaussPtsDetJac()(gg);
	  }

	  ierr = getDivergenceMatrixOperato_Hdiv(
	    col_side,col_type,col_data,gg,div_vec); CHKERRQ(ierr);
      
	  //FIXME this multiplication should be done in blas or ublas
	  for(int rr = 0;rr<nb_row;rr++) {
	    for(int cc = 0;cc<nb_col;cc++) {
	      NN(rr,cc) += w*(row_data.getN(gg)[rr]*div_vec[cc]);
	    }
	  }

	}

	/*cerr << row_side << " " << col_side << endl;
	cerr << row_type << " " << col_type << endl;
	cerr << row_data.getIndices() << endl;
	cerr << col_data.getIndices() << endl;
	cerr << row_data.getFieldData() << endl;
	cerr << col_data.getFieldData() << endl;
	cerr << NN << endl << endl;*/

	ierr = MatSetValues(
	  Aij,
	  nb_row,&row_data.getIndices()[0],
	  nb_col,&col_data.getIndices()[0],
	  &NN(0,0),ADD_VALUES); CHKERRQ(ierr);

	transNN.resize(nb_col,nb_row);
	ublas::noalias(transNN) = trans(NN);
	ierr = MatSetValues(
	  Aij,
	  nb_col,&col_data.getIndices()[0],
	  nb_row,&row_data.getIndices()[0],
	  &transNN(0,0),ADD_VALUES); CHKERRQ(ierr);


      } catch (const std::exception& ex) {
	ostringstream ss;
	ss << "throw in method: " << ex.what() << endl;
	SETERRQ(PETSC_COMM_SELF,1,ss.str().c_str());
      }

      PetscFunctionReturn(0);
    }

  };

  /** \brief calulate source therms
    */
  struct OpL2Source: public TetElementForcesAndSourcesCore::UserDataOperator {

    UltraWeakTransportElement &cTx;
    Vec F;

    OpL2Source(
      UltraWeakTransportElement &ctx,
      const string field_name,Vec _F):
      TetElementForcesAndSourcesCore::UserDataOperator(field_name),
      cTx(ctx),F(_F) {}
      
    ublas::vector<FieldData> Nf;
    PetscErrorCode doWork(
      int side,EntityType type,DataForcesAndSurcesCore::EntData &data) {
      PetscFunctionBegin;

      PetscErrorCode ierr;

      try {

	if(data.getFieldData().size()==0) PetscFunctionReturn(0);

	int nb_row = data.getFieldData().size();
	Nf.resize(nb_row);
	bzero(&Nf[0],Nf.data().size()*sizeof(FieldData));

	int nb_gauss_pts = data.getHdivN().size1();
	int gg = 0;
	for(;gg<nb_gauss_pts;gg++) {

	  double w = getGaussPts()(3,gg)*getVolume();
	  if(getHoGaussPtsDetJac().size()>0) {
	    w *= getHoGaussPtsDetJac()(gg);
	  }
      

	  const double x = getCoordsAtGaussPts()(gg,0);
	  const double y = getCoordsAtGaussPts()(gg,1);
	  const double z = getCoordsAtGaussPts()(gg,2);
	  double flux;
	  ierr = cTx.getFlux(x,y,z,flux); CHKERRQ(ierr);

	  for(int rr = 0;rr<nb_row;rr++) {
	    Nf[rr] += w*data.getN(gg)[rr]*flux;
	  }
	
	}

	ierr = VecSetValues(
	  F,nb_row,&data.getIndices()[0],
	  &Nf[0],ADD_VALUES); CHKERRQ(ierr);

      } catch (const std::exception& ex) {
	ostringstream ss;
	ss << "throw in method: " << ex.what() << endl;
	SETERRQ(PETSC_COMM_SELF,1,ss.str().c_str());
      }


      PetscFunctionReturn(0);
    }

  };

  ublas::vector<ublas::vector<FieldData> > valuesGradientAtGaussPts;
  ublas::vector<FieldData> divergenceAtGaussPts;
  ublas::vector<ublas::vector<FieldData> > fluxesAtGaussPts;

  struct OpValuesGradientAtGaussPts: public TetElementForcesAndSourcesCore::UserDataOperator {

    UltraWeakTransportElement &cTx;

    OpValuesGradientAtGaussPts(
      UltraWeakTransportElement &ctx,
      const string field_name):
      TetElementForcesAndSourcesCore::UserDataOperator(field_name),
      cTx(ctx) {}

    PetscErrorCode doWork(
      int side,EntityType type,DataForcesAndSurcesCore::EntData &data) {
      PetscFunctionBegin;

      try {

      if(data.getFieldData().size() == 0)  PetscFunctionReturn(0);
  
      int nb_gauss_pts = data.getDiffN().size1();

      cTx.valuesGradientAtGaussPts.resize(nb_gauss_pts);

      for(int gg = 0;gg<nb_gauss_pts;gg++) {

	cTx.valuesGradientAtGaussPts[gg].resize(3);
	noalias(cTx.valuesGradientAtGaussPts[gg]) = prod( trans(data.getDiffN(gg)), data.getFieldData() );
      }
      
      } catch (const std::exception& ex) {
	ostringstream ss;
	ss << "throw in method: " << ex.what() << endl;
	SETERRQ(PETSC_COMM_SELF,1,ss.str().c_str());
      }

      PetscFunctionReturn(0);
    }

  };

  struct OpFluxDivergenceAtGaussPts: public TetElementForcesAndSourcesCore::UserDataOperator {

    UltraWeakTransportElement &cTx;

    OpFluxDivergenceAtGaussPts(
      UltraWeakTransportElement &ctx,
      const string field_name):
      TetElementForcesAndSourcesCore::UserDataOperator(field_name),
      cTx(ctx) {}

    ublas::vector<FieldData> div_vec;
    PetscErrorCode doWork(
      int side,EntityType type,DataForcesAndSurcesCore::EntData &data) {
      PetscFunctionBegin;
      PetscErrorCode ierr;

      try {

      if(data.getFieldData().size() == 0)  PetscFunctionReturn(0);

      int nb_gauss_pts = data.getDiffN().size1();
      int nb_dofs = data.getFieldData().size();

      cTx.fluxesAtGaussPts.resize(nb_gauss_pts);
      cTx.divergenceAtGaussPts.resize(nb_gauss_pts);
      if(type == MBTRI && side == 0) {
	for(int gg = 0;gg<nb_gauss_pts;gg++) {
	  cTx.fluxesAtGaussPts[gg].resize(3);
	  bzero(&(cTx.fluxesAtGaussPts[gg][0]),3*sizeof(FieldData));
	}
	bzero(&cTx.divergenceAtGaussPts[0],nb_gauss_pts*sizeof(FieldData));
      }

      div_vec.resize(nb_dofs);

      for(int gg = 0;gg<nb_gauss_pts;gg++) {

	ierr = getDivergenceMatrixOperato_Hdiv(
	  side,type,data,gg,div_vec); CHKERRQ(ierr);

	cTx.divergenceAtGaussPts[gg] += inner_prod(div_vec,data.getFieldData());
	noalias(cTx.fluxesAtGaussPts[gg]) += prod( trans(data.getHdivN(gg)), data.getFieldData());

      }

      } catch (const std::exception& ex) {
	ostringstream ss;
	ss << "throw in method: " << ex.what() << endl;
	SETERRQ(PETSC_COMM_SELF,1,ss.str().c_str());
      }

      PetscFunctionReturn(0);
    }

  };


  /** \brief calulate error evaluator
    */
  struct OpError_L2Norm: public TetElementForcesAndSourcesCore::UserDataOperator {

    UltraWeakTransportElement &cTx;

    OpError_L2Norm(
      UltraWeakTransportElement &ctx,
      const string field_name):
      TetElementForcesAndSourcesCore::UserDataOperator(field_name),
      cTx(ctx) {}

    ublas::vector<FieldData> deltaFlux;
    PetscErrorCode doWork(
      int side,EntityType type,DataForcesAndSurcesCore::EntData &data) {
      PetscFunctionBegin;

      PetscErrorCode ierr;
      ErrorCode rval;

      try {

      if(data.getFieldData().size() == 0)  PetscFunctionReturn(0);
      int nb_gauss_pts = data.getN().size1();
      EntityHandle fe_ent = getMoFEMFEPtr()->get_ent();
    
      double def_VAL = 0;
      Tag th_error_div_sigma;
      rval = cTx.mField.get_moab().tag_get_handle("ERRORL2_DIVSIGMA_F",1,MB_TYPE_DOUBLE,th_error_div_sigma,MB_TAG_CREAT|MB_TAG_SPARSE,&def_VAL); CHKERR_PETSC(rval);
      Tag th_error_flux;
      rval = cTx.mField.get_moab().tag_get_handle("ERRORL2_FLUX",1,MB_TYPE_DOUBLE,th_error_flux,MB_TAG_CREAT|MB_TAG_SPARSE,&def_VAL); CHKERR_PETSC(rval);
      if(type == MBTRI && side == 0) {
	cTx.mField.get_moab().tag_set_data(th_error_div_sigma,&fe_ent,1,&def_VAL);
      }
      
      double* error_div_ptr;
      rval = cTx.mField.get_moab().tag_get_by_ptr(th_error_div_sigma,&fe_ent,1,(const void**)&error_div_ptr); CHKERR_PETSC(rval);
      double* error_flux_ptr;
      rval = cTx.mField.get_moab().tag_get_by_ptr(th_error_flux,&fe_ent,1,(const void**)&error_flux_ptr); CHKERR_PETSC(rval);

      deltaFlux.resize(3);
      for(int gg = 0;gg<nb_gauss_pts;gg++) {
	double w = getGaussPts()(3,gg)*getVolume();
	if(getHoGaussPtsDetJac().size()>0) {
	  w *= getHoGaussPtsDetJac()(gg);
	}

	const double x = getCoordsAtGaussPts()(gg,0);
	const double y = getCoordsAtGaussPts()(gg,1);
	const double z = getCoordsAtGaussPts()(gg,2);
	double flux;
	ierr = cTx.getFlux(x,y,z,flux); CHKERRQ(ierr);

	*error_div_ptr += w*( pow(cTx.divergenceAtGaussPts[gg] - flux,2) );
	noalias(deltaFlux) = cTx.fluxesAtGaussPts[gg]-cTx.valuesGradientAtGaussPts[gg];
	*error_flux_ptr += w*( inner_prod(deltaFlux,deltaFlux) );

      }

      } catch (const std::exception& ex) {
	ostringstream ss;
	ss << "throw in method: " << ex.what() << endl;
	SETERRQ(PETSC_COMM_SELF,1,ss.str().c_str());
      }

      PetscFunctionReturn(0);
    }

  };

};


}

#endif //__ULTRAWEAK_TARNSPORT_ELEMENT_HPP

/***************************************************************************//**
 * \defgroup mofem_ultra_weak_transport_elem Ultra weak transport element
 * \ingroup mofem_forces_and_sources 
 ******************************************************************************/



