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

  struct virtual PetscErrorCode getFlux(const double x,const double y,const double z,double &flux) = 0; 

  /// \brief  definition of volume element
  struct MyVolumeFE: public TetElementForcesAndSourcesCore {
    MyVolumeFE(FieldInterface &_mField): TetElementForcesAndSourcesCore(_mField) {}
    int getRule(int order) { return order-1; };
  };

  MyVolumeFE fe;   

  /** \brief define surface element
    *
    */
  struct MyTriFE: public TriElementForcesAndSurcesCore {
    MyTriFE(FieldInterface &_mField): TriElementForcesAndSurcesCore(_mField) {}
    int getRule(int order) { return ceil(order/2); };
  };


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

      try {

	if(row_data.getFieldData().size()==0) PetscFunctionReturn(0);
	if(col_data.getFieldData().size()==0) PetscFunctionReturn(0);

	int nb_row = row_data.getFieldData().size();
	int nb_col = col_data.getFieldData().size();
	NN.resize(nb_row,nb_col);
	bzero(&NN(0,0),NN.data().size()*sizeof(FieldData));

	int nb_gauss_pts = data_row.getHdivN().size1();
	int gg = 0;
	for(;gg<nb_gauss_pts;gg++) {

	  double w = getGaussPts(gg,3)*getVolume();
	  if(getHoGaussPtsDetJac().size()>0) {
	    w *= getHoGaussPtsDetJac();
	  }
      
	  for(int rr = 0;rr<nb_row;rr++) {
	    for(int cc = 0;cc<nb_col;cc++) {
	      NN(rr,cc) += w*(
		data_row.getHdivN(gg)(rr,0)*data_col.getHdivN(gg)(rr,0)+
		data_row.getHdivN(gg)(rr,1)*data_col.getHdivN(gg)(rr,1)+
		data_row.getHdivN(gg)(rr,2)*data_col.getHdivN(gg)(rr,2) );
	    }
	  }

	}

	ierr = MatSetValues(
	  Aij,
	  nb_row,&row_data.getIndices()[0],
	  nb_col,&col_data.getIndices()[0],
	  &NN(0,0),ADD_VALUES); CHKERRQ(ierr);
	if(row_side != col_side || row_type != col_type) {
	  transK.resize(nb_col,nb_row);
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
  struct OpLhsDivTauDotU_HdivL2: public TetElementForcesAndSourcesCore::UserDataOperator {

    UltraWeakTransportElement &cTx;
    Mat Aij;

    OpLhsDivTauDotU_HdivL2(
      UltraWeakTransportElement &ctx,
      const string field_name_row,string field_name_col,Mat _Aij):
      TetElementForcesAndSourcesCore::UserDataOperator(field_name_row,field_name_col),
      cTx(ctx),Aij(_Aij) {
  
      //this operator is not symmetric settig this varible makes element
      //operator to loop over element entities (subsimplicies) without
      //assumption that off-diagonal matrices are symmetric.
      symm = false;

    }

    ublas::matrix<FieldData> NN;
    ublas::vector<FieldData> div_vec;

    PetscErrorCode doWork(
      int row_side,int col_side,
      EntityType row_type,EntityType col_type,
      DataForcesAndSurcesCore::EntData &row_data,
      DataForcesAndSurcesCore::EntData &col_data) {
      PetscFunctionBegin;

      try {

	if(row_data.getFieldData().size()==0) PetscFunctionReturn(0);
	if(col_data.getFieldData().size()==0) PetscFunctionReturn(0);

	int nb_row = row_data.getFieldData().size();
	int nb_col = col_data.getFieldData().size();
	NN.resize(nb_row,nb_col);
	bzero(&NN(0,0),NN.data().size()*sizeof(FieldData));

	div_vec.resize(nb_row,0);

	int nb_gauss_pts = data_row.getHdivN().size1();
	int gg = 0;
	for(;gg<nb_gauss_pts;gg++) {

	  double w = getGaussPts(gg,3)*getVolume();
	  if(getHoGaussPtsDetJac().size()>0) {
	    w *= getHoGaussPtsDetJac();
	  }

	  ierr = getDivergenceMatrixOperato_Hdiv(
	    row_side,row_type,row_data,gg,div_vec); CHKERRQ(ierr);
      
	  for(int rr = 0;rr<nb_row;rr++) {
	    for(int cc = 0;cc<nb_col;cc++) {
	      NN(rr,cc) += w*(div_vec[rr]*getN(gg)[cc]);
	    }
	  }

	}

	ierr = MatSetValues(
	  Aij,
	  nb_row,&row_data.getIndices()[0],
	  nb_col,&col_data.getIndices()[0],
	  &NN(0,0),ADD_VALUES); CHKERRQ(ierr);

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
  struct OpLhsDivTauDotU_L2Hdiv: public TetElementForcesAndSourcesCore::UserDataOperator {

    UltraWeakTransportElement &cTx;
    Mat Aij;

    OpLhsDivTauDotU_L2Hdiv(
      UltraWeakTransportElement &ctx,
      const string field_name_row,string field_name_col,Mat _Aij):
      TetElementForcesAndSourcesCore::UserDataOperator(field_name_row,field_name_col),
      cTx(ctx),Aij(_Aij) {
  
      //this operator is not symmetric settig this varible makes element
      //operator to loop over element entities (subsimplicies) without
      //assumption that off-diagonal matrices are symmetric.
      symm = false;

    }

    ublas::matrix<FieldData> NN;
    ublas::vector<FieldData> div_vec;

    PetscErrorCode doWork(
      int row_side,int col_side,
      EntityType row_type,EntityType col_type,
      DataForcesAndSurcesCore::EntData &row_data,
      DataForcesAndSurcesCore::EntData &col_data) {
      PetscFunctionBegin;

      try {

	if(row_data.getFieldData().size()==0) PetscFunctionReturn(0);
	if(col_data.getFieldData().size()==0) PetscFunctionReturn(0);

	int nb_row = row_data.getFieldData().size();
	int nb_col = col_data.getFieldData().size();
	NN.resize(nb_row,nb_col);
	bzero(&NN(0,0),NN.data().size()*sizeof(FieldData));

	div_vec.resize(nb_col,0);

	int nb_gauss_pts = data_row.getHdivN().size1();
	int gg = 0;
	for(;gg<nb_gauss_pts;gg++) {

	  double w = getGaussPts(gg,3)*getVolume();
	  if(getHoGaussPtsDetJac().size()>0) {
	    w *= getHoGaussPtsDetJac();
	  }

	  ierr = getDivergenceMatrixOperato_Hdiv(
	    col_side,col_type,col_data,gg,div_vec); CHKERRQ(ierr);
      
	  for(int rr = 0;rr<nb_row;rr++) {
	    for(int cc = 0;cc<nb_col;cc++) {
	      NN(rr,cc) += w*(getN(gg)[rr]*div_vec[cc]);
	    }
	  }

	}

	ierr = MatSetValues(
	  Aij,
	  nb_row,&row_data.getIndices()[0],
	  nb_col,&col_data.getIndices()[0],
	  &NN(0,0),ADD_VALUES); CHKERRQ(ierr);

      } catch (const std::exception& ex) {
	ostringstream ss;
	ss << "throw in method: " << ex.what() << endl;
	SETERRQ(PETSC_COMM_SELF,1,ss.str().c_str());
      }

      PetscFunctionReturn(0);
    }

  };

  struct OpL2Source public TetElementForcesAndSourcesCore::UserDataOperator {

    UltraWeakTransportElement &cTx;
    Vec F;

    OpL2Source(
      UltraWeakTransportElement &ctx,
      const string field_name,Mat _F):
      TetElementForcesAndSourcesCore::UserDataOperator(field_name),
      cTx(ctx),F(_F) {}
      
    ublas::vector<FieldData> Nf;
    PetscErrorCode doWork(
      int side,EntityType type,DataForcesAndSurcesCore::EntData &data) {
      PetscFunctionBegin;

      try {

	if(data.getFieldData().size()==0) PetscFunctionReturn(0);

	int nb_row = data.getFieldData().size();
	Nf.resize(nb_row);
	bzero(&Nf(0,0),Nf.data().size()*sizeof(FieldData));

	int nb_gauss_pts = data_row.getHdivN().size1();
	int gg = 0;
	for(;gg<nb_gauss_pts;gg++) {

	  double w = getGaussPts(gg,3)*getVolume();
	  if(getHoGaussPtsDetJac().size()>0) {
	    w *= getHoGaussPtsDetJac();
	  }
      
	  for(int rr = 0;rr<nb_row;rr++) {
	    const double x = getCoordsAtGaussPts()(gg,0);
	    const double y = getCoordsAtGaussPts()(gg,1);
	    const double z = getCoordsAtGaussPts()(gg,2);
	    double flux;
	    ierr = cTx.getFlux(x,y,z,flux); CHKERRQ(ierr);
	    Nf[rr] += w*data.getN(gg)[rr]*flux;
	  }
	
	}

	ierr = MatSetValues(
	  Aij,
	  nb_row,&row_data.getIndices()[0],
	  nb_col,&col_data.getIndices()[0],
	  &NN(0,0),ADD_VALUES); CHKERRQ(ierr);

      } catch (const std::exception& ex) {
	ostringstream ss;
	ss << "throw in method: " << ex.what() << endl;
	SETERRQ(PETSC_COMM_SELF,1,ss.str().c_str());
      }


      PetscFunctionReturn(0);
    }

  };



}


}

#endif //__ULTRAWEAK_TARNSPORT_ELEMENT_HPP

/***************************************************************************//**
 * \defgroup mofem_ultra_weak_transport_elem Ultra weak transport element
 * \ingroup mofem_forces_and_sources 
 ******************************************************************************/



