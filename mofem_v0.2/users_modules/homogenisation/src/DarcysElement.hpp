/** \file DarcysElement.hpp
 * \brief Operators and data structures for thermal analys
 *
 * Implementation of Darceys element for unsteady and steady case.
 *
 */

/* Copyright (C) 2013, Zahur Ullah (Zahur.Ullah@glasgow.ac.uk)
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

#ifndef __THERMAL_ELEMENT_HPP
#define __THERMAL_ELEMENT_HPP

//#include<moab/Skinner.hpp>

namespace MoFEM {

/** \brief struture grouping operators and data used for thermal problems
  * \ingroup mofem_thermal_elem 
  *
  * In order to assemble matrices and right hand vectors, the loops over
  * elements, enetities over that elememnts and finally loop over intergration
  * points are executed.
  *
  * Following implementation separte those three cegories of loops and to eeach
  * loop attach operator.
  *
  */
struct DarcysElement {

  /// \brief  definition of volume element
  struct MyVolumeFE: public TetElementForcesAndSourcesCore {
    MyVolumeFE(FieldInterface &_mField): TetElementForcesAndSourcesCore(_mField) {}
    
    /** \brief it is used to calculate nb. of Gauss integartion points
     *
     * for more details pleas look 
     *   Reference:
     *
     * Albert Nijenhuis, Herbert Wilf,
     * Combinatorial Algorithms for Computers and Calculators,
     * Second Edition,
     * Academic Press, 1978,
     * ISBN: 0-12-519260-6,
     * LC: QA164.N54.
     *
     * More details about algorithm 
     * http://people.sc.fsu.edu/~jburkardt/cpp_src/gm_rule/gm_rule.html
    **/
    int getRule(int order) { return order-1; };
  };
  
  MyVolumeFE feRhs; ///< cauclate right hand side for tetrahedral elements
  MyVolumeFE& getLoopFeRhs() { return feRhs; } ///< get rhs volume element 
  MyVolumeFE feLhs; //< calculate left hand side for tetrahedral elements
  MyVolumeFE& getLoopFeLhs() { return feLhs; } ///< get lhs volume element

  FieldInterface &mField;
  DarcysElement(
    FieldInterface &m_field):
    feRhs(m_field),feLhs(m_field),mField(m_field) {}

  /** \brief data for calulation het conductivity and heat capacity elements
    * \infroup mofem_thermal_elem 
    */
  struct BlockData {
    double vIscosity;
    double pErmeability;
    Range tEts; ///< constatins elements in block set
  }; 
  map<int,BlockData> setOfBlocks; ///< maps block set id with appropiate BlockData


  struct CommonData {
    ublas::vector<double> temperatureAtGaussPts;
    ublas::vector<double> temperatureRateAtGaussPts;
    ublas::matrix<double> gradAtGaussPts;
    inline ublas::matrix_row<ublas::matrix<double> > getGradAtGaussPts(const int gg) { 
      return ublas::matrix_row<ublas::matrix<double> >(gradAtGaussPts,gg); 
    }
  };
  CommonData commonData;

  /// \brief operator to calulete temeperature gradient at Gauss points
  struct OpGetGradAtGaussPts: public TetElementForcesAndSourcesCore::UserDataOperator {

    CommonData &commonData;
    OpGetGradAtGaussPts(const string field_name,CommonData &common_data):
      TetElementForcesAndSourcesCore::UserDataOperator(field_name),
      commonData(common_data) {}

    /** \brief operator calulating temeratire gradients
      *
      * temerature gradient is calculated multiplying direvatives of shape functions by degrees of freedom
      */
    PetscErrorCode doWork(
      int side,EntityType type,DataForcesAndSurcesCore::EntData &data) {
      PetscFunctionBegin;
      try {

        if(data.getIndices().size()==0) PetscFunctionReturn(0);
 	int nb_dofs = data.getFieldData().size();
	int nb_gauss_pts = data.getN().size1();

	//initialize
	commonData.gradAtGaussPts.resize(nb_gauss_pts,3);
	if(type == MBVERTEX) {
	  fill(commonData.gradAtGaussPts.data().begin(),commonData.gradAtGaussPts.data().end(),0);
	}	

	for(int gg = 0;gg<nb_gauss_pts;gg++) {
	  ublas::noalias(commonData.getGradAtGaussPts(gg)) += prod( trans(data.getDiffN(gg,nb_dofs)), data.getFieldData() );
	}
  
      } catch (const std::exception& ex) {
	ostringstream ss;
	ss << "throw in method: " << ex.what() << endl;
	SETERRQ(PETSC_COMM_SELF,1,ss.str().c_str());
      }

      PetscFunctionReturn(0);
    }

  };


  /** \biref operator to calculate right hand side of het conductivity terms
    * \infroup mofem_thermal_elem 
    */
  struct OpDarceysRhs: public TetElementForcesAndSourcesCore::UserDataOperator {

    BlockData &dAta;
    CommonData &commonData;
    bool useTsF;
    OpDarceysRhs(const string field_name,BlockData &data,CommonData &common_data):
      TetElementForcesAndSourcesCore::UserDataOperator(field_name),
      dAta(data),commonData(common_data),useTsF(true) {}

    Vec F;
    OpDarceysRhs(const string field_name,Vec _F,BlockData &data,CommonData &common_data):
      TetElementForcesAndSourcesCore::UserDataOperator(field_name),
      dAta(data),commonData(common_data),useTsF(false),F(_F) { }

    ublas::vector<double> Nf;

    /** \brief calculate thermal conductivity matrix
      *
      * F = int diffN^T k gard_T dOmega^2
      *
      */
    PetscErrorCode doWork(
      int side,EntityType type,DataForcesAndSurcesCore::EntData &data) {
      PetscFunctionBegin;

      if(dAta.tEts.find(getMoFEMFEPtr()->get_ent()) == dAta.tEts.end()) {
	PetscFunctionReturn(0);
      }

      try {

      if(data.getIndices().size()==0) PetscFunctionReturn(0);
      if(dAta.tEts.find(getMoFEMFEPtr()->get_ent())==dAta.tEts.end()) PetscFunctionReturn(0);

      PetscErrorCode ierr;

      int nb_row_dofs = data.getIndices().size();
      Nf.resize(nb_row_dofs);
      bzero(&*Nf.data().begin(),data.getIndices().size()*sizeof(FieldData));

      //cerr << data.getIndices() << endl;
      //cerr << data.getDiffN() << endl;

      for(unsigned int gg = 0;gg<data.getN().size1();gg++) {

	double val = (dAta.pErmeability/dAta.vIscosity)*getVolume()*getGaussPts()(3,gg);

	if(getHoGaussPtsDetJac().size()>0) {
	  val *= getHoGaussPtsDetJac()[gg]; ///< higher order geometry
	}

	//cerr << val << endl;
	//cerr << data.getDiffN() << endl;
	//cerr << data.getIndices() << endl;
	//cerr << commonData.gradAtGaussPts << endl;
	//cblas
	//cblas_dgemv(CblasRowMajor,CblasNoTrans,nb_row_dofs,3,val,
	  //&data.getDiffN()(gg,0),3,&commonData.gradAtGaussPts(gg,0),1,
	  //1.,&Nf[0],1);

	//ublas
	ublas::noalias(Nf) += val*prod(data.getDiffN(gg,nb_row_dofs),commonData.getGradAtGaussPts(gg));

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

  /** \biref operator to calculate left hand side of het conductivity terms
    * \infroup mofem_thermal_elem 
    */
  struct OpDarceysLhs: public TetElementForcesAndSourcesCore::UserDataOperator {

    BlockData &dAta;
    CommonData &commonData;
    bool useTsB;
    OpDarceysLhs(const string field_name,BlockData &data,CommonData &common_data):
      TetElementForcesAndSourcesCore::UserDataOperator(field_name),
      dAta(data),commonData(common_data),useTsB(true) { }

    Mat A;
    OpDarceysLhs(const string field_name,Mat _A,BlockData &data,CommonData &common_data):
      TetElementForcesAndSourcesCore::UserDataOperator(field_name),
      dAta(data),commonData(common_data),useTsB(false),A(_A) {}

    ublas::matrix<double> K,transK;

    /** \brief calculate thermal conductivity matrix
      *
      * K = int diffN^T k diffN^T dOmega^2
      *
      */
    PetscErrorCode doWork(
      int row_side,int col_side,
      EntityType row_type,EntityType col_type,
      DataForcesAndSurcesCore::EntData &row_data,
      DataForcesAndSurcesCore::EntData &col_data) {
      PetscFunctionBegin;

      if(dAta.tEts.find(getMoFEMFEPtr()->get_ent()) == dAta.tEts.end()) {
	PetscFunctionReturn(0);
      }

      try {
  
	if(row_data.getIndices().size()==0) PetscFunctionReturn(0);
	if(col_data.getIndices().size()==0) PetscFunctionReturn(0);
  
	int nb_row = row_data.getN().size2();
	int nb_col = col_data.getN().size2();
	K.resize(nb_row,nb_col);
	bzero(&*K.data().begin(),nb_row*nb_col*sizeof(double));
  
	for(unsigned int gg = 0;gg<row_data.getN().size1();gg++) {

	  double val = (dAta.pErmeability/dAta.vIscosity)*getVolume()*getGaussPts()(3,gg);
	  if(getHoGaussPtsDetJac().size()>0) {
	    val *= getHoGaussPtsDetJac()[gg]; ///< higher order geometry
	  }

	  //cblas
	  //double *diff_N_row,*diff_N_col;
	  //diff_N_row = &row_data.getDiffN()(gg,0);
	  //diff_N_col = &col_data.getDiffN()(gg,0);
	  //cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasTrans,
	    //nb_row,nb_col,3,
	    //val,diff_N_row,3,diff_N_col,3,1.,&K(0,0),nb_col);

	  //ublas
	  noalias(K) += val*prod(row_data.getDiffN(gg,nb_row),trans(col_data.getDiffN(gg,nb_col)));

	}

	PetscErrorCode ierr;
	if(!useTsB) {
	  const_cast<FEMethod*>(getFEMethod())->ts_B = A;
	}
	ierr = MatSetValues(
	  (getFEMethod()->ts_B),
	  nb_row,&row_data.getIndices()[0],
	  nb_col,&col_data.getIndices()[0],
	  &K(0,0),ADD_VALUES); CHKERRQ(ierr);
	if(row_side != col_side || row_type != col_type) {
	  transK.resize(nb_col,nb_row);
	  noalias(transK) = trans( K );
	  ierr = MatSetValues(
	    (getFEMethod()->ts_B),
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


  /** \brief add thermal element on tets
    * \infroup mofem_thermal_elem 
    *
    * It get data from block set and define elemenet in moab
    *
    * \param problem name
    * \param field name
    * \param name of mesh nodal positions (if not defined nodal coordinates are used)
    */
  PetscErrorCode addDarceysElements(
    const string problem_name,const string field_name,const string mesh_nodals_positions = "MESH_NODE_POSITIONS") {
    PetscFunctionBegin;

    PetscErrorCode ierr;
    ErrorCode rval;

    ierr = mField.add_finite_element("DARCEYS_FE",MF_ZERO); CHKERRQ(ierr);
    ierr = mField.modify_finite_element_add_field_row("DARCEYS_FE",field_name); CHKERRQ(ierr);
    ierr = mField.modify_finite_element_add_field_col("DARCEYS_FE",field_name); CHKERRQ(ierr);
    ierr = mField.modify_finite_element_add_field_data("DARCEYS_FE",field_name); CHKERRQ(ierr);
    if(mField.check_field(mesh_nodals_positions)) {
      ierr = mField.modify_finite_element_add_field_data("DARCEYS_FE",mesh_nodals_positions); CHKERRQ(ierr);
    }
    ierr = mField.modify_problem_add_finite_element(problem_name,"DARCEYS_FE"); CHKERRQ(ierr);

    //takes skin of block of entities
    //Skinner skin(&mField.get_moab());
    // loop over all blocksets and get data which name is FluidPressure
    for(_IT_CUBITMESHSETS_BY_BCDATA_TYPE_FOR_LOOP_(mField,BLOCKSET|MAT_MOISTURESET,it)) {
      
      Mat_Moisture moisture_data;
      ierr = it->get_attribute_data_structure(moisture_data); CHKERRQ(ierr);
      setOfBlocks[it->get_msId()].vIscosity = moisture_data.data.Viscosity;
      setOfBlocks[it->get_msId()].pErmeability = moisture_data.data.Permeability;
      rval = mField.get_moab().get_entities_by_type(it->meshset,MBTET,setOfBlocks[it->get_msId()].tEts,true); CHKERR_PETSC(rval);
      ierr = mField.add_ents_to_finite_element_by_TETs(setOfBlocks[it->get_msId()].tEts,"DARCEYS_FE"); CHKERRQ(ierr);

    }

    PetscFunctionReturn(0);
  }

  /** \brief this function is used in case of stationary problem to set elements for rhs
    * \infroup mofem_thermal_elem 
    */
  PetscErrorCode setDarceysFiniteElementRhsOperators(string field_name,Vec &F) {
    PetscFunctionBegin;
    map<int,BlockData>::iterator sit = setOfBlocks.begin();
    feRhs.get_op_to_do_Rhs().push_back(new OpGetGradAtGaussPts(field_name,commonData));
    for(;sit!=setOfBlocks.end();sit++) {
      //add finite element
      feRhs.get_op_to_do_Rhs().push_back(new OpDarceysRhs(field_name,F,sit->second,commonData));
    }
    PetscFunctionReturn(0);
  }

  /** \brief this fucntion is used in case of stationary heat conductivity problem for lhs
    * \infroup mofem_thermal_elem 
    */
  PetscErrorCode setDarceysFiniteElementLhsOperators(string field_name,Mat A) {
    PetscFunctionBegin;
    map<int,BlockData>::iterator sit = setOfBlocks.begin();
    for(;sit!=setOfBlocks.end();sit++) {
      //add finite elemen
      feLhs.get_op_to_do_Lhs().push_back(new OpDarceysLhs(field_name,A,sit->second,commonData));
    }
    PetscFunctionReturn(0);
  }


};

}

#endif //__THERMAL_ELEMENT_HPP

/***************************************************************************//**
 * \defgroup mofem_thermal_elem Thermal element
 * \ingroup mofem_forces_and_sources 
 ******************************************************************************/



