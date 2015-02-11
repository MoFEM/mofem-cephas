/** 
 * \brief Operators and data structures for thermal analys
 *
 * Implementation of nonliear eleastic element.
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

#ifndef __NONLINEAR_ELASTIC_HPP
#define __NONLINEAR_ELASTIC_HPP

#ifndef WITH_ADOL_C
  #error "MoFEM need to be compiled with ADOL-C"
#endif 

namespace MoFEM {

/** \brief struture grouping operators and data used for calulation of mass (convective) element
  * \ingroup nonlinear_eleastic_elem
  *
  * In order to assemble matrices and right hand vectors, the loops over
  * elements, enetities over that elememnts and finally loop over intergration
  * points are executed.
  *
  * Following implementation separte those three cegories of loops and to eeach
  * loop attach operator.
  *
  */
struct NonlinearElasticElement {

  /// \brief  definition of volume element
  struct MyVolumeFE: public TetElementForcesAndSourcesCore {
    MyVolumeFE(FieldInterface &_mField);
    
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
    int getRule(int order);

    PetscErrorCode preProcess(); 

  };
  
  MyVolumeFE feRhs; ///< cauclate right hand side for tetrahedral elements
  MyVolumeFE& getLoopFeRhs() { return feRhs; } ///< get rhs volume element 
  MyVolumeFE feLhs; //< calculate left hand side for tetrahedral elements
  MyVolumeFE& getLoopFeLhs() { return feLhs; } ///< get lhs volume element

  FieldInterface &mField;
  short int tAg;

  NonlinearElasticElement(
    FieldInterface &m_field,short int tag);

  /** \brief data for calulation het conductivity and heat capacity elements
    * \infroup mofem_forces_and_sources 
    */
  struct BlockData {
    int iD;
    double E;
    double PoissonRatio;
    Range tEts; ///< constatins elements in block set
  }; 
  map<int,BlockData> setOfBlocks; ///< maps block set id with appropiate BlockData

  /** \brief common data used by volume elements
    * \infroup mofem_forces_and_sources 
    */
  struct CommonData {
    map<string,vector<ublas::vector<double> > > dataAtGaussPts;
    map<string,vector<ublas::matrix<double> > > gradAtGaussPts;
    string spatialPositions;
    string meshPositions;
    vector<ublas::matrix<double> > P; ///< this need to be I Piola-Kirchoff stress tensor
    vector<vector<double*> > jacStressRowPtr;
    vector<ublas::matrix<double> > jacStress; ///< this is simply material tangent operator
  };
  CommonData commonData;

  struct OpGetDataAtGaussPts: public TetElementForcesAndSourcesCore::UserDataOperator {

    vector<ublas::vector<double> > &valuesAtGaussPts;
    vector<ublas::matrix<double> > &gradientAtGaussPts;
    const EntityType zeroAtType;

    OpGetDataAtGaussPts(const string field_name,
      vector<ublas::vector<double> > &values_at_gauss_pts,
      vector<ublas::matrix<double> > &gardient_at_gauss_pts);

    /** \brief operator calulating deformation gradient
      *
      * temerature gradient is calculated multiplying direvatives of shape functions by degrees of freedom
      */
    PetscErrorCode doWork(
      int side,EntityType type,DataForcesAndSurcesCore::EntData &data); 

  };

  struct OpGetCommonDataAtGaussPts: public OpGetDataAtGaussPts {
    OpGetCommonDataAtGaussPts(const string field_name,CommonData &common_data);

  };
 
  template<typename TYPE> 
  struct FunctionsToCalulatePiolaKirchhoffI {

    PetscErrorCode dEterminatnt(ublas::matrix<TYPE> a,TYPE &det) {
      PetscFunctionBegin;
      // a11a22a33
      //+a21a32a13
      //+a31a12a23
      //-a11a32a23
      //-a31a22a13
      //-a21a12a33
      //http://www.cg.info.hiroshima-cu.ac.jp/~miyazaki/knowledge/teche23.html
      //http://mathworld.wolfram.com/MatrixInverse.html
      det = a(0,0)*a(1,1)*a(2,2)
        +a(1,0)*a(2,1)*a(0,2)
        +a(2,0)*a(0,1)*a(1,2)
        -a(0,0)*a(2,1)*a(1,2)
        -a(2,0)*a(1,1)*a(0,2)
        -a(1,0)*a(0,1)*a(2,2);
      PetscFunctionReturn(0);
    }
  
    PetscErrorCode iNvert(TYPE det,ublas::matrix<TYPE> a,ublas::matrix<TYPE> &inv_a) {
      PetscFunctionBegin;
      //PetscErrorCode ierr;
      inv_a.resize(3,3);
      //http://www.cg.info.hiroshima-cu.ac.jp/~miyazaki/knowledge/teche23.html
      //http://mathworld.wolfram.com/MatrixInverse.html
      inv_a(0,0) = a(1,1)*a(2,2)-a(1,2)*a(2,1);
      inv_a(0,1) = a(0,2)*a(2,1)-a(0,1)*a(2,2);
      inv_a(0,2) = a(0,1)*a(1,2)-a(0,2)*a(1,1);
      inv_a(1,0) = a(1,2)*a(2,0)-a(1,0)*a(2,2);
      inv_a(1,1) = a(0,0)*a(2,2)-a(0,2)*a(2,0);
      inv_a(1,2) = a(0,2)*a(1,0)-a(0,0)*a(1,2);
      inv_a(2,0) = a(1,0)*a(2,1)-a(1,1)*a(2,0);
      inv_a(2,1) = a(0,1)*a(2,0)-a(0,0)*a(2,1);
      inv_a(2,2) = a(0,0)*a(1,1)-a(0,1)*a(1,0);
      inv_a /= det;
      PetscFunctionReturn(0);
    }

    double lambda,mu;
    ublas::matrix<TYPE> F,C,E,S,invF,P;
    TYPE J;

    int gG;
    CommonData *commonData_ptr;

    PetscErrorCode CalulateC_CauchyDefromationTensor() {
      PetscFunctionBegin;
      C.resize(3,3);
      noalias(C) = prod(trans(F),F);
      PetscFunctionReturn(0);
    }

    PetscErrorCode CalulateE_GreenStrain() {
      PetscFunctionBegin;
      E.resize(3,3);
      noalias(E) = C;
      for(int dd = 0;dd<3;dd++) {
	E(dd,dd) -= 1;
      }
      E *= 0.5;
      PetscFunctionReturn(0);
    }

    //St. Venant–Kirchhoff Material
    PetscErrorCode CalculateS_PiolaKirchhoffII() {
      PetscFunctionBegin;
      TYPE trE = 0;
      for(int dd = 0;dd<3;dd++) {
	trE += E(dd,dd);
      }
      S.resize(3,3);
      S.clear();
      for(int dd = 0;dd<3;dd++) {
	S(dd,dd) = trE*lambda;
      }
      S += 2*mu*E;
      PetscFunctionReturn(0);
    }

    virtual PetscErrorCode CalualteP_PiolaKirchhoffI(
      const BlockData block_data,
      const NumeredMoFEMFiniteElement *fe_ptr) {
      PetscFunctionBegin;
      PetscErrorCode ierr;
      lambda = LAMBDA(block_data.E,block_data.PoissonRatio);
      mu = MU(block_data.E,block_data.PoissonRatio);
      ierr = CalulateC_CauchyDefromationTensor(); CHKERRQ(ierr);
      ierr = CalulateE_GreenStrain(); CHKERRQ(ierr);
      ierr = CalculateS_PiolaKirchhoffII(); CHKERRQ(ierr);
      P.resize(3,3);
      noalias(P) = prod(F,S);
      //cerr << "P: " << P << endl;
      PetscFunctionReturn(0);
    }

    virtual PetscErrorCode SetUserActiveVariables(
      int &nb_active_variables) {
      PetscFunctionBegin;
      PetscFunctionReturn(0);
    }

    virtual PetscErrorCode SetUserActiveVariables(
      ublas::vector<double> &active_varibles) {
      PetscFunctionBegin;
      PetscFunctionReturn(0);
    }


  };

  struct OpJacobian: public TetElementForcesAndSourcesCore::UserDataOperator {

    BlockData &dAta;
    CommonData &commonData;
    FunctionsToCalulatePiolaKirchhoffI<adouble> &fUn;
    int tAg;//,lastId;
    bool jAcobian;
    bool fieldDisp;

    OpJacobian(
      const string field_name,
      BlockData &data,
      CommonData &common_data,
      FunctionsToCalulatePiolaKirchhoffI<adouble> &fun,
      int tag,bool jacobian,bool field_disp);

    ublas::vector<double> active_varibles;
    int nb_active_variables;

    PetscErrorCode doWork(
      int row_side,EntityType row_type,DataForcesAndSurcesCore::EntData &row_data); 
  };

  struct OpRhs: public TetElementForcesAndSourcesCore::UserDataOperator {

    BlockData &dAta;
    CommonData &commonData;

    OpRhs(const string field_name,BlockData &data,CommonData &common_data);

    ublas::vector<double> nf;
    PetscErrorCode doWork(
      int row_side,EntityType row_type,DataForcesAndSurcesCore::EntData &row_data);

  };

  struct OpLhs_dx: public TetElementForcesAndSourcesCore::UserDataOperator {

    BlockData &dAta;
    CommonData &commonData;
    int tAg;

    OpLhs_dx(const string vel_field,const string field_name,BlockData &data,CommonData &common_data);

    ublas::matrix<double> k,trans_k,jac,F;
    virtual PetscErrorCode getJac(DataForcesAndSurcesCore::EntData &col_data,int gg);

    PetscErrorCode doWork(
      int row_side,int col_side,
      EntityType row_type,EntityType col_type,
      DataForcesAndSurcesCore::EntData &row_data,
      DataForcesAndSurcesCore::EntData &col_data);

  };

  PetscErrorCode setBlocks();

  PetscErrorCode addElement(string element_name,
    string spatial_position_field_name,
    string material_position_field_name = "MESH_NODE_POSITIONS",bool ale = false);

  /** \brief Set opperators to calculate left hand tangent matrix and right hand residual
    *
    * \param fun class needed to calulate Piola Kirchoff I Sterss tensor
    * \param spatial_position_field_name name of appraximation field
    * \param material_position_field_name name of field to define geometry
    * \param ale true if arbitray lagrangian eulerian formulation
    * \param field_disp true if approximation field represents displacements otherwise it is field of spatial positions
    */
  PetscErrorCode setOperators(
    FunctionsToCalulatePiolaKirchhoffI<adouble> &fun,
    string spatial_position_field_name,
    string material_position_field_name = "MESH_NODE_POSITIONS",
    bool ale = false,bool field_disp = false);

};

}

#endif //__NONLINEAR_ELASTIC_HPP

/***************************************************************************//**
 * \defgroup nonlinear_eleastic_elem Non-Linear Elastic Element 
 * \ingroup mofem_forces_and_sources 
 ******************************************************************************/



