/** \file NonLienarElasticElement.hpp
 * \ingroup nonlinear_elastic_elem
 * \brief Operators and data structures for non-linear elastic analysis
 *
 * Implementation of nonlinear elastic element.
 */

/*
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

/** \brief structure grouping operators and data used for calculation of nonlinear elastic element
  * \ingroup nonlinear_elastic_elem
  *
  * In order to assemble matrices and right hand vectors, the loops over
  * elements, entities over that elements and finally loop over integration
  * points are executed.
  *
  * Following implementation separate those three categories of loops and to each
  * loop attach operator.
  *
  */
struct NonlinearElasticElement {

  /// \brief  definition of volume element
  struct MyVolumeFE: public VolumeElementForcesAndSourcesCore {

    Mat A;
    Vec F;

    int addToRule;

    MyVolumeFE(FieldInterface &m_field);

    /** \brief it is used to calculate nb. of Gauss integration points
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

    Vec V;
    double eNergy;

    PetscErrorCode preProcess();
    PetscErrorCode postProcess();

  };

  MyVolumeFE feRhs; ///< calculate right hand side for tetrahedral elements
  MyVolumeFE& getLoopFeRhs() { return feRhs; } ///< get rhs volume element
  MyVolumeFE feLhs; //< calculate left hand side for tetrahedral elements
  MyVolumeFE& getLoopFeLhs() { return feLhs; } ///< get lhs volume element

  MyVolumeFE feEnergy; ///< calculate elastic energy
  MyVolumeFE& getLoopFeEnergy() { return feEnergy; } ///< get energy fe

  FieldInterface &mField;
  short int tAg;

  NonlinearElasticElement(
    FieldInterface &m_field,short int tag);

  template<typename TYPE>
  struct FunctionsToCalulatePiolaKirchhoffI;

  /** \brief data for calculation het conductivity and heat capacity elements
    * \ingroup nonlinear_elastic_elem
    */
  struct BlockData {
    int iD;
    double E;
    double PoissonRatio;
    Range tEts; ///< constatins elements in block set
    FunctionsToCalulatePiolaKirchhoffI<adouble> *materialAdoublePtr;
    FunctionsToCalulatePiolaKirchhoffI<double> *materialDoublePtr;
    Range forcesOnlyOnEntitiesRow;
    Range forcesOnlyOnEntitiesCol;
  };
  map<int,BlockData> setOfBlocks; ///< maps block set id with appropriate BlockData

  /** \brief common data used by volume elements
    * \ingroup nonlinear_elastic_elem
    */
  struct CommonData {
    map<string,vector<VectorDouble > > dataAtGaussPts;
    map<string,vector<MatrixDouble > > gradAtGaussPts;
    string spatialPositions;
    string meshPositions;
    vector<MatrixDouble > sTress;
    vector<vector<double*> > jacStressRowPtr;
    vector<MatrixDouble > jacStress; ///< this is simply material tangent operator
   };
  CommonData commonData;

  /** \brief Implementation of elastic (non-linear) St. Kirchoff equation
    * \ingroup nonlinear_elastic_elem
    */
  template<typename TYPE>
  struct FunctionsToCalulatePiolaKirchhoffI {

    /** \brief Calculate determinant of 3x3 matrix
      */
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


    /** \brief Calculate inverse of 3x3 matrix
      */
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
    ublas::matrix<TYPE> F,C,E,S,invF,P,SiGma;
    TYPE J,eNergy;

    int gG;	///< Gauss point number
    CommonData *commonDataPtr; ///< common data shared between entities (f.e. field values at Gauss pts.)
    VolumeElementForcesAndSourcesCore::UserDataOperator *opPtr; ///< pointer to finite element tetrahedral operator

    PetscErrorCode calculateC_CauchyDefromationTensor() {
      PetscFunctionBegin;
      C.resize(3,3);
      noalias(C) = prod(trans(F),F);
      PetscFunctionReturn(0);
    }

    PetscErrorCode calculateE_GreenStrain() {
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
    PetscErrorCode calculateS_PiolaKirchhoffII() {
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

    /** \brief Function overload to implement user material
      *

      * Calculation of Piola Kirchhoff I is implemented by user. Tangent matrix
      * user implemented physical equation is calculated using automatic
      * differentiation.

      * \f$\mathbf{S} = \lambda\textrm{tr}[\mathbf{E}]\mathbf{I}+2\mu\mathbf{E}\f$

      * Notes: <br>
      * Number of actual Gauss point is accessed from variable gG. <br>
      * Access to operator data structures is available by variable opPtr. <br>
      * Access to common data is by commonDataPtr. <br>

      * \param block_data used to give access to material parameters
      * \param fe_ptr pointer to element data structures

      For details look to: <br>
      NONLINEAR CONTINUUM MECHANICS FOR FINITE ELEMENT ANALYSIS, Javier Bonet,
      Richard D. Wood

      */
    virtual PetscErrorCode calculateP_PiolaKirchhoffI(
      const BlockData block_data,
      const NumeredMoFEMFiniteElement *fe_ptr
    ) {
      PetscFunctionBegin;
      PetscErrorCode ierr;
      lambda = LAMBDA(block_data.E,block_data.PoissonRatio);
      mu = MU(block_data.E,block_data.PoissonRatio);
      ierr = calculateC_CauchyDefromationTensor(); CHKERRQ(ierr);
      ierr = calculateE_GreenStrain(); CHKERRQ(ierr);
      ierr = calculateS_PiolaKirchhoffII(); CHKERRQ(ierr);
      P.resize(3,3);
      noalias(P) = prod(F,S);
      //cerr << "P: " << P << endl;
      PetscFunctionReturn(0);
    }

    virtual PetscErrorCode setUserActiveVariables(
      int &nb_active_variables) {
      PetscFunctionBegin;
      PetscFunctionReturn(0);
    }

    virtual PetscErrorCode setUserActiveVariables(
      VectorDouble &active_varibles) {
      PetscFunctionBegin;
      PetscFunctionReturn(0);
    }

    /** \brief Calculate elastic energy density
      *
      * \f[\Psi = \frac{1}{2}\lambda(\textrm{tr}[\mathbf{E}])^2+\mu\mathbf{E}:\mathbf{E}\f]
      */
    virtual PetscErrorCode calculateElasticEnergy(
      const BlockData block_data,
      const NumeredMoFEMFiniteElement *fe_ptr
    ) {
      PetscFunctionBegin;
      PetscErrorCode ierr;
      lambda = LAMBDA(block_data.E,block_data.PoissonRatio);
      mu = MU(block_data.E,block_data.PoissonRatio);
      ierr = calculateC_CauchyDefromationTensor(); CHKERRQ(ierr);
      ierr = calculateE_GreenStrain(); CHKERRQ(ierr);
      TYPE trace = 0;
      eNergy = 0;
      for(int ii = 0;ii<3;ii++) {
        trace += E(ii,ii);
        for(int jj = 0;jj<3;jj++) {
          TYPE e = E(ii,jj);
          eNergy += mu*e*e;
        }
      }
      eNergy += 0.5*lambda*trace*trace;
      PetscFunctionReturn(0);
    }

    /** \brief Calculate Eshelby stress
    */
    virtual PetscErrorCode calculateSiGma_EshelbyStress(
      const BlockData block_data,
      const NumeredMoFEMFiniteElement *fe_ptr
    ) {
      PetscFunctionBegin;
      PetscErrorCode ierr;
      ierr = calculateP_PiolaKirchhoffI(block_data,fe_ptr); CHKERRQ(ierr);
      ierr = calculateElasticEnergy(block_data,fe_ptr); CHKERRQ(ierr);
      SiGma.resize(3,3,false);
      noalias(SiGma) = -prod(trans(F),P);
      for(int dd = 0;dd<3;dd++) {
        SiGma(dd,dd) += eNergy;
      }
      PetscFunctionReturn(0);
    }

  };

  struct OpGetDataAtGaussPts: public VolumeElementForcesAndSourcesCore::UserDataOperator {

    vector<VectorDouble > &valuesAtGaussPts;
    vector<MatrixDouble > &gradientAtGaussPts;
    const EntityType zeroAtType;

    OpGetDataAtGaussPts(const string field_name,
      vector<VectorDouble > &values_at_gauss_pts,
      vector<MatrixDouble > &gardient_at_gauss_pts);

    /** \brief operator calculating deformation gradient
      *
      * temperature gradient is calculated multiplying derivatives of shape functions by degrees of freedom
      */
    PetscErrorCode doWork(
      int side,EntityType type,DataForcesAndSurcesCore::EntData &data);

  };

  struct OpGetCommonDataAtGaussPts: public OpGetDataAtGaussPts {
    OpGetCommonDataAtGaussPts(const string field_name,CommonData &common_data);
  };

  struct OpJacobianPiolaKirchhoffStress: public VolumeElementForcesAndSourcesCore::UserDataOperator {

    BlockData &dAta;
    CommonData &commonData;
    int tAg;//,lastId;
    int adlocReturnValue;
    bool jAcobian;
    bool fUnction;
    bool aLe;
    bool fieldDisp;

    /**
      \brief Construct operator to calculate Pilo Kirchhoff stress or its derivatives over gradient deformation

      \param field_name approximation field name of spatial positions or displacements
      \param data reference to block data (what is Young modulus, Poisson ratio or what elements are part of the block)
      \param tag adol-c unique tag of the tape
      \param jacobian if true derivative of Piola Stress is calculated otherwise just stress is calculated
      \param field_disp if true approximation field keeps displacements not spatial positions

    */
    OpJacobianPiolaKirchhoffStress(
      const string field_name,
      BlockData &data,
      CommonData &common_data,
      int tag,
      bool jacobian,
      bool ale,
      bool field_disp
    );

    VectorDouble active_varibles;
    int nb_active_variables;

    vector<MatrixDouble > *ptrh;
    vector<MatrixDouble > *ptrH;

    adouble detH;
    ublas::matrix<adouble> h;
    ublas::matrix<adouble> H;
    ublas::matrix<adouble> invH;

    virtual PetscErrorCode calculateStress();

    PetscErrorCode doWork(int row_side,EntityType row_type,DataForcesAndSurcesCore::EntData &row_data);

  };

  struct OpRhsPiolaKirchhoff: public VolumeElementForcesAndSourcesCore::UserDataOperator {

    BlockData &dAta;
    CommonData &commonData;
    bool fieldDisp;
    bool aLe;

    ublas::vector<int> iNdices;
    OpRhsPiolaKirchhoff(const string field_name,BlockData &data,CommonData &common_data);

    VectorDouble nf;
    PetscErrorCode doWork(
      int row_side,EntityType row_type,DataForcesAndSurcesCore::EntData &row_data
    );

  };

  struct OpEnergy: public VolumeElementForcesAndSourcesCore::UserDataOperator {

    BlockData &dAta;
    CommonData &commonData;
    Vec *Vptr;
    bool fieldDisp;

    OpEnergy(const string field_name,BlockData &data,CommonData &common_data,Vec *v_ptr,bool field_disp);

    PetscErrorCode doWork(
      int row_side,EntityType row_type,DataForcesAndSurcesCore::EntData &row_data);

  };

  struct OpLhsPiolaKirchhoff_dx: public VolumeElementForcesAndSourcesCore::UserDataOperator {

    BlockData &dAta;
    CommonData &commonData;
    int tAg;
    bool aLe;

    ublas::vector<int> rowIndices;
    ublas::vector<int> colIndices;

    OpLhsPiolaKirchhoff_dx(const string vel_field,const string field_name,BlockData &data,CommonData &common_data);

    MatrixDouble k,trans_k,jac,F;

    /**
      \brief Directive of Piola Kirchhoff stress over spatial DOFs

      This project derivative \f$\frac{\partial P}{\partial F}\f$, that is
      \f[
      \frac{\partial P}{\partial x_\textrm{DOF}} =  \frac{\partial P}{\partial F}\frac{\partial F}{\partial x_\textrm{DOF}},
      \f]
      where second therm \f$\frac{\partial F}{\partial x_\textrm{DOF}}\f$ is derivative of shape function

    */
    virtual PetscErrorCode getJac(DataForcesAndSurcesCore::EntData &col_data,int gg);

    virtual PetscErrorCode aSemble(
      int row_side,int col_side,
      EntityType row_type,EntityType col_type,
      DataForcesAndSurcesCore::EntData &row_data,
      DataForcesAndSurcesCore::EntData &col_data
    );

    PetscErrorCode doWork(
      int row_side,int col_side,
      EntityType row_type,EntityType col_type,
      DataForcesAndSurcesCore::EntData &row_data,
      DataForcesAndSurcesCore::EntData &col_data);

  };


  struct OpLhsPiolaKirchhoff_dX: public OpLhsPiolaKirchhoff_dx {

    OpLhsPiolaKirchhoff_dX(const string vel_field,const string field_name,BlockData &data,CommonData &common_data);

    /// \berief Derivative of Piola Kirchhoff stress over material DOFs
    PetscErrorCode getJac(DataForcesAndSurcesCore::EntData &col_data,int gg);

    PetscErrorCode aSemble(
      int row_side,int col_side,
      EntityType row_type,EntityType col_type,
      DataForcesAndSurcesCore::EntData &row_data,
      DataForcesAndSurcesCore::EntData &col_data
    );

  };

  struct OpJacobianEshelbyStress: public OpJacobianPiolaKirchhoffStress {

    OpJacobianEshelbyStress(
      const string field_name,
      BlockData &data,
      CommonData &common_data,
      int tag,
      bool jacobian,
      bool ale
    );

    PetscErrorCode calculateStress();

  };

  struct OpRhsEshelbyStrees: public OpRhsPiolaKirchhoff {

    OpRhsEshelbyStrees(
      const string field_name,BlockData &data,CommonData &common_data
    );

  };

  struct OpLhsEshelby_dx: public OpLhsPiolaKirchhoff_dX {

    OpLhsEshelby_dx(
      const string vel_field,const string field_name,BlockData &data,CommonData &common_data
    );

    PetscErrorCode getJac(DataForcesAndSurcesCore::EntData &col_data,int gg);

  };

  struct OpLhsEshelby_dX: public OpLhsPiolaKirchhoff_dx {

    OpLhsEshelby_dX(
      const string vel_field,const string field_name,BlockData &data,CommonData &common_data
    );

    PetscErrorCode getJac(DataForcesAndSurcesCore::EntData &col_data,int gg);

  };

  PetscErrorCode setBlocks(
    FunctionsToCalulatePiolaKirchhoffI<double> *materialDoublePtr,
    FunctionsToCalulatePiolaKirchhoffI<adouble> *materialAdoublePtr);

  PetscErrorCode addElement(string element_name,
    string spatial_position_field_name,
    string material_position_field_name = "MESH_NODE_POSITIONS",bool ale = false);

  /** \brief Set operators to calculate left hand tangent matrix and right hand residual
    *
    * \param fun class needed to calculate Piola Kirchoff I Stress tensor
    * \param spatial_position_field_name name of approximation field
    * \param material_position_field_name name of field to define geometry
    * \param ale true if arbitrary Lagrangian Eulerian formulation
    * \param field_disp true if approximation field represents displacements otherwise it is field of spatial positions
    */
  PetscErrorCode setOperators(
    string spatial_position_field_name,
    string material_position_field_name = "MESH_NODE_POSITIONS",
    bool ale = false,bool field_disp = false);

};

}

#endif //__NONLINEAR_ELASTIC_HPP

/***************************************************************************//**
 * \defgroup nonlinear_elastic_elem NonLinear Elastic Element
 * \ingroup mofem_forces_and_sources
 ******************************************************************************/
