/** \file ThermalElement.hpp
 \ingroup mofem_thermal_elem

 \brief Operators and data structures for thermal analysis

 Implementation of thermal element for unsteady and steady case.
 Radiation and convection blocks implemented by Xuan Meng

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

#ifndef __THERMAL_ELEMENT_HPP
#define __THERMAL_ELEMENT_HPP

namespace MoFEM {

/** \brief structure grouping operators and data used for thermal problems
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
struct ThermalElement {

  /// \brief  definition of volume element
  struct MyVolumeFE: public VolumeElementForcesAndSourcesCore {
    MyVolumeFE(FieldInterface &_mField): VolumeElementForcesAndSourcesCore(_mField) {}

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
    int getRule(int order) { return order-1; };
  };
  MyVolumeFE feRhs; ///< cauclate right hand side for tetrahedral elements
  MyVolumeFE& getLoopFeRhs() { return feRhs; } ///< get rhs volume element
  MyVolumeFE feLhs; //< calculate left hand side for tetrahedral elements
  MyVolumeFE& getLoopFeLhs() { return feLhs; } ///< get lhs volume element


  /** \brief define surface element
    *
    * This element is used to integrate heat fluxes; convection and radiation
    */
  struct MyTriFE: public FaceElementForcesAndSourcesCore {
    MyTriFE(FieldInterface &_mField): FaceElementForcesAndSourcesCore(_mField) {}
    int getRule(int order) { return order; };
  };

  MyTriFE feFlux; //< heat flux element
  MyTriFE& getLoopFeFlux() { return feFlux; } //< get heat flux element

  MyTriFE feConvectionRhs; //< convection element
  MyTriFE feConvectionLhs;
  MyTriFE& getLoopFeConvectionRhs() { return feConvectionRhs; } //< get convection element
  MyTriFE& getLoopFeConvectionLhs() { return feConvectionLhs; }

  MyTriFE feRadiationRhs; //< radiation element
  MyTriFE feRadiationLhs;
  MyTriFE& getLoopFeRadiationRhs() { return feRadiationRhs; } //< get radiation element
  MyTriFE& getLoopFeRadiationLhs() { return feRadiationLhs; }

  FieldInterface &mField;
  ThermalElement(FieldInterface &m_field):
    feRhs(m_field),feLhs(m_field),
    feFlux(m_field),
    feConvectionRhs(m_field),feConvectionLhs(m_field),
    feRadiationRhs(m_field),feRadiationLhs(m_field),
    mField(m_field) {}

  /** \brief data for calculation heat conductivity and heat capacity elements
    * \infroup mofem_thermal_elem
    */
  struct BlockData {
    //double cOnductivity;
    ublas::matrix<double> cOnductivity_mat;  //This is (3x3) conductivity matrix
    double cApacity;   // rou * c_p == material density multiple heat capacity
    Range tEts; ///< contains elements in block set
  };
  map<int,BlockData> setOfBlocks; ///< maps block set id with appropriate BlockData

  /** \brief data for calulation heat flux
    * \infroup mofem_thermal_elem
    */
  struct FluxData {
    HeatfluxCubitBcData dAta; ///< for more details look to BCMultiIndices.hpp to see details of HeatfluxCubitBcData
    Range tRis; ///< suraface triangles where hate flux is applied
  };
  map<int,FluxData> setOfFluxes; ///< maps side set id with appropriate FluxData


  /** \brief data for convection
    * \infroup mofem_thermal_elem
    */
  struct ConvectionData {
    double cOnvection; /*The summation of Convection coefficients*/
    double tEmperature; /*Ambient temperature of the area contains the black body */
    Range tRis; ///< those will be on body skin, except this with contact with other body where temperature is applied
  };
  map<int,ConvectionData> setOfConvection; //< maps block set id with appropriate data

  /** \brief data for radiation
    * \infroup mofem_thermal_elem
    */
  struct RadiationData {
        double sIgma; /* The Stefan-Boltzmann constant*/
        double eMissivity; /* The surface emissivity coefficients range = [0,1] */
        //double aBsorption; /* The surface absorption coefficients */
        double aMbienttEmp; /* The incident radiant heat flow per unit surface area; or the ambient temperature of space*/
        Range tRis; ///< those will be on body skin, except this with contact with other body where temperature is applied
  };
  map<int,RadiationData> setOfRadiation; //< maps block set id with appropriate data

  /** \brief common data used by volume elements
    * \infroup mofem_thermal_elem
    */
  struct CommonData {
    ublas::vector<double> temperatureAtGaussPts;
    ublas::vector<double> temperatureRateAtGaussPts;
    ublas::matrix<double> gradAtGaussPts;
    inline ublas::matrix_row<ublas::matrix<double> > getGradAtGaussPts(const int gg) {
      return ublas::matrix_row<ublas::matrix<double> >(gradAtGaussPts,gg);
    }
  };
  CommonData commonData;

  /// \brief operator to calculate temperature gradient at Gauss points
  struct OpGetGradAtGaussPts: public VolumeElementForcesAndSourcesCore::UserDataOperator {

    CommonData &commonData;
    OpGetGradAtGaussPts(const string field_name,CommonData &common_data):
      VolumeElementForcesAndSourcesCore::UserDataOperator(field_name),
      commonData(common_data) {}

    /** \brief operator calculating temperature gradients
      *
      * temperature gradient is calculated multiplying derivatives of shape functions by degrees of freedom.
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

  /** \brief operator to calculate temperature  and rate of temperature at Gauss points
    * \infroup mofem_thermal_elem
    */
  template<typename OP>
  struct OpGetFieldAtGaussPts: public OP::UserDataOperator {

    ublas::vector<double> &fieldAtGaussPts;
    OpGetFieldAtGaussPts(const string field_name,ublas::vector<double> &field_at_gauss_pts):
      OP::UserDataOperator(field_name),
      fieldAtGaussPts(field_at_gauss_pts) {}

    /** \brief operator calculating temperature and rate of temperature
      *
      * temperature temperature or rate of temperature is calculated multiplying shape functions by degrees of freedom
      */
    PetscErrorCode doWork(
      int side,EntityType type,DataForcesAndSurcesCore::EntData &data) {
      PetscFunctionBegin;
      try {

        if(data.getFieldData().size()==0) PetscFunctionReturn(0);
        int nb_dofs = data.getFieldData().size();
        int nb_gauss_pts = data.getN().size1();

        //initialize
        fieldAtGaussPts.resize(nb_gauss_pts);
        if(type == MBVERTEX) {
          //loop over shape functions on entities always start from
          //vertices, so if nodal shape functions are processed, vector of
          //field values is zero at initialization
          fill(fieldAtGaussPts.begin(),fieldAtGaussPts.end(),0);
        }

        for(int gg = 0;gg<nb_gauss_pts;gg++) {
          fieldAtGaussPts[gg] += inner_prod(data.getN(gg,nb_dofs),data.getFieldData());

        }

      } catch (const std::exception& ex) {
        ostringstream ss;
        ss << "throw in method: " << ex.what() << endl;
        SETERRQ(PETSC_COMM_SELF,MOFEM_STD_EXCEPTION_THROW,ss.str().c_str());
      }

      PetscFunctionReturn(0);
    }

  };

  /** \brief operator to calculate temperature at Gauss pts
    * \infroup mofem_thermal_elem
    */
  struct OpGetTetTemperatureAtGaussPts: public OpGetFieldAtGaussPts<VolumeElementForcesAndSourcesCore> {
    OpGetTetTemperatureAtGaussPts(const string field_name,CommonData &common_data):
      OpGetFieldAtGaussPts<VolumeElementForcesAndSourcesCore>(field_name,common_data.temperatureAtGaussPts) {}
  };

  /** \brief operator to calculate temperature at Gauss pts
    * \infroup mofem_thermal_elem
    */
  struct OpGetTriTemperatureAtGaussPts: public OpGetFieldAtGaussPts<FaceElementForcesAndSourcesCore> {
    OpGetTriTemperatureAtGaussPts(const string field_name,CommonData &common_data):
      OpGetFieldAtGaussPts<FaceElementForcesAndSourcesCore>(field_name,common_data.temperatureAtGaussPts) {}
  };

  /** \brief operator to calculate temperature rate at Gauss pts
    * \infroup mofem_thermal_elem
    */
  struct OpGetTetRateAtGaussPts: public OpGetFieldAtGaussPts<VolumeElementForcesAndSourcesCore> {
    OpGetTetRateAtGaussPts(const string field_name,CommonData &common_data):
      OpGetFieldAtGaussPts<VolumeElementForcesAndSourcesCore>(field_name,common_data.temperatureRateAtGaussPts) {}
  };

  /** \biref operator to calculate right hand side of heat conductivity terms
    * \infroup mofem_thermal_elem
    */
  struct OpThermalRhs: public VolumeElementForcesAndSourcesCore::UserDataOperator {

    BlockData &dAta;
    CommonData &commonData;
    bool useTsF;
    OpThermalRhs(const string field_name,BlockData &data,CommonData &common_data):
      VolumeElementForcesAndSourcesCore::UserDataOperator(field_name),
      dAta(data),commonData(common_data),useTsF(true) {}

    Vec F;
    OpThermalRhs(const string field_name,Vec _F,BlockData &data,CommonData &common_data):
      VolumeElementForcesAndSourcesCore::UserDataOperator(field_name),
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
        Nf.clear();
        //cerr << data.getIndices() << endl;
        //cerr << data.getDiffN() << endl;

        for(unsigned int gg = 0;gg<data.getN().size1();gg++) {

          ublas::matrix<double>  val = dAta.cOnductivity_mat*getVolume()*getGaussPts()(3,gg);

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
          ublas::noalias(Nf) += prod(prod(data.getDiffN(gg,nb_row_dofs),val), commonData.getGradAtGaussPts(gg));

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

  /** \biref operator to calculate left hand side of heat conductivity terms
    * \infroup mofem_thermal_elem
    */
  struct OpThermalLhs: public VolumeElementForcesAndSourcesCore::UserDataOperator {

    BlockData &dAta;
    CommonData &commonData;
    bool useTsB;
    OpThermalLhs(const string field_name,BlockData &data,CommonData &common_data):
      VolumeElementForcesAndSourcesCore::UserDataOperator(field_name),
      dAta(data),commonData(common_data),useTsB(true) { }

    Mat A;
    OpThermalLhs(const string field_name,Mat _A,BlockData &data,CommonData &common_data):
      VolumeElementForcesAndSourcesCore::UserDataOperator(field_name),
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
	K.clear();
        for(unsigned int gg = 0;gg<row_data.getN().size1();gg++) {

          ublas::matrix<double>  val = dAta.cOnductivity_mat*getVolume()*getGaussPts()(3,gg);
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
          ublas::matrix<double> K1=prod(row_data.getDiffN(gg,nb_row),val);
          noalias(K) += prod(K1,trans(col_data.getDiffN(gg,nb_col)));
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

  /** \brief operator to calculate right hand side of heat capacity terms
    * \infroup mofem_thermal_elem
    */
  struct OpHeatCapacityRhs: public VolumeElementForcesAndSourcesCore::UserDataOperator {

    BlockData &dAta;
    CommonData &commonData;
    OpHeatCapacityRhs(const string field_name,BlockData &data,CommonData &common_data):
      VolumeElementForcesAndSourcesCore::UserDataOperator(field_name),
      dAta(data),commonData(common_data) {}

    ublas::vector<double> Nf;

    /** \brief calculate thermal conductivity matrix
      *
      * F = int N^T c (dT/dt) dOmega^2
      *
      */
    PetscErrorCode doWork(
      int side,EntityType type,DataForcesAndSurcesCore::EntData &data) {
      PetscFunctionBegin;

      try {

        if(data.getIndices().size()==0) PetscFunctionReturn(0);
        int nb_row = data.getN().size2();
        Nf.resize(nb_row);
	Nf.clear();
        for(unsigned int gg = 0;gg<data.getN().size1();gg++) {
          double val = getGaussPts()(3,gg);
          if(getHoGaussPtsDetJac().size()>0) {
            val *= getHoGaussPtsDetJac()[gg]; ///< higher order geometry
          }
          val *= commonData.temperatureRateAtGaussPts[gg];
          ////////////
          //cblas
          //cblas_daxpy(nb_row,val,&data.getN()(gg,0),1,&*Nf.data().begin(),1);
          //ublas
          ublas::noalias(Nf) += val*data.getN(gg);
        }
	Nf *= getVolume()*dAta.cApacity;
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

  /** \brief operator to calculate left hand side of heat capacity terms
    * \infroup mofem_thermal_elem
    */
  struct OpHeatCapacityLsh: public VolumeElementForcesAndSourcesCore::UserDataOperator {

    BlockData &dAta;
    CommonData &commonData;
    OpHeatCapacityLsh(const string field_name,BlockData &data,CommonData &common_data):
      VolumeElementForcesAndSourcesCore::UserDataOperator(field_name),
      dAta(data),commonData(common_data) {}

    ublas::matrix<double> M,transM;

    /** \brief calculate heat capacity matrix
      *
      * M = int N^T c N dOmega^2
      *
      */
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
	M.clear();

        for(unsigned int gg = 0;gg<row_data.getN().size1();gg++) {

          double val = getGaussPts()(3,gg);
          if(getHoGaussPtsDetJac().size()>0) {
            val *= getHoGaussPtsDetJac()[gg]; ///< higher order geometry
          }

          //cblas
          //double *N_row,*N_col;
          //N_row = &row_data.getN()(gg,0);
          //N_col = &col_data.getN()(gg,0);
          //cblas_dger(CblasRowMajor,
          //  nb_row,nb_col,val,N_row,1,N_col,1,&M(0,0),nb_col);
          //ublas
          noalias(M) += val*outer_prod( row_data.getN(gg,nb_row),col_data.getN(gg,nb_col) );

        }

	M *= getVolume()*dAta.cApacity*getFEMethod()->ts_a;

        PetscErrorCode ierr;
        ierr = MatSetValues(
               (getFEMethod()->ts_B),
               nb_row,&row_data.getIndices()[0],
               nb_col,&col_data.getIndices()[0],
               &M(0,0),ADD_VALUES); CHKERRQ(ierr);
        if(row_side != col_side || row_type != col_type) {
          transM.resize(nb_col,nb_row);
          noalias(transM) = trans(M);
          ierr = MatSetValues(
                 (getFEMethod()->ts_B),
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



  /** \brief operator for calculate heat flux and assemble to right hand side
    * \infroup mofem_thermal_elem
    */
  struct OpHeatFlux:public FaceElementForcesAndSourcesCore::UserDataOperator {

    FluxData &dAta;
    bool ho_geometry;
    bool useTsF;
    OpHeatFlux(const string field_name,FluxData &data,bool _ho_geometry = false):
      FaceElementForcesAndSourcesCore::UserDataOperator(field_name),
      dAta(data),ho_geometry(_ho_geometry),useTsF(true) { }

    Vec F;
    OpHeatFlux(const string field_name,Vec _F,
           FluxData &data,bool _ho_geometry = false):
      FaceElementForcesAndSourcesCore::UserDataOperator(field_name),
      dAta(data),ho_geometry(_ho_geometry),useTsF(false),F(_F) { }

    ublas::vector<FieldData> Nf;

    /** \brief calculate heat flux
      *
      * F = int_S N^T * flux dS
      *
      */
    PetscErrorCode doWork(
      int side,EntityType type,DataForcesAndSurcesCore::EntData &data) {
      PetscFunctionBegin;

      if(data.getIndices().size()==0) PetscFunctionReturn(0);
      if(dAta.tRis.find(getMoFEMFEPtr()->get_ent())==dAta.tRis.end()) PetscFunctionReturn(0);

      PetscErrorCode ierr;

      const FENumeredDofMoFEMEntity *dof_ptr;
      ierr = getMoFEMFEPtr()->get_row_dofs_by_petsc_gloabl_dof_idx(data.getIndices()[0],&dof_ptr); CHKERRQ(ierr);
      int rank = dof_ptr->get_max_rank();

      int nb_dofs = data.getIndices().size()/rank;

      Nf.resize(data.getIndices().size());
      Nf.clear();
      //cerr << getNormal() << endl;
      //cerr << getNormals_at_GaussPt() << endl;

      for(unsigned int gg = 0;gg<data.getN().size1();gg++) {

        double val = getGaussPts()(2,gg);
        double flux;
        if(ho_geometry) {
          double area = norm_2(getNormals_at_GaussPt(gg))*0.5; //cblas_dnrm2(3,&getNormals_at_GaussPt()(gg,0),1);
          flux = dAta.dAta.data.value1*area;  //FluxData.HeatfluxCubitBcData.data.value1 * area
        } else {
          flux = dAta.dAta.data.value1*getArea();
        }
        //cblas_daxpy(nb_row_dofs,val*flux,&data.getN()(gg,0),1,&*Nf.data().begin(),1);
        ublas::noalias(Nf) += val*flux*data.getN(gg,nb_dofs);

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


  /**
    * operator to calculate radiation therms on body surface and assemble to lhs of equations
    * for the jocabian Matrix of Picard Linearization
    * \infroup mofem_thermal_elem
    */
  struct OpRadiationLhs:public FaceElementForcesAndSourcesCore::UserDataOperator {
    CommonData &commonData; //get the temperature or temperature Rate from CommonData
    RadiationData &dAta;
    bool ho_geometry;
    bool useTsB;

    OpRadiationLhs(const string field_name,RadiationData &data,CommonData &common_data,bool _ho_geometry = false):
      FaceElementForcesAndSourcesCore::UserDataOperator(field_name),
      commonData(common_data),dAta(data),ho_geometry(_ho_geometry),useTsB(true) { }

    Mat A;
    OpRadiationLhs(const string field_name,Mat _A,RadiationData &data,CommonData &common_data,bool _ho_geometry = false):
      FaceElementForcesAndSourcesCore::UserDataOperator(field_name),
      commonData(common_data),dAta(data),ho_geometry(_ho_geometry),useTsB(false),A(_A) { }

    ublas::matrix<double> N,transN;

    /** \brief calculate thermal radiation term in the lhs of equations(Tangent Matrix) for transient Thermal Problem
      *
      * K = intS 4* N^T* sIgma* eMissivity* N*  T^3 dS (Reference _ see Finite Element Simulation of Heat Transfer
      * by jean-Michel Bergheau)
    */
    PetscErrorCode doWork(
      int row_side,int col_side,
      EntityType row_type,EntityType col_type,
      DataForcesAndSurcesCore::EntData &row_data,
      DataForcesAndSurcesCore::EntData &col_data) {
      PetscFunctionBegin;

      PetscErrorCode ierr;

      try {

	if(row_data.getIndices().size()==0) PetscFunctionReturn(0);
        if(col_data.getIndices().size()==0) PetscFunctionReturn(0);

        int nb_row = row_data.getN().size2();
        int nb_col = col_data.getN().size2();

        N.resize(nb_row,nb_col);
	N.clear();

        for(unsigned int gg = 0;gg<row_data.getN().size1();gg++) {
	  double T3_at_Gauss_pt = pow(commonData.temperatureAtGaussPts[gg],3.0);

	  double radiationConst;
	  if(ho_geometry) {
	    double area = norm_2(getNormals_at_GaussPt(gg))*0.5;
	    radiationConst = dAta.sIgma*dAta.eMissivity*area;
	  } else {
	    radiationConst = dAta.sIgma*dAta.eMissivity*getArea();
	  }
	  const double fOur = 4.0;
	  double val = fOur*getGaussPts()(2,gg)*radiationConst*T3_at_Gauss_pt;
	  noalias(N) += val*outer_prod( row_data.getN(gg,nb_row),col_data.getN(gg,nb_col) );
        }

	if(!useTsB) {
          const_cast<FEMethod*>(getFEMethod())->ts_B = A;
	}
	ierr = MatSetValues(
	  (getFEMethod()->ts_B),
	  nb_row,&row_data.getIndices()[0],
	  nb_col,&col_data.getIndices()[0],
	  &N(0,0),ADD_VALUES); CHKERRQ(ierr);
        if(row_side != col_side || row_type != col_type) {
          transN.resize(nb_col,nb_row);
          noalias(transN) = trans( N );
          ierr = MatSetValues(
                 (getFEMethod()->ts_B),
                 nb_col,&col_data.getIndices()[0],
                 nb_row,&row_data.getIndices()[0],
                 &transN(0,0),ADD_VALUES); CHKERRQ(ierr);
        }

        } catch (const std::exception& ex) {
          ostringstream ss;
          ss << "throw in method: " << ex.what() << endl;
          SETERRQ(PETSC_COMM_SELF,1,ss.str().c_str());
        }

        PetscFunctionReturn(0);
      }

  };

  /** \brief operator to calculate radiation therms on body surface and assemble to rhs of transient equations(Residual Vector)
    * \infroup mofem_thermal_elem
    */
  struct OpRadiationRhs:public FaceElementForcesAndSourcesCore::UserDataOperator {

    CommonData &commonData; //get the temperature or temperature Rate from CommonData
    RadiationData &dAta;
    bool ho_geometry;
    bool useTsF;
    OpRadiationRhs(const string field_name,RadiationData &data,CommonData &common_data,bool _ho_geometry = false):
      FaceElementForcesAndSourcesCore::UserDataOperator(field_name),
      commonData(common_data),dAta(data),ho_geometry(_ho_geometry),useTsF(true) {}

    Vec F;
    OpRadiationRhs(const string field_name,Vec _F,RadiationData &data,CommonData &common_data,bool _ho_geometry = false):
      FaceElementForcesAndSourcesCore::UserDataOperator(field_name),
      commonData(common_data),dAta(data),ho_geometry(_ho_geometry),useTsF(false),F(_F) {}

    ublas::vector<FieldData> Nf;

    /** \brief calculate Transient Radiation condition on the right hand side residual
      *
      *  R=int_S N^T * sIgma * eMissivity * (Ta^4 -Ts^4) dS
     **/
    PetscErrorCode doWork(int side,EntityType type,DataForcesAndSurcesCore::EntData &data) {
      PetscFunctionBegin;

      if(data.getIndices().size()==0) PetscFunctionReturn(0);
      if(dAta.tRis.find(getMoFEMFEPtr()->get_ent())==dAta.tRis.end()) PetscFunctionReturn(0);

      PetscErrorCode ierr;

      const FENumeredDofMoFEMEntity *dof_ptr;
      ierr = getMoFEMFEPtr()->get_row_dofs_by_petsc_gloabl_dof_idx(data.getIndices()[0],&dof_ptr); CHKERRQ(ierr);
      int rank = dof_ptr->get_max_rank();
      int nb_row_dofs = data.getIndices().size()/rank;

      Nf.resize(data.getIndices().size());
      Nf.clear();
      //cerr << getNormal() << endl;
      //cerr << getNormals_at_GaussPt() << endl;

      for(unsigned int gg = 0;gg<data.getN().size1();gg++) {

	double T4_at_Gauss_pt = pow(commonData.temperatureAtGaussPts[gg],4.0);
	double ambientTemp = pow(dAta.aMbienttEmp,4.0);
	double tEmp = 0;

	if(ambientTemp > 0) {
	  tEmp = -ambientTemp + T4_at_Gauss_pt;
	}

	double val = getGaussPts()(2,gg);
	double radiationConst;

	if(ho_geometry) {
	  double area = norm_2(getNormals_at_GaussPt(gg))*0.5;
          radiationConst = dAta.sIgma*dAta.eMissivity*tEmp*area;
        } else {
          radiationConst = dAta.sIgma*dAta.eMissivity*tEmp*getArea();
        }
        ublas::noalias(Nf) += val*radiationConst*data.getN(gg,nb_row_dofs);

      }

      //cerr << "VecSetValues\n";
      //cerr << Nf << endl;
      //cerr << data.getIndices() << endl;

      if(useTsF) {
	ierr = VecSetValues(getFEMethod()->ts_F,data.getIndices().size(),&data.getIndices()[0],&Nf[0],ADD_VALUES); CHKERRQ(ierr);
      } else {
        ierr = VecSetValues(F,data.getIndices().size(),&data.getIndices()[0],&Nf[0],ADD_VALUES); CHKERRQ(ierr);
      }

      PetscFunctionReturn(0);
    }

  };

  /** \brief operator to calculate convection therms on body surface and assemble to rhs of equations
    * \infroup mofem_thermal_elem
    */
  struct OpConvectionRhs:public FaceElementForcesAndSourcesCore::UserDataOperator {

    CommonData &commonData; //get the temperature or temperature Rate from CommonData
    ConvectionData &dAta;
    bool ho_geometry;
    bool useTsF;
    OpConvectionRhs(const string field_name,ConvectionData &data,CommonData &common_data,bool _ho_geometry = false):
      FaceElementForcesAndSourcesCore::UserDataOperator(field_name),
      commonData(common_data),dAta(data),ho_geometry(_ho_geometry),useTsF(true) {}

    Vec F;
    OpConvectionRhs(const string field_name,Vec _F,ConvectionData &data,CommonData &common_data,bool _ho_geometry = false):
      FaceElementForcesAndSourcesCore::UserDataOperator(field_name),
      commonData(common_data),dAta(data),ho_geometry(_ho_geometry),useTsF(false),F(_F) {}

    ublas::vector<FieldData> Nf;

    /** brief calculate Convection condition on the right hand side
      *  R=int_S N^T*alpha*N_f  dS **/

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
      Nf.clear();
      //fill(Nf.begin(),Nf.end(),0);
      //cerr << getNormal() << endl;
      //cerr << getNormals_at_GaussPt() << endl;

      for(unsigned int gg = 0;gg<data.getN().size1();gg++) {

	double T_at_Gauss_pt = commonData.temperatureAtGaussPts[gg];
	double convectionConst;
	if(ho_geometry) {
	  double area = norm_2(getNormals_at_GaussPt(gg))*0.5;
	  convectionConst = dAta.cOnvection*area*(T_at_Gauss_pt-dAta.tEmperature);
	} else {
	  convectionConst = dAta.cOnvection*getArea()*(T_at_Gauss_pt-dAta.tEmperature);
	}
	double val = getGaussPts()(2,gg)*convectionConst;
	ublas::noalias(Nf) += val*data.getN(gg,nb_row_dofs);

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


  /// \biref operator to calculate convection therms on body surface and assemble to lhs of equations
  struct OpConvectionLhs:public FaceElementForcesAndSourcesCore::UserDataOperator {

    ConvectionData &dAta;
    bool ho_geometry;
    bool useTsB;

    OpConvectionLhs(const string field_name,
            ConvectionData &data,bool _ho_geometry = false):
      FaceElementForcesAndSourcesCore::UserDataOperator(field_name),
      dAta(data),ho_geometry(_ho_geometry),useTsB(true) {}

    Mat A;
    OpConvectionLhs(const string field_name,Mat _A,
            ConvectionData &data,bool _ho_geometry = false):
      FaceElementForcesAndSourcesCore::UserDataOperator(field_name),
      dAta(data),ho_geometry(_ho_geometry),useTsB(false),A(_A) {}

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

        if(row_data.getIndices().size()==0) PetscFunctionReturn(0);
        if(col_data.getIndices().size()==0) PetscFunctionReturn(0);

        int nb_row = row_data.getN().size2();
        int nb_col = col_data.getN().size2();
        K.resize(nb_row,nb_col);
	K.clear();

        for(unsigned int gg = 0;gg<row_data.getN().size1();gg++) {

          double convectionConst;
          if(ho_geometry) {
            double area = norm_2(getNormals_at_GaussPt(gg))*0.5;
            convectionConst = dAta.cOnvection*area;
          }   else {
            convectionConst = dAta.cOnvection*getArea();
          }
          double val = getGaussPts()(2,gg)*convectionConst;
          noalias(K) += val*outer_prod( row_data.getN(gg,nb_row),col_data.getN(gg,nb_col) );

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

  /** \brief this calas is to control time stepping
    * \infroup mofem_thermal_elem
    *
    * It is used to save data for temperature rate vector to MoFEM field.
    */
  struct UpdateAndControl: public FEMethod {

    FieldInterface& mField;
    const string tempName;
    const string rateName;

    UpdateAndControl(FieldInterface& _mField,
      const string temp_name,const string rate_name):
      mField(_mField),tempName(temp_name),rateName(rate_name) {}

    PetscErrorCode preProcess() {
      PetscFunctionBegin;
      PetscErrorCode ierr;
      ierr = mField.set_other_local_ghost_vector(
        problemPtr,tempName,rateName,ROW,ts_u_t,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
      PetscFunctionReturn(0);
    }

    PetscErrorCode postProcess() {
      PetscFunctionBegin;
      PetscFunctionReturn(0);
    }

  };

  /** \brief TS monitore it records temperature at time steps
    * \infroup mofem_thermal_elem
    */
  struct TimeSeriesMonitor: public FEMethod {

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

      ierr = mField.set_global_ghost_vector(
             problemPtr,ROW,ts_u,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);

      BitRefLevel proble_bit_level = problemPtr->get_BitRefLevel();

      SeriesRecorder *recorder_ptr = NULL;
      ierr = mField.query_interface(recorder_ptr); CHKERRQ(ierr);
      ierr = recorder_ptr->record_begin(seriesName); CHKERRQ(ierr);
      ierr = recorder_ptr->record_field(seriesName,tempName,proble_bit_level,mask); CHKERRQ(ierr);
      ierr = recorder_ptr->record_end(seriesName,ts_t); CHKERRQ(ierr);

      PetscFunctionReturn(0);
    }

  };

  /** \brief add thermal element on tets
    * \infroup mofem_thermal_elem
    *
    * It get data from block set and define element in moab
    *w
    * \param field name
    * \param name of mesh nodal positions (if not defined nodal coordinates are used)
    */
  PetscErrorCode addThermalElements(const string field_name,const string mesh_nodals_positions = "MESH_NODE_POSITIONS");

  /** \brief add heat flux element
    * \infroup mofem_thermal_elem
    *
    * It get data from heat flux set and define element in moab. Alternatively
    * uses block set with name HEAT_FLUX.
    *
    * \param field name
    * \param name of mesh nodal positions (if not defined nodal coordinates are used)
    */
  PetscErrorCode addThermalFluxElement(const string field_name,const string mesh_nodals_positions = "MESH_NODE_POSITIONS");


  /** \brief add convection element
  * \infroup mofem_thermal_elem
  *
  * It get data from convection set and define element in moab. Alternatively
  * uses block set with name CONVECTION.
  *
  * \param field name
  * \param name of mesh nodal positions (if not defined nodal coordinates are used)
  */
  PetscErrorCode addThermalConvectionElement(const string field_name,const string mesh_nodals_positions = "MESH_NODE_POSITIONS");

  /** \brief add Non-linear Radiation element
  * \infroup mofem_thermal_elem
  *
  * It get data from Radiation set and define element in moab. Alternatively
  * uses block set with name RADIATION.
  *
  * \param field name
  * \param name of mesh nodal positions (if not defined nodal coordinates are used)
  */
  PetscErrorCode addThermalRadiationElement(const string field_name,const string mesh_nodals_positions = "MESH_NODE_POSITIONS");

  /** \brief this function is used in case of stationary problem to set elements for rhs
    * \infroup mofem_thermal_elem
    */
  PetscErrorCode setThermalFiniteElementRhsOperators(string field_name,Vec &F);

  /** \brief this function is used in case of stationary heat conductivity problem for lhs
    * \infroup mofem_thermal_elem
    */
  PetscErrorCode setThermalFiniteElementLhsOperators(string field_name,Mat A);

  /** \brief this function is used in case of stationary problem for heat flux terms
    * \infroup mofem_thermal_elem
    */
  PetscErrorCode setThermalFluxFiniteElementRhsOperators(string field_name,Vec &F,const string mesh_nodals_positions = "MESH_NODE_POSITIONS");

  /* \brief linear Steady convection terms in lhs
   */
  PetscErrorCode setThermalConvectionFiniteElementRhsOperators(string field_name,Vec &F,const string mesh_nodals_positions = "MESH_NODE_POSITIONS");

  /* \brief linear Steady convection terms in rhs
   */
  PetscErrorCode setThermalConvectionFiniteElementLhsOperators(string field_name,Mat A,const string mesh_nodals_positions = "MESH_NODE_POSITIONS");

  /** \brief set up operators for unsteady heat flux; convection; radiation problem
    * \infroup mofem_thermal_elem
    */
  PetscErrorCode setTimeSteppingProblem(string field_name,string rate_name,const string mesh_nodals_positions = "MESH_NODE_POSITIONS");

  /** \brief set up operators for unsteady heat flux; convection; radiation problem
    * \infroup mofem_thermal_elem
    */
  PetscErrorCode setTimeSteppingProblem(TsCtx &ts_ctx,string field_name,string rate_name,const string mesh_nodals_positions = "MESH_NODE_POSITIONS");

};

}

#endif //__THERMAL_ELEMENT_HPP

/***************************************************************************//**
 * \defgroup mofem_thermal_elem Thermal element
 * \ingroup mofem_forces_and_sources
 ******************************************************************************/
