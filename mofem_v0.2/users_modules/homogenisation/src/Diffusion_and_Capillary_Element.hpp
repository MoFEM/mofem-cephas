/** \file MoistureElement.hpp
 * \brief Operators and data structures for moisture analys
 *
 * Implementation of moisture element for unsteady and steady case.
 *
 */

/* Copyright (C) 2014, Zahur Ullah (Zahur.Ullah AT glasgow.ac.uk)
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

#ifndef __MOISTURE_ELEMENT_HPP
#define __MOISTURE_ELEMENT_HPP

//#include<moab/Skinner.hpp>

namespace MoFEM {
  
  /** \brief struture grouping operators and data used for moisture problems
   * \ingroup mofem_moisture_elem
   *
   * In order to assemble matrices and right hand vectors, the loops over
   * elements, enetities over that elememnts and finally loop over intergration
   * points are executed.
   *
   * Following implementation separte those three cegories of loops and to eeach
   * loop attach operator.
   *
   */
  
  struct MoistureElement {

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
    
    
    MyVolumeFE feRhsDiffusion; ///< cauclate right hand side for tetrahedral elements
    MyVolumeFE& getLoopFeRhs() { return feRhsDiffusion; } ///< get rhs volume element
    
    MyVolumeFE feLhsDiffusion, feLhsCapillary; //< calculate left hand side for tetrahedral elements
    MyVolumeFE& getLoopFeLhsDiffusion() { return feLhsDiffusion; } ///< get lhs volume element
    MyVolumeFE& getLoopFeLhsCapillary() { return feLhsCapillary; } ///< get lhs volume element

    /** \brief define surface element
     *
     * This element is used to integrate het fluxes and radiation
     */
    struct MyTriFE: public TriElementForcesAndSurcesCore {
      MyTriFE(FieldInterface &_mField): TriElementForcesAndSurcesCore(_mField) {}
      int getRule(int order) { return ceil(order/2); };
    };
    
    MyTriFE feFluxDiffusion; //< heat flux element
    MyTriFE& getLoopFeFluxDiffusion() { return feFluxDiffusion; } //< get heat flux element

    
    FieldInterface &mField;
    MoistureElement(FieldInterface &m_field):
    feRhsDiffusion(m_field), feLhsDiffusion(m_field),feLhsCapillary(m_field), feFluxDiffusion(m_field), mField(m_field){}
    
     /** \brief data for calulation het conductivity and heat capacity elements
     * \infroup mofem_moisture_elem
     */
    struct BlockData {
      double dIffusivity;
      double vIscosity;
      double pErmeability;
      Range tEts; ///< constatins elements in block set
    };
    map<int,BlockData> setOfBlocks; ///< maps block set id with appropiate BlockData
    
    
    /** \brief data for calulation heat flux
     * \infroup mofem_thermal_elem
     */
    struct FluxData {
      double dAta_mass_flux;
      Range tRis; ///< suraface triangles where hate flux is applied
    };
    map<int,FluxData> setOfFluxes; ///< maps side set id with appropiate FluxData

    
    /** \brief add moisture element on tets
     * \infroup mofem_moisture_elem
     *
     * It get data from block set and define elemenet in moab
     *
     * \param problem name
     * \param field name
     * \param name of mesh nodal positions (if not defined nodal coordinates are used)
     */
    PetscErrorCode addDiffusionElements(const string problem_name,const string field_name,const string mesh_nodals_positions = "MESH_NODE_POSITIONS") {
      PetscFunctionBegin;
      PetscErrorCode ierr;
      ErrorCode rval;
//      cout<<"insides the addDiffusionElements = "<<endl;
      ierr = mField.add_finite_element("DIFFUSION_FE",MF_ZERO); CHKERRQ(ierr);
      ierr = mField.modify_finite_element_add_field_row("DIFFUSION_FE",field_name); CHKERRQ(ierr);
      ierr = mField.modify_finite_element_add_field_col("DIFFUSION_FE",field_name); CHKERRQ(ierr);
      ierr = mField.modify_finite_element_add_field_data("DIFFUSION_FE",field_name); CHKERRQ(ierr);
      if(mField.check_field(mesh_nodals_positions)) {
        ierr = mField.modify_finite_element_add_field_data("DIFFUSION_FE",mesh_nodals_positions); CHKERRQ(ierr);
      }
      ierr = mField.modify_problem_add_finite_element(problem_name,"DIFFUSION_FE"); CHKERRQ(ierr);
      
        // loop over all blocksets
      for(_IT_CUBITMESHSETS_BY_BCDATA_TYPE_FOR_LOOP_(mField,BLOCKSET|MAT_MOISTURESET,it)) {
        if(it->get_Cubit_name().compare(0,19,"MAT_MOISTURE") == 0){
          Mat_Moisture diffusion_data;
          ierr = it->get_attribute_data_structure(diffusion_data); CHKERRQ(ierr);
//          cout<<"diffusion_data.data.Diffusivity = "<<diffusion_data.data.Diffusivity<<endl;
//          cout<<"it->get_msId() = "<<it->get_msId()<<endl;
          setOfBlocks[it->get_msId()].dIffusivity = diffusion_data.data.Diffusivity;
          rval = mField.get_moab().get_entities_by_type(it->meshset,MBTET,setOfBlocks[it->get_msId()].tEts,true); CHKERR_PETSC(rval);
//          cout<<"setOfBlocks[it->get_msIdx()].tEts.size() = "<<setOfBlocks[it->get_msId()].tEts.size()<<endl;
          ierr = mField.add_ents_to_finite_element_by_TETs(setOfBlocks[it->get_msId()].tEts,"DIFFUSION_FE"); CHKERRQ(ierr);
        }
      }
      
      PetscFunctionReturn(0);
    }
    
    
    PetscErrorCode addCapillaryElements(const string problem_name,const string field_name,const string mesh_nodals_positions = "MESH_NODE_POSITIONS") {
      PetscFunctionBegin;
      PetscErrorCode ierr;
      ErrorCode rval;
      cout<<"insides the addCapillaryElements = "<<endl;
      ierr = mField.add_finite_element("CAPILLARY_FE",MF_ZERO); CHKERRQ(ierr);
      ierr = mField.modify_finite_element_add_field_row("CAPILLARY_FE",field_name); CHKERRQ(ierr);
      ierr = mField.modify_finite_element_add_field_col("CAPILLARY_FE",field_name); CHKERRQ(ierr);
      ierr = mField.modify_finite_element_add_field_data("CAPILLARY_FE",field_name); CHKERRQ(ierr);
      if(mField.check_field(mesh_nodals_positions)) {
        ierr = mField.modify_finite_element_add_field_data("CAPILLARY_FE",mesh_nodals_positions); CHKERRQ(ierr);
      }
      ierr = mField.modify_problem_add_finite_element(problem_name,"CAPILLARY_FE"); CHKERRQ(ierr);
      
      // loop over all blocksets
      for(_IT_CUBITMESHSETS_BY_BCDATA_TYPE_FOR_LOOP_(mField,BLOCKSET|MAT_MOISTURESET,it)) {
        if(it->get_Cubit_name().compare(0,18,"MAT_MOISTURE_FIBRE") == 0){
          Mat_Moisture capillary_data;
          ierr = it->get_attribute_data_structure(capillary_data); CHKERRQ(ierr);
          setOfBlocks[it->get_msId()].vIscosity = capillary_data.data.Viscosity;
          setOfBlocks[it->get_msId()].pErmeability = capillary_data.data.Permeability;
          cout<<"capillary_data.data.Viscosity = "<<capillary_data.data.Viscosity<<endl;
          cout<<"capillary_data.data.Permeability = "<<capillary_data.data.Permeability<<endl;
          rval = mField.get_moab().get_entities_by_type(it->meshset,MBTET,setOfBlocks[it->get_msId()].tEts,true); CHKERR_PETSC(rval);
          cout<<"setOfBlocks[it->get_msId()].tEts.size() = "<<setOfBlocks[it->get_msId()].tEts.size()<<endl;
          ierr = mField.add_ents_to_finite_element_by_TETs(setOfBlocks[it->get_msId()].tEts,"CAPILLARY_FE"); CHKERRQ(ierr);
        }
      }
      
      PetscFunctionReturn(0);
    }

    
    /** \brief add heat flux element
     * \infroup mofem_thermal_elem
     *
     * It get data from het flux set and define elemenet in moab. Aletrantively
     * uses block set with name HET_FLUX.
     *
     * \param problem name
     * \param field name
     * \param name of mesh nodal positions (if not defined nodal coordinates are used)
     */
    PetscErrorCode addDiffusionFluxElement(
                                         const string problem_name,const string field_name,const string mesh_nodals_positions = "MESH_NODE_POSITIONS") {
      PetscFunctionBegin;
      
      PetscErrorCode ierr;
      ErrorCode rval;
      
      ierr = mField.add_finite_element("DIFFUSION_FLUX_FE",MF_ZERO); CHKERRQ(ierr);
      ierr = mField.modify_finite_element_add_field_row("DIFFUSION_FLUX_FE",field_name); CHKERRQ(ierr);
      ierr = mField.modify_finite_element_add_field_col("DIFFUSION_FLUX_FE",field_name); CHKERRQ(ierr);
      ierr = mField.modify_finite_element_add_field_data("DIFFUSION_FLUX_FE",field_name); CHKERRQ(ierr);
      if(mField.check_field(mesh_nodals_positions)) {
        ierr = mField.modify_finite_element_add_field_data("DIFFUSION_FLUX_FE",mesh_nodals_positions); CHKERRQ(ierr);
      }
      ierr = mField.modify_problem_add_finite_element(problem_name,"DIFFUSION_FLUX_FE"); CHKERRQ(ierr);
      
      //this is alternative method for setting boundary conditions, to bypass bu in cubit file reader.
      //not elegant, but good enough
      for(_IT_CUBITMESHSETS_BY_SET_TYPE_FOR_LOOP_(mField,BLOCKSET,it)) {
        if(it->get_Cubit_name().compare(0,9,"MASS_FLUX") == 0) {
          vector<double> data;
          ierr = it->get_Cubit_attributes(data); CHKERRQ(ierr);
          if(data.size()!=1) {
            SETERRQ(PETSC_COMM_SELF,1,"Data inconsistency");
          }
//          cout<<"data[0]   "<<data[0]<<endl;
//          std::string wait;
//          std::cin >> wait;
          
          setOfFluxes[it->get_msId()].dAta_mass_flux = data[0];
          //cerr << setOfFluxes[it->get_msId()].dAta << endl;
          rval = mField.get_moab().get_entities_by_type(it->meshset,MBTRI,setOfFluxes[it->get_msId()].tRis,true); CHKERR_PETSC(rval);
          ierr = mField.add_ents_to_finite_element_by_TRIs(setOfFluxes[it->get_msId()].tRis,"DIFFUSION_FLUX_FE"); CHKERRQ(ierr);
          
        }
      }

      PetscFunctionReturn(0);
    }
   
    
    
    /** \brief this calas is to control time stepping
     * \infroup mofem_thermal_elem
     *
     * It is used to save data for temerature rate vectot to MoFEM field.
     */
    struct UpdateAndControl: public FEMethod {
      
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
                                                     problemPtr,tempName,rateName,ROW,ts_u_t,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
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
        
        ierr = mField.set_global_VecCreateGhost(
                                                problemPtr,ROW,ts_u,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
        
        BitRefLevel proble_bit_level = problemPtr->get_BitRefLevel();
        
        SeriesRecorder *recorder_ptr;
        ierr = mField.query_interface(recorder_ptr); CHKERRQ(ierr);
        ierr = recorder_ptr->record_begin(seriesName); CHKERRQ(ierr);
        ierr = recorder_ptr->record_field(seriesName,tempName,proble_bit_level,mask); CHKERRQ(ierr);
        ierr = recorder_ptr->record_end(seriesName); CHKERRQ(ierr);
        
        PetscFunctionReturn(0);
      }
      
    };

    
    
    
    
    
    
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
    
    
    /** \brief this fucntion is used in case of stationary heat conductivity problem for lhs
     * \infroup mofem_moisture_elem
     */
    PetscErrorCode setDiffusionElementLhsOperators(string field_name,Mat A) {
      PetscFunctionBegin;
      map<int,BlockData>::iterator sit = setOfBlocks.begin();
      for(;sit!=setOfBlocks.end();sit++) {
        //add finite elemen
        feLhsDiffusion.get_op_to_do_Lhs().push_back(new OpDiffusionLhs(field_name,A,sit->second,commonData));
      }
      PetscFunctionReturn(0);
    }
    
    
    /** \biref operator to calculate left hand side of het conductivity terms
     * \infroup mofem_moisture_elem
     */
    struct OpDiffusionLhs: public TetElementForcesAndSourcesCore::UserDataOperator {
      
      BlockData &dAta;
      CommonData &commonData;
      bool useTsB;
      OpDiffusionLhs(const string field_name,BlockData &data,CommonData &common_data):
      TetElementForcesAndSourcesCore::UserDataOperator(field_name),
      dAta(data),commonData(common_data),useTsB(true) { }
      
      Mat A;
      OpDiffusionLhs(const string field_name,Mat _A,BlockData &data,CommonData &common_data):
      TetElementForcesAndSourcesCore::UserDataOperator(field_name),
      dAta(data),commonData(common_data),useTsB(false),A(_A) {}
      
      ublas::matrix<double> K,transK;
      
      /** \brief calculate moisture diffusivity matrix
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
        
//        cout<<"OpDiffusionLhs  start "<<endl;
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
          
//          cout<<"dAta.dIffusivity = "<<dAta.dIffusivity<<endl;

          for(unsigned int gg = 0;gg<row_data.getN().size1();gg++) {
            //            cout<<"gg = "<<gg<<endl;
            //            cout<<"dAta.dIffusivity = "<<dAta.dIffusivity<<endl;
            //            cout<<"getGaussPts()(3,gg) = "<<getGaussPts()(3,gg)<<endl;
            //            cout<<"row_data.getN().size1() = "<<row_data.getN().size1()<<endl;
            double val = dAta.dIffusivity*getVolume()*getGaussPts()(3,gg);
            if(getHoGaussPtsDetJac().size()>0) {
              val *= getHoGaussPtsDetJac()[gg]; ///< higher order geometry
            }
            
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
        
//        cout<<"OpDiffusionLhs  End "<<endl;
        PetscFunctionReturn(0);
      }
    
    };
    
    
    /** \brief this fucntion is used in case of stationary heat conductivity problem for lhs
     * \infroup mofem_moisture_elem
     */
    PetscErrorCode setCapillaryElementLhsOperators(string field_name,Mat A) {
      PetscFunctionBegin;
      map<int,BlockData>::iterator sit = setOfBlocks.begin();
      for(;sit!=setOfBlocks.end();sit++) {
        //add finite elemen
        feLhsCapillary.get_op_to_do_Lhs().push_back(new OpCapillaryLhs(field_name,A,sit->second,commonData));
      }
      PetscFunctionReturn(0);
    }

    /** \biref operator to calculate left hand side of het conductivity terms
     * \infroup mofem_thermal_elem
     */
    struct OpCapillaryLhs: public TetElementForcesAndSourcesCore::UserDataOperator {
      
      BlockData &dAta;
      CommonData &commonData;
      bool useTsB;
      OpCapillaryLhs(const string field_name,BlockData &data,CommonData &common_data):
      TetElementForcesAndSourcesCore::UserDataOperator(field_name),
      dAta(data),commonData(common_data),useTsB(true) { }
      
      Mat A;
      OpCapillaryLhs(const string field_name,Mat _A,BlockData &data,CommonData &common_data):
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

    
    
    /** \biref operator to calculate left hand side of miisture capacity terms
     * \infroup mofem_thermal_elem
     */
    struct OpMoistureCapacityLhs: public TetElementForcesAndSourcesCore::UserDataOperator {
      
      BlockData &dAta;
      CommonData &commonData;
      OpMoistureCapacityLhs(const string field_name,BlockData &data,CommonData &common_data):
      TetElementForcesAndSourcesCore::UserDataOperator(field_name),
      dAta(data),commonData(common_data) {}

      ublas::matrix<double> M,transM;
      
      /** \brief calculate heat capacity matrix
       *
       * M = int N^T  N dOmega
       *
       */
      PetscErrorCode doWork(
                            int row_side,int col_side,
                            EntityType row_type,EntityType col_type,
                            DataForcesAndSurcesCore::EntData &row_data,
                            DataForcesAndSurcesCore::EntData &col_data) {
        PetscFunctionBegin;
//        cout<<"OpMoistureCapacityLhs  start "<<endl;


        try {
          if(row_data.getIndices().size()==0) PetscFunctionReturn(0);
          if(col_data.getIndices().size()==0) PetscFunctionReturn(0);
          
          int nb_row = row_data.getN().size2();
          int nb_col = col_data.getN().size2();
          
//          cout<<"nb_row  =  "<<nb_row<<endl;
//          cout<<"nb_col  =  "<<nb_col<<endl;
//          //  std::string wait;
//          //  std::cin >> wait;

          M.resize(nb_row,nb_col);
          bzero(&*M.data().begin(),nb_row*nb_col*sizeof(double));

          for(unsigned int gg = 0;gg<row_data.getN().size1();gg++) {
            double val = getVolume()*getGaussPts()(3,gg);
            if(getHoGaussPtsDetJac().size()>0) {
              val *= getHoGaussPtsDetJac()[gg]; ///< higher order geometry
            }
            val *= getFEMethod()->ts_a;

            //cblas
            //double *N_row,*N_col;
            //N_row = &row_data.getN()(gg,0);
            //N_col = &col_data.getN()(gg,0);
            //cblas_dger(CblasRowMajor,
            //  nb_row,nb_col,val,N_row,1,N_col,1,&M(0,0),nb_col);
            
            //ublas
            noalias(M) += val*outer_prod( row_data.getN(gg,nb_row),col_data.getN(gg,nb_col) );
            
          }

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
        
//        cout<<"OpMoistureCapacityLhs  end "<<endl;

        PetscFunctionReturn(0);
      }
      
    };

    
    
    
    
    /** \brief opearator to caulate concentration and rate of concentration at Gauss points
     * \infroup mofem_thermal_elem
     */
    struct OpGetFieldAtGaussPts: public TetElementForcesAndSourcesCore::UserDataOperator {
      
      ublas::vector<double> &fieldAtGaussPts;
      OpGetFieldAtGaussPts(const string field_name,ublas::vector<double> &field_at_gauss_pts):
      TetElementForcesAndSourcesCore::UserDataOperator(field_name),
      fieldAtGaussPts(field_at_gauss_pts) {}
      
      /** \brief operator calulating temererature and rate of temperature
       *
       * temerature temerature or rate of temperature is calculated multiplyingshape functions by degrees of freedom
       */
      PetscErrorCode doWork(
                            int side,EntityType type,DataForcesAndSurcesCore::EntData &data) {
        PetscFunctionBegin;
//        cout<<"OpGetFieldAtGaussPts  start "<<endl;

        try {
          
          if(data.getFieldData().size()==0) PetscFunctionReturn(0);
          int nb_dofs = data.getFieldData().size();
          int nb_gauss_pts = data.getN().size1();
          
          //initialize
          fieldAtGaussPts.resize(nb_gauss_pts);
          if(type == MBVERTEX) {
            //loop over shape functions on entities allways start form
            //vertices, so if nodal shape functions are processed, vector of
            //field values is zeroad at initialization
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
        
//        cout<<"OpGetFieldAtGaussPts  start "<<endl;
        PetscFunctionReturn(0);
      }
      
    };

    
    
    /** \brief operator to calculate tempereature rate at Gauss pts
     * \infroup mofem_thermal_elem
     */
    struct OpGetRateAtGaussPts: public OpGetFieldAtGaussPts {
      OpGetRateAtGaussPts(const string field_name,CommonData &common_data):
      OpGetFieldAtGaussPts(field_name,common_data.temperatureRateAtGaussPts) {}
    };

    
    
    
    /// \brief operator to calulete concentration gradient at Gauss points
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
//        cout<<"OpGetGradAtGaussPts  start "<<endl;

        try {
          
          if(data.getIndices().size()==0) PetscFunctionReturn(0);
          int nb_dofs = data.getFieldData().size();
          int nb_gauss_pts = data.getN().size1();
          
          //initialize
          commonData.gradAtGaussPts.resize(nb_gauss_pts,3);
          if(type == MBVERTEX) {
            fill(commonData.gradAtGaussPts.data().begin(),commonData.gradAtGaussPts.data().end(),0);
          }
          
//          cout<<"nb_gauss_pts ="<< nb_gauss_pts<<endl;
          for(int gg = 0;gg<nb_gauss_pts;gg++) {
            ublas::noalias(commonData.getGradAtGaussPts(gg)) += prod( trans(data.getDiffN(gg,nb_dofs)), data.getFieldData() );
//            cout<<"Hi from OpGetGradAtGaussPts "<<endl;
          }
          
        } catch (const std::exception& ex) {
          ostringstream ss;
          ss << "throw in method: " << ex.what() << endl;
          SETERRQ(PETSC_COMM_SELF,1,ss.str().c_str());
        }
//        cout<<"OpGetGradAtGaussPts  end "<<endl;
        PetscFunctionReturn(0);
      }
      
    };
    

    
    
    
    
    /** \biref operator to calculate right hand side of Moisture conductivity terms
     * \infroup mofem_thermal_elem
     */
    struct OpMoistureRhs: public TetElementForcesAndSourcesCore::UserDataOperator {
      
      BlockData &dAta;
      CommonData &commonData;
      bool useTsF;
      OpMoistureRhs(const string field_name,BlockData &data,CommonData &common_data):
      TetElementForcesAndSourcesCore::UserDataOperator(field_name),
      dAta(data),commonData(common_data),useTsF(true) {}
      
      Vec F;
      OpMoistureRhs(const string field_name,Vec _F,BlockData &data,CommonData &common_data):
      TetElementForcesAndSourcesCore::UserDataOperator(field_name),
      dAta(data),commonData(common_data),useTsF(false),F(_F) { }
      
      ublas::vector<double> Nf;
      
      /** \brief calculate thermal conductivity matrix
       *
       * F = int diffN^T k gard_T dOmega  (internal force)
       *
       */
      PetscErrorCode doWork(
                            int side,EntityType type,DataForcesAndSurcesCore::EntData &data) {
        PetscFunctionBegin;
//        cout<<"OpMoistureRhs  start "<<endl;

        
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
            
            double val = dAta.dIffusivity*getVolume()*getGaussPts()(3,gg);
            
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
        
//        cout<<"OpMoistureRhs  end "<<endl;
        PetscFunctionReturn(0);
      }
      
    };

    
    
    
    /** \biref operator to calculate right hand side of het capacity terms
     * \infroup mofem_thermal_elem
     */
    struct OpMoistureCapacityRhs: public TetElementForcesAndSourcesCore::UserDataOperator {
      
      BlockData &dAta;
      CommonData &commonData;
      OpMoistureCapacityRhs(const string field_name,BlockData &data,CommonData &common_data):
      TetElementForcesAndSourcesCore::UserDataOperator(field_name),
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
//        cout<<"OpMoistureCapacityRhs  start "<<endl;

        try {
          
          if(data.getIndices().size()==0) PetscFunctionReturn(0);
          int nb_row = data.getN().size2();
          Nf.resize(nb_row);
          bzero(&Nf[0],nb_row*sizeof(double));
          
//          cout<<"data.getN().size1() "<<data.getN().size1()<<endl;

          for(unsigned int gg = 0;gg<data.getN().size1();gg++) {
            double val = getVolume()*getGaussPts()(3,gg);
            if(getHoGaussPtsDetJac().size()>0) {
              val *= getHoGaussPtsDetJac()[gg]; ///< higher order geometry
            }
//            cout<<"val "<<val<<endl;
//            cout<<"commonData.temperatureRateAtGaussPts[gg] "<<val<<endl;
            
            val *= commonData.temperatureRateAtGaussPts[gg];
            //cblas
            //cblas_daxpy(nb_row,val,&data.getN()(gg,0),1,&*Nf.data().begin(),1);
            //ublas
            ublas::noalias(Nf) += val*data.getN(gg);
//            cout<<"Nf "<<Nf<<endl;

          }
          PetscErrorCode ierr;
          ierr = VecSetValues(getFEMethod()->ts_F,data.getIndices().size(),
                              &data.getIndices()[0],&Nf[0],ADD_VALUES); CHKERRQ(ierr);
          
          
        } catch (const std::exception& ex) {
          ostringstream ss;
          ss << "throw in method: " << ex.what() << endl;
          SETERRQ(PETSC_COMM_SELF,1,ss.str().c_str());
        }
        
//        cout<<"OpMoistureCapacityRhs  end "<<endl;

        PetscFunctionReturn(0);
      }
      
    };

    
    
    /** \brief operator for calculate heat flux and assemble to right hand side
     * \infroup mofem_thermal_elem
     */
    struct OpMoistureFlux:public TriElementForcesAndSurcesCore::UserDataOperator {
      
      FluxData &dAta;
      bool ho_geometry;
      bool useTsF;
      OpMoistureFlux(const string field_name,FluxData &data,bool _ho_geometry = false):
      TriElementForcesAndSurcesCore::UserDataOperator(field_name),
      dAta(data),ho_geometry(_ho_geometry),useTsF(true) { }

      Vec F;
      OpMoistureFlux(const string field_name,Vec _F,
                 FluxData &data,bool _ho_geometry = false):
      TriElementForcesAndSurcesCore::UserDataOperator(field_name),
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
//        cout<<"OpMoistureFlux start "<<endl;

        if(data.getIndices().size()==0) PetscFunctionReturn(0);
        if(dAta.tRis.find(getMoFEMFEPtr()->get_ent())==dAta.tRis.end()) PetscFunctionReturn(0);
        PetscErrorCode ierr;
        const FENumeredDofMoFEMEntity *dof_ptr;
        ierr = getMoFEMFEPtr()->get_row_dofs_by_petsc_gloabl_dof_idx(data.getIndices()[0],&dof_ptr); CHKERRQ(ierr);
        int rank = dof_ptr->get_max_rank();
        
        int nb_dofs = data.getIndices().size()/rank;
        
        Nf.resize(data.getIndices().size());
        //bzero(&*Nf.data().begin(),data.getIndices().size()*sizeof(FieldData));
        fill(Nf.begin(),Nf.end(),0);
        
        //cerr << getNormal() << endl;
        //cerr << getNormals_at_GaussPt() << endl;
        
        for(unsigned int gg = 0;gg<data.getN().size1();gg++) {
          
          double val = getGaussPts()(2,gg);
          double flux;
          if(ho_geometry) {
            double area = norm_2(getNormals_at_GaussPt(gg)); //cblas_dnrm2(3,&getNormals_at_GaussPt()(gg,0),1);
            flux = dAta.dAta_mass_flux*area;
          } else {
            flux = dAta.dAta_mass_flux*getArea();
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
        
//        cout<<"OpMoistureFlux end "<<endl;

        PetscFunctionReturn(0);
      }
      
    };

    
    
    
    
    
    /** \brief set up operators for unsedy heat flux problem
     * \infroup mofem_thermal_elem
     */
    PetscErrorCode setTimeSteppingProblem(TsCtx &ts_ctx,string field_name,string rate_name,const string mesh_nodals_positions = "MESH_NODE_POSITIONS") {
      PetscFunctionBegin;
      
      {
        map<int,BlockData>::iterator sit = setOfBlocks.begin();
        for(;sit!=setOfBlocks.end();sit++) {
          //add finite element
          feLhsDiffusion.get_op_to_do_Lhs().push_back(new OpDiffusionLhs(field_name,sit->second,commonData));
          feLhsDiffusion.get_op_to_do_Lhs().push_back(new OpMoistureCapacityLhs(field_name,sit->second,commonData));
          
          feRhsDiffusion.get_op_to_do_Rhs().push_back(new OpGetRateAtGaussPts(rate_name,commonData));
          feRhsDiffusion.get_op_to_do_Rhs().push_back(new OpGetGradAtGaussPts(field_name,commonData));
          feRhsDiffusion.get_op_to_do_Rhs().push_back(new OpMoistureRhs(field_name,sit->second,commonData));
          feRhsDiffusion.get_op_to_do_Rhs().push_back(new OpMoistureCapacityRhs(field_name,sit->second,commonData));
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
          feFluxDiffusion.get_op_to_do_Rhs().push_back(new OpMoistureFlux(field_name,sit->second,ho_geometry));
        }
      }
      
      //rhs
      TsCtx::loops_to_do_type& loops_to_do_Rhs = ts_ctx.get_loops_to_do_IFunction();
      loops_to_do_Rhs.push_back(TsCtx::loop_pair_type("DIFFUSION_FE",&feRhsDiffusion));
      loops_to_do_Rhs.push_back(TsCtx::loop_pair_type("DIFFUSION_FLUX_FE",&feFluxDiffusion));
      
      //lhs
      TsCtx::loops_to_do_type& loops_to_do_Mat = ts_ctx.get_loops_to_do_IJacobian();
      loops_to_do_Mat.push_back(TsCtx::loop_pair_type("DIFFUSION_FE",&feLhsDiffusion));
      
      PetscFunctionReturn(0);
    }

    
    
  };
  
}

#endif //__MOISTURE_ELEMENT_HPP

/***************************************************************************//**
                                                                              * \defgroup mofem_moisture_elem Moisture element
                                                                              * \ingroup mofem_forces_and_sources
                                                                              ******************************************************************************/

