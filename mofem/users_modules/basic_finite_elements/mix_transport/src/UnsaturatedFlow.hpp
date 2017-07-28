/** \file UnsaturatedFlow.hpp
 * \brief Mix implementation of transport element
 *
 * \ingroup mofem_mix_transport_elem
 *
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

#ifndef __UNSATURATD_FLOW_HPP__
#define __UNSATURATD_FLOW_HPP__

namespace MixTransport {

  struct GenericMaterial {

    virtual ~GenericMaterial() {
    }

    double h;           ///< hydraulic head

    double K;           ///< Hydraulic conductivity [L/s]
    double diffK;       ///< Derivative of hydraulic conductivity [L/s * L^2/F]
    double C;           ///< Capacity [S^2/L^2]
    double diffC;       ///< Derivative of capacity [S^2/L^2 * L^2/F ]

    Range tEts;         ///< Elements with this material

    double x,y,z;       ///< in meters (L)


    /**
     * \brief Initialize head
     * @return value of head
     */
    virtual double initalPcEval() const = 0;

    virtual PetscErrorCode calK() {
      PetscFunctionBegin;
      SETERRQ(
        PETSC_COMM_SELF,MOFEM_NOT_IMPLEMENTED,
        "Not implemented how to calculate hydraulic conductivity"
      );
      PetscFunctionReturn(0);
    }

    virtual PetscErrorCode calDiffK() {
      PetscFunctionBegin;
      SETERRQ(
        PETSC_COMM_SELF,MOFEM_NOT_IMPLEMENTED,
        "Not implemented how to calculate derivative of hydraulic conductivity"
      );
      PetscFunctionReturn(0);
    }

    virtual PetscErrorCode calC() {
      PetscFunctionBegin;
      SETERRQ(
        PETSC_COMM_SELF,MOFEM_NOT_IMPLEMENTED,
        "Not implemented how to calculate capacity"
      );
      PetscFunctionReturn(0);
    }

    virtual PetscErrorCode calDiffC() {
      PetscFunctionBegin;
      SETERRQ(
        PETSC_COMM_SELF,MOFEM_NOT_IMPLEMENTED,
        "Not implemented how to calculate capacity"
      );
      PetscFunctionReturn(0);
    }

  };

  struct UnsaturatedFlowElement: public MixTransportElement {

    DM dM;

    UnsaturatedFlowElement(MoFEM::Interface &m_field):
    MixTransportElement(m_field),
    dM(PETSC_NULL),
    lastEvalBcValEnt(0),
    lastEvalBcBlockValId(-1),
    lastEvalBcFluxEnt(0),
    lastEvalBcBlockFluxId(-1) {
    }

    ~UnsaturatedFlowElement() {

      if(dM!=PETSC_NULL) {
        ierr = DMDestroy(&dM); CHKERRABORT(PETSC_COMM_WORLD,ierr);
      }
    }

    typedef std::map<int,boost::shared_ptr<GenericMaterial> > MaterialsDoubleMap;
    MaterialsDoubleMap dMatMap; ///< materials

    boost::shared_ptr<VectorDouble> headRateAtGaussPts;

    virtual PetscErrorCode getMaterial(
      const EntityHandle ent,int &block_id
    ) const {
      PetscFunctionBegin;
      for(
        MaterialsDoubleMap::const_iterator mit = dMatMap.begin();
        mit!=dMatMap.end();mit++
      ) {
        if(mit->second->tEts.find(ent)!=mit->second->tEts.end()) {
          block_id = mit->first;
          PetscFunctionReturn(0);
        }
      }
      SETERRQ(mField.get_comm(),MOFEM_DATA_INCONSISTENCY,"Element not found, no material data");
      PetscFunctionReturn(0);
    }

    struct BcData {
      Range eNts;
      double fixValue;
      boost::function<
      double (const double x,const double y,const double z)
      > hookFun;
      BcData():
      hookFun(NULL) {
      }
    };
    typedef map<int,boost::shared_ptr<BcData> > BcMap;
    BcMap bcValueMap;

    EntityHandle lastEvalBcValEnt;
    int lastEvalBcBlockValId;

    PetscErrorCode getBcOnValues(
      const EntityHandle ent,const int gg,
      const double x,const double y,const double z,
      double &value
    ) {
      PetscFunctionBegin;
      int block_id = -1;
      if(lastEvalBcValEnt==ent) {
        block_id = lastEvalBcBlockValId;
      } else {
        for(BcMap::iterator it = bcValueMap.begin();it!=bcValueMap.end();it++) {
          if(it->second->eNts.find(ent)!=it->second->eNts.end()) {
            block_id = it->first;
          }
        }
        lastEvalBcValEnt = ent;
        lastEvalBcBlockValId = block_id;
      }
      if(block_id>=0) {
        if(bcValueMap.at(block_id)->hookFun) {
          value = bcValueMap.at(block_id)->hookFun(x,y,z);
        } else {
          value = bcValueMap.at(block_id)->fixValue;
        }
      } else {
        value = 0;
      }
      PetscFunctionReturn(0);
    }

    BcMap bcFluxMap;
    EntityHandle lastEvalBcFluxEnt;
    int lastEvalBcBlockFluxId;

    /**
     * \brief essential (Neumann) boundary condition (set fluxes)
     * @param  ent  handle to finite element entity
     * @param  x    coord
     * @param  y    coord
     * @param  z    coord
     * @param  flux reference to flux which is set by function
     * @return      [description]
     */
    PetscErrorCode getBcOnFluxes(
      const EntityHandle ent,
      const double x,const double y,const double z,
      double &flux
    ) {
      PetscFunctionBegin;
      int block_id = -1;
      if(lastEvalBcFluxEnt==ent) {
        block_id = lastEvalBcBlockFluxId;
      } else {
        for(BcMap::iterator it = bcFluxMap.begin();it!=bcFluxMap.end();it++) {
          if(it->second->eNts.find(ent)!=it->second->eNts.end()) {
            block_id = it->first;
          }
        }
        lastEvalBcFluxEnt = ent;
        lastEvalBcBlockFluxId = block_id;
      }
      if(block_id>=0) {
        if(bcFluxMap.at(block_id)->hookFun) {
          flux = bcFluxMap.at(block_id)->hookFun(x,y,z);
        } else {
          flux = bcFluxMap.at(block_id)->fixValue;
        }
      } else {
        flux = 0;
      }
      PetscFunctionReturn(0);
    }

    struct OpRhsBcOnValues: public MoFEM::FaceElementForcesAndSourcesCore::UserDataOperator {

      UnsaturatedFlowElement &cTx;
      boost::shared_ptr<MethodForForceScaling> valueScale;

      /**
       * \brief Constructor
       */
      OpRhsBcOnValues(
        UnsaturatedFlowElement &ctx,const std::string fluxes_name,
        boost::shared_ptr<MethodForForceScaling>& value_scale
      ):
      MoFEM::FaceElementForcesAndSourcesCore::UserDataOperator(fluxes_name,UserDataOperator::OPROW),
      cTx(ctx),
      valueScale(value_scale) {
      }

      VectorDouble nF;

      /**
       * \brief Integrate boundary condition
       * @param  side local index of entity
       * @param  type type of entity
       * @param  data data on entity
       * @return      error code
       */
      PetscErrorCode doWork(
        int side,EntityType type,DataForcesAndSurcesCore::EntData &data
      ) {
        PetscFunctionBegin;
        try {
          if(data.getFieldData().size()==0) PetscFunctionReturn(0);
          EntityHandle fe_ent = getNumeredEntFiniteElementPtr()->getEnt();
          nF.resize(data.getIndices().size());
          nF.clear();
          int nb_gauss_pts = data.getHdivN().size1();
          for(int gg = 0;gg<nb_gauss_pts;gg++) {
            double x,y,z;
            x = getCoordsAtGaussPts()(gg,0);
            y = getCoordsAtGaussPts()(gg,1);
            z = getCoordsAtGaussPts()(gg,2);
            double value;
            ierr = cTx.getBcOnValues(fe_ent,gg,x,y,z,value); CHKERRQ(ierr);
            value = -z;
            double w = getGaussPts()(2,gg)*0.5;
            noalias(nF) += w*prod(data.getHdivN(gg),getNormal())*value;
          }
          Vec f = getFEMethod()->ts_F;
          if(valueScale) {
            ierr = valueScale->scaleNf(getFEMethod(),nF); CHKERRQ(ierr);
          }
          ierr = VecSetValues(
            f,data.getIndices().size(),&data.getIndices()[0],&nF[0],ADD_VALUES
          ); CHKERRQ(ierr);
        } catch (const std::exception& ex) {
          std::ostringstream ss;
          ss << "throw in method: " << ex.what() << std::endl;
          SETERRQ(PETSC_COMM_SELF,MOFEM_STD_EXCEPTION_THROW,ss.str().c_str());
        }
        PetscFunctionReturn(0);
      }
    };

    struct OpResidualFlux: public MoFEM::VolumeElementForcesAndSourcesCore::UserDataOperator {

      UnsaturatedFlowElement &cTx;

      OpResidualFlux(UnsaturatedFlowElement &ctx,const std::string& flux_name):
      MoFEM::VolumeElementForcesAndSourcesCore::UserDataOperator(flux_name,UserDataOperator::OPROW),
      cTx(ctx) {}

      VectorDouble divVec,nF;
      FTensor::Index<'i',3> i;

      PetscErrorCode doWork(int side,EntityType type,DataForcesAndSurcesCore::EntData &data) {
        PetscFunctionBegin;
        const int nb_dofs = data.getIndices().size();
        if(nb_dofs==0) PetscFunctionReturn(0);
        nF.resize(nb_dofs,false);
        nF.clear();
        // Get EntityHandle of the finite element
        EntityHandle fe_ent = getNumeredEntFiniteElementPtr()->getEnt();
        // Get material block id
        int block_id;
        ierr = cTx.getMaterial(fe_ent,block_id); CHKERRQ(ierr);
        // Get material block
        boost::shared_ptr<GenericMaterial>& block_data = cTx.dMatMap.at(block_id);
        // Get base function
        FTensor::Tensor1<double*,3> t_n_hdiv = data.getFTensor1HdivN<3>();
        // Get pressure
        FTensor::Tensor0<double*> t_h = getTensor0FormData(cTx.valuesAtGaussPts);
        // Get flux
        FTensor::Tensor1<double*,3> t_flux = getTensor1FormData<3>(cTx.fluxesAtGaussPts);
        // Coords at integration points
        FTensor::Tensor1<double*,3> t_coords = getTensor1CoordsAtGaussPts();
        // Get integration weight
        FTensor::Tensor0<double*> t_w = getFTensor0IntegrationWeight();
        // Get volume
        double vol = getVolume();
        // Get material parameters
        int nb_gauss_pts = data.getHdivN().size1();
        for(int gg = 0;gg!=nb_gauss_pts;gg++) {
          const double alpha = t_w*vol;
          block_data->h = t_h;
          block_data->x = t_coords(0);
          block_data->y = t_coords(1);
          block_data->z = t_coords(2);
          ierr = block_data->calK(); CHKERRQ(ierr);
          const double K = block_data->K;
          const double z = t_coords(2); /// z-coordinate at Gauss pt
          // Get divergence
          ierr = getDivergenceOfHDivBaseFunctions(side,type,data,gg,divVec); CHKERRQ(ierr);
          noalias(nF) -= alpha*(t_h-z)*divVec;
          FTensor::Tensor0<double*> t_nf(&*nF.begin());
          for(int rr = 0;rr!=nb_dofs;rr++) {
            t_nf += alpha*(1/K)*(t_n_hdiv(i)*t_flux(i));
            ++t_n_hdiv;
            ++t_nf;
          }
          ++t_h;
          ++t_flux;
          ++t_coords;
          ++t_w;
        }
        ierr = VecSetValues(
          getFEMethod()->ts_F,nb_dofs,
          &*data.getIndices().begin(),&*nF.begin(),ADD_VALUES
        ); CHKERRQ(ierr);
        PetscFunctionReturn(0);
      }

    };

    struct OpResidualMass: public MoFEM::VolumeElementForcesAndSourcesCore::UserDataOperator {

      UnsaturatedFlowElement &cTx;

      OpResidualMass(UnsaturatedFlowElement &ctx,const std::string& val_name):
      MoFEM::VolumeElementForcesAndSourcesCore::UserDataOperator(val_name,UserDataOperator::OPROW),
      cTx(ctx) {
      }

      VectorDouble nF;

      PetscErrorCode doWork(int side,EntityType type,DataForcesAndSurcesCore::EntData &data) {
        PetscFunctionBegin;
        PetscFunctionBegin;
        const int nb_dofs = data.getIndices().size();
        if(nb_dofs==0) PetscFunctionReturn(0);
        nF.resize(nb_dofs,false);
        nF.clear();
        // Get EntityHandle of the finite element
        EntityHandle fe_ent = getNumeredEntFiniteElementPtr()->getEnt();
        // Get material block id
        int block_id;
        ierr = cTx.getMaterial(fe_ent,block_id); CHKERRQ(ierr);
        // Get material block
        boost::shared_ptr<GenericMaterial>& block_data = cTx.dMatMap.at(block_id);
        // Get pressure
        FTensor::Tensor0<double*> t_h = getTensor0FormData(cTx.valuesAtGaussPts);
        // Get pressure rate
        FTensor::Tensor0<double*> t_h_t = getTensor0FormData(*cTx.headRateAtGaussPts);
        // Flux divergence
        FTensor::Tensor0<double*> t_div_flux = getTensor0FormData(cTx.divergenceAtGaussPts);
        // Get integration weight
        FTensor::Tensor0<double*> t_w = getFTensor0IntegrationWeight();
        // Coords at integration points
        FTensor::Tensor1<double*,3> t_coords = getTensor1CoordsAtGaussPts();
        // Get volume
        const double vol = getVolume();
        // Get number of integration points
        int nb_gauss_pts = data.getN().size1();
        for(int gg = 0;gg!=nb_gauss_pts;gg++) {
          const double alpha = t_w*vol;
          block_data->h = t_h;
          block_data->x = t_coords(0);
          block_data->y = t_coords(1);
          block_data->z = t_coords(2);
          ierr = block_data->calC(); CHKERRQ(ierr);
          const double C = block_data->C;
          noalias(nF) += (alpha*(t_div_flux+C*t_h_t))*data.getN(gg);
          ++t_h;
          ++t_h_t;
          ++t_div_flux;
          ++t_coords;
          ++t_w;
        }
        Vec f = getFEMethod()->ts_F;
        ierr = VecSetValues(
          f,nb_dofs,&*data.getIndices().begin(),&*nF.begin(),ADD_VALUES
        ); CHKERRQ(ierr);
        PetscFunctionReturn(0);
      }

    };

    struct OpTauDotSigma_HdivHdiv: public MoFEM::VolumeElementForcesAndSourcesCore::UserDataOperator {

      UnsaturatedFlowElement &cTx;

      OpTauDotSigma_HdivHdiv(UnsaturatedFlowElement &ctx,const std::string flux_name):
      MoFEM::VolumeElementForcesAndSourcesCore::UserDataOperator(
        flux_name,flux_name,UserDataOperator::OPROWCOL
      ),
      cTx(ctx) {
        sYmm = true;
      }

      MatrixDouble nN,transNN;

      FTensor::Index<'j',3> j;

      /**
       * \brief Assemble matrix
       * @param  row_side local index of row entity on element
       * @param  col_side local index of col entity on element
       * @param  row_type type of row entity, f.e. MBVERTEX, MBEDGE, or MBTET
       * @param  col_type type of col entity, f.e. MBVERTEX, MBEDGE, or MBTET
       * @param  row_data data for row
       * @param  col_data data for col
       * @return          error code
       */
      PetscErrorCode doWork(
        int row_side,int col_side,
        EntityType row_type,EntityType col_type,
        DataForcesAndSurcesCore::EntData &row_data,
        DataForcesAndSurcesCore::EntData &col_data
      ) {
        PetscFunctionBegin;
        try {
          const int nb_row = row_data.getIndices().size();
          const int nb_col = col_data.getIndices().size();
          if(nb_row==0) PetscFunctionReturn(0);
          if(nb_col==0) PetscFunctionReturn(0);
          nN.resize(nb_row,nb_col,false);
          nN.clear();
          // Get EntityHandle of the finite element
          EntityHandle fe_ent = getNumeredEntFiniteElementPtr()->getEnt();
          // Get material block id
          int block_id;
          ierr = cTx.getMaterial(fe_ent,block_id); CHKERRQ(ierr);
          // Get material block
          boost::shared_ptr<GenericMaterial>& block_data = cTx.dMatMap.at(block_id);
          // Get pressure
          FTensor::Tensor0<double*> t_h = getTensor0FormData(cTx.valuesAtGaussPts);
          // Coords at integration points
          FTensor::Tensor1<double*,3> t_coords = getTensor1CoordsAtGaussPts();
          // Get base functions
          FTensor::Tensor1<double*,3> t_n_hdiv_row = row_data.getFTensor1HdivN<3>();
          // Get integration weight
          FTensor::Tensor0<double*> t_w = getFTensor0IntegrationWeight();
          // Get volume
          const double vol = getVolume();
          int nb_gauss_pts = row_data.getHdivN().size1();
          for(int gg = 0;gg!=nb_gauss_pts;gg++) {
            block_data->h = t_h;
            block_data->x = t_coords(0);
            block_data->y = t_coords(1);
            block_data->z = t_coords(2);
            ierr = block_data->calK(); CHKERRQ(ierr);
            const double K = block_data->K;
            // get integration weight and multiply by element volume
            double alpha = t_w*vol;
            FTensor::Tensor0<double*> t_a(&*nN.data().begin());
            for(int kk = 0;kk!=nb_row;kk++) {
              FTensor::Tensor1<double*,3> t_n_hdiv_col = col_data.getFTensor1HdivN<3>(gg,0);
              for(int ll = 0;ll!=nb_col;ll++) {
                t_a += alpha*(1/K)*(t_n_hdiv_row(j)*t_n_hdiv_col(j));
                ++t_n_hdiv_col;
                ++t_a;
              }
              ++t_n_hdiv_row;
            }
            ++t_coords;
            ++t_h;
            ++t_w;
          }
          Mat a = getFEMethod()->ts_B;
          ierr = MatSetValues(
            a,
            nb_row,&*row_data.getIndices().begin(),
            nb_col,&*col_data.getIndices().begin(),
            &*nN.data().begin(),ADD_VALUES
          ); CHKERRQ(ierr);
          // matrix is symmetric, assemble other part
          if(row_side != col_side || row_type != col_type) {
            transNN.resize(nb_col,nb_row);
            noalias(transNN) = trans(nN);
            ierr = MatSetValues(
              a,
              nb_col,&*col_data.getIndices().begin(),
              nb_row,&*row_data.getIndices().begin(),
              &*transNN.data().begin(),ADD_VALUES
            ); CHKERRQ(ierr);
          }
        } catch (const std::exception& ex) {
          std::ostringstream ss;
          ss << "throw in method: " << ex.what() << std::endl;
          SETERRQ(PETSC_COMM_SELF,MOFEM_STD_EXCEPTION_THROW,ss.str().c_str());
        }

        PetscFunctionReturn(0);
      }

    };

    struct OpVU_L2L2: public MoFEM::VolumeElementForcesAndSourcesCore::UserDataOperator {

      UnsaturatedFlowElement &cTx;

      OpVU_L2L2(UnsaturatedFlowElement &ctx,const std::string value_name):
      MoFEM::VolumeElementForcesAndSourcesCore::UserDataOperator(
        value_name,value_name,UserDataOperator::OPROWCOL
      ),
      cTx(ctx) {
        sYmm = true;
      }

      MatrixDouble nN;

      /**
       * \brief Assemble matrix
       * @param  row_side local index of row entity on element
       * @param  col_side local index of col entity on element
       * @param  row_type type of row entity, f.e. MBVERTEX, MBEDGE, or MBTET
       * @param  col_type type of col entity, f.e. MBVERTEX, MBEDGE, or MBTET
       * @param  row_data data for row
       * @param  col_data data for col
       * @return          error code
       */
      PetscErrorCode doWork(
        int row_side,int col_side,
        EntityType row_type,EntityType col_type,
        DataForcesAndSurcesCore::EntData &row_data,
        DataForcesAndSurcesCore::EntData &col_data
      ) {
        PetscFunctionBegin;
        try {
          int nb_row = row_data.getIndices().size();
          int nb_col = col_data.getIndices().size();
          if(nb_row==0) PetscFunctionReturn(0);
          if(nb_col==0) PetscFunctionReturn(0);
          nN.resize(nb_row,nb_col,false);
          nN.clear();
          // Get EntityHandle of the finite element
          EntityHandle fe_ent = getNumeredEntFiniteElementPtr()->getEnt();
          // Get material block id
          int block_id;
          ierr = cTx.getMaterial(fe_ent,block_id); CHKERRQ(ierr);
          // Get material block
          boost::shared_ptr<GenericMaterial>& block_data = cTx.dMatMap.at(block_id);
          // Get pressure
          FTensor::Tensor0<double*> t_h = getTensor0FormData(cTx.valuesAtGaussPts);
          // Get pressure rate
          FTensor::Tensor0<double*> t_h_t = getTensor0FormData(*cTx.headRateAtGaussPts);
          // Get integration weight
          FTensor::Tensor0<double*> t_w = getFTensor0IntegrationWeight();
          // Coords at integration points
          FTensor::Tensor1<double*,3> t_coords = getTensor1CoordsAtGaussPts();
          // Time step factor
          double ts_a = getFEMethod()->ts_a;
          // Get volume
          const double vol = getVolume();
          int nb_gauss_pts = row_data.getN().size1();
          // Get base functions
          FTensor::Tensor0<double*> t_n_row = row_data.getFTensor0N();
          for(int gg = 0;gg!=nb_gauss_pts;gg++) {
            // get integration weight and multiply by element volume
            double alpha = t_w*vol;
            block_data->h = t_h;
            block_data->x = t_coords(0);
            block_data->y = t_coords(1);
            block_data->z = t_coords(2);
            ierr = block_data->calC(); CHKERRQ(ierr);
            ierr = block_data->calDiffC(); CHKERRQ(ierr);
            const double C = block_data->C;
            const double diffC = block_data->diffC;
            FTensor::Tensor0<double*> t_a(&*nN.data().begin());
            for(int kk = 0;kk!=nb_row;kk++) {
              FTensor::Tensor0<double*> t_n_col = col_data.getFTensor0N(gg,0);
              for(int ll = 0;ll!=nb_col;ll++) {
                t_a += (alpha*(C*ts_a+diffC*t_h_t))*t_n_row*t_n_col;
                ++t_n_col;
                ++t_a;
              }
              ++t_n_row;
            }
            ++t_w;
            ++t_coords;
            ++t_h;
            ++t_h_t;
          }
          Mat a = getFEMethod()->ts_B;
          ierr = MatSetValues(
            a,
            nb_row,&row_data.getIndices()[0],
            nb_col,&col_data.getIndices()[0],
            &*nN.data().begin(),ADD_VALUES
          ); CHKERRQ(ierr);
        } catch (const std::exception& ex) {
          std::ostringstream ss;
          ss << "throw in method: " << ex.what() << std::endl;
          SETERRQ(PETSC_COMM_SELF,MOFEM_STD_EXCEPTION_THROW,ss.str().c_str());
        }

        PetscFunctionReturn(0);
      }

    };

    struct OpVDivSigma_L2Hdiv: public MoFEM::VolumeElementForcesAndSourcesCore::UserDataOperator {

      UnsaturatedFlowElement &cTx;

      /**
       * \brief Constructor
       */
      OpVDivSigma_L2Hdiv(
        UnsaturatedFlowElement &ctx,
        const std::string& val_name_row,
        const std::string& flux_name_col
      ):
      MoFEM::VolumeElementForcesAndSourcesCore::UserDataOperator(
        val_name_row,flux_name_col,UserDataOperator::OPROWCOL,false
      ),
      cTx(ctx) {
      }

      MatrixDouble nN;
      VectorDouble divVec;

      /**
       * \brief Do calculations
       * @param  row_side local index of entity on row
       * @param  col_side local index of entity on column
       * @param  row_type type of row entity
       * @param  col_type type of col entity
       * @param  row_data row data structure carrying information about base functions, DOFs indices, etc.
       * @param  col_data column data structure carrying information about base functions, DOFs indices, etc.
       * @return          error code
       */
      PetscErrorCode doWork(
        int row_side,int col_side,
        EntityType row_type,EntityType col_type,
        DataForcesAndSurcesCore::EntData &row_data,
        DataForcesAndSurcesCore::EntData &col_data
      ) {
        PetscFunctionBegin;
        try {
          int nb_row = row_data.getFieldData().size();
          int nb_col = col_data.getFieldData().size();
          if(nb_row==0) PetscFunctionReturn(0);
          if(nb_col==0) PetscFunctionReturn(0);
          nN.resize(nb_row,nb_col,false);
          divVec.resize(nb_col,false);
          nN.clear();
          int nb_gauss_pts = row_data.getHdivN().size1();
          for(int gg = 0;gg<nb_gauss_pts;gg++) {
            double alpha = getGaussPts()(3,gg)*getVolume();
            ierr = getDivergenceOfHDivBaseFunctions(
              col_side,col_type,col_data,gg,divVec
            ); CHKERRQ(ierr);
            noalias(nN) += alpha*outer_prod(row_data.getN(gg),divVec);
          }
          ierr = MatSetValues(
            getFEMethod()->ts_B,
            nb_row,&row_data.getIndices()[0],
            nb_col,&col_data.getIndices()[0],
            &nN(0,0),ADD_VALUES
          ); CHKERRQ(ierr);
        } catch (const std::exception& ex) {
          std::ostringstream ss;
          ss << "throw in method: " << ex.what() << std::endl;
          SETERRQ(PETSC_COMM_SELF,MOFEM_STD_EXCEPTION_THROW,ss.str().c_str());
        }
        PetscFunctionReturn(0);
      }

    };

    struct OpDivTauU_HdivL2: public MoFEM::VolumeElementForcesAndSourcesCore::UserDataOperator {

      UnsaturatedFlowElement &cTx;

      /**
       * \brief Constructor
       */
      OpDivTauU_HdivL2(
        UnsaturatedFlowElement &ctx,
        const std::string& flux_name_col,
        const std::string& val_name_row
      ):
      MoFEM::VolumeElementForcesAndSourcesCore::UserDataOperator(
        flux_name_col,val_name_row,UserDataOperator::OPROWCOL,false
      ),
      cTx(ctx) {
      }

      MatrixDouble nN;
      VectorDouble divVec;
      FTensor::Index<'i',3> i;

      /**
       * \brief Do calculations
       * @param  row_side local index of entity on row
       * @param  col_side local index of entity on column
       * @param  row_type type of row entity
       * @param  col_type type of col entity
       * @param  row_data row data structure carrying information about base functions, DOFs indices, etc.
       * @param  col_data column data structure carrying information about base functions, DOFs indices, etc.
       * @return          error code
       */
      PetscErrorCode doWork(
        int row_side,int col_side,
        EntityType row_type,EntityType col_type,
        DataForcesAndSurcesCore::EntData &row_data,
        DataForcesAndSurcesCore::EntData &col_data
      ) {
        PetscFunctionBegin;
        try {
          int nb_row = row_data.getFieldData().size();
          int nb_col = col_data.getFieldData().size();
          if(nb_row==0) PetscFunctionReturn(0);
          if(nb_col==0) PetscFunctionReturn(0);
          nN.resize(nb_row,nb_col,false);
          divVec.resize(nb_row,false);
          nN.clear();
          // Get EntityHandle of the finite element
          EntityHandle fe_ent = getNumeredEntFiniteElementPtr()->getEnt();
          // Get material block id
          int block_id;
          ierr = cTx.getMaterial(fe_ent,block_id); CHKERRQ(ierr);
          // Get material block
          boost::shared_ptr<GenericMaterial>& block_data = cTx.dMatMap.at(block_id);
          // Get pressure
          FTensor::Tensor0<double*> t_h = getTensor0FormData(cTx.valuesAtGaussPts);
          // Get flux
          FTensor::Tensor1<double*,3> t_flux = getTensor1FormData<3>(cTx.fluxesAtGaussPts);
          // Coords at integration points
          FTensor::Tensor1<double*,3> t_coords = getTensor1CoordsAtGaussPts();
          // Get integration weight
          FTensor::Tensor0<double*> t_w = getFTensor0IntegrationWeight();
          // Get base function
          FTensor::Tensor1<double*,3> t_n_hdiv_row = row_data.getFTensor1HdivN<3>();
          // Get volume
          double vol = getVolume();
          int nb_gauss_pts = row_data.getHdivN().size1();
          for(int gg = 0;gg<nb_gauss_pts;gg++) {
            block_data->h = t_h;
            block_data->x = t_coords(0);
            block_data->y = t_coords(1);
            block_data->z = t_coords(2);
            ierr = block_data->calK(); CHKERRQ(ierr);
            ierr = block_data->calDiffK(); CHKERRQ(ierr);
            const double K = block_data->K;
            const double KK = K*K;
            const double diffK = block_data->diffK;
            double alpha = t_w*vol;
            ierr = getDivergenceOfHDivBaseFunctions(
              row_side,row_type,row_data,gg,divVec
            ); CHKERRQ(ierr);
            noalias(nN) -= alpha*outer_prod(divVec,col_data.getN(gg));
            FTensor::Tensor0<double*> t_a(&*nN.data().begin());
            for(int rr = 0;rr!=nb_row;rr++) {
              double beta = alpha*(-diffK/KK)*(t_n_hdiv_row(i)*t_flux(i));
              FTensor::Tensor0<double*> t_n_col = col_data.getFTensor0N(gg,0);
              for(int cc =0;cc!=nb_col;cc++) {
                t_a += beta*t_n_col;
                ++t_n_col;
                ++t_a;
              }
              ++t_n_hdiv_row;
            }
            ++t_w;
            ++t_coords;
            ++t_h;
            ++t_flux;
          }
          ierr = MatSetValues(
            getFEMethod()->ts_B,
            nb_row,&row_data.getIndices()[0],
            nb_col,&col_data.getIndices()[0],
            &nN(0,0),ADD_VALUES
          ); CHKERRQ(ierr);
        } catch (const std::exception& ex) {
          std::ostringstream ss;
          ss << "throw in method: " << ex.what() << std::endl;
          SETERRQ(PETSC_COMM_SELF,MOFEM_STD_EXCEPTION_THROW,ss.str().c_str());
        }
        PetscFunctionReturn(0);
      }

    };

    struct OpEvaluateInitiallHead: public MoFEM::VolumeElementForcesAndSourcesCore::UserDataOperator {
      UnsaturatedFlowElement &cTx;
      OpEvaluateInitiallHead(UnsaturatedFlowElement &ctx,const std::string& val_name):
      MoFEM::VolumeElementForcesAndSourcesCore::UserDataOperator(val_name,UserDataOperator::OPROW),
      cTx(ctx) {
      }

      MatrixDouble nN;
      VectorDouble nF;

      PetscErrorCode doWork(int side,EntityType type,DataForcesAndSurcesCore::EntData &data) {
        PetscFunctionBegin;
        try {
          if(data.getFieldData().size()==0) PetscFunctionReturn(0);
          int nb_dofs = data.getFieldData().size();
          int nb_gauss_pts = data.getN().size1();
          if(nb_dofs!=data.getN().size2()) {
            SETERRQ(PETSC_COMM_WORLD,MOFEM_DATA_INCONSISTENCY,"wrong number of dofs");
          }
          nN.resize(nb_dofs,nb_dofs,false);
          nF.resize(nb_dofs,false);
          nN.clear();
          nF.clear();

          // Get EntityHandle of the finite element
          EntityHandle fe_ent = getFEEntityHandle();
          // Get material block id
          int block_id;
          ierr = cTx.getMaterial(fe_ent,block_id); CHKERRQ(ierr);
          // Get material block
          boost::shared_ptr<GenericMaterial>& block_data = cTx.dMatMap.at(block_id);

          // loop over integration points
          for(int gg = 0;gg<nb_gauss_pts;gg++) {
            // get coordinates at integration point
            block_data->x = getCoordsAtGaussPts()(gg,0);
            block_data->y = getCoordsAtGaussPts()(gg,1);
            block_data->z = getCoordsAtGaussPts()(gg,2);
            // get weight for integration rule
            double alpha = getGaussPts()(2,gg)*getVolume();
            nN += alpha*outer_prod(data.getN(gg),data.getN(gg));
            nF += alpha*block_data->initalPcEval()*data.getN(gg);
          }

          // factor matrix
          cholesky_decompose(nN);
          // solve local problem
          cholesky_solve(nN,nF,ublas::lower());

          // set solution to vector
          ierr = VecSetValues(
            cTx.D1,nb_dofs,&*data.getIndices().begin(),
            &*nF.begin(),INSERT_VALUES
          ); CHKERRQ(ierr);

        } catch (const std::exception& ex) {
          std::ostringstream ss;
          ss << "throw in method: " << ex.what() << std::endl;
          SETERRQ(PETSC_COMM_SELF,1,ss.str().c_str());
        }

        PetscFunctionReturn(0);
      }

    };

    struct MonitorPostProc: public FEMethod {

        UnsaturatedFlowElement &cTx;
        boost::shared_ptr<PostProcVolumeOnRefinedMesh> postProc;
        const int fRequency;

        MonitorPostProc(
          UnsaturatedFlowElement &ctx,
          boost::shared_ptr<PostProcVolumeOnRefinedMesh>& post_proc,
          const int frequency
        ):
        cTx(ctx),
        postProc(post_proc),
        fRequency(frequency) {
        }

        PetscErrorCode preProcess() {
          PetscFunctionBegin;
          PetscFunctionReturn(0);
        }

        PetscErrorCode operator()() {
          PetscFunctionBegin;
          PetscFunctionReturn(0);
        }

        PetscErrorCode postProcess() {
          PetscFunctionBegin;
          int step;
          ierr = TSGetTimeStepNumber(ts,&step); CHKERRQ(ierr);
          if((step)%fRequency==0) {
            ierr = DMoFEMLoopFiniteElements(cTx.dM,"MIX",postProc); CHKERRQ(ierr);
            ierr = postProc->writeFile(
              string("out_")+boost::lexical_cast<std::string>(step)+".h5m"
            ); CHKERRQ(ierr);
          }
          PetscFunctionReturn(0);
        }

      };


    /// \brief add fields
    PetscErrorCode addFields(const std::string &values,const std::string &fluxes,const int order) {

      PetscFunctionBegin;
      //Fields
      ierr = mField.add_field(fluxes,HDIV,DEMKOWICZ_JACOBI_BASE,1); CHKERRQ(ierr);
      ierr = mField.add_field(values,L2,AINSWORTH_LEGENDRE_BASE,1); CHKERRQ(ierr);
      ierr = mField.add_field(values+"_t",L2,AINSWORTH_LEGENDRE_BASE,1); CHKERRQ(ierr);

      //meshset consisting all entities in mesh
      EntityHandle root_set = mField.get_moab().get_root_set();
      //add entities to field

      for(_IT_CUBITMESHSETS_BY_NAME_FOR_LOOP_(mField,"SOIL",it)) {
        ierr = mField.add_ents_to_field_by_type(
          dMatMap[it->getMeshsetId()]->tEts,MBTET,fluxes
        ); CHKERRQ(ierr);
        ierr = mField.add_ents_to_field_by_type(
          dMatMap[it->getMeshsetId()]->tEts,MBTET,values
        ); CHKERRQ(ierr);
        ierr = mField.add_ents_to_field_by_type(
          dMatMap[it->getMeshsetId()]->tEts,MBTET,values+"_t"
        ); CHKERRQ(ierr);
      }

      ierr = mField.set_field_order(root_set,MBTET,fluxes,order+1); CHKERRQ(ierr);
      ierr = mField.set_field_order(root_set,MBTRI,fluxes,order+1); CHKERRQ(ierr);
      ierr = mField.set_field_order(root_set,MBTET,values,order); CHKERRQ(ierr);
      ierr = mField.set_field_order(root_set,MBTET,values+"_t",order); CHKERRQ(ierr);

      PetscFunctionReturn(0);
    }

    /// \brief add finite elements
    PetscErrorCode addFiniteElements(
      const std::string &fluxes_name,const std::string &values_name
    ) {
      PetscFunctionBegin;

      // Define element "MIX". Note that this element will work with fluxes_name and
      // values_name. This reflect bilinear form for the problem
      ierr = mField.add_finite_element("MIX",MF_ZERO); CHKERRQ(ierr);
      ierr = mField.modify_finite_element_add_field_row("MIX",fluxes_name); CHKERRQ(ierr);
      ierr = mField.modify_finite_element_add_field_col("MIX",fluxes_name); CHKERRQ(ierr);
      ierr = mField.modify_finite_element_add_field_row("MIX",values_name); CHKERRQ(ierr);
      ierr = mField.modify_finite_element_add_field_col("MIX",values_name); CHKERRQ(ierr);
      ierr = mField.modify_finite_element_add_field_data("MIX",fluxes_name); CHKERRQ(ierr);
      ierr = mField.modify_finite_element_add_field_data("MIX",values_name); CHKERRQ(ierr);
      ierr = mField.modify_finite_element_add_field_data("MIX",values_name+"_t"); CHKERRQ(ierr);

      for(_IT_CUBITMESHSETS_BY_NAME_FOR_LOOP_(mField,"SOIL",it)) {
        ierr = mField.add_ents_to_finite_element_by_type(
          dMatMap[it->getMeshsetId()]->tEts,MBTET,"MIX"
        ); CHKERRQ(ierr);
      }

      // Define element to integrate natural boundary conditions, i.e. set values.
      ierr = mField.add_finite_element("MIX_BCVALUE",MF_ZERO); CHKERRQ(ierr);
      ierr = mField.modify_finite_element_add_field_row("MIX_BCVALUE",fluxes_name); CHKERRQ(ierr);
      ierr = mField.modify_finite_element_add_field_col("MIX_BCVALUE",fluxes_name); CHKERRQ(ierr);
      ierr = mField.modify_finite_element_add_field_data("MIX_BCVALUE",fluxes_name); CHKERRQ(ierr);
      ierr = mField.modify_finite_element_add_field_data("MIX_BCVALUE",values_name); CHKERRQ(ierr);

      for(_IT_CUBITMESHSETS_BY_NAME_FOR_LOOP_(mField,"CAPILLARY_PRESSURE",it)) {
        ierr = mField.add_ents_to_finite_element_by_type(
          bcValueMap[it->getMeshsetId()]->eNts,MBTRI,"MIX_BCVALUE"
        ); CHKERRQ(ierr);
      }

      // Define element to apply essential boundary conditions.
      ierr = mField.add_finite_element("MIX_BCFLUX",MF_ZERO); CHKERRQ(ierr);
      ierr = mField.modify_finite_element_add_field_row("MIX_BCFLUX",fluxes_name); CHKERRQ(ierr);
      ierr = mField.modify_finite_element_add_field_col("MIX_BCFLUX",fluxes_name); CHKERRQ(ierr);
      ierr = mField.modify_finite_element_add_field_data("MIX_BCFLUX",fluxes_name); CHKERRQ(ierr);
      ierr = mField.modify_finite_element_add_field_data("MIX_BCFLUX",values_name); CHKERRQ(ierr);

      for(_IT_CUBITMESHSETS_BY_NAME_FOR_LOOP_(mField,"WATER_FLUX",it)) {
        ierr = mField.add_ents_to_finite_element_by_type(
          bcFluxMap[it->getMeshsetId()]->eNts,MBTRI,"MIX_BCFLUX"
        ); CHKERRQ(ierr);
      }

      PetscFunctionReturn(0);
    }

    /**
     * \brief Build problem
     * @param  ref_level mesh refinement on which mesh problem you like to built.
     * @return           error code
     */
    PetscErrorCode buildProblem(BitRefLevel ref_level = BitRefLevel().set(0)) {
      PetscFunctionBegin;

      // Build fields
      ierr = mField.build_fields(); CHKERRQ(ierr);
      // Build finite elements
      ierr = mField.build_finite_elements("MIX"); CHKERRQ(ierr);
      ierr = mField.build_finite_elements("MIX_BCFLUX"); CHKERRQ(ierr);
      ierr = mField.build_finite_elements("MIX_BCVALUE"); CHKERRQ(ierr);
      //Build adjacencies of degrees of freedom and elements
      ierr = mField.build_adjacencies(ref_level); CHKERRQ(ierr);


      ierr = DMCreate(PETSC_COMM_WORLD,&dM);CHKERRQ(ierr);
      ierr = DMSetType(dM,"DMMOFEM");CHKERRQ(ierr);
      ierr = DMMoFEMSetIsPartitioned(dM,PETSC_TRUE);
      //set dM data structure which created mofem data structures
      ierr = DMMoFEMCreateMoFEM(dM,&mField,"MIX",ref_level); CHKERRQ(ierr);
      ierr = DMMoFEMSetSquareProblem(dM,PETSC_TRUE); CHKERRQ(ierr);
      ierr = DMMoFEMSetIsPartitioned(dM,PETSC_TRUE); CHKERRQ(ierr);
      ierr = DMSetFromOptions(dM); CHKERRQ(ierr);
      ierr = DMMoFEMAddElement(dM,"MIX"); CHKERRQ(ierr);
      ierr = DMMoFEMAddElement(dM,"MIX_BCFLUX"); CHKERRQ(ierr);
      ierr = DMMoFEMAddElement(dM,"MIX_BCVALUE"); CHKERRQ(ierr);
      ierr = DMSetUp(dM); CHKERRQ(ierr);

      PetscFunctionReturn(0);
    }

    boost::shared_ptr<ForcesAndSurcesCore> feFaceBc;
    boost::shared_ptr<ForcesAndSurcesCore> feFaceRhs;
    boost::shared_ptr<ForcesAndSurcesCore> feVolInitialPc;
    boost::shared_ptr<ForcesAndSurcesCore> feVolRhs;
    boost::shared_ptr<ForcesAndSurcesCore> feVolLhs;
    boost::shared_ptr<MethodForForceScaling> scaleMethodFlux;
    boost::shared_ptr<MethodForForceScaling> scaleMethodValue;
    boost::shared_ptr<FEMethod> tsMonitor;

    /**
    * \brief Set integration rule to volume elements
    *
    */
    struct VolRule {
      int operator()(int,int,int p_data) const {
        return 2*p_data;
      }
    };
    /**
    * \brief Set integration rule to boundary elements
    *
    */
    struct FaceRule {
      int operator()(int p_row,int p_col,int p_data) const {
        return 2*p_data;
      }
    };

    std::vector<int> bcVecIds;
    VectorDouble bcVecVals,vecValsOnBc;

    struct preProcessVol {
      UnsaturatedFlowElement& cTx;
      boost::shared_ptr<ForcesAndSurcesCore> fePtr;
      //std::string mArk;

      preProcessVol(
        UnsaturatedFlowElement& ctx,
        boost::shared_ptr<ForcesAndSurcesCore>& fe_ptr/*,std::string mark*/
      ):
      cTx(ctx),
      fePtr(fe_ptr)/*,mArk(mark)*/ {
      }
      PetscErrorCode operator()() {
        PetscFunctionBegin;
        // Update pressure rates
        ierr = fePtr->mField.set_other_local_ghost_vector(
          fePtr->problemPtr,"VALUES",string("VALUES")+"_t",
          ROW,fePtr->ts_u_t,INSERT_VALUES,SCATTER_REVERSE
        ); CHKERRQ(ierr);
        switch (fePtr->ts_ctx) {
          case TSMethod::CTX_TSSETIFUNCTION:
          if(!cTx.bcIndices.empty()) {
            double scale;
            ierr = cTx.scaleMethodFlux->getForceScale(fePtr->ts_t,scale); CHKERRQ(ierr);
            if(cTx.bcVecIds.size()!=cTx.bcIndices.size()) {
              cTx.bcVecIds.insert(cTx.bcVecIds.begin(),cTx.bcIndices.begin(),cTx.bcIndices.end());
              cTx.bcVecVals.resize(cTx.bcVecIds.size(),false);
              cTx.vecValsOnBc.resize(cTx.bcVecIds.size(),false);
            }
            ierr = VecGetValues(
              cTx.D0,cTx.bcVecIds.size(),&*cTx.bcVecIds.begin(),&*cTx.bcVecVals.begin()
            ); CHKERRQ(ierr);
            ierr = VecGetValues(
              fePtr->ts_u,cTx.bcVecIds.size(),&*cTx.bcVecIds.begin(),&*cTx.vecValsOnBc.begin()
            ); CHKERRQ(ierr);
            cTx.bcVecVals *= scale;
            // cerr << mArk << endl;
            // cerr << "v " << cTx.vecValsOnBc << endl;
            // cerr << "v " << cTx.bcVecVals << endl;
            VectorDouble::iterator vit = cTx.bcVecVals.begin();
            const NumeredDofEntity *dof_ptr;
            for(
              std::vector<int>::iterator it = cTx.bcVecIds.begin();
              it!=cTx.bcVecIds.end();it++,vit++
            ) {
              ierr = fePtr->problemPtr->getColDofsByPetscGlobalDofIdx(*it,&dof_ptr); CHKERRQ(ierr);
              dof_ptr->getFieldData() = *vit;
            }
          } else {
            cTx.bcVecIds.resize(0);
            cTx.bcVecVals.resize(0);
            cTx.vecValsOnBc.resize(0);
          }
          break;
          default:
          // don nothing
          break;
        }
        PetscFunctionReturn(0);
      }
    };

    struct postProcessVol {
      UnsaturatedFlowElement& cTx;
      boost::shared_ptr<ForcesAndSurcesCore> fePtr;
      // std::string mArk;
      postProcessVol(
        UnsaturatedFlowElement& ctx,
        boost::shared_ptr<ForcesAndSurcesCore>& fe_ptr//,std::string mark
      ):
      cTx(ctx),
      fePtr(fe_ptr)/*,mArk(mark)*/ {
      }
      PetscErrorCode operator()() {
        PetscFunctionBegin;
        switch (fePtr->ts_ctx) {
          case TSMethod::CTX_TSSETIJACOBIAN: {
            ierr = MatAssemblyBegin(fePtr->ts_B,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
            ierr = MatAssemblyEnd(fePtr->ts_B,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
            ierr = MatZeroRowsColumns(
              fePtr->ts_B,cTx.bcVecIds.size(),&*cTx.bcVecIds.begin(),1,PETSC_NULL,PETSC_NULL
            ); CHKERRQ(ierr);
            ierr = MatAssemblyBegin(fePtr->ts_B,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
            ierr = MatAssemblyEnd(fePtr->ts_B,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
            // MatView(fePtr->ts_B,PETSC_VIEWER_DRAW_WORLD);
            // std::string wait;
            // std::cin >> wait;
          }
          break;
          case TSMethod::CTX_TSSETIFUNCTION: {
            ierr = VecAssemblyBegin(fePtr->ts_F); CHKERRQ(ierr);
            ierr = VecAssemblyEnd(fePtr->ts_F); CHKERRQ(ierr);
            if(!cTx.bcVecIds.empty()) {
              cTx.vecValsOnBc -= cTx.bcVecVals;
              // cerr << mArk << endl;
              // cerr << "a " << cTx.vecValsOnBc << endl;
              // cerr << "a " << cTx.bcVecVals << endl;
              ierr = VecSetValues(
                fePtr->ts_F,cTx.bcVecIds.size(),
                &*cTx.bcVecIds.begin(),
                &*cTx.vecValsOnBc.begin(),INSERT_VALUES
              ); CHKERRQ(ierr);
            }
            ierr = VecAssemblyBegin(fePtr->ts_F); CHKERRQ(ierr);
            ierr = VecAssemblyEnd(fePtr->ts_F); CHKERRQ(ierr);
          }
          break;
          default:
          // don nothing
          break;
        }
        PetscFunctionReturn(0);
      }
    };

    PetscErrorCode setFiniteElements(
      ForcesAndSurcesCore::RuleHookFun vol_rule = VolRule(),
      ForcesAndSurcesCore::RuleHookFun face_rule = FaceRule()
    ) {
      PetscFunctionBegin;
      // create finite element instances
      feFaceBc = boost::shared_ptr<ForcesAndSurcesCore>(new FaceElementForcesAndSourcesCore(mField));
      feFaceRhs = boost::shared_ptr<ForcesAndSurcesCore>(new FaceElementForcesAndSourcesCore(mField));
      feVolInitialPc = boost::shared_ptr<ForcesAndSurcesCore>(new VolumeElementForcesAndSourcesCore(mField));
      feVolLhs = boost::shared_ptr<ForcesAndSurcesCore>(new VolumeElementForcesAndSourcesCore(mField));
      feVolRhs = boost::shared_ptr<ForcesAndSurcesCore>(new VolumeElementForcesAndSourcesCore(mField));
      // set integration rule to elements
      feFaceBc->getRuleHook = face_rule;
      feFaceRhs->getRuleHook = face_rule;
      feVolInitialPc->getRuleHook = vol_rule;
      feVolLhs->getRuleHook = vol_rule;
      feVolRhs->getRuleHook = vol_rule;
      feVolRhs->preProcessHook = preProcessVol(*this,feVolRhs);
      feVolLhs->preProcessHook = preProcessVol(*this,feVolLhs);
      feVolRhs->postProcessHook = postProcessVol(*this,feVolRhs);
      feVolLhs->postProcessHook = postProcessVol(*this,feVolLhs);

      scaleMethodFlux= boost::shared_ptr<MethodForForceScaling>(
        new TimeForceScale("-flux_history",false)
      );

      scaleMethodValue = boost::shared_ptr<MethodForForceScaling>(
        new TimeForceScale("-pressure_history",false)
      );

      // Set operator to calculate essential boundary conditions
      feFaceBc->getOpPtrVector().push_back(new OpEvaluateBcOnFluxes(*this,"FLUXES"));

      // Set operator to calculate initial capillary pressure
      feVolInitialPc->getOpPtrVector().push_back(new OpEvaluateInitiallHead(*this,"VALUES"));

      // set residual face from Neumann terms, i.e. applied pressure
      feFaceRhs->getOpPtrVector().push_back(new OpRhsBcOnValues(*this,"FLUXES",scaleMethodValue));
      // set residual
      headRateAtGaussPts = boost::make_shared<VectorDouble>();
      feVolRhs->getOpPtrVector().push_back(
        new OpCalculateScalarFieldValues(string("VALUES")+"_t",headRateAtGaussPts,MBTET)
      );
      feVolRhs->getOpPtrVector().push_back(new OpValuesAtGaussPts(*this,"VALUES"));
      feVolRhs->getOpPtrVector().push_back(new OpFluxDivergenceAtGaussPts(*this,"FLUXES"));
      feVolRhs->getOpPtrVector().push_back(new OpResidualFlux(*this,"FLUXES"));
      feVolRhs->getOpPtrVector().push_back(new OpResidualMass(*this,"VALUES"));
      feVolRhs->getOpPtrVector().back().opType = ForcesAndSurcesCore::UserDataOperator::OPROW;
      // set tangent matrix
      feVolLhs->getOpPtrVector().push_back(
        new OpCalculateScalarFieldValues(string("VALUES")+"_t",headRateAtGaussPts,MBTET)
      );
      feVolLhs->getOpPtrVector().push_back(new OpValuesAtGaussPts(*this,"VALUES"));
      feVolLhs->getOpPtrVector().push_back(new OpFluxDivergenceAtGaussPts(*this,"FLUXES"));
      feVolLhs->getOpPtrVector().push_back(new OpTauDotSigma_HdivHdiv(*this,"FLUXES"));
      feVolLhs->getOpPtrVector().push_back(new OpVU_L2L2(*this,"VALUES"));
      feVolLhs->getOpPtrVector().push_back(new OpVDivSigma_L2Hdiv(*this,"VALUES","FLUXES"));
      feVolLhs->getOpPtrVector().push_back(new OpDivTauU_HdivL2(*this,"FLUXES","VALUES"));

      boost::shared_ptr<FEMethod> null;
      ierr = DMMoFEMTSSetIFunction(dM,"MIX_BCVALUE",feFaceRhs,null,null); CHKERRQ(ierr);
      ierr = DMMoFEMTSSetIFunction(dM,"MIX",feVolRhs,null,null); CHKERRQ(ierr);
      ierr = DMMoFEMTSSetIJacobian(dM,"MIX",feVolLhs,null,null); CHKERRQ(ierr);

      // setting up post-processing
      boost::shared_ptr<PostProcVolumeOnRefinedMesh> post_process(new PostProcVolumeOnRefinedMesh(mField));
      ierr = post_process->generateReferenceElementMesh(); CHKERRQ(ierr);
      ierr = post_process->addFieldValuesPostProc("VALUES"); CHKERRQ(ierr);
      ierr = post_process->addFieldValuesPostProc("VALUES_t"); CHKERRQ(ierr);
      ierr = post_process->addFieldValuesPostProc("FLUXES"); CHKERRQ(ierr);
      int frequency = 1;
      ierr = PetscOptionsBegin(
        PETSC_COMM_WORLD,"Monitor post proc","","none"
      ); CHKERRQ(ierr);
      ierr = PetscOptionsInt(
        "-how_often_output",
        "frequency how often results are dumped on hard disk","",
        frequency,&frequency,NULL
      ); CHKERRQ(ierr);
      ierr = PetscOptionsEnd(); CHKERRQ(ierr);
      tsMonitor = boost::shared_ptr<FEMethod>(new MonitorPostProc(*this,post_process,frequency));
      TsCtx *ts_ctx;
      DMMoFEMGetTsCtx(dM,&ts_ctx);
      ts_ctx->get_postProcess_to_do_Monitor().push_back(tsMonitor);
      PetscFunctionReturn(0);
    }

    Vec D1;

    PetscErrorCode createMatrices() {
      PetscFunctionBegin;
      ierr = DMCreateMatrix(dM,&Aij); CHKERRQ(ierr);
      ierr = DMCreateGlobalVector(dM,&D0); CHKERRQ(ierr);
      ierr = VecDuplicate(D0,&D1); CHKERRQ(ierr);
      ierr = VecDuplicate(D0,&D); CHKERRQ(ierr);
      ierr = VecDuplicate(D0,&F); CHKERRQ(ierr);
      PetscFunctionReturn(0);
    }

    /// \brief destroy matrices
    PetscErrorCode destroyMatrices() {
      PetscFunctionBegin;
      ierr = MatDestroy(&Aij); CHKERRQ(ierr);
      ierr = VecDestroy(&D); CHKERRQ(ierr);
      ierr = VecDestroy(&D0); CHKERRQ(ierr);
      ierr = VecDestroy(&D1); CHKERRQ(ierr);
      ierr = VecDestroy(&F); CHKERRQ(ierr);
      PetscFunctionReturn(0);
    }

    PetscErrorCode calculateEssentialBc() {
      PetscFunctionBegin;
      // clear vectors
      ierr = VecZeroEntries(D0); CHKERRQ(ierr);
      ierr = VecGhostUpdateBegin(D0,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
      ierr = VecGhostUpdateEnd(D0,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
      // clear essential bc indices, it could have dofs from other mesh refinement
      bcIndices.clear();
      // set operator to calculate essential boundary conditions
      ierr = DMoFEMLoopFiniteElements(dM,"MIX_BCFLUX",feFaceBc); CHKERRQ(ierr);
      ierr = VecGhostUpdateBegin(D0,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
      ierr = VecGhostUpdateEnd(D0,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
      ierr = VecAssemblyBegin(D0); CHKERRQ(ierr);
      ierr = VecAssemblyEnd(D0); CHKERRQ(ierr);
      double norm2D0;
      ierr = VecNorm(D0,NORM_2,&norm2D0); CHKERRQ(ierr);
      // ierr = VecView(D0,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);
      PetscPrintf(PETSC_COMM_WORLD,"norm2D0 = %6.4e\n",norm2D0);
      PetscFunctionReturn(0);
    }

    PetscErrorCode calculateInitialPc() {
      PetscFunctionBegin;
      // clear vectors
      ierr = VecZeroEntries(D1); CHKERRQ(ierr);
      ierr = VecGhostUpdateBegin(D1,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
      ierr = VecGhostUpdateEnd(D1,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
      // set operator to calculate essential boundary conditions
      ierr = DMoFEMLoopFiniteElements(dM,"MIX",feVolInitialPc); CHKERRQ(ierr);
      ierr = VecGhostUpdateBegin(D1,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
      ierr = VecGhostUpdateEnd(D1,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
      ierr = VecAssemblyBegin(D1); CHKERRQ(ierr);
      ierr = VecAssemblyEnd(D1); CHKERRQ(ierr);
      double norm2D1;
      ierr = VecNorm(D1,NORM_2,&norm2D1); CHKERRQ(ierr);
      // ierr = VecView(D0,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);
      PetscPrintf(PETSC_COMM_WORLD,"norm2D1 = %6.4e\n",norm2D1);
      PetscFunctionReturn(0);
    }

    /**
     * \brief solve problem
     * @return error code
     */
    PetscErrorCode solveProblem(bool set_initial_pc = true) {
      PetscFunctionBegin;
      if(set_initial_pc) {
        // Set initial head
        ierr = DMoFEMMeshToLocalVector(dM,D1,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
      }

      // Initiate vector from data on the mesh
      ierr = DMoFEMMeshToLocalVector(dM,D,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);

      TS ts;
      ierr = TSCreate(PETSC_COMM_WORLD,&ts); CHKERRQ(ierr);
      ierr = TSSetType(ts,TSBEULER); CHKERRQ(ierr);
      ierr = TSSetIFunction(ts,F,PETSC_NULL,PETSC_NULL); CHKERRQ(ierr);
      ierr = TSSetIJacobian(ts,Aij,Aij,PETSC_NULL,PETSC_NULL); CHKERRQ(ierr);
      double ftime = 1;
      ierr = TSSetDuration(ts,PETSC_DEFAULT,ftime); CHKERRQ(ierr);
      ierr = TSSetFromOptions(ts); CHKERRQ(ierr);
      ierr = TSSetDM(ts,dM); CHKERRQ(ierr);
      #if PETSC_VERSION_GE(3,7,0)
      ierr = TSSetExactFinalTime(ts,TS_EXACTFINALTIME_STEPOVER); CHKERRQ(ierr);
      #endif
      // Set monitor
      TsCtx *ts_ctx;
      DMMoFEMGetTsCtx(dM,&ts_ctx);
      ierr = TSMonitorSet(ts,f_TSMonitorSet,ts_ctx,PETSC_NULL); CHKERRQ(ierr);
      ierr = TSSolve(ts,D); CHKERRQ(ierr);
      ierr = TSGetTime(ts,&ftime); CHKERRQ(ierr);
      PetscInt steps,snesfails,rejects,nonlinits,linits;
      ierr = TSGetTimeStepNumber(ts,&steps); CHKERRQ(ierr);
      ierr = TSGetSNESFailures(ts,&snesfails); CHKERRQ(ierr);
      ierr = TSGetStepRejections(ts,&rejects); CHKERRQ(ierr);
      ierr = TSGetSNESIterations(ts,&nonlinits); CHKERRQ(ierr);
      ierr = TSGetKSPIterations(ts,&linits); CHKERRQ(ierr);
      PetscPrintf(PETSC_COMM_WORLD,
        "steps %D (%D rejected, %D SNES fails), ftime %g, nonlinits %D, linits %D\n",
        steps,rejects,snesfails,ftime,nonlinits,linits
      );
      PetscFunctionReturn(0);
    }

  };

}

#endif //  __UNSATURATD_FLOW_HPP__
