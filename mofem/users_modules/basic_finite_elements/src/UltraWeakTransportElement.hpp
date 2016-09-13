/** \file UltraWeakTransportElement.hpp
 * \brief Ultra weak implementation of transport element
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

#ifndef __ULTRAWEAK_TARNSPORT_ELEMENT_HPP
#define __ULTRAWEAK_TARNSPORT_ELEMENT_HPP

/** \brief Ultra weak transport problem
  \ingroup mofem_ultra_weak_transport_elem

  Note to solve this system you need to use direct solver or proper preconditioner
  for saddle problem.

  It is based on \cite arnold2006differential \cite arnold2012mixed
  <https://www.researchgate.net/profile/Richard_Falk/publication/226454406_Differential_Complexes_and_Stability_of_Finite_Element_Methods_I._The_de_Rham_Complex/links/02e7e5214f0426ff77000000.pdf>

*/
struct UltraWeakTransportElement {

  MoFEM::Interface &mField;

  /// \brief  definition of volume element
  struct MyVolumeFE: public MoFEM::VolumeElementForcesAndSourcesCore {
    MyVolumeFE(MoFEM::Interface &m_field): MoFEM::VolumeElementForcesAndSourcesCore(m_field) {}
    int getRule(int order) { return 2*order+1; };
  };

  MyVolumeFE feVol;

  /** \brief define surface element
    *
    */
  struct MyTriFE: public MoFEM::FaceElementForcesAndSourcesCore {
    MyTriFE(MoFEM::Interface &m_field): MoFEM::FaceElementForcesAndSourcesCore(m_field) {}
    int getRule(int order) { return 2*order+1; };
  };

  MyTriFE feTriFluxValue;

  UltraWeakTransportElement(MoFEM::Interface &m_field):
  mField(m_field),
  feVol(m_field),feTriFluxValue(m_field) {};

  virtual ~UltraWeakTransportElement() {}

  VectorDouble valuesAtGaussPts;
  ublas::vector<VectorDouble > valuesGradientAtGaussPts;
  VectorDouble divergenceAtGaussPts;
  ublas::vector<VectorDouble > fluxesAtGaussPts;

  set<PetscInt> bcIndices;
  PetscErrorCode getDirichletBCIndices(IS *is) {
    PetscFunctionBegin;
    std::vector<PetscInt> ids;
    ids.insert(ids.begin(),bcIndices.begin(),bcIndices.end());
    PetscErrorCode ierr;
    IS is_local;
    ierr = ISCreateGeneral(
      mField.get_comm(),ids.size(),ids.empty()?PETSC_NULL:&ids[0],PETSC_COPY_VALUES,&is_local
    ); CHKERRQ(ierr);
    ierr = ISAllGather(is_local,is); CHKERRQ(ierr);
    ierr = ISDestroy(&is_local); CHKERRQ(ierr);
    PetscFunctionReturn(0);
  }

  virtual PetscErrorCode getFlux(
    const EntityHandle ent,
    const double x,const double y,const double z,
    double &flux) {
    PetscFunctionBegin;
    flux  = 0;
    PetscFunctionReturn(0);
  }

  virtual PetscErrorCode getResistivity(
    const EntityHandle ent,
    const double x,const double y,const double z,
    ublas::matrix<FieldData> &inv_k) {
    PetscFunctionBegin;
    inv_k.resize(3,3);
    bzero(&*inv_k.data().begin(),9*sizeof(FieldData));
    for(int dd = 0;dd<3;dd++) {
      inv_k(dd,dd) = 1;
    }
    PetscFunctionReturn(0);
  }

  virtual PetscErrorCode getBcOnValues(
    const EntityHandle ent,
    const double x,const double y,const double z,
    double &value) {
    PetscFunctionBegin;
    value = 0;
    PetscFunctionReturn(0);
  }

  virtual PetscErrorCode getBcOnFluxes(
    const EntityHandle ent,
    const double x,const double y,const double z,
    double &flux) {
    PetscFunctionBegin;
    flux = 0;
    PetscFunctionReturn(0);
  }

  /** \brief data for evaluation of het conductivity and heat capacity elements
    * \infroup mofem_thermal_elem
    */
  struct BlockData {
    double cOnductivity;
    double cApacity;
    Range tEts; ///< constatins elements in block set
  };
  std::map<int,BlockData> setOfBlocks; ///< maps block set id with appropriate BlockData

  /// \brief add finite elements
  PetscErrorCode addFiniteElements(
    const std::string &fluxes_name,const std::string &values_name,
    const std::string &error_name,const std::string mesh_nodals_positions = "MESH_NODE_POSITIONS") {
    PetscFunctionBegin;

    PetscErrorCode ierr;
    ErrorCode rval;

    ierr = mField.add_finite_element("ULTRAWEAK",MF_ZERO); CHKERRQ(ierr);
    ierr = mField.modify_finite_element_add_field_row("ULTRAWEAK",fluxes_name); CHKERRQ(ierr);
    ierr = mField.modify_finite_element_add_field_col("ULTRAWEAK",fluxes_name); CHKERRQ(ierr);
    ierr = mField.modify_finite_element_add_field_row("ULTRAWEAK",values_name); CHKERRQ(ierr);
    ierr = mField.modify_finite_element_add_field_col("ULTRAWEAK",values_name); CHKERRQ(ierr);
    ierr = mField.modify_finite_element_add_field_data("ULTRAWEAK",fluxes_name); CHKERRQ(ierr);
    ierr = mField.modify_finite_element_add_field_data("ULTRAWEAK",values_name); CHKERRQ(ierr);
    if(mField.check_field(mesh_nodals_positions)) {
      ierr = mField.modify_finite_element_add_field_data("ULTRAWEAK",mesh_nodals_positions); CHKERRQ(ierr);
    }

    ierr = mField.add_finite_element("ULTRAWEAK_ERROR",MF_ZERO); CHKERRQ(ierr);
    ierr = mField.modify_finite_element_add_field_row("ULTRAWEAK_ERROR",error_name); CHKERRQ(ierr);
    ierr = mField.modify_finite_element_add_field_col("ULTRAWEAK_ERROR",error_name); CHKERRQ(ierr);
    ierr = mField.modify_finite_element_add_field_data("ULTRAWEAK_ERROR",fluxes_name); CHKERRQ(ierr);
    ierr = mField.modify_finite_element_add_field_data("ULTRAWEAK_ERROR",values_name); CHKERRQ(ierr);
    ierr = mField.modify_finite_element_add_field_data("ULTRAWEAK_ERROR",error_name); CHKERRQ(ierr);
    if(mField.check_field(mesh_nodals_positions)) {
      ierr = mField.modify_finite_element_add_field_data("ULTRAWEAK_ERROR",mesh_nodals_positions); CHKERRQ(ierr);
    }

    for(_IT_CUBITMESHSETS_BY_BCDATA_TYPE_FOR_LOOP_(mField,BLOCKSET|MAT_THERMALSET,it)) {

      Mat_Thermal temp_data;
      ierr = it->getAttributeDataStructure(temp_data); CHKERRQ(ierr);
      setOfBlocks[it->getMeshsetId()].cOnductivity = temp_data.data.Conductivity;
      setOfBlocks[it->getMeshsetId()].cApacity = temp_data.data.HeatCapacity;
      rval = mField.get_moab().get_entities_by_type(it->meshset,MBTET,setOfBlocks[it->getMeshsetId()].tEts,true); CHKERRQ_MOAB(rval);
      ierr = mField.add_ents_to_finite_element_by_TETs(setOfBlocks[it->getMeshsetId()].tEts,"ULTRAWEAK"); CHKERRQ(ierr);
      ierr = mField.add_ents_to_finite_element_by_TETs(setOfBlocks[it->getMeshsetId()].tEts,"ULTRAWEAK_ERROR"); CHKERRQ(ierr);

    }

    ierr = mField.add_finite_element("ULTRAWEAK_FLUXNEUMANN",MF_ZERO); CHKERRQ(ierr);
    ierr = mField.modify_finite_element_add_field_row("ULTRAWEAK_FLUXNEUMANN",fluxes_name); CHKERRQ(ierr);
    ierr = mField.modify_finite_element_add_field_col("ULTRAWEAK_FLUXNEUMANN",fluxes_name); CHKERRQ(ierr);
    ierr = mField.modify_finite_element_add_field_data("ULTRAWEAK_FLUXNEUMANN",fluxes_name); CHKERRQ(ierr);
    ierr = mField.modify_finite_element_add_field_data("ULTRAWEAK_FLUXNEUMANN",values_name); CHKERRQ(ierr);
    if(mField.check_field(mesh_nodals_positions)) {
      ierr = mField.modify_finite_element_add_field_data("ULTRAWEAK_FLUXNEUMANN",mesh_nodals_positions); CHKERRQ(ierr);
    }

    ierr = mField.add_finite_element("ULTRAWEAK_FLUXDIRICHLET",MF_ZERO); CHKERRQ(ierr);
    ierr = mField.modify_finite_element_add_field_row("ULTRAWEAK_FLUXDIRICHLET",fluxes_name); CHKERRQ(ierr);
    ierr = mField.modify_finite_element_add_field_col("ULTRAWEAK_FLUXDIRICHLET",fluxes_name); CHKERRQ(ierr);
    ierr = mField.modify_finite_element_add_field_data("ULTRAWEAK_FLUXDIRICHLET",fluxes_name); CHKERRQ(ierr);
    ierr = mField.modify_finite_element_add_field_data("ULTRAWEAK_FLUXDIRICHLET",values_name); CHKERRQ(ierr);
    if(mField.check_field(mesh_nodals_positions)) {
      ierr = mField.modify_finite_element_add_field_data("ULTRAWEAK_FLUXDIRICHLET",mesh_nodals_positions); CHKERRQ(ierr);
    }

    PetscFunctionReturn(0);
  }

  /** \brief tau,sigma in Hdiv, calculates Aij = Asemble int sigma_dot_tau dTet
  */
  struct OpTauDotSigma_HdivHdiv: public MoFEM::VolumeElementForcesAndSourcesCore::UserDataOperator {

    UltraWeakTransportElement &cTx;
    Mat Aij;
    Vec F;

    OpTauDotSigma_HdivHdiv(
      UltraWeakTransportElement &ctx,
      const std::string field_name,Mat _Aij,Vec _F
    ):
    MoFEM::VolumeElementForcesAndSourcesCore::UserDataOperator(
      field_name,
      UserDataOperator::OPROW|UserDataOperator::OPROWCOL
    ),
    cTx(ctx),Aij(_Aij),F(_F) {}
    virtual ~OpTauDotSigma_HdivHdiv() {}

    ublas::matrix<FieldData> NN,transNN;
    ublas::matrix<FieldData> invK,invKN;
    VectorDouble Nf;
    VectorDouble invKFlux;

    PetscErrorCode doWork(
      int row_side,int col_side,
      EntityType row_type,EntityType col_type,
      DataForcesAndSurcesCore::EntData &row_data,
      DataForcesAndSurcesCore::EntData &col_data
    ) {
      PetscFunctionBegin;

      PetscErrorCode ierr;

      try {

        if(row_data.getFieldData().size()==0) PetscFunctionReturn(0);
        if(col_data.getFieldData().size()==0) PetscFunctionReturn(0);

        EntityHandle fe_ent = getNumeredEntFiniteElementPtr()->getEnt();

        int nb_row = row_data.getFieldData().size();
        int nb_col = col_data.getFieldData().size();
        NN.resize(nb_row,nb_col);
        bzero(&NN(0,0),NN.data().size()*sizeof(FieldData));

        invKN.resize(3,nb_col);
        int nb_gauss_pts = row_data.getHdivN().size1();
        int gg = 0;
        for(;gg<nb_gauss_pts;gg++) {

          double w = getGaussPts()(3,gg)*getVolume();
          if(getHoGaussPtsDetJac().size()>0) {
            w *= getHoGaussPtsDetJac()(gg);
          }

          const double x = getCoordsAtGaussPts()(gg,0);
          const double y = getCoordsAtGaussPts()(gg,1);
          const double z = getCoordsAtGaussPts()(gg,2);
          ierr = cTx.getResistivity(fe_ent,x,y,z,invK); CHKERRQ(ierr);
          noalias(invKN) = prod(invK,trans(col_data.getHdivN(gg)));
          noalias(NN) += w*prod(row_data.getHdivN(gg),invKN);

        }

        ierr = MatSetValues(
          Aij,
          nb_row,&row_data.getIndices()[0],
          nb_col,&col_data.getIndices()[0],
          &NN(0,0),ADD_VALUES
        ); CHKERRQ(ierr);

        if(row_side != col_side || row_type != col_type) {
          transNN.resize(nb_col,nb_row);
          noalias(transNN) = trans(NN);
          ierr = MatSetValues(
            Aij,
            nb_col,&col_data.getIndices()[0],
            nb_row,&row_data.getIndices()[0],
            &transNN(0,0),ADD_VALUES
          ); CHKERRQ(ierr);
        }


      } catch (const std::exception& ex) {
        std::ostringstream ss;
        ss << "throw in method: " << ex.what() << std::endl;
        SETERRQ(PETSC_COMM_SELF,MOFEM_STD_EXCEPTION_THROW,ss.str().c_str());
      }

      PetscFunctionReturn(0);
    }


    PetscErrorCode doWork(
      int side,EntityType type,DataForcesAndSurcesCore::EntData &data
    ) {
      PetscFunctionBegin;

      PetscErrorCode ierr;

      try {

        if(data.getFieldData().size()==0) PetscFunctionReturn(0);

        EntityHandle fe_ent = getNumeredEntFiniteElementPtr()->getEnt();

        int nb_row = data.getFieldData().size();
        Nf.resize(nb_row);
        bzero(&*Nf.data().begin(),Nf.data().size()*sizeof(FieldData));

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
          ierr = cTx.getResistivity(fe_ent,x,y,z,invK); CHKERRQ(ierr);

          invKFlux.resize(3);
          noalias(invKFlux) = prod(invK,cTx.fluxesAtGaussPts[gg]);
          noalias(Nf) += w*prod(data.getHdivN(gg),invKFlux);

        }

        ierr = VecSetValues(
          F,nb_row,&data.getIndices()[0],
          &Nf[0],ADD_VALUES
        ); CHKERRQ(ierr);

      } catch (const std::exception& ex) {
        std::ostringstream ss;
        ss << "throw in method: " << ex.what() << std::endl;
        SETERRQ(PETSC_COMM_SELF,MOFEM_STD_EXCEPTION_THROW,ss.str().c_str());
      }

      PetscFunctionReturn(0);
    }

  };

  /** \brief u in L2 and tau in Hdiv, calculates Aij = Asemble int u * div(tau) dTet
    */
  struct OpDivTauU_HdivL2: public MoFEM::VolumeElementForcesAndSourcesCore::UserDataOperator {

    UltraWeakTransportElement &cTx;
    Mat Aij;
    Vec F;

    OpDivTauU_HdivL2(
      UltraWeakTransportElement &ctx,
      const std::string field_name_row,string field_name_col,Mat _Aij,Vec _F
    ):
    MoFEM::VolumeElementForcesAndSourcesCore::UserDataOperator(
      field_name_row,field_name_col,
      UserDataOperator::OPROW
    ),
    cTx(ctx),
    Aij(_Aij),
    F(_F) {
      //this operator is not symmetric settig this varible makes element
      //operator to loop over element entities (subsimplicies) without
      //assumption that off-diagonal matrices are symmetric.
      sYmm = false;
    }
    virtual ~OpDivTauU_HdivL2() {}

    VectorDouble div_vec,Nf;

    PetscErrorCode doWork(
      int side,EntityType type,
      DataForcesAndSurcesCore::EntData &data
    ) {
      PetscFunctionBegin;

      PetscErrorCode ierr;

      try {

        if(data.getFieldData().size()==0) PetscFunctionReturn(0);

        int nb_row = data.getIndices().size();
        Nf.resize(nb_row);
        bzero(&*Nf.data().begin(),Nf.data().size()*sizeof(FieldData));

        div_vec.resize(data.getHdivN().size2()/3,0);
        if(div_vec.size()!=data.getIndices().size()) {
          SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"data inconsistency");
        }

        int nb_gauss_pts = data.getN().size1();
        int gg = 0;
        for(;gg<nb_gauss_pts;gg++) {

          double w = getGaussPts()(3,gg)*getVolume();
          if(getHoGaussPtsDetJac().size()>0) {
            w *= getHoGaussPtsDetJac()(gg);
          }

          ierr = getDivergenceOfHDivBaseFunctions(side,type,data,gg,div_vec); CHKERRQ(ierr);

          noalias(Nf) -= w*div_vec*cTx.valuesAtGaussPts[gg];

        }

        ierr = VecSetValues(
          F,nb_row,&data.getIndices()[0],
          &Nf[0],ADD_VALUES
        ); CHKERRQ(ierr);

      } catch (const std::exception& ex) {
        std::ostringstream ss;
        ss << "throw in method: " << ex.what() << std::endl;
        SETERRQ(PETSC_COMM_SELF,MOFEM_STD_EXCEPTION_THROW,ss.str().c_str());
      }

      PetscFunctionReturn(0);
    }

  };

  /** \brief V in L2 and sigma in Hdiv, calculates Aij = Asemble int V * div(sigma) dTet
    */
  struct OpVDotDivSigma_L2Hdiv: public MoFEM::VolumeElementForcesAndSourcesCore::UserDataOperator {

    UltraWeakTransportElement &cTx;
    Mat Aij;
    Vec F;

    OpVDotDivSigma_L2Hdiv(
      UltraWeakTransportElement &ctx,
      const std::string field_name_row,string field_name_col,Mat _Aij,Vec _F
    ):
    MoFEM::VolumeElementForcesAndSourcesCore::UserDataOperator(
      field_name_row,field_name_col,
      UserDataOperator::OPROW|UserDataOperator::OPROWCOL
    ),
    cTx(ctx),Aij(_Aij),F(_F) {

      //this operator is not symmetric setting this variable makes element
      //operator to loop over element entities without
      //assumption that off-diagonal matrices are symmetric.
      sYmm = false;

    }
    virtual ~OpVDotDivSigma_L2Hdiv() {}

    ublas::matrix<FieldData> NN,transNN;
    VectorDouble div_vec,Nf;

    PetscErrorCode doWork(
      int row_side,int col_side,
      EntityType row_type,EntityType col_type,
      DataForcesAndSurcesCore::EntData &row_data,
      DataForcesAndSurcesCore::EntData &col_data
    ) {
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

          ierr = getDivergenceOfHDivBaseFunctions(
            col_side,col_type,col_data,gg,div_vec
          ); CHKERRQ(ierr);

          //FIXME this multiplication should be done in blas or ublas
          for(int rr = 0;rr<nb_row;rr++) {
            for(int cc = 0;cc<nb_col;cc++) {
              NN(rr,cc) -= w*(row_data.getN(gg)[rr]*div_vec[cc]);
            }
          }

        }

        ierr = MatSetValues(
          Aij,
          nb_row,&row_data.getIndices()[0],
          nb_col,&col_data.getIndices()[0],
          &NN(0,0),ADD_VALUES
        ); CHKERRQ(ierr);

        transNN.resize(nb_col,nb_row);
        ublas::noalias(transNN) = trans(NN);
        ierr = MatSetValues(
          Aij,
          nb_col,&col_data.getIndices()[0],
          nb_row,&row_data.getIndices()[0],
          &transNN(0,0),ADD_VALUES
        ); CHKERRQ(ierr);


      } catch (const std::exception& ex) {
        std::ostringstream ss;
        ss << "throw in method: " << ex.what() << std::endl;
        SETERRQ(PETSC_COMM_SELF,MOFEM_STD_EXCEPTION_THROW,ss.str().c_str());
      }

      PetscFunctionReturn(0);
    }

    PetscErrorCode doWork(
      int side,EntityType type,
      DataForcesAndSurcesCore::EntData &data
    ) {
      PetscFunctionBegin;

      PetscErrorCode ierr;

      try {

        if(data.getIndices().size()==0) PetscFunctionReturn(0);
        if(data.getIndices().size()!=data.getN().size2()) {
          SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"data inconsistency");
        }

        int nb_row = data.getIndices().size();
        Nf.resize(nb_row);
        bzero(&*Nf.data().begin(),Nf.data().size()*sizeof(FieldData));

        int nb_gauss_pts = data.getN().size1();
        int gg = 0;
        for(;gg<nb_gauss_pts;gg++) {

          double w = getGaussPts()(3,gg)*getVolume();
          if(getHoGaussPtsDetJac().size()>0) {
            w *= getHoGaussPtsDetJac()(gg);
          }

          noalias(Nf) -= w*data.getN(gg)*cTx.divergenceAtGaussPts[gg];

        }

        ierr = VecSetValues(
          F,nb_row,&data.getIndices()[0],
          &Nf[0],ADD_VALUES
        ); CHKERRQ(ierr);

      } catch (const std::exception& ex) {
        std::ostringstream ss;
        ss << "throw in method: " << ex.what() << std::endl;
        SETERRQ(PETSC_COMM_SELF,1,ss.str().c_str());
      }

      PetscFunctionReturn(0);
    }

  };

  /** \brief calculate source therms
  */
  struct OpL2Source: public MoFEM::VolumeElementForcesAndSourcesCore::UserDataOperator {

    UltraWeakTransportElement &cTx;
    Vec F;

    OpL2Source(
      UltraWeakTransportElement &ctx,
      const std::string field_name,Vec _F
    ):
    MoFEM::VolumeElementForcesAndSourcesCore::UserDataOperator(field_name,UserDataOperator::OPROW),
    cTx(ctx),
    F(_F) {}

    virtual ~OpL2Source() {}

    VectorDouble Nf;
    PetscErrorCode doWork(
      int side,EntityType type,DataForcesAndSurcesCore::EntData &data) {
        PetscFunctionBegin;

        PetscErrorCode ierr;

        try {

          if(data.getFieldData().size()==0) PetscFunctionReturn(0);
          EntityHandle fe_ent = getNumeredEntFiniteElementPtr()->getEnt();

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
            ierr = cTx.getFlux(fe_ent,x,y,z,flux); CHKERRQ(ierr);

            for(int rr = 0;rr<nb_row;rr++) {
              Nf[rr] += w*data.getN(gg)[rr]*flux;
            }

          }

          ierr = VecSetValues(
            F,nb_row,&data.getIndices()[0],
            &Nf[0],ADD_VALUES
          ); CHKERRQ(ierr);

        } catch (const std::exception& ex) {
          std::ostringstream ss;
          ss << "throw in method: " << ex.what() << std::endl;
          SETERRQ(PETSC_COMM_SELF,1,ss.str().c_str());
        }

        PetscFunctionReturn(0);
      }

    };

  /** \brief calualte F = int_\gamma tau*n u_bar d d\Gamma
    */
  struct OpRhsBcOnValues: public MoFEM::FaceElementForcesAndSourcesCore::UserDataOperator {

    UltraWeakTransportElement &cTx;
    Vec F;

    OpRhsBcOnValues(
      UltraWeakTransportElement &ctx,const std::string field_name,Vec _F):
      MoFEM::FaceElementForcesAndSourcesCore::UserDataOperator(field_name,UserDataOperator::OPROW),
      cTx(ctx),F(_F) {}

    VectorDouble Nf;
    PetscErrorCode doWork(
      int side,EntityType type,DataForcesAndSurcesCore::EntData &data
    ) {
      PetscFunctionBegin;

      PetscErrorCode ierr;

      try {

        if(data.getFieldData().size()==0) PetscFunctionReturn(0);
        EntityHandle fe_ent = getNumeredEntFiniteElementPtr()->getEnt();

        Nf.resize(data.getIndices().size());
        bzero(&*Nf.data().begin(),Nf.size()*sizeof(FieldData));

        int nb_gauss_pts = data.getHdivN().size1();
        for(int gg = 0;gg<nb_gauss_pts;gg++) {

          const double x = getCoordsAtGaussPts()(gg,0);
          const double y = getCoordsAtGaussPts()(gg,1);
          const double z = getCoordsAtGaussPts()(gg,2);

          double value;
          ierr = cTx.getBcOnValues(fe_ent,x,y,z,value); CHKERRQ(ierr);

          double w = getGaussPts()(2,gg)*0.5;
          if(getNormalsAtGaussPt().size1() == (unsigned int)nb_gauss_pts) {
            noalias(Nf) += w*prod(data.getHdivN(gg),getNormalsAtGaussPt(gg))*value;
          } else {
            noalias(Nf) += w*prod(data.getHdivN(gg),getNormal())*value;
          }

        }

        //std::cerr << Nf << std::endl << std::endl;

        ierr = VecSetValues(F,data.getIndices().size(),
        &data.getIndices()[0],&Nf[0],ADD_VALUES); CHKERRQ(ierr);

      } catch (const std::exception& ex) {
        std::ostringstream ss;
        ss << "throw in method: " << ex.what() << std::endl;
        SETERRQ(PETSC_COMM_SELF,MOFEM_STD_EXCEPTION_THROW,ss.str().c_str());
      }

      PetscFunctionReturn(0);
    }

  };

  struct OpEvaluateBcOnFluxes: public MoFEM::FaceElementForcesAndSourcesCore::UserDataOperator {

    UltraWeakTransportElement &cTx;
    Vec X;

    OpEvaluateBcOnFluxes(
      UltraWeakTransportElement &ctx,const std::string field_name,Vec _X):
      MoFEM::FaceElementForcesAndSourcesCore::UserDataOperator(field_name,UserDataOperator::OPROW),
      cTx(ctx),X(_X) {}
    virtual ~OpEvaluateBcOnFluxes() {}

    ublas::matrix<FieldData> NN,L;
    VectorDouble Nf,normalN,x;
    PetscErrorCode doWork(int side,EntityType type,DataForcesAndSurcesCore::EntData &data) {
      PetscFunctionBegin;
      PetscErrorCode ierr;

      try {

        if(data.getFieldData().size()==0) PetscFunctionReturn(0);
        EntityHandle fe_ent = getNumeredEntFiniteElementPtr()->getEnt();

        int nb_dofs = data.getFieldData().size();
        NN.resize(nb_dofs,nb_dofs);
        L.resize(nb_dofs,nb_dofs);
        Nf.resize(nb_dofs);
        normalN.resize(nb_dofs);

        bzero(&*NN.data().begin(),NN.data().size()*sizeof(FieldData));
        bzero(&*Nf.data().begin(),Nf.data().size()*sizeof(FieldData));
        L.resize(NN.size1(),NN.size2());

        int nb_gauss_pts = data.getHdivN().size1();
        for(int gg = 0;gg<nb_gauss_pts;gg++) {

          const double x = getCoordsAtGaussPts()(gg,0);
          const double y = getCoordsAtGaussPts()(gg,1);
          const double z = getCoordsAtGaussPts()(gg,2);

          double flux;
          ierr = cTx.getBcOnFluxes(fe_ent,x,y,z,flux); CHKERRQ(ierr);

          //std::cerr << data.getHdivN() << std::endl;

          double area;
          if(getNormalsAtGaussPt().size1() == (unsigned int)nb_gauss_pts) {
            area = 2.*norm_2(getNormalsAtGaussPt(gg));
            noalias(normalN) = prod(data.getHdivN(gg),getNormalsAtGaussPt(gg))/area;
          } else {
            area = 2.*getArea();
            noalias(normalN) = prod(data.getHdivN(gg),getNormal())/area;
          }

          double w = getGaussPts()(2,gg);
          noalias(NN) += w*outer_prod(normalN,normalN);
          noalias(Nf) -= w*normalN*flux;

        }

        cTx.bcIndices.insert(data.getIndices().begin(),data.getIndices().end());

        cholesky_decompose(NN,L);
        cholesky_solve(L,Nf,ublas::lower());

        ierr = VecSetValues(X,data.getIndices().size(),&data.getIndices()[0],&Nf[0],INSERT_VALUES); CHKERRQ(ierr);

      } catch (const std::exception& ex) {
        std::ostringstream ss;
        ss << "throw in method: " << ex.what() << std::endl;
        SETERRQ(PETSC_COMM_SELF,1,ss.str().c_str());
      }

    PetscFunctionReturn(0);
  }

  };

  struct OpValuesAtGaussPts: public MoFEM::VolumeElementForcesAndSourcesCore::UserDataOperator {

    UltraWeakTransportElement &cTx;

    OpValuesAtGaussPts(
      UltraWeakTransportElement &ctx,
      const std::string field_name
    ):
    MoFEM::VolumeElementForcesAndSourcesCore::UserDataOperator(field_name,UserDataOperator::OPROW),
    cTx(ctx) {}

    virtual ~OpValuesAtGaussPts() {}

    PetscErrorCode doWork(int side,EntityType type,DataForcesAndSurcesCore::EntData &data) {
      PetscFunctionBegin;

      try {

        if(data.getFieldData().size() == 0)  PetscFunctionReturn(0);

        int nb_gauss_pts = data.getN().size1();
        cTx.valuesAtGaussPts.resize(nb_gauss_pts);
        for(int gg = 0;gg<nb_gauss_pts;gg++) {
          cTx.valuesAtGaussPts[gg] = inner_prod( trans(data.getN(gg)), data.getFieldData() );
        }

      } catch (const std::exception& ex) {
        std::ostringstream ss;
        ss << "throw in method: " << ex.what() << std::endl;
        SETERRQ(PETSC_COMM_SELF,1,ss.str().c_str());
      }

      PetscFunctionReturn(0);
    }

  };

  struct OpValuesGradientAtGaussPts: public MoFEM::VolumeElementForcesAndSourcesCore::UserDataOperator {

    UltraWeakTransportElement &cTx;

    OpValuesGradientAtGaussPts(
      UltraWeakTransportElement &ctx,
      const std::string field_name
    ):
    MoFEM::VolumeElementForcesAndSourcesCore::UserDataOperator(field_name,UserDataOperator::OPROW),
    cTx(ctx) {}
    virtual ~OpValuesGradientAtGaussPts() {}

    PetscErrorCode doWork(int side,EntityType type,DataForcesAndSurcesCore::EntData &data) {
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
        std::ostringstream ss;
        ss << "throw in method: " << ex.what() << std::endl;
        SETERRQ(PETSC_COMM_SELF,1,ss.str().c_str());
      }

      PetscFunctionReturn(0);
    }

  };

  struct OpFluxDivergenceAtGaussPts: public MoFEM::VolumeElementForcesAndSourcesCore::UserDataOperator {

    UltraWeakTransportElement &cTx;

    OpFluxDivergenceAtGaussPts(
      UltraWeakTransportElement &ctx,
      const std::string field_name
    ):
    MoFEM::VolumeElementForcesAndSourcesCore::UserDataOperator(field_name,UserDataOperator::OPROW),
    cTx(ctx) {}

    virtual ~OpFluxDivergenceAtGaussPts() {}

    VectorDouble div_vec;
    PetscErrorCode doWork(int side,EntityType type,DataForcesAndSurcesCore::EntData &data) {
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
          cTx.fluxesAtGaussPts[gg].clear();
        }
        cTx.divergenceAtGaussPts.clear();
      }

      div_vec.resize(nb_dofs);
      for(int gg = 0;gg<nb_gauss_pts;gg++) {

        ierr = getDivergenceOfHDivBaseFunctions(side,type,data,gg,div_vec); CHKERRQ(ierr);
        cTx.divergenceAtGaussPts[gg] += inner_prod(div_vec,data.getFieldData());
        noalias(cTx.fluxesAtGaussPts[gg]) += prod( trans(data.getHdivN(gg)), data.getFieldData());

      }

    } catch (const std::exception& ex) {
      std::ostringstream ss;
      ss << "throw in method: " << ex.what() << std::endl;
      SETERRQ(PETSC_COMM_SELF,1,ss.str().c_str());
    }

    PetscFunctionReturn(0);
  }

};


  /** \brief calculate error evaluator
    */
  struct OpError_L2Norm: public MoFEM::VolumeElementForcesAndSourcesCore::UserDataOperator {

    UltraWeakTransportElement &cTx;

    OpError_L2Norm(
      UltraWeakTransportElement &ctx,
      const std::string field_name):
      MoFEM::VolumeElementForcesAndSourcesCore::UserDataOperator(field_name,UserDataOperator::OPROW),
      cTx(ctx) {}
    virtual ~OpError_L2Norm() {}

    VectorDouble deltaFlux;
    PetscErrorCode doWork(
      int side,EntityType type,DataForcesAndSurcesCore::EntData &data) {
      PetscFunctionBegin;

      PetscErrorCode ierr;
      ErrorCode rval;

      try {

        if(data.getFieldData().size() == 0)  PetscFunctionReturn(0);
        int nb_gauss_pts = data.getN().size1();
        EntityHandle fe_ent = getNumeredEntFiniteElementPtr()->getEnt();

        double def_VAL = 0;
        Tag th_error_div_sigma;
        rval = cTx.mField.get_moab().tag_get_handle("ERRORL2_DIVSIGMA_F",1,MB_TYPE_DOUBLE,th_error_div_sigma,MB_TAG_CREAT|MB_TAG_SPARSE,&def_VAL); CHKERRQ_MOAB(rval);
        Tag th_error_flux;
        rval = cTx.mField.get_moab().tag_get_handle("ERRORL2_FLUX",1,MB_TYPE_DOUBLE,th_error_flux,MB_TAG_CREAT|MB_TAG_SPARSE,&def_VAL); CHKERRQ_MOAB(rval);
        if(type == MBTRI && side == 0) {
          cTx.mField.get_moab().tag_set_data(th_error_div_sigma,&fe_ent,1,&def_VAL);
        }

        double* error_div_ptr;
        rval = cTx.mField.get_moab().tag_get_by_ptr(th_error_div_sigma,&fe_ent,1,(const void**)&error_div_ptr); CHKERRQ_MOAB(rval);
        double* error_flux_ptr;
        rval = cTx.mField.get_moab().tag_get_by_ptr(th_error_flux,&fe_ent,1,(const void**)&error_flux_ptr); CHKERRQ_MOAB(rval);

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
          ierr = cTx.getFlux(fe_ent,x,y,z,flux); CHKERRQ(ierr);

          if(gg == 0) {
            *error_div_ptr = 0;
            *error_flux_ptr = 0;
          }

          *error_div_ptr += w*( pow(cTx.divergenceAtGaussPts[gg] + flux,2) );
          noalias(deltaFlux) = cTx.fluxesAtGaussPts[gg]+cTx.valuesGradientAtGaussPts[gg];
          *error_flux_ptr += w*( inner_prod(deltaFlux,deltaFlux) );

          //std::cerr << cTx.fluxesAtGaussPts[gg] << " " << cTx.valuesGradientAtGaussPts[gg] << std::endl;

        }

        if(type == MBTET) {
          *error_div_ptr = sqrt(*error_div_ptr);
          *error_flux_ptr = sqrt(*error_flux_ptr);
        }

        const FENumeredDofEntity *dof_ptr;
        ierr = getNumeredEntFiniteElementPtr()->getRowDofsByPetscGlobalDofIdx(
          data.getIndices()[0],&dof_ptr
        ); CHKERRQ(ierr);
        dof_ptr->getFieldData() = *error_flux_ptr;

      } catch (const std::exception& ex) {
        std::ostringstream ss;
        ss << "throw in method: " << ex.what() << std::endl;
        SETERRQ(PETSC_COMM_SELF,1,ss.str().c_str());
      }

      PetscFunctionReturn(0);
    }

  };

};

#endif //__ULTRAWEAK_TARNSPORT_ELEMENT_HPP

/***************************************************************************//**
 * \defgroup mofem_ultra_weak_transport_elem Ultra weak transport element
 * \ingroup user_modules
 ******************************************************************************/
