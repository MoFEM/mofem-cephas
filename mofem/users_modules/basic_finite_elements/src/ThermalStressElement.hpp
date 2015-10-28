/** \file ThermalStressElement.hpp
  \ingroup mofem_thermal_elem
  \brief Implementation of thermal stresses element

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

#ifndef __THERMALSTRESSELEMENT_HPP
#define __THERMALSTRESSELEMENT_HPP

/** \brief Implentation of thermal stress element
\ingroup mofem_thermal_elem
*/
struct ThermalStressElement {

  struct MyVolumeFE: public VolumeElementForcesAndSourcesCore {
    MyVolumeFE(FieldInterface &m_field): VolumeElementForcesAndSourcesCore(m_field) {}
    int getRule(int order) { return order-1; };
  };

  MyVolumeFE feThermalStressRhs;
  MyVolumeFE& getLoopThermalStressRhs() { return feThermalStressRhs; }

  FieldInterface &mField;
  ThermalStressElement( FieldInterface &m_field):
  feThermalStressRhs(m_field),
  mField(m_field) {

  }

  struct BlockData {
    double youngModulus;
    double poissonRatio;
    double thermalExpansion;
    double refTemperature;
    BlockData(): refTemperature(0) {}
    Range tEts;
  };
  map<int,BlockData> setOfBlocks;

  struct CommonData {
    ublas::vector<double> temperatureAtGaussPts;
  };
  CommonData commonData;

  struct OpGetTemperatureAtGaussPts: public VolumeElementForcesAndSourcesCore::UserDataOperator {

    CommonData &commonData;
    int verb;
    OpGetTemperatureAtGaussPts(const string field_name,CommonData &common_data,int _verb = 0):
    VolumeElementForcesAndSourcesCore::UserDataOperator(field_name,UserDataOperator::OPROW),
    commonData(common_data),
    verb(_verb) {}

    PetscErrorCode doWork(int side,EntityType type,DataForcesAndSurcesCore::EntData &data) {
      PetscFunctionBegin;
      try {
        if(data.getFieldData().size()==0) PetscFunctionReturn(0);
        int nb_dofs = data.getFieldData().size();
        int nb_gauss_pts = data.getN().size1();
        //initialize
        commonData.temperatureAtGaussPts.resize(nb_gauss_pts);
        if(type == MBVERTEX) {
          fill(commonData.temperatureAtGaussPts.begin(),commonData.temperatureAtGaussPts.end(),0);
        }
        for(int gg = 0;gg<nb_gauss_pts;gg++) {
          commonData.temperatureAtGaussPts[gg] += inner_prod(data.getN(gg,nb_dofs),data.getFieldData());
        }
      } catch (const std::exception& ex) {
        ostringstream ss;
        ss << "throw in method: " << ex.what() << endl;
        SETERRQ(PETSC_COMM_SELF,1,ss.str().c_str());
      }
      PetscFunctionReturn(0);
    }

  };


  struct OpThermalStressRhs: public VolumeElementForcesAndSourcesCore::UserDataOperator {

    Vec F;
    BlockData &dAta;
    CommonData &commonData;
    int verb;
    OpThermalStressRhs(const string field_name,Vec _F,BlockData &data,CommonData &common_data,int _verb = 0):
    VolumeElementForcesAndSourcesCore::UserDataOperator(field_name,UserDataOperator::OPROW),
    F(_F),
    dAta(data),
    commonData(common_data),
    verb(_verb) { }

    ublas::vector<double> Nf;
    PetscErrorCode doWork(
      int side,EntityType type,DataForcesAndSurcesCore::EntData &data) {
        PetscFunctionBegin;

        try {

          if(data.getIndices().size()==0) PetscFunctionReturn(0);
          if(dAta.tEts.find(getMoFEMFEPtr()->get_ent())==dAta.tEts.end()) PetscFunctionReturn(0);

          PetscErrorCode ierr;

          const FENumeredDofMoFEMEntity *dof_ptr;
          ierr = getMoFEMFEPtr()->get_row_dofs_by_petsc_gloabl_dof_idx(data.getIndices()[0],&dof_ptr); CHKERRQ(ierr);
          int rank = dof_ptr->get_max_rank();
          if(rank != 3) {
            SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"data inconsistency");
          }

          unsigned int nb_dofs = data.getIndices().size();
          if(nb_dofs % rank != 0) {
            SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"data inconsistency");
          }
          if(data.getN().size2() < nb_dofs/rank) {
            SETERRQ3(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,
              "data inconsistency N.size2 %d nb_dof %d rank %d",data.getN().size2(),nb_dofs,rank);
            }

            Nf.resize(nb_dofs);
            bzero(&*Nf.data().begin(),nb_dofs*sizeof(FieldData));

            if(verb > 0) {
              if(type == MBVERTEX) {
                cout << commonData.temperatureAtGaussPts << endl;
                cout << "thermal expansion " << dAta.thermalExpansion << endl;
              }
            }

            for(unsigned int gg = 0;gg<data.getN().size1();gg++) {

              //eps_thermal = (phi-phi_ref)*alpha
              //sig_thernal = - (E/1-2mu) * eps_thermal
              //var_eps = [ diff_N[0], diffN[1], diffN[2] ]

              if(dAta.refTemperature != dAta.refTemperature) {
                SETERRQ(PETSC_COMM_SELF,MOFEM_INVALID_DATA ,"invalid data");
              }

              double phi = (commonData.temperatureAtGaussPts[gg]-dAta.refTemperature);
              double eps_thermal = dAta.thermalExpansion*phi;
              double sig_thermal = -eps_thermal*(dAta.youngModulus/(1.-2*dAta.poissonRatio));

              double val = sig_thermal*getVolume()*getGaussPts()(3,gg);

              double *diff_N;
              diff_N = &data.getDiffN()(gg,0);
              cblas_daxpy(nb_dofs,val,diff_N,1,&Nf[0],1);

            }

            /*for(unsigned int ii = 0;ii<Nf.size();ii++) {
            if(Nf[ii] != Nf[ii]) {
            SETERRQ(PETSC_COMM_SELF,MOFEM_INVALID_DATA ,"invalid data");
          }
        }*/

        ierr = VecSetValues(F,data.getIndices().size(),
        &data.getIndices()[0],&Nf[0],ADD_VALUES); CHKERRQ(ierr);

      } catch (const std::exception& ex) {
        ostringstream ss;
        ss << "throw in method: " << ex.what() << endl;
        SETERRQ(PETSC_COMM_SELF,MOFEM_STD_EXCEPTION_THROW,ss.str().c_str());
      }

      PetscFunctionReturn(0);
    }

  };

  PetscErrorCode addThermalSterssElement(
    const string fe_name,const string field_name,const string thermal_field_name,
    const string mesh_nodals_positions = "MESH_NODE_POSITIONS"
  ) {
    PetscFunctionBegin;
    if(mField.check_field(thermal_field_name)) {
      PetscErrorCode ierr;
      ErrorCode rval;
      ierr = mField.add_finite_element(fe_name,MF_ZERO); CHKERRQ(ierr);
      ierr = mField.modify_finite_element_add_field_row(fe_name,field_name); CHKERRQ(ierr);
      ierr = mField.modify_finite_element_add_field_col(fe_name,field_name); CHKERRQ(ierr);
      ierr = mField.modify_finite_element_add_field_data(fe_name,field_name); CHKERRQ(ierr);
      ierr = mField.modify_finite_element_add_field_data(fe_name,thermal_field_name); CHKERRQ(ierr);
      if(mField.check_field(mesh_nodals_positions)) {
        ierr = mField.modify_finite_element_add_field_data(fe_name,mesh_nodals_positions); CHKERRQ(ierr);
      }
      for(_IT_CUBITMESHSETS_BY_BCDATA_TYPE_FOR_LOOP_(mField,BLOCKSET|MAT_ELASTICSET,it)) {
        Mat_Elastic mydata;
        ierr = it->get_attribute_data_structure(mydata); CHKERRQ(ierr);
        setOfBlocks[it->get_msId()].youngModulus = mydata.data.Young;
        setOfBlocks[it->get_msId()].poissonRatio = mydata.data.Poisson;
        setOfBlocks[it->get_msId()].thermalExpansion = mydata.data.ThermalExpansion;
        rval = mField.get_moab().get_entities_by_type(it->meshset,MBTET,setOfBlocks[it->get_msId()].tEts,true); CHKERR_PETSC(rval);
        ierr = mField.add_ents_to_finite_element_by_TETs(setOfBlocks[it->get_msId()].tEts,fe_name); CHKERRQ(ierr);
        double ref_temp;
        PetscBool flg;
        ierr = PetscOptionsGetReal(PETSC_NULL,"-my_ref_temp",&ref_temp,&flg); CHKERRQ(ierr);
        if(flg == PETSC_TRUE) {
          PetscPrintf(mField.get_comm(),"set refernce temperature %3.2f\n",ref_temp);
          setOfBlocks[it->get_msId()].refTemperature = ref_temp;
        }
      }
    }
    PetscFunctionReturn(0);
  }

  PetscErrorCode setThermalStressRhsOperators(string field_name,string thermal_field_name,Vec &F,int verb = 0) {
    PetscFunctionBegin;
    if(mField.check_field(thermal_field_name)) {
      map<int,BlockData>::iterator sit = setOfBlocks.begin();
      for(;sit!=setOfBlocks.end();sit++) {
        //add finite elemen
        feThermalStressRhs.getOpPtrVector().push_back(new OpGetTemperatureAtGaussPts(thermal_field_name,commonData,verb));
        feThermalStressRhs.getOpPtrVector().push_back(new OpThermalStressRhs(field_name,F,sit->second,commonData,verb));
      }
    }
    PetscFunctionReturn(0);
  }


  };

#endif //__THERMALSTRESSELEMENT_HPP
