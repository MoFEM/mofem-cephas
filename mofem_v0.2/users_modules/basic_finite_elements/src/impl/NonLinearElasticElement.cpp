/**
 * \brief Operators and data structures for thermal analyse
 *
 * Implementation of nonlinear elastic element.
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

#include <MoFEM.hpp>

#include <Projection10NodeCoordsOnField.hpp>

#include <boost/numeric/ublas/vector_proxy.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/vector.hpp>

#include <adolc/adolc.h>
#include <NonLinearElasticElement.hpp>

namespace MoFEM {

NonlinearElasticElement::MyVolumeFE::MyVolumeFE(FieldInterface &m_field):
  VolumeElementForcesAndSourcesCore(m_field),
  A(PETSC_NULL),
  F(PETSC_NULL),
  addToRule(1) {}

int NonlinearElasticElement::MyVolumeFE::getRule(int order) { return (order-1+addToRule); };

PetscErrorCode NonlinearElasticElement::MyVolumeFE::preProcess() {
  PetscFunctionBegin;
  PetscErrorCode ierr;

  ierr = VolumeElementForcesAndSourcesCore::preProcess(); CHKERRQ(ierr);

  if(A != PETSC_NULL) {
    snes_B = A;
  }

  if(F != PETSC_NULL) {
    snes_f = F;
  }

  switch (ts_ctx) {
    case CTX_TSSETIFUNCTION: {
      if(!F) {
        snes_ctx = CTX_SNESSETFUNCTION;
        snes_f = ts_F;
      }
      break;
    }
    case CTX_TSSETIJACOBIAN: {
      if(!A) {
        snes_ctx = CTX_SNESSETJACOBIAN;
        snes_B = ts_B;
      }
      break;
    }
    default:
    break;
  }

  int ghosts[] = { 0 };
  int rank;
  MPI_Comm_rank(mField.get_comm(),&rank);

  switch (snes_ctx) {
    case CTX_SNESNONE:
    if(rank == 0) {
      ierr = VecCreateGhost(mField.get_comm(),1,1,1,ghosts,&V); CHKERRQ(ierr);
    } else {
      ierr = VecCreateGhost(mField.get_comm(),0,1,1,ghosts,&V); CHKERRQ(ierr);
    }
    ierr = VecZeroEntries(V); CHKERRQ(ierr);
    ierr = VecGhostUpdateBegin(V,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
    ierr = VecGhostUpdateEnd(V,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
    break;
    default:
    break;
  }

  PetscFunctionReturn(0);
}


PetscErrorCode NonlinearElasticElement::MyVolumeFE::postProcess() {
  PetscFunctionBegin;

  double *array;

  switch (snes_ctx) {
    case CTX_SNESNONE:
      ierr = VecAssemblyBegin(V); CHKERRQ(ierr);
      ierr = VecAssemblyEnd(V); CHKERRQ(ierr);
      ierr = VecGhostUpdateBegin(V,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
      ierr = VecGhostUpdateEnd(V,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
      ierr = VecGhostUpdateBegin(V,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
      ierr = VecGhostUpdateEnd(V,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
      ierr = VecGetArray(V,&array); CHKERRQ(ierr);
      eNergy = array[0];
      ierr = VecRestoreArray(V,&array); CHKERRQ(ierr);
      ierr = VecDestroy(&V); CHKERRQ(ierr);
      break;
    default:
      break;
  }

  ierr = VolumeElementForcesAndSourcesCore::postProcess(); CHKERRQ(ierr);


  PetscFunctionReturn(0);
}

NonlinearElasticElement::NonlinearElasticElement(
  FieldInterface &m_field,short int tag):
  feRhs(m_field),feLhs(m_field),
  feEnergy(m_field),
  mField(m_field),tAg(tag) {}

NonlinearElasticElement::OpGetDataAtGaussPts::OpGetDataAtGaussPts(const string field_name,
  vector<VectorDouble > &values_at_gauss_pts,
  vector<MatrixDouble > &gardient_at_gauss_pts):
  VolumeElementForcesAndSourcesCore::UserDataOperator(field_name,UserDataOperator::OPROW),
  valuesAtGaussPts(values_at_gauss_pts),gradientAtGaussPts(gardient_at_gauss_pts),
  zeroAtType(MBVERTEX) {}

PetscErrorCode NonlinearElasticElement::OpGetDataAtGaussPts::doWork(
  int side,EntityType type,DataForcesAndSurcesCore::EntData &data) {
  PetscFunctionBegin;
  try {

    int nb_dofs = data.getFieldData().size();
    if(nb_dofs == 0) {
      PetscFunctionReturn(0);
    }
    int nb_gauss_pts = data.getN().size1();

    //initialize
    VectorDouble& values = data.getFieldData();
    valuesAtGaussPts.resize(nb_gauss_pts);
    gradientAtGaussPts.resize(nb_gauss_pts);
    for(int gg = 0;gg<nb_gauss_pts;gg++) {
      valuesAtGaussPts[gg].resize(3,false);
      gradientAtGaussPts[gg].resize(3,3,false);
    }

    if(type == zeroAtType) {
      for(int gg = 0;gg<nb_gauss_pts;gg++) {
        valuesAtGaussPts[gg].clear();
        gradientAtGaussPts[gg].clear();
      }
    }

    //cerr << valuesAtGaussPts[0] << " : ";
    for(int gg = 0;gg<nb_gauss_pts;gg++) {
      VectorDouble N = data.getN(gg,nb_dofs/3);
      MatrixDouble diffN = data.getDiffN(gg,nb_dofs/3);
      for(int dd = 0;dd<nb_dofs/3;dd++) {
        for(int rr1 = 0;rr1<3;rr1++) {
          valuesAtGaussPts[gg][rr1] += N[dd]*values[3*dd+rr1];
          for(int rr2 = 0;rr2<3;rr2++) {
            gradientAtGaussPts[gg](rr1,rr2) += diffN(dd,rr2)*values[3*dd+rr1];
          }
        }
      }
    }

    //cerr << row_field_name << " " << col_field_name << endl;
    //cerr << side << " " << type << endl;
    //cerr << values << endl;
    //cerr << valuesAtGaussPts[0] << endl;

  } catch (const std::exception& ex) {
    ostringstream ss;
    ss << "throw in method: " << ex.what() << endl;
    SETERRQ(PETSC_COMM_SELF,1,ss.str().c_str());
  }

  PetscFunctionReturn(0);
}

NonlinearElasticElement::OpGetCommonDataAtGaussPts::OpGetCommonDataAtGaussPts(const string field_name,CommonData &common_data):
  OpGetDataAtGaussPts(field_name,
  common_data.dataAtGaussPts[field_name],
  common_data.gradAtGaussPts[field_name]) {}

NonlinearElasticElement::OpJacobianPiolaKirchhoffStress::OpJacobianPiolaKirchhoffStress(
  const string field_name,
  BlockData &data,
  CommonData &common_data,
  int tag,
  bool jacobian,
  bool ale,
  bool field_disp):
  VolumeElementForcesAndSourcesCore::UserDataOperator(field_name,UserDataOperator::OPROW),
  dAta(data),
  commonData(common_data),
  tAg(tag),
  adlocReturnValue(0),
  jAcobian(jacobian),
  fUnction(!jacobian),
  aLe(ale),
  fieldDisp(field_disp)
  {}

PetscErrorCode NonlinearElasticElement::OpJacobianPiolaKirchhoffStress::calculateStress() {
  PetscFunctionBegin;

  try {

    PetscErrorCode ierr;
    ierr = dAta.materialAdoublePtr->calculateP_PiolaKirchhoffI(dAta,getMoFEMFEPtr()); CHKERRQ(ierr);
    if(aLe) {
      dAta.materialAdoublePtr->P =
        dAta.materialAdoublePtr->detH*prod(dAta.materialAdoublePtr->P,trans(dAta.materialAdoublePtr->invH));
    }
    commonData.sTress[0].resize(3,3,false);
    for(int dd1 = 0;dd1<3;dd1++) {
      for(int dd2 = 0;dd2<3;dd2++) {
        dAta.materialAdoublePtr->P(dd1,dd2) >>= (commonData.sTress[0])(dd1,dd2);
      }
    }

  } catch (const std::exception& ex) {
    ostringstream ss;
    ss << "throw in method: " << ex.what() << endl;
    SETERRQ(PETSC_COMM_SELF,1,ss.str().c_str());
  }

  PetscFunctionReturn(0);
}

PetscErrorCode NonlinearElasticElement::OpJacobianPiolaKirchhoffStress::doWork(
  int row_side,EntityType row_type,DataForcesAndSurcesCore::EntData &row_data
) {
  PetscFunctionBegin;

  //do it only once, no need to repeat this for edges,faces or tets
  if(row_type != MBVERTEX) PetscFunctionReturn(0);

  PetscErrorCode ierr;
  if(dAta.tEts.find(getMoFEMFEPtr()->get_ent()) == dAta.tEts.end()) {
    PetscFunctionReturn(0);
  }

  int nb_dofs = row_data.getFieldData().size();
  if(nb_dofs==0) PetscFunctionReturn(0);
  dAta.materialAdoublePtr->commonDataPtr = &commonData;
  dAta.materialAdoublePtr->opPtr = this;

  try {

    int nb_gauss_pts = row_data.getN().size1();
    commonData.sTress.resize(nb_gauss_pts);
    commonData.jacStressRowPtr.resize(nb_gauss_pts);
    commonData.jacStress.resize(nb_gauss_pts);

    ptrh = &(commonData.gradAtGaussPts[commonData.spatialPositions]);
    if(aLe) {
      ptrH = &(commonData.gradAtGaussPts[commonData.meshPositions]);
    }

    for(int gg = 0;gg<nb_gauss_pts;gg++) {

      dAta.materialAdoublePtr->gG = gg;

      if(gg == 0) {

        trace_on(tAg);

        dAta.materialAdoublePtr->F.resize(3,3,false);

        if(!aLe) {

          nb_active_variables = 0;
          for(int dd1 = 0;dd1<3;dd1++) {
            for(int dd2 = 0;dd2<3;dd2++) {
              dAta.materialAdoublePtr->F(dd1,dd2) <<= (*ptrh)[gg](dd1,dd2);
              if(fieldDisp) {
                if(dd1 == dd2) {
                  dAta.materialAdoublePtr->F(dd1,dd2) += 1;
                }
              }
              nb_active_variables++;
            }
          }


        } else {

          nb_active_variables = 0;

          dAta.materialAdoublePtr->h.resize(3,3,false);
          for(int dd1 = 0;dd1<3;dd1++) {
            for(int dd2 = 0;dd2<3;dd2++) {
              dAta.materialAdoublePtr->h(dd1,dd2) <<= (*ptrh)[gg](dd1,dd2);
              nb_active_variables++;
            }
          }

          dAta.materialAdoublePtr->H.resize(3,3,false);
          for(int dd1 = 0;dd1<3;dd1++) {
            for(int dd2 = 0;dd2<3;dd2++) {
              dAta.materialAdoublePtr->H(dd1,dd2) <<= (*ptrH)[gg](dd1,dd2);
              nb_active_variables++;
            }
          }

          ierr = dAta.materialAdoublePtr->dEterminatnt(
            dAta.materialAdoublePtr->H,dAta.materialAdoublePtr->detH
          ); CHKERRQ(ierr);
          dAta.materialAdoublePtr->invH.resize(3,3,false);
          ierr = dAta.materialAdoublePtr->iNvert(
            dAta.materialAdoublePtr->detH,dAta.materialAdoublePtr->H,dAta.materialAdoublePtr->invH
          ); CHKERRQ(ierr);
          noalias(dAta.materialAdoublePtr->F) = prod(
            dAta.materialAdoublePtr->h,dAta.materialAdoublePtr->invH
          );

        }

        ierr = dAta.materialAdoublePtr->setUserActiveVariables(nb_active_variables); CHKERRQ(ierr);
        ierr = calculateStress(); CHKERRQ(ierr);

        trace_off();

      }

      activeVariables.resize(nb_active_variables,false);

      if(!aLe) {
        for(int dd1 = 0;dd1<3;dd1++) {
          for(int dd2 = 0;dd2<3;dd2++) {
            activeVariables(dd1*3+dd2) = (*ptrh)[gg](dd1,dd2);
            if(fieldDisp) {
              if(dd1 == dd2) {
                activeVariables(dd1*3+dd2) += 1;
              }
            }
          }
        }
      } else {
        for(int dd1 = 0;dd1<3;dd1++) {
          for(int dd2 = 0;dd2<3;dd2++) {
            activeVariables(dd1*3+dd2) = (*ptrh)[gg](dd1,dd2);
          }
        }
        for(int dd1 = 0;dd1<3;dd1++) {
          for(int dd2 = 0;dd2<3;dd2++) {
            activeVariables(9+dd1*3+dd2) = (*ptrH)[gg](dd1,dd2);
          }
        }
      }
      ierr = dAta.materialAdoublePtr->setUserActiveVariables(activeVariables); CHKERRQ(ierr);

      if(fUnction) {
        commonData.sTress[gg].resize(3,3,false);
        int r;
        //play recorder for values
        r = ::function(tAg,9,nb_active_variables,&activeVariables[0],&commonData.sTress[gg](0,0));
        if(r<adlocReturnValue) { // function is locally analytic
          SETERRQ1(PETSC_COMM_SELF,MOFEM_OPERATION_UNSUCCESSFUL,"ADOL-C function evaluation with error r = %d",r);
        }
      }

      if(jAcobian){
        if(commonData.jacStressRowPtr[gg].size()!=9) {
          commonData.jacStressRowPtr[gg].resize(9);
          commonData.jacStress[gg].resize(9,nb_active_variables,false);
          for(int dd = 0;dd<9;dd++) {
            (commonData.jacStressRowPtr[gg])[dd] = &(commonData.jacStress[gg](dd,0));
          }
        }
        int r;
        //play recorder for jacobians
        r = jacobian(
          tAg,9,nb_active_variables,
          &activeVariables[0],
          &(commonData.jacStressRowPtr[gg])[0]
        );
        if(r<adlocReturnValue) {
          SETERRQ(PETSC_COMM_SELF,MOFEM_OPERATION_UNSUCCESSFUL,"ADOL-C function evaluation with error");
        }
      }

      }

    } catch (const std::exception& ex) {
      ostringstream ss;
      ss << "throw in method: " << ex.what() << endl;
      SETERRQ(PETSC_COMM_SELF,1,ss.str().c_str());
    }

  PetscFunctionReturn(0);
}

NonlinearElasticElement::OpRhsPiolaKirchhoff::OpRhsPiolaKirchhoff(const string field_name,BlockData &data,CommonData &common_data):
  VolumeElementForcesAndSourcesCore::UserDataOperator(field_name,UserDataOperator::OPROW),
  dAta(data),
  commonData(common_data),
  aLe(false) {}

PetscErrorCode NonlinearElasticElement::OpRhsPiolaKirchhoff::aSemble(
  int row_side,EntityType row_type,DataForcesAndSurcesCore::EntData &row_data
) {
  PetscFunctionBegin;
  PetscErrorCode ierr;

  int nb_dofs = row_data.getIndices().size();

  int *indices_ptr = &row_data.getIndices()[0];
  if(!dAta.forcesOnlyOnEntitiesRow.empty()) {
    iNdices.resize(nb_dofs,false);
    noalias(iNdices) = row_data.getIndices();
    indices_ptr = &iNdices[0];
    ublas::vector<const FEDofMoFEMEntity*>& dofs = row_data.getFieldDofs();
    ublas::vector<const FEDofMoFEMEntity*>::iterator dit = dofs.begin();
    for(int ii = 0;dit!=dofs.end();dit++,ii++) {
      if(dAta.forcesOnlyOnEntitiesRow.find((*dit)->get_ent())==dAta.forcesOnlyOnEntitiesRow.end()) {
        iNdices[ii] = -1;
      }
    }
  }

  ierr = VecSetValues(
    getFEMethod()->snes_f,
    nb_dofs,
    indices_ptr,&nf[0],
    ADD_VALUES
  ); CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

PetscErrorCode NonlinearElasticElement::OpRhsPiolaKirchhoff::doWork(
  int row_side,EntityType row_type,DataForcesAndSurcesCore::EntData &row_data
) {
  PetscFunctionBegin;

  PetscErrorCode ierr;
  if(dAta.tEts.find(getMoFEMFEPtr()->get_ent()) == dAta.tEts.end()) {
    PetscFunctionReturn(0);
  }
  if(row_data.getIndices().size()==0) PetscFunctionReturn(0);
  int nb_dofs = row_data.getIndices().size();

  try {

    nf.resize(nb_dofs,false);
    nf.clear();

    for(unsigned int gg = 0;gg<row_data.getN().size1();gg++) {
      //diffN - on rows has degrees of freedom
      //diffN - on columns has rerevatives direvatives of shape functin
      const DataForcesAndSurcesCore::MatrixAdaptor &diffN = row_data.getDiffN(gg,nb_dofs/3);
      const MatrixDouble& stress = commonData.sTress[gg];

      double val = getVolume()*getGaussPts()(3,gg);
      if((!aLe)&&getHoGaussPtsDetJac().size()>0) {
        val *= getHoGaussPtsDetJac()[gg]; ///< higher order geometry
      }
      for(int dd = 0;dd<nb_dofs/3;dd++) {
        for(int rr = 0;rr<3;rr++) {
          for(int nn = 0;nn<3;nn++) {
            nf[3*dd+rr] += val*diffN(dd,nn)*stress(rr,nn);
          }
        }
      }

    }

    if((unsigned int)nb_dofs > 3*row_data.getN().size2()) {
      SETERRQ(PETSC_COMM_SELF,1,"data inconsistency");
    }

    ierr = aSemble(row_side,row_type,row_data); CHKERRQ(ierr);

  } catch (const std::exception& ex) {
    ostringstream ss;
    ss << "throw in method: " << ex.what() << endl;
    SETERRQ(PETSC_COMM_SELF,1,ss.str().c_str());
  }

  PetscFunctionReturn(0);
}

NonlinearElasticElement::OpEnergy::OpEnergy(const string field_name,BlockData &data,CommonData &common_data,Vec *v_ptr,bool field_disp):
  VolumeElementForcesAndSourcesCore::UserDataOperator(field_name,UserDataOperator::OPROW),
  dAta(data),commonData(common_data),
  Vptr(v_ptr),
  fieldDisp(field_disp) { }

PetscErrorCode NonlinearElasticElement::OpEnergy::doWork(
  int row_side,EntityType row_type,DataForcesAndSurcesCore::EntData &row_data) {
  PetscFunctionBegin;

  PetscErrorCode ierr;
  if(row_type != MBVERTEX) PetscFunctionReturn(0);
  if(dAta.tEts.find(getMoFEMFEPtr()->get_ent()) == dAta.tEts.end()) {
    PetscFunctionReturn(0);
  }

  try {

    vector<MatrixDouble > &F = (commonData.gradAtGaussPts[commonData.spatialPositions]);

    for(unsigned int gg = 0;gg<row_data.getN().size1();gg++) {
      double val = getVolume()*getGaussPts()(3,gg);
      if(getHoGaussPtsDetJac().size()>0) {
        val *= getHoGaussPtsDetJac()[gg]; ///< higher order geometry
      }

      dAta.materialDoublePtr->F.resize(3,3,false);
      noalias(dAta.materialDoublePtr->F) = F[gg];
      if(fieldDisp) {
        for(int dd = 0;dd<3;dd++) {
          dAta.materialAdoublePtr->F(dd,dd) += 1;
        }
      }
      ierr = dAta.materialDoublePtr->calculateElasticEnergy(dAta,getMoFEMFEPtr()); CHKERRQ(ierr);
      ierr = VecSetValue(*Vptr,0,val*dAta.materialDoublePtr->eNergy,ADD_VALUES); CHKERRQ(ierr);

    }

  } catch (const std::exception& ex) {
    ostringstream ss;
    ss << "throw in method: " << ex.what() << endl;
    SETERRQ(PETSC_COMM_SELF,1,ss.str().c_str());
  }

  PetscFunctionReturn(0);
}


NonlinearElasticElement::OpLhsPiolaKirchhoff_dx::OpLhsPiolaKirchhoff_dx(
  const string vel_field,const string field_name,BlockData &data,CommonData &common_data):
  VolumeElementForcesAndSourcesCore::UserDataOperator(vel_field,field_name,UserDataOperator::OPROWCOL),
  dAta(data),
  commonData(common_data),
  aLe(false) { }

PetscErrorCode NonlinearElasticElement::OpLhsPiolaKirchhoff_dx::getJac(DataForcesAndSurcesCore::EntData &col_data,int gg) {
  PetscFunctionBegin;
  jac.clear();

  int nb_col = col_data.getFieldData().size();
  const DataForcesAndSurcesCore::MatrixAdaptor diffN = col_data.getDiffN(gg,nb_col/3);
  ublas::matrix<double> &jac_stress = commonData.jacStress[gg];
  // FIXME: this is efficiency bottle neck
  for(int dd = 0;dd<nb_col/3;dd++) {
    for(int rr = 0;rr<3;rr++) {
      for(int ii = 0;ii<9;ii++) {
        for(int jj = 0;jj<3;jj++) {
          jac(ii,3*dd+rr) += jac_stress(ii,3*rr+jj)*diffN(dd,jj);
        }
      }
    }
  }

  PetscFunctionReturn(0);
}


PetscErrorCode NonlinearElasticElement::OpLhsPiolaKirchhoff_dx::aSemble(
  int row_side,int col_side,
  EntityType row_type,EntityType col_type,
  DataForcesAndSurcesCore::EntData &row_data,
  DataForcesAndSurcesCore::EntData &col_data
) {
  PetscFunctionBegin;

  PetscErrorCode ierr;
  int nb_row = row_data.getIndices().size();
  int nb_col = col_data.getIndices().size();

  int *row_indices_ptr = &row_data.getIndices()[0];
  int *col_indices_ptr = &col_data.getIndices()[0];

  /*for(int dd1 = 0;dd1<k.size1();dd1++) {
    for(int dd2 = 0;dd2<k.size2();dd2++) {
      if(k(dd1,dd2)!=k(dd1,dd2)) {
        SETERRQ(PETSC_COMM_SELF,1,"Wrong result");
      }
    }
  }*/

  if(!dAta.forcesOnlyOnEntitiesRow.empty()) {
    rowIndices.resize(nb_row,false);
    noalias(rowIndices) = row_data.getIndices();
    row_indices_ptr = &rowIndices[0];
    ublas::vector<const FEDofMoFEMEntity*>& dofs = row_data.getFieldDofs();
    ublas::vector<const FEDofMoFEMEntity*>::iterator dit = dofs.begin();
    for(int ii = 0;dit!=dofs.end();dit++,ii++) {
      if(dAta.forcesOnlyOnEntitiesRow.find((*dit)->get_ent())==dAta.forcesOnlyOnEntitiesRow.end()) {
        rowIndices[ii] = -1;
      }
    }
  }

  if(!dAta.forcesOnlyOnEntitiesCol.empty()) {
    colIndices.resize(nb_col,false);
    noalias(colIndices) = col_data.getIndices();
    col_indices_ptr = &colIndices[0];
    ublas::vector<const FEDofMoFEMEntity*>& dofs = col_data.getFieldDofs();
    ublas::vector<const FEDofMoFEMEntity*>::iterator dit = dofs.begin();
    for(int ii = 0;dit!=dofs.end();dit++,ii++) {
      if(dAta.forcesOnlyOnEntitiesCol.find((*dit)->get_ent())==dAta.forcesOnlyOnEntitiesCol.end()) {
        colIndices[ii] = -1;
      }
    }
  }

  ierr = MatSetValues(getFEMethod()->snes_B,
    nb_row,row_indices_ptr,
    nb_col,col_indices_ptr,
    &k(0,0),ADD_VALUES
  ); CHKERRQ(ierr);

  //is symmetric
  if(row_side != col_side || row_type != col_type) {

    row_indices_ptr = &row_data.getIndices()[0];
    col_indices_ptr = &col_data.getIndices()[0];

    if(!dAta.forcesOnlyOnEntitiesCol.empty()) {
      rowIndices.resize(nb_row,false);
      noalias(rowIndices) = row_data.getIndices();
      row_indices_ptr = &rowIndices[0];
      ublas::vector<const FEDofMoFEMEntity*>& dofs = row_data.getFieldDofs();
      ublas::vector<const FEDofMoFEMEntity*>::iterator dit = dofs.begin();
      for(int ii = 0;dit!=dofs.end();dit++,ii++) {
        if(dAta.forcesOnlyOnEntitiesCol.find((*dit)->get_ent())==dAta.forcesOnlyOnEntitiesCol.end()) {
          rowIndices[ii] = -1;
        }
      }
    }

    if(!dAta.forcesOnlyOnEntitiesRow.empty()) {
      colIndices.resize(nb_col,false);
      noalias(colIndices) = col_data.getIndices();
      col_indices_ptr = &colIndices[0];
      ublas::vector<const FEDofMoFEMEntity*>& dofs = col_data.getFieldDofs();
      ublas::vector<const FEDofMoFEMEntity*>::iterator dit = dofs.begin();
      for(int ii = 0;dit!=dofs.end();dit++,ii++) {
        if(dAta.forcesOnlyOnEntitiesRow.find((*dit)->get_ent())==dAta.forcesOnlyOnEntitiesRow.end()) {
          colIndices[ii] = -1;
        }
      }
    }

    trans_k.resize(nb_col,nb_row,false);
    noalias(trans_k) = trans(k);
    ierr = MatSetValues(
      getFEMethod()->snes_B,
      nb_col,col_indices_ptr,
      nb_row,row_indices_ptr,
      &trans_k(0,0),ADD_VALUES
    ); CHKERRQ(ierr);
  }

  PetscFunctionReturn(0);
}

PetscErrorCode NonlinearElasticElement::OpLhsPiolaKirchhoff_dx::doWork(
  int row_side,int col_side,
  EntityType row_type,EntityType col_type,
  DataForcesAndSurcesCore::EntData &row_data,
  DataForcesAndSurcesCore::EntData &col_data
) {
  PetscFunctionBegin;

  PetscErrorCode ierr;

  int nb_row = row_data.getIndices().size();
  int nb_col = col_data.getIndices().size();
  if(nb_row == 0) PetscFunctionReturn(0);
  if(nb_col == 0) PetscFunctionReturn(0);

  if(dAta.tEts.find(getMoFEMFEPtr()->get_ent()) == dAta.tEts.end()) {
    PetscFunctionReturn(0);
  }

  try {

    k.resize(nb_row,nb_col,false);
    k.clear();
    jac.resize(9,nb_col,false);

    for(unsigned int gg = 0;gg<row_data.getN().size1();gg++) {

      ierr = getJac(col_data,gg); CHKERRQ(ierr);
      double val = getVolume()*getGaussPts()(3,gg);
      if((!aLe)&&(getHoGaussPtsDetJac().size()>0)) {
        val *= getHoGaussPtsDetJac()[gg]; ///< higher order geometry
      }
      jac *= val;

      const DataForcesAndSurcesCore::MatrixAdaptor &diffN = row_data.getDiffN(gg,nb_row/3);

      { //integrate element stiffness matrix
        for(int dd1 = 0;dd1<nb_row/3;dd1++) {
          for(int rr1 = 0;rr1<3;rr1++) {
            for(int dd2 = 0;dd2<nb_col/3;dd2++) {
              for(int rr2 = 0;rr2<3;rr2++) {
                k(3*dd1+rr1,3*dd2+rr2) +=
                diffN(dd1,0)*jac(3*rr1+0,3*dd2+rr2)+
                diffN(dd1,1)*jac(3*rr1+1,3*dd2+rr2)+
                diffN(dd1,2)*jac(3*rr1+2,3*dd2+rr2);
              }
            }
          }
        }
      }

    }

    ierr = aSemble(row_side,col_side,row_type,col_type,row_data,col_data); CHKERRQ(ierr);

  } catch (const std::exception& ex) {
    ostringstream ss;
    ss << "throw in method: " << ex.what() << endl;
    SETERRQ(PETSC_COMM_SELF,1,ss.str().c_str());
  }

  PetscFunctionReturn(0);
}

NonlinearElasticElement::OpLhsPiolaKirchhoff_dX::OpLhsPiolaKirchhoff_dX(
  const string vel_field,const string field_name,BlockData &data,CommonData &common_data):
  OpLhsPiolaKirchhoff_dx(vel_field,field_name,data,common_data)
  { sYmm = false; }

PetscErrorCode NonlinearElasticElement::OpLhsPiolaKirchhoff_dX::getJac(DataForcesAndSurcesCore::EntData &col_data,int gg) {
  PetscFunctionBegin;
  jac.clear();
  int nb_col = col_data.getFieldData().size();
  const DataForcesAndSurcesCore::MatrixAdaptor diffN = col_data.getDiffN(gg,nb_col/3);
  for(int dd = 0;dd<nb_col/3;dd++) {
    for(int rr = 0;rr<3;rr++) {
      for(int ii = 0;ii<9;ii++) {
        for(int jj = 0;jj<3;jj++) {
          jac(ii,3*dd+rr) += commonData.jacStress[gg](ii,9+3*rr+jj)*diffN(dd,jj);
        }
      }
    }
  }
  PetscFunctionReturn(0);
}

PetscErrorCode NonlinearElasticElement::OpLhsPiolaKirchhoff_dX::aSemble(
  int row_side,int col_side,
  EntityType row_type,EntityType col_type,
  DataForcesAndSurcesCore::EntData &row_data,
  DataForcesAndSurcesCore::EntData &col_data
) {
  PetscFunctionBegin;

  PetscErrorCode ierr;
  int nb_row = row_data.getIndices().size();
  int nb_col = col_data.getIndices().size();

  int *row_indices_ptr = &row_data.getIndices()[0];
  if(!dAta.forcesOnlyOnEntitiesRow.empty()) {
    rowIndices.resize(nb_row,false);
    noalias(rowIndices) = row_data.getIndices();
    row_indices_ptr = &rowIndices[0];
    ublas::vector<const FEDofMoFEMEntity*>& dofs = row_data.getFieldDofs();
    ublas::vector<const FEDofMoFEMEntity*>::iterator dit = dofs.begin();
    for(int ii = 0;dit!=dofs.end();dit++,ii++) {
      if(dAta.forcesOnlyOnEntitiesRow.find((*dit)->get_ent())==dAta.forcesOnlyOnEntitiesRow.end()) {
        rowIndices[ii] = -1;
      }
    }
  }

  int *col_indices_ptr = &col_data.getIndices()[0];
  if(!dAta.forcesOnlyOnEntitiesCol.empty()) {
    colIndices.resize(nb_col,false);
    noalias(colIndices) = col_data.getIndices();
    col_indices_ptr = &colIndices[0];
    ublas::vector<const FEDofMoFEMEntity*>& dofs = col_data.getFieldDofs();
    ublas::vector<const FEDofMoFEMEntity*>::iterator dit = dofs.begin();
    for(int ii = 0;dit!=dofs.end();dit++,ii++) {
      if(dAta.forcesOnlyOnEntitiesCol.find((*dit)->get_ent())==dAta.forcesOnlyOnEntitiesCol.end()) {
        colIndices[ii] = -1;
      }
    }
  }

  /*for(int dd1 = 0;dd1<k.size1();dd1++) {
    for(int dd2 = 0;dd2<k.size2();dd2++) {
      if(k(dd1,dd2)!=k(dd1,dd2)) {
        SETERRQ(PETSC_COMM_SELF,1,"Wrong result");
      }
    }
  }*/

  ierr = MatSetValues(
    getFEMethod()->snes_B,
    nb_row,row_indices_ptr,
    nb_col,col_indices_ptr,
    &k(0,0),ADD_VALUES
  ); CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

NonlinearElasticElement::OpJacobianEshelbyStress::OpJacobianEshelbyStress(
  const string field_name,
  BlockData &data,
  CommonData &common_data,
  int tag,
  bool jacobian,
  bool ale):
  OpJacobianPiolaKirchhoffStress(field_name,data,common_data,tag,jacobian,ale,false)
  {}

PetscErrorCode NonlinearElasticElement::OpJacobianEshelbyStress::calculateStress() {
  PetscFunctionBegin;

  try {

    PetscErrorCode ierr;
    ierr = dAta.materialAdoublePtr->calculateSiGma_EshelbyStress(dAta,getMoFEMFEPtr()); CHKERRQ(ierr);
    if(aLe) {
      dAta.materialAdoublePtr->SiGma =
      dAta.materialAdoublePtr->detH*prod(dAta.materialAdoublePtr->SiGma,trans(dAta.materialAdoublePtr->invH));
    }
    commonData.sTress[0].resize(3,3,false);
    for(int dd1 = 0;dd1<3;dd1++) {
      for(int dd2 = 0;dd2<3;dd2++) {
        dAta.materialAdoublePtr->SiGma(dd1,dd2) >>= (commonData.sTress[0])(dd1,dd2);
      }
    }

  } catch (const std::exception& ex) {
    ostringstream ss;
    ss << "throw in method: " << ex.what() << endl;
    SETERRQ(PETSC_COMM_SELF,1,ss.str().c_str());
  }

  PetscFunctionReturn(0);
}

NonlinearElasticElement::OpRhsEshelbyStrees::OpRhsEshelbyStrees(
  const string field_name,BlockData &data,CommonData &common_data
):
OpRhsPiolaKirchhoff(field_name,data,common_data)
{}

NonlinearElasticElement::OpLhsEshelby_dx::OpLhsEshelby_dx(
  const string vel_field,const string field_name,BlockData &data,CommonData &common_data
):
OpLhsPiolaKirchhoff_dX(vel_field,field_name,data,common_data) {}

PetscErrorCode NonlinearElasticElement::OpLhsEshelby_dx::getJac(DataForcesAndSurcesCore::EntData &col_data,int gg) {
  PetscFunctionBegin;
  jac.clear();
  int nb_col = col_data.getFieldData().size();
  const DataForcesAndSurcesCore::MatrixAdaptor diffN = col_data.getDiffN(gg,nb_col/3);
  for(int dd = 0;dd<nb_col/3;dd++) {
    for(int rr = 0;rr<3;rr++) {
      for(int ii = 0;ii<9;ii++) {
        for(int jj = 0;jj<3;jj++) {
          jac(ii,3*dd+rr) += commonData.jacStress[gg](ii,3*rr+jj)*diffN(dd,jj);
        }
      }
    }
  }
  PetscFunctionReturn(0);
}

NonlinearElasticElement::OpLhsEshelby_dX::OpLhsEshelby_dX(
  const string vel_field,const string field_name,BlockData &data,CommonData &common_data
):
OpLhsPiolaKirchhoff_dx(vel_field,field_name,data,common_data)
{}

PetscErrorCode NonlinearElasticElement::OpLhsEshelby_dX::getJac(DataForcesAndSurcesCore::EntData &col_data,int gg) {
  PetscFunctionBegin;
  jac.clear();
  int nb_col = col_data.getFieldData().size();
  const DataForcesAndSurcesCore::MatrixAdaptor diffN = col_data.getDiffN(gg,nb_col/3);
  for(int dd = 0;dd<nb_col/3;dd++) {
    for(int rr = 0;rr<3;rr++) {
      for(int ii = 0;ii<9;ii++) {
        for(int jj = 0;jj<3;jj++) {
          jac(ii,3*dd+rr) += commonData.jacStress[gg](ii,9+3*rr+jj)*diffN(dd,jj);
        }
      }
    }
  }
  PetscFunctionReturn(0);
}

PetscErrorCode NonlinearElasticElement::setBlocks(
  FunctionsToCalulatePiolaKirchhoffI<double> *materialDoublePtr,
  FunctionsToCalulatePiolaKirchhoffI<adouble> *materialAdoublePtr) {
  PetscFunctionBegin;
  ErrorCode rval;
  PetscErrorCode ierr;

  for(_IT_CUBITMESHSETS_BY_BCDATA_TYPE_FOR_LOOP_(mField,BLOCKSET|MAT_ELASTICSET,it)) {
    Mat_Elastic mydata;
    ierr = it->get_attribute_data_structure(mydata); CHKERRQ(ierr);
    int id = it->get_msId();
    EntityHandle meshset = it->get_meshset();
    rval = mField.get_moab().get_entities_by_type(meshset,MBTET,setOfBlocks[id].tEts,true); CHKERR_PETSC(rval);
    setOfBlocks[id].iD = id;
    setOfBlocks[id].E = mydata.data.Young;
    setOfBlocks[id].PoissonRatio = mydata.data.Poisson;
    setOfBlocks[id].materialDoublePtr = materialDoublePtr;
    setOfBlocks[id].materialAdoublePtr = materialAdoublePtr;
    //cerr << setOfBlocks[id].tEts << endl;
  }

  PetscFunctionReturn(0);
}

PetscErrorCode NonlinearElasticElement::addElement(string element_name,
  string spatial_position_field_name,
  string material_position_field_name,bool ale) {
  PetscFunctionBegin;

  PetscErrorCode ierr;
  //ErrorCode rval;

  ierr = mField.add_finite_element(element_name,MF_ZERO); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_row(element_name,spatial_position_field_name); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_col(element_name,spatial_position_field_name); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_data(element_name,spatial_position_field_name); CHKERRQ(ierr);
  if(mField.check_field(material_position_field_name)) {
    if(ale) {
      ierr = mField.modify_finite_element_add_field_row(element_name,material_position_field_name); CHKERRQ(ierr);
      ierr = mField.modify_finite_element_add_field_col(element_name,material_position_field_name); CHKERRQ(ierr);
    }
    ierr = mField.modify_finite_element_add_field_data(element_name,material_position_field_name); CHKERRQ(ierr);
  }

  map<int,BlockData>::iterator sit = setOfBlocks.begin();
  for(;sit!=setOfBlocks.end();sit++) {
    ierr = mField.add_ents_to_finite_element_by_TETs(sit->second.tEts,element_name); CHKERRQ(ierr);
  }

  PetscFunctionReturn(0);
}

PetscErrorCode NonlinearElasticElement::setOperators(
  string spatial_position_field_name,
  string material_position_field_name,
  bool ale,bool field_disp) {
  PetscFunctionBegin;

  commonData.spatialPositions = spatial_position_field_name;
  commonData.meshPositions = material_position_field_name;

  //Rhs
  feRhs.getOpPtrVector().push_back(new OpGetCommonDataAtGaussPts(spatial_position_field_name,commonData));
  if(mField.check_field(material_position_field_name)) {
    feRhs.getOpPtrVector().push_back(new OpGetCommonDataAtGaussPts(material_position_field_name,commonData));
  }
  map<int,BlockData>::iterator sit = setOfBlocks.begin();
  for(;sit!=setOfBlocks.end();sit++) {
    feRhs.getOpPtrVector().push_back(
      new OpJacobianPiolaKirchhoffStress(spatial_position_field_name,sit->second,commonData,tAg,false,ale,field_disp)
    );
    feRhs.getOpPtrVector().push_back(new OpRhsPiolaKirchhoff(spatial_position_field_name,sit->second,commonData));
  }

  //Energy
  feEnergy.getOpPtrVector().push_back(new OpGetCommonDataAtGaussPts(spatial_position_field_name,commonData));
  if(mField.check_field(material_position_field_name)) {
    feEnergy.getOpPtrVector().push_back(new OpGetCommonDataAtGaussPts(material_position_field_name,commonData));
  }
  sit = setOfBlocks.begin();
  for(;sit!=setOfBlocks.end();sit++) {
    feEnergy.getOpPtrVector().push_back(new OpEnergy(spatial_position_field_name,sit->second,commonData,&feEnergy.V,field_disp));
  }

  //Lhs
  feLhs.getOpPtrVector().push_back(new OpGetCommonDataAtGaussPts(spatial_position_field_name,commonData));
  if(mField.check_field(material_position_field_name)) {
    feLhs.getOpPtrVector().push_back(new OpGetCommonDataAtGaussPts(material_position_field_name,commonData));
  }
  sit = setOfBlocks.begin();
  for(;sit!=setOfBlocks.end();sit++) {
    feLhs.getOpPtrVector().push_back(
      new OpJacobianPiolaKirchhoffStress(spatial_position_field_name,sit->second,commonData,tAg,true,ale,field_disp)
    );
    feLhs.getOpPtrVector().push_back(new OpLhsPiolaKirchhoff_dx(spatial_position_field_name,spatial_position_field_name,sit->second,commonData));
  }

  PetscFunctionReturn(0);
}

}
