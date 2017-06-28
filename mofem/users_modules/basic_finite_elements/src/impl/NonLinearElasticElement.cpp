/**
 * \brief Operators and data structures for nonlinear elastic material
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
using namespace MoFEM;
#include <Projection10NodeCoordsOnField.hpp>

#include <boost/numeric/ublas/vector_proxy.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/vector.hpp>

#include <adolc/adolc.h>
#include <NonLinearElasticElement.hpp>

NonlinearElasticElement::MyVolumeFE::MyVolumeFE(MoFEM::Interface &m_field):
  VolumeElementForcesAndSourcesCore(m_field),
  A(PETSC_NULL),
  F(PETSC_NULL),
  addToRule(1) {}

int NonlinearElasticElement::MyVolumeFE::getRule(int order) { return 2*(order-1)+addToRule; };

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
  MoFEM::Interface &m_field,short int tag):
  feRhs(m_field),feLhs(m_field),
  feEnergy(m_field),
  mField(m_field),tAg(tag) {}

NonlinearElasticElement::OpGetDataAtGaussPts::OpGetDataAtGaussPts(const std::string field_name,
  std::vector<VectorDouble > &values_at_gauss_pts,
  std::vector<MatrixDouble > &gardient_at_gauss_pts):
  VolumeElementForcesAndSourcesCore::UserDataOperator(field_name,UserDataOperator::OPROW),
  valuesAtGaussPts(values_at_gauss_pts),
  gradientAtGaussPts(gardient_at_gauss_pts),
  zeroAtType(MBVERTEX) {}

PetscErrorCode NonlinearElasticElement::OpGetDataAtGaussPts::doWork(
  int side,EntityType type,DataForcesAndSurcesCore::EntData &data
) {
  PetscFunctionBegin;
  try {

    const int nb_dofs = data.getFieldData().size();
    const int nb_base_functions = data.getN().size2();
    if(nb_dofs == 0) {
      PetscFunctionReturn(0);
    }
    const int nb_gauss_pts = data.getN().size1();
    const int rank = data.getFieldDofs()[0]->getNbOfCoeffs();

    //initialize
    if(type == zeroAtType) {
      valuesAtGaussPts.resize(nb_gauss_pts);
      gradientAtGaussPts.resize(nb_gauss_pts);
      for(int gg = 0;gg!=nb_gauss_pts;gg++) {
        valuesAtGaussPts[gg].resize(rank,false);
        gradientAtGaussPts[gg].resize(rank,3,false);
      }
      for(int gg = 0;gg!=nb_gauss_pts;gg++) {
        valuesAtGaussPts[gg].clear();
        gradientAtGaussPts[gg].clear();
      }
    }

    FTensor::Tensor0<double*> base_function = data.getFTensor0N();
    FTensor::Tensor1<double*,3> diff_base_functions = data.getFTensor1DiffN<3>();
    FTensor::Index<'i',3> i;
    FTensor::Index<'j',3> j;

    if(rank==1) {

      for(int gg = 0;gg!=nb_gauss_pts;gg++) {
        FTensor::Tensor0<double*> field_data = data.getFTensor0FieldData();
        double &val = valuesAtGaussPts[gg][0];
        FTensor::Tensor1<double*,3> grad(
          &gradientAtGaussPts[gg](0,0),
          &gradientAtGaussPts[gg](0,1),
          &gradientAtGaussPts[gg](0,2)
        );
        int bb = 0;
        for(;bb!=nb_dofs;bb++) {
          val += base_function*field_data;
          grad(i) += diff_base_functions(i)*field_data;
          ++diff_base_functions;
          ++base_function;
          ++field_data;
        }
        for(;bb!=nb_base_functions;bb++) {
          ++diff_base_functions;
          ++base_function;
        }
      }

    } else if(rank==3) {

      for(int gg = 0;gg!=nb_gauss_pts;gg++) {
        FTensor::Tensor1<double*,3> field_data = data.getFTensor1FieldData<3>();
        FTensor::Tensor1<double*,3> values(
          &valuesAtGaussPts[gg][0],
          &valuesAtGaussPts[gg][1],
          &valuesAtGaussPts[gg][2]
        );
        FTensor::Tensor2<double*,3,3> gradient(
          &gradientAtGaussPts[gg](0,0),&gradientAtGaussPts[gg](0,1),&gradientAtGaussPts[gg](0,2),
          &gradientAtGaussPts[gg](1,0),&gradientAtGaussPts[gg](1,1),&gradientAtGaussPts[gg](1,2),
          &gradientAtGaussPts[gg](2,0),&gradientAtGaussPts[gg](2,1),&gradientAtGaussPts[gg](2,2)
        );
        int bb = 0;
        for(;bb!=nb_dofs/3;bb++) {
          values(i) += base_function*field_data(i);
          gradient(i,j) += field_data(i)*diff_base_functions(j);
          ++diff_base_functions;
          ++base_function;
          ++field_data;
        }
        for(;bb!=nb_base_functions;bb++) {
          ++diff_base_functions;
          ++base_function;
        }
      }

    } else {
      // FIXME: THat part is inefficient
      VectorDouble& values = data.getFieldData();
      //std::cerr << valuesAtGaussPts[0] << " : ";
      for(int gg = 0;gg<nb_gauss_pts;gg++) {
        VectorAdaptor N = data.getN(gg,nb_dofs/rank);
        MatrixAdaptor diffN = data.getDiffN(gg,nb_dofs/rank);
        for(int dd = 0;dd<nb_dofs/rank;dd++) {
          for(int rr1 = 0;rr1<rank;rr1++) {
            valuesAtGaussPts[gg][rr1] += N[dd]*values[rank*dd+rr1];
            for(int rr2 = 0;rr2<3;rr2++) {
              gradientAtGaussPts[gg](rr1,rr2) += diffN(dd,rr2)*values[rank*dd+rr1];
            }
          }
        }
      }
    }

    //std::cerr << row_field_name << " " << col_field_name << std::endl;
    //std::cerr << side << " " << type << std::endl;
    //std::cerr << values << std::endl;
    //std::cerr << valuesAtGaussPts[0] << std::endl;

  } catch (MoFEMException const &e) {
    SETERRQ(PETSC_COMM_SELF,e.errorCode,e.errorMessage);
  }

  PetscFunctionReturn(0);
}

NonlinearElasticElement::OpGetCommonDataAtGaussPts::OpGetCommonDataAtGaussPts(const std::string field_name,CommonData &common_data):
  OpGetDataAtGaussPts(field_name,
  common_data.dataAtGaussPts[field_name],
  common_data.gradAtGaussPts[field_name]) {}

NonlinearElasticElement::OpJacobianPiolaKirchhoffStress::OpJacobianPiolaKirchhoffStress(
  const std::string field_name,
  BlockData &data,
  CommonData &common_data,
  int tag,
  bool jacobian,
  bool ale,
  bool field_disp
):
VolumeElementForcesAndSourcesCore::UserDataOperator(field_name,UserDataOperator::OPROW),
dAta(data),
commonData(common_data),
tAg(tag),
adlocReturnValue(0),
jAcobian(jacobian),
fUnction(!jacobian),
aLe(ale),
fieldDisp(field_disp) {

}

PetscErrorCode NonlinearElasticElement::OpJacobianPiolaKirchhoffStress::calculateStress(const int gg) {
  PetscFunctionBegin;
  try {
    PetscErrorCode ierr;
    ierr = dAta.materialAdoublePtr->calculateP_PiolaKirchhoffI(
      dAta,getNumeredEntFiniteElementPtr()
    ); CHKERRQ(ierr);
    if(aLe) {
      dAta.materialAdoublePtr->P =
        dAta.materialAdoublePtr->detH*prod(dAta.materialAdoublePtr->P,trans(dAta.materialAdoublePtr->invH));
    }
    commonData.sTress[gg].resize(3,3,false);
    for(int dd1 = 0;dd1<3;dd1++) {
      for(int dd2 = 0;dd2<3;dd2++) {
        dAta.materialAdoublePtr->P(dd1,dd2) >>= (commonData.sTress[gg])(dd1,dd2);
      }
    }
    //std::cerr << "P " << dAta.materialAdoublePtr->P << std::endl;
  } catch (MoFEMException const &e) {
    SETERRQ(PETSC_COMM_SELF,e.errorCode,e.errorMessage);
  }

  PetscFunctionReturn(0);
}

PetscErrorCode NonlinearElasticElement::OpJacobianPiolaKirchhoffStress::recordTag(const int gg) {
  PetscFunctionBegin;

  trace_on(tAg);

  try {

    dAta.materialAdoublePtr->F.resize(3,3,false);

    if(!aLe) {

      nbActiveVariables = 0;
      for(int dd1 = 0;dd1<3;dd1++) {
        for(int dd2 = 0;dd2<3;dd2++) {
          dAta.materialAdoublePtr->F(dd1,dd2) <<= (*ptrh)[gg](dd1,dd2);
          if(fieldDisp) {
            if(dd1 == dd2) {
              dAta.materialAdoublePtr->F(dd1,dd2) += 1;
            }
          }
          nbActiveVariables++;
        }
      }

    } else {

      nbActiveVariables = 0;

      dAta.materialAdoublePtr->h.resize(3,3,false);
      for(int dd1 = 0;dd1<3;dd1++) {
        for(int dd2 = 0;dd2<3;dd2++) {
          dAta.materialAdoublePtr->h(dd1,dd2) <<= (*ptrh)[gg](dd1,dd2);
          nbActiveVariables++;
        }
      }

      dAta.materialAdoublePtr->H.resize(3,3,false);
      for(int dd1 = 0;dd1<3;dd1++) {
        for(int dd2 = 0;dd2<3;dd2++) {
          dAta.materialAdoublePtr->H(dd1,dd2) <<= (*ptrH)[gg](dd1,dd2);
          nbActiveVariables++;
        }
      }

      ierr = dAta.materialAdoublePtr->dEterminatnt(
        dAta.materialAdoublePtr->H,dAta.materialAdoublePtr->detH
      ); CHKERRQ(ierr);
      dAta.materialAdoublePtr->invH.resize(3,3,false);
      ierr = dAta.materialAdoublePtr->iNvert(
        dAta.materialAdoublePtr->detH,dAta.materialAdoublePtr->H,
        dAta.materialAdoublePtr->invH
      ); CHKERRQ(ierr);
      noalias(dAta.materialAdoublePtr->F) = prod(
        dAta.materialAdoublePtr->h,dAta.materialAdoublePtr->invH
      );

    }

    ierr = dAta.materialAdoublePtr->setUserActiveVariables(nbActiveVariables); CHKERRQ(ierr);
    ierr = calculateStress(gg); CHKERRQ(ierr);

    trace_off();

  } catch (MoFEMException const &e) {
    SETERRQ(PETSC_COMM_SELF,e.errorCode,e.errorMessage);
  }

  PetscFunctionReturn(0);
}

PetscErrorCode NonlinearElasticElement::OpJacobianPiolaKirchhoffStress::playTag(const int gg) {
  PetscFunctionBegin;

  int r;

  if(fUnction) {
    commonData.sTress[gg].resize(3,3,false);
    //play recorder for values
    r = ::function(tAg,9,nbActiveVariables,&activeVariables[0],&commonData.sTress[gg](0,0));
    if(r<adlocReturnValue) { // function is locally analytic
      SETERRQ1(PETSC_COMM_SELF,MOFEM_OPERATION_UNSUCCESSFUL,"ADOL-C function evaluation with error r = %d",r);
    }
  }

  if(jAcobian) {
    commonData.jacStress[gg].resize(9,nbActiveVariables,false);
    double *jac_ptr[] = {
      &(commonData.jacStress[gg](0,0)),&(commonData.jacStress[gg](1,0)),&(commonData.jacStress[gg](2,0)),
      &(commonData.jacStress[gg](3,0)),&(commonData.jacStress[gg](4,0)),&(commonData.jacStress[gg](5,0)),
      &(commonData.jacStress[gg](6,0)),&(commonData.jacStress[gg](7,0)),&(commonData.jacStress[gg](8,0))
    };
    //play recorder for jacobians
    r = jacobian(
      tAg,9,nbActiveVariables,&activeVariables[0],jac_ptr
    );
    if(r<adlocReturnValue) {
      SETERRQ(PETSC_COMM_SELF,MOFEM_OPERATION_UNSUCCESSFUL,"ADOL-C function evaluation with error");
    }
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
  if(dAta.tEts.find(getNumeredEntFiniteElementPtr()->getEnt()) == dAta.tEts.end()) {
    PetscFunctionReturn(0);
  }

  int nb_dofs = row_data.getFieldData().size();
  if(nb_dofs==0) PetscFunctionReturn(0);
  dAta.materialAdoublePtr->commonDataPtr = &commonData;
  dAta.materialAdoublePtr->opPtr = this;

  try {

    int nb_gauss_pts = row_data.getN().size1();
    commonData.sTress.resize(nb_gauss_pts);
    commonData.jacStress.resize(nb_gauss_pts);

    ptrh = &(commonData.gradAtGaussPts[commonData.spatialPositions]);
    if(aLe) {
      ptrH = &(commonData.gradAtGaussPts[commonData.meshPositions]);
    }

    for(int gg = 0;gg!=nb_gauss_pts;gg++) {

      dAta.materialAdoublePtr->gG = gg;

      // Record tag and calualte stress
      if(recordTagForIntegrationPoint(gg)) {
        ierr = recordTag(gg); CHKERRQ(ierr);
      }

      // Set active variables vector
      if(jAcobian||(!recordTagForIntegrationPoint(gg))) {
        activeVariables.resize(nbActiveVariables,false);
        if(!aLe) {
          for(int dd1 = 0;dd1<3;dd1++) {
            for(int dd2 = 0;dd2<3;dd2++) {
              activeVariables(dd1*3+dd2) = (*ptrh)[gg](dd1,dd2);
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

        // Play tag and calculate stress or tangent
        if(jAcobian||(!recordTagForIntegrationPoint(gg))) {
          ierr = playTag(gg); CHKERRQ(ierr);
        }

      }

    }

  } catch (MoFEMException const &e) {
    SETERRQ(PETSC_COMM_SELF,e.errorCode,e.errorMessage);
  }

  PetscFunctionReturn(0);
}

NonlinearElasticElement::OpJacobianEnergy::OpJacobianEnergy(
  const std::string field_name, ///< field name for spatial positions or displacements
  BlockData &data,
  CommonData &common_data,
  int tag,
  bool gradient,
  bool hessian,
  bool ale,
  bool field_disp
):
VolumeElementForcesAndSourcesCore::UserDataOperator(field_name,UserDataOperator::OPROW),
dAta(data),
commonData(common_data),
tAg(tag),
gRadient(gradient),
hEssian(hessian),
aLe(ale),
fieldDisp(field_disp) {
}

PetscErrorCode NonlinearElasticElement::OpJacobianEnergy::calculateEnergy(const int gg) {
  PetscErrorCode ierr;
  PetscFunctionBegin;
  try {
    ierr = dAta.materialAdoublePtr->calculateElasticEnergy(dAta,getNumeredEntFiniteElementPtr()); CHKERRQ(ierr);
    dAta.materialAdoublePtr->eNergy >>= commonData.eNergy[gg];
  } catch (MoFEMException const &e) {
    SETERRQ(PETSC_COMM_SELF,e.errorCode,e.errorMessage);
  }
  PetscFunctionReturn(0);
}

PetscErrorCode NonlinearElasticElement::OpJacobianEnergy::recordTag(const int gg) {
  PetscErrorCode ierr;
  PetscFunctionBegin;

  trace_on(tAg);

  try {

    dAta.materialAdoublePtr->F.resize(3,3,false);

    if(!aLe) {

      nbActiveVariables = 0;
      for(int dd1 = 0;dd1<3;dd1++) {
        for(int dd2 = 0;dd2<3;dd2++) {
          dAta.materialAdoublePtr->F(dd1,dd2) <<= (*ptrh)[gg](dd1,dd2);
          if(fieldDisp) {
            if(dd1 == dd2) {
              dAta.materialAdoublePtr->F(dd1,dd2) += 1;
            }
          }
          nbActiveVariables++;
        }
      }

    } else {

      nbActiveVariables = 0;

      dAta.materialAdoublePtr->h.resize(3,3,false);
      for(int dd1 = 0;dd1<3;dd1++) {
        for(int dd2 = 0;dd2<3;dd2++) {
          dAta.materialAdoublePtr->h(dd1,dd2) <<= (*ptrh)[gg](dd1,dd2);
          nbActiveVariables++;
        }
      }

      dAta.materialAdoublePtr->H.resize(3,3,false);
      for(int dd1 = 0;dd1<3;dd1++) {
        for(int dd2 = 0;dd2<3;dd2++) {
          dAta.materialAdoublePtr->H(dd1,dd2) <<= (*ptrH)[gg](dd1,dd2);
          nbActiveVariables++;
        }
      }

      ierr = dAta.materialAdoublePtr->dEterminatnt(
        dAta.materialAdoublePtr->H,dAta.materialAdoublePtr->detH
      ); CHKERRQ(ierr);
      dAta.materialAdoublePtr->invH.resize(3,3,false);
      ierr = dAta.materialAdoublePtr->iNvert(
        dAta.materialAdoublePtr->detH,dAta.materialAdoublePtr->H,
        dAta.materialAdoublePtr->invH
      ); CHKERRQ(ierr);
      noalias(dAta.materialAdoublePtr->F) = prod(
        dAta.materialAdoublePtr->h,dAta.materialAdoublePtr->invH
      );

    }

    ierr = dAta.materialAdoublePtr->setUserActiveVariables(nbActiveVariables); CHKERRQ(ierr);
    ierr = calculateEnergy(gg); CHKERRQ(ierr);

    trace_off();

  } catch (MoFEMException const &e) {
    SETERRQ(PETSC_COMM_SELF,e.errorCode,e.errorMessage);
  }

  PetscFunctionReturn(0);
}

PetscErrorCode NonlinearElasticElement::OpJacobianEnergy::playTag(const int gg) {
  PetscFunctionBegin;

  int r;

  if(gRadient) {
    commonData.jacEnergy[gg].resize(nbActiveVariables,false);
    int r = ::gradient(
      tAg,nbActiveVariables,&activeVariables[0],&commonData.jacEnergy[gg][0]
    );
    if(r<0) {
      // That means that energy function is not smooth and derivative
      // can not be calculated,
      SETERRQ(
        PETSC_COMM_SELF,
        MOFEM_OPERATION_UNSUCCESSFUL,
        "ADOL-C function evaluation with error"
      );
    }
  }

  if(hEssian) {
    SETERRQ(
      PETSC_COMM_SELF,
      MOFEM_NOT_IMPLEMENTED,
      "Not yet implemented"
    );
  }

  PetscFunctionReturn(0);
}

PetscErrorCode NonlinearElasticElement::OpJacobianEnergy::doWork(
  int row_side,EntityType row_type,DataForcesAndSurcesCore::EntData &row_data
) {
  PetscFunctionBegin;
  PetscFunctionBegin;

  //do it only once, no need to repeat this for edges,faces or tets
  if(row_type != MBVERTEX) PetscFunctionReturn(0);

  PetscErrorCode ierr;
  if(dAta.tEts.find(getNumeredEntFiniteElementPtr()->getEnt()) == dAta.tEts.end()) {
    PetscFunctionReturn(0);
  }

  int nb_dofs = row_data.getFieldData().size();
  if(nb_dofs==0) PetscFunctionReturn(0);
  dAta.materialAdoublePtr->commonDataPtr = &commonData;
  dAta.materialAdoublePtr->opPtr = this;

  try {

    int nb_gauss_pts = row_data.getN().size1();
    commonData.eNergy.resize(nb_gauss_pts);
    commonData.jacEnergy.resize(nb_gauss_pts);

    ptrh = &(commonData.gradAtGaussPts[commonData.spatialPositions]);
    if(aLe) {
      ptrH = &(commonData.gradAtGaussPts[commonData.meshPositions]);
    }

    for(int gg = 0;gg!=nb_gauss_pts;gg++) {

      dAta.materialAdoublePtr->gG = gg;

      // Record tag and calualte stress
      if(recordTagForIntegrationPoint(gg)) {
        ierr = recordTag(gg); CHKERRQ(ierr);
      }

      activeVariables.resize(nbActiveVariables,false);
      if(!aLe) {
        for(int dd1 = 0;dd1<3;dd1++) {
          for(int dd2 = 0;dd2<3;dd2++) {
            activeVariables(dd1*3+dd2) = (*ptrh)[gg](dd1,dd2);
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

      // Play tag and calualte stress or tannget
      ierr = playTag(gg); CHKERRQ(ierr);

    }

  } catch (MoFEMException const &e) {
    SETERRQ(PETSC_COMM_SELF,e.errorCode,e.errorMessage);
  }

  PetscFunctionReturn(0);
}


NonlinearElasticElement::OpRhsPiolaKirchhoff::OpRhsPiolaKirchhoff(const std::string field_name,BlockData &data,CommonData &common_data):
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
    VectorDofs& dofs = row_data.getFieldDofs();
    VectorDofs::iterator dit = dofs.begin();
    for(int ii = 0;dit!=dofs.end();dit++,ii++) {
      if(dAta.forcesOnlyOnEntitiesRow.find((*dit)->getEnt())==dAta.forcesOnlyOnEntitiesRow.end()) {
        iNdices[ii] = -1;
      }
    }
  }
  ierr = VecSetValues(
    getFEMethod()->snes_f,
    nb_dofs,
    indices_ptr,
    &nf[0],
    ADD_VALUES
  ); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode NonlinearElasticElement::OpRhsPiolaKirchhoff::doWork(
  int row_side,EntityType row_type,DataForcesAndSurcesCore::EntData &row_data
) {
  PetscFunctionBegin;

  PetscErrorCode ierr;
  if(dAta.tEts.find(getNumeredEntFiniteElementPtr()->getEnt()) == dAta.tEts.end()) {
    PetscFunctionReturn(0);
  }

  const int nb_dofs = row_data.getIndices().size();
  if(nb_dofs==0) PetscFunctionReturn(0);
  if((unsigned int)nb_dofs > 3*row_data.getN().size2()) {
    SETERRQ(PETSC_COMM_SELF,1,"data inconsistency");
  }
  const int nb_base_functions = row_data.getN().size2();
  const int nb_gauss_pts = row_data.getN().size1();

  try {

    nf.resize(nb_dofs,false);
    nf.clear();

    FTensor::Tensor1<double*,3> diff_base_functions = row_data.getFTensor1DiffN<3>();
    FTensor::Index<'i',3> i;
    FTensor::Index<'j',3> j;

    for(int gg = 0;gg!=nb_gauss_pts;gg++) {
      double val = getVolume()*getGaussPts()(3,gg);
      if((!aLe)&&getHoGaussPtsDetJac().size()>0) {
        val *= getHoGaussPtsDetJac()[gg]; ///< higher order geometry
      }
      const MatrixDouble& stress = commonData.sTress[gg];
      FTensor::Tensor2<const double *,3,3> t3(
        &stress(0,0),&stress(0,1),&stress(0,2),
        &stress(1,0),&stress(1,1),&stress(1,2),
        &stress(2,0),&stress(2,1),&stress(2,2)
      );
      FTensor::Tensor1<double*,3> rhs(&nf[0],&nf[1],&nf[2],3);
      int bb = 0;
      for(;bb!=nb_dofs/3;bb++) {
        rhs(i) += val*t3(i,j)*diff_base_functions(j);
        ++rhs;
        ++diff_base_functions;
      }
      for(;bb!=nb_base_functions;bb++) {
        ++diff_base_functions;
      }
    }

    //std::cerr << "nf : " << nf << std::endl;
    ierr = aSemble(row_side,row_type,row_data); CHKERRQ(ierr);

  } catch (MoFEMException const &e) {
    SETERRQ(PETSC_COMM_SELF,e.errorCode,e.errorMessage);
  }

  PetscFunctionReturn(0);
}

NonlinearElasticElement::OpEnergy::OpEnergy(
  const std::string field_name,BlockData &data,CommonData &common_data,Vec *v_ptr,bool field_disp
):
VolumeElementForcesAndSourcesCore::UserDataOperator(field_name,UserDataOperator::OPROW),
dAta(data),commonData(common_data),
Vptr(v_ptr),
fieldDisp(field_disp) {
}

PetscErrorCode NonlinearElasticElement::OpEnergy::doWork(
  int row_side,EntityType row_type,DataForcesAndSurcesCore::EntData &row_data
) {
  PetscFunctionBegin;

  PetscErrorCode ierr;
  if(row_type != MBVERTEX) PetscFunctionReturn(0);
  if(dAta.tEts.find(getNumeredEntFiniteElementPtr()->getEnt()) == dAta.tEts.end()) {
    PetscFunctionReturn(0);
  }

  try {

    std::vector<MatrixDouble > &F = (commonData.gradAtGaussPts[commonData.spatialPositions]);

    for(unsigned int gg = 0;gg<row_data.getN().size1();gg++) {
      double val = getVolume()*getGaussPts()(3,gg);
      if(getHoGaussPtsDetJac().size()>0) {
        val *= getHoGaussPtsDetJac()[gg]; ///< higher order geometry
      }

      dAta.materialDoublePtr->F.resize(3,3,false);
      noalias(dAta.materialDoublePtr->F) = F[gg];
      if(fieldDisp) {
        for(int dd = 0;dd<3;dd++) {
          dAta.materialDoublePtr->F(dd,dd) += 1;
        }
      }

      int nb_active_variables = 0;
      ierr = dAta.materialDoublePtr->setUserActiveVariables(nb_active_variables); CHKERRQ(ierr);
      ierr = dAta.materialDoublePtr->calculateElasticEnergy(dAta,getNumeredEntFiniteElementPtr()); CHKERRQ(ierr);
      ierr = VecSetValue(*Vptr,0,val*dAta.materialDoublePtr->eNergy,ADD_VALUES); CHKERRQ(ierr);

    }

  } catch (MoFEMException const &e) {
    SETERRQ(PETSC_COMM_SELF,e.errorCode,e.errorMessage);
  }

  PetscFunctionReturn(0);
}


NonlinearElasticElement::OpLhsPiolaKirchhoff_dx::OpLhsPiolaKirchhoff_dx(
  const std::string vel_field,const std::string field_name,BlockData &data,CommonData &common_data
):
VolumeElementForcesAndSourcesCore::UserDataOperator(vel_field,field_name,UserDataOperator::OPROWCOL),
dAta(data),
commonData(common_data),
aLe(false) {
}

PetscErrorCode NonlinearElasticElement::OpLhsPiolaKirchhoff_dx::getJac(
  DataForcesAndSurcesCore::EntData &col_data,int gg
) {
  PetscFunctionBegin;
  jac.clear();
  FTensor::Index<'i',3> i;
  FTensor::Index<'j',3> j;
  FTensor::Index<'k',3> k;
  MatrixDouble& jac_stress = commonData.jacStress[gg];
  int nb_col = col_data.getFieldData().size();
  double *diff_ptr = const_cast<double*>(&(col_data.getDiffN(gg,nb_col/3)(0,0)));
  // First two indices 'i','j' derivatives of 1st Piola-stress, third index 'k' is
  // displacement component
  FTensor::Tensor3<double*,3,3,3> t3_1(
    &jac_stress(3*0+0,0),&jac_stress(3*0+0,1),&jac_stress(3*0+0,2),
    &jac_stress(3*0+1,0),&jac_stress(3*0+1,1),&jac_stress(3*0+1,2),
    &jac_stress(3*0+2,0),&jac_stress(3*0+2,1),&jac_stress(3*0+2,2),
    &jac_stress(3*1+0,0),&jac_stress(3*1+0,1),&jac_stress(3*1+0,2),
    &jac_stress(3*1+1,0),&jac_stress(3*1+1,1),&jac_stress(3*1+1,2),
    &jac_stress(3*1+2,0),&jac_stress(3*1+2,1),&jac_stress(3*1+2,2),
    &jac_stress(3*2+0,0),&jac_stress(3*2+0,1),&jac_stress(3*2+0,2),
    &jac_stress(3*2+1,0),&jac_stress(3*2+1,1),&jac_stress(3*2+1,2),
    &jac_stress(3*2+2,0),&jac_stress(3*2+2,1),&jac_stress(3*2+2,2),3
  );
  for(int rr = 0;rr!=3;rr++) {
    // Derivate of 1st Piola-stress multiplied by gradient of defamation for
    // base function (dd) and displacement component (rr)
    FTensor::Tensor2<double*,3,3> t2_1(
      &jac(0,rr),&jac(1,rr),&jac(2,rr),
      &jac(3,rr),&jac(4,rr),&jac(5,rr),
      &jac(6,rr),&jac(7,rr),&jac(8,rr),3
    );
    FTensor::Tensor1<double*,3> diff(diff_ptr,&diff_ptr[1],&diff_ptr[2],3);
    for(int dd = 0;dd!=nb_col/3;dd++) {
      t2_1(i,j) += t3_1(i,j,k)*diff(k);
      ++t2_1;
      ++diff;
    }
    ++t3_1;
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
    VectorDofs& dofs = row_data.getFieldDofs();
    VectorDofs::iterator dit = dofs.begin();
    for(int ii = 0;dit!=dofs.end();dit++,ii++) {
      if(dAta.forcesOnlyOnEntitiesRow.find((*dit)->getEnt())==dAta.forcesOnlyOnEntitiesRow.end()) {
        rowIndices[ii] = -1;
      }
    }
  }

  if(!dAta.forcesOnlyOnEntitiesCol.empty()) {
    colIndices.resize(nb_col,false);
    noalias(colIndices) = col_data.getIndices();
    col_indices_ptr = &colIndices[0];
    VectorDofs& dofs = col_data.getFieldDofs();
    VectorDofs::iterator dit = dofs.begin();
    for(int ii = 0;dit!=dofs.end();dit++,ii++) {
      if(dAta.forcesOnlyOnEntitiesCol.find((*dit)->getEnt())==dAta.forcesOnlyOnEntitiesCol.end()) {
        colIndices[ii] = -1;
      }
    }
  }

  ierr = MatSetValues(
    getFEMethod()->snes_B,
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
      VectorDofs& dofs = row_data.getFieldDofs();
      VectorDofs::iterator dit = dofs.begin();
      for(int ii = 0;dit!=dofs.end();dit++,ii++) {
        if(dAta.forcesOnlyOnEntitiesCol.find((*dit)->getEnt())==dAta.forcesOnlyOnEntitiesCol.end()) {
          rowIndices[ii] = -1;
        }
      }
    }

    if(!dAta.forcesOnlyOnEntitiesRow.empty()) {
      colIndices.resize(nb_col,false);
      noalias(colIndices) = col_data.getIndices();
      col_indices_ptr = &colIndices[0];
      VectorDofs& dofs = col_data.getFieldDofs();
      VectorDofs::iterator dit = dofs.begin();
      for(int ii = 0;dit!=dofs.end();dit++,ii++) {
        if(dAta.forcesOnlyOnEntitiesRow.find((*dit)->getEnt())==dAta.forcesOnlyOnEntitiesRow.end()) {
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

  if(dAta.tEts.find(getNumeredEntFiniteElementPtr()->getEnt()) == dAta.tEts.end()) {
    PetscFunctionReturn(0);
  }

  // const int nb_base_functions = row_data.getN().size2();
  const int nb_gauss_pts = row_data.getN().size1();

  FTensor::Index<'i',3> i;
  FTensor::Index<'j',3> j;
  FTensor::Index<'m',3> m;

  try {

    k.resize(nb_row,nb_col,false);
    k.clear();
    jac.resize(9,nb_col,false);

    for(unsigned int gg = 0;gg!=nb_gauss_pts;gg++) {
      ierr = getJac(col_data,gg); CHKERRQ(ierr);
      double val = getVolume()*getGaussPts()(3,gg);
      if((!aLe)&&(getHoGaussPtsDetJac().size()>0)) {
        val *= getHoGaussPtsDetJac()[gg]; ///< higher order geometry
      }
      FTensor::Tensor3<double*,3,3,3> t3_1(
        &jac(3*0+0,0),&jac(3*0+0,1),&jac(3*0+0,2),
        &jac(3*0+1,0),&jac(3*0+1,1),&jac(3*0+1,2),
        &jac(3*0+2,0),&jac(3*0+2,1),&jac(3*0+2,2),
        &jac(3*1+0,0),&jac(3*1+0,1),&jac(3*1+0,2),
        &jac(3*1+1,0),&jac(3*1+1,1),&jac(3*1+1,2),
        &jac(3*1+2,0),&jac(3*1+2,1),&jac(3*1+2,2),
        &jac(3*2+0,0),&jac(3*2+0,1),&jac(3*2+0,2),
        &jac(3*2+1,0),&jac(3*2+1,1),&jac(3*2+1,2),
        &jac(3*2+2,0),&jac(3*2+2,1),&jac(3*2+2,2),3
      );
      for(int rr = 0;rr!=nb_col/3;rr++) {
        FTensor::Tensor1<double*,3> diff_base_functions = row_data.getFTensor1DiffN<3>(gg,0);
        int bb = 0;
        for(;bb!=nb_row/3;bb++) {
          FTensor::Tensor2<double*,3,3> lhs(
            &k(3*bb+0,3*rr+0),&k(3*bb+0,3*rr+1),&k(3*bb+0,3*rr+2),
            &k(3*bb+1,3*rr+0),&k(3*bb+1,3*rr+1),&k(3*bb+1,3*rr+2),
            &k(3*bb+2,3*rr+0),&k(3*bb+2,3*rr+1),&k(3*bb+2,3*rr+2)
          );
          lhs(i,j) += val*t3_1(i,m,j)*diff_base_functions(m);
          ++diff_base_functions;
        }
        ++t3_1;
      }
    }

    //std::cerr << "N " << getNumeredEntFiniteElementPtr()->getRefEnt() << std::endl << k << std::endl;
    ierr = aSemble(row_side,col_side,row_type,col_type,row_data,col_data); CHKERRQ(ierr);

  } catch (MoFEMException const &e) {
    SETERRQ(PETSC_COMM_SELF,e.errorCode,e.errorMessage);
  }

  PetscFunctionReturn(0);
}

NonlinearElasticElement::OpLhsPiolaKirchhoff_dX::OpLhsPiolaKirchhoff_dX(
  const std::string vel_field,const std::string field_name,BlockData &data,CommonData &common_data):
  OpLhsPiolaKirchhoff_dx(vel_field,field_name,data,common_data)
  { sYmm = false; }

PetscErrorCode NonlinearElasticElement::OpLhsPiolaKirchhoff_dX::getJac(DataForcesAndSurcesCore::EntData &col_data,int gg) {
  PetscFunctionBegin;
  jac.clear();
  int nb_col = col_data.getFieldData().size();
  const MatrixAdaptor diffN = col_data.getDiffN(gg,nb_col/3);
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
    VectorDofs& dofs = row_data.getFieldDofs();
    VectorDofs::iterator dit = dofs.begin();
    for(int ii = 0;dit!=dofs.end();dit++,ii++) {
      if(dAta.forcesOnlyOnEntitiesRow.find((*dit)->getEnt())==dAta.forcesOnlyOnEntitiesRow.end()) {
        rowIndices[ii] = -1;
      }
    }
  }

  int *col_indices_ptr = &col_data.getIndices()[0];
  if(!dAta.forcesOnlyOnEntitiesCol.empty()) {
    colIndices.resize(nb_col,false);
    noalias(colIndices) = col_data.getIndices();
    col_indices_ptr = &colIndices[0];
    VectorDofs& dofs = col_data.getFieldDofs();
    VectorDofs::iterator dit = dofs.begin();
    for(int ii = 0;dit!=dofs.end();dit++,ii++) {
      if(dAta.forcesOnlyOnEntitiesCol.find((*dit)->getEnt())==dAta.forcesOnlyOnEntitiesCol.end()) {
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
  const std::string field_name,
  BlockData &data,
  CommonData &common_data,
  int tag,
  bool jacobian,
  bool ale
):
OpJacobianPiolaKirchhoffStress(field_name,data,common_data,tag,jacobian,ale,false) {
}

PetscErrorCode NonlinearElasticElement::OpJacobianEshelbyStress::calculateStress(const int gg) {
  PetscFunctionBegin;
  try {
    PetscErrorCode ierr;
    ierr = dAta.materialAdoublePtr->calculateSiGma_EshelbyStress(dAta,getNumeredEntFiniteElementPtr()); CHKERRQ(ierr);
    if(aLe) {
      dAta.materialAdoublePtr->SiGma =
      dAta.materialAdoublePtr->detH*prod(dAta.materialAdoublePtr->SiGma,trans(dAta.materialAdoublePtr->invH));
    }
    commonData.sTress[gg].resize(3,3,false);
    for(int dd1 = 0;dd1<3;dd1++) {
      for(int dd2 = 0;dd2<3;dd2++) {
        dAta.materialAdoublePtr->SiGma(dd1,dd2) >>= (commonData.sTress[gg])(dd1,dd2);
      }
    }
  } catch (MoFEMException const &e) {
    SETERRQ(PETSC_COMM_SELF,e.errorCode,e.errorMessage);
  }
  PetscFunctionReturn(0);
}

NonlinearElasticElement::OpRhsEshelbyStrees::OpRhsEshelbyStrees(
  const std::string field_name,BlockData &data,CommonData &common_data
):
OpRhsPiolaKirchhoff(field_name,data,common_data)
{}

NonlinearElasticElement::OpLhsEshelby_dx::OpLhsEshelby_dx(
  const std::string vel_field,const std::string field_name,BlockData &data,CommonData &common_data
):
OpLhsPiolaKirchhoff_dX(vel_field,field_name,data,common_data) {}

PetscErrorCode NonlinearElasticElement::OpLhsEshelby_dx::getJac(DataForcesAndSurcesCore::EntData &col_data,int gg) {
  PetscFunctionBegin;
  jac.clear();
  FTensor::Index<'i',3> i;
  FTensor::Index<'j',3> j;
  FTensor::Index<'k',3> k;
  MatrixDouble &jac_stress = commonData.jacStress[gg];
  int nb_col = col_data.getFieldData().size();
  double *diff_ptr = const_cast<double*>(&(col_data.getDiffN(gg,nb_col/3)(0,0)));
  FTensor::Tensor1<double*,3> diff(diff_ptr,&diff_ptr[1],&diff_ptr[2],3);
  for(int dd = 0;dd!=nb_col/3;dd++) {
    for(int rr = 0;rr!=3;rr++) {
      // Derivate of 1st Piola-stress multiplied by gradient of defamation for
      // base function (dd) and displacement component (rr)
      FTensor::Tensor2<double*,3,3> t2_1(
        &jac(0,3*dd+rr),&jac(1,3*dd+rr),&jac(2,3*dd+rr),
        &jac(3,3*dd+rr),&jac(4,3*dd+rr),&jac(5,3*dd+rr),
        &jac(6,3*dd+rr),&jac(7,3*dd+rr),&jac(8,3*dd+rr)
      );
      // First two indices 'i','j' derivatives of 1st Piola-stress, third index 'k' is
      // displacement component
      FTensor::Tensor3<double*,3,3,3> t3_1(
        &jac_stress(3*0+0,3*rr+0),&jac_stress(3*0+0,3*rr+1),&jac_stress(3*0+0,3*rr+2),
        &jac_stress(3*0+1,3*rr+0),&jac_stress(3*0+1,3*rr+1),&jac_stress(3*0+1,3*rr+2),
        &jac_stress(3*0+2,3*rr+0),&jac_stress(3*0+2,3*rr+1),&jac_stress(3*0+2,3*rr+2),
        &jac_stress(3*1+0,3*rr+0),&jac_stress(3*1+0,3*rr+1),&jac_stress(3*1+0,3*rr+2),
        &jac_stress(3*1+1,3*rr+0),&jac_stress(3*1+1,3*rr+1),&jac_stress(3*1+1,3*rr+2),
        &jac_stress(3*1+2,3*rr+0),&jac_stress(3*1+2,3*rr+1),&jac_stress(3*1+2,3*rr+2),
        &jac_stress(3*2+0,3*rr+0),&jac_stress(3*2+0,3*rr+1),&jac_stress(3*2+0,3*rr+2),
        &jac_stress(3*2+1,3*rr+0),&jac_stress(3*2+1,3*rr+1),&jac_stress(3*2+1,3*rr+2),
        &jac_stress(3*2+2,3*rr+0),&jac_stress(3*2+2,3*rr+1),&jac_stress(3*2+2,3*rr+2)
      );
      t2_1(i,j) += t3_1(i,j,k)*diff(k);
    }
    ++diff;
  }
  // int nb_col = col_data.getFieldData().size();
  // const MatrixAdaptor diffN = col_data.getDiffN(gg,nb_col/3);
  // for(int dd = 0;dd<nb_col/3;dd++) {
  //   for(int rr = 0;rr<3;rr++) {
  //     for(int ii = 0;ii<9;ii++) {
  //       for(int jj = 0;jj<3;jj++) {
  //         jac(ii,3*dd+rr) += commonData.jacStress[gg](ii,3*rr+jj)*diffN(dd,jj);
  //       }
  //     }
  //   }
  // }
  PetscFunctionReturn(0);
}

NonlinearElasticElement::OpLhsEshelby_dX::OpLhsEshelby_dX(
  const std::string vel_field,const std::string field_name,BlockData &data,CommonData &common_data
):
OpLhsPiolaKirchhoff_dx(vel_field,field_name,data,common_data)
{}

PetscErrorCode NonlinearElasticElement::OpLhsEshelby_dX::getJac(DataForcesAndSurcesCore::EntData &col_data,int gg) {
  PetscFunctionBegin;
  jac.clear();
  FTensor::Index<'i',3> i;
  FTensor::Index<'j',3> j;
  FTensor::Index<'k',3> k;
  MatrixDouble &jac_stress = commonData.jacStress[gg];
  int nb_col = col_data.getFieldData().size();
  double *diff_ptr = const_cast<double*>(&(col_data.getDiffN(gg,nb_col/3)(0,0)));
  // First two indices 'i','j' derivatives of 1st Piola-stress, third index 'k' is
  // displacement component
  FTensor::Tensor3<double*,3,3,3> t3_1(
    &jac_stress(3*0+0,9),&jac_stress(3*0+0,9+1),&jac_stress(3*0+0,9+2),
    &jac_stress(3*0+1,9),&jac_stress(3*0+1,9+1),&jac_stress(3*0+1,9+2),
    &jac_stress(3*0+2,9),&jac_stress(3*0+2,9+1),&jac_stress(3*0+2,9+2),
    &jac_stress(3*1+0,9),&jac_stress(3*1+0,9+1),&jac_stress(3*1+0,9+2),
    &jac_stress(3*1+1,9),&jac_stress(3*1+1,9+1),&jac_stress(3*1+1,9+2),
    &jac_stress(3*1+2,9),&jac_stress(3*1+2,9+1),&jac_stress(3*1+2,9+2),
    &jac_stress(3*2+0,9),&jac_stress(3*2+0,9+1),&jac_stress(3*2+0,9+2),
    &jac_stress(3*2+1,9),&jac_stress(3*2+1,9+1),&jac_stress(3*2+1,9+2),
    &jac_stress(3*2+2,9),&jac_stress(3*2+2,9+1),&jac_stress(3*2+2,9+2),3
  );
  for(int rr = 0;rr!=3;rr++) {
    // Derivate of 1st Piola-stress multiplied by gradient of defamation for
    // base function (dd) and displacement component (rr)
    FTensor::Tensor2<double*,3,3> t2_1(
      &jac(0,rr),&jac(1,rr),&jac(2,rr),
      &jac(3,rr),&jac(4,rr),&jac(5,rr),
      &jac(6,rr),&jac(7,rr),&jac(8,rr),3
    );
    FTensor::Tensor1<double*,3> diff(diff_ptr,&diff_ptr[1],&diff_ptr[2],3);
    for(int dd = 0;dd!=nb_col/3;dd++) {
      t2_1(i,j) += t3_1(i,j,k)*diff(k);
      ++t2_1;
      ++diff;
    }
    ++t3_1;
  }
  // MatrixDouble &jac_stress = commonData.jacStress[gg];
  // int nb_col = col_data.getFieldData().size();
  // double *diff_ptr = const_cast<double*>(&(col_data.getDiffN(gg,nb_col/3)(0,0)));
  // FTensor::Tensor1<double*,3> diff(diff_ptr,&diff_ptr[1],&diff_ptr[2],3);
  // for(int dd = 0;dd!=nb_col/3;dd++) {
  //   for(int rr = 0;rr!=3;rr++) {
  //     // Derivate of 1st Piola-stress multiplied by gradient of defamation for
  //     // base function (dd) and displacement component (rr)
  //     FTensor::Tensor2<double*,3,3> t2_1(
  //       &jac(0,3*dd+rr),&jac(1,3*dd+rr),&jac(2,3*dd+rr),
  //       &jac(3,3*dd+rr),&jac(4,3*dd+rr),&jac(5,3*dd+rr),
  //       &jac(6,3*dd+rr),&jac(7,3*dd+rr),&jac(8,3*dd+rr)
  //     );
  //     // First two indices 'i','j' derivatives of 1st Piola-stress, third index 'k' is
  //     // displacement component
  //     FTensor::Tensor3<double*,3,3,3> t3_1(
  //       &jac_stress(3*0+0,9+3*rr+0),&jac_stress(3*0+0,9+3*rr+1),&jac_stress(3*0+0,9+3*rr+2),
  //       &jac_stress(3*0+1,9+3*rr+0),&jac_stress(3*0+1,9+3*rr+1),&jac_stress(3*0+1,9+3*rr+2),
  //       &jac_stress(3*0+2,9+3*rr+0),&jac_stress(3*0+2,9+3*rr+1),&jac_stress(3*0+2,9+3*rr+2),
  //       &jac_stress(3*1+0,9+3*rr+0),&jac_stress(3*1+0,9+3*rr+1),&jac_stress(3*1+0,9+3*rr+2),
  //       &jac_stress(3*1+1,9+3*rr+0),&jac_stress(3*1+1,9+3*rr+1),&jac_stress(3*1+1,9+3*rr+2),
  //       &jac_stress(3*1+2,9+3*rr+0),&jac_stress(3*1+2,9+3*rr+1),&jac_stress(3*1+2,9+3*rr+2),
  //       &jac_stress(3*2+0,9+3*rr+0),&jac_stress(3*2+0,9+3*rr+1),&jac_stress(3*2+0,9+3*rr+2),
  //       &jac_stress(3*2+1,9+3*rr+0),&jac_stress(3*2+1,9+3*rr+1),&jac_stress(3*2+1,9+3*rr+2),
  //       &jac_stress(3*2+2,9+3*rr+0),&jac_stress(3*2+2,9+3*rr+1),&jac_stress(3*2+2,9+3*rr+2)
  //     );
  //     t2_1(i,j) += t3_1(i,j,k)*diff(k);
  //   }
  //   ++diff;
  // }
  // int nb_col = col_data.getFieldData().size();
  // const MatrixAdaptor diffN = col_data.getDiffN(gg,nb_col/3);
  // for(int dd = 0;dd<nb_col/3;dd++) {
  //   for(int rr = 0;rr<3;rr++) {
  //     for(int ii = 0;ii<9;ii++) {
  //       for(int jj = 0;jj<3;jj++) {
  //         jac(ii,3*dd+rr) += commonData.jacStress[gg](ii,9+3*rr+jj)*diffN(dd,jj);
  //       }
  //     }
  //   }
  // }
  PetscFunctionReturn(0);
}

PetscErrorCode NonlinearElasticElement::setBlocks(
  boost::shared_ptr<FunctionsToCalculatePiolaKirchhoffI<double> > materialDoublePtr,
  boost::shared_ptr<FunctionsToCalculatePiolaKirchhoffI<adouble> > materialAdoublePtr
) {
  PetscFunctionBegin;
  ErrorCode rval;
  PetscErrorCode ierr;

  if(!materialDoublePtr) {
    SETERRQ(
      mField.get_comm(),
      MOFEM_DATA_INCONSISTENCY,
      "Pointer for materialDoublePtr not allocated"
    );
  }
  if(!materialAdoublePtr) {
    SETERRQ(
      mField.get_comm(),
      MOFEM_DATA_INCONSISTENCY,
      "Pointer for materialAdoublePtr not allocated"
    );
  }

  for(_IT_CUBITMESHSETS_BY_BCDATA_TYPE_FOR_LOOP_(mField,BLOCKSET|MAT_ELASTICSET,it)) {
    Mat_Elastic mydata;
    ierr = it->getAttributeDataStructure(mydata); CHKERRQ(ierr);
    int id = it->getMeshsetId();
    EntityHandle meshset = it->getMeshset();
    rval = mField.get_moab().get_entities_by_type(meshset,MBTET,setOfBlocks[id].tEts,true); CHKERRQ_MOAB(rval);
    setOfBlocks[id].iD = id;
    setOfBlocks[id].E = mydata.data.Young;
    setOfBlocks[id].PoissonRatio = mydata.data.Poisson;
    setOfBlocks[id].materialDoublePtr = materialDoublePtr;
    setOfBlocks[id].materialAdoublePtr = materialAdoublePtr;
    //std::cerr << setOfBlocks[id].tEts << std::endl;
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

  std::map<int,BlockData>::iterator sit = setOfBlocks.begin();
  for(;sit!=setOfBlocks.end();sit++) {
    ierr = mField.add_ents_to_finite_element_by_type(sit->second.tEts,MBTET,element_name); CHKERRQ(ierr);
  }

  PetscFunctionReturn(0);
}

PetscErrorCode NonlinearElasticElement::setOperators(
  string spatial_position_field_name,
  string material_position_field_name,
  bool ale,bool field_disp
) {
  PetscFunctionBegin;

  commonData.spatialPositions = spatial_position_field_name;
  commonData.meshPositions = material_position_field_name;

  //Rhs
  feRhs.getOpPtrVector().push_back(new OpGetCommonDataAtGaussPts(spatial_position_field_name,commonData));
  if(mField.check_field(material_position_field_name)) {
    feRhs.getOpPtrVector().push_back(new OpGetCommonDataAtGaussPts(material_position_field_name,commonData));
  }
  std::map<int,BlockData>::iterator sit = setOfBlocks.begin();
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
