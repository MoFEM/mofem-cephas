/** \file VolumeElementForcesAndSourcesCore.cpp

\brief Implementation of volume element

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

#include <Includes.hpp>
#include <version.h>
#include <definitions.h>
#include <Common.hpp>

#include <base_functions.h>
#include <h1_hdiv_hcurl_l2.h>
#include <fem_tools.h>

#include <UnknownInterface.hpp>

#include <MaterialBlocks.hpp>
#include <BCData.hpp>
#include <TagMultiIndices.hpp>
#include <CoordSysMultiIndices.hpp>
#include <FieldMultiIndices.hpp>
#include <EntsMultiIndices.hpp>
#include <DofsMultiIndices.hpp>
#include <FEMMultiIndices.hpp>
#include <ProblemsMultiIndices.hpp>
#include <AdjacencyMultiIndices.hpp>
#include <BCMultiIndices.hpp>
#include <CoreDataStructures.hpp>
#include <SeriesMultiIndices.hpp>

#include <LoopMethods.hpp>
#include <Interface.hpp>
#include <MeshRefinment.hpp>
#include <PrismInterface.hpp>
#include <SeriesRecorder.hpp>
#include <Core.hpp>

#include <BaseFunction.hpp>
#include <LegendrePolynomial.hpp>
#include <LobattoPolynomial.hpp>

#include <FTensor.hpp>
#include <DataStructures.hpp>
#include <EntPolynomialBaseCtx.hpp>
#include <TetPolynomialBase.hpp> // Base functions on tet
#include <DataOperators.hpp>
#include <ElementsOnEntities.hpp>
#include <VolumeElementForcesAndSourcesCore.hpp>
#include <FaceElementForcesAndSourcesCore.hpp>

#ifdef __cplusplus
extern "C" {
#endif
  #include <cblas.h>
  #include <lapack_wrap.h>
  // #include <gm_rule.h>
  #include <quad.h>
#ifdef __cplusplus
}
#endif

namespace MoFEM {

VolumeElementForcesAndSourcesCore::VolumeElementForcesAndSourcesCore(
  Interface &m_field,const EntityType type
):
ForcesAndSurcesCore(m_field),
coords(12),
jAc(3,3),
invJac(3,3),
dataH1(type),
derivedDataH1(dataH1),
dataL2(type),
derivedDataL2(dataL2),
dataHdiv(type),
derivedDataHdiv(dataHdiv),
dataHcurl(type),
derivedDataHcurl(dataHcurl),
dataNoField(type),
dataNoFieldCol(type),
opSetInvJacH1(invJac),
opContravariantPiolaTransform(vOlume,jAc),
opCovariantPiolaTransform(invJac),
opSetInvJacHdivAndHcurl(invJac),
meshPositionsFieldName("MESH_NODE_POSITIONS"),
opHOatGaussPoints(hoCoordsAtGaussPts,hoGaussPtsJac,3,3),
opSetHoInvJacH1(hoGaussPtsInvJac),
opHoContravariantTransform(hoGaussPtsDetJac,hoGaussPtsJac),
opHoCovariantTransform(hoGaussPtsInvJac),
opSetHoInvJacHdivAndHcurl(hoGaussPtsInvJac),
tJac(
  &jAc(0,0),&jAc(0,1),&jAc(0,2),
  &jAc(1,0),&jAc(1,1),&jAc(1,2),
  &jAc(2,0),&jAc(2,1),&jAc(2,2)
),
tInvJac(
  &invJac(0,0),&invJac(0,1),&invJac(0,2),
  &invJac(1,0),&invJac(1,1),&invJac(1,2),
  &invJac(2,0),&invJac(2,1),&invJac(2,2)
) {
}

PetscErrorCode VolumeElementForcesAndSourcesCore::setIntegartionPts() {
  PetscFunctionBegin;
  int order_data = getMaxDataOrder();
  int order_row = getMaxRowOrder();
  int order_col = getMaxColOrder();
  int rule = getRule(order_row,order_col,order_data);
  // std::cerr << order_data << " " << order_row << " " << order_col << " " << rule << std::endl;
  if(rule >= 0) {
    if(rule<QUAD_3D_TABLE_SIZE) {
      if(QUAD_3D_TABLE[rule]->dim!=3) {
        SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"wrong dimension");
      }
      if(QUAD_3D_TABLE[rule]->order<rule) {
        SETERRQ2(
          PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"wrong order %d != %d",
          QUAD_3D_TABLE[rule]->order,rule
        );
      }
      nbGaussPts = QUAD_3D_TABLE[rule]->npoints;
      gaussPts.resize(4,nbGaussPts,false);
      cblas_dcopy(
        nbGaussPts,&QUAD_3D_TABLE[rule]->points[1],4,&gaussPts(0,0),1
      );
      cblas_dcopy(
        nbGaussPts,&QUAD_3D_TABLE[rule]->points[2],4,&gaussPts(1,0),1
      );
      cblas_dcopy(
        nbGaussPts,&QUAD_3D_TABLE[rule]->points[3],4,&gaussPts(2,0),1
      );
      cblas_dcopy(
        nbGaussPts,QUAD_3D_TABLE[rule]->weights,1,&gaussPts(3,0),1
      );
      dataH1.dataOnEntities[MBVERTEX][0].getN(NOBASE).resize(nbGaussPts,4,false);
      double *shape_ptr = dataH1.dataOnEntities[MBVERTEX][0].getN(NOBASE).data().begin();
      cblas_dcopy(
        4*nbGaussPts,QUAD_3D_TABLE[rule]->points,1,shape_ptr,1
      );
    } else {
      SETERRQ2(
        PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"rule > quadrature order %d < %d",
        rule,QUAD_3D_TABLE_SIZE
      );
      nbGaussPts = 0;
    }
  } else {
    ierr = setGaussPts(order_row,order_col,order_data); CHKERRQ(ierr);
    nbGaussPts = gaussPts.size2();
    dataH1.dataOnEntities[MBVERTEX][0].getN(NOBASE).resize(nbGaussPts,4,false);
    if(nbGaussPts>0) {
      ierr = ShapeMBTET(
        &*dataH1.dataOnEntities[MBVERTEX][0].getN(NOBASE).data().begin(),
        &gaussPts(0,0),&gaussPts(1,0),&gaussPts(2,0),nbGaussPts
      ); CHKERRQ(ierr);
    }
  }
  PetscFunctionReturn(0);
}

PetscErrorCode VolumeElementForcesAndSourcesCore::calculateVolumeAndJacobian() {
  PetscFunctionBegin;
  EntityHandle ent = numeredEntFiniteElementPtr->getEnt();
  rval = mField.get_moab().get_connectivity(ent,conn,num_nodes,true); CHKERRQ_MOAB(rval);
  rval = mField.get_moab().get_coords(conn,num_nodes,&*coords.data().begin()); CHKERRQ_MOAB(rval);
  double diff_n[12];
  ierr = ShapeDiffMBTET(diff_n); CHKERRQ(ierr);
  FTensor::Tensor1<double*,3> t_diff_n(&diff_n[0],&diff_n[1],&diff_n[2],3);
  FTensor::Tensor1<double*,3> t_coords(&coords[0],&coords[1],&coords[2],3);
  FTensor::Index<'i',3> i;
  FTensor::Index<'j',3> j;
  jAc.clear();
  for(int nn = 0;nn!=4;nn++) {
    tJac(i,j) += t_coords(i)*t_diff_n(j);
    ++t_coords;
    ++t_diff_n;
  }
  ierr = determinantTensor3by3(tJac,vOlume); CHKERRQ(ierr);
  ierr = invertTensor3by3(tJac,vOlume,tInvJac); CHKERRQ(ierr);
  vOlume *= G_TET_W1[0]/6.;
  PetscFunctionReturn(0);
}

PetscErrorCode VolumeElementForcesAndSourcesCore::calculateCoordinatesAtGaussPts() {
  PetscFunctionBegin;
  try {
    // Get coords at Gauss points
    FTensor::Index<'i',3> i;
    double *shape_functions_ptr = &*dataH1.dataOnEntities[MBVERTEX][0].getN(NOBASE).data().begin();
    coordsAtGaussPts.resize(nbGaussPts,3,false);
    coordsAtGaussPts.clear();
    FTensor::Tensor1<double*,3> t_coords_at_gauss_ptr(
      &coordsAtGaussPts(0,0),&coordsAtGaussPts(0,1),&coordsAtGaussPts(0,2),3
    );
    FTensor::Tensor0<double*> t_shape_functions(shape_functions_ptr);
    for(int gg = 0;gg<nbGaussPts;gg++) {
      FTensor::Tensor1<double*,3> t_coords(&coords[0],&coords[1],&coords[2],3);
      for(int bb = 0;bb<4;bb++) {
        t_coords_at_gauss_ptr(i) += t_coords(i)*t_shape_functions;
        ++t_coords;
        ++t_shape_functions;
      };
      ++t_coords_at_gauss_ptr;
    }

  } catch (std::exception& ex) {
    std::ostringstream ss;
    ss << "thorw in method: " << ex.what() << " at line " << __LINE__ << " in file " << __FILE__;
    SETERRQ(PETSC_COMM_SELF,MOFEM_STD_EXCEPTION_THROW,ss.str().c_str());
  }
  PetscFunctionReturn(0);
}

PetscErrorCode VolumeElementForcesAndSourcesCore::getSpaceBaseAndOrderOnElement() {
  PetscFunctionBegin;
  try {
    ierr = getSpacesAndBaseOnEntities(dataH1); CHKERRQ(ierr);
    ierr = getFaceTriNodes(dataH1); CHKERRQ(ierr);
    //H1
    if((dataH1.spacesOnEntities[MBEDGE]).test(H1)) {
      ierr = getEdgesSense(dataH1); CHKERRQ(ierr);
      ierr = getEdgesDataOrder(dataH1,H1); CHKERRQ(ierr);
    }
    if((dataH1.spacesOnEntities[MBTRI]).test(H1)) {
      ierr = getTrisSense(dataH1); CHKERRQ(ierr);
      ierr = getTrisDataOrder(dataH1,H1); CHKERRQ(ierr);
    }
    if((dataH1.spacesOnEntities[MBTET]).test(H1)) {
      ierr = getTetDataOrder(dataH1,H1); CHKERRQ(ierr);
    }
    //Hcurl
    if((dataH1.spacesOnEntities[MBEDGE]).test(HCURL)) {
      ierr = getEdgesSense(dataHcurl); CHKERRQ(ierr);
      ierr = getEdgesDataOrder(dataHcurl,HCURL); CHKERRQ(ierr);
      dataHcurl.spacesOnEntities[MBEDGE].set(HCURL);
    }
    if((dataH1.spacesOnEntities[MBTRI]).test(HCURL)) {
      ierr = getTrisSense(dataHcurl); CHKERRQ(ierr);
      ierr = getFaceTriNodes(dataHcurl); CHKERRQ(ierr);
      ierr = getTrisDataOrder(dataHcurl,HCURL); CHKERRQ(ierr);
      dataHcurl.spacesOnEntities[MBTRI].set(HCURL);
    }
    if((dataH1.spacesOnEntities[MBTET]).test(HCURL)) {
      ierr = getTetDataOrder(dataHcurl,HCURL); CHKERRQ(ierr);
      dataHcurl.spacesOnEntities[MBTET].set(HCURL);
    }
    //Hdiv
    if((dataH1.spacesOnEntities[MBTRI]).test(HDIV)) {
      ierr = getTrisSense(dataHdiv); CHKERRQ(ierr);
      ierr = getFaceTriNodes(dataHdiv); CHKERRQ(ierr);
      ierr = getTrisDataOrder(dataHdiv,HDIV); CHKERRQ(ierr);
      dataHdiv.spacesOnEntities[MBTRI].set(HDIV);
    }
    if((dataH1.spacesOnEntities[MBTET]).test(HDIV)) {
      ierr = getTetDataOrder(dataHdiv,HDIV); CHKERRQ(ierr);
      dataHdiv.spacesOnEntities[MBTET].set(HDIV);
    }
    //L2
    if((dataH1.spacesOnEntities[MBTET]).test(L2)) {
      ierr = getTetDataOrder(dataL2,L2); CHKERRQ(ierr);
      dataL2.spacesOnEntities[MBTET].set(L2);
    }
  } catch (std::exception& ex) {
    std::ostringstream ss;
    ss << "thorw in method: " << ex.what() << " at line " << __LINE__ << " in file " << __FILE__;
    SETERRQ(PETSC_COMM_SELF,MOFEM_STD_EXCEPTION_THROW,ss.str().c_str());
  }
  PetscFunctionReturn(0);
}

PetscErrorCode VolumeElementForcesAndSourcesCore::calculateBaseFunctionsOnElement() {
  PetscFunctionBegin;
  try {
    /// Use the some node base
    dataHdiv.dataOnEntities[MBVERTEX][0].getNSharedPtr(NOBASE) = dataH1.dataOnEntities[MBVERTEX][0].getNSharedPtr(NOBASE);
    dataHcurl.dataOnEntities[MBVERTEX][0].getNSharedPtr(NOBASE) = dataH1.dataOnEntities[MBVERTEX][0].getNSharedPtr(NOBASE);
    dataL2.dataOnEntities[MBVERTEX][0].getNSharedPtr(NOBASE) = dataH1.dataOnEntities[MBVERTEX][0].getNSharedPtr(NOBASE);
    std::vector<FieldApproximationBase> shape_functions_for_bases;
    for(int b = AINSWORTH_COLE_BASE;b!=LASTBASE;b++) {
      if(dataH1.bAse.test(b)) {
        switch (ApproximationBaseArray[b]) {
          case AINSWORTH_COLE_BASE:
          case LOBATTO_BASE:
          if(dataH1.spacesOnEntities[MBVERTEX].test(H1)) {
            ierr = TetPolynomialBase().getValue(
              gaussPts,
              boost::shared_ptr<BaseFunctionCtx>(
                new EntPolynomialBaseCtx(dataH1,H1,ApproximationBaseArray[b],NOBASE)
              )
            ); CHKERRQ(ierr);
          }
          if(dataH1.spacesOnEntities[MBTRI].test(HDIV)) {
            ierr = TetPolynomialBase().getValue(
              gaussPts,
              boost::shared_ptr<BaseFunctionCtx>(
                new EntPolynomialBaseCtx(dataHdiv,HDIV,ApproximationBaseArray[b],NOBASE)
              )
            ); CHKERRQ(ierr);
          }
          if(dataH1.spacesOnEntities[MBEDGE].test(HCURL)) {
            ierr = TetPolynomialBase().getValue(
              gaussPts,
              boost::shared_ptr<BaseFunctionCtx>(
                new EntPolynomialBaseCtx(dataHcurl,HCURL,ApproximationBaseArray[b],NOBASE)
              )
            ); CHKERRQ(ierr);
          }
          if(dataH1.spacesOnEntities[MBTET].test(L2)) {
            ierr = TetPolynomialBase().getValue(
              gaussPts,
              boost::shared_ptr<BaseFunctionCtx>(
                new EntPolynomialBaseCtx(dataL2,L2,ApproximationBaseArray[b],NOBASE)
              )
            ); CHKERRQ(ierr);
          }
          break;
          default:
          SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"Not yet implemented");
        }
      }
    }
  } catch (std::exception& ex) {
    std::ostringstream ss;
    ss << "thorw in method: " << ex.what() << " at line " << __LINE__ << " in file " << __FILE__;
    SETERRQ(PETSC_COMM_SELF,MOFEM_STD_EXCEPTION_THROW,ss.str().c_str());
  }
  PetscFunctionReturn(0);
}

PetscErrorCode VolumeElementForcesAndSourcesCore::operator()() {
  PetscFunctionBegin;

  try {

    if(numeredEntFiniteElementPtr->getEntType() != MBTET) PetscFunctionReturn(0);

    ierr = calculateVolumeAndJacobian(); CHKERRQ(ierr);
    ierr = getSpaceBaseAndOrderOnElement(); CHKERRQ(ierr);
    ierr = setIntegartionPts(); CHKERRQ(ierr);
    if(nbGaussPts == 0) PetscFunctionReturn(0);
    ierr = calculateCoordinatesAtGaussPts(); CHKERRQ(ierr);
    ierr = calculateBaseFunctionsOnElement(); CHKERRQ(ierr);

    try {
      ierr = opSetInvJacH1.opRhs(dataH1); CHKERRQ(ierr);
      if(dataH1.spacesOnEntities[MBEDGE].test(HCURL)) {
        ierr = opCovariantPiolaTransform.opRhs(dataHcurl); CHKERRQ(ierr);
        ierr = opSetInvJacHdivAndHcurl.opRhs(dataHcurl); CHKERRQ(ierr);
      }
      if(dataH1.spacesOnEntities[MBTRI].test(HDIV)) {
        ierr = opContravariantPiolaTransform.opRhs(dataHdiv); CHKERRQ(ierr);
        ierr = opSetInvJacHdivAndHcurl.opRhs(dataHdiv); CHKERRQ(ierr);
      }
      if(dataH1.spacesOnEntities[MBTET].test(L2)) {
        ierr = opSetInvJacH1.opRhs(dataL2); CHKERRQ(ierr);
      }
    } catch (std::exception& ex) {
      std::ostringstream ss;
      ss << "thorw in method: " << ex.what() << " at line " << __LINE__ << " in file " << __FILE__;
      SETERRQ(PETSC_COMM_SELF,MOFEM_STD_EXCEPTION_THROW,ss.str().c_str());
    }

    if(
      dataPtr->get<FieldName_mi_tag>().find(meshPositionsFieldName)!=
      dataPtr->get<FieldName_mi_tag>().end()
    ) {
      const Field* field_struture = mField.get_field_structure(meshPositionsFieldName);
      BitFieldId id = field_struture->getId();

      if((numeredEntFiniteElementPtr->getBitFieldIdData()&id).none()) {
        SETERRQ(PETSC_COMM_SELF,MOFEM_NOT_FOUND,"no MESH_NODE_POSITIONS in element data");
      }

      ierr = getEdgesDataOrderSpaceAndBase(dataH1,meshPositionsFieldName); CHKERRQ(ierr);
      ierr = getTrisDataOrderSpaceAndBase(dataH1,meshPositionsFieldName); CHKERRQ(ierr);
      ierr = getTetDataOrderSpaceAndBase(dataH1,meshPositionsFieldName); CHKERRQ(ierr);
      ierr = getNodesFieldData(dataH1,meshPositionsFieldName); CHKERRQ(ierr);
      if(dataH1.dataOnEntities[MBVERTEX][0].getFieldData().size()!=12) {
        SETERRQ(PETSC_COMM_SELF,MOFEM_NOT_FOUND,"no MESH_NODE_POSITIONS in element data");
      }
      ierr = getEdgesFieldData(dataH1,meshPositionsFieldName); CHKERRQ(ierr);
      ierr = getTrisFieldData(dataH1,meshPositionsFieldName); CHKERRQ(ierr);
      ierr = getTetsFieldData(dataH1,meshPositionsFieldName); CHKERRQ(ierr);
      try {
        ierr = opHOatGaussPoints.opRhs(dataH1); CHKERRQ(ierr);
        hoGaussPtsInvJac.resize(hoGaussPtsJac.size1(),hoGaussPtsJac.size2(),false);
        FTensor::Tensor2<double*,3,3> jac(
          &hoGaussPtsJac(0,0),&hoGaussPtsJac(0,1),&hoGaussPtsJac(0,2),
          &hoGaussPtsJac(0,3),&hoGaussPtsJac(0,4),&hoGaussPtsJac(0,5),
          &hoGaussPtsJac(0,6),&hoGaussPtsJac(0,7),&hoGaussPtsJac(0,8),9
        );
        FTensor::Tensor2<double*,3,3> inv_jac(
          &hoGaussPtsInvJac(0,0),&hoGaussPtsInvJac(0,1),&hoGaussPtsInvJac(0,2),
          &hoGaussPtsInvJac(0,3),&hoGaussPtsInvJac(0,4),&hoGaussPtsInvJac(0,5),
          &hoGaussPtsInvJac(0,6),&hoGaussPtsInvJac(0,7),&hoGaussPtsInvJac(0,8),9
        );
        hoGaussPtsDetJac.resize(nbGaussPts,false);
        FTensor::Tensor0<double*> det(&hoGaussPtsDetJac[0]);
        for(int gg = 0;gg!=nbGaussPts;gg++) {
          ierr = determinantTensor3by3(jac,det); CHKERRQ(ierr);
          // if(det<0) {
          //   SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"Negative volume");
          // }
          ierr = invertTensor3by3(jac,det,inv_jac); CHKERRQ(ierr);
          ++jac;
          ++inv_jac;
          ++det;
        }
        // MatrixDouble jac(3,3);
        // for(int gg = 0;gg<nbGaussPts;gg++) {
          // cblas_dcopy(9,&hoGaussPtsJac(gg,0),1,&jac(0,0),1);
        //   hoGaussPtsDetJac[gg] = ShapeDetJacVolume(&jac(0,0));
        //   if(hoGaussPtsDetJac[gg]<0) {
        //     SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"Negative volume");
        //   }
        //   ierr = ShapeInvJacVolume(&hoGaussPtsInvJac(gg,0)); CHKERRQ(ierr);
        // }
        ierr = opSetHoInvJacH1.opRhs(dataH1); CHKERRQ(ierr);
        if(dataH1.spacesOnEntities[MBTET].test(L2)) {
          ierr = opSetHoInvJacH1.opRhs(dataL2); CHKERRQ(ierr);
        }
        if(dataH1.spacesOnEntities[MBTRI].test(HDIV)) {
          ierr = opHoContravariantTransform.opRhs(dataHdiv); CHKERRQ(ierr);
          ierr = opSetHoInvJacHdivAndHcurl.opRhs(dataHdiv); CHKERRQ(ierr);
        }
        if(dataH1.spacesOnEntities[MBEDGE].test(HCURL)) {
          ierr = opHoCovariantTransform.opRhs(dataHcurl); CHKERRQ(ierr);
          ierr = opSetHoInvJacHdivAndHcurl.opRhs(dataHcurl); CHKERRQ(ierr);
        }

      } catch (std::exception& ex) {
        std::ostringstream ss;
        ss << "problem with indices in method: " << ex.what() << " at line " << __LINE__ << " in file " << __FILE__;
        SETERRQ(PETSC_COMM_SELF,MOFEM_STD_EXCEPTION_THROW,ss.str().c_str());
      }
    } else {
      hoCoordsAtGaussPts.resize(0,0,false);
      hoGaussPtsInvJac.resize(0,0,false);
      hoGaussPtsDetJac.resize(0,false);
      try {
        for(int b = AINSWORTH_COLE_BASE;b!=LASTBASE;b++) {
          if(dataH1.dataOnEntities[MBVERTEX][0].getDiffN(ApproximationBaseArray[b]).size1()!=4) continue;
          if(dataH1.dataOnEntities[MBVERTEX][0].getDiffN(ApproximationBaseArray[b]).size2()!=3) continue;
          MatrixDouble diffN(nbGaussPts,12);
          for(int gg = 0;gg<nbGaussPts;gg++) {
            for(int nn = 0;nn<4;nn++) {
              for(int dd = 0;dd<3;dd++) {
                diffN(gg,nn*3+dd) = dataH1.dataOnEntities[MBVERTEX][0].getDiffN(ApproximationBaseArray[b])(nn,dd);
              }
            }
          }
          dataH1.dataOnEntities[MBVERTEX][0].getDiffN(ApproximationBaseArray[b]).resize(diffN.size1(),diffN.size2(),false);
          dataH1.dataOnEntities[MBVERTEX][0].getDiffN(ApproximationBaseArray[b]).data().swap(diffN.data());
        }
      } catch (std::exception& ex) {
        std::ostringstream ss;
        ss << "thorw in method: " << ex.what() << " at line " << __LINE__ << " in file " << __FILE__;
        SETERRQ(PETSC_COMM_SELF,MOFEM_STD_EXCEPTION_THROW,ss.str().c_str());
      }
    }

    const UserDataOperator::OpType types[2] = {
      UserDataOperator::OPROW, UserDataOperator::OPCOL
    };
    std::vector<std::string> last_eval_field_name(2);
    DataForcesAndSurcesCore *op_data[2];
    FieldSpace space[2];
    FieldApproximationBase base[2];

    boost::ptr_vector<UserDataOperator>::iterator oit,hi_oit;
    oit = opPtrVector.begin();
    hi_oit = opPtrVector.end();

    for(;oit!=hi_oit;oit++) {

      oit->setPtrFE(this);

      for(int ss = 0;ss!=2;ss++) {

        std::string field_name = !ss ? oit->rowFieldName : oit->colFieldName;
        const Field* field_struture = mField.get_field_structure(field_name);
        BitFieldId data_id = field_struture->getId();

        if((oit->getNumeredEntFiniteElementPtr()->getBitFieldIdData()&data_id).none()) {
          SETERRQ2(
            PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"no data field < %s > on finite element < %s >",
            field_name.c_str(),feName.c_str()
          );
        }

        if(oit->getOpType()&types[ss] || oit->getOpType()&UserDataOperator::OPROWCOL) {

          space[ss] = field_struture->getSpace();
          switch(space[ss]) {
            case NOSPACE:
            SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"unknown space");
            break;
            case H1:
            op_data[ss] = !ss ? &dataH1 : &derivedDataH1;
            break;
            case HCURL:
            op_data[ss] = !ss ? &dataHcurl : &derivedDataHcurl;
            break;
            case HDIV:
            op_data[ss] = !ss ? &dataHdiv : &derivedDataHdiv;
            break;
            case L2:
            op_data[ss] = !ss ? &dataL2 : &derivedDataL2;
            break;
            case NOFIELD:
            op_data[ss] = !ss ? &dataNoField : &dataNoFieldCol;
            break;
            case LASTSPACE:
            SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"unknown space");
            break;
          }

          base[ss] = field_struture->getApproxBase();
          switch(base[ss]) {
            case AINSWORTH_COLE_BASE:
            break;
            case LOBATTO_BASE:
            break;
            default:
            SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"unknown or not implemented base");
            break;
          }

          if(last_eval_field_name[ss]!=field_name) {

            switch(space[ss]) {
              case NOSPACE:
              SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"unknown space");
              break;
              case H1:
              if(!ss) {
                ierr = getRowNodesIndices(*op_data[ss],field_name); CHKERRQ(ierr);
              } else {
                ierr = getColNodesIndices(*op_data[ss],field_name); CHKERRQ(ierr);
              }
              ierr = getNodesFieldData(*op_data[ss],field_name); CHKERRQ(ierr);
              case HCURL:
              if(!ss) {
                ierr = getEdgesRowIndices(*op_data[ss],field_name); CHKERRQ(ierr);
              } else {
                ierr = getEdgesColIndices(*op_data[ss],field_name); CHKERRQ(ierr);
              }
              ierr = getEdgesDataOrderSpaceAndBase(*op_data[ss],field_name); CHKERRQ(ierr);
              ierr = getEdgesFieldData(*op_data[ss],field_name); CHKERRQ(ierr);
              case HDIV:
              if(!ss) {
                ierr = getTrisRowIndices(*op_data[ss],field_name); CHKERRQ(ierr);
              } else {
                ierr = getTrisColIndices(*op_data[ss],field_name); CHKERRQ(ierr);
              }
              ierr = getTrisDataOrderSpaceAndBase(*op_data[ss],field_name); CHKERRQ(ierr);
              ierr = getTrisFieldData(*op_data[ss],field_name); CHKERRQ(ierr);
              case L2:
              if(!ss) {
                ierr = getTetsRowIndices(*op_data[ss],field_name); CHKERRQ(ierr);
              } else {
                ierr = getTetsColIndices(*op_data[ss],field_name); CHKERRQ(ierr);
              }
              ierr = getTetDataOrderSpaceAndBase(*op_data[ss],field_name); CHKERRQ(ierr);
              ierr = getTetsFieldData(*op_data[ss],field_name); CHKERRQ(ierr);
              break;
              case NOFIELD:
              if(!getNinTheLoop()) {
                // NOFIELD data are the same for each element, can be retrieved only once
                if(!ss) {
                  ierr = getNoFieldRowIndices(*op_data[ss],field_name); CHKERRQ(ierr);
                } else {
                  ierr = getNoFieldColIndices(*op_data[ss],field_name); CHKERRQ(ierr);
                }
                ierr = getNoFieldFieldData(*op_data[ss],field_name); CHKERRQ(ierr);
              }
              break;
              case LASTSPACE:
              SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"unknown space");
              break;
            }
            last_eval_field_name[ss]=field_name;

          }
        }
      }

      if(oit->getOpType()&UserDataOperator::OPROW) {
        try {
          ierr = oit->opRhs(
            *op_data[0],
            oit->doVerticesRow,
            oit->doEdgesRow,
            oit->doQuadsRow,
            oit->doTrisRow,
            oit->doTetsRow,
            false
          ); CHKERRQ(ierr);
        } catch (std::exception& ex) {
          std::ostringstream ss;
          ss << "thorw in method: " << ex.what() << " at line " << __LINE__ << " in file " << __FILE__;
          ss << " operator number " << std::distance<boost::ptr_vector<UserDataOperator>::iterator>(opPtrVector.begin(),oit);
          SETERRQ(PETSC_COMM_SELF,MOFEM_STD_EXCEPTION_THROW,ss.str().c_str());
        }
      }

      if(oit->getOpType()&UserDataOperator::OPCOL) {
        try {
          ierr = oit->opRhs(
            *op_data[1],
            oit->doVerticesCol,
            oit->doEdgesCol,
            oit->doQuadsCol,
            oit->doTrisCol,
            oit->doTetsCol,
            false
          ); CHKERRQ(ierr);
        } catch (std::exception& ex) {
          std::ostringstream ss;
          ss << "thorw in method: " << ex.what() << " at line " << __LINE__ << " in file " << __FILE__;
          ss << " operator number " << std::distance<boost::ptr_vector<UserDataOperator>::iterator>(opPtrVector.begin(),oit);
          SETERRQ(PETSC_COMM_SELF,MOFEM_STD_EXCEPTION_THROW,ss.str().c_str());
        }
      }

      if(oit->getOpType()&UserDataOperator::OPROWCOL) {
        try {
          ierr = oit->opLhs(*op_data[0],*op_data[1],oit->sYmm); CHKERRQ(ierr);
        } catch (std::exception& ex) {
          std::ostringstream ss;
          ss << "thorw in method: " << ex.what() << " at line " << __LINE__ << " in file " << __FILE__;
          ss << " operator number " << std::distance<boost::ptr_vector<UserDataOperator>::iterator>(opPtrVector.begin(),oit);
          SETERRQ(PETSC_COMM_SELF,MOFEM_STD_EXCEPTION_THROW,ss.str().c_str());
        }
      }

    }

  } catch (std::exception& ex) {
    std::ostringstream ss;
    ss << "thorw in method: " << ex.what() << " at line " << __LINE__ << " in file " << __FILE__;
    SETERRQ(PETSC_COMM_SELF,MOFEM_STD_EXCEPTION_THROW,ss.str().c_str());
  }

  PetscFunctionReturn(0);
}

PetscErrorCode VolumeElementForcesAndSourcesCore::UserDataOperator::getDivergenceOfHDivBaseFunctions(
  int side,EntityType type,
  DataForcesAndSurcesCore::EntData &data,
  int gg,
  VectorDouble &div
) {
  PetscFunctionBegin;

  try {

    int nb_dofs = data.getFieldData().size();
    if(nb_dofs==0) PetscFunctionReturn(0);

    if(data.getSpace()!=HDIV) {
      SETERRQ1(
        PETSC_COMM_SELF,
        MOFEM_DATA_INCONSISTENCY,
        "This function should be used for HDIV used but is used with %s",
        FieldSpaceNames[data.getSpace()]
      );
    }

    if((unsigned int)nb_dofs != data.getDiffHdivN().size2()/9) {
      std::cerr << "side " << side << " type " << type << std::endl;
      SETERRQ3(
        PETSC_COMM_SELF,
        MOFEM_DATA_INCONSISTENCY,
        "Data inositency, wrong number of dofs  = %s "
        "%d != %d/9",
        FieldSpaceNames[data.getSpace()],
        nb_dofs,data.getDiffHdivN().size2()
      );
    }

    div.resize(nb_dofs,false);

    FTensor::Tensor0<double*> t_div(&*div.data().begin());
    const double *grad_ptr = &data.getDiffHdivN()(gg,0);
    FTensor::Tensor1<const double*,3> t_grad_base(
      grad_ptr,
      &grad_ptr[HDIV1_1],
      &grad_ptr[HDIV2_2],9
    );
    FTensor::Index<'i',3> i;
    for(int dd = 0;dd<nb_dofs;dd++) {
      t_div = t_grad_base(0)+t_grad_base(1)+t_grad_base(2);
      ++t_div;
      ++t_grad_base;
    }

  } catch (std::exception& ex) {
    std::ostringstream ss;
    ss << "thorw in method: " << ex.what() << " at line " << __LINE__ << " in file " << __FILE__;
    SETERRQ(PETSC_COMM_SELF,1,ss.str().c_str());
  }

  PetscFunctionReturn(0);
}

PetscErrorCode VolumeElementForcesAndSourcesCore::UserDataOperator::getCurlOfHCurlBaseFunctions(
  int side,
  EntityType type,
  DataForcesAndSurcesCore::EntData &data,
  int gg,
  MatrixDouble &curl
) {
  PetscFunctionBegin;

  try {

    int nb_dofs = data.getFieldData().size();
    if(nb_dofs==0) PetscFunctionReturn(0);

    if(data.getSpace()!=HCURL) {
      SETERRQ1(
        PETSC_COMM_SELF,
        MOFEM_DATA_INCONSISTENCY,
        "This function should be used for HCURL used but is used with %s",
        FieldSpaceNames[data.getSpace()]
      );
    }

    if((unsigned int)nb_dofs != data.getDiffHcurlN().size2()/9) {
      std::cerr << "side " << side << " type " << type << std::endl;
      SETERRQ3(
        PETSC_COMM_SELF,
        MOFEM_DATA_INCONSISTENCY,
        "Data inositency, wrong number of dofs  = %s "
        "%d != %d/9",
        FieldSpaceNames[data.getSpace()],
        nb_dofs,data.getDiffHcurlN().size2()
      );
    }

    curl.resize(nb_dofs,3,false);
    FTensor::Tensor1<double*,3> t_curl(&curl(0,0),&curl(0,1),&curl(0,2),3);
    const double *grad_ptr = &data.getDiffHcurlN()(gg,0);

    FTensor::Tensor2<const double*,3,3> t_grad_base(
      grad_ptr,&grad_ptr[HCURL0_1],&grad_ptr[HCURL0_2],
      &grad_ptr[HCURL1_0],&grad_ptr[HCURL1_1],&grad_ptr[HCURL1_2],
      &grad_ptr[HCURL2_0],&grad_ptr[HCURL2_1],&grad_ptr[HCURL2_2],9
    );
    FTensor::Index<'i',3> i;
    for(int dd = 0;dd<nb_dofs;dd++) {
      t_curl(0) = t_grad_base(2,1)-t_grad_base(1,2);
      t_curl(1) = t_grad_base(0,2)-t_grad_base(2,0);
      t_curl(2) = t_grad_base(1,0)-t_grad_base(0,1);
      ++t_curl;
      ++t_grad_base;
    }

  } catch (std::exception& ex) {
    std::ostringstream ss;
    ss << "thorw in method: " << ex.what() << " at line " << __LINE__ << " in file " << __FILE__;
    SETERRQ(PETSC_COMM_SELF,1,ss.str().c_str());
  }

  PetscFunctionReturn(0);
}

PetscErrorCode VolumeElementForcesAndSourcesCoreOnSide::setGaussPts(int order) {
  PetscFunctionBegin;
  if(faceFEPtr==NULL) {
    SETERRQ(PETSC_COMM_WORLD,MOFEM_DATA_INCONSISTENCY,"Pointer to face element is not set");
  }
  const EntityHandle face_entity = faceFEPtr->numeredEntFiniteElementPtr->getEnt();
  SideNumber_multiIndex& side_table = const_cast<SideNumber_multiIndex&>(
    numeredEntFiniteElementPtr->getSideNumberTable()
  );
  SideNumber_multiIndex::nth_index<0>::type::iterator sit = side_table.get<0>().find(face_entity);
  if(sit==side_table.get<0>().end()) {
    SETERRQ(PETSC_COMM_WORLD,MOFEM_DATA_INCONSISTENCY,"Face can not be found on volume element");
  }
  faceSense = (*sit)->sense;
  faceSideNumber = (*sit)->side_number;
  fill(tetConnMap,&tetConnMap[4],-1);
  for(int nn = 0;nn!=3;nn++) {
    faceConnMap[nn] = std::distance(conn,find(conn,&conn[4],faceFEPtr->conn[nn]));
    tetConnMap[faceConnMap[nn]] = nn;
    if(faceConnMap[nn]>3) {
      SETERRQ(PETSC_COMM_WORLD,MOFEM_DATA_INCONSISTENCY,"No common node on face and element can not be found");
    }
  }
  oppositeNode = std::distance(tetConnMap,find(tetConnMap,&tetConnMap[4],-1));
  const int nb_gauss_pts = faceFEPtr->gaussPts.size2();
  gaussPts.resize(4,nb_gauss_pts,false);
  gaussPts.clear();
  const MatrixDouble &face_shape_funtions = faceFEPtr->dataH1.dataOnEntities[MBVERTEX][0].getN(NOBASE);
  const double tet_coords[] = { 0,0,0, 1,0,0, 0,1,0, 0,0,1 };
  for(int gg = 0;gg!=nb_gauss_pts;gg++) {
    gaussPts(0,gg) =
    face_shape_funtions(gg,0)*tet_coords[3*faceConnMap[0]+0]+
    face_shape_funtions(gg,1)*tet_coords[3*faceConnMap[1]+0]+
    face_shape_funtions(gg,2)*tet_coords[3*faceConnMap[2]+0];
    gaussPts(1,gg) =
    face_shape_funtions(gg,0)*tet_coords[3*faceConnMap[0]+1]+
    face_shape_funtions(gg,1)*tet_coords[3*faceConnMap[1]+1]+
    face_shape_funtions(gg,2)*tet_coords[3*faceConnMap[2]+1];
    gaussPts(2,gg) =
    face_shape_funtions(gg,0)*tet_coords[3*faceConnMap[0]+2]+
    face_shape_funtions(gg,1)*tet_coords[3*faceConnMap[1]+2]+
    face_shape_funtions(gg,2)*tet_coords[3*faceConnMap[2]+2];
    gaussPts(3,gg) = faceFEPtr->gaussPts(2,gg);
  }
  PetscFunctionReturn(0);
}

VectorDouble& VolumeElementForcesAndSourcesCoreOnSide::UserDataOperator::getNormal() {
  return getFaceFEPtr()->normal;
}

MatrixDouble& VolumeElementForcesAndSourcesCoreOnSide::UserDataOperator::getNormalsAtGaussPt() {
  return getFaceFEPtr()->nOrmals_at_GaussPt;
}

MatrixDouble& VolumeElementForcesAndSourcesCoreOnSide::UserDataOperator::getFaceCoordsAtGaussPts() {
  return getFaceFEPtr()->coordsAtGaussPts;
}

/** \brief if higher order geometry return normals at Gauss pts.
  *
  * \param gg gauss point number
  */
ublas::matrix_row<MatrixDouble > VolumeElementForcesAndSourcesCoreOnSide::UserDataOperator::getNormalsAtGaussPt(const int gg) {
  return ublas::matrix_row<MatrixDouble >(getNormalsAtGaussPt(),gg);
}


}
