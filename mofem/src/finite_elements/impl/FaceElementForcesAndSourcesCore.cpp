/** \file FaceElementForcesAndSourcesCore.cpp

\brief Implementation of face element

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
#include <FEMultiIndices.hpp>
#include <ProblemsMultiIndices.hpp>
#include <AdjacencyMultiIndices.hpp>
#include <BCMultiIndices.hpp>
#include <CoreDataStructures.hpp>
#include <SeriesMultiIndices.hpp>

#include <LoopMethods.hpp>
#include <Interface.hpp>
#include <MeshRefinement.hpp>
#include <PrismInterface.hpp>
#include <SeriesRecorder.hpp>
#include <Core.hpp>

#include <FTensor.hpp>
#include <DataStructures.hpp>
#include <DataOperators.hpp>
#include <ForcesAndSurcesCore.hpp>
#include <VolumeElementForcesAndSourcesCore.hpp>
#include <FaceElementForcesAndSourcesCore.hpp>

#include <BaseFunction.hpp>
#include <EntPolynomialBaseCtx.hpp>
#include <TriPolynomialBase.hpp>

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

PetscErrorCode FaceElementForcesAndSourcesCore::UserDataOperator::loopSideVolumes(
  const string &fe_name,VolumeElementForcesAndSourcesCoreOnSide &method
) {
  
  PetscFunctionBegin;

  const EntityHandle ent = getNumeredEntFiniteElementPtr()->getEnt();
  const Problem *problem_ptr = getFEMethod()->problemPtr;
  Range adjacent_volumes;
  ierr = getTriFE()->mField.get_adjacencies_any(ent,3,adjacent_volumes); CHKERRQ(ierr);
  typedef NumeredEntFiniteElement_multiIndex::index<Composite_Name_And_Ent_mi_tag>::type FEByComposite;
  FEByComposite &numered_fe =
  (const_cast<NumeredEntFiniteElement_multiIndex&>(
    problem_ptr->numeredFiniteElements)
  ).get<Composite_Name_And_Ent_mi_tag>();

  method.feName = fe_name;

  ierr = method.setFaceFEPtr(getTriFE()); CHKERRQ(ierr);
  ierr = method.copyBasicMethod(*getFEMethod()); CHKERRQ(ierr);
  ierr = method.copyKsp(*getFEMethod()); CHKERRQ(ierr);
  ierr = method.copySnes(*getFEMethod()); CHKERRQ(ierr);
  ierr = method.copyTs(*getFEMethod()); CHKERRQ(ierr);

  try {
    ierr = method.preProcess(); CHKERRQ(ierr);
  } catch (const std::exception& ex) {
    std::ostringstream ss;
    ss << "throw in method: " << ex.what() << " at line " << __LINE__ << " in file " << __FILE__ << std::endl;
    SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,ss.str().c_str());
  }

  int nn = 0;
  method.loopSize = adjacent_volumes.size();
  for(Range::iterator vit = adjacent_volumes.begin();vit!=adjacent_volumes.end();vit++) {
    FEByComposite::iterator miit = numered_fe.find(boost::make_tuple(fe_name,*vit));
    if(miit!=numered_fe.end()) {
      // cerr << **miit << endl;
      // cerr << &(**miit) << endl;
      // cerr << (*miit)->getEnt() << endl;
      method.nInTheLoop = nn++;
      method.numeredEntFiniteElementPtr = *miit;
      method.dataPtr = (*miit)->sPtr->data_dofs;
      method.rowPtr = (*miit)->rows_dofs;
      method.colPtr = (*miit)->cols_dofs;

      try {
        ierr = method(); CHKERRQ(ierr);
      } catch (const std::exception& ex) {
        std::ostringstream ss;
        ss << "throw in method: " << ex.what() << " at line " << __LINE__ << " in file " << __FILE__ << std::endl;
        SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,ss.str().c_str());
      }
    }
  }

  try {
    ierr = method.postProcess(); CHKERRQ(ierr);
  } catch (const std::exception& ex) {
    std::ostringstream ss;
    ss << "throw in method: " << ex.what() << " at line " << __LINE__ << " in file " << __FILE__ << std::endl;
    SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,ss.str().c_str());
  }

  PetscFunctionReturn(0);
}

PetscErrorCode FaceElementForcesAndSourcesCore::calculateAreaAndNormal() {
  PetscFunctionBegin;
  EntityHandle ent = numeredEntFiniteElementPtr->getEnt();
  rval = mField.get_moab().get_connectivity(ent,conn,num_nodes,true); CHKERRQ_MOAB(rval);
  coords.resize(num_nodes*3,false);
  rval = mField.get_moab().get_coords(conn,num_nodes,&*coords.data().begin()); CHKERRQ_MOAB(rval);
  double diff_n[6];
  ierr = ShapeDiffMBTRI(diff_n); CHKERRQ(ierr);
  nOrmal.resize(3,false);
  ierr = ShapeFaceNormalMBTRI(
    diff_n,&*coords.data().begin(),&*nOrmal.data().begin()
  ); CHKERRQ(ierr);
  aRea = cblas_dnrm2(3,&*nOrmal.data().begin(),1)*0.5;
  tangentOne.resize(3,false);
  tangentTwo.resize(3,false);
  for(int dd = 0;dd!=3;dd++) {
    tangentOne[dd] = cblas_ddot(3,&diff_n[0],2,&coords[dd],3);
    tangentTwo[dd] = cblas_ddot(3,&diff_n[1],2,&coords[dd],3);
  }
  PetscFunctionReturn(0);
}

PetscErrorCode FaceElementForcesAndSourcesCore::setIntegartionPts() {
  PetscFunctionBegin;
  // Set integration points
  int order_data = getMaxDataOrder();
  int order_row = getMaxRowOrder();
  int order_col = getMaxColOrder();
  int rule = getRule(order_row,order_col,order_data);
  if(rule >= 0) {
    if(rule<QUAD_2D_TABLE_SIZE) {
      if(QUAD_2D_TABLE[rule]->dim!=2) {
        SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"wrong dimension");
      }
      if(QUAD_2D_TABLE[rule]->order<rule) {
        SETERRQ2(
          PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"wrong order %d != %d",
          QUAD_2D_TABLE[rule]->order,rule
        );
      }
      nbGaussPts = QUAD_2D_TABLE[rule]->npoints;
      gaussPts.resize(3,nbGaussPts,false);
      cblas_dcopy(
        nbGaussPts,&QUAD_2D_TABLE[rule]->points[1],3,&gaussPts(0,0),1
      );
      cblas_dcopy(
        nbGaussPts,&QUAD_2D_TABLE[rule]->points[2],3,&gaussPts(1,0),1
      );
      cblas_dcopy(
        nbGaussPts,QUAD_2D_TABLE[rule]->weights,1,&gaussPts(2,0),1
      );
      dataH1.dataOnEntities[MBVERTEX][0].getN(NOBASE).resize(nbGaussPts,3,false);
      double *shape_ptr = &*dataH1.dataOnEntities[MBVERTEX][0].getN(NOBASE).data().begin();
      cblas_dcopy(
        3*nbGaussPts,QUAD_2D_TABLE[rule]->points,1,shape_ptr,1
      );
    } else {
      SETERRQ2(
        PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"rule > quadrature order %d < %d",
        rule,QUAD_2D_TABLE_SIZE
      );
      nbGaussPts = 0;
    }
  } else {
    // If rule is negative, set user defined integration points
    ierr = setGaussPts(order_row,order_col,order_data); CHKERRQ(ierr);
    nbGaussPts = gaussPts.size2();
    dataH1.dataOnEntities[MBVERTEX][0].getN(NOBASE).resize(nbGaussPts,3,false);
    if(nbGaussPts) {
      ierr = ShapeMBTRI(
        &*dataH1.dataOnEntities[MBVERTEX][0].getN(NOBASE).data().begin(),
        &gaussPts(0,0),
        &gaussPts(1,0),
        nbGaussPts
      ); CHKERRQ(ierr);
    }
  }
  PetscFunctionReturn(0);
}

PetscErrorCode FaceElementForcesAndSourcesCore::getSpaceBaseAndOrderOnElement() {
  PetscFunctionBegin;
  // Get spaces order/base and sense of entities.

  ierr = getSpacesAndBaseOnEntities(dataH1); CHKERRQ(ierr);

  //H1
  if(dataH1.spacesOnEntities[MBEDGE].test(H1)) {
    ierr = getEdgesSense(dataH1); CHKERRQ(ierr);
    ierr = getEdgesDataOrder(dataH1,H1); CHKERRQ(ierr);
  }
  if(dataH1.spacesOnEntities[MBTRI].test(H1)) {
    ierr = getTrisSense(dataH1); CHKERRQ(ierr);
    ierr = getTrisDataOrder(dataH1,H1); CHKERRQ(ierr);
  }

  //Hcurl
  if(dataH1.spacesOnEntities[MBEDGE].test(HCURL)) {
    ierr = getEdgesSense(dataHcurl); CHKERRQ(ierr);
    ierr = getEdgesDataOrder(dataHcurl,HCURL); CHKERRQ(ierr);
    dataHcurl.spacesOnEntities[MBEDGE].set(HCURL);
  }
  if(dataH1.spacesOnEntities[MBTRI].test(HCURL)) {
    ierr = getTrisSense(dataHcurl); CHKERRQ(ierr);
    ierr = getTrisDataOrder(dataHcurl,HCURL); CHKERRQ(ierr);
    dataHcurl.spacesOnEntities[MBTRI].set(HCURL);
  }

  //Hdiv
  if(dataH1.spacesOnEntities[MBTRI].test(HDIV)) {
    ierr = getTrisSense(dataHdiv); CHKERRQ(ierr);
    ierr = getTrisDataOrder(dataHdiv,HDIV); CHKERRQ(ierr);
    dataHcurl.spacesOnEntities[MBTRI].set(HDIV);
  }

  //L2
  if(dataH1.spacesOnEntities[MBTRI].test(L2)) {
    ierr = getTrisSense(dataL2); CHKERRQ(ierr);
    ierr = getTrisDataOrder(dataL2,L2); CHKERRQ(ierr);
    dataHcurl.spacesOnEntities[MBTRI].set(L2);
  }

  PetscFunctionReturn(0);
}

PetscErrorCode FaceElementForcesAndSourcesCore::calculateCoordinatesAtGaussPts() {
  PetscFunctionBegin;
  double *shape_functions = &*dataH1.dataOnEntities[MBVERTEX][0].getN(NOBASE).data().begin();
  coordsAtGaussPts.resize(nbGaussPts,3,false);
  for(int gg = 0;gg<nbGaussPts;gg++) {
    for(int dd = 0;dd<3;dd++) {
      coordsAtGaussPts(gg,dd) = cblas_ddot(3,&shape_functions[3*gg],1,&coords[dd],3);
    }
  }
  PetscFunctionReturn(0);
}

PetscErrorCode FaceElementForcesAndSourcesCore::calculateBaseFunctionsOnElement() {
  PetscFunctionBegin;
  // Calculate base base functions for faces.
  try {

    for(int b = AINSWORTH_LEGENDRE_BASE;b!=LASTBASE;b++) {
      if(dataH1.bAse.test(b)) {
        switch (ApproximationBaseArray[b]) {
          case NOBASE:
          break;
          case AINSWORTH_LEGENDRE_BASE:
          case AINSWORTH_LOBBATO_BASE:
          if(
            dataH1.spacesOnEntities[MBVERTEX].test(H1)&&
            dataH1.basesOnEntities[MBVERTEX].test(b)
          ) {
            ierr = TriPolynomialBase().getValue(
              gaussPts,
              boost::shared_ptr<BaseFunctionCtx>(
                new EntPolynomialBaseCtx(dataH1,H1,ApproximationBaseArray[b],NOBASE)
              )
            ); CHKERRQ(ierr);
          }
          if(
            dataH1.spacesOnEntities[MBEDGE].test(HCURL)&&
            dataH1.basesOnEntities[MBEDGE].test(b)
          ) {
            ierr = TriPolynomialBase().getValue(
              gaussPts,
              boost::shared_ptr<BaseFunctionCtx>(
                new EntPolynomialBaseCtx(dataHcurl,HCURL,ApproximationBaseArray[b],NOBASE)
              )
            ); CHKERRQ(ierr);
          }
          if(
            dataH1.spacesOnEntities[MBTRI].test(HDIV)&&
            dataH1.basesOnEntities[MBTRI].test(b)
          ) {
            ierr = TriPolynomialBase().getValue(
              gaussPts,
              boost::shared_ptr<BaseFunctionCtx>(
                new EntPolynomialBaseCtx(dataHdiv,HDIV,ApproximationBaseArray[b],NOBASE)
              )
            ); CHKERRQ(ierr);
          }
          if(
            dataH1.spacesOnEntities[MBTRI].test(L2)&&
            dataH1.basesOnEntities[MBTRI].test(b)
          ) {
            ierr = TriPolynomialBase().getValue(
              gaussPts,
              boost::shared_ptr<BaseFunctionCtx>(
                new EntPolynomialBaseCtx(dataL2,L2,ApproximationBaseArray[b],NOBASE)
              )
            ); CHKERRQ(ierr);
          }
          break;
          case DEMKOWICZ_JACOBI_BASE:
          if(
            dataH1.spacesOnEntities[MBTRI].test(HDIV)&&
            dataH1.basesOnEntities[MBTRI].test(b)
          ) {
            ierr = TriPolynomialBase().getValue(
              gaussPts,
              boost::shared_ptr<BaseFunctionCtx>(
                new EntPolynomialBaseCtx(dataHdiv,HDIV,ApproximationBaseArray[b],NOBASE)
              )
            ); CHKERRQ(ierr);
          }
          break;
          default:
          SETERRQ1(
            PETSC_COMM_SELF,
            MOFEM_DATA_INCONSISTENCY,
            "Base <%s> not yet implemented",
            ApproximationBaseNames[ApproximationBaseArray[b]]
          );
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

PetscErrorCode FaceElementForcesAndSourcesCore::calculateHoNormal() {
  PetscFunctionBegin;
  // Check if field for high-order geometry is set and if it is set calculate
  // higher-order normals and face tangent vectors.
  if(
    dataPtr->get<FieldName_mi_tag>().find(meshPositionsFieldName)!=
    dataPtr->get<FieldName_mi_tag>().end()
  ) {
    const Field* field_struture = mField.get_field_structure(meshPositionsFieldName);
    BitFieldId id = field_struture->getId();

    if((numeredEntFiniteElementPtr->getBitFieldIdData()&id).none()) {
      SETERRQ(PETSC_COMM_SELF,MOFEM_NOT_FOUND,"no MESH_NODE_POSITIONS in element data");
    }

    // Calculate normal for high-order geometry
    ierr = getEdgesDataOrderSpaceAndBase(dataH1,meshPositionsFieldName); CHKERRQ(ierr);
    ierr = getTrisDataOrderSpaceAndBase(dataH1,meshPositionsFieldName); CHKERRQ(ierr);
    ierr = getNodesFieldData(dataH1,meshPositionsFieldName); CHKERRQ(ierr);
    ierr = getEdgesFieldData(dataH1,meshPositionsFieldName); CHKERRQ(ierr);
    ierr = getTrisFieldData(dataH1,meshPositionsFieldName); CHKERRQ(ierr);
    try {
      ierr = opHOCoordsAndNormals.opRhs(dataH1); CHKERRQ(ierr);
      ierr = opHOCoordsAndNormals.calculateNormals(); CHKERRQ(ierr);
    } catch (std::exception& ex) {
      std::ostringstream ss;
      ss << "thorw in method: " << ex.what() << " at line " << __LINE__ << " in file " << __FILE__;
      SETERRQ(PETSC_COMM_SELF,1,ss.str().c_str());
    }
  } else {
    hoCoordsAtGaussPts.resize(0,0,false);
    normalsAtGaussPt.resize(0,0,false);
    tangentOneAtGaussPt.resize(0,0,false);
    tangentTwoAtGaussPt.resize(0,0,false);
  }
  PetscFunctionReturn(0);
}


PetscErrorCode FaceElementForcesAndSourcesCore::operator()() {
  PetscFunctionBegin;

  if(numeredEntFiniteElementPtr->getEntType() != MBTRI) PetscFunctionReturn(0);

  // Calculate normal and tangent vectors for face geometry given by 3 nodes.
  ierr = calculateAreaAndNormal(); CHKERRQ(ierr);
  ierr = getSpaceBaseAndOrderOnElement(); CHKERRQ(ierr);

  ierr = setIntegartionPts(); CHKERRQ(ierr);
  if(nbGaussPts == 0) PetscFunctionReturn(0);

  dataH1.dataOnEntities[MBVERTEX][0].getDiffN(NOBASE).resize(3,2,false);
  ierr = ShapeDiffMBTRI(
    &*dataH1.dataOnEntities[MBVERTEX][0].getDiffN(NOBASE).data().begin()
  ); CHKERRQ(ierr);

  /// Use the some node base

  ierr = calculateCoordinatesAtGaussPts(); CHKERRQ(ierr);

  // Share base shape functions between spaces
  dataHdiv.dataOnEntities[MBVERTEX][0].getNSharedPtr(NOBASE) = dataH1.dataOnEntities[MBVERTEX][0].getNSharedPtr(NOBASE);
  dataHcurl.dataOnEntities[MBVERTEX][0].getNSharedPtr(NOBASE) = dataH1.dataOnEntities[MBVERTEX][0].getNSharedPtr(NOBASE);
  dataL2.dataOnEntities[MBVERTEX][0].getNSharedPtr(NOBASE) = dataH1.dataOnEntities[MBVERTEX][0].getNSharedPtr(NOBASE);
  dataHdiv.dataOnEntities[MBVERTEX][0].getDiffNSharedPtr(NOBASE) = dataH1.dataOnEntities[MBVERTEX][0].getDiffNSharedPtr(NOBASE);
  dataHcurl.dataOnEntities[MBVERTEX][0].getDiffNSharedPtr(NOBASE) = dataH1.dataOnEntities[MBVERTEX][0].getDiffNSharedPtr(NOBASE);
  dataL2.dataOnEntities[MBVERTEX][0].getDiffNSharedPtr(NOBASE) = dataH1.dataOnEntities[MBVERTEX][0].getDiffNSharedPtr(NOBASE);

  ierr = calculateBaseFunctionsOnElement(); CHKERRQ(ierr);
  ierr = calculateHoNormal(); CHKERRQ(ierr);

  // Apply Piola transform to HDiv and HCurl spaces, uses previously calculated
  // faces normal and tangent vectors.
  if(dataH1.spacesOnEntities[MBTRI].test(HDIV)) {
    ierr = opContravariantTransoform.opRhs(dataHdiv); CHKERRQ(ierr);
  }
  if(dataH1.spacesOnEntities[MBEDGE].test(HCURL)) {
    // cerr << dataHcurl.dataOnEntities[MBEDGE][0].getN(AINSWORTH_LEGENDRE_BASE) << endl;
    ierr = opCovariantTransoform.opRhs(dataHcurl); CHKERRQ(ierr);
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

    if(oit->sPace!=LASTSPACE) {

      // Set field
      switch(oit->sPace) {
        case NOSPACE:
        SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"unknown space");
        case H1:
        op_data[0] = &dataH1;
        break;
        case HCURL:
        op_data[0] = &dataHcurl;
        break;
        case HDIV:
        op_data[0] = &dataHdiv;
        break;
        case L2:
        op_data[0] = &dataL2;
        break;
        case NOFIELD:
        op_data[0] = &dataNoField;
        break;
        case LASTSPACE:
        SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"unknown space");
      }

      // Reseat all data which all field dependent
      op_data[0]->resetFieldDepenentData();
      last_eval_field_name[0] = "";
      last_eval_field_name[1] = "";

      // Run operator
      ierr = oit->opRhs(
        *op_data[0],
        oit->doVertices,
        oit->doEdges,
        oit->doQuads,
        oit->doTris,
        false,
        false
      ); CHKERRQ(ierr);

    } else {

      for(int ss = 0;ss!=2;ss++) {

        std::string field_name = !ss ? oit->rowFieldName : oit->colFieldName;
        const Field* field_struture = mField.get_field_structure(field_name);
        BitFieldId data_id = field_struture->getId();

        if((oit->getNumeredEntFiniteElementPtr()->getBitFieldIdData()&data_id).none()) {
          SETERRQ2(
            PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,
            "no data field < %s > on finite element < %s >",
            field_name.c_str(),feName.c_str()
          );
        }

        if(oit->getOpType()&types[ss] || oit->getOpType()&UserDataOperator::OPROWCOL) {

          space[ss] = field_struture->getSpace();
          switch(space[ss]) {
            case NOSPACE:
            SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"unknown space");
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
          }

          base[ss] = field_struture->getApproxBase();
          switch(base[ss]) {
            case NOBASE:
            case AINSWORTH_LEGENDRE_BASE:
            case AINSWORTH_LOBBATO_BASE:
            case DEMKOWICZ_JACOBI_BASE:
            break;
            default:
            SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"unknown or not implemented base");
            break;
          }

          if(last_eval_field_name[ss]!=field_name) {

            switch(space[ss]) {
              case NOSPACE:
              SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"unknown space");
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
              case L2:
              if(!ss) {
                ierr = getTrisRowIndices(*op_data[ss],field_name); CHKERRQ(ierr);
              } else {
                ierr = getTrisColIndices(*op_data[ss],field_name); CHKERRQ(ierr);
              }
              ierr = getTrisDataOrderSpaceAndBase(*op_data[ss],field_name); CHKERRQ(ierr);
              ierr = getTrisFieldData(*op_data[ss],field_name); CHKERRQ(ierr);
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
            }
            last_eval_field_name[ss]=field_name;

          }
        }
      }

      if(oit->getOpType()&UserDataOperator::OPROW) {
        try {
          ierr = oit->opRhs(
            *op_data[0],
            oit->doVertices,
            oit->doEdges,
            oit->doQuads,
            oit->doTris,
            false,
            false
          ); CHKERRQ(ierr);
        } catch (std::exception& ex) {
          std::ostringstream ss;
          ss << "Operator " << typeid(*oit).name() //boost::core::demangle(typeid(*oit).name())
          << " operator number " << std::distance<boost::ptr_vector<UserDataOperator>::iterator>(opPtrVector.begin(),oit)
          << " thorw in method: " << ex.what()
          << " at line " << __LINE__
          << " in file " << __FILE__;
          SETERRQ(PETSC_COMM_SELF,MOFEM_STD_EXCEPTION_THROW,ss.str().c_str());
        }
      }

      if(oit->getOpType()&UserDataOperator::OPCOL) {
        try {
          ierr = oit->opRhs(
            *op_data[1],
            oit->doVertices,
            oit->doEdges,
            oit->doQuads,
            oit->doTris,
            false,
            false
          ); CHKERRQ(ierr);
        } catch (std::exception& ex) {
          std::ostringstream ss;
          ss << "Operator " << typeid(*oit).name() //boost::core::demangle(typeid(*oit).name())
          << " operator number " << std::distance<boost::ptr_vector<UserDataOperator>::iterator>(opPtrVector.begin(),oit)
          << " thorw in method: " << ex.what()
          << " at line " << __LINE__
          << " in file " << __FILE__;
          SETERRQ(PETSC_COMM_SELF,MOFEM_STD_EXCEPTION_THROW,ss.str().c_str());
        }
      }

      if(oit->getOpType()&UserDataOperator::OPROWCOL) {
        try {
          ierr = oit->opLhs(*op_data[0],*op_data[1],oit->sYmm); CHKERRQ(ierr);
        } catch (std::exception& ex) {
          std::ostringstream ss;
          ss << "Operator " << typeid(*oit).name() //boost::core::demangle(typeid(*oit).name())
          << " operator number " << std::distance<boost::ptr_vector<UserDataOperator>::iterator>(opPtrVector.begin(),oit)
          << " thorw in method: " << ex.what()
          << " at line " << __LINE__
          << " in file " << __FILE__;
          SETERRQ(PETSC_COMM_SELF,MOFEM_STD_EXCEPTION_THROW,ss.str().c_str());
        }
      }

    }

  }

  PetscFunctionReturn(0);
}

PetscErrorCode OpCalculateInvJacForFace::doWork(
  int side,
  EntityType type,
  DataForcesAndSurcesCore::EntData &data
) {
  
  PetscFunctionBegin;

  if(getNumeredEntFiniteElementPtr()->getEntType()!=MBTRI) {
    SETERRQ(
      PETSC_COMM_SELF,
      MOFEM_DATA_INCONSISTENCY,
      "This operator can be used only with element which is triangle"
    );
  }

  try {

    if(type == MBVERTEX) {
      VectorDouble &coords = getCoords();
      double *coords_ptr = &*coords.data().begin();
      double diff_n[6];
      ierr = ShapeDiffMBTRI(diff_n); CHKERRQ(ierr);
      double j00,j01,j10,j11;
      for(int gg = 0;gg<1;gg++) {
        // this is triangle, derivative of nodal shape functions is constant.
        // So only need to do one node.
        j00 = cblas_ddot(3,&coords_ptr[0],3,&diff_n[0],2);
        j01 = cblas_ddot(3,&coords_ptr[0],3,&diff_n[1],2);
        j10 = cblas_ddot(3,&coords_ptr[1],3,&diff_n[0],2);
        j11 = cblas_ddot(3,&coords_ptr[1],3,&diff_n[1],2);
      }
      double det = j00*j11-j01*j10;
      invJac.resize(2,2,false);
      invJac(0,0) = j11/det;
      invJac(0,1) = -j01/det;
      invJac(1,0) = -j10/det;
      invJac(1,1) = j00/det;
    }
  } catch (std::exception& ex) {
    std::ostringstream ss;
    ss << "thorw in method: " << ex.what() << " at line " << __LINE__ << " in file " << __FILE__;
    SETERRQ(PETSC_COMM_SELF,MOFEM_STD_EXCEPTION_THROW,ss.str().c_str());
  }

  doVertices = true;
  doEdges = false;
  doQuads = false;
  doTris = false;
  doTets= false;
  doPrisms = false;

  PetscFunctionReturn(0);
}

PetscErrorCode OpSetInvJacH1ForFace::doWork(
  int side,
  EntityType type,
  DataForcesAndSurcesCore::EntData &data
) {
  PetscFunctionBegin;
  // 

  if(
    getNumeredEntFiniteElementPtr()->getEntType()!=MBTRI &&
    getNumeredEntFiniteElementPtr()->getEntType()!=MBQUAD
  ) {
    SETERRQ(
      PETSC_COMM_SELF,
      MOFEM_DATA_INCONSISTENCY,
      "This operator can be used only with element which is triangle"
    );
  }


  for(int b = AINSWORTH_LEGENDRE_BASE; b!=USER_BASE; b++) {

    FieldApproximationBase base = ApproximationBaseArray[b];

    try {

      unsigned int nb_dofs = data.getN(base).size2();
      if(nb_dofs==0) PetscFunctionReturn(0);
      unsigned int nb_gauss_pts = data.getN(base).size1();
      diffNinvJac.resize(nb_gauss_pts,2*nb_dofs,false);

      if(type!=MBVERTEX) {
        if(nb_dofs != data.getDiffN(base).size2()/2) {
          SETERRQ2(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,
            "data inconsistency nb_dofs != data.diffN.size2()/2 ( %u != %u/2 )",
            nb_dofs,data.getDiffN(base).size2()
          );
        }
      }

      //std::cerr << type << std::endl;
      //std::cerr << data.getDiffN() << std::endl;
      //std::cerr << std::endl;
      switch (type) {
        case MBVERTEX:
        case MBEDGE:
        case MBTRI: {
          for(unsigned int gg = 0;gg<nb_gauss_pts;gg++) {
            for(unsigned int dd = 0;dd<nb_dofs;dd++) {
              cblas_dgemv(
                CblasRowMajor,CblasTrans,2,2,1,
                &*invJac.data().begin(),2,
                &data.getDiffN(base)(gg,2*dd),1,
                0,&diffNinvJac(gg,2*dd),1
              );
            }
          }
          data.getDiffN(base).data().swap(diffNinvJac.data());
        }
        break;
        default:
        SETERRQ(PETSC_COMM_SELF,MOFEM_NOT_IMPLEMENTED,"not implemented");
      }

    } catch (std::exception& ex) {
      std::ostringstream ss;
      ss << "Error in OpSetInvJacH1ForFace "
      << "thorw in method: "
      << ex.what()
      << " at line " << __LINE__ << " in file " << __FILE__;
      SETERRQ(PETSC_COMM_SELF,MOFEM_STD_EXCEPTION_THROW,ss.str().c_str());
    }

  }

  PetscFunctionReturn(0);
}

PetscErrorCode OpSetInvJacHcurlFace::doWork(
  int side,EntityType type,DataForcesAndSurcesCore::EntData &data
) {
  PetscFunctionBegin;

  if(type != MBEDGE && type != MBTRI) PetscFunctionReturn(0);

  if(
    getNumeredEntFiniteElementPtr()->getEntType()!=MBTRI
  ) {
    SETERRQ(
      PETSC_COMM_SELF,
      MOFEM_DATA_INCONSISTENCY,
      "This operator can be used only with element which is triangle"
    );
  }

  FTensor::Tensor2<double*,2,2> t_inv_jac = FTensor::Tensor2<double*,2,2>(
    &invJac(0,0),&invJac(0,1),&invJac(1,0),&invJac(1,1)
  );

  FTensor::Index<'i',3> i;
  FTensor::Index<'j',2> j;
  FTensor::Index<'k',2> k;

  for(int b = AINSWORTH_LEGENDRE_BASE; b!=USER_BASE; b++) {

    FieldApproximationBase base = ApproximationBaseArray[b];

    try {

      const unsigned int nb_base_functions = data.getDiffHcurlN(base).size2()/6;
      if(!nb_base_functions) continue;
      const unsigned int nb_gauss_pts = data.getDiffHcurlN(base).size1();

      diffHcurlInvJac.resize(
        nb_gauss_pts,data.getDiffHcurlN(base).size2(),false
      );

      // cerr << data.getDiffHcurlN(base) << endl;

      FTensor::Tensor2<double*,3,2> t_diff_n = data.getFTensor2DiffHcurlN<3,2>(base);
      double *t_inv_diff_n_ptr = &*diffHcurlInvJac.data().begin();
      FTensor::Tensor2<double*,3,2> t_inv_diff_n(
        t_inv_diff_n_ptr,
        &t_inv_diff_n_ptr[HCURL0_1],
        &t_inv_diff_n_ptr[HCURL1_0],
        &t_inv_diff_n_ptr[HCURL1_1],
        &t_inv_diff_n_ptr[HCURL2_0],
        &t_inv_diff_n_ptr[HCURL2_1],6
      );

      for(unsigned int gg = 0;gg!=nb_gauss_pts;gg++) {
        for(unsigned int bb = 0;bb!=nb_base_functions;bb++) {
          t_inv_diff_n(i,j) = t_diff_n(i,k)*t_inv_jac(k,j);
          ++t_diff_n;
          ++t_inv_diff_n;
        }
      }

      data.getDiffHcurlN(base).data().swap(diffHcurlInvJac.data());


    } catch (std::exception& ex) {
      std::ostringstream ss;
      ss << "thorw in method: " << ex.what() << " at line " << __LINE__ << " in file " << __FILE__;
      SETERRQ(PETSC_COMM_SELF,MOFEM_STD_EXCEPTION_THROW,ss.str().c_str());
    }

  }

  PetscFunctionReturn(0);
}

}
