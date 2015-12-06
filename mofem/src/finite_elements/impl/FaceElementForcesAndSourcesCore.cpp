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
// #include <version.h>
#include <definitions.h>
#include <Common.hpp>

#include <h1_hdiv_hcurl_l2.h>
#include <fem_tools.h>

#include <MaterialBlocks.hpp>
#include <CubitBCData.hpp>
#include <TagMultiIndices.hpp>
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
#include <FieldInterface.hpp>
#include <MeshRefinment.hpp>
#include <PrismInterface.hpp>
#include <SeriesRecorder.hpp>
#include <Core.hpp>

#include <DataStructures.hpp>
#include <DataOperators.hpp>
#include <ElementsOnEntities.hpp>
#include <FaceElementForcesAndSourcesCore.hpp>

#ifdef __cplusplus
extern "C" {
#endif
  #include <cblas.h>
  #include <lapack_wrap.h>
  #include <gm_rule.h>
#ifdef __cplusplus
}
#endif

namespace MoFEM {

PetscErrorCode FaceElementForcesAndSourcesCore::operator()() {
  PetscFunctionBegin;

  if(fePtr->get_ent_type() != MBTRI) PetscFunctionReturn(0);

  ierr = getSpacesOnEntities(dataH1); CHKERRQ(ierr);

  //H1
  if((dataH1.spacesOnEntities[MBEDGE]).test(H1)) {
    ierr = getEdgesSense(dataH1); CHKERRQ(ierr);
    ierr = getEdgesOrder(dataH1,H1); CHKERRQ(ierr);
  }
  if((dataH1.spacesOnEntities[MBTRI]).test(H1)) {
    ierr = getTrisSense(dataH1); CHKERRQ(ierr);
    ierr = getTrisOrder(dataH1,H1); CHKERRQ(ierr);
  }

  //Hdiv
  if((dataH1.spacesOnEntities[MBTRI]).test(HDIV)) {
    ierr = getTrisSense(dataHdiv); CHKERRQ(ierr);
    ierr = getTrisOrder(dataHdiv,HDIV); CHKERRQ(ierr);
  }

  int order = 1;
  for(unsigned int ee = 0;ee<dataH1.dataOnEntities[MBEDGE].size();ee++) {
    order = max(order,dataH1.dataOnEntities[MBEDGE][ee].getOrder());
  }
  for(unsigned int ff = 0;ff<dataH1.dataOnEntities[MBTRI].size();ff++) {
    order = max(order,dataH1.dataOnEntities[MBTRI][ff].getOrder());
  }
  for(unsigned int ff = 0;ff<dataHdiv.dataOnEntities[MBTRI].size();ff++) {
    order = max(order,dataHdiv.dataOnEntities[MBTRI][ff].getOrder());
  }

  int nb_gauss_pts;
  int rule = getRule(order);
  if(rule >= 0) {
    //if(mField.check_field(meshPositionsFieldName)) {
    //rule += 1;
    //}
    nb_gauss_pts = gm_rule_size(rule,2);
    gaussPts.resize(3,nb_gauss_pts,false);
    ierr = Grundmann_Moeller_integration_points_2D_TRI(
      rule,&gaussPts(0,0),&gaussPts(1,0),&gaussPts(2,0)
    ); CHKERRQ(ierr);
  } else {
    ierr = setGaussPts(order); CHKERRQ(ierr);
    nb_gauss_pts = gaussPts.size2();
  }
  if(nb_gauss_pts == 0) PetscFunctionReturn(0);

  ierr = shapeTRIFunctions_H1(dataH1,&gaussPts(0,0),&gaussPts(1,0),nb_gauss_pts); CHKERRQ(ierr);
  if(dataH1.spacesOnEntities[MBTRI].test(HDIV)) {
    ierr = shapeTRIFunctions_Hdiv(
      dataHdiv,&gaussPts(0,0),&gaussPts(1,0),nb_gauss_pts); CHKERRQ(ierr
      ); CHKERRQ(ierr);
  }

  EntityHandle ent = fePtr->get_ent();
  rval = mField.get_moab().get_connectivity(ent,conn,num_nodes,true); CHKERR_PETSC(rval);
  coords.resize(num_nodes*3,false);
  rval = mField.get_moab().get_coords(conn,num_nodes,&*coords.data().begin()); CHKERR_PETSC(rval);

  normal.resize(3,false);
  ierr = ShapeFaceNormalMBTRI(
    &*dataH1.dataOnEntities[MBVERTEX][0].getDiffN().data().begin(),
    &*coords.data().begin(),&*normal.data().begin()
  ); CHKERRQ(ierr);
  aRea = cblas_dnrm2(3,&*normal.data().begin(),1)*0.5;

  coordsAtGaussPts.resize(nb_gauss_pts,3,false);
  for(int gg = 0;gg<nb_gauss_pts;gg++) {
    for(int dd = 0;dd<3;dd++) {
      coordsAtGaussPts(gg,dd) = cblas_ddot(3,&dataH1.dataOnEntities[MBVERTEX][0].getN()(gg,0),1,&coords[dd],3);
    }
  }

  // In linear geometry derivatives are constant,
  // this in expense of efficiency makes implementation
  // constant between vertices and other types of entities
  MatrixDouble diffN(nb_gauss_pts,6);
  for(int gg = 0;gg<nb_gauss_pts;gg++) {
    for(int nn = 0;nn<3;nn++) {
      for(int dd = 0;dd<2;dd++) {
        diffN(gg,nn*2+dd) = dataH1.dataOnEntities[MBVERTEX][0].getDiffN()(nn,dd);
      }
    }
  }
  dataH1.dataOnEntities[MBVERTEX][0].getDiffN().resize(diffN.size1(),diffN.size2(),false);
  dataH1.dataOnEntities[MBVERTEX][0].getDiffN().data().swap(diffN.data());

  if(
    dataPtr->get<FieldName_mi_tag>().find(meshPositionsFieldName)!=
    dataPtr->get<FieldName_mi_tag>().end()
  ) {
    ierr = getEdgesOrder(dataH1,meshPositionsFieldName); CHKERRQ(ierr);
    ierr = getTrisOrder(dataH1,meshPositionsFieldName); CHKERRQ(ierr);
    ierr = getNodesFieldData(dataH1,meshPositionsFieldName); CHKERRQ(ierr);
    ierr = getEdgesFieldData(dataH1,meshPositionsFieldName); CHKERRQ(ierr);
    ierr = getTrisFieldData(dataH1,meshPositionsFieldName); CHKERRQ(ierr);
    ierr = getNodesFieldDofs(dataH1,meshPositionsFieldName); CHKERRQ(ierr);
    ierr = getEdgesFieldDofs(dataH1,meshPositionsFieldName); CHKERRQ(ierr);
    ierr = getTrisFieldDofs(dataH1,meshPositionsFieldName); CHKERRQ(ierr);
    try {
      ierr = opHOCoordsAndNormals.opRhs(dataH1); CHKERRQ(ierr);
      ierr = opHOCoordsAndNormals.calculateNormals(); CHKERRQ(ierr);
    } catch (exception& ex) {
      ostringstream ss;
      ss << "thorw in method: " << ex.what() << " at line " << __LINE__ << " in file " << __FILE__;
      SETERRQ(PETSC_COMM_SELF,1,ss.str().c_str());
    }
  } else {
    hoCoordsAtGaussPts.resize(0,0,false);
    nOrmals_at_GaussPt.resize(0,0,false);
    tAngent1_at_GaussPt.resize(0,0,false);
    tAngent2_at_GaussPt.resize(0,0,false);
  }

  if(dataH1.spacesOnEntities[MBTRI].test(HDIV)) {
    ierr = opSetPiolaTransoformOnTriangle.opRhs(dataHdiv); CHKERRQ(ierr);
  }

  const UserDataOperator::OpType types[2] = {
    UserDataOperator::OPROW, UserDataOperator::OPCOL
  };
  vector<string> last_eval_field_name(2);
  DataForcesAndSurcesCore *op_data[2];
  FieldSpace space[2];

  boost::ptr_vector<UserDataOperator>::iterator oit,hi_oit;
  oit = opPtrVector.begin();
  hi_oit = opPtrVector.end();

  for(;oit!=hi_oit;oit++) {

    oit->setPtrFE(this);

    for(int ss = 0;ss!=2;ss++) {

      string field_name = !ss ? oit->rowFieldName : oit->colFieldName;
      BitFieldId data_id = mField.get_field_structure(field_name)->get_id();
      if((oit->getMoFEMFEPtr()->get_BitFieldId_data()&data_id).none()) {
        SETERRQ2(
          PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,
          "no data field < %s > on finite element < %s >",
          field_name.c_str(),feName.c_str()
        );
      }

      if(oit->getOpType()&types[ss] || oit->getOpType()&UserDataOperator::OPROWCOL) {

        space[ss] = mField.get_field_structure(field_name)->get_space();

        switch(space[ss]) {
          case H1:
          op_data[ss] = !ss ? &dataH1 : &derivedDataH1;
          break;
          case HCURL:
          SETERRQ(PETSC_COMM_SELF,MOFEM_NOT_IMPLEMENTED,"not implemented yet");
          break;
          case HDIV:
          op_data[ss] = !ss ? &dataHdiv : &derivedDataHdiv;
          break;
          case L2:
          SETERRQ(PETSC_COMM_SELF,MOFEM_NOT_IMPLEMENTED,"not make sanes on face");
          break;
          case NOFIELD:
          op_data[ss] = !ss ? &dataNoField : &dataNoFieldCol;
          break;
          case LASTSPACE:
          SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"unknown space");
          break;
        }

        if(last_eval_field_name[ss]!=field_name) {

          switch(space[ss]) {
            case H1:
            if(!ss) {
              ierr = getRowNodesIndices(*op_data[ss],field_name); CHKERRQ(ierr);
            } else {
              ierr = getColNodesIndices(*op_data[ss],field_name); CHKERRQ(ierr);
            }
            ierr = getNodesFieldData(*op_data[ss],field_name); CHKERRQ(ierr);
            ierr = getNodesFieldDofs(*op_data[ss],field_name); CHKERRQ(ierr);
            case HCURL:
            if(!ss) {
              ierr = getEdgesRowIndices(*op_data[ss],field_name); CHKERRQ(ierr);
            } else {
              ierr = getEdgesColIndices(*op_data[ss],field_name); CHKERRQ(ierr);
            }
            ierr = getEdgesOrder(*op_data[ss],field_name); CHKERRQ(ierr);
            ierr = getEdgesFieldData(*op_data[ss],field_name); CHKERRQ(ierr);
            ierr = getEdgesFieldDofs(*op_data[ss],field_name); CHKERRQ(ierr);
            case HDIV:
            if(!ss) {
              ierr = getTrisRowIndices(*op_data[ss],field_name); CHKERRQ(ierr);
            } else {
              ierr = getTrisColIndices(*op_data[ss],field_name); CHKERRQ(ierr);
            }
            ierr = getTrisOrder(*op_data[ss],field_name); CHKERRQ(ierr);
            ierr = getTrisFieldData(*op_data[ss],field_name); CHKERRQ(ierr);
            ierr = getTrisFieldDofs(*op_data[ss],field_name); CHKERRQ(ierr);
            break;
            case L2:
            SETERRQ(PETSC_COMM_SELF,MOFEM_NOT_IMPLEMENTED,"not make sanes on face");
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
        ierr = oit->opRhs(*op_data[0]); CHKERRQ(ierr);
      } catch (exception& ex) {
        ostringstream ss;
        ss << "thorw in method: " << ex.what() << " at line " << __LINE__ << " in file " << __FILE__;
        SETERRQ(PETSC_COMM_SELF,MOFEM_STD_EXCEPTION_THROW,ss.str().c_str());
      }
    }


    if(oit->getOpType()&UserDataOperator::OPCOL) {
      try {
        ierr = oit->opRhs(*op_data[1]); CHKERRQ(ierr);
      } catch (exception& ex) {
        ostringstream ss;
        ss << "thorw in method: " << ex.what() << " at line " << __LINE__ << " in file " << __FILE__;
        SETERRQ(PETSC_COMM_SELF,MOFEM_STD_EXCEPTION_THROW,ss.str().c_str());
      }
    }


    if(oit->getOpType()&UserDataOperator::OPROWCOL) {
      try {
        ierr = oit->opLhs(*op_data[0],*op_data[1],oit->sYmm); CHKERRQ(ierr);
      } catch (exception& ex) {
        ostringstream ss;
        ss << "thorw in method: " << ex.what() << " at line " << __LINE__ << " in file " << __FILE__;
        SETERRQ(PETSC_COMM_SELF,MOFEM_STD_EXCEPTION_THROW,ss.str().c_str());
      }
    }

  }

  PetscFunctionReturn(0);
}

}
