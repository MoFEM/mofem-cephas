/** \file EdgeElementForcesAndSurcesCore.cpp

\brief Implementation of edge element

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
#include <CubitBCData.hpp>
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
#include <EdgePolynomialBase.hpp> // Base functions on tet
#include <DataOperators.hpp>
#include <ElementsOnEntities.hpp>
#include <EdgeElementForcesAndSurcesCore.hpp>

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

PetscErrorCode EdgeElementForcesAndSurcesCore::operator()() {
  PetscFunctionBegin;

  if(numeredEntFiniteElementPtr->getEntType() != MBEDGE) PetscFunctionReturn(0);

  EntityHandle ent = numeredEntFiniteElementPtr->getEnt();
  {
    int num_nodes;
    const EntityHandle* conn;
    rval = mField.get_moab().get_connectivity(ent,conn,num_nodes,true); CHKERRQ_MOAB(rval);
    cOords.resize(num_nodes*3,false);
    rval = mField.get_moab().get_coords(conn,num_nodes,&*cOords.data().begin()); CHKERRQ_MOAB(rval);
    dIrection.resize(3,false);
    cblas_dcopy(3,&cOords[3],1,&*dIrection.data().begin(),1);
    cblas_daxpy(3,-1.,&cOords[0],1,&*dIrection.data().begin(),1);
    lEngth = cblas_dnrm2(3,&*dIrection.data().begin(),1);
  }

  //PetscAttachDebugger();
  ierr = getSpacesAndBaseOnEntities(dataH1); CHKERRQ(ierr);
  ierr = getEdgesDataOrder(dataH1,H1); CHKERRQ(ierr);
  ierr = getEdgesDataOrder(dataH1,H1); CHKERRQ(ierr);

  dataH1.dataOnEntities[MBEDGE][0].getSense() = 1; // set sense to 1, this is this entity

  int order_data = getMaxDataOrder();
  int order_row = getMaxRowOrder();
  int order_col = getMaxColOrder();
  int rule = getRule(order_row,order_col,order_data);
  int nb_gauss_pts;
  {
    if(rule<QUAD_1D_TABLE_SIZE) {
      if(QUAD_1D_TABLE[rule]->dim!=1) {
        SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"wrong dimension");
      }
      if(QUAD_1D_TABLE[rule]->order<rule) {
        SETERRQ2(
          PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"wrong order %d != %d",
          QUAD_1D_TABLE[rule]->order,rule
        );
      }
      nb_gauss_pts = QUAD_1D_TABLE[rule]->npoints;
      gaussPts.resize(2,nb_gauss_pts,false);
      cblas_dcopy(
        nb_gauss_pts,&QUAD_1D_TABLE[rule]->points[1],2,&gaussPts(0,0),1
      );
      cblas_dcopy(
        nb_gauss_pts,QUAD_1D_TABLE[rule]->weights,1,&gaussPts(1,0),1
      );
      dataH1.dataOnEntities[MBVERTEX][0].getN(NOBASE).resize(nb_gauss_pts,2,false);
      double *shape_ptr = &*dataH1.dataOnEntities[MBVERTEX][0].getN(NOBASE).data().begin();
      cblas_dcopy(
        2*nb_gauss_pts,QUAD_1D_TABLE[rule]->points,1,shape_ptr,1
      );
    } else {
      SETERRQ2(
        PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"rule > quadrature order %d < %d",
        rule,QUAD_1D_TABLE_SIZE
      );
      nb_gauss_pts = 0;
    }
  }

  coordsAtGaussPts.resize(nb_gauss_pts,3,false);
  for(int gg = 0;gg<nb_gauss_pts;gg++) {
    for(int dd = 0;dd<3;dd++) {
      coordsAtGaussPts(gg,dd) = N_MBEDGE0(gaussPts(0,gg))*cOords[dd] + N_MBEDGE1(gaussPts(0,gg))*cOords[3+dd];
    }
  }

  try {

    for(int b = AINSWORTH_COLE_BASE;b!=LASTBASE;b++) {
      if(dataH1.bAse.test(b)) {
        switch (ApproximationBaseArray[b]) {
          case AINSWORTH_COLE_BASE:
          case LOBATTO_BASE:
          if(dataH1.spacesOnEntities[MBVERTEX].test(H1)) {
            ierr = EdgePolynomialBase().getValue(
              gaussPts,
              boost::shared_ptr<BaseFunctionCtx>(
                new EntPolynomialBaseCtx(dataH1,H1,ApproximationBaseArray[b],NOBASE)
              )
            ); CHKERRQ(ierr);
          }
          if(dataH1.spacesOnEntities[MBEDGE].test(HCURL)) {
            ierr = EdgePolynomialBase().getValue(
              gaussPts,
              boost::shared_ptr<BaseFunctionCtx>(
                new EntPolynomialBaseCtx(dataH1,HCURL,ApproximationBaseArray[b],NOBASE)
              )
            ); CHKERRQ(ierr);
          }
          if(dataH1.spacesOnEntities[MBEDGE].test(HDIV)) {
            ierr = EdgePolynomialBase().getValue(
              gaussPts,
              boost::shared_ptr<BaseFunctionCtx>(
                new EntPolynomialBaseCtx(dataH1,HDIV,ApproximationBaseArray[b],NOBASE)
              )
            ); CHKERRQ(ierr);
          }
          if(dataH1.spacesOnEntities[MBEDGE].test(L2)) {
            ierr = EdgePolynomialBase().getValue(
              gaussPts,
              boost::shared_ptr<BaseFunctionCtx>(
                new EntPolynomialBaseCtx(dataH1,L2,ApproximationBaseArray[b],NOBASE)
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

  if(
    dataPtr->get<FieldName_mi_tag>().find(meshPositionsFieldName)!=
    dataPtr->get<FieldName_mi_tag>().end()
  ) {

    ierr = getEdgesDataOrderSpaceAndBase(dataH1,meshPositionsFieldName); CHKERRQ(ierr);
    ierr = getEdgesFieldData(dataH1,meshPositionsFieldName); CHKERRQ(ierr);
    ierr = getNodesFieldData(dataH1,meshPositionsFieldName); CHKERRQ(ierr);
    try {
      ierr = opGetHoTangentOnEdge.opRhs(dataH1); CHKERRQ(ierr);
    } catch (std::exception& ex) {
      std::ostringstream ss;
      ss << "thorw in method: " << ex.what() << " at line " << __LINE__ << " in file " << __FILE__;
      SETERRQ(PETSC_COMM_SELF,1,ss.str().c_str());
    }
  } else {
    tAngent_at_GaussPt.resize(0,3,false);
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
          case H1:
          op_data[ss] = !ss ? &dataH1 : &derivedDataH1;
          break;
          case HCURL:
          SETERRQ(PETSC_COMM_SELF,MOFEM_NOT_IMPLEMENTED,"not implemented yet");
          break;
          case HDIV:
          SETERRQ(PETSC_COMM_SELF,MOFEM_NOT_IMPLEMENTED,"not make sanes on edge");
          break;
          case L2:
          SETERRQ(PETSC_COMM_SELF,MOFEM_NOT_IMPLEMENTED,"not make sanes on edge");
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
            break;
            case HDIV:
            SETERRQ(PETSC_COMM_SELF,MOFEM_NOT_IMPLEMENTED,"not make sanes on edge");
            break;
            case L2:
            SETERRQ(PETSC_COMM_SELF,MOFEM_NOT_IMPLEMENTED,"not make sanes on edge");
            break;
            case NOFIELD:
            if(!getNinTheLoop()) {
              // NOFIELD data are the same for each element, can be retreived only once
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
          false,
          false,
          false,
          false
        ); CHKERRQ(ierr);
      } catch (std::exception& ex) {
        std::ostringstream ss;
        ss << "thorw in method: " << ex.what() << " at line " << __LINE__ << " in file " << __FILE__;
        SETERRQ(PETSC_COMM_SELF,MOFEM_STD_EXCEPTION_THROW,ss.str().c_str());
      }
    }


    if(oit->getOpType()&UserDataOperator::OPCOL) {
      try {
        ierr = oit->opRhs(
          *op_data[1],
          oit->doVerticesCol,
          oit->doEdgesCol,
          false,
          false,
          false,
          false
        ); CHKERRQ(ierr);
      } catch (std::exception& ex) {
        std::ostringstream ss;
        ss << "thorw in method: " << ex.what() << " at line " << __LINE__ << " in file " << __FILE__;
        SETERRQ(PETSC_COMM_SELF,MOFEM_STD_EXCEPTION_THROW,ss.str().c_str());
      }
    }


    if(oit->getOpType()&UserDataOperator::OPROWCOL) {
      try {
        ierr = oit->opLhs(*op_data[0],*op_data[1],oit->sYmm); CHKERRQ(ierr);
      } catch (std::exception& ex) {
        std::ostringstream ss;
        ss << "thorw in method: " << ex.what() << " at line " << __LINE__ << " in file " << __FILE__;
        SETERRQ(PETSC_COMM_SELF,MOFEM_STD_EXCEPTION_THROW,ss.str().c_str());
      }
    }

  }

  PetscFunctionReturn(0);
}

}
