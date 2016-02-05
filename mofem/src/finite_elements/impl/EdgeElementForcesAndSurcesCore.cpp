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
// #include <version.h>
#include <definitions.h>
#include <Common.hpp>

#include <h1_hdiv_hcurl_l2.h>
#include <fem_tools.h>

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
#include <FieldInterface.hpp>
#include <MeshRefinment.hpp>
#include <PrismInterface.hpp>
#include <SeriesRecorder.hpp>
#include <Core.hpp>

#include <DataStructures.hpp>
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

  if(fePtr->get_ent_type() != MBEDGE) PetscFunctionReturn(0);

  //PetscAttachDebugger();

  ierr = getEdgesOrder(dataH1,H1); CHKERRQ(ierr);

  dataH1.dataOnEntities[MBEDGE][0].getSense() = 1; // set sense to 1, this is this entity
  int order = dataH1.dataOnEntities[MBEDGE][0].getOrder();
  int rule = getRule(order);
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
    } else {
      SETERRQ2(
        PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"rule > quadrature order %d < %d",
        rule,QUAD_1D_TABLE_SIZE
      );
      nb_gauss_pts = 0;
    }
  }
  // {
  //   nb_gauss_pts = gm_rule_size(rule,1);
  //   gaussPts.resize(2,nb_gauss_pts,false);
  //   ierr = Grundmann_Moeller_integration_points_1D_EDGE(
  //     rule,&gaussPts(0,0),&gaussPts(1,0)
  //   ); CHKERRQ(ierr);
  // }

  ierr = shapeEDGEFunctions_H1(dataH1,0,&gaussPts(0,0),nb_gauss_pts); CHKERRQ(ierr);

  EntityHandle ent = fePtr->get_ent();
  int num_nodes;
  const EntityHandle* conn;
  rval = mField.get_moab().get_connectivity(ent,conn,num_nodes,true); CHKERR_PETSC(rval);
  cOords.resize(num_nodes*3,false);
  rval = mField.get_moab().get_coords(conn,num_nodes,&*cOords.data().begin()); CHKERR_PETSC(rval);

  dIrection.resize(3,false);
  cblas_dcopy(3,&cOords[3],1,&*dIrection.data().begin(),1);
  cblas_daxpy(3,-1.,&cOords[0],1,&*dIrection.data().begin(),1);
  lEngth = cblas_dnrm2(3,&*dIrection.data().begin(),1);

  coordsAtGaussPts.resize(nb_gauss_pts,3,false);
  for(int gg = 0;gg<nb_gauss_pts;gg++) {
    for(int dd = 0;dd<3;dd++) {
      coordsAtGaussPts(gg,dd) = N_MBEDGE0(gaussPts(0,gg))*cOords[dd] + N_MBEDGE1(gaussPts(0,gg))*cOords[3+dd];
    }
  }
  //cerr << coordsAtGaussPts << endl;

  if(
    dataPtr->get<FieldName_mi_tag>().find(meshPositionsFieldName)!=
    dataPtr->get<FieldName_mi_tag>().end()
  ) {

    ierr = getEdgesOrder(dataH1,meshPositionsFieldName); CHKERRQ(ierr);
    ierr = getNodesFieldData(dataH1,meshPositionsFieldName); CHKERRQ(ierr);
    ierr = getEdgesFieldData(dataH1,meshPositionsFieldName); CHKERRQ(ierr);
    try {
      ierr = opGetHoTangentOnEdge.opRhs(dataH1); CHKERRQ(ierr);
    } catch (exception& ex) {
      ostringstream ss;
      ss << "thorw in method: " << ex.what() << " at line " << __LINE__ << " in file " << __FILE__;
      SETERRQ(PETSC_COMM_SELF,1,ss.str().c_str());
    }
  } else {
    tAngent_at_GaussPt.resize(0,3,false);
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
          PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"no data field < %s > on finite element < %s >",
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
      } catch (exception& ex) {
        ostringstream ss;
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
