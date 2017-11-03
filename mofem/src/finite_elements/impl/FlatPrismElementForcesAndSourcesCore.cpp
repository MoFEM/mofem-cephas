/** \file FlatPrismElementForcesAndSourcesCore.cpp

\brief Implementation of flat prism element

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
#include <ForcesAndSourcesCore.hpp>
#include <BaseFunction.hpp>
#include <EntPolynomialBaseCtx.hpp>
#include <FlatPrismPolynomialBase.hpp>
#include <FlatPrismElementForcesAndSourcesCore.hpp>

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

MoFEMErrorCode FlatPrismElementForcesAndSourcesCore::operator()() {
  MoFEMFunctionBeginHot;

  if(numeredEntFiniteElementPtr->getEntType() != MBPRISM) MoFEMFunctionReturnHot(0);

  try {

    EntityHandle ent = numeredEntFiniteElementPtr->getEnt();
    int num_nodes;
    const EntityHandle* conn;
    rval = mField.get_moab().get_connectivity(ent,conn,num_nodes,true); CHKERRQ_MOAB(rval);
    {
      coords.resize(num_nodes*3,false);
      rval = mField.get_moab().get_coords(conn,num_nodes,&*coords.data().begin()); CHKERRQ_MOAB(rval);

      double diff_n[6];
      ierr = ShapeDiffMBTRI(diff_n); CHKERRQ(ierr);
      normal.resize(6,false);
      ierr = ShapeFaceNormalMBTRI(diff_n,&coords[0],&normal[0]); CHKERRQ(ierr);
      ierr = ShapeFaceNormalMBTRI(diff_n, &coords[9],&normal[3]); CHKERRQ(ierr);
      aRea[0] = cblas_dnrm2(3,&normal[0],1)*0.5;
      aRea[1] = cblas_dnrm2(3,&normal[3],1)*0.5;
    }

    ierr = getSpacesAndBaseOnEntities(dataH1); CHKERRQ(ierr);

    //H1
    if((dataH1.spacesOnEntities[MBEDGE]).test(H1)) {
      ierr = getEdgesSense(dataH1); CHKERRQ(ierr);
      ierr = getEdgesDataOrder(dataH1,H1); CHKERRQ(ierr);
      ierr = getTrisSense(dataH1); CHKERRQ(ierr);
      ierr = getTrisDataOrder(dataH1,H1); CHKERRQ(ierr);
    }

    //H1
    if((dataH1.spacesOnEntities[MBEDGE]).test(HCURL)) {
      ierr = getEdgesSense(dataHcurl); CHKERRQ(ierr);
      ierr = getEdgesDataOrder(dataHcurl,HCURL); CHKERRQ(ierr);
      ierr = getTrisSense(dataHcurl); CHKERRQ(ierr);
      ierr = getTrisDataOrder(dataHcurl,HCURL); CHKERRQ(ierr);
    }

    //Hdiv
    if((dataH1.spacesOnEntities[MBTRI]).test(HDIV)) {
      ierr = getTrisSense(dataHdiv); CHKERRQ(ierr);
      ierr = getTrisDataOrder(dataHdiv,HDIV); CHKERRQ(ierr);
    }

    //Hdiv
    if((dataH1.spacesOnEntities[MBTRI]).test(L2)) {
      ierr = getTrisSense(dataL2); CHKERRQ(ierr);
      ierr = getTrisDataOrder(dataL2,L2); CHKERRQ(ierr);
    }

    int order_data = getMaxDataOrder();
    int order_row = getMaxRowOrder();
    int order_col = getMaxColOrder();
    int rule = getRule(order_row,order_col,order_data);
    int nb_gauss_pts;
    // int rule = getRule(order);
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
        nb_gauss_pts = QUAD_2D_TABLE[rule]->npoints;
        gaussPts.resize(3,nb_gauss_pts,false);
        cblas_dcopy(
          nb_gauss_pts,&QUAD_2D_TABLE[rule]->points[1],3,&gaussPts(0,0),1
        );
        cblas_dcopy(
          nb_gauss_pts,&QUAD_2D_TABLE[rule]->points[2],3,&gaussPts(1,0),1
        );
        cblas_dcopy(
          nb_gauss_pts,QUAD_2D_TABLE[rule]->weights,1,&gaussPts(2,0),1
        );
        dataH1.dataOnEntities[MBVERTEX][0].getN(NOBASE).resize(nb_gauss_pts,3,false);
        double *shape_ptr = &*dataH1.dataOnEntities[MBVERTEX][0].getN(NOBASE).data().begin();
        cblas_dcopy(
          3*nb_gauss_pts,QUAD_2D_TABLE[rule]->points,1,shape_ptr,1
        );
      } else {
        SETERRQ2(
          PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"rule > quadrature order %d < %d",
          rule,QUAD_2D_TABLE_SIZE
        );
        nb_gauss_pts = 0;
      }
    } else {
      ierr = setGaussPts(order_row,order_col,order_data); CHKERRQ(ierr);
      nb_gauss_pts = gaussPts.size2();
      dataH1.dataOnEntities[MBVERTEX][0].getN(NOBASE).resize(nb_gauss_pts,3,false);
      if(nb_gauss_pts) {
        ierr = ShapeMBTRI(
          &*dataH1.dataOnEntities[MBVERTEX][0].getN(NOBASE).data().begin(),
          &gaussPts(0,0),
          &gaussPts(1,0),
          nb_gauss_pts
        ); CHKERRQ(ierr);
      }
    }
    if(nb_gauss_pts == 0) MoFEMFunctionReturnHot(0);

    dataHdiv.dataOnEntities[MBVERTEX][0].getNSharedPtr(NOBASE) = dataH1.dataOnEntities[MBVERTEX][0].getNSharedPtr(NOBASE);
    dataHcurl.dataOnEntities[MBVERTEX][0].getNSharedPtr(NOBASE) = dataH1.dataOnEntities[MBVERTEX][0].getNSharedPtr(NOBASE);
    dataL2.dataOnEntities[MBVERTEX][0].getNSharedPtr(NOBASE) = dataH1.dataOnEntities[MBVERTEX][0].getNSharedPtr(NOBASE);
    {
      coordsAtGaussPts.resize(nb_gauss_pts,6,false);
      for(int gg = 0;gg<nb_gauss_pts;gg++) {
        for(int dd = 0;dd<3;dd++) {
          coordsAtGaussPts(gg,dd) = cblas_ddot(3,&dataH1.dataOnEntities[MBVERTEX][0].getN(NOBASE)(gg,0),1,&coords[dd],3);
          coordsAtGaussPts(gg,3+dd) = cblas_ddot(3,&dataH1.dataOnEntities[MBVERTEX][0].getN(NOBASE)(gg,0),1,&coords[9+dd],3);
        }
      }
    }

    try {

      std::vector<FieldApproximationBase> shape_functions_for_bases;
      for(int b = AINSWORTH_LEGENDRE_BASE;b!=LASTBASE;b++) {
        if(dataH1.bAse.test(b)) {
          switch (ApproximationBaseArray[b]) {
            case AINSWORTH_LEGENDRE_BASE:
            case AINSWORTH_LOBATTO_BASE:
            if(dataH1.spacesOnEntities[MBVERTEX].test(H1)) {
              ierr = FlatPrismPolynomialBase().getValue(
                gaussPts,
                boost::shared_ptr<BaseFunctionCtx>(
                  new FlatPrismPolynomialBaseCtx(
                    dataH1,mField.get_moab(),numeredEntFiniteElementPtr.get(),H1,ApproximationBaseArray[b],NOBASE
                  )
                )
              ); CHKERRQ(ierr);
            }
            if(dataH1.spacesOnEntities[MBTRI].test(HDIV)) {
              SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"Not yet implemented");
            }
            if(dataH1.spacesOnEntities[MBEDGE].test(HCURL)) {
              SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"Not yet implemented");
            }
            if(dataH1.spacesOnEntities[MBTET].test(L2)) {
              SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"Not yet implemented");
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

    try {

      if(
        dataPtr->get<FieldName_mi_tag>().find(meshPositionsFieldName)!=
        dataPtr->get<FieldName_mi_tag>().end()
      ) {
        hoCoordsAtGaussPtsF3.resize(nb_gauss_pts,3,false);
        nOrmals_at_GaussPtF3.resize(nb_gauss_pts,3,false);
        tAngent1_at_GaussPtF3.resize(nb_gauss_pts,3,false);
        tAngent2_at_GaussPtF3.resize(nb_gauss_pts,3,false);
        hoCoordsAtGaussPtsF4.resize(nb_gauss_pts,3,false);
        nOrmals_at_GaussPtF4.resize(nb_gauss_pts,3,false);
        tAngent1_at_GaussPtF4.resize(nb_gauss_pts,3,false);
        tAngent2_at_GaussPtF4.resize(nb_gauss_pts,3,false);
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
        hoCoordsAtGaussPtsF3.resize(0,0,false);
        nOrmals_at_GaussPtF3.resize(0,0,false);
        tAngent1_at_GaussPtF3.resize(0,0,false);
        tAngent2_at_GaussPtF3.resize(0,0,false);
        hoCoordsAtGaussPtsF4.resize(0,0,false);
        nOrmals_at_GaussPtF4.resize(0,0,false);
        tAngent1_at_GaussPtF4.resize(0,0,false);
        tAngent2_at_GaussPtF4.resize(0,0,false);
      }
    } catch (std::exception& ex) {
      std::ostringstream ss;
      ss << "thorw in method: " << ex.what() << " at line " << __LINE__ << " in file " << __FILE__;
      SETERRQ(PETSC_COMM_SELF,MOFEM_STD_EXCEPTION_THROW,ss.str().c_str());
    }

    if(dataH1.spacesOnEntities[MBTRI].test(HDIV)) {
      SETERRQ(PETSC_COMM_SELF,MOFEM_NOT_IMPLEMENTED,"not implemented yet");
    }

    const UserDataOperator::OpType types[2] = {
      UserDataOperator::OPROW, UserDataOperator::OPCOL
    };
    std::vector<std::string> last_eval_field_name(2);
    DataForcesAndSourcesCore *op_data[2];
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
            case AINSWORTH_LEGENDRE_BASE:
            case AINSWORTH_LOBATTO_BASE:
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
            oit->doVertices,
            oit->doEdges,
            oit->doQuads,
            oit->doTris,
            false,
            false
          ); CHKERRQ(ierr);
        } catch (std::exception& ex) {
          std::ostringstream ss;
          ss << "Operator "
             << boost::typeindex::type_id_runtime(*oit).pretty_name()
             << " operator number "
             << std::distance<boost::ptr_vector<UserDataOperator>::iterator>(
                    opPtrVector.begin(), oit)
             << " thorw in method: " << ex.what() << " at line " << __LINE__
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
          ss << "Operator "
             << boost::typeindex::type_id_runtime(*oit).pretty_name()
            << " operator number "
             << std::distance<boost::ptr_vector<UserDataOperator>::iterator>(
                    opPtrVector.begin(), oit)
             << " thorw in method: " << ex.what() << " at line " << __LINE__
             << " in file " << __FILE__;
          SETERRQ(PETSC_COMM_SELF,MOFEM_STD_EXCEPTION_THROW,ss.str().c_str());
        }
      }


      if(oit->getOpType()&UserDataOperator::OPROWCOL) {
        try {
          ierr = oit->opLhs(*op_data[0],*op_data[1],oit->sYmm); CHKERRQ(ierr);
        } catch (std::exception& ex) {
          std::ostringstream ss;
          ss << "Operator " 
          << boost::typeindex::type_id_runtime(*oit).pretty_name()
          << " operator number " << std::distance<boost::ptr_vector<UserDataOperator>::iterator>(opPtrVector.begin(),oit)
          << " thorw in method: " << ex.what()
          << " at line " << __LINE__
          << " in file " << __FILE__;
          SETERRQ(PETSC_COMM_SELF,MOFEM_STD_EXCEPTION_THROW,ss.str().c_str());
        }
      }

    }

  } catch (std::exception& ex) {
    std::ostringstream ss;
    ss << "thorw in method: " << ex.what() << " at line " << __LINE__ << " in file " << __FILE__;
    SETERRQ(PETSC_COMM_SELF,MOFEM_STD_EXCEPTION_THROW,ss.str().c_str());
  }

  MoFEMFunctionReturnHot(0);
}

MoFEMErrorCode OpCalculateInvJacForFlatPrism::doWork(
  int side,
  EntityType type,
  DataForcesAndSourcesCore::EntData &data
) {
  
  MoFEMFunctionBeginHot;

  try {

    if(type == MBVERTEX) {

      VectorDouble &coords = getCoords();
      double *coords_ptr = &*coords.data().begin();
      double diff_n[6];
      ierr = ShapeDiffMBTRI(diff_n); CHKERRQ(ierr);
      double j00_f3,j01_f3,j10_f3,j11_f3;
      for(int gg = 0;gg<1;gg++) {
        // this is triangle, derivative of nodal shape functions is constant.
        // So only need to do one node.
        j00_f3 = cblas_ddot(3,&coords_ptr[0],3,&diff_n[0],2);
        j01_f3 = cblas_ddot(3,&coords_ptr[0],3,&diff_n[1],2);
        j10_f3 = cblas_ddot(3,&coords_ptr[1],3,&diff_n[0],2);
        j11_f3 = cblas_ddot(3,&coords_ptr[1],3,&diff_n[1],2);
      }
      double det_f3 = j00_f3*j11_f3-j01_f3*j10_f3;
      invJacF3.resize(2,2,false);
      invJacF3(0,0) = j11_f3/det_f3;
      invJacF3(0,1) = -j01_f3/det_f3;
      invJacF3(1,0) = -j10_f3/det_f3;
      invJacF3(1,1) = j00_f3/det_f3;

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
  doTets = false;
  doPrisms = false;

  MoFEMFunctionReturnHot(0);
}

MoFEMErrorCode OpSetInvJacH1ForFlatPrism::doWork(
  int side,
  EntityType type,
  DataForcesAndSourcesCore::EntData &data
) {
  MoFEMFunctionBeginHot;
  // 

  for(int b = AINSWORTH_LEGENDRE_BASE; b!=USER_BASE; b++) {

    FieldApproximationBase base = ApproximationBaseArray[b];

    try {
      unsigned int nb_dofs = data.getN(base).size2();
      if(nb_dofs==0) MoFEMFunctionReturnHot(0);
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
      //std::cerr << data.getDiffN(base) << std::endl;
      //std::cerr << std::endl;
      switch (type) {
        case MBVERTEX:
        case MBEDGE:
        case MBTRI: {
          for(unsigned int gg = 0;gg<nb_gauss_pts;gg++) {
            for(unsigned int dd = 0;dd<nb_dofs;dd++) {
              cblas_dgemv(
                CblasRowMajor,CblasTrans,2,2,1,
                &*invJacF3.data().begin(),2,
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
      ss << "thorw in method: " << ex.what() << " at line " << __LINE__ << " in file " << __FILE__;
      SETERRQ(PETSC_COMM_SELF,MOFEM_STD_EXCEPTION_THROW,ss.str().c_str());
    }

  }

  MoFEMFunctionReturnHot(0);
}

}
