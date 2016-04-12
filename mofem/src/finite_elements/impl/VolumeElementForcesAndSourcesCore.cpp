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
#include <FieldInterface.hpp>
#include <MeshRefinment.hpp>
#include <PrismInterface.hpp>
#include <SeriesRecorder.hpp>
#include <Core.hpp>

#include <DataStructures.hpp>
#include <DataOperators.hpp>
#include <ElementsOnEntities.hpp>
#include <VolumeElementForcesAndSourcesCore.hpp>

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

PetscErrorCode VolumeElementForcesAndSourcesCore::operator()() {
  PetscFunctionBegin;

  try {

    if(fePtr->get_ent_type() != MBTET) PetscFunctionReturn(0);

    ierr = getSpacesAndBaseOnEntities(dataH1); CHKERRQ(ierr);
    ierr = getFaceTriNodes(dataH1); CHKERRQ(ierr);
    //H1
    if((dataH1.spacesOnEntities[MBEDGE]).test(H1)) {
      ierr = getEdgesSense(dataH1); CHKERRQ(ierr);
    }
    if((dataH1.spacesOnEntities[MBTRI]).test(H1)) {
      ierr = getTrisSense(dataH1); CHKERRQ(ierr);
    }
    //Hdiv
    if((dataH1.spacesOnEntities[MBTRI]).test(HDIV)) {
      ierr = getTrisSense(dataHdiv); CHKERRQ(ierr);
      ierr = getFaceTriNodes(dataHdiv); CHKERRQ(ierr);
    }

    //H1
    if((dataH1.spacesOnEntities[MBEDGE]).test(H1)) {
      ierr = getEdgesDataOrder(dataH1,H1); CHKERRQ(ierr);
    }
    if((dataH1.spacesOnEntities[MBTRI]).test(H1)) {
      ierr = getTrisDataOrder(dataH1,H1); CHKERRQ(ierr);
    }
    if((dataH1.spacesOnEntities[MBTET]).test(H1)) {
      ierr = getTetsDataOrder(dataH1,H1); CHKERRQ(ierr);
    }

    //Hdiv
    if((dataH1.spacesOnEntities[MBTRI]).test(HDIV)) {
      ierr = getTrisSense(dataHdiv); CHKERRQ(ierr);
      ierr = getTrisDataOrder(dataHdiv,HDIV); CHKERRQ(ierr);
      ierr = getTetsDataOrder(dataHdiv,HDIV); CHKERRQ(ierr);
      ierr = getFaceTriNodes(dataHdiv); CHKERRQ(ierr);
    }

    //L2
    if((dataH1.spacesOnEntities[MBTET]).test(L2)) {
      ierr = getTetsDataOrder(dataL2,L2); CHKERRQ(ierr);
    }

    // int order = 1;
    // for(unsigned int ee = 0;ee<dataH1.dataOnEntities[MBEDGE].size();ee++) {
    //   order = max(order,dataH1.dataOnEntities[MBEDGE][ee].getDataOrder());
    // }
    // for(unsigned int ff = 0;ff<dataH1.dataOnEntities[MBTRI].size();ff++) {
    //   order = max(order,dataH1.dataOnEntities[MBTRI][ff].getDataOrder());
    // }
    // order = max(order,dataH1.dataOnEntities[MBTET][0].getDataOrder());
    // for(unsigned int ff = 0;ff<dataHdiv.dataOnEntities[MBTRI].size();ff++) {
    //   order = max(order,dataHdiv.dataOnEntities[MBTRI][ff].getDataOrder());
    // }
    // order = max(order,dataL2.dataOnEntities[MBTET][0].getDataOrder());

    int order_data = getMaxDataOrder();
    int order_row = getMaxRowOrder();
    int order_col = getMaxColOrder();
    int rule = getRule(order_row,order_col,order_data);
    // cerr << order_data << " " << order_row << " " << order_col << " " << rule << endl;
    int nb_gauss_pts;
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
        nb_gauss_pts = QUAD_3D_TABLE[rule]->npoints;
        gaussPts.resize(4,nb_gauss_pts,false);
        cblas_dcopy(
          nb_gauss_pts,&QUAD_3D_TABLE[rule]->points[1],4,&gaussPts(0,0),1
        );
        cblas_dcopy(
          nb_gauss_pts,&QUAD_3D_TABLE[rule]->points[2],4,&gaussPts(1,0),1
        );
        cblas_dcopy(
          nb_gauss_pts,&QUAD_3D_TABLE[rule]->points[3],4,&gaussPts(2,0),1
        );
        cblas_dcopy(
          nb_gauss_pts,QUAD_3D_TABLE[rule]->weights,1,&gaussPts(3,0),1
        );
        dataH1.dataOnEntities[MBVERTEX][0].getN().resize(nb_gauss_pts,4,false);
        double *shape_ptr = dataH1.dataOnEntities[MBVERTEX][0].getN().data().begin();
        cblas_dcopy(
          4*nb_gauss_pts,QUAD_3D_TABLE[rule]->points,1,shape_ptr,1
        );
      } else {
        SETERRQ2(
          PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"rule > quadrature order %d < %d",
          rule,QUAD_3D_TABLE_SIZE
        );
        nb_gauss_pts = 0;
      }
      // nb_gauss_pts = gm_rule_size(rule,3);
      // gaussPts.resize(4,nb_gauss_pts,false);
      // ierr = Grundmann_Moeller_integration_points_3D_TET(
      //   rule,&gaussPts(0,0),&gaussPts(1,0),&gaussPts(2,0),&gaussPts(3,0)
      // ); CHKERRQ(ierr);
    } else {
      ierr = setGaussPts(order_row,order_col,order_data); CHKERRQ(ierr);
      nb_gauss_pts = gaussPts.size2();
      dataH1.dataOnEntities[MBVERTEX][0].getN().resize(nb_gauss_pts,4,false);
      if(nb_gauss_pts>0) {
        ierr = ShapeMBTET(
          &*dataH1.dataOnEntities[MBVERTEX][0].getN().data().begin(),
          &gaussPts(0,0),&gaussPts(1,0),&gaussPts(2,0),nb_gauss_pts
        ); CHKERRQ(ierr);
      }
    }
    if(nb_gauss_pts == 0) PetscFunctionReturn(0);

    cerr << dataH1.basesOnEntities[MBVERTEX] << endl;
    cerr << dataH1.basesOnEntities[MBVERTEX].test(AINSWORTH_COLE_BASE) << endl;
    cerr << dataH1.basesOnEntities[MBEDGE] << endl;
    cerr << dataH1.basesOnEntities[MBEDGE].test(AINSWORTH_COLE_BASE) << endl;

    if(dataH1.bAse.test(AINSWORTH_COLE_BASE)) {
      ierr = shapeTETFunctions_H1(
        dataH1,&gaussPts(0,0),&gaussPts(1,0),&gaussPts(2,0),nb_gauss_pts,Legendre_polynomials
      ); CHKERRQ(ierr);
      if((dataH1.spacesOnEntities[MBTRI]).test(HDIV)) {
        ierr = shapeTETFunctions_Hdiv(
          dataHdiv,&gaussPts(0,0),&gaussPts(1,0),&gaussPts(2,0),nb_gauss_pts,Legendre_polynomials
        ); CHKERRQ(ierr);
      }
      if((dataH1.spacesOnEntities[MBTET]).test(L2)) {
        ierr = shapeTETFunctions_L2(
          dataL2,&gaussPts(0,0),&gaussPts(1,0),&gaussPts(2,0),nb_gauss_pts,Legendre_polynomials
        ); CHKERRQ(ierr);
      }
    }
    if(dataH1.basesOnEntities[MBVERTEX].test(LOBATTO_BASE)) {
    }


    EntityHandle ent = fePtr->get_ent();
    rval = mField.get_moab().get_connectivity(ent,conn,num_nodes,true); CHKERRQ_MOAB(rval);
    coords.resize(num_nodes*3,false);
    rval = mField.get_moab().get_coords(conn,num_nodes,&*coords.data().begin()); CHKERRQ_MOAB(rval);
    vOlume = ShapeVolumeMBTET(&*dataH1.dataOnEntities[MBVERTEX][0].getDiffN().data().begin(),&*coords.data().begin());
    Jac.resize(3,3,false);
    invJac.resize(3,3,false);
    ierr = ShapeJacMBTET(
      &*dataH1.dataOnEntities[MBVERTEX][0].getDiffN().data().begin(),&*coords.begin(),&*Jac.data().begin()
    ); CHKERRQ(ierr);
    noalias(invJac) = Jac;
    ierr = ShapeInvJacVolume(&*invJac.data().begin()); CHKERRQ(ierr);

    coordsAtGaussPts.resize(nb_gauss_pts,3,false);
    for(int gg = 0;gg<nb_gauss_pts;gg++) {
      for(int dd = 0;dd<3;dd++) {
        coordsAtGaussPts(gg,dd) = cblas_ddot(4,&dataH1.dataOnEntities[MBVERTEX][0].getN()(gg,0),1,&coords[dd],3);
      }
    }

    try {
      ierr = opSetInvJacH1.opRhs(dataH1); CHKERRQ(ierr);
      if((dataH1.spacesOnEntities[MBTET]).test(L2)) {
        ierr = opSetInvJacH1.opRhs(dataL2); CHKERRQ(ierr);
      }
      if((dataH1.spacesOnEntities[MBTRI]).test(HDIV)) {
        ierr = opPiolaTransform.opRhs(dataHdiv); CHKERRQ(ierr);
        ierr = opSetInvJacHdiv.opRhs(dataHdiv); CHKERRQ(ierr);
      }
    } catch (exception& ex) {
      ostringstream ss;
      ss << "thorw in method: " << ex.what() << " at line " << __LINE__ << " in file " << __FILE__;
      SETERRQ(PETSC_COMM_SELF,MOFEM_STD_EXCEPTION_THROW,ss.str().c_str());
    }

    if(
      dataPtr->get<FieldName_mi_tag>().find(meshPositionsFieldName)!=
      dataPtr->get<FieldName_mi_tag>().end()
    ) {
      BitFieldId id = mField.get_field_structure(meshPositionsFieldName)->get_id();
      if((fePtr->get_BitFieldId_data()&id).none()) {
        SETERRQ(PETSC_COMM_SELF,MOFEM_NOT_FOUND,"no MESH_NODE_POSITIONS in element data");
      }
      ierr = getEdgesDataOrder(dataH1,meshPositionsFieldName); CHKERRQ(ierr);
      ierr = getTrisDataOrder(dataH1,meshPositionsFieldName); CHKERRQ(ierr);
      ierr = getTetsDataOrder(dataH1,meshPositionsFieldName); CHKERRQ(ierr);
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
        ublas::noalias(hoGaussPtsInvJac) = hoGaussPtsJac;
        MatrixDouble jac(3,3);
        hoGaussPtsDetJac.resize(nb_gauss_pts,false);
        for(int gg = 0;gg<nb_gauss_pts;gg++) {
          cblas_dcopy(9,&hoGaussPtsJac(gg,0),1,&jac(0,0),1);
          hoGaussPtsDetJac[gg] = ShapeDetJacVolume(&jac(0,0));
          if(hoGaussPtsDetJac[gg]<0) {
            SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"Negative volume");
          }
          ierr = ShapeInvJacVolume(&hoGaussPtsInvJac(gg,0)); CHKERRQ(ierr);
        }
        ierr = opSetHoInvJacH1.opRhs(dataH1); CHKERRQ(ierr);
        if((dataH1.spacesOnEntities[MBTET]).test(L2)) {
          ierr = opSetHoInvJacH1.opRhs(dataL2); CHKERRQ(ierr);
        }
        if((dataH1.spacesOnEntities[MBTRI]).test(HDIV)) {
          ierr = opSetHoPiolaTransform.opRhs(dataHdiv); CHKERRQ(ierr);
          ierr = opSetHoInvJacHdiv.opRhs(dataHdiv); CHKERRQ(ierr);
        }
      } catch (exception& ex) {
        ostringstream ss;
        ss << "problem with indices in method: " << ex.what() << " at line " << __LINE__ << " in file " << __FILE__;
        SETERRQ(PETSC_COMM_SELF,MOFEM_STD_EXCEPTION_THROW,ss.str().c_str());
      }
    } else {

      hoCoordsAtGaussPts.resize(0,0,false);
      hoGaussPtsInvJac.resize(0,0,false);
      hoGaussPtsDetJac.resize(0,false);

      MatrixDouble diffN(nb_gauss_pts,12);
      for(int gg = 0;gg<nb_gauss_pts;gg++) {
        for(int nn = 0;nn<4;nn++) {
          for(int dd = 0;dd<3;dd++) {
            diffN(gg,nn*3+dd) = dataH1.dataOnEntities[MBVERTEX][0].getDiffN()(nn,dd);
          }
        }
      }
      dataH1.dataOnEntities[MBVERTEX][0].getDiffN().resize(diffN.size1(),diffN.size2(),false);
      dataH1.dataOnEntities[MBVERTEX][0].getDiffN().data().swap(diffN.data());

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

          if(last_eval_field_name[ss]!=field_name) {

            switch(space[ss]) {
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
              ierr = getEdgesDataOrder(*op_data[ss],field_name); CHKERRQ(ierr);
              ierr = getEdgesFieldData(*op_data[ss],field_name); CHKERRQ(ierr);
              // ierr = getEdgesFieldDofs(*op_data[ss],field_name); CHKERRQ(ierr);
              case HDIV:
              if(!ss) {
                ierr = getTrisRowIndices(*op_data[ss],field_name); CHKERRQ(ierr);
              } else {
                ierr = getTrisColIndices(*op_data[ss],field_name); CHKERRQ(ierr);
              }
              ierr = getTrisDataOrder(*op_data[ss],field_name); CHKERRQ(ierr);
              ierr = getTrisFieldData(*op_data[ss],field_name); CHKERRQ(ierr);
              case L2:
              if(!ss) {
                ierr = getTetsRowIndices(*op_data[ss],field_name); CHKERRQ(ierr);
              } else {

                ierr = getTetsColIndices(*op_data[ss],field_name); CHKERRQ(ierr);
              }
              ierr = getTetsDataOrder(*op_data[ss],field_name); CHKERRQ(ierr);
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
            oit->doQuadsCol,
            oit->doTrisCol,
            oit->doTetsCol,
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

  } catch (exception& ex) {
    ostringstream ss;
    ss << "thorw in method: " << ex.what() << " at line " << __LINE__ << " in file " << __FILE__;
    SETERRQ(PETSC_COMM_SELF,MOFEM_STD_EXCEPTION_THROW,ss.str().c_str());
  }

  PetscFunctionReturn(0);
}

PetscErrorCode VolumeElementForcesAndSourcesCore::UserDataOperator::getDivergenceMatrixOperato_Hdiv(
  int side,EntityType type,DataForcesAndSurcesCore::EntData &data,
  int gg,VectorDouble &div) {
    PetscFunctionBegin;

    try {

      int nb_dofs = data.getFieldData().size();
      if((unsigned int)nb_dofs != data.getDiffHdivN().size2()/9) {
        SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"data inconsistency");
      }

      if(nb_dofs == 0) PetscFunctionReturn(0);

      int dd = 0;
      for(;dd<nb_dofs;dd++) {
        div[dd] =
        (data.getDiffHdivN(dd,gg))(0,0)+
        (data.getDiffHdivN(dd,gg))(1,1)+
        (data.getDiffHdivN(dd,gg))(2,2);
      }

    } catch (exception& ex) {
      ostringstream ss;
      ss << "thorw in method: " << ex.what() << " at line " << __LINE__ << " in file " << __FILE__;
      SETERRQ(PETSC_COMM_SELF,1,ss.str().c_str());
    }

    PetscFunctionReturn(0);
  }

}
