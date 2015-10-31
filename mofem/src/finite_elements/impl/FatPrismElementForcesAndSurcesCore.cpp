/** \file FatPrismElementForcesAndSurcesCore.cpp

\brief Implementation of fat prism element

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

#include <petscsys.h>
#include <petscvec.h>
#include <petscmat.h>
#include <petscsnes.h>
#include <petscts.h>

#include <definitions.h>
#include <h1_hdiv_hcurl_l2.h>
#include <fem_tools.h>

#include <Common.hpp>

#include <LoopMethods.hpp>
#include <FieldInterface.hpp>

#define BOOST_UBLAS_SHALLOW_ARRAY_ADAPTOR

#include <boost/numeric/ublas/storage.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/io.hpp>

#include <boost/ptr_container/ptr_vector.hpp>
#include <boost/ptr_container/ptr_map.hpp>

#include <DataStructures.hpp>
#include <DataOperators.hpp>
#include <ElementsOnEntities.hpp>
#include <FatPrismElementForcesAndSurcesCore.hpp>

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

PetscErrorCode FatPrismElementForcesAndSurcesCore::operator()() {
  PetscFunctionBegin;


  if(fePtr->get_ent_type() != MBPRISM) PetscFunctionReturn(0);

  try {

    ierr = getSpacesOnEntities(dataH1); CHKERRQ(ierr);
    ierr = getSpacesOnEntities(dataH1TrianglesOnly); CHKERRQ(ierr);

    //H1
    if((dataH1.spacesOnEntities[MBEDGE]).test(H1)) {
      ierr = getEdgesSense(dataH1); CHKERRQ(ierr);
      ierr = getTrisSense(dataH1); CHKERRQ(ierr);
      ierr = getQuadSense(dataH1); CHKERRQ(ierr);
      ierr = getEdgesOrder(dataH1,H1); CHKERRQ(ierr);
      ierr = getTrisOrder(dataH1,H1); CHKERRQ(ierr);
      ierr = getQuadOrder(dataH1,H1); CHKERRQ(ierr);
      ierr = getPrismOrder(dataH1,H1); CHKERRQ(ierr);
      //Triangles only
      ierr = getEdgesSense(dataH1TrianglesOnly); CHKERRQ(ierr);
      ierr = getTrisSense(dataH1TrianglesOnly); CHKERRQ(ierr);
      ierr = getEdgesOrder(dataH1TrianglesOnly,H1); CHKERRQ(ierr);
      ierr = getTrisOrder(dataH1TrianglesOnly,H1); CHKERRQ(ierr);
    }

    //Hdiv
    if((dataH1.spacesOnEntities[MBTRI]).test(HDIV)) {
      SETERRQ(PETSC_COMM_SELF,MOFEM_NOT_IMPLEMENTED,"not implemented yet");
    }

    //Hcurl
    if((dataH1.spacesOnEntities[MBEDGE]).test(HCURL)) {
      SETERRQ(PETSC_COMM_SELF,MOFEM_NOT_IMPLEMENTED,"not implemented yet");
    }

    //L2
    if((dataH1.spacesOnEntities[MBPRISM]).test(L2)) {
      SETERRQ(PETSC_COMM_SELF,MOFEM_NOT_IMPLEMENTED,"not implemented yet");
    }

    int order_triangles_only = 1;
    for(unsigned int ee = 0;ee<dataH1TrianglesOnly.dataOnEntities[MBEDGE].size();ee++) {
      order_triangles_only = max(
        order_triangles_only,dataH1TrianglesOnly.dataOnEntities[MBEDGE][ee].getOrder()
      );
    }
    for(unsigned int ff = 0;ff<dataH1TrianglesOnly.dataOnEntities[MBTRI].size();ff++) {
      order_triangles_only = max(
        order_triangles_only,dataH1TrianglesOnly.dataOnEntities[MBTRI][ff].getOrder()
      );
    }

    int nb_gauss_pts_on_face;
    int rule = getRule(order_triangles_only);
    if(rule >= 0) {
      nb_gauss_pts_on_face = gm_rule_size(rule,2);
      gaussPtsTrianglesOnly.resize(3,nb_gauss_pts_on_face,false);
      ierr = Grundmann_Moeller_integration_points_2D_TRI(
        rule,
        &gaussPtsTrianglesOnly(0,0),
        &gaussPtsTrianglesOnly(1,0),
        &gaussPtsTrianglesOnly(2,0)
      ); CHKERRQ(ierr);
    } else {
      ierr = setGaussPts(order_triangles_only); CHKERRQ(ierr);
      nb_gauss_pts_on_face = gaussPtsTrianglesOnly.size2();
    }
    if(nb_gauss_pts_on_face == 0) PetscFunctionReturn(0);

    ierr = shapeFlatPRISMFunctions_H1(
      dataH1TrianglesOnly,
      &gaussPtsTrianglesOnly(0,0),
      &gaussPtsTrianglesOnly(1,0),
      nb_gauss_pts_on_face
    ); CHKERRQ(ierr);

    EntityHandle ent = fePtr->get_ent();
    int num_nodes;
    const EntityHandle* conn;
    rval = mField.get_moab().get_connectivity(ent,conn,num_nodes,true); CHKERR_PETSC(rval);
    coords.resize(num_nodes*3,false);
    rval = mField.get_moab().get_coords(conn,num_nodes,&*coords.data().begin()); CHKERR_PETSC(rval);

    normal.resize(6,false);
    ierr = ShapeFaceNormalMBTRI(
      &*dataH1TrianglesOnly.dataOnEntities[MBVERTEX][0].getDiffN().data().begin(),
      &coords[0],
      &normal[0]
    ); CHKERRQ(ierr);
    ierr = ShapeFaceNormalMBTRI(
      &*dataH1TrianglesOnly.dataOnEntities[MBVERTEX][0].getDiffN().data().begin(),
      &coords[9],
      &normal[3]
    ); CHKERRQ(ierr);
    aRea[0] = cblas_dnrm2(3,&normal[0],1)*0.5;
    aRea[1] = cblas_dnrm2(3,&normal[3],1)*0.5;

    coordsAtGaussPtsTrianglesOnly.resize(nb_gauss_pts_on_face,6,false);
    for(int gg = 0;gg<nb_gauss_pts_on_face;gg++) {
      for(int dd = 0;dd<3;dd++) {
        coordsAtGaussPtsTrianglesOnly(gg,dd) = cblas_ddot(
          3,&dataH1TrianglesOnly.dataOnEntities[MBVERTEX][0].getN()(gg,0),1,&coords[dd],3
        );
        coordsAtGaussPtsTrianglesOnly(gg,3+dd) = cblas_ddot(
          3,&dataH1TrianglesOnly.dataOnEntities[MBVERTEX][0].getN()(gg,0),1,&coords[9+dd],3
        );
      }
    }

    try {

      if(
        dataPtr->get<FieldName_mi_tag>().find(meshPositionsFieldName)!=
        dataPtr->get<FieldName_mi_tag>().end()
      ) {
        hoCoordsAtGaussPtsF3.resize(nb_gauss_pts_on_face,3,false);
        nOrmals_at_GaussPtF3.resize(nb_gauss_pts_on_face,3,false);
        tAngent1_at_GaussPtF3.resize(nb_gauss_pts_on_face,3,false);
        tAngent2_at_GaussPtF3.resize(nb_gauss_pts_on_face,3,false);
        hoCoordsAtGaussPtsF4.resize(nb_gauss_pts_on_face,3,false);
        nOrmals_at_GaussPtF4.resize(nb_gauss_pts_on_face,3,false);
        tAngent1_at_GaussPtF4.resize(nb_gauss_pts_on_face,3,false);
        tAngent2_at_GaussPtF4.resize(nb_gauss_pts_on_face,3,false);
        ierr = getEdgesOrder(dataH1TrianglesOnly,meshPositionsFieldName); CHKERRQ(ierr);
        ierr = getTrisOrder(dataH1TrianglesOnly,meshPositionsFieldName); CHKERRQ(ierr);
        ierr = getNodesFieldData(dataH1TrianglesOnly,meshPositionsFieldName); CHKERRQ(ierr);
        ierr = getEdgesFieldData(dataH1TrianglesOnly,meshPositionsFieldName); CHKERRQ(ierr);
        ierr = getTrisFieldData(dataH1TrianglesOnly,meshPositionsFieldName); CHKERRQ(ierr);
        ierr = getNodesFieldDofs(dataH1TrianglesOnly,meshPositionsFieldName); CHKERRQ(ierr);
        ierr = getEdgesFieldDofs(dataH1TrianglesOnly,meshPositionsFieldName); CHKERRQ(ierr);
        ierr = getTrisFieldDofs(dataH1TrianglesOnly,meshPositionsFieldName); CHKERRQ(ierr);
        try {
          ierr = opHOCoordsAndNormals.opRhs(dataH1TrianglesOnly); CHKERRQ(ierr);
          ierr = opHOCoordsAndNormals.calculateNormals(); CHKERRQ(ierr);
        } catch (exception& ex) {
          ostringstream ss;
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
    } catch (exception& ex) {
      ostringstream ss;
      ss << "thorw in method: " << ex.what() << " at line " << __LINE__ << " in file " << __FILE__;
      SETERRQ(PETSC_COMM_SELF,MOFEM_STD_EXCEPTION_THROW,ss.str().c_str());
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
            SETERRQ(PETSC_COMM_SELF,MOFEM_NOT_IMPLEMENTED,"not implemented yet");
            break;
            case L2:
            SETERRQ(PETSC_COMM_SELF,MOFEM_NOT_IMPLEMENTED,"not implemented yet");
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
              if(!ss) {
                ierr = getQuadRowIndices(*op_data[ss],field_name); CHKERRQ(ierr);
              } else {
                ierr = getQuadColIndices(*op_data[ss],field_name); CHKERRQ(ierr);
              }
              ierr = getQuadOrder(*op_data[ss],field_name); CHKERRQ(ierr);
              ierr = getQuadFieldData(*op_data[ss],field_name); CHKERRQ(ierr);
              ierr = getQuadFieldDofs(*op_data[ss],field_name); CHKERRQ(ierr);
              case L2:
              if(!ss) {
                ierr = getPrismRowIndices(*op_data[ss],field_name); CHKERRQ(ierr);
              } else {
                ierr = getPrismColIndices(*op_data[ss],field_name); CHKERRQ(ierr);
              }
              ierr = getPrismOrder(*op_data[ss],field_name); CHKERRQ(ierr);
              ierr = getPrismFieldData(*op_data[ss],field_name); CHKERRQ(ierr);
              ierr = getPrismFieldDofs(*op_data[ss],field_name); CHKERRQ(ierr);
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

  } catch (MoFEMException const &e) {
    SETERRQ(PETSC_COMM_SELF,e.errorCode,e.errorMessage);
  } catch (exception& ex) {
    ostringstream ss;
    ss << "thorw in method: " << ex.what() << " at line " << __LINE__ << " in file " << __FILE__;
    SETERRQ(PETSC_COMM_SELF,MOFEM_STD_EXCEPTION_THROW,ss.str().c_str());
  }

  PetscFunctionReturn(0);
}

}
