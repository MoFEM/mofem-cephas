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
  FieldInterface &m_field,const EntityType type
):
ForcesAndSurcesCore(m_field),
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
opPiolaTransform(vOlume,jAc),
opSetInvJacHdiv(invJac),
meshPositionsFieldName("MESH_NODE_POSITIONS"),
opHOatGaussPoints(hoCoordsAtGaussPts,hoGaussPtsJac,3,3),
opSetHoInvJacH1(hoGaussPtsInvJac),
opSetHoPiolaTransform(hoGaussPtsDetJac,hoGaussPtsJac),
opSetHoInvJacHdiv(hoGaussPtsInvJac),
coords(12),
jAc(3,3),
invJac(3,3),
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
};

PetscErrorCode VolumeElementForcesAndSourcesCore::operator()() {
  PetscFunctionBegin;

  try {

    if(fePtr->get_ent_type() != MBTET) PetscFunctionReturn(0);

    // Calculate volume and inverse jacobian
    {
      EntityHandle ent = fePtr->get_ent();
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
      ierr = determinantTensor2<3,double*,double>(tJac,vOlume); CHKERRQ(ierr);
      ierr = invertTensor2<3,double*,double>(tJac,vOlume,tInvJac); CHKERRQ(ierr);
      vOlume *= G_TET_W1[0]/6.;
    }

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
      ierr = getTetDataOrder(dataH1,H1); CHKERRQ(ierr);
    }

    //Hdiv
    if((dataH1.spacesOnEntities[MBTRI]).test(HDIV)) {
      ierr = getTrisSense(dataHdiv); CHKERRQ(ierr);
      ierr = getTrisDataOrder(dataHdiv,HDIV); CHKERRQ(ierr);
      ierr = getTetDataOrder(dataHdiv,HDIV); CHKERRQ(ierr);
      ierr = getFaceTriNodes(dataHdiv); CHKERRQ(ierr);
    }

    //L2
    if((dataH1.spacesOnEntities[MBTET]).test(L2)) {
      ierr = getTetDataOrder(dataL2,L2); CHKERRQ(ierr);
    }

    int order_data = getMaxDataOrder();
    int order_row = getMaxRowOrder();
    int order_col = getMaxColOrder();
    int rule = getRule(order_row,order_col,order_data);
    // std::cerr << order_data << " " << order_row << " " << order_col << " " << rule << std::endl;
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
        dataH1.dataOnEntities[MBVERTEX][0].getN(NOBASE).resize(nb_gauss_pts,4,false);
        double *shape_ptr = dataH1.dataOnEntities[MBVERTEX][0].getN(NOBASE).data().begin();
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
    } else {
      ierr = setGaussPts(order_row,order_col,order_data); CHKERRQ(ierr);
      nb_gauss_pts = gaussPts.size2();
      dataH1.dataOnEntities[MBVERTEX][0].getN(NOBASE).resize(nb_gauss_pts,4,false);
      if(nb_gauss_pts>0) {
        ierr = ShapeMBTET(
          &*dataH1.dataOnEntities[MBVERTEX][0].getN(NOBASE).data().begin(),
          &gaussPts(0,0),&gaussPts(1,0),&gaussPts(2,0),nb_gauss_pts
        ); CHKERRQ(ierr);
      }
    }
    if(nb_gauss_pts == 0) PetscFunctionReturn(0);

    /// Use the some node base
    dataHdiv.dataOnEntities[MBVERTEX][0].getNSharedPtr(NOBASE) = dataH1.dataOnEntities[MBVERTEX][0].getNSharedPtr(NOBASE);
    dataHcurl.dataOnEntities[MBVERTEX][0].getNSharedPtr(NOBASE) = dataH1.dataOnEntities[MBVERTEX][0].getNSharedPtr(NOBASE);
    dataL2.dataOnEntities[MBVERTEX][0].getNSharedPtr(NOBASE) = dataH1.dataOnEntities[MBVERTEX][0].getNSharedPtr(NOBASE);

    // Get coords at Gauss points
    {
      FTensor::Index<'i',3> i;
      double *shape_functions_ptr = &*dataH1.dataOnEntities[MBVERTEX][0].getN(NOBASE).data().begin();
      coordsAtGaussPts.resize(nb_gauss_pts,3,false);
      coordsAtGaussPts.clear();
      FTensor::Tensor1<double*,3> t_coords_at_gauss_ptr(
        &coordsAtGaussPts(0,0),&coordsAtGaussPts(0,1),&coordsAtGaussPts(0,2),3
      );
      FTensor::Tensor0<double*> t_shape_functions(shape_functions_ptr);
      for(int gg = 0;gg<nb_gauss_pts;gg++) {
        FTensor::Tensor1<double*,3> t_coords(&coords[0],&coords[1],&coords[2],3);
        for(int bb = 0;bb<4;bb++) {
          t_coords_at_gauss_ptr(i) += t_coords(i)*t_shape_functions;
          ++t_coords;
          ++t_shape_functions;
        };
        ++t_coords_at_gauss_ptr;
      }
    }

    try {

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

    try {
      ierr = opSetInvJacH1.opRhs(dataH1); CHKERRQ(ierr);
      if((dataH1.spacesOnEntities[MBTET]).test(L2)) {
        ierr = opSetInvJacH1.opRhs(dataL2); CHKERRQ(ierr);
      }
      if((dataH1.spacesOnEntities[MBTRI]).test(HDIV)) {
        ierr = opPiolaTransform.opRhs(dataHdiv); CHKERRQ(ierr);
        ierr = opSetInvJacHdiv.opRhs(dataHdiv); CHKERRQ(ierr);
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
      BitFieldId id = field_struture->get_id();

      if((fePtr->get_BitFieldId_data()&id).none()) {
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
          MatrixDouble diffN(nb_gauss_pts,12);
          for(int gg = 0;gg<nb_gauss_pts;gg++) {
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
        BitFieldId data_id = field_struture->get_id();

        if((oit->getNumeredEntFiniteElementPtr()->get_BitFieldId_data()&data_id).none()) {
          SETERRQ2(
            PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"no data field < %s > on finite element < %s >",
            field_name.c_str(),feName.c_str()
          );
        }

        if(oit->getOpType()&types[ss] || oit->getOpType()&UserDataOperator::OPROWCOL) {

          space[ss] = field_struture->get_space();
          switch(space[ss]) {
            case NOSPACE:
            SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"unknown space");
            break;
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

          base[ss] = field_struture->get_approx_base();
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

  } catch (std::exception& ex) {
    std::ostringstream ss;
    ss << "thorw in method: " << ex.what() << " at line " << __LINE__ << " in file " << __FILE__;
    SETERRQ(PETSC_COMM_SELF,MOFEM_STD_EXCEPTION_THROW,ss.str().c_str());
  }

  PetscFunctionReturn(0);
}

PetscErrorCode VolumeElementForcesAndSourcesCore::UserDataOperator::getDivergenceMatrixOperator_Hdiv(
  int side,EntityType type,
  DataForcesAndSurcesCore::EntData &data,
  int gg,
  VectorDouble &div
) {
  PetscFunctionBegin;

  try {

    int nb_dofs = data.getFieldData().size();
    if((unsigned int)nb_dofs != data.getDiffHdivN().size2()/9) {
      std::cerr << "side " << side << " type " << type << std::endl;
      SETERRQ3(
        PETSC_COMM_SELF,
        MOFEM_DATA_INCONSISTENCY,
        "data inconsistency base = %s "
        "%d != %d/9",
        FieldSpaceNames[data.getSpace()],
        nb_dofs,data.getDiffHdivN().size2()
      );
    }

    if(nb_dofs == 0) PetscFunctionReturn(0);

    int dd = 0;
    for(;dd<nb_dofs;dd++) {
      div[dd] =
      (data.getDiffHdivN(dd,gg))(0,0)+
      (data.getDiffHdivN(dd,gg))(1,1)+
      (data.getDiffHdivN(dd,gg))(2,2);
    }

  } catch (std::exception& ex) {
    std::ostringstream ss;
    ss << "thorw in method: " << ex.what() << " at line " << __LINE__ << " in file " << __FILE__;
    SETERRQ(PETSC_COMM_SELF,1,ss.str().c_str());
  }

  PetscFunctionReturn(0);
}

}
