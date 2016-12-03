/** \file VertexElementForcesAndSourcesCore.cpp

\brief Implementation of vertex element

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
#include <MeshRefinement.hpp>
#include <PrismInterface.hpp>
#include <SeriesRecorder.hpp>
#include <Core.hpp>

#include <FTensor.hpp>
#include <DataStructures.hpp>
#include <DataOperators.hpp>
#include <ForcesAndSurcesCore.hpp>
#include <VertexElementForcesAndSourcesCore.hpp>

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

PetscErrorCode VertexElementForcesAndSourcesCore::operator()() {
  PetscFunctionBegin;

  if(numeredEntFiniteElementPtr->getEntType() != MBVERTEX) PetscFunctionReturn(0);

  EntityHandle ent = numeredEntFiniteElementPtr->getEnt();
  coords.resize(3,false);
  rval = mField.get_moab().get_coords(&ent,1,&*coords.data().begin()); CHKERRQ_MOAB(rval);

  const UserDataOperator::OpType types[2] = {
    UserDataOperator::OPROW, UserDataOperator::OPCOL
  };
  std::vector<std::string> last_eval_field_name(2);
  DataForcesAndSurcesCore *op_data[2];
  FieldSpace space[2];

  boost::ptr_vector<UserDataOperator>::iterator oit,hi_oit;
  oit = opPtrVector.begin();
  hi_oit = opPtrVector.end();

  for(;oit!=hi_oit;oit++) {

    oit->setPtrFE(this);

    for(int ss = 0;ss!=2;ss++) {

      std::string field_name = !ss ? oit->rowFieldName : oit->colFieldName;
      BitFieldId data_id = mField.get_field_structure(field_name)->getId();
      if((oit->getNumeredEntFiniteElementPtr()->getBitFieldIdData()&data_id).none()) {
        SETERRQ2(
          PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"no data field < %s > on finite element < %s >",
          field_name.c_str(),feName.c_str()
        );
      }

      if(oit->getOpType()&types[ss] || oit->getOpType()&UserDataOperator::OPROWCOL) {

        space[ss] = mField.get_field_structure(field_name)->getSpace();

        switch(space[ss]) {
          case NOSPACE:
          SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"unknown space");
          case H1:
          op_data[ss] = !ss ? &data : &derivedData;
          break;
          case HCURL:
          SETERRQ(PETSC_COMM_SELF,MOFEM_NOT_IMPLEMENTED,"not make sense for vertex");
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
            case NOSPACE:
            SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"unknown space");
            case H1:
            if(!ss) {
              ierr = getRowNodesIndices(*op_data[ss],field_name); CHKERRQ(ierr);
            } else {
              ierr = getColNodesIndices(*op_data[ss],field_name); CHKERRQ(ierr);
            }
            ierr = getNodesFieldData(*op_data[ss],field_name); CHKERRQ(ierr);
            break;
            case HCURL:
            SETERRQ(PETSC_COMM_SELF,MOFEM_NOT_IMPLEMENTED,"not make sanes on edge");
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
          false,
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
          false,
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
