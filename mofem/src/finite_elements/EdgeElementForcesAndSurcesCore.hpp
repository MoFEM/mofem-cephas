/** \file ElementsOnEntities.hpp

  \brief Implementation of elements on entities.

  Those element are inherited by user to implement specific implementation of
  particular problem.

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

#ifndef __EDGEELEMENTFORCESANDSURCESCORE_HPP__
#define __EDGEELEMENTFORCESANDSURCESCORE_HPP__

using namespace boost::numeric;

namespace MoFEM {

/** \brief Edge finite element
 * \ingroup mofem_forces_and_sources_edge_element
 *
 * User is implementing own operator at Gauss points level, by own object
 * derived from EdgeElementForcesAndSurcesCoreL::UserDataOperator.  Arbitrary
 * number of operator added pushing objects to rowOpPtrVector and
 * rowColOpPtrVector.
 *
 */
struct EdgeElementForcesAndSurcesCore: public ForcesAndSurcesCore {

  DataForcesAndSurcesCore dataH1;
  DerivedDataForcesAndSurcesCore derivedDataH1;
  DataForcesAndSurcesCore dataNoField,dataNoFieldCol;
  std::string meshPositionsFieldName;

  MatrixDouble tAngent_at_GaussPt;
  OpGetHoTangentOnEdge opGetHoTangentOnEdge;

  EdgeElementForcesAndSurcesCore(FieldInterface &m_field):
    ForcesAndSurcesCore(m_field),
    dataH1(MBEDGE),
    derivedDataH1(dataH1),
    dataNoField(MBEDGE),
    dataNoFieldCol(MBEDGE),
    meshPositionsFieldName("MESH_NODE_POSITIONS"),
    opGetHoTangentOnEdge(tAngent_at_GaussPt)
  {};

  MoABErrorCode rval;
  double lEngth;;
  VectorDouble dIrection;
  VectorDouble cOords;
  MatrixDouble gaussPts;
  MatrixDouble coordsAtGaussPts;

  /** \brief default operator for EDGE element
    \ingroup mofem_forces_and_sources_edge_element
    */
  struct UserDataOperator: public ForcesAndSurcesCore::UserDataOperator {

    UserDataOperator(
      const std::string &field_name,const char type):
      ForcesAndSurcesCore::UserDataOperator(field_name,type) {}

    UserDataOperator(
      const std::string &row_field_name,const std::string &col_field_name,const char type):
      ForcesAndSurcesCore::UserDataOperator(row_field_name,col_field_name,type) {}

    inline double getLength() {
      return static_cast<EdgeElementForcesAndSurcesCore*>(ptrFE)->lEngth;
    }
    inline VectorDouble& getDirection() {
      return static_cast<EdgeElementForcesAndSurcesCore*>(ptrFE)->dIrection;
    }
    inline VectorDouble& getCoords() {
      return static_cast<EdgeElementForcesAndSurcesCore*>(ptrFE)->cOords;
    }
    inline MatrixDouble& getGaussPts() {
      return static_cast<EdgeElementForcesAndSurcesCore*>(ptrFE)->gaussPts;
    }
    inline MatrixDouble& getCoordsAtGaussPts() {
      return static_cast<EdgeElementForcesAndSurcesCore*>(ptrFE)->coordsAtGaussPts;
    }
    inline MatrixDouble& getTangetAtGaussPtrs() {
      return static_cast<EdgeElementForcesAndSurcesCore*>(ptrFE)->tAngent_at_GaussPt;
    }
    inline const EdgeElementForcesAndSurcesCore* getEdgeFE() {
      return static_cast<EdgeElementForcesAndSurcesCore*>(ptrFE);
    }

  };

  PetscErrorCode preProcess() {
    PetscFunctionBegin;
    PetscFunctionReturn(0);
  }
  PetscErrorCode operator()();
  PetscErrorCode postProcess() {
    PetscFunctionBegin;
    PetscFunctionReturn(0);
  }

};

}

#endif //__EDGEELEMENTFORCESANDSURCESCORE_HPP__

/***************************************************************************//**
 * \defgroup mofem_forces_and_sources_edge_element Edge Element
 *
 * \brief Implementation of general edge element.
 *
 * \ingroup mofem_forces_and_sources
 ******************************************************************************/
