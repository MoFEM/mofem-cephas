/** \file ForcesAndSurcesCore.hpp

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


#ifndef __VERTEXELEMENTFORCESANDSOURCESCORE_HPP__
#define __VERTEXELEMENTFORCESANDSOURCESCORE_HPP__

using namespace boost::numeric;

namespace MoFEM {

/** \brief Vertex finite element
 * \ingroup mofem_forces_and_sources_vertex_element

 User is implementing own operator at Gauss points level, by own object
 derived from VertexElementForcesAndSourcesCoreL::UserDataOperator.  Arbitrary
 number of operator added pushing objects to rowOpPtrVector and
 rowColOpPtrVector.

 */
struct VertexElementForcesAndSourcesCore: public ForcesAndSurcesCore {

  DataForcesAndSurcesCore data;
  DerivedDataForcesAndSurcesCore derivedData;
  DataForcesAndSurcesCore dataNoField,dataNoFieldCol;
  std::string meshPositionsFieldName;

  VertexElementForcesAndSourcesCore(Interface &m_field):
    ForcesAndSurcesCore(m_field),
    data(MBVERTEX),
    derivedData(data),
    dataNoField(MBVERTEX),
    dataNoFieldCol(MBVERTEX)
  {};

  MoABErrorCode rval;
  VectorDouble coords;

  /** \brief default operator for VERTEX element
    \ingroup mofem_forces_and_sources_vertex_element
    */
  struct UserDataOperator: public ForcesAndSurcesCore::UserDataOperator {

    UserDataOperator(
      const std::string &field_name,const char type):
      ForcesAndSurcesCore::UserDataOperator(field_name,type) {}

    UserDataOperator(
      const std::string &row_field_name,const std::string &col_field_name,const char type):
      ForcesAndSurcesCore::UserDataOperator(row_field_name,col_field_name,type) {}

    inline VectorDouble& getCoords() {
      return static_cast<VertexElementForcesAndSourcesCore*>(ptrFE)->coords;
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

#endif //__VERTEXELEMENTFORCESANDSOURCESCORE_HPP__

/***************************************************************************//**
 * \defgroup mofem_forces_and_sources_vertex_element Vertex Element
 * \brief Finite element and operators for vertex entity
 *
 * \ingroup mofem_forces_and_sources
 ******************************************************************************/
