/** \file ForcesAndSourcesCore.hpp

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
struct VertexElementForcesAndSourcesCore : public ForcesAndSourcesCore {

  std::string meshPositionsFieldName;

  VertexElementForcesAndSourcesCore(Interface &m_field)
      : ForcesAndSourcesCore(m_field){};

  VectorDouble coords;

  /** \brief default operator for VERTEX element
    \ingroup mofem_forces_and_sources_vertex_element
    */
  struct UserDataOperator : public ForcesAndSourcesCore::UserDataOperator {

    using ForcesAndSourcesCore::UserDataOperator::UserDataOperator;

    inline VectorDouble &getCoords() {
      return static_cast<VertexElementForcesAndSourcesCore *>(ptrFE)->coords;
    }
  };

  MoFEMErrorCode operator()();
};

} // namespace MoFEM

#endif //__VERTEXELEMENTFORCESANDSOURCESCORE_HPP__

/**
 * \defgroup mofem_forces_and_sources_vertex_element Vertex Element
 * \brief Finite element and operators for vertex entity
 *
 * \ingroup mofem_forces_and_sources
 **/
