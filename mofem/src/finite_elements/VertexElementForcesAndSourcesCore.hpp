/** \file ForcesAndSourcesCore.hpp

  \brief Implementation of elements on entities.

  Those element are inherited by user to implement specific implementation of
  particular problem.

*/

/* MIT License
 *
 * Copyright (c) 2022
 * 
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 * 
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 * 
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */

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

  VertexElementForcesAndSourcesCore(Interface &m_field);

  /** \brief default operator for VERTEX element
    \ingroup mofem_forces_and_sources_vertex_element
    */
  struct UserDataOperator : public ForcesAndSourcesCore::UserDataOperator {

    using ForcesAndSourcesCore::UserDataOperator::UserDataOperator;


    inline VectorDouble3 &getCoords();

 protected:
    MoFEMErrorCode setPtrFE(ForcesAndSourcesCore *ptr);

  };

  MoFEMErrorCode operator()();

protected:
  VectorDouble3 coords;
  friend class UserDataOperator;
};

VectorDouble3 &
VertexElementForcesAndSourcesCore::UserDataOperator::getCoords() {
  return static_cast<VertexElementForcesAndSourcesCore *>(ptrFE)->coords;
}

} // namespace MoFEM

#endif //__VERTEXELEMENTFORCESANDSOURCESCORE_HPP__

/**
 * \defgroup mofem_forces_and_sources_vertex_element Vertex Element
 * \brief Finite element and operators for vertex entity
 *
 * \ingroup mofem_forces_and_sources
 **/
