/** \file ForcesAndSourcesCore.hpp

  \brief Implementation of elements on entities.

  Those element are inherited by user to implement specific implementation of
  particular problem.

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
