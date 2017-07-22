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
  DataForcesAndSurcesCore dataHcurl;
  DerivedDataForcesAndSurcesCore derivedDataHcurl;
  DataForcesAndSurcesCore dataNoField,dataNoFieldCol;
  std::string meshPositionsFieldName;

  MatrixDouble tAngent_at_GaussPt;
  OpGetHoTangentOnEdge opGetHoTangentOnEdge;
  OpSetCovariantPiolaTransoformOnEdge opCovariantTransoform;

  EdgeElementForcesAndSurcesCore(Interface &m_field):
  ForcesAndSurcesCore(m_field),
  dataH1(MBEDGE),
  derivedDataH1(dataH1),
  dataHcurl(MBEDGE),
  derivedDataHcurl(dataHcurl),
  dataNoField(MBEDGE),
  dataNoFieldCol(MBEDGE),
  meshPositionsFieldName("MESH_NODE_POSITIONS"),
  opGetHoTangentOnEdge(tAngent_at_GaussPt),
  opCovariantTransoform(dIrection,tAngent_at_GaussPt) {
  }

  
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

    /**
     * \brief get edge length
     */
    inline double getLength() {
      return static_cast<EdgeElementForcesAndSurcesCore*>(ptrFE)->lEngth;
    }

    /**
     * \brief get measure of element
     * @return length of face
     */
    inline double getMeasure() {
      return getLength();
    }

    /**
     * \brief get edge dIrection
     */
    inline VectorDouble& getDirection() {
      return static_cast<EdgeElementForcesAndSurcesCore*>(ptrFE)->dIrection;
    }

    /**
     * \brief get edge node coordinates
     */
    inline VectorDouble& getCoords() {
      return static_cast<EdgeElementForcesAndSurcesCore*>(ptrFE)->cOords;
    }

    /**
     * \brief get matrix of integration (Gauss) points on the edge
     *  where column 0 is a coordinate X and column 1 is a value of weight
     * for example getGaussPts()(0,13) returns 0 coordinate of 13th Gauss point on particular edge element
     */
    inline MatrixDouble& getGaussPts() {
      return static_cast<EdgeElementForcesAndSurcesCore*>(ptrFE)->gaussPts;
    }

    inline FTensor::Tensor0<double*> getFTensor0IntegrationWeight() {
      return FTensor::Tensor0<double*>(&(getGaussPts()(1,0)),1);
    }

    /**
     * \brief get coordinate at integration point
     */
    inline MatrixDouble& getCoordsAtGaussPts() {
      return static_cast<EdgeElementForcesAndSurcesCore*>(ptrFE)->coordsAtGaussPts;
    }

    /**
     * \brief get tangent vector to edge curve at integration points
     */
    inline MatrixDouble& getTangetAtGaussPts() {
      return static_cast<EdgeElementForcesAndSurcesCore*>(ptrFE)->tAngent_at_GaussPt;
    }

    /**
     * \brief get pointer to this finite element
     */
    inline const EdgeElementForcesAndSurcesCore* getEdgeFE() {
      return static_cast<EdgeElementForcesAndSurcesCore*>(ptrFE);
    }

    inline FTensor::Tensor1<double*,3> getTensor1Direction() {
      double *ptr = &*getDirection().data().begin();
      return FTensor::Tensor1<double*,3>(ptr,&ptr[1],&ptr[2]);
    }

    /**
     * \brief get get coords at gauss points

     \code
     FTensor::Index<'i',3> i;
     FTensor::Tensor1<double,3> t_center;
     FTensor::Tensor1<double*,3> t_coords = getTensor1Coords();
     t_center(i) = 0;
     for(int nn = 0;nn!=2;nn++) {
        t_center(i) += t_coords(i);
        ++t_coords;
      }
      t_center(i) /= 2;
    \endcode

     */
    inline FTensor::Tensor1<double*,3> getTensor1Coords() {
      double *ptr = &*getCoords().data().begin();
      return FTensor::Tensor1<double*,3>(ptr,&ptr[1],&ptr[2],3);
    }

    inline FTensor::Tensor1<double*,3> getTensor1TangentAtGaussPts() {
      double *ptr = &*getTangetAtGaussPts().data().begin();
      return FTensor::Tensor1<double*,3>(ptr,&ptr[1],&ptr[2],3);
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
