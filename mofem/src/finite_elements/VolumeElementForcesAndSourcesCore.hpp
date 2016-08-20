/** \file VolumeElementForcesAndSourcesCore.hpp

  \brief Volume element.

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

#ifndef __VOLUMEELEMENTFORCESANDSOURCESCORE_HPP__
#define __VOLUMEELEMENTFORCESANDSOURCESCORE_HPP__

using namespace boost::numeric;

namespace MoFEM {

/** \brief Volume finite element
 \ingroup mofem_forces_and_sources_volume_element

 User is implementing own operator at Gauss point level, by own object
 derived from VolumeElementForcesAndSourcesCore::UserDataOperator.  Arbitrary
 number of operator added pushing objects to rowOpPtrVector and
 rowColOpPtrVector.

 */
struct VolumeElementForcesAndSourcesCore: public ForcesAndSurcesCore {

  VectorDouble coords;
  MatrixDouble3by3 jAc;
  MatrixDouble3by3 invJac;

  DataForcesAndSurcesCore dataH1;
  DerivedDataForcesAndSurcesCore derivedDataH1;
  DataForcesAndSurcesCore dataL2;
  DerivedDataForcesAndSurcesCore derivedDataL2;
  DataForcesAndSurcesCore dataHdiv;
  DerivedDataForcesAndSurcesCore derivedDataHdiv;
  DataForcesAndSurcesCore dataHcurl;
  DerivedDataForcesAndSurcesCore derivedDataHcurl;
  DataForcesAndSurcesCore dataNoField;
  DataForcesAndSurcesCore dataNoFieldCol;

  OpSetInvJacH1 opSetInvJacH1;
  OpSetContravariantPiolaTransform opContravariantPiolaTransform;
  OpSetCovariantPiolaTransform opCovariantPiolaTransform;
  OpSetInvJacHdivAndHcurl opSetInvJacHdivAndHcurl;

  std::string meshPositionsFieldName;
  MatrixDouble hoCoordsAtGaussPts;
  MatrixDouble hoGaussPtsJac;
  MatrixDouble hoGaussPtsInvJac;
  VectorDouble hoGaussPtsDetJac;

  OpGetDataAndGradient opHOatGaussPoints; ///< higher order geometry data at Gauss pts
  OpSetHoInvJacH1 opSetHoInvJacH1;
  OpSetHoContravariantPiolaTransform opHoContravariantTransform;
  OpSetHoCovariantPiolaTransform opHoCovariantTransform;
  OpSetHoInvJacHdivAndHcurl opSetHoInvJacHdivAndHcurl;

  VolumeElementForcesAndSourcesCore(Interface &m_field,const EntityType type = MBTET);
  virtual ~VolumeElementForcesAndSourcesCore() {}

  MoABErrorCode rval;
  double vOlume;

  int num_nodes;
  const EntityHandle* conn;
  FTensor::Tensor2<double*,3,3> tJac;
  FTensor::Tensor2<double*,3,3> tInvJac;

  MatrixDouble gaussPts;
  MatrixDouble coordsAtGaussPts;

  /** \brief default operator for TET element
    * \ingroup mofem_forces_and_sources_volume_element
    */
  struct UserDataOperator: public ForcesAndSurcesCore::UserDataOperator {

    UserDataOperator(
      const std::string &field_name,const char type
    ):
    ForcesAndSurcesCore::UserDataOperator(field_name,type) {
    }

    UserDataOperator(
      const std::string &row_field_name,const std::string &col_field_name,const char type
    ):
    ForcesAndSurcesCore::UserDataOperator(row_field_name,col_field_name,type) {
    }

    /** \brief get element number of nodes
    */
    inline int getNumNodes() {
      return static_cast<VolumeElementForcesAndSourcesCore*>(ptrFE)->num_nodes;
    }

    /** \brief get element connectivity
     */
    inline const EntityHandle* getConn() {
      return static_cast<VolumeElementForcesAndSourcesCore*>(ptrFE)->conn;
    }

    /** \brief element volume (linear geometry)
      */
    inline double getVolume() {
      return static_cast<VolumeElementForcesAndSourcesCore*>(ptrFE)->vOlume;
    }

    /** \brief nodal coordinates
      */
    inline VectorDouble& getCoords() { return
      static_cast<VolumeElementForcesAndSourcesCore*>(ptrFE)->coords;
    }

    /** \brief matrix of Gauss pts
      */
    inline MatrixDouble& getGaussPts() {
      return static_cast<VolumeElementForcesAndSourcesCore*>(ptrFE)->gaussPts;
    }

    /** \brief Gauss points and weight, matrix (nb. of points x 4)

      Column 0-3 and 4 represents Gauss pts coordinate and weight, respectively.

      */
    inline MatrixDouble& getCoordsAtGaussPts() {
      return static_cast<VolumeElementForcesAndSourcesCore*>(ptrFE)->coordsAtGaussPts;
    }

    /** \brief coordinate at Gauss points (if hierarchical approximation of element geometry)
      */
    inline MatrixDouble& getHoCoordsAtGaussPts() {
      return static_cast<VolumeElementForcesAndSourcesCore*>(ptrFE)->hoCoordsAtGaussPts;
    }

    inline MatrixDouble& getHoGaussPtsInvJac() {
      return static_cast<VolumeElementForcesAndSourcesCore*>(ptrFE)->hoGaussPtsInvJac;
    }
    inline VectorDouble& getHoGaussPtsDetJac() {
      return static_cast<VolumeElementForcesAndSourcesCore*>(ptrFE)->hoGaussPtsDetJac;
    }

    /** \brief return pointer to Generic Volume Finite Element object
     */
    inline const VolumeElementForcesAndSourcesCore* getVolumeFE() {
      return static_cast<VolumeElementForcesAndSourcesCore*>(ptrFE);
    }

    //differential operators
    PetscErrorCode getDivergenceMatrixOperator_Hdiv(
      int side,EntityType type,DataForcesAndSurcesCore::EntData &data,
      int gg,VectorDouble &div
    );

  };

  int nbGaussPts;
  virtual PetscErrorCode setIntegartionPts();
  virtual PetscErrorCode calculateVolumeAndJacobian();
  virtual PetscErrorCode calculateCoordinatesAtGaussPts();
  virtual PetscErrorCode getSpaceBaseAndOrderOnElement();
  virtual PetscErrorCode calculateBaseFunctionsOnElement();

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

#endif //__VOLUMEELEMENTFORCESANDSOURCESCORE_HPP__

/***************************************************************************//**
 * \defgroup mofem_forces_and_sources_volume_element Volume Element
 * \brief Implementation of general volume element.
 *
 * \ingroup mofem_forces_and_sources
 ******************************************************************************/
