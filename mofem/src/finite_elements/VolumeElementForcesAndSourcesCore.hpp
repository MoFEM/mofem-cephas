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
 number of operator added pushing objects to OpPtrVector

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

    UserDataOperator(const FieldSpace space):
    ForcesAndSurcesCore::UserDataOperator(space) {}

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
    inline VectorDouble& getCoords() {
      return static_cast<VolumeElementForcesAndSourcesCore*>(ptrFE)->coords;
    }

    /** \brief matrix of integration (Gauss) points for Volume Element
      *  where columns 0,1,2 are x,y,z coordinates respectively and column 3 is a weight value
      * for example getGaussPts()(1,13) returns y coordinate of 13th Gauss point on particular volume element
      */
    inline MatrixDouble& getGaussPts() {
      return static_cast<VolumeElementForcesAndSourcesCore*>(ptrFE)->gaussPts;
    }

    /** \brief Gauss points and weight, matrix (nb. of points x 3)

    Column 0-2 integration points coordinate x, y and z, respectively. At rows are
    integration points.

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

    /**
     * \brief Get divergence of base functions at integration point
     *
     * Works only for H-div space.
     *
     * @param  side side (local) number of entity on element
     * @param  type type of entity
     * @param  data data structure
     * @param  gg   gauss pts
     * @param  div  divergence vector, size of vector is equal to number of base functions
     * @return      error code
     */
    PetscErrorCode getDivergenceOfHDivBaseFunctions(
      int side,
      EntityType type,
      DataForcesAndSurcesCore::EntData &data,
      int gg,
      VectorDouble &div
    );

    /**
     * \brief Get curl of base functions at integration point
     *
     * Works only for H-curl space.
     *
     * @param  side side (local) number of entity on element
     * @param  type type of entity
     * @param  data data structure
     * @param  gg   gauss pts
     * @param  curl curl matrix, nb. of rows is equal to number of base functions, columns are curl of base vector
     * @return      error code
     */
    PetscErrorCode getCurlOfHCurlBaseFunctions(
      int side,
      EntityType type,
      DataForcesAndSurcesCore::EntData &data,
      int gg,
      MatrixDouble &curl
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

struct FaceElementForcesAndSourcesCore;

/**
 * \brief Volume element used to integrate on skeleton
 * \ingroup mofem_forces_and_sources_volume_element
 */
struct VolumeElementForcesAndSourcesCoreOnSide: public VolumeElementForcesAndSourcesCore {

  FaceElementForcesAndSourcesCore *faceFEPtr;
  VolumeElementForcesAndSourcesCoreOnSide(
    Interface &m_field,const EntityType type = MBTET
  ):
  VolumeElementForcesAndSourcesCore(m_field,type),
  faceFEPtr(NULL) {
  }
  ~VolumeElementForcesAndSourcesCoreOnSide() {}

  inline PetscErrorCode setFaceFEPtr(const FaceElementForcesAndSourcesCore *face_fe_ptr) {
    PetscFunctionBegin;
    faceFEPtr = const_cast<FaceElementForcesAndSourcesCore*>(face_fe_ptr);
    PetscFunctionReturn(0);
  }

  int getRule(int order) { return -1; };

  int faceSense;       ///< Sense of face, could be 1 or -1
  int faceSideNumber;  ///< Face side number
  int faceConnMap[3];
  int tetConnMap[4];
  int oppositeNode;

  PetscErrorCode setGaussPts(int order);

  /** \brief default operator for TET element
    * \ingroup mofem_forces_and_sources_volume_element
    */
  struct UserDataOperator: public VolumeElementForcesAndSourcesCore::UserDataOperator {

    UserDataOperator(
      const std::string &field_name,const char type
    ):
    VolumeElementForcesAndSourcesCore::UserDataOperator(field_name,type) {
    }

    UserDataOperator(
      const std::string &row_field_name,const std::string &col_field_name,const char type
    ):
    VolumeElementForcesAndSourcesCore::UserDataOperator(row_field_name,col_field_name,type) {
    }

    /** \brief return pointer to Generic Volume Finite Element object
     */
    inline const VolumeElementForcesAndSourcesCoreOnSide* getVolumeFE() const {
      return static_cast<VolumeElementForcesAndSourcesCoreOnSide*>(ptrFE);
    }

    inline FaceElementForcesAndSourcesCore* getFaceFEPtr() const {
      return getVolumeFE()->faceFEPtr;
    }

    /**
     * \brief get face sense in respect to volume
     * @return error code
     */
    inline int getFaceSense() const {
      return getVolumeFE()->faceSense;
    }

    /**
     * \brief get face side number in respect to volume
     * @return error code
     */
    inline int getFaceSideNumber() const {
      return getVolumeFE()->faceSideNumber;
    }

    inline bool getEdgeFace(const int ee) const {
      const bool edges_on_faces[6][4] = {
        { true, false, false, true }, // e0
        { false, true, false, true }, // e1
        { false, false, true, true }, // e2
        { true, false, true, false }, // e3
        { true, true, false, false }, // e4
        { false, true, true, false }
      };
      return edges_on_faces[ee][getFaceSideNumber()];
    }

    /**
     * get face normal on side which is this element
     * @return face normal
     */
    VectorDouble& getNormal();

    /** \brief get normal as tensor
    */
    inline FTensor::Tensor1<double*,3> getTensor1Normal() {
      double *ptr = &*getNormal().data().begin();
      return FTensor::Tensor1<double*,3>(ptr,&ptr[1],&ptr[2]);
    }

    /** \brief if higher order geometry return normals at Gauss pts.

    Note: returned matrix has size 0 in rows and columns if no HO approximation
    of geometry is available.

     */
    MatrixDouble& getNormalsAtGaussPt();

    /** \brief if higher order geometry return normals at Gauss pts.
      *
      * \param gg gauss point number
      */
    ublas::matrix_row<MatrixDouble > getNormalsAtGaussPt(const int gg);

    /** \brief get normal at integration points

      Example:
      \code
      double nrm2;
      FTensor::Index<'i',3> i;
      FTensor::Tensor1<double*,3> t_normal = getTensor1NormalsAtGaussPt();
      for(int gg = gg!=data.getN().size1();gg++) {
        nrm2 = sqrt(t_normal(i)*t_normal(i));
        ++t_normal;
      }
      \endcode

    */
    inline FTensor::Tensor1<double*,3> getTensor1NormalsAtGaussPt() {
      double *ptr = &*getNormalsAtGaussPt().data().begin();
      return FTensor::Tensor1<double*,3>(ptr,&ptr[1],&ptr[2],3);
    }

    /** \brief get face coordinates at Gauss pts.

    \note Coordinates should be the same what function getCoordsAtGaussPts
    on tets is returning. If both coordinates are different it is error, or you
    do something very unusual.

     */
    MatrixDouble& getFaceCoordsAtGaussPts();

  };


};

}

#endif //__VOLUMEELEMENTFORCESANDSOURCESCORE_HPP__

/***************************************************************************//**
 * \defgroup mofem_forces_and_sources_volume_element Volume Element
 * \brief Implementation of general volume element.
 *
 * \ingroup mofem_forces_and_sources
 ******************************************************************************/
