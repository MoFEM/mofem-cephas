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

 User is implementing own operator at Gauss point level, by class
 derived from VolumeElementForcesAndSourcesCore::UserDataOperator. Arbitrary
 number of operator can be added by pushing objects to OpPtrVector

 */
struct VolumeElementForcesAndSourcesCore : public ForcesAndSourcesCore {

  VectorDouble coords;
  MatrixDouble3by3 jAc;
  MatrixDouble3by3 invJac;

  DataForcesAndSourcesCore dataH1;
  DerivedDataForcesAndSourcesCore derivedDataH1;
  DataForcesAndSourcesCore dataL2;
  DerivedDataForcesAndSourcesCore derivedDataL2;
  DataForcesAndSourcesCore dataHdiv;
  DerivedDataForcesAndSourcesCore derivedDataHdiv;
  DataForcesAndSourcesCore dataHcurl;
  DerivedDataForcesAndSourcesCore derivedDataHcurl;
  DataForcesAndSourcesCore dataNoField;
  DataForcesAndSourcesCore dataNoFieldCol;

  OpSetInvJacH1 opSetInvJacH1;
  OpSetContravariantPiolaTransform opContravariantPiolaTransform;
  OpSetCovariantPiolaTransform opCovariantPiolaTransform;
  OpSetInvJacHdivAndHcurl opSetInvJacHdivAndHcurl;

  std::string meshPositionsFieldName;
  MatrixDouble hoCoordsAtGaussPts;
  MatrixDouble hoGaussPtsJac;
  MatrixDouble hoGaussPtsInvJac;
  VectorDouble hoGaussPtsDetJac;

  OpGetDataAndGradient<3, 3>
      opHOatGaussPoints; ///< higher order geometry data at Gauss pts
  OpSetHoInvJacH1 opSetHoInvJacH1;
  OpSetHoContravariantPiolaTransform opHoContravariantTransform;
  OpSetHoCovariantPiolaTransform opHoCovariantTransform;
  OpSetHoInvJacHdivAndHcurl opSetHoInvJacHdivAndHcurl;

  VolumeElementForcesAndSourcesCore(Interface &m_field,
                                    const EntityType type = MBTET);
  virtual ~VolumeElementForcesAndSourcesCore() {}

  double vOlume;

  int num_nodes;
  const EntityHandle *conn;
  FTensor::Tensor2<double *, 3, 3> tJac;
  FTensor::Tensor2<double *, 3, 3> tInvJac;

  MatrixDouble gaussPts;
  MatrixDouble coordsAtGaussPts;

  /** \brief default operator for TET element
   * \ingroup mofem_forces_and_sources_volume_element
   */
  struct UserDataOperator : public ForcesAndSourcesCore::UserDataOperator {

    UserDataOperator(const FieldSpace space)
        : ForcesAndSourcesCore::UserDataOperator(space) {}

    UserDataOperator(const std::string &field_name, const char type)
        : ForcesAndSourcesCore::UserDataOperator(field_name, type) {}

    UserDataOperator(const std::string &row_field_name,
                     const std::string &col_field_name, const char type,
                     const bool symm = true)
        : ForcesAndSourcesCore::UserDataOperator(row_field_name, col_field_name,
                                                 type, symm) {}

    /** \brief get element number of nodes
     */
    inline int getNumNodes() {
      return static_cast<VolumeElementForcesAndSourcesCore *>(ptrFE)->num_nodes;
    }

    /** \brief get element connectivity
     */
    inline const EntityHandle *getConn() {
      return static_cast<VolumeElementForcesAndSourcesCore *>(ptrFE)->conn;
    }

    /** \brief element volume (linear geometry)
     */
    inline double getVolume() {
      return static_cast<VolumeElementForcesAndSourcesCore *>(ptrFE)->vOlume;
    }

    /**
     * \brief get element Jacobian
     */
    inline FTensor::Tensor2<double *, 3, 3> &getJac() {
      return static_cast<VolumeElementForcesAndSourcesCore *>(ptrFE)->tJac;
    }

    /**
     * \brief get element inverse Jacobian
     */
    inline FTensor::Tensor2<double *, 3, 3> &getInvJac() {
      return static_cast<VolumeElementForcesAndSourcesCore *>(ptrFE)->tInvJac;
    }

    /**
     * \brief get measure of element
     * @return area of face
     */
    inline double getMeasure() { return getVolume(); }

    /** \brief nodal coordinates
     */
    inline VectorDouble &getCoords() {
      return static_cast<VolumeElementForcesAndSourcesCore *>(ptrFE)->coords;
    }

    /** \brief matrix of integration (Gauss) points for Volume Element
     *  where columns 0,1,2 are x,y,z coordinates respectively and column 3 is a
     * weight value for example getGaussPts()(1,13) returns y coordinate of 13th
     * Gauss point on particular volume element
     */
    inline MatrixDouble &getGaussPts() {
      return static_cast<VolumeElementForcesAndSourcesCore *>(ptrFE)->gaussPts;
    }

    inline FTensor::Tensor0<double *> getFTensor0IntegrationWeight() {
      return FTensor::Tensor0<double *>(&(getGaussPts()(3, 0)), 1);
    }

    /** \brief Gauss points and weight, matrix (nb. of points x 3)

    Column 0-2 integration points coordinate x, y and z, respectively. At rows
    are integration points.

    */
    inline MatrixDouble &getCoordsAtGaussPts() {
      return static_cast<VolumeElementForcesAndSourcesCore *>(ptrFE)
          ->coordsAtGaussPts;
    }

    /** \brief coordinate at Gauss points (if hierarchical approximation of
     * element geometry)
     */
    inline MatrixDouble &getHoCoordsAtGaussPts() {
      return static_cast<VolumeElementForcesAndSourcesCore *>(ptrFE)
          ->hoCoordsAtGaussPts;
    }

    inline MatrixDouble &getHoGaussPtsInvJac() {
      return static_cast<VolumeElementForcesAndSourcesCore *>(ptrFE)
          ->hoGaussPtsInvJac;
    }

    inline VectorDouble &getHoGaussPtsDetJac() {
      return static_cast<VolumeElementForcesAndSourcesCore *>(ptrFE)
          ->hoGaussPtsDetJac;
    }

    inline FTensor::Tensor0<double *> getFTenosr0HoMeasure() {
      return FTensor::Tensor0<double *>(&*getHoGaussPtsDetJac().data().begin());
    }

    /**
     * \brief Get coordinates at integration points assuming linear geometry
     *
     * \code
     * FTensor::Tensor1<double*> t_coords = getTensor1CoordsAtGaussPts();
     * for(int gg = 0;gg!=nb_int_ptrs;gg++) {
     *   // do something
     *   ++t_coords;
     * }
     * \endcode
     *
     */
    inline FTensor::Tensor1<double *, 3> getTensor1CoordsAtGaussPts() {
      return FTensor::Tensor1<double *, 3>(&getCoordsAtGaussPts()(0, 0),
                                           &getCoordsAtGaussPts()(0, 1),
                                           &getCoordsAtGaussPts()(0, 2), 3);
    }

    /**
     * \brief Get coordinates at integration points taking geometry from field
     *
     * This is HO geometry given by arbitrary order polynomial
     * \code
     * FTensor::Tensor1<double*> t_coords = getTensor1HoCoordsAtGaussPts();
     * for(int gg = 0;gg!=nb_int_ptrs;gg++) {
     *   // do something
     *   ++t_coords;
     * }
     * \endcode
     *
     */
    inline FTensor::Tensor1<double *, 3> getFTensor1HoCoordsAtGaussPts() {
      return FTensor::Tensor1<double *, 3>(&getHoCoordsAtGaussPts()(0, 0),
                                           &getHoCoordsAtGaussPts()(0, 1),
                                           &getHoCoordsAtGaussPts()(0, 2), 3);
    }

    /** \brief return pointer to Generic Volume Finite Element object
     */
    inline const VolumeElementForcesAndSourcesCore *getVolumeFE() {
      return static_cast<VolumeElementForcesAndSourcesCore *>(ptrFE);
    }

    /**
     * \brief Get divergence of base functions at integration point
     *
     * Works only for H-div space.
     *
     * How to use it:
     * \code
     * VectorDouble div_vec(data.getFieldData().size());
     * for(int gg = 0;gg<nb_gauss_pts;gg++) {
     *  ierr = getDivergenceOfHDivBaseFunctions(side,type,data,gg,div_vec);
     * CHKERRG(ierr);
     *  // do something with vec_div
     * }
     * \endcode
     * where vec_div has size of nb. of dofs, at each element represents
     * divergence of base functions.
     *
     * @param  side side (local) number of entity on element
     * @param  type type of entity
     * @param  data data structure
     * @param  gg   gauss pts
     * @param  div  divergence vector, size of vector is equal to number of base
     * functions
     * @return      error code
     */
    MoFEMErrorCode
    getDivergenceOfHDivBaseFunctions(int side, EntityType type,
                                     DataForcesAndSourcesCore::EntData &data,
                                     int gg, VectorDouble &div);

    /**
     * \brief Get curl of base functions at integration point
     *
     * \f[
     * \nabla \times \mathbf{\phi}
     * \f]
     *
     * Works only for H-curl and H-div space.
     *
     * How to use it:
     * \code
     * MatrixDouble curl_mat(data.getFieldData().size(),3);
     * for(int gg = 0;gg<nb_gauss_pts;gg++) {
     *  ierr = getCurlOfHCurlBaseFunctions(side,type,data,gg,curl_mat);
     * CHKERRG(ierr);
     *  // do something with vec_div
     * }
     * \endcode
     * where curl_mat is matrix which number of rows is equal to nb. of base
     * functions. Number of columns is 3, since we work in 3d here. Rows
     * represents curl of vase functions.
     *
     * @param  side side (local) number of entity on element
     * @param  type type of entity
     * @param  data data structure
     * @param  gg   gauss pts
     * @param  curl curl matrix, nb. of rows is equal to number of base
     * functions, columns are curl of base vector
     * @return      error code
     */
    MoFEMErrorCode
    getCurlOfHCurlBaseFunctions(int side, EntityType type,
                                DataForcesAndSourcesCore::EntData &data, int gg,
                                MatrixDouble &curl);
  };

  unsigned int nbGaussPts; ///< Number of integration points

  // Note that functions below could be overloaded by user to change default
  // behavior of the element.

  /**
   * \brief Set integration points
   * @return Error code
   */
  virtual MoFEMErrorCode setIntegrationPts();

  /// \deprecated function with spelling mistake, use setIntegrationPts
  DEPRECATED virtual MoFEMErrorCode setIntegartionPts() {
    return setIntegrationPts();
  }

  /**
   * \brief Calculate element volume and Jacobian
   *
   * Note that at that point is assumed that geometry is exclusively defined by
   * corner nodes.
   *
   * @return Error code
   */
  virtual MoFEMErrorCode calculateVolumeAndJacobian();

  /**
   * \brief Calculate coordinate at integration points
   * @return Error code
   */
  virtual MoFEMErrorCode calculateCoordinatesAtGaussPts();

  /**
   * \brief Determine approximation space and order of base functions
   * @return Error code
   */
  virtual MoFEMErrorCode getSpaceBaseAndOrderOnElement();

  /**
   * \brief Calculate base functions
   * @return Error code
   */
  virtual MoFEMErrorCode calculateBaseFunctionsOnElement(const int b);

  /**
   * \brief Calculate base functions
   * @return Error code
   */
  virtual MoFEMErrorCode calculateBaseFunctionsOnElement();

  /**
   * \brief Transform base functions based on geometric element Jacobian.
   *
   * This function apply transformation to base functions and its derivatives.
   * For example when base functions for H-div are present the
   * Piola-Transformarion is applied to base functions and their derivatives.
   *
   * @return Error code
   */
  virtual MoFEMErrorCode transformBaseFunctions();

  /** \brief Calculate Jacobian for HO geometry
   *
   * MoFEM use hierarchical approximate base to describe geometry of the body.
   * This function transform derivatives of base functions when HO geometry is
   * set and calculate Jacobian, inverse of Jacobian and determinant of
   * transformation.
   *
   */
  virtual MoFEMErrorCode calculateHoJacobian();

  /**
   * \brief Transform base functions based on ho-geometry element Jacobian.
   *
   * This function apply transformation to base functions and its derivatives.
   * For example when base functions for H-div are present the
   * Piola-Transformarion is applied to base functions and their derivatives.
   *
   * @return Error code
   */
  virtual MoFEMErrorCode transformHoBaseFunctions();

  MoFEMErrorCode operator()();
};

struct FaceElementForcesAndSourcesCore;

/**
 * \brief Volume element used to integrate on skeleton
 * \ingroup mofem_forces_and_sources_volume_element
 */
struct VolumeElementForcesAndSourcesCoreOnSide
    : public VolumeElementForcesAndSourcesCore {

  FaceElementForcesAndSourcesCore *faceFEPtr;
  VolumeElementForcesAndSourcesCoreOnSide(Interface &m_field,
                                          const EntityType type = MBTET)
      : VolumeElementForcesAndSourcesCore(m_field, type), faceFEPtr(NULL) {}
  ~VolumeElementForcesAndSourcesCoreOnSide() {}

  inline MoFEMErrorCode
  setFaceFEPtr(const FaceElementForcesAndSourcesCore *face_fe_ptr) {
    MoFEMFunctionBeginHot;
    faceFEPtr = const_cast<FaceElementForcesAndSourcesCore *>(face_fe_ptr);
    MoFEMFunctionReturnHot(0);
  }

  int getRule(int order) { return -1; };

  int faceSense;      ///< Sense of face, could be 1 or -1
  int faceSideNumber; ///< Face side number
  int faceConnMap[3];
  int tetConnMap[4];
  int oppositeNode;

  MoFEMErrorCode setGaussPts(int order);

  /** \brief default operator for TET element
   * \ingroup mofem_forces_and_sources_volume_element
   */
  struct UserDataOperator
      : public VolumeElementForcesAndSourcesCore::UserDataOperator {

    UserDataOperator(const std::string &field_name, const char type)
        : VolumeElementForcesAndSourcesCore::UserDataOperator(field_name,
                                                              type) {}

    UserDataOperator(const std::string &row_field_name,
                     const std::string &col_field_name, const char type)
        : VolumeElementForcesAndSourcesCore::UserDataOperator(
              row_field_name, col_field_name, type) {}

    /** \brief return pointer to Generic Volume Finite Element object
     */
    inline const VolumeElementForcesAndSourcesCoreOnSide *getVolumeFE() const {
      return static_cast<VolumeElementForcesAndSourcesCoreOnSide *>(ptrFE);
    }

    inline FaceElementForcesAndSourcesCore *getFaceFEPtr() const {
      return getVolumeFE()->faceFEPtr;
    }

    /**
     * \brief get face sense in respect to volume
     * @return error code
     */
    inline int getFaceSense() const { return getVolumeFE()->faceSense; }

    /**
     * \brief get face side number in respect to volume
     * @return error code
     */
    inline int getFaceSideNumber() const {
      return getVolumeFE()->faceSideNumber;
    }

    inline bool getEdgeFace(const int ee) const {
      const bool edges_on_faces[6][4] = {{true, false, false, true}, // e0
                                         {false, true, false, true}, // e1
                                         {false, false, true, true}, // e2
                                         {true, false, true, false}, // e3
                                         {true, true, false, false}, // e4
                                         {false, true, true, false}};
      return edges_on_faces[ee][getFaceSideNumber()];
    }

    /**
     * get face normal on side which is this element
     * @return face normal
     */
    VectorDouble &getNormal();

    /** \brief get normal as tensor
     */
    inline FTensor::Tensor1<double *, 3> getTensor1Normal() {
      double *ptr = &*getNormal().data().begin();
      return FTensor::Tensor1<double *, 3>(ptr, &ptr[1], &ptr[2]);
    }

    /** \brief if higher order geometry return normals at Gauss pts.

    Note: returned matrix has size 0 in rows and columns if no HO approximation
    of geometry is available.

     */
    MatrixDouble &getNormalsAtGaussPt();

    /** \brief if higher order geometry return normals at Gauss pts.
     *
     * \param gg gauss point number
     */
    ublas::matrix_row<MatrixDouble> getNormalsAtGaussPt(const int gg);

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
    inline FTensor::Tensor1<double *, 3> getTensor1NormalsAtGaussPt() {
      double *ptr = &*getNormalsAtGaussPt().data().begin();
      return FTensor::Tensor1<double *, 3>(ptr, &ptr[1], &ptr[2], 3);
    }

    /** \brief get face coordinates at Gauss pts.

    \note Coordinates should be the same what function getCoordsAtGaussPts
    on tets is returning. If both coordinates are different it is error, or you
    do something very unusual.

     */
    MatrixDouble &getFaceCoordsAtGaussPts();
  };
};

/// \deprecated Use VolumeElementForcesAndSourcesCore
DEPRECATED typedef VolumeElementForcesAndSourcesCore
    VolumeElementForcesAndSurcesCore;

} // namespace MoFEM

#endif //__VOLUMEELEMENTFORCESANDSOURCESCORE_HPP__

/**
 * \defgroup mofem_forces_and_sources_volume_element Volume Element
 * \brief Implementation of general volume element.
 *
 * \ingroup mofem_forces_and_sources
 **/
