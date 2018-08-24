/** \file FaceElementForcesAndSourcesCore.hpp
  \brief Implementation of face element.

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

#ifndef __FACEELEMENTFORCESANDSOURCESCORE_HPP__
#define __FACEELEMENTFORCESANDSOURCESCORE_HPP__

using namespace boost::numeric;

namespace MoFEM {

/** \brief Face finite element
 \ingroup mofem_forces_and_sources_tri_element

 User is implementing own operator at Gauss point level, by own object
 derived from FaceElementForcesAndSourcesCoreL::UserDataOperator.  Arbitrary
 number of operator added pushing objects to OpPtrVector

 */
struct FaceElementForcesAndSourcesCore : public ForcesAndSourcesCore {

  double aRea;
  ;
  int num_nodes;
  const EntityHandle *conn;
  VectorDouble nOrmal, tangentOne, tangentTwo;
  VectorDouble coords;
  MatrixDouble gaussPts;
  MatrixDouble coordsAtGaussPts;

  DataForcesAndSourcesCore dataH1;
  DerivedDataForcesAndSourcesCore derivedDataH1;
  DataForcesAndSourcesCore dataHdiv;
  DerivedDataForcesAndSourcesCore derivedDataHdiv;
  DataForcesAndSourcesCore dataHcurl;
  DerivedDataForcesAndSourcesCore derivedDataHcurl;
  DataForcesAndSourcesCore dataL2;
  DerivedDataForcesAndSourcesCore derivedDataL2;
  DataForcesAndSourcesCore dataNoField, dataNoFieldCol;

  std::string meshPositionsFieldName; ///< Name of the field with geometry

  MatrixDouble hoCoordsAtGaussPts;
  MatrixDouble normalsAtGaussPts;
  MatrixDouble tangentOneAtGaussPts;
  MatrixDouble tangentTwoAtGaussPts;
  OpGetCoordsAndNormalsOnFace opHOCoordsAndNormals;
  OpSetContravariantPiolaTransformOnTriangle opContravariantTransform;
  OpSetCovariantPiolaTransformOnTriangle opCovariantTransform;

  FaceElementForcesAndSourcesCore(Interface &m_field)
      : ForcesAndSourcesCore(m_field), dataH1(MBTRI), derivedDataH1(dataH1),
        dataHdiv(MBTRI), derivedDataHdiv(dataHdiv), dataHcurl(MBTRI),
        derivedDataHcurl(dataHcurl), dataL2(MBTRI), derivedDataL2(dataL2),
        dataNoField(MBTRI), dataNoFieldCol(MBTRI),
        meshPositionsFieldName("MESH_NODE_POSITIONS"),
        opHOCoordsAndNormals(hoCoordsAtGaussPts, normalsAtGaussPts,
                             tangentOneAtGaussPts, tangentTwoAtGaussPts),
        opContravariantTransform(nOrmal, normalsAtGaussPts),
        opCovariantTransform(nOrmal, normalsAtGaussPts, tangentOne,
                             tangentOneAtGaussPts, tangentTwo,
                             tangentTwoAtGaussPts) {}

  /** \brief default operator for TRI element
   * \ingroup mofem_forces_and_sources_tri_element
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
                                                 type, symm){};

    /**
     * \brief get area of face
     * @return area of face
     */
    inline double getArea() {
      return static_cast<FaceElementForcesAndSourcesCore *>(ptrFE)->aRea;
    }

    /**
     * \brief get measure of element
     * @return area of face
     */
    inline double getMeasure() { return getArea(); }

    /** \brief get triangle normal
     */
    inline VectorDouble &getNormal() {
      return static_cast<FaceElementForcesAndSourcesCore *>(ptrFE)->nOrmal;
    }

    /** \brief get triangle tangent 1
     */
    inline VectorDouble &getTangent1() {
      return static_cast<FaceElementForcesAndSourcesCore *>(ptrFE)->tangentOne;
    }

    /** \brief get triangle tangent 2
     */
    inline VectorDouble &getTangent2() {
      return static_cast<FaceElementForcesAndSourcesCore *>(ptrFE)->tangentTwo;
    }

    /** \brief get normal as tensor
     */
    inline auto getFTensor1Normal() {
      double *ptr = &*getNormal().data().begin();
      return FTensor::Tensor1<double *, 3>(ptr, &ptr[1], &ptr[2]);
    }
    
    /// \deprecated use getTensor1Normal()
    DEPRECATED inline auto getTensor1Normal() { return getFTensor1Normal(); }

    /** \brief get tangentOne as tensor
     */
    inline auto getFTensor1Tangent1() {
      double *ptr = &*getTangent1().data().begin();
      return FTensor::Tensor1<double *, 3>(ptr, &ptr[1], &ptr[2]);
    }

    /// \deprecated use getFTensor1Tangent1
    DEPRECATED inline auto getTensor1Tangent1() {
      return getFTensor1Tangent1();
    }

    /** \brief get tangentTwo as tensor
     */
    inline auto getFTensor2Tangent1() {
      double *ptr = &*getTangent2().data().begin();
      return FTensor::Tensor1<double *, 3>(ptr, &ptr[1], &ptr[2]);
    }

    /// \deprecated use getFTensor2Tangent1
    DEPRECATED inline auto getTensor2Tangent1() {
      return getFTensor2Tangent1();
    }

    /** \brief get element number of nodes
     */
    inline int getNumNodes() {
      return static_cast<FaceElementForcesAndSourcesCore *>(ptrFE)->num_nodes;
    }

    /** \brief get element connectivity
     */
    inline const EntityHandle *getConn() {
      return static_cast<FaceElementForcesAndSourcesCore *>(ptrFE)->conn;
    }

    /** \brief get triangle coordinates
     */
    inline VectorDouble &getCoords() {
      return static_cast<FaceElementForcesAndSourcesCore *>(ptrFE)->coords;
    }

    /**
     * \brief get get coords at gauss points

     \code
     FTensor::Index<'i',3> i;
     FTensor::Tensor1<double,3> t_center;
     auto t_coords = getFTensor1Coords();
     t_center(i) = 0;
     for(int nn = 0;nn!=3;nn++) {
        t_center(i) += t_coords(i);
        ++t_coords;
      }
      t_center(i) /= 3;
    \endcode

     */
    inline auto getFTensor1Coords() {
      double *ptr = &*getCoords().data().begin();
      return FTensor::Tensor1<FTensor::PackPtr<double *, 3>, 3>(ptr, &ptr[1],
                                                                &ptr[2]);
    }

    /// \deprecated use getFTensor1Coords
    DEPRECATED inline auto getTensor1Coords() {
      return getFTensor1Coords(); }

    /** \brief get matrix of integration (Gauss) points on Face Element
     *  where columns 0,1 are x,y coordinates respectively and column 2 is a
     * value of weight for example getGaussPts()(1,13) returns y coordinate of
     * 13th Gauss point on particular face element
     */
    inline MatrixDouble &getGaussPts() {
      return static_cast<FaceElementForcesAndSourcesCore *>(ptrFE)->gaussPts;
    }

    /**
     * @brief Get integration weights
     *
     * \code
     * auto t_w = getFTensor0IntegrationWeight();
     * for(int gg = 0; gg!=getGaussPts.size2(); ++gg) {
     *  // integrate something
     *  ++t_w;
     * }
     * \endcode
     *
     * @return FTensor::Tensor0<FTensor::PackPtr<double *, 1>>
     */
    inline auto getFTensor0IntegrationWeight() {
      return FTensor::Tensor0<FTensor::PackPtr<double *, 1>>(
          &(getGaussPts()(2, 0)));
    }

    /** \brief Gauss points and weight, matrix (nb. of points x 3)

    Column 0-2 integration points coordinate x and y, respectively. At rows are
    integration points.

    */
    inline MatrixDouble &getCoordsAtGaussPts() {
      return static_cast<FaceElementForcesAndSourcesCore *>(ptrFE)
          ->coordsAtGaussPts;
    }

    /** \brief get coordinates at Gauss pts.
     */
    inline auto getFTensor1CoordsAtGaussPts() {
      double *ptr = &*getCoordsAtGaussPts().data().begin();
      return FTensor::Tensor1<FTensor::PackPtr<double *, 3>, 3>(ptr, &ptr[1],
                                                                &ptr[2]);
    }

    /// \deprecated use getFTensor1CoordsAtGaussPts
    DEPRECATED inline auto getTensor1CoordsAtGaussPts() {
      return getFTensor1CoordsAtGaussPts();
    }

    /** \brief coordinate at Gauss points (if hierarchical approximation of
    element geometry)

    Note: returned matrix has size 0 in rows and columns if no HO approximation
    of geometry is available.

      */
    inline MatrixDouble &getHoCoordsAtGaussPts() {
      return static_cast<FaceElementForcesAndSourcesCore *>(ptrFE)
          ->hoCoordsAtGaussPts;
    }

    /** \brief get coordinates at Gauss pts (takes in account ho approx. of
     * geometry)
     */
    inline auto getFTensor1HoCoordsAtGaussPts() {
      if (getHoCoordsAtGaussPts().size1() == 0 &&
          getHoCoordsAtGaussPts().size2() != 3) {
        return getFTensor1Coords();
      }
      double *ptr = &*getHoCoordsAtGaussPts().data().begin();
      return FTensor::Tensor1<FTensor::PackPtr<double *, 3>, 3>(ptr, &ptr[1],
                                                                &ptr[2]);
    }

    /// \deprecated use getTensor1HoCoordsAtGaussPts
    DEPRECATED inline auto getTensor1HoCoordsAtGaussPts() {
      return getFTensor1HoCoordsAtGaussPts();
    }

    /** \brief if higher order geometry return normals at Gauss pts.

    Note: returned matrix has size 0 in rows and columns if no HO approximation
    of geometry is available.

     */
    inline MatrixDouble &getNormalsAtGaussPts() {
      return static_cast<FaceElementForcesAndSourcesCore *>(ptrFE)
          ->normalsAtGaussPts;
    }

    /// \deprecated use getNormalsAtGaussPts
    DEPRECATED inline MatrixDouble &getNormalsAtGaussPt() {
      return getNormalsAtGaussPts();
    }

    /** \brief if higher order geometry return normals at Gauss pts.
     *
     * \param gg gauss point number
     */
    inline ublas::matrix_row<MatrixDouble> getNormalsAtGaussPts(const int gg) {
      return ublas::matrix_row<MatrixDouble>(
          static_cast<FaceElementForcesAndSourcesCore *>(ptrFE)
              ->normalsAtGaussPts,
          gg);
    }

    /// \deprecated use getNormalsAtGaussPts
    DEPRECATED inline auto getNormalsAtGaussPt(const int gg) {
      return getNormalsAtGaussPts(gg);
    }

    /** \brief if higher order geometry return tangent vector to triangle at
    Gauss pts.

    Note: returned matrix has size 0 in rows and columns if no HO approximation
    of geometry is avaliable.

     */
    inline MatrixDouble &getTangent1AtGaussPts() {
      return static_cast<FaceElementForcesAndSourcesCore *>(ptrFE)
          ->tangentOneAtGaussPts;
    }

    /// \deprecated use getTangent1AtGaussPts
    DEPRECATED inline MatrixDouble &getTangent1AtGaussPt() {
      return getTangent1AtGaussPts();
    }

    /** \brief if higher order geometry return tangent vector to triangle at
    Gauss pts.

    Note: returned matrix has size 0 in rows and columns if no HO approximation
    of geometry is avaliable.

     */
    inline MatrixDouble &getTangent2AtGaussPts() {
      return static_cast<FaceElementForcesAndSourcesCore *>(ptrFE)
          ->tangentTwoAtGaussPts;
    }

    /// \deprecated use getTangent2AtGaussPts
    DEPRECATED inline MatrixDouble &getTangent2AtGaussPt() {
      return getTangent2AtGaussPts();
    }

    /** \brief get normal at integration points

      Example:
      \code
      double nrm2;
      FTensor::Index<'i',3> i;
      FTensor::Tensor1<double*,3> t_normal = getFTensor1NormalsAtGaussPts();
      for(int gg = gg!=data.getN().size1();gg++) {
        nrm2 = sqrt(t_normal(i)*t_normal(i));
        ++t_normal;
      }
      \endcode

    */
    inline auto getFTensor1NormalsAtGaussPts() {
      double *ptr = &*getNormalsAtGaussPts().data().begin();
      return FTensor::Tensor1<FTensor::PackPtr<double *, 3>, 3>(ptr, &ptr[1],
                                                                &ptr[2]);
    }

    /// \deprecated use getFTensor1NormalsAtGaussPt
    DEPRECATED inline auto getTensor1NormalsAtGaussPt() {
      return getFTensor1NormalsAtGaussPts();
    }

    /** \brief get tangent 1 at integration points

    */
    inline auto getFTensor1Tangent1AtGaussPts() {
      double *ptr = &*getTangent1AtGaussPts().data().begin();
      return FTensor::Tensor1<FTensor::PackPtr<double *, 3>, 3>(ptr, &ptr[1],
                                                                &ptr[2]);
    }

    /// \deprecated use getFTensor1Tangent1AtGaussPt
    DEPRECATED inline auto getTensor1Tangent1AtGaussPt() {
      return getFTensor1Tangent1AtGaussPts();
    }

    /** \brief get tangent 2 at integration points

    */
    inline auto getFTensor1Tangent2AtGaussPts() {
      double *ptr = &*getTangent2AtGaussPts().data().begin();
      return FTensor::Tensor1<FTensor::PackPtr<double *, 3>, 3>(ptr, &ptr[1],
                                                                &ptr[2]);
    }

    /// \deprecated use getFTensor1Tangent2AtGaussPt
    DEPRECATED inline auto getTensor1Tangent2AtGaussPt() {
      return getFTensor1Tangent2AtGaussPts();
    }

    /** \deprecated use getFaceFE instead
     */
    DEPRECATED inline const FaceElementForcesAndSourcesCore *
    getFaceElementForcesAndSourcesCore() {
      return getFaceFE();
    }

    /** \deprecated use getFaceFE instead
     */
    DEPRECATED inline const FaceElementForcesAndSourcesCore *getTriFE() {
      return getFaceFE();
    }

    /** \brief return pointer to Generic Triangle Finite Element object
     */
    inline const FaceElementForcesAndSourcesCore *getFaceFE() {
      return static_cast<FaceElementForcesAndSourcesCore *>(ptrFE);
    }

    /**
     *
     * User call this function to loop over elements on the side of face. This
     * function calls MoFEM::VolumeElementForcesAndSourcesCoreOnSide with is
     * operator to do calculations.
     *
     * @param  fe_name Name of the element
     * @param  method  Finite element object
     * @return         error code
     */
    MoFEMErrorCode
    loopSideVolumes(const string &fe_name,
                    VolumeElementForcesAndSourcesCoreOnSide &method);
  };

  /**
   * \brief Calculate element area and normal of the face
   *
   * Note that at that point is assumed that geometry is exclusively defined by
   * corner nodes.
   *
   * @return Error code
   */
  virtual MoFEMErrorCode calculateAreaAndNormal();

  int nbGaussPts; ///< Number of integration points

  /**
   * \brief Set integration points
   * @return Error code
   */
  virtual MoFEMErrorCode setIntegrationPts();

  /// \deprecated method with spelling mistake, use setIntegrationPts
  DEPRECATED MoFEMErrorCode setIntegartionPts() { return setIntegrationPts(); }

  /**
   * \brief Determine approximation space and order of base functions
   * @return Error code
   */
  virtual MoFEMErrorCode getSpaceBaseAndOrderOnElement();

  /**
   * \brief Calculate coordinate at integration points
   * @return Error code
   */
  virtual MoFEMErrorCode calculateCoordinatesAtGaussPts();

  /**
   * \brief Calculate base functions
   * @return Error code
   */
  virtual MoFEMErrorCode calculateBaseFunctionsOnElement();

  /**
   * \brief Calculate normal on curved elements
   *
   *  Geometry is given by other field.
   *
   * @return error code
   */
  virtual MoFEMErrorCode calculateHoNormal();

  MoFEMErrorCode operator()();
};

/** \brief Calculate inverse of jacobian for face element

  It is assumed that face element is XY plane. Applied
  only for 2d problems.

  \todo Generalize function for arbitrary face orientation in 3d space

  \ingroup mofem_forces_and_sources_tri_element

*/
struct OpCalculateInvJacForFace
    : public FaceElementForcesAndSourcesCore::UserDataOperator {
  MatrixDouble &invJac;

  // /**
  //  * \deprecated Field name do not needed to construct class, change v0.5.17.
  //  */
  // DEPRECATED OpCalculateInvJacForFace(const std::string
  // &field_name,MatrixDouble &inv_jac):
  // FaceElementForcesAndSourcesCore::UserDataOperator(H1),
  // invJac(inv_jac) {}

  OpCalculateInvJacForFace(MatrixDouble &inv_jac)
      : FaceElementForcesAndSourcesCore::UserDataOperator(H1), invJac(inv_jac) {
  }

  MoFEMErrorCode doWork(int side, EntityType type,
                        DataForcesAndSourcesCore::EntData &data);
};

/** \brief Transform local reference derivatives of shape functions to global
derivatives

\ingroup mofem_forces_and_sources_tri_element

*/
struct OpSetInvJacH1ForFace
    : public FaceElementForcesAndSourcesCore::UserDataOperator {
  MatrixDouble &invJac;

  // /**
  //  * \deprecated Field name do not needed to construct class, change v0.5.17.
  //  */
  // DEPRECATED OpSetInvJacH1ForFace(const std::string &field_name,MatrixDouble
  // &inv_jac): FaceElementForcesAndSourcesCore::UserDataOperator(H1),
  // invJac(inv_jac) {}

  OpSetInvJacH1ForFace(MatrixDouble &inv_jac)
      : FaceElementForcesAndSourcesCore::UserDataOperator(H1), invJac(inv_jac) {
  }

  MatrixDouble diffNinvJac;
  MoFEMErrorCode doWork(int side, EntityType type,
                        DataForcesAndSourcesCore::EntData &data);
};

/**
 * \brief brief Transform local reference derivatives of shape function to
 global derivatives for face

 * \ingroup mofem_forces_and_sources_tri_element
 */
struct OpSetInvJacHcurlFace
    : public FaceElementForcesAndSourcesCore::UserDataOperator {

  MatrixDouble &invJac;

  // /**
  //  * \deprecated Field name do not needed to construct class, change v0.5.17.
  //  */
  // DEPRECATED OpSetInvJacHcurlFace(const std::string &field_name,MatrixDouble
  // &inv_jac): FaceElementForcesAndSourcesCore::UserDataOperator(HCURL),
  // invJac(inv_jac) {
  // }

  OpSetInvJacHcurlFace(MatrixDouble &inv_jac)
      : FaceElementForcesAndSourcesCore::UserDataOperator(HCURL),
        invJac(inv_jac) {}

  MatrixDouble diffHcurlInvJac;

  MoFEMErrorCode doWork(int side, EntityType type,
                        DataForcesAndSourcesCore::EntData &data);
};

/// \deprecated Use FaceElementForcesAndSourcesCore
DEPRECATED typedef FaceElementForcesAndSourcesCore
    FaceElementForcesAndSurcesCore;

} // namespace MoFEM

#endif //__FACEELEMENTFORCESANDSOURCESCORE_HPP__

/**
 * \defgroup mofem_forces_and_sources_tri_element Face Element
 * \brief Implementation of face element
 * 
 * \ingroup mofem_forces_and_sources
 ******************************************************************************/
