/** \file FaceElementForcesAndSourcesCore.hpp
  \brief Implementation of face element.

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
struct FaceElementForcesAndSourcesCoreBase : public ForcesAndSourcesCore {

  std::string meshPositionsFieldName; ///< Name of the field with geometry

  FaceElementForcesAndSourcesCoreBase(Interface &m_field);

  /** \brief default operator for TRI element
   * \ingroup mofem_forces_and_sources_tri_element
   */
  struct UserDataOperator : public ForcesAndSourcesCore::UserDataOperator {

    using ForcesAndSourcesCore::UserDataOperator::UserDataOperator;

    /**
     * \brief get area of face
     * @return area of face
     */
    inline double getArea() {
      return static_cast<FaceElementForcesAndSourcesCoreBase *>(ptrFE)->aRea;
    }

    /**
     * \brief get measure of element
     * @return area of face
     */
    inline double getMeasure() { return getArea(); }

    /** \brief get triangle normal
     */
    inline VectorDouble &getNormal() {
      return static_cast<FaceElementForcesAndSourcesCoreBase *>(ptrFE)->nOrmal;
    }

    /** \brief get triangle tangent 1
     */
    inline VectorDouble &getTangent1() {
      return static_cast<FaceElementForcesAndSourcesCoreBase *>(ptrFE)
          ->tangentOne;
    }

    /** \brief get triangle tangent 2
     */
    inline VectorDouble &getTangent2() {
      return static_cast<FaceElementForcesAndSourcesCoreBase *>(ptrFE)
          ->tangentTwo;
    }

    /** \brief get normal as tensor
     */
    inline auto getFTensor1Normal() {
      double *ptr = &*getNormal().data().begin();
      return FTensor::Tensor1<double *, 3>(ptr, &ptr[1], &ptr[2]);
    }

    /** \brief get tangentOne as tensor
     */
    inline auto getFTensor1Tangent1() {
      double *ptr = &*getTangent1().data().begin();
      return FTensor::Tensor1<double *, 3>(ptr, &ptr[1], &ptr[2]);
    }

    /** \brief get tangentTwo as tensor
     */
    inline auto getFTensor2Tangent1() {
      double *ptr = &*getTangent2().data().begin();
      return FTensor::Tensor1<double *, 3>(ptr, &ptr[1], &ptr[2]);
    }

    /** \brief get element number of nodes
     */
    inline int getNumNodes() {
      return static_cast<FaceElementForcesAndSourcesCoreBase *>(ptrFE)
          ->num_nodes;
    }

    /** \brief get element connectivity
     */
    inline const EntityHandle *getConn() {
      return static_cast<FaceElementForcesAndSourcesCoreBase *>(ptrFE)->conn;
    }

    /** \brief get triangle coordinates
     */
    inline VectorDouble &getCoords() {
      return static_cast<FaceElementForcesAndSourcesCoreBase *>(ptrFE)->coords;
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

    /** \brief Gauss points and weight, matrix (nb. of points x 3)

    Column 0-2 integration points coordinate x and y, respectively. At rows are
    integration points.

    */
    inline MatrixDouble &getCoordsAtGaussPts() {
      return static_cast<FaceElementForcesAndSourcesCoreBase *>(ptrFE)
          ->coordsAtGaussPts;
    }

    /** \brief get coordinates at Gauss pts.
     */
    inline auto getFTensor1CoordsAtGaussPts() {
      double *ptr = &*getCoordsAtGaussPts().data().begin();
      return FTensor::Tensor1<FTensor::PackPtr<double *, 3>, 3>(ptr, &ptr[1],
                                                                &ptr[2]);
    }

    /** \brief coordinate at Gauss points (if hierarchical approximation of
    element geometry)

    Note: returned matrix has size 0 in rows and columns if no HO approximation
    of geometry is available.

      */
    inline MatrixDouble &getHoCoordsAtGaussPts() {
      return static_cast<FaceElementForcesAndSourcesCoreBase *>(ptrFE)
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

    /** \brief if higher order geometry return normals at Gauss pts.

    Note: returned matrix has size 0 in rows and columns if no HO approximation
    of geometry is available.

     */
    inline MatrixDouble &getNormalsAtGaussPts() {
      return static_cast<FaceElementForcesAndSourcesCoreBase *>(ptrFE)
          ->normalsAtGaussPts;
    }

    /**
     * @deprecated Use getNormalsAtGaussPts
     */
    DEPRECATED inline MatrixDouble &getNormalsAtGaussPt() {
      return getNormalsAtGaussPts();
    }

    /** \brief if higher order geometry return normals at Gauss pts.
     *
     * \param gg gauss point number
     */
    inline ublas::matrix_row<MatrixDouble> getNormalsAtGaussPts(const int gg) {
      return ublas::matrix_row<MatrixDouble>(
          static_cast<FaceElementForcesAndSourcesCoreBase *>(ptrFE)
              ->normalsAtGaussPts,
          gg);
    }

    /**
     * @deprecated Cotrect name is getNormalsAtGaussPts
     */
    DEPRECATED inline ublas::matrix_row<MatrixDouble>
    getNormalsAtGaussPt(const int gg) {
      return getNormalsAtGaussPts(gg);
    }

    /** \brief if higher order geometry return tangent vector to triangle at
    Gauss pts.

    Note: returned matrix has size 0 in rows and columns if no HO approximation
    of geometry is avaliable.

     */
    inline MatrixDouble &getTangent1AtGaussPts() {
      return static_cast<FaceElementForcesAndSourcesCoreBase *>(ptrFE)
          ->tangentOneAtGaussPts;
    }

    /** \brief if higher order geometry return tangent vector to triangle at
    Gauss pts.

    Note: returned matrix has size 0 in rows and columns if no HO approximation
    of geometry is avaliable.

     */
    inline MatrixDouble &getTangent2AtGaussPts() {
      return static_cast<FaceElementForcesAndSourcesCoreBase *>(ptrFE)
          ->tangentTwoAtGaussPts;
    }

    /** \brief get normal at integration points

      Example:
      \code
      double nrm2;
      FTensor::Index<'i',3> i;
      auto t_normal = getFTensor1NormalsAtGaussPts();
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

    /** \brief get tangent 1 at integration points

    */
    inline auto getFTensor1Tangent1AtGaussPts() {
      double *ptr = &*getTangent1AtGaussPts().data().begin();
      return FTensor::Tensor1<FTensor::PackPtr<double *, 3>, 3>(ptr, &ptr[1],
                                                                &ptr[2]);
    }

    /** \brief get tangent 2 at integration points

    */
    inline auto getFTensor1Tangent2AtGaussPts() {
      double *ptr = &*getTangent2AtGaussPts().data().begin();
      return FTensor::Tensor1<FTensor::PackPtr<double *, 3>, 3>(ptr, &ptr[1],
                                                                &ptr[2]);
    }

    /** \brief return pointer to Generic Triangle Finite Element object
     */
    inline const FaceElementForcesAndSourcesCoreBase *getFaceFE() {
      return static_cast<FaceElementForcesAndSourcesCoreBase *>(ptrFE);
    }

    /**
     * @deprecated Use getFaceFE
     */
    DEPRECATED inline const FaceElementForcesAndSourcesCoreBase *
    getFaceElementForcesAndSourcesCore() {
      return getFaceFE();
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

  enum Switches {
    NO_HO_GEOMETRY = 1 << 0,
    NO_CONTRAVARIANT_TRANSFORM_HDIV = 1 << 1,
    NO_COVARIANT_TRANSFORM_HCURL = 1 << 2,
  };

  template <int SWITCH> MoFEMErrorCode OpSwitch() {
    MoFEMFunctionBegin;

    const EntityType type = numeredEntFiniteElementPtr->getEntType();
    if (type != lastEvaluatedElementEntityType) {
      switch (type) {
      case MBTRI:
        getElementPolynomialBase() =
            boost::shared_ptr<BaseFunction>(new TriPolynomialBase());
        break;
      case MBQUAD:
        getElementPolynomialBase() =
            boost::shared_ptr<BaseFunction>(new QuadPolynomialBase());
        break;
      default:
        MoFEMFunctionReturnHot(0);
      }
      CHKERR createDataOnElement();    
    }

    // Calculate normal and tangent vectors for face geometry
    CHKERR calculateAreaAndNormal();
    CHKERR getSpaceBaseAndOrderOnElement();

    CHKERR setIntegrationPts();
    if (gaussPts.size2() == 0)
      MoFEMFunctionReturnHot(0);

    DataForcesAndSourcesCore &data_curl = *dataOnElement[HCURL];
    DataForcesAndSourcesCore &data_div = *dataOnElement[HDIV];

    /// Use the some node base
    CHKERR calculateCoordinatesAtGaussPts();
    CHKERR calculateBaseFunctionsOnElement();

    // FIXME for testing only
    switch (numeredEntFiniteElementPtr->getEntType()) {
    case MBTRI:
      if (!(NO_HO_GEOMETRY & SWITCH))
        CHKERR calculateHoNormal();
      break;
    case MBQUAD:
      CHKERR calculateAreaAndNormalAtIntegrationPts();
      break;
    default:
      SETERRQ(PETSC_COMM_SELF, MOFEM_NOT_IMPLEMENTED,
              "Element type not implemented");
    }

    // Apply Piola transform to HDiv and HCurl spaces, uses previously
    // calculated faces normal and tangent vectors.
    if (!(NO_CONTRAVARIANT_TRANSFORM_HDIV & SWITCH))
      if (dataH1.spacesOnEntities[MBTRI].test(HDIV))
        CHKERR opContravariantTransform.opRhs(data_div);

    if (!(NO_COVARIANT_TRANSFORM_HCURL & SWITCH))
      if (dataH1.spacesOnEntities[MBEDGE].test(HCURL))
        CHKERR opCovariantTransform.opRhs(data_curl);

    // Iterate over operators
    CHKERR loopOverOperators();

    MoFEMFunctionReturn(0);
  }

protected:
  MoFEMErrorCode getNumberOfNodes(int &num_nodes) const;

  /**
   * \brief Calculate element area and normal of the face at integration points
   *
   * @return Error code
   */
  virtual MoFEMErrorCode calculateAreaAndNormalAtIntegrationPts();

  /**
   * \brief Calculate element area and normal of the face
   *
   * Note that at that point is assumed that geometry is exclusively defined by
   * corner nodes.
   *
   * @return Error code
   */
  virtual MoFEMErrorCode calculateAreaAndNormal();

  /**
   * \brief Set integration points
   * @return Error code
   */
  virtual MoFEMErrorCode setIntegrationPts();

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
   * \brief Calculate normal on curved elements
   *
   *  Geometry is given by other field.
   *
   * @return error code
   */
  virtual MoFEMErrorCode calculateHoNormal();

  double aRea;
  int num_nodes;
  const EntityHandle *conn;
  VectorDouble nOrmal, tangentOne, tangentTwo;
  VectorDouble coords;
  MatrixDouble coordsAtGaussPts;

  MatrixDouble hoCoordsAtGaussPts;
  MatrixDouble normalsAtGaussPts;
  MatrixDouble tangentOneAtGaussPts;
  MatrixDouble tangentTwoAtGaussPts;
  OpGetCoordsAndNormalsOnFace opHOCoordsAndNormals;
  OpSetContravariantPiolaTransformOnFace opContravariantTransform;
  OpSetCovariantPiolaTransformOnFace opCovariantTransform;

  friend class UserDataOperator;
  friend class VolumeElementForcesAndSourcesCoreOnSideBase;
};

/** \brief Face finite element switched
 \ingroup mofem_forces_and_sources_tri_element

 */
template <int SWITCH>
struct FaceElementForcesAndSourcesCoreSwitch
    : public FaceElementForcesAndSourcesCoreBase {

  using FaceElementForcesAndSourcesCoreBase::
      FaceElementForcesAndSourcesCoreBase;

  using UserDataOperator =
      FaceElementForcesAndSourcesCoreBase::UserDataOperator;

  MoFEMErrorCode operator()() { return OpSwitch<SWITCH>(); }
};

/** \brief Face finite element default
 \ingroup mofem_forces_and_sources_tri_element

 */
using FaceElementForcesAndSourcesCore =
    FaceElementForcesAndSourcesCoreSwitch<0>;

} // namespace MoFEM

#endif //__FACEELEMENTFORCESANDSOURCESCORE_HPP__

/**
 * \defgroup mofem_forces_and_sources_tri_element Face Element
 * \brief Implementation of face element
 *
 * \ingroup mofem_forces_and_sources
 **/
