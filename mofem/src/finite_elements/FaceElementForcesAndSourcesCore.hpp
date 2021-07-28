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

template <int SWITCH> struct VolumeElementForcesAndSourcesCoreOnSideSwitch;

/** \brief Face finite element
 \ingroup mofem_forces_and_sources_tri_element

 User is implementing own operator at Gauss point level, by own object
 derived from FaceElementForcesAndSourcesCoreL::UserDataOperator.  Arbitrary
 number of operator added pushing objects to OpPtrVector

 */
struct FaceElementForcesAndSourcesCoreBase : public ForcesAndSourcesCore {

  std::string meshPositionsFieldName; ///< Name of the field with geometry

  /** \brief default operator for TRI element
   * \ingroup mofem_forces_and_sources_tri_element
   */
  struct UserDataOperator;
  
  enum Switches {
    NO_HO_GEOMETRY = 1 << 0,
    NO_CONTRAVARIANT_TRANSFORM_HDIV = 1 << 1,
    NO_COVARIANT_TRANSFORM_HCURL = 1 << 2,
  };

  template <int SWITCH> MoFEMErrorCode opSwitch();

protected:
  FaceElementForcesAndSourcesCoreBase(Interface &m_field);

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
   * Note that at that point is assumed that geometry is exclusively defined
   * by corner nodes.
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

  MatrixDouble normalsAtGaussPts;
  MatrixDouble tangentOneAtGaussPts;
  MatrixDouble tangentTwoAtGaussPts;
  OpGetCoordsAndNormalsOnFace opHOCoordsAndNormals;
  OpSetContravariantPiolaTransformOnFace opContravariantTransform;
  OpSetCovariantPiolaTransformOnFace opCovariantTransform;

  friend class UserDataOperator;
  friend class VolumeElementForcesAndSourcesCoreOnSideBase;
};

/** \brief default operator for TRI element
 * \ingroup mofem_forces_and_sources_tri_element
 */
struct FaceElementForcesAndSourcesCoreBase::UserDataOperator
    : public ForcesAndSourcesCore::UserDataOperator {

  using ForcesAndSourcesCore::UserDataOperator::UserDataOperator;

  /**
   * \brief get area of face
   * @return area of face
   */
  inline double getArea();

  /**
   * \brief get measure of element
   * @return area of face
   */
  inline double getMeasure();

  /** \brief get triangle normal
   */
  inline VectorDouble &getNormal();

  /** \brief get triangle tangent 1
   */
  inline VectorDouble &getTangent1();

  /** \brief get triangle tangent 2
   */
  inline VectorDouble &getTangent2();

  /** \brief get normal as tensor
   */
  inline auto getFTensor1Normal();

  /** \brief get tangentOne as tensor
   */
  inline auto getFTensor1Tangent1();

  /** \brief get tangentTwo as tensor
   */
  inline auto getFTensor1Tangent2();

  /** \brief get element number of nodes
   */
  inline int getNumNodes();

  /** \brief get element connectivity
   */
  inline const EntityHandle *getConn();

  /** \brief get triangle coordinates
   */
  inline VectorDouble &getCoords();

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
  inline auto getFTensor1Coords();

  /** \brief if higher order geometry return normals at Gauss pts.

  Note: returned matrix has size 0 in rows and columns if no HO approximation
  of geometry is available.

   */
  inline MatrixDouble &getNormalsAtGaussPts();

  /** \brief if higher order geometry return normals at Gauss pts.
   *
   * \param gg gauss point number
   */
  inline ublas::matrix_row<MatrixDouble> getNormalsAtGaussPts(const int gg);

  /** \brief if higher order geometry return tangent vector to triangle at
  Gauss pts.

  Note: returned matrix has size 0 in rows and columns if no HO approximation
  of geometry is avaliable.

   */
  inline MatrixDouble &getTangent1AtGaussPts();

  /** \brief if higher order geometry return tangent vector to triangle at
  Gauss pts.

  Note: returned matrix has size 0 in rows and columns if no HO approximation
  of geometry is avaliable.

   */
  inline MatrixDouble &getTangent2AtGaussPts();

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
  inline auto getFTensor1NormalsAtGaussPts();

  /** \brief get tangent 1 at integration points

  */
  inline auto getFTensor1Tangent1AtGaussPts();

  /** \brief get tangent 2 at integration points

  */
  inline auto getFTensor1Tangent2AtGaussPts();

  /** \brief return pointer to Generic Triangle Finite Element object
   */
  inline FaceElementForcesAndSourcesCoreBase *getFaceFE();


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
  template <int SWITCH>
  MoFEMErrorCode loopSideVolumes(
      const string &fe_name,
      VolumeElementForcesAndSourcesCoreOnSideSwitch<SWITCH> &fe_method);

private:
  MoFEMErrorCode setPtrFE(ForcesAndSourcesCore *ptr);
};

/** \brief Face finite element switched
 \ingroup mofem_forces_and_sources_tri_element

 */
template <int SWITCH>
struct FaceElementForcesAndSourcesCoreSwitch
    : public FaceElementForcesAndSourcesCoreBase {

  FaceElementForcesAndSourcesCoreSwitch(Interface &m_field)
      : FaceElementForcesAndSourcesCoreBase(m_field) {}

  using UserDataOperator =
      FaceElementForcesAndSourcesCoreBase::UserDataOperator;

  MoFEMErrorCode operator()();
};

/** \brief Face finite element default
 \ingroup mofem_forces_and_sources_tri_element

 */
using FaceElementForcesAndSourcesCore =
    FaceElementForcesAndSourcesCoreSwitch<0>;

template <int SWITCH>
MoFEMErrorCode FaceElementForcesAndSourcesCoreBase::opSwitch() {
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

  CHKERR calculateCoordinatesAtGaussPts();
  CHKERR calHierarchicalBaseFunctionsOnElement();
  CHKERR calBernsteinBezierBaseFunctionsOnElement();
  CHKERR calculateAreaAndNormalAtIntegrationPts();

  if (!(NO_HO_GEOMETRY & SWITCH)) {
    CHKERR calculateHoNormal();
  }

  // Apply Piola transform to HDiv and HCurl spaces, uses previously
  // calculated faces normal and tangent vectors.
  // if (!(NO_CONTRAVARIANT_TRANSFORM_HDIV & SWITCH)) {
  //   if (dataH1.spacesOnEntities[MBTRI].test(HDIV))
  //     CHKERR opContravariantTransform.opRhs(data_div);
  //   if (dataH1.spacesOnEntities[MBQUAD].test(HDIV))
  //     CHKERR opContravariantTransform.opRhs(data_div);
  // }

  // if (!(NO_COVARIANT_TRANSFORM_HCURL & SWITCH)) {
  //   if (dataH1.spacesOnEntities[MBEDGE].test(HCURL))
  //     CHKERR opCovariantTransform.opRhs(data_curl);
  // }

  // Iterate over operators
  CHKERR loopOverOperators();

  MoFEMFunctionReturn(0);
}

template <int SWITCH>
MoFEMErrorCode FaceElementForcesAndSourcesCoreSwitch<SWITCH>::operator()() {
  return opSwitch<SWITCH>();
}

double FaceElementForcesAndSourcesCoreBase::UserDataOperator::getArea() {
  return static_cast<FaceElementForcesAndSourcesCoreBase *>(ptrFE)->aRea;
}

double FaceElementForcesAndSourcesCoreBase::UserDataOperator::getMeasure() {
  return getArea();
}

VectorDouble &
FaceElementForcesAndSourcesCoreBase::UserDataOperator::getNormal() {
  return static_cast<FaceElementForcesAndSourcesCoreBase *>(ptrFE)->nOrmal;
}

VectorDouble &
FaceElementForcesAndSourcesCoreBase::UserDataOperator::getTangent1() {
  return static_cast<FaceElementForcesAndSourcesCoreBase *>(ptrFE)->tangentOne;
}

VectorDouble &
FaceElementForcesAndSourcesCoreBase::UserDataOperator::getTangent2() {
  return static_cast<FaceElementForcesAndSourcesCoreBase *>(ptrFE)->tangentTwo;
}

auto FaceElementForcesAndSourcesCoreBase::UserDataOperator::
    getFTensor1Normal() {
  double *ptr = &*getNormal().data().begin();
  return FTensor::Tensor1<double *, 3>(ptr, &ptr[1], &ptr[2]);
}

auto FaceElementForcesAndSourcesCoreBase::UserDataOperator::
    getFTensor1Tangent1() {
  double *ptr = &*getTangent1().data().begin();
  return FTensor::Tensor1<double *, 3>(ptr, &ptr[1], &ptr[2]);
}

auto FaceElementForcesAndSourcesCoreBase::UserDataOperator::
    getFTensor1Tangent2() {
  double *ptr = &*getTangent2().data().begin();
  return FTensor::Tensor1<double *, 3>(ptr, &ptr[1], &ptr[2]);
}

int FaceElementForcesAndSourcesCoreBase::UserDataOperator::getNumNodes() {
  return static_cast<FaceElementForcesAndSourcesCoreBase *>(ptrFE)->num_nodes;
}

const EntityHandle *
FaceElementForcesAndSourcesCoreBase::UserDataOperator::getConn() {
  return static_cast<FaceElementForcesAndSourcesCoreBase *>(ptrFE)->conn;
}

VectorDouble &
FaceElementForcesAndSourcesCoreBase::UserDataOperator::getCoords() {
  return static_cast<FaceElementForcesAndSourcesCoreBase *>(ptrFE)->coords;
}

auto FaceElementForcesAndSourcesCoreBase::UserDataOperator::
    getFTensor1Coords() {
  double *ptr = &*getCoords().data().begin();
  return FTensor::Tensor1<FTensor::PackPtr<double *, 3>, 3>(ptr, &ptr[1],
                                                            &ptr[2]);
}

MatrixDouble &
FaceElementForcesAndSourcesCoreBase::UserDataOperator::getNormalsAtGaussPts() {
  return static_cast<FaceElementForcesAndSourcesCoreBase *>(ptrFE)
      ->normalsAtGaussPts;
}

ublas::matrix_row<MatrixDouble>
FaceElementForcesAndSourcesCoreBase::UserDataOperator::getNormalsAtGaussPts(
    const int gg) {
  return ublas::matrix_row<MatrixDouble>(
      static_cast<FaceElementForcesAndSourcesCoreBase *>(ptrFE)
          ->normalsAtGaussPts,
      gg);
}

MatrixDouble &
FaceElementForcesAndSourcesCoreBase::UserDataOperator::getTangent1AtGaussPts() {
  return static_cast<FaceElementForcesAndSourcesCoreBase *>(ptrFE)
      ->tangentOneAtGaussPts;
}

MatrixDouble &
FaceElementForcesAndSourcesCoreBase::UserDataOperator::getTangent2AtGaussPts() {
  return static_cast<FaceElementForcesAndSourcesCoreBase *>(ptrFE)
      ->tangentTwoAtGaussPts;
}

auto FaceElementForcesAndSourcesCoreBase::UserDataOperator::
    getFTensor1NormalsAtGaussPts() {
  double *ptr = &*getNormalsAtGaussPts().data().begin();
  return FTensor::Tensor1<FTensor::PackPtr<double *, 3>, 3>(ptr, &ptr[1],
                                                            &ptr[2]);
}

auto FaceElementForcesAndSourcesCoreBase::UserDataOperator::
    getFTensor1Tangent1AtGaussPts() {
  double *ptr = &*getTangent1AtGaussPts().data().begin();
  return FTensor::Tensor1<FTensor::PackPtr<double *, 3>, 3>(ptr, &ptr[1],
                                                            &ptr[2]);
}

auto FaceElementForcesAndSourcesCoreBase::UserDataOperator::
    getFTensor1Tangent2AtGaussPts() {
  double *ptr = &*getTangent2AtGaussPts().data().begin();
  return FTensor::Tensor1<FTensor::PackPtr<double *, 3>, 3>(ptr, &ptr[1],
                                                            &ptr[2]);
}

FaceElementForcesAndSourcesCoreBase *
FaceElementForcesAndSourcesCoreBase::UserDataOperator::getFaceFE() {
  return static_cast<FaceElementForcesAndSourcesCoreBase *>(ptrFE);
}



template <int SWITCH>
MoFEMErrorCode
FaceElementForcesAndSourcesCoreBase::UserDataOperator::loopSideVolumes(
    const string &fe_name,
    VolumeElementForcesAndSourcesCoreOnSideSwitch<SWITCH> &fe_method) {
  return loopSide(fe_name, &fe_method, 3);
}

} // namespace MoFEM

#endif //__FACEELEMENTFORCESANDSOURCESCORE_HPP__

/**
 * \defgroup mofem_forces_and_sources_tri_element Face Element
 * \brief Implementation of face element
 *
 * \ingroup mofem_forces_and_sources
 **/
