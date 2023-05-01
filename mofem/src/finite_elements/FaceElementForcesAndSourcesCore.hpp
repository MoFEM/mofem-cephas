/** \file FaceElementForcesAndSourcesCore.hpp
  \brief Implementation of face element.

*/

#ifndef __FACEELEMENTFORCESANDSOURCESCORE_HPP__
#define __FACEELEMENTFORCESANDSOURCESCORE_HPP__

using namespace boost::numeric;

namespace MoFEM {

struct VolumeElementForcesAndSourcesCoreOnSide;

/** \brief Face finite element
 \ingroup mofem_forces_and_sources_tri_element

 User is implementing own operator at Gauss point level, by own object
 derived from FaceElementForcesAndSourcesCoreL::UserDataOperator.  Arbitrary
 number of operator added pushing objects to OpPtrVector

 */
struct FaceElementForcesAndSourcesCore : public ForcesAndSourcesCore {

  /**
   * @deprecated not used anumore, will be removed in next versions
   *
   */
  std::string meshPositionsFieldName;

  /** \brief default operator for TRI element
   * \ingroup mofem_forces_and_sources_tri_element
   */
  struct UserDataOperator;

  MoFEMErrorCode operator()();

  FaceElementForcesAndSourcesCore(Interface &m_field);

protected:
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

  double aRea;
  double elementCircumDiam;
  int num_nodes;
  const EntityHandle *conn;
  VectorDouble nOrmal, tangentOne, tangentTwo;
  VectorDouble coords;

  MatrixDouble normalsAtGaussPts;
  MatrixDouble tangentOneAtGaussPts;
  MatrixDouble tangentTwoAtGaussPts;

  friend class UserDataOperator;
  friend class VolumeElementForcesAndSourcesCoreOnSide;
};

/** \brief default operator for TRI element
 * \ingroup mofem_forces_and_sources_tri_element
 */
struct FaceElementForcesAndSourcesCore::UserDataOperator
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

  /**
   * \brief get circumference of element (or largest diagonal)
   * @return of face
   */
  inline double getCircumDiam();

  /**
   * \brief get measure of element
   * @return area of face
   */
  inline double getElementCharacteristicLength();

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
  inline FaceElementForcesAndSourcesCore *getFaceFE();

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
  loopSideVolumes(const string fe_name,
                  VolumeElementForcesAndSourcesCoreOnSide &fe_method);

private:
  MoFEMErrorCode setPtrFE(ForcesAndSourcesCore *ptr);
};

double FaceElementForcesAndSourcesCore::UserDataOperator::getArea() {
  return static_cast<FaceElementForcesAndSourcesCore *>(ptrFE)->aRea;
}

double FaceElementForcesAndSourcesCore::UserDataOperator::getMeasure() {
  return getArea();
}

double FaceElementForcesAndSourcesCore::UserDataOperator::getCircumDiam() {
  return static_cast<FaceElementForcesAndSourcesCore *>(ptrFE)
      ->elementCircumDiam;
}

double FaceElementForcesAndSourcesCore::UserDataOperator::
    getElementCharacteristicLength() {
  return getCircumDiam();
}

VectorDouble &FaceElementForcesAndSourcesCore::UserDataOperator::getNormal() {
  return static_cast<FaceElementForcesAndSourcesCore *>(ptrFE)->nOrmal;
}

VectorDouble &FaceElementForcesAndSourcesCore::UserDataOperator::getTangent1() {
  return static_cast<FaceElementForcesAndSourcesCore *>(ptrFE)->tangentOne;
}

VectorDouble &FaceElementForcesAndSourcesCore::UserDataOperator::getTangent2() {
  return static_cast<FaceElementForcesAndSourcesCore *>(ptrFE)->tangentTwo;
}

auto FaceElementForcesAndSourcesCore::UserDataOperator::getFTensor1Normal() {
  double *ptr = &*getNormal().data().begin();
  return FTensor::Tensor1<double *, 3>(ptr, &ptr[1], &ptr[2]);
}

auto FaceElementForcesAndSourcesCore::UserDataOperator::getFTensor1Tangent1() {
  double *ptr = &*getTangent1().data().begin();
  return FTensor::Tensor1<double *, 3>(ptr, &ptr[1], &ptr[2]);
}

auto FaceElementForcesAndSourcesCore::UserDataOperator::getFTensor1Tangent2() {
  double *ptr = &*getTangent2().data().begin();
  return FTensor::Tensor1<double *, 3>(ptr, &ptr[1], &ptr[2]);
}

int FaceElementForcesAndSourcesCore::UserDataOperator::getNumNodes() {
  return static_cast<FaceElementForcesAndSourcesCore *>(ptrFE)->num_nodes;
}

const EntityHandle *
FaceElementForcesAndSourcesCore::UserDataOperator::getConn() {
  return static_cast<FaceElementForcesAndSourcesCore *>(ptrFE)->conn;
}

VectorDouble &FaceElementForcesAndSourcesCore::UserDataOperator::getCoords() {
  return static_cast<FaceElementForcesAndSourcesCore *>(ptrFE)->coords;
}

auto FaceElementForcesAndSourcesCore::UserDataOperator::getFTensor1Coords() {
  double *ptr = &*getCoords().data().begin();
  return FTensor::Tensor1<FTensor::PackPtr<double *, 3>, 3>(ptr, &ptr[1],
                                                            &ptr[2]);
}

MatrixDouble &
FaceElementForcesAndSourcesCore::UserDataOperator::getNormalsAtGaussPts() {
  return static_cast<FaceElementForcesAndSourcesCore *>(ptrFE)
      ->normalsAtGaussPts;
}

ublas::matrix_row<MatrixDouble>
FaceElementForcesAndSourcesCore::UserDataOperator::getNormalsAtGaussPts(
    const int gg) {
  return ublas::matrix_row<MatrixDouble>(
      static_cast<FaceElementForcesAndSourcesCore *>(ptrFE)->normalsAtGaussPts,
      gg);
}

MatrixDouble &
FaceElementForcesAndSourcesCore::UserDataOperator::getTangent1AtGaussPts() {
  return static_cast<FaceElementForcesAndSourcesCore *>(ptrFE)
      ->tangentOneAtGaussPts;
}

MatrixDouble &
FaceElementForcesAndSourcesCore::UserDataOperator::getTangent2AtGaussPts() {
  return static_cast<FaceElementForcesAndSourcesCore *>(ptrFE)
      ->tangentTwoAtGaussPts;
}

auto FaceElementForcesAndSourcesCore::UserDataOperator::
    getFTensor1NormalsAtGaussPts() {
  double *ptr = &*getNormalsAtGaussPts().data().begin();
  return FTensor::Tensor1<FTensor::PackPtr<double *, 3>, 3>(ptr, &ptr[1],
                                                            &ptr[2]);
}

auto FaceElementForcesAndSourcesCore::UserDataOperator::
    getFTensor1Tangent1AtGaussPts() {
  double *ptr = &*getTangent1AtGaussPts().data().begin();
  return FTensor::Tensor1<FTensor::PackPtr<double *, 3>, 3>(ptr, &ptr[1],
                                                            &ptr[2]);
}

auto FaceElementForcesAndSourcesCore::UserDataOperator::
    getFTensor1Tangent2AtGaussPts() {
  double *ptr = &*getTangent2AtGaussPts().data().begin();
  return FTensor::Tensor1<FTensor::PackPtr<double *, 3>, 3>(ptr, &ptr[1],
                                                            &ptr[2]);
}

FaceElementForcesAndSourcesCore *
FaceElementForcesAndSourcesCore::UserDataOperator::getFaceFE() {
  return static_cast<FaceElementForcesAndSourcesCore *>(ptrFE);
}

/** \deprecated use FaceElementForcesAndSourcesCore
 */
DEPRECATED typedef FaceElementForcesAndSourcesCore
    FaceElementForcesAndSourcesCoreBase;

} // namespace MoFEM

#endif //__FACEELEMENTFORCESANDSOURCESCORE_HPP__

/**
 * \defgroup mofem_forces_and_sources_tri_element Face Element
 * \brief Implementation of face element
 *
 * \ingroup mofem_forces_and_sources
 **/
