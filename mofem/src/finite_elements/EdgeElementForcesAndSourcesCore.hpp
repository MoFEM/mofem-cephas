/** \file EdgeElementForcesAndSourcesCore.hpp

  \brief Implementation of elements on entities.

  Those element are inherited by user to implement specific implementation of
  particular problem.

*/



#ifndef __EDGEELEMENTFORCESANDSURCESCORE_HPP__
#define __EDGEELEMENTFORCESANDSURCESCORE_HPP__

using namespace boost::numeric;

namespace MoFEM {

struct FaceElementForcesAndSourcesCoreOnSide;

/** \brief Edge finite element
 * \ingroup mofem_forces_and_sources_edge_element
 *
 * User is implementing own operator at Gauss points level, by own object
 * derived from EdgeElementForcesAndSourcesCoreL::UserDataOperator.  Arbitrary
 * number of operator added pushing objects to rowOpPtrVector and
 * rowColOpPtrVector.
 *
 */
struct EdgeElementForcesAndSourcesCore : public ForcesAndSourcesCore {
  EdgeElementForcesAndSourcesCore(Interface &m_field);

  std::string meshPositionsFieldName;

  /** \brief default operator for EDGE element
    \ingroup mofem_forces_and_sources_edge_element
    */
  struct UserDataOperator;

  enum Switches {};

  MoFEMErrorCode operator()();

  static FTensor::Tensor1<double, 3> tFaceOrientation;

protected:


  MatrixDouble tangentAtGaussPts;

  double lEngth;

  int numNodes;
  const EntityHandle *cOnn;
  VectorDouble dIrection;
  VectorDouble cOords;

  MoFEMErrorCode calculateEdgeDirection();
  MoFEMErrorCode setIntegrationPts();
  MoFEMErrorCode calculateCoordsAtIntegrationPts();

  friend class FaceElementForcesAndSourcesCoreOnSide;
};

/** \brief default operator for EDGE element
  \ingroup mofem_forces_and_sources_edge_element
  */
struct EdgeElementForcesAndSourcesCore::UserDataOperator
    : public ForcesAndSourcesCore::UserDataOperator {

  using ForcesAndSourcesCore::UserDataOperator::UserDataOperator;

  /** \brief get element connectivity
   */
  inline const EntityHandle *getConn();

  /**
   * \brief get edge length
   */
  inline double getLength();

  /**
   * \brief get measure of element
   * @return length of face
   */
  inline double getMeasure();

  /**
   * \brief get edge direction
   */
  inline VectorDouble &getDirection();

  /**
   * \brief get edge normal
   * NOTE: it should be used only in 2D analysis
   */
  inline auto getFTensor1Normal();

  /**
   * @brief get ftensor1 edge normal
   *
   * @param vec vector in third direction
   * @return auto
   */
  inline auto getFTensor1Normal(const FTensor::Tensor1<double, 3> &vec);

  /**
   * \brief get edge node coordinates
   */
  inline VectorDouble &getCoords();

  /**
   * \brief get tangent vector to edge curve at integration points
   */
  inline MatrixDouble &getTangentAtGaussPts();

  DEPRECATED inline MatrixDouble &getTangetAtGaussPts() {
    return getTangentAtGaussPts();
  }

  /**
   * \brief get pointer to this finite element
   */
  inline const EdgeElementForcesAndSourcesCore *getEdgeFE();

  inline FTensor::Tensor1<double, 3> getFTensor1Direction();

  /**
   * \brief get get coords at gauss points

   \code
   FTensor::Index<'i',3> i;
   auto t_center;
   auto t_coords = getFTensor1Coords();
   t_center(i) = 0;
   for(int nn = 0;nn!=2;nn++) {
      t_center(i) += t_coords(i);
      ++t_coords;
    }
    t_center(i) /= 2;
  \endcode

   */
  inline FTensor::Tensor1<FTensor::PackPtr<double *, 3>, 3> getFTensor1Coords();

  /**
   * @deprecated use getFTensor1Coords
   */
  DEPRECATED FTensor::Tensor1<FTensor::PackPtr<double *, 3>, 3>
  getTensor1Coords() {
    return getFTensor1Coords();
  }

  /**
   * @brief Get tangent vector at integration points
   *
   * @return FTensor::Tensor1<FTensor::PackPtr<double *, 3>, DIM>
   */
  template <int DIM = 3>
  inline FTensor::Tensor1<FTensor::PackPtr<double *, 3>, DIM>
  getFTensor1TangentAtGaussPts();

  MoFEMErrorCode
  loopSideFaces(const string fe_name,
                FaceElementForcesAndSourcesCoreOnSide &fe_side);

protected:
  MoFEMErrorCode setPtrFE(ForcesAndSourcesCore *ptr);
};

const EntityHandle *
EdgeElementForcesAndSourcesCore::UserDataOperator::getConn() {
  return static_cast<EdgeElementForcesAndSourcesCore *>(ptrFE)->cOnn;
}

double EdgeElementForcesAndSourcesCore::UserDataOperator::getLength() {
  return static_cast<EdgeElementForcesAndSourcesCore *>(ptrFE)->lEngth;
}

double EdgeElementForcesAndSourcesCore::UserDataOperator::getMeasure() {
  return getLength();
}

VectorDouble &
EdgeElementForcesAndSourcesCore::UserDataOperator::getDirection() {
  return static_cast<EdgeElementForcesAndSourcesCore *>(ptrFE)->dIrection;
}

VectorDouble &
EdgeElementForcesAndSourcesCore::UserDataOperator::getCoords() {
  return static_cast<EdgeElementForcesAndSourcesCore *>(ptrFE)->cOords;
}

MatrixDouble &
EdgeElementForcesAndSourcesCore::UserDataOperator::getTangentAtGaussPts() {
  return static_cast<EdgeElementForcesAndSourcesCore *>(ptrFE)
      ->tangentAtGaussPts;
}

const EdgeElementForcesAndSourcesCore *
EdgeElementForcesAndSourcesCore::UserDataOperator::getEdgeFE() {
  return static_cast<EdgeElementForcesAndSourcesCore *>(ptrFE);
}

inline FTensor::Tensor1<double, 3>
EdgeElementForcesAndSourcesCore::UserDataOperator::getFTensor1Direction() {
  return FTensor::Tensor1<double, 3>(getDirection()[0], getDirection()[1],
                                     getDirection()[2]);
}

FTensor::Tensor1<FTensor::PackPtr<double *, 3>, 3>
EdgeElementForcesAndSourcesCore::UserDataOperator::getFTensor1Coords() {
  double *ptr = &*getCoords().data().begin();
  return FTensor::Tensor1<FTensor::PackPtr<double *, 3>, 3>(ptr, &ptr[1],
                                                            &ptr[2]);
}

template <>
inline FTensor::Tensor1<FTensor::PackPtr<double *, 3>, 3>
EdgeElementForcesAndSourcesCore::UserDataOperator::getFTensor1TangentAtGaussPts<
    3>() {
  double *ptr = &*getTangentAtGaussPts().data().begin();
  return FTensor::Tensor1<FTensor::PackPtr<double *, 3>, 3>(ptr, &ptr[1],
                                                            &ptr[2]);
}

template <>
inline FTensor::Tensor1<FTensor::PackPtr<double *, 3>, 2>
EdgeElementForcesAndSourcesCore::UserDataOperator::getFTensor1TangentAtGaussPts<
    2>() {
  double *ptr = &*getTangentAtGaussPts().data().begin();
  return FTensor::Tensor1<FTensor::PackPtr<double *, 3>, 2>(ptr, &ptr[1]);
}

auto EdgeElementForcesAndSourcesCore::UserDataOperator::getFTensor1Normal(
    const FTensor::Tensor1<double, 3> &vec) {
  FTensor::Tensor1<double, 3> t_normal;
  FTensor::Index<'i', 3> i;
  FTensor::Index<'j', 3> j;
  FTensor::Index<'k', 3> k;
  auto t_dir = getFTensor1Direction();
  t_normal(i) = FTensor::levi_civita(i, j, k) * t_dir(j) * vec(k);
  return t_normal;
}

auto EdgeElementForcesAndSourcesCore::UserDataOperator::
    getFTensor1Normal() {
  return getFTensor1Normal(tFaceOrientation);
}

/**
 * @deprecated use EdgeElementForcesAndSourcesCore
 */
DEPRECATED typedef EdgeElementForcesAndSourcesCore
    EdgeElementForcesAndSourcesCoreBase;

} // namespace MoFEM

#endif //__EDGEELEMENTFORCESANDSURCESCORE_HPP__

/**
 * \defgroup mofem_forces_and_sources_edge_element Edge Element
 *
 * \brief Implementation of edge element.
 *
 * \ingroup mofem_forces_and_sources
 */
