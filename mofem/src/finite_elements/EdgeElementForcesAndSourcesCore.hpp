/** \file EdgeElementForcesAndSourcesCore.hpp

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

struct FaceElementForcesAndSourcesCoreOnSideBase;

/** \brief Edge finite element
 * \ingroup mofem_forces_and_sources_edge_element
 *
 * User is implementing own operator at Gauss points level, by own object
 * derived from EdgeElementForcesAndSourcesCoreL::UserDataOperator.  Arbitrary
 * number of operator added pushing objects to rowOpPtrVector and
 * rowColOpPtrVector.
 *
 */
struct EdgeElementForcesAndSourcesCoreBase : public ForcesAndSourcesCore {

  std::string meshPositionsFieldName;

  /** \brief default operator for EDGE element
    \ingroup mofem_forces_and_sources_edge_element
    */
  struct UserDataOperator;

  enum Switches {};

  template <int SWITCH> MoFEMErrorCode opSwitch();

protected:
  EdgeElementForcesAndSourcesCoreBase(Interface &m_field);

  MatrixDouble tangentAtGaussPts;

  double lEngth;

  int numNodes;
  const EntityHandle *cOnn;
  VectorDouble dIrection;
  VectorDouble cOords;

  MoFEMErrorCode calculateEdgeDirection();
  MoFEMErrorCode setIntegrationPts();
  MoFEMErrorCode calculateCoordsAtIntegrationPts();

  friend class FaceElementForcesAndSourcesCoreOnSideBase;
};

/** \brief default operator for EDGE element
  \ingroup mofem_forces_and_sources_edge_element
  */
struct EdgeElementForcesAndSourcesCoreBase::UserDataOperator
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
  inline MatrixDouble &getTangetAtGaussPts();

  /**
   * \brief get pointer to this finite element
   */
  inline const EdgeElementForcesAndSourcesCoreBase *getEdgeFE();

  inline FTensor::Tensor1<double, 3> getFTensor1Direction();

  /**
   * \brief get get coords at gauss points

   \code
   FTensor::Index<'i',3> i;
   auto t_center;
   auto t_coords = getTensor1Coords();
   t_center(i) = 0;
   for(int nn = 0;nn!=2;nn++) {
      t_center(i) += t_coords(i);
      ++t_coords;
    }
    t_center(i) /= 2;
  \endcode

   */
  inline FTensor::Tensor1<FTensor::PackPtr<double *, 3>, 3> getTensor1Coords();

  inline FTensor::Tensor1<FTensor::PackPtr<double *, 3>, 3>
  getFTensor1TangentAtGaussPts();

  MoFEMErrorCode
  loopSideFaces(const string fe_name,
                FaceElementForcesAndSourcesCoreOnSideBase &fe_side);

protected:
  MoFEMErrorCode setPtrFE(ForcesAndSourcesCore *ptr);
};

/** \brief Edge finite element
 * \ingroup mofem_forces_and_sources_edge_element
 *
 * User is implementing own operator at Gauss points level, by own object
 * derived from EdgeElementForcesAndSourcesCoreL::UserDataOperator.  Arbitrary
 * number of operator added pushing objects to rowOpPtrVector and
 * rowColOpPtrVector.
 *
 */
template <int SWITCH>
struct EdgeElementForcesAndSourcesCoreSwitch
    : public EdgeElementForcesAndSourcesCoreBase {

  EdgeElementForcesAndSourcesCoreSwitch(Interface &m_field)
      : EdgeElementForcesAndSourcesCoreBase(m_field) {}

  MoFEMErrorCode operator()();
};

/** \brief Edge finite element default
 \ingroup mofem_forces_and_sources_edge_element

 */
using EdgeElementForcesAndSourcesCore =
    EdgeElementForcesAndSourcesCoreSwitch<0>;

template <int SWITCH>
MoFEMErrorCode EdgeElementForcesAndSourcesCoreBase::opSwitch() {
  MoFEMFunctionBegin;

  if (numeredEntFiniteElementPtr->getEntType() != MBEDGE)
    MoFEMFunctionReturnHot(0);

  CHKERR calculateEdgeDirection();
  CHKERR getSpacesAndBaseOnEntities(dataH1);
  CHKERR getEntityDataOrder<MBEDGE>(dataH1, H1);
  dataH1.dataOnEntities[MBEDGE][0].getSense() =
      1; // set sense to 1, this is this entity

  // L2
  if (dataH1.spacesOnEntities[MBEDGE].test(L2)) {
    auto &data_l2 = *dataOnElement[L2];
    CHKERR getEntityDataOrder<MBEDGE>(data_l2, L2);
    data_l2.dataOnEntities[MBEDGE][0].getSense() =
        1; // set sense to 1, this is this entity
    data_l2.spacesOnEntities[MBEDGE].set(L2);
  }

  // Hcurl
  if (dataH1.spacesOnEntities[MBEDGE].test(HCURL)) {
    auto &data_curl = *dataOnElement[HCURL];
    CHKERR getEntityDataOrder<MBEDGE>(data_curl, HCURL);
    data_curl.dataOnEntities[MBEDGE][0].getSense() =
        1; // set sense to 1, this is this entity
    data_curl.spacesOnEntities[MBEDGE].set(HCURL);
  }

  CHKERR setIntegrationPts();
  CHKERR calculateCoordsAtIntegrationPts();
  CHKERR calHierarchicalBaseFunctionsOnElement();
  CHKERR calBernsteinBezierBaseFunctionsOnElement();

  // Iterate over operators
  CHKERR loopOverOperators();

  MoFEMFunctionReturn(0);
}

template <int SWITCH>
MoFEMErrorCode EdgeElementForcesAndSourcesCoreSwitch<SWITCH>::operator()() {
  return opSwitch<SWITCH>();
}

const EntityHandle *
EdgeElementForcesAndSourcesCoreBase::UserDataOperator::getConn() {
  return static_cast<EdgeElementForcesAndSourcesCoreBase *>(ptrFE)->cOnn;
}

double EdgeElementForcesAndSourcesCoreBase::UserDataOperator::getLength() {
  return static_cast<EdgeElementForcesAndSourcesCoreBase *>(ptrFE)->lEngth;
}

double EdgeElementForcesAndSourcesCoreBase::UserDataOperator::getMeasure() {
  return getLength();
}

VectorDouble &
EdgeElementForcesAndSourcesCoreBase::UserDataOperator::getDirection() {
  return static_cast<EdgeElementForcesAndSourcesCoreBase *>(ptrFE)->dIrection;
}

VectorDouble &
EdgeElementForcesAndSourcesCoreBase::UserDataOperator::getCoords() {
  return static_cast<EdgeElementForcesAndSourcesCoreBase *>(ptrFE)->cOords;
}

MatrixDouble &
EdgeElementForcesAndSourcesCoreBase::UserDataOperator::getTangetAtGaussPts() {
  return static_cast<EdgeElementForcesAndSourcesCoreBase *>(ptrFE)
      ->tangentAtGaussPts;
}

const EdgeElementForcesAndSourcesCoreBase *
EdgeElementForcesAndSourcesCoreBase::UserDataOperator::getEdgeFE() {
  return static_cast<EdgeElementForcesAndSourcesCoreBase *>(ptrFE);
}

inline FTensor::Tensor1<double, 3>
EdgeElementForcesAndSourcesCoreBase::UserDataOperator::getFTensor1Direction() {
  return FTensor::Tensor1<double, 3>(getDirection()[0], getDirection()[1],
                                     getDirection()[2]);
}

FTensor::Tensor1<FTensor::PackPtr<double *, 3>, 3>
EdgeElementForcesAndSourcesCoreBase::UserDataOperator::getTensor1Coords() {
  double *ptr = &*getCoords().data().begin();
  return FTensor::Tensor1<FTensor::PackPtr<double *, 3>, 3>(ptr, &ptr[1],
                                                            &ptr[2]);
}

FTensor::Tensor1<FTensor::PackPtr<double *, 3>, 3>
EdgeElementForcesAndSourcesCoreBase::UserDataOperator::
    getFTensor1TangentAtGaussPts() {
  double *ptr = &*getTangetAtGaussPts().data().begin();
  return FTensor::Tensor1<FTensor::PackPtr<double *, 3>, 3>(ptr, &ptr[1],
                                                            &ptr[2]);
}

auto EdgeElementForcesAndSourcesCoreBase::UserDataOperator::getFTensor1Normal(
    const FTensor::Tensor1<double, 3> &vec) {
  FTensor::Tensor1<double, 3> t_normal;
  FTensor::Index<'i', 3> i;
  FTensor::Index<'j', 3> j;
  FTensor::Index<'k', 3> k;
  auto t_dir = getFTensor1Direction();
  t_normal(i) = FTensor::levi_civita(i, j, k) * t_dir(j) * vec(k);
  return t_normal;
}

auto EdgeElementForcesAndSourcesCoreBase::UserDataOperator::
    getFTensor1Normal() {
  FTensor::Tensor1<double, 3> t_normal{0., 0., 1.};
  return getFTensor1Normal(t_normal);
}

} // namespace MoFEM

#endif //__EDGEELEMENTFORCESANDSURCESCORE_HPP__

/**
 * \defgroup mofem_forces_and_sources_edge_element Edge Element
 *
 * \brief Implementation of edge element.
 *
 * \ingroup mofem_forces_and_sources
 */
