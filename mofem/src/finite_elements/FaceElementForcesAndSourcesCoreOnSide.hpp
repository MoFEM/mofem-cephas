/** \file FaceElementForcesAndSourcesCoreOnSide.hpp
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

#ifndef __FACEELEMENTFORCESANDSOURCESCORE_ONSIDE__HPP__
#define __FACEELEMENTFORCESANDSOURCESCORE_ONSIDE__HPP__

using namespace boost::numeric;

namespace MoFEM {

/**
 * \brief Base face element used to integrate on skeleton
 * \ingroup mofem_forces_and_sources_volume_element
 */
struct FaceElementForcesAndSourcesCoreOnSideBase
    : public FaceElementForcesAndSourcesCoreBase {

  using FaceElementForcesAndSourcesCoreBase::
      FaceElementForcesAndSourcesCoreBase;

  int getRule(int order);

  /**
   * @brief Get the face nodes mapped on volume element
   *
   * \todo That this is not general, e.g., for quad number of nodes is 4.
   *
   * @return const std::array<int, 3>&
   */
  inline const std::array<int, 2> &getEdgeConnMap() const;

  /**
   * @brief Get face nodes maped on volume
   *
   * @return const sdt::array<int, 3>&
   */
  inline const std::array<int, 3> &getFaceConnMap() const;

  /**
   * @brief Get node on volume opposite to volume element
   *
   * @return int
   */
  inline int getOppositeNode() const;

  /**
   * @brief Sense face on volume
   *
   * @return int
   */
  inline int getEdgeSense() const;

  /**
   * @brief Face number on the volume
   *
   * @return int
   */
  inline int getEdgeSideNumber() const;

  /** \brief default operator for Face element
   * \ingroup mofem_forces_and_sources_volume_element
   */
  struct UserDataOperator
      : public FaceElementForcesAndSourcesCoreBase::UserDataOperator {

    using FaceElementForcesAndSourcesCoreBase::UserDataOperator::
        UserDataOperator;

    /** \brief return pointer to Generic Volume Finite Element object
     */
    inline const FaceElementForcesAndSourcesCoreOnSideBase *getFaceFE() const;

    /**
     * @brief Get the edge side finite element
     *
     * @return EdgeElementForcesAndSourcesCoreBase*
     */
    inline EdgeElementForcesAndSourcesCoreBase *getEdgeFE() const;

    /**
     * \brief get face sense in respect to volume
     * @return error code
     */
    inline int getEdgeSense() const;

    /**
     * \brief get face side number in respect to volume
     * @return error code
     */
    inline int getEdgeSideNumber() const;

    /**
     * get face normal on side which is this element
     * @return face normal
     */
    inline VectorDouble &getDirection();

    /** \brief get normal as tensor
     */
    inline auto getFTensor1Direction();

    /** \brief get face coordinates at Gauss pts.

    \note Coordinates should be the same what function getCoordsAtGaussPts
    on face is returning. If both coordinates are different it is error, or you
    do something very unusual.

     */
    inline MatrixDouble &getEdgeCoordsAtGaussPts();

  protected:

    MoFEMErrorCode setPtrFE(ForcesAndSourcesCore *ptr);

  };

protected:

  MoFEMErrorCode setGaussPts(int order);

private:
  int edgeSense;      ///< Sense of edge, could be 1 or -1
  int edgeSideNumber; ///< Edge side number
  std::array<int, 2> edgeConnMap;
  std::array<int, 3> faceConnMap;
  int oppositeNode;
};

/**
 * @brief Face side finite element with switches
 *
 * Using SWITCH to off functions
 *
 * @tparam SWITCH
 */
template <int SWITCH>
struct FaceElementForcesAndSourcesCoreOnSideSwitch
    : public FaceElementForcesAndSourcesCoreOnSideBase {

  FaceElementForcesAndSourcesCoreOnSideSwitch(MoFEM::Interface &m_field)
      : FaceElementForcesAndSourcesCoreOnSideBase(m_field) {}

  using UserDataOperator =
      FaceElementForcesAndSourcesCoreOnSideBase::UserDataOperator;

  MoFEMErrorCode operator()();
};

const std::array<int, 2> &
FaceElementForcesAndSourcesCoreOnSideBase::getEdgeConnMap() const {
  return edgeConnMap;
}

const std::array<int, 3> &
FaceElementForcesAndSourcesCoreOnSideBase::getFaceConnMap() const {
  return faceConnMap;
}

int FaceElementForcesAndSourcesCoreOnSideBase::getOppositeNode() const {
  return oppositeNode;
}

int FaceElementForcesAndSourcesCoreOnSideBase::getEdgeSense() const {
  return edgeSense;
}

int FaceElementForcesAndSourcesCoreOnSideBase::getEdgeSideNumber() const {
  return edgeSideNumber;
}

const FaceElementForcesAndSourcesCoreOnSideBase *
FaceElementForcesAndSourcesCoreOnSideBase::UserDataOperator::getFaceFE() const {
  return static_cast<FaceElementForcesAndSourcesCoreOnSideBase *>(ptrFE);
}

EdgeElementForcesAndSourcesCoreBase *
FaceElementForcesAndSourcesCoreOnSideBase::UserDataOperator::getEdgeFE() const {
  return static_cast<EdgeElementForcesAndSourcesCoreBase *>(ptrFE->sidePtrFE);
}

int FaceElementForcesAndSourcesCoreOnSideBase::UserDataOperator::getEdgeSense()
    const {
  return getFaceFE()->edgeSense;
}

int FaceElementForcesAndSourcesCoreOnSideBase::UserDataOperator::
    getEdgeSideNumber() const {
  return getFaceFE()->edgeSideNumber;
}

VectorDouble &
FaceElementForcesAndSourcesCoreOnSideBase::UserDataOperator::getDirection() {
  return getEdgeFE()->dIrection;
}

auto FaceElementForcesAndSourcesCoreOnSideBase::UserDataOperator::
    getFTensor1Direction() {
  double *ptr = &*getDirection().data().begin();
  return FTensor::Tensor1<double *, 3>(ptr, &ptr[1], &ptr[2]);
}

MatrixDouble &FaceElementForcesAndSourcesCoreOnSideBase::UserDataOperator::
    getEdgeCoordsAtGaussPts() {
  return getEdgeFE()->coordsAtGaussPts;
}

template <int SWITCH>
MoFEMErrorCode FaceElementForcesAndSourcesCoreOnSideSwitch<SWITCH>::
operator()() {
  return OpSwitch<SWITCH>();
}

} // namespace MoFEM

#endif //__FACEELEMENTFORCESANDSOURCESCORE_ONSIDE___HPP__
