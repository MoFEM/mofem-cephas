/** \file FaceElementForcesAndSourcesCoreOnSide.hpp
  \brief Implementation of face element.

*/

/* MIT License
 *
 * Copyright (c) 2022
 * 
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 * 
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 * 
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */

#ifndef __FACEELEMENTFORCESANDSOURCESCORE_ONSIDE__HPP__
#define __FACEELEMENTFORCESANDSOURCESCORE_ONSIDE__HPP__

using namespace boost::numeric;

namespace MoFEM {

/**
 * \brief Base face element used to integrate on skeleton
 * \ingroup mofem_forces_and_sources_volume_element
 */
struct FaceElementForcesAndSourcesCoreOnSide
    : public FaceElementForcesAndSourcesCore {

  using FaceElementForcesAndSourcesCore::
      FaceElementForcesAndSourcesCore;

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
   * @return const sdt::array<int, 4>&
   */
  inline const std::array<int, 4> &getFaceConnMap() const;

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
  struct UserDataOperator;

protected:
  MoFEMErrorCode setGaussPts(int order);

private:
  int edgeSense;      ///< Sense of edge, could be 1 or -1
  int edgeSideNumber; ///< Edge side number
  std::array<int, 2> edgeConnMap;
  std::array<int, 4> faceConnMap;
  int oppositeNode;
};

/** \brief default operator for Face element
 * \ingroup mofem_forces_and_sources_volume_element
 */
struct FaceElementForcesAndSourcesCoreOnSide::UserDataOperator
    : public FaceElementForcesAndSourcesCore::UserDataOperator {

  using FaceElementForcesAndSourcesCore::UserDataOperator::UserDataOperator;

  /** \brief return pointer to Generic Volume Finite Element object
   */
  inline const FaceElementForcesAndSourcesCoreOnSide *getFaceFE() const;

  /**
   * @brief Get the edge side finite element
   *
   * @return EdgeElementForcesAndSourcesCore*
   */
  inline EdgeElementForcesAndSourcesCore *getEdgeFE() const;

  /**
   * \brief get face sense in respect to volume
   * @return edge sense
   */
  inline int getEdgeSense() const;

  /**
   * \brief get face side number in respect to volume
   * @return edge side number
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

const std::array<int, 2> &
FaceElementForcesAndSourcesCoreOnSide::getEdgeConnMap() const {
  return edgeConnMap;
}

const std::array<int, 4> &
FaceElementForcesAndSourcesCoreOnSide::getFaceConnMap() const {
  return faceConnMap;
}

int FaceElementForcesAndSourcesCoreOnSide::getOppositeNode() const {
  return oppositeNode;
}

int FaceElementForcesAndSourcesCoreOnSide::getEdgeSense() const {
  return edgeSense;
}

int FaceElementForcesAndSourcesCoreOnSide::getEdgeSideNumber() const {
  return edgeSideNumber;
}

const FaceElementForcesAndSourcesCoreOnSide *
FaceElementForcesAndSourcesCoreOnSide::UserDataOperator::getFaceFE() const {
  return static_cast<FaceElementForcesAndSourcesCoreOnSide *>(ptrFE);
}

EdgeElementForcesAndSourcesCore *
FaceElementForcesAndSourcesCoreOnSide::UserDataOperator::getEdgeFE() const {
  return static_cast<EdgeElementForcesAndSourcesCore *>(ptrFE->sidePtrFE);
}

int FaceElementForcesAndSourcesCoreOnSide::UserDataOperator::getEdgeSense()
    const {
  return getFaceFE()->edgeSense;
}

int FaceElementForcesAndSourcesCoreOnSide::UserDataOperator::
    getEdgeSideNumber() const {
  return getFaceFE()->edgeSideNumber;
}

VectorDouble &
FaceElementForcesAndSourcesCoreOnSide::UserDataOperator::getDirection() {
  return getEdgeFE()->dIrection;
}

auto FaceElementForcesAndSourcesCoreOnSide::UserDataOperator::
    getFTensor1Direction() {
  double *ptr = &*getDirection().data().begin();
  return FTensor::Tensor1<double *, 3>(ptr, &ptr[1], &ptr[2]);
}

MatrixDouble &FaceElementForcesAndSourcesCoreOnSide::UserDataOperator::
    getEdgeCoordsAtGaussPts() {
  return getEdgeFE()->coordsAtGaussPts;
}

/**
 * @deprecated do not use needed for back compatibility
 */
template <int SWITCH>
struct FaceElementForcesAndSourcesCoreOnSideSwitch
    : public FaceElementForcesAndSourcesCoreOnSide {
  using FaceElementForcesAndSourcesCoreOnSide::
      FaceElementForcesAndSourcesCoreOnSide;
  using UserDataOperator =
      FaceElementForcesAndSourcesCoreOnSide::UserDataOperator;
};

} // namespace MoFEM

#endif //__FACEELEMENTFORCESANDSOURCESCORE_ONSIDE___HPP__
