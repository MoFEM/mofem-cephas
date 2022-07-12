/** \file VolumeElementForcesAndSourcesCoreOnSide.hpp
  \brief Volume element.

  Those element are inherited by user to implement specific implementation of
  particular problem.

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

#ifndef __VOLUMEELEMENTFORCESANDSOURCESCORE_ONSIDE_HPP__
#define __VOLUMEELEMENTFORCESANDSOURCESCORE_ONSIDE_HPP__

using namespace boost::numeric;

namespace MoFEM {

/**
 * \brief Base volume element used to integrate on skeleton
 * \ingroup mofem_forces_and_sources_volume_element
 */
struct VolumeElementForcesAndSourcesCoreOnSide
    : public VolumeElementForcesAndSourcesCore {

  using VolumeElementForcesAndSourcesCore::VolumeElementForcesAndSourcesCore;

  /**
   * @brief Get the face nodes mapped on volume element
   *
   * \todo That this is not general, e.g., for quad number of nodes is 4.
   *
   * @return const std::array<int, 4>&
   */
  inline const std::array<int, 4> &getFaceConnMap() const;

  /**
   * @brief Get face nodes maped on volume
   *
   * \todo That this is not general, e.g., for prism or hex, size of fixed array
   * is wrong.
   *
   * @return const sdt::array<int, 4>&
   */
  inline const std::array<int, 8> &getTetConnMap() const;

  /**
   * @brief Get node on volume opposite to volume element
   *
   * \todo That this is not general, e.g., for prism or hex, opoosite node is
   * not unique.
   *
   * @return int
   */
  inline int getOppositeNode() const;

  /**
   * @brief Sense face on volume
   *
   * @return int
   */
  inline int getFaceSense() const;

  /**
   * @brief Face number on the volume
   *
   * @return int
   */
  inline int getFaceSideNumber() const;

  /** \brief default operator for TET element
   * \ingroup mofem_forces_and_sources_volume_element
   */
  struct UserDataOperator;

  int getRule(int order);
  MoFEMErrorCode setGaussPts(int order);

private:
  int faceSense;      ///< Sense of face, could be 1 or -1
  int faceSideNumber; ///< Face side number
  std::array<int, 4> faceConnMap;
  std::array<int, 8> tetConnMap;
  int oppositeNode;
};

/** \brief default operator for TET element
 * \ingroup mofem_forces_and_sources_volume_element
 */
struct VolumeElementForcesAndSourcesCoreOnSide::UserDataOperator
    : public VolumeElementForcesAndSourcesCore::UserDataOperator {

  using VolumeElementForcesAndSourcesCore::UserDataOperator::UserDataOperator;

  inline VolumeElementForcesAndSourcesCoreOnSide *getVolumeFE() const;

  inline FaceElementForcesAndSourcesCore *getFaceFE() const;

  /**
   * \brief get face sense in respect to volume
   * @return error code
   */
  inline int getFaceSense() const;

  /**
   * \brief get face side number in respect to volume
   * @return error code
   */
  inline int getFaceSideNumber() const;

  inline bool getEdgeFace(const int ee) const;

  /**
   * get face normal on side which is this element
   * @return face normal
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

  /** \brief get face coordinates at Gauss pts.

  \note Coordinates should be the same what function getCoordsAtGaussPts
  on tets is returning. If both coordinates are different it is error, or you
  do something very unusual.

   */
  inline MatrixDouble &getFaceCoordsAtGaussPts();

protected:
  MoFEMErrorCode setPtrFE(ForcesAndSourcesCore *ptr);
};

const std::array<int, 4> &
VolumeElementForcesAndSourcesCoreOnSide::getFaceConnMap() const {
  return faceConnMap;
}

const std::array<int, 8> &
VolumeElementForcesAndSourcesCoreOnSide::getTetConnMap() const {
  return tetConnMap;
}

int VolumeElementForcesAndSourcesCoreOnSide::getOppositeNode() const {
  return oppositeNode;
}

int VolumeElementForcesAndSourcesCoreOnSide::getFaceSense() const {
  return faceSense;
}

int VolumeElementForcesAndSourcesCoreOnSide::getFaceSideNumber() const {
  return faceSideNumber;
}

FaceElementForcesAndSourcesCore *
VolumeElementForcesAndSourcesCoreOnSide::UserDataOperator::getFaceFE()
    const {
  return static_cast<FaceElementForcesAndSourcesCore *>(
      getVolumeFE()->sidePtrFE);
}

VectorDouble &
VolumeElementForcesAndSourcesCoreOnSide::UserDataOperator::getNormal() {
  return getFaceFE()->nOrmal;
}

VectorDouble &
VolumeElementForcesAndSourcesCoreOnSide::UserDataOperator::getTangent1() {
  return getFaceFE()->tangentOne;
}

VectorDouble &
VolumeElementForcesAndSourcesCoreOnSide::UserDataOperator::getTangent2() {
  return getFaceFE()->tangentTwo;
}

VolumeElementForcesAndSourcesCoreOnSide *
VolumeElementForcesAndSourcesCoreOnSide::UserDataOperator::getVolumeFE()
    const {
  return static_cast<VolumeElementForcesAndSourcesCoreOnSide *>(ptrFE);
}

int VolumeElementForcesAndSourcesCoreOnSide::UserDataOperator::
    getFaceSense() const {
  return getVolumeFE()->faceSense;
}

auto VolumeElementForcesAndSourcesCoreOnSide::UserDataOperator::
    getFTensor1Normal() {
  double *ptr = &*getNormal().data().begin();
  return FTensor::Tensor1<double *, 3>(ptr, &ptr[1], &ptr[2]);
}

auto VolumeElementForcesAndSourcesCoreOnSide::UserDataOperator::
    getFTensor1Tangent1() {
  double *ptr = &*getTangent1().data().begin();
  return FTensor::Tensor1<double *, 3>(ptr, &ptr[1], &ptr[2]);
}

auto VolumeElementForcesAndSourcesCoreOnSide::UserDataOperator::
    getFTensor1Tangent2() {
  double *ptr = &*getTangent2().data().begin();
  return FTensor::Tensor1<double *, 3>(ptr, &ptr[1], &ptr[2]);
}

inline auto VolumeElementForcesAndSourcesCoreOnSide::UserDataOperator::
    getFTensor1NormalsAtGaussPts() {
  double *ptr = &*getNormalsAtGaussPts().data().begin();
  return FTensor::Tensor1<FTensor::PackPtr<double *, 3>, 3>(ptr, &ptr[1],
                                                            &ptr[2]);
}

int VolumeElementForcesAndSourcesCoreOnSide::UserDataOperator::
    getFaceSideNumber() const {
  return getVolumeFE()->faceSideNumber;
}

MatrixDouble &VolumeElementForcesAndSourcesCoreOnSide::UserDataOperator::
    getNormalsAtGaussPts() {
  return getFaceFE()->normalsAtGaussPts;
}

MatrixDouble &VolumeElementForcesAndSourcesCoreOnSide::UserDataOperator::
    getFaceCoordsAtGaussPts() {
  return getFaceFE()->coordsAtGaussPts;
}

bool VolumeElementForcesAndSourcesCoreOnSide::UserDataOperator::getEdgeFace(
    const int ee) const {
  constexpr bool edges_on_faces[6][4] = {{true, false, false, true}, // e0
                                         {false, true, false, true}, // e1
                                         {false, false, true, true}, // e2
                                         {true, false, true, false}, // e3
                                         {true, true, false, false}, // e4
                                         {false, true, true, false}};
  return edges_on_faces[ee][getFaceSideNumber()];
}

ublas::matrix_row<MatrixDouble> VolumeElementForcesAndSourcesCoreOnSide::
    UserDataOperator::getNormalsAtGaussPts(const int gg) {
  return ublas::matrix_row<MatrixDouble>(getNormalsAtGaussPts(), gg);
}

/**
 * @deprecated use VolumeElementForcesAndSourcesCore
 * 
 */
DEPRECATED typedef VolumeElementForcesAndSourcesCoreOnSide
    VolumeElementForcesAndSourcesCoreOnSideBase;

/**
 * @deprecated do not use this template, it is obsolete
 */
template <int SWITCH>
struct VolumeElementForcesAndSourcesCoreOnSideSwitch
    : public VolumeElementForcesAndSourcesCoreOnSide {
  using VolumeElementForcesAndSourcesCoreOnSide::
      VolumeElementForcesAndSourcesCoreOnSide;
  using UserDataOperator =
      VolumeElementForcesAndSourcesCoreOnSide::UserDataOperator;
};



} // namespace MoFEM

#endif //__VOLUMEELEMENTFORCESANDSOURCESCORE_ONSIDE_HPP__
