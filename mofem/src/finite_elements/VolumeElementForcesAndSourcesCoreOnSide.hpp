/** \file VolumeElementForcesAndSourcesCoreOnSide.hpp
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

#ifndef __VOLUMEELEMENTFORCESANDSOURCESCORE_ONSIDE_HPP__
#define __VOLUMEELEMENTFORCESANDSOURCESCORE_ONSIDE_HPP__

using namespace boost::numeric;

namespace MoFEM {

/**
 * \brief Base volume element used to integrate on skeleton
 * \ingroup mofem_forces_and_sources_volume_element
 */
struct VolumeElementForcesAndSourcesCoreOnSideBase
    : public VolumeElementForcesAndSourcesCore {

  using VolumeElementForcesAndSourcesCore::VolumeElementForcesAndSourcesCore;

  /**
   * @brief Get the face nodes mapped on volume element
   *
   * \todo That this is not general, e.g., for quad number of nodes is 4.
   *
   * @return const std::array<int, 3>&
   */
  inline const std::array<int, 3> &getFaceConnMap() const;

  /**
   * @brief Get face nodes maped on volume
   *
   * \todo That this is not general, e.g., for prism or hex, size of fixed array
   * is wrong.
   *
   * @return const sdt::array<int, 4>&
   */
  inline const std::array<int, 4> &getTetConnMap() const;

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
  struct UserDataOperator
      : public VolumeElementForcesAndSourcesCore::UserDataOperator {

    using VolumeElementForcesAndSourcesCore::UserDataOperator::UserDataOperator;

    inline VolumeElementForcesAndSourcesCoreOnSideBase *getVolumeFE() const;

    inline FaceElementForcesAndSourcesCoreBase *getFaceFE() const;

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

  int getRule(int order);
  MoFEMErrorCode setGaussPts(int order);

private:
  int faceSense;      ///< Sense of face, could be 1 or -1
  int faceSideNumber; ///< Face side number
  std::array<int, 3> faceConnMap;
  std::array<int, 4> tetConnMap;
  int oppositeNode;
};

/**
 * @brief Volume side finite element with switches
 *
 * Using SWITCH to off functions
 *
 * @tparam SWITCH
 */
template <int SWITCH>
struct VolumeElementForcesAndSourcesCoreOnSideSwitch
    : public VolumeElementForcesAndSourcesCoreOnSideBase {

  using VolumeElementForcesAndSourcesCoreOnSideBase::
      VolumeElementForcesAndSourcesCoreOnSideBase;

  using UserDataOperator =
      VolumeElementForcesAndSourcesCoreOnSideBase::UserDataOperator;

  MoFEMErrorCode operator()();
};

/** \brief Volume element used to integrate on skeleton
 \ingroup mofem_forces_and_sources_volume_element

 */
using VolumeElementForcesAndSourcesCoreOnSide =
    VolumeElementForcesAndSourcesCoreOnSideSwitch<0>;

const std::array<int, 3> &
VolumeElementForcesAndSourcesCoreOnSideBase::getFaceConnMap() const {
  return faceConnMap;
}

const std::array<int, 4> &
VolumeElementForcesAndSourcesCoreOnSideBase::getTetConnMap() const {
  return tetConnMap;
}

int VolumeElementForcesAndSourcesCoreOnSideBase::getOppositeNode() const {
  return oppositeNode;
}

int VolumeElementForcesAndSourcesCoreOnSideBase::getFaceSense() const {
  return faceSense;
}

int VolumeElementForcesAndSourcesCoreOnSideBase::getFaceSideNumber() const {
  return faceSideNumber;
}

FaceElementForcesAndSourcesCoreBase *
VolumeElementForcesAndSourcesCoreOnSideBase::UserDataOperator::getFaceFE()
    const {
  return static_cast<FaceElementForcesAndSourcesCoreBase *>(
      getVolumeFE()->sidePtrFE);
}

VectorDouble &
VolumeElementForcesAndSourcesCoreOnSideBase::UserDataOperator::getNormal() {
  return getFaceFE()->nOrmal;
}

VectorDouble &
VolumeElementForcesAndSourcesCoreOnSideBase::UserDataOperator::getTangent1() {
  return getFaceFE()->tangentOne;
}

VectorDouble &
VolumeElementForcesAndSourcesCoreOnSideBase::UserDataOperator::getTangent2() {
  return getFaceFE()->tangentTwo;
}

VolumeElementForcesAndSourcesCoreOnSideBase *
VolumeElementForcesAndSourcesCoreOnSideBase::UserDataOperator::getVolumeFE()
    const {
  return static_cast<VolumeElementForcesAndSourcesCoreOnSideBase *>(ptrFE);
}

int VolumeElementForcesAndSourcesCoreOnSideBase::UserDataOperator::
    getFaceSense() const {
  return getVolumeFE()->faceSense;
}

auto VolumeElementForcesAndSourcesCoreOnSideBase::UserDataOperator::
    getFTensor1Normal() {
  double *ptr = &*getNormal().data().begin();
  return FTensor::Tensor1<double *, 3>(ptr, &ptr[1], &ptr[2]);
}

auto VolumeElementForcesAndSourcesCoreOnSideBase::UserDataOperator::
    getFTensor1Tangent1() {
  double *ptr = &*getTangent1().data().begin();
  return FTensor::Tensor1<double *, 3>(ptr, &ptr[1], &ptr[2]);
}

auto VolumeElementForcesAndSourcesCoreOnSideBase::UserDataOperator::
    getFTensor1Tangent2() {
  double *ptr = &*getTangent2().data().begin();
  return FTensor::Tensor1<double *, 3>(ptr, &ptr[1], &ptr[2]);
}

inline auto VolumeElementForcesAndSourcesCoreOnSideBase::UserDataOperator::
    getFTensor1NormalsAtGaussPts() {
  double *ptr = &*getNormalsAtGaussPts().data().begin();
  return FTensor::Tensor1<FTensor::PackPtr<double *, 3>, 3>(ptr, &ptr[1],
                                                            &ptr[2]);
}

int VolumeElementForcesAndSourcesCoreOnSideBase::UserDataOperator::
    getFaceSideNumber() const {
  return getVolumeFE()->faceSideNumber;
}

MatrixDouble &VolumeElementForcesAndSourcesCoreOnSideBase::UserDataOperator::
    getNormalsAtGaussPts() {
  return getFaceFE()->normalsAtGaussPts;
}

MatrixDouble &VolumeElementForcesAndSourcesCoreOnSideBase::UserDataOperator::
    getFaceCoordsAtGaussPts() {
  return getFaceFE()->coordsAtGaussPts;
}

bool VolumeElementForcesAndSourcesCoreOnSideBase::UserDataOperator::getEdgeFace(
    const int ee) const {
  constexpr bool edges_on_faces[6][4] = {{true, false, false, true}, // e0
                                         {false, true, false, true}, // e1
                                         {false, false, true, true}, // e2
                                         {true, false, true, false}, // e3
                                         {true, true, false, false}, // e4
                                         {false, true, true, false}};
  return edges_on_faces[ee][getFaceSideNumber()];
}

ublas::matrix_row<MatrixDouble> VolumeElementForcesAndSourcesCoreOnSideBase::
    UserDataOperator::getNormalsAtGaussPts(const int gg) {
  return ublas::matrix_row<MatrixDouble>(getNormalsAtGaussPts(), gg);
}

template <int SWITCH>
MoFEMErrorCode VolumeElementForcesAndSourcesCoreOnSideSwitch<SWITCH>::
operator()() {
  return opSwitch<SWITCH>();
}

} // namespace MoFEM

#endif //__VOLUMEELEMENTFORCESANDSOURCESCORE_ONSIDE_HPP__
