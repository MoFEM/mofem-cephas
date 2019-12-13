/** \file FaceElementForcesAndSourcesCoreOnVolumeSide.hpp
  \brief Face element.

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

#ifndef __FACEELEMENTFORCESANDSOURCESCORE_ONVOLUMESIDE_HPP__
#define __FACEELEMENTFORCESANDSOURCESCORE_ONVOLUMESIDE_HPP__

using namespace boost::numeric;

namespace MoFEM {

/**
 * \brief Base volume element used to integrate on contact surface (could be
 * extended to other volume elements) \ingroup
 * mofem_forces_and_sources_volume_element
 */
struct FaceElementForcesAndSourcesCoreOnVolumeSideBase
    : public FaceElementForcesAndSourcesCore {

  using FaceElementForcesAndSourcesCore::FaceElementForcesAndSourcesCore;

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
   * \todo That this is not general, e.g., for prism or hex, opposite node is
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
      : public FaceElementForcesAndSourcesCore::UserDataOperator {

    using FaceElementForcesAndSourcesCore::UserDataOperator::UserDataOperator;

    inline FaceElementForcesAndSourcesCoreOnVolumeSideBase *
    getFaceFE() const;

    inline ContactPrismElementForcesAndSourcesCore *getContactFE() const;

    /** \brief get face coordinates at Gauss pts.

    \note Coordinates should be the same what function getMasterCoordsAtGaussPts
    on tets is returning. If both coordinates are different it is error, or you
    do something very unusual.

     */
    inline MatrixDouble &getMasterCoordsAtGaussPts();

    /** \brief get face coordinates at Gauss pts.

    \note Coordinates should be the same what function getSlaveCoordsAtGaussPts
    on tets is returning. If both coordinates are different it is error, or you
    do something very unusual.

    */
    inline MatrixDouble &getSlaveCoordsAtGaussPts();

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
  };

  int getRule(int order);
  MoFEMErrorCode setGaussPts(int order);

private:
  int faceSense;      ///< Sense of face, could be 1 or -1
  int faceSideNumber; ///< Face side number
  std::array<int, 3> faceConnMap;
};

/**
 * @brief Volume side finite element with switches
 *
 * Using SWITCH to off functions
 *
 * @tparam SWITCH
 */
template <int SWITCH>
struct FaceElementForcesAndSourcesCoreOnVolumeSideSwitch
    : public FaceElementForcesAndSourcesCoreOnVolumeSideBase {

  using FaceElementForcesAndSourcesCoreOnVolumeSideBase::
      FaceElementForcesAndSourcesCoreOnVolumeSideBase;

  using UserDataOperator =
      FaceElementForcesAndSourcesCoreOnVolumeSideBase::UserDataOperator;

  MoFEMErrorCode operator()();
};

/** \brief Volume element used to integrate on contact element (could be
 extended for other volume elements) \ingroup
 mofem_forces_and_sources_volume_element

 */
using FaceElementForcesAndSourcesCoreOnVolumeSide =
    FaceElementForcesAndSourcesCoreOnVolumeSideSwitch<0>;

const std::array<int, 3> &
FaceElementForcesAndSourcesCoreOnVolumeSideBase::getFaceConnMap() const {
  return faceConnMap;
}


int FaceElementForcesAndSourcesCoreOnVolumeSideBase::getFaceSense() const {
  return faceSense;
}

int FaceElementForcesAndSourcesCoreOnVolumeSideBase::getFaceSideNumber()
    const {
  return faceSideNumber;
}

FaceElementForcesAndSourcesCoreOnVolumeSideBase *
FaceElementForcesAndSourcesCoreOnVolumeSideBase::UserDataOperator::
    getFaceFE() const {
  return static_cast<FaceElementForcesAndSourcesCoreOnVolumeSideBase *>(
      ptrFE);
}

ContactPrismElementForcesAndSourcesCore *
FaceElementForcesAndSourcesCoreOnVolumeSideBase::UserDataOperator::
    getContactFE() const {
  return static_cast<ContactPrismElementForcesAndSourcesCore *>(
      getFaceFE()->sidePtrFE);
}

MatrixDouble &FaceElementForcesAndSourcesCoreOnVolumeSideBase::
    UserDataOperator::getMasterCoordsAtGaussPts() {
  return getContactFE()->getGaussPtsMasterFromEleSide();
}

MatrixDouble &FaceElementForcesAndSourcesCoreOnVolumeSideBase::
    UserDataOperator::getSlaveCoordsAtGaussPts() {
  return getContactFE()->getGaussPtsSlaveFromEleSide();
}

int FaceElementForcesAndSourcesCoreOnVolumeSideBase::UserDataOperator::
    getFaceSense() const {
  return getFaceFE()->faceSense;
}

int FaceElementForcesAndSourcesCoreOnVolumeSideBase::UserDataOperator::
    getFaceSideNumber() const {
  return getFaceFE()->faceSideNumber;
}

template <int SWITCH>
MoFEMErrorCode FaceElementForcesAndSourcesCoreOnVolumeSideSwitch<SWITCH>::
operator()() {
  return OpSwitch<SWITCH>();
}

} // namespace MoFEM

#endif //__FACEELEMENTFORCESANDSOURCESCORE_ONVOLUMESIDE_HPP__
