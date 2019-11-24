/** \file VolumeElementForcesAndSourcesCoreOnVolumeSide.hpp
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

#ifndef __VOLUMEELEMENTFORCESANDSOURCESCORE_ONVOLUMESIDE_HPP__
#define __VOLUMEELEMENTFORCESANDSOURCESCORE_ONVOLUMESIDE_HPP__

using namespace boost::numeric;

namespace MoFEM {

/**
 * \brief Base volume element used to integrate on contact surface (could be
 * extended to other volume elements) \ingroup
 * mofem_forces_and_sources_volume_element
 */
struct VolumeElementForcesAndSourcesCoreOnVolumeSideBase
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
   * \todo That this is not general, e.g., for prism or hex, opposite node is
   * not unique.
   *
   * @return int
   */
  inline int getOppositeNode() const;

  /** \brief default operator for TET element
   * \ingroup mofem_forces_and_sources_volume_element
   */
  struct UserDataOperator
      : public VolumeElementForcesAndSourcesCore::UserDataOperator {

    using VolumeElementForcesAndSourcesCore::UserDataOperator::UserDataOperator;

    inline VolumeElementForcesAndSourcesCoreOnVolumeSideBase *getVolumeFE() const;

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
struct VolumeElementForcesAndSourcesCoreOnVolumeSideSwitch
    : public VolumeElementForcesAndSourcesCoreOnVolumeSideBase {

  using VolumeElementForcesAndSourcesCoreOnVolumeSideBase::
      VolumeElementForcesAndSourcesCoreOnVolumeSideBase;

  using UserDataOperator =
      VolumeElementForcesAndSourcesCoreOnVolumeSideBase::UserDataOperator;

  MoFEMErrorCode operator()();
};

/** \brief Volume element used to integrate on contact element (could be extended for other volume elements)
 \ingroup mofem_forces_and_sources_volume_element

 */
using VolumeElementForcesAndSourcesCoreOnVolumeSide =
    VolumeElementForcesAndSourcesCoreOnVolumeSideSwitch<0>;

const std::array<int, 3> &
VolumeElementForcesAndSourcesCoreOnVolumeSideBase::getFaceConnMap() const {
  return faceConnMap;
}

const std::array<int, 4> &
VolumeElementForcesAndSourcesCoreOnVolumeSideBase::getTetConnMap() const {
  return tetConnMap;
}

int VolumeElementForcesAndSourcesCoreOnVolumeSideBase::getOppositeNode() const {
  return oppositeNode;
}

VolumeElementForcesAndSourcesCoreOnVolumeSideBase *
VolumeElementForcesAndSourcesCoreOnVolumeSideBase::UserDataOperator::getVolumeFE()
    const {
  return static_cast<VolumeElementForcesAndSourcesCoreOnVolumeSideBase *>(ptrFE);
}

template <int SWITCH>
MoFEMErrorCode VolumeElementForcesAndSourcesCoreOnVolumeSideSwitch<SWITCH>::
operator()() {
  return OpSwitch<SWITCH>();
}

} // namespace MoFEM

#endif //__VOLUMEELEMENTFORCESANDSOURCESCORE_ONVOLUMESIDE_HPP__
