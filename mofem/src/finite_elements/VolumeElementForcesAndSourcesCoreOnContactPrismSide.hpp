/** \file VolumeElementForcesAndSourcesCoreOnContactPrismSide.hpp
  \brief Volume element.

  Those element are inherited by user to implement specific implementation of
  particular problem.

*/



#ifndef __VOLUMEELEMENTFORCESANDSOURCESCORE_ONVOLUMESIDE_HPP__
#define __VOLUMEELEMENTFORCESANDSOURCESCORE_ONVOLUMESIDE_HPP__

using namespace boost::numeric;

namespace MoFEM {

/**
 * \brief Base volume element used to integrate on contact surface (could be
 * extended to other volume elements) \ingroup
 * mofem_forces_and_sources_volume_element
 */
struct VolumeElementForcesAndSourcesCoreOnContactPrismSide
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
      : public VolumeElementForcesAndSourcesCore::UserDataOperator {

    using VolumeElementForcesAndSourcesCore::UserDataOperator::UserDataOperator;

    inline VolumeElementForcesAndSourcesCoreOnContactPrismSide *
    getVolumeFE() const;

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
  std::array<int, 4> faceConnMap;
  std::array<int, 8> tetConnMap;
  int oppositeNode;
};

const std::array<int, 4> &
VolumeElementForcesAndSourcesCoreOnContactPrismSide::getFaceConnMap() const {
  return faceConnMap;
}

const std::array<int, 8> &
VolumeElementForcesAndSourcesCoreOnContactPrismSide::getTetConnMap() const {
  return tetConnMap;
}

int VolumeElementForcesAndSourcesCoreOnContactPrismSide::getOppositeNode() const {
  return oppositeNode;
}

int VolumeElementForcesAndSourcesCoreOnContactPrismSide::getFaceSense() const {
  return faceSense;
}

int VolumeElementForcesAndSourcesCoreOnContactPrismSide::getFaceSideNumber()
    const {
  return faceSideNumber;
}

VolumeElementForcesAndSourcesCoreOnContactPrismSide *
VolumeElementForcesAndSourcesCoreOnContactPrismSide::UserDataOperator::
    getVolumeFE() const {
  return static_cast<VolumeElementForcesAndSourcesCoreOnContactPrismSide *>(
      ptrFE);
}

ContactPrismElementForcesAndSourcesCore *
VolumeElementForcesAndSourcesCoreOnContactPrismSide::UserDataOperator::
    getContactFE() const {
  return static_cast<ContactPrismElementForcesAndSourcesCore *>(
      getVolumeFE()->sidePtrFE);
}

MatrixDouble &VolumeElementForcesAndSourcesCoreOnContactPrismSide::
    UserDataOperator::getMasterCoordsAtGaussPts() {
  return getContactFE()->getGaussPtsMasterFromEleSide();
}

MatrixDouble &VolumeElementForcesAndSourcesCoreOnContactPrismSide::
    UserDataOperator::getSlaveCoordsAtGaussPts() {
  return getContactFE()->getGaussPtsSlaveFromEleSide();
}

int VolumeElementForcesAndSourcesCoreOnContactPrismSide::UserDataOperator::
    getFaceSense() const {
  return getVolumeFE()->faceSense;
}

int VolumeElementForcesAndSourcesCoreOnContactPrismSide::UserDataOperator::
    getFaceSideNumber() const {
  return getVolumeFE()->faceSideNumber;
}

/**
 * @deprecated use VolumeElementForcesAndSourcesCore
 * 
 */
DEPRECATED typedef VolumeElementForcesAndSourcesCoreOnContactPrismSide
    VolumeElementForcesAndSourcesCoreOnContactPrismSideBase;

/**
 * @deprecated do not use this template, it is obsolete
 */
template <int SWITCH>
struct VolumeElementForcesAndSourcesCoreOnContactPrismSideSwitch
    : public VolumeElementForcesAndSourcesCoreOnContactPrismSide {
  using VolumeElementForcesAndSourcesCoreOnContactPrismSide::
      VolumeElementForcesAndSourcesCoreOnContactPrismSide;
  using UserDataOperator =
      VolumeElementForcesAndSourcesCoreOnSide::UserDataOperator;
};

} // namespace MoFEM

#endif //__VOLUMEELEMENTFORCESANDSOURCESCORE_ONVOLUMESIDE_HPP__
