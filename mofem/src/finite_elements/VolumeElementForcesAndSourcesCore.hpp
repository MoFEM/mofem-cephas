/** \file VolumeElementForcesAndSourcesCore.hpp
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

#ifndef __VOLUMEELEMENTFORCESANDSOURCESCORE_HPP__
#define __VOLUMEELEMENTFORCESANDSOURCESCORE_HPP__

using namespace boost::numeric;

namespace MoFEM {

/** \brief Volume finite element base
 \ingroup mofem_forces_and_sources_volume_element

 User is implementing own operator at Gauss point level, by class
 derived from VolumeElementForcesAndSourcesCore::UserDataOperator. Arbitrary
 number of operator can be added by pushing objects to OpPtrVector

 */
struct VolumeElementForcesAndSourcesCore : public ForcesAndSourcesCore {

  VolumeElementForcesAndSourcesCore(Interface &m_field,
                                        const EntityType type = MBTET);


  std::string meshPositionsFieldName; ///< \deprecated DO NOT USE!

  /** \brief default operator for TET element
   * \ingroup mofem_forces_and_sources_volume_element
   */
  struct UserDataOperator;

  enum Switches {
  };

  MoFEMErrorCode operator()();

protected:

  // Note that functions below could be overloaded by user to change default
  // behavior of the element.

  /**
   * \brief Set integration points
   * @return Error code
   */
  virtual MoFEMErrorCode setIntegrationPts();

  /**
   * \brief Calculate element volume and Jacobian
   *
   * Note that at that point is assumed that geometry is exclusively defined by
   * corner nodes.
   *
   * @return Error code
   */
  virtual MoFEMErrorCode calculateVolumeAndJacobian();

  /**
   * \brief Calculate coordinate at integration points
   * @return Error code
   */
  virtual MoFEMErrorCode calculateCoordinatesAtGaussPts();

  /**
   * \brief Determine approximation space and order of base functions
   * @return Error code
   */
  virtual MoFEMErrorCode getSpaceBaseAndOrderOnElement();

  /**
   * \brief Transform base functions based on geometric element Jacobian.
   *
   * This function apply transformation to base functions and its derivatives.
   * For example when base functions for H-div are present the
   * Piola-Transformarion is applied to base functions and their derivatives.
   *
   * @return Error code
   */
  virtual MoFEMErrorCode transformBaseFunctions();


  VectorDouble coords;
  MatrixDouble3by3 jAc;
  MatrixDouble3by3 invJac;

  OpSetInvJacH1 opSetInvJacH1;
  OpSetContravariantPiolaTransform opContravariantPiolaTransform;
  OpSetCovariantPiolaTransform opCovariantPiolaTransform;
  OpSetInvJacHdivAndHcurl opSetInvJacHdivAndHcurl;

  double &vOlume;

  int num_nodes;
  const EntityHandle *conn;
  FTensor::Tensor2<double *, 3, 3> tJac;
  FTensor::Tensor2<double *, 3, 3> tInvJac;

  friend class UserDataOperator;
};

struct VolumeElementForcesAndSourcesCore::UserDataOperator
    : public ForcesAndSourcesCore::UserDataOperator {

  using ForcesAndSourcesCore::UserDataOperator::UserDataOperator;

  /** \brief get element number of nodes
   */
  inline int getNumNodes();

  /** \brief get element connectivity
   */
  inline const EntityHandle *getConn();

  /** \brief element volume (linear geometry)
   */
  inline double getVolume() const;

  /** \brief element volume (linear geometry)
   */
  inline double &getVolume();

  /**
   * \brief get element Jacobian
   */
  inline FTensor::Tensor2<double *, 3, 3> &getJac();

  /**
   * \brief get element inverse Jacobian
   */
  inline FTensor::Tensor2<double *, 3, 3> &getInvJac();

  /** \brief nodal coordinates
   */
  inline VectorDouble &getCoords();

  /** \brief return pointer to Generic Volume Finite Element object
   */
  inline VolumeElementForcesAndSourcesCore *getVolumeFE() const;

protected:
  MoFEMErrorCode setPtrFE(ForcesAndSourcesCore *ptr);
};

int VolumeElementForcesAndSourcesCore::UserDataOperator::getNumNodes() {
  return static_cast<VolumeElementForcesAndSourcesCore *>(ptrFE)->num_nodes;
}

const EntityHandle *
VolumeElementForcesAndSourcesCore::UserDataOperator::getConn() {
  return static_cast<VolumeElementForcesAndSourcesCore *>(ptrFE)->conn;
}

double
VolumeElementForcesAndSourcesCore::UserDataOperator::getVolume() const {
  return static_cast<VolumeElementForcesAndSourcesCore *>(ptrFE)->vOlume;
}

double &VolumeElementForcesAndSourcesCore::UserDataOperator::getVolume() {
  return static_cast<VolumeElementForcesAndSourcesCore *>(ptrFE)->vOlume;
}

FTensor::Tensor2<double *, 3, 3> &
VolumeElementForcesAndSourcesCore::UserDataOperator::getJac() {
  return static_cast<VolumeElementForcesAndSourcesCore *>(ptrFE)->tJac;
}

FTensor::Tensor2<double *, 3, 3> &
VolumeElementForcesAndSourcesCore::UserDataOperator::getInvJac() {
  return static_cast<VolumeElementForcesAndSourcesCore *>(ptrFE)->tInvJac;
}

VectorDouble &
VolumeElementForcesAndSourcesCore::UserDataOperator::getCoords() {
  return static_cast<VolumeElementForcesAndSourcesCore *>(ptrFE)->coords;
}

VolumeElementForcesAndSourcesCore *
VolumeElementForcesAndSourcesCore::UserDataOperator::getVolumeFE() const {
  return static_cast<VolumeElementForcesAndSourcesCore *>(ptrFE);
}

/**
 * @deprecated use VolumeElementForcesAndSourcesCore
 * 
 */
DEPRECATED typedef VolumeElementForcesAndSourcesCore
    VolumeElementForcesAndSourcesCoreBase;

/**
 * @deprecated obsolete template. it will be removed in next versions
 */
template <int SWITCH>
struct VolumeElementForcesAndSourcesCoreSwitch
    : public VolumeElementForcesAndSourcesCore {

  using VolumeElementForcesAndSourcesCore::VolumeElementForcesAndSourcesCore;
  using UserDataOperator = VolumeElementForcesAndSourcesCore::UserDataOperator;
};

} // namespace MoFEM

#endif //__VOLUMEELEMENTFORCESANDSOURCESCORE_HPP__

/**
 * \defgroup mofem_forces_and_sources_volume_element Volume Element
 * \brief Implementation of general volume element.
 *
 * \ingroup mofem_forces_and_sources
 **/
