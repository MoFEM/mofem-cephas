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
struct VolumeElementForcesAndSourcesCoreBase : public ForcesAndSourcesCore {

  std::string meshPositionsFieldName;

  /** \brief default operator for TET element
   * \ingroup mofem_forces_and_sources_volume_element
   */
  struct UserDataOperator : public ForcesAndSourcesCore::UserDataOperator {

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

    /**
     * \brief get measure of element
     * @return volume
     */
    inline double getMeasure() const;

    /**
     * \brief get measure of element
     * @return volume
     */
    inline double &getMeasure();

    /** \brief nodal coordinates
     */
    inline VectorDouble &getCoords();

    /** \brief coordinate at Gauss points (if hierarchical approximation of
     * element geometry)
     */
    inline MatrixDouble &getHoCoordsAtGaussPts();

    inline MatrixDouble &getHoGaussPtsJac();

    inline MatrixDouble &getHoGaussPtsInvJac();

    inline VectorDouble &getHoGaussPtsDetJac();

    inline auto getFTenosr0HoMeasure();

    /**
     * \brief Get coordinates at integration points assuming linear geometry
     *
     * \code
     * auto t_coords = getFTensor1CoordsAtGaussPts();
     * for(int gg = 0;gg!=nb_int_ptrs;gg++) {
     *   // do something
     *   ++t_coords;
     * }
     * \endcode
     *
     */
    inline auto getFTensor1CoordsAtGaussPts();

    /**
     * \brief Get coordinates at integration points taking geometry from field
     *
     * This is HO geometry given by arbitrary order polynomial
     * \code
     * auto t_coords = getFTensor1HoCoordsAtGaussPts();
     * for(int gg = 0;gg!=nb_int_ptrs;gg++) {
     *   // do something
     *   ++t_coords;
     * }
     * \endcode
     *
     */
    inline auto getFTensor1HoCoordsAtGaussPts();

    inline auto getFTensor2HoGaussPtsJac();

    inline auto getFTensor2HoGaussPtsInvJac();

    /** \brief return pointer to Generic Volume Finite Element object
     */
    inline VolumeElementForcesAndSourcesCoreBase *getVolumeFE() const;

  protected:
    MoFEMErrorCode setPtrFE(ForcesAndSourcesCore *ptr);
  };

  enum Switches {
    NO_HO_GEOMETRY = 1 << 0 | 1 << 2,
    NO_TRANSFORM = 1 << 1 | 1 << 2,
    NO_HO_TRANSFORM = 1 << 2
  };

  template <int SWITCH> MoFEMErrorCode OpSwitch();

protected:
  VolumeElementForcesAndSourcesCoreBase(Interface &m_field,
                                        const EntityType type = MBTET);

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

  /** \brief Calculate Jacobian for HO geometry
   *
   * MoFEM use hierarchical approximate base to describe geometry of the body.
   * This function transform derivatives of base functions when HO geometry is
   * set and calculate Jacobian, inverse of Jacobian and determinant of
   * transformation.
   *
   */
  virtual MoFEMErrorCode calculateHoJacobian();

  /**
   * \brief Transform base functions based on ho-geometry element Jacobian.
   *
   * This function apply transformation to base functions and its derivatives.
   * For example when base functions for H-div are present the
   * Piola-Transformarion is applied to base functions and their derivatives.
   *
   * @return Error code
   */
  virtual MoFEMErrorCode transformHoBaseFunctions();

  VectorDouble coords;
  MatrixDouble3by3 jAc;
  MatrixDouble3by3 invJac;

  OpSetInvJacH1 opSetInvJacH1;
  OpSetContravariantPiolaTransform opContravariantPiolaTransform;
  OpSetCovariantPiolaTransform opCovariantPiolaTransform;
  OpSetInvJacHdivAndHcurl opSetInvJacHdivAndHcurl;

  MatrixDouble hoCoordsAtGaussPts;
  MatrixDouble hoGaussPtsJac;
  MatrixDouble hoGaussPtsInvJac;
  VectorDouble hoGaussPtsDetJac;

  OpGetDataAndGradient<3, 3>
      opHOatGaussPoints; ///< higher order geometry data at Gauss pts
  OpSetHoInvJacH1 opSetHoInvJacH1;
  OpSetHoContravariantPiolaTransform opHoContravariantTransform;
  OpSetHoCovariantPiolaTransform opHoCovariantTransform;
  OpSetHoInvJacHdivAndHcurl opSetHoInvJacHdivAndHcurl;

  double vOlume;

  int num_nodes;
  const EntityHandle *conn;
  FTensor::Tensor2<double *, 3, 3> tJac;
  FTensor::Tensor2<double *, 3, 3> tInvJac;

  friend class UserDataOperator;
};

/**
 * @brief Volume finite element with switches
 *
 * Using SWITCH to off functions
 *
 * @tparam SWITCH
 */
template <int SWITCH>
struct VolumeElementForcesAndSourcesCoreSwitch
    : public VolumeElementForcesAndSourcesCoreBase {

  VolumeElementForcesAndSourcesCoreSwitch(Interface &m_field,
                                          const EntityType type = MBTET)
      : VolumeElementForcesAndSourcesCoreBase(m_field, MBTET) {}
  using UserDataOperator =
      VolumeElementForcesAndSourcesCoreBase::UserDataOperator;

  MoFEMErrorCode operator()();
};

/** \brief Volume finite element default
 \ingroup mofem_forces_and_sources_volume_element

 */
using VolumeElementForcesAndSourcesCore =
    VolumeElementForcesAndSourcesCoreSwitch<0>;

template <int SWITCH>
MoFEMErrorCode VolumeElementForcesAndSourcesCoreBase::OpSwitch() {
  MoFEMFunctionBegin;

  if (numeredEntFiniteElementPtr->getEntType() != MBTET)
    MoFEMFunctionReturnHot(0);
  CHKERR createDataOnElement();

  CHKERR calculateVolumeAndJacobian();
  CHKERR getSpaceBaseAndOrderOnElement();
  CHKERR setIntegrationPts();
  if (gaussPts.size2() == 0)
    MoFEMFunctionReturnHot(0);
  CHKERR calculateCoordinatesAtGaussPts();
  CHKERR calHierarchicalBaseFunctionsOnElement();
  CHKERR calBernsteinBezierBaseFunctionsOnElement();

  if (!(NO_TRANSFORM & SWITCH))
    CHKERR transformBaseFunctions();

  if (!(NO_HO_GEOMETRY & SWITCH))
    CHKERR calculateHoJacobian();

  if (!(NO_HO_TRANSFORM & SWITCH))
    CHKERR transformHoBaseFunctions();

  // Iterate over operators
  CHKERR loopOverOperators();

  MoFEMFunctionReturn(0);
}

template <int SWITCH>
MoFEMErrorCode VolumeElementForcesAndSourcesCoreSwitch<SWITCH>::operator()() {
  return OpSwitch<SWITCH>();
}

int VolumeElementForcesAndSourcesCoreBase::UserDataOperator::getNumNodes() {
  return static_cast<VolumeElementForcesAndSourcesCoreBase *>(ptrFE)->num_nodes;
}

const EntityHandle *
VolumeElementForcesAndSourcesCoreBase::UserDataOperator::getConn() {
  return static_cast<VolumeElementForcesAndSourcesCoreBase *>(ptrFE)->conn;
}

double
VolumeElementForcesAndSourcesCoreBase::UserDataOperator::getVolume() const {
  return static_cast<VolumeElementForcesAndSourcesCoreBase *>(ptrFE)->vOlume;
}

double &VolumeElementForcesAndSourcesCoreBase::UserDataOperator::getVolume() {
  return static_cast<VolumeElementForcesAndSourcesCoreBase *>(ptrFE)->vOlume;
}

FTensor::Tensor2<double *, 3, 3> &
VolumeElementForcesAndSourcesCoreBase::UserDataOperator::getJac() {
  return static_cast<VolumeElementForcesAndSourcesCoreBase *>(ptrFE)->tJac;
}

FTensor::Tensor2<double *, 3, 3> &
VolumeElementForcesAndSourcesCoreBase::UserDataOperator::getInvJac() {
  return static_cast<VolumeElementForcesAndSourcesCoreBase *>(ptrFE)->tInvJac;
}

double
VolumeElementForcesAndSourcesCoreBase::UserDataOperator::getMeasure() const {
  return getVolume();
}

double &VolumeElementForcesAndSourcesCoreBase::UserDataOperator::getMeasure() {
  return getVolume();
}

VectorDouble &
VolumeElementForcesAndSourcesCoreBase::UserDataOperator::getCoords() {
  return static_cast<VolumeElementForcesAndSourcesCoreBase *>(ptrFE)->coords;
}

MatrixDouble &VolumeElementForcesAndSourcesCoreBase::UserDataOperator::
    getHoCoordsAtGaussPts() {
  return static_cast<VolumeElementForcesAndSourcesCoreBase *>(ptrFE)
      ->hoCoordsAtGaussPts;
}

MatrixDouble &
VolumeElementForcesAndSourcesCoreBase::UserDataOperator::getHoGaussPtsJac() {
  return static_cast<VolumeElementForcesAndSourcesCoreBase *>(ptrFE)
      ->hoGaussPtsJac;
}

MatrixDouble &
VolumeElementForcesAndSourcesCoreBase::UserDataOperator::getHoGaussPtsInvJac() {
  return static_cast<VolumeElementForcesAndSourcesCoreBase *>(ptrFE)
      ->hoGaussPtsInvJac;
}

VectorDouble &
VolumeElementForcesAndSourcesCoreBase::UserDataOperator::getHoGaussPtsDetJac() {
  return static_cast<VolumeElementForcesAndSourcesCoreBase *>(ptrFE)
      ->hoGaussPtsDetJac;
}

auto VolumeElementForcesAndSourcesCoreBase::UserDataOperator::
    getFTenosr0HoMeasure() {
  return FTensor::Tensor0<FTensor::PackPtr<double *, 1>>(
      &*getHoGaussPtsDetJac().data().begin());
}

auto VolumeElementForcesAndSourcesCoreBase::UserDataOperator::
    getFTensor1CoordsAtGaussPts() {
  return FTensor::Tensor1<FTensor::PackPtr<double *, 3>, 3>(
      &getCoordsAtGaussPts()(0, 0), &getCoordsAtGaussPts()(0, 1),
      &getCoordsAtGaussPts()(0, 2));
}

auto VolumeElementForcesAndSourcesCoreBase::UserDataOperator::
    getFTensor1HoCoordsAtGaussPts() {
  double *ptr = &*getHoCoordsAtGaussPts().data().begin();
  return FTensor::Tensor1<FTensor::PackPtr<double *, 3>, 3>(ptr, ptr + 1,
                                                            ptr + 2);
}

auto VolumeElementForcesAndSourcesCoreBase::UserDataOperator::
    getFTensor2HoGaussPtsJac() {
  double *ptr = &*getHoGaussPtsJac().data().begin();
  FTensor::Tensor2<FTensor::PackPtr<double *, 9>, 3, 3> jac(
      ptr, ptr + 1, ptr + 2, ptr + 3, ptr + 4, ptr + 5, ptr + 6, ptr + 7,
      ptr + 8);
}

auto VolumeElementForcesAndSourcesCoreBase::UserDataOperator::
    getFTensor2HoGaussPtsInvJac() {
  double *ptr = &*getHoGaussPtsInvJac().data().begin();
  FTensor::Tensor2<FTensor::PackPtr<double *, 9>, 3, 3> jac(
      ptr, ptr + 1, ptr + 2, ptr + 3, ptr + 4, ptr + 5, ptr + 6, ptr + 7,
      ptr + 8);
}

VolumeElementForcesAndSourcesCoreBase *
VolumeElementForcesAndSourcesCoreBase::UserDataOperator::getVolumeFE() const {
  return static_cast<VolumeElementForcesAndSourcesCoreBase *>(ptrFE);
}

} // namespace MoFEM

#endif //__VOLUMEELEMENTFORCESANDSOURCESCORE_HPP__

/**
 * \defgroup mofem_forces_and_sources_volume_element Volume Element
 * \brief Implementation of general volume element.
 *
 * \ingroup mofem_forces_and_sources
 **/
