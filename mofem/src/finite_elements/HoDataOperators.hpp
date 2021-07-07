/** \file HoDataOperators.hpp
  * \brief Operators managing HO geometry

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

#ifndef __HO_DATA_OPERATORS_HPP__
#define __HO_DATA_OPERATORS_HPP__

namespace MoFEM {

/**
 * @brief Calculate HO coordinates at gauss points
 *
 */
struct OpCalculateHoCoords : public ForcesAndSourcesCore::UserDataOperator {

  OpCalculateHoCoords(const std::string field_name)
      : ForcesAndSourcesCore::UserDataOperator(field_name, OPROW) {}

  MoFEMErrorCode doWork(int side, EntityType type,
                        DataForcesAndSourcesCore::EntData &data);
};

/**
 * \brief Set inverse jacobian to base functions
 *
 */
struct OpSetHOInvJacToScalarBases
    : public ForcesAndSourcesCore::UserDataOperator {

  OpSetHOInvJacToScalarBases(const FieldSpace space,
                             boost::shared_ptr<MatrixDouble> inv_jac_ptr)
      : ForcesAndSourcesCore::UserDataOperator(space), invJacPtr(inv_jac_ptr) {}

  MoFEMErrorCode doWork(int side, EntityType type,
                        DataForcesAndSourcesCore::EntData &data);

private:
  boost::shared_ptr<MatrixDouble> invJacPtr;
  MatrixDouble diffNinvJac;
};

/**
 * \brief transform local reference derivatives of shape function to global
 derivatives if higher order geometry is given
 *

 * \ingroup mofem_forces_and_sources
*/
struct OpSetHOInvJacVectorBase : public ForcesAndSourcesCore::UserDataOperator {

  OpSetHOInvJacVectorBase(const FieldSpace space,
                          boost::shared_ptr<MatrixDouble> inv_jac_ptr)
      : ForcesAndSourcesCore::UserDataOperator(space), invJacPtr(inv_jac_ptr) {
    doVertices = false;
    if (space == HDIV)
      doEdges = false;
  }

  MoFEMErrorCode doWork(int side, EntityType type,
                        DataForcesAndSourcesCore::EntData &data);

private:
  boost::shared_ptr<MatrixDouble> invJacPtr;
  MatrixDouble diffHdivInvJac;
  MatrixDouble diffNinvJac;
};

/**
 * @brief Modify integration weights on face to take in account higher-order
 * geometry
 * @ingroup mofem_forces_and_sources_tri_element
 *
 */
struct OpMakeHighOrderGeometryWeightsOnFace
    : public FaceElementForcesAndSourcesCoreBase::UserDataOperator {
  OpMakeHighOrderGeometryWeightsOnFace()
      : FaceElementForcesAndSourcesCoreBase::UserDataOperator(NOSPACE) {}
  MoFEMErrorCode doWork(int side, EntityType type,
                        DataForcesAndSourcesCore::EntData &data);
};

/**
 * \brief Set inverse jacobian to base functions
 *
 */
struct OpSetHOWeights : public ForcesAndSourcesCore::UserDataOperator {

  OpSetHOWeights(boost::shared_ptr<VectorDouble> det_ptr)
      : ForcesAndSourcesCore::UserDataOperator(NOSPACE, OPLAST),
        detPtr(det_ptr) {}

  MoFEMErrorCode doWork(int side, EntityType type,
                        DataForcesAndSourcesCore::EntData &data);

private:
  boost::shared_ptr<VectorDouble> detPtr;
};

/** \brief Apply contravariant (Piola) transfer to Hdiv space for HO geometr

* \ingroup mofem_forces_and_sources
*/
struct OpSetHOContravariantPiolaTransform
    : public ForcesAndSourcesCore::UserDataOperator {

  OpSetHOContravariantPiolaTransform(const FieldSpace space,
                                     boost::shared_ptr<VectorDouble> det_ptr,
                                     boost::shared_ptr<MatrixDouble> jac_ptr)
      : ForcesAndSourcesCore::UserDataOperator(space, OPLAST), detPtr(det_ptr),
        jacPtr(jac_ptr) {
    doVertices = false;
    if (space == HDIV)
      doEdges = false;
  }

  MoFEMErrorCode doWork(int side, EntityType type,
                        DataForcesAndSourcesCore::EntData &data);

private:
  boost::shared_ptr<VectorDouble> detPtr;
  boost::shared_ptr<MatrixDouble> jacPtr;

  MatrixDouble piolaN;
  MatrixDouble piolaDiffN;
};

/** \brief Apply covariant (Piola) transfer to Hcurl space for HO geometry
 */
struct OpSetHOCovariantPiolaTransform
    : public ForcesAndSourcesCore::UserDataOperator {

  OpSetHOCovariantPiolaTransform(const FieldSpace space,
                                 boost::shared_ptr<MatrixDouble> jac_inv_ptr)
      : ForcesAndSourcesCore::UserDataOperator(space, OPLAST),
        jacInvPtr(jac_inv_ptr) {
    doVertices = false;
    if (space == HDIV)
      doEdges = false;
  }

  MoFEMErrorCode doWork(int side, EntityType type,
                        DataForcesAndSourcesCore::EntData &data);

private:
  boost::shared_ptr<MatrixDouble> jacInvPtr;

  MatrixDouble piolaN;
  MatrixDouble piolaDiffN;
};

}; // namespace MoFEM

#endif