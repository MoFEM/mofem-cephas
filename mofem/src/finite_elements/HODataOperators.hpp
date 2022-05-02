/** \file HODataOperators.hpp
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
 * @brief Calculate jacobian on Hex or other volume which is not simplex.
 *
 * Using not a field data, but geometrical coordinates of element.
 *
 * \note Use OpCalculateVectorFieldGradient to calculate Jacobian from field
 * data.
 *
 */
struct OpCalculateHOJacVolume
    : public VolumeElementForcesAndSourcesCoreBase::UserDataOperator {

  OpCalculateHOJacVolume(boost::shared_ptr<MatrixDouble> jac_ptr);

  MoFEMErrorCode doWork(int side, EntityType type,
                        EntitiesFieldData::EntData &data);

private:
  boost::shared_ptr<MatrixDouble> jacPtr;
};

/**
 * @brief Calculate HO coordinates at gauss points
 *
 */
struct OpCalculateHOCoords : public ForcesAndSourcesCore::UserDataOperator {

  OpCalculateHOCoords(const std::string field_name)
      : ForcesAndSourcesCore::UserDataOperator(field_name, OPROW) {}

  MoFEMErrorCode doWork(int side, EntityType type,
                        EntitiesFieldData::EntData &data);
};

/**
 * \brief Set inverse jacobian to base functions
 *
 */
struct OpSetHOInvJacToScalarBases
    : public ForcesAndSourcesCore::UserDataOperator {

  OpSetHOInvJacToScalarBases(const FieldSpace space,
                             boost::shared_ptr<MatrixDouble> inv_jac_ptr)
      : ForcesAndSourcesCore::UserDataOperator(space), invJacPtr(inv_jac_ptr) {

    if (!inv_jac_ptr)
      CHK_THROW_MESSAGE(MOFEM_DATA_INCONSISTENCY, "invJacPtr not allocated");

    if (space == L2) {
      doVertices = false;
    }
  }

  MoFEMErrorCode doWork(int side, EntityType type,
                        EntitiesFieldData::EntData &data);

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

    if (!invJacPtr)
      CHK_THROW_MESSAGE(MOFEM_DATA_INCONSISTENCY,
                        "Pointer for invJacPtr not allocated");

    doVertices = false;
    if (space == HDIV)
      doEdges = false;
  }

  MoFEMErrorCode doWork(int side, EntityType type,
                        EntitiesFieldData::EntData &data);

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
struct OpSetHOWeightsOnFace
    : public FaceElementForcesAndSourcesCoreBase::UserDataOperator {
  OpSetHOWeightsOnFace()
      : FaceElementForcesAndSourcesCoreBase::UserDataOperator(NOSPACE) {}
  MoFEMErrorCode doWork(int side, EntityType type,
                        EntitiesFieldData::EntData &data);
};

/**
 * \brief Set inverse jacobian to base functions
 *
 */
struct OpSetHOWeights : public ForcesAndSourcesCore::UserDataOperator {

  OpSetHOWeights(boost::shared_ptr<VectorDouble> det_ptr)
      : ForcesAndSourcesCore::UserDataOperator(NOSPACE, OPSPACE),
        detPtr(det_ptr) {
    if (!detPtr)
      CHK_THROW_MESSAGE(MOFEM_DATA_INCONSISTENCY,
                        "Pointer for detPtr not allocated");
  }

  MoFEMErrorCode doWork(int side, EntityType type,
                        EntitiesFieldData::EntData &data);

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
      : ForcesAndSourcesCore::UserDataOperator(space, OPSPACE), detPtr(det_ptr),
        jacPtr(jac_ptr) {
    doVertices = false;
    if (space == HDIV)
      doEdges = false;
  }

  MoFEMErrorCode doWork(int side, EntityType type,
                        EntitiesFieldData::EntData &data);

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
      : ForcesAndSourcesCore::UserDataOperator(space, OPSPACE),
        jacInvPtr(jac_inv_ptr) {
    doVertices = false;
    if (space == HDIV)
      doEdges = false;
    if (!jacInvPtr)
      CHK_THROW_MESSAGE(MOFEM_DATA_INCONSISTENCY,
                        "Pointer for jacPtr not allocated");
  }

  MoFEMErrorCode doWork(int side, EntityType type,
                        EntitiesFieldData::EntData &data);

private:
  boost::shared_ptr<MatrixDouble> jacInvPtr;

  MatrixDouble piolaN;
  MatrixDouble piolaDiffN;
};

/** \brief Calculate jacobian for face element

  It is assumed that face element is XY plane. Applied
  only for 2d problems and 2d problems embedded in 3d space

  \note If you push operators for HO normal befor this operator, HO geometry is
  taken into account when you calculate jacobian.

*/
template <int DIM> struct OpCalculateHOJacForFaceImpl;

template <>
struct OpCalculateHOJacForFaceImpl<2>
    : public FaceElementForcesAndSourcesCoreBase::UserDataOperator {

  OpCalculateHOJacForFaceImpl(boost::shared_ptr<MatrixDouble> jac_ptr);

  MoFEMErrorCode doWork(int side, EntityType type,
                        EntitiesFieldData::EntData &data);

protected:
  boost::shared_ptr<MatrixDouble> jacPtr;
};

template <>
struct OpCalculateHOJacForFaceImpl<3> : public OpCalculateHOJacForFaceImpl<2> {

  using OpCalculateHOJacForFaceImpl<2>::OpCalculateHOJacForFaceImpl;

  MoFEMErrorCode doWork(int side, EntityType type,
                        EntitiesFieldData::EntData &data);
};

using OpCalculateHOJacForFace = OpCalculateHOJacForFaceImpl<2>;
using OpCalculateHOJacForFaceEmbeddedIn3DSpace = OpCalculateHOJacForFaceImpl<3>;

/** \brief Calculate normals at Gauss points of triangle element
 * \ingroup mofem_forces_and_source
 */
struct OpGetHONormalsOnFace
    : public FaceElementForcesAndSourcesCoreBase::UserDataOperator {

  OpGetHONormalsOnFace(std::string field_name);

  MoFEMErrorCode doWork(int side, EntityType type,
                        EntitiesFieldData::EntData &data);
};

/** \brief transform Hdiv base fluxes from reference element to physical
 * triangle \ingroup mofem_forces_and_sources
 *
 * \note Matrix which keeps normal is assumed to have three columns, and number
 * of rows should be equal to number of integration points.
 *
 */
struct OpHOSetContravariantPiolaTransformOnFace3D
    : public FaceElementForcesAndSourcesCoreBase::UserDataOperator {

  OpHOSetContravariantPiolaTransformOnFace3D(
      const FieldSpace space,
      boost::shared_ptr<MatrixDouble> normals_at_gauss_pts = nullptr)
      : FaceElementForcesAndSourcesCoreBase::UserDataOperator(space, OPSPACE),
        normalsAtGaussPts(normals_at_gauss_pts) {}

  MoFEMErrorCode doWork(int side, EntityType type,
                        EntitiesFieldData::EntData &data);

private:
  boost::shared_ptr<MatrixDouble> normalsAtGaussPts;
};

/** \brief transform Hcurl base fluxes from reference element to physical edge
 * \ingroup mofem_forces_and_sources
 */
struct OpHOSetContravariantPiolaTransformOnEdge3D
    : public EdgeElementForcesAndSourcesCoreBase::UserDataOperator {

  OpHOSetContravariantPiolaTransformOnEdge3D(
      const FieldSpace space = HCURL,
      boost::shared_ptr<MatrixDouble> tangent_at_pts = nullptr)
      : EdgeElementForcesAndSourcesCoreBase::UserDataOperator(space),
        tangentAtGaussPts(tangent_at_pts) {}

  MoFEMErrorCode doWork(int side, EntityType type,
                        EntitiesFieldData::EntData &data);

private:
  boost::shared_ptr<MatrixDouble> tangentAtGaussPts;
};

/** \brief transform Hcurl base fluxes from reference element to physical
 * triangle \ingroup mofem_forces_and_sources
 */
struct OpHOSetCovariantPiolaTransformOnFace3D
    : public FaceElementForcesAndSourcesCoreBase::UserDataOperator {

  OpHOSetCovariantPiolaTransformOnFace3D(
      const FieldSpace space,
      boost::shared_ptr<MatrixDouble> normals_at_pts = nullptr,
      boost::shared_ptr<MatrixDouble> tangent1_at_pts = nullptr,
      boost::shared_ptr<MatrixDouble> tangent2_at_pts = nullptr)
      : FaceElementForcesAndSourcesCoreBase::UserDataOperator(space, OPSPACE),
        normalsAtPts(normals_at_pts), tangent1AtPts(tangent1_at_pts),
        tangent2AtPts(tangent2_at_pts) {
    if (normals_at_pts || tangent1_at_pts || tangent2_at_pts)
      if (normals_at_pts && tangent1_at_pts && tangent2_at_pts)
        CHK_THROW_MESSAGE(MOFEM_DATA_INCONSISTENCY,
                          "All elements in constructor have to set pointer");
  }

  MoFEMErrorCode doWork(int side, EntityType type,
                        EntitiesFieldData::EntData &data);

private:
  boost::shared_ptr<MatrixDouble> normalsAtPts;
  boost::shared_ptr<MatrixDouble> tangent1AtPts;
  boost::shared_ptr<MatrixDouble> tangent2AtPts;

  MatrixDouble piolaN;
  MatrixDouble diffPiolaN;
};

/** \brief Calculate tangent vector on edge form HO geometry approximation
 * \ingroup mofem_forces_and_sources
 */
struct OpGetHOTangentsOnEdge
    : public EdgeElementForcesAndSourcesCoreBase::UserDataOperator {

  OpGetHOTangentsOnEdge(
      std::string field_name,
      boost::shared_ptr<MatrixDouble> tangents_at_pts = nullptr)
      : EdgeElementForcesAndSourcesCoreBase::UserDataOperator(field_name,
                                                              OPROW),
        tangentsAtPts(tangents_at_pts) {}

  MoFEMErrorCode doWork(int side, EntityType type,
                        EntitiesFieldData::EntData &data);

private:
  boost::shared_ptr<MatrixDouble> tangentsAtPts;
};

/**
 * @brief Scale base functions by inverses of measure of element
 *
 * @tparam OP
 */
template <typename OP> struct OpScaleBaseBySpaceInverseOfMeasure : public OP {

  OpScaleBaseBySpaceInverseOfMeasure(
      boost::shared_ptr<VectorDouble> det_jac_ptr, const FieldSpace space = L2)
      : OP(space), fieldSpace(space), detJacPtr(det_jac_ptr) {}

  MoFEMErrorCode doWork(int side, EntityType type,
                        EntitiesFieldData::EntData &data) {
    MoFEMFunctionBegin;

    if (!detJacPtr) {
      SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY, "detJacPtr not set");
    }

    auto scale = [&]() {
      for (int b = AINSWORTH_LEGENDRE_BASE; b != USER_BASE; b++) {
        FieldApproximationBase base = static_cast<FieldApproximationBase>(b);
        auto &base_fun = data.getN(base);
        auto &diff_base_fun = data.getDiffN(base);
        if (detJacPtr) {

          auto &det_vec = *detJacPtr;
          const auto nb_base_fun = base_fun.size2();
          const auto nb_int_pts = base_fun.size1();

          if (nb_int_pts != det_vec.size())
            SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
                    "Number of integration pts in detJacPtr does not mush "
                    "number of integration pts in base function");

          auto get_row = [](auto &m, auto gg) {
            return ublas::matrix_row<decltype(m)>(m, gg);
          };

          for (auto gg = 0; gg != nb_int_pts; ++gg)
            get_row(base_fun) /= det_vec[gg];

          if (diff_base_fun.size1() == nb_int_pts) {
            for (auto gg = 0; gg != nb_int_pts; ++gg)
              get_row(diff_base_fun) /= det_vec[gg];
          }
        }
      }
    };

    if (this->getFEDim() == 3) {
      switch (fieldSpace) {
      case H1:
        scale();
        break;
      case HCURL:
        if (type >= MBEDGE)
          scale();
        break;
      case HDIV:
        if (type >= MBTRI)
          scale();
        break;
      case L2:
        if (type >= MBTET)
          scale();
        break;
      default:
        SETERRQ(PETSC_COMM_SELF, MOFEM_IMPOSIBLE_CASE, "impossible case");
      }
    } else if (this->getFEDim() == 2) {
      switch (fieldSpace) {
      case H1:
        scale();
        break;
      case HCURL:
        if (type >= MBEDGE)
          scale();
        break;
      case HDIV:
      case L2:
        if (type >= MBTRI)
          scale();
        break;
      default:
        SETERRQ(PETSC_COMM_SELF, MOFEM_IMPOSIBLE_CASE, "impossible case");
      }
    } else if (this->getFEDim() == 1) {
      switch (fieldSpace) {
      case H1:
        scale();
        break;
      case HCURL:
      case L2:
        if (type >= MBEDGE)
          scale();
        break;
      default:
        SETERRQ(PETSC_COMM_SELF, MOFEM_IMPOSIBLE_CASE, "impossible case");
      }
    } else {
      SETERRQ(PETSC_COMM_SELF, MOFEM_IMPOSIBLE_CASE, "impossible case");
    }

    MoFEMFunctionReturn(0);
  }

private:
  FieldSpace fieldSpace;
  boost::shared_ptr<VectorDouble> detJacPtr;
};

template <typename E>
MoFEMErrorCode addHOOpsVol(const std::string field, E &e, bool h1, bool hcurl,
                           bool hdiv, bool l2) {
  MoFEMFunctionBegin;
  auto material_grad_mat = boost::make_shared<MatrixDouble>();
  auto material_det_vec = boost::make_shared<VectorDouble>();
  auto material_inv_grad_mat = boost::make_shared<MatrixDouble>();
  e.getOpPtrVector().push_back(
      new OpCalculateVectorFieldGradient<3, 3>(field, material_grad_mat));
  e.getOpPtrVector().push_back(new OpInvertMatrix<3>(
      material_grad_mat, material_det_vec, material_inv_grad_mat));
  e.getOpPtrVector().push_back(new OpSetHOWeights(material_det_vec));
  if (h1)
    e.getOpPtrVector().push_back(
        new OpSetHOInvJacToScalarBases(H1, material_inv_grad_mat));
  if (l2)
    e.getOpPtrVector().push_back(
        new OpSetHOInvJacToScalarBases(L2, material_inv_grad_mat));
  if (hdiv) {
    e.getOpPtrVector().push_back(new OpSetHOContravariantPiolaTransform(
        HDIV, material_det_vec, material_grad_mat));
    e.getOpPtrVector().push_back(
        new OpSetHOInvJacVectorBase(HDIV, material_inv_grad_mat));
  }
  if (hcurl) {
    e.getOpPtrVector().push_back(
        new OpSetHOCovariantPiolaTransform(HCURL, material_inv_grad_mat));
    e.getOpPtrVector().push_back(
        new OpSetHOInvJacVectorBase(HCURL, material_inv_grad_mat));
  }
  MoFEMFunctionReturn(0);
}

template <typename E>
MoFEMErrorCode addHOOpsFace3D(const std::string field, E &e, bool hcurl,
                              bool hdiv) {
  MoFEMFunctionBegin;
  e.meshPositionsFieldName = "none";
  e.getOpPtrVector().push_back(new OpGetHONormalsOnFace(field));
  e.getOpPtrVector().push_back(new OpCalculateHOCoords(field));
  if (hcurl) {
    e.getOpPtrVector().push_back(
        new OpHOSetContravariantPiolaTransformOnFace3D(HDIV));
  }
  if (hdiv) {
    e.getOpPtrVector().push_back(
        new OpHOSetCovariantPiolaTransformOnFace3D(HDIV));
  }
  MoFEMFunctionReturn(0);
}

}; // namespace MoFEM

#endif