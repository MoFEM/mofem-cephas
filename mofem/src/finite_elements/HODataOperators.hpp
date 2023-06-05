/** \file HODataOperators.hpp
  * \brief Operators managing HO geometry

*/

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
struct OpCalculateHOJacForVolume
    : public VolumeElementForcesAndSourcesCore::UserDataOperator {

  OpCalculateHOJacForVolume(boost::shared_ptr<MatrixDouble> jac_ptr);

  MoFEMErrorCode doWork(int side, EntityType type,
                        EntitiesFieldData::EntData &data);

private:
  boost::shared_ptr<MatrixDouble> jacPtr;
};

/** \deprecated use OpCalculateHOJacForVolume
 */
using OpCalculateHOJacVolume = OpCalculateHOJacForVolume;

/**
 * @brief Calculate HO coordinates at gauss points
 *
 */
template <int FIELD_DIM = 3>
struct OpCalculateHOCoords : public ForcesAndSourcesCore::UserDataOperator {

  OpCalculateHOCoords(const std::string field_name)
      : ForcesAndSourcesCore::UserDataOperator(field_name, OPROW) {}

  MoFEMErrorCode doWork(int side, EntityType type,
                        EntitiesFieldData::EntData &data);
};

/**
 * @deprecated This class should be DIM = 3 specialization when default
 * parameter is removed
 *
 */
struct OpSetHOInvJacToScalarBasesImpl : public OpSetInvJacToScalarBasesBasic {

  using OpSetInvJacToScalarBasesBasic::OpSetInvJacToScalarBasesBasic;

  MoFEMErrorCode doWork(int side, EntityType type,
                        EntitiesFieldData::EntData &data);
};

/**
 * \brief Set inverse jacobian to base functions
 *
 * \deprecated  Version with default variant DIM = 3 is deprecated, keep for
 * back compatibility. Should be removed in future. Use of default variant make
 * code implicit, what can be source of some hidden error. Explict interface is
 * better, when user see and control outcome, and is aware of existing variants.
 *
 */
template <int DIM = 3>
struct OpSetHOInvJacToScalarBases : public OpSetHOInvJacToScalarBasesImpl {
  using OpSetHOInvJacToScalarBasesImpl::OpSetHOInvJacToScalarBasesImpl;
};

template <>
struct OpSetHOInvJacToScalarBases<2>
    : public OpSetInvJacSpaceForFaceImpl<2, 1> {

  using OpSetInvJacSpaceForFaceImpl<2, 1>::OpSetInvJacSpaceForFaceImpl;
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
    : public FaceElementForcesAndSourcesCore::UserDataOperator {
  OpSetHOWeightsOnFace()
      : FaceElementForcesAndSourcesCore::UserDataOperator(NOSPACE) {}
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

/** \brief Apply contravariant (Piola) transfer to Hdiv space for HO geometry

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
    : public FaceElementForcesAndSourcesCore::UserDataOperator {

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

template <int DIM> struct OpCalculateHOJac;

template <> struct OpCalculateHOJac<3> : public OpCalculateHOJacForVolume {
  using OpCalculateHOJacForVolume::OpCalculateHOJacForVolume;
};

template <> struct OpCalculateHOJac<2> : public OpCalculateHOJacForFaceImpl<2> {
  using OpCalculateHOJacForFaceImpl<2>::OpCalculateHOJacForFaceImpl;
};

/** \brief Calculate normals at Gauss points of triangle element
 * \ingroup mofem_forces_and_source
 */
template <int FIELD_DIM = 3>
struct OpGetHONormalsOnFace
    : public FaceElementForcesAndSourcesCore::UserDataOperator {

  OpGetHONormalsOnFace(std::string field_name)
      : FaceElementForcesAndSourcesCore::UserDataOperator(field_name, OPROW) {}

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
    : public FaceElementForcesAndSourcesCore::UserDataOperator {

  OpHOSetContravariantPiolaTransformOnFace3D(
      const FieldSpace space,
      boost::shared_ptr<MatrixDouble> normals_at_gauss_pts = nullptr)
      : FaceElementForcesAndSourcesCore::UserDataOperator(space, OPSPACE),
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
    : public EdgeElementForcesAndSourcesCore::UserDataOperator {

  OpHOSetContravariantPiolaTransformOnEdge3D(
      const FieldSpace space = HCURL,
      boost::shared_ptr<MatrixDouble> tangent_at_pts = nullptr)
      : EdgeElementForcesAndSourcesCore::UserDataOperator(space),
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
    : public FaceElementForcesAndSourcesCore::UserDataOperator {

  OpHOSetCovariantPiolaTransformOnFace3D(
      const FieldSpace space,
      boost::shared_ptr<MatrixDouble> normals_at_pts = nullptr,
      boost::shared_ptr<MatrixDouble> tangent1_at_pts = nullptr,
      boost::shared_ptr<MatrixDouble> tangent2_at_pts = nullptr)
      : FaceElementForcesAndSourcesCore::UserDataOperator(space, OPSPACE),
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
template <int FIELD_DIM = 3>
struct OpGetHOTangentsOnEdge
    : public EdgeElementForcesAndSourcesCore::UserDataOperator {

  OpGetHOTangentsOnEdge(
      std::string field_name,
      boost::shared_ptr<MatrixDouble> tangents_at_pts = nullptr)
      : EdgeElementForcesAndSourcesCore::UserDataOperator(field_name, OPROW),
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
        SETERRQ(PETSC_COMM_SELF, MOFEM_IMPOSSIBLE_CASE, "impossible case");
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
        SETERRQ(PETSC_COMM_SELF, MOFEM_IMPOSSIBLE_CASE, "impossible case");
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
        SETERRQ(PETSC_COMM_SELF, MOFEM_IMPOSSIBLE_CASE, "impossible case");
      }
    } else {
      SETERRQ(PETSC_COMM_SELF, MOFEM_IMPOSSIBLE_CASE, "impossible case");
    }

    MoFEMFunctionReturn(0);
  }

private:
  FieldSpace fieldSpace;
  boost::shared_ptr<VectorDouble> detJacPtr;
};

/**
 * @brief Add operators pushing bases from local to physical configuration
 *
 * @tparam FE_DIM dimension of element
 * @tparam PROBLEM_DIM problem dimension
 * @tparam SPACE_DIM space dimension
 */
template <int FE_DIM, int PROBLEM_DIM, int SPACE_DIM> struct AddHOOps;

template <> struct AddHOOps<2, 2, 2> {
  AddHOOps() = delete;
  static MoFEMErrorCode
  add(boost::ptr_deque<ForcesAndSourcesCore::UserDataOperator> &pipeline,
      std::vector<FieldSpace> spaces, std::string geom_field_name = "",
      boost::shared_ptr<MatrixDouble> jac = nullptr,
      boost::shared_ptr<VectorDouble> det = nullptr,
      boost::shared_ptr<MatrixDouble> inv_jac = nullptr);
};

template <> struct AddHOOps<1, 2, 2> {
  AddHOOps() = delete;
  static MoFEMErrorCode
  add(boost::ptr_deque<ForcesAndSourcesCore::UserDataOperator> &pipeline,
      std::vector<FieldSpace> spaces, std::string geom_field_name = "");
};
template <> struct AddHOOps<1, 3, 3> {
  AddHOOps() = delete;
  static MoFEMErrorCode
  add(boost::ptr_deque<ForcesAndSourcesCore::UserDataOperator> &pipeline,
      std::vector<FieldSpace> space, std::string geom_field_name = "");
};

template <> struct AddHOOps<2, 3, 3> {
  AddHOOps() = delete;
  static MoFEMErrorCode
  add(boost::ptr_deque<ForcesAndSourcesCore::UserDataOperator> &pipeline,
      std::vector<FieldSpace> space, std::string geom_field_name = "");
};

template <> struct AddHOOps<3, 3, 3> {
  AddHOOps() = delete;
  static MoFEMErrorCode
  add(boost::ptr_deque<ForcesAndSourcesCore::UserDataOperator> &pipeline,
      std::vector<FieldSpace> space, std::string geom_field_name = "",
      boost::shared_ptr<MatrixDouble> jac = nullptr,
      boost::shared_ptr<VectorDouble> det = nullptr,
      boost::shared_ptr<MatrixDouble> inv_jac = nullptr);
};

template <int FIELD_DIM>
MoFEMErrorCode
OpCalculateHOCoords<FIELD_DIM>::doWork(int side, EntityType type,
                                       EntitiesFieldData::EntData &data) {
  FTensor::Index<'i', FIELD_DIM> i;
  MoFEMFunctionBegin;
  const auto nb_dofs = data.getFieldData().size() / FIELD_DIM;
  if (nb_dofs) {
    if (type == MBVERTEX)
      getCoordsAtGaussPts().clear();
    auto t_base = data.getFTensor0N();
    auto t_coords = getFTensor1CoordsAtGaussPts();
    const auto nb_integration_pts = data.getN().size1();
    const auto nb_base_functions = data.getN().size2();
    for (auto gg = 0; gg != nb_integration_pts; ++gg) {
      auto t_dof = data.getFTensor1FieldData<FIELD_DIM>();
      size_t bb = 0;
      for (; bb != nb_dofs; ++bb) {
        t_coords(i) += t_base * t_dof(i);
        ++t_dof;
        ++t_base;
      }
      for (; bb < nb_base_functions; ++bb)
        ++t_base;
      ++t_coords;
    }
  }
  MoFEMFunctionReturn(0);
}

template <int FIELD_DIM>
MoFEMErrorCode
OpGetHONormalsOnFace<FIELD_DIM>::doWork(int side, EntityType type,
                                        EntitiesFieldData::EntData &data) {
  MoFEMFunctionBeginHot;

  FTensor::Number<0> N0;
  FTensor::Number<1> N1;

  auto get_ftensor1 = [](MatrixDouble &m) {
    return FTensor::Tensor1<FTensor::PackPtr<double *, 3>, 3>(
        &m(0, 0), &m(0, 1), &m(0, 2));
  };

  unsigned int nb_dofs = data.getFieldData().size();
  if (nb_dofs != 0) {

    int nb_gauss_pts = data.getN().size1();
    auto &tangent1_at_gauss_pts = getTangent1AtGaussPts();
    auto &tangent2_at_gauss_pts = getTangent2AtGaussPts();
    tangent1_at_gauss_pts.resize(nb_gauss_pts, 3, false);
    tangent2_at_gauss_pts.resize(nb_gauss_pts, 3, false);

    switch (type) {
    case MBVERTEX: {
      tangent1_at_gauss_pts.clear();
      tangent2_at_gauss_pts.clear();
    }
    case MBEDGE:
    case MBTRI:
    case MBQUAD: {

#ifndef NDEBUG
      if (2 * data.getN().size2() != data.getDiffN().size2()) {
        SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
                "expected two direcatives in local element coordinates");
      }
      if (nb_dofs % FIELD_DIM != 0) {
        SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
                "expected that number of dofs is multiplicative of field "
                "dimension");
      }
#endif

      if (nb_dofs > FIELD_DIM * data.getN().size2()) {
        unsigned int nn = 0;
        for (; nn != nb_dofs; nn++) {
          if (!data.getFieldDofs()[nn]->getActive())
            break;
        }
        if (nn > FIELD_DIM * data.getN().size2()) {
          SETERRQ1(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
                   "Data inconsistency for base %s",
                   ApproximationBaseNames[data.getBase()]);
        } else {
          nb_dofs = nn;
          if (!nb_dofs)
            MoFEMFunctionReturnHot(0);
        }
      }
      const int nb_base_functions = data.getN().size2();
      auto t_base = data.getFTensor0N();
      auto t_diff_base = data.getFTensor1DiffN<2>();
      auto t_t1 = get_ftensor1(tangent1_at_gauss_pts);
      auto t_t2 = get_ftensor1(tangent2_at_gauss_pts);
      for (int gg = 0; gg != nb_gauss_pts; ++gg) {
        auto t_data = data.getFTensor1FieldData<FIELD_DIM>();
        int bb = 0;
        for (; bb != nb_dofs / FIELD_DIM; ++bb) {
          FTensor::Index<'i', FIELD_DIM> i;
          t_t1(i) += t_data(i) * t_diff_base(N0);
          t_t2(i) += t_data(i) * t_diff_base(N1);
          ++t_data;
          ++t_base;
          ++t_diff_base;
        }
        for (; bb != nb_base_functions; ++bb) {
          ++t_base;
          ++t_diff_base;
        }
        ++t_t1;
        ++t_t2;
      }
    } break;
    default:
      SETERRQ(PETSC_COMM_SELF, MOFEM_NOT_IMPLEMENTED, "not implemented");
    }
  }

  if (moab::CN::Dimension(type) == 2) {

    auto &normal_at_gauss_pts = getNormalsAtGaussPts();
    auto &tangent1_at_gauss_pts = getTangent1AtGaussPts();
    auto &tangent2_at_gauss_pts = getTangent2AtGaussPts();

    const auto nb_gauss_pts = tangent1_at_gauss_pts.size1();
    normal_at_gauss_pts.resize(nb_gauss_pts, 3, false);

    auto t_normal = get_ftensor1(normal_at_gauss_pts);
    auto t_t1 = get_ftensor1(tangent1_at_gauss_pts);
    auto t_t2 = get_ftensor1(tangent2_at_gauss_pts);
    for (unsigned int gg = 0; gg != nb_gauss_pts; ++gg) {
      FTensor::Index<'i', 3> i;
      FTensor::Index<'j', 3> j;
      FTensor::Index<'k', 3> k;
      t_normal(j) = FTensor::levi_civita(i, j, k) * t_t1(k) * t_t2(i);
      ++t_normal;
      ++t_t1;
      ++t_t2;
    }
  }

  MoFEMFunctionReturnHot(0);
}

template <int FIELD_DIM>
MoFEMErrorCode
OpGetHOTangentsOnEdge<FIELD_DIM>::doWork(int side, EntityType type,
                                         EntitiesFieldData::EntData &data) {
  MoFEMFunctionBegin;

  auto get_tangent = [&]() -> MatrixDouble & {
    if (tangentsAtPts)
      return *tangentsAtPts;
    else
      return getTangentAtGaussPts();
  };

  auto &tangent = get_tangent();
  int nb_gauss_pts = getGaussPts().size2();

  if (type == MBVERTEX) {
    tangent.resize(nb_gauss_pts, 3, false);
    tangent.clear();
  }

  int nb_dofs = data.getFieldData().size();
  if (nb_dofs != 0) {
    const int nb_base_functions = data.getN().size2();
    double *diff_base_function = &*data.getDiffN().data().begin();
    auto tangent_at_gauss_pts =
        getFTensor1FromPtr<FIELD_DIM, 3>(&*tangent.data().begin());

    FTensor::Index<'i', FIELD_DIM> i;
    int size = nb_dofs / FIELD_DIM;
    if (nb_dofs % FIELD_DIM) {
      SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
              "Inconsistent number of dofs and filed dimension");
    }
    for (int gg = 0; gg != nb_gauss_pts; ++gg) {
      auto field_data = data.getFTensor1FieldData<FIELD_DIM>();
      int bb = 0;
      for (; bb < size; ++bb) {
        tangent_at_gauss_pts(i) += field_data(i) * (*diff_base_function);
        ++field_data;
        ++diff_base_function;
      }
      for (; bb != nb_base_functions; ++bb) {
        ++diff_base_function;
      }
      ++tangent_at_gauss_pts;
    }
  }
  MoFEMFunctionReturn(0);
}

/**
 * @deprecated do not use this function, instead use AddHOOps.
 *
 * @tparam E
 * @param field
 * @param e
 * @param h1
 * @param hcurl
 * @param hdiv
 * @param l2
 * @return MoFEMErrorCode
 */
template <typename E>
MoFEMErrorCode addHOOpsVol(const std::string field, E &e, bool h1,
                                      bool hcurl, bool hdiv, bool l2) {
  std::vector<FieldSpace> spaces;
  if (h1)
    spaces.push_back(H1);
  if (hcurl)
    spaces.push_back(HCURL);
  if (hdiv)
    spaces.push_back(HDIV);
  if (l2)
    spaces.push_back(L2);
  return AddHOOps<3, 3, 3>::add(e.getOpPtrVector(), spaces, field);
}

/**
 * @deprecated do not use this function, instead use AddHOOps.
 *
 * @tparam E
 * @param field
 * @param e
 * @param hcurl
 * @param hdiv
 * @return MoFEMErrorCode
 */
template <typename E>
MoFEMErrorCode addHOOpsFace3D(const std::string field, E &e,
                                         bool hcurl, bool hdiv) {
  std::vector<FieldSpace> spaces;
  if (hcurl)
    spaces.push_back(HCURL);
  if (hdiv)
    spaces.push_back(HDIV);
  return AddHOOps<2, 3, 3>::add(e.getOpPtrVector(), spaces, field);
}

}; // namespace MoFEM

#endif