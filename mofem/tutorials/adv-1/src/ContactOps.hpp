

/**
 * \file ContactOps.hpp
 * \example ContactOps.hpp
 */

#ifndef __CONTACTOPS_HPP__
#define __CONTACTOPS_HPP__

namespace ContactOps {

//! [Common data]
struct CommonData : public boost::enable_shared_from_this<CommonData> {
  // MatrixDouble contactStress;
  MatrixDouble contactTraction;
  MatrixDouble contactDisp;
  MatrixDouble contactDispGrad;

  VectorDouble sdfVals;  ///< size is equal to number of gauss points on element
  MatrixDouble gradsSdf; ///< nb of rows is equals to dimension, and nb of cols
                         ///< is equals to number of gauss points on element
  MatrixDouble hessSdf;  ///< nb of rows is equals to nb of element of symmetric
                        ///< matrix, and nb of cols is equals to number of gauss
                        ///< points on element
  VectorDouble constraintVals;

  static SmartPetscObj<Vec>
      totalTraction; // User have to release and create vector when appropiate.

  static auto createTotalTraction(MoFEM::Interface &m_field) {
    constexpr int ghosts[] = {0, 1, 2, 3, 4};
    totalTraction =
        createGhostVector(m_field.get_comm(),

                          (m_field.get_comm_rank() == 0) ? 5 : 0, 5,

                          (m_field.get_comm_rank() == 0) ? 0 : 5, ghosts);
    return totalTraction;
  }

  static auto getFTensor1TotalTraction() {
    if (CommonData::totalTraction) {
      const double *t_ptr;
      CHK_THROW_MESSAGE(VecGetArrayRead(CommonData::totalTraction, &t_ptr),
                        "get array");
      FTensor::Tensor1<double, 5> t{t_ptr[0], t_ptr[1], t_ptr[2], t_ptr[3],
                                    t_ptr[4]};
      CHK_THROW_MESSAGE(VecRestoreArrayRead(CommonData::totalTraction, &t_ptr),
                        "restore array");
      return t;
    } else {
      return FTensor::Tensor1<double, 5>{0., 0., 0., 0., 0.};
    }
  }

  inline auto contactTractionPtr() {
    return boost::shared_ptr<MatrixDouble>(shared_from_this(),
                                           &contactTraction);
  }

  inline auto contactDispPtr() {
    return boost::shared_ptr<MatrixDouble>(shared_from_this(), &contactDisp);
  }

  inline auto contactDispGradPtr() {
    return boost::shared_ptr<MatrixDouble>(shared_from_this(),
                                           &contactDispGrad);
  }

  inline auto sdfPtr() {
    return boost::shared_ptr<VectorDouble>(shared_from_this(), &sdfVals);
  }

  inline auto gradSdfPtr() {
    return boost::shared_ptr<MatrixDouble>(shared_from_this(), &gradsSdf);
  }

  inline auto hessSdfPtr() {
    return boost::shared_ptr<MatrixDouble>(shared_from_this(), &hessSdf);
  }

  inline auto constraintPtr() {
    return boost::shared_ptr<VectorDouble>(shared_from_this(), &constraintVals);
  }
};

SmartPetscObj<Vec> CommonData::totalTraction;

//! [Common data]

//! [Surface distance function from python]
#ifdef PYTHON_SDF
struct SDFPython {
  SDFPython() = default;
  virtual ~SDFPython() = default;

  MoFEMErrorCode sdfInit(const std::string py_file) {
    MoFEMFunctionBegin;
    try {

      // create main module
      auto main_module = bp::import("__main__");
      mainNamespace = main_module.attr("__dict__");
      bp::exec_file(py_file.c_str(), mainNamespace, mainNamespace);
      // create a reference to python function
      sdfFun = mainNamespace["sdf"];
      sdfGradFun = mainNamespace["grad_sdf"];
      sdfHessFun = mainNamespace["hess_sdf"];

    } catch (bp::error_already_set const &) {
      // print all other errors to stderr
      PyErr_Print();
      CHK_THROW_MESSAGE(MOFEM_OPERATION_UNSUCCESSFUL, "Python error");
    }
    MoFEMFunctionReturn(0);
  };

  template <typename T>
  inline std::vector<T>
  py_list_to_std_vector(const boost::python::object &iterable) {
    return std::vector<T>(boost::python::stl_input_iterator<T>(iterable),
                          boost::python::stl_input_iterator<T>());
  }

  MoFEMErrorCode evalSdf(

      double delta_t, double t, np::ndarray x, np::ndarray y, np::ndarray z,
      np::ndarray tx, np::ndarray ty, np::ndarray tz, int block_id,
      np::ndarray &sdf

  ) {
    MoFEMFunctionBegin;
    try {

      // call python function
      sdf = bp::extract<np::ndarray>(
          sdfFun(delta_t, t, x, y, z, tx, ty, tz, block_id));

    } catch (bp::error_already_set const &) {
      // print all other errors to stderr
      PyErr_Print();
      CHK_THROW_MESSAGE(MOFEM_OPERATION_UNSUCCESSFUL, "Python error");
    }
    MoFEMFunctionReturn(0);
  }

  MoFEMErrorCode evalGradSdf(

      double delta_t, double t, np::ndarray x, np::ndarray y, np::ndarray z,
      np::ndarray tx, np::ndarray ty, np::ndarray tz, int block_id,
      np::ndarray &grad_sdf

  ) {
    MoFEMFunctionBegin;
    try {

      // call python function
      grad_sdf = bp::extract<np::ndarray>(
          sdfGradFun(delta_t, t, x, y, z, tx, ty, tz, block_id));

    } catch (bp::error_already_set const &) {
      // print all other errors to stderr
      PyErr_Print();
      CHK_THROW_MESSAGE(MOFEM_OPERATION_UNSUCCESSFUL, "Python error");
    }
    MoFEMFunctionReturn(0);
  }

  MoFEMErrorCode evalHessSdf(

      double delta_t, double t, np::ndarray x, np::ndarray y, np::ndarray z,
      np::ndarray tx, np::ndarray ty, np::ndarray tz, int block_id,
      np::ndarray &hess_sdf

  ) {
    MoFEMFunctionBegin;
    try {

      // call python function
      hess_sdf = bp::extract<np::ndarray>(
          sdfHessFun(delta_t, t, x, y, z, tx, ty, tz, block_id));

    } catch (bp::error_already_set const &) {
      // print all other errors to stderr
      PyErr_Print();
      CHK_THROW_MESSAGE(MOFEM_OPERATION_UNSUCCESSFUL, "Python error");
    }
    MoFEMFunctionReturn(0);
  }

private:
  bp::object mainNamespace;
  bp::object sdfFun;
  bp::object sdfGradFun;
  bp::object sdfHessFun;
};

static boost::weak_ptr<SDFPython> sdfPythonWeakPtr;

inline np::ndarray convert_to_numpy(VectorDouble &data, int nb_gauss_pts,
                                    int id) {
  auto dtype = np::dtype::get_builtin<double>();
  auto size = bp::make_tuple(nb_gauss_pts);
  auto stride = bp::make_tuple(3 * sizeof(double));
  return (np::from_data(&data[id], dtype, size, stride, bp::object()));
};
#endif
//! [Surface distance function from python]

using SurfaceDistanceFunction = boost::function<VectorDouble(
    double delta_t, double t, int nb_gauss_pts, MatrixDouble &spatial_coords,
    MatrixDouble &normals_at_pts, int block_id)>;

using GradSurfaceDistanceFunction = boost::function<MatrixDouble(
    double delta_t, double t, int nb_gauss_pts, MatrixDouble &spatial_coords,
    MatrixDouble &normals_at_pts, int block_id)>;

using HessSurfaceDistanceFunction = boost::function<MatrixDouble(
    double delta_t, double t, int nb_gauss_pts, MatrixDouble &spatial_coords,
    MatrixDouble &normals_at_pts, int block_id)>;

inline VectorDouble surface_distance_function(double delta_t, double t,
                                              int nb_gauss_pts,
                                              MatrixDouble &m_spatial_coords,
                                              MatrixDouble &m_normals_at_pts,
                                              int block_id) {

#ifdef PYTHON_SDF
  if (auto sdf_ptr = sdfPythonWeakPtr.lock()) {

    VectorDouble v_spatial_coords = m_spatial_coords.data();
    VectorDouble v_normal_at_pts = m_normals_at_pts.data();

    bp::list python_coords;
    bp::list python_normals;

    for (int idx = 0; idx < 3; ++idx) {
      python_coords.append(
          convert_to_numpy(v_spatial_coords, nb_gauss_pts, idx));
      python_normals.append(
          convert_to_numpy(v_normal_at_pts, nb_gauss_pts, idx));
    }

    np::ndarray np_sdf = np::empty(bp::make_tuple(nb_gauss_pts),
                                   np::dtype::get_builtin<double>());
    CHK_MOAB_THROW(sdf_ptr->evalSdf(delta_t, t,
                                    bp::extract<np::ndarray>(python_coords[0]),
                                    bp::extract<np::ndarray>(python_coords[1]),
                                    bp::extract<np::ndarray>(python_coords[2]),
                                    bp::extract<np::ndarray>(python_normals[0]),
                                    bp::extract<np::ndarray>(python_normals[1]),
                                    bp::extract<np::ndarray>(python_normals[2]),
                                    block_id, np_sdf),
                   "Failed python call");

    double *sdf_val_ptr = reinterpret_cast<double *>(np_sdf.get_data());

    VectorDouble v_sdf;
    v_sdf.resize(nb_gauss_pts, false);

    for (size_t gg = 0; gg < nb_gauss_pts; ++gg)
      v_sdf[gg] = *(sdf_val_ptr + gg);

    return v_sdf;
  }
#endif
  VectorDouble v_sdf;
  v_sdf.resize(nb_gauss_pts, false);
  auto t_coords = getFTensor1FromPtr<3>(&m_spatial_coords(0, 0));

  for (size_t gg = 0; gg < nb_gauss_pts; ++gg) {
    v_sdf[gg] = -t_coords(2) - 0.1;
    ++t_coords;
  }

  return v_sdf;
}

inline MatrixDouble
grad_surface_distance_function(double delta_t, double t, int nb_gauss_pts,
                               MatrixDouble &m_spatial_coords,
                               MatrixDouble &m_normals_at_pts, int block_id) {
#ifdef PYTHON_SDF
  if (auto sdf_ptr = sdfPythonWeakPtr.lock()) {

    VectorDouble v_spatial_coords = m_spatial_coords.data();
    VectorDouble v_normal_at_pts = m_normals_at_pts.data();

    bp::list python_coords;
    bp::list python_normals;

    for (int idx = 0; idx < 3; ++idx) {
      python_coords.append(
          convert_to_numpy(v_spatial_coords, nb_gauss_pts, idx));
      python_normals.append(
          convert_to_numpy(v_normal_at_pts, nb_gauss_pts, idx));
    }

    np::ndarray np_grad_sdf = np::empty(bp::make_tuple(nb_gauss_pts, 3),
                                        np::dtype::get_builtin<double>());
    CHK_MOAB_THROW(sdf_ptr->evalGradSdf(
                       delta_t, t, bp::extract<np::ndarray>(python_coords[0]),
                       bp::extract<np::ndarray>(python_coords[1]),
                       bp::extract<np::ndarray>(python_coords[2]),
                       bp::extract<np::ndarray>(python_normals[0]),
                       bp::extract<np::ndarray>(python_normals[1]),
                       bp::extract<np::ndarray>(python_normals[2]), block_id,
                       np_grad_sdf),
                   "Failed python call");

    double *grad_ptr = reinterpret_cast<double *>(np_grad_sdf.get_data());

    MatrixDouble m_grad_sdf;
    m_grad_sdf.resize(3, nb_gauss_pts, false);
    for (size_t gg = 0; gg < nb_gauss_pts; ++gg) {
      for (int idx = 0; idx < 3; ++idx)
        m_grad_sdf(idx, gg) = *(grad_ptr + (3 * gg + idx));
    }
    return m_grad_sdf;
  }
#endif
  MatrixDouble m_grad_sdf;
  m_grad_sdf.resize(3, nb_gauss_pts, false);
  FTensor::Index<'i', 3> i;
  FTensor::Tensor1<double, 3> t_grad_sdf_set{0.0, 0.0, -1.0};
  auto t_grad_sdf = getFTensor1FromMat<3>(m_grad_sdf);

  for (size_t gg = 0; gg < nb_gauss_pts; ++gg) {
    t_grad_sdf(i) = t_grad_sdf_set(i);
    ++t_grad_sdf;
  }

  return m_grad_sdf;
}

inline MatrixDouble
hess_surface_distance_function(double delta_t, double t, int nb_gauss_pts,
                               MatrixDouble &m_spatial_coords,
                               MatrixDouble &m_normals_at_pts, int block_id) {
#ifdef PYTHON_SDF
  if (auto sdf_ptr = sdfPythonWeakPtr.lock()) {

    VectorDouble v_spatial_coords = m_spatial_coords.data();
    VectorDouble v_normal_at_pts = m_normals_at_pts.data();

    bp::list python_coords;
    bp::list python_normals;

    for (int idx = 0; idx < 3; ++idx) {
      python_coords.append(
          convert_to_numpy(v_spatial_coords, nb_gauss_pts, idx));
      python_normals.append(
          convert_to_numpy(v_normal_at_pts, nb_gauss_pts, idx));
    };

    np::ndarray np_hess_sdf = np::empty(bp::make_tuple(nb_gauss_pts, 6),
                                        np::dtype::get_builtin<double>());
    CHK_MOAB_THROW(sdf_ptr->evalHessSdf(
                       delta_t, t, bp::extract<np::ndarray>(python_coords[0]),
                       bp::extract<np::ndarray>(python_coords[1]),
                       bp::extract<np::ndarray>(python_coords[2]),
                       bp::extract<np::ndarray>(python_normals[0]),
                       bp::extract<np::ndarray>(python_normals[1]),
                       bp::extract<np::ndarray>(python_normals[2]), block_id,
                       np_hess_sdf),
                   "Failed python call");

    double *hess_ptr = reinterpret_cast<double *>(np_hess_sdf.get_data());

    MatrixDouble m_hess_sdf;
    m_hess_sdf.resize(6, nb_gauss_pts, false);
    for (size_t gg = 0; gg < nb_gauss_pts; ++gg) {
      for (int idx = 0; idx < 6; ++idx)
        m_hess_sdf(idx, gg) =
            *(hess_ptr + (6 * gg + idx));
    }
    return m_hess_sdf;
  }
#endif
  MatrixDouble m_hess_sdf;
  m_hess_sdf.resize(6, nb_gauss_pts, false);
  FTensor::Index<'i', 3> i;
  FTensor::Index<'j', 3> j;
  FTensor::Tensor2_symmetric<double, 3> t_hess_sdf_set{0., 0., 0., 0., 0., 0.};
  auto t_hess_sdf = getFTensor2SymmetricFromMat<3>(m_hess_sdf);

  for (size_t gg = 0; gg < nb_gauss_pts; ++gg) {
    t_hess_sdf(i, j) = t_hess_sdf_set(i, j);
    ++t_hess_sdf;
  }
  return m_hess_sdf;
}

template <int DIM, IntegrationType I, typename BoundaryEleOp>
struct OpAssembleTotalContactTractionImpl;

template <int DIM, IntegrationType I, typename BoundaryEleOp>
struct OpAssembleTotalContactAreaImpl;

template <int DIM, IntegrationType I, typename BoundaryEleOp>
struct OpEvaluateSDFImpl;

template <int DIM, IntegrationType I, typename AssemblyBoundaryEleOp>
struct OpConstrainBoundaryRhsImpl;

template <int DIM, IntegrationType I, typename AssemblyBoundaryEleOp>
struct OpConstrainBoundaryLhs_dUImpl;

template <int DIM, IntegrationType I, typename AssemblyBoundaryEleOp>
struct OpConstrainBoundaryLhs_dTractionImpl;

template <typename T1, typename T2, int DIM1, int DIM2>
inline auto get_spatial_coords(FTensor::Tensor1<T1, DIM1> &&t_coords,
                               FTensor::Tensor1<T2, DIM2> &&t_disp,
                               size_t nb_gauss_pts) {
  MatrixDouble m_spatial_coords(nb_gauss_pts, 3);
  m_spatial_coords.clear();
  auto t_spatial_coords = getFTensor1FromPtr<3>(&m_spatial_coords(0, 0));
  FTensor::Index<'i', DIM2> i;
  for (auto gg = 0; gg != nb_gauss_pts; ++gg) {
    t_spatial_coords(i) = t_coords(i) + t_disp(i);
    ++t_spatial_coords;
    ++t_coords;
    ++t_disp;
  }
  return m_spatial_coords;
}

template <typename T1, int DIM1>
inline auto get_normalize_normals(FTensor::Tensor1<T1, DIM1> &&t_normal_at_pts,
                                  size_t nb_gauss_pts) {
  MatrixDouble m_normals_at_pts(3, nb_gauss_pts);
  m_normals_at_pts.clear();
  FTensor::Index<'i', DIM1> i;
  auto t_set_normal = getFTensor1FromMat<3>(m_normals_at_pts);
  for (auto gg = 0; gg != nb_gauss_pts; ++gg) {
    t_set_normal(i) = t_normal_at_pts(i) / t_normal_at_pts.l2();
    ++t_set_normal;
    ++t_normal_at_pts;
  }
  return m_normals_at_pts;
}

template <int DIM, typename BoundaryEleOp>
struct OpAssembleTotalContactTractionImpl<DIM, GAUSS, BoundaryEleOp>
    : public BoundaryEleOp {
  OpAssembleTotalContactTractionImpl(
      boost::shared_ptr<CommonData> common_data_ptr, double scale = 1,
      bool is_axisymmetric = false);
  MoFEMErrorCode doWork(int side, EntityType type, EntData &data);

private:
  boost::shared_ptr<CommonData> commonDataPtr;
  const double scaleTraction;
  bool isAxisymmetric;
};

template <int DIM, typename BoundaryEleOp>
struct OpAssembleTotalContactAreaImpl<DIM, GAUSS, BoundaryEleOp>
    : public BoundaryEleOp {
  OpAssembleTotalContactAreaImpl(
      boost::shared_ptr<CommonData> common_data_ptr,
      bool is_axisymmetric = false,
      boost::shared_ptr<Range> contact_range_ptr = nullptr);
  MoFEMErrorCode doWork(int side, EntityType type, EntData &data);

  SurfaceDistanceFunction surfaceDistanceFunction = surface_distance_function;
  GradSurfaceDistanceFunction gradSurfaceDistanceFunction =
      grad_surface_distance_function;

private:
  boost::shared_ptr<CommonData> commonDataPtr;
  bool isAxisymmetric;
  boost::shared_ptr<Range> contactRange;
};

template <int DIM, typename BoundaryEleOp>
struct OpEvaluateSDFImpl<DIM, GAUSS, BoundaryEleOp> : public BoundaryEleOp {
  OpEvaluateSDFImpl(boost::shared_ptr<CommonData> common_data_ptr);
  MoFEMErrorCode doWork(int side, EntityType type, EntData &data);

private:
  boost::shared_ptr<CommonData> commonDataPtr;

  SurfaceDistanceFunction surfaceDistanceFunction = surface_distance_function;
  GradSurfaceDistanceFunction gradSurfaceDistanceFunction =
      grad_surface_distance_function;
  HessSurfaceDistanceFunction hessSurfaceDistanceFunction =
      hess_surface_distance_function;
};

template <int DIM, typename AssemblyBoundaryEleOp>
struct OpConstrainBoundaryRhsImpl<DIM, GAUSS, AssemblyBoundaryEleOp>
    : public AssemblyBoundaryEleOp {
  OpConstrainBoundaryRhsImpl(const std::string field_name,
                             boost::shared_ptr<CommonData> common_data_ptr,
                             bool is_axisymmetric = false);
  MoFEMErrorCode iNtegrate(EntitiesFieldData::EntData &data);

  SurfaceDistanceFunction surfaceDistanceFunction = surface_distance_function;
  GradSurfaceDistanceFunction gradSurfaceDistanceFunction =
      grad_surface_distance_function;

private:
  boost::shared_ptr<CommonData> commonDataPtr;
  bool isAxisymmetric;
};

template <int DIM, typename AssemblyBoundaryEleOp>
struct OpConstrainBoundaryLhs_dUImpl<DIM, GAUSS, AssemblyBoundaryEleOp>
    : public AssemblyBoundaryEleOp {
  OpConstrainBoundaryLhs_dUImpl(const std::string row_field_name,
                                const std::string col_field_name,
                                boost::shared_ptr<CommonData> common_data_ptr,
                                bool is_axisymmetric = false);
  MoFEMErrorCode iNtegrate(EntitiesFieldData::EntData &row_data,
                           EntitiesFieldData::EntData &col_data);

  SurfaceDistanceFunction surfaceDistanceFunction = surface_distance_function;
  GradSurfaceDistanceFunction gradSurfaceDistanceFunction =
      grad_surface_distance_function;
  HessSurfaceDistanceFunction hessSurfaceDistanceFunction =
      hess_surface_distance_function;

  boost::shared_ptr<CommonData> commonDataPtr;
  bool isAxisymmetric;
};

template <int DIM, typename AssemblyBoundaryEleOp>
struct OpConstrainBoundaryLhs_dTractionImpl<DIM, GAUSS, AssemblyBoundaryEleOp>
    : public AssemblyBoundaryEleOp {
  OpConstrainBoundaryLhs_dTractionImpl(
      const std::string row_field_name, const std::string col_field_name,
      boost::shared_ptr<CommonData> common_data_ptr,
      bool is_axisymmetric = false);
  MoFEMErrorCode iNtegrate(EntitiesFieldData::EntData &row_data,
                           EntitiesFieldData::EntData &col_data);

  SurfaceDistanceFunction surfaceDistanceFunction = surface_distance_function;
  GradSurfaceDistanceFunction gradSurfaceDistanceFunction =
      grad_surface_distance_function;

private:
  boost::shared_ptr<CommonData> commonDataPtr;
  bool isAxisymmetric;
};

template <typename BoundaryEleOp> struct ContactIntegrators {
  template <int DIM, IntegrationType I>
  using OpAssembleTotalContactTraction =
      OpAssembleTotalContactTractionImpl<DIM, I, BoundaryEleOp>;

  template <int DIM, IntegrationType I>
  using OpAssembleTotalContactArea =
      OpAssembleTotalContactAreaImpl<DIM, I, BoundaryEleOp>;

  template <int DIM, IntegrationType I>
  using OpEvaluateSDF = OpEvaluateSDFImpl<DIM, I, BoundaryEleOp>;

  template <AssemblyType A> struct Assembly {

    using AssemblyBoundaryEleOp =
        typename FormsIntegrators<BoundaryEleOp>::template Assembly<A>::OpBase;

    template <int DIM, IntegrationType I>
    using OpConstrainBoundaryRhs =
        OpConstrainBoundaryRhsImpl<DIM, I, AssemblyBoundaryEleOp>;

    template <int DIM, IntegrationType I>
    using OpConstrainBoundaryLhs_dU =
        OpConstrainBoundaryLhs_dUImpl<DIM, I, AssemblyBoundaryEleOp>;

    template <int DIM, IntegrationType I>
    using OpConstrainBoundaryLhs_dTraction =
        OpConstrainBoundaryLhs_dTractionImpl<DIM, I, AssemblyBoundaryEleOp>;
  };
};

inline double sign(double x) {
  constexpr auto eps = std::numeric_limits<float>::epsilon();
  if (std::abs(x) < eps)
    return 0;
  else if (x > eps)
    return 1;
  else
    return -1;
};

inline double w(const double sdf, const double tn) {
  return sdf - cn_contact * tn;
}

/**
 * @brief constrain function
 *
 * return 1 if negative sdf or positive tn
 *
 * @param sdf signed distance
 * @param tn traction
 * @return double
 */
inline double constrain(double sdf, double tn) {
  const auto s = sign(w(sdf, tn));
  return (1 - s) / 2;
}

template <int DIM, typename BoundaryEleOp>
OpAssembleTotalContactTractionImpl<DIM, GAUSS, BoundaryEleOp>::
    OpAssembleTotalContactTractionImpl(
        boost::shared_ptr<CommonData> common_data_ptr, double scale,
        bool is_axisymmetric)
    : BoundaryEleOp(NOSPACE, BoundaryEleOp::OPSPACE),
      commonDataPtr(common_data_ptr), scaleTraction(scale),
      isAxisymmetric(is_axisymmetric) {}

template <int DIM, typename BoundaryEleOp>
MoFEMErrorCode
OpAssembleTotalContactTractionImpl<DIM, GAUSS, BoundaryEleOp>::doWork(
    int side, EntityType type, EntData &data) {
  MoFEMFunctionBegin;

  FTensor::Index<'i', DIM> i;
  FTensor::Tensor1<double, 3> t_sum_t{0., 0., 0.};

  auto t_w = BoundaryEleOp::getFTensor0IntegrationWeight();
  auto t_traction = getFTensor1FromMat<DIM>(commonDataPtr->contactTraction);
  auto t_coords = BoundaryEleOp::getFTensor1CoordsAtGaussPts();

  const auto nb_gauss_pts = BoundaryEleOp::getGaussPts().size2();
  for (auto gg = 0; gg != nb_gauss_pts; ++gg) {
    double jacobian = 1.;
    if (isAxisymmetric) {
      jacobian = 2. * M_PI * t_coords(0);
    }
    const double alpha = t_w * jacobian * BoundaryEleOp::getMeasure();
    t_sum_t(i) += alpha * t_traction(i);
    ++t_w;
    ++t_traction;
    ++t_coords;
  }

  t_sum_t(i) *= scaleTraction;

  constexpr int ind[] = {0, 1, 2};
  CHKERR VecSetValues(commonDataPtr->totalTraction, 3, ind, &t_sum_t(0),
                      ADD_VALUES);

  MoFEMFunctionReturn(0);
}
template <int DIM, typename BoundaryEleOp>
OpAssembleTotalContactAreaImpl<DIM, GAUSS, BoundaryEleOp>::
    OpAssembleTotalContactAreaImpl(
        boost::shared_ptr<CommonData> common_data_ptr, bool is_axisymmetric,
        boost::shared_ptr<Range> contact_range_ptr)
    : BoundaryEleOp(NOSPACE, BoundaryEleOp::OPSPACE),
      commonDataPtr(common_data_ptr), isAxisymmetric(is_axisymmetric),
      contactRange(contact_range_ptr) {}

template <int DIM, typename BoundaryEleOp>
MoFEMErrorCode
OpAssembleTotalContactAreaImpl<DIM, GAUSS, BoundaryEleOp>::doWork(
    int side, EntityType type, EntData &data) {
  MoFEMFunctionBegin;

  auto fe_type = BoundaryEleOp::getFEType();

  const auto fe_ent = BoundaryEleOp::getFEEntityHandle();

  if (contactRange->find(fe_ent) != contactRange->end()) {
    FTensor::Index<'i', DIM> i;
    FTensor::Index<'j', DIM> j;
    FTensor::Tensor1<double, 2> t_sum_a{0., 0.};

    auto t_w = BoundaryEleOp::getFTensor0IntegrationWeight();
    auto t_traction = getFTensor1FromMat<DIM>(commonDataPtr->contactTraction);
    auto t_coords = BoundaryEleOp::getFTensor1CoordsAtGaussPts();

    auto t_grad = getFTensor2FromMat<DIM, DIM>(commonDataPtr->contactDispGrad);
    auto t_normal_at_pts = BoundaryEleOp::getFTensor1NormalsAtGaussPts();

    const auto nb_gauss_pts = BoundaryEleOp::getGaussPts().size2();
    auto m_spatial_coords = get_spatial_coords(
        BoundaryEleOp::getFTensor1CoordsAtGaussPts(),
        getFTensor1FromMat<DIM>(commonDataPtr->contactDisp), nb_gauss_pts);
    auto m_normals_at_pts = get_normalize_normals(
        BoundaryEleOp::getFTensor1NormalsAtGaussPts(), nb_gauss_pts);

    auto t_normal = getFTensor1FromMat<3>(m_normals_at_pts);
    auto ts_time = BoundaryEleOp::getTStime();
    auto ts_time_step = BoundaryEleOp::getTStimeStep();
    int block_id = 0;
    auto v_sdf =
        surfaceDistanceFunction(ts_time_step, ts_time, nb_gauss_pts,
                                m_spatial_coords, m_normals_at_pts, block_id);
    auto m_grad_sdf = gradSurfaceDistanceFunction(
        ts_time_step, ts_time, nb_gauss_pts, m_spatial_coords, m_normals_at_pts,
        block_id);
    auto t_sdf = getFTensor0FromVec(v_sdf);
    auto t_grad_sdf = getFTensor1FromMat<3>(m_grad_sdf);
    for (auto gg = 0; gg != nb_gauss_pts; ++gg) {
      double jacobian = 1.;
      if (isAxisymmetric) {
        jacobian = 2. * M_PI * t_coords(0); // Axisymmetric Jacobian
      }
      auto tn = -t_traction(i) * t_grad_sdf(i);
      auto c = constrain(t_sdf, tn);
      double alpha = t_w * jacobian;

      FTensor::Tensor2<double, DIM, DIM> F;
      FTensor::Tensor2<double, DIM, DIM> invF;
      FTensor::Tensor1<double, DIM> t_normal_current;

      F(i, j) = t_grad(i, j) + kronecker_delta(i, j);
      auto det = determinantTensor(F);
      CHKERR invertTensor(F, det, invF);
      t_normal_current(i) = det * (invF(j, i) * t_normal_at_pts(j));

      alpha *= sqrt(t_normal_current(i) * t_normal_current(i));

      if (fe_type == MBTRI) {
        alpha /= 2;
      }
      if (c > 1e-12) {
        t_sum_a(0) += alpha; // real area
      }
      t_sum_a(1) += alpha; // Potential area
      ++t_w;
      ++t_traction;
      ++t_coords;
      ++t_sdf;
      ++t_grad_sdf;

      ++t_grad;
      ++t_normal_at_pts;
    }
    constexpr int ind[] = {3, 4};
    CHKERR VecSetValues(commonDataPtr->totalTraction, 2, ind, &t_sum_a(0),
                        ADD_VALUES);
  }
  MoFEMFunctionReturn(0);
}

template <int DIM, typename BoundaryEleOp>
OpEvaluateSDFImpl<DIM, GAUSS, BoundaryEleOp>::OpEvaluateSDFImpl(
    boost::shared_ptr<CommonData> common_data_ptr)
    : BoundaryEleOp(NOSPACE, BoundaryEleOp::OPSPACE),
      commonDataPtr(common_data_ptr) {}

template <int DIM, typename BoundaryEleOp>
MoFEMErrorCode
OpEvaluateSDFImpl<DIM, GAUSS, BoundaryEleOp>::doWork(int side, EntityType type,
                                                     EntData &data) {
  MoFEMFunctionBegin;

  const auto nb_gauss_pts = BoundaryEleOp::getGaussPts().size2();
  auto &sdf_vec = commonDataPtr->sdfVals;
  auto &grad_mat = commonDataPtr->gradsSdf;
  auto &hess_mat = commonDataPtr->hessSdf;
  auto &constraint_vec = commonDataPtr->constraintVals;
  auto &contactTraction_mat = commonDataPtr->contactTraction;

  sdf_vec.resize(nb_gauss_pts, false);
  grad_mat.resize(DIM, nb_gauss_pts, false);
  hess_mat.resize((DIM * (DIM + 1)) / 2, nb_gauss_pts, false);
  constraint_vec.resize(nb_gauss_pts, false);

  auto t_traction = getFTensor1FromMat<DIM>(contactTraction_mat);

  auto t_sdf = getFTensor0FromVec(sdf_vec);
  auto t_grad_sdf = getFTensor1FromMat<DIM>(grad_mat);
  auto t_hess_sdf = getFTensor2SymmetricFromMat<DIM>(hess_mat);
  auto t_constraint = getFTensor0FromVec(constraint_vec);

  auto t_disp = getFTensor1FromMat<DIM>(commonDataPtr->contactDisp);
  auto t_coords = BoundaryEleOp::getFTensor1CoordsAtGaussPts();
  auto t_normal_at_pts = BoundaryEleOp::getFTensor1NormalsAtGaussPts();

  FTensor::Index<'i', DIM> i;
  FTensor::Index<'j', DIM> j;

  auto ts_time = BoundaryEleOp::getTStime();
  auto ts_time_step = BoundaryEleOp::getTStimeStep();

  auto m_spatial_coords = get_spatial_coords(
      BoundaryEleOp::getFTensor1CoordsAtGaussPts(),
      getFTensor1FromMat<DIM>(commonDataPtr->contactDisp), nb_gauss_pts);
  auto m_normals_at_pts = get_normalize_normals(
      BoundaryEleOp::getFTensor1NormalsAtGaussPts(), nb_gauss_pts);

  // placeholder to pass boundary block id to python
  int block_id = 0;

  auto v_sdf =
      surfaceDistanceFunction(ts_time_step, ts_time, nb_gauss_pts,
                              m_spatial_coords, m_normals_at_pts, block_id);

  auto m_grad_sdf =
      gradSurfaceDistanceFunction(ts_time_step, ts_time, nb_gauss_pts,
                                  m_spatial_coords, m_normals_at_pts, block_id);

  auto m_hess_sdf =
      hessSurfaceDistanceFunction(ts_time_step, ts_time, nb_gauss_pts,
                                  m_spatial_coords, m_normals_at_pts, block_id);

  auto t_sdf_v = getFTensor0FromVec(v_sdf);
  auto t_grad_sdf_v = getFTensor1FromMat<3>(m_grad_sdf);
  auto t_hess_sdf_v = getFTensor2SymmetricFromMat<3>(m_hess_sdf);

  auto next = [&]() {
    ++t_sdf;
    ++t_sdf_v;
    ++t_grad_sdf;
    ++t_grad_sdf_v;
    ++t_hess_sdf;
    ++t_hess_sdf_v;
    ++t_disp;
    ++t_traction;
    ++t_constraint;
  };

  for (auto gg = 0; gg != nb_gauss_pts; ++gg) {

    auto tn = -t_traction(i) * t_grad_sdf_v(i);
    auto c = constrain(t_sdf_v, tn);

    t_sdf = t_sdf_v;
    t_grad_sdf(i) = t_grad_sdf_v(i);
    t_hess_sdf(i, j) = t_hess_sdf_v(i, j);
    t_constraint = c;

    next();
  }

  MoFEMFunctionReturn(0);
}

template <int DIM, typename AssemblyBoundaryEleOp>
OpConstrainBoundaryRhsImpl<DIM, GAUSS, AssemblyBoundaryEleOp>::
    OpConstrainBoundaryRhsImpl(const std::string field_name,
                               boost::shared_ptr<CommonData> common_data_ptr,
                               bool is_axisymmetric)
    : AssemblyBoundaryEleOp(field_name, field_name,
                            AssemblyBoundaryEleOp::OPROW),
      commonDataPtr(common_data_ptr), isAxisymmetric(is_axisymmetric) {}

template <int DIM, typename AssemblyBoundaryEleOp>
MoFEMErrorCode
OpConstrainBoundaryRhsImpl<DIM, GAUSS, AssemblyBoundaryEleOp>::iNtegrate(
    EntitiesFieldData::EntData &data) {
  MoFEMFunctionBegin;

  FTensor::Index<'i', DIM> i;
  FTensor::Index<'j', DIM> j;
  FTensor::Index<'k', DIM> k;
  FTensor::Index<'l', DIM> l;

  const size_t nb_gauss_pts = AssemblyBoundaryEleOp::getGaussPts().size2();

  auto &nf = AssemblyBoundaryEleOp::locF;

  auto t_normal_at_pts = AssemblyBoundaryEleOp::getFTensor1NormalsAtGaussPts();

  auto t_w = AssemblyBoundaryEleOp::getFTensor0IntegrationWeight();
  auto t_disp = getFTensor1FromMat<DIM>(commonDataPtr->contactDisp);
  auto t_traction = getFTensor1FromMat<DIM>(commonDataPtr->contactTraction);
  auto t_coords = AssemblyBoundaryEleOp::getFTensor1CoordsAtGaussPts();

  size_t nb_base_functions = data.getN().size2() / 3;
  auto t_base = data.getFTensor1N<3>();

  auto m_spatial_coords = get_spatial_coords(
      BoundaryEleOp::getFTensor1CoordsAtGaussPts(),
      getFTensor1FromMat<DIM>(commonDataPtr->contactDisp), nb_gauss_pts);
  auto m_normals_at_pts = get_normalize_normals(
      BoundaryEleOp::getFTensor1NormalsAtGaussPts(), nb_gauss_pts);

  auto t_normal = getFTensor1FromMat<3>(m_normals_at_pts);

  auto ts_time = AssemblyBoundaryEleOp::getTStime();
  auto ts_time_step = AssemblyBoundaryEleOp::getTStimeStep();

  // placeholder to pass boundary block id to python
  int block_id = 0;

  auto v_sdf =
      surfaceDistanceFunction(ts_time_step, ts_time, nb_gauss_pts,
                              m_spatial_coords, m_normals_at_pts, block_id);

  auto m_grad_sdf =
      gradSurfaceDistanceFunction(ts_time_step, ts_time, nb_gauss_pts,
                                  m_spatial_coords, m_normals_at_pts, block_id);

  auto t_sdf = getFTensor0FromVec(v_sdf);
  auto t_grad_sdf = getFTensor1FromMat<3>(m_grad_sdf);

  for (size_t gg = 0; gg != nb_gauss_pts; ++gg) {

    auto t_nf = getFTensor1FromPtr<DIM>(&nf[0]);

    double jacobian = 1.;
    if (isAxisymmetric) {
      jacobian = 2. * M_PI * t_coords(0);
    }
    const double alpha = t_w * jacobian * AssemblyBoundaryEleOp::getMeasure();

    auto tn = -t_traction(i) * t_grad_sdf(i);
    auto c = constrain(t_sdf, tn);

    FTensor::Tensor2<double, DIM, DIM> t_cP;
    t_cP(i, j) = (c * t_grad_sdf(i)) * t_grad_sdf(j);
    FTensor::Tensor2<double, DIM, DIM> t_cQ;
    t_cQ(i, j) = kronecker_delta(i, j) - t_cP(i, j);

    FTensor::Tensor1<double, DIM> t_rhs;
    t_rhs(i) =

        t_cQ(i, j) * (t_disp(j) - cn_contact * t_traction(j))

        +

        t_cP(i, j) * t_disp(j) +
        c * (t_sdf * t_grad_sdf(i)); // add gap0 displacements

    size_t bb = 0;
    for (; bb != AssemblyBoundaryEleOp::nbRows / DIM; ++bb) {
      const double beta = alpha * (t_base(i) * t_normal(i));
      t_nf(i) -= beta * t_rhs(i);

      ++t_nf;
      ++t_base;
    }
    for (; bb < nb_base_functions; ++bb)
      ++t_base;

    ++t_disp;
    ++t_traction;
    ++t_coords;
    ++t_w;
    ++t_normal;
    ++t_sdf;
    ++t_grad_sdf;
  }

  MoFEMFunctionReturn(0);
}

template <int DIM, typename AssemblyBoundaryEleOp>
OpConstrainBoundaryLhs_dUImpl<DIM, GAUSS, AssemblyBoundaryEleOp>::
    OpConstrainBoundaryLhs_dUImpl(const std::string row_field_name,
                                  const std::string col_field_name,
                                  boost::shared_ptr<CommonData> common_data_ptr,
                                  bool is_axisymmetric)
    : AssemblyBoundaryEleOp(row_field_name, col_field_name,
                            AssemblyBoundaryEleOp::OPROWCOL),
      commonDataPtr(common_data_ptr), isAxisymmetric(is_axisymmetric) {
  AssemblyBoundaryEleOp::sYmm = false;
}

template <int DIM, typename AssemblyBoundaryEleOp>
MoFEMErrorCode
OpConstrainBoundaryLhs_dUImpl<DIM, GAUSS, AssemblyBoundaryEleOp>::iNtegrate(
    EntitiesFieldData::EntData &row_data,
    EntitiesFieldData::EntData &col_data) {
  MoFEMFunctionBegin;

  FTensor::Index<'i', DIM> i;
  FTensor::Index<'j', DIM> j;
  FTensor::Index<'k', DIM> k;

  const size_t nb_gauss_pts = AssemblyBoundaryEleOp::getGaussPts().size2();
  auto &locMat = AssemblyBoundaryEleOp::locMat;

  auto t_normal_at_pts = AssemblyBoundaryEleOp::getFTensor1NormalsAtGaussPts();
  auto t_traction = getFTensor1FromMat<DIM>(commonDataPtr->contactTraction);
  auto t_coords = AssemblyBoundaryEleOp::getFTensor1CoordsAtGaussPts();

  auto t_w = AssemblyBoundaryEleOp::getFTensor0IntegrationWeight();
  auto t_row_base = row_data.getFTensor1N<3>();
  size_t nb_face_functions = row_data.getN().size2() / 3;

  auto m_spatial_coords = get_spatial_coords(
      BoundaryEleOp::getFTensor1CoordsAtGaussPts(),
      getFTensor1FromMat<DIM>(commonDataPtr->contactDisp), nb_gauss_pts);
  auto m_normals_at_pts = get_normalize_normals(
      BoundaryEleOp::getFTensor1NormalsAtGaussPts(), nb_gauss_pts);

  auto t_normal = getFTensor1FromMat<3>(m_normals_at_pts);

  auto ts_time = AssemblyBoundaryEleOp::getTStime();
  auto ts_time_step = AssemblyBoundaryEleOp::getTStimeStep();

  // placeholder to pass boundary block id to python
  int block_id = 0;

  auto v_sdf =
      surfaceDistanceFunction(ts_time_step, ts_time, nb_gauss_pts,
                              m_spatial_coords, m_normals_at_pts, block_id);

  auto m_grad_sdf =
      gradSurfaceDistanceFunction(ts_time_step, ts_time, nb_gauss_pts,
                                  m_spatial_coords, m_normals_at_pts, block_id);

  auto m_hess_sdf =
      hessSurfaceDistanceFunction(ts_time_step, ts_time, nb_gauss_pts,
                                  m_spatial_coords, m_normals_at_pts, block_id);

  auto t_sdf = getFTensor0FromVec(v_sdf);
  auto t_grad_sdf = getFTensor1FromMat<3>(m_grad_sdf);
  auto t_hess_sdf = getFTensor2SymmetricFromMat<3>(m_hess_sdf);

  for (size_t gg = 0; gg != nb_gauss_pts; ++gg) {

    double jacobian = 1.;
    if (isAxisymmetric) {
      jacobian = 2. * M_PI * t_coords(0);
    }
    const double alpha = t_w * jacobian * AssemblyBoundaryEleOp::getMeasure();

    auto tn = -t_traction(i) * t_grad_sdf(i);
    auto c = constrain(t_sdf, tn);

    FTensor::Tensor2<double, DIM, DIM> t_cP;
    t_cP(i, j) = (c * t_grad_sdf(i)) * t_grad_sdf(j);
    FTensor::Tensor2<double, DIM, DIM> t_cQ;
    t_cQ(i, j) = kronecker_delta(i, j) - t_cP(i, j);

    FTensor::Tensor2<double, DIM, DIM> t_res_dU;
    t_res_dU(i, j) = kronecker_delta(i, j) + t_cP(i, j);

    if (c > 0) {
      t_res_dU(i, j) +=
          (c * cn_contact) *
              (t_hess_sdf(i, j) * (t_grad_sdf(k) * t_traction(k)) +
               t_grad_sdf(i) * t_hess_sdf(k, j) * t_traction(k)) +
          c * t_sdf * t_hess_sdf(i, j);
    }

    size_t rr = 0;
    for (; rr != AssemblyBoundaryEleOp::nbRows / DIM; ++rr) {

      auto t_mat = getFTensor2FromArray<DIM, DIM, DIM>(locMat, DIM * rr);

      const double row_base = t_row_base(i) * t_normal(i);

      auto t_col_base = col_data.getFTensor0N(gg, 0);
      for (size_t cc = 0; cc != AssemblyBoundaryEleOp::nbCols / DIM; ++cc) {
        const double beta = alpha * row_base * t_col_base;

        t_mat(i, j) -= beta * t_res_dU(i, j);

        ++t_col_base;
        ++t_mat;
      }

      ++t_row_base;
    }
    for (; rr < nb_face_functions; ++rr)
      ++t_row_base;

    ++t_traction;
    ++t_coords;
    ++t_w;
    ++t_normal;
    ++t_sdf;
    ++t_grad_sdf;
    ++t_hess_sdf;
  }

  MoFEMFunctionReturn(0);
}

template <int DIM, typename AssemblyBoundaryEleOp>
OpConstrainBoundaryLhs_dTractionImpl<DIM, GAUSS, AssemblyBoundaryEleOp>::
    OpConstrainBoundaryLhs_dTractionImpl(
        const std::string row_field_name, const std::string col_field_name,
        boost::shared_ptr<CommonData> common_data_ptr, bool is_axisymmetric)
    : AssemblyBoundaryEleOp(row_field_name, col_field_name,
                            AssemblyBoundaryEleOp::OPROWCOL),
      commonDataPtr(common_data_ptr), isAxisymmetric(is_axisymmetric) {
  AssemblyBoundaryEleOp::sYmm = false;
}

template <int DIM, typename AssemblyBoundaryEleOp>
MoFEMErrorCode
OpConstrainBoundaryLhs_dTractionImpl<DIM, GAUSS, AssemblyBoundaryEleOp>::
    iNtegrate(EntitiesFieldData::EntData &row_data,
              EntitiesFieldData::EntData &col_data) {
  MoFEMFunctionBegin;

  FTensor::Index<'i', DIM> i;
  FTensor::Index<'j', DIM> j;
  FTensor::Index<'k', DIM> k;

  const size_t nb_gauss_pts = AssemblyBoundaryEleOp::getGaussPts().size2();
  auto &locMat = AssemblyBoundaryEleOp::locMat;

  auto t_normal_at_pts = AssemblyBoundaryEleOp::getFTensor1NormalsAtGaussPts();
  auto t_traction = getFTensor1FromMat<DIM>(commonDataPtr->contactTraction);
  auto t_coords = AssemblyBoundaryEleOp::getFTensor1CoordsAtGaussPts();

  auto t_w = AssemblyBoundaryEleOp::getFTensor0IntegrationWeight();
  auto t_row_base = row_data.getFTensor1N<3>();
  size_t nb_face_functions = row_data.getN().size2() / 3;

  auto m_spatial_coords = get_spatial_coords(
      BoundaryEleOp::getFTensor1CoordsAtGaussPts(),
      getFTensor1FromMat<DIM>(commonDataPtr->contactDisp), nb_gauss_pts);
  auto m_normals_at_pts = get_normalize_normals(
      BoundaryEleOp::getFTensor1NormalsAtGaussPts(), nb_gauss_pts);

  auto t_normal = getFTensor1FromMat<3>(m_normals_at_pts);

  auto ts_time = AssemblyBoundaryEleOp::getTStime();
  auto ts_time_step = AssemblyBoundaryEleOp::getTStimeStep();

  // placeholder to pass boundary block id to python
  int block_id = 0;

  auto v_sdf =
      surfaceDistanceFunction(ts_time_step, ts_time, nb_gauss_pts,
                              m_spatial_coords, m_normals_at_pts, block_id);

  auto m_grad_sdf =
      gradSurfaceDistanceFunction(ts_time_step, ts_time, nb_gauss_pts,
                                  m_spatial_coords, m_normals_at_pts, block_id);

  auto t_sdf = getFTensor0FromVec(v_sdf);
  auto t_grad_sdf = getFTensor1FromMat<3>(m_grad_sdf);

  for (size_t gg = 0; gg != nb_gauss_pts; ++gg) {

    double jacobian = 1.;
    if (isAxisymmetric) {
      jacobian = 2. * M_PI * t_coords(0);
    }
    const double alpha = t_w * jacobian * AssemblyBoundaryEleOp::getMeasure();

    auto tn = -t_traction(i) * t_grad_sdf(i);
    auto c = constrain(t_sdf, tn);

    FTensor::Tensor2<double, DIM, DIM> t_cP;
    t_cP(i, j) = (c * t_grad_sdf(i)) * t_grad_sdf(j);
    FTensor::Tensor2<double, DIM, DIM> t_cQ;
    t_cQ(i, j) = kronecker_delta(i, j) - t_cP(i, j);

    FTensor::Tensor2<double, DIM, DIM> t_res_dt;
    t_res_dt(i, j) = -cn_contact * t_cQ(i, j);

    size_t rr = 0;
    for (; rr != AssemblyBoundaryEleOp::nbRows / DIM; ++rr) {

      auto t_mat = getFTensor2FromArray<DIM, DIM, DIM>(locMat, DIM * rr);
      const double row_base = t_row_base(i) * t_normal(i);

      auto t_col_base = col_data.getFTensor1N<3>(gg, 0);
      for (size_t cc = 0; cc != AssemblyBoundaryEleOp::nbCols / DIM; ++cc) {
        const double col_base = t_col_base(i) * t_normal(i);
        const double beta = alpha * row_base * col_base;

        t_mat(i, j) -= beta * t_res_dt(i, j);

        ++t_col_base;
        ++t_mat;
      }

      ++t_row_base;
    }
    for (; rr < nb_face_functions; ++rr)
      ++t_row_base;

    ++t_traction;
    ++t_coords;
    ++t_w;
    ++t_normal;
    ++t_sdf;
    ++t_grad_sdf;
  }

  MoFEMFunctionReturn(0);
}

template <int DIM, AssemblyType A, IntegrationType I, typename DomainEleOp>
MoFEMErrorCode opFactoryDomainRhs(
    boost::ptr_deque<ForcesAndSourcesCore::UserDataOperator> &pip,
    std::string sigma, std::string u, bool is_axisymmetric = false) {
  MoFEMFunctionBegin;

  using B = typename FormsIntegrators<DomainEleOp>::template Assembly<
      A>::template LinearForm<I>;
  using OpMixDivURhs = typename B::template OpMixDivTimesU<3, DIM, DIM>;
  using OpMixDivUCylRhs =
      typename B::template OpMixDivTimesU<3, DIM, DIM, CYLINDRICAL>;

  using OpMixLambdaGradURhs = typename B::template OpMixTensorTimesGradU<DIM>;
  using OpMixUTimesDivLambdaRhs =
      typename B::template OpMixVecTimesDivLambda<SPACE_DIM>;
  using OpMixUTimesLambdaRhs =
      typename B::template OpGradTimesTensor<1, DIM, DIM>;

  auto common_data_ptr = boost::make_shared<ContactOps::CommonData>();
  auto mat_grad_ptr = boost::make_shared<MatrixDouble>();
  auto div_stress_ptr = boost::make_shared<MatrixDouble>();
  auto contact_stress_ptr = boost::make_shared<MatrixDouble>();

  auto jacobian = [is_axisymmetric](const double r, const double,
                                    const double) {
    if (is_axisymmetric)
      return 2. * M_PI * r;
    else
      return 1.;
  };

  pip.push_back(new OpCalculateVectorFieldValues<DIM>(
      u, common_data_ptr->contactDispPtr()));
  pip.push_back(
      new OpCalculateHVecTensorField<DIM, DIM>(sigma, contact_stress_ptr));

  if (!is_axisymmetric) {
    pip.push_back(
        new OpCalculateHVecTensorDivergence<DIM, DIM>(sigma, div_stress_ptr));
  } else {
    pip.push_back(new OpCalculateHVecTensorDivergence<DIM, DIM, CYLINDRICAL>(
        sigma, div_stress_ptr));
  }

  pip.push_back(new OpCalculateVectorFieldGradient<DIM, DIM>(u, mat_grad_ptr));

  if (!is_axisymmetric) {
    pip.push_back(
        new OpMixDivURhs(sigma, common_data_ptr->contactDispPtr(), jacobian));
  } else {
    pip.push_back(new OpMixDivUCylRhs(sigma, common_data_ptr->contactDispPtr(),
                                      jacobian));
  }

  pip.push_back(new OpMixLambdaGradURhs(sigma, mat_grad_ptr, jacobian));
  pip.push_back(new OpMixUTimesDivLambdaRhs(u, div_stress_ptr, jacobian));
  pip.push_back(new OpMixUTimesLambdaRhs(u, contact_stress_ptr, jacobian));

  MoFEMFunctionReturn(0);
}

template <typename OpMixLhs> struct OpMixLhsSide : public OpMixLhs {
  using OpMixLhs::OpMixLhs;
  MoFEMErrorCode doWork(int row_side, int col_side, EntityType row_type,
                        EntityType col_type,
                        EntitiesFieldData::EntData &row_data,
                        EntitiesFieldData::EntData &col_data) {
    MoFEMFunctionBegin;
    auto side_fe_entity = OpMixLhs::getSidePtrFE()->getFEEntityHandle();
    auto side_fe_data = OpMixLhs::getSideEntity(row_side, row_type);
    // Only assemble side which correspond to edge entity on boundary
    if (side_fe_entity == side_fe_data) {
      CHKERR OpMixLhs::doWork(row_side, col_side, row_type, col_type, row_data,
                              col_data);
    }
    MoFEMFunctionReturn(0);
  }
};

template <int DIM, AssemblyType A, IntegrationType I, typename DomainEle>
MoFEMErrorCode opFactoryBoundaryToDomainLhs(
    MoFEM::Interface &m_field,
    boost::ptr_deque<ForcesAndSourcesCore::UserDataOperator> &pip,
    std::string fe_domain_name, std::string sigma, std::string u,
    std::string geom, ForcesAndSourcesCore::RuleHookFun rule,
    bool is_axisymmetric = false) {
  MoFEMFunctionBegin;

  using DomainEleOp = typename DomainEle::UserDataOperator;

  auto op_loop_side = new OpLoopSide<DomainEle>(
      m_field, fe_domain_name, DIM, Sev::noisy,
      boost::make_shared<ForcesAndSourcesCore::UserDataOperator::AdjCache>());
  pip.push_back(op_loop_side);

  CHKERR AddHOOps<DIM, DIM, DIM>::add(op_loop_side->getOpPtrVector(),
                                      {H1, HDIV}, geom);

  using B = typename FormsIntegrators<DomainEleOp>::template Assembly<
      A>::template BiLinearForm<I>;

  using OpMixDivULhs = typename B::template OpMixDivTimesVec<DIM>;
  using OpMixDivUCylLhs =
      typename B::template OpMixDivTimesVec<DIM, CYLINDRICAL>;
  using OpLambdaGraULhs = typename B::template OpMixTensorTimesGrad<DIM>;

  using OpMixDivULhsSide = OpMixLhsSide<OpMixDivULhs>;
  using OpMixDivUCylLhsSide = OpMixLhsSide<OpMixDivUCylLhs>;
  using OpLambdaGraULhsSide = OpMixLhsSide<OpLambdaGraULhs>;

  auto unity = []() { return 1; };
  auto jacobian = [is_axisymmetric](const double r, const double,
                                    const double) {
    if (is_axisymmetric)
      return 2. * M_PI * r;
    else
      return 1.;
  };

  if (!is_axisymmetric) {
    op_loop_side->getOpPtrVector().push_back(
        new OpMixDivULhsSide(sigma, u, unity, jacobian, true));
  } else {
    op_loop_side->getOpPtrVector().push_back(
        new OpMixDivUCylLhsSide(sigma, u, unity, jacobian, true));
  }
  op_loop_side->getOpPtrVector().push_back(
      new OpLambdaGraULhsSide(sigma, u, unity, jacobian, true));

  op_loop_side->getSideFEPtr()->getRuleHook = rule;
  MoFEMFunctionReturn(0);
}

template <int DIM, AssemblyType A, IntegrationType I, typename BoundaryEleOp>
MoFEMErrorCode opFactoryBoundaryLhs(
    boost::ptr_deque<ForcesAndSourcesCore::UserDataOperator> &pip,
    std::string sigma, std::string u, bool is_axisymmetric = false) {
  MoFEMFunctionBegin;

  using C = ContactIntegrators<BoundaryEleOp>;

  auto common_data_ptr = boost::make_shared<ContactOps::CommonData>();

  pip.push_back(new OpCalculateVectorFieldValues<DIM>(
      u, common_data_ptr->contactDispPtr()));
  pip.push_back(new OpCalculateHVecTensorTrace<DIM, BoundaryEleOp>(
      sigma, common_data_ptr->contactTractionPtr()));
  pip.push_back(
      new typename C::template Assembly<A>::template OpConstrainBoundaryLhs_dU<
          DIM, GAUSS>(sigma, u, common_data_ptr, is_axisymmetric));
  pip.push_back(new typename C::template Assembly<A>::
                    template OpConstrainBoundaryLhs_dTraction<DIM, GAUSS>(
                        sigma, sigma, common_data_ptr, is_axisymmetric));

  MoFEMFunctionReturn(0);
}

template <int DIM, AssemblyType A, IntegrationType I, typename BoundaryEleOp>
MoFEMErrorCode opFactoryBoundaryRhs(
    boost::ptr_deque<ForcesAndSourcesCore::UserDataOperator> &pip,
    std::string sigma, std::string u, bool is_axisymmetric = false) {
  MoFEMFunctionBegin;

  using C = ContactIntegrators<BoundaryEleOp>;

  auto common_data_ptr = boost::make_shared<ContactOps::CommonData>();

  pip.push_back(new OpCalculateVectorFieldValues<DIM>(
      u, common_data_ptr->contactDispPtr()));
  pip.push_back(new OpCalculateHVecTensorTrace<DIM, BoundaryEleOp>(
      sigma, common_data_ptr->contactTractionPtr()));
  pip.push_back(
      new typename C::template Assembly<A>::template OpConstrainBoundaryRhs<
          DIM, GAUSS>(sigma, common_data_ptr, is_axisymmetric));

  MoFEMFunctionReturn(0);
}

template <int DIM, IntegrationType I, typename BoundaryEleOp>
MoFEMErrorCode opFactoryCalculateTraction(
    boost::ptr_deque<ForcesAndSourcesCore::UserDataOperator> &pip,
    std::string sigma, bool is_axisymmetric = false) {
  MoFEMFunctionBegin;

  using C = ContactIntegrators<BoundaryEleOp>;

  auto common_data_ptr = boost::make_shared<ContactOps::CommonData>();
  pip.push_back(new OpCalculateHVecTensorTrace<DIM, BoundaryEleOp>(
      sigma, common_data_ptr->contactTractionPtr()));
  pip.push_back(new typename C::template OpAssembleTotalContactTraction<DIM, I>(
      common_data_ptr, 1. / scale, is_axisymmetric));

  MoFEMFunctionReturn(0);
}

template <int DIM, IntegrationType I, typename BoundaryEleOp>
MoFEMErrorCode opFactoryCalculateArea(
    boost::ptr_deque<ForcesAndSourcesCore::UserDataOperator> &pip,
    OpLoopSide<SideEle> *op_loop_side, std::string sigma, std::string u,
    bool is_axisymmetric = false,
    boost::shared_ptr<Range> contact_range_ptr = nullptr) {
  MoFEMFunctionBegin;
  using C = ContactIntegrators<BoundaryEleOp>;

  auto common_data_ptr = boost::make_shared<ContactOps::CommonData>();

  op_loop_side->getOpPtrVector().push_back(
      new OpCalculateVectorFieldGradient<SPACE_DIM, SPACE_DIM>(
          "U", common_data_ptr->contactDispGradPtr()));

  if (contact_range_ptr) {
    pip.push_back(new OpCalculateVectorFieldValues<DIM>(
        u, common_data_ptr->contactDispPtr()));
    pip.push_back(new OpCalculateHVecTensorTrace<DIM, BoundaryEleOp>(
        sigma, common_data_ptr->contactTractionPtr()));
    pip.push_back(op_loop_side);
    pip.push_back(new typename C::template OpAssembleTotalContactArea<DIM, I>(
        common_data_ptr, is_axisymmetric, contact_range_ptr));
  }
  MoFEMFunctionReturn(0);
}

}; // namespace ContactOps

#endif // __CONTACTOPS_HPP__
