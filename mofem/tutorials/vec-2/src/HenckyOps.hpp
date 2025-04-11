/**
 * \file HenckyOps.hpp
 * \example HenckyOps.hpp
 *
 * @copyright Copyright (c) 2023
 */

#ifndef __HENCKY_OPS_HPP__
#define __HENCKY_OPS_HPP__

namespace HenckyOps {

constexpr double eps = std::numeric_limits<float>::epsilon();

auto f = [](double v) { return 0.5 * std::log(v); };
auto d_f = [](double v) { return 0.5 / v; };
auto dd_f = [](double v) { return -0.5 / (v * v); };

struct isEq {
  static inline auto check(const double &a, const double &b) {
    return std::abs(a - b) / absMax < eps;
  }
  static double absMax;
};

double isEq::absMax = 1;

inline auto is_eq(const double &a, const double &b) {
  return isEq::check(a, b);
};

template <int DIM> inline auto get_uniq_nb(double *ptr) {
  std::array<double, DIM> tmp;
  std::copy(ptr, &ptr[DIM], tmp.begin());
  std::sort(tmp.begin(), tmp.end());
  isEq::absMax = std::max(std::abs(tmp[0]), std::abs(tmp[DIM - 1]));
  return std::distance(tmp.begin(), std::unique(tmp.begin(), tmp.end(), is_eq));
};

template <int DIM>
inline auto sort_eigen_vals(FTensor::Tensor1<double, DIM> &eig,
                            FTensor::Tensor2<double, DIM, DIM> &eigen_vec) {

  isEq::absMax =
      std::max(std::max(std::abs(eig(0)), std::abs(eig(1))), std::abs(eig(2)));

  int i = 0, j = 1, k = 2;

  if (is_eq(eig(0), eig(1))) {
    i = 0;
    j = 2;
    k = 1;
  } else if (is_eq(eig(0), eig(2))) {
    i = 0;
    j = 1;
    k = 2;
  } else if (is_eq(eig(1), eig(2))) {
    i = 1;
    j = 0;
    k = 2;
  }

  FTensor::Tensor2<double, 3, 3> eigen_vec_c{
      eigen_vec(i, 0), eigen_vec(i, 1), eigen_vec(i, 2),

      eigen_vec(j, 0), eigen_vec(j, 1), eigen_vec(j, 2),

      eigen_vec(k, 0), eigen_vec(k, 1), eigen_vec(k, 2)};

  FTensor::Tensor1<double, 3> eig_c{eig(i), eig(j), eig(k)};

  {
    FTensor::Index<'i', DIM> i;
    FTensor::Index<'j', DIM> j;
    eig(i) = eig_c(i);
    eigen_vec(i, j) = eigen_vec_c(i, j);
  }
};

struct CommonData : public boost::enable_shared_from_this<CommonData> {
  boost::shared_ptr<MatrixDouble> matGradPtr;
  boost::shared_ptr<MatrixDouble> matDPtr;
  boost::shared_ptr<MatrixDouble> matLogCPlastic;

  MatrixDouble matEigVal;
  MatrixDouble matEigVec;
  MatrixDouble matLogC;
  MatrixDouble matLogCdC;
  MatrixDouble matFirstPiolaStress;
  MatrixDouble matSecondPiolaStress;
  MatrixDouble matHenckyStress;
  MatrixDouble matTangent;

  inline auto getMatFirstPiolaStress() {
    return boost::shared_ptr<MatrixDouble>(shared_from_this(),
                                           &matFirstPiolaStress);
  }

  inline auto getMatHenckyStress() {
    return boost::shared_ptr<MatrixDouble>(shared_from_this(),
                                           &matHenckyStress);
  }

  inline auto getMatLogC() {
    return boost::shared_ptr<MatrixDouble>(shared_from_this(), &matLogC);
  }

  inline auto getMatTangent() {
    return boost::shared_ptr<MatrixDouble>(shared_from_this(), &matTangent);
  }
};

template <int DIM, IntegrationType I, typename DomainEleOp>
struct OpCalculateEigenValsImpl;

template <int DIM, IntegrationType I, typename DomainEleOp>
struct OpCalculateLogCImpl;

template <int DIM, IntegrationType I, typename DomainEleOp>
struct OpCalculateLogC_dCImpl;

template <int DIM, IntegrationType I, typename DomainEleOp, int S>
struct OpCalculateHenckyStressImpl;

template <int DIM, IntegrationType I, typename DomainEleOp, int S>
struct OpCalculateHenckyThermalStressImpl;

template <int DIM, IntegrationType I, typename DomainEleOp, int S>
struct OpCalculateHenckyThermalStressImpl;

template <int DIM, IntegrationType I, typename DomainEleOp, int S>
struct OpCalculateHenckyPlasticStressImpl;

template <int DIM, IntegrationType I, typename DomainEleOp, int S>
struct OpCalculatePiolaStressImpl;

template <int DIM, IntegrationType I, typename DomainEleOp, int S>
struct OpHenckyTangentImpl;

template <int DIM, IntegrationType I, typename DomainEleOp, int S>
struct OpCalculateHenckyThermalStressdTImpl;

template <int DIM, typename DomainEleOp>
struct OpCalculateEigenValsImpl<DIM, GAUSS, DomainEleOp> : public DomainEleOp {

  OpCalculateEigenValsImpl(const std::string field_name,
                           boost::shared_ptr<CommonData> common_data)
      : DomainEleOp(field_name, DomainEleOp::OPROW),
        commonDataPtr(common_data) {
    std::fill(&DomainEleOp::doEntities[MBEDGE],
              &DomainEleOp::doEntities[MBMAXTYPE], false);
  }

  MoFEMErrorCode doWork(int side, EntityType type, EntData &data) {
    MoFEMFunctionBegin;

    FTensor::Index<'i', DIM> i;
    FTensor::Index<'j', DIM> j;
    FTensor::Index<'k', DIM> k;

    constexpr auto t_kd = FTensor::Kronecker_Delta<int>();
    // const size_t nb_gauss_pts = matGradPtr->size2();
    const size_t nb_gauss_pts = DomainEleOp::getGaussPts().size2();

    auto t_grad = getFTensor2FromMat<DIM, DIM>(*(commonDataPtr->matGradPtr));

    commonDataPtr->matEigVal.resize(DIM, nb_gauss_pts, false);
    commonDataPtr->matEigVec.resize(DIM * DIM, nb_gauss_pts, false);
    auto t_eig_val = getFTensor1FromMat<DIM>(commonDataPtr->matEigVal);
    auto t_eig_vec = getFTensor2FromMat<DIM, DIM>(commonDataPtr->matEigVec);

    for (size_t gg = 0; gg != nb_gauss_pts; ++gg) {

      FTensor::Tensor2<double, DIM, DIM> F;
      FTensor::Tensor2_symmetric<double, DIM> C;
      FTensor::Tensor1<double, DIM> eig;
      FTensor::Tensor2<double, DIM, DIM> eigen_vec;

      F(i, j) = t_grad(i, j) + t_kd(i, j);
      C(i, j) = F(k, i) ^ F(k, j);

      for (int ii = 0; ii != DIM; ii++)
        for (int jj = 0; jj != DIM; jj++)
          eigen_vec(ii, jj) = C(ii, jj);

      CHKERR computeEigenValuesSymmetric(eigen_vec, eig);
      for (auto ii = 0; ii != DIM; ++ii)
        eig(ii) = std::max(std::numeric_limits<double>::epsilon(), eig(ii));

      // rare case when two eigen values are equal
      auto nb_uniq = get_uniq_nb<DIM>(&eig(0));
      if constexpr (DIM == 3) {
        if (nb_uniq == 2) {
          sort_eigen_vals<DIM>(eig, eigen_vec);
        }
      }

      t_eig_val(i) = eig(i);
      t_eig_vec(i, j) = eigen_vec(i, j);

      ++t_grad;
      ++t_eig_val;
      ++t_eig_vec;
    }

    MoFEMFunctionReturn(0);
  }

private:
  boost::shared_ptr<CommonData> commonDataPtr;
};

template <int DIM, typename DomainEleOp>
struct OpCalculateLogCImpl<DIM, GAUSS, DomainEleOp> : public DomainEleOp {

  OpCalculateLogCImpl(const std::string field_name,
                      boost::shared_ptr<CommonData> common_data)
      : DomainEleOp(field_name, DomainEleOp::OPROW),
        commonDataPtr(common_data) {
    std::fill(&DomainEleOp::doEntities[MBEDGE],
              &DomainEleOp::doEntities[MBMAXTYPE], false);
  }

  MoFEMErrorCode doWork(int side, EntityType type, EntData &data) {
    MoFEMFunctionBegin;

    FTensor::Index<'i', DIM> i;
    FTensor::Index<'j', DIM> j;

    // const size_t nb_gauss_pts = matGradPtr->size2();
    const size_t nb_gauss_pts = DomainEleOp::getGaussPts().size2();
    constexpr auto size_symm = (DIM * (DIM + 1)) / 2;
    commonDataPtr->matLogC.resize(size_symm, nb_gauss_pts, false);

    auto t_eig_val = getFTensor1FromMat<DIM>(commonDataPtr->matEigVal);
    auto t_eig_vec = getFTensor2FromMat<DIM, DIM>(commonDataPtr->matEigVec);

    auto t_logC = getFTensor2SymmetricFromMat<DIM>(commonDataPtr->matLogC);

    for (size_t gg = 0; gg != nb_gauss_pts; ++gg) {
      t_logC(i, j) = EigenMatrix::getMat(t_eig_val, t_eig_vec, f)(i, j);
      t_logC(i, j) = EigenMatrix::getMat(t_eig_val, t_eig_vec, f)(i, j);
      ++t_eig_val;
      ++t_eig_vec;
      ++t_logC;
    }

    MoFEMFunctionReturn(0);
  }

private:
  boost::shared_ptr<CommonData> commonDataPtr;
};

template <int DIM, typename DomainEleOp>
struct OpCalculateLogC_dCImpl<DIM, GAUSS, DomainEleOp> : public DomainEleOp {

  OpCalculateLogC_dCImpl(const std::string field_name,
                         boost::shared_ptr<CommonData> common_data)
      : DomainEleOp(field_name, DomainEleOp::OPROW),
        commonDataPtr(common_data) {
    std::fill(&DomainEleOp::doEntities[MBEDGE],
              &DomainEleOp::doEntities[MBMAXTYPE], false);
  }

  MoFEMErrorCode doWork(int side, EntityType type, EntData &data) {
    MoFEMFunctionBegin;

    FTensor::Index<'i', DIM> i;
    FTensor::Index<'j', DIM> j;
    FTensor::Index<'k', DIM> k;
    FTensor::Index<'l', DIM> l;

    const size_t nb_gauss_pts = DomainEleOp::getGaussPts().size2();
    constexpr auto size_symm = (DIM * (DIM + 1)) / 2;
    commonDataPtr->matLogCdC.resize(size_symm * size_symm, nb_gauss_pts, false);
    auto t_logC_dC = getFTensor4DdgFromMat<DIM, DIM>(commonDataPtr->matLogCdC);
    auto t_eig_val = getFTensor1FromMat<DIM>(commonDataPtr->matEigVal);
    auto t_eig_vec = getFTensor2FromMat<DIM, DIM>(commonDataPtr->matEigVec);

    for (size_t gg = 0; gg != nb_gauss_pts; ++gg) {
      // rare case when two eigen values are equal
      auto nb_uniq = get_uniq_nb<DIM>(&t_eig_val(0));
      t_logC_dC(i, j, k, l) =
          2 * EigenMatrix::getDiffMat(t_eig_val, t_eig_vec, f, d_f,
                                      nb_uniq)(i, j, k, l);
      auto nb_uniq = get_uniq_nb<DIM>(&t_eig_val(0));
      t_logC_dC(i, j, k, l) =
          2 * EigenMatrix::getDiffMat(t_eig_val, t_eig_vec, f, d_f,
                                      nb_uniq)(i, j, k, l);

      ++t_logC_dC;
      ++t_eig_val;
      ++t_eig_vec;
    }

    MoFEMFunctionReturn(0);
  }

private:
  boost::shared_ptr<CommonData> commonDataPtr;
};

template <int DIM, typename DomainEleOp, int S>
struct OpCalculateHenckyStressImpl<DIM, GAUSS, DomainEleOp, S>
    : public DomainEleOp {

  OpCalculateHenckyStressImpl(const std::string field_name,
                              boost::shared_ptr<CommonData> common_data)
      : DomainEleOp(field_name, DomainEleOp::OPROW),
        commonDataPtr(common_data) {
    std::fill(&DomainEleOp::doEntities[MBEDGE],
              &DomainEleOp::doEntities[MBMAXTYPE], false);
  }

  MoFEMErrorCode doWork(int side, EntityType type, EntData &data) {
    MoFEMFunctionBegin;

    FTensor::Index<'i', DIM> i;
    FTensor::Index<'j', DIM> j;
    FTensor::Index<'k', DIM> k;
    FTensor::Index<'l', DIM> l;

    // const size_t nb_gauss_pts = matGradPtr->size2();
    const size_t nb_gauss_pts = DomainEleOp::getGaussPts().size2();
    auto t_D = getFTensor4DdgFromMat<DIM, DIM, S>(*commonDataPtr->matDPtr);
    auto t_logC = getFTensor2SymmetricFromMat<DIM>(commonDataPtr->matLogC);
    constexpr auto size_symm = (DIM * (DIM + 1)) / 2;
    commonDataPtr->matHenckyStress.resize(size_symm, nb_gauss_pts, false);
    auto t_T = getFTensor2SymmetricFromMat<DIM>(commonDataPtr->matHenckyStress);

    for (size_t gg = 0; gg != nb_gauss_pts; ++gg) {
      t_T(i, j) = t_D(i, j, k, l) * t_logC(k, l);
      ++t_logC;
      ++t_T;
      ++t_D;
    }

    MoFEMFunctionReturn(0);
  }

private:
  boost::shared_ptr<CommonData> commonDataPtr;
};

template <int DIM, typename DomainEleOp, int S>
struct OpCalculateHenckyThermalStressImpl<DIM, GAUSS, DomainEleOp, S>
    : public DomainEleOp {

  OpCalculateHenckyThermalStressImpl(
      const std::string field_name, boost::shared_ptr<VectorDouble> temperature,
      boost::shared_ptr<CommonData> common_data,
      boost::shared_ptr<VectorDouble> coeff_expansion_ptr,
      boost::shared_ptr<double> ref_temp_ptr)
      : DomainEleOp(field_name, DomainEleOp::OPROW), tempPtr(temperature),
        commonDataPtr(common_data), coeffExpansionPtr(coeff_expansion_ptr),
        refTempPtr(ref_temp_ptr) {
    std::fill(&DomainEleOp::doEntities[MBEDGE],
              &DomainEleOp::doEntities[MBMAXTYPE], false);
  }

  MoFEMErrorCode doWork(int side, EntityType type, EntData &data) {
    MoFEMFunctionBegin;

    FTensor::Index<'i', DIM> i;
    FTensor::Index<'j', DIM> j;
    FTensor::Index<'k', DIM> k;
    FTensor::Index<'l', DIM> l;

    // const size_t nb_gauss_pts = matGradPtr->size2();
    const size_t nb_gauss_pts = DomainEleOp::getGaussPts().size2();
    auto t_D = getFTensor4DdgFromMat<DIM, DIM, S>(*commonDataPtr->matDPtr);
    auto t_logC = getFTensor2SymmetricFromMat<DIM>(commonDataPtr->matLogC);
    auto t_logC_dC = getFTensor4DdgFromMat<DIM, DIM>(commonDataPtr->matLogCdC);
    constexpr auto size_symm = (DIM * (DIM + 1)) / 2;
    commonDataPtr->matHenckyStress.resize(size_symm, nb_gauss_pts, false);
    commonDataPtr->matFirstPiolaStress.resize(DIM * DIM, nb_gauss_pts, false);
    commonDataPtr->matSecondPiolaStress.resize(size_symm, nb_gauss_pts, false);
    auto t_T = getFTensor2SymmetricFromMat<DIM>(commonDataPtr->matHenckyStress);
    auto t_P = getFTensor2FromMat<DIM, DIM>(commonDataPtr->matFirstPiolaStress);
    auto t_S =
        getFTensor2SymmetricFromMat<DIM>(commonDataPtr->matSecondPiolaStress);
    auto t_grad = getFTensor2FromMat<DIM, DIM>(*(commonDataPtr->matGradPtr));
    auto t_temp = getFTensor0FromVec(*tempPtr);

    auto t_coeff_exp = FTensor::Tensor2<double, SPACE_DIM, SPACE_DIM>();
    t_coeff_exp(i, j) = 0;
    for (auto d = 0; d != SPACE_DIM; ++d) {
      t_coeff_exp(d, d) = (*coeffExpansionPtr)[d];
    }

    for (size_t gg = 0; gg != nb_gauss_pts; ++gg) {
#ifdef HENCKY_SMALL_STRAIN
      t_P(i, j) = t_D(i, j, k, l) *
                  (t_grad(k, l) - t_coeff_exp(k, l) * (t_temp - (*refTempPtr)));
#else
      t_T(i, j) = t_D(i, j, k, l) *
                  (t_logC(k, l) - t_coeff_exp(k, l) * (t_temp - (*refTempPtr)));
      FTensor::Tensor2<double, DIM, DIM> t_F;
      t_F(i, j) = t_grad(i, j) + t_kd(i, j);
      t_S(k, l) = t_T(i, j) * t_logC_dC(i, j, k, l);
      t_P(i, l) = t_F(i, k) * t_S(k, l);
#endif
      ++t_grad;
      ++t_logC;
      ++t_logC_dC;
      ++t_P;
      ++t_T;
      ++t_S;
      ++t_D;
      ++t_temp;
    }

    MoFEMFunctionReturn(0);
  }

private:
  boost::shared_ptr<CommonData> commonDataPtr;
  boost::shared_ptr<VectorDouble> tempPtr;
  boost::shared_ptr<VectorDouble> coeffExpansionPtr;
  boost::shared_ptr<double> refTempPtr;
};

template <int DIM, typename DomainEleOp, int S>
struct OpCalculateHenckyThermalStressImpl<DIM, GAUSS, DomainEleOp, S>
    : public DomainEleOp {

  OpCalculateHenckyThermalStressImpl(
      const std::string field_name, boost::shared_ptr<VectorDouble> temperature,
      boost::shared_ptr<CommonData> common_data,
      boost::shared_ptr<VectorDouble> coeff_expansion_ptr,
      boost::shared_ptr<double> ref_temp_ptr)
      : DomainEleOp(field_name, DomainEleOp::OPROW), tempPtr(temperature),
        commonDataPtr(common_data), coeffExpansionPtr(coeff_expansion_ptr),
        refTempPtr(ref_temp_ptr) {
    std::fill(&DomainEleOp::doEntities[MBEDGE],
              &DomainEleOp::doEntities[MBMAXTYPE], false);
  }

  MoFEMErrorCode doWork(int side, EntityType type, EntData &data) {
    MoFEMFunctionBegin;

    FTensor::Index<'i', DIM> i;
    FTensor::Index<'j', DIM> j;
    FTensor::Index<'k', DIM> k;
    FTensor::Index<'l', DIM> l;

    // const size_t nb_gauss_pts = matGradPtr->size2();
    const size_t nb_gauss_pts = DomainEleOp::getGaussPts().size2();
    auto t_D = getFTensor4DdgFromMat<DIM, DIM, S>(*commonDataPtr->matDPtr);
    auto t_logC = getFTensor2SymmetricFromMat<DIM>(commonDataPtr->matLogC);
    auto t_logC_dC = getFTensor4DdgFromMat<DIM, DIM>(commonDataPtr->matLogCdC);
    constexpr auto size_symm = (DIM * (DIM + 1)) / 2;
    commonDataPtr->matHenckyStress.resize(size_symm, nb_gauss_pts, false);
    commonDataPtr->matFirstPiolaStress.resize(DIM * DIM, nb_gauss_pts, false);
    commonDataPtr->matSecondPiolaStress.resize(size_symm, nb_gauss_pts, false);
    auto t_T = getFTensor2SymmetricFromMat<DIM>(commonDataPtr->matHenckyStress);
    auto t_P = getFTensor2FromMat<DIM, DIM>(commonDataPtr->matFirstPiolaStress);
    auto t_S =
        getFTensor2SymmetricFromMat<DIM>(commonDataPtr->matSecondPiolaStress);
    auto t_grad = getFTensor2FromMat<DIM, DIM>(*(commonDataPtr->matGradPtr));
    auto t_temp = getFTensor0FromVec(*tempPtr);

    auto t_coeff_exp = FTensor::Tensor2<double, SPACE_DIM, SPACE_DIM>();
    t_coeff_exp(i, j) = 0;
    for (auto d = 0; d != SPACE_DIM; ++d) {
      t_coeff_exp(d, d) = (*coeffExpansionPtr)[d];
    }

    for (size_t gg = 0; gg != nb_gauss_pts; ++gg) {
#ifdef HENCKY_SMALL_STRAIN
      t_P(i, j) = t_D(i, j, k, l) *
                  (t_grad(k, l) - t_coeff_exp(k, l) * (t_temp - (*refTempPtr)));
#else
      t_T(i, j) = t_D(i, j, k, l) *
                  (t_logC(k, l) - t_coeff_exp(k, l) * (t_temp - (*refTempPtr)));
      FTensor::Tensor2<double, DIM, DIM> t_F;
      t_F(i, j) = t_grad(i, j) + t_kd(i, j);
      t_S(k, l) = t_T(i, j) * t_logC_dC(i, j, k, l);
      t_P(i, l) = t_F(i, k) * t_S(k, l);
#endif
      ++t_grad;
      ++t_logC;
      ++t_logC_dC;
      ++t_P;
      ++t_T;
      ++t_S;
      ++t_D;
      ++t_temp;
    }

    MoFEMFunctionReturn(0);
  }

private:
  boost::shared_ptr<CommonData> commonDataPtr;
  boost::shared_ptr<VectorDouble> tempPtr;
  boost::shared_ptr<VectorDouble> coeffExpansionPtr;
  boost::shared_ptr<double> refTempPtr;
};

template <int DIM, typename DomainEleOp, int S>
struct OpCalculateHenckyPlasticStressImpl<DIM, GAUSS, DomainEleOp, S>
    : public DomainEleOp {

  OpCalculateHenckyPlasticStressImpl(const std::string field_name,
                                     boost::shared_ptr<CommonData> common_data,
                                     boost::shared_ptr<MatrixDouble> mat_D_ptr,
                                     const double scale = 1)
      : DomainEleOp(field_name, DomainEleOp::OPROW), commonDataPtr(common_data),
        scaleStress(scale), matDPtr(mat_D_ptr) {
    std::fill(&DomainEleOp::doEntities[MBEDGE],
              &DomainEleOp::doEntities[MBMAXTYPE], false);

    matLogCPlastic = commonDataPtr->matLogCPlastic;
  }

  MoFEMErrorCode doWork(int side, EntityType type, EntData &data) {
    MoFEMFunctionBegin;

    FTensor::Index<'i', DIM> i;
    FTensor::Index<'j', DIM> j;
    FTensor::Index<'k', DIM> k;
    FTensor::Index<'l', DIM> l;

    // const size_t nb_gauss_pts = matGradPtr->size2();
    const size_t nb_gauss_pts = DomainEleOp::getGaussPts().size2();
    auto t_D = getFTensor4DdgFromMat<DIM, DIM, S>(*matDPtr);
    auto t_logC = getFTensor2SymmetricFromMat<DIM>(commonDataPtr->matLogC);
    auto t_logCPlastic = getFTensor2SymmetricFromMat<DIM>(*matLogCPlastic);
    constexpr auto size_symm = (DIM * (DIM + 1)) / 2;
    commonDataPtr->matHenckyStress.resize(size_symm, nb_gauss_pts, false);
    auto t_T = getFTensor2SymmetricFromMat<DIM>(commonDataPtr->matHenckyStress);

    for (size_t gg = 0; gg != nb_gauss_pts; ++gg) {
      t_T(i, j) = t_D(i, j, k, l) * (t_logC(k, l) - t_logCPlastic(k, l));
      t_T(i, j) /= scaleStress;
      ++t_logC;
      ++t_T;
      ++t_D;
      ++t_logCPlastic;
    }

    MoFEMFunctionReturn(0);
  }

private:
  boost::shared_ptr<CommonData> commonDataPtr;
  boost::shared_ptr<MatrixDouble> matDPtr;
  boost::shared_ptr<MatrixDouble> matLogCPlastic;
  const double scaleStress;
};

template <int DIM, typename DomainEleOp, int S>
struct OpCalculatePiolaStressImpl<DIM, GAUSS, DomainEleOp, S>
    : public DomainEleOp {

  OpCalculatePiolaStressImpl(const std::string field_name,
                             boost::shared_ptr<CommonData> common_data)
      : DomainEleOp(field_name, DomainEleOp::OPROW),
        commonDataPtr(common_data) {
    std::fill(&DomainEleOp::doEntities[MBEDGE],
              &DomainEleOp::doEntities[MBMAXTYPE], false);
  }

  MoFEMErrorCode doWork(int side, EntityType type, EntData &data) {
    MoFEMFunctionBegin;

    FTensor::Index<'i', DIM> i;
    FTensor::Index<'j', DIM> j;
    FTensor::Index<'k', DIM> k;
    FTensor::Index<'l', DIM> l;

    constexpr auto t_kd = FTensor::Kronecker_Delta<int>();

    // const size_t nb_gauss_pts = matGradPtr->size2();
    const size_t nb_gauss_pts = DomainEleOp::getGaussPts().size2();
#ifdef HENCKY_SMALL_STRAIN
    auto t_D = getFTensor4DdgFromMat<DIM, DIM, S>(*commonDataPtr->matDPtr);
#endif
    auto t_logC = getFTensor2SymmetricFromMat<DIM>(commonDataPtr->matLogC);
    auto t_logC_dC = getFTensor4DdgFromMat<DIM, DIM>(commonDataPtr->matLogCdC);
    constexpr auto size_symm = (DIM * (DIM + 1)) / 2;
    commonDataPtr->matFirstPiolaStress.resize(DIM * DIM, nb_gauss_pts, false);
    commonDataPtr->matSecondPiolaStress.resize(size_symm, nb_gauss_pts, false);
    auto t_P = getFTensor2FromMat<DIM, DIM>(commonDataPtr->matFirstPiolaStress);
    auto t_T = getFTensor2SymmetricFromMat<DIM>(commonDataPtr->matHenckyStress);
    auto t_S =
        getFTensor2SymmetricFromMat<DIM>(commonDataPtr->matSecondPiolaStress);
    auto t_grad = getFTensor2FromMat<DIM, DIM>(*(commonDataPtr->matGradPtr));

    for (size_t gg = 0; gg != nb_gauss_pts; ++gg) {

#ifdef HENCKY_SMALL_STRAIN
      t_P(i, j) = t_D(i, j, k, l) * t_grad(k, l);
#else
      FTensor::Tensor2<double, DIM, DIM> t_F;
      t_F(i, j) = t_grad(i, j) + t_kd(i, j);
      t_S(k, l) = t_T(i, j) * t_logC_dC(i, j, k, l);
      t_P(i, l) = t_F(i, k) * t_S(k, l);
#endif

      ++t_grad;
      ++t_logC;
      ++t_logC_dC;
      ++t_P;
      ++t_T;
      ++t_S;
#ifdef HENCKY_SMALL_STRAIN
      ++t_D;
#endif
    }

    MoFEMFunctionReturn(0);
  }

private:
  boost::shared_ptr<CommonData> commonDataPtr;
};

template <int DIM, typename DomainEleOp, int S>
struct OpHenckyTangentImpl<DIM, GAUSS, DomainEleOp, S> : public DomainEleOp {
  OpHenckyTangentImpl(const std::string field_name,
                      boost::shared_ptr<CommonData> common_data,
                      boost::shared_ptr<MatrixDouble> mat_D_ptr = nullptr)
      : DomainEleOp(field_name, DomainEleOp::OPROW),
        commonDataPtr(common_data) {
    std::fill(&DomainEleOp::doEntities[MBEDGE],
              &DomainEleOp::doEntities[MBMAXTYPE], false);
    if (mat_D_ptr)
      matDPtr = mat_D_ptr;
    else
      matDPtr = commonDataPtr->matDPtr;
  }

  MoFEMErrorCode doWork(int side, EntityType type, EntData &data) {
    MoFEMFunctionBegin;

    FTensor::Index<'i', DIM> i;
    FTensor::Index<'j', DIM> j;
    FTensor::Index<'k', DIM> k;
    FTensor::Index<'l', DIM> l;
    FTensor::Index<'m', DIM> m;
    FTensor::Index<'n', DIM> n;
    FTensor::Index<'o', DIM> o;
    FTensor::Index<'p', DIM> p;

    constexpr auto t_kd = FTensor::Kronecker_Delta<int>();
    // const size_t nb_gauss_pts = matGradPtr->size2();
    const size_t nb_gauss_pts = DomainEleOp::getGaussPts().size2();
    commonDataPtr->matTangent.resize(DIM * DIM * DIM * DIM, nb_gauss_pts);
    auto dP_dF =
        getFTensor4FromMat<DIM, DIM, DIM, DIM, 1>(commonDataPtr->matTangent);

    auto t_D = getFTensor4DdgFromMat<DIM, DIM, S>(*matDPtr);
    auto t_eig_val = getFTensor1FromMat<DIM>(commonDataPtr->matEigVal);
    auto t_eig_vec = getFTensor2FromMat<DIM, DIM>(commonDataPtr->matEigVec);
    auto t_T = getFTensor2SymmetricFromMat<DIM>(commonDataPtr->matHenckyStress);
    auto t_S =
        getFTensor2SymmetricFromMat<DIM>(commonDataPtr->matSecondPiolaStress);
    auto t_grad = getFTensor2FromMat<DIM, DIM>(*(commonDataPtr->matGradPtr));
    auto t_logC_dC = getFTensor4DdgFromMat<DIM, DIM>(commonDataPtr->matLogCdC);

    for (size_t gg = 0; gg != nb_gauss_pts; ++gg) {

#ifdef HENCKY_SMALL_STRAIN
      dP_dF(i, j, k, l) = t_D(i, j, k, l);
#else

      FTensor::Tensor2<double, DIM, DIM> t_F;
      t_F(i, j) = t_grad(i, j) + t_kd(i, j);

      // rare case when two eigen values are equal
      auto nb_uniq = get_uniq_nb<DIM>(&t_eig_val(0));
      FTensor::Tensor4<double, DIM, DIM, DIM, DIM> dC_dF;
      dC_dF(i, j, k, l) = (t_kd(i, l) * t_F(k, j)) + (t_kd(j, l) * t_F(k, i));

      auto TL = EigenMatrix::getDiffDiffMat(t_eig_val, t_eig_vec, f, d_f, dd_f,
                                            t_T, nb_uniq);
      TL(i, j, k, l) *= 4;
      FTensor::Ddg<double, DIM, DIM> P_D_P_plus_TL;
      P_D_P_plus_TL(i, j, k, l) =
          TL(i, j, k, l) +
          (t_logC_dC(i, j, o, p) * t_D(o, p, m, n)) * t_logC_dC(m, n, k, l);
      P_D_P_plus_TL(i, j, k, l) *= 0.5;
      dP_dF(i, j, m, n) = t_kd(i, m) * (t_kd(k, n) * t_S(k, j));
      dP_dF(i, j, m, n) +=
          t_F(i, k) * (P_D_P_plus_TL(k, j, o, p) * dC_dF(o, p, m, n));

#endif

      ++dP_dF;

      ++t_grad;
      ++t_eig_val;
      ++t_eig_vec;
      ++t_logC_dC;
      ++t_S;
      ++t_T;
      ++t_D;
    }

    MoFEMFunctionReturn(0);
  }

private:
  boost::shared_ptr<CommonData> commonDataPtr;
  boost::shared_ptr<MatrixDouble> matDPtr;
};

template <int DIM, typename AssemblyDomainEleOp, int S>
struct OpCalculateHenckyThermalStressdTImpl<DIM, GAUSS, AssemblyDomainEleOp, S>
    : public AssemblyDomainEleOp {
  OpCalculateHenckyThermalStressdTImpl(
      const std::string row_field_name, const std::string col_field_name,
      boost::shared_ptr<CommonData> elastic_common_data_ptr,
      boost::shared_ptr<VectorDouble> coeff_expansion_ptr);

  MoFEMErrorCode iNtegrate(EntitiesFieldData::EntData &row_data,
                           EntitiesFieldData::EntData &col_data);

private:
  boost::shared_ptr<CommonData> elasticCommonDataPtr;
  boost::shared_ptr<VectorDouble> coeffExpansionPtr;
};

template <int DIM, typename AssemblyDomainEleOp, int S>
OpCalculateHenckyThermalStressdTImpl<DIM, GAUSS, AssemblyDomainEleOp, S>::
    OpCalculateHenckyThermalStressdTImpl(
        const std::string row_field_name, const std::string col_field_name,
        boost::shared_ptr<CommonData> elastic_common_data_ptr,
        boost::shared_ptr<VectorDouble> coeff_expansion_ptr)
    : AssemblyDomainEleOp(row_field_name, col_field_name,
                          AssemblyDomainEleOp::OPROWCOL),
      elasticCommonDataPtr(elastic_common_data_ptr),
      coeffExpansionPtr(coeff_expansion_ptr) {
  this->sYmm = false;
}

template <int DIM, typename AssemblyDomainEleOp, int S>
MoFEMErrorCode
OpCalculateHenckyThermalStressdTImpl<DIM, GAUSS, AssemblyDomainEleOp, S>::
    iNtegrate(EntitiesFieldData::EntData &row_data,
              EntitiesFieldData::EntData &col_data) {
  MoFEMFunctionBegin;

  auto &locMat = AssemblyDomainEleOp::locMat;

  const auto nb_integration_pts = row_data.getN().size1();
  const auto nb_row_base_functions = row_data.getN().size2();
  auto t_w = this->getFTensor0IntegrationWeight();

  constexpr auto t_kd = FTensor::Kronecker_Delta<int>();
  auto t_row_diff_base = row_data.getFTensor1DiffN<DIM>();
  auto t_D = getFTensor4DdgFromMat<DIM, DIM, 0>(*elasticCommonDataPtr->matDPtr);
  auto t_grad =
      getFTensor2FromMat<DIM, DIM>(*(elasticCommonDataPtr->matGradPtr));
  auto t_logC_dC =
      getFTensor4DdgFromMat<DIM, DIM>(elasticCommonDataPtr->matLogCdC);
  FTensor::Tensor2_symmetric<double, DIM> t_eigen_strain;
  FTensor::Tensor2<double, DIM, DIM> t_F;

  FTensor::Index<'i', DIM> i;
  FTensor::Index<'j', DIM> j;
  FTensor::Index<'k', DIM> k;
  FTensor::Index<'l', DIM> l;
  FTensor::Index<'m', DIM> m;
  FTensor::Index<'n', DIM> n;
  FTensor::Index<'o', DIM> o;

  auto t_coeff_exp = FTensor::Tensor2<double, SPACE_DIM, SPACE_DIM>();
  t_coeff_exp(i, j) = 0;
  for (auto d = 0; d != SPACE_DIM; ++d) {
    t_coeff_exp(d, d) = (*coeffExpansionPtr)[d];
  }

  t_eigen_strain(i, j) = (t_D(i, j, k, l) * t_coeff_exp(k, l));

  for (auto gg = 0; gg != nb_integration_pts; ++gg) {

    t_F(i, j) = t_grad(i, j) + t_kd(i, j);

    double alpha = this->getMeasure() * t_w;
    auto rr = 0;
    for (; rr != AssemblyDomainEleOp::nbRows / DIM; ++rr) {
      auto t_mat = getFTensor1FromMat<DIM, 1>(locMat, rr * DIM);
      auto t_col_base = col_data.getFTensor0N(gg, 0);
      for (auto cc = 0; cc != AssemblyDomainEleOp::nbCols; cc++) {
#ifdef HENCKY_SMALL_STRAIN
        t_mat(i) -=
            (t_row_diff_base(j) * t_eigen_strain(i, j)) * (t_col_base * alpha);
#else
        t_mat(i) -= (t_row_diff_base(j) *
                     (t_F(i, o) * ((t_D(m, n, k, l) * t_coeff_exp(k, l)) *
                                   t_logC_dC(m, n, o, j)))) *
                    (t_col_base * alpha);
#endif

        ++t_mat;
        ++t_col_base;
      }

      ++t_row_diff_base;
    }
    for (; rr != nb_row_base_functions; ++rr)
      ++t_row_diff_base;

    ++t_w;
    ++t_grad;
    ++t_logC_dC;
    ++t_D;
  }

  MoFEMFunctionReturn(0);
}

template <typename DomainEleOp> struct HenckyIntegrators {
  template <int DIM, IntegrationType I>
  using OpCalculateEigenVals = OpCalculateEigenValsImpl<DIM, I, DomainEleOp>;

  template <int DIM, IntegrationType I>
  using OpCalculateLogC = OpCalculateLogCImpl<DIM, I, DomainEleOp>;

  template <int DIM, IntegrationType I>
  using OpCalculateLogC_dC = OpCalculateLogC_dCImpl<DIM, I, DomainEleOp>;

  template <int DIM, IntegrationType I, int S>
  using OpCalculateHenckyStress =
      OpCalculateHenckyStressImpl<DIM, I, DomainEleOp, S>;

  template <int DIM, IntegrationType I, int S>
  using OpCalculateHenckyThermalStress =
      OpCalculateHenckyThermalStressImpl<DIM, I, DomainEleOp, S>;

  template <int DIM, IntegrationType I, int S>
  using OpCalculateHenckyPlasticStress =
      OpCalculateHenckyPlasticStressImpl<DIM, I, DomainEleOp, S>;

  template <int DIM, IntegrationType I, int S>
  using OpCalculatePiolaStress =
      OpCalculatePiolaStressImpl<DIM, GAUSS, DomainEleOp, S>;

  template <int DIM, IntegrationType I, int S>
  using OpHenckyTangent = OpHenckyTangentImpl<DIM, GAUSS, DomainEleOp, S>;

  template <int DIM, IntegrationType I, typename AssemblyDomainEleOp, int S>
  using OpCalculateHenckyThermalStressdT =
      OpCalculateHenckyThermalStressdTImpl<DIM, GAUSS, AssemblyDomainEleOp, S>;
};

template <int DIM>
MoFEMErrorCode
addMatBlockOps(MoFEM::Interface &m_field,
               boost::ptr_deque<ForcesAndSourcesCore::UserDataOperator> &pip,
               std::string block_name,
               boost::shared_ptr<MatrixDouble> mat_D_Ptr, Sev sev,
               double scale = 1) {
  MoFEMFunctionBegin;

  PetscBool plane_strain_flag = PETSC_FALSE;
  CHKERR PetscOptionsGetBool(PETSC_NULL, "", "-plane_strain",
                             &plane_strain_flag, PETSC_NULL);

  struct OpMatBlocks : public DomainEleOp {
    OpMatBlocks(boost::shared_ptr<MatrixDouble> m, double bulk_modulus_K,
                double shear_modulus_G, MoFEM::Interface &m_field, Sev sev,
                std::vector<const CubitMeshSets *> meshset_vec_ptr,
                double scale, PetscBool plane_strain_flag)
        : DomainEleOp(NOSPACE, DomainEleOp::OPSPACE), matDPtr(m),
          bulkModulusKDefault(bulk_modulus_K),
          shearModulusGDefault(shear_modulus_G), scaleYoungModulus(scale),
          planeStrainFlag(plane_strain_flag) {
      CHK_THROW_MESSAGE(extractBlockData(m_field, meshset_vec_ptr, sev),
                        "Can not get data from block");
    }

    MoFEMErrorCode doWork(int side, EntityType type,
                          EntitiesFieldData::EntData &data) {
      MoFEMFunctionBegin;

      for (auto &b : blockData) {

        if (b.blockEnts.find(getFEEntityHandle()) != b.blockEnts.end()) {
          CHKERR getMatDPtr(matDPtr, b.bulkModulusK * scaleYoungModulus,
                            b.shearModulusG * scaleYoungModulus,
                            planeStrainFlag);
          MoFEMFunctionReturnHot(0);
        }
      }

      CHKERR getMatDPtr(matDPtr, bulkModulusKDefault * scaleYoungModulus,
                        shearModulusGDefault * scaleYoungModulus,
                        planeStrainFlag);
      MoFEMFunctionReturn(0);
    }

  private:
    boost::shared_ptr<MatrixDouble> matDPtr;
    const double scaleYoungModulus;
    const PetscBool planeStrainFlag;

    struct BlockData {
      double bulkModulusK;
      double shearModulusG;
      Range blockEnts;
    };

    double bulkModulusKDefault;
    double shearModulusGDefault;
    std::vector<BlockData> blockData;

    MoFEMErrorCode
    extractBlockData(MoFEM::Interface &m_field,
                     std::vector<const CubitMeshSets *> meshset_vec_ptr,
                     Sev sev) {
      MoFEMFunctionBegin;

      for (auto m : meshset_vec_ptr) {
        MOFEM_TAG_AND_LOG("WORLD", sev, "MatBlock") << *m;
        std::vector<double> block_data;
        CHKERR m->getAttributes(block_data);
        if (block_data.size() != 2) {
          SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
                  "Expected that block has two attribute");
        }
        auto get_block_ents = [&]() {
          Range ents;
          CHKERR
          m_field.get_moab().get_entities_by_handle(m->meshset, ents, true);
          return ents;
        };

        double young_modulus = block_data[0];
        double poisson_ratio = block_data[1];
        double bulk_modulus_K = young_modulus / (3 * (1 - 2 * poisson_ratio));
        double shear_modulus_G = young_modulus / (2 * (1 + poisson_ratio));

        MOFEM_TAG_AND_LOG("WORLD", sev, "MatBlock")
            << "E = " << young_modulus << " nu = " << poisson_ratio;

        blockData.push_back(
            {bulk_modulus_K, shear_modulus_G, get_block_ents()});
      }
      MOFEM_LOG_CHANNEL("WORLD");
      MoFEMFunctionReturn(0);
    }

    MoFEMErrorCode getMatDPtr(boost::shared_ptr<MatrixDouble> mat_D_ptr,
                              double bulk_modulus_K, double shear_modulus_G,
                              PetscBool is_plane_strain) {
      MoFEMFunctionBegin;
      //! [Calculate elasticity tensor]
      auto set_material_stiffness = [&]() {
        FTensor::Index<'i', DIM> i;
        FTensor::Index<'j', DIM> j;
        FTensor::Index<'k', DIM> k;
        FTensor::Index<'l', DIM> l;
        constexpr auto t_kd = FTensor::Kronecker_Delta_symmetric<int>();
        double A = (SPACE_DIM == 2 && !is_plane_strain)
                       ? 2 * shear_modulus_G /
                             (bulk_modulus_K + (4. / 3.) * shear_modulus_G)
                       : 1;
        auto t_D = getFTensor4DdgFromMat<DIM, DIM, 0>(*mat_D_ptr);
        t_D(i, j, k, l) =
            2 * shear_modulus_G * ((t_kd(i, k) ^ t_kd(j, l)) / 4.) +
            A * (bulk_modulus_K - (2. / 3.) * shear_modulus_G) * t_kd(i, j) *
                t_kd(k, l);
      };
      //! [Calculate elasticity tensor]
      constexpr auto size_symm = (DIM * (DIM + 1)) / 2;
      mat_D_ptr->resize(size_symm * size_symm, 1);
      set_material_stiffness();
      MoFEMFunctionReturn(0);
    }
  };

  double E = 1.0;
  double nu = 0.3;

  CHKERR PetscOptionsBegin(PETSC_COMM_WORLD, "", "", "none");
  CHKERR PetscOptionsScalar("-young_modulus", "Young modulus", "", E, &E,
                            PETSC_NULL);
  CHKERR PetscOptionsScalar("-poisson_ratio", "poisson ratio", "", nu, &nu,
                            PETSC_NULL);
  ierr = PetscOptionsEnd();

  double bulk_modulus_K = E / (3 * (1 - 2 * nu));
  double shear_modulus_G = E / (2 * (1 + nu));
  pip.push_back(new OpMatBlocks(
      mat_D_Ptr, bulk_modulus_K, shear_modulus_G, m_field, sev,

      // Get blockset using regular expression
      m_field.getInterface<MeshsetsManager>()->getCubitMeshsetPtr(std::regex(

          (boost::format("%s(.*)") % block_name).str()

              )),
      scale, plane_strain_flag

      ));

  MoFEMFunctionReturn(0);
}

template <int DIM, IntegrationType I, typename DomainEleOp>
auto commonDataFactory(
    MoFEM::Interface &m_field,
    boost::ptr_deque<ForcesAndSourcesCore::UserDataOperator> &pip,
    std::string field_name, std::string block_name, Sev sev, double scale = 1) {

  auto common_ptr = boost::make_shared<HenckyOps::CommonData>();
  common_ptr->matDPtr = boost::make_shared<MatrixDouble>();
  common_ptr->matGradPtr = boost::make_shared<MatrixDouble>();

  CHK_THROW_MESSAGE(addMatBlockOps<DIM>(m_field, pip, block_name,
                                        common_ptr->matDPtr, sev, scale),
                    "addMatBlockOps");

  using H = HenckyIntegrators<DomainEleOp>;

  pip.push_back(new OpCalculateVectorFieldGradient<DIM, DIM>(
      field_name, common_ptr->matGradPtr));
  pip.push_back(new typename H::template OpCalculateEigenVals<DIM, I>(
      field_name, common_ptr));
  pip.push_back(
      new typename H::template OpCalculateLogC<DIM, I>(field_name, common_ptr));
  pip.push_back(new typename H::template OpCalculateLogC_dC<DIM, I>(
      field_name, common_ptr));
  // Assumes constant D matrix per entity
  pip.push_back(new typename H::template OpCalculateHenckyStress<DIM, I, 0>(
      field_name, common_ptr));
  pip.push_back(new typename H::template OpCalculatePiolaStress<DIM, I, 0>(
      field_name, common_ptr));

  return common_ptr;
}

template <int DIM, AssemblyType A, IntegrationType I, typename DomainEleOp>
MoFEMErrorCode opFactoryDomainRhs(
    MoFEM::Interface &m_field,
    boost::ptr_deque<ForcesAndSourcesCore::UserDataOperator> &pip,
    std::string field_name, boost::shared_ptr<HenckyOps::CommonData> common_ptr,
    Sev sev) {
  MoFEMFunctionBegin;

  using B = typename FormsIntegrators<DomainEleOp>::template Assembly<
      A>::template LinearForm<I>;
  using OpInternalForcePiola =
      typename B::template OpGradTimesTensor<1, DIM, DIM>;
  pip.push_back(
      new OpInternalForcePiola("U", common_ptr->getMatFirstPiolaStress()));

  MoFEMFunctionReturn(0);
}

template <int DIM, AssemblyType A, IntegrationType I, typename DomainEleOp>
MoFEMErrorCode opFactoryDomainRhs(
    MoFEM::Interface &m_field,
    boost::ptr_deque<ForcesAndSourcesCore::UserDataOperator> &pip,
    std::string field_name, std::string block_name, Sev sev, double scale = 1) {
  MoFEMFunctionBegin;

  auto common_ptr = commonDataFactory<DIM, I, DomainEleOp>(
      m_field, pip, field_name, block_name, sev, scale);
  CHKERR opFactoryDomainRhs<DIM, A, I, DomainEleOp>(m_field, pip, field_name,
                                                    common_ptr, sev);

  MoFEMFunctionReturn(0);
}

template <int DIM, AssemblyType A, IntegrationType I, typename DomainEleOp>
MoFEMErrorCode opFactoryDomainLhs(
    MoFEM::Interface &m_field,
    boost::ptr_deque<ForcesAndSourcesCore::UserDataOperator> &pip,
    std::string field_name, boost::shared_ptr<HenckyOps::CommonData> common_ptr,
    Sev sev) {
  MoFEMFunctionBegin;

  using B = typename FormsIntegrators<DomainEleOp>::template Assembly<
      A>::template BiLinearForm<I>;
  using OpKPiola = typename B::template OpGradTensorGrad<1, DIM, DIM, 1>;

  using H = HenckyIntegrators<DomainEleOp>;
  // Assumes constant D matrix per entity
  pip.push_back(new typename H::template OpHenckyTangent<DIM, I, 0>(
      field_name, common_ptr));
  pip.push_back(
      new OpKPiola(field_name, field_name, common_ptr->getMatTangent()));

  MoFEMFunctionReturn(0);
}

template <int DIM, AssemblyType A, IntegrationType I, typename DomainEleOp>
MoFEMErrorCode opFactoryDomainLhs(
    MoFEM::Interface &m_field,
    boost::ptr_deque<ForcesAndSourcesCore::UserDataOperator> &pip,
    std::string field_name, std::string block_name, Sev sev, double scale = 1) {
  MoFEMFunctionBegin;

  auto common_ptr = commonDataFactory<DIM, I, DomainEleOp>(
      m_field, pip, field_name, block_name, sev, scale);
  CHKERR opFactoryDomainLhs<DIM, A, I, DomainEleOp>(m_field, pip, field_name,
                                                    common_ptr, sev);

  MoFEMFunctionReturn(0);
}
} // namespace HenckyOps

#endif // __HENCKY_OPS_HPP__