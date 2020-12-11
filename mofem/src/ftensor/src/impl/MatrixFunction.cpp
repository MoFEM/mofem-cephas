#include <type_traits>
#include <FTensor.hpp>
#include <MoFEM.hpp>
#include <Includes.hpp>
#include <MatrixFunction.hpp>
#include <MatrixFunctionTemplate.hpp>

namespace EigenMatrix {

template <typename T, int Dim>
FTensor::Tensor2_symmetric<double, Dim>
getMatImpl(Val<T, Dim> &t_val, Vec<T, Dim> &t_vec, Fun<double> f) {
  return EigenMatrixImp<T, T, 3, Dim>(t_val, t_vec).getMat(f);
}

template <typename T, int Dim>
FTensor::Ddg<double, Dim, Dim> getDiffMatImpl(Val<T, Dim> &t_val,
                                              Vec<T, Dim> &t_vec, Fun<double> f,
                                              Fun<double> d_f, const int nb) {
  switch (nb) {
  case 1:
    return EigenMatrixImp<T, T, 1, Dim>(t_val, t_vec).getDiffMat(f, d_f);
  case 2:
    return EigenMatrixImp<T, T, 2, Dim>(t_val, t_vec).getDiffMat(f, d_f);
  case 3:
    break;
  default:
    THROW_MESSAGE("third parameter should be 1,2 or 3");
  }
  return EigenMatrixImp<T, T, 3, Dim>(t_val, t_vec).getDiffMat(f, d_f);
};

template <typename T, typename S, int Dim>
FTensor::Ddg<double, Dim, Dim>
getDiffDiffMatImpl(Val<T, Dim> &t_val, Vec<T, Dim> &t_vec, Fun<double> f,
                   Fun<double> d_f, Fun<double> dd_f, S &t_S, const int nb) {
  switch (nb) {
  case 1:
    return EigenMatrixImp<T, T, 1, Dim>(t_val, t_vec)
        .getDiffDiffMat(f, d_f, dd_f, t_S);
  case 2:
    return EigenMatrixImp<T, T, 2, Dim>(t_val, t_vec)
        .getDiffDiffMat(f, d_f, dd_f, t_S);
  case 3:
    break;
  default:
    THROW_MESSAGE("third parameter should be 1,2 or 3");
  }
  return EigenMatrixImp<T, T, 3, Dim>(t_val, t_vec)
      .getDiffDiffMat(f, d_f, dd_f, t_S);
};

FTensor::Tensor2_symmetric<double, 3> getMat(Val<double, 3> &t_val,
                                             Vec<double, 3> &t_vec, Fun<double> f) {
  return getMatImpl<double, 3>(t_val, t_vec, f);
}

FTensor::Tensor2_symmetric<double, 3>
getMat(Val<FTensor::PackPtr<double *, 1>, 3> &t_val,
       Vec<FTensor::PackPtr<double *, 1>, 3> &t_vec, Fun<double> f) {
  return getMatImpl<FTensor::PackPtr<double *, 1>, 3>(t_val, t_vec, f);
}

FTensor::Ddg<double, 3, 3> getDiffMat(Val<double, 3> &t_val,
                                      Vec<double, 3> &t_vec, Fun<double> f,
                                      Fun<double> d_f, const int nb) {
  return getDiffMatImpl<double, 3>(t_val, t_vec, f, d_f, nb);
}

FTensor::Ddg<double, 3, 3>
getDiffMat(Val<FTensor::PackPtr<double *, 1>, 3> &t_val,
           Vec<FTensor::PackPtr<double *, 1>, 3> &t_vec, Fun<double> f,
           Fun<double> d_f, const int nb) {
  return getDiffMatImpl<FTensor::PackPtr<double *, 1>, 3>(t_val, t_vec, f, d_f,
                                                          nb);
}

FTensor::Ddg<double, 3, 3> getDiffDiffMat(Val<double, 3> &t_val,
                                          Vec<double, 3> &t_vec, Fun<double> f,
                                          Fun<double> d_f, Fun<double> dd_f,
                                          FTensor::Tensor2<double, 3, 3> &t_S,
                                          const int nb) {
  return getDiffDiffMatImpl<double, FTensor::Tensor2<double, 3, 3>, 3>(
      t_val, t_vec, f, d_f, dd_f, t_S, nb);
}

FTensor::Ddg<double, 3, 3>
getDiffDiffMat(Val<FTensor::PackPtr<double *, 1>, 3> &t_val,
               Vec<FTensor::PackPtr<double *, 1>, 3> &t_vec, Fun<double> f,
               Fun<double> d_f, Fun<double> dd_f,
               FTensor::Tensor2<double, 3, 3> &t_S, const int nb) {
  return getDiffDiffMatImpl<FTensor::PackPtr<double *, 1>,
                            FTensor::Tensor2<double, 3, 3>, 3>(
      t_val, t_vec, f, d_f, dd_f, t_S, nb);
}

FTensor::Ddg<double, 3, 3>
getDiffDiffMat(Val<double, 3> &t_val, Vec<double, 3> &t_vec, Fun<double> f,
               Fun<double> d_f, Fun<double> dd_f,
               FTensor::Tensor2_symmetric<double, 3> &t_S, const int nb) {
  return getDiffDiffMatImpl<double, FTensor::Tensor2_symmetric<double, 3>, 3>(
      t_val, t_vec, f, d_f, dd_f, t_S, nb);
}

FTensor::Tensor2_symmetric<double, 2> getMat(Val<double, 2> &t_val,
                                             Vec<double, 2> &t_vec,
                                             Fun<double> f) {
  return getMatImpl<double, 2>(t_val, t_vec, f);
}

FTensor::Ddg<double, 2, 2> getDiffMat(Val<double, 2> &t_val,
                                      Vec<double, 2> &t_vec, Fun<double> f,
                                      Fun<double> d_f, const int nb) {
  return getDiffMatImpl<double, 2>(t_val, t_vec, f, d_f, nb);
}

FTensor::Ddg<double, 2, 2> getDiffDiffMat(Val<double, 2> &t_val,
                                          Vec<double, 2> &t_vec, Fun<double> f,
                                          Fun<double> d_f, Fun<double> dd_f,
                                          FTensor::Tensor2<double, 2, 2> &t_S,
                                          const int nb) {
  return getDiffDiffMatImpl<double, FTensor::Tensor2<double, 2, 2>, 2>(
      t_val, t_vec, f, d_f, dd_f, t_S, nb);
}

FTensor::Ddg<double, 2, 2>
getDiffDiffMat(Val<double, 2> &t_val, Vec<double, 2> &t_vec, Fun<double> f,
               Fun<double> d_f, Fun<double> dd_f,
               FTensor::Tensor2_symmetric<double, 2> &t_S, const int nb) {
  return getDiffDiffMatImpl<double, FTensor::Tensor2_symmetric<double, 2>, 2>(
      t_val, t_vec, f, d_f, dd_f, t_S, nb);
}

} // namespace EigenMatrix