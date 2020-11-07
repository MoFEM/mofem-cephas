#include <type_traits>
#include <FTensor.hpp>
#include <MoFEM.hpp>
#include <Includes.hpp>
#include <MatrixFunction.hpp>
#include <MatrixFunctionTemplate.hpp>

namespace EigenMatrix {

template <typename T, int Dim>
FTensor::Tensor2_symmetric<T, 3>
getMatImpl(Val<T, Dim> &t_val, Vec<T, Dim> &t_vec, Fun<T> f, const int nb) {
  switch (nb) {
  case 1:
    return EigenMatrixImp<T, T, 1, Dim>(t_val, t_vec).getMat(f);
  case 2:
    return EigenMatrixImp<T, T, 2, Dim>(t_val, t_vec).getMat(f);
  case 3:
    break;
  default:
    THROW_MESSAGE("third parameter should be 1,2 or 3");
  }
  return EigenMatrixImp<T, T, 3, Dim>(t_val, t_vec).getMat(f);
}

template <typename T, int Dim>
FTensor::Ddg<T, 3, 3> getDiffMatImpl(Val<T, Dim> &t_val, Vec<T, Dim> &t_vec,
                                     Fun<T> f, Fun<T> d_f, const int nb) {
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
FTensor::Ddg<T, 3, 3> getDiffDiffMatImpl(Val<T, Dim> &t_val, Vec<T, Dim> &t_vec,
                                         Fun<T> f, Fun<T> d_f, Fun<T> dd_f,
                                         S &t_S, const int nb) {
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
                                             Vec<double, 3> &t_vec,
                                             Fun<double> f, const int nb) {
  return getMatImpl<double, 3>(t_val, t_vec, f, nb);
}

FTensor::Ddg<double, 3, 3> getDiffMat(Val<double, 3> &t_val,
                                      Vec<double, 3> &t_vec, Fun<double> f,
                                      Fun<double> d_f, const int nb) {
  return getDiffMatImpl<double, 3>(t_val, t_vec, f, d_f, nb);
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
getDiffDiffMat(Val<double, 3> &t_val, Vec<double, 3> &t_vec, Fun<double> f,
               Fun<double> d_f, Fun<double> dd_f,
               FTensor::Tensor2_symmetric<double, 3> &t_S, const int nb) {
  return getDiffDiffMatImpl<double, FTensor::Tensor2_symmetric<double, 3>, 3>(
      t_val, t_vec, f, d_f, dd_f, t_S, nb);
}

} // namespace EigenMatrix