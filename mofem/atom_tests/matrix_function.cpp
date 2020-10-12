/**
 * @file matrix_function.cpp
 * @example matrix_function.cpp
 * @brief Test and example for matrix function
 *
 */

#define FTENSOR_DEBUG
#include <FTensor.hpp>
#include <MatrixFunction.hpp>

#include <MoFEM.hpp>
using namespace MoFEM;

static char help[] = "...\n\n";

int main(int argc, char *argv[]) {

  MoFEM::Core::Initialize(&argc, &argv, (char *)0, help);

  auto core_log = logging::core::get();
  core_log->add_sink(
      LogManager::createSink(LogManager::getStrmSelf(), "ATOM_TEST"));
  LogManager::setLog("ATOM_TEST");

  try {

    FTensor::Index<'i', 3> i;
    FTensor::Index<'j', 3> j;

    // Test matrix
    {
      FTensor::Tensor2<double, 3, 3> t_A{

          1.,   0.1, -0.5,

          0.1,  2.,  0.,

          -0.5, 0.,  3.};

      FTensor::Tensor2<double, 3, 3> t_N{

          0.969837,  -0.0860972, 0.228042,

          0.0790574, 0.996073,   0.0398449,

          -0.230577, -0.0206147, 0.972836};

      FTensor::Tensor1<double, 3> t_L{0.873555, 2.00794, 3.11851};

      MatrixDouble a_sym_mat(6, 1);
      auto t_sym_a = getFTensor2SymmetricFromMat<3>(a_sym_mat);
      t_sym_a(i, j) = (t_A(i, j) || t_A(i, j)) / 2;

      MatrixDouble n_mat(9, 1);
      auto t_n = getFTensor2FromMat<3, 3>(n_mat);
      t_n(i, j) = t_N(i, j);

      MatrixDouble l_mat(3, 1);
      auto t_l = getFTensor1FromMat<3>(n_mat);
      t_l(i) = t_L(i);

      auto f = [](double) { return 1; };

      auto t_b = EigenProjection<double, 3, 1>::getMat<3>(t_l, t_n, f);
  }
}
CATCH_ERRORS;

CHKERR MoFEM::Core::Finalize();
}
