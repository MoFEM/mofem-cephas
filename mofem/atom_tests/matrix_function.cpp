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

      auto f = [](double v) { return v; };

      auto t_b = EigenProjection<double, double>::getMat<3>(t_L, t_N, f);

      cout << t_A << endl;
      cout << t_b << endl;

      auto t_d = EigenProjection<double, double>::getDiffMat<3>(t_L, t_N, f, f);

      for(int ii = 0; ii!=3; ++ii)
        for (int jj = 0; jj <= ii; ++jj)
          for (int kk = 0; kk != 3; ++kk)
            for (int ll = 0; ll <=kk; ++ll)
              cout << ii << " " << jj << " " << kk << " " << ll << " : "
                   << t_d(ii, jj, kk, ll) << endl;
    }
  }
  CATCH_ERRORS;

  CHKERR MoFEM::Core::Finalize();
}
