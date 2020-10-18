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

      auto print_ddg = [](auto &t) {
        for (int ii = 0; ii != 3; ++ii)
          for (int jj = 0; jj <= ii; ++jj)
            for (int kk = 0; kk != 3; ++kk)
              for (int ll = 0; ll <= kk; ++ll)
                MOFEM_LOG("ATOM_TEST", Sev::inform)
                    << ii + 1 << " " << jj + 1 << " " << kk + 1 << " " << ll + 1
                    << " : " << t(ii, jj, kk, ll);
      };

      auto print_mat = [](auto &t) {
        for (int ii = 0; ii != 3; ++ii)
          for (int jj = 0; jj <= ii; ++jj)
            MOFEM_LOG("ATOM_TEST", Sev::inform)
                << ii + 1 << " " << jj + 1 << " : " << t(ii, jj);
      };

      MOFEM_LOG("ATOM_TEST", Sev::inform) << "Diff A";
      print_mat(t_A);

      auto t_d2m_0 =
          EigenProjection<double, double>::getD2M<0, 3, 0, 2>(t_L, t_N);
      MOFEM_LOG("ATOM_TEST", Sev::inform) << "Diff d2m 0";
      print_mat(t_d2m_0);

      auto t_d2m_1 =
          EigenProjection<double, double>::getD2M<1, 3, 0, 2>(t_L, t_N);
      MOFEM_LOG("ATOM_TEST", Sev::inform) << "Diff d2m 1";
      print_mat(t_d2m_1);

      auto t_d2m_2 =
          EigenProjection<double, double>::getD2M<2, 3, 0, 2>(t_L, t_N);
      MOFEM_LOG("ATOM_TEST", Sev::inform) << "Diff d2m 2";
      print_mat(t_d2m_2);

      auto f = [](double v) { return v; };

      auto t_b = EigenProjection<double, double>::getMat<3>(t_L, t_N, f);
      MOFEM_LOG("ATOM_TEST", Sev::inform) << "Reconstruct mat";
      print_mat(t_b);

      auto t_d = EigenProjection<double, double>::getDiffMat<3>(t_L, t_N, f, f);

      MOFEM_LOG("ATOM_TEST", Sev::inform) << "Diff";
      print_ddg(t_d);

      FTensor::Tensor2<double, 3, 3> t_S{

          1., 0., 0.,

          0., 1., 0.,

          0., 0., 1.};

      auto t_dd =
          EigenProjection<double, double>::getDiffDiffMat<decltype(t_S), 3>(
              t_L, t_N, f, f, f, t_S);

      MOFEM_LOG("ATOM_TEST", Sev::inform) << "Diff Diff";
      print_ddg(t_dd);
    }
  }
  CATCH_ERRORS;

  CHKERR MoFEM::Core::Finalize();
}
