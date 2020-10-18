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
    FTensor::Index<'k', 3> k;
    FTensor::Index<'l', 3> l;

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
                MOFEM_LOG("ATOM_TEST", Sev::verbose)
                    << ii + 1 << " " << jj + 1 << " " << kk + 1 << " " << ll + 1
                    << " : " << t(ii, jj, kk, ll);
      };

      auto print_ddg_direction = [](auto &t, auto kk, int ll) {
        for (int ii = 0; ii != 3; ++ii)
          for (int jj = 0; jj <= ii; ++jj)
            MOFEM_LOG("ATOM_TEST", Sev::verbose)
                << ii + 1 << " " << jj + 1 << " " << kk + 1 << " " << ll + 1
                << " : " << t(ii, jj, kk, ll);
      };

      auto print_mat = [](auto &t) {
        for (int ii = 0; ii != 3; ++ii)
          for (int jj = 0; jj <= ii; ++jj)
            MOFEM_LOG("ATOM_TEST", Sev::verbose)
                << ii + 1 << " " << jj + 1 << " : " << t(ii, jj);
      };

      MOFEM_LOG("ATOM_TEST", Sev::verbose) << "Diff A";
      print_mat(t_A);

      auto t_d2m_0 =
          EigenProjection<double, double>::getD2M<3, 0, 0, 2>(t_L, t_N);
      MOFEM_LOG("ATOM_TEST", Sev::verbose) << "Diff d2m 0";
      print_mat(t_d2m_0);

      auto t_d2m_1 =
          EigenProjection<double, double>::getD2M<3, 1, 0, 2>(t_L, t_N);
      MOFEM_LOG("ATOM_TEST", Sev::verbose) << "Diff d2m 1";
      print_mat(t_d2m_1);

      auto t_d2m_2 =
          EigenProjection<double, double>::getD2M<3, 2, 0, 2>(t_L, t_N);
      MOFEM_LOG("ATOM_TEST", Sev::verbose) << "Diff d2m 2";
      print_mat(t_d2m_2);

      auto t_dd4m_0 =
          EigenProjection<double, double>::getDD4M<3,0, 0, 2, 0, 2>(t_L, t_N);
      MOFEM_LOG("ATOM_TEST", Sev::verbose) << "Diff dd4m 0";
      print_mat(t_dd4m_0);

      auto t_dd4m_1 =
          EigenProjection<double, double>::getDD4M<3, 1, 0, 2, 0, 2>(t_L, t_N);
      MOFEM_LOG("ATOM_TEST", Sev::verbose) << "Diff dd4m 1";
      print_mat(t_dd4m_1);

      auto t_dd4m_2 =
          EigenProjection<double, double>::getDD4M<3, 2, 0, 2, 0, 2>(t_L, t_N);
      MOFEM_LOG("ATOM_TEST", Sev::verbose) << "Diff dd4m 2";
      print_mat(t_dd4m_2);

      // auto f = [](double v) { return v*v; };
      // auto d_f = [](double v) { return 2*v; };
      // auto dd_f = [](double v) { return 2; };

      auto f = [](double v) { return exp(v); };
      auto d_f = [](double v) { return exp(v); };
      auto dd_f = [](double v) { return exp(v); };

      // auto f = [](double v) { return v; };
      // auto d_f = [](double v) { return 1; };
      // auto dd_f = [](double v) { return 0; };

      auto t_b = EigenProjection<double, double>::getMat<3>(t_L, t_N, f);
      MOFEM_LOG("ATOM_TEST", Sev::verbose) << "Reconstruct mat";
      print_mat(t_b);

      auto t_d =
          EigenProjection<double, double>::getDiffMat<3>(t_L, t_N, f, d_f);

      MOFEM_LOG("ATOM_TEST", Sev::verbose) << "Diff";
      print_ddg_direction(t_d, 0, 2);

      FTensor::Tensor2<double, 3, 3> t_S{

          1., 0., 0.,

          0., 1., 0.,

          0., 0., 1.};

      auto t_dd =
          EigenProjection<double, double>::getDiffDiffMat<decltype(t_S), 3>(
              t_L, t_N, f, d_f, dd_f, t_S);

      MOFEM_LOG("ATOM_TEST", Sev::verbose) << "Diff Diff";
      print_ddg_direction(t_dd, 0, 2);

      auto norm2_t_b = t_b(i, j) * t_b(i, j);
      MOFEM_LOG("ATOM_TEST", Sev::inform) << "norm2_t_b " << norm2_t_b;

      auto get_norm_t4 = [](auto t) {
        double r = 0;
        for (int ii = 0; ii != 3; ++ii)
          for (int jj = 0; jj != 3; ++jj)
            for (int kk = 0; kk != 3; ++kk)
              for (int ll = 0; ll != 3; ++ll)
                r += t(ii, jj, kk, ll) * t(ii, jj, kk, ll);
        return r;
      };

      auto norm2_t_d = get_norm_t4(t_d);
      MOFEM_LOG("ATOM_TEST", Sev::inform) << "norm2_t_d " << norm2_t_d;

      auto norm2_t_dd = get_norm_t4(t_dd);
      MOFEM_LOG("ATOM_TEST", Sev::inform) << "norm2_t_dd " << norm2_t_dd;

      constexpr double regression_t_b = 572.543;
      constexpr double regression_t_d = 859.939;
      constexpr double regression_t_dd = 859.939;

      constexpr double eps = 1e-2;
      if(std::abs(norm2_t_b - regression_t_b) > eps)
        SETERRQ(PETSC_COMM_WORLD, MOFEM_ATOM_TEST_INVALID, "Wrong t_b");
      if (std::abs(norm2_t_d - regression_t_d) > eps)
        SETERRQ(PETSC_COMM_WORLD, MOFEM_ATOM_TEST_INVALID, "Wrong t_d");
      if (std::abs(norm2_t_dd - regression_t_dd) > eps)
        SETERRQ(PETSC_COMM_WORLD, MOFEM_ATOM_TEST_INVALID, "Wrong t_dd");

    }
  }
  CATCH_ERRORS;

  CHKERR MoFEM::Core::Finalize();
}
