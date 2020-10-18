/**
 * @file matrix_function.cpp
 * @example matrix_function.cpp
 * @brief Test and example for matrix function
 *
 * For reference see \cite miehe2001algorithms
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

      FTensor::Tensor2<double, 3, 3> t_S{

          1., 0., 0.,

          0., 1., 0.,

          0., 0., 1.};

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

      auto get_norm_t4 = [](auto t) {
        double r = 0;
        for (int ii = 0; ii != 3; ++ii)
          for (int jj = 0; jj != 3; ++jj)
            for (int kk = 0; kk != 3; ++kk)
              for (int ll = 0; ll != 3; ++ll)
                r += t(ii, jj, kk, ll) * t(ii, jj, kk, ll);
        return r;
      };

      MOFEM_LOG("ATOM_TEST", Sev::verbose) << "Diff A";
      print_mat(t_A);

      // Testing against values in mathematica for 0,2 directive
      {
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
            EigenProjection<double, double>::getDD4M<3, 0, 0, 2, 0, 2>(t_L,
                                                                       t_N);
        MOFEM_LOG("ATOM_TEST", Sev::verbose) << "Diff dd4m 0";
        print_mat(t_dd4m_0);

        auto t_dd4m_1 =
            EigenProjection<double, double>::getDD4M<3, 1, 0, 2, 0, 2>(t_L,
                                                                       t_N);
        MOFEM_LOG("ATOM_TEST", Sev::verbose) << "Diff dd4m 1";
        print_mat(t_dd4m_1);

        auto t_dd4m_2 =
            EigenProjection<double, double>::getDD4M<3, 2, 0, 2, 0, 2>(t_L,
                                                                       t_N);
        MOFEM_LOG("ATOM_TEST", Sev::verbose) << "Diff dd4m 2";
        print_mat(t_dd4m_2);

        auto f = [](double v) { return exp(v); };
        auto d_f = [](double v) { return exp(v); };
        auto dd_f = [](double v) { return exp(v); };

        auto t_b = EigenProjection<double, double>::getMat<3>(t_L, t_N, f);
        MOFEM_LOG("ATOM_TEST", Sev::verbose) << "Reconstruct mat";
        print_mat(t_b);

        auto t_d =
            EigenProjection<double, double>::getDiffMat<3>(t_L, t_N, f, d_f);

        MOFEM_LOG("ATOM_TEST", Sev::verbose) << "Diff";
        print_ddg_direction(t_d, 0, 2);

        auto t_dd =
            EigenProjection<double, double>::getDiffDiffMat<decltype(t_S), 3>(
                t_L, t_N, f, d_f, dd_f, t_S);

        MOFEM_LOG("ATOM_TEST", Sev::verbose) << "Diff Diff";
        print_ddg_direction(t_dd, 0, 2);

        auto norm2_t_b = t_b(i, j) * t_b(i, j);
        MOFEM_LOG("ATOM_TEST", Sev::inform) << "norm2_t_b " << norm2_t_b;

        auto norm2_t_d = get_norm_t4(t_d);
        MOFEM_LOG("ATOM_TEST", Sev::inform) << "norm2_t_d " << norm2_t_d;

        auto norm2_t_dd = get_norm_t4(t_dd);
        MOFEM_LOG("ATOM_TEST", Sev::inform) << "norm2_t_dd " << norm2_t_dd;

        constexpr double regression_t_b = 572.543;
        constexpr double regression_t_d = 859.939;
        constexpr double regression_t_dd = 824.683;

        constexpr double eps = 1e-2;
        if (std::abs(norm2_t_b - regression_t_b) > eps)
          SETERRQ(PETSC_COMM_WORLD, MOFEM_ATOM_TEST_INVALID, "Wrong t_b");
        if (std::abs(norm2_t_d - regression_t_d) > eps)
          SETERRQ(PETSC_COMM_WORLD, MOFEM_ATOM_TEST_INVALID, "Wrong t_d");
        if (std::abs(norm2_t_dd - regression_t_dd) > eps)
          SETERRQ(PETSC_COMM_WORLD, MOFEM_ATOM_TEST_INVALID, "Wrong t_dd");
      }

      // Comparing with lapack calculated eigen values. Note resulst should be
      // invarinat to the direction of eiegn vector. Eigen vector can be
      // multiplied by -1 and result should be unchanged
      {
        int info;
        double wkopt;
        double w[3];

        std::array<double, 9> a{1.,   0.1, -0.5,

                                0.1,  2.,  0.,

                                -0.5, 0.,  3.};

        /* Query and allocate the optimal workspace */
        int lwork = -1;
        info = lapack_dsyev('V', 'U', 3, a.data(), 3, w, &wkopt, lwork);
        if (info > 0)
          SETERRQ(PETSC_COMM_SELF, MOFEM_ATOM_TEST_INVALID,
                  "The algorithm failed to compute eigenvalues.");
        lwork = (int)wkopt;
        std::vector<double> work(lwork);
        /* Solve eigenproblem */
        info = lapack_dsyev('V', 'U', 3, a.data(), 3, w, &*work.begin(), lwork);
        if (info > 0)
          SETERRQ(PETSC_COMM_SELF, MOFEM_ATOM_TEST_INVALID,
                  "The algorithm failed to compute eigenvalues.");

        FTensor::Tensor2<double, 3, 3> t_eig_vec{

            a[0 * 3 + 0], a[0 * 3 + 1], a[0 * 3 + 2],

            a[1 * 3 + 0], a[1 * 3 + 1], a[1 * 3 + 2],

            a[2 * 3 + 0], a[2 * 3 + 1], a[2 * 3 + 2]};

        FTensor::Tensor1<double, 3> t_eig_vals{w[0], w[1], w[2]};

        auto t_eig_val_diff =
            (t_eig_vals(i) - t_L(i)) * (t_eig_vals(i) - t_L(i));
        MOFEM_LOG("ATOM_TEST", Sev::inform)
            << "t_eig_val_diff " << t_eig_val_diff;

        auto f = [](double v) { return exp(v); };
        auto d_f = [](double v) { return exp(v); };
        auto dd_f = [](double v) { return exp(v); };

        auto t_b = EigenProjection<double, double>::getMat<3>(t_L, t_N, f);
        auto t_c = EigenProjection<double, double>::getMat<3>(t_eig_vals,
                                                              t_eig_vec, f);
        t_c(i, j) -= t_b(i, j);
        print_mat(t_c);

        auto norm2_t_c = t_c(i, j) * t_c(i, j);
        MOFEM_LOG("ATOM_TEST", Sev::verbose)
            << "Reconstruct mat difference with lapack eigen valyes and "
               "vectors "
            << norm2_t_c;

        constexpr double eps = 1e-8;
        if (fabs(norm2_t_c) > eps)
          SETERRQ(PETSC_COMM_SELF, MOFEM_ATOM_TEST_INVALID,
                  "Matrix not reeconstructed");
      }

      // constexpr auto t_kd = FTensor::Kronecker_Delta_symmetric<int>();
      // FTensor::Ddg<double, 3, 3> t_one;
      // t_one(i, j, k, l) = (t_kd(i, k) ^ t_kd(j, l)) / 4;
      // print_ddg(t_one);

      // Teestsimg linear function second second direvarive zero
      {

        auto f = [](double v) { return v; };
        auto d_f = [](double v) { return 1; };
        auto dd_f = [](double v) { return 0; };

        auto t_b = EigenProjection<double, double>::getMat<3>(t_L, t_N, f);

        auto t_d =
            EigenProjection<double, double>::getDiffMat<3>(t_L, t_N, f, d_f);
        // print_ddg(t_d);

        FTensor::Tensor2_symmetric<double, 3> t_c;
        t_c(i, j) = t_d(i, j, k, l) * t_b(k, l) - t_b(i, j);
        print_mat(t_c);
        auto norm2_t_c = t_c(i, j) * t_c(i, j);
        MOFEM_LOG("ATOM_TEST", Sev::verbose)
            << "Directive of symmetric matrix multiplied by matrix itself, "
               "should be equal to matrix itself "
            << norm2_t_c;
        constexpr double eps = 1e-10;
        if (norm2_t_c > eps)
          SETERRQ(PETSC_COMM_SELF, MOFEM_ATOM_TEST_INVALID,
                  "This norm should be zero");

        FTensor::Tensor2<double, 3, 3> t_S{

            1., 0., 0.,

            0., 1., 0.,

            0., 0., 1.};

        auto t_dd =
            EigenProjection<double, double>::getDiffDiffMat<decltype(t_S), 3>(
                t_L, t_N, f, d_f, dd_f, t_S);

        auto norm2_t_dd = get_norm_t4(t_dd);
        MOFEM_LOG("ATOM_TEST", Sev::inform) << "norm2_t_dd " << norm2_t_dd;
        if (norm2_t_dd > eps)
          SETERRQ(PETSC_COMM_SELF, MOFEM_ATOM_TEST_INVALID,
                  "This norm should be zero");
      }
    }
  }
  CATCH_ERRORS;

  CHKERR MoFEM::Core::Finalize();
}
