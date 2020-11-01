/**
 * @file matrix_function.cpp
 * @example matrix_function.cpp
 * @brief Test and example for matrix function
 *
 * For reference see \cite miehe2001algorithms
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

    auto print_ddg = [](auto &t, auto str = "") {
      constexpr double eps = 1e-6;
      for (int ii = 0; ii != 3; ++ii)
        for (int jj = 0; jj != 3; ++jj)
          for (int kk = 0; kk != 3; ++kk)
            for (int ll = 0; ll != 3; ++ll) {
              double v = t(ii, jj, kk, ll);
              double w = std::abs(v) < eps ? 0 : v;
              MOFEM_LOG("ATOM_TEST", Sev::verbose)
                  << str << std::fixed << std::setprecision(3) << std::showpos
                  << ii + 1 << " " << jj + 1 << " " << kk + 1 << " " << ll + 1
                  << " : " << w;
            }
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
        for (int jj = 0; jj != 3; ++jj)
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

    auto run_lapack = [](auto &a) {
      int info;
      double wkopt;
      double w[3];

      FTensor::Tensor2<double, 3, 3> t_a{

          a[0], a[1], a[2],

          a[3], a[4], a[5],

          a[6], a[7], a[8]};

      /* Query and allocate the optimal workspace */
      int lwork = -1;
      info = lapack_dsyev('V', 'U', 3, a.data(), 3, w, &wkopt, lwork);
      if (info > 0)
        THROW_MESSAGE("The algorithm failed to compute eigenvalues.");
      lwork = (int)wkopt;
      std::vector<double> work(lwork);
      /* Solve eigenproblem */
      info = lapack_dsyev('V', 'U', 3, a.data(), 3, w, &*work.begin(), lwork);
      if (info > 0)
        THROW_MESSAGE("The algorithm failed to compute eigenvalues.");

      FTensor::Tensor2<double, 3, 3> t_eig_vec{

          a[0 * 3 + 0], a[0 * 3 + 1], a[0 * 3 + 2],

          a[2 * 3 + 0], a[2 * 3 + 1], a[2 * 3 + 2],

          a[1 * 3 + 0], a[1 * 3 + 1], a[1 * 3 + 2]};

      FTensor::Tensor1<double, 3> t_eig_vals{w[0], w[2], w[1]};

      return std::make_tuple(t_a, t_eig_vec, t_eig_vals);
    };

    auto diff_ddg = [i, j, k, l](auto &t_1, auto &t_2) {
      constexpr double eps = 1e-4;
      for (int ii = 0; ii != 3; ++ii)
        for (int jj = 0; jj != 3; ++jj)
          for (int kk = 0; kk != 3; ++kk)
            for (int ll = 0; ll != 3; ++ll) {

              if (std::abs(t_1(ii, jj, kk, ll) - t_2(ii, jj, kk, ll)) > eps)
                MOFEM_LOG("ATOM_TEST", Sev::error)
                    << "Error " << ii << " " << jj << " " << kk << " " << ll
                    << " " << t_1(ii, jj, kk, ll) << " " << t_2(ii, jj, kk, ll);
            }

      for (int ii = 0; ii != 3; ++ii)
        for (int jj = 0; jj != 3; ++jj)
          for (int kk = 0; kk != 3; ++kk)
            for (int ll = 0; ll != 3; ++ll)
              t_1(ii, jj, kk, ll) -= t_2(ii, jj, kk, ll);
    };

    auto get_diff_matrix = [i, j, k, l, diff_ddg](auto &t_d) {
      constexpr auto t_kd = FTensor::Kronecker_Delta<double>();
      FTensor::Tensor4<double, 3, 3, 3, 3> t_d_a;
      t_d_a(i, j, k, l) = 0;

      for (int ii = 0; ii != 3; ++ii)
        for (int jj = 0; jj != 3; ++jj)
          for (int kk = 0; kk != 3; ++kk)
            for (int ll = 0; ll != 3; ++ll) {

              auto diff = [&](auto ii, auto jj, auto kk, auto ll) {
                return t_kd(ii, kk) * t_kd(jj, ll);
              };

              t_d_a(ii, jj, kk, ll) =
                  (diff(ii, jj, kk, ll) + diff(ii, jj, ll, kk)) / 2.;
            }

      diff_ddg(t_d_a, t_d);

      return t_d_a;
    };

    auto get_diff_matrix2 = [i, j, k, l, diff_ddg](auto &t_a, auto &t_d) {
      constexpr auto t_kd = FTensor::Kronecker_Delta<double>();
      FTensor::Tensor4<double, 3, 3, 3, 3> t_d_a;
      t_d_a(i, j, k, l) = 0;

      for (int ii = 0; ii != 3; ++ii)
        for (int jj = 0; jj != 3; ++jj)
          for (int kk = 0; kk != 3; ++kk)
            for (int ll = 0; ll != 3; ++ll)
              for (int zz = 0; zz != 3; ++zz) {

                auto diff = [&](auto ii, auto jj, auto kk, auto ll, int zz) {
                  return

                      t_a(ii, zz) * t_kd(zz, kk) * t_kd(jj, ll)

                      +

                      t_kd(ii, kk) * t_kd(zz, ll) * t_a(zz, jj);
                };

                t_d_a(ii, jj, kk, ll) +=
                    (diff(ii, jj, kk, ll, zz) + diff(ii, jj, ll, kk, zz)) / 2.;
              }

      diff_ddg(t_d_a, t_d);

      return t_d_a;
    };

    auto get_diff2_matrix2 = [i, j, k, l, diff_ddg](auto &t_s, auto &t_dd) {
      constexpr auto t_kd = FTensor::Kronecker_Delta<double>();
      FTensor::Tensor4<double, 3, 3, 3, 3> t_dd_a;
      t_dd_a(i, j, k, l) = 0;

      for (int ii = 0; ii != 3; ++ii)
        for (int jj = 0; jj != 3; ++jj)
          for (int kk = 0; kk != 3; ++kk)
            for (int ll = 0; ll != 3; ++ll)
              for (int mm = 0; mm != 3; ++mm)
                for (int nn = 0; nn != 3; ++nn)
                  for (int zz = 0; zz != 3; ++zz) {

                    auto diff = [&](auto ii, auto jj, auto kk, auto ll, int mm,
                                    int nn, int zz) {
                      return

                          t_s(ii, jj) * (t_kd(ii, mm) * t_kd(zz, nn) *
                                             t_kd(zz, kk) * t_kd(jj, ll)

                                         +

                                         t_kd(ii, kk) * t_kd(zz, ll) *
                                             t_kd(zz, mm) * t_kd(jj, nn));
                    };

                    t_dd_a(kk, ll, mm, nn) +=
                        (

                            diff(ii, jj, kk, ll, mm, nn, zz)

                            +

                            diff(ii, jj, ll, kk, mm, nn, zz)

                            +

                            diff(ii, jj, kk, ll, nn, mm, zz)

                            +

                            diff(ii, jj, ll, kk, nn, mm, zz)

                                ) /
                        4.;
                  }

      diff_ddg(t_dd_a, t_dd);

      return t_dd_a;
    };

    // Test matrix againsst mathematica results
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

      MOFEM_LOG("ATOM_TEST", Sev::verbose) << "Diff A";
      print_mat(t_A);

      // Testing against values in mathematica for 0,2 directive
      {
        auto f = [](double v) { return exp(v); };
        auto d_f = [](double v) { return exp(v); };
        auto dd_f = [](double v) { return exp(v); };

        auto t_b = EigenProjection<double, double, 3>::getMat(t_L, t_N, f);
        MOFEM_LOG("ATOM_TEST", Sev::verbose) << "Reconstruct mat";
        print_mat(t_b);

        auto t_d =
            EigenProjection<double, double, 3>::getDiffMat(t_L, t_N, f, d_f);

        MOFEM_LOG("ATOM_TEST", Sev::verbose) << "Diff";
        print_ddg_direction(t_d, 0, 2);

        auto t_dd = EigenProjection<double, double, 3>::getDiffDiffMat(
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
        constexpr double regression_t_dd = 859.939;

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

        std::array<double, 9> a{1.,   0.1, -0.5,

                                0.1,  2.,  0.,

                                -0.5, 0.,  3.};

        auto tuple = run_lapack(a);
        auto &t_a = std::get<0>(tuple);
        auto &t_eig_vec = std::get<1>(tuple);
        auto &t_eig_vals = std::get<2>(tuple);

        auto t_eig_val_diff =
            (t_eig_vals(i) - t_L(i)) * (t_eig_vals(i) - t_L(i));
        MOFEM_LOG("ATOM_TEST", Sev::inform)
            << "t_eig_val_diff " << t_eig_val_diff;

        auto f = [](double v) { return exp(v); };
        auto d_f = [](double v) { return exp(v); };
        auto dd_f = [](double v) { return exp(v); };

        auto t_b = EigenProjection<double, double, 3>::getMat(t_L, t_N, f);
        auto t_c = EigenProjection<double, double, 3>::getMat(t_eig_vals,
                                                              t_eig_vec, f);
        t_c(i, j) -= t_b(i, j);
        print_mat(t_c);

        auto norm2_t_c = t_c(i, j) * t_c(i, j);
        MOFEM_LOG("ATOM_TEST", Sev::inform)
            << "Reconstruct mat difference with lapack eigen valyes and "
               "vectors "
            << norm2_t_c;

        constexpr double eps = 1e-8;
        if (fabs(norm2_t_c) > eps)
          SETERRQ(PETSC_COMM_SELF, MOFEM_ATOM_TEST_INVALID,
                  "Matrix not reeconstructed");
      }

      constexpr auto t_kd = FTensor::Kronecker_Delta_symmetric<int>();
      FTensor::Ddg<double, 3, 3> t_one;
      t_one(i, j, k, l) = (t_kd(i, k) ^ t_kd(j, l)) / 4.;

      // Testsing linear function second second direvarive zero
      {
        auto f = [](double v) { return v; };
        auto d_f = [](double v) { return 1; };
        auto dd_f = [](double v) { return 0; };

        constexpr double eps = 1e-10;
        {
          auto t_b = EigenProjection<double, double, 3>::getMat(t_L, t_N, f);
          t_b(i, j) -= (t_A(i, j) || t_A(j, i)) / 2;
          auto norm2_t_b = t_b(i, j) * t_b(i, j);
          MOFEM_LOG("ATOM_TEST", Sev::inform)
              << "Result should be matrix itself " << norm2_t_b;
          if (norm2_t_b > eps)
            SETERRQ(PETSC_COMM_SELF, MOFEM_ATOM_TEST_INVALID,
                    "This norm should be zero");
        }

        {

          auto t_d =
              EigenProjection<double, double, 3>::getDiffMat(t_L, t_N, f, d_f);
          auto t_d_a = get_diff_matrix(t_d);

          MOFEM_LOG("ATOM_TEST", Sev::verbose) << "t_d_a";
          print_ddg(t_d_a, "hand ");
          MOFEM_LOG("ATOM_TEST", Sev::verbose) << "t_d";
          print_ddg(t_d, "code ");

          double nrm2_t_d_a = get_norm_t4(t_d_a);
          MOFEM_LOG("ATOM_TEST", Sev::inform)
              << "Direvarive hand calculation minus code " << nrm2_t_d_a;
          if (nrm2_t_d_a > eps)
            SETERRQ(PETSC_COMM_SELF, MOFEM_ATOM_TEST_INVALID,
                    "This norm should be zero");
        }

        {
          FTensor::Tensor2<double, 3, 3> t_S{

              1., 0., 0.,

              0., 1., 0.,

              0., 0., 1.};

          auto t_dd = EigenProjection<double, double, 3>::getDiffDiffMat(
              t_L, t_N, f, d_f, dd_f, t_S);

          auto norm2_t_dd = get_norm_t4(t_dd);
          MOFEM_LOG("ATOM_TEST", Sev::inform) << "norm2_t_dd " << norm2_t_dd;
          if (norm2_t_dd > eps)
            SETERRQ(PETSC_COMM_SELF, MOFEM_ATOM_TEST_INVALID,
                    "This norm should be zero");
        }
      }

      // Testsing quadratic function second second direvarive zero
      {
        auto f = [](double v) { return v * v; };
        auto d_f = [](double v) { return 2 * v; };
        auto dd_f = [](double v) { return 2; };

        constexpr double eps = 1e-9;

        // check if multiplication gives right value
        {
          auto t_b = EigenProjection<double, double, 3>::getMat(t_L, t_N, f);
          FTensor::Tensor2<double, 3, 3> t_a;
          t_a(i, j) = t_b(i, j) - t_A(i, k) * t_A(k, j);
          print_mat(t_a);
          auto norm2_t_a = t_a(i, j) * t_a(i, j);
          MOFEM_LOG("ATOM_TEST", Sev::inform)
              << "Result should be matrix times matrix " << norm2_t_a;
          if (norm2_t_a > eps)
            SETERRQ(PETSC_COMM_SELF, MOFEM_ATOM_TEST_INVALID,
                    "This norm should be zero");
        }

        // check first directive
        {
          auto t_d =
              EigenProjection<double, double, 3>::getDiffMat(t_L, t_N, f, d_f);
          print_ddg_direction(t_d, 0, 2);
          auto t_d_a = get_diff_matrix2(t_A, t_d);
          MOFEM_LOG("ATOM_TEST", Sev::verbose) << "t_d_a";
          print_ddg(t_d_a, "hand ");
          MOFEM_LOG("ATOM_TEST", Sev::verbose) << "t_d";
          print_ddg(t_d, "code ");
          double nrm2_t_d_a = get_norm_t4(t_d_a);
          MOFEM_LOG("ATOM_TEST", Sev::inform)
              << "Direvarive hand calculation minus code " << nrm2_t_d_a;
          if (nrm2_t_d_a > eps)
            SETERRQ(PETSC_COMM_SELF, MOFEM_ATOM_TEST_INVALID,
                    "This norm should be zero");
        }

        // check second directive
        {
          FTensor::Tensor2<double, 3, 3> t_S{

              1.,      1. / 2., 1. / 3.,

              2. / 2., 1.,      2. / 3.,

              3. / 2., 1.,      3. / 3.};

          auto t_dd = EigenProjection<double, double, 3>::getDiffDiffMat(
              t_L, t_N, f, d_f, dd_f, t_S);
          auto t_dd_a = get_diff2_matrix2(t_S, t_dd);

          MOFEM_LOG("ATOM_TEST", Sev::verbose) << "t_dd_a";
          print_ddg(t_dd_a, "hand ");
          MOFEM_LOG("ATOM_TEST", Sev::verbose) << "t_dd";
          print_ddg(t_dd, "code ");

          double nrm2_t_dd_a = get_norm_t4(t_dd_a);
          MOFEM_LOG("ATOM_TEST", Sev::inform)
              << "Direvarive hand calculation minus code " << nrm2_t_dd_a;
          if (nrm2_t_dd_a > eps)
            SETERRQ(PETSC_COMM_SELF, MOFEM_ATOM_TEST_INVALID,
                    "This norm should be zero");
        }
      }
    }

    // Testing two same eigen values
    {

      std::array<double, 9> a{0.1, 0.,  0,

                              0.,  0.1, 0.,

                              0.0, 0.,  4.};

      auto tuple = run_lapack(a);
      auto &t_a = std::get<0>(tuple);
      auto &t_eig_vecs = std::get<1>(tuple);
      auto &t_eig_vals = std::get<2>(tuple);

      auto f = [](double v) { return v; };
      auto d_f = [](double v) { return 1; };
      auto dd_f = [](double v) { return 0; };

      constexpr double eps = 1e-10;
      {
        auto t_b = EigenProjection<double, double, 3>::getMat(t_eig_vals,
                                                              t_eig_vecs, f);
        t_b(i, j) -= (t_a(i, j) || t_a(j, i)) / 2;
        auto norm2_t_b = t_b(i, j) * t_b(i, j);
        MOFEM_LOG("ATOM_TEST", Sev::inform)
            << "Result should be matrix itself " << norm2_t_b;
        if (norm2_t_b > eps)
          SETERRQ(PETSC_COMM_SELF, MOFEM_ATOM_TEST_INVALID,
                  "This norm should be zero");
      }

      {
        auto t_d = EigenProjection<double, double, 2>::getDiffMat(
            t_eig_vals, t_eig_vecs, f, d_f);
        auto t_d_a = get_diff_matrix(t_d);
        MOFEM_LOG("ATOM_TEST", Sev::verbose) << "t_d_a";
        print_ddg(t_d_a, "hand ");
        MOFEM_LOG("ATOM_TEST", Sev::verbose) << "t_d";
        print_ddg(t_d, "code ");
        double nrm2_t_d_a = get_norm_t4(t_d_a);
        MOFEM_LOG("ATOM_TEST", Sev::inform)
            << "Direvarive hand calculation minus code " << nrm2_t_d_a;
        if (nrm2_t_d_a > eps)
          SETERRQ(PETSC_COMM_SELF, MOFEM_ATOM_TEST_INVALID,
                  "This norm should be zero");
      }

      {
        auto f = [](double v) { return v * v; };
        auto d_f = [](double v) { return 2 * v; };
        auto dd_f = [](double v) { return 2; };
        auto t_d = EigenProjection<double, double, 2>::getDiffMat(
            t_eig_vals, t_eig_vecs, f, d_f);
        auto t_d_a = get_diff_matrix2(t_a, t_d);
        MOFEM_LOG("ATOM_TEST", Sev::verbose) << "t_d_a";
        print_ddg(t_d_a, "hand ");
        MOFEM_LOG("ATOM_TEST", Sev::verbose) << "t_d";
        print_ddg(t_d, "code ");
        double nrm2_t_d_a = get_norm_t4(t_d_a);
        MOFEM_LOG("ATOM_TEST", Sev::inform)
            << "Direvarive hand calculation minus code " << nrm2_t_d_a;
        if (nrm2_t_d_a > eps)
          SETERRQ(PETSC_COMM_SELF, MOFEM_ATOM_TEST_INVALID,
                  "This norm should be zero");
      }
    }

    // Testing three same eigen values
    {

      std::array<double, 9> a{4.,  0., 0,

                              0.,  4., 0.,

                              0.0, 0., 4.};

      auto f = [](double v) { return v; };
      auto d_f = [](double v) { return 1; };
      auto dd_f = [](double v) { return 0; };

      auto tuple = run_lapack(a);
      auto &t_a = std::get<0>(tuple);
      auto &t_eig_vecs = std::get<1>(tuple);
      auto &t_eig_vals = std::get<2>(tuple);

      // cerr << t_eig_vecs << endl;
      // cerr << t_eig_vals << endl;

      constexpr double eps = 1e-10;
      {
        auto t_b = EigenProjection<double, double, 3>::getMat(t_eig_vals,
                                                              t_eig_vecs, f);
        t_b(i, j) -= (t_a(i, j) || t_a(j, i)) / 2;
        auto norm2_t_b = t_b(i, j) * t_b(i, j);
        MOFEM_LOG("ATOM_TEST", Sev::inform)
            << "Result should be matrix itself " << norm2_t_b;
        if (norm2_t_b > eps)
          SETERRQ(PETSC_COMM_SELF, MOFEM_ATOM_TEST_INVALID,
                  "This norm should be zero");
      }

      {
        auto t_d = EigenProjection<double, double, 1>::getDiffMat(
            t_eig_vals, t_eig_vecs, f, d_f);
        auto t_d_a = get_diff_matrix(t_d);
        MOFEM_LOG("ATOM_TEST", Sev::verbose) << "t_d_a";
        print_ddg(t_d_a, "hand ");
        MOFEM_LOG("ATOM_TEST", Sev::verbose) << "t_d";
        print_ddg(t_d, "code ");
        double nrm2_t_d_a = get_norm_t4(t_d_a);
        MOFEM_LOG("ATOM_TEST", Sev::inform)
            << "Direvarive hand calculation minus code " << nrm2_t_d_a;
        if (nrm2_t_d_a > eps)
          SETERRQ(PETSC_COMM_SELF, MOFEM_ATOM_TEST_INVALID,
                  "This norm should be zero");
      }

      {
        auto f = [](double v) { return v * v; };
        auto d_f = [](double v) { return 2 * v; };
        auto dd_f = [](double v) { return 2; };
        auto t_d = EigenProjection<double, double, 1>::getDiffMat(
            t_eig_vals, t_eig_vecs, f, d_f);
        auto t_d_a = get_diff_matrix2(t_a, t_d);
        MOFEM_LOG("ATOM_TEST", Sev::verbose) << "t_d_a";
        print_ddg(t_d_a, "hand ");
        MOFEM_LOG("ATOM_TEST", Sev::verbose) << "t_d";
        print_ddg(t_d, "code ");
        double nrm2_t_d_a = get_norm_t4(t_d_a);
        MOFEM_LOG("ATOM_TEST", Sev::inform)
            << "Direvarive hand calculation minus code " << nrm2_t_d_a;
        if (nrm2_t_d_a > eps)
          SETERRQ(PETSC_COMM_SELF, MOFEM_ATOM_TEST_INVALID,
                  "This norm should be zero");
      }
    }

    // check second directive
    {

      std::array<double, 9> a{0.1, 0.,  0.,

                              0.,  0.1, 0.,

                              0.,  0.,  0.1};

      auto tuple = run_lapack(a);
      auto &t_a = std::get<0>(tuple);
      auto &t_eig_vecs = std::get<1>(tuple);
      auto &t_eig_vals = std::get<2>(tuple);

      t_eig_vals(0) -= 1e-4;
      t_eig_vals(2) += 1e-4;

      constexpr double eps = 1e-10;

      auto f = [](double v) { return v; };
      auto d_f = [](double v) { return 1; };
      auto dd_f = [](double v) { return 0; };

      FTensor::Tensor2<double, 3, 3> t_S{

          0., 0., 0.,

          0., 1., 0.,

          0., 0., 0.};

      // FTensor::Tensor2<double, 3, 3> t_S{

      //     1.,      1. / 2., 1. / 3.,

      //     2. / 2., 1.,      2. / 3.,

      //     3. / 2., 1.,      3. / 3.};

      auto t_dd = EigenProjection<double, double, 1>::getDiffDiffMat(
          t_eig_vals, t_eig_vecs, f, d_f, dd_f, t_S);

      MOFEM_LOG("ATOM_TEST", Sev::verbose) << "t_dd";
      print_ddg(t_dd, "test ");

      double nrm2_t_dd = get_norm_t4(t_dd);
      MOFEM_LOG("ATOM_TEST", Sev::inform)
          << "Direvarive hand calculation minus code " << nrm2_t_dd;
      if (nrm2_t_dd > eps)
        SETERRQ(PETSC_COMM_SELF, MOFEM_ATOM_TEST_INVALID,
                "This norm should be zero");
    }

    // check second directive
    {

      std::array<double, 9> a{2,  0., 0.,

                              0., 2,  0.,

                              0., 0., 2};

      auto tuple = run_lapack(a);
      auto &t_a = std::get<0>(tuple);
      auto &t_eig_vecs = std::get<1>(tuple);
      auto &t_eig_vals = std::get<2>(tuple);

      constexpr double eps = 1e-10;

      auto f = [](double v) { return v * v; };
      auto d_f = [](double v) { return 2 * v; };
      auto dd_f = [](double v) { return 2; };
      FTensor::Tensor2<double, 3, 3> t_S{

          1.,      1. / 2., 1. / 3.,

          2. / 1., 1.,      2. / 3.,

          3. / 1., 3. / 1., 1.};

      auto t_dd = EigenProjection<double, double, 1>::getDiffDiffMat(
          t_eig_vals, t_eig_vecs, f, d_f, dd_f, t_S);
      // print_ddg(t_dd, "test ");

      auto t_dd_a = get_diff2_matrix2(t_S, t_dd);

      double nrm2_t_dd_a = get_norm_t4(t_dd_a);
      MOFEM_LOG("ATOM_TEST", Sev::inform)
          << "Direvarive hand calculation minus code " << nrm2_t_dd_a;
      if (nrm2_t_dd_a > eps)
        SETERRQ(PETSC_COMM_SELF, MOFEM_ATOM_TEST_INVALID,
                "This norm should be zero");
    }

    // check second directive
    {

      std::array<double, 9> a{0.2, 0.,  0.,

                              0.,  0.2, 0.,

                              0.,  0.,  2};

      auto tuple = run_lapack(a);
      auto &t_a = std::get<0>(tuple);
      auto &t_eig_vecs = std::get<1>(tuple);
      auto &t_eig_vals = std::get<2>(tuple);

      constexpr double eps = 1e-10;

      auto f = [](double v) { return v * v; };
      auto d_f = [](double v) { return 2 * v; };
      auto dd_f = [](double v) { return 2; };
      FTensor::Tensor2<double, 3, 3> t_S{

          1.,      1. / 2., 1. / 3.,

          2. / 1., 1.,      2. / 3.,

          3. / 1., 3. / 1., 1.};

      auto t_dd = EigenProjection<double, double, 2>::getDiffDiffMat(
          t_eig_vals, t_eig_vecs, f, d_f, dd_f, t_S);
      // print_ddg(t_dd, "test ");

      auto t_dd_a = get_diff2_matrix2(t_S, t_dd);

      double nrm2_t_dd_a = get_norm_t4(t_dd_a);
      MOFEM_LOG("ATOM_TEST", Sev::inform)
          << "Direvarive hand calculation minus code " << nrm2_t_dd_a;
      if (nrm2_t_dd_a > eps)
        SETERRQ(PETSC_COMM_SELF, MOFEM_ATOM_TEST_INVALID,
                "This norm should be zero");
    }

    // check second directive exponent
    {

      std::array<double, 9> a{2,  0., 0.,

                              0., 2,  0.,

                              0., 0., 2};

      auto tuple = run_lapack(a);
      auto &t_a = std::get<0>(tuple);
      auto &t_eig_vecs = std::get<1>(tuple);
      auto &t_eig_vals = std::get<2>(tuple);

      t_eig_vals(0) -= 1e-5;
      t_eig_vals(2) += 1e-5;

      constexpr double eps = 1e-7;

      auto f = [](double v) { return exp(v); };
      auto d_f = [](double v) { return exp(v); };
      auto dd_f = [](double v) { return exp(v); };
      FTensor::Tensor2<double, 3, 3> t_S{

          1.,      1. / 2., 1. / 3.,

          2. / 1., 1.,      2. / 3.,

          3. / 1., 3. / 1., 1.};

      auto t_dd_1 = EigenProjection<double, double, 3>::getDiffDiffMat(
          t_eig_vals, t_eig_vecs, f, d_f, dd_f, t_S);
      auto t_dd_2 = EigenProjection<double, double, 1>::getDiffDiffMat(
          t_eig_vals, t_eig_vecs, f, d_f, dd_f, t_S);

      double nrm2_t_dd_t1 = get_norm_t4(t_dd_1);
      MOFEM_LOG("ATOM_TEST", Sev::inform)
          << "Direvarive nor t_dd_1 " << nrm2_t_dd_t1;

      double nrm2_t_dd_t2 = get_norm_t4(t_dd_2);
      MOFEM_LOG("ATOM_TEST", Sev::inform)
          << "Direvarive norm t_dd_2 " << nrm2_t_dd_t2;

      print_ddg(t_dd_1, "t_dd_1 ");
      print_ddg(t_dd_2, "t_dd_2 ");

      FTensor::Ddg<double,3,3> t_dd_3;
      t_dd_3(i, j, k, l) = t_dd_1(i, j, k, l) - t_dd_2(i, j, k, l);

      for (int ii = 0; ii != 3; ++ii)
        for (int jj = 0; jj != 3; ++jj)
          for (int kk = 0; kk != 3; ++kk)
            for (int ll = 0; ll != 3; ++ll) {
              constexpr double eps = 1e-4;
              if (std::abs(t_dd_3(ii, jj, kk, ll)) > eps)
                MOFEM_LOG("ATOM_TEST", Sev::error)
                    << "Error " << ii << " " << jj << " " << kk << " " << ll
                    << " " << t_dd_1(ii, jj, kk, ll) << " "
                    << t_dd_2(ii, jj, kk, ll);
            }

      double nrm2_t_dd_3 = get_norm_t4(t_dd_3);
      MOFEM_LOG("ATOM_TEST", Sev::inform)
          << "Direvarive approx. calculation minus code " << nrm2_t_dd_3;
      if (nrm2_t_dd_3 > eps)
        SETERRQ(PETSC_COMM_SELF, MOFEM_ATOM_TEST_INVALID,
                "This norm should be zero");
    }

  }
  CATCH_ERRORS;

  CHKERR MoFEM::Core::Finalize();
}
