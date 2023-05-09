/**
 * @file matrix_function.cpp
 * @example matrix_function.cpp
 * @brief Test and example for matrix function
 *
 * For reference see \cite miehe2001algorithms
 *
 */

//#define FTENSOR_DEBUG
#include <FTensor.hpp>

#include <MoFEM.hpp>
using namespace MoFEM;

#include <MatrixFunction.hpp>

FTensor::Index<'i', 3> i;
FTensor::Index<'j', 3> j;
FTensor::Index<'k', 3> k;
FTensor::Index<'l', 3> l;

template <typename T1, typename T2, int DIM>
void diff_ddg(T1 &t_1, T2 &t_2, const FTensor::Number<DIM> &) {
  constexpr double eps = 1e-4;
  for (int ii = 0; ii != DIM; ++ii)
    for (int jj = 0; jj != DIM; ++jj)
      for (int kk = 0; kk != DIM; ++kk)
        for (int ll = 0; ll != DIM; ++ll) {

          if (std::abs(t_1(ii, jj, kk, ll) - t_2(ii, jj, kk, ll)) > eps)
            MOFEM_LOG("ATOM_TEST", Sev::error)
                << "Error " << ii << " " << jj << " " << kk << " " << ll << " "
                << t_1(ii, jj, kk, ll) << " " << t_2(ii, jj, kk, ll);
        }

  for (int ii = 0; ii != DIM; ++ii)
    for (int jj = 0; jj != DIM; ++jj)
      for (int kk = 0; kk != DIM; ++kk)
        for (int ll = 0; ll != DIM; ++ll)
          t_1(ii, jj, kk, ll) -= t_2(ii, jj, kk, ll);
};

template <typename T1, typename T2, int DIM>
auto get_diff_matrix2(T1 &t_a, T2 &t_d, const FTensor::Number<DIM> &) {
  constexpr auto t_kd = FTensor::Kronecker_Delta<double>();
  FTensor::Tensor4<double, DIM, DIM, DIM, DIM> t_d_a;

  FTensor::Index<'i', DIM> i;
  FTensor::Index<'j', DIM> j;
  FTensor::Index<'k', DIM> k;
  FTensor::Index<'l', DIM> l;

  t_d_a(i, j, k, l) = 0;

  for (int ii = 0; ii != DIM; ++ii)
    for (int jj = 0; jj != DIM; ++jj)
      for (int kk = 0; kk != DIM; ++kk)
        for (int ll = 0; ll != DIM; ++ll)
          for (int zz = 0; zz != DIM; ++zz) {

            auto diff = [&](auto ii, auto jj, auto kk, auto ll, int zz) {
              return

                  t_a(ii, zz) * t_kd(zz, kk) * t_kd(jj, ll)

                  +

                  t_kd(ii, kk) * t_kd(zz, ll) * t_a(zz, jj);
            };

            t_d_a(ii, jj, kk, ll) +=
                (diff(ii, jj, kk, ll, zz) + diff(ii, jj, ll, kk, zz)) / 2.;
          }

  diff_ddg(t_d_a, t_d, FTensor::Number<DIM>());

  return t_d_a;
};

template <typename T1, typename T2, int DIM>
auto get_diff2_matrix2(T1 &t_s, T2 &t_dd, const FTensor::Number<DIM> &) {
  constexpr auto t_kd = FTensor::Kronecker_Delta<double>();
  FTensor::Tensor4<double, DIM, DIM, DIM, DIM> t_dd_a;

  FTensor::Index<'i', DIM> i;
  FTensor::Index<'j', DIM> j;
  FTensor::Index<'k', DIM> k;
  FTensor::Index<'l', DIM> l;

  t_dd_a(i, j, k, l) = 0;

  for (int ii = 0; ii != DIM; ++ii)
    for (int jj = 0; jj != DIM; ++jj)
      for (int kk = 0; kk != DIM; ++kk)
        for (int ll = 0; ll != DIM; ++ll)
          for (int mm = 0; mm != DIM; ++mm)
            for (int nn = 0; nn != DIM; ++nn)
              for (int zz = 0; zz != DIM; ++zz) {

                auto diff = [&](auto ii, auto jj, auto kk, auto ll, int mm,
                                int nn, int zz) {
                  return

                      t_s(ii, jj) *
                      (t_kd(ii, mm) * t_kd(zz, nn) * t_kd(zz, kk) * t_kd(jj, ll)

                       +

                       t_kd(ii, kk) * t_kd(zz, ll) * t_kd(zz, mm) *
                           t_kd(jj, nn));
                };

                t_dd_a(kk, ll, mm, nn) += (

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

  diff_ddg(t_dd_a, t_dd, FTensor::Number<DIM>());

  return t_dd_a;
};

template <typename T1, int DIM>
auto get_diff_matrix(T1 &t_d, const FTensor::Number<DIM> &) {
  constexpr auto t_kd = FTensor::Kronecker_Delta<double>();
  FTensor::Tensor4<double, DIM, DIM, DIM, DIM> t_d_a;

  FTensor::Index<'i', DIM> i;
  FTensor::Index<'j', DIM> j;
  FTensor::Index<'k', DIM> k;
  FTensor::Index<'l', DIM> l;

  t_d_a(i, j, k, l) = 0;

  for (int ii = 0; ii != DIM; ++ii)
    for (int jj = 0; jj != DIM; ++jj)
      for (int kk = 0; kk != DIM; ++kk)
        for (int ll = 0; ll != DIM; ++ll) {

          auto diff = [&](auto ii, auto jj, auto kk, auto ll) {
            return t_kd(ii, kk) * t_kd(jj, ll);
          };

          t_d_a(ii, jj, kk, ll) =
              (diff(ii, jj, kk, ll) + diff(ii, jj, ll, kk)) / 2.;
        }

  diff_ddg(t_d_a, t_d, FTensor::Number<3>());

  return t_d_a;
};

template <typename T, int DIM>
auto get_norm_t4(T &t, const FTensor::Number<DIM> &) {
  double r = 0;
  for (int ii = 0; ii != DIM; ++ii)
    for (int jj = 0; jj != DIM; ++jj)
      for (int kk = 0; kk != DIM; ++kk)
        for (int ll = 0; ll != DIM; ++ll)
          r += t(ii, jj, kk, ll) * t(ii, jj, kk, ll);
  return r;
};

static char help[] = "...\n\n";

int main(int argc, char *argv[]) {

  MoFEM::Core::Initialize(&argc, &argv, (char *)0, help);

  auto core_log = logging::core::get();
  core_log->add_sink(
      LogManager::createSink(LogManager::getStrmSelf(), "ATOM_TEST"));
  LogManager::setLog("ATOM_TEST");
  BOOST_LOG_SCOPED_THREAD_ATTR("Timeline", attrs::timer());
  MOFEM_LOG_ATTRIBUTES("ATOM_TEST", LogManager::BitLineID);

  try {

    auto print_ddg = [](auto &t, auto str = "") {
      constexpr double eps = 1e-6;
      for (int ii = 0; ii != 3; ++ii)
        for (int jj = 0; jj != 3; ++jj)
          for (int kk = 0; kk != 3; ++kk)
            for (int ll = 0; ll != 3; ++ll) {
              double v = t(ii, jj, kk, ll);
              double w = std::abs(v) < eps ? 0 : v;
              MOFEM_LOG("ATOM_TEST", Sev::noisy)
                  << str << std::fixed << std::setprecision(3) << std::showpos
                  << ii + 1 << " " << jj + 1 << " " << kk + 1 << " " << ll + 1
                  << " : " << w;
            }
    };

    auto print_ddg_direction = [](auto &t, auto kk, int ll) {
      for (int ii = 0; ii != 3; ++ii)
        for (int jj = 0; jj <= ii; ++jj)
          MOFEM_LOG("ATOM_TEST", Sev::noisy)
              << ii + 1 << " " << jj + 1 << " " << kk + 1 << " " << ll + 1
              << " : " << t(ii, jj, kk, ll);
    };

    auto print_mat = [](auto &t) {
      for (int ii = 0; ii != 3; ++ii)
        for (int jj = 0; jj != 3; ++jj)
          MOFEM_LOG("ATOM_TEST", Sev::noisy)
              << ii + 1 << " " << jj + 1 << " : " << t(ii, jj);
    };

    enum swap { swap12, swap01 };
    auto run_lapack = [](auto &a, swap s = swap12) {
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

      if (s == swap12) {
        FTensor::Tensor2<double, 3, 3> t_eig_vec{

            a[0 * 3 + 0], a[0 * 3 + 1], a[0 * 3 + 2],

            a[2 * 3 + 0], a[2 * 3 + 1], a[2 * 3 + 2],

            a[1 * 3 + 0], a[1 * 3 + 1], a[1 * 3 + 2]};

        FTensor::Tensor1<double, 3> t_eig_vals{w[0], w[2], w[1]};
        return std::make_tuple(t_a, t_eig_vec, t_eig_vals);
      }

      FTensor::Tensor2<double, 3, 3> t_eig_vec{

          a[1 * 3 + 0], a[1 * 3 + 1], a[1 * 3 + 2],

          a[0 * 3 + 0], a[0 * 3 + 1], a[0 * 3 + 2],

          a[2 * 3 + 0], a[2 * 3 + 1], a[2 * 3 + 2],

      };

      FTensor::Tensor1<double, 3> t_eig_vals{w[1], w[0], w[2]};
      return std::make_tuple(t_a, t_eig_vec, t_eig_vals);
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

      FTensor::Tensor2_symmetric<double, 3> t_S_sym;
      t_S_sym(i, j) = (t_S(i, j) || t_S(j, i)) / 2.;

      MOFEM_LOG("ATOM_TEST", Sev::verbose) << "Diff A";
      print_mat(t_A);

      // Testing against values in mathematica for 0,2 directive
      {
        auto f = [](double v) { return exp(v); };
        auto d_f = [](double v) { return exp(v); };
        auto dd_f = [](double v) { return exp(v); };

        auto t_b = EigenMatrix::getMat(t_L, t_N, f);
        MOFEM_LOG("ATOM_TEST", Sev::verbose) << "Reconstruct mat";
        print_mat(t_b);

        auto t_d = EigenMatrix::getDiffMat(t_L, t_N, f, d_f, 3);

        MOFEM_LOG("ATOM_TEST", Sev::verbose) << "Diff";
        print_ddg_direction(t_d, 0, 2);

        auto t_dd =
            EigenMatrix::getDiffDiffMat(t_L, t_N, f, d_f, dd_f, t_S_sym, 3);

        MOFEM_LOG("ATOM_TEST", Sev::verbose) << "Diff Diff";
        print_ddg_direction(t_dd, 0, 2);

        auto norm2_t_b = t_b(i, j) * t_b(i, j);
        MOFEM_LOG("ATOM_TEST", Sev::inform) << "norm2_t_b " << norm2_t_b;

        auto norm2_t_d = get_norm_t4(t_d, FTensor::Number<3>());
        MOFEM_LOG("ATOM_TEST", Sev::inform) << "norm2_t_d " << norm2_t_d;

        auto norm2_t_dd = get_norm_t4(t_dd, FTensor::Number<3>());
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
        // auto &t_a = std::get<0>(tuple);
        auto &t_eig_vec = std::get<1>(tuple);
        auto &t_eig_vals = std::get<2>(tuple);

        auto t_eig_val_diff =
            (t_eig_vals(i) - t_L(i)) * (t_eig_vals(i) - t_L(i));
        MOFEM_LOG("ATOM_TEST", Sev::inform)
            << "t_eig_val_diff " << t_eig_val_diff;

        auto f = [](double v) { return exp(v); };
        // auto d_f = [](double v) { return exp(v); };
        // auto dd_f = [](double v) { return exp(v); };

        auto t_b = EigenMatrix::getMat(t_L, t_N, f);
        auto t_c = EigenMatrix::getMat(t_eig_vals, t_eig_vec, f);
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
          auto t_b = EigenMatrix::getMat(t_L, t_N, f);
          t_b(i, j) -= (t_A(i, j) || t_A(j, i)) / 2;
          auto norm2_t_b = t_b(i, j) * t_b(i, j);
          MOFEM_LOG("ATOM_TEST", Sev::inform)
              << "Result should be matrix itself " << norm2_t_b;
          if (norm2_t_b > eps)
            SETERRQ(PETSC_COMM_SELF, MOFEM_ATOM_TEST_INVALID,
                    "This norm should be zero");
        }

        {

          auto t_d = EigenMatrix::getDiffMat(t_L, t_N, f, d_f, 3);
          auto t_d_a = get_diff_matrix(t_d, FTensor::Number<3>());

          MOFEM_LOG("ATOM_TEST", Sev::verbose) << "t_d_a";
          print_ddg(t_d_a, "hand ");
          MOFEM_LOG("ATOM_TEST", Sev::verbose) << "t_d";
          print_ddg(t_d, "code ");

          double nrm2_t_d_a = get_norm_t4(t_d_a, FTensor::Number<3>());
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

          FTensor::Tensor2_symmetric<double, 3> t_S_sym;
          t_S_sym(i, j) = t_S(i, j) || t_S(j, i);

          auto t_dd =
              EigenMatrix::getDiffDiffMat(t_L, t_N, f, d_f, dd_f, t_S_sym, 3);

          auto norm2_t_dd = get_norm_t4(t_dd, FTensor::Number<3>());
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
          auto t_b = EigenMatrix::getMat(t_L, t_N, f);
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
          auto t_d = EigenMatrix::getDiffMat(t_L, t_N, f, d_f, 3);
          print_ddg_direction(t_d, 0, 2);
          auto t_d_a = get_diff_matrix2(t_A, t_d, FTensor::Number<3>());
          MOFEM_LOG("ATOM_TEST", Sev::verbose) << "t_d_a";
          print_ddg(t_d_a, "hand ");
          MOFEM_LOG("ATOM_TEST", Sev::verbose) << "t_d";
          print_ddg(t_d, "code ");
          double nrm2_t_d_a = get_norm_t4(t_d_a, FTensor::Number<3>());
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

          FTensor::Tensor2_symmetric<double, 3> t_S_sym;
          t_S_sym(i, j) = t_S(i, j) || t_S(j, i);

          auto t_dd =
              EigenMatrix::getDiffDiffMat(t_L, t_N, f, d_f, dd_f, t_S_sym, 3);
          auto t_dd_a = get_diff2_matrix2(t_S_sym, t_dd, FTensor::Number<3>());

          MOFEM_LOG("ATOM_TEST", Sev::verbose) << "t_dd_a";
          print_ddg(t_dd_a, "hand ");
          MOFEM_LOG("ATOM_TEST", Sev::verbose) << "t_dd";
          print_ddg(t_dd, "code ");

          double nrm2_t_dd_a = get_norm_t4(t_dd_a, FTensor::Number<3>());
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

      std::array<double, 9> a{5.,  4., 0,

                              4.,  5,  0.,

                              0.0, 0., 9};

      auto tuple = run_lapack(a, swap01);
      auto &t_a = std::get<0>(tuple);
      auto &t_eig_vecs = std::get<1>(tuple);
      auto &t_eig_vals = std::get<2>(tuple);

      auto f = [](double v) { return v; };
      auto d_f = [](double v) { return 1; };
      // auto dd_f = [](double v) { return 0; };

      constexpr double eps = 1e-10;
      {
        auto t_b = EigenMatrix::getMat(t_eig_vals, t_eig_vecs, f);
        t_b(i, j) -= (t_a(i, j) || t_a(j, i)) / 2;
        auto norm2_t_b = t_b(i, j) * t_b(i, j);
        MOFEM_LOG("ATOM_TEST", Sev::inform)
            << "Result should be matrix itself " << norm2_t_b;
        if (norm2_t_b > eps)
          SETERRQ(PETSC_COMM_SELF, MOFEM_ATOM_TEST_INVALID,
                  "This norm should be zero");
      }

      {
        auto t_d = EigenMatrix::getDiffMat(t_eig_vals, t_eig_vecs, f, d_f, 2);
        auto t_d_a = get_diff_matrix(t_d, FTensor::Number<3>());
        MOFEM_LOG("ATOM_TEST", Sev::verbose) << "t_d_a";
        print_ddg(t_d_a, "hand ");
        MOFEM_LOG("ATOM_TEST", Sev::verbose) << "t_d";
        print_ddg(t_d, "code ");
        double nrm2_t_d_a = get_norm_t4(t_d_a, FTensor::Number<3>());
        MOFEM_LOG("ATOM_TEST", Sev::inform)
            << "Direvarive hand calculation minus code " << nrm2_t_d_a;
        if (nrm2_t_d_a > eps)
          SETERRQ(PETSC_COMM_SELF, MOFEM_ATOM_TEST_INVALID,
                  "This norm should be zero");
      }

      {
        auto f = [](double v) { return v * v; };
        auto d_f = [](double v) { return 2 * v; };
        // auto dd_f = [](double v) { return 2; };
        auto t_d = EigenMatrix::getDiffMat(t_eig_vals, t_eig_vecs, f, d_f, 2);
        auto t_d_a = get_diff_matrix2(t_a, t_d, FTensor::Number<3>());
        MOFEM_LOG("ATOM_TEST", Sev::verbose) << "t_d_a";
        print_ddg(t_d_a, "hand ");
        MOFEM_LOG("ATOM_TEST", Sev::verbose) << "t_d";
        print_ddg(t_d, "code ");
        double nrm2_t_d_a = get_norm_t4(t_d_a, FTensor::Number<3>());
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
      // auto dd_f = [](double v) { return 0; };

      auto tuple = run_lapack(a);
      auto &t_a = std::get<0>(tuple);
      auto &t_eig_vecs = std::get<1>(tuple);
      auto &t_eig_vals = std::get<2>(tuple);

      constexpr double eps = 1e-10;
      {
        auto t_b = EigenMatrix::getMat(t_eig_vals, t_eig_vecs, f);
        t_b(i, j) -= (t_a(i, j) || t_a(j, i)) / 2;
        auto norm2_t_b = t_b(i, j) * t_b(i, j);
        MOFEM_LOG("ATOM_TEST", Sev::inform)
            << "Result should be matrix itself " << norm2_t_b;
        if (norm2_t_b > eps)
          SETERRQ(PETSC_COMM_SELF, MOFEM_ATOM_TEST_INVALID,
                  "This norm should be zero");
      }

      {
        auto t_d = EigenMatrix::getDiffMat(t_eig_vals, t_eig_vecs, f, d_f, 1);
        auto t_d_a = get_diff_matrix(t_d, FTensor::Number<3>());
        MOFEM_LOG("ATOM_TEST", Sev::verbose) << "t_d_a";
        print_ddg(t_d_a, "hand ");
        MOFEM_LOG("ATOM_TEST", Sev::verbose) << "t_d";
        print_ddg(t_d, "code ");
        double nrm2_t_d_a = get_norm_t4(t_d_a, FTensor::Number<3>());
        MOFEM_LOG("ATOM_TEST", Sev::inform)
            << "Direvarive hand calculation minus code " << nrm2_t_d_a;
        if (nrm2_t_d_a > eps)
          SETERRQ(PETSC_COMM_SELF, MOFEM_ATOM_TEST_INVALID,
                  "This norm should be zero");
      }

      {
        auto f = [](double v) { return v * v; };
        auto d_f = [](double v) { return 2 * v; };
        // auto dd_f = [](double v) { return 2; };
        auto t_d = EigenMatrix::getDiffMat(t_eig_vals, t_eig_vecs, f, d_f, 1);
        auto t_d_a = get_diff_matrix2(t_a, t_d, FTensor::Number<3>());
        MOFEM_LOG("ATOM_TEST", Sev::verbose) << "t_d_a";
        print_ddg(t_d_a, "hand ");
        MOFEM_LOG("ATOM_TEST", Sev::verbose) << "t_d";
        print_ddg(t_d, "code ");
        double nrm2_t_d_a = get_norm_t4(t_d_a, FTensor::Number<3>());
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
      // auto &t_a = std::get<0>(tuple);
      auto &t_eig_vecs = std::get<1>(tuple);
      auto &t_eig_vals = std::get<2>(tuple);

      t_eig_vals(0) -= 1e-4;
      t_eig_vals(2) += 1e-4;

      constexpr double eps = 1e-10;

      auto f = [](double v) { return v; };
      auto d_f = [](double v) { return 1; };
      auto dd_f = [](double v) { return 0; };

      FTensor::Tensor2<double, 3, 3> t_S{

          1.,      1. / 2., 1. / 3.,

          2. / 2., 1.,      2. / 3.,

          3. / 2., 1.,      3. / 3.};

      FTensor::Tensor2_symmetric<double, 3> t_S_sym;
      t_S_sym(i, j) = t_S(i, j) || t_S(j, i);

      auto t_dd = EigenMatrix::getDiffDiffMat(t_eig_vals, t_eig_vecs, f, d_f,
                                              dd_f, t_S_sym, 1);

      MOFEM_LOG("ATOM_TEST", Sev::verbose) << "t_dd";
      print_ddg(t_dd, "test ");

      double nrm2_t_dd = get_norm_t4(t_dd, FTensor::Number<3>());
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
      // auto &t_a = std::get<0>(tuple);
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

      FTensor::Tensor2_symmetric<double, 3> t_S_sym;
      t_S_sym(i, j) = t_S(i, j) || t_S(j, i);

      auto t_dd = EigenMatrix::getDiffDiffMat(t_eig_vals, t_eig_vecs, f, d_f,
                                              dd_f, t_S_sym, 1);
      // print_ddg(t_dd, "test ");

      auto t_dd_a = get_diff2_matrix2(t_S_sym, t_dd, FTensor::Number<3>());

      double nrm2_t_dd_a = get_norm_t4(t_dd_a, FTensor::Number<3>());
      MOFEM_LOG("ATOM_TEST", Sev::inform)
          << "Direvarive hand calculation minus code " << nrm2_t_dd_a;
      if (nrm2_t_dd_a > eps)
        SETERRQ(PETSC_COMM_SELF, MOFEM_ATOM_TEST_INVALID,
                "This norm should be zero");
    }

    // check second directive two reapeating eiegn values
    {

      std::array<double, 9> a{5., 4., 0.,

                              4., 5., 0.,

                              0., 0., 9};

      auto tuple = run_lapack(a, swap01);
      // auto &t_a = std::get<0>(tuple);
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

      FTensor::Tensor2_symmetric<double, 3> t_S_sym;
      t_S_sym(i, j) = t_S(i, j) || t_S(j, i);

      auto t_dd = EigenMatrix::getDiffDiffMat(t_eig_vals, t_eig_vecs, f, d_f,
                                              dd_f, t_S_sym, 2);
      print_ddg(t_dd, "test ");

      auto t_dd_a = get_diff2_matrix2(t_S_sym, t_dd, FTensor::Number<3>());

      double nrm2_t_dd_a = get_norm_t4(t_dd_a, FTensor::Number<3>());
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
      // auto &t_a = std::get<0>(tuple);
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

      FTensor::Tensor2_symmetric<double, 3> t_S_sym;
      t_S_sym(i, j) = t_S(i, j) || t_S(j, i);

      auto t_dd_1 = EigenMatrix::getDiffDiffMat(t_eig_vals, t_eig_vecs, f, d_f,
                                                dd_f, t_S_sym, 3);
      auto t_dd_2 = EigenMatrix::getDiffDiffMat(t_eig_vals, t_eig_vecs, f, d_f,
                                                dd_f, t_S_sym, 1);

      double nrm2_t_dd_t1 = get_norm_t4(t_dd_1, FTensor::Number<3>());
      MOFEM_LOG("ATOM_TEST", Sev::verbose)
          << "Direvarive nor t_dd_1 " << nrm2_t_dd_t1;

      double nrm2_t_dd_t2 = get_norm_t4(t_dd_2, FTensor::Number<3>());
      MOFEM_LOG("ATOM_TEST", Sev::verbose)
          << "Direvarive norm t_dd_2 " << nrm2_t_dd_t2;

      print_ddg(t_dd_1, "t_dd_1 ");
      print_ddg(t_dd_2, "t_dd_2 ");

      FTensor::Ddg<double, 3, 3> t_dd_3;
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

      double nrm2_t_dd_3 = get_norm_t4(t_dd_3, FTensor::Number<3>());
      MOFEM_LOG("ATOM_TEST", Sev::inform)
          << "Direvarive approx. calculation minus code " << nrm2_t_dd_3;
      if (nrm2_t_dd_3 > eps)
        SETERRQ(PETSC_COMM_SELF, MOFEM_ATOM_TEST_INVALID,
                "This norm should be zero");
    }

    // check second directive exponent agains perturned
    {

      std::array<double, 9> a{5., 4., 0.,

                              4., 5., 0.,

                              0., 0., 9};

      auto tuple = run_lapack(a, swap01);
      // auto &t_a = std::get<0>(tuple);
      auto &t_eig_vecs = std::get<1>(tuple);
      auto &t_eig_vals = std::get<2>(tuple);

      t_eig_vals(0) -= 1e-4;
      t_eig_vals(2) += 1e-4;

      constexpr double eps = 1e-4;

      auto f = [](double v) { return v * v; };
      auto d_f = [](double v) { return 2 * v; };
      auto dd_f = [](double v) { return 2; };

      FTensor::Tensor2<double, 3, 3> t_S{

          1.,      1. / 2., 1. / 3.,

          2. / 1., 1.,      2. / 3.,

          3. / 1., 3. / 1., 1.};

      FTensor::Tensor2_symmetric<double, 3> t_S_sym;
      t_S_sym(i, j) = t_S(i, j) || t_S(j, i);

      auto t_dd_1 = EigenMatrix::getDiffDiffMat(t_eig_vals, t_eig_vecs, f, d_f,
                                                dd_f, t_S_sym, 3);
      auto t_dd_2 = EigenMatrix::getDiffDiffMat(t_eig_vals, t_eig_vecs, f, d_f,
                                                dd_f, t_S_sym, 2);

      double nrm2_t_dd_t1 = get_norm_t4(t_dd_1, FTensor::Number<3>());
      MOFEM_LOG("ATOM_TEST", Sev::verbose)
          << "Direvarive nor t_dd_1 " << nrm2_t_dd_t1;

      double nrm2_t_dd_t2 = get_norm_t4(t_dd_2, FTensor::Number<3>());
      MOFEM_LOG("ATOM_TEST", Sev::verbose)
          << "Direvarive norm t_dd_2 " << nrm2_t_dd_t2;

      print_ddg(t_dd_1, "t_dd_1 ");
      print_ddg(t_dd_2, "t_dd_2 ");

      FTensor::Ddg<double, 3, 3> t_dd_3;
      t_dd_3(i, j, k, l) = t_dd_1(i, j, k, l) - t_dd_2(i, j, k, l);

      for (int ii = 0; ii != 3; ++ii)
        for (int jj = 0; jj != 3; ++jj)
          for (int kk = 0; kk != 3; ++kk)
            for (int ll = 0; ll != 3; ++ll) {
              constexpr double eps = 1e-3;
              if (std::abs(t_dd_3(ii, jj, kk, ll)) > eps)
                MOFEM_LOG("ATOM_TEST", Sev::error)
                    << "Error " << ii << " " << jj << " " << kk << " " << ll
                    << " " << t_dd_1(ii, jj, kk, ll) << " "
                    << t_dd_2(ii, jj, kk, ll);
            }

      double nrm2_t_dd_3 = get_norm_t4(t_dd_3, FTensor::Number<3>());
      MOFEM_LOG("ATOM_TEST", Sev::inform)
          << "Direvarive approx. calculation minus code " << nrm2_t_dd_3;
      if (nrm2_t_dd_3 > eps)
        SETERRQ(PETSC_COMM_SELF, MOFEM_ATOM_TEST_INVALID,
                "This norm should be zero");
    }

    // check second directive exponent
    {

      std::array<double, 9> a{5., 4., 0.,

                              4., 5., 0.,

                              0., 0., 9};

      auto tuple = run_lapack(a, swap01);
      // auto &t_a = std::get<0>(tuple);
      auto &t_eig_vecs = std::get<1>(tuple);
      auto &t_eig_vals = std::get<2>(tuple);

      constexpr double eps = 1e-4;
      constexpr int p = 3;

      auto f = [](double v) { return pow(v, p); };
      auto d_f = [](double v) { return p * pow(v, p - 1); };
      auto dd_f = [](double v) {
        return p * (p - 1) * pow(v, std::max(0, p - 2));
      };

      FTensor::Tensor2<double, 3, 3> t_S{

          1.,      1. / 2., 1. / 3.,

          2. / 1., 1.,      2. / 3.,

          3. / 1., 3. / 1., 1.};

      // FTensor::Tensor2_symmetric<double, 3> t_S_sym;
      // t_S_sym(i, j) = t_S(i, j) || t_S(j, i);

      t_eig_vals(0) += 2e-5;
      t_eig_vals(2) -= 2e-5;
      auto t_dd_1 = EigenMatrix::getDiffDiffMat(t_eig_vals, t_eig_vecs, f, d_f,
                                                dd_f, t_S, 3);
      auto t_dd_2 = EigenMatrix::getDiffDiffMat(t_eig_vals, t_eig_vecs, f, d_f,
                                                dd_f, t_S, 2);

      double nrm2_t_dd_t1 = get_norm_t4(t_dd_1, FTensor::Number<3>());
      MOFEM_LOG("ATOM_TEST", Sev::verbose)
          << "Direvarive nor t_dd_1 " << nrm2_t_dd_t1;

      double nrm2_t_dd_t2 = get_norm_t4(t_dd_2, FTensor::Number<3>());
      MOFEM_LOG("ATOM_TEST", Sev::verbose)
          << "Direvarive norm t_dd_2 " << nrm2_t_dd_t2;

      print_ddg(t_dd_1, "t_dd_1 ");
      print_ddg(t_dd_2, "t_dd_2 ");

      FTensor::Ddg<double, 3, 3> t_dd_3;
      t_dd_3(i, j, k, l) = t_dd_1(i, j, k, l) - t_dd_2(i, j, k, l);

      for (int ii = 0; ii != 3; ++ii)
        for (int jj = 0; jj != 3; ++jj)
          for (int kk = 0; kk != 3; ++kk)
            for (int ll = 0; ll != 3; ++ll) {
              constexpr double eps = 1e-3;
              if (std::abs(t_dd_3(ii, jj, kk, ll)) > eps)
                MOFEM_LOG("ATOM_TEST", Sev::error)
                    << "Error " << ii << " " << jj << " " << kk << " " << ll
                    << " " << t_dd_1(ii, jj, kk, ll) << " "
                    << t_dd_2(ii, jj, kk, ll) << " " << t_dd_3(ii, jj, kk, ll);
            }

      double nrm2_t_dd_3 = get_norm_t4(t_dd_3, FTensor::Number<3>());
      MOFEM_LOG("ATOM_TEST", Sev::inform)
          << "Direvarive approx. calculation minus code " << nrm2_t_dd_3;
      if (nrm2_t_dd_3 > eps)
        SETERRQ(PETSC_COMM_SELF, MOFEM_ATOM_TEST_INVALID,
                "This norm should be zero");
    }

    // Speed
    {

      std::array<double, 9> a{1.,   0.1, -0.5,

                              0.1,  2.,  0.,

                              -0.5, 0.,  3.};

      auto tuple = run_lapack(a);
      // auto &t_a = std::get<0>(tuple);
      auto &t_eig_vecs = std::get<1>(tuple);
      auto &t_eig_vals = std::get<2>(tuple);

      auto f = [](double v) { return exp(v); };
      auto d_f = [](double v) { return exp(v); };
      auto dd_f = [](double v) { return exp(v); };

      FTensor::Tensor2<double, 3, 3> t_S{

          1.,      1. / 2., 1. / 3.,

          2. / 1., 1.,      2. / 3.,

          3. / 1., 3. / 1., 1.};

      FTensor::Tensor2_symmetric<double, 3> t_S_sym;
      t_S_sym(i, j) = t_S(i, j) || t_S(j, i);

      MOFEM_LOG("ATOM_TEST", Sev::inform) << "Start";
      for (int ii = 0; ii != 1000; ++ii) {
        auto t_d = EigenMatrix::getDiffMat(t_eig_vals, t_eig_vecs, f, d_f, 3);
        auto t_dd = EigenMatrix::getDiffDiffMat(t_eig_vals, t_eig_vecs, f, d_f,
                                                dd_f, t_S_sym, 3);
        std::ignore = t_d;
        std::ignore = t_dd;
      }
      MOFEM_LOG("ATOM_TEST", Sev::inform) << "End";
    }

    // 2d case

    auto run_lapack_2d = [](auto &a) {
      int info;
      double wkopt;
      double w[2];

      FTensor::Tensor2<double, 2, 2> t_a{

          a[0], a[1],

          a[2], a[3]};

      /* Query and allocate the optimal workspace */
      int lwork = -1;
      info = lapack_dsyev('V', 'U', 2, a.data(), 2, w, &wkopt, lwork);
      if (info > 0)
        THROW_MESSAGE("The algorithm failed to compute eigenvalues.");
      lwork = (int)wkopt;
      std::vector<double> work(lwork);
      /* Solve eigenproblem */
      info = lapack_dsyev('V', 'U', 2, a.data(), 2, w, &*work.begin(), lwork);
      if (info > 0)
        THROW_MESSAGE("The algorithm failed to compute eigenvalues.");

      FTensor::Tensor2<double, 2, 2> t_eig_vecs{

          a[0 * 2 + 0], a[0 * 2 + 1],

          a[1 * 2 + 0], a[1 * 2 + 1]};

      FTensor::Tensor1<double, 2> t_eig_vals{w[0], w[1]};

      return std::make_tuple(t_a, t_eig_vecs, t_eig_vals);
    };

    // Testsing quadratic function for 2d
    {

      std::array<double, 9> a{1., 0.1,

                              0.1, 2.};

      auto tuple = run_lapack_2d(a);
      auto &t_A = std::get<0>(tuple);
      auto &t_eig_vecs = std::get<1>(tuple);
      auto &t_eig_vals = std::get<2>(tuple);

      auto f = [](double v) { return v * v; };
      auto d_f = [](double v) { return 2 * v; };
      auto dd_f = [](double v) { return 2; };

      constexpr double eps = 1e-6;

      FTensor::Index<'i', 2> i;
      FTensor::Index<'j', 2> j;
      FTensor::Index<'k', 2> k;
      FTensor::Index<'l', 2> l;

      // check if multiplication gives right value
      {
        auto t_b = EigenMatrix::getMat(t_eig_vals, t_eig_vecs, f);
        FTensor::Tensor2<double, 2, 2> t_a;
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
        auto t_d = EigenMatrix::getDiffMat(t_eig_vals, t_eig_vecs, f, d_f, 2);
        auto t_d_a = get_diff_matrix2(t_A, t_d, FTensor::Number<2>());
        double nrm2_t_d_a = get_norm_t4(t_d_a, FTensor::Number<2>());
        MOFEM_LOG("ATOM_TEST", Sev::inform)
            << "Direvarive hand calculation minus code " << nrm2_t_d_a;
        if (nrm2_t_d_a > eps)
          SETERRQ(PETSC_COMM_SELF, MOFEM_ATOM_TEST_INVALID,
                  "This norm should be zero");
      }

      // check second directive
      {
        FTensor::Tensor2<double, 2, 2> t_S{

            1., 1. / 2,

            2. / 2., 1.};

        // FTensor::Index<'i', 2> i;
        // FTensor::Index<'j', 2> j;
        // FTensor::Tensor2_symmetric<double, 2> t_S_sym;
        // t_S_sym(i, j) = t_S(i, j) || t_S(j, i);

        auto t_dd = EigenMatrix::getDiffDiffMat(t_eig_vals, t_eig_vecs, f, d_f,
                                                dd_f, t_S, 2);
        auto t_dd_a = get_diff2_matrix2(t_S, t_dd, FTensor::Number<2>());

        MOFEM_LOG("ATOM_TEST", Sev::verbose) << "t_dd_a";
        MOFEM_LOG("ATOM_TEST", Sev::verbose) << "t_dd";

        double nrm2_t_dd_a = get_norm_t4(t_dd_a, FTensor::Number<2>());
        MOFEM_LOG("ATOM_TEST", Sev::inform)
            << "Direvarive hand calculation minus code " << nrm2_t_dd_a;
        if (nrm2_t_dd_a > eps)
          SETERRQ(PETSC_COMM_SELF, MOFEM_ATOM_TEST_INVALID,
                  "This norm should be zero");
      }
    }

    // Testsing quadratic function for repeating eigen valsues
    {

      std::array<double, 9> a{2., 0,

                              0, 2.};

      auto tuple = run_lapack_2d(a);
      auto &t_A = std::get<0>(tuple);
      auto &t_eig_vecs = std::get<1>(tuple);
      auto &t_eig_vals = std::get<2>(tuple);

      auto f = [](double v) { return v * v; };
      auto d_f = [](double v) { return 2 * v; };
      auto dd_f = [](double v) { return 2; };

      constexpr double eps = 1e-6;

      FTensor::Index<'i', 2> i;
      FTensor::Index<'j', 2> j;
      FTensor::Index<'k', 2> k;
      FTensor::Index<'l', 2> l;

      // check if multiplication gives right value
      {
        auto t_b = EigenMatrix::getMat(t_eig_vals, t_eig_vecs, f);
        FTensor::Tensor2<double, 2, 2> t_a;
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
        auto t_d = EigenMatrix::getDiffMat(t_eig_vals, t_eig_vecs, f, d_f, 1);
        auto t_d_a = get_diff_matrix2(t_A, t_d, FTensor::Number<2>());
        double nrm2_t_d_a = get_norm_t4(t_d_a, FTensor::Number<2>());
        MOFEM_LOG("ATOM_TEST", Sev::inform)
            << "Direvarive hand calculation minus code " << nrm2_t_d_a;
        if (nrm2_t_d_a > eps)
          SETERRQ(PETSC_COMM_SELF, MOFEM_ATOM_TEST_INVALID,
                  "This norm should be zero");
      }

      // check second directive
      {
        FTensor::Tensor2<double, 2, 2> t_S{

            1., 1. / 2,

            2. / 2., 1.};

        FTensor::Index<'i', 2> i;
        FTensor::Index<'j', 2> j;
        FTensor::Tensor2_symmetric<double, 2> t_S_sym;
        t_S_sym(i, j) = t_S(i, j) || t_S(j, i);

        auto t_dd = EigenMatrix::getDiffDiffMat(t_eig_vals, t_eig_vecs, f, d_f,
                                                dd_f, t_S_sym, 1);
        auto t_dd_a = get_diff2_matrix2(t_S_sym, t_dd, FTensor::Number<2>());

        MOFEM_LOG("ATOM_TEST", Sev::verbose) << "t_dd_a";
        MOFEM_LOG("ATOM_TEST", Sev::verbose) << "t_dd";

        double nrm2_t_dd_a = get_norm_t4(t_dd_a, FTensor::Number<2>());
        MOFEM_LOG("ATOM_TEST", Sev::inform)
            << "Direvarive hand calculation minus code " << nrm2_t_dd_a;
        if (nrm2_t_dd_a > eps)
          SETERRQ(PETSC_COMM_SELF, MOFEM_ATOM_TEST_INVALID,
                  "This norm should be zero");
      }
    }
  }
  CATCH_ERRORS;

  CHKERR MoFEM::Core::Finalize();
}
