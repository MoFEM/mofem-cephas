#include "../../../src/FTensor.hpp"
#include "../test_for_zero.hpp"
#include <iostream>
using namespace FTensor;
using namespace std;

void test_T4ddg_006(const Tensor2_symmetric<double, 3> &t2s_2,
                    const Tensor2_symmetric<double, 3> &t2s_3) {
  Index<'i', 3> i;
  Index<'j', 3> j;
  Index<'k', 3> k;
  Index<'l', 3> l;

  {
    Ddg<double, 3, 3> t4ddg_3_1, t4ddg_3_2, t4ddg_3_3;
    t4ddg_3_1(i, j, k, l) = t2s_2(i, j) * t2s_3(k, l);
    t4ddg_3_2(i, j, k, l) = t2s_3(i, j) * t2s_2(k, l);
    t4ddg_3_3(i, j, k, l) = t4ddg_3_1(i, j, k, l);
    t4ddg_3_3(i, j, k, l) += t4ddg_3_2(i, j, k, l);

    for (int ii = 0; ii != 3; ++ii)
      for (int jj = 0; jj != 3; ++jj)
        for (int kk = 0; kk != 3; ++kk)
          for (int ll = 0; ll != 3; ++ll) {
            test_for_zero(t4ddg_3_3(ii, jj, kk, ll) -
                              t4ddg_3_1(ii, jj, kk, ll) -
                              t4ddg_3_2(ii, jj, kk, ll),
                          "T4ddg(i,j,k,l)+=T4ddg(i,j,k,l)(" + to_string(ii) +
                              "," + to_string(jj) + "," + to_string(kk) + "," +
                              to_string(ll) + ")");
          }

    Ddg<double, 3, 3> t4ddg_3_4;
    t4ddg_3_4(i, j, k, l) = t4ddg_3_3(i, j, k, l);
    t4ddg_3_4(i, j, k, l) -= t4ddg_3_3(i, j, k, l);
    for (int ii = 0; ii != 3; ++ii)
      for (int jj = 0; jj != 3; ++jj)
        for (int kk = 0; kk != 3; ++kk)
          for (int ll = 0; ll != 3; ++ll) {
            test_for_zero(t4ddg_3_4(ii, jj, kk, ll),
                          "T4ddg(i,j,k,l)-=T4ddg(i,j,k,l)(" + to_string(ii) +
                              "," + to_string(jj) + "," + to_string(kk) + "," +
                              to_string(ll) + ")");
          }
  }

  {
    Ddg<double, 3, 3> t4ddg_3_1, t4ddg_3_2;
    t4ddg_3_1(i, j, k, l) = t2s_2(i, j) * t2s_3(k, l);
    t4ddg_3_2(i, j, k, l) = t4ddg_3_1(i, j, k, l);
    t4ddg_3_1(i, j, k, l) *= 2.;
    for (int ii = 0; ii != 3; ++ii)
      for (int jj = 0; jj != 3; ++jj)
        for (int kk = 0; kk != 3; ++kk)
          for (int ll = 0; ll != 3; ++ll) {
            test_for_zero(
                t4ddg_3_1(ii, jj, kk, ll) - 2 * t4ddg_3_2(ii, jj, kk, ll),
                "T4ddg(i,j,k,l)*=2(" + to_string(ii) + "," + to_string(jj) +
                    "," + to_string(kk) + "," + to_string(ll) + ")");
          }
  }

  {
    Ddg<double, 3, 3> t4ddg_3_1;
    t4ddg_3_1(i, j, k, l) = t2s_2(i, j) * t2s_3(k, l);

    Tensor2<double, 3, 3> t2;
    for (int ii = 0; ii != 3; ++ii)
      for (int jj = 0; jj != 3; ++jj)
        t2(ii, jj) = 1 + ii + 10 * jj;

    Tensor2<double, 3, 3> t2_1, t2_2;
    t2_1(i, k) = t4ddg_3_1(i, j, k, l) * t2(j, l);
    t2_2(i, k) = t2s_2(i, j) * (t2s_3(k, l) * t2(j, l));

    for (int ii = 0; ii != 3; ++ii)
      for (int jj = 0; jj != 3; ++jj) {
        test_for_zero(t2_1(ii, jj) - t2_2(ii, jj), "T4ddg(i,j,k,l)*T2(j,l)(" +
                                                       to_string(ii) + "," +
                                                       to_string(jj) + ")");
      }
  }

  {
    Ddg<double, 3, 3> t4ddg_3_1;
    t4ddg_3_1(i, j, k, l) = t2s_2(i, j) * t2s_3(k, l);

    Tensor2<double, 3, 3> t2;
    for (int ii = 0; ii != 3; ++ii)
      for (int jj = 0; jj != 3; ++jj)
        t2(ii, jj) = 1 + ii + 10 * jj;

    Tensor2_symmetric<double, 3> t2s_1, t2s_4;
    t2s_1(k, l) = t4ddg_3_1(i, j, k, l) * t2(i, j);

    t2s_4(k, l) = 0;
    for (int kk = 0; kk != 3; ++kk)
      for (int ll = kk; ll != 3; ++ll)
        for (int ii = 0; ii != 3; ++ii)
          for (int jj = 0; jj != 3; ++jj)
            t2s_4(kk, ll) += t2s_3(kk, ll) * (t2s_2(ii, jj) * t2(ii, jj));

    for (int ii = 0; ii != 3; ++ii)
      for (int jj = 0; jj != 3; ++jj) {
        test_for_zero(t2s_1(ii, jj) - t2s_4(ii, jj), "T4ddg(i,j,k,l)*T2(i,j)(" +
                                                         to_string(ii) + "," +
                                                         to_string(jj) + ")");
      }
  }
}