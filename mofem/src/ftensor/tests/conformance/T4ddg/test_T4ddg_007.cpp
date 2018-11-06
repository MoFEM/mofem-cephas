#include "../../../src/FTensor.hpp"
#include "../test_for_zero.hpp"
#include <iostream>
using namespace FTensor;
using namespace std;

void test_T4ddg_007(const Tensor2_symmetric<double, 3> &t2s_2,
                    const Tensor2_symmetric<double, 3> &t2s_3) {
  Index<'i', 3> i;
  Index<'j', 3> j;
  Index<'k', 3> k;
  Index<'l', 3> l;

  Ddg<double, 3, 3> t4ddg_3_1;
  t4ddg_3_1(i, j, k, l) = t2s_2(i, j) * t2s_3(k, l);

  Tensor1<double, 3> t1_1;
  t1_1(0) = 1;
  t1_1(1) = 10;
  t1_1(2) = 100;

  {
    Dg<double, 3, 3> t3dg_3_1, t3dg_3_2;
    t3dg_3_1(i, j, l) = t4ddg_3_1(i, j, k, l) * t1_1(k);
    t3dg_3_2(i, j, l) = t2s_2(i, j) * (t2s_3(k, l) * t1_1(k));

    for (int ii = 0; ii != 3; ++ii)
      for (int jj = 0; jj != 3; ++jj)
        for (int ll = 0; ll != 3; ++ll) {
          test_for_zero(t3dg_3_1(ii, jj, ll) - t3dg_3_2(ii, jj, ll),
                        "t4ddg_3_1(i, j, k, l) * t1_1(k)(" + to_string(ii) +
                            "," + to_string(jj) + "," + to_string(ll) + ")");
        }
  }

  {
    Christof<double, 3, 3> t3ch_3_1, t3ch_3_2;
    t3ch_3_1(l, j, k) = t4ddg_3_1(i, j, k, l) * t1_1(i);
    t3ch_3_2(i, j, k) = 0;
    for (int jj = 0; jj != 3; ++jj)
      for (int kk = 0; kk != 3; ++kk)
        for (int ll = kk; ll != 3; ++ll) {
          auto &v = t3ch_3_2(jj, kk, ll);
          for (int ii = 0; ii != 3; ++ii) {
            v += (t2s_2(ii, jj) * t1_1(ii)) * t2s_3(kk, ll);
          }
        }

    for (int jj = 0; jj != 3; ++jj)
      for (int kk = 0; kk != 3; ++kk)
        for (int ll = 0; ll != 3; ++ll) {
          test_for_zero(t3ch_3_1(jj, kk, ll) - t3ch_3_2(jj, kk, ll),
                        "t4ddg_3_1(i, j, k, l) * t1_1(i)(" + to_string(jj) +
                            "," + to_string(kk) + "," + to_string(ll) + ")");
        }
  }
}