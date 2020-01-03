#include "../../../src/FTensor.hpp"
#include "../test_for_zero.hpp"
#include <iostream>
using namespace FTensor;
using namespace std;

void test_T4ddg_008(const Tensor2<double, 3, 3> &t2_1,
                    const Tensor2<double, 3, 3> &t2_2,
                    const Tensor2_symmetric<double, 3> &t2s_2,
                    const Tensor2_symmetric<double, 3> &t2s_3) {

  Index<'i', 3> i;
  Index<'j', 3> j;
  Index<'k', 3> k;
  Index<'l', 3> l;
  Index<'m', 3> m;
  Index<'n', 3> n;

  Ddg<double, 3, 3> t4ddg_1;
  t4ddg_1(i, j, k, l) = t2s_2(i, j) * t2s_3(k, l);
  Tensor4<double, 3, 3, 3, 3> t4_1;
  t4_1(i, j, k, l) = t2_1(i, j) * t2_2(k, l);

  {

    Tensor4<double, 3, 3, 3, 3> t4_2, t4_3, t4_4;

    t4_2(i, j, k, l) = t4ddg_1(i, j, m, n) * t4_1(m, n, k, l);

    t4_3(i, j, k, l) = 0;
    for (int ii = 0; ii != 3; ++ii) {
      for (int jj = 0; jj != 3; ++jj) {
        for (int kk = 0; kk != 3; ++kk) {
          for (int ll = 0; ll != 3; ++ll) {
            for (int mm = 0; mm != 3; ++mm) {
              for (int nn = 0; nn != 3; ++nn) {
                t4_3(ii, jj, kk, ll) +=
                    t4ddg_1(ii, jj, mm, nn) * t4_1(mm, nn, kk, ll);
              }
            }
          }
        }
      }
    }

    for (int ii = 0; ii != 3; ++ii)
      for (int jj = 0; jj != 3; ++jj)
        for (int kk = 0; kk != 3; ++kk)
          for (int ll = 0; ll != 3; ++ll) {
            test_for_zero(t4_2(ii, jj, kk, ll) - t4_3(ii, jj, kk, ll),
                          "t4_2(i, j, k, l) = t4ddg_1(i, j, m, n) * "
                          "t4_1(m, n, k, l)");
          }
  
      t4_4(i, j, k, l) = t4_1(m, n, k, l) * t4ddg_1(i, j, m, n);


    for (int ii = 0; ii != 3; ++ii)
      for (int jj = 0; jj != 3; ++jj)
        for (int kk = 0; kk != 3; ++kk)
          for (int ll = 0; ll != 3; ++ll) {
            test_for_zero(t4_4(ii, jj, kk, ll) - t4_3(ii, jj, kk, ll),
                          "t4_2(i, j, k, l) = t4ddg_1(i, j, m, n) * "
                          "t4_1(m, n, k, l)");
          }

  
  }

};