#include "../../../src/FTensor.hpp"
#include "../test_for_zero.hpp"
#include <iostream>
using namespace FTensor;
using namespace std;

void test_T4_008() {

  Index<'i', 3> i;
  Index<'j', 3> j;
  Index<'k', 3> k;
  Index<'l', 3> l;

  Tensor4<double, 3, 3, 3, 3> t_4;

  for (int ii = 0; ii != 3;++ii)
    for (int jj = 0; jj != 3;++jj)
      for (int kk = 0; kk != 3;++kk)
        for (int ll = 0; ll != 3;++ll)
          t_4(ii, jj, kk, ll) = 1 + ii + 10 * jj + 100 * kk + 1000 * ll;

  Tensor2<double, 3, 3> t_2;
  Tensor2_symmetric<double, 3> ts_2;
  for (int ii = 0; ii != 3;++ii)
    for (int jj = ii; jj != 3;++jj) {
      ts_2(ii, jj) = 1 + ii + 10 * jj;
      t_2(ii, jj) = 1 + ii + 10 * jj;
      t_2(jj, ii) = 1 + ii + 10 * jj;
    }

  Tensor2<double, 3, 3> t_2_2;
  t_2_2(k, l) = t_4(i, j, k, l) * ts_2(i, j);
  t_2_2(k, l) -= t_4(i, j, k, l) * t_2(i, j);
  for (int ii = 0; ii != 3; ++ii)
    for (int jj = 0; jj != 3; ++jj) {
      test_for_zero(t_2_2(ii, jj), "T4(i,j,k,l)*Ts2(i,j)(" + to_string(ii) +
                                       "," + to_string(jj) + ")");
    }


}