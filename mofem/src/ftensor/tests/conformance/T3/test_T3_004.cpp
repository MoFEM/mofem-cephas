#include "../../../src/FTensor.hpp"
#include "../test_for_zero.hpp"
#include <iostream>
using namespace FTensor;
using namespace std;

void test_T3_004(const Tensor3<double, 3, 3, 3> &t3_1) {

  Index<'i', 3> i;
  Index<'j', 3> j;
  Index<'k', 3> k;

  Dg<double, 3, 3> t_dg;
  t_dg(i, j, k) = t3_1(i, j, k) || t3_1(j, i, k);
  FTensor::Tensor3<double, 3, 3, 3> t_sym_3;
  t_sym_3(i, j, k) = t3_1(i, j, k) + t3_1(j, i, k);

  for (int ii = 0; ii != 3; ++ii)
    for (int jj = 0; jj != 3; ++jj)
      for (int kk = 0; kk != 3; ++kk) {
        test_for_zero(t_dg(ii, jj, kk) - t_sym_3(ii, jj, kk),
                      "T3(i,j,k)||T3(j,i,k)(" + to_string(ii) + "," +
                          to_string(jj) + "," + to_string(kk) + ")");
      }
}