#include "../../../src/FTensor.hpp"
#include "../test_for_zero.hpp"
#include <iostream>
using namespace FTensor;
using namespace std;

void test_T3_004(const Tensor3<double, 3, 3, 3> &t3_1) {

  Index<'i', 3> i;
  Index<'j', 3> j;
  Index<'k', 3> k;

  {
    // Symmetrization of first two indices yields dg
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

  // Tensor 3 times symmetric tensor 2 yields tensor 3 
  {
    Tensor3<double, 1, 3, 2> t_3_1;
    Tensor2<double, 3, 3> t_2_1;
    Tensor2_symmetric<double, 3> t_2s_1;
    Tensor2<double, 3, 3> t_2_2;
    Tensor3<double, 1, 3, 2> t_3_2;
    for (int ii = 0; ii != 1; ++ii)
      for (int jj = 0; jj != 3; ++jj) {
        for (int kk = 0; kk != 2; ++kk) {
          t_3_1(ii, jj, kk) = ii + 10. * jj + 100. * kk;
        }
      }
    for (int jj = 0; jj != 3; ++jj)
      for (int ll = 0; ll != 3; ++ll) {
        t_2_1(jj, ll) = jj + 10. * ll;
      }
    Index<'i', 1> i;
    Index<'j', 3> j;
    Index<'k', 2> k;
    Index<'l', 3> l;
    t_2s_1(j, l) = t_2_1(j, l) || t_2_1(l, j);
    t_2_2(j, l) = t_2_1(j, l) + t_2_1(l, j);
    t_3_2(i, l, k) =
        t_3_1(i, j, k) * t_2s_1(j, l) -  t_2s_1(j, l)*t_3_1(i, j, k);
    for (int ii = 0; ii != 1; ++ii)
      for (int ll = 0; ll != 3; ++ll)
        for (int kk = 0; kk != 2; ++kk) {
          test_for_zero(t_3_2(ii, ll, kk),
                        "T3(i,j,k)||T2s(j,l)(" + to_string(ii) + "," +
                            to_string(ll) + "," + to_string(kk) + ")");
        }
    t_3_2(i, l, k) =
        t_3_1(i, j, k) * t_2s_1(j, l) - t_2_2(j, l) * t_3_1(i, j, k);
    for (int ii = 0; ii != 1; ++ii)
      for (int ll = 0; ll != 3; ++ll)
        for (int kk = 0; kk != 2; ++kk) {
          test_for_zero(t_3_2(ii, ll, kk),
                        "T3(i,j,k)||T2s(j,l)(" + to_string(ii) + "," +
                            to_string(ll) + "," + to_string(kk) + ")");
        }
 }
}