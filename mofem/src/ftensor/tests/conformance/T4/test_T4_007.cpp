#include "../../../src/FTensor.hpp"
#include "../test_for_zero.hpp"
#include <iostream>
using namespace FTensor;
using namespace std;

void test_T4_007(const Tensor4<double, 1, 2, 3, 4> &t4,
                 const Tensor3<double, 2, 3, 4> &t3_2) {
  Index<'i', 1> i;
  Index<'j', 2> j;
  Index<'k', 3> k;
  Index<'l', 4> l;
  Index<'n', 3> n;

  Number<0> N0;
  Number<1> N1;
  Number<2> N2;

  Tensor3<double, 2, 3, 4> t32;
  t32(j, k, l) = 0.;

  // Yield tensor 4 setting 4th slot
  {
    Tensor4<double, 2, 3, 4, 3> t_4;
    t_4(j, k, l, N0) = t3_2(j, k, l);
    t_4(j, k, l, N1) = t32(j, k, l);
    t_4(j, k, l, N2) = t32(j, k, l);
    for (int ii = 0; ii != 2; ++ii)
      for (int jj = 0; jj != 3; ++jj)
        for (int kk = 0; kk != 4; ++kk) {
          test_for_zero(t_4(ii, jj, kk, 0) - t3_2(ii, jj, kk),
                        "T4(i,j,k,N0)(" + to_string(ii) + "," + to_string(jj) +
                            "," + to_string(kk) + ")");
          for (int ll : {1, 2}) {
            test_for_zero(t_4(ii, jj, kk, ll), "T4(i,j,k,ll)(" + to_string(ii) +
                                                   "," + to_string(jj) + "," +
                                                   to_string(kk) + "," +
                                                   to_string(ll) + ")");
          }
        }
  }
  {
    Tensor4<double, 2, 3, 4, 3> t_4;
    t_4(j, k, l, N0) = t32(j, k, l);
    t_4(j, k, l, N1) = t3_2(j, k, l);
    t_4(j, k, l, N2) = t32(j, k, l);
    for (int ii = 0; ii != 2; ++ii)
      for (int jj = 0; jj != 3; ++jj)
        for (int kk = 0; kk != 4; ++kk) {
          test_for_zero(t_4(ii, jj, kk, 1) - t3_2(ii, jj, kk),
                        "T4(i,j,k,N1)(" + to_string(ii) + "," + to_string(jj) +
                            "," + to_string(kk) + ")");
          for (int ll : {0, 2}) {
            test_for_zero(t_4(ii, jj, kk, ll), "T4(i,j,k,ll)(" + to_string(ii) +
                                                   "," + to_string(jj) + "," +
                                                   to_string(kk) + "," +
                                                   to_string(ll) + ")");
          }
        }
  }
  {
    Tensor4<double, 2, 3, 4, 3> t_4;
    t_4(j, k, l, N0) = t32(j, k, l);
    t_4(j, k, l, N1) = t32(j, k, l);
    t_4(j, k, l, N2) = t3_2(j, k, l);
    for (int ii = 0; ii != 2; ++ii)
      for (int jj = 0; jj != 3; ++jj)
        for (int kk = 0; kk != 4; ++kk) {
          test_for_zero(t_4(ii, jj, kk, 2) - t3_2(ii, jj, kk),
                        "T4(i,j,k,N2)(" + to_string(ii) + "," + to_string(jj) +
                            "," + to_string(kk) + ")");
          for (int ll : {0, 1}) {
            test_for_zero(t_4(ii, jj, kk, ll), "T4(i,j,k,ll)(" + to_string(ii) +
                                                   "," + to_string(jj) + "," +
                                                   to_string(kk) + "," +
                                                   to_string(ll) + ")");
          }
        }
  }

  // Yield tensor 4 setting 3th slot
  {
    Tensor4<double, 2, 3, 3, 4> t_4;
    t_4(j, k, N0, l) = t3_2(j, k, l);
    t_4(j, k, N1, l) = t32(j, k, l);
    t_4(j, k, N2, l) = t32(j, k, l);
    for (int ii = 0; ii != 2; ++ii)
      for (int jj = 0; jj != 3; ++jj)
        for (int kk = 0; kk != 4; ++kk) {
          test_for_zero(t_4(ii, jj, 0, kk) - t3_2(ii, jj, kk),
                        "T4(i,j,N0,k)(" + to_string(ii) + "," + to_string(jj) +
                            "," + to_string(kk) + ")");
          for (int ll : {1, 2}) {
            test_for_zero(t_4(ii, jj, ll, kk), "T4(i,j,l, k)(" + to_string(ii) +
                                                   "," + to_string(jj) + "," +
                                                   to_string(ll) + "," +
                                                   to_string(kk) + ")");
          }
        }
  }
  {
    Tensor4<double, 2, 3, 3, 4> t_4;
    t_4(j, k, N0, l) = t32(j, k, l);
    t_4(j, k, N1, l) = t3_2(j, k, l);
    t_4(j, k, N2, l) = t32(j, k, l);
    for (int ii = 0; ii != 2; ++ii)
      for (int jj = 0; jj != 3; ++jj)
        for (int kk = 0; kk != 4; ++kk) {
          test_for_zero(t_4(ii, jj, 1, kk) - t3_2(ii, jj, kk),
                        "T4(i,j,N1,k)(" + to_string(ii) + "," + to_string(jj) +
                            "," + to_string(kk) + ")");
          for (int ll : {0, 2}) {
            test_for_zero(t_4(ii, jj, ll, kk), "T4(i,j,k,ll)(" + to_string(ii) +
                                                   "," + to_string(jj) + "," +
                                                   to_string(ll) + "," +
                                                   to_string(kk) + ")");
          }
        }
  }
  {
    Tensor4<double, 2, 3, 3, 4> t_4;
    t_4(j, k, N0, l) = t32(j, k, l);
    t_4(j, k, N1, l) = t32(j, k, l);
    t_4(j, k, N2, l) = t3_2(j, k, l);
    for (int ii = 0; ii != 2; ++ii)
      for (int jj = 0; jj != 3; ++jj)
        for (int kk = 0; kk != 4; ++kk) {
          test_for_zero(t_4(ii, jj, 2, kk) - t3_2(ii, jj, kk),
                        "T4(i,j,N2,k)(" + to_string(ii) + "," + to_string(jj) +
                            "," + to_string(kk) + ")");
          for (int ll : {0, 1}) {
            test_for_zero(t_4(ii, jj, ll, kk), "T4(i,j,k,ll)(" + to_string(ii) +
                                                   "," + to_string(jj) + "," +
                                                   to_string(ll) + "," +
                                                   to_string(kk) + ")");
          }
        }
  }
}