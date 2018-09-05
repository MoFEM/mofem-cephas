#include "../../../src/FTensor.hpp"
#include "../test_for_zero.hpp"
#include <iostream>
using namespace FTensor;
using namespace std;

void test_T3dg_038(const Dg<double, 3, 3> &t3dg_2,
                  const Dg<double, 3, 3> &t3dg_3) {
  Index<'i', 3> i;
  Index<'j', 3> j;
  Index<'k', 3> k;

  {
    Dg<double, 3, 3> t3dg_3_1;
    t3dg_3_1(i, j, k) = t3dg_3(i, j, k);
    t3dg_3_1(i, j, k) += t3dg_2(i, j, k);
    for (int ii = 0; ii != 3;++ii)
      for (int jj = 0; jj != 3;++jj)
        for (int kk = 0; kk != 3;++kk) {
          test_for_zero(t3dg_3_1(ii, jj, kk) - t3dg_2(ii, jj, kk) -
                            t3dg_3(ii, jj, kk),
                        "T3(i,j,k)+=T3(i,j,k)(" + to_string(ii) + "," +
                            to_string(jj) + "," + to_string(kk) + ")");
        }
  }
}