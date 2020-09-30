#include "../../../src/FTensor.hpp"
#include "../test_for_zero.hpp"
#include <iostream>
using namespace FTensor;
using namespace std;

void test_T2s_47(const Tensor2_symmetric<double, 3> &t2s_1,
                 const Tensor2_symmetric<double, 3> &t2s_2) {

  Index<'i', 3> i;
  Index<'j', 3> j;
  Index<'k', 3> k;
  Index<'l', 3> l;

  Ddg<double, 3, 3> t4s;
  t4s(i, j, k, l) = t2s_1(i, k) ^ t2s_2(j, l);

  for (int ii = 0; ii != 3; ++ii)
    for (int jj = 0; jj != 3; ++jj)
      for (int kk = 0; kk != 3; ++kk)
        for (int ll = 0; ll != 3; ++ll) {

          auto small_eval = [&](auto n1, auto n2, auto n3, auto n4) {
            return t2s_1(n1, n3) * t2s_2(n2, n4);
          };

          const auto n1 = ii;
          const auto n2 = jj;
          const auto n3 = kk;
          const auto n4 = ll;

          const double t =
              small_eval(n1, n2, n3, n4) + small_eval(n2, n1, n3, n4) +
              small_eval(n1, n2, n4, n3) + small_eval(n2, n1, n4, n3);

          test_for_zero(t4s(ii, jj, kk, ll) - t, "t2s_1(i, k) ^ t2s_2(j, k)");

        }
}
