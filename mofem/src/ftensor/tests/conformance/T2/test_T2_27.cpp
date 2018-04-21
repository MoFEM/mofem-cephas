#include <iostream>
#include "../../../src/FTensor.hpp"
#include "../test_for_zero.hpp"
using namespace FTensor;
using namespace std;

void test_T2_27(const Tensor2<double, 3, 3> &t2_2,
                const Tensor2<double, 3, 3> &t2_3)
{
  Index<'i', 3> i;
  Index<'j', 3> j;
  Index<'k', 3> k;
  Index<'l', 3> l;
  Index<'m', 3> m;
  Index<'n', 3> n;

  Number<0> N0;
  Number<1> N1;
  Number<2> N2;

  /* Tensor2 tests */

  /* Tensor2*Tensor2 */
  Tensor2<double, 3, 3> t2;

  t2(i, k) = t2_2(i, j) * t2_3(k, j);
  test_for_zero(t2(0, 0)
                  - (t2_2(0, 0) * t2_3(0, 0) + t2_2(0, 1) * t2_3(0, 1)
                     + t2_2(0, 2) * t2_3(0, 2)),
                "T2(i,j)*T2(k,j)(0,0)");
  test_for_zero(t2(0, 1)
                  - (t2_2(0, 0) * t2_3(1, 0) + t2_2(0, 1) * t2_3(1, 1)
                     + t2_2(0, 2) * t2_3(1, 2)),
                "T2(i,j)*T2(k,j)(0,1)");
  test_for_zero(t2(0, 2)
                  - (t2_2(0, 0) * t2_3(2, 0) + t2_2(0, 1) * t2_3(2, 1)
                     + t2_2(0, 2) * t2_3(2, 2)),
                "T2(i,j)*T2(k,j)(0,2)");
  test_for_zero(t2(1, 0)
                  - (t2_2(1, 0) * t2_3(0, 0) + t2_2(1, 1) * t2_3(0, 1)
                     + t2_2(1, 2) * t2_3(0, 2)),
                "T2(i,j)*T2(k,j)(1,0)");
  test_for_zero(t2(1, 1)
                  - (t2_2(1, 0) * t2_3(1, 0) + t2_2(1, 1) * t2_3(1, 1)
                     + t2_2(1, 2) * t2_3(1, 2)),
                "T2(i,j)*T2(k,j)(1,1)");
  test_for_zero(t2(1, 2)
                  - (t2_2(1, 0) * t2_3(2, 0) + t2_2(1, 1) * t2_3(2, 1)
                     + t2_2(1, 2) * t2_3(2, 2)),
                "T2(i,j)*T2(k,j)(1,2)");
  test_for_zero(t2(2, 0)
                  - (t2_2(2, 0) * t2_3(0, 0) + t2_2(2, 1) * t2_3(0, 1)
                     + t2_2(2, 2) * t2_3(0, 2)),
                "T2(i,j)*T2(k,j)(2,0)");
  test_for_zero(t2(2, 1)
                  - (t2_2(2, 0) * t2_3(1, 0) + t2_2(2, 1) * t2_3(1, 1)
                     + t2_2(2, 2) * t2_3(1, 2)),
                "T2(i,j)*T2(k,j)(2,1)");
  test_for_zero(t2(2, 2)
                  - (t2_2(2, 0) * t2_3(2, 0) + t2_2(2, 1) * t2_3(2, 1)
                     + t2_2(2, 2) * t2_3(2, 2)),
                "T2(i,j)*T2(k,j)(2,2)");
}
