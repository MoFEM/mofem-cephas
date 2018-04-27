#include "../../../src/FTensor.hpp"
#include "../test_for_zero.hpp"
#include <iostream>
using namespace FTensor;
using namespace std;

void test_T2_02(const Tensor1<double, 3> &t1_1)
{
  Index<'i', 3> i;

  Number<0> N0;
  Number<1> N1;
  Number<2> N2;

  /* Tensor2 tests */
  Tensor2<double, 3, 3> t2_1;

  t2_1(N0, i) = t1_1(i);
  t2_1(N1, i) = t1_1(i);
  t2_1(N2, i) = t1_1(i);
  test_for_zero(t2_1(0, 0) - t1_1(0), "T2(N,i)=T1(i)(0,0)");
  test_for_zero(t2_1(0, 1) - t1_1(1), "T2(N,i)=T1(i)(0,1)");
  test_for_zero(t2_1(0, 2) - t1_1(2), "T2(N,i)=T1(i)(0,2)");
  test_for_zero(t2_1(1, 0) - t1_1(0), "T2(N,i)=T1(i)(1,0)");
  test_for_zero(t2_1(1, 1) - t1_1(1), "T2(N,i)=T1(i)(1,1)");
  test_for_zero(t2_1(1, 2) - t1_1(2), "T2(N,i)=T1(i)(1,2)");
  test_for_zero(t2_1(2, 0) - t1_1(0), "T2(N,i)=T1(i)(2,0)");
  test_for_zero(t2_1(2, 1) - t1_1(1), "T2(N,i)=T1(i)(2,1)");
  test_for_zero(t2_1(2, 2) - t1_1(2), "T2(N,i)=T1(i)(2,2)");
}
