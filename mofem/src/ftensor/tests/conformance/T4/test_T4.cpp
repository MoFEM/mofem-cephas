#include <iostream>
#include "../../../src/FTensor.hpp"
#include "../test_for_zero.hpp"
using namespace FTensor;
using namespace std;

void test_T4_01(const Tensor4<double, 3, 3, 3, 3> &t4_1);

void test_T4_02(const Tensor4<double, 3, 3, 3, 3> &t4_1,
                const Tensor4<double, 3, 3, 3, 3> &t4_2);

void test_T4_03(const Tensor4<double, 3, 3, 3, 3> &t4_1,
                const Tensor4<double, 3, 3, 3, 3> &t4_2);

void test_T4_04(const Tensor4<double, 3, 3, 3, 3> &t4_1);
void test_T4_iostream();

void test_T4(const Tensor4<double, 3, 3, 3, 3> &t4_1,
             const Tensor4<double, 3, 3, 3, 3> &t4_2)
{
  test_T4_01(t4_1);
  test_T4_02(t4_1, t4_2);
  test_T4_03(t4_1, t4_2);
  test_T4_04(t4_1);
  test_T4_iostream();
}
