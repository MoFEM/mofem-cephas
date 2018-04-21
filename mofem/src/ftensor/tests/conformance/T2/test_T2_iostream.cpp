#include "../../../src/FTensor.hpp"
#include "../test_ostream.hpp"
#include "../test_for_zero.hpp"

using namespace FTensor;
using namespace std;
void test_T2_iostream()
{
  test_ostream(
    Tensor2<double, 3, 3, column_major>(3, 7, 11, 13, 17, 23, 27, 31, 37),
    "[[3,7,11],[13,17,23],[27,31,37]]", "operator<<(T2<3,3,c>)");
  test_ostream(
    Tensor2<double, 3, 3, row_major>(3, 7, 11, 13, 17, 23, 27, 31, 37),
    "[[3,7,11],[13,17,23],[27,31,37]]", "operator<<(T2<3,3,r>)");
  test_ostream(Tensor2<double, 1, 2, column_major>(13, 17), "[[13,17]]",
               "operator<<(T2<1,2,c>)");
  test_ostream(Tensor2<double, 1, 2, row_major>(13, 17), "[[13,17]]",
               "operator<<(T2<1,2,r>)");
  test_ostream(Tensor2<double, 2, 1, column_major>(13, 17), "[[13],[17]]",
               "operator<<(T2<2,1,c>)");
  test_ostream(Tensor2<double, 2, 1, row_major>(13, 17), "[[13],[17]]",
               "operator<<(T2<2,1,r>)");

  Tensor2<double, 3, 3> t2_1;
  std::stringstream ss("[[3,7,11],[13,17,23],[27,31,37]]");
  ss >> t2_1;
  test_for_zero(t2_1(0, 0) - 3, "operator>>(T2)(0,0)");
  test_for_zero(t2_1(0, 1) - 7, "operator>>(T2)(0,1)");
  test_for_zero(t2_1(0, 2) - 11, "operator>>(T2)(0,2)");
  test_for_zero(t2_1(1, 0) - 13, "operator>>(T2)(1,0)");
  test_for_zero(t2_1(1, 1) - 17, "operator>>(T2)(1,1)");
  test_for_zero(t2_1(1, 2) - 23, "operator>>(T2)(1,2)");
  test_for_zero(t2_1(2, 0) - 27, "operator>>(T2)(2,0)");
  test_for_zero(t2_1(2, 1) - 31, "operator>>(T2)(2,1)");
  test_for_zero(t2_1(2, 2) - 37, "operator>>(T2)(2,2)");
}
