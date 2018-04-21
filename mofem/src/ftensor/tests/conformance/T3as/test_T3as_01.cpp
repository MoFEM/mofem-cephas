#include "../../../src/FTensor.hpp"
#include "../test_for_zero.hpp"
#include <iostream>
using namespace FTensor;
using namespace std;

void test_T3as_01(const Dg<double, 3, 3> &t3dg_2,
                  Tensor3_antisymmetric<double, 3, 3> &t3as_1)
{
  Index<'i', 3> i;
  Index<'j', 3> j;
  Index<'k', 3> k;

  Number<0> N0;
  Number<1> N1;
  Number<2> N2;

  t3as_1(i, j, k) = (t3dg_2(i, j, k) && t3dg_2(i, k, j));

  test_for_zero(t3as_1(0, 0, 0) - (t3dg_2(0, 0, 0) - t3dg_2(0, 0, 0)),
                "T3as(i,j,k)=T3dg(i,j,k) && T3dg(i,k,j)(0,0,0)");
  test_for_zero(t3as_1(0, 0, 1) - (t3dg_2(0, 0, 1) - t3dg_2(0, 1, 0)),
                "T3as(i,j,k)=T3dg(i,j,k) && T3dg(i,k,j)(0,0,1)");
  test_for_zero(t3as_1(0, 0, 2) - (t3dg_2(0, 0, 2) - t3dg_2(0, 2, 0)),
                "T3as(i,j,k)=T3dg(i,j,k) && T3dg(i,k,j)(0,0,2)");
  test_for_zero(t3as_1(0, 1, 0) - (t3dg_2(0, 1, 0) - t3dg_2(0, 0, 1)),
                "T3as(i,j,k)=T3dg(i,j,k) && T3dg(i,k,j)(0,1,0)");
  test_for_zero(t3as_1(0, 1, 1) - (t3dg_2(0, 1, 1) - t3dg_2(0, 1, 1)),
                "T3as(i,j,k)=T3dg(i,j,k) && T3dg(i,k,j)(0,1,1)");
  test_for_zero(t3as_1(0, 1, 2) - (t3dg_2(0, 1, 2) - t3dg_2(0, 2, 1)),
                "T3as(i,j,k)=T3dg(i,j,k) && T3dg(i,k,j)(0,1,2)");
  test_for_zero(t3as_1(0, 2, 0) - (t3dg_2(0, 2, 0) - t3dg_2(0, 0, 2)),
                "T3as(i,j,k)=T3dg(i,j,k) && T3dg(i,k,j)(0,2,0)");
  test_for_zero(t3as_1(0, 2, 1) - (t3dg_2(0, 2, 1) - t3dg_2(0, 1, 2)),
                "T3as(i,j,k)=T3dg(i,j,k) && T3dg(i,k,j)(0,2,1)");
  test_for_zero(t3as_1(0, 2, 2) - (t3dg_2(0, 2, 2) - t3dg_2(0, 2, 2)),
                "T3as(i,j,k)=T3dg(i,j,k) && T3dg(i,k,j)(0,2,2)");
  test_for_zero(t3as_1(1, 0, 0) - (t3dg_2(1, 0, 0) - t3dg_2(1, 0, 0)),
                "T3as(i,j,k)=T3dg(i,j,k) && T3dg(i,k,j)(1,0,0)");
  test_for_zero(t3as_1(1, 0, 1) - (t3dg_2(1, 0, 1) - t3dg_2(1, 1, 0)),
                "T3as(i,j,k)=T3dg(i,j,k) && T3dg(i,k,j)(1,0,1)");
  test_for_zero(t3as_1(1, 0, 2) - (t3dg_2(1, 0, 2) - t3dg_2(1, 2, 0)),
                "T3as(i,j,k)=T3dg(i,j,k) && T3dg(i,k,j)(1,0,2)");
  test_for_zero(t3as_1(1, 1, 0) - (t3dg_2(1, 1, 0) - t3dg_2(1, 0, 1)),
                "T3as(i,j,k)=T3dg(i,j,k) && T3dg(i,k,j)(1,1,0)");
  test_for_zero(t3as_1(1, 1, 1) - (t3dg_2(1, 1, 1) - t3dg_2(1, 1, 1)),
                "T3as(i,j,k)=T3dg(i,j,k) && T3dg(i,k,j)(1,1,1)");
  test_for_zero(t3as_1(1, 1, 2) - (t3dg_2(1, 1, 2) - t3dg_2(1, 2, 1)),
                "T3as(i,j,k)=T3dg(i,j,k) && T3dg(i,k,j)(1,1,2)");
  test_for_zero(t3as_1(1, 2, 0) - (t3dg_2(1, 2, 0) - t3dg_2(1, 0, 2)),
                "T3as(i,j,k)=T3dg(i,j,k) && T3dg(i,k,j)(1,2,0)");
  test_for_zero(t3as_1(1, 2, 1) - (t3dg_2(1, 2, 1) - t3dg_2(1, 1, 2)),
                "T3as(i,j,k)=T3dg(i,j,k) && T3dg(i,k,j)(1,2,1)");
  test_for_zero(t3as_1(1, 2, 2) - (t3dg_2(1, 2, 2) - t3dg_2(1, 2, 2)),
                "T3as(i,j,k)=T3dg(i,j,k) && T3dg(i,k,j)(1,2,2)");
  test_for_zero(t3as_1(2, 0, 0) - (t3dg_2(2, 0, 0) - t3dg_2(2, 0, 0)),
                "T3as(i,j,k)=T3dg(i,j,k) && T3dg(i,k,j)(2,0,0)");
  test_for_zero(t3as_1(2, 0, 1) - (t3dg_2(2, 0, 1) - t3dg_2(2, 1, 0)),
                "T3as(i,j,k)=T3dg(i,j,k) && T3dg(i,k,j)(2,0,1)");
  test_for_zero(t3as_1(2, 0, 2) - (t3dg_2(2, 0, 2) - t3dg_2(2, 2, 0)),
                "T3as(i,j,k)=T3dg(i,j,k) && T3dg(i,k,j)(2,0,2)");
  test_for_zero(t3as_1(2, 1, 0) - (t3dg_2(2, 1, 0) - t3dg_2(2, 0, 1)),
                "T3as(i,j,k)=T3dg(i,j,k) && T3dg(i,k,j)(2,1,0)");
  test_for_zero(t3as_1(2, 1, 1) - (t3dg_2(2, 1, 1) - t3dg_2(2, 1, 1)),
                "T3as(i,j,k)=T3dg(i,j,k) && T3dg(i,k,j)(2,1,1)");
  test_for_zero(t3as_1(2, 1, 2) - (t3dg_2(2, 1, 2) - t3dg_2(2, 2, 1)),
                "T3as(i,j,k)=T3dg(i,j,k) && T3dg(i,k,j)(2,1,2)");
  test_for_zero(t3as_1(2, 2, 0) - (t3dg_2(2, 2, 0) - t3dg_2(2, 0, 2)),
                "T3as(i,j,k)=T3dg(i,j,k) && T3dg(i,k,j)(2,2,0)");
  test_for_zero(t3as_1(2, 2, 1) - (t3dg_2(2, 2, 1) - t3dg_2(2, 1, 2)),
                "T3as(i,j,k)=T3dg(i,j,k) && T3dg(i,k,j)(2,2,1)");
  test_for_zero(t3as_1(2, 2, 2) - (t3dg_2(2, 2, 2) - t3dg_2(2, 2, 2)),
                "T3as(i,j,k)=T3dg(i,j,k) && T3dg(i,k,j)(2,2,2)");
}
