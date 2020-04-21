#include "../../../src/FTensor.hpp"
#include "../test_for_zero.hpp"
#include <iostream>
using namespace FTensor;
using namespace std;

void test_Kornecker_Delta_01() {

  Index<'i', 2> i;
  Index<'j', 2> j;

  Index<'k', 3> k;
  Index<'l', 3> l;

  FTensor::Tensor2<int, 2, 2> t_one2{1, 0, 0, 1};
  FTensor::Tensor2<int, 3, 3> t_one3{1, 0, 0, 0, 1, 0, 0, 0, 1};

  FTensor::Tensor2_symmetric<int, 2> t_one2_symmetric;
  t_one2_symmetric(i, j) = t_one2(i, j) || t_one2(j, i);
  FTensor::Tensor2_symmetric<int, 3> t_one3_symmetric;
  t_one3_symmetric(k, l) = t_one3(k, l) || t_one3(l, k);

  t_one2(i, j) -= kronecker_delta(i, j);
  for (auto ii : {0, 1})
    for (auto jj : {0, 1})
      test_for_zero(t_one2(ii, jj), "kronecker_delta 2 by 2");

  t_one3(k, l) -= kronecker_delta(k, l);
  for (auto ii : {0, 1, 2})
    for (auto jj : {0, 1, 2})
      test_for_zero(t_one3(ii, jj), "kronecker_delta 3 by 3");

  t_one2_symmetric(i, j) -= kronecker_delta_symmetric(i, j);
  t_one2_symmetric(i, j) -= kronecker_delta_symmetric(j, i);
  for (auto ii : {0, 1})
    for (auto jj : {0, 1})
      test_for_zero(t_one2_symmetric(ii, jj),
                    "kronecker_delta_symmetric 2 by 2");

  t_one3_symmetric(k, l) -= kronecker_delta_symmetric(k, l);
  t_one3_symmetric(k, l) -= kronecker_delta_symmetric(l, k);
  for (auto ii : {0, 1, 2})
    for (auto jj : {0, 1, 2})
      test_for_zero(t_one3_symmetric(ii, jj),
                    "kronecker_delta_symmetric 3 by 3");
}