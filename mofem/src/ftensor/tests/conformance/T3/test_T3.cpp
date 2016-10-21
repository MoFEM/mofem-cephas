#include <iostream>
#include "../../../src/FTensor.hpp"
#include "../test_for_zero.hpp"
using namespace FTensor;
using namespace std;

void test_T3_01(const Tensor1<double,3> &t1_1, const Tensor1<double,3> &t1_2,
                const Tensor2<double,3,3> &t2_2, const Tensor2<double,3,3> &t2_3,
                const Tensor2_symmetric<double,3> &t2s_2,
                const Tensor2_symmetric<double,3> &t2s_3,
                const Tensor3_dg<double,3,3> &t3dg_2);

void test_T3_02(const Tensor1<double,3> &t1_1, const Tensor1<double,3> &t1_2,
                const Tensor2<double,3,3> &t2_2,
                const Tensor2<double,3,3> &t2_3);

void test_T3_03(const Tensor3<double,3,3,3> &t3_1);

void test_T3(const Tensor1<double,3> &t1_1, const Tensor1<double,3> &t1_2,
             const Tensor2<double,3,3> &t2_2,
             const Tensor2<double,3,3> &t2_3,
             const Tensor2_symmetric<double,3> &t2s_2,
             const Tensor2_symmetric<double,3> &t2s_3,
             const Tensor3<double,3,3,3> &t3_1,
             const Tensor3_dg<double,3,3> &t3dg_2)
{
  test_T3_01(t1_1,t1_2,t2_2,t2_3,t2s_2,t2s_3,t3dg_2);
  test_T3_02(t1_1,t1_2,t2_2,t2_3);
  test_T3_03(t3_1);

  Index<'i',3> i;
  Index<'j',3> j;
  Index<'k',3> k;

  Tensor3<double,3,3,3> t3_2;
  t3_2(i,j,k) = t2_2(i,j)*t1_1(k);
  Tensor2<double,3,3> t2_1,t2_4;
  t2_1(i,j) = t3_2(i,j,k)*t1_2(k);
  t2_4(i,j) = t2_2(i,j)*(t1_1(k)*t1_2(k));

  for(int ii = 0;ii!=3;ii++) {
    for(int jj = 0;jj!=3;jj++) {
      test_for_zero(t2_1(ii,jj)-t2_4(ii,jj),"Tensor3_times_Tensor1_2");
    }
  }

}
