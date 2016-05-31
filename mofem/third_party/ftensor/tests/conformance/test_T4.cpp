#include <iostream>
#include "../../src/FTensor.hpp"
#include "test_for_zero.hpp"
using namespace FTensor;
using namespace std;

void test_T4(
  const Tensor1<double,3> &t1_1,
  const Tensor2<double,3,3> &t2_1,
  const Tensor4<double,3,3,3,3> &t4_1
) {
  Index<'i',3> i;
  Index<'j',3> j;
  Index<'k',3> k;
  Index<'l',3> l;
  Index<'m',3> m;
  Index<'n',3> n;

  Number<0> N0;
  Number<1> N1;
  Number<2> N2;

  // for(int ii = 0;ii!=3;ii++) {
  //   for(int jj = 0;jj!=3;jj++) {
  //     for(int kk = 0;kk!=3;kk++) {
  //       for(int ll = 0;ll!=3;ll++) {
  //         std::cout << t4_1(ii,jj,kk,ll) << endl;
  //       }
  //     }
  //   }
  // }

  Tensor4<double,3,3,3,3> t4;
  t4(i,j,k,l) = t4_1(i,j,k,l);
  for(int ii = 0;ii!=3;ii++) {
    for(int jj = 0;jj!=3;jj++) {
      for(int kk = 0;kk!=3;kk++) {
        for(int ll = 0;ll!=3;ll++) {
          double a = (ii+1)*1000+(jj+1)*100+(kk+1)*10+(ll+1)*1;
          // std::cout << t4(ii,jj,kk,ll) << " " << a << endl;
          test_for_zero(t4(ii,jj,kk,ll) - a,"T4_equals_T4");
        }
      }
    }
  }

  FTensor::Tensor2<double,3,3> t2;
  t2(i,j) = t1_1(i)*t1_1(j);
  FTensor::Tensor3<double,3,3,3> t3;
  t3(i,j,k) = t1_1(i)*t1_1(j)*t1_1(k);

  FTensor::Tensor4<double,3,3,3,3> t4_31;
  t4_31(i,j,k,l) = t3(i,j,k)*t1_1(l);
  FTensor::Tensor4<double,3,3,3,3> t4_22;
  t4_22(i,j,k,l) = t2(i,j)*t2(k,l);

  FTensor::Tensor4<double,3,3,3,3> t4_1111;
  t4_1111(i,j,k,l) = t1_1(i)*t1_1(j)*t1_1(k)*t1_1(l);

  // t4(i,j,k,l) = t2_1(i,j)*t2_1(k,l);
  for(int ii = 0;ii!=3;ii++) {
    for(int jj = 0;jj!=3;jj++) {
      for(int kk = 0;kk!=3;kk++) {
        for(int ll = 0;ll!=3;ll++) {
          // double a = (ii+1)*1000+(jj+1)*100+(kk+1)*10+(ll+1)*1;
          // std::cout << t4(ii,jj,kk,ll) << endl;
          test_for_zero(t4_22(ii,jj,kk,ll) - t4_31(ii,jj,kk,ll),"T2_times_T2 - T3_times_T1 -> T4");
          test_for_zero(t4_22(ii,jj,kk,ll) - t4_1111(ii,jj,kk,ll),"T2_times_T2 - T1_T1_T1_T1 -> T4");
        }
      }
    }
  }


}
