#include <iostream>
#include "../../src/FTensor.hpp"
#include "test_for_zero.hpp"
using namespace FTensor;
using namespace std;

void test_T4(
  const Tensor1<double,3> &t1_1,
  const Tensor2<double,3,3> &t2_1,
  const Tensor2_symmetric<double,3>  &t2s_1,
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

  #define TESTING_ASSIGMENT(I,J,K,L,II,JJ,KK,LL) \
  t4(i,j,k,l) = t4_1(I,J,K,L); \
  for(ii = 0;ii!=3;ii++) { \
    for(jj = 0;jj!=3;jj++) { \
      for(kk = 0;kk!=3;kk++) { \
        for(int ll = 0;ll!=3;ll++) { \
          test_for_zero(t4(ii,jj,kk,ll) - t4_1(II,JJ,KK,LL), \
          "T4_equals_T4 Assignment_" # I # J # K # L \
        ); \
        } \
      } \
    } \
  }

  int ii,jj,kk,ll;

  // jikl
  TESTING_ASSIGMENT(j,i,k,l, jj,ii,kk,ll);

  // jkil
  TESTING_ASSIGMENT(j,k,i,l, jj,kk,ii,ll);

  // jkli
  TESTING_ASSIGMENT(j,k,l,i, jj,kk,ll,ii);

  // kjli
  TESTING_ASSIGMENT(k,j,l,i, kk,jj,ll,ii);

  // klji
  TESTING_ASSIGMENT(k,l,j,i, kk,ll,jj,ii);

  // klij
  TESTING_ASSIGMENT(k,l,i,j, kk,ll,ii,jj);

  // lkij
  TESTING_ASSIGMENT(l,k,i,j, ll,kk,ii,jj);

  // likj
  TESTING_ASSIGMENT(l,i,k,j, ll,ii,kk,jj);

  // lijk
  TESTING_ASSIGMENT(l,i,j,k, ll,ii,jj,kk);

  // iljk
  TESTING_ASSIGMENT(i,l,j,k, ii,ll,jj,kk);

  // ijlk
  TESTING_ASSIGMENT(i,j,l,k, ii,jj,ll,kk);

  // lkji
  TESTING_ASSIGMENT(l,k,j,i, ll,kk,jj,ii);

  // ikjl
  TESTING_ASSIGMENT(i,k,j,l, ii,kk,jj,ll);

  // iklj
  TESTING_ASSIGMENT(i,k,l,j, ii,kk,ll,jj);

  // ilkj
  TESTING_ASSIGMENT(i,l,k,j, ii,ll,kk,jj);

  // jilk
  TESTING_ASSIGMENT(j,i,l,k, jj,ii,ll,kk);

  // kijl
  TESTING_ASSIGMENT(k,i,j,l, kk,ii,jj,ll);

  // kilj
  TESTING_ASSIGMENT(k,i,j,l, kk,ii,jj,ll);

  // jlik
  TESTING_ASSIGMENT(j,l,i,k, jj,ll,ii,kk);

  // jlki
  TESTING_ASSIGMENT(j,l,k,i, jj,ll,kk,ii);

  // kjil
  TESTING_ASSIGMENT(k,j,i,l, kk,jj,ii,ll);

  // ljik
  TESTING_ASSIGMENT(l,j,i,k, ll,jj,ii,kk);

  // ljki
  TESTING_ASSIGMENT(l,j,k,i, ll,jj,kk,ii);


  #undef TESTING_ASSIGMENT

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

  t4_31(i,j,k,l) = t3(i,j,k)*t1_1(l);
  FTensor::Tensor4<double,3,3,3,3> t4_13;
  t4_13(i,j,k,l) = t1_1(i)*t3(j,k,l);

  t4(i,j,k,l) = t4_1(i,j,k,l);
  t4(i,j,k,l) = 1;
  for(int ii = 0;ii!=3;ii++) {
    for(int jj = 0;jj!=3;jj++) {
      for(int kk = 0;kk!=3;kk++) {
        for(int ll = 0;ll!=3;ll++) {
          test_for_zero(t4(ii,jj,kk,ll) - 1,"T4_equals_generic");
        }
      }
    }
  }

  Tensor4_ddg<double,3,3> t4ddg_1;
  Tensor4<double,3,3,3,3> t4_222;
  Tensor2<double,3,3> t2_cpy;
  for(int ii = 0;ii!=3;ii++) {
    for(int jj = 0;jj!=3;jj++) {
      t2_cpy(ii,jj) = t2s_1(ii,jj);
    }
  }

  t4ddg_1(i,j,k,l)=t2s_1(i,j)*t2s_1(k,l);
  t4(i,j,k,m) = t4ddg_1(i,j,k,l)*t2_cpy(l,m);
  t2(k,m) = t2s_1(k,l)*t2s_1(l,m);
  t4_222(i,j,k,m) = t2_cpy(i,j)*t2(k,m);

  for(int ii = 0;ii!=3;ii++) {
    for(int jj = 0;jj!=3;jj++) {
      for(int kk = 0;kk!=3;kk++) {
        for(int ll = 0;ll!=3;ll++) {
          std::cerr << t4(ii,jj,kk,ll) << " " << t4_222(ii,jj,kk,ll) << std::endl;
          test_for_zero(t4(ii,jj,kk,ll) - t4_222(ii,jj,kk,ll),"Tensor4_ddg_times_Tensor2_3");
        }
      }
    }
  }


}
