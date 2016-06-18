#include <iostream>
#include "../../../src/FTensor.hpp"
#include "../test_for_zero.hpp"
using namespace FTensor;
using namespace std;

void test_T4ddg_01(const Tensor2_symmetric<double,3> &t2s_2,
                   const Tensor2_symmetric<double,3> &t2s_3,
                   const Tensor3_dg<double,3,3> &t3dg_2,
                   const Tensor3_dg<double,3,3> &t3dg_3);
void test_T4ddg_02(const Tensor2_symmetric<double,3> &t2s_2);
void test_T4ddg_03(Tensor2<double,3,3> &t2_1,
                   Tensor2_symmetric<double,3> &t2s_1,
                   const Tensor2_symmetric<double,3> &t2s_2,
                   const Tensor2_symmetric<double,3> &t2s_3);
void test_T4ddg_04(const Tensor1<double,3> &t1_2,
                   const Tensor2<double,3,3> &t2_2,
                   const Tensor2<double,3,3> &t2_3,
                   Tensor2_symmetric<double,3> &t2s_1,
                   const Tensor2_symmetric<double,3> &t2s_2,
                   const Tensor2_symmetric<double,3> &t2s_3,
                   Tensor3_dg<double,3,3> &t3dg_1);
void test_T4ddg_05(Tensor1<double,3> &t1_1, const Tensor1<double,3> &t1_2,
                   Tensor2<double,3,3> &t2_1,
                   Tensor2_symmetric<double,3> &t2s_1,
                   const Tensor2_symmetric<double,3> &t2s_2,
                   const Tensor2_symmetric<double,3> &t2s_3,
                   Tensor3_dg<double,3,3> &t3dg_1,
                   const Tensor3_dg<double,3,3> &t3dg_2,
                   const Tensor3_dg<double,3,3> &t3dg_3);


void test_T4ddg(Tensor1<double,3> &t1_1, const Tensor1<double,3> &t1_2,
                Tensor2<double,3,3> &t2_1, const Tensor2<double,3,3> &t2_2,
                const Tensor2<double,3,3> &t2_3,
                Tensor2_symmetric<double,3> &t2s_1,
                const Tensor2_symmetric<double,3> &t2s_2,
                const Tensor2_symmetric<double,3> &t2s_3,
                Tensor3_dg<double,3,3> &t3dg_1,
                const Tensor3_dg<double,3,3> &t3dg_2,
                const Tensor3_dg<double,3,3> &t3dg_3)
{
 test_T4ddg_01(t2s_2,t2s_3,t3dg_2,t3dg_3);
 test_T4ddg_02(t2s_2);
 test_T4ddg_03(t2_1,t2s_1,t2s_2,t2s_3);
 test_T4ddg_04(t1_2,t2_2,t2_3,t2s_1,t2s_2,t2s_3,t3dg_1);
 test_T4ddg_05(t1_1,t1_2,t2_1,t2s_1,t2s_2,t2s_3,t3dg_1,t3dg_2,t3dg_3);

 /* Generic assignment and operators that works with doubles, ints etc. */

 Index<'i',3> i;
 Index<'j',3> j;
 Index<'k',3> k;
 Index<'l',3> l;

 Tensor4_ddg<double,3,3> t4ddg_1;
 t4ddg_1(i,j,k,l)=1;

 for(int ii = 0;ii!=3;ii++) {
   for(int jj = 0;jj!=3;jj++) {
     for(int kk = 0;kk!=3;kk++) {
       for(int ll = 0;ll!=3;ll++) {
         // std::cerr << t4(ii,jj,kk,ll) << " " << t4_222(ii,jj,kk,ll) << std::endl;
         test_for_zero(t4ddg_1(ii,jj,kk,ll) - 1,"T4_ddg_equals_generic");
       }
     }
   }
 }

 Tensor4_ddg<double,3,3> t4ddg_2;

 t4ddg_2(i,k,j,l)=t2s_2(i,k)*t2s_2(j,l);
 t4ddg_1(i,j,k,l)=t4ddg_2(i,j,k,l);
 t4ddg_2(i,j,k,l)*=2;

 for(int ii = 0;ii!=3;ii++) {
   for(int jj = 0;jj!=3;jj++) {
     for(int kk = 0;kk!=3;kk++) {
       for(int ll = 0;ll!=3;ll++) {
         std::cout << t4ddg_2(ii,jj,kk,ll) << " " << t4ddg_1(ii,jj,kk,ll) << std::endl;
         test_for_zero(t4ddg_2(ii,jj,kk,ll) - 2*t4ddg_1(ii,jj,kk,ll),"T4_ddg_equals_generic");
       }
     }
   }
 }

}
