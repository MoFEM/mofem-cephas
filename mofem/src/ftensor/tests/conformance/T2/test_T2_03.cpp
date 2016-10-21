#include <iostream>
#include "../../../src/FTensor.hpp"
#include "../test_for_zero.hpp"
using namespace FTensor;
using namespace std;

void test_T2_03()
{
  Index<'i',3> i;
  Index<'j',3> j;
  Index<'k',3> k;
  Index<'l',3> l;
  Index<'m',3> m;
  Index<'n',3> n;
  
  Number<0> N0;
  Number<1> N1;
  Number<2> N2;

  /* Tensor2 tests */
  Tensor2<double,3,3> t2_1;

  t2_1(N0,i)=10;
  t2_1(N1,i)=11;
  t2_1(N2,i)=12;
  test_for_zero(t2_1(0,0)-10,"T2(N,i)=T(0,0)");
  test_for_zero(t2_1(0,1)-10,"T2(N,i)=T(0,1)");
  test_for_zero(t2_1(0,2)-10,"T2(N,i)=T(0,2)");
  test_for_zero(t2_1(1,0)-11,"T2(N,i)=T(1,0)");
  test_for_zero(t2_1(1,1)-11,"T2(N,i)=T(1,1)");
  test_for_zero(t2_1(1,2)-11,"T2(N,i)=T(1,2)");
  test_for_zero(t2_1(2,0)-12,"T2(N,i)=T(2,0)");
  test_for_zero(t2_1(2,1)-12,"T2(N,i)=T(2,1)");
  test_for_zero(t2_1(2,2)-12,"T2(N,i)=T(2,2)");
}
