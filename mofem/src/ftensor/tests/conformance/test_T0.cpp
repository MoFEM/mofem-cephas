#include <iostream>
#include "../../src/FTensor.hpp"
#include "test_for_zero.hpp"
using namespace FTensor;
using namespace std;

void test_T0(const int &T, Tensor0<double*> &t0_1,
	     const Tensor0<double*> &t0_2)
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

  /* Tensor0 test cases. */

  t0_1=T+t0_2;
  test_for_zero(t0_1-T-t0_2,"T+T0");
  t0_1+=10+t0_2;
  test_for_zero(t0_1-(10+T+2*t0_2),"T0+=T");
  t0_1-=5+3*t0_2;
  test_for_zero(t0_1-(5+T-t0_2),"T0-=T");
  t0_1*=2+t0_2;
  test_for_zero(t0_1-(5+T-t0_2)*(2+t0_2),"T0*=T");
  t0_1/=7.0+t0_2;
  test_for_zero(t0_1-(5+T-t0_2)*(2+t0_2)/(7.0+t0_2),"T0/=T");
}
