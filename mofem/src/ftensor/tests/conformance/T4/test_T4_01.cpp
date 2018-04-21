#include "../../../src/FTensor.hpp"
#include "../test_ostream.hpp"
#include "../test_for_zero.hpp"

using namespace FTensor;
using namespace std;

void test_T4_01(const Tensor4<double, 3, 3, 3, 3> &t4_1)
{
  Index<'i', 3> i;
  Index<'j', 3> j;
  Index<'k', 3> k;
  Index<'l', 3> l;
  Index<'m', 3> m;
  Index<'n', 3> n;

  Number<0> N0;
  Number<1> N1;
  Number<2> N2;

  Tensor4<double, 3, 3, 3, 3> t4;

  /*----------------------------------------------------------------------------*/
  t4(i, j, k, l) = t4_1(i, j, k, l);
  /*----------------------------------------------------------------------------*/

  test_for_zero(t4(0, 0, 0, 0) - t4_1(0, 0, 0, 0),
                "T4(i,j,k,l)=T4(i,j,k,l)(0,0,0,0)");
  test_for_zero(t4(0, 0, 0, 1) - t4_1(0, 0, 0, 1),
                "T4(i,j,k,l)=T4(i,j,k,l)(0,0,0,1)");
  test_for_zero(t4(0, 0, 0, 2) - t4_1(0, 0, 0, 2),
                "T4(i,j,k,l)=T4(i,j,k,l)(0,0,0,2)");
  test_for_zero(t4(0, 0, 1, 0) - t4_1(0, 0, 1, 0),
                "T4(i,j,k,l)=T4(i,j,k,l)(0,0,1,0)");
  test_for_zero(t4(0, 0, 1, 1) - t4_1(0, 0, 1, 1),
                "T4(i,j,k,l)=T4(i,j,k,l)(0,0,1,1)");
  test_for_zero(t4(0, 0, 1, 2) - t4_1(0, 0, 1, 2),
                "T4(i,j,k,l)=T4(i,j,k,l)(0,0,1,2)");
  test_for_zero(t4(0, 0, 2, 0) - t4_1(0, 0, 2, 0),
                "T4(i,j,k,l)=T4(i,j,k,l)(0,0,2,0)");
  test_for_zero(t4(0, 0, 2, 1) - t4_1(0, 0, 2, 1),
                "T4(i,j,k,l)=T4(i,j,k,l)(0,0,2,1)");
  test_for_zero(t4(0, 0, 2, 2) - t4_1(0, 0, 2, 2),
                "T4(i,j,k,l)=T4(i,j,k,l)(0,0,2,2)");
  test_for_zero(t4(0, 1, 0, 0) - t4_1(0, 1, 0, 0),
                "T4(i,j,k,l)=T4(i,j,k,l)(0,1,0,0)");
  test_for_zero(t4(0, 1, 0, 1) - t4_1(0, 1, 0, 1),
                "T4(i,j,k,l)=T4(i,j,k,l)(0,1,0,1)");
  test_for_zero(t4(0, 1, 0, 2) - t4_1(0, 1, 0, 2),
                "T4(i,j,k,l)=T4(i,j,k,l)(0,1,0,2)");
  test_for_zero(t4(0, 1, 1, 0) - t4_1(0, 1, 1, 0),
                "T4(i,j,k,l)=T4(i,j,k,l)(0,1,1,0)");
  test_for_zero(t4(0, 1, 1, 1) - t4_1(0, 1, 1, 1),
                "T4(i,j,k,l)=T4(i,j,k,l)(0,1,1,1)");
  test_for_zero(t4(0, 1, 1, 2) - t4_1(0, 1, 1, 2),
                "T4(i,j,k,l)=T4(i,j,k,l)(0,1,1,2)");
  test_for_zero(t4(0, 1, 2, 0) - t4_1(0, 1, 2, 0),
                "T4(i,j,k,l)=T4(i,j,k,l)(0,1,2,0)");
  test_for_zero(t4(0, 1, 2, 1) - t4_1(0, 1, 2, 1),
                "T4(i,j,k,l)=T4(i,j,k,l)(0,1,2,1)");
  test_for_zero(t4(0, 1, 2, 2) - t4_1(0, 1, 2, 2),
                "T4(i,j,k,l)=T4(i,j,k,l)(0,1,2,2)");
  test_for_zero(t4(0, 2, 0, 0) - t4_1(0, 2, 0, 0),
                "T4(i,j,k,l)=T4(i,j,k,l)(0,2,0,0)");
  test_for_zero(t4(0, 2, 0, 1) - t4_1(0, 2, 0, 1),
                "T4(i,j,k,l)=T4(i,j,k,l)(0,2,0,1)");
  test_for_zero(t4(0, 2, 0, 2) - t4_1(0, 2, 0, 2),
                "T4(i,j,k,l)=T4(i,j,k,l)(0,2,0,2)");
  test_for_zero(t4(0, 2, 1, 0) - t4_1(0, 2, 1, 0),
                "T4(i,j,k,l)=T4(i,j,k,l)(0,2,1,0)");
  test_for_zero(t4(0, 2, 1, 1) - t4_1(0, 2, 1, 1),
                "T4(i,j,k,l)=T4(i,j,k,l)(0,2,1,1)");
  test_for_zero(t4(0, 2, 1, 2) - t4_1(0, 2, 1, 2),
                "T4(i,j,k,l)=T4(i,j,k,l)(0,2,1,2)");
  test_for_zero(t4(0, 2, 2, 0) - t4_1(0, 2, 2, 0),
                "T4(i,j,k,l)=T4(i,j,k,l)(0,2,2,0)");
  test_for_zero(t4(0, 2, 2, 1) - t4_1(0, 2, 2, 1),
                "T4(i,j,k,l)=T4(i,j,k,l)(0,2,2,1)");
  test_for_zero(t4(0, 2, 2, 2) - t4_1(0, 2, 2, 2),
                "T4(i,j,k,l)=T4(i,j,k,l)(0,2,2,2)");

  test_for_zero(t4(1, 0, 0, 0) - t4_1(1, 0, 0, 0),
                "T4(i,j,k,l)=T4(i,j,k,l)(1,0,0,0)");
  test_for_zero(t4(1, 0, 0, 1) - t4_1(1, 0, 0, 1),
                "T4(i,j,k,l)=T4(i,j,k,l)(1,0,0,1)");
  test_for_zero(t4(1, 0, 0, 2) - t4_1(1, 0, 0, 2),
                "T4(i,j,k,l)=T4(i,j,k,l)(1,0,0,2)");
  test_for_zero(t4(1, 0, 1, 0) - t4_1(1, 0, 1, 0),
                "T4(i,j,k,l)=T4(i,j,k,l)(1,0,1,0)");
  test_for_zero(t4(1, 0, 1, 1) - t4_1(1, 0, 1, 1),
                "T4(i,j,k,l)=T4(i,j,k,l)(1,0,1,1)");
  test_for_zero(t4(1, 0, 1, 2) - t4_1(1, 0, 1, 2),
                "T4(i,j,k,l)=T4(i,j,k,l)(1,0,1,2)");
  test_for_zero(t4(1, 0, 2, 0) - t4_1(1, 0, 2, 0),
                "T4(i,j,k,l)=T4(i,j,k,l)(1,0,2,0)");
  test_for_zero(t4(1, 0, 2, 1) - t4_1(1, 0, 2, 1),
                "T4(i,j,k,l)=T4(i,j,k,l)(1,0,2,1)");
  test_for_zero(t4(1, 0, 2, 2) - t4_1(1, 0, 2, 2),
                "T4(i,j,k,l)=T4(i,j,k,l)(1,0,2,2)");
  test_for_zero(t4(1, 1, 0, 0) - t4_1(1, 1, 0, 0),
                "T4(i,j,k,l)=T4(i,j,k,l)(1,1,0,0)");
  test_for_zero(t4(1, 1, 0, 1) - t4_1(1, 1, 0, 1),
                "T4(i,j,k,l)=T4(i,j,k,l)(1,1,0,1)");
  test_for_zero(t4(1, 1, 0, 2) - t4_1(1, 1, 0, 2),
                "T4(i,j,k,l)=T4(i,j,k,l)(1,1,0,2)");
  test_for_zero(t4(1, 1, 1, 0) - t4_1(1, 1, 1, 0),
                "T4(i,j,k,l)=T4(i,j,k,l)(1,1,1,0)");
  test_for_zero(t4(1, 1, 1, 1) - t4_1(1, 1, 1, 1),
                "T4(i,j,k,l)=T4(i,j,k,l)(1,1,1,1)");
  test_for_zero(t4(1, 1, 1, 2) - t4_1(1, 1, 1, 2),
                "T4(i,j,k,l)=T4(i,j,k,l)(1,1,1,2)");
  test_for_zero(t4(1, 1, 2, 0) - t4_1(1, 1, 2, 0),
                "T4(i,j,k,l)=T4(i,j,k,l)(1,1,2,0)");
  test_for_zero(t4(1, 1, 2, 1) - t4_1(1, 1, 2, 1),
                "T4(i,j,k,l)=T4(i,j,k,l)(1,1,2,1)");
  test_for_zero(t4(1, 1, 2, 2) - t4_1(1, 1, 2, 2),
                "T4(i,j,k,l)=T4(i,j,k,l)(1,1,2,2)");
  test_for_zero(t4(1, 2, 0, 0) - t4_1(1, 2, 0, 0),
                "T4(i,j,k,l)=T4(i,j,k,l)(1,2,0,0)");
  test_for_zero(t4(1, 2, 0, 1) - t4_1(1, 2, 0, 1),
                "T4(i,j,k,l)=T4(i,j,k,l)(1,2,0,1)");
  test_for_zero(t4(1, 2, 0, 2) - t4_1(1, 2, 0, 2),
                "T4(i,j,k,l)=T4(i,j,k,l)(1,2,0,2)");
  test_for_zero(t4(1, 2, 1, 0) - t4_1(1, 2, 1, 0),
                "T4(i,j,k,l)=T4(i,j,k,l)(1,2,1,0)");
  test_for_zero(t4(1, 2, 1, 1) - t4_1(1, 2, 1, 1),
                "T4(i,j,k,l)=T4(i,j,k,l)(1,2,1,1)");
  test_for_zero(t4(1, 2, 1, 2) - t4_1(1, 2, 1, 2),
                "T4(i,j,k,l)=T4(i,j,k,l)(1,2,1,2)");
  test_for_zero(t4(1, 2, 2, 0) - t4_1(1, 2, 2, 0),
                "T4(i,j,k,l)=T4(i,j,k,l)(1,2,2,0)");
  test_for_zero(t4(1, 2, 2, 1) - t4_1(1, 2, 2, 1),
                "T4(i,j,k,l)=T4(i,j,k,l)(1,2,2,1)");
  test_for_zero(t4(1, 2, 2, 2) - t4_1(1, 2, 2, 2),
                "T4(i,j,k,l)=T4(i,j,k,l)(1,2,2,2)");

  test_for_zero(t4(2, 0, 0, 0) - t4_1(2, 0, 0, 0),
                "T4(i,j,k,l)=T4(i,j,k,l)(2,0,0,0)");
  test_for_zero(t4(2, 0, 0, 1) - t4_1(2, 0, 0, 1),
                "T4(i,j,k,l)=T4(i,j,k,l)(2,0,0,1)");
  test_for_zero(t4(2, 0, 0, 2) - t4_1(2, 0, 0, 2),
                "T4(i,j,k,l)=T4(i,j,k,l)(2,0,0,2)");
  test_for_zero(t4(2, 0, 1, 0) - t4_1(2, 0, 1, 0),
                "T4(i,j,k,l)=T4(i,j,k,l)(2,0,1,0)");
  test_for_zero(t4(2, 0, 1, 1) - t4_1(2, 0, 1, 1),
                "T4(i,j,k,l)=T4(i,j,k,l)(2,0,1,1)");
  test_for_zero(t4(2, 0, 1, 2) - t4_1(2, 0, 1, 2),
                "T4(i,j,k,l)=T4(i,j,k,l)(2,0,1,2)");
  test_for_zero(t4(2, 0, 2, 0) - t4_1(2, 0, 2, 0),
                "T4(i,j,k,l)=T4(i,j,k,l)(2,0,2,0)");
  test_for_zero(t4(2, 0, 2, 1) - t4_1(2, 0, 2, 1),
                "T4(i,j,k,l)=T4(i,j,k,l)(2,0,2,1)");
  test_for_zero(t4(2, 0, 2, 2) - t4_1(2, 0, 2, 2),
                "T4(i,j,k,l)=T4(i,j,k,l)(2,0,2,2)");
  test_for_zero(t4(2, 1, 0, 0) - t4_1(2, 1, 0, 0),
                "T4(i,j,k,l)=T4(i,j,k,l)(2,1,0,0)");
  test_for_zero(t4(2, 1, 0, 1) - t4_1(2, 1, 0, 1),
                "T4(i,j,k,l)=T4(i,j,k,l)(2,1,0,1)");
  test_for_zero(t4(2, 1, 0, 2) - t4_1(2, 1, 0, 2),
                "T4(i,j,k,l)=T4(i,j,k,l)(2,1,0,2)");
  test_for_zero(t4(2, 1, 1, 0) - t4_1(2, 1, 1, 0),
                "T4(i,j,k,l)=T4(i,j,k,l)(2,1,1,0)");
  test_for_zero(t4(2, 1, 1, 1) - t4_1(2, 1, 1, 1),
                "T4(i,j,k,l)=T4(i,j,k,l)(2,1,1,1)");
  test_for_zero(t4(2, 1, 1, 2) - t4_1(2, 1, 1, 2),
                "T4(i,j,k,l)=T4(i,j,k,l)(2,1,1,2)");
  test_for_zero(t4(2, 1, 2, 0) - t4_1(2, 1, 2, 0),
                "T4(i,j,k,l)=T4(i,j,k,l)(2,1,2,0)");
  test_for_zero(t4(2, 1, 2, 1) - t4_1(2, 1, 2, 1),
                "T4(i,j,k,l)=T4(i,j,k,l)(2,1,2,1)");
  test_for_zero(t4(2, 1, 2, 2) - t4_1(2, 1, 2, 2),
                "T4(i,j,k,l)=T4(i,j,k,l)(2,1,2,2)");
  test_for_zero(t4(2, 2, 0, 0) - t4_1(2, 2, 0, 0),
                "T4(i,j,k,l)=T4(i,j,k,l)(2,2,0,0)");
  test_for_zero(t4(2, 2, 0, 1) - t4_1(2, 2, 0, 1),
                "T4(i,j,k,l)=T4(i,j,k,l)(2,2,0,1)");
  test_for_zero(t4(2, 2, 0, 2) - t4_1(2, 2, 0, 2),
                "T4(i,j,k,l)=T4(i,j,k,l)(2,2,0,2)");
  test_for_zero(t4(2, 2, 1, 0) - t4_1(2, 2, 1, 0),
                "T4(i,j,k,l)=T4(i,j,k,l)(2,2,1,0)");
  test_for_zero(t4(2, 2, 1, 1) - t4_1(2, 2, 1, 1),
                "T4(i,j,k,l)=T4(i,j,k,l)(2,2,1,1)");
  test_for_zero(t4(2, 2, 1, 2) - t4_1(2, 2, 1, 2),
                "T4(i,j,k,l)=T4(i,j,k,l)(2,2,1,2)");
  test_for_zero(t4(2, 2, 2, 0) - t4_1(2, 2, 2, 0),
                "T4(i,j,k,l)=T4(i,j,k,l)(2,2,2,0)");
  test_for_zero(t4(2, 2, 2, 1) - t4_1(2, 2, 2, 1),
                "T4(i,j,k,l)=T4(i,j,k,l)(2,2,2,1)");
  test_for_zero(t4(2, 2, 2, 2) - t4_1(2, 2, 2, 2),
                "T4(i,j,k,l)=T4(i,j,k,l)(2,2,2,2)");

  /*---------------------------------------------------------------------------
   */
  t4(i, j, k, l) = t4_1(i, l, j, k);
  /*---------------------------------------------------------------------------
   */

  test_for_zero(t4(0, 0, 0, 0) - t4_1(0, 0, 0, 0),
                "T4(i,j,k,l)=T4(i,l,j,k)(0,0,0,0)");
  test_for_zero(t4(0, 0, 0, 1) - t4_1(0, 1, 0, 0),
                "T4(i,j,k,l)=T4(i,l,j,k)(0,0,0,1)");
  test_for_zero(t4(0, 0, 0, 2) - t4_1(0, 2, 0, 0),
                "T4(i,j,k,l)=T4(i,l,j,k)(0,0,0,2)");
  test_for_zero(t4(0, 0, 1, 0) - t4_1(0, 0, 0, 1),
                "T4(i,j,k,l)=T4(i,l,j,k)(0,0,1,0)");
  test_for_zero(t4(0, 0, 1, 1) - t4_1(0, 1, 0, 1),
                "T4(i,j,k,l)=T4(i,l,j,k)(0,0,1,1)");
  test_for_zero(t4(0, 0, 1, 2) - t4_1(0, 2, 0, 1),
                "T4(i,j,k,l)=T4(i,l,j,k)(0,0,1,2)");
  test_for_zero(t4(0, 0, 2, 0) - t4_1(0, 0, 0, 2),
                "T4(i,j,k,l)=T4(i,l,j,k)(0,0,2,0)");
  test_for_zero(t4(0, 0, 2, 1) - t4_1(0, 1, 0, 2),
                "T4(i,j,k,l)=T4(i,l,j,k)(0,0,2,1)");
  test_for_zero(t4(0, 0, 2, 2) - t4_1(0, 2, 0, 2),
                "T4(i,j,k,l)=T4(i,l,j,k)(0,0,2,2)");
  test_for_zero(t4(0, 1, 0, 0) - t4_1(0, 0, 1, 0),
                "T4(i,j,k,l)=T4(i,l,j,k)(0,1,0,0)");
  test_for_zero(t4(0, 1, 0, 1) - t4_1(0, 1, 1, 0),
                "T4(i,j,k,l)=T4(i,l,j,k)(0,1,0,1)");
  test_for_zero(t4(0, 1, 0, 2) - t4_1(0, 2, 1, 0),
                "T4(i,j,k,l)=T4(i,l,j,k)(0,1,0,2)");
  test_for_zero(t4(0, 1, 1, 0) - t4_1(0, 0, 1, 1),
                "T4(i,j,k,l)=T4(i,l,j,k)(0,1,1,0)");
  test_for_zero(t4(0, 1, 1, 1) - t4_1(0, 1, 1, 1),
                "T4(i,j,k,l)=T4(i,l,j,k)(0,1,1,1)");
  test_for_zero(t4(0, 1, 1, 2) - t4_1(0, 2, 1, 1),
                "T4(i,j,k,l)=T4(i,l,j,k)(0,1,1,2)");
  test_for_zero(t4(0, 1, 2, 0) - t4_1(0, 0, 1, 2),
                "T4(i,j,k,l)=T4(i,l,j,k)(0,1,2,0)");
  test_for_zero(t4(0, 1, 2, 1) - t4_1(0, 1, 1, 2),
                "T4(i,j,k,l)=T4(i,l,j,k)(0,1,2,1)");
  test_for_zero(t4(0, 1, 2, 2) - t4_1(0, 2, 1, 2),
                "T4(i,j,k,l)=T4(i,l,j,k)(0,1,2,2)");
  test_for_zero(t4(0, 2, 0, 0) - t4_1(0, 0, 2, 0),
                "T4(i,j,k,l)=T4(i,l,j,k)(0,2,0,0)");
  test_for_zero(t4(0, 2, 0, 1) - t4_1(0, 1, 2, 0),
                "T4(i,j,k,l)=T4(i,l,j,k)(0,2,0,1)");
  test_for_zero(t4(0, 2, 0, 2) - t4_1(0, 2, 2, 0),
                "T4(i,j,k,l)=T4(i,l,j,k)(0,2,0,2)");
  test_for_zero(t4(0, 2, 1, 0) - t4_1(0, 0, 2, 1),
                "T4(i,j,k,l)=T4(i,l,j,k)(0,2,1,0)");
  test_for_zero(t4(0, 2, 1, 1) - t4_1(0, 1, 2, 1),
                "T4(i,j,k,l)=T4(i,l,j,k)(0,2,1,1)");
  test_for_zero(t4(0, 2, 1, 2) - t4_1(0, 2, 2, 1),
                "T4(i,j,k,l)=T4(i,l,j,k)(0,2,1,2)");
  test_for_zero(t4(0, 2, 2, 0) - t4_1(0, 0, 2, 2),
                "T4(i,j,k,l)=T4(i,l,j,k)(0,2,2,0)");
  test_for_zero(t4(0, 2, 2, 1) - t4_1(0, 1, 2, 2),
                "T4(i,j,k,l)=T4(i,l,j,k)(0,2,2,1)");
  test_for_zero(t4(0, 2, 2, 2) - t4_1(0, 2, 2, 2),
                "T4(i,j,k,l)=T4(i,l,j,k)(0,2,2,2)");

  test_for_zero(t4(1, 0, 0, 0) - t4_1(1, 0, 0, 0),
                "T4(i,j,k,l)=T4(i,l,j,k)(1,0,0,0)");
  test_for_zero(t4(1, 0, 0, 1) - t4_1(1, 1, 0, 0),
                "T4(i,j,k,l)=T4(i,l,j,k)(1,0,0,1)");
  test_for_zero(t4(1, 0, 0, 2) - t4_1(1, 2, 0, 0),
                "T4(i,j,k,l)=T4(i,l,j,k)(1,0,0,2)");
  test_for_zero(t4(1, 0, 1, 0) - t4_1(1, 0, 0, 1),
                "T4(i,j,k,l)=T4(i,l,j,k)(1,0,1,0)");
  test_for_zero(t4(1, 0, 1, 1) - t4_1(1, 1, 0, 1),
                "T4(i,j,k,l)=T4(i,l,j,k)(1,0,1,1)");
  test_for_zero(t4(1, 0, 1, 2) - t4_1(1, 2, 0, 1),
                "T4(i,j,k,l)=T4(i,l,j,k)(1,0,1,2)");
  test_for_zero(t4(1, 0, 2, 0) - t4_1(1, 0, 0, 2),
                "T4(i,j,k,l)=T4(i,l,j,k)(1,0,2,0)");
  test_for_zero(t4(1, 0, 2, 1) - t4_1(1, 1, 0, 2),
                "T4(i,j,k,l)=T4(i,l,j,k)(1,0,2,1)");
  test_for_zero(t4(1, 0, 2, 2) - t4_1(1, 2, 0, 2),
                "T4(i,j,k,l)=T4(i,l,j,k)(1,0,2,2)");
  test_for_zero(t4(1, 1, 0, 0) - t4_1(1, 0, 1, 0),
                "T4(i,j,k,l)=T4(i,l,j,k)(1,1,0,0)");
  test_for_zero(t4(1, 1, 0, 1) - t4_1(1, 1, 1, 0),
                "T4(i,j,k,l)=T4(i,l,j,k)(1,1,0,1)");
  test_for_zero(t4(1, 1, 0, 2) - t4_1(1, 2, 1, 0),
                "T4(i,j,k,l)=T4(i,l,j,k)(1,1,0,2)");
  test_for_zero(t4(1, 1, 1, 0) - t4_1(1, 0, 1, 1),
                "T4(i,j,k,l)=T4(i,l,j,k)(1,1,1,0)");
  test_for_zero(t4(1, 1, 1, 1) - t4_1(1, 1, 1, 1),
                "T4(i,j,k,l)=T4(i,l,j,k)(1,1,1,1)");
  test_for_zero(t4(1, 1, 1, 2) - t4_1(1, 2, 1, 1),
                "T4(i,j,k,l)=T4(i,l,j,k)(1,1,1,2)");
  test_for_zero(t4(1, 1, 2, 0) - t4_1(1, 0, 1, 2),
                "T4(i,j,k,l)=T4(i,l,j,k)(1,1,2,0)");
  test_for_zero(t4(1, 1, 2, 1) - t4_1(1, 1, 1, 2),
                "T4(i,j,k,l)=T4(i,l,j,k)(1,1,2,1)");
  test_for_zero(t4(1, 1, 2, 2) - t4_1(1, 2, 1, 2),
                "T4(i,j,k,l)=T4(i,l,j,k)(1,1,2,2)");
  test_for_zero(t4(1, 2, 0, 0) - t4_1(1, 0, 2, 0),
                "T4(i,j,k,l)=T4(i,l,j,k)(1,2,0,0)");
  test_for_zero(t4(1, 2, 0, 1) - t4_1(1, 1, 2, 0),
                "T4(i,j,k,l)=T4(i,l,j,k)(1,2,0,1)");
  test_for_zero(t4(1, 2, 0, 2) - t4_1(1, 2, 2, 0),
                "T4(i,j,k,l)=T4(i,l,j,k)(1,2,0,2)");
  test_for_zero(t4(1, 2, 1, 0) - t4_1(1, 0, 2, 1),
                "T4(i,j,k,l)=T4(i,l,j,k)(1,2,1,0)");
  test_for_zero(t4(1, 2, 1, 1) - t4_1(1, 1, 2, 1),
                "T4(i,j,k,l)=T4(i,l,j,k)(1,2,1,1)");
  test_for_zero(t4(1, 2, 1, 2) - t4_1(1, 2, 2, 1),
                "T4(i,j,k,l)=T4(i,l,j,k)(1,2,1,2)");
  test_for_zero(t4(1, 2, 2, 0) - t4_1(1, 0, 2, 2),
                "T4(i,j,k,l)=T4(i,l,j,k)(1,2,2,0)");
  test_for_zero(t4(1, 2, 2, 1) - t4_1(1, 1, 2, 2),
                "T4(i,j,k,l)=T4(i,l,j,k)(1,2,2,1)");
  test_for_zero(t4(1, 2, 2, 2) - t4_1(1, 2, 2, 2),
                "T4(i,j,k,l)=T4(i,l,j,k)(1,2,2,2)");

  test_for_zero(t4(2, 0, 0, 0) - t4_1(2, 0, 0, 0),
                "T4(i,j,k,l)=T4(i,l,j,k)(2,0,0,0)");
  test_for_zero(t4(2, 0, 0, 1) - t4_1(2, 1, 0, 0),
                "T4(i,j,k,l)=T4(i,l,j,k)(2,0,0,1)");
  test_for_zero(t4(2, 0, 0, 2) - t4_1(2, 2, 0, 0),
                "T4(i,j,k,l)=T4(i,l,j,k)(2,0,0,2)");
  test_for_zero(t4(2, 0, 1, 0) - t4_1(2, 0, 0, 1),
                "T4(i,j,k,l)=T4(i,l,j,k)(2,0,1,0)");
  test_for_zero(t4(2, 0, 1, 1) - t4_1(2, 1, 0, 1),
                "T4(i,j,k,l)=T4(i,l,j,k)(2,0,1,1)");
  test_for_zero(t4(2, 0, 1, 2) - t4_1(2, 2, 0, 1),
                "T4(i,j,k,l)=T4(i,l,j,k)(2,0,1,2)");
  test_for_zero(t4(2, 0, 2, 0) - t4_1(2, 0, 0, 2),
                "T4(i,j,k,l)=T4(i,l,j,k)(2,0,2,0)");
  test_for_zero(t4(2, 0, 2, 1) - t4_1(2, 1, 0, 2),
                "T4(i,j,k,l)=T4(i,l,j,k)(2,0,2,1)");
  test_for_zero(t4(2, 0, 2, 2) - t4_1(2, 2, 0, 2),
                "T4(i,j,k,l)=T4(i,l,j,k)(2,0,2,2)");
  test_for_zero(t4(2, 1, 0, 0) - t4_1(2, 0, 1, 0),
                "T4(i,j,k,l)=T4(i,l,j,k)(2,1,0,0)");
  test_for_zero(t4(2, 1, 0, 1) - t4_1(2, 1, 1, 0),
                "T4(i,j,k,l)=T4(i,l,j,k)(2,1,0,1)");
  test_for_zero(t4(2, 1, 0, 2) - t4_1(2, 2, 1, 0),
                "T4(i,j,k,l)=T4(i,l,j,k)(2,1,0,2)");
  test_for_zero(t4(2, 1, 1, 0) - t4_1(2, 0, 1, 1),
                "T4(i,j,k,l)=T4(i,l,j,k)(2,1,1,0)");
  test_for_zero(t4(2, 1, 1, 1) - t4_1(2, 1, 1, 1),
                "T4(i,j,k,l)=T4(i,l,j,k)(2,1,1,1)");
  test_for_zero(t4(2, 1, 1, 2) - t4_1(2, 2, 1, 1),
                "T4(i,j,k,l)=T4(i,l,j,k)(2,1,1,2)");
  test_for_zero(t4(2, 1, 2, 0) - t4_1(2, 0, 1, 2),
                "T4(i,j,k,l)=T4(i,l,j,k)(2,1,2,0)");
  test_for_zero(t4(2, 1, 2, 1) - t4_1(2, 1, 1, 2),
                "T4(i,j,k,l)=T4(i,l,j,k)(2,1,2,1)");
  test_for_zero(t4(2, 1, 2, 2) - t4_1(2, 2, 1, 2),
                "T4(i,j,k,l)=T4(i,l,j,k)(2,1,2,2)");
  test_for_zero(t4(2, 2, 0, 0) - t4_1(2, 0, 2, 0),
                "T4(i,j,k,l)=T4(i,l,j,k)(2,2,0,0)");
  test_for_zero(t4(2, 2, 0, 1) - t4_1(2, 1, 2, 0),
                "T4(i,j,k,l)=T4(i,l,j,k)(2,2,0,1)");
  test_for_zero(t4(2, 2, 0, 2) - t4_1(2, 2, 2, 0),
                "T4(i,j,k,l)=T4(i,l,j,k)(2,2,0,2)");
  test_for_zero(t4(2, 2, 1, 0) - t4_1(2, 0, 2, 1),
                "T4(i,j,k,l)=T4(i,l,j,k)(2,2,1,0)");
  test_for_zero(t4(2, 2, 1, 1) - t4_1(2, 1, 2, 1),
                "T4(i,j,k,l)=T4(i,l,j,k)(2,2,1,1)");
  test_for_zero(t4(2, 2, 1, 2) - t4_1(2, 2, 2, 1),
                "T4(i,j,k,l)=T4(i,l,j,k)(2,2,1,2)");
  test_for_zero(t4(2, 2, 2, 0) - t4_1(2, 0, 2, 2),
                "T4(i,j,k,l)=T4(i,l,j,k)(2,2,2,0)");
  test_for_zero(t4(2, 2, 2, 1) - t4_1(2, 1, 2, 2),
                "T4(i,j,k,l)=T4(i,l,j,k)(2,2,2,1)");
  test_for_zero(t4(2, 2, 2, 2) - t4_1(2, 2, 2, 2),
                "T4(i,j,k,l)=T4(i,l,j,k)(2,2,2,2)");
}
