#include "../../../src/FTensor.hpp"
#include "../test_for_zero.hpp"
#include <iostream>
using namespace FTensor;
using namespace std;

void test_T3_004(const Tensor3<double, 3, 3, 3> &t3_1) {

  Index<'i', 3> i;
  Index<'j', 3> j;
  Index<'k', 3> k;

  {
    // Symmetrize of first two indices yields dg
    Dg<double, 3, 3> t_dg;
    t_dg(i, j, k) = t3_1(i, j, k) || t3_1(j, i, k);
    FTensor::Tensor3<double, 3, 3, 3> t_sym_3;
    t_sym_3(i, j, k) = t3_1(i, j, k) + t3_1(j, i, k);

    for (int ii = 0; ii != 3; ++ii)
      for (int jj = 0; jj != 3; ++jj)
        for (int kk = 0; kk != 3; ++kk) {
          test_for_zero(t_dg(ii, jj, kk) - t_sym_3(ii, jj, kk),
                        "T3(i,j,k)||T3(j,i,k)(" + to_string(ii) + "," +
                            to_string(jj) + "," + to_string(kk) + ")");
        }

  }

  // Tensor 3 times symmetric tensor 2 yields tensor 3 
  {
    Tensor3<double, 1, 3, 2> t_3_1;
    Tensor2<double, 3, 3> t_2_1;
    Tensor2_symmetric<double, 3> t_2s_1;
    Tensor2<double, 3, 3> t_2_2;
    Tensor3<double, 1, 3, 2> t_3_2;
    for (int ii = 0; ii != 1; ++ii)
      for (int jj = 0; jj != 3; ++jj) {
        for (int kk = 0; kk != 2; ++kk) {
          t_3_1(ii, jj, kk) = 1 + ii + 10. * jj + 100. * kk;
        }
      }
    for (int jj = 0; jj != 3; ++jj)
      for (int ll = 0; ll != 3; ++ll) {
        t_2_1(jj, ll) = 1 + jj + 10. * ll;
      }
    Index<'i', 1> i;
    Index<'j', 3> j;
    Index<'k', 2> k;
    Index<'l', 3> l;
    t_2s_1(j, l) = t_2_1(j, l) || t_2_1(l, j);
    t_2_2(j, l) = t_2_1(j, l) + t_2_1(l, j);
    t_3_2(i, l, k) =
        t_3_1(i, j, k) * t_2s_1(j, l) -  t_2s_1(j, l)*t_3_1(i, j, k);
    for (int ii = 0; ii != 1; ++ii)
      for (int ll = 0; ll != 3; ++ll)
        for (int kk = 0; kk != 2; ++kk) {
          test_for_zero(t_3_2(ii, ll, kk),
                        "T3(i,j,k)||T2s(j,l)(" + to_string(ii) + "," +
                            to_string(ll) + "," + to_string(kk) + ")");
        }
    t_3_2(i, l, k) =
        t_3_1(i, j, k) * t_2s_1(j, l) - t_2_2(j, l) * t_3_1(i, j, k);
    for (int ii = 0; ii != 1; ++ii)
      for (int ll = 0; ll != 3; ++ll)
        for (int kk = 0; kk != 2; ++kk) {
          test_for_zero(t_3_2(ii, ll, kk),
                        "T3(i,j,k)||T2s(j,l)(" + to_string(ii) + "," +
                            to_string(ll) + "," + to_string(kk) + ")");
        }
 }

   // Tensor 3 times symmetric tensor 2 yields tensor 3
 {
   Tensor3<double, 3, 3, 3> t_3_1;
   Tensor2<double, 3, 3> t_2_1;
   Tensor2_symmetric<double, 3> t_2s_1;
   Tensor2<double, 3, 3> t_2_2;
   Tensor3<double, 3, 3, 3> t_3_2;
   for (int ii = 0; ii != 3; ++ii)
     for (int jj = 0; jj != 3; ++jj) {
       for (int kk = 0; kk != 3; ++kk) {
         t_3_1(ii, jj, kk) = 1 + ii + 10. * jj + 100. * kk;
       }
     }
   for (int jj = 0; jj != 3; ++jj)
     for (int ll = 0; ll != 3; ++ll) {
       t_2_1(jj, ll) = 1 + jj + 10. * ll;
     }
   Index<'i', 3> i;
   Index<'j', 3> j;
   Index<'k', 3> k;
   Index<'l', 3> l;
   t_2s_1(i, j) = t_2_1(i, j) || t_2_1(j, i);
   t_2_2(i, j) = t_2_1(i, j) + t_2_1(j, i);
   t_3_2(l, j, k) =
       t_3_1(i, j, k) * t_2s_1(i, l) - t_2s_1(i, l) * t_3_1(i, j, k);
   for (int ll = 0; ll != 3; ++ll)
     for (int jj = 0; jj != 3; ++jj)
       for (int kk = 0; kk != 3; ++kk) {
         test_for_zero(t_3_2(ll, jj, kk),
                       "T3(i,j,k)||T2s(i,l)(" + to_string(ll) + "," +
                           to_string(jj) + "," + to_string(kk) + ")");
       }
   t_3_2(l, j, k) =
       t_3_1(i, j, k) * t_2s_1(i, l) - t_3_1(i, j, k) * t_2_2(i, l);

   for (int ll = 0; ll != 3; ++ll)
     for (int jj = 0; jj != 3; ++jj)
       for (int kk = 0; kk != 3; ++kk) {
         test_for_zero(t_3_2(ll, jj, kk),
                       "T3(i,j,k)||T2s(i,l)(" + to_string(ll) + "," +
                           to_string(jj) + "," + to_string(kk) + ")");
       }
 }

 // Tensor 3 times tensor 3 yield tensor 2 A(i,j,k)*B(j,i,l)
 {
   Tensor1<double, 3> t_1;
   Tensor2<double, 3, 3> t_2;
   for (int ii = 0; ii != 3; ++ii) {
     t_1(ii) = 1 + ii;
     for (int jj = 0; jj != 3; ++jj) {
       t_2(ii, jj) = 1 + ii + 10. * jj;
     }
   }
   Tensor3<double, 3, 3, 3> t_3_1; 
   Index<'i', 3> i;
   Index<'j', 3> j;
   Index<'k', 3> k;
   Index<'l', 3> l;
   t_3_1(i, j, k) = t_2(i, j) * t_1(k);
   Tensor2<double, 3, 3> t_2_1;
   t_2_1(k, l) = t_3_1(i, j, k) * t_3_1(j, i, l);
   Tensor2<double, 3, 3> t_2_2;
   t_2_2(k, l) = (t_2(i, j) * t_2(j, i)) * (t_1(k) * t_1(l));
   for (int ii = 0; ii != 3; ++ii)
     for (int jj = 0; jj != 3; ++jj) {
       test_for_zero(t_2_1(ii, jj) - t_2_2(ii, jj), "T3(i,j,k)*T3(j,i,l)(" +
                                                        to_string(ii) + "," +
                                                        to_string(jj) + ")");
     }
 }

 // Tensor 3 times tensor 3 yield tensor 2 A(j,l,k)*B(i,k,l)
 {
   Tensor1<double, 3> t_1;
   Tensor2<double, 3, 3> t_2;
   for (int ii = 0; ii != 3; ++ii) {
     t_1(ii) = 1 + ii;
     for (int jj = 0; jj != 3; ++jj) {
       t_2(ii, jj) = 1 + ii + 10. * jj;
     }
   }
   Tensor3<double, 3, 3, 3> t_3_1; 
   Index<'i', 3> i;
   Index<'j', 3> j;
   Index<'k', 3> k;
   Index<'l', 3> l;
   t_3_1(i, j, k) = t_1(i) * t_2(j, k);
   Tensor2<double, 3, 3> t_2_1;
   t_2_1(i, j) = t_3_1(j, k, l) * t_3_1(i, l, k);
   Tensor2<double, 3, 3> t_2_2;
   t_2_2(i, j) = (t_2(k, l) * t_2(l, k)) * (t_1(i) * t_1(j));
   for (int ii = 0; ii != 3; ++ii)
     for (int jj = 0; jj != 3; ++jj) {
       test_for_zero(t_2_1(ii, jj) - t_2_2(ii, jj), "T3(j,k,l)*T3(i,l,k)(" +
                                                        to_string(ii) + "," +
                                                        to_string(jj) + ")");
     }
 }

 // Tensor 3 times tensor 3 yield tensor 2 A(k,i,j)*B(l,i,j)
 {
   Tensor1<double, 3> t_1;
   Tensor2<double, 3, 3> t_2;
   for (int ii = 0; ii != 3; ++ii) {
     t_1(ii) = 1 + ii;
     for (int jj = 0; jj != 3; ++jj) {
       t_2(ii, jj) = 1 + ii + 10. * jj;
     }
   }
   Tensor3<double, 3, 3, 3> t_3_1; 
   Index<'i', 3> i;
   Index<'j', 3> j;
   Index<'k', 3> k;
   Index<'l', 3> l;
   t_3_1(i, j, k) = t_1(i) * t_2(j, k);
   Tensor2<double, 3, 3> t_2_1;
   t_2_1(k, l) = t_3_1(k, i, j) * t_3_1(l, i, j);
   Tensor2<double, 3, 3> t_2_2;
   t_2_2(k, l) = (t_2(i, j) * t_2(i, j)) * (t_1(k) * t_1(l));
   for (int kk = 0; kk != 3; ++kk)
     for (int ll = 0; ll != 3; ++ll) {
       test_for_zero(t_2_1(kk, ll) - t_2_2(kk, ll), "T3(k,i,j)*T3(l,i,j)(" +
                                                        to_string(kk) + "," +
                                                        to_string(ll) + ")");
     }
 }

 // Tensor 3 times tensor 3 yield tensor 2 A(i,j,k)*B(i,j,k)
 {
   Tensor1<double, 3> t_1;
   Tensor2<double, 3, 3> t_2;
   for (int ii = 0; ii != 3; ++ii) {
     t_1(ii) = 1 + ii;
     for (int jj = 0; jj != 3; ++jj) {
       t_2(ii, jj) = 1 + ii + 10. * jj;
     }
   }
   Tensor3<double, 3, 3, 3> t_3_1; 
   Index<'i', 3> i;
   Index<'j', 3> j;
   Index<'k', 3> k;
   Index<'l', 3> l;
   t_3_1(i, j, k) = t_2(i, j) * t_1(k);
   Tensor2<double, 3, 3> t_2_1;
   t_2_1(k, l) = t_3_1(i, j, k) * t_3_1(i, j, l);
   Tensor2<double, 3, 3> t_2_2;
   t_2_2(k, l) = (t_2(i, j) * t_2(i, j)) * (t_1(k) * t_1(l));
   for (int kk = 0; kk != 3; ++kk)
     for (int ll = 0; ll != 3; ++ll) {
       test_for_zero(t_2_1(kk, ll) - t_2_2(kk, ll), "T3(i,j,k)*T3(i,j,l)(" +
                                                        to_string(kk) + "," +
                                                        to_string(ll) + ")");
     }
 }


}