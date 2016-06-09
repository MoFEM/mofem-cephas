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

  //  t4(i,j,k,l) = t2_1(i,j)*t2_1(k,l);
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

  /* A(i,j,k,l)*B(l,m) */

  t4ddg_1(i,j,k,l)=t2s_1(i,j)*t2s_1(k,l);
  t4(i,j,k,m) = t4ddg_1(i,j,k,l)*t2_cpy(l,m);
  t2(k,m) = t2s_1(k,l)*t2s_1(l,m);
  t4_222(i,j,k,m) = t2_cpy(i,j)*t2(k,m);

  for(int ii = 0;ii!=3;ii++) {
    for(int jj = 0;jj!=3;jj++) {
      for(int kk = 0;kk!=3;kk++) {
        for(int ll = 0;ll!=3;ll++) {
          // std::cerr << t4(ii,jj,kk,ll) << " " << t4_222(ii,jj,kk,ll) << std::endl;
          test_for_zero(t4(ii,jj,kk,ll) - t4_222(ii,jj,kk,ll),"Tensor4_ddg_times_Tensor2_3_1");
        }
      }
    }
  }

  /* B(l,m)*A(i,j,k,l) */

  t4(i,j,k,m) = t2_cpy(l,m)*t4ddg_1(i,j,k,l);

  for(int ii = 0;ii!=3;ii++) {
    for(int jj = 0;jj!=3;jj++) {
      for(int kk = 0;kk!=3;kk++) {
        for(int ll = 0;ll!=3;ll++) {
          // std::cerr << t4(ii,jj,kk,ll) << " " << t4_222(ii,jj,kk,ll) << std::endl;
          test_for_zero(t4(ii,jj,kk,ll) - t4_222(ii,jj,kk,ll),"Tensor4_ddg_times_Tensor2_3_1");
        }
      }
    }
  }

  /* A(i,j,k,l)*B(m,l) */

  t4(i,j,k,m) = t4ddg_1(i,j,k,l)*t2_cpy(m,l);
  t2(k,m) = t2s_1(k,l)*t2s_1(m,l);
  t4_222(i,j,k,m) = t2_cpy(i,j)*t2(k,m);

  for(int ii = 0;ii!=3;ii++) {
    for(int jj = 0;jj!=3;jj++) {
      for(int kk = 0;kk!=3;kk++) {
        for(int ll = 0;ll!=3;ll++) {
          // std::cerr << t4(ii,jj,kk,ll) << " " << t4_222(ii,jj,kk,ll) << std::endl;
          test_for_zero(t4(ii,jj,kk,ll) - t4_222(ii,jj,kk,ll),"Tensor4_ddg_times_Tensor2_3_0");
        }
      }
    }
  }

  /* A(i,j,k,l)*B(j,m) */

  t4(i,m,k,l) = t4ddg_1(i,j,k,l)*t2_cpy(j,m);
  t2(i,m) = t2s_1(i,j)*t2s_1(j,m);
  t4_222(i,m,k,l) = t2(i,m)*t2_cpy(k,l);

  for(int ii = 0;ii!=3;ii++) {
    for(int jj = 0;jj!=3;jj++) {
      for(int kk = 0;kk!=3;kk++) {
        for(int ll = 0;ll!=3;ll++) {
          // std::cerr << t4(ii,jj,kk,ll) << " " << t4_222(ii,jj,kk,ll) << std::endl;
          test_for_zero(t4(ii,jj,kk,ll) - t4_222(ii,jj,kk,ll),"Tensor4_ddg_times_Tensor2_1_0");
        }
      }
    }
  }

  /* B(j,m)*A(i,j,k,l)*/

  t4(i,m,k,l) = t2_cpy(j,m)*t4ddg_1(i,j,k,l);

  for(int ii = 0;ii!=3;ii++) {
    for(int jj = 0;jj!=3;jj++) {
      for(int kk = 0;kk!=3;kk++) {
        for(int ll = 0;ll!=3;ll++) {
          // std::cerr << t4(ii,jj,kk,ll) << " " << t4_222(ii,jj,kk,ll) << std::endl;
          test_for_zero(t4(ii,jj,kk,ll) - t4_222(ii,jj,kk,ll),"Tensor4_ddg_times_Tensor2_1_0");
        }
      }
    }
  }

  /* A(i,j,k,l)*B(m,j) */

  t4(i,m,k,l) = t4ddg_1(i,j,k,l)*t2_cpy(m,j);
  t2(i,m) = t2s_1(i,j)*t2s_1(m,j);
  t4_222(i,m,k,l) = t2(i,m)*t2_cpy(k,l);

  for(int ii = 0;ii!=3;ii++) {
    for(int jj = 0;jj!=3;jj++) {
      for(int kk = 0;kk!=3;kk++) {
        for(int ll = 0;ll!=3;ll++) {
          // std::cout << t4(ii,jj,kk,ll) << " " << t4_222(ii,jj,kk,ll) << std::endl;
          test_for_zero(t4(ii,jj,kk,ll) - t4_222(ii,jj,kk,ll),"Tensor4_ddg_times_Tensor2_1_1");
        }
      }
    }
  }

  /* B(m,m)*A(i,j,k,l)*/

  t4(i,m,k,l) = t2_cpy(m,j)*t4ddg_1(i,j,k,l);

  for(int ii = 0;ii!=3;ii++) {
    for(int jj = 0;jj!=3;jj++) {
      for(int kk = 0;kk!=3;kk++) {
        for(int ll = 0;ll!=3;ll++) {
          // std::cerr << t4(ii,jj,kk,ll) << " " << t4_222(ii,jj,kk,ll) << std::endl;
          test_for_zero(t4(ii,jj,kk,ll) - t4_222(ii,jj,kk,ll),"Tensor4_ddg_times_Tensor2_1_1");
        }
      }
    }
  }

  /* A(i,j,k,l)*B(i,m) */

  t4(m,j,k,l) = t4ddg_1(i,j,k,l)*t2_cpy(i,m);
  t2(m,j) = t2s_1(i,j)*t2s_1(i,m);
  t4_222(m,j,k,l) = t2(m,j)*t2_cpy(k,l);

  for(int ii = 0;ii!=3;ii++) {
    for(int jj = 0;jj!=3;jj++) {
      for(int kk = 0;kk!=3;kk++) {
        for(int ll = 0;ll!=3;ll++) {
          // std::cout << t4(ii,jj,kk,ll) << " " << t4_222(ii,jj,kk,ll) << std::endl;
          test_for_zero(t4(ii,jj,kk,ll) - t4_222(ii,jj,kk,ll),"Tensor4_ddg_times_Tensor2_0_0");
        }
      }
    }
  }

  /* B(i,m)*A(i,j,k,l) */

  t4(m,j,k,l) = t2_cpy(i,m)*t4ddg_1(i,j,k,l);

  for(int ii = 0;ii!=3;ii++) {
    for(int jj = 0;jj!=3;jj++) {
      for(int kk = 0;kk!=3;kk++) {
        for(int ll = 0;ll!=3;ll++) {
          // std::cout << t4(ii,jj,kk,ll) << " " << t4_222(ii,jj,kk,ll) << std::endl;
          test_for_zero(t4(ii,jj,kk,ll) - t4_222(ii,jj,kk,ll),"Tensor4_ddg_times_Tensor2_0_0");
        }
      }
    }
  }

  /* A(i,j,k,l)*B(m,i) */

  t4(m,j,k,l) = t4ddg_1(i,j,k,l)*t2_cpy(m,i);
  t2(m,j) = t2s_1(i,j)*t2s_1(m,i);
  t4_222(m,j,k,l) = t2(m,j)*t2_cpy(k,l);

  for(int ii = 0;ii!=3;ii++) {
    for(int jj = 0;jj!=3;jj++) {
      for(int kk = 0;kk!=3;kk++) {
        for(int ll = 0;ll!=3;ll++) {
          // std::cout << t4(ii,jj,kk,ll) << " " << t4_222(ii,jj,kk,ll) << std::endl;
          test_for_zero(t4(ii,jj,kk,ll) - t4_222(ii,jj,kk,ll),"Tensor4_ddg_times_Tensor2_0_0");
        }
      }
    }
  }

  /* B(m,i)*A(i,j,k,l) */

  t4(m,j,k,l) = t2_cpy(m,i)*t4ddg_1(i,j,k,l);

  for(int ii = 0;ii!=3;ii++) {
    for(int jj = 0;jj!=3;jj++) {
      for(int kk = 0;kk!=3;kk++) {
        for(int ll = 0;ll!=3;ll++) {
          // std::cout << t4(ii,jj,kk,ll) << " " << t4_222(ii,jj,kk,ll) << std::endl;
          test_for_zero(t4(ii,jj,kk,ll) - t4_222(ii,jj,kk,ll),"Tensor4_ddg_times_Tensor2_0_0");
        }
      }
    }
  }

  /* Testing Tensor4 pointer */

  double t4d[81];
  Tensor4<double*,3,3,3,3> t4_pointer(
    &t4d[ 0],&t4d[ 1],&t4d[ 2],&t4d[ 3],&t4d[ 4],&t4d[ 5],&t4d[ 6],&t4d[ 7],&t4d[ 8],
    &t4d[ 9],&t4d[10],&t4d[11],&t4d[12],&t4d[13],&t4d[14],&t4d[15],&t4d[16],&t4d[17],
    &t4d[18],&t4d[19],&t4d[20],&t4d[21],&t4d[22],&t4d[23],&t4d[24],&t4d[25],&t4d[26],
    &t4d[27],&t4d[28],&t4d[29],&t4d[30],&t4d[31],&t4d[32],&t4d[33],&t4d[34],&t4d[35],
    &t4d[36],&t4d[37],&t4d[38],&t4d[39],&t4d[40],&t4d[41],&t4d[42],&t4d[43],&t4d[44],
    &t4d[45],&t4d[46],&t4d[47],&t4d[48],&t4d[49],&t4d[50],&t4d[51],&t4d[52],&t4d[53],
    &t4d[54],&t4d[55],&t4d[56],&t4d[57],&t4d[58],&t4d[59],&t4d[60],&t4d[61],&t4d[62],
    &t4d[63],&t4d[64],&t4d[65],&t4d[66],&t4d[67],&t4d[68],&t4d[69],&t4d[70],&t4d[71],
    &t4d[72],&t4d[73],&t4d[74],&t4d[75],&t4d[76],&t4d[77],&t4d[78],&t4d[79],&t4d[80]
  );

  t4_pointer(i,j,k,l) = t4(i,j,k,l);

  for(int ii = 0;ii!=3;ii++) {
    for(int jj = 0;jj!=3;jj++) {
      for(int kk = 0;kk!=3;kk++) {
        for(int ll = 0;ll!=3;ll++) {
          // std::cout << ii << jj << kk << ll << std::endl;
          test_for_zero(t4_pointer(ii,jj,kk,ll) - t4(ii,jj,kk,ll),"Tensor4_pointer");
        }
      }
    }
  }

  Tensor4<double*,3,3,3,3> t4_pointer_2(t4d,1);
  t4_pointer_2(i,j,k,l) = t4(i,j,k,l);

  for(int ii = 0;ii!=3;ii++) {
    for(int jj = 0;jj!=3;jj++) {
      for(int kk = 0;kk!=3;kk++) {
        for(int ll = 0;ll!=3;ll++) {
          // std::cout << ii << jj << kk << ll << std::endl;
          test_for_zero(t4_pointer(ii,jj,kk,ll) - t4(ii,jj,kk,ll),"Tensor4_pointer");
        }
      }
    }
  }

  /* T4 times T2 yielding T4 */

  /* A(i,j,k,l)*B(l,m) */

  t4(i,j,k,m) = ( t2_1(i,j)*t2_1(k,l) )*t2_1(l,m);
  t4_222(i,j,k,m) = t2_1(i,j)*( t2_1(k,l)*t2_1(l,m) );

  for(int ii = 0;ii!=3;ii++) {
    for(int jj = 0;jj!=3;jj++) {
      for(int kk = 0;kk!=3;kk++) {
        for(int ll = 0;ll!=3;ll++) {
          // std::cerr << t4(ii,jj,kk,ll) << " " << t4_222(ii,jj,kk,ll) << std::endl;
          test_for_zero(t4(ii,jj,kk,ll) - t4_222(ii,jj,kk,ll),"Tensor4_times_Tensor2_3_1");
        }
      }
    }
  }

  /* B(l,m)*A(i,j,k,l) */

  t4(i,j,k,m) = t2_1(l,m)*( t2_1(i,j)*t2_1(k,l) );

  for(int ii = 0;ii!=3;ii++) {
    for(int jj = 0;jj!=3;jj++) {
      for(int kk = 0;kk!=3;kk++) {
        for(int ll = 0;ll!=3;ll++) {
          // std::cerr << t4(ii,jj,kk,ll) << " " << t4_222(ii,jj,kk,ll) << std::endl;
          test_for_zero(t4(ii,jj,kk,ll) - t4_222(ii,jj,kk,ll),"Tensor4_times_Tensor2_3_1");
        }
      }
    }
  }

  /* A(i,j,k,l)*B(m,l) */

  t4(i,j,k,m) = ( t2_1(i,j)*t2_1(k,l) )*t2_1(m,l);
  t4_222(i,j,k,m) = t2_1(i,j)*( t2_1(k,l)*t2_1(m,l) );

  for(int ii = 0;ii!=3;ii++) {
    for(int jj = 0;jj!=3;jj++) {
      for(int kk = 0;kk!=3;kk++) {
        for(int ll = 0;ll!=3;ll++) {
          // std::cerr << t4(ii,jj,kk,ll) << " " << t4_222(ii,jj,kk,ll) << std::endl;
          test_for_zero(t4(ii,jj,kk,ll) - t4_222(ii,jj,kk,ll),"Tensor4_times_Tensor2_3_0");
        }
      }
    }
  }

  /* B(m,l)*A(i,j,k,l) */

  t4(i,j,k,m) = t2_1(m,l)*( t2_1(i,j)*t2_1(k,l) );

  for(int ii = 0;ii!=3;ii++) {
    for(int jj = 0;jj!=3;jj++) {
      for(int kk = 0;kk!=3;kk++) {
        for(int ll = 0;ll!=3;ll++) {
          // std::cerr << t4(ii,jj,kk,ll) << " " << t4_222(ii,jj,kk,ll) << std::endl;
          test_for_zero(t4(ii,jj,kk,ll) - t4_222(ii,jj,kk,ll),"Tensor4_times_Tensor2_3_0");
        }
      }
    }
  }

  /* A(i,j,k,l)*B(j,m) */

  t4_222(i,j,k,l) = t2_1(i,j)*t2_1(k,l);
  t4(i,m,k,l) = t4_222(i,j,k,l)*t2_1(j,m);
  t4_222(i,m,k,l) = (t2_1(i,j)*t2_1(j,m))*t2_1(k,l);

  for(int ii = 0;ii!=3;ii++) {
    for(int jj = 0;jj!=3;jj++) {
      for(int kk = 0;kk!=3;kk++) {
        for(int ll = 0;ll!=3;ll++) {
          // std::cerr << t4(ii,jj,kk,ll) << " " << t4_222(ii,jj,kk,ll) << std::endl;
          test_for_zero(t4(ii,jj,kk,ll) - t4_222(ii,jj,kk,ll),"Tensor4_times_Tensor2_1_0");
        }
      }
    }
  }

  /* A(i,j,k,l)*B(m,j) */

  t4_222(i,j,k,l) = t2_1(i,j)*t2_1(k,l);
  t4(i,m,k,l) = t4_222(i,j,k,l)*t2_1(m,j);
  t4_222(i,m,k,l) = (t2_1(i,j)*t2_1(m,j))*t2_1(k,l);

  for(int ii = 0;ii!=3;ii++) {
    for(int jj = 0;jj!=3;jj++) {
      for(int kk = 0;kk!=3;kk++) {
        for(int ll = 0;ll!=3;ll++) {
          // std::cerr << t4(ii,jj,kk,ll) << " " << t4_222(ii,jj,kk,ll) << std::endl;
          test_for_zero(t4(ii,jj,kk,ll) - t4_222(ii,jj,kk,ll),"Tensor4_times_Tensor2_1_1");
        }
      }
    }
  }

  /* A(i,j,k,l)*B(i,m) */

  t4_222(i,j,k,l) = t2_1(i,j)*t2_1(k,l);
  t4(m,j,k,l) = t4_222(i,j,k,l)*t2_1(i,m);
  t4_222(m,j,k,l) = (t2_1(i,j)*t2_1(i,m))*t2_1(k,l);

  for(int ii = 0;ii!=3;ii++) {
    for(int jj = 0;jj!=3;jj++) {
      for(int kk = 0;kk!=3;kk++) {
        for(int ll = 0;ll!=3;ll++) {
          // std::cerr << t4(ii,jj,kk,ll) << " " << t4_222(ii,jj,kk,ll) << std::endl;
          test_for_zero(t4(ii,jj,kk,ll) - t4_222(ii,jj,kk,ll),"Tensor4_times_Tensor2_0_0");
        }
      }
    }
  }

  /* A(i,j,k,l)*B(m,i) */

  t4_222(i,j,k,l) = t2_1(i,j)*t2_1(k,l);
  t4(m,j,k,l) = t4_222(i,j,k,l)*t2_1(m,i);
  t4_222(m,j,k,l) = (t2_1(i,j)*t2_1(m,i))*t2_1(k,l);

  for(int ii = 0;ii!=3;ii++) {
    for(int jj = 0;jj!=3;jj++) {
      for(int kk = 0;kk!=3;kk++) {
        for(int ll = 0;ll!=3;ll++) {
          // std::cerr << t4(ii,jj,kk,ll) << " " << t4_222(ii,jj,kk,ll) << std::endl;
          test_for_zero(t4(ii,jj,kk,ll) - t4_222(ii,jj,kk,ll),"Tensor4_times_Tensor2_0_1");
        }
      }
    }
  }

  /* A(i,j,k,l)*B(k,m) */

  t4_222(i,j,k,l) = t2_1(i,j)*t2_1(k,l);
  t4(i,j,m,l) = t4_222(i,j,k,l)*t2_1(k,m);
  t4_222(i,j,m,l) = t2_1(i,j)*(t2_1(k,l)*t2_1(k,m));

  for(int ii = 0;ii!=3;ii++) {
    for(int jj = 0;jj!=3;jj++) {
      for(int kk = 0;kk!=3;kk++) {
        for(int ll = 0;ll!=3;ll++) {
          // std::cerr << t4(ii,jj,kk,ll) << " " << t4_222(ii,jj,kk,ll) << std::endl;
          test_for_zero(t4(ii,jj,kk,ll) - t4_222(ii,jj,kk,ll),"Tensor4_times_Tensor2_2_0");
        }
      }
    }
  }

  /* B(m,k)*A(i,j,k,l) */

  t4_222(i,j,k,l) = t2_1(i,j)*t2_1(k,l);
  t4(i,j,m,l) = t4_222(i,j,k,l)*t2_1(k,m);
  t4_222(i,j,m,l) = t2_1(i,j)*(t2_1(k,l)*t2_1(k,m));

  for(int ii = 0;ii!=3;ii++) {
    for(int jj = 0;jj!=3;jj++) {
      for(int kk = 0;kk!=3;kk++) {
        for(int ll = 0;ll!=3;ll++) {
          // std::cerr << t4(ii,jj,kk,ll) << " " << t4_222(ii,jj,kk,ll) << std::endl;
          test_for_zero(t4(ii,jj,kk,ll) - t4_222(ii,jj,kk,ll),"Tensor4_times_Tensor2_2_1");
        }
      }
    }
  }

  /* T4*=U */

  t4(i,j,m,l) *= 4;

  for(int ii = 0;ii!=3;ii++) {
    for(int jj = 0;jj!=3;jj++) {
      for(int kk = 0;kk!=3;kk++) {
        for(int ll = 0;ll!=3;ll++) {
          // std::cerr << t4(ii,jj,kk,ll) << " " << t4_222(ii,jj,kk,ll) << std::endl;
          test_for_zero(t4(ii,jj,kk,ll) - 4*t4_222(ii,jj,kk,ll),"Tensor4_times_Tensor2_2_1");
        }
      }
    }
  }

}
