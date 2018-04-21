
#ifndef WITH_ADOL_C
#error "MoFEM need to be compiled with ADOL-C"
#endif

#include <ostream>

#define ADOLC_TAPELESS
#define NUMBER_DIRECTIONS 6

// #include <adolc/adolc.h>
#include <adolc/adtl.h>
// typedef adtl::adouble adouble;

#define FTENSOR_DEBUG
#include <FTensor.hpp>

int main(int argc, char *argv[])
{
  adtl::setNumDir(6);
  FTensor::Index<'I', 2> I;
  FTensor::Index<'J', 2> J;

  FTensor::Tensor1<double, 2> t1(1, 1);
  FTensor::Tensor2<double, 2, 2> t2(1, 0, 0, 1);

  FTensor::Tensor1<adtl::adouble, 2> a_t1;
  FTensor::Tensor2<adtl::adouble, 2, 2> a_t2;

  int dd0 = 0;
  for(int ii = 0; ii != 2; ii++, dd0++)
    {
      a_t1(ii).setValue(t1(ii));
      for(int kk = 0; kk != 6; kk++)
        {
          if(dd0 == kk)
            {
              a_t1(ii).setADValue(kk, 1);
            }
          else
            {
              a_t1(ii).setADValue(kk, 0);
            }
        }
    }
  for(int ii = 0; ii != 2; ii++)
    {
      for(int jj = 0; jj != 2; jj++, dd0++)
        {
          a_t2(ii, jj).setValue(t2(ii, jj));
          for(int kk = 0; kk != 6; kk++)
            {
              if(kk == dd0)
                {
                  a_t2(ii, jj).setADValue(kk, 1);
                }
              else
                {
                  a_t2(ii, jj).setADValue(kk, 0);
                }
            }
        }
    }
  for(int ii = 0; ii != 2; ii++)
    {
      std::cout << "a_t1 ( " << ii << " ) = " << a_t1(ii) << std::endl;
    }
  for(int ii = 0; ii != 2; ii++)
    {
      for(int jj = 0; jj != 2; jj++)
        {
          std::cout << "a_t1 ( " << ii << "," << jj << " ) = " << a_t2(ii, jj)
                    << std::endl;
        }
    }

  adtl::adouble a_t0 = a_t1(I) * a_t2(I, J) * a_t1(J);

  double t0;
  t0 = a_t0.getValue();
  std::cout << "Value: " << t0 << " ( " << a_t0 << " ) " << std::endl;

  FTensor::Tensor1<double, 2> ad_t0_t1;
  FTensor::Tensor2<double, 2, 2> ad_t0_t2;

  int dd = 0;
  for(int ii = 0; ii != 2; ii++, dd++)
    {
      ad_t0_t1(ii) = a_t0.getADValue(dd);
    }
  for(int ii = 0; ii != 2; ii++)
    {
      for(int jj = 0; jj != 2; jj++, dd++)
        {
          ad_t0_t2(ii, jj) = a_t0.getADValue(dd);
        }
    }

  // 2nd derivative
  FTensor::Tensor2<double, 2, 2> ad_t0_t1_t1;

  std::cout << "Derivatives t0_t1" << std::endl;
  for(int ii = 0; ii != 2; ii++, dd++)
    {
      std::cout << ad_t0_t1(ii) << std::endl;
      if(ad_t0_t1(ii) != 2)
        {
          std::cerr << "Wrong result, should be 2" << std::endl;
          exit(-1);
        }
    }
  std::cout << "Derivatives t0_t2" << std::endl;
  for(int ii = 0; ii != 2; ii++)
    {
      for(int jj = 0; jj != 2; jj++, dd++)
        {
          std::cout << ad_t0_t2(ii, jj) << " ";
          if(ad_t0_t2(ii, ii) != 1)
            {
              std::cerr << "Wrong result, should be 1" << std::endl;
              exit(-1);
            }
        }
      std::cout << "\n";
    }

  return 0;
}
