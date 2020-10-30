
#ifndef WITH_ADOL_C
#error "MoFEM need to be compiled with ADOL-C"
#endif

#include <ostream>

#include <adolc/adolc.h>

#define FTENSOR_DEBUG
#include <FTensor.hpp>

#define BOOST_UBLAS_SHALLOW_ARRAY_ADAPTOR

#include <boost/numeric/ublas/storage.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/io.hpp>

int main(int argc, char *argv[]) {
  FTensor::Index<'I', 2> I;
  FTensor::Index<'J', 2> J;

  boost::numeric::ublas::vector<double> input(3);
  input.clear();
  FTensor::Tensor2_symmetric<double *, 2> t2(&input[0], &input[1], &input[2]);
  t2(0, 0) = 1;
  t2(0, 1) = 0.5;
  t2(1, 1) = 1;

  FTensor::Tensor2_symmetric<double, 2> one;
  one(0, 0) = one(1, 1) = 1;
  one(1, 0) = one(0, 1) = 1;

  for (int ii = 0; ii != 2; ii++) {
    for (int jj = 0; jj != 2; jj++) {
      std::cout << t2(ii, jj) << " ";
    }
    std::cout << std::endl;
  }
  std::cout << std::endl;

  double r;

  int tag = 1;
  int keep = 1;

  trace_on(tag, keep);
  FTensor::Tensor2_symmetric<adouble, 2> a_t2;
  // FTensor::Tensor2_symmetric<adouble,2> a_t2_tmp;
  a_t2(I, J) <<= t2(I, J);
  // a_t2_tmp = a_t2;
  adouble a_r = a_t2(I, J) * a_t2(I, J);
  // a_t2(0,0)*a_t2(0,0)+a_t2(1,1)*a_t2(1,1)+a_t2(0,1)*a_t2(0,1);
  a_r >>= r;
  trace_off();

  for (int ii = 0; ii != 2; ii++) {
    for (int jj = 0; jj != 2; jj++) {
      std::cout << a_t2(ii, jj) << " ";
    }
    std::cout << std::endl;
  }
  std::cout << std::endl;

  std::cout << a_r << " " << r << std::endl << std::endl;

  boost::numeric::ublas::vector<double> output(3);
  FTensor::Tensor2_symmetric<double *, 2> jac_t2(&output[0], &output[1],
                                                 &output[2]);

  /*int ret_jac =*/::gradient(tag, 3, &input[0], &output[0]);

  std::cout << input << std::endl;
  std::cout << output << std::endl;

  for (int ii = 0; ii != 2; ii++) {
    for (int jj = 0; jj != 2; jj++) {
      std::cout << jac_t2(ii, jj) << " ";
      if (jac_t2(ii, jj) != 2) {
        std::cerr << "Wrong value should be 2\n" << std::endl;
        exit(-1);
      }
    }
    std::cout << std::endl;
  }
  std::cout << std::endl;

  boost::numeric::ublas::matrix<double> Hessian(3, 3);
  Hessian.clear();
  double *H[3];
  for (int nn = 0; nn != 3; nn++) {
    H[nn] = &Hessian(nn, 0);
  }
  /*int ret_val_hessian =*/::hessian(tag, 3, &input[0], H);
  std::cout << Hessian << std::endl;

  FTensor::Ddg<double *, 2, 2> t4(
      &Hessian(0, 0), &Hessian(0, 1), &Hessian(0, 2), &Hessian(1, 0),
      &Hessian(1, 1), &Hessian(1, 2), &Hessian(2, 0), &Hessian(2, 1),
      &Hessian(2, 2));
  for (int ii = 0; ii != 2; ii++) {
    for (int jj = 0; jj != 2; jj++) {
      for (int II = 0; II != 2; II++) {
        for (int JJ = 0; JJ != 2; JJ++) {
          std::cout << "(" << ii << "," << jj << "," << II << "," << JJ << ") "
                    << t4(ii, jj, II, JJ) << std::endl;
        }
      }
    }
  }
  std::cout << std::endl;

  if (t4(0, 0, 0, 0) != 2 || t4(1, 1, 1, 1) != 2 || t4(1, 0, 1, 0) != 4 ||
      t4(0, 1, 0, 1) != 4 || t4(0, 1, 1, 0) != 4 || t4(1, 0, 0, 1) != 4) {
    std::cerr << "Wrong value\n" << std::endl;
    exit(-1);
  }

  return 0;
}
