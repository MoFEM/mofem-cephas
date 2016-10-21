#include <stdlib.h>
#include <iostream>
#include <cmath>

// #ifdef FTENOR_MOFEMTESTING

#ifndef __TEST_FOR_ZERO_HPP__
#define __TEST_FOR_ZERO_HPP__

template<class T>
void test_for_zero(const T &t, const char *s)
{

  if(std::abs(t)<1e-14) {
    std::cout << "PASS: " << s << std::endl;
  } else {
    std::cout << "FAIL: " << s << " " << std::abs(t) << std::endl;
    exit(-1);
  }
}

#endif // __TEST_FOR_ZERO_HPP__
