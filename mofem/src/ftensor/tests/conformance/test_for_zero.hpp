#include <stdlib.h>
#include <iostream>
#include <cmath>

// #ifdef FTENOR_MOFEMTESTING

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

// #else
//
// template<class T>
// void test_for_zero(const T &t, const char *s)
// {
//   if(std::abs(t)<1e-14)
//   std::cout << "PASS: " << s << std::endl;
//   else
//   std::cout << "FAIL: " << s << " " << std::abs(t) << std::endl;
// }
//
// #endif
