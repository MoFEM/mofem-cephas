#include <iostream>
#include <cmath>

#pragma once

#ifndef __TEST_FOR_ZERO_HPP__
#define __TEST_FOR_ZERO_HPP__

template<class T>
static void test_for_zero(const T &t, const char *s)
{
  if(std::abs(t)<1e-14)
    std::cout << "PASS: " << s << "\n";
  else
    std::cout << "FAIL: " << s << " " << std::abs(t) << "\n";
}

#endif // __TEST_FOR_ZERO_HPP__