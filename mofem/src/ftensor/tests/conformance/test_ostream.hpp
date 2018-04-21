#include <iostream>
#include <sstream>

#ifndef __TES_OSTREAM_HPP__
#define __TES_OSTREAM_HPP__

template <class T>
static void test_ostream(const T &t, const std::string &expected,
                         const std::string &test_name)
{
  std::stringstream ss;
  ss << t;
  if(ss.str() == expected)
    {
      std::cout << "PASS: " << test_name << "\n";
    }
  else
    {
      std::cout << "FAIL: " << test_name << ": " << ss.str()
                << "!=" << expected << "\n";
    }
}

#endif //__TES_OSTREAM_HPP__