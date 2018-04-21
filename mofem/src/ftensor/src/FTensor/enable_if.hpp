/* Workaround to pre c++11 compilers for std::enable_if */

#pragma once

namespace std
{
  template <bool Condition, typename T = void> struct enable_if
  {};

  template <typename T> struct enable_if<true, T>
  {
    typedef T type;
  };
}
