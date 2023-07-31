/* A helper class that allows simple initialization of the Tensor1,
   but only if it has the correct number of elements. */

#pragma once

namespace FTensor
{
  template <class T, int Tensor_Dim> class Tensor1_constructor;

  template <class T> class Tensor1_constructor<T, 1>
  {
  public:
    Tensor1_constructor(T data[], T d0) { data[0] = d0; }
  };

  template <class T> class Tensor1_constructor<T, 2>
  {
  public:
    Tensor1_constructor(T data[], T d0, T d1)
    {
      data[0] = d0;
      data[1] = d1;
    }
  };

  template <class T> class Tensor1_constructor<T, 3>
  {
  public:
    Tensor1_constructor(T data[], T d0, T d1, T d2)
    {
      data[0] = d0;
      data[1] = d1;
      data[2] = d2;
    }
  };

  template <class T> class Tensor1_constructor<T, 4>
  {
  public:
    Tensor1_constructor(T data[], T d0, T d1, T d2, T d3)
    {
      data[0] = d0;
      data[1] = d1;
      data[2] = d2;
      data[3] = d3;
    }
  };
  template <class T> class Tensor1_constructor<T, 6>
  {
  public:
    Tensor1_constructor(T data[], T d0, T d1, T d2, T d3, T d4, T d5)
    {
      data[0] = d0;
      data[1] = d1;
      data[2] = d2;
      data[3] = d3;
      data[4] = d4;
      data[5] = d5;
    }
  };
  template <class T> class Tensor1_constructor<T, 9> {
  public:
    Tensor1_constructor(T data[], T d0, T d1, T d2, T d3, T d4, T d5, T d6,
                        T d7, T d8) {
      data[0] = d0;
      data[1] = d1;
      data[2] = d2;
      data[3] = d3;
      data[4] = d4;
      data[5] = d5;
      data[6] = d6;
      data[7] = d7;
      data[8] = d8;
    }
  };
}
