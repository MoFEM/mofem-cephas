/* A helper class that allows simple initialization of the Tensor2_antisymmetric,
   but only if it has the correct number of elements. */

template<class T, int Tensor_Dim>
class Tensor2_antisymmetric_constructor;

template<class T>
class Tensor2_antisymmetric_constructor<T,2>
{
public:
  Tensor2_antisymmetric_constructor(T data[(1*2)/2], T d01)
  {
    data[0]=d01;
  }
};

template<class T>
class Tensor2_antisymmetric_constructor<T,3>
{
public:
  Tensor2_antisymmetric_constructor(T data[(2*3)/2], T d01, T d02, T d12)
  {
    data[0]=d01;
    data[1]=d02;
    data[2]=d12;
  }
};

template<class T>
class Tensor2_antisymmetric_constructor<T,4>
{
public:
  Tensor2_antisymmetric_constructor(T data[(3*4)/2], T d01, T d02, T d03,
                                    T d12, T d13, T d23)
  {
    data[0]=d01;
    data[1]=d02;
    data[2]=d03;
    data[3]=d12;
    data[4]=d13;
    data[5]=d23;
  }
};
