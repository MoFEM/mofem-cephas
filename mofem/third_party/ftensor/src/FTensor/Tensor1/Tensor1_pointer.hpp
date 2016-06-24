/* A version for pointers, useful for previously
   constructed arrays. */

template<class T,int Dim,int Current_Dim>
inline void T1_increment(
  const Tensor1<T,Dim> &iter,const Number<Current_Dim> &
) {
  iter.increment(Number<Current_Dim>());
  T1_increment(iter,Number<Current_Dim-1>());
}

template<class T,int Dim>
inline void T1_increment(
  const Tensor1<T,Dim> &iter,const Number<1> &
) {
  iter.increment(Number<1>());
}

template <class T, int Tensor_Dim>
class Tensor1<T*,Tensor_Dim>
{
protected:
  /* Note that the T *'s are mutable, so the pointer can change,
     allowing iterating over a array. */
  const int inc;
  mutable T * restrict data[Tensor_Dim];
public:
  /* Initializations for varying numbers of elements, with each one
     defined for a particular Tensor_Dim.  To initialize a different
     dimension, just add the appropriate constructor and call to
     the Tensor1_constructor constructor. */
  Tensor1(T *d0, T *d1,const int i = 1): inc(i)
  {
    Tensor1_constructor<T* restrict,Tensor_Dim>(data,d0,d1);
  }
  Tensor1(T *d0, T *d1, T *d2,const int i = 1): inc(i)
  {
    Tensor1_constructor<T* restrict,Tensor_Dim>(data,d0,d1,d2);
  }
  Tensor1(T *d0, T *d1, T *d2, T *d3,const int i = 1): inc(i)
  {
    Tensor1_constructor<T* restrict,Tensor_Dim>(data,d0,d1,d2,d3);
  }

  /* There are two operator(int)'s, one for non-consts that lets you
     change the value, and one for consts that doesn't. */

  T & operator()(const int N) {
    #ifdef FTENSOR_DEBUG
    if(N>=Tensor_Dim || N<0)
    {
      std::stringstream s;
      s << "Bad index in Tensor1<T*," << Tensor_Dim
      << ">.operator(" << N << ")" << std::endl;
      throw std::runtime_error(s.str());
    }
    if(!data[N]) {
      std::stringstream s;
      s << "Null pointer in Tensor1<T*," << Tensor_Dim
      << ">.operator(" << N << ")" << std::endl;
      throw std::runtime_error(s.str());
    }
    #endif
    return *data[N];
  }
  T operator()(const int N) const {
    #ifdef FTENSOR_DEBUG
    if(N>=Tensor_Dim || N<0)
    {
      std::stringstream s;
      s << "Bad index in Tensor1<T*," << Tensor_Dim
      << ">.operator(" << N << ") const" << std::endl;
      throw std::runtime_error(s.str());
    }
    #endif
    return data[N] ? *data[N] : 0;
  }
  T * ptr(const int N) const {
    #ifdef FTENSOR_DEBUG
    if(N>=Tensor_Dim || N<0)
    {
      std::stringstream s;
      s << "Bad index in Tensor1<T*," << Tensor_Dim
      << ">.ptr(" << N << ")" << std::endl;
      throw std::runtime_error(s.str());
    }
    #endif
    return data[N];
  }

  /* These operator()'s are the first part in constructing template
     expressions.  They can be used to slice off lower dimensional
     parts. They are not entirely safe, since you can accidently use a
     higher dimension than what is really allowed (like Dim=5). */

  template<char i, int Dim>
  Tensor1_Expr<Tensor1<T*,Tensor_Dim>,T,Dim,i>
  operator()(const Index<i,Dim> &index)
  {
    return Tensor1_Expr<Tensor1<T*,Tensor_Dim>,T,Dim,i>(*this);
  }

  template<char i, int Dim>
  Tensor1_Expr<const Tensor1<T*,Tensor_Dim>,T,Dim,i>
  operator()(const Index<i,Dim> &index) const
  {
    return Tensor1_Expr<const Tensor1<T*,Tensor_Dim>,T,Dim,i>(*this);
  }

  /* The ++ operator increments the pointer, not the number that the
     pointer points to.  This allows iterating over a grid. */

  template<int Current_Dim>
  inline void increment(const Number<Current_Dim> &) const {
    if(data[Current_Dim-1]) data[Current_Dim-1]+=inc;
  }

  const Tensor1 & operator++() const
  {
    T1_increment(*this,Number<Tensor_Dim>());
    return *this;
  }

};

class AAA {
	void A() {
		double d[9];
		Tensor1<double*,3> t2(
			&d[0],&d[1],&d[2]
		);
		++t2;
	}
};
