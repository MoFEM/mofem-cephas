/* A version for pointers, useful for previously
   constructed arrays. */

#pragma once

namespace FTensor
{

template <class T, int Tensor_Dim> class Tensor1<T *, Tensor_Dim> {
protected:

  const int inc;

  /* Note that the T *'s are mutable, so the pointer can change,
     allowing iterating over a array. */

  mutable T *restrict data[Tensor_Dim];

public:
  /* Initializations for varying numbers of elements, with each one
     defined for a particular Tensor_Dim.  To initialize a different
     dimension, just add the appropriate constructor and call to
     the Tensor1_constructor constructor. */
  Tensor1(T *d0, const int i = 1) : inc(i) {
    Tensor1_constructor<T * restrict, Tensor_Dim>(data, d0);
  }
  Tensor1(T *d0, T *d1, const int i = 1) : inc(i) {
    Tensor1_constructor<T * restrict, Tensor_Dim>(data, d0, d1);
  }
  Tensor1(T *d0, T *d1, T *d2, const int i = 1) : inc(i) {
    Tensor1_constructor<T * restrict, Tensor_Dim>(data, d0, d1, d2);
  }
  Tensor1(T *d0, T *d1, T *d2, T *d3, const int i = 1) : inc(i) {
    Tensor1_constructor<T * restrict, Tensor_Dim>(data, d0, d1, d2, d3);
  }
  Tensor1(T *d0, T *d1, T *d2, T *d3, T *d4, T *d5, const int i = 1) : inc(i) {
    Tensor1_constructor<T * restrict, Tensor_Dim>(data, d0, d1, d2, d3, d4, d5);
  }
  /* Initializations for varying numbers of elements. */
  template <class... U> Tensor1(U *... d) : data(d...), inc(1) {}

  /* Initialization from array */
  template <class U>
  Tensor1(std::array<U *, Tensor_Dim> &a, const int i = 1) : inc(i) {
    std::copy(a.begin(), a.end(), data);
  }

  Tensor1(const int i = 1) : inc(i) {}

  Tensor1<T, Tensor_Dim> normalize() {
    const Index<'a', Tensor_Dim> a;
    (*this)(a) /= l2();
    return *this;
  }

  T l2() const { return sqrt(l2_squared(Number<Tensor_Dim>())); }

  template <int Current_Dim> T l2_squared(const Number<Current_Dim> &) const {
    return (*data[Current_Dim - 1]) * (*data[Current_Dim - 1]) +
           l2_squared(Number<Current_Dim - 1>());
  }
  T l2_squared(const Number<1> &) const { return (*data[0]) * (*data[0]); }

  /* There are two operator(int)'s, one for non-consts that lets you
     change the value, and one for consts that doesn't. */

  T &operator()(const int N) {
#ifdef FTENSOR_DEBUG
      if(N >= Tensor_Dim || N < 0)
        {
          std::stringstream s;
          s << "Bad index in Tensor1<T*," << Tensor_Dim << ">.operator(" << N
            << ")" << std::endl;
          throw std::out_of_range(s.str());
        }
#endif
      return *data[N];
    }
    T operator()(const int N) const
    {
#ifdef FTENSOR_DEBUG
      if(N >= Tensor_Dim || N < 0)
        {
          std::stringstream s;
          s << "Bad index in Tensor1<T*," << Tensor_Dim << ">.operator(" << N
            << ") const" << std::endl;
          throw std::out_of_range(s.str());
        }
#endif
      return *data[N];
    }
    T *ptr(const int N) const
    {
#ifdef FTENSOR_DEBUG
      if(N >= Tensor_Dim || N < 0)
        {
          std::stringstream s;
          s << "Bad index in Tensor1<T*," << Tensor_Dim << ">.ptr(" << N << ")"
            << std::endl;
          throw std::out_of_range(s.str());
        }
#endif
      return data[N];
    }

    T *&ptr(const int N) {
#ifdef FTENSOR_DEBUG
      if (N >= Tensor_Dim || N < 0) {
        std::stringstream s;
        s << "Bad index in Tensor1<T*," << Tensor_Dim << ">.ptr(" << N << ")"
          << std::endl;
        throw std::out_of_range(s.str());
      }
#endif
      return data[N];
    }

    /* These operator()'s are the first part in constructing template
       expressions.  They can be used to slice off lower dimensional
       parts. They are not entirely safe, since you can accidentaly use a
       higher dimension than what is really allowed (like Dim=5). */

    template <char i, int Dim>
    Tensor1_Expr<Tensor1<T *, Tensor_Dim>, T, Dim, i>
    operator()(const Index<i, Dim> &index)
    {
      return Tensor1_Expr<Tensor1<T *, Tensor_Dim>, T, Dim, i>(*this);
    }

    template <char i, int Dim>
    Tensor1_Expr<const Tensor1<T *, Tensor_Dim>, T, Dim, i>
    operator()(const Index<i, Dim> &index) const
    {
      return Tensor1_Expr<const Tensor1<T *, Tensor_Dim>, T, Dim, i>(*this);
    }

    /* The ++ operator increments the pointer, not the number that the
       pointer points to.  This allows iterating over a grid. */

    const Tensor1 &operator++() const
    {
      for(int i = 0; i < Tensor_Dim; ++i)
        data[i] += inc;
      return *this;
    }

    private:

      /**
       * @brief Preventing casting on derived class
       * 
       * That can be source of errors
       * 
       */
      template <int I>
      Tensor1(const Tensor1<PackPtr<T *, I>, Tensor_Dim> &) = delete;
  };

  template <class T, int Tensor_Dim, int I>
  class Tensor1<PackPtr<T *, I>, Tensor_Dim> : public Tensor1<T *, Tensor_Dim> {

  public:
  
    /* Initializations for varying numbers of elements. */
    template <class... U> Tensor1(U *... d) : Tensor1<T *, Tensor_Dim>(d...) {}

    template <class U>
    Tensor1(std::array<U *, Tensor_Dim> &a) : Tensor1<T *, Tensor_Dim>(a) {}

    Tensor1(): Tensor1<T *, Tensor_Dim>() {}

    /* The ++ operator increments the pointer, not the number that the
       pointer points to.  This allows iterating over a grid. */

    const Tensor1 &operator++() const
    {
      for(int i = 0; i < Tensor_Dim; ++i)
        Tensor1<T *, Tensor_Dim>::data[i] += I;
      return *this;
    }
  };

  
}
