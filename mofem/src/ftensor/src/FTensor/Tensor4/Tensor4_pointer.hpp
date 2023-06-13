/* A version for pointers. */

#pragma once

namespace FTensor
{
template <class T, int Dim0, int Dim1, int Dim2, int Dim3, int Current_Dim0,
          int Current_Dim1, int Current_Dim2, int Current_Dim3>
inline void
T4_increment(const Tensor4<T, Dim0, Dim1, Dim2, Dim3> &iter,
             const Number<Current_Dim0> &, const Number<Current_Dim1> &,
             const Number<Current_Dim2> &, const Number<Current_Dim3> &)
{
  iter.increment(Number<Current_Dim0>(), Number<Current_Dim1>(),
                 Number<Current_Dim2>(), Number<Current_Dim3>());
  T4_increment(iter, Number<Current_Dim0 - 1>(), Number<Current_Dim1>(),
               Number<Current_Dim2>(), Number<Current_Dim3>());
}

template <class T, int Dim0, int Dim1, int Dim2, int Dim3, int Current_Dim1,
          int Current_Dim2, int Current_Dim3>
inline void
T4_increment(const Tensor4<T, Dim0, Dim1, Dim2, Dim3> &iter,
             const Number<1> &, const Number<Current_Dim1> &,
             const Number<Current_Dim2> &, const Number<Current_Dim3> &)
{
  iter.increment(Number<1>(), Number<Current_Dim1>(), Number<Current_Dim2>(),
                 Number<Current_Dim3>());
  T4_increment(iter, Number<Dim0>(), Number<Current_Dim1 - 1>(),
               Number<Current_Dim2>(), Number<Current_Dim3>());
}

template <class T, int Dim0, int Dim1, int Dim2, int Dim3, int Current_Dim2,
          int Current_Dim3>
inline void
T4_increment(const Tensor4<T, Dim0, Dim1, Dim2, Dim3> &iter,
             const Number<1> &, const Number<1> &,
             const Number<Current_Dim2> &, const Number<Current_Dim3> &)
{
  iter.increment(Number<1>(), Number<1>(), Number<Current_Dim2>(),
                 Number<Current_Dim3>());
  T4_increment(iter, Number<Dim0>(), Number<Dim1>(),
               Number<Current_Dim2 - 1>(), Number<Current_Dim3>());
}

template <class T, int Dim0, int Dim1, int Dim2, int Dim3, int Current_Dim3>
inline void T4_increment(const Tensor4<T, Dim0, Dim1, Dim2, Dim3> &iter,
                         const Number<1> &, const Number<1> &,
                         const Number<1> &, const Number<Current_Dim3> &)
{
  iter.increment(Number<1>(), Number<1>(), Number<1>(),
                 Number<Current_Dim3>());
  T4_increment(iter, Number<Dim0>(), Number<Dim1>(), Number<Dim2>(),
               Number<Current_Dim3 - 1>());
}

template <class T, int Dim0, int Dim1, int Dim2, int Dim3>
inline void T4_increment(const Tensor4<T, Dim0, Dim1, Dim2, Dim3> &iter,
                         const Number<1> &, const Number<1> &,
                         const Number<1> &, const Number<1> &)
{
  iter.increment(Number<1>(), Number<1>(), Number<1>(), Number<1>());
}

template <class T, int Tensor_Dim0, int Tensor_Dim1, int Tensor_Dim2,
          int Tensor_Dim3>
class Tensor4<T *, Tensor_Dim0, Tensor_Dim1, Tensor_Dim2, Tensor_Dim3>
{
  const int inc;

protected:  
  mutable T *restrict
      data[Tensor_Dim0][Tensor_Dim1][Tensor_Dim2][Tensor_Dim3];

public:
  Tensor4(T *d0000, T *d0001, T *d0010, T *d0011, T *d0100, T *d0101,
          T *d0110, T *d0111, T *d1000, T *d1001, T *d1010, T *d1011,
          T *d1100, T *d1101, T *d1110, T *d1111, const int i = 1)
      : inc(i)
  {
    Tensor4_constructor<T * restrict, Tensor_Dim0, Tensor_Dim1, Tensor_Dim2,
                        Tensor_Dim3>(
        data, d0000, d0001, d0010, d0011, d0100, d0101, d0110, d0111, d1000,
        d1001, d1010, d1011, d1100, d1101, d1110, d1111);
  }

  Tensor4(T *d0000, T *d0001, T *d0002, T *d0010, T *d0011, T *d0012,
          T *d0020, T *d0021, T *d0022, T *d0100, T *d0101, T *d0102,
          T *d0110, T *d0111, T *d0112, T *d0120, T *d0121, T *d0122,
          T *d0200, T *d0201, T *d0202, T *d0210, T *d0211, T *d0212,
          T *d0220, T *d0221, T *d0222, T *d1000, T *d1001, T *d1002,
          T *d1010, T *d1011, T *d1012, T *d1020, T *d1021, T *d1022,
          T *d1100, T *d1101, T *d1102, T *d1110, T *d1111, T *d1112,
          T *d1120, T *d1121, T *d1122, T *d1200, T *d1201, T *d1202,
          T *d1210, T *d1211, T *d1212, T *d1220, T *d1221, T *d1222,
          T *d2000, T *d2001, T *d2002, T *d2010, T *d2011, T *d2012,
          T *d2020, T *d2021, T *d2022, T *d2100, T *d2101, T *d2102,
          T *d2110, T *d2111, T *d2112, T *d2120, T *d2121, T *d2122,
          T *d2200, T *d2201, T *d2202, T *d2210, T *d2211, T *d2212,
          T *d2220, T *d2221, T *d2222, const int i = 1)
      : inc(i)
  {
    Tensor4_constructor<T * restrict, Tensor_Dim0, Tensor_Dim1, Tensor_Dim2,
                        Tensor_Dim3>(
        data, d0000, d0001, d0002, d0010, d0011, d0012, d0020, d0021, d0022,
        d0100, d0101, d0102, d0110, d0111, d0112, d0120, d0121, d0122, d0200,
        d0201, d0202, d0210, d0211, d0212, d0220, d0221, d0222, d1000, d1001,
        d1002, d1010, d1011, d1012, d1020, d1021, d1022, d1100, d1101, d1102,
        d1110, d1111, d1112, d1120, d1121, d1122, d1200, d1201, d1202, d1210,
        d1211, d1212, d1220, d1221, d1222, d2000, d2001, d2002, d2010, d2011,
        d2012, d2020, d2021, d2022, d2100, d2101, d2102, d2110, d2111, d2112,
        d2120, d2121, d2122, d2200, d2201, d2202, d2210, d2211, d2212, d2220,
        d2221, d2222);
  }

  Tensor4(T *d, const int shift, const int i = 1) : inc(i)
  {
    int s = 0;
    for (int i = 0; i != Tensor_Dim0; ++i)
    {
      for (int j = 0; j != Tensor_Dim1; ++j)
      {
        for (int k = 0; k != Tensor_Dim2; ++k)
        {
          for (int l = 0; l != Tensor_Dim3; ++l)
          {
            data[i][j][k][l] = &d[s];
            s += shift;
          }
        }
      }
    }
  }

  /* Initializations for varying numbers of elements. */
  template <class... U> Tensor4(U *... d) : inc(1), data{d...} {}

  /* There are two operator(int,int)'s, one for non-consts that lets you
       change the value, and one for consts that doesn't. */

  T &operator()(const int N1, const int N2, const int N3, const int N4)
  {
#ifdef FTENSOR_DEBUG
    if (N1 >= Tensor_Dim0 || N1 < 0 || N2 >= Tensor_Dim1 || N2 < 0 || N3 >= Tensor_Dim2 || N3 < 0 || N4 >= Tensor_Dim3 || N4 < 0)
    {
      std::stringstream s;
      s << "Bad index in Tensor4<T," << Tensor_Dim0 << "," << Tensor_Dim1
        << "," << Tensor_Dim2 << "," << Tensor_Dim3 << ">.operator(" << N1
        << "," << N2 << "," << N3 << "," << N4 << ") const" << std::endl;
      throw std::runtime_error(s.str());
    }
#endif
    return *data[N1][N2][N3][N4];
  }

  T operator()(const int N1, const int N2, const int N3, const int N4) const
  {
#ifdef FTENSOR_DEBUG
    if (N1 >= Tensor_Dim0 || N1 < 0 || N2 >= Tensor_Dim1 || N2 < 0 || N3 >= Tensor_Dim2 || N3 < 0 || N4 >= Tensor_Dim3 || N4 < 0)
    {
      std::stringstream s;
      s << "Bad index in Tensor4<T," << Tensor_Dim0 << "," << Tensor_Dim1
        << "," << Tensor_Dim2 << "," << Tensor_Dim3 << ">.operator(" << N1
        << "," << N2 << "," << N3 << "," << N4 << ") const" << std::endl;
      throw std::runtime_error(s.str());
    }
#endif
    return *data[N1][N2][N3][N4];
  }

  /* These operator()'s are the first part in constructing template
       expressions.  They can be used to slice off lower dimensional
       parts. They are not entirely safe, since you can accidently use a
       higher dimension than what is really allowed (like Dim=5). */

  template <char i, char j, char k, char l, int Dim0, int Dim1, int Dim2,
            int Dim3>
  Tensor4_Expr<
      Tensor4<T *, Tensor_Dim0, Tensor_Dim1, Tensor_Dim2, Tensor_Dim3>, T,
      Dim0, Dim1, Dim2, Dim3, i, j, k, l>
  operator()(const Index<i, Dim0>, const Index<j, Dim1>,
             const Index<k, Dim2>, const Index<l, Dim3>)
  {
    return Tensor4_Expr<
        Tensor4<T *, Tensor_Dim0, Tensor_Dim1, Tensor_Dim2, Tensor_Dim3>, T,
        Dim0, Dim1, Dim2, Dim3, i, j, k, l>(*this);
  }

  template <char i, char j, char k, char l, int Dim0, int Dim1, int Dim2,
            int Dim3>
  Tensor4_Expr<
      const Tensor4<T *, Tensor_Dim0, Tensor_Dim1, Tensor_Dim2, Tensor_Dim3>,
      T, Dim0, Dim1, Dim2, Dim3, i, j, k, l>
  operator()(const Index<i, Dim0>, const Index<j, Dim1>,
             const Index<k, Dim2>, const Index<l, Dim3>) const
  {
    return Tensor4_Expr<
        const Tensor4<T *, Tensor_Dim0, Tensor_Dim1, Tensor_Dim2, Tensor_Dim3>,
        T, Dim0, Dim1, Dim2, Dim3, i, j, k, l>(*this);
  }

  /* The ++ operator increments the pointer, not the number that the
       pointer points to.  This allows iterating over a grid. */

  template <int Current_Dim0, int Current_Dim1, int Current_Dim2,
            int Current_Dim3>
  inline void
  increment(const Number<Current_Dim0> &, const Number<Current_Dim1> &,
            const Number<Current_Dim2> &, const Number<Current_Dim3> &) const
  {
    data[Current_Dim0 - 1][Current_Dim1 - 1][Current_Dim2 - 1]
        [Current_Dim3 - 1] += inc;
  }

  const Tensor4 &operator++() const
  {
    T4_increment(*this, Number<Tensor_Dim0>(), Number<Tensor_Dim1>(),
                 Number<Tensor_Dim2>(), Number<Tensor_Dim3>());
    return *this;
  }

  T *ptr(const int N1, const int N2, const int N3, const int N4) const
  {
#ifdef FTENSOR_DEBUG
    if (N1 >= Tensor_Dim0 || N1 < 0 || N2 >= Tensor_Dim1 || N2 < 0 ||
        N3 >= Tensor_Dim2 || N3 < 0 || N4 >= Tensor_Dim3 || N4 < 0)
    {
      std::stringstream s;
      s << "Bad index in Tensor4<T*," << Tensor_Dim0 << "," << Tensor_Dim1
        << "," << Tensor_Dim2 << "," << Tensor_Dim3 << ">.ptr(" << N1 << ","
        << N2 << "," << N3 << "," << N4 << ")" << std::endl;
      throw std::out_of_range(s.str());
    }
#endif
    return data[N1][N2][N3][N4];
  }

  T *restrict &ptr(const int N1, const int N2, const int N3, const int N4)
  {
#ifdef FTENSOR_DEBUG
    if (N1 >= Tensor_Dim0 || N1 < 0 || N2 >= Tensor_Dim1 || N2 < 0 ||
        N3 >= Tensor_Dim2 || N3 < 0 || N4 >= Tensor_Dim3 || N4 < 0)
    {
      std::stringstream s;
      s << "Bad index in Tensor4<T*," << Tensor_Dim0 << "," << Tensor_Dim1
        << "," << Tensor_Dim2 << "," << Tensor_Dim3 << ">.ptr(" << N1 << ","
        << N2 << "," << N3 << "," << N4 << ")" << std::endl;
      throw std::out_of_range(s.str());
    }
#endif
    return data[N1][N2][N3][N4];
  }

  private:
    template <int I>
    Tensor4(const Tensor4<PackPtr<T *, I>, Tensor_Dim0, Tensor_Dim1,
                          Tensor_Dim2, Tensor_Dim3> &) = delete;
};

template <class T, int Tensor_Dim0, int Tensor_Dim1, int Tensor_Dim2,
          int Tensor_Dim3, int I>
class Tensor4<PackPtr<T *, I>, Tensor_Dim0, Tensor_Dim1, Tensor_Dim2,
              Tensor_Dim3> : public Tensor4<T *, Tensor_Dim0, Tensor_Dim1,
                                            Tensor_Dim2, Tensor_Dim3>
{

public:
  Tensor4(T *d, const int shift)
      : Tensor4<T *, Tensor_Dim0, Tensor_Dim1, Tensor_Dim2, Tensor_Dim3>(
            d, shift) {}

  template <class... U>
  Tensor4(U *... d)
      : Tensor4<T *, Tensor_Dim0, Tensor_Dim1, Tensor_Dim2, Tensor_Dim3>(
            d...) {}

  /* The ++ operator increments the pointer, not the number that the
       pointer points to.  This allows iterating over a grid. */

  template <int Current_Dim0, int Current_Dim1, int Current_Dim2,
            int Current_Dim3>
  inline void increment(const Number<Current_Dim0> &,
                        const Number<Current_Dim1> &,
                        const Number<Current_Dim2> &,
                        const Number<Current_Dim3> &) const
  {
    Tensor4<T *, Tensor_Dim0, Tensor_Dim1, Tensor_Dim2,
            Tensor_Dim3>::data[Current_Dim0 - 1][Current_Dim1 - 1]
                              [Current_Dim2 - 1][Current_Dim3 - 1] += I;
  }

  const Tensor4 &operator++() const
  {
    T4_increment(*this, Number<Tensor_Dim0>(), Number<Tensor_Dim1>(),
                 Number<Tensor_Dim2>(), Number<Tensor_Dim3>());
    return *this;
  }
};
} // namespace FTensor