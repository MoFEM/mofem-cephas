/* A version for pointers. */

#pragma once

namespace FTensor {
template <class T, int Dim0, int Dim1, int Dim2, int Current_Dim0,
          int Current_Dim1, int Current_Dim2>
inline void T3_increment(const Tensor3<T, Dim0, Dim1, Dim2> &iter,
                         const Number<Current_Dim0> &,
                         const Number<Current_Dim1> &,
                         const Number<Current_Dim2> &) {
  iter.increment(Number<Current_Dim0>(), Number<Current_Dim1>(),
                 Number<Current_Dim2>());
  T3_increment(iter, Number<Current_Dim0 - 1>(), Number<Current_Dim1>(),
               Number<Current_Dim2>());
}

template <class T, int Dim0, int Dim1, int Dim2, int Current_Dim1,
          int Current_Dim2>
inline void T3_increment(const Tensor3<T, Dim0, Dim1, Dim2> &iter,
                         const Number<1> &, const Number<Current_Dim1> &,
                         const Number<Current_Dim2> &) {
  iter.increment(Number<1>(), Number<Current_Dim1>(), Number<Current_Dim2>());
  T3_increment(iter, Number<Dim0>(), Number<Current_Dim1 - 1>(),
               Number<Current_Dim2>());
}

template <class T, int Dim0, int Dim1, int Dim2, int Current_Dim2>
inline void T3_increment(const Tensor3<T, Dim0, Dim1, Dim2> &iter,
                         const Number<1> &, const Number<1> &,
                         const Number<Current_Dim2> &) {
  iter.increment(Number<1>(), Number<1>(), Number<Current_Dim2>());
  T3_increment(iter, Number<Dim0>(), Number<Dim1>(),
               Number<Current_Dim2 - 1>());
}

template <class T, int Dim0, int Dim1, int Dim2>
inline void T3_increment(const Tensor3<T, Dim0, Dim1, Dim2> &iter,
                         const Number<1> &, const Number<1> &,
                         const Number<1> &) {
  iter.increment(Number<1>(), Number<1>(), Number<1>());
}

template <class T, int Tensor_Dim0, int Tensor_Dim1, int Tensor_Dim2>
class Tensor3<T *, Tensor_Dim0, Tensor_Dim1, Tensor_Dim2> {
  const int inc;

protected:
  mutable T *restrict data[Tensor_Dim0][Tensor_Dim1][Tensor_Dim2];

public:
  Tensor3() {}

  /* Tensor_Dim0=2, Tensor_Dim1=2, Tensor_Dim2=2 */
  Tensor3(T *d000, T *d001, T *d010, T *d011, T *d100, T *d101, T *d110,
          T *d111, const int i)
      : inc(i) {
    Tensor3_constructor<T *restrict, Tensor_Dim0, Tensor_Dim1, Tensor_Dim2>(
        data, d000, d001, d010, d011, d100, d101, d110, d111);
  }

  /* Tensor_Dim0=3, Tensor_Dim1=3, Tensor_Dim2=3 */
  Tensor3(T *d000, T *d001, T *d002, T *d010, T *d011, T *d012, T *d020,
          T *d021, T *d022, T *d100, T *d101, T *d102, T *d110, T *d111,
          T *d112, T *d120, T *d121, T *d122, T *d200, T *d201, T *d202,
          T *d210, T *d211, T *d212, T *d220, T *d221, T *d222, const int i)
      : inc(i) {
    Tensor3_constructor<T *restrict, Tensor_Dim0, Tensor_Dim1, Tensor_Dim2>(
        data, d000, d001, d002, d010, d011, d012, d020, d021, d022, d100, d101,
        d102, d110, d111, d112, d120, d121, d122, d200, d201, d202, d210, d211,
        d212, d220, d221, d222);
  }

  /* Tensor_Dim0=4, Tensor_Dim1=4, Tensor_Dim2=4 */
  Tensor3(T *d000, T *d001, T *d002, T *d003, T *d010, T *d011, T *d012,
          T *d013, T *d020, T *d021, T *d022, T *d023, T *d030, T *d031,
          T *d032, T *d033,

          T *d100, T *d101, T *d102, T *d103, T *d110, T *d111, T *d112,
          T *d113, T *d120, T *d121, T *d122, T *d123, T *d130, T *d131,
          T *d132, T *d133,

          T *d200, T *d201, T *d202, T *d203, T *d210, T *d211, T *d212,
          T *d213, T *d220, T *d221, T *d222, T *d223, T *d230, T *d231,
          T *d232, T *d233,

          T *d300, T *d301, T *d302, T *d303, T *d310, T *d311, T *d312,
          T *d313, T *d320, T *d321, T *d322, T *d323, T *d330, T *d331,
          T *d332, T *d333, const int i)
      : inc(i) {
    Tensor3_constructor<T *restrict, Tensor_Dim0, Tensor_Dim1, Tensor_Dim2>(
        data, d000, d001, d002, d003, d010, d011, d012, d013, d020, d021, d022,
        d023, d030, d031, d032, d033,

        d100, d101, d102, d103, d110, d111, d112, d113, d120, d121, d122, d123,
        d130, d131, d132, d133,

        d200, d201, d202, d203, d210, d211, d212, d213, d220, d221, d222, d223,
        d230, d231, d232, d233,

        d300, d301, d302, d303, d310, d311, d312, d313, d320, d321, d322, d323,
        d330, d331, d332, d333);
  }

  /* Initializations for varying numbers of elements. */
  template <class... U> Tensor3(U *...d) : data{d...}, inc(1) {}

  template <class U>
  Tensor3(std::array<U *, Tensor_Dim0 * Tensor_Dim1 * Tensor_Dim2> &a,
          const int i = 1)
      : inc(i) {
    int l = 0;
    for (int i = 0; i != Tensor_Dim0; ++i) {
      for (int j = 0; j != Tensor_Dim1; ++j) {
        for (int k = 0; k != Tensor_Dim2; ++k, ++l) {
          data[i][j][k] = a[l];
        }
      }
    }
  }

  /* There are two operator(int,int,int)'s, one for non-consts that lets you
     change the value, and one for consts that doesn't. */

  T &operator()(const int N1, const int N2, const int N3) {
#ifdef FTENSOR_DEBUG
    if (N1 >= Tensor_Dim0 || N1 < 0 || N2 >= Tensor_Dim1 || N2 < 0 ||
        N3 >= Tensor_Dim2 || N3 < 0) {
      std::stringstream s;
      s << "Bad index in Tensor3<T," << Tensor_Dim0 << "," << Tensor_Dim1 << ","
        << Tensor_Dim2 << ">.operator(" << N1 << "," << N2 << "," << N3 << ")"
        << std::endl;
      throw std::runtime_error(s.str());
    }
#endif
    return *data[N1][N2][N3];
  }

  T operator()(const int N1, const int N2, const int N3) const {
#ifdef FTENSOR_DEBUG
    if (N1 >= Tensor_Dim0 || N1 < 0 || N2 >= Tensor_Dim1 || N2 < 0 ||
        N3 >= Tensor_Dim2 || N3 < 0) {
      std::stringstream s;
      s << "Bad index in Tensor3<T," << Tensor_Dim0 << "," << Tensor_Dim1 << ","
        << Tensor_Dim2 << ">.operator(" << N1 << "," << N2 << "," << N3
        << ") const" << std::endl;
      throw std::runtime_error(s.str());
    }
#endif
    return *data[N1][N2][N3];
  }

  /* These operator()'s are the first part in constructing template
     expressions. */

  template <char i, char j, char k, int Dim0, int Dim1, int Dim2>
  Tensor3_Expr<Tensor3<T *, Tensor_Dim0, Tensor_Dim1, Tensor_Dim2>, T, Dim0,
               Dim1, Dim2, i, j, k>
  operator()(const Index<i, Dim0>, const Index<j, Dim1>, const Index<k, Dim2>) {
    return Tensor3_Expr<Tensor3<T *, Tensor_Dim0, Tensor_Dim1, Tensor_Dim2>, T,
                        Dim0, Dim1, Dim2, i, j, k>(*this);
  }

  template <char i, char j, char k, int Dim0, int Dim1, int Dim2>
  Tensor3_Expr<const Tensor3<T *, Tensor_Dim0, Tensor_Dim1, Tensor_Dim2>, T,
               Dim0, Dim1, Dim2, i, j, k>
  operator()(const Index<i, Dim0>, const Index<j, Dim1>,
             const Index<k, Dim2>) const {
    return Tensor3_Expr<
        const Tensor3<T *, Tensor_Dim0, Tensor_Dim1, Tensor_Dim2>, T, Dim0,
        Dim1, Dim2, i, j, k>(*this);
  }

  /* The ++ operator increments the pointer, not the number that the
     pointer points to.  This allows iterating over a grid. */

  template <int Current_Dim0, int Current_Dim1, int Current_Dim2>
  inline void increment(const Number<Current_Dim0> &,
                        const Number<Current_Dim1> &,
                        const Number<Current_Dim2> &) const {
    data[Current_Dim0 - 1][Current_Dim1 - 1][Current_Dim2 - 1] += inc;
  }

  const Tensor3 &operator++() const {
    T3_increment(*this, Number<Tensor_Dim0>(), Number<Tensor_Dim1>(),
                 Number<Tensor_Dim2>());
    return *this;
  }

private:
  template <int I>
  Tensor3(const Tensor3<PackPtr<T *, I>, Tensor_Dim0, Tensor_Dim1, Tensor_Dim2>
              &) = delete;
};

template <class T, int Tensor_Dim0, int Tensor_Dim1, int Tensor_Dim2, int I>
class Tensor3<PackPtr<T *, I>, Tensor_Dim0, Tensor_Dim1, Tensor_Dim2>
    : public Tensor3<T *, Tensor_Dim0, Tensor_Dim1, Tensor_Dim2> {

public:
  Tensor3() {}

  /* Initializations for varying numbers of elements. */
  template <class... U>
  Tensor3(U *...d)
      : Tensor3<T *, Tensor_Dim0, Tensor_Dim1, Tensor_Dim2>{d...} {}

  template <class U>
  Tensor3(std::array<U *, Tensor_Dim0 * Tensor_Dim1 * Tensor_Dim2> &a)
      : Tensor3<T *, Tensor_Dim0, Tensor_Dim1, Tensor_Dim2>(a) {}

  /* The ++ operator increments the pointer, not the number that the
     pointer points to.  This allows iterating over a grid. */

  template <int Current_Dim0, int Current_Dim1, int Current_Dim2>
  inline void increment(const Number<Current_Dim0> &,
                        const Number<Current_Dim1> &,
                        const Number<Current_Dim2> &) const {
    Tensor3<T *, Tensor_Dim0, Tensor_Dim1,
            Tensor_Dim2>::data[Current_Dim0 - 1][Current_Dim1 - 1]
                              [Current_Dim2 - 1] += I;
  }

  const Tensor3 &operator++() const {
    T3_increment(*this, Number<Tensor_Dim0>(), Number<Tensor_Dim1>(),
                 Number<Tensor_Dim2>());
    return *this;
  }
};
} // namespace FTensor