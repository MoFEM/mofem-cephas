/* A version for pointers. */

#pragma once

namespace FTensor
{
  template <class T, int Tensor_Dim01, int Tensor_Dim23>
  class Ddg<T *, Tensor_Dim01, Tensor_Dim23>
  {
    const int inc;

  protected:

    mutable T *restrict data[(Tensor_Dim01 * (Tensor_Dim01 + 1)) / 2]
                            [(Tensor_Dim23 * (Tensor_Dim23 + 1)) / 2];

  public:
    Ddg(T *d0000, T *d0010, T *d0011, T *d1000, T *d1010, T *d1011, T *d1100,
        T *d1110, T *d1111, const int i = 1)
        : inc(i)
    {
      ptr(0, 0, 0, 0) = d0000;
      ptr(0, 0, 1, 0) = d0010;
      ptr(0, 0, 1, 1) = d0011;
      ptr(1, 0, 0, 0) = d1000;
      ptr(1, 0, 1, 0) = d1010;
      ptr(1, 0, 1, 1) = d1011;
      ptr(1, 1, 0, 0) = d1100;
      ptr(1, 1, 1, 0) = d1110;
      ptr(1, 1, 1, 1) = d1111;
    }

    Ddg(T *d0000, T *d0001, T *d0002, T *d0011, T *d0012, T *d0022, T *d0100,
        T *d0101, T *d0102, T *d0111, T *d0112, T *d0122, T *d0200, T *d0201,
        T *d0202, T *d0211, T *d0212, T *d0222, T *d1100, T *d1101, T *d1102,
        T *d1111, T *d1112, T *d1122, T *d1200, T *d1201, T *d1202, T *d1211,
        T *d1212, T *d1222, T *d2200, T *d2201, T *d2202, T *d2211, T *d2212,
        T *d2222, const int i = 1)
        : inc(i)
    {
      ptr(0, 0, 0, 0) = d0000;
      ptr(0, 0, 0, 1) = d0001;
      ptr(0, 0, 0, 2) = d0002;
      ptr(0, 0, 1, 1) = d0011;
      ptr(0, 0, 1, 2) = d0012;
      ptr(0, 0, 2, 2) = d0022;
      ptr(0, 1, 0, 0) = d0100;
      ptr(0, 1, 0, 1) = d0101;
      ptr(0, 1, 0, 2) = d0102;
      ptr(0, 1, 1, 1) = d0111;
      ptr(0, 1, 1, 2) = d0112;
      ptr(0, 1, 2, 2) = d0122;
      ptr(0, 2, 0, 0) = d0200;
      ptr(0, 2, 0, 1) = d0201;
      ptr(0, 2, 0, 2) = d0202;
      ptr(0, 2, 1, 1) = d0211;
      ptr(0, 2, 1, 2) = d0212;
      ptr(0, 2, 2, 2) = d0222;
      ptr(1, 1, 0, 0) = d1100;
      ptr(1, 1, 0, 1) = d1101;
      ptr(1, 1, 0, 2) = d1102;
      ptr(1, 1, 1, 1) = d1111;
      ptr(1, 1, 1, 2) = d1112;
      ptr(1, 1, 2, 2) = d1122;
      ptr(1, 2, 0, 0) = d1200;
      ptr(1, 2, 0, 1) = d1201;
      ptr(1, 2, 0, 2) = d1202;
      ptr(1, 2, 1, 1) = d1211;
      ptr(1, 2, 1, 2) = d1212;
      ptr(1, 2, 2, 2) = d1222;
      ptr(2, 2, 0, 0) = d2200;
      ptr(2, 2, 0, 1) = d2201;
      ptr(2, 2, 0, 2) = d2202;
      ptr(2, 2, 1, 1) = d2211;
      ptr(2, 2, 1, 2) = d2212;
      ptr(2, 2, 2, 2) = d2222;
    }

    /* Initializations for varying numbers of elements. */
    template <class... U> Ddg(U *... d) : data{d...}, inc(1) {}

    /* There are two operator(int,int,int,int)'s, one for non-consts
       that lets you change the value, and one for consts that
       doesn't. */

    T &operator()(const int N1, const int N2, const int N3, const int N4)
    {
#ifdef FTENSOR_DEBUG
      if(N1 >= Tensor_Dim01 || N1 < 0 || N2 >= Tensor_Dim01 || N2 < 0
         || N3 >= Tensor_Dim23 || N3 < 0 || N4 >= Tensor_Dim23 || N4 < 0)
        {
          std::stringstream s;
          s << "Bad index in Dg<T*," << Tensor_Dim01 << "," << Tensor_Dim23
            << ">.operator(" << N1 << "," << N2 << "," << N3 << "," << N4
            << ")" << std::endl;
          throw std::out_of_range(s.str());
        }
#endif
      return N1 > N2
               ? (N3 > N4 ? *data[N1 + (N2 * (2 * Tensor_Dim01 - N2 - 1)) / 2]
                                 [N3 + (N4 * (2 * Tensor_Dim23 - N4 - 1)) / 2]
                          : *data[N1 + (N2 * (2 * Tensor_Dim01 - N2 - 1)) / 2]
                                 [N4 + (N3 * (2 * Tensor_Dim23 - N3 - 1)) / 2])
               : (N3 > N4
                    ? *data[N2 + (N1 * (2 * Tensor_Dim01 - N1 - 1)) / 2]
                           [N3 + (N4 * (2 * Tensor_Dim23 - N4 - 1)) / 2]
                    : *data[N2 + (N1 * (2 * Tensor_Dim01 - N1 - 1)) / 2]
                           [N4 + (N3 * (2 * Tensor_Dim23 - N3 - 1)) / 2]);
    }

    T operator()(const int N1, const int N2, const int N3, const int N4) const
    {
#ifdef FTENSOR_DEBUG
      if(N1 >= Tensor_Dim01 || N1 < 0 || N2 >= Tensor_Dim01 || N2 < 0
         || N3 >= Tensor_Dim23 || N3 < 0 || N4 >= Tensor_Dim23 || N4 < 0)
        {
          std::stringstream s;
          s << "Bad index in Dg<T*," << Tensor_Dim01 << "," << Tensor_Dim23
            << ">.operator(" << N1 << "," << N2 << "," << N3 << "," << N4
            << ") const" << std::endl;
          throw std::out_of_range(s.str());
        }
#endif
      return N1 > N2
               ? (N3 > N4 ? *data[N1 + (N2 * (2 * Tensor_Dim01 - N2 - 1)) / 2]
                                 [N3 + (N4 * (2 * Tensor_Dim23 - N4 - 1)) / 2]
                          : *data[N1 + (N2 * (2 * Tensor_Dim01 - N2 - 1)) / 2]
                                 [N4 + (N3 * (2 * Tensor_Dim23 - N3 - 1)) / 2])
               : (N3 > N4
                    ? *data[N2 + (N1 * (2 * Tensor_Dim01 - N1 - 1)) / 2]
                           [N3 + (N4 * (2 * Tensor_Dim23 - N4 - 1)) / 2]
                    : *data[N2 + (N1 * (2 * Tensor_Dim01 - N1 - 1)) / 2]
                           [N4 + (N3 * (2 * Tensor_Dim23 - N3 - 1)) / 2]);
    }

    T *restrict &ptr(const int N1, const int N2, const int N3, const int N4) const
    {
#ifdef FTENSOR_DEBUG
      if(N1 >= Tensor_Dim01 || N1 < 0 || N2 >= Tensor_Dim01 || N2 < 0
         || N3 >= Tensor_Dim23 || N3 < 0 || N4 >= Tensor_Dim23 || N4 < 0)
        {
          std::stringstream s;
          s << "Bad index in Dg<T," << Tensor_Dim01 << "," << Tensor_Dim23
            << ">.ptr(" << N1 << "," << N2 << "," << N3 << "," << N4 << ")"
            << std::endl;
          throw std::out_of_range(s.str());
        }
#endif
      return N1 > N2
               ? (N3 > N4 ? data[N1 + (N2 * (2 * Tensor_Dim01 - N2 - 1)) / 2]
                                [N3 + (N4 * (2 * Tensor_Dim23 - N4 - 1)) / 2]
                          : data[N1 + (N2 * (2 * Tensor_Dim01 - N2 - 1)) / 2]
                                [N4 + (N3 * (2 * Tensor_Dim23 - N3 - 1)) / 2])
               : (N3 > N4 ? data[N2 + (N1 * (2 * Tensor_Dim01 - N1 - 1)) / 2]
                                [N3 + (N4 * (2 * Tensor_Dim23 - N4 - 1)) / 2]
                          : data[N2 + (N1 * (2 * Tensor_Dim01 - N1 - 1)) / 2]
                                [N4 + (N3 * (2 * Tensor_Dim23 - N3 - 1)) / 2]);
    }

    /* These operator()'s are the first part in constructing template
       expressions.  They can be used to slice off lower dimensional
       parts. They are not entirely safe, since you can accidently use a
       higher dimension than what is really allowed (like Dim=5). */

    template <char i, char j, char k, char l, int Dim01, int Dim23>
    Ddg_Expr<Ddg<T *, Tensor_Dim01, Tensor_Dim23>, T, Dim01, Dim23, i, j, k, l>
    operator()(const Index<i, Dim01> index1, const Index<j, Dim01> index2,
               const Index<k, Dim23> index3, const Index<l, Dim23> index4)
    {
      return Ddg_Expr<Ddg<T *, Tensor_Dim01, Tensor_Dim23>, T, Dim01, Dim23, i,
                      j, k, l>(*this);
    }

    template <char i, char j, char k, char l, int Dim01, int Dim23>
    Ddg_Expr<const Ddg<T *, Tensor_Dim01, Tensor_Dim23>, T, Dim01, Dim23, i, j,
             k, l>
    operator()(const Index<i, Dim01> index1, const Index<j, Dim01> index2,
               const Index<k, Dim23> index3,
               const Index<l, Dim23> index4) const
    {
      return Ddg_Expr<const Ddg<T *, Tensor_Dim01, Tensor_Dim23>, T, Dim01,
                      Dim23, i, j, k, l>(*this);
    }

    /* This is for expressions where a number is used for two slots, and
       an index for the other two, yielding a Tensor2_symmetric_Expr. */

    template <char i, char j, int N0, int N1, int Dim>
    Tensor2_symmetric_Expr<
      Ddg_number_01<const Ddg<T *, Tensor_Dim01, Tensor_Dim23>, T, N0, N1>, T,
      Dim, i, j>
    operator()(const Number<N0> n1, const Number<N1> n2,
               const Index<i, Dim> index1, const Index<j, Dim> index2) const
    {
      using TensorExpr
        = Ddg_number_01<const Ddg<T *, Tensor_Dim01, Tensor_Dim23>, T, N0, N1>;
      return Tensor2_symmetric_Expr<TensorExpr, T, Dim, i, j>(
        TensorExpr(*this));
    }

    template <char i, char j, int N0, int N1, int Dim>
    Tensor2_symmetric_Expr<
      Ddg_number_rhs_01<Ddg<T *, Tensor_Dim01, Tensor_Dim23>, T, N0, N1>, T,
      Dim, i, j>
    operator()(const Number<N0> n1, const Number<N1> n2,
               const Index<i, Dim> index1, const Index<j, Dim> index2)
    {
      using TensorExpr
        = Ddg_number_rhs_01<Ddg<T *, Tensor_Dim01, Tensor_Dim23>, T, N0, N1>;
      return Tensor2_symmetric_Expr<TensorExpr, T, Dim, i, j>(*this);
    }

    /* This is for expressions where a number is used for one slot, and
       an index for the other three, yielding a Dg_Expr. */

    template <char i, char j, char k, int N0, int Dim1, int Dim23>
    Dg_Expr<Ddg_number_0<const Ddg<T *, Tensor_Dim01, Tensor_Dim23>, T, N0>, T,
            Dim23, Dim1, i, j, k>
    operator()(const Number<N0> n1, const Index<k, Dim1> index3,
               const Index<i, Dim23> index1,
               const Index<j, Dim23> index2) const
    {
      using TensorExpr
        = Ddg_number_0<const Ddg<T *, Tensor_Dim01, Tensor_Dim23>, T, N0>;
      return Dg_Expr<TensorExpr, T, Dim23, Dim1, i, j, k>(TensorExpr(*this));
    }

    template <char i, char j, char k, int N0, int Dim1, int Dim23>
    Dg_Expr<Ddg_number_rhs_0<Ddg<T *, Tensor_Dim01, Tensor_Dim23>, T, N0>, T,
            Dim23, Dim1, i, j, k>
    operator()(const Number<N0> n1, const Index<k, Dim1> index3,
               const Index<i, Dim23> index1, const Index<j, Dim23> index2)
    {
      using TensorExpr
        = Ddg_number_rhs_0<Ddg<T *, Tensor_Dim01, Tensor_Dim23>, T, N0>;
      return Dg_Expr<TensorExpr, T, Dim23, Dim1, i, j, k>(*this);
    }

    /* This is for expressions where an int (not a Number) is used for
       two slots, and an index for the other two, yielding a
       Tensor2_symmetric_Expr. */

    template <char i, char j, int Dim>
    Tensor2_symmetric_Expr<
      Ddg_numeral_01<const Ddg<T *, Tensor_Dim01, Tensor_Dim23>, T>, T, Dim, i,
      j>
    operator()(const int N0, const int N1, const Index<i, Dim> index1,
               const Index<j, Dim> index2) const
    {
      using TensorExpr
        = Ddg_numeral_01<const Ddg<T *, Tensor_Dim01, Tensor_Dim23>, T>;
      return Tensor2_symmetric_Expr<TensorExpr, T, Dim, i, j>(
        TensorExpr(*this, N0, N1));
    }

    template <char i, char j, int Dim>
    Tensor2_symmetric_Expr<
      Ddg_numeral_23<Ddg<T *, Tensor_Dim01, Tensor_Dim23>, T>, T, Dim, i, j>
    operator()(const Index<i, Dim> index1, const Index<j, Dim> index2,
               const int N2, const int N3) const
    {
      using TensorExpr
        = Ddg_numeral_23<Ddg<T *, Tensor_Dim01, Tensor_Dim23>, T>;
      return Tensor2_symmetric_Expr<TensorExpr, T, Dim, i, j>(
        TensorExpr(*this, N2, N3));
    }

    /* int in three slots, Index in the other yielding a Tensor1_Expr. */

    template <char i, int Dim>
    Tensor1_Expr<Ddg_numeral_123<Ddg<T *, Tensor_Dim01, Tensor_Dim23>, T>, T,
                 Dim, i>
    operator()(const Index<i, Dim> index1, const int N1, const int N2,
               const int N3)
    {
      using TensorExpr
        = Ddg_numeral_123<Ddg<T *, Tensor_Dim01, Tensor_Dim23>, T>;
      return Tensor1_Expr<TensorExpr, T, Dim, i>(
        TensorExpr(*this, N1, N2, N3));
    }

    template <char i, int Dim>
    Tensor1_Expr<Ddg_numeral_123<Ddg<T *, Tensor_Dim01, Tensor_Dim23>, T>, T,
                 Dim, i>
    operator()(const int N1, const Index<i, Dim> index1, const int N2,
               const int N3)
    {
      using TensorExpr
        = Ddg_numeral_123<Ddg<T *, Tensor_Dim01, Tensor_Dim23>, T>;
      return Tensor1_Expr<TensorExpr, T, Dim, i>(
        TensorExpr(*this, N1, N2, N3));
    }

    /* This is for expressions where an int (not a Number) is used for
       one slot, and an index for the other three, yielding a
       Dg_Expr. */

    template <char i, char j, char k, int Dim1, int Dim23>
    Dg_Expr<Ddg_numeral_0<const Ddg<T *, Tensor_Dim01, Tensor_Dim23>, T>, T,
            Dim23, Dim1, i, j, k>
    operator()(const int N0, const Index<k, Dim1> index3,
               const Index<i, Dim23> index1,
               const Index<j, Dim23> index2) const
    {
      using TensorExpr
        = Ddg_numeral_0<const Ddg<T *, Tensor_Dim01, Tensor_Dim23>, T>;
      return Dg_Expr<TensorExpr, T, Dim23, Dim1, i, j, k>(
        TensorExpr(*this, N0));
    }

    /* The ++ operator increments the pointer, not the number that the
       pointer points to.  This allows iterating over a grid. */

    const Ddg &operator++() const {
      for(int i = 0; i < (Tensor_Dim01 * (Tensor_Dim01 + 1)) / 2; ++i)
        for(int j = 0; j < (Tensor_Dim23 * (Tensor_Dim23 + 1)) / 2; ++j)
          data[i][j] += inc;
      return *this;
    }

  private:
    template <int I>
    Ddg(const Ddg<PackPtr<T *, I>, Tensor_Dim01, Tensor_Dim23> &) = delete;
  };

  template <class T, int Tensor_Dim01, int Tensor_Dim23, int I>
  class Ddg<PackPtr<T *, I>, Tensor_Dim01, Tensor_Dim23>
      : public Ddg<T *, Tensor_Dim01, Tensor_Dim23> {

  public:
    /* Initializations for varying numbers of elements. */
    template <class... U>
    Ddg(U *... d) : Ddg<T *, Tensor_Dim01, Tensor_Dim23>(d...) {}

    /* The ++ operator increments the pointer, not the number that the
       pointer points to.  This allows iterating over a grid. */

    const Ddg &operator++() const {
      for (int i = 0; i < (Tensor_Dim01 * (Tensor_Dim01 + 1)) / 2; ++i)
        for (int j = 0; j < (Tensor_Dim23 * (Tensor_Dim23 + 1)) / 2; ++j)
          Ddg<T *, Tensor_Dim23, Tensor_Dim23>::data[i][j] += I;
      return *this;
    }
  };
}
