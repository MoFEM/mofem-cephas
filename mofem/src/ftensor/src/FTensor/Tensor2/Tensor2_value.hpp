#pragma once

/* A general version, not for pointers. */

#include <ostream>

#ifdef FTENSOR_DEBUG
#include <sstream>
#include <stdexcept>
#endif

#pragma once

namespace FTensor
{
  template <class T, int Tensor_Dim0, int Tensor_Dim1, Layout layout>
  class Tensor2
  {
    T data[(layout == column_major) ? Tensor_Dim0 : Tensor_Dim1]
          [(layout == column_major) ? Tensor_Dim1 : Tensor_Dim0];

  public:
    /* Initializations for varying numbers of elements, with each one
       defined for a particular Tensor_Dim.  To initialize a different
       dimension, just add the appropriate constructor and call to
       the Tensor2_constructor constructor. */
    Tensor2(T d00, T d01)
    {
      Tensor2_constructor<T, Tensor_Dim0, Tensor_Dim1, layout>(data, d00, d01);
    }
    Tensor2(T d00, T d01, T d10, T d11)
    {
      Tensor2_constructor<T, Tensor_Dim0, Tensor_Dim1, layout>(data, d00, d01,
                                                               d10, d11);
    }
    Tensor2(T d00, T d01, T d02, T d10, T d11, T d12, T d20, T d21, T d22)
    {
      Tensor2_constructor<T, Tensor_Dim0, Tensor_Dim1, layout>(
        data, d00, d01, d02, d10, d11, d12, d20, d21, d22);
    }
    Tensor2(T d00, T d01, T d10, T d11, T d20, T d21)
    {
      Tensor2_constructor<T, Tensor_Dim0, Tensor_Dim1, layout>(
        data, d00, d01, d10, d11, d20, d21);
    }
    Tensor2(T d00, T d01, T d02, T d03, T d10, T d11, T d12, T d13, T d20,
            T d21, T d22, T d23, T d30, T d31, T d32, T d33)
    {
      Tensor2_constructor<T, Tensor_Dim0, Tensor_Dim1, layout>(
        data, d00, d01, d02, d03, d10, d11, d12, d13, d20, d21, d22, d23, d30,
        d31, d32, d33);
    }
    Tensor2() {}

    /* There are two operator(int,int)'s, one for non-consts that lets you
       change the value, and one for consts that doesn't. */

    T &operator()(const int N1, const int N2)
    {
#ifdef FTENSOR_DEBUG
      if(N1 >= Tensor_Dim0 || N1 < 0 || N2 >= Tensor_Dim1 || N2 < 0)
        {
          std::stringstream s;
          s << "Bad index in Tensor2<T," << Tensor_Dim0 << "," << Tensor_Dim1
            << ">.operator(" << N1 << "," << N2 << ")" << std::endl;
          throw std::runtime_error(s.str());
        }
#endif
      return ((layout == column_major) ? data[N1][N2] : data[N2][N1]);
    }

    T operator()(const int N1, const int N2) const
    {
#ifdef FTENSOR_DEBUG
      if(N1 >= Tensor_Dim0 || N1 < 0 || N2 >= Tensor_Dim1 || N2 < 0)
        {
          std::stringstream s;
          s << "Bad index in Tensor2<T," << Tensor_Dim0 << "," << Tensor_Dim1
            << ">.operator(" << N1 << "," << N2 << ") const" << std::endl;
          throw std::runtime_error(s.str());
        }
#endif
      return ((layout == column_major) ? data[N1][N2] : data[N2][N1]);
    }

    /* These operator()'s are the first part in constructing template
       expressions.  They can be used to slice off lower dimensional
       parts. They are not entirely safe, since you can accidently use a
       higher dimension than what is really allowed (like Dim=5). */

    template <char i, char j, int Dim0, int Dim1>
    typename std::enable_if<
      (Tensor_Dim0 >= Dim0 && Tensor_Dim1 >= Dim1),
      Tensor2_Expr<Tensor2<T, Tensor_Dim0, Tensor_Dim1, layout>, T, Dim0, Dim1,
                   i, j>>::type
    operator()(const Index<i, Dim0>, const Index<j, Dim1>)
    {
      return Tensor2_Expr<Tensor2<T, Tensor_Dim0, Tensor_Dim1, layout>, T,
                          Dim0, Dim1, i, j>(*this);
    }

    template <char i, char j, int Dim0, int Dim1>
    typename std::enable_if<
      (Tensor_Dim0 >= Dim0 && Tensor_Dim1 >= Dim1),
      Tensor2_Expr<const Tensor2<T, Tensor_Dim0, Tensor_Dim1, layout>, T, Dim0,
                   Dim1, i, j>>::type
    operator()(const Index<i, Dim0>, const Index<j, Dim1>) const
    {
      return Tensor2_Expr<const Tensor2<T, Tensor_Dim0, Tensor_Dim1, layout>,
                          T, Dim0, Dim1, i, j>(*this);
    }

    /* This is for expressions where a number is used for one slot, and
       an index for another, yielding a Tensor1_Expr.  The non-const
       versions don't actually create a Tensor2_number_rhs_[01] object.
       They create a Tensor1_Expr directly, which provides the
       appropriate indexing operators.  The const versions do create a
       Tensor2_number_[01]. */

    template <char i, int Dim, int N>
    typename std::enable_if<
      (Tensor_Dim0 >= Dim && Tensor_Dim1 > N),
      Tensor1_Expr<Tensor2_number_rhs_1<
                     Tensor2<T, Tensor_Dim0, Tensor_Dim1, layout>, T, N>,
                   T, Dim, i>>::type
    operator()(const Index<i, Dim>, const Number<N>)
    {
      typedef Tensor2_number_rhs_1<
        Tensor2<T, Tensor_Dim0, Tensor_Dim1, layout>, T, N>
        TensorExpr;
      return Tensor1_Expr<TensorExpr, T, Dim, i>(*this);
    }

    template <char i, int Dim, int N>
    typename std::enable_if<
      (Tensor_Dim0 > N && Tensor_Dim1 >= Dim),
      Tensor1_Expr<Tensor2_number_rhs_0<
                     Tensor2<T, Tensor_Dim0, Tensor_Dim1, layout>, T, N>,
                   T, Dim, i>>::type
    operator()(const Number<N>, const Index<i, Dim>)
    {
      typedef Tensor2_number_rhs_0<
        Tensor2<T, Tensor_Dim0, Tensor_Dim1, layout>, T, N>
        TensorExpr;
      return Tensor1_Expr<TensorExpr, T, Dim, i>(*this);
    }

    template <char i, int Dim, int N>
    typename std::enable_if<
      (Tensor_Dim0 >= Dim && Tensor_Dim1 > N),
      Tensor1_Expr<const Tensor2_number_1<
                     const Tensor2<T, Tensor_Dim0, Tensor_Dim1, layout>, T, N>,
                   T, Dim, i>>::type
    operator()(const Index<i, Dim>, const Number<N>) const
    {
      typedef const Tensor2_number_1<
        const Tensor2<T, Tensor_Dim0, Tensor_Dim1, layout>, T, N>
        TensorExpr;
      return Tensor1_Expr<TensorExpr, T, Dim, i>(TensorExpr(*this));
    }

    template <char i, int Dim, int N>
    typename std::enable_if<
      (Tensor_Dim0 > N && Tensor_Dim1 >= Dim),
      Tensor1_Expr<const Tensor2_number_0<
                     const Tensor2<T, Tensor_Dim0, Tensor_Dim1, layout>, T, N>,
                   T, Dim, i>>::type
    operator()(const Number<N>, const Index<i, Dim>) const
    {
      typedef const Tensor2_number_0<
        const Tensor2<T, Tensor_Dim0, Tensor_Dim1, layout>, T, N>
        TensorExpr;
      return Tensor1_Expr<TensorExpr, T, Dim, i>(TensorExpr(*this));
    }

    /* This is for expressions where an actual number (not a Number<>)
       is used for one slot, and an index for another, yielding a
       Tensor1_Expr. */

    /* Unfortunately since this integers can only be checked at run time
       i can only partially protect this expressions. We should suggest
       when ever possible to use Number<>. */

    template <char i, int Dim>
    typename std::enable_if<
      (Tensor_Dim0 >= Dim),
      Tensor1_Expr<
        const Tensor2_numeral_1<const Tensor2<T, Tensor_Dim0, Tensor_Dim1>, T>,
        T, Dim, i>>::type
    operator()(const Index<i, Dim>, const int N) const
    {
      typedef const Tensor2_numeral_1<
        const Tensor2<T, Tensor_Dim0, Tensor_Dim1, layout>, T>
        TensorExpr;
      return Tensor1_Expr<TensorExpr, T, Dim, i>(TensorExpr(*this, N));
    }

    template <char i, int Dim>
    typename std::enable_if<
      (Tensor_Dim1 >= Dim),
      Tensor1_Expr<
        const Tensor2_numeral_0<const Tensor2<T, Tensor_Dim0, Tensor_Dim1>, T>,
        T, Dim, i>>::type
    operator()(const int N, const Index<i, Dim>) const
    {
      typedef const Tensor2_numeral_0<
        const Tensor2<T, Tensor_Dim0, Tensor_Dim1, layout>, T>
        TensorExpr;
      return Tensor1_Expr<TensorExpr, T, Dim, i>(TensorExpr(*this, N));
    }

    /* These two operator()'s return the Tensor2 with internal
       contractions, yielding a T.  I have to specify one for both
       const and non-const because otherwise they compiler will use the
       operator() which gives a Tensor2_Expr<>. */

    template <char i, int Dim>
    typename std::enable_if<(Tensor_Dim0 >= Dim && Tensor_Dim1 >= Dim), T>::type
    operator()(const Index<i, Dim>, const Index<i, Dim>)
    {
      return internal_contract(Number<Dim>());
    }

    template <char i, int Dim>
    typename std::enable_if<(Tensor_Dim0 >= Dim && Tensor_Dim1 >= Dim), T>::type
    operator()(const Index<i, Dim>, const Index<i, Dim>) const
    {
      return internal_contract(Number<Dim>());
    }

  private:
    template <int N> T internal_contract(Number<N>) const
    {
      return data[N - 1][N - 1] + internal_contract(Number<N - 1>());
    }

    T internal_contract(Number<1>) const { return data[0][0]; }
  };
}

/// JSON compatible output

namespace FTensor
{
  template <class T, int Tensor_Dim0, int Tensor_Dim1, FTensor::Layout layout>
  std::ostream &Tensor2_ostream_row(
    std::ostream &os,
    const FTensor::Tensor2<T, Tensor_Dim0, Tensor_Dim1, layout> &t,
    const int &i)
  {
    os << '[';
    for(int j = 0; j + 1 < Tensor_Dim1; ++j)
      {
        os << t(i, j) << ',';
      }
    if(Tensor_Dim1 > 0)
      {
        os << t(i, Tensor_Dim1 - 1);
      }
    os << ']';
    return os;
  }
}

template <class T, int Tensor_Dim0, int Tensor_Dim1, FTensor::Layout layout>
std::ostream &
operator<<(std::ostream &os,
           const FTensor::Tensor2<T, Tensor_Dim0, Tensor_Dim1, layout> &t)
{
  os << '[';
  for(int i = 0; i + 1 < Tensor_Dim0; ++i)
    {
      FTensor::Tensor2_ostream_row(os, t, i);
      os << ',';
    }
  if(Tensor_Dim0 > 0)
    {
      FTensor::Tensor2_ostream_row(os, t, Tensor_Dim0 - 1);
    }
  os << ']';
  return os;
}

namespace FTensor
{
  template <class T, int Tensor_Dim0, int Tensor_Dim1, FTensor::Layout layout>
  std::istream &
  Tensor2_istream_row(std::istream &is,
                      FTensor::Tensor2<T, Tensor_Dim0, Tensor_Dim1, layout> &t,
                      const int &i)
  {
    char c;
    is >> c;
    for(int j = 0; j + 1 < Tensor_Dim1; ++j)
      {
        is >> t(i, j) >> c;
      }
    if(Tensor_Dim1 > 0)
      {
        is >> t(i, Tensor_Dim1 - 1);
      }
    is >> c;
    return is;
  }
}

template <class T, int Tensor_Dim0, int Tensor_Dim1, FTensor::Layout layout>
std::istream &
operator>>(std::istream &is,
           FTensor::Tensor2<T, Tensor_Dim0, Tensor_Dim1, layout> &t)
{
  char c;
  is >> c;
  for(int i = 0; i + 1 < Tensor_Dim0; ++i)
    {
      FTensor::Tensor2_istream_row(is, t, i);
      is >> c;
    }
  if(Tensor_Dim0 > 0)
    {
      FTensor::Tensor2_istream_row(is, t, Tensor_Dim0 - 1);
    }
  is >> c;
  return is;
}
