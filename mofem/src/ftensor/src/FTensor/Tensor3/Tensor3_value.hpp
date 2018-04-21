/* A general version, not for pointers. */

#pragma once

#include <ostream>

namespace FTensor
{
  template <class T, int Tensor_Dim0, int Tensor_Dim1, int Tensor_Dim2>
  class Tensor3
  {
    T data[Tensor_Dim1][Tensor_Dim1][Tensor_Dim2];

  public:
    Tensor3() {}

    /* Tensor_Dim0=2, Tensor_Dim1=2, Tensor_Dim2=2 */
    Tensor3(T d000, T d001, T d010, T d011, T d100, T d101, T d110, T d111)
    {
      Tensor3_constructor<T, Tensor_Dim0, Tensor_Dim1, Tensor_Dim2>(
        data, d000, d001, d010, d011, d100, d101, d110, d111);
    }

    /* Tensor_Dim0=3, Tensor_Dim1=3, Tensor_Dim2=3 */
    Tensor3(T d000, T d001, T d002, T d010, T d011, T d012, T d020, T d021,
            T d022, T d100, T d101, T d102, T d110, T d111, T d112, T d120,
            T d121, T d122, T d200, T d201, T d202, T d210, T d211, T d212,
            T d220, T d221, T d222)
    {
      Tensor3_constructor<T, Tensor_Dim0, Tensor_Dim1, Tensor_Dim2>(
        data, d000, d001, d002, d010, d011, d012, d020, d021, d022, d100, d101,
        d102, d110, d111, d112, d120, d121, d122, d200, d201, d202, d210, d211,
        d212, d220, d221, d222);
    }

    /* Tensor_Dim0=4, Tensor_Dim1=4, Tensor_Dim2=4 */
    Tensor3(T d000, T d001, T d002, T d003, T d010, T d011, T d012, T d013,
            T d020, T d021, T d022, T d023, T d030, T d031, T d032, T d033,

            T d100, T d101, T d102, T d103, T d110, T d111, T d112, T d113,
            T d120, T d121, T d122, T d123, T d130, T d131, T d132, T d133,

            T d200, T d201, T d202, T d203, T d210, T d211, T d212, T d213,
            T d220, T d221, T d222, T d223, T d230, T d231, T d232, T d233,

            T d300, T d301, T d302, T d303, T d310, T d311, T d312, T d313,
            T d320, T d321, T d322, T d323, T d330, T d331, T d332, T d333)
    {
      Tensor3_constructor<T, Tensor_Dim0, Tensor_Dim1, Tensor_Dim2>(
        data, d000, d001, d002, d003, d010, d011, d012, d013, d020, d021, d022,
        d023, d030, d031, d032, d033,

        d100, d101, d102, d103, d110, d111, d112, d113, d120, d121, d122, d123,
        d130, d131, d132, d133,

        d200, d201, d202, d203, d210, d211, d212, d213, d220, d221, d222, d223,
        d230, d231, d232, d233,

        d300, d301, d302, d303, d310, d311, d312, d313, d320, d321, d322, d323,
        d330, d331, d332, d333);
    }

    /* There are two operator(int,int,int)'s, one for non-consts that lets you
       change the value, and one for consts that doesn't. */

    T &operator()(const int N1, const int N2, const int N3)
    {
#ifdef FTENSOR_DEBUG
      if(N1 >= Tensor_Dim0 || N1 < 0 || N2 >= Tensor_Dim1 || N2 < 0
         || N3 >= Tensor_Dim2 || N3 < 0)
        {
          std::stringstream s;
          s << "Bad index in Tensor3<T," << Tensor_Dim0 << "," << Tensor_Dim1
            << "," << Tensor_Dim2 << ">.operator(" << N1 << "," << N2 << ","
            << N3 << ")" << std::endl;
          throw std::runtime_error(s.str());
        }
#endif
      return data[N1][N2][N3];
    }

    T operator()(const int N1, const int N2, const int N3) const
    {
#ifdef FTENSOR_DEBUG
      if(N1 >= Tensor_Dim0 || N1 < 0 || N2 >= Tensor_Dim1 || N2 < 0
         || N3 >= Tensor_Dim2 || N3 < 0)
        {
          std::stringstream s;
          s << "Bad index in Tensor3<T," << Tensor_Dim0 << "," << Tensor_Dim1
            << "," << Tensor_Dim2 << ">.operator(" << N1 << "," << N2 << ","
            << N3 << ") const" << std::endl;
          throw std::runtime_error(s.str());
        }
#endif
      return data[N1][N2][N3];
    }

    /* These operator()'s are the first part in constructing template
       expressions. */

    template <char i, char j, char k, int Dim0, int Dim1, int Dim2>
    typename std::enable_if<
      (Tensor_Dim0 >= Dim0 && Tensor_Dim1 >= Dim1 && Tensor_Dim2 >= Dim2),
      Tensor3_Expr<Tensor3<T, Tensor_Dim0, Tensor_Dim1, Tensor_Dim2>, T, Dim0,
                   Dim1, Dim2, i, j, k>>::type
    operator()(const Index<i, Dim0>, const Index<j, Dim1>,
               const Index<k, Dim2>)
    {
      return Tensor3_Expr<Tensor3<T, Tensor_Dim0, Tensor_Dim1, Tensor_Dim2>, T,
                          Dim0, Dim1, Dim2, i, j, k>(*this);
    }

    template <char i, char j, char k, int Dim0, int Dim1, int Dim2>
    typename std::enable_if<
      (Tensor_Dim0 >= Dim0 && Tensor_Dim1 >= Dim1 && Tensor_Dim2 >= Dim2),
      Tensor3_Expr<const Tensor3<T, Tensor_Dim0, Tensor_Dim1, Tensor_Dim2>, T,
                   Dim0, Dim1, Dim2, i, j, k>>::type
    operator()(const Index<i, Dim0>, const Index<j, Dim1>,
               const Index<k, Dim2>) const
    {
      return Tensor3_Expr<
        const Tensor3<T, Tensor_Dim0, Tensor_Dim1, Tensor_Dim2>, T, Dim0, Dim1,
        Dim2, i, j, k>(*this);
    }

    /* These operators are for internal contractions. We can not do
       something like A(i,j,j) where i and j have different dimensions,
       because it becomes ambiguous. */
    /*TODO I dont know how i and j having different dimensions can be
     * ambigous?? Since you do Use that in A(j,j,i)*/
    /* A(i,j,j) */

    template <char i, char j, int Dim>
    inline typename std::enable_if<
      (Tensor_Dim0 >= Dim && Tensor_Dim1 >= Dim && Tensor_Dim2 >= Dim),
      Tensor1_Expr<const Tensor3_contracted_12<
                     Tensor3<T, Tensor_Dim0, Tensor_Dim1, Tensor_Dim2>, T, Dim>,
                   T, Dim, i>>::type
    operator()(const Index<i, Dim>, const Index<j, Dim>,
               const Index<j, Dim>) const
    {
      typedef const Tensor3_contracted_12<
        Tensor3<T, Tensor_Dim0, Tensor_Dim1, Tensor_Dim2>, T, Dim>
        TensorExpr;
      return Tensor1_Expr<TensorExpr, T, Dim, i>(TensorExpr(*this));
    }

    /* A(j,i,j) */

    template <char i, char j, int Dim>
    inline typename std::enable_if<
      (Tensor_Dim0 >= Dim && Tensor_Dim1 >= Dim && Tensor_Dim2 >= Dim),
      Tensor1_Expr<const Tensor3_contracted_02<
                     Tensor3<T, Tensor_Dim0, Tensor_Dim1, Tensor_Dim2>, T, Dim>,
                   T, Dim, i>>::type
    operator()(const Index<j, Dim>, const Index<i, Dim>,
               const Index<j, Dim>) const
    {
      typedef const Tensor3_contracted_02<
        Tensor3<T, Tensor_Dim0, Tensor_Dim1, Tensor_Dim2>, T, Dim>
        TensorExpr;
      return Tensor1_Expr<TensorExpr, T, Dim, i>(TensorExpr(*this));
    }

    /* A(j,j,i) */

    template <char i, char j, int Dim, int Dim2>
    inline typename std::enable_if<
      (Tensor_Dim0 >= Dim && Tensor_Dim1 >= Dim && Tensor_Dim2 >= Dim2),
      Tensor1_Expr<const Tensor3_contracted_01<
                     Tensor3<T, Tensor_Dim0, Tensor_Dim1, Tensor_Dim2>, T, Dim>,
                   T, Dim2, i>>::type
    operator()(const Index<j, Dim>, const Index<j, Dim>,
               const Index<i, Dim2>) const
    {
      typedef const Tensor3_contracted_01<
        Tensor3<T, Tensor_Dim0, Tensor_Dim1, Tensor_Dim2>, T, Dim2>
        TensorExpr;
      return Tensor1_Expr<TensorExpr, T, Dim2, i>(TensorExpr(*this));
    }

    /* This is for expressions where a number is used for one slot, and
       indices for the others, yielding a Tensor2_Expr or
       Tensor2_symmetric_Expr.  The non-const versions don't actually
       create a Tensor3_number_rhs_* object, while the const versions
       do create a Tensor3_number_*. */

    /* First slot. */

    template <char i, char j, int N, int Dim1, int Dim2>
    typename std::enable_if<
      (Tensor_Dim0 > N && Tensor_Dim1 >= Dim1 && Tensor_Dim2 >= Dim2),
      Tensor2_Expr<Tensor3_number_rhs_0<
                     Tensor3<T, Tensor_Dim0, Tensor_Dim1, Tensor_Dim2>, T, N>,
                   T, Dim1, Dim2, i, j>>::type
    operator()(const Number<N>, const Index<i, Dim1>, const Index<j, Dim2>)
    {
      typedef Tensor3_number_rhs_0<
        Tensor3<T, Tensor_Dim0, Tensor_Dim1, Tensor_Dim2>, T, N>
        TensorExpr;
      return Tensor2_Expr<TensorExpr, T, Dim1, Dim2, i, j>(*this);
    }

    template <char i, char j, int N, int Dim1, int Dim2>
    typename std::enable_if<
      (Tensor_Dim0 > N && Tensor_Dim1 >= Dim1 && Tensor_Dim2 >= Dim2),
      const Tensor2_Expr<
        const Tensor3_number_0<
          const Tensor3<T, Tensor_Dim0, Tensor_Dim1, Tensor_Dim2>, T, N>,
        T, Dim1, Dim2, i, j>>::type
    operator()(const Number<N>, const Index<i, Dim1>,
               const Index<j, Dim2>) const
    {
      typedef const Tensor3_number_0<
        const Tensor3<T, Tensor_Dim0, Tensor_Dim1, Tensor_Dim2>, T, N>
        TensorExpr;
      return Tensor2_Expr<TensorExpr, T, Dim1, Dim2, i, j>(TensorExpr(*this));
    }

    /* Second slot. */

    template <char i, char j, int N, int Dim0, int Dim2>
    typename std::enable_if<
      (Tensor_Dim0 >= Dim0 && Tensor_Dim1 > N && Tensor_Dim2 >= Dim2),
      Tensor2_Expr<Tensor3_number_rhs_0<
                     Tensor3<T, Tensor_Dim0, Tensor_Dim1, Tensor_Dim2>, T, N>,
                   T, Dim0, Dim2, i, j>>::type
    operator()(const Index<i, Dim0>, const Number<N>, const Index<j, Dim2>)
    {
      typedef Tensor3_number_rhs_0<
        Tensor3<T, Tensor_Dim0, Tensor_Dim1, Tensor_Dim2>, T, N>
        TensorExpr;
      return Tensor2_Expr<TensorExpr, T, Dim0, Dim2, i, j>(*this);
    }

    template <char i, char j, int N, int Dim0, int Dim2>
    typename std::enable_if<
      (Tensor_Dim0 >= Dim0 && Tensor_Dim1 > N && Tensor_Dim2 >= Dim2),
      const Tensor2_Expr<
        const Tensor3_number_0<
          const Tensor3<T, Tensor_Dim0, Tensor_Dim1, Tensor_Dim2>, T, N>,
        T, Dim0, Dim2, i, j>>::type
    operator()(const Index<i, Dim0>, const Number<N>,
               const Index<j, Dim2>) const
    {
      typedef const Tensor3_number_0<
        const Tensor3<T, Tensor_Dim0, Tensor_Dim1, Tensor_Dim2>, T, N>
        TensorExpr;
      return Tensor2_Expr<TensorExpr, T, Dim0, Dim2, i, j>(TensorExpr(*this));
    }

    /* Third slot. */
    // TODO allow two different dimensions here
    template <char i, char j, int N, int Dim>
    typename std::enable_if<
      (Tensor_Dim0 >= Dim && Tensor_Dim1 >= Dim && Tensor_Dim2 > N),
      Tensor2_symmetric_Expr<
        Tensor3_number_rhs_2<Tensor3<T, Tensor_Dim0, Tensor_Dim1, Tensor_Dim2>,
                             T, N>,
        T, Dim, i, j>>::type
    operator()(const Index<i, Dim>, const Index<j, Dim>, const Number<N>)
    {
      typedef Tensor3_number_rhs_2<
        Tensor3<T, Tensor_Dim0, Tensor_Dim1, Tensor_Dim2>, T, N>
        TensorExpr;
      return Tensor2_symmetric_Expr<TensorExpr, T, Dim, i, j>(*this);
    }
    // TODO allow two different dimensions here
    template <char i, char j, int N, int Dim>
    typename std::enable_if<
      (Tensor_Dim0 >= Dim && Tensor_Dim1 >= Dim && Tensor_Dim2 > N),
      Tensor2_symmetric_Expr<
        const Tensor3_number_2<
          const Tensor3<T, Tensor_Dim0, Tensor_Dim1, Tensor_Dim2>, T, N>,
        T, Dim, i, j>>::type
    operator()(const Index<i, Dim>, const Index<j, Dim>, const Number<N>) const
    {
      typedef const Tensor3_number_2<
        const Tensor3<T, Tensor_Dim0, Tensor_Dim1, Tensor_Dim2>, T, N>
        TensorExpr;
      return Tensor2_symmetric_Expr<TensorExpr, T, Dim, i, j>(
        TensorExpr(*this));
    }

    /* This is for expressions where a number is used for two slots, and
       an Index for the other, yielding a Tensor1_Expr.  The non-const
       versions don't actually create a Tensor3_number_rhs_* object,
       while the const versions do create a Tensor3_number_*. */

    /* Index in first slot. */

    template <char i, int N1, int N2, int Dim>
    typename std::enable_if<
      (Tensor_Dim0 >= Dim && Tensor_Dim1 > N1 && Tensor_Dim2 > N2),
      Tensor1_Expr<
        Tensor3_number_rhs_12<
          Tensor3<T, Tensor_Dim0, Tensor_Dim1, Tensor_Dim2>, T, N1, N2>,
        T, Dim, i>>::type
    operator()(const Index<i, Dim> index, const Number<N1>, const Number<N2>)
    {
      typedef Tensor3_number_rhs_12<
        Tensor3<T, Tensor_Dim0, Tensor_Dim1, Tensor_Dim2>, T, N1, N2>
        TensorExpr;
      return Tensor1_Expr<TensorExpr, T, Dim, i>(*this);
    }

    template <char i, int N1, int N2, int Dim>
    typename std::enable_if<
      (Tensor_Dim0 >= Dim && Tensor_Dim1 > N1 && Tensor_Dim2 > N2),
      const Tensor1_Expr<
        const Tensor3_number_12<
          const Tensor3<T, Tensor_Dim0, Tensor_Dim1, Tensor_Dim2>, T, N1, N2>,
        T, Dim, i>>::type
    operator()(const Index<i, Dim> index, const Number<N1>,
               const Number<N2>) const
    {
      typedef const Tensor3_number_12<
        const Tensor3<T, Tensor_Dim0, Tensor_Dim1, Tensor_Dim2>, T, N1, N2>
        TensorExpr;
      return Tensor1_Expr<TensorExpr, T, Dim, i>(TensorExpr(*this));
    }

    /* Index in second slot.  I use the same structures as for the Index
       in the first slot since the tensor is symmetric on the first two
       indices. */

    template <char i, int N0, int N2, int Dim>
    typename std::enable_if<
      (Tensor_Dim0 > N0 && Tensor_Dim1 >= Dim && Tensor_Dim2 > N2),
      Tensor1_Expr<
        Tensor3_number_rhs_12<
          Tensor3<T, Tensor_Dim0, Tensor_Dim1, Tensor_Dim2>, T, N0, N2>,
        T, Dim, i>>::type
    operator()(const Number<N0>, const Index<i, Dim> index, const Number<N2>)
    {
      typedef Tensor3_number_rhs_12<
        Tensor3<T, Tensor_Dim0, Tensor_Dim1, Tensor_Dim2>, T, N0, N2>
        TensorExpr;
      return Tensor1_Expr<TensorExpr, T, Dim, i>(*this);
    }

    template <char i, int N0, int N2, int Dim>
    typename std::enable_if<
      (Tensor_Dim0 > N0 && Tensor_Dim1 >= Dim && Tensor_Dim2 > N2),
      Tensor1_Expr<
        const Tensor3_number_12<
          const Tensor3<T, Tensor_Dim0, Tensor_Dim1, Tensor_Dim2>, T, N0, N2>,
        T, Dim, i>>::type
    operator()(const Number<N0>, const Index<i, Dim> index,
               const Number<N2>) const
    {
      typedef const Tensor3_number_12<
        const Tensor3<T, Tensor_Dim0, Tensor_Dim1, Tensor_Dim2>, T, N0, N2>
        TensorExpr;
      return Tensor1_Expr<TensorExpr, T, Dim, i>(TensorExpr(*this));
    }

    /* Index in third slot. */

    template <char i, int N0, int N1, int Dim>
    typename std::enable_if<
      (Tensor_Dim0 > N0 && Tensor_Dim1 > N1 && Tensor_Dim2 >= Dim),
      Tensor1_Expr<
        Tensor3_number_rhs_01<
          Tensor3<T, Tensor_Dim0, Tensor_Dim1, Tensor_Dim2>, T, N0, N1>,
        T, Dim, i>>::type
    operator()(const Number<N0>, const Number<N1>, const Index<i, Dim> index)
    {
      typedef Tensor3_number_rhs_01<
        Tensor3<T, Tensor_Dim0, Tensor_Dim1, Tensor_Dim2>, T, N0, N1>
        TensorExpr;
      return Tensor1_Expr<TensorExpr, T, Dim, i>(*this);
    }

    template <char i, int N0, int N1, int Dim>
    typename std::enable_if<
      (Tensor_Dim0 > N0 && Tensor_Dim1 > N1 && Tensor_Dim2 >= Dim),
      const Tensor1_Expr<
        const Tensor3_number_01<
          const Tensor3<T, Tensor_Dim0, Tensor_Dim1, Tensor_Dim2>, T, N0, N1>,
        T, Dim, i>>::type
    operator()(const Number<N0>, const Number<N1>,
               const Index<i, Dim> index) const
    {
      typedef const Tensor3_number_01<
        const Tensor3<T, Tensor_Dim0, Tensor_Dim1, Tensor_Dim2>, T, N0, N1>
        TensorExpr;
      return Tensor1_Expr<TensorExpr, T, Dim, i>(TensorExpr(*this));
    }

    /* Specializations for using actual numbers instead of Number<>.
       This is for expressions where an actual number is used for one
       slot, and indices for the others, yielding a Tensor2_Expr or
       Tensor2_symmetric_Expr. I only define the const versions. */

    /* First slot. */

    template <char i, char j, int Dim1, int Dim2>
    typename std::enable_if<
      (Tensor_Dim1 >= Dim1 && Tensor_Dim2 >= Dim2),
      Tensor2_Expr<const Tensor3_numeral_0<
                     const Tensor3<T, Tensor_Dim0, Tensor_Dim1, Tensor_Dim2>, T>,
                   T, Dim1, Dim2, i, j>>::type
    operator()(const int N, const Index<i, Dim1>, const Index<j, Dim2>) const
    {
      typedef const Tensor3_numeral_0<
        const Tensor3<T, Tensor_Dim0, Tensor_Dim1, Tensor_Dim2>, T>
        TensorExpr;
      return Tensor2_Expr<TensorExpr, T, Dim1, Dim2, i, j>(
        TensorExpr(*this, N));
    }

    /* Second slot. */

    template <char i, char j, int Dim0, int Dim2>
    typename std::enable_if<
      (Tensor_Dim0 >= Dim0 && Tensor_Dim2 >= Dim2),
      Tensor2_Expr<const Tensor3_numeral_0<
                     const Tensor3<T, Tensor_Dim0, Tensor_Dim1, Tensor_Dim2>, T>,
                   T, Dim0, Dim2, i, j>>::type
    operator()(const Index<i, Dim0>, const int N, const Index<j, Dim2>) const
    {
      typedef const Tensor3_numeral_0<
        const Tensor3<T, Tensor_Dim0, Tensor_Dim1, Tensor_Dim2>, T>
        TensorExpr;
      return Tensor2_Expr<TensorExpr, T, Dim0, Dim2, i, j>(
        TensorExpr(*this, N));
    }

    /* Third slot. */
    // TODO Allow multiple dimensions here
    template <char i, char j, int Dim>
    typename std::enable_if<
      (Tensor_Dim0 >= Dim && Tensor_Dim1 >= Dim),
      const Tensor2_symmetric_Expr<
        const Tensor3_numeral_2<
          const Tensor3<T, Tensor_Dim0, Tensor_Dim1, Tensor_Dim2>, T>,
        T, Dim, i, j>>::type
    operator()(const Index<i, Dim>, const Index<j, Dim>, const int N) const
    {
      typedef const Tensor3_numeral_2<
        const Tensor3<T, Tensor_Dim0, Tensor_Dim1, Tensor_Dim2>, T>
        TensorExpr;
      return Tensor2_symmetric_Expr<TensorExpr, T, Dim, i, j>(
        TensorExpr(*this, N));
    }

    /* This is for expressions where a numeral is used for two slots, and
       an Index for the other, yielding a Tensor1_Expr. */

    /* Index in first slot. */

    template <char i, int Dim>
    typename std::enable_if<
      (Tensor_Dim0 >= Dim),
      const Tensor1_Expr<
        const Tensor3_numeral_12<
          const Tensor3<T, Tensor_Dim0, Tensor_Dim1, Tensor_Dim2>, T>,
        T, Dim, i>>::type
    operator()(const Index<i, Dim> index, const int N1, const int N2) const
    {
      typedef const Tensor3_numeral_12<
        const Tensor3<T, Tensor_Dim0, Tensor_Dim1, Tensor_Dim2>, T>
        TensorExpr;
      return Tensor1_Expr<TensorExpr, T, Dim, i>(TensorExpr(*this, N1, N2));
    }

    /* Index in second slot.  I use the same structures as for the Index
       in the first slot since the tensor is symmetric on the first two
       indices. */

    template <char i, int Dim>
    typename std::enable_if<
      (Tensor_Dim1 >= Dim),
      const Tensor1_Expr<
        const Tensor3_numeral_12<
          const Tensor3<T, Tensor_Dim0, Tensor_Dim1, Tensor_Dim2>, T>,
        T, Dim, i>>::type
    operator()(const int N1, const Index<i, Dim> index, const int N2) const
    {
      typedef const Tensor3_numeral_12<
        const Tensor3<T, Tensor_Dim0, Tensor_Dim1, Tensor_Dim2>, T>
        TensorExpr;
      return Tensor1_Expr<TensorExpr, T, Dim, i>(TensorExpr(*this, N1, N2));
    }

    /* Index in third slot. */

    template <char i, int Dim>
    typename std::enable_if<
      (Tensor_Dim2 >= Dim),
      const Tensor1_Expr<
        const Tensor3_numeral_01<
          const Tensor3<T, Tensor_Dim0, Tensor_Dim1, Tensor_Dim2>, T>,
        T, Dim, i>>::type
    operator()(const int N1, const int N2, const Index<i, Dim> index) const
    {
      typedef const Tensor3_numeral_01<
        const Tensor3<T, Tensor_Dim0, Tensor_Dim1, Tensor_Dim2>, T>
        TensorExpr;
      return Tensor1_Expr<TensorExpr, T, Dim, i>(TensorExpr(*this, N1, N2));
    }
  };
}

/// JSON compatible output for numbers.  If T is something funky
/// (strings, arrays, etc.) then you are going to have a bad time.

namespace FTensor
{
  template <class T, int Tensor_Dim0, int Tensor_Dim1, int Tensor_Dim2>
  std::ostream &Tensor3_ostream_row(
    std::ostream &os,
    const FTensor::Tensor3<T, Tensor_Dim0, Tensor_Dim1, Tensor_Dim2> &t,
    const int &i, const int &j)
  {
    os << '[';
    for(int k = 0; k + 1 < Tensor_Dim2; ++k)
      {
        os << t(i, j, k) << ',';
      }
    if(Tensor_Dim2 > 0)
      {
        os << t(i, j, Tensor_Dim2 - 1);
      }
    os << ']';
    return os;
  }

  template <class T, int Tensor_Dim0, int Tensor_Dim1, int Tensor_Dim2>
  std::ostream &Tensor3_ostream_block(
    std::ostream &os,
    const FTensor::Tensor3<T, Tensor_Dim0, Tensor_Dim1, Tensor_Dim2> &t,
    const int &i)
  {
    os << '[';
    for(int j = 0; j + 1 < Tensor_Dim1; ++j)
      {
        FTensor::Tensor3_ostream_row(os, t, i, j);
        os << ',';
      }
    if(Tensor_Dim1 > 0)
      {
        FTensor::Tensor3_ostream_row(os, t, i, Tensor_Dim1 - 1);
      }
    os << ']';
    return os;
  }
}

template <class T, int Tensor_Dim0, int Tensor_Dim1, int Tensor_Dim2>
std::ostream &
operator<<(std::ostream &os,
           const FTensor::Tensor3<T, Tensor_Dim0, Tensor_Dim1, Tensor_Dim2> &t)
{
  os << '[';
  for(int i = 0; i + 1 < Tensor_Dim0; ++i)
    {
      FTensor::Tensor3_ostream_block(os, t, i);
      os << ',';
    }
  if(Tensor_Dim0 > 0)
    {
      FTensor::Tensor3_ostream_block(os, t, Tensor_Dim0 - 1);
    }
  os << ']';
  return os;
}

/// JSON compatible input.  It does not validate that separators are
/// actually braces '[' or commas ','.  It also ignores errors from
/// missing trailing characters.  So you could do something like
///
///   Tensor3<double,2,2,2> t3_1;
///   std::stringstream ss(":::3:4:::7:8:::::11:12:::13:14");
///   ss >> t3_1;

namespace FTensor
{
  template <class T, int Tensor_Dim0, int Tensor_Dim1, int Tensor_Dim2>
  std::istream &Tensor3_istream_row(
    std::istream &is,
    FTensor::Tensor3<T, Tensor_Dim0, Tensor_Dim1, Tensor_Dim2> &t,
    const int &i, const int &j)
  {
    char c;
    is >> c;
    for(int k = 0; k + 1 < Tensor_Dim2; ++k)
      {
        is >> t(i, j, k) >> c;
      }
    if(Tensor_Dim2 > 0)
      {
        is >> t(i, j, Tensor_Dim2 - 1);
      }
    is >> c;
    return is;
  }

  template <class T, int Tensor_Dim0, int Tensor_Dim1, int Tensor_Dim2>
  std::istream &Tensor3_istream_block(
    std::istream &is,
    FTensor::Tensor3<T, Tensor_Dim0, Tensor_Dim1, Tensor_Dim2> &t,
    const int &i)
  {
    char c;
    is >> c;
    for(int j = 0; j + 1 < Tensor_Dim1; ++j)
      {
        Tensor3_istream_row(is, t, i, j);
        is >> c;
      }
    if(Tensor_Dim1 > 0)
      {
        Tensor3_istream_row(is, t, i, Tensor_Dim1 - 1);
      }
    is >> c;
    return is;
  }
}

template <class T, int Tensor_Dim0, int Tensor_Dim1, int Tensor_Dim2>
std::istream &
operator>>(std::istream &is,
           FTensor::Tensor3<T, Tensor_Dim0, Tensor_Dim1, Tensor_Dim2> &t)
{
  char c;
  is >> c;
  for(int i = 0; i + 1 < Tensor_Dim0; ++i)
    {
      FTensor::Tensor3_istream_block(is, t, i);
      is >> c;
    }
  if(Tensor_Dim0 > 0)
    {
      FTensor::Tensor3_istream_block(is, t, Tensor_Dim0 - 1);
    }
  is >> c;
  return is;
}
