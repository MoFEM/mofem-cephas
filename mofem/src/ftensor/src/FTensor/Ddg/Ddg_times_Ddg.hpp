/* Fully contract a Ddg with a Ddg. */

#pragma once

namespace FTensor
{
  /* A(i,j,k,l)*B(i,k,j,l) */

  template <class A, class B, class T, class U, int Dim, char i, char j,
            char k, char l, int Current_Dim0, int Current_Dim1,
            int Current_Dim2, int Current_Dim3>
  typename promote<T, U>::V
  T4ddg_times_T4ddg_0213(const Ddg_Expr<A, T, Dim, Dim, i, j, k, l> &a,
                         const Ddg_Expr<B, U, Dim, Dim, i, k, j, l> &b,
                         const Number<Current_Dim0> &,
                         const Number<Current_Dim1> &,
                         const Number<Current_Dim2> &,
                         const Number<Current_Dim3> &)
  {
    return a(Current_Dim0 - 1, Current_Dim1 - 1, Current_Dim2 - 1,
             Current_Dim3 - 1)
             * b(Current_Dim0 - 1, Current_Dim2 - 1, Current_Dim1 - 1,
                 Current_Dim3 - 1)
           + T4ddg_times_T4ddg_0213(
               a, b, Number<Current_Dim0 - 1>(), Number<Current_Dim1>(),
               Number<Current_Dim2>(), Number<Current_Dim3>());
  }

  template <class A, class B, class T, class U, int Dim, char i, char j, char k,
            char l, int Current_Dim1, int Current_Dim2, int Current_Dim3>
  typename promote<T, U>::V
  T4ddg_times_T4ddg_0213(const Ddg_Expr<A, T, Dim, Dim, i, j, k, l> &a,
                         const Ddg_Expr<B, U, Dim, Dim, i, k, j, l> &b,
                         const Number<1> &, const Number<Current_Dim1> &,
                         const Number<Current_Dim2> &,
                         const Number<Current_Dim3> &)
  {
    return a(0, Current_Dim1 - 1, Current_Dim2 - 1, Current_Dim3 - 1)
             * b(0, Current_Dim2 - 1, Current_Dim1 - 1, Current_Dim3 - 1)
           + T4ddg_times_T4ddg_0213(
               a, b, Number<Dim>(), Number<Current_Dim1 - 1>(),
               Number<Current_Dim2>(), Number<Current_Dim3>());
  }

  template <class A, class B, class T, class U, int Dim, char i, char j,
            char k, char l, int Current_Dim2, int Current_Dim3>
  typename promote<T, U>::V
  T4ddg_times_T4ddg_0213(const Ddg_Expr<A, T, Dim, Dim, i, j, k, l> &a,
                         const Ddg_Expr<B, U, Dim, Dim, i, k, j, l> &b,
                         const Number<1> &, const Number<1> &,
                         const Number<Current_Dim2> &,
                         const Number<Current_Dim3> &)
  {
    return a(0, 0, Current_Dim2 - 1, Current_Dim3 - 1)
             * b(0, Current_Dim2 - 1, 0, Current_Dim3 - 1)
           + T4ddg_times_T4ddg_0213(a, b, Number<Dim>(), Number<Dim>(),
                                    Number<Current_Dim2 - 1>(),
                                    Number<Current_Dim3>());
  }

  template <class A, class B, class T, class U, int Dim, char i, char j,
            char k, char l, int Current_Dim3>
  typename promote<T, U>::V
  T4ddg_times_T4ddg_0213(const Ddg_Expr<A, T, Dim, Dim, i, j, k, l> &a,
                         const Ddg_Expr<B, U, Dim, Dim, i, k, j, l> &b,
                         const Number<1> &, const Number<1> &,
                         const Number<1> &, const Number<Current_Dim3> &)
  {
    return a(0, 0, 0, Current_Dim3 - 1) * b(0, 0, 0, Current_Dim3 - 1)
           + T4ddg_times_T4ddg_0213(a, b, Number<Dim>(), Number<Dim>(),
                                    Number<Dim>(), Number<Current_Dim3 - 1>());
  }

  template <class A, class B, class T, class U, int Dim, char i, char j,
            char k, char l>
  typename promote<T, U>::V
  T4ddg_times_T4ddg_0213(const Ddg_Expr<A, T, Dim, Dim, i, j, k, l> &a,
                         const Ddg_Expr<B, U, Dim, Dim, i, k, j, l> &b,
                         const Number<1> &, const Number<1> &,
                         const Number<1> &, const Number<1> &)
  {
    return a(0, 0, 0, 0) * b(0, 0, 0, 0);
  }

  template <class A, class B, class T, class U, int Dim, char i, char j,
            char k, char l>
  typename promote<T, U>::V
  operator*(const Ddg_Expr<A, T, Dim, Dim, i, j, k, l> &a,
            const Ddg_Expr<B, U, Dim, Dim, i, k, j, l> &b)
  {
    return T4ddg_times_T4ddg_0213(a, b, Number<Dim>(), Number<Dim>(),
                                  Number<Dim>(), Number<Dim>());
  }

  /* A(m,m,i,j)*B(m,nk,l) -> Ddg */

  template <class A, class B, class T, class U, int Dim01, int Dim23, int Dim45,
            char i, char j, char k, char l, char m, char n>
  class Ddg_times_Ddg_0101 {

    using IterA = Ddg_Expr<A, T, Dim45, Dim01, m, n, i, j>;
    using IterB = Ddg_Expr<B, U, Dim45, Dim23, m, n, k, l>;

    IterA iterA;
    IterB iterB;

  public:
    Ddg_times_Ddg_0101(const IterA &a, const IterB &b) : iterA(a), iterB(b) {}

    typename promote<T, U>::V operator()(const int N0, const int N1,
                                         const int N2, const int N3) const {
      typename promote<T, U>::V ret_val = 0;
      auto index_sequence = std::make_index_sequence<Dim45>();

      auto outer = [&](auto N4) {
        auto inner = [&](auto N5) {
          ret_val += iterA(N4, N5, N0, N1) * iterA(N4, N5, N2, N3);
        };
        boost::hana::for_each(index_sequence, inner);
      };
      boost::hana::for_each(index_sequence, outer);

      return ret_val;
    }
  };

  template <class A, class B, class T, class U, int Dim01, int Dim23, int Dim45,
            char i, char j, char k, char l, char m, char n>
  Ddg_Expr<
      Ddg_times_Ddg_0101<A, B, T, U, Dim01, Dim23, Dim45, i, j, k, l, m, n>,
      typename promote<T, U>::V, Dim01, Dim23, i, j, k, l>
  operator*(const Ddg_Expr<A, T, Dim45, Dim01, m, n, i, j> &a,
            const Ddg_Expr<B, U, Dim45, Dim23, m, n, k, l> &b) {
    using TensorExpr =
        Ddg_times_Ddg_0101<A, B, T, U, Dim01, Dim23, Dim45, i, j, k, l, m, n>;
    return Ddg_Expr<TensorExpr, typename promote<T, U>::V, Dim01, Dim23, i, j,
                    k, l>(TensorExpr(a, b));
  }
}
