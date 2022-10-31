/* Fully contract a Ddg with a Ddg. */

#pragma once

namespace FTensor {
/* A(i,j,k,l)*B(i,k,j,l) */

template <class A, class B, class T, class U, int Dim, char i, char j, char k,
          char l, int Current_Dim0, int Current_Dim1, int Current_Dim2,
          int Current_Dim3>
typename promote<T, U>::V T4ddg_times_T4ddg_0213(
    const Ddg_Expr<A, T, Dim, Dim, i, j, k, l> &a,
    const Ddg_Expr<B, U, Dim, Dim, i, k, j, l> &b, const Number<Current_Dim0> &,
    const Number<Current_Dim1> &, const Number<Current_Dim2> &,
    const Number<Current_Dim3> &) {
  return a(Current_Dim0 - 1, Current_Dim1 - 1, Current_Dim2 - 1,
           Current_Dim3 - 1) *
             b(Current_Dim0 - 1, Current_Dim2 - 1, Current_Dim1 - 1,
               Current_Dim3 - 1) +
         T4ddg_times_T4ddg_0213(a, b, Number<Current_Dim0 - 1>(),
                                Number<Current_Dim1>(), Number<Current_Dim2>(),
                                Number<Current_Dim3>());
}

template <class A, class B, class T, class U, int Dim, char i, char j, char k,
          char l, int Current_Dim1, int Current_Dim2, int Current_Dim3>
typename promote<T, U>::V
T4ddg_times_T4ddg_0213(const Ddg_Expr<A, T, Dim, Dim, i, j, k, l> &a,
                       const Ddg_Expr<B, U, Dim, Dim, i, k, j, l> &b,
                       const Number<1> &, const Number<Current_Dim1> &,
                       const Number<Current_Dim2> &,
                       const Number<Current_Dim3> &) {
  return a(0, Current_Dim1 - 1, Current_Dim2 - 1, Current_Dim3 - 1) *
             b(0, Current_Dim2 - 1, Current_Dim1 - 1, Current_Dim3 - 1) +
         T4ddg_times_T4ddg_0213(a, b, Number<Dim>(), Number<Current_Dim1 - 1>(),
                                Number<Current_Dim2>(), Number<Current_Dim3>());
}

template <class A, class B, class T, class U, int Dim, char i, char j, char k,
          char l, int Current_Dim2, int Current_Dim3>
typename promote<T, U>::V
T4ddg_times_T4ddg_0213(const Ddg_Expr<A, T, Dim, Dim, i, j, k, l> &a,
                       const Ddg_Expr<B, U, Dim, Dim, i, k, j, l> &b,
                       const Number<1> &, const Number<1> &,
                       const Number<Current_Dim2> &,
                       const Number<Current_Dim3> &) {
  return a(0, 0, Current_Dim2 - 1, Current_Dim3 - 1) *
             b(0, Current_Dim2 - 1, 0, Current_Dim3 - 1) +
         T4ddg_times_T4ddg_0213(a, b, Number<Dim>(), Number<Dim>(),
                                Number<Current_Dim2 - 1>(),
                                Number<Current_Dim3>());
}

template <class A, class B, class T, class U, int Dim, char i, char j, char k,
          char l, int Current_Dim3>
typename promote<T, U>::V
T4ddg_times_T4ddg_0213(const Ddg_Expr<A, T, Dim, Dim, i, j, k, l> &a,
                       const Ddg_Expr<B, U, Dim, Dim, i, k, j, l> &b,
                       const Number<1> &, const Number<1> &, const Number<1> &,
                       const Number<Current_Dim3> &) {
  return a(0, 0, 0, Current_Dim3 - 1) * b(0, 0, 0, Current_Dim3 - 1) +
         T4ddg_times_T4ddg_0213(a, b, Number<Dim>(), Number<Dim>(),
                                Number<Dim>(), Number<Current_Dim3 - 1>());
}

template <class A, class B, class T, class U, int Dim, char i, char j, char k,
          char l>
typename promote<T, U>::V
T4ddg_times_T4ddg_0213(const Ddg_Expr<A, T, Dim, Dim, i, j, k, l> &a,
                       const Ddg_Expr<B, U, Dim, Dim, i, k, j, l> &b,
                       const Number<1> &, const Number<1> &, const Number<1> &,
                       const Number<1> &) {
  return a(0, 0, 0, 0) * b(0, 0, 0, 0);
}

template <class A, class B, class T, class U, int Dim, char i, char j, char k,
          char l>
typename promote<T, U>::V
operator*(const Ddg_Expr<A, T, Dim, Dim, i, j, k, l> &a,
          const Ddg_Expr<B, U, Dim, Dim, i, k, j, l> &b) {
  return T4ddg_times_T4ddg_0213(a, b, Number<Dim>(), Number<Dim>(),
                                Number<Dim>(), Number<Dim>());
}

/* A(m,n,i,j)*B(m,n,k,l) -> Ddg */

template <class A, class B, class T, class U, int Dim01, int Dim23, int Dim45,
          char i, char j, char k, char l, char m, char n>
class Ddg_times_Ddg_0101 {

  using IterA = Ddg_Expr<A, T, Dim45, Dim01, m, n, i, j>;
  using IterB = Ddg_Expr<B, U, Dim45, Dim23, m, n, k, l>;

  IterA iterA;
  IterB iterB;

public:
  Ddg_times_Ddg_0101(const IterA &a, const IterB &b) : iterA(a), iterB(b) {}

  template <int Current_Dim0, int Current_Dim1>
  inline typename promote<T, U>::V eval(const int N1, const int N2, const int N3,
                                 const int N4, const Number<Current_Dim0> &,
                                 const Number<Current_Dim1> &) const {
    return iterA(Current_Dim0 - 1, Current_Dim1 - 1, N1, N2) *
               iterB(Current_Dim0 - 1, Current_Dim1 - 1, N3, N4) +
           eval(N1, N2, N3, N4, Number<Current_Dim0 - 1>(),
                Number<Current_Dim1>());
  }
  template <int Current_Dim1>
  inline typename promote<T, U>::V eval(const int N1, const int N2, const int N3,
                                 const int N4, const Number<1> &,
                                 const Number<Current_Dim1> &) const {
    return iterA(0, Current_Dim1 - 1, N1, N2) *
               iterB(0, Current_Dim1 - 1, N3, N4) +
           eval(N1, N2, N3, N4, Number<Dim45>(),
                Number<Current_Dim1 - 1>());
  }
  inline typename promote<T, U>::V eval(const int N1, const int N2, const int N3,
                                 const int N4, const Number<1> &,
                                 const Number<1> &) const {
    return iterA(0, 0, N1, N2) * iterB(0, 0, N3, N4);
  }

  typename promote<T, U>::V operator()(const int N1, const int N2, const int N3,
                                       const int N4) const {
    return eval(N1, N2, N3, N4, Number<Dim45>(), Number<Dim45>());
  }
};

template <class A, class B, class T, class U, int Dim01, int Dim23, int Dim45,
          char i, char j, char k, char l, char m, char n>
inline auto operator*(const Ddg_Expr<A, T, Dim45, Dim01, m, n, i, j> &a,
                      const Ddg_Expr<B, U, Dim45, Dim23, m, n, k, l> &b) {
  using TensorExpr =
      Ddg_times_Ddg_0101<A, B, T, U, Dim01, Dim23, Dim45, i, j, k, l, m, n>;
  return Ddg_Expr<TensorExpr, typename promote<T, U>::V, Dim01, Dim23, i, j, k,
                  l>(TensorExpr(a, b));
}

/* B(m,n,k,l)* A(m,n,i,j) -> Ddg */

template <class A, class B, class T, class U, int Dim01, int Dim23, int Dim45,
          char i, char j, char k, char l, char m, char n>
inline auto operator*(const Ddg_Expr<B, U, Dim45, Dim23, m, n, k, l> &b,
                      const Ddg_Expr<A, T, Dim45, Dim01, m, n, i, j> &a) {
  using TensorExpr =
      Ddg_times_Ddg_0101<A, B, T, U, Dim01, Dim23, Dim45, i, j, k, l, m, n>;
  return Ddg_Expr<TensorExpr, typename promote<T, U>::V, Dim01, Dim23, i, j, k,
                  l>(TensorExpr(a, b));
}

/* A(i,j, m, n)*B(m,n,k,l) -> Ddg */

template <class A, class B, class T, class U, int Dim01, int Dim23, int Dim45,
          char i, char j, char k, char l, char m, char n>
class Ddg_times_Ddg_2301 {

  using IterA = Ddg_Expr<A, T, Dim45, Dim01, i, j, m, n>;
  using IterB = Ddg_Expr<B, U, Dim45, Dim23, m, n, k, l>;

  IterA iterA;
  IterB iterB;

public:
  Ddg_times_Ddg_2301(const IterA &a, const IterB &b) : iterA(a), iterB(b) {}

  template <int Current_Dim0, int Current_Dim1>
  inline typename promote<T, U>::V eval(const int N1, const int N2, const int N3,
                                 const int N4, const Number<Current_Dim0> &,
                                 const Number<Current_Dim1> &) const {
    return iterA(N1, N2, Current_Dim0 - 1, Current_Dim1 - 1) *
               iterB(Current_Dim0 - 1, Current_Dim1 - 1, N3, N4) +
           eval(N1, N2, N3, N4, Number<Current_Dim0 - 1>(),
                Number<Current_Dim1>());
  }
  template <int Current_Dim1>
  inline typename promote<T, U>::V eval(const int N1, const int N2, const int N3,
                                 const int N4, const Number<1> &,
                                 const Number<Current_Dim1> &) const {
    return iterA(N1, N2, 0, Current_Dim1 - 1) *
               iterB(0, Current_Dim1 - 1, N3, N4) +
           eval(N1, N2, N3, N4, Number<Dim45>(), Number<Current_Dim1 - 1>());
  }
  inline typename promote<T, U>::V eval(const int N1, const int N2, const int N3,
                                 const int N4, const Number<1> &,
                                 const Number<1> &) const {
    return iterA(N1, N2, 0, 0) * iterB(0, 0, N3, N4);
  }

  typename promote<T, U>::V operator()(const int N1, const int N2, const int N3,
                                       const int N4) const {
    return eval(N1, N2, N3, N4, Number<Dim45>(), Number<Dim45>());
  }

};

template <class A, class B, class T, class U, int Dim01, int Dim23, int Dim45,
          char i, char j, char k, char l, char m, char n>
Ddg_Expr<Ddg_times_Ddg_2301<A, B, T, U, Dim01, Dim23, Dim45, i, j, k, l, m, n>,
         typename promote<T, U>::V, Dim01, Dim23, i, j, k, l>
operator*(const Ddg_Expr<A, T, Dim01, Dim45, i, j, m, n> &a,
          const Ddg_Expr<B, U, Dim45, Dim23, m, n, k, l> &b) {
  using TensorExpr =
      Ddg_times_Ddg_2301<A, B, T, U, Dim01, Dim23, Dim45, i, j, k, l, m, n>;
  return Ddg_Expr<TensorExpr, typename promote<T, U>::V, Dim01, Dim23, i, j, k,
                  l>(TensorExpr(a, b));
}

/* A(m,n,i,j)*B(k,l,m,n) -> Ddg */

template <class A, class B, class T, class U, int Dim01, int Dim23, int Dim45,
          char i, char j, char k, char l, char m, char n>
class Ddg_times_Ddg_0123 {

  using IterA = Ddg_Expr<A, T, Dim45, Dim01, m, n, i, j>;
  using IterB = Ddg_Expr<B, U, Dim23, Dim45, k, l, m, n>;

  IterA iterA;
  IterB iterB;

public:
  Ddg_times_Ddg_0123(const IterA &a, const IterB &b) : iterA(a), iterB(b) {}

  template <int Current_Dim0, int Current_Dim1>
  inline typename promote<T, U>::V eval(const int N1, const int N2, const int N3,
                                 const int N4, const Number<Current_Dim0> &,
                                 const Number<Current_Dim1> &) const {
    return iterA(Current_Dim0 - 1, Current_Dim1 - 1, N1, N2) *
               iterB(N3, N4, Current_Dim0 - 1, Current_Dim1 - 1) +
           eval(N1, N2, N3, N4, Number<Current_Dim0 - 1>(),
                Number<Current_Dim1>());
  }
  template <int Current_Dim1>
  inline typename promote<T, U>::V eval(const int N1, const int N2, const int N3,
                                 const int N4, const Number<1> &,
                                 const Number<Current_Dim1> &) const {
    return iterA(0, Current_Dim1 - 1, N1, N2) *
               iterB(N3, N4, 0, Current_Dim1 - 1) +
           eval(N1, N2, N3, N4, Number<Dim45>(), Number<Current_Dim1 - 1>());
  }
  inline typename promote<T, U>::V eval(const int N1, const int N2, const int N3,
                                 const int N4, const Number<1> &,
                                 const Number<1> &) const {
    return iterA(0, 0, N1, N2) * iterB(N3, N4, 0, 0);
  }

  typename promote<T, U>::V operator()(const int N1, const int N2, const int N3,
                                       const int N4) const {
    return eval(N1, N2, N3, N4, Number<Dim45>(), Number<Dim45>());
  }


};

template <class A, class B, class T, class U, int Dim01, int Dim23, int Dim45,
          char i, char j, char k, char l, char m, char n>
Ddg_Expr<Ddg_times_Ddg_0123<A, B, T, U, Dim01, Dim23, Dim45, i, j, k, l, m, n>,
         typename promote<T, U>::V, Dim01, Dim23, i, j, k, l>
operator*(const Ddg_Expr<A, T, Dim45, Dim01, m, n, i, j> &a,
          const Ddg_Expr<B, U, Dim23, Dim45, k, l, m, n> &b) {
  using TensorExpr =
      Ddg_times_Ddg_0123<A, B, T, U, Dim01, Dim23, Dim45, i, j, k, l, m, n>;
  return Ddg_Expr<TensorExpr, typename promote<T, U>::V, Dim01, Dim23, i, j, k,
                  l>(TensorExpr(a, b));
}

/* A(i,j,m,n)*B(k,l,m,n) -> Ddg */

template <class A, class B, class T, class U, int Dim01, int Dim23, int Dim45,
          char i, char j, char k, char l, char m, char n>
class Ddg_times_Ddg_2323 {

  using IterA = Ddg_Expr<A, T, Dim45, Dim01, i, j, m, n>;
  using IterB = Ddg_Expr<B, U, Dim23, Dim45, k, l, m, n>;

  IterA iterA;
  IterB iterB;

public:
  Ddg_times_Ddg_2323(const IterA &a, const IterB &b) : iterA(a), iterB(b) {}

  template <int Current_Dim0, int Current_Dim1>
  inline typename promote<T, U>::V
  eval(const int N1, const int N2, const int N3, const int N4,
       const Number<Current_Dim0> &, const Number<Current_Dim1> &) const {
    return iterA(N1, N2, Current_Dim0 - 1, Current_Dim1 - 1) *
               iterB(N3, N4, Current_Dim0 - 1, Current_Dim1 - 1) +
           eval(N1, N2, N3, N4, Number<Current_Dim0 - 1>(),
                Number<Current_Dim1>());
  }
  template <int Current_Dim1>
  inline typename promote<T, U>::V
  eval(const int N1, const int N2, const int N3, const int N4,
       const Number<1> &, const Number<Current_Dim1> &) const {
    return iterA(N1, N2, 0, Current_Dim1 - 1) *
               iterB(N3, N4, 0, Current_Dim1 - 1) +
           eval(N1, N2, N3, N4, Number<Dim45>(), Number<Current_Dim1 - 1>());
  }
  inline typename promote<T, U>::V eval(const int N1, const int N2,
                                        const int N3, const int N4,
                                        const Number<1> &,
                                        const Number<1> &) const {
    return iterA(N1, N2, 0, 0) * iterB(N3, N4, 0, 0);
  }

  typename promote<T, U>::V operator()(const int N1, const int N2, const int N3,
                                       const int N4) const {
    return eval(N1, N2, N3, N4, Number<Dim45>(), Number<Dim45>());
  }

};

template <class A, class B, class T, class U, int Dim01, int Dim23, int Dim45,
          char i, char j, char k, char l, char m, char n>
Ddg_Expr<Ddg_times_Ddg_2323<A, B, T, U, Dim01, Dim23, Dim45, i, j, k, l, m, n>,
         typename promote<T, U>::V, Dim01, Dim23, i, j, k, l>
operator*(const Ddg_Expr<A, T, Dim01, Dim45, i, j, m, n> &a,
          const Ddg_Expr<B, U, Dim23, Dim45, k, l, m, n> &b) {
  using TensorExpr =
      Ddg_times_Ddg_2323<A, B, T, U, Dim01, Dim23, Dim45, i, j, k, l, m, n>;
  return Ddg_Expr<TensorExpr, typename promote<T, U>::V, Dim01, Dim23, i, j, k,
                  l>(TensorExpr(a, b));
}

} // namespace FTensor
