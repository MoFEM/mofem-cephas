/* Various assignments operators */

/* T4=T4 */

#pragma once

namespace FTensor
{
  /* T4 = T4 */

  template<class A,class B,class U,int Dim0,int Dim1,int Dim2,int Dim3,
          char i,char j,char k,char l,
          int Current_Dim0,int Current_Dim1,int Current_Dim2,int Current_Dim3>
  inline void T4_equals_T4(A &iter,const Tensor4_Expr<B,U,Dim0,Dim1,Dim2,Dim3,i,j,k,l> &result,
                           const Number<Current_Dim0> &,
                           const Number<Current_Dim1> &,
                           const Number<Current_Dim2> &,
                           const Number<Current_Dim3> &)
  {
    iter(Current_Dim0-1,Current_Dim1-1,Current_Dim2-1,Current_Dim3-1) =
            result(Current_Dim0-1,Current_Dim1-1,Current_Dim2-1,Current_Dim3-1);
    T4_equals_T4(iter,result,
                 Number<Current_Dim0-1>(),
                 Number<Current_Dim1>(),
                 Number<Current_Dim2>(),
                 Number<Current_Dim3>());
  }

  template<class A,class B,class U,int Dim0,int Dim1,int Dim2,int Dim3,
          char i,char j,char k,char l,
          int Current_Dim1,int Current_Dim2,int Current_Dim3>
  inline void T4_equals_T4(A &iter,const Tensor4_Expr<B,U,Dim0,Dim1,Dim2,Dim3,i,j,k,l> &result,
                           const Number<1> &,
                           const Number<Current_Dim1> &,
                           const Number<Current_Dim2> &,
                           const Number<Current_Dim3> &)
  {
    iter(0,Current_Dim1-1,Current_Dim2-1,Current_Dim3-1) =
            result(0,Current_Dim1-1,Current_Dim2-1,Current_Dim3-1);
    T4_equals_T4(iter,result,
                 Number<Dim0>(),
                 Number<Current_Dim1-1>(),
                 Number<Current_Dim2>(),
                 Number<Current_Dim3>());
  }

  template<class A,class B,class U,int Dim0,int Dim1,int Dim2,int Dim3,
          char i,char j,char k,char l,
          int Current_Dim2,int Current_Dim3>
  inline void T4_equals_T4(A &iter,const Tensor4_Expr<B,U,Dim0,Dim1,Dim2,Dim3,i,j,k,l> &result,
                           const Number<1> &,
                           const Number<1> &,
                           const Number<Current_Dim2> &,
                           const Number<Current_Dim3> &)
  {
    iter(0,0,Current_Dim2-1,Current_Dim3-1) =
            result(0,0,Current_Dim2-1,Current_Dim3-1);
    T4_equals_T4(iter,result,
                 Number<Dim0>(),
                 Number<Dim1>(),
                 Number<Current_Dim2-1>(),
                 Number<Current_Dim3>());
  }

  template<class A,class B,class U,int Dim0,int Dim1,int Dim2,int Dim3,
          char i,char j,char k,char l,
          int Current_Dim3>
  inline void T4_equals_T4(A &iter,const Tensor4_Expr<B,U,Dim0,Dim1,Dim2,Dim3,i,j,k,l> &result,
                           const Number<1> &,
                           const Number<1> &,
                           const Number<1> &,
                           const Number<Current_Dim3> &)
  {
    iter(0,0,0,Current_Dim3-1) =
            result(0,0,0,Current_Dim3-1);
    T4_equals_T4(iter,result,
                 Number<Dim0>(),
                 Number<Dim1>(),
                 Number<Dim2>(),
                 Number<Current_Dim3-1>());
  }

  template<class A,class B,class U,int Dim0,int Dim1,int Dim2,int Dim3,
          char i,char j,char k,char l>
  inline void T4_equals_T4(A &iter,const Tensor4_Expr<B,U,Dim0,Dim1,Dim2,Dim3,i,j,k,l> &result,
                           const Number<1> &,
                           const Number<1> &,
                           const Number<1> &,
                           const Number<1> &)
  {
    iter(0,0,0,0) = result(0,0,0,0);
  }

  template<class A, class T, int Dim0, int Dim1, int Dim2, int Dim3,
          char i, char j, char k, char l>
  template<class B, class U> inline
  const Tensor4_Expr<Tensor4<A,Dim0,Dim1,Dim2,Dim3>,T,Dim0,Dim1,Dim2,Dim3,i,j,k,l> &
  Tensor4_Expr<Tensor4<A,Dim0,Dim1,Dim2,Dim3>,T,Dim0,Dim1,Dim2,Dim3,i,j,k,l>::
  operator=(const Tensor4_Expr<B,U,Dim0,Dim1,Dim2,Dim3,i,j,k,l> &result)
  {
    T4_equals_T4(iter,result,Number<Dim0>(),Number<Dim1>(),
                           Number<Dim2>(),Number<Dim3>());
    return *this;
  }

  /* T2=T2_Expr(T2) */

  template<class A, class T, int Dim0, int Dim1, int Dim2, int Dim3,
          char i, char j, char k, char l>
  inline const Tensor4_Expr<Tensor4<A,Dim0,Dim1,Dim2,Dim3>,T,Dim0,Dim1,Dim2,Dim3,i,j,k,l> &
          Tensor4_Expr<Tensor4<A,Dim0,Dim1,Dim2,Dim3>,T,Dim0,Dim1,Dim2,Dim3,i,j,k,l>::
  operator =(const Tensor4_Expr<Tensor4<A,Dim0,Dim1,Dim2,Dim3>,T,Dim0,Dim1,Dim2,Dim3,i,j,k,l> &result)
  {
    return operator=<Tensor4<A,Dim0,Dim1,Dim2,Dim3>,T>(result);
  }

  /* T4(i,j,k,l) = T4(i,l,j,k) */

  template<class A,class B,class U,int Dim0,int Dim1,int Dim2,int Dim3,
           char i,char j,char k,char l,
           int Current_Dim0,int Current_Dim1,int Current_Dim2,int Current_Dim3>
  inline void T4_equals_T4_0312(A &iter,const Tensor4_Expr<B,U,Dim0,Dim1,Dim2,Dim3,i,j,k,l> &result,
                                const Number<Current_Dim0> &,
                                const Number<Current_Dim1> &,
                                const Number<Current_Dim2> &,
                                const Number<Current_Dim3> &)
  {
    iter(Current_Dim0-1,Current_Dim1-1,Current_Dim2-1,Current_Dim3-1) =
            result(Current_Dim0-1,Current_Dim3-1,Current_Dim1-1,Current_Dim2-1);
    T4_equals_T4_0312(iter,result,
                      Number<Current_Dim0-1>(),
                      Number<Current_Dim1>(),
                      Number<Current_Dim2>(),
                      Number<Current_Dim3>());
  }

  template<class A,class B,class U,int Dim0,int Dim1,int Dim2,int Dim3,
           char i,char j,char k,char l,
           int Current_Dim1,int Current_Dim2,int Current_Dim3>
  inline void T4_equals_T4_0312(A &iter,const Tensor4_Expr<B,U,Dim0,Dim1,Dim2,Dim3,i,j,k,l> &result,
                                const Number<1> &,
                                const Number<Current_Dim1> &,
                                const Number<Current_Dim2> &,
                                const Number<Current_Dim3> &)
  {
    iter(0,Current_Dim1-1,Current_Dim2-1,Current_Dim3-1) =
            result(0,Current_Dim3-1,Current_Dim1-1,Current_Dim2-1);
    T4_equals_T4_0312(iter,result,
                      Number<Dim0>(),
                      Number<Current_Dim1-1>(),
                      Number<Current_Dim2>(),
                      Number<Current_Dim3>());
  }

  template<class A,class B,class U,int Dim0,int Dim1,int Dim2,int Dim3,
           char i,char j,char k,char l,
           int Current_Dim2,int Current_Dim3>
  inline void T4_equals_T4_0312(A &iter,const Tensor4_Expr<B,U,Dim0,Dim1,Dim2,Dim3,i,j,k,l> &result,
                                const Number<1> &,
                                const Number<1> &,
                                const Number<Current_Dim2> &,
                                const Number<Current_Dim3> &)
  {
    iter(0,0,Current_Dim2-1,Current_Dim3-1) =
            result(0,Current_Dim3-1,0,Current_Dim2-1);
    T4_equals_T4_0312(iter,result,
                      Number<Dim0>(),
                      Number<Dim1>(),
                      Number<Current_Dim2-1>(),
                      Number<Current_Dim3>());
  }

  template<class A,class B,class U,int Dim0,int Dim1,int Dim2,int Dim3,
           char i,char j,char k,char l,
           int Current_Dim3>
  inline void T4_equals_T4_0312(A &iter,const Tensor4_Expr<B,U,Dim0,Dim1,Dim2,Dim3,i,j,k,l> &result,
                                const Number<1> &,
                                const Number<1> &,
                                const Number<1> &,
                                const Number<Current_Dim3> &)
  {
    iter(0,0,0,Current_Dim3-1) =
            result(0,Current_Dim3-1,0,0);
    T4_equals_T4_0312(iter,result,
                      Number<Dim0>(),
                      Number<Dim1>(),
                      Number<Dim2>(),
                      Number<Current_Dim3-1>());
  }

  template<class A,class B,class U,int Dim0,int Dim1,int Dim2,int Dim3,
           char i,char j,char k,char l>
  inline void T4_equals_T4_0312(A &iter,const Tensor4_Expr<B,U,Dim0,Dim1,Dim2,Dim3,i,j,k,l> &result,
                                const Number<1> &,
                                const Number<1> &,
                                const Number<1> &,
                                const Number<1> &)
  {
    iter(0,0,0,0) = result(0,0,0,0);
  }

  template<class A, class T, int Dim0, int Dim1, int Dim2, int Dim3,
           char i, char j, char k, char l>
  template<class B, class U> inline
  const Tensor4_Expr<Tensor4<A,Dim0,Dim1,Dim2,Dim3>,T,Dim0,Dim1,Dim2,Dim3,i,j,k,l> &
  Tensor4_Expr<Tensor4<A,Dim0,Dim1,Dim2,Dim3>,T,Dim0,Dim1,Dim2,Dim3,i,j,k,l>::
  operator=(const Tensor4_Expr<B,U,Dim0,Dim3,Dim1,Dim2,i,l,j,k> &result)
  {
    T4_equals_T4_0312(iter,result,Number<Dim0>(),Number<Dim1>(),
                             Number<Dim2>(),Number<Dim3>());
    return *this;
  }

    // Assignment variants

#define EQUAL(T4_equals_T4_IJKL, D0, D1, D2, D3, I, J, K, L, S10, S11, S12,    \
              S13, S20, S21, S22, S23, S30, S31, S32, S33, S40, S41, S42, S43) \
  template <class A, class B, class U, int Dim0, int Dim1, int Dim2, int Dim3, \
            char i, char j, char k, char l, int Current_Dim0,                  \
            int Current_Dim1, int Current_Dim2, int Current_Dim3>              \
  inline void T4_equals_T4_IJKL(                                               \
      A &iter,                                                                 \
      const Tensor4_Expr<B, U, Dim0, Dim1, Dim2, Dim3, I, J, K, L> &result,    \
      const Number<Current_Dim0> &, const Number<Current_Dim1> &,              \
      const Number<Current_Dim2> &, const Number<Current_Dim3> &) {            \
    iter(Current_Dim0 - 1, Current_Dim1 - 1, Current_Dim2 - 1,                 \
         Current_Dim3 - 1) = result(S10, S11, S12, S13);                       \
    T4_equals_T4_IJKL(iter, result, Number<Current_Dim0 - 1>(),                \
                      Number<Current_Dim1>(), Number<Current_Dim2>(),          \
                      Number<Current_Dim3>());                                 \
  }                                                                            \
  template <class A, class B, class U, int Dim0, int Dim1, int Dim2, int Dim3, \
            char i, char j, char k, char l, int Current_Dim1,                  \
            int Current_Dim2, int Current_Dim3>                                \
  inline void T4_equals_T4_IJKL(                                               \
      A &iter,                                                                 \
      const Tensor4_Expr<B, U, Dim0, Dim1, Dim2, Dim3, I, J, K, L> &result,    \
      const Number<1> &, const Number<Current_Dim1> &,                         \
      const Number<Current_Dim2> &, const Number<Current_Dim3> &) {            \
    iter(0, Current_Dim1 - 1, Current_Dim2 - 1, Current_Dim3 - 1) =            \
        result(S20, S21, S22, S23);                                            \
    T4_equals_T4_IJKL(iter, result, Number<Dim0>(),                            \
                      Number<Current_Dim1 - 1>(), Number<Current_Dim2>(),      \
                      Number<Current_Dim3>());                                 \
  }                                                                            \
  template <class A, class B, class U, int Dim0, int Dim1, int Dim2, int Dim3, \
            char i, char j, char k, char l, int Current_Dim2,                  \
            int Current_Dim3>                                                  \
  inline void T4_equals_T4_IJKL(                                               \
      A &iter,                                                                 \
      const Tensor4_Expr<B, U, Dim0, Dim1, Dim2, Dim3, I, J, K, L> &result,    \
      const Number<1> &, const Number<1> &, const Number<Current_Dim2> &,      \
      const Number<Current_Dim3> &) {                                          \
    iter(0, 0, Current_Dim2 - 1, Current_Dim3 - 1) =                           \
        result(S30, S31, S32, S33);                                            \
    T4_equals_T4_IJKL(iter, result, Number<Dim0>(), Number<Dim1>(),            \
                      Number<Current_Dim2 - 1>(), Number<Current_Dim3>());     \
  }                                                                            \
  template <class A, class B, class U, int Dim0, int Dim1, int Dim2, int Dim3, \
            char i, char j, char k, char l, int Current_Dim3>                  \
  inline void T4_equals_T4_IJKL(                                               \
      A &iter,                                                                 \
      const Tensor4_Expr<B, U, Dim0, Dim1, Dim2, Dim3, I, J, K, L> &result,    \
      const Number<1> &, const Number<1> &, const Number<1> &,                 \
      const Number<Current_Dim3> &) {                                          \
    iter(0, 0, 0, Current_Dim3 - 1) = result(S40, S41, S42, S43);              \
    T4_equals_T4_IJKL(iter, result, Number<Dim0>(), Number<Dim1>(),            \
                      Number<Dim2>(), Number<Current_Dim3 - 1>());             \
  }                                                                            \
  template <class A, class B, class U, int Dim0, int Dim1, int Dim2, int Dim3, \
            char i, char j, char k, char l>                                    \
  inline void T4_equals_T4_IJKL(                                               \
      A &iter,                                                                 \
      const Tensor4_Expr<B, U, Dim0, Dim1, Dim2, Dim3, I, J, K, L> &result,    \
      const Number<1> &, const Number<1> &, const Number<1> &,                 \
      const Number<1> &) {                                                     \
    iter(0, 0, 0, 0) = result(0, 0, 0, 0);                                     \
  }                                                                            \
  template <class A, class T, int Dim0, int Dim1, int Dim2, int Dim3, char i,  \
            char j, char k, char l>                                            \
  template <class B, class U>                                                  \
  inline const Tensor4_Expr<Tensor4<A, Dim0, Dim1, Dim2, Dim3>, T, Dim0, Dim1, \
                            Dim2, Dim3, i, j, k, l> &                          \
  Tensor4_Expr<Tensor4<A, Dim0, Dim1, Dim2, Dim3>, T, Dim0, Dim1, Dim2, Dim3,  \
               i, j, k, l>::                                                   \
  operator=(const Tensor4_Expr<B, U, D0, D1, D2, D3, I, J, K, L> &result) {    \
    T4_equals_T4_IJKL(iter, result, Number<Dim0>(), Number<Dim1>(),            \
                      Number<Dim2>(), Number<Dim3>());                         \
    return *this;                                                              \
  }

  // kijl
  EQUAL(T4_equals_T4_2013, Dim2, Dim0, Dim1, Dim3, k, i, j, l, Current_Dim2 - 1,
        Current_Dim0 - 1, Current_Dim1 - 1, Current_Dim3 - 1, Current_Dim2 - 1,
        0, Current_Dim1 - 1, Current_Dim3 - 1, Current_Dim2 - 1, 0, 0,
        Current_Dim3 - 1, 0, 0, 0, Current_Dim3 - 1);

#undef EQUAL
}
