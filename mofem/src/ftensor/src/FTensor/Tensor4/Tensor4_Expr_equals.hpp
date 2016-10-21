/* Various assignment operators. */

/* T4=T4 */

template<class A, class B, class U, int Dim0, int Dim1,int Dim2, int Dim3, char i, char j, char k, char l,
  int Current_Dim0, int Current_Dim1,int Current_Dim2, int Current_Dim3>
inline void T4_equals_T4(
  A &iter,
  const Tensor4_Expr<B,U,Dim0,Dim1,Dim2,Dim3,i,j,k,l> &result,
  const Number<Current_Dim0> &,
  const Number<Current_Dim1> &,
  const Number<Current_Dim2> &,
  const Number<Current_Dim3> &
) {
  iter(Current_Dim0-1,Current_Dim1-1,Current_Dim2-1,Current_Dim3-1)=result(Current_Dim0-1,Current_Dim1-1,Current_Dim2-1,Current_Dim3-1);
  T4_equals_T4(iter,result,Number<Current_Dim0-1>(),Number<Current_Dim1>(),Number<Current_Dim2>(),Number<Current_Dim3>());
}

template<class A, class B, class U, int Dim0, int Dim1,int Dim2, int Dim3, char i, char j, char k, char l,
  int Current_Dim1,int Current_Dim2,int Current_Dim3>
inline void T4_equals_T4(A &iter,
			 const Tensor4_Expr<B,U,Dim0,Dim1,Dim2,Dim3,i,j,k,l> &result,
			 const Number<1> &,
			 const Number<Current_Dim1> &,
       const Number<Current_Dim2> &,
       const Number<Current_Dim3> &
     )
{
  iter(0,Current_Dim1-1,Current_Dim2-1,Current_Dim3-1)=result(0,Current_Dim1-1,Current_Dim2-1,Current_Dim3-1);
  T4_equals_T4(iter,result,Number<Dim0>(),Number<Current_Dim1-1>(),Number<Current_Dim2>(),Number<Current_Dim3>());
}

template<class A, class B, class U, int Dim0, int Dim1,int Dim2, int Dim3, char i, char j, char k, char l,
  int Current_Dim2,int Current_Dim3>
inline void T4_equals_T4(A &iter,
			 const Tensor4_Expr<B,U,Dim0,Dim1,Dim2,Dim3,i,j,k,l> &result,
			 const Number<1> &,
			 const Number<1> &,
       const Number<Current_Dim2> &,
       const Number<Current_Dim3> &
     )
{
  iter(0,0,Current_Dim2-1,Current_Dim3-1)=result(0,0,Current_Dim2-1,Current_Dim3-1);
  T4_equals_T4(iter,result,Number<Dim0>(),Number<Dim1>(),Number<Current_Dim2-1>(),Number<Current_Dim3>());
}

template<class A, class B, class U, int Dim0, int Dim1,int Dim2, int Dim3, char i, char j, char k, char l,
  int Current_Dim3>
inline void T4_equals_T4(A &iter,
			 const Tensor4_Expr<B,U,Dim0,Dim1,Dim2,Dim3,i,j,k,l> &result,
			 const Number<1> &,
			 const Number<1> &,
       const Number<1> &,
       const Number<Current_Dim3> &
     )
{
  iter(0,0,0,Current_Dim3-1)=result(0,0,0,Current_Dim3-1);
  T4_equals_T4(iter,result,Number<Dim0>(),Number<Dim1>(),Number<Dim2>(),Number<Current_Dim3-1>());
}

template<class A, class B, class U, int Dim0, int Dim1,int Dim2,int Dim3, char i, char j,char k,char l>
inline void T4_equals_T4(A &iter,
  const Tensor4_Expr<B,U,Dim0,Dim1,Dim2,Dim3,i,j,k,l> &result,
  const Number<1> &,
  const Number<1> &,
  const Number<1> &,
  const Number<1> &
)
{
  iter(0,0,0,0)=result(0,0,0,0);
}

template<class A, class T, int Dim0, int Dim1,int Dim2,int Dim3, char i, char j, char k, char l>
template<class B, class U> inline
const Tensor4_Expr<Tensor4<A,Dim0,Dim1,Dim2,Dim3>,T,Dim0,Dim1,Dim2,Dim3,i,j,k,l> &
Tensor4_Expr<Tensor4<A,Dim0,Dim1,Dim2,Dim3>,T,Dim0,Dim1,Dim2,Dim3,i,j,k,l>::
operator=(const Tensor4_Expr<B,U,Dim0,Dim1,Dim2,Dim3,i,j,k,l> &result)
{
  T4_equals_T4(iter,result,Number<Dim0>(),Number<Dim1>(),Number<Dim2>(),Number<Dim3>());
  return *this;
}

/* T4=T4_Expr(T4): I have to explicitly declare this operator= because
   otherwise the compiler will generate its own and not use the
   template code. */

  template<class A, class T, int Dim0, int Dim1,int Dim2,int Dim3, char i, char j, char k, char l> inline
  const Tensor4_Expr<Tensor4<A,Dim0,Dim1,Dim2,Dim3>,T,Dim0,Dim1,Dim2,Dim3,i,j,k,l> &
  Tensor4_Expr<Tensor4<A,Dim0,Dim1,Dim2,Dim3>,T,Dim0,Dim1,Dim2,Dim3,i,j,k,l>::
  operator=(const Tensor4_Expr<Tensor4<A,Dim0,Dim1,Dim2,Dim3>,T,Dim0,Dim1,Dim2,Dim3,i,j,k,l> &result) {
    return operator=<Tensor4<A,Dim0,Dim1,Dim2,Dim3>,T>(result);
  }

  // Assignment variants

  #define EQUAL(T4_equals_T4_IJKL,D0,D1,D2,D3,I,J,K,L,S10,S11,S12,S13,S20,S21,S22,S23,S30,S31,S32,S33,S40,S41,S42,S43) \
  template<class A, class B, class U, int Dim0, int Dim1,int Dim2, int Dim3, char i, char j, char k, char l, \
    int Current_Dim0, int Current_Dim1,int Current_Dim2, int Current_Dim3> \
  inline void T4_equals_T4_IJKL( \
    A &iter, \
    const Tensor4_Expr<B,U,Dim0,Dim1,Dim2,Dim3,I,J,K,L> &result, \
    const Number<Current_Dim0> &, \
    const Number<Current_Dim1> &, \
    const Number<Current_Dim2> &, \
    const Number<Current_Dim3> & \
  ) { \
    iter(Current_Dim0-1,Current_Dim1-1,Current_Dim2-1,Current_Dim3-1)=result(S10,S11,S12,S13); \
    T4_equals_T4_IJKL(iter,result,Number<Current_Dim0-1>(),Number<Current_Dim1>(),Number<Current_Dim2>(),Number<Current_Dim3>()); \
  } \
  template<class A, class B, class U, int Dim0, int Dim1,int Dim2, int Dim3, char i, char j, char k, char l, \
    int Current_Dim1,int Current_Dim2,int Current_Dim3> \
  inline void T4_equals_T4_IJKL(A &iter, \
  			 const Tensor4_Expr<B,U,Dim0,Dim1,Dim2,Dim3,I,J,K,L> &result, \
  			 const Number<1> &, \
  			 const Number<Current_Dim1> &, \
         const Number<Current_Dim2> &, \
         const Number<Current_Dim3> & \
       ) \
  { \
    iter(0,Current_Dim1-1,Current_Dim2-1,Current_Dim3-1)=result(S20,S21,S22,S23); \
    T4_equals_T4_IJKL(iter,result,Number<Dim0>(),Number<Current_Dim1-1>(),Number<Current_Dim2>(),Number<Current_Dim3>()); \
  } \
  template<class A, class B, class U, int Dim0, int Dim1,int Dim2, int Dim3, char i, char j, char k, char l, \
    int Current_Dim2,int Current_Dim3> \
  inline void T4_equals_T4_IJKL(A &iter, \
  			 const Tensor4_Expr<B,U,Dim0,Dim1,Dim2,Dim3,I,J,K,L> &result, \
  			 const Number<1> &, \
  			 const Number<1> &, \
         const Number<Current_Dim2> &, \
         const Number<Current_Dim3> & \
       ) \
  { \
    iter(0,0,Current_Dim2-1,Current_Dim3-1)=result(S30,S31,S32,S33); \
    T4_equals_T4_IJKL(iter,result,Number<Dim0>(),Number<Dim1>(),Number<Current_Dim2-1>(),Number<Current_Dim3>()); \
  } \
  template<class A, class B, class U, int Dim0, int Dim1,int Dim2, int Dim3, char i, char j, char k, char l, \
    int Current_Dim3> \
  inline void T4_equals_T4_IJKL(A &iter, \
  			 const Tensor4_Expr<B,U,Dim0,Dim1,Dim2,Dim3,I,J,K,L> &result, \
  			 const Number<1> &, \
  			 const Number<1> &, \
         const Number<1> &, \
         const Number<Current_Dim3> & \
       ) \
  { \
    iter(0,0,0,Current_Dim3-1)=result(S40,S41,S42,S43); \
    T4_equals_T4_IJKL(iter,result,Number<Dim0>(),Number<Dim1>(),Number<Dim2>(),Number<Current_Dim3-1>()); \
  } \
  template<class A, class B, class U, int Dim0, int Dim1,int Dim2,int Dim3, char i, char j,char k,char l> \
  inline void T4_equals_T4_IJKL(A &iter, \
    const Tensor4_Expr<B,U,Dim0,Dim1,Dim2,Dim3,I,J,K,L> &result, \
    const Number<1> &, \
    const Number<1> &, \
    const Number<1> &, \
    const Number<1> & \
  ) \
  { \
    iter(0,0,0,0)=result(0,0,0,0); \
  } \
  template<class A, class T, int Dim0, int Dim1,int Dim2,int Dim3, char i, char j, char k, char l> \
  template<class B, class U> inline \
  const Tensor4_Expr<Tensor4<A,Dim0,Dim1,Dim2,Dim3>,T,Dim0,Dim1,Dim2,Dim3,i,j,k,l> & \
  Tensor4_Expr<Tensor4<A,Dim0,Dim1,Dim2,Dim3>,T,Dim0,Dim1,Dim2,Dim3,i,j,k,l>:: \
  operator=(const Tensor4_Expr<B,U,D0,D1,D2,D3,I,J,K,L> &result) \
  { \
    T4_equals_T4_IJKL(iter,result,Number<Dim0>(),Number<Dim1>(),Number<Dim2>(),Number<Dim3>()); \
    return *this; \
  }

  // jikl
  EQUAL(
    T4_equals_T4_1023,
    Dim1,Dim0,Dim2,Dim3,
    j,i,k,l,
    Current_Dim1-1,Current_Dim0-1,Current_Dim2-1,Current_Dim3-1,
    Current_Dim1-1,0,Current_Dim2-1,Current_Dim3-1,
    0,0,Current_Dim2-1,Current_Dim3-1,
    0,0,0,Current_Dim3-1
  );


  // jkil
  EQUAL(
    T4_equals_T4_1203,
    Dim1,Dim2,Dim0,Dim3,
    j,k,i,l,
    Current_Dim1-1,Current_Dim2-1,Current_Dim0-1,Current_Dim3-1,
    Current_Dim1-1,Current_Dim2-1,0,Current_Dim3-1,
    0,Current_Dim2-1,0,Current_Dim3-1,
    0,0,0,Current_Dim3-1
  );

  // jkli
  EQUAL(
    T4_equals_T4_1230,
    Dim1,Dim2,Dim3,Dim0,
    j,k,l,i,
    Current_Dim1-1,Current_Dim2-1,Current_Dim3-1,Current_Dim0-1,
    Current_Dim1-1,Current_Dim2-1,Current_Dim3-1,0,
    0,Current_Dim2-1,Current_Dim3-1,0,
    0,0,Current_Dim3-1,0
  );

  // kjli
  EQUAL(
    T4_equals_T4_2130,
    Dim2,Dim1,Dim3,Dim0,
    k,j,l,i,
    Current_Dim2-1,Current_Dim1-1,Current_Dim3-1,Current_Dim0-1,
    Current_Dim2-1,Current_Dim1-1,Current_Dim3-1,0,
    Current_Dim2-1,0,Current_Dim3-1,0,
    0,0,Current_Dim3-1,0
  );

  // klji
  EQUAL(
    T4_equals_T4_2310,
    Dim2,Dim3,Dim1,Dim0,
    k,l,j,i,
    Current_Dim2-1,Current_Dim3-1,Current_Dim1-1,Current_Dim0-1,
    Current_Dim2-1,Current_Dim3-1,Current_Dim1-1,0,
    Current_Dim2-1,Current_Dim3-1,0,0,
    0,Current_Dim3-1,0,0
  );

  // klij
  EQUAL(
    T4_equals_T4_2301,
    Dim2,Dim3,Dim0,Dim1,
    k,l,i,j,
    Current_Dim2-1,Current_Dim3-1,Current_Dim0-1,Current_Dim1-1,
    Current_Dim2-1,Current_Dim3-1,0,Current_Dim1-1,
    Current_Dim2-1,Current_Dim3-1,0,0,
    0,Current_Dim3-1,0,0
  );

  // lkij
  EQUAL(
    T4_equals_T4_3201,
    Dim3,Dim2,Dim0,Dim1,
    l,k,i,j,
    Current_Dim3-1,Current_Dim2-1,Current_Dim0-1,Current_Dim1-1,
    Current_Dim3-1,Current_Dim2-1,0,Current_Dim1-1,
    Current_Dim3-1,Current_Dim2-1,0,0,
    Current_Dim3-1,0,0,0
  );

  // likj
  EQUAL(
    T4_equals_T4_3021,
    Dim3,Dim0,Dim2,Dim1,
    l,i,k,j,
    Current_Dim3-1,Current_Dim0-1,Current_Dim2-1,Current_Dim1-1,
    Current_Dim3-1,0,Current_Dim2-1,Current_Dim1-1,
    Current_Dim3-1,0,Current_Dim2-1,0,
    Current_Dim3-1,0,0,0
  );

  // lijk
  EQUAL(
    T4_equals_T4_3012,
    Dim3,Dim0,Dim1,Dim2,
    l,i,j,k,
    Current_Dim3-1,Current_Dim0-1,Current_Dim1-1,Current_Dim2-1,
    Current_Dim3-1,0,Current_Dim1-1,Current_Dim2-1,
    Current_Dim3-1,0,0,Current_Dim2-1,
    Current_Dim3-1,0,0,0
  );

  // iljk
  EQUAL(
    T4_equals_T4_0312,
    Dim0,Dim3,Dim1,Dim2,
    i,l,j,k,
    Current_Dim0-1,Current_Dim3-1,Current_Dim1-1,Current_Dim2-1,
    0,Current_Dim3-1,Current_Dim1-1,Current_Dim2-1,
    0,Current_Dim3-1,0,Current_Dim2-1,
    0,Current_Dim3-1,0,0
  );

  // ijlk
  EQUAL(
    T4_equals_T4_0132,
    Dim0,Dim1,Dim3,Dim2,
    i,j,l,k,
    Current_Dim0-1,Current_Dim1-1,Current_Dim3-1,Current_Dim2-1,
    0,Current_Dim1-1,Current_Dim3-1,Current_Dim2-1,
    0,0,Current_Dim3-1,Current_Dim2-1,
    0,0,Current_Dim3-1,0
  );

  // lkji
  EQUAL(
    T4_equals_T4_3210,
    Dim3,Dim2,Dim1,Dim0,
    l,k,j,i,
    Current_Dim3-1,Current_Dim2-1,Current_Dim1-1,Current_Dim0-1,
    Current_Dim3-1,Current_Dim2-1,Current_Dim1-1,0,
    Current_Dim3-1,Current_Dim2-1,0,0,
    Current_Dim3-1,0,0,0
  );

  // ikjl
  EQUAL(
    T4_equals_T4_0213,
    Dim0,Dim2,Dim1,Dim3,
    i,k,j,l,
    Current_Dim0-1,Current_Dim2-1,Current_Dim1-1,Current_Dim3-1,
    0,Current_Dim2-1,Current_Dim1-1,Current_Dim3-1,
    0,Current_Dim2-1,0,Current_Dim3-1,
    0,0,0,Current_Dim3-1
  );

  // iklj
  EQUAL(
    T4_equals_T4_0231,
    Dim0,Dim2,Dim3,Dim1,
    i,k,l,j,
    Current_Dim0-1,Current_Dim2-1,Current_Dim3-1,Current_Dim1-1,
    0,Current_Dim2-1,Current_Dim3-1,Current_Dim1-1,
    0,Current_Dim2-1,Current_Dim3-1,0,
    0,0,Current_Dim3-1,0
  );

  // ilkj
  EQUAL(
    T4_equals_T4_0321,
    Dim0,Dim3,Dim2,Dim1,
    i,l,k,j,
    Current_Dim0-1,Current_Dim3-1,Current_Dim2-1,Current_Dim1-1,
    0,Current_Dim3-1,Current_Dim2-1,Current_Dim1-1,
    0,Current_Dim3-1,Current_Dim2-1,0,
    0,Current_Dim3-1,0,0
  );

  // jilk
  EQUAL(
    T4_equals_T4_1032,
    Dim1,Dim0,Dim3,Dim2,
    j,i,l,k,
    Current_Dim1-1,Current_Dim0-1,Current_Dim3-1,Current_Dim2-1,
    Current_Dim1-1,0,Current_Dim3-1,Current_Dim2-1,
    0,0,Current_Dim3-1,Current_Dim2-1,
    0,0,Current_Dim3-1,0
  );

  // kijl
  EQUAL(
    T4_equals_T4_2013,
    Dim2,Dim0,Dim1,Dim3,
    k,i,j,l,
    Current_Dim2-1,Current_Dim0-1,Current_Dim1-1,Current_Dim3-1,
    Current_Dim2-1,0,Current_Dim1-1,Current_Dim3-1,
    Current_Dim2-1,0,0,Current_Dim3-1,
    0,0,0,Current_Dim3-1
  );

  // kilj
  EQUAL(
    T4_equals_T4_2031,
    Dim2,Dim0,Dim3,Dim1,
    k,i,l,j,
    Current_Dim2-1,Current_Dim0-1,Current_Dim3-1,Current_Dim1-1,
    Current_Dim2-1,0,Current_Dim3-1,Current_Dim1-1,
    Current_Dim2-1,0,Current_Dim3-1,0,
    0,0,Current_Dim3-1,0
  );

  // jlik
  EQUAL(
    T4_equals_T4_1302,
    Dim1,Dim3,Dim0,Dim2,
    j,l,i,k,
    Current_Dim1-1,Current_Dim3-1,Current_Dim0-1,Current_Dim2-1,
    Current_Dim1-1,Current_Dim3-1,0,Current_Dim2-1,
    0,Current_Dim3-1,0,Current_Dim2-1,
    0,Current_Dim3-1,0,0
  );

  // jlki
  EQUAL(
    T4_equals_T4_1320,
    Dim1,Dim3,Dim2,Dim0,
    j,l,k,i,
    Current_Dim1-1,Current_Dim3-1,Current_Dim2-1,Current_Dim0-1,
    Current_Dim1-1,Current_Dim3-1,Current_Dim2-1,0,
    0,Current_Dim3-1,Current_Dim2-1,0,
    0,Current_Dim3-1,0,0
  );

  // kjil
  EQUAL(
    T4_equals_T4_2103,
    Dim2,Dim1,Dim0,Dim3,
    k,j,i,l,
    Current_Dim2-1,Current_Dim1-1,Current_Dim0-1,Current_Dim3-1,
    Current_Dim2-1,Current_Dim1-1,0,Current_Dim3-1,
    Current_Dim2-1,0,0,Current_Dim3-1,
    0,0,0,Current_Dim3-1
  );

  // ljik
  EQUAL(
    T4_equals_T4_3102,
    Dim3,Dim1,Dim0,Dim2,
    l,j,i,k,
    Current_Dim3-1,Current_Dim1-1,Current_Dim0-1,Current_Dim2-1,
    Current_Dim3-1,Current_Dim1-1,0,Current_Dim2-1,
    Current_Dim3-1,0,0,Current_Dim2-1,
    Current_Dim3-1,0,0,0
  );

  // ljki
  EQUAL(
    T4_equals_T4_3120,
    Dim3,Dim1,Dim2,Dim0,
    l,j,k,i,
    Current_Dim3-1,Current_Dim1-1,Current_Dim2-1,Current_Dim0-1,
    Current_Dim3-1,Current_Dim1-1,Current_Dim2-1,0,
    Current_Dim3-1,0,Current_Dim2-1,0,
    Current_Dim3-1,0,0,0
  );


  #undef EQUAL

  /* T3=U */

  template<class A, class U, int Dim0, int Dim1, int Dim2, int Dim3,
    int Current_Dim0, int Current_Dim1, int Current_Dim2, int Current_Dim3>
  inline void T4_equals_generic(A &iter, const U &u,
  				const Number<Current_Dim0> &,
  				const Number<Current_Dim1> &,
  				const Number<Current_Dim2> &,
          const Number<Current_Dim3> &,
  				const Number<Dim0> &,
  				const Number<Dim1> &,
  				const Number<Dim2> &,
          const Number<Dim3> &)
  {
    iter(Current_Dim0-1,Current_Dim1-1,Current_Dim2-1,Current_Dim3-1)=u;
    T4_equals_generic(
      iter,u,
      Number<Current_Dim0-1>(),Number<Current_Dim1>(),
      Number<Current_Dim2>(),Number<Current_Dim3>(),
  		Number<Dim0>(),Number<Dim1>(),Number<Dim2>(),Number<Dim3>()
    );
  }

  template<class A, class U, int Dim0, int Dim1, int Dim2, int Dim3,
    int Current_Dim1, int Current_Dim2,int Current_Dim3>
  inline void T4_equals_generic(A &iter, const U &u,
  				const Number<1> &,
  				const Number<Current_Dim1> &,
  				const Number<Current_Dim2> &,
          const Number<Current_Dim3> &,
  				const Number<Dim0> &,
  				const Number<Dim1> &,
  				const Number<Dim2> &,
          const Number<Dim3> &)
  {
    iter(0,Current_Dim1-1,Current_Dim2-1,Current_Dim3-1)=u;
    T4_equals_generic(
      iter,u,
      Number<Dim0>(),Number<Current_Dim1-1>(),
      Number<Current_Dim2>(),Number<Current_Dim3>(),
  		Number<Dim0>(),Number<Dim1>(),Number<Dim2>(),Number<Dim3>()
    );
  }

  template<class A, class U, int Dim0, int Dim1, int Dim2,int Dim3,int Current_Dim2,int Current_Dim3>
  inline void T4_equals_generic(A &iter, const U &u,
  				const Number<1> &, const Number<1> &,
  				const Number<Current_Dim2> &, const Number<Current_Dim3> &,
  				const Number<Dim0> &,
  				const Number<Dim1> &,
  				const Number<Dim2> &,
          const Number<Dim3> &)
  {
    iter(0,0,Current_Dim2-1,Current_Dim3-1)=u;
    T4_equals_generic(iter,u,
        Number<Dim0>(),Number<Dim1>(),
        Number<Current_Dim2-1>(),Number<Current_Dim3>(),
  		  Number<Dim0>(),Number<Dim1>(),Number<Dim2>(),Number<Dim3>()
      );
  }

  template<class A, class U, int Dim0, int Dim1, int Dim2,int Dim3 ,int Current_Dim3>
  inline void T4_equals_generic(A &iter, const U &u,
          const Number<1> &, const Number<1> &,
          const Number<1> &, const Number<Current_Dim3> &,
          const Number<Dim0> &,
          const Number<Dim1> &,
          const Number<Dim2> &,
          const Number<Dim3> &)
  {
    iter(0,0,0,Current_Dim3-1)=u;
    T4_equals_generic(iter,u,
            Number<Dim0>(),Number<Dim1>(),
            Number<Dim2>(),Number<Current_Dim3-1>(),
            Number<Dim0>(),Number<Dim1>(),
            Number<Dim2>(),Number<Dim3>());
  }

  template<class A, class U, int Dim0, int Dim1, int Dim2,int Dim3>
  inline void T4_equals_generic(A &iter, const U &u,
  				const Number<1> &, const Number<1> &,
  				const Number<1> &, const Number<1> &,
  				const Number<Dim0> &,
  				const Number<Dim1> &,
  				const Number<Dim2> &,
          const Number<Dim3> &)
  {
    iter(0,0,0,0)=u;
  }

  template<class A, class T, int Dim0, int Dim1,int Dim2,int Dim3, char i, char j, char k, char l>
  template <class U> inline
  const Tensor4_Expr<Tensor4<A,Dim0,Dim1,Dim2,Dim3>,T,Dim0,Dim1,Dim2,Dim3,i,j,k,l> &
  Tensor4_Expr<Tensor4<A,Dim0,Dim1,Dim2,Dim3>,T,Dim0,Dim1,Dim2,Dim3,i,j,k,l>::
  operator=(const U &u)
  {
    T4_equals_generic(
      iter,u,
      Number<Dim0>(),Number<Dim1>(),Number<Dim2>(),Number<Dim3>(),
  		Number<Dim0>(),Number<Dim1>(),Number<Dim2>(),Number<Dim3>()
    );
    return *this;
  }

  /* T3*=U */

  template<class A, class U, int Dim0, int Dim1, int Dim2, int Dim3,
    int Current_Dim0, int Current_Dim1, int Current_Dim2, int Current_Dim3>
  inline void T4_times_equals_generic(A &iter, const U &u,
          const Number<Current_Dim0> &,
          const Number<Current_Dim1> &,
          const Number<Current_Dim2> &,
          const Number<Current_Dim3> &,
          const Number<Dim0> &,
          const Number<Dim1> &,
          const Number<Dim2> &,
          const Number<Dim3> &)
  {
    iter(Current_Dim0-1,Current_Dim1-1,Current_Dim2-1,Current_Dim3-1)*=u;
    T4_times_equals_generic(
      iter,u,
      Number<Current_Dim0-1>(),Number<Current_Dim1>(),
      Number<Current_Dim2>(),Number<Current_Dim3>(),
      Number<Dim0>(),Number<Dim1>(),Number<Dim2>(),Number<Dim3>()
    );
  }

  template<class A, class U, int Dim0, int Dim1, int Dim2, int Dim3,
    int Current_Dim1, int Current_Dim2,int Current_Dim3>
  inline void T4_times_equals_generic(A &iter, const U &u,
          const Number<1> &,
          const Number<Current_Dim1> &,
          const Number<Current_Dim2> &,
          const Number<Current_Dim3> &,
          const Number<Dim0> &,
          const Number<Dim1> &,
          const Number<Dim2> &,
          const Number<Dim3> &)
  {
    iter(0,Current_Dim1-1,Current_Dim2-1,Current_Dim3-1)*=u;
    T4_times_equals_generic(
      iter,u,
      Number<Dim0>(),Number<Current_Dim1-1>(),
      Number<Current_Dim2>(),Number<Current_Dim3>(),
      Number<Dim0>(),Number<Dim1>(),Number<Dim2>(),Number<Dim3>()
    );
  }

  template<class A, class U, int Dim0, int Dim1, int Dim2,int Dim3,int Current_Dim2,int Current_Dim3>
  inline void T4_times_equals_generic(A &iter, const U &u,
          const Number<1> &, const Number<1> &,
          const Number<Current_Dim2> &, const Number<Current_Dim3> &,
          const Number<Dim0> &,
          const Number<Dim1> &,
          const Number<Dim2> &,
          const Number<Dim3> &)
  {
    iter(0,0,Current_Dim2-1,Current_Dim3-1)*=u;
    T4_times_equals_generic(iter,u,
        Number<Dim0>(),Number<Dim1>(),
        Number<Current_Dim2-1>(),Number<Current_Dim3>(),
        Number<Dim0>(),Number<Dim1>(),Number<Dim2>(),Number<Dim3>()
      );
  }

  template<class A, class U, int Dim0, int Dim1, int Dim2,int Dim3 ,int Current_Dim3>
  inline void T4_times_equals_generic(A &iter, const U &u,
          const Number<1> &, const Number<1> &,
          const Number<1> &, const Number<Current_Dim3> &,
          const Number<Dim0> &,
          const Number<Dim1> &,
          const Number<Dim2> &,
          const Number<Dim3> &)
  {
    iter(0,0,0,Current_Dim3-1)*=u;
    T4_times_equals_generic(iter,u,
            Number<Dim0>(),Number<Dim1>(),
            Number<Dim2>(),Number<Current_Dim3-1>(),
            Number<Dim0>(),Number<Dim1>(),
            Number<Dim2>(),Number<Dim3>());
  }

  template<class A, class U, int Dim0, int Dim1, int Dim2,int Dim3>
  inline void T4_times_equals_generic(A &iter, const U &u,
          const Number<1> &, const Number<1> &,
          const Number<1> &, const Number<1> &,
          const Number<Dim0> &,
          const Number<Dim1> &,
          const Number<Dim2> &,
          const Number<Dim3> &)
  {
    iter(0,0,0,0)*=u;
  }

  template<class A, class T, int Dim0, int Dim1,int Dim2,int Dim3, char i, char j, char k, char l>
  template <class U> inline
  const Tensor4_Expr<Tensor4<A,Dim0,Dim1,Dim2,Dim3>,T,Dim0,Dim1,Dim2,Dim3,i,j,k,l> &
  Tensor4_Expr<Tensor4<A,Dim0,Dim1,Dim2,Dim3>,T,Dim0,Dim1,Dim2,Dim3,i,j,k,l>::
  operator*=(const U &u)
  {
    T4_times_equals_generic(
      iter,u,
      Number<Dim0>(),Number<Dim1>(),Number<Dim2>(),Number<Dim3>(),
      Number<Dim0>(),Number<Dim1>(),Number<Dim2>(),Number<Dim3>()
    );
    return *this;
  }
