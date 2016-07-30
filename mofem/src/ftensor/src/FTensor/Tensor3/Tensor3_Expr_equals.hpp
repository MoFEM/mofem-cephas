/* Various assignment operators, */

/* T3=T3 */

template<class A, class B, class U, int Dim0, int Dim1, int Dim2,
  char i, char j, char k, int Current_Dim0, int Current_Dim1, int Current_Dim2>
inline void T3_equals_T3
(A &iter, const Tensor3_Expr<B,U,Dim0,Dim1,Dim2,i,j,k> &result,
 const Number<Current_Dim0> &, const Number<Current_Dim1> &,
 const Number<Current_Dim2> &)
{
  iter(Current_Dim0-1,Current_Dim1-1,Current_Dim2-1)=
    result(Current_Dim0-1,Current_Dim1-1,Current_Dim2-1);
  T3_equals_T3(iter,result,Number<Current_Dim0-1>(),
		   Number<Current_Dim1>(),Number<Current_Dim2>());
}

template<class A, class B, class U, int Dim0, int Dim1, int Dim2,
  char i, char j, char k, int Current_Dim1, int Current_Dim2>
inline void T3_equals_T3
(A &iter, const Tensor3_Expr<B,U,Dim0,Dim1,Dim2,i,j,k> &result,
 const Number<1> &, const Number<Current_Dim1> &,
 const Number<Current_Dim2> &)
{
  iter(0,Current_Dim1-1,Current_Dim2-1)=
    result(0,Current_Dim1-1,Current_Dim2-1);
  T3_equals_T3(iter,result,Number<Dim0>(),
		   Number<Current_Dim1-1>(),Number<Current_Dim2>());
}

template<class A, class B, class U, int Dim0, int Dim1, int Dim2,
  char i, char j, char k, int Current_Dim2>
inline void T3_equals_T3
(A &iter, const Tensor3_Expr<B,U,Dim0,Dim1,Dim2,i,j,k> &result,
 const Number<1> &, const Number<1> &,
 const Number<Current_Dim2> &)
{
  iter(0,0,Current_Dim2-1)=result(0,0,Current_Dim2-1);
  T3_equals_T3(iter,result,Number<Dim0>(),
		   Number<Dim1>(),Number<Current_Dim2-1>());
}

template<class A, class B, class U, int Dim0, int Dim1, int Dim2,
  char i, char j, char k>
inline void T3_equals_T3
(A &iter, const Tensor3_Expr<B,U,Dim0,Dim1,Dim2,i,j,k> &result,
 const Number<1> &, const Number<1> &, const Number<1> &)
{
  iter(0,0,0)=result(0,0,0);
}

template<class A, class T, int Tensor_Dim0, int Tensor_Dim1, int Tensor_Dim2,
  int Dim0, int Dim1, int Dim2, char i, char j, char k>
template<class B, class U> inline
const Tensor3_Expr<Tensor3<A,Tensor_Dim0,Tensor_Dim1,Tensor_Dim2>,T,
                   Dim0,Dim1,Dim2,i,j,k> &
Tensor3_Expr<Tensor3<A,Tensor_Dim0,Tensor_Dim1,Tensor_Dim2>,T,
             Dim0,Dim1,Dim2,i,j,k>::
operator=(const Tensor3_Expr<B,U,Dim0,Dim1,Dim2,i,j,k> &result)
{
  T3_equals_T3(iter,result,Number<Dim0>(),Number<Dim1>(),Number<Dim2>());
  return *this;
}

/* T3=T3_Expr(T3): I have to explicitly declare this operator= because
   otherwise the compiler will generate its own and not use the
   template code. */

template<class A, class T, int Tensor_Dim0,int Tensor_Dim1, int Tensor_Dim2,
  int Dim0, int Dim1, int Dim2, char i, char j, char k> inline
const Tensor3_Expr<Tensor3<A,Tensor_Dim0,Tensor_Dim1,Tensor_Dim2>,T,
                   Dim0,Dim1,Dim2,i,j,k> &
Tensor3_Expr<Tensor3<A,Tensor_Dim0,Tensor_Dim1,Tensor_Dim2>,T,
             Dim0,Dim1,Dim2,i,j,k>::
operator=(const Tensor3_Expr<Tensor3<A,Tensor_Dim0,Tensor_Dim1,Tensor_Dim2>,T,
          Dim0,Dim1,Dim2,i,j,k> &result)
{
  return operator=<Tensor3<A,Tensor_Dim0,Tensor_Dim1,Tensor_Dim2>,T>(result);
}

/* T3(i,j,k)=T3(i,k,j) */

template<class A, class B, class U, int Dim0, int Dim1, int Dim2,
  char i, char j, char k, int Current_Dim0, int Current_Dim1, int Current_Dim2>
inline void T3_equals_T3_021
(A &iter, const Tensor3_Expr<B,U,Dim0,Dim2,Dim1,i,j,k> &result,
 const Number<Current_Dim0> &, const Number<Current_Dim1> &,
 const Number<Current_Dim2> &)
{
  iter(Current_Dim0-1,Current_Dim1-1,Current_Dim2-1)=
    result(Current_Dim0-1,Current_Dim2-1,Current_Dim1-1);
  T3_equals_T3_021(iter,result,Number<Current_Dim0-1>(),
		   Number<Current_Dim1>(),Number<Current_Dim2>());
}

template<class A, class B, class U, int Dim0, int Dim1, int Dim2,
  char i, char j, char k, int Current_Dim1, int Current_Dim2>
inline void T3_equals_T3_021
(A &iter, const Tensor3_Expr<B,U,Dim0,Dim2,Dim1,i,j,k> &result,
 const Number<1> &, const Number<Current_Dim1> &,
 const Number<Current_Dim2> &)
{
  iter(0,Current_Dim1-1,Current_Dim2-1)=
    result(0,Current_Dim2-1,Current_Dim1-1);
  T3_equals_T3_021(iter,result,Number<Dim0>(),
		   Number<Current_Dim1-1>(),Number<Current_Dim2>());
}

template<class A, class B, class U, int Dim0, int Dim1, int Dim2,
  char i, char j, char k, int Current_Dim2>
inline void T3_equals_T3_021
(A &iter, const Tensor3_Expr<B,U,Dim0,Dim2,Dim1,i,j,k> &result,
 const Number<1> &, const Number<1> &,
 const Number<Current_Dim2> &)
{
  iter(0,0,Current_Dim2-1)=result(0,Current_Dim2-1,0);
  T3_equals_T3_021(iter,result,Number<Dim0>(),
		   Number<Dim1>(),Number<Current_Dim2-1>());
}

template<class A, class B, class U, int Dim0, int Dim1, int Dim2,
  char i, char j, char k>
inline void T3_equals_T3_021
(A &iter, const Tensor3_Expr<B,U,Dim0,Dim2,Dim1,i,j,k> &result,
 const Number<1> &, const Number<1> &, const Number<1> &)
{
  iter(0,0,0)=result(0,0,0);
}

template<class A, class T, int Tensor_Dim0, int Tensor_Dim1, int Tensor_Dim2,
  int Dim0, int Dim1, int Dim2, char i, char j, char k>
template<class B, class U> inline
const Tensor3_Expr<Tensor3<A,Tensor_Dim0,Tensor_Dim1,Tensor_Dim2>,T,
                   Dim0,Dim1,Dim2,i,j,k> &
Tensor3_Expr<Tensor3<A,Tensor_Dim0,Tensor_Dim1,Tensor_Dim2>,T,
             Dim0,Dim1,Dim2,i,j,k>::
operator=(const Tensor3_Expr<B,U,Dim0,Dim2,Dim1,i,k,j> &result)
{
  T3_equals_T3_021(iter,result,Number<Dim0>(),Number<Dim1>(),Number<Dim2>());
  return *this;
}

/* T3(i,j,k)=T3(j,i,k) */

template<class A, class B, class U, int Dim0, int Dim1, int Dim2,
  char i, char j, char k, int Current_Dim0, int Current_Dim1, int Current_Dim2>
inline void T3_equals_T3_102
(A &iter, const Tensor3_Expr<B,U,Dim1,Dim0,Dim2,i,j,k> &result,
 const Number<Current_Dim0> &, const Number<Current_Dim1> &,
 const Number<Current_Dim2> &)
{
  iter(Current_Dim0-1,Current_Dim1-1,Current_Dim2-1)=
    result(Current_Dim1-1,Current_Dim0-1,Current_Dim2-1);
  T3_equals_T3_102(iter,result,Number<Current_Dim0-1>(),
		   Number<Current_Dim1>(),Number<Current_Dim2>());
}

template<class A, class B, class U, int Dim0, int Dim1, int Dim2,
  char i, char j, char k, int Current_Dim1, int Current_Dim2>
inline void T3_equals_T3_102
(A &iter, const Tensor3_Expr<B,U,Dim1,Dim0,Dim2,i,j,k> &result,
 const Number<1> &, const Number<Current_Dim1> &,
 const Number<Current_Dim2> &)
{
  iter(0,Current_Dim1-1,Current_Dim2-1)=
    result(Current_Dim1-1,0,Current_Dim2-1);
  T3_equals_T3_102(iter,result,Number<Dim0>(),
		   Number<Current_Dim1-1>(),Number<Current_Dim2>());
}

template<class A, class B, class U, int Dim0, int Dim1, int Dim2,
  char i, char j, char k, int Current_Dim2>
inline void T3_equals_T3_102
(A &iter, const Tensor3_Expr<B,U,Dim1,Dim0,Dim2,i,j,k> &result,
 const Number<1> &, const Number<1> &,
 const Number<Current_Dim2> &)
{
  iter(0,0,Current_Dim2-1)=result(0,0,Current_Dim2-1);
  T3_equals_T3_102(iter,result,Number<Dim0>(),
		   Number<Dim1>(),Number<Current_Dim2-1>());
}

template<class A, class B, class U, int Dim0, int Dim1, int Dim2,
  char i, char j, char k>
inline void T3_equals_T3_102
(A &iter, const Tensor3_Expr<B,U,Dim1,Dim0,Dim2,i,j,k> &result,
 const Number<1> &, const Number<1> &, const Number<1> &)
{
  iter(0,0,0)=result(0,0,0);
}

template<class A, class T, int Tensor_Dim0, int Tensor_Dim1, int Tensor_Dim2,
  int Dim0, int Dim1, int Dim2, char i, char j, char k>
template<class B, class U> inline
const Tensor3_Expr<Tensor3<A,Tensor_Dim0,Tensor_Dim1,Tensor_Dim2>,T,
                   Dim0,Dim1,Dim2,i,j,k> &
Tensor3_Expr<Tensor3<A,Tensor_Dim0,Tensor_Dim1,Tensor_Dim2>,T,
             Dim0,Dim1,Dim2,i,j,k>::
operator=(const Tensor3_Expr<B,U,Dim1,Dim0,Dim2,j,i,k> &result)
{
  T3_equals_T3_102(iter,result,Number<Dim0>(),Number<Dim1>(),Number<Dim2>());
  return *this;
}

/* T3(i,j,k)=T3(j,k,i) */

template<class A, class B, class U, int Dim0, int Dim1, int Dim2,
  char i, char j, char k, int Current_Dim0, int Current_Dim1, int Current_Dim2>
inline void T3_equals_T3_120
(A &iter, const Tensor3_Expr<B,U,Dim1,Dim2,Dim0,i,j,k> &result,
 const Number<Current_Dim0> &, const Number<Current_Dim1> &,
 const Number<Current_Dim2> &)
{
  iter(Current_Dim0-1,Current_Dim1-1,Current_Dim2-1)=
    result(Current_Dim1-1,Current_Dim2-1,Current_Dim0-1);
  T3_equals_T3_120(iter,result,Number<Current_Dim0-1>(),
		   Number<Current_Dim1>(),Number<Current_Dim2>());
}

template<class A, class B, class U, int Dim0, int Dim1, int Dim2,
  char i, char j, char k, int Current_Dim1, int Current_Dim2>
inline void T3_equals_T3_120
(A &iter, const Tensor3_Expr<B,U,Dim1,Dim2,Dim0,i,j,k> &result,
 const Number<1> &, const Number<Current_Dim1> &,
 const Number<Current_Dim2> &)
{
  iter(0,Current_Dim1-1,Current_Dim2-1)=
    result(Current_Dim1-1,Current_Dim2-1,0);
  T3_equals_T3_120(iter,result,Number<Dim0>(),
		   Number<Current_Dim1-1>(),Number<Current_Dim2>());
}

template<class A, class B, class U, int Dim0, int Dim1, int Dim2,
  char i, char j, char k, int Current_Dim2>
inline void T3_equals_T3_120
(A &iter, const Tensor3_Expr<B,U,Dim1,Dim2,Dim0,i,j,k> &result,
 const Number<1> &, const Number<1> &,
 const Number<Current_Dim2> &)
{
  iter(0,0,Current_Dim2-1)=result(0,Current_Dim2-1,0);
  T3_equals_T3_120(iter,result,Number<Dim0>(),
		   Number<Dim1>(),Number<Current_Dim2-1>());
}

template<class A, class B, class U, int Dim0, int Dim1, int Dim2,
  char i, char j, char k>
inline void T3_equals_T3_120
(A &iter, const Tensor3_Expr<B,U,Dim1,Dim2,Dim0,i,j,k> &result,
 const Number<1> &, const Number<1> &, const Number<1> &)
{
  iter(0,0,0)=result(0,0,0);
}

template<class A, class T, int Tensor_Dim0, int Tensor_Dim1, int Tensor_Dim2,
  int Dim0, int Dim1, int Dim2, char i, char j, char k>
template<class B, class U> inline
const Tensor3_Expr<Tensor3<A,Tensor_Dim0,Tensor_Dim1,Tensor_Dim2>,T,
                   Dim0,Dim1,Dim2,i,j,k> &
Tensor3_Expr<Tensor3<A,Tensor_Dim0,Tensor_Dim1,Tensor_Dim2>,T,
             Dim0,Dim1,Dim2,i,j,k>::
operator=(const Tensor3_Expr<B,U,Dim1,Dim2,Dim0,j,k,i> &result)
{
  T3_equals_T3_120(iter,result,Number<Dim0>(),Number<Dim1>(),Number<Dim2>());
  return *this;
}

/* T3(i,j,k)=T3(k,i,j) */

template<class A, class B, class U, int Dim0, int Dim1, int Dim2,
  char i, char j, char k, int Current_Dim0, int Current_Dim1, int Current_Dim2>
inline void T3_equals_T3_201
(A &iter, const Tensor3_Expr<B,U,Dim2,Dim0,Dim1,i,j,k> &result,
 const Number<Current_Dim0> &, const Number<Current_Dim1> &,
 const Number<Current_Dim2> &)
{
  iter(Current_Dim0-1,Current_Dim1-1,Current_Dim2-1)=
    result(Current_Dim2-1,Current_Dim0-1,Current_Dim1-1);
  T3_equals_T3_201(iter,result,Number<Current_Dim0-1>(),
		   Number<Current_Dim1>(),Number<Current_Dim2>());
}

template<class A, class B, class U, int Dim0, int Dim1, int Dim2,
  char i, char j, char k, int Current_Dim1, int Current_Dim2>
inline void T3_equals_T3_201
(A &iter, const Tensor3_Expr<B,U,Dim2,Dim0,Dim1,i,j,k> &result,
 const Number<1> &, const Number<Current_Dim1> &,
 const Number<Current_Dim2> &)
{
  iter(0,Current_Dim1-1,Current_Dim2-1)=
    result(Current_Dim2-1,0,Current_Dim1-1);
  T3_equals_T3_201(iter,result,Number<Dim0>(),
		   Number<Current_Dim1-1>(),Number<Current_Dim2>());
}

template<class A, class B, class U, int Dim0, int Dim1, int Dim2,
  char i, char j, char k, int Current_Dim2>
inline void T3_equals_T3_201
(A &iter, const Tensor3_Expr<B,U,Dim2,Dim0,Dim1,i,j,k> &result,
 const Number<1> &, const Number<1> &,
 const Number<Current_Dim2> &)
{
  iter(0,0,Current_Dim2-1)=result(Current_Dim2-1,0,0);
  T3_equals_T3_201(iter,result,Number<Dim0>(),
		   Number<Dim1>(),Number<Current_Dim2-1>());
}

template<class A, class B, class U, int Dim0, int Dim1, int Dim2,
  char i, char j, char k>
inline void T3_equals_T3_201
(A &iter, const Tensor3_Expr<B,U,Dim2,Dim0,Dim1,i,j,k> &result,
 const Number<1> &, const Number<1> &, const Number<1> &)
{
  iter(0,0,0)=result(0,0,0);
}

template<class A, class T, int Tensor_Dim0, int Tensor_Dim1, int Tensor_Dim2,
  int Dim0, int Dim1, int Dim2, char i, char j, char k>
template<class B, class U> inline
const Tensor3_Expr<Tensor3<A,Tensor_Dim0,Tensor_Dim1,Tensor_Dim2>,T,
                   Dim0,Dim1,Dim2,i,j,k> &
Tensor3_Expr<Tensor3<A,Tensor_Dim0,Tensor_Dim1,Tensor_Dim2>,T,
             Dim0,Dim1,Dim2,i,j,k>::
operator=(const Tensor3_Expr<B,U,Dim2,Dim0,Dim1,k,i,j> &result)
{
  T3_equals_T3_201(iter,result,Number<Dim0>(),Number<Dim1>(),Number<Dim2>());
  return *this;
}

/* T3(i,j,k)=T3(k,j,i) */

template<class A, class B, class U, int Dim0, int Dim1, int Dim2,
  char i, char j, char k, int Current_Dim0, int Current_Dim1, int Current_Dim2>
inline void T3_equals_T3_210
(A &iter, const Tensor3_Expr<B,U,Dim2,Dim1,Dim0,i,j,k> &result,
 const Number<Current_Dim0> &, const Number<Current_Dim1> &,
 const Number<Current_Dim2> &)
{
  iter(Current_Dim0-1,Current_Dim1-1,Current_Dim2-1)=
    result(Current_Dim2-1,Current_Dim1-1,Current_Dim0-1);
  T3_equals_T3_210(iter,result,Number<Current_Dim0-1>(),
		   Number<Current_Dim1>(),Number<Current_Dim2>());
}

template<class A, class B, class U, int Dim0, int Dim1, int Dim2,
  char i, char j, char k, int Current_Dim1, int Current_Dim2>
inline void T3_equals_T3_210
(A &iter, const Tensor3_Expr<B,U,Dim2,Dim1,Dim0,i,j,k> &result,
 const Number<1> &, const Number<Current_Dim1> &,
 const Number<Current_Dim2> &)
{
  iter(0,Current_Dim1-1,Current_Dim2-1)=
    result(Current_Dim2-1,Current_Dim1-1,0);
  T3_equals_T3_210(iter,result,Number<Dim0>(),
		   Number<Current_Dim1-1>(),Number<Current_Dim2>());
}

template<class A, class B, class U, int Dim0, int Dim1, int Dim2,
  char i, char j, char k, int Current_Dim2>
inline void T3_equals_T3_210
(A &iter, const Tensor3_Expr<B,U,Dim2,Dim1,Dim0,i,j,k> &result,
 const Number<1> &, const Number<1> &,
 const Number<Current_Dim2> &)
{
  iter(0,0,Current_Dim2-1)=result(Current_Dim2-1,0,0);
  T3_equals_T3_210(iter,result,Number<Dim0>(),
		   Number<Dim1>(),Number<Current_Dim2-1>());
}

template<class A, class B, class U, int Dim0, int Dim1, int Dim2,
  char i, char j, char k>
inline void T3_equals_T3_210
(A &iter, const Tensor3_Expr<B,U,Dim2,Dim1,Dim0,i,j,k> &result,
 const Number<1> &, const Number<1> &, const Number<1> &)
{
  iter(0,0,0)=result(0,0,0);
}

template<class A, class T, int Tensor_Dim0, int Tensor_Dim1, int Tensor_Dim2,
  int Dim0, int Dim1, int Dim2, char i, char j, char k>
template<class B, class U> inline
const Tensor3_Expr<Tensor3<A,Tensor_Dim0,Tensor_Dim1,Tensor_Dim2>,T,
                   Dim0,Dim1,Dim2,i,j,k> &
Tensor3_Expr<Tensor3<A,Tensor_Dim0,Tensor_Dim1,Tensor_Dim2>,T,
             Dim0,Dim1,Dim2,i,j,k>::
operator=(const Tensor3_Expr<B,U,Dim2,Dim1,Dim0,k,j,i> &result)
{
  T3_equals_T3_210(iter,result,Number<Dim0>(),Number<Dim1>(),Number<Dim2>());
  return *this;
}

/* T3=U */

template<class A, class U, int Dim0, int Dim1, int Dim2,
  int Current_Dim0, int Current_Dim1, int Current_Dim2>
inline void T3_equals_generic(A &iter, const U &u,
				const Number<Current_Dim0> &,
				const Number<Current_Dim1> &,
				const Number<Current_Dim2> &,
				const Number<Dim0> &,
				const Number<Dim1> &,
				const Number<Dim2> &)
{
  iter(Current_Dim0-1,Current_Dim1-1,Current_Dim2-1)=u;
  T3_equals_generic(iter,u,Number<Current_Dim0-1>(),
		      Number<Current_Dim1>(),Number<Current_Dim2>(),
		      Number<Dim0>(),Number<Dim1>(),Number<Dim2>());
}

template<class A, class U, int Dim0, int Dim1, int Dim2,
  int Current_Dim1, int Current_Dim2>
inline void T3_equals_generic(A &iter, const U &u,
				const Number<1> &,
				const Number<Current_Dim1> &,
				const Number<Current_Dim2> &,
				const Number<Dim0> &,
				const Number<Dim1> &,
				const Number<Dim2> &)
{
  iter(0,Current_Dim1-1,Current_Dim2-1)=u;
  T3_equals_generic(iter,u,Number<Dim0>(),
		      Number<Current_Dim1-1>(),Number<Current_Dim2>(),
		      Number<Dim0>(),Number<Dim1>(),Number<Dim2>());
}

template<class A, class U, int Dim0, int Dim1, int Dim2, int Current_Dim2>
inline void T3_equals_generic(A &iter, const U &u,
				const Number<1> &, const Number<1> &,
				const Number<Current_Dim2> &,
				const Number<Dim0> &,
				const Number<Dim1> &,
				const Number<Dim2> &)
{
  iter(0,0,Current_Dim2-1)=u;
  T3_equals_generic(iter,u,Number<Dim0>(),
		      Number<Dim1>(),Number<Current_Dim2-1>(),
		      Number<Dim0>(),Number<Dim1>(),Number<Dim2>());
}

template<class A, class U, int Dim0, int Dim1, int Dim2>
inline void T3_equals_generic(A &iter, const U &u,
				const Number<1> &, const Number<1> &,
				const Number<1> &,
				const Number<Dim0> &,
				const Number<Dim1> &,
				const Number<Dim2> &)
{
  iter(0,0,0)=u;
}

template<class A, class T, int Tensor_Dim0, int Tensor_Dim1, int Tensor_Dim2,
  int Dim0, int Dim1, int Dim2, char i, char j, char k>
template<class U> inline
const Tensor3_Expr<Tensor3<A,Tensor_Dim0,Tensor_Dim1,Tensor_Dim2>,T,
                   Dim0,Dim1,Dim2,i,j,k> &
Tensor3_Expr<Tensor3<A,Tensor_Dim0,Tensor_Dim1,Tensor_Dim2>,T,
             Dim0,Dim1,Dim2,i,j,k>::
operator=(const U &u)
{
  T3_equals_generic(iter,u,Number<Dim0>(),Number<Dim1>(),Number<Dim2>(),
		      Number<Dim0>(),Number<Dim1>(),Number<Dim2>());
  return *this;
}

