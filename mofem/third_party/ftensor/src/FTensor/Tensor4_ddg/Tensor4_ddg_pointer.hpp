/* A version for pointers. */

template<
  class T,
  int Dim01,int Dim23,
  int Current_Position0,int Current_Position1
>
inline void T4_ddg_increment(
  const Tensor4_ddg<T,Dim01,Dim23> &iter,
  const Number<Current_Position0> &,
  const Number<Current_Position1> &
) {
  iter.increment(
    Number<Current_Position0>(),Number<Current_Position1>()
  );
  T4_ddg_increment(
    iter,Number<Current_Position0-1>(),Number<Current_Position1>()
  );
}

template<class T,int Dim01,int Dim23,int Current_Position1>
inline void T4_ddg_increment(
  const Tensor4_ddg<T,Dim01,Dim23> &iter,
  const Number<1> &,
  const Number<Current_Position1> &
) {
	iter.increment(
		Number<1>(),Number<Current_Position1>()
	);
	T4_ddg_increment(
		iter,Number<(Dim01*(Dim01+1))/2>(),Number<Current_Position1-1>()
	);
}

template<class T,int Dim01,int Dim23>
inline void T4_ddg_increment(
  const Tensor4_ddg<T,Dim01,Dim23> &iter,
  const Number<1> &,
  const Number<1> &
) {
	iter.increment(
		Number<1>(),Number<1>()
	);
}

template <class T, int Tensor_Dim01, int Tensor_Dim23>
class Tensor4_ddg<T*,Tensor_Dim01,Tensor_Dim23>
{
  const int inc;
  mutable T * restrict data[(Tensor_Dim01*(Tensor_Dim01+1))/2][(Tensor_Dim23*(Tensor_Dim23+1))/2];
public:
  /* There are two operator(int,int,int,int)'s, one for non-consts
     that lets you change the value, and one for consts that
     doesn't. */

  Tensor4_ddg(
    T* d0000, T* d0010, T* d0011,
    T* d1000, T* d1010, T* d1011,
    T* d1100, T* d1110, T* d1111,
    const int i = 1
  ):
  inc(i) {
    ptr(0,0,0,0) = d0000; ptr(0,0,1,0) = d0010; ptr(0,0,1,1) = d0011;
    ptr(1,0,0,0) = d1000; ptr(1,0,1,0) = d1010; ptr(1,0,1,1) = d1011;
    ptr(1,1,0,0) = d1100; ptr(1,1,1,0) = d1110; ptr(1,1,1,1) = d1111;
  }

  Tensor4_ddg(
    T* d0000, T* d0001, T* d0002, T* d0011, T* d0012, T* d0022,
    T* d0100, T* d0101, T* d0102, T* d0111, T* d0112, T* d0122,
    T* d0200, T* d0201, T* d0202, T* d0211, T* d0212, T* d0222,
    T* d1100, T* d1101, T* d1102, T* d1111, T* d1112, T* d1122,
    T* d1200, T* d1201, T* d1202, T* d1211, T* d1212, T* d1222,
    T* d2200, T* d2201, T* d2202, T* d2211, T* d2212, T* d2222,
    const int i = 1
  ):
  inc(i) {
    ptr(0,0,0,0) = d0000; ptr(0,0,0,1) = d0001; ptr(0,0,0,2) = d0002; ptr(0,0,1,1) = d0011; ptr(0,0,1,2) = d0012; ptr(0,0,2,2) = d0022;
    ptr(0,1,0,0) = d0100; ptr(0,1,0,1) = d0101; ptr(0,1,0,2) = d0102; ptr(0,1,1,1) = d0111; ptr(0,1,1,2) = d0112; ptr(0,1,2,2) = d0122;
    ptr(0,2,0,0) = d0200; ptr(0,2,0,1) = d0201; ptr(0,2,0,2) = d0202; ptr(0,2,1,1) = d0211; ptr(0,2,1,2) = d0212; ptr(0,2,2,2) = d0222;
    ptr(1,1,0,0) = d1100; ptr(1,1,0,1) = d1101; ptr(1,1,0,2) = d1102; ptr(1,1,1,1) = d1111; ptr(1,1,1,2) = d1112; ptr(1,1,2,2) = d1122;
    ptr(1,2,0,0) = d1200; ptr(1,2,0,1) = d1201; ptr(1,2,0,2) = d1202; ptr(1,2,1,1) = d1211; ptr(1,2,1,2) = d1212; ptr(1,2,2,2) = d1222;
    ptr(2,2,0,0) = d2200; ptr(2,2,0,1) = d2201; ptr(2,2,0,2) = d2202; ptr(2,2,1,1) = d2211; ptr(2,2,1,2) = d2212; ptr(2,2,2,2) = d2222;
  }

  T & operator()(const int N1, const int N2, const int N3, const int N4)
  {
#ifdef FTENSOR_DEBUG
    if(N1>=Tensor_Dim01 || N1<0 || N2>=Tensor_Dim01 || N2<0
       || N3>=Tensor_Dim23 || N3<0 || N4>=Tensor_Dim23 || N4<0)
      {
        std::stringstream s;
        s << "Bad index in Tensor3_dg<T*,"
          << Tensor_Dim01 << "," << Tensor_Dim23
          << ">.operator("
          << N1 << "," << N2 << "," << N3 << "," << N4 << ")"
          << std::endl;
        throw std::runtime_error(s.str());
      }
#endif
    return
    N1>N2 ? (
      N3>N4
      ? *data[N1+(N2*(2*Tensor_Dim01-N2-1))/2][N3+(N4*(2*Tensor_Dim23-N4-1))/2]
      : *data[N1+(N2*(2*Tensor_Dim01-N2-1))/2][N4+(N3*(2*Tensor_Dim23-N3-1))/2]
    )
    : (
      N3>N4
      ? *data[N2+(N1*(2*Tensor_Dim01-N1-1))/2][N3+(N4*(2*Tensor_Dim23-N4-1))/2]
      : *data[N2+(N1*(2*Tensor_Dim01-N1-1))/2][N4+(N3*(2*Tensor_Dim23-N3-1))/2]
    );
  }

  T operator()(const int N1, const int N2, const int N3, const int N4)
    const
  {
#ifdef FTENSOR_DEBUG
    if(N1>=Tensor_Dim01 || N1<0 || N2>=Tensor_Dim01 || N2<0
       || N3>=Tensor_Dim23 || N3<0 || N4>=Tensor_Dim23 || N4<0)
      {
        std::stringstream s;
        s << "Bad index in Tensor3_dg<T*,"
          << Tensor_Dim01 << "," << Tensor_Dim23
          << ">.operator("
          << N1 << "," << N2 << "," << N3 << "," << N4
          << ") const"
          << std::endl;
        throw std::runtime_error(s.str());
      }
#endif
    return
    N1>N2 ?
    (
      N3>N4
      ? *data[N1+(N2*(2*Tensor_Dim01-N2-1))/2] [N3+(N4*(2*Tensor_Dim23-N4-1))/2]
		  : *data[N1+(N2*(2*Tensor_Dim01-N2-1))/2] [N4+(N3*(2*Tensor_Dim23-N3-1))/2]
    )
    : (
      N3>N4
      ? *data[N2+(N1*(2*Tensor_Dim01-N1-1))/2] [N3+(N4*(2*Tensor_Dim23-N4-1))/2]
	    : *data[N2+(N1*(2*Tensor_Dim01-N1-1))/2] [N4+(N3*(2*Tensor_Dim23-N3-1))/2]
    );
  }

  T* ptr(const int N1, const int N2, const int N3, const int N4) const
  {
#ifdef FTENSOR_DEBUG
    if(N1>=Tensor_Dim01 || N1<0 || N2>=Tensor_Dim01 || N2<0
       || N3>=Tensor_Dim23 || N3<0 || N4>=Tensor_Dim23 || N4<0)
      {
        std::stringstream s;
        s << "Bad index in Tensor3_dg<T,"
          << Tensor_Dim01 << "," << Tensor_Dim23 << ">.ptr("
          << N1 << "," << N2 << "," << N3 << "," << N4 << ")"
          << std::endl;
        throw std::runtime_error(s.str());
      }
#endif
    return
    N1>N2 ?
    (
      N3>N4
      ? data[N1+(N2*(2*Tensor_Dim01-N2-1))/2] [N3+(N4*(2*Tensor_Dim23-N4-1))/2]
      : data[N1+(N2*(2*Tensor_Dim01-N2-1))/2] [N4+(N3*(2*Tensor_Dim23-N3-1))/2]
    )
    : (
      N3>N4
      ? data[N2+(N1*(2*Tensor_Dim01-N1-1))/2] [N3+(N4*(2*Tensor_Dim23-N4-1))/2]
      : data[N2+(N1*(2*Tensor_Dim01-N1-1))/2] [N4+(N3*(2*Tensor_Dim23-N3-1))/2]
    );
  }

  /* Return reference to pointer, that will allow to set internal
  data pointer from storage structure.
  */
  T* restrict & ptr(const int N1, const int N2, const int N3, const int N4)
  {
    #ifdef FTENSOR_DEBUG
    if(
      N1>=Tensor_Dim01 || N1<0 || N2>=Tensor_Dim01 || N2<0
      || N3>=Tensor_Dim23 || N3<0 || N4>=Tensor_Dim23 || N4<0
    )
    {
      std::stringstream s;
      s << "Bad index in Tensor3_dg<T,"
      << Tensor_Dim01 << "," << Tensor_Dim23 << ">.ptr("
      << N1 << "," << N2 << "," << N3 << "," << N4 << ")"
      << std::endl;
      throw std::runtime_error(s.str());
    }
    #endif
    return
    N1>N2 ?
    (
      N3>N4
      ? data[N1+(N2*(2*Tensor_Dim01-N2-1))/2][N3+(N4*(2*Tensor_Dim23-N4-1))/2]
      : data[N1+(N2*(2*Tensor_Dim01-N2-1))/2][N4+(N3*(2*Tensor_Dim23-N3-1))/2]
    )
    : (
      N3>N4
      ? data[N2+(N1*(2*Tensor_Dim01-N1-1))/2][N3+(N4*(2*Tensor_Dim23-N4-1))/2]
      : data[N2+(N1*(2*Tensor_Dim01-N1-1))/2][N4+(N3*(2*Tensor_Dim23-N3-1))/2]
    );
  }

  /* These operator()'s are the first part in constructing template
     expressions.  They can be used to slice off lower dimensional
     parts. They are not entirely safe, since you can accidently use a
     higher dimension than what is really allowed (like Dim=5). */

  template<char i, char j, char k, char l, int Dim01, int Dim23>
  Tensor4_ddg_Expr<Tensor4_ddg<T*,Tensor_Dim01,Tensor_Dim23>,
    T,Dim01,Dim23,i,j,k,l> operator()
    (const Index<i,Dim01> index1, const Index<j,Dim01> index2,
     const Index<k,Dim23> index3, const Index<l,Dim23> index4)
  {
    return Tensor4_ddg_Expr<Tensor4_ddg<T*,Tensor_Dim01,Tensor_Dim23>,
      T,Dim01,Dim23,i,j,k,l> (*this);
  }


  template<char i, char j, char k, char l, int Dim01, int Dim23>
  const Tensor4_ddg_Expr<const Tensor4_ddg<T*,Tensor_Dim01,Tensor_Dim23>,
    T,Dim01,Dim23,i,j,k,l> operator()
    (const Index<i,Dim01> index1, const Index<j,Dim01> index2,
     const Index<k,Dim23> index3, const Index<l,Dim23> index4) const
  {
    return Tensor4_ddg_Expr<const Tensor4_ddg<T*,Tensor_Dim01,Tensor_Dim23>,
      T,Dim01,Dim23,i,j,k,l> (*this);
  }

  /* This is for expressions where a number is used for two slots, and
     an index for the other two, yielding a Tensor2_symmetric_Expr. */

  template<char i, char j, int N0, int N1, int Dim>
  const Tensor2_symmetric_Expr<const Tensor4_ddg_number_01
  <const Tensor4_ddg<T*,Tensor_Dim01,Tensor_Dim23>,T,N0,N1>,T,Dim,i,j>
  operator()(const Number<N0> n1, const Number<N1> n2,
	     const Index<i,Dim> index1, const Index<j,Dim> index2) const
  {
    typedef const Tensor4_ddg_number_01<const Tensor4_ddg
      <T*,Tensor_Dim01,Tensor_Dim23>,T,N0,N1> TensorExpr;
    return Tensor2_symmetric_Expr<TensorExpr,T,Dim,i,j>(TensorExpr(*this));
  }

  template<char i, char j, int N0, int N1, int Dim>
  Tensor2_symmetric_Expr<Tensor4_ddg_number_rhs_01
  <Tensor4_ddg<T*,Tensor_Dim01,Tensor_Dim23>,T,N0,N1>,T,Dim,i,j>
  operator()(const Number<N0> n1, const Number<N1> n2,
	     const Index<i,Dim> index1, const Index<j,Dim> index2)
  {
    typedef Tensor4_ddg_number_rhs_01<Tensor4_ddg<T*,Tensor_Dim01,Tensor_Dim23>,
      T,N0,N1> TensorExpr;
    return Tensor2_symmetric_Expr<TensorExpr,T,Dim,i,j>(*this);
  }

  /* This is for expressions where a number is used for one slot, and
     an index for the other three, yielding a Tensor3_dg_Expr. */

  template<char i, char j, char k, int N0, int Dim1, int Dim23>
  const Tensor3_dg_Expr<const Tensor4_ddg_number_0
  <const Tensor4_ddg<T*,Tensor_Dim01,Tensor_Dim23>,T,N0>,T,Dim23,Dim1,i,j,k>
  operator()(const Number<N0> n1, const Index<k,Dim1> index3,
	     const Index<i,Dim23> index1, const Index<j,Dim23> index2) const
  {
    typedef const Tensor4_ddg_number_0<const Tensor4_ddg
      <T*,Tensor_Dim01,Tensor_Dim23>,T,N0> TensorExpr;
    return Tensor3_dg_Expr<TensorExpr,T,Dim23,Dim1,i,j,k>(TensorExpr(*this));
  }

  template<char i, char j, char k, int N0, int Dim1, int Dim23>
  Tensor3_dg_Expr<Tensor4_ddg_number_rhs_0
  <Tensor4_ddg<T*,Tensor_Dim01,Tensor_Dim23>,T,N0>,T,Dim23,Dim1,i,j,k>
  operator()(const Number<N0> n1, const Index<k,Dim1> index3,
	     const Index<i,Dim23> index1, const Index<j,Dim23> index2)
  {
    typedef Tensor4_ddg_number_rhs_0<Tensor4_ddg<T*,Tensor_Dim01,Tensor_Dim23>,
      T,N0> TensorExpr;
    return Tensor3_dg_Expr<TensorExpr,T,Dim23,Dim1,i,j,k>(*this);
  }

  /* This is for expressions where an int (not a Number) is used for
     two slots, and an index for the other two, yielding a
     Tensor2_symmetric_Expr. */

  template<char i, char j, int Dim>
  const Tensor2_symmetric_Expr<const Tensor4_ddg_numeral_01
  <const Tensor4_ddg<T*,Tensor_Dim01,Tensor_Dim23>,T>,T,Dim,i,j>
  operator()(const int N0, const int N1,
	     const Index<i,Dim> index1, const Index<j,Dim> index2) const
  {
    typedef const Tensor4_ddg_numeral_01<const Tensor4_ddg
      <T*,Tensor_Dim01,Tensor_Dim23>,T> TensorExpr;
    return Tensor2_symmetric_Expr<TensorExpr,T,Dim,i,j>
      (TensorExpr(*this,N0,N1));
  }

  template<char i, char j, int Dim>
  const Tensor2_symmetric_Expr<const Tensor4_ddg_numeral_23
  <const Tensor4_ddg<T*,Tensor_Dim01,Tensor_Dim23>,T>,T,Dim,i,j>
  operator()(const Index<i,Dim> index1, const Index<j,Dim> index2,
	     const int N2, const int N3) const
  {
    typedef const Tensor4_ddg_numeral_23<const Tensor4_ddg
      <T*,Tensor_Dim01,Tensor_Dim23>,T> TensorExpr;
    return Tensor2_symmetric_Expr<TensorExpr,T,Dim,i,j>
      (TensorExpr(*this,N2,N3));
  }

  /* int in three slots, Index in the other yielding a Tensor1_Expr. */

  template<char i, int Dim>
  const Tensor1_Expr<const Tensor4_ddg_numeral_123<const Tensor4_ddg
  <T*,Tensor_Dim01,Tensor_Dim23>,T>,T,Dim,i>
  operator()(const Index<i,Dim> index1, const int N1, const int N2,
	     const int N3)
  {
    typedef const Tensor4_ddg_numeral_123<const Tensor4_ddg
      <T*,Tensor_Dim01,Tensor_Dim23>,T> TensorExpr;
    return Tensor1_Expr<TensorExpr,T,Dim,i>(TensorExpr(*this,N1,N2,N3));
  }

  template<char i, int Dim>
  const Tensor1_Expr<const Tensor4_ddg_numeral_123<const Tensor4_ddg
  <T*,Tensor_Dim01,Tensor_Dim23>,T>,T,Dim,i>
  operator()(const int N1, const Index<i,Dim> index1, const int N2,
	     const int N3)
  {
    typedef const Tensor4_ddg_numeral_123<const Tensor4_ddg
      <T*,Tensor_Dim01,Tensor_Dim23>,T> TensorExpr;
    return Tensor1_Expr<TensorExpr,T,Dim,i>(TensorExpr(*this,N1,N2,N3));
  }

  /* This is for expressions where an int (not a Number) is used for
     one slot, and an index for the other three, yielding a
     Tensor3_dg_Expr. */

  template<char i, char j, char k, int Dim1, int Dim23>
  const Tensor3_dg_Expr<const Tensor4_ddg_numeral_0
  <const Tensor4_ddg<T*,Tensor_Dim01,Tensor_Dim23>,T>,T,Dim23,Dim1,i,j,k>
  operator()(const int N0, const Index<k,Dim1> index3,
	     const Index<i,Dim23> index1, const Index<j,Dim23> index2) const
  {
    typedef const Tensor4_ddg_numeral_0<const Tensor4_ddg
      <T*,Tensor_Dim01,Tensor_Dim23>,T> TensorExpr;
    return Tensor3_dg_Expr<TensorExpr,T,Dim23,Dim1,i,j,k>
      (TensorExpr(*this,N0));
  }

  /* The ++ operator increments the pointer, not the number that the
     pointer points to.  This allows iterating over a grid. */

  template<int Current_Position0,int Current_Position1>
  inline void increment(
    const Number<Current_Position0> &,
    const Number<Current_Position1> &
  ) const {
    data[Current_Position0-1][Current_Position1-1]+=inc;
  }

  const Tensor4_ddg<T*,Tensor_Dim01,Tensor_Dim23> & operator++() const
  {
    T4_ddg_increment(*this,Number<(Tensor_Dim01*(Tensor_Dim01+1))/2>(),Number<(Tensor_Dim23*(Tensor_Dim23+1))/2>());
    return *this;
  }

};

class AAA {
	void A() {
		double d[9];
		Tensor4_ddg<double*,3,3> t2(
      &d[0],&d[1],&d[2],
      &d[1],&d[2],&d[3],
      &d[4],&d[5],&d[6]
		);
		++t2;
	}
};
