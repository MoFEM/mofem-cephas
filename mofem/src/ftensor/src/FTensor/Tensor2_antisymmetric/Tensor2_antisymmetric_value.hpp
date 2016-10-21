/* A general version, not for pointers. */

template <class T, int Tensor_Dim>
class Tensor2_antisymmetric
{
  T data[(Tensor_Dim*(Tensor_Dim-1))/2];
public:
  Tensor2_antisymmetric() {}

  /* Tensor_Dim=2 */
  Tensor2_antisymmetric(T d01)
  {
    Tensor2_antisymmetric_constructor<T,Tensor_Dim>(data,d01);
  }

  /* Tensor_Dim=3 */
  Tensor2_antisymmetric(T d01, T d02, T d12)
  {
    Tensor2_antisymmetric_constructor<T,Tensor_Dim>(data,d01,d02,d12);
  }

  /* Tensor_Dim=4 */
  Tensor2_antisymmetric(T d01, T d02, T d03, T d12, T d13, T d23)
  {
    Tensor2_antisymmetric_constructor<T,Tensor_Dim>
      (data,d01,d02,d03,d12,d13,d23);
  }

  /* There are two ways of accessing the values inside,
     unsafe(int,int) and operator(int,int).  unsafe(int,int) will give
     you a wrong answer if you aren't careful.  The problem is that we
     only store the minimal set of components, but some have different
     signs.  We can't return the negative of a component, and assign
     something to it, because that would assign something to a
     temporary.  To get the correct answer if you don't want to change
     the value, just use operator(int,int). */

  T & unsafe(const int N1, const int N2)
  {
#ifdef FTENSOR_DEBUG
    if(N1>=Tensor_Dim || N1<0 || N2>=Tensor_Dim || N2<0 || N1>=N2)
      {
        std::stringstream s;
        s << "Bad index in Tensor2_antisymmetric<T," << Tensor_Dim
          << ">.operator(" << N1 << "," << N2 << ")"
          << std::endl;
        throw std::runtime_error(s.str());
      }
#endif
    return data[(N2-1)+(N1*(2*(Tensor_Dim-1)-N1-1))/2];
  }

  T operator()(const int N1, const int N2) const
  {
#ifdef FTENSOR_DEBUG
    if(N1>=Tensor_Dim || N1<0 || N2>=Tensor_Dim || N2<0)
      {
        std::stringstream s;
        s << "Bad index in Tensor2_antisymmetric<T," << Tensor_Dim
          << ">.operator(" << N1 << "," << N2 << ") const"
          << std::endl;
        throw std::runtime_error(s.str());
      }
#endif
    return N1==N2 ? 0 : (N1<N2 ? data[(N2-1)+(N1*(2*(Tensor_Dim-1)-N1-1))/2]
                         : -data[(N1-1)+(N2*(2*(Tensor_Dim-1)-N2-1))/2]);
  }

  /* These operator()'s are the first part in constructing template
     expressions.  They can be used to slice off lower dimensional
     parts. They are not entirely safe, since you can accidently use a
     higher dimension than what is really allowed (like Dim=5). */

  /* This returns a Tensor2_Expr, since the indices are not really
     antisymmetric anymore since they cover different dimensions. */

  template<char i, char j, int Dim0, int Dim1>
  Tensor2_Expr<Tensor2_antisymmetric<T,Tensor_Dim>,T,Dim0,Dim1,i,j>
  operator()(const Index<i,Dim0> index1, const Index<j,Dim1> index2)
  {
    return Tensor2_Expr<Tensor2_antisymmetric<T,Tensor_Dim>,T,Dim0,Dim1,i,j>
      (*this);
  }

  template<char i, char j, int Dim0, int Dim1>
  Tensor2_Expr<const Tensor2_antisymmetric<T,Tensor_Dim>,T,Dim0,Dim1,i,j>
  operator()(const Index<i,Dim0> index1, const Index<j,Dim1> index2) const
  {
    return Tensor2_Expr<const Tensor2_antisymmetric<T,Tensor_Dim>,T,Dim0,Dim1,i,j>
      (*this);
  }

  /* This returns a Tensor2_antisymmetric_Expr, since the indices are still
     antisymmetric on the lower dimensions. */

  template<char i, char j, int Dim>
  Tensor2_antisymmetric_Expr<Tensor2_antisymmetric<T,Tensor_Dim>,T,Dim,i,j>
  operator()(const Index<i,Dim> index1, const Index<j,Dim> index2)
  {
    return Tensor2_antisymmetric_Expr<Tensor2_antisymmetric<T,Tensor_Dim>,T,Dim,i,j>
      (*this);
  }

  template<char i, char j, int Dim>
  Tensor2_antisymmetric_Expr<const Tensor2_antisymmetric<T,Tensor_Dim>,T,Dim,i,j>
  operator()(const Index<i,Dim> index1, const Index<j,Dim> index2) const
  {
    return Tensor2_antisymmetric_Expr<const Tensor2_antisymmetric<T,Tensor_Dim>,
      T,Dim,i,j>(*this);
  }

  /* This is for expressions where a number is used for one slot, and
     an index for another, yielding a Tensor1_Expr.  The non-const
     versions don't actually create a Tensor2_number_rhs_[01] object.
     They create a Tensor1_Expr directly, which provides the
     appropriate indexing operators.  The const versions do create a
     Tensor2_number_[01]. */

  template<char i, int N, int Dim>
  Tensor1_Expr<Tensor2_number_rhs_1<Tensor2_antisymmetric<T,Tensor_Dim>,T,N>,
    T,Dim,i>
  operator()(const Index<i,Dim> index1, const Number<N> &n1)
  {
    typedef Tensor2_number_rhs_1<Tensor2_antisymmetric<T,Tensor_Dim>,T,N>
      TensorExpr;
    return Tensor1_Expr<TensorExpr,T,Dim,i>(*this);
  }

  template<char i, int N, int Dim>
  Tensor1_Expr<const Tensor2_number_1<const Tensor2_antisymmetric<T,Tensor_Dim>,
    T,N>,T,Dim,i>
  operator()(const Index<i,Dim> index1, const Number<N> &n1) const
  {
    typedef const Tensor2_number_1<const Tensor2_antisymmetric<T,Tensor_Dim>,T,N>
      TensorExpr;
    return Tensor1_Expr<TensorExpr,T,Dim,i>(TensorExpr(*this));
  }

  template<char i, int N, int Dim>
  Tensor1_Expr<Tensor2_number_rhs_0<Tensor2_antisymmetric<T,Tensor_Dim>,T,N>,
    T,Dim,i>
  operator()(const Number<N> &n1, const Index<i,Dim> index1)
  {
    typedef Tensor2_number_rhs_0<Tensor2_antisymmetric<T,Tensor_Dim>,T,N>
      TensorExpr;
    return Tensor1_Expr<TensorExpr,T,Dim,i>(*this);
  }

  template<char i, int N, int Dim>
  Tensor1_Expr<const Tensor2_number_0<const Tensor2_antisymmetric<T,Tensor_Dim>,
    T,N>,T,Dim,i>
  operator()(const Number<N> &n1, const Index<i,Dim> index1) const
  {
    typedef const Tensor2_number_0<const Tensor2_antisymmetric<T,Tensor_Dim>,T,N>
      TensorExpr;
    return Tensor1_Expr<TensorExpr,T,Dim,i>(TensorExpr(*this));
  }

  /* Specializations for using actual numbers instead of Number<> */

  template<char i, int Dim>
  Tensor1_Expr<const Tensor2_numeral_1<const Tensor2_antisymmetric<T,Tensor_Dim>,
    T>,T,Dim,i>
  operator()(const Index<i,Dim> index1, const int N) const
  {
    typedef const Tensor2_numeral_1<const Tensor2_antisymmetric<T,Tensor_Dim>,T>
      TensorExpr;
    return Tensor1_Expr<TensorExpr,T,Dim,i>(TensorExpr(*this,N));
  }

  template<char i, int Dim>
  Tensor1_Expr<const Tensor2_numeral_0<const Tensor2_antisymmetric<T,Tensor_Dim>,
    T>,T,Dim,i>
  operator()(const int N, const Index<i,Dim> index1) const
  {
    typedef const Tensor2_numeral_0<const Tensor2_antisymmetric<T,Tensor_Dim>,T>
      TensorExpr;
    return Tensor1_Expr<TensorExpr,T,Dim,i>(TensorExpr(*this,N));
  }
  
  /* These two operator()'s return the Tensor2 with internal
     contractions, yielding a T.  I have to specify one for both
     const and non-const because otherwise the compiler will use the
     operator() which gives a Tensor2_Expr<>. */
  
  template<char i, int Dim>
  T operator()(const Index<i,Dim> index1, const Index<i,Dim> index2)
  {
    return 0;
  }

  template<char i, int Dim>
  T operator()(const Index<i,Dim> index1, const Index<i,Dim> index2) const
  {
    return 0;
  }
};
